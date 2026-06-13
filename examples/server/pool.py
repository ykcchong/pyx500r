"""Thread-safe, path-keyed pool of pyx500r readers.

pyx500r readers wrap an **in-memory SQLite connection** created with the default
``check_same_thread=True``. Such a connection may only be touched by the thread
that created it — a per-reader ``Lock`` is *not* sufficient, because serializing
access still lets *different* threads touch the same connection.

The robust pattern (and the one used here) is to give every reader its own
dedicated single-worker thread. Both the ``open_*`` call and every subsequent
method call are submitted to that thread, so the connection is always used by
its creator. Callers never touch the reader directly; they pass a function::

    pool = ReaderPool(max_open=16)
    samples = pool.with_wiff(path, lambda r: r.list_samples())

The submitted callable runs in the reader's thread; its return value is passed
back to the caller. Combine with ``fastapi.concurrency.run_in_threadpool`` so
the event loop itself never blocks.

See ../../docs/GUI_INTEGRATION.md §4.
"""

from __future__ import annotations

import threading
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Callable, TypeVar

from pyx500r import QSessionReader, WiffReader, open_qsession, open_wiff2

T = TypeVar("T")


class _Entry:
    """A reader bound to its own dedicated single-worker thread."""

    def __init__(self, opener: Callable[[], object], key: str) -> None:
        self._executor = ThreadPoolExecutor(
            max_workers=1, thread_name_prefix=f"pyx500r:{Path(key).name}"
        )
        # Open the reader *inside* the dedicated thread so the SQLite
        # connection is owned by that thread.
        self.reader = self._executor.submit(opener).result()

    def call(self, fn: Callable[[object], T]) -> T:
        """Run ``fn(reader)`` in the reader's dedicated thread and return result."""
        return self._executor.submit(fn, self.reader).result()

    def close(self) -> None:
        def _close(r: object) -> None:
            close = getattr(r, "close", None)
            if callable(close):
                close()

        try:
            self._executor.submit(_close, self.reader).result()
        finally:
            self._executor.shutdown(wait=True)


class ReaderPool:
    """LRU cache of dedicated-thread readers, keyed by resolved file path."""

    def __init__(self, max_open: int = 16) -> None:
        self._lock = threading.Lock()
        self._wiff: "OrderedDict[str, _Entry]" = OrderedDict()
        self._qs: "OrderedDict[str, _Entry]" = OrderedDict()
        self._max = max_open

    # -- public API ----------------------------------------------------------
    def with_wiff(self, path: str | Path, fn: Callable[[WiffReader], T]) -> T:
        """Run ``fn(reader)`` against the WiffReader for *path*."""
        return self._entry(self._wiff, path, lambda: open_wiff2(path)).call(fn)  # type: ignore[arg-type]

    def with_qsession(self, path: str | Path, fn: Callable[[QSessionReader], T]) -> T:
        """Run ``fn(reader)`` against the QSessionReader for *path*."""
        return self._entry(self._qs, path, lambda: open_qsession(path)).call(fn)  # type: ignore[arg-type]

    # -- internals -----------------------------------------------------------
    def _entry(
        self,
        od: "OrderedDict[str, _Entry]",
        path: str | Path,
        opener: Callable[[], object],
    ) -> _Entry:
        key = str(Path(path).resolve())
        with self._lock:
            entry = od.get(key)
            if entry is None:
                entry = _Entry(opener, key)
                od[key] = entry
                self._evict(od)
            od.move_to_end(key)
            return entry

    def _evict(self, od: "OrderedDict[str, _Entry]") -> None:
        """Evict least-recently-used readers. Caller must hold ``self._lock``."""
        while len(od) > self._max:
            _key, entry = od.popitem(last=False)
            try:
                entry.close()
            except Exception:
                pass

    def close_all(self) -> None:
        with self._lock:
            for od in (self._wiff, self._qs):
                for _key, entry in od.items():
                    try:
                        entry.close()
                    except Exception:
                        pass
                od.clear()
