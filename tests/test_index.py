"""Data-free unit tests for the precursor index.

These build synthetic ``PrecursorIndex`` objects in memory, so they run without
any mass-spectrometry data files. They pin the bug fixes ported from the
upstream cleanup:

* ``find(..., ms_level != 2)`` returns ``[]`` (only MS2 is indexed).
* ``CrossFilePrecursorIndex.find`` populates ``precursor_mz`` / ``retention_time``
  (previously left as ``None``) and sorts results.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from pyx500r.index import (
    CrossFilePrecursorIndex,
    PrecursorIndex,
)


def _make_index(file_path: str = "synthetic.wiff2") -> PrecursorIndex:
    return PrecursorIndex(
        file_path=file_path,
        precursor_mz=np.array([100.0, 200.0, 200.05, 300.0], dtype=np.float64),
        indices=[(0, 1, 0), (0, 1, 1), (0, 1, 2), (0, 1, 3)],
        retention_times=np.array([1.0, 2.0, 2.1, 3.0], dtype=np.float64),
        n_ms2=4,
    )


class TestPrecursorIndex(unittest.TestCase):
    def setUp(self) -> None:
        self.index = _make_index()

    def test_find_within_tolerance(self) -> None:
        hits = self.index.find(200.0, tolerance_da=0.1)
        self.assertEqual(set(hits), {(0, 1, 1), (0, 1, 2)})

    def test_find_no_match(self) -> None:
        self.assertEqual(self.index.find(999.0, tolerance_da=0.1), [])

    def test_find_non_ms2_returns_empty(self) -> None:
        self.assertEqual(self.index.find(200.0, tolerance_da=0.1, ms_level=1), [])
        self.assertGreater(len(self.index.find(200.0, tolerance_da=0.1, ms_level=2)), 0)

    def test_to_from_dict_round_trip(self) -> None:
        restored = PrecursorIndex.from_dict(self.index.to_dict())
        self.assertEqual(restored.file_path, self.index.file_path)
        self.assertEqual(restored.n_ms2, self.index.n_ms2)
        self.assertTrue(np.allclose(restored.precursor_mz, self.index.precursor_mz))
        self.assertEqual(restored.indices, self.index.indices)


class TestCrossFilePrecursorIndex(unittest.TestCase):
    def setUp(self) -> None:
        self.cross = CrossFilePrecursorIndex()
        self.cross.add(_make_index("a.wiff2"))
        self.cross.add(_make_index("b.wiff2"))

    def test_len_and_repr(self) -> None:
        self.assertEqual(len(self.cross), 2)
        self.assertIn("files=2", repr(self.cross))

    def test_find_populates_all_fields(self) -> None:
        results = self.cross.find(200.0, precursor_tolerance_da=0.1)
        self.assertGreater(len(results), 0)
        for r in results:
            # Regression: these were previously left as None.
            self.assertIsNotNone(r["precursor_mz"])
            self.assertIsNotNone(r["retention_time"])
            self.assertIsInstance(r["precursor_mz"], float)
            self.assertIsInstance(r["retention_time"], float)
            self.assertIn(r["file_path"], ("a.wiff2", "b.wiff2"))

    def test_find_sorted(self) -> None:
        results = self.cross.find(250.0, precursor_tolerance_da=200.0)
        keys = [(r["file_path"], r["precursor_mz"]) for r in results]
        self.assertEqual(keys, sorted(keys))

    def test_find_non_ms2_empty(self) -> None:
        self.assertEqual(
            self.cross.find(200.0, precursor_tolerance_da=0.1, ms_level=1), []
        )

    def test_save_load_round_trip(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "idx.json"
            self.cross.save(p)
            restored = CrossFilePrecursorIndex.load(p)
            self.assertEqual(len(restored), len(self.cross))


if __name__ == "__main__":
    unittest.main(verbosity=2)
