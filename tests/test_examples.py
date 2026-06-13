"""Data-free tests for the example apps under ``examples/``.

These assert the apps and helpers import and build without any
mass-spectrometry data files. The data-driven endpoint tests from upstream are
intentionally omitted because they require ``.wiff2`` / ``.qsession`` /
LibraryView ``.sqlite`` files, which are never shipped with this package.
"""

from __future__ import annotations

import os
import unittest
from datetime import datetime

try:
    import fastapi  # noqa: F401

    _HAVE_FASTAPI = True
except Exception:
    _HAVE_FASTAPI = False


class TestSerializers(unittest.TestCase):
    """The serialization helpers must be data-independent and JSON-safe."""

    def test_jsonify_datetime_and_arrays(self) -> None:
        import numpy as np

        from examples.server.serializers import _jsonify, _to_list

        self.assertEqual(_jsonify(datetime(2026, 1, 2, 3, 4, 5)), "2026-01-02T03:04:05")
        self.assertEqual(_to_list(np.array([1.0, 2.0])), [1.0, 2.0])
        self.assertEqual(_to_list([1.0, 2.0]), [1.0, 2.0])
        self.assertEqual(_to_list(None), [])

    def test_jsonify_nested_and_fallback(self) -> None:
        from examples.server.serializers import _jsonify

        out = _jsonify({"a": [1, 2], "b": {"c": True}})
        self.assertEqual(out, {"a": [1, 2], "b": {"c": True}})
        # unknown object types are stringified rather than raising
        class Weird:
            def __str__(self) -> str:
                return "weird"

        self.assertEqual(_jsonify(Weird()), "weird")


@unittest.skipUnless(_HAVE_FASTAPI, "fastapi not installed")
class TestExampleAppsBuild(unittest.TestCase):
    """Both FastAPI apps must construct without opening any data file."""

    def test_server_app_builds(self) -> None:
        os.environ["PYX500R_DATA_ROOT"] = "/tmp/pyx500r-nonexistent"
        import importlib

        from examples.server import app as app_module

        importlib.reload(app_module)
        self.assertTrue(app_module.app.title)

    def test_library_browser_app_builds(self) -> None:
        # DB is opened lazily, so building the app must not require the file.
        os.environ["PYX500R_LIBRARY_DB"] = "/tmp/pyx500r-nonexistent.sqlite"
        import importlib

        from examples.library_browser import app as app_module

        importlib.reload(app_module)
        self.assertTrue(app_module.app.title)


if __name__ == "__main__":
    unittest.main(verbosity=2)
