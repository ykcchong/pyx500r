"""Data-free CLI contract tests.

These do not touch any mass-spectrometry data files — they only assert the
importable/argv contract of the CLI entry points (every ``main`` accepts an
``argv`` list and ``--help`` exits 0). The data-driven CLI behaviour tests from
upstream are intentionally omitted because they require ``.wiff2`` files.
"""

from __future__ import annotations

import inspect
import io
import tomllib
import unittest
from contextlib import redirect_stdout
from pathlib import Path


class TestCliArgvContracts(unittest.TestCase):
    def test_main_functions_accept_argv(self) -> None:
        from pyx500r import cli, cli_bridge, cli_parallel, libsearch_cli, w2searcher

        for mod in (cli, cli_bridge, cli_parallel, libsearch_cli, w2searcher):
            with self.subTest(module=mod.__name__):
                sig = inspect.signature(mod.main)
                self.assertIn("argv", sig.parameters)

    def test_cli_help_exits_zero(self) -> None:
        from pyx500r import cli

        buf = io.StringIO()
        with self.assertRaises(SystemExit) as ctx, redirect_stdout(buf):
            cli.main(["--help"])
        self.assertEqual(ctx.exception.code, 0)

    def test_libsearch_help_exits_zero(self) -> None:
        from pyx500r import libsearch_cli

        buf = io.StringIO()
        with self.assertRaises(SystemExit) as ctx, redirect_stdout(buf):
            libsearch_cli.main(["--help"])
        self.assertEqual(ctx.exception.code, 0)

    def test_short_console_script_aliases_are_declared(self) -> None:
        pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
        scripts = tomllib.loads(pyproject.read_text())["project"]["scripts"]

        expected = {
            "x500r": "pyx500r.cli:main",
            "x500rp": "pyx500r.cli_parallel:main",
            "x500rqsession": "pyx500r.cli_bridge:main",
            "x500rindex": "pyx500r.index_builder:main",
            "x500rlibsearch": "pyx500r.libsearch_cli:main",
            "x500rsearch": "pyx500r.w2searcher:main",
            "x500rgui": "pyx500r.wiff_gui:main",
        }

        for name, target in expected.items():
            with self.subTest(script=name):
                self.assertEqual(scripts.get(name), target)


if __name__ == "__main__":
    unittest.main(verbosity=2)
