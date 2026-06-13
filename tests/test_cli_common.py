"""Unit tests for the shared CLI helpers in ``pyx500r._cli_common``.

These also guard against regressions from consolidating the previously
duplicated copies that lived in ``cli.py`` and ``w2searcher.py``.
"""

from __future__ import annotations

import unittest

from pyx500r._cli_common import parse_transitions, ppm_tolerance


class TestPpmTolerance(unittest.TestCase):
    def test_basic_conversion(self) -> None:
        # 20 ppm of 500 m/z == 0.01 Da
        self.assertAlmostEqual(ppm_tolerance(500.0, 20.0), 0.01, places=12)

    def test_zero_ppm(self) -> None:
        self.assertEqual(ppm_tolerance(500.0, 0.0), 0.0)

    def test_scales_with_mz(self) -> None:
        self.assertAlmostEqual(
            ppm_tolerance(1000.0, 10.0), 2 * ppm_tolerance(500.0, 10.0), places=12
        )


class TestParseTransitions(unittest.TestCase):
    def test_none_returns_empty(self) -> None:
        self.assertEqual(parse_transitions(None), [])

    def test_empty_list_returns_empty(self) -> None:
        self.assertEqual(parse_transitions([]), [])

    def test_single_transition(self) -> None:
        result = parse_transitions(["250.1587:191.0857,163.0907,109.0443"])
        self.assertEqual(len(result), 1)
        prec, prods = result[0]
        self.assertAlmostEqual(prec, 250.1587)
        self.assertEqual(prods, [191.0857, 163.0907, 109.0443])

    def test_multiple_transitions(self) -> None:
        result = parse_transitions(["100.0:50.0", "200.0:75.0,80.0"])
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0], (100.0, [50.0]))
        self.assertEqual(result[1], (200.0, [75.0, 80.0]))

    def test_whitespace_is_stripped(self) -> None:
        result = parse_transitions(["250.0: 191.0 , 163.0 "])
        self.assertEqual(result[0], (250.0, [191.0, 163.0]))

    def test_missing_colon_raises(self) -> None:
        with self.assertRaises(ValueError):
            parse_transitions(["250.0,191.0"])

    def test_no_products_raises(self) -> None:
        with self.assertRaises(ValueError):
            parse_transitions(["250.0:"])

    def test_cli_and_w2searcher_share_impl(self) -> None:
        """Both CLIs must reference the exact same helper objects."""
        from pyx500r import cli, w2searcher

        self.assertIs(cli._parse_transitions, parse_transitions)
        self.assertIs(cli._ppm_tolerance, ppm_tolerance)
        self.assertIs(w2searcher._parse_transitions, parse_transitions)
        self.assertIs(w2searcher._ppm_tolerance, ppm_tolerance)


if __name__ == "__main__":
    unittest.main(verbosity=2)
