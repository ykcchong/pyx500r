"""Baseline and regression tests for pyx500r.centroid.

These tests establish expected outputs for the current centroiding
implementation.  After refactoring, they must continue to pass with
bit-identical (or near-identical) results.
"""

from __future__ import annotations

import unittest

from pyx500r.centroid import (
    Peak,
    add_framing_zeros,
    centroid_spectrum,
    moving_average_smooth,
)


# ---------------------------------------------------------------------------
# 1. moving_average_smooth
# ---------------------------------------------------------------------------

class TestMovingAverageSmooth(unittest.TestCase):
    """Unit tests for the moving-average smoothing filter."""

    def test_empty_returns_empty(self) -> None:
        self.assertEqual(moving_average_smooth([], 1), [])

    def test_half_window_zero_is_identity(self) -> None:
        y = [1.0, 2.0, 3.0, 4.0]
        self.assertEqual(moving_average_smooth(y, 0), y)

    def test_single_element(self) -> None:
        self.assertEqual(moving_average_smooth([5.0], 1), [5.0])

    def test_no_smoothing_when_n_lt_window(self) -> None:
        y = [1.0, 2.0]
        self.assertEqual(moving_average_smooth(y, 2), y)

    def test_triangle_wave(self) -> None:
        y = [1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0]
        result = moving_average_smooth(y, 1)
        expected = [
            1.3333333333333333,
            2.0,
            3.0,
            4.0,
            4.333333333333333,
            4.0,
            3.0,
            2.0,
            1.3333333333333333,
        ]
        self.assertEqual(len(result), len(y))
        for a, b in zip(result, expected):
            self.assertAlmostEqual(a, b, places=12)

    def test_constant_signal_unchanged(self) -> None:
        y = [5.0] * 20
        result = moving_average_smooth(y, 2)
        self.assertTrue(all(abs(v - 5.0) < 1e-12 for v in result))

    def test_larger_half_window(self) -> None:
        y = [1.0, 0.0, 0.0, 0.0, 1.0]
        result = moving_average_smooth(y, 2)
        # window=5; edges clamped, middle = average of all 5 = 0.4
        self.assertAlmostEqual(result[2], 0.4, places=12)


# ---------------------------------------------------------------------------
# 2. add_framing_zeros
# ---------------------------------------------------------------------------

class TestAddFramingZeros(unittest.TestCase):
    """Unit tests for zero-point insertion."""

    def test_empty_returns_empty(self) -> None:
        mz, ints = add_framing_zeros([], [], lambda x: 1.0)
        self.assertEqual(mz, [])
        self.assertEqual(ints, [])

    def test_single_point(self) -> None:
        mz, ints = add_framing_zeros([100.0], [10.0], lambda x: 1.0)
        self.assertEqual(mz, [99.0, 100.0, 101.0])
        self.assertEqual(ints, [0.0, 10.0, 0.0])

    def test_gap_insertion(self) -> None:
        mz = [100.0, 101.0, 103.0]
        ints = [10.0, 20.0, 30.0]
        mz_out, ints_out = add_framing_zeros(mz, ints, lambda x: 1.0, half_insert=1)
        self.assertEqual(mz_out, [99.0, 100.0, 101.0, 102.0, 103.0, 104.0])
        self.assertEqual(ints_out, [0.0, 10.0, 20.0, 0.0, 30.0, 0.0])

    def test_no_insert_before_first(self) -> None:
        mz = [100.0, 101.0]
        ints = [10.0, 20.0]
        mz_out, ints_out = add_framing_zeros(
            mz, ints, lambda x: 1.0, insert_before_first=False, insert_after_last=False
        )
        self.assertEqual(mz_out, [100.0, 101.0])
        self.assertEqual(ints_out, [10.0, 20.0])


# ---------------------------------------------------------------------------
# 3. centroid_spectrum — synthetic / small data
# ---------------------------------------------------------------------------

class TestCentroidSpectrumSynthetic(unittest.TestCase):
    """Tests against hand-crafted spectra."""

    def test_too_few_points_returns_nonzero(self) -> None:
        mz = [100.0, 101.0]
        ints = [0.0, 5.0]
        cmz, cint = centroid_spectrum(mz, ints)
        self.assertEqual(cmz, [101.0])
        self.assertEqual(cint, [5.0])

    def test_all_zeros_returns_empty(self) -> None:
        mz = [100.0, 101.0, 102.0]
        ints = [0.0, 0.0, 0.0]
        cmz, cint = centroid_spectrum(mz, ints)
        self.assertEqual(cmz, [])
        self.assertEqual(cint, [])

    def test_single_peak(self) -> None:
        # A simple Gaussian-like peak centered at 100.5
        mz = [100.0, 100.5, 101.0]
        ints = [1.0, 10.0, 1.0]
        cmz, cint = centroid_spectrum(mz, ints)
        self.assertEqual(len(cmz), 1)
        self.assertEqual(len(cint), 1)
        # Centroid should be close to 100.5
        self.assertAlmostEqual(cmz[0], 100.5, places=1)
        # Area should be sum of intensities above threshold
        self.assertGreater(cint[0], 0.0)

    def test_two_separate_peaks(self) -> None:
        # Dense points with a deep valley between two peaks
        mz = [100.0, 100.1, 100.2, 100.3, 100.4, 100.5,
              101.0, 101.1, 101.2, 101.3, 101.4, 101.5]
        ints = [1.0, 5.0, 20.0, 50.0, 20.0, 5.0,
                5.0, 20.0, 50.0, 20.0, 5.0, 1.0]
        cmz, cint = centroid_spectrum(mz, ints)
        self.assertEqual(len(cmz), 2)
        self.assertTrue(cmz[0] < cmz[1])

    def test_centroid_percentage_100(self) -> None:
        mz = [100.0, 100.5, 101.0]
        ints = [1.0, 10.0, 1.0]
        cmz, cint = centroid_spectrum(mz, ints, centroid_percentage=100.0)
        # Only the apex should contribute at 100%
        self.assertEqual(len(cmz), 1)
        self.assertAlmostEqual(cmz[0], 100.5, places=3)

    def test_output_sorted_by_mz(self) -> None:
        mz = [100.0, 101.0, 102.0, 103.0, 104.0]
        ints = [10.0, 1.0, 10.0, 1.0, 10.0]
        cmz, cint = centroid_spectrum(mz, ints)
        self.assertTrue(all(a < b for a, b in zip(cmz, cmz[1:])))


# ---------------------------------------------------------------------------
# 4. Peak dataclass
# ---------------------------------------------------------------------------

class TestPeakDataclass(unittest.TestCase):
    def test_defaults(self) -> None:
        p = Peak(x_value=100.0, area=50.0, height=10.0,
                 apex_x=100.0, apex_y=10.0,
                 start_x=99.0, start_y=1.0,
                 end_x=101.0, end_y=1.0)
        self.assertEqual(p.noise, 0.0)
        self.assertEqual(p.start_x_half_height, 0.0)
        self.assertEqual(p.end_x_half_height, 0.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
