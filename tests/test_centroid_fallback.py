"""Tests for the pure-Python fallback centroiding path.

These ensure that pyx500r.centroid_fallback works correctly even when
numba is not installed.
"""

from __future__ import annotations

import unittest

from pyx500r.centroid_fallback import (
    Peak,
    add_framing_zeros,
    centroid_spectrum,
    moving_average_smooth,
)


class TestFallbackMovingAverageSmooth(unittest.TestCase):
    def test_empty_returns_empty(self) -> None:
        self.assertEqual(moving_average_smooth([], 1), [])

    def test_triangle_wave(self) -> None:
        y = [1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0]
        result = moving_average_smooth(y, 1)
        expected = [
            1.3333333333333333, 2.0, 3.0, 4.0, 4.333333333333333,
            4.0, 3.0, 2.0, 1.3333333333333333,
        ]
        self.assertEqual(len(result), len(y))
        for a, b in zip(result, expected):
            self.assertAlmostEqual(a, b, places=12)


class TestFallbackAddFramingZeros(unittest.TestCase):
    def test_single_point(self) -> None:
        mz, ints = add_framing_zeros([100.0], [10.0], lambda x: 1.0)
        self.assertEqual(mz, [99.0, 100.0, 101.0])
        self.assertEqual(ints, [0.0, 10.0, 0.0])


class TestFallbackCentroidSpectrumSynthetic(unittest.TestCase):
    def test_single_peak(self) -> None:
        mz = [100.0, 100.5, 101.0]
        ints = [1.0, 10.0, 1.0]
        cmz, cint = centroid_spectrum(mz, ints)
        self.assertEqual(len(cmz), 1)
        self.assertAlmostEqual(cmz[0], 100.5, places=1)

    def test_output_sorted_by_mz(self) -> None:
        mz = [100.0, 101.0, 102.0, 103.0, 104.0]
        ints = [10.0, 1.0, 10.0, 1.0, 10.0]
        cmz, cint = centroid_spectrum(mz, ints)
        self.assertTrue(all(a < b for a, b in zip(cmz, cmz[1:])))

    def test_return_arrays(self) -> None:
        import numpy as np
        mz = [100.0, 100.5, 101.0]
        ints = [1.0, 10.0, 1.0]
        cmz, cint = centroid_spectrum(mz, ints, return_arrays=True)
        self.assertIsInstance(cmz, np.ndarray)
        self.assertIsInstance(cint, np.ndarray)


if __name__ == "__main__":
    unittest.main(verbosity=2)
