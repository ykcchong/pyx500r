"""Pure-Python fallback centroiding and spectrum processing for SCIEX WIFF2 data.

This module has **no numba dependency**. It provides a pure-Python + NumPy
implementation of the centroiding algorithm that is bit-exact with the
numba-accelerated fast path but 5-10× slower.

Public API
----------
* ``centroid_spectrum`` — convert profile to centroid format
* ``add_framing_zeros`` — insert zero-intensity boundary points
* ``moving_average_smooth`` — simple moving-average filter
* ``Peak`` — dataclass for a single centroided peak
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from bisect import bisect_left


# --------------------------------------------------------------------------- #
# Data structures
# --------------------------------------------------------------------------- #

@dataclass
class Peak:
    """A single centroided peak.

    Attributes mirror ``Clearcore2.RawXYProcessing.PeakClass``.
    """
    x_value: float   # centroid m/z
    area: float
    height: float
    apex_x: float
    apex_y: float
    start_x: float
    start_y: float
    end_x: float
    end_y: float
    start_x_half_height: float = 0.0
    end_x_half_height: float = 0.0
    noise: float = 0.0


# --------------------------------------------------------------------------- #
# Moving average smooth (MovingAverageSmooth)
# --------------------------------------------------------------------------- #

def moving_average_smooth(y, half_window: int):
    """Simple moving-average smoothing.

    Parameters
    ----------
    y: intensity values (list or 1-D array).
    half_window: radius of the smoothing window (total width = 2*half_window+1).

    Returns
    -------
    Smoothed intensity values, same type and length as input.
    """
    arr = np.asarray(y, dtype=np.float64)
    n = arr.size
    if n == 0 or half_window == 0:
        result = arr.copy()
        return result.tolist() if isinstance(y, list) else result

    window = 2 * half_window + 1
    if n < window:
        result = arr.copy()
        return result.tolist() if isinstance(y, list) else result

    out = np.empty(n, dtype=np.float64)

    for j in range(half_window + 1):
        s = 0.0
        for k in range(j - half_window, j + half_window + 1):
            idx = 0 if k < 0 else k
            s += arr[idx]
        out[j] = s / window

    for j in range(half_window + 1, n - half_window):
        s = 0.0
        for k in range(j - half_window, j + half_window + 1):
            s += arr[k]
        out[j] = s / window

    for j in range(n - half_window, n):
        s = 0.0
        for k in range(j - half_window, j + half_window + 1):
            idx = n - 1 if k >= n else k
            s += arr[idx]
        out[j] = s / window

    return out.tolist() if isinstance(y, list) else out


# --------------------------------------------------------------------------- #
# Add zeros (AddZeros)
# --------------------------------------------------------------------------- #

def add_framing_zeros(
    mz: list[float],
    intensities: list[float],
    step_size_fn,
    half_insert: int = 1,
    insert_before_first: bool = True,
    insert_after_last: bool = True,
) -> tuple[list[float], list[float]]:
    """Insert zero-intensity points at spectrum boundaries.

    Port of ``AddZeros.AddSomeMissingZeros``.
    """
    n = len(mz)
    if n == 0:
        return [], []

    x_out: list[float] = []
    y_out: list[float] = []

    if half_insert == 1:
        _add_one_missing_zero(
            mz, intensities, step_size_fn,
            insert_before_first, insert_after_last,
            x_out, y_out,
        )
    else:
        for i in range(n):
            _half_insert_left(mz, intensities, step_size_fn, i,
                              insert_before_first, x_out, y_out)
            x_out.append(mz[i])
            y_out.append(intensities[i])
            _half_insert_right(mz, intensities, step_size_fn, i,
                               insert_after_last, x_out, y_out)

    return x_out, y_out


def _add_one_missing_zero(
    mz: list[float], intensities: list[float],
    step_size_fn,
    insert_before_first: bool, insert_after_last: bool,
    x_out: list[float], y_out: list[float],
) -> None:
    """Fast path for half_insert == 1."""
    n = len(mz)
    if n == 0:
        return

    if insert_before_first:
        step = step_size_fn(mz[0])
        x_out.append(mz[0] - step)
        y_out.append(0.0)

    prev_x = mz[0]
    num_steps = 1
    next_insert_at = mz[0] - 1.0

    for i in range(n - 1):
        cur_x = mz[i]
        if cur_x >= next_insert_at:
            step = step_size_fn(cur_x)
            next_insert_at = cur_x + num_steps * step

        if cur_x - step_size_fn(cur_x) > prev_x + step_size_fn(prev_x) / 2.0:
            insert_x = cur_x - step_size_fn(cur_x)
            x_out.append(insert_x)
            y_out.append(0.0)

        x_out.append(cur_x)
        y_out.append(intensities[i])
        prev_x = cur_x

        cur_x_next = mz[i + 1]
        if cur_x_next - step_size_fn(cur_x_next) > prev_x + step_size_fn(prev_x) / 2.0:
            insert_x = prev_x + step_size_fn(prev_x)
            x_out.append(insert_x)
            y_out.append(0.0)
            prev_x = insert_x

    last_x = mz[n - 1]
    step = step_size_fn(last_x)
    if last_x - step > prev_x + step_size_fn(prev_x) / 2.0:
        x_out.append(last_x - step)
        y_out.append(0.0)
    x_out.append(last_x)
    y_out.append(intensities[n - 1])

    if insert_after_last:
        step = step_size_fn(last_x)
        x_out.append(last_x + step)
        y_out.append(0.0)


def _half_insert_left(
    mz: list[float], intensities: list[float],
    step_size_fn, old_index: int, insert_before_first: bool,
    x_out: list[float], y_out: list[float],
) -> None:
    if old_index == 0 and not insert_before_first:
        return

    prev_x = x_out[-1] if x_out else 0.0
    cur_x = mz[old_index]
    step = step_size_fn(cur_x)

    for _ in range(1):
        cur_x -= step
        step = step_size_fn(cur_x)
        if x_out and cur_x < prev_x + step / 2.0:
            cur_x += step
            break
        x_out.append(cur_x)
        y_out.append(0.0)


def _half_insert_right(
    mz: list[float], intensities: list[float],
    step_size_fn, old_index: int, insert_after_last: bool,
    x_out: list[float], y_out: list[float],
) -> None:
    n = len(mz)
    if old_index == n - 1 and not insert_after_last:
        return

    next_x = mz[old_index + 1] if old_index < n - 1 else None
    cur_x = mz[old_index]
    step = step_size_fn(cur_x)

    for _ in range(1):
        cur_x += step
        step = step_size_fn(cur_x)
        if next_x is not None and cur_x > next_x - step / 2.0:
            cur_x -= step
            break
        x_out.append(cur_x)
        y_out.append(0.0)


# --------------------------------------------------------------------------- #
# Pure-Python fallback helpers (list-based)
# --------------------------------------------------------------------------- #

def _find_point_at_height(
    x: list[float], y: list[float], index: int,
    height_fraction: float, look_right: bool,
    peak_indices: list[int],
) -> int:
    """Walk left/right from a peak apex until 85% height or assigned territory."""
    step = _estimate_step(x, index)
    threshold = y[index] * height_fraction
    direction = 1 if look_right else -1
    num_points = len(x)
    num = index

    while (0 < num < num_points - 1
           and peak_indices[num] == -1
           and y[num] > threshold
           and abs(x[num + direction] - x[num]) < step * 1.5):
        num += direction

    if peak_indices[num] != -1:
        if y[num] <= threshold:
            num -= direction
        else:
            peak_id = peak_indices[num]
            max_idx = num
            max_val = y[num]
            scan = num + direction
            while 0 <= scan < num_points and peak_indices[scan] == peak_id:
                if y[scan] > max_val:
                    max_val = y[scan]
                    max_idx = scan
                scan += direction
            if abs(x[max_idx] - x[index]) > 0.75:
                num -= direction

    return num


def _get_peak_centre_range(
    x: list[float], y: list[float], index: int,
    look_right: bool, peak_indices: list[int],
    centroid_fraction: float, step_size_fn,
) -> int:
    """Find the centre-of-range boundary for an unassigned peak."""
    if centroid_fraction < 0.4999:
        height_frac = 0.25
    elif centroid_fraction <= 0.75:
        height_frac = 0.5
    else:
        height_frac = 0.7

    rough_width = _guess_rough_peak_width(step_size_fn, x[index])
    step = step_size_fn(x[index])
    x_center = x[index]
    y_peak = y[index]
    threshold = height_frac * y_peak
    direction = 1 if look_right else -1
    num_points = len(x)

    num = -1
    candidate = index + direction

    while (0 < candidate < num_points - 1
           and y[candidate] > threshold
           and peak_indices[candidate] == -1
           and y[candidate] <= y_peak
           and abs(x[candidate] - x_center) < rough_width
           and abs(x[candidate - direction] - x[candidate]) < 1.5 * step):
        if y[candidate] < threshold * 1.33:
            y_peak *= 0.75
        if num == -1 or y[candidate] < y[num]:
            num = candidate
        candidate += direction

    if y[candidate] > threshold or abs(x[candidate - direction] - x[candidate]) > 1.5 * step:
        if num > 0:
            return num
        return candidate - direction

    return candidate


def _expand_peak_range(
    x: list[float], y: list[float],
    left: int, right: int,
    peak_indices: list[int],
    step_size_fn,
) -> tuple[int, int]:
    """Expand peak boundaries to include trailing/leading ramps."""
    num_points = len(x)

    step = step_size_fn(x[left])
    left -= 1
    while (left >= 0 and peak_indices[left] == -1
           and y[left] <= y[left + 1]
           and abs(x[left + 1] - x[left] - step) < step / 2.0):
        left -= 1
    left += 1
    while left < num_points and y[left] == 0.0 and peak_indices[left] == -1:
        left += 1

    step = step_size_fn(x[right])
    right += 1
    while (right < num_points and peak_indices[right] == -1
           and y[right - 1] >= y[right]
           and abs(x[right] - x[right - 1] - step) < step / 2.0):
        right += 1
    right -= 1
    while right >= 0 and y[right] == 0.0 and peak_indices[right] == -1:
        right -= 1

    return left, right


def _guess_rough_peak_width(step_size_fn, x_val: float) -> float:
    """Estimate peak width based on step size."""
    step = step_size_fn(x_val)
    if step < 0.031:
        return 0.2
    if step < 0.061:
        return 0.3
    if step < 0.21:
        return 0.8
    if step < 0.51:
        return 1.0
    return 4.0 * step


def _estimate_step(x: list[float], index: int) -> float:
    """Estimate local step size around a point."""
    n = len(x)
    if n < 2:
        return 1.0
    if index > 0 and index < n - 1:
        return (x[index + 1] - x[index - 1]) / 2.0
    return x[1] - x[0]


def _set_peak_range(peak_indices: list[int], start: int, end: int, peak_id: int) -> None:
    for i in range(start, end + 1):
        peak_indices[i] = peak_id


def _closest_point(mz: list[float], x_value: float, mode: str) -> int:
    """Binary search for closest point in sorted m/z array."""
    n = len(mz)
    if n == 0:
        return -1
    if n == 1:
        if mode in ("closest", "lower_or_first", "higher_or_last"):
            return 0
        if mode == "lower" and mz[0] <= x_value:
            return 0
        if mode == "higher" and mz[0] >= x_value:
            return 0
        return -1

    lo, hi = -1, n
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if x_value > mz[mid]:
            lo = mid
        else:
            hi = mid

    if lo == -1:
        if mode != "lower":
            return 0
        return 0 if mz[0] == x_value else -1

    if lo >= n - 1:
        if mode != "higher":
            return n - 1
        return n - 1 if mz[n - 1] == x_value else -1

    if x_value == mz[lo]:
        return lo
    if x_value == mz[lo + 1]:
        return lo + 1

    if mode == "closest":
        if abs(mz[lo] - x_value) > abs(mz[lo + 1] - x_value):
            return lo + 1
        return lo
    if mode in ("higher", "higher_or_last"):
        return lo + 1
    return lo


# --------------------------------------------------------------------------- #
# Peak range assignment (pure Python)
# --------------------------------------------------------------------------- #

def _assign_peak_ranges_py(
    x_smooth: np.ndarray,
    y_smooth: np.ndarray,
    local_max_indices: np.ndarray,
    sorted_order: np.ndarray,
    peak_indices: np.ndarray,
    centroid_fraction: float,
    step_size_fn,
) -> None:
    """Pure-Python peak range assignment."""
    peak_indices_list = peak_indices.tolist()
    x_list = x_smooth.tolist()
    y_list = y_smooth.tolist()

    for order_idx in sorted_order:
        apex = int(local_max_indices[order_idx])
        if peak_indices_list[apex] != -1:
            continue

        left = _find_point_at_height(x_list, y_list, apex, 0.85, False, peak_indices_list)
        right = _find_point_at_height(x_list, y_list, apex, 0.85, True, peak_indices_list)

        left_assigned = peak_indices_list[left]
        right_assigned = peak_indices_list[right]

        if left_assigned != -1 and right_assigned != -1:
            if (x_list[apex] - x_list[left]) < (x_list[right] - x_list[apex]):
                _set_peak_range(peak_indices_list, left, right - 1, left_assigned)
            else:
                _set_peak_range(peak_indices_list, left + 1, right, right_assigned)
        elif left_assigned != -1:
            _set_peak_range(peak_indices_list, left, right, left_assigned)
        elif right_assigned != -1:
            _set_peak_range(peak_indices_list, left, right, right_assigned)
        else:
            left2 = _get_peak_centre_range(x_list, y_list, apex, False,
                                           peak_indices_list, centroid_fraction, step_size_fn)
            right2 = _get_peak_centre_range(x_list, y_list, apex, True,
                                            peak_indices_list, centroid_fraction, step_size_fn)
            left3, right3 = _expand_peak_range(x_list, y_list, left2, right2,
                                               peak_indices_list, step_size_fn)
            _set_peak_range(peak_indices_list, left3, right3, apex)

    peak_indices[:] = np.asarray(peak_indices_list, dtype=np.int64)


# --------------------------------------------------------------------------- #
# Centroid extraction (numpy vectorised, no numba)
# --------------------------------------------------------------------------- #

def _extract_centroids_numpy(
    mz_arr: np.ndarray,
    int_arr: np.ndarray,
    x_smooth: np.ndarray,
    y_smooth: np.ndarray,
    peak_indices: np.ndarray,
    centroid_fraction: float,
    allow_non_zero_baseline: bool,
) -> list[Peak]:
    """Vectorised centroid extraction from assigned peak ranges."""
    peaks: list[Peak] = []
    num_points = len(peak_indices)
    if num_points == 0:
        return peaks

    valid = peak_indices != -1
    if not np.any(valid):
        return peaks

    is_start = np.zeros(num_points, dtype=bool)
    is_start[0] = valid[0]
    is_start[1:] = valid[1:] & ((peak_indices[1:] != peak_indices[:-1]) | (~valid[:-1]))

    is_end = np.zeros(num_points, dtype=bool)
    is_end[-1] = valid[-1]
    is_end[:-1] = valid[:-1] & ((peak_indices[:-1] != peak_indices[1:]) | (~valid[1:]))

    region_starts = np.where(is_start)[0]
    region_ends = np.where(is_end)[0]

    for start, end in zip(region_starts, region_ends):
        apex_idx = int(peak_indices[start])
        if apex_idx < 0:
            continue

        orig_start = int(np.searchsorted(mz_arr, x_smooth[start], side="right"))
        if orig_start > 0 and mz_arr[orig_start - 1] == x_smooth[start]:
            orig_start -= 1

        orig_end = int(np.searchsorted(mz_arr, x_smooth[end], side="left"))
        if orig_end >= len(mz_arr):
            orig_end = len(mz_arr) - 1

        if orig_end <= orig_start:
            continue

        threshold = float(y_smooth[apex_idx]) * centroid_fraction

        if allow_non_zero_baseline:
            baseline = min(
                float(int_arr[orig_start]),
                float(int_arr[orig_end]),
                float(y_smooth[start]),
                float(y_smooth[end]),
            )
        else:
            baseline = 0.0

        mz_slice = mz_arr[orig_start:orig_end + 1]
        int_slice = int_arr[orig_start:orig_end + 1]
        contrib = int_slice - baseline
        mask = contrib > threshold

        if not np.any(mask):
            continue

        weights = contrib[mask]
        weight_sum = float(weights.sum())
        if weight_sum <= 0:
            continue

        centroid_mz_val = float(np.dot(mz_slice[mask], weights)) / weight_sum
        area = weight_sum

        peak_height = float(int_slice.max())
        apex_local = int(int_slice.argmax())
        peak_apex_x = float(mz_slice[apex_local])
        peak_apex_y = float(int_slice[apex_local])

        peaks.append(Peak(
            x_value=centroid_mz_val,
            area=area,
            height=peak_height,
            apex_x=peak_apex_x,
            apex_y=peak_apex_y,
            start_x=float(mz_arr[orig_start]),
            start_y=float(int_arr[orig_start]),
            end_x=float(mz_arr[orig_end]),
            end_y=float(int_arr[orig_end]),
        ))

    return peaks


# --------------------------------------------------------------------------- #
# Centroiding fallback
# --------------------------------------------------------------------------- #

def _centroid_spectrum_fallback(mz, intensities, centroid_percentage, allow_non_zero_baseline, return_arrays=False):
    """Pure-Python fallback when numba is absent."""
    mz_arr = np.asarray(mz, dtype=np.float64)
    int_arr = np.asarray(intensities, dtype=np.float64)
    centroid_fraction = centroid_percentage / 100.0

    def step_size_fn(x_val):
        idx = bisect_left(mz_arr, x_val)
        n = len(mz_arr)
        if 0 < idx < n:
            return float(mz_arr[idx] - mz_arr[idx - 1])
        return float(mz_arr[1] - mz_arr[0]) if n > 1 else 1.0

    x_smooth_list, y_smooth_list = add_framing_zeros(
        mz_arr.tolist(), int_arr.tolist(), step_size_fn,
        half_insert=1, insert_before_first=True, insert_after_last=True,
    )
    x_smooth = np.asarray(x_smooth_list, dtype=np.float64)
    y_smooth_raw = np.asarray(y_smooth_list, dtype=np.float64)

    limit = min(100, len(mz_arr))
    if limit > 1:
        avg_step = float(np.mean(np.diff(mz_arr[:limit])))
    else:
        avg_step = 1.0
    half_window = 2 if (avg_step < 0.01 and centroid_fraction <= 0.5) else 1
    y_smooth = moving_average_smooth(y_smooth_raw, half_window)

    num_points = len(x_smooth)
    if num_points < 3:
        if return_arrays:
            return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float64)
        return [], []

    is_max = (y_smooth[1:-1] > y_smooth[:-2]) & (y_smooth[1:-1] >= y_smooth[2:])
    local_max_indices = np.where(is_max)[0] + 1
    local_max_values = y_smooth[local_max_indices]

    if local_max_indices.size == 0:
        if return_arrays:
            return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float64)
        return [], []

    sorted_order = np.argsort(local_max_values, kind="mergesort")[::-1]

    peak_indices = np.full(num_points, -1, dtype=np.int64)
    _assign_peak_ranges_py(
        x_smooth, y_smooth, local_max_indices, sorted_order,
        peak_indices, centroid_fraction, step_size_fn,
    )

    peaks = _extract_centroids_numpy(
        mz_arr, int_arr, x_smooth, y_smooth, peak_indices,
        centroid_fraction, allow_non_zero_baseline,
    )

    if not peaks:
        if return_arrays:
            return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float64)
        return [], []

    cmz = np.array([p.x_value for p in peaks], dtype=np.float64)
    cint = np.array([p.area for p in peaks], dtype=np.float64)
    order = np.argsort(cmz)

    if return_arrays:
        return cmz[order], cint[order]
    return cmz[order].tolist(), cint[order].tolist()


def centroid_spectrum(
    mz,
    intensities,
    centroid_percentage: float = 50.0,
    allow_non_zero_baseline: bool = False,
    return_arrays: bool = False,
):
    """Convert a profile spectrum to centroid format (pure-Python fallback)."""
    n = len(mz)
    if n < 3:
        mz_out = [m for m, v in zip(mz, intensities) if v > 0]
        int_out = [v for v in intensities if v > 0]
        if return_arrays:
            return np.asarray(mz_out, dtype=np.float64), np.asarray(int_out, dtype=np.float64)
        return mz_out, int_out

    return _centroid_spectrum_fallback(
        mz, intensities, centroid_percentage, allow_non_zero_baseline, return_arrays=return_arrays
    )
