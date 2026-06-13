"""Trial module: monolithic numba kernel with inline='never'.

If this produces correct results, the approach can be merged into centroid.py.
"""

from __future__ import annotations

import numpy as np

try:
    from numba import njit
    _HAS_NUMBA = True
except ImportError:
    _HAS_NUMBA = False

    def njit(*args, **kwargs):  # noqa: F811
        def _identity(fn):
            return fn
        return _identity


# --------------------------------------------------------------------------- #
# Helper kernels  — inline='never' to prevent inlining into the monolith
# --------------------------------------------------------------------------- #

@njit(cache=True, inline='never')
def _bisect_left_njit(arr, x):
    lo, hi = 0, len(arr)
    while lo < hi:
        mid = (lo + hi) // 2
        if arr[mid] < x:
            lo = mid + 1
        else:
            hi = mid
    return lo


@njit(cache=True, inline='never')
def _step_size_njit(mz, x):
    idx = _bisect_left_njit(mz, x)
    n = len(mz)
    if 0 < idx < n:
        return mz[idx] - mz[idx - 1]
    return mz[1] - mz[0] if n > 1 else 1.0


@njit(cache=True, inline='never')
def _precompute_x_steps_linear_njit(mz, x_smooth, x_steps):
    """Two-pointer linear scan — O(n_smooth + n_mz) instead of O(n_smooth * log n_mz)."""
    n_mz = len(mz)
    n_smooth = len(x_smooth)
    mz_idx = 0
    for i in range(n_smooth):
        x = x_smooth[i]
        while mz_idx < n_mz and mz[mz_idx] < x:
            mz_idx += 1
        if 0 < mz_idx < n_mz:
            x_steps[i] = mz[mz_idx] - mz[mz_idx - 1]
        elif n_mz > 1:
            x_steps[i] = mz[1] - mz[0]
        else:
            x_steps[i] = 1.0


@njit(cache=True, inline='never')
def _add_zeros_njit(mz, intensities, insert_before_first, insert_after_last, x_out, y_out):
    n = len(mz)
    if n == 0:
        return 0
    pos = 0

    if insert_before_first:
        step = _step_size_njit(mz, mz[0])
        x_out[pos] = mz[0] - step
        y_out[pos] = 0.0
        pos += 1

    prev_x = mz[0]
    num_steps = 1
    next_insert_at = mz[0] - 1.0

    for i in range(n - 1):
        cur_x = mz[i]
        if cur_x >= next_insert_at:
            step = _step_size_njit(mz, cur_x)
            next_insert_at = cur_x + num_steps * step

        step_cur = _step_size_njit(mz, cur_x)
        step_prev = _step_size_njit(mz, prev_x)
        if cur_x - step_cur > prev_x + step_prev / 2.0:
            insert_x = cur_x - step_cur
            x_out[pos] = insert_x
            y_out[pos] = 0.0
            pos += 1

        x_out[pos] = cur_x
        y_out[pos] = intensities[i]
        pos += 1
        prev_x = cur_x

        cur_x_next = mz[i + 1]
        step_next = _step_size_njit(mz, cur_x_next)
        step_prev = _step_size_njit(mz, prev_x)
        if cur_x_next - step_next > prev_x + step_prev / 2.0:
            insert_x = prev_x + step_prev
            x_out[pos] = insert_x
            y_out[pos] = 0.0
            pos += 1
            prev_x = insert_x

    last_x = mz[n - 1]
    step_last = _step_size_njit(mz, last_x)
    step_prev = _step_size_njit(mz, prev_x)
    if last_x - step_last > prev_x + step_prev / 2.0:
        x_out[pos] = last_x - step_last
        y_out[pos] = 0.0
        pos += 1

    x_out[pos] = last_x
    y_out[pos] = intensities[n - 1]
    pos += 1

    if insert_after_last:
        step = _step_size_njit(mz, last_x)
        x_out[pos] = last_x + step
        y_out[pos] = 0.0
        pos += 1

    return pos


@njit(cache=True, inline='never')
def _moving_average_smooth_njit(y, half_window, out):
    n = len(y)
    window = 2 * half_window + 1

    for j in range(half_window + 1):
        s = 0.0
        for k in range(j - half_window, j + half_window + 1):
            idx = 0 if k < 0 else k
            s += y[idx]
        out[j] = s / window

    for j in range(half_window + 1, n - half_window):
        s = 0.0
        for k in range(j - half_window, j + half_window + 1):
            s += y[k]
        out[j] = s / window

    for j in range(n - half_window, n):
        s = 0.0
        for k in range(j - half_window, j + half_window + 1):
            idx = n - 1 if k >= n else k
            s += y[idx]
        out[j] = s / window


@njit(cache=True, inline='never')
def _find_local_maxima_njit(y, out_indices, out_values):
    n = len(y)
    count = 0
    for i in range(1, n - 1):
        if y[i] > y[i - 1] and y[i] >= y[i + 1]:
            out_indices[count] = i
            out_values[count] = y[i]
            count += 1
    return count


@njit(cache=True, inline='never')
def _estimate_step_njit(x, index):
    n = len(x)
    if n < 2:
        return 1.0
    if index > 0 and index < n - 1:
        return (x[index + 1] - x[index - 1]) / 2.0
    return x[1] - x[0]


@njit(cache=True, inline='never')
def _guess_rough_peak_width_njit(step):
    if step < 0.031:
        return 0.2
    if step < 0.061:
        return 0.3
    if step < 0.21:
        return 0.8
    if step < 0.51:
        return 1.0
    return 4.0 * step


@njit(cache=True, inline='never')
def _find_point_at_height_njit(x, y, index, height_fraction, look_right, peak_indices):
    step = _estimate_step_njit(x, index)
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


@njit(cache=True, inline='never')
def _get_peak_centre_range_njit(x, y, index, look_right, peak_indices, centroid_fraction, step):
    if centroid_fraction < 0.4999:
        height_frac = 0.25
    elif centroid_fraction <= 0.75:
        height_frac = 0.5
    else:
        height_frac = 0.7

    rough_width = _guess_rough_peak_width_njit(step)
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


@njit(cache=True, inline='never')
def _expand_peak_range_njit(x, y, left, right, peak_indices, step_left, step_right):
    num_points = len(x)

    left -= 1
    while (left >= 0 and peak_indices[left] == -1
           and y[left] <= y[left + 1]
           and abs(x[left + 1] - x[left] - step_left) < step_left / 2.0):
        left -= 1
    left += 1
    while left < num_points and y[left] == 0.0 and peak_indices[left] == -1:
        left += 1

    right += 1
    while (right < num_points and peak_indices[right] == -1
           and y[right - 1] >= y[right]
           and abs(x[right] - x[right - 1] - step_right) < step_right / 2.0):
        right += 1
    right -= 1
    while right >= 0 and y[right] == 0.0 and peak_indices[right] == -1:
        right -= 1

    return left, right


@njit(cache=True, inline='never')
def _searchsorted_right_njit(arr, x):
    lo, hi = 0, len(arr)
    while lo < hi:
        mid = (lo + hi) // 2
        if arr[mid] <= x:
            lo = mid + 1
        else:
            hi = mid
    return lo


@njit(cache=True, inline='never')
def _searchsorted_left_njit(arr, x):
    lo, hi = 0, len(arr)
    while lo < hi:
        mid = (lo + hi) // 2
        if arr[mid] < x:
            lo = mid + 1
        else:
            hi = mid
    return lo


@njit(cache=True, inline='never')
def _extract_centroids_njit(mz_orig, int_orig, x_smooth, y_smooth, peak_indices, centroid_fraction, allow_non_zero_baseline):
    n = len(peak_indices)
    max_peaks = n
    cmz = np.empty(max_peaks, dtype=np.float64)
    cint = np.empty(max_peaks, dtype=np.float64)
    num_peaks = 0

    i = 0
    while i < n:
        if peak_indices[i] == -1:
            i += 1
            continue

        peak_id = peak_indices[i]
        start = i
        while i < n and peak_indices[i] == peak_id:
            i += 1
        end = i - 1

        apex_idx = int(peak_id)

        orig_start = _searchsorted_right_njit(mz_orig, x_smooth[start])
        if orig_start > 0 and mz_orig[orig_start - 1] == x_smooth[start]:
            orig_start -= 1

        orig_end = _searchsorted_left_njit(mz_orig, x_smooth[end])
        if orig_end >= len(mz_orig):
            orig_end = len(mz_orig) - 1

        if orig_end <= orig_start:
            continue

        threshold = y_smooth[apex_idx] * centroid_fraction

        if allow_non_zero_baseline:
            baseline = min(int_orig[orig_start], int_orig[orig_end], y_smooth[start], y_smooth[end])
        else:
            baseline = 0.0

        weight_sum = 0.0
        weighted_mz_sum = 0.0

        for k in range(orig_start, orig_end + 1):
            contrib = int_orig[k] - baseline
            if contrib > threshold:
                weight_sum += contrib
                weighted_mz_sum += mz_orig[k] * contrib

        if weight_sum <= 0.0:
            continue

        cmz[num_peaks] = weighted_mz_sum / weight_sum
        cint[num_peaks] = weight_sum
        num_peaks += 1

    return cmz[:num_peaks], cint[:num_peaks]


# --------------------------------------------------------------------------- #
# Monolithic kernel  — calls the helpers above
# --------------------------------------------------------------------------- #

@njit(cache=True)
def _centroid_spectrum_fast_njit(mz, intensities, centroid_fraction, allow_non_zero_baseline):
    n = len(mz)
    if n < 3:
        max_out = n
        cmz = np.empty(max_out, dtype=np.float64)
        cint = np.empty(max_out, dtype=np.float64)
        count = 0
        for i in range(n):
            if intensities[i] > 0:
                cmz[count] = mz[i]
                cint[count] = intensities[i]
                count += 1
        return cmz[:count], cint[:count]

    # 1. Add zeros
    max_smooth = 3 * n + 2
    x_smooth = np.empty(max_smooth, dtype=np.float64)
    y_smooth_raw = np.empty(max_smooth, dtype=np.float64)
    num_smooth = _add_zeros_njit(mz, intensities, True, True, x_smooth, y_smooth_raw)

    if num_smooth < 3:
        return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float64)

    x_smooth = x_smooth[:num_smooth]
    y_smooth_raw = y_smooth_raw[:num_smooth]

    # 2. Determine half_window
    limit = min(100, n)
    if limit > 1:
        total = 0.0
        for i in range(1, limit):
            total += mz[i] - mz[i - 1]
        avg_step = total / (limit - 1)
    else:
        avg_step = 1.0
    half_window = 2 if (avg_step < 0.01 and centroid_fraction <= 0.5) else 1

    # 3. Smooth
    y_smooth = np.empty(num_smooth, dtype=np.float64)
    _moving_average_smooth_njit(y_smooth_raw, half_window, y_smooth)

    # 4. Local maxima
    max_indices = np.empty(num_smooth, dtype=np.int64)
    max_values = np.empty(num_smooth, dtype=np.float64)
    num_maxima = _find_local_maxima_njit(y_smooth, max_indices, max_values)

    if num_maxima == 0:
        return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float64)

    # Sort by intensity descending
    sorted_order = np.argsort(max_values[:num_maxima], kind="mergesort")
    for i in range(num_maxima // 2):
        j = num_maxima - 1 - i
        sorted_order[i], sorted_order[j] = sorted_order[j], sorted_order[i]

    # Pre-compute step sizes for x_smooth
    x_steps = np.empty(num_smooth, dtype=np.float64)
    _precompute_x_steps_linear_njit(mz, x_smooth, x_steps)

    # 5. Assign peak ranges
    peak_indices = np.full(num_smooth, -1, dtype=np.int64)

    for order_idx in range(num_maxima):
        apex = int(max_indices[sorted_order[order_idx]])
        if peak_indices[apex] != -1:
            continue

        left = _find_point_at_height_njit(x_smooth, y_smooth, apex, 0.85, False, peak_indices)
        right = _find_point_at_height_njit(x_smooth, y_smooth, apex, 0.85, True, peak_indices)

        left_assigned = peak_indices[left]
        right_assigned = peak_indices[right]

        if left_assigned != -1 and right_assigned != -1:
            if (x_smooth[apex] - x_smooth[left]) < (x_smooth[right] - x_smooth[apex]):
                peak_indices[left:right] = left_assigned
            else:
                peak_indices[left + 1:right + 1] = right_assigned
        elif left_assigned != -1:
            peak_indices[left:right + 1] = left_assigned
        elif right_assigned != -1:
            peak_indices[left:right + 1] = right_assigned
        else:
            left2 = _get_peak_centre_range_njit(x_smooth, y_smooth, apex, False, peak_indices, centroid_fraction, x_steps[apex])
            right2 = _get_peak_centre_range_njit(x_smooth, y_smooth, apex, True, peak_indices, centroid_fraction, x_steps[apex])
            left3, right3 = _expand_peak_range_njit(x_smooth, y_smooth, left2, right2, peak_indices, x_steps[left2], x_steps[right2])
            peak_indices[left3:right3 + 1] = apex

    # 6. Extract centroids
    return _extract_centroids_njit(mz, intensities, x_smooth, y_smooth, peak_indices, centroid_fraction, allow_non_zero_baseline)


# --------------------------------------------------------------------------- #
# Warm-up
# --------------------------------------------------------------------------- #

def _warmup():
    x = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    y = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    _centroid_spectrum_fast_njit(x, y, 0.5, False)


if _HAS_NUMBA:
    _warmup()
    del _warmup


# --------------------------------------------------------------------------- #
# Public API
# --------------------------------------------------------------------------- #

def centroid_spectrum(mz, intensities, centroid_percentage=50.0, allow_non_zero_baseline=False, return_arrays=False):
    n = len(mz)
    if n < 3:
        mz_out = [m for m, v in zip(mz, intensities) if v > 0]
        int_out = [v for v in intensities if v > 0]
        if return_arrays:
            return np.asarray(mz_out, dtype=np.float64), np.asarray(int_out, dtype=np.float64)
        return mz_out, int_out

    if not _HAS_NUMBA:
        raise RuntimeError("numba is required for centroid_new")

    mz_arr = np.asarray(mz, dtype=np.float64)
    int_arr = np.asarray(intensities, dtype=np.float64)
    centroid_fraction = centroid_percentage / 100.0
    cmz, cint = _centroid_spectrum_fast_njit(mz_arr, int_arr, centroid_fraction, allow_non_zero_baseline)

    if len(cmz) == 0:
        if return_arrays:
            return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float64)
        return [], []

    order = np.argsort(cmz)
    if return_arrays:
        return cmz[order], cint[order]
    return cmz[order].tolist(), cint[order].tolist()


# --------------------------------------------------------------------------- #
# Debug helper — returns peak_indices for comparison
# --------------------------------------------------------------------------- #

@njit(cache=False)
def _debug_peak_indices_njit(mz, intensities, centroid_fraction):
    n = len(mz)
    max_smooth = 3 * n + 2
    x_smooth = np.empty(max_smooth, dtype=np.float64)
    y_smooth_raw = np.empty(max_smooth, dtype=np.float64)
    num_smooth = _add_zeros_njit(mz, intensities, True, True, x_smooth, y_smooth_raw)
    x_smooth = x_smooth[:num_smooth]
    y_smooth_raw = y_smooth_raw[:num_smooth]

    limit = min(100, n)
    if limit > 1:
        total = 0.0
        for i in range(1, limit):
            total += mz[i] - mz[i - 1]
        avg_step = total / (limit - 1)
    else:
        avg_step = 1.0
    half_window = 2 if (avg_step < 0.01 and centroid_fraction <= 0.5) else 1

    y_smooth = np.empty(num_smooth, dtype=np.float64)
    _moving_average_smooth_njit(y_smooth_raw, half_window, y_smooth)

    max_indices = np.empty(num_smooth, dtype=np.int64)
    max_values = np.empty(num_smooth, dtype=np.float64)
    num_maxima = _find_local_maxima_njit(y_smooth, max_indices, max_values)

    sorted_order = np.argsort(max_values[:num_maxima], kind="mergesort")
    for i in range(num_maxima // 2):
        j = num_maxima - 1 - i
        sorted_order[i], sorted_order[j] = sorted_order[j], sorted_order[i]

    x_steps = np.empty(num_smooth, dtype=np.float64)
    for i in range(num_smooth):
        x_steps[i] = _step_size_njit(mz, x_smooth[i])

    peak_indices = np.full(num_smooth, -1, dtype=np.int64)

    for order_idx in range(num_maxima):
        apex = int(max_indices[sorted_order[order_idx]])
        if peak_indices[apex] != -1:
            continue
        left = _find_point_at_height_njit(x_smooth, y_smooth, apex, 0.85, False, peak_indices)
        right = _find_point_at_height_njit(x_smooth, y_smooth, apex, 0.85, True, peak_indices)
        left_assigned = peak_indices[left]
        right_assigned = peak_indices[right]
        if left_assigned != -1 and right_assigned != -1:
            if (x_smooth[apex] - x_smooth[left]) < (x_smooth[right] - x_smooth[apex]):
                peak_indices[left:right] = left_assigned
            else:
                peak_indices[left + 1:right + 1] = right_assigned
        elif left_assigned != -1:
            peak_indices[left:right + 1] = left_assigned
        elif right_assigned != -1:
            peak_indices[left:right + 1] = right_assigned
        else:
            left2 = _get_peak_centre_range_njit(x_smooth, y_smooth, apex, False, peak_indices, centroid_fraction, x_steps[apex])
            right2 = _get_peak_centre_range_njit(x_smooth, y_smooth, apex, True, peak_indices, centroid_fraction, x_steps[apex])
            left3, right3 = _expand_peak_range_njit(x_smooth, y_smooth, left2, right2, peak_indices, x_steps[left2], x_steps[right2])
            peak_indices[left3:right3 + 1] = apex

    return peak_indices


# --------------------------------------------------------------------------- #
# Debug: return assignment log (apex, left, right, assigned_to)
# --------------------------------------------------------------------------- #

@njit(cache=False)
def _debug_assignment_log_njit(mz, intensities, centroid_fraction):
    n = len(mz)
    max_smooth = 3 * n + 2
    x_smooth_big = np.empty(max_smooth, dtype=np.float64)
    y_smooth_raw_big = np.empty(max_smooth, dtype=np.float64)
    num_smooth = _add_zeros_njit(mz, intensities, True, True, x_smooth_big, y_smooth_raw_big)

    x_smooth = np.empty(num_smooth, dtype=np.float64)
    y_smooth_raw = np.empty(num_smooth, dtype=np.float64)
    for i in range(num_smooth):
        x_smooth[i] = x_smooth_big[i]
        y_smooth_raw[i] = y_smooth_raw_big[i]

    limit = min(100, n)
    if limit > 1:
        total = 0.0
        for i in range(1, limit):
            total += mz[i] - mz[i - 1]
        avg_step = total / (limit - 1)
    else:
        avg_step = 1.0
    half_window = 2 if (avg_step < 0.01 and centroid_fraction <= 0.5) else 1

    y_smooth = np.empty(num_smooth, dtype=np.float64)
    _moving_average_smooth_njit(y_smooth_raw, half_window, y_smooth)

    max_indices = np.empty(num_smooth, dtype=np.int64)
    max_values = np.empty(num_smooth, dtype=np.float64)
    num_maxima = _find_local_maxima_njit(y_smooth, max_indices, max_values)

    sorted_order = np.argsort(max_values[:num_maxima], kind="mergesort")
    for i in range(num_maxima // 2):
        j = num_maxima - 1 - i
        sorted_order[i], sorted_order[j] = sorted_order[j], sorted_order[i]

    x_steps = np.empty(num_smooth, dtype=np.float64)
    for i in range(num_smooth):
        x_steps[i] = _step_size_njit(mz, x_smooth[i])

    peak_indices = np.full(num_smooth, -1, dtype=np.int64)
    log = np.empty((num_maxima, 5), dtype=np.int64)
    log_idx = 0

    for order_idx in range(num_maxima):
        apex = int(max_indices[sorted_order[order_idx]])
        if peak_indices[apex] != -1:
            continue
        left = _find_point_at_height_njit(x_smooth, y_smooth, apex, 0.85, False, peak_indices)
        right = _find_point_at_height_njit(x_smooth, y_smooth, apex, 0.85, True, peak_indices)
        left_assigned = peak_indices[left]
        right_assigned = peak_indices[right]
        action = 0
        if left_assigned != -1 and right_assigned != -1:
            action = 1
            if (x_smooth[apex] - x_smooth[left]) < (x_smooth[right] - x_smooth[apex]):
                peak_indices[left:right] = left_assigned
            else:
                peak_indices[left + 1:right + 1] = right_assigned
        elif left_assigned != -1:
            action = 2
            peak_indices[left:right + 1] = left_assigned
        elif right_assigned != -1:
            action = 3
            peak_indices[left:right + 1] = right_assigned
        else:
            action = 4
            left2 = _get_peak_centre_range_njit(x_smooth, y_smooth, apex, False, peak_indices, centroid_fraction, x_steps[apex])
            right2 = _get_peak_centre_range_njit(x_smooth, y_smooth, apex, True, peak_indices, centroid_fraction, x_steps[apex])
            left3, right3 = _expand_peak_range_njit(x_smooth, y_smooth, left2, right2, peak_indices, x_steps[left2], x_steps[right2])
            peak_indices[left3:right3 + 1] = apex
            left = left3
            right = right3
        log[log_idx, 0] = apex
        log[log_idx, 1] = left
        log[log_idx, 2] = right
        log[log_idx, 3] = action
        log[log_idx, 4] = sorted_order[order_idx]
        log_idx += 1

    return log[:log_idx]


# --------------------------------------------------------------------------- #
# Debug: return intermediate arrays for comparison
# --------------------------------------------------------------------------- #

@njit(cache=False)
def _debug_intermediates_njit(mz, intensities, centroid_fraction):
    n = len(mz)
    max_smooth = 3 * n + 2
    x_smooth_big = np.empty(max_smooth, dtype=np.float64)
    y_smooth_raw_big = np.empty(max_smooth, dtype=np.float64)
    num_smooth = _add_zeros_njit(mz, intensities, True, True, x_smooth_big, y_smooth_raw_big)

    x_smooth = np.empty(num_smooth, dtype=np.float64)
    y_smooth_raw = np.empty(num_smooth, dtype=np.float64)
    for i in range(num_smooth):
        x_smooth[i] = x_smooth_big[i]
        y_smooth_raw[i] = y_smooth_raw_big[i]

    limit = min(100, n)
    if limit > 1:
        total = 0.0
        for i in range(1, limit):
            total += mz[i] - mz[i - 1]
        avg_step = total / (limit - 1)
    else:
        avg_step = 1.0
    half_window = 2 if (avg_step < 0.01 and centroid_fraction <= 0.5) else 1

    y_smooth = np.empty(num_smooth, dtype=np.float64)
    _moving_average_smooth_njit(y_smooth_raw, half_window, y_smooth)

    max_indices = np.empty(num_smooth, dtype=np.int64)
    max_values = np.empty(num_smooth, dtype=np.float64)
    num_maxima = _find_local_maxima_njit(y_smooth, max_indices, max_values)

    sorted_order = np.argsort(max_values[:num_maxima], kind="mergesort")
    for i in range(num_maxima // 2):
        j = num_maxima - 1 - i
        sorted_order[i], sorted_order[j] = sorted_order[j], sorted_order[i]

    return x_smooth, y_smooth, max_indices[:num_maxima], max_values[:num_maxima], sorted_order
