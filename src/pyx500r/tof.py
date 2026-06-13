"""Pure-Python SCIEX TOF spectrum codec.

Reverse-engineered, byte-exact reimplementation of
``Clearcore2.Compression.DecompressionAlgorithmTof`` and
``Sciex.FMan.DefaultTofCalibration`` with **no .NET / DLL dependency**.

A WIFF2 acquisition stores each TOF spectrum as a compressed run-length stream
of (time-bin, intensity) pairs inside the companion ``.wiff.scan`` file. The
m/z axis is recovered from a per-scan quadratic calibration.

Performance-critical paths are JIT-compiled with **numba** (when available) for
near-C speed. A pure-Python fallback is always available.
"""

from __future__ import annotations

import struct
from math import sqrt

try:
    from numba import njit
    _HAS_NUMBA = True
except ImportError:
    _HAS_NUMBA = False

    # no-op decorator when numba is absent
    def njit(*args, **kwargs):  # noqa: F811
        def _identity(fn):
            return fn
        return _identity


# DecompressionAlgorithmTof control constants (decompiled).
_IS_ZEROS_MASK = 0x80
_ONE_BYTE_VAL = 124
_TWO_BYTE_VAL = 125
_FOUR_BYTE_VAL = 126
_MAX_ONE_BYTE_VAL = 123
_STOP_MARKER = 0xFF
_FIXED_BIN_MARKER = b"\xff\xff\xff\xff"


# --------------------------------------------------------------------------- #
# Core decode loop (pure Python — always available)
# --------------------------------------------------------------------------- #
def _decode_loop(
    stream: bytes,
    pos_start: int,
    start_bin: int,
    step: int,
    min_bin: int,
) -> tuple[list[int], list[int]]:
    """Pure-Python RLE decode loop. Returns (bins, intensities)."""
    pos = pos_start
    time_bin = start_bin
    bins: list[int] = []
    ints: list[int] = []

    while pos < len(stream):
        byte = stream[pos]
        pos += 1
        if byte == _STOP_MARKER:
            break
        token = byte & 0x7F
        if token <= _MAX_ONE_BYTE_VAL:
            value = token
        elif token == _ONE_BYTE_VAL:
            value = stream[pos]
            pos += 1
        elif token == _TWO_BYTE_VAL:
            value = stream[pos] | (stream[pos + 1] << 8)
            pos += 2
        elif token == _FOUR_BYTE_VAL:
            value = stream[pos] | (stream[pos + 1] << 8) | (stream[pos + 2] << 16) | (stream[pos + 3] << 24)
            pos += 4
        else:
            value = 0

        if byte & _IS_ZEROS_MASK:
            time_bin += value * step
        else:
            if value != 0 and time_bin >= min_bin:
                bins.append(time_bin)
                ints.append(value)
            time_bin += step

    return bins, ints


# --------------------------------------------------------------------------- #
# Numba JIT kernels (conditionally compiled)
# --------------------------------------------------------------------------- #
if _HAS_NUMBA:
    import numpy as np  # numba always pulls in numpy

    @njit(cache=True)
    def _decompress_tof_kernel_fast(
        stream, pos_start, start_bin, step, min_bin,
        bins_out, ints_out, max_out,
    ) -> int:
        """Core JIT loop — returns the number of written points."""
        pos = pos_start
        time_bin = start_bin
        count = 0
        while pos < len(stream):
            byte = stream[pos]
            pos += 1
            if byte == _STOP_MARKER:
                break
            token = byte & 0x7F
            if token <= _MAX_ONE_BYTE_VAL:
                value = token
            elif token == _ONE_BYTE_VAL:
                value = stream[pos]; pos += 1
            elif token == _TWO_BYTE_VAL:
                value = stream[pos] | (stream[pos + 1] << 8); pos += 2
            elif token == _FOUR_BYTE_VAL:
                value = stream[pos] | (stream[pos + 1] << 8) | (stream[pos + 2] << 16) | (stream[pos + 3] << 24); pos += 4
            else:
                value = 0
            if byte & _IS_ZEROS_MASK:
                time_bin += value * step
            else:
                if value != 0 and time_bin >= min_bin:
                    if count < max_out:
                        bins_out[count] = time_bin
                        ints_out[count] = value
                        count += 1
                time_bin += step
        return count

    @njit(cache=True)
    def _decompress_tof_kernel_calibrated(
        stream, cal_a, cal_t0, time_resolution,
        pos_start, start_bin, step, min_bin,
        mz_out, ints_out, max_out,
    ) -> int:
        """JIT loop with m/z calibration fused into the kernel."""
        pos = pos_start
        time_bin = start_bin
        count = 0
        while pos < len(stream):
            byte = stream[pos]
            pos += 1
            if byte == _STOP_MARKER:
                break
            token = byte & 0x7F
            if token <= _MAX_ONE_BYTE_VAL:
                value = token
            elif token == _ONE_BYTE_VAL:
                value = stream[pos]; pos += 1
            elif token == _TWO_BYTE_VAL:
                value = stream[pos] | (stream[pos + 1] << 8); pos += 2
            elif token == _FOUR_BYTE_VAL:
                value = stream[pos] | (stream[pos + 1] << 8) | (stream[pos + 2] << 16) | (stream[pos + 3] << 24); pos += 4
            else:
                value = 0
            if byte & _IS_ZEROS_MASK:
                time_bin += value * step
            else:
                if value != 0 and time_bin >= min_bin:
                    if count < max_out:
                        mz_out[count] = (cal_a * time_resolution * time_bin - cal_a * cal_t0) ** 2
                        ints_out[count] = value
                        count += 1
                time_bin += step
        return count

    def _decompress_tof_numba(stream_bytes, number_of_time_bins_to_sum, min_bin, return_arrays=False):
        """Numba-accelerated decompression (returns raw time bins)."""
        stream = np.frombuffer(stream_bytes, dtype=np.uint8)
        length = len(stream)
        has_marker = length >= 8 and stream[0] == 0xFF and stream[1] == 0xFF and stream[2] == 0xFF and stream[3] == 0xFF
        if has_marker:
            start_bin = struct.unpack_from("<I", stream_bytes, 4)[0]
            pos_start = 8
            step = number_of_time_bins_to_sum
        else:
            start_bin = 1; pos_start = 0; step = 1
        max_out = length
        bins_arr = np.empty(max_out, dtype=np.int64)
        ints_arr = np.empty(max_out, dtype=np.int64)
        count = _decompress_tof_kernel_fast(stream, pos_start, start_bin, step, min_bin, bins_arr, ints_arr, max_out)
        if return_arrays:
            return bins_arr[:count], ints_arr[:count]
        return bins_arr[:count].tolist(), ints_arr[:count].tolist()

    def _decompress_tof_numba_calibrated(
        stream_bytes, cal_a, cal_t0, time_resolution,
        number_of_time_bins_to_sum, min_bin, return_arrays=False,
    ):
        """Numba-accelerated decompression with fused m/z calibration."""
        stream = np.frombuffer(stream_bytes, dtype=np.uint8)
        length = len(stream)
        has_marker = length >= 8 and stream[0] == 0xFF and stream[1] == 0xFF and stream[2] == 0xFF and stream[3] == 0xFF
        if has_marker:
            start_bin = struct.unpack_from("<I", stream_bytes, 4)[0]
            pos_start = 8
            step = number_of_time_bins_to_sum
        else:
            start_bin = 1; pos_start = 0; step = 1
        max_out = length
        mz_arr = np.empty(max_out, dtype=np.float64)
        ints_arr = np.empty(max_out, dtype=np.int64)
        count = _decompress_tof_kernel_calibrated(
            stream, cal_a, cal_t0, time_resolution,
            pos_start, start_bin, step, min_bin,
            mz_arr, ints_arr, max_out,
        )
        if return_arrays:
            return mz_arr[:count], ints_arr[:count]
        return mz_arr[:count].tolist(), ints_arr[:count].tolist()

    # JIT warm-up at import time
    def _warmup():
        e = np.empty(0, dtype=np.int64)
        f = np.empty(0, dtype=np.float64)
        _decompress_tof_kernel_fast(np.array([], dtype=np.uint8), 0, 0, 1, 0, e, e, 0)
        _decompress_tof_kernel_calibrated(np.array([], dtype=np.uint8), 1.0, 0.0, 1.0, 0, 0, 1, 0, f, e, 0)
    _warmup()
    del _warmup

else:
    def _decompress_tof_numba(stream_bytes, number_of_time_bins_to_sum, min_bin):
        raise NotImplementedError("numba not installed — use the pure-Python fallback")
    def _decompress_tof_numba_calibrated(stream_bytes, cal_a, cal_t0, time_resolution, number_of_time_bins_to_sum, min_bin):
        raise NotImplementedError("numba not installed — use the pure-Python fallback")


# --------------------------------------------------------------------------- #
# Public API
# --------------------------------------------------------------------------- #
def decompress_tof(
    stream: bytes,
    number_of_time_bins_to_sum: int = 1,
    min_bin: int = 0,
    cal_a: float | None = None,
    cal_t0: float | None = None,
    time_resolution: float | None = None,
    return_arrays: bool = False,
) -> tuple:
    """Decode a compressed TOF stream into ``(mz_or_bins, intensities)``.

    ``stream`` must begin at the packed-data sentinel. When the stream is
    prefixed with the ``FF FF FF FF`` fixed-bin marker, the starting time bin
    is read from the following ``uint32`` and bins advance in steps of
    ``number_of_time_bins_to_sum``; otherwise bins start at 1 and advance by 1.

    When **cal_a**, **cal_t0**, and **time_resolution** are provided, the m/z
    axis is computed and returned in place of raw time bins.

    If **numba** is installed, JIT-compiled kernels give near-C performance.
    Otherwise, a pure-Python loop is used (5-10x slower but always available).

    Parameters
    ----------
    return_arrays: if True, return numpy arrays instead of Python lists.
        This avoids a list→numpy round-trip when the output is fed directly
        into centroiding or other array-based consumers.
    """
    has_marker = len(stream) >= 8 and stream[0:4] == _FIXED_BIN_MARKER
    if has_marker:
        start_bin = struct.unpack_from("<I", stream, 4)[0]
        pos_start = 8
        step = number_of_time_bins_to_sum
    else:
        start_bin = 1
        pos_start = 0
        step = 1

    # --- fast path: numba JIT ---
    if _HAS_NUMBA:
        if cal_a is not None and cal_t0 is not None and time_resolution is not None:
            return _decompress_tof_numba_calibrated(
                stream, cal_a, cal_t0, time_resolution,
                number_of_time_bins_to_sum, min_bin, return_arrays,
            )
        return _decompress_tof_numba(stream, number_of_time_bins_to_sum, min_bin, return_arrays)

    # --- pure-Python fallback ---
    bins, ints = _decode_loop(stream, pos_start, start_bin, step, min_bin)
    if cal_a is not None and cal_t0 is not None and time_resolution is not None:
        cal = TofCalibration(cal_a=cal_a, cal_t0=cal_t0, time_resolution=time_resolution)
        mz = cal.bins_to_masses(bins)
        if return_arrays:
            import numpy as np
            return np.array(mz, dtype=np.float64), np.array(ints, dtype=np.int64)
        return mz, ints
    if return_arrays:
        import numpy as np
        return np.array(bins, dtype=np.int64), np.array(ints, dtype=np.int64)
    return bins, ints


class TofCalibration:
    """Quadratic TOF m/z calibration (``Sciex.FMan.DefaultTofCalibration``).

    ``m/z = (cal_a * time_resolution * bin - cal_a * cal_t0) ** 2``
    """

    __slots__ = ("cal_a", "cal_t0", "time_resolution", "_a_res", "_a_t0")

    def __init__(self, cal_a: float, cal_t0: float, time_resolution: float):
        if cal_a <= 0.0:
            raise ValueError("cal_a (slope) must be positive")
        if time_resolution <= 0.0:
            raise ValueError("time_resolution must be positive")
        self.cal_a = cal_a
        self.cal_t0 = cal_t0
        self.time_resolution = time_resolution
        self._a_res = cal_a * time_resolution
        self._a_t0 = cal_a * cal_t0

    def bin_to_mass(self, time_bin: float) -> float:
        root = self._a_res * time_bin - self._a_t0
        return root * root

    def bins_to_masses(self, time_bins: list[int]) -> list[float]:
        a_res = self._a_res
        a_t0 = self._a_t0
        return [(a_res * b - a_t0) ** 2 for b in time_bins]

    def mass_to_bin(self, mass: float) -> float:
        if mass <= 0.0:
            return 0.0
        return (sqrt(mass) / self.cal_a + self.cal_t0) / self.time_resolution


# --------------------------------------------------------------------------- #
# TOF compression (CompressionAlgorithmTof)
# --------------------------------------------------------------------------- #
# Compression uses the same token scheme as decompression but with separate
# tokens for zero-run lengths (high bit set) vs intensity values (high bit
# clear).  Zero-run tokens: 128+value (<=123), 252+1byte, 253+2byte, 254+4byte.
# Intensity tokens: same as decompression (<=123, 124+1byte, 125+2byte, 126+4byte).

_COMPRESS_ZERO_ONE_BYTE = 252   # 0xFC — 1-byte zero-run follows
_COMPRESS_ZERO_TWO_BYTE = 253   # 0xFD — 2-byte zero-run follows
_COMPRESS_ZERO_FOUR_BYTE = 254  # 0xFE — 4-byte zero-run follows


def compress_tof(
    bins: list[int],
    intensities: list[int],
    number_of_time_bins_to_sum: int = 1,
    threshold: int = 0,
) -> bytes:
    """Compress (bin, intensity) pairs into the TOF RLE format.

    This is the inverse of ``decompress_tof``.  The output begins with the
    fixed-bin marker (``FF FF FF FF``) followed by the starting bin as a
    little-endian uint32, then the RLE stream padded to a 4-byte boundary
    with ``0xFF`` stop markers.

    Parameters
    ----------
    bins: time-bin indices (must be strictly increasing).
    intensities: intensity values aligned with *bins*.
    number_of_time_bins_to_sum: bin step size (1 for no binning).
    threshold: intensities <= *threshold* are silently dropped.

    Returns
    -------
    bytes
        The compressed TOF record suitable for writing into a ``.wiff.scan``
        file (prefixed with the ``0x4C`` data-offset header if needed).
    """
    if not bins:
        return b""

    # Filter by threshold
    pairs = [(b, v) for b, v in zip(bins, intensities) if v > threshold]
    if not pairs:
        return b""

    out: list[int] = []

    # Fixed-bin marker (-1 as little-endian uint32 = FF FF FF FF)
    out.extend([0xFF, 0xFF, 0xFF, 0xFF])

    # Starting bin: first bin minus one step (the DLL stores bin - binSize)
    start_bin = pairs[0][0] - number_of_time_bins_to_sum
    out.extend(struct.pack("<I", start_bin & 0xFFFFFFFF))

    # Leading zero marker (the DLL always adds a 0x00 byte = "zero zeros")
    out.append(0x00)

    prev_bin = pairs[0][0]

    for bin_num, intensity in pairs:
        # Emit zero-run if there's a gap
        if bin_num != prev_bin:
            gaps = (bin_num - prev_bin) // number_of_time_bins_to_sum - 1
            if gaps > 0:
                _emit_value(out, gaps, is_zero_run=True)

        prev_bin = bin_num
        _emit_value(out, intensity, is_zero_run=False)

    # Pad to 4-byte boundary with 0xFF stop markers
    padding = (4 - (len(out) % 4)) % 4
    out.extend([0xFF] * padding)

    return bytes(out)


def _emit_value(out: list[int], value: int, is_zero_run: bool) -> None:
    """Append an RLE token to the output byte list."""
    if is_zero_run:
        # Zero-run: high bit set on the token
        if value <= 123:
            out.append(value | 0x80)
        elif value <= 255:
            out.append(_COMPRESS_ZERO_ONE_BYTE)
            out.append(value & 0xFF)
        elif value <= 65535:
            out.append(_COMPRESS_ZERO_TWO_BYTE)
            out.extend(struct.pack("<H", value))
        else:
            out.append(_COMPRESS_ZERO_FOUR_BYTE)
            out.extend(struct.pack("<I", value))
    else:
        # Intensity value: high bit clear
        if value <= 123:
            out.append(value)
        elif value <= 255:
            out.append(_ONE_BYTE_VAL)
            out.append(value & 0xFF)
        elif value <= 65535:
            out.append(_TWO_BYTE_VAL)
            out.extend(struct.pack("<H", value))
        else:
            out.append(_FOUR_BYTE_VAL)
            out.extend(struct.pack("<I", value))


# --------------------------------------------------------------------------- #
# Quadrupole decompression (DecompressionAlgorithmQuad)
# --------------------------------------------------------------------------- #
# Quad data has a different format: mass ranges with start/stop/step/scale,
# followed by RLE-encoded intensity values per mass bin.


class MassRange:
    """Quadrupole mass range descriptor.

    Each range has a start mass, stop mass, step mass (all in millim/z),
    and a scale factor to convert raw counts to intensity.
    """

    __slots__ = ("start_mass", "stop_mass", "step_mass", "scale_factor")

    def __init__(
        self,
        start_mass: float,
        stop_mass: float,
        step_mass: float,
        scale_factor: float,
    ):
        self.start_mass = int(start_mass * 1000.0)
        self.stop_mass = int(stop_mass * 1000.0)
        self.step_mass = int(step_mass * 1000.0)
        self.scale_factor = scale_factor


def decompress_quad(stream: bytes) -> tuple[list[float], list[float]]:
    """Decompress a quadrupole spectrum from its compressed byte stream.

    Returns ``(mz, intensities)`` with zero-intensity points removed.

    The format is:
    1. uint32: number of mass ranges
    2. For each range: 4 x double (start, stop, step, scale)
    3. uint32: total number of data points
    4. RLE-encoded intensity values with bin-skip control bytes
    """
    pos = 0

    def read_int() -> int:
        nonlocal pos
        val = struct.unpack_from("<I", stream, pos)[0]
        pos += 4
        return val

    def read_double() -> float:
        nonlocal pos
        val = struct.unpack_from("<d", stream, pos)[0]
        pos += 8
        return val

    # Read mass ranges
    num_ranges = read_int()
    ranges: list[MassRange] = []
    for _ in range(num_ranges):
        ranges.append(MassRange(read_double(), read_double(), read_double(), read_double()))

    if not ranges:
        return [], []

    # Read total point count
    total_points = read_int()

    # Decompress data points
    mz_out: list[float] = []
    int_out: list[float] = []

    current_range_idx = 0
    bin_number = 0
    points_written = 0
    flag = True  # alternates between reading header byte and data byte

    while pos < len(stream) and points_written < total_points:
        if flag:
            # Header byte: bits 7-5 = count (0-7), bit 7 indicates sign
            header = stream[pos] & 0xFF
            count = (header & 0x70) >> 4
            is_data = ((header & 0x80) == 0)  # bit 7 clear = data point
            pos += 1

            if count == 0:
                # Need to read more — this is a continuation
                continue
            flag = False
        else:
            # Data byte
            control = stream[pos] & 0xFF
            count = control & 0x07
            is_data = ((control & 0x08) == 0)
            pos += 1
            flag = True

        if pos >= len(stream) or points_written >= total_points:
            break

        # Read value
        if count == 0:
            value = stream[pos] & 0xFF
            pos += 1
            if value == 0:
                value = read_int()
        else:
            value = count

        if pos >= len(stream) or points_written >= total_points:
            break

        if is_data:
            step_mass = ranges[0].step_mass
            start_mass = ranges[0].start_mass
            mass_mz = (bin_number * step_mass + start_mass) / 1000.0

            if value == 0:
                # Zero intensity — skip
                pass
            else:
                if current_range_idx < len(ranges) - 1 and (bin_number * step_mass + start_mass) > ranges[current_range_idx + 1].start_mass:
                    current_range_idx += 1
                scale = ranges[current_range_idx].scale_factor
                intensity = float(value) * scale
                mz_out.append(mass_mz)
                int_out.append(intensity)
                points_written += 1

            bin_number += 1
        else:
            # Skip bins
            bin_number += count

    return mz_out, int_out


# --------------------------------------------------------------------------- #
# Zero-width (SRM/MRM) decompression (DecompressionAlgorithmZeroWidth)
# --------------------------------------------------------------------------- #
# Zero-width data is stored as an array of IEEE 754 floats.
# Positive float = intensity at current transition index.
# Negative float = skip N transitions (N = abs(int(value))).
# The first float is a header and is skipped.


def decompress_zero_width(
    stream: bytes,
    number_of_transitions: int,
    include_zero_intensity: bool = False,
) -> tuple[list[float], list[float]]:
    """Decompress a zero-width (SRM/MRM) spectrum.

    Parameters
    ----------
    stream: raw bytes (must be a multiple of 4).
    number_of_transitions: expected number of transition channels.
    include_zero_intensity: if True, return all transitions (zeros included).

    Returns
    ``(transition_index, intensities)``.
    """
    if not stream or len(stream) < 4:
        return ([i for i in range(number_of_transitions)] if include_zero_intensity else [],
                [0.0] * number_of_transitions if include_zero_intensity else [])

    import numpy as np

    vals = np.frombuffer(stream, dtype=np.float32)[1:]  # skip header float

    if include_zero_intensity:
        y_values = [0.0] * number_of_transitions
        x_values = [float(i) for i in range(number_of_transitions)]
        idx = 0
        for val in vals:
            if val < 0.0:
                idx += int(-val)
            else:
                if idx < number_of_transitions:
                    y_values[idx] = float(val)
                    idx += 1
        return x_values, y_values
    else:
        x_values: list[float] = []
        y_values: list[float] = []
        idx = 0
        for val in vals:
            if val < 0.0:
                idx += int(-val)
            else:
                x_values.append(float(idx))
                y_values.append(float(val))
                idx += 1
        return x_values, y_values


# --------------------------------------------------------------------------- #
# Mass/time conversion utilities (CompressionAlgorithmTof static methods)
# --------------------------------------------------------------------------- #


def mass_to_time(
    mass: float,
    slope: float,
    delay: float,
    resolution: float,
) -> int:
    """Convert m/z to TOF time bin.

    ``bin = round((sqrt(mass) / slope + delay) / resolution)``

    This is the inverse of :func:`time_to_mass`.
    """
    return int((sqrt(mass) / slope + delay) / resolution + 0.5)


def time_to_mass(
    time_bin: int,
    slope: float,
    delay: float,
    resolution: float,
) -> float:
    """Convert TOF time bin to m/z.

    ``m/z = (slope * (bin * resolution - delay)) ** 2``

    This is the inverse of :func:`mass_to_time`.
    """
    val = slope * (time_bin * resolution - delay)
    return val * val
