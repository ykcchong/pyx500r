"""JSON serialization helpers for the pyx500r reference server.

Converts pyx500r's frozen dataclasses and the ``UnifiedPeak`` view into
JSON-native dicts, handling the two non-native cases:

* numpy arrays  -> Python lists (via ``.tolist()``)
* ``datetime``  -> ISO-8601 strings

See ../../docs/GUI_INTEGRATION.md §3 for the rationale.
"""

from __future__ import annotations

import dataclasses
from datetime import datetime
from typing import Any


def _to_list(seq: Any) -> list:
    """Coerce a Sequence/np.ndarray of numbers into a JSON-native list[float]."""
    if seq is None:
        return []
    tolist = getattr(seq, "tolist", None)
    if callable(tolist):
        return tolist()
    return [float(x) for x in seq]


def _jsonify(value: Any) -> Any:
    """Recursively coerce a value into something json.dumps can handle.

    Anything that is not a JSON-native scalar/list/dict (e.g. a stray
    BinaryFormatter reference placeholder or an enum) is stringified as a last
    resort, so an endpoint never 500s on an unexpected nested type.
    """
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, datetime):
        return value.isoformat()
    if hasattr(value, "tolist") and not isinstance(value, (str, bytes)):
        return value.tolist()
    if isinstance(value, dict):
        return {str(k): _jsonify(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_jsonify(v) for v in value]
    return str(value)


def dataclass_to_dict(obj: Any) -> dict[str, Any]:
    """Serialize any pyx500r frozen dataclass to a JSON-native dict.

    Handles ``datetime`` and numpy arrays anywhere in the structure.
    """
    raw = dataclasses.asdict(obj)
    return {k: _jsonify(v) for k, v in raw.items()}


def serialize_spectrum(spec: Any) -> dict[str, Any]:
    """Serialize a SpectrumData, forcing mz/intensities to lists."""
    return {
        "sample_index": spec.sample_index,
        "experiment_index": spec.experiment_index,
        "cycle_index": spec.cycle_index,
        "scan_time": spec.scan_time,
        "mz": _to_list(spec.mz),
        "intensities": _to_list(spec.intensities),
        "centroided": spec.centroided,
        "precursor_mz": spec.precursor_mz,
        "isolation_target_mz": spec.isolation_target_mz,
        "isolation_lower_offset": spec.isolation_lower_offset,
        "isolation_upper_offset": spec.isolation_upper_offset,
    }


def serialize_chromatogram(chrom: Any) -> dict[str, Any]:
    """Serialize a Chromatogram or XicChromatogram (both have times/intensities)."""
    out: dict[str, Any] = {
        "times": _to_list(chrom.times),
        "intensities": _to_list(chrom.intensities),
    }
    # XicChromatogram-specific fields
    for attr in ("xic_id", "sample_key", "mz_lower", "mz_upper", "status",
                 "experiment_index", "ms_level"):
        if hasattr(chrom, attr):
            out[attr] = getattr(chrom, attr)
    return out


def unified_to_dict(up: Any) -> dict[str, Any]:
    """Project a UnifiedPeak (a property view, not a dataclass) to a flat row.

    This defines the front-end's "results row" shape (``UnifiedPeakDTO`` in
    docs/types.ts).
    """
    mass_error = up.mass_error
    return {
        "name": up.name,
        "formula": up.formula,
        "sample_index": up.sample_index,
        "compound_index": up.compound_index,
        "area": up.area,
        "retention_time": up.retention_time,
        "height": up.height,
        "signal_to_noise": up.signal_to_noise,
        "found_mass": up.found_mass,
        "found_rt": up.found_rt,
        "mass_error_ppm": (mass_error * 1e6) if mass_error is not None else None,
        "isotope_diff": up.isotope_diff,
        "rt_diff": up.rt_diff,
        "contains_msms": up.contains_msms,
        "has_been_calculated": up.has_been_calculated,
        "valid_integration": up.valid_integration,
        "is_valid": up.is_valid(),
        "library_hits": [dataclass_to_dict(h) for h in up.library_hits],
    }
