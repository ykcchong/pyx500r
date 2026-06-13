"""Pure-Python parser for RTParts (compound definitions) in .qsession files.

The RTParts table contains a custom binary serialization produced by
``IterativeSerializer`` (Clearcore2.Data).  Primitive fields are raw
binary; complex objects (``double[]``, ``IntegrationParameters``) are
wrapped in .NET ``BinaryFormatter`` blobs.

This module uses a pure-Python BinaryFormatter reimplementation
(:mod:`pyx500r.binaryformatter`) — no .NET / pythonnet dependency.
"""
from __future__ import annotations

import io
import json
import os
import sqlite3
import struct
import subprocess
import tempfile
from pathlib import Path
import re
from typing import Any

def _find_gap_parser() -> str | None:
    """Locate the ``bf_gap_parser2.exe`` C# batch binary."""
    module_dir = Path(__file__).resolve().parent
    candidates = [
        module_dir / "bf_gap_parser2.exe",
        module_dir.parent / "bf_gap_parser2.exe",
        module_dir.parent.parent / "bf_gap_parser2.exe",
        Path("bf_gap_parser2.exe"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    for env_path in os.environ.get("PATH", "").split(os.pathsep):
        p = Path(env_path) / "bf_gap_parser2.exe"
        if p.exists():
            return str(p)
    return None


def _parse_gap(raw: bytes, start: int, end: int) -> dict[tuple[int, int], dict[str, Any]] | None:
    """Parse XicManagerXic objects from the gap using the C# batch parser."""
    exe = _find_gap_parser()
    if exe is None:
        return None

    with tempfile.NamedTemporaryFile(delete=False, suffix=".gap") as tmp:
        tmp.write(raw[start:end])
        tmp_path = tmp.name

    try:
        result = subprocess.run(
            ["mono", exe, tmp_path, "0"],
            capture_output=True,
            text=True,
            timeout=120,
        )
        if result.returncode != 0:
            return None
        objs = json.loads(result.stdout)
        lookup: dict[tuple[int, int], dict[str, Any]] = {}
        for obj in objs:
            lookup[(obj["sample"], obj["compound"])] = obj["data"]
        return lookup
    except Exception:
        return None
    finally:
        os.unlink(tmp_path)

from .binaryformatter import (
    parse_bf,
    parse_bf_with_consumed,
    parse_datetime,
    parse_double_array,
    parse_float_array,
    parse_hashtable,
    parse_integration_parameters,
    parse_int_array,
)
# Empirically determined for this file-format version.
# After the "fIntegrationParameters" name tag the stream contains
#   <1-byte null-flag> + <BinaryFormatter data>
# The total is always 1009 bytes before the next name tag.
_INTEGRATION_PARAMETERS_BLOB_SIZE = 1009

# ``double[]`` extraction-value arrays always consume 36 bytes in the
# observed format (1 null-flag + 35 bytes of BinaryFormatter data).
_EXTRACTION_ARRAY_SIZE = 36


# ---------------------------------------------------------------------------
# Low-level stream helpers
# ---------------------------------------------------------------------------
def _read_name(stream: io.BytesIO) -> tuple[str, bool]:
    b = stream.read(1)[0]
    is_null = b > 127
    name_len = b & 127
    return stream.read(name_len).decode("ascii"), is_null

def _read_tagged_double(stream: io.BytesIO) -> float:
    """Read a double value, consuming the preceding name tag (no verify)."""
    b = stream.read(1)[0]
    name_len = b & 127
    stream.seek(stream.tell() + name_len)  # skip name bytes
    if b > 127:  # null
        return 0.0
    return struct.unpack("<d", stream.read(8))[0]


def _read_tagged_bool(stream: io.BytesIO) -> bool:
    """Read a bool value, consuming the preceding name tag (no verify)."""
    b = stream.read(1)[0]
    name_len = b & 127
    stream.seek(stream.tell() + name_len)  # skip name bytes
    if b > 127:  # null
        return False
    return stream.read(1)[0] != 0


def _read_string(stream: io.BytesIO) -> str:
    length = struct.unpack("<H", stream.read(2))[0]
    return stream.read(length).decode("utf-8") if length else ""


def _read_int(stream: io.BytesIO) -> int:
    return struct.unpack("<i", stream.read(4))[0]


def _read_double(stream: io.BytesIO) -> float:
    return struct.unpack("<d", stream.read(8))[0]


def _read_bool(stream: io.BytesIO) -> bool:
    return stream.read(1)[0] != 0


def _read_version(stream: io.BytesIO) -> int:
    return struct.unpack("<h", stream.read(2))[0]


# ---------------------------------------------------------------------------
# BinaryFormatter helpers
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# BinaryFormatter helpers (pure-Python, no pythonnet required)
# ---------------------------------------------------------------------------
def _parse_bf_double_array(stream: io.BytesIO) -> list[float] | None:
    """Parse a BinaryFormatter ``double[]`` from the stream."""
    start = stream.tell()
    raw = stream.getvalue()
    arr, consumed = parse_double_array(raw, start)
    if arr is not None:
        stream.seek(start + consumed)
        return arr
    stream.read(36)
    return None


def _parse_bf_int_array(stream: io.BytesIO) -> list[int] | None:
    """Parse a BinaryFormatter ``int[]`` from the stream."""
    start = stream.tell()
    raw = stream.getvalue()
    arr, consumed = parse_int_array(raw, start)
    if arr is not None:
        stream.seek(start + consumed)
        return arr
    stream.read(36)
    return None


def _parse_bf_integration_parameters(stream: io.BytesIO) -> dict[str, Any] | None:
    """Parse a BinaryFormatter ``IntegrationParameters`` from the stream."""
    start = stream.tell()
    raw = stream.getvalue()
    result, consumed = parse_integration_parameters(raw[start:start + 1024])
    if result:
        stream.seek(start + consumed)
        return result
    stream.read(1009)
    return None


def _hashtable_to_dict(obj: Any) -> dict[Any, Any] | None:
    """Convert a parsed ``System.Collections.Hashtable`` object into a dict.

    The BinaryFormatter representation stores parallel ``Keys`` / ``Values``
    arrays (SerializationInfo form). Zip them into a plain ``{key: value}``
    mapping, dropping the internal ``__class``/bucket bookkeeping. Returns
    ``None`` if *obj* is not a recognisable hashtable.
    """
    if not isinstance(obj, dict):
        return None
    keys = obj.get("Keys")
    values = obj.get("Values")
    if isinstance(keys, list) and isinstance(values, list):
        return {k: v for k, v in zip(keys, values) if k is not None}
    # Older bucket-based / already-flat dict: strip bookkeeping keys.
    cls = obj.get("__class", "")
    if "Hashtable" in str(cls):
        return {
            k: v
            for k, v in obj.items()
            if k not in ("__class", "LoadFactor", "Version", "Comparer",
                         "HashCodeProvider", "HashSize", "Keys", "Values")
        }
    return None


def _parse_bf_hashtable(stream: io.BytesIO) -> dict[Any, Any] | None:
    """Parse a BinaryFormatter ``Hashtable`` from the stream.

    Parses the *full* object graph (so that ``MemberReference`` placeholders to
    the ``Keys``/``Values`` arrays are resolved) and zips it into a plain dict.
    """
    start = stream.tell()
    raw = stream.getvalue()
    # Quick check that a hashtable starts here, using a small slice.
    probe = parse_hashtable(raw[start : start + 500])
    if not probe:
        stream.read(50)
        return None
    # Re-parse the full graph from *start* so Keys/Values refs resolve.
    obj, consumed = parse_bf_with_consumed(raw[start:])
    stream.seek(start + consumed)
    result = _hashtable_to_dict(obj)
    if result is not None:
        return result
    # Fall back to the (possibly ref-laden) probe result.
    return probe
def _parse_bf_float_array(stream: io.BytesIO) -> list[float] | None:
    """Parse a BinaryFormatter ``float[]`` from the stream."""
    start = stream.tell()
    raw = stream.getvalue()
    arr, consumed = parse_float_array(raw, start)
    if arr is not None:
        stream.seek(start + consumed)
        return arr
    stream.read(36)
    return None


def _parse_bf_datetime(stream: io.BytesIO) -> Any:
    """Parse a BinaryFormatter ``DateTime`` from the stream."""
    start = stream.tell()
    raw = stream.getvalue()
    dt = parse_datetime(raw, start)
    if dt is not None:
        _, consumed = parse_bf_with_consumed(raw[start:start + 20])
        stream.seek(start + consumed)
        return dt
    stream.read(9)
    return None


def _parse_bf_object(stream: io.BytesIO) -> Any:
    """Parse an arbitrary BinaryFormatter object from the stream."""
    start = stream.tell()
    raw = stream.getvalue()
    obj, consumed = parse_bf_with_consumed(raw, start)
    stream.seek(start + consumed)
    return obj


# ---------------------------------------------------------------------------
def _decode_compound(stream: io.BytesIO) -> dict[str, Any]:
    c: dict[str, Any] = {}
    c["version"] = _read_version(stream)
    c["is_analyte"] = stream.read(1)[0] == 1

    n, nl = _read_name(stream)
    assert n == "fName", f"expected fName, got {n!r}"
    c["name"] = None if nl else _read_string(stream)

    n, nl = _read_name(stream)
    assert n == "fGroupName"
    c["group_name"] = None if nl else _read_string(stream)

    n, nl = _read_name(stream)
    assert n == "fUnits"
    c["units"] = None if nl else _read_string(stream)

    n, nl = _read_name(stream)
    assert n == "fComment"
    c["comment"] = None if nl else _read_string(stream)

    n, nl = _read_name(stream)
    assert n == "fIsReportable"
    c["is_reportable"] = None if nl else _read_bool(stream)

    n, nl = _read_name(stream)
    assert n == "fExtractionType"
    c["extraction_type"] = None if nl else _read_int(stream)

    n, nl = _read_name(stream)
    assert n == "fPeriod"
    c["period"] = None if nl else _read_int(stream)

    n, nl = _read_name(stream)
    assert n == "fExperiment"
    c["experiment"] = None if nl else _read_int(stream)

    n, nl = _read_name(stream)
    assert n == "fAcquisitionIndices"
    c["acquisition_indices"] = None if nl else _parse_bf_int_array(stream)

    n, nl = _read_name(stream)
    assert n == "fExtractionValues1"
    c["extraction_values1"] = None if nl else _parse_bf_double_array(stream)

    n, nl = _read_name(stream)
    assert n == "fExtractionValues2"
    c["extraction_values2"] = None if nl else _parse_bf_double_array(stream)

    n, nl = _read_name(stream)
    assert n == "fITraqTagMass"
    c["itraq_tag_mass"] = None if nl else _read_double(stream)

    n, nl = _read_name(stream)
    assert n == "fIntegrationParameters"
    c["integration_parameters"] = None if nl else _parse_bf_integration_parameters(stream)

    n, nl = _read_name(stream)
    assert n == "fADCChannels"
    c["adc_channels"] = None if nl else _read_string(stream)

    n, nl = _read_name(stream)
    assert n == "fIsNonTargeted"
    c["is_non_targeted"] = None if nl else _read_bool(stream)

    n, nl = _read_name(stream)
    assert n == "fFormula"
    c["formula"] = None if nl else _read_string(stream)

    n, nl = _read_name(stream)
    assert n == "fIsotopeIndex"
    c["isotope_index"] = None if nl else _read_int(stream)

    n, nl = _read_name(stream)
    assert n == "fChargeFormula"
    c["charge_formula"] = None if nl else _read_string(stream)

    n, nl = _read_name(stream)
    assert n == "fAdductFormula"
    c["adduct_formula"] = None if nl else _read_string(stream)

    n, nl = _read_name(stream)
    assert n == "fPrecursorMass"
    c["precursor_mass"] = None if nl else _read_double(stream)

    n, nl = _read_name(stream)
    assert n == "fSummedCompounds"
    c["summed_compounds"] = None if nl else _parse_bf_int_array(stream)

    n, nl = _read_name(stream)
    assert n == "fIsSummed"
    c["is_summed"] = None if nl else _read_bool(stream)

    n, nl = _read_name(stream)
    assert n == "fFragmentMass"
    c["fragment_mass"] = None if nl else _read_double(stream)

    n, nl = _read_name(stream)
    assert n == "fStartStopMass"
    c["start_stop_mass"] = None if nl else _read_string(stream)

    n, nl = _read_name(stream)
    assert n == "fIsFromMultiPeriodData"
    c["is_from_multi_period_data"] = None if nl else _read_bool(stream)

    n, nl = _read_name(stream)
    assert n == "_fAutoSampleIndex"
    c["auto_sample_index"] = None if nl else _read_int(stream)

    n, nl = _read_name(stream)
    assert n == "fExpectedMW"
    c["expected_mw"] = None if nl else _read_double(stream)

    # --- Analyte / ISTD specific extension ---
    sub_ver = _read_version(stream)
    c["sub_version"] = sub_ver
    if c["is_analyte"]:
        n, nl = _read_name(stream)
        assert n == "fInternalStdName"
        c["internal_std_name"] = None if nl else _read_string(stream)

        n, nl = _read_name(stream)
        assert n == "fRegressionArea"
        c["regression_area"] = None if nl else _read_bool(stream)

        n, nl = _read_name(stream)
        assert n == "fRegressionType"
        c["regression_type"] = None if nl else _read_int(stream)

        n, nl = _read_name(stream)
        assert n == "fRegressionWeighting"
        c["regression_weighting"] = None if nl else _read_int(stream)

        if sub_ver >= 2:
            n, nl = _read_name(stream)
            assert n == "_useAutoRegression"
            c["use_auto_regression"] = None if nl else _read_bool(stream)

    return c


# ---------------------------------------------------------------------------
# BinaryFormatter skip helpers
# ---------------------------------------------------------------------------
_KNOWN_FIELD_NAMES: set[str] = {
    # MultiPeak fields
    "_peakIndex", "_integrationParameters", "_use", "_peakComment",
    "_actualConcentration", "_failedQuery", "_customFields", "_customPeakFields",
    "_validIntegration", "_modified", "_retentionTime", "_area", "_correctedArea",
    "_height", "_correctedHeight", "_startRT", "_startY", "_endRT", "_endY",
    "_halfHeightStartRT", "_halfHeightEndRT", "_noise", "_profile", "_profileType",
    "_peakType", "_apexRT", "_apexY", "_regionArea", "_regionHeight",
    "_sMrmRetentionTimeShift", "_rowColour", "_rowHidden", "_startX5PctHeight",
    "_endX5PctHeight", "_startX10PctHeight", "_endX10PctHeight",
    "_pointsAcrossBaseline", "_pointsAcrossHalfHeight", "_overrideExperimentIndex",
    "fStdAddnActualConcentration", "fExtractedMsMs", "_reportable", "_molecularWeight",
    "_originalArea", "_superGroupId",
    # MultiSample fields
    "_multiPeaks", "_samples", "_fManSample", "_sampleOffset", "_sampleName", "_sampleId",
    "_rack", "_plate", "_vial", "_dateTime", "_sampleType", "_sampleComment",
    "_dilutionFactor", "_injectionVolume", "_userName", "_acqMethodName",
    "_instrumentName", "_instrumentSerialNumber", "_batchName", "_barcode",
    "_scannedBarcode", "_autosamplerMethodSupportsBarcode", "_sampleComparison",
    "_msMethod", "_lcMethod", "_sampleSignature", "_isTriggeredMsMs",
    "_assayInfo", "_experimentPolarities", "_timeSinceFirstSampleMin",
    "_timeSinceLastSampleSec", "_differenceFromAverageSampleTimeSecs", "_transferData",
    # MultiData tail fields
    "_customFieldNames", "_externalCal", "_isotopePatterns", "_comment",
    "_projectPath", "IsSmrmDataInSession", "_customFieldEditables",
    "_customFieldVisibles", "_customFieldTypes", "_customFieldFormula",
    "_customFieldFormatted", "_combinedRuleColumns",
}

_NAME_TAGS: dict[str, tuple[bytes, bytes]] = {}
for _name in _KNOWN_FIELD_NAMES:
    _nb = _name.encode("ascii")
    _nlen = len(_nb)
    _NAME_TAGS[_name] = (bytes([_nlen]) + _nb, bytes([_nlen + 128]) + _nb)


def _is_valid_name_tag(raw: bytes, pos: int) -> str | None:
    """Return the field name if *pos* points to a valid name tag, else None."""
    if pos >= len(raw):
        return None
    b = raw[pos]
    nlen = b & 0x7F
    if nlen == 0 or pos + 1 + nlen > len(raw):
        return None
    try:
        name = raw[pos + 1 : pos + 1 + nlen].decode("ascii")
    except UnicodeDecodeError:
        return None
    if name not in _KNOWN_FIELD_NAMES:
        return None
    return name


def _skip_bf_object(
    raw: bytes,
    pos: int,
    next_expected: list[str],
    max_scan: int = 50_000,
) -> int:
    """Return the byte position of the next expected field name tag.

    Uses fast ``bytes.find()`` instead of regex.
    """
    end = min(pos + max_scan, len(raw))
    for name in next_expected:
        tag, tag_null = _NAME_TAGS[name]
        idx = raw.find(tag, pos, end)
        if idx != -1:
            return idx
        idx = raw.find(tag_null, pos, end)
        if idx != -1:
            return idx
    raise ValueError(
        f"Could not find any of {next_expected!r} after position {pos}"
    )


# ---------------------------------------------------------------------------
# MultiSample / MultiPeak decoders
# ---------------------------------------------------------------------------
def _decode_multipeak(
    stream: io.BytesIO,
    sample_index: int,
    compound_index: int,
) -> dict[str, Any]:
    """Decode one MultiPeak from the iterative serialization stream."""
    raw = stream.getvalue()
    version = _read_version(stream)
    # Local bindings – avoid repeated global/attribute lookups
    read_name = _read_name
    read_double = _read_double
    read_bool = _read_bool
    read_int = _read_int
    skip_bf = _skip_bf_object
    p: dict[str, Any] = {
        "sample_index": sample_index,
        "compound_index": compound_index,
        "version": version,
    }

    # --- Core fields (present in all versions) ---
    name, is_null = _read_name(stream)
    if name == "_peakIndex" and not is_null:
        p["peak_index"] = _read_int(stream)
    elif name == "_peakIndex":
        p["peak_index"] = -1
    else:
        raise ValueError(f"Expected _peakIndex, got {name!r}")

    name, is_null = _read_name(stream)
    if name == "_integrationParameters":
        if not is_null:
            stream.seek(_skip_bf_object(raw, stream.tell(), ["_use"]))
    else:
        raise ValueError(f"Expected _integrationParameters, got {name!r}")

    name, is_null = _read_name(stream)
    if name == "_use" and not is_null:
        p["use_for_calibration"] = _read_bool(stream)
    elif name == "_use":
        p["use_for_calibration"] = False
    else:
        raise ValueError(f"Expected _use, got {name!r}")

    _read_opt_str = _read_optional_string
    _read_opt_bool = _read_optional_bool
    _read_opt_double = _read_optional_double

    p["peak_comment"] = _read_opt_str(stream, "_peakComment")
    p["actual_concentration"] = _read_tagged_double(stream)
    p["failed_query"] = _read_tagged_bool(stream)

    # _customFields (BinaryFormatter Hashtable)
    n, nl = _read_name(stream)
    if n == "_customFields" and not nl:
        p["custom_fields"] = _parse_bf_hashtable(stream)
    elif n != "_customFields":
        raise ValueError(f"Expected _customFields, got {n!r}")
    # _customPeakFields (BinaryFormatter Hashtable, may be absent in old versions)
    pos = stream.tell()
    b = stream.read(1)[0]
    stream.seek(pos)
    is_null_tag = b >= 128
    nlen = b & 0x7F
    if nlen == len("_customPeakFields") and raw[pos + 1 : pos + 1 + nlen] == b"_customPeakFields":
        n, nl = _read_name(stream)
        if not nl:
            p["custom_peak_fields"] = _parse_bf_hashtable(stream)
    # else: _customPeakFields is not present, stream stays at current position

    p["valid_integration"] = _read_tagged_bool(stream)
    p["modified"] = _read_tagged_bool(stream)
    p["retention_time"] = _read_tagged_double(stream)
    p["area"] = _read_tagged_double(stream)
    p["corrected_area"] = _read_tagged_double(stream)
    p["height"] = _read_tagged_double(stream)
    p["corrected_height"] = _read_tagged_double(stream)
    p["start_rt"] = _read_tagged_double(stream)
    p["start_y"] = _read_tagged_double(stream)
    p["end_rt"] = _read_tagged_double(stream)
    p["end_y"] = _read_tagged_double(stream)
    p["half_height_start_rt"] = _read_tagged_double(stream)
    p["half_height_end_rt"] = _read_tagged_double(stream)
    p["noise"] = _read_tagged_double(stream)

    # _profile (BinaryFormatter float[])
    n, nl = _read_name(stream)
    if n == "_profile" and not nl:
        stream.seek(_skip_bf_object(raw, stream.tell(), ["_profileType"]))
    elif n != "_profile":
        raise ValueError(f"Expected _profile, got {n!r}")

    p["profile_type"] = _read_optional_int(stream, "_profileType")
    p["peak_type"] = _read_optional_int(stream, "_peakType")
    p["apex_rt"] = _read_tagged_double(stream)
    p["apex_y"] = _read_tagged_double(stream)
    p["region_area"] = _read_tagged_double(stream)
    p["region_height"] = _read_tagged_double(stream)
    p["s_mrm_retention_time_shift"] = _read_tagged_bool(stream)

    # _rowColour (BinaryFormatter Color)
    n, nl = _read_name(stream)
    if n == "_rowColour" and not nl:
        stream.seek(_skip_bf_object(raw, stream.tell(), ["_rowHidden", "_startX5PctHeight"]))
    elif n != "_rowColour":
        raise ValueError(f"Expected _rowColour, got {n!r}")

    p["row_hidden"] = _read_tagged_bool(stream)

    # --- Version-dependent tail fields ---
    if version >= 15:
        p["start_x5_pct_height"] = _read_tagged_double(stream)
        p["end_x5_pct_height"] = _read_tagged_double(stream)
        p["start_x10_pct_height"] = _read_tagged_double(stream)
        p["end_x10_pct_height"] = _read_tagged_double(stream)
        p["points_across_baseline"] = _read_optional_int(stream, "_pointsAcrossBaseline")
        p["points_across_half_height"] = _read_optional_int(stream, "_pointsAcrossHalfHeight")

    if version >= 16:
        p["override_experiment_index"] = _read_optional_int(stream, "_overrideExperimentIndex")

    if version >= 17:
        n, nl = _read_name(stream)
        if n == "fStdAddnActualConcentration" and not nl:
            p["std_addn_actual_concentration"] = _read_double(stream)
        elif n == "fStdAddnActualConcentration":
            p["std_addn_actual_concentration"] = 0.0
        else:
            stream.seek(stream.tell() - 1 - len(n))

    if version >= 18:
        n, nl = _read_name(stream)
        if n == "fExtractedMsMs" and not nl:
            p["extracted_ms_ms"] = _read_double(stream)
        elif n == "fExtractedMsMs":
            p["extracted_ms_ms"] = None
        else:
            stream.seek(stream.tell() - 1 - len(n))

    if version >= 19:
        p["reportable"] = _read_tagged_bool(stream)

    if version >= 20:
        p["molecular_weight"] = _read_tagged_double(stream)

    if version >= 21:
        p["original_area"] = _read_tagged_double(stream)

    if version >= 22:
        p["super_group_id"] = _read_opt_str(stream, "_superGroupId")

    # Compute signal-to-noise from height and noise (matches .NET logic)
    noise = p.get("noise", -1.0)
    height = p.get("height", 0.0)
    p["signal_to_noise"] = height / noise if noise > 0.0 else -1.0

    return p


def _read_optional_int(stream: io.BytesIO, field: str) -> int:
    n, nl = _read_name(stream)
    if n != field:
        raise ValueError(f"Expected {field!r}, got {n!r}")
    return 0 if nl else _read_int(stream)

def _read_optional_string(stream: io.BytesIO, field: str) -> str | None:
    n, nl = _read_name(stream)
    if n != field:
        raise ValueError(f"Expected {field!r}, got {n!r}")
    return None if nl else _read_string(stream)


def _read_optional_bool(stream: io.BytesIO, field: str) -> bool:
    n, nl = _read_name(stream)
    if n != field:
        raise ValueError(f"Expected {field!r}, got {n!r}")
    return False if nl else _read_bool(stream)


def _read_optional_double(stream: io.BytesIO, field: str) -> float:
    n, nl = _read_name(stream)
    if n != field:
        raise ValueError(f"Expected {field!r}, got {n!r}")
    return 0.0 if nl else _read_double(stream)


def _read_string_field(stream: io.BytesIO, field: str) -> str | None:
    n, nl = _read_name(stream)
    if n != field:
        raise ValueError(f"Expected {field!r}, got {n!r}")
    return None if nl else _read_string(stream)


def _read_bool_field(stream: io.BytesIO, field: str) -> bool:
    n, nl = _read_name(stream)
    if n != field:
        raise ValueError(f"Expected {field!r}, got {n!r}")
    return False if nl else _read_bool(stream)


def _read_double_field(stream: io.BytesIO, field: str) -> float:
    n, nl = _read_name(stream)
    if n != field:
        raise ValueError(f"Expected {field!r}, got {n!r}")
    return 0.0 if nl else _read_double(stream)


def _decode_multisample(
    stream: io.BytesIO,
    sample_index: int,
    num_compounds: int,
) -> dict[str, Any]:
    """Decode one MultiSample from the iterative serialization stream."""
    raw = stream.getvalue()
    version = _read_version(stream)
    s: dict[str, Any] = {"index": sample_index, "version": version}
    # _multiPeaks array header
    n, nl = _read_name(stream)
    if n != "_multiPeaks":
        raise ValueError(f"Expected _multiPeaks, got {n!r}")
    if nl:
        s["peaks"] = []
        return s

    num_peaks = _read_int(stream)
    peaks: list[dict[str, Any]] = []
    for ci in range(num_peaks):
        peaks.append(_decode_multipeak(stream, sample_index, ci))
    s["peaks"] = peaks

    # Sample metadata fields – bind to module-level helpers for speed
    _read_str_fld = _read_string_field
    _read_bool_fld = _read_bool_field
    _read_dbl_fld = _read_double_field

    # _fManSample (SampleLocator) – BinaryFormatter, skip
    n, nl = _read_name(stream)
    if n == "_fManSample" and not nl:
        stream.seek(_skip_bf_object(raw, stream.tell(), ["_sampleOffset", "_sampleName"]))
    elif n != "_fManSample":
        raise ValueError(f"Expected _fManSample, got {n!r}")

    s["sample_offset"] = _read_dbl_fld(stream, "_sampleOffset")
    s["sample_name"] = _read_str_fld(stream, "_sampleName")
    s["sample_id"] = _read_str_fld(stream, "_sampleId")
    s["rack"] = _read_str_fld(stream, "_rack")
    s["plate"] = _read_str_fld(stream, "_plate")
    s["vial"] = _read_str_fld(stream, "_vial")

    # _dateTime (BinaryFormatter DateTime)
    n, nl = _read_name(stream)
    if n == "_dateTime" and not nl:
        stream.seek(_skip_bf_object(raw, stream.tell(), ["_sampleType"]))
    elif n != "_dateTime":
        raise ValueError(f"Expected _dateTime, got {n!r}")

    s["sample_type"] = _read_int_field(stream, "_sampleType")
    s["sample_comment"] = _read_str_fld(stream, "_sampleComment")
    s["dilution_factor"] = _read_dbl_fld(stream, "_dilutionFactor")
    s["injection_volume"] = _read_dbl_fld(stream, "_injectionVolume")
    s["user_name"] = _read_str_fld(stream, "_userName")
    s["acq_method_name"] = _read_str_fld(stream, "_acqMethodName")
    s["instrument_name"] = _read_str_fld(stream, "_instrumentName")
    s["instrument_serial_number"] = _read_str_fld(stream, "_instrumentSerialNumber")

    if version >= 6:
        s["batch_name"] = _read_str_fld(stream, "_batchName")
        s["barcode"] = _read_str_fld(stream, "_barcode")
        s["scanned_barcode"] = _read_str_fld(stream, "_scannedBarcode")
    if version >= 7:
        s["autosampler_method_supports_barcode"] = _read_bool_fld(stream, "_autosamplerMethodSupportsBarcode")
    if version >= 8:
        s["sample_comparison"] = _read_bool_fld(stream, "_sampleComparison")
    if version >= 9:
        s["ms_method"] = _read_str_fld(stream, "_msMethod")
        s["lc_method"] = _read_str_fld(stream, "_lcMethod")
    if version >= 10:
        s["sample_signature"] = _read_str_fld(stream, "_sampleSignature")
    if version >= 11:
        s["is_triggered_ms_ms"] = _read_bool_fld(stream, "_isTriggeredMsMs")
    if version >= 12:
        n, nl = _read_name(stream)
        if n == "_assayInfo" and not nl:
            stream.seek(_skip_bf_object(raw, stream.tell(), ["_experimentPolarities", "_timeSinceFirstSampleMin", "_msMethod"]))
        elif n != "_assayInfo":
            raise ValueError(f"Expected _assayInfo, got {n!r}")
    if version >= 13:
        n, nl = _read_name(stream)
        if n == "_experimentPolarities" and not nl:
            stream.seek(_skip_bf_object(raw, stream.tell(), ["_timeSinceFirstSampleMin", "_msMethod"]))
        elif n != "_experimentPolarities":
            raise ValueError(f"Expected _experimentPolarities, got {n!r}")
    if version >= 14:
        for tf in ("_timeSinceFirstSampleMin", "_timeSinceLastSampleSec", "_differenceFromAverageSampleTimeSecs"):
            n, nl = _read_name(stream)
            if n == tf and not nl:
                stream.seek(_skip_bf_object(raw, stream.tell(), ["_transferData", "_msMethod"]))
            elif n != tf:
                raise ValueError(f"Expected {tf!r}, got {n!r}")
    if version >= 15:
        n, nl = _read_name(stream)
        if n == "_transferData" and not nl:
            stream.seek(_skip_bf_object(raw, stream.tell(), ["_msMethod"]))
        elif n != "_transferData":
            raise ValueError(f"Expected _transferData, got {n!r}")

    return s


def _read_int_field(stream: io.BytesIO, field: str) -> int:
    n, nl = _read_name(stream)
    if n != field:
        raise ValueError(f"Expected {field!r}, got {n!r}")
    return 0 if nl else _read_int(stream)

def read_multidata(stream: io.BytesIO) -> dict[str, Any]:
    """Parse the full MultiData object graph from an RTParts stream.

    Returns a dictionary with keys:
    * ``md_version`` – MultiData format version
    * ``qm_version`` – QuantMethod format version
    * ``compounds`` – list of compound dicts
    * ``samples`` – list of sample dicts (each with a ``peaks`` list)
    """
    md_version = _read_version(stream)
    qm_version = _read_version(stream)

    # fCompounds
    n, nl = _read_name(stream)
    if n != "fCompounds":
        raise ValueError(f"Expected fCompounds, got {n!r}")
    num_compounds = 0 if nl else _read_int(stream)
    compounds = [_decode_compound(stream) for _ in range(num_compounds)]

    # Remaining QuantMethod fields — parse the gap (XicManagerXic objects)
    raw = stream.getvalue()
    start = stream.tell()
    tag, tag_null = _NAME_TAGS["_samples"]
    idx = raw.find(tag, start)
    if idx == -1:
        idx = raw.find(tag_null, start)
    if idx == -1:
        raise ValueError("Could not find _samples in stream")

    xic_lookup = _parse_gap(raw, start, idx)
    stream.seek(idx)

    # _samples
    n, nl = _read_name(stream)
    if n != "_samples":
        raise ValueError(f"Expected _samples, got {n!r}")
    num_samples = 0 if nl else _read_int(stream)
    samples = [_decode_multisample(stream, i, num_compounds) for i in range(num_samples)]

    return {
        "md_version": md_version,
        "qm_version": qm_version,
        "compounds": compounds,
        "samples": samples,
        "xic_lookup": xic_lookup,
    }


# ---------------------------------------------------------------------------
def load_rtparts_stream(
    conn: sqlite3.Connection,
    at_record_timestamp: int | None = None,
) -> io.BytesIO:
    """Reassemble the concatenated RTParts blob for a given timestamp.

    If *at_record_timestamp* is ``None`` the most recent record is used.
    """
    c = conn.cursor()
    if at_record_timestamp is None:
        c.execute("SELECT ATRecordTimeStamp FROM RTParts ORDER BY ATRecordTimeStamp DESC LIMIT 1")
        row = c.fetchone()
        if row is None:
            raise ValueError("RTParts table is empty")
        at_record_timestamp = row[0]

    c.execute(
        "SELECT PartContent FROM RTParts WHERE ATRecordTimeStamp = ? ORDER BY PartId",
        (at_record_timestamp,),
    )
    all_bytes = b"".join(row[0] for row in c)
    return io.BytesIO(all_bytes)


def read_compounds(stream: io.BytesIO) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Parse *all* compounds from an RTParts stream.

    Returns ``(compounds, metadata)`` where *metadata* contains the
    MultiData version, QuantMethod version and compound count.
    """
    md_version = _read_version(stream)
    qm_version = _read_version(stream)
    n, _ = _read_name(stream)
    assert n == "fCompounds"
    num_compounds = _read_int(stream)

    compounds = [_decode_compound(stream) for _ in range(num_compounds)]
    meta = {"md_version": md_version, "qm_version": qm_version, "count": num_compounds}
    return compounds, meta
