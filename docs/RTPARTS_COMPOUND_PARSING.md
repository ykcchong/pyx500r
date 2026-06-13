# RTParts Compound Parsing — Implementation Notes

## What exists in the file

A `.qsession` file stores compound definitions in the `RTParts` table as a custom
binary stream produced by `IterativeSerializer` (Clearcore2.Data).  The stream is
**not** plain BinaryFormatter — it is a text-tag protocol where primitives are raw
binary and complex objects are wrapped in BinaryFormatter sub-blobs.

## What we decode today

The parser in `src/pyx500r/rtparts.py` successfully reads **all 2471 compounds**
from the sample file, including:

| Field | Type | Notes |
|-------|------|-------|
| `name` | string | |
| `group_name` | string | |
| `formula` | string \| None | |
| `charge_formula` | string \| None | e.g. `"[M+H]+"` |
| `adduct_formula` | string \| None | often `None` |
| `precursor_mass` | double | `-1.0` when not set |
| `fragment_mass` | double | `-1.0` when not set |
| `extraction_type` | int | `1` = XIC extraction in sample file |
| `period` | int | |
| `experiment` | int | |
| `extraction_values1` | double[] | m/z window lower bound (length 1) |
| `extraction_values2` | double[] | m/z window upper bound (length 1) |
| `mz_lower` / `mz_upper` | float | derived from `ev1[0]` / `ev2[0]` |
| `is_analyte` | bool | `True` for all 2471 compounds in this file |
| `is_reportable` | bool | |
| `is_non_targeted` | bool | |
| `is_summed` | bool | |
| `is_from_multi_period_data` | bool | |
| `isotope_index` | int | |
| `expected_mw` | double | |
| `units` | string | often empty |
| `comment` | string \| None | |
| `internal_std_name` | string \| None | analyte-only |
| `regression_area` | bool \| None | analyte-only |
| `regression_type` | int \| None | analyte-only |
| `regression_weighting` | int \| None | analyte-only |
| `use_auto_regression` | bool \| None | analyte-only, version ≥ 2 |

## What is skipped

### 1. `fIntegrationParameters`

A 1009-byte BinaryFormatter blob containing ~30 integration settings
(smoothing, baseline, peak detection, S/N threshold, noise regions, etc.).

**Why skipped:** it is an `ISerializable` object whose type lives in
`Clearcore2.QuantLibrary.dll`.  The assembly is not loadable on macOS (missing
WPF/PresentationFramework dependencies), so `BinaryFormatter.Deserialize()` fails
with `SerializationException: Unable to find assembly`.

**C# source (decompiled):** `decompiled/SciexAnalytics/Clearcore2.QuantLibrary/Integration/IntegrationParameters.cs`

Key fields inside the blob (version 17):
- `fSmoothingHalfWindow`, `fSmoothingType` (Gaussian/Savitzky-Golay)
- `fExpectedRT`, `fRTHalfWindow`
- `fMinimumHeight`, `fMinimumWidth`
- `fReportLargest`, `fUpdateRT`
- `fXicWidth` (observed default 0.02)
- `fSignalToNoiseThreshold`
- `fAlgorithmInputType`, `fRetentionTimeMode`
- `fIsParameterLockingApplied` + 20+ min/max enable flags
- `fIsIntactQuant` + intact-quant sub-fields (resolution, peak width, etc.)
- `signalToNoiseAlgorithm` + noise regions (`noiseRegion{N}Begin/End/IsManuallySet`)

### 2. `fAcquisitionIndices`

Always `null` in the sample file.  Uses `SerializeOneObject` → BinaryFormatter array.

### 3. `fSummedCompounds`

Always `null` in the sample file.  Uses `SerializeOneObject` → BinaryFormatter string array.

## The "green / yellow / red" question

There is **no dedicated color/status enum** in the qsession binary format.
MultiQuant's traffic-light indicators are computed at runtime from:
1. Integration results (area, height, RT deviation vs. expected)
2. Qualifier ion ratios
3. Concentration vs. calibration curve

The only place "green" appears in the sample file is inside a compound **name**:
```
Fluoxetine_low CE cfm if mass and RT is green
```
This is user-entered text, not a machine-readable flag.

## BinaryFormatter array parsing trick

`double[]` arrays (`fExtractionValues1`, `fExtractionValues2`) are parsed with
pythonnet's `BinaryFormatter` by reading a 100-byte probe chunk and checking
`MemoryStream.Position` after deserialization.  Each array consumes exactly
**36 bytes** in this file format (1 null-flag + 35 bytes of BinaryFormatter data).

## File-format constants (empirical)

| Constant | Value | Meaning |
|----------|-------|---------|
| `INTEGRATION_PARAMETERS_BLOB_SIZE` | 1009 | Bytes to skip after `fIntegrationParameters` name tag |
| `EXTRACTION_ARRAY_SIZE` | 36 | Bytes per `double[]` extraction-value array |
| Compound base version | 18 | `QuantCompound.SerializeIteratively()` version |
| Analyte sub-version | 2 | `QuantAnalyte` extension version |
| MultiData version | 17 | Overall stream version |
| QuantMethod version | 31 | Method version |
