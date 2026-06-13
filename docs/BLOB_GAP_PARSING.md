# RTParts Blob Gap Parsing

> **Date**: 2026-06-10
> **Status**: ✅ Complete — 4,942 XicManagerXic objects parsed end-to-end.

## File context

* `data/qsession/example_session.qsession` — the SQLite database that stores the decrypted session.
* `RTParts` table — contains the serialized MultiData / QuantMethod / sample / peak object graph split into 484 `PartId` rows.
* `load_rtparts_stream()` reassembles the parts into one contiguous 49 MB `io.BytesIO` stream.

## Overall stream layout

| Section | Start offset | Size | Notes |
|---------|-------------|------|-------|
| `md_version` + `qm_version` | 0 | 4 bytes | Two int16 version fields |
| `fCompounds` header | 4 | variable | Name tag + count (2471) |
| **2471 compound definitions** | variable | ~4.1 MB | Parsed by `_decode_compound()` |
| `fDefaultIntegrationParams` | 4,147,379 | 1,009 bytes | BinaryFormatter `DefaultIntegrationParams` |
| **GAP — 30.9 MB** | 4,148,388 | 30,933,856 bytes | See below |
| `_samples` tag | 35,081,235 | — | Start of `MultiSample` / `MultiPeak` data |
| **2 samples × 2471 peaks** | 35,081,235 | ~14 MB | Parsed by `_decode_multisample()` / `_decode_multipeak()` |

> **Key finding**: the Python parser currently skips the entire 30.9 MB gap.  It scans for the `_samples` tag after reading compounds, discarding everything in between.

## What's in the gap

### Object arrangement

The gap contains **4,942** `XicManagerXic` objects = 2471 compounds × 2 samples.

- The objects are laid out **grouped by sample**:
  - **Sample 0**: compounds 0 … 2470 (2471 objects)
  - **Sample 1**: compounds 0 … 2470 (2471 objects)
- Each object carries its own `_index` field matching the compound index.
- Before the first `XicManagerXic` there is a small **78-byte `DateTime`** object (value 0/0).

### Size distribution

| Size cluster | Count | Reason |
|-------------|-------|--------|
| ~5,600–6,000 bytes | ~2,000 | `_hasBeenCalculated = true` — has `_foundAtMass`, `_foundAtRt`, `_area`, etc. |
| ~6,500–6,600 bytes | ~2,900 | `_hasBeenCalculated = false` — fewer populated fields |

> Individual sizes vary because string fields (`_compoundName`, `_formula`, `_adductDescription`, etc.) change length per compound.

### Wire format

Each `XicManagerXic` is a **BinaryFormatter** blob that starts with:

```
00 00 00                      ← BinaryFormatter Header (0x00) + RootId + HeaderId
0C 02 00 00 00 4E ...         ← BinaryLibrary #2  (Sciex.XicManager.Engine)
0C 03 00 00 00 4E ...         ← BinaryLibrary #3  (Clearcore2.QuantLibrary)
05 01 00 00 00 25 ...         ← ClassWithMembersAndTypes for XicManagerXic
```

The class uses a **non-standard `BinaryTypeEnum` mapping** on the wire.  The MS-NRBF mapping (`0=Primitive, 1=String, 2=Object, 3=SystemClass, 4=Class, …`) does **not** match the runtime types.  For example:

- `_isQualifier` is a `Boolean` but is tagged as `SystemClass`
- `_foundAtMass` is a `Double` but is tagged as `String`
- `_librarySearchResults` is an `ArrayList` but is tagged as `Primitive`

Because the pure-Python `_BfReader` relies on the wire tags to read member values, it reads garbage and fails with `UnicodeDecodeError` or `struct.error`.

## C# fallback

The C# `BinaryFormatter` with a `GenericBinder` can parse every object successfully because it uses **runtime type metadata** (reflection over the `GenericObject` class) rather than the wire-format tags.

### Tools

| Tool | Purpose |
|------|---------|
| `bf_cli.exe` | Single-object parser; reads a byte slice, returns `{"success": true, "consumed": N, "data": {...}}` |
| `scan_gap_objects.exe` | Scans a byte range, finds valid BF objects, returns JSON array with offset + consumed + data |

### Example — Famotidine (compound 616, sample 0)

Extracted via `bf_cli` at offset **23,737,029** (within `rtparts_blob.bin`):

```json
{
  "_compoundName": "Famotidine",
  "_baseMass": 337.0449362768,
  "_extractionMass": 338.052212729121,
  "_foundAtMass": 338.052936825565,
  "_foundAtRt": 2.37642436482406,
  "_foundAtRtApex": 2.363,
  "_rt": 2.32,
  "_area": 111209.815178138,
  "_intensity": 20126.9958859162,
  "_index": 616,
  "_hasBeenCalculated": true,
  "_librarySearchResults": {
    "_items": [
      {"_fit": 0.98712201185893, "_reverseFit": 0.999357224681658, "_purity": 0.986515348908501},
      {"_fit": 0.986788488636293, "_reverseFit": 1.0, "_purity": 0.986788488636293}
    ]
  }
}
```

## Fields observed on `XicManagerXic`

### Core compound / mass fields
| Field | Type | Notes |
|-------|------|-------|
| `_version` | int | Usually 7 |
| `_disposed` | bool | Always false |
| `_index` | int | Compound index (0 … 2470) |
| `_compoundName` | string | Same as `CompoundInfo.name` |
| `_formula` | string | Molecular formula |
| `_baseMass` | double | Monoisotopic mass |
| `_extractionMass` | double | `[M+H]+` or adduct mass |
| `_foundAtMass` | double | Observed centroid mass (0 if not calculated) |
| `_foundAtRt` | double | Observed RT at peak start (0 if not calculated) |
| `_foundAtRtApex` | double | Observed RT at apex |
| `_foundAtRtStart` | double | Chromatogram start RT |
| `_foundAtRtEnd` | double | Chromatogram end RT |
| `_rt` | double | Expected / theoretical RT |
| `_area` | double | Integrated area |
| `_intensity` | double | Peak intensity |
| `_charge` | int | Charge state |
| `_adductDescription` | string | e.g. `"[M+H]+"` |
| `_fragmentMass` | double? | MS/MS fragment mass |
| `_foundAtFragmentMass` | double? | Observed fragment mass |

### Flags
| Field | Type |
|-------|------|
| `_hasBeenCalculated` | bool |
| `_include` | bool |
| `_isInternalStandard` | bool |
| `_isNonTargetPeak` | bool |
| `_isQualifier` | bool |
| `_quantifierIndex` | int |
| `_containsMSMS` | bool |
| `_isDeconvolvedUsedForSearch` | bool |
| `_isIDAMsMs` | bool |

### Isotope / extraction
| Field | Type | Notes |
|-------|------|-------|
| `_extractionType` | int enum | 0 = Standard, 1 = Scheduled, 2 = Targeted |
| `_extractionWidth` | double | Mass extraction window (Da) |
| `_expectedRtWidth` | double | RT window (min) |
| `_isotopeOneBasedIndex` | int | 1-based isotope index |
| `_retTimeToIsotopeDifference` | Hashtable | Maps RT → isotope ratio diff |

### MS/MS
| Field | Type | Notes |
|-------|------|-------|
| `_msMsCycle` | int | MS/MS cycle number |
| `_msMsPeriod` | int | MS/MS period |
| `_msMsSpectrum` | object | Raw MS/MS spectrum data |
| `_msMsSpectrumDeconvolved` | object | Deconvolved MS/MS |

### Library search
| Field | Type | Notes |
|-------|------|-------|
| `_librarySearchResults` | ArrayList | `LibrarySearchResult` objects |
| `_formulaFindResults` | ArrayList | `FormulaFindResult` objects |
| `_librarySearchCriteriaStartRt` | double | RT window start for search |
| `_librarySearchCriteriaEndRt` | double | RT window end for search |

### Result detail fields (`LibrarySearchResult`)
| Field | Type |
|-------|------|
| `_fit` | double |
| `_reverseFit` | double |
| `_purity` | double |
| `_librarySearchResultId` | int |
| `_ionIndex` | int |

### Identification
| Field | Type |
|-------|------|
| `_id` | Guid (as `_a` … `_k` byte fields) |

## Mapping to Python objects

The intended mapping is:

```
(sample_index, compound_index) → XicManagerXic
```

Where:
- `sample_index` = 0 or 1 (matches `QuantSampleInfo.sample_index`)
- `compound_index` = 0 … 2470 (matches `CompoundInfo.compound_index`)

Each `QuantPeakInfo` (from `_decode_multipeak`) should be enriched with its corresponding `XicManagerXic` fields:
- `measured_mass` → `_foundAtMass`
- `measured_rt` → `_foundAtRt`
- `library_search_results` → `_librarySearchResults._items`
- `isotope_pattern` → `_retTimeToIsotopeDifference`
- `msms_info` → `_msMsCycle`, `_msMsPeriod`, `_containsMSMS`
- etc.

## Implementation (2026-06-10)

The gap is now fully parsed using a C# helper (`parse_xics.exe`) that reads the
entire RTParts stream and extracts all 4,942 `XicManagerXic` blobs in ~0.5 s.

### Architecture

1. **`parse_xics.cs` / `parse_xics.exe`** (mono): Scans the raw byte stream for
   ``XicManagerXic`` class-name strings, backtracks to find the containing
   BinaryFormatter header, deserialises via ``GenericBinder``→
   ``GenericObject`` (ISerializable shim), and emits one compact JSON line per
   blob (keys: ``o`` offset, ``c`` consumed, ``n`` name, ``x`` index, ``m``
   found mass, ``r`` found RT, ``a`` area, ``I`` intensity, ``h``
   has_been_calculated, ``q`` is_qualifier, ``s`` is_internal_standard, ``A``
   apex RT, ``b`` base mass, ``l`` library results).

2. **`src/pyx500r/xic_gap.py`**: Python wrapper that calls ``parse_xics.exe``
   via subprocess, converts the compact JSON to .NET-style field names, and
   builds a `(sample_index, compound_index)` lookup dict.

3. **`QSessionReader._load_multidata`** now calls ``parse_xic_blobs`` +
   ``build_xic_index`` after ``read_multidata``, storing the result as
   ``data["xic_lookup"]``.  If mono is unavailable the gap is silently skipped.

4. **``QuantPeakInfo.xic_result``** now carries the full ``dict`` of
   ``XicManagerXic`` fields (``_foundAtMass``, ``_foundAtRt``, ``_area``,
   ``_intensity``, ``_librarySearchResults``, etc.) for every peak.

### Performance

| Step | Time |
|------|------|
| Decrypt + load | 0.8 s |
| Parse compounds (2,471) | 0.1 s |
| Parse samples + peaks (2 × 2,471) | 2.5 s |
| Parse XIC gap (4,942 blobs via mono) | 1.0 s |
| **Total** | **4.8 s** |

4,942 / 4,942 peaks have XIC data.  57 peaks have library-search results.
