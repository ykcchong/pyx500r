# pyx500r — API Reference

> Accurate as of `pyx500r` **0.2.0**. Every signature and field in this document
> is generated from and verified against the source in `src/pyx500r/`.

`pyx500r` reads SCIEX `.wiff2` raw acquisition files and `.qsession`
quantitation-results files in **pure Python** (no .NET / SCIEX DLL dependency).
It is designed to be embedded in services and GUIs — see
[`GUI_INTEGRATION.md`](./GUI_INTEGRATION.md) for a FastAPI/Uvicorn + TypeScript
recipe and [`types.ts`](./types.ts) for ready-made TypeScript interfaces.

## Table of contents

- [Installation & feature flags](#installation--feature-flags)
- [Architecture at a glance](#architecture-at-a-glance)
- [WIFF2 reader (`WiffReader`)](#wiff2-reader-wiffreader)
- [QSession reader (`QSessionReader`)](#qsession-reader-qsessionreader)
- [Bridge (`WiffQSessionBridge`)](#bridge-wiffqsessionbridge)
- [UnifiedPeak](#unifiedpeak)
- [Precursor index](#precursor-index)
- [Library search](#library-search)
- [TOF codec & calibration](#tof-codec--calibration)
- [Centroiding](#centroiding)
- [Decryption primitives](#decryption-primitives)
- [Data model reference](#data-model-reference)
- [Serialization notes (for web APIs)](#serialization-notes-for-web-apis)
- [Concurrency & lifecycle](#concurrency--lifecycle)
- [Performance](#performance)
- [CLI tools](#cli-tools)
- [Exceptions](#exceptions)

---

## Installation & feature flags

```bash
pip install -e .            # pure-Python: crypto + models + readers + slow TOF
pip install -e .[numba]     # + numba/numpy JIT — ~10× faster TOF & centroiding
pip install -e .[cli]       # + tqdm for CLI progress bars
pip install -e .[all]       # everything
pip install -e .[dev]       # + pytest, ruff
```

| Capability | Base install | Needs `numpy` | Needs `numba` |
|------------|:---:|:---:|:---:|
| Decrypt `.wiff2` / `.qsession` | ✅ | | |
| Read metadata / compounds / peaks | ✅ | | |
| Decode TOF spectra (pure Python) | ✅ | | |
| Fast TOF decode / `return_arrays=True` | | ✅ | ✅ (best) |
| `PrecursorIndex`, `build_precursor_index` | | ✅ | |
| `centroid_spectrum` (fast path) | falls back to NumPy/Python | ✅ | ✅ (best) |

`numpy` is pulled in by both the `numba` extra and transitively used by the
index/centroiding code. The base install degrades gracefully when it is absent.

---

## Architecture at a glance

```
                ┌─────────────────────────────────────────────┐
                │                 pyx500r                       │
                ├───────────────┬───────────────┬───────────────┤
  .wiff2  ─────▶│  WiffReader   │               │ crypto.py     │  AES-128-OFB
  .wiff.scan    │  (reader.py)  │               │ (decrypt_*)   │  SQLite SEE
                ├───────────────┤    Bridge     ├───────────────┤
 .qsession ────▶│ QSessionReader│ (bridge.py)   │ tof.py        │  RLE codec
                │ (qsession.py) │               │ centroid.py   │  + calibration
                ├───────────────┴───────────────┴───────────────┤
                │  models.py — frozen dataclasses (the wire API) │
                └───────────────────────────────────────────────┘
```

All public read methods return **frozen dataclasses** from `models.py`. These
are the stable surface you serialize to JSON for a GUI. The byte-format parsers
(`tof`, `crypto`, `rtparts`, `binaryformatter`, `xic_gap`, `lbp_reader`) are
internal but documented for completeness.

---

## WIFF2 reader (`WiffReader`)

A `.wiff2` is an encrypted SQLite database holding acquisition metadata; the
spectra themselves live in a companion `.wiff.scan` binary next to it. Both
files must be present.

### Opening

```python
from pyx500r import open_wiff2, WiffReader

with open_wiff2("acquisition.wiff2") as reader:   # recommended
    ...

reader = WiffReader("acquisition.wiff2")          # manual
try:
    ...
finally:
    reader.close()
```

`open_wiff2(path, password=WIFF2_PASSWORD) -> WiffReader`
`WiffReader(wiff_path, password=WIFF2_PASSWORD)`

Raises `ValueError` if the suffix is not `.wiff2`, `FileNotFoundError` if the
`.wiff2` or its `.wiff.scan` companion is missing.

### Method reference

| Method | Returns | Notes |
|--------|---------|-------|
| `sample_count` *(property)* | `int` | Number of samples in the file |
| `list_samples()` | `list[SampleInfo]` | Ordered by sample index |
| `sample_start_time(sample_index=0)` | `datetime \| None` | Parsed acquisition timestamp |
| `list_instruments(sample_index=0)` | `list[InstrumentInfo]` | All devices |
| `get_ms_instrument(sample_index=0)` | `InstrumentInfo \| None` | First MS device |
| `get_experiments(sample_index=0)` | `list[ExperimentInfo]` | Includes `cycle_count` |
| `get_experiment_tic(sample_index=0, experiment_index=0)` | `Chromatogram` | Total-ion chromatogram |
| `get_sample_tic(sample_index=0)` | `Chromatogram` | Merged across experiments |
| `get_cycle_times(sample_index=0, experiment_index=0)` | `list[float]` | RT per cycle (cached) |
| `get_spectrum(...)` | `SpectrumData` | See below |
| `get_spectrum_metadata(...)` | `SpectrumMetadata` | Cheap — no m/z calibration |
| `iter_spectra(...)` | `Iterator[SpectrumData]` | Auto-prefetches |
| `prefetch_experiment(sample_index=0, experiment_index=None)` | `None` | Bulk-load scan rows |
| `clear_prefetch()` | `None` | Free the prefetch cache |
| `close()` | `None` | Closes the in-memory DB |

### Samples, experiments, instruments

```python
with open_wiff2("acquisition.wiff2") as reader:
    print(reader.sample_count)                       # 1

    for s in reader.list_samples():
        print(s.index, s.name, s.source)             # SampleInfo

    for e in reader.get_experiments(sample_index=0):
        print(e.index, e.scan_type, e.ms_level,      # ExperimentInfo
              e.polarity, e.cycle_count)
        # 0 TOFMS 1 negative 293

    ms = reader.get_ms_instrument()                   # InstrumentInfo | None
    print(ms.model_name, ms.serial_number)            # 'X500 QTOF' 'DM230252004'

    print(reader.sample_start_time())                 # datetime | None
```

### Spectra

```python
spec = reader.get_spectrum(
    sample_index=0,
    experiment_index=0,
    cycle_index=51,
    centroid=False,        # True → centroid before returning
    return_arrays=False,   # True → numpy arrays; False → Python lists
)
# SpectrumData:
spec.mz             # Sequence[float]  (list, or np.ndarray if return_arrays)
spec.intensities    # Sequence[float]  ← note the plural; there is NO .intensity
spec.scan_time      # float (minutes)
spec.centroided     # bool
spec.precursor_mz   # float | None       (MS2 only)
spec.isolation_target_mz / isolation_lower_offset / isolation_upper_offset
```

> ⚠️ The field is **`intensities`** (plural). There is no `.intensity` and no
> `.mz is None` sentinel — an empty spectrum yields a length-0 sequence.

Metadata-only access skips the per-point quadratic m/z calibration and is
markedly cheaper for large profile scans:

```python
meta = reader.get_spectrum_metadata(0, 0, 51)   # SpectrumMetadata
meta.point_count    # 156460 — identical to len(get_spectrum(...).mz)
meta.ms_level       # 1
meta.polarity       # 'negative'
```

### Iterating efficiently

`iter_spectra` auto-prefetches all scan rows for the experiment on first use to
avoid per-cycle SQLite round-trips. For tight loops, request numpy arrays:

```python
for spec in reader.iter_spectra(0, experiment_index=1,
                                limit=None, return_arrays=True):
    if spec.precursor_mz is None:
        continue
    # spec.mz / spec.intensities are np.ndarray here
```

Manual prefetch control (e.g. iterating the same experiment repeatedly):

```python
reader.prefetch_experiment(sample_index=0, experiment_index=1)
... # many get_spectrum() calls hit the in-memory cache
reader.clear_prefetch()
```

---

## QSession reader (`QSessionReader`)

A `.qsession` is an encrypted SQLite quantitation session (MultiQuant). Two
on-disk layouts exist and are auto-detected on open:

- **v1** — results live in a custom `RTParts` binary blob (parsed by
  `rtparts.py` + `binaryformatter.py` + `xic_gap.py`).
- **v2** — results live in regular `MultiSample` / `MultiPeak` / `QualPeak`
  tables.

The public API is identical for both.

### Opening

```python
from pyx500r import open_qsession, QSessionReader

with open_qsession("results.qsession") as qs:
    ...
```

`open_qsession(path, password=QSESSION_PASSWORD) -> QSessionReader`

### Method reference

| Method | Returns | Notes |
|--------|---------|-------|
| `version` *(property)* | `str \| None` | Creating software version |
| `qmap_version` *(property)* | `str \| None` | Quant-method version |
| `locked` *(property)* | `bool` | Read-only session flag |
| `sample_keys` *(property)* | `list[str]` | Distinct XIC sample keys |
| `xic_count` *(property)* | `int` | Rows in `XicRawTable` |
| `list_samples()` / `iter_samples()` | `list/Iterator[QuantSampleInfo]` | |
| `get_sample(index)` | `QuantSampleInfo` | Raises `IndexError` |
| `list_compounds()` / `iter_compounds()` | `list/Iterator[CompoundInfo]` | |
| `iter_peaks(sample_index=None)` | `Iterator[QuantPeakInfo]` | All, or one sample |
| `get_peak(sample_index, compound_index)` | `QuantPeakInfo` | |
| `results_matrix()` | `list[list[QuantPeakInfo]]` | `samples × compounds` |
| `list_xics(sample_key=None)` | `list[XicInfo]` | XIC metadata |
| `get_xic(xic_id)` | `XicChromatogram` | Raises `KeyError` |
| `get_xic_by_mz_range(sample_key, mz_lower, mz_upper)` | `XicChromatogram \| None` | |
| `iter_xics(sample_key=None, return_arrays=False)` | `Iterator[XicChromatogram]` | |
| `get_chromatogram(sample_index, compound_index)` | `XicChromatogram \| None` | Cached |
| `list_audit_records(limit=100)` | `list[dict]` | Audit trail |
| `resolve_library_hits(library_db)` | `int` | Enrich hits with names |
| `search_library(reader, library_db, ...)` | `dict[(int,int), list[LibraryHit]]` | See below |
| `close()` | `None` | |

### Compounds, samples, peaks

```python
with open_qsession("results.qsession") as qs:
    compounds = qs.list_compounds()       # list[CompoundInfo]
    c = compounds[0]
    c.name, c.formula, c.mz_lower, c.mz_upper, c.is_analyte

    samples = qs.list_samples()           # list[QuantSampleInfo]
    s = samples[0]
    s.sample_name, s.sample_id, s.sample_signature

    peak = qs.get_peak(sample_index=0, compound_index=616)   # QuantPeakInfo
    peak.area, peak.retention_time, peak.height, peak.valid_integration

    matrix = qs.results_matrix()          # list[list[QuantPeakInfo | None]]
    #   matrix[si][ci] is None where a sample has no peak for that compound
```

### XIC chromatograms

There are two kinds:

1. **Pre-computed XICs** stored in the session (`XicRawTable`) — fast, no raw
   file needed:

   ```python
   xic = qs.get_chromatogram(sample_index=1, compound_index=616)  # XicChromatogram | None
   xic.times          # np.ndarray (fast frombuffer path)
   xic.intensities    # np.ndarray
   ```

2. **Re-extracted XICs** from the raw `.wiff2` — see `Bridge.extract_xic`.

### Per-peak XIC detail (`xic_result`)

Every `QuantPeakInfo.xic_result` is a dict of the underlying `XicManagerXic`
fields (measured mass/RT, MS/MS linkage, library hits). Prefer accessing these
through [`UnifiedPeak`](#unifiedpeak), which exposes them with clean names.

```python
x = peak.xic_result          # dict | None
x["_foundAtMass"], x["_foundAtRt"], x["_area"], x["_containsMSMS"]
```

---

## Bridge (`WiffQSessionBridge`)

Pairs a `.qsession` with its `.wiff2` acquisition(s) and routes sample lookups
automatically. This is the recommended entry point for a results-viewer GUI.

```python
from pyx500r import WiffQSessionBridge

with WiffQSessionBridge(
    qsession="results.qsession",
    wiffs=["sample_N.wiff2", "sample_P.wiff2"],
    match_by="name",                       # or "position"
    library_db="libview.sqlite",           # optional, enriches library hits
) as bridge:
    if bridge.library_hits_resolved is not None:
        print(f"resolved {bridge.library_hits_resolved} library names")
    ...
```

Constructor: `WiffQSessionBridge(qsession, wiffs, *, match_by="name", library_db=None)`
- `match_by="name"` — pair qsession samples to wiff samples by sample name
  (case-insensitive).
- `match_by="position"` — pair by index order.
- `library_db` — optional `libview_*.sqlite`; if provided and present, library
  hit GUIDs are resolved to names on `open()`. The count is stored in
  `bridge.library_hits_resolved` (`None` if no DB was supplied).

### Methods

| Method | Returns |
|--------|---------|
| `match_samples()` | `list[dict]` — routing: `qsession_index`, `qsession_sample`, `wiff_index`, `wiff_sample` |
| `get_extraction_window(compound_index)` | `ExtractionWindow \| None` |
| `extract_xic(sample_index, compound_index, rt_start=None, rt_end=None)` | `XicChromatogram \| None` |
| `unified_results()` | `list[list[UnifiedPeak]]` |
| `compounds` *(property)* | `list[CompoundInfo]` |
| `samples` *(property)* | `list[QuantSampleInfo]` |
| `peaks` *(property)* | `list[QuantPeakInfo]` |
| `get_peak(si, ci)` | `QuantPeakInfo` |
| `results_matrix()` | `list[list[QuantPeakInfo]]` |
| `get_chromatogram(si, ci)` | `XicChromatogram \| None` |
| `qsession` *(property)* | `QSessionReader` |
| `wiffs` *(property)* | `list[WiffReader]` |
| `library_hits_resolved` *(attr)* | `int \| None` |

```python
window = bridge.get_extraction_window(42)
# ExtractionWindow(period=0, experiment=0, mz_center=310.13, mz_half_window=0.01)

xic = bridge.extract_xic(0, 42, rt_start=5.0, rt_end=7.0)   # re-extracted from raw
matrix = bridge.unified_results()
up = matrix[1][616]                                          # UnifiedPeak
```

`extract_xic` re-extracts a chromatogram from raw scan data by summing every
profile point in `[mz_center ± mz_half_window]` across the cycles whose RT falls
in `[rt_start, rt_end]`. Returns `None` if no extraction window exists, the
sample can't be routed to a wiff file, or no cycles match. Raises `IndexError`
for an out-of-range `compound_index`.

---

## UnifiedPeak

`UnifiedPeak` flattens **compound + peak + XIC** into one read-only namespace —
ideal as a single GUI row model. Construct directly, or get a full matrix from
`bridge.unified_results()`.

```python
from pyx500r import UnifiedPeak
up = UnifiedPeak(peak, compound)        # both may be None-tolerant
```

Selected properties (all are read-only):

```python
# identity / compound
up.name, up.formula, up.mz_lower, up.mz_upper, up.period, up.experiment
up.is_analyte, up.is_reportable, up.internal_std_name

# integration (peak)
up.area, up.corrected_area, up.height, up.corrected_height
up.retention_time, up.start_rt, up.end_rt, up.apex_rt, up.apex_y
up.noise, up.signal_to_noise, up.valid_integration

# measured values (XIC, no underscore prefix)
up.found_mass, up.found_rt, up.found_rt_apex, up.found_rt_start, up.found_rt_end
up.xic_area, up.xic_intensity, up.base_mass, up.extraction_mass, up.extraction_width
up.mass_error          # (found_mass - extraction_mass) / extraction_mass, or None
up.isotope_diff        # isotope ratio diff × 100, or None
up.rt_diff             # found_rt - expected_rt, or None
up.contains_msms, up.msms_experiment, up.msms_retention_time
up.ms1_cycle, up.msms_cycle, up.charge, up.modification_text
up.has_been_calculated, up.has_library_been_searched, up.is_qualifier

# library matches
up.library_hits        # list[LibraryHit]

# raw escape hatches
up.peak                # QuantPeakInfo
up.compound            # CompoundInfo | None
up.xic                 # dict | None  (raw XicManagerXic)

# convenience
up.is_valid()          # valid_integration AND has_been_calculated
```

Any property returns a safe default (`0.0`, `""`, `None`, `[]`) when the backing
peak/compound/xic is missing, so it never raises `AttributeError` on partial
data — convenient for serializing sparse matrices.

---

## Precursor index

For batch MS2 precursor lookup across many files without decoding spectra.
Requires `numpy`.

```python
from pyx500r import build_precursor_index, PrecursorIndex, CrossFilePrecursorIndex

idx = build_precursor_index(Path("a.wiff2"))   # PrecursorIndex
hits = idx.find(target_mz=500.0, tolerance_da=0.01)   # list[(sample, exp, cycle)]

# Aggregate across many files:
cross = CrossFilePrecursorIndex()
for p in wiff_paths:
    cross.add(build_precursor_index(p))
results = cross.find(500.0, precursor_tolerance_da=0.01)
# list[dict]: file_path, sample_index, experiment_index, cycle_index,
#             precursor_mz, retention_time  (sorted by file_path, precursor_mz)
cross.save(Path("index.json"))
cross = CrossFilePrecursorIndex.load(Path("index.json"))
```

Only MS2 scans are indexed, so `find(..., ms_level != 2)` returns `[]`.

---

## Library search

`QSessionReader.search_library` matches each session MS/MS peak against an
external LibraryView SQLite DB by extracting the linked spectrum from a wiff
reader and dot-product scoring.

```python
results = qs.search_library(
    reader,                      # a WiffReader covering the session's wiff files
    "libview.sqlite",
    ppm_tol=50.0,
    dot_product_ppm=20.0,
    prescreen_n=200,
    top_n=5,
    polarity_override=None,      # "POS"/"NEG" to override experiment metadata
)
# dict[(sample_index, compound_index)] -> list[LibraryHit]
```

The lower-level `pyx500r.libsearch.LibrarySearcher` can be used directly to
search arbitrary peak lists; see `pyx500r/libsearch.py` and the
`pyx500r-libsearch` CLI.

---

## TOF codec & calibration

```python
from pyx500r import decompress_tof, TofCalibration

# Decode a compressed TOF stream. With cal_* params the m/z axis is computed
# inside the (optionally numba-JIT) kernel.
mz, intensities = decompress_tof(
    stream,                              # bytes, starting at the FF FF FF FF sentinel
    number_of_time_bins_to_sum=1,
    min_bin=0,
    cal_a=slope, cal_t0=delay, time_resolution=tdc_res,
    return_arrays=True,                  # numpy arrays vs Python lists
)

cal = TofCalibration(cal_a=slope, cal_t0=delay, time_resolution=tdc_res)
cal.bin_to_mass(time_bin)        # float
cal.bins_to_masses([...])        # list[float]
cal.mass_to_bin(mass)            # float
```

> ⚠️ `TofCalibration` is constructed with explicit `cal_a`, `cal_t0`,
> `time_resolution` — there is **no** `TofCalibration.from_database(...)`.
> Inside `WiffReader` the calibration is read per-scan from the `scanItems`
> columns `slope`, `delay`, `tdcResolution`.

Also available: `compress_tof`, `decompress_quad`, `decompress_zero_width`,
`mass_to_time`, `time_to_mass`, `MassRange`.

---

## Centroiding

```python
from pyx500r import centroid_spectrum

cmz, cints = centroid_spectrum(
    mz, intensities,
    centroid_percentage=50.0,
    return_arrays=True,
)
```

`centroid_spectrum` dispatches to a numba kernel when available, else a
NumPy/pure-Python fallback that is bit-compatible. Helpers `add_framing_zeros`,
`moving_average_smooth`, and the `Peak` dataclass are also exported.

---

## Decryption primitives

For tooling that needs the raw SQLite bytes:

```python
from pyx500r import (
    decrypt_database, decrypt_page,
    WIFF2_PASSWORD, QSESSION_PASSWORD,
    PAGE_SIZE, QSESSION_PAGE_SIZE, RESERVED_BYTES,
)

raw_sqlite = decrypt_database("a.wiff2", WIFF2_PASSWORD)            # bytes
raw_sqlite = decrypt_database("r.qsession", QSESSION_PASSWORD,
                              page_size=QSESSION_PAGE_SIZE)
```

---

## Data model reference

All read methods return these frozen dataclasses (`pyx500r.models`). `slots=True`
and `frozen=True`, so they are immutable and serializable with
`dataclasses.asdict`.

### SampleInfo
| Field | Type |
|-------|------|
| `index` | `int` |
| `sample_id` | `str` |
| `name` | `str` |
| `source` | `str` |
| `start_timestamp` | `str \| None` |

### ExperimentInfo
| Field | Type |
|-------|------|
| `index` | `int` |
| `experiment_id` | `str` |
| `scan_type` | `str` |
| `ms_level` | `int` |
| `polarity` | `str` (`"positive"`/`"negative"`/`"unknown"`) |
| `cycle_count` | `int` |

### InstrumentInfo
| Field | Type |
|-------|------|
| `sample_index` | `int` |
| `instrument_index` | `int` |
| `device_type` | `int \| None` |
| `device_name` | `str \| None` |
| `model_name` | `str \| None` |
| `serial_number` | `str \| None` |
| `is_mass_spectrometer` | `bool` |

### SpectrumMetadata
| Field | Type |
|-------|------|
| `sample_index`, `experiment_index`, `cycle_index` | `int` |
| `scan_time` | `float` |
| `scan_type` | `str` |
| `ms_level` | `int` |
| `polarity` | `str` |
| `point_count` | `int` |
| `precursor_mz` | `float \| None` |
| `isolation_target_mz` / `isolation_lower_offset` / `isolation_upper_offset` | `float \| None` |

### SpectrumData
| Field | Type |
|-------|------|
| `sample_index`, `experiment_index`, `cycle_index` | `int` |
| `scan_time` | `float` |
| `mz` | `Sequence[float]` — list or `np.ndarray` |
| `intensities` | `Sequence[float]` — **plural** |
| `centroided` | `bool` |
| `precursor_mz` | `float \| None` |
| `isolation_target_mz` / `isolation_lower_offset` / `isolation_upper_offset` | `float \| None` |

### Chromatogram
| Field | Type |
|-------|------|
| `times` | `list[float]` |
| `intensities` | `list[float]` |
| `experiment_index` | `int \| None` |
| `ms_level` | `int \| None` |

### XicInfo
| Field | Type |
|-------|------|
| `xic_id` | `str` |
| `sample_key` | `str` |
| `mz_lower` / `mz_upper` | `float` |
| `group_index` / `replicate_index` | `int` |
| `status` | `int` |

### XicChromatogram
| Field | Type |
|-------|------|
| `xic_id` | `str` |
| `sample_key` | `str` |
| `times` | `Sequence[float]` |
| `intensities` | `Sequence[float]` |
| `mz_lower` / `mz_upper` | `float` |
| `status` | `int` |

### CompoundInfo (selected; 32 fields total)
| Field | Type | Description |
|-------|------|-------------|
| `name` | `str` | |
| `group_name` | `str` | |
| `formula` | `str \| None` | |
| `charge_formula` / `adduct_formula` | `str \| None` | |
| `precursor_mass` / `fragment_mass` | `float` | |
| `extraction_type` | `int` | 0=signal, 1=MS/MS |
| `period` / `experiment` | `int` | |
| `mz_lower` / `mz_upper` | `float` | extraction window |
| `is_analyte` / `is_reportable` / `is_non_targeted` / `is_summed` / `is_from_multi_period_data` | `bool` | |
| `isotope_index` | `int` | |
| `expected_mw` | `float` | |
| `units` | `str` | |
| `comment` / `internal_std_name` | `str \| None` | |
| `regression_area` | `bool \| None` | |
| `regression_type` / `regression_weighting` | `int \| None` | |
| `use_auto_regression` | `bool \| None` | |
| `integration_parameters` | `dict \| None` | |
| `acquisition_indices` / `summed_compounds` | `list[int] \| None` | |
| `extraction_values1` / `extraction_values2` | `list[float] \| None` | raw m/z windows |

### QuantSampleInfo (selected; 24 fields total)
| Field | Type |
|-------|------|
| `index` | `int` |
| `sample_name` / `sample_id` | `str` |
| `sample_type` | `int` |
| `sample_comment` | `str \| None` |
| `dilution_factor` / `injection_volume` | `float` |
| `user_name` / `acq_method_name` / `instrument_name` / `instrument_serial_number` | `str \| None` |
| `batch_name` / `barcode` / `scanned_barcode` | `str \| None` |
| `autosampler_method_supports_barcode` / `sample_comparison` | `bool` |
| `ms_method` / `lc_method` / `sample_signature` | `str \| None` |
| `rack` / `plate` / `vial` | `str \| None` |
| `acquisition_date` | `datetime \| None` |

### QuantPeakInfo (selected; 48 fields total)
| Field | Type | Description |
|-------|------|-------------|
| `sample_index` / `compound_index` / `peak_index` | `int` | |
| `area` / `corrected_area` / `original_area` / `region_area` | `float` | |
| `height` / `corrected_height` / `region_height` | `float` | |
| `retention_time` | `float` | |
| `start_rt` / `start_y` / `end_rt` / `end_y` | `float` | peak bounds |
| `apex_rt` / `apex_y` | `float` | |
| `half_height_start_rt` / `half_height_end_rt` | `float` | |
| `noise` / `signal_to_noise` | `float` | |
| `valid_integration` / `reportable` / `failed_query` / `modified` / `row_hidden` | `bool` | |
| `use_for_calibration` | `bool` | |
| `actual_concentration` / `std_addn_actual_concentration` | `float` | |
| `profile_type` / `peak_type` | `int` | |
| `molecular_weight` | `float` | |
| `points_across_baseline` / `points_across_half_height` | `int` | |
| `integration_parameters` | `dict \| None` | |
| `profile` | `list[float] \| None` | raw chromatogram |
| `custom_fields` / `custom_peak_fields` | `dict \| None` | |
| `xic_result` | `dict \| None` | full XicManagerXic |
| `super_group_id` | `str \| None` | |

### LibraryHit
| Field | Type | Source |
|-------|------|--------|
| `fit` / `reverse_fit` / `purity` | `float` | built-in MultiQuant search |
| `is_smart_confirmation` | `bool` | |
| `library_entry_id` | `str` | resolved GUID |
| `name` / `formula` / `cas` | `str` | from external library |
| `score` | `float` | external dot-product |
| `precursor_mz` / `collision_energy` | `float` | |
| `num_peaks` | `int` | |
| `spectrum_id` / `compound_id` | `str` | |

---

## Serialization notes (for web APIs)

- **Dataclasses → JSON**: every model is a frozen dataclass.
  `dataclasses.asdict(obj)` yields a JSON-native dict (FastAPI/Pydantic accept
  these directly). `acquisition_date` is the only `datetime` field — encode as
  ISO-8601.
- **`mz`/`intensities`/`times`** are Python `list[float]` unless you pass
  `return_arrays=True`, in which case they are `np.ndarray`. **Call
  `.tolist()` before JSON encoding numpy arrays**, or simply omit
  `return_arrays` when feeding a serializer.
- **`xic_result`** dicts contain only JSON-native scalar/list/dict values
  (verified on real data) — safe to pass through, though the keys use .NET names
  (`_foundAtMass`, etc.). Prefer projecting through `UnifiedPeak` for clean keys.
- **`UnifiedPeak` is not a dataclass** — it's a property-backed view. Build an
  explicit dict for the wire (see `examples/server/` for a `unified_to_dict`
  helper).
- **Large arrays**: a dense TOF profile is ~156k points (~2.5 MB as JSON).
  Prefer centroided spectra, server-side downsampling, or a binary endpoint for
  raw profiles (see [`GUI_INTEGRATION.md`](./GUI_INTEGRATION.md)).

---

## Concurrency & lifecycle

- Each reader wraps an **in-memory SQLite connection** created with the default
  `check_same_thread=True`. **Do not share a single reader across threads.**
- For a multi-worker server, open one reader **per request** (cheap: ~5 ms to
  decrypt+open the example file) or maintain a **per-file pool keyed by path**
  with a lock, or run CPU-bound decode work in a process pool. See
  [`GUI_INTEGRATION.md`](./GUI_INTEGRATION.md#concurrency).
- Always `close()` (or use the context manager) to release the in-memory DB and
  the memory-mapped `.wiff.scan` bytes.
- Caches (`prefetch`, cycle-times, experiment-info, chromatogram) live on the
  reader instance and are freed on `close()`.

---

## Performance

Measured on the bundled 156k-point reference acquisition (`numba` installed):

| Operation | Time | Notes |
|-----------|------|-------|
| `open_wiff2` (decrypt + in-memory SQLite) | ~5 ms | per file |
| `get_spectrum` (dense, 156k points) | ~6 ms | JIT TOF decode + calibration |
| `get_spectrum_metadata` (point count only) | ~3 ms | skips m/z calibration |
| Open `.qsession` + full load (2,471 cmpds / 4,942 peaks) | ~4 s | one-time; cached after |
| Parse XIC gap (4,942 blobs) | ~0.7 s | via `nrbf` |

Without `numba`, TOF decode and centroiding are ~5–10× slower but byte-identical.

---

## CLI tools

Installed as console scripts (see `pyproject.toml [project.scripts]`):

| Command | Module | Purpose |
|---------|--------|---------|
| `pyx500r` | `pyx500r.cli` | List samples / find precursor→product transitions |
| `ppyx500r` | `pyx500r.cli_parallel` | Same, multiprocess (`-j N`) |
| `pyx500r-index` | `pyx500r.index_builder` | Build an `.npz` MS2 product-ion index |
| `pyx500r-search` | `pyx500r.w2searcher` | Query an `.npz` index (TUI + one-shot) |
| `pyx500r-qsession` | `pyx500r.cli_bridge` | Interactive qsession results explorer |
| `pyx500r-libsearch` | `pyx500r.libsearch_cli` | Search MGF/peaks against a LibraryView DB |

```bash
pyx500r list acquisition.wiff2
pyx500r transitions "*.wiff2" --precursor-mz 456.2 --tolerance-ppm 20
pyx500r transitions file.wiff2 -t "250.1587:191.0857,163.0907,109.0443" --json
pyx500r-qsession results.qsession ./wiff_dir --library-db libview.sqlite
```

Every CLI `main(argv=None) -> int` is importable and testable.

---

## Exceptions

| Raised by | Exception | When |
|-----------|-----------|------|
| `WiffReader.__init__` | `ValueError` | non-`.wiff2` suffix |
| `WiffReader.__init__` | `FileNotFoundError` | missing `.wiff2` or `.wiff.scan` |
| `WiffReader.get_spectrum` / `_scan_item` | `IndexError` | sample/experiment/cycle out of range |
| `QSessionReader.__init__` | `ValueError` | non-`.qsession` suffix |
| `QSessionReader.__init__` | `OSError` | cannot decrypt with any known page size |
| `QSessionReader.get_xic` | `KeyError` | unknown XIC id |
| `QSessionReader.get_peak` / `get_sample` | `IndexError` | out-of-range indices |
| `WiffQSessionBridge.__init__` | `ValueError` | bad `match_by` or wrong file suffix |
| `WiffQSessionBridge.extract_xic` | `IndexError` | out-of-range `compound_index` |
| `decrypt_database` | `ValueError` | file size not a multiple of page size |

---

## See also

- [`GUI_INTEGRATION.md`](./GUI_INTEGRATION.md) — FastAPI/Uvicorn + TypeScript guide
- [`types.ts`](./types.ts) — TypeScript interfaces for all models
- [`examples/server/`](../examples/server/) — runnable reference FastAPI app
- [`QSESSION_DATA_MODEL.md`](./QSESSION_DATA_MODEL.md) — qsession internals
- [`LBP_FORMAT.md`](./LBP_FORMAT.md) — LibraryView `.lbp` format
