# pyx500r

Pure-Python reader for SCIEX X500R QTOF acquisition files (`.wiff2` + `.wiff.scan`)
and MultiQuant qsession result tables (`.qsession`) — **no .NET / SCIEX DLL
dependency**.

Built for **small-molecule screening**, **toxicology** and **forensic
toxicology** workflows where you need programmatic access to:

- Full-scan TOF spectra and MS/MS precursor data from X500R acquisitions
- Extracted ion chromatograms (XICs) and quantitation peak results from MultiQuant sessions
- A `WiffQSessionBridge` that pairs acquisitions with quant results for unified compound × sample access
- Batch precursor-indexed search across large file cohorts
- LibraryView (`.lbp` / `.sqlite`) reference-spectrum reading and library search

Works on macOS, Linux and Windows.

## Installation

```bash
pip install -e .            # core reader (cryptography + numba + numpy + nrbf)
pip install -e .[cli]       # + tqdm progress bars
pip install -e .[server]    # + FastAPI/Uvicorn for the example web apps
pip install -e .[dev]       # + pytest
```

## Quick start — X500R acquisition files

```python
from pyx500r import open_wiff2

with open_wiff2("toxicology_run.wiff2") as reader:
    print(reader.list_samples())
    print(reader.get_experiments())

    tic = reader.get_experiment_tic(experiment_index=0)
    spectrum = reader.get_spectrum(experiment_index=0, cycle_index=51)
    print(f"{len(spectrum.mz)} points, m/z[0] = {spectrum.mz[0]:.7f}")
```

## Quick start — MultiQuant qsession result tables

```python
from pyx500r import open_qsession

with open_qsession("quant_results.qsession") as qs:
    print(qs.version, len(qs.list_compounds()), "compounds")
    matrix = qs.results_matrix()                 # samples × compounds
    peak = qs.get_peak(sample_index=0, compound_index=0)
    print(peak.area, peak.retention_time)
```

## Quick start — bridge acquisitions + results

```python
from pyx500r import WiffQSessionBridge

with WiffQSessionBridge("quant_results.qsession",
                        ["sample_N.wiff2", "sample_P.wiff2"],
                        match_by="name") as bridge:
    for row in bridge.unified_results():
        for up in row:
            if up.is_valid():
                print(up.name, up.area, up.found_mass, up.contains_msms)
```

## Documentation

| Document | What it covers |
|----------|----------------|
| [`docs/API.md`](docs/API.md) | Full API reference — readers, bridge, `UnifiedPeak`, TOF codec, every dataclass, exceptions, performance |
| [`docs/GUI_INTEGRATION.md`](docs/GUI_INTEGRATION.md) | Building a **FastAPI/Uvicorn + TypeScript** GUI: serialization, concurrency, large-spectrum strategy, endpoint design |
| [`docs/types.ts`](docs/types.ts) | TypeScript interfaces mirroring every model |
| [`examples/server/`](examples/server/) | Reference FastAPI JSON API over `.wiff2` / `.qsession` |
| [`examples/library_browser/`](examples/library_browser/) | GUI app to browse a LibraryView `.sqlite`: compounds, settings, reference-spectrum plots |
| [`docs/QSESSION_DATA_MODEL.md`](docs/QSESSION_DATA_MODEL.md) · [`docs/QSESSION_FORMAT.md`](docs/QSESSION_FORMAT.md) | qsession internals & container format |
| [`docs/LBP_FORMAT.md`](docs/LBP_FORMAT.md) · [`docs/RTPARTS_COMPOUND_PARSING.md`](docs/RTPARTS_COMPOUND_PARSING.md) · [`docs/BLOB_GAP_PARSING.md`](docs/BLOB_GAP_PARSING.md) | Binary format internals |

## Command-line tools

| Command | Purpose |
|---------|---------|
| `pyx500r` | List samples / find precursor→product transitions |
| `ppyx500r` | Same, multiprocess (`-j N`) |
| `pyx500r-index` | Build an `.npz` MS2 product-ion index |
| `pyx500r-search` | Query an `.npz` index (TUI + one-shot) |
| `pyx500r-qsession` | Interactive qsession results explorer |
| `pyx500r-libsearch` | Search MGF/peaks against a LibraryView DB |

## Running tests

```bash
python -m pytest -q
```

The committed test suite is **data-free** and runs without any mass-spectrometry
files. Tests that exercise real `.wiff2` / `.qsession` / LibraryView data are not
part of this repository — no MS data is shipped here. If you have your own files,
place them under a local (git-ignored) `data/` directory; the example apps and a
local test harness can then point at them via environment variables
(`PYX500R_DATA_ROOT`, `PYX500R_LIBRARY_DB`).

## License

GPL-3.0-or-later
