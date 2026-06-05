# pyx500r

Pure-Python reader for SCIEX X500R QTOF acquisition files (`.wiff2` + `.wiff.scan`) and MultiQuant qsession result tables (`.qsession`).

Built for **small-molecule screening**, **toxicology** and **forensic toxicology** workflows where you need programmatic access to:

- Full-scan TOF spectra and MS/MS precursor data from X500R acquisitions
- Extracted ion chromatograms (XICs) and quantitation peak results from MultiQuant sessions
- Batch precursor-indexed search across large file cohorts

Works on macOS, Linux and Windows.

## Installation

```bash
# Install the package
pip install -e .

# Development
pip install -e .[dev]
```

## Quick start — X500R acquisition files

```python
from pyx500r import open_wiff2

with open_wiff2("toxicology_run.wiff2") as reader:
    print(reader.list_samples())
    print(reader.get_experiments())

    # Total ion chromatogram for experiment 0
    tic = reader.get_experiment_tic(experiment_index=0)

    # Spectrum at cycle 51 (TOF-MS or MS/MS)
    spectrum = reader.get_spectrum(experiment_index=0, cycle_index=51)
    print(f"{len(spectrum.mz)} points, m/z[0] = {spectrum.mz[0]:.7f}")
```

## Quick start — MultiQuant qsession result tables

```python
from pyx500r import open_qsession

with open_qsession("quant_results.qsession") as qs:
    # Iterate extracted-ion chromatograms
    for xic in qs.iter_xics():
        print(f"{xic.xic_id}: {len(xic.times)} points")

    # Pull the full quantitation results matrix
    matrix = qs.results_matrix()
    print(f"{len(matrix)} samples × {len(matrix[0])} compounds")
```

## Reader paths

| Install   | Reader class   | Requirements         | OS  |
|-----------|----------------|----------------------|-----|
| `pyx500r` | `WiffReader`   | `cryptography`, `numba`, `numpy` | Any |

`numba` is included by default, so there is no separate slow-path install.

## Running tests

```bash
python -m pytest tests/ -v
```

Test data should be placed in `data/encoded/` (`.wiff2` + `.wiff.scan`) and `data/qsession/` (`.qsession`).

## License

GPL-3.0-or-later
