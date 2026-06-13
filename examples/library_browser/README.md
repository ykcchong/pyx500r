# LibraryView Browser

An interactive browser for a SCIEX **LibraryView** `.sqlite` database (the kind
under `data/`, e.g. `libview.sqlite`). Search compounds, inspect their
settings/metadata, and view their **reference spectra** — including compounds
that have *multiple* spectra (different polarities, collision energies,
instruments) — rendered as interactive plots in the browser.

```
┌─────────────┬────────────────────┬──────────────────────────┐
│ compound    │ compound settings  │ reference spectra          │
│ search list │ (names, formula,   │  [POS CE35] [NEG CE20] …   │
│             │  thresholds, libs, │  ┌──────────────────────┐  │
│ Caffeine    │  retention times)  │  │  ▍ ▍▍  ▍ ▍   spectrum │  │
│ Cocaine     │                    │  │ ▍▍▍▍▍▍▍▍▍▍▍   plot     │  │
│ …           │                    │  └──────────────────────┘  │
└─────────────┴────────────────────┴──────────────────────────┘
```

## Run

```bash
# from the repo root
pip install -e ".[numba]"
pip install -r examples/library_browser/requirements.txt

export PYX500R_LIBRARY_DB=data/libview.sqlite
uvicorn examples.library_browser.app:app --reload --port 8001
```

Open **http://localhost:8001/** for the GUI, or **/docs** for the Swagger API.

## Configuration

| Variable | Default | Meaning |
|----------|---------|---------|
| `PYX500R_LIBRARY_DB` | `data/libview.sqlite` | Path to the LibraryView `.sqlite` file (opened read-only) |

## Features

- **Search & filter** compounds by name / formula / CAS / identifier, and by
  library (collection). Server-side paginated, so the 7k-compound DB stays snappy.
- **Compound settings** panel: identifier, formula, MW / monoisotopic mass,
  purity / yellow-flag / red-flag thresholds, all alternative names (with
  regions), library memberships, and retention times per instrument.
- **Reference spectra**: each compound's spectra are shown as selectable tabs
  (polarity, collision energy, instrument, precursor m/z). Selecting one renders
  it as a plot.
  - **Centroid** view → labelled stick plot (top peaks annotated with m/z).
  - **Raw / profile** view → continuous trace.
  - Spectra whose XY blobs are stored under an extra encryption layer
    (`Encryption` like `HRAIO|2.0`) are flagged and dimmed — the reader can't
    decode those.

## Architecture

| File | Role |
|------|------|
| `library_db.py` | Read-only data-access layer over the LibraryView schema. Thread-local SQLite connections; decodes XY blobs via `pyx500r.libsearch._extract_double_arrays_from_blob`. Importable/usable on its own. |
| `app.py` | FastAPI backend + static file serving. |
| `static/index.html`, `styles.css`, `app.js` | Zero-dependency single-page GUI (vanilla JS + Canvas — no build step). |

## API

| Endpoint | Returns |
|----------|---------|
| `GET /api/stats` | compound / spectrum / library counts |
| `GET /api/libraries` | libraries with compound counts |
| `GET /api/compounds?search=&library_id=&offset=&limit=` | paged compound list |
| `GET /api/compounds/{id}` | full compound settings + spectrum metadata |
| `GET /api/compounds/{id}/spectra` | spectrum metadata list only |
| `GET /api/spectra/{id}?kind=centroid\|raw` | spectrum metadata + decoded `mz[]`/`intensity[]` |

## Using the data layer directly

```python
from examples.library_browser.library_db import LibraryDB

db = LibraryDB("data/libview.sqlite")
hit = db.list_compounds(search="caffeine")["items"][0]
detail = db.get_compound(hit["id"])
for s in detail["spectra"]:
    spec = db.get_spectrum(s["id"], kind="centroid")
    print(s["polarity"], s["collision_energy"], "→", spec["num_peaks"], "peaks")
```

## Notes

- The database is opened **read-only** (`mode=ro`); the browser never writes.
- Decode is CPU-bound and runs in a thread pool so the event loop stays free.
- For a production deployment add auth and a response cache (spectra are
  immutable). This example favours clarity over completeness.
