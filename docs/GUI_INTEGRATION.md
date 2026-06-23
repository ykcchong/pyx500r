# Building a Web GUI on pyx500r (FastAPI / Uvicorn + TypeScript)

This guide shows how to expose `pyx500r` as a JSON HTTP API suitable for a
single-page TypeScript front-end (React/Vue/Svelte/etc.) that browses
`.wiff2` acquisitions and `.qsession` quantitation results.

A complete, runnable reference implementation lives in
[`examples/server/`](../examples/server/). This document explains the *why*
behind it and the patterns you should follow when extending it.

> Read [`API.md`](./API.md) first for the underlying Python API, and pair this
> guide with [`types.ts`](./types.ts) for the TypeScript model definitions.

---

## 1. Design overview

```
 ┌──────────────┐     HTTP/JSON      ┌─────────────────────┐     pyx500r     ┌──────────────┐
 │  TypeScript  │ ◀───────────────▶ │  FastAPI (Uvicorn)   │ ──────────────▶ │  .wiff2 /    │
 │  SPA (React) │   /api/...         │  + ReaderPool        │   dataclasses   │  .qsession   │
 └──────────────┘                    └─────────────────────┘                 └──────────────┘
```

Key principles:

1. **The dataclasses in `models.py` are your wire contract.** They are frozen,
   slotted, and JSON-native (except `datetime`). Serialize them with
   `dataclasses.asdict` and you get stable, typed payloads that map 1:1 onto the
   interfaces in [`types.ts`](./types.ts).
2. **Readers are stateful and single-threaded.** Wrap them in a pool keyed by
   file path; never share one across threads (§4).
3. **Decode work is CPU-bound.** Offload large-spectrum and bulk operations off
   the event loop with `run_in_threadpool` / a process pool (§5).
4. **Spectra are big.** A dense profile scan is ~156k points. Centroid,
   downsample, or use a binary endpoint for raw data (§6).

---

## 2. Project layout

```
examples/server/
├── README.md
├── requirements.txt        # fastapi, uvicorn[standard]
├── app.py                  # FastAPI app + endpoints
├── pool.py                 # path-keyed reader pool (thread-safe)
└── serializers.py          # dataclass / UnifiedPeak → JSON-native dict
```

Install and run:

```bash
pip install -e ".[server]"
uvicorn examples.server.app:app --reload --port 8000
# Interactive docs at http://localhost:8000/docs
```

---

## 3. Serialization

### 3.1 Dataclasses

Every reader return value is a frozen dataclass. Convert to a dict for the wire:

```python
import dataclasses
from pyx500r import open_wiff2

with open_wiff2(path) as r:
    payload = [dataclasses.asdict(e) for e in r.get_experiments()]
# payload is JSON-native; FastAPI returns it directly
```

FastAPI can also return the dataclass instances directly (it knows how to encode
them), but going through `asdict` gives you a place to fix up the two gotchas
below.

### 3.2 The two gotchas

1. **numpy arrays.** If you pass `return_arrays=True`, `mz`/`intensities`/`times`
   become `np.ndarray`, which is **not** JSON-serializable by the stdlib encoder.
   Either omit `return_arrays` (you get lists), or call `.tolist()`:

   ```python
   mz = spec.mz.tolist() if hasattr(spec.mz, "tolist") else list(spec.mz)
   ```

2. **`datetime`.** `QuantSampleInfo.acquisition_date` is the only `datetime`
   field. Encode as ISO-8601:

   ```python
   d = s.acquisition_date
   out["acquisition_date"] = d.isoformat() if d else None
   ```

The reference `serializers.py` handles both centrally.

### 3.3 UnifiedPeak

`UnifiedPeak` is a property-backed view, **not** a dataclass, so `asdict` won't
work. Build an explicit projection (this also defines your front-end's "results
row" shape):

```python
def unified_to_dict(up) -> dict:
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
        "mass_error_ppm": (up.mass_error * 1e6) if up.mass_error is not None else None,
        "contains_msms": up.contains_msms,
        "valid_integration": up.valid_integration,
        "is_valid": up.is_valid(),
        "library_hits": [dataclasses.asdict(h) for h in up.library_hits],
    }
```

### 3.4 `xic_result` raw keys

`QuantPeakInfo.xic_result` is a dict using .NET names (`_foundAtMass`, …). Its
values are JSON-native. Expose it raw only for "advanced/debug" views; for the
main UI, project through `UnifiedPeak` for clean camel/snake keys.

---

## 4. Concurrency: the reader pool

`WiffReader` / `QSessionReader` wrap an in-memory SQLite connection created with
the default `check_same_thread=True`. **A connection must be used from one
thread at a time.** Three viable strategies:

| Strategy | When | Tradeoff |
|----------|------|----------|
| **Open per request** | Low traffic, many distinct files | Simple; ~5 ms open cost each time |
| **Path-keyed pool + lock** | Repeated access to the same files | Fast; must serialize access per file |
| **Process pool** | Heavy CPU decode under load | True parallelism; serialization overhead |

The reference implementation uses a **path-keyed pool guarded by per-file
locks**, with LRU eviction. Sketch:

```python
import threading
from pathlib import Path
from collections import OrderedDict
from pyx500r import open_wiff2, open_qsession

class ReaderPool:
    def __init__(self, max_open: int = 16):
        self._lock = threading.Lock()
        self._wiff: OrderedDict[str, tuple] = OrderedDict()   # path -> (reader, rlock)
        self._qs: OrderedDict[str, tuple] = OrderedDict()
        self._max = max_open

    def wiff(self, path: str):
        key = str(Path(path).resolve())
        with self._lock:
            if key not in self._wiff:
                self._wiff[key] = (open_wiff2(key), threading.Lock())
                self._evict(self._wiff)
            self._wiff.move_to_end(key)
            return self._wiff[key]   # (reader, rlock)

    def _evict(self, od):
        while len(od) > self._max:
            _, (reader, _) = od.popitem(last=False)
            reader.close()
```

Endpoints then do:

```python
reader, rlock = pool.wiff(path)
with rlock:                       # serialize access to this reader
    spec = reader.get_spectrum(...)
```

Because the work inside `with rlock` is synchronous + CPU-bound, run it off the
event loop (next section) so one slow decode doesn't block all requests.

> **Uvicorn workers:** with `--workers N` each worker process gets its own
> `ReaderPool` (no shared memory). That's fine — the pool is a per-process cache.
> Just size `max_open` with total memory in mind (each open file holds its
> decrypted DB + memory-mapped `.wiff.scan`).

---

## 5. Keeping the event loop responsive

`async def` endpoints run on the event loop; a 6 ms decode is fine, but bulk
operations (a whole experiment, a 4 s qsession load, building an index) will
stall every concurrent request. Offload them:

```python
from fastapi.concurrency import run_in_threadpool

@app.get("/api/wiff/{file}/spectrum")
async def spectrum(file: str, exp: int, cycle: int, centroid: bool = False):
    reader, rlock = pool.wiff(resolve(file))

    def work():
        with rlock:
            return reader.get_spectrum(0, exp, cycle, centroid=centroid)

    spec = await run_in_threadpool(work)
    return serialize_spectrum(spec)
```

For genuinely parallel CPU work across files (e.g. batch transitions over
hundreds of files), use a `ProcessPoolExecutor` — `pyx500r`'s own
`cli_parallel` and `index_builder` already do this with `multiprocessing.Pool`.

---

## 6. Spectra are large — strategies

A dense TOF profile scan is **~156,000 points** → ~2.5 MB of JSON. Options, in
order of preference for a GUI:

1. **Centroid server-side** (`centroid=True`) — turns 156k profile points into a
   few thousand peaks. Best default for MS display.
2. **Downsample for overview** — for a zoomed-out trace, bin to the pixel width
   of the chart (e.g. ≤4000 points) server-side; fetch full resolution only on
   zoom.
3. **m/z window slicing** — accept `mz_min`/`mz_max` query params and slice
   before serializing.
4. **Binary endpoint** — for raw profiles, return `Float32Array`-friendly bytes
   instead of JSON (much smaller + faster to parse in TS):

   ```python
   from fastapi import Response
   import numpy as np

   @app.get("/api/wiff/{file}/spectrum.bin")
   async def spectrum_bin(file: str, exp: int, cycle: int):
       reader, rlock = pool.wiff(resolve(file))
       def work():
           with rlock:
               return reader.get_spectrum(0, exp, cycle, return_arrays=True)
       s = await run_in_threadpool(work)
       mz = np.ascontiguousarray(s.mz, dtype=np.float32)
       it = np.ascontiguousarray(s.intensities, dtype=np.float32)
       n = len(mz).to_bytes(4, "little")
       return Response(n + mz.tobytes() + it.tobytes(),
                       media_type="application/octet-stream")
   ```

   TypeScript decode:
   ```ts
   const buf = await (await fetch(url)).arrayBuffer();
   const n = new DataView(buf).getUint32(0, true);
   const mz = new Float32Array(buf, 4, n);
   const intensity = new Float32Array(buf, 4 + n * 4, n);
   ```

5. **Always serve metadata cheaply** via `get_spectrum_metadata` (`point_count`,
   `precursor_mz`, …) so the UI can decide before fetching a heavy payload.

---

## 7. Recommended endpoint map

A pragmatic REST surface for a viewer. (The reference app implements the
starred ★ ones.)

### WIFF2

| Method & path | Backed by | Returns |
|---------------|-----------|---------|
| ★ `GET /api/wiff/{file}/samples` | `list_samples` | `SampleInfo[]` |
| ★ `GET /api/wiff/{file}/experiments?sample=0` | `get_experiments` | `ExperimentInfo[]` |
| `GET /api/wiff/{file}/instruments?sample=0` | `list_instruments` | `InstrumentInfo[]` |
| ★ `GET /api/wiff/{file}/tic?sample=0&experiment=0` | `get_experiment_tic` | `Chromatogram` |
| ★ `GET /api/wiff/{file}/spectrum?...&centroid=` | `get_spectrum` | `SpectrumData` |
| `GET /api/wiff/{file}/spectrum/meta?...` | `get_spectrum_metadata` | `SpectrumMetadata` |
| `GET /api/wiff/{file}/spectrum.bin?...` | `get_spectrum` | `application/octet-stream` |

### QSession

| Method & path | Backed by | Returns |
|---------------|-----------|---------|
| ★ `GET /api/qsession/{file}/info` | `version`/`locked`/counts | summary object |
| ★ `GET /api/qsession/{file}/compounds` | `list_compounds` | `CompoundInfo[]` |
| ★ `GET /api/qsession/{file}/samples` | `list_samples` | `QuantSampleInfo[]` |
| ★ `GET /api/qsession/{file}/peak?sample=&compound=` | `get_peak` | `QuantPeakInfo` |
| `GET /api/qsession/{file}/chromatogram?sample=&compound=` | `get_chromatogram` | `XicChromatogram` |
| ★ `GET /api/qsession/{file}/results?page=&size=` | `results_matrix` | paginated `UnifiedPeak`-like rows |

### Bridge

| Method & path | Backed by | Returns |
|---------------|-----------|---------|
| `GET /api/bridge/results?qs=&wiff=` | `unified_results` | unified rows + sample routing |
| `GET /api/bridge/extract-xic?...` | `extract_xic` | `XicChromatogram` |

### Pagination

`results_matrix()` is `samples × compounds` (e.g. 2 × 2471 ≈ 5k cells). For a
results grid, flatten and paginate **on the server**:

```python
rows = [unified_to_dict(cell) for row in bridge.unified_results() for cell in row]
return {"total": len(rows), "items": rows[offset:offset + size]}
```

---

## 8. Path safety & file resolution

Never let the client pass an arbitrary filesystem path. Constrain to a
configured data root and reject traversal:

```python
from pathlib import Path
DATA_ROOT = Path(os.environ.get("PYWIFF2_DATA_ROOT", "data")).resolve()

def resolve(rel: str) -> Path:
    p = (DATA_ROOT / rel).resolve()
    if DATA_ROOT not in p.parents and p != DATA_ROOT:
        raise HTTPException(400, "path escapes data root")
    if not p.exists():
        raise HTTPException(404, "file not found")
    return p
```

Expose a discovery endpoint so the UI lists available files rather than guessing:

```python
@app.get("/api/files")
def files():
    return {
        "wiff2": [str(p.relative_to(DATA_ROOT)) for p in DATA_ROOT.rglob("*.wiff2")],
        "qsession": [str(p.relative_to(DATA_ROOT)) for p in DATA_ROOT.rglob("*.qsession")],
    }
```

---

## 9. CORS for local dev

When the SPA runs on a different port (e.g. Vite on 5173):

```python
from fastapi.middleware.cors import CORSMiddleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173"],
    allow_methods=["GET"],
    allow_headers=["*"],
)
```

---

## 10. Error handling

Map pyx500r's exceptions (see [API.md §Exceptions](./API.md#exceptions)) to HTTP
status codes once, via exception handlers:

```python
from fastapi import Request
from fastapi.responses import JSONResponse

@app.exception_handler(IndexError)
async def _index_err(_: Request, exc: IndexError):
    return JSONResponse(status_code=404, content={"detail": str(exc)})

@app.exception_handler(KeyError)
async def _key_err(_: Request, exc: KeyError):
    return JSONResponse(status_code=404, content={"detail": str(exc)})

@app.exception_handler(ValueError)
async def _value_err(_: Request, exc: ValueError):
    return JSONResponse(status_code=400, content={"detail": str(exc)})
```

`FileNotFoundError` from a reader → 404; `OSError` from a qsession that can't be
decrypted → 422/500 depending on whether it's user input.

---

## 11. TypeScript front-end

Use the generated [`types.ts`](./types.ts) verbatim. A minimal typed client:

```ts
import type {
  SampleInfo, ExperimentInfo, SpectrumData, Chromatogram,
  CompoundInfo, QuantPeakInfo, UnifiedPeakDTO,
} from "./types";

const API = "http://localhost:8000";

async function getJSON<T>(path: string): Promise<T> {
  const res = await fetch(`${API}${path}`);
  if (!res.ok) throw new Error(`${res.status} ${await res.text()}`);
  return res.json() as Promise<T>;
}

export const api = {
  samples: (f: string) =>
    getJSON<SampleInfo[]>(`/api/wiff/${encodeURIComponent(f)}/samples`),
  experiments: (f: string, sample = 0) =>
    getJSON<ExperimentInfo[]>(`/api/wiff/${encodeURIComponent(f)}/experiments?sample=${sample}`),
  tic: (f: string, sample = 0, experiment = 0) =>
    getJSON<Chromatogram>(`/api/wiff/${encodeURIComponent(f)}/tic?sample=${sample}&experiment=${experiment}`),
  spectrum: (f: string, experiment: number, cycle: number, centroid = false) =>
    getJSON<SpectrumData>(
      `/api/wiff/${encodeURIComponent(f)}/spectrum?experiment=${experiment}&cycle=${cycle}&centroid=${centroid}`),
  compounds: (f: string) =>
    getJSON<CompoundInfo[]>(`/api/qsession/${encodeURIComponent(f)}/compounds`),
  results: (f: string, page = 0, size = 100) =>
    getJSON<{ total: number; items: UnifiedPeakDTO[] }>(
      `/api/qsession/${encodeURIComponent(f)}/results?page=${page}&size=${size}`),
};
```

Plotting tips:
- Use a WebGL chart (e.g. `uPlot`, `plotly.js` scattergl, or `regl`) for spectra
  — 100k+ points will choke SVG/canvas-DOM charts.
- For the TIC/XIC traces, `times`/`intensities` are parallel arrays ready for
  most charting libs.

---

## 12. Production checklist

- [ ] `uvicorn app:app --workers N` (each worker = its own `ReaderPool`).
- [ ] Bound `ReaderPool.max_open` to fit memory (decrypted DB + scan mmap per file).
- [ ] Offload decode with `run_in_threadpool`; consider a `ProcessPoolExecutor`
      for batch endpoints.
- [ ] Centroid/downsample large spectra; offer the binary endpoint for raw data.
- [ ] Lock down file paths to a data root; never trust client paths.
- [ ] Add response caching (e.g. `Cache-Control`) — spectra/TICs are immutable.
- [ ] Install the base package in the server image; `numba` is included for fast decode.
- [ ] Pre-warm numba JIT at startup (import `pyx500r.tof` triggers warm-up) so
      the first request isn't slow.

---

## See also

- [`API.md`](./API.md) — full Python API reference
- [`types.ts`](./types.ts) — TypeScript interfaces
- [`examples/server/`](../examples/server/) — runnable FastAPI app
