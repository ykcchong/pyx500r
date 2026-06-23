# pyx500r reference server

A minimal, production-shaped FastAPI app that exposes `pyx500r` over JSON HTTP
for a TypeScript GUI. See [`../../docs/GUI_INTEGRATION.md`](../../docs/GUI_INTEGRATION.md)
for the full design guide and [`../../docs/types.ts`](../../docs/types.ts) for
the matching TypeScript interfaces.

## Install & run

```bash
# from the repo root
pip install -e ".[server]"

export PYX500R_DATA_ROOT=./data          # folder containing your .wiff2/.qsession
uvicorn examples.server.app:app --reload --port 8000
```

Open http://localhost:8000/docs for interactive Swagger UI, or
http://localhost:8000/api/health to verify it's up.

## Configuration (environment variables)

| Variable | Default | Meaning |
|----------|---------|---------|
| `PYX500R_DATA_ROOT` | `data` | Directory searched for files; all paths are constrained inside it |
| `PYX500R_MAX_OPEN` | `16` | Max readers kept open in the LRU pool |
| `PYX500R_CORS_ORIGINS` | `http://localhost:5173,http://localhost:3000` | Comma-separated allowed origins |

## Endpoints

Discovery:
- `GET /api/files` — list available `.wiff2` / `.qsession` paths
- `GET /api/health`

WIFF2 (`{file}` is a path relative to the data root):
- `GET /api/wiff/{file}/samples`
- `GET /api/wiff/{file}/experiments?sample=0`
- `GET /api/wiff/{file}/instruments?sample=0`
- `GET /api/wiff/{file}/tic?sample=0&experiment=0`
- `GET /api/wiff/{file}/spectrum/meta?sample=0&experiment=0&cycle=0`
- `GET /api/wiff/{file}/spectrum?sample=0&experiment=0&cycle=0&centroid=false&mz_min=&mz_max=`
- `GET /api/wiff/{file}/spectrum.bin?...` — binary `uint32 n + float32 mz[n] + float32 int[n]`

QSession:
- `GET /api/qsession/{file}/info`
- `GET /api/qsession/{file}/compounds`
- `GET /api/qsession/{file}/samples`
- `GET /api/qsession/{file}/peak?sample=&compound=`
- `GET /api/qsession/{file}/chromatogram?sample=&compound=`

Bridge:
- `GET /api/bridge/results?qs=...&wiff=...&wiff=...&page=0&size=100`

## Example

```bash
curl "http://localhost:8000/api/files"
curl "http://localhost:8000/api/wiff/YOUR_FILE.wiff2/samples"
curl "http://localhost:8000/api/wiff/YOUR_FILE.wiff2/spectrum?experiment=0&cycle=51&centroid=true"
```

> `{file}` is URL-encoded; `/` becomes `%2F`. The route uses a `:path`
> converter so sub-directories work.

## Notes

- Readers are pooled per file path and guarded by per-reader locks, because each
  wraps a single-thread SQLite connection. CPU-bound decode runs in a thread
  pool so it doesn't block the event loop.
- For large/raw spectra prefer `centroid=true`, an `mz_min`/`mz_max` window, or
  the `.bin` endpoint.
- This app is a starting point — add auth, caching, and rate limiting for real
  deployments.
