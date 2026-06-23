"""LibraryView browser — FastAPI backend + bundled single-page GUI.

Browse a SCIEX LibraryView ``.sqlite`` database: search compounds, inspect
their settings/metadata, and view their (possibly multiple) reference spectra
rendered as interactive plots.

Run::

    pip install -e ".[server]"
    export PYX500R_LIBRARY_DB=path/to/your/libview.sqlite
    uvicorn examples.library_browser.app:app --reload --port 8001

Then open http://localhost:8001/ for the GUI, or /docs for the API.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

try:
    from fastapi import FastAPI, HTTPException, Query, Request
    from fastapi.concurrency import run_in_threadpool
    from fastapi.responses import FileResponse, JSONResponse
    from fastapi.staticfiles import StaticFiles
except ImportError as exc:  # pragma: no cover - dependency hint
    raise SystemExit(
        "FastAPI is not installed. Run:\n"
        "    pip install -e \".[server]\""
    ) from exc

from .library_db import LibraryDB

# ── configuration ───────────────────────────────────────────────────────────
DB_PATH = os.environ.get("PYX500R_LIBRARY_DB", "data/libview.sqlite")
_STATIC = Path(__file__).resolve().parent / "static"

_db: LibraryDB | None = None


def get_db() -> LibraryDB:
    global _db
    if _db is None:
        _db = LibraryDB(DB_PATH)
    return _db


app = FastAPI(
    title="pyx500r LibraryView browser",
    version="0.1.0",
    description="Browse compounds and reference spectra from a SCIEX LibraryView DB.",
)


@app.exception_handler(FileNotFoundError)
async def _fnf(_: Request, exc: FileNotFoundError) -> JSONResponse:
    return JSONResponse(status_code=404, content={"detail": str(exc)})


@app.exception_handler(ValueError)
async def _ve(_: Request, exc: ValueError) -> JSONResponse:
    return JSONResponse(status_code=400, content={"detail": str(exc)})


# ── API ──────────────────────────────────────────────────────────────────────
@app.get("/api/stats")
async def stats() -> dict[str, Any]:
    return await run_in_threadpool(get_db().stats)


@app.get("/api/libraries")
async def libraries() -> list[dict[str, Any]]:
    return await run_in_threadpool(get_db().libraries)


@app.get("/api/compounds")
async def compounds(
    search: str | None = None,
    library_id: str | None = None,
    offset: int = 0,
    limit: int = Query(default=100, le=500),
) -> dict[str, Any]:
    return await run_in_threadpool(
        lambda: get_db().list_compounds(
            search=search, library_id=library_id, offset=offset, limit=limit
        )
    )


@app.get("/api/compounds/{compound_id}")
async def compound_detail(compound_id: str) -> dict[str, Any]:
    result = await run_in_threadpool(get_db().get_compound, compound_id)
    if result is None:
        raise HTTPException(status_code=404, detail="compound not found")
    return result


@app.get("/api/compounds/{compound_id}/spectra")
async def compound_spectra(compound_id: str) -> list[dict[str, Any]]:
    return await run_in_threadpool(get_db().list_spectra, compound_id)


@app.get("/api/spectra/{spectrum_id}")
async def spectrum(spectrum_id: str, kind: str = "centroid") -> dict[str, Any]:
    result = await run_in_threadpool(get_db().get_spectrum, spectrum_id, kind)
    if result is None:
        raise HTTPException(status_code=404, detail="spectrum not found")
    return result


# ── GUI (served last so /api/* wins) ─────────────────────────────────────────
@app.get("/")
async def index() -> FileResponse:
    return FileResponse(_STATIC / "index.html")


app.mount("/", StaticFiles(directory=str(_STATIC)), name="static")
