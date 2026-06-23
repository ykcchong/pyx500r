"""Reference FastAPI server exposing pyx500r over JSON HTTP.

A minimal but production-shaped API for a TypeScript GUI that browses
``.wiff2`` acquisitions and ``.qsession`` quantitation results.

Run::

    pip install -e ".[server]"
    export PYX500R_DATA_ROOT=./data        # directory containing your files
    uvicorn examples.server.app:app --reload --port 8000

Then open http://localhost:8000/docs for interactive Swagger UI.

This module is intentionally dependency-light and self-contained; see
../../docs/GUI_INTEGRATION.md for the design rationale and extension points.
"""

from __future__ import annotations

import os
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Any

import numpy as np

try:
    from fastapi import FastAPI, HTTPException, Query, Request, Response
    from fastapi.concurrency import run_in_threadpool
    from fastapi.middleware.cors import CORSMiddleware
    from fastapi.responses import JSONResponse
except ImportError as exc:  # pragma: no cover - dependency hint
    raise SystemExit(
        "FastAPI is not installed. Run:\n"
        "    pip install -e \".[server]\""
    ) from exc

from pyx500r import WiffQSessionBridge

from .pool import ReaderPool
from .serializers import (
    dataclass_to_dict,
    serialize_chromatogram,
    unified_to_dict,
)

# ── configuration ──────────────────────────────────────────────────────────
DATA_ROOT = Path(os.environ.get("PYX500R_DATA_ROOT", "data")).resolve()
MAX_OPEN = int(os.environ.get("PYX500R_MAX_OPEN", "16"))
CORS_ORIGINS = os.environ.get(
    "PYX500R_CORS_ORIGINS", "http://localhost:5173,http://localhost:3000"
).split(",")

pool = ReaderPool(max_open=MAX_OPEN)


@asynccontextmanager
async def _lifespan(_: "FastAPI"):
    # startup: nothing to pre-warm (readers are opened lazily, per file)
    yield
    # shutdown: release all pooled readers + their dedicated threads
    pool.close_all()


app = FastAPI(
    title="pyx500r API",
    version="0.2.0",
    description="JSON access to SCIEX .wiff2 / .qsession files.",
    lifespan=_lifespan,
)
app.add_middleware(
    CORSMiddleware,
    allow_origins=CORS_ORIGINS,
    allow_methods=["GET"],
    allow_headers=["*"],
)


# ── path safety ─────────────────────────────────────────────────────────────
def resolve(rel: str, suffix: str) -> Path:
    """Resolve a client-supplied relative path inside DATA_ROOT, safely."""
    p = (DATA_ROOT / rel).resolve()
    if p != DATA_ROOT and DATA_ROOT not in p.parents:
        raise HTTPException(status_code=400, detail="path escapes data root")
    if p.suffix.lower() != suffix:
        raise HTTPException(status_code=400, detail=f"expected a {suffix} file")
    if not p.exists():
        raise HTTPException(status_code=404, detail=f"file not found: {rel}")
    return p


# ── exception mapping ────────────────────────────────────────────────────────
@app.exception_handler(IndexError)
async def _index_err(_: Request, exc: IndexError) -> JSONResponse:
    return JSONResponse(status_code=404, content={"detail": str(exc)})


@app.exception_handler(KeyError)
async def _key_err(_: Request, exc: KeyError) -> JSONResponse:
    return JSONResponse(status_code=404, content={"detail": str(exc)})


@app.exception_handler(FileNotFoundError)
async def _fnf_err(_: Request, exc: FileNotFoundError) -> JSONResponse:
    return JSONResponse(status_code=404, content={"detail": str(exc)})


@app.exception_handler(ValueError)
async def _value_err(_: Request, exc: ValueError) -> JSONResponse:
    return JSONResponse(status_code=400, content={"detail": str(exc)})


# ── discovery ────────────────────────────────────────────────────────────────
@app.get("/api/files")
def list_files() -> dict[str, list[str]]:
    """List the .wiff2 and .qsession files available under the data root."""
    return {
        "wiff2": sorted(str(p.relative_to(DATA_ROOT)) for p in DATA_ROOT.rglob("*.wiff2")),
        "qsession": sorted(str(p.relative_to(DATA_ROOT)) for p in DATA_ROOT.rglob("*.qsession")),
    }


# ── WIFF2 endpoints ──────────────────────────────────────────────────────────
@app.get("/api/wiff/{file:path}/samples")
async def wiff_samples(file: str) -> list[dict[str, Any]]:
    path = resolve(file, ".wiff2")
    return await run_in_threadpool(
        pool.with_wiff, path,
        lambda r: [dataclass_to_dict(s) for s in r.list_samples()],
    )


@app.get("/api/wiff/{file:path}/experiments")
async def wiff_experiments(file: str, sample: int = 0) -> list[dict[str, Any]]:
    path = resolve(file, ".wiff2")
    return await run_in_threadpool(
        pool.with_wiff, path,
        lambda r: [dataclass_to_dict(e) for e in r.get_experiments(sample)],
    )


@app.get("/api/wiff/{file:path}/instruments")
async def wiff_instruments(file: str, sample: int = 0) -> list[dict[str, Any]]:
    path = resolve(file, ".wiff2")
    return await run_in_threadpool(
        pool.with_wiff, path,
        lambda r: [dataclass_to_dict(i) for i in r.list_instruments(sample)],
    )


@app.get("/api/wiff/{file:path}/tic")
async def wiff_tic(file: str, sample: int = 0, experiment: int = 0) -> dict[str, Any]:
    path = resolve(file, ".wiff2")
    return await run_in_threadpool(
        pool.with_wiff, path,
        lambda r: serialize_chromatogram(
            r.get_experiment_tic(sample_index=sample, experiment_index=experiment)
        ),
    )


@app.get("/api/wiff/{file:path}/spectrum/meta")
async def wiff_spectrum_meta(
    file: str, sample: int = 0, experiment: int = 0, cycle: int = 0
) -> dict[str, Any]:
    path = resolve(file, ".wiff2")
    return await run_in_threadpool(
        pool.with_wiff, path,
        lambda r: dataclass_to_dict(r.get_spectrum_metadata(sample, experiment, cycle)),
    )


def _spectrum_payload(
    reader, sample: int, experiment: int, cycle: int, centroid: bool,
    mz_min: float | None, mz_max: float | None,
) -> dict[str, Any]:
    spec = reader.get_spectrum(
        sample, experiment, cycle, centroid=centroid, return_arrays=True
    )
    mz = np.asarray(spec.mz)
    it = np.asarray(spec.intensities)
    if mz_min is not None or mz_max is not None:
        lo = -np.inf if mz_min is None else mz_min
        hi = np.inf if mz_max is None else mz_max
        mask = (mz >= lo) & (mz <= hi)
        mz, it = mz[mask], it[mask]
    return {
        "sample_index": spec.sample_index,
        "experiment_index": spec.experiment_index,
        "cycle_index": spec.cycle_index,
        "scan_time": spec.scan_time,
        "mz": mz.tolist(),
        "intensities": it.tolist(),
        "centroided": spec.centroided,
        "precursor_mz": spec.precursor_mz,
        "isolation_target_mz": spec.isolation_target_mz,
        "isolation_lower_offset": spec.isolation_lower_offset,
        "isolation_upper_offset": spec.isolation_upper_offset,
    }


@app.get("/api/wiff/{file:path}/spectrum")
async def wiff_spectrum(
    file: str,
    sample: int = 0,
    experiment: int = 0,
    cycle: int = 0,
    centroid: bool = False,
    mz_min: float | None = Query(default=None),
    mz_max: float | None = Query(default=None),
) -> dict[str, Any]:
    """Return one spectrum as JSON.

    Tip: pass ``centroid=true`` (or an ``mz_min``/``mz_max`` window) to keep the
    payload small. A dense profile scan can be ~156k points.
    """
    path = resolve(file, ".wiff2")
    return await run_in_threadpool(
        pool.with_wiff, path,
        lambda r: _spectrum_payload(r, sample, experiment, cycle, centroid, mz_min, mz_max),
    )


def _spectrum_bin(reader, sample: int, experiment: int, cycle: int, centroid: bool) -> bytes:
    spec = reader.get_spectrum(
        sample, experiment, cycle, centroid=centroid, return_arrays=True
    )
    mz = np.ascontiguousarray(spec.mz, dtype=np.float32)
    it = np.ascontiguousarray(spec.intensities, dtype=np.float32)
    return len(mz).to_bytes(4, "little") + mz.tobytes() + it.tobytes()


@app.get("/api/wiff/{file:path}/spectrum.bin")
async def wiff_spectrum_bin(
    file: str, sample: int = 0, experiment: int = 0, cycle: int = 0, centroid: bool = False
) -> Response:
    """Binary spectrum: uint32 count + float32 mz[] + float32 intensity[].

    Much smaller and faster to parse in TypeScript than JSON for raw profiles.
    See docs/GUI_INTEGRATION.md §6 for the decode snippet.
    """
    path = resolve(file, ".wiff2")
    payload = await run_in_threadpool(
        pool.with_wiff, path,
        lambda r: _spectrum_bin(r, sample, experiment, cycle, centroid),
    )
    return Response(content=payload, media_type="application/octet-stream")


# ── QSession endpoints ───────────────────────────────────────────────────────
@app.get("/api/qsession/{file:path}/info")
async def qsession_info(file: str) -> dict[str, Any]:
    path = resolve(file, ".qsession")

    def work(r) -> dict[str, Any]:
        return {
            "version": r.version,
            "qmap_version": r.qmap_version,
            "locked": r.locked,
            "sample_count": len(r.list_samples()),
            "compound_count": len(r.list_compounds()),
            "xic_count": r.xic_count,
        }

    return await run_in_threadpool(pool.with_qsession, path, work)


@app.get("/api/qsession/{file:path}/compounds")
async def qsession_compounds(file: str) -> list[dict[str, Any]]:
    path = resolve(file, ".qsession")
    return await run_in_threadpool(
        pool.with_qsession, path,
        lambda r: [dataclass_to_dict(c) for c in r.list_compounds()],
    )


@app.get("/api/qsession/{file:path}/samples")
async def qsession_samples(file: str) -> list[dict[str, Any]]:
    path = resolve(file, ".qsession")
    return await run_in_threadpool(
        pool.with_qsession, path,
        lambda r: [dataclass_to_dict(s) for s in r.list_samples()],
    )


@app.get("/api/qsession/{file:path}/peak")
async def qsession_peak(file: str, sample: int, compound: int) -> dict[str, Any]:
    path = resolve(file, ".qsession")
    return await run_in_threadpool(
        pool.with_qsession, path,
        lambda r: dataclass_to_dict(r.get_peak(sample, compound)),
    )


@app.get("/api/qsession/{file:path}/chromatogram")
async def qsession_chromatogram(file: str, sample: int, compound: int) -> dict[str, Any]:
    path = resolve(file, ".qsession")

    def work(r) -> dict[str, Any] | None:
        chrom = r.get_chromatogram(sample, compound)
        return serialize_chromatogram(chrom) if chrom is not None else None

    result = await run_in_threadpool(pool.with_qsession, path, work)
    if result is None:
        raise HTTPException(status_code=404, detail="no chromatogram for this cell")
    return result


# ── Bridge endpoint: paginated unified results ───────────────────────────────
@app.get("/api/bridge/results")
async def bridge_results(
    qs: str,
    wiff: list[str] = Query(default=[]),
    match_by: str = "name",
    page: int = 0,
    size: int = Query(default=100, le=1000),
) -> dict[str, Any]:
    """Flattened, paginated compound×sample results across a qsession+wiffs.

    ``qs`` and each ``wiff`` are paths relative to the data root. The bridge is
    opened per request (it owns its own readers); for heavy use, cache results.
    """
    qs_path = resolve(qs, ".qsession")
    wiff_paths = [resolve(w, ".wiff2") for w in wiff]

    def work() -> dict[str, Any]:
        with WiffQSessionBridge(qs_path, wiff_paths, match_by=match_by) as br:
            rows = [
                unified_to_dict(cell)
                for row in br.unified_results()
                for cell in row
            ]
        total = len(rows)
        start = page * size
        return {
            "total": total,
            "page": page,
            "size": size,
            "items": rows[start : start + size],
        }

    return await run_in_threadpool(work)


@app.get("/api/health")
def health() -> dict[str, Any]:
    return {"status": "ok", "data_root": str(DATA_ROOT)}
