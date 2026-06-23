"""WIFF2 GUI — native Windows viewer for SCIEX X500R acquisitions.

Provides:
* Open ``.wiff2`` files and browse samples/experiments.
* Extract and display extracted ion chromatograms (XIC) from TOF-MS data
  by specifying a target m/z and tolerance.  The MS1 experiment is
  auto-selected — no manual experiment picking required.
* Input a chemical formula and adduct (H+/Na+/NH4+) to calculate the
  monoisotopic m/z automatically (requires ``pyteomics``).
* Automatically find all MS/MS (TOFMSMS) spectra whose precursor m/z
  matches the target within tolerance, and display their stick spectra.

Usage::

    x500rgui
    # or
    python -m pyx500r.wiff_gui

Dependencies: ``matplotlib``.  Optional: ``pyteomics`` for formula mass calculation.
"""

from __future__ import annotations

import csv
import json
import queue
import threading
import tkinter as tk
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from typing import Any

import numpy as np

# ── matplotlib with TkAgg backend ──────────────────────────────────────────
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.transforms import Bbox
from matplotlib.ticker import FuncFormatter, MaxNLocator

from pyx500r.reader import WiffReader
from pyx500r.models import Chromatogram, ExperimentInfo, SampleInfo, SpectrumData

# ── optional: formula mass calculation ───────────────────────────────────
try:
    from pyteomics import mass as _pyteomics_mass
    _HAS_PYTEOMICS = True
except ImportError:
    _pyteomics_mass = None  # type: ignore[assignment]
    _HAS_PYTEOMICS = False

# Known adducts for positive-ion mode
_ADDUCTS = {
    "H+":     {"formula": "H+",    "charge": 1},
    "Na+":    {"formula": "Na+",   "charge": 1},
    "NH4+":   {"formula": "NH4+",  "charge": 1},
}

_DATA_DIR = Path(__file__).resolve().parents[2] / "data"
_SPECTRA_DIR = _DATA_DIR / "spectra"
_DEFAULT_CSV_LIBRARY = _SPECTRA_DIR / "highresnps april 2026 shimadzu.csv"
_DEFAULT_JSON_LIBRARY = _DATA_DIR / "library.json"

# ── helpers ────────────────────────────────────────────────────────────────


def _calculate_formula_mass(formula: str, adduct: str) -> float | None:
    """Calculate monoisotopic [M+adduct]+ mass using pyteomics.

    Returns None if pyteomics is not installed or the formula is invalid.
    """
    if not _HAS_PYTEOMICS:
        return None
    try:
        return _pyteomics_mass.calculate_mass(  # type: ignore[union-attr]
            formula=formula,
            charge=_ADDUCTS[adduct]["charge"],
            charge_carrier=_ADDUCTS[adduct]["formula"],
        )
    except Exception:
        return None


def _compute_isotopic_distribution(
    formula: str, adduct_key: str,
) -> list[tuple[float, float]]:
    """Compute theoretical isotopic distribution for formula + adduct.

    Returns list of (m/z, relative_abundance_pct) sorted by m/z, where the
    monoisotopic peak is normalised to 100 %.  Peaks below 0.1 % are dropped.
    """
    if not _HAS_PYTEOMICS:
        return []
    try:
        adduct_info = _ADDUCTS[adduct_key]
        charge_carrier = adduct_info["formula"]
        charge = adduct_info["charge"]

        isotopologues = _pyteomics_mass.isotopologues(formula=formula)  # type: ignore[union-attr]
        result: list[tuple[float, float]] = []
        for composition in isotopologues:
            mz_val = _pyteomics_mass.calculate_mass(  # type: ignore[union-attr]
                composition=composition, charge=charge, charge_carrier=charge_carrier,
            )
            abund = _pyteomics_mass.isotopic_composition_abundance(  # type: ignore[union-attr]
                formula=formula, composition=composition,
            )
            if abund >= 0.001:  # drop < 0.1 %
                result.append((mz_val, abund * 100.0))
        # Normalise so the most abundant isotopologue = 100 %
        if result:
            max_abund = max(a for _, a in result)
            result = [(mz, a / max_abund * 100.0) for mz, a in result]
        return sorted(result, key=lambda x: x[0])
    except Exception:
        return []


def _ppm_to_da(mz: float, ppm: float) -> float:
    """Convert ppm tolerance at a given m/z to a ± Da half-window."""
    return mz * ppm * 1e-6


def _compact_sci(value: float) -> str:
    """Return compact scientific notation such as 3e5 for axis labels."""
    if not np.isfinite(value) or value == 0:
        return "0"
    mantissa, exponent = f"{value:.1e}".split("e")
    mantissa = mantissa.rstrip("0").rstrip(".")
    return f"{mantissa}e{int(exponent)}"


def _style_dense_axes(ax: Any, *, grid: bool = True) -> None:
    """Use compact margins/ticks suitable for embedded Tk plots."""
    ax.tick_params(axis="both", which="major", labelsize=8, pad=1)
    ax.title.set_size(9)
    if grid:
        ax.grid(True, axis="both", color="#d0d0d0", linewidth=0.45, alpha=0.8)
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)


def _remove_artists(artists: list[Any]) -> None:
    for artist in artists:
        try:
            artist.remove()
        except Exception:
            pass


def _label_bbox(
    ax: Any,
    x: float,
    y: float,
    text: str,
    fontsize: float,
    rotation: float,
    *,
    ha: str = "center",
    va: str = "bottom",
) -> Bbox:
    """Estimate a text label bbox in display pixels for collision checks."""
    x_px, y_px = ax.transData.transform((x, y))
    dpi_scale = ax.figure.dpi / 72.0
    char_w = fontsize * 0.58 * dpi_scale
    char_h = fontsize * 1.15 * dpi_scale
    if rotation % 180:
        width = char_h * 1.4
        height = max(char_w * len(text), char_h)
    else:
        width = max(char_w * len(text), char_h)
        height = char_h * 1.4

    if ha == "left":
        x0 = x_px
    elif ha == "right":
        x0 = x_px - width
    else:
        x0 = x_px - width / 2

    if va == "center":
        y0 = y_px - height / 2
    elif va == "top":
        y0 = y_px - height
    else:
        y0 = y_px

    return Bbox.from_bounds(x0, y0, width, height)


def _place_peak_labels(
    ax: Any,
    mz: np.ndarray,
    rel: np.ndarray,
    *,
    label_for: Any,
    y_values: np.ndarray | None = None,
    threshold_pct: float = 2.0,
    max_labels: int | None = None,
    color: str = "#333333",
    fontsize: int = 7,
    rotation: int = 90,
) -> list[Any]:
    """Place peak labels above visible peaks, skipping labels that would overlap."""
    if mz.size == 0 or rel.size == 0:
        return []
    y = rel if y_values is None else y_values

    x_lo, x_hi = ax.get_xlim()
    y_lo, y_hi = ax.get_ylim()
    visible = (mz >= x_lo) & (mz <= x_hi) & (rel >= threshold_pct) & (y > max(y_lo, 0))
    if not visible.any():
        return []

    candidates = np.flatnonzero(visible)
    candidates = candidates[np.argsort(rel[candidates])[::-1]]
    if max_labels is not None:
        candidates = candidates[:max_labels]

    y_span = max(y_hi - y_lo, 1.0)
    y_pad = max(y_span * 0.025, 1.0)
    axes_bbox = ax.bbox.padded(-2)
    accepted: list[Bbox] = []
    artists: list[Any] = []

    for idx in candidates:
        label = label_for(idx)
        x = float(mz[idx])
        peak_y = float(y[idx])
        y_text = min(peak_y + y_pad, y_hi - y_pad)

        x_px, y_px = ax.transData.transform((x, min(max(peak_y, y_lo), y_hi)))
        right_x = ax.transData.inverted().transform((x_px + 5, y_px))[0]
        left_x = ax.transData.inverted().transform((x_px - 5, y_px))[0]
        attempts = [
            (x, y_text, rotation, "center", "bottom"),
            (right_x, min(max(peak_y, y_lo), y_hi), 0, "left", "center"),
            (left_x, min(max(peak_y, y_lo), y_hi), 0, "right", "center"),
        ]

        for lx, ly, lrot, ha, va in attempts:
            bbox = _label_bbox(ax, lx, ly, label, fontsize, lrot, ha=ha, va=va)
            if not axes_bbox.contains(bbox.x0, bbox.y0) or not axes_bbox.contains(bbox.x1, bbox.y1):
                continue
            padded = bbox.padded(2)
            if any(padded.overlaps(existing) for existing in accepted):
                continue

            artist = ax.text(
                lx, ly, label,
                ha=ha, va=va, fontsize=fontsize, rotation=lrot,
                color=color, clip_on=True,
            )
            artists.append(artist)
            accepted.append(padded)
            break

    return artists


class _ExtractionCancelled(Exception):
    """Raised internally when the user cancels a running extraction."""


def _sum_sorted_window(mz_arr: np.ndarray, int_arr: np.ndarray, mz_lo: float, mz_hi: float) -> float:
    """Sum intensities in a sorted m/z half-open index window."""
    if mz_arr.size == 0:
        return 0.0
    lo = int(np.searchsorted(mz_arr, mz_lo, side="left"))
    hi = int(np.searchsorted(mz_arr, mz_hi, side="right"))
    if hi <= lo:
        return 0.0
    return float(int_arr[lo:hi].sum())


def _extract_xic_from_ms1(
    reader: WiffReader,
    sample_index: int,
    ms1_experiment_index: int,
    target_mz: float,
    mz_tolerance_da: float,
    rt_start: float | None = None,
    rt_end: float | None = None,
    progress_cb: Callable[[int, int], None] | None = None,
    cancel_event: threading.Event | None = None,
) -> Chromatogram | None:
    """Extract an XIC by summing intensities in the m/z window for every MS1 cycle."""
    mz_lo = target_mz - mz_tolerance_da
    mz_hi = target_mz + mz_tolerance_da

    try:
        cycle_times = reader.get_cycle_times(sample_index, ms1_experiment_index)
    except Exception:
        return None
    if not cycle_times:
        return None

    reader.prefetch_experiment(sample_index, ms1_experiment_index)

    times: list[float] = []
    intensities: list[float] = []
    total = len(cycle_times)

    for cycle in range(total):
        if cancel_event is not None and cancel_event.is_set():
            raise _ExtractionCancelled
        t = cycle_times[cycle]
        if rt_start is not None and t < rt_start:
            continue
        if rt_end is not None and t > rt_end:
            continue
        try:
            spec = reader.get_spectrum(
                sample_index, ms1_experiment_index, cycle,
                centroid=False, return_arrays=True,
            )
        except Exception:
            continue
        mz_arr = np.asarray(spec.mz)
        int_arr = np.asarray(spec.intensities)
        total_int = _sum_sorted_window(mz_arr, int_arr, mz_lo, mz_hi)
        times.append(t)
        intensities.append(total_int)
        if progress_cb is not None:
            progress_cb(cycle + 1, total)

    if not times:
        return None
    return Chromatogram(
        times=times,
        intensities=intensities,
        experiment_index=ms1_experiment_index,
        ms_level=1,
    )


@dataclass
class _MsMsMatch:
    """A single MS/MS spectrum whose precursor matches the target m/z."""
    experiment_index: int
    cycle_index: int
    scan_time: float
    precursor_mz: float
    tic_50_parent: float  # total ion current in [50, precursor_mz]


@dataclass
class _ExtractionResult:
    """Completed XIC extraction and MS/MS match payload."""
    xic: Chromatogram | None
    msms_matches: list[_MsMsMatch]
    target_mz: float
    tol_da: float
    ms1_idx: int


@dataclass
class _LibraryHit:
    """One matchms library search hit."""
    rank: int
    name: str
    formula: str
    database_id: str
    precursor_mz: float
    score: float
    matches: int
    num_peaks: int
    mz: np.ndarray
    intensity: np.ndarray


def _float_or_none(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _csv_row_to_matchms_json(row: dict[str, str]) -> dict[str, Any] | None:
    peaks: list[list[float]] = []
    for idx in range(1, 8):
        mz = _float_or_none(row.get(f"m/z {idx}"))
        intensity = _float_or_none(row.get(f"Intensity {idx}"))
        if mz is None or intensity is None or intensity <= 0:
            continue
        peaks.append([mz, intensity])
    if not peaks:
        return None
    peaks.sort(key=lambda item: item[0])

    precursor_mz = _float_or_none(row.get("Precursor m/z"))
    metadata: dict[str, Any] = {
        "compound_name": row.get("Compound Name", "").strip(),
        "formula": row.get("Formula", "").strip(),
        "precursor_mz": precursor_mz or 0.0,
        "ionmode": row.get("Polarity", "").strip().lower(),
        "adduct": row.get("Precursor Ion", "").strip(),
        "smiles": row.get("SMILES", "").strip(),
        "inchi": row.get("InChI", "").strip(),
        "database_id": row.get("Comment", "").strip(),
        "compound_class": row.get("Class", "").strip(),
        "theory_mw": _float_or_none(row.get("Theory MW")) or 0.0,
        "retention_time": _float_or_none(row.get("RT")) or 0.0,
        "collision_gas_voltage": _float_or_none(row.get("Collision Gas Vol.")) or 0.0,
        "peaks_json": peaks,
    }
    if row.get("CAS #"):
        metadata["cas"] = row["CAS #"].strip()
    return metadata


def _convert_csv_to_matchms_json(csv_path: Path, json_path: Path) -> int:
    """Convert the Shimadzu CSV library into MatchMS JSON."""
    spectra: list[dict[str, Any]] = []
    with csv_path.open("r", newline="", encoding="utf-8-sig", errors="replace") as handle:
        for row in csv.DictReader(handle):
            item = _csv_row_to_matchms_json(row)
            if item is not None:
                spectra.append(item)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    with json_path.open("w", encoding="utf-8") as handle:
        json.dump(spectra, handle, ensure_ascii=False, separators=(",", ":"))
    return len(spectra)


def _ensure_matchms_json_library(csv_path: Path, json_path: Path) -> None:
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV library not found: {csv_path}")
    if (
        not json_path.exists()
        or json_path.stat().st_mtime < csv_path.stat().st_mtime
    ):
        _convert_csv_to_matchms_json(csv_path, json_path)


def _load_matchms_json_library(path: Path) -> list[dict[str, Any]]:
    with path.open("r", encoding="utf-8") as handle:
        raw = json.load(handle)
    library: list[dict[str, Any]] = []
    for item in raw:
        peaks = item.get("peaks_json", [])
        if not peaks:
            continue
        mz = np.array([float(p[0]) for p in peaks], dtype=np.float64)
        intensity = np.array([float(p[1]) for p in peaks], dtype=np.float64)
        library.append({
            "name": str(item.get("compound_name") or item.get("name") or ""),
            "formula": str(item.get("formula") or ""),
            "database_id": str(item.get("database_id") or ""),
            "precursor_mz": float(item.get("precursor_mz") or 0.0),
            "ionmode": str(item.get("ionmode") or ""),
            "mz": mz,
            "intensity": intensity,
            "num_peaks": int(len(mz)),
        })
    return library


def _search_json_library_with_matchms(
    query_mz: np.ndarray,
    query_intensity: np.ndarray,
    library: list[dict[str, Any]],
    *,
    precursor_mz: float | None = None,
    precursor_ppm: float = 50.0,
    tolerance_da: float = 0.02,
    top_n: int = 10,
) -> list[_LibraryHit]:
    """Search a centroided MS/MS peak list against MatchMS JSON library entries."""
    try:
        from matchms import Spectrum
        from matchms.similarity import CosineGreedy
    except ImportError as exc:
        raise RuntimeError("matchms is not installed. Install with: pip install -e '.[gui]'") from exc

    if query_mz.size == 0 or query_intensity.size == 0:
        return []

    query = Spectrum(
        mz=np.asarray(query_mz, dtype=np.float64),
        intensities=np.asarray(query_intensity, dtype=np.float64),
        metadata={"precursor_mz": 0.0},
    )
    similarity = CosineGreedy(tolerance=tolerance_da)
    scored: list[tuple[float, int, dict[str, Any]]] = []

    for entry in library:
        ref_precursor = float(entry.get("precursor_mz") or 0.0)
        if precursor_mz and ref_precursor > 0:
            tol = precursor_mz * precursor_ppm * 1e-6
            if abs(ref_precursor - precursor_mz) > tol:
                continue
        ref = Spectrum(
            mz=entry["mz"],
            intensities=entry["intensity"],
            metadata={
                "compound_name": entry.get("name", ""),
                "formula": entry.get("formula", ""),
                "database_id": entry.get("database_id", ""),
                "precursor_mz": 0.0,
            },
        )
        pair = similarity.pair(query, ref)
        score = float(pair["score"])
        matches = int(pair["matches"])
        if score > 0 and matches > 0:
            scored.append((score, matches, entry))

    scored.sort(key=lambda item: (item[0], item[1]), reverse=True)
    hits: list[_LibraryHit] = []
    for rank, (score, matches, entry) in enumerate(scored[:top_n], start=1):
        hits.append(_LibraryHit(
            rank=rank,
            name=str(entry.get("name", "")),
            formula=str(entry.get("formula", "")),
            database_id=str(entry.get("database_id", "")),
            precursor_mz=float(entry.get("precursor_mz") or 0.0),
            score=score,
            matches=matches,
            num_peaks=int(entry.get("num_peaks", len(entry.get("mz", [])))),
            mz=np.asarray(entry.get("mz", []), dtype=np.float64),
            intensity=np.asarray(entry.get("intensity", []), dtype=np.float64),
        ))
    return hits


def _find_matching_msms(
    reader: WiffReader,
    sample_index: int,
    experiments: list[ExperimentInfo],
    target_mz: float,
    tol_da: float,
    cancel_event: threading.Event | None = None,
) -> list[_MsMsMatch]:
    """Find all MS/MS spectra whose precursor m/z falls within the window."""
    mz_lo = target_mz - tol_da
    mz_hi = target_mz + tol_da
    matches: list[_MsMsMatch] = []

    for exp in experiments:
        if exp.ms_level < 2:
            continue
        # Prefetch scan-item metadata in one batch (avoids per-cycle SQLite queries)
        try:
            reader.prefetch_experiment(sample_index, exp.index)
        except Exception:
            continue
        for ci in range(exp.cycle_count):
            if cancel_event is not None and cancel_event.is_set():
                raise _ExtractionCancelled
            row = reader._prefetch_cache.get((sample_index, exp.index, ci))
            if row is None:
                continue
            precursor_mz = reader._precursor_mz(row)
            if precursor_mz is not None and mz_lo <= precursor_mz <= mz_hi:
                # Compute TIC in [50, precursor_mz] range
                tic = 0.0
                try:
                    spec = reader.get_spectrum(
                        sample_index, exp.index, ci,
                        centroid=False, return_arrays=True,
                    )
                    mz_arr = np.asarray(spec.mz)
                    int_arr = np.asarray(spec.intensities)
                    tic = _sum_sorted_window(mz_arr, int_arr, 50.0, precursor_mz)
                except Exception:
                    pass
                matches.append(_MsMsMatch(
                    experiment_index=exp.index,
                    cycle_index=ci,
                    scan_time=float(row["retentionTime"]),
                    precursor_mz=precursor_mz,
                    tic_50_parent=tic,
                ))
    # Sort by retention time (ascending)
    matches.sort(key=lambda m: m.scan_time)
    return matches


# ── main window ────────────────────────────────────────────────────────────


class WiffGuiApp(tk.Tk):
    """Native Windows GUI for browsing WIFF2 TOF-MS and MS/MS data."""

    def __init__(self) -> None:
        super().__init__()
        self.title("pyx500r — WIFF2 Viewer")
        self.geometry("1280x860")
        self.minsize(960, 700)

        # ── state ──
        self._reader: WiffReader | None = None
        self._current_path: Path | None = None
        self._file_label_var = tk.StringVar(value="No file opened")
        self._samples: list[SampleInfo] = []
        self._experiments: list[ExperimentInfo] = []
        self._current_sample_idx = 0
        self._current_xic: Chromatogram | None = None
        self._current_msms_spectrum: SpectrumData | None = None
        self._current_target_mz: float = 0.0
        self._current_ms1_idx: int = 0
        self._current_formula: str = ""
        self._current_adduct: str = "H+"
        self._xic_click_cid: int | None = None  # matplotlib event connection id
        self._iso_xlim_cid: int | None = None
        self._iso_ylim_cid: int | None = None
        self._iso_label_artists: list[Any] = []
        self._iso_label_data: tuple[np.ndarray, np.ndarray] | None = None
        self._msms_xlim_cid: int | None = None
        self._msms_ylim_cid: int | None = None
        self._msms_label_artists: list[Any] = []
        self._msms_label_data: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None
        self._extraction_queue: queue.Queue[tuple[str, Any]] | None = None
        self._extraction_thread: threading.Thread | None = None
        self._extraction_cancel_event: threading.Event | None = None
        self._library_queue: queue.Queue[tuple[str, Any]] | None = None
        self._library_thread: threading.Thread | None = None
        self._json_library_cache: list[dict[str, Any]] | None = None
        self._library_hits_by_iid: dict[str, _LibraryHit] = {}
        self._current_library_hit: _LibraryHit | None = None

        # ── build UI ──
        self._build_menu()
        self._build_ui()

        # Keyboard shortcuts
        self.bind("<Control-o>", lambda _e: self._open_file())
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    # ── menu ──────────────────────────────────────────────────────────────
    def _build_menu(self) -> None:
        menubar = tk.Menu(self)
        self.config(menu=menubar)

        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Open WIFF2…", command=self._open_file, accelerator="Ctrl+O")
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self._on_close)
        menubar.add_cascade(label="File", menu=file_menu)

        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="About", command=self._show_about)
        menubar.add_cascade(label="Help", menu=help_menu)

    # ── main layout ───────────────────────────────────────────────────────
    def _build_ui(self) -> None:
        # ── top toolbar ──
        toolbar = ttk.Frame(self, padding=4)
        toolbar.pack(side=tk.TOP, fill=tk.X)

        ttk.Button(toolbar, text="Open WIFF2…", command=self._open_file).pack(
            side=tk.LEFT, padx=(0, 8),
        )
        ttk.Label(toolbar, textvariable=self._file_label_var, foreground="#555").pack(
            side=tk.LEFT,
        )

        # ── main paned window ──
        paned = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        paned.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=4, pady=(0, 4))

        # Left panel
        left = ttk.Frame(paned, width=340)
        paned.add(left, weight=0)

        # Right panel — vertical stack (XIC top, MS/MS bottom)
        right = ttk.Frame(paned)
        paned.add(right, weight=1)

        self._build_left_panel(left)
        self._build_right_panel(right)

        # ── status bar ──
        self._status_var = tk.StringVar(value="Ready")
        status = ttk.Label(
            self, textvariable=self._status_var,
            relief=tk.SUNKEN, anchor=tk.W, padding=(4, 2),
        )
        status.pack(side=tk.BOTTOM, fill=tk.X)

    # ── left panel ────────────────────────────────────────────────────────
    def _build_left_panel(self, parent: ttk.Frame) -> None:
        # Formula → mass calculator
        formula_frame = ttk.LabelFrame(parent, text="Formula → m/z (optional, needs pyteomics)", padding=4)
        formula_frame.pack(fill=tk.X, padx=2, pady=(2, 4))

        frow1 = ttk.Frame(formula_frame)
        frow1.pack(fill=tk.X, pady=(0, 2))
        ttk.Label(frow1, text="Formula:").pack(side=tk.LEFT)
        self._formula_var = tk.StringVar(value="")
        self._formula_entry = ttk.Entry(frow1, textvariable=self._formula_var, width=14)
        self._formula_entry.pack(side=tk.LEFT, padx=(4, 8))
        ttk.Label(frow1, text="Adduct:").pack(side=tk.LEFT)
        self._adduct_var = tk.StringVar(value="H+")
        self._adduct_combo = ttk.Combobox(
            frow1, textvariable=self._adduct_var, state="readonly",
            values=list(_ADDUCTS.keys()), width=5,
        )
        self._adduct_combo.pack(side=tk.LEFT, padx=(4, 8))
        self._calc_btn = ttk.Button(
            frow1, text="→ m/z", command=self._on_calculate_formula,
        )
        self._calc_btn.pack(side=tk.LEFT)

        # Mass search
        search_frame = ttk.LabelFrame(parent, text="Extract XIC (auto TOF-MS → find MS/MS)", padding=4)
        search_frame.pack(fill=tk.X, padx=2, pady=(0, 4))

        row1 = ttk.Frame(search_frame)
        row1.pack(fill=tk.X, pady=(0, 2))
        ttk.Label(row1, text="Target m/z:").pack(side=tk.LEFT)
        self._mass_var = tk.StringVar(value="")
        self._mass_entry = ttk.Entry(row1, textvariable=self._mass_var, width=12)
        self._mass_entry.pack(side=tk.LEFT, padx=(4, 8))
        ttk.Label(row1, text="Tol (±ppm):").pack(side=tk.LEFT)
        self._tol_var = tk.StringVar(value="10")
        self._tol_entry = ttk.Entry(row1, textvariable=self._tol_var, width=6)
        self._tol_entry.pack(side=tk.LEFT, padx=(4, 0))

        row2 = ttk.Frame(search_frame)
        row2.pack(fill=tk.X)
        self._extract_btn = ttk.Button(
            row2, text="Extract XIC + Find MS/MS",
            command=self._extract_xic, state=tk.DISABLED,
        )
        self._extract_btn.pack(side=tk.LEFT)
        self._cancel_btn = ttk.Button(
            row2, text="Cancel", command=self._cancel_extraction, state=tk.DISABLED,
        )
        self._cancel_btn.pack(side=tk.LEFT, padx=(4, 0))

        # Progress bar
        self._progress = ttk.Progressbar(search_frame, mode="determinate")
        self._progress.pack(fill=tk.X, pady=(4, 0))

        # MS/MS match list
        msms_frame = ttk.LabelFrame(parent, text="Matching MS/MS Spectra", padding=4)
        msms_frame.pack(fill=tk.BOTH, expand=True, padx=2, pady=(0, 4))

        self._msms_tree = ttk.Treeview(
            msms_frame, columns=("time", "precursor", "tic"),
            show="tree headings", selectmode="browse", height=8,
        )
        self._msms_tree.heading("#0", text="Exp")
        self._msms_tree.heading("time", text="RT (min)")
        self._msms_tree.heading("precursor", text="Precursor m/z")
        self._msms_tree.heading("tic", text="TIC (50→prec)")
        self._msms_tree.column("#0", width=35, anchor=tk.CENTER)
        self._msms_tree.column("time", width=65, anchor=tk.CENTER)
        self._msms_tree.column("precursor", width=100, anchor=tk.CENTER)
        self._msms_tree.column("tic", width=90, anchor=tk.CENTER)
        self._msms_tree.pack(fill=tk.BOTH, expand=True)
        self._msms_tree.bind("<<TreeviewSelect>>", self._on_msms_selected)

        lib_frame = ttk.LabelFrame(parent, text="HighResNPS Library Search", padding=4)
        lib_frame.pack(fill=tk.BOTH, expand=True, padx=2, pady=(0, 4))

        lib_row = ttk.Frame(lib_frame)
        lib_row.pack(fill=tk.X, pady=(0, 4))
        self._library_btn = ttk.Button(
            lib_row, text="Search HighResNPS JSON",
            command=self._search_current_msms_library, state=tk.DISABLED,
        )
        self._library_btn.pack(side=tk.LEFT)

        self._library_tree = ttk.Treeview(
            lib_frame, columns=("score", "matches", "precursor", "formula"),
            show="tree headings", selectmode="browse", height=7,
        )
        self._library_tree.heading("#0", text="Hit")
        self._library_tree.heading("score", text="Score")
        self._library_tree.heading("matches", text="Peaks")
        self._library_tree.heading("precursor", text="Prec.")
        self._library_tree.heading("formula", text="Formula")
        self._library_tree.column("#0", width=130, minwidth=110)
        self._library_tree.column("score", width=55, anchor=tk.CENTER)
        self._library_tree.column("matches", width=45, anchor=tk.CENTER)
        self._library_tree.column("precursor", width=70, anchor=tk.CENTER)
        self._library_tree.column("formula", width=75, anchor=tk.CENTER)
        self._library_tree.pack(fill=tk.BOTH, expand=True)
        self._library_tree.bind("<<TreeviewSelect>>", self._on_library_hit_selected)

    # ── right panel (vertical stack: XIC top, MS/MS bottom) ───────────────
    def _build_right_panel(self, parent: ttk.Frame) -> None:
        # Split the right panel vertically with a PanedWindow
        self._right_pane = ttk.PanedWindow(parent, orient=tk.VERTICAL)
        self._right_pane.pack(fill=tk.BOTH, expand=True)

        # ── XIC + isotope (top, split 70/30) ──
        xic_frame = ttk.Frame(self._right_pane)
        self._right_pane.add(xic_frame, weight=1)

        xic_hpane = ttk.PanedWindow(xic_frame, orient=tk.HORIZONTAL)
        xic_hpane.pack(fill=tk.BOTH, expand=True)

        # XIC plot (70%)
        xic_plot_frame = ttk.Frame(xic_hpane)
        xic_hpane.add(xic_plot_frame, weight=7)

        self._xic_fig = Figure(figsize=(6, 3.5), dpi=100)
        self._xic_ax = self._xic_fig.add_subplot(111)
        self._xic_ax.set_title("Extracted Ion Chromatogram (TOF-MS)")
        _style_dense_axes(self._xic_ax)
        self._xic_fig.subplots_adjust(left=0.055, right=0.985, bottom=0.08, top=0.9)

        self._xic_canvas = FigureCanvasTkAgg(self._xic_fig, xic_plot_frame)
        self._xic_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self._xic_toolbar = NavigationToolbar2Tk(self._xic_canvas, xic_plot_frame)
        self._xic_toolbar.update()

        # Isotope pattern at clicked RT (30%)
        iso_frame = ttk.Frame(xic_hpane, padding=2)
        xic_hpane.add(iso_frame, weight=3)

        self._iso_fig = Figure(figsize=(3, 3.5), dpi=100)
        self._iso_ax = self._iso_fig.add_subplot(111)
        self._iso_ax.set_xlabel("m/z")
        self._iso_ax.set_title("Click XIC to view isotope pattern")
        _style_dense_axes(self._iso_ax)
        self._iso_fig.subplots_adjust(left=0.12, right=0.985, bottom=0.12, top=0.9)

        self._iso_canvas = FigureCanvasTkAgg(self._iso_fig, iso_frame)
        self._iso_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self._iso_toolbar = NavigationToolbar2Tk(self._iso_canvas, iso_frame)
        self._iso_toolbar.update()

        # ── MS/MS (bottom) ──
        msms_frame = ttk.Frame(self._right_pane)
        self._right_pane.add(msms_frame, weight=1)

        self._msms_fig = Figure(figsize=(8, 3.5), dpi=100)
        self._msms_ax = self._msms_fig.add_subplot(111)
        self._msms_ax.set_xlabel("m/z")
        self._msms_ax.set_title("MS/MS Spectrum (select from left panel)")
        _style_dense_axes(self._msms_ax)
        self._msms_fig.subplots_adjust(left=0.045, right=0.992, bottom=0.11, top=0.9)

        self._msms_canvas = FigureCanvasTkAgg(self._msms_fig, msms_frame)
        self._msms_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self._msms_toolbar = NavigationToolbar2Tk(self._msms_canvas, msms_frame)
        self._msms_toolbar.update()

        # Force the sash to split the space evenly (prevents collapsed pane)
        self.update_idletasks()
        self._right_pane.update_idletasks()
        self.after(100, self._fix_pane_sash)

    # ── file open ─────────────────────────────────────────────────────────
    def _open_file(self) -> None:
        path = filedialog.askopenfilename(
            title="Open WIFF2 file",
            filetypes=[("WIFF2 files", "*.wiff2"), ("All files", "*.*")],
        )
        if not path:
            return
        self._load_file(Path(path))

    def _load_file(self, path: Path) -> None:
        self._set_status(f"Loading {path.name}…")
        self._close_reader()
        try:
            self._reader = WiffReader(path)
        except Exception as exc:
            messagebox.showerror("Error", f"Failed to open file:\n{exc}")
            self._set_status("Ready")
            return

        self._current_path = path
        self._file_label_var.set(path.name)
        self._samples = self._reader.list_samples()
        if self._samples:
            self._current_sample_idx = 0
            self._experiments = self._reader.get_experiments(self._current_sample_idx)
            self._extract_btn.config(state=tk.NORMAL if self._experiments else tk.DISABLED)
        else:
            self._experiments = []
            self._extract_btn.config(state=tk.DISABLED)
        self._set_status(f"Opened {path.name}")

    def _close_reader(self) -> None:
        if self._extraction_cancel_event is not None:
            self._extraction_cancel_event.set()
        if self._reader is not None:
            try:
                self._reader.close()
            except Exception:
                pass
            self._reader = None
        self._current_path = None
        self._samples = []
        self._experiments = []
        self._current_xic = None
        self._current_msms_spectrum = None
        self._clear_xic_plot()
        self._clear_msms_plot()
        self._msms_tree.delete(*self._msms_tree.get_children())
        self._library_tree.delete(*self._library_tree.get_children())
        self._library_hits_by_iid = {}
        self._current_library_hit = None
        self._library_btn.config(state=tk.DISABLED)
        self._extract_btn.config(state=tk.DISABLED)
        self._file_label_var.set("No file opened")
        self._current_formula = ""

    def _on_close(self) -> None:
        self._close_reader()
        self.destroy()

    def _get_ms1_experiment(self) -> ExperimentInfo | None:
        """Return the first TOF-MS (MS1) experiment, or None."""
        for exp in self._experiments:
            if exp.ms_level == 1:
                return exp
        return None

    # ── formula calculation ───────────────────────────────────────────────
    def _on_calculate_formula(self) -> None:
        formula = self._formula_var.get().strip()
        if not formula:
            return
        adduct = self._adduct_var.get()
        mass_val = _calculate_formula_mass(formula, adduct)
        if mass_val is None:
            if not _HAS_PYTEOMICS:
                messagebox.showwarning(
                    "pyteomics not installed",
                    "Formula calculation requires the 'pyteomics' package.\n\n"
                    "Install it with:  pip install pyteomics",
                )
            else:
                messagebox.showwarning(
                    "Invalid formula",
                    f"Could not parse formula '{formula}'.\n"
                    "Use standard Hill notation, e.g. C12H17NO2",
                )
            return
        self._mass_var.set(f"{mass_val:.6f}")
        self._current_formula = formula
        self._current_adduct = adduct
        self._set_status(f"Formula: {formula} + {adduct} → m/z {mass_val:.6f}")

    # ── XIC extraction + auto MS/MS matching ──────────────────────────────
    def _extract_xic(self) -> None:
        if self._reader is None or self._current_path is None:
            return

        if self._extraction_thread is not None and self._extraction_thread.is_alive():
            return

        ms1 = self._get_ms1_experiment()
        if ms1 is None:
            messagebox.showwarning("No TOF-MS", "No TOF-MS (MS1) experiment found in this file.")
            return

        try:
            target_mz = float(self._mass_var.get().strip())
        except ValueError:
            messagebox.showwarning("Invalid m/z", "Please enter a valid numeric m/z value.")
            return
        try:
            tol_ppm = float(self._tol_var.get().strip())
        except ValueError:
            messagebox.showwarning("Invalid tolerance", "Please enter a valid numeric ppm tolerance.")
            return

        if target_mz <= 0 or tol_ppm <= 0:
            messagebox.showwarning("Invalid value", "m/z and tolerance must be positive.")
            return

        tol_da = _ppm_to_da(target_mz, tol_ppm)

        result_queue: queue.Queue[tuple[str, Any]] = queue.Queue()
        cancel_event = threading.Event()
        self._extraction_queue = result_queue
        self._extraction_cancel_event = cancel_event
        self._extract_btn.config(state=tk.DISABLED)
        self._cancel_btn.config(state=tk.NORMAL)
        self._progress.config(maximum=100, value=0)
        self._set_status(f"Extracting XIC: m/z {target_mz:.4f} ± {tol_da:.4f} Da…")

        worker_args = (
            self._current_path,
            self._current_sample_idx,
            list(self._experiments),
            ms1.index,
            target_mz,
            tol_da,
            result_queue,
            cancel_event,
        )
        self._extraction_thread = threading.Thread(
            target=self._run_extraction_worker,
            args=worker_args,
            daemon=True,
        )
        self._extraction_thread.start()
        self.after(50, self._poll_extraction_queue)

    @staticmethod
    def _run_extraction_worker(
        path: Path,
        sample_idx: int,
        experiments: list[ExperimentInfo],
        ms1_idx: int,
        target_mz: float,
        tol_da: float,
        result_queue: queue.Queue[tuple[str, Any]],
        cancel_event: threading.Event,
    ) -> None:
        try:
            def _progress_cb(current: int, total: int) -> None:
                result_queue.put(("progress", (current, total)))

            with WiffReader(path) as reader:
                xic = _extract_xic_from_ms1(
                    reader,
                    sample_idx,
                    ms1_idx,
                    target_mz,
                    tol_da,
                    progress_cb=_progress_cb,
                    cancel_event=cancel_event,
                )
                result_queue.put(("status", "Finding matching MS/MS spectra…"))
                msms_matches = _find_matching_msms(
                    reader,
                    sample_idx,
                    experiments,
                    target_mz,
                    tol_da,
                    cancel_event=cancel_event,
                )
            result_queue.put(("result", _ExtractionResult(
                xic=xic,
                msms_matches=msms_matches,
                target_mz=target_mz,
                tol_da=tol_da,
                ms1_idx=ms1_idx,
            )))
        except _ExtractionCancelled:
            result_queue.put(("cancelled", None))
        except Exception as exc:
            result_queue.put(("error", str(exc)))

    def _poll_extraction_queue(self) -> None:
        active_queue = self._extraction_queue
        if active_queue is None:
            return
        keep_polling = True
        while True:
            try:
                kind, payload = active_queue.get_nowait()
            except queue.Empty:
                break
            if kind == "progress":
                current, total = payload
                self._progress.config(maximum=total, value=current)
            elif kind == "status":
                self._set_status(str(payload))
            elif kind == "result":
                self._finish_extraction()
                self._handle_extraction_result(payload)
                keep_polling = False
            elif kind == "cancelled":
                self._finish_extraction()
                self._set_status("Extraction cancelled.")
                keep_polling = False
            elif kind == "error":
                self._finish_extraction()
                messagebox.showerror("Error", f"XIC extraction failed:\n{payload}")
                keep_polling = False

        if keep_polling and self._extraction_queue is active_queue:
            self.after(50, self._poll_extraction_queue)

    def _handle_extraction_result(self, result: _ExtractionResult) -> None:
        # ── Display XIC ──
        if result.xic is None or not result.xic.times:
            messagebox.showinfo("No data", "No intensities found in the specified m/z window.")
            return

        self._current_xic = result.xic
        self._current_target_mz = result.target_mz
        self._current_ms1_idx = result.ms1_idx
        self._plot_xic(result.xic, result.target_mz, result.tol_da)

        # ── Populate MS/MS match list ──
        self._populate_msms_matches(result.msms_matches)

        # ── Auto-select first MS/MS match if any ──
        if result.msms_matches:
            first = self._msms_tree.get_children()[0]
            self._msms_tree.selection_set(first)
            self._on_msms_selected()

        msg = (
            f"XIC: m/z {result.target_mz:.4f} ± {result.tol_da:.4f} Da — "
            f"{len(result.xic.times)} points, max={max(result.xic.intensities):.0f}"
        )
        if result.msms_matches:
            msg += f" | {len(result.msms_matches)} MS/MS match(es)"
        self._set_status(msg)

    def _cancel_extraction(self) -> None:
        if self._extraction_cancel_event is not None:
            self._extraction_cancel_event.set()
        self._cancel_btn.config(state=tk.DISABLED)
        self._set_status("Cancelling extraction…")

    def _finish_extraction(self) -> None:
        self._extract_btn.config(state=tk.NORMAL)
        self._cancel_btn.config(state=tk.DISABLED)
        self._extraction_cancel_event = None
        self._extraction_queue = None
        self._extraction_thread = None

    # ── XIC plotting ──────────────────────────────────────────────────────
    def _format_xic_axis(self, intensities: np.ndarray | None = None) -> None:
        self._xic_ax.set_xlabel("")
        self._xic_ax.set_ylabel("")
        self._xic_ax.xaxis.set_major_locator(MaxNLocator(nbins=9, prune=None))
        if intensities is not None and intensities.size:
            y_max = float(np.nanmax(intensities))
            if y_max > 0:
                y_top = y_max * 1.04
                self._xic_ax.set_ylim(0, y_top)
                self._xic_ax.set_yticks([y_top])
                self._xic_ax.set_yticklabels([_compact_sci(y_max)])
            else:
                self._xic_ax.set_ylim(0, 1)
                self._xic_ax.set_yticks([1])
                self._xic_ax.set_yticklabels(["0"])
        _style_dense_axes(self._xic_ax)

    def _plot_xic(self, xic: Chromatogram, target_mz: float, tol_da: float) -> None:
        self._xic_ax.clear()
        t = np.array(xic.times)
        y = np.array(xic.intensities)
        self._xic_ax.plot(t, y, linewidth=0.8, color="#1f77b4", picker=True, pickradius=3)
        self._xic_ax.set_title(
            f"XIC: m/z {target_mz:.4f} ± {tol_da:.4f} Da (TOF-MS)"
        )
        self._xic_ax.relim()
        self._xic_ax.autoscale_view()
        self._format_xic_axis(y)
        self._xic_fig.subplots_adjust(left=0.055, right=0.985, bottom=0.08, top=0.9)

        # Connect click handler (disconnect previous if any)
        if self._xic_click_cid is not None:
            self._xic_fig.canvas.mpl_disconnect(self._xic_click_cid)
        self._xic_click_cid = self._xic_fig.canvas.mpl_connect(
            "button_press_event", self._on_xic_click,
        )

        self._xic_canvas.draw_idle()

    def _clear_xic_plot(self) -> None:
        self._xic_ax.clear()
        self._xic_ax.set_title("Extracted Ion Chromatogram (TOF-MS)")
        self._format_xic_axis()
        self._xic_fig.subplots_adjust(left=0.055, right=0.985, bottom=0.08, top=0.9)
        self._xic_canvas.draw_idle()
        self._clear_isotope_plot()

    # ── XIC click → isotope pattern ──────────────────────────────────────
    def _on_xic_click(self, event: Any) -> None:
        """Handle click on XIC: read MS1 spectrum at nearest RT and show isotope pattern."""
        if event.inaxes != self._xic_ax:
            return
        if self._reader is None or self._current_xic is None:
            return
        if not self._current_xic.times:
            return

        # Find nearest retention time
        times_arr = np.array(self._current_xic.times)
        idx = int(np.argmin(np.abs(times_arr - event.xdata)))
        rt_clicked = times_arr[idx]
        cycle = idx  # XIC index maps 1:1 to cycle index

        self._set_status(f"Loading spectrum at RT={rt_clicked:.2f} min (cycle {cycle})…")

        try:
            spec = self._reader.get_spectrum(
                self._current_sample_idx, self._current_ms1_idx, cycle,
                centroid=True, return_arrays=True,
            )
        except Exception as exc:
            self._set_status(f"Failed to load spectrum: {exc}")
            return

        self._plot_isotope(spec, rt_clicked)
        self._set_status(
            f"Isotope pattern at RT={rt_clicked:.2f} min (cycle {cycle}), "
            f"{len(spec.mz)} points"
        )

    def _clear_isotope_labels(self) -> None:
        _remove_artists(self._iso_label_artists)
        self._iso_label_artists = []

    def _refresh_isotope_labels(self, _ax: Any = None) -> None:
        if self._iso_label_data is None:
            return
        mz, rel = self._iso_label_data
        if mz.size == 0:
            return

        self._clear_isotope_labels()
        self._iso_label_artists = _place_peak_labels(
            self._iso_ax,
            mz,
            rel,
            label_for=lambda idx: f"{mz[idx]:.4f}",
            threshold_pct=2.0,
            max_labels=18,
            color="#165c26",
            rotation=90,
        )
        self._iso_canvas.draw_idle()

    def _plot_isotope(self, spec: SpectrumData, rt: float) -> None:
        """Plot isotope pattern: m/z in [target-0.5, target+2.5], relative abundance."""
        self._iso_ax.clear()
        self._clear_isotope_labels()
        self._iso_label_data = None

        target = self._current_target_mz
        mz_lo = target - 0.5
        mz_hi = target + 2.5

        mz_all = np.asarray(spec.mz)
        int_all = np.asarray(spec.intensities)
        mask = (mz_all >= mz_lo) & (mz_all <= mz_hi)
        mz = mz_all[mask]
        intensities = int_all[mask]

        if len(mz) == 0:
            self._iso_ax.set_title(f"No peaks in [{mz_lo:.2f}, {mz_hi:.2f}]")
            self._iso_ax.set_xlabel("m/z")
            self._iso_ax.set_ylabel("")
            _style_dense_axes(self._iso_ax)
            self._iso_fig.subplots_adjust(left=0.12, right=0.985, bottom=0.12, top=0.9)
            self._iso_canvas.draw_idle()
            return

        # Normalize to max = 100%
        max_int = float(intensities.max())
        if max_int > 0:
            rel = intensities / max_int * 100.0
        else:
            rel = np.zeros_like(intensities)

        self._iso_ax.vlines(mz, 0, rel, color="#2ca02c", linewidth=0.8)

        self._iso_ax.set_xlabel("m/z")
        self._iso_ax.set_ylabel("")
        self._iso_ax.set_title(f"Isotope pattern @ RT={rt:.2f} min")
        self._iso_ax.set_xlim(mz_lo, mz_hi)
        self._iso_ax.set_ylim(0, 105)
        self._iso_ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _pos: f"{v:.0f}%"))
        self._iso_ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        _style_dense_axes(self._iso_ax)
        self._iso_fig.subplots_adjust(left=0.12, right=0.985, bottom=0.12, top=0.9)

        # ── Overlay expected isotope distribution from formula ─────────
        if self._current_formula:
            expected = _compute_isotopic_distribution(
                self._current_formula, self._current_adduct,
            )
            for exp_mz, exp_abund in expected:
                if exp_abund < 0.05:
                    continue
                half_width_da = exp_mz * 10e-6  # ± 10 ppm half-width → 20 ppm total
                # Rectangle height matches expected abundance (ymax in axes coords 0-1)
                self._iso_ax.axvspan(
                    exp_mz - half_width_da, exp_mz + half_width_da,
                    ymin=0, ymax=exp_abund / 105.0,
                    alpha=0.25, color="#d62728", zorder=0,
                )
                # Label the monoisotopic peak
                if exp_abund > 99.0:
                    self._iso_ax.annotate(
                        f"{exp_mz:.4f}",
                        xy=(exp_mz, 102), ha="center", va="bottom",
                        fontsize=7, color="#d62728",
                    )

        self._iso_label_data = (mz, rel)
        if self._iso_xlim_cid is None:
            self._iso_xlim_cid = self._iso_ax.callbacks.connect(
                "xlim_changed", self._refresh_isotope_labels,
            )
        if self._iso_ylim_cid is None:
            self._iso_ylim_cid = self._iso_ax.callbacks.connect(
                "ylim_changed", self._refresh_isotope_labels,
            )
        self._refresh_isotope_labels()
        self._iso_canvas.draw_idle()

    def _clear_isotope_plot(self) -> None:
        self._iso_ax.clear()
        self._clear_isotope_labels()
        self._iso_label_data = None
        self._iso_ax.set_xlabel("m/z")
        self._iso_ax.set_ylabel("")
        self._iso_ax.set_title("Select MS/MS to view isotope pattern")
        _style_dense_axes(self._iso_ax)
        self._iso_fig.subplots_adjust(left=0.12, right=0.985, bottom=0.12, top=0.9)
        self._iso_canvas.draw_idle()

    # ── MS/MS matching ────────────────────────────────────────────────────
    def _populate_msms_matches(self, matches: list[_MsMsMatch]) -> None:
        self._msms_tree.delete(*self._msms_tree.get_children())
        for m in matches:
            iid = f"{m.experiment_index}:{m.cycle_index}"
            tic_str = f"{m.tic_50_parent:.0f}" if m.tic_50_parent > 0 else "0"
            self._msms_tree.insert(
                "", tk.END, iid=iid,
                text=str(m.experiment_index),
                values=(f"{m.scan_time:.2f}", f"{m.precursor_mz:.4f}", tic_str),
            )

    def _on_msms_selected(self, _event: Any = None) -> None:
        sel = self._msms_tree.selection()
        if not sel or self._reader is None:
            return
        iid = sel[0]
        exp_idx_str, cycle_str = iid.split(":")
        exp_idx = int(exp_idx_str)
        cycle = int(cycle_str)

        try:
            spec = self._reader.get_spectrum(
                self._current_sample_idx, exp_idx, cycle,
                centroid=True, return_arrays=True,
            )
        except Exception as exc:
            messagebox.showerror("Error", f"Failed to load spectrum:\n{exc}")
            return

        self._current_msms_spectrum = spec
        self._current_library_hit = None
        self._plot_msms(spec)
        self._library_tree.delete(*self._library_tree.get_children())
        self._library_hits_by_iid = {}
        self._library_btn.config(state=tk.NORMAL)

        # ── Auto-show isotope pattern from MS1 at this RT ────────────
        if self._current_target_mz > 0:
            self._show_isotope_at_rt(spec.scan_time)

        self._set_status(
            f"MS/MS: exp={exp_idx}, cycle={cycle}, "
            f"RT={spec.scan_time:.2f} min, "
            f"precursor={spec.precursor_mz or '?'}, {len(spec.mz)} points"
        )

    def _search_current_msms_library(self) -> None:
        if self._current_msms_spectrum is None:
            messagebox.showinfo("No MS/MS selected", "Select a matching MS/MS spectrum first.")
            return
        if self._library_thread is not None and self._library_thread.is_alive():
            return
        if not _DEFAULT_CSV_LIBRARY.exists():
            messagebox.showerror("Library not found", f"CSV file not found:\n{_DEFAULT_CSV_LIBRARY}")
            return

        spec = self._current_msms_spectrum
        query_mz = np.asarray(spec.mz, dtype=np.float64)
        query_intensity = np.asarray(spec.intensities, dtype=np.float64)
        precursor_mz = float(spec.precursor_mz or 0.0)
        result_queue: queue.Queue[tuple[str, Any]] = queue.Queue()
        self._library_queue = result_queue
        self._library_btn.config(state=tk.DISABLED)
        self._library_tree.delete(*self._library_tree.get_children())
        self._set_status("Searching HighResNPS MatchMS JSON library...")

        def _worker() -> None:
            try:
                _ensure_matchms_json_library(_DEFAULT_CSV_LIBRARY, _DEFAULT_JSON_LIBRARY)
                if self._json_library_cache is None:
                    library = _load_matchms_json_library(_DEFAULT_JSON_LIBRARY)
                    self._json_library_cache = library
                else:
                    library = self._json_library_cache
                hits = _search_json_library_with_matchms(
                    query_mz,
                    query_intensity,
                    library,
                    precursor_mz=precursor_mz or None,
                    tolerance_da=0.02,
                    top_n=10,
                )
                result_queue.put(("result", hits))
            except Exception as exc:
                result_queue.put(("error", str(exc)))

        self._library_thread = threading.Thread(target=_worker, daemon=True)
        self._library_thread.start()
        self.after(50, self._poll_library_queue)

    def _poll_library_queue(self) -> None:
        active_queue = self._library_queue
        if active_queue is None:
            return
        keep_polling = True
        while True:
            try:
                kind, payload = active_queue.get_nowait()
            except queue.Empty:
                break
            if kind == "result":
                self._finish_library_search()
                self._populate_library_hits(payload)
                keep_polling = False
            elif kind == "error":
                self._finish_library_search()
                messagebox.showerror("Library search failed", str(payload))
                keep_polling = False

        if keep_polling and self._library_queue is active_queue:
            self.after(50, self._poll_library_queue)

    def _finish_library_search(self) -> None:
        self._library_queue = None
        self._library_thread = None
        self._library_btn.config(state=tk.NORMAL if self._current_msms_spectrum is not None else tk.DISABLED)

    def _populate_library_hits(self, hits: list[_LibraryHit]) -> None:
        self._library_tree.delete(*self._library_tree.get_children())
        self._library_hits_by_iid = {}
        self._current_library_hit = None
        if not hits:
            self._set_status("HighResNPS library search: no hits.")
            return
        for hit in hits:
            display = f"{hit.rank}. {hit.name}"
            iid = f"hit:{hit.rank}"
            self._library_hits_by_iid[iid] = hit
            self._library_tree.insert(
                "", tk.END, iid=iid, text=display,
                values=(
                    f"{hit.score:.3f}",
                    str(hit.matches),
                    f"{hit.precursor_mz:.4f}" if hit.precursor_mz else "",
                    hit.formula,
                ),
            )
        best = hits[0]
        self._set_status(
            f"HighResNPS library search: best {best.name} "
            f"(score={best.score:.3f}, matched peaks={best.matches}, precursor={best.precursor_mz:.4f})"
        )

    def _on_library_hit_selected(self, _event: Any = None) -> None:
        sel = self._library_tree.selection()
        if not sel or self._current_msms_spectrum is None:
            return
        hit = self._library_hits_by_iid.get(sel[0])
        if hit is None:
            return
        self._current_library_hit = hit
        self._plot_msms(self._current_msms_spectrum, library_hit=hit)
        self._set_status(
            f"Library comparison: {hit.name} "
            f"(score={hit.score:.3f}, matched peaks={hit.matches})"
        )

    def _show_isotope_at_rt(self, rt: float) -> None:
        """Load MS1 spectrum at the given retention time and display isotope pattern."""
        if self._reader is None or self._current_xic is None:
            return
        if not self._current_xic.times:
            return
        # Find nearest RT in the XIC → cycle index
        times_arr = np.array(self._current_xic.times)
        idx = int(np.argmin(np.abs(times_arr - rt)))
        cycle = idx  # XIC index maps 1:1 to cycle index
        try:
            spec = self._reader.get_spectrum(
                self._current_sample_idx, self._current_ms1_idx, cycle,
                centroid=True, return_arrays=True,
            )
        except Exception:
            return
        self._plot_isotope(spec, self._current_xic.times[idx])

    # ── MS/MS plotting ────────────────────────────────────────────────────
    def _clear_msms_labels(self) -> None:
        _remove_artists(self._msms_label_artists)
        self._msms_label_artists = []

    def _refresh_msms_labels(self, _ax: Any = None) -> None:
        if self._msms_label_data is None:
            return
        mz, intensities, rel = self._msms_label_data
        if mz.size == 0:
            return

        self._clear_msms_labels()
        self._msms_label_artists = _place_peak_labels(
            self._msms_ax,
            mz,
            rel,
            y_values=intensities,
            label_for=lambda idx: f"[{mz[idx]:.4f}, {rel[idx]:.0f}%]",
            threshold_pct=2.0,
            max_labels=10,
            color="#8c1d1d",
            rotation=90,
        )
        self._msms_canvas.draw_idle()

    def _plot_msms(self, spec: SpectrumData, library_hit: _LibraryHit | None = None) -> None:
        self._msms_ax.clear()
        self._clear_msms_labels()
        self._msms_label_data = None
        mz_all = np.asarray(spec.mz)
        intensities_all = np.asarray(spec.intensities)

        # Trim: only show m/z from 50 up to precursor+0.5
        precursor = spec.precursor_mz or 0.0
        mz_hi_data = max(precursor + 0.5, 50.0)
        mask = (mz_all >= 50.0) & (mz_all <= mz_hi_data)
        mz = mz_all[mask]
        intensities = intensities_all[mask]

        if len(mz) > 0:
            self._msms_ax.vlines(mz, 0, intensities, color="#d62728", linewidth=0.7)

        lib_mz = np.array([], dtype=np.float64)
        lib_intensities = np.array([], dtype=np.float64)
        if library_hit is not None and library_hit.mz.size and library_hit.intensity.size:
            lib_mask = (library_hit.mz >= 50.0) & (library_hit.mz <= mz_hi_data)
            lib_mz = library_hit.mz[lib_mask]
            lib_intensities = library_hit.intensity[lib_mask]
            if lib_intensities.size:
                lib_max = float(np.nanmax(lib_intensities))
                if lib_max > 0:
                    scale_to = float(np.nanmax(intensities)) if len(intensities) else lib_max
                    lib_plot = -(lib_intensities / lib_max * scale_to)
                    self._msms_ax.vlines(
                        lib_mz, 0, lib_plot,
                        color="#1f77b4", linewidth=0.8, alpha=0.85,
                    )

        self._msms_ax.set_xlabel("m/z")
        self._msms_ax.set_ylabel("")

        # Fixed x-axis for uniform display: 50 to precursor+10
        x_hi = max(precursor + 10.0, 60.0)
        self._msms_ax.set_xlim(50.0, x_hi)
        y_max = 0.0
        if len(intensities) > 0:
            y_max = float(np.nanmax(intensities))
            if y_max > 0:
                y_bottom = -y_max * 1.12 if library_hit is not None and lib_intensities.size else 0
                self._msms_ax.set_ylim(y_bottom, y_max * 1.18)
                rel_pct = intensities / y_max * 100.0
                self._msms_label_data = (mz, intensities, rel_pct)

        title = f"MS/MS Spectrum — RT={spec.scan_time:.2f} min"
        if spec.precursor_mz:
            title += f", precursor m/z={spec.precursor_mz:.4f}"
        if library_hit is not None:
            title += f" | library: {library_hit.name[:36]}"
        self._msms_ax.set_title(title)
        if library_hit is not None and lib_intensities.size and y_max > 0:
            self._msms_ax.axhline(0, color="#444444", linewidth=0.6)
        self._msms_ax.xaxis.set_major_locator(MaxNLocator(nbins=12))
        self._msms_ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _pos: _compact_sci(abs(v))))
        _style_dense_axes(self._msms_ax)
        self._msms_fig.subplots_adjust(left=0.045, right=0.992, bottom=0.11, top=0.9)
        if self._msms_xlim_cid is None:
            self._msms_xlim_cid = self._msms_ax.callbacks.connect(
                "xlim_changed", self._refresh_msms_labels,
            )
        if self._msms_ylim_cid is None:
            self._msms_ylim_cid = self._msms_ax.callbacks.connect(
                "ylim_changed", self._refresh_msms_labels,
            )
        self._refresh_msms_labels()
        self._msms_canvas.draw_idle()

    def _fix_pane_sash(self) -> None:
        """Ensure the vertical PanedWindow sash splits space evenly."""
        try:
            height = self._right_pane.winfo_height()
            if height > 50:
                self._right_pane.sashpos(0, height // 2)
        except Exception:
            pass

    def _clear_msms_plot(self) -> None:
        self._msms_ax.clear()
        self._clear_msms_labels()
        self._msms_label_data = None
        self._msms_ax.set_xlabel("m/z")
        self._msms_ax.set_ylabel("")
        self._msms_ax.set_title("MS/MS Spectrum (select from left panel)")
        _style_dense_axes(self._msms_ax)
        self._msms_fig.subplots_adjust(left=0.045, right=0.992, bottom=0.11, top=0.9)
        self._msms_canvas.draw_idle()

    # ── utils ─────────────────────────────────────────────────────────────
    def _set_status(self, text: str) -> None:
        self._status_var.set(text)

    def _show_about(self) -> None:
        messagebox.showinfo(
            "About pyx500r WIFF2 Viewer",
            "pyx500r WIFF2 Viewer\n\n"
            "Pure-Python SCIEX X500R acquisition browser.\n"
            "https://github.com/ykcchong/pyx500r",
        )


# ── entry point ────────────────────────────────────────────────────────────


def main() -> None:
    """Launch the WIFF2 GUI application."""
    app = WiffGuiApp()
    app.mainloop()


if __name__ == "__main__":
    main()
