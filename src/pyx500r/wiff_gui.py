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

    pyx500r-gui
    # or
    python -m pyx500r.wiff_gui

Dependencies: ``matplotlib``.  Optional: ``pyteomics`` for formula mass calculation.
"""

from __future__ import annotations

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


def _extract_xic_from_ms1(
    reader: WiffReader,
    sample_index: int,
    ms1_experiment_index: int,
    target_mz: float,
    mz_tolerance_da: float,
    rt_start: float | None = None,
    rt_end: float | None = None,
    progress_cb: Callable[[int, int], None] | None = None,
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
        if mz_arr.size == 0:
            total_int = 0.0
        else:
            total_int = float(int_arr[(mz_arr >= mz_lo) & (mz_arr <= mz_hi)].sum())
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


def _find_matching_msms(
    reader: WiffReader,
    sample_index: int,
    experiments: list[ExperimentInfo],
    target_mz: float,
    tol_da: float,
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
            times = reader.get_cycle_times(sample_index, exp.index)
        except Exception:
            continue
        for ci in range(min(len(times), exp.cycle_count)):
            try:
                meta = reader.get_spectrum_metadata(sample_index, exp.index, ci)
            except Exception:
                continue
            if meta.precursor_mz is not None and mz_lo <= meta.precursor_mz <= mz_hi:
                # Compute TIC in [50, precursor_mz] range
                tic = 0.0
                try:
                    spec = reader.get_spectrum(
                        sample_index, exp.index, ci,
                        centroid=False, return_arrays=True,
                    )
                    mz_arr = np.asarray(spec.mz)
                    int_arr = np.asarray(spec.intensities)
                    mask = (mz_arr >= 50.0) & (mz_arr <= meta.precursor_mz)
                    tic = float(int_arr[mask].sum()) if mask.any() else 0.0
                except Exception:
                    pass
                matches.append(_MsMsMatch(
                    experiment_index=exp.index,
                    cycle_index=ci,
                    scan_time=meta.scan_time,
                    precursor_mz=meta.precursor_mz,
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
        self._extraction_cancelled = False

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
        # Sample selector
        sample_frame = ttk.LabelFrame(parent, text="Sample", padding=4)
        sample_frame.pack(fill=tk.X, padx=2, pady=(2, 4))

        self._sample_var = tk.StringVar(value="(no file)")
        self._sample_combo = ttk.Combobox(
            sample_frame, textvariable=self._sample_var, state="readonly",
        )
        self._sample_combo.pack(fill=tk.X)
        self._sample_combo.bind("<<ComboboxSelected>>", self._on_sample_changed)

        # Experiment tree (read-only info)
        exp_frame = ttk.LabelFrame(parent, text="Experiments", padding=4)
        exp_frame.pack(fill=tk.BOTH, expand=True, padx=2, pady=(0, 4))

        self._exp_tree = ttk.Treeview(
            exp_frame, columns=("cycles", "polarity"),
            show="tree headings", selectmode="browse", height=6,
        )
        self._exp_tree.heading("#0", text="Experiment")
        self._exp_tree.heading("cycles", text="Cycles")
        self._exp_tree.heading("polarity", text="Polarity")
        self._exp_tree.column("#0", width=140, minwidth=80)
        self._exp_tree.column("cycles", width=50, anchor=tk.CENTER)
        self._exp_tree.column("polarity", width=60, anchor=tk.CENTER)
        self._exp_tree.pack(fill=tk.BOTH, expand=True)

        # Formula → mass calculator
        formula_frame = ttk.LabelFrame(parent, text="Formula → m/z (optional, needs pyteomics)", padding=4)
        formula_frame.pack(fill=tk.X, padx=2, pady=(0, 4))

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
        self._xic_ax.set_xlabel("Retention Time (min)")
        self._xic_ax.set_ylabel("Intensity")
        self._xic_ax.set_title("Extracted Ion Chromatogram (TOF-MS)")
        self._xic_fig.tight_layout()

        self._xic_canvas = FigureCanvasTkAgg(self._xic_fig, xic_plot_frame)
        self._xic_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self._xic_toolbar = NavigationToolbar2Tk(self._xic_canvas, xic_plot_frame)
        self._xic_toolbar.update()

        # Isotope pattern at clicked RT (30%)
        iso_frame = ttk.LabelFrame(xic_hpane, text="Isotope @ clicked RT", padding=2)
        xic_hpane.add(iso_frame, weight=3)

        self._iso_fig = Figure(figsize=(3, 3.5), dpi=100)
        self._iso_ax = self._iso_fig.add_subplot(111)
        self._iso_ax.set_xlabel("m/z")
        self._iso_ax.set_ylabel("Rel. Abundance (%)")
        self._iso_ax.set_title("Click XIC to view isotope pattern")
        self._iso_fig.tight_layout()

        self._iso_canvas = FigureCanvasTkAgg(self._iso_fig, iso_frame)
        self._iso_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # ── MS/MS (bottom) ──
        msms_frame = ttk.Frame(self._right_pane)
        self._right_pane.add(msms_frame, weight=1)

        self._msms_fig = Figure(figsize=(8, 3.5), dpi=100)
        self._msms_ax = self._msms_fig.add_subplot(111)
        self._msms_ax.set_xlabel("m/z")
        self._msms_ax.set_ylabel("Intensity")
        self._msms_ax.set_title("MS/MS Spectrum (select from left panel)")
        self._msms_fig.tight_layout()

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

        self._file_label_var.set(path.name)
        self._samples = self._reader.list_samples()
        self._populate_samples()
        self._set_status(f"Opened {path.name} — {len(self._samples)} sample(s)")

    def _close_reader(self) -> None:
        if self._reader is not None:
            try:
                self._reader.close()
            except Exception:
                pass
            self._reader = None
        self._samples = []
        self._experiments = []
        self._current_xic = None
        self._current_msms_spectrum = None
        self._clear_xic_plot()
        self._clear_msms_plot()
        self._sample_var.set("(no file)")
        self._sample_combo["values"] = []
        self._exp_tree.delete(*self._exp_tree.get_children())
        self._msms_tree.delete(*self._msms_tree.get_children())
        self._extract_btn.config(state=tk.DISABLED)
        self._file_label_var.set("No file opened")
        self._current_formula = ""

    def _on_close(self) -> None:
        self._close_reader()
        self.destroy()

    # ── sample selection ──────────────────────────────────────────────────
    def _populate_samples(self) -> None:
        names = [s.name or f"Sample {s.index}" for s in self._samples]
        self._sample_combo["values"] = names
        if names:
            self._sample_combo.current(0)
            self._on_sample_changed()
        else:
            self._sample_var.set("(no samples)")

    def _on_sample_changed(self, _event: Any = None) -> None:
        idx = self._sample_combo.current()
        if idx < 0 or self._reader is None:
            return
        self._current_sample_idx = idx
        self._experiments = self._reader.get_experiments(idx)
        self._populate_experiments()

        # Enable extract button as soon as a file + sample are available
        if self._experiments:
            self._extract_btn.config(state=tk.NORMAL)

    # ── experiment tree (read-only info) ───────────────────────────────────
    def _populate_experiments(self) -> None:
        self._exp_tree.delete(*self._exp_tree.get_children())
        for exp in self._experiments:
            ms_label = f"MS{exp.ms_level}" if exp.ms_level else "?"
            display = f"{exp.scan_type or 'Exp ' + str(exp.index)} ({ms_label})"
            self._exp_tree.insert(
                "", tk.END,
                iid=str(exp.index),
                text=display,
                values=(str(exp.cycle_count), exp.polarity),
            )

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
        if self._reader is None:
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

        self._extraction_cancelled = False
        self._extract_btn.config(state=tk.DISABLED)
        self._cancel_btn.config(state=tk.NORMAL)
        self._progress.config(value=0)
        self._set_status(f"Extracting XIC: m/z {target_mz:.4f} ± {tol_da:.4f} Da…")

        def _progress_cb(current: int, total: int) -> None:
            if self._extraction_cancelled:
                return
            self._progress.config(maximum=total, value=current)
            self.update_idletasks()

        self.after(50, lambda: self._do_extract_and_match(
            ms1.index, target_mz, tol_da, _progress_cb,
        ))

    def _do_extract_and_match(
        self, ms1_idx: int, target_mz: float, tol_da: float,
        progress_cb: Callable[[int, int], None],
    ) -> None:
        # ── Step 1: extract XIC ──
        try:
            xic = _extract_xic_from_ms1(
                self._reader,  # type: ignore[arg-type]
                self._current_sample_idx,
                ms1_idx,
                target_mz,
                tol_da,
                progress_cb=progress_cb,
            )
        except Exception as exc:
            self._finish_extraction()
            messagebox.showerror("Error", f"XIC extraction failed:\n{exc}")
            return

        # ── Step 2: find matching MS/MS spectra ──
        msms_matches: list[_MsMsMatch] = []
        try:
            msms_matches = _find_matching_msms(
                self._reader,  # type: ignore[arg-type]
                self._current_sample_idx,
                self._experiments,
                target_mz,
                tol_da,
            )
        except Exception:
            pass  # non-fatal — still show the XIC

        self._finish_extraction()

        # ── Display XIC ──
        if xic is None or not xic.times:
            messagebox.showinfo("No data", "No intensities found in the specified m/z window.")
            return

        self._current_xic = xic
        self._current_target_mz = target_mz
        self._current_ms1_idx = ms1_idx
        self._plot_xic(xic, target_mz, tol_da)

        # ── Populate MS/MS match list ──
        self._populate_msms_matches(msms_matches)

        # ── Auto-select first MS/MS match if any ──
        if msms_matches:
            first = self._msms_tree.get_children()[0]
            self._msms_tree.selection_set(first)
            self._on_msms_selected()

        msg = (
            f"XIC: m/z {target_mz:.4f} ± {tol_da:.4f} Da — "
            f"{len(xic.times)} points, max={max(xic.intensities):.0f}"
        )
        if msms_matches:
            msg += f" | {len(msms_matches)} MS/MS match(es)"
        self._set_status(msg)

    def _cancel_extraction(self) -> None:
        self._extraction_cancelled = True
        self._finish_extraction()
        self._set_status("Extraction cancelled.")

    def _finish_extraction(self) -> None:
        self._extract_btn.config(state=tk.NORMAL)
        self._cancel_btn.config(state=tk.DISABLED)

    # ── XIC plotting ──────────────────────────────────────────────────────
    def _plot_xic(self, xic: Chromatogram, target_mz: float, tol_da: float) -> None:
        self._xic_ax.clear()
        t = np.array(xic.times)
        y = np.array(xic.intensities)
        self._xic_ax.plot(t, y, linewidth=0.8, color="#1f77b4", picker=True, pickradius=3)
        self._xic_ax.set_xlabel("Retention Time (min)")
        self._xic_ax.set_ylabel("Intensity")
        self._xic_ax.set_title(
            f"XIC: m/z {target_mz:.4f} ± {tol_da:.4f} Da (TOF-MS)"
        )
        self._xic_ax.relim()
        self._xic_ax.autoscale_view()
        self._xic_fig.tight_layout()

        # Connect click handler (disconnect previous if any)
        if self._xic_click_cid is not None:
            self._xic_fig.canvas.mpl_disconnect(self._xic_click_cid)
        self._xic_click_cid = self._xic_fig.canvas.mpl_connect(
            "button_press_event", self._on_xic_click,
        )

        self._xic_canvas.draw_idle()

    def _clear_xic_plot(self) -> None:
        self._xic_ax.clear()
        self._xic_ax.set_xlabel("Retention Time (min)")
        self._xic_ax.set_ylabel("Intensity")
        self._xic_ax.set_title("Extracted Ion Chromatogram (TOF-MS)")
        self._xic_fig.tight_layout()
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
                centroid=False, return_arrays=True,
            )
        except Exception as exc:
            self._set_status(f"Failed to load spectrum: {exc}")
            return

        self._plot_isotope(spec, rt_clicked)
        self._set_status(
            f"Isotope pattern at RT={rt_clicked:.2f} min (cycle {cycle}), "
            f"{len(spec.mz)} points"
        )

    def _plot_isotope(self, spec: SpectrumData, rt: float) -> None:
        """Plot isotope pattern: m/z in [target-0.5, target+2.5], relative abundance."""
        self._iso_ax.clear()

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
            self._iso_ax.set_ylabel("Rel. Abundance (%)")
            self._iso_fig.tight_layout()
            self._iso_canvas.draw_idle()
            return

        # Normalize to max = 100%
        max_int = float(intensities.max())
        if max_int > 0:
            rel = intensities / max_int * 100.0
        else:
            rel = np.zeros_like(intensities)

        # Bar chart
        if len(mz) > 1:
            bar_width = max((mz[-1] - mz[0]) * 0.001, 0.0005)
        else:
            bar_width = 0.01
        self._iso_ax.bar(mz, rel, width=bar_width, color="#2ca02c",
                         edgecolor="#2ca02c", linewidth=0)

        self._iso_ax.set_xlabel("m/z")
        self._iso_ax.set_ylabel("Rel. Abundance (%)")
        self._iso_ax.set_title(f"Isotope pattern @ RT={rt:.2f} min")
        self._iso_ax.set_xlim(mz_lo, mz_hi)
        self._iso_ax.set_ylim(0, 105)
        self._iso_fig.tight_layout()

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

        self._iso_canvas.draw_idle()

    def _clear_isotope_plot(self) -> None:
        self._iso_ax.clear()
        self._iso_ax.set_xlabel("m/z")
        self._iso_ax.set_ylabel("Rel. Abundance (%)")
        self._iso_ax.set_title("Select MS/MS to view isotope pattern")
        self._iso_fig.tight_layout()
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
                centroid=False, return_arrays=True,
            )
        except Exception as exc:
            messagebox.showerror("Error", f"Failed to load spectrum:\n{exc}")
            return

        self._current_msms_spectrum = spec
        self._plot_msms(spec)

        # ── Auto-show isotope pattern from MS1 at this RT ────────────
        if self._current_target_mz > 0:
            self._show_isotope_at_rt(spec.scan_time)

        self._set_status(
            f"MS/MS: exp={exp_idx}, cycle={cycle}, "
            f"RT={spec.scan_time:.2f} min, "
            f"precursor={spec.precursor_mz or '?'}, {len(spec.mz)} points"
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
                centroid=False, return_arrays=True,
            )
        except Exception:
            return
        self._plot_isotope(spec, self._current_xic.times[idx])

    # ── MS/MS plotting ────────────────────────────────────────────────────
    def _plot_msms(self, spec: SpectrumData) -> None:
        self._msms_ax.clear()
        mz_all = np.asarray(spec.mz)
        intensities_all = np.asarray(spec.intensities)

        # Trim: only show m/z from 50 up to precursor+0.5
        precursor = spec.precursor_mz or 0.0
        mz_hi_data = max(precursor + 0.5, 50.0)
        mask = (mz_all >= 50.0) & (mz_all <= mz_hi_data)
        mz = mz_all[mask]
        intensities = intensities_all[mask]

        if len(mz) > 0:
            # Bar chart (stick spectrum).  Use a small relative bar width.
            if len(mz) > 1:
                bar_width = max((mz[-1] - mz[0]) * 0.001, 0.0001)
            else:
                bar_width = 0.01
            self._msms_ax.bar(mz, intensities, width=bar_width, color="#d62728",
                              edgecolor="#d62728", linewidth=0)

        self._msms_ax.set_xlabel("m/z")
        self._msms_ax.set_ylabel("Intensity")

        # Fixed x-axis for uniform display: 50 to precursor+10
        x_hi = max(precursor + 10.0, 60.0)
        self._msms_ax.set_xlim(50.0, x_hi)

        title = f"MS/MS Spectrum — RT={spec.scan_time:.2f} min"
        if spec.precursor_mz:
            title += f", precursor m/z={spec.precursor_mz:.4f}"
        self._msms_ax.set_title(title)
        self._msms_ax.relim()
        self._msms_ax.autoscale_view(scaley=True)
        self._msms_fig.tight_layout()
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
        self._msms_ax.set_xlabel("m/z")
        self._msms_ax.set_ylabel("Intensity")
        self._msms_ax.set_title("MS/MS Spectrum (select from left panel)")
        self._msms_fig.tight_layout()
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
