"""Bridge between ``.wiff2`` acquisitions and ``.qsession`` results."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .models import UnifiedPeak, XicChromatogram
from .qsession import QSessionReader
from .reader import WiffReader


@dataclass
class ExtractionWindow:
    period: int
    experiment: int
    mz_center: float
    mz_half_window: float
    rt_start: float | None = None
    rt_end: float | None = None


class WiffQSessionBridge:
    #: Number of library hits resolved during :meth:`open` (``None`` until opened).
    library_hits_resolved: int | None

    def __init__(
        self, qsession: str | Path, wiffs: list[str | Path],
        *, match_by: str = "name",
        library_db: str | Path | None = None,
    ):
        if match_by not in ("name", "position"):
            raise ValueError(f"match_by must be 'name' or 'position', got {match_by!r}")
        qsession = Path(qsession).resolve()
        if qsession.suffix.lower() != ".qsession":
            raise ValueError(f"Expected .qsession, got: {qsession.name}")
        self._qsession_path = qsession
        self._wiff_paths = [Path(w).resolve() for w in wiffs]
        for p in self._wiff_paths:
            if p.suffix.lower() != ".wiff2":
                raise ValueError(f"Expected .wiff2, got: {p.name}")
        self._match_by = match_by
        self._qsession: QSessionReader | None = None
        self._wiffs: list[WiffReader] | None = None
        self._sample_map: dict[int, tuple[int, int]] | None = None
        self._library_db: Path | None = Path(library_db).resolve() if library_db else None
        self.library_hits_resolved = None

    def __enter__(self): self.open(); return self
    def __exit__(self, *a): self.close()

    def open(self):
        self._qsession = QSessionReader(self._qsession_path)
        self._wiffs = [WiffReader(p) for p in self._wiff_paths]
        self._sample_map = self._build_sample_map()
        if self._library_db is not None and self._library_db.exists():
            self.library_hits_resolved = self._qsession.resolve_library_hits(self._library_db)

    def close(self):
        if self._qsession: self._qsession.close(); self._qsession = None
        if self._wiffs:
            for w in self._wiffs: w.close()
            self._wiffs = None
        self._sample_map = None

    @property
    def _qs(self) -> QSessionReader:
        if self._qsession is None: raise RuntimeError("Bridge not opened.")
        return self._qsession

    @property
    def _w(self) -> list[WiffReader]:
        if self._wiffs is None: raise RuntimeError("Bridge not opened.")
        return self._wiffs

    def _build_sample_map(self) -> dict[int, tuple[int, int]]:
        result: dict[int, tuple[int, int]] = {}
        qs_samples = self._qs.list_samples()
        if self._match_by == "name":
            name_lookup: dict[str, tuple[int, int]] = {}
            for wi, wiff in enumerate(self._wiffs):
                for s in wiff.list_samples():
                    key = (s.name or "").lower()
                    if key: name_lookup[key] = (wi, s.index)
            for qi, qs in enumerate(qs_samples):
                name = (qs.sample_name or "").lower()
                match = name_lookup.get(name)
                if match is not None: result[qi] = match
        else:
            for qi, qs in enumerate(qs_samples):
                if qi < len(self._wiffs):
                    wiff = self._wiffs[qi]
                    samples = wiff.list_samples()
                    if samples: result[qi] = (qi, samples[0].index)
        return result

    def _resolve_wiff(self, qs_sample_index: int) -> tuple[WiffReader, int] | None:
        if self._sample_map is None: return None
        entry = self._sample_map.get(qs_sample_index)
        if entry is None: return None
        wi, si = entry
        return self._wiffs[wi], si

    def match_samples(self) -> list[dict[str, Any]]:
        result: list[dict[str, Any]] = []
        for qi, qs in enumerate(self._qs.list_samples()):
            resolved = self._sample_map.get(qi) if self._sample_map else None
            wiff_idx, wiff_sample = None, None
            if resolved is not None:
                wiff_idx, si = resolved
                wiff_sample = self._wiffs[wiff_idx].list_samples()[si]
            result.append({"qsession_index": qi, "qsession_sample": qs,
                           "wiff_index": wiff_idx, "wiff_sample": wiff_sample})
        return result

    def get_extraction_window(self, compound_index: int) -> ExtractionWindow | None:
        compounds = self._qs.list_compounds()
        if compound_index < 0 or compound_index >= len(compounds):
            raise IndexError(f"Compound index {compound_index} out of range")
        c = compounds[compound_index]
        ev1 = c.extraction_values1 or []
        ev2 = c.extraction_values2 or []
        if not ev1 or not ev2: return None
        return ExtractionWindow(
            period=c.period, experiment=c.experiment,
            mz_center=float(ev1[0]), mz_half_window=float(ev2[0] - ev1[0]) / 2.0,
        )

    def extract_xic(self, sample_index: int, compound_index: int,
                    rt_start: float | None = None, rt_end: float | None = None,
                    ) -> XicChromatogram | None:
        window = self.get_extraction_window(compound_index)
        if window is None: return None
        resolved = self._resolve_wiff(sample_index)
        if resolved is None: return None
        wiff, wiff_sample_idx = resolved
        try:
            cycle_times = wiff.get_cycle_times(wiff_sample_idx, window.experiment)
        except Exception:
            return None
        if not cycle_times: return None
        mz_lo = window.mz_center - window.mz_half_window
        mz_hi = window.mz_center + window.mz_half_window
        sample_key = ""
        try:
            sample_key = self._qs.list_samples()[sample_index].sample_name
        except (IndexError, RuntimeError):
            pass
        times, intensities = [], []
        for cycle in range(len(cycle_times)):
            t = cycle_times[cycle]
            if rt_start is not None and t < rt_start: continue
            if rt_end is not None and t > rt_end: continue
            try:
                spec = wiff.get_spectrum(wiff_sample_idx, window.experiment, cycle,
                                         centroid=False, return_arrays=True)
            except Exception:
                continue
            if spec is None or spec.mz is None or spec.intensities is None: continue
            mz_arr = np.asarray(spec.mz)
            int_arr = np.asarray(spec.intensities)
            if mz_arr.size == 0:
                total = 0.0
            else:
                total = float(int_arr[(mz_arr >= mz_lo) & (mz_arr <= mz_hi)].sum())
            times.append(t); intensities.append(total)
        if not times: return None
        return XicChromatogram(
            xic_id=f"extract:{sample_index}:{compound_index}",
            sample_key=sample_key,
            times=times,
            intensities=intensities,
            mz_lower=mz_lo,
            mz_upper=mz_hi,
            status=0,
        )

    # ── unified access ──────────────────────────────────────────────

    def unified_results(self) -> list[list[UnifiedPeak]]:
        """Return the matrix with each cell as a ``UnifiedPeak``.

        ``UnifiedPeak`` merges peak + compound + XIC fields into one
        flat namespace::

            up = bridge.unified_results()[0][42]
            print(up.name, up.area, up.retention_time, up.found_mass)
            for hit in up.library_hits:
                print(hit.fit, hit.reverse_fit, hit.purity)
        """
        matrix = self.results_matrix()
        compounds = self.compounds
        result: list[list[UnifiedPeak]] = []
        for si, row in enumerate(matrix):
            new_row: list[UnifiedPeak] = []
            for ci, peak in enumerate(row):
                c = compounds[ci] if ci < len(compounds) else None
                new_row.append(UnifiedPeak(peak, c))
            result.append(new_row)
        return result

    # ── pass-through ────────────────────────────────────────────────

    @property
    def compounds(self): return self._qs.list_compounds()
    @property
    def samples(self): return self._qs.list_samples()
    @property
    def peaks(self): return list(self._qs.iter_peaks())
    def get_peak(self, si, ci): return self._qs.get_peak(si, ci)
    def results_matrix(self): return self._qs.results_matrix()
    def get_chromatogram(self, si, ci): return self._qs.get_chromatogram(si, ci)
    @property
    def qsession(self): return self._qs
    @property
    def wiffs(self): return self._w
