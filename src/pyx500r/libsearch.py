"""Library search identification using the LibraryView SQLite database.

Two-tier search strategy:

1. **Fast pre-screening** — 16-mass ``MassSpectrumSignature`` fingerprint
   with precursor m/z tolerance.
2. **Full spectral matching** — dot-product similarity on centroided peak
   lists extracted from ``CentroidedXYData`` blobs.

Usage::

    from pyx500r.libsearch import LibrarySearcher

    searcher = LibrarySearcher("data/libview.sqlite")
    results = searcher.search(mz_array, intensity_array, precursor_mz, polarity="POS")
    for r in results[:5]:
        print(f"{r['name']} | score={r['score']:.3f} | prec={r['precursor_mz']:.4f}")
"""
from __future__ import annotations

import sqlite3
import struct
from pathlib import Path
from typing import Any, Optional

import numpy as np

# ── helpers ───────────────────────────────────────────────────────

def _extract_double_arrays_from_blob(blob: bytes) -> tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Extract _xValues (m/z) and _yValues (intensity) directly from a
    BinaryFormatter blob by scanning for ``ArraySinglePrimitive`` (0x0F)
    records.

    Returns ``(mz_array, intensity_array)`` as float64 numpy arrays,
    or ``(None, None)`` on failure.
    """
    mz: Optional[np.ndarray] = None
    intensity: Optional[np.ndarray] = None

    pos = 0
    while pos < len(blob) - 15:
        if blob[pos] == 0x0F:  # ArraySinglePrimitive
            try:
                obj_id = struct.unpack_from("<i", blob, pos + 1)[0]
                arr_len = struct.unpack_from("<i", blob, pos + 5)[0]
                prim_type = blob[pos + 9]

                if prim_type == 6 and 0 < arr_len < 1_000_000:  # double
                    vals = np.frombuffer(
                        blob, dtype=np.float64, count=arr_len, offset=pos + 10
                    )
                    if obj_id == 3 and mz is None:
                        mz = vals
                    elif obj_id == 4 and intensity is None:
                        intensity = vals

                pos += 10 + arr_len * 8
                continue
            except (struct.error, ValueError):
                pass
        pos += 1

    return mz, intensity


def _dot_product(
    query_mz: np.ndarray,
    query_intensity: np.ndarray,
    ref_mz: np.ndarray,
    ref_intensity: np.ndarray,
    ppm_tol: float = 20.0,
) -> float:
    """Compute normalised dot-product similarity between two spectra.

    Uses bidirectional peak matching: for each peak in both spectra,
    the closest matching peak in the other spectrum (within ``ppm_tol``)
    contributes.  The geometric mean of forward and reverse scores is
    returned, bounded to [0, 1].

    .. math::

        \\text{forward}  = \\frac{\\sum q_i \\cdot r_{\\text{match}(i)}}
                                {\\|q\\| \\cdot \\|r\\|}

        \\text{reverse}  = \\frac{\\sum r_j \\cdot q_{\\text{match}(j)}}
                                {\\|q\\| \\cdot \\|r\\|}

        \\text{score} = \\sqrt{\\text{forward} \\cdot \\text{reverse}}
    """
    if len(query_mz) == 0 or len(ref_mz) == 0:
        return 0.0

    q_norm = float(np.linalg.norm(query_intensity))
    r_norm = float(np.linalg.norm(ref_intensity))
    if q_norm == 0 or r_norm == 0:
        return 0.0

    q_int_n = query_intensity / q_norm
    r_int_n = ref_intensity / r_norm

    norm_factor = q_norm * r_norm

    def _match(a_mz, a_int, b_mz, b_int, tol_da):
        """Sum of a_int[i] * best_matching_b_intensity for each a peak."""
        total = 0.0
        indices = np.searchsorted(b_mz, a_mz)
        for i, (mz, a_i) in enumerate(zip(a_mz, a_int)):
            lo = max(0, indices[i] - 2)
            hi = min(len(b_mz), indices[i] + 3)
            best = 0.0
            td = tol_da[i] if hasattr(tol_da, '__len__') else tol_da
            for j in range(lo, hi):
                if abs(b_mz[j] - mz) <= td:
                    if b_int[j] > best:
                        best = b_int[j]
            total += a_i * best
        return total

    # Forward: query → reference
    q_tol = query_mz * ppm_tol * 1e-6
    fwd = _match(query_mz, q_int_n * q_norm, ref_mz, r_int_n, q_tol) / norm_factor

    # Reverse: reference → query
    r_tol = ref_mz * ppm_tol * 1e-6
    rev = _match(ref_mz, r_int_n * r_norm, query_mz, q_int_n, r_tol) / norm_factor

    score = float(np.sqrt(max(fwd, 0.0) * max(rev, 0.0)))
    return min(score, 1.0)


def _signature_similarity(
    query_sig: np.ndarray, ref_sig: np.ndarray, tol_da: float = 0.05
) -> float:
    """Fast pre-screening using 16-mass signatures.

    Counts how many signature masses match within absolute tolerance.
    Returns fraction of query masses that match (0.0 – 1.0).
    """
    if len(query_sig) == 0 or len(ref_sig) == 0:
        return 0.0

    matches = 0
    # Only compare non-zero masses
    q_mask = query_sig > 0
    r_mask = ref_sig > 0

    q_nonzero = query_sig[q_mask]
    r_nonzero = ref_sig[r_mask]

    if len(q_nonzero) == 0:
        return 0.0

    for qm in q_nonzero:
        # Check if any reference mass is within tolerance
        if np.any(np.abs(r_nonzero - qm) <= tol_da):
            matches += 1

    return matches / len(q_nonzero)


# ── main search class ─────────────────────────────────────────────

class LibrarySearcher:
    """Search a LibraryView SQLite database for matching MS/MS spectra.

    Parameters
    ----------
    db_path:
        Path to the ``.sqlite`` database (converted from LBP SqlCE).
    """

    def __init__(self, db_path: str | Path) -> None:
        self._conn = sqlite3.connect(str(db_path))
        self._conn.row_factory = sqlite3.Row
        self._ensure_indexes()

    def _ensure_indexes(self) -> None:
        """Create helper indexes if not present."""
        self._conn.executescript("""
            CREATE INDEX IF NOT EXISTS idx_libsearch_prec
                ON MassSpectrum(PrecursorMass1, PositivePolarity);
            CREATE INDEX IF NOT EXISTS idx_libsearch_sig
                ON MassSpectrumSignature(MassSpectrumID);
        """)

    def close(self) -> None:
        self._conn.close()

    def __enter__(self) -> "LibrarySearcher":
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()

    # ── pre-screening ──────────────────────────────────────────

    def prescreen(
        self,
        query_signature: np.ndarray,
        precursor_mz: float,
        polarity: str = "POS",
        ppm_tol: float = 50.0,
        sig_tol_da: float = 0.05,
        top_n: int = 200,
    ) -> list[sqlite3.Row]:
        """Fast pre-screening using mass signatures and precursor m/z.

        Parameters
        ----------
        query_signature:
            16-element array of signature masses (like in
            ``MassSpectrumSignature``). Non-zero entries only.
        precursor_mz:
            Precursor m/z of the query spectrum.
        polarity:
            ``"POS"`` or ``"NEG"``.
        ppm_tol:
            Precursor m/z tolerance in ppm.
        sig_tol_da:
            Signature mass matching tolerance in Da.
        top_n:
            Maximum number of candidates to return.

        Returns
        -------
        List of sqlite3.Row objects ordered by signature similarity
        (best first).  Each row has columns from ``MassSpectrum``,
        ``MassSpectrumSignature``, ``Compound``, and ``CompoundName``.
        """
        pol = 1 if polarity.upper() == "POS" else 0
        mz_lo = precursor_mz * (1 - ppm_tol * 1e-6)
        mz_hi = precursor_mz * (1 + ppm_tol * 1e-6)

        # Fetch candidates within precursor m/z window
        cur = self._conn.execute(
            """
            SELECT ms.*, sig.*, cp.Identifier, cp.Formula,
                   cp.MolecularWeight, cp.MonoIsotopicMass, cp.CAS,
                   cn.Name as CompoundName
            FROM MassSpectrum ms
            JOIN MassSpectrumSignature sig ON ms.Id = sig.MassSpectrumID
            JOIN Compound cp ON ms.CompoundId = cp.Id
            JOIN CompoundName cn ON cp.Id = cn.CompoundId AND cn.IsDefault = 1
            WHERE ms.PrecursorMass1 BETWEEN ? AND ?
              AND ms.PositivePolarity = ?
            """,
            (mz_lo, mz_hi, pol),
        )

        candidates: list[tuple[float, sqlite3.Row]] = []

        for row in cur:
            # Build reference signature
            ref_sig = np.array(
                [row[f"Mass{i}"] for i in range(16)], dtype=np.float64
            )
            sim = _signature_similarity(query_signature, ref_sig, sig_tol_da)
            if sim > 0:
                candidates.append((sim, row))

        # Sort by similarity descending, take top_n
        candidates.sort(key=lambda x: x[0], reverse=True)
        return [row for _, row in candidates[:top_n]]

    # ── full spectral search ───────────────────────────────────

    def search(
        self,
        query_mz: np.ndarray,
        query_intensity: np.ndarray,
        precursor_mz: float,
        polarity: str = "POS",
        ppm_tol: float = 50.0,
        sig_tol_da: float = 0.05,
        dot_product_ppm: float = 20.0,
        prescreen_n: int = 200,
        top_n: int = 10,
    ) -> list[dict[str, Any]]:
        """Search the library for matching compounds.

        Two-tier search:
        1. Pre-screen using mass signatures + precursor m/z
        2. Full dot-product similarity on centroided spectra

        Parameters
        ----------
        query_mz, query_intensity:
            Centroided peak list from the unknown spectrum.
        precursor_mz:
            Precursor m/z of the unknown.
        polarity:
            ``"POS"`` or ``"NEG"``.
        ppm_tol:
            Precursor m/z window for pre-screening (ppm).
        sig_tol_da:
            Signature mass matching tolerance (Da).
        dot_product_ppm:
            PPM tolerance for dot-product peak matching.
        prescreen_n:
            Max candidates from pre-screening.
        top_n:
            Number of results to return.

        Returns
        -------
        List of dicts with keys: ``name``, ``formula``, ``cas``,
        ``molecular_weight``, ``precursor_mz``, ``collision_energy``,
        ``score`` (0–1), ``num_peaks``, ``spectrum_id``, ``compound_id``.
        """
        # Build query signature from the query spectrum
        query_sig = self._build_signature(query_mz, query_intensity)

        # Pre-screen
        candidates = self.prescreen(
            query_sig, precursor_mz, polarity, ppm_tol, sig_tol_da, prescreen_n
        )

        # Full spectral matching
        results: list[dict[str, Any]] = []

        for row in candidates:
            blob = row["CentroidedXYData"]
            if blob is None:
                continue

            try:
                ref_mz, ref_intensity = _extract_double_arrays_from_blob(blob)
            except Exception:
                continue

            if ref_mz is None or ref_intensity is None:
                continue
            if len(ref_mz) == 0 or len(ref_intensity) == 0:
                continue

            score = _dot_product(
                query_mz, query_intensity, ref_mz, ref_intensity, dot_product_ppm
            )

            results.append(
                {
                    "name": row["CompoundName"] or row["Identifier"] or "",
                    "formula": row["Formula"] or "",
                    "cas": row["CAS"] or "",
                    "molecular_weight": row["MolecularWeight"] or 0.0,
                    "monoisotopic_mass": row["MonoIsotopicMass"] or 0.0,
                    "precursor_mz": row["PrecursorMass1"],
                    "collision_energy": row["CollisionEnergy"],
                    "score": score,
                    "num_peaks": len(ref_mz),
                    "spectrum_id": row["Id"],
                    "compound_id": row["CompoundId"],
                }
            )

        results.sort(key=lambda x: x["score"], reverse=True)
        return results[:top_n]

    # ── signature builder ──────────────────────────────────────

    @staticmethod
    def _build_signature(
        mz: np.ndarray, intensity: np.ndarray, n: int = 16
    ) -> np.ndarray:
        """Build a 16-mass signature from a centroided spectrum.

        Selects the top *n* most intense peaks and returns their m/z
        values sorted by mass.  Pads with zeros if fewer than *n* peaks.
        """
        if len(mz) == 0:
            return np.zeros(n, dtype=np.float64)

        # Take top N by intensity
        order = np.argsort(intensity)[::-1][:n]
        top_mz = mz[order]
        top_mz = np.sort(top_mz)

        sig = np.zeros(n, dtype=np.float64)
        sig[: len(top_mz)] = top_mz
        return sig

    # ── database statistics ────────────────────────────────────

    @property
    def stats(self) -> dict[str, Any]:
        """Return summary statistics about the library."""
        return {
            "compounds": self._conn.execute(
                "SELECT COUNT(*) FROM Compound"
            ).fetchone()[0],
            "spectra": self._conn.execute(
                "SELECT COUNT(*) FROM MassSpectrum"
            ).fetchone()[0],
            "positive": self._conn.execute(
                "SELECT COUNT(*) FROM MassSpectrum WHERE PositivePolarity=1"
            ).fetchone()[0],
            "negative": self._conn.execute(
                "SELECT COUNT(*) FROM MassSpectrum WHERE PositivePolarity=0"
            ).fetchone()[0],
            "precursor_range": self._conn.execute(
                "SELECT MIN(PrecursorMass1), MAX(PrecursorMass1) FROM MassSpectrum"
            ).fetchone(),
            "libraries": [
                row["Name"]
                for row in self._conn.execute("SELECT Name FROM Library")
            ],
        }

    def get_spectrum(self, spectrum_id: str) -> Optional[dict[str, Any]]:
        """Retrieve a single library spectrum by ID.

        Returns a dict with ``mz``, ``intensity``, and metadata, or
        ``None`` if not found.
        """
        row = self._conn.execute(
            """
            SELECT ms.*, cp.Identifier, cp.Formula, cp.CAS,
                   cp.MolecularWeight, cp.MonoIsotopicMass,
                   cn.Name as CompoundName
            FROM MassSpectrum ms
            JOIN Compound cp ON ms.CompoundId = cp.Id
            JOIN CompoundName cn ON cp.Id = cn.CompoundId AND cn.IsDefault = 1
            WHERE ms.Id = ?
            """,
            (spectrum_id,),
        ).fetchone()

        if row is None:
            return None

        blob = row["CentroidedXYData"]
        mz, intensity = None, None
        if blob:
            mz, intensity = _extract_double_arrays_from_blob(blob)

        return {
            "spectrum_id": row["Id"],
            "compound_id": row["CompoundId"],
            "name": row["CompoundName"] or row["Identifier"] or "",
            "formula": row["Formula"] or "",
            "cas": row["CAS"] or "",
            "molecular_weight": row["MolecularWeight"],
            "monoisotopic_mass": row["MonoIsotopicMass"],
            "precursor_mz": row["PrecursorMass1"],
            "collision_energy": row["CollisionEnergy"],
            "polarity": "POS" if row["PositivePolarity"] else "NEG",
            "mz": mz,
            "intensity": intensity,
        }
