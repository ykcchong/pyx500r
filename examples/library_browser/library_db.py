"""Read-only data access layer over a SCIEX LibraryView ``.sqlite`` database.

This wraps the converted LibraryView schema (see docs/LBP_FORMAT.md) and exposes
a small, GUI-friendly API:

* :meth:`LibraryDB.libraries`            — the library folders/collections
* :meth:`LibraryDB.list_compounds`       — paged/filtered compound list
* :meth:`LibraryDB.get_compound`         — full compound "settings" detail
* :meth:`LibraryDB.list_spectra`         — reference-spectrum metadata for a compound
* :meth:`LibraryDB.get_spectrum`         — one spectrum's metadata + decoded peaks

Spectral peak arrays are decoded from the ``CentroidedXYData`` / ``RawXYData``
BinaryFormatter blobs using :func:`pyx500r.libsearch._extract_double_arrays_from_blob`.

All connections are read-only and opened per-thread (see :class:`LibraryDB`).
"""

from __future__ import annotations

import sqlite3
import threading
from pathlib import Path
from typing import Any

from pyx500r.libsearch import _extract_double_arrays_from_blob


class LibraryDB:
    """Thread-safe read-only accessor for a LibraryView ``.sqlite`` file.

    A separate SQLite connection is created per thread (cached in thread-local
    storage), so a single ``LibraryDB`` instance can be shared across a FastAPI
    thread pool safely.
    """

    def __init__(self, db_path: str | Path) -> None:
        self.db_path = Path(db_path).resolve()
        if not self.db_path.exists():
            raise FileNotFoundError(f"library DB not found: {self.db_path}")
        self._local = threading.local()
        # Validate schema up front + warm a connection on the main thread.
        conn = self._conn()
        names = {r[0] for r in conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        )}
        missing = {"Compound", "CompoundName", "MassSpectrum", "Library"} - names
        if missing:
            raise ValueError(f"not a LibraryView database (missing tables: {missing})")

    # ------------------------------------------------------------------ #
    # connection handling
    # ------------------------------------------------------------------ #
    def _conn(self) -> sqlite3.Connection:
        conn = getattr(self._local, "conn", None)
        if conn is None:
            # read-only URI connection; check_same_thread=False is safe because
            # each thread gets its *own* connection via thread-local storage.
            conn = sqlite3.connect(
                f"file:{self.db_path}?mode=ro", uri=True, check_same_thread=False
            )
            conn.row_factory = sqlite3.Row
            self._local.conn = conn
        return conn

    def close(self) -> None:
        conn = getattr(self._local, "conn", None)
        if conn is not None:
            conn.close()
            self._local.conn = None

    # ------------------------------------------------------------------ #
    # libraries
    # ------------------------------------------------------------------ #
    def libraries(self) -> list[dict[str, Any]]:
        """Return all libraries (collections) with a compound count each."""
        rows = self._conn().execute(
            """
            SELECT l.Id, l.Name, l.IsFolder,
                   (SELECT COUNT(*) FROM CompoundLibrary cl WHERE cl.LibraryId = l.Id) AS compound_count
            FROM Library l
            ORDER BY l.Name
            """
        ).fetchall()
        return [
            {
                "id": r["Id"],
                "name": r["Name"],
                "is_folder": bool(r["IsFolder"]),
                "compound_count": r["compound_count"],
            }
            for r in rows
        ]

    def stats(self) -> dict[str, Any]:
        c = self._conn()
        return {
            "compounds": c.execute("SELECT COUNT(*) FROM Compound").fetchone()[0],
            "spectra": c.execute("SELECT COUNT(*) FROM MassSpectrum").fetchone()[0],
            "libraries": c.execute("SELECT COUNT(*) FROM Library").fetchone()[0],
            "db_path": str(self.db_path),
        }

    # ------------------------------------------------------------------ #
    # compounds
    # ------------------------------------------------------------------ #
    def list_compounds(
        self,
        *,
        search: str | None = None,
        library_id: str | None = None,
        offset: int = 0,
        limit: int = 100,
    ) -> dict[str, Any]:
        """Return a paged compound list, optionally filtered.

        Parameters
        ----------
        search:
            Case-insensitive substring matched against the default compound
            name, the identifier, the formula, and the CAS number.
        library_id:
            Restrict to compounds belonging to this library.
        offset, limit:
            Pagination window (``limit`` is capped at 500 by the caller).

        Returns ``{"total": int, "offset": int, "limit": int, "items": [...]}``.
        """
        where: list[str] = []
        params: list[Any] = []

        join = ""
        if library_id:
            join = "JOIN CompoundLibrary cl ON cl.CompoundId = c.Id AND cl.LibraryId = ?"
            params.append(library_id)

        if search:
            like = f"%{search.strip()}%"
            where.append(
                "(cn.Name LIKE ? OR c.Identifier LIKE ? OR c.Formula LIKE ? OR c.CAS LIKE ?)"
            )
            params.extend([like, like, like, like])

        where_sql = ("WHERE " + " AND ".join(where)) if where else ""

        base = f"""
            FROM Compound c
            LEFT JOIN CompoundName cn ON cn.CompoundId = c.Id AND cn.IsDefault = 1
            {join}
            {where_sql}
        """

        total = self._conn().execute(
            f"SELECT COUNT(*) {base}", params
        ).fetchone()[0]

        rows = self._conn().execute(
            f"""
            SELECT c.Id, cn.Name AS name, c.Identifier, c.Formula, c.CAS,
                   c.MolecularWeight, c.MonoIsotopicMass,
                   (SELECT COUNT(*) FROM MassSpectrum ms WHERE ms.CompoundId = c.Id) AS spectrum_count
            {base}
            ORDER BY cn.Name COLLATE NOCASE
            LIMIT ? OFFSET ?
            """,
            [*params, limit, offset],
        ).fetchall()

        items = [
            {
                "id": r["Id"],
                "name": r["name"] or r["Identifier"] or "(unnamed)",
                "identifier": r["Identifier"],
                "formula": r["Formula"],
                "cas": r["CAS"],
                "molecular_weight": r["MolecularWeight"],
                "monoisotopic_mass": r["MonoIsotopicMass"],
                "spectrum_count": r["spectrum_count"],
            }
            for r in rows
        ]
        return {"total": total, "offset": offset, "limit": limit, "items": items}

    def get_compound(self, compound_id: str) -> dict[str, Any] | None:
        """Return full compound "settings": identity, thresholds, all names,
        library memberships, retention times, and spectrum metadata."""
        c = self._conn()
        row = c.execute(
            "SELECT * FROM Compound WHERE Id = ?", (compound_id,)
        ).fetchone()
        if row is None:
            return None

        names = c.execute(
            """
            SELECT cn.Name, cn.IsDefault, r.Name AS region
            FROM CompoundName cn
            LEFT JOIN Region r ON r.Id = cn.RegionId
            WHERE cn.CompoundId = ?
            ORDER BY cn.IsDefault DESC, cn.Name
            """,
            (compound_id,),
        ).fetchall()

        libraries = c.execute(
            """
            SELECT l.Id, l.Name
            FROM CompoundLibrary cl JOIN Library l ON l.Id = cl.LibraryId
            WHERE cl.CompoundId = ?
            ORDER BY l.Name
            """,
            (compound_id,),
        ).fetchall()

        retention_times = c.execute(
            """
            SELECT rt.Value, rt.IsDefault, it.ModelName AS instrument
            FROM RetentionTime rt
            LEFT JOIN Instrument i ON i.Id = rt.InstrumentId
            LEFT JOIN InstrumentType it ON it.Id = i.InstrumentTypeId
            WHERE rt.CompoundId = ?
            ORDER BY rt.IsDefault DESC
            """,
            (compound_id,),
        ).fetchall()

        default_name = next(
            (n["Name"] for n in names if n["IsDefault"]), row["Identifier"]
        )

        return {
            "id": row["Id"],
            "name": default_name or "(unnamed)",
            "identifier": row["Identifier"],
            "formula": row["Formula"],
            "cas": row["CAS"],
            "molecular_weight": row["MolecularWeight"],
            "monoisotopic_mass": row["MonoIsotopicMass"],
            "molecular_structure_source": row["MolecularStructureSource"],
            "purity_threshold": row["PurityThreshold"],
            "red_flag_threshold": row["RedFlagThreshold"],
            "yellow_flag_threshold": row["YellowFlagThreshold"],
            "comment": row["Comment"],
            "created_date": row["CreatedDate"],
            "last_updated": row["LastUpdated"],
            "active": bool(row["Active"]) if row["Active"] is not None else None,
            "names": [
                {"name": n["Name"], "is_default": bool(n["IsDefault"]), "region": n["region"]}
                for n in names
            ],
            "libraries": [{"id": l["Id"], "name": l["Name"]} for l in libraries],
            "retention_times": [
                {"value": rt["Value"], "is_default": bool(rt["IsDefault"]),
                 "instrument": rt["instrument"]}
                for rt in retention_times
            ],
            "spectra": self.list_spectra(compound_id),
        }

    # ------------------------------------------------------------------ #
    # spectra
    # ------------------------------------------------------------------ #
    def list_spectra(self, compound_id: str) -> list[dict[str, Any]]:
        """Return metadata for every reference spectrum of a compound.

        A compound may have many spectra (different polarities, collision
        energies, instruments). Peak arrays are *not* decoded here — call
        :meth:`get_spectrum` for that.
        """
        rows = self._conn().execute(
            """
            SELECT ms.Id, ms.PrecursorMass1, ms.PrecursorMass2,
                   ms.PrecursorChargeState1, ms.PositivePolarity,
                   ms.CollisionEnergy, ms.CollisionEnergySpread,
                   ms.Type, ms.IonSource, ms.CADGasType, ms.Encryption,
                   ms.StartRT, ms.EndRT, ms.CreatedDate,
                   ms.RawXYData IS NOT NULL AS has_raw,
                   ms.CentroidedXYData IS NOT NULL AS has_centroid,
                   it.ModelName AS instrument,
                   st.ScanTypeName AS scan_type
            FROM MassSpectrum ms
            LEFT JOIN Instrument i ON i.Id = ms.InstrumentId
            LEFT JOIN InstrumentType it ON it.Id = i.InstrumentTypeId
            LEFT JOIN ScanType st ON st.Id = ms.ScanTypeId
            WHERE ms.CompoundId = ?
            ORDER BY ms.PositivePolarity DESC, ms.CollisionEnergy
            """,
            (compound_id,),
        ).fetchall()
        return [self._spectrum_meta(r) for r in rows]

    @staticmethod
    def _spectrum_meta(r: sqlite3.Row) -> dict[str, Any]:
        keys = r.keys()
        # ``has_raw``/``has_centroid`` are present in the list query; in the
        # detail query (SELECT ms.*) derive them from the blob columns instead.
        if "has_raw" in keys:
            has_raw = bool(r["has_raw"])
            has_centroid = bool(r["has_centroid"])
        else:
            has_raw = r["RawXYData"] is not None
            has_centroid = r["CentroidedXYData"] is not None
        # Some legacy spectra store their XY blobs under an extra encryption
        # layer (Encryption like "HRAIO|2.0") that this reader cannot decode.
        encryption = (r["Encryption"] if "Encryption" in keys else "") or ""
        return {
            "id": r["Id"],
            "precursor_mz": r["PrecursorMass1"],
            "precursor_mz2": r["PrecursorMass2"],
            "charge_state": r["PrecursorChargeState1"],
            "polarity": "POS" if r["PositivePolarity"] else "NEG",
            "collision_energy": r["CollisionEnergy"],
            "collision_energy_spread": r["CollisionEnergySpread"],
            "type": r["Type"] or None,
            "ion_source": r["IonSource"] or None,
            "cad_gas_type": r["CADGasType"] or None,
            "scan_type": r["scan_type"],
            "instrument": r["instrument"],
            "start_rt": r["StartRT"],
            "end_rt": r["EndRT"],
            "created_date": r["CreatedDate"],
            "has_raw": has_raw,
            "has_centroid": has_centroid,
            "encrypted": bool(encryption.strip()),
        }

    def get_spectrum(
        self, spectrum_id: str, kind: str = "centroid"
    ) -> dict[str, Any] | None:
        """Return one spectrum's metadata plus decoded peak arrays.

        Parameters
        ----------
        kind:
            ``"centroid"`` (default) decodes ``CentroidedXYData``;
            ``"raw"`` decodes ``RawXYData`` (profile). Falls back to the other
            if the requested blob is absent.

        Returns a dict with ``mz``, ``intensity`` (parallel ``list[float]``),
        the resolved ``kind``, and all spectrum metadata, or ``None`` if the
        spectrum id is unknown.
        """
        if kind not in ("centroid", "raw"):
            raise ValueError("kind must be 'centroid' or 'raw'")

        c = self._conn()
        r = c.execute(
            """
            SELECT ms.*, it.ModelName AS instrument, st.ScanTypeName AS scan_type,
                   cn.Name AS compound_name
            FROM MassSpectrum ms
            LEFT JOIN Instrument i ON i.Id = ms.InstrumentId
            LEFT JOIN InstrumentType it ON it.Id = i.InstrumentTypeId
            LEFT JOIN ScanType st ON st.Id = ms.ScanTypeId
            LEFT JOIN CompoundName cn ON cn.CompoundId = ms.CompoundId AND cn.IsDefault = 1
            WHERE ms.Id = ?
            """,
            (spectrum_id,),
        ).fetchone()
        if r is None:
            return None

        # Pick the blob, falling back to whichever is present.
        order = ["CentroidedXYData", "RawXYData"] if kind == "centroid" else ["RawXYData", "CentroidedXYData"]
        blob = None
        resolved = kind
        for col in order:
            if r[col] is not None:
                blob = r[col]
                resolved = "centroid" if col == "CentroidedXYData" else "raw"
                break

        mz_list: list[float] = []
        int_list: list[float] = []
        if blob is not None and len(blob) >= 16:
            try:
                mz, inten = _extract_double_arrays_from_blob(bytes(blob))
            except Exception:
                # Malformed / unexpectedly-encrypted blob — return no peaks
                # rather than failing the whole request.
                mz, inten = None, None
            if mz is not None and inten is not None:
                mz_list = [float(x) for x in mz]
                int_list = [float(x) for x in inten]

        meta = self._spectrum_meta(r)
        meta.update(
            {
                "compound_id": r["CompoundId"],
                "compound_name": r["compound_name"],
                "kind": resolved,
                "num_peaks": len(mz_list),
                "base_peak_intensity": max(int_list) if int_list else 0.0,
                "mz": mz_list,
                "intensity": int_list,
            }
        )
        return meta
