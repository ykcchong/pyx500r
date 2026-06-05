"""Pure-Python reader for SCIEX MultiQuant qsession (quantitation session) files.

A ``.qsession`` file is an encrypted SQLite database produced by SCIEX
MultiQuant software after processing X500R QTOF acquisition data. It stores
extracted ion chromatograms (XICs), peak integration results, compound
definitions, audit trails, and quantitation metadata for small-molecule
and toxicology panels.

The encryption scheme is the same AES-128-OFB SEE cipher used for
``.wiff2`` files, but with a different password and page size:

* Page size: **1024** bytes (wiff2 uses 4096)
* Password: ``"PQS1 is not Sirius"``

Usage::

    from pyx500r import open_qsession

    with open_qsession("quant_results.qsession") as qs:
        print(qs.version)
        for xic in qs.iter_xics():
            print(f"{xic.xic_id}: {len(xic.times)} points")
"""

from __future__ import annotations

import re
import sqlite3
import struct
from collections.abc import Iterator, Sequence
from pathlib import Path

import numpy as np

from .crypto import QSESSION_PAGE_SIZE, QSESSION_PASSWORD, decrypt_database
from .models import (
    CompoundInfo,
    QuantPeakInfo,
    QuantSampleInfo,
    XicChromatogram,
    XicInfo,
)
from .rtparts import load_rtparts_stream, read_compounds, read_multidata

_XIC_ID_RE = re.compile(
    r"^\[(?P<group>\d+):(?P<replicate>\d+)\]_(?P<type>.) "
    r"XicScan (?P<sample>.+) "
    r"(?P<mz_lower>[\d.]+)-(?P<mz_upper>[\d.]+)$"
)


def open_qsession(path: str | Path, password: str = QSESSION_PASSWORD) -> "QSessionReader":
    """Open a ``.qsession`` file and return a :class:`QSessionReader`."""
    return QSessionReader(path, password=password)


class QSessionReader:
    """Pure-Python reader over a single ``.qsession`` quantitation session."""

    def __init__(self, path: str | Path, password: str = QSESSION_PASSWORD):
        self.path = Path(path).resolve()
        if self.path.suffix.lower() != ".qsession":
            raise ValueError(f"Only .qsession is supported, got: {self.path.name}")
        if not self.path.exists():
            raise FileNotFoundError(f"QSession file does not exist: {self.path}")

        raw_db = decrypt_database(
            self.path, password, page_size=QSESSION_PAGE_SIZE
        )
        self._conn = self._connect(raw_db)
        self._conn.row_factory = sqlite3.Row

        # caches
        self._xic_info_cache: dict[str, XicInfo] = {}
        self._xic_data_cache: dict[str, XicChromatogram] = {}
        self._chromatogram_cache: dict[tuple[int, int], XicChromatogram | None] = {}
        self._multidata_cache: dict[str, Any] | None = None

    # ------------------------------------------------------------------ #
    # lifecycle
    # ------------------------------------------------------------------ #
    def __enter__(self) -> "QSessionReader":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def close(self) -> None:
        if self._conn is not None:
            self._conn.close()
            self._conn = None  # type: ignore[assignment]

    @staticmethod
    def _connect(raw_db: bytes) -> sqlite3.Connection:
        conn = sqlite3.connect(":memory:")
        if hasattr(conn, "deserialize"):
            conn.deserialize(raw_db)
            return conn
        # Fallback for Python builds without serialize/deserialize support.
        conn.close()
        import tempfile

        tmp = tempfile.NamedTemporaryFile(suffix=".db", delete=False)
        try:
            tmp.write(raw_db)
            tmp.flush()
        finally:
            tmp.close()
        return sqlite3.connect(tmp.name)

    # ------------------------------------------------------------------ #
    # metadata
    # ------------------------------------------------------------------ #
    @property
    def version(self) -> str | None:
        """Software version that created the session (e.g. ``'MultiQuant MD'``)."""
        row = self._conn.execute(
            "SELECT version FROM VersionInformation LIMIT 1"
        ).fetchone()
        return str(row[0]) if row and row[0] else None

    @property
    def qmap_version(self) -> str | None:
        """Quantitation method version (e.g. ``'MultiQuant 2.0'``)."""
        row = self._conn.execute(
            "SELECT version FROM QMapVersionInformation LIMIT 1"
        ).fetchone()
        return str(row[0]) if row and row[0] else None

    @property
    def locked(self) -> bool:
        """Whether the session is locked (read-only)."""
        row = self._conn.execute(
            "SELECT Locked FROM LockInformation LIMIT 1"
        ).fetchone()
        if row is None:
            return False
        val = row[0]
        if isinstance(val, str):
            return val.lower() in ("true", "1", "yes")
        return bool(val)

    @property
    def sample_keys(self) -> list[str]:
        """All distinct sample keys referenced by XICs."""
        rows = self._conn.execute(
            "SELECT DISTINCT SampleKey FROM XicRawTable ORDER BY SampleKey"
        ).fetchall()
        return [str(r[0]) for r in rows if r[0]]

    @property
    def xic_count(self) -> int:
        """Total number of XIC rows in the session."""
        return self._conn.execute("SELECT COUNT(*) FROM XicRawTable").fetchone()[0]

    # ------------------------------------------------------------------ #
    # XIC queries
    # ------------------------------------------------------------------ #
    def list_xics(self, sample_key: str | None = None) -> list[XicInfo]:
        """Return metadata for all (or filtered) XICs.

        Parameters
        ----------
        sample_key :
            If given, only XICs belonging to this sample key are returned.
        """
        if sample_key is not None:
            rows = self._conn.execute(
                "SELECT ID, SampleKey, status FROM XicRawTable "
                "WHERE SampleKey = ? ORDER BY ID",
                (sample_key,),
            ).fetchall()
        else:
            rows = self._conn.execute(
                "SELECT ID, SampleKey, status FROM XicRawTable ORDER BY ID"
            ).fetchall()

        result: list[XicInfo] = []
        for row in rows:
            info = self._parse_xic_id(row["ID"], row["SampleKey"], row["status"])
            if info is not None:
                result.append(info)
                self._xic_info_cache[info.xic_id] = info
        return result

    def get_xic(self, xic_id: str) -> XicChromatogram:
        """Retrieve a single XIC by its ID."""
        cached = self._xic_data_cache.get(xic_id)
        if cached is not None:
            return cached

        row = self._conn.execute(
            "SELECT ID, SampleKey, Xdata, Ydata, status FROM XicRawTable WHERE ID = ?",
            (xic_id,),
        ).fetchone()
        if row is None:
            raise KeyError(f"XIC not found: {xic_id!r}")

        return self._build_xic(row)

    def get_xic_by_mz_range(
        self,
        sample_key: str,
        mz_lower: float,
        mz_upper: float,
    ) -> XicChromatogram | None:
        """Retrieve the first XIC matching a sample key and m/z window."""
        # Use a tolerance-based match to account for floating-point drift
        tol = 0.001
        rows = self._conn.execute(
            "SELECT ID, SampleKey, Xdata, Ydata, status FROM XicRawTable "
            "WHERE SampleKey = ? AND ID LIKE ?",
            (sample_key, f"% {mz_lower:.6f}-{mz_upper:.6f}"),
        ).fetchall()

        # Fallback: fuzzy match
        if not rows:
            all_rows = self._conn.execute(
                "SELECT ID, SampleKey, Xdata, Ydata, status FROM XicRawTable "
                "WHERE SampleKey = ?",
                (sample_key,),
            ).fetchall()
            for r in all_rows:
                info = self._parse_xic_id(r["ID"], r["SampleKey"], r["status"])
                if info and abs(info.mz_lower - mz_lower) < tol and abs(info.mz_upper - mz_upper) < tol:
                    rows = [r]
                    break

        if not rows:
            return None
        return self._build_xic(rows[0])

    def iter_xics(
        self,
        sample_key: str | None = None,
        return_arrays: bool = False,
    ) -> Iterator[XicChromatogram]:
        """Iterate over all XIC chromatograms.

        Parameters
        ----------
        sample_key :
            If given, only XICs for this sample key are yielded.
        return_arrays :
            If ``True``, ``times`` and ``intensities`` are numpy arrays;
            otherwise Python lists.
        """
        if sample_key is not None:
            rows = self._conn.execute(
                "SELECT ID, SampleKey, Xdata, Ydata, status FROM XicRawTable "
                "WHERE SampleKey = ? ORDER BY ID",
                (sample_key,),
            )
        else:
            rows = self._conn.execute(
                "SELECT ID, SampleKey, Xdata, Ydata, status FROM XicRawTable ORDER BY ID"
            )

        for row in rows:
            xic = self._build_xic(row, return_arrays=return_arrays)
            yield xic

    # ------------------------------------------------------------------ #
    # audit trail
    # ------------------------------------------------------------------ #
    def list_audit_records(self, limit: int = 100) -> list[dict]:
        """Return recent audit-trail records."""
        rows = self._conn.execute(
            "SELECT EventName, EventTypeID, ATRecordTimeStamp, UserName, "
            "FullUserName, Reason, ESig, ATSaveRecordTimeStamp "
            "FROM AuditTrailRecords ORDER BY ATRecordTimeStamp DESC LIMIT ?",
            (limit,),
        ).fetchall()
        return [
            {
                "event_name": r["EventName"],
                "event_type_id": r["EventTypeID"],
                "timestamp": r["ATRecordTimeStamp"],
                "username": r["UserName"],
                "full_username": r["FullUserName"],
                "reason": r["Reason"],
                "esig": bool(r["ESig"]),
                "save_timestamp": r["ATSaveRecordTimeStamp"],
            }
            for r in rows
        ]

    # ------------------------------------------------------------------ #
    # internals
    # ------------------------------------------------------------------ #
    @staticmethod
    def _parse_xic_id(xic_id: str, sample_key: str, status: int) -> XicInfo | None:
        m = _XIC_ID_RE.match(xic_id)
        if not m:
            return None
        return XicInfo(
            xic_id=xic_id,
            sample_key=sample_key,
            mz_lower=float(m.group("mz_lower")),
            mz_upper=float(m.group("mz_upper")),
            group_index=int(m.group("group")),
            replicate_index=int(m.group("replicate")),
            status=int(status),
        )

    @staticmethod
    def _sample_info_from_dict(s: dict[str, Any]) -> QuantSampleInfo:
        """Build a QuantSampleInfo from a MultiData sample dict."""
        return QuantSampleInfo(
            index=s["index"],
            sample_name=s["sample_name"] or "",
            sample_id=s["sample_id"] or "",
            sample_type=s.get("sample_type", 0),
            sample_comment=s.get("sample_comment"),
            dilution_factor=s.get("dilution_factor", 0.0),
            injection_volume=s.get("injection_volume", 0.0),
            user_name=s.get("user_name"),
            acq_method_name=s.get("acq_method_name"),
            instrument_name=s.get("instrument_name"),
            instrument_serial_number=s.get("instrument_serial_number"),
            batch_name=s.get("batch_name"),
            barcode=s.get("barcode"),
            scanned_barcode=s.get("scanned_barcode"),
            autosampler_method_supports_barcode=s.get("autosampler_method_supports_barcode", False),
            sample_comparison=s.get("sample_comparison", False),
            ms_method=s.get("ms_method"),
            lc_method=s.get("lc_method"),
            sample_signature=s.get("sample_signature"),
            rack=s.get("rack"),
            plate=s.get("plate"),
            vial=s.get("vial"),
            acquisition_date=s.get("acquisition_date"),
        )

    @staticmethod
    def _peak_info_from_dict(p: dict[str, Any], sample_index: int) -> QuantPeakInfo:
        """Build a QuantPeakInfo from a MultiData peak dict."""
        return QuantPeakInfo(
            sample_index=sample_index,
            compound_index=p["compound_index"],
            peak_index=p.get("peak_index", -1),
            use_for_calibration=p.get("use_for_calibration", False),
            peak_comment=p.get("peak_comment"),
            actual_concentration=p.get("actual_concentration", 0.0),
            failed_query=p.get("failed_query", False),
            valid_integration=p.get("valid_integration", False),
            modified=p.get("modified", False),
            retention_time=p.get("retention_time", 0.0),
            area=p.get("area", 0.0),
            corrected_area=p.get("corrected_area", 0.0),
            height=p.get("height", 0.0),
            corrected_height=p.get("corrected_height", 0.0),
            start_rt=p.get("start_rt", 0.0),
            start_y=p.get("start_y", 0.0),
            end_rt=p.get("end_rt", 0.0),
            end_y=p.get("end_y", 0.0),
            half_height_start_rt=p.get("half_height_start_rt", 0.0),
            half_height_end_rt=p.get("half_height_end_rt", 0.0),
            noise=p.get("noise", -1.0),
            signal_to_noise=p.get("signal_to_noise", -1.0),
            profile_type=p.get("profile_type", 0),
            peak_type=p.get("peak_type", 0),
            apex_rt=p.get("apex_rt", 0.0),
            apex_y=p.get("apex_y", 0.0),
            region_area=p.get("region_area", 0.0),
            region_height=p.get("region_height", 0.0),
            s_mrm_retention_time_shift=p.get("s_mrm_retention_time_shift", False),
            row_hidden=p.get("row_hidden", False),
            reportable=p.get("reportable", False),
            molecular_weight=p.get("molecular_weight", 0.0),
            original_area=p.get("original_area", 0.0),
            override_experiment_index=p.get("override_experiment_index", 0),
            points_across_baseline=p.get("points_across_baseline", 0),
            points_across_half_height=p.get("points_across_half_height", 0),
            integration_parameters=p.get("integration_parameters"),
            profile=p.get("profile"),
            custom_fields=p.get("custom_fields"),
            custom_peak_fields=p.get("custom_peak_fields"),
            start_x5_pct_height=p.get("start_x5_pct_height", 0.0),
            end_x5_pct_height=p.get("end_x5_pct_height", 0.0),
            start_x10_pct_height=p.get("start_x10_pct_height", 0.0),
            end_x10_pct_height=p.get("end_x10_pct_height", 0.0),
            std_addn_actual_concentration=p.get("std_addn_actual_concentration", 0.0),
            extracted_ms_ms=p.get("extracted_ms_ms"),
            super_group_id=p.get("super_group_id"),
        )

    def _build_xic(
        self,
        row: sqlite3.Row,
        return_arrays: bool = False,
    ) -> XicChromatogram:
        xic_id = row["ID"]
        cached = self._xic_data_cache.get(xic_id)
        if cached is not None:
            if not return_arrays:
                return cached
            # caller wants arrays — convert from the cached list version
            times = np.asarray(cached.times, dtype=np.float64)
            intensities = np.asarray(cached.intensities, dtype=np.float64)
            return XicChromatogram(
                xic_id=xic_id,
                sample_key=cached.sample_key,
                times=times,
                intensities=intensities,
                mz_lower=cached.mz_lower,
                mz_upper=cached.mz_upper,
                status=cached.status,
            )

        xdata: bytes = row["Xdata"]
        ydata: bytes = row["Ydata"]

        times = self._unpack_doubles(xdata)
        intensities = self._unpack_doubles(ydata)

        info = self._parse_xic_id(xic_id, row["SampleKey"], row["status"])
        if info is None:
            info = XicInfo(
                xic_id=xic_id,
                sample_key=row["SampleKey"],
                mz_lower=0.0,
                mz_upper=0.0,
                group_index=0,
                replicate_index=0,
                status=row["status"],
            )

        if not return_arrays:
            times = times.tolist() if isinstance(times, np.ndarray) else list(times)
            intensities = (
                intensities.tolist()
                if isinstance(intensities, np.ndarray)
                else list(intensities)
            )

        xic = XicChromatogram(
            xic_id=xic_id,
            sample_key=info.sample_key,
            times=times,
            intensities=intensities,
            mz_lower=info.mz_lower,
            mz_upper=info.mz_upper,
            status=info.status,
        )
        self._xic_data_cache[xic_id] = xic
        return xic

    # ------------------------------------------------------------------ #
    # compounds (RTParts)
    # ------------------------------------------------------------------ #
    def list_compounds(self) -> list[CompoundInfo]:
        """Return all compound definitions from the RTParts table."""
        return list(self.iter_compounds())

    def iter_compounds(self) -> Iterator[CompoundInfo]:
        """Iterate over compound definitions from the RTParts table."""
        stream = load_rtparts_stream(self._conn)
        compounds, _meta = read_compounds(stream)
        for c in compounds:
            ev1 = c.get("extraction_values1") or []
            ev2 = c.get("extraction_values2") or []
            yield CompoundInfo(
                name=c["name"] or "",
                group_name=c["group_name"] or "",
                formula=c.get("formula"),
                charge_formula=c.get("charge_formula"),
                adduct_formula=c.get("adduct_formula"),
                precursor_mass=c.get("precursor_mass") or 0.0,
                fragment_mass=c.get("fragment_mass") or 0.0,
                extraction_type=c.get("extraction_type") or 0,
                period=c.get("period") or 0,
                experiment=c.get("experiment") or 0,
                mz_lower=ev1[0] if ev1 else 0.0,
                mz_upper=ev2[0] if ev2 else 0.0,
                is_analyte=bool(c.get("is_analyte")),
                is_reportable=bool(c.get("is_reportable")),
                is_non_targeted=bool(c.get("is_non_targeted")),
                is_summed=bool(c.get("is_summed")),
                is_from_multi_period_data=bool(c.get("is_from_multi_period_data")),
                isotope_index=c.get("isotope_index") or 0,
                expected_mw=c.get("expected_mw") or 0.0,
                units=c.get("units") or "",
                comment=c.get("comment"),
                internal_std_name=c.get("internal_std_name"),
                regression_area=c.get("regression_area"),
                regression_type=c.get("regression_type"),
                regression_weighting=c.get("regression_weighting"),
                use_auto_regression=c.get("use_auto_regression"),
                integration_parameters=c.get("integration_parameters"),
                acquisition_indices=c.get("acquisition_indices"),
                summed_compounds=c.get("summed_compounds"),
                extraction_values1=c.get("extraction_values1"),
                extraction_values2=c.get("extraction_values2"),
            )

    # ------------------------------------------------------------------ #
    # full RTParts (samples + peaks)
    # ------------------------------------------------------------------ #
    def _load_multidata(self) -> dict[str, Any]:
        """Lazy-load the full MultiData object graph from RTParts."""
        if self._multidata_cache is None:
            stream = load_rtparts_stream(self._conn)
            self._multidata_cache = read_multidata(stream)
        return self._multidata_cache

    def list_samples(self) -> list[QuantSampleInfo]:
        """Return all sample definitions from the RTParts table."""
        return list(self.iter_samples())

    def iter_samples(self) -> Iterator[QuantSampleInfo]:
        """Iterate over sample definitions from the RTParts table."""
        data = self._load_multidata()
        for s in data["samples"]:
            yield self._sample_info_from_dict(s)

    def iter_peaks(self, sample_index: int | None = None) -> Iterator[QuantPeakInfo]:
        """Iterate over peak results.

        Parameters
        ----------
        sample_index :
            If given, only peaks for this sample are yielded.
            Otherwise all peaks across all samples are yielded.
        """
        data = self._load_multidata()
        samples = data["samples"]
        if sample_index is not None:
            samples = [samples[sample_index]]
        for s in samples:
            si = s["index"]
            for p in s["peaks"]:
                yield self._peak_info_from_dict(p, si)

    def get_sample(self, index: int) -> QuantSampleInfo:
        """Return a single sample by index."""
        data = self._load_multidata()
        if index < 0 or index >= len(data["samples"]):
            raise IndexError(f"Sample index {index} out of range ({len(data['samples'])} samples)")
        return self._sample_info_from_dict(data["samples"][index])

    def get_peak(self, sample_index: int, compound_index: int) -> QuantPeakInfo:
        """Return the peak for a specific sample and compound."""
        data = self._load_multidata()
        if sample_index < 0 or sample_index >= len(data["samples"]):
            raise IndexError(f"Sample index {sample_index} out of range")
        sample = data["samples"][sample_index]
        if compound_index < 0 or compound_index >= len(sample["peaks"]):
            raise IndexError(f"Compound index {compound_index} out of range")
        return self._peak_info_from_dict(sample["peaks"][compound_index], sample_index)

    def get_chromatogram(
        self,
        sample_index: int,
        compound_index: int,
    ) -> XicChromatogram | None:
        """Return the full extracted-ion chromatogram for a compound × sample.

        The chromatogram is looked up in ``XicRawTable`` using the XIC ID
        constructed from the compound's period / experiment / extraction
        window and the sample's signature. Results are cached per
        *(sample_index, compound_index)* pair.

        Parameters
        ----------
        sample_index :
            Index of the sample (row in the results table).
        compound_index :
            Index of the compound (column in the results table).

        Returns
        -------
        XicChromatogram | None
            The chromatogram if found in the cache, otherwise ``None``.
        """
        key = (sample_index, compound_index)
        if key in self._chromatogram_cache:
            return self._chromatogram_cache[key]

        data = self._load_multidata()
        if sample_index < 0 or sample_index >= len(data["samples"]):
            raise IndexError(f"Sample index {sample_index} out of range")
        sample = data["samples"][sample_index]
        if compound_index < 0 or compound_index >= len(data["compounds"]):
            raise IndexError(f"Compound index {compound_index} out of range")
        compound = data["compounds"][compound_index]

        ev1 = compound.get("extraction_values1") or []
        ev2 = compound.get("extraction_values2") or []
        if not ev1 or not ev2:
            self._chromatogram_cache[key] = None
            return None

        period = compound.get("period", 0)
        experiment = compound.get("experiment", 0)
        sample_signature = sample.get("sample_signature") or ""
        xic_id = (
            f"[{period}:{experiment}]_X "
            f"XicScan {sample_signature} "
            f"{ev1[0]}-{ev2[0]}"
        )

        row = self._conn.execute(
            "SELECT Xdata, Ydata, status FROM XicRawTable WHERE ID = ?",
            (xic_id,),
        ).fetchone()
        if row is None:
            self._chromatogram_cache[key] = None
            return None

        xblob, yblob, status = row
        if not xblob or not yblob:
            self._chromatogram_cache[key] = None
            return None

        times = struct.unpack(f"<{len(xblob) // 8}d", xblob)
        intensities = struct.unpack(f"<{len(yblob) // 8}d", yblob)
        chrom = XicChromatogram(
            xic_id=xic_id,
            sample_key=sample_signature,
            times=times,
            intensities=intensities,
            mz_lower=ev1[0],
            mz_upper=ev2[0],
            status=status,
        )
        self._chromatogram_cache[key] = chrom
        return chrom

    def results_matrix(self) -> list[list[QuantPeakInfo]]:
        """Return the full results table as a 2-D matrix.

        ``results[sample_index][compound_index]`` gives the
        :class:`QuantPeakInfo` for that cell.
        """
        data = self._load_multidata()
        return [
            [
                self._peak_info_from_dict(p, s["index"])
                for p in s["peaks"]
            ]
            for s in data["samples"]
        ]

    @staticmethod
    def _unpack_doubles(blob: bytes) -> np.ndarray | Sequence[float]:
        if len(blob) % 8 != 0:
            return []
        # read-only view is safe — the sqlite blob buffer is immutable
        return np.frombuffer(blob, dtype=np.float64)
