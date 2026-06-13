"""Pure-Python reader for SCIEX qsession (quantitation session) files.

A ``.qsession`` file is an encrypted SQLite database produced by SCIEX
quantitation software (e.g. MultiQuant). It stores extracted ion
chromatograms (XICs), audit trails, and result metadata.

The encryption scheme is the same AES-128-OFB SEE cipher used for
``.wiff2`` files, but with a different password and page size:

* Page size: **1024** bytes (wiff2 uses 4096)
* Password: ``"PQS1 is not Sirius"``

Usage::

    from pyx500r import open_qsession

    with open_qsession("analysis.qsession") as qs:
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
    LibraryHit,
    QuantPeakInfo,
    QuantSampleInfo,
    XicChromatogram,
    XicInfo,
)
from .rtparts import load_rtparts_stream, read_compounds, read_multidata
from .xic_gap import build_xic_index, parse_xic_blobs

_XIC_ID_RE = re.compile(
    r"^\[(?P<group>\d+):(?P<replicate>\d+)\]_(?P<type>.) "
    r"XicScan (?P<sample>.+) "
    r"(?P<mz_lower>[\d.]+)-(?P<mz_upper>[\d.]+)$"
)


def open_qsession(path: str | Path, password: str = QSESSION_PASSWORD) -> "QSessionReader":
    """Open a ``.qsession`` file and return a :class:`QSessionReader`."""
    return QSessionReader(path, password=password)


class QSessionReader:
    """Pure-Python reader over a single ``.qsession`` quantitation session.

    Supports both v1 (1024-byte pages) and v2 (4096-byte pages) formats.
    Page size is auto-detected on open.
    """

    # Both known page sizes, tried in order
    _PAGE_SIZES = (4096, 1024)

    def __init__(self, path: str | Path, password: str = QSESSION_PASSWORD):
        self.path = Path(path).resolve()
        if self.path.suffix.lower() != ".qsession":
            raise ValueError(f"Only .qsession is supported, got: {self.path.name}")
        if not self.path.exists():
            raise FileNotFoundError(f"QSession file does not exist: {self.path}")

        raw_db = None
        last_err = None
        for page_size in self._PAGE_SIZES:
            try:
                raw_db = decrypt_database(
                    self.path, password, page_size=page_size
                )
                self._conn = self._connect(raw_db)
                break
            except Exception as e:
                last_err = e
                continue

        if self._conn is None:
            raise OSError(
                f"Cannot open qsession (tried page sizes {self._PAGE_SIZES}): {last_err}"
            ) from last_err

        self._conn.row_factory = sqlite3.Row

        # Detect format version
        cur = self._conn.cursor()
        cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
        table_names = {r[0] for r in cur.fetchall()}
        self._is_v2 = 'MultiData' in table_names

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
        """Connect to an in-memory SQLite database from raw bytes.

        Tries ``deserialize()`` first (fast, in-memory), falls back to a
        temp file if the SQLite version doesn't support the page size.
        """
        # Try deserialize first
        if hasattr(sqlite3.Connection, "deserialize"):
            try:
                conn = sqlite3.connect(":memory:")
                conn.deserialize(raw_db)
                return conn
            except sqlite3.Error:
                pass  # fall through to temp file

        # Temp file fallback
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

    def _peak_info_from_dict(self, p: dict[str, Any], sample_index: int, xic_result: dict[str, Any] | None = None) -> QuantPeakInfo:
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
            xic_result=xic_result,
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
    # compounds
    # ------------------------------------------------------------------ #
    def list_compounds(self) -> list[CompoundInfo]:
        """Return all compound definitions."""
        return list(self.iter_compounds())

    def iter_compounds(self) -> Iterator[CompoundInfo]:
        """Iterate over compound definitions.

        V1: from RTParts BinaryFormatter blobs.
        V2: from MultiData.customColumns blob (BF-serialized QuantCompounds).
        """
        if self._is_v2:
            yield from self._iter_compounds_v2()
        else:
            yield from self._iter_compounds_v1()

    def _iter_compounds_v1(self) -> Iterator[CompoundInfo]:
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

    def _iter_compounds_v2(self) -> Iterator[CompoundInfo]:
        """V2: extract compound names from QualPeak (latest QuantResultId per PeakIndex)."""
        cur = self._conn.cursor()
        cur.execute(
            "SELECT q.PeakIndex, q.CompoundName, q.Formula"
            " FROM QualPeak q"
            " INNER JOIN ("
            "  SELECT PeakIndex, MAX(QuantResultId) AS MaxResult"
            "  FROM QualPeak GROUP BY PeakIndex"
            " ) latest ON q.PeakIndex = latest.PeakIndex"
            "  AND q.QuantResultId = latest.MaxResult"
            " ORDER BY q.PeakIndex"
        )
        for r in cur.fetchall():
            ci = r["PeakIndex"]
            yield CompoundInfo(
                name=r["CompoundName"] or f"Compound_{ci}",
                group_name="",
                formula=r["Formula"],
                charge_formula=None, adduct_formula=None,
                precursor_mass=0.0, fragment_mass=0.0,
                extraction_type=0, period=0, experiment=0,
                mz_lower=0.0, mz_upper=0.0,
                is_analyte=True, is_reportable=True,
                is_non_targeted=False, is_summed=False,
                is_from_multi_period_data=False,
                isotope_index=0, expected_mw=0.0, units="",
                comment=None, internal_std_name=None,
                regression_area=None, regression_type=None,
                regression_weighting=None, use_auto_regression=None,
            )

    # ------------------------------------------------------------------ #
    # full RTParts (samples + peaks) — V1 only
    # ------------------------------------------------------------------ #
    def _load_multidata(self) -> dict[str, Any]:
        """Lazy-load the full MultiData object graph from RTParts (V1 only)."""
        if self._multidata_cache is None:
            stream = load_rtparts_stream(self._conn)
            data = read_multidata(stream)
            raw = stream.getvalue()
            blobs = parse_xic_blobs(raw)
            data["xic_lookup"] = build_xic_index(blobs, len(data["samples"]))
            data["xic_blobs"] = blobs
            self._multidata_cache = data
        return self._multidata_cache

    def list_samples(self) -> list[QuantSampleInfo]:
        """Return all sample definitions."""
        return list(self.iter_samples())

    def iter_samples(self) -> Iterator[QuantSampleInfo]:
        """Iterate over sample definitions."""
        if self._is_v2:
            yield from self._iter_samples_v2()
        else:
            yield from self._iter_samples_v1()

    def _iter_samples_v1(self) -> Iterator[QuantSampleInfo]:
        data = self._load_multidata()
        for s in data["samples"]:
            yield self._sample_info_from_dict(s)

    def _iter_samples_v2(self) -> Iterator[QuantSampleInfo]:
        cur = self._conn.cursor()
        cur.execute("SELECT * FROM MultiSample ORDER BY Id")
        for row in cur.fetchall():
            yield QuantSampleInfo(
                index=row["Id"] - 1,  # 0-based
                sample_name=row["sampleName"] or "",
                sample_id=row["sampleId"] or "",
                sample_type=row["sampleType"] or 0,
                sample_comment=row["sampleComment"],
                dilution_factor=row["dilutionFactor"] or 1.0,
                injection_volume=row["injectionVolume"] or 0.0,
                user_name=row["userName"],
                acq_method_name=row["acqMethodName"],
                instrument_name=row["instrumentName"],
                instrument_serial_number=row["instrumentSerialNumber"],
                batch_name=row["batchName"],
                barcode=row["barcode"],
                scanned_barcode=row["scannedBarcode"],
                autosampler_method_supports_barcode=bool(row["autosamplerMethodSupportsBarcode"]),
                sample_comparison=bool(row["sampleComparison"]),
                ms_method=row["msMethod"],
                lc_method=row["lcMethod"],
                sample_signature=row["sampleSignature"],
                rack=row["rack"],
                plate=row["plate"],
                vial=row["vial"],
                acquisition_date=row["dateTime"],
            )

    # ------------------------------------------------------------------ #
    # peaks
    # ------------------------------------------------------------------ #
    def iter_peaks(self, sample_index: int | None = None) -> Iterator[QuantPeakInfo]:
        """Iterate over peak results.

        Parameters
        ----------
        sample_index :
            If given, only peaks for this sample are yielded.
            Otherwise all peaks across all samples are yielded.
        """
        if self._is_v2:
            yield from self._iter_peaks_v2(sample_index)
        else:
            yield from self._iter_peaks_v1(sample_index)

    def _iter_peaks_v1(self, sample_index: int | None = None) -> Iterator[QuantPeakInfo]:
        data = self._load_multidata()
        samples = data["samples"]
        if sample_index is not None:
            samples = [samples[sample_index]]
        xic = data.get("xic_lookup")
        for s in samples:
            si = s["index"]
            for p in s["peaks"]:
                xic_result = xic.get((si, p["compound_index"])) if xic else None
                yield self._peak_info_from_dict(p, si, xic_result)

    def _iter_peaks_v2(self, sample_index: int | None = None) -> Iterator[QuantPeakInfo]:
        cur = self._conn.cursor()
        if sample_index is not None:
            cur.execute(
                "SELECT * FROM MultiPeak WHERE SampleId = ? ORDER BY Id",
                (sample_index + 1,),
            )
        else:
            cur.execute("SELECT * FROM MultiPeak ORDER BY Id")
        for row in cur.fetchall():
            yield self._peak_info_from_v2_row(row)

    def get_sample(self, index: int) -> QuantSampleInfo:
        """Return a single sample by index."""
        if self._is_v2:
            cur = self._conn.cursor()
            cur.execute("SELECT * FROM MultiSample WHERE Id = ?", (index + 1,))
            row = cur.fetchone()
            if row is None:
                raise IndexError(f"Sample index {index} out of range")
            return next(self._iter_samples_v2_for_row(row))
        data = self._load_multidata()
        if index < 0 or index >= len(data["samples"]):
            raise IndexError(f"Sample index {index} out of range ({len(data['samples'])} samples)")
        return self._sample_info_from_dict(data["samples"][index])

    def _iter_samples_v2_for_row(self, row) -> Iterator[QuantSampleInfo]:
        yield QuantSampleInfo(
            index=row["Id"] - 1,
            sample_name=row["sampleName"] or "",
            sample_id=row["sampleId"] or "",
            sample_type=row["sampleType"] or 0,
            sample_comment=row["sampleComment"],
            dilution_factor=row["dilutionFactor"] or 1.0,
            injection_volume=row["injectionVolume"] or 0.0,
            user_name=row["userName"],
            acq_method_name=row["acqMethodName"],
            instrument_name=row["instrumentName"],
            instrument_serial_number=row["instrumentSerialNumber"],
            batch_name=row["batchName"],
            barcode=row["barcode"],
            scanned_barcode=row["scannedBarcode"],
            autosampler_method_supports_barcode=bool(row["autosamplerMethodSupportsBarcode"]),
            sample_comparison=bool(row["sampleComparison"]),
            ms_method=row["msMethod"],
            lc_method=row["lcMethod"],
            sample_signature=row["sampleSignature"],
            rack=row["rack"],
            plate=row["plate"],
            vial=row["vial"],
            acquisition_date=row["dateTime"],
        )

    def get_peak(self, sample_index: int, compound_index: int) -> QuantPeakInfo:
        """Return the peak for a specific sample and compound."""
        if self._is_v2:
            cur = self._conn.cursor()
            # In v2, PeakIndex within a sample corresponds to compound_index
            # MultiPeak.Id = sample_id * N + compound_index + 1 (roughly)
            # Better: use SampleId and PeakIndex
            cur.execute(
                "SELECT * FROM MultiPeak WHERE SampleId = ? AND PeakIndex = ?",
                (sample_index + 1, compound_index),
            )
            row = cur.fetchone()
            if row is None:
                raise IndexError(
                    f"No peak for sample {sample_index}, compound {compound_index}"
                )
            return self._peak_info_from_v2_row(row)

        data = self._load_multidata()
        if sample_index < 0 or sample_index >= len(data["samples"]):
            raise IndexError(f"Sample index {sample_index} out of range")
        sample = data["samples"][sample_index]
        if compound_index < 0 or compound_index >= len(sample["peaks"]):
            raise IndexError(f"Compound index {compound_index} out of range")
        xic = data.get("xic_lookup")
        xic_result = xic.get((sample_index, compound_index)) if xic else None
        return self._peak_info_from_dict(sample["peaks"][compound_index], sample_index, xic_result)

    def _peak_info_from_v2_row(self, row: sqlite3.Row) -> QuantPeakInfo:
        """Convert a v2 MultiPeak row to QuantPeakInfo,
        enriched with XIC data from QualPeak."""
        ci = row["PeakIndex"]
        si = row["SampleId"] - 1
        pi = row["Id"]
        xic_result = self._v2_xic_result(si, ci)

        return QuantPeakInfo(
            sample_index=si,
            compound_index=ci,
            peak_index=row["Id"],
            use_for_calibration=bool(row["Use"]),
            peak_comment=row["PeakComment"],
            actual_concentration=row["ActualConcentration"] or 0.0,
            failed_query=bool(row["FailedQuery"]),
            valid_integration=bool(row["ValidIntegration"]),
            modified=bool(row["Modified"]),
            retention_time=row["RetentionTime"] or 0.0,
            area=row["Area"] or 0.0,
            corrected_area=row["CorrectedArea"] or 0.0,
            height=row["Height"] or 0.0,
            corrected_height=row["CorrectedHeight"] or 0.0,
            start_rt=row["StartRt"] or 0.0,
            start_y=row["StartY"] or 0.0,
            end_rt=row["EndRt"] or 0.0,
            end_y=row["EndY"] or 0.0,
            half_height_start_rt=row["HalfHeightStartRt"] or 0.0,
            half_height_end_rt=row["HalfHeightEndRt"] or 0.0,
            noise=row["Noise"] or 0.0,
            signal_to_noise=xic_result.get("_signalToNoise") or 0.0,
            profile_type=row["ProfileType"] or 0,
            peak_type=row["PeakType"] or 0,
            apex_rt=row["ApexRt"] or 0.0,
            apex_y=row["ApexY"] or 0.0,
            region_area=row["RegionArea"] or 0.0,
            region_height=row["RegionHeight"] or 0.0,
            s_mrm_retention_time_shift=bool(row["SMrmRetentionTimeShift"]),
            row_hidden=bool(row["RowHidden"]),
            reportable=bool(row["Reportable"]),
            molecular_weight=row["MolecularWeight"] or 0.0,
            original_area=row["OriginalArea"] or 0.0,
            override_experiment_index=row["OverrideExperimentIndex"] or 0,
            points_across_baseline=row["PointsAcrossBaseline"] or 0,
            points_across_half_height=row["PointsAcrossHalfHeight"] or 0,
            integration_parameters=None,
            profile=None,
            custom_fields=None,
            custom_peak_fields=None,
            start_x5_pct_height=row["StartX5PctHeight"] or 0.0,
            end_x5_pct_height=row["EndX5PctHeight"] or 0.0,
            start_x10_pct_height=row["StartX10PctHeight"] or 0.0,
            end_x10_pct_height=row["EndX10PctHeight"] or 0.0,
            std_addn_actual_concentration=row["StdAddnActualConcentration"] or 0.0,
            extracted_ms_ms=row["ExtractedMsms"],
            xic_result=xic_result,
            super_group_id=row["SuperGroupId"],
        )

    def _v2_xic_result(self, sample_index: int, compound_index: int) -> dict[str, Any]:
        """Build XIC result dict from QualPeak (latest QuantResultId per peak)."""
        if not hasattr(self, '_v2_qualpeak_index'):
            cur = self._conn.cursor()
            cur.execute(
                "SELECT q.QuantSampleId, q.PeakIndex,"
                " q.ContainsMSMS, q.HasBeenCalculated, q.HasLibraryBeenSearched,"
                " q.FoundAtMass, q.FoundAtRt, q.FoundAtRtApex,"
                " q.FoundAtRtStart, q.FoundAtRtEnd,"
                " q.ExtractionMass, q.ExtractionWidth,"
                " q.IsotopeRatioDiffFromExpected,"
                " q.LibrarySearchResult,"
                " q.Formula, q.BaseMass, q.SignalToNoise,"
                " q.Area, q.Intensity, q.Charge,"
                " q.ExpectedRt, q.ExpectedRtWidth,"
                " q.IsQualifier, q.IsInternalStandard,"
                " q.MsMsExperimentIndex, q.MsMsRetentionTime,"
                " q.Ms1Cycle, q.MsMsCycle, q.ModificationText"
                " FROM QualPeak q"
                " INNER JOIN ("
                "  SELECT QuantSampleId, PeakIndex, MAX(QuantResultId) AS MaxResult"
                "  FROM QualPeak GROUP BY QuantSampleId, PeakIndex"
                " ) latest ON q.QuantSampleId = latest.QuantSampleId"
                "  AND q.PeakIndex = latest.PeakIndex"
                "  AND q.QuantResultId = latest.MaxResult"
            )
            self._v2_qualpeak_index = {}
            for r in cur.fetchall():
                key = (r["QuantSampleId"] - 1, r["PeakIndex"])
                entry: dict[str, Any] = {
                    "_containsMSMS": bool(r["ContainsMSMS"]),
                    "_hasBeenCalculated": bool(r["HasBeenCalculated"]),
                    "<HasLibraryBeenSearched>k__BackingField": bool(r["HasLibraryBeenSearched"]),
                    "_foundAtMass": r["FoundAtMass"],
                    "_foundAtRt": r["FoundAtRt"],
                    "_foundAtRtApex": r["FoundAtRtApex"],
                    "<FoundAtRtStart>k__BackingField": r["FoundAtRtStart"],
                    "<FoundAtRtEnd>k__BackingField": r["FoundAtRtEnd"],
                    "_extractionMass": r["ExtractionMass"],
                    "_extractionWidth": r["ExtractionWidth"],
                    "<IsoptopeRatioDiffFromExpected>k__BackingField": r["IsotopeRatioDiffFromExpected"],
                    "_formula": r["Formula"],
                    "_baseMass": r["BaseMass"],
                    "_signalToNoise": r["SignalToNoise"],
                    "_area": r["Area"],
                    "_intensity": r["Intensity"],
                    "_charge": r["Charge"],
                    "_rt": r["ExpectedRt"],
                    "_expectedRtWidth": r["ExpectedRtWidth"],
                    "_isQualifier": bool(r["IsQualifier"]),
                    "_isInternalStandard": bool(r["IsInternalStandard"]),
                    "_msMsExperimentIndex": r["MsMsExperimentIndex"],
                    "_msMsRetentionTime": r["MsMsRetentionTime"],
                    "_ms1Cycle": r["Ms1Cycle"],
                    "_msMsCycle": r["MsMsCycle"],
                    "_modificationText": r["ModificationText"],
                }
                lib_raw = r["LibrarySearchResult"]
                if lib_raw and lib_raw.strip():
                    entry["_librarySearchResults"] = self._v2_parse_library_result(lib_raw)
                self._v2_qualpeak_index[key] = entry

        return self._v2_qualpeak_index.get(
            (sample_index, compound_index),
            {"_containsMSMS": False, "_hasBeenCalculated": False},
        )

    @staticmethod
    def _v2_parse_library_result(raw: str) -> list[dict[str, Any]]:
        """Parse QualPeak.LibrarySearchResult pipe-delimited format.

        Format variants (pipe-delimited per line):
          5-field: count|fit|reverse_fit|guid|is_smart
          5-field: fit|reverse_fit|count|guid|is_smart
          3-field: fit|reverse_fit|purity
        """
        results = []
        for line in raw.strip().split("\n"):
            parts = line.split("|")
            if len(parts) < 3:
                continue
            # Determine which parts are floats vs GUIDs
            floats = []
            guid = ""
            is_smart = False
            for p in parts:
                p = p.strip()
                if not p:
                    continue
                try:
                    floats.append(float(p))
                except ValueError:
                    if "-" in p and len(p) > 30:
                        guid = p
                    elif p.lower() in ("true", "false"):
                        is_smart = p.lower() == "true"

            if len(floats) >= 2:
                results.append({
                    "_fit": floats[-3] if len(floats) >= 3 else floats[0],
                    "_reverseFit": floats[-2] if len(floats) >= 3 else floats[1],
                    "_purity": floats[-1] if len(floats) >= 3 else 0.0,
                    "_librarySearchResultId": guid,
                    "_isSmartConfirmation": is_smart,
                })
        return results

    # ── library GUID resolution ──────────────────────────────────

    @staticmethod
    def guid_dict_to_str(guid_dict: Any) -> str:
        """Convert a .NET-serialized ``System.Guid`` dict to hex string.

        Handles both the dict form (v1 RTParts) and a plain string (v2 QualPeak).
        """
        if isinstance(guid_dict, str):
            return guid_dict
        if isinstance(guid_dict, dict):
            a = guid_dict.get("_a", 0)
            b = guid_dict.get("_b", 0)
            c = guid_dict.get("_c", 0)
            rest = bytes(guid_dict.get(f"_{chr(ord('d') + i)}", 0) for i in range(8))
            return f"{a:08x}-{b:04x}-{c:04x}-{rest[:2].hex()}-{rest[2:].hex()}"
        return str(guid_dict) if guid_dict else ""

    def resolve_library_hits(
        self,
        library_db: str | Path,
    ) -> int:
        """Resolve ``_librarySearchResultId`` GUIDs to library entry names.

        Opens *library_db* (a ``libview_*.sqlite``) and populates
        ``_librarySearchResults`` entries with ``resolved_name``,
        ``resolved_formula``, ``resolved_cas`` keys looked up via::

            MassSpectrum.Id → CompoundId → CompoundName.Name

        Returns the number of hits resolved.
        """
        import sqlite3

        # Collect all unique GUID strings from cached XIC data
        guid_set: set[str] = set()

        # Check v2 cache
        xic_cache = getattr(self, "_v2_qualpeak_index", None) or {}
        for entry in xic_cache.values():
            for lr in (entry.get("_librarySearchResults") or []):
                g = self.guid_dict_to_str(lr.get("_librarySearchResultId", ""))
                if g:
                    lr["_librarySearchResultId"] = g
                    guid_set.add(g)

        # Check v1 multidata cache
        try:
            md = self._load_multidata()
            xic_lookup = md.get("xic_lookup") or {}
            for key, entry in xic_lookup.items():
                for lr in (entry.get("_librarySearchResults") or []):
                    g = self.guid_dict_to_str(lr.get("_librarySearchResultId", ""))
                    if g:
                        lr["_librarySearchResultId"] = g
                        guid_set.add(g)
        except Exception:
            pass

        if not guid_set:
            return 0

        lib = sqlite3.connect(str(library_db))
        try:
            cur = lib.cursor()
            # Build lookup: GUID → (name, formula, cas)
            lookup: dict[str, tuple[str, str, str]] = {}
            placeholders = ",".join("?" for _ in guid_set)
            cur.execute(
                f"""SELECT DISTINCT ms.Id, cn.Name, c.Formula, c.CAS
                    FROM MassSpectrum ms
                    JOIN CompoundName cn ON cn.CompoundId = ms.CompoundId
                    JOIN Compound c ON c.Id = ms.CompoundId
                    WHERE ms.Id IN ({placeholders})""",
                list(guid_set),
            )
            for row in cur.fetchall():
                lookup[row[0]] = (row[1] or "", row[2] or "", row[3] or "")

            # Apply to v2 cache
            resolved = 0
            for entry in xic_cache.values():
                for lr in (entry.get("_librarySearchResults") or []):
                    g = lr.get("_librarySearchResultId", "")
                    if g in lookup:
                        lr["resolved_name"] = lookup[g][0]
                        lr["resolved_formula"] = lookup[g][1]
                        lr["resolved_cas"] = lookup[g][2]
                        resolved += 1

            # Apply to v1 multidata cache
            try:
                for key, entry in xic_lookup.items():
                    for lr in (entry.get("_librarySearchResults") or []):
                        g = lr.get("_librarySearchResultId", "")
                        if g in lookup:
                            lr["resolved_name"] = lookup[g][0]
                            lr["resolved_formula"] = lookup[g][1]
                            lr["resolved_cas"] = lookup[g][2]
                            resolved += 1
            except Exception:
                pass
        finally:
            lib.close()
        return resolved



    def search_library(
        self,
        reader: Any,
        library_db: str | Path,
        *,
        ppm_tol: float = 50.0,
        dot_product_ppm: float = 20.0,
        prescreen_n: int = 200,
        top_n: int = 5,
        polarity_override: str | None = None,
    ) -> dict[tuple[int, int], list["LibraryHit"]]:
        """Search an external library for MS/MS spectra linked to this session.

        For every compound × sample peak that has an associated MS/MS
        spectrum, the spectrum is extracted from the wiff2 *reader* and
        searched against the *library_db* (a ``libview_*.sqlite`` file).

        Parameters
        ----------
        reader:
            A :class:`~pyx500r.reader.Wiff2Reader` (or compatible object
            with ``get_spectrum(sample_index, experiment_index, cycle_index)``
            and ``iter_spectra(sample_index, experiment_index)``) covering
            the wiff2 files referenced by this session.
        library_db:
            Path to a LibraryView ``.sqlite`` database (converted from LBP).
        ppm_tol:
            Precursor m/z window in ppm for library pre-screening.
        dot_product_ppm:
            Peak-matching tolerance in ppm for dot-product scoring.
        prescreen_n:
            Maximum candidates from signature-based pre-screening.
        top_n:
            Number of top-scoring library matches to return per peak.
        polarity_override:
            If set (``"POS"`` or ``"NEG"``), overrides the polarity
            detected from the MS/MS experiment. Useful when the wiff2
            experiment metadata is unreliable.

        Returns
        -------
        dict mapping ``(sample_index, compound_index)`` → ``list[LibraryHit]``.
        Each :class:`LibraryHit` carries ``score``, ``name``, ``formula``,
        ``cas`` in addition to the standard ``fit``/``reverse_fit`` fields.
        """
        from .libsearch import (
            LibrarySearcher,
        )
        import numpy as np

        # Open library once
        if isinstance(library_db, LibrarySearcher):
            searcher = library_db
            _own_searcher = False
        else:
            searcher = LibrarySearcher(str(library_db))
            _own_searcher = True

        try:
            results: dict[tuple[int, int], list[LibraryHit]] = {}

            # Build compound lookup: compound_index -> CompoundInfo
            # list_compounds returns compounds in index order
            compounds_list = self.list_compounds()
            compound_by_index: dict[int, Any] = {
                i: c for i, c in enumerate(compounds_list)
            }

            # Collect all peaks with MS/MS, grouped by sample
            # Each entry: (QuantPeakInfo, CompoundInfo, xic_result)
            msms_peaks: dict[
                int, dict[int, tuple[Any, Any, dict[str, Any]]]
            ] = {}
            for peak in self.iter_peaks():
                if not peak.xic_result:
                    continue
                xic = peak.xic_result
                if not xic.get("_containsMSMS"):
                    continue
                msms_exp = xic.get("_msMsExperimentIndex")
                if msms_exp is None:
                    continue
                ci = peak.compound_index
                if ci in compound_by_index:
                    msms_peaks.setdefault(peak.sample_index, {})[ci] = (
                        peak, compound_by_index[ci], xic
                    )

            if not msms_peaks:
                return results

            # For each sample, scan MS/MS experiments and find matching spectra
            for sample_index, peaks in msms_peaks.items():
                # Get all MS/MS experiment indices used in this sample
                msms_experiments: set[int] = set()
                for _peak, _compound, xic in peaks.values():
                    exp_idx = xic.get("_msMsExperimentIndex")
                    if exp_idx is not None:
                        msms_experiments.add(exp_idx)

                # For each MS/MS experiment, collect spectra near the
                # expected retention times
                spectrum_cache: dict[
                    int, dict[float, tuple[np.ndarray, np.ndarray, str]]
                ] = {}

                for exp_idx in msms_experiments:
                    spectrum_cache[exp_idx] = {}
                    try:
                        for spec in reader.iter_spectra(
                            sample_index=sample_index,
                            experiment_index=exp_idx,
                            return_arrays=True,
                        ):
                            if spec.precursor_mz is None:
                                continue
                            rt = spec.scan_time
                            spectrum_cache[exp_idx][rt] = (
                                np.asarray(spec.mz, dtype=np.float64),
                                np.asarray(spec.intensities, dtype=np.float64),
                                spec.precursor_mz,
                            )
                    except Exception:
                        continue

                # Now search each compound against the library
                for compound_index, (peak, compound, xic) in peaks.items():
                    exp_idx = xic.get("_msMsExperimentIndex")
                    if exp_idx is None or exp_idx not in spectrum_cache:
                        continue

                    target_rt = xic.get("_msMsRetentionTime") or xic.get("_foundAtRt") or 0.0
                    cache = spectrum_cache[exp_idx]

                    # Find closest MS/MS scan by retention time
                    if not cache:
                        continue
                    rts = sorted(cache.keys())
                    idx = np.searchsorted(rts, target_rt)
                    best_rt = target_rt
                    if idx == 0:
                        best_rt = rts[0]
                    elif idx >= len(rts):
                        best_rt = rts[-1]
                    else:
                        if abs(rts[idx] - target_rt) < abs(
                            rts[idx - 1] - target_rt
                        ):
                            best_rt = rts[idx]
                        else:
                            best_rt = rts[idx - 1]

                    mz_arr, int_arr, prec_mz = cache[best_rt]

                    # Determine polarity
                    if polarity_override:
                        pol = polarity_override
                    else:
                        try:
                            exp = reader._experiment_info(
                                sample_index, exp_idx
                            )
                            pol = exp.polarity
                        except Exception:
                            pol = "POS"

                    # Use precursor mass from compound if available
                    search_prec = float(prec_mz)
                    if search_prec <= 0 and hasattr(compound, 'precursor_mass'):
                        search_prec = float(compound.precursor_mass)

                    # Search
                    hits_raw = searcher.search(
                        mz_arr,
                        int_arr,
                        precursor_mz=search_prec,
                        polarity=pol,
                        ppm_tol=ppm_tol,
                        dot_product_ppm=dot_product_ppm,
                        prescreen_n=prescreen_n,
                        top_n=top_n,
                    )

                    results[(sample_index, compound_index)] = [
                        LibraryHit(
                            fit=h["score"],
                            reverse_fit=h["score"],
                            purity=h["score"],
                            name=h["name"],
                            formula=h["formula"],
                            cas=h["cas"],
                            score=h["score"],
                            precursor_mz=h["precursor_mz"],
                            collision_energy=h["collision_energy"],
                            num_peaks=h["num_peaks"],
                            spectrum_id=h["spectrum_id"],
                            compound_id=h["compound_id"],
                        )
                        for h in hits_raw
                    ]

            return results
        finally:
            if _own_searcher:
                searcher.close()

    def results_matrix(self) -> list[list[QuantPeakInfo]]:
        """Return the full samples × compounds peak matrix."""
        if self._is_v2:
            return self._results_matrix_v2()
        return self._results_matrix_v1()

    def _results_matrix_v2(self) -> list[list[QuantPeakInfo]]:
        samples = self.list_samples()
        compounds = self.list_compounds()
        ns = len(samples)
        nc = len(compounds)
        matrix: list[list[QuantPeakInfo | None]] = [
            [None] * nc for _ in range(ns)
        ]
        for peak in self.iter_peaks():
            si = peak.sample_index
            ci = peak.compound_index
            if 0 <= si < ns and 0 <= ci < nc:
                matrix[si][ci] = peak
        return matrix

    def _results_matrix_v1(self) -> list[list[QuantPeakInfo]]:
        data = self._load_multidata()
        xic = data.get("xic_lookup")
        matrix: list[list[QuantPeakInfo]] = []
        for s in data["samples"]:
            si = s["index"]
            row: list[QuantPeakInfo] = []
            for p in s["peaks"]:
                xr = xic.get((si, p["compound_index"])) if xic else None
                row.append(self._peak_info_from_dict(p, si, xr))
            matrix.append(row)
        return matrix

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

    @staticmethod
    def _unpack_doubles(blob: bytes) -> np.ndarray | Sequence[float]:
        if len(blob) % 8 != 0:
            return []
        # read-only view is safe — the sqlite blob buffer is immutable
        return np.frombuffer(blob, dtype=np.float64)
