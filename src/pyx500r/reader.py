"""Pure-Python reader for SCIEX X500R QTOF acquisitions.

Reads ``.wiff2`` + ``.wiff.scan`` files produced by the SCIEX X500R QTOF
system (and compatible SCIEX QTOF platforms). The reader decrypts the
embedded SQLite database (see :mod:`pyx500r.crypto`), parses acquisition
metadata, and decodes TOF spectra via :mod:`pyx500r.tof`.

Designed for **small-molecule screening**, **toxicology** and **forensic
toxicology** workflows where automated extraction of TOF-MS and MS/MS
spectra from X500R acquisitions is required.

Centroiding is provided by :mod:`pyx500r.centroid`.
"""

from __future__ import annotations

import sqlite3
from collections.abc import Iterator, Sequence
from datetime import datetime
from functools import lru_cache
from pathlib import Path

import numpy as np

from .crypto import WIFF2_PASSWORD, decrypt_database
from .models import (
    Chromatogram,
    ExperimentInfo,
    InstrumentInfo,
    SampleInfo,
    SpectrumData,
    SpectrumMetadata,
)
from .tof import decompress_tof

from .centroid import centroid_spectrum as _centroid_impl

# Offset of the contiguous record section inside a ``.wiff.scan`` file.
_SCAN_DATA_OFFSET = 44
_SCAN_SENTINEL = b"\xff\xff\xff\xff"

_POLARITY = {1: "positive", 2: "negative"}


def open_wiff2(path: str | Path, password: str = WIFF2_PASSWORD) -> "WiffReader":
    """Open a ``.wiff2`` file and return a :class:`WiffReader`."""
    return WiffReader(path, password=password)


class WiffReader:
    """Pure-Python reader over a single ``.wiff2`` acquisition."""

    def __init__(self, wiff_path: str | Path, password: str = WIFF2_PASSWORD):
        self.wiff_path = Path(wiff_path).resolve()
        if self.wiff_path.suffix.lower() != ".wiff2":
            raise ValueError(f"Only .wiff2 is supported, got: {self.wiff_path.name}")
        if not self.wiff_path.exists():
            raise FileNotFoundError(f"WIFF2 file does not exist: {self.wiff_path}")

        self._scan_path = self.wiff_path.with_suffix(".wiff.scan")
        if not self._scan_path.exists():
            raise FileNotFoundError(f"Companion scan file is missing: {self._scan_path}")

        raw_db = decrypt_database(self.wiff_path, password)
        self._conn = self._connect(raw_db)
        self._conn.row_factory = sqlite3.Row
        self._scan_bytes: bytes | None = None  # lazily memory-mapped/read
        self._prefetch_cache: dict[tuple[int, int, int], sqlite3.Row] = {}
        self._time_bins_cache: dict[tuple[int, int], int] = {}
        self._cycle_times_cache: dict[tuple[int, int], list[float]] = {}
        self._experiment_info_cache: dict[tuple[int, int], ExperimentInfo] = {}

    # ------------------------------------------------------------------ #
    # lifecycle
    # ------------------------------------------------------------------ #
    def __enter__(self) -> "WiffReader":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def close(self) -> None:
        if self._conn is not None:
            self._conn.close()
            self._conn = None  # type: ignore[assignment]
        self._scan_bytes = None

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

    @property
    def _scan(self) -> bytes:
        if self._scan_bytes is None:
            self._scan_bytes = self._scan_path.read_bytes()
        return self._scan_bytes

    # ------------------------------------------------------------------ #
    # samples
    # ------------------------------------------------------------------ #
    @property
    def sample_count(self) -> int:
        return self._conn.execute("SELECT COUNT(*) FROM sample").fetchone()[0]

    def list_samples(self) -> list[SampleInfo]:
        rows = self._conn.execute(
            "SELECT sample_id, sample_index, guid, name, file_path, created_date "
            "FROM sample ORDER BY sample_index"
        ).fetchall()
        return [
            SampleInfo(
                index=row["sample_index"],
                sample_id=str(row["guid"] or row["sample_id"]),
                name=str(row["name"]),
                source=str(row["file_path"] or ""),
                start_timestamp=str(row["created_date"]) if row["created_date"] else None,
            )
            for row in rows
        ]

    def sample_start_time(self, sample_index: int = 0) -> datetime | None:
        created = self._sample_row(sample_index)["created_date"]
        if not created:
            return None
        text = str(created).replace("Z", "+00:00")
        try:
            return datetime.fromisoformat(text)
        except ValueError:
            return None

    # ------------------------------------------------------------------ #
    # instruments
    # ------------------------------------------------------------------ #
    def list_instruments(self, sample_index: int = 0) -> list[InstrumentInfo]:
        rows = self._conn.execute(
            "SELECT device_type, device_model, device_model_name, serial_number "
            "FROM device_identifier ORDER BY device_identifier_id"
        ).fetchall()
        result: list[InstrumentInfo] = []
        for instrument_index, row in enumerate(rows):
            device_type = row["device_type"]
            result.append(
                InstrumentInfo(
                    sample_index=sample_index,
                    instrument_index=instrument_index,
                    device_type=device_type,
                    device_name=str(row["device_model_name"]) if row["device_model_name"] else None,
                    model_name=str(row["device_model"]) if row["device_model"] else None,
                    serial_number=str(row["serial_number"]) if row["serial_number"] else None,
                    is_mass_spectrometer=device_type == 0,
                )
            )
        return result

    def get_ms_instrument(self, sample_index: int = 0) -> InstrumentInfo | None:
        for instrument in self.list_instruments(sample_index):
            if instrument.is_mass_spectrometer:
                return instrument
        return None

    # ------------------------------------------------------------------ #
    # experiments
    # ------------------------------------------------------------------ #
    def get_experiments(self, sample_index: int = 0) -> list[ExperimentInfo]:
        sample_id = self._sample_row(sample_index)["sample_id"]
        rows = self._conn.execute(
            "SELECT experiment_id, name, polarity, idaType, start_cycle, stop_cycle "
            "FROM experimentRunInfo WHERE deviceRunInfo_id = ("
            "  SELECT deviceRunInfo_id FROM sample WHERE sample_id = ?) "
            "ORDER BY experiment_id",
            (sample_id,),
        ).fetchall()

        # Batch-fetch cycle counts to avoid N+1 queries
        exp_ids = [row["experiment_id"] for row in rows]
        count_map: dict[int, int] = {}
        if exp_ids:
            placeholders = ",".join("?" * len(exp_ids))
            count_rows = self._conn.execute(
                f"SELECT experiment_index, COUNT(*) as cnt FROM scanItems "
                f"WHERE sample_id = ? AND experiment_index IN ({placeholders}) "
                f"GROUP BY experiment_index",
                (sample_id, *exp_ids),
            ).fetchall()
            count_map = {r["experiment_index"]: r["cnt"] for r in count_rows}

        result: list[ExperimentInfo] = []
        for row in rows:
            exp_index = row["experiment_id"]
            result.append(
                ExperimentInfo(
                    index=exp_index,
                    experiment_id=str(exp_index),
                    scan_type=str(row["name"]),
                    ms_level=int(row["idaType"]) if row["idaType"] else 1,
                    polarity=_POLARITY.get(row["polarity"], "unknown"),
                    cycle_count=count_map.get(exp_index, 0),
                )
            )
        return result

    # ------------------------------------------------------------------ #
    # chromatograms
    # ------------------------------------------------------------------ #
    def get_experiment_tic(self, sample_index: int = 0, experiment_index: int = 0) -> Chromatogram:
        sample_id = self._sample_row(sample_index)["sample_id"]
        rows = self._conn.execute(
            "SELECT retentionTime, tic FROM scanItems "
            "WHERE sample_id = ? AND experiment_index = ? ORDER BY cycleIndex",
            (sample_id, experiment_index),
        ).fetchall()
        ms_level = self._experiment_ms_level(sample_index, experiment_index)
        return Chromatogram(
            times=[float(r["retentionTime"]) for r in rows],
            intensities=[float(r["tic"]) for r in rows],
            experiment_index=experiment_index,
            ms_level=ms_level,
        )

    def get_sample_tic(self, sample_index: int = 0) -> Chromatogram:
        from heapq import merge
        sequences = []
        for experiment in self.get_experiments(sample_index):
            tic = self.get_experiment_tic(sample_index=sample_index, experiment_index=experiment.index)
            sequences.append(zip(tic.times, tic.intensities))
        merged = list(merge(*sequences, key=lambda item: item[0]))
        return Chromatogram(
            times=[item[0] for item in merged],
            intensities=[item[1] for item in merged],
        )

    def prefetch_experiment(
        self,
        sample_index: int = 0,
        experiment_index: int | None = None,
    ) -> None:
        """Pre-fetch all scanItems for an experiment into memory.

        This eliminates per-spectrum SQLite round-trips when iterating over
        many cycles.  Call this before ``iter_spectra`` or repeated
        ``get_spectrum`` calls for the same experiment.

        If **experiment_index** is None, all experiments for the sample are
        pre-fetched.
        """
        sample_id = self._sample_row(sample_index)["sample_id"]

        if experiment_index is not None:
            rows = self._conn.execute(
                "SELECT * FROM scanItems "
                "WHERE sample_id = ? AND experiment_index = ? "
                "ORDER BY cycleIndex",
                (sample_id, experiment_index),
            ).fetchall()
            for row in rows:
                key = (sample_index, row["experiment_index"], row["cycleIndex"])
                self._prefetch_cache[key] = row
        else:
            rows = self._conn.execute(
                "SELECT * FROM scanItems WHERE sample_id = ? ORDER BY experiment_index, cycleIndex",
                (sample_id,),
            ).fetchall()
            for row in rows:
                key = (sample_index, row["experiment_index"], row["cycleIndex"])
                self._prefetch_cache[key] = row

    def clear_prefetch(self) -> None:
        """Clear the pre-fetch cache."""
        self._prefetch_cache.clear()

    def get_cycle_times(self, sample_index: int = 0, experiment_index: int = 0) -> list[float]:
        key = (sample_index, experiment_index)
        if key in self._cycle_times_cache:
            return self._cycle_times_cache[key]
        times = self.get_experiment_tic(
            sample_index=sample_index, experiment_index=experiment_index
        ).times
        self._cycle_times_cache[key] = times
        return times

    # ------------------------------------------------------------------ #
    # spectra
    # ------------------------------------------------------------------ #
    def get_spectrum(
        self,
        sample_index: int = 0,
        experiment_index: int = 0,
        cycle_index: int = 0,
        centroid: bool = False,
        return_arrays: bool = False,
    ) -> SpectrumData:
        row = self._scan_item(sample_index, experiment_index, cycle_index)
        mz, intensities = self._decode_spectrum(sample_index, experiment_index, row)

        if centroid:
            mz, intensities = _centroid_impl(mz, intensities, return_arrays=True)
        elif not return_arrays:
            # Ensure float dtype for consistency before converting to list
            if intensities.dtype != np.float64:
                intensities = intensities.astype(np.float64)

        if not return_arrays:
            # Convert to Python lists at the API boundary for backward compatibility
            mz = mz.tolist() if isinstance(mz, np.ndarray) else list(mz)
            intensities = intensities.tolist() if isinstance(intensities, np.ndarray) else list(intensities)

        precursor_mz = self._precursor_mz(row)
        target, lower, upper = self._isolation_window(row)
        return SpectrumData(
            sample_index=sample_index,
            experiment_index=experiment_index,
            cycle_index=cycle_index,
            scan_time=float(row["retentionTime"]),
            mz=mz,
            intensities=intensities,
            centroided=centroid,
            precursor_mz=precursor_mz,
            isolation_target_mz=target,
            isolation_lower_offset=lower,
            isolation_upper_offset=upper,
        )

    def get_spectrum_metadata(
        self,
        sample_index: int = 0,
        experiment_index: int = 0,
        cycle_index: int = 0,
    ) -> SpectrumMetadata:
        row = self._scan_item(sample_index, experiment_index, cycle_index)
        mz, _ = self._decode_spectrum(sample_index, experiment_index, row)
        experiment = self._experiment_info(sample_index, experiment_index)
        precursor_mz = self._precursor_mz(row)
        target, lower, upper = self._isolation_window(row)
        return SpectrumMetadata(
            sample_index=sample_index,
            experiment_index=experiment_index,
            cycle_index=cycle_index,
            scan_time=float(row["retentionTime"]),
            scan_type=experiment.scan_type,
            ms_level=experiment.ms_level,
            polarity=experiment.polarity,
            point_count=len(mz),
            precursor_mz=precursor_mz,
            isolation_target_mz=target,
            isolation_lower_offset=lower,
            isolation_upper_offset=upper,
        )

    def iter_spectra(
        self,
        sample_index: int = 0,
        experiment_index: int = 0,
        limit: int | None = None,
        return_arrays: bool = False,
    ) -> Iterator[SpectrumData]:
        # Auto-prefetch on first use to avoid per-cycle SQLite round-trips
        if (sample_index, experiment_index, 0) not in self._prefetch_cache:
            self.prefetch_experiment(sample_index, experiment_index)

        cycle_count = len(
            self.get_cycle_times(sample_index=sample_index, experiment_index=experiment_index)
        )
        max_cycles = cycle_count if limit is None else min(cycle_count, limit)
        for cycle_index in range(max_cycles):
            yield self.get_spectrum(
                sample_index=sample_index,
                experiment_index=experiment_index,
                cycle_index=cycle_index,
                return_arrays=return_arrays,
            )

    # ------------------------------------------------------------------ #
    # internals
    # ------------------------------------------------------------------ #
    @lru_cache(maxsize=8)
    def _sample_row(self, sample_index: int) -> sqlite3.Row:
        row = self._conn.execute(
            "SELECT * FROM sample WHERE sample_index = ?", (sample_index,)
        ).fetchone()
        if row is None:
            raise IndexError(f"sample_index out of range: {sample_index}")
        return row

    def _experiment_info(self, sample_index: int, experiment_index: int) -> ExperimentInfo:
        key = (sample_index, experiment_index)
        if key in self._experiment_info_cache:
            return self._experiment_info_cache[key]
        for experiment in self.get_experiments(sample_index):
            self._experiment_info_cache[(sample_index, experiment.index)] = experiment
            if experiment.index == experiment_index:
                return experiment
        raise IndexError(f"experiment_index out of range: {experiment_index}")

    def _experiment_ms_level(self, sample_index: int, experiment_index: int) -> int:
        return self._experiment_info(sample_index, experiment_index).ms_level

    def _time_bins_to_sum(self, sample_index: int, experiment_index: int) -> int:
        key = (sample_index, experiment_index)
        if key in self._time_bins_cache:
            return self._time_bins_cache[key]
        sample_id = self._sample_row(sample_index)["sample_id"]
        row = self._conn.execute(
            "SELECT timeBinsToSum FROM experimentRunInfo "
            "WHERE experiment_id = ? AND deviceRunInfo_id = ("
            "  SELECT deviceRunInfo_id FROM sample WHERE sample_id = ?)",
            (experiment_index, sample_id),
        ).fetchone()
        result = 1 if (row is None or row["timeBinsToSum"] is None) else (int(row["timeBinsToSum"]) or 1)
        self._time_bins_cache[key] = result
        return result

    def _scan_item(self, sample_index: int, experiment_index: int, cycle_index: int) -> sqlite3.Row:
        cache_key = (sample_index, experiment_index, cycle_index)
        if cache_key in self._prefetch_cache:
            row = self._prefetch_cache[cache_key]
            if row is None:
                raise IndexError(
                    f"No scan for sample={sample_index} experiment={experiment_index} cycle={cycle_index}"
                )
            return row
        sample_id = self._sample_row(sample_index)["sample_id"]
        row = self._conn.execute(
            "SELECT * FROM scanItems "
            "WHERE sample_id = ? AND experiment_index = ? AND cycleIndex = ?",
            (sample_id, experiment_index, cycle_index),
        ).fetchone()
        if row is None:
            raise IndexError(
                f"No scan for sample={sample_index} experiment={experiment_index} cycle={cycle_index}"
            )
        return row

    def _decode_spectrum(
        self, sample_index: int, experiment_index: int, row: sqlite3.Row
    ) -> tuple[np.ndarray, np.ndarray]:
        size = int(row["size"])
        if size <= 0:
            return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float64)
        start = _SCAN_DATA_OFFSET + int(row["offset"])
        record = self._scan[start : start + size]
        sentinel = record.find(_SCAN_SENTINEL)
        if sentinel < 0:
            return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float64)
        stream = record[sentinel:]
        n = self._time_bins_to_sum(sample_index, experiment_index)
        # Use calibrated decompression — m/z computed inside the JIT kernel
        cal_a = float(row["slope"])
        cal_t0 = float(row["delay"])
        tdc_res = float(row["tdcResolution"])
        return decompress_tof(
            stream,
            number_of_time_bins_to_sum=n,
            cal_a=cal_a,
            cal_t0=cal_t0,
            time_resolution=tdc_res,
            return_arrays=True,
        )

    @staticmethod
    def _precursor_mz(row: sqlite3.Row) -> float | None:
        value = row["precursorMass"]
        if value is None or float(value) <= 0.0:
            return None
        return float(value)

    @staticmethod
    def _isolation_window(row: sqlite3.Row) -> tuple[float | None, float | None, float | None]:
        center = row["centerMass"] if "centerMass" in row else None
        width = row["centerMassWidth"] if "centerMassWidth" in row else None
        if center is None or float(center) <= 0.0:
            return None, None, None
        if width is None or float(width) <= 0.0:
            return float(center), None, None
        half = float(width) / 2.0
        return float(center), half, half
