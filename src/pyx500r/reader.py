"""Pure-Python reader for SCIEX WIFF2 acquisitions.

This reader has **no .NET / DLL dependency**. It decrypts the ``.wiff2`` SQLite
database (see :mod:`pyx500r.crypto`), reads acquisition metadata from the
embedded tables, and decodes TOF spectra directly from the companion
``.wiff.scan`` file using the reverse-engineered codec in :mod:`pyx500r.tof`.

Centroiding is provided by :mod:`pyx500r.centroid` (ported from
``Clearcore2.RawXYProcessing.SpectralPeakFinder``).
"""

from __future__ import annotations

import sqlite3
from collections.abc import Iterator
from dataclasses import dataclass
from datetime import datetime
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


def _row_get(row: sqlite3.Row, key: str, default=None):
    return row[key] if key in row.keys() else default


def _positive_float(value) -> float | None:
    if value is None:
        return None
    result = float(value)
    return result if result > 0.0 else None


@dataclass(frozen=True, slots=True)
class _SampleRow:
    sample_id: int
    sample_index: int
    guid: str | None
    name: str
    file_path: str | None
    created_date: str | None
    device_run_info_id: int

    def __getitem__(self, key: str):
        mapping = {
            "sample_id": self.sample_id,
            "sample_index": self.sample_index,
            "guid": self.guid,
            "name": self.name,
            "file_path": self.file_path,
            "created_date": self.created_date,
            "deviceRunInfo_id": self.device_run_info_id,
        }
        return mapping[key]


@dataclass(frozen=True, slots=True)
class _DeviceRow:
    device_type: int | None
    device_model: str | None
    device_model_name: str | None
    serial_number: str | None


@dataclass(frozen=True, slots=True)
class _ScanRow:
    sample_index: int
    sample_id: int
    experiment_index: int
    cycle_index: int
    retention_time: float
    tic: float
    size: int
    offset: int
    slope: float
    delay: float
    tdc_resolution: float
    precursor_mass: float | None
    center_mass: float | None
    center_mass_width: float | None

    def __contains__(self, key: str) -> bool:
        return key in {
            "sample_id",
            "experiment_index",
            "cycleIndex",
            "retentionTime",
            "tic",
            "size",
            "offset",
            "slope",
            "delay",
            "tdcResolution",
            "precursorMass",
            "centerMass",
            "centerMassWidth",
        }

    def __getitem__(self, key: str):
        mapping = {
            "sample_id": self.sample_id,
            "experiment_index": self.experiment_index,
            "cycleIndex": self.cycle_index,
            "retentionTime": self.retention_time,
            "tic": self.tic,
            "size": self.size,
            "offset": self.offset,
            "slope": self.slope,
            "delay": self.delay,
            "tdcResolution": self.tdc_resolution,
            "precursorMass": self.precursor_mass,
            "centerMass": self.center_mass,
            "centerMassWidth": self.center_mass_width,
        }
        return mapping[key]


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

        self._scan_bytes: bytes | None = self._scan_path.read_bytes()
        self._scan_view: memoryview | None = memoryview(self._scan_bytes)
        self._time_bins_cache: dict[tuple[int, int], int] = {}
        self._cycle_times_cache: dict[tuple[int, int], list[float]] = {}
        self._experiment_info_cache: dict[tuple[int, int], ExperimentInfo] = {}
        self._point_count_cache: dict[tuple[int, int, int], int] = {}
        self._sample_rows: dict[int, _SampleRow] = {}
        self._sample_infos: list[SampleInfo] = []
        self._device_rows: list[_DeviceRow] = []
        self._instrument_cache: dict[int, list[InstrumentInfo]] = {}
        self._experiments_by_sample: dict[int, list[ExperimentInfo]] = {}
        self._scan_cache: dict[tuple[int, int, int], _ScanRow] = {}
        self._prefetch_cache: dict[tuple[int, int, int], _ScanRow] = {}
        self._scan_by_experiment: dict[tuple[int, int], list[_ScanRow]] = {}
        self._tic_cache: dict[tuple[int, int], Chromatogram] = {}

        raw_db = decrypt_database(self.wiff_path, password)
        conn = self._connect(raw_db)
        conn.row_factory = sqlite3.Row
        self._conn: sqlite3.Connection | None = conn
        try:
            self._load_native_metadata()
        finally:
            conn.close()
            self._conn = None

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
        self._scan_view = None
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
            self._scan_view = memoryview(self._scan_bytes)
        return self._scan_bytes

    @property
    def _scan_buffer(self) -> memoryview:
        if self._scan_view is None:
            self._scan_view = memoryview(self._scan)
        return self._scan_view

    def _load_native_metadata(self) -> None:
        """Materialize WIFF2 SQLite metadata into native Python indexes."""
        if self._conn is None:
            raise RuntimeError("WIFF2 metadata connection is not available")

        sample_rows = self._conn.execute(
            "SELECT sample_id, sample_index, guid, name, file_path, created_date, deviceRunInfo_id "
            "FROM sample ORDER BY sample_index"
        ).fetchall()
        sample_id_to_index: dict[int, int] = {}
        for row in sample_rows:
            sample = _SampleRow(
                sample_id=int(row["sample_id"]),
                sample_index=int(row["sample_index"]),
                guid=str(row["guid"]) if row["guid"] else None,
                name=str(row["name"]),
                file_path=str(row["file_path"]) if row["file_path"] else None,
                created_date=str(row["created_date"]) if row["created_date"] else None,
                device_run_info_id=int(row["deviceRunInfo_id"]),
            )
            self._sample_rows[sample.sample_index] = sample
            sample_id_to_index[sample.sample_id] = sample.sample_index
            self._sample_infos.append(
                SampleInfo(
                    index=sample.sample_index,
                    sample_id=str(sample.guid or sample.sample_id),
                    name=sample.name,
                    source=str(sample.file_path or ""),
                    start_timestamp=sample.created_date,
                )
            )

        self._device_rows = [
            _DeviceRow(
                device_type=row["device_type"],
                device_model=str(row["device_model"]) if row["device_model"] else None,
                device_model_name=str(row["device_model_name"]) if row["device_model_name"] else None,
                serial_number=str(row["serial_number"]) if row["serial_number"] else None,
            )
            for row in self._conn.execute(
                "SELECT device_type, device_model, device_model_name, serial_number "
                "FROM device_identifier ORDER BY device_identifier_id"
            ).fetchall()
        ]

        scan_rows = self._conn.execute(
            "SELECT * FROM scanItems ORDER BY sample_id, experiment_index, cycleIndex"
        ).fetchall()
        for row in scan_rows:
            sample_id = int(row["sample_id"])
            sample_index = sample_id_to_index.get(sample_id)
            if sample_index is None:
                continue
            scan = _ScanRow(
                sample_index=sample_index,
                sample_id=sample_id,
                experiment_index=int(row["experiment_index"]),
                cycle_index=int(row["cycleIndex"]),
                retention_time=float(row["retentionTime"]),
                tic=float(row["tic"] or 0.0),
                size=int(row["size"] or 0),
                offset=int(row["offset"] or 0),
                slope=float(row["slope"] or 0.0),
                delay=float(row["delay"] or 0.0),
                tdc_resolution=float(row["tdcResolution"] or 0.0),
                precursor_mass=_positive_float(_row_get(row, "precursorMass")),
                center_mass=_positive_float(_row_get(row, "centerMass")),
                center_mass_width=_positive_float(_row_get(row, "centerMassWidth")),
            )
            key = (scan.sample_index, scan.experiment_index, scan.cycle_index)
            self._scan_cache[key] = scan
            self._scan_by_experiment.setdefault(
                (scan.sample_index, scan.experiment_index), []
            ).append(scan)

        count_map = {
            key: len(rows)
            for key, rows in self._scan_by_experiment.items()
        }

        for sample in self._sample_rows.values():
            rows = self._conn.execute(
                "SELECT experiment_id, name, polarity, idaType, timeBinsToSum "
                "FROM experimentRunInfo WHERE deviceRunInfo_id = ? "
                "ORDER BY experiment_id",
                (sample.device_run_info_id,),
            ).fetchall()
            experiments: list[ExperimentInfo] = []
            for row in rows:
                exp_index = int(row["experiment_id"])
                time_bins = 1 if row["timeBinsToSum"] is None else (int(row["timeBinsToSum"]) or 1)
                self._time_bins_cache[(sample.sample_index, exp_index)] = time_bins
                info = ExperimentInfo(
                    index=exp_index,
                    experiment_id=str(exp_index),
                    scan_type=str(row["name"]),
                    ms_level=int(row["idaType"]) if row["idaType"] else 1,
                    polarity=_POLARITY.get(row["polarity"], "unknown"),
                    cycle_count=count_map.get((sample.sample_index, exp_index), 0),
                )
                experiments.append(info)
                self._experiment_info_cache[(sample.sample_index, exp_index)] = info
            self._experiments_by_sample[sample.sample_index] = experiments

    # ------------------------------------------------------------------ #
    # samples
    # ------------------------------------------------------------------ #
    @property
    def sample_count(self) -> int:
        return len(self._sample_infos)

    def list_samples(self) -> list[SampleInfo]:
        return list(self._sample_infos)

    def sample_start_time(self, sample_index: int = 0) -> datetime | None:
        created = self._sample_row(sample_index).created_date
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
        self._sample_row(sample_index)
        if sample_index not in self._instrument_cache:
            self._instrument_cache[sample_index] = [
                InstrumentInfo(
                    sample_index=sample_index,
                    instrument_index=instrument_index,
                    device_type=row.device_type,
                    device_name=row.device_model_name,
                    model_name=row.device_model,
                    serial_number=row.serial_number,
                    is_mass_spectrometer=row.device_type == 0,
                )
                for instrument_index, row in enumerate(self._device_rows)
            ]
        return list(self._instrument_cache[sample_index])

    def get_ms_instrument(self, sample_index: int = 0) -> InstrumentInfo | None:
        for instrument in self.list_instruments(sample_index):
            if instrument.is_mass_spectrometer:
                return instrument
        return None

    # ------------------------------------------------------------------ #
    # experiments
    # ------------------------------------------------------------------ #
    def get_experiments(self, sample_index: int = 0) -> list[ExperimentInfo]:
        self._sample_row(sample_index)
        return list(self._experiments_by_sample.get(sample_index, []))

    # ------------------------------------------------------------------ #
    # chromatograms
    # ------------------------------------------------------------------ #
    def get_experiment_tic(self, sample_index: int = 0, experiment_index: int = 0) -> Chromatogram:
        key = (sample_index, experiment_index)
        if key in self._tic_cache:
            return self._tic_cache[key]
        rows = self._scan_by_experiment.get(key, [])
        ms_level = self._experiment_ms_level(sample_index, experiment_index)
        tic = Chromatogram(
            times=[r.retention_time for r in rows],
            intensities=[r.tic for r in rows],
            experiment_index=experiment_index,
            ms_level=ms_level,
        )
        self._tic_cache[key] = tic
        return tic

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
        self._sample_row(sample_index)
        if experiment_index is not None:
            rows = self._scan_by_experiment.get((sample_index, experiment_index), [])
            for row in rows:
                key = (sample_index, row.experiment_index, row.cycle_index)
                self._prefetch_cache[key] = row
        else:
            for (scan_sample_index, _), rows in self._scan_by_experiment.items():
                if scan_sample_index != sample_index:
                    continue
                for row in rows:
                    key = (sample_index, row.experiment_index, row.cycle_index)
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
            scan_time=row.retention_time,
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
        point_count = self._decode_point_count(sample_index, experiment_index, row)
        experiment = self._experiment_info(sample_index, experiment_index)
        precursor_mz = self._precursor_mz(row)
        target, lower, upper = self._isolation_window(row)
        return SpectrumMetadata(
            sample_index=sample_index,
            experiment_index=experiment_index,
            cycle_index=cycle_index,
            scan_time=row.retention_time,
            scan_type=experiment.scan_type,
            ms_level=experiment.ms_level,
            polarity=experiment.polarity,
            point_count=point_count,
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
    def _sample_row(self, sample_index: int) -> _SampleRow:
        row = self._sample_rows.get(sample_index)
        if row is None:
            raise IndexError(f"sample_index out of range: {sample_index}")
        return row

    def _experiment_info(self, sample_index: int, experiment_index: int) -> ExperimentInfo:
        key = (sample_index, experiment_index)
        if key in self._experiment_info_cache:
            return self._experiment_info_cache[key]
        raise IndexError(f"experiment_index out of range: {experiment_index}")

    def _experiment_ms_level(self, sample_index: int, experiment_index: int) -> int:
        return self._experiment_info(sample_index, experiment_index).ms_level

    def _time_bins_to_sum(self, sample_index: int, experiment_index: int) -> int:
        key = (sample_index, experiment_index)
        if key in self._time_bins_cache:
            return self._time_bins_cache[key]
        self._experiment_info(sample_index, experiment_index)
        return 1

    def _scan_item(self, sample_index: int, experiment_index: int, cycle_index: int) -> _ScanRow:
        cache_key = (sample_index, experiment_index, cycle_index)
        if cache_key in self._prefetch_cache:
            return self._prefetch_cache[cache_key]
        row = self._scan_cache.get(cache_key)
        if row is None:
            raise IndexError(
                f"No scan for sample={sample_index} experiment={experiment_index} cycle={cycle_index}"
            )
        return row

    def _decode_point_count(
        self, sample_index: int, experiment_index: int, row: _ScanRow
    ) -> int:
        """Return the number of decoded points without computing the m/z axis.

        Skips the per-point quadratic m/z calibration done in
        :meth:`_decode_spectrum`, which is unnecessary for metadata.
        """
        key = (sample_index, experiment_index, row.cycle_index)
        if key in self._point_count_cache:
            return self._point_count_cache[key]
        if row.size <= 0:
            self._point_count_cache[key] = 0
            return 0
        start = _SCAN_DATA_OFFSET + row.offset
        stop = start + row.size
        sentinel = self._scan.find(_SCAN_SENTINEL, start, stop)
        if sentinel < 0:
            self._point_count_cache[key] = 0
            return 0
        stream = self._scan_buffer[sentinel:stop]
        n = self._time_bins_to_sum(sample_index, experiment_index)
        bins, _ = decompress_tof(
            stream,
            number_of_time_bins_to_sum=n,
            return_arrays=True,
        )
        point_count = len(bins)
        self._point_count_cache[key] = point_count
        return point_count

    def _decode_spectrum(
        self, sample_index: int, experiment_index: int, row: _ScanRow
    ) -> tuple[np.ndarray, np.ndarray]:
        if row.size <= 0:
            return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float64)
        start = _SCAN_DATA_OFFSET + row.offset
        stop = start + row.size
        sentinel = self._scan.find(_SCAN_SENTINEL, start, stop)
        if sentinel < 0:
            return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float64)
        stream = self._scan_buffer[sentinel:stop]
        n = self._time_bins_to_sum(sample_index, experiment_index)
        # Use calibrated decompression — m/z computed inside the JIT kernel
        return decompress_tof(
            stream,
            number_of_time_bins_to_sum=n,
            cal_a=row.slope,
            cal_t0=row.delay,
            time_resolution=row.tdc_resolution,
            return_arrays=True,
        )

    @staticmethod
    def _precursor_mz(row: _ScanRow) -> float | None:
        return row.precursor_mass

    @staticmethod
    def _isolation_window(row: _ScanRow) -> tuple[float | None, float | None, float | None]:
        center = row.center_mass
        width = row.center_mass_width
        if center is None:
            return None, None, None
        if width is None:
            return center, None, None
        half = width / 2.0
        return center, half, half
