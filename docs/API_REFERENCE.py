"""API Reference — all public functions and classes in sciex_wiff.

This document covers the **three reader paths** available:

1. **Pure-Python reader** (`WiffReader`) — no DLL dependency
2. **DLL-backed reader** (`Wiff2File` via `Wiff2Api`) — requires pythonnet + SCIEX DLLs
3. **Precursor index** (`PrecursorIndex`) — bulk MS2 precursor search

Both readers (paths 1 & 2) expose the **same public surface** with identical
dataclass return types, so they are drop-in compatible.

Data files
----------
Reference acquisitions in ``data/encoded/``:

    * ``example_neg.wiff2`` — negative polarity, 1 sample, 11 experiments
    * ``example_pos.wiff2`` — positive polarity, 1 sample, 11 experiments

Each ``.wiff2`` has companion files:
    * ``*.wiff.scan`` — unencrypted protobuf TOF spectral data
    * ``*.timeseries.data`` — encrypted SQLite (timeseries, not yet decoded)

---

*Auto-generated:* ``2026-06-02``
"""

# ===========================================================================
# 1. MODULE: sciex_wiff.crypto
# ===========================================================================
"""
``sciex_wiff.crypto`` — AES-128-OFB page decryption for SCIEX SQLite containers.
================================================================================

Constants
---------

.. py:data:: WIFF2_PASSWORD
   :type: str

   The SCIEX wiff2 connection password (hardcoded UUID):
   ``"F90CA3B4-CC7B-4439-A479-2097CB8AE246"``.
   The AES-128 key is ``WIFF2_PASSWORD.encode("utf-8")[:16]``.

.. py:data:: PAGE_SIZE
   :type: int

   SQLite page size: ``4096`` bytes. Every page in ``.wiff2`` and
   ``.timeseries.data`` is this size.

.. py:data:: RESERVED_BYTES
   :type: int

   Reserved trailer at the end of each page: ``12`` bytes at offset
   ``4084:4096``. These hold the per-page random nonce in clear and
   are NOT encrypted.

Functions
---------

.. py:function:: decrypt_page(page_number: int, page: bytes, password: str = WIFF2_PASSWORD) -> bytes

   Decrypt a **single** 4096-byte SQLite page.

   :param page_number: 1-based SQLite page number.
   :param page: Exactly 4096 bytes of encrypted page data.
   :param password: Connection password (default: ``WIFF2_PASSWORD``).
   :returns: 4096 bytes of plaintext page data.
   :raises ValueError: if ``page`` is not exactly 4096 bytes.

   Special handling for page 1:
       * Bytes ``[0:16]`` are replaced with ``b"SQLite format 3\\x00"``.
       * Bytes ``[16:24]`` (page-size/reserved fields) are kept from ciphertext.

.. py:function:: decrypt_database(source: str | Path | bytes, password: str = WIFF2_PASSWORD) -> bytes

   Decrypt a **whole** encrypted SQLite database into plaintext bytes.

   :param source: File path or raw encrypted bytes.
   :param password: Connection password.
   :returns: Complete plaintext SQLite database bytes, openable with
             ``sqlite3.connect(":memory:").deserialize(result)``.
   :raises ValueError: if data size is not a multiple of ``PAGE_SIZE``.
"""

# ===========================================================================
# 2. MODULE: sciex_wiff.tof
# ===========================================================================
"""
``sciex_wiff.tof`` — TOF spectrum decompression & m/z calibration.
==================================================================

Functions
---------

.. py:function:: decompress_tof(stream: bytes, number_of_time_bins_to_sum: int = 1, min_bin: int = 0, cal_a: float | None = None, cal_t0: float | None = None, time_resolution: float | None = None) -> tuple[list[int | float], list[int]]

   Decode a compressed TOF stream into ``(mz_or_bins, intensities)``.

   :param stream: Compressed stream starting at the ``FF FF FF FF`` sentinel.
   :param number_of_time_bins_to_sum: Step size ``n`` (= ``timeBinsToSum`` from DB).
       1 for plain; typically 4 for instrument data.
   :param min_bin: Skip points below this time bin (default 0).
   :param cal_a: Slope from ``scanItems.slope``. If provided with `cal_t0` and
       `time_resolution`, returns **m/z values** instead of raw bins.
   :param cal_t0: Delay from ``scanItems.delay``.
   :param time_resolution: TDC resolution from ``scanItems.tdcResolution``.
   :returns:
       * If calibration params provided: ``(mz: list[float], intensities: list[int])``
       * Otherwise: ``(bins: list[int], intensities: list[int])``

   **Stream format** (reverse-engineered from ``Clearcore2.Compression.DecompressionAlgorithmTof``):

   * Optional ``FF FF FF FF`` fixed-bin-marker + ``uint32`` start bin.
   * Then variable-length RLE tokens until ``0xFF`` stop marker.
   * Each token: ``byte & 0x80`` → zero-run (gap); else → intensity.
   * Value encoding: ``≤123`` = literal; ``124`` = u8; ``125`` = u16 LE; ``126`` = u32 LE.

   Uses **numba** JIT for near-C performance; kernels are pre-compiled at import time.

Classes
-------

.. py:class:: TofCalibration(cal_a: float, cal_t0: float, time_resolution: float)

   Quadratic TOF m/z calibration (``Sciex.FMan.DefaultTofCalibration``).

   :param cal_a: Slope parameter (must be > 0).
   :param cal_t0: Delay parameter.
   :param time_resolution: TDC time resolution (must be > 0).

   **Formula**: ``m/z = (cal_a × time_resolution × bin − cal_a × cal_t0)²``

   .. py:method:: bin_to_mass(time_bin: float) -> float

       Convert a single time bin to m/z.

   .. py:method:: bins_to_masses(time_bins: list[int]) -> list[float]

       Convert a list of time bins to m/z values.

   .. py:method:: mass_to_bin(mass: float) -> float

       Convert m/z back to an approximate time bin (inverse calibration).
       Returns ``0.0`` for non-positive mass.
"""

# ===========================================================================
# 3. MODULE: sciex_wiff.models
# ===========================================================================
"""
``sciex_wiff.models`` — shared immutable dataclasses.
======================================================

All six dataclasses are ``@dataclass(frozen=True)`` and shared between
both reader implementations.

.. py:class:: SampleInfo

   :param index:           Zero-based sample index.
   :param sample_id:       GUID or sample identifier string.
   :param name:            Human-readable sample name.
   :param source:          Source file path.
   :param start_timestamp: ISO-8601 acquisition start time (or ``None``).

.. py:class:: ExperimentInfo

   :param index:        Zero-based experiment index.
   :param experiment_id: Experiment identifier string.
   :param scan_type:    Scan type name (e.g. ``"TOFMS"``).
   :param ms_level:     MS level (1 = MS1, 2 = MS2).
   :param polarity:     ``"positive"`` or ``"negative"``.
   :param cycle_count:  Number of cycles (scans) in this experiment.

.. py:class:: InstrumentInfo

   :param sample_index:       Owning sample index.
   :param instrument_index:   Zero-based instrument index.
   :param device_type:        Device type code (0 = mass spectrometer).
   :param device_name:        Human-readable device name.
   :param model_name:         Instrument model (e.g. ``"X500 QTOF"``).
   :param serial_number:      Instrument serial number.
   :param is_mass_spectrometer: ``True`` if ``device_type == 0``.

.. py:class:: Chromatogram

   :param times:             List of retention times (seconds).
   :param intensities:       List of corresponding intensities.
   :param experiment_index:  Source experiment index (optional).
   :param ms_level:          MS level (optional).

.. py:class:: SpectrumMetadata

   :param sample_index:          Owning sample index.
   :param experiment_index:      Owning experiment index.
   :param cycle_index:           Cycle (scan) number.
   :param scan_time:             Retention time (seconds).
   :param scan_type:             Scan type string.
   :param ms_level:              MS level.
   :param polarity:              ``"positive"`` or ``"negative"``.
   :param point_count:           Number of (m/z, intensity) pairs.
   :param precursor_mz:          Precursor m/z (or ``None``).
   :param isolation_target_mz:   Isolation window center (or ``None``).
   :param isolation_lower_offset: Lower half-width (or ``None``).
   :param isolation_upper_offset: Upper half-width (or ``None``).

.. py:class:: SpectrumData

   :param sample_index:          Owning sample index.
   :param experiment_index:      Owning experiment index.
   :param cycle_index:           Cycle (scan) number.
   :param scan_time:             Retention time (seconds).
   :param mz:                    List of m/z values.
   :param intensities:           List of intensity values.
   :param centroided:            ``True`` if centroid mode (always ``False`` for pure-Python).
   :param precursor_mz:          Precursor m/z (or ``None``).
   :param isolation_target_mz:   Isolation window center (or ``None``).
   :param isolation_lower_offset: Lower half-width (or ``None``).
   :param isolation_upper_offset: Upper half-width (or ``None``).
"""

# ===========================================================================
# 4. MODULE: sciex_wiff.reader (Pure-Python path)
# ===========================================================================
"""
``sciex_wiff.reader`` — pure-Python WIFF2 reader (NO .NET dependency).
======================================================================

Factory
-------

.. py:function:: open_wiff2(path: str | Path, password: str = WIFF2_PASSWORD) -> WiffReader

   Open a ``.wiff2`` file and return a fully-initialized :class:`WiffReader`.

Class
-----

.. py:class:: WiffReader(wiff_path: str | Path, password: str = WIFF2_PASSWORD)

   Pure-Python reader over a single ``.wiff2`` acquisition.
   Supports the context-manager protocol.

   **Properties**

   .. py:attribute:: sample_count
      :type: int

      Number of samples in the file (read-only).

   **Methods — lifecycle**

   .. py:method:: close() -> None

      Close the internal SQLite connection and release memory.

   **Methods — metadata**

   .. py:method:: list_samples() -> list[SampleInfo]

      Return metadata for all samples.

   .. py:method:: sample_start_time(sample_index: int = 0) -> datetime | None

      Return the acquisition start time for a sample.

   .. py:method:: list_instruments(sample_index: int = 0) -> list[InstrumentInfo]

      Return instrument details for a sample.

   .. py:method:: get_ms_instrument(sample_index: int = 0) -> InstrumentInfo | None

      Return the mass spectrometer instrument (``device_type == 0``), or ``None``.

   .. py:method:: get_experiments(sample_index: int = 0) -> list[ExperimentInfo]

      Return all experiments for a sample.

   **Methods — chromatograms**

   .. py:method:: get_experiment_tic(sample_index: int = 0, experiment_index: int = 0) -> Chromatogram

      Return the TIC (total ion chromatogram) for one experiment.

   .. py:method:: get_sample_tic(sample_index: int = 0) -> Chromatogram

      Return the merged TIC across all experiments for a sample.

   .. py:method:: get_cycle_times(sample_index: int = 0, experiment_index: int = 0) -> list[float]

      Return retention times for all cycles in an experiment.

   **Methods — spectra**

   .. py:method:: get_spectrum(sample_index: int = 0, experiment_index: int = 0, cycle_index: int = 0, centroid: bool = False) -> SpectrumData

      Return a single spectrum. ``centroid`` must be ``False`` (centroiding is not supported).

   .. py:method:: get_spectrum_metadata(sample_index: int = 0, experiment_index: int = 0, cycle_index: int = 0) -> SpectrumMetadata

      Return metadata for a single spectrum (no m/z or intensity data).

   .. py:method:: iter_spectra(sample_index: int = 0, experiment_index: int = 0, limit: int | None = None) -> Iterator[SpectrumData]

      Iterate over spectra in an experiment, optionally capped at ``limit``.
"""

# ===========================================================================
# 5. MODULE: sciex_wiff.wiff2 (DLL-backed path)
# ===========================================================================
"""
``sciex_wiff.wiff2`` — DLL-backed reader via SCIEX Data API.
=============================================================

This path requires **pythonnet** and the **SCIEX vendor DLLs** (``app-x64/`` or
``app/`` directory). It uses ``SCIEX.Apis.Data.v1.dll`` via .NET interop.

Factory
-------

.. py:class:: Wiff2Api(vendor_dir: str | Path | None = None, license_key: str | None = None)

   Initialize the SCIEX Data API runtime.

   :param vendor_dir: Directory containing ``SCIEX.Apis.Data.v1.dll``.
       Defaults to the ``proprietary/vendor/`` directory.
   :param license_key: SCIEX API license XML string. If ``None``, reads
       ``SCIEX_API_LICENSE_KEY`` from environment.

   **Methods**

   .. py:method:: open(wiff_path: str | Path) -> Wiff2File

      Open a ``.wiff2`` file and return a :class:`Wiff2File`.

Class
-----

.. py:class:: Wiff2File(api: Wiff2Api, wiff_path: str | Path)

   DLL-backed reader over a ``.wiff2`` file. Supports the context-manager protocol.

   **Public surface is identical to** :class:`WiffReader`:

   * :py:attr:`sample_count` (property)
   * :py:meth:`list_samples()`
   * :py:meth:`sample_start_time()`
   * :py:meth:`list_instruments()`
   * :py:meth:`get_ms_instrument()`
   * :py:meth:`get_experiments()`
   * :py:meth:`get_experiment_tic()`
   * :py:meth:`get_sample_tic()`
   * :py:meth:`get_cycle_times()`
   * :py:meth:`get_spectrum()` — **additionally supports** ``centroid=True`` and ``add_framing_zeros=True``
   * :py:meth:`get_spectrum_metadata()`
   * :py:meth:`iter_spectra()`
   * :py:meth:`close()`
"""

# ===========================================================================
# 6. MODULE: sciex_wiff.index (Precursor index path)
# ===========================================================================
"""
``sciex_wiff.index`` — precursor m/z index for fast MS2 lookup.
================================================================

Builder
-------

.. py:function:: build_precursor_index(wiff_path: Path) -> PrecursorIndex

   Build a sorted precursor m/z index from a single ``.wiff2`` file.
   Extracts all MS2 scan items with a valid precursor mass.

Class
-----

.. py:class:: PrecursorIndex

   Sorted precursor m/z index for O(log n) binary-search lookup.

   :param file_path:      Path to the source ``.wiff2`` file.
   :param precursor_mz:   ``np.ndarray[float64]`` sorted precursor m/z values.
   :param indices:        Aligned ``(sample_index, experiment_index, cycle_index)`` tuples.
   :param retention_times: Aligned retention times.
   :param n_ms2:          Number of MS2 scan items indexed.

   .. py:method:: find(target_mz: float, tolerance_da: float, ms_level: int = 2) -> list[tuple[int, int, int]]

      Binary search for precursors within ``tolerance_da`` of ``target_mz``.
      Returns ``(sample_index, experiment_index, cycle_index)`` tuples.

   .. py:method:: to_dict() -> dict

      Serialize the index to a JSON-compatible dict.

   .. py:method:: from_dict(d: dict) -> PrecursorIndex

      Deserialize an index from a dict (classmethod).
"""

# ===========================================================================
# 7. Complete API surface summary
# ===========================================================================
"""
Summary — all public symbols
=============================

.. code-block:: python

    # Data models (shared across all paths)
    from sciex_wiff import (
        SampleInfo,          # dataclass: sample metadata
        ExperimentInfo,      # dataclass: experiment metadata
        InstrumentInfo,      # dataclass: instrument metadata
        Chromatogram,        # dataclass: time/intensity trace
        SpectrumMetadata,    # dataclass: spectrum header
        SpectrumData,        # dataclass: m/z + intensity arrays
    )

    # Pure-Python reader (path 1 — NO DLL)
    from sciex_wiff import (
        WiffReader,          # class: pure-Python .wiff2 reader
        open_wiff2,          # function: factory for WiffReader
    )

    # DLL-backed reader (path 2 — requires pythonnet + SCIEX DLLs)
    from sciex_wiff import (
        Wiff2Api,            # class: SCIEX Data API factory
        Wiff2File,           # class: DLL-backed .wiff2 reader
    )

    # TOF codec & calibration
    from sciex_wiff import (
        TofCalibration,      # class: quadratic m/z calibration
        decompress_tof,      # function: decompress TOF stream → (mz, ints)
    )

    # Crypto primitives
    from sciex_wiff import (
        WIFF2_PASSWORD,      # str: "F90CA3B4-CC7B-4439-A479-2097CB8AE246"
        PAGE_SIZE,           # int: 4096
        RESERVED_BYTES,      # int: 12
        decrypt_page,        # function: decrypt single SQLite page
        decrypt_database,    # function: decrypt whole .wiff2 → bytes
    )

    # Precursor index (path 3 — bulk MS2 lookup)
    from sciex_wiff.index import (
        PrecursorIndex,      # class: sorted precursor m/z index
        build_precursor_index,  # function: build index from .wiff2
    )


Quick start — each path
------------------------

**Path 1 — Pure-Python (no DLL):**

.. code-block:: python

    from sciex_wiff import open_wiff2

    with open_wiff2("example_neg.wiff2") as reader:
        print(reader.sample_count)          # 1
        print(reader.list_samples()[0].name) # "example_neg"

        for exp in reader.get_experiments():
            print(f"  {exp.scan_type} (MS{exp.ms_level}, {exp.polarity})")

        tic = reader.get_experiment_tic(0, 0)
        spectrum = reader.get_spectrum(0, 0, 51)
        print(f"  spectrum: {len(spectrum.mz)} points")
        print(f"  first m/z: {spectrum.mz[0]:.7f}")


**Path 2 — DLL-backed (requires pythonnet + SCIEX license):**

.. code-block:: python

    from sciex_wiff import Wiff2Api

    api = Wiff2Api(vendor_dir="app-x64", license_key="<license_key>...</license_key>")
    with api.open("example_neg.wiff2") as reader:
        # Same API as path 1, plus centroid support
        spectrum = reader.get_spectrum(0, 0, 51, centroid=True)


**Path 3 — Precursor index:**

.. code-block:: python

    from sciex_wiff.index import build_precursor_index

    index = build_precursor_index(Path("example_neg.wiff2"))
    matches = index.find(target_mz=500.0, tolerance_da=0.01)
    for sample_idx, exp_idx, cycle_idx in matches:
        print(f"  sample={sample_idx} exp={exp_idx} cycle={cycle_idx}")
"""
