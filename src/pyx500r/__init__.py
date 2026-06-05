"""pyx500r — Pure-Python reader for SCIEX X500R QTOF data.

* ``WiffReader`` reads X500R ``.wiff2`` + ``.wiff.scan`` acquisitions
  (TOF-MS and MS/MS) with zero .NET dependency.
* ``QSessionReader`` reads MultiQuant ``.qsession`` quantitation result
  tables (XICs, peak integration results, compound libraries).
* ``PrecursorIndex`` / ``build_precursor_index`` provide O(log n) binary
  search for batch precursor lookup across large toxicology/small-molecule
  file cohorts.

Quick start — X500R acquisition::

    from pyx500r import open_wiff2

    with open_wiff2("toxicology_run.wiff2") as reader:
        print(reader.list_samples())
        tic = reader.get_experiment_tic(experiment_index=0)
        spectrum = reader.get_spectrum(experiment_index=0, cycle_index=51)
        print(f"{len(spectrum.mz)} points, first m/z={spectrum.mz[0]:.7f}")

Quick start — MultiQuant qsession::

    from pyx500r import open_qsession

    with open_qsession("quant_results.qsession") as qs:
        for xic in qs.iter_xics():
            print(f"{xic.xic_id}: {len(xic.times)} points")
"""

from __future__ import annotations

# ── Always available (pure-Python) ──
from .centroid import Peak, add_framing_zeros, centroid_spectrum, moving_average_smooth
from .crypto import (
    PAGE_SIZE,
    QSESSION_PAGE_SIZE,
    QSESSION_PASSWORD,
    RESERVED_BYTES,
    WIFF2_PASSWORD,
    decrypt_database,
    decrypt_page,
)
from .models import (
    Chromatogram,
    CompoundInfo,
    ExperimentInfo,
    InstrumentInfo,
    QuantPeakInfo,
    QuantSampleInfo,
    SampleInfo,
    SpectrumData,
    SpectrumMetadata,
    XicChromatogram,
    XicInfo,
)
from .qsession import QSessionReader, open_qsession
from .reader import WiffReader, open_wiff2
from .rtparts import load_rtparts_stream, read_compounds
from .tof import (
    MassRange,
    TofCalibration,
    compress_tof,
    decompress_quad,
    decompress_tof,
    decompress_zero_width,
    mass_to_time,
    time_to_mass,
)

# ── Precursor index ──
from .index import PrecursorIndex, build_precursor_index

__all__ = [
    # data models
    "Chromatogram",
    "CompoundInfo",
    "ExperimentInfo",
    "InstrumentInfo",
    "QuantPeakInfo",
    "QuantSampleInfo",
    "SampleInfo",
    "SpectrumData",
    "SpectrumMetadata",
    "XicChromatogram",
    "XicInfo",
    "WiffReader",
    "open_wiff2",
    "QSessionReader",
    "open_qsession",
    "load_rtparts_stream",
    "read_compounds",
    # TOF decompression / compression
    "TofCalibration",
    "decompress_tof",
    "compress_tof",
    "decompress_quad",
    "decompress_zero_width",
    "mass_to_time",
    "time_to_mass",
    "MassRange",
    # centroiding
    "centroid_spectrum",
    "add_framing_zeros",
    "moving_average_smooth",
    "Peak",
    # decryption primitives
    "WIFF2_PASSWORD",
    "QSESSION_PASSWORD",
    "PAGE_SIZE",
    "QSESSION_PAGE_SIZE",
    "RESERVED_BYTES",
    "decrypt_database",
    "decrypt_page",
    # precursor index (requires numpy)
    "PrecursorIndex",
    "build_precursor_index",
]
