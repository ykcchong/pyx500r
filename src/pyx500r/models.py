"""Shared, dependency-free data models for SCIEX X500R QTOF and MultiQuant readers.

These dataclasses are shared by the acquisition and qsession readers so the
two implementations expose an identical surface.
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime
from typing import Any

@dataclass(frozen=True, slots=True)
class SampleInfo:
    index: int
    sample_id: str
    name: str
    source: str
    start_timestamp: str | None


@dataclass(frozen=True, slots=True)
class ExperimentInfo:
    index: int
    experiment_id: str
    scan_type: str
    ms_level: int
    polarity: str
    cycle_count: int


@dataclass(frozen=True, slots=True)
class Chromatogram:
    times: list[float]
    intensities: list[float]
    experiment_index: int | None = None
    ms_level: int | None = None


@dataclass(frozen=True, slots=True)
class InstrumentInfo:
    sample_index: int
    instrument_index: int
    device_type: int | None
    device_name: str | None
    model_name: str | None
    serial_number: str | None
    is_mass_spectrometer: bool


@dataclass(frozen=True, slots=True)
class SpectrumMetadata:
    sample_index: int
    experiment_index: int
    cycle_index: int
    scan_time: float
    scan_type: str
    ms_level: int
    polarity: str
    point_count: int
    precursor_mz: float | None = None
    isolation_target_mz: float | None = None
    isolation_lower_offset: float | None = None
    isolation_upper_offset: float | None = None


@dataclass(frozen=True, slots=True)
class SpectrumData:
    sample_index: int
    experiment_index: int
    cycle_index: int
    scan_time: float
    mz: Sequence[float]
    intensities: Sequence[float]
    centroided: bool
    precursor_mz: float | None = None
    isolation_target_mz: float | None = None
    isolation_lower_offset: float | None = None
    isolation_upper_offset: float | None = None


@dataclass(frozen=True, slots=True)
class XicInfo:
    """Metadata for an extracted ion chromatogram (XIC) in a qsession file."""

    xic_id: str
    sample_key: str
    mz_lower: float
    mz_upper: float
    group_index: int
    replicate_index: int
    status: int


@dataclass(frozen=True, slots=True)
class XicChromatogram:
    """Extracted ion chromatogram (XIC) from a qsession file."""

    xic_id: str
    sample_key: str
    times: Sequence[float]
    intensities: Sequence[float]
    mz_lower: float
    mz_upper: float
    status: int


@dataclass(frozen=True, slots=True)
class CompoundInfo:
    """Compound definition from a qsession RTParts table."""

    name: str
    group_name: str
    formula: str | None
    charge_formula: str | None
    adduct_formula: str | None
    precursor_mass: float
    fragment_mass: float
    extraction_type: int
    period: int
    experiment: int
    mz_lower: float
    mz_upper: float
    is_analyte: bool
    is_reportable: bool
    is_non_targeted: bool
    is_summed: bool
    is_from_multi_period_data: bool
    isotope_index: int
    expected_mw: float
    units: str
    comment: str | None
    internal_std_name: str | None
    regression_area: bool | None
    regression_type: int | None
    regression_weighting: int | None
    use_auto_regression: bool | None
    integration_parameters: dict[str, Any] | None = None
    acquisition_indices: list[int] | None = None
    summed_compounds: list[int] | None = None
    extraction_values1: list[float] | None = None
    extraction_values2: list[float] | None = None
@dataclass(frozen=True, slots=True)
class QuantSampleInfo:
    """Sample metadata from a qsession RTParts quantitation results table."""

    index: int
    sample_name: str
    sample_id: str
    sample_type: int
    sample_comment: str | None
    dilution_factor: float
    injection_volume: float
    user_name: str | None
    acq_method_name: str | None
    instrument_name: str | None
    instrument_serial_number: str | None
    batch_name: str | None
    barcode: str | None
    scanned_barcode: str | None
    autosampler_method_supports_barcode: bool
    sample_comparison: bool
    ms_method: str | None
    lc_method: str | None
    sample_signature: str | None
    rack: str | None
    plate: str | None
    vial: str | None
    acquisition_date: datetime | None = None

@dataclass(frozen=True, slots=True)
class QuantPeakInfo:
    """Peak integration results for one compound in one sample."""

    sample_index: int
    compound_index: int
    peak_index: int
    use_for_calibration: bool
    peak_comment: str | None
    actual_concentration: float
    failed_query: bool
    valid_integration: bool
    modified: bool
    retention_time: float
    area: float
    corrected_area: float
    height: float
    corrected_height: float
    start_rt: float
    start_y: float
    end_rt: float
    end_y: float
    half_height_start_rt: float
    half_height_end_rt: float
    noise: float
    signal_to_noise: float
    profile_type: int
    peak_type: int
    apex_rt: float
    apex_y: float
    region_area: float
    region_height: float
    s_mrm_retention_time_shift: bool
    row_hidden: bool
    reportable: bool
    molecular_weight: float
    original_area: float
    override_experiment_index: int
    points_across_baseline: int
    points_across_half_height: int
    integration_parameters: dict[str, Any] | None = None
    profile: list[float] | None = None
    custom_fields: dict[str, str] | None = None
    custom_peak_fields: dict[str, str] | None = None
    start_x5_pct_height: float = 0.0
    end_x5_pct_height: float = 0.0
    start_x10_pct_height: float = 0.0
    end_x10_pct_height: float = 0.0
    std_addn_actual_concentration: float = 0.0
    extracted_ms_ms: float | None = None
    super_group_id: str | None = None
