"""Shared, dependency-free data models for SCIEX WIFF2 readers.

These dataclasses are used by both the DLL-backed reader (``wiff2``) and the
pure-Python reader (``reader``) so the two implementations expose an identical
surface.
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
    xic_result: dict[str, Any] | None = None
    super_group_id: str | None = None


class UnifiedPeak:
    """A single compound result — peak + compound + XIC data in one flat namespace.

    Attribute lookup order:
    1. Peak fields (area, retention_time, …)
    2. Compound fields (name, formula, mz_lower, …)
    3. XIC fields from ``xic_result`` (found_mass, found_rt, library_hits, …)

    XIC fields are accessible without the ``_`` prefix
    (e.g. ``up.found_mass`` instead of ``up.xic_result[\"_foundAtMass\"]``).
    ``library_hits`` is a ``list[LibraryHit]`` instead of the raw
    ``_librarySearchResults`` list-of-dicts.
    """

    __slots__ = ("_peak", "_compound", "_xic")

    def __init__(
        self,
        peak: QuantPeakInfo | None,
        compound: CompoundInfo | None = None,
    ) -> None:
        self._peak = peak
        self._compound = compound
        self._xic: dict[str, Any] | None = peak.xic_result if peak is not None else None

    # -- direct access -------------------------------------------------

    @property
    def peak(self) -> QuantPeakInfo:
        return self._peak

    @property
    def compound(self) -> CompoundInfo | None:
        return self._compound

    @property
    def xic(self) -> dict[str, Any] | None:
        return self._xic

    # -- compound pass-through -----------------------------------------

    @property
    def name(self) -> str:
        return self._compound.name if self._compound else ""

    @property
    def formula(self) -> str | None:
        return self._compound.formula if self._compound else None

    @property
    def mz_lower(self) -> float:
        return self._compound.mz_lower if self._compound else 0.0

    @property
    def mz_upper(self) -> float:
        return self._compound.mz_upper if self._compound else 0.0

    @property
    def period(self) -> int:
        return self._compound.period if self._compound else 0

    @property
    def experiment(self) -> int:
        return self._compound.experiment if self._compound else 0

    @property
    def is_analyte(self) -> bool:
        return self._compound.is_analyte if self._compound else False

    @property
    def is_reportable(self) -> bool:
        return self._compound.is_reportable if self._compound else False

    @property
    def internal_std_name(self) -> str | None:
        return self._compound.internal_std_name if self._compound else None

    # -- peak pass-through ---------------------------------------------

    @property
    def sample_index(self) -> int:
        return self._peak.sample_index if self._peak is not None else -1

    @property
    def compound_index(self) -> int:
        return self._peak.compound_index if self._peak is not None else -1

    @property
    def area(self) -> float:
        return self._peak.area if self._peak is not None else 0.0

    @property
    def retention_time(self) -> float:
        return self._peak.retention_time if self._peak is not None else 0.0

    @property
    def height(self) -> float:
        return self._peak.height if self._peak is not None else 0.0

    @property
    def corrected_area(self) -> float:
        return self._peak.corrected_area if self._peak is not None else 0.0

    @property
    def corrected_height(self) -> float:
        return self._peak.corrected_height if self._peak is not None else 0.0

    @property
    def start_rt(self) -> float:
        return self._peak.start_rt if self._peak is not None else 0.0

    @property
    def end_rt(self) -> float:
        return self._peak.end_rt if self._peak is not None else 0.0

    @property
    def apex_rt(self) -> float:
        return self._peak.apex_rt if self._peak is not None else 0.0

    @property
    def apex_y(self) -> float:
        return self._peak.apex_y if self._peak is not None else 0.0

    @property
    def noise(self) -> float:
        return self._peak.noise if self._peak is not None else 0.0

    @property
    def signal_to_noise(self) -> float:
        return self._peak.signal_to_noise if self._peak is not None else 0.0

    @property
    def valid_integration(self) -> bool:
        return self._peak.valid_integration if self._peak is not None else False

    @property
    def failed_query(self) -> bool:
        return self._peak.failed_query if self._peak is not None else False

    @property
    def modified(self) -> bool:
        return self._peak.modified if self._peak is not None else False

    @property
    def row_hidden(self) -> bool:
        return self._peak.row_hidden if self._peak is not None else False

    @property
    def region_area(self) -> float:
        return self._peak.region_area if self._peak is not None else 0.0

    @property
    def region_height(self) -> float:
        return self._peak.region_height if self._peak is not None else 0.0

    @property
    def actual_concentration(self) -> float:
        return self._peak.actual_concentration if self._peak is not None else 0.0

    @property
    def profile_type(self) -> int:
        return self._peak.profile_type if self._peak is not None else 0

    @property
    def peak_type(self) -> int:
        return self._peak.peak_type if self._peak is not None else 0

    @property
    def molecular_weight(self) -> float:
        return self._peak.molecular_weight if self._peak is not None else 0.0

    @property
    def original_area(self) -> float:
        return self._peak.original_area if self._peak is not None else 0.0

    @property
    def half_height_start_rt(self) -> float:
        return self._peak.half_height_start_rt if self._peak is not None else 0.0

    @property
    def half_height_end_rt(self) -> float:
        return self._peak.half_height_end_rt if self._peak is not None else 0.0

    @property
    def use_for_calibration(self) -> bool:
        return self._peak.use_for_calibration if self._peak is not None else False

    @property
    def peak_comment(self) -> str | None:
        return self._peak.peak_comment if self._peak is not None else None

    # -- XIC pass-through (without _ prefix) ---------------------------

    @property
    def found_mass(self) -> float | None:
        return self._xic.get("_foundAtMass") if self._xic else None

    @property
    def found_rt(self) -> float | None:
        return self._xic.get("_foundAtRt") if self._xic else None

    @property
    def found_rt_apex(self) -> float | None:
        return self._xic.get("_foundAtRtApex") if self._xic else None

    @property
    def xic_area(self) -> float | None:
        return self._xic.get("_area") if self._xic else None

    @property
    def xic_intensity(self) -> float | None:
        return self._xic.get("_intensity") if self._xic else None

    @property
    def base_mass(self) -> float | None:
        return self._xic.get("_baseMass") if self._xic else None

    @property
    def extraction_mass(self) -> float | None:
        return self._xic.get("_extractionMass") if self._xic else None

    @property
    def extraction_width(self) -> float | None:
        return self._xic.get("_extractionWidth") if self._xic else None

    @property
    def mass_error(self) -> float | None:
        """Relative mass error: (found_mass - extraction_mass) / extraction_mass."""
        if not self._xic:
            return None
        fm = self._xic.get("_foundAtMass")
        em = self._xic.get("_extractionMass")
        if fm is None or em is None or em == 0:
            return None
        return (fm - em) / em

    @property
    def isotope_diff(self) -> float | None:
        """Isotope ratio difference from expected (raw value × 100)."""
        if not self._xic:
            return None
        val = self._xic.get("<IsoptopeRatioDiffFromExpected>k__BackingField")
        if val is None:
            return None
        return val * 100

    @property
    def rt_diff(self) -> float | None:
        """Retention time difference: found_rt - expected_rt."""
        if not self._xic:
            return None
        fr = self._xic.get("_foundAtRt")
        er = self._xic.get("_rt")
        if fr is None or er is None:
            return None
        return fr - er

    @property
    def has_been_calculated(self) -> bool:
        return bool(self._xic.get("_hasBeenCalculated")) if self._xic else False

    @property
    def is_qualifier(self) -> bool:
        return bool(self._xic.get("_isQualifier")) if self._xic else False

    @property
    def is_internal_standard(self) -> bool:
        return bool(self._xic.get("_isInternalStandard")) if self._xic else False

    @property
    def contains_msms(self) -> bool:
        return bool(self._xic.get("_containsMSMS")) if self._xic else False

    @property
    def msms_experiment(self) -> int | None:
        return self._xic.get("_msMsExperimentIndex") if self._xic else None

    @property
    def msms_retention_time(self) -> float | None:
        return self._xic.get("_msMsRetentionTime") if self._xic else None

    @property
    def ms1_cycle(self) -> int | None:
        return self._xic.get("_ms1Cycle") if self._xic else None

    @property
    def msms_cycle(self) -> int | None:
        return self._xic.get("_msMsCycle") if self._xic else None

    @property
    def expected_rt(self) -> float | None:
        return self._xic.get("_rt") if self._xic else None

    @property
    def expected_rt_width(self) -> float | None:
        return self._xic.get("_expectedRtWidth") if self._xic else None

    @property
    def charge(self) -> int | None:
        return self._xic.get("_charge") if self._xic else None

    @property
    def modification_text(self) -> str | None:
        return self._xic.get("_modificationText") if self._xic else None

    @property
    def found_rt_start(self) -> float | None:
        return self._xic.get("<FoundAtRtStart>k__BackingField") if self._xic else None

    @property
    def found_rt_end(self) -> float | None:
        return self._xic.get("<FoundAtRtEnd>k__BackingField") if self._xic else None

    @property
    def library_hits(self) -> list["LibraryHit"]:
        items = self._xic.get("_librarySearchResults") if self._xic else None
        if not items or not isinstance(items, list):
            return []
        hits: list[LibraryHit] = []
        for h in items:
            if h is None:
                continue
            # Resolve library_entry_id from both forms
            raw_id = h.get("_librarySearchResultId", "")
            if isinstance(raw_id, dict):
                from .qsession import QSessionReader
                entry_id = QSessionReader.guid_dict_to_str(raw_id)
            else:
                entry_id = str(raw_id) if raw_id else ""
            lib_name = h.get("resolved_name") or ""
            lib_formula = h.get("resolved_formula") or ""
            lib_cas = h.get("resolved_cas") or ""
            hits.append(LibraryHit(
                fit=h.get("_fit", 0.0),
                reverse_fit=h.get("_reverseFit", 0.0),
                purity=h.get("_purity", 0.0),
                is_smart_confirmation=bool(h.get("_isSmartConfirmation", False)),
                library_entry_id=entry_id,
                name=lib_name,
                formula=lib_formula,
                cas=lib_cas,
            ))
        return hits

    @property
    def has_library_been_searched(self) -> bool:
        return bool(
            self._xic.get("<HasLibraryBeenSearched>k__BackingField")
        ) if self._xic else False

    # -- convenience --------------------------------------------------

    def is_valid(self) -> bool:
        """True if this compound has a real integration result."""
        return self._peak.valid_integration and self.has_been_calculated

    def __repr__(self) -> str:
        name = self.name or "?"
        area = self.area or 0
        rt = self.retention_time or 0
        return (
            f"UnifiedPeak({name!r}, area={area:.1f}, rt={rt:.3f}, "
            f"calc={self.has_been_calculated})"
        )


@dataclass(frozen=True, slots=True)
class LibraryHit:
    """One library search match from the XIC result.

    Fields from the built-in MultiQuant search:
        ``fit``, ``reverse_fit``, ``purity``, ``is_smart_confirmation``.

    Fields added by :meth:`QSessionReader.search_library` (external library):
        ``name``, ``formula``, ``cas``, ``score``,
        ``precursor_mz``, ``collision_energy``, ``num_peaks``,
        ``spectrum_id``, ``compound_id``.
    """

    fit: float
    reverse_fit: float
    purity: float
    is_smart_confirmation: bool = False
    library_entry_id: str = ""

    # ── external library search fields ──
    name: str = ""
    formula: str = ""
    cas: str = ""
    score: float = 0.0
    precursor_mz: float = 0.0
    collision_energy: float = 0.0
    num_peaks: int = 0
    spectrum_id: str = ""
    compound_id: str = ""