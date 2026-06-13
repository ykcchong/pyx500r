/**
 * TypeScript interfaces for the pyx500r JSON API.
 *
 * These mirror the frozen dataclasses in `src/pyx500r/models.py` 1:1 and are
 * intended to be the wire contract between a pyx500r-backed FastAPI service
 * (see ../examples/server/) and a TypeScript front-end.
 *
 * Conventions
 * -----------
 * - Python `float | None`  -> `number | null`
 * - Python `str | None`    -> `string | null`
 * - Python `datetime`      -> ISO-8601 `string | null` (encode server-side)
 * - Python `Sequence[float]` (list or numpy) -> `number[]`
 *   (call `.tolist()` on numpy arrays before JSON encoding; see GUI_INTEGRATION.md)
 * - `dict[str, Any] | None` -> `Record<string, unknown> | null`
 *
 * Generated to match pyx500r 0.2.0.
 */

/* ────────────────────────────────────────────────────────────────────────
 * WIFF2 — acquisition metadata & spectra
 * ──────────────────────────────────────────────────────────────────────── */

/** A single sample within a .wiff2 acquisition. (models.SampleInfo) */
export interface SampleInfo {
  index: number;
  sample_id: string;
  name: string;
  source: string;
  start_timestamp: string | null;
}

/** One acquisition experiment (scan function). (models.ExperimentInfo) */
export interface ExperimentInfo {
  index: number;
  experiment_id: string;
  scan_type: string;
  ms_level: number;
  /** "positive" | "negative" | "unknown" */
  polarity: string;
  cycle_count: number;
}

/** A device/instrument record. (models.InstrumentInfo) */
export interface InstrumentInfo {
  sample_index: number;
  instrument_index: number;
  device_type: number | null;
  device_name: string | null;
  model_name: string | null;
  serial_number: string | null;
  is_mass_spectrometer: boolean;
}

/** Lightweight spectrum metadata — cheap to fetch. (models.SpectrumMetadata) */
export interface SpectrumMetadata {
  sample_index: number;
  experiment_index: number;
  cycle_index: number;
  scan_time: number;
  scan_type: string;
  ms_level: number;
  polarity: string;
  point_count: number;
  precursor_mz: number | null;
  isolation_target_mz: number | null;
  isolation_lower_offset: number | null;
  isolation_upper_offset: number | null;
}

/**
 * A full spectrum. (models.SpectrumData)
 * NOTE: the intensity field is `intensities` (plural). `mz` and `intensities`
 * are parallel arrays of equal length.
 */
export interface SpectrumData {
  sample_index: number;
  experiment_index: number;
  cycle_index: number;
  scan_time: number;
  mz: number[];
  intensities: number[];
  centroided: boolean;
  precursor_mz: number | null;
  isolation_target_mz: number | null;
  isolation_lower_offset: number | null;
  isolation_upper_offset: number | null;
}

/** A chromatogram (TIC/BPC). (models.Chromatogram) */
export interface Chromatogram {
  times: number[];
  intensities: number[];
  experiment_index: number | null;
  ms_level: number | null;
}

/* ────────────────────────────────────────────────────────────────────────
 * QSession — quantitation results
 * ──────────────────────────────────────────────────────────────────────── */

/** XIC metadata row. (models.XicInfo) */
export interface XicInfo {
  xic_id: string;
  sample_key: string;
  mz_lower: number;
  mz_upper: number;
  group_index: number;
  replicate_index: number;
  status: number;
}

/** An extracted-ion chromatogram. (models.XicChromatogram) */
export interface XicChromatogram {
  xic_id: string;
  sample_key: string;
  times: number[];
  intensities: number[];
  mz_lower: number;
  mz_upper: number;
  status: number;
}

/** A compound (target) definition. (models.CompoundInfo) */
export interface CompoundInfo {
  name: string;
  group_name: string;
  formula: string | null;
  charge_formula: string | null;
  adduct_formula: string | null;
  precursor_mass: number;
  fragment_mass: number;
  /** 0 = signal, 1 = MS/MS */
  extraction_type: number;
  period: number;
  experiment: number;
  mz_lower: number;
  mz_upper: number;
  is_analyte: boolean;
  is_reportable: boolean;
  is_non_targeted: boolean;
  is_summed: boolean;
  is_from_multi_period_data: boolean;
  isotope_index: number;
  expected_mw: number;
  units: string;
  comment: string | null;
  internal_std_name: string | null;
  regression_area: boolean | null;
  regression_type: number | null;
  regression_weighting: number | null;
  use_auto_regression: boolean | null;
  integration_parameters: Record<string, unknown> | null;
  acquisition_indices: number[] | null;
  summed_compounds: number[] | null;
  extraction_values1: number[] | null;
  extraction_values2: number[] | null;
}

/** A sample within a quantitation session. (models.QuantSampleInfo) */
export interface QuantSampleInfo {
  index: number;
  sample_name: string;
  sample_id: string;
  sample_type: number;
  sample_comment: string | null;
  dilution_factor: number;
  injection_volume: number;
  user_name: string | null;
  acq_method_name: string | null;
  instrument_name: string | null;
  instrument_serial_number: string | null;
  batch_name: string | null;
  barcode: string | null;
  scanned_barcode: string | null;
  autosampler_method_supports_barcode: boolean;
  sample_comparison: boolean;
  ms_method: string | null;
  lc_method: string | null;
  sample_signature: string | null;
  rack: string | null;
  plate: string | null;
  vial: string | null;
  /** ISO-8601 timestamp, or null. */
  acquisition_date: string | null;
}

/**
 * Raw per-peak XIC detail (models.QuantPeakInfo.xic_result).
 * Keys use the original .NET names. Values are JSON-native. Optional in the
 * UI — prefer the projected fields on UnifiedPeakDTO.
 */
export interface XicResult {
  _foundAtMass?: number | null;
  _foundAtRt?: number | null;
  _foundAtRtApex?: number | null;
  _area?: number | null;
  _intensity?: number | null;
  _baseMass?: number | null;
  _extractionMass?: number | null;
  _extractionWidth?: number | null;
  _hasBeenCalculated?: boolean;
  _containsMSMS?: boolean;
  _msMsExperimentIndex?: number | null;
  _msMsRetentionTime?: number | null;
  _ms1Cycle?: number | null;
  _msMsCycle?: number | null;
  _charge?: number | null;
  _rt?: number | null;
  _expectedRtWidth?: number | null;
  _isQualifier?: boolean;
  _isInternalStandard?: boolean;
  _modificationText?: string | null;
  _librarySearchResults?: Array<Record<string, unknown>>;
  [key: string]: unknown;
}

/** A peak integration result. (models.QuantPeakInfo) */
export interface QuantPeakInfo {
  sample_index: number;
  compound_index: number;
  peak_index: number;
  use_for_calibration: boolean;
  peak_comment: string | null;
  actual_concentration: number;
  failed_query: boolean;
  valid_integration: boolean;
  modified: boolean;
  retention_time: number;
  area: number;
  corrected_area: number;
  height: number;
  corrected_height: number;
  start_rt: number;
  start_y: number;
  end_rt: number;
  end_y: number;
  half_height_start_rt: number;
  half_height_end_rt: number;
  noise: number;
  signal_to_noise: number;
  profile_type: number;
  peak_type: number;
  apex_rt: number;
  apex_y: number;
  region_area: number;
  region_height: number;
  s_mrm_retention_time_shift: boolean;
  row_hidden: boolean;
  reportable: boolean;
  molecular_weight: number;
  original_area: number;
  override_experiment_index: number;
  points_across_baseline: number;
  points_across_half_height: number;
  integration_parameters: Record<string, unknown> | null;
  profile: number[] | null;
  custom_fields: Record<string, string> | null;
  custom_peak_fields: Record<string, string> | null;
  start_x5_pct_height: number;
  end_x5_pct_height: number;
  start_x10_pct_height: number;
  end_x10_pct_height: number;
  std_addn_actual_concentration: number;
  extracted_ms_ms: number | null;
  xic_result: XicResult | null;
  super_group_id: string | null;
}

/** One library-search match. (models.LibraryHit) */
export interface LibraryHit {
  fit: number;
  reverse_fit: number;
  purity: number;
  is_smart_confirmation: boolean;
  library_entry_id: string;
  name: string;
  formula: string;
  cas: string;
  score: number;
  precursor_mz: number;
  collision_energy: number;
  num_peaks: number;
  spectrum_id: string;
  compound_id: string;
}

/* ────────────────────────────────────────────────────────────────────────
 * Bridge / derived DTOs (not 1:1 dataclasses — defined by the server)
 * ──────────────────────────────────────────────────────────────────────── */

/** m/z extraction window for a compound. (bridge.ExtractionWindow) */
export interface ExtractionWindow {
  period: number;
  experiment: number;
  mz_center: number;
  mz_half_window: number;
  rt_start: number | null;
  rt_end: number | null;
}

/**
 * Flattened compound+peak+XIC row. Projected from `UnifiedPeak`
 * (which is a property view, not a dataclass). The exact shape is defined by
 * your server's `unified_to_dict`; this matches examples/server/serializers.py.
 */
export interface UnifiedPeakDTO {
  name: string;
  formula: string | null;
  sample_index: number;
  compound_index: number;
  area: number;
  retention_time: number;
  height: number;
  signal_to_noise: number;
  found_mass: number | null;
  found_rt: number | null;
  /** mass_error * 1e6, or null. */
  mass_error_ppm: number | null;
  isotope_diff: number | null;
  rt_diff: number | null;
  contains_msms: boolean;
  has_been_calculated: boolean;
  valid_integration: boolean;
  is_valid: boolean;
  library_hits: LibraryHit[];
}

/** Sample routing entry from bridge.match_samples(). */
export interface SampleRouting {
  qsession_index: number;
  qsession_sample: QuantSampleInfo;
  wiff_index: number | null;
  wiff_sample: SampleInfo | null;
}

/* ────────────────────────────────────────────────────────────────────────
 * Server-defined response envelopes (match examples/server/app.py)
 * ──────────────────────────────────────────────────────────────────────── */

/** GET /api/files */
export interface FileListing {
  wiff2: string[];
  qsession: string[];
}

/** GET /api/qsession/{file}/info */
export interface QSessionInfo {
  version: string | null;
  qmap_version: string | null;
  locked: boolean;
  sample_count: number;
  compound_count: number;
  xic_count: number;
}

/** Generic server-side pagination envelope. */
export interface Page<T> {
  total: number;
  page: number;
  size: number;
  items: T[];
}
