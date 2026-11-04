use mzdata::io::SpectrumSource;
use mzdata::params::ParamDescribed;
use mzdata::prelude::SpectrumLike;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashSet};

use crate::envelope::evaluate_isotope_envelope_peaks;
use crate::mzml::{open_reader, scan_token};
use crate::select::header_total_ion_current;
use crate::{
    query_candidate_comparison_mzml, query_fragment_match_mzml, query_isotope_envelope_mzml,
    query_select_mzml, CandidateComparisonQuery, CandidateComparisonResult, FragmentMatchQuery,
    FragmentMatchResult, IsotopeEnvelopeQuery, IsotopeEnvelopeResult, MassSpecError,
    PrecursorMetadata, SelectQueryResult, SelectedSpectrum,
};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", tag = "kind", content = "value")]
pub enum SpectrumSelector {
    Index(usize),
    Id(String),
}

/// Metadata selection shared by extraction, pooling and fragment measurements.
/// Explicit ids are resolved to native ids before intersections are applied.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SpectrumSelection {
    #[serde(default)]
    pub scans: Vec<String>,
    pub ms_level: Option<u8>,
    pub rt_min: Option<f64>,
    pub rt_max: Option<f64>,
    pub scan_min: Option<u64>,
    pub scan_max: Option<u64>,
}

pub fn select_spectrum_ids_reader(
    reader: &mut crate::Reader,
    selection: &SpectrumSelection,
) -> Result<Vec<String>, MassSpecError> {
    if selection.ms_level == Some(0)
        || selection.rt_min.is_some_and(|v| !v.is_finite())
        || selection.rt_max.is_some_and(|v| !v.is_finite())
        || selection
            .rt_min
            .zip(selection.rt_max)
            .is_some_and(|(a, b)| a > b)
        || selection
            .scan_min
            .zip(selection.scan_max)
            .is_some_and(|(a, b)| a > b)
    {
        return Err(MassSpecError::InvalidQuery(
            "invalid spectrum selection bounds".into(),
        ));
    }
    let matches = |id: &str, level: u8, rt: f64| {
        if selection.ms_level.is_some_and(|expected| expected != level) {
            return false;
        }
        if (selection.rt_min.is_some() || selection.rt_max.is_some())
            && (!rt.is_finite()
                || selection.rt_min.is_some_and(|min| rt < min)
                || selection.rt_max.is_some_and(|max| rt > max))
        {
            return false;
        }
        if selection.scan_min.is_some() || selection.scan_max.is_some() {
            let scan = scan_token(id).and_then(|v| v.parse::<u64>().ok());
            if scan.is_none_or(|scan| {
                selection.scan_min.is_some_and(|min| scan < min)
                    || selection.scan_max.is_some_and(|max| scan > max)
            }) {
                return false;
            }
        }
        true
    };
    let mut requested = HashSet::new();
    let mut selected = Vec::new();
    if !selection.scans.is_empty() {
        for scan in &selection.scans {
            let id = if scan.chars().all(|c| c.is_ascii_digit()) && !scan.is_empty() {
                format!("scan={scan}")
            } else {
                scan.clone()
            };
            let selector =
                crate::fragment::resolve_selector_from_reader(reader, &SpectrumSelector::Id(id))?;
            let spectrum = match &selector {
                SpectrumSelector::Index(i) => reader.get_spectrum_by_index(*i),
                SpectrumSelector::Id(id) => reader.get_spectrum_by_id(id),
            }
            .ok_or_else(|| MassSpecError::SpectrumNotFound(scan.clone()))?;
            if !requested.insert(spectrum.id().to_owned()) {
                return Err(MassSpecError::InvalidQuery(
                    "duplicate scan in selection".into(),
                ));
            }
            if matches(spectrum.id(), spectrum.ms_level(), spectrum.start_time()) {
                selected.push(spectrum.id().to_owned());
            }
        }
    } else {
        for spectrum in reader.iter() {
            if matches(spectrum.id(), spectrum.ms_level(), spectrum.start_time()) {
                selected.push(spectrum.id().to_owned());
            }
        }
    }
    if selected.is_empty() {
        return Err(MassSpecError::SpectrumNotFound(
            "no spectra satisfy the selection".into(),
        ));
    }
    Ok(selected)
}

#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum PeakOrder {
    #[default]
    MzAscending,
    IntensityDescending,
    IntensityCappedMzAscending,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ScanQuery {
    pub selector: SpectrumSelector,
    pub mz_min: Option<f64>,
    pub mz_max: Option<f64>,
    pub min_intensity: Option<f32>,
    pub max_peaks: Option<usize>,
    #[serde(default)]
    pub peak_order: PeakOrder,
    /// MS level gate for the scans mode: None matches any level, Some(level)
    /// sends scans that exist only at another level to missingScans.
    #[serde(default)]
    pub ms_level: Option<u8>,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase", tag = "unit", content = "value")]
pub enum MassTolerance {
    Da(f64),
    Ppm(f64),
}

impl MassTolerance {
    fn contains(self, query: f64, observed: f64) -> bool {
        let width = match self {
            Self::Da(value) => value,
            Self::Ppm(value) => query * value / 1_000_000.0,
        };
        (observed - query).abs() <= width
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ProductQuery {
    pub mz: f64,
    pub tolerance: MassTolerance,
    pub max_matches: Option<usize>,
    /// Restrict the census to these scan tokens (bare numbers or full
    /// nativeIDs). None censuses every MS2 spectrum in the file.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub scans: Option<Vec<String>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum MassSpecQuery {
    Index {
        max_spectra: Option<usize>,
    },
    Scan {
        selector: SpectrumSelector,
        mz_min: Option<f64>,
        mz_max: Option<f64>,
        min_intensity: Option<f32>,
        max_peaks: Option<usize>,
        #[serde(default)]
        peak_order: PeakOrder,
    },
    Scans {
        scans: Vec<String>,
        mz_min: Option<f64>,
        mz_max: Option<f64>,
        min_intensity: Option<f32>,
        max_peaks: Option<usize>,
        #[serde(default)]
        peak_order: PeakOrder,
    },
    Sum {
        scans: Vec<String>,
        bin_ppm: f64,
        mz_min: Option<f64>,
        mz_max: Option<f64>,
        min_intensity: Option<f32>,
        max_peaks: Option<usize>,
    },
    Window {
        ms_level: u8,
        rt_min: f64,
        rt_max: f64,
        max_scans: usize,
        bin_ppm: f64,
        mz_min: Option<f64>,
        mz_max: Option<f64>,
        min_intensity: Option<f32>,
        max_peaks: Option<usize>,
    },
    Xic {
        mz: f64,
        tolerance_ppm: f64,
        ms_level: u8,
        rt_min: Option<f64>,
        rt_max: Option<f64>,
    },
    Product {
        mz: f64,
        tolerance: MassTolerance,
        max_matches: Option<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        scans: Option<Vec<String>>,
    },
    Precursor {
        mz: f64,
        tolerance: MassTolerance,
        rt_min: Option<f64>,
        rt_max: Option<f64>,
        charge: Option<i32>,
        max_matches: Option<usize>,
    },
    FragmentMatch {
        selector: SpectrumSelector,
        sequence: String,
        #[serde(default)]
        modifications: Vec<crate::MassModification>,
        ion_types: Vec<crate::IonType>,
        max_charge: u8,
        tolerance_ppm: f64,
        min_intensity: Option<f32>,
        max_matches: Option<usize>,
        #[serde(default)]
        neutral_losses: Vec<crate::NamedNeutralLoss>,
        #[serde(default)]
        engine_neutral_losses: bool,
    },
    CandidateComparison {
        selector: SpectrumSelector,
        sequence: String,
        candidate_models: Vec<crate::CandidateModel>,
        ion_types: Vec<crate::IonType>,
        max_charge: u8,
        tolerance_ppm: f64,
        min_intensity: Option<f32>,
    },
    IsotopeEnvelope {
        mass: Option<f64>,
        proforma: Option<String>,
        #[serde(default)]
        external_mods: Vec<crate::proforma::ExternalModification>,
        charge: u8,
        rt_min: f64,
        rt_max: f64,
        tolerance_ppm: f64,
        isotopes: usize,
        max_offset: i32,
    },
    Select {
        statement: String,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Peak {
    pub mz: f64,
    pub intensity: f32,
}

impl From<(f64, f32)> for Peak {
    fn from((mz, intensity): (f64, f32)) -> Self {
        Self { mz, intensity }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ScanQueryResult {
    pub id: String,
    pub index: usize,
    pub ms_level: u8,
    pub centroided: bool,
    pub retention_time: f64,
    pub retention_time_unit: String,
    pub precursor: Option<ScanPrecursorMetadata>,
    pub total_peak_count: usize,
    pub matched_peak_count: usize,
    pub truncated: bool,
    pub peaks: Vec<Peak>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ScanPrecursorMetadata {
    pub selected_mz: f64,
    pub intensity: f32,
    pub charge: Option<i32>,
    pub isolation_target: f32,
    pub isolation_lower_bound: f32,
    pub isolation_upper_bound: f32,
    pub activation: Option<String>,
    pub collision_energy: Option<f32>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct IndexQueryResult {
    pub spectra_in_file: usize,
    pub spectra_inspected: usize,
    pub complete: bool,
    pub index_source: String,
    pub spectra_per_ms_level: BTreeMap<String, usize>,
    pub activation_types: BTreeMap<String, usize>,
    pub retention_time_min: Option<f64>,
    pub retention_time_max: Option<f64>,
    pub retention_time_unit: String,
    pub native_id_example: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct BatchSpectrum {
    pub scan: String,
    pub spectrum: ScanQueryResult,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ScansQueryResult {
    pub requested: usize,
    pub returned: usize,
    pub missing_scans: Vec<String>,
    pub spectra: Vec<BatchSpectrum>,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct MergedPeak {
    pub mz: f64,
    pub intensity: f32,
    pub merged_from: usize,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SumQueryResult {
    pub spectra_metadata: Vec<crate::SpectrumMetadata>,
    pub bin_ppm: f64,
    pub truncated: bool,
    pub scan_numbers: Vec<String>,
    /// Per-scan precursor metadata, positional with `scan_numbers`. Carries the
    /// isolation window so callers (e.g. the fragment-spectrum renderer) can
    /// exclude the unfragmented precursor region without a second mzML read —
    /// the precursor data is already produced by the scan traversal the sum
    /// runs internally.
    #[serde(default)]
    pub precursors: Vec<Option<ScanPrecursorMetadata>>,
    pub peak_count: usize,
    pub peaks: Vec<MergedPeak>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct WindowQueryResult {
    pub spectra_metadata: Vec<crate::SpectrumMetadata>,
    pub bin_ppm: f64,
    pub scans_summed: usize,
    pub scan_numbers: Vec<String>,
    pub retention_time_unit: String,
    pub truncated: bool,
    pub peak_count: usize,
    pub peaks: Vec<MergedPeak>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct XicPoint {
    pub retention_time: f64,
    pub scan: String,
    pub intensity: f32,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct XicQueryResult {
    pub point_count: usize,
    pub retention_time_unit: String,
    pub apex: Option<XicPoint>,
    pub points: Vec<XicPoint>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ProductMatch {
    pub id: String,
    pub index: usize,
    pub ms_level: u8,
    pub mz: f64,
    pub intensity: f32,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ProductQueryResult {
    pub query: ProductQuerySummary,
    pub scanned_spectrum_count: usize,
    pub matches: Vec<ProductMatch>,
    pub truncated: bool,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PrecursorQuery {
    pub mz: f64,
    pub tolerance: MassTolerance,
    pub rt_min: Option<f64>,
    pub rt_max: Option<f64>,
    pub charge: Option<i32>,
    pub max_matches: Option<usize>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PrecursorQueryResult {
    pub query: PrecursorQuerySummary,
    pub scanned_spectrum_count: usize,
    pub matches: Vec<SelectedSpectrum>,
    pub truncated: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PrecursorQuerySummary {
    pub mz: f64,
    pub tolerance: MassTolerance,
    pub rt_min: Option<f64>,
    pub rt_max: Option<f64>,
    pub charge: Option<i32>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ProductQuerySummary {
    pub mz: f64,
    pub tolerance: MassTolerance,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub scans: Option<Vec<String>>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum MassSpecQueryResult {
    Index(IndexQueryResult),
    Scan(ScanQueryResult),
    Scans(ScansQueryResult),
    Sum(SumQueryResult),
    Window(WindowQueryResult),
    Xic(XicQueryResult),
    Product(ProductQueryResult),
    Precursor(PrecursorQueryResult),
    FragmentMatch(FragmentMatchResult),
    CandidateComparison(CandidateComparisonResult),
    IsotopeEnvelope(IsotopeEnvelopeResult),
    Select(SelectQueryResult),
}

pub fn execute_query_mzml(
    bytes: &[u8],
    query: MassSpecQuery,
) -> Result<MassSpecQueryResult, MassSpecError> {
    match query {
        MassSpecQuery::Index { max_spectra } => {
            query_index_mzml(bytes, max_spectra).map(MassSpecQueryResult::Index)
        }
        MassSpecQuery::Scan {
            selector,
            mz_min,
            mz_max,
            min_intensity,
            max_peaks,
            peak_order,
        } => query_mzml(
            bytes,
            &ScanQuery {
                selector,
                mz_min,
                mz_max,
                min_intensity,
                max_peaks,
                peak_order,
                ms_level: None,
            },
        )
        .map(MassSpecQueryResult::Scan),
        MassSpecQuery::Scans {
            scans,
            mz_min,
            mz_max,
            min_intensity,
            max_peaks,
            peak_order,
        } => query_scans_mzml(
            bytes,
            &scans,
            &ScanQuery {
                selector: SpectrumSelector::Index(0),
                mz_min,
                mz_max,
                min_intensity,
                max_peaks,
                peak_order,
                ms_level: None,
            },
        )
        .map(MassSpecQueryResult::Scans),
        MassSpecQuery::Sum {
            scans,
            bin_ppm,
            mz_min,
            mz_max,
            min_intensity,
            max_peaks,
        } => query_sum_mzml(
            bytes,
            &scans,
            bin_ppm,
            mz_min,
            mz_max,
            min_intensity,
            max_peaks,
        )
        .map(MassSpecQueryResult::Sum),
        MassSpecQuery::Window {
            ms_level,
            rt_min,
            rt_max,
            max_scans,
            bin_ppm,
            mz_min,
            mz_max,
            min_intensity,
            max_peaks,
        } => query_window_mzml(
            bytes,
            ms_level,
            rt_min,
            rt_max,
            max_scans,
            bin_ppm,
            mz_min,
            mz_max,
            min_intensity,
            max_peaks,
        )
        .map(MassSpecQueryResult::Window),
        MassSpecQuery::Xic {
            mz,
            tolerance_ppm,
            ms_level,
            rt_min,
            rt_max,
        } => query_xic_mzml(bytes, mz, tolerance_ppm, ms_level, rt_min, rt_max)
            .map(MassSpecQueryResult::Xic),
        MassSpecQuery::Product {
            mz,
            tolerance,
            max_matches,
            scans,
        } => query_product_mzml(
            bytes,
            &ProductQuery {
                mz,
                tolerance,
                max_matches,
                scans,
            },
        )
        .map(MassSpecQueryResult::Product),
        MassSpecQuery::Precursor {
            mz,
            tolerance,
            rt_min,
            rt_max,
            charge,
            max_matches,
        } => query_precursor_mzml(
            bytes,
            &PrecursorQuery {
                mz,
                tolerance,
                rt_min,
                rt_max,
                charge,
                max_matches,
            },
        )
        .map(MassSpecQueryResult::Precursor),
        MassSpecQuery::FragmentMatch {
            selector,
            sequence,
            modifications,
            ion_types,
            max_charge,
            tolerance_ppm,
            min_intensity,
            max_matches,
            neutral_losses,
            engine_neutral_losses,
        } => query_fragment_match_mzml(
            bytes,
            &FragmentMatchQuery {
                selector,
                sequence,
                modifications,
                ion_types,
                max_charge,
                tolerance_ppm,
                min_intensity,
                max_matches,
                neutral_losses,
                engine_neutral_losses,
                min_correlation: 0.7,
            },
        )
        .map(MassSpecQueryResult::FragmentMatch),
        MassSpecQuery::CandidateComparison {
            selector,
            sequence,
            candidate_models,
            ion_types,
            max_charge,
            tolerance_ppm,
            min_intensity,
        } => query_candidate_comparison_mzml(
            bytes,
            &CandidateComparisonQuery {
                max_matches: None,
                selector,
                sequence,
                candidate_models,
                neutral_losses: Vec::new(),
                engine_neutral_losses: false,
                ion_types,
                max_charge,
                tolerance_ppm,
                min_intensity,
                min_correlation: 0.7,
            },
        )
        .map(MassSpecQueryResult::CandidateComparison),
        MassSpecQuery::IsotopeEnvelope {
            mass,
            proforma,
            external_mods,
            charge,
            rt_min,
            rt_max,
            tolerance_ppm,
            isotopes,
            max_offset,
        } => query_isotope_envelope_mzml(
            bytes,
            &IsotopeEnvelopeQuery {
                mass,
                proforma,
                external_mods,
                charge,
                rt_min,
                rt_max,
                tolerance_ppm,
                isotopes,
                max_offset,
                theory_locked: false,
            },
        )
        .map(MassSpecQueryResult::IsotopeEnvelope),
        MassSpecQuery::Select { statement } => {
            query_select_mzml(bytes, &statement).map(MassSpecQueryResult::Select)
        }
    }
}

pub fn query_index_mzml(
    bytes: &[u8],
    max_spectra: Option<usize>,
) -> Result<IndexQueryResult, MassSpecError> {
    if max_spectra == Some(0) {
        return Err(MassSpecError::InvalidQuery(
            "maxSpectra must be greater than zero".to_owned(),
        ));
    }
    let mut reader = open_reader(bytes)?;
    let spectra_in_file = reader.len();
    let mut spectra_per_ms_level = BTreeMap::new();
    let mut activation_types = BTreeMap::new();
    let mut retention_time_min: Option<f64> = None;
    let mut retention_time_max: Option<f64> = None;
    let mut native_id_example = None;
    let mut spectra_inspected = 0;

    for spectrum in reader.iter() {
        if max_spectra.is_some_and(|limit| spectra_inspected >= limit) {
            break;
        }
        spectra_inspected += 1;
        native_id_example.get_or_insert_with(|| spectrum.id().to_owned());
        *spectra_per_ms_level
            .entry(spectrum.ms_level().to_string())
            .or_insert(0) += 1;
        if let Some(activation) = spectrum
            .precursor()
            .and_then(|value| value.activation.method())
            .map(|method| method.name().to_owned())
        {
            *activation_types.entry(activation).or_insert(0) += 1;
        }
        let rt = spectrum.start_time();
        if rt.is_finite() {
            retention_time_min = Some(retention_time_min.map_or(rt, |value| value.min(rt)));
            retention_time_max = Some(retention_time_max.map_or(rt, |value| value.max(rt)));
        }
    }

    Ok(IndexQueryResult {
        spectra_in_file,
        spectra_inspected,
        complete: spectra_inspected >= spectra_in_file,
        index_source: "mzdata".to_owned(),
        spectra_per_ms_level,
        activation_types,
        retention_time_min,
        retention_time_max,
        retention_time_unit: "minute".to_owned(),
        native_id_example,
    })
}

pub fn query_scans_mzml(
    bytes: &[u8],
    scans: &[String],
    filters: &ScanQuery,
) -> Result<ScansQueryResult, MassSpecError> {
    if scans.is_empty() {
        return Err(MassSpecError::InvalidQuery(
            "scans must contain at least one scan number".to_owned(),
        ));
    }
    validate_query(filters)?;
    let mut reader = open_reader(bytes)?;
    let mut identities = HashSet::new();
    let mut spectra = Vec::new();
    let mut missing_scans = Vec::new();
    for scan in scans {
        let id = if scan.chars().all(|c| c.is_ascii_digit()) && !scan.is_empty() {
            format!("scan={scan}")
        } else {
            scan.clone()
        };
        let selector =
            crate::fragment::resolve_selector_from_reader(&mut reader, &SpectrumSelector::Id(id))?;
        let query = ScanQuery {
            selector,
            ..filters.clone()
        };
        match query_scan_reader(&mut reader, &query) {
            Ok(spectrum) => {
                if !identities.insert(spectrum.id.clone()) {
                    return Err(MassSpecError::InvalidQuery(
                        "duplicate scan in scan list".into(),
                    ));
                }
                spectra.push(BatchSpectrum {
                    scan: scan.clone(),
                    spectrum,
                });
            }
            Err(MassSpecError::SpectrumNotFound(_)) => missing_scans.push(scan.clone()),
            Err(error) => return Err(error),
        }
    }
    if spectra.is_empty() {
        return Err(MassSpecError::SpectrumNotFound(format!(
            "scans {missing_scans:?}"
        )));
    }
    Ok(ScansQueryResult {
        requested: scans.len(),
        returned: spectra.len(),
        missing_scans,
        spectra,
    })
}

pub fn query_sum_mzml(
    bytes: &[u8],
    scans: &[String],
    bin_ppm: f64,
    mz_min: Option<f64>,
    mz_max: Option<f64>,
    min_intensity: Option<f32>,
    max_peaks: Option<usize>,
) -> Result<SumQueryResult, MassSpecError> {
    if !bin_ppm.is_finite() || bin_ppm < 0.0 {
        return Err(MassSpecError::InvalidQuery(
            "binPpm must be a non-negative finite number".to_owned(),
        ));
    }
    if max_peaks == Some(0) {
        return Err(MassSpecError::InvalidQuery(
            "maxPeaks must be greater than zero".to_owned(),
        ));
    }

    let result = query_scans_mzml(
        bytes,
        scans,
        &ScanQuery {
            selector: SpectrumSelector::Index(0),
            mz_min,
            mz_max,
            min_intensity,
            max_peaks: None,
            peak_order: PeakOrder::MzAscending,
            ms_level: None,
        },
    )?;
    if !result.missing_scans.is_empty() {
        return Err(MassSpecError::SpectrumNotFound(format!(
            "scans {:?}",
            result.missing_scans
        )));
    }
    let levels: HashSet<u8> = result
        .spectra
        .iter()
        .map(|item| item.spectrum.ms_level)
        .collect();
    let flat: Vec<Peak> = result
        .spectra
        .iter()
        .flat_map(|item| item.spectrum.peaks.iter().copied())
        .collect();
    let mut peaks = merge_peaks(flat, bin_ppm);
    let flat_peak_count = peaks.len();
    cap_merged_peaks(&mut peaks, max_peaks);

    let mut warnings = Vec::new();
    if levels.len() > 1 {
        let mut levels: Vec<u8> = levels.into_iter().collect();
        levels.sort_unstable();
        warnings.push(format!("summed spectra span MS levels {levels:?}"));
    }

    let spectra_metadata = result
        .spectra
        .iter()
        .map(|item| crate::SpectrumMetadata {
            id: item.spectrum.id.clone(),
            ms_level: item.spectrum.ms_level,
            retention_time: item.spectrum.retention_time,
            centroided: item.spectrum.centroided,
            precursor: item.spectrum.precursor.clone(),
        })
        .collect();
    let truncated = peaks.len() < flat_peak_count;
    let (scan_numbers, precursors) = result
        .spectra
        .into_iter()
        .map(|item| (item.scan, item.spectrum.precursor))
        .unzip();
    Ok(SumQueryResult {
        spectra_metadata,
        bin_ppm,
        truncated,
        scan_numbers,
        precursors,
        peak_count: peaks.len(),
        peaks,
        warnings,
    })
}

#[allow(clippy::too_many_arguments)]
pub fn query_window_mzml(
    bytes: &[u8],
    ms_level: u8,
    rt_min: f64,
    rt_max: f64,
    max_scans: usize,
    bin_ppm: f64,
    mz_min: Option<f64>,
    mz_max: Option<f64>,
    min_intensity: Option<f32>,
    max_peaks: Option<usize>,
) -> Result<WindowQueryResult, MassSpecError> {
    if ms_level == 0 || !rt_min.is_finite() || !rt_max.is_finite() || rt_min > rt_max {
        return Err(MassSpecError::InvalidQuery(
            "window requires a positive msLevel and ordered finite retention-time range".to_owned(),
        ));
    }
    if max_scans == 0 || !bin_ppm.is_finite() || bin_ppm < 0.0 || max_peaks == Some(0) {
        return Err(MassSpecError::InvalidQuery(
            "window limits must be positive and binPpm must be non-negative".to_owned(),
        ));
    }
    let mut reader = open_reader(bytes)?;
    let mut scans = select_spectrum_ids_reader(
        &mut reader,
        &SpectrumSelection {
            ms_level: Some(ms_level),
            rt_min: Some(rt_min),
            rt_max: Some(rt_max),
            ..Default::default()
        },
    )?;
    let scan_truncated = scans.len() > max_scans;
    scans.truncate(max_scans);
    let sum = query_sum_mzml(
        bytes,
        &scans,
        bin_ppm,
        mz_min,
        mz_max,
        min_intensity,
        max_peaks,
    )?;
    Ok(WindowQueryResult {
        spectra_metadata: sum.spectra_metadata,
        bin_ppm,
        scans_summed: sum.scan_numbers.len(),
        scan_numbers: sum.scan_numbers,
        retention_time_unit: "minute".into(),
        truncated: scan_truncated || sum.truncated,
        peak_count: sum.peak_count,
        peaks: sum.peaks,
    })
}

pub fn query_xic_mzml(
    bytes: &[u8],
    mz: f64,
    tolerance_ppm: f64,
    ms_level: u8,
    rt_min: Option<f64>,
    rt_max: Option<f64>,
) -> Result<XicQueryResult, MassSpecError> {
    if !mz.is_finite()
        || mz <= 0.0
        || !tolerance_ppm.is_finite()
        || tolerance_ppm <= 0.0
        || ms_level == 0
    {
        return Err(MassSpecError::InvalidQuery(
            "xic requires positive mz, tolerancePpm, and msLevel".to_owned(),
        ));
    }
    if rt_min.zip(rt_max).is_some_and(|(min, max)| min > max) {
        return Err(MassSpecError::InvalidQuery(
            "xic retention-time range must be ordered".to_owned(),
        ));
    }
    let window = mz * tolerance_ppm / 1_000_000.0;
    let mut reader = open_reader(bytes)?;
    let mut points = Vec::new();
    for mut spectrum in reader.iter() {
        if spectrum.ms_level() != ms_level {
            continue;
        }
        let rt = spectrum.start_time();
        if !rt.is_finite()
            || rt_min.is_some_and(|min| rt < min)
            || rt_max.is_some_and(|max| rt > max)
        {
            continue;
        }
        spectrum
            .try_build_peaks()
            .map_err(|error| MassSpecError::PeakDecode(error.to_string()))?;
        let intensity = spectrum
            .peaks()
            .iter()
            .filter(|peak| (peak.mz - mz).abs() <= window)
            .map(|peak| peak.intensity)
            .sum();
        points.push(XicPoint {
            retention_time: rt,
            scan: scan_token(spectrum.id()).unwrap_or_else(|| spectrum.id().to_owned()),
            intensity,
        });
    }
    let apex = points
        .iter()
        .max_by(|left, right| left.intensity.total_cmp(&right.intensity))
        .cloned();
    Ok(XicQueryResult {
        point_count: points.len(),
        retention_time_unit: if points.is_empty() {
            String::new()
        } else {
            "minute".to_owned()
        },
        apex,
        points,
    })
}

pub(crate) fn merge_peaks(mut peaks: Vec<Peak>, bin_ppm: f64) -> Vec<MergedPeak> {
    peaks.sort_unstable_by(|left, right| left.mz.total_cmp(&right.mz));
    let mut merged: Vec<MergedPeak> = Vec::new();
    let mut cluster_start = 0.0;
    for peak in peaks {
        if let Some(last) = merged.last_mut() {
            // Bound the entire cluster span; a drifting centroid permits chaining.
            if peak.mz - cluster_start <= cluster_start * bin_ppm / 1_000_000.0 {
                let total = last.intensity + peak.intensity;
                if total != 0.0 {
                    last.mz = (last.mz * f64::from(last.intensity)
                        + peak.mz * f64::from(peak.intensity))
                        / f64::from(total);
                }
                last.intensity = total;
                last.merged_from += 1;
                continue;
            }
        }
        cluster_start = peak.mz;
        merged.push(MergedPeak {
            mz: peak.mz,
            intensity: peak.intensity,
            merged_from: 1,
        });
    }
    merged
}

#[cfg(test)]
mod aggregation_tests {
    use super::*;

    #[test]
    fn merging_tolerates_drift_without_chaining_distinct_peaks() {
        let peaks = vec![
            Peak {
                mz: 1000.000,
                intensity: 100.0,
            },
            Peak {
                mz: 1000.009,
                intensity: 100.0,
            },
            Peak {
                mz: 1000.014,
                intensity: 100.0,
            },
        ];
        let merged = merge_peaks(peaks.clone(), 10.0);
        assert_eq!(
            merged.len(),
            2,
            "the full cluster span must stay within 10 ppm"
        );
        assert!((merged[0].mz - 1000.0045).abs() < 1e-9);
        assert_eq!(merged[0].intensity, 200.0);
        assert_eq!(merged.iter().map(|p| p.intensity).sum::<f32>(), 300.0);
        let reverse = merge_peaks(peaks.into_iter().rev().collect(), 10.0);
        assert_eq!(merged, reverse, "scan ordering cannot change the sum");
    }
}

fn cap_merged_peaks(peaks: &mut Vec<MergedPeak>, max_peaks: Option<usize>) {
    if let Some(limit) = max_peaks {
        peaks.sort_unstable_by(|left, right| {
            right
                .intensity
                .total_cmp(&left.intensity)
                .then_with(|| left.mz.total_cmp(&right.mz))
        });
        peaks.truncate(limit);
        peaks.sort_unstable_by(|left, right| left.mz.total_cmp(&right.mz));
    }
}

fn order_and_cap(peaks: &mut Vec<Peak>, order: PeakOrder, max_peaks: Option<usize>) {
    match order {
        PeakOrder::MzAscending => {
            peaks.sort_unstable_by(|left, right| left.mz.total_cmp(&right.mz));
        }
        PeakOrder::IntensityDescending => peaks.sort_unstable_by(|left, right| {
            right
                .intensity
                .total_cmp(&left.intensity)
                .then_with(|| left.mz.total_cmp(&right.mz))
        }),
        PeakOrder::IntensityCappedMzAscending => {
            peaks.sort_unstable_by(|left, right| {
                right
                    .intensity
                    .total_cmp(&left.intensity)
                    .then_with(|| left.mz.total_cmp(&right.mz))
            });
            if let Some(limit) = max_peaks {
                peaks.truncate(limit);
            }
            peaks.sort_unstable_by(|left, right| left.mz.total_cmp(&right.mz));
            return;
        }
    }
    if let Some(limit) = max_peaks {
        peaks.truncate(limit);
    }
}

pub fn query_mzml(bytes: &[u8], query: &ScanQuery) -> Result<ScanQueryResult, MassSpecError> {
    let mut reader = open_reader(bytes)?;
    let selector = crate::fragment::resolve_selector_from_reader(&mut reader, &query.selector)?;
    query_scan_reader(
        &mut reader,
        &ScanQuery {
            selector,
            ..query.clone()
        },
    )
}

pub(crate) fn query_scan_reader(
    reader: &mut crate::Reader,
    query: &ScanQuery,
) -> Result<ScanQueryResult, MassSpecError> {
    validate_query(query)?;
    let selector_label = match &query.selector {
        SpectrumSelector::Index(index) => format!("index {index}"),
        SpectrumSelector::Id(id) => format!("id {id}"),
    };
    let mut spectrum = match &query.selector {
        SpectrumSelector::Index(index) => reader.get_spectrum_by_index(*index),
        SpectrumSelector::Id(id) => reader.get_spectrum_by_id(id),
    }
    .ok_or(MassSpecError::SpectrumNotFound(selector_label))?;
    if query
        .ms_level
        .is_some_and(|level| spectrum.ms_level() != level)
    {
        return Err(MassSpecError::SpectrumNotFound(spectrum.id().to_owned()));
    }

    spectrum
        .try_build_peaks()
        .map_err(|error| MassSpecError::PeakDecode(error.to_string()))?;

    let total_peak_count = spectrum.peaks().len();
    let mut peaks: Vec<Peak> = spectrum
        .peaks()
        .iter()
        .filter(|peak| query.mz_min.is_none_or(|minimum| peak.mz >= minimum))
        .filter(|peak| query.mz_max.is_none_or(|maximum| peak.mz <= maximum))
        .filter(|peak| {
            query
                .min_intensity
                .is_none_or(|minimum| peak.intensity >= minimum)
        })
        .map(|peak| Peak {
            mz: peak.mz,
            intensity: peak.intensity,
        })
        .collect();
    let matched_peak_count = peaks.len();

    order_and_cap(&mut peaks, query.peak_order, query.max_peaks);

    let precursor = spectrum.precursor().and_then(|value| {
        let ion = value.ions.first()?;
        Some(ScanPrecursorMetadata {
            selected_mz: ion.mz,
            intensity: ion.intensity,
            charge: ion.charge,
            isolation_target: value.isolation_window.target,
            isolation_lower_bound: value.isolation_window.lower_bound,
            isolation_upper_bound: value.isolation_window.upper_bound,
            activation: value
                .activation
                .method()
                .map(|method| method.name().to_owned()),
            collision_energy: (value.activation.energy != 0.0).then_some(value.activation.energy),
        })
    });

    Ok(ScanQueryResult {
        id: spectrum.id().to_owned(),
        index: spectrum.index(),
        ms_level: spectrum.ms_level(),
        centroided: spectrum.signal_continuity() == mzdata::spectrum::SignalContinuity::Centroid,
        retention_time: spectrum.start_time(),
        retention_time_unit: "minute".to_owned(),
        precursor,
        total_peak_count,
        matched_peak_count,
        truncated: peaks.len() < matched_peak_count,
        peaks,
    })
}

pub fn query_product_mzml(
    bytes: &[u8],
    query: &ProductQuery,
) -> Result<ProductQueryResult, MassSpecError> {
    if !query.mz.is_finite() || query.mz <= 0.0 {
        return Err(MassSpecError::InvalidQuery(
            "product mz must be positive".to_owned(),
        ));
    }
    let limit = query.max_matches.unwrap_or(100);
    if limit == 0 {
        return Err(MassSpecError::InvalidQuery(
            "maxMatches must be greater than zero".to_owned(),
        ));
    }
    let mut reader = open_reader(bytes)?;
    let mut scanned_spectrum_count = 0;
    let mut matches = Vec::new();
    let mut truncated = false;

    for mut spectrum in reader.iter() {
        if spectrum.ms_level() != 2 {
            continue;
        }
        // A scan restriction turns the whole-file census into evidence about
        // the named spectra only; unfiltered stays the file-wide survey.
        if let Some(wanted) = &query.scans {
            let token = scan_token(spectrum.id());
            let id = spectrum.id();
            if !wanted.iter().any(|selector| {
                selector == id
                    || token.as_deref().is_some_and(|token| {
                        selector.strip_prefix("scan=").unwrap_or(selector) == token
                    })
            }) {
                continue;
            }
        }
        scanned_spectrum_count += 1;
        spectrum
            .try_build_peaks()
            .map_err(|error| MassSpecError::PeakDecode(error.to_string()))?;
        let matched = spectrum
            .peaks()
            .iter()
            .filter(|peak| query.tolerance.contains(query.mz, peak.mz))
            .max_by(|left, right| left.intensity.total_cmp(&right.intensity));
        if let Some(peak) = matched {
            if matches.len() == limit {
                truncated = true;
                break;
            }
            matches.push(ProductMatch {
                id: spectrum.id().to_owned(),
                index: spectrum.index(),
                ms_level: spectrum.ms_level(),
                mz: peak.mz,
                intensity: peak.intensity,
            });
        }
    }

    Ok(ProductQueryResult {
        query: ProductQuerySummary {
            mz: query.mz,
            tolerance: query.tolerance,
            scans: query.scans.clone(),
        },
        scanned_spectrum_count,
        matches,
        truncated,
    })
}

pub fn query_precursor_mzml(
    bytes: &[u8],
    query: &PrecursorQuery,
) -> Result<PrecursorQueryResult, MassSpecError> {
    if !query.mz.is_finite() || query.mz <= 0.0 {
        return Err(MassSpecError::InvalidQuery(
            "precursor mz must be positive".to_owned(),
        ));
    }
    let limit = query.max_matches.unwrap_or(100);
    if limit == 0 {
        return Err(MassSpecError::InvalidQuery(
            "maxMatches must be greater than zero".to_owned(),
        ));
    }
    let mut reader = open_reader(bytes)?;
    let mut scanned_spectrum_count = 0;
    let mut matches = Vec::new();
    let mut truncated = false;

    if query
        .rt_min
        .zip(query.rt_max)
        .is_some_and(|(minimum, maximum)| minimum > maximum)
    {
        return Err(MassSpecError::InvalidQuery(
            "precursor retention-time range must be ordered".to_owned(),
        ));
    }

    for mut spectrum in reader.iter() {
        if spectrum.ms_level() != 2 {
            continue;
        }
        scanned_spectrum_count += 1;
        let retention_time = spectrum.start_time();
        if query.rt_min.is_some_and(|minimum| retention_time < minimum)
            || query.rt_max.is_some_and(|maximum| retention_time > maximum)
        {
            continue;
        }
        let precursor = spectrum.precursor().and_then(|value| {
            let ion = value.ions.first()?;
            if !query.tolerance.contains(query.mz, ion.mz)
                || query
                    .charge
                    .is_some_and(|charge| ion.charge.is_some_and(|value| value != charge))
            {
                return None;
            }
            let activation = value
                .activation
                .method()
                .map(|method| method.name().to_owned());
            let fragment_ion_types =
                crate::fragment::ion_types_for_activation(activation.as_deref()).ok();
            Some(PrecursorMetadata {
                selected_mz: ion.mz,
                intensity: ion.intensity,
                charge: ion.charge,
                isolation_target: value.isolation_window.target,
                isolation_lower_bound: value.isolation_window.lower_bound,
                isolation_upper_bound: value.isolation_window.upper_bound,
                activation,
                fragment_ion_types,
                collision_energy: (value.activation.energy != 0.0)
                    .then_some(value.activation.energy),
            })
        });
        if precursor.is_none() {
            continue;
        }
        if matches.len() == limit {
            truncated = true;
            break;
        }
        let header_tic = header_total_ion_current(spectrum.params());
        spectrum
            .try_build_peaks()
            .map_err(|error| MassSpecError::PeakDecode(error.to_string()))?;
        matches.push(SelectedSpectrum {
            id: spectrum.id().to_owned(),
            scan: scan_token(spectrum.id()).unwrap_or_else(|| spectrum.id().to_owned()),
            index: spectrum.index(),
            ms_level: spectrum.ms_level(),
            retention_time,
            retention_time_unit: "minute".to_owned(),
            total_ion_current: header_tic.unwrap_or_else(|| {
                spectrum
                    .peaks()
                    .iter()
                    .map(|peak| peak.intensity as f64)
                    .sum()
            }),
            peak_count: Some(spectrum.peaks().len()),
            centroided: spectrum.signal_continuity()
                == mzdata::spectrum::SignalContinuity::Centroid,
            precursor,
            matched_products: Vec::new(),
        });
    }

    Ok(PrecursorQueryResult {
        query: PrecursorQuerySummary {
            mz: query.mz,
            tolerance: query.tolerance,
            rt_min: query.rt_min,
            rt_max: query.rt_max,
            charge: query.charge,
        },
        scanned_spectrum_count,
        matches,
        truncated,
    })
}

fn validate_query(query: &ScanQuery) -> Result<(), MassSpecError> {
    if query
        .mz_min
        .zip(query.mz_max)
        .is_some_and(|(minimum, maximum)| minimum > maximum)
    {
        return Err(MassSpecError::InvalidQuery(
            "mzMin must be less than or equal to mzMax".to_owned(),
        ));
    }
    if query.max_peaks == Some(0) {
        return Err(MassSpecError::InvalidQuery(
            "maxPeaks must be greater than zero".to_owned(),
        ));
    }
    Ok(())
}

// ── MS1 isotope-based precursor validation (IsoForma / Degnan 2024) ─────────

/// Validate an MS2 scan's precursor against its nearest MS1 isotope envelope.
///
/// Faithful to the IsoForma scan-detection step: for a given MS2 scan, the
/// precursor's theoretical isotopic distribution is matched to the MS1 spectra
/// nearest in retention time; if at least `min_matched_peaks` (typically 3)
/// isotope peaks match, the precursor is accepted. Uses the same
/// `query_isotope_envelope_mzml` kernel as the `envelope --model` path so the
/// theoretical shape and the observed aggregation come from one source.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ValidatePrecursorQuery {
    pub proforma: String,
    pub external_mods: Vec<crate::proforma::ExternalModification>,
    /// Charge used to build the theoretical envelope and aggregate the MS1.
    /// When `None`, the MS2 scan's own precursor charge is read from the
    /// spectrum header (IsoForma uses the acquisition-assigned charge; a
    /// caller can override when that charge is known to be wrong).
    pub charge: Option<i32>,
    /// MS2 scan token (e.g. "4337" or "scan=4337").
    pub scan: String,
    pub tolerance_ppm: f64,
    pub min_matched_peaks: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ValidatePrecursorResult {
    pub ms2_scan: String,
    pub ms2_rt_min: f64,
    pub ms1_scan_count: usize,
    pub charge: i32,
    pub matched_peaks: usize,
    pub matched_peaks_in_isolation_window: usize,
    pub theoretical_peak_count: usize,
    pub minimum_profile_matches: usize,
    pub min_matched_peaks: usize,
    pub isolation_window: [f64; 2],
    pub is_valid: bool,
    pub pearson: Option<f64>,
    pub refined_mass: f64,
    pub envelope_evaluations: Vec<crate::envelope::EnvelopeEvaluation>,
}

pub fn validate_precursor_mzml(
    bytes: &[u8],
    query: &ValidatePrecursorQuery,
) -> Result<ValidatePrecursorResult, MassSpecError> {
    validate_precursors_mzml(bytes, query, std::slice::from_ref(&query.scan))?
        .into_iter()
        .next()
        .ok_or_else(|| MassSpecError::InvalidQuery("no precursor scan requested".to_owned()))
}

/// Batch form of [`validate_precursor_mzml`]. The mzML and its parent MS1
/// spectra are traversed once for all requested MS2 scans.
pub fn validate_precursors_mzml(
    bytes: &[u8],
    query: &ValidatePrecursorQuery,
    scans: &[String],
) -> Result<Vec<ValidatePrecursorResult>, MassSpecError> {
    if scans.is_empty() {
        return Err(MassSpecError::InvalidQuery(
            "at least one MS2 scan is required".to_owned(),
        ));
    }
    let targets: HashSet<String> = scans
        .iter()
        .map(|scan| scan.strip_prefix("scan=").unwrap_or(scan).to_owned())
        .collect();
    #[derive(Clone)]
    struct GateMetadata {
        rt: f64,
        charge: Option<i32>,
        precursor_id: Option<String>,
        isolation_window: [f64; 2],
    }
    let mut reader = open_reader(bytes)?;
    let mut ms1 = Vec::new();
    let mut metadata = BTreeMap::new();
    for spectrum in reader.iter() {
        if spectrum.ms_level() == 1 {
            ms1.push((spectrum.id().to_owned(), spectrum.start_time()));
            continue;
        }
        if spectrum.ms_level() != 2 {
            continue;
        }
        let Some(token) = scan_token(spectrum.id()) else {
            continue;
        };
        if targets.contains(&token) {
            let mut gate = GateMetadata {
                rt: spectrum.start_time(),
                charge: None,
                precursor_id: None,
                isolation_window: [f64::NAN, f64::NAN],
            };
            if let Some(precursor) = spectrum.precursor() {
                gate.charge = precursor.ions.first().and_then(|ion| ion.charge);
                gate.precursor_id = precursor.precursor_id.clone();
                gate.isolation_window = [
                    f64::from(precursor.isolation_window.lower_bound),
                    f64::from(precursor.isolation_window.upper_bound),
                ];
            }
            metadata.insert(token, gate);
        }
    }
    let mut parent_reader = open_reader(bytes)?;
    let mut output = Vec::with_capacity(scans.len());
    for requested_scan in scans {
        let token = requested_scan
            .strip_prefix("scan=")
            .unwrap_or(requested_scan);
        let gate = metadata
            .get(token)
            .ok_or_else(|| MassSpecError::InvalidQuery(format!("MS2 scan {token} not found")))?;
        let charge = query.charge.or(gate.charge).ok_or_else(|| {
            MassSpecError::InvalidQuery(format!(
                "no precursor charge for scan {token}: pass --charge or rely on the MS2 header"
            ))
        })?;
        let charge_u8 = u8::try_from(charge).map_err(|_| {
            MassSpecError::InvalidQuery(format!("invalid precursor charge {charge}"))
        })?;
        let parent_id = gate.precursor_id.clone().or_else(|| {
            ms1.iter()
                .min_by(|left, right| {
                    (left.1 - gate.rt)
                        .abs()
                        .total_cmp(&(right.1 - gate.rt).abs())
                })
                .map(|(id, _)| id.clone())
        });
        let mut parent = parent_id
            .as_deref()
            .and_then(|id| parent_reader.get_spectrum_by_id(id))
            .ok_or_else(|| {
                MassSpecError::InvalidQuery("no parent MS1 spectrum found".to_owned())
            })?;
        parent
            .try_build_peaks()
            .map_err(|error| MassSpecError::PeakDecode(error.to_string()))?;
        let peaks: Vec<crate::Peak> = parent
            .peaks()
            .iter()
            .map(|peak| crate::Peak {
                mz: peak.mz,
                intensity: peak.intensity,
            })
            .collect();
        let envelope_query = IsotopeEnvelopeQuery {
            mass: None,
            proforma: Some(query.proforma.clone()),
            external_mods: query.external_mods.clone(),
            charge: charge_u8,
            rt_min: gate.rt,
            rt_max: gate.rt,
            tolerance_ppm: query.tolerance_ppm,
            isotopes: 0,
            max_offset: 0,
            theory_locked: false,
        };
        let best = evaluate_isotope_envelope_peaks(&peaks, &envelope_query)?;
        let theoretical_peak_count = best.isotopes_considered;
        let minimum_profile_matches = (theoretical_peak_count as f64 * 0.4).floor() as usize;
        let matched_peaks_in_isolation_window = best
            .matched_isotopes
            .iter()
            .filter(|peak| {
                peak.theoretical_mz >= gate.isolation_window[0]
                    && peak.theoretical_mz <= gate.isolation_window[1]
            })
            .count();
        output.push(ValidatePrecursorResult {
            ms2_scan: requested_scan.clone(),
            ms2_rt_min: gate.rt,
            ms1_scan_count: 1,
            charge,
            matched_peaks: best.isotopes_matched,
            matched_peaks_in_isolation_window,
            theoretical_peak_count,
            minimum_profile_matches,
            min_matched_peaks: query.min_matched_peaks,
            isolation_window: gate.isolation_window,
            is_valid: matched_peaks_in_isolation_window >= query.min_matched_peaks
                && best.isotopes_matched >= minimum_profile_matches,
            pearson: best.pearson,
            refined_mass: best.monoisotopic_mass,
            envelope_evaluations: vec![best],
        });
    }
    Ok(output)
}
