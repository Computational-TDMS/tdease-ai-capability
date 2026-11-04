//! Reusable centroid observations. Selection and pooling are independent of
//! fragment theory; every matcher sees the same ordered, indexed peak list.

use std::collections::HashSet;

use serde::{Deserialize, Serialize};

use crate::fragment::resolve_selector_from_reader;
use crate::{MassSpecError, Peak, Reader, ScanPrecursorMetadata, SpectrumSelector};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SpectrumMetadata {
    pub id: String,
    pub ms_level: u8,
    pub retention_time: f64,
    pub centroided: bool,
    pub precursor: Option<ScanPrecursorMetadata>,
}

/// Also reads the result of `spectrum --mode sum|window`: extra projection
/// fields (peakCount, mergedFrom, etc.) do not change the observation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SpectrumObservation {
    pub spectra_metadata: Vec<SpectrumMetadata>,
    pub peaks: Vec<Peak>,
    #[serde(default)]
    pub truncated: bool,
    pub bin_ppm: Option<f64>,
    #[serde(default)]
    pub minimum_abundance_pct: f64,
}

impl SpectrumObservation {
    pub(crate) fn selected_peaks(&self, minimum: Option<f32>) -> std::borrow::Cow<'_, [Peak]> {
        if minimum.is_none_or(|min| self.peaks.iter().all(|peak| peak.intensity >= min)) {
            std::borrow::Cow::Borrowed(&self.peaks)
        } else {
            std::borrow::Cow::Owned(
                self.peaks
                    .iter()
                    .copied()
                    .filter(|peak| minimum.is_none_or(|min| peak.intensity >= min))
                    .collect(),
            )
        }
    }
    pub fn validate(&self) -> Result<(), MassSpecError> {
        if self.spectra_metadata.is_empty() || self.truncated {
            return Err(MassSpecError::InvalidQuery(
                "matching requires a complete observation with source spectrum metadata".into(),
            ));
        }
        let mut ids = HashSet::new();
        for metadata in &self.spectra_metadata {
            if metadata.id.is_empty() || !ids.insert(&metadata.id) {
                return Err(MassSpecError::InvalidQuery(
                    "observation contains an empty or duplicate spectrum id".into(),
                ));
            }
            if !metadata.centroided || metadata.ms_level == 0 {
                return Err(MassSpecError::InvalidQuery(
                    "observation requires centroided spectra with a declared MS level".into(),
                ));
            }
        }
        validate_peaks(&self.peaks)?;
        if (self.spectra_metadata.len() > 1 && self.bin_ppm.is_none())
            || self
                .bin_ppm
                .is_some_and(|value| !value.is_finite() || value < 0.0)
            || !self.minimum_abundance_pct.is_finite()
            || !(0.0..=100.0).contains(&self.minimum_abundance_pct)
        {
            return Err(MassSpecError::InvalidQuery(
                "invalid observation aggregation parameters".into(),
            ));
        }
        Ok(())
    }

    /// Pooling arithmetic accepts arbitrary scans. A fragment measurement
    /// additionally needs a coherent activation method. An explicit method
    /// can fill absent metadata, but cannot silently override a known conflict.
    pub fn fragment_activation(
        &self,
        explicit: Option<&str>,
    ) -> Result<Option<String>, MassSpecError> {
        self.validate()?;
        if self.spectra_metadata.iter().any(|m| m.ms_level != 2) {
            return Err(MassSpecError::InvalidQuery(
                "fragment matching requires MS2 observations".into(),
            ));
        }
        let mut known: Option<String> = None;
        let mut missing = false;
        for metadata in &self.spectra_metadata {
            let method = metadata
                .precursor
                .as_ref()
                .and_then(|p| p.activation.as_ref());
            match method {
                Some(method) => {
                    if known.as_ref().is_some_and(|previous| previous != method) {
                        return Err(MassSpecError::InvalidQuery(
                            "fragment pools cannot mix activation methods".into(),
                        ));
                    }
                    known = Some(method.clone());
                }
                None => missing = true,
            }
        }
        if self.spectra_metadata.len() > 1 && missing && explicit.is_none() {
            return Err(MassSpecError::InvalidQuery(
                "fragment pool has unknown activation; supply --activation explicitly".into(),
            ));
        }
        if let (Some(method), Some(declared)) = (&known, explicit) {
            if crate::ion_types_for_activation(Some(method))?
                != crate::ion_types_for_activation(Some(declared))?
            {
                return Err(MassSpecError::InvalidQuery(
                    "explicit activation conflicts with pool metadata".into(),
                ));
            }
        }
        Ok(explicit.map(str::to_owned).or(known))
    }
}

pub(crate) fn validate_peaks(peaks: &[Peak]) -> Result<(), MassSpecError> {
    if peaks
        .iter()
        .any(|p| !p.mz.is_finite() || p.mz <= 0.0 || !p.intensity.is_finite() || p.intensity <= 0.0)
        || peaks.windows(2).any(|pair| pair[0].mz > pair[1].mz)
    {
        return Err(MassSpecError::InvalidQuery(
            "centroid peaks must have positive finite coordinates/intensities and be sorted by m/z"
                .into(),
        ));
    }
    Ok(())
}

/// Extract once, optionally pool, then reuse for any number of model queries.
pub fn extract_spectrum_observation_reader(
    reader: &mut Reader,
    selectors: &[SpectrumSelector],
    bin_ppm: Option<f64>,
    minimum_abundance_pct: f64,
) -> Result<SpectrumObservation, MassSpecError> {
    if selectors.is_empty()
        || bin_ppm.is_some_and(|v| !v.is_finite() || v < 0.0)
        || !minimum_abundance_pct.is_finite()
        || !(0.0..=100.0).contains(&minimum_abundance_pct)
        || (selectors.len() > 1 && bin_ppm.is_none())
    {
        return Err(MassSpecError::InvalidQuery(
            "invalid spectrum pool selection or aggregation parameters".into(),
        ));
    }
    let mut seen = HashSet::new();
    let mut spectra_metadata = Vec::with_capacity(selectors.len());
    let mut peaks = Vec::new();
    for requested in selectors {
        let selector = resolve_selector_from_reader(reader, requested)?;
        let spectrum = crate::query::query_scan_reader(
            reader,
            &crate::ScanQuery {
                selector,
                mz_min: None,
                mz_max: None,
                min_intensity: None,
                max_peaks: None,
                peak_order: crate::PeakOrder::MzAscending,
                ms_level: None,
            },
        )?;
        if !seen.insert(spectrum.id.clone()) {
            return Err(MassSpecError::InvalidQuery(
                "duplicate scan in spectrum pool".into(),
            ));
        }
        if !spectrum.centroided {
            return Err(MassSpecError::InvalidQuery(
                "spectrum pooling requires centroided spectra".into(),
            ));
        }
        spectra_metadata.push(SpectrumMetadata {
            id: spectrum.id,
            ms_level: spectrum.ms_level,
            retention_time: spectrum.retention_time,
            centroided: spectrum.centroided,
            precursor: spectrum.precursor,
        });
        peaks.extend(spectrum.peaks.into_iter().filter(|p| {
            p.mz.is_finite() && p.mz > 0.0 && p.intensity.is_finite() && p.intensity > 0.0
        }));
    }
    let minimum = peaks.iter().map(|p| p.intensity).fold(0.0_f32, f32::max)
        * (minimum_abundance_pct as f32 / 100.0);
    peaks.retain(|p| p.intensity > minimum);
    let peaks = if let Some(bin) = bin_ppm {
        crate::query::merge_peaks(peaks, bin)
            .into_iter()
            .map(|p| Peak {
                mz: p.mz,
                intensity: p.intensity,
            })
            .collect()
    } else {
        peaks.sort_unstable_by(|a, b| a.mz.total_cmp(&b.mz));
        peaks
    };
    validate_peaks(&peaks)?;
    Ok(SpectrumObservation {
        spectra_metadata,
        peaks,
        truncated: false,
        bin_ppm,
        minimum_abundance_pct,
    })
}
