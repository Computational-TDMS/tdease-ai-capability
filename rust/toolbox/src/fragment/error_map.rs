use serde::{Deserialize, Serialize};

use crate::constants::{NEUTRON_SHIFT, PROTON};
use crate::Peak;

use super::{FragmentMatch, FragmentQuantification, IonType, TheoreticalFragment};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum FragmentMatchState {
    Matched,
    MissingEnvelope,
    UnreasonableCharge,
    MissingApex,
    CorrelationRejected,
    /// Deconvoluted (msalign) mode: no peak of this charge exists at this
    /// neutral mass. Raw-mzML statuses never use it.
    MissingPeak,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct FragmentMatchStatus {
    pub ion: IonType,
    pub number: usize,
    pub charge: u8,
    pub cleavage_after: usize,
    pub theoretical_mz: f64,
    pub theoretical_apex_mz: f64,
    pub apex_offset: usize,
    pub state: FragmentMatchState,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub envelope_correlation: Option<f64>,
}

/// Raw-MS2 error-map matching: chemical isotope envelopes (intact-shape
/// approximations for numeric named losses) combined with
/// Informed-Proteomics' `GetAllIsotopePeaks` traversal and LCMS-Spectator's
/// `LabeledIonViewModel.GetPeakDataPoints` gates:
/// 1. the most-abundant theoretical isotope must find a peak;
/// 2. neighboring isotope positions are visited with a monotone peak cursor;
/// 3. missing positions enter the full-envelope Pearson as zero;
/// 4. charge must be at most `round(max(4, neutral_mass / 1000 * 3))`;
/// 5. the error-map cell uses the most intense observed isotope.
pub(super) fn match_fragments_envelope_with_status(
    fragments: &[TheoreticalFragment],
    peaks: &[Peak],
    tolerance_ppm: f64,
    min_correlation: f64,
) -> (Vec<FragmentMatch>, Vec<FragmentMatchStatus>) {
    const REL_INT_THRESHOLD: f64 = 0.1;
    let find_index = |target: f64| -> Option<usize> {
        let window = target * tolerance_ppm / 1_000_000.0;
        let lo = peaks.partition_point(|peak| peak.mz <= target - window);
        let hi = peaks.partition_point(|peak| peak.mz < target + window);
        (lo..hi).max_by(|left, right| peaks[*left].intensity.total_cmp(&peaks[*right].intensity))
    };
    let mut matches = Vec::new();
    let mut statuses = Vec::with_capacity(fragments.len());
    for fragment in fragments {
        // Spectator's matching coordinate is nominal C13 spacing, not the
        // BRAIN centroid coordinate exported by general fragment theory.
        let isotope_mz = |offset: usize| {
            (fragment.neutral_mass + offset as f64 * NEUTRON_SHIFT) / f64::from(fragment.charge)
                + PROTON
        };
        let status = |state, envelope_correlation| FragmentMatchStatus {
            ion: fragment.ion,
            number: fragment.number,
            charge: fragment.charge,
            cleavage_after: fragment.cleavage_after,
            theoretical_mz: fragment.mz,
            theoretical_apex_mz: isotope_mz(fragment.apex_offset),
            apex_offset: fragment.apex_offset,
            state,
            envelope_correlation,
        };
        if fragment.isotope_envelope.is_empty() {
            statuses.push(status(FragmentMatchState::MissingEnvelope, None));
            continue;
        }
        let max_reasonable_charge = 4.0f64.max(fragment.neutral_mass / 1000.0 * 3.0).round() as u8;
        if fragment.charge > max_reasonable_charge {
            statuses.push(status(FragmentMatchState::UnreasonableCharge, None));
            continue;
        }
        let envelope = &fragment.isotope_envelope;
        let apex_index = match envelope.iter().position(|position| {
            position.neutron_offset == fragment.apex_offset
                && position.relative_intensity >= REL_INT_THRESHOLD
        }) {
            Some(index) => index,
            _ => {
                statuses.push(status(FragmentMatchState::MissingEnvelope, None));
                continue;
            }
        };
        let apex = &envelope[apex_index];
        let apex_peak_index = match find_index(isotope_mz(apex.neutron_offset)) {
            Some(index) => index,
            None => {
                statuses.push(status(FragmentMatchState::MissingApex, None));
                continue;
            }
        };

        let mut observed: Vec<Option<usize>> = vec![None; envelope.len()];
        observed[apex_index] = Some(apex_peak_index);
        let mut peak_index = apex_peak_index as isize - 1;
        for offset in (0..apex_index).rev() {
            let position = match envelope.get(offset) {
                Some(position) if position.relative_intensity >= REL_INT_THRESHOLD => position,
                _ => break,
            };
            let mz = isotope_mz(position.neutron_offset);
            let window = mz * tolerance_ppm / 1_000_000.0;
            let (min_mz, max_mz) = (mz - window, mz + window);
            let mut i = peak_index;
            while i >= 0 {
                let index = i as usize;
                let mz = peaks[index].mz;
                if mz < min_mz {
                    peak_index = i;
                    break;
                }
                if mz <= max_mz {
                    let slot = &mut observed[offset];
                    if slot.is_none_or(|current| peaks[index].intensity > peaks[current].intensity)
                    {
                        *slot = Some(index);
                    }
                }
                i -= 1;
            }
            if i < 0 {
                peak_index = -1;
            }
        }
        let mut peak_index = apex_peak_index + 1;
        for (offset, slot) in observed.iter_mut().enumerate().skip(apex_index + 1) {
            let position = match envelope.get(offset) {
                Some(position) if position.relative_intensity >= REL_INT_THRESHOLD => position,
                _ => break,
            };
            let mz = isotope_mz(position.neutron_offset);
            let window = mz * tolerance_ppm / 1_000_000.0;
            let (min_mz, max_mz) = (mz - window, mz + window);
            let mut i = peak_index;
            while i < peaks.len() {
                let mz = peaks[i].mz;
                if mz > max_mz {
                    peak_index = i;
                    break;
                }
                if mz >= min_mz {
                    if slot.is_none_or(|current| peaks[i].intensity > peaks[current].intensity) {
                        *slot = Some(i);
                    }
                }
                i += 1;
            }
            if i >= peaks.len() {
                peak_index = peaks.len();
            }
        }

        let observed_intensities: Vec<f64> = observed
            .iter()
            .map(|slot| {
                slot.map(|index| f64::from(peaks[index].intensity))
                    .unwrap_or(0.0)
            })
            .collect();
        let theoretical_intensities: Vec<f64> = envelope
            .iter()
            .map(|peak| peak.relative_intensity)
            .collect();
        let Some(pearson) =
            crate::statistics::pearson(&observed_intensities, &theoretical_intensities)
        else {
            statuses.push(status(FragmentMatchState::CorrelationRejected, None));
            continue;
        };
        if pearson < min_correlation {
            statuses.push(status(
                FragmentMatchState::CorrelationRejected,
                Some(pearson),
            ));
            continue;
        }
        let (observed_isotope_offset, observed_peak_index) = observed
            .iter()
            .enumerate()
            .filter(|(offset, _)| envelope[*offset].relative_intensity >= REL_INT_THRESHOLD)
            .filter_map(|(offset, index)| index.map(|index| (offset, index)))
            .max_by(|(_, left), (_, right)| {
                peaks[*left].intensity.total_cmp(&peaks[*right].intensity)
            })
            .unwrap_or((apex_index, apex_peak_index));
        let observed_isotope_offset = envelope[observed_isotope_offset].neutron_offset;
        let display_theoretical_mz = isotope_mz(observed_isotope_offset);
        let observed_monoisotopic_mass = (peaks[observed_peak_index].mz - PROTON)
            * f64::from(fragment.charge)
            - observed_isotope_offset as f64 * NEUTRON_SHIFT;
        matches.push(FragmentMatch {
            ion: fragment.ion,
            number: fragment.number,
            charge: fragment.charge,
            cleavage_after: fragment.cleavage_after,
            neutral_mass: fragment.neutral_mass,
            theoretical_mz: fragment.mz,
            theoretical_apex_mz: isotope_mz(apex.neutron_offset),
            apex_offset: fragment.apex_offset,
            observed_isotope_offset,
            observed_mz: peaks[observed_peak_index].mz,
            observed_neutral_mass: observed_monoisotopic_mass,
            observed_peak_index,
            intensity: peaks[observed_peak_index].intensity,
            quantification: Some(FragmentQuantification {
                observed_intensity: observed_intensities.iter().sum(),
                theoretical_fraction: observed
                    .iter()
                    .zip(&theoretical_intensities)
                    .filter(|(index, _)| index.is_some())
                    .map(|(_, intensity)| intensity)
                    .sum::<f64>()
                    / theoretical_intensities.iter().sum::<f64>(),
                observed_peak_indices: observed.iter().flatten().copied().collect(),
            }),
            envelope_correlation: Some(pearson),
            error_ppm: (peaks[observed_peak_index].mz - display_theoretical_mz)
                / display_theoretical_mz
                * 1_000_000.0,
            neutral_loss: fragment.neutral_loss.clone(),
        });
        statuses.push(status(FragmentMatchState::Matched, Some(pearson)));
    }
    (matches, statuses)
}
