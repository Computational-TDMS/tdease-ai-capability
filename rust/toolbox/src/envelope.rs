/*!
Observed-vs-theoretical isotope envelope matching over an mzML MS1 window.

This module now delegates its chemistry and its theoretical distribution to the
crate's single sources: `proforma` (exact composition) and `theory` (BRAIN
envelope). There is no second convolution here — a given ProForma produces the
identical theoretical distribution through `isotopeEnvelope` and through
`calibrateMassShift`.
*/

use mzdata::io::SpectrumSource;
use mzdata::prelude::SpectrumLike;
use serde::{Deserialize, Serialize};

use crate::constants::{NEUTRON_SHIFT, PROTON};
use crate::mzml::open_reader;
use crate::proforma::{chemical_from_counts, parse_proforma, ExternalModification};
use crate::theory::{anchored_envelope, intensity_by_neutron_index, theoretical_envelope};
use crate::{MassSpecError, Peak};

const MAX_ISOTOPES: usize = 64;
const AUTO_RELATIVE_APEX_CUTOFF: f64 = 0.01;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct IsotopeEnvelopeQuery {
    pub mass: Option<f64>,
    pub proforma: Option<String>,
    /// Named-modification injection for the ProForma parser (name -> net formula).
    #[serde(default)]
    pub external_mods: Vec<ExternalModification>,
    pub charge: u8,
    pub rt_min: f64,
    pub rt_max: f64,
    pub tolerance_ppm: f64,
    pub isotopes: usize,
    pub max_offset: i32,
    /// Strict theory-locked collection: match observed peaks against the
    /// theoretical ladder positions only — no per-scan seed-residual search.
    /// The seed freedom lets an off-ladder noise peak re-anchor a whole scan
    /// and mixes misaligned picks into the aggregate; locking to theory
    /// keeps noise out and lets slightly-offset peaks still land in their
    /// own bin (they stay within the tolerance window of THEIR target).
    /// Measured on real top-down data (f198): pearson 0.671 -> 0.847/0.884
    /// (±10/±3 ppm lock). Default false keeps the interactive `envelope`
    /// behavior unchanged; the refine/calibration paths set it.
    #[serde(default)]
    pub theory_locked: bool,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct EnvelopeModel {
    pub kind: String,
    pub neutral_mass: f64,
    /// Exact composition from the ProForma parser (or averagine for a bare mass).
    pub composition: Vec<crate::isotopic::ElementCount>,
    /// BRAIN-derived relative intensity per neutron index (legacy shape).
    pub theoretical_pattern: Vec<f64>,
    /// Absolute m/z per neutron index, anchored to the model's neutral mass —
    /// the same coordinate frame as the observed envelope.
    pub theoretical_mz: Vec<f64>,
    pub warnings: Vec<String>,
}

/// Legacy CHNOSPBr composition shape, kept for backward compatibility. New
/// code should use `Vec<ElementCount>` (via `proforma`/`theory`), which carries
/// every element rather than just C H N O S P Br.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "UPPERCASE")]
pub struct ElementalComposition {
    pub c: i32,
    pub h: i32,
    pub n: i32,
    pub o: i32,
    pub s: i32,
    pub p: i32,
    pub br: i32,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RtAggregation {
    pub scan_count: usize,
    pub retention_time_unit: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct MatchedIsotope {
    pub isotope: usize,
    pub theoretical_mz: f64,
    pub observed_mz: f64,
    pub error_ppm: f64,
    pub intensity: f32,
    pub relative_theoretical: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct EnvelopeEvaluation {
    pub charge: u8,
    pub offset: i32,
    pub monoisotopic_mass: f64,
    pub isotopes_matched: usize,
    pub isotopes_considered: usize,
    pub pearson: Option<f64>,
    pub summed_intensity: f64,
    pub matched_isotopes: Vec<MatchedIsotope>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct IsotopeEnvelopeResult {
    pub model: EnvelopeModel,
    pub rt_aggregation: RtAggregation,
    pub evaluations: Vec<EnvelopeEvaluation>,
    pub best: EnvelopeEvaluation,
}

pub fn query_isotope_envelope_mzml(
    bytes: &[u8],
    query: &IsotopeEnvelopeQuery,
) -> Result<IsotopeEnvelopeResult, MassSpecError> {
    validate_query(query)?;
    let mass = query
        .mass
        .or_else(|| {
            query
                .proforma
                .as_ref()
                .and_then(|_| build_model(query).ok().map(|m| m.neutral_mass))
        })
        .ok_or_else(|| {
            MassSpecError::InvalidQuery("isotopeEnvelope mass unavailable".to_owned())
        })?;
    let mut resolved = query.clone();
    resolved.isotopes = resolve_isotope_count(query, mass)?;
    let model = build_model(&resolved)?;
    let (peaks, scan_count) = aggregate_window(bytes, &resolved, model.neutral_mass)?;
    let mut evaluations: Vec<_> = (-resolved.max_offset..=resolved.max_offset)
        .map(|offset| evaluate(&peaks, &resolved, &model, offset))
        .collect();
    evaluations.sort_by_key(|left| left.offset);
    let best = evaluations
        .iter()
        .max_by(|left, right| {
            left.isotopes_matched
                .cmp(&right.isotopes_matched)
                .then_with(|| {
                    left.pearson
                        .unwrap_or(f64::NEG_INFINITY)
                        .total_cmp(&right.pearson.unwrap_or(f64::NEG_INFINITY))
                })
                .then_with(|| left.summed_intensity.total_cmp(&right.summed_intensity))
                .then_with(|| right.offset.abs().cmp(&left.offset.abs()))
        })
        .cloned()
        .ok_or_else(|| MassSpecError::InvalidQuery("no envelope offsets requested".to_owned()))?;
    Ok(IsotopeEnvelopeResult {
        model,
        rt_aggregation: RtAggregation {
            scan_count,
            retention_time_unit: "minute".to_owned(),
        },
        evaluations,
        best,
    })
}

/// Evaluate one already-selected MS1 spectrum against the same theoretical
/// envelope used by the window query. This is the IsoForma precursor gate:
/// the MS2 header points to one parent MS1, so no chromatographic aggregation
/// is appropriate here.
pub(crate) fn evaluate_isotope_envelope_peaks(
    peaks: &[Peak],
    query: &IsotopeEnvelopeQuery,
) -> Result<EnvelopeEvaluation, MassSpecError> {
    validate_query(query)?;
    let mass = query
        .mass
        .or_else(|| {
            query
                .proforma
                .as_ref()
                .and_then(|_| build_model(query).ok().map(|model| model.neutral_mass))
        })
        .ok_or_else(|| {
            MassSpecError::InvalidQuery("isotopeEnvelope mass unavailable".to_owned())
        })?;
    let mut resolved = query.clone();
    resolved.isotopes = resolve_isotope_count(query, mass)?;
    let model = build_model(&resolved)?;
    Ok(evaluate(peaks, &resolved, &model, 0))
}

fn validate_query(query: &IsotopeEnvelopeQuery) -> Result<(), MassSpecError> {
    if query.mass.is_none() == query.proforma.is_none() {
        return Err(MassSpecError::InvalidQuery(
            "isotopeEnvelope requires exactly one of mass or proforma".to_owned(),
        ));
    }
    if query.charge == 0
        || !query.rt_min.is_finite()
        || !query.rt_max.is_finite()
        || query.rt_min > query.rt_max
        || !query.tolerance_ppm.is_finite()
        || query.tolerance_ppm <= 0.0
        || query.isotopes > MAX_ISOTOPES
        || query.max_offset < 0
    {
        return Err(MassSpecError::InvalidQuery(
            "invalid isotopeEnvelope charge, RT range, tolerance, isotope count, or offset"
                .to_owned(),
        ));
    }
    Ok(())
}

fn build_model(query: &IsotopeEnvelopeQuery) -> Result<EnvelopeModel, MassSpecError> {
    let (kind, neutral_mass, composition, warnings) = if let Some(mass) = query.mass {
        if !mass.is_finite() || mass <= 0.0 {
            return Err(MassSpecError::InvalidQuery(
                "isotopeEnvelope mass must be positive".to_owned(),
            ));
        }
        (
            "averagine".to_owned(),
            mass,
            crate::theory::averagine_composition(mass),
            Vec::new(),
        )
    } else {
        let parsed = parse_proforma(
            query.proforma.as_deref().unwrap_or_default(),
            &query.external_mods,
        )?;
        (
            "proformaComposition".to_owned(),
            parsed.monoisotopic_mass,
            chemical_from_counts(&parsed.composition),
            parsed.warnings,
        )
    };

    // Single BRAIN kernel: exact composition -> anchored theoretical envelope.
    let env = theoretical_envelope(&composition, i32::from(query.charge), 0.99, 0.001);
    let anchored = anchored_envelope(&env, neutral_mass);
    let (theoretical_pattern, theoretical_mz) =
        intensity_by_neutron_index(&anchored, query.isotopes);

    Ok(EnvelopeModel {
        kind,
        neutral_mass,
        composition: crate::proforma::element_counts_from_chemical(&composition),
        theoretical_pattern,
        theoretical_mz,
        warnings,
    })
}

fn aggregate_window(
    bytes: &[u8],
    query: &IsotopeEnvelopeQuery,
    mass: f64,
) -> Result<(Vec<Peak>, usize), MassSpecError> {
    // Exposed pub(crate) so the mass-shift refiner can reuse the same MS1
    // window aggregation instead of duplicating it.
    aggregate_window_inner(bytes, query, mass)
}

pub fn aggregate_window_inner(
    bytes: &[u8],
    query: &IsotopeEnvelopeQuery,
    mass: f64,
) -> Result<(Vec<Peak>, usize), MassSpecError> {
    aggregate_windows_inner(bytes, &[(query, mass)])
        .into_iter()
        .next()
        .unwrap_or_else(|| {
            Err(MassSpecError::InvalidQuery(
                "empty envelope batch result".into(),
            ))
        })
}

/// Aggregate several independent isotope windows in one mzML traversal.
/// Each spectrum is decoded once, then tested only against requests whose RT
/// windows contain it. The per-request matching and apex selection are exactly
/// the same as [`aggregate_window_inner`].
pub fn aggregate_windows_inner(
    bytes: &[u8],
    requests: &[(&IsotopeEnvelopeQuery, f64)],
) -> Vec<Result<(Vec<Peak>, usize), MassSpecError>> {
    let prepared: Vec<Result<Vec<f64>, MassSpecError>> = requests
        .iter()
        .map(|(query, mass)| envelope_targets(query, *mass))
        .collect();
    let mut matches: Vec<Vec<(f64, CoherentScanMatch)>> =
        (0..requests.len()).map(|_| Vec::new()).collect();
    let mut reader = match open_reader(bytes) {
        Ok(reader) => reader,
        Err(error) => {
            let message = error.to_string();
            return requests
                .iter()
                .map(|_| Err(MassSpecError::InvalidQuery(message.clone())))
                .collect();
        }
    };
    for mut spectrum in reader.iter() {
        if spectrum.ms_level() != 1 {
            continue;
        }
        let rt = spectrum.start_time();
        let active: Vec<usize> = requests
            .iter()
            .enumerate()
            .filter(|(index, (query, _))| {
                prepared[*index].is_ok() && rt >= query.rt_min && rt <= query.rt_max
            })
            .map(|(index, _)| index)
            .collect();
        if active.is_empty() {
            continue;
        }
        if let Err(error) = spectrum.try_build_peaks() {
            let message = error.to_string();
            return requests
                .iter()
                .map(|_| Err(MassSpecError::PeakDecode(message.clone())))
                .collect();
        }
        let centroid_peaks: Vec<Peak> = spectrum
            .peaks()
            .iter()
            .map(|peak| Peak {
                mz: peak.mz,
                intensity: peak.intensity,
            })
            .collect();
        for index in active {
            let (query, _) = requests[index];
            if let Ok(targets) = &prepared[index] {
                if let Some(matched) = coherent_scan_match(
                    &centroid_peaks,
                    targets,
                    query.tolerance_ppm,
                    query.theory_locked,
                ) {
                    matches[index].push((rt, matched));
                }
            }
        }
    }
    prepared
        .into_iter()
        .zip(matches)
        .map(|(targets, scan_matches)| match targets {
            Ok(targets) => aggregate_scan_matches(scan_matches, targets.len()),
            Err(error) => Err(error),
        })
        .collect()
}

fn envelope_targets(query: &IsotopeEnvelopeQuery, mass: f64) -> Result<Vec<f64>, MassSpecError> {
    let isotope_count = resolve_isotope_count(query, mass)?;
    Ok((-query.max_offset..isotope_count as i32 + query.max_offset)
        .map(|isotope| {
            (mass + f64::from(isotope) * NEUTRON_SHIFT + f64::from(query.charge) * PROTON)
                / f64::from(query.charge)
        })
        .collect())
}

fn aggregate_scan_matches(
    mut scan_matches: Vec<(f64, CoherentScanMatch)>,
    target_count: usize,
) -> Result<(Vec<Peak>, usize), MassSpecError> {
    if scan_matches.is_empty() {
        return Err(MassSpecError::InvalidQuery(
            "no coherent MS1 isotope chain in isotopeEnvelope RT window".to_owned(),
        ));
    }
    let apex_index = scan_matches
        .iter()
        .enumerate()
        .max_by(|(_, left), (_, right)| left.1.total_intensity.total_cmp(&right.1.total_intensity))
        .map(|(index, _)| index)
        .unwrap_or(0);
    let apex_intensity = scan_matches[apex_index].1.total_intensity;
    let min_intensity = apex_intensity * 0.05;
    let mut left = apex_index;
    while left > 0 && scan_matches[left - 1].1.total_intensity >= min_intensity {
        left -= 1;
    }
    let mut right = apex_index + 1;
    while right < scan_matches.len() && scan_matches[right].1.total_intensity >= min_intensity {
        right += 1;
    }
    let selected: Vec<_> = scan_matches.drain(left..right).collect();
    let scan_count = selected.len();
    let mut accumulated = vec![(0.0_f64, 0.0_f32); target_count];
    for (_, matched) in selected {
        for (index, peak) in matched.peaks.into_iter().enumerate() {
            if let Some(peak) = peak {
                accumulated[index].0 += peak.mz * f64::from(peak.intensity);
                accumulated[index].1 += peak.intensity;
            }
        }
    }
    let peaks = accumulated
        .into_iter()
        .filter_map(|(weighted_mz, intensity)| {
            (intensity > 0.0).then_some(Peak {
                mz: weighted_mz / f64::from(intensity),
                intensity,
            })
        })
        .collect();
    Ok((peaks, scan_count))
}

fn resolve_isotope_count(query: &IsotopeEnvelopeQuery, mass: f64) -> Result<usize, MassSpecError> {
    if query.isotopes > 0 {
        return Ok(query.isotopes);
    }
    let composition = if let Some(proforma) = query.proforma.as_deref() {
        let parsed = parse_proforma(proforma, &query.external_mods)?;
        chemical_from_counts(&parsed.composition)
    } else {
        crate::theory::averagine_composition(mass)
    };
    // Use the same finite support as `theory`. Passing (1.0, 0.0) through
    // mzcore's truncation API retains only M0, which silently reduces the
    // automatic query to one isotope and makes a >=3-peak validation
    // impossible for large proteoforms.
    let envelope = theoretical_envelope(&composition, i32::from(query.charge), 0.99, 0.001);
    let anchored = anchored_envelope(&envelope, mass);
    Ok(relative_apex_isotope_count(
        &anchored,
        AUTO_RELATIVE_APEX_CUTOFF,
    ))
}

fn relative_apex_isotope_count(
    anchored: &crate::theory::AnchoredEnvelope,
    relative_cutoff: f64,
) -> usize {
    let apex = anchored
        .peaks
        .iter()
        .map(|peak| peak.intensity)
        .fold(0.0_f64, f64::max);
    if apex <= 0.0 {
        return 1;
    }
    let threshold = apex * relative_cutoff;
    anchored
        .peaks
        .iter()
        .zip(&anchored.center_masses)
        .filter(|(peak, _)| peak.intensity >= threshold)
        .map(|(_, center_mass)| {
            ((center_mass - anchored.anchor_mass) / NEUTRON_SHIFT)
                .round()
                .max(0.0) as usize
        })
        .max()
        .map_or(1, |last| (last + 1).min(MAX_ISOTOPES))
}

#[derive(Debug, Clone)]
struct CoherentScanMatch {
    peaks: Vec<Option<Peak>>,
    total_intensity: f32,
    matched_count: usize,
    residual: f64,
}

fn coherent_scan_match(
    peaks: &[Peak],
    targets: &[f64],
    tolerance_ppm: f64,
    theory_locked: bool,
) -> Option<CoherentScanMatch> {
    // Theory-locked collection tries only the theoretical positions; the
    // seed-residual search below is the per-scan re-anchoring freedom that
    // admits off-ladder noise and misaligned picks.
    let mut seeds = vec![0.0_f64];
    if !theory_locked {
        for target in targets {
            let window = target * tolerance_ppm / 1_000_000.0;
            let mut nearby: Vec<_> = peaks
                .iter()
                .filter(|peak| (peak.mz - target).abs() <= window)
                .collect();
            nearby.sort_by(|a, b| b.intensity.total_cmp(&a.intensity));
            seeds.extend(nearby.into_iter().take(3).map(|peak| peak.mz - target));
        }
    }

    let local_ppm = tolerance_ppm.clamp(5.0, 20.0);
    let mut best: Option<CoherentScanMatch> = None;
    for residual in seeds {
        let mut used = vec![false; peaks.len()];
        let mut matched = vec![None; targets.len()];
        let mut total_intensity = 0.0_f32;
        let mut matched_count = 0usize;
        for (target_index, target) in targets.iter().enumerate() {
            let shifted = target + residual;
            let window = shifted * local_ppm / 1_000_000.0;
            let candidate = peaks
                .iter()
                .enumerate()
                .filter(|(index, peak)| !used[*index] && (peak.mz - shifted).abs() <= window)
                .min_by(|(_, a), (_, b)| (a.mz - shifted).abs().total_cmp(&(b.mz - shifted).abs()));
            if let Some((index, peak)) = candidate {
                used[index] = true;
                let value = *peak;
                total_intensity += value.intensity.sqrt();
                matched[target_index] = Some(value);
                matched_count += 1;
            }
        }
        let candidate = CoherentScanMatch {
            peaks: matched,
            total_intensity,
            matched_count,
            residual,
        };
        let replace = best.as_ref().is_none_or(|current| {
            candidate.matched_count > current.matched_count
                || (candidate.matched_count == current.matched_count
                    && candidate.total_intensity > current.total_intensity)
                || (candidate.matched_count == current.matched_count
                    && candidate.total_intensity == current.total_intensity
                    && candidate.residual.abs() < current.residual.abs())
        });
        if replace {
            best = Some(candidate);
        }
    }
    best.filter(|matched| matched.matched_count >= 2 || targets.len() == 1)
}

fn evaluate(
    peaks: &[Peak],
    query: &IsotopeEnvelopeQuery,
    model: &EnvelopeModel,
    offset: i32,
) -> EnvelopeEvaluation {
    let monoisotopic_mass = model.neutral_mass + f64::from(offset) * NEUTRON_SHIFT;
    let mut observed = Vec::with_capacity(query.isotopes);
    let mut matched_isotopes = Vec::new();
    for isotope in 0..query.isotopes {
        // BRAIN-anchored absolute m/z (not fixed spacing) + the integer offset.
        let base_mz = model
            .theoretical_mz
            .get(isotope)
            .copied()
            .unwrap_or_else(|| {
                (monoisotopic_mass
                    + isotope as f64 * NEUTRON_SHIFT
                    + f64::from(query.charge) * PROTON)
                    / f64::from(query.charge)
            });
        let target = base_mz + f64::from(offset) * NEUTRON_SHIFT / f64::from(query.charge);
        let window = target * query.tolerance_ppm / 1_000_000.0;
        let matched = peaks
            .iter()
            .filter(|peak| (peak.mz - target).abs() <= window)
            .min_by(|left, right| {
                (left.mz - target)
                    .abs()
                    .total_cmp(&(right.mz - target).abs())
            });
        let relative_theoretical = model
            .theoretical_pattern
            .get(isotope)
            .copied()
            .unwrap_or(0.0);
        observed.push(matched.map_or(0.0, |peak| f64::from(peak.intensity)));
        if let Some(peak) = matched {
            matched_isotopes.push(MatchedIsotope {
                isotope,
                theoretical_mz: target,
                observed_mz: peak.mz,
                error_ppm: (peak.mz - target) / target * 1_000_000.0,
                intensity: peak.intensity,
                relative_theoretical,
            });
        }
    }
    // Pearson over matched isotope pairs only. Padding unmatched isotopologue
    // positions with 0 pulled the correlation toward a false "missing shape"
    // agreement; only the peaks actually matched to the candidate are compared.
    let pearson = if matched_isotopes.len() >= 3 {
        let obs_m: Vec<f64> = matched_isotopes
            .iter()
            .map(|m| f64::from(m.intensity))
            .collect();
        let theo_m: Vec<f64> = matched_isotopes
            .iter()
            .map(|m| m.relative_theoretical)
            .collect();
        crate::statistics::pearson(&obs_m, &theo_m)
    } else {
        None
    };
    EnvelopeEvaluation {
        charge: query.charge,
        offset,
        monoisotopic_mass,
        isotopes_matched: matched_isotopes.len(),
        isotopes_considered: query.isotopes,
        pearson,
        summed_intensity: observed.iter().sum(),
        matched_isotopes,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_percent_theoretical_support_stops_before_the_noise_tail() {
        let anchored = crate::theory::AnchoredEnvelope {
            peaks: vec![
                crate::isotopic::TheoreticalPeakPoint {
                    mz: 100.0,
                    intensity: 0.001,
                },
                crate::isotopic::TheoreticalPeakPoint {
                    mz: 101.0,
                    intensity: 0.02,
                },
                crate::isotopic::TheoreticalPeakPoint {
                    mz: 102.0,
                    intensity: 1.0,
                },
                crate::isotopic::TheoreticalPeakPoint {
                    mz: 103.0,
                    intensity: 0.01,
                },
                crate::isotopic::TheoreticalPeakPoint {
                    mz: 104.0,
                    intensity: 0.009,
                },
            ],
            center_masses: vec![100.0, 101.0, 102.0, 103.0, 104.0],
            charge: 1,
            anchor_mass: 100.0,
        };
        assert_eq!(relative_apex_isotope_count(&anchored, 0.01), 4);
    }

    #[test]
    fn automatic_support_keeps_large_proteoform_isotope_ladder() {
        let query = IsotopeEnvelopeQuery {
            mass: None,
            proforma: Some("[Acetyl]-SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRTVLKTFLENVIRDSVTYTEHARRKTVTAMDVVYALKRQGRTLYGFGG".to_owned()),
            external_mods: Vec::new(),
            charge: 14,
            rt_min: 0.0,
            rt_max: 1.0,
            tolerance_ppm: 10.0,
            isotopes: 0,
            max_offset: 0,
            theory_locked: false,
        };
        assert!(resolve_isotope_count(&query, 11_391.0).unwrap() >= 10);
    }

    #[test]
    fn coherent_chain_beats_nearer_isolated_noise() {
        let targets = [100.0, 100.5, 101.0];
        let peaks = [
            Peak {
                mz: 100.0,
                intensity: 1.0,
            },
            Peak {
                mz: 100.5,
                intensity: 1.0,
            },
            Peak {
                mz: 101.0,
                intensity: 1.0,
            },
            Peak {
                mz: 100.01,
                intensity: 100.0,
            },
            Peak {
                mz: 100.51,
                intensity: 80.0,
            },
            Peak {
                mz: 101.01,
                intensity: 60.0,
            },
        ];
        let matched = coherent_scan_match(&peaks, &targets, 150.0, false).unwrap();
        let mz: Vec<_> = matched
            .peaks
            .into_iter()
            .flatten()
            .map(|peak| peak.mz)
            .collect();
        assert_eq!(mz, vec![100.01, 100.51, 101.01]);
    }

    #[test]
    fn proforma_formula_preserves_exact_composition_and_mass() {
        let q = IsotopeEnvelopeQuery {
            mass: None,
            proforma: Some("[Formula:C2H2O]-PEPTIDE".to_owned()),
            external_mods: Vec::new(),
            charge: 1,
            rt_min: -1.0,
            rt_max: 1.0,
            tolerance_ppm: 10.0,
            isotopes: 8,
            max_offset: 2,
            theory_locked: false,
        };
        let model = build_model(&q).unwrap();
        assert_eq!(model.kind, "proformaComposition");
        // PEPTIDE carbons (34) + acetyl (2) = 36.
        let c = model
            .composition
            .iter()
            .find(|e| e.element == "C")
            .map(|e| e.count)
            .unwrap();
        assert_eq!(c, 36);
    }

    #[test]
    fn bromine_formula_preserves_the_near_equal_m_plus_two_signature() {
        let q = IsotopeEnvelopeQuery {
            mass: None,
            proforma: Some("[Formula:Br-H]-PEPTIDE".to_owned()),
            external_mods: Vec::new(),
            charge: 1,
            rt_min: -1.0,
            rt_max: 1.0,
            tolerance_ppm: 10.0,
            isotopes: 8,
            max_offset: 2,
            theory_locked: false,
        };
        let model = build_model(&q).unwrap();
        // M+2 (81Br) must be prominent relative to M0. On a peptide of several
        // dozen carbons the M+2 carries both the 81Br contribution and the
        // natural C13 M+2, so it can exceed 1.0; the point is it is far above
        // a pure-C13 M+2 (~0.1) — i.e. the chemical_elements 0.7.0 Br fix is
        // preserved.
        let m0 = model.theoretical_pattern[0];
        let m2 = model.theoretical_pattern.get(2).copied().unwrap_or(0.0);
        let ratio = m2 / m0.max(1e-12);
        assert!(ratio > 0.6 && ratio < 1.5, "M+2/M0 {ratio}");
        // The M+2 m/z must sit near the Br doublet spacing (~2.008/1), not 2*NEUTRON.
        let mz2 = model.theoretical_mz[2];
        let mz0 = model.theoretical_mz[0];
        let spacing = mz2 - mz0;
        assert!(
            (spacing - 2.0 * 1.003_354_837_8).abs() < 0.02,
            "Br M+2 spacing {spacing}"
        );
    }
}
