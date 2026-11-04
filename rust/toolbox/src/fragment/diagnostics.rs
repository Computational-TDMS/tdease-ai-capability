use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};

use crate::constants::NEUTRON_SHIFT;

use super::{CandidateModel, FragmentMatch, IonType};

/// Signed fragment-error statistics for one ion series.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SeriesErrorStats {
    pub n: usize,
    pub mean_ppm: f64,
    pub sd_ppm: f64,
}

/// N-terminal (a/b/c) vs C-terminal (x/y/z/z•/z′) signed-error divergence.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SeriesDivergence {
    pub n_minus_c_ppm: f64,
    pub threshold_ppm: f64,
    pub flagged: bool,
}

/// One cross-candidate isotope-borrowing pair.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct IsotopeBorrowPair {
    pub ion: IonType,
    pub number: usize,
    pub charge: u8,
    /// Candidate whose theoretical fragment sits higher; its match is the
    /// borrowed isotope peak.
    pub borrower: String,
    pub lender: String,
    /// Observed neutral-mass gap between the two matches, in neutron units.
    pub neutron_offset: f64,
    /// Gap residual after removing the integer neutron, in ppm of the
    /// borrower's theoretical fragment mass.
    pub residual_ppm: f64,
}

const SERIES_DIVERGENCE_FLOOR_PPM: f64 = 2.0;
const SERIES_DIVERGENCE_MIN_N: usize = 5;

fn error_stats(errors: &[f64]) -> Option<SeriesErrorStats> {
    let n = errors.len();
    if n == 0 {
        return None;
    }
    let mean = errors.iter().sum::<f64>() / n as f64;
    let sd = if n > 1 {
        let variance = errors.iter().map(|e| (e - mean).powi(2)).sum::<f64>() / (n - 1) as f64;
        variance.sqrt()
    } else {
        0.0
    };
    Some(SeriesErrorStats {
        n,
        mean_ppm: mean,
        sd_ppm: sd,
    })
}

pub(super) fn series_error_stats(
    matches: &[FragmentMatch],
) -> (
    BTreeMap<IonType, SeriesErrorStats>,
    Option<SeriesDivergence>,
) {
    let mut per_series = BTreeMap::new();
    for ion in matches
        .iter()
        .map(|matched| matched.ion)
        .collect::<BTreeSet<_>>()
    {
        if let Some(stats) = error_stats(
            &matches
                .iter()
                .filter(|matched| matched.ion == ion)
                .map(|matched| matched.error_ppm)
                .collect::<Vec<_>>(),
        ) {
            per_series.insert(ion, stats);
        }
    }
    let side = |n_terminal: bool| {
        error_stats(
            &matches
                .iter()
                .filter(|matched| {
                    let is_n = matches!(matched.ion, IonType::A | IonType::B | IonType::C);
                    is_n == n_terminal
                })
                .map(|matched| matched.error_ppm)
                .collect::<Vec<_>>(),
        )
    };
    let (Some(n_side), Some(c_side)) = (side(true), side(false)) else {
        return (per_series, None);
    };
    if n_side.n < SERIES_DIVERGENCE_MIN_N || c_side.n < SERIES_DIVERGENCE_MIN_N {
        return (per_series, None);
    }
    let divergence = SeriesDivergence {
        n_minus_c_ppm: n_side.mean_ppm - c_side.mean_ppm,
        threshold_ppm: SERIES_DIVERGENCE_FLOOR_PPM.max(
            2.0 * ((n_side.sd_ppm.powi(2) / n_side.n as f64)
                + (c_side.sd_ppm.powi(2) / c_side.n as f64))
                .sqrt(),
        ),
        flagged: false,
    };
    let flagged = divergence.n_minus_c_ppm.abs() > divergence.threshold_ppm;
    (
        per_series,
        Some(SeriesDivergence {
            flagged,
            ..divergence
        }),
    )
}

/// Detect a candidate claiming an isotope peak already explained by the same
/// ion from another candidate.
pub(super) fn detect_isotope_borrow_pairs<T: AsRef<[FragmentMatch]>>(
    models: &[CandidateModel],
    per_candidate: &[T],
    tolerance_ppm: f64,
) -> Vec<IsotopeBorrowPair> {
    let mut pairs = Vec::new();
    let maps: Vec<BTreeMap<_, _>> = per_candidate
        .iter()
        .map(|matches| keyed_matches(matches.as_ref()).collect())
        .collect();
    for (ai, amap) in maps.iter().enumerate() {
        for (bi, bmap) in maps.iter().enumerate().skip(ai + 1) {
            for (key, &mb) in bmap {
                let Some(&ma) = amap.get(key) else {
                    continue;
                };
                if ma.observed_peak_index == mb.observed_peak_index {
                    continue;
                }
                let (hi, lo, hi_model, lo_model) = if ma.theoretical_mz > mb.theoretical_mz {
                    (ma, mb, ai, bi)
                } else {
                    (mb, ma, bi, ai)
                };
                let gap_neutral = (hi.observed_mz - lo.observed_mz) * f64::from(key.2);
                let offset = gap_neutral / NEUTRON_SHIFT;
                let k = offset.round();
                if k < 1.0 || (offset - k).abs() > 0.25 {
                    continue;
                }
                let residual_ppm =
                    (gap_neutral - k * NEUTRON_SHIFT) / hi.neutral_mass.abs().max(1.0) * 1e6;
                if residual_ppm.abs() > tolerance_ppm {
                    continue;
                }
                pairs.push(IsotopeBorrowPair {
                    ion: key.0,
                    number: key.1,
                    charge: key.2,
                    borrower: models[hi_model].id.clone(),
                    lender: models[lo_model].id.clone(),
                    neutron_offset: offset,
                    residual_ppm,
                });
            }
        }
    }
    pairs.sort_by(|a, b| {
        (a.ion, a.number, a.charge, &a.borrower).cmp(&(b.ion, b.number, b.charge, &b.borrower))
    });
    pairs
}

fn keyed_matches(
    matches: &[FragmentMatch],
) -> impl Iterator<Item = ((IonType, usize, u8), &FragmentMatch)> {
    matches
        .iter()
        .map(|matched| ((matched.ion, matched.number, matched.charge), matched))
}
