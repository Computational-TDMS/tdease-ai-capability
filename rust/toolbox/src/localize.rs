/*!
PrSM mass-shift localization by carrier-presence window sweep over TopFD's
deconvoluted MS2 spectra (msalign).

TopPIC's open search reports an unexpected mass shift as a residue *range*
(its diagonal alignment only bounds the gap between matched segments, and the
MIScore refinement only localizes shifts already in the PTM list). This
module closes that gap on the same spectra TopPIC searched -- the `_ms2.msalign`
TopFD output -- where every peak is already a deconvoluted neutral monoisotopic
mass, so matching happens directly in neutral mass space with no charge
enumeration layer.

The physics on deconvoluted peaks is carrier presence: a fragment whose
cleavage spans the true site must be observed at its theoretical mass PLUS the
shift, and a fragment that stops short of the site must be observed at its
plain mass. One fragment cannot be both (the shift moves it by whole Da), so
every observable fragment casts exactly one vote. Sweeping a unit window over
every cleavage position turns the votes into a localization gradient read by
two independent ladders:

- the N-terminal ladder (b/c; a cleavage `cut` spans the site at residue
  `p` iff `cut >= p`) saturates only from the site upward (below it,
  non-carriers observed at their plain masses dilute it),
- the C-terminal ladder (y/z; spans iff `cut <= p`) holds only up to the
  site (above it, plain-matched non-carriers dilute it).

Their crossover is the site; each ladder alone would plateau outward, but
their agreement cannot. (TopFD absorbs the 81Br isotope partner into the same
deconvoluted envelope, so on msalign the partner does not exist as a separate
peak and carrier presence is the honest observable.)

Self-diagnosis is built in: fragments below the window (`pre` class) share the
parameterization, so a healthy control pool beside a collapsed spanning class
means the run is UNINFORMATIVE, never negative.
*/

use serde::{Deserialize, Serialize};

use crate::constants::H2O;
use crate::fragment::{residue_mass, MassModification};
use crate::MassSpecError;

/// Fragment-type offsets over the plain residue sums (Da, neutral frame).
const C_OFFSET: f64 = 17.026_549;
const ZDOT_OFFSET: f64 = 1.991_840_6; // z_k = suffix_residues + ZDOT_OFFSET
const ZP1_OFFSET: f64 = 2.999_665_6; // z+1_k = suffix_residues + ZP1_OFFSET

/// Peaks below this fraction of the base peak never participate (both classes).
/// TopFD's own S/N convention, not a tuned verdict rule.
const INTENSITY_FLOOR: f64 = 0.01;
/// Laplace smoothing for the per-window noise-rate estimate (standard
/// statistical practice, not a case-tuned threshold).
const NOISE_SMOOTHING: f64 = 0.5;
/// A ladder's maximum must reach this -log10 p-value (p < 0.01) before any
/// site is called from it.
const EVIDENCE_BAR: f64 = 2.0;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PrsmLocalizeQuery {
    /// Scan numbers (the msalign SCANS= field). Several scans of the same
    /// precursor combine: per-position log-likelihoods sum across scans, so
    /// an observation gap in one spectrum is closed by another.
    pub scans: Vec<String>,
    pub sequence: String,
    /// The unknown shift to localize (Da, monoisotopic).
    pub delta_mass: f64,
    /// Pinned known modifications (1-based positions or "nTerm"/"cTerm").
    /// Pinned events only -- a ranged family must be decomposed first (see
    /// `_physics._parse_extra` for the measured false-negative a range causes).
    #[serde(default)]
    pub modifications: Vec<MassModification>,
    /// Match tolerance in ppm (default 8).
    #[serde(default = "default_ppm")]
    pub tolerance_ppm: f64,
    #[serde(default = "default_mass_lo")]
    pub mass_min: f64,
    #[serde(default = "default_mass_hi")]
    pub mass_max: f64,
    /// Restrict the sweep to a TopPIC-style 1-based inclusive range
    /// (e.g. 33-91). Default: the whole sequence -- every cleavage position
    /// is a candidate and the site must be *searched out*, never assumed.
    #[serde(default)]
    pub range: Option<(usize, usize)>,
}

fn default_ppm() -> f64 {
    8.0
}
fn default_mass_lo() -> f64 {
    50.0
}
fn default_mass_hi() -> f64 {
    100_000.0
}

/// Evidence summary over all fragments (the model's measured inputs).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct LocalizeEvidence {
    /// Fragments observed at the shifted mass (carriers under the hypothesis
    /// that the site lies inside them).
    pub shifted: usize,
    /// Fragments observed only at the plain mass (non-carriers under the
    /// same hypothesis).
    pub plain_only: usize,
    /// Fragments observed at BOTH masses -- ambiguous, and the basis of the
    /// measured noise rate.
    pub both: usize,
    /// Fragments observed nowhere (uninformative, excluded from the model).
    pub unobserved: usize,
    /// Measured probability that a non-carrier fragment spuriously matches
    /// at the shifted mass: both / (plain_only + both), Laplace-smoothed.
    pub noise_rate: f64,
}

/// Log-likelihood of one candidate residue position. The score is relative
/// to the best position, so the maximum is 0 and every other position is
/// negative; a flat top is an unresolved interval, not a bug.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PositionStats {
    /// 1-based residue position.
    pub position: usize,
    pub residue: char,
    /// log10 likelihood relative to the best position (<= 0).
    pub relative_log10_likelihood: f64,
    /// Carrier / non-carrier split among fragments that contain this
    /// position and are observable.
    pub containing_observable: usize,
    pub containing_shifted: usize,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct LocalizeRow {
    /// The msalign scan this row was matched on.
    pub scan: String,
    pub label: String,
    /// The deconvoluted peak's charge (as TopFD assigned it).
    pub charge: Option<u8>,
    /// `span` = observed at the shifted mass (a carrier vote);
    /// `spanPlain` = spanning but observed at the plain mass (a non-carrier
    /// vote); `preShifted` = does not contain the site yet matches shifted
    /// (noise, visible for audit).
    pub class: String,
    pub theoretical_mass: f64,
    pub observed_mass: f64,
    pub error_ppm: f64,
    pub intensity: f32,
    /// TopFD's deconvolution quality for the matched peak (msalign 4th
    /// column): TopFD's own per-peak confidence, NOT an envelope
    /// correlation — provenance stays labeled.
    pub topfd_score: Option<f32>,
}

/// One input scan's contribution to the combined localization.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ScanLocalizeSummary {
    pub scan: String,
    /// (measured precursor mass, theoretical sequence+mods+delta mass, ppm).
    pub precursor_check: Option<(f64, f64, f64)>,
    pub shifted: usize,
    pub plain_only: usize,
    pub both: usize,
    pub unobserved: usize,
    pub noise_rate: f64,
    pub evidence_log10: Option<f64>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PrsmLocalizeResult {
    pub sequence: String,
    pub delta_mass: f64,
    /// Per-scan details; the position curve and site call below are the
    /// combined (summed-likelihood) result over all scans.
    pub scans: Vec<ScanLocalizeSummary>,
    pub evidence: LocalizeEvidence,
    /// The per-position likelihood curve, best position first-normalized.
    pub positions: Vec<PositionStats>,
    /// Residue interval the likelihood pins the shift into (1-based
    /// inclusive): the flat top of the curve.
    pub site_interval: (usize, usize),
    /// Single residue when the likelihood top is unique.
    pub best_site: Option<usize>,
    /// log10 likelihood ratio of the best position against the no-shift null
    /// (all shifted matches being noise) -- the localization's evidence.
    pub evidence_log10: Option<f64>,
    /// ProForma with the localized shift `[+mass@pos]` / `[+mass|a-b]`.
    pub proforma: String,
    /// Matched rows classified at the best position, for manual inspection.
    pub best_rows: Vec<LocalizeRow>,
    pub notes: Vec<String>,
}

/// One TopFD deconvoluted peak — every field of the msalign peak line is
/// kept: neutral monoisotopic mass, intensity, charge, and TopFD's
/// deconvolution quality score (4th column). TopFD is the deconvoluter in
/// this pipeline, so this is our `DeconvolutedPeak` (IP's carries Corr/Dist
/// from ITS deconvoluter; ours carries TopFD's score — same role, different
/// provenance, never silently interchanged).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DeconvolutedPeak {
    pub mass: f64,
    pub intensity: f32,
    pub charge: Option<u8>,
    pub score: Option<f32>,
}

/// One deconvoluted msalign spectrum.
#[derive(Debug, Clone)]
pub struct MsAlignSpectrum {
    pub scans: String,
    pub level: u8,
    pub activation: Option<String>,
    pub precursor_mass: Option<f64>,
    pub precursor_charge: Option<u8>,
    /// Deconvoluted peaks (all four msalign columns).
    pub peaks: Vec<DeconvolutedPeak>,
}

/// Extract one spectrum by scan number from TopFD's `_ms2.msalign` text.
/// Peak lines are TAB-separated `mass intensity charge [score]`; lines with
/// `=` are header fields; `#` starts a comment.
pub fn parse_msalign_spectrum(bytes: &[u8], scan: &str) -> Result<MsAlignSpectrum, MassSpecError> {
    let text = std::str::from_utf8(bytes)
        .map_err(|error| MassSpecError::InvalidQuery(format!("msalign is not utf-8: {error}")))?;
    let mut current: Option<MsAlignSpectrum> = None;
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if line == "BEGIN IONS" {
            current = Some(MsAlignSpectrum {
                scans: String::new(),
                level: 2,
                activation: None,
                precursor_mass: None,
                precursor_charge: None,
                peaks: Vec::new(),
            });
            continue;
        }
        if line == "END IONS" {
            if let Some(spectrum) = current.take() {
                if spectrum.scans == scan {
                    return Ok(spectrum);
                }
            }
            continue;
        }
        let Some(spectrum) = current.as_mut() else {
            continue;
        };
        if let Some((key, value)) = line.split_once('=') {
            match key.trim() {
                "SCANS" => spectrum.scans = value.trim().to_owned(),
                "LEVEL" => {
                    spectrum.level = value.trim().parse().unwrap_or(2);
                }
                "ACTIVATION" => spectrum.activation = Some(value.trim().to_owned()),
                "PRECURSOR_MASS" => {
                    spectrum.precursor_mass = value.trim().parse().ok();
                }
                "PRECURSOR_CHARGE" => {
                    spectrum.precursor_charge = value
                        .trim()
                        .parse::<i64>()
                        .ok()
                        .and_then(|c| u8::try_from(c.unsigned_abs()).ok());
                }
                _ => {}
            }
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() >= 2 {
            let (Ok(mass), Ok(intensity)) = (
                columns[0].trim().parse::<f64>(),
                columns[1].trim().parse::<f32>(),
            ) else {
                continue;
            };
            let charge = columns
                .get(2)
                .and_then(|c| c.trim().parse::<i64>().ok())
                .and_then(|c| u8::try_from(c.unsigned_abs()).ok());
            let score = columns.get(3).and_then(|c| c.trim().parse::<f32>().ok());
            spectrum.peaks.push(DeconvolutedPeak {
                mass,
                intensity,
                charge,
                score,
            });
        }
    }
    Err(MassSpecError::SpectrumNotFound(format!(
        "msalign scan {scan}"
    )))
}

pub fn query_prsm_localize_msalign(
    bytes: &[u8],
    query: &PrsmLocalizeQuery,
) -> Result<PrsmLocalizeResult, MassSpecError> {
    let mut spectra = Vec::with_capacity(query.scans.len());
    for scan in &query.scans {
        let spectrum = parse_msalign_spectrum(bytes, scan)?;
        if spectrum.level != 2 {
            return Err(MassSpecError::InvalidQuery(format!(
                "prsmLocalize needs an MS2 msalign spectrum; scan {scan} is level {}",
                spectrum.level
            )));
        }
        spectra.push((scan.clone(), spectrum));
    }
    query_prsm_localize_spectra(&spectra, query)
}

/// Evaluate localization against already parsed spectra. Batch callers use
/// this to avoid reparsing the same MSAlign scan for every library composition.
pub fn query_prsm_localize_spectra(
    spectra: &[(String, MsAlignSpectrum)],
    query: &PrsmLocalizeQuery,
) -> Result<PrsmLocalizeResult, MassSpecError> {
    query_prsm_localize_spectra_restricted(spectra, query, None)
}

/// Internal variant used when a searched modification has residue-specific
/// chemistry. The public localizer remains whole-sequence by default; the
/// sequence-graph component sweep supplies the same legal positions that its
/// `LibrarySpec` used during DP placement.
pub(crate) fn query_prsm_localize_spectra_restricted(
    spectra: &[(String, MsAlignSpectrum)],
    query: &PrsmLocalizeQuery,
    allowed_positions: Option<&[usize]>,
) -> Result<PrsmLocalizeResult, MassSpecError> {
    let sequence = query.sequence.trim().to_ascii_uppercase();
    let n = sequence.chars().count();
    if n < 3 {
        return Err(MassSpecError::InvalidQuery(
            "prsmLocalize sequence must contain at least three residues".to_owned(),
        ));
    }
    if !query.delta_mass.is_finite() || query.delta_mass == 0.0 {
        return Err(MassSpecError::InvalidQuery(
            "prsmLocalize deltaMass must be a finite non-zero value".to_owned(),
        ));
    }
    let (lo, hi) = match query.range {
        Some((start, end)) => {
            if start == 0 || end < start || end > n {
                return Err(MassSpecError::InvalidQuery(format!(
                    "prsmLocalize range must be 1..{n} inclusive, got {start}-{end}"
                )));
            }
            (start.max(1), end.min(n))
        }
        None => (1, n),
    };
    let candidate_positions: Vec<usize> = match allowed_positions {
        Some(allowed) => {
            let mut positions: Vec<usize> = allowed
                .iter()
                .copied()
                .filter(|position| (lo..=hi).contains(position))
                .collect();
            positions.sort_unstable();
            positions.dedup();
            if positions.is_empty() {
                return Err(MassSpecError::InvalidQuery(
                    "prsmLocalize has no legal positions for this modification".to_owned(),
                ));
            }
            positions
        }
        None => (lo..=hi).collect(),
    };

    if query.scans.is_empty() {
        return Err(MassSpecError::InvalidQuery(
            "prsmLocalize requires at least one scan".to_owned(),
        ));
    }
    if spectra.is_empty() {
        return Err(MassSpecError::InvalidQuery(
            "prsmLocalize requires at least one parsed spectrum".to_owned(),
        ));
    }

    // Pinned modification masses folded into the fragment ladder.
    let mut residue_deltas = vec![0.0_f64; n + 1];
    let mut n_term = 0.0;
    let mut c_term = 0.0;
    for modification in &query.modifications {
        if !modification.delta_mass.is_finite() {
            return Err(MassSpecError::InvalidQuery(
                "prsmLocalize modification deltaMass must be finite".to_owned(),
            ));
        }
        match modification.position.as_str() {
            "nTerm" => n_term += modification.delta_mass,
            "cTerm" => c_term += modification.delta_mass,
            value => {
                let position = value.parse::<usize>().map_err(|_| {
                    MassSpecError::InvalidQuery(format!("unknown modification position {value}"))
                })?;
                if position == 0 || position > n {
                    return Err(MassSpecError::InvalidQuery(format!(
                        "modification position {position} outside 1..{n}"
                    )));
                }
                residue_deltas[position] += modification.delta_mass;
            }
        }
    }
    // extra_before[i] = pinned mass carried by residues 1..=i.
    let mut extra_before = vec![0.0_f64; n + 1];
    for i in 1..=n {
        extra_before[i] = extra_before[i - 1] + residue_deltas[i];
    }

    let mut residue_prefix = vec![0.0_f64; n + 1];
    for (i, residue) in sequence.chars().enumerate() {
        let mass = residue_mass(residue).ok_or_else(|| {
            MassSpecError::InvalidQuery(format!(
                "unknown residue {residue}; B, J, X and Z are not resolvable"
            ))
        })?;
        residue_prefix[i + 1] = residue_prefix[i] + mass;
    }
    let total = residue_prefix[n] + H2O + n_term + c_term + extra_before[n];

    // Fragment ladder: label, neutral base mass (pinned extras included, test
    // mass NOT), cleavage cut (residues 1..cut on the prefix side), side.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    enum Side {
        Prefix,
        Suffix,
    }
    let mut fragments: Vec<(String, f64, usize, Side)> = Vec::with_capacity(4 * n);
    for cut in 1..n {
        let extra_n = extra_before[cut] + n_term;
        fragments.push((
            format!("c{cut}"),
            residue_prefix[cut] + C_OFFSET + extra_n,
            cut,
            Side::Prefix,
        ));
        fragments.push((
            format!("b{cut}"),
            residue_prefix[cut] + extra_n,
            cut,
            Side::Prefix,
        ));
        if cut >= 2 {
            // At cut=1 the suffix is the whole sequence: y_n would pose as the
            // intact precursor, so the C-terminal family starts at y_{n-1}.
            let k = n - cut + 1;
            let base_c = total - residue_prefix[cut - 1] - extra_before[cut - 1] - n_term;
            fragments.push((format!("y{k}"), base_c, cut, Side::Suffix));
            fragments.push((
                format!("z{k}"),
                base_c - H2O + ZDOT_OFFSET,
                cut,
                Side::Suffix,
            ));
            fragments.push((
                format!("z+{k}"),
                base_c - H2O + ZP1_OFFSET,
                cut,
                Side::Suffix,
            ));
        }
    }

    // -- The localization model ------------------------------------------
    //
    // One parameter per candidate residue position p. A fragment that
    // CONTAINS p should be observed at its shifted mass; one that does not
    // should be observed at its plain mass. Observations that contradict the
    // hypothesis cost log10(noise_rate), where the noise rate is MEASURED
    // per spectrum: plain-matched fragments that also carry a shifted match
    // are exactly "a non-carrier matching at the shifted mass".
    //
    //   shifted & contains p    -> 0            (expected)
    //   plain   & !contains p   -> 0            (expected)
    //   shifted & !contains p   -> log10(noise) (noise match)
    //   plain   & contains p    -> log10(noise) (carrier showing plain)
    //   unobserved              -> 0            (excluded)
    //
    // No thresholds anywhere: the site is the argmax, the interval is the
    // flat top, and the reported evidence is the likelihood ratio against
    // the no-shift null. Scans combine by summing per-position
    // log-likelihoods, so an observation gap in one spectrum (flat top) is
    // closed by another scan whose cleavages cover the gap.
    struct Matched {
        label: String,
        plain: Option<usize>,
        shifted: Option<usize>,
        base_mass: f64,
        test_mass: f64,
    }
    let contains = |cut: usize, side: Side, p: usize| -> bool {
        match side {
            Side::Prefix => cut >= p,
            Side::Suffix => cut <= p,
        }
    };
    let theoretical_precursor = total + query.delta_mass;
    let mut scan_summaries = Vec::with_capacity(spectra.len());
    let mut combined_likelihoods = vec![0.0_f64; candidate_positions.len()];
    let mut combined_containing_observable = vec![0usize; candidate_positions.len()];
    let mut combined_containing_shifted = vec![0usize; candidate_positions.len()];
    let mut summed_shifted = 0usize;
    let mut summed_plain_only = 0usize;
    let mut summed_both = 0usize;
    let mut summed_unobserved = 0usize;
    let mut scan_models: Vec<(String, Vec<DeconvolutedPeak>, Vec<(Matched, usize, Side)>)> =
        Vec::with_capacity(spectra.len());
    for (scan_label, spectrum) in spectra {
        let precursor_check = spectrum.precursor_mass.map(|measured| {
            let ppm = (measured - theoretical_precursor) / theoretical_precursor * 1_000_000.0;
            (measured, theoretical_precursor, ppm)
        });
        let mut peaks = spectrum.peaks.clone();
        peaks.sort_by(|left, right| left.mass.total_cmp(&right.mass));
        let base = peaks
            .iter()
            .map(|peak| peak.intensity)
            .fold(0.0_f32, f32::max);
        let floor = INTENSITY_FLOOR as f32 * base;
        let ppm = query.tolerance_ppm;
        let find = |target: f64| -> Option<usize> {
            let tol = target * ppm / 1_000_000.0;
            let mut i = partition_point(&peaks, target - tol);
            let mut best: Option<usize> = None;
            while i < peaks.len() && peaks[i].mass <= target + tol {
                if peaks[i].intensity >= floor
                    && best.is_none_or(|b| peaks[i].intensity > peaks[b].intensity)
                {
                    best = Some(i);
                }
                i += 1;
            }
            best
        };
        // Per fragment precomputed matches: `plain` (at the base mass) and
        // `shifted` (at base + delta). A fragment that matches both is a
        // carrier (the shifted match is the specific one).
        let mut matched: Vec<(Matched, usize, Side)> = Vec::with_capacity(fragments.len());
        for (label, base_mass, cut, side) in &fragments {
            if *base_mass < query.mass_min || *base_mass > query.mass_max {
                continue;
            }
            let test_mass = base_mass + query.delta_mass;
            let shifted = if test_mass >= query.mass_min && test_mass <= query.mass_max {
                find(test_mass)
            } else {
                None
            };
            matched.push((
                Matched {
                    label: label.clone(),
                    plain: find(*base_mass),
                    shifted,
                    base_mass: *base_mass,
                    test_mass,
                },
                *cut,
                *side,
            ));
        }
        let shifted_count = matched
            .iter()
            .filter(|(m, _, _)| m.shifted.is_some())
            .count();
        let plain_only_count = matched
            .iter()
            .filter(|(m, _, _)| m.plain.is_some() && m.shifted.is_none())
            .count();
        let both_count = matched
            .iter()
            .filter(|(m, _, _)| m.plain.is_some() && m.shifted.is_some())
            .count();
        let unobserved_count = matched
            .iter()
            .filter(|(m, _, _)| m.plain.is_none() && m.shifted.is_none())
            .count();
        let noise_rate = ((both_count as f64 + NOISE_SMOOTHING)
            / (plain_only_count as f64 + both_count as f64 + 2.0 * NOISE_SMOOTHING))
            .clamp(1e-6, 1.0 - 1e-6);
        let log10_noise = noise_rate.log10();
        let mut scan_likelihoods = Vec::with_capacity(candidate_positions.len());
        for (offset, p) in candidate_positions.iter().copied().enumerate() {
            let mut log10_l = 0.0;
            for (m, cut, side) in &matched {
                let inside = contains(*cut, *side, p);
                if inside && (m.shifted.is_some() || m.plain.is_some()) {
                    combined_containing_observable[offset] += 1;
                }
                match (m.shifted.is_some(), m.plain.is_some(), inside) {
                    (true, _, true) => {
                        combined_containing_shifted[offset] += 1;
                    }
                    (true, _, false) => log10_l += log10_noise, // noise match
                    (false, true, true) => log10_l += log10_noise, // carrier plain
                    (false, true, false) => {}                  // expected
                    _ => {}                                     // unobserved
                }
            }
            scan_likelihoods.push(log10_l);
        }
        // The no-shift null: every shifted match is noise.
        let null_likelihood = shifted_count as f64 * log10_noise;
        let scan_max = scan_likelihoods
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        for (offset, likelihood) in scan_likelihoods.iter().enumerate() {
            combined_likelihoods[offset] += likelihood;
        }
        summed_shifted += shifted_count;
        summed_plain_only += plain_only_count;
        summed_both += both_count;
        summed_unobserved += unobserved_count;
        scan_summaries.push(ScanLocalizeSummary {
            scan: scan_label.clone(),
            precursor_check,
            shifted: shifted_count,
            plain_only: plain_only_count,
            both: both_count,
            unobserved: unobserved_count,
            noise_rate,
            evidence_log10: (shifted_count > 0).then_some(scan_max - null_likelihood),
        });
        scan_models.push((scan_label.clone(), peaks, matched));
    }
    let likelihoods = combined_likelihoods;
    let shifted_count = summed_shifted;
    let evidence = LocalizeEvidence {
        shifted: summed_shifted,
        plain_only: summed_plain_only,
        both: summed_both,
        unobserved: summed_unobserved,
        noise_rate: ((summed_both as f64 + NOISE_SMOOTHING)
            / (summed_plain_only as f64 + summed_both as f64 + 2.0 * NOISE_SMOOTHING))
            .clamp(1e-6, 1.0 - 1e-6),
    };
    let mut positions: Vec<PositionStats> = Vec::with_capacity(candidate_positions.len());
    let max_likelihood = likelihoods
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    for (offset, p) in candidate_positions.iter().copied().enumerate() {
        positions.push(PositionStats {
            position: p,
            residue: sequence.chars().nth(p - 1).unwrap_or('?'),
            relative_log10_likelihood: likelihoods[offset] - max_likelihood,
            containing_observable: combined_containing_observable[offset],
            containing_shifted: combined_containing_shifted[offset],
        });
    }
    let summed_null: f64 = scan_summaries
        .iter()
        .map(|summary| summary.shifted as f64 * summary.noise_rate.log10())
        .sum();
    let evidence_log10 = (shifted_count > 0).then_some(max_likelihood - summed_null);

    // The called region: the likelihood's flat top (within floating-point
    // equality of the maximum). A unique top is a site call; a flat top is an
    // honestly unresolved interval.
    let mut top: Vec<usize> = positions
        .iter()
        .zip(&likelihoods)
        .filter(|(_, likelihood)| **likelihood >= max_likelihood - 1e-9)
        .map(|(position, _)| position.position)
        .collect();
    top.sort_unstable();
    let mut notes = Vec::new();
    let (site_interval, best_site) = if top.is_empty() || shifted_count == 0 {
        notes.push(
            "no fragment was observed at the shifted mass -- the shift is not \
             localized by this measurement"
                .to_owned(),
        );
        ((0, 0), None)
    } else {
        let lower = top[0];
        let upper = *top.last().unwrap();
        match evidence_log10 {
            Some(value) if value < EVIDENCE_BAR => {
                notes.push(format!(
                    "best position only {value:.1} log10 units above the no-shift null \
                     (below the p<0.01 bar) -- not localized by this measurement"
                ));
                ((0, 0), None)
            }
            _ => {
                if lower == upper {
                    notes.push(format!(
                        "unique likelihood top at residue {lower}; evidence {:+.1} log10 \
                         over the no-shift null at measured noise rate {:.4}",
                        evidence_log10.unwrap_or(0.0),
                        evidence.noise_rate
                    ));
                    ((lower, lower), Some(lower))
                } else {
                    notes.push(format!(
                        "likelihood flat over residues {lower}-{upper} (an observation gap: \
                         every observable fragment either contains the whole interval or \
                         misses it); evidence {:+.1} log10 over the no-shift null at \
                         measured noise rate {:.4}",
                        evidence_log10.unwrap_or(0.0),
                        evidence.noise_rate
                    ));
                    ((lower, upper), None)
                }
            }
        }
    };

    // Rows classified at the called region's best position (or the argmax
    // when no call survived the evidence bar) for manual inspection.
    let rows_position = best_site
        .or_else(|| {
            likelihoods
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.total_cmp(b))
                .map(|(index, _)| candidate_positions[index])
        })
        .unwrap_or(candidate_positions[0]);
    let mut best_rows: Vec<LocalizeRow> = Vec::new();
    for (scan_label, peaks, matched) in &scan_models {
        for (m, cut, side) in matched {
            let inside = contains(*cut, *side, rows_position);
            let (index, class, theoretical) = match (m.shifted, m.plain, inside) {
                (Some(index), _, true) => (index, "span", m.test_mass),
                (Some(index), _, false) => (index, "preShifted", m.test_mass),
                (None, Some(index), true) => (index, "spanPlain", m.base_mass),
                _ => continue,
            };
            let peak = peaks[index];
            best_rows.push(LocalizeRow {
                scan: scan_label.clone(),
                label: m.label.clone(),
                charge: peak.charge,
                class: class.to_owned(),
                theoretical_mass: theoretical,
                observed_mass: peak.mass,
                error_ppm: (peak.mass - theoretical) / theoretical * 1_000_000.0,
                intensity: peak.intensity,
                topfd_score: peak.score,
            });
        }
    }

    let (annotation, annotate_at) = match best_site {
        Some(site) => (format!("[{:+.4}@{site}]", query.delta_mass), Some(site)),
        None if site_interval.0 > 0 => (
            format!(
                "[{:+.4}|{}-{}]",
                query.delta_mass, site_interval.0, site_interval.1
            ),
            None,
        ),
        None => (String::new(), None),
    };
    let mut proforma = String::with_capacity(sequence.len() + annotation.len() + 16);
    if n_term != 0.0 {
        proforma.push_str(&format!("[{n_term:+.12}]-"));
    }
    for (index, residue) in sequence.chars().enumerate() {
        proforma.push(residue);
        if residue_deltas[index + 1] != 0.0 {
            proforma.push_str(&format!("[{:+.12}]", residue_deltas[index + 1]));
        }
        if Some(index + 1) == annotate_at {
            proforma.push_str(&annotation);
        } else if annotate_at.is_none() && !annotation.is_empty() && index + 1 == site_interval.0 {
            // Range annotation opens at the interval start.
            proforma.push_str(&annotation);
        }
    }
    if c_term != 0.0 {
        proforma.push_str(&format!("-[{c_term:+.12}]"));
    }

    Ok(PrsmLocalizeResult {
        sequence,
        delta_mass: query.delta_mass,
        scans: scan_summaries,
        evidence,
        positions,
        site_interval,
        best_site,
        evidence_log10,
        proforma,
        best_rows,
        notes,
    })
}

fn partition_point(peaks: &[DeconvolutedPeak], target: f64) -> usize {
    let mut lo = 0usize;
    let mut hi = peaks.len();
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if peaks[mid].mass < target {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}

#[cfg(test)]
mod tests {
    use super::*;

    fn query() -> PrsmLocalizeQuery {
        PrsmLocalizeQuery {
            scans: vec!["1".to_owned()],
            sequence: "AAAAAYAAAAK".to_owned(), // Y6 carries the shift
            delta_mass: 77.926_9,
            modifications: Vec::new(),
            tolerance_ppm: 8.0,
            mass_min: 50.0,
            mass_max: 100_000.0,
            range: None,
        }
    }

    /// msalign text with every fragment spanning Y6 placed at its shifted
    /// mass, fragments below the site at their plain mass.
    fn synthetic_msalign() -> String {
        let sequence = "AAAAAYAAAAK";
        let n = sequence.chars().count();
        let mut prefix = vec![0.0_f64; n + 1];
        for (i, residue) in sequence.chars().enumerate() {
            prefix[i + 1] = prefix[i] + residue_mass(residue).unwrap();
        }
        let total = prefix[n] + H2O;
        let mut lines: Vec<String> = vec![
            "BEGIN IONS".into(),
            "SCANS=1".into(),
            "LEVEL=2".into(),
            "ACTIVATION=ETD".into(),
            format!("PRECURSOR_MASS={:.6}", total + 77.926_9),
            "PRECURSOR_CHARGE=3".into(),
        ];
        let peak = |lines: &mut Vec<String>, mass: f64, charge: u8| {
            lines.push(format!("{mass:.6}\t1000.0\t{charge}\t1.00"));
        };
        for cut in 1..n {
            let carries = cut >= 6;
            peak(
                &mut lines,
                prefix[cut] + if carries { 77.926_9 } else { 0.0 },
                2,
            );
            peak(
                &mut lines,
                prefix[cut] + C_OFFSET + if carries { 77.926_9 } else { 0.0 },
                3,
            );
            if cut >= 2 {
                let carries_s = cut <= 6;
                // total - prefix[cut-1] is the y-ion mass (suffix + H2O).
                peak(
                    &mut lines,
                    total - prefix[cut - 1] + if carries_s { 77.926_9 } else { 0.0 },
                    4,
                );
            }
        }
        lines.push("END IONS".into());
        lines.join("\n")
    }

    #[test]
    fn sweep_pins_the_site_that_carries_the_shift() {
        let result = query_prsm_localize_msalign(synthetic_msalign().as_bytes(), &query()).unwrap();
        assert_eq!(result.best_site, Some(6), "site must be searched out as Y6");
        assert!(result.proforma.contains("@6"), "{}", result.proforma);
        // The precursor check closes to zero ppm on this self-consistent input.
        let (_, _, ppm) = result.scans[0].precursor_check.unwrap();
        assert!(ppm.abs() < 1.0, "precursor ppm {ppm}");
        // The likelihood top is unique and far above the no-shift null.
        assert!(result.evidence_log10.unwrap_or(0.0) > 2.0);
        // Neighbouring positions sit strictly below the top.
        let at = |p: usize| {
            result
                .positions
                .iter()
                .find(|position| position.position == p)
                .unwrap()
                .relative_log10_likelihood
        };
        assert!(at(5) < -1e-9 && at(7) < -1e-9);
    }

    #[test]
    fn missing_scan_is_an_honest_error() {
        let error = query_prsm_localize_msalign(
            synthetic_msalign().as_bytes(),
            &PrsmLocalizeQuery {
                scans: vec!["999".to_owned()],
                ..query()
            },
        )
        .unwrap_err();
        assert!(error.to_string().contains("999"), "{error}");
    }
}
