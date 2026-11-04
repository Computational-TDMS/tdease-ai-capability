/*!
Mass-shift calibration from an MS1 isotope envelope.

This is the **convenience workflow** that composes the crate's three atomic
capabilities — `proforma::parse_proforma` (exact composition), `theory`
(BRAIN theoretical envelope), and `fit_envelope` (fixed-frame fit). Each of
those is independently callable; this module orchestrates them into the
"given a candidate proteoform and an observed precursor, what is the refined
mass shift?" answer.

The three first-principles fixes live here:
- **P0-a exact composition**: when a `proforma` is supplied the composition
  comes from real residues + exact modification net formulas (incl.
  substitutions like `Bromo = Br - H`), not from averagine.
- **P0-b anchored coordinates**: `theoretical` peaks are the *anchored* BRAIN
  envelope (absolute m/z with M0 at the candidate's exact mass) — the same
  coordinate frame the fit uses, so a downstream plot aligns observed vs
  theoretical.
- **P0-c interface honesty**: declared base + prior positions extraction; the fit's
  match window is an explicit `match_window_ppm` parameter, never silently
  overridden by the extraction tolerance.
*/

use serde::{Deserialize, Serialize};

use crate::constants::{neutral_from_mz, NEUTRON_SHIFT, PROTON};
use crate::envelope::{aggregate_window_inner, IsotopeEnvelopeQuery};
use crate::error::MassSpecError;
use crate::isotopic::TheoreticalPeakPoint;
use crate::proforma::{
    apply_net_formula, chemical_from_counts, parse_proforma, ExternalModification,
};
use crate::theory::{anchored_envelope, averagine_composition, theoretical_envelope};
use crate::Peak;

/// Input for calibrating a mass shift from an MS1 isotope envelope.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RefineMassShiftQuery {
    /// Exact baseline proteoform, e.g. `"[Acetyl]-SGRGK..."`. When a separate
    /// `modification_formula` is supplied it is added to this composition;
    /// otherwise this is treated as the complete candidate proteoform.
    pub proforma: Option<String>,
    /// Named-modification injection for the ProForma parser (name -> net formula).
    #[serde(default)]
    pub external_mods: Vec<ExternalModification>,
    /// The monoisotopic neutral mass of the *unmodified* sequence. Used when
    /// `proforma` is None (legacy path; composition is averagine-approximate).
    #[serde(default)]
    pub sequence_mass: Option<f64>,
    /// The neutral mass contribution of the candidate modification
    /// (monoisotopic). Provide one of this or `modification_formula`.
    pub modification_mass: Option<f64>,
    /// The elemental formula of the candidate modification, e.g. "PO3H".
    pub modification_formula: Option<String>,
    /// Measured precursor m/z, which may correspond to a later isotope.
    /// Extraction uses the parsed declared base + prior_shift, not this isotope.
    pub mz: f64,
    /// Charge state of the ion (must not be zero).
    pub charge: i32,
    /// Retention-time window (minutes) to aggregate MS1 over.
    pub rt_min: f64,
    pub rt_max: f64,
    /// Aggregation tolerance (ppm) for the extraction bins.
    pub tolerance_ppm: f64,
    /// The prior mass shift reported by the search engine, if any.
    #[serde(default)]
    pub prior_shift: f64,
    /// The fit's per-peak match window (ppm), default 30, floored at 0.015 Da.
    /// Independent of `tolerance_ppm` (P0-c).
    #[serde(default = "default_match_window_ppm")]
    pub match_window_ppm: f64,
}

fn default_match_window_ppm() -> f64 {
    30.0
}

/// Fixed isotope-index neighbourhood used only by the ingest calibration
/// workflow. Interactive candidate fitting is always anchored at offset 0.
const CALIBRATION_MAX_OFFSET: i32 = 2;

/// One fitted candidate evaluation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RefinedShiftEvaluation {
    /// Fixed C13 isotope-index correction selected during ingest calibration.
    pub offset: i32,
    /// Number of observed peaks matched to theoretical peaks.
    pub matched_peaks: usize,
    /// Sum of the theoretical intensities of the matched peaks (normalized,
    /// ~0..1). Weighting matches by the *theoretical* intensity stops a weak
    /// low-abundance isotopologue from swinging the offset selection.
    pub coverage: f64,
    /// Pearson correlation between observed and theoretical intensity shapes.
    /// Shape QC only — NOT an identity probability.
    pub pearson: Option<f64>,
    /// The refined monoisotopic neutral mass.
    pub refined_mass: f64,
    /// The shared sub-Da offset absorbed from the weighted mean of per-peak residuals.
    pub fitted_da: f64,
    /// Per-peak residuals `observed - theoretical_center_mass` in Da.
    pub residuals: Vec<f64>,
}

/// The result of a refinement pass.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RefineMassShiftResult {
    /// The refined monoisotopic neutral mass.
    pub refined_mass: f64,
    /// The refined mass shift in the **bare-sequence frame**
    /// (refined mass - unmodified sequence mass): every declared
    /// modification — including terminal backgrounds such as N-terminal
    /// acetylation — is INSIDE this number. Kept for compatibility; prefer
    /// `refined_shift_after_declared` when comparing against a search
    /// engine's residual.
    pub refined_shift: f64,
    /// Sum of the proforma's declared modification deltas (0 on the legacy
    /// averagine path). Tool-derived; consumers never hand-compute the
    /// background.
    #[serde(default)]
    pub declared_delta_mass: f64,
    /// The refined shift in the **post-declared-forms frame**
    /// (refined mass - bare sequence - declared deltas) — directly
    /// comparable to TopPIC's `Unexpected modifications` residual, which is
    /// also measured after the assigned fixed forms.
    #[serde(default)]
    pub refined_shift_after_declared: f64,
    /// The selected C13 isotope-index correction.
    pub offset: i32,
    /// Number of matched peaks at the best offset.
    pub matched_peaks: usize,
    /// Pearson correlation at the best offset (shape QC).
    pub pearson: Option<f64>,
    /// The theoretical envelope used — **anchored** absolute m/z, the same
    /// coordinate frame as `observed` (P0-b).
    pub theoretical: Vec<TheoreticalPeakPoint>,
    /// The aggregated observed envelope peaks.
    pub observed: Vec<TheoreticalPeakPoint>,
    /// Per-offset evaluations (all candidates, not just best).
    pub evaluations: Vec<RefinedShiftEvaluation>,
    /// Scan count aggregated over.
    pub scan_count: usize,
    /// Warnings (e.g. mass-only path used averagine composition).
    pub warnings: Vec<String>,
    /// Boundary adjudication: when the top two offset candidates are not
    /// separated (equal matched peaks and fitted_da within 0.01 Da), this is
    /// set instead of the best being silently trusted. `None` means the best
    /// was separated from the runner-up. The floor rule: where the evidence
    /// does not decide, the output says so.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ambiguity: Option<String>,
}

/// Refine a mass shift using the MS1 isotope envelope of the candidate.
///
/// This is the convenience workflow composing the three atomic capabilities.
/// `refine_mass_shift` is the public entry point; a `calibrate_mass_shift`
/// alias is also provided for callers who prefer the workflow name.
pub fn refine_mass_shift(
    bytes: &[u8],
    query: &RefineMassShiftQuery,
) -> Result<RefineMassShiftResult, MassSpecError> {
    calibrate_mass_shift(bytes, query)
}

/// Convenience alias matching the workflow's intended name.
pub fn calibrate_mass_shift(
    bytes: &[u8],
    query: &RefineMassShiftQuery,
) -> Result<RefineMassShiftResult, MassSpecError> {
    validate(query)?;
    let (composition, anchor_mass, bare_mass, declared_base, declared_delta, mut warnings) =
        candidate_chemistry(query)?;
    let (theoretical, anchored) = theoretical_for(query, &composition, anchor_mass);
    let (expected_center, heal) = resolve_expected_center(anchor_mass, declared_base, query.prior_shift);
    if let Some(note) = heal {
        warnings.push(note);
    }
    let (observed_peaks, scan_count) = extract_observed(bytes, query, expected_center)?;
    finish_calibrate(
        query,
        observed_peaks,
        scan_count,
        &theoretical,
        &anchored,
        bare_mass,
        declared_base,
        declared_delta,
        expected_center,
        warnings,
    )
}

/// Calibrate a mass shift from a **pre-extracted** observed envelope.
///
/// This is the `--observed` reuse path: the caller aggregates the MS1 envelope
/// once (its anchor depends only on `mz`/`prior_shift`, not the hypothesis) and
/// then evaluates N candidate hypotheses without re-reading the mzML. Shares
/// the exact same composition/theory/fit pipeline as [`calibrate_mass_shift`]
/// (D5 — one algorithm, two entry points).
pub fn calibrate_mass_shift_from_observed(
    observed: &[Peak],
    scan_count: usize,
    query: &RefineMassShiftQuery,
) -> Result<RefineMassShiftResult, MassSpecError> {
    validate(query)?;
    let (composition, anchor_mass, bare_mass, declared_base, declared_delta, mut warnings) =
        candidate_chemistry(query)?;
    let (theoretical, anchored) = theoretical_for(query, &composition, anchor_mass);
    let (expected_center, heal) = resolve_expected_center(anchor_mass, declared_base, query.prior_shift);
    if let Some(note) = heal {
        warnings.push(note);
    }
    finish_calibrate(
        query,
        observed.to_vec(),
        scan_count,
        &theoretical,
        &anchored,
        bare_mass,
        declared_base,
        declared_delta,
        expected_center,
        warnings,
    )
}

/// Calibrate independent shifts while parsing the mzML only once.
pub fn calibrate_mass_shifts(
    bytes: &[u8],
    queries: &[RefineMassShiftQuery],
) -> Vec<Result<RefineMassShiftResult, MassSpecError>> {
    struct Prepared {
        theoretical: Vec<TheoreticalPeakPoint>,
        anchored: crate::theory::AnchoredEnvelope,
        bare_mass: f64,
        declared_base: f64,
        declared_delta: f64,
        expected_center: f64,
        warnings: Vec<String>,
        envelope_query: IsotopeEnvelopeQuery,
        window_center: f64,
    }
    let prepared: Vec<Result<Prepared, MassSpecError>> = queries
        .iter()
        .map(|query| {
            validate(query)?;
            let (composition, anchor_mass, bare_mass, declared_base, declared_delta, mut warnings) =
                candidate_chemistry(query)?;
            let (theoretical, anchored) = theoretical_for(query, &composition, anchor_mass);
            let (expected_center, heal) = resolve_expected_center(anchor_mass, declared_base, query.prior_shift);
            if let Some(note) = heal {
                warnings.push(note);
            }
            let (envelope_query, window_center) = extraction_query(query, expected_center);
            Ok(Prepared {
                theoretical,
                anchored,
                bare_mass,
                declared_base,
                declared_delta,
                expected_center,
                warnings,
                envelope_query,
                window_center,
            })
        })
        .collect();
    let extraction_inputs: Vec<_> = prepared
        .iter()
        .filter_map(|item| item.as_ref().ok())
        .map(|item| (&item.envelope_query, item.window_center))
        .collect();
    let mut extracted =
        crate::envelope::aggregate_windows_inner(bytes, &extraction_inputs).into_iter();
    prepared
        .into_iter()
        .zip(queries)
        .map(|(item, query)| match item {
            Err(error) => Err(error),
            Ok(item) => {
                let (peaks, scan_count) = extracted
                    .next()
                    .expect("one extraction per prepared query")?;
                finish_calibrate(
                    query,
                    peaks,
                    scan_count,
                    &item.theoretical,
                    &item.anchored,
                    item.bare_mass,
                    item.declared_base,
                    item.declared_delta,
                    item.expected_center,
                    item.warnings,
                )
            }
        })
        .collect()
}

/// Resolve the candidate's chemistry. Returns
/// `(composition, anchor_mass, bare_mass, declared_base_mass, declared_delta_mass, warnings)`:
/// - `bare_mass` — the unmodified sequence (all declared deltas stripped);
/// - `declared_base_mass` — bare + declared deltas (the fixed base the
///   residual is measured against);
/// - `declared_delta_mass` — the tool-summed declared deltas (terminal
///   backgrounds included). The refined shift is reported in BOTH frames so
///   no consumer ever hand-computes a background.
fn candidate_chemistry(
    query: &RefineMassShiftQuery,
) -> Result<
    (
        chemical_elements::ChemicalComposition<'static>,
        f64,
        f64,
        f64,
        f64,
        Vec<String>,
    ),
    MassSpecError,
> {
    // Prefer a real ProFormA; fall back to the legacy sequence_mass path.
    if let Some(proforma) = &query.proforma {
        let parsed = parse_proforma(proforma, &query.external_mods)?;
        let mut composition = chemical_from_counts(&parsed.composition);
        let declared_delta: f64 = parsed.modifications.iter().map(|m| m.delta_mass).sum();
        let declared_base = parsed.monoisotopic_mass;
        let bare_mass = declared_base - declared_delta;
        if let Some(formula) = &query.modification_formula {
            apply_net_formula(&mut composition, formula)?;
            let candidate_mass = query.modification_mass.map_or_else(
                || composition.mass(),
                |delta| parsed.monoisotopic_mass + delta,
            );
            return Ok((
                composition,
                candidate_mass,
                bare_mass,
                declared_base,
                declared_delta,
                parsed.warnings,
            ));
        }
        if let Some(residual) = query.modification_mass {
            // Declared named tags are the FIXED BASE — the same doctrine as
            // resolve-shift's fixed forms — and `modification_mass` is the
            // UNKNOWN residual delta measured in the post-declared frame
            // (TopPIC's unexpected-shift frame). No identity is claimed for
            // the residual, so this is not a floor-rule double source. The
            // protein's own isotope pattern dominates a small-residual
            // envelope, so the parsed composition (declared mods included)
            // drives the shape.
            let anchor_mass = parsed.monoisotopic_mass + residual;
            let mut warnings = parsed.warnings;
            warnings.push(
                "mass-only residual: isotope shape uses the declared baseline composition; calibration is not composition or identity evidence"
                    .to_owned(),
            );
            return Ok((
                composition,
                anchor_mass,
                bare_mass,
                declared_base,
                declared_delta,
                warnings,
            ));
        }
        // Backward-compatible complete-candidate path: the proforma itself
        // closes the model; anchor at its own mass.
        Ok((
            composition,
            parsed.monoisotopic_mass,
            bare_mass,
            declared_base,
            declared_delta,
            parsed.warnings,
        ))
    } else {
        let sequence_mass = query.sequence_mass.ok_or_else(|| {
            MassSpecError::InvalidQuery(
                "refineMassShift requires proforma or sequence_mass".to_owned(),
            )
        })?;
        let sequence_comp = averagine_composition(sequence_mass);
        let mod_comp = match &query.modification_formula {
            Some(formula) => {
                let mut c = chemical_elements::ChemicalComposition::new();
                apply_formula(&mut c, formula)?;
                Some(c)
            }
            None => {
                if let Some(mass) = query.modification_mass {
                    Some(averagine_composition(mass))
                } else {
                    return Err(MassSpecError::InvalidQuery(
                        "refineMassShift requires modification_formula or modification_mass"
                            .to_owned(),
                    ));
                }
            }
        };
        let mut composition = sequence_comp;
        if let Some(mc) = &mod_comp {
            add_composition(&mut composition, mc);
        }
        // Anchor: prefer the caller's exact modification_mass (the
        // UniMod-style delta). A formula alone is NOT the delta for a
        // substitution (Bromo = 77.9105, but "Br" = 78.9183).
        let mod_mono = query
            .modification_mass
            .unwrap_or_else(|| mod_comp.as_ref().map(|c| c.mass()).unwrap_or(0.0));
        let anchor_mass = sequence_mass + mod_mono;
        let warnings = vec![
            "legacy path used averagine composition; provide a proforma for exact chemistry"
                .to_owned(),
        ];
        Ok((
            composition,
            anchor_mass,
            sequence_mass,
            sequence_mass,
            0.0,
            warnings,
        ))
    }
}

/// P0-b: the anchored BRAIN theoretical envelope and its absolute-m/z peaks.
fn theoretical_for(
    query: &RefineMassShiftQuery,
    composition: &chemical_elements::ChemicalComposition<'static>,
    anchor_mass: f64,
) -> (Vec<TheoreticalPeakPoint>, crate::theory::AnchoredEnvelope) {
    let env = theoretical_envelope(composition, query.charge, 0.99, 0.001);
    let anchored = anchored_envelope(&env, anchor_mass);
    let theoretical: Vec<TheoreticalPeakPoint> = anchored
        .peaks
        .iter()
        .map(|p| TheoreticalPeakPoint {
            mz: p.mz,
            intensity: p.intensity,
        })
        .collect();
    (theoretical, anchored)
}

/// One extraction basis for single and batch calibration. The window centers
/// on the resolved expected position (see `resolve_expected_center`): the
/// engine prior when it is coherent with the modification anchor, the anchor
/// itself when the caller left the prior at zero while declaring a residual.
/// The raw selected m/z may be a later isotope, so it cannot replace this M0
/// prior.
fn extraction_query(
    query: &RefineMassShiftQuery,
    expected_center: f64,
) -> (IsotopeEnvelopeQuery, f64) {
    let window_center = expected_center;
    let envelope_query = IsotopeEnvelopeQuery {
        mass: Some(window_center),
        proforma: None,
        external_mods: Vec::new(),
        charge: query.charge as u8,
        rt_min: query.rt_min,
        rt_max: query.rt_max,
        tolerance_ppm: query.tolerance_ppm,
        isotopes: 0,
        max_offset: CALIBRATION_MAX_OFFSET,
        theory_locked: true,
    };
    (envelope_query, window_center)
}

/// Resolve where the envelope is expected: the engine prior
/// (`declared_base + prior_shift`) when one was given, the modification
/// anchor itself when the prior is absent. A caller that declares a residual
/// via `modification_mass` while leaving `prior_shift` at the zero default is
/// internally inconsistent — the theoretical comb is anchored at the
/// residual, so a window at the bare base collects peaks the comb can never
/// match. Healing to the anchor keeps such calls working and says so. A
/// NON-zero prior that disagrees with the anchor is NOT healed: multi-
/// candidate scans deliberately point the prior at the observed feature
/// while anchoring theories at different candidate masses, bridged by the
/// isotope-offset grid.
fn resolve_expected_center(
    anchor_mass: f64,
    declared_base: f64,
    prior_shift: f64,
) -> (f64, Option<String>) {
    if prior_shift == 0.0 {
        if (anchor_mass - declared_base).abs() > 0.5 {
            return (
                anchor_mass,
                Some(format!(
                    "no prior_shift given while the modification anchor sits \
                     {:.4} Da from the declared base; using the anchor as the \
                     expected center — pass the search engine's reported shift \
                     as prior-shift when available",
                    anchor_mass - declared_base
                )),
            );
        }
        return (declared_base, None);
    }
    (declared_base + prior_shift, None)
}

/// Aggregate MS1 around the resolved expected position.
fn extract_observed(
    bytes: &[u8],
    query: &RefineMassShiftQuery,
    expected_center: f64,
) -> Result<(Vec<Peak>, usize), MassSpecError> {
    let (envelope_query, window_center) = extraction_query(query, expected_center);
    aggregate_window_inner(bytes, &envelope_query, window_center)
}

/// The shared workflow tail: per-offset fit with the residual model, the
/// boundary adjudication, and the result.
fn finish_calibrate(
    query: &RefineMassShiftQuery,
    observed_peaks: Vec<Peak>,
    scan_count: usize,
    theoretical: &[TheoreticalPeakPoint],
    anchored: &crate::theory::AnchoredEnvelope,
    bare_mass: f64,
    declared_base: f64,
    declared_delta: f64,
    expected_center: f64,
    mut warnings: Vec<String>,
) -> Result<RefineMassShiftResult, MassSpecError> {
    if observed_peaks.is_empty() {
        return Err(MassSpecError::InvalidQuery(
            "calibrateMassShift: no MS1 peaks matched in the aggregation window".to_owned(),
        ));
    }
    let observed: Vec<TheoreticalPeakPoint> = observed_peaks
        .iter()
        .map(|p| TheoreticalPeakPoint {
            mz: p.mz,
            intensity: f64::from(p.intensity),
        })
        .collect();

    // ---- fit: per-offset evaluation with residual model --------------------
    // The scoring prior sits on the same resolved expected center the
    // extraction used. Anchoring it on the bare declared base instead made
    // every neutron step toward that base worth ~100 score points — enough
    // to drag the offset selection off the correct registration whenever a
    // residual is in play.
    let prior_center = expected_center;
    let mut evaluations = Vec::new();
    let mut best_score = f64::NEG_INFINITY;
    let mut best_index = 0usize;
    // Calibration may inspect a bounded C13 index neighbourhood, but this is
    // fixed scientific correction logic, never a caller-controlled fit knob.
    for offset in -2..=2 {
        let eval = evaluate_calibration_offset(
            &observed_peaks,
            anchored,
            query.charge,
            offset,
            query.match_window_ppm,
        );
        // Score: residual model dominates (P0-b), shape is QC, prior keeps the
        // refinement a correction rather than a whole-neutron re-assignment.
        let score = eval.matched_peaks as f64 * 5.0
            + eval.coverage * 10.0
            + eval.pearson.unwrap_or(-1.0) * 5.0
            - (eval.refined_mass - prior_center).abs() * 100.0
            + (1.0 - (offset as f64).abs() * 0.001);
        if score > best_score {
            best_score = score;
            best_index = evaluations.len();
        }
        evaluations.push(eval);
    }
    let best_eval = evaluations[best_index].clone();
    if best_eval.matched_peaks == 0 {
        warnings.push(
            "calibrateMassShift: no peaks matched at any offset; result is not trustworthy"
                .to_owned(),
        );
    }
    let ambiguity = runner_up_ambiguity(&evaluations, best_index);

    Ok(RefineMassShiftResult {
        refined_mass: best_eval.refined_mass,
        refined_shift: best_eval.refined_mass - bare_mass,
        declared_delta_mass: declared_delta,
        refined_shift_after_declared: best_eval.refined_mass - declared_base,
        offset: best_eval.offset,
        matched_peaks: best_eval.matched_peaks,
        pearson: best_eval.pearson,
        theoretical: theoretical.to_vec(),
        observed,
        evaluations,
        scan_count,
        warnings,
        ambiguity,
    })
}

/// Boundary adjudication: when a non-best candidate matches as many peaks and
/// its refined mass comes within **half a neutron** (0.5 Da) of the best, the
/// envelope cannot pick an offset on its own — the two are not separated, and
/// the result says so (the floor rule: report it, don't silently trust the top
/// rank).
///
/// Note on `fitted_da`: it is the shared sub-Da calibration offset and is
/// nearly constant across integer offsets (empirically 0.0376 vs 0.0385 Da), so
/// it is *not* a separator. The refined masses, which include the integer
/// shift, differ by a whole neutron (1.003 Da) between adjacent offsets — that
/// is the quantity that separates them.
fn runner_up_ambiguity(
    evaluations: &[RefinedShiftEvaluation],
    best_index: usize,
) -> Option<String> {
    let best = &evaluations[best_index];
    let runner = evaluations
        .iter()
        .enumerate()
        .filter(|(index, _)| *index != best_index)
        .map(|(_, evaluation)| evaluation)
        .max_by(|a, b| {
            a.matched_peaks
                .cmp(&b.matched_peaks)
                .then_with(|| a.coverage.total_cmp(&b.coverage))
                .then_with(|| {
                    a.pearson
                        .unwrap_or(-1.0)
                        .total_cmp(&b.pearson.unwrap_or(-1.0))
                })
        });
    let runner = runner?;
    let mass_delta = (runner.refined_mass - best.refined_mass).abs();
    if runner.matched_peaks == best.matched_peaks && mass_delta < 0.5 {
        Some(format!(
            "not separated: offset {} and offset {} refine to within {:.3} Da; \
             the envelope cannot pick an offset without the prior",
            runner.offset, best.offset, mass_delta
        ))
    } else {
        None
    }
}

fn validate(query: &RefineMassShiftQuery) -> Result<(), MassSpecError> {
    if query.charge <= 0 || query.charge > u8::MAX as i32 {
        return Err(MassSpecError::InvalidQuery(
            "refineMassShift charge must be positive and fit in an 8-bit charge state".to_owned(),
        ));
    }
    if query.proforma.is_none() && query.sequence_mass.is_none() {
        return Err(MassSpecError::InvalidQuery(
            "refineMassShift requires proforma or sequence_mass".to_owned(),
        ));
    }
    if !query.mz.is_finite() || query.mz <= 0.0 {
        return Err(MassSpecError::InvalidQuery(
            "refineMassShift mz must be positive".to_owned(),
        ));
    }
    if query.proforma.is_none()
        && query.modification_formula.is_none()
        && query.modification_mass.is_none()
    {
        return Err(MassSpecError::InvalidQuery(
            "refineMassShift requires modification_formula or modification_mass".to_owned(),
        ));
    }
    // candidate_chemistry treats numeric residuals as positional calibration,
    // including residuals on a modified baseline; it does not assign chemistry.
    if !query.rt_min.is_finite() || !query.rt_max.is_finite() || query.rt_min > query.rt_max {
        return Err(MassSpecError::InvalidQuery(
            "refineMassShift RT window must be ordered and finite".to_owned(),
        ));
    }
    if !query.tolerance_ppm.is_finite() || query.tolerance_ppm <= 0.0 {
        return Err(MassSpecError::InvalidQuery(
            "refineMassShift tolerance_ppm must be positive".to_owned(),
        ));
    }
    Ok(())
}

struct EnvelopeMatches<'a> {
    observed: Vec<&'a Peak>,
    theoretical: Vec<usize>,
    residuals: Vec<f64>,
}

pub const DEFAULT_ENVELOPE_WINDOW_MZ: f64 = 0.015;

/// Match one anchored theoretical envelope to observed peaks. `shift` is used
/// only by the ingest calibration workflow; interactive `fit_envelope` passes
/// zero and never exposes a shift control.
fn match_envelope<'a>(
    observed: &'a [Peak],
    anchored: &crate::theory::AnchoredEnvelope,
    charge: i32,
    shift: f64,
    match_window_ppm: f64,
    minimum_window_mz: f64,
) -> EnvelopeMatches<'a> {
    let mut matched_obs = Vec::new();
    let mut matched_theo = Vec::new();
    let mut residuals = Vec::new();
    // Each observed peak may be matched by at most one theoretical peak.
    // `find` (the old behaviour) returned the first in-m/z-order peak within
    // the window and never marked a peak as taken, so a dense or co-eluting
    // window could pair two theoretical peaks with the same observed peak.
    // The nearest unused peak is used instead (TDEaseNext align_envelope).
    let mut used: Vec<bool> = vec![false; observed.len()];
    for (i, cm) in anchored.center_masses.iter().enumerate() {
        let neutral = cm + shift;
        let target = (neutral + (charge as f64 * PROTON)) / charge.abs() as f64;
        // Match window: explicit comparison tolerance, floored so small
        // measurement error does not silently drop low-isotopologue peaks.
        let window = (target * match_window_ppm / 1_000_000.0).max(minimum_window_mz);
        let mut best: Option<usize> = None;
        let mut best_delta = f64::INFINITY;
        for (j, peak) in observed.iter().enumerate() {
            if used[j] {
                continue;
            }
            let delta = (peak.mz - target).abs();
            if delta <= window && delta < best_delta {
                best_delta = delta;
                best = Some(j);
            }
        }
        if let Some(j) = best {
            used[j] = true;
            let p = &observed[j];
            let obs_neutral = neutral_from_mz(p.mz, charge, PROTON);
            residuals.push(obs_neutral - neutral);
            matched_obs.push(p);
            matched_theo.push(i);
        }
    }
    EnvelopeMatches {
        observed: matched_obs,
        theoretical: matched_theo,
        residuals,
    }
}

/// Evaluate one integer isotope-index offset for ingest calibration.
fn evaluate_calibration_offset(
    observed: &[Peak],
    anchored: &crate::theory::AnchoredEnvelope,
    charge: i32,
    offset: i32,
    match_window_ppm: f64,
) -> RefinedShiftEvaluation {
    let shift = offset as f64 * NEUTRON_SHIFT;
    let matches = match_envelope(
        observed,
        anchored,
        charge,
        shift,
        match_window_ppm,
        DEFAULT_ENVELOPE_WINDOW_MZ,
    );
    // Weighted mean of residuals by observed intensity (the calibration offset).
    let (fitted_da, total_w) = matches.observed.iter().zip(matches.residuals.iter()).fold(
        (0.0, 0.0),
        |(acc, w), (p, r)| {
            let wi = f64::from(p.intensity);
            (acc + wi * r, w + wi)
        },
    );
    let fitted_da = if total_w > 0.0 {
        fitted_da / total_w
    } else {
        0.0
    };

    let refined_mass = anchored.center_masses.first().copied().unwrap_or(0.0) + shift + fitted_da;

    let pearson = pearson_match(&matches.observed, &matches.theoretical, anchored);
    let coverage = matches
        .theoretical
        .iter()
        .map(|i| anchored.peaks[*i].intensity)
        .sum::<f64>();
    RefinedShiftEvaluation {
        offset,
        matched_peaks: matches.observed.len(),
        coverage,
        pearson,
        refined_mass,
        fitted_da,
        residuals: matches.residuals,
    }
}

fn pearson_match(
    observed: &[&Peak],
    theo_indices: &[usize],
    anchored: &crate::theory::AnchoredEnvelope,
) -> Option<f64> {
    if observed.len() < 2 {
        return None;
    }
    let obs: Vec<f64> = observed.iter().map(|p| f64::from(p.intensity)).collect();
    let theo: Vec<f64> = theo_indices
        .iter()
        .map(|i| anchored.peaks[*i].intensity)
        .collect();
    crate::statistics::pearson(&obs, &theo)
}

/// Add a modification formula into a composition (signed net formulas allowed).
/// Routes through the ProForma net-formula parser so `Br-H`, `O-N` etc. work —
/// `chemical_elements::ChemicalComposition::parse` cannot handle negative
/// counts.
fn apply_formula(
    composition: &mut chemical_elements::ChemicalComposition<'_>,
    formula: &str,
) -> Result<(), MassSpecError> {
    apply_net_formula(composition, formula)
}

/// Add `add` into `base` element-wise, in place.
fn add_composition(
    base: &mut chemical_elements::ChemicalComposition<'_>,
    add: &chemical_elements::ChemicalComposition<'_>,
) {
    for (elt, count) in add.iter() {
        base.inc_str(&elt.element.symbol, *count);
    }
}

// -------------------------------------------------------------------------
// Atomic capability: fitEnvelope
//
// Observed + theoretical (anchored) -> one fixed-zero fit.
// Pure data flow, independent of mzML extraction: the caller supplies the
// observed peaks and an already-anchored theoretical envelope. Integer C13
// registration belongs only to ingest calibration above, never this operation.
// -------------------------------------------------------------------------

/// The fit's match window (ppm), floored so small measurement error does not
/// silently drop low-isotopologue peaks. This is a comparison tolerance, not a
/// mass translation.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct FitEnvelopeQuery {
    #[serde(default = "default_match_window_ppm")]
    pub match_window_ppm: f64,
    /// Absolute m/z floor; zero requests a strict ppm-only window.
    #[serde(default = "default_minimum_window_mz")]
    pub minimum_window_mz: f64,
}

fn default_minimum_window_mz() -> f64 {
    DEFAULT_ENVELOPE_WINDOW_MZ
}

/// One fixed-frame MS1 comparison. It deliberately has no isotope-offset or
/// fitted-mass field: candidate models differ by their complete chemistry.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct FitEnvelopeEvaluation {
    /// Number of observed peaks matched to the anchored theoretical envelope.
    pub matched_peaks: usize,
    /// Sum of the theoretical intensities of the matched peaks (~0..1).
    pub coverage: f64,
    /// Correlation of observed and theoretical intensities for matched peaks.
    /// Shape QC only — not an identity probability.
    pub pearson: Option<f64>,
    /// Per-peak residuals `observed - theoretical_center_mass` in Da.
    pub residuals: Vec<f64>,
}

/// Fit an observed envelope against an anchored theoretical envelope.
///
/// Returns one fixed-frame evaluation. Candidate identity comparisons must
/// differ in their complete ProForma, not in a hidden mass offset.
pub fn fit_envelope(
    observed: &[Peak],
    theoretical: &crate::theory::AnchoredEnvelope,
    query: &FitEnvelopeQuery,
) -> Result<FitEnvelopeEvaluation, MassSpecError> {
    validate_envelope_fit(observed, theoretical, query)?;
    let matches = match_envelope(
        observed,
        theoretical,
        theoretical.charge,
        0.0,
        query.match_window_ppm,
        query.minimum_window_mz,
    );
    let coverage = matches
        .theoretical
        .iter()
        .map(|index| theoretical.peaks[*index].intensity)
        .sum();
    Ok(FitEnvelopeEvaluation {
        matched_peaks: matches.observed.len(),
        coverage,
        pearson: pearson_match(&matches.observed, &matches.theoretical, theoretical),
        residuals: matches.residuals,
    })
}

pub fn validate_envelope_fit(
    observed: &[Peak],
    theoretical: &crate::theory::AnchoredEnvelope,
    query: &FitEnvelopeQuery,
) -> Result<(), MassSpecError> {
    if !query.match_window_ppm.is_finite()
        || query.match_window_ppm <= 0.0
        || theoretical.charge <= 0
        || theoretical.peaks.len() != theoretical.center_masses.len()
    {
        return Err(MassSpecError::InvalidQuery(
            "fitEnvelope requires positive finite tolerance/charge and aligned theoretical arrays"
                .into(),
        ));
    }
    if theoretical.center_masses.is_empty() {
        return Err(MassSpecError::InvalidQuery(
            "fitEnvelope: theoretical envelope is empty".to_owned(),
        ));
    }
    if !query.minimum_window_mz.is_finite()
        || query.minimum_window_mz < 0.0
        || observed.iter().any(|p| {
            !p.mz.is_finite() || p.mz <= 0.0 || !p.intensity.is_finite() || p.intensity < 0.0
        })
        || theoretical
            .center_masses
            .iter()
            .any(|mass| !mass.is_finite() || *mass <= 0.0)
        || theoretical
            .peaks
            .iter()
            .any(|p| !p.intensity.is_finite() || p.intensity < 0.0)
    {
        return Err(MassSpecError::InvalidQuery("fitEnvelope requires finite valid coordinates, intensities and a non-negative minimumWindowMz".into()));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::fs;

    use chemical_elements::isotopic_pattern::baffling::NumPeaksSpec;
    use chemical_elements::isotopic_pattern::BafflingRecursiveIsotopicPatternGenerator;
    use chemical_elements::ChemicalComposition;

    use super::*;

    fn run_refine(xml: &[u8], query: &RefineMassShiftQuery) -> RefineMassShiftResult {
        refine_mass_shift(xml, query).expect("refine should succeed")
    }

    fn brain_ints(formula: &str) -> Vec<f64> {
        let comp = ChemicalComposition::parse(formula).unwrap();
        let mut gen = BafflingRecursiveIsotopicPatternGenerator::new();
        let peaks = gen.isotopic_variants(comp, NumPeaksSpec::Guess, 1, PROTON);
        peaks.iter().map(|p| p.intensity).collect()
    }

    const REAL_FIXTURE: &str = "fixtures/chlamy-h4-scan4113-z11.mzML";

    #[test]
    fn batch_refine_matches_individual_results() {
        let xml = fs::read(REAL_FIXTURE).unwrap();
        let query = RefineMassShiftQuery {
            proforma: None,
            external_mods: Vec::new(),
            sequence_mass: Some(11319.3476),
            modification_formula: None,
            modification_mass: Some(91.9674),
            mz: (11411.315 + 11.0 * PROTON) / 11.0,
            charge: 11,
            rt_min: -1.0,
            rt_max: 1.0,
            tolerance_ppm: 100.0,
            prior_shift: 92.0074,
            match_window_ppm: 30.0,
        };
        let individual = run_refine(&xml, &query);
        let batch = calibrate_mass_shifts(&xml, &[query.clone(), query])
            .into_iter()
            .map(Result::unwrap)
            .collect::<Vec<_>>();
        for result in batch {
            assert_eq!(result.offset, individual.offset);
            assert_eq!(result.matched_peaks, individual.matched_peaks);
            assert_eq!(result.scan_count, individual.scan_count);
            assert!((result.refined_shift - individual.refined_shift).abs() < 1e-9);
            assert!((result.pearson.unwrap() - individual.pearson.unwrap()).abs() < 1e-12);
        }
    }

    #[test]
    fn validate_requires_proforma_or_mass_and_positive_mz() {
        let mut q = RefineMassShiftQuery {
            proforma: None,
            external_mods: Vec::new(),
            sequence_mass: Some(1000.0),
            modification_formula: Some("PO3H".to_owned()),
            modification_mass: None,
            mz: 500.0,
            charge: 2,
            rt_min: 0.0,
            rt_max: 10.0,
            tolerance_ppm: 10.0,
            prior_shift: 0.0,
            match_window_ppm: 30.0,
        };
        assert!(validate(&q).is_ok());
        q.proforma = None;
        q.sequence_mass = None;
        assert!(validate(&q).is_err());
        q.sequence_mass = Some(1000.0);
        q.charge = 0;
        assert!(validate(&q).is_err());
    }

    #[test]
    fn proforma_path_yields_exact_composition_and_mass() {
        let q = RefineMassShiftQuery {
            proforma: Some("[Acetyl]-SGRGK".to_owned()),
            external_mods: Vec::new(),
            sequence_mass: None,
            modification_formula: None,
            modification_mass: None,
            mz: 1000.0,
            charge: 2,
            rt_min: 0.0,
            rt_max: 1.0,
            tolerance_ppm: 100.0,
            prior_shift: 0.0,
            match_window_ppm: 30.0,
        };
        let parsed = parse_proforma(q.proforma.as_deref().unwrap(), &q.external_mods).unwrap();
        let c = parsed
            .composition
            .iter()
            .find(|e| e.element == "C")
            .map(|e| e.count)
            .unwrap();
        assert_eq!(c, 19 + 2);
    }

    #[test]
    fn bare_proforma_plus_shift_mass_uses_exact_base_and_pure_delta() {
        let q = RefineMassShiftQuery {
            proforma: Some("SGRGK".to_owned()),
            external_mods: Vec::new(),
            sequence_mass: None,
            modification_formula: None,
            modification_mass: Some(42.0106),
            mz: 1000.0,
            charge: 2,
            rt_min: 0.0,
            rt_max: 1.0,
            tolerance_ppm: 100.0,
            prior_shift: 42.0106,
            match_window_ppm: 30.0,
        };
        let (composition, anchor_mass, bare_mass, declared_base, declared_delta, warnings) =
            candidate_chemistry(&q).unwrap();
        let parsed = parse_proforma("SGRGK", &[]).unwrap();
        assert!((bare_mass - parsed.monoisotopic_mass).abs() < 1e-9);
        assert!((anchor_mass - bare_mass - 42.0106).abs() < 1e-9);
        assert!(declared_delta == 0.0 && (declared_base - bare_mass).abs() < 1e-12);
        // Exact base: element counts match the protein's own atoms (averagine
        // would give different C/N/O/H ratios for the same mass).
        for element in ["C", "H", "N", "O", "S"] {
            let parsed_count = parsed
                .composition
                .iter()
                .find(|e| e.element == element)
                .map(|e| e.count)
                .unwrap_or(0);
            assert_eq!(composition[element], parsed_count, "element {element}");
        }
        assert!(!warnings.iter().any(|w| w.contains("averagine")));
    }

    #[test]
    fn declared_mod_proforma_plus_residual_mass_uses_declared_base() {
        // Named tags are the fixed base; modification_mass is the unknown
        // residual in the post-declared frame (TopPIC unexpected-shift frame).
        // The anchor must sit on declared_base + residual — never on the bare
        // sequence, which would mis-place the window by the background.
        let q = RefineMassShiftQuery {
            proforma: Some("[Acetyl]-SGRGK".to_owned()),
            external_mods: Vec::new(),
            sequence_mass: None,
            modification_formula: None,
            modification_mass: Some(121.972),
            mz: 1000.0,
            charge: 2,
            rt_min: 0.0,
            rt_max: 1.0,
            tolerance_ppm: 100.0,
            prior_shift: 121.972,
            match_window_ppm: 30.0,
        };
        let (_, anchor_mass, bare_mass, declared_base, declared_delta, _) =
            candidate_chemistry(&q).unwrap();
        let parsed = parse_proforma("SGRGK", &[]).unwrap();
        assert!(
            (declared_delta - 42.010565).abs() < 1e-6,
            "{declared_delta}"
        );
        assert!((bare_mass - parsed.monoisotopic_mass).abs() < 1e-9);
        assert!((declared_base - bare_mass - 42.010565).abs() < 1e-6);
        assert!((anchor_mass - declared_base - 121.972).abs() < 1e-9);
    }

    #[test]
    fn proforma_baseline_accepts_separate_candidate_chemistry() {
        let q = RefineMassShiftQuery {
            proforma: Some("[Acetyl]-SGRGK".to_owned()),
            external_mods: Vec::new(),
            sequence_mass: None,
            modification_formula: Some("H(-1) Br".to_owned()),
            modification_mass: Some(77.910511),
            mz: 1000.0,
            charge: 2,
            rt_min: 0.0,
            rt_max: 1.0,
            tolerance_ppm: 100.0,
            prior_shift: 77.910511,
            match_window_ppm: 30.0,
        };
        let (_, candidate_mass, bare_mass, declared_base, declared_delta, _) =
            candidate_chemistry(&q).unwrap();
        let parsed = parse_proforma(q.proforma.as_deref().unwrap(), &[]).unwrap();
        let bare = parse_proforma("SGRGK", &[]).unwrap();
        assert!((declared_base - parsed.monoisotopic_mass).abs() < 1e-9);
        assert!((bare_mass - bare.monoisotopic_mass).abs() < 1e-9);
        assert!((declared_base - bare_mass - declared_delta).abs() < 1e-9);
        assert!((candidate_mass - declared_base - 77.910511).abs() < 1e-9);
    }

    #[test]
    fn theoretical_peaks_are_anchored_absolute_mz() {
        let xml = fs::read(REAL_FIXTURE).unwrap();
        let q = RefineMassShiftQuery {
            proforma: None,
            external_mods: Vec::new(),
            sequence_mass: Some(11319.3476),
            modification_formula: None,
            modification_mass: Some(91.9674),
            mz: (11411.315 + 11.0 * PROTON) / 11.0,
            charge: 11,
            rt_min: -1.0,
            rt_max: 1.0,
            tolerance_ppm: 100.0,
            prior_shift: 92.0074,
            match_window_ppm: 30.0,
        };
        let res = run_refine(&xml, &q);
        let m0_mz = res.theoretical.first().unwrap().mz;
        let expected_m0 = (11411.315 + 11.0 * PROTON) / 11.0;
        assert!(
            (m0_mz - expected_m0).abs() < 0.01,
            "theoretical M0 {m0_mz} must be anchored to the candidate mass ({expected_m0})"
        );
        let obs_m0 = res.observed.first().unwrap().mz;
        assert!(
            (obs_m0 - expected_m0).abs() < 0.5,
            "observed M0 {obs_m0} must be near the anchored theoretical ({expected_m0})"
        );
    }

    #[test]
    fn real_chlamy_h4_envelope_converges_to_toppic_feature_mass() {
        let xml = fs::read(REAL_FIXTURE).unwrap();
        let true_shift = 11411.315 - 11319.3476;
        let q = RefineMassShiftQuery {
            proforma: None,
            external_mods: Vec::new(),
            sequence_mass: Some(11319.3476),
            modification_formula: None,
            modification_mass: Some(true_shift),
            mz: (11411.315 + 11.0 * PROTON) / 11.0,
            charge: 11,
            rt_min: -1.0,
            rt_max: 1.0,
            tolerance_ppm: 100.0,
            prior_shift: true_shift + 0.04,
            match_window_ppm: 30.0,
        };
        let res = run_refine(&xml, &q);
        assert!(
            (res.refined_mass - 11411.315).abs() < 0.05,
            "converge to TopPIC feature mass 11411.315 (got {}, delta {})",
            res.refined_mass,
            (res.refined_mass - 11411.315).abs()
        );
        assert!(
            res.matched_peaks >= 3,
            "got {} matched peaks",
            res.matched_peaks
        );
    }

    #[test]
    fn real_bromine_envelope_refines_to_feature_mass() {
        let xml = fs::read("fixtures/histone48-scan1157-z14.mzML").unwrap();
        let seq_mass = 11315.3778;
        let observed_mono = 11393.314;
        let prior = observed_mono - seq_mass - 0.06;
        let q = RefineMassShiftQuery {
            proforma: None,
            external_mods: Vec::new(),
            sequence_mass: Some(seq_mass),
            modification_formula: Some("Br-H".to_owned()),
            modification_mass: Some(77.9105),
            mz: (observed_mono + 14.0 * PROTON) / 14.0,
            charge: 14,
            rt_min: -1.0,
            rt_max: 1.0,
            tolerance_ppm: 100.0,
            prior_shift: prior,
            match_window_ppm: 30.0,
        };
        let res = run_refine(&xml, &q);
        assert!(
            (res.refined_mass - observed_mono).abs() < 0.05,
            "brominated refined mass {} should match feature 67 ({})",
            res.refined_mass,
            observed_mono
        );
        assert_eq!(
            res.offset, 0,
            "a 0.06 Da prior error must not be corrected by a whole-neutron offset (got {})",
            res.offset
        );
        assert!(
            res.matched_peaks >= 10,
            "got {} matched peaks",
            res.matched_peaks
        );
    }

    #[test]
    fn real_data_ms1_envelope_cannot_shape_discriminate_br_phospho() {
        let xml = fs::read("fixtures/histone48-scan1157-z14.mzML").unwrap();
        let seq_mass = 11315.3778;
        let observed_mono = 11393.314;
        let mut pearsons = Vec::new();
        for (label, formula, mod_mass) in [
            ("Br", Some("Br-H"), 77.9105),
            ("PO3H", Some("PO3H"), 79.9663),
            ("Sulfo", Some("SO3"), 79.9568),
        ] {
            let q = RefineMassShiftQuery {
                proforma: None,
                external_mods: Vec::new(),
                sequence_mass: Some(seq_mass),
                modification_formula: formula.map(|f| f.to_owned()),
                modification_mass: Some(mod_mass),
                mz: (observed_mono + 14.0 * PROTON) / 14.0,
                charge: 14,
                rt_min: -1.0,
                rt_max: 1.0,
                tolerance_ppm: 100.0,
                prior_shift: observed_mono - seq_mass,
                match_window_ppm: 30.0,
            };
            let res = run_refine(&xml, &q);
            pearsons.push((label, res.pearson.unwrap_or(-1.0)));
        }
        pearsons.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        let spread = pearsons[0].1 - pearsons[2].1;
        eprintln!(
            "MS1 shape discrimination spread: {spread:.3} ({:?})",
            pearsons
        );
        assert!(
            spread < 0.05,
            "expected MS1 shape spread to be tiny on a large protein (got {spread:.3}, {:?})",
            pearsons
        );
    }

    #[test]
    fn brain_bromine_mplus2_is_prominent() {
        let br_ints = brain_ints("C20H40O4Br");
        let m0 = br_ints.first().copied().unwrap_or(f64::NAN);
        let m2 = br_ints.get(2).copied().unwrap_or(0.0);
        let ratio = m2 / m0;
        assert!(
            ratio > 0.5,
            "BRAIN should produce a strong Br M+2 (got {ratio})"
        );
    }
}
