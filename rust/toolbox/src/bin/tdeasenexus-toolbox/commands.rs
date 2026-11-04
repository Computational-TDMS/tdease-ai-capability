//! Command implementations. Thin mappings over the core: the only logic here is
//! argument marshalling and JSON/CSV shaping (D2 — the CLI is an adapter, not a
//! second algorithm).

use std::fs;
use std::io::Read;
use std::path::Path;

use serde::Deserialize;
use serde_json::{json, Value};

use tdeasenexus_toolbox::{
    aggregate_window_inner, aggregate_windows_inner, anchored_envelope, calibrate_mass_shift,
    calibrate_mass_shift_from_observed, calibrate_mass_shifts, chemical_from_counts, fit_envelope,
    ion_types_for_activation, parse_proforma, query_index_mzml, query_precursor_mzml,
    query_product_mzml, query_scans_mzml, query_select_mzml, query_sum_mzml, query_xic_mzml,
    search_by_mass, search_by_name, solve_combinations, theoretical_envelope,
    theoretical_fragment_list_proforma_variants, theoretical_fragments_mzpaf_proforma_variants,
    validate_precursor_mzml, validate_precursors_mzml, AnchoredEnvelope, ApprovedFilter,
    CandidateComparisonQuery, CandidateModel, EnvelopeMode, ExternalModification, FitEnvelopeQuery,
    FixedMod, FragmentMatchQuery, IonType, IsotopeBorrowPair, IsotopeEnvelopeQuery,
    MassModification, MassTolerance, MzTabWriter, NamedNeutralLoss, ParsedProForma, Peak,
    PeakOrder, PrecursorQuery, ProductQuery, RefineMassShiftQuery, ScanQuery, SearchQuery,
    SolveQuery, SolveSite, SpectrumSelector, Tol, TolUnit, UniModFilters, UniModNameEntry,
    UniModRecord, ValidatePrecursorQuery, USER_OFFSET,
};

use crate::args::{
    EnvelopeArgs, EnvelopeWindowsArgs, FitEnvelopeArgs, FragmentArgs, FragmentMode,
    FragmentTheoreticalArgs, InspectArgs, IsotopicBatchArgs, ParseProformaArgs, RefineShiftArgs,
    RefineShiftBatchArgs, SelectArgs, SpectrumArgs, SpectrumMode, TheoryArgs, ToleranceUnit,
    UnimodAddArgs, UnimodBatchSearchArgs, UnimodBatchSolveArgs, UnimodGetArgs, UnimodResiduesArgs,
    UnimodSearchArgs, UnimodSolveArgs, ValidatePrecursorArgs,
};
use crate::output::{full_digest, head_digest, OutputEnvelope, SourceDigest, SourceInfo};

/// A command's stdout is the standard JSON envelope, except for interchange
/// formats for external tools: (`--out-format mgf`) a Mascot Generic Format
/// spectrum list and (`--out-format mztab|mzpaf`) a tab-separated standard
/// document. Tabular materialization for ourselves goes through
/// `spectra.export`, never the kernel.
pub enum DispatchOutput {
    Json(OutputEnvelope),
    Mgf(String),
    /// A tab-separated standard document: the mzTab-M molecular table emitted
    /// by `fragment-theoretical --out-format mztab`.
    Tab(String),
}

/// A runtime error that maps to the stable stderr JSON contract.
pub enum CommandError {
    FileRead {
        path: String,
        error: String,
    },
    JsonParse {
        path: String,
        error: String,
    },
    ProformaParse(String),
    InvalidQuery(String),
    MissingArgument(String),
    MassSpec {
        code: String,
        message: String,
        path: Option<String>,
    },
}

impl CommandError {
    /// `{code, message, details?}` — the stable error shape (D4).
    pub fn to_json(&self) -> Value {
        let (code, message, path) = match self {
            Self::FileRead { path, error } => (
                "FILE_READ".to_owned(),
                format!("{path}: {error}"),
                Some(path.clone()),
            ),
            Self::JsonParse { path, error } => (
                "JSON_PARSE".to_owned(),
                format!("{path}: {error}"),
                Some(path.clone()),
            ),
            Self::ProformaParse(message) => ("PROFORMA_PARSE".to_owned(), message.clone(), None),
            Self::InvalidQuery(message) => {
                ("MASS_SPEC_INVALID_QUERY".to_owned(), message.clone(), None)
            }
            Self::MissingArgument(message) => ("USAGE".to_owned(), message.clone(), None),
            Self::MassSpec {
                code,
                message,
                path,
            } => (code.clone(), message.clone(), path.clone()),
        };
        match path {
            Some(path) => json!({
                "error": { "code": code, "message": message, "details": { "path": path } }
            }),
            None => json!({ "error": { "code": code, "message": message } }),
        }
    }
}

fn mass_spec_error(path: Option<&str>, error: tdeasenexus_toolbox::MassSpecError) -> CommandError {
    CommandError::MassSpec {
        code: error.code().to_owned(),
        message: error.to_string(),
        path: path.map(str::to_owned),
    }
}

pub fn dispatch(args: &crate::args::Args) -> Result<DispatchOutput, CommandError> {
    match &args.command {
        crate::args::Command::Inspect(a) => inspect(a).map(DispatchOutput::Json),
        crate::args::Command::ParseProforma(a) => parse_proforma_cmd(a).map(DispatchOutput::Json),
        crate::args::Command::Theory(a) => theory(a).map(DispatchOutput::Json),
        crate::args::Command::Envelope(a) => envelope(a).map(DispatchOutput::Json),
        crate::args::Command::EnvelopeWindows(a) => envelope_windows(a).map(DispatchOutput::Json),
        crate::args::Command::FitEnvelope(a) => fit_envelope_cmd(a).map(DispatchOutput::Json),
        crate::args::Command::RefineShift(a) => refine_shift(a).map(DispatchOutput::Json),
        crate::args::Command::RefineShiftBatch(a) => {
            refine_shift_batch(a).map(DispatchOutput::Json)
        }
        crate::args::Command::Spectrum(a) => spectrum(a),
        crate::args::Command::Select(a) => select_cmd(a).map(DispatchOutput::Json),
        crate::args::Command::Fragment(a) => fragment(a).map(DispatchOutput::Json),
        crate::args::Command::EnvelopeRegions(a) => envelope_regions(a).map(DispatchOutput::Json),
        crate::args::Command::FragmentTheoretical(a) => fragment_theoretical(a),
        crate::args::Command::UnimodGet(a) => unimod_get(a).map(DispatchOutput::Json),
        crate::args::Command::UnimodSearch(a) => unimod_search(a).map(DispatchOutput::Json),
        crate::args::Command::UnimodBatchSearch(a) => {
            unimod_batch_search(a).map(DispatchOutput::Json)
        }
        crate::args::Command::UnimodBatchSolve(a) => {
            unimod_batch_solve(a).map(DispatchOutput::Json)
        }
        crate::args::Command::UnimodSolve(a) => unimod_solve(a).map(DispatchOutput::Json),
        crate::args::Command::UnimodResidues(a) => unimod_residues(a).map(DispatchOutput::Json),
        crate::args::Command::UnimodAdd(a) => unimod_add(a).map(DispatchOutput::Json),
        crate::args::Command::IsotopicBatch(a) => isotopic_batch(a).map(DispatchOutput::Json),
        crate::args::Command::ValidatePrecursor(a) => {
            validate_precursor_cmd(a).map(DispatchOutput::Json)
        }
        crate::args::Command::ResolveShift(a) => resolve_shift(a).map(DispatchOutput::Json),
        crate::args::Command::Schema(a) => schema_cmd(a).map(DispatchOutput::Json),
    }
}

/// Every command's output schema, keyed by command name. The envelope wrapper
/// (schemaVersion/tool/version/command/params/source/elapsedMs/warnings/result)
/// is fixed; this is the shape of the `result` payload. Documentation schema —
/// `additionalProperties` stays open so adding a field never breaks a consumer
/// that fetched an older schema (the contract's `schemaVersion` is the gate).
fn result_schema(command: &str) -> Option<Value> {
    let fields: &[&str] = match command {
        "inspect" => &[
            "spectraInFile",
            "spectraInspected",
            "complete",
            "indexSource",
            "spectraPerMsLevel",
            "activationTypes",
            "retentionTimeMin",
            "retentionTimeMax",
            "retentionTimeUnit",
            "nativeIdExample",
        ],
        "parse-proforma" => &[
            "sequence",
            "monoisotopicMass",
            "composition",
            "modifications",
            "warnings",
        ],
        "theory" => &[
            "models",
            "count",
            "proforma",
            "sequence",
            "monoisotopicMass",
            "chemicalFormula",
            "composition",
            "modifications",
            "warnings",
            "envelope",
            "anchored",
        ],
        "envelope" => &["peaks", "scanCount", "charge", "mass"],
        "envelope-windows" => &["windows"],
        "fit-envelope" => &["evaluation", "charge", "fits", "count"],
        "refine-shift" => &[
            "refinedMass",
            "refinedShift",
            "declaredDeltaMass",
            "refinedShiftAfterDeclared",
            "offset",
            "matchedPeaks",
            "pearson",
            "theoretical",
            "observed",
            "evaluations",
            "scanCount",
            "warnings",
            "ambiguity",
        ],
        "refine-shift-batch" => &["results"],
        // spectrum's result depends on --mode (window/sum/xic/scans/precursor/product).
        "spectrum" => &[
            "spectra",
            "spectraMetadata",
            "scanNumbers",
            "peaks",
            "peakCount",
            "binPpm",
            "truncated",
        ],
        "select" => &[
            "plan",
            "scannedSpectrumCount",
            "matchedSpectrumCount",
            "truncated",
            "spectra",
        ],
        // fragment's result depends on --mode (match/compare).
        "fragment" => &["<mode-dependent>"],
        "fragment-theoretical" => &[
            "proforma",
            "sequence",
            "monoisotopicMass",
            "chemicalFormula",
            "modifications",
            "ionTypes",
            "maxCharge",
            "fragmentCount",
            "fragments",
            "warnings",
        ],
        "unimod-get" => &["source", "record"],
        "unimod-search" => &["source", "matches", "massMatches"],
        "unimod-batch-search" => &["source", "queries"],
        "unimod-batch-solve" => &["source", "queries"],
        "unimod-solve" => &[
            "source",
            "candidates",
            "matchCount",
            "returnedCount",
            "truncated",
            "completeness",
            "scannedCombinationCount",
            "prunedCombinationCount",
        ],
        "unimod-residues" => &["source", "aminoAcids"],
        "unimod-add" => &["source", "added", "userId"],
        "isotopic-batch" => &["variants"],
        "envelope-regions" => &["scan", "proforma", "ions", "charge", "pad"],
        "validate-precursor" => &[
            "ms2Scan",
            "ms2RtMin",
            "ms1ScanCount",
            "charge",
            "matchedPeaks",
            "minMatchedPeaks",
            "isValid",
            "pearson",
            "refinedMass",
            "envelopeEvaluations",
        ],
        "resolve-shift" => &[
            "sequence",
            "observedShift",
            "shiftFrame",
            "scans",
            "ions",
            "candidates",
            "enumeratedCombinations",
            "evaluatedCombinations",
            "notes",
        ],
        "schema" => &["$schema", "title", "type", "properties", "required"],
        _ => return None,
    };
    Some(json!({
        "type": "object",
        "properties": fields.iter().map(|name| (name.to_string(), json!({ "description": "see --help" }))).collect::<serde_json::Map<_, _>>(),
        "additionalProperties": true,
    }))
}

/// `resolve-shift` -- the sequence-graph engine (Informed-Proteomics core)
/// resolving an observed shift into positioned ProForma components: library
/// components plus an optional free residual SHIFT, one DP over
/// (position x combination). Without `--library` the library is empty and
/// the DP localizes the free SHIFT alone (single-event localization); with
/// `--library builtin|<json>|none` the shift is decomposed over the library
/// (`builtin` is the default), keeping the unexplained residual as a
/// positioned free SHIFT unless `--no-free-shift` forces closure.
/// Output ranks candidate models; it is a measurement
/// (`candidate_model_evaluation`, notIdentity), never a verdict.
fn resolve_shift(a: &crate::args::ResolveShiftArgs) -> Result<OutputEnvelope, CommandError> {
    // "none" is the explicit label-free single-event localization (the old
    // omit-`--library` behavior); the default is now builtin chemistry with
    // the residual kept open.
    let label_free = matches!(a.library.as_deref(), Some("none") | None);
    let library = match a.library.as_deref() {
        None | Some("builtin") => None,
        Some("none") => Some("[]".to_string()),
        Some(text) => Some(text.to_string()),
    };
    let bundle = crate::args::GraphSearchArgs {
        file: a.file.clone(),
        scan: a.scan.clone(),
        proforma: a.proforma.clone(),
        sequence: a.sequence.clone(),
        observed_shift: a.shift,
        library,
        max_components: a.max_components,
        no_free_shift: a.no_free_shift,
        mods: a.mods.clone(),
        ion_types: a.ion_types.clone(),
        tolerance_ppm: a.tolerance_ppm,
        mass_min: a.mass_min,
        mass_max: a.mass_max,
        scorer: a.scorer.clone(),
        overshoot_window_da: a.overshoot_window_da,
        max_candidates: a.max_candidates,
    };
    let mut envelope = run_sequence_graph(&bundle)?;
    if let Some(params) = envelope.params.as_object_mut() {
        params.insert("details".to_string(), json!(a.details));
    }
    if label_free {
        if let Some(params) = envelope.params.as_object_mut() {
            params.insert("librarySource".to_string(), json!("empty-free-shift"));
        }
    }
    if !a.details {
        compact_resolve_shift_output(&mut envelope);
    }
    Ok(envelope)
}

/// Keep stdout small enough to be consumed directly by an interactive shell
/// or Agent turn. Complete fields for the returned candidates remain available
/// with `resolve-shift --details`; this projection changes no search or ranking.
fn compact_resolve_shift_output(envelope: &mut OutputEnvelope) {
    if let Some(params) = envelope.params.as_object_mut() {
        params.remove("library");
        params.retain(|_, value| !value.is_null());
    }

    let Some(result) = envelope.result.as_object_mut() else {
        return;
    };
    let compact_warnings = result
        .get("notes")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .filter_map(Value::as_str)
        .filter_map(|note| {
            if note.contains("truncated below a") {
                Some("candidate view pruned at an evidence cliff")
            } else if note.contains("truncated by maxCandidates") {
                Some("candidate view capped by --max-candidates")
            } else {
                None
            }
        })
        .map(str::to_owned)
        .collect::<Vec<_>>();
    result.remove("sequence");
    result.remove("notes");

    let Some(candidates) = result.get_mut("candidates").and_then(Value::as_array_mut) else {
        return;
    };
    for candidate in candidates {
        let Some(candidate) = candidate.as_object_mut() else {
            continue;
        };
        candidate.remove("alternatives");
        candidate.remove("modifications");
        candidate.remove("numPlacements");
        candidate.remove("spectralEvalue");
        candidate.remove("scoreProbability");
        if let Some(cliffs) = candidate
            .get_mut("componentCliffs")
            .and_then(Value::as_array_mut)
        {
            for cliff in cliffs {
                if let Some(cliff) = cliff.as_object_mut() {
                    cliff.remove("dpSite");
                }
            }
        }
        candidate.retain(|_, value| !value.is_null());
    }
    envelope.warnings.extend(compact_warnings);
}

#[cfg(test)]
mod resolve_shift_output_tests {
    use super::*;

    #[test]
    fn compact_view_keeps_candidate_decision_fields_only() {
        let mut envelope = OutputEnvelope::new(
            "resolve-shift",
            json!({"library": [{}]}),
            Vec::new(),
            json!({
                "notes": ["1 candidates truncated by maxCandidates=1"],
                "candidates": [{
                    "proforma": "PEP[Phospho]TIDE", "score": 12.0,
                    "componentCliffs": [{"name": "Phospho", "dpSite": "3", "siteInterval": [3, 3]}],
                    "alternatives": [["large"]], "spectralEvalue": 9.0
                }]
            }),
        );

        compact_resolve_shift_output(&mut envelope);
        assert!(envelope.params.get("library").is_none());
        assert!(envelope.result.get("notes").is_none());
        let candidate = &envelope.result["candidates"][0];
        assert_eq!(candidate["proforma"], "PEP[Phospho]TIDE");
        assert!(candidate["componentCliffs"][0].get("dpSite").is_none());
        assert!(candidate.get("alternatives").is_none());
        assert!(candidate.get("spectralEvalue").is_none());
        assert_eq!(envelope.warnings.len(), 1);
    }
}

/// The sequence-graph engine behind `resolve-shift`: decompose an observed
/// mass shift over a modification library and localize every component
/// (plus the free residual SHIFT) with one DP. Output ranks candidate
/// models; it is a measurement (`candidate_model_evaluation`,
/// notIdentity), never a verdict.
/// Parse `--mods` positional fixed forms, name-first. An entry carries a
/// modification `"name"` resolved through the same ProForma/UniMod chemistry
/// authority the record path uses — callers never memorize or hand-write a
/// numeric mass. A literal `"deltaMass"` remains accepted for chemistry
/// absent from every database, but it is the exception, not the pattern.
fn parse_positional_mods(
    text: &str,
) -> Result<Vec<tdeasenexus_toolbox::PositionalModification>, CommandError> {
    let entries: Vec<Value> = parse_json_arg("--mods", text)?;
    let mut out = Vec::with_capacity(entries.len());
    for entry in entries {
        let position = entry
            .get("position")
            .and_then(Value::as_str)
            .ok_or_else(|| {
                CommandError::InvalidQuery(
                    "--mods entries require a position (nTerm, cTerm, or a 1-based residue index)"
                        .to_owned(),
                )
            })?
            .to_owned();
        let delta_mass = match entry.get("deltaMass").and_then(Value::as_f64) {
            Some(mass) => mass,
            None => {
                let name = entry.get("name").and_then(Value::as_str).ok_or_else(|| {
                    CommandError::InvalidQuery(
                        "--mods entries carry a modification name (resolved by the chemistry \
                         authority, e.g. {\"position\":\"nTerm\",\"name\":\"Acetyl\"}) — a \
                         hand-written deltaMass only for chemistry absent from all databases"
                            .to_owned(),
                    )
                })?;
                tdeasenexus_toolbox::mod_delta_mass(name)
                    .map_err(|error| mass_spec_error(None, error))?
            }
        };
        out.push(tdeasenexus_toolbox::PositionalModification {
            position,
            delta_mass,
        });
    }
    Ok(out)
}

fn run_sequence_graph(a: &crate::args::GraphSearchArgs) -> Result<OutputEnvelope, CommandError> {
    use tdeasenexus_toolbox::{graph_search_msalign, GraphSearchQuery};

    let (bytes, source) = read_file(&a.file)?;
    let (library, library_source) = match a.library.as_deref() {
        Some(text) => (parse_graph_library(text, &a.file)?, "explicit"),
        None => (common_graph_library()?, "builtin-common"),
    };
    let fixed: Vec<tdeasenexus_toolbox::PositionalModification> = match &a.mods {
        Some(text) => parse_positional_mods(text)?,
        None => Vec::new(),
    };
    // Frame gate: TopPIC unexpected shifts are residuals measured AFTER the
    // assigned N-terminal form (NONE/NME/NME_ACETYLATION/M_ACETYLATION) —
    // the terminal background is hidden inside that residual. This kernel
    // never injects or defaults it: the caller must capture the source
    // row's terminal form back into an EXPLICIT annotation. A bare
    // sequence+shift call (no --proforma, no --mods) is frame-degenerate
    // and refused. The preferred form is the complete ProForma from the
    // record (--proforma: named tags become fixed forms, a numeric bracket
    // is the observed shift); --mods accepts the same named chemistry.
    if a.proforma.is_none() && a.mods.is_none() && a.observed_shift.is_some() {
        return Err(CommandError::InvalidQuery(
            "resolve-shift requires an explicit frame annotation for --shift: TopPIC \
             unexpected shifts are residuals after the assigned N-terminal form, and this \
             kernel never injects that hidden background. Capture the source row's \
             terminal form yourself — preferred: pass the complete ProForma (--proforma), \
             where named tags such as [Acetyl]- become fixed forms and a numeric bracket \
             is the observed shift. Alternative: --mods with the terminal form as a NAMED \
             modification ({\"position\":\"nTerm\",\"name\":\"Acetyl\"}); use [] to \
             acknowledge a NONE terminal form. Never hand-write a numeric mass"
                .to_owned(),
        ));
    }
    let query = GraphSearchQuery {
        file: a.file.clone(),
        proforma: a.proforma.clone(),
        scans: a.scan.clone(),
        sequence: a.sequence.clone().unwrap_or_default(),
        fixed_modifications: fixed,
        library: library.clone(),
        max_components: a.max_components,
        observed_shift: a.observed_shift.unwrap_or(0.0),
        no_free_shift: a.no_free_shift,
        tolerance_ppm: a.tolerance_ppm,
        ion_types: a.ion_types.clone(),
        mass_min: a.mass_min,
        mass_max: a.mass_max,
        scorer: a.scorer.clone(),
        overshoot_window_da: a.overshoot_window_da,
        max_candidates: a.max_candidates,
    };
    let result = graph_search_msalign(&bytes, &query)
        .map_err(|error| mass_spec_error(Some(a.file.as_str()), error))?;
    Ok(OutputEnvelope::new(
        "resolve-shift",
        json!({
            "file": a.file, "scan": a.scan, "proforma": a.proforma,
            "sequence": a.sequence, "observedShift": a.observed_shift,
            "library": library, "librarySource": library_source,
            "mods": a.mods, "maxComponents": a.max_components,
            "noFreeShift": a.no_free_shift, "ionTypes": a.ion_types,
            "tolerancePpm": a.tolerance_ppm, "scorer": a.scorer,
        }),
        vec![source],
        serde_json::to_value(result).unwrap(),
    ))
}

/// Compact OpenSearch PTM candidates. This is intentionally a small, explicit
/// `LibrarySpec` set rather than a substitute for UniMod or a chemistry claim.
/// Its JSON shape is the public `--library` template: callers can copy it,
/// remove entries, add named modifications, or add `positions` constraints.
fn common_graph_library() -> Result<Vec<tdeasenexus_toolbox::LibrarySpec>, CommandError> {
    use tdeasenexus_toolbox::{mod_delta_mass, LibrarySpec};
    let entry = |name: &str, chemistry: &str, residues: &[char], location: Option<&str>| {
        Ok(LibrarySpec {
            name: name.to_owned(),
            delta_mass: mod_delta_mass(chemistry).map_err(|error| mass_spec_error(None, error))?,
            residues: residues.to_vec(),
            location: location.map(str::to_owned),
            positions: None,
        })
    };
    [
        entry("Acetyl@K", "Acetyl", &['K'], None),
        entry("Acetyl@nTerm", "Acetyl", &['*'], Some("proteinNTerm")),
        entry("Methyl@K/R", "Methyl", &['K', 'R'], None),
        entry("Dimethyl@K/R", "Dimethyl", &['K', 'R'], None),
        entry("Trimethyl@K", "Trimethyl", &['K'], None),
        entry("Phospho@S/T/Y", "Phospho", &['S', 'T', 'Y'], None),
        entry("Bromo@F/H/W/Y", "Bromo", &['F', 'H', 'W', 'Y'], None),
        entry("Oxidation@M", "Oxidation", &['M'], None),
        entry("Deamidation@N/Q", "Deamidation", &['N', 'Q'], None),
    ]
    .into_iter()
    .collect()
}

fn parse_graph_library(
    text: &str,
    file: &str,
) -> Result<Vec<tdeasenexus_toolbox::LibrarySpec>, CommandError> {
    let mut library_value: serde_json::Value = parse_json_arg("--library", text)?;
    let entries = library_value.as_array_mut().ok_or_else(|| {
        CommandError::InvalidQuery(
            "--library must be a JSON array of LibrarySpec entries".to_owned(),
        )
    })?;
    for entry in entries {
        if let Some(residues) = entry.get_mut("residues") {
            // bare-string shorthand: the string's characters are the targets
            if let Some(text) = residues.as_str() {
                let chars: Vec<Value> =
                    text.chars().map(|c| Value::String(c.to_string())).collect();
                *residues = Value::Array(chars);
            }
        }
        if entry.get("deltaMass").is_none() {
            if let Some(name) = entry.get("name").and_then(|v| v.as_str()) {
                let mass = tdeasenexus_toolbox::mod_delta_mass(name)
                    .map_err(|error| mass_spec_error(Some(file), error))?;
                entry["deltaMass"] = json!(mass);
            }
        }
    }
    serde_json::from_value(library_value).map_err(|error| CommandError::JsonParse {
        path: "--library".to_owned(),
        error: error.to_string(),
    })
}

/// `schema <command>` — print the JSON Schema (draft-07) for a command's output
/// envelope, so a caller (the Python shim, a skill) can check a response shape
/// without guessing. `--help` is the prose interface; this is the machine form.
/// One variant in an `isotopic-batch` request.
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct BatchVariant {
    formula: String,
    charge: i32,
    #[serde(default = "default_batch_npeaks")]
    npeaks: usize,
}

fn default_batch_npeaks() -> usize {
    4
}

fn isotopic_batch(a: &IsotopicBatchArgs) -> Result<OutputEnvelope, CommandError> {
    let (bytes, source) = read_file(&a.file)?;
    let variants: Vec<BatchVariant> =
        serde_json::from_slice(&bytes).map_err(|error| CommandError::JsonParse {
            path: a.file.clone(),
            error: error.to_string(),
        })?;
    if !a.truncate_after.is_finite()
        || !(0.0..=1.0).contains(&a.truncate_after)
        || a.truncate_after == 0.0
        || !a.ignore_below.is_finite()
        || !(0.0..=1.0).contains(&a.ignore_below)
    {
        return Err(CommandError::InvalidQuery(
            "isotope thresholds require truncateAfter in (0,1] and ignoreBelow in [0,1]".into(),
        ));
    }
    if variants.is_empty() {
        return Err(CommandError::InvalidQuery(
            "isotopic-batch requires nonempty variants".into(),
        ));
    }
    let count = variants.len();
    let prepared: Vec<_> = variants
        .into_iter()
        .enumerate()
        .map(|(index, variant)| {
            if variant.charge <= 0 || variant.npeaks == 0 {
                return Err(CommandError::InvalidQuery(format!(
                    "variant {index} requires positive charge and npeaks"
                )));
            }
            let composition = chemical_elements::ChemicalComposition::parse(&variant.formula)
                .map_err(|error| {
                    CommandError::InvalidQuery(format!(
                        "variant {index}: invalid elemental formula {}: {error}",
                        variant.formula
                    ))
                })?;
            Ok((variant, composition))
        })
        .collect::<Result<_, CommandError>>()?;
    let mut out_variants = Vec::with_capacity(count);
    for (variant, composition) in prepared {
        // BRAIN envelope (M0 always kept by the theory M0-preservation fix),
        // truncated by the caller's thresholds, then take the first npeaks
        // isotopologues (M0..M_{n-1}) — the shape a flat-`npeaks` caller
        // expects (matches brainpy's npeaks semantics, with exact m/z + true
        // relative intensity).
        let env = theoretical_envelope(
            &composition,
            variant.charge,
            a.truncate_after,
            a.ignore_below,
        );
        let peaks: Vec<Value> = env
            .peaks
            .iter()
            .take(variant.npeaks)
            .map(|p| json!({"mz": p.mz, "intensity": p.intensity}))
            .collect();
        out_variants.push(json!({
            "formula": variant.formula,
            "charge": variant.charge,
            "npeaks": variant.npeaks,
            "peaks": peaks,
        }));
    }
    let result = json!({ "variants": out_variants });
    let out = OutputEnvelope::new(
        "isotopic-batch",
        json!({
            "file": a.file,
            "truncateAfter": a.truncate_after,
            "ignoreBelow": a.ignore_below,
            "count": count,
        }),
        vec![source],
        result,
    );
    Ok(out)
}

fn schema_cmd(a: &crate::args::SchemaArgs) -> Result<OutputEnvelope, CommandError> {
    let result = result_schema(&a.command).ok_or_else(|| {
        CommandError::InvalidQuery(format!(
            "unknown command {:?}: run `tdeasenexus-toolbox --help` for the public command list",
            a.command
        ))
    })?;
    let schema = json!({
        "$schema": "http://json-schema.org/draft-07/schema#",
        "title": format!("tdeasenexus-toolbox {}", a.command),
        "x-batch-inputs": {
            "theory": "--models JSON array; each row has proforma, optional id/charge/anchor; row values override shared flags",
            "fit-envelope": "--pairs JSON array; each row has observed/theoretical paths and optional id/matchWindowPpm; files are reused",
            "fragment": "--candidate-proformas JSON models and --scans-input JSON scan ids; one model and scan lists also work; --sum-scans is explicit",
            "envelope": "envelope-windows --windows JSON",
            "calibration": "refine-shift-batch --input JSON",
            "isotopes": "isotopic-batch JSON file",
            "fragment-theoretical": "--models JSON models with shared ion-types/max-charge; batch output is json",
            "json": "Array flags accept inline JSON or file paths; positional isotopic-batch input is a file; scan arrays accept integers or native-id strings",
        },
        "x-observation-contract": {
            "fragmentInputs": ["centroided mzML", "TopFD msalign", "spectrum sum/window JSON"],
            "selection": "scan ids, native scan-number bounds and RT bounds intersect; RT unit is minute",
            "pooling": "explicit sumScans; binPpm is distinct from matching tolerancePpm",
            "rawMatching": "isotope envelope; full-slot Pearson, missing slots zero, theoretical apex required",
            "deconvolutedMatching": "neutral mass and charge; no isotope-shape measurement",
            "comparison": "each model owns its backbone; positional ion comparison requires sharedBackbone and ionComparisonAvailable",
            "observation": "spectraMetadata and ordered peaks are reusable; truncated observations cannot enter fragment matching",
        },
        "x-mass-constants": {
            "protonMass": tdeasenexus_toolbox::PROTON,
            "carbon13Spacing": tdeasenexus_toolbox::NEUTRON_SHIFT,
        },
        "type": "object",
        "properties": {
            "schemaVersion": { "const": 1 },
            "tool": { "const": "tdeasenexus-toolbox" },
            "version": { "type": "string" },
            "command": { "const": a.command },
            "params": { "type": "object" },
            "source": { "type": "array" },
            "elapsedMs": { "type": "integer" },
            "warnings": { "type": "array", "items": { "type": "string" } },
            "result": result,
        },
        "required": ["schemaVersion", "tool", "version", "command", "result"],
    });
    Ok(OutputEnvelope::new(
        "schema",
        json!({ "command": a.command }),
        Vec::new(),
        schema,
    ))
}

/// Read a file into memory with provenance (source path, size, head digest).
fn read_file(path: &str) -> Result<(Vec<u8>, SourceInfo), CommandError> {
    let bytes = fs::read(path).map_err(|error| CommandError::FileRead {
        path: path.to_owned(),
        error: error.to_string(),
    })?;
    let source = SourceInfo {
        path: path.to_owned(),
        size: bytes.len() as u64,
        digest: SourceDigest {
            algo: "sha256:head1MiB".to_owned(),
            value: head_digest(&bytes),
        },
    };
    Ok((bytes, source))
}

/// Read a JSON file into a typed value with provenance.
fn read_json<T: serde::de::DeserializeOwned>(path: &str) -> Result<(T, SourceInfo), CommandError> {
    let (bytes, source) = read_file(path)?;
    let value = serde_json::from_slice(&bytes).map_err(|error| CommandError::JsonParse {
        path: path.to_owned(),
        error: error.to_string(),
    })?;
    Ok((value, source))
}

/// Load an external-mods JSON file: `[{"name": "...", "formula": "..."}]`.
fn load_mods(
    path: Option<&str>,
) -> Result<(Vec<ExternalModification>, Vec<SourceInfo>), CommandError> {
    let Some(path) = path else {
        return Ok((Vec::new(), Vec::new()));
    };
    let (mods, source) = read_json::<Vec<ExternalModification>>(path)?;
    Ok((mods, vec![source]))
}

/// Read an anchored theoretical envelope: a bare `AnchoredEnvelope` JSON, a
/// `theory`-command object `{envelope, anchored}`, or a full CLI envelope whose
/// `result` holds `{envelope, anchored}`.
fn read_anchored(path: &str) -> Result<(AnchoredEnvelope, SourceInfo), CommandError> {
    let (bytes, source) = read_file(path)?;
    let mut value: Value =
        serde_json::from_slice(&bytes).map_err(|error| CommandError::JsonParse {
            path: path.to_owned(),
            error: error.to_string(),
        })?;
    if let Value::Object(map) = &mut value {
        if let Some(result) = map.get_mut("result") {
            value = result.take();
        }
        if let Value::Object(map) = &mut value {
            if let Some(anchored) = map.get_mut("anchored") {
                value = anchored.take();
            }
        }
    }
    let anchored: AnchoredEnvelope =
        serde_json::from_value(value).map_err(|error| CommandError::JsonParse {
            path: path.to_owned(),
            error: error.to_string(),
        })?;
    Ok((anchored, source))
}

/// An `envelope`-command output reused via `refine-shift --observed`.
#[derive(Deserialize)]
struct ObservedEnvelopeFile {
    peaks: Vec<Peak>,
    #[serde(rename = "scanCount", default)]
    scan_count: usize,
}

/// Read an observed envelope: a bare `[{"mz":..,"intensity":..}]` array, an
/// `envelope`-command object `{peaks, scanCount, ...}`, or a full CLI envelope
/// whose `result` holds `{peaks, scanCount, ...}`.
fn read_observed(path: &str) -> Result<(Vec<Peak>, usize, SourceInfo), CommandError> {
    let (bytes, source) = read_file(path)?;
    let mut value: Value =
        serde_json::from_slice(&bytes).map_err(|error| CommandError::JsonParse {
            path: path.to_owned(),
            error: error.to_string(),
        })?;
    match &mut value {
        Value::Array(_) => {
            let peaks: Vec<Peak> =
                serde_json::from_value(value).map_err(|error| CommandError::JsonParse {
                    path: path.to_owned(),
                    error: error.to_string(),
                })?;
            Ok((peaks, 0, source))
        }
        Value::Object(map) => {
            // A full CLI envelope nests the command result under `result`.
            if let Some(result) = map.get_mut("result") {
                value = result.take();
            }
            let envelope: ObservedEnvelopeFile =
                serde_json::from_value(value).map_err(|error| CommandError::JsonParse {
                    path: path.to_owned(),
                    error: error.to_string(),
                })?;
            Ok((envelope.peaks, envelope.scan_count, source))
        }
        _ => Err(CommandError::JsonParse {
            path: path.to_owned(),
            error: "expected an array of peaks or an envelope {peaks, scanCount}".to_owned(),
        }),
    }
}

fn inspect(a: &InspectArgs) -> Result<OutputEnvelope, CommandError> {
    let (bytes, _) = read_file(&a.file)?;
    let result = query_index_mzml(&bytes, a.max_spectra)
        .map_err(|error| mass_spec_error(Some(&a.file), error))?;
    let source = SourceInfo {
        path: a.file.clone(),
        size: bytes.len() as u64,
        digest: SourceDigest {
            algo: if a.digest_full {
                "sha256:full"
            } else {
                "sha256:head1MiB"
            }
            .to_owned(),
            value: if a.digest_full {
                full_digest(&bytes)
            } else {
                head_digest(&bytes)
            },
        },
    };
    Ok(OutputEnvelope::new(
        "inspect",
        json!({ "file": a.file, "maxSpectra": a.max_spectra, "digestFull": a.digest_full }),
        vec![source],
        serde_json::to_value(result).unwrap(),
    ))
}

fn parse_proforma_cmd(a: &ParseProformaArgs) -> Result<OutputEnvelope, CommandError> {
    // Default: the single core parser (mzcore backend, embedded UniMod
    // ontology, TopPIC flanking format). The CLI never resolves tags itself.
    // An explicit --db/--user-db switches to the core-exported dual-source
    // resolver ([USER:<id>] tags, name-collision ambiguity floor).
    if a.db.is_some() || a.user_db.is_some() {
        let db_path = crate::unimod_db::resolve_db_path(a.db.as_deref())?;
        let (records, names, _specificities) = load_merged_sources(&db_path, a.user_db.as_deref())?;
        let db = crate::unimod_db::UnimodDb::open(&db_path)?;
        let meta = db.meta()?;
        let amino_acids = db.amino_acids()?;
        let resolver =
            tdeasenexus_toolbox::UniModProFormaResolver::new(&records, &names, &amino_acids);
        let result = resolver
            .resolve(&a.proforma)
            .map_err(|error| CommandError::ProformaParse(error.to_string()))?;
        let mut result = serde_json::to_value(result).unwrap();
        if let Value::Object(map) = &mut result {
            map.insert("source".to_owned(), db_source(&meta));
        }
        return Ok(OutputEnvelope::new(
            "parse-proforma",
            json!({ "proforma": a.proforma, "mods": a.mods, "db": a.db, "userDb": a.user_db }),
            vec![db_source_info(&db_path)?],
            result,
        ));
    }
    let (mods, source) = load_mods(a.mods.as_deref())?;
    let parsed = parse_proforma(&a.proforma, &mods)
        .map_err(|error| CommandError::ProformaParse(error.to_string()))?;
    Ok(OutputEnvelope::new(
        "parse-proforma",
        json!({ "proforma": a.proforma, "mods": a.mods }),
        source,
        serde_json::to_value(parsed).unwrap(),
    ))
}

#[derive(serde::Serialize, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
struct TheoryModel {
    id: Option<String>,
    label: Option<String>,
    proforma: String,
    charge: Option<i32>,
    anchor: Option<f64>,
}

fn theory(a: &TheoryArgs) -> Result<OutputEnvelope, CommandError> {
    if !a.truncate_after.is_finite()
        || a.truncate_after <= 0.0
        || a.truncate_after > 1.0
        || !a.ignore_below.is_finite()
        || !(0.0..=1.0).contains(&a.ignore_below)
    {
        return Err(CommandError::InvalidQuery(
            "theory thresholds require truncateAfter in (0,1] and ignoreBelow in [0,1]".into(),
        ));
    }
    let models: Vec<TheoryModel> = if let Some(value) = &a.models {
        parse_json_arg("--models", value)?
    } else {
        vec![TheoryModel {
            id: None,
            label: None,
            proforma: a.proforma.clone().ok_or_else(|| {
                CommandError::MissingArgument("theory requires --proforma or --models".into())
            })?,
            charge: None,
            anchor: None,
        }]
    };
    if models.is_empty() {
        return Err(CommandError::InvalidQuery(
            "theory requires nonempty models".into(),
        ));
    }
    let (mods, source) = load_mods(a.mods.as_deref())?;
    let mut ids = std::collections::HashSet::new();
    let prepared: Vec<_> = models.iter().enumerate().map(|(index, model)| {
        let id = model.id.clone().unwrap_or_else(|| format!("model{}", index + 1));
        if id.trim().is_empty() || !ids.insert(id.clone()) { return Err(CommandError::InvalidQuery("theory model ids must be nonempty and unique".into())); }
        let charge = model.charge.unwrap_or(a.charge);
        let anchor = model.anchor.or(a.anchor);
        if charge <= 0 || anchor.is_some_and(|v| !v.is_finite() || v <= 0.0) {
            return Err(CommandError::InvalidQuery(format!("model {id}: charge and any anchor mass must be positive and finite")));
        }
        let parsed = parse_proforma(&model.proforma, &mods)
            .map_err(|error| CommandError::ProformaParse(format!("model {id}: {error}")))?;
        if !parsed.warnings.is_empty() {
            return Err(CommandError::InvalidQuery(format!("model {id}: theory asserts an isotope shape, which a numeric mass-only tag cannot define; use named or Formula: tags")));
        }
        Ok((id, charge, anchor, parsed))
    }).collect::<Result<_, CommandError>>()?;
    let generate = |proforma: &str, charge: i32, anchor: Option<f64>, parsed: &ParsedProForma| {
        let composition = chemical_from_counts(&parsed.composition);
        let envelope = theoretical_envelope(&composition, charge, a.truncate_after, a.ignore_below);
        let anchored = anchored_envelope(&envelope, anchor.unwrap_or(parsed.monoisotopic_mass));
        json!({"proforma": proforma, "sequence": parsed.sequence,
            "monoisotopicMass": parsed.monoisotopic_mass, "chemicalFormula": net_formula_string(&parsed.composition),
            "composition": parsed.composition, "modifications": parsed.modifications, "warnings": parsed.warnings,
            "envelope": envelope, "anchored": anchored})
    };
    let result = if a.models.is_some() {
        let rows: Vec<_> = models
            .iter()
            .zip(&prepared)
            .map(|(model, (id, charge, anchor, parsed))| {
                json!({"id": id, "label": model.label, "charge": charge,
                "measurement": generate(&model.proforma, *charge, *anchor, parsed)})
            })
            .collect();
        json!({"models": rows, "count": rows.len()})
    } else {
        let (_, charge, anchor, parsed) = &prepared[0];
        generate(&models[0].proforma, *charge, *anchor, parsed)
    };
    Ok(OutputEnvelope::new(
        "theory",
        json!({"proforma": a.proforma, "models": if a.models.is_some() { Some(&models) } else { None },
        "mods": a.mods, "charge": a.charge, "anchor": a.anchor,
        "truncateAfter": a.truncate_after, "ignoreBelow": a.ignore_below}),
        source,
        result,
    ))
}

fn envelope(a: &EnvelopeArgs) -> Result<OutputEnvelope, CommandError> {
    let (bytes, source) = read_file(&a.file)?;
    let query = IsotopeEnvelopeQuery {
        mass: Some(a.mass),
        proforma: None,
        external_mods: Vec::new(),
        charge: a.charge,
        rt_min: a.rt_min,
        rt_max: a.rt_max,
        tolerance_ppm: a.tolerance_ppm,
        isotopes: a.isotopes,
        max_offset: 0,
        theory_locked: false,
    };
    let params = json!({
        "file": a.file, "mass": a.mass, "charge": a.charge, "rtMin": a.rt_min,
        "rtMax": a.rt_max, "tolerancePpm": a.tolerance_ppm, "isotopes": a.isotopes,
    });
    let (peaks, scan_count) = aggregate_window_inner(&bytes, &query, a.mass)
        .map_err(|error| mass_spec_error(Some(&a.file), error))?;
    let result = json!({
        "peaks": peaks, "scanCount": scan_count, "charge": a.charge, "mass": a.mass,
    });
    Ok(OutputEnvelope::new(
        "envelope",
        params,
        vec![source],
        result,
    ))
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
struct EnvelopeWindowArg {
    label: Option<String>,
    mass: f64,
    charge: u8,
    rt_min: f64,
    rt_max: f64,
}

/// Aggregate several observed MS1 isotope envelopes in one mzML traversal —
/// the multiplicity-ladder surface (one file load for the whole ladder;
/// per-request matching is exactly the single-window command's).
fn envelope_windows(a: &EnvelopeWindowsArgs) -> Result<OutputEnvelope, CommandError> {
    let windows: Vec<EnvelopeWindowArg> = parse_json_arg("--windows", &a.windows)?;
    if windows.is_empty() {
        return Err(CommandError::InvalidQuery(
            "envelope-windows requires at least one window".to_owned(),
        ));
    }
    let (bytes, source) = read_file(&a.file)?;
    let queries: Vec<IsotopeEnvelopeQuery> = windows
        .iter()
        .map(|window| IsotopeEnvelopeQuery {
            mass: Some(window.mass),
            proforma: None,
            external_mods: Vec::new(),
            charge: window.charge,
            rt_min: window.rt_min,
            rt_max: window.rt_max,
            tolerance_ppm: a.tolerance_ppm,
            isotopes: a.isotopes,
            max_offset: 0,
            theory_locked: false,
        })
        .collect();
    let requests: Vec<(&IsotopeEnvelopeQuery, f64)> = queries
        .iter()
        .zip(windows.iter())
        .map(|(query, window)| (query, window.mass))
        .collect();
    let outcomes = aggregate_windows_inner(&bytes, &requests);
    let results: Vec<serde_json::Value> = windows
        .iter()
        .zip(outcomes)
        .map(|(window, outcome)| match outcome {
            Ok((peaks, scan_count)) => json!({
                "label": window.label, "mass": window.mass, "charge": window.charge,
                "rtMin": window.rt_min, "rtMax": window.rt_max,
                "peaks": peaks, "scanCount": scan_count,
            }),
            Err(error) => json!({
                "label": window.label, "mass": window.mass, "charge": window.charge,
                "rtMin": window.rt_min, "rtMax": window.rt_max,
                "error": error.to_string(),
            }),
        })
        .collect();
    Ok(OutputEnvelope::new(
        "envelope-windows",
        json!({
            "file": a.file, "windows": a.windows,
            "tolerancePpm": a.tolerance_ppm, "isotopes": a.isotopes,
        }),
        vec![source],
        json!({ "windows": results }),
    ))
}

/// Map a CLI peak-order token onto the kernel's `PeakOrder` (case-insensitive).
fn parse_peak_order(text: &str) -> Result<PeakOrder, CommandError> {
    match text.to_ascii_lowercase().as_str() {
        "mzascending" => Ok(PeakOrder::MzAscending),
        "intensitydescending" => Ok(PeakOrder::IntensityDescending),
        "intensitycappedmzascending" => Ok(PeakOrder::IntensityCappedMzAscending),
        _ => Err(CommandError::InvalidQuery(format!(
            "unknown peak order {text:?}; expected mzAscending | intensityDescending | intensityCappedMzAscending"
        ))),
    }
}

#[derive(serde::Serialize, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
struct EnvelopeFitPair {
    id: Option<String>,
    observed: String,
    theoretical: String,
    match_window_ppm: Option<f64>,
    minimum_window_mz: Option<f64>,
}

fn fit_envelope_cmd(a: &FitEnvelopeArgs) -> Result<OutputEnvelope, CommandError> {
    let pairs: Vec<EnvelopeFitPair> = if let Some(value) = &a.pairs {
        parse_json_arg("--pairs", value)?
    } else {
        vec![EnvelopeFitPair {
            id: None,
            observed: a.observed.clone().ok_or_else(|| {
                CommandError::MissingArgument("fit-envelope requires --observed or --pairs".into())
            })?,
            theoretical: a.theoretical.clone().ok_or_else(|| {
                CommandError::MissingArgument(
                    "fit-envelope requires --theoretical or --pairs".into(),
                )
            })?,
            match_window_ppm: None,
            minimum_window_mz: None,
        }]
    };
    if pairs.is_empty() {
        return Err(CommandError::InvalidQuery(
            "fit-envelope requires nonempty pairs".into(),
        ));
    }
    let mut ids = std::collections::HashSet::new();
    let requests: Vec<_> = pairs
        .iter()
        .enumerate()
        .map(|(index, pair)| {
            let id = pair
                .id
                .clone()
                .unwrap_or_else(|| format!("fit{}", index + 1));
            let tolerance = pair.match_window_ppm.unwrap_or(a.match_window_ppm);
            let minimum_window_mz = pair.minimum_window_mz.unwrap_or(a.minimum_window_mz);
            if id.trim().is_empty()
                || !ids.insert(id.clone())
                || !tolerance.is_finite()
                || tolerance <= 0.0
                || !minimum_window_mz.is_finite()
                || minimum_window_mz < 0.0
            {
                return Err(CommandError::InvalidQuery(format!(
                    "fit {id}: unique nonempty id and positive finite matchWindowPpm required"
                )));
            }
            Ok((
                id,
                FitEnvelopeQuery {
                    match_window_ppm: tolerance,
                    minimum_window_mz,
                },
            ))
        })
        .collect::<Result<_, CommandError>>()?;
    let mut observed_cache = std::collections::HashMap::new();
    let mut theory_cache = std::collections::HashMap::new();
    let mut sources = Vec::new();
    // Read every input before computing the first fit; repeated observations
    // and theories share their immutable data throughout the batch.
    for pair in &pairs {
        if !observed_cache.contains_key(&pair.observed) {
            let (peaks, _, source) = read_observed(&pair.observed)?;
            observed_cache.insert(pair.observed.clone(), peaks);
            if !sources
                .iter()
                .any(|existing: &SourceInfo| existing.path == source.path)
            {
                sources.push(source);
            }
        }
        if !theory_cache.contains_key(&pair.theoretical) {
            let (theory, source) = read_anchored(&pair.theoretical)?;
            theory_cache.insert(pair.theoretical.clone(), theory);
            if !sources
                .iter()
                .any(|existing: &SourceInfo| existing.path == source.path)
            {
                sources.push(source);
            }
        }
    }
    for (pair, (_, query)) in pairs.iter().zip(&requests) {
        tdeasenexus_toolbox::validate_envelope_fit(
            &observed_cache[&pair.observed],
            &theory_cache[&pair.theoretical],
            query,
        )
        .map_err(|error| mass_spec_error(None, error))?;
    }
    let evaluations: Vec<_> = pairs.iter().zip(&requests).map(|(pair, (id, query))| {
        let theoretical = &theory_cache[&pair.theoretical];
        let evaluation = fit_envelope(&observed_cache[&pair.observed], theoretical, query)
            .map_err(|error| mass_spec_error(None, error))?;
        Ok(json!({"id": id, "observed": pair.observed, "theoretical": pair.theoretical,
            "matchWindowPpm": query.match_window_ppm, "minimumWindowMz": query.minimum_window_mz, "evaluation": evaluation, "charge": theoretical.charge}))
    }).collect::<Result<_, CommandError>>()?;
    let result = if a.pairs.is_some() {
        json!({"fits": evaluations, "count": evaluations.len()})
    } else {
        json!({"evaluation": evaluations[0]["evaluation"], "charge": evaluations[0]["charge"]})
    };
    Ok(OutputEnvelope::new(
        "fit-envelope",
        json!({"observed": a.observed, "theoretical": a.theoretical,
        "pairs": if a.pairs.is_some() { Some(&pairs) } else { None }, "matchWindowPpm": a.match_window_ppm, "minimumWindowMz": a.minimum_window_mz}),
        sources,
        result,
    ))
}

fn validate_precursor_cmd(a: &ValidatePrecursorArgs) -> Result<OutputEnvelope, CommandError> {
    let (external_mods, mut sources) = load_mods(a.mods.as_deref())?;
    // Evidence lattice: the precursor validation is an isotope-envelope test,
    // so a mass-only delta — which cannot define the shape — is rejected here
    // instead of silently matching against a base-composition envelope.
    let parsed = parse_proforma(&a.proforma, &external_mods)
        .map_err(|error| CommandError::ProformaParse(error.to_string()))?;
    if !parsed.warnings.is_empty() {
        return Err(CommandError::InvalidQuery(
            "validate-precursor is an isotope-envelope test; a numeric mass-only tag cannot \
             define the shape — use a named or Formula: tag"
                .to_owned(),
        ));
    }
    let (bytes, source) = read_file(&a.file)?;
    sources.insert(0, source);
    let query = ValidatePrecursorQuery {
        proforma: a.proforma.clone(),
        external_mods,
        charge: a.charge,
        scan: a.scan.clone(),
        tolerance_ppm: a.tolerance_ppm,
        min_matched_peaks: a.min_matched_peaks,
    };
    let scans: Vec<String> = a.scan.split(',').map(str::to_owned).collect();
    let result = if scans.len() == 1 {
        serde_json::to_value(
            validate_precursor_mzml(&bytes, &query)
                .map_err(|error| mass_spec_error(Some(&a.file), error))?,
        )
        .unwrap()
    } else {
        serde_json::to_value(
            validate_precursors_mzml(&bytes, &query, &scans)
                .map_err(|error| mass_spec_error(Some(&a.file), error))?,
        )
        .unwrap()
    };
    Ok(OutputEnvelope::new(
        "validate-precursor",
        json!({
            "file": a.file, "proforma": a.proforma, "mods": a.mods, "charge": a.charge,
            "scan": a.scan, "tolerancePpm": a.tolerance_ppm,
            "minMatchedPeaks": a.min_matched_peaks,
        }),
        sources,
        result,
    ))
}

fn refine_shift(a: &RefineShiftArgs) -> Result<OutputEnvelope, CommandError> {
    let (external_mods, mut sources) = load_mods(a.mods.as_deref())?;
    let query = RefineMassShiftQuery {
        proforma: a.proforma.clone(),
        external_mods,
        sequence_mass: a.sequence_mass,
        modification_mass: a.mod_mass,
        modification_formula: a.mod_formula.clone(),
        mz: a.mz,
        charge: a.charge,
        rt_min: a.rt_min,
        rt_max: a.rt_max,
        tolerance_ppm: a.tolerance_ppm,
        prior_shift: a.prior_shift,
        match_window_ppm: a.match_window_ppm,
    };
    let result = match &a.observed {
        // Reuse a pre-extracted envelope: evaluate N hypotheses without
        // re-reading the mzML (D5).
        Some(path) => {
            let (peaks, scan_count, source) = read_observed(path)?;
            sources.push(source);
            calibrate_mass_shift_from_observed(&peaks, scan_count, &query)
                .map_err(|error| mass_spec_error(Some(path), error))?
        }
        None => {
            let (bytes, source) = read_file(&a.file)?;
            sources.insert(0, source);
            calibrate_mass_shift(&bytes, &query)
                .map_err(|error| mass_spec_error(Some(&a.file), error))?
        }
    };
    Ok(OutputEnvelope::new(
        "refine-shift",
        json!({
            "file": a.file, "observed": a.observed, "proforma": a.proforma, "mods": a.mods,
            "sequenceMass": a.sequence_mass, "modMass": a.mod_mass, "modFormula": a.mod_formula,
            "mz": a.mz, "charge": a.charge, "rtMin": a.rt_min, "rtMax": a.rt_max,
            "tolerancePpm": a.tolerance_ppm, "priorShift": a.prior_shift,
            "matchWindowPpm": a.match_window_ppm,
        }),
        sources,
        serde_json::to_value(result).unwrap(),
    ))
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
struct RefineShiftBatchRequest {
    id: String,
    sequence: Option<String>,
    proforma: Option<String>,
    sequence_mass: Option<f64>,
    modification_mass: Option<f64>,
    modification_formula: Option<String>,
    mz: f64,
    charge: i32,
    rt_min: f64,
    rt_max: f64,
    #[serde(default = "default_refine_tolerance")]
    tolerance_ppm: f64,
    #[serde(default)]
    prior_shift: f64,
    #[serde(default = "default_refine_match_window")]
    match_window_ppm: f64,
}

fn default_refine_tolerance() -> f64 {
    10.0
}
fn default_refine_match_window() -> f64 {
    30.0
}

fn refine_shift_batch(a: &RefineShiftBatchArgs) -> Result<OutputEnvelope, CommandError> {
    let (requests, mut sources): (Vec<RefineShiftBatchRequest>, _) =
        load_json_arg("--input", &a.input)?;
    validate_batch_ids(requests.iter().map(|request| request.id.as_str()))?;
    let (bytes, source) = read_file(&a.file)?;
    sources.insert(0, source);
    let mut ids = Vec::with_capacity(requests.len());
    let mut queries = Vec::with_capacity(requests.len());
    let mut early_errors = Vec::with_capacity(requests.len());
    for request in requests {
        ids.push(request.id.clone());
        let sequence_mass = match (&request.sequence, request.sequence_mass) {
            (Some(sequence), _) => match parse_proforma(sequence, &[]) {
                Ok(parsed) => Some(parsed.monoisotopic_mass),
                Err(error) => {
                    early_errors.push(Some((error.code().to_owned(), error.to_string())));
                    queries.push(None);
                    continue;
                }
            },
            (None, mass) => mass,
        };
        early_errors.push(None);
        queries.push(Some(RefineMassShiftQuery {
            // Carry the exact sequence into the query as a bare proforma so
            // `candidate_chemistry` uses the protein's real composition instead
            // of the averagine fallback. `parse_proforma` above already proved
            // it parses; a bare sequence declares no modifications, so the
            // shift mass is applied as a pure delta (refine.rs F3 path).
            proforma: request.proforma.or_else(|| request.sequence.clone()),
            external_mods: Vec::new(),
            sequence_mass,
            modification_mass: request.modification_mass,
            modification_formula: request.modification_formula,
            mz: request.mz,
            charge: request.charge,
            rt_min: request.rt_min,
            rt_max: request.rt_max,
            tolerance_ppm: request.tolerance_ppm,
            prior_shift: request.prior_shift,
            match_window_ppm: request.match_window_ppm,
        }));
    }
    let valid_queries: Vec<_> = queries.iter().filter_map(Clone::clone).collect();
    let mut calibrated = calibrate_mass_shifts(&bytes, &valid_queries).into_iter();
    let results: Vec<Value> = ids
        .into_iter()
        .zip(queries)
        .zip(early_errors)
        .map(|((id, query), early_error)| {
            if let Some((code, message)) = early_error {
                return json!({"id": id, "ok": false, "error": {
                    "code": code, "message": message
                }});
            }
            debug_assert!(query.is_some());
            match calibrated.next().expect("one result per valid query") {
                Ok(result) => json!({"id": id, "ok": true, "result": result}),
                Err(error) => json!({"id": id, "ok": false, "error": {
                    "code": error.code(), "message": error.to_string()
                }}),
            }
        })
        .collect();
    Ok(OutputEnvelope::new(
        "refine-shift-batch",
        json!({"file": a.file, "input": a.input}),
        sources,
        json!({"results": results}),
    ))
}

fn spectrum(a: &SpectrumArgs) -> Result<DispatchOutput, CommandError> {
    if a.out_format != "json" && !(matches!(a.mode, SpectrumMode::Scans) && a.out_format == "mgf") {
        return Err(CommandError::InvalidQuery(format!(
            "unknown --out-format {:?}; expected json{} (tabular exports go through spectra.export, not the kernel)",
            a.out_format,
            if matches!(a.mode, SpectrumMode::Scans) {
                " or mgf"
            } else {
                ""
            }
        )));
    }
    let (bytes, source) = read_file(&a.file)?;
    let tolerance = match (a.tolerance, a.tolerance_unit) {
        (Some(value), ToleranceUnit::Da) => Some(MassTolerance::Da(value)),
        (Some(value), ToleranceUnit::Ppm) => Some(MassTolerance::Ppm(value)),
        (None, _) => None,
    };
    let mut scans: Vec<String> = a
        .scans
        .as_deref()
        .map(|text| text.split(',').map(str::trim).map(str::to_owned).collect())
        .unwrap_or_default();
    if let Some(value) = &a.scans_input {
        scans = parse_scan_list(value)?;
    }

    if a.max_scans.is_some() && !matches!(a.mode, SpectrumMode::Window) {
        return Err(CommandError::InvalidQuery(
            "--max-scans applies to window mode".into(),
        ));
    }
    let scoped = matches!(a.mode, SpectrumMode::Window)
        || a.scan_min.is_some()
        || a.scan_max.is_some()
        || (matches!(
            a.mode,
            SpectrumMode::Sum | SpectrumMode::Scans | SpectrumMode::Product
        ) && (a.rt_min.is_some() || a.rt_max.is_some()))
        || (matches!(a.mode, SpectrumMode::Sum) && a.ms_level > 0)
        || (matches!(a.mode, SpectrumMode::Window) && !scans.is_empty());
    if scoped {
        if !matches!(
            a.mode,
            SpectrumMode::Sum | SpectrumMode::Window | SpectrumMode::Scans | SpectrumMode::Product
        ) {
            return Err(CommandError::InvalidQuery(
                "scan ranges are supported by sum/window/scans/product".into(),
            ));
        }
        let mut reader = tdeasenexus_toolbox::open_reader(&bytes)
            .map_err(|error| mass_spec_error(Some(&a.file), error))?;
        scans = tdeasenexus_toolbox::select_spectrum_ids_reader(
            &mut reader,
            &tdeasenexus_toolbox::SpectrumSelection {
                scans,
                ms_level: if a.ms_level > 0 {
                    Some(a.ms_level)
                } else if matches!(a.mode, SpectrumMode::Window) {
                    Some(1)
                } else if matches!(a.mode, SpectrumMode::Product) {
                    Some(2)
                } else {
                    None
                },
                rt_min: a.rt_min,
                rt_max: a.rt_max,
                scan_min: a.scan_min,
                scan_max: a.scan_max,
            },
        )
        .map_err(|error| mass_spec_error(Some(&a.file), error))?;
    }
    let result: Value = match a.mode {
        SpectrumMode::Window => {
            let truncated = a.max_scans.is_some_and(|limit| scans.len() > limit);
            if a.max_scans == Some(0) {
                return Err(CommandError::InvalidQuery(
                    "maxScans must be positive".into(),
                ));
            }
            if let Some(limit) = a.max_scans {
                scans.truncate(limit);
            }
            let sum = query_sum_mzml(
                &bytes,
                &scans,
                a.bin_ppm,
                a.mz_min,
                a.mz_max,
                a.min_intensity,
                a.max_peaks,
            )
            .map_err(|error| mass_spec_error(Some(&a.file), error))?;
            let mut result = serde_json::to_value(&sum).unwrap();
            result["scansSummed"] = json!(sum.scan_numbers.len());
            result["retentionTimeUnit"] = json!("minute");
            result["truncated"] = json!(truncated || sum.truncated);
            result
        }
        SpectrumMode::Sum => {
            let sum = query_sum_mzml(
                &bytes,
                &scans,
                a.bin_ppm,
                a.mz_min,
                a.mz_max,
                a.min_intensity,
                a.max_peaks,
            )
            .map_err(|error| mass_spec_error(Some(&a.file), error))?;
            json!(sum)
        }
        SpectrumMode::Xic => {
            let mz = a.mz.ok_or_else(|| {
                CommandError::MissingArgument("spectrum --mode xic requires --mz".to_owned())
            })?;
            let tolerance_ppm = a.tolerance.unwrap_or(20.0);
            json!(query_xic_mzml(
                &bytes,
                mz,
                tolerance_ppm,
                a.ms_level.max(1),
                a.rt_min,
                a.rt_max
            )
            .map_err(|error| mass_spec_error(Some(&a.file), error))?)
        }
        SpectrumMode::Scans => {
            let peak_order = match a.peak_order.as_deref() {
                Some(text) => parse_peak_order(text)?,
                None => PeakOrder::MzAscending,
            };
            let filters = ScanQuery {
                selector: SpectrumSelector::Index(0),
                mz_min: a.mz_min,
                mz_max: a.mz_max,
                min_intensity: a.min_intensity,
                max_peaks: a.max_peaks,
                peak_order,
                // 0 = any level (the level-agnostic default); an explicit
                // --ms-level filters and mismatches land in missingScans.
                ms_level: (a.ms_level != 0).then_some(a.ms_level),
            };
            let scans_result = query_scans_mzml(&bytes, &scans, &filters)
                .map_err(|error| mass_spec_error(Some(&a.file), error))?;
            if a.out_format == "mgf" {
                return scans_to_mgf(&a.file, &scans_result, a.charge);
            }
            json!(scans_result)
        }
        SpectrumMode::Precursor => {
            let mz = a.mz.ok_or_else(|| {
                CommandError::MissingArgument("spectrum --mode precursor requires --mz".to_owned())
            })?;
            let tolerance = tolerance.ok_or_else(|| {
                CommandError::MissingArgument(
                    "spectrum --mode precursor requires --tolerance".to_owned(),
                )
            })?;
            json!(query_precursor_mzml(
                &bytes,
                &PrecursorQuery {
                    mz,
                    tolerance,
                    rt_min: a.rt_min,
                    rt_max: a.rt_max,
                    charge: a.charge,
                    max_matches: a.max_peaks,
                },
            )
            .map_err(|error| mass_spec_error(Some(&a.file), error))?)
        }
        SpectrumMode::Product => {
            let mz = a.mz.ok_or_else(|| {
                CommandError::MissingArgument("spectrum --mode product requires --mz".to_owned())
            })?;
            let tolerance = tolerance.ok_or_else(|| {
                CommandError::MissingArgument(
                    "spectrum --mode product requires --tolerance".to_owned(),
                )
            })?;
            // --scans restricts the census to the named spectra; without it
            // the query stays the whole-file survey.
            let scan_filter: Option<Vec<String>> = a.scans.as_deref().map(|text| {
                text.split(',')
                    .map(str::trim)
                    .filter(|scan| !scan.is_empty())
                    .map(str::to_owned)
                    .collect::<Vec<String>>()
            });
            json!(query_product_mzml(
                &bytes,
                &ProductQuery {
                    mz,
                    tolerance,
                    max_matches: a.max_peaks,
                    scans: scan_filter,
                },
            )
            .map_err(|error| mass_spec_error(Some(&a.file), error))?)
        }
    };

    Ok(DispatchOutput::Json(OutputEnvelope::new(
        "spectrum",
        json!({
            "file": a.file, "mode": format!("{:?}", a.mode).to_lowercase(),
            "msLevel": a.ms_level, "rtMin": a.rt_min, "rtMax": a.rt_max,
            "maxScans": a.max_scans, "scanMin": a.scan_min, "scanMax": a.scan_max, "binPpm": a.bin_ppm, "mz": a.mz,
            "tolerance": a.tolerance, "toleranceUnit": a.tolerance_unit.as_str(), "charge": a.charge,
            "scans": a.scans, "mzMin": a.mz_min, "mzMax": a.mz_max,
            "minIntensity": a.min_intensity, "maxPeaks": a.max_peaks,
            "outFormat": a.out_format,
        }),
        vec![source],
        result,
    )))
}

fn select_cmd(a: &SelectArgs) -> Result<OutputEnvelope, CommandError> {
    let (bytes, _) = read_file(&a.file)?;
    let result = query_select_mzml(&bytes, &a.statement)
        .map_err(|error| mass_spec_error(Some(&a.file), error))?;
    let source = SourceInfo {
        path: a.file.clone(),
        size: bytes.len() as u64,
        digest: SourceDigest {
            algo: "sha256:head1MiB".to_owned(),
            value: head_digest(&bytes),
        },
    };
    Ok(OutputEnvelope::new(
        "select",
        json!({ "file": a.file, "statement": a.statement }),
        vec![source],
        serde_json::to_value(result).unwrap(),
    ))
}

/// Parse a comma-separated ion-type list (`B,Y,C,Z`) into core `IonType`s.
fn parse_ion_types(text: &str) -> Result<Vec<IonType>, CommandError> {
    text.split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|token| match token.to_ascii_uppercase().as_str() {
            "A" => Ok(IonType::A),
            "B" => Ok(IonType::B),
            "C" => Ok(IonType::C),
            "X" => Ok(IonType::X),
            "Y" => Ok(IonType::Y),
            "Z" => Ok(IonType::Z),
            "ZDOT" | "ZPRIME" => Err(CommandError::InvalidQuery(
                "ZDot/ZPrime are not generated: theoretical fragment series are decided \
                 by the fragmentation method (ETD -> c/z, HCD/CID -> b/y); use Z"
                    .to_owned(),
            )),
            _ => Err(CommandError::InvalidQuery(format!(
                "unknown ion type {token:?}; expected A/B/C/X/Y/Z"
            ))),
        })
        .collect()
}

/// Load a JSON argument once, retaining provenance when it names a file.
fn load_json_arg<T: serde::de::DeserializeOwned>(
    flag: &str,
    text: &str,
) -> Result<(T, Vec<SourceInfo>), CommandError> {
    if Path::new(text).is_file() {
        let (bytes, source) = read_file(text)?;
        let value = serde_json::from_slice(&bytes).map_err(|error| CommandError::JsonParse {
            path: text.into(),
            error: error.to_string(),
        })?;
        Ok((value, vec![source]))
    } else {
        let value = serde_json::from_str(text).map_err(|error| CommandError::JsonParse {
            path: flag.into(),
            error: error.to_string(),
        })?;
        Ok((value, Vec::new()))
    }
}

fn parse_json_arg<T: serde::de::DeserializeOwned>(
    flag: &str,
    text: &str,
) -> Result<T, CommandError> {
    load_json_arg(flag, text).map(|(value, _)| value)
}

fn validate_batch_ids<'a>(ids: impl IntoIterator<Item = &'a str>) -> Result<(), CommandError> {
    let mut seen = std::collections::HashSet::new();
    for id in ids {
        if id.trim().is_empty() || !seen.insert(id) {
            return Err(CommandError::InvalidQuery(
                "batch ids must be nonempty and unique".into(),
            ));
        }
    }
    if seen.is_empty() {
        return Err(CommandError::InvalidQuery(
            "batch input cannot be empty".into(),
        ));
    }
    Ok(())
}

/// A file or inline JSON array; integers and native-id strings share one
/// normalization path. Reject empty, fractional and misspelled structures.
fn parse_scan_list(value: &str) -> Result<Vec<String>, CommandError> {
    let values: Vec<Value> = parse_json_arg("--scans-input", value)?;
    if values.is_empty() {
        return Err(CommandError::InvalidQuery(
            "--scans-input cannot be empty".into(),
        ));
    }
    values
        .into_iter()
        .enumerate()
        .map(|(index, value)| match value {
            Value::String(id) if !id.trim().is_empty() => Ok(id.trim().to_owned()),
            Value::Number(number) if number.as_u64().is_some() => Ok(number.to_string()),
            _ => Err(CommandError::InvalidQuery(format!(
                "--scans-input[{index}] must be a scan integer or native-id string"
            ))),
        })
        .collect()
}

/// `--neutral-losses` resolution: the sentinel `engine` turns on the
/// fragmentation engine's own loss variants; anything else is the JSON array
/// of site-scoped named losses.
fn parse_neutral_losses(text: &str) -> Result<(Vec<NamedNeutralLoss>, bool), CommandError> {
    if text.trim() == "engine" {
        Ok((Vec::new(), true))
    } else {
        let named: Vec<NamedNeutralLoss> = parse_json_arg("--neutral-losses", text)?;
        Ok((named, false))
    }
}

#[derive(serde::Serialize, Deserialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
struct CandidateProforma {
    id: String,
    proforma: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    label: Option<String>,
}

/// Evidence-lattice rung of a parsed hypothesis: which class of claim its
/// tags can support. `massDelta` (any numeric tag) supports ladder/position
/// measurements only; `composition` (formula tags) adds isotope-shape
/// claims; `named` tags additionally carry a candidate identity. The rung
/// travels with every fragment-family result so downstream adjudication
/// knows exactly what was tested.
fn hypothesis_rung(parsed: &ParsedProForma) -> &'static str {
    if !parsed.warnings.is_empty() {
        "massDelta"
    } else if parsed
        .modifications
        .iter()
        .any(|m| m.tag.starts_with("Formula:"))
    {
        "composition"
    } else if parsed.modifications.is_empty() {
        "unmodified"
    } else {
        "named"
    }
}

fn rung_rank(rung: &str) -> u8 {
    match rung {
        "massDelta" => 0,
        "composition" => 1,
        _ => 2,
    }
}

/// Diagnostic for a candidate whose parsed backbone differs from the expected
/// sequence: residue counts and the first divergent position with local
/// context, so an inserted/deleted residue in a hand-built model is visible
/// immediately instead of silently invalidating the comparison.
fn backbone_diff_error(candidate_id: &str, actual: &str, expected: &str) -> String {
    let first_diff = actual
        .char_indices()
        .zip(expected.char_indices())
        .find(|((_, a), (_, b))| a != b)
        .map(|((i, a), (_, b))| (i, a, b));
    let context = |sequence: &str, at: usize| -> String {
        let start = at.saturating_sub(6);
        let end = (at + 7).min(sequence.len());
        sequence[start..end].to_owned()
    };
    let detail = match first_diff {
        Some((index, a, b)) => format!(
            "first divergence at position {} ('{}' vs '{}', context '{}' vs '{}')",
            index + 1,
            a,
            b,
            context(actual, index),
            context(expected, index)
        ),
        None => format!(
            "one sequence is a prefix of the other (lengths {} vs {})",
            actual.len(),
            expected.len()
        ),
    };
    format!(
        "candidate '{candidate_id}' does not match the expected backbone: {} residues vs {}; {detail}. \
         Rebuild the model against the source-declared sequence and pass it via --sequence",
        actual.chars().count(),
        expected.chars().count()
    )
}

fn fragment_model(
    proforma: &str,
) -> Result<(String, Vec<MassModification>, &'static str), CommandError> {
    let parsed = parse_proforma(proforma, &[])
        .map_err(|error| CommandError::ProformaParse(error.to_string()))?;
    // Evidence lattice: numeric mass-only tags are legal here — fragment m/z
    // depends only on the summed deltas, never on the modification's
    // elements. The rung labels the output; matching never claims shape.
    let rung = hypothesis_rung(&parsed);
    let modifications = parsed
        .modifications
        .into_iter()
        .map(|modification| {
            let position = modification.position.ok_or_else(|| {
                CommandError::InvalidQuery(format!(
                    "ProForma modification [{}] has no exact site",
                    modification.tag
                ))
            })?;
            if position.starts_with("site:") {
                return Err(CommandError::InvalidQuery(format!(
                    "ProForma modification [{}] has residue-only, not positional, localization",
                    modification.tag
                )));
            }
            if modification.localization_start != modification.localization_end {
                return Err(CommandError::InvalidQuery(format!(
                    "ProForma modification [{}] has an unresolved localization range",
                    modification.tag
                )));
            }
            let position = position
                .strip_prefix("residue:")
                .unwrap_or(&position)
                .to_owned();
            Ok(MassModification {
                position,
                delta_mass: modification.delta_mass,
                formula: modification.formula,
            })
        })
        .collect::<Result<Vec<_>, CommandError>>()?;
    Ok((parsed.sequence, modifications, rung))
}

// ── envelope-regions: selective observed extraction inside theoretical envelopes ──

/// One ion selection: `y50`, `y50:6`, or a range `b4-b7` (charge-less).
fn parse_ion_selection(spec: &str) -> Result<Vec<(char, usize, Option<u8>)>, CommandError> {
    let mut out = Vec::new();
    for token in spec.split(',').map(str::trim).filter(|t| !t.is_empty()) {
        let bad = || {
            CommandError::InvalidQuery(format!(
                "cannot parse ion selector {token:?}: expected y50, y50:6, or b4-b7"
            ))
        };
        let (range, charge) = match token.split_once(':') {
            Some((r, z)) => (r, Some(z.parse::<u8>().map_err(|_| bad())?)),
            None => (token, None),
        };
        let (head, numbers) = range.split_at(1);
        let ion = head.to_ascii_lowercase().chars().next().ok_or_else(bad)?;
        if !matches!(ion, 'a' | 'b' | 'c' | 'x' | 'y' | 'z') {
            return Err(bad());
        }
        let nums: Vec<usize> = if let Some((lo, hi)) = numbers.split_once('-') {
            let (lo, hi): (usize, usize) = (
                lo.parse().map_err(|_| bad())?,
                hi.parse().map_err(|_| bad())?,
            );
            if hi < lo || hi - lo > 200 {
                return Err(bad());
            }
            (lo..=hi).collect()
        } else {
            vec![numbers.parse().map_err(|_| bad())?]
        };
        for number in nums {
            out.push((ion, number, charge));
        }
    }
    if out.is_empty() {
        return Err(CommandError::InvalidQuery("--ions is empty".to_owned()));
    }
    Ok(out)
}

fn envelope_regions(a: &crate::args::EnvelopeRegionsArgs) -> Result<OutputEnvelope, CommandError> {
    let (bytes, mzml_source) = read_file(&a.file)?;
    let selection = parse_ion_selection(&a.ions)?;
    let scan_selector = |scan: &str| {
        if !scan.is_empty() && scan.chars().all(|c| c.is_ascii_digit()) {
            tdeasenexus_toolbox::SpectrumSelector::Id(format!("scan={scan}"))
        } else {
            tdeasenexus_toolbox::SpectrumSelector::Id(scan.to_owned())
        }
    };

    // Charge ceiling mirrors `fragment`: explicit --charge wins, else the
    // precursor's declared charge (the physical ceiling).
    let mut reader = tdeasenexus_toolbox::open_reader(&bytes)
        .map_err(|error| mass_spec_error(Some(&a.file), error))?;
    let (max_charge, max_charge_source) = match a.charge {
        Some(value) => (value, "explicit"),
        None => match tdeasenexus_toolbox::fragment_precursor_charge(
            &mut reader,
            &scan_selector(&a.scan),
        )
        .map_err(|error| mass_spec_error(Some(&a.file), error))?
        {
            Some(z) => (z, "precursor"),
            None => {
                eprintln!(
                    "WARNING envelope-regions: no precursor charge declared for {} scan {}; charge ceiling falls back to 1 — pass --charge explicitly.",
                    a.file, a.scan
                );
                (1, "fallback")
            }
        },
    };

    let ion_types: Vec<tdeasenexus_toolbox::IonType> = {
        let mut chars: Vec<char> = selection.iter().map(|(ion, _, _)| *ion).collect();
        chars.sort();
        chars.dedup();
        chars
            .into_iter()
            .map(|c| match c {
                'a' => tdeasenexus_toolbox::IonType::A,
                'b' => tdeasenexus_toolbox::IonType::B,
                'c' => tdeasenexus_toolbox::IonType::C,
                'x' => tdeasenexus_toolbox::IonType::X,
                'z' => tdeasenexus_toolbox::IonType::Z,
                _ => tdeasenexus_toolbox::IonType::Y,
            })
            .collect()
    };
    let fragments = tdeasenexus_toolbox::theoretical_fragment_list_proforma_variants(
        std::slice::from_ref(&a.proforma),
        &[],
        &ion_types,
        max_charge,
        &[],
        false,
        tdeasenexus_toolbox::EnvelopeMode::Exact,
    )
    .map_err(|error| CommandError::InvalidQuery(error.to_string()))?;

    // One full MS2 peak list for the scan; each region filters its own window.
    let filters = tdeasenexus_toolbox::ScanQuery {
        selector: tdeasenexus_toolbox::SpectrumSelector::Index(0),
        mz_min: None,
        mz_max: None,
        min_intensity: None,
        max_peaks: None,
        peak_order: tdeasenexus_toolbox::PeakOrder::MzAscending,
        ms_level: Some(2),
    };
    let scans =
        tdeasenexus_toolbox::query_scans_mzml(&bytes, std::slice::from_ref(&a.scan), &filters)
            .map_err(|error| mass_spec_error(Some(&a.file), error))?;
    let spectrum = scans
        .spectra
        .first()
        .ok_or_else(|| CommandError::InvalidQuery(format!("scan {} not found", a.scan)))?;
    let all_peaks: Vec<(f64, f32)> = spectrum
        .spectrum
        .peaks
        .iter()
        .map(|peak| (peak.mz, peak.intensity))
        .collect();
    drop(reader);

    let mut regions: Vec<serde_json::Value> = Vec::new();
    let mut empty_regions: Vec<String> = Vec::new();
    // A requested selection that yields no region (no theoretical fragment at
    // this charge, or an empty isotope envelope) must be stated, never
    // silently dropped: exit-0-with-nothing looks like "the ion does not
    // exist" when it may be a one-off envelope-generation failure.
    let mut missing_regions: Vec<String> = Vec::new();
    for (ion, number, wanted_charge) in &selection {
        let before = regions.len();
        for fragment in fragments.iter().filter(|f| {
            let ion_char = match f.ion {
                tdeasenexus_toolbox::IonType::A => 'a',
                tdeasenexus_toolbox::IonType::B => 'b',
                tdeasenexus_toolbox::IonType::C => 'c',
                tdeasenexus_toolbox::IonType::X => 'x',
                tdeasenexus_toolbox::IonType::Y => 'y',
                tdeasenexus_toolbox::IonType::Z => 'z',
            };
            ion_char == *ion && f.number == *number
        }) {
            let charge = fragment.charge;
            if let Some(want) = wanted_charge {
                if charge != *want {
                    continue;
                }
            }
            if fragment.isotope_envelope.is_empty() {
                continue;
            }
            let dz = tdeasenexus_toolbox::NEUTRON_SHIFT / f64::from(charge);
            let pad = f64::from(a.pad) * dz;
            let first = fragment
                .isotope_envelope
                .first()
                .map(|p| p.mz)
                .unwrap_or(fragment.mz);
            let last = fragment
                .isotope_envelope
                .last()
                .map(|p| p.mz)
                .unwrap_or(fragment.mz);
            let lo = first - pad;
            let hi = last + pad;
            // Pure extraction: every peak inside the window is returned,
            // annotated with its nearest theoretical slot. Whether a peak
            // "matches" is the caller's prior judgment (fragment match), not
            // this command's — no gating, no selection.
            let mut observed: Vec<serde_json::Value> = Vec::new();
            for (mz, intensity) in all_peaks.iter().filter(|(mz, _)| *mz >= lo && *mz <= hi) {
                let mut nearest: Option<(usize, f64)> = None;
                for (slot, stick) in fragment.isotope_envelope.iter().enumerate() {
                    let err = (mz - stick.mz) * 1.0e6 / stick.mz;
                    if nearest.is_none_or(|(_, best)| err.abs() < best.abs()) {
                        nearest = Some((slot, err));
                    }
                }
                if let Some((slot, err)) = nearest {
                    observed.push(serde_json::json!({
                        "mz": mz, "intensity": intensity,
                        "nearestNeutronOffset": fragment.isotope_envelope[slot].neutron_offset,
                        "errPpm": err,
                    }));
                }
            }
            if observed.is_empty() {
                empty_regions.push(format!("{ion}{number}:z{charge}"));
            }
            regions.push(serde_json::json!({
                "ion": ion.to_string(), "number": number, "charge": charge,
                "neutralMass": fragment.neutral_mass,
                "monoMz": fragment.mz, "apexMz": fragment.apex_mz, "apexOffset": fragment.apex_offset,
                "theoretical": fragment.isotope_envelope.iter().map(|stick| serde_json::json!({
                    "neutronOffset": stick.neutron_offset,
                    "mz": stick.mz,
                    "relativeIntensity": stick.relative_intensity,
                })).collect::<Vec<_>>(),
                "observed": observed,
                "window": {"mzMin": lo, "mzMax": hi, "padSlots": a.pad},
            }));
        }
        if regions.len() == before {
            let wanted = wanted_charge.map(|z| format!(":z{z}")).unwrap_or_default();
            missing_regions.push(format!("{ion}{number}{wanted}"));
        }
    }
    if !missing_regions.is_empty() {
        eprintln!(
            "WARNING envelope-regions: requested ion(s) produced no region (no theoretical fragment at this charge, or an empty isotope envelope): {}",
            missing_regions.join(", ")
        );
    }

    Ok(OutputEnvelope::new(
        "envelope-regions",
        serde_json::json!({
            "file": a.file, "scan": a.scan, "proforma": a.proforma,
            "ions": a.ions, "charge": a.charge,
            "maxCharge": max_charge, "maxChargeSource": max_charge_source,
            "pad": a.pad,
        }),
        vec![mzml_source],
        serde_json::json!({
            "regions": regions,
            "requested": selection.len(),
            "extracted": regions.len(),
            "emptyRegions": empty_regions,
            "missingRegions": missing_regions,
        }),
    ))
}

fn fragment(a: &FragmentArgs) -> Result<OutputEnvelope, CommandError> {
    let requested_ions = a.ion_types.as_deref().map(parse_ion_types).transpose()?;
    tdeasenexus_toolbox::validate_fragment_parameters(
        a.max_charge,
        a.tolerance_ppm,
        a.min_correlation,
        a.min_intensity,
        a.max_matches,
    )
    .map_err(|error| mass_spec_error(Some(&a.file), error))?;
    if !a.bin_ppm.is_finite()
        || a.bin_ppm < 0.0
        || !a.minimum_abundance_pct.is_finite()
        || !(0.0..=100.0).contains(&a.minimum_abundance_pct)
    {
        return Err(CommandError::InvalidQuery(
            "binPpm must be non-negative and minimumAbundancePct in 0..100".into(),
        ));
    }

    let mode = a.mode.unwrap_or(if a.proforma.is_some() {
        FragmentMode::Match
    } else {
        FragmentMode::Compare
    });
    let models: Vec<CandidateProforma> = match mode {
        FragmentMode::Match => {
            if a.candidate_proformas.is_some() {
                return Err(CommandError::InvalidQuery(
                    "match accepts --proforma, not --candidate-proformas".into(),
                ));
            }
            vec![CandidateProforma {
                id: "model".into(),
                label: None,
                proforma: a.proforma.clone().ok_or_else(|| {
                    CommandError::MissingArgument(
                        "fragment --mode match requires --proforma".into(),
                    )
                })?,
            }]
        }
        FragmentMode::Compare => {
            if a.proforma.is_some() {
                return Err(CommandError::InvalidQuery(
                    "compare accepts --candidate-proformas, not --proforma".into(),
                ));
            }
            parse_json_arg(
                "--candidate-proformas",
                a.candidate_proformas.as_deref().ok_or_else(|| {
                    CommandError::MissingArgument(
                        "fragment requires --proforma or --candidate-proformas".into(),
                    )
                })?,
            )?
        }
    };
    if models.is_empty() {
        return Err(CommandError::InvalidQuery(
            "at least one complete candidate is required".into(),
        ));
    }
    let mut parsed = Vec::with_capacity(models.len());
    let mut hypothesis_type = "unmodified";
    let mut ids = std::collections::HashSet::new();
    for model in &models {
        if model.id.trim().is_empty() || !ids.insert(&model.id) {
            return Err(CommandError::InvalidQuery(
                "model ids must be nonempty and unique".into(),
            ));
        }
        let (sequence, modifications, rung) = fragment_model(&model.proforma)?;
        if let Some(expected) = &a.sequence {
            let expected = expected.trim().to_ascii_uppercase();
            if sequence != expected {
                return Err(CommandError::InvalidQuery(backbone_diff_error(
                    &model.id, &sequence, &expected,
                )));
            }
        }
        if rung_rank(rung) < rung_rank(hypothesis_type) || hypothesis_type == "unmodified" {
            hypothesis_type = rung;
        }
        parsed.push((sequence, modifications));
    }
    let shared_backbone = parsed.iter().all(|model| model.0 == parsed[0].0);
    let candidate_models: Vec<_> = models
        .iter()
        .zip(&parsed)
        .map(|(model, (sequence, modifications))| CandidateModel {
            id: model.id.clone(),
            sequence: (!shared_backbone).then(|| sequence.clone()),
            modifications: modifications.clone(),
        })
        .collect();
    let (neutral_losses, engine_losses) = a
        .neutral_losses
        .as_deref()
        .map(parse_neutral_losses)
        .transpose()?
        .unwrap_or_default();
    let evidence_path = a.evidence_out.as_deref().ok_or_else(|| {
        CommandError::MissingArgument(
            "fragment matching requires --evidence-out; full match tables are file-backed".into(),
        )
    })?;
    if let Ok(destination) = Path::new(evidence_path).canonicalize() {
        for input in [Some(a.file.as_str()), a.candidate_proformas.as_deref()]
            .into_iter()
            .flatten()
        {
            if Path::new(input)
                .canonicalize()
                .is_ok_and(|source| source == destination)
            {
                return Err(CommandError::InvalidQuery(
                    "evidence output must not overwrite an input file".into(),
                ));
            }
        }
    }
    let (bytes, source) = read_file(&a.file)?;
    let sources = vec![source];
    let is_json = bytes.iter().copied().find(|b| !b.is_ascii_whitespace()) == Some(b'{');
    let mut observation_source = Value::Null;
    let mut observation_parameters = Value::Null;
    let mut saved: Option<tdeasenexus_toolbox::SpectrumObservation> = if is_json {
        let doc: Value =
            serde_json::from_slice(&bytes).map_err(|error| CommandError::JsonParse {
                path: a.file.clone(),
                error: error.to_string(),
            })?;
        let payload = if doc.get("schemaVersion").is_some() {
            if doc["schemaVersion"] != 1
                || doc["tool"] != "tdeasenexus-toolbox"
                || doc["command"] != "spectrum"
            {
                return Err(CommandError::InvalidQuery(
                    "expected a schemaVersion 1 spectrum observation artifact".into(),
                ));
            }
            observation_source = doc.get("source").cloned().unwrap_or(Value::Null);
            observation_parameters = doc.get("params").cloned().unwrap_or(Value::Null);
            doc["result"].clone()
        } else {
            doc
        };
        Some(serde_json::from_value(payload).map_err(|error| {
            CommandError::InvalidQuery(format!(
                "invalid spectrum observation: {error}; use spectrum --mode sum/window"
            ))
        })?)
    } else {
        None
    };
    let has_saved = saved.is_some();
    // Content determines the observation domain; a renamed msalign remains
    // deconvoluted and can never enter the raw envelope matcher.
    let is_msalign = !is_json
        && bytes
            .windows(b"BEGIN IONS".len())
            .any(|w| w == b"BEGIN IONS")
        && bytes
            .windows(b"PRECURSOR_MASS=".len())
            .any(|w| w == b"PRECURSOR_MASS=");
    let mut scans: Vec<String> = if let Some(value) = &a.scans_input {
        parse_scan_list(value)?
    } else {
        a.scan
            .split(',')
            .map(str::trim)
            .filter(|id| !id.is_empty())
            .map(str::to_owned)
            .collect()
    };
    let scoped =
        a.rt_min.is_some() || a.rt_max.is_some() || a.scan_min.is_some() || a.scan_max.is_some();
    let mut reader = if saved.is_none() && !is_msalign {
        Some(
            tdeasenexus_toolbox::open_reader(&bytes)
                .map_err(|error| mass_spec_error(Some(&a.file), error))?,
        )
    } else {
        None
    };
    if let Some(observation) = &saved {
        if !scans.is_empty() || scoped || a.sum_scans {
            return Err(CommandError::InvalidQuery(
            "saved observations already declare their scans and aggregation; selection applies to source mzML".into()));
        }
        scans = observation
            .spectra_metadata
            .iter()
            .map(|m| m.id.clone())
            .collect();
    } else if is_msalign {
        if scoped || a.sum_scans {
            return Err(CommandError::InvalidQuery(
            "msalign supports scan enumeration; raw-spectrum pooling and RT selection require mzML".into()));
        }
        if scans.is_empty() {
            return Err(CommandError::MissingArgument(
                "msalign requires --scan".into(),
            ));
        }
        let mut unique = std::collections::HashSet::new();
        if scans
            .iter()
            .any(|s| !unique.insert(s.strip_prefix("scan=").unwrap_or(s)))
        {
            return Err(CommandError::InvalidQuery(
                "duplicate scan in selection".into(),
            ));
        }
    } else {
        if scans.is_empty() && !scoped {
            return Err(CommandError::MissingArgument(
                "fragment requires scan ids, RT bounds, or scan-number bounds".into(),
            ));
        }
        let selected = tdeasenexus_toolbox::select_spectrum_ids_reader(
            reader.as_mut().unwrap(),
            &tdeasenexus_toolbox::SpectrumSelection {
                scans: scans.clone(),
                ms_level: Some(2),
                rt_min: a.rt_min,
                rt_max: a.rt_max,
                scan_min: a.scan_min,
                scan_max: a.scan_max,
            },
        )
        .map_err(|error| mass_spec_error(Some(&a.file), error))?;
        if scoped {
            scans = selected;
        }
    }
    let single_match =
        matches!(mode, FragmentMode::Match) && (scans.len() == 1 || a.sum_scans || has_saved);
    let mut params = json!({
        "file": a.file, "mode": if matches!(mode, FragmentMode::Match) { "match" } else { "compare" },
        "scan": a.scan, "scansInput": a.scans_input, "scans": scans, "rtMin": a.rt_min, "rtMax": a.rt_max,
        "scanMin": a.scan_min, "scanMax": a.scan_max,
        "proforma": a.proforma, "candidateProformas": a.candidate_proformas, "models": models,
        "evidenceOut": evidence_path, "activation": a.activation,
        "maxCharge": a.max_charge, "tolerancePpm": a.tolerance_ppm,
        "minCorrelation": a.min_correlation, "minIntensity": a.min_intensity,
        "maxMatches": a.max_matches, "neutralLosses": a.neutral_losses,
        "observationDomain": if is_msalign { "deconvolutedNeutralMass" } else { "centroidMz" },
        "observationSource": observation_source, "observationParameters": observation_parameters,
        "observationAggregation": saved.as_ref().map(|o| json!({"binPpm": o.bin_ppm, "minimumAbundancePct": o.minimum_abundance_pct})),
        "sumScans": a.sum_scans, "binPpm": a.bin_ppm, "minimumAbundancePct": a.minimum_abundance_pct,
    });
    let mut series_warnings = Vec::new();
    let mut evidence_comparisons = Vec::new();
    let mut stdout_comparisons = Vec::new();
    let mut single = None;
    let groups = if a.sum_scans || has_saved {
        vec![scans.clone()]
    } else {
        scans.iter().map(|scan| vec![scan.clone()]).collect()
    };
    let mut prepared_matches =
        std::collections::HashMap::<(String, u8), tdeasenexus_toolbox::PreparedFragmentMatch>::new(
        );
    let mut prepared_comparisons = std::collections::HashMap::<
        (String, u8),
        tdeasenexus_toolbox::PreparedCandidateComparison,
    >::new();
    for group in groups {
        let scan = &group[0];
        let selector = SpectrumSelector::Id(if scan.chars().all(|c| c.is_ascii_digit()) {
            format!("scan={scan}")
        } else {
            scan.clone()
        });
        let mut observation = if has_saved {
            saved.take()
        } else if is_msalign {
            None
        } else {
            let selectors: Vec<_> = group
                .iter()
                .map(|scan| {
                    SpectrumSelector::Id(if scan.chars().all(|c| c.is_ascii_digit()) {
                        format!("scan={scan}")
                    } else {
                        scan.clone()
                    })
                })
                .collect();
            Some(
                tdeasenexus_toolbox::extract_spectrum_observation_reader(
                    reader.as_mut().unwrap(),
                    &selectors,
                    a.sum_scans.then_some(a.bin_ppm),
                    if a.sum_scans {
                        a.minimum_abundance_pct
                    } else {
                        0.0
                    },
                )
                .map_err(|error| mass_spec_error(Some(&a.file), error))?,
            )
        };
        let deconvoluted = if is_msalign {
            Some(
                tdeasenexus_toolbox::parse_msalign_spectrum(
                    &bytes,
                    scan.strip_prefix("scan=").unwrap_or(scan),
                )
                .map_err(|error| mass_spec_error(Some(&a.file), error))?,
            )
        } else {
            None
        };
        let activation = if let Some(observation) = &observation {
            observation
                .fragment_activation(a.activation.as_deref())
                .map_err(|error| mass_spec_error(Some(&a.file), error))?
        } else {
            a.activation
                .clone()
                .or_else(|| deconvoluted.as_ref().and_then(|s| s.activation.clone()))
        };
        let ion_types = if let Some(ions) = &requested_ions {
            ions.clone()
        } else {
            ion_types_for_activation(activation.as_deref())
                .map_err(|error| mass_spec_error(Some(&a.file), error))?
        };
        let ion_types_text = ion_types
            .iter()
            .map(|ion| ion_label(*ion).to_uppercase())
            .collect::<Vec<_>>()
            .join(",");
        let declared = observation
            .as_ref()
            .and_then(|o| {
                o.spectra_metadata
                    .iter()
                    .filter_map(|m| {
                        m.precursor
                            .as_ref()
                            .and_then(|p| p.charge)
                            .and_then(|z| u8::try_from(z).ok())
                            .filter(|z| *z > 0)
                    })
                    .max()
            })
            .or_else(|| deconvoluted.as_ref().and_then(|s| s.precursor_charge));
        let (max_charge, max_charge_source) = match (a.max_charge, declared) {
            (Some(z), _) => (z, "explicit"),
            (None, Some(z)) => (z, "precursor"),
            _ => {
                series_warnings.push(format!("scan {scan}: no precursor charge; maxCharge falls back to 1; supply --max-charge for another ceiling"));
                (1, "fallback")
            }
        };
        if let Some(observation) = &mut observation {
            observation
                .peaks
                .retain(|p| a.min_intensity.is_none_or(|min| p.intensity >= min));
        }
        if scans.len() == 1 || a.sum_scans || has_saved {
            params["resolvedActivation"] = json!(activation);
            params["ionTypes"] = json!(ion_types_text);
            params["maxCharge"] = json!(max_charge);
            params["maxChargeSource"] = json!(max_charge_source);
        }
        let is_pool = observation
            .as_ref()
            .is_some_and(|o| o.spectra_metadata.len() > 1)
            || a.sum_scans;
        let metadata = json!({"scan": scan, "activation": activation, "ionTypes": ion_types,
            "maxCharge": max_charge, "maxChargeSource": max_charge_source,
            "spectraMetadata": observation.as_ref().map(|o| &o.spectra_metadata),
            "pooledScans": if is_pool { Some(&group) } else { None },
            "binPpm": observation.as_ref().and_then(|o| o.bin_ppm),
            "minimumAbundancePct": observation.as_ref().map(|o| o.minimum_abundance_pct),
            "pooledPeaks": if is_pool { observation.as_ref().map(|o| &o.peaks) } else { None },
            "observedPeaks": if !is_pool { observation.as_ref().map(|o| &o.peaks) } else { None },
        });
        if single_match {
            let key = (ion_types_text.clone(), max_charge);
            let prepared = match prepared_matches.entry(key) {
                std::collections::hash_map::Entry::Occupied(entry) => entry.into_mut(),
                std::collections::hash_map::Entry::Vacant(entry) => {
                    let query = FragmentMatchQuery {
                        selector,
                        sequence: parsed[0].0.clone(),
                        modifications: parsed[0].1.clone(),
                        ion_types,
                        max_charge,
                        tolerance_ppm: a.tolerance_ppm,
                        min_intensity: a.min_intensity,
                        max_matches: a.max_matches,
                        neutral_losses: neutral_losses.clone(),
                        engine_neutral_losses: engine_losses,
                        min_correlation: a.min_correlation,
                    };
                    entry.insert(
                        tdeasenexus_toolbox::PreparedFragmentMatch::new(&query)
                            .map_err(|error| mass_spec_error(Some(&a.file), error))?,
                    )
                }
            };
            let matched = if let Some(observation) = &observation {
                prepared.match_observation(observation, a.activation.as_deref())
            } else {
                prepared.match_deconvoluted(deconvoluted.as_ref().unwrap())
            }
            .map_err(|error| mass_spec_error(Some(&a.file), error))?;
            if let Some(d) = matched.series_divergence.as_ref().filter(|d| d.flagged) {
                series_warnings.push(format!("series error divergence: N-terminal minus C-terminal signed mean {:+.2} ppm exceeds {:.2} ppm", d.n_minus_c_ppm, d.threshold_ppm));
            }
            let matched = serde_json::to_value(matched).unwrap();
            single = Some(
                json!({"measurementRole": "candidate_model_evaluation", "notIdentity": true,
                "hypothesisType": hypothesis_type, "query": {"proforma": models[0].proforma,
                    "sequence": parsed[0].0, "modifications": parsed[0].1, "ionTypes": ion_types_text,
                    "maxCharge": max_charge, "maxChargeSource": max_charge_source, "tolerancePpm": a.tolerance_ppm,
                    "neutralLosses": a.neutral_losses}, "observation": metadata, "measurement": matched}),
            );
        } else {
            let key = (ion_types_text.clone(), max_charge);
            let prepared = match prepared_comparisons.entry(key) {
                std::collections::hash_map::Entry::Occupied(entry) => entry.into_mut(),
                std::collections::hash_map::Entry::Vacant(entry) => {
                    let query = CandidateComparisonQuery {
                        max_matches: a.max_matches,
                        selector,
                        sequence: parsed[0].0.clone(),
                        candidate_models: candidate_models.clone(),
                        neutral_losses: neutral_losses.clone(),
                        engine_neutral_losses: engine_losses,
                        ion_types: ion_types.clone(),
                        max_charge,
                        tolerance_ppm: a.tolerance_ppm,
                        min_intensity: a.min_intensity,
                        min_correlation: a.min_correlation,
                    };
                    entry.insert(
                        tdeasenexus_toolbox::PreparedCandidateComparison::new(&query)
                            .map_err(|error| mass_spec_error(Some(&a.file), error))?,
                    )
                }
            };
            let compared = if let Some(observation) = &observation {
                prepared.match_observation(observation, a.activation.as_deref())
            } else {
                prepared.match_deconvoluted(deconvoluted.as_ref().unwrap())
            }
            .map_err(|error| mass_spec_error(Some(&a.file), error))?;
            let flagged: Vec<String> = compared
                .candidates
                .iter()
                .filter_map(|candidate| {
                    candidate
                        .series_divergence
                        .as_ref()
                        .filter(|d| d.flagged)
                        .map(|d| {
                            format!(
                                "{} {:+.2} ppm (threshold {:.2})",
                                candidate.candidate_id, d.n_minus_c_ppm, d.threshold_ppm
                            )
                        })
                })
                .collect();
            if !flagged.is_empty() {
                series_warnings.push(format!(
                        "series error divergence (N-minus-C signed mean over threshold — model-side mass/site misassignment or isotope borrowing): {}",
                        flagged.join("; ")
                    ));
            }
            if !compared.isotope_borrow_pairs.is_empty() {
                let mut agg: std::collections::BTreeMap<
                    (String, String),
                    (usize, Option<&IsotopeBorrowPair>),
                > = Default::default();
                for pair in &compared.isotope_borrow_pairs {
                    let entry = agg
                        .entry((pair.borrower.clone(), pair.lender.clone()))
                        .or_default();
                    entry.0 += 1;
                    if entry.1.is_none() {
                        entry.1 = Some(pair);
                    }
                }
                let parts: Vec<String> = agg
                        .iter()
                        .map(|((borrower, lender), (count, example))| {
                            let ex = example.expect("count > 0 implies an example");
                            let ion = serde_json::to_value(ex.ion).unwrap_or_default();
                            format!(
                                "{borrower} borrows {count} isotope peak(s) from {lender} (e.g. {}{} q{}: {:+.2} neutrons, {:+.1} ppm residual)",
                                ion.as_str().unwrap_or("?"),
                                ex.number,
                                ex.charge,
                                ex.neutron_offset,
                                ex.residual_ppm
                            )
                        })
                        .collect();
                series_warnings.push(format!(
                        "isotope-borrowing signature: same ion matched at integer-neutron-spaced observed peaks — the borrower's match counts include the other candidate's fragment isotope peaks, not independent support: {}",
                        parts.join("; ")
                    ));
            }

            if !compared.ion_comparison_available {
                series_warnings.push("positional ion comparison unavailable: different backbones or neutral-loss variants; independent model measurements remain available".into());
            }
            let compared = serde_json::to_value(compared).unwrap();
            let compact = compact_fragment_comparison(&compared);
            let mut entry = metadata;
            entry["measurement"] = compared;
            evidence_comparisons.push(entry);
            stdout_comparisons.push(
                json!({"scan": scan, "activation": activation, "ionTypes": ion_types,
                "maxCharge": max_charge, "maxChargeSource": max_charge_source,
                "measurement": compact}),
            );
        }
    }
    params["preparedTheorySets"] = json!(prepared_matches.len() + prepared_comparisons.len());
    let full = single.unwrap_or_else(|| json!({"measurementRole": "candidate_model_evaluation", "notIdentity": true,
        "hypothesisType": hypothesis_type, "scans": scans,
        "aggregation": {"sumScans": a.sum_scans, "binPpm": a.bin_ppm, "minimumAbundancePct": a.minimum_abundance_pct},
        "comparisons": evidence_comparisons}));
    // The artifact records the invocation that produced the measurement —
    // not where the artifact itself was filed. Dropping evidenceOut keeps a
    // measurement written to two destinations byte-identical.
    let mut filed_params = params.clone();
    if let Some(object) = filed_params.as_object_mut() {
        object.remove("evidenceOut");
    }
    let mut artifact = OutputEnvelope::new("fragment", filed_params, sources.clone(), full);
    artifact.warnings = series_warnings.clone();
    let evidence_ref = write_evidence_artifact(evidence_path, &artifact)?;
    let result = if single_match {
        json!({"measurementRole": "candidate_model_evaluation", "notIdentity": true, "hypothesisType": hypothesis_type,
            "measurement": compact_fragment_match(&artifact.result["measurement"]), "evidenceRef": evidence_ref})
    } else {
        json!({"measurementRole": "candidate_model_evaluation", "notIdentity": true, "hypothesisType": hypothesis_type,
        "scans": scans, "comparisons": stdout_comparisons, "evidenceRef": evidence_ref})
    };
    let mut output = OutputEnvelope::new("fragment", params, sources, result);
    output.warnings = series_warnings;
    Ok(output)
}

fn write_evidence_artifact(
    path: &str,
    payload: &impl serde::Serialize,
) -> Result<Value, CommandError> {
    let content = serde_json::to_vec_pretty(payload)
        .map_err(|error| CommandError::InvalidQuery(error.to_string()))?;
    if let Some(parent) = Path::new(path)
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent).map_err(|error| {
            CommandError::InvalidQuery(format!(
                "cannot create evidence directory {parent:?}: {error}"
            ))
        })?;
    }
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_err(|error| CommandError::InvalidQuery(error.to_string()))?
        .as_nanos();
    let temporary = Path::new(path).with_extension(format!("tmp-{}-{nonce}", std::process::id()));
    let written = (|| -> std::io::Result<()> {
        use std::io::Write;
        let mut file = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)?;
        file.write_all(&content)?;
        file.sync_all()?;
        drop(file);
        fs::rename(&temporary, path)
    })();
    if let Err(error) = written {
        let _ = fs::remove_file(&temporary);
        return Err(CommandError::InvalidQuery(format!(
            "cannot write evidence artifact {path}: {error}"
        )));
    }
    Ok(json!({
        "path": path,
        "sha256": full_digest(&content),
        "size": content.len(),
        "mediaType": "application/json"
    }))
}

fn compact_object(value: &Value, excluded: &[&str]) -> serde_json::Map<String, Value> {
    value
        .as_object()
        .into_iter()
        .flat_map(|object| object.iter())
        .filter(|(name, _)| !excluded.contains(&name.as_str()))
        .map(|(name, value)| (name.clone(), value.clone()))
        .collect()
}

fn table_count(value: &Value, name: &str) -> usize {
    value
        .get(name)
        .and_then(Value::as_array)
        .map_or(0, Vec::len)
}

fn compact_fragment_comparison(value: &Value) -> Value {
    let mut summary = compact_object(
        value,
        &["candidates", "isotopeBorrowPairs", "discriminatingIons"],
    );
    let discriminating_count = table_count(value, "discriminatingIons");
    summary.insert("discriminatingIonCount".into(), json!(discriminating_count));
    let preview: Vec<_> = value["discriminatingIons"]
        .as_array()
        .into_iter()
        .flatten()
        .take(32)
        .cloned()
        .collect();
    summary.insert("discriminatingIons".into(), json!(preview));
    summary.insert(
        "discriminatingIonsTruncated".into(),
        json!(discriminating_count > 32),
    );
    let mut by_borrower = std::collections::BTreeMap::<&str, usize>::new();
    for pair in value["isotopeBorrowPairs"].as_array().into_iter().flatten() {
        if let Some(borrower) = pair["borrower"].as_str() {
            *by_borrower.entry(borrower).or_default() += 1;
        }
    }
    summary.insert(
        "isotopeBorrowPairCount".into(),
        json!(table_count(value, "isotopeBorrowPairs")),
    );
    summary.insert("isotopeBorrowBy".into(), json!(by_borrower));
    let candidates: Vec<_> = value["candidates"]
        .as_array()
        .into_iter()
        .flatten()
        .map(|candidate| {
            let mut row = compact_object(candidate, &["matches", "fragmentStatuses"]);
            row.insert(
                "matchTableCount".into(),
                json!(table_count(candidate, "matches")),
            );
            row.insert(
                "fragmentStatusCount".into(),
                json!(table_count(candidate, "fragmentStatuses")),
            );
            Value::Object(row)
        })
        .collect();
    summary.insert("candidates".into(), json!(candidates));
    Value::Object(summary)
}

fn compact_fragment_match(value: &Value) -> Value {
    let mut summary = compact_object(value, &["matches", "fragmentStatuses", "topUnmatchedPeaks"]);
    summary.insert(
        "matchTableCount".into(),
        json!(table_count(value, "matches")),
    );
    summary.insert(
        "fragmentStatusCount".into(),
        json!(table_count(value, "fragmentStatuses")),
    );
    summary.insert(
        "unmatchedPreviewCount".into(),
        json!(table_count(value, "topUnmatchedPeaks")),
    );
    Value::Object(summary)
}

/// Fragment-ion label (spectrum_utils-style lowercase: c, y, z·, z').
fn ion_label(ion: IonType) -> &'static str {
    match ion {
        IonType::A => "a",
        IonType::B => "b",
        IonType::C => "c",
        IonType::X => "x",
        IonType::Y => "y",
        IonType::Z => "z",
    }
}

/// The residue sequence carried by a fragment ion, with each site-specific
/// modification annotated inline at exactly the residue it modifies
/// (`SGR[Acetyl]GK`). N-terminal ions (a/b/c) carry residues `1..=cleavage_after`;
/// C-terminal ions (x/y/z) carry `cleavage_after+1..=len`. A modification at
/// position `p` is rendered at residue `p` — one rule, no special cases.
fn fragment_sequence(
    sequence: &str,
    parsed: &ParsedProForma,
    ion: IonType,
    cleavage_after: usize,
) -> String {
    let n_terminal = matches!(ion, IonType::A | IonType::B | IonType::C);
    let len = sequence.len();
    let (lo, hi) = if n_terminal {
        (0, cleavage_after)
    } else {
        (cleavage_after, len)
    };
    // The fragment's residues, with the terminal modification prepended/appended
    // once (nTerm on N-terminal ions, cTerm on C-terminal ions).
    let prefix = if n_terminal {
        terminal_tag(parsed, "nTerm")
    } else {
        String::new()
    };
    let suffix = if n_terminal {
        String::new()
    } else {
        terminal_tag(parsed, "cTerm")
    };
    let mut out = String::new();
    out.push_str(&prefix);
    for (i, ch) in sequence.chars().enumerate().take(hi).skip(lo) {
        out.push(ch);
        // Site modifications render at exactly the residue they modify.
        let pos = i + 1;
        for m in &parsed.modifications {
            let site: Option<usize> = m
                .position
                .as_deref()
                .and_then(|p| p.strip_prefix("residue:"))
                .and_then(|s| s.parse().ok());
            if site == Some(pos) {
                out.push('[');
                out.push_str(&m.tag);
                out.push(']');
            }
        }
    }
    out.push_str(&suffix);
    out
}

/// The `[Tag]` for a terminal modification of the given position ("" if none).
fn terminal_tag(parsed: &ParsedProForma, position: &str) -> String {
    let mut out = String::new();
    for m in &parsed.modifications {
        if m.position.as_deref() == Some(position) {
            out.push('[');
            out.push_str(&m.tag);
            out.push(']');
        }
    }
    out
}

/// Render the parsed composition as a canonical net-formula string
/// (`C35H63N15O12`, `C9H8O2Br`), sorted by element symbol.
fn net_formula_string(counts: &[tdeasenexus_toolbox::ElementCount]) -> String {
    let mut parts: Vec<String> = counts
        .iter()
        .filter(|c| c.count != 0)
        .map(|c| {
            let symbol = if c.isotope != 0 {
                format!("{}{}", c.isotope, c.element)
            } else {
                c.element.clone()
            };
            if c.count == 1 {
                symbol
            } else {
                format!("{symbol}{}", c.count)
            }
        })
        .collect();
    parts.sort();
    parts.join("")
}

/// `fragment-theoretical`: generate the theoretical fragment ions for a
/// ProFormA with no observed spectrum — the spectrum_utils
/// `get_theoretical_fragments` analogue. JSON by default (full fragment
/// detail + per-fragment sequence annotation); `--out-format mztab` emits a
/// standard mzTab-M molecular table, and `--out-format mzpaf` emits standard
/// HUPO-PSI mzPAF fragment annotations (one per line) for external consumers.
fn fragment_theoretical(a: &FragmentTheoreticalArgs) -> Result<DispatchOutput, CommandError> {
    let out_format = a.out_format.as_deref().unwrap_or("json");
    if !matches!(out_format, "json" | "mztab" | "mzpaf")
        || (a.models.is_some() && out_format != "json")
    {
        return Err(CommandError::InvalidQuery(
            "fragment-theoretical accepts json/mztab/mzpaf; batch models use json".into(),
        ));
    }
    if a.max_charge == 0 {
        return Err(CommandError::InvalidQuery(
            "maxCharge must be positive".into(),
        ));
    }
    let ion_types = parse_ion_types(&a.ion_types)?;
    let (neutral_losses, engine_losses) = a
        .neutral_losses
        .as_deref()
        .map(parse_neutral_losses)
        .transpose()?
        .unwrap_or_default();
    let models: Vec<CandidateProforma> = if let Some(value) = &a.models {
        parse_json_arg("--models", value)?
    } else {
        vec![CandidateProforma {
            id: "model".into(),
            label: None,
            proforma: a.proforma.clone().ok_or_else(|| {
                CommandError::MissingArgument(
                    "fragment-theoretical requires --proforma or --models".into(),
                )
            })?,
        }]
    };
    if models.is_empty() {
        return Err(CommandError::InvalidQuery(
            "nonempty theoretical models required".into(),
        ));
    }
    let (mods, source) = load_mods(a.mods.as_deref())?;
    let mut ids = std::collections::HashSet::new();
    let prepared: Vec<_> = models
        .iter()
        .map(|model| {
            if model.id.trim().is_empty() || !ids.insert(&model.id) {
                return Err(CommandError::InvalidQuery(
                    "model ids must be nonempty and unique".into(),
                ));
            }
            let parsed = parse_proforma(&model.proforma, &mods).map_err(|error| {
                CommandError::ProformaParse(format!("model {}: {error}", model.id))
            })?;
            let normalized =
                tdeasenexus_toolbox::normalize_for_fragmentation(&model.proforma, &mods).map_err(
                    |error| CommandError::ProformaParse(format!("model {}: {error}", model.id)),
                )?;
            Ok((parsed, normalized))
        })
        .collect::<Result<_, CommandError>>()?;
    let generate = |proforma: &str,
                    parsed: &ParsedProForma,
                    fragmentation: &tdeasenexus_toolbox::FragmentationProForma|
     -> Result<DispatchOutput, CommandError> {
        let hypothesis_type = hypothesis_rung(parsed);
        // The apex anchor + full isotope envelope are computed only for the JSON
        // detail output — the interchange formats only need the fragment count /
        // mzPAF lines and skip the BRAIN pass entirely. `--no-envelope` drops the
        // envelope stick (but keeps apexMz) for bulk batch output.
        let out_format = a.out_format.as_deref().unwrap_or("json");
        let envelope_mode = if out_format != "json" || a.no_envelope {
            EnvelopeMode::None
        } else if hypothesis_type == "massDelta" {
            EnvelopeMode::Averagine
        } else {
            EnvelopeMode::Exact
        };
        let fragments = theoretical_fragment_list_proforma_variants(
            &fragmentation.variants,
            &fragmentation.ambiguous_ranges,
            &ion_types,
            a.max_charge,
            &neutral_losses,
            engine_losses,
            envelope_mode,
        )
        .map_err(|error| CommandError::InvalidQuery(error.to_string()))?;

        match out_format {
            "json" => {
                let fragment_values: Vec<Value> = fragments
                .iter()
                .map(|f| {
                    let envelope: Vec<Value> = f
                        .isotope_envelope
                        .iter()
                        .map(|peak| {
                            json!({
                                "neutronOffset": peak.neutron_offset,
                                "mz": peak.mz,
                                "relativeIntensity": peak.relative_intensity,
                            })
                        })
                        .collect();
                    json!({
                        "ion": ion_label(f.ion),
                        "number": f.number,
                        "charge": f.charge,
                        "cleavageAfter": f.cleavage_after,
                        "neutralMass": f.neutral_mass,
                        "mz": f.mz,
                        // Most-abundant isotopologue m/z + its neutron index
                        // above M0 (the fragment apex anchor for top-down
                        // matching; Predator MCP 2024 rationale).
                        "apexMz": f.apex_mz,
                        "apexOffset": f.apex_offset,
                        // Full BRAIN isotope stick: per-isotopologue m/z and
                        // apex-normalized relative intensity. Omitted (empty
                        // array) under `--no-envelope` or in the interchange
                        // formats.
                        "isotopeEnvelope": envelope,
                        "sequence": fragment_sequence(&parsed.sequence, parsed, f.ion, f.cleavage_after),
                        "neutralLoss": f.neutral_loss,
                    })
                })
                .collect();
                let mut warnings = parsed.warnings.clone();
                if envelope_mode == EnvelopeMode::Averagine {
                    warnings.push(
                        "isotope envelopes use the averagine model of each fragment's mass \
                     (mass-only hypothesis): neutron positions are exact, intensities are \
                     an approximation and never a shape claim"
                            .to_owned(),
                    );
                }
                let payload = json!({
                    "proforma": proforma,
                    "hypothesisType": hypothesis_type,
                    "sequence": parsed.sequence,
                    "monoisotopicMass": parsed.monoisotopic_mass,
                    "chemicalFormula": net_formula_string(&parsed.composition),
                    "modifications": parsed.modifications,
                    "warnings": warnings,
                    "ionTypes": a.ion_types,
                    "maxCharge": a.max_charge,
                    "fragmentCount": fragments.len(),
                    "fragments": fragment_values,
                });
                Ok(DispatchOutput::Json(OutputEnvelope::new(
                    "fragment-theoretical",
                    json!({ "proforma": proforma, "mods": a.mods, "outFormat": out_format }),
                    source.clone(),
                    payload,
                )))
            }
            "mztab" => {
                let mods_label = parsed
                    .modifications
                    .iter()
                    .map(|m| {
                        let site = m.position.as_deref().unwrap_or("?");
                        format!("{}({})", m.tag, site)
                    })
                    .collect::<Vec<_>>()
                    .join(",");
                let mut mz = MzTabWriter::new();
                mz.meta("mzTab-version", "2.0.0-M");
                mz.meta("mzTab-mode", "Complete");
                mz.meta("mzTab-type", "Identification");
                mz.meta("mzTab-ID", "TDEase-fragment-theoretical");
                mz.meta(
                    "description",
                    format!("theoretical fragments for {}", proforma),
                );
                mz.set_molecular_columns(vec![
                    "mzTab_id",
                    "SMILES",
                    "chemical_formula",
                    "unique",
                    "opt_mzTab_proforma",
                    "opt_mzTab_sequence",
                    "opt_mzTab_modifications",
                    "opt_mzTab_theoretical_mass",
                    "opt_mzTab_fragment_count",
                    "opt_mzTab_ion_types",
                ]);
                mz.add_molecular_row(vec![
                    "1".to_string(),
                    String::new(), // SMILES — no stereo chemistry for a protein form
                    net_formula_string(&parsed.composition),
                    "0".to_string(),
                    proforma.to_owned(),
                    parsed.sequence.clone(),
                    mods_label,
                    format!("{:.6}", parsed.monoisotopic_mass),
                    fragments.len().to_string(),
                    a.ion_types.clone(),
                ]);
                Ok(DispatchOutput::Tab(mz.render()))
            }
            "mzpaf" => {
                let mzpaf = theoretical_fragments_mzpaf_proforma_variants(
                    &fragmentation.variants,
                    &fragmentation.ambiguous_ranges,
                    &ion_types,
                    a.max_charge,
                    &neutral_losses,
                )
                .map_err(|error| CommandError::InvalidQuery(error.to_string()))?;
                Ok(DispatchOutput::Tab(mzpaf))
            }
            other => Err(CommandError::InvalidQuery(format!(
                "unknown --out-format {other:?}; expected json, mztab or mzpaf"
            ))),
        }
    };
    if a.models.is_none() {
        return generate(&models[0].proforma, &prepared[0].0, &prepared[0].1);
    }
    let rows: Vec<_> = models
        .iter()
        .zip(&prepared)
        .map(
            |(model, (parsed, normalized))| match generate(&model.proforma, parsed, normalized)? {
                DispatchOutput::Json(output) => {
                    Ok(json!({"id": model.id, "label": model.label, "measurement": output.result}))
                }
                DispatchOutput::Tab(_) | DispatchOutput::Mgf(_) => Err(CommandError::InvalidQuery(
                    "batch fragment theory requires json".into(),
                )),
            },
        )
        .collect::<Result<_, CommandError>>()?;
    Ok(DispatchOutput::Json(OutputEnvelope::new(
        "fragment-theoretical",
        json!({"models": models, "ionTypes": a.ion_types, "maxCharge": a.max_charge,
            "mods": a.mods, "noEnvelope": a.no_envelope, "neutralLosses": a.neutral_losses}),
        source,
        json!({"models": rows, "count": rows.len()}),
    )))
}

/// `--out-format mgf`: render the selected scans as a Mascot Generic Format
/// spectrum list. MGF is the community-standard raw-m/z peak exchange format:
/// each scan is one `BEGIN IONS` block with `PEPMASS`/`CHARGE`/`RTINSECONDS`
/// and the (m/z, intensity) pairs. quantify.py (isoforma) consumes this
/// directly — the raw-m/z domain, unlike MSAlign's neutral-mass domain, keeps
/// the discriminating fragments.
fn scans_to_mgf(
    file: &str,
    result: &tdeasenexus_toolbox::ScansQueryResult,
    charge_override: Option<i32>,
) -> Result<DispatchOutput, CommandError> {
    if result.spectra.is_empty() {
        return Err(CommandError::InvalidQuery(
            "mgf export requires at least one matching scan".to_owned(),
        ));
    }
    let mut out = String::new();
    out.push_str(&format!(
        "MGF from tdeasenexus-toolbox spectrum (scans) on {file}; {} spectrum(es)\n",
        result.returned
    ));
    for batch in &result.spectra {
        let spectrum = &batch.spectrum;
        out.push_str("BEGIN IONS\n");
        // Activation rides in TITLE (freeform MGF practice) so MGF consumers can
        // separate complementary series — ETD and HCD scans of one precursor
        // must not be blindly summed. The `TITLE=scan=` prefix is stable;
        // isoforma's parser reads PEPMASS/CHARGE and ignores TITLE.
        let mut title = format!("TITLE=scan={}", batch.scan);
        if let Some(activation) = spectrum
            .precursor
            .as_ref()
            .and_then(|p| p.activation.as_deref())
        {
            title.push_str(&format!(" activation={activation}"));
        }
        out.push_str(&title);
        out.push('\n');
        if let Some(precursor) = &spectrum.precursor {
            // PEPMASS keeps the raw acquisition fact; CHARGE prefers an
            // explicit analyst override — acquisition headers are known to
            // carry wrong charges (survey-frame charges on chlamy H4), and
            // downstream MGF filters compute neutrals from this header.
            out.push_str(&format!("PEPMASS={:.6}\n", precursor.selected_mz));
            let charge = charge_override.or(precursor.charge);
            if let Some(charge) = charge {
                let sign = if charge < 0 { "-" } else { "+" };
                out.push_str(&format!("CHARGE={}{}\n", charge.abs(), sign));
            }
        }
        if spectrum.retention_time_unit == "minute" {
            out.push_str(&format!(
                "RTINSECONDS={:.3}\n",
                spectrum.retention_time * 60.0
            ));
        }
        for peak in &spectrum.peaks {
            out.push_str(&format!("{:.6} {:.4}\n", peak.mz, peak.intensity));
        }
        out.push_str("END IONS\n");
    }
    Ok(DispatchOutput::Mgf(out))
}

// ── UniMod reference queries (D12-R) ────────────────────────────────────────

/// The read-only source provenance attached to every unimod answer.
fn db_source(meta: &tdeasenexus_toolbox::UniModSourceMeta) -> Value {
    json!({
        "database": "unimod",
        "digest": meta.source_sha256,
        "version": format!("schema-{}", meta.schema_version),
        "release": format!("{}.{}", meta.release_major, meta.release_minor),
    })
}

/// Provenance for the database file itself (path, size, head digest).
fn db_source_info(path: &str) -> Result<SourceInfo, CommandError> {
    let metadata = fs::metadata(path).map_err(|error| CommandError::FileRead {
        path: path.to_owned(),
        error: error.to_string(),
    })?;
    let mut bytes = Vec::new();
    fs::File::open(path)
        .map_err(|error| CommandError::FileRead {
            path: path.to_owned(),
            error: error.to_string(),
        })?
        .take(1 << 20)
        .read_to_end(&mut bytes)
        .map_err(|error| CommandError::FileRead {
            path: path.to_owned(),
            error: error.to_string(),
        })?;
    Ok(SourceInfo {
        path: path.to_owned(),
        size: metadata.len(),
        digest: SourceDigest {
            algo: "sha256:head1MiB".to_owned(),
            value: head_digest(&bytes),
        },
    })
}

fn split_comma(text: &Option<String>) -> Vec<String> {
    text.as_deref()
        .map(|value| {
            value
                .split(',')
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .map(str::to_owned)
                .collect()
        })
        .unwrap_or_default()
}

fn approved_filter(value: &str) -> Option<ApprovedFilter> {
    match value {
        "only" => Some(ApprovedFilter::Only),
        "prefer" => Some(ApprovedFilter::Prefer),
        _ => None,
    }
}

fn tol(unit: ToleranceUnit, value: f64) -> Tol {
    Tol {
        unit: match unit {
            ToleranceUnit::Da => TolUnit::Da,
            ToleranceUnit::Ppm => TolUnit::Ppm,
        },
        value,
    }
}

/// Load the official + user stores and merge records, names and specificities.
/// A missing user store contributes nothing (dual-source is additive, D17).
#[allow(clippy::type_complexity)]
fn load_merged_sources(
    db_path: &str,
    user_db: Option<&str>,
) -> Result<
    (
        Vec<tdeasenexus_toolbox::UniModRecord>,
        Vec<tdeasenexus_toolbox::UniModNameEntry>,
        Vec<tdeasenexus_toolbox::UniModSpecificity>,
    ),
    CommandError,
> {
    let db = crate::unimod_db::UnimodDb::open(db_path)?;
    let mut records = db.all_records()?;
    let mut names = db.names()?;
    let mut specificities = db.all_specificities()?;
    if let Some(user) = user_db {
        records.extend(crate::unimod_db::load_user_records(user)?);
        names.extend(crate::unimod_db::load_user_names(user)?);
        specificities.extend(crate::unimod_db::load_user_specificities(user)?);
    }
    Ok((records, names, specificities))
}

fn unimod_get(a: &UnimodGetArgs) -> Result<OutputEnvelope, CommandError> {
    let db_path = crate::unimod_db::resolve_db_path(a.db.as_deref())?;
    let (records, _names, specificities) = load_merged_sources(&db_path, a.user_db.as_deref())?;
    let db = crate::unimod_db::UnimodDb::open(&db_path)?;
    let meta = db.meta()?;
    let source = db_source(&meta);
    let record = if let Some(id) = a.id {
        // An explicit --id targets the official `UNIMOD:n` space only.
        records.iter().find(|r| r.record_id == id).cloned()
    } else if let Some(title) = a.title.as_deref() {
        let matches: Vec<&tdeasenexus_toolbox::UniModRecord> = records
            .iter()
            .filter(|r| r.title.to_lowercase() == title.trim().to_lowercase())
            .collect();
        match matches.len() {
            0 => None,
            1 => Some(matches[0].clone()),
            _ => {
                return Err(CommandError::InvalidQuery(format!(
                    "PROFORMA_MODIFICATION_AMBIGUOUS: {title} matches both UNIMOD and a user record; use an explicit [UNIMOD:<id>] or [USER:<id>]"
                )));
            }
        }
    } else {
        return Err(CommandError::MissingArgument(
            "unimod-get requires --id or --title".to_owned(),
        ));
    };
    let result = match record {
        Some(record) => {
            let rec_specificities: Vec<tdeasenexus_toolbox::UniModSpecificity> = specificities
                .iter()
                .filter(|s| s.modification_id == record.record_id)
                .cloned()
                .collect();
            let mut value = serde_json::to_value(&record).unwrap();
            if let Value::Object(map) = &mut value {
                map.insert(
                    "unimodId".to_owned(),
                    json!(tdeasenexus_toolbox::display_id(
                        record.record_id,
                        &record.source
                    )),
                );
                map.insert("custom".to_owned(), json!(record.is_custom()));
                map.insert(
                    "specificities".to_owned(),
                    serde_json::to_value(rec_specificities).unwrap(),
                );
            }
            json!({ "source": source, "record": value })
        }
        None => json!({ "source": source, "record": null }),
    };
    Ok(OutputEnvelope::new(
        "unimod-get",
        json!({ "db": a.db, "userDb": a.user_db, "id": a.id, "title": a.title }),
        vec![db_source_info(&db_path)?],
        result,
    ))
}

fn unimod_search(a: &UnimodSearchArgs) -> Result<OutputEnvelope, CommandError> {
    let db_path = crate::unimod_db::resolve_db_path(a.db.as_deref())?;
    let (records, names, specificities) = load_merged_sources(&db_path, a.user_db.as_deref())?;
    let db = crate::unimod_db::UnimodDb::open(&db_path)?;
    let meta = db.meta()?;
    let query = SearchQuery {
        text: a.text.clone(),
        mass: a.mass,
        tolerance: a.tolerance.map(|value| tol(a.tolerance_unit, value)),
        filters: UniModFilters {
            residues: split_comma(&a.residues),
            approved: approved_filter(&a.approved),
            contains: split_comma(&a.contains),
            excludes: split_comma(&a.excludes),
            keep_classification: split_comma(&a.keep_classification),
        },
        limit: a.limit,
    };
    let matches = search_by_name(&records, &names, &specificities, &query);
    let mass_matches = search_by_mass(&records, &specificities, &query);
    let result = json!({
        "source": db_source(&meta),
        "matches": matches,
        "massMatches": mass_matches,
    });
    Ok(OutputEnvelope::new(
        "unimod-search",
        json!({
            "db": a.db, "userDb": a.user_db, "text": a.text, "mass": a.mass,
            "tolerance": a.tolerance, "toleranceUnit": a.tolerance_unit.as_str(),
            "contains": a.contains, "excludes": a.excludes, "residues": a.residues,
            "approved": a.approved, "keepClassification": a.keep_classification, "limit": a.limit,
        }),
        vec![db_source_info(&db_path)?],
        result,
    ))
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct UnimodBatchQuery {
    id: String,
    mass: f64,
    tolerance: f64,
    /// Per-row fixed modifications, same grammar as `--fixed-mod`.
    #[serde(default, rename = "fixedMods")]
    fixed_mods: Vec<String>,
}

fn unimod_batch_search(a: &UnimodBatchSearchArgs) -> Result<OutputEnvelope, CommandError> {
    let (queries, mut input_sources): (Vec<UnimodBatchQuery>, _) =
        load_json_arg("--input", &a.input)?;
    validate_batch_ids(queries.iter().map(|query| query.id.as_str()))?;
    if queries
        .iter()
        .any(|query| query.id.is_empty() || query.tolerance < 0.0)
    {
        return Err(CommandError::InvalidQuery(
            "unimod-batch-search requires non-empty ids and non-negative tolerances".to_owned(),
        ));
    }

    let db_path = crate::unimod_db::resolve_db_path(a.db.as_deref())?;
    let (records, _names, specificities) = load_merged_sources(&db_path, a.user_db.as_deref())?;
    let db = crate::unimod_db::UnimodDb::open(&db_path)?;
    let meta = db.meta()?;
    let filters = UniModFilters {
        residues: split_comma(&a.residues),
        approved: approved_filter(&a.approved),
        contains: split_comma(&a.contains),
        excludes: split_comma(&a.excludes),
        keep_classification: split_comma(&a.keep_classification),
    };
    let results: Vec<Value> = queries
        .iter()
        .map(|item| {
            let query = SearchQuery {
                text: None,
                mass: Some(item.mass),
                tolerance: Some(tol(a.tolerance_unit, item.tolerance)),
                filters: filters.clone(),
                limit: a.limit,
            };
            json!({
                "id": item.id,
                "mass": item.mass,
                "tolerance": item.tolerance,
                "massMatches": search_by_mass(&records, &specificities, &query),
            })
        })
        .collect();
    Ok(OutputEnvelope::new(
        "unimod-batch-search",
        json!({
            "input": a.input,
            "db": a.db,
            "userDb": a.user_db,
            "toleranceUnit": a.tolerance_unit.as_str(),
            "contains": a.contains,
            "excludes": a.excludes,
            "residues": a.residues,
            "approved": a.approved,
            "keepClassification": a.keep_classification,
            "limit": a.limit,
            "count": results.len(),
        }),
        {
            input_sources.push(db_source_info(&db_path)?);
            input_sources
        },
        json!({ "source": db_source(&meta), "queries": results }),
    ))
}

fn unimod_batch_solve(a: &UnimodBatchSolveArgs) -> Result<OutputEnvelope, CommandError> {
    let (queries, mut input_sources): (Vec<UnimodBatchQuery>, _) =
        load_json_arg("--input", &a.input)?;
    validate_batch_ids(queries.iter().map(|query| query.id.as_str()))?;
    if queries
        .iter()
        .any(|query| query.id.is_empty() || query.tolerance < 0.0)
    {
        return Err(CommandError::InvalidQuery(
            "unimod-batch-solve requires non-empty ids and non-negative tolerances".to_owned(),
        ));
    }

    let db_path = crate::unimod_db::resolve_db_path(a.db.as_deref())?;
    let (records, names, specificities) = load_merged_sources(&db_path, a.user_db.as_deref())?;
    let db = crate::unimod_db::UnimodDb::open(&db_path)?;
    let meta = db.meta()?;
    let neutral_losses = db.all_neutral_losses()?;
    let (records_unfiltered, specificities_unfiltered) = (records.clone(), specificities.clone());
    // Same default artefact exclusions as unimod-solve — both commands solve
    // unknown composition gaps over the same pool (see unimod_solve).
    let default_exclusions = "Isotopic label,AA substitution,Artefact";
    let exclude_classification = a
        .exclude_classification
        .clone()
        .or_else(|| Some(default_exclusions.to_owned()));
    let (records, specificities) = filter_solve_records(
        records,
        specificities,
        &a.keep_classification,
        &exclude_classification,
        &a.element_whitelist,
    );
    let filters_applied = a
        .keep_classification
        .as_deref()
        .is_some_and(|v| !v.is_empty())
        || exclude_classification
            .as_deref()
            .is_some_and(|v| !v.is_empty())
        || !a.element_whitelist.is_empty()
        || !a.allowed_sites.is_empty();
    if filters_applied {
        eprintln!(
            "WARNING unimod-batch-solve: classification/site filters applied — an unconstrainedBaseline (same query, no filters) is attached to every solution; chemistry present only there was removed by the filters."
        );
    }
    let allowed_sites = parse_solve_sites(&a.allowed_sites)?;
    let row_fixed_entries: Vec<Vec<String>> = queries
        .iter()
        .map(|item| {
            a.fixed_mods
                .iter()
                .chain(&item.fixed_mods)
                .cloned()
                .collect()
        })
        .collect();
    let row_fixed_mods: Vec<Vec<FixedMod>> = row_fixed_entries
        .iter()
        .map(|entries| resolve_fixed_mods(entries, &records_unfiltered, &records, &names))
        .collect::<Result<_, _>>()?;
    let results: Vec<Value> = queries
        .iter()
        .zip(row_fixed_mods)
        .zip(&row_fixed_entries)
        .map(|((item, fixed_mods), fixed_entries)| {
            let query = SolveQuery {
                target_mass: item.mass,
                tolerance: tol(a.tolerance_unit, item.tolerance),
                max_modifications: a.max_modifications,
                max_distinct_modifications: a.max_distinct,
                allow_repeated_modification: a.allow_repeated,
                approved: approved_filter(&a.approved),
                allowed_sites: allowed_sites.clone(),
                fixed_modifications: fixed_mods,
                limit: a.limit,
            };
            let solution = solve_with_unconstrained_baseline(
                &records,
                &specificities,
                &records_unfiltered,
                &specificities_unfiltered,
                &neutral_losses,
                &query,
                filters_applied,
            )?;
            Ok(json!({
                "id": item.id,
                "mass": item.mass,
                "tolerance": item.tolerance,
                "fixedMods": fixed_entries,
                "solution": solution,
            }))
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(OutputEnvelope::new(
        "unimod-batch-solve",
        json!({
            "input": a.input,
            "db": a.db,
            "userDb": a.user_db,
            "toleranceUnit": a.tolerance_unit.as_str(),
            "maxModifications": a.max_modifications,
            "maxDistinct": a.max_distinct,
            "allowRepeated": a.allow_repeated,
            "approved": a.approved,
            "keepClassification": a.keep_classification,
            "excludeClassification": exclude_classification,
            "elementWhitelist": a.element_whitelist,
            "allowedSites": a.allowed_sites,
            "fixedMods": a.fixed_mods,
            "limit": a.limit,
            "count": results.len(),
        }),
        {
            input_sources.push(db_source_info(&db_path)?);
            input_sources
        },
        json!({ "source": db_source(&meta), "queries": results }),
    ))
}
/// Extract element names from a UniMod composition string
/// ("H(-7) O Fe(3)", "Se", "C(6) H(12) N(2) O(2)"). Isotope prefixes
/// (15N) still yield the base element (N).
fn parse_composition_elements(composition: &str) -> std::collections::HashSet<String> {
    let mut out = std::collections::HashSet::new();
    let mut chars = composition.chars().peekable();
    while let Some(ch) = chars.next() {
        if ch.is_ascii_uppercase() {
            let mut el = ch.to_string();
            if let Some(&next) = chars.peek() {
                if next.is_ascii_lowercase() {
                    el.push(next);
                    chars.next();
                }
            }
            out.insert(el);
        }
    }
    out
}

fn filter_solve_records(
    mut records: Vec<tdeasenexus_toolbox::UniModRecord>,
    mut specificities: Vec<tdeasenexus_toolbox::UniModSpecificity>,
    keep_classification: &Option<String>,
    exclude_classification: &Option<String>,
    element_whitelist: &str,
) -> (
    Vec<tdeasenexus_toolbox::UniModRecord>,
    Vec<tdeasenexus_toolbox::UniModSpecificity>,
) {
    let keep = split_comma(keep_classification);
    if !keep.is_empty() {
        let kept_ids: std::collections::HashSet<u32> = specificities
            .iter()
            .filter(|item| keep.iter().any(|class| class == &item.classification))
            .map(|item| item.modification_id)
            .collect();
        records.retain(|item| kept_ids.contains(&item.record_id));
        specificities.retain(|item| kept_ids.contains(&item.modification_id));
    }

    let excluded: std::collections::HashSet<String> =
        split_comma(exclude_classification).into_iter().collect();
    if !excluded.is_empty() {
        let class_by_mod: std::collections::HashMap<u32, Vec<&str>> =
            specificities
                .iter()
                .fold(std::collections::HashMap::new(), |mut map, item| {
                    map.entry(item.modification_id)
                        .or_default()
                        .push(&item.classification);
                    map
                });
        records.retain(|item| {
            class_by_mod.get(&item.record_id).is_none_or(|classes| {
                classes.is_empty() || !classes.iter().all(|class| excluded.contains(*class))
            })
        });
    }

    let whitelist: std::collections::HashSet<String> =
        split_comma(&Some(element_whitelist.to_owned()))
            .into_iter()
            .collect();
    if !whitelist.is_empty() {
        records.retain(|item| parse_composition_elements(&item.composition).is_subset(&whitelist));
    }
    let retained: std::collections::HashSet<u32> =
        records.iter().map(|item| item.record_id).collect();
    specificities.retain(|item| retained.contains(&item.modification_id));
    (records, specificities)
}

/// Resolve one `--fixed-mod` reference: `UNIMOD:n` / `USER:n` / bare record id,
/// or an exact (case-insensitive) title from the merged name table.
fn resolve_fixed_mod_reference(reference: &str, names: &[UniModNameEntry]) -> Option<u32> {
    if let Some(rest) = reference.strip_prefix("UNIMOD:") {
        return rest.parse::<u32>().ok();
    }
    if let Some(rest) = reference.strip_prefix("USER:") {
        return rest.parse::<u32>().ok().map(|n| n + USER_OFFSET);
    }
    if let Ok(n) = reference.parse::<u32>() {
        return Some(n);
    }
    let needle = reference.trim().to_lowercase();
    names
        .iter()
        .find(|entry| entry.name.to_lowercase() == needle)
        .map(|entry| entry.modification_id)
}

/// Parse one `--fixed-mod TITLE|ID[:COUNT]` entry into a pinned (record, count).
fn parse_fixed_mod_entry(
    entry: &str,
    names: &[UniModNameEntry],
) -> Result<(u32, usize), CommandError> {
    let invalid = |message: String| {
        CommandError::InvalidQuery(format!("invalid --fixed-mod \"{entry}\": {message}"))
    };
    if let Some(record_id) = resolve_fixed_mod_reference(entry, names) {
        return Ok((record_id, 1));
    }
    let (prefix, suffix) = entry
        .rsplit_once(':')
        .ok_or_else(|| invalid("unknown modification title or id".to_owned()))?;
    let count = suffix
        .parse::<usize>()
        .ok()
        .filter(|count| *count >= 1)
        .ok_or_else(|| invalid("count must be a positive integer".to_owned()))?;
    let record_id = resolve_fixed_mod_reference(prefix, names)
        .ok_or_else(|| invalid(format!("unknown modification title or id \"{prefix}\"")))?;
    Ok((record_id, count))
}

/// Resolve every pinned entry against the merged record table. A pinned record
/// must exist in the reference database AND survive this query's filters —
/// pinning chemistry the filters removed would silently change the answer.
fn resolve_fixed_mods(
    entries: &[String],
    records_unfiltered: &[UniModRecord],
    records_filtered: &[UniModRecord],
    names: &[UniModNameEntry],
) -> Result<Vec<FixedMod>, CommandError> {
    let present =
        |list: &[UniModRecord], record_id: u32| list.iter().any(|r| r.record_id == record_id);
    entries
        .iter()
        .map(|entry| {
            let (record_id, count) = parse_fixed_mod_entry(entry, names)?;
            if !present(records_unfiltered, record_id) {
                return Err(CommandError::InvalidQuery(format!(
                    "--fixed-mod \"{entry}\": record {record_id} does not exist in the reference database"
                )));
            }
            if !present(records_filtered, record_id) {
                return Err(CommandError::InvalidQuery(format!(
                    "--fixed-mod \"{entry}\": record {record_id} was removed by this query's classification/element filters and cannot be pinned"
                )));
            }
            Ok(FixedMod { record_id, count })
        })
        .collect()
}

/// Filtered solve plus the unconstrained baseline when query-side filters
/// (classification/element/site) were applied. Shared by `unimod-solve` and
/// `unimod-batch-solve` so both modes return the same shape — the cli
/// contract requires batch to preserve single-query constraints, baseline
/// included. Chemistry present only in the baseline was removed by the
/// filters and must enter the comparison manifest.
#[allow(clippy::too_many_arguments)]
fn solve_with_unconstrained_baseline(
    records: &[tdeasenexus_toolbox::UniModRecord],
    specificities: &[tdeasenexus_toolbox::UniModSpecificity],
    unfiltered: &[tdeasenexus_toolbox::UniModRecord],
    unfiltered_specificities: &[tdeasenexus_toolbox::UniModSpecificity],
    neutral_losses: &[tdeasenexus_toolbox::UniModNeutralLoss],
    query: &SolveQuery,
    filters_applied: bool,
) -> Result<Value, CommandError> {
    let mut result = serde_json::to_value(
        solve_combinations(records, specificities, neutral_losses, query)
            .map_err(|error| mass_spec_error(None, error))?,
    )
    .unwrap();
    if !filters_applied {
        return Ok(result);
    }
    let baseline_query = SolveQuery {
        allowed_sites: Vec::new(),
        ..query.clone()
    };
    let baseline = solve_combinations(
        unfiltered,
        unfiltered_specificities,
        neutral_losses,
        &baseline_query,
    )
    .map_err(|error| mass_spec_error(None, error))?;
    let mut baseline_object = match serde_json::to_value(&baseline).unwrap() {
        Value::Object(map) => map,
        _ => serde_json::Map::new(),
    };
    let baseline_candidates = match baseline_object.remove("candidates") {
        Some(Value::Array(items)) => Value::Array(items.into_iter().take(15).collect()),
        _ => Value::Array(Vec::new()),
    };
    baseline_object.insert("candidates".to_owned(), baseline_candidates);
    baseline_object.insert(
        "note".to_owned(),
        json!("unconstrained re-run of the same query (no classification/site filters); top 15 shown — any chemistry present here but absent from the filtered candidates was removed by this query's filters and must enter the comparison manifest"),
    );
    if let Value::Object(map) = &mut result {
        map.insert(
            "unconstrainedBaseline".to_owned(),
            Value::Object(baseline_object),
        );
    }
    Ok(result)
}

fn unimod_solve(a: &UnimodSolveArgs) -> Result<OutputEnvelope, CommandError> {
    let db_path = crate::unimod_db::resolve_db_path(a.db.as_deref())?;
    let (records, names, specificities) = load_merged_sources(&db_path, a.user_db.as_deref())?;
    let db = crate::unimod_db::UnimodDb::open(&db_path)?;
    let meta = db.meta()?;
    // Default artefact exclusions for an unknown-gap solver: isotopic labels
    // (the label is pinned when a labelling study is being solved, never
    // searched for), amino-acid substitutions (the sequence-variant
    // hypothesis belongs to its own lane, not to composition search) and the
    // explicitly artefactual class. Without these, label/substitution records
    // dominate every error-ordered window and bury biological answers.
    // `--exclude-classification ""` takes manual control and disables the
    // default; the unconstrained baseline keeps the excluded chemistry
    // inspectable either way.
    let default_exclusions = "Isotopic label,AA substitution,Artefact";
    let exclude_classification = a
        .exclude_classification
        .clone()
        .or_else(|| Some(default_exclusions.to_owned()));
    let filters_applied = a
        .keep_classification
        .as_deref()
        .is_some_and(|v| !v.is_empty())
        || exclude_classification
            .as_deref()
            .is_some_and(|v| !v.is_empty())
        || !a.element_whitelist.is_empty()
        || !a.allowed_sites.is_empty();
    let (records_unfiltered, specificities_unfiltered) = (records.clone(), specificities.clone());
    let (records, specificities) = filter_solve_records(
        records,
        specificities,
        &a.keep_classification,
        &exclude_classification,
        &a.element_whitelist,
    );
    let fixed_mods = resolve_fixed_mods(&a.fixed_mods, &records_unfiltered, &records, &names)?;
    let neutral_losses = db.all_neutral_losses()?;
    let query = SolveQuery {
        target_mass: a.target_mass,
        tolerance: tol(a.tolerance_unit, a.tolerance),
        max_modifications: a.max_modifications,
        max_distinct_modifications: a.max_distinct,
        allow_repeated_modification: a.allow_repeated,
        approved: approved_filter(&a.approved),
        allowed_sites: parse_solve_sites(&a.allowed_sites)?,
        fixed_modifications: fixed_mods,
        limit: a.limit,
    };
    let mut result = solve_with_unconstrained_baseline(
        &records,
        &specificities,
        &records_unfiltered,
        &specificities_unfiltered,
        &neutral_losses,
        &query,
        filters_applied,
    )?;
    if let Value::Object(map) = &mut result {
        map.insert("source".to_owned(), db_source(&meta));
    }
    if filters_applied {
        eprintln!(
            "WARNING unimod-solve: classification/site filters applied — an unconstrainedBaseline (same query, no filters) is attached; chemistry present only there was removed by the filters."
        );
    }
    if result
        .get("truncated")
        .and_then(Value::as_bool)
        .unwrap_or(false)
    {
        eprintln!(
            "WARNING unimod-solve: result TRUNCATED — candidates beyond the returned window were dropped. The returned set is NOT exhaustive: narrow --tolerance, add --allowed-site constraints, or raise --limit before ranking or claiming no viable candidate."
        );
    }
    Ok(OutputEnvelope::new(
        "unimod-solve",
        json!({
            "db": a.db, "userDb": a.user_db, "targetMass": a.target_mass,
            "tolerance": a.tolerance, "toleranceUnit": a.tolerance_unit.as_str(),
            "maxModifications": a.max_modifications, "maxDistinct": a.max_distinct,
            "allowRepeated": a.allow_repeated, "approved": a.approved,
            "keepClassification": a.keep_classification,
            "excludeClassification": exclude_classification,
            "elementWhitelist": a.element_whitelist,
            "allowedSites": a.allowed_sites, "fixedMods": a.fixed_mods, "limit": a.limit,
        }),
        vec![db_source_info(&db_path)?],
        result,
    ))
}

fn parse_solve_sites(values: &[String]) -> Result<Vec<SolveSite>, CommandError> {
    values
        .iter()
        .map(|value| {
            let (residue, position) = value.split_once(':').unwrap_or((value.as_str(), ""));
            let residue = residue.trim().to_ascii_uppercase();
            if residue.len() != 1 || !residue.as_bytes()[0].is_ascii_alphabetic() {
                return Err(CommandError::InvalidQuery(format!(
                    "invalid allowed site {value}"
                )));
            }
            let position = if position.trim().is_empty() {
                None
            } else {
                Some(position.trim().parse::<usize>().map_err(|_| {
                    CommandError::InvalidQuery(format!("invalid allowed-site position in {value}"))
                })?)
            };
            Ok(SolveSite { residue, position })
        })
        .collect()
}

fn unimod_residues(a: &UnimodResiduesArgs) -> Result<OutputEnvelope, CommandError> {
    let db_path = crate::unimod_db::resolve_db_path(a.db.as_deref())?;
    let db = crate::unimod_db::UnimodDb::open(&db_path)?;
    let meta = db.meta()?;
    let amino_acids = db.amino_acids()?;
    let result = json!({ "source": db_source(&meta), "aminoAcids": amino_acids });
    Ok(OutputEnvelope::new(
        "unimod-residues",
        json!({ "db": a.db, "userDb": a.user_db }),
        vec![db_source_info(&db_path)?],
        result,
    ))
}

fn unimod_add(a: &UnimodAddArgs) -> Result<OutputEnvelope, CommandError> {
    if a.title.trim().is_empty() {
        return Err(CommandError::MissingArgument(
            "unimod-add requires --title".to_owned(),
        ));
    }
    if a.formula.trim().is_empty() {
        return Err(CommandError::MissingArgument(
            "unimod-add requires --formula".to_owned(),
        ));
    }
    let sites = split_comma(&a.sites);
    let user_id = crate::unimod_db::add_user_modification(
        &a.user_db,
        a.title.trim(),
        a.title.trim(),
        a.formula.trim(),
        a.mono_mass,
        a.average_mass.unwrap_or(a.mono_mass),
        &sites,
    )?;
    let result = json!({
        "added": true,
        "userId": user_id,
        "unimodId": format!("USER:{}", user_id),
    });
    Ok(OutputEnvelope::new(
        "unimod-add",
        json!({
            "userDb": a.user_db, "title": a.title, "formula": a.formula,
            "monoMass": a.mono_mass, "averageMass": a.average_mass, "sites": a.sites,
        }),
        Vec::new(),
        result,
    ))
}
