//! UniMod candidate search and combination solving — the reference-data query
//! kernel (D12-R). This is pure logic over in-memory records; SQLite lives in
//! the CLI (`unimod_db.rs`), so the same functions back both the CLI commands
//! and the divergence tests. TS and Python reach UniMod through the CLI.
//!
//! Ported from the former TS layer (`unimod-search.ts`, `unimod-combination.ts`)
//! so the ranking tiers, composition filter and combination ordering are a
//! single Rust kernel instead of a TS implementation that duplicated the CLI.

use std::collections::HashMap;

use serde::Serialize;

use crate::error::MassSpecError;
use crate::proforma::split_net_formula;

/// Offset applied to user modification ids (`USER:n`) so they never collide
/// with official `UNIMOD:n` when both sources are merged (D17). Official ids
/// are ~1..1560; the offset keeps the spaces disjoint.
pub const USER_OFFSET: u32 = 1_000_000;

/// The display id for a record: `UNIMOD:<n>` for the official source,
/// `USER:<n>` (id minus the offset) for a user record.
pub fn display_id(record_id: u32, source: &str) -> String {
    if source == "user" {
        format!("USER:{}", record_id.saturating_sub(USER_OFFSET))
    } else {
        format!("UNIMOD:{record_id}")
    }
}

// ── Reference data shapes (mirror the DB projection) ────────────────────────

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct UniModRecord {
    pub record_id: u32,
    pub title: String,
    pub full_name: String,
    pub approved: bool,
    pub delta_mono: f64,
    pub delta_average: f64,
    pub composition: String,
    /// The modification source: `"unimod"` (official read-only) or `"user"`
    /// (writable `user-mods.sqlite`). `custom` derives from this (D17).
    pub source: String,
}

impl UniModRecord {
    pub fn is_custom(&self) -> bool {
        self.source == "user"
    }
}

#[derive(Debug, Clone)]
pub struct UniModNameEntry {
    pub modification_id: u32,
    pub name: String,
    pub kind: NameKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NameKind {
    Title,
    FullName,
    AltName,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct UniModSpecificity {
    pub id: u32,
    pub modification_id: u32,
    pub site: String,
    pub position: String,
    pub classification: String,
    pub spec_group_id: Option<u32>,
    pub hidden: bool,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct UniModNeutralLoss {
    pub specificity_id: u32,
    pub mono_mass: f64,
    pub average_mass: f64,
    pub composition: String,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct UniModAminoAcid {
    pub title: String,
    pub three_letter: String,
    pub full_name: String,
    pub mono_mass: f64,
    pub average_mass: f64,
}

/// The read-only source provenance carried with every answer.
#[derive(Debug, Clone, Serialize)]
pub struct UniModSourceMeta {
    pub source_sha256: String,
    pub schema_version: u32,
    pub release_major: String,
    pub release_minor: String,
}

// ── Queries ─────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct Tol {
    pub unit: TolUnit,
    pub value: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TolUnit {
    Da,
    Ppm,
}

#[derive(Debug, Clone, Default)]
pub struct UniModFilters {
    pub residues: Vec<String>,
    pub approved: Option<ApprovedFilter>,
    pub contains: Vec<String>,
    pub excludes: Vec<String>,
    /// Record-level classification gate. When non-empty, a record is kept only
    /// if at least one of its specificities has a `classification` in this set.
    /// This is the *correct* PTM filter: UniMod's `approved` flag is unreliable
    /// (Oxidation — the most common PTM — is `approved=false`). Keep the
    /// real-PTM classes (Post-translational, Artefact, Multiple, glycosylation)
    /// and drop AA substitutions / isotopic labels / pure chemical derivatives.
    /// The check is **record-level, not per-site**: a mod like Oxidation, whose
    /// Met site is "Artefact" but other sites are "Chemical derivative", is kept
    /// — filtering per-site would invert the biology (drop the canonical site).
    pub keep_classification: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApprovedFilter {
    Only,
    Prefer,
}

#[derive(Debug, Clone, Default)]
pub struct SearchQuery {
    pub text: Option<String>,
    pub mass: Option<f64>,
    pub tolerance: Option<Tol>,
    pub filters: UniModFilters,
    pub limit: usize,
}

#[derive(Debug, Clone)]
pub struct SolveSite {
    pub residue: String,
    pub position: Option<usize>,
}

/// `count` copies of a modification pinned into every candidate combination
/// before the variable search starts — the "fixed modification" of
/// database-search engines (e.g. Carbamidomethyl on C). The solver enumerates
/// only the residual mass `target − Σ(delta×count)` with the remaining
/// modification budget; pinned instances never count against
/// `allow_repeated_modification` (they are not searched, only reported).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FixedMod {
    pub record_id: u32,
    pub count: usize,
}

#[derive(Debug, Clone)]
pub struct SolveQuery {
    pub target_mass: f64,
    pub tolerance: Tol,
    pub max_modifications: usize,
    pub max_distinct_modifications: Option<usize>,
    pub allow_repeated_modification: bool,
    pub approved: Option<ApprovedFilter>,
    pub allowed_sites: Vec<SolveSite>,
    pub fixed_modifications: Vec<FixedMod>,
    pub limit: usize,
}

// ── Search output ───────────────────────────────────────────────────────────

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct NameMatch {
    pub record_id: u32,
    pub title: String,
    pub full_name: String,
    pub approved: bool,
    pub delta_mono: f64,
    pub composition: String,
    pub source: String,
    pub custom: bool,
    /// 1 (strongest) .. 8 (weakest).
    pub text_rank: u32,
    pub matched_field: &'static str,
    pub matched_text: String,
    pub match_reason: String,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct MassMatch {
    pub record_id: u32,
    pub title: String,
    pub approved: bool,
    pub delta_mono: f64,
    pub theoretical_delta_mass: f64,
    pub signed_error_da: f64,
    pub absolute_error_da: f64,
    pub signed_error_ppm: f64,
    pub absolute_error_ppm: f64,
    pub feasible: bool,
    pub source: String,
    pub custom: bool,
    pub specificities: Vec<MassSpecificity>,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct MassSpecificity {
    pub site: String,
    pub position: String,
    pub classification: String,
}

/// PPM is defined as (theoretical − observed) / |observed| × 10⁶. Near-zero
/// mass differences make ppm unstable, so a Da tolerance is always honored on
/// its own merits even when the relative error is large (mirrors the TS rule).
#[derive(Debug, Clone, Copy)]
pub struct MassError {
    pub signed_da: f64,
    pub abs_da: f64,
    pub signed_ppm: f64,
    pub abs_ppm: f64,
    pub within: bool,
}

pub fn mass_error(observed: f64, theoretical: f64, tol: &Tol) -> MassError {
    let signed_da = theoretical - observed;
    let abs_da = signed_da.abs();
    let denom = observed.abs();
    let signed_ppm = if denom > 0.0 {
        (signed_da / denom) * 1e6
    } else {
        0.0
    };
    let abs_ppm = signed_ppm.abs();
    let within = match tol.unit {
        TolUnit::Da => abs_da <= tol.value,
        TolUnit::Ppm => denom > 0.0 && abs_ppm <= tol.value,
    };
    MassError {
        signed_da,
        abs_da,
        signed_ppm,
        abs_ppm,
        within,
    }
}

/// Element-membership filter on a UniMod net composition. The grammar is the
/// shared §5.1 superset (`split_net_formula`); two filter semantics layer on
/// top, matching the TS `parseComposition`: stable-isotope prefixes (`2H`,
/// `13C`) fold into their base element, and a lowercase-led token (`dHex`,
/// `HexNAc`) is a carbohydrate shorthand with no elements.
pub fn composition_matches(composition: &str, contains: &[String], excludes: &[String]) -> bool {
    if contains.is_empty() && excludes.is_empty() {
        return true;
    }
    let counts = composition_counts(composition);
    contains
        .iter()
        .all(|element| counts.get(element).copied().unwrap_or(0) >= 1)
        && excludes
            .iter()
            .all(|element| counts.get(element).copied().unwrap_or(0) == 0)
}

fn composition_counts(composition: &str) -> HashMap<String, i32> {
    let mut counts = HashMap::new();
    // `dHex`, `HexNAc` … start lowercase and carry no element counts.
    if composition
        .chars()
        .next()
        .is_some_and(|c| c.is_ascii_lowercase())
    {
        return counts;
    }
    for (sym, count) in split_net_formula(composition).unwrap_or_default() {
        let base = sym.trim_start_matches(|c: char| c.is_ascii_digit());
        *counts.entry(base.to_string()).or_insert(0) += count;
    }
    counts
}

// ── Name search: explainable 8-tier ranking ─────────────────────────────────

pub fn search_by_name(
    records: &[UniModRecord],
    names: &[UniModNameEntry],
    specificities: &[UniModSpecificity],
    query: &SearchQuery,
) -> Vec<NameMatch> {
    let Some(text) = query.text.as_deref() else {
        return Vec::new();
    };
    let needle = text.trim().to_lowercase();
    if needle.is_empty() {
        return Vec::new();
    }
    let record_by_id: HashMap<u32, &UniModRecord> = records
        .iter()
        .map(|record| (record.record_id, record))
        .collect();
    let specs_by_mod = specificities_by_mod(specificities);

    let mut matches: Vec<NameMatch> = Vec::new();
    let mut seen = std::collections::HashSet::new();
    for entry in names {
        let Some(record) = record_by_id.get(&entry.modification_id) else {
            continue;
        };
        if !passes_filters(record, query.filters.approved) {
            continue;
        }
        let specs = specs_by_mod
            .get(&record.record_id)
            .cloned()
            .unwrap_or_default();
        if !classification_ok(&specs, &query.filters.keep_classification) {
            continue;
        }
        let Some(rank) = rank_name(&needle, entry) else {
            continue;
        };
        if seen.contains(&record.record_id) {
            // Keep the strongest reason for a record that has several names.
            if let Some(existing) = matches.iter_mut().find(|m| m.record_id == record.record_id) {
                if rank.tier < existing.text_rank {
                    existing.text_rank = rank.tier;
                    existing.matched_field = rank.field;
                    existing.matched_text = rank.matched_text;
                    existing.match_reason = rank.reason;
                }
            }
            continue;
        }
        seen.insert(record.record_id);
        matches.push(NameMatch {
            record_id: record.record_id,
            title: record.title.clone(),
            full_name: record.full_name.clone(),
            approved: record.approved,
            delta_mono: record.delta_mono,
            composition: record.composition.clone(),
            source: record.source.clone(),
            custom: record.is_custom(),
            text_rank: rank.tier,
            matched_field: rank.field,
            matched_text: rank.matched_text,
            match_reason: rank.reason,
        });
    }

    // Sort by tier, then title length (shorter first), then title for stability.
    matches.sort_by(|a, b| {
        a.text_rank
            .cmp(&b.text_rank)
            .then(a.title.len().cmp(&b.title.len()))
            .then(a.title.cmp(&b.title))
    });
    matches.truncate(query.limit.clamp(1, 1000));
    matches
}

struct NameRank {
    tier: u32,
    field: &'static str,
    matched_text: String,
    reason: String,
}

fn rank_name(needle: &str, entry: &UniModNameEntry) -> Option<NameRank> {
    let hay = entry.name.to_lowercase();
    let is_title = entry.kind == NameKind::Title;
    let is_alt = entry.kind == NameKind::FullName || entry.kind == NameKind::AltName;
    // 1 — exact accession (numeric id as text); 2 — exact title. The TS code
    // checked both under the same `is_title` branch, so an exact title ranks 1
    // with field "accession"; preserved verbatim.
    if is_title && hay == *needle {
        return Some(NameRank {
            tier: 1,
            field: "accession",
            matched_text: entry.name.clone(),
            reason: format!("exact accession \"{}\"", entry.name),
        });
    }
    if is_title && hay == *needle {
        return Some(NameRank {
            tier: 2,
            field: "title",
            matched_text: entry.name.clone(),
            reason: format!("exact title \"{}\"", entry.name),
        });
    }
    // 3 — exact full/alt name.
    if is_alt && hay == *needle {
        return Some(NameRank {
            tier: 3,
            field: if entry.kind == NameKind::FullName {
                "full_name"
            } else {
                "alt_name"
            },
            matched_text: entry.name.clone(),
            reason: format!(
                "exact {} \"{}\"",
                if entry.kind == NameKind::FullName {
                    "full_name"
                } else {
                    "alt_name"
                },
                entry.name
            ),
        });
    }
    // 4 — title prefix.
    if is_title && hay.starts_with(needle) {
        return Some(NameRank {
            tier: 4,
            field: "title",
            matched_text: entry.name.clone(),
            reason: format!("title prefix \"{}\"", entry.name),
        });
    }
    // 5 — token prefix (split on non-alphanumeric).
    let tokens: Vec<&str> = hay
        .split(|c: char| !c.is_ascii_alphanumeric())
        .filter(|t| !t.is_empty())
        .collect();
    if tokens.iter().any(|t| t.starts_with(needle)) {
        return Some(NameRank {
            tier: 5,
            field: if entry.kind == NameKind::Title {
                "title"
            } else if entry.kind == NameKind::FullName {
                "full_name"
            } else {
                "alt_name"
            },
            matched_text: entry.name.clone(),
            reason: format!("token prefix in \"{}\"", entry.name),
        });
    }
    // 6 — substring.
    if hay.contains(needle) {
        return Some(NameRank {
            tier: 6,
            field: if entry.kind == NameKind::Title {
                "title"
            } else if entry.kind == NameKind::FullName {
                "full_name"
            } else {
                "alt_name"
            },
            matched_text: entry.name.clone(),
            reason: format!("substring of \"{}\"", entry.name),
        });
    }
    // 7 — trigram similarity.
    if trigram_similarity(needle, &hay) >= 0.4 {
        return Some(NameRank {
            tier: 7,
            field: if entry.kind == NameKind::Title {
                "title"
            } else if entry.kind == NameKind::FullName {
                "full_name"
            } else {
                "alt_name"
            },
            matched_text: entry.name.clone(),
            reason: format!("trigram similarity to \"{}\"", entry.name),
        });
    }
    // 8 — edit distance (cheap; only short needles).
    if needle.len() <= 32 && levenshtein(needle, &hay) <= std::cmp::max(1, needle.len() / 3) {
        return Some(NameRank {
            tier: 8,
            field: if entry.kind == NameKind::Title {
                "title"
            } else if entry.kind == NameKind::FullName {
                "full_name"
            } else {
                "alt_name"
            },
            matched_text: entry.name.clone(),
            reason: format!("edit distance to \"{}\"", entry.name),
        });
    }
    None
}

fn passes_filters(record: &UniModRecord, approved: Option<ApprovedFilter>) -> bool {
    !matches!(approved, Some(ApprovedFilter::Only) if !record.approved)
}

fn trigram_similarity(a: &str, b: &str) -> f64 {
    if a.is_empty() || b.is_empty() {
        return if a == b { 1.0 } else { 0.0 };
    }
    let ta = trigrams(a);
    let tb = trigrams(b);
    if ta.is_empty() || tb.is_empty() {
        return if a == b { 1.0 } else { 0.0 };
    }
    let shared = ta.iter().filter(|g| tb.contains(*g)).count();
    shared as f64 / std::cmp::max(ta.len(), tb.len()) as f64
}

fn trigrams(s: &str) -> std::collections::HashSet<String> {
    let padded = format!(" {s} ");
    let chars: Vec<char> = padded.chars().collect();
    let mut set = std::collections::HashSet::new();
    for window in chars.windows(3) {
        set.insert(window.iter().collect::<String>());
    }
    set
}

fn levenshtein(a: &str, b: &str) -> usize {
    let a: Vec<char> = a.chars().collect();
    let b: Vec<char> = b.chars().collect();
    let (m, n) = (a.len(), b.len());
    if m == 0 {
        return n;
    }
    if n == 0 {
        return m;
    }
    let mut prev: Vec<usize> = (0..=n).collect();
    let mut curr = vec![0usize; n + 1];
    for i in 1..=m {
        curr[0] = i;
        for j in 1..=n {
            let cost = if a[i - 1] == b[j - 1] { 0 } else { 1 };
            curr[j] = (curr[j - 1] + 1).min(prev[j] + 1).min(prev[j - 1] + cost);
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    prev[n]
}

// ── Mass search ─────────────────────────────────────────────────────────────

pub fn search_by_mass(
    records: &[UniModRecord],
    specificities: &[UniModSpecificity],
    query: &SearchQuery,
) -> Vec<MassMatch> {
    let Some(mass) = query.mass else {
        return Vec::new();
    };
    let Some(tol) = query.tolerance.as_ref() else {
        return Vec::new();
    };
    let specs_by_mod = specificities_by_mod(specificities);
    let mut matches: Vec<MassMatch> = Vec::new();
    for record in records {
        if !passes_filters(record, query.filters.approved) {
            continue;
        }
        if !composition_matches(
            &record.composition,
            &query.filters.contains,
            &query.filters.excludes,
        ) {
            continue;
        }
        let specs = specs_by_mod
            .get(&record.record_id)
            .cloned()
            .unwrap_or_default();
        if !classification_ok(&specs, &query.filters.keep_classification) {
            continue;
        }
        let theoretical = record.delta_mono;
        let error = mass_error(mass, theoretical, tol);
        if !error.within {
            continue;
        }
        let feasible = site_feasible(&specs, &query.filters.residues);
        matches.push(MassMatch {
            record_id: record.record_id,
            title: record.title.clone(),
            approved: record.approved,
            delta_mono: record.delta_mono,
            theoretical_delta_mass: theoretical,
            signed_error_da: error.signed_da,
            absolute_error_da: error.abs_da,
            signed_error_ppm: error.signed_ppm,
            absolute_error_ppm: error.abs_ppm,
            feasible,
            source: record.source.clone(),
            custom: record.is_custom(),
            specificities: specs
                .iter()
                .map(|s| MassSpecificity {
                    site: s.site.clone(),
                    position: s.position.clone(),
                    classification: s.classification.clone(),
                })
                .collect(),
        });
    }
    // Physical fit first; UniMod's editorial approved flag is returned as data
    // but must not steer blind candidate ranking.
    matches.sort_by(|a, b| {
        a.absolute_error_ppm
            .partial_cmp(&b.absolute_error_ppm)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(b.feasible.cmp(&a.feasible))
            .then(a.record_id.cmp(&b.record_id))
    });
    matches.truncate(query.limit.clamp(1, 1000));
    matches
}

fn site_feasible(specs: &[UniModSpecificity], residues: &[String]) -> bool {
    if residues.is_empty() {
        return true;
    }
    let wanted: std::collections::HashSet<String> =
        residues.iter().map(|r| r.to_uppercase()).collect();
    specs
        .iter()
        .any(|s| wanted.contains(&s.site.to_uppercase()))
}

/// Record-level classification gate: keep iff some specificity's classification
/// is in `keep` (empty `keep` = no gate). See `UniModFilters.keep_classification`.
fn classification_ok(specs: &[UniModSpecificity], keep: &[String]) -> bool {
    if keep.is_empty() {
        return true;
    }
    specs
        .iter()
        .any(|s| keep.iter().any(|c| &s.classification == c))
}

fn specificities_by_mod(
    specificities: &[UniModSpecificity],
) -> HashMap<u32, Vec<UniModSpecificity>> {
    let mut map: HashMap<u32, Vec<UniModSpecificity>> = HashMap::new();
    for spec in specificities {
        map.entry(spec.modification_id)
            .or_default()
            .push(spec.clone());
    }
    map
}

// ── Combination solving ─────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct CombinationModification {
    pub unimod_id: String,
    pub record_id: u32,
    pub title: String,
    pub count: usize,
    /// At least one instance of this modification is pinned (`fixed`
    /// modification) — mandatory in every candidate, not found by the search.
    pub fixed: bool,
    pub source: String,
    pub custom: bool,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct CandidateNeutralLoss {
    pub unimod_id: String,
    pub mono_mass: f64,
    pub composition: String,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct CombinationCandidate {
    pub rank: usize,
    pub modifications: Vec<CombinationModification>,
    pub theoretical_delta_mass: f64,
    pub signed_error_da: f64,
    pub absolute_error_da: f64,
    pub signed_error_ppm: f64,
    pub absolute_error_ppm: f64,
    pub modification_count: usize,
    pub distinct_modification_count: usize,
    pub approved_count: usize,
    pub site_feasible: bool,
    pub neutral_losses: Vec<CandidateNeutralLoss>,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct SolveMassResult {
    pub candidates: Vec<CombinationCandidate>,
    pub match_count: usize,
    pub returned_count: usize,
    pub truncated: bool,
    pub completeness: &'static str,
    pub scanned_combination_count: usize,
    pub pruned_combination_count: usize,
}

const MAX_SCAN_BUDGET: usize = 2_000_000;

struct PoolEntry {
    record_id: u32,
    title: String,
    approved: bool,
    delta: f64,
    source: String,
}

struct Half {
    indices: Vec<usize>,
    sum: f64,
}

struct Solved {
    multiset: Vec<usize>,
    theoretical: f64,
}

/// Find all combinations of up to `max_modifications` whose summed mass matches
/// `target_mass` within `tolerance`. The rank is absolute ppm error first (the
/// legacy `(terms, abs(delta_da))` sort let "fewer terms" outrank a closer
/// mass — the corrected order is ppm first, fewer total mods, fewer distinct,
/// then canonical member for stability). The returned window is
/// level-stratified: the best candidate of every distinct-modification level
/// present leads the window (fewest components first), and the remaining
/// slots follow the ppm order — deeper combinations have more free parameters
/// and fit any mass more tightly, so a pure ppm window can bury the most
/// parsimonious viable answer beyond `limit`. `rank` stays the pure ppm
/// order. When sites are supplied, combinations
/// that cannot be jointly assigned to them are excluded before ranking.
///
/// Pinned (`fixed`) modifications are mandatory members of every candidate:
/// the search runs on the residual mass `target − Σ(delta×count)` with the
/// remaining budget, and candidates report those members with `fixed: true`.
/// An unknown pinned id is a broken query and errors; a pinned record pushed
/// out of the search pool by `approved`/site filters keeps its metadata from
/// `records` and simply cannot be placed on a filtered-away site.
pub fn solve_combinations(
    records: &[UniModRecord],
    specificities: &[UniModSpecificity],
    neutral_losses: &[UniModNeutralLoss],
    query: &SolveQuery,
) -> Result<SolveMassResult, MassSpecError> {
    let target = query.target_mass;
    let max_mods = query.max_modifications.clamp(1, 4);
    let limit = query.limit.clamp(1, 1000);
    let record_by_id: HashMap<u32, &UniModRecord> =
        records.iter().map(|r| (r.record_id, r)).collect();
    let mut fixed_counts: HashMap<u32, usize> = HashMap::new();
    for fixed in &query.fixed_modifications {
        if fixed.count == 0 {
            return Err(MassSpecError::InvalidQuery(format!(
                "unimod solve: fixed modification {} has count 0",
                fixed.record_id
            )));
        }
        if !record_by_id.contains_key(&fixed.record_id) {
            return Err(MassSpecError::InvalidQuery(format!(
                "unimod solve: fixed modification {} is not in the record pool",
                fixed.record_id
            )));
        }
        *fixed_counts.entry(fixed.record_id).or_insert(0) += fixed.count;
    }
    let fixed_sum: f64 = fixed_counts
        .iter()
        .map(|(&id, &count)| record_by_id[&id].delta_mono * count as f64)
        .sum();
    let fixed_total: usize = fixed_counts.values().sum();
    let fixed_distinct = fixed_counts.len();
    // Pinned mass leaves the target; pinned counts leave the budget and the
    // distinct-modification cap.
    let residual_target = target - fixed_sum;
    let variable_budget = max_mods.saturating_sub(fixed_total);
    let max_distinct = query
        .max_distinct_modifications
        .unwrap_or(max_mods)
        .min(max_mods)
        .saturating_sub(fixed_distinct);
    let allow_repeat = query.allow_repeated_modification;
    let specs_by_mod = specificities_by_mod(specificities);
    let allowed_residues = query
        .allowed_sites
        .iter()
        .map(|site| site.residue.as_str())
        .collect::<std::collections::HashSet<_>>();

    let mut pool: Vec<PoolEntry> = Vec::new();
    for record in records {
        if matches!(query.approved, Some(ApprovedFilter::Only) if !record.approved) {
            continue;
        }
        if !allowed_residues.is_empty()
            && !specs_by_mod
                .get(&record.record_id)
                .into_iter()
                .flatten()
                .any(|spec| allowed_residues.contains(spec.site.as_str()))
        {
            continue;
        }
        pool.push(PoolEntry {
            record_id: record.record_id,
            title: record.title.clone(),
            approved: record.approved,
            delta: record.delta_mono,
            source: record.source.clone(),
        });
    }
    // delta ascending — required by the pruning and two-pointer joins.
    pool.sort_by(|a, b| {
        a.delta
            .partial_cmp(&b.delta)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let abs_tol_da = match query.tolerance.unit {
        TolUnit::Da => query.tolerance.value,
        TolUnit::Ppm => {
            if target.abs() > 0.0 {
                (query.tolerance.value / 1e6) * target.abs()
            } else {
                query.tolerance.value
            }
        }
    };
    let min_delta = pool.first().map(|p| p.delta).unwrap_or(0.0);
    let max_delta = pool.last().map(|p| p.delta).unwrap_or(0.0);

    let mut state = SolveState {
        pool: &pool,
        target: residual_target,
        abs_tol_da,
        min_delta,
        max_delta,
        allow_repeat,
        max_distinct,
        candidates: Vec::new(),
        scanned: 0,
        pruned: 0,
        budget_exhausted: false,
    };

    // The variable part may be empty: a pinned multiset that already matches
    // the full target is the single candidate. Without pinned modifications
    // this never fires — the empty combination stays excluded (legacy).
    if fixed_total > 0 {
        let error = mass_error(
            residual_target,
            0.0,
            &Tol {
                unit: TolUnit::Da,
                value: abs_tol_da,
            },
        );
        if error.within {
            state.candidates.push(Solved {
                multiset: Vec::new(),
                theoretical: 0.0,
            });
        }
    }

    for modification_count in 1..=variable_budget {
        if modification_count <= 2 {
            state.enumerate_depth(0, modification_count, residual_target, &mut Vec::new());
        } else {
            state.meet_in_the_middle(modification_count);
        }
        if state.budget_exhausted {
            break;
        }
    }

    // Deduplicate identical multisets (allowRepeat produces many index orderings
    // that resolve to the same mod multiset).
    let mut canonical: Vec<&Solved> = Vec::new();
    let mut seen_keys = std::collections::HashSet::new();
    for c in &state.candidates {
        let key = c
            .multiset
            .iter()
            .map(|idx| idx.to_string())
            .collect::<Vec<_>>()
            .join(",");
        if seen_keys.insert(key) {
            canonical.push(c);
        }
    }
    let budget_exhausted = state.budget_exhausted;
    let nl_by_spec = neutral_losses_by_spec(neutral_losses);

    // Build ranked candidates. When the caller supplied sites, they are a
    // query constraint: an unassignable multiset is not a candidate and must
    // not consume the bounded return window.
    let mut ranked: Vec<CombinationCandidate> = canonical
        .iter()
        .map(|c| {
            build_candidate(
                &pool,
                &record_by_id,
                &specs_by_mod,
                &nl_by_spec,
                c,
                &query.allowed_sites,
                &fixed_counts,
                fixed_sum,
                target,
                abs_tol_da,
            )
        })
        .collect();
    if !query.allowed_sites.is_empty() {
        ranked.retain(|candidate| candidate.site_feasible);
    }
    let match_count = ranked.len();
    ranked.sort_by(|a, b| {
        a.absolute_error_ppm
            .partial_cmp(&b.absolute_error_ppm)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.modification_count.cmp(&b.modification_count))
            .then(
                a.distinct_modification_count
                    .cmp(&b.distinct_modification_count),
            )
            .then(
                a.modifications
                    .iter()
                    .map(|m| m.record_id)
                    .collect::<Vec<_>>()
                    .cmp(
                        &b.modifications
                            .iter()
                            .map(|m| m.record_id)
                            .collect::<Vec<_>>(),
                    ),
            )
    });
    for (index, candidate) in ranked.iter_mut().enumerate() {
        candidate.rank = index + 1;
    }

    // Level-stratified window: one guaranteed slot per distinct-modification
    // level (that level's ppm-best, fewest components first), then the rest
    // in ppm order. Deep combinations overfit — with more free parameters
    // they fit any target more tightly — so the ppm tail alone can push the
    // parsimonious viable answer past `limit` and out of the reader's view.
    let mut level_first: Vec<usize> = Vec::new();
    let mut seen_levels = std::collections::HashSet::new();
    for (idx, c) in ranked.iter().enumerate() {
        if seen_levels.insert(c.distinct_modification_count) {
            level_first.push(idx);
        }
    }
    level_first.sort_by_key(|&idx| ranked[idx].distinct_modification_count);
    let mut window: Vec<usize> = Vec::with_capacity(limit);
    for idx in &level_first {
        if window.len() < limit {
            window.push(*idx);
        }
    }
    for idx in 0..ranked.len() {
        if window.len() >= limit {
            break;
        }
        if !level_first.contains(&idx) {
            window.push(idx);
        }
    }
    let returned = window
        .into_iter()
        .map(|i| ranked[i].clone())
        .collect::<Vec<_>>();
    Ok(SolveMassResult {
        candidates: returned.clone(),
        match_count,
        returned_count: returned.len(),
        truncated: returned.len() < match_count || budget_exhausted,
        completeness: if budget_exhausted {
            "partial"
        } else {
            "complete"
        },
        scanned_combination_count: state.scanned,
        pruned_combination_count: state.pruned,
    })
}

struct SolveState<'a> {
    pool: &'a [PoolEntry],
    target: f64,
    abs_tol_da: f64,
    min_delta: f64,
    max_delta: f64,
    allow_repeat: bool,
    max_distinct: usize,
    candidates: Vec<Solved>,
    scanned: usize,
    pruned: usize,
    budget_exhausted: bool,
}

impl<'a> SolveState<'a> {
    fn record_combo(&mut self, indices: &[usize]) {
        self.scanned += 1;
        if self.scanned > MAX_SCAN_BUDGET {
            self.budget_exhausted = true;
            return;
        }
        let theoretical: f64 = indices.iter().map(|&idx| self.pool[idx].delta).sum();
        let error = mass_error(
            self.target,
            theoretical,
            &Tol {
                unit: TolUnit::Da,
                value: self.abs_tol_da,
            },
        );
        if !error.within {
            self.pruned += 1;
            return;
        }
        let mut multiset: Vec<usize> = indices.to_vec();
        multiset.sort_unstable();
        self.candidates.push(Solved {
            multiset,
            theoretical,
        });
    }

    fn enumerate_depth(
        &mut self,
        start_idx: usize,
        terms_left: usize,
        remaining: f64,
        chosen: &mut Vec<usize>,
    ) {
        if self.budget_exhausted {
            return;
        }
        if terms_left == 0 {
            if !chosen.is_empty() {
                self.record_combo(chosen);
            }
            return;
        }
        for i in start_idx..self.pool.len() {
            let delta = self.pool[i].delta;
            let remaining_min = delta + (terms_left - 1) as f64 * self.min_delta;
            let remaining_max = delta + (terms_left - 1) as f64 * self.max_delta;
            if remaining_min > remaining + self.abs_tol_da {
                self.pruned += 1;
                break;
            }
            if remaining_max < remaining - self.abs_tol_da {
                self.pruned += 1;
                continue;
            }
            if !self.allow_repeat && chosen.contains(&i) {
                continue;
            }
            chosen.push(i);
            let next_start = if self.allow_repeat { i } else { i + 1 };
            self.enumerate_depth(next_start, terms_left - 1, remaining - delta, chosen);
            chosen.pop();
            if self.budget_exhausted {
                return;
            }
        }
    }

    fn meet_in_the_middle(&mut self, k: usize) {
        let half = k.div_ceil(2);
        let other = k - half;
        let mut lefts: Vec<Half> = Vec::new();
        let mut rights: Vec<Half> = vec![Half {
            indices: Vec::new(),
            sum: 0.0,
        }];
        self.build_halves(half, 0, 0.0, &mut Vec::new(), &mut lefts);
        self.build_halves(other, 0, 0.0, &mut Vec::new(), &mut rights);

        rights.sort_by(|a, b| {
            a.sum
                .partial_cmp(&b.sum)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let right_sums: Vec<f64> = rights.iter().map(|r| r.sum).collect();

        for left in &lefts {
            if self.budget_exhausted {
                return;
            }
            let need = self.target - left.sum;
            let lo = need - self.abs_tol_da;
            let hi = need + self.abs_tol_da;
            let start = lower_bound(&right_sums, lo);
            for j in start..rights.len() {
                if right_sums[j] > hi {
                    break;
                }
                if self.budget_exhausted {
                    return;
                }
                let mut combined = left.indices.clone();
                combined.extend_from_slice(&rights[j].indices);
                if !self.allow_repeat {
                    let distinct: std::collections::HashSet<usize> =
                        combined.iter().copied().collect();
                    if distinct.len() != combined.len() {
                        continue;
                    }
                }
                let distinct = combined
                    .iter()
                    .collect::<std::collections::HashSet<_>>()
                    .len();
                if distinct > self.max_distinct {
                    self.pruned += 1;
                    continue;
                }
                self.record_combo(&combined);
            }
        }
    }

    fn build_halves(
        &mut self,
        terms_left: usize,
        start_idx: usize,
        sum: f64,
        chosen: &mut Vec<usize>,
        emit: &mut Vec<Half>,
    ) {
        if self.budget_exhausted {
            return;
        }
        if terms_left == 0 {
            emit.push(Half {
                indices: chosen.clone(),
                sum,
            });
            return;
        }
        for i in start_idx..self.pool.len() {
            let delta = self.pool[i].delta;
            if terms_left > 1
                && sum + delta + (terms_left - 1) as f64 * self.min_delta
                    > self.target + self.abs_tol_da
            {
                self.pruned += 1;
                break;
            }
            if !self.allow_repeat && chosen.contains(&i) {
                continue;
            }
            chosen.push(i);
            let next_start = if self.allow_repeat { i } else { i + 1 };
            self.build_halves(terms_left - 1, next_start, sum + delta, chosen, emit);
            chosen.pop();
            if self.budget_exhausted {
                return;
            }
        }
    }
}

/// Merge the variable multiset with the pinned counts and build the reported
/// candidate. `solved.theoretical` is the variable sum only — the pinned mass
/// is added here and the error is recomputed against the full target, so
/// every candidate (pinned members or not) reports the same frame.
fn build_candidate(
    pool: &[PoolEntry],
    record_by_id: &HashMap<u32, &UniModRecord>,
    specs_by_mod: &HashMap<u32, Vec<UniModSpecificity>>,
    nl_by_spec: &HashMap<u32, Vec<&UniModNeutralLoss>>,
    solved: &Solved,
    allowed_sites: &[SolveSite],
    fixed_counts: &HashMap<u32, usize>,
    fixed_sum: f64,
    target: f64,
    abs_tol_da: f64,
) -> CombinationCandidate {
    let mut counts: HashMap<u32, usize> = HashMap::new();
    for &idx in &solved.multiset {
        *counts.entry(pool[idx].record_id).or_insert(0) += 1;
    }
    for (&record_id, &count) in fixed_counts {
        *counts.entry(record_id).or_insert(0) += count;
    }
    let mut modifications: Vec<CombinationModification> = counts
        .iter()
        .map(|(&record_id, &count)| {
            let (title, source) = match pool.iter().find(|p| p.record_id == record_id) {
                Some(p) => (p.title.clone(), p.source.clone()),
                // A pinned id can sit outside the search pool (approved=only,
                // site filters) while still being mandatory; the record table
                // is the metadata authority.
                None => {
                    let r = record_by_id
                        .get(&record_id)
                        .expect("pinned and multiset ids exist in the record table");
                    (r.title.clone(), r.source.clone())
                }
            };
            CombinationModification {
                unimod_id: display_id(record_id, &source),
                record_id,
                title,
                count,
                fixed: fixed_counts.contains_key(&record_id),
                custom: source == "user",
                source,
            }
        })
        .collect();
    modifications.sort_by_key(|m| m.record_id);
    let approved_count = counts
        .keys()
        .filter(
            |record_id| match pool.iter().find(|p| p.record_id == **record_id) {
                Some(p) => p.approved,
                None => record_by_id.get(record_id).is_some_and(|r| r.approved),
            },
        )
        .count();
    // Per-modification neutral losses, resolved through specificities (the port
    // filled these in TS; the kernel owns them now). One entry per specificity's
    // neutral loss, keyed by the mod's UNIMOD id — `count` does not multiply.
    let mut neutral_losses: Vec<CandidateNeutralLoss> = Vec::new();
    for &record_id in counts.keys() {
        let unimod_id = format!("UNIMOD:{record_id}");
        for spec in specs_by_mod.get(&record_id).into_iter().flatten() {
            for nl in nl_by_spec.get(&spec.id).into_iter().flatten() {
                neutral_losses.push(CandidateNeutralLoss {
                    unimod_id: unimod_id.clone(),
                    mono_mass: nl.mono_mass,
                    composition: nl.composition.clone(),
                });
            }
        }
    }
    let theoretical = solved.theoretical + fixed_sum;
    let error = mass_error(
        target,
        theoretical,
        &Tol {
            unit: TolUnit::Da,
            value: abs_tol_da,
        },
    );
    CombinationCandidate {
        rank: 0,
        modifications,
        theoretical_delta_mass: theoretical,
        signed_error_da: error.signed_da,
        absolute_error_da: error.abs_da,
        signed_error_ppm: error.signed_ppm,
        absolute_error_ppm: error.abs_ppm,
        modification_count: counts.values().sum(),
        distinct_modification_count: counts.len(),
        approved_count,
        site_feasible: combination_site_feasible(&counts, specs_by_mod, allowed_sites),
        neutral_losses,
    }
}

fn combination_site_feasible(
    counts: &HashMap<u32, usize>,
    specs_by_mod: &HashMap<u32, Vec<UniModSpecificity>>,
    allowed_sites: &[SolveSite],
) -> bool {
    if allowed_sites.is_empty() {
        return true;
    }
    let mut instances = Vec::new();
    for (&record_id, &count) in counts {
        for _ in 0..count {
            instances.push(record_id);
        }
    }
    if instances.len() > allowed_sites.len() {
        return false;
    }
    instances.sort_by_key(|record_id| {
        specs_by_mod
            .get(record_id)
            .map(|specs| {
                specs
                    .iter()
                    .filter(|spec| allowed_sites.iter().any(|site| site.residue == spec.site))
                    .count()
            })
            .unwrap_or(0)
    });

    fn assign(
        index: usize,
        instances: &[u32],
        specs_by_mod: &HashMap<u32, Vec<UniModSpecificity>>,
        allowed_sites: &[SolveSite],
        used: &mut [bool],
    ) -> bool {
        if index == instances.len() {
            return true;
        }
        let Some(specs) = specs_by_mod.get(&instances[index]) else {
            return false;
        };
        for (site_index, site) in allowed_sites.iter().enumerate() {
            if used[site_index] || !specs.iter().any(|spec| spec.site == site.residue) {
                continue;
            }
            used[site_index] = true;
            if assign(index + 1, instances, specs_by_mod, allowed_sites, used) {
                return true;
            }
            used[site_index] = false;
        }
        false
    }

    assign(
        0,
        &instances,
        specs_by_mod,
        allowed_sites,
        &mut vec![false; allowed_sites.len()],
    )
}

fn neutral_losses_by_spec(
    neutral_losses: &[UniModNeutralLoss],
) -> HashMap<u32, Vec<&UniModNeutralLoss>> {
    let mut map: HashMap<u32, Vec<&UniModNeutralLoss>> = HashMap::new();
    for nl in neutral_losses {
        map.entry(nl.specificity_id).or_default().push(nl);
    }
    map
}

fn lower_bound(arr: &[f64], x: f64) -> usize {
    let mut lo = 0;
    let mut hi = arr.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if arr[mid] < x {
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

    fn records() -> Vec<UniModRecord> {
        vec![
            UniModRecord {
                record_id: 1,
                title: "Acetyl".to_owned(),
                full_name: "Acetylation".to_owned(),
                approved: true,
                delta_mono: 42.010565,
                delta_average: 42.0367,
                composition: "C2 H2 O1".to_owned(),
                source: "unimod".to_owned(),
            },
            UniModRecord {
                record_id: 21,
                title: "Phospho".to_owned(),
                full_name: "Phosphorylation".to_owned(),
                approved: true,
                delta_mono: 79.966331,
                delta_average: 79.9799,
                composition: "H O4 P".to_owned(),
                source: "unimod".to_owned(),
            },
            UniModRecord {
                record_id: 340,
                title: "Bromo".to_owned(),
                full_name: "Bromination".to_owned(),
                approved: true,
                delta_mono: 77.9105,
                delta_average: 77.9309,
                composition: "Br H -1".to_owned(),
                source: "unimod".to_owned(),
            },
        ]
    }

    fn names() -> Vec<UniModNameEntry> {
        let mut names = Vec::new();
        for record in records() {
            names.push(UniModNameEntry {
                modification_id: record.record_id,
                name: record.title.clone(),
                kind: NameKind::Title,
            });
            names.push(UniModNameEntry {
                modification_id: record.record_id,
                name: record.full_name.clone(),
                kind: NameKind::FullName,
            });
        }
        names
    }

    fn spec(mod_id: u32, site: &str) -> UniModSpecificity {
        UniModSpecificity {
            id: mod_id * 10,
            modification_id: mod_id,
            site: site.to_owned(),
            position: String::new(),
            classification: "Anywhere".to_owned(),
            spec_group_id: None,
            hidden: false,
        }
    }

    #[test]
    fn name_search_phospho_ranks_as_title_prefix() {
        let query = SearchQuery {
            text: Some("phosph".to_owned()),
            ..Default::default()
        };
        let matches = search_by_name(&records(), &names(), &[], &query);
        let phospho = matches
            .iter()
            .find(|m| m.record_id == 21)
            .expect("Phospho found");
        assert_eq!(phospho.text_rank, 4, "title prefix tier");
        assert_eq!(phospho.matched_field, "title");
        assert!(phospho.match_reason.contains("prefix"));
    }

    #[test]
    fn exact_title_outranks_substring() {
        let query = SearchQuery {
            text: Some("phosphorylation".to_owned()),
            ..Default::default()
        };
        let matches = search_by_name(&records(), &names(), &[], &query);
        assert_eq!(
            matches.first().unwrap().record_id,
            21,
            "exact full-name ranks first"
        );
        assert!(matches.first().unwrap().text_rank <= 3);
    }

    #[test]
    fn mass_search_with_composition_filter_finds_bromo() {
        let query = SearchQuery {
            mass: Some(77.9105),
            tolerance: Some(Tol {
                unit: TolUnit::Da,
                value: 0.01,
            }),
            filters: UniModFilters {
                contains: vec!["Br".to_owned()],
                excludes: vec!["P".to_owned()],
                ..Default::default()
            },
            ..Default::default()
        };
        let matches = search_by_mass(&records(), &[], &query);
        assert_eq!(matches.len(), 1, "only Bromo carries Br and no P");
        assert_eq!(matches[0].record_id, 340);
        assert!(matches[0].absolute_error_da < 0.01);
    }

    #[test]
    fn composition_matches_folds_isotopes_and_treats_carbs_as_empty() {
        assert!(composition_matches(
            "Br H -1",
            &["Br".to_owned()],
            &["P".to_owned()]
        ));
        assert!(!composition_matches("Br H -1", &["P".to_owned()], &[]));
        // Stable-isotope prefix folds into the base element for the filter.
        assert!(composition_matches("2H(8)", &["H".to_owned()], &[]));
        // A carbohydrate shorthand carries no elements.
        assert!(!composition_matches("dHex", &["H".to_owned()], &[]));
        assert!(composition_matches("dHex", &[], &["P".to_owned()]));
    }

    #[test]
    fn solve_two_phospho_ranks_by_ppm_first() {
        let query = SolveQuery {
            target_mass: 159.932662, // 2 × Phospho (79.966331) exactly
            tolerance: Tol {
                unit: TolUnit::Ppm,
                value: 20.0,
            },
            max_modifications: 3,
            max_distinct_modifications: None,
            allow_repeated_modification: true,
            approved: None,
            allowed_sites: Vec::new(),
            fixed_modifications: Vec::new(),
            limit: 10,
        };
        let specificities = vec![spec(21, "S"), spec(21, "T")];
        let result =
            solve_combinations(&records(), &specificities, &[], &query).expect("solve succeeds");
        assert!(
            !result.candidates.is_empty(),
            "at least one combination matches"
        );
        let best = &result.candidates[0];
        assert_eq!(best.rank, 1);
        assert!(best.absolute_error_ppm < 20.0);
        assert_eq!(
            best.modifications,
            vec![CombinationModification {
                unimod_id: "UNIMOD:21".to_owned(),
                record_id: 21,
                title: "Phospho".to_owned(),
                count: 2,
                fixed: false,
                source: "unimod".to_owned(),
                custom: false,
            }]
        );
        assert!(result.completeness == "complete");
    }

    #[test]
    fn solve_up_to_two_includes_single_and_double_modifications() {
        let synthetic = vec![
            UniModRecord {
                record_id: 1,
                title: "Single80".to_owned(),
                full_name: "Single80".to_owned(),
                approved: true,
                delta_mono: 80.0,
                delta_average: 80.0,
                composition: "".to_owned(),
                source: "unimod".to_owned(),
            },
            UniModRecord {
                record_id: 2,
                title: "Half40".to_owned(),
                full_name: "Half40".to_owned(),
                approved: true,
                delta_mono: 40.0,
                delta_average: 40.0,
                composition: "".to_owned(),
                source: "unimod".to_owned(),
            },
        ];
        let query = SolveQuery {
            target_mass: 80.0,
            tolerance: Tol {
                unit: TolUnit::Da,
                value: 0.001,
            },
            max_modifications: 2,
            max_distinct_modifications: None,
            allow_repeated_modification: true,
            approved: None,
            allowed_sites: Vec::new(),
            fixed_modifications: Vec::new(),
            limit: 100,
        };
        let result = solve_combinations(&synthetic, &[], &[], &query).expect("solve succeeds");
        let counts = result
            .candidates
            .iter()
            .map(|candidate| candidate.modification_count)
            .collect::<std::collections::HashSet<_>>();
        assert!(counts.contains(&1), "single modifications are enumerated");
        assert!(
            counts.contains(&2),
            "two-modification combinations are enumerated"
        );
    }

    #[test]
    fn solve_sites_exclude_combinations_that_cannot_be_placed() {
        let synthetic = vec![
            UniModRecord {
                record_id: 1,
                title: "Whole".to_owned(),
                full_name: "Whole".to_owned(),
                approved: true,
                delta_mono: 80.0,
                delta_average: 80.0,
                composition: "".to_owned(),
                source: "unimod".to_owned(),
            },
            UniModRecord {
                record_id: 2,
                title: "Half".to_owned(),
                full_name: "Half".to_owned(),
                approved: true,
                delta_mono: 40.0,
                delta_average: 40.0,
                composition: "".to_owned(),
                source: "unimod".to_owned(),
            },
        ];
        let query = SolveQuery {
            target_mass: 80.0,
            tolerance: Tol {
                unit: TolUnit::Da,
                value: 0.001,
            },
            max_modifications: 2,
            max_distinct_modifications: None,
            allow_repeated_modification: true,
            approved: None,
            allowed_sites: vec![SolveSite {
                residue: "K".to_owned(),
                position: Some(10),
            }],
            fixed_modifications: Vec::new(),
            limit: 1,
        };
        let specificities = vec![spec(1, "K"), spec(2, "K")];

        let result =
            solve_combinations(&synthetic, &specificities, &[], &query).expect("solve succeeds");

        assert_eq!(result.match_count, 1);
        assert!(!result.truncated);
        assert_eq!(result.candidates.len(), 1);
        assert_eq!(result.candidates[0].modifications[0].title, "Whole");
        assert!(result.candidates[0].site_feasible);
    }

    fn cam_record() -> UniModRecord {
        UniModRecord {
            record_id: 4,
            title: "Carbamidomethyl".to_owned(),
            full_name: "Carbamidomethylation".to_owned(),
            approved: true,
            delta_mono: 57.021464,
            delta_average: 57.0519,
            composition: "C2 H3 N1 O1".to_owned(),
            source: "unimod".to_owned(),
        }
    }

    #[test]
    fn solve_fixed_mod_pins_and_solves_the_residual() {
        let mut pool = records();
        pool.push(cam_record());
        // Decoy matching the FULL target alone: under pin semantics the
        // variable search only sees the residual, so it must be excluded.
        pool.push(UniModRecord {
            record_id: 900,
            title: "Whole".to_owned(),
            full_name: "Whole".to_owned(),
            approved: true,
            delta_mono: 57.021464 + 79.966331,
            delta_average: 137.0,
            composition: String::new(),
            source: "unimod".to_owned(),
        });
        let query = SolveQuery {
            target_mass: 57.021464 + 79.966331, // CAM + Phospho
            tolerance: Tol {
                unit: TolUnit::Da,
                value: 0.001,
            },
            max_modifications: 2,
            max_distinct_modifications: None,
            allow_repeated_modification: true,
            approved: None,
            allowed_sites: Vec::new(),
            fixed_modifications: vec![FixedMod {
                record_id: 4,
                count: 1,
            }],
            limit: 10,
        };
        let result = solve_combinations(&pool, &[], &[], &query).expect("solve succeeds");
        assert_eq!(
            result.match_count, 1,
            "only CAM + Phospho survives the pin: {:?}",
            result.candidates
        );
        let best = &result.candidates[0];
        assert_eq!(best.modification_count, 2);
        let cam = best
            .modifications
            .iter()
            .find(|m| m.record_id == 4)
            .unwrap();
        assert!(cam.fixed);
        assert_eq!(cam.count, 1);
        let phospho = best
            .modifications
            .iter()
            .find(|m| m.record_id == 21)
            .unwrap();
        assert!(!phospho.fixed);
        assert!((best.theoretical_delta_mass - (57.021464 + 79.966331)).abs() < 1e-9);
    }

    #[test]
    fn solve_fixed_only_budget_matches_the_pinned_multiset_alone() {
        let mut pool = records();
        pool.push(cam_record());
        let query = SolveQuery {
            target_mass: 2.0 * 57.021464,
            tolerance: Tol {
                unit: TolUnit::Da,
                value: 0.001,
            },
            max_modifications: 2,
            max_distinct_modifications: None,
            allow_repeated_modification: true,
            approved: None,
            allowed_sites: Vec::new(),
            fixed_modifications: vec![FixedMod {
                record_id: 4,
                count: 2,
            }],
            limit: 10,
        };
        let result = solve_combinations(&pool, &[], &[], &query).expect("solve succeeds");
        assert_eq!(result.match_count, 1);
        assert_eq!(result.completeness, "complete");
        assert_eq!(result.candidates[0].modifications.len(), 1);
        assert_eq!(result.candidates[0].modification_count, 2);
        assert!(result.candidates[0].modifications[0].fixed);

        // A target the pinned multiset cannot cover within tolerance yields
        // zero candidates, not an error.
        let unmatched = SolveQuery {
            target_mass: 2.0 * 57.021464 + 0.5,
            ..query.clone()
        };
        let result = solve_combinations(&pool, &[], &[], &unmatched).expect("solve succeeds");
        assert_eq!(result.match_count, 0);
        assert!(result.candidates.is_empty());
    }

    #[test]
    fn solve_fixed_unknown_record_is_a_query_error() {
        let query = SolveQuery {
            target_mass: 100.0,
            tolerance: Tol {
                unit: TolUnit::Da,
                value: 0.01,
            },
            max_modifications: 2,
            max_distinct_modifications: None,
            allow_repeated_modification: true,
            approved: None,
            allowed_sites: Vec::new(),
            fixed_modifications: vec![FixedMod {
                record_id: 4242,
                count: 1,
            }],
            limit: 10,
        };
        let result = solve_combinations(&records(), &[], &[], &query);
        assert!(
            result.is_err(),
            "a pinned id absent from the pool is a broken query, not an empty answer"
        );
    }

    #[test]
    fn solve_fixed_does_not_consume_allow_repeated() {
        let mut pool = records();
        pool.push(cam_record());
        let query = SolveQuery {
            target_mass: 2.0 * 57.021464,
            tolerance: Tol {
                unit: TolUnit::Da,
                value: 0.001,
            },
            max_modifications: 2,
            max_distinct_modifications: None,
            allow_repeated_modification: false,
            approved: None,
            allowed_sites: Vec::new(),
            fixed_modifications: vec![FixedMod {
                record_id: 4,
                count: 1,
            }],
            limit: 10,
        };
        let result = solve_combinations(&pool, &[], &[], &query).expect("solve succeeds");
        assert_eq!(
            result.match_count, 1,
            "a pinned instance is not a searched repeat"
        );
        let mods = &result.candidates[0].modifications;
        assert_eq!(mods.len(), 1);
        assert_eq!(
            mods[0].count, 2,
            "one pinned + one searched instance of CAM"
        );
        assert!(mods[0].fixed);
    }

    #[test]
    fn solve_fixed_instances_consume_allowed_sites() {
        let mut pool = records();
        pool.push(cam_record());
        let specificities = vec![spec(4, "K"), spec(21, "S"), spec(21, "T")];
        let base = SolveQuery {
            target_mass: 57.021464 + 79.966331, // CAM + Phospho
            tolerance: Tol {
                unit: TolUnit::Da,
                value: 0.001,
            },
            max_modifications: 2,
            max_distinct_modifications: None,
            allow_repeated_modification: true,
            approved: None,
            allowed_sites: vec![SolveSite {
                residue: "K".to_owned(),
                position: None,
            }],
            fixed_modifications: vec![
                FixedMod {
                    record_id: 4,
                    count: 1,
                },
                FixedMod {
                    record_id: 21,
                    count: 1,
                },
            ],
            limit: 10,
        };
        // One K site: the pinned CAM takes it, the pinned Phospho (S/T only)
        // has nowhere to go.
        let result = solve_combinations(&pool, &specificities, &[], &base).expect("solve succeeds");
        assert_eq!(
            result.match_count, 0,
            "pinned instances occupy allowed sites"
        );

        let result = solve_combinations(
            &pool,
            &specificities,
            &[],
            &SolveQuery {
                allowed_sites: vec![
                    SolveSite {
                        residue: "K".to_owned(),
                        position: None,
                    },
                    SolveSite {
                        residue: "S".to_owned(),
                        position: None,
                    },
                ],
                ..base.clone()
            },
        )
        .expect("solve succeeds");
        assert_eq!(result.match_count, 1);
        assert!(result.candidates[0].site_feasible);
    }

    #[test]
    fn site_filter_on_mass_search_requires_a_specificity_site() {
        let query = SearchQuery {
            mass: Some(79.966331),
            tolerance: Some(Tol {
                unit: TolUnit::Da,
                value: 0.01,
            }),
            filters: UniModFilters {
                residues: vec!["Y".to_owned()], // Phospho sites are S/T → not feasible
                ..Default::default()
            },
            ..Default::default()
        };
        let matches = search_by_mass(&records(), &[spec(21, "S"), spec(21, "T")], &query);
        assert_eq!(matches.len(), 1);
        assert!(!matches[0].feasible, "Y not among S/T sites");
    }

    #[test]
    fn classification_gate_is_record_level_keeps_oxidation() {
        // Oxidation mirrors real UniMod: approved=false, with a Met=Artefact
        // specificity alongside Chemical-derivative sites. A record-level gate
        // keep=[Artefact, Post-translational] must KEEP it (a per-site filter
        // would drop it on the Chemical-derivative sites — the inversion bug).
        // A pure AA-substitution record at the same mass must be DROPPED.
        let mut recs = records();
        recs.push(UniModRecord {
            record_id: 35,
            title: "Oxidation".to_owned(),
            full_name: "Oxidation".to_owned(),
            approved: false,
            delta_mono: 15.9949,
            delta_average: 15.999,
            composition: "O1".to_owned(),
            source: "unimod".to_owned(),
        });
        recs.push(UniModRecord {
            record_id: 999,
            title: "X->Y".to_owned(),
            full_name: "substitution".to_owned(),
            approved: true,
            delta_mono: 15.9949,
            delta_average: 15.999,
            composition: "O1".to_owned(),
            source: "unimod".to_owned(),
        });
        let mut ox_m = spec(35, "M");
        ox_m.classification = "Artefact".to_owned();
        let mut ox_t = spec(35, "T");
        ox_t.classification = "Chemical derivative".to_owned();
        let mut subst = spec(999, "A");
        subst.classification = "AA substitution".to_owned();
        let query = SearchQuery {
            mass: Some(15.9949),
            tolerance: Some(Tol {
                unit: TolUnit::Da,
                value: 0.01,
            }),
            filters: UniModFilters {
                keep_classification: vec!["Artefact".to_owned(), "Post-translational".to_owned()],
                ..Default::default()
            },
            ..Default::default()
        };
        let matches = search_by_mass(&recs, &[ox_m, ox_t, subst], &query);
        let titles: Vec<&str> = matches.iter().map(|m| m.title.as_str()).collect();
        assert!(
            titles.contains(&"Oxidation"),
            "Oxidation kept (has an Artefact site): {titles:?}"
        );
        assert!(
            !titles.contains(&"X->Y"),
            "pure AA-substitution dropped: {titles:?}"
        );
    }
}
