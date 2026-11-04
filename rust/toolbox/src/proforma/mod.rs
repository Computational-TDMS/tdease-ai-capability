/*!
ProForma — the single chemical-truth source for the crate.

Every mass-spec API that needs an exact elemental composition (isotope
envelope, mass-shift calibration, fragment matching) resolves its chemistry
through this module. There is deliberately no second parser: a ProForma must
produce the same exact composition and monoisotopic mass everywhere in the
crate.

Supported syntax (a deliberate, explicit subset):
  - single-letter amino-acid residues
  - `[Formula:CHNOSPBr...]` — an explicit net formula; element counts may be
    signed (`Br-H` = Br replaces an aromatic H), so substitutions are
    expressed as net changes, not as `+Br`.
  - `[+42.0106]` / `[-18.0106]` — a numeric monoisotopic mass delta. Because a
    bare number carries no element information, the composition cannot be
    exact; a warning is emitted and the delta mass is used directly.
  - named modifications and ontology IDs resolved by mzcore's embedded
    ontologies; there is no local name-to-formula override table.
  - caller-supplied extra names via [`parse_proforma`]'s `extra_mods` slice.
    These explicit formulas override ontology lookup. Unknown names hard-error.
  - `[Mod@site]` — an optional `@site` suffix (`[Acetyl@K4]`, `[UNIMOD:7@K]`),
    recorded as the modification's position.
  - net formulas accept **two real-world spellings** with the same meaning:
    `Br-H` / `O-N` / `C2H2O` (net syntax) and `Br H -1` / `H -1 N -1 O 1`
    (UniMod composition syntax). Both reduce to the same signed (symbol, count)
    terms, so a UniMod record's composition feeds the exact-composition path
    with zero normalization.
*/

use std::num::NonZeroU16;

use chemical_elements::ChemicalComposition;
use mzcore::chemistry::{Element, MolecularFormula};
use serde::{Deserialize, Serialize};

use crate::error::MassSpecError;
use crate::fragment::residue_mass;
use crate::isotopic::ElementCount;

mod mzcore_adapter;
mod preprocess;
pub(crate) use preprocess::normalize_formula;

/// An externally supplied modification name -> net formula, for names the
/// embedded ontology does not know (e.g. from a user database).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ExternalModification {
    pub name: String,
    /// Net formula with signed counts, e.g. "Br-H" or "C2H2O".
    pub formula: String,
}

/// One parsed modification tag on the sequence.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ParsedModification {
    /// The modification tag as it appeared in the ProForma (name only; an
    /// `@site` suffix is parsed into `position`), e.g. `Acetyl` or `Br-H`.
    pub tag: String,
    /// The monoisotopic delta mass this modification contributes (Da).
    pub delta_mass: f64,
    /// The exact net formula, when the tag carried one (either `Formula:`,
    /// an ontology name, or an injected extra). `None` for a bare numeric tag.
    pub formula: Option<String>,
    /// The modification's site: `"nTerm"`, `"cTerm"`, `"residue:<N>"` (1-based),
    /// or `"site:<AA>"` (a bare residue letter, position unknown). Inferred from
    /// placement, or from an explicit `[Mod@site]` suffix when present.
    pub position: Option<String>,
    /// The residue letter this modification sits on, when the position pins a
    /// single 1-based residue (`residue:<N>`). Every numeric position is thus
    /// self-describing — an off-by-one is mechanically visible at the parse
    /// boundary (e.g. `87:V` when the manifest intended Y88) instead of a
    /// silent convention mismatch between artifact layers.
    pub residue: Option<char>,
    /// Inclusive 1-based TopPIC localization interval. Equal bounds are an
    /// exact residue; wider bounds retain the search engine's ambiguity.
    pub localization_start: Option<usize>,
    pub localization_end: Option<usize>,
}

/// The parsed truth: exact composition + exact monoisotopic mass + provenance.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ParsedProForma {
    /// The residue sequence as parsed (tags removed).
    pub sequence: String,
    /// Exact elemental composition, including water and all net changes.
    /// Counts may be negative for substitutions (e.g. `Bromo` -> Br + H -1).
    pub composition: Vec<ElementCount>,
    /// The exact monoisotopic neutral mass (residues + water + deltas).
    pub monoisotopic_mass: f64,
    /// The parsed modifications, in order of appearance.
    pub modifications: Vec<ParsedModification>,
    /// Warnings (e.g. a bare numeric delta has no exact composition).
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FragmentationProForma {
    pub variants: Vec<String>,
    pub ambiguous_ranges: Vec<(usize, usize)>,
}

/// Parse a ProFormA into the exact composition and monoisotopic mass.
///
/// The parsing backend is mzcore (full ProFormA 2.0/2.1, embedded
/// UniMod/PSI-MOD/GNOme ontologies). This crate layers on top:
///
/// - TopPIC compatibility: `M.(S)[Acetyl]...K.K` flanking dots are stripped by
///   [`preprocess`].
/// - The chemistry truth is mzcore's whole-peptide `MolecularFormula`: both
///   `composition` and `monoisotopic_mass` are projected verbatim from it
///   ([`mzcore_adapter`]), and every modification's `delta_mass` is mzcore's
///   formula mass — there is no second element/mass table here.
/// - `external_mods` supplies explicit `(name, net_formula)` overrides.
///   Otherwise mzcore resolves names and IDs; unknown names are hard errors.
pub fn parse_proforma(
    text: &str,
    external_mods: &[ExternalModification],
) -> Result<ParsedProForma, MassSpecError> {
    // 1. TopPIC -> ProFormA normalization: strip flanking dots, rewrite tags.
    let mut normalized = preprocess::normalize(text, external_mods)?;
    let concrete = preprocess::concretize_mod_ranges(&normalized.mzcore_text)?;
    normalized.mzcore_text = concrete.into_iter().next().ok_or_else(|| {
        MassSpecError::InvalidQuery("ProForma produced no concrete localization".to_owned())
    })?;

    // 2. mzcore parse + chemistry reduction. mzcore's `MolecularFormula` is the
    //    single chemical truth: the whole-peptide formula supplies both the
    //    composition and the monoisotopic mass, and each modification carries
    //    its own mzcore-computed delta mass.
    let (chemistries, whole_formula) = mzcore_adapter::parse_with_mzcore(&normalized)?;

    // 3. Project mzcore's chemistry onto our flat model.
    let composition = formula_to_element_counts(&whole_formula);
    let monoisotopic_mass = whole_formula.monoisotopic_mass().value;
    let mut modifications = Vec::new();
    let mut warnings = Vec::new();

    for (tag, chemistry) in normalized.tags.iter().zip(&chemistries) {
        if chemistry.counts.is_empty() {
            // Mass-only (or INFO) — no element contribution.
            warnings.push(format!(
                "mass-only modification [{}] has no exact elemental composition",
                tag.tag
            ));
        }
        modifications.push(ParsedModification {
            residue: tag.position.as_deref().and_then(|position| {
                let index = position.strip_prefix("residue:").unwrap_or(position);
                index
                    .parse::<usize>()
                    .ok()
                    .and_then(|index| normalized.sequence.chars().nth(index.checked_sub(1)?))
            }),
            tag: tag.tag.clone(),
            delta_mass: chemistry.delta_mass,
            formula: if chemistry.counts.is_empty() {
                None
            } else {
                Some(counts_to_net_formula(&chemistry.counts))
            },
            position: tag.position.clone(),
            localization_start: tag.range.map(|range| range.0),
            localization_end: tag.range.map(|range| range.1),
        });
    }

    Ok(ParsedProForma {
        sequence: normalized.sequence,
        composition,
        monoisotopic_mass,
        modifications,
        warnings,
    })
}

/// The exact monoisotopic delta of a named modification, computed from its
/// formula by the same chemistry as [`parse_proforma`] (builtin table or
/// mzcore's ontologies). Interfaces accept names so that chemistry-known
/// masses are never hand-typed (E64): `Bromo` = Br−H = 77.910513,
/// `Phospho` = HO3P = 79.966331, `Dimethyl` = C2H4 = 28.031300.
pub fn mod_delta_mass(name: &str) -> Result<f64, MassSpecError> {
    let parsed = parse_proforma(&format!("G[{name}]"), &[])?;
    parsed
        .modifications
        .first()
        .filter(|m| m.tag.eq_ignore_ascii_case(name))
        .map(|m| m.delta_mass)
        .ok_or_else(|| {
            MassSpecError::InvalidQuery(format!(
                "modification name {name:?} resolved to no delta mass"
            ))
        })
}

/// Normalize TopPIC flanks and external modification names into the exact
/// HUPO-PSI ProForma expression consumed by mzcore/mzannotate.
pub fn normalize_for_fragmentation(
    text: &str,
    external_mods: &[ExternalModification],
) -> Result<FragmentationProForma, MassSpecError> {
    let normalized = preprocess::normalize(text, external_mods)?;
    let ambiguous_ranges = normalized
        .tags
        .iter()
        .filter_map(|tag| tag.range)
        .filter(|(start, end)| start < end)
        .collect();
    Ok(FragmentationProForma {
        variants: preprocess::concretize_mod_ranges(&normalized.mzcore_text)?,
        ambiguous_ranges,
    })
}

/// Flatten an mzcore `MolecularFormula` into the crate's `ElementCount` list.
///
/// mzcore represents a charged formula by a negative `Element::Electron` entry
/// (the charge carrier). Those are skipped here: the charge is carried
/// separately by the caller (`theoretical_envelope`'s `charge` parameter), and
/// `chemical_elements` has no "e" element — feeding it through `inc_str` would
/// panic on `UnknownElement`.
pub(crate) fn formula_to_element_counts(formula: &MolecularFormula) -> Vec<ElementCount> {
    formula
        .elements()
        .iter()
        .filter(|(element, _, _)| *element != Element::Electron)
        .map(|(element, isotope, count)| ElementCount {
            element: element.symbol().to_owned(),
            count: *count,
            isotope: isotope.map(NonZeroU16::get).unwrap_or(0),
        })
        .collect()
}

/// Render an `ElementCount` list back into a canonical net-formula string,
/// e.g. `C2H2O`, `BrH-1`. Positive counts are `SYM{n}` (n omitted for 1);
/// negative counts are `SYM-{n}`. This is the `formula` field of a parsed
/// modification.
fn counts_to_net_formula(counts: &[ElementCount]) -> String {
    let mut sorted = counts.to_vec();
    sorted.sort_by(|a, b| a.element.cmp(&b.element));
    let mut out = String::new();
    for entry in sorted {
        match entry.count {
            1 => out.push_str(&entry.element),
            n => {
                out.push_str(&entry.element);
                out.push_str(&n.to_string());
            }
        }
    }
    if out.is_empty() {
        // No element contributes (shouldn't happen for a non-empty counts list,
        // but keep the output well-formed).
        "0".to_owned()
    } else {
        out
    }
}

/// Parse an `@site` suffix into a position string.
///
/// Accepts `K4`/`4` (→ `residue:<N>`), `K` (→ `site:<AA>`), `N-term`/`C-term`
/// (→ `nTerm`/`cTerm`), and returns `unknown` for anything else.
pub(crate) fn parse_site(site: &str) -> String {
    let s = site.trim();
    if s.eq_ignore_ascii_case("nterm") || s.eq_ignore_ascii_case("n-term") {
        return "nTerm".to_owned();
    }
    if s.eq_ignore_ascii_case("cterm") || s.eq_ignore_ascii_case("c-term") {
        return "cTerm".to_owned();
    }
    let mut letters = String::new();
    let mut digits = String::new();
    let mut known = true;
    for c in s.chars() {
        if c.is_ascii_alphabetic() && digits.is_empty() {
            letters.push(c);
        } else if c.is_ascii_digit() {
            digits.push(c);
        } else {
            known = false;
        }
    }
    if !known || (letters.is_empty() && digits.is_empty()) {
        "unknown".to_owned()
    } else if !digits.is_empty() {
        format!("residue:{}", digits)
    } else {
        format!("site:{}", letters.to_uppercase())
    }
}

/// Compute the monoisotopic mass of a net formula like `C2H2O` or `Br-H`.
#[cfg(test)]
fn net_formula_mass(formula: &str) -> Result<f64, MassSpecError> {
    let terms = split_net_formula(formula)?;
    let mut mass = 0.0;
    for (symbol, count) in terms {
        let elt_mass = element_monoisotopic_mass(&symbol).ok_or_else(|| {
            MassSpecError::InvalidQuery(format!(
                "unsupported element {symbol} in net formula [{formula}]"
            ))
        })?;
        mass += elt_mass * f64::from(count);
    }
    Ok(mass)
}

/// Split a net formula into signed `(symbol, count)` terms.
///
/// Accepts a **superset of the two real-world spellings**, reducing both to the
/// same terms:
///
/// | input | terms | source |
/// |---|---|---|
/// | `C2H2O` | C2 H2 O1 | both |
/// | `Br-H` / `O-N` | Br1 H-1 / O1 N-1 | net |
/// | `-H+Br` | H-1 Br1 | net |
/// | `Br H -1` | Br1 H-1 | UniMod |
/// | `H -1 N -1 O 1` | H-1 N-1 O1 | UniMod |
/// | `S -1 Se 1` | S-1 Se1 | UniMod (substitution) |
/// | `H(-1) Br` | H-1 Br1 | UniMod composition (signed count in parens) |
/// | `S(1) O(3)` | S1 O3 | UniMod composition |
///
/// Whitespace is allowed between any tokens. The sign rules: a `+`/`-` followed
/// by a digit is the signed count of the *current* element; a `+`/`-` followed
/// by an element starts a new element with that sign; a digit with no element
/// is an error; a trailing sign is an error.
pub fn split_net_formula(formula: &str) -> Result<Vec<(String, i32)>, MassSpecError> {
    let chars: Vec<char> = formula.chars().collect();
    let mut i = 0;
    let mut terms: Vec<(String, i32)> = Vec::new();
    // The open element: (symbol, 0) means its count is not yet read.
    let mut cur: Option<(String, i32)> = None;
    // The sign to apply to the next count or next element (default +1).
    let mut cur_sign = 1i32;

    while i < chars.len() {
        let c = chars[i];
        match c {
            '+' | '-' => {
                let sign = if c == '+' { 1 } else { -1 };
                // Find the next non-whitespace token (a sign may be separated
                // from its count or element by spaces).
                let mut j = i + 1;
                while j < chars.len() && chars[j].is_ascii_whitespace() {
                    j += 1;
                }
                match chars.get(j).copied() {
                    Some(next) if next.is_ascii_digit() => {
                        // Signed count of the current element: `H -1`, `H-2`.
                        let (sym, _) = cur.take().ok_or_else(|| {
                            MassSpecError::InvalidQuery(format!(
                                "malformed net formula [{formula}]: signed count with no element"
                            ))
                        })?;
                        let (count, consumed) = read_count(&chars, j, formula)?;
                        terms.push((sym, sign * count));
                        cur_sign = 1;
                        i = j + consumed;
                    }
                    _ => {
                        // The sign applies to the next element: `-H`, `Br-H`, `H - N`.
                        if let Some((sym, _)) = cur.take() {
                            terms.push((sym, cur_sign));
                        }
                        cur_sign = sign;
                        i += 1;
                        // A sign must be followed by something, else it dangles.
                        if !chars[i..].iter().any(|c| !c.is_ascii_whitespace()) {
                            return Err(MassSpecError::InvalidQuery(format!(
                                "malformed net formula [{formula}]: trailing sign"
                            )));
                        }
                    }
                }
            }
            c if c.is_ascii_uppercase() => {
                // A new element begins: flush the open one with its implicit count.
                if let Some((sym, _)) = cur.take() {
                    terms.push((sym, cur_sign));
                }
                let (sym, consumed) = read_element(&chars, i);
                cur = Some((sym, 0));
                // cur_sign is deliberately left untouched: it carries the sign
                // for this new element when a `+`/`-` preceded it.
                i += consumed;
            }
            c if c.is_ascii_digit() => {
                if let Some((sym, _)) = cur.take() {
                    // Count of the current element: `C2`, `H -1`.
                    let (count, consumed) = read_count(&chars, i, formula)?;
                    terms.push((sym, cur_sign * count));
                    cur_sign = 1;
                    i += consumed;
                } else {
                    // Stable-isotope prefix: `2H`, `13C`, `15N`, `18O` — a digit
                    // sequence directly followed by an element symbol. This is
                    // UniMod's labelling-modification spelling.
                    let (digits, consumed) = read_count(&chars, i, formula)?;
                    match chars.get(i + consumed).copied() {
                        Some(next) if next.is_ascii_uppercase() => {
                            let (elem, elem_consumed) = read_element(&chars, i + consumed);
                            let sym = format!("{digits}{elem}");
                            cur = Some((sym, 0));
                            i += consumed + elem_consumed;
                        }
                        _ => {
                            return Err(MassSpecError::InvalidQuery(format!(
                                "malformed net formula [{formula}]: digit before element"
                            )));
                        }
                    }
                }
            }
            '(' => {
                // UniMod's signed-count-in-parens spelling: `H(-1)` = H × -1,
                // `Br(1)` = Br × +1. The sign inside the parens is independent.
                let (sym, _) = cur.take().ok_or_else(|| {
                    MassSpecError::InvalidQuery(format!(
                        "malformed net formula [{formula}]: parens count with no element"
                    ))
                })?;
                let mut j = i + 1;
                let inner_sign = match chars.get(j).copied() {
                    Some('+') => {
                        j += 1;
                        1
                    }
                    Some('-') => {
                        j += 1;
                        -1
                    }
                    _ => 1,
                };
                let start = j;
                while j < chars.len() && chars[j].is_ascii_digit() {
                    j += 1;
                }
                if start == j || chars.get(j).copied() != Some(')') {
                    return Err(MassSpecError::InvalidQuery(format!(
                        "malformed net formula [{formula}]: expected (signed count) after {sym}"
                    )));
                }
                let count: i32 =
                    chars[start..j]
                        .iter()
                        .collect::<String>()
                        .parse()
                        .map_err(|_| {
                            MassSpecError::InvalidQuery(format!(
                                "invalid count in net formula [{formula}]"
                            ))
                        })?;
                terms.push((sym, inner_sign * count));
                cur_sign = 1;
                i = j + 1;
            }
            c if c.is_ascii_whitespace() => i += 1,
            c => {
                return Err(MassSpecError::InvalidQuery(format!(
                    "unsupported character {c} in net formula [{formula}]"
                )));
            }
        }
    }
    if let Some((sym, _)) = cur.take() {
        terms.push((sym, cur_sign));
    }
    if terms.is_empty() {
        return Err(MassSpecError::InvalidQuery(format!(
            "empty net formula [{formula}]"
        )));
    }
    Ok(terms)
}

/// Read an element symbol: an uppercase letter plus following lowercase letters
/// (`Br`, `Cl`, `Se`, ...). Returns `(symbol, chars_consumed)`.
fn read_element(chars: &[char], start: usize) -> (String, usize) {
    let mut i = start;
    let mut symbol = String::new();
    symbol.push(chars[i]);
    i += 1;
    while i < chars.len() && chars[i].is_ascii_lowercase() {
        symbol.push(chars[i]);
        i += 1;
    }
    (symbol, i - start)
}

/// Read a count: one or more digits. Returns `(count, chars_consumed)`.
fn read_count(chars: &[char], start: usize, formula: &str) -> Result<(i32, usize), MassSpecError> {
    let mut i = start;
    while i < chars.len() && chars[i].is_ascii_digit() {
        i += 1;
    }
    let text: String = chars[start..i].iter().collect();
    let count = text.parse::<i32>().map_err(|_| {
        MassSpecError::InvalidQuery(format!("invalid count in net formula [{formula}]"))
    })?;
    Ok((count, i - start))
}

/// Monoisotopic mass of supported elements (Da), used only to validate and
/// resolve `apply_net_formula` inputs (external-mod formula strings). The
/// parsed-composition masses no longer come through here — they are mzcore's
/// `MolecularFormula` masses.
///
/// The core proteomics elements use explicit IUPAC monoisotopic masses; any
/// other element resolves through `chemical_elements`' periodic table, where
/// the most-abundant-isotope mass is the monoisotopic mass (covering the full
/// UniMod element set: Se, I, Cl, Fe, Cu, Zn, Na, K, Mg, Ca, ...). An unknown
/// symbol returns None and is a hard error upstream (never a silent guess).
fn element_monoisotopic_mass(symbol: &str) -> Option<f64> {
    const EXPLICIT: [(&str, f64); 6] = [
        ("C", 12.0),
        ("H", 1.007_825_032_23),
        ("N", 14.003_074_004_43),
        ("O", 15.994_914_619_57),
        ("S", 31.972_071_174_4),
        ("P", 30.973_761_998_42),
    ];
    if let Some((_, mass)) = EXPLICIT.iter().find(|(symbol_, _)| *symbol_ == symbol) {
        return Some(*mass);
    }
    // Stable-isotope symbols that UniMod uses in labelling modifications. Their
    // masses are the isotope masses, not the natural-abundance element mass.
    const ISOTOPIC: [(&str, f64); 4] = [
        ("2H", 2.014_101_777_8),
        ("13C", 13.003_354_837_8),
        ("15N", 15.000_108_898_2),
        ("18O", 17.999_159_612_9),
    ];
    if let Some((_, mass)) = ISOTOPIC.iter().find(|(symbol_, _)| *symbol_ == symbol) {
        return Some(*mass);
    }
    chemical_elements::PERIODIC_TABLE
        .get(symbol)
        .map(|element| element.mass())
}

/// Apply a net formula (`C2H2O`, `Br-H`, `Br H -1`, `O-N`, `2H(8)`) to a
/// running composition.
///
/// Elements are validated before incrementing: `chemical_elements`' `inc_str`
/// **panics** on an unknown element (via `ElementSpecification::parse().unwrap()`),
/// so this checks the element table first and turns an unknown symbol into a
/// query error rather than a panic. Stable-isotope symbols (`2H`, `13C`) are
/// bridged to `chemical_elements`' bracket-isotope spelling (`H[2]`, `C[13]`).
pub(crate) fn apply_net_formula(
    composition: &mut ChemicalComposition<'_>,
    formula: &str,
) -> Result<(), MassSpecError> {
    let terms = split_net_formula(formula)?;
    for (symbol, count) in terms {
        if element_monoisotopic_mass(&symbol).is_none() {
            return Err(MassSpecError::InvalidQuery(format!(
                "unsupported element {symbol} in net formula [{formula}]"
            )));
        }
        // inc_str with a negative count subtracts.
        increment(&mut *composition, &symbol, count)?;
    }
    Ok(())
}

/// Increment one (possibly isotope-prefixed) element symbol into a composition.
fn increment(
    composition: &mut ChemicalComposition<'_>,
    symbol: &str,
    count: i32,
) -> Result<(), MassSpecError> {
    if let Some((base, isotope)) = split_isotope_symbol(symbol) {
        let spec = chemical_elements::ElementSpecification::parse(&format!("{base}[{isotope}]"))
            .map_err(|error| {
                MassSpecError::InvalidQuery(format!(
                    "unsupported isotope element {symbol}: {error:?}"
                ))
            })?;
        composition.inc(spec, count);
        Ok(())
    } else {
        composition.inc_str(symbol, count);
        Ok(())
    }
}

/// Split an isotope-prefixed symbol (`2H`, `13C`) into `(base, mass_number)`.
fn split_isotope_symbol(symbol: &str) -> Option<(&str, u16)> {
    let digits: String = symbol.chars().take_while(|c| c.is_ascii_digit()).collect();
    if digits.is_empty() {
        return None;
    }
    let rest = &symbol[digits.len()..];
    let mass: u16 = digits.parse().ok()?;
    if rest.is_empty() {
        return None;
    }
    Some((rest, mass))
}

/// Build a `ChemicalComposition` from `ElementCount` entries (signed counts
/// allowed — substitutions like `Bromo` = +Br -H are representable).
///
/// Public because native CLI adapters compose `parse_proforma` output
/// into `theoretical_envelope`'s input.
pub fn chemical_from_counts(counts: &[ElementCount]) -> ChemicalComposition<'static> {
    let mut composition = ChemicalComposition::new();
    for entry in counts {
        if entry.isotope != 0 {
            // Re-express the isotope as chemical_elements' bracket spelling so
            // the composition carries the isotope's own mass.
            let spec = chemical_elements::ElementSpecification::parse(&format!(
                "{}[{}]",
                entry.element, entry.isotope
            ))
            .expect("isotope counts come from parse_proforma and are valid");
            composition.inc(spec, entry.count);
        } else {
            let _ = increment(&mut composition, &entry.element, entry.count);
        }
    }
    composition
}

/// Flatten a `ChemicalComposition` into sorted `ElementCount` entries.
pub(crate) fn element_counts_from_chemical(
    composition: &ChemicalComposition<'_>,
) -> Vec<ElementCount> {
    element_counts(composition)
}

/// Flatten a `ChemicalComposition` into sorted `ElementCount` entries.
fn element_counts(composition: &ChemicalComposition<'_>) -> Vec<ElementCount> {
    let mut counts: Vec<ElementCount> = composition
        .iter()
        .filter(|(_, count)| **count != 0)
        .map(|(element, count)| ElementCount {
            element: element.element.symbol.clone(),
            count: *count,
            isotope: element.isotope,
        })
        .collect();
    counts.sort_by(|a, b| a.element.cmp(&b.element));
    counts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residue_echo_pins_each_named_modification_to_its_letter() {
        // The incident shape, miniature: intent "Br on the Y", placement one
        // residue early. The echo makes the mismatch mechanically visible at
        // the parse boundary (87:V vs 88:Y) instead of a silent convention
        // drift between artifact layers.
        let early = parse_proforma("[Acetyl]-SGRGKDVV[Bromo]YALKR", &[]).unwrap();
        let br = early
            .modifications
            .iter()
            .find(|m| m.tag == "Bromo")
            .expect("Bromo tag parsed");
        assert_eq!(br.position.as_deref(), Some("residue:8"));
        assert_eq!(br.residue, Some('V'), "the early placement echoes V, not Y");
        let acetyl = early
            .modifications
            .iter()
            .find(|m| m.tag == "Acetyl")
            .expect("Acetyl tag parsed");
        assert_eq!(acetyl.position.as_deref(), Some("nTerm"));
        assert_eq!(acetyl.residue, None, "terminal mods carry no residue echo");

        // One residue later — the intended shape — echoes Y.
        let intended = parse_proforma("[Acetyl]-SGRGKDVVY[Bromo]ALKR", &[]).unwrap();
        let br = intended
            .modifications
            .iter()
            .find(|m| m.tag == "Bromo")
            .unwrap();
        assert_eq!(br.position.as_deref(), Some("residue:9"));
        assert_eq!(br.residue, Some('Y'));
    }

    fn count_of(composition: &[ElementCount], element: &str) -> i32 {
        composition
            .iter()
            .find(|e| e.element == element)
            .map(|e| e.count)
            .unwrap_or(0)
    }

    // ── TopPIC flanking compatibility ────────────────────────────────────────

    #[test]
    fn toppic_flanking_sequence_is_stripped() {
        // TopPIC proteoform format: `M.` N-terminal leader and `.K` C-terminal
        // trailer separated by dots. The flanks carry no mass and are stripped;
        // the `(S)[Acetyl]` modRange is standard ProFormA.
        let p = parse_proforma("M.(S)[Acetyl]ESSSKSSQPLASKQEKDGTEKRGRGRPRKQPPK.K", &[]).unwrap();
        assert_eq!(
            p.sequence, "SESSSKSSQPLASKQEKDGTEKRGRGRPRKQPPK",
            "sequence without flanks, with (S) residue kept"
        );
        let m = &p.modifications[0];
        assert_eq!(m.tag, "Acetyl");
        assert!(
            (m.delta_mass - 42.0106).abs() < 1e-3,
            "delta {}",
            m.delta_mass
        );
        assert_eq!(m.position.as_deref(), Some("site:S"));
        assert!(p.monoisotopic_mass > 0.0);
    }

    #[test]
    fn toppic_leading_only_flank() {
        let p = parse_proforma("M.PEPTIDE", &[]).unwrap();
        assert_eq!(p.sequence, "PEPTIDE");
    }

    #[test]
    fn toppic_trailing_only_flank() {
        let p = parse_proforma("PEPTIDE.K", &[]).unwrap();
        assert_eq!(p.sequence, "PEPTIDE");
    }

    #[test]
    fn toppic_multiple_ranges() {
        let p = parse_proforma("M.ESS(S)[Acetyl]K(T)[Phospho]PK.K", &[]).unwrap();
        assert_eq!(p.sequence, "ESSSKTPK");
        assert_eq!(p.modifications.len(), 2);
        assert_eq!(p.modifications[0].tag, "Acetyl");
        assert_eq!(p.modifications[1].tag, "Phospho");
    }

    #[test]
    fn middle_dot_is_rejected() {
        // A dot inside the sequence (not a flank) must error, not be silently
        // dropped.
        assert!(parse_proforma("PE.PTIDE", &[]).is_err());
    }

    // ── Embedded-ontology resolution (new capability) ──────────────────────

    #[test]
    fn bare_name_resolves_through_embedded_ontology() {
        // Carboxymethyl was previously a hard error; now it resolves through
        // UniMod.
        let p = parse_proforma("[Carboxymethyl]-PEPTIDE", &[]).unwrap();
        assert_eq!(p.modifications[0].tag, "Carboxymethyl");
        assert!(
            (p.modifications[0].delta_mass - 58.005).abs() < 0.01,
            "delta {}",
            p.modifications[0].delta_mass
        );
    }

    // ── Glycan support ──────────────────────────────────────────────────────

    #[test]
    fn glycan_composition_resolves() {
        // [Glycan:HexNAc1Hex2] — one N-acetylhexosamine + two hexoses.
        let p = parse_proforma("GAG[Glycan:HexNAc1Hex2]AGF", &[]).unwrap();
        let m = &p.modifications[0];
        assert_eq!(m.tag, "Glycan:HexNAc1Hex2");
        assert!(m.delta_mass > 0.0, "glycan delta {}", m.delta_mass);
        // HexNAc = C8H15NO6, Hex = C6H12O6 → total C20 H39 NO18 minus water
        // per monosaccharide; just assert the glycan adds carbon/nitrogen.
        assert!(
            count_of(&p.composition, "N")
                > count_of(&parse_proforma("GAGAGF", &[]).unwrap().composition, "N")
        );
    }

    #[test]
    fn info_tag_is_zero_mass() {
        // [INFO:...] contributes no mass or composition.
        let p = parse_proforma("[INFO:flank]-PEPTIDE", &[]).unwrap();
        let base = parse_proforma("PEPTIDE", &[]).unwrap();
        assert!((p.monoisotopic_mass - base.monoisotopic_mass).abs() < 1e-6);
    }

    #[test]
    fn named_acetyl_resolves_to_exact_net_formula() {
        let p = parse_proforma("[Acetyl]-PEPTIDE", &[]).unwrap();
        // PEPTIDE = P + E + P + T + I + D + E (7 residues) + water.
        // Acetyl = C2H2O must appear in the composition.
        assert!(p.modifications.iter().any(|m| m.tag == "Acetyl"));
        // PEPTIDE residue carbons: P=5 E=5 P=5 T=4 I=6 D=4 E=5 = 34 + Acetyl C2.
        assert_eq!(count_of(&p.composition, "C"), 34 + 2);
        // monoisotopic mass: residues + H2O + acetyl delta.
        let expected_seq = 97.052763875
            + 129.042593135
            + 97.052763875
            + 101.047678505
            + 113.084063975
            + 115.026943065
            + 129.042593135;
        assert!(
            (p.monoisotopic_mass - (expected_seq + crate::constants::H2O + 42.010565)).abs() < 1e-3
        );
    }

    #[test]
    fn bromo_is_net_change_br_minus_h() {
        let p = parse_proforma("[Bromo]-PEPTIDE", &[]).unwrap();
        let br = count_of(&p.composition, "Br");
        let h = count_of(&p.composition, "H");
        assert_eq!(br, 1, "Bromo must add exactly one Br");
        // H count: residues' H + water H2 + (Bromo replaces one H -> -1).
        // We just assert Br is present and H reflects the substitution.
        let p_plain = parse_proforma("PEPTIDE", &[]).unwrap();
        let h_plain = count_of(&p_plain.composition, "H");
        assert_eq!(h, h_plain - 1, "Bromo must replace one H (net Br-H)");
    }

    #[test]
    fn deamidation_expresses_oxygen_substituting_nitrogen() {
        let p = parse_proforma("[Deamidation]-PEPTIDE", &[]).unwrap();
        let p_plain = parse_proforma("PEPTIDE", &[]).unwrap();
        let o = count_of(&p.composition, "O");
        let n = count_of(&p.composition, "N");
        let o_plain = count_of(&p_plain.composition, "O");
        let n_plain = count_of(&p_plain.composition, "N");
        assert_eq!(o, o_plain + 1, "Deamidation adds one O");
        assert_eq!(n, n_plain - 1, "Deamidation removes one N (O-N)");
        assert_eq!(
            count_of(&p.composition, "H"),
            count_of(&p_plain.composition, "H") - 1
        );
        let canonical = parse_proforma("[UNIMOD:7]-PEPTIDE", &[]).unwrap();
        assert_eq!(p.composition, canonical.composition);
        assert_eq!(p.monoisotopic_mass, canonical.monoisotopic_mass);
    }

    #[test]
    fn numeric_delta_warns_but_keeps_mass() {
        let p = parse_proforma("[+42.0106]-PEPTIDE", &[]).unwrap();
        assert!(p.warnings.iter().any(|w| w.contains("mass-only")));
        assert!((p.modifications[0].delta_mass - 42.0106).abs() < 1e-4);
    }

    #[test]
    fn formula_tag_applies_net_formula() {
        let p = parse_proforma("[Formula:C2H2O]-PEPTIDE", &[]).unwrap();
        assert!(
            count_of(&p.composition, "C")
                > count_of(&parse_proforma("PEPTIDE", &[]).unwrap().composition, "C")
        );
        assert_eq!(p.modifications[0].formula.as_deref(), Some("C2H2O"));
    }

    #[test]
    fn unknown_name_hard_errors_unless_injected() {
        let err = parse_proforma("[MyCustom]-PEPTIDE", &[]).unwrap_err();
        assert!(err.to_string().contains("unknown ProForma modification"));

        let p = parse_proforma(
            "[MyCustom]-PEPTIDE",
            &[ExternalModification {
                name: "MyCustom".to_owned(),
                formula: "C2H2O".to_owned(),
            }],
        )
        .unwrap();
        assert!(p.modifications.iter().any(|m| m.tag == "MyCustom"));
        assert_eq!(
            count_of(&p.composition, "C"),
            count_of(
                &parse_proforma("[Acetyl]-PEPTIDE", &[]).unwrap().composition,
                "C"
            )
        );
    }

    #[test]
    fn injected_override_wins_over_builtin() {
        let p = parse_proforma(
            "[Acetyl]-PEPTIDE",
            &[ExternalModification {
                name: "Acetyl".to_owned(),
                formula: "C2H4".to_owned(),
            }],
        )
        .unwrap();
        let m = p.modifications.iter().find(|m| m.tag == "Acetyl").unwrap();
        assert_eq!(m.formula.as_deref(), Some("C2H4"));
        assert!(
            (m.delta_mass - 28.0313).abs() < 1e-3,
            "overridden Acetyl delta {}",
            m.delta_mass
        );
    }

    #[test]
    fn sequence_mass_matches_residue_sum() {
        let p = parse_proforma("PEPTIDE", &[]).unwrap();
        let expected = 97.052763875
            + 129.042593135
            + 97.052763875
            + 101.047678505
            + 113.084063975
            + 115.026943065
            + 129.042593135
            + crate::constants::H2O;
        assert!((p.monoisotopic_mass - expected).abs() < 1e-3);
    }

    #[test]
    fn monoisotopic_mass_matches_mzcore_whole_formula() {
        // The single chemical truth is mzcore's whole-peptide `MolecularFormula`:
        // `parse_proforma` must project its mass and composition verbatim, with
        // no second mass source of its own. Re-deriving the formula directly from
        // mzcore and comparing is the guard against a future dual-source drift.
        use mzcore::chemistry::{AmbiguousMolecule, OutputMolecularFormula};
        use mzcore::ontology::STATIC_ONTOLOGIES;
        use mzcore::prelude::SequencePosition;
        use mzcore::sequence::Peptidoform;
        for proforma in ["PEPTIDE", "[Acetyl]-PEPTIDE", "[Bromo]-PEPTIDE"] {
            let p = parse_proforma(proforma, &[]).unwrap();
            let (pf, _) = Peptidoform::pro_forma(proforma, &STATIC_ONTOLOGIES).unwrap();
            let linear = pf.into_linear().unwrap();
            let whole: MolecularFormula = linear
                .calculate_masses_inner::<OutputMolecularFormula>(SequencePosition::default(), 0, 0)
                .single()
                .expect("one whole-peptide formula");
            assert!(
                (p.monoisotopic_mass - whole.monoisotopic_mass().value).abs() < 1e-6,
                "{proforma}: reported {} vs mzcore {}",
                p.monoisotopic_mass,
                whole.monoisotopic_mass().value
            );
            assert_eq!(
                formula_to_element_counts(&whole),
                p.composition,
                "{proforma}"
            );
        }
    }

    // ---- unified net-formula syntax ----------------------------------------

    #[test]
    fn net_formula_spellings_agree() {
        // The three spellings of a bromine substitution must be identical.
        let net = split_net_formula("Br-H").unwrap();
        let unimod = split_net_formula("Br H -1").unwrap();
        let sign_first = split_net_formula("-H+Br").unwrap();
        // Order may differ across spellings (the sign can precede the element);
        // the *multiset* of terms is the invariant.
        let mut net_sorted = net.clone();
        let mut unimod_sorted = unimod.clone();
        let mut sign_first_sorted = sign_first.clone();
        net_sorted.sort();
        unimod_sorted.sort();
        sign_first_sorted.sort();
        assert_eq!(net_sorted, unimod_sorted);
        assert_eq!(net_sorted, sign_first_sorted);
        assert_eq!(net, vec![("Br".to_owned(), 1), ("H".to_owned(), -1)]);
        // And the masses agree to the Da the mass table carries.
        let m_net = net_formula_mass("Br-H").unwrap();
        let m_unimod = net_formula_mass("Br H -1").unwrap();
        let m_sign_first = net_formula_mass("-H+Br").unwrap();
        assert!((m_net - m_unimod).abs() < 1e-9);
        assert!((m_net - m_sign_first).abs() < 1e-9);
        assert!((m_net - 77.91051).abs() < 1e-4, "Br-H mono {}", m_net);
    }

    #[test]
    fn net_formula_unimod_deamidation() {
        let terms = split_net_formula("H -1 N -1 O 1").unwrap();
        assert_eq!(
            terms,
            vec![
                ("H".to_owned(), -1),
                ("N".to_owned(), -1),
                ("O".to_owned(), 1),
            ]
        );
    }

    #[test]
    fn net_formula_unimod_parens_syntax() {
        // UniMod stores composition as `H(-1) Br` (signed count in parens).
        let terms = split_net_formula("H(-1) Br").unwrap();
        assert_eq!(terms, vec![("H".to_owned(), -1), ("Br".to_owned(), 1)]);
        let s = split_net_formula("S(1) O(3)").unwrap();
        assert_eq!(s, vec![("S".to_owned(), 1), ("O".to_owned(), 3)]);
        // And it must agree with the spaced spelling and the net spelling.
        let a = split_net_formula("H(-1) Br").unwrap();
        let b = split_net_formula("H -1 Br 1").unwrap();
        let c = split_net_formula("-H+Br").unwrap();
        let mut s1 = a.clone();
        s1.sort();
        let mut s2 = b.clone();
        s2.sort();
        let mut s3 = c.clone();
        s3.sort();
        assert_eq!(s1, s2);
        assert_eq!(s1, s3);
    }

    #[test]
    fn net_formula_unimod_corpus() {
        // Representative UniMod composition strings (signed counts, spaces).
        let corpus: &[(&str, &[(&str, i32)])] = &[
            ("C2 H2 O 1", &[("C", 2), ("H", 2), ("O", 1)]),
            ("H1 O3 P1", &[("H", 1), ("O", 3), ("P", 1)]),
            ("H -1 N -1 O 1", &[("H", -1), ("N", -1), ("O", 1)]),
            ("O 1", &[("O", 1)]),
            ("C1 H2", &[("C", 1), ("H", 2)]),
            ("C2 H3 N1 O1", &[("C", 2), ("H", 3), ("N", 1), ("O", 1)]),
            ("H1 O3 S1", &[("H", 1), ("O", 3), ("S", 1)]),
            ("S -1 Se 1", &[("S", -1), ("Se", 1)]),
            ("I1", &[("I", 1)]),
            ("Cl1", &[("Cl", 1)]),
            ("Fe1", &[("Fe", 1)]),
            (
                "H -3 C -1 N -1 O 1",
                &[("H", -3), ("C", -1), ("N", -1), ("O", 1)],
            ),
            ("C -1 H -2", &[("C", -1), ("H", -2)]),
        ];
        for (formula, expected) in corpus {
            let got =
                split_net_formula(formula).unwrap_or_else(|error| panic!("{formula}: {error}"));
            let expected: Vec<(String, i32)> = expected
                .iter()
                .map(|(symbol, count)| (symbol.to_string(), *count))
                .collect();
            assert_eq!(got, expected, "composition {formula}");
        }
    }

    #[test]
    fn net_formula_unimod_isotope_prefix() {
        // UniMod labelling compositions carry stable-isotope symbols.
        let terms = split_net_formula("H(30) 2H(8) C(22) N(4) O(6) S").unwrap();
        assert!(terms.contains(&("2H".to_owned(), 8)), "{terms:?}");
        assert!(terms.contains(&("C".to_owned(), 22)));
        // `C2H2O` must NOT treat the `2` after C as an isotope prefix.
        let normal = split_net_formula("C2H2O").unwrap();
        assert_eq!(
            normal,
            vec![
                ("C".to_owned(), 2),
                ("H".to_owned(), 2),
                ("O".to_owned(), 1)
            ]
        );
        // The deuterium mass is the isotope mass, not hydrogen's.
        let m = net_formula_mass("2H(1)").unwrap();
        assert!((m - 2.014_101_777_8).abs() < 1e-9, "2H mono {m}");
    }

    #[test]
    fn net_formula_rejects_malformed() {
        assert!(split_net_formula("").is_err(), "empty");
        assert!(split_net_formula("2").is_err(), "bare digit");
        assert!(split_net_formula("+2").is_err(), "signed count first");
        assert!(split_net_formula("Br H -").is_err(), "trailing sign");
        assert!(split_net_formula("H ? 2").is_err(), "bad char");
        // A digit prefix followed by a lowercase letter is not an element.
        assert!(split_net_formula("2h").is_err(), "lowercase after digit");
    }

    #[test]
    fn isotope_label_keeps_exact_mass_and_identity() {
        // A deuterium label must be representable in the composition with its
        // isotope identity, so the mass uses 2.0141, not 1.0078. Empirically
        // (UNIMOD delta_mono == composition element-mass sum), the label's
        // delta is the isotope's own mass — it is an addition, not an H swap.
        let p = parse_proforma(
            "[D2]-PEPTIDE",
            &[ExternalModification {
                name: "D2".to_owned(),
                formula: "2H(2)".to_owned(),
            }],
        )
        .unwrap();
        let plain = parse_proforma("PEPTIDE", &[]).unwrap();
        let d2 = p
            .composition
            .iter()
            .find(|e| e.element == "H" && e.isotope == 2)
            .expect("deuterium must be its own isotope entry");
        assert_eq!(d2.count, 2);
        assert!(
            (p.monoisotopic_mass - plain.monoisotopic_mass - 2.0 * 2.014_101_777_8).abs() < 1e-9,
            "2H(2) delta is two deuterium masses (UNIMOD semantics)"
        );
        // Round-trip through the composition must preserve the isotope.
        let comp = chemical_from_counts(&p.composition);
        let has_deuterium = comp
            .iter()
            .any(|(elt, count)| elt.isotope == 2 && *count == 2);
        assert!(has_deuterium, "chemical_from_counts keeps isotope identity");
    }

    #[test]
    fn injected_unimod_composition_string_parses() {
        // A UniMod record's raw composition feeds the exact path directly.
        let p = parse_proforma(
            "[UniModStyle]-PEPTIDE",
            &[ExternalModification {
                name: "UniModStyle".to_owned(),
                formula: "S -1 Se 1".to_owned(),
            }],
        )
        .unwrap();
        let s = count_of(&p.composition, "S");
        let se = count_of(&p.composition, "Se");
        assert_eq!(
            s,
            count_of(&parse_proforma("PEPTIDE", &[]).unwrap().composition, "S") - 1
        );
        assert_eq!(se, 1, "Se must be substituted in");
        assert!((p.modifications[0].delta_mass - 47.9440).abs() < 1e-3);
    }

    // ---- position tracking + @site ----------------------------------------

    #[test]
    fn position_inferred_from_placement() {
        let p = parse_proforma("[Acetyl]-SGRGK-[Bromo]", &[]).unwrap();
        let ac = p.modifications.iter().find(|m| m.tag == "Acetyl").unwrap();
        assert_eq!(ac.position.as_deref(), Some("nTerm"));
        let br = p.modifications.iter().find(|m| m.tag == "Bromo").unwrap();
        assert_eq!(br.position.as_deref(), Some("cTerm"));
    }

    #[test]
    fn site_suffix_overrides_placement() {
        // Inline tag after residue 2 (SG), explicit @K4 and @S forms.
        let p = parse_proforma("SG[Phospho@S2]RGK[Acetyl@K]", &[]).unwrap();
        let ph = p.modifications.iter().find(|m| m.tag == "Phospho").unwrap();
        assert_eq!(ph.position.as_deref(), Some("residue:2"));
        let ac = p.modifications.iter().find(|m| m.tag == "Acetyl").unwrap();
        assert_eq!(ac.position.as_deref(), Some("site:K"));
    }

    #[test]
    fn unimod_tag_resolves_via_embedded_ontology() {
        // [UNIMOD:7] = Deamidated (O-N). The mzcore backend resolves it through
        // the embedded UniMod ontology, so this now parses instead of erroring.
        let p = parse_proforma("[UNIMOD:7]-PEPTIDE", &[]).unwrap();
        let m = p
            .modifications
            .iter()
            .find(|m| m.tag == "UNIMOD:7")
            .unwrap();
        assert!(
            (m.delta_mass - 0.984).abs() < 1e-3,
            "delta {}",
            m.delta_mass
        );
        // Deamidated is a net O-N substitution.
        let n = count_of(&p.composition, "N");
        let o = count_of(&p.composition, "O");
        let base = parse_proforma("PEPTIDE", &[]).unwrap();
        assert_eq!(n, count_of(&base.composition, "N") - 1);
        assert_eq!(o, count_of(&base.composition, "O") + 1);
    }
}
