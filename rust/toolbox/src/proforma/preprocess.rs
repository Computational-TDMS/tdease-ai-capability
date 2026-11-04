//! TopPIC → ProFormA preprocessing layer.
//!
//! TopPIC proteoform strings use a flanking notation that ProFormA does not:
//! `M.(S)[Acetyl]ESSSKSSQPLASKQEKDGTEKRGRGRPRKQPPK.K` — a single flanking
//! residue outside a leading/trailing `.`, with the N-terminal leader residue
//! (M) and the C-terminal trailer residue (K) separated by dots. ProFormA has
//! no dot token at all (`parse_linear_peptide` errors "not a valid amino
//! acid"), so this layer strips the flanks and, optionally, maps each bracketed
//! modification to a form the parser accepts:
//!
//! - `[+15.9949]` / `[-18.0106]` — numeric delta, pass through (mass-only).
//! - `[Formula:...]` — pass through after normalization via
//!   [`crate::split_net_formula`] into a compact `C2H2O`-style string the
//!   downstream parser accepts.
//! - a name that matches an `external_mods` entry — rewritten to
//!   `[Formula:<formula>]` (explicit caller override of ontology chemistry).
//! - `UNIMOD:` / `GNO:` / `Glycan:` / `glycanstructure:` — passed through
//!   verbatim for the downstream parser's embedded ontologies.
//! - any other bare name — passed through; the downstream parser resolves it
//!   against UniMod/PSI-MOD. If it resolves to nothing, the error surfaces with
//!   the "unknown ProFormA modification" prefix the CLI contract requires.
//!
//! `@site` suffixes (`[Acetyl@K4]`) are stripped here and the position is
//! computed by placement; the downstream parser never sees an `@`.
//!
//! A `(X)` range before a tag is ProFormA 2.1's `modRange` — the range's
//! residues are recorded and the *single* following tag's position is inferred
//! from them (`site:AA`). The range produces exactly one modification:
//! downstream mass reconstruction sums `modifications[].delta_mass`, so
//! expanding per residue would over-count.

use crate::proforma::{parse_site, ExternalModification};
use crate::MassSpecError;

/// One bracket tag after scanning, before the downstream parse.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RawTag {
    /// The bracket body without any `@site` suffix, e.g. `Acetyl`, `UNIMOD:7`.
    pub(crate) tag: String,
    /// The explicit `@site` suffix, if any.
    pub(crate) site: Option<String>,
    /// The inferred position (`nTerm`/`cTerm`/`residue:N`/`site:AA`/`unknown`).
    pub(crate) position: Option<String>,
    /// Exact TopPIC localization interval, when the tag follows `(SUBSEQ)`.
    pub(crate) range: Option<(usize, usize)>,
    /// The rewritten bracket body fed to the downstream parser.
    pub(crate) rewritten: String,
}

/// The normalized, downstream-safe form of the input.
#[derive(Debug)]
pub(crate) struct Normalized {
    /// The residue sequence without flanks or brackets.
    pub(crate) sequence: String,
    /// Leading flank residue, if any (e.g. `M` in `M.PEPTIDE`). Recorded for
    /// completeness; flanking residues carry no mass and are not currently
    /// consumed downstream.
    #[allow(dead_code)]
    pub(crate) n_flank: Option<char>,
    /// Trailing flank residue, if any (e.g. `K` in `PEPTIDE.K`). Same note as
    /// `n_flank`.
    #[allow(dead_code)]
    pub(crate) c_flank: Option<char>,
    /// The rewritten string ready for the downstream parser.
    pub(crate) mzcore_text: String,
    /// Tags, in order of appearance, aligned with downstream modifications.
    pub(crate) tags: Vec<RawTag>,
}

/// Replace each multi-residue modRange with the two concrete endpoint
/// placements that span every possible fragment-mass state. A cleavage inside
/// the range can then be identified as ambiguous and excluded; cleavages
/// outside the range collapse to one deterministic mzcore result.
pub(crate) fn concretize_mod_ranges(text: &str) -> Result<Vec<String>, MassSpecError> {
    let mut variants = vec![text.to_owned()];
    loop {
        let mut expanded = false;
        let mut next = Vec::new();
        for variant in variants {
            let Some((start, close, tag_end)) = first_mod_range(&variant)? else {
                next.push(variant);
                continue;
            };
            expanded = true;
            let residues = &variant[start + 1..close];
            let tag = &variant[close + 1..=tag_end];
            let mut chars = residues.chars();
            let first = chars.next().ok_or_else(|| {
                MassSpecError::InvalidQuery("empty ProForma localization range".to_owned())
            })?;
            let last = residues.chars().last().unwrap_or(first);
            let first_len = first.len_utf8();
            let last_start = residues.len() - last.len_utf8();
            let prefix = &variant[..start];
            let suffix = &variant[tag_end + 1..];
            next.push(format!(
                "{prefix}{first}{tag}{}{suffix}",
                &residues[first_len..]
            ));
            next.push(format!(
                "{prefix}{}{last}{tag}{suffix}",
                &residues[..last_start]
            ));
        }
        next.sort();
        next.dedup();
        variants = next;
        if !expanded {
            return Ok(variants);
        }
    }
}

fn first_mod_range(text: &str) -> Result<Option<(usize, usize, usize)>, MassSpecError> {
    let bytes = text.as_bytes();
    let Some(start) = bytes.iter().position(|byte| *byte == b'(') else {
        return Ok(None);
    };
    let close = bytes[start + 1..]
        .iter()
        .position(|byte| *byte == b')')
        .map(|offset| start + 1 + offset)
        .ok_or_else(|| MassSpecError::InvalidQuery("unclosed ProForma modRange".to_owned()))?;
    if bytes.get(close + 1) != Some(&b'[') {
        return Err(MassSpecError::InvalidQuery(
            "ProForma modRange must be followed by a modification".to_owned(),
        ));
    }
    let mut depth = 0usize;
    for (index, byte) in bytes.iter().enumerate().skip(close + 1) {
        match byte {
            b'[' => depth += 1,
            b']' => {
                depth -= 1;
                if depth == 0 {
                    return Ok(Some((start, close, index)));
                }
            }
            _ => {}
        }
    }
    Err(MassSpecError::InvalidQuery(
        "unclosed modification after ProForma modRange".to_owned(),
    ))
}

/// Strip TopPIC flanking dots and rewrite tags, producing a downstream-safe
/// string plus per-tag metadata.
pub(crate) fn normalize(
    input: &str,
    external_mods: &[ExternalModification],
) -> Result<Normalized, MassSpecError> {
    // ── Flank stripping ─────────────────────────────────────────────────────
    let mut body = input;
    let mut n_flank = None;
    let mut c_flank = None;

    let bytes = body.as_bytes();
    // Leading flank: `M.` — an uppercase letter then a dot.
    if bytes.len() >= 2 && bytes[0].is_ascii_uppercase() && bytes[1] == b'.' {
        n_flank = Some(bytes[0] as char);
        body = &body[2..];
    }
    // Trailing flank: `.K` — a dot then an uppercase letter at the very end.
    let b = body.as_bytes();
    if b.len() >= 2 && b[b.len() - 2] == b'.' && b[b.len() - 1].is_ascii_uppercase() {
        c_flank = Some(b[b.len() - 1] as char);
        body = &body[..b.len() - 2];
    }

    if body.is_empty() {
        return Err(MassSpecError::InvalidQuery(
            "ProFormA contains no amino-acid sequence".to_owned(),
        ));
    }

    // ── Tag scan ────────────────────────────────────────────────────────────
    let mut out = String::new();
    let mut sequence = String::new();
    let mut tags = Vec::new();
    let mut i = 0usize;
    let chars: Vec<char> = body.chars().collect();
    // Position inference state.
    let mut placement = Place::Before;
    let mut residue_index = 0usize;
    // A pending `(X)` range: its residue letters (for site: inference).
    let mut pending_range: Option<(String, usize, usize)> = None;

    while i < chars.len() {
        let c = chars[i];
        match c {
            '[' => {
                let end = chars[i + 1..]
                    .iter()
                    .position(|ch| *ch == ']')
                    .map(|off| i + 1 + off)
                    .ok_or_else(|| {
                        MassSpecError::InvalidQuery("unclosed ProFormA tag".to_owned())
                    })?;
                let raw: String = chars[i + 1..end].iter().collect();
                let tag = classify_tag(&raw, external_mods)?;
                // Position: explicit @site wins; else pending range; else
                // placement.
                let position = tag
                    .site
                    .as_deref()
                    .map(parse_site)
                    .or_else(|| {
                        pending_range
                            .as_ref()
                            .map(|(letters, _, _)| format!("site:{}", letters.to_uppercase()))
                    })
                    .or_else(|| match placement {
                        Place::Before => Some("nTerm".to_owned()),
                        Place::AfterResidue(n) => Some(format!("residue:{n}")),
                        Place::AfterDash => Some("cTerm".to_owned()),
                    });
                out.push_str(&format!("[{}]", tag.rewritten));
                tags.push(RawTag {
                    position,
                    range: pending_range.take().map(|(_, start, end)| (start, end)),
                    site: tag.site,
                    tag: tag.tag,
                    rewritten: tag.rewritten,
                });
                i = end + 1;
            }
            '(' => {
                // A `(X)` range — collect its letters for site inference. The
                // downstream parser handles the range itself; it stays in `out`
                // verbatim so the modRange semantics are preserved.
                let close = chars[i + 1..]
                    .iter()
                    .position(|ch| *ch == ')')
                    .map(|off| i + 1 + off)
                    .ok_or_else(|| {
                        MassSpecError::InvalidQuery("unclosed '(' range in ProFormA".to_owned())
                    })?;
                let inner: String = chars[i + 1..close].iter().collect();
                let letters: String = inner
                    .chars()
                    .filter(|ch| ch.is_ascii_alphabetic())
                    .collect();
                let start = residue_index + 1;
                let end = residue_index + letters.len();
                pending_range = if letters.is_empty() {
                    None
                } else {
                    Some((letters.clone(), start, end))
                };
                // A one-residue TopPIC range is exact localization, not
                // ambiguity. Canonicalize it to ordinary residue placement;
                // mzcore 0.2 otherwise duplicates an N-terminal ambiguous mod.
                if letters.len() == 1 {
                    out.push_str(&letters);
                } else {
                    out.push('(');
                    out.push_str(&inner);
                    out.push(')');
                }
                // Range residues are not part of the top-level sequence string
                // but ARE part of the chemical mass (mzcore counts them). Track
                // them in `sequence` for composition reconstruction below.
                for ch in inner.chars() {
                    if ch.is_ascii_alphabetic() {
                        sequence.push(ch);
                        residue_index += 1;
                    }
                }
                i = close + 1;
                continue;
            }
            '-' => {
                placement = Place::AfterDash;
                out.push(c);
                i += 1;
            }
            ch if ch.is_ascii_uppercase() => {
                // A residue. Validate it.
                if crate::proforma::residue_mass(ch).is_none() {
                    return Err(MassSpecError::InvalidQuery(format!(
                        "unknown ProForma residue {ch}"
                    )));
                }
                out.push(ch);
                sequence.push(ch);
                residue_index += 1;
                placement = Place::AfterResidue(residue_index);
                i += 1;
            }
            ch if ch.is_ascii_whitespace() => {
                out.push(ch);
                i += 1;
            }
            '.' => {
                // A dot in the middle (not a flank) is a downstream error.
                return Err(MassSpecError::InvalidQuery(format!(
                    "unsupported ProForma token '.' at position {i}"
                )));
            }
            other => {
                return Err(MassSpecError::InvalidQuery(format!(
                    "unsupported ProForma token {other}"
                )));
            }
        }
    }

    // ── N-terminal `-` normalization ────────────────────────────────────────
    // ProFormA requires a `-` between an N-terminal modification and the
    // sequence (`[Acetyl]-PEPTIDE`). If the output starts with an N-terminal
    // `[` tag but the first *outer* `]` is not followed by a `-`, insert one so
    // the downstream parser treats it as N-term. This matches the historical
    // parser's acceptance of `[Acetyl]PEPTIDE`. Nested brackets (an isotope
    // formula like `[Formula:[2H2]]`) are matched by depth so their inner `]`
    // does not fool the check.
    let mut needs_dash = false;
    if out.starts_with('[') {
        needs_dash = true;
        // Find the closing `]` of the first (outer) tag, skipping nested `[...]`.
        let mut depth = 0i32;
        let chars = out.char_indices();
        for (idx, ch) in chars {
            match ch {
                '[' => depth += 1,
                ']' => {
                    depth -= 1;
                    if depth == 0 {
                        let after = out[idx + 1..].trim_start();
                        if after.starts_with('-') {
                            needs_dash = false;
                        }
                        break;
                    }
                }
                _ => {}
            }
        }
    }
    if needs_dash {
        out.insert(0, '-');
    }

    Ok(Normalized {
        sequence,
        n_flank,
        c_flank,
        mzcore_text: out,
        tags,
    })
}

/// Where the scanner currently sits, for position inference.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Place {
    Before,
    AfterResidue(usize),
    AfterDash,
}

struct Classified {
    tag: String,
    site: Option<String>,
    rewritten: String,
}

impl Classified {
    fn new(tag: &str, site: Option<String>, rewritten: String) -> Self {
        Self {
            tag: tag.to_owned(),
            site,
            rewritten,
        }
    }
}

/// Classify a single bracketed body (no brackets) into a rewritten form.
fn classify_tag(
    raw: &str,
    external_mods: &[ExternalModification],
) -> Result<Classified, MassSpecError> {
    // Split an optional `@site` suffix.
    let (body, site) = match raw.find('@') {
        Some(pos) => (&raw[..pos], Some(raw[pos + 1..].to_owned())),
        None => (raw, None),
    };

    // 1. Numeric delta.
    if body.parse::<f64>().is_ok() {
        let rewritten = if body.starts_with('+') || body.starts_with('-') {
            body.to_owned()
        } else {
            format!("+{body}")
        };
        return Ok(Classified::new(body, site, rewritten));
    }

    // 2. `Formula:` prefix — normalize the net formula.
    if let Some(formula) = body.strip_prefix("Formula:") {
        let normalized = normalize_formula(formula)?;
        return Ok(Classified::new(body, site, format!("Formula:{normalized}")));
    }

    // Explicit caller-supplied formulas override ontology lookup.
    if let Some(mod_entry) = external_mods
        .iter()
        .find(|m| m.name.eq_ignore_ascii_case(body))
    {
        let normalized = normalize_formula(&mod_entry.formula)?;
        return Ok(Classified::new(body, site, format!("Formula:{normalized}")));
    }

    // Names and ontology tags pass through; the downstream parser resolves them
    //    against UniMod/PSI-MOD, or errors with the required prefix.
    Ok(Classified::new(body, site, body.to_owned()))
}

/// Normalize a net formula (`C2H2O`, `Br-H`, `Br H -1`, `O-N`, `H(-1) Br`) into
/// a compact downstream-safe form. Uses `split_net_formula` (already accepting
/// all real-world spellings) then re-emits a canonical `C2H2O`-style string.
///
/// Stable-isotope symbols (`2H`, `13C`) are rendered in the downstream parser's
/// bracket spelling `[2H2]` (mzcore reads isotopes as `[<mass><element><count>]`,
/// e.g. `[13C2]`), so an injected UniMod labelling formula round-trips exactly.
pub(crate) fn normalize_formula(formula: &str) -> Result<String, MassSpecError> {
    let terms = crate::split_net_formula(formula)?;
    let mut out = String::new();
    for (sym, count) in terms {
        let is_isotope = sym.chars().next().is_some_and(|c| c.is_ascii_digit());
        if is_isotope {
            // Bracket spelling: [2H2], [13C2].
            out.push('[');
            out.push_str(&sym);
            out.push_str(&count.to_string());
            out.push(']');
        } else {
            match count {
                1 => out.push_str(&sym),
                n => {
                    out.push_str(&sym);
                    out.push_str(&n.to_string());
                }
            }
        }
    }
    if out.is_empty() {
        return Err(MassSpecError::InvalidQuery(format!(
            "empty net formula [{formula}]"
        )));
    }
    Ok(out)
}
