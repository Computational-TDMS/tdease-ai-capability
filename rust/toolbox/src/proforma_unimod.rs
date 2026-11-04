//! UniMod-aware ProFormA resolution (D12-R completion). The name resolver that
//! used to live in TS (`proforma-mass.ts` + a minimal SQLite index) now runs in
//! the Rust kernel: given the in-memory UniMod records, names and amino-acid
//! table (loaded from the CLI's `--db`), it resolves `[UNIMOD:n]` and exact
//! names to records, computes the average mass from reference data, rewrites
//! resolved tags to `[M<id>]`, and hands the chemistry to `parse_proforma` —
//! the single exact-composition kernel. TS reaches this through the CLI, so it
//! holds no UniMod SQLite reader at all.

use std::collections::HashMap;

use serde::Serialize;

use crate::proforma::{parse_proforma, ExternalModification};
use crate::unimod::{display_id, UniModAminoAcid, UniModNameEntry, UniModRecord, USER_OFFSET};
use crate::MassSpecError;

const WATER_AVERAGE: f64 = 18.01528;

/// One entry in the modification trace, serialized with only the keys its
/// variant carries (resolved record / pass-through / numeric delta).
#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ResolvedModification {
    pub location: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub unimod_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub title: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub delta_mono: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub delta_average: Option<f64>,
    /// The record source: `"unimod"` or `"user"` (D17). Absent for
    /// pass-through/numeric tokens that resolve no record.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub custom: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub token: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mass_by: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub average_mass_available: Option<bool>,
}

/// The resolved ProFormA: exact sequence/mass from `parse_proforma`, plus the
/// reference-data average mass and the UniMod modification trace.
#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ProFormaResolution {
    pub proforma: String,
    pub sequence: String,
    pub monoisotopic_mass: f64,
    /// Serialized as `null` when a numeric/Formula/unknown mod carries no
    /// average-mass definition (the TS port consumed a `null`, not an omission).
    pub average_mass: Option<f64>,
    pub modifications: Vec<ResolvedModification>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub average_mass_unavailable_reason: Option<String>,
}

/// The UniMod reference data the resolver reads (records, names, amino acids).
/// Loaded once by the caller; the resolver performs no I/O.
pub struct UniModProFormaResolver<'a> {
    record_by_id: HashMap<u32, &'a UniModRecord>,
    name_to_records: HashMap<String, Vec<u32>>,
    amino_acids: HashMap<String, &'a UniModAminoAcid>,
}

impl<'a> UniModProFormaResolver<'a> {
    pub fn new(
        records: &'a [UniModRecord],
        names: &'a [UniModNameEntry],
        amino_acids: &'a [UniModAminoAcid],
    ) -> Self {
        let mut record_by_id = HashMap::new();
        for record in records {
            record_by_id.insert(record.record_id, record);
        }
        let mut name_to_records: HashMap<String, Vec<u32>> = HashMap::new();
        for entry in names {
            name_to_records
                .entry(entry.name.to_lowercase())
                .or_default()
                .push(entry.modification_id);
        }
        // Deterministic ordering within a name → the lowest record id wins when
        // a name matches several records (the TS sort kept ids ascending).
        for ids in name_to_records.values_mut() {
            ids.sort_unstable();
            ids.dedup();
        }
        let mut amino = HashMap::new();
        for aa in amino_acids {
            amino.insert(aa.title.to_uppercase(), aa);
        }
        Self {
            record_by_id,
            name_to_records,
            amino_acids: amino,
        }
    }

    fn get_record(&self, id: u32) -> Option<&'a UniModRecord> {
        self.record_by_id.get(&id).copied()
    }

    fn exact_names(&self, name: &str) -> Vec<&'a UniModRecord> {
        let Some(ids) = self.name_to_records.get(&name.to_lowercase()) else {
            return Vec::new();
        };
        ids.iter()
            .filter_map(|id| self.record_by_id.get(id).copied())
            .collect()
    }

    /// Resolve one ProFormA against the reference data and hand chemistry to
    /// the exact-composition kernel.
    pub fn resolve(&self, proforma: &str) -> Result<ProFormaResolution, MassSpecError> {
        let text = proforma.trim();
        if text.is_empty() {
            return Err(MassSpecError::InvalidQuery(
                "PROFORMA_EMPTY: proforma must contain a sequence".to_owned(),
            ));
        }
        let mut modifications: Vec<ResolvedModification> = Vec::new();
        let mut external_mods: Vec<ExternalModification> = Vec::new();
        let mut sequence = String::new();
        let mut rewritten = String::new();
        let mut average: Option<f64> = Some(WATER_AVERAGE);
        let mut cursor = 0usize;

        if text.starts_with('[') {
            self.consume_bracket(
                text,
                &mut cursor,
                "nTerm",
                &mut rewritten,
                &mut average,
                &mut modifications,
                &mut external_mods,
            )?;
            if text[cursor..].starts_with('-') {
                cursor += 1;
                rewritten.push('-');
            } else {
                return Err(MassSpecError::InvalidQuery(
                    "PROFORMA_SYNTAX: N-terminal modification must be followed by -".to_owned(),
                ));
            }
        }

        let chars: Vec<char> = text.chars().collect();
        while cursor < chars.len() {
            let token = chars[cursor];
            if token == '-' {
                cursor += 1;
                self.consume_bracket(
                    text,
                    &mut cursor,
                    "cTerm",
                    &mut rewritten,
                    &mut average,
                    &mut modifications,
                    &mut external_mods,
                )?;
                if cursor != chars.len() {
                    return Err(MassSpecError::InvalidQuery(
                        "PROFORMA_SYNTAX: C-terminal modification must end the ProFormA string"
                            .to_owned(),
                    ));
                }
                break;
            }
            if !token.is_ascii_uppercase() {
                return Err(MassSpecError::InvalidQuery(format!(
                    "PROFORMA_SEQUENCE_INVALID: unexpected token {token:?}"
                )));
            }
            let residue = self
                .amino_acids
                .get(&token.to_string())
                .copied()
                .ok_or_else(|| {
                    MassSpecError::InvalidQuery(format!("PROFORMA_RESIDUE_UNKNOWN: {token}"))
                })?;
            sequence.push(token);
            rewritten.push(token);
            if let Some(a) = average.as_mut() {
                *a += residue.average_mass;
            }
            cursor += 1;
            while cursor < chars.len() && chars[cursor] == '[' {
                let location = format!("residue:{}", sequence.len());
                self.consume_bracket(
                    text,
                    &mut cursor,
                    &location,
                    &mut rewritten,
                    &mut average,
                    &mut modifications,
                    &mut external_mods,
                )?;
            }
        }

        if sequence.is_empty() {
            return Err(MassSpecError::InvalidQuery(
                "PROFORMA_EMPTY_SEQUENCE: proforma must contain at least one residue".to_owned(),
            ));
        }

        let has_numeric = modifications
            .iter()
            .any(|m| m.delta_mono.is_some() && m.average_mass_available == Some(false));
        let average_mass_unavailable_reason = if average.is_none() {
            Some(if has_numeric {
                "A signed numeric modification delta has no elemental composition or average-mass definition.".to_owned()
            } else {
                "A modification without an elemental composition (Formula: or unknown name) has no average-mass definition.".to_owned()
            })
        } else {
            None
        };

        let parsed = parse_proforma(&rewritten, &external_mods)?;
        Ok(ProFormaResolution {
            proforma: text.to_owned(),
            sequence: parsed.sequence,
            monoisotopic_mass: parsed.monoisotopic_mass,
            average_mass: average,
            modifications,
            average_mass_unavailable_reason,
        })
    }
}

/// The TS regex was `/^[+-](?:\d+(?:\.\d*)?|\.\d+)$/`: an optional sign, then
/// either digits with an optional (possibly empty) fraction, or a leading dot
/// with at least one digit.
fn is_numeric(token: &str) -> bool {
    let rest = token.strip_prefix(['+', '-']).unwrap_or(token);
    if rest.is_empty() {
        return false;
    }
    let mut chars = rest.chars().peekable();
    let mut seen_digit = false;
    let mut seen_dot = false;
    if chars.peek() == Some(&'.') {
        chars.next();
        seen_dot = true;
        if chars.peek().is_none_or(|c| !c.is_ascii_digit()) {
            return false;
        }
    }
    for c in chars {
        if c.is_ascii_digit() {
            seen_digit = true;
        } else if c == '.' && !seen_dot {
            seen_dot = true;
        } else {
            return false;
        }
    }
    seen_digit
}

impl<'a> UniModProFormaResolver<'a> {
    /// Append a resolved UniMod record's contribution: average mass, provenance
    /// trace, and — when the record carries a composition — an external mod so
    /// the exact-composition kernel resolves it; mass-only records rewrite to a
    /// numeric delta.
    fn add_resolved(
        &self,
        record: &UniModRecord,
        location: &str,
        rewritten: &mut String,
        average: &mut Option<f64>,
        modifications: &mut Vec<ResolvedModification>,
        external_mods: &mut Vec<ExternalModification>,
    ) {
        if let Some(a) = average.as_mut() {
            *a += record.delta_average;
        }
        modifications.push(ResolvedModification {
            location: location.to_owned(),
            unimod_id: Some(display_id(record.record_id, &record.source)),
            title: Some(record.title.clone()),
            delta_mono: Some(record.delta_mono),
            delta_average: Some(record.delta_average),
            source: Some(record.source.clone()),
            custom: Some(record.source == "user"),
            token: None,
            mass_by: None,
            average_mass_available: None,
        });
        if !record.composition.is_empty() {
            let name = format!("M{}", record.record_id);
            external_mods.push(ExternalModification {
                name: name.clone(),
                formula: record.composition.clone(),
            });
            rewritten.push_str(&format!("[{name}]"));
        } else {
            *average = None;
            rewritten.push_str(&format!("[+{}]", record.delta_mono));
        }
    }

    /// Consume one `[...]` bracket in place, resolving its tag (UNIMOD id,
    /// exact name, numeric delta, or Formula:/unknown pass-through).
    fn consume_bracket(
        &self,
        text: &str,
        cursor: &mut usize,
        location: &str,
        rewritten: &mut String,
        average: &mut Option<f64>,
        modifications: &mut Vec<ResolvedModification>,
        external_mods: &mut Vec<ExternalModification>,
    ) -> Result<(), MassSpecError> {
        let open = *cursor + 1;
        let Some(close) = text[*cursor + 1..].find(']').map(|i| *cursor + 1 + i) else {
            return Err(MassSpecError::InvalidQuery(
                "PROFORMA_SYNTAX: unclosed modification bracket".to_owned(),
            ));
        };
        let token = text[open..close].trim();
        if token.is_empty() {
            return Err(MassSpecError::InvalidQuery(
                "PROFORMA_SYNTAX: empty modification bracket".to_owned(),
            ));
        }
        // Explicit `UNIMOD:<n>` / `USER:<n>` ids resolve directly (D17). A
        // `USER:n` maps to the offset user record and is verified to be one.
        if let Some(id) = token
            .strip_prefix("UNIMOD:")
            .or_else(|| token.strip_prefix("unimod:"))
            .and_then(|rest| rest.parse::<u32>().ok())
        {
            let record = self
                .get_record(id)
                .filter(|r| r.source == "unimod")
                .ok_or_else(|| {
                    MassSpecError::InvalidQuery(format!("PROFORMA_UNIMOD_UNKNOWN: {token}"))
                })?;
            self.add_resolved(
                record,
                location,
                rewritten,
                average,
                modifications,
                external_mods,
            );
        } else if let Some(id) = token
            .strip_prefix("USER:")
            .and_then(|rest| rest.parse::<u32>().ok())
        {
            let record = self
                .get_record(USER_OFFSET + id)
                .filter(|r| r.source == "user")
                .ok_or_else(|| {
                    MassSpecError::InvalidQuery(format!("PROFORMA_USER_UNKNOWN: {token}"))
                })?;
            self.add_resolved(
                record,
                location,
                rewritten,
                average,
                modifications,
                external_mods,
            );
        } else if is_numeric(token) {
            *average = None;
            modifications.push(ResolvedModification {
                location: location.to_owned(),
                unimod_id: None,
                title: None,
                delta_mono: token.parse::<f64>().ok(),
                delta_average: None,
                source: None,
                custom: None,
                token: None,
                mass_by: None,
                average_mass_available: Some(false),
            });
            rewritten.push_str(&format!("[{token}]"));
        } else if token.to_ascii_lowercase().starts_with("formula:") {
            *average = None;
            modifications.push(ResolvedModification {
                location: location.to_owned(),
                unimod_id: None,
                title: None,
                delta_mono: None,
                delta_average: None,
                source: None,
                custom: None,
                token: Some(token.to_owned()),
                mass_by: Some("core".to_owned()),
                average_mass_available: None,
            });
            rewritten.push_str(&format!("[{token}]"));
        } else {
            let named = self.exact_names(token);
            match named.len() {
                1 => {
                    self.add_resolved(
                        named[0],
                        location,
                        rewritten,
                        average,
                        modifications,
                        external_mods,
                    );
                }
                n if n > 1 => {
                    let list = named
                        .iter()
                        .map(|r| format!("{} ({})", display_id(r.record_id, &r.source), r.title))
                        .collect::<Vec<_>>()
                        .join(", ");
                    return Err(MassSpecError::InvalidQuery(format!(
                        "PROFORMA_MODIFICATION_AMBIGUOUS: {token} matches {list}. Use an explicit UNIMOD:<id>."
                    )));
                }
                _ => {
                    // Unknown to UniMod: maybe a Rust built-in name (e.g.
                    // PyroGlu) or an error — let the core decide.
                    *average = None;
                    modifications.push(ResolvedModification {
                        location: location.to_owned(),
                        unimod_id: None,
                        title: None,
                        delta_mono: None,
                        delta_average: None,
                        source: None,
                        custom: None,
                        token: Some(token.to_owned()),
                        mass_by: Some("core".to_owned()),
                        average_mass_available: None,
                    });
                    rewritten.push_str(&format!("[{token}]"));
                }
            }
        }
        *cursor = close + 1;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::unimod::{UniModAminoAcid, UniModNameEntry, UniModRecord};

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
        let mut out = Vec::new();
        for record in records() {
            out.push(UniModNameEntry {
                modification_id: record.record_id,
                name: record.title.clone(),
                kind: crate::unimod::NameKind::Title,
            });
        }
        out
    }

    fn amino_acids() -> Vec<UniModAminoAcid> {
        vec![
            UniModAminoAcid {
                title: "A".to_owned(),
                three_letter: "Ala".to_owned(),
                full_name: "Alanine".to_owned(),
                mono_mass: 71.037114,
                average_mass: 71.0779,
            },
            UniModAminoAcid {
                title: "C".to_owned(),
                three_letter: "Cys".to_owned(),
                full_name: "Cysteine".to_owned(),
                mono_mass: 103.009185,
                average_mass: 103.1388,
            },
        ]
    }

    #[test]
    fn resolves_acetyl_name_to_record_with_average_mass() {
        let (recs, names, aas) = (records(), names(), amino_acids());
        let resolver = UniModProFormaResolver::new(&recs, &names, &aas);
        let result = resolver.resolve("AC[Acetyl]").unwrap();
        assert_eq!(result.sequence, "AC");
        assert!(result.monoisotopic_mass > 0.0);
        assert!(result.average_mass.is_some());
        let first = &result.modifications[0];
        assert_eq!(first.title.as_deref(), Some("Acetyl"));
        assert_eq!(first.unimod_id.as_deref(), Some("UNIMOD:1"));
        assert_eq!(first.location, "residue:2");
    }

    #[test]
    fn bromo_resolves_to_same_exact_mass_as_formula() {
        let (recs, names, aas) = (records(), names(), amino_acids());
        let resolver = UniModProFormaResolver::new(&recs, &names, &aas);
        let named = resolver.resolve("[Bromo]-AC").unwrap();
        let formula = resolver.resolve("[Formula:Br-H]-AC").unwrap();
        assert!((named.monoisotopic_mass - formula.monoisotopic_mass).abs() < 1e-9);
        assert!(formula.average_mass.is_none());
        assert!(formula
            .average_mass_unavailable_reason
            .as_deref()
            .is_some_and(|r| r.contains("no average-mass definition")));
        assert_eq!(
            formula.modifications[0].token.as_deref(),
            Some("Formula:Br-H")
        );
    }

    #[test]
    fn numeric_delta_has_no_average_mass() {
        let (recs, names, aas) = (records(), names(), amino_acids());
        let resolver = UniModProFormaResolver::new(&recs, &names, &aas);
        let result = resolver.resolve("AC[+15.9949]").unwrap();
        assert!(result.average_mass.is_none());
        assert!(result
            .average_mass_unavailable_reason
            .as_deref()
            .is_some_and(|r| r.contains("numeric modification delta")));
        assert_eq!(result.modifications[0].delta_mono, Some(15.9949));
        assert_eq!(result.modifications[0].average_mass_available, Some(false));
    }

    #[test]
    fn unknown_name_passes_through_to_the_core() {
        let (recs, names, aas) = (records(), names(), amino_acids());
        let resolver = UniModProFormaResolver::new(&recs, &names, &aas);
        let result = resolver.resolve("AC[Carboxymethyl]").unwrap();
        assert_eq!(result.modifications[0].mass_by.as_deref(), Some("core"));
        assert!(result.average_mass.is_none());
    }

    #[test]
    fn unknown_unimod_id_is_an_honest_error() {
        let (recs, names, aas) = (records(), names(), amino_acids());
        let resolver = UniModProFormaResolver::new(&recs, &names, &aas);
        let err = resolver.resolve("AC[UNIMOD:9999]").unwrap_err();
        assert!(err.to_string().contains("PROFORMA_UNIMOD_UNKNOWN"));
    }
}
