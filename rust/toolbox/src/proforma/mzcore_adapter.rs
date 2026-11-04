//! mzcore-backed ProFormA parsing and chemistry conversion.
//!
//! This is the only module that touches `mzcore` types. The public
//! [`crate::parse_proforma`] stays f64/`ElementCount`-based; the uom `Mass`
//! types never leak past this boundary. Every modification the parser returns
//! is reduced to its `MolecularFormula` (extracted directly from the
//! `SimpleModificationInner` variants) and flattened into [`crate::ElementCount`]
//! plus a `delta_mass` — mzcore's formula mass, the single mass source. This
//! module also produces the whole-peptide `MolecularFormula` that backs
//! `ParsedProForma`'s `composition` and `monoisotopic_mass`.
//!
//! Named modifications are resolved against the embedded UniMod/PSI-MOD
//! ontologies ([`mzcore::ontology::STATIC_ONTOLOGIES`]) so `[UNIMOD:n]`, bare
//! names like `[Carboxymethyl]`, `[GNO:...]` and `[Glycan:...]` all parse
//! without an external database.

use std::num::NonZeroU16;

use mzcore::chemistry::{AmbiguousMolecule, MolecularFormula, Molecule, OutputMolecularFormula};
use mzcore::ontology::STATIC_ONTOLOGIES;
use mzcore::prelude::SequencePosition;
use mzcore::sequence::{Linked, Peptidoform, SimpleModificationInner};

use crate::proforma::preprocess::Normalized;
use crate::MassSpecError;

/// The chemistry of one parsed modification, aligned by index with the input
/// tags.
#[derive(Debug)]
pub(crate) struct TagChemistry {
    /// The exact element counts contributed by this modification (empty for a
    /// mass-only or `INFO` tag). Negative counts are allowed (substitutions).
    pub(crate) counts: Vec<crate::ElementCount>,
    /// The monoisotopic delta mass (Da), computed by mzcore from the
    /// modification's `MolecularFormula` — the single mass source. For a
    /// mass-only tag this is the bare numeric delta.
    pub(crate) delta_mass: f64,
}

/// Parse the normalized text with mzcore and reduce every modification to its
/// chemistry, plus the whole-peptide `MolecularFormula` — the single source of
/// truth for the final composition and monoisotopic mass. The caller is
/// responsible for ensuring `normalized.tags.len()` matches the number of
/// modifications mzcore reports (guarded in [`crate::parse_proforma`]).
pub(crate) fn parse_with_mzcore(
    normalized: &Normalized,
) -> Result<(Vec<TagChemistry>, MolecularFormula), MassSpecError> {
    let (peptidoform, _warnings) =
        Peptidoform::pro_forma(&normalized.mzcore_text, &STATIC_ONTOLOGIES).map_err(|errors| {
            let detail = errors
                .iter()
                .map(|e| e.to_string())
                .collect::<Vec<_>>()
                .join("; ");
            // The CLI contract (cli_contract.rs) expects failures for
            // unresolvable tags to carry this prefix; mzcore's own error
            // text is appended as the detail.
            MassSpecError::InvalidQuery(format!("unknown ProForma modification: {detail}"))
        })?;

    let mut chemistries = Vec::new();
    collect_modifications(&peptidoform, &mut chemistries)?;

    // The number of modifications must match the number of tags we classified.
    // A mismatch means our scanning and mzcore disagree about the string, which
    // would silently produce wrong masses — refuse rather than guess.
    if chemistries.len() != normalized.tags.len() {
        return Err(MassSpecError::InvalidQuery(format!(
            "ProFormA modification count mismatch: expected {}, got {}",
            normalized.tags.len(),
            chemistries.len()
        )));
    }

    // The whole-peptide formula (residues + water + all modifications) is the
    // single source of truth for the final composition and monoisotopic mass.
    // A TopPIC proteoform is always linear; a non-linear peptidoform is an
    // error the caller cannot represent in the flat `ElementCount` model anyway.
    let whole_formula: MolecularFormula = peptidoform
        .into_linear()
        .ok_or_else(|| MassSpecError::InvalidQuery("ProFormA must be a linear peptide".to_owned()))?
        .calculate_masses_inner::<OutputMolecularFormula>(SequencePosition::default(), 0, 0)
        .single()
        .ok_or_else(|| {
            MassSpecError::InvalidQuery(
                "ambiguous ProFormA composition is not representable".to_owned(),
            )
        })?;

    Ok((chemistries, whole_formula))
}

/// Collect every `SimpleModificationInner` from N-term, per-residue, and C-term
/// slots, in order.
///
/// mzcore represents a `(X)[mod]` range that touches the N-terminus by placing
/// the *same* `Ambiguous` modification in both the N-term slot and the first
/// residue's slot. Such duplicates share a `(group, id)` pair and are counted
/// once, so the tag count matches what the preprocessor scanned.
fn collect_modifications(
    peptidoform: &Peptidoform<Linked>,
    out: &mut Vec<TagChemistry>,
) -> Result<(), MassSpecError> {
    use std::collections::HashSet;
    let mut seen_ambiguous = HashSet::new();
    for modification in peptidoform.get_n_term() {
        push_modification_dedup(modification, &mut seen_ambiguous, out)?;
    }
    for element in peptidoform.sequence() {
        for modification in &element.modifications {
            push_modification_dedup(modification, &mut seen_ambiguous, out)?;
        }
    }
    for modification in peptidoform.get_c_term() {
        push_modification_dedup(modification, &mut seen_ambiguous, out)?;
    }
    Ok(())
}

/// Push a modification, skipping a duplicate `Ambiguous` already seen.
fn push_modification_dedup(
    modification: &mzcore::sequence::Modification,
    seen: &mut std::collections::HashSet<(String, usize)>,
    out: &mut Vec<TagChemistry>,
) -> Result<(), MassSpecError> {
    if let mzcore::sequence::Modification::Ambiguous { group, id, .. } = modification {
        if !seen.insert((group.clone(), *id)) {
            // Already counted at another slot (N-term / residue duplicate).
            return Ok(());
        }
    }
    push_modification(modification, out)
}

/// Flatten a `MolecularFormula` into the crate's `ElementCount` list.
fn formula_counts(formula: &MolecularFormula) -> Vec<crate::ElementCount> {
    formula
        .elements()
        .iter()
        .map(|(element, isotope, count)| crate::ElementCount {
            element: element.symbol().to_owned(),
            count: *count,
            isotope: isotope.map(NonZeroU16::get).unwrap_or(0),
        })
        .collect::<Vec<_>>()
}

/// Reduce one `Modification` to its chemistry. Cross-links are not
/// representable in our flat `ElementCount` model and are rejected explicitly
/// rather than silently dropped. `Ambiguous` modifications (which a `(X)` range
/// produces) are unwrapped to their underlying simple modification — the mass
/// and formula are the same regardless of placement.
fn push_modification(
    modification: &mzcore::sequence::Modification,
    out: &mut Vec<TagChemistry>,
) -> Result<(), MassSpecError> {
    let inner = match modification {
        mzcore::sequence::Modification::Simple(sm) => &**sm,
        mzcore::sequence::Modification::Ambiguous { modification, .. } => &**modification,
        other => {
            return Err(MassSpecError::InvalidQuery(format!(
                "unsupported ProFormA modification kind: {other:?}"
            )));
        }
    };

    // mzcore 0.2 removed the `Chemical` trait; extract the `MolecularFormula`
    // directly from the `SimpleModificationInner` variant.
    let formula = match inner {
        SimpleModificationInner::Formula(f) => f,
        SimpleModificationInner::Database { formula, .. } => formula,
        SimpleModificationInner::Mass(_, mass, _) => {
            // A bare numeric delta `[+42.0106]` carries only a mass shift.
            // OrderedMass -> Mass (uom Quantity); read its `.value` field (Da).
            let delta_mass = mass.into_inner().value;
            out.push(TagChemistry {
                counts: Vec::new(),
                delta_mass,
            });
            return Ok(());
        }
        SimpleModificationInner::Info(_) => {
            out.push(TagChemistry {
                counts: Vec::new(),
                delta_mass: 0.0,
            });
            return Ok(());
        }
        SimpleModificationInner::Glycan(monosaccharides) => {
            // Sum each monosaccharide × count into one formula via the
            // `Molecule` trait (OutputMolecularFormula mode).
            let mut formula = MolecularFormula::default();
            for (saccharide, count) in monosaccharides {
                let f: MolecularFormula = saccharide.calculate_mass::<OutputMolecularFormula>();
                formula += &(f * count);
            }
            let counts = formula_counts(&formula);
            out.push(TagChemistry {
                counts,
                delta_mass: formula.monoisotopic_mass().value,
            });
            return Ok(());
        }
        other => {
            return Err(MassSpecError::InvalidQuery(format!(
                "unsupported ProFormA modification kind: {other:?}"
            )));
        }
    };

    let counts = formula_counts(formula);

    out.push(TagChemistry {
        counts,
        // mzcore's molecular-formula mass is the single mass source; it already
        // includes the formula's `additional_mass`.
        delta_mass: formula.monoisotopic_mass().value,
    });
    Ok(())
}
