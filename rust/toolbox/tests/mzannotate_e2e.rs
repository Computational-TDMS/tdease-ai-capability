//! Golden regression: mzannotate fragments a REAL TopPIC proteoform through the
//! public API, and the key m/z is exact.
//!
//! This pins two facts that matter:
//! 1. mzcore 0.2.0's `MassOutput::contains_negative_amount` is inverted
//!    (`mass.value >= 0.0`), so the `OutputMonoIsotopic` mode `continue`s every
//!    normal positive-mass fragment and yields an empty series. We therefore
//!    fragment in `OutputMolecularFormula` mode, which implements the correct
//!    element-count semantics. If someone flips this back to the mono mode, this
//!    test turns red instead of silently returning an empty fragment set.
//! 2. mzannotate's fragmentation of a real (modified) TopPIC proteoform gives
//!    the expected c1 m/z: Acetyl + Glu + NH3 + proton = 42.010565 +
//!    129.042593 + 17.026549 + 1.007276 = 189.086983.
use mzannotate::fragment::FragmentType;
use mzannotate::prelude::*;
use mzcore::chemistry::{MassMode, OutputMolecularFormula};
use mzcore::prelude::*;
use mzcore::system::e;
use mzcore::system::isize::Charge;

#[test]
fn real_toppic_proteoform_fragments_golden_mz() {
    // A real TopPIC proteoform, N-term acetylated, in the standard ProFormA
    // form our preprocess layer emits (TopPIC `M.(S)[Acetyl]...K.K` → this).
    let proforma = "[Formula:C2H2O]-ESSSKSSQPLASKQEKDGTEKRGRGRPRKQPPK";
    let (peptidoform, warnings) =
        Peptidoform::pro_forma(proforma, &mzcore::ontology::STATIC_ONTOLOGIES).expect("parses");
    assert!(warnings.is_empty(), "no parse warnings: {warnings:?}");
    let peptide = peptidoform.into_linear().expect("linear");

    let fragments: Vec<Fragment<OutputMolecularFormula>> = peptide
        .generate_theoretical_fragments::<OutputMolecularFormula>(
            Charge::new::<e>(1),
            FragmentationModel::all(),
        );
    assert!(!fragments.is_empty(), "mzannotate produced zero fragments");

    // Find c1: ion type `c`, series number 1, charge 1.
    let c1 = fragments
        .iter()
        .find(|f| {
            matches!(&f.ion, FragmentType::c(pos, 0) if pos.series_number == 1)
                && f.charge.value == 1
        })
        .unwrap_or_else(|| panic!("no c1 among {} fragments", fragments.len()));
    let c1_mz = c1.mz(MassMode::Monoisotopic).expect("c1 has a mass").value;
    let expected = 189.086_983;
    assert!(
        (c1_mz - expected).abs() < 1e-3,
        "c1 m/z {c1_mz} != {expected}"
    );
}
