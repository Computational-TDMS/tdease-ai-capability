/*!
One physical-constant truth for the whole crate.

Previously the proton mass, neutron spacing, and water mass each existed in
several copies that drifted in the 8th decimal place (`envelope.rs` used
`1.007_276_466_77`, `refine.rs`/`isotopic.rs` used `chemical_elements::PROTON`).
A sub-Da mass shift can be obscured by exactly that kind of drift, so there is
now a single source.
*/

use chemical_elements::PROTON as _CHEM_PROTON;

/// Mass of H+ (proton), used to convert between neutral mass and m/z.
pub const PROTON: f64 = _CHEM_PROTON;
/// 13C–12C mass difference; one neutron's contribution to an isotopologue.
pub const NEUTRON_SHIFT: f64 = 1.003_354_837_8;
/// Neutral water mass; added to residue sums to get an intact peptide mass.
pub const H2O: f64 = 18.010_564_686_3;

/// Convert a measured m/z to a neutral mass, handling positive and negative charge.
pub fn neutral_from_mz(mz: f64, charge: i32, carrier: f64) -> f64 {
    let z = charge.abs() as f64;
    (mz * z) - (charge as f64 * carrier)
}
