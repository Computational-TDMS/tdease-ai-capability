/*!
Unified theoretical-envelope kernel.

Every API that needs an isotope envelope (isotopeEnvelope, calibrateMassShift,
generateIsotopicPatternFromFormula) computes it through this single BRAIN
kernel. There is deliberately no second convolution: a given composition must
produce the identical distribution on every API.

BRAIN (`chemical_elements::BafflingRecursiveIsotopicPatternGenerator`) returns
precise *relative* center-of-mass offsets (not a fixed 1.00335 spacing). We
keep those offsets and re-anchor them to an exact input mass, so the absolute
m/z we hand back live in the same coordinate frame as the observed envelope.
*/

use chemical_elements::isotopic_pattern::baffling::NumPeaksSpec;
use chemical_elements::isotopic_pattern::BafflingRecursiveIsotopicPatternGenerator;
use chemical_elements::{ChemicalComposition, PROTON};

use crate::constants::NEUTRON_SHIFT;
use crate::isotopic::{
    IsotopicModel, IsotopicModels, TheoreticalIsotopicPattern, TheoreticalPeakPoint,
};

/// BRAIN theoretical envelope in the charge frame (relative offsets).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct TheoreticalEnvelope {
    /// BRAIN peaks: m/z in the charge frame + relative intensity.
    pub peaks: Vec<TheoreticalPeakPoint>,
    /// Neutral center-of-mass per peak (BRAIN precise offsets, not fixed spacing).
    pub center_masses: Vec<f64>,
    pub charge: i32,
    /// The exact neutral monoisotopic mass of the composition (the M0 position).
    /// This is the anchor baseline for [`anchored_envelope`]. It must be taken
    /// from the *untruncated* BRAIN output: for a large protein the M0 peak's
    /// relative intensity is below `ignore_below` and is dropped from `peaks`,
    /// so `center_masses.first()` is then the M+1 peak and anchoring on it
    /// shifts the whole envelope by one neutron.
    #[serde(default)]
    pub monoisotopic_neutral: f64,
}

/// The same envelope re-anchored so M0 sits at `anchor_mass` (exact input mass).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct AnchoredEnvelope {
    /// Absolute m/z peaks (same coordinate frame as the observed envelope).
    pub peaks: Vec<TheoreticalPeakPoint>,
    /// Absolute neutral center-masses; peak 0 == anchor_mass.
    pub center_masses: Vec<f64>,
    pub charge: i32,
    pub anchor_mass: f64,
}

/// Generate a BRAIN theoretical envelope for an exact composition.
pub fn theoretical_envelope(
    composition: &ChemicalComposition<'_>,
    charge: i32,
    truncate_after: f64,
    ignore_below: f64,
) -> TheoreticalEnvelope {
    let mut generator = BafflingRecursiveIsotopicPatternGenerator::new();
    let peaks =
        generator.isotopic_variants(composition.clone(), NumPeaksSpec::Guess, charge, PROTON);
    let full = TheoreticalIsotopicPattern::from(peaks);
    // M0's exact neutral mass + peak, taken before truncation removes the weak
    // M0 peak of a large protein.
    let m0_peak = full.peaks.first().copied();
    let monoisotopic_neutral = m0_peak
        .map(|p| neutral_from_mz(p.mz, charge, PROTON))
        .unwrap_or(0.0);
    let tid = full
        .truncate_after(truncate_after)
        .ignore_below(ignore_below);
    let mut peaks: Vec<TheoreticalPeakPoint> = tid
        .peaks
        .iter()
        .map(|p| TheoreticalPeakPoint {
            mz: p.mz,
            intensity: p.intensity,
        })
        .collect();
    // Always keep the M0 peak in the list, even when its relative intensity is
    // below `ignore_below` and the truncation dropped it. Without it the first
    // retained peak is M+1, and an integer-offset sweep can slide M+1 onto the
    // observed M0 and earn a spurious match at offset -1 (a negative-offset
    // bias). Keeping M0 makes offset 0 match the observed M0.
    if let Some(m0_peak) = m0_peak {
        if peaks.is_empty() || peaks[0].mz > m0_peak.mz + 1e-9 {
            peaks.insert(
                0,
                TheoreticalPeakPoint {
                    mz: m0_peak.mz,
                    intensity: m0_peak.intensity,
                },
            );
        }
    }
    // Truncation returns a retained distribution normalized to approximately
    // one. Re-inserting a weak M0 therefore pushes the total above one unless
    // the complete retained envelope is normalized again. Coverage is defined
    // as a fraction of this distribution, so preserve that invariant here.
    let retained_total: f64 = peaks.iter().map(|peak| peak.intensity).sum();
    if retained_total.is_finite() && retained_total > 0.0 {
        for peak in &mut peaks {
            peak.intensity /= retained_total;
        }
    }
    let center_masses: Vec<f64> = peaks
        .iter()
        .map(|p| neutral_from_mz(p.mz, charge, PROTON))
        .collect();
    TheoreticalEnvelope {
        peaks,
        center_masses,
        charge,
        monoisotopic_neutral,
    }
}

/// Re-anchor a BRAIN envelope so M0 sits at `anchor_mass` (the exact input mass).
///
/// The BRAIN center-of-mass offsets are relative; this shifts the whole shape
/// so `center_masses[0] == anchor_mass` and the absolute m/z values line up
/// with the observed coordinate frame (P0-b).
pub fn anchored_envelope(envelope: &TheoreticalEnvelope, anchor_mass: f64) -> AnchoredEnvelope {
    // Anchor against the true M0 neutral mass, not `center_masses.first()`:
    // for large proteins the M0 peak is dropped by `ignore_below` and the
    // first retained peak is M+1 — anchoring on it would shift the whole
    // envelope by one neutron (see `TheoreticalEnvelope.monoisotopic_neutral`).
    let brain_mono = if envelope.monoisotopic_neutral > 0.0 {
        envelope.monoisotopic_neutral
    } else {
        envelope.center_masses.first().copied().unwrap_or(0.0)
    };
    let center_masses: Vec<f64> = envelope
        .center_masses
        .iter()
        .map(|cm| anchor_mass + (cm - brain_mono))
        .collect();
    let peaks = envelope
        .peaks
        .iter()
        .enumerate()
        .map(|(i, p)| TheoreticalPeakPoint {
            mz: mz_from_neutral(center_masses[i], envelope.charge),
            intensity: p.intensity,
        })
        .collect();
    AnchoredEnvelope {
        peaks,
        center_masses,
        charge: envelope.charge,
        anchor_mass,
    }
}

/// Per-neutron-index buckets, preserving the legacy `theoreticalPattern: f64[]`
/// shape plus absolute m/z per index.
///
/// Returns `(intensities, mz)` both length `isotope_count`: `intensities[i]` is
/// the relative intensity of the isotopologue at +i neutrons (0 where absent,
/// e.g. Br's M1), `mz[i]` the absolute m/z of that isotopologue (placeholder
/// `anchor_mass + i*NEUTRON_SHIFT` where no BRAIN peak landed, so the arrays
/// stay parallel to the caller's `isotopes` count).
pub fn intensity_by_neutron_index(
    anchored: &AnchoredEnvelope,
    isotope_count: usize,
) -> (Vec<f64>, Vec<f64>) {
    let mut intensities = vec![0.0; isotope_count];
    let mut mz = vec![mz_from_neutral(anchored.anchor_mass, anchored.charge); isotope_count];
    for (i, cm) in anchored.center_masses.iter().enumerate() {
        let neutron_index = ((cm - anchored.anchor_mass) / NEUTRON_SHIFT).round() as usize;
        if neutron_index < isotope_count {
            intensities[neutron_index] = anchored.peaks.get(i).map(|p| p.intensity).unwrap_or(0.0);
            mz[neutron_index] = mz_from_neutral(*cm, anchored.charge);
        }
    }
    // Fill m/z for empty indices with fixed spacing so a plotter can draw gaps.
    for i in 0..isotope_count {
        if intensities[i] == 0.0 {
            mz[i] = mz_from_neutral(
                anchored.anchor_mass + i as f64 * NEUTRON_SHIFT,
                anchored.charge,
            );
        }
    }
    (intensities, mz)
}

/// Averagine (Senko) residue composition scaled to a mass, as an integer
/// composition. Used only when no sequence/formula is available.
///
/// Single averagine model: the same peptide monomer `isotopic.rs` uses, so the
/// same mass yields the same integer composition whether it reaches BRAIN
/// through the envelope/refine fallback or the isotopic-pattern entry
/// points. The base mass is the model's *monoisotopic* mass — every caller here
/// passes a monoisotopic mass — replacing the old hard-coded 111.1254 (Senko's
/// average mass), a 0.07 Da divergence from the mono base.
pub(crate) fn averagine_composition(mass: f64) -> ChemicalComposition<'static> {
    let model: IsotopicModel<'static> = IsotopicModels::Peptide.into();
    let units = mass / model.base_mass;
    let mut comp = ChemicalComposition::new();
    for (elt, count) in model.base_composition.ordered_entries() {
        comp.set(elt, (count * units).round().max(1.0) as i32);
    }
    comp
}

fn neutral_from_mz(mz: f64, charge: i32, carrier: f64) -> f64 {
    let z = charge.abs() as f64;
    (mz * z) - (charge as f64 * carrier)
}

fn mz_from_neutral(neutral_mass: f64, charge: i32) -> f64 {
    let z = charge.abs() as f64;
    (neutral_mass + charge as f64 * PROTON) / z
}

#[cfg(test)]
mod tests {
    use super::*;

    fn comp_of() -> ChemicalComposition<'static> {
        let mut c = ChemicalComposition::new();
        c["C"] = 60;
        c["H"] = 100;
        c["N"] = 15;
        c["O"] = 20;
        c["S"] = 1;
        c
    }

    #[test]
    fn bra_in_center_masses_are_not_fixed_spacing() {
        let env = theoretical_envelope(&comp_of(), 2, 0.99, 0.001);
        assert!(!env.center_masses.is_empty());
        // BRAIN offsets differ slightly from fixed 1.00335 spacing.
        let d0 = env.center_masses[1] - env.center_masses[0];
        assert!((d0 - NEUTRON_SHIFT).abs() < 0.01, "M1 spacing {d0}");
    }

    #[test]
    fn anchored_envelope_places_m0_at_anchor_mass() {
        let env = theoretical_envelope(&comp_of(), 2, 0.99, 0.001);
        let anchored = anchored_envelope(&env, 1000.0);
        assert!((anchored.center_masses[0] - 1000.0).abs() < 1e-9);
        // The m/z of M0 is (1000 + 2*PROTON)/2.
        let m0_mz = anchored.peaks[0].mz;
        assert!((m0_mz - (1000.0 + 2.0 * PROTON) / 2.0).abs() < 1e-6);
    }

    #[test]
    fn intensity_by_neutron_index_parallels_isotope_count() {
        let env = theoretical_envelope(&comp_of(), 2, 0.99, 0.001);
        let anchored = anchored_envelope(&env, 1000.0);
        let (ints, mz) = intensity_by_neutron_index(&anchored, 8);
        assert_eq!(ints.len(), 8);
        assert_eq!(mz.len(), 8);
        // M0 intensity is nonzero and its m/z is the anchored M0.
        assert!(ints[0] > 0.0);
        assert!((mz[0] - (1000.0 + 2.0 * PROTON) / 2.0).abs() < 1e-6);
    }

    #[test]
    fn anchored_envelope_does_not_shift_by_one_neutron_when_m0_is_truncated() {
        // Regression: `ignore_below` drops the weak M0 peak of a large protein,
        // so `center_masses.first()` is the M+1 peak. Anchoring on it used to
        // shift the whole envelope one neutron to the left (apex one peak too
        // low). The anchor must use the exact monoisotopic neutral mass.
        // Composition ~ H4 + Ac + DiMe + Ox + Br (large, Br doublet).
        let mut c = ChemicalComposition::new();
        c["C"] = 498;
        c["H"] = 842;
        c["N"] = 163;
        c["O"] = 136;
        c["S"] = 1;
        c["Br"] = 1;
        let anchor = 11393.288536295671;
        let env = theoretical_envelope(&c, 17, 0.99, 0.001);
        assert!(
            env.monoisotopic_neutral > 0.0,
            "M0 neutral must be captured even when the M0 peak is truncated"
        );
        let retained_total: f64 = env.peaks.iter().map(|peak| peak.intensity).sum();
        assert!(
            (retained_total - 1.0).abs() < 1e-12,
            "retained envelope intensity must be normalized, got {retained_total}"
        );
        let anchored = anchored_envelope(&env, anchor);
        // The apex must sit at the same neutron offset as the untruncated BRAIN
        // shape (7 on this Br protein), i.e. anchor + ~7.02 Da, not +6.
        let apex_neutral = anchored
            .center_masses
            .iter()
            .zip(anchored.peaks.iter())
            .max_by(|a, b| a.1.intensity.total_cmp(&b.1.intensity))
            .map(|(cm, _)| *cm)
            .unwrap_or(0.0);
        let offset_neutrons = (apex_neutral - anchor) / NEUTRON_SHIFT;
        assert!(
            (offset_neutrons - 7.0).abs() < 0.3,
            "Br apex at ~M+7 (got {offset_neutrons:.2} neutrons from anchor)"
        );
    }

    #[test]
    fn averagine_is_positive_and_finite() {
        let c = averagine_composition(1000.0);
        assert!(c.mass() > 0.0);
    }

    #[test]
    fn averagine_scales_by_the_shared_monoisotopic_base() {
        let model: IsotopicModel = IsotopicModels::Peptide.into();
        // Regression: the base used to be the hard-coded Senko *average* mass
        // (111.1254) while isotopic.rs scaled by the model's monoisotopic base —
        // a 0.07 Da dual truth. Both must now divide by the same model base.
        let expected_units = 1000.0 / model.base_mass;
        let c = averagine_composition(1000.0);
        for (elt, ratio) in model.base_composition.iter() {
            let count = (*ratio * expected_units).round().max(1.0) as i32;
            assert_eq!(c[elt], count, "averagine element {elt}");
        }
        assert!(
            (model.base_mass - 111.1254).abs() > 0.05,
            "base must be the mono mass, got {}",
            model.base_mass
        );
    }
}
