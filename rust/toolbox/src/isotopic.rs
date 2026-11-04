/*!
Theoretical isotopic pattern generation.

Port of `mzdeisotope`'s `isotopic_model.rs` (Apache-2.0), which layers
Senko's *Averagine* model over the *BRAIN* recursive isotopic distribution
engine from the `chemical_elements` crate.

Pipeline: `scale(mz, charge)` → averaged integer composition → BRAIN →
exact center-mass peaks → truncate/ignore → shift to m/z.

- [`IsotopicModel`]: generates an isotopic pattern from a mass and charge by
  extrapolating a fractional "average monomer" composition (Senko, 1995).
- [`IsotopicModels`]: a set of named average monomer models for biomolecules.
- [`CachingIsotopicModel`]: same as [`IsotopicModel`] but caches patterns by
  rounded m/z so nearby queries reuse the computed BRAIN distribution.
- [`TheoreticalIsotopicDistributionScalingMethod`]: scales a theoretical
  pattern to observed peak intensities (sum / max / top-3).
- [`BasePeakToMonoisotopicOffsetEstimator`]: estimates how many isotopic
  peaks the base peak sits above the monoisotopic peak for a given mass.
*/

use std::cmp::Ordering;
use std::collections::btree_map::{self, BTreeMap, Entry as BEntry};
use std::collections::hash_map::{self, HashMap};
use std::hash;

use chemical_elements::isotopic_pattern::baffling::NumPeaksSpec;
#[doc(hidden)]
pub use chemical_elements::isotopic_pattern::{
    BafflingRecursiveIsotopicPatternGenerator, TheoreticalIsotopicPattern,
};
use chemical_elements::{
    mass_charge_ratio, neutral_mass, ChemicalComposition, ElementSpecification,
};
use mzpeaks::{CentroidLike, IntensityMeasurement};
use serde::{Deserialize, Serialize};

use crate::error::MassSpecError;

fn isclose(a: f64, b: f64, delta: f64) -> bool {
    (a - b).abs() < delta
}

/// The mass of H+, a hydrogen atom minus an electron.
/// Single source of truth: `crate::constants::PROTON`.
pub use crate::constants::PROTON;

/// A fractional elemental composition with non-ordinal element counts used to represent
/// "averaged" chemical compositions.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct FractionalComposition<'a>(HashMap<ElementSpecification<'a>, f64>);

impl<'a> Extend<(ElementSpecification<'a>, f64)> for FractionalComposition<'a> {
    fn extend<T: IntoIterator<Item = (ElementSpecification<'a>, f64)>>(&mut self, iter: T) {
        <HashMap<ElementSpecification<'a>, f64> as Extend<(ElementSpecification<'a>, f64)>>::extend(
            &mut self.0,
            iter,
        )
    }
}

impl<'a> FromIterator<(ElementSpecification<'a>, f64)> for FractionalComposition<'a> {
    fn from_iter<T: IntoIterator<Item = (ElementSpecification<'a>, f64)>>(iter: T) -> Self {
        let mut this = Self::default();
        this.extend(iter);
        this
    }
}

impl<'a> From<Vec<(ElementSpecification<'a>, f64)>> for FractionalComposition<'a> {
    fn from(value: Vec<(ElementSpecification<'a>, f64)>) -> Self {
        value.into_iter().collect()
    }
}

impl<'a> FractionalComposition<'a> {
    #[inline]
    pub fn new(composition: HashMap<ElementSpecification<'a>, f64>) -> Self {
        Self(composition)
    }

    #[inline]
    pub fn get<Q>(&self, k: &Q) -> Option<&f64>
    where
        ElementSpecification<'a>: std::borrow::Borrow<Q>,
        Q: hash::Hash + Eq + ?Sized,
    {
        self.0.get(k)
    }

    #[inline]
    pub fn iter(&self) -> hash_map::Iter<'_, ElementSpecification<'a>, f64> {
        self.0.iter()
    }

    pub(crate) fn ordered_entries(&self) -> Vec<(ElementSpecification<'a>, f64)> {
        let mut entries = self
            .0
            .iter()
            .map(|(element, count)| (*element, *count))
            .collect::<Vec<_>>();
        entries.sort_by_key(|(element, _)| (element.element.element_number, element.isotope));
        entries
    }

    #[inline]
    pub fn mass(&self) -> f64 {
        self.ordered_entries()
            .iter()
            .map(|(e, c)| e.element.most_abundant_mass * c)
            .sum()
    }

    pub fn insert(&mut self, k: ElementSpecification<'a>, v: f64) -> Option<f64> {
        self.0.insert(k, v)
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn contains_key<Q>(&self, k: &Q) -> bool
    where
        ElementSpecification<'a>: std::borrow::Borrow<Q>,
        Q: hash::Hash + Eq + ?Sized,
    {
        self.0.contains_key(k)
    }
}

/// The capability to generate theoretical isotopic patterns from a given mass
pub trait IsotopicPatternGenerator {
    /// Generate a theoretical isotopic pattern for a given m/z and charge.
    ///
    /// # Arguments
    /// - `mz`: The m/z to use as a point of reference when generating the isotopic pattern
    /// - `charge`: The charge state to use to convert `mz` to a mass
    /// - `charge_carrier`: The mass of the charge carrier, usually a proton
    /// - `truncate_after`: The cumulative abundance percentage of isotopic signal to retain
    /// - `ignore_below`: The minimum abundance percentage of isotopic signal a peak must have to be kept
    fn isotopic_cluster(
        &mut self,
        mz: f64,
        charge: i32,
        charge_carrier: f64,
        truncate_after: f64,
        ignore_below: f64,
    ) -> TheoreticalIsotopicPattern;

    /// Populate a cache, if available, pre-computing isotopic patterns and storing
    /// them away for future use.
    #[allow(unused, clippy::too_many_arguments)]
    fn populate_cache(
        &mut self,
        min_mz: f64,
        max_mz: f64,
        min_charge: i32,
        max_charge: i32,
        charge_carrier: f64,
        truncate_after: f64,
        ignore_below: f64,
    ) {
        let _ = (
            min_mz,
            max_mz,
            min_charge,
            max_charge,
            charge_carrier,
            truncate_after,
            ignore_below,
        );
    }

    /// Get the largest width isotopic pattern this generator has created so far.
    fn largest_isotopic_width(&self) -> f64 {
        f64::INFINITY
    }
}

/// The mass difference between isotopes `C[13]` and `C[12]`. Not precisely universal, but the
/// majority of expected applications are carbon-based.
/// Single source of truth: `crate::constants::NEUTRON_SHIFT`.
pub use crate::constants::NEUTRON_SHIFT;

#[allow(dead_code)] // kept from the mzdeisotope port for the deconvolution stack
const ISOTOPIC_SHIFT: [f64; 10] = [
    NEUTRON_SHIFT / 1.0,
    NEUTRON_SHIFT / 2.0,
    NEUTRON_SHIFT / 3.0,
    NEUTRON_SHIFT / 4.0,
    NEUTRON_SHIFT / 5.0,
    NEUTRON_SHIFT / 6.0,
    NEUTRON_SHIFT / 7.0,
    NEUTRON_SHIFT / 8.0,
    NEUTRON_SHIFT / 9.0,
    NEUTRON_SHIFT / 10.0,
];

/// Get the m/z difference between isotopic peaks at a given charge state
#[allow(dead_code)] // kept from the mzdeisotope port for the deconvolution stack
#[inline(always)]
pub fn isotopic_shift(charge: i32) -> f64 {
    if charge > 0 && charge < 11 {
        ISOTOPIC_SHIFT[(charge - 1) as usize]
    } else {
        NEUTRON_SHIFT / charge as f64
    }
}

/// A package of parameters used to generate isotopic patterns
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct IsotopicPatternParams {
    /// The cumulative abundance percentage of isotopic signal to retain
    pub truncate_after: f64,
    /// The minimum abundance percentage of isotopic signal a peak must have
    /// to be kept
    pub ignore_below: f64,
    /// Whether or not to use the incremental truncation algorithm, and how far
    /// down to truncate
    pub incremental_truncation: Option<f64>,
    /// The mass of the charge carrier, e.g. proton mass
    pub charge_carrier: f64,
}

impl Default for IsotopicPatternParams {
    fn default() -> Self {
        Self {
            truncate_after: 0.95,
            ignore_below: 0.001,
            incremental_truncation: None,
            charge_carrier: PROTON,
        }
    }
}

impl IsotopicPatternParams {
    pub fn new(
        truncate_after: f64,
        ignore_below: f64,
        incremental_truncation: Option<f64>,
        charge_carrier: f64,
    ) -> Self {
        Self {
            truncate_after,
            ignore_below,
            incremental_truncation,
            charge_carrier,
        }
    }
}

/// A model for converting an m/z and a theoretical charge state into a theoretical
/// isotopic pattern based upon an "average monomer" and linear extension.
///
/// This is an implementation of Senko's Averagine [^1]
///
/// # References
/// [^1]: Senko M, Beu S, McLafferty F: Determination of Monoisotopic Masses and Ion
///       Populations for Large Biomolecules from Resolved Isotopic Distributions.
///       Journal of the American Society for Mass Spectrometry 1995, 6:229-233
///       <https://doi.org/10.1016/1044-0305(95)00017-8>
#[derive(Debug, Clone)]
pub struct IsotopicModel<'lifespan> {
    /// The "average" monomer composition
    pub base_composition: FractionalComposition<'lifespan>,
    /// The mass of the average monomer to interpolate with
    pub base_mass: f64,
    hydrogen: ElementSpecification<'lifespan>,
    generator: BafflingRecursiveIsotopicPatternGenerator<'lifespan>,
}

impl PartialEq for IsotopicModel<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.base_composition == other.base_composition
    }
}

impl<'lifespan: 'transient, 'transient> IsotopicModel<'lifespan> {
    /// Create a new [`IsotopicModel`] from a fractional composition
    pub fn new<C: Into<FractionalComposition<'lifespan>>>(base_composition: C) -> Self {
        let base_composition = base_composition.into();
        Self {
            base_mass: base_composition.mass(),
            base_composition,
            hydrogen: ElementSpecification::parse("H").expect("H is a valid element"),
            generator: BafflingRecursiveIsotopicPatternGenerator::new(),
        }
    }

    /// Extrapolate the average monomer composition to an integer composition
    /// for the neutral mass implied by `mz` at `charge`.
    pub fn scale(
        &self,
        mz: f64,
        charge: i32,
        charge_carrier: f64,
    ) -> ChemicalComposition<'transient> {
        let neutral = neutral_mass(mz, charge, charge_carrier);
        let scale = neutral / self.base_mass;

        let mut scaled = ChemicalComposition::new();
        for (elt, count) in self.base_composition.ordered_entries() {
            scaled.set(elt, (count * scale).round() as i32);
        }
        let scaled_mass = scaled.mass();
        let delta = (scaled_mass - neutral).round() as i32;
        let hydrogens = scaled[&self.hydrogen];
        if hydrogens > delta {
            scaled[&self.hydrogen] -= delta;
        } else {
            scaled[&self.hydrogen] = 0;
        }
        scaled
    }
}

impl IsotopicPatternGenerator for IsotopicModel<'_> {
    fn isotopic_cluster(
        &mut self,
        mz: f64,
        charge: i32,
        charge_carrier: f64,
        truncate_after: f64,
        ignore_below: f64,
    ) -> TheoreticalIsotopicPattern {
        let composition = self.scale(mz, charge, charge_carrier);
        let peaks = self.generator.isotopic_variants(
            composition,
            NumPeaksSpec::Guess,
            charge,
            charge_carrier,
        );
        let peaks = TheoreticalIsotopicPattern::from(peaks);
        let diff = mz - peaks.origin;
        peaks
            .truncate_after(truncate_after)
            .ignore_below(ignore_below)
            .shift(diff)
    }
}

impl<T: IntoIterator<Item = (&'static str, f64)>> From<T> for IsotopicModel<'_> {
    fn from(iter: T) -> Self {
        let mut f = FractionalComposition::default();
        for (e, c) in iter {
            f.insert(e.parse().expect("Failed to parse element specification"), c);
        }
        IsotopicModel::new(f)
    }
}

impl<'a> From<IsotopicModel<'a>> for FractionalComposition<'a> {
    fn from(value: IsotopicModel<'a>) -> Self {
        value.base_composition
    }
}

/// A set of named average monomer isotopic models for biomolecules.
///
/// Variants convert to [`IsotopicModel`] and [`CachingIsotopicModel`].
#[derive(Debug, Clone, Copy)]
pub enum IsotopicModels {
    Peptide,
    Glycan,
    Glycopeptide,
    PermethylatedGlycan,
    Heparin,
    HeparanSulfate,
}

impl IsotopicModels {
    pub fn name(self) -> &'static str {
        match self {
            Self::Peptide => "peptide",
            Self::Glycan => "glycan",
            Self::Glycopeptide => "glycopeptide",
            Self::PermethylatedGlycan => "permethylated-glycan",
            Self::Heparin => "heparin",
            Self::HeparanSulfate => "heparan-sulfate",
        }
    }
}

impl From<IsotopicModels> for IsotopicModel<'_> {
    fn from(source: IsotopicModels) -> IsotopicModel<'static> {
        match source {
            IsotopicModels::Peptide => vec![
                ("H", 7.7583f64),
                ("C", 4.9384),
                ("S", 0.0417),
                ("O", 1.4773),
                ("N", 1.3577),
            ],
            IsotopicModels::Glycan => vec![("H", 11.8333), ("C", 7.0), ("O", 5.16666), ("N", 0.5)],
            IsotopicModels::Glycopeptide => vec![
                ("H", 15.75),
                ("C", 10.93),
                ("S", 0.02054),
                ("O", 6.4773),
                ("N", 1.6577),
            ],
            IsotopicModels::PermethylatedGlycan => {
                vec![("C", 12.0), ("H", 21.8333), ("N", 0.5), ("O", 5.16666)]
            }
            IsotopicModels::Heparin => {
                vec![("H", 10.5), ("C", 6.0), ("S", 0.5), ("O", 5.5), ("N", 0.5)]
            }
            IsotopicModels::HeparanSulfate => vec![
                ("H", 10.667),
                ("C", 6.0),
                ("S", 1.333),
                ("O", 9.0),
                ("N", 0.667),
            ],
        }
        .into()
    }
}

impl From<IsotopicModels> for CachingIsotopicModel<'_> {
    fn from(source: IsotopicModels) -> CachingIsotopicModel<'static> {
        let model: IsotopicModel = source.into();
        model.into()
    }
}

/// Strategies for scaling a theoretical isotopic pattern to pair with an experimental
/// isotopic pattern expected to follow the same distribution.
#[derive(Debug, Clone, Copy, Default)]
pub enum TheoreticalIsotopicDistributionScalingMethod {
    #[default]
    /// Assume that the sum of all experimental peaks' intensities are contributed by the same isotopic pattern
    Sum,
    /// Assume the most abundant theoretical peak is the most reliable peak to determine the scaling factor for
    /// the entire theoretical distribution
    Max,
    Top3,
}

impl TheoreticalIsotopicDistributionScalingMethod {
    fn find_top_3_peaks<T: IntensityMeasurement>(
        &self,
        peaks: &[T],
    ) -> ((usize, usize, usize), (f32, f32, f32)) {
        let mut t1_index: usize = 0;
        let mut t2_index: usize = 0;
        let mut t3_index: usize = 0;
        let mut t1 = 0.0f32;
        let mut t2 = 0.0f32;
        let mut t3 = 0.0f32;

        for (i, p) in peaks.iter().enumerate() {
            let y = p.intensity();
            if y > t1 {
                t3 = t2;
                t2 = t1;
                t1 = y;
                t3_index = t2_index;
                t2_index = t1_index;
                t1_index = i;
            } else if y > t2 {
                t3_index = t2_index;
                t3 = t2;
                t2 = y;
                t2_index = i;
            } else if y > t3 {
                t3_index = i;
                t3 = y;
            }
        }

        ((t1_index, t2_index, t3_index), (t1, t2, t3))
    }

    pub fn scale<C: CentroidLike>(
        &self,
        experimental: &[C],
        theoretical: &mut TheoreticalIsotopicPattern,
    ) {
        if theoretical.is_empty() {
            return;
        }
        match self {
            Self::Sum => {
                let total: f32 = experimental.iter().map(|p| p.intensity()).sum();
                theoretical
                    .iter_mut()
                    .for_each(|p| p.intensity *= total as f64);
            }
            Self::Max => {
                let (index, theo_max) = theoretical
                    .iter()
                    .enumerate()
                    .max_by(|a, b| a.1.intensity().total_cmp(&b.1.intensity()))
                    .map(|(i, p)| (i, p.intensity()))
                    .unwrap();
                let scale = experimental[index].intensity() / theo_max;
                theoretical
                    .iter_mut()
                    .for_each(|p| p.intensity *= scale as f64);
            }
            Self::Top3 => {
                let ((t1_index, t2_index, t3_index), (t1, t2, t3)) =
                    self.find_top_3_peaks(&theoretical.peaks);
                // The top-3 indices come from `theoretical`, which `experimental` is
                // matched against one-to-one, so the indices are always in range.
                let scale = {
                    let mut scale = experimental[t1_index].intensity() / t1;
                    scale += experimental[t2_index].intensity() / t2;
                    scale += experimental[t3_index].intensity() / t3;
                    scale / 3.0
                };
                theoretical
                    .iter_mut()
                    .for_each(|p| p.intensity *= scale as f64);
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct IsotopicPatternSpec {
    pub mz: f64,
    pub charge: i32,
    pub charge_carrier: f64,
    pub truncate_after: f64,
    pub ignore_below: f64,
}

impl PartialOrd for IsotopicPatternSpec {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for IsotopicPatternSpec {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.charge.cmp(&other.charge) {
            Ordering::Equal => match self.mz.total_cmp(&other.mz) {
                Ordering::Equal => match self.truncate_after.total_cmp(&other.truncate_after) {
                    Ordering::Equal => match self.ignore_below.total_cmp(&other.ignore_below) {
                        Ordering::Equal => self.charge_carrier.total_cmp(&other.charge_carrier),
                        x => x,
                    },
                    x => x,
                },
                x => x,
            },
            x => x,
        }
    }
}

impl PartialEq for IsotopicPatternSpec {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        isclose(self.mz, other.mz, 1e-6)
            && self.charge == other.charge
            && isclose(self.charge_carrier, other.charge_carrier, 1e-6)
            && isclose(self.truncate_after, other.truncate_after, 1e-6)
            && isclose(self.ignore_below, other.ignore_below, 1e-6)
    }
}

impl Eq for IsotopicPatternSpec {}

#[derive(Debug, Clone, Copy)]
struct FloatRange {
    start: f64,
    end: f64,
    step: f64,
    index: usize,
}

impl FloatRange {
    fn new(start: f64, end: f64, step: f64) -> Self {
        Self {
            start,
            end,
            step,
            index: 0,
        }
    }
}

impl Iterator for FloatRange {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        let val = self.start + self.step * (self.index as f64);
        if val < self.end {
            self.index += 1;
            Some(val)
        } else {
            None
        }
    }
}

/// A wrapper around [`IsotopicModel`] which includes a cache over isotopic patterns,
/// mapping similar m/z values with the same charge and parameters to a previously calculated
/// pattern if one exists, otherwise computing a new pattern and saving it in the cache.
#[derive(Debug, Clone)]
pub struct CachingIsotopicModel<'lifespan> {
    cache_truncation: f64,
    inner: IsotopicModel<'lifespan>,
    cache: BTreeMap<IsotopicPatternSpec, TheoreticalIsotopicPattern>,
    largest_isotopic_width: Option<f64>,
}

fn isotopic_pattern_width(tid: &TheoreticalIsotopicPattern) -> f64 {
    let peaks = &tid.peaks;
    let n = peaks.len();
    if n < 2 {
        return 0.0;
    }
    let first = &peaks[0];
    let last = &peaks[n - 1];
    last.mz - first.mz
}

impl<'lifespan> CachingIsotopicModel<'lifespan> {
    pub fn new<C: Into<FractionalComposition<'lifespan>>>(
        base_composition: C,
        cache_truncation: f64,
    ) -> Self {
        Self {
            inner: IsotopicModel::new(base_composition),
            cache: BTreeMap::new(),
            cache_truncation,
            largest_isotopic_width: Some(0.0),
        }
    }

    pub fn len(&self) -> usize {
        self.cache.len()
    }

    pub fn is_empty(&self) -> bool {
        self.cache.is_empty()
    }

    pub fn iter(&self) -> btree_map::Iter<'_, IsotopicPatternSpec, TheoreticalIsotopicPattern> {
        self.cache.iter()
    }

    pub fn clear(&mut self) {
        self.cache.clear();
    }

    fn get_largest_isotopic_width(&self) -> f64 {
        self.cache
            .values()
            .map(isotopic_pattern_width)
            .max_by(|a, b| a.total_cmp(b))
            .unwrap_or_default()
    }

    pub fn largest_isotopic_width(&self) -> f64 {
        self.largest_isotopic_width
            .unwrap_or_else(|| self.get_largest_isotopic_width())
    }

    #[inline(always)]
    fn truncate_mz(&self, mz: f64) -> f64 {
        (mz / self.cache_truncation).round() * self.cache_truncation
    }

    pub fn make_cache_key(
        &self,
        mz: f64,
        charge: i32,
        charge_carrier: f64,
        truncate_after: f64,
        ignore_below: f64,
    ) -> IsotopicPatternSpec {
        IsotopicPatternSpec {
            mz: self.truncate_mz(mz),
            charge,
            charge_carrier,
            truncate_after,
            ignore_below,
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub fn populate_cache(
        &mut self,
        min_mz: f64,
        max_mz: f64,
        min_charge: i32,
        max_charge: i32,
        charge_carrier: f64,
        truncate_after: f64,
        ignore_below: f64,
    ) {
        let sign = min_charge / min_charge.abs();
        for mz in FloatRange::new(min_mz, max_mz, 0.1) {
            for charge in min_charge.abs()..max_charge.abs() {
                self.isotopic_cluster(
                    mz,
                    charge * sign,
                    charge_carrier,
                    truncate_after,
                    ignore_below,
                );
            }
        }
    }
}

impl IsotopicPatternGenerator for CachingIsotopicModel<'_> {
    fn largest_isotopic_width(&self) -> f64 {
        self.largest_isotopic_width()
    }

    fn isotopic_cluster(
        &mut self,
        mz: f64,
        charge: i32,
        charge_carrier: f64,
        truncate_after: f64,
        ignore_below: f64,
    ) -> TheoreticalIsotopicPattern {
        let key = self.make_cache_key(mz, charge, charge_carrier, truncate_after, ignore_below);
        match self.cache.entry(key) {
            BEntry::Occupied(ent) => {
                let res = ent.get();
                let offset = mz - res.origin;
                res.clone_shifted(offset)
            }
            BEntry::Vacant(ent) => {
                let res = self.inner.isotopic_cluster(
                    key.mz,
                    charge,
                    charge_carrier,
                    truncate_after,
                    ignore_below,
                );
                let offset = mz - res.origin;
                let out = ent.insert(res).clone_shifted(offset);
                let width = isotopic_pattern_width(&out);
                match &mut self.largest_isotopic_width {
                    Some(w) => {
                        *w = w.max(width);
                    }
                    None => {
                        self.largest_isotopic_width = Some(width);
                    }
                }
                out
            }
        }
    }
}

impl<'a> From<IsotopicModel<'a>> for CachingIsotopicModel<'a> {
    fn from(inst: IsotopicModel<'a>) -> CachingIsotopicModel<'a> {
        CachingIsotopicModel::new(inst.base_composition, 1.0)
    }
}

/// An estimator that maps a neutral mass to the number of isotopic peaks by which
/// the base peak exceeds the monoisotopic peak, using only the isotopic model.
///
/// Results are binned by mass to avoid recomputing for every query.
#[derive(Debug, Clone)]
pub struct BasePeakToMonoisotopicOffsetEstimator<I: IsotopicPatternGenerator> {
    isotopic_model: I,
    bins: Vec<usize>,
    step_size: f64,
}

impl<I: IsotopicPatternGenerator> BasePeakToMonoisotopicOffsetEstimator<I> {
    pub fn new(isotopic_model: I, step_size: f64) -> Self {
        Self {
            isotopic_model,
            bins: Vec::new(),
            step_size,
        }
    }

    pub fn len(&self) -> usize {
        self.bins.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bins.is_empty()
    }

    fn max_bin_mass(&self) -> f64 {
        self.len() as f64 * self.step_size
    }

    fn bin_for(&self, mass: f64) -> usize {
        let offset = mass / self.step_size;
        offset as usize
    }

    fn estimate_for_peak_offset(&mut self, mass: f64) -> usize {
        let tid = self
            .isotopic_model
            .isotopic_cluster(mass, 1, PROTON, 0.99, 0.0);
        tid.into_iter()
            .enumerate()
            .max_by(|(_, pa), (_, pb)| pa.intensity.total_cmp(&pb.intensity))
            .map(|(i, _)| i)
            .unwrap()
    }

    fn populate_bins(&mut self, max_mass: f64) {
        let mut current_bin = self.max_bin_mass();

        while max_mass >= current_bin {
            let next_bin_mass = current_bin + self.step_size;
            let delta = self.estimate_for_peak_offset(next_bin_mass);
            self.bins.push(delta);
            current_bin = next_bin_mass;
        }
    }

    pub fn get_peak_offset(&mut self, mass: f64) -> usize {
        let index = self.bin_for(mass);
        if self.len() <= index {
            self.populate_bins(mass);
        }
        self.bins[index]
    }

    pub fn get_peak_offset_direct(&mut self, mass: f64) -> usize {
        self.estimate_for_peak_offset(mass)
    }
}

// ---------------------------------------------------------------------------
// Isotopic query types and entry points
// ---------------------------------------------------------------------------

/// A single theoretical isotopic peak: exact m/z and relative intensity.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct TheoreticalPeakPoint {
    pub mz: f64,
    pub intensity: f64,
}

/// Named averagine models accepted by the native query entry points.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
#[derive(Default)]
pub enum IsotopicModelName {
    #[default]
    Peptide,
    Glycan,
    Glycopeptide,
    #[serde(rename = "permethylated-glycan")]
    PermethylatedGlycan,
    Heparin,
    #[serde(rename = "heparan-sulfate")]
    HeparanSulfate,
}

impl From<IsotopicModelName> for IsotopicModels {
    fn from(value: IsotopicModelName) -> Self {
        match value {
            IsotopicModelName::Peptide => IsotopicModels::Peptide,
            IsotopicModelName::Glycan => IsotopicModels::Glycan,
            IsotopicModelName::Glycopeptide => IsotopicModels::Glycopeptide,
            IsotopicModelName::PermethylatedGlycan => IsotopicModels::PermethylatedGlycan,
            IsotopicModelName::Heparin => IsotopicModels::Heparin,
            IsotopicModelName::HeparanSulfate => IsotopicModels::HeparanSulfate,
        }
    }
}

impl IsotopicModelName {
    pub fn name(self) -> &'static str {
        IsotopicModels::from(self).name()
    }
}

/// Input for generating a theoretical isotopic pattern.
///
/// Provide exactly one of `mz` (an observed m/z, which fixes the monoisotopic
/// peak position) or `mass` (a neutral monoisotopic mass).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct IsotopicPatternQuery {
    /// The m/z whose monoisotopic peak the generated pattern is aligned to.
    pub mz: Option<f64>,
    /// A neutral monoisotopic mass; converted to m/z using `charge` and the proton mass.
    pub mass: Option<f64>,
    /// The charge state of the ion (positive or negative, must not be zero).
    pub charge: i32,
    /// The averagine model to use. Defaults to `peptide`.
    #[serde(default)]
    pub model: IsotopicModelName,
    /// The cumulative abundance of isotopic signal to retain. Defaults to 0.95.
    #[serde(default = "default_truncate_after")]
    pub truncate_after: f64,
    /// The minimum relative abundance a peak must have to be kept. Defaults to 0.001.
    #[serde(default = "default_ignore_below")]
    pub ignore_below: f64,
}

fn default_truncate_after() -> f64 {
    0.95
}

fn default_ignore_below() -> f64 {
    0.001
}

/// The generated theoretical isotopic pattern.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct IsotopicPatternResult {
    /// The theoretical peaks, in ascending m/z order. Intensities are relative
    /// (sum to 1.0 over the retained pattern).
    pub peaks: Vec<TheoreticalPeakPoint>,
    /// The averagine-derived integer composition used to compute the pattern.
    pub composition: Vec<ElementCount>,
    /// The model name used to generate the pattern.
    pub model: String,
    /// The m/z the pattern was generated at.
    pub mz: f64,
    /// The charge state used.
    pub charge: i32,
}

/// A single element's count in a composition.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ElementCount {
    pub element: String,
    pub count: i32,
    /// Mass number when this count is a specific isotope (0 = the element's
    /// natural/monoisotopic form). Stable-isotope labels from UniMod (`2H`,
    /// `13C`) carry this so the exact composition keeps its isotope identity —
    /// mass and envelope computations must use the isotope mass, not the
    /// natural-abundance one.
    #[serde(default)]
    pub isotope: u16,
}

/// Generate a theoretical isotopic pattern from an m/z (or neutral mass) and charge
/// using the named averagine model.
pub fn generate_isotopic_pattern(
    query: &IsotopicPatternQuery,
) -> Result<IsotopicPatternResult, MassSpecError> {
    validate_query(query)?;
    let mz = resolve_mz(query)?;

    let mut model: IsotopicModel = IsotopicModels::from(query.model).into();
    let tid = model.isotopic_cluster(
        mz,
        query.charge,
        PROTON,
        query.truncate_after,
        query.ignore_below,
    );

    let composition = model.scale(mz, query.charge, PROTON);
    let peaks: Vec<TheoreticalPeakPoint> = tid
        .peaks
        .iter()
        .map(|p| TheoreticalPeakPoint {
            mz: p.mz,
            intensity: p.intensity,
        })
        .collect();
    let composition: Vec<ElementCount> = composition
        .iter()
        .map(|(elt, count)| ElementCount {
            element: elt.element.symbol.clone(),
            count: *count,
            isotope: elt.isotope,
        })
        .collect();

    Ok(IsotopicPatternResult {
        peaks,
        composition,
        model: query.model.name().to_owned(),
        mz,
        charge: query.charge,
    })
}

/// Generate a theoretical isotopic pattern from an explicit elemental formula.
///
/// Unlike [`generate_isotopic_pattern`], this does not use an averagine model;
/// the pattern is computed exactly for the given composition.
pub fn generate_isotopic_pattern_from_formula(
    formula: &str,
    charge: i32,
    truncate_after: f64,
    ignore_below: f64,
) -> Result<IsotopicPatternResult, MassSpecError> {
    if charge <= 0 {
        return Err(MassSpecError::InvalidQuery(
            "isotopicPattern charge must be positive".to_owned(),
        ));
    }
    let composition = ChemicalComposition::parse(formula).map_err(|error| {
        MassSpecError::InvalidQuery(format!("invalid elemental formula {formula}: {error}"))
    })?;

    let mut generator = BafflingRecursiveIsotopicPatternGenerator::new();
    let peaks =
        generator.isotopic_variants(composition.clone(), NumPeaksSpec::Guess, charge, PROTON);
    let tid = TheoreticalIsotopicPattern::from(peaks)
        .truncate_after(truncate_after)
        .ignore_below(ignore_below);

    let mz = tid.peaks.first().map(|p| p.mz).unwrap_or_default();
    let peaks: Vec<TheoreticalPeakPoint> = tid
        .peaks
        .iter()
        .map(|p| TheoreticalPeakPoint {
            mz: p.mz,
            intensity: p.intensity,
        })
        .collect();
    let composition: Vec<ElementCount> = composition
        .iter()
        .map(|(elt, count)| ElementCount {
            element: elt.element.symbol.clone(),
            count: *count,
            isotope: elt.isotope,
        })
        .collect();

    Ok(IsotopicPatternResult {
        peaks,
        composition,
        model: "formula".to_owned(),
        mz,
        charge,
    })
}

fn validate_query(query: &IsotopicPatternQuery) -> Result<(), MassSpecError> {
    if query.charge <= 0 {
        return Err(MassSpecError::InvalidQuery(
            "isotopicPattern charge must be positive".to_owned(),
        ));
    }
    if query.mz.is_some() == query.mass.is_some() {
        return Err(MassSpecError::InvalidQuery(
            "isotopicPattern requires exactly one of mz or mass".to_owned(),
        ));
    }
    if !(0.0..=1.0).contains(&query.truncate_after) {
        return Err(MassSpecError::InvalidQuery(
            "isotopicPattern truncateAfter must be within (0, 1]".to_owned(),
        ));
    }
    if !(0.0..=1.0).contains(&query.ignore_below) {
        return Err(MassSpecError::InvalidQuery(
            "isotopicPattern ignoreBelow must be within [0, 1)".to_owned(),
        ));
    }
    Ok(())
}

fn resolve_mz(query: &IsotopicPatternQuery) -> Result<f64, MassSpecError> {
    if let Some(mz) = query.mz {
        if !mz.is_finite() || mz <= 0.0 {
            return Err(MassSpecError::InvalidQuery(
                "isotopicPattern mz must be positive".to_owned(),
            ));
        }
        Ok(mz)
    } else if let Some(mass) = query.mass {
        if !mass.is_finite() || mass <= 0.0 {
            return Err(MassSpecError::InvalidQuery(
                "isotopicPattern mass must be positive".to_owned(),
            ));
        }
        Ok(mass_charge_ratio(mass, query.charge, PROTON))
    } else {
        Err(MassSpecError::InvalidQuery(
            "isotopicPattern requires exactly one of mz or mass".to_owned(),
        ))
    }
}

#[cfg(test)]
mod test {
    use mzpeaks::{CentroidPeak, MZPeakSetType};

    use super::*;

    fn assert_is_close(t1: f64, t2: f64, tol: f64, label: &str) {
        assert!(
            (t1 - t2).abs() < tol,
            "Observed {label} {t1}, expected {t2}, difference {}",
            t1 - t2,
        );
    }

    #[test]
    fn test_fc() {
        let model: IsotopicModel = IsotopicModels::Peptide.into();
        assert!(model.base_composition.contains_key("C"));
        assert_eq!(model.base_composition.get("C").copied(), Some(4.9384));
        assert_eq!(model.base_composition.len(), 5);
    }

    #[test]
    fn test_tid() {
        let mut model: IsotopicModel = IsotopicModels::Peptide.into();
        let tid = model
            .isotopic_cluster(1000.0, 2, PROTON, 0.95, 0.001)
            .scale_by(100.0);

        assert!((tid.total() - 100.0).abs() < 1e-6);
        assert!((tid[0].mz() - 1000.0).abs() < 1e-6);

        let diff = tid[1].mz() - tid[0].mz();
        let neutron2 = 0.5014313;
        assert!((diff - neutron2).abs() < 1e-6);

        // Values reflect the BRAIN isotope-indexing fix in `chemical_elements`
        // (isotope_by_shift), which corrected the C/O/N/S isotope placement.
        let expected = 31.2953548;
        assert_is_close(
            f64::from(tid[0].intensity()),
            expected,
            1e-6,
            "M0 intensity",
        );

        let tid = model
            .isotopic_cluster(1000.5, 2, PROTON, 0.95, 0.001)
            .scale_by(100.0);
        assert!((tid.total() - 100.0).abs() < 1e-6);
        assert!((tid[0].mz() - 1000.5).abs() < 1e-6);

        let expected = 31.2918873;
        assert_is_close(
            f64::from(tid[0].intensity()),
            expected,
            1e-6,
            "M0 intensity",
        );
    }

    #[test]
    fn test_scale() {
        let mut isomod: IsotopicModel = IsotopicModels::Heparin.into();
        let mut ref_tid = isomod.isotopic_cluster(1500.0, 2, PROTON, 0.95, 0.1);
        ref_tid = ref_tid.scale_by(1000.0);
        let eid: MZPeakSetType<_> = ref_tid
            .iter()
            .map(|p| CentroidPeak::new(p.mz(), p.intensity(), 0))
            .collect();

        let mut tid = isomod.isotopic_cluster(1500.0, 2, PROTON, 0.95, 0.1);
        TheoreticalIsotopicDistributionScalingMethod::Sum.scale(eid.as_slice(), &mut tid);

        let tic_ref: f32 = ref_tid.iter().map(|p| p.intensity()).sum();
        let tic_sum: f32 = tid.iter().map(|p| p.intensity()).sum();
        assert!((tic_ref - tic_sum).abs() < 1e-3);

        let mut tid = isomod.isotopic_cluster(1500.0, 2, PROTON, 0.95, 0.1);
        TheoreticalIsotopicDistributionScalingMethod::Max.scale(eid.as_slice(), &mut tid);

        let tic_max: f32 = tid.iter().map(|p| p.intensity()).sum();
        assert!((tic_ref - tic_max).abs() < 1e-3);

        let mut tid = isomod.isotopic_cluster(1500.0, 2, PROTON, 0.95, 0.1);
        TheoreticalIsotopicDistributionScalingMethod::Top3.scale(eid.as_slice(), &mut tid);

        let tic_top3: f32 = tid.iter().map(|p| p.intensity()).sum();
        assert!(
            (tic_ref - tic_top3).abs() < 1e-3,
            "Observed TIC {tic_top3}, expected TIC {tic_ref}, delta: {}",
            tic_ref - tic_top3
        );
    }

    #[test]
    fn test_cache() {
        let mut model: CachingIsotopicModel = IsotopicModels::Peptide.into();
        assert_eq!(model.get_largest_isotopic_width(), 0.0);
        assert!(model.is_empty());

        model.populate_cache(200.0, 400.0, 2, 4, PROTON, 0.95, 0.001);
        assert_is_close(
            model.get_largest_isotopic_width(),
            1.0027695470485014,
            1e-6,
            "width",
        );

        model.clear();
        assert_eq!(model.get_largest_isotopic_width(), 0.0);
        assert!(model.is_empty());

        let mut model: IsotopicModel = IsotopicModels::Peptide.into();
        assert_eq!(model.largest_isotopic_width(), f64::INFINITY);
        model.populate_cache(200.0, 400.0, 2, 4, PROTON, 0.95, 0.001);
        assert_eq!(model.largest_isotopic_width(), f64::INFINITY);
    }

    #[test]
    fn test_tid_cached() {
        let mut model: CachingIsotopicModel = IsotopicModels::Peptide.into();
        let tid = model
            .isotopic_cluster(1000.0, 2, PROTON, 0.95, 0.001)
            .scale_by(100.0);

        assert!((tid.total() - 100.0).abs() < 1e-6);
        assert!((tid[0].mz() - 1000.0).abs() < 1e-6);

        let expected = 31.2953548;
        assert_is_close(
            f64::from(tid[0].intensity()),
            expected,
            1e-6,
            "M0 intensity",
        );

        let tid = model
            .isotopic_cluster(1000.001, 2, PROTON, 0.95, 0.001)
            .scale_by(100.0);

        assert!((tid.total() - 100.0).abs() < 1e-6);
        assert!((tid[0].mz() - 1000.001).abs() < 1e-6);
        assert_is_close(
            f64::from(tid[0].intensity()),
            expected,
            1e-6,
            "M0 intensity",
        );
        assert_eq!(model.len(), 1);
    }

    #[test]
    fn test_mono_to_base() {
        let mut model: BasePeakToMonoisotopicOffsetEstimator<IsotopicModel> =
            BasePeakToMonoisotopicOffsetEstimator::new(IsotopicModels::Peptide.into(), 10.0);

        assert_eq!(model.get_peak_offset(1000.0), 0);
        assert_eq!(model.get_peak_offset(2500.0), 1);
        assert_eq!(model.get_peak_offset(10000.0), 6);
    }

    #[test]
    fn test_isotopic_pattern_overflow() {
        let mut model: IsotopicModel = IsotopicModels::Peptide.into();
        let tid = model.isotopic_cluster(40000.0, 1, PROTON, 0.95, 0.001);
        assert!(tid.len() > 5, "{} should be longer than 5 peaks", tid.len())
    }

    #[test]
    fn test_generate_isotopic_pattern_query() {
        let query = IsotopicPatternQuery {
            mz: Some(1000.0),
            mass: None,
            charge: 2,
            model: IsotopicModelName::Peptide,
            truncate_after: 0.95,
            ignore_below: 0.001,
        };
        let result = generate_isotopic_pattern(&query).unwrap();
        assert!(!result.peaks.is_empty());
        assert!((result.peaks[0].mz - 1000.0).abs() < 1e-6);
        assert_eq!(result.model, "peptide");
        let total: f64 = result.peaks.iter().map(|p| p.intensity).sum();
        assert!((total - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_generate_isotopic_pattern_from_formula() {
        let result = generate_isotopic_pattern_from_formula("C6H12O6", 1, 0.95, 0.001).unwrap();
        assert!(!result.peaks.is_empty());
        // Glucose monoisotopic m/z ~ 180.063
        assert!((result.peaks[0].mz - 181.07).abs() < 0.01);
        assert_eq!(result.model, "formula");
    }

    #[test]
    fn brominated_fragment_doublet_is_discriminating_against_phospho() {
        // Real-element discrimination (T1 ground truth, J. Proteome Res.
        // 2017): a fragment carrying Y88-Br vs the same fragment carrying
        // phospho. The whole-protein MS1 envelope cannot separate them (see
        // refine.rs `real_data_ms1_envelope_cannot_shape_discriminate_br_phospho`),
        // but a small fragment can: 81Br is ~49% abundant, so the +2 peak of a
        // Br-carrying fragment is ~as strong as M0, while a phospho fragment's
        // +2 is C13-only (~10%). This is the fragment-level discrimination the
        // published ground truth uses ("isotope pattern of fragments containing
        // Y88 shows a distinctive doublet absent from fragments that do not").
        //
        // A Br-carrying tyrosine fragment (Tyr residue C9H9NO2, Br replaces an
        // aromatic H -> C9H8NO2Br) — the +2 relative intensity must clearly
        // exceed the phospho form. Phospho on the residue = C9H9NO2 + PO3H
        // (the phosphate on the phenolic OH) = C9H10NO6P.
        let br = generate_isotopic_pattern_from_formula("C9H8NO2Br", 1, 0.99, 0.001).unwrap();
        let po3h = generate_isotopic_pattern_from_formula("C9H10NO6P", 1, 0.99, 0.001).unwrap();

        let br_m0 = br.peaks[0].intensity;
        // +2 peak of the Br fragment (isotope index 2).
        let br_m2 = br.peaks.get(2).map(|p| p.intensity).unwrap_or(0.0);
        let po3h_m0 = po3h.peaks[0].intensity;
        let po3h_m2 = po3h.peaks.get(2).map(|p| p.intensity).unwrap_or(0.0);

        let br_ratio = br_m2 / br_m0;
        let po3h_ratio = po3h_m2 / po3h_m0;
        eprintln!("fragment M+2/M0: Br = {br_ratio:.3}, PO3H = {po3h_ratio:.3}");
        // Br doublet: 81Br gives an ~equal-intensity +2; phospho's +2 is C13
        // only and stays well below 0.5.
        assert!(
            br_ratio > 0.5,
            "Br-carrying fragment must show a prominent +2 doublet (got {br_ratio:.3})",
        );
        assert!(
            po3h_ratio < br_ratio - 0.3,
            "phospho fragment +2 ({po3h_ratio:.3}) must be clearly weaker than Br ({br_ratio:.3})",
        );
    }

    #[test]
    fn brain_brainpy_isotopic_reference_parity() {
        // Relative intensity ratios (M+1/M0, M+2/M0) computed from
        // brain-isotopic-distribution 1.5.19 (brainpy), the Python engine used
        // by isoforma/scripts/quantify.py. Ratios are scale-invariant, so the
        // two engines' different truncation/normalisation (brainpy fixes the
        // peak count; BRAIN truncates by cumulative abundance) cannot hide a
        // real isotope-shape mismatch. Mirrored in tests/test_isoforma_parity.py.
        for (formula, charge, m1_m0, m2_m0) in [
            ("C9H8NO2Br", 1, 0.102678, 0.981601),
            ("C9H10NO6P", 1, 0.104429, 0.017249),
            ("C35H63N15O12", 1, 0.445163, 0.121599),
            ("C100H150N30O35S1", 2, 1.229660, 0.866611),
        ] {
            let result =
                generate_isotopic_pattern_from_formula(formula, charge, 0.9999, 0.0001).unwrap();
            assert!(
                result.peaks.len() >= 3,
                "{formula}: expected >=3 isotopologues, got {}",
                result.peaks.len(),
            );
            let m0 = result.peaks[0].intensity;
            let got_m1 = result.peaks[1].intensity / m0;
            let got_m2 = result.peaks[2].intensity / m0;
            assert!(
                (got_m1 - m1_m0).abs() < 0.01,
                "{formula} M+1/M0: BRAIN {got_m1:.5} != brainpy reference {m1_m0}",
            );
            assert!(
                (got_m2 - m2_m0).abs() < 0.01,
                "{formula} M+2/M0: BRAIN {got_m2:.5} != brainpy reference {m2_m0}",
            );
        }
    }
}
