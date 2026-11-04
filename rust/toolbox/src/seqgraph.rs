/*!
Sequence-graph modification search — the Informed-Proteomics core, ported.

Source of record: Park et al., *Informed-Proteomics*, Nat. Methods 14, 909
(2017); C# sources `Backend/Data/Sequence/{SequenceGraph,ShiftedSequenceGraph,
ModificationParams,AminoAcidSet}.cs`, `Scoring/TopDown/FilteredProteinMass
Binning.cs`, `TopDown/Scoring/CompositeScorer*.cs` (repo snapshot
2024). The port is line-faithful where the original is load-bearing, and
departs in exactly two documented places (bottom of this comment).

The idea (why this replaces enumeration search): the sequence is built into a
layered graph whose nodes are (sequence position × modification-combination
state), deduplicated BY COMBINATION — never by site arrangement. One memoized
max-DP over the graph then yields, in a single pass, the best placement of
every modification in a combination plus its traceback (which modification at
which position), with prefix compositions shared by all arrangements. The
combination space itself is precompiled (`ModificationParams`): every
multiset of library modifications up to `max_mods`, with an O(1) transition
table (prev combination, added mod) -> new combination.

Integer mass space (the "atoms replace residue masses" trick):
- residues and modifications are `GraphComposition`s — integer atom counts
  (C/H/N/O/S/P plus an additional-element list), mass being a lazy projection;
- `RESCALING_HIGH_PRECISION = 274.335215` maps masses onto fine integer bins
  (~3.645 mDa), and peak lookups are integer-bin hash hits, never float
  window scans.

Departures from the C# original (both deliberate, both cited):
1. Free-shift handling. IP has TWO graphs: `SequenceGraph` (closed
   combinations only; the precursor mass filter runs over sink compositions)
   and `ShiftedSequenceGraph`, whose layer 0 is a pseudo-residue "Shift"
   carrying `CompositionWithDeltaMass(shiftedMass)` — used exclusively by
   `TagMatchFinder` to hang one free mass at a FIXED anchor (a tag boundary).
   Our problem — a large TopPIC-observed shift decomposed over a library with
   an unexplained residual placed anywhere — needs the free mass to float
   across the sequence TOGETHER with pinned library mods, so we express it
   as ONE designated virtual modification `SHIFT` inside the plain
   `SequenceGraph`: its delta mass is resolved per combination state as
   (observed_shift − Σ other members), which generalizes the anchored
   pseudo-residue. One graph, one DP, closed and unclosed candidates ranked
   in the same pass.
2. Fixed modifications are positional (TopPIC ProForma positions: "nTerm",
   "cTerm", or 1-based residue), pinned into the residue composition at
   build time, instead of IP's residue-class baking in `AminoAcidSet`
   (they come from a mod file there; ours come from the PrSM row).

Nothing here assigns chemistry: scores rank candidate models; verdicts stay
with the adjudication layer.
*/

use std::collections::{HashMap, HashSet};

use crate::fragment::MassModification;
use crate::localize::{
    parse_msalign_spectrum, query_prsm_localize_spectra_restricted, DeconvolutedPeak,
    MsAlignSpectrum, PrsmLocalizeQuery, PrsmLocalizeResult,
};
use crate::MassSpecError;

// ── integer mass space ─────────────────────────────────────────────────────

/// Nominal-mass binning (IP `Constants.RescalingConstant`).
pub const RESCALING_NOMINAL: f64 = 0.9995;
/// High-precision binning (IP `Constants.RescalingConstantHighPrecision`):
/// bin width ≈ 3.645 mDa, fine enough to separate every distinct
/// residue/modification composition sum.
pub const RESCALING_HIGH_PRECISION: f64 = 274.335_215;

pub fn nominal_bin(mass: f64) -> i64 {
    (mass * RESCALING_NOMINAL).round() as i64
}

pub fn fine_bin(mass: f64) -> i64 {
    (mass * RESCALING_HIGH_PRECISION).round() as i64
}

/// IP `Atom.AtomArr` (UniMod atom masses; compound "atoms" included so glycan
/// modifications are expressible). Index into `ELEMENTS` is the element id.
const ELEMENTS: &[(&str, f64, i32)] = &[
    ("H", 1.007_825_035, 1),
    ("2H", 2.014_101_779, 2),
    ("D", 2.014_101_779, 2),
    ("Li", 7.016_003, 7),
    ("C", 12.0, 12),
    ("13C", 13.003_354_83, 13),
    ("N", 14.003_074, 14),
    ("15N", 15.000_108_97, 15),
    ("O", 15.994_914_63, 16),
    ("18O", 17.999_160_3, 18),
    ("F", 18.998_403_22, 19),
    ("Na", 22.989_767_7, 23),
    ("P", 30.973_762, 31),
    ("S", 31.972_070_7, 32),
    ("Cl", 34.968_852_72, 35),
    ("K", 38.963_707_4, 39),
    ("Ca", 39.962_590_6, 40),
    ("Fe", 55.934_939_3, 56),
    ("Ni", 57.935_346_2, 58),
    ("Cu", 62.929_598_9, 63),
    ("Zn", 63.929_144_8, 64),
    ("Br", 78.918_336_1, 79),
    ("Se", 79.916_519_6, 80),
    ("Mo", 97.905_407_3, 98),
    ("Ag", 106.905_092, 107),
    ("I", 126.904_473, 127),
    ("Au", 196.966_543, 197),
    ("Hg", 201.970_617, 202),
    ("B", 11.009_305_54, 11),
    ("As", 74.921_59, 75),
    ("Mg", 23.985_043, 24),
    ("Hex", 162.052_824, 162),
    ("HexNAc", 203.079_373, 203),
    ("dHex", 146.057_909, 146),
    ("NeuAc", 291.095_417, 291),
    ("NeuGc", 307.090_331, 307),
    ("Hep", 192.063_388, 192),
    ("Pent", 132.042_257, 85),
];

/// Integer atom-count composition (IP `Composition` + its
/// `CompositionWithDeltaMass` extension folded in: an optional free delta
/// mass carries masses no element formula spells).
#[derive(Clone, Debug, PartialEq, Default)]
pub struct GraphComposition {
    pub c: i16,
    pub h: i16,
    pub n: i16,
    pub o: i16,
    pub s: i16,
    pub p: i16,
    /// Additional elements (Br, Se, glycans, metals), sorted by element id.
    pub extra: Vec<(u16, i16)>,
    /// Free delta mass (Da) when the composition carries one.
    pub delta_mass: Option<f64>,
}

impl GraphComposition {
    pub fn zero() -> Self {
        Self::default()
    }

    pub fn from_counts(c: i16, h: i16, n: i16, o: i16, s: i16, p: i16) -> Self {
        Self {
            c,
            h,
            n,
            o,
            s,
            p,
            extra: Vec::new(),
            delta_mass: None,
        }
    }

    pub fn from_delta(mass: f64) -> Self {
        Self {
            delta_mass: Some(mass),
            ..Self::default()
        }
    }

    pub fn is_zero(&self) -> bool {
        self.c == 0
            && self.h == 0
            && self.n == 0
            && self.o == 0
            && self.s == 0
            && self.p == 0
            && self.extra.is_empty()
            && self.delta_mass.is_none_or(|m| m == 0.0)
    }

    /// Count arithmetic addition (IP `Composition.AddComposition`): counts
    /// add element-wise, free deltas add when both carry one.
    pub fn add(&self, other: &Self) -> Self {
        let mut extra = self.extra.clone();
        for &(el, count) in &other.extra {
            if let Some(slot) = extra.iter_mut().find(|(e, _)| *e == el) {
                slot.1 += count;
            } else {
                extra.push((el, count));
            }
        }
        extra.retain(|&(_, count)| count != 0);
        extra.sort_unstable();
        let delta = match (self.delta_mass, other.delta_mass) {
            (Some(a), Some(b)) => Some(a + b),
            (Some(a), None) => Some(a),
            (None, Some(b)) => Some(b),
            (None, None) => None,
        };
        Self {
            c: self.c + other.c,
            h: self.h + other.h,
            n: self.n + other.n,
            o: self.o + other.o,
            s: self.s + other.s,
            p: self.p + other.p,
            extra,
            delta_mass: delta.filter(|m| *m != 0.0),
        }
    }

    pub fn sub(&self, other: &Self) -> Self {
        self.add(&other.negate())
    }

    pub fn negate(&self) -> Self {
        let extra = self.extra.iter().map(|&(el, count)| (el, -count)).collect();
        Self {
            c: -self.c,
            h: -self.h,
            n: -self.n,
            o: -self.o,
            s: -self.s,
            p: -self.p,
            extra,
            delta_mass: self.delta_mass.map(|m| -m),
        }
    }

    /// Monoisotopic mass — the lazy double projection of the integer counts.
    pub fn mass(&self) -> f64 {
        let mut mass = self.c as f64 * 12.0
            + self.h as f64 * 1.007_825_035
            + self.n as f64 * 14.003_074
            + self.o as f64 * 15.994_914_63
            + self.s as f64 * 31.972_070_7
            + self.p as f64 * 30.973_762;
        for &(el, count) in &self.extra {
            mass += ELEMENTS[el as usize].1 * count as f64;
        }
        mass + self.delta_mass.unwrap_or(0.0)
    }

    /// Integer nominal mass (IP `Composition.NominalMass`).
    pub fn nominal_mass(&self) -> i32 {
        let mut nominal = self.c as i32 * 12
            + self.h as i32
            + self.n as i32 * 14
            + self.o as i32 * 16
            + self.s as i32 * 32
            + self.p as i32 * 31;
        for &(el, count) in &self.extra {
            nominal += ELEMENTS[el as usize].2 * count as i32;
        }
        nominal + self.delta_mass.map(|m| nominal_bin(m) as i32).unwrap_or(0)
    }
}

/// IP `AminoAcid.StandardAminoAcidArr` — residue formulas (water already
/// removed): residue composition per one-letter code.
pub fn residue_composition(residue: char) -> Option<GraphComposition> {
    let (c, h, n, o, s) = match residue {
        'G' => (2, 3, 1, 1, 0),
        'A' => (3, 5, 1, 1, 0),
        'S' => (3, 5, 1, 2, 0),
        'P' => (5, 7, 1, 1, 0),
        'V' => (5, 9, 1, 1, 0),
        'T' => (4, 7, 1, 2, 0),
        'C' => (3, 5, 1, 1, 1),
        'L' | 'I' => (6, 11, 1, 1, 0),
        'N' => (4, 6, 2, 2, 0),
        'D' => (4, 5, 1, 3, 0),
        'Q' => (5, 8, 2, 2, 0),
        'K' => (6, 12, 2, 1, 0),
        'E' => (5, 7, 1, 3, 0),
        'M' => (5, 9, 1, 1, 1),
        'H' => (6, 7, 3, 1, 0),
        'F' => (9, 9, 1, 1, 0),
        'R' => (6, 12, 4, 1, 0),
        'Y' => (9, 9, 1, 2, 0),
        'W' => (11, 10, 2, 1, 0),
        _ => return None,
    };
    Some(GraphComposition::from_counts(c, h, n, o, s, 0))
}

// ── modification space ─────────────────────────────────────────────────────

/// IP `SequenceLocation` (with the `AffectedLocations` hierarchy from
/// `AminoAcidSet`'s static constructor: Everywhere fans out to all sites;
/// PeptideNTerm also covers ProteinNTerm, and so on).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SequenceLocation {
    Everywhere,
    PeptideNTerm,
    PeptideCTerm,
    ProteinNTerm,
    ProteinCTerm,
}

impl SequenceLocation {
    pub fn affected(self) -> &'static [SequenceLocation] {
        use SequenceLocation::*;
        match self {
            Everywhere => &[
                Everywhere,
                PeptideNTerm,
                PeptideCTerm,
                ProteinNTerm,
                ProteinCTerm,
            ],
            PeptideNTerm => &[PeptideNTerm, ProteinNTerm],
            PeptideCTerm => &[PeptideCTerm, ProteinCTerm],
            ProteinNTerm => &[ProteinNTerm],
            ProteinCTerm => &[ProteinCTerm],
        }
    }

    pub fn parse(text: &str) -> Option<Self> {
        use SequenceLocation::*;
        match text
            .to_ascii_lowercase()
            .replace(['-', '_', ' '], "")
            .as_str()
        {
            "everywhere" => Some(Everywhere),
            "peptidenterm" | "peptident" => Some(PeptideNTerm),
            "peptidecterm" | "peptidecterminus" | "peptidec" => Some(PeptideCTerm),
            "proteinnterm" | "proteinnterminus" | "proteinn" => Some(ProteinNTerm),
            "proteincterm" | "proteincterminus" | "proteinc" => Some(ProteinCTerm),
            _ => None,
        }
    }
}

/// One searchable modification (IP `SearchModification`): composition +
/// target residues (`'*'` = any residue — our extension so the virtual SHIFT
/// can float; in IP `'*'` only names the terminal marker) + location.
#[derive(Clone, Debug)]
pub struct SearchModification {
    pub name: String,
    pub composition: GraphComposition,
    pub residues: Vec<char>,
    pub location: SequenceLocation,
    /// Optional positional prior (1-based protein positions, "nTerm"/"cTerm"
    /// allowed): the search may place this modification ONLY inside it.
    /// Distinct from fixed mods (pinned, unsearched) — this RESTRICTS the
    /// placement space (TopPIC interval priors, surgical re-search).
    pub positions: Option<Vec<String>>,
}

/// Designated free-mass modification: at most one per combination; its
/// contribution is resolved per combination state (see `ModificationParams`).
pub const SHIFT_NAME: &str = "SHIFT";

/// Precompiled combination space (IP `ModificationParams`): every multiset of
/// library modifications up to `max_mods`, indexed, with the transition table
/// `(prev combination, added mod) -> new combination` and the reverse map
/// used by the traceback.
pub struct ModificationParams {
    modifications: Vec<SearchModification>,
    max_mods: usize,
    /// Sorted multisets of mod indices (the "no mod" member removed).
    combinations: Vec<Vec<u16>>,
    index_of: HashMap<Vec<u16>, u32>,
    transition: HashMap<(u32, u16), u32>,
    reverse: HashMap<(u32, u32), u16>,
    shift_index: Option<u16>,
    /// Observed shift the SHIFT delta resolves against (see
    /// `combination_composition`).
    shift_observed: Option<f64>,
}

/// Cap matching IP's practical library sizes; larger spaces are a caller
/// error, not a slower run.
const MAX_COMBINATIONS: usize = 200_000;

impl ModificationParams {
    pub fn new(
        modifications: Vec<SearchModification>,
        max_mods: usize,
        shift_observed: Option<f64>,
    ) -> Result<Self, MassSpecError> {
        let shift_index = modifications
            .iter()
            .position(|m| m.name == SHIFT_NAME)
            .map(|i| i as u16);
        if max_mods == 0 {
            return Err(MassSpecError::InvalidQuery(
                "max_mods must be at least 1".into(),
            ));
        }
        if modifications.len() > u16::MAX as usize - 1 {
            return Err(MassSpecError::InvalidQuery("too many modifications".into()));
        }
        let mut combinations: Vec<Vec<u16>> = Vec::new();
        let mut current: Vec<u16> = Vec::new();
        enumerate_combinations(
            &(0..modifications.len() as u16).collect::<Vec<_>>(),
            0,
            max_mods,
            &mut current,
            &mut combinations,
            shift_index,
        );
        if combinations.len() > MAX_COMBINATIONS {
            return Err(MassSpecError::InvalidQuery(format!(
                "modification combination space ({} entries) exceeds {} — narrow the library or max_mods",
                combinations.len(),
                MAX_COMBINATIONS
            )));
        }
        let index_of: HashMap<Vec<u16>, u32> = combinations
            .iter()
            .enumerate()
            .map(|(i, comb)| (comb.clone(), i as u32))
            .collect();
        let mut transition = HashMap::new();
        let mut reverse = HashMap::new();
        for (i, comb) in combinations.iter().enumerate() {
            if comb.len() >= max_mods {
                continue;
            }
            for m in 0..modifications.len() as u16 {
                if Some(m) == shift_index && comb.contains(&m) {
                    continue; // SHIFT appears at most once per combination
                }
                let mut next = comb.clone();
                next.push(m);
                next.sort_unstable();
                if let Some(&j) = index_of.get(&next) {
                    transition.insert((i as u32, m), j);
                    reverse.insert((i as u32, j), m);
                }
            }
        }
        Ok(Self {
            modifications,
            max_mods,
            combinations,
            index_of,
            transition,
            reverse,
            shift_index,
            shift_observed,
        })
    }

    pub fn num_modifications(&self) -> usize {
        self.modifications.len()
    }

    pub fn max_mods(&self) -> usize {
        self.max_mods
    }

    pub fn num_combinations(&self) -> usize {
        self.combinations.len()
    }

    pub fn modification(&self, index: u16) -> &SearchModification {
        &self.modifications[index as usize]
    }

    pub fn combination(&self, index: u32) -> &[u16] {
        &self.combinations[index as usize]
    }

    /// The combination's total composition. The designated SHIFT member
    /// resolves dynamically: its contribution is
    /// `shift_observed − Σ(other members)`, so that SHIFT + pinned library
    /// members always sum to the observed shift (the generalization of
    /// `ShiftedSequenceGraph`'s anchored pseudo-residue).
    pub fn combination_composition(&self, index: u32) -> GraphComposition {
        let comb = self.combination(index);
        let mut total = GraphComposition::zero();
        let mut shift_present = false;
        for &m in comb {
            if Some(m) == self.shift_index {
                shift_present = true;
                continue;
            }
            total = total.add(&self.modifications[m as usize].composition);
        }
        if shift_present {
            let observed = self.shift_observed.unwrap_or(0.0);
            let free = observed - total.mass();
            total = total.add(&GraphComposition::from_delta(free));
        }
        total
    }

    pub fn combination_mass(&self, index: u32) -> f64 {
        self.combination_composition(index).mass()
    }

    pub fn transition_index(&self, prev: u32, added: u16) -> Option<u32> {
        self.transition.get(&(prev, added)).copied()
    }

    /// Index of the (sorted-independent) combination key, for tests and
    /// callers that address combinations by member multiset.
    pub fn index_of_key(&self, key: &[u16]) -> Option<u32> {
        let mut sorted = key.to_vec();
        sorted.sort_unstable();
        self.index_of.get(&sorted).copied()
    }

    pub fn modification_between(&self, prev: u32, cur: u32) -> Option<u16> {
        self.reverse.get(&(prev, cur)).copied()
    }

    pub fn shift_index(&self) -> Option<u16> {
        self.shift_index
    }

    /// Number of DISTINCT placements (site arrangements) of combination
    /// `index` — informational (the search itself never enumerates them).
    pub fn combination_label(&self, index: u32) -> String {
        let comb = self.combination(index);
        if comb.is_empty() {
            return "unmodified".to_string();
        }
        comb.iter()
            .map(|&m| self.modifications[m as usize].name.as_str())
            .collect::<Vec<_>>()
            .join("+")
    }
}

fn enumerate_combinations(
    pool: &[u16],
    start: usize,
    max: usize,
    current: &mut Vec<u16>,
    out: &mut Vec<Vec<u16>>,
    shift_index: Option<u16>,
) {
    out.push(current.clone());
    if current.len() >= max {
        return;
    }
    for (k, &m) in pool.iter().enumerate().skip(start) {
        if Some(m) == shift_index && current.contains(&m) {
            continue;
        }
        current.push(m);
        enumerate_combinations(pool, k, max, current, out, shift_index);
        current.pop();
    }
}

// ── the sequence graph ─────────────────────────────────────────────────────

/// One positional fixed modification (TopPIC ProForma position space:
/// "nTerm", "cTerm", or a 1-based residue number).
#[derive(Clone, Debug, serde::Deserialize, serde::Serialize)]
#[serde(rename_all = "camelCase")]
pub struct PositionalModification {
    pub position: String,
    pub delta_mass: f64,
}

#[derive(Clone)]
pub struct Node {
    pub mod_comb_index: u32,
    /// Previous-layer node indices in insertion order (IP tie-breaking on
    /// paths keeps this order: strict `>` while scanning).
    pub prev: Vec<u32>,
    composition: Option<GraphComposition>,
    complementary: Option<GraphComposition>,
}

/// Layered sequence graph (IP `SequenceGraph`, reversed construction:
/// layer 0 is the source, layer 1 the C-terminal marker, residues follow
/// C→N, the last layer is the N-terminal marker).
pub struct SequenceGraph<'a> {
    params: &'a ModificationParams,
    sequence: Vec<char>,
    max_seq_index: usize,
    index: usize,
    num_nterm_cleavages: usize,
    /// Per layer: the residue char actually placed ('[' / ']' are the
    /// protein terminal markers) and its composition including any
    /// positional fixed modification baked in.
    residues: Vec<(char, GraphComposition)>,
    fixed_instances: Vec<Vec<(String, String)>>,
    fixed_positions: HashSet<String>,
    suffix_composition: Vec<GraphComposition>,
    graph: Vec<Vec<Node>>,
    candidates: ModCandidates,
    sink: Option<(u32, GraphComposition)>,
}

/// Residue classes that may carry variable modifications at a layer, per
/// location — built from the library once (IP `AminoAcidSet`'s
/// location-specific maps).
struct ModCandidates {
    by_location: HashMap<(SequenceLocation, char), Vec<u16>>,
}

impl ModCandidates {
    fn build(params: &ModificationParams) -> Self {
        let mut by_location: HashMap<(SequenceLocation, char), Vec<u16>> = HashMap::new();
        for (i, m) in params
            .modifications
            .iter()
            .enumerate()
            .take(params.num_modifications())
        {
            for &loc in m.location.affected() {
                for &res in &m.residues {
                    by_location.entry((loc, res)).or_default().push(i as u16);
                }
            }
        }
        Self { by_location }
    }

    fn get(&self, location: SequenceLocation, residue: char) -> &[u16] {
        static EMPTY: [u16; 0] = [];
        self.by_location
            .get(&(location, residue))
            .map_or(&EMPTY, |v| v.as_slice())
    }
}

impl<'a> SequenceGraph<'a> {
    /// Build the graph over `sequence` (one-letter, no termini) with protein
    /// terminal markers and positional fixed modifications pinned in.
    pub fn new(
        params: &'a ModificationParams,
        sequence: &str,
        fixed: &[PositionalModification],
    ) -> Result<Self, MassSpecError> {
        let chars: Vec<char> = sequence.chars().filter(|c| !c.is_whitespace()).collect();
        if chars.len() < 3 {
            return Err(MassSpecError::InvalidQuery(
                "sequence must contain at least three residues".into(),
            ));
        }
        let len = chars.len();
        let max_seq_index = len + 3;
        let candidates = ModCandidates::build(params);
        let mut graph = Self {
            params,
            sequence: chars.clone(),
            max_seq_index,
            index: 0,
            num_nterm_cleavages: 0,
            residues: vec![('\0', GraphComposition::zero()); max_seq_index],
            fixed_instances: vec![Vec::new(); max_seq_index],
            fixed_positions: fixed
                .iter()
                .map(|item| item.position.trim().to_ascii_lowercase())
                .collect(),
            suffix_composition: Vec::with_capacity(max_seq_index),
            graph: Vec::with_capacity(max_seq_index),
            candidates,
            sink: None,
        };
        graph.suffix_composition.push(GraphComposition::zero());
        graph.graph.push(vec![Node {
            mod_comb_index: 0,
            prev: Vec::new(),
            composition: None,
            complementary: None,
        }]);

        let fixed_at = |position: &str| -> Option<usize> {
            let lower = position.to_ascii_lowercase();
            if lower == "nterm" {
                Some(max_seq_index - 1)
            } else if lower == "cterm" {
                Some(1)
            } else {
                lower.parse::<usize>().ok().map(|p| len + 2 - p)
            }
        };
        let mut pinned: HashMap<usize, Vec<GraphComposition>> = HashMap::new();
        for m in fixed {
            let layer = fixed_at(&m.position).ok_or_else(|| {
                MassSpecError::InvalidQuery(format!(
                    "fixed modification position {:?} is not nTerm/cTerm/1-based",
                    m.position
                ))
            })?;
            if layer == 1 || layer == max_seq_index - 1 {
                // Terminal markers themselves carry the composition.
            } else if !(2..=len + 1).contains(&layer) {
                return Err(MassSpecError::InvalidQuery(format!(
                    "fixed modification position {} outside sequence 1..{}",
                    m.position, len
                )));
            }
            pinned
                .entry(layer)
                .or_default()
                .push(GraphComposition::from_delta(m.delta_mass));
            let name = format!("fixed({:.4})", m.delta_mass);
            graph.fixed_instances[layer].push((name, m.position.clone()));
        }

        graph.put_marker(']', 1, pinned.get(&1).map(|v| v.as_slice()));
        for i in (0..len).rev() {
            let layer = len + 2 - (i + 1);
            if !graph.put_residue(chars[i], layer, pinned.get(&layer).map(|v| v.as_slice()))? {
                return Err(MassSpecError::InvalidQuery(format!(
                    "non-standard residue {:?} in sequence",
                    chars[i]
                )));
            }
        }
        graph.put_marker(
            '[',
            max_seq_index - 1,
            pinned.get(&(max_seq_index - 1)).map(|v| v.as_slice()),
        );
        Ok(graph)
    }

    pub fn is_valid(&self) -> bool {
        !self.graph.is_empty()
    }

    pub fn index(&self) -> usize {
        self.index
    }

    fn location_for(&self, layer: usize, residue: char) -> SequenceLocation {
        use SequenceLocation::*;
        let len = self.sequence.len();
        let n_marker_slot = self.max_seq_index - 1 - self.num_nterm_cleavages;
        if layer == 1 {
            return if residue == ')' {
                PeptideCTerm
            } else {
                ProteinCTerm
            };
        }
        if layer == n_marker_slot {
            return if residue == '(' {
                PeptideNTerm
            } else {
                ProteinNTerm
            };
        }
        if layer == 2 {
            return if self.residues[1].0 == ')' {
                PeptideCTerm
            } else {
                ProteinCTerm
            };
        }
        if layer == n_marker_slot.saturating_sub(1) {
            return if self.residues[n_marker_slot].0 == '(' {
                PeptideNTerm
            } else {
                ProteinNTerm
            };
        }
        let _ = len;
        Everywhere
    }

    /// IP `PutAminoAcid`: resolve the location, accumulate the suffix
    /// composition, then expand the layer — each previous node keeps an
    /// unmodified edge and gains one edge per applicable modification;
    /// nodes are deduplicated BY combination index.
    fn put_residue(
        &mut self,
        residue: char,
        layer: usize,
        pinned: Option<&[GraphComposition]>,
    ) -> Result<bool, MassSpecError> {
        let base = residue_composition(residue)
            .ok_or_else(|| MassSpecError::InvalidQuery(format!("unknown residue {residue}")))?;
        let mut composition = base;
        if let Some(list) = pinned {
            for extra in list {
                composition = composition.add(extra);
            }
        }
        self.put(residue, composition, layer, true)
    }

    fn put_marker(&mut self, marker: char, layer: usize, pinned: Option<&[GraphComposition]>) {
        let mut composition = GraphComposition::zero();
        if let Some(list) = pinned {
            for extra in list {
                composition = composition.add(extra);
            }
        }
        // Marker placement cannot fail (zero composition).
        let _ = self.put(marker, composition, layer, true);
    }

    fn put(
        &mut self,
        residue: char,
        composition: GraphComposition,
        layer: usize,
        _is_residue: bool,
    ) -> Result<bool, MassSpecError> {
        self.index = layer;
        self.residues[layer] = (residue, composition.clone());
        let suffix = self.suffix_composition[layer - 1].add(&composition);
        while self.suffix_composition.len() <= layer {
            self.suffix_composition.push(GraphComposition::zero());
        }
        self.suffix_composition[layer] = suffix;

        let location = self.location_for(layer, residue);
        let mut mod_indices: Vec<u16> =
            if residue == '[' || residue == ']' || residue == '(' || residue == ')' {
                self.candidates.get(location, residue).to_vec()
            } else {
                let mut merged: Vec<u16> = self
                    .candidates
                    .get(location, residue)
                    .iter()
                    .copied()
                    .chain(self.candidates.get(location, '*').iter().copied())
                    .collect();
                merged.sort_unstable();
                merged.dedup();
                merged
            };
        // positional priors: a mod with `positions` may only occupy the
        // listed protein positions (or nTerm/cTerm on the marker layers)
        if !mod_indices.is_empty() {
            let len = self.sequence.len();
            let position_label: String = if residue == '[' {
                "nTerm".to_string()
            } else if residue == ']' {
                "cTerm".to_string()
            } else {
                (len + 2 - layer).to_string()
            };
            let fixed_here = self.fixed_positions.contains(&position_label);
            let fixed_nterm = self.fixed_positions.contains("nterm");
            let fixed_cterm = self.fixed_positions.contains("cterm");
            mod_indices.retain(|&m| {
                let modification = self.params.modification(m);
                let occupied = fixed_here
                    || (fixed_nterm
                        && (residue == '['
                            || matches!(
                                modification.location,
                                SequenceLocation::PeptideNTerm | SequenceLocation::ProteinNTerm
                            )))
                    || (fixed_cterm
                        && (residue == ']'
                            || matches!(
                                modification.location,
                                SequenceLocation::PeptideCTerm | SequenceLocation::ProteinCTerm
                            )));
                !occupied
                    && modification
                        .positions
                        .as_ref()
                        .is_none_or(|allowed| allowed.contains(&position_label))
            });
        }

        while self.graph.len() <= layer {
            self.graph.push(Vec::new());
        }
        let prev_layer = self.graph[layer - 1].clone();
        if mod_indices.is_empty() {
            self.graph[layer] = prev_layer
                .iter()
                .map(|n| Node {
                    mod_comb_index: n.mod_comb_index,
                    prev: vec![],
                    composition: None,
                    complementary: None,
                })
                .collect();
            for (i, node) in self.graph[layer].iter_mut().enumerate() {
                node.prev.push(i as u32);
            }
        } else {
            let mut by_comb: HashMap<u32, Node> = HashMap::new();
            for (i, prev_node) in prev_layer.iter().enumerate() {
                let prev_comb = prev_node.mod_comb_index;
                // unmodified edge
                let entry = by_comb.entry(prev_comb).or_insert_with(|| Node {
                    mod_comb_index: prev_comb,
                    prev: Vec::new(),
                    composition: None,
                    complementary: None,
                });
                entry.prev.push(i as u32);
                // modified edges
                for &m in mod_indices.iter() {
                    if let Some(new_comb) = self.params.transition_index(prev_comb, m) {
                        let entry = by_comb.entry(new_comb).or_insert_with(|| Node {
                            mod_comb_index: new_comb,
                            prev: Vec::new(),
                            composition: None,
                            complementary: None,
                        });
                        entry.prev.push(i as u32);
                    }
                }
            }
            let mut nodes: Vec<Node> = by_comb.into_values().collect();
            // Deterministic layer order (the combination table order).
            nodes.sort_unstable_by_key(|n| n.mod_comb_index);
            self.graph[layer] = nodes;
        }
        Ok(true)
    }

    /// Layer-resolved suffix composition + combination composition,
    /// cached on the node (IP `GetComposition`).
    pub fn node_composition(&mut self, layer: usize, node: usize) -> GraphComposition {
        let comb = self.graph[layer][node].mod_comb_index;
        if let Some(cached) = &self.graph[layer][node].composition {
            return cached.clone();
        }
        let composition =
            self.suffix_composition[layer].add(&self.params.combination_composition(comb));
        self.graph[layer][node].composition = Some(composition.clone());
        composition
    }

    /// sink composition − node composition (IP `GetComplementaryComposition`);
    /// cleared by `set_sink` (IP `SetSink` clears the whole cache).
    pub fn complementary_composition(
        &mut self,
        layer: usize,
        node: usize,
    ) -> Option<GraphComposition> {
        // clone the sink composition out first: node_composition needs &mut
        let sink = self.sink.as_ref()?.1.clone();
        if let Some(cached) = &self.graph[layer][node].complementary {
            return Some(cached.clone());
        }
        let node_comp = self.node_composition(layer, node);
        let comp = sink.sub(&node_comp);
        self.graph[layer][node].complementary = Some(comp.clone());
        Some(comp)
    }

    /// Sink-layer nodes and their combination indices — the candidate
    /// proteoform set (IP `GetNumProteoformCompositions` /
    /// `GetModificationCombinations`).
    /// (node index, combination index) pairs of the sink layer.
    pub fn top_nodes(&self) -> Vec<(u32, u32)> {
        self.graph[self.index]
            .iter()
            .enumerate()
            .map(|(i, n)| (i as u32, n.mod_comb_index))
            .collect()
    }

    /// Fix the sink: the chosen modification combination whose best
    /// arrangement the DP will score.
    pub fn set_sink(&mut self, node: usize) -> Result<(), MassSpecError> {
        if node >= self.graph[self.index].len() {
            return Err(MassSpecError::InvalidQuery(format!(
                "sink node {node} outside layer {}",
                self.index
            )));
        }
        let composition = self.node_composition(self.index, node);
        for layer in self.graph.iter_mut() {
            for n in layer.iter_mut() {
                n.complementary = None;
            }
        }
        self.sink = Some((node as u32, composition));
        Ok(())
    }

    pub fn sink_composition_with_h2o(&self) -> Option<GraphComposition> {
        self.sink
            .as_ref()
            .map(|(_, c)| c.add(&GraphComposition::from_counts(0, 2, 0, 1, 0, 0)))
    }

    /// The one-pass DP (IP `GetFragmentScore`): forward max-DP over layers,
    /// node score from the scorer seeing BOTH sides of the cleavage
    /// (prefix = complementary composition, suffix = node composition).
    /// Returns (score, per-layer best-prev) for traceback reuse.
    fn run_dp(
        &mut self,
        scorer: &mut dyn GraphScorer,
    ) -> (f64, Vec<Vec<Option<u32>>>, Vec<Vec<Option<f64>>>) {
        let sink_node = self.sink.as_ref().expect("set_sink before scoring").0 as usize;
        let mut node_score: Vec<Vec<Option<f64>>> =
            self.graph.iter().map(|l| vec![None; l.len()]).collect();
        let mut max_score: Vec<Vec<Option<f64>>> =
            self.graph.iter().map(|l| vec![None; l.len()]).collect();
        let mut best_prev: Vec<Vec<Option<u32>>> =
            self.graph.iter().map(|l| vec![None; l.len()]).collect();
        max_score[0][0] = Some(0.0);
        for layer in 1..=self.index {
            let nodes = if layer == self.index {
                vec![sink_node]
            } else {
                (0..self.graph[layer].len()).collect()
            };
            for node in nodes {
                // sink node contributes no cleavage of its own (IP pre-seeds
                // nodeScore[sink] = 0)
                let score = if layer == self.index {
                    0.0
                } else {
                    let suffix = self.node_composition(layer, node);
                    let Some(prefix) = self.complementary_composition(layer, node) else {
                        continue;
                    };
                    let s = scorer.fragment_score(&prefix, &suffix);
                    node_score[layer][node] = Some(s);
                    s
                };
                let mut best: Option<(f64, u32)> = None;
                for &prev in &self.graph[layer][node].prev {
                    if let Some(ps) = max_score[layer - 1][prev as usize] {
                        // strict `>` keeps the first prev on ties (IP order)
                        if best.is_none_or(|(bs, _)| ps > bs) {
                            best = Some((ps, prev));
                        }
                    }
                }
                let total = score + best.map_or(0.0, |(b, _)| b);
                max_score[layer][node] = Some(total);
                best_prev[layer][node] = best.map(|(_, p)| p);
            }
        }
        (
            max_score[self.index][sink_node].unwrap_or(0.0),
            best_prev,
            max_score,
        )
    }

    pub fn get_fragment_score(&mut self, scorer: &mut dyn GraphScorer) -> f64 {
        self.run_dp(scorer).0
    }

    /// All OPTIMAL placements (IP tie semantics made explicit): every path
    /// achieving the sink's max score, as modification lists. A prev node is
    /// on an optimal path iff its max score equals the best prev max — the
    /// same criterion the DP's strict-`>` scan already encodes. Capped at
    /// `cap`; `capped` flags that more exist. Mods come out in ascending
    /// protein position order (the walk descends N->C).
    pub fn optimal_placements(
        &mut self,
        scorer: &mut dyn GraphScorer,
        cap: usize,
    ) -> (f64, Vec<Vec<(String, String)>>, bool) {
        const EPS: f64 = 1e-9;
        let (score, _best_prev, max_score) = self.run_dp(scorer);
        let len = self.sequence.len();
        let sink_node = self.sink.as_ref().expect("set_sink first").0 as usize;
        let label = |residue: char, layer: usize| -> String {
            if residue == '[' {
                "nTerm".to_string()
            } else if residue == ']' {
                "cTerm".to_string()
            } else {
                (len + 2 - layer).to_string()
            }
        };
        let mut out: Vec<Vec<(String, String)>> = Vec::new();
        let mut capped = false;
        let mut stack: Vec<(usize, usize, Vec<(String, String)>)> =
            vec![(self.index, sink_node, Vec::new())];
        while let Some((layer, node, mods)) = stack.pop() {
            if layer == 0 {
                out.push(mods);
                if out.len() >= cap {
                    capped = true;
                    break;
                }
                continue;
            }
            let prevs = &self.graph[layer][node].prev;
            let best_prev_max = prevs
                .iter()
                .filter_map(|p| max_score[layer - 1][*p as usize])
                .fold(f64::NEG_INFINITY, f64::max);
            for &prev in prevs {
                let Some(ps) = max_score[layer - 1][prev as usize] else {
                    continue;
                };
                if (ps - best_prev_max).abs() > EPS {
                    continue;
                }
                let mut next_mods = mods.clone();
                let (residue, _) = self.residues[layer];
                let cur_comb = self.graph[layer][node].mod_comb_index;
                let prev_comb = self.graph[layer - 1][prev as usize].mod_comb_index;
                if prev_comb != cur_comb {
                    if let Some(m) = self.params.modification_between(prev_comb, cur_comb) {
                        next_mods.push((
                            self.params.modification(m).name.clone(),
                            label(residue, layer),
                        ));
                    }
                }
                stack.push((layer - 1, prev as usize, next_mods));
            }
        }
        (score, out, capped)
    }

    /// DP + traceback (IP `GetFragmentScoreAndModifications`): reconstructs
    /// which modification sits at which position along the best path.
    /// Positions are TopPIC-style: 1-based residue numbers, or
    /// "nTerm"/"cTerm" on the marker layers.
    pub fn get_fragment_score_and_modifications(
        &mut self,
        scorer: &mut dyn GraphScorer,
    ) -> (f64, Vec<(String, String)>) {
        let (score, best_prev, _max_score) = self.run_dp(scorer);
        let len = self.sequence.len();
        let mut mods: Vec<(String, String)> = Vec::new();
        // fixed instances along the walk, then dynamic ones
        let mut layer = self.index;
        let mut node = self.sink.as_ref().expect("set_sink first").0 as usize;
        let mut fixed_seen: Vec<(String, String)> = Vec::new();
        while layer > 0 {
            let (residue, _) = self.residues[layer];
            let position: String = if residue == '[' {
                "nTerm".to_string()
            } else if residue == ']' {
                "cTerm".to_string()
            } else {
                (len + 2 - layer).to_string()
            };
            fixed_seen.extend(self.fixed_instances[layer].iter().cloned());
            let Some(prev) = best_prev[layer][node] else {
                break;
            };
            let cur_comb = self.graph[layer][node].mod_comb_index;
            let prev_comb = self.graph[layer - 1][prev as usize].mod_comb_index;
            if prev_comb != cur_comb {
                if let Some(m) = self.params.modification_between(prev_comb, cur_comb) {
                    mods.push((self.params.modification(m).name.clone(), position));
                }
            }
            node = prev as usize;
            layer -= 1;
        }
        fixed_seen.extend(mods);
        (score, fixed_seen)
    }

    /// Number of distinct placements leading to each top-layer node
    /// (IP `GetNumProteoformSequences`' count table, computed on demand).
    pub fn num_placements(&self) -> Vec<f64> {
        let mut counts: Vec<Vec<f64>> = self.graph.iter().map(|l| vec![0.0; l.len()]).collect();
        counts[0][0] = 1.0;
        for layer in 1..=self.index {
            for node in 0..self.graph[layer].len() {
                let mut total = 0.0;
                for prev in &self.graph[layer][node].prev {
                    total += counts[layer - 1][*prev as usize];
                }
                counts[layer][node] = total;
            }
        }
        counts[self.index].clone()
    }

    /// N-terminal truncation (IP `CleaveNTerm`): drop the leading residue,
    /// re-purpose the two top layers, keep the lower layers untouched (their
    /// suffix compositions are C-side accumulations and unchanged).
    pub fn cleave_n_term(&mut self) -> Result<(), MassSpecError> {
        if self.index < 4 {
            return Err(MassSpecError::InvalidQuery(
                "sequence too short to cleave".into(),
            ));
        }
        self.index -= 3;
        for layer in self.index + 1..=self.index + 2 {
            for node in self.graph[layer].iter_mut() {
                node.composition = None;
                node.complementary = None;
            }
        }
        self.num_nterm_cleavages += 1;
        let marker_slot = self.max_seq_index - 1 - self.num_nterm_cleavages;
        self.residues[marker_slot] = ('[', GraphComposition::zero());
        let next = self.sequence[self.num_nterm_cleavages];
        // layer numbers must be captured BEFORE put_* advances self.index
        let residue_layer = self.index + 1;
        let marker_layer = self.index + 2;
        self.put_residue(next, residue_layer, None)?;
        self.put_marker('[', marker_layer, None);
        Ok(())
    }

    pub fn num_nterm_cleavages(&self) -> usize {
        self.num_nterm_cleavages
    }
}

// ── scoring ────────────────────────────────────────────────────────────────

/// Edge cost of the DP (IP `IScorer`): one call per (position, combination)
/// node, seeing both sides of the cleavage.
pub trait GraphScorer {
    fn fragment_score(&mut self, prefix: &GraphComposition, suffix: &GraphComposition) -> f64;
}

pub type BaseIon = crate::fragment::IonType;

impl BaseIon {
    pub fn parse_set(text: &str) -> Result<Vec<Self>, MassSpecError> {
        let ions: Result<Vec<_>, _> = text
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(|s| match s.to_ascii_uppercase().as_str() {
                "B" => Ok(Self::B),
                "C" => Ok(Self::C),
                "Y" => Ok(Self::Y),
                "Z" => Ok(Self::Z),
                _ => Err(MassSpecError::InvalidQuery(format!(
                    "unknown graph ion {s:?}; expected B/C/Y/Z"
                ))),
            })
            .collect();
        let ions = ions?;
        if ions.is_empty() {
            return Err(MassSpecError::InvalidQuery("empty ion set".into()));
        }
        Ok(ions)
    }

    /// DP uses residue-sum coordinates. Project offsets once from the same
    /// mzannotate-generated ions as spectral matching, not a second mass table.
    pub fn offset(self) -> (bool, f64) {
        static OFFSETS: std::sync::LazyLock<std::collections::BTreeMap<BaseIon, f64>> =
            std::sync::LazyLock::new(|| {
                crate::theoretical_fragment_list(
                    "AA",
                    &[],
                    &[
                        BaseIon::A,
                        BaseIon::B,
                        BaseIon::C,
                        BaseIon::X,
                        BaseIon::Y,
                        BaseIon::Z,
                    ],
                    1,
                    &[],
                    false,
                )
                .expect("mzannotate supports the canonical AA ion ladder")
                .into_iter()
                .map(|fragment| {
                    (
                        fragment.ion,
                        fragment.neutral_mass - crate::fragment::residue_mass('A').unwrap(),
                    )
                })
                .collect()
            });
        (matches!(self, Self::A | Self::B | Self::C), OFFSETS[&self])
    }
}

/// Counting scorer over the fine integer bins (IP
/// `ProductScorerBasedOnDeconvolutedSpectra.DeconvScorer`): deconvoluted
/// neutral masses are smeared over every fine bin their ± tolerance covers;
/// a fragment scores +1 per ion type whose bin is present.
pub struct BinCountScorer {
    bins: HashSet<i64>,
    ions: Vec<BaseIon>,
}

impl BinCountScorer {
    pub fn new(
        peaks: &[DeconvolutedPeak],
        ions: Vec<BaseIon>,
        tolerance_ppm: f64,
        mass_min: f64,
        mass_max: f64,
    ) -> Self {
        let mut bins = HashSet::new();
        for &peak in peaks {
            let mass = peak.mass;
            if mass < mass_min || mass > mass_max {
                continue;
            }
            let tol = mass * tolerance_ppm * 1e-6;
            let lo = fine_bin(mass - tol);
            let hi = fine_bin(mass + tol);
            for b in lo..=hi {
                bins.insert(b);
            }
        }
        Self { bins, ions }
    }
}

impl GraphScorer for BinCountScorer {
    fn fragment_score(&mut self, prefix: &GraphComposition, suffix: &GraphComposition) -> f64 {
        let mut score = 0.0;
        for &ion in &self.ions {
            let (is_prefix, offset) = ion.offset();
            let mass = if is_prefix {
                prefix.mass()
            } else {
                suffix.mass()
            } + offset;
            if self.bins.contains(&fine_bin(mass)) {
                score += 1.0;
            }
        }
        score
    }
}

/// IP `CompositeScorer`'s trained linear discriminant (default weight set,
/// `WeightScaleFactor = 4.0`) adapted to TopFD's deconvoluted peak contract.
/// Count, intensity and mass error retain their IP meaning. TopFD's generic
/// deconvolution score supplies the 0..1 peak-quality term; it is not exposed
/// as IP's Pearson Corr. IP Dist has no msalign field and is omitted.
pub struct WeightedBinScorer {
    map: HashMap<i64, (f64, f32, f32)>, // fine bin -> (mass, intensity, TopFD quality) of best peak
    reference_intensity: f64,
    ions: Vec<BaseIon>,
}

pub struct IonScoreWeight {
    pub count: f64,
    pub consecutive_match: f64,
    pub intensity: f64,
    /// IP Corr coefficient applied to TopFD's 0..1 deconvolution-quality
    /// proxy. The coefficient is retained, but the observable is never
    /// reported as a Pearson correlation.
    pub deconvolution_quality: f64,
    pub mass_error: f64,
}

pub struct ScoreWeights {
    pub beta0: f64,
    pub cutoff: f64,
    pub complementary_ion_count: f64,
    pub prefix: IonScoreWeight,
    pub suffix: IonScoreWeight,
}

/// IP `CompositeScorer` static constructor, default (non-flip) set,
/// `× WeightScaleFactor (4.0)` applied.
pub const fn ip_default_weights() -> ScoreWeights {
    let s = 4.0;
    ScoreWeights {
        beta0: -5.393_357_744_174_74 * s,
        cutoff: 3.0 * s,
        complementary_ion_count: -0.378_768_218_613_132 * s,
        // Dist (prefix 0.180, suffix -0.284) has no msalign column and stays
        // a documented gap rather than fabricated evidence.
        prefix: IonScoreWeight {
            count: -1.125_486_077 * s,
            consecutive_match: 0.125_523_007 * s,
            deconvolution_quality: 1.473_096_149 * s,
            intensity: 0.323_923_971 * s,
            mass_error: -0.032_959_575 * s,
        },
        suffix: IonScoreWeight {
            count: -0.886_170_429 * s,
            consecutive_match: 0.083_704_707 * s,
            deconvolution_quality: 1.189_293_38 * s,
            intensity: 0.372_673_677 * s,
            mass_error: -0.026_141_488 * s,
        },
    }
}

/// IP `CompositeScorer.GetProbability`: logistic on `Beta0 + score` — the
/// per-match probability IcTopDownLauncher reports alongside the score.
/// `Cutoff` (also ported) is the launcher's reporting filter
/// (`match.Score > Cutoff`).
/// IP `CompositeScorer.GetRefIntensity`: the 87.5th-percentile intensity of
/// the (sorted) peak list — `intensities[(int)(0.8750 * n)]` — NOT the max.
pub fn get_ref_intensity(mut intensities: Vec<f64>) -> f64 {
    if intensities.is_empty() {
        return 0.0;
    }
    intensities.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    intensities[(0.8750 * intensities.len() as f64) as usize]
        .min(intensities[intensities.len() - 1])
}

pub fn get_probability(score: f64) -> f64 {
    let eta = (ip_default_weights().beta0 + score).exp();
    eta / (eta + 1.0)
}

impl WeightedBinScorer {
    pub fn new(
        peaks: &[DeconvolutedPeak],
        ions: Vec<BaseIon>,
        tolerance_ppm: f64,
        mass_min: f64,
        mass_max: f64,
    ) -> Self {
        // IP's `UpdateDeconvPeak` keeps the most intense peak per bin. The
        // msalign fourth column is TopFD's generic deconvolution score: a
        // useful 0..1 quality proxy, but not IP's Pearson Corr.
        let mut map: HashMap<i64, (f64, f32, f32)> = HashMap::new();
        let mut intensities: Vec<f64> = Vec::new();
        for peak in peaks {
            let (mass, intensity) = (peak.mass, peak.intensity);
            if mass < mass_min || mass > mass_max {
                continue;
            }
            intensities.push(f64::from(intensity));
            let topfd_quality = peak
                .score
                .filter(|score| score.is_finite())
                .unwrap_or(0.0)
                .clamp(0.0, 1.0);
            let tol = mass * tolerance_ppm * 1e-6;
            let lo = fine_bin(mass - tol);
            let hi = fine_bin(mass + tol);
            for b in lo..=hi {
                match map.get(&b) {
                    Some(&(_, best_i, _)) if best_i >= intensity => {}
                    _ => {
                        map.insert(b, (mass, intensity, topfd_quality));
                    }
                }
            }
        }
        Self {
            map,
            reference_intensity: get_ref_intensity(intensities).max(f64::MIN_POSITIVE),
            ions,
        }
    }
}

impl GraphScorer for WeightedBinScorer {
    fn fragment_score(&mut self, prefix: &GraphComposition, suffix: &GraphComposition) -> f64 {
        let w = ip_default_weights();
        let mut score = 0.0;
        let mut prefix_hit = false;
        let mut suffix_hit = false;
        for &ion in &self.ions {
            let (is_prefix, offset) = ion.offset();
            let (comp, param) = if is_prefix {
                (prefix, &w.prefix)
            } else {
                (suffix, &w.suffix)
            };
            let mass = comp.mass() + offset;
            if let Some(&(observed, intensity, topfd_quality)) = self.map.get(&fine_bin(mass)) {
                let mass_error_ppm = ((observed - mass).abs() / mass) * 1e6;
                score += param.count;
                score +=
                    param.intensity * (f64::from(intensity) / self.reference_intensity).min(1.0);
                // TopFD quality occupies IP's peak-quality slot without being
                // relabeled as Pearson Corr. Dist is the one unavailable IP
                // term and contributes zero.
                score += param.deconvolution_quality * f64::from(topfd_quality);
                score += param.mass_error * mass_error_ppm;
                if is_prefix {
                    prefix_hit = true;
                } else {
                    suffix_hit = true;
                }
            }
        }
        if prefix_hit && suffix_hit {
            score += w.complementary_ion_count;
        }
        score
    }
}

// ── spectral E-value: mass-bin scoring graph + generating function ────────

/// f64 top-bits truncation bin (IP `MzComparerWithBinning`): keeping the
/// top `num_bits` bits of the IEEE double makes bin equality a ppm-scale
/// mass comparison with no float arithmetic (27 bits ~ 16 ppm half-width).
pub fn truncated_bin(mass: f64, num_bits: u32) -> i64 {
    (mass.to_bits() >> (64 - num_bits)) as i64
}

/// IP `FilteredProteinMassBinning`: fine bins (`x274.335215`) restricted to
/// masses REACHABLE as residue sums (forward closure from bin 0), then
/// compressed onto dense 27-bit truncation bins. Unreachable masses have no
/// bin at all.
pub struct MassBinning {
    coarse: Vec<i64>,     // sorted distinct truncation bins (dense index -> coarse)
    dense_mass: Vec<f64>, // representative mass per dense bin
    num_bits: u32,
}

impl MassBinning {
    pub fn build(residue_masses: &[f64], mass_max: f64, num_bits: u32) -> Self {
        let max_fine = fine_bin(mass_max).max(0) as usize;
        let mut reachable = vec![false; max_fine + 1];
        if !reachable.is_empty() {
            reachable[0] = true;
        }
        for b in 0..=max_fine {
            if !reachable[b] {
                continue;
            }
            let m = b as f64 / RESCALING_HIGH_PRECISION;
            for &r in residue_masses {
                let j = fine_bin(m + r);
                if j >= 0 && (j as usize) <= max_fine {
                    reachable[j as usize] = true;
                }
            }
        }
        use std::collections::BTreeMap;
        let mut coarse_sums: BTreeMap<i64, (f64, usize)> = BTreeMap::new();
        for b in 0..=max_fine {
            if !reachable[b] {
                continue;
            }
            let m = b as f64 / RESCALING_HIGH_PRECISION;
            let c = truncated_bin(m, num_bits);
            let entry = coarse_sums.entry(c).or_insert((0.0, 0));
            entry.0 += m;
            entry.1 += 1;
        }
        let coarse: Vec<i64> = coarse_sums.keys().copied().collect();
        let dense_mass: Vec<f64> = coarse_sums
            .values()
            .map(|(sum, n)| sum / *n as f64)
            .collect();
        Self {
            coarse,
            dense_mass,
            num_bits,
        }
    }

    pub fn num_bins(&self) -> usize {
        self.coarse.len()
    }

    /// Dense bin index of a mass, or None when no reachable bin covers it.
    pub fn bin_number(&self, mass: f64) -> Option<usize> {
        if mass <= 0.0 {
            return Some(0);
        }
        let c = truncated_bin(mass, self.num_bits);
        self.coarse.binary_search(&c).ok()
    }

    pub fn bin_mass(&self, index: usize) -> f64 {
        self.dense_mass[index]
    }
}

/// IP `ProteinScoringGraph` + factory: nodes are dense mass bins; adjacency
/// `j <- i` when one residue mass spans bin i to bin j; node score is the
/// matched-ion score at that prefix mass (prefix + suffix), edge score a
/// consecutive-match bonus. Sequence-independent — the null space over
/// which the generating function integrates.
pub struct MassBinScoringGraph {
    pub adj: Vec<Vec<usize>>, // adj[j] lists prev bins
    pub prefix_scores: Vec<Option<f64>>,
    pub suffix_scores: Vec<Option<f64>>,
}

impl MassBinScoringGraph {
    pub fn build(
        binning: &MassBinning,
        residue_masses: &[f64],
        prefix_scores: Vec<Option<f64>>,
        suffix_scores: Vec<Option<f64>>,
    ) -> Self {
        let n = binning.num_bins();
        let mut adj = vec![Vec::new(); n];
        for i in 0..n {
            let m = binning.bin_mass(i);
            for &r in residue_masses {
                if let Some(j) = binning.bin_number(m + r) {
                    adj[j].push(i);
                }
            }
        }
        Self {
            adj,
            prefix_scores,
            suffix_scores,
        }
    }

    fn node_score(&self, node: usize) -> f64 {
        self.prefix_scores[node].unwrap_or(0.0) + self.suffix_scores[node].unwrap_or(0.0)
    }

    fn edge_score(&self, prev: usize, node: usize) -> f64 {
        // IP CompositeScorer consecutive-match terms (x WeightScaleFactor)
        let w = ip_default_weights();
        let mut score = 0.0;
        if self.prefix_scores[prev].is_some() && self.prefix_scores[node].is_some() {
            score += w.prefix.consecutive_match;
        }
        if self.suffix_scores[prev].is_some() && self.suffix_scores[node].is_some() {
            score += w.suffix.consecutive_match;
        }
        score
    }
}

/// IP `ScoringGraphEdge.Weight` (constant): the per-edge probability
/// discount every convolution multiplies — the amino-acid prior IP bakes
/// into spectral E-values.
const GF_EDGE_WEIGHT: f64 = 0.05;

/// IP `ScoreDistribution`: unnormalized spectral probability mass over
/// integer scores (path counts through the mass-bin graph).
#[derive(Clone, Debug)]
pub struct ScoreDistribution {
    min_score: i32,
    probs: Vec<f64>,
}

impl ScoreDistribution {
    fn new(min_score: i32, max_score: i32) -> Self {
        Self {
            min_score,
            probs: vec![0.0; (max_score - min_score + 1).max(1) as usize],
        }
    }

    fn get(&self, score: i32) -> f64 {
        if score < self.min_score || score >= self.min_score + self.probs.len() as i32 {
            return 0.0;
        }
        self.probs[(score - self.min_score) as usize]
    }

    fn add_shifted(&mut self, prev: &ScoreDistribution, shift: i32, weight: f64) {
        for (k, &p) in prev.probs.iter().enumerate() {
            if p == 0.0 {
                continue;
            }
            let target = prev.min_score + k as i32 + shift;
            if target < self.min_score {
                continue;
            }
            let idx = (target - self.min_score) as usize;
            if idx < self.probs.len() {
                self.probs[idx] += p * weight;
            }
        }
    }

    fn max_score(&self) -> i32 {
        self.min_score + self.probs.len() as i32 - 1
    }
}

/// IP `GeneratingFunction`: one DP pass over the mass-bin graph produces the
/// score distribution of ALL residue sequences consistent with the mass;
/// `spectral_e_value` is the expected number achieving a score >= observed.
pub struct GeneratingFunction<'g> {
    graph: &'g MassBinScoringGraph,
    distribution: Option<ScoreDistribution>,
}

impl<'g> GeneratingFunction<'g> {
    pub fn new(graph: &'g MassBinScoringGraph) -> Self {
        Self {
            graph,
            distribution: None,
        }
    }

    pub fn compute(&mut self) {
        if self.distribution.is_some() {
            return;
        }
        let n = self.graph.adj.len();
        let mut table: Vec<Option<ScoreDistribution>> = vec![None; n];
        let mut source = ScoreDistribution::new(0, 0);
        source.probs[0] = 1.0;
        if n > 0 {
            table[0] = Some(source);
        }
        for node in 1..n {
            let own = self.graph.node_score(node);
            // IP tracks the min/max bounds from the UNROUNDED combined
            // score and floors/ceils them at construction (`GetScoreDistribution`).
            let mut min_f = f64::INFINITY;
            let mut max_f = f64::NEG_INFINITY;
            let mut has_edge = false;
            for &prev in &self.graph.adj[node] {
                let Some(prev_dist) = &table[prev] else {
                    continue;
                };
                let combined = own + self.graph.edge_score(prev, node);
                min_f = min_f.min(prev_dist.min_score as f64 + combined);
                max_f = max_f.max(prev_dist.max_score() as f64 + combined);
                has_edge = true;
            }
            if !has_edge {
                continue;
            }
            let mut dist = ScoreDistribution::new(min_f.floor() as i32, max_f.ceil() as i32);
            for &prev in &self.graph.adj[node] {
                let Some(prev_dist) = &table[prev] else {
                    continue;
                };
                let combined = own + self.graph.edge_score(prev, node);
                dist.add_shifted(prev_dist, combined.round() as i32, GF_EDGE_WEIGHT);
            }
            table[node] = Some(dist);
        }
        self.distribution = table[n - 1].clone();
    }

    pub fn spectral_e_value(&mut self, score: f64) -> Option<f64> {
        self.compute();
        let dist = self.distribution.as_ref()?;
        let s = score.round() as i32;
        let mut tail = 0.0;
        let lo = s.max(dist.min_score);
        for k in lo..=dist.max_score() {
            tail += dist.get(k);
        }
        Some(if tail <= 0.0 { f64::MIN_POSITIVE } else { tail })
    }
}

// ── top-level search over a msalign file ───────────────────────────────────

#[derive(Clone, Debug, serde::Deserialize, serde::Serialize)]
#[serde(rename_all = "camelCase")]
pub struct LibrarySpec {
    pub name: String,
    pub delta_mass: f64,
    /// Target residues; '*' (or missing) = any residue.
    #[serde(default)]
    pub residues: Vec<char>,
    /// "everywhere" (default) / "proteinNTerm" / "peptideNTerm" / ...
    #[serde(default)]
    pub location: Option<String>,
    /// Optional positional prior: 1-based protein positions (or
    /// "nTerm"/"cTerm") where this modification may be placed.
    #[serde(default)]
    pub positions: Option<Vec<String>>,
}

#[derive(Clone, Debug, serde::Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct GraphSearchQuery {
    /// TopFD `_ms2.msalign` file (read-only).
    pub file: String,
    /// TopPIC ProFormA (e.g. `[Acetyl]-SGRGK...VVYALK)[+121.9569]RQ...`).
    /// When given, `sequence`, `observed_shift` and `fixed_modifications`
    /// are DERIVED from it — the mass-only bracket is the observed shift,
    /// the named tags with formulas are the fixed forms — and setting them
    /// separately is an error. This is the parameter-honesty route: every
    /// session derives the frame from the TopPIC record itself instead of
    /// transcribing it by hand (a dropped S1 shifted every site by one).
    #[serde(default)]
    pub proforma: Option<String>,
    /// MS2 scan numbers (msalign SCANS= field), comma-separated; sibling
    /// scans of one precursor pool their peaks.
    #[serde(default)]
    pub scans: Vec<String>,
    /// Protein sequence, one-letter uppercase. Required unless `proforma`
    /// derives it.
    #[serde(default)]
    pub sequence: String,
    /// Positional fixed modifications (TopPIC ProForma positions).
    #[serde(default)]
    pub fixed_modifications: Vec<PositionalModification>,
    /// Modification library (name, delta mass, residues, location).
    pub library: Vec<LibrarySpec>,
    /// Max total components (library + SHIFT) per combination.
    pub max_components: usize,
    /// The observed unexpected shift (Da) the library + free mass must
    /// explain. Required unless `proforma` derives it (mass-only bracket).
    #[serde(default)]
    pub observed_shift: f64,
    /// Disable the free SHIFT (require the library to close the shift).
    #[serde(default)]
    pub no_free_shift: bool,
    /// Fragment match tolerance (ppm, default 10 like IP's product tolerance).
    #[serde(default = "default_tolerance_ppm")]
    pub tolerance_ppm: f64,
    /// Fragment ion series (default "c,z" — the ETD ladder our data uses).
    #[serde(default = "default_ion_types")]
    pub ion_types: String,
    /// Fragment neutral-mass bounds.
    #[serde(default = "default_mass_min")]
    pub mass_min: f64,
    #[serde(default = "default_mass_max")]
    pub mass_max: f64,
    /// Scorer: "weighted" (default; IP `CompositeScorer` weights over
    /// Count/Intensity/TopFD-quality/MassError, with unavailable Dist=0) or
    /// "count" (diagnostic IP `DeconvScorer`, +1 per matched ion bin).
    #[serde(default = "default_scorer")]
    pub scorer: String,
    /// Skip combinations whose library mass overshoots the observed shift
    /// by more than this window (Da); co-eluting mixtures calibrate between
    /// forms, so a small overshoot stays evaluable.
    #[serde(default = "default_overshoot")]
    pub overshoot_window_da: f64,
    /// Output cap on the ranked candidate list. The search itself is
    /// complete; this only limits the returned view (ties are counted on
    /// the FULL list before truncation, and the note records how many were
    /// cut).
    #[serde(default = "default_max_candidates")]
    pub max_candidates: usize,
}

fn default_tolerance_ppm() -> f64 {
    10.0
}
fn default_ion_types() -> String {
    "c,z".to_string()
}
fn default_mass_min() -> f64 {
    50.0
}
fn default_mass_max() -> f64 {
    100_000.0
}
fn default_scorer() -> String {
    "weighted".to_string()
}
fn default_overshoot() -> f64 {
    0.06
}
fn default_max_candidates() -> usize {
    5
}

#[derive(Clone, Debug, serde::Serialize)]
#[serde(rename_all = "camelCase")]
pub struct GraphSearchCandidate {
    pub combination: String,
    pub components: Vec<String>,
    pub combination_mass: f64,
    pub free_mass: Option<f64>,
    pub mass_error_da: f64,
    pub score: f64,
    pub modifications: Vec<(String, String)>,
    /// Complete, exact-chemistry model ready for `theory` and `fragment`.
    /// Present when the sequence-graph engine received a TopPIC ProForma frame. The
    /// observed mass-only SHIFT is replaced by the candidate's named library
    /// components; a remaining free SHIFT is deliberately omitted. Rendered from
    /// MEASUREMENT pins (component-cliff resolutions) whenever every component
    /// localized — never from the DP's tie-canonical placements, which stay
    /// audit-only in `modifications` and `dpSite`. Round-trip validated at the
    /// parse boundary (see `render_measured_proforma`).
    pub proforma: Option<String>,
    pub num_placements: f64,
    /// Expected number of mass-consistent random sequences scoring at least
    /// `score` (IP spectral E-value; None when the GF layer is skipped).
    pub spectral_evalue: Option<f64>,
    /// IP `CompositeScorer.GetProbability` (logistic on Beta0 + score);
    /// weighted scorer only — count scores have no IP probability analogue.
    pub score_probability: Option<f64>,
    /// Distinct optimal arrangements (tie group size actually enumerated,
    /// capped) and up to five of them, when more than one is optimal.
    pub alternative_placements: Option<usize>,
    pub alternatives: Vec<Vec<(String, String)>>,
    /// One-flow free-shift localization (MSPathFinder semantics: closed
    /// search leaves a shift, the shift gets an approximate localization —
    /// E56 coordinate-ascent protocol): likelihood-ratio evidence that the
    /// free mass is one real modification at a specific site, vs the
    /// no-shift null. A measurement, never a sort key — the composition
    /// `score` above stays the ranking.
    pub shift_evidence_log10: Option<f64>,
    pub shift_best_site: Option<usize>,
    pub shift_site_interval: Option<(usize, usize)>,
    /// Per-component cliff profile from the ascent round (each library
    /// member re-localized with the SHIFT pinned at its own site): a real
    /// model's components each carry a cliff, imposters are flat.
    pub component_cliffs: Vec<ComponentCliff>,
    /// The library closes the mass to calibration scale (free residual
    /// < 0.05 Da): there is no free mass to localize — the discovery
    /// measurement for this candidate is its componentCliffs. Such
    /// candidates form their own view band, never cut by the evidence
    /// cliff (a closed model is the enumeration's answer, not a weak one).
    pub closed_by_library: bool,
}

fn candidate_proforma(
    source: Option<&str>,
    modifications: &[(String, String)],
) -> Result<Option<String>, MassSpecError> {
    let Some(source) = source else {
        return Ok(None);
    };
    let parsed = crate::proforma::parse_proforma(source, &[])?;
    let mut at_site: std::collections::BTreeMap<usize, Vec<String>> =
        std::collections::BTreeMap::new();
    let mut nterm = Vec::new();

    for modification in parsed.modifications.iter().filter(|m| m.formula.is_some()) {
        match modification.position.as_deref() {
            Some("nTerm") => nterm.push(modification.tag.clone()),
            Some(position) => {
                let position = position.strip_prefix("residue:").unwrap_or(position);
                if let Ok(position) = position.parse::<usize>() {
                    at_site
                        .entry(position)
                        .or_default()
                        .push(modification.tag.clone());
                }
            }
            None => {}
        }
    }
    for (name, position) in modifications {
        if name == "SHIFT" || name.starts_with("fixed(") {
            continue;
        }
        let tag = name.split('@').next().unwrap_or(name).to_owned();
        if position == "nTerm" {
            nterm.push(tag);
        } else if let Ok(position) = position.parse::<usize>() {
            at_site.entry(position).or_default().push(tag);
        }
    }

    let mut rendered = String::new();
    for tag in nterm {
        rendered.push('[');
        rendered.push_str(&tag);
        rendered.push_str("]-");
    }
    for (index, residue) in parsed.sequence.chars().enumerate() {
        rendered.push(residue);
        if let Some(tags) = at_site.get(&(index + 1)) {
            for tag in tags {
                rendered.push('[');
                rendered.push_str(tag);
                rendered.push(']');
            }
        }
    }
    Ok(Some(rendered))
}

/// Re-render a candidate ProForma from MEASUREMENT pins (component-cliff
/// resolutions). The rendered string is the consumable model ("ready for
/// `theory` and `fragment`"), so it must carry the measurement — never the
/// DP's tie-canonical placements, which stay audit-only in `modifications`
/// and `dpSite`. The result is round-trip validated through `parse_proforma`:
/// every rendered tag must land back exactly at its pinned 1-based residue,
/// so a placement-convention bug is a hard contract error, not a silent
/// artifact-layer drift. Returns `Ok(None)` when any component site is
/// unresolved (`*`-free numeric positions only) — the caller then keeps the
/// DP-arrangement rendering and has already noted it.
fn render_measured_proforma(
    source: Option<&str>,
    measured: &[(String, String)],
    label: &str,
    notes: &mut Vec<String>,
) -> Result<Option<String>, MassSpecError> {
    let all_resolved = measured
        .iter()
        .all(|(_, position)| position.parse::<usize>().is_ok());
    if !all_resolved {
        notes.push(format!(
            "candidate {label}: not every component localized — its proforma keeps the \
             DP tie-canonical placements (audit only; read componentCliffs for the \
             measurement)"
        ));
        return Ok(None);
    }
    let Some(rendered) = candidate_proforma(source, measured)? else {
        return Ok(None);
    };
    let Ok(parsed) = crate::proforma::parse_proforma(&rendered, &[]) else {
        // Tag names outside the chemistry tables: placement is still
        // index-based and contract-tested; note the skipped validation.
        notes.push(format!(
            "candidate {label}: measured proforma round-trip skipped (tag names outside \
             the chemistry tables)"
        ));
        return Ok(Some(rendered));
    };
    for (name, position) in measured {
        let index: usize = position.parse().expect("all_resolved checked above");
        let expected_tag = name.split('@').next().unwrap_or(name);
        let residue = index
            .checked_sub(1)
            .and_then(|offset| parsed.sequence.chars().nth(offset));
        let Some(residue) = residue else {
            // The renderer silently drops out-of-range pins — exactly the
            // silent-drift class this contract exists to kill.
            return Err(MassSpecError::InvalidQuery(format!(
                "measured proforma round-trip failed: [{name}] pin {index} is outside \
                 the sequence — rendering/measurement convention mismatch"
            )));
        };
        let landed = parsed.modifications.iter().any(|m| {
            m.tag == expected_tag
                && m.position.as_deref() == Some(format!("residue:{index}").as_str())
                && m.residue == Some(residue)
        });
        if !landed {
            return Err(MassSpecError::InvalidQuery(format!(
                "measured proforma round-trip failed: [{name}] did not land back at \
                 residue {index} ({residue}) — rendering/measurement convention mismatch"
            )));
        }
    }
    Ok(Some(rendered))
}

/// One library component's leave-one-out cliff in a free-shift candidate.
#[derive(Clone, Debug, serde::Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ComponentCliff {
    pub name: String,
    /// DP placement (tie-canonical first-best) — audit trail only, the
    /// measurement never pins to it.
    pub dp_site: String,
    /// The component's own cliff: unique site when the likelihood top is a
    /// single residue, else the flat-top interval (an observation gap —
    /// every observable fragment either contains the whole interval or
    /// misses it). (0,0) = not localized by this measurement.
    pub site_interval: (usize, usize),
    pub evidence_log10: Option<f64>,
}

#[derive(Debug, serde::Serialize)]
#[serde(rename_all = "camelCase")]
pub struct GraphSearchResult {
    pub sequence: String,
    pub observed_shift: f64,
    pub shift_frame: String,
    pub scans: Vec<String>,
    pub ions: Vec<String>,
    pub candidates: Vec<GraphSearchCandidate>,
    pub enumerated_combinations: usize,
    pub evaluated_combinations: usize,
    pub max_candidates: usize,
    pub notes: Vec<String>,
}

/// Derive the search frame (sequence, fixed forms, observed shift) from a
/// TopPIC ProFormA. The mass-only bracket (no elemental formula) is THE
/// observed unexpected shift; every named tag with a formula becomes a
/// pinned fixed modification. Refusing explicit duplicates of derived
/// fields keeps exactly one authority for the frame.
/// The frame a TopPIC record hands a measurement: the sequence plus the
/// pinned fixed forms (every named tag with a formula). The mass-only
/// bracket is NOT consumed here — the caller decides what it means
/// (the sequence-graph engine consumes it as the observed shift; callers
/// comparing complete hypotheses must construct those models explicitly).
pub fn proforma_frame(text: &str) -> Result<(String, Vec<MassModification>), MassSpecError> {
    let parsed = crate::proforma::parse_proforma(text, &[])?;
    let mut pins = Vec::new();
    for m in &parsed.modifications {
        if m.formula.is_none() {
            continue;
        }
        let position = match m.position.as_deref() {
            Some("nTerm") | Some("cTerm") => m.position.clone().unwrap(),
            Some(pos) if pos.starts_with("residue:") => {
                pos.trim_start_matches("residue:").to_string()
            }
            other => {
                return Err(MassSpecError::InvalidQuery(format!(
                    "proforma fixed modification [{}] has no pinnable site ({other:?}) — \
                     resolve the site first",
                    m.tag
                )))
            }
        };
        pins.push(MassModification {
            position,
            delta_mass: m.delta_mass,
            formula: m.formula.clone(),
        });
    }
    Ok((parsed.sequence, pins))
}

/// Derive the search frame (sequence, fixed forms, observed shift) from a
/// TopPIC ProFormA. The mass-only bracket (no elemental formula) is THE
/// observed unexpected shift. Refusing explicit duplicates of derived
/// fields keeps exactly one authority for the frame.
fn derive_frame_from_proforma(
    query: &mut GraphSearchQuery,
    proforma: &str,
) -> Result<(), MassSpecError> {
    let parsed = crate::proforma::parse_proforma(proforma, &[])?;
    let mass_only: Vec<&crate::proforma::ParsedModification> = parsed
        .modifications
        .iter()
        .filter(|m| m.formula.is_none())
        .collect();
    if mass_only.len() != 1 {
        return Err(MassSpecError::InvalidQuery(format!(
            "proforma must carry exactly one mass-only bracket (the observed shift); got {} — \
             decompose per bracket or run per PrSM",
            mass_only.len()
        )));
    }
    if !query.sequence.trim().is_empty()
        || query.observed_shift != 0.0
        || !query.fixed_modifications.is_empty()
    {
        return Err(MassSpecError::InvalidQuery(
            "--sequence/--shift/--mods conflict with --proforma: the proforma is the \
             single authority for the frame"
                .into(),
        ));
    }
    let (sequence, pins) = proforma_frame(proforma)?;
    query.sequence = sequence;
    query.observed_shift = mass_only[0].delta_mass;
    query.fixed_modifications = pins
        .into_iter()
        .map(|m| PositionalModification {
            position: m.position,
            delta_mass: m.delta_mass,
        })
        .collect();
    Ok(())
}

pub fn graph_search_msalign(
    bytes: &[u8],
    query: &GraphSearchQuery,
) -> Result<GraphSearchResult, MassSpecError> {
    if query.library.is_empty() && query.no_free_shift {
        return Err(MassSpecError::InvalidQuery(
            "an empty library cannot close the shift alone; allow the free SHIFT or provide --library"
                .into(),
        ));
    }
    let mut query = query.clone();
    if let Some(proforma) = query.proforma.clone() {
        derive_frame_from_proforma(&mut query, &proforma)?;
    } else if query.sequence.trim().is_empty() || query.observed_shift == 0.0 {
        return Err(MassSpecError::InvalidQuery(
            "resolve-shift needs --sequence and --shift, or --proforma to derive them".into(),
        ));
    }
    let mut spectra = Vec::new();
    for scan in &query.scans {
        let spectrum = parse_msalign_spectrum(bytes, scan)?;
        if spectrum.level != 2 {
            return Err(MassSpecError::InvalidQuery(format!(
                "resolve-shift needs an MS2 spectrum; scan {scan} is level {}",
                spectrum.level
            )));
        }
        spectra.push(spectrum);
    }
    if has_mixed_activation(&spectra) {
        return Err(MassSpecError::InvalidQuery(
            "resolve-shift cannot pool ETD/ECD with HCD/CID because one search uses one fragment-ion series; run each activation separately"
                .into(),
        ));
    }
    query.ion_types = resolve_ion_types(
        &query.ion_types,
        spectra.first().and_then(|s| s.activation.as_deref()),
    )?;
    let ions = BaseIon::parse_set(&query.ion_types)?;
    let peaks: Vec<DeconvolutedPeak> = spectra
        .iter()
        .flat_map(|s| s.peaks.iter().copied())
        .collect();
    // The GF null model is defined over the counting observable, so spectral
    // E-values run only with the explicit counting scorer.
    let weighted = match query.scorer.trim().to_ascii_lowercase().as_str() {
        "weighted" => true,
        "count" => false,
        name => {
            return Err(MassSpecError::InvalidQuery(format!(
                "unknown resolve-shift scorer {name:?}; expected weighted or count"
            )))
        }
    };
    let mut scorer: Box<dyn GraphScorer> = if weighted {
        Box::new(WeightedBinScorer::new(
            &peaks,
            ions.clone(),
            query.tolerance_ppm,
            query.mass_min,
            query.mass_max,
        ))
    } else {
        Box::new(BinCountScorer::new(
            &peaks,
            ions.clone(),
            query.tolerance_ppm,
            query.mass_min,
            query.mass_max,
        ))
    };
    let mut result = graph_search_run(&query, scorer.as_mut(), &ions, &peaks, !weighted, weighted)?;
    // One-flow integration: the search has spoken (enumeration + gating);
    // each free-shift candidate now gets its SHIFT localized — closed
    // search, leftover shift, approximate localization, one command.
    let scan_spectra: Vec<(String, MsAlignSpectrum)> =
        query.scans.iter().cloned().zip(spectra).collect();
    localize_free_shifts(&mut result, &scan_spectra, &query)?;
    rank_by_discovery(&mut result);
    // view cap AFTER the discovery ranking — never between measurement
    // and view (the E53 lesson: the truth sat at composition rank 14
    // under a cap of 10 and died unmeasured). Library-closed candidates
    // lead the list, so a plain tail truncation can never cut them
    // unless the closed band itself exceeds the cap.
    if result.candidates.len() > query.max_candidates {
        let truncated = result.candidates.len() - query.max_candidates;
        result.candidates.truncate(query.max_candidates);
        result.notes.push(format!(
            "{truncated} candidates truncated by maxCandidates={} — the search and \
             the discovery ranking were complete; request a higher cap to see the tail",
            query.max_candidates
        ));
    }
    Ok(result)
}

fn has_mixed_activation(spectra: &[MsAlignSpectrum]) -> bool {
    let activations: Vec<String> = spectra
        .iter()
        .filter_map(|s| s.activation.as_deref())
        .map(|a| a.to_ascii_uppercase())
        .collect();
    let has_etd = activations
        .iter()
        .any(|a| a.contains("ETD") || a.contains("ECD"));
    let has_hcd = activations
        .iter()
        .any(|a| a.contains("HCD") || a.contains("CID"));
    has_etd && has_hcd
}

/// Discovery-ordered view (E59/E61): LIBRARY-CLOSED candidates lead
/// (composition order) — a model that explains the whole shift with known
/// chemistry is the preferred explanation; the open-residual route is the
/// discovery channel for when the library CANNOT close (MSPathFinder's
/// closed-search-first philosophy). Then the free-shift candidates carrying
/// a localization measurement, strongest evidence at the top; candidates
/// without the measurement keep composition order after them. When the
/// largest adjacent evidence gap reaches DISCOVERY_CLIFF_LOG10, everything
/// below the cliff is truncated from the view — the cut, its size and the
/// gap are recorded in the notes; the composition evaluation counts are
/// untouched. Order differences smaller than the cliff INSIDE the kept
/// head are mass degeneracy, not a verdict.
fn rank_by_discovery(result: &mut GraphSearchResult) {
    if result
        .candidates
        .iter()
        .all(|c| c.shift_evidence_log10.is_none())
    {
        return; // nothing measured (closed search / all flat): composition order stands
    }
    let mut closed: Vec<GraphSearchCandidate> = Vec::new();
    let mut with_evidence: Vec<GraphSearchCandidate> = Vec::new();
    let mut without: Vec<GraphSearchCandidate> = Vec::new();
    for candidate in result.candidates.drain(..) {
        if candidate.free_mass.is_none() || candidate.closed_by_library {
            // SHIFT-less combinations and calibration-scale residuals are
            // both "the library closes the mass"
            closed.push(candidate);
        } else if candidate.shift_evidence_log10.is_some() {
            with_evidence.push(candidate);
        } else {
            without.push(candidate);
        }
    }
    // closed band keeps composition order (the sort in graph_search_run);
    // the discovery axis orders the open candidates
    with_evidence.sort_by(|a, b| {
        b.shift_evidence_log10
            .partial_cmp(&a.shift_evidence_log10)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    // cliff truncation at the LARGEST adjacent gap when it reaches the
    // significance convention (E2: a gap < 2 log10 is a tie, not a verdict)
    let mut cut_after = None;
    let mut largest_gap = DISCOVERY_CLIFF_LOG10;
    for i in 0..with_evidence.len().saturating_sub(1) {
        // Candidates without an evidence score sort last; treat them as
        // -inf instead of panicking if one slips into the gap window.
        let gap = (with_evidence[i]
            .shift_evidence_log10
            .unwrap_or(f64::NEG_INFINITY)
            - with_evidence[i + 1]
                .shift_evidence_log10
                .unwrap_or(f64::NEG_INFINITY))
        .abs();
        if gap >= largest_gap {
            largest_gap = gap;
            cut_after = Some(i);
        }
    }
    if let Some(i) = cut_after {
        let dropped = with_evidence.len() - (i + 1) + without.len();
        let gap = largest_gap;
        with_evidence.truncate(i + 1);
        without.clear();
        result.notes.push(format!(
            "{dropped} candidates truncated below a {gap:.1} log10 evidence cliff — \
             the kept head is separated by an order of magnitude or more; inside the head, \
             order differences smaller than {DISCOVERY_CLIFF_LOG10:.0} log10 are mass \
             degeneracy (the duel's job, not the ranking's); library-closed candidates are \
             kept below the head in their own band — they are the enumeration's answers, not \
             weak evidence"
        ));
    }
    result.notes.push(
        "candidate order: library-closed models first (composition order — the preferred \
         explanation when the library closes the mass), then the discovery measurement \
         (free-shift localization evidence), then the rest; the composition score stays \
         the gating axis — read both"
            .to_string(),
    );
    result.candidates = closed;
    result.candidates.extend(with_evidence);
    result.candidates.extend(without);
}

/// Adjacent-gap threshold for the discovery-view truncation: the log's
/// significance convention (E2) — a runner-up gap below this is a tie.
const DISCOVERY_CLIFF_LOG10: f64 = 2.0;

/// A localize result counts as resolved when its likelihood top is
/// non-empty (`site_interval != (0,0)`): a unique site OR an honest
/// observation-gap interval. Pinning anywhere inside a flat top is
/// mass-equivalent for every observable fragment — pin at its start.
fn resolved_pin(r: &PrsmLocalizeResult) -> Option<usize> {
    (r.site_interval != (0, 0)).then_some(r.site_interval.0)
}

/// Derive the component sweep's legal residue positions from the same
/// `LibrarySpec` used to build the sequence graph. Terminal-only locations
/// intentionally return an empty set: the localizer measures residue sites,
/// so an unresolved terminal component keeps its valid DP placement instead
/// of being reinterpreted as a residue modification.
fn component_allowed_positions(query: &GraphSearchQuery) -> HashMap<String, Vec<usize>> {
    query
        .library
        .iter()
        .map(|spec| {
            let location = spec
                .location
                .as_deref()
                .and_then(SequenceLocation::parse)
                .unwrap_or(SequenceLocation::Everywhere);
            let positions = if location != SequenceLocation::Everywhere {
                Vec::new()
            } else {
                query
                    .sequence
                    .chars()
                    .enumerate()
                    .filter_map(|(offset, residue)| {
                        let position = offset + 1;
                        let residue_allowed = spec.residues.is_empty()
                            || spec.residues.contains(&'*')
                            || spec.residues.contains(&residue);
                        let position_allowed = spec.positions.as_ref().is_none_or(|allowed| {
                            allowed
                                .iter()
                                .any(|value| value.parse::<usize>().ok() == Some(position))
                        });
                        (residue_allowed && position_allowed).then_some(position)
                    })
                    .collect()
            };
            (spec.name.clone(), positions)
        })
        .collect()
}

/// Two-pass per-component localization sweep. Pass 0 measures every
/// component alone (others unknown ⇒ excluded from the ladder — fragments
/// spanning them become unobserved rather than mis-predicted); pass 1
/// re-measures each with every OTHER resolved component pinned (full
/// leave-one-out). A component that never resolves stays EXCLUDED —
/// pinning it anywhere would corrupt the voters it wrongly spans;
/// resolutions are sticky (a weaker pass-1 never drops a pass-0 site).
/// `free_pin` = (site, mass) pins the free mass during the sweep; None =
/// closed model (library closes the mass, nothing to pin).
fn sweep_components(
    run: &mut dyn FnMut(
        f64,
        &[MassModification],
        Option<&[usize]>,
        &mut Vec<String>,
    ) -> Option<PrsmLocalizeResult>,
    fixed_pins: &[MassModification],
    dp_members: &[(String, String)],
    delta_of: &HashMap<&str, f64>,
    allowed_positions: &HashMap<String, Vec<usize>>,
    free_pin: Option<(usize, f64)>,
    notes: &mut Vec<String>,
) -> (Vec<ComponentCliff>, Vec<Option<MassModification>>) {
    let mut component_cliffs = Vec::with_capacity(dp_members.len());
    let mut comp_pins: Vec<Option<MassModification>> = vec![None; dp_members.len()];
    for pass in 0..2 {
        let mut pass_pins = Vec::with_capacity(dp_members.len());
        let mut pass_cliffs = Vec::with_capacity(dp_members.len());
        for (i, (name, dp_pos)) in dp_members.iter().enumerate() {
            let delta = delta_of[name.as_str()];
            let mut pins = fixed_pins.to_vec();
            if pass == 1 {
                for (j, other) in comp_pins.iter().enumerate() {
                    if j != i {
                        if let Some(p) = other {
                            pins.push(p.clone());
                        }
                    }
                }
            }
            if let Some((site, mass)) = free_pin {
                pins.push(MassModification {
                    position: site.to_string(),
                    delta_mass: mass,
                    formula: None,
                });
            }
            let allowed = allowed_positions.get(name).map(Vec::as_slice);
            let rc = if allowed.is_some_and(|positions| positions.is_empty()) {
                None
            } else {
                run(delta, &pins, allowed, notes)
            };
            let pin = rc
                .as_ref()
                .and_then(resolved_pin)
                .map(|s| MassModification {
                    position: s.to_string(),
                    delta_mass: delta,
                    formula: None,
                });
            pass_cliffs.push(ComponentCliff {
                name: name.clone(),
                dp_site: dp_pos.clone(),
                site_interval: rc.as_ref().map(|r| r.site_interval).unwrap_or((0, 0)),
                evidence_log10: rc.as_ref().and_then(|r| r.evidence_log10),
            });
            pass_pins.push(pin.or_else(|| comp_pins[i].clone()));
        }
        comp_pins = pass_pins;
        component_cliffs = pass_cliffs;
    }
    (component_cliffs, comp_pins)
}

/// Re-localize every free-shift candidate's SHIFT. Protocol (three steps,
/// no dependence on DP tie placement):
/// 1. the free mass discovers its own site on the bare ladder (fixed mods
///    only) — fragments that span the site and no library region vote with
///    no pinning assumptions;
/// 2. the free mass pinned at that site, each library component is
///    localized in turn (already-measured components pinned at their own
///    sites) — the per-component cliff profile;
/// 3. the free mass is re-measured with every component at its own cliff
///    site — the headline evidence.
///
/// Fallback when step 1 is flat: the DP's composition-tied arrangements
/// (canonical + alternatives) are measured as multi-starts and the
/// localization likelihood picks the reported one (E56).
fn localize_free_shifts(
    result: &mut GraphSearchResult,
    spectra: &[(String, MsAlignSpectrum)],
    query: &GraphSearchQuery,
) -> Result<(), MassSpecError> {
    if spectra.is_empty() {
        return Ok(());
    }
    let delta_of: HashMap<&str, f64> = query
        .library
        .iter()
        .map(|spec| (spec.name.as_str(), spec.delta_mass))
        .collect();
    let fixed_pins: Vec<MassModification> = query
        .fixed_modifications
        .iter()
        .map(|m| MassModification {
            position: m.position.clone(),
            delta_mass: m.delta_mass,
            formula: None,
        })
        .collect();
    let allowed_positions = component_allowed_positions(query);
    let mut localize_error_noted = false;
    let mut run = |delta: f64,
                   pinned: &[MassModification],
                   allowed: Option<&[usize]>,
                   notes: &mut Vec<String>|
     -> Option<PrsmLocalizeResult> {
        let lq = PrsmLocalizeQuery {
            scans: query.scans.clone(),
            sequence: query.sequence.clone(),
            delta_mass: delta,
            modifications: pinned.to_vec(),
            tolerance_ppm: query.tolerance_ppm,
            mass_min: query.mass_min,
            mass_max: query.mass_max,
            range: None,
        };
        match query_prsm_localize_spectra_restricted(spectra, &lq, allowed) {
            Ok(r) => {
                if std::env::var_os("TDEASE_DEBUG_LOCALIZE").is_some() {
                    eprintln!(
                        "[localize] delta={} pins={:?} -> ev={} best={} ivl={}-{}",
                        lq.delta_mass,
                        lq.modifications,
                        r.evidence_log10
                            .map(|v| format!("{v:.1}"))
                            .unwrap_or_else(|| "-".into()),
                        r.best_site
                            .map(|v| v.to_string())
                            .unwrap_or_else(|| "-".into()),
                        r.site_interval.0,
                        r.site_interval.1
                    );
                }
                Some(r)
            }
            Err(e) => {
                if !localize_error_noted {
                    localize_error_noted = true;
                    notes.push(format!(
                        "free-shift localization skipped for at least one candidate: {e}"
                    ));
                }
                None
            }
        }
    };
    for candidate in result.candidates.iter_mut() {
        let Some(free) = candidate.free_mass else {
            continue;
        };
        if !free.is_finite() || free.abs() < 1e-4 {
            continue;
        }
        // Library members of this candidate (fixed("...") instances pin from
        // the query spec; the SHIFT is the mass being swept).
        let dp_members: Vec<(String, String)> = candidate
            .modifications
            .iter()
            .filter(|(name, _)| name.as_str() != SHIFT_NAME && delta_of.contains_key(name.as_str()))
            .cloned()
            .collect();
        if free.abs() < 0.05 {
            // A residual this small is calibration scale (TopPIC's own
            // precursor window is ±0.4 Da), not a modification hypothesis:
            // localizing a 0.02 Da shift measures nothing. The model is
            // CLOSED by the library — the discovery measurement is the
            // per-component cliffs (no free mass to pin).
            let members = dp_members.clone();
            let (component_cliffs, comp_pins) = sweep_components(
                &mut run,
                &fixed_pins,
                &members,
                &delta_of,
                &allowed_positions,
                None,
                &mut result.notes,
            );
            candidate.component_cliffs = component_cliffs;
            candidate.closed_by_library = true;
            // The consumable model carries the measurement: re-render from
            // the cliff pins (falling back to the DP arrangement, noted,
            // when some component never localized).
            // Unresolved components keep a non-numeric marker so the
            // renderer falls back to the DP arrangement instead of silently
            // producing a mass-incomplete model.
            let measured: Vec<(String, String)> = candidate
                .component_cliffs
                .iter()
                .zip(comp_pins.iter())
                .map(|(cliff, pin)| {
                    let position = pin
                        .as_ref()
                        .map(|p| p.position.clone())
                        .unwrap_or_else(|| "unresolved".to_owned());
                    (cliff.name.clone(), position)
                })
                .collect();
            if let Some(rendered) = render_measured_proforma(
                query.proforma.as_deref(),
                &measured,
                &candidate.combination,
                &mut result.notes,
            )? {
                candidate.proforma = Some(rendered);
            }
            candidate.combination = candidate.combination.trim_end_matches("+SHIFT").to_owned();
            candidate.components.retain(|name| name != SHIFT_NAME);
            candidate
                .modifications
                .retain(|(name, _)| name != SHIFT_NAME);
            candidate.combination_mass -= free;
            candidate.mass_error_da = -free;
            candidate.free_mass = None;
            result.notes.push(format!(
                "candidate {}: signed mass error {:.4} Da is calibration scale — the library \
                 closes the mass with the complete candidate ProForma; the discovery measurement \
                 is the component localization above",
                candidate.combination, free
            ));
            continue;
        }
        // Step 1: the free mass discovers its own site on the bare ladder.
        let bare = run(free, &fixed_pins, None, &mut result.notes);
        if let Some(site0) = bare.as_ref().and_then(resolved_pin) {
            // Step 2: each component localized with the free mass pinned at
            // its discovered site (two-pass leave-one-out — see
            // sweep_components). A component that never resolves stays
            // EXCLUDED — pinning it anywhere would corrupt the voters it
            // wrongly spans.
            let (component_cliffs, comp_pins) = sweep_components(
                &mut run,
                &fixed_pins,
                &dp_members,
                &delta_of,
                &allowed_positions,
                Some((site0, free)),
                &mut result.notes,
            );
            // Step 3: headline measurement — free re-localized with every
            // RESOLVED component at its own cliff site (the free pin is
            // dropped: it is the mass being swept; unresolved components
            // stay excluded).
            let mut final_pins = fixed_pins.clone();
            final_pins.extend(comp_pins.iter().flatten().cloned());
            let final_r = if dp_members.is_empty() {
                bare
            } else {
                run(free, &final_pins, None, &mut result.notes)
            };
            let site = final_r.as_ref().and_then(|r| r.best_site);
            let interval = final_r.as_ref().map(|r| r.site_interval).unwrap_or((0, 0));
            candidate.shift_evidence_log10 = final_r.as_ref().and_then(|r| r.evidence_log10);
            candidate.shift_best_site = site;
            candidate.shift_site_interval = (interval != (0, 0)).then_some(interval);
            candidate.component_cliffs = component_cliffs;
            // The consumable model carries the measurement: components at
            // their cliff pins (a remaining free SHIFT stays deliberately
            // omitted — its localized site lives in the fields above).
            // Unresolved components keep a non-numeric marker so the
            // renderer falls back to the DP arrangement instead of silently
            // producing a mass-incomplete model.
            let measured: Vec<(String, String)> = candidate
                .component_cliffs
                .iter()
                .zip(comp_pins.iter())
                .map(|(cliff, pin)| {
                    let position = pin
                        .as_ref()
                        .map(|p| p.position.clone())
                        .unwrap_or_else(|| "unresolved".to_owned());
                    (cliff.name.clone(), position)
                })
                .collect();
            if let Some(rendered) = render_measured_proforma(
                query.proforma.as_deref(),
                &measured,
                &candidate.combination,
                &mut result.notes,
            )? {
                candidate.proforma = Some(rendered);
            }
            continue;
        }
        // Fallback (bare curve flat): multi-start over the DP's distinct
        // composition-tied library arrangements; the localization
        // likelihood — not canonical order — picks the reported one.
        let mut starts: Vec<Vec<(String, String)>> = Vec::new();
        let mut push_start = |arr: &[(String, String)]| {
            let mut lib: Vec<(String, String)> = arr
                .iter()
                .filter(|(name, _)| {
                    name.as_str() != SHIFT_NAME && delta_of.contains_key(name.as_str())
                })
                .cloned()
                .collect();
            lib.sort();
            if !lib.is_empty() && !starts.contains(&lib) {
                starts.push(lib);
            }
        };
        push_start(&candidate.modifications);
        for alt in &candidate.alternatives {
            push_start(alt);
        }
        if starts.len() > 6 {
            starts.truncate(6);
        }
        let mut best: Option<(f64, usize, PrsmLocalizeResult)> = None;
        for (si, lib) in starts.iter().enumerate() {
            let pins: Vec<MassModification> = fixed_pins
                .iter()
                .cloned()
                .chain(lib.iter().map(|(name, pos)| MassModification {
                    position: pos.clone(),
                    delta_mass: delta_of[name.as_str()],
                    formula: None,
                }))
                .collect();
            if let Some(r) = run(free, &pins, None, &mut result.notes) {
                let ev = r.evidence_log10.unwrap_or(f64::NEG_INFINITY);
                if best.as_ref().is_none_or(|(b, _, _)| ev > *b) {
                    best = Some((ev, si, r));
                }
            }
        }
        let Some((_, si, r0)) = best else {
            candidate.shift_evidence_log10 = bare.as_ref().and_then(|r| r.evidence_log10);
            continue;
        };
        let dp_pins = starts[si].clone();
        let Some(site0) = resolved_pin(&r0) else {
            candidate.shift_evidence_log10 = r0.evidence_log10;
            continue;
        };
        let mut component_cliffs = Vec::with_capacity(dp_pins.len());
        let mut own_pins = fixed_pins.clone();
        let mut pin_sites: Vec<Option<usize>> = Vec::with_capacity(dp_pins.len());
        for i in 0..dp_pins.len() {
            let (name, dp_pos) = &dp_pins[i];
            let delta = delta_of[name.as_str()];
            let mut pins = fixed_pins.clone();
            for (j, (other, other_pos)) in dp_pins.iter().enumerate() {
                if j != i {
                    pins.push(MassModification {
                        position: other_pos.clone(),
                        delta_mass: delta_of[other.as_str()],
                        formula: None,
                    });
                }
            }
            pins.push(MassModification {
                position: site0.to_string(),
                delta_mass: free,
                formula: None,
            });
            let allowed = allowed_positions.get(name).map(Vec::as_slice);
            let rc = if allowed.is_some_and(|positions| positions.is_empty()) {
                None
            } else {
                run(delta, &pins, allowed, &mut result.notes)
            };
            let pin_site = rc.as_ref().and_then(resolved_pin);
            component_cliffs.push(ComponentCliff {
                name: name.clone(),
                dp_site: dp_pos.clone(),
                site_interval: rc.as_ref().map(|r| r.site_interval).unwrap_or((0, 0)),
                evidence_log10: rc.as_ref().and_then(|r| r.evidence_log10),
            });
            // unresolved ⇒ excluded (pinning it anywhere would corrupt the
            // voters it wrongly spans)
            if let Some(s) = pin_site {
                own_pins.push(MassModification {
                    position: s.to_string(),
                    delta_mass: delta,
                    formula: None,
                });
            }
            pin_sites.push(pin_site);
        }
        // All components must have localized; otherwise keep the DP
        // arrangement rendering (the renderer notes the fallback).
        let measured: Vec<(String, String)> = dp_pins
            .iter()
            .zip(pin_sites.iter())
            .map(|((name, _), site)| {
                (
                    name.clone(),
                    site.map(|s| s.to_string())
                        .unwrap_or_else(|| "unresolved".to_owned()),
                )
            })
            .collect();
        let final_r = if dp_pins.is_empty() {
            Some(r0)
        } else {
            run(free, &own_pins, None, &mut result.notes)
        };
        let site = final_r.as_ref().and_then(|r| r.best_site);
        let interval = final_r.as_ref().map(|r| r.site_interval).unwrap_or((0, 0));
        candidate.shift_evidence_log10 = final_r.as_ref().and_then(|r| r.evidence_log10);
        candidate.shift_best_site = site;
        candidate.shift_site_interval = (interval != (0, 0)).then_some(interval);
        candidate.component_cliffs = component_cliffs;
        if let Some(rendered) = render_measured_proforma(
            query.proforma.as_deref(),
            &measured,
            &candidate.combination,
            &mut result.notes,
        )? {
            candidate.proforma = Some(rendered);
        }
    }
    Ok(())
}

/// Resolve auto ion types from the spectra's activation method (ETD/ECD
/// -> c/z, HCD/CID -> b/y). Unknown activation must be explicit.
/// Explicit ion lists pass through unchanged.
fn resolve_ion_types(requested: &str, activation: Option<&str>) -> Result<String, MassSpecError> {
    if !requested.trim().eq_ignore_ascii_case("auto") {
        return Ok(requested.to_owned());
    }
    match activation.map(str::to_ascii_uppercase) {
        Some(tag) if tag.contains("ETD") || tag.contains("ECD") => Ok("c,z".to_owned()),
        Some(tag) if tag.contains("HCD") || tag.contains("CID") => Ok("b,y".to_owned()),
        value => Err(MassSpecError::InvalidQuery(format!(
            "resolve-shift could not infer fragment ions from ACTIVATION={value:?}; pass --ion-types explicitly"
        ))),
    }
}

/// The shared engine core: graph construction, top-node loop, ranking, and
/// (count-only) spectral E-values. `peaks` feeds the GF null model and is
/// meaningful only on the deconvoluted msalign path.
fn graph_search_run(
    query: &GraphSearchQuery,
    scorer: &mut dyn GraphScorer,
    ions: &[BaseIon],
    peaks: &[DeconvolutedPeak],
    run_evalue: bool,
    weighted: bool,
) -> Result<GraphSearchResult, MassSpecError> {
    let mut modifications: Vec<SearchModification> = query
        .library
        .iter()
        .map(|spec| {
            let residues = if spec.residues.is_empty() {
                vec!['*']
            } else {
                spec.residues.clone()
            };
            let location = spec
                .location
                .as_deref()
                .and_then(SequenceLocation::parse)
                .unwrap_or(SequenceLocation::Everywhere);
            SearchModification {
                name: spec.name.clone(),
                composition: GraphComposition::from_delta(spec.delta_mass),
                residues,
                location,
                positions: spec.positions.clone(),
            }
        })
        .collect();
    if !query.no_free_shift {
        modifications.push(SearchModification {
            name: SHIFT_NAME.to_string(),
            composition: GraphComposition::from_delta(query.observed_shift),
            residues: vec!['*'],
            location: SequenceLocation::Everywhere,
            positions: None,
        });
    }
    // The free residual is a discovery coordinate, not a PTM component. Give
    // the graph one extra slot when SHIFT is enabled, then enforce the public
    // max_components limit on library members below.
    let graph_component_limit = query.max_components.max(1) + usize::from(!query.no_free_shift);
    let params = ModificationParams::new(
        modifications,
        graph_component_limit,
        Some(query.observed_shift),
    )?;
    let mut graph = SequenceGraph::new(&params, &query.sequence, &query.fixed_modifications)?;

    let mut candidates: Vec<GraphSearchCandidate> = Vec::new();
    let mut evaluated = 0usize;
    let mut notes = Vec::new();
    for (node, comb) in graph.top_nodes() {
        let members = params.combination(comb).to_vec();
        let library_member_count = members
            .iter()
            .filter(|&&m| Some(m) != params.shift_index())
            .count();
        if library_member_count > query.max_components.max(1) {
            continue;
        }
        let library_mass: f64 = members
            .iter()
            .filter(|&&m| Some(m) != params.shift_index())
            .map(|&m| params.modification(m).composition.mass())
            .sum();
        if members.iter().any(|&m| Some(m) == params.shift_index()) {
            // free-shift combination: SHIFT + pinned library members must
            // not overshoot the observed shift beyond the mixture window
            // (co-eluting mixtures calibrate between their forms).
            if library_mass > query.observed_shift + query.overshoot_window_da {
                continue;
            }
        } else {
            // closed mode gates at the fragment tolerance applied to the
            // SHIFT mass (IP gates sinks by precursor tolerance; the 60 mDa
            // mixture window belongs to the free-shift path only). A closed
            // model with a large mass error would otherwise win count
            // scoring on ions that never span the misplaced components —
            // its discriminating ions are exactly the ones that fail.
            let closed_tol = query.observed_shift.abs() * query.tolerance_ppm * 1e-6;
            if (library_mass - query.observed_shift).abs() > closed_tol {
                continue;
            }
        }
        graph.set_sink(node as usize)?;
        let total_mass = params.combination_mass(comb);
        let free = query.observed_shift - library_mass;
        let (score, mods) = graph.get_fragment_score_and_modifications(scorer);
        let (_s2, alternatives, _capped) = graph.optimal_placements(scorer, 6);
        evaluated += 1;
        let placements = graph.num_placements()[node as usize];
        let (alternative_placements, shown_alternatives) = if alternatives.len() > 1 {
            (
                Some(alternatives.len()),
                alternatives.into_iter().take(5).collect(),
            )
        } else {
            (None, Vec::new())
        };
        let proforma = candidate_proforma(query.proforma.as_deref(), &mods)?;
        candidates.push(GraphSearchCandidate {
            combination: params.combination_label(comb),
            components: members
                .iter()
                .map(|&m| params.modification(m).name.clone())
                .collect(),
            combination_mass: total_mass,
            free_mass: if members.iter().any(|&m| Some(m) == params.shift_index()) {
                Some(free)
            } else {
                None
            },
            mass_error_da: total_mass - query.observed_shift,
            score,
            modifications: mods,
            proforma,
            num_placements: placements,
            spectral_evalue: None,
            score_probability: weighted.then(|| get_probability(score)),
            alternative_placements,
            alternatives: shown_alternatives,
            shift_evidence_log10: None,
            shift_best_site: None,
            shift_site_interval: None,
            component_cliffs: Vec::new(),
            closed_by_library: false,
        });
    }
    if candidates.is_empty() {
        notes.push(format!(
            "no combination within overshoot window {:.3} Da of the observed shift; \
             widen the library, max_components, or overshoot_window_da",
            query.overshoot_window_da
        ));
    }
    candidates.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    // Spectral E-values (IP RunGeneratingFunction): one mass-bin scoring
    // graph per DISTINCT candidate protein mass, cached like IP's per-scan
    // distribution cache.
    if weighted {
        notes.push(
            "spectral E-value is skipped for the weighted scorer: the current generating-function null is defined over count scores only".to_string(),
        );
        let cutoff = ip_default_weights().cutoff;
        notes.push(format!(
            "IP ScoreParam.Cutoff = {cutoff:.1}: IcTopDownLauncher only reports matches with Score > Cutoff; candidates at or below it are not IP-reportable"
        ));
    }
    if run_evalue && !candidates.is_empty() {
        // Node scores follow the counting scorer (+1 per matched
        // prefix/suffix ion at a bin) — the GF null space is defined over
        // the same observable the ranking used.
        let residue_masses: Vec<f64> = query
            .sequence
            .chars()
            .filter_map(|c| residue_composition(c).map(|comp| comp.mass()))
            .collect();
        let mut base_mass = 0.0f64;
        for c in query.sequence.chars() {
            if let Some(comp) = residue_composition(c) {
                base_mass += comp.mass();
            }
        }
        base_mass += crate::constants::H2O;
        let fixed_total: f64 = query.fixed_modifications.iter().map(|m| m.delta_mass).sum();
        // the binning must cover the largest PROTEIN mass (base + shift +
        // fixed), not the shift alone — nodes span every prefix mass
        let mass_max = base_mass
            + candidates
                .iter()
                .map(|c| c.combination_mass)
                .fold(0.0f64, f64::max)
            + fixed_total;
        let (prefix_offset, suffix_offset) = {
            let p = ions
                .iter()
                .find(|i| i.offset().0)
                .map(|i| i.offset().1)
                .unwrap_or(0.0);
            let q = ions
                .iter()
                .find(|i| !i.offset().0)
                .map(|i| i.offset().1)
                .unwrap_or(crate::constants::H2O);
            (p, q)
        };
        if mass_max <= 50_000.0 {
            let binning = MassBinning::build(&residue_masses, mass_max + 1.0, 27);
            let mut graph_cache: HashMap<i64, MassBinScoringGraph> = HashMap::new();
            for candidate in candidates.iter_mut() {
                let protein_mass = base_mass + candidate.combination_mass + fixed_total;
                let key = fine_bin(protein_mass);
                let scoring_graph = graph_cache.entry(key).or_insert_with(|| {
                    let n = binning.num_bins();
                    let mut prefix = vec![None; n];
                    let mut suffix = vec![None; n];
                    for &peak in peaks.iter() {
                        let mass = peak.mass;
                        if mass < query.mass_min || mass > query.mass_max {
                            continue;
                        }
                        if let Some(b) = binning.bin_number(mass - prefix_offset) {
                            prefix[b] = Some(1.0);
                        }
                        if let Some(b) = binning.bin_number(protein_mass - (mass - suffix_offset)) {
                            suffix[b] = Some(1.0);
                        }
                    }
                    MassBinScoringGraph::build(&binning, &residue_masses, prefix, suffix)
                });
                let mut gf = GeneratingFunction::new(scoring_graph);
                candidate.spectral_evalue = gf.spectral_e_value(candidate.score);
            }
        } else {
            notes.push("spectral E-value skipped: candidate mass above 50 kDa".to_string());
        }
    }
    if candidates
        .first()
        .and_then(|c| c.spectral_evalue)
        .is_some_and(|ev| ev > 1.0)
    {
        notes.push(
            "spectral E-values exceed 1 (non-discriminating): at intact scale the              counting observable cannot separate the true ladder from the mass-              consistent null space — IP's discriminative power comes from the              envelope-fit terms (Corr/Dist), absent on msalign. Treat the              ranking as the signal; E-values become meaningful once raw-mzML              envelope fits supply those terms.".to_string(),
        );
    }
    let tied_groups = candidates
        .windows(2)
        .filter(|w| (w[0].score - w[1].score).abs() < 1e-6)
        .count();
    if tied_groups > 0 {
        notes.push(format!(
            "{tied_groups} adjacent score ties: positions inside tied groups are              the canonical first-best path, not evidence-discriminated              placements (IP tie semantics); treat them as ambiguous until              discriminating ions separate the candidates"
        ));
    }
    // NOTE: the maxCandidates VIEW cap is applied by the caller AFTER the
    // discovery ranking — cutting before the localization measurement is
    // the E53 reproduction failure (the truth model sat at composition
    // rank 14 under a cap of 10 and never got measured).
    Ok(GraphSearchResult {
        sequence: query.sequence.clone(),
        observed_shift: query.observed_shift,
        shift_frame: "toppic-unexpected-after-fixed-forms".to_string(),
        scans: query.scans.clone(),
        ions: ions.iter().map(|i| format!("{i:?}")).collect(),
        candidates,
        enumerated_combinations: params.num_combinations(),
        evaluated_combinations: evaluated,
        max_candidates: query.max_candidates,
        notes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn auto_ion_types_follow_activation_and_fail_when_unknown() {
        assert_eq!(resolve_ion_types("auto", Some("ETD")).unwrap(), "c,z");
        assert_eq!(resolve_ion_types("auto", Some("HCD")).unwrap(), "b,y");
        assert_eq!(resolve_ion_types("C,Z", None).unwrap(), "C,Z");
        assert!(resolve_ion_types("auto", None).is_err());
    }

    #[test]
    fn measured_proforma_carries_the_pins_and_replaces_the_shift_bracket() {
        let source = "[Acetyl]-PEPTIDEK[+121.972]";
        let measured = vec![
            ("Phospho".to_owned(), "3".to_owned()),
            ("Methyl".to_owned(), "8".to_owned()),
        ];
        let mut notes = Vec::new();
        let rendered = render_measured_proforma(Some(source), &measured, "cand", &mut notes)
            .unwrap()
            .expect("all pins resolved");
        assert_eq!(rendered, "[Acetyl]-PEP[Phospho]TIDEK[Methyl]");
        assert!(!rendered.contains("+121.972"), "shift bracket replaced");
        // Round-trip: the parsed echo must land each tag on its own residue.
        let parsed = crate::proforma::parse_proforma(&rendered, &[]).unwrap();
        let phospho = parsed
            .modifications
            .iter()
            .find(|m| m.tag == "Phospho")
            .unwrap();
        assert_eq!(phospho.residue, Some('P'));
        let methyl = parsed
            .modifications
            .iter()
            .find(|m| m.tag == "Methyl")
            .unwrap();
        assert_eq!(methyl.residue, Some('K'));
    }

    #[test]
    fn measured_proforma_falls_back_when_any_component_is_unresolved() {
        let source = "[Acetyl]-PEPTIDEK[+121.972]";
        let measured = vec![
            ("Phospho".to_owned(), "3".to_owned()),
            ("Methyl".to_owned(), "unresolved".to_owned()),
        ];
        let mut notes = Vec::new();
        let rendered =
            render_measured_proforma(Some(source), &measured, "cand", &mut notes).unwrap();
        assert!(
            rendered.is_none(),
            "no mass-incomplete model may be rendered"
        );
        assert!(
            notes
                .iter()
                .any(|n| n.contains("keeps the DP tie-canonical")),
            "fallback is noted: {notes:?}"
        );
    }

    #[test]
    fn measured_proforma_round_trip_failure_is_a_hard_contract_error() {
        // A tag name outside the chemistry tables is noted, not fatal; but a
        // placement that does not survive parsing is a contract breach. Force
        // the latter by pinning past the sequence end (renderer would place
        // nothing) — the rendered string then lacks the tag and the round
        // trip must fail loudly.
        let source = "[Acetyl]-PEPTIDEK[+121.972]";
        let measured = vec![("Phospho".to_owned(), "99".to_owned())];
        let mut notes = Vec::new();
        let result = render_measured_proforma(Some(source), &measured, "cand", &mut notes);
        assert!(result.is_err(), "out-of-range pin must be a hard error");
    }

    fn methyl_k() -> SearchModification {
        SearchModification {
            name: "Methyl".to_string(),
            composition: GraphComposition::from_delta(14.015_650_064),
            residues: vec!['K'],
            location: SequenceLocation::Everywhere,
            positions: None,
        }
    }

    #[test]
    fn residue_composition_k_is_c6h12n2o() {
        let k = residue_composition('K').unwrap();
        assert_eq!((k.c, k.h, k.n, k.o, k.s), (6, 12, 2, 1, 0));
        assert!((k.mass() - 128.094_963).abs() < 1e-5);
    }

    #[test]
    fn combination_space_and_transitions() {
        let params = ModificationParams::new(vec![methyl_k()], 2, None).unwrap();
        // multisets of size <= 2 over one mod: {}, {Me}, {Me,Me}
        assert_eq!(params.num_combinations(), 3);
        let empty = params.combination(0);
        assert!(empty.is_empty());
        let one = params.index_of_key(&[0]).unwrap();
        let two = params.index_of_key(&[0, 0]).unwrap();
        assert_eq!(params.transition_index(0, 0), Some(one));
        assert_eq!(params.transition_index(one, 0), Some(two));
        assert_eq!(params.transition_index(two, 0), None); // max_mods reached
        assert_eq!(params.modification_between(0, one), Some(0));
    }

    #[test]
    fn shift_appears_at_most_once_and_resolves_against_observed() {
        let params =
            ModificationParams::new(vec![methyl_k(), shift_mod(121.98)], 3, Some(121.98)).unwrap();
        // {} {Me} {Me,Me} {Me,Me,Me} {SHIFT} {Me,SHIFT} {Me,Me,SHIFT}
        assert_eq!(params.num_combinations(), 7);
        let me_shift = params
            .index_of_key(&{
                let mut v = vec![0u16, 1u16];
                v.sort_unstable();
                v
            })
            .unwrap();
        // SHIFT contribution = observed − Σ(other members)
        let comp = params.combination_composition(me_shift);
        let expected = 14.015_650_064 + (121.98 - 14.015_650_064);
        assert!((comp.mass() - expected).abs() < 1e-9, "got {}", comp.mass());
        // transition adding SHIFT twice is absent
        let shift = params.index_of_key(&[1]).unwrap();
        assert_eq!(params.transition_index(shift, 1), None);
    }

    fn shift_mod(delta: f64) -> SearchModification {
        SearchModification {
            name: SHIFT_NAME.to_string(),
            composition: GraphComposition::from_delta(delta),
            residues: vec!['*'],
            location: SequenceLocation::Everywhere,
            positions: None,
        }
    }

    #[test]
    fn graph_layers_collapse_arrangements_into_combination_nodes() {
        // 2 K sites, max 1 Methyl: top layer holds {} and {Me} — the two
        // PLACEMENTS collapse into one {Me} node with two paths.
        let params = ModificationParams::new(vec![methyl_k()], 1, None).unwrap();
        let graph = SequenceGraph::new(&params, "GGKGGK", &[]).unwrap();
        let top = graph.top_nodes();
        assert_eq!(top.len(), 2);
        let counts = graph.num_placements();
        // node order = combination table order: {} first, {Me} second
        assert_eq!(counts[0], 1.0);
        assert_eq!(counts[1], 2.0);
    }

    #[test]
    fn dp_places_modification_where_the_peak_demands() {
        // Sequence GGGKG; peaks engineered so the c3 ion matches ONLY when
        // residue 4 (K) carries the modification.
        let params = ModificationParams::new(vec![methyl_k()], 1, None).unwrap();
        let mut graph = SequenceGraph::new(&params, "GGGKG", &[]).unwrap();
        // find the {Me} sink node
        let top = graph.top_nodes();
        let me_node = top
            .iter()
            .find(|(_, comb)| !params.combination(*comb).is_empty())
            .map(|(node, _)| *node as usize)
            .unwrap();
        graph.set_sink(me_node).unwrap();
        // prefix(GGG) = G+G+G residues; c3 = prefix + BaseIon::C.offset().1 (+ Methyl if
        // placed among the first 3 — impossible, K is residue 4)
        let ggg = residue_composition('G').unwrap();
        let prefix3 = ggg.add(&ggg).add(&ggg);
        // c3 covers residues 1..3; a modification at residue 4 must NOT move
        // it. c4 (residues 1..4 + methyl) is the discriminating ion. Peaks
        // are deconvoluted NEUTRAL masses (msalign convention).
        let k = residue_composition('K').unwrap();
        let c4_me = prefix3
            .add(&k)
            .add(&GraphComposition::from_delta(14.015_650_064));
        let peaks = vec![
            peak(prefix3.mass() + BaseIon::C.offset().1, 100.0), // plain c3
            peak(c4_me.mass() + BaseIon::C.offset().1, 500.0),   // c4 + Methyl
        ];
        let mut scorer = BinCountScorer::new(&peaks, vec![BaseIon::C], 5.0, 50.0, 100_000.0);
        let (score, mods) = graph.get_fragment_score_and_modifications(&mut scorer);
        assert!(score > 0.0);
        assert_eq!(
            mods.iter()
                .find(|(name, _)| name == "Methyl")
                .map(|(_, pos)| pos.clone()),
            Some("4".to_string()),
            "Methyl must land on K4, got {mods:?}"
        );
    }

    #[test]
    fn free_shift_lands_on_the_discriminating_position() {
        // Same construction with the modification expressed as a free SHIFT.
        let params =
            ModificationParams::new(vec![shift_mod(14.015_650_064)], 1, Some(14.015_650_064))
                .unwrap();
        let mut graph = SequenceGraph::new(&params, "GGGKG", &[]).unwrap();
        let top = graph.top_nodes();
        let shift_node = top
            .iter()
            .find(|(_, comb)| !params.combination(*comb).is_empty())
            .map(|(node, _)| *node as usize)
            .unwrap();
        graph.set_sink(shift_node).unwrap();
        let ggg = residue_composition('G').unwrap();
        let prefix3 = ggg.add(&ggg).add(&ggg);
        let k = residue_composition('K').unwrap();
        let c4_me = prefix3
            .add(&k)
            .add(&GraphComposition::from_delta(14.015_650_064));
        let peaks = vec![
            peak(prefix3.mass() + BaseIon::C.offset().1, 100.0),
            peak(c4_me.mass() + BaseIon::C.offset().1, 500.0),
        ];
        let mut scorer = BinCountScorer::new(&peaks, vec![BaseIon::C], 5.0, 50.0, 100_000.0);
        let (_score, mods) = graph.get_fragment_score_and_modifications(&mut scorer);
        assert_eq!(
            mods.iter()
                .find(|(name, _)| name == SHIFT_NAME)
                .map(|(_, pos)| pos.clone()),
            Some("4".to_string()),
            "SHIFT must land on residue 4, got {mods:?}"
        );
    }

    #[test]
    fn positional_fixed_mods_bake_into_the_layer() {
        let params = ModificationParams::new(vec![], 1, None).unwrap();
        let fixed = vec![PositionalModification {
            position: "nTerm".to_string(),
            delta_mass: 42.010_565,
        }];
        let mut graph = SequenceGraph::new(&params, "GGK", &fixed).unwrap();
        let node = graph.top_nodes()[0].0 as usize;
        graph.set_sink(node).unwrap();
        let mass = graph.sink_composition_with_h2o().unwrap().mass();
        let base = residue_composition('G').unwrap();
        let expected = base
            .add(&base)
            .add(&residue_composition('K').unwrap())
            .mass()
            + 42.010_565
            + crate::constants::H2O;
        // the graph accumulates in reverse residue order; association
        // rounding differs at ~1e-8 Da (bin width is 3.6 mDa)
        assert!(
            (mass - expected).abs() < 1e-6,
            "mass {mass} vs expected {expected}"
        );
        let (_score, mods) = {
            struct Zero;
            impl GraphScorer for Zero {
                fn fragment_score(&mut self, _: &GraphComposition, _: &GraphComposition) -> f64 {
                    0.0
                }
            }
            let mut z = Zero;
            graph.get_fragment_score_and_modifications(&mut z)
        };
        assert!(mods
            .iter()
            .any(|(name, pos)| name.starts_with("fixed(") && pos == "nTerm"));
    }

    #[test]
    fn cleave_n_term_drops_the_leading_residue() {
        let params = ModificationParams::new(vec![], 1, None).unwrap();
        let mut graph = SequenceGraph::new(&params, "GGKGG", &[]).unwrap();
        graph.cleave_n_term().unwrap();
        assert_eq!(graph.num_nterm_cleavages(), 1);
        let node = graph.top_nodes()[0].0 as usize;
        graph.set_sink(node).unwrap();
        let mass = graph.sink_composition_with_h2o().unwrap().mass();
        let g = residue_composition('G').unwrap();
        // residues 2..5 of "GGKGG" are G,K,G,G
        let expected = g
            .add(&residue_composition('K').unwrap())
            .add(&g)
            .add(&g)
            .mass()
            + crate::constants::H2O;
        assert!(
            (mass - expected).abs() < 1e-6,
            "mass {mass} vs expected {expected}"
        );
    }

    #[test]
    fn weighted_scorer_counts_negatively_and_rebukes_lone_hits() {
        // One prefix hit only, weakly supported by TopFD (score 0.3):
        // the negative Count weight must dominate — IP's discriminant
        // penalizes lone non-complementary matches that quality cannot repay.
        let g = residue_composition('G').unwrap();
        let prefix = g.add(&g);
        let peaks = vec![peak_with_score(
            prefix.mass() + BaseIon::C.offset().1,
            100.0,
            0.3,
        )];
        let mut scorer =
            WeightedBinScorer::new(&peaks, vec![BaseIon::C, BaseIon::Z], 10.0, 50.0, 1e5);
        let score = scorer.fragment_score(&prefix, &GraphComposition::zero());
        assert!(
            score < 0.0,
            "weakly-backed lone hit must score below zero, got {score}"
        );
        // Strong TopFD deconvolution quality (score 1.0) repays the Count
        // price through the IP peak-quality coefficient.
        let peaks = vec![peak(prefix.mass() + BaseIon::C.offset().1, 100.0)];
        let mut scorer =
            WeightedBinScorer::new(&peaks, vec![BaseIon::C, BaseIon::Z], 101.0, 50.0, 1e5);
        let score = scorer.fragment_score(&prefix, &GraphComposition::zero());
        assert!(
            score > 0.0,
            "envelope-backed hit must be repaid, got {score}"
        );
    }

    #[test]
    fn mass_binning_reaches_residue_sums() {
        let g = residue_composition('G').unwrap().mass();
        let k = residue_composition('K').unwrap().mass();
        let binning = MassBinning::build(&[g, k], 300.0, 27);
        assert!(binning.num_bins() > 1);
        // a residue sum is reachable and maps to a bin
        assert!(binning.bin_number(g + g + k).is_some());
    }

    #[test]
    fn generating_function_e_value_is_monotone() {
        let g = residue_composition('G').unwrap().mass();
        let binning = MassBinning::build(&[g], 300.0, 27);
        // one peak landing on a prefix ion mass inside the graph
        let mut prefix = vec![None; binning.num_bins()];
        let suffix = vec![None; binning.num_bins()];
        for (i, m) in binning.dense_mass.iter().enumerate() {
            if (*m - g * 2.0).abs() < 0.01 {
                prefix[i] = Some(1.0);
            }
        }
        let graph = MassBinScoringGraph::build(&binning, &[g], prefix, suffix);
        let mut gf = GeneratingFunction::new(&graph);
        let low = gf.spectral_e_value(-100.0).expect("distribution exists");
        let mut gf2 = GeneratingFunction::new(&graph);
        let high = gf2.spectral_e_value(0.5).expect("distribution exists");
        assert!(high <= low, "E(high)={high} must be <= E(low)={low}");
        // IP `ScoringGraphEdge.Weight = 0.05`: every convolution discounts a
        // path by 0.05 per edge, so the total mass is strictly below 1 for
        // any nonempty graph (the old `>= 1.0` pin held only under the
        // erroneous weight-1.0 port).
        assert!(
            low > 0.0 && low < 1.0,
            "0.05/edge discount applied, got {low}"
        );
    }

    #[test]
    fn alternatives_report_tied_placements() {
        // No discriminating peaks: every placement of Methyl is optimal, so
        // both K3 and K6 arrangements must be enumerated.
        let params = ModificationParams::new(vec![methyl_k()], 1, None).unwrap();
        let mut graph = SequenceGraph::new(&params, "GGKGGK", &[]).unwrap();
        let top = graph.top_nodes();
        let me_node = top
            .iter()
            .find(|(_, comb)| !params.combination(*comb).is_empty())
            .map(|(node, _)| *node as usize)
            .unwrap();
        graph.set_sink(me_node).unwrap();
        let peaks: Vec<DeconvolutedPeak> = Vec::new();
        let mut scorer = BinCountScorer::new(&peaks, vec![BaseIon::C], 5.0, 50.0, 1e5);
        let (_score, alts, capped) = graph.optimal_placements(&mut scorer, 6);
        assert!(!capped);
        assert_eq!(alts.len(), 2, "both K placements must appear: {alts:?}");
        let positions: Vec<Vec<String>> = alts
            .iter()
            .map(|a| a.iter().map(|(_, p)| p.clone()).collect())
            .collect();
        assert!(positions.contains(&vec!["3".to_string()]));
        assert!(positions.contains(&vec!["6".to_string()]));
    }

    #[test]
    fn graph_search_end_to_end_on_synthetic_msalign() {
        // c3 plain + c4+Methyl peaks in msalign text; the search must rank
        // the Methyl@K candidate above the free SHIFT at another position.
        let g = residue_composition('G').unwrap();
        let prefix3 = g.add(&g).add(&g);
        let k = residue_composition('K').unwrap();
        let c4_me = prefix3
            .add(&k)
            .add(&GraphComposition::from_delta(14.015_650_064));
        let text = format!(
            "BEGIN IONS\nSCANS=1\nPRECURSOR_MASS=600.0\n{} {}\n{} {}\nEND IONS\n",
            prefix3.mass() + BaseIon::C.offset().1,
            100.0,
            c4_me.mass() + BaseIon::C.offset().1,
            500.0
        );
        let query = GraphSearchQuery {
            proforma: None,
            file: "synthetic.msalign".to_string(),
            scans: vec!["1".to_string()],
            sequence: "GGGKG".to_string(),
            fixed_modifications: vec![],
            library: vec![LibrarySpec {
                name: "Methyl".to_string(),
                delta_mass: 14.015_650_064,
                residues: vec!['K'],
                location: None,
                positions: None,
            }],
            max_components: 1,
            observed_shift: 14.015_650_064,
            no_free_shift: true,
            tolerance_ppm: 5.0,
            ion_types: "c".to_string(),
            mass_min: 50.0,
            mass_max: 100_000.0,
            scorer: "count".to_string(),
            overshoot_window_da: 0.06,
            max_candidates: 10,
        };
        let result = graph_search_msalign(text.as_bytes(), &query).unwrap();
        assert_eq!(result.candidates.len(), 1);
        let top = &result.candidates[0];
        assert_eq!(top.combination, "Methyl");
        assert!(top
            .modifications
            .iter()
            .any(|(name, pos)| name == "Methyl" && pos == "4"));
    }

    #[test]
    fn graph_search_keeps_unclosed_residual_as_ambiguous_shift() {
        // The library explains the methyl component at K4. No fragments
        // distinguish the remaining +1 Da residual between K4 and K6, so it
        // must survive as a named SHIFT with both placements, not be silently
        // absorbed into Methyl or forced to one site for mass closure.
        let g = residue_composition('G').unwrap();
        let k = residue_composition('K').unwrap();
        let c3 = g.add(&g).add(&g);
        let c4_me = c3
            .add(&k)
            .add(&GraphComposition::from_delta(14.015_650_064));
        let text = format!(
            "BEGIN IONS\nSCANS=1\nACTIVATION=ETD\n{} 100\n{} 500\nEND IONS\n",
            c3.mass() + BaseIon::C.offset().1,
            c4_me.mass() + BaseIon::C.offset().1,
        );
        let query = GraphSearchQuery {
            proforma: None,
            file: "synthetic.msalign".into(),
            scans: vec!["1".into()],
            sequence: "GGGKGK".into(),
            fixed_modifications: vec![],
            library: vec![LibrarySpec {
                name: "Methyl".into(),
                delta_mass: 14.015_650_064,
                residues: vec!['K'],
                location: None,
                positions: None,
            }],
            max_components: 2,
            observed_shift: 15.015_650_064,
            no_free_shift: false,
            tolerance_ppm: 5.0,
            ion_types: "auto".into(),
            mass_min: 50.0,
            mass_max: 100_000.0,
            scorer: "count".into(),
            overshoot_window_da: 0.06,
            max_candidates: 10,
        };
        let result = graph_search_msalign(text.as_bytes(), &query).unwrap();
        assert_eq!(result.ions, vec!["C", "Z"]);
        let model = result
            .candidates
            .iter()
            .find(|candidate| {
                candidate.components.iter().any(|m| m == "Methyl")
                    && candidate.components.iter().any(|m| m == SHIFT_NAME)
            })
            .expect("Methyl + residual SHIFT candidate");
        assert!((model.free_mass.unwrap() - 1.0).abs() < 1e-9);
        assert!(model.alternative_placements.is_some());
        assert!(model.alternatives.iter().any(|placement| placement
            .iter()
            .any(|(name, position)| name == SHIFT_NAME && position == "4")));
        assert!(model.alternatives.iter().any(|placement| placement
            .iter()
            .any(|(name, position)| name == SHIFT_NAME && position == "6")));
    }

    #[test]
    fn free_shift_candidate_carries_its_own_localization() {
        // One-flow contract: the search output must carry the free mass's
        // own localization (site + evidence + per-component cliffs), not
        // just the DP placement. Plain c1-c3 plus c4/c5/c6 and two z ions
        // all shifted by exactly the free mass: it must pin uniquely to
        // residue 4 (the K) with evidence over the p<0.01 bar.
        let g = residue_composition('G').unwrap();
        let k = residue_composition('K').unwrap();
        let free = 100.0;
        // sequence GGGKGGGG: prefix masses for cuts 1..6, suffix z at cuts 2,3
        let mut prefix = vec![GraphComposition::zero()];
        for comp in [&g, &g, &g, &k, &g, &g, &g, &g] {
            prefix.push(prefix.last().unwrap().add(comp));
        }
        let total = prefix[8].clone();
        let suffix = |cut: usize| total.sub(&prefix[cut]);
        let peaks = [
            (prefix[1].mass() + BaseIon::C.offset().1, 100.0),
            (prefix[2].mass() + BaseIon::C.offset().1, 100.0),
            (prefix[3].mass() + BaseIon::C.offset().1, 100.0),
            (prefix[4].mass() + BaseIon::C.offset().1 + free, 500.0),
            (prefix[5].mass() + BaseIon::C.offset().1 + free, 400.0),
            (prefix[6].mass() + BaseIon::C.offset().1 + free, 300.0),
            (suffix(2).mass() + BaseIon::Z.offset().1 + free, 350.0),
            (suffix(3).mass() + BaseIon::Z.offset().1 + free, 250.0),
        ];
        let body = peaks
            .iter()
            .map(|(m, i)| format!("{m}\t{i}\t1\t1.0"))
            .collect::<Vec<_>>()
            .join("\n");
        let text = format!("BEGIN IONS\nSCANS=1\nACTIVATION=ETD\n{body}\nEND IONS\n");
        let query = GraphSearchQuery {
            proforma: None,
            file: "synthetic.msalign".into(),
            scans: vec!["1".into()],
            sequence: "GGGKGGGG".into(),
            fixed_modifications: vec![],
            library: vec![LibrarySpec {
                name: "Methyl".into(),
                delta_mass: 14.015_650_064,
                residues: vec!['K'],
                location: None,
                positions: None,
            }],
            max_components: 2,
            observed_shift: 114.015_650_064,
            no_free_shift: false,
            tolerance_ppm: 5.0,
            ion_types: "c,z".into(),
            mass_min: 50.0,
            mass_max: 100_000.0,
            scorer: "count".into(),
            overshoot_window_da: 0.06,
            max_candidates: 10,
        };
        let result = graph_search_msalign(text.as_bytes(), &query).unwrap();
        let model = result
            .candidates
            .iter()
            .find(|candidate| {
                candidate.components.iter().any(|m| m == "Methyl")
                    && candidate.components.iter().any(|m| m == SHIFT_NAME)
            })
            .expect("Methyl + free SHIFT candidate");
        assert!((model.free_mass.unwrap() - free).abs() < 1e-6);
        let evidence = model
            .shift_evidence_log10
            .expect("free SHIFT must carry localization evidence");
        assert!(
            evidence > 2.0,
            "evidence {evidence} must clear the p<0.01 bar"
        );
        assert_eq!(
            model.shift_best_site,
            Some(4),
            "free mass pins to residue 4"
        );
        assert_eq!(model.shift_site_interval, Some((4, 4)));
        // the Methyl component has no peaks of its own: its cliff must be
        // honestly unresolved (interval (0,0)), never pinned to the DP tie
        assert!(model
            .component_cliffs
            .iter()
            .any(|c| c.name == "Methyl" && c.site_interval == (0, 0)));
    }

    #[test]
    fn component_localization_cannot_escape_library_residue_specificity() {
        // The observed c ladder favors a shift at T2, but Bromo is allowed
        // only on Y4. The component sweep must stay inside the same candidate
        // space as the graph instead of rendering the impossible Bromo@T2.
        let sequence = "GTGYGG";
        let bromo = 77.910_512_542_102;
        let residual = 0.01;
        let mut prefix = GraphComposition::zero();
        let mut peaks = Vec::new();
        for (index, residue) in sequence.chars().enumerate() {
            prefix = prefix.add(&residue_composition(residue).unwrap());
            let cut = index + 1;
            if cut < sequence.len() {
                let mass =
                    prefix.mass() + BaseIon::C.offset().1 + if cut >= 2 { bromo } else { 0.0 };
                peaks.push(format!("{mass}\t{}\t1\t1.0", 100 + cut * 10));
            }
        }
        let text = format!(
            "BEGIN IONS\nSCANS=1\nACTIVATION=ETD\n{}\nEND IONS\n",
            peaks.join("\n")
        );
        let query = GraphSearchQuery {
            proforma: Some(format!("{sequence}[+{}]", bromo + residual)),
            file: "synthetic.msalign".into(),
            scans: vec!["1".into()],
            sequence: String::new(),
            fixed_modifications: vec![],
            library: vec![LibrarySpec {
                name: "Bromo".into(),
                delta_mass: bromo,
                residues: vec!['Y'],
                location: None,
                positions: None,
            }],
            max_components: 1,
            observed_shift: 0.0,
            no_free_shift: false,
            tolerance_ppm: 5.0,
            ion_types: "c,z".into(),
            mass_min: 50.0,
            mass_max: 100_000.0,
            scorer: "count".into(),
            overshoot_window_da: 0.06,
            max_candidates: 10,
        };

        let result = graph_search_msalign(text.as_bytes(), &query)
            .expect("an impossible component pin must not abort the whole search");
        let model = result
            .candidates
            .iter()
            .find(|candidate| candidate.components.iter().any(|name| name == "Bromo"))
            .expect("Bromo candidate");
        let cliff = model
            .component_cliffs
            .iter()
            .find(|cliff| cliff.name == "Bromo")
            .expect("Bromo component localization");
        assert!(
            cliff.site_interval == (0, 0) || cliff.site_interval == (4, 4),
            "Bromo localization escaped Y specificity: {:?}",
            cliff.site_interval
        );
    }

    #[test]
    fn max_candidates_caps_the_output_view_not_the_search() {
        // A free-shift search produces many candidates; the cap must cut
        // the OUTPUT only, with the truncation count in the notes and the
        // full evaluation count still reported.
        let g = residue_composition('G').unwrap();
        let prefix3 = g.add(&gg()).add(&g);
        let text = format!(
            "BEGIN IONS
SCANS=1
{}	100.0	1	1.0
END IONS
",
            prefix3.mass() + BaseIon::C.offset().1
        );
        let query = GraphSearchQuery {
            proforma: None,
            file: "x.msalign".into(),
            scans: vec!["1".into()],
            sequence: "GGGKG".into(),
            fixed_modifications: vec![],
            library: vec![LibrarySpec {
                name: "Methyl".into(),
                delta_mass: 14.015_650_064,
                residues: vec!['K'],
                location: None,
                positions: None,
            }],
            max_components: 1,
            observed_shift: 14.015_650_064,
            no_free_shift: false,
            tolerance_ppm: 5.0,
            ion_types: "c".into(),
            mass_min: 50.0,
            mass_max: 100_000.0,
            scorer: "count".into(),
            overshoot_window_da: 0.06,
            max_candidates: 1,
        };
        let result = graph_search_msalign(text.as_bytes(), &query).unwrap();
        assert_eq!(result.candidates.len(), 1);
        assert!(result.evaluated_combinations >= 2, "search stays complete");
        assert!(result
            .notes
            .iter()
            .any(|n| n.contains("truncated by maxCandidates=1")));
    }

    #[test]
    fn positional_prior_restricts_placement() {
        // Methyl(K) constrained to position 4: the DP must place it there
        // even when another K exists, and the constrained-out arrangement
        // must not appear among alternatives.
        let params = ModificationParams::new(
            vec![SearchModification {
                name: "Methyl".into(),
                composition: GraphComposition::from_delta(14.015_650_064),
                residues: vec!['K'],
                location: SequenceLocation::Everywhere,
                positions: Some(vec!["3".to_string()]),
            }],
            1,
            None,
        )
        .unwrap();
        let mut graph = SequenceGraph::new(&params, "GGKGGK", &[]).unwrap();
        let top = graph.top_nodes();
        let me_node = top
            .iter()
            .find(|(_, comb)| !params.combination(*comb).is_empty())
            .map(|(node, _)| *node as usize)
            .unwrap();
        graph.set_sink(me_node).unwrap();
        let peaks: Vec<DeconvolutedPeak> = Vec::new();
        let mut scorer = BinCountScorer::new(&peaks, vec![BaseIon::C], 5.0, 50.0, 1e5);
        let (_s, alts, _c) = graph.optimal_placements(&mut scorer, 6);
        assert!(!alts.is_empty());
        for alt in &alts {
            assert!(
                alt.iter().all(|(_, p)| p == "3"),
                "placement outside the prior (K6 must be excluded) leaked: {alt:?}"
            );
        }
    }

    fn gg() -> GraphComposition {
        residue_composition('G').unwrap()
    }

    fn peak(mass: f64, intensity: f32) -> DeconvolutedPeak {
        peak_with_score(mass, intensity, 1.0)
    }

    fn peak_with_score(mass: f64, intensity: f32, topfd_score: f32) -> DeconvolutedPeak {
        DeconvolutedPeak {
            mass,
            intensity,
            charge: Some(1),
            score: Some(topfd_score),
        }
    }
}
