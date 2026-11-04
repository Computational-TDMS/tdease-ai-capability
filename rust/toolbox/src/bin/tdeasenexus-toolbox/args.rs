//! Command-line surface. `--help` renders as the interface documentation.

use clap::{Args as ClapArgs, Parser, Subcommand};

#[derive(Debug, Clone, Copy, clap::ValueEnum)]
pub enum ToleranceUnit {
    Ppm,
    Da,
}

impl ToleranceUnit {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Ppm => "ppm",
            Self::Da => "da",
        }
    }
}

#[derive(Parser)]
#[command(
    name = "tdeasenexus-toolbox",
    version,
    about = "Native adapter for the TDEaseNexus mass-spec core (mzML, isotopes, ProFormA).",
    long_about = "One JSON object on stdout per call; code/message errors on stderr.\n\
                  Input files are read-only; output defaults to stdout. See\n\
                  docs/design/tdeasenexus-toolbox-native-adapter.md for the contract."
)]
pub struct Args {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Subcommand)]
pub enum Command {
    /// Index an mzML and report a file summary.
    Inspect(InspectArgs),
    /// Parse a ProFormA into exact composition, mass, and modification trace.
    ParseProforma(ParseProformaArgs),
    /// Generate the BRAIN theoretical isotope envelope for a composition.
    Theory(TheoryArgs),
    /// Aggregate the observed MS1 isotope envelope across an RT window.
    Envelope(EnvelopeArgs),
    /// Aggregate several observed MS1 isotope envelopes across RT windows in
    /// one mzML load (the multiplicity-ladder surface).
    EnvelopeWindows(EnvelopeWindowsArgs),
    /// Fit observed peaks against an anchored theoretical envelope.
    FitEnvelope(FitEnvelopeArgs),
    /// Calibrate a mass shift from the MS1 isotope envelope of a candidate.
    RefineShift(RefineShiftArgs),
    /// Calibrate many mass shifts from one mzML load.
    RefineShiftBatch(RefineShiftBatchArgs),
    /// Query the mzML in one of several modes.
    Spectrum(SpectrumArgs),
    /// Run a SELECT statement over an mzML — the condition+range query
    /// surface (TIC/XIC extraction, survey listings, fragment evidence).
    Select(SelectArgs),
    /// Measure complete proteoform models on selected or pooled MS2
    /// observations (raw mzML, deconvoluted msalign, or saved spectrum JSON).
    Fragment(FragmentArgs),
    /// Selectively extract observed raw peaks inside the theoretical isotope
    /// envelope of chosen fragments (`--ions y3,y2` or `b4-b7`, optional
    /// per-ion charge `y3:2`), paired with the BRAIN theoretical sticks —
    /// the compact per-request counterpart of a full envelope export, for
    /// mirror-comparison visualization.
    EnvelopeRegions(EnvelopeRegionsArgs),
    /// Generate the theoretical fragment ions for a ProFormA (the spectrum_utils
    /// `get_theoretical_fragments` analogue) — a pure-theory command that needs
    /// no observed spectrum. Output: JSON envelope by default, or a standard
    /// mzTab-M molecular table with `--out-format mztab`.
    FragmentTheoretical(FragmentTheoreticalArgs),
    /// Get one UniMod modification by id or exact title (dual-source read).
    UnimodGet(UnimodGetArgs),
    /// Search UniMod by fuzzy name and/or by mass with composition filters.
    UnimodSearch(UnimodSearchArgs),
    /// Search many UniMod delta masses while loading the database once.
    UnimodBatchSearch(UnimodBatchSearchArgs),
    /// Solve many UniMod combination masses while loading the database once.
    UnimodBatchSolve(UnimodBatchSolveArgs),
    /// Solve multi-modification mass combinations with ppm-ranked candidates.
    UnimodSolve(UnimodSolveArgs),
    /// List the amino-acid residue table.
    UnimodResidues(UnimodResiduesArgs),
    /// Add a user modification (`USER:n`) to the writable user store.
    UnimodAdd(UnimodAddArgs),
    /// Batch BRAIN isotope envelopes for many compositions (one call).
    IsotopicBatch(IsotopicBatchArgs),
    /// Validate an MS2 precursor against its nearest MS1 isotope envelope
    /// (IsoForma/Degnan 2024: match the precursor's theoretical isotope
    /// distribution to the closest MS1, require >= N matched peaks).
    ValidatePrecursor(ValidatePrecursorArgs),
    /// Resolve an observed PrSM mass shift into positioned ProForma
    /// components on deconvoluted TopFD msalign scans with the
    /// sequence-graph engine (the Informed-Proteomics core): one DP over
    /// (position x combination) states. The default library is builtin;
    /// --library none localizes a free residual alone, while builtin|<json>
    /// also proposes decomposition over the modification library.
    ResolveShift(ResolveShiftArgs),
    /// Print the JSON Schema for a command's output envelope.
    Schema(SchemaArgs),
}

/// Sequence-graph engine parameter bundle of `resolve-shift` (constructed
/// by its dispatcher; not a CLI subcommand).
pub struct GraphSearchArgs {
    /// TopFD `_ms2.msalign` (deconvoluted spectra — the sole deep-dive input).
    pub file: String,
    /// MS2 scan numbers (msalign SCANS= field), comma-separated; sibling
    /// scans of one precursor pool their peaks.
    pub scan: Vec<String>,
    /// TopPIC ProFormA from the source record:
    /// derives sequence/observed-shift and mods from the record itself
    /// (mass-only bracket = the shift; named tags = fixed forms).
    pub proforma: Option<String>,
    /// Protein sequence, one-letter uppercase (omit with --proforma).
    pub sequence: Option<String>,
    /// The observed unexpected shift (Da) the library + free mass must
    /// explain (TopPIC shift frame: after fixed forms; omit with
    /// --proforma).
    pub observed_shift: Option<f64>,
    /// Modification library JSON. Omit for the built-in common PTM library;
    /// provide this only to replace that library. Entries may omit deltaMass
    /// when `name` is a built-in ProForma modification name.
    pub library: Option<String>,
    /// Max total components (library + SHIFT) per combination.
    pub max_components: usize,
    /// Disable the free SHIFT (require the library alone to close the shift).
    pub no_free_shift: bool,
    /// Positional fixed modifications JSON (TopPIC ProForma positions),
    /// name-first: [{"position":"nTerm","name":"Acetyl"}]. The name is
    /// resolved by the chemistry authority (ProForma/UniMod); a literal
    /// deltaMass is accepted only for chemistry absent from all databases.
    pub mods: Option<String>,
    /// Fragment ion series; auto derives ETD/ECD c/z or HCD/CID b/y from
    /// the msalign ACTIVATION field.
    pub ion_types: String,
    /// Fragment match tolerance (ppm).
    pub tolerance_ppm: f64,
    /// Fragment neutral-mass bounds.
    pub mass_min: f64,
    pub mass_max: f64,
    /// Scorer: "weighted" (default; IP weights over count, intensity, TopFD
    /// quality and mass error; unavailable Dist=0) or "count" (diagnostic
    /// IP DeconvScorer, +1 per matched ion).
    pub scorer: String,
    /// Library mass may overshoot the observed shift by at most this window
    /// (Da) and stay evaluable.
    pub overshoot_window_da: f64,
    /// Output cap on the ranked candidate list (search stays complete).
    pub max_candidates: usize,
}

/// Shift-resolution surface: the sequence-graph engine (Informed-Proteomics
/// core) resolving an observed shift into positioned components. The
/// dispatcher constructs the engine bundle (`GraphSearchArgs`); the builtin
/// library is the default, and `none` localizes only the free residual.
#[derive(ClapArgs)]
pub struct ResolveShiftArgs {
    /// TopFD `_ms2.msalign` (read-only) -- the deconvoluted spectra TopPIC
    /// searched; the sole deep-dive input.
    #[arg(value_name = "FILE")]
    pub file: String,
    /// MS2 scan numbers (msalign SCANS= field), comma-separated. Scans may
    /// pool evidence only when they use the same activation method.
    #[arg(long, value_name = "SCANS", value_delimiter = ',')]
    pub scan: Vec<String>,
    /// TopPIC ProForma: derives --sequence, --shift and the fixed forms from
    /// the record itself (mass-only bracket = the shift; named tags = fixed
    /// forms).
    #[arg(long, value_name = "PF")]
    pub proforma: Option<String>,
    /// Protein sequence, one-letter uppercase; alternative to --proforma.
    #[arg(long, value_name = "SEQ")]
    pub sequence: Option<String>,
    /// The observed unexpected delta mass to resolve (Da).
    #[arg(long, value_name = "MASS", allow_hyphen_values = true)]
    pub shift: Option<f64>,
    /// Modification library: `builtin` for the common PTM library (the
    /// default — known chemistry is labeled, the unexplained residual stays
    /// as a positioned free SHIFT), a JSON array supplied as a file path or
    /// inline for an explicit replacement, or `none` for the label-free
    /// single-event localization (free residual SHIFT alone).
    #[arg(long, value_name = "LIB", default_value = "builtin")]
    pub library: Option<String>,
    /// Max total components (library + free shift) per combination
    /// (default 3).
    #[arg(long, default_value_t = 3, value_name = "N")]
    pub max_components: usize,
    /// Forbid the free residual shift: the library alone must close the
    /// observed mass.
    #[arg(long)]
    pub no_free_shift: bool,
    /// Positional fixed modifications JSON, inline or file, name-first
    /// (ProForma/UniMod-resolved): [{"position":"nTerm","name":"Acetyl"}];
    /// a literal deltaMass only for chemistry absent from all databases.
    #[arg(long, value_name = "JSON")]
    pub mods: Option<String>,
    /// Fragment ion series (default auto: derive from the msalign ACTIVATION
    /// field).
    #[arg(long, default_value = "auto", value_name = "IONS")]
    pub ion_types: String,
    /// Fragment match tolerance in ppm (default 10).
    #[arg(long, default_value_t = 10.0, value_name = "PPM")]
    pub tolerance_ppm: f64,
    /// Fragment neutral-mass lower bound.
    #[arg(long, default_value_t = 50.0)]
    pub mass_min: f64,
    /// Fragment neutral-mass upper bound.
    #[arg(long, default_value_t = 100_000.0)]
    pub mass_max: f64,
    /// Scorer: "weighted" (default; IP weights over available msalign peak
    /// evidence) or "count" (diagnostic, +1 per matched ion).
    #[arg(long, default_value = "weighted", value_name = "NAME")]
    pub scorer: String,
    /// Library mass may overshoot the observed shift by at most this
    /// window (Da) and stay evaluable.
    #[arg(long, default_value_t = 0.06, value_name = "DA")]
    pub overshoot_window_da: f64,
    /// Output cap on the ranked candidate list (the search stays complete).
    #[arg(long, default_value_t = 5, value_name = "N")]
    pub max_candidates: usize,
    /// Expand diagnostics for the candidates selected by --max-candidates,
    /// including alternative placements, spectral fields, and notes.
    /// Default stdout is the compact candidate view for interactive use.
    #[arg(long)]
    pub details: bool,
}

#[derive(ClapArgs)]
pub struct SchemaArgs {
    /// The command whose output schema to print (e.g. refine-shift).
    #[arg(value_name = "COMMAND")]
    pub command: String,
}

#[derive(ClapArgs)]
pub struct InspectArgs {
    /// mzML file to index (read-only).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Compute a full-file sha256 instead of the fast head digest.
    #[arg(long)]
    pub digest_full: bool,
    /// Cap the number of spectra indexed.
    #[arg(long)]
    pub max_spectra: Option<usize>,
}

#[derive(ClapArgs)]
pub struct ParseProformaArgs {
    /// ProFormA string, e.g. "ACDM[Oxidation]K".
    #[arg(long)]
    pub proforma: String,
    /// JSON file of external mods: [{"name": "...", "formula": "..."}].
    #[arg(long)]
    pub mods: Option<String>,
    /// Explicit UniMod SQLite path: switches to the dual-source resolver
    /// (same records, explicit source control). Default path is the single
    /// core parser with the embedded ontology.
    #[arg(long)]
    pub db: Option<String>,
    /// User modification store (user-mods.sqlite). A name matching both
    /// sources is an ambiguity error demanding an explicit [USER:<id>].
    #[arg(long)]
    pub user_db: Option<String>,
}

#[derive(ClapArgs)]
pub struct TheoryArgs {
    /// Complete ProFormA seeding the envelope. Every modification must be a
    /// named or formula-bearing tag; bare numeric deltas are rejected because
    /// they do not define an isotope composition.
    #[arg(long, required_unless_present = "models", conflicts_with = "models")]
    pub proforma: Option<String>,
    /// Batch complete models, inline JSON or a file:
    /// [{"id":"a","proforma":"PEPTIDE","charge":2,"anchor":...}, ...].
    /// Row charge/anchor override the shared flags; every row is checked first.
    #[arg(long, value_name = "JSON")]
    pub models: Option<String>,
    /// JSON file of external mods (see parse-proforma).
    #[arg(long)]
    pub mods: Option<String>,
    /// Charge state of the ion.
    #[arg(long, default_value_t = 1)]
    pub charge: i32,
    /// Truncate the envelope after this cumulative intensity (default 0.99).
    #[arg(long, default_value_t = 0.99)]
    pub truncate_after: f64,
    /// Drop peaks below this relative intensity (default 0.001).
    #[arg(long, default_value_t = 0.001)]
    pub ignore_below: f64,
    /// Anchor M0 to this neutral mass. Defaults to the proforma's exact
    /// monoisotopic mass.
    #[arg(long)]
    pub anchor: Option<f64>,
}

#[derive(ClapArgs)]
pub struct EnvelopeArgs {
    /// mzML file to aggregate over (read-only).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Anchor neutral mass for the isotope targets (the candidate mass).
    #[arg(long)]
    pub mass: f64,
    /// Charge state of the ion.
    #[arg(long)]
    pub charge: u8,
    /// RT window start (minutes).
    #[arg(long, default_value_t = -1.0, allow_negative_numbers = true)]
    pub rt_min: f64,
    /// RT window end (minutes).
    #[arg(long, default_value_t = 1.0, allow_negative_numbers = true)]
    pub rt_max: f64,
    /// Extraction tolerance per isotope (ppm).
    #[arg(long, default_value_t = 10.0)]
    pub tolerance_ppm: f64,
    /// Number of isotope targets to aggregate. Zero (default) derives a
    /// bounded extraction window from the ordinary isotope support.
    #[arg(long, default_value_t = 0)]
    pub isotopes: usize,
}

#[derive(ClapArgs)]
pub struct EnvelopeWindowsArgs {
    /// mzML file to aggregate over (read-only).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// JSON array of extraction windows:
    /// [{"label","mass","charge","rtMin","rtMax"}, ...] (or a JSON file path).
    #[arg(long, value_name = "JSON")]
    pub windows: String,
    /// Extraction tolerance per isotope (ppm).
    #[arg(long, default_value_t = 10.0)]
    pub tolerance_ppm: f64,
    /// Number of isotope targets to aggregate. Zero (default) derives a
    /// bounded extraction window from the ordinary isotope support.
    #[arg(long, default_value_t = 0)]
    pub isotopes: usize,
}

#[derive(ClapArgs)]
pub struct FitEnvelopeArgs {
    /// JSON file of observed peaks: [{"mz": ..., "intensity": ...}].
    #[arg(long, required_unless_present = "pairs", conflicts_with = "pairs")]
    pub observed: Option<String>,
    /// JSON file of an anchored theoretical envelope (theory --anchor output).
    #[arg(long, required_unless_present = "pairs", conflicts_with = "pairs")]
    pub theoretical: Option<String>,
    /// Batch fits, inline JSON or a file:
    /// [{"id":"a","observed":"obs.json","theoretical":"theory.json",
    ///   "matchWindowPpm":30}, ...]. Repeated files are loaded once.
    #[arg(long, value_name = "JSON")]
    pub pairs: Option<String>,
    /// The fit's per-peak match window (ppm).
    #[arg(long, default_value_t = 30.0)]
    pub match_window_ppm: f64,
    /// Absolute m/z window floor. Set to 0 for a strict ppm-only fit.
    #[arg(long, default_value_t = tdeasenexus_toolbox::DEFAULT_ENVELOPE_WINDOW_MZ)]
    pub minimum_window_mz: f64,
}

#[derive(ClapArgs)]
pub struct RefineShiftArgs {
    /// mzML file (read-only).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Reuse a previously extracted observed envelope (JSON file) instead of
    /// aggregating over the mzML again.
    #[arg(long)]
    pub observed: Option<String>,
    /// Complete candidate ProFormA, or declared baseline with --mod-mass/--mod-formula.
    #[arg(long)]
    pub proforma: Option<String>,
    /// JSON file of external mods.
    #[arg(long)]
    pub mods: Option<String>,
    /// Unmodified sequence mass (legacy path, averagine composition).
    #[arg(long)]
    pub sequence_mass: Option<f64>,
    /// Unassigned residual mass added to the declared baseline (calibration only).
    #[arg(long)]
    pub mod_mass: Option<f64>,
    /// Candidate modification net formula added to the declared baseline.
    #[arg(long)]
    pub mod_formula: Option<String>,
    /// Measured precursor m/z; the M0 extraction basis is baseline + --prior-shift.
    #[arg(long)]
    pub mz: f64,
    /// Charge state of the ion.
    #[arg(long)]
    pub charge: i32,
    /// RT window start (minutes).
    #[arg(long, default_value_t = -1.0, allow_negative_numbers = true)]
    pub rt_min: f64,
    /// RT window end (minutes).
    #[arg(long, default_value_t = 1.0, allow_negative_numbers = true)]
    pub rt_max: f64,
    /// Aggregation tolerance per isotope (ppm).
    #[arg(long, default_value_t = 10.0)]
    pub tolerance_ppm: f64,
    /// The prior mass shift reported by the search engine, in the
    /// post-declared-forms frame (same frame as refinedShiftAfterDeclared;
    /// TopPIC's unexpected shift is already in this frame).
    #[arg(long, default_value_t = 0.0, allow_negative_numbers = true)]
    pub prior_shift: f64,
    /// The fit's per-peak match window (ppm).
    #[arg(long, default_value_t = 30.0)]
    pub match_window_ppm: f64,
}

#[derive(ClapArgs)]
pub struct SpectrumArgs {
    /// mzML file (read-only).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Query mode.
    #[arg(long, value_enum)]
    pub mode: SpectrumMode,
    /// MS level filter. 0 (the default) matches any level in scans mode —
    /// `--scans N` alone stays level-agnostic; an explicit value filters, and
    /// a scan that exists only at another level lands in missingScans. The
    /// window and XIC modes treat 0 as MS1 (their historical default).
    #[arg(long, default_value_t = 0)]
    pub ms_level: u8,
    /// RT lower bound in minutes (window/sum/scans/product/precursor/xic).
    #[arg(long, allow_negative_numbers = true)]
    pub rt_min: Option<f64>,
    /// RT upper bound in minutes (window/sum/scans/product/precursor/xic).
    #[arg(long, allow_negative_numbers = true)]
    pub rt_max: Option<f64>,
    /// Optional scan limit for window mode; omitted means the full selection.
    #[arg(long)]
    pub max_scans: Option<usize>,
    /// Inclusive native scan-number bounds, intersected with RT and scan ids.
    #[arg(long)]
    pub scan_min: Option<u64>,
    #[arg(long)]
    pub scan_max: Option<u64>,
    /// Binning tolerance in ppm (window/sum).
    #[arg(long, default_value_t = 20.0)]
    pub bin_ppm: f64,
    /// m/z of interest (xic/product/precursor).
    #[arg(long)]
    pub mz: Option<f64>,
    /// Tolerance value (product/precursor).
    #[arg(long)]
    pub tolerance: Option<f64>,
    /// Tolerance unit (product/precursor): ppm (default) or da.
    #[arg(long, value_enum, default_value_t = ToleranceUnit::Ppm)]
    pub tolerance_unit: ToleranceUnit,
    /// Charge filter (precursor).
    #[arg(long)]
    pub charge: Option<i32>,
    /// Comma-separated scan numbers (sum/scans; product restricts its
    /// otherwise whole-file census to these scans).
    #[arg(long)]
    pub scans: Option<String>,
    /// JSON array of scan ids for extraction/aggregation. Avoids command-line length
    /// limits for large batches.
    #[arg(long, conflicts_with = "scans")]
    pub scans_input: Option<String>,
    /// m/z lower bound.
    #[arg(long)]
    pub mz_min: Option<f64>,
    /// m/z upper bound.
    #[arg(long)]
    pub mz_max: Option<f64>,
    /// Minimum peak intensity.
    #[arg(long)]
    pub min_intensity: Option<f32>,
    /// Max peaks to return.
    #[arg(long)]
    pub max_peaks: Option<usize>,
    /// Peak ordering for scans mode: mzAscending (default) | intensityDescending
    /// | intensityCappedMzAscending.
    #[arg(long)]
    pub peak_order: Option<String>,
    /// Output format: json (default) or mgf (scans, interchange for external
    /// tools). Tabular exports go through `spectra.export`, not the kernel.
    #[arg(long, default_value = "json")]
    pub out_format: String,
}

#[derive(Debug, Clone, Copy, clap::ValueEnum)]
pub enum SpectrumMode {
    Window,
    Sum,
    Xic,
    Scans,
    Precursor,
    Product,
}

#[derive(ClapArgs)]
pub struct SelectArgs {
    /// mzML file (read-only).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// The SELECT statement: SELECT spectra FROM mzml WHERE <predicates>
    /// [LIMIT n];
    #[arg(long, value_name = "SQL")]
    pub statement: String,
}

#[derive(Debug, Clone, Copy, clap::ValueEnum)]
pub enum FragmentMode {
    /// Match one hypothesis: the spectrum's observed peaks against the
    /// theoretical fragments of one complete `--proforma` model.
    Match,
    /// Compare one or more complete `--candidate-proformas` models and report
    /// which matched ions distinguish those models.
    Compare,
}

#[derive(ClapArgs)]
pub struct FragmentArgs {
    /// Sum selected centroided MS2 scans before measuring any model set.
    #[arg(long)]
    pub sum_scans: bool,
    /// Ppm tolerance used to merge centroid peaks in a summed MS2 spectrum.
    #[arg(long, default_value_t = 10.0)]
    pub bin_ppm: f64,
    /// Drop pooled MS2 peaks below this percent of the strongest input peak.
    /// IsoForma uses 0.01 percent.
    #[arg(long, default_value_t = 0.01)]
    pub minimum_abundance_pct: f64,
    /// Raw centroided mzML, TopFD msalign, or a saved spectrum observation JSON.
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Query mode: match (one `--proforma`) or compare
    /// (`--candidate-proformas`). One matching engine either way — compare is
    /// the same pass fanned out over candidates plus the diff. Omit --mode
    /// and it is derived from which of the two is present.
    #[arg(long, value_enum)]
    pub mode: Option<FragmentMode>,
    /// Scan ids to measure (comma-separated), a scan number (resolved to nativeID scan=N)
    /// or a full nativeID (e.g. scan=17).
    #[arg(
        long,
        value_name = "SCAN",
        default_value = "",
        conflicts_with = "scans_input"
    )]
    pub scan: String,
    /// JSON array of scan ids/numbers, inline or a file path. Same selection
    /// as --scan, without command-line length limits.
    #[arg(long, value_name = "JSON")]
    pub scans_input: Option<String>,
    /// RT bounds in minutes, intersected with scan ids/range (mzML).
    #[arg(long)]
    pub rt_min: Option<f64>,
    #[arg(long)]
    pub rt_max: Option<f64>,
    /// Inclusive native scan-number bounds (mzML).
    #[arg(long)]
    pub scan_min: Option<u64>,
    #[arg(long)]
    pub scan_max: Option<u64>,
    /// Fragmentation method (ETD/ECD/HCD/CID). When supplied, the core derives
    /// the appropriate ion series; otherwise it reads source metadata.
    #[arg(long, value_name = "METHOD")]
    pub activation: Option<String>,
    /// One complete ProFormA candidate for `--mode match`. It contains the
    /// sequence, terminal processing and every known/proposed modification.
    #[arg(long, value_name = "PROFORMA")]
    pub proforma: Option<String>,
    /// JSON array for `--mode compare`:
    /// [{"id":"phospho","proforma":"[Acetyl]-S[Phospho]GRGK"}, ...].
    /// Every model is independently parsed to exact chemistry before matching.
    #[arg(long, value_name = "JSON")]
    pub candidate_proformas: Option<String>,
    /// Evidence artifact for fragment matching or candidate comparison. Full
    /// match tables are written here; stdout contains only a bounded summary.
    #[arg(long, value_name = "FILE")]
    pub evidence_out: Option<String>,
    /// Comma-separated fragment ion types: A,B,C,X,Y,Z. When omitted, the
    /// kernel infers ETD/ECD -> C,Z and HCD/CID -> B,Y from source metadata;
    /// unknown activation must be explicit.
    #[arg(long, value_name = "IONS")]
    pub ion_types: Option<String>,
    /// Pearson threshold for raw-mzML envelope matching (Spectator default
    /// 0.7). mzML matching is always envelope-level: anchor on the
    /// most-abundant isotope, walk the cluster with a monotone cursor, and
    /// require observed-vs-theoretical Pearson >= this threshold. A
    /// mass-shifted wrong model fails coherently on every fragment crossing
    /// the misassigned site instead of riding satellites. (msalign input is
    /// deconvoluted and matches in neutral-mass space — no envelope stage.)
    #[arg(long, default_value_t = 0.7)]
    pub min_correlation: f64,
    /// Max charge state of fragment ions. Default: the precursor's declared
    /// charge — the physical ceiling, since a fragment cannot carry more
    /// charge than its precursor (multi-scan compare uses the largest
    /// declared charge in the batch). Override only with a physical reason,
    /// and state the cap whenever match counts are reported.
    #[arg(long, value_name = "N")]
    pub max_charge: Option<u8>,
    /// Fragment match tolerance in ppm (default 10).
    #[arg(long, default_value_t = 10.0, value_name = "PPM")]
    pub tolerance_ppm: f64,
    /// Minimum peak intensity to consider.
    #[arg(long, value_name = "F")]
    pub min_intensity: Option<f32>,
    /// Max matches to return.
    #[arg(long, value_name = "N")]
    pub max_matches: Option<usize>,
    /// JSON array of model-local site-scoped named neutral losses:
    /// [{"label":"H3PO4","deltaMass":97.9769,"position":88}]. The loss is
    /// emitted only for fragments whose residue span contains `position`
    /// (1-based index, nTerm, or cTerm). The sentinel value "engine" emits
    /// the fragmentation engine's own chemically-gated loss variants.
    #[arg(long, value_name = "JSON")]
    pub neutral_losses: Option<String>,
    /// Expected backbone sequence (compare mode): the carrier protein's
    /// residues as declared by the source PrSM/UniProt record, plain
    /// one-letter codes. Every candidate proforma must resolve to exactly
    /// this sequence — a mismatch aborts with the first divergent position.
    /// Omit this optional assertion when comparing alternative backbones.
    /// Each model retains its own sequence; positional comparison is then
    /// available only when the backbones are identical.
    #[arg(long, value_name = "SEQ")]
    pub sequence: Option<String>,
}

#[derive(Debug, clap::Args)]
pub struct EnvelopeRegionsArgs {
    /// Centroided raw mzML holding the MS2 scan.
    #[arg(value_name = "MZML")]
    pub file: String,
    /// Scan number (or nativeID).
    #[arg(long, value_name = "N")]
    pub scan: String,
    /// One complete site-resolved ProFormA model.
    #[arg(long, value_name = "PROFORMA")]
    pub proforma: String,
    /// Fragments to extract: `y3,y2`, ranges `b4-b7`, optional per-ion
    /// charge `y3:2`. A bare ion extracts every charge variant that has
    /// observed support inside its theoretical envelope.
    #[arg(long, value_name = "IONS")]
    pub ions: String,
    /// Restrict every extraction to one charge (default: each ion's own
    /// supported charges, capped by the precursor charge).
    #[arg(long, value_name = "N")]
    pub charge: Option<u8>,
    /// Isotope slots of margin on each side of the theoretical envelope.
    #[arg(long, default_value_t = 2, value_name = "N")]
    pub pad: u8,
}

/// Pure-theory fragment generation (the `get_theoretical_fragments` analogue).
#[derive(Debug, clap::Args)]
pub struct FragmentTheoreticalArgs {
    /// ProFormA string; TopPIC-style site tags and named/formula tags resolve
    /// to exact composition and delta mass. Numeric tags carry mass only.
    #[arg(
        long,
        value_name = "PROFORMA",
        required_unless_present = "models",
        conflicts_with = "models"
    )]
    pub proforma: Option<String>,
    /// Batch complete models [{"id":"a","proforma":"PEPTIDE"}, ...],
    /// inline or a JSON file. Batch output is JSON with one result per id.
    #[arg(long, value_name = "JSON")]
    pub models: Option<String>,
    /// JSON file of external mods: [{"name": "...", "formula": "..."}].
    #[arg(long, value_name = "JSON")]
    pub mods: Option<String>,
    /// Comma-separated fragment ion types: A,B,C,X,Y,Z (ETD: C,Z; HCD/CID:
    /// B,Y — the series are decided by the fragmentation method).
    #[arg(long, default_value = "C,Z", value_name = "IONS")]
    pub ion_types: String,
    /// Max charge state of fragment ions.
    #[arg(long, default_value_t = 1, value_name = "N")]
    pub max_charge: u8,
    /// JSON array of site-scoped named neutral losses:
    /// [{"label":"H3PO4","deltaMass":97.9769,"position":88}]. The loss is
    /// emitted only for fragments whose residue span contains `position`
    /// (1-based index, nTerm, or cTerm). The sentinel value "engine" emits
    /// the fragmentation engine's own chemically-gated loss variants.
    #[arg(long, value_name = "JSON")]
    pub neutral_losses: Option<String>,
    /// Output format: json (default, full fragment detail), mztab
    /// (standard mzTab-M molecular table), or mzpaf (HUPO-PSI mzPAF fragment
    /// annotations, one per line).
    #[arg(long, value_name = "FMT")]
    pub out_format: Option<String>,
    /// Omit each fragment's isotopeEnvelope from the json output (the BRAIN
    /// isotope stick the spectrum viewer draws). On by default for json —
    /// the envelope shares its kernel pass with apexMz, so it costs nothing
    /// extra; pass this flag to shrink bulk batch output.
    #[arg(long)]
    pub no_envelope: bool,
}

#[derive(Debug, clap::Args)]
pub struct UnimodGetArgs {
    /// UniMod SQLite path (default: $TDEASE_SYSTEM_ASSETS/reference/unimod.sqlite).
    #[arg(long, value_name = "DB")]
    pub db: Option<String>,
    /// User modification store (user-mods.sqlite, optional).
    #[arg(long, value_name = "DB")]
    pub user_db: Option<String>,
    /// Numeric UniMod record id.
    #[arg(long, value_name = "ID")]
    pub id: Option<u32>,
    /// Exact modification title (case-insensitive).
    #[arg(long, value_name = "TITLE")]
    pub title: Option<String>,
}

#[derive(Debug, clap::Args)]
pub struct UnimodSearchArgs {
    /// UniMod SQLite path (default: $TDEASE_SYSTEM_ASSETS/reference/unimod.sqlite).
    #[arg(long, value_name = "DB")]
    pub db: Option<String>,
    /// User modification store (user-mods.sqlite, optional).
    #[arg(long, value_name = "DB")]
    pub user_db: Option<String>,
    /// Fuzzy name text (8-tier explainable ranking).
    #[arg(long, value_name = "TEXT")]
    pub text: Option<String>,
    /// Observed mass to match by delta.
    #[arg(long, value_name = "MASS", allow_hyphen_values = true)]
    pub mass: Option<f64>,
    /// Tolerance value for --mass (absolute Da by default — we match PTM
    /// delta masses like +15.99, not whole proteins, so ppm is meaningless).
    #[arg(long, value_name = "TOL")]
    pub tolerance: Option<f64>,
    /// Tolerance unit: da (default) or ppm.
    #[arg(long, value_enum, default_value_t = ToleranceUnit::Da, value_name = "UNIT")]
    pub tolerance_unit: ToleranceUnit,
    /// Elements that must have a positive count (comma list).
    #[arg(long, value_name = "ELEMS")]
    pub contains: Option<String>,
    /// Elements that must have a zero count (comma list).
    #[arg(long, value_name = "ELEMS")]
    pub excludes: Option<String>,
    /// Restrict to residues a specificity must expose (comma list).
    #[arg(long, value_name = "RESIDUES")]
    pub residues: Option<String>,
    /// approved=only filters to approved records.
    #[arg(long, default_value = "any", value_name = "FILTER")]
    pub approved: String,
    /// Record-level classification gate (comma list). When set, a record is
    /// kept only if some specificity classification is in this set — the
    /// correct PTM filter (`approved` is unreliable: Oxidation is approved=0).
    /// Typical: "Post-translational,Artefact,Multiple,N-linked glycosylation,
    /// O-linked glycosylation,Other glycosylation". Check is record-level, so
    /// a mod with one biological site + chemical-derivative sites is kept.
    #[arg(long, value_name = "CLASSES")]
    pub keep_classification: Option<String>,
    /// Max matches to return.
    #[arg(long, default_value_t = 10, value_name = "N")]
    pub limit: usize,
}

#[derive(ClapArgs)]
pub struct RefineShiftBatchArgs {
    /// mzML file (read-only).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Inline JSON or a file of refine-shift requests. Each object requires `id` plus the
    /// camelCase equivalents of the single-command parameters.
    #[arg(long)]
    pub input: String,
}

#[derive(Debug, clap::Args)]
pub struct UnimodBatchSearchArgs {
    /// Inline JSON or a file of {"id": string, "mass": number, "tolerance": number}.
    #[arg(long, value_name = "JSON")]
    pub input: String,
    /// UniMod SQLite path (default: $TDEASE_SYSTEM_ASSETS/reference/unimod.sqlite).
    #[arg(long, value_name = "DB")]
    pub db: Option<String>,
    /// User modification store (user-mods.sqlite, optional).
    #[arg(long, value_name = "DB")]
    pub user_db: Option<String>,
    /// Tolerance unit for every input row: da (default) or ppm.
    #[arg(long, value_enum, default_value_t = ToleranceUnit::Da, value_name = "UNIT")]
    pub tolerance_unit: ToleranceUnit,
    /// Elements that must have a positive count (comma list).
    #[arg(long, value_name = "ELEMS")]
    pub contains: Option<String>,
    /// Elements that must have a zero count (comma list).
    #[arg(long, value_name = "ELEMS")]
    pub excludes: Option<String>,
    /// Restrict to residues a specificity must expose (comma list).
    #[arg(long, value_name = "RESIDUES")]
    pub residues: Option<String>,
    /// approved=only filters to approved records.
    #[arg(long, default_value = "any", value_name = "FILTER")]
    pub approved: String,
    /// Record-level classification gate (comma list).
    #[arg(long, value_name = "CLASSES")]
    pub keep_classification: Option<String>,
    /// Max matches returned per input row.
    #[arg(long, default_value_t = 10, value_name = "N")]
    pub limit: usize,
}

#[derive(Debug, clap::Args)]
pub struct UnimodBatchSolveArgs {
    /// JSON array of {"id": string, "mass": number, "tolerance": number,
    /// "fixedMods"?: [string]}. Row pins are added to global --fixed-mod pins;
    /// repeated entries add counts. Mass includes all pinned components.
    #[arg(long, value_name = "FILE")]
    pub input: String,
    /// UniMod SQLite path (default: $TDEASE_SYSTEM_ASSETS/reference/unimod.sqlite).
    #[arg(long, value_name = "DB")]
    pub db: Option<String>,
    /// User modification store (user-mods.sqlite, optional).
    #[arg(long, value_name = "DB")]
    pub user_db: Option<String>,
    /// Tolerance unit for every input row: da (default) or ppm.
    #[arg(long, value_enum, default_value_t = ToleranceUnit::Da, value_name = "UNIT")]
    pub tolerance_unit: ToleranceUnit,
    /// Maximum modifications in a combination.
    #[arg(long, default_value_t = 2, value_name = "N")]
    pub max_modifications: usize,
    /// Maximum distinct modifications.
    #[arg(long, value_name = "N")]
    pub max_distinct: Option<usize>,
    /// Disallow the same modification multiple times (default: allowed).
    #[arg(long = "no-allow-repeated", action = clap::ArgAction::SetFalse, default_value_t = true)]
    pub allow_repeated: bool,
    /// approved=only filters to approved records.
    #[arg(long, default_value = "any", value_name = "FILTER")]
    pub approved: String,
    /// Keep only records having at least one specificity in these classes.
    #[arg(long, value_name = "CLASSES")]
    pub keep_classification: Option<String>,
    /// Exclude records whose specificities are all in these classes.
    #[arg(long, value_name = "CLASSES")]
    pub exclude_classification: Option<String>,
    /// Element whitelist for record compositions.
    #[arg(long, default_value = "Br,C,Cl,H,N,O,P,S", value_name = "ELEMENTS")]
    pub element_whitelist: String,
    /// Sequence sites available to every combination.
    #[arg(long = "allowed-site", value_name = "RESIDUE[:POSITION]")]
    pub allowed_sites: Vec<String>,
    /// Fixed modifications pinned into every candidate of every row:
    /// TITLE|UNIMOD:n|USER:n|n, optionally :COUNT (default 1). Repeatable.
    #[arg(long = "fixed-mod", value_name = "TITLE|ID[:COUNT]")]
    pub fixed_mods: Vec<String>,
    /// Max candidates returned per input row.
    #[arg(long, default_value_t = 20, value_name = "N")]
    pub limit: usize,
}

#[derive(Debug, clap::Args)]
pub struct UnimodSolveArgs {
    /// UniMod SQLite path (default: $TDEASE_SYSTEM_ASSETS/reference/unimod.sqlite).
    #[arg(long, value_name = "DB")]
    pub db: Option<String>,
    /// User modification store (user-mods.sqlite, optional).
    #[arg(long, value_name = "DB")]
    pub user_db: Option<String>,
    /// Target neutral delta mass to match with a combination.
    #[arg(long, value_name = "MASS", allow_hyphen_values = true)]
    pub target_mass: f64,
    /// Tolerance value (absolute Da by default — the target is a small PTM
    /// delta, not a whole-protein mass, so ppm is meaningless).
    #[arg(long, value_name = "TOL")]
    pub tolerance: f64,
    /// Tolerance unit: da (default) or ppm.
    #[arg(long, value_enum, default_value_t = ToleranceUnit::Da, value_name = "UNIT")]
    pub tolerance_unit: ToleranceUnit,
    /// Max modifications in a combination (1..4).
    #[arg(long, default_value_t = 2, value_name = "N")]
    pub max_modifications: usize,
    /// Max distinct modifications.
    #[arg(long, value_name = "N")]
    pub max_distinct: Option<usize>,
    /// Disallow the same modification multiple times (default: allowed).
    #[arg(long = "no-allow-repeated", action = clap::ArgAction::SetFalse, default_value_t = true)]
    pub allow_repeated: bool,
    /// approved=only filters to approved records.
    #[arg(long, default_value = "any", value_name = "FILTER")]
    pub approved: String,
    /// Keep only records having at least one specificity in these classes.
    #[arg(long, value_name = "CLASSES")]
    pub keep_classification: Option<String>,
    /// Sequence sites available to the combination, repeated as RESIDUE[:POSITION].
    #[arg(long = "allowed-site", value_name = "RESIDUE[:POSITION]")]
    pub allowed_sites: Vec<String>,
    /// Fixed modifications pinned into every candidate: TITLE,
    /// UNIMOD:n / USER:n / bare id, optionally :COUNT (default 1).
    /// Repeatable. The pinned mass leaves the target before the search;
    /// every candidate reports those members with "fixed": true.
    #[arg(long = "fixed-mod", value_name = "TITLE|ID[:COUNT]")]
    pub fixed_mods: Vec<String>,
    /// Max candidates to return.
    #[arg(long, default_value_t = 20, value_name = "N")]
    pub limit: usize,

    /// Exclude records whose specificities are ALL in these classes
    /// (experimental labels, cross-linkers, chemical derivatives).
    #[arg(long, value_name = "CLASSES")]
    pub exclude_classification: Option<String>,
    /// Element whitelist for record compositions (comma list); records
    /// containing any element outside it are chemistry-impossible in a
    /// routine proteomics sample (Se/Fe/Ni/Al/Cu/Zn...) and are dropped
    /// BEFORE combination solving.
    #[arg(long, default_value = "Br,C,Cl,H,N,O,P,S", value_name = "ELEMENTS")]
    pub element_whitelist: String,
}

#[derive(Debug, clap::Args)]
pub struct UnimodResiduesArgs {
    /// UniMod SQLite path (default: $TDEASE_SYSTEM_ASSETS/reference/unimod.sqlite).
    #[arg(long, value_name = "DB")]
    pub db: Option<String>,
    /// User modification store (user-mods.sqlite, optional).
    #[arg(long, value_name = "DB")]
    pub user_db: Option<String>,
}

#[derive(Debug, clap::Args)]
pub struct UnimodAddArgs {
    /// User modification store to write (user-mods.sqlite).
    #[arg(long, value_name = "DB")]
    pub user_db: String,
    /// Modification title (the name ProFormA `[Title]` resolves).
    #[arg(long, value_name = "TITLE")]
    pub title: String,
    /// Exact net formula, e.g. "Br-H" or "C2H2O" (composition stays exact).
    #[arg(long, value_name = "FORMULA")]
    pub formula: String,
    /// Monoisotopic delta mass (Da).
    #[arg(long, value_name = "MASS")]
    pub mono_mass: f64,
    /// Average delta mass (Da); defaults to mono when omitted.
    #[arg(long, value_name = "MASS")]
    pub average_mass: Option<f64>,
    /// Comma-separated specificity sites (e.g. "S,T").
    #[arg(long, value_name = "SITES")]
    pub sites: Option<String>,
}

#[derive(ClapArgs)]
pub struct IsotopicBatchArgs {
    /// JSON file of variants: [{"formula":"C35H63N15O12","charge":1,"npeaks":4}, ...]
    #[arg(value_name = "FILE")]
    pub file: String,
    /// Truncate the envelope after this cumulative intensity (default 0.99).
    #[arg(long, default_value_t = 0.99)]
    pub truncate_after: f64,
    /// Drop peaks below this relative intensity (default 0.001).
    #[arg(long, default_value_t = 0.001)]
    pub ignore_below: f64,
}

#[derive(ClapArgs)]
pub struct ValidatePrecursorArgs {
    /// mzML file (read-only).
    #[arg(value_name = "FILE")]
    pub file: String,
    /// ProFormA string seeding the precursor theoretical envelope.
    #[arg(long)]
    pub proforma: String,
    /// JSON file of external mods (see parse-proforma).
    #[arg(long)]
    pub mods: Option<String>,
    /// Charge state of the precursor ion. When omitted, the MS2 scan's own
    /// precursor charge is read from the spectrum header (IsoForma uses the
    /// acquisition-assigned charge; override only when it is known to be wrong).
    #[arg(long)]
    pub charge: Option<i32>,
    /// The MS2 scan to validate (scan number or nativeID).
    #[arg(long, value_name = "SCAN")]
    pub scan: String,
    /// Aggregation tolerance per isotope (ppm).
    #[arg(long, default_value_t = 10.0)]
    pub tolerance_ppm: f64,
    /// Minimum matched MS1 isotope peaks for the precursor to be valid.
    #[arg(long, default_value_t = 3)]
    pub min_matched_peaks: usize,
}
