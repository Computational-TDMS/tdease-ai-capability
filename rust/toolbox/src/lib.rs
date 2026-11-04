mod constants;
mod envelope;
mod error;
mod fragment;
mod isotopic;
pub mod localize;
mod mzml;
pub use mzml::{open_reader, Reader};
mod mztab;
mod observation;
mod proforma;
mod query;
mod refine;
mod select;
mod seqgraph;
mod statistics;
mod theory;
// The UniMod query kernel (D12-R) is pure logic with no SQLite. TS reaches the
// UniMod query support is part of the native CLI and compiles unconditionally.
mod proforma_unimod;
mod unimod;

pub use constants::{neutral_from_mz, H2O, NEUTRON_SHIFT, PROTON};
pub use envelope::{
    aggregate_window_inner, aggregate_windows_inner, query_isotope_envelope_mzml,
    ElementalComposition, EnvelopeEvaluation, EnvelopeModel, IsotopeEnvelopeQuery,
    IsotopeEnvelopeResult, MatchedIsotope, RtAggregation,
};
pub use error::MassSpecError;
pub use fragment::{
    fragment_activation, fragment_precursor_charge, ion_types_for_activation, modification_deltas,
    query_candidate_comparison_deconvoluted, query_candidate_comparison_msalign,
    query_candidate_comparison_mzml, query_candidate_comparison_observation,
    query_candidate_comparison_reader, query_candidate_comparison_sum_reader,
    query_fragment_match_deconvoluted, query_fragment_match_msalign, query_fragment_match_mzml,
    query_fragment_match_observation, query_fragment_match_reader, theoretical_fragment_list,
    theoretical_fragment_list_proforma, theoretical_fragment_list_proforma_variants,
    theoretical_fragments_mzpaf, theoretical_fragments_mzpaf_proforma,
    theoretical_fragments_mzpaf_proforma_variants, CandidateComparisonQuery,
    CandidateComparisonResult, CandidateFragmentSummary, CandidateModel, DiscriminatingIon,
    EnvelopeMode, FragmentMatch, FragmentMatchQuery, FragmentMatchResult, FragmentMatchState,
    FragmentMatchStatus, FragmentSummary, IonSeriesSummary, IonType, IsotopeBorrowPair,
    MassModification, NamedNeutralLoss, SeriesDivergence, SeriesErrorStats, SiteWindowEvidence,
    TheoreticalFragment,
};
pub use fragment::{
    validate_fragment_parameters, PreparedCandidateComparison, PreparedFragmentMatch,
};
pub use isotopic::{
    generate_isotopic_pattern, generate_isotopic_pattern_from_formula,
    BasePeakToMonoisotopicOffsetEstimator, CachingIsotopicModel, ElementCount,
    FractionalComposition, IsotopicModel, IsotopicModelName, IsotopicModels,
    IsotopicPatternGenerator, IsotopicPatternParams, IsotopicPatternQuery, IsotopicPatternResult,
    IsotopicPatternSpec, TheoreticalIsotopicDistributionScalingMethod, TheoreticalPeakPoint,
};
pub use localize::{
    parse_msalign_spectrum, query_prsm_localize_msalign, query_prsm_localize_spectra,
    DeconvolutedPeak, LocalizeEvidence, LocalizeRow, MsAlignSpectrum, PositionStats,
    PrsmLocalizeQuery, PrsmLocalizeResult, ScanLocalizeSummary,
};
pub use mztab::MzTabWriter;
pub use observation::{extract_spectrum_observation_reader, SpectrumMetadata, SpectrumObservation};
pub use proforma::{
    chemical_from_counts, mod_delta_mass, normalize_for_fragmentation, parse_proforma,
    split_net_formula, ExternalModification, FragmentationProForma, ParsedModification,
    ParsedProForma,
};
pub use proforma_unimod::{ProFormaResolution, ResolvedModification, UniModProFormaResolver};
pub use query::{
    execute_query_mzml, query_index_mzml, query_mzml, query_precursor_mzml, query_product_mzml,
    query_scans_mzml, query_sum_mzml, query_window_mzml, query_xic_mzml,
    select_spectrum_ids_reader, validate_precursor_mzml, validate_precursors_mzml, BatchSpectrum,
    IndexQueryResult, MassSpecQuery, MassSpecQueryResult, MassTolerance, MergedPeak, Peak,
    PeakOrder, PrecursorQuery, PrecursorQueryResult, PrecursorQuerySummary, ProductMatch,
    ProductQuery, ProductQueryResult, ScanPrecursorMetadata, ScanQuery, ScanQueryResult,
    ScansQueryResult, SpectrumSelection, SpectrumSelector, SumQueryResult, ValidatePrecursorQuery,
    ValidatePrecursorResult, WindowQueryResult, XicPoint, XicQueryResult,
};
pub use refine::{
    calibrate_mass_shift, calibrate_mass_shift_from_observed, calibrate_mass_shifts, fit_envelope,
    refine_mass_shift, validate_envelope_fit, FitEnvelopeEvaluation, FitEnvelopeQuery,
    RefineMassShiftQuery, RefineMassShiftResult, RefinedShiftEvaluation,
    DEFAULT_ENVELOPE_WINDOW_MZ,
};
pub use select::{
    query_select_mzml, PrecursorMetadata, SelectPlan, SelectPredicate, SelectQueryResult,
    SelectedSpectrum,
};
pub use seqgraph::{
    fine_bin, graph_search_msalign, nominal_bin, proforma_frame, BaseIon, GraphComposition,
    GraphSearchCandidate, GraphSearchQuery, GraphSearchResult, LibrarySpec, ModificationParams,
    PositionalModification, SearchModification, SequenceGraph, RESCALING_HIGH_PRECISION,
    RESCALING_NOMINAL, SHIFT_NAME,
};
pub use theory::{
    anchored_envelope, intensity_by_neutron_index, theoretical_envelope, AnchoredEnvelope,
    TheoreticalEnvelope,
};
pub use unimod::{
    composition_matches, display_id, mass_error, search_by_mass, search_by_name,
    solve_combinations, ApprovedFilter, CombinationCandidate, CombinationModification, FixedMod,
    MassError, MassMatch, MassSpecificity, NameKind, NameMatch, SearchQuery, SolveMassResult,
    SolveQuery, SolveSite, Tol, TolUnit, UniModAminoAcid, UniModFilters, UniModNameEntry,
    UniModNeutralLoss, UniModRecord, UniModSourceMeta, UniModSpecificity, USER_OFFSET,
};
