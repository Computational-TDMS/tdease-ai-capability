# Proteoform Annotation Contract

OpenSearch is Copilot's source-linked discovery layer. Its
feature merge, mass-shift inventory, UniMod candidate generation, and sequence
interval records are consumed by downstream analysis and visualization.

## One observation record, small overlays

`opensearch.json` is the canonical compact inventory for mass-shift, family views,
HTML and focused Copilot work. It has one observation record per assignment,
plus typed family/routing collections, and preserves the complete source chain
without copying raw tables:

```text
assignment -> feature -> proteoform -> PrSM -> scan
```

Each record retains the exact observed shift and mass basis, feature mass/RT,
source ProForma, TopPIC fixed/variable PTMs, terminal forms,
unexpected-modification intervals, precursor/proteoform masses, activation,
quality and fragment counts. The assigned-protein facts are the single owner of source
assigned/declared forms; a missing field is unknown, and downstream code
must not reconstruct it from a rounded shift.

The compact bundle keeps routing next to its observation, so no assignment row
is copied into a second queue. Its default UniMod pass stores only a bounded
top-three summary; `--include-unimod-candidates` retains full candidate arrays.
`--skip-unimod-candidates` is available for pure recurrence/family triage.
The embedded `modificationLibrary` array can be extracted as `resolve-shift`
library input; the whole inventory bundle is not a library JSON array. This
projection contains returned single-modification proposals, not every possible
combination or mechanism. Empty or bounded contents do not exhaust chemistry.

## Two different outputs

`observations` preserve source measurements and reported assignments:

- source `Feature ID`, `Proteoform ID`, and `PrSM ID`;
- raw mass values and separately recorded derived/calibrated values with their basis;
- protein accession, sequence interval, residue text, RT, and recurrence;
- PrSM/feature counts.

Preserve the distinction between a reported assignment and a verified sample
identity. Calibration and interpretation do not overwrite the original source
measurement. Use accession and explicitly parsed gene/name fields with their
provenance; arbitrary description tokens are not inferred gene or family names.

`annotations` are projections over those facts:

They retain the observed mass, proposed label, role/stage, basis and source
references. Automatic overlays have `notIdentity: true`; their display label
is a reference-derived proposal, not an assignment to the sample.

An annotation never changes an inventory bin or feature identity. Automatic UniMod
matches are first-pass proposals, not final chemistry. AI/deep review may
support, revise or refute a proposal with measurements while retaining its
origin. An explicit human verdict is recorded through Copilot's human-decision
route; an overlay role or a model score does not create that verdict.

## TopPICR-compatible projection

`toppicr-py` consumes the same TopPIC proteoform rows for its proteoform
map: bars are sequence intervals, fill is spectral count, and modification
markers use the TopPICR mass-centroid/UniMod matching algorithm. Unmatched
masses remain displayed as Da values.

OpenSearch consumes the same annotation shape for mass-shift inventory figures.
Its sticks and counts remain source-derived; labels are overlays linked to
`basis` and `sources`. This makes the inventory useful to the batch Copilot
without turning a prior into an identity claim.

`ingest build` embeds these overlays in `opensearch.json`. Full
candidate evidence is embedded only when the optional candidate pass is enabled.
