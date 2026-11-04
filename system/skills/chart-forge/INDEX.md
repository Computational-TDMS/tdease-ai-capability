# ChartForge routes

These routes produce reproducible scientific figures per [SKILL.md](SKILL.md).

ChartForge is the shared presentation layer for source-backed scientific
figures. A route has four parts:

```text
source artifact -> domain adapter -> chart contract -> example package
```

The domain adapter owns parsing, filtering, grouping, statistics, and
provenance. The example owns visual encoding and interaction. Shared report
panels under `system/skills/_shared/panels/` are adapters for product reports;
they are not additional template packages.

## Choose by scientific object

| Object or question | Route | Contract |
|---|---|---|
| TopFD MS1 feature intervals | [`feature-map`](examples/feature-map/) | flat rows with `rt0`, `rt1`, `mz`, `logI` |
| MS1 isotope families | [`ms1-family-map`](examples/ms1-family-map/) | one row per family |
| Observed versus theoretical isotope envelope | [`envelope-compare`](examples/envelope-compare/) | envelope panels |
| Fragment-ion isotope-window mirror (multi-source cards) | [`envelope-mirror`](examples/envelope-mirror/) | observed up, theoretical comb down, one card per source |
| Accumulated MS2 scans | [`mgf-accumulated`](examples/mgf-accumulated/) | scans, peaks, annotations |
| Fragment matches on a sequence | [`fragment-annotation`](examples/fragment-annotation/) | peaks, residues, annotations |
| Candidate discrimination for one MS2 | [`candidate-discrimination`](examples/candidate-discrimination/) | spectrum and candidate maps |
| Fragment support across a model | [`proteoform-coverage`](examples/proteoform-coverage/) | sequence and model rows |
| OpenSearch mass-shift census | [`opensearch-census`](examples/opensearch-census/) | clusters, marks, limits |
| Proteoform definition levels | [`grading-levels`](examples/grading-levels/) | `data.levels[]` |
| One-accession proteoform map | [`proteoform-map`](examples/proteoform-map/) | sequence lanes and PTM boxes |
| PCA or sample separation | `system/skills/_shared/panels/pca.js` | point rows with group and size |
| Differential abundance | `system/skills/_shared/panels/volcano.js` | point rows with log2FC and BH FDR |
| Categorical values | [`categorical-bars`](examples/categorical-bars/) | `data.rows[]` |
| Attention inspection | [`attention-bars`](examples/attention-bars/) | `data.panels[].rows[]` |
| MassQL traces | [`massql-traces`](examples/massql-traces/) | `data.traces[]` |
| PTM structure context | [`structure-ptm`](examples/structure-ptm/) | structure plus site metadata |

The route is chosen from the object represented by the rows. A family is not
an envelope, a feature interval is not a proteoform centroid, and a statistical
point is not a raw spectrum peak.

## The `annotations` field

Every spectrum route accepts an optional `annotations` array — call-time
custom flags the agent attaches without touching the adapter's science.
Both kinds render through the shared
`forgeMs2` atoms (`components/ms2-annotation.js`), so every spectrum figure
gains them without template changes:

- `{kind: "callout", mz, label, color?, dash?}` — dashed flag at one m/z
- `{kind: "band", from, to, label, color?}` — tinted m/z region

Keep annotations sparse: they point attention, they do not carry data.

## Package contract

Each directory under `examples/` is a stable kebab-case route ID and contains
`chart.js`, `data.json`, and `build.json`. Optional `gen-data.py` files show an
adapter for that source format. The input to `build.py` is exactly:

```json
{"data": {}, "meta": {}}
```

The authoritative data contract for a route is its `examples/<route>/data.json`
(structure) and `build.json` (build inputs) — read those before constructing
data for a new figure; `gen-data.py` is a worked adapter example, not the
contract.

Use an existing route when its contract matches. New package maintenance uses
`references/engineering-contract.md`; it is not a required step during analysis.

Before handing built figures or reports to the user, verify self-containment
in one shot — no unreplaced tokens, no external network references, a data
payload present — with `scripts/check-figures.py <files-or-dirs>`; it prints
one bounded JSON verdict. Do not re-derive these checks by hand per figure.

## Shared language

- Tokens: [`templates/tokens.css`](templates/tokens.css), journal presets, and
  [`templates/micro.js`](templates/micro.js).
- Interaction rules: [`paradigms/interaction.md`](paradigms/interaction.md).
- Data contracts and provenance: [`paradigms/data-shaping.md`](paradigms/data-shaping.md).
- HTML and export pipeline: [`paradigms/html-pipeline.md`](paradigms/html-pipeline.md).
- Scientific safeguards: [`references/scientific-integrity.md`](references/scientific-integrity.md).
- Compatibility inventory: [`references/template-catalog.md`](references/template-catalog.md).

The runtime entry point is [`SKILL.md`](SKILL.md). This is the maintained route
table; the old catalog path points here for compatibility.
