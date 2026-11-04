---
name: toppicr-py
description: Quantitative proteoform matrices, QC reports, differential analysis and derived tables/figures from TopPIC results.
---

# TopPICR-Py

Use completed TopPIC PrSM tables, the matching FASTA and optional TopFD MS1
features to construct a quantitative matrix.

    uv run --script "$TDEASE_SYSTEM_ASSETS/skills/toppicr-py/scripts/toppic_run.py" \
      --ms2 <run>/sample_ms2_toppic_prsm.tsv --fasta <run>/database.fasta \
      --ms1 <run>/sample_ms1.feature --output-dir <analysis-dir>

Omit --ms1 when unavailable. --help covers multi-file and --project-dir inputs.
Long builds use tdease.job.start and tdease.job.wait. The runner refuses a
nonempty output directory; --overwrite explicitly replaces its previous
outputs. Inspect failed/partial results before deciding to replace them.

Stdout links the manifest, report and section statuses. Resolve saved table
paths from run_manifest.json; expression.csv, features.csv and samples.csv
are the core outputs. Read warnings and stage metrics in relation to recorded
parameters. Zero in the matrix means not detected/recovered, not biological
absence. Without MS1, quantification uses MS2 and has no MS1-based back-fill.

## Downstream comparison

For a compatible independent two-group design:

    uv run --script "$TDEASE_SYSTEM_ASSETS/skills/toppicr-py/scripts/toppic_plot.py" \
      compare --expression <expression.csv> --features <features.csv> \
      --group-a <regex> --group-b <regex> --out <comparison.json>

Tier the result before judging it — a fold change resting on two detections in
the weaker group is a different finding from one seen in every sample. Tiers
annotate evidence strength; they do not retire findings. The largest
between-group separations and recurring unnamed-shift families stay visible at
every tier (references/interpreting-results.md, Prioritization semantics):

    uv run --script "$TDEASE_SYSTEM_ASSETS/skills/toppicr-py/scripts/toppic_shortlist.py" \
      --comparison <comparison.json> --occupancy <features_with_occupancy.csv> \
      --groups <A,B> [--out-shortlist tiers.csv --out-genes gene_concordance.csv]

Check the selected samples against the experimental design. The default is a
Python port of the reference normalization and empirical-Bayes moderated
t-test with BH correction. The result's method and definitions blocks record
transforms, parameters, exclusions and field meanings. Cite the Python port
as such, rather than claiming it ran R limma.

--normalize none --stats welch is a sensitivity analysis on detected log2
values. For pairing, repeated measurements, covariates, more groups or other
estimands, choose an appropriate analysis directly. Matrix construction does
not validate downstream test assumptions.

## Reports and interpretation

The build emits qc_report.html. To redraw completed saved results:

    uv run --script "$TDEASE_SYSTEM_ASSETS/skills/toppicr-py/scripts/toppic_plot.py" \
      report --project <analysis-dir>

The report reuses recorded data, including calibration, without refitting.
Absent inputs yield skipped sections. Inspect failed/skipped statuses and
present the relevant views; additional figures are useful when they answer
the question, not to complete a fixed tour.

- references/interpreting-results.md: table semantics and stage changes.
- references/interpreting-figures.md: section-specific interpretation.

Feature modification labels inherit the source search's assumptions: the
terminal form and mass residual in a source row are that search's model, not
facts about the sample. Read them in the context of the search that produced
them before interpreting their chemistry. Additional figures can be rendered
from the saved tables.
