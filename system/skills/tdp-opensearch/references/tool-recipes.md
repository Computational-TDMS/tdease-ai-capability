# OpenSearch invocation

OpenSearch is a script, not a PATH command or registered tool:

    uv run --script "$TDEASE_SYSTEM_ASSETS/skills/tdp-opensearch/main.py" <command> ...

Commands include inspect, ingest, families and figures. Use command-specific --help for
unfamiliar arguments. Select source-linked records with inspect or build a
reusable projection with ingest build.

## Optional mass-reference work

| Profile | Effect |
|---|---|
| --skip-unimod-candidates | Inventory without UniMod candidate generation |
| default | Bounded search with compact candidate summaries |
| --include-unimod-candidates | Retain full candidate arrays |

These flags apply to ingest build; inspect has its own read-side schema and
does not run this candidate pass. For initial unknown-shift discovery, use
--skip-unimod-candidates and hand the
selected context to Copilot's initial reference-prior pass. Subsequent focused
mass searches can use those priors; rebuilding the whole inventory is
unnecessary. Existing bundles remain reusable. These profiles control cost
and payload, not confidence.

## Family and refinement behavior

families finds MS1 mass ladders under spacing, ppm and shared-RT constraints.
Supply --spacing-min and --spacing-max explicitly in neutral Da; choose their
bounds from the question or observed range. MS1-only ladders do not establish
protein identity or sequence interval. The identified-family projections in
ingest also use source protein/interval relations. Preserve these distinct
grains; ladder gaps and recurrence are observations, not assigned chemistry.

With a same-stem mzML beside --feature, or explicit --mzml, the tool performs
RT-stacked refine-shift by default; --no-refine skips it. Failed shape QC keeps
raw mass. Read the recorded basis, warnings and calibration evidence; successful
averagine refinement is not an external accuracy standard. Refinement can
take long enough to need a job.

ingest build instead controls its embedded calibration with --skip-calibration.
Use the option belonging to the invoked command; refinement and inventory
remain distinct operations.

The bundle retains modificationLibrary for focused follow-up. It does not run
resolve-shift over every residual.

## Figures

The report renders supported views whose inputs exist (census, families,
grading). --marks adds curated overlays. Ranking fields are triage data, not
a separate attention report. --census / --families support existing single-view
artifacts. Add the raw feature table/accession for a proteoform view only when
that view answers the request. Read returned paths and section status.
