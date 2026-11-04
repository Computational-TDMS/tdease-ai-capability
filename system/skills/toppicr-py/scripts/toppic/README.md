# TopPICPy

TopPICPy is a Python downstream analysis pipeline for processing the output of
[TopPIC](https://www.toppic.org/) top-down proteomics. It is a Python rewrite of
[TopPICR](https://github.com/PNNL-Comp-Mass-Spec/TopPICR) that converts TopPIC PrSM
files, optional TopFD MS1 feature files, and a FASTA database into a proteoform
expression matrix of "features × samples".

The package works both for direct analysis and for embedding into a larger
data-processing platform. The Python API, the command-line interface, and the
optional Streamlit interface all drive the same `TopPICPipeline` service, so there
is only one analysis code path to maintain.

Exports include a `calibration` table containing the actual post-FDR,
pre-clustering MS2 observations: raw/aligned RT, raw/adjusted/recalibrated
mass, source dataset/scan, and the mass-reference membership used by the
estimator. `run_manifest.json` records its path and reference dataset.
Reports consume these saved observations without rerunning calibration;
older exports without them cannot supply calibration panels.

## Current status

Version `0.1.0` supports the complete core workflow:

```text
Parsing → FASTA mapping → filtering/FDR → protein inference
→ retention-time alignment → mass recalibration → clustering
→ match-between-runs → expression matrix
```

This project has been validated against output from TopPIC 1.5.4 and 1.8.1+.
TopPIC 1.7.0 is handled by the same version-detection layer, but has not yet
reached the same validation coverage as the reference data.

## Installation

Install the pinned official GitHub release:

```bash
pip install "toppic-py @ git+https://github.com/11111tao/toppic-py.git@v0.1.0"
```

Install the optional Streamlit and plotting interfaces:

```bash
pip install "toppic-py[web] @ git+https://github.com/11111tao/toppic-py.git@v0.1.0"
```

Set up a development environment:

```bash
git clone https://github.com/11111tao/toppic-py.git
cd toppic-py
python -m venv .venv
.venv/bin/python -m pip install -e ".[web,parquet,dev]"
```

## Python API in five minutes

Assume the following project directory layout:

```text
project/
├── fasta/
│   └── database.fasta.gz
├── toppic/
│   ├── sample_a_ms2_toppic_prsm.tsv
│   └── sample_b_ms2_toppic_prsm.tsv
└── topfd/
    ├── sample_a_ms1.feature
    └── sample_b_ms1.feature
```

Run:

```python
from toppic import PipelineConfig, PipelineInputs, TopPICPipeline

inputs = PipelineInputs.from_project_directory("project")
config = PipelineConfig(
    fdr_threshold=0.01,
    occurrence_threshold=2,
    n_jobs=4,
)

result = TopPICPipeline(config).run(
    inputs,
    output_dir="results",
)

print(result.expression.shape)
print(result.features.head())
print(result.warnings)
```

For uploaded files, the host platform should first write them to a local temporary
directory, then hand them to TopPICPy as a project directory or explicit file
paths. The lifecycle of HTTP, S3, or browser uploads is not managed by TopPICPy.

## CLI in five minutes

```bash
toppic run \
  --project-dir project \
  --output-dir results \
  --fdr-threshold 0.01 \
  --occurrence-threshold 2 \
  --n-jobs 4
```

File paths can also be given explicitly; repeat the corresponding flags for
multiple samples:

```bash
toppic run \
  --fasta project/fasta/database.fasta.gz \
  --ms2 project/toppic/sample_a_ms2_toppic_prsm.tsv \
  --ms2 project/toppic/sample_b_ms2_toppic_prsm.tsv \
  --ms1 project/topfd/sample_a_ms1.feature \
  --ms1 project/topfd/sample_b_ms1.feature \
  --output-dir results
```

## Output

The Python API returns a `PipelineResult`:

- `expression`: the feature × sample intensity matrix;
- `features`: feature annotations;
- `samples`: sample annotations;
- `stage_metrics`: row-count statistics per processing stage;
- `warnings`: warnings about recoverable problems;
- `manifest`: inputs, parameters, version, runtime, and matrix dimensions.

With `output_dir` set, TopPICPy writes:

```text
results/
├── expression.csv
├── features.csv
├── samples.csv
└── run_manifest.json
```

The default output format is CSV. Parquet is also available once the `parquet`
optional dependency is installed.

## Validation summary

We previously evaluated this Python rewrite against the 10-sample JPR TopPICR
reference workflow. The goal of the validation was comparable workflow behavior
and scientific conclusions, not line-by-line, value-by-value identity.

| Metric | Published R results | Python reference-run results |
|---|---:|---:|
| PrSM count | 23,136 | 22,636 |
| Proteoform count | 6,878 | 6,560 |
| Amino acid sequence count | 3,679 | 3,411 |
| PrSM FDR | 0.079% | 0.037% |
| Proteoform FDR | 0.129% | 0.056% |
| Sequence FDR | 0.315% | 0.136% |
| Gene FDR | 0.784% | 0.492% |
| Missing rate after match-between-runs | 22.8% | 25.8% |

Validation also covered retention-time alignment, mass recalibration, protein
inference, clustering, match-between-runs, expression-matrix construction, PTM
visualization, and reproduction of paper Figures 2–6.

## Optional web interface

After installing the `web` optional dependency and cloning the repository, run:

```bash
TOPPIC_PROJECTS_ROOT=/path/to/projects streamlit run app.py
```

The interface is a demo and interactive-analysis layer; it still calls the same
`TopPICPipeline` service used by the Python API and the CLI.

## Known differences

- FASTA mapping uses Python joins on recognizable accessions instead of TopPICR's
  parallel substring-search implementation, which produces the already-documented
  count differences.
- For ambiguous cases, the tie-breaking rule in gene selection differs from the R
  workflow, so the unique gene count can be lower than in the R workflow.
- Missing expression values in the output matrix are represented as zeros, not
  `NA` as in R.
- The MSnSet log2, global, and loess normalization steps used by the paper's
  downstream statistical analysis are not part of the core pipeline.
- When no MS1 files are provided, the pipeline explicitly falls back to aggregated
  MS2 feature intensities and records a warning.

## Testing

Run the default package and integration tests:

```bash
python -m pytest -m "not jpr"
```

The full regression tests on external data are optional:

```bash
TOPPIC_JPR_DATA_DIR=/path/to/JPR \
  python -m pytest -m jpr
```

CI does not execute the R implementation.

## Citation and license

This project is a Python rewrite of TopPICR `v0.0.7` (commit `0acc421`). When
using this workflow, please cite TopPIC/TopPICR and the related scientific
publications.

TopPICPy is distributed under the GPL-3.0-or-later license. See [LICENSE](LICENSE)
and [NOTICE](NOTICE) for details.
