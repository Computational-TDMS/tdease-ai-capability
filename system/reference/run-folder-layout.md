# Workspace Run flat layout

A TDEase Workspace Run is flat: data, active workflow files and returned scientific
outputs sit directly at `workspace/runs/<run>/`. No `raw/`, `fasta/`,
`results/` or per-sample subdirectories are required.

## Why flat

- TopFD and related tools naturally write outputs next to their inputs.
- Relative paths stay identical through upload, Desktop link optimization, Backend
  execution and selective output synchronization.
- Every Run directory remains understandable with ordinary filesystem tools.
- The `workspace/runs/` parent is only a stable discovery boundary.

## Canonical Client layout

```text
workspace/runs/<run>/
  run.json
  Snakefile
  config.yaml
  <sample>.mzML
  database.fasta
  <sample>_ms1.feature
  <sample>_ms1.msalign
  <sample>_ms2.msalign
  <sample>_ms2_toppic_proteoform.tsv
  analysis.svg
  logs/
```

Backend-private engine persistence, Node payloads, caches and temporary files do not
belong in this directory and are never Client identities.

## Conventions

- Sample names and workflow input/output roles come from `config.yaml` and the
  Snakefile, not from Backend paths.
- A workflow may use conventional names such as `database.fasta` and
  `<sample>.mzML`, but TDEase does not impose extra directory grouping.
- Returned outputs are selected from workflow declarations, sample configuration and
  user choice. Backend does not mirror every technical artifact.
- Logs and events are synchronized into the Run for inspection.

## Submission

1. Create or select `workspace/runs/<run>/`.
2. Instantiate a System Snakefile template into it, or author one directly.
3. Configure samples and parameters and put required source files in the Run root.
4. Dry-run and submit the exact relative-path/digest manifest.
5. Backend requests missing bytes and materializes a private execution directory.
6. Desktop may satisfy transfer with verified direct reads, links or copies; Server
   uploads bytes. Both are the same protocol operation.
7. Synchronize execution logs and selected declared outputs back into the Workspace
   Run.
