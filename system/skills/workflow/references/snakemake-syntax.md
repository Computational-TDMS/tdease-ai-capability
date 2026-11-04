# Snakemake syntax (project subset)

This is **not** the Snakemake manual. It is the subset this project's Run Snakefiles
actually use, calibrated against the shipped `toppicr-py` template and real Node
`rules.smk`. For anything outside this subset, read the Node's `rules.smk` first
(`tdease.nodes.get`) — the Node is the authority.

## Imports and Node include

```python
import os
from pathlib import Path
```

Node rules are pulled in through the config-injected binding: the Backend
injects the package's absolute catalog path as `config["nodes"][<id>]`, and
the Snakefile includes from it (`references/node-binding.md`). A root
Snakefile has exactly one binding block per Node it composes.

## config

`config` is the parsed `config.yaml`, a plain dict. Required scientific inputs
should fail clearly when absent; optional settings may use documented defaults:

```python
SAMPLES = config["samples"]  # required, source-grounded sample list
TOPFD_ACTIVATION = config.get("topfd", {}).get("activation", "FILE")
```

Use direct indexing or an explicit validation error for required keys. Do not
invent sample names or scientific inputs to make a missing configuration parse.

## expand()

Turns a wildcard pattern into a concrete file list over samples. This is how `rule all`
declares the Run's final outputs:

```python
rule all:
    input:
        expand("{sample}_ms2_toppic_proteoform.tsv", sample=SAMPLES)
```

`{sample}` is a wildcard; `SAMPLES` supplies its values. The pattern string is also the
output path convention — Node outputs land at the run root as `{sample}_...`.

## rule

A rule names one step. The fields this project uses:

| field | purpose |
|-------|---------|
| `input` | upstream files (Node outputs, staged inputs); may be Python expressions |
| `output` | files this rule produces; drives the DAG and `rule all` |
| `params` | non-file values passed to the Node (often `config.get(...)`) |
| `threads` | CPU threads; declared, enforced by Backend |
| `resources` | memory/time ceilings (`mem_mb`, `time_min`); enforced by Backend |
| `shell` | the command; usually calls the Node binary with `params`/`input`/`output` |
| `benchmark` | optional; writes a TSV of runtime stats |
| `log` | optional; captures stderr to a file |

Most Run-level Snakefiles do **not** write `shell`/`threads`/`resources` themselves —
those live in the Node's `rules.smk`, which the root Snakefile only `include:`s. The
root Snakefile typically defines only `rule all` and any orchestration glue (sample
expansion, serialization guards).

## Wildcards

`{sample}` in `output` is a wildcard, inferred from requested targets. Keep wildcard
expansion sample-driven via `expand(..., sample=SAMPLES)` so a Run's scope is exactly
`config.yaml`'s `samples:` list.

## Serialization

Some Nodes cannot run concurrently (e.g. TopPIC jobs sharing one FASTA index). The
Node's `rules.smk` declares this (often a `resources` lock or an explicit dependency on
an index-marker output). Do not add ad-hoc locks in the root Snakefile; if you believe
jobs collide, read the Node — serialization is the Node's responsibility.

## Shell portability

The engine executes `shell:` rules through cmd.exe on Windows. POSIX-only
commands (`cp`, `mv`, `rm`, `ln`) do not exist there — and a dry-run cannot
catch this, because a dry-run never executes shell commands. The failure
surfaces only at execution time, after the Run's compute has already run.
Snakemake tutorials use `cp` freely; that idiom is wrong here.

For any file operation in a Run-level glue rule, use the established
`python -c` one-liner pattern (same as the topindex Node's `rules.smk`):

```python
import sys

rule export_inputs:
    input:
        src="database.fasta",
    output:
        dst="analysis/inp/database.fasta",
    params:
        python_bin=sys.executable,
    shell:
        '"{params.python_bin}" -c "import shutil,sys; shutil.copyfile(sys.argv[1],sys.argv[2])" '
        '"{input.src}" "{output.dst}"'
```

`&&` chaining is safe (cmd.exe has it); `mkdir` without `-p` is safe. When in
doubt, do the whole operation in one `python -c` call — `shutil`,
`pathlib.Path`, and `os` cover everything `cp`/`mv`/`rm`/`ln` would.

## What not to do

- Do not invent `shell:` blocks that call Node binaries directly with hand-written
  flags — bind through the Node's `rules.smk` so its calibrated flags win.
- Do not write a generic command/Shell rule. There is no arbitrary-command workflow
  block (AGENT-CONTRACT.md §6).
- Do not infer Run terminal state from log text. The Backend reports status; read that.
