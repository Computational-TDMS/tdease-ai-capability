---
name: toppic-suite
description: Configure and diagnose TopPIC tools or inspect PrSM results; use workflow for Run submission and control.
---

# TopPIC Suite

Select the required conversion, deconvolution, indexing or search stage from
the objective. The connected Backend Node owns executable paths, parameters
and output declarations; examples do not override its installed contract.
Use workflow for authoring, submission, revisions and recovery.

Parameter choices depend on acquisition, database, search space and intended
error control. Recover facts from metadata and protocols; use suitable Node
defaults for routine choices. Ask only for unresolved experimental facts or
choices that materially change the intended analysis.

## Read results

TopPIC tables can contain preambles and version-specific columns. Use the
source reader for bounded schema, summary or row selection:

    uv run --script "$TDEASE_SYSTEM_ASSETS/skills/toppic-suite/scripts/prsm_read.py" \
      rows --path <table.tsv> --scan <scan> --limit <n>

Use --help for filters and schema/summary commands. Preserve source version,
terminal processing, localization intervals, multiple hits and missing fields.
Table Spectrum IDs, raw scan IDs and MSAlign records are distinct identifiers.

E/Q-values concern search-match error; MIScore concerns localization conditional
on a mass. Neither establishes chemical identity. Hand the question to
spectral-evidence, tdp-copilot or toppicr-py as appropriate.

## References on demand

- references/suite-and-files.md: stages and file roles.
- references/parameter-decisions.md: acquisition and search choices.
- references/database-resources.md: FASTA and modification inputs.
- references/failure-playbook.md: TopPIC-specific faults.
- system:tdp-opensearch, resource references/toppic-output.md: result fields
  and table grammar, shared with downstream analysis.
- references/identification-quality.md: scores and definition levels.
- references/toppic-workflow-execution.md: Wilson & Zhou 2022 reproduction
  protocol; use only when that protocol is relevant, not as default settings.
