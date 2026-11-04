---
name: ptm-database
description: Query UniProt, dbPTM, iPTMnet, ELM motifs, cohort, UniRef and BLAST for PTM candidates and sequence context.
---

# PTM and sequence references

Return reference records and query coverage. Gather them as chemical priors
for judgment — modification families, motif classes and the chemistry of their
linkages — not as evidence that the sample carries them. A focused lookup can
finish with its records; connecting priors to explanations belongs to the
investigation, not to this skill. UniMod chemical searches use the kernel
tooling; this skill does not measure spectra or assign chemistry.

## Query interface

Run from the Workspace, with an existing analysis output directory:

```bash
uv run --script "$TDEASE_SYSTEM_ASSETS/skills/ptm-database/scripts/query.py" \
  <command> <query-arguments> --out <analysis-dir>/reference.json
```

Stdout carries the matched records themselves — site rows, motif classes with
their matches, cohort exemplars, top homolog hits — so a query's evidence
lands in the CLI return, never only in a file. Every command also archives the
full structured query report via `--out` (current directory if omitted; keep
it in the analysis directory instead), and `--full-stdout` prints that report.
Stdout previews may be bounded; compare counts and truncation with displayed
records and read the report when omitted records could change the explanation.
UniProt reports retain all returned feature annotations and their evidence,
not just modified-residue rows. Inspect relevant sequence/processing records
when they could change the model or its coordinate mapping.

| Command | Input | Scientific question | Access |
|---|---|---|---|
| `dbptm` | `--accession`, optional `--type` | Which sites/types have been annotated on this protein? | Local mirror |
| `iptmnet` | `--accession` | Which site, enzyme and provenance records exist? | Local mirror |
| `elm` | `--sequence`, optional `--class-prefix` | Which motif-class patterns match this sequence? | Local mirror |
| `cohort` | `--regex`, optional `--type` | Which annotated sites share a grounded sequence pattern? | Local mirror |
| `uniprot` | `--accession` | What sequence, processing and entry annotations are available? | Live REST |
| `uniref` | `--accession` | Which related sequences belong to the cluster? | Live REST |
| `blast` | `--sequence`, optional `--email` | Which homologs align to this sequence/window? Reviewed Swiss-Prot by default. | Live REST |

Each source answers its stated question within its coverage. An empty result
does not answer a different source's question. Choose inputs from the observed
sequence, accession or justified window and preserve their coordinates and
provenance. Candidate compatibility with a complete model is an investigation
step, not a source property.

Use `--help` for command-specific options. Type aliases resolve through the
mirror's vocabulary table. dbPTM sequence windows are 21 residues, with the
modified residue at index 10; iPTMnet counts represent deduplicated sites,
not publication counts. Motif matches do not include the ELM website's
structural or compartment filters.

## Combined lookup

`flow --accession <accession> [--cohort-regex <regex>]` combines UniProt,
UniRef, BLAST, direct and homolog mirror lookups, and motif queries.
`--skip-blast` skips BLAST only; UniProt and UniRef still need network access.
Use the atomic local commands for offline queries.

Use atomic calls for focused questions or `flow` for the initial prior pass.
Reuse the resulting reports across mass targets and subsequent measurements.
Interpret coverage from the returned sources and records. Atomic calls cover
only their named source; flow covers the lanes it actually completed. Keep a
relevant unqueried source visible as pending or explicitly deferred, rather
than treating one empty lookup as completion of the reference investigation.
For a shift under active investigation, deferring a lane needs a scope reason
—an expectation of an empty result is not one— and the sequence-motif lane
(`elm`) cannot be deferred before an "uncharacterised" verdict.

Source failures remain in `lane_errors`. `source_status` reports each source;
`lane_status` summarizes groups and can be `partial` when only some sources
completed. Distinguish `hit`, `empty`, `unavailable`, `skipped` and `pending`;
an unavailable source is not a negative result. `complete` means the command
settled, not that all sources were available or the explanation is complete.
With sequence-only input, use the applicable sequence-based
commands. A cohort query needs a pattern grounded in the supplied sequence or
reference context; an absent pattern is not a reason to invent one.

## Source and interpretation

The local dbPTM, iPTMnet and ELM mirror is built by `scripts/fetch.py` during
maintenance. Read `references/databases.md` for source scope and access details;
use `--db` only for an explicitly supplied mirror. PhosphoSitePlus data must be
user-supplied and are not redistributed. Other listed resources may be access
pointers rather than bundled data.

A missing mirror makes its lanes unavailable; continue with available reference
sources and measurements. Mirror provisioning is maintenance, not a prerequisite
for completing the investigation. flow runs independent local queries before
slow homology requests and updates the same report while running
(`complete: false`). Read that report to reuse completed results if a job is
interrupted; pending sources remain unfinished. Supply the actual target
sequence with an accession when available, so motif queries use that sequence
even if UniProt cannot be reached. Keep its coordinate mapping to the entry
explicit. Atomic queries remain available for focused follow-up or recovery.

Return the source identity/date when available, the relevant records and the
report path. Separate exact-protein annotations, homolog annotations and motif
matches. An annotation supports a candidate; it does not establish that the
sample carries it. Empty results mean no matching knowledge in the queried
source, not absence of the modification.

When a record could explain part of a mass gap, the gap is not closed here:
return the record, its applicability and the remaining question. Whether to
develop the model, seek further context or measure a discriminating prediction
is an investigation decision; lookup completion does not settle it.
