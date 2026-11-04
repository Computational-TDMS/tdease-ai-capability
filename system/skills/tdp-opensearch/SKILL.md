---
name: tdp-opensearch
description: Inspect existing TDP results and prioritize source-linked observations; does not execute a TopPIC search.
---

# TDP OpenSearch

Turn a broad results question into traceable observations for follow-up.
Rankings and reference matches prioritize investigation, not chemistry.

For a named protein, feature, PrSM, scan or family, inspect it directly.
For a run-wide question, build an inventory and select groups relevant to the
request. A batch report can finish there; an investigation continues through
tdp-copilot. Supplied model comparisons use spectral-evidence directly;
TopPIC execution uses toppic-suite and workflow.

## Selection

Use inspect --summary-out for a disposable pass or ingest build for a reusable
inventory. Read references/tool-recipes.md for invocation. Preserve source IDs,
object grain, mass frame, uncertainty, terminal processing and exclusions.

For TopPIC-derived observations, read references/toppic-output.md before
interpreting source fields or extracting chemical candidates. An unexpected
shift is relative to the row's assigned sequence/processing and PTMs; its
N-terminal form is a search assignment. The reference supplies the source-row
reader and field contract without loading the search-execution Skill.

Associate observations through supported source relations and compatible
precursor/RT/mass facts. Do not merge physical carriers because their shifts
are similar. Rank by the scientific question, observed support, recurrence,
ambiguity and the value of the next measurement. Use a manageable attention
set and state its cap and coverage; it does not limit underlying accounting.

Return selected groups with their source chain, relevant observations, priors,
selection reason and next unresolved question. Existing discovery artifacts
can be reused without rebuilding the inventory.

references/attention-and-discovery.md explains grouping and reported grades.
The result reader and qualification parameters own the implemented filters;
describe those used instead of treating a ranking as a scientific verdict.
