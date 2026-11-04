---
name: tdp-execution
description: Select TDP processing stages and input dependencies; use workflow to execute and toolchain Skills for parameters.
---

# TDP pipeline selection

Choose the stages needed by the objective and existing data. Use workflow for
Run authoring/control and the deployed Node contract for executable parameters.

| Objective / available input | Route |
|---|---|
| Raw acquisition to open-search identification | toppic-raw-pipeline: msconvert → TopFD → TopPIC |
| Targeted known-modification search and/or MS1 features | informed-proteomics-pipeline: PbfGen, ProMex, MSPathFinderT as needed |
| Native MS1 deconvolution or isotope prediction | unidec |
| Existing result investigation | tdp-opensearch / tdp-copilot |
| Existing result matrices | toppicr-py |

Templates under templates/ are starting points, not obligatory pipelines.
They mirror authoritative system/templates/; adapt stages to the request.
The Informed-Proteomics template supports a discovery-to-cross-check route,
but another search is not required to interpret suitable existing evidence.

## Decision-changing constraints

- Inspect acquisition metadata when it determines flags. An MS2-only input
  needs paired TopFD -o / TopPIC -x behavior; see toppic-suite.
- A cross-check can use a focused FASTA from identified accessions and a bounded
  modification model. State this restricted search space; agreement is not
  independent confirmation of the original hypothesis.
- Modification definitions need grounded masses/compositions. Use supplied
  artifacts or the native reference library; unknown shifts can remain
  hypotheses for spectral investigation.
- Preserve feature artifacts required downstream. ProMex produces .ms1ft;
  a method accepting TopFD features need not run ProMex for a redundant map.
- MSPathFinderT and PbfGen can reuse existing outputs. Check inputs and resolved
  overwrite behavior so a new request does not silently consume stale results.

Installed Nodes own binaries, parameters and outputs. Use toppic-suite for
TopFD/TopPIC specifics and references/node-manuals.md for Informed-Proteomics.
Verify outputs before interpretation; process exit alone does not establish
a usable scientific result.
