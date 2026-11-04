---
name: workflow
description: Author, validate, submit and repair Run workflows; deployed Nodes own executable contracts.
---

# Scientific workflows

Author one contained Run for the requested objective. Backend Nodes supply
versioned rules, binaries and declared outputs. Use matching templates as
drafts, adapting to real inputs rather than inventing missing data.

## Author

Inspect the relevant Node contracts and source inputs. Separate external files
from derived intermediates: external inputs need to exist and match their
format; derived files need valid producing rules/dependencies. Resolve missing
samples or experimental choices with the user. Fetch needed public references
with provenance when the objective defines them; see toppic-suite's
references/database-resources.md for FASTA access.

Keep one root Snakefile/config and Node-declared layouts. Bind through
config["nodes"] and copy the Node's explicit executable path/parameters.
A bare command name may select unrelated host software. Snakemake shell rules
run in cmd.exe on Windows; use portable Python for file operations rather than
POSIX commands. Read references/node-binding.md for binding and
references/snakemake-syntax.md for unfamiliar syntax.

## Execute

- tdease.workflow.run validates, freezes and submits a definition atomically.
- tdease.workflow.validate checks without execution.
- tdease.workflow.inspect inspects structure/Node availability.
- tdease.runs.resume continues an unchanged definition and inputs.
- Changed executable/scientific inputs need a new revision.
- tdease.runs.overwrite is an explicit full recompute.

Node digests record provenance, not a resume gate. After submission the app
supervises execution; wait as needed to fulfill the request. Use published
result manifests for outputs. Distinguish validation, submission, execution
and verified results; process success alone is not scientific validity.

## Repair

Use tdease.runs.diagnose, the failed rule/command, inputs, exit code and logs
to locate the failing boundary. Worker stage lines identify failures before
engine output. Correct the cause in its owner instead of compensating with
unrelated parameters or Run edits.

RUN_PLAN_STALE means frozen text differs from disk: read current text and
submit it. A recurring error with no new corrective information calls for its
code/message and blocker, not unchanged retries.
