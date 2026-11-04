# Failure Playbook

## Diagnose the boundary

Collect, in order:

1. Run status and attempt/revision
2. failed DAG rule and sample/wildcards
3. actual executable command and Node digest
4. exit code and bounded log tail
5. expected inputs and outputs
6. surviving child processes
7. CPU affinity, priority, memory, free disk, and index-directory growth

## Failure classes

### Definition

Symptoms: Snakefile syntax, wildcard resolution, missing config key, missing declared
input. Edit the Run, validate, then create a revision if a submitted definition changed.

### Protocol or permission

Symptoms: Workspace containment, forbidden path, missing capability, upload mismatch.
Do not edit scientific files. Follow the public-protocol remedy.

### Node environment

Symptoms: executable missing, dependency DLL missing, unsupported platform/version.
Inspect the Backend Node. Do not search Workspace paths for a substitute binary.

### Data

Symptoms: malformed FASTA, duplicate accession, incompatible mzML/feature pairing,
empty or corrupt spectrum file. Validate the smallest affected input before rerunning.

TopFD exiting with `No spectrum to read in mzML file!` means the acquisition has
no MS1 scans (MS2-only method), not a corrupt file: probe the mzML's scan levels,
then rerun with `topfd -o` (missing-level-one) and `toppic -x` (no-topfd-feature).
The two flags are one data fact — set both. The `toppic-raw-pipeline` template
carries this wiring in its config.yaml comments.

### Resource

Symptoms: disk growth, out-of-memory, service starvation, excessive thread use. Cancel
the whole process tree, preserve the diagnostic facts, and change one resource or method
choice at a time.

### Concurrency

Symptoms: index/output "in use", multiple TopPIC writers, overlapping Runs. Confirm the
process list and scheduler lease. Serialize shared writers; do not reduce per-job threads
as a substitute for serialization.

### Orphaned execution

Symptoms: Backend restarted while run.json still says RUNNING. Project the terminal state
from events/process ownership. Resume only if the frozen definition is unchanged.

## Retry decision

```text
same Snakefile + config + inputs
  -> resume original Backend Run

scientific or executable definition, or inputs changed
  -> create a new revision of the same Workspace Run
```

Node digests record provenance; catalog redeployment alone is not a resume gate.
The workflow Skill owns revision and recovery rules. A diagnostic probe should
answer one explicit unknown rather than repeat an unchanged failure.
