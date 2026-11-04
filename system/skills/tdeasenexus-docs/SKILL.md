---
name: tdeasenexus-docs
description: Diagnose TDEase Desktop, Backend and installation faults; route workflow and scientific failures to their owning Skills.
---

# TDEase system diagnosis

Use this Skill for installation, startup, connection and runtime faults.
Workflow authoring belongs to `workflow`, TopPIC toolchain failures to
`toppic-suite`, and result interpretation to the scientific Skills.

## Diagnose before repairing

1. Read the actual error and identify its owner: Desktop, Agent, Backend,
   Scientific System or Workspace.
2. For a Run, start with `tdease.runs.diagnose` and its status/log evidence.
   A diagnostic classification guides the investigation; inspect the failed
   stage rather than treating the label as proof of the root cause.
3. Read `references/repair-playbook.md` before changing installation or runtime
   state. Use `references/system-map.md` to distinguish durable records from
   rebuildable projections. Restart only when it addresses the diagnosed fault.
4. Verify the affected operation after repair. Report what changed, the observed
   result and anything still unresolved.

Never edit Backend records to make a Run appear successful. Local Workspace
files remain usable offline, but remote-only results and execution history
still depend on Backend storage. Do not delete that storage as a repair shortcut.

## Read-only health report

```bash
uv run --script "$TDEASE_SYSTEM_ASSETS/skills/tdeasenexus-docs/scripts/doctor.py" \
  --workspace "$TDEASE_WORKSPACE" --json
```

Add `--data-root <known-path>` when diagnosing the installation or Backend.
The doctor checks the roots it is given, reports findings and next actions,
and exits non-zero on failed checks. It does not repair them. Use it when the
fault is not already localized; preserve exact diagnostic codes in the report.

## Run lifecycle

Submission records a frozen definition and queues execution. The worker
materializes that definition, binds Nodes, invokes Snakemake and publishes a
result manifest. Journals record lifecycle facts; definitions record what was
submitted; result revisions identify published files. Read
`references/run-mechanism.md` for the execution and recovery details.

## Access boundary

Use registered tools and known paths. File tools can reach the bound Workspace
and authorized product-data locations; provider credentials remain forbidden.
`system-wiki/` is a virtual file-tool prefix, not a shell directory. Ask for a
missing installation path rather than guessing it. Workspace cwd and command
checks do not make arbitrary bash an operating-system filesystem sandbox.
