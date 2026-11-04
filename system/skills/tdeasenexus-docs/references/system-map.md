# Where TDEaseNexus keeps its files

Every root, what owns it, and whether deleting it is a repair or a disaster.
Read this before touching anything outside a Workspace.

## Four layers, one boundary

| Layer | What it is | Where it lives |
|---|---|---|
| Desktop shell | Electron main + renderer: windows, IPC, Backend supervision, Agent runtime | install directory (`TDEaseNexus/`), no user data |
| Scientific System | read-only shipped assets: skills, prompts, templates, reference data, native tools | install directory `resources/app/assets/` (dev: `tdease-ai-capability/system/`) |
| Backend | Python/FastAPI service + embedded Snakemake engine, Node packages, Run storage | writable data root `backend/` |
| Workspace | user-selected project folder: Runs, sessions, wiki, evidence | anywhere the user chose |

Local Workspace files remain usable independently of Backend storage. The
Backend owns execution history and results that have not been pulled locally;
deleting its storage loses those records and bytes even if local files survive.

## The writable data root

One root holds both `system/` and `backend/`. Where it is:

| Install | Data root |
|---|---|
| Packaged, default | `TDEaseNexus-Data/` **beside the install directory** (e.g. `D:\Apps\TDEaseNexus-Data` next to `D:\Apps\TDEaseNexus`) |
| Packaged, install dir not writable | drive root `X:\TDEaseNexus-Data\`, else the per-user profile |
| Development | `backend/` is the repo's `data/backend`; the writable `system/` is `%APPDATA%\TDEase\system` |

Interactive reinstall wipes `TDEaseNexus-Data/` (and legacy `TDEase-Data/`)
beside the install directory. Silent updates replace the application while
preserving the data root. Launch does not relocate old data folders.

A fresh root is deliberately kept space-free: TopPIC 1.8.x hard-exits when its
own directory contains a space. An existing root is never relocated for that.

```text
<data root>/
  system/                  # Client-owned runtime state (writable)
    settings.json          # backend host/port, providers (API keys), runtime, mcp
    workspaces.json        # workspace registry: id, local_path, name
    models-cache.json      # provider model metadata cache
    logs/                  # backend.log, electron-main.log, crash-dumps/
    agent/                 # agent session host state
      tools/               # provisioned tdeasenexus-toolbox (the bash PATH entry)
    cache/uv/              # PEP 723 skill environments
    cache/results/         # result cache
    wiki/                  # Client-wide knowledge pages (user-owned)
    skills/ prompts/ templates/ reference/            # created at startup;
                           # shipped assets are read from the product assets root
  backend/                 # Backend-private storage (derived + execution state)
    config.json            # deployment mode, dataRoot ".", host/port, logging
    database/db.sqlite3    # rebuildable index (runs, run_files, run_graph, ...)
    nodes/                 # Node catalog: one directory per package + _template/
    runs/<run_id>/         # one directory per Run (see below)
    store/                 # content-addressed blob store (2-hex prefix dirs)
    cache/  trash/         # trash/ holds deleted Runs pending retention
    .venv/                 # dev only: uv environment (packaged: runtime/venv/)
```

`nodes/` is user-manageable: a package the user edits is never overwritten by an
app update, and a deleted package is re-seeded at next start from the shipped
`resources/backend/nodes/`. A directory starting with `_` is not a package
(`_template/` is the copyable starting point).

## A Backend Run directory

```text
<data root>/backend/runs/<run_id>/
  run.json                 # frozen definition: workspace_files (path/blob/size),
                           # snakefile, config, nodes, definition_revision_id,
                           # local_root_id + local_run_path + workspace_run_path
  events.jsonl             # append-only lifecycle journal — the Run's truth
  revisions/rev-<seq>.json # one result manifest per finished execution range
  results.manifest.json    # latest manifest (files + workflow_files + digests)
  results/                 # engine-declared outputs, as produced
  logs/tasks/, logs/attempts/
  .tdease/
    engine/                # Snakemake state (the `.snakemake` junction points here)
    engine-request.json    # exactly what the engine was asked to run
    job-events.jsonl       # per-job engine events
    local-inputs.json      # each materialized input: path, sha256, source, strategy
  .runtime/{locks,incomplete,task-records}/
  <materialized inputs>    # hardlinked from the Workspace, never copied
                           # (node packages never enter the Run: the engine
                           #  reaches them via config["nodes"] catalog paths)
  Snakefile, config.yaml, scripts/   # written control files
```

## A Workspace

```text
<workspace>/
  runs/<run>/              # flat: sources, workflow files, returned outputs
    run.json               # client-visible identity + output policy
    Snakefile, config.yaml, <sample>.mzML, database.fasta, ...
    results/               # optional grouping the workflow itself declares
  analysis/<slug>/         # ad-hoc analysis intermediates (scratch until cited)
  runs/<run>/analysis/<slug>/   # Run-scoped analysis intermediates
  runs/<run>/evidence/, evidence/   # cited evidence for a Run / for the Workspace
  exports/                 # conversation exports
  wiki/                    # workspace-scoped knowledge pages
  .tdease/
    sessions/<session_id>/ # session.json, events.jsonl, turns.jsonl,
                           # snapshots/, content/, model-session/, subagents/, jobs/
    skills/                # workspace skill overrides (scope `workspace:`)
    attachments/ tmp/ quarantine/
    workspace.json permissions.json settings.json user-mods.sqlite
```

Run directories are flat siblings under `runs/`; the Backend never creates a
sibling Client Run directory for a Backend run id.

## Durable vs derived

**Durable — never delete to fix a fault.** Deleting these loses work:

| Path | Why |
|---|---|
| `<workspace>/**` | user files, Runs, sessions, evidence, wiki |
| `<data root>/system/settings.json` | provider keys and app configuration |
| `<data root>/system/workspaces.json` | which folders the app may open |
| `<data root>/system/wiki/` | Client-wide knowledge pages |
| `<data root>/backend/runs/<id>/run.json` + `events.jsonl` | a Run's definition and history |
| `<data root>/backend/store/` | result blobs referenced by manifests |
| User-edited or user-installed `<data root>/backend/nodes/<pkg>/` | package customizations not recoverable from the shipped catalog |

**Derived — safe to delete; the app re-creates it.** Deleting is a legitimate repair:

| Path | Re-created by |
|---|---|
| `<data root>/backend/.venv/` (or `runtime/venv/`) | next launch: `uv sync --locked --no-dev` (needs network or a warm uv cache) |
| `<data root>/backend/database/db.sqlite3` | next Backend start: index reconcile from `runs/*/events.jsonl` |
| Unmodified shipped `<data root>/backend/nodes/<pkg>/` | next launch: re-seeded from the shipped catalog; verify it has no user changes before removal |
| `<data root>/system/agent/tools/` | next launch: copied from the shipped native tool |
| `<data root>/backend/cache/`, `system/cache/` | next use |
| `runs/<id>/.tdease/`, materialized inputs | next execution of that Run (Node
packages are never materialized into a Run) |
| `<workspace>/.tdease/tmp/` | next use |

Deleting a *whole* Backend run directory removes the Run's history too — it is
not a derived-state repair, it is a Run deletion. Use the app's trash action.
