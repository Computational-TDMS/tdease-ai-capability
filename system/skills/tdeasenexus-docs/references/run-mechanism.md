# How a Run actually executes

The mechanism from "user presses Run" to "results in the Workspace", and the
vocabulary for saying where it stopped. Use this to locate a fault, not to
re-derive science.

## The chain

1. **Plan (Client).** The Client scans the Workspace Run, hashes files, builds
   an opaque `planRef`. The model never authors or sees a raw commit plan.
2. **Dry-run (Backend).** The real Snakemake dry-run decides whether the
   workflow is valid. A dry-run is not execution.
3. **Submit (Backend).** `POST /api/v2/runs` validates identity, containment,
   blobs and the node bindings, appends `run.submitted`, then `run.queue_requested`.
   Submission persists the definition and queue request; it does not yet
   materialize the workflow's working files.
4. **Queue.** One Backend deployment has one global execution slot: Runs start
   in submission order, at most one executing. Snakemake alone decides
   parallelism inside that Run.
5. **Materialize (worker, immediately before execution).** The worker exposes
   the revision it is about to start: inputs are **hardlinked** from sealed,
   read-only Backend CAS blobs, control files (`Snakefile`,
   `*.smk`, `config.yaml`, `scripts/`) are written, and engine state is created
   under `.tdease/engine` with a `.snakemake` junction pointing at it. Node
   packages are never materialized: the engine request's config carries each
   bound package's absolute catalog path (`config["nodes"]`). A link the
   platform refuses **fails the Run** with a structured error — the Backend
   never silently copies multi-GB inputs.
6. **Execute.** A pinned Snakemake child process runs the Snakefile from
   `.tdease/engine-request.json`. The engine is supervised; the worker writes a
   process identity (`pid:create_time`) so a recycled pid can never be mistaken
   for a live engine.
7. **Journal.** Every lifecycle fact is appended to `events.jsonl`
   (`schema_version: 3`, monotonic `seq`, `type`, `data`). The status projection
   is a pure fold of that journal; SQLite is only a rebuildable index.
8. **Results.** When an execution range ends, the worker writes
   `revisions/rev-<seq>.json` and refreshes `results.manifest.json`: declared
   outputs first, then other produced files, each with size and SHA-256. Blobs
   land in the Backend store. The Client pulls a manifest entry into the
   Workspace Run on request (`runs.pullFile`).

## Lifecycle vocabulary

Statuses: `QUEUED`, `RUNNING`, `SUCCEEDED`, `FAILED`, `CANCELLED` — declared in
`src/contracts/run_lifecycle.json` (Backend) and mirrored by the Client. Events:
`run.submitted`, `run.queue_requested`, `run.started`, `run.succeeded`,
`run.failed`, `run.cancelled`, plus the non-lifecycle `run.action`,
`workflow.graph` and `task.started`/`progress`/`succeeded`/`failed`. Task output
is a log file, not an event.

Two rules that surprise people:

- **No status locks a definition.** A Run is always editable and always
  continuable. Submitting again appends a new definition revision and a queue
  request from any status. If an execution is in flight it keeps its own journal
  range and manifest, and the new definition is queued behind it: the projected
  status stays `RUNNING` until that range terminates, then becomes `QUEUED`.
  Reproducibility is **recorded, not enforced** — the revision id is a content
  hash, every execution range names the revision it ran, and the manifest names
  it too.
- **Execution is deliberately serial.** "Runs should run in parallel" is a
  regression, not an improvement.

## Where it stopped — the failure classes

The Client's `runs.diagnose` classifies from status, error text and log tail:

| Class | Meaning | First action |
|---|---|---|
| `IN_PROGRESS` | QUEUED or RUNNING | `runs.wait`, `runs.tailLog` |
| `BACKEND_ORPHANED` | the worker/engine died with the app or machine | `runs.resume` |
| `PERMISSION_PROTOCOL` | path escape, forbidden path, containment refusal | fix the path or the Workspace root; never patch Backend permissions |
| `RESOURCE_EXHAUSTION` | out of memory, no space, budget refused | reduce the request or free space |
| `CONCURRENCY` | file in use, locked, another process | stop the other process; do not "unlock" by deleting state |
| `INPUT_DATA` | missing or corrupt input, bad FASTA/mzML | fix the input in the Workspace |
| `NODE_RUNTIME` | node binary/DLL missing or wrong platform | `nodes.get`, re-seed the package |
| `WORKFLOW_DEFINITION` | Snakefile/config/syntax/wildcard error | edit the Run definition |
| `EXECUTION_UNKNOWN` | FAILED with no classified signature | read the log tail and `runs.diagnose` facts |

Backend errors arrive in one envelope with a derived `domain`
(`workspace`, `run`, `engine`, `agent`, `backend`, `environment`) and a closed
`code`. A workflow that fails is **not** an API error: the Run reaches `FAILED`
with an `engine` code and HTTP 200. HTTP 4xx/5xx means the *request* was wrong.

## Reading a Run's evidence without guessing

- `runs.status` → projected status, last event seq, error.
- `runs.tailLog` → bounded engine log tail.
- `runs.listFiles` → the results manifest (Backend truth for what was produced).
- `runs.diagnose` → graph + failed/running jobs + log tail + classified next action.
- On disk: `events.jsonl` is authoritative; `results.manifest.json` and
  `revisions/*.json` record what each range produced; `.tdease/engine-request.json`
  records exactly what the engine was asked to do.

Never reconstruct a DAG or a status by parsing a Snakefile or a log; ask the
Backend, or fold the journal.
