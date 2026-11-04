# Repair playbook

One principle: **repair means making derived state match durable truth.**
Never delete truth or edit Backend state to fake a result. Preserve durable
records and use their owning application operations for authorized corrections;
a derived-state repair does not justify deleting user data.

## Why "reinstall" appears to work

A reinstall re-provisions environments, Node packages and native tools. An
interactive reinstall also wipes the sibling data root; silent updates preserve
it. Reinstalling is therefore not a data-preserving repair shortcut. Prefer the
diagnosed component repair below. If product assets are damaged and reinstall
is necessary, account for durable configuration, history and unpulled results
before recommending that reset.

## First run (bootstrap)

A fresh install creates its writable data root and fills it on first launch: the
`system/` tree, `backend/config.json`, a `settings.json` (a packaged install
seeds the provider key once — user edits always win afterwards), the Node
catalog seeded from the shipped packages, the native tool provisioned into
`system/agent/tools/`, and the Python environment synced from the network into
`backend/runtime/venv`. Since 1.0.2 the installer ships neither a wheelhouse nor
a bundled interpreter, so **the first launch needs network**; skill scripts
fetch their own PEP 723 dependencies on first use.

What the user must supply: a folder to open as a Workspace, and a provider key
if the install did not seed one. What proves it worked: the window reaches the
main UI, the Backend answers, and `doctor.py` reports no fail.

## The ladder — stop at the first step that fixes it

1. **Observe before repairing.** Reproduce once, then read: the splash message,
   `system/logs/backend.log`, `system/logs/electron-main.log`, `crash-dumps/`,
   the Run's `events.jsonl` tail, `runs.diagnose`.
2. **Restart the app.** At launch the Desktop probes the Backend environment,
   re-syncs the Python environment when the probe fails, re-seeds missing Node
   packages, re-provisions native tools, and the Backend reconciles the index
   from `events.jsonl` and recovers orphaned Runs. Most faults are fixed here.
3. **Delete one named derived path, then restart.** Only from the derived table
   in `system-map.md`, and only after step 1 named the path.
4. **Collect the log evidence** and report exact facts (paths, event seq, error
   code). A repair without a diagnosis is a guess.

## Symptoms

**App never reaches the main window; splash shows an error.**
Read the splash text and `system/logs/backend.log`. If it names `uv`, `trampoline`
or a Python path, the install was moved/copied: the app rebuilds
`<data root>/backend/runtime/venv` (dev: `data/backend/.venv`) on the next
launch and needs network for that. If it still fails, delete that environment
directory and relaunch — it is derived. If the message says the data root is not
writable, install to a writable directory; the app otherwise falls back to the
per-user profile.

**UI says the Backend is offline / reconnecting.**
The Desktop picks the Backend port at spawn (configured port `0` means
OS-assigned); never assume a fixed port. Check `backend.log` for the actual bind
and the last exception.

**A Run stays `QUEUED` for a long time.**
Execution is serial: another Run is executing. Confirm with the Run list before
treating it as a fault. If nothing is executing, the worker died — restart.

**A Run stays `RUNNING` with no new events.**
The engine process died without a terminal event. The Backend's identity check
(`pid:create_time`) prevents a recycled pid from being mistaken for a live
engine; at next start the orphan is recovered, and `runs.resume` continues from
the journal. Do not hand-write events.

**A Run reached `FAILED`.**
Use the class table in `run-mechanism.md`. The fix is in the Run definition or
its inputs — editing Backend storage, permissions or the engine's state
directory to make it pass is a fabrication, not a repair.

**Results are missing from the Workspace (cloud badge).**
They are in the Backend manifest. Pull the file through the app
(`runs.pullFile`); never copy blobs out of `backend/store/` by hand — the
manifest's digest is the identity.

**A Node package is missing or a binding digest cannot be resolved.**
Check `backend/nodes/<id>/rules.smk`. A **deleted shipped package** is re-seeded
at next launch. A **user-edited package** is authoritative and is never
overwritten — repair it by hand or restore the shipped version deliberately.
Directories starting with `_` are not packages.

**The app was moved to another drive, machine, or folder name.**
Python environments are not relocatable: uv's Windows console-script trampolines
embed the old absolute path. The launch probe catches this and rebuilds. There
is no data-folder migration: the installer force-reinstalls and wipes any
`TDEaseNexus-Data/` or legacy `TDEase-Data/` beside the install directory, and a
moved install simply provisions a fresh data root on first launch.

**A Workspace cannot be opened / "not registered".**
`system/workspaces.json` maps ids to paths. Re-opening the folder re-registers
it. If the folder genuinely moved, remove the stale entry through the app, never
by editing another Workspace's Runs.

**No model / provider errors.**
Provider keys live in `system/settings.json` and never leave the main process.
Ask the user to enter the key in Settings; never read, echo or store a key.

**Disk pressure.**
Reclaimable: `backend/cache/`, `system/cache/`, `backend/trash/` (deleted Runs
pending retention). Not reclaimable: `backend/store/` holds blobs referenced by
result manifests.

**`tdeasenexus-toolbox: command not found` inside bash.**
The managed bin directory `system/agent/tools/` is provisioned at launch from
the shipped tool. Restart; if the shipped tool is absent, the install assets are
damaged.

**A session will not load.**
Session truth is `.tdease/sessions/<id>/events.jsonl`; the app trims a truncated
last line on load. Never delete a session directory to clear a UI glitch — that
destroys the record.

## Reporting a fault you cannot repair

State, in order: what was observed (message + where), which layer it belongs to,
what the derived/durable split says, what you tried, and the exact next action
for the user. A precise report is worth more than a speculative repair.

## What the Agent can and cannot reach

The shell receives controlled environment variables, including
`TDEASE_SYSTEM_ASSETS` (installed assets) and `TDEASE_WORKSPACE` (the bound
Workspace). Workspace cwd and command checks are not an OS filesystem sandbox.
Use only known, authorized paths; provider credentials remain forbidden.
`system-wiki/` resolves through file tools, not as a shell path. Prefer the
app's repair operations when they own the affected state.

To run `scripts/doctor.py` against the Desktop side, ask the user where the app
is installed, then pass the data root explicitly. It is `TDEaseNexus-Data` beside
the install directory; if that directory is not writable the app falls back to
`%APPDATA%\TDEaseNexus`. In a development tree the Backend data root is
`<repo>/data/backend` and the writable System root is `%APPDATA%\TDEase\system`
(pass `--system-root` for that one).

```bash
uv run --script "$TDEASE_SYSTEM_ASSETS/skills/tdeasenexus-docs/scripts/doctor.py" \
  --workspace "$TDEASE_WORKSPACE" --data-root "<the user's data root>"
```
