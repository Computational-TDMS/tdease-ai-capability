# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""TDEaseNexus doctor — a read-only health report for the installation, the
Backend data root and a Workspace.

It never writes, never deletes and never talks to the network. Every finding is
a verdict plus the next action, so a fault is reported instead of guessed.

    uv run --script doctor.py --workspace "$TDEASE_WORKSPACE"
    uv run --script doctor.py --workspace C:\\proj --data-root D:\\Apps\\TDEaseNexus-Data
    uv run --script doctor.py --workspace C:\\proj --data-root D:\\Apps\\TDEaseNexus-Data --json

Roots are never guessed from product environment variables: the Workspace and
the read-only System assets default to the two sandbox pointers
(``TDEASE_WORKSPACE`` / ``TDEASE_SYSTEM_ASSETS``); the data root must be passed
explicitly (a small set of conventional locations is probed and reported when
it is not). Exit code is 1 when any check fails.
"""

from __future__ import annotations

import argparse
import json
import os
import sqlite3
import sys
from dataclasses import dataclass, field
from pathlib import Path

OK, WARN, FAIL, SKIP = "ok", "warn", "fail", "skip"
RANK = {OK: 0, SKIP: 1, WARN: 2, FAIL: 3}
MARK = {OK: "[ok]  ", WARN: "[warn]", FAIL: "[fail]", SKIP: "[skip]"}

# A journal is read line by line; a corrupt tail must not hide the rest.
_MAX_RUNS_DETAIL = 50
_LOG_TAIL_BYTES = 64 * 1024
_ERROR_TOKENS = ("error", "critical", "traceback", "exception")


@dataclass
class Check:
    id: str
    title: str
    status: str
    detail: str
    fix: str | None = None


@dataclass
class Report:
    roots: dict[str, str] = field(default_factory=dict)
    checks: list[Check] = field(default_factory=list)

    def add(self, id: str, title: str, status: str, detail: str, fix: str | None = None) -> None:
        self.checks.append(Check(id, title, status, detail, fix))

    def worse(self, current: str, candidate: str) -> str:
        return candidate if RANK[candidate] > RANK[current] else current

    @property
    def exit_code(self) -> int:
        return 1 if any(check.status == FAIL for check in self.checks) else 0


def _read_json(path: Path) -> tuple[object | None, str | None]:
    """Parse JSON, returning (value, error). Never raises."""
    try:
        return json.loads(path.read_text(encoding="utf-8-sig")), None
    except FileNotFoundError:
        return None, "missing"
    except (OSError, UnicodeDecodeError) as exc:
        return None, f"unreadable: {exc}"
    except json.JSONDecodeError as exc:
        return None, f"invalid JSON: {exc}"


def _pid_alive(pid: int) -> bool:
    """Whether a live process with this pid exists. Identity is not verified."""
    if pid <= 0:
        return False
    if os.name == "nt":
        import ctypes

        PROCESS_QUERY_LIMITED_INFORMATION, STILL_ACTIVE = 0x1000, 259
        kernel32 = ctypes.windll.kernel32
        handle = kernel32.OpenProcess(PROCESS_QUERY_LIMITED_INFORMATION, False, pid)
        if not handle:
            return False
        try:
            code = ctypes.c_ulong()
            if kernel32.GetExitCodeProcess(handle, ctypes.byref(code)):
                return code.value == STILL_ACTIVE
            return False
        finally:
            kernel32.CloseHandle(handle)
    try:
        os.kill(pid, 0)
        return True
    except (OSError, ProcessLookupError):
        return False


def _tail_lines(path: Path, limit: int) -> list[str]:
    try:
        with path.open("rb") as stream:
            stream.seek(0, os.SEEK_END)
            size = stream.tell()
            stream.seek(max(0, size - _LOG_TAIL_BYTES))
            text = stream.read().decode("utf-8", errors="replace")
    except OSError:
        return []
    return [line.rstrip() for line in text.splitlines()[-limit:]]


def _read_journal(path: Path) -> tuple[list[dict], int]:
    """Return (parsed events, count of unparseable lines)."""
    events: list[dict] = []
    bad = 0
    try:
        with path.open("r", encoding="utf-8", errors="replace") as stream:
            for line in stream:
                line = line.strip()
                if not line:
                    continue
                try:
                    value = json.loads(line)
                except json.JSONDecodeError:
                    bad += 1
                    continue
                if isinstance(value, dict):
                    events.append(value)
                else:
                    bad += 1
    except OSError:
        return [], 0
    return events, bad


def _size(path: Path) -> str:
    try:
        value = float(path.stat().st_size)
    except OSError:
        return "unknown size"
    for unit in ("B", "KB", "MB", "GB"):
        if value < 1024 or unit == "GB":
            return f"{value:.0f} {unit}" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024
    return "?"


def _native_tool(directory: Path) -> Path | None:
    for name in ("tdeasenexus-toolbox.exe", "tdeasenexus-toolbox"):
        candidate = directory / name
        if candidate.is_file():
            return candidate
    return None


def _shipped_native_tool(assets: Path) -> Path | None:
    """The shipped native CLI: `tools/` in a development System root, the
    `native/` resource beside the packaged assets otherwise."""
    for directory in (assets / "tools", assets.parent / "native", assets.parent.parent / "native"):
        tool = _native_tool(directory)
        if tool:
            return tool
    return None


def _list_dirs(directory: Path) -> list[Path]:
    try:
        return sorted(p for p in directory.iterdir() if p.is_dir())
    except OSError:
        return []


# ── checks ────────────────────────────────────────────────────────────────────


def check_workspace(report: Report, workspace: Path | None) -> None:
    if workspace is None:
        report.add("workspace", "Workspace", SKIP, "none given (pass --workspace or set TDEASE_WORKSPACE)")
        return
    if not workspace.is_dir():
        report.add("workspace", "Workspace", FAIL, f"not a directory: {workspace}",
                   "Open an existing folder as a Workspace; do not create it by hand.")
        return
    has_layout = (workspace / ".tdease").is_dir()
    runs = _list_dirs(workspace / "runs")
    sessions = _list_dirs(workspace / ".tdease" / "sessions")
    # A Run the app created carries run.json; a hand-made data folder under
    # runs/ is plain files and not a fault. Only an app-shaped folder missing
    # its definition is a finding.
    app_created = [p for p in runs if (p / "run.json").is_file() or (p / ".tdease").is_dir()]
    broken = [p.name for p in app_created if not (p / "run.json").is_file()]
    report.add("workspace.layout", "Workspace layout",
               WARN if broken else OK,
               f"{len(runs)} Run folder(s) ({len(app_created)} app-created), "
               f"{len(sessions)} Session(s), .tdease {'present' if has_layout else 'absent'}"
               + (f"; missing run.json: {', '.join(broken[:5])}" if broken else ""),
               "An app-created Run folder without run.json lost its identity; "
               "recover it from the Backend journal or re-create the Run." if broken else None)


def check_system_assets(report: Report, assets: Path | None) -> None:
    if assets is None:
        report.add("system-assets", "System assets", SKIP, "not located (pass --system-assets)")
        return
    if not assets.is_dir():
        report.add("system-assets", "System assets", FAIL, f"not a directory: {assets}",
                   "The install assets are missing or the path is wrong. A reinstall is "
                   "justified only when resources/app/assets is genuinely damaged.")
        return
    categories = [name for name in ("skills", "prompts", "templates", "reference") if (assets / name).is_dir()]
    skills = sorted(p.parent.name for p in (assets / "skills").glob("*/SKILL.md")) if (assets / "skills").is_dir() else []
    report.add("system-assets", "System assets", OK if len(categories) == 4 else WARN,
               f"{', '.join(categories) or 'no categories'}; {len(skills)} skill(s)")
    tool = _shipped_native_tool(assets)
    report.add("system-assets.tool", "Shipped native tool",
               OK if tool else FAIL,
               f"{tool} ({_size(tool)})" if tool else
               f"missing under {assets / 'tools'} and the sibling native/ resource",
               None if tool else "The shipped tdeasenexus-toolbox is absent — the install is incomplete.")


def check_system_data(report: Report, system: Path | None) -> None:
    if system is None or not system.is_dir():
        report.add("system-data", "System data", WARN,
                   f"no System root at {system}" if system else "not located",
                   "Pass --system-root, or --data-root for an install whose data root holds "
                   "system/ (development trees keep it in %APPDATA%\\TDEase\\system).")
        return

    settings, error = _read_json(system / "settings.json")
    if error:
        report.add("system-data.settings", "settings.json", FAIL, f"{system / 'settings.json'}: {error}",
                   "Delete settings.json and relaunch: the app writes a default file. "
                   "Provider keys must be re-entered afterwards.")
    else:
        providers = (settings or {}).get("providers") or {}
        configured = [name for name, value in providers.items()
                      if isinstance(value, dict) and str(value.get("api_key") or "").strip()]
        report.add("system-data.settings", "settings.json", OK,
                   f"{len(configured)} provider(s) with a key: {', '.join(configured) or 'none'}")

    registry, error = _read_json(system / "workspaces.json")
    if error:
        report.add("system-data.registry", "workspaces.json", FAIL, f"unreadable: {error}",
                   "Re-open each Workspace folder in the app; do not hand-edit Runs.")
    else:
        entries = (registry or {}).get("workspaces") if isinstance(registry, dict) else None
        entries = entries if isinstance(entries, list) else []
        missing = [str(e.get("local_path")) for e in entries
                   if isinstance(e, dict) and not Path(str(e.get("local_path", ""))).is_dir()]
        report.add("system-data.registry", "workspaces.json",
                   WARN if missing else OK,
                   f"{len(entries)} registered Workspace(s)"
                   + (f"; missing on disk: {', '.join(missing[:5])}" if missing else ""),
                   "Remove the stale entry through the app and re-open the folder." if missing else None)

    tool = _native_tool(system / "agent" / "tools")
    report.add("system-data.tools", "Provisioned native tool",
               OK if tool else WARN,
               f"{tool} ({_size(tool)})" if tool else f"missing under {system / 'agent' / 'tools'}",
               None if tool else "Relaunch the app: it re-provisions the managed bin directory.")

    log = system / "logs" / "backend.log"
    if log.is_file():
        lines = _tail_lines(log, 400)
        hits = [line for line in lines if any(token in line.lower() for token in _ERROR_TOKENS)]
        report.add("system-data.logs", "Backend log",
                   WARN if hits else OK,
                   f"{_size(log)}; {len(hits)} error-ish line(s) in the last 400"
                   + (f"; last: {hits[-1][-200:]}" if hits else ""))
    else:
        report.add("system-data.logs", "Backend log", SKIP, f"no {log} yet")


def _resolve_project_root(backend: Path | None, explicit: Path | None, assets: Path | None) -> Path | None:
    """Where pyproject.toml / uv.lock / src live.

    A packaged install keeps the project in the read-only `resources/backend`
    beside the assets while the data root holds the runtime state; a
    development tree keeps it in the sibling `tdease-backend` checkout while the
    data root is `<repo>/data/backend`.
    """
    candidates = [explicit] if explicit else []
    if backend is not None:
        candidates += [
            backend,
            backend.parent / "tdease-backend",
            backend.parent.parent / "tdease-backend",
        ]
    if assets is not None:
        candidates.append(assets.parent.parent / "backend")
    for candidate in candidates:
        if candidate and (candidate / "pyproject.toml").is_file():
            return candidate
    return None


def check_backend(report: Report, backend: Path | None, project: Path | None,
                  lifecycle: Path | None = None) -> None:
    if backend is None:
        report.add("backend", "Backend data root", SKIP, "not located")
        return
    if not backend.is_dir():
        report.add("backend", "Backend data root", FAIL, f"not a directory: {backend}",
                   "The Backend creates it at first start; a wrong --data-root is the usual cause.")
        return
    config, error = _read_json(backend / "config.json")
    report.add("backend.config", "config.json", WARN if error else OK,
               f"unreadable: {error}" if error else
               f"mode={(config or {}).get('deploymentMode')} env={(config or {}).get('environment')}",
               "Delete config.json and relaunch to regenerate it." if error else None)

    present = [name for name in ("database", "nodes", "runs", "store") if (backend / name).is_dir()]
    report.add("backend.layout", "Backend layout", OK if len(present) == 4 else WARN,
               f"{', '.join(present) or 'no expected subdirectory'}",
               None if len(present) == 4 else
               "Restart the app; a data root missing store/ or runs/ is usually a wrong --data-root.")

    for label, root in (("project .venv", project / ".venv" if project else None),
                        ("runtime/venv", backend / "runtime" / "venv")):
        if root is None or not root.is_dir():
            continue
        command = root / ("Scripts" if os.name == "nt" else "bin") / (
            "tdease-backend.exe" if os.name == "nt" else "tdease-backend")
        report.add("backend.env", "Backend environment",
                   OK if command.is_file() else WARN,
                   f"{root} ({label}); console script "
                   + ("present" if command.is_file() else "MISSING"),
                   None if command.is_file() else
                   "Relaunch the app: it re-runs `uv sync --locked --no-dev`. If that "
                   "fails, delete this environment directory and relaunch.")
        break
    else:
        report.add("backend.env", "Backend environment", WARN,
                   "no .venv or runtime/venv found",
                   "Relaunch the app: it provisions the environment on first start "
                   "(needs network or a warm uv cache).")

    check_index(report, backend / "database" / "db.sqlite3")
    check_nodes(report, backend / "nodes")
    _check_stale_catalog_sibling(report, backend.parent / "nodes")
    check_runs(report, backend / "runs", project, lifecycle)


def check_index(report: Report, db: Path) -> None:
    if not db.is_file():
        report.add("backend.index", "Run index", WARN, f"no {db}",
                   "The Backend rebuilds it from runs/*/events.jsonl at next start.")
        return
    try:
        connection = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
        try:
            integrity = connection.execute("PRAGMA quick_check").fetchone()
            tables = {row[0] for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type='table'")}
            rows = connection.execute("SELECT COUNT(*) FROM runs").fetchone()[0] if "runs" in tables else 0
        finally:
            connection.close()
    except sqlite3.Error as exc:
        report.add("backend.index", "Run index", WARN, f"unreadable: {exc}",
                   "Delete db.sqlite3 and restart the app: it is a rebuildable projection "
                   "of events.jsonl, never the truth.")
        return
    healthy = bool(integrity) and integrity[0] == "ok"
    report.add("backend.index", "Run index", OK if healthy else WARN,
               f"quick_check={integrity[0] if integrity else '?'}, {len(tables)} table(s), {rows} run row(s)",
               None if healthy else "Delete db.sqlite3 and restart the app; it rebuilds from the journals.")


def _check_stale_catalog_sibling(report: Report, stale: Path) -> None:
    """A catalog left at <dataRoot>/nodes (beside backend/) is from the retired
    sibling layout. It is never read by the Backend and silently confuses
    repairs; surface it instead of merging anything automatically."""
    if stale.is_dir():
        report.add("backend.nodes.stale-sibling", "Stale node catalog location", WARN,
                   f"a second catalog exists at {stale}; the Backend reads backend/nodes only",
                   f"Delete the stale folder (or move packages you edited into the live catalog) after "
                   f"checking its contents: {stale}")


def check_nodes(report: Report, nodes: Path) -> None:
    if not nodes.is_dir():
        report.add("backend.nodes", "Node catalog", FAIL, f"no {nodes}",
                   "Restart the app: the shipped Node packages are re-seeded into it.")
        return
    packages = [p for p in _list_dirs(nodes) if not p.name.startswith("_")]
    without_rules = [p.name for p in packages if not (p / "rules.smk").is_file()]
    template = (nodes / "_template").is_dir()
    status = OK if packages and not without_rules and template else WARN
    report.add("backend.nodes", "Node catalog", status,
               f"{len(packages)} package(s); _template {'present' if template else 'MISSING'}"
               + (f"; no rules.smk: {', '.join(without_rules[:5])}" if without_rules else ""),
               None if status == OK else
               "Restart the app to re-seed deleted shipped packages; a user-edited package "
               "is never overwritten, so repair it deliberately.")


def check_runs(report: Report, runs_dir: Path, project: Path | None, explicit: Path | None) -> None:
    if not runs_dir.is_dir():
        report.add("backend.runs", "Runs", SKIP, f"no {runs_dir} yet")
        return
    lifecycle_path = explicit or (project / "src" / "contracts" / "run_lifecycle.json" if project else None)
    lifecycle, _ = _read_json(lifecycle_path) if lifecycle_path else (None, "missing")
    statuses = (lifecycle or {}).get("statuses") or {}
    terminals = set((lifecycle or {}).get("terminal_events") or ())
    note = "" if statuses else " (lifecycle declaration not found; statuses are raw event names)"
    run_dirs = _list_dirs(runs_dir)
    findings: list[str] = []
    worst = OK
    has_open_range = False
    for run_dir in run_dirs[:_MAX_RUNS_DETAIL]:
        _, error = _read_json(run_dir / "run.json")
        if error:
            findings.append(f"{run_dir.name}: run.json {error}")
            worst = report.worse(worst, FAIL)
            continue
        events, bad = _read_journal(run_dir / "events.jsonl")
        if not events:
            findings.append(f"{run_dir.name}: no readable journal")
            worst = report.worse(worst, FAIL)
            continue
        status = None
        open_range = False
        for event in events:
            event_type = str(event.get("type") or "")
            if event_type in statuses:
                status = statuses[event_type]
            if event_type == "run.started":
                open_range = True
            elif event_type in terminals:
                open_range = False
        notes = [f"{status or events[-1].get('type') or 'unknown'} @seq {events[-1].get('seq')}"]
        if bad:
            notes.append(f"{bad} unparseable line(s)")
        if open_range:
            has_open_range = True
            process, _ = _read_json(run_dir / ".runtime" / "process.json")
            pid = int((process or {}).get("engine_pid") or 0)
            alive = _pid_alive(pid)
            notes.append(f"OPEN range, engine pid {pid} " + ("alive" if alive else "dead"))
            if not alive:
                worst = report.worse(worst, WARN)
        if status == "SUCCEEDED" and not (run_dir / "results.manifest.json").is_file():
            notes.append("SUCCEEDED without a results manifest")
            worst = report.worse(worst, WARN)
        findings.append(f"{run_dir.name}: " + "; ".join(notes))
    if len(run_dirs) > _MAX_RUNS_DETAIL:
        findings.append(f"… and {len(run_dirs) - _MAX_RUNS_DETAIL} more Run(s) not inspected")
    report.add("backend.runs", "Run journals", worst if findings else OK,
               f"{len(run_dirs)} Run(s){note}\n      " + "\n      ".join(findings) if findings
               else f"none{note}",
               "An OPEN range with a dead engine is recovered at the next Backend start; "
               "`runs.resume` continues it. Never hand-write events." if has_open_range else
               ("A SUCCEEDED Run without a results manifest produced no returned files."
                if worst != OK else None))


# ── roots ─────────────────────────────────────────────────────────────────────


def _looks_like_data_root(path: Path) -> bool:
    return (path / "backend").is_dir() or (path / "system").is_dir()


def _probe_data_root(workspace: Path | None) -> tuple[Path | None, list[str]]:
    candidates: list[Path] = []
    if workspace is not None:
        candidates.append(workspace.parent)
    appdata = os.environ.get("APPDATA")
    if appdata:
        candidates.extend([Path(appdata) / "TDEaseNexus", Path(appdata) / "TDEase"])
    local = os.environ.get("LOCALAPPDATA")
    if local:
        candidates.append(Path(local) / "TDEaseNexus")
    probed = [str(p) for p in candidates]
    for candidate in candidates:
        if _looks_like_data_root(candidate):
            return candidate, probed
    return None, probed


def _probe_system_root(data_root: Path | None) -> Path | None:
    candidates: list[Path] = []
    if data_root is not None:
        candidates.append(data_root / "system")
    appdata = os.environ.get("APPDATA")
    if appdata:
        candidates.extend([Path(appdata) / "TDEase" / "system", Path(appdata) / "TDEaseNexus" / "system"])
    for candidate in candidates:
        if candidate.is_dir():
            return candidate
    return None


def _probe_system_assets(workspace: Path | None, data_root: Path | None) -> Path | None:
    candidates: list[Path] = []
    if data_root is not None:
        candidates.append(data_root.parent / "app" / "assets")
        candidates.append(data_root.parent / "resources" / "app" / "assets")
    if workspace is not None:
        # Walk up a few levels for a development checkout.
        directory = workspace.resolve()
        for _ in range(4):
            candidates.append(directory / "tdease-ai-capability" / "system")
            if directory.parent == directory:
                break
            directory = directory.parent
    for candidate in candidates:
        if (candidate / "skills").is_dir():
            return candidate
    return None


def build_report(args: argparse.Namespace) -> Report:
    report = Report()
    workspace = Path(args.workspace).expanduser() if args.workspace else None
    if workspace is None and os.environ.get("TDEASE_WORKSPACE"):
        workspace = Path(os.environ["TDEASE_WORKSPACE"])
    if workspace is None:
        workspace = Path.cwd()

    data_root = Path(args.data_root).expanduser() if args.data_root else None
    probed: list[str] = []
    if data_root is None:
        data_root, probed = _probe_data_root(workspace)
    backend = Path(args.backend_root).expanduser() if args.backend_root else (
        data_root / "backend" if data_root else None)

    system = Path(args.system_root).expanduser() if args.system_root else _probe_system_root(data_root)

    assets = Path(args.system_assets).expanduser() if args.system_assets else None
    if assets is None and os.environ.get("TDEASE_SYSTEM_ASSETS"):
        assets = Path(os.environ["TDEASE_SYSTEM_ASSETS"])
    if assets is None:
        assets = _probe_system_assets(workspace, data_root)

    project = _resolve_project_root(backend, Path(args.backend_project).expanduser()
                                    if args.backend_project else None, assets)
    lifecycle = Path(args.lifecycle).expanduser() if args.lifecycle else None

    # Resolve once, here: every check reports absolute paths, so a report is
    # readable regardless of the directory it was produced from.
    workspace = workspace.resolve() if workspace else None
    data_root = data_root.resolve() if data_root else None
    system = system.resolve() if system else None
    backend = backend.resolve() if backend else None
    project = project.resolve() if project else None
    assets = assets.resolve() if assets else None

    report.roots = {
        "workspace": str(workspace) if workspace else "",
        "data_root": str(data_root) if data_root else "",
        "system_root": str(system) if system else "",
        "backend_root": str(backend) if backend else "",
        "backend_project": str(project) if project else "",
        "system_assets": str(assets) if assets else "",
    }
    if data_root is None:
        report.add("roots.data", "Data root", WARN,
                   "not found; probed: " + (", ".join(probed) or "nothing"),
                   "Pass --data-root <the TDEaseNexus-Data folder beside the install, or the "
                   "folder holding system/ and backend/ in a development tree>.")

    check_workspace(report, workspace)
    check_system_assets(report, assets)
    check_system_data(report, system)
    check_backend(report, backend, project, lifecycle)
    return report


# ── output ────────────────────────────────────────────────────────────────────


def render(report: Report) -> str:
    lines = ["TDEaseNexus doctor", ""]
    for label, value in report.roots.items():
        lines.append(f"  {label.replace('_', ' '):<14} {value or '-'}")
    lines.append("")
    for check in report.checks:
        lines.append(f"{MARK[check.status]} {check.id:<22} {check.detail}")
        if check.fix:
            lines.append(f"       {'':<22} fix: {check.fix}")
    counts = {status: sum(1 for c in report.checks if c.status == status) for status in RANK}
    lines.append("")
    lines.append(f"{counts[FAIL]} fail, {counts[WARN]} warn, {counts[OK]} ok, {counts[SKIP]} skipped")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workspace", help="Workspace root (default: $TDEASE_WORKSPACE, else cwd)")
    parser.add_argument("--data-root", help="Writable data root holding system/ and backend/")
    parser.add_argument("--system-root", help="Writable System root (default: <data-root>/system, else the profile)")
    parser.add_argument("--backend-root", help="Backend root (default: <data-root>/backend)")
    parser.add_argument("--backend-project", help="Backend project root holding pyproject.toml/uv.lock/src")
    parser.add_argument("--system-assets", help="Read-only shipped assets (default: $TDEASE_SYSTEM_ASSETS)")
    parser.add_argument("--lifecycle", help="Path to run_lifecycle.json (default: <backend-project>/src/contracts/)")
    parser.add_argument("--json", action="store_true", help="Machine-readable report")
    args = parser.parse_args(argv)
    report = build_report(args)
    if args.json:
        print(json.dumps({
            "roots": report.roots,
            "checks": [check.__dict__ for check in report.checks],
            "exit_code": report.exit_code,
        }, ensure_ascii=False, indent=2))
    else:
        print(render(report))
    return report.exit_code


if __name__ == "__main__":
    raise SystemExit(main())
