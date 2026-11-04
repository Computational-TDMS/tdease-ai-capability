#!/usr/bin/env python3
"""chart-forge build: inject tokens + d3 + data + chart code into base.html.

Produces a self-contained HTML scientific asset (offline-loadable, archival).
Usage:
  python scripts/build.py --example categorical-bars --data data.json -o fig.html

Data contract: JSON {data: ..., meta: {...}} or CSV (rows -> array of objects,
numeric strings auto-typed; meta from --meta key=value flags).
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

FORGE = Path(__file__).resolve().parent.parent


def fail(message: str) -> None:
    raise SystemExit(message)


def run_export(out: Path, export_arg: str) -> None:
    """--export passthrough: DIR resolves against the output file's directory
    ('export' stays the repo convention), absolute paths honoured as-is."""
    import shutil
    node = shutil.which("node")
    if not node:
        fail("export requested but node is not on PATH")
    p = Path(export_arg)
    d = p if p.is_absolute() else out.parent / export_arg
    subprocess.run([node, str(FORGE / "scripts" / "export.mjs"),
                    "--in", str(out), "--out", str(d)], check=True)

# Capability-package layout: skills/<skill>/scripts/build.py sits beside skills/_shared/.
# The complete runtime bundle must provide that shared provenance contract.
sys.path.insert(0, str(FORGE.parent))
try:
    from _shared.manifest import write_manifest
except ImportError:
    write_manifest = None


def load_data(path: Path, meta_overrides: dict) -> dict:
    if path.suffix.lower() == ".json":
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(payload, dict) or "data" not in payload:
            fail(f"{path}: JSON contract must be an object containing 'data'")
    elif path.suffix.lower() == ".csv":
        with path.open(encoding="utf-8-sig", newline="") as fh:
            rows = list(csv.DictReader(fh))
        typed = []
        for row in rows:
            out = {}
            for k, v in row.items():
                if v is None or v == "":
                    out[k] = None
                elif re.fullmatch(r"-?\d+", v):
                    out[k] = int(v)
                elif re.fullmatch(r"-?\d*\.\d+([eE][+-]?\d+)?", v):
                    out[k] = float(v)
                else:
                    out[k] = v
            typed.append(out)
        payload = {"data": typed, "meta": {}}
    else:
        fail(f"unsupported data format: {path.suffix} (use .json or .csv)")

    meta = payload.setdefault("meta", {})
    if not isinstance(meta, dict):
        fail(f"{path}: 'meta' must be an object")
    meta.update(meta_overrides)
    meta.setdefault("source", path.name)
    meta.setdefault("generated", datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"))
    meta.setdefault("tool", "chart-forge")
    if (meta.get("source") == path.name
            and isinstance(payload["data"], list)
            and payload["data"]
            and not meta.get("n")):
        meta["n"] = len(payload["data"])
    return payload


def check_grid(payload: dict, path: Path) -> None:
    """Family-8/11 contract (data = {marks, conds, cells}): cells must cover the
    declared marks×conds grid exactly once. Structural absence must be explicit
    (occ: null); an omitted cell is indistinguishable from a typo'd key and once
    crashed the figure at render time ("Cannot read properties of undefined")."""
    data = payload["data"]
    if not isinstance(data, dict):
        return
    if not all(isinstance(data.get(k), list) for k in ("marks", "conds", "cells")):
        return
    declared = {(m, c) for m in data["marks"] for c in data["conds"]}
    seen: dict[tuple[str, str], int] = {}
    for cell in data["cells"]:
        if not isinstance(cell, dict) or not {"mark", "cond"} <= set(cell):
            fail(f"{path}: cell without mark/cond: {cell!r}")
        key = (cell["mark"], cell["cond"])
        if key not in declared:
            fail(f"{path}: cell {key} references a mark/cond outside the declared "
                 f"marks×conds grid")
        seen[key] = seen.get(key, 0) + 1
    missing = sorted(declared - set(seen))
    if missing:
        sample = ", ".join(f"{{'mark': {m!r}, 'cond': {c!r}, 'occ': None}}"
                           for m, c in missing[:3])
        fail(f"{path}: cells cover {len(seen)}/{len(declared)} of the marks×conds "
             f"grid; missing {len(missing)} (e.g. {sample}). Emit structural gaps "
             f"explicitly with occ: null")
    dupes = sorted(k for k, v in seen.items() if v > 1)
    if dupes:
        fail(f"{path}: duplicate cells for {dupes[:3]}")


def load_example(name: str) -> tuple[Path, dict]:
    if not re.fullmatch(r"[a-z0-9][a-z0-9-]*", name):
        fail(f"invalid example name: {name!r}")
    examples_root = (FORGE / "examples").resolve()
    exdir = (examples_root / name).resolve()
    if exdir.parent != examples_root or not exdir.is_dir():
        fail(f"unknown example: {name}")

    spec_path = exdir / "build.json"
    if not spec_path.is_file():
        fail(f"example {name!r} is missing required build.json")
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    if not isinstance(spec, dict):
        fail(f"{spec_path}: build spec must be an object")
    unknown = sorted(set(spec) - {"bundle", "journal", "scripts"})
    if unknown:
        fail(f"{spec_path}: unknown field(s): {', '.join(unknown)}")
    if spec.get("bundle", "d3") not in {"d3", "3dmol", "none"}:
        fail(f"{spec_path}: invalid bundle {spec.get('bundle')!r}")
    if spec.get("journal") not in {None, "nature", "cell", "analytical-chemistry"}:
        fail(f"{spec_path}: invalid journal {spec.get('journal')!r}")
    scripts = spec.get("scripts", [])
    if not isinstance(scripts, list) or not all(isinstance(item, str) for item in scripts):
        fail(f"{spec_path}: scripts must be a string array")
    forge_root = FORGE.resolve()
    for item in scripts:
        resolved = (exdir / item).resolve()
        if forge_root not in resolved.parents or not resolved.is_file():
            fail(f"{spec_path}: script must resolve inside chart-forge: {item}")
    return exdir, spec


def json_for_html(payload: dict) -> str:
    """Serialize JSON without allowing source text to close the script element."""
    return (json.dumps(payload, ensure_ascii=False)
            .replace("&", "\\u0026")
            .replace("<", "\\u003c")
            .replace(">", "\\u003e"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--template", default=str(FORGE / "templates" / "base.html"))
    ap.add_argument("--example", metavar="NAME",
                    help="gallery template: resolve chart.js and build.json "
                         "(scripts/bundle/journal defaults) from examples/NAME — "
                         "the one-line form; explicit flags override the spec")
    ap.add_argument("--data", required=True, help="data.json ({data, meta}) or data.csv")
    ap.add_argument("--chart", help="chart .js file (reads window.DATA); "
                                    "implied by --example")
    ap.add_argument("--tokens", default=str(FORGE / "templates" / "tokens.css"))
    ap.add_argument("--bundle", default=None, choices=["d3", "3dmol", "none"],
                    help="JS library to inline at the /*__D3__*/ placeholder "
                         "(default: example build.json, else d3)")
    ap.add_argument("--journal", default=None,
                    choices=["nature", "cell", "analytical-chemistry"],
                    help="journal style preset inlined as window.JOURNAL before the chart")
    ap.add_argument("--lib", default=None, help="explicit lib file path (overrides --bundle)")
    ap.add_argument("--no-micro", action="store_true",
                    help="skip the templates/micro.js scaffolding injection")
    ap.add_argument("--export", nargs="?", const="export", metavar="DIR",
                    help="after building, export PNG/SVG via export.mjs into DIR "
                         "(default: <out-dir>/export)")
    ap.add_argument("--out", "-o", required=True)
    ap.add_argument("--title"), ap.add_argument("--subtitle")
    ap.add_argument("--meta", action="append", default=[],
                    help="meta key=value, repeatable (e.g. --meta thresholds='{\"padj\":0.05}')")
    ap.add_argument("--script", action="append", default=[],
                    help="extra .js file(s) inlined before the chart script; 'export const X' rewritten to 'const X'")
    ap.add_argument("--no-manifest", action="store_true",
                    help="skip the _shared provenance sidecar even when available")
    args = ap.parse_args()
    if write_manifest is None and not args.no_manifest:
        fail("_shared.manifest is unavailable; use the complete capability bundle "
             "or pass --no-manifest when a domain adapter writes the sidecar")

    # --example: resolve the gallery spec first; explicit flags win over it.
    example_scripts: list[str] = []
    example_spec_path: str | None = None
    if args.example:
        exdir, spec = load_example(args.example)
        if not args.chart:
            args.chart = str(exdir / "chart.js")
        args.bundle = args.bundle or spec.get("bundle", "d3")
        args.journal = args.journal or spec.get("journal")
        example_scripts = [str((exdir / item).resolve()) for item in spec.get("scripts", [])]
        example_spec_path = str(exdir / "build.json")
    args.bundle = args.bundle or "d3"
    if not args.chart:
        fail("--chart (or --example) is required")

    overrides = {}
    if args.title:
        overrides["title"] = args.title
    if args.subtitle:
        overrides["subtitle"] = args.subtitle
    for kv in args.meta:
        if "=" not in kv:
            fail(f"--meta requires key=value, got {kv!r}")
        k, _, v = kv.partition("=")
        if not k:
            fail("--meta key cannot be empty")
        try:
            overrides[k] = json.loads(v)
        except json.JSONDecodeError:
            overrides[k] = v

    data_path = Path(args.data)
    chart_path = Path(args.chart)
    for label, path in (("data", data_path), ("chart", chart_path),
                        ("template", Path(args.template)), ("tokens", Path(args.tokens))):
        if not path.is_file():
            fail(f"{label} file not found: {path}")
    payload = load_data(data_path, overrides)
    check_grid(payload, data_path)
    html = Path(args.template).read_text(encoding="utf-8")
    tokens = Path(args.tokens).read_text(encoding="utf-8")
    lib_path = Path(args.lib) if args.lib else {
        "d3": FORGE / "assets" / "d3" / "vendor" / "d3.min.js",
        "3dmol": FORGE / "assets" / "vendor" / "3Dmol.min.js",
        "none": None,
    }[args.bundle]
    if lib_path and not lib_path.is_file():
        fail(f"library file not found: {lib_path}")
    lib = lib_path.read_text(encoding="utf-8") if lib_path else ""
    extras = []
    if not args.no_micro:
        extras.append((FORGE / "templates" / "micro.js").read_text(encoding="utf-8"))
    if args.journal:
        code = (FORGE / "templates" / "journals" / f"{args.journal}.js").read_text(encoding="utf-8")
        extras.append(re.sub(r"^export const (\w+)", r"const \1", code, flags=re.M))
    for s in [*example_scripts, *args.script]:
        script_path = Path(s)
        if not script_path.is_file():
            fail(f"script file not found: {script_path}")
        code = script_path.read_text(encoding="utf-8")
        code = re.sub(r"^export const (\w+)", r"const \1", code, flags=re.M)
        extras.append(code)
    chart = "\n;\n".join(extras + [Path(args.chart).read_text(encoding="utf-8")])

    # Missing-provider guard: a chart that references a known companion global
    # (component/ramp) without the providing script injected produces a figure
    # that dies on open ("forgeMs2 is not defined"). Fail here, loudly, instead.
    def _strip_js_comments(code: str) -> str:
        code = re.sub(r"/\*.*?\*/", "", code, flags=re.S)
        return re.sub(r"^\s*//.*$", "", code, flags=re.M)

    PROVIDERS = {
        "forgeMs2": "components/ms2-annotation.js",
        "batlow": "templates/crameri-batlow.js",
        "vik": "templates/crameri-vik.js",
        "roma": "templates/crameri-roma.js",
        "lajolla": "templates/crameri-lajolla.js",
        "bilbao": "templates/crameri-bilbao.js",
    }
    provided = "\n".join(extras)
    code_only = _strip_js_comments(Path(args.chart).read_text(encoding="utf-8"))
    if args.bundle == "d3" and 'd3.select("#chart")' not in code_only:
        # Report-mount convention (engineering-contract: the chart "reads
        # window.DATA and renders into #chart"), checked at build time so
        # every compiled D3 figure is report-mountable by construction.
        fail("chart never mounts into #chart; the report compiler cannot embed it")
    for ident, provider in PROVIDERS.items():
        if re.search(rf"\b{ident}\b", code_only) and not re.search(rf"\b{ident}\b", provided):
            fail(f"chart references '{ident}' but its provider was not injected; "
                 f"add {provider} to the example build.json scripts array")

    for placeholder, content in [
        ("/*__TOKENS__*/", tokens),
        ("/*__D3__*/", lib),
        ("/*__DATA__*/", json_for_html(payload)),
        ("/*__CHART__*/", chart),
    ]:
        if placeholder not in html:
            fail(f"template missing placeholder {placeholder}")
        html = html.replace(placeholder, content)

    leftovers = re.findall(r"/\*__[A-Z]+__\*/", html)
    if leftovers:
        fail(f"unreplaced placeholders remain: {leftovers}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    data = payload["data"]
    shape = (f"list[{len(data)}]" if isinstance(data, list)
             else f"object[{len(data)}]" if isinstance(data, dict)
             else type(data).__name__)
    print(f"OK {out} ({out.stat().st_size / 1024:.0f} KB, data={shape}, "
          f"meta={sorted(k for k in payload['meta'])})")

    if write_manifest is not None and not args.no_manifest:
        rows = payload["data"] if isinstance(payload["data"], list) and payload["data"] else None
        mpath = write_manifest(
            out,
            figure=out.name,
            script="chart-forge/scripts/build.py",
            params={k: payload["meta"][k]
                    for k in ("source", "thresholds", "colour_scale") if k in payload["meta"]},
            inputs=[args.data, args.chart, *example_scripts, *args.script,
                    *([example_spec_path] if example_spec_path else [])],
            caption={"title": payload["meta"].get("title"),
                     "subtitle": payload["meta"].get("subtitle"),
                     "n": len(payload["data"]) if isinstance(payload["data"], list) else None},
            plot_data=rows,
        )
        print(f"manifest {mpath}")

    if args.export:
        run_export(out, args.export)


if __name__ == "__main__":
    main()
