#!/usr/bin/env python3
"""Acceptance audit for every shipped chart-forge template."""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

FORGE = Path(__file__).resolve().parent.parent
REQUIRED = ("chart.js", "data.json", "build.json")


def run(command: list[str], *, cwd: Path | None = None) -> None:
    completed = subprocess.run(command, cwd=cwd, text=True, capture_output=True)
    if completed.returncode:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise SystemExit(f"command failed ({completed.returncode}): {' '.join(command)}\n{detail}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--browser", action="store_true",
                        help="open every built figure in Chromium and export it")
    parser.add_argument("--keep", type=Path,
                        help="copy built HTML and exports to this directory for inspection")
    args = parser.parse_args()

    examples = sorted(path for path in (FORGE / "examples").iterdir() if path.is_dir())
    failures = {
        path.name: [name for name in REQUIRED if not (path / name).is_file()]
        for path in examples
    }
    failures = {name: missing for name, missing in failures.items() if missing}
    if failures:
        detail = "; ".join(f"{name}: {', '.join(items)}" for name, items in failures.items())
        raise SystemExit(f"incomplete example package(s): {detail}")

    with tempfile.TemporaryDirectory(prefix="chart-forge-audit-") as temp_name:
        temp = Path(temp_name)
        outputs: list[Path] = []
        for example in examples:
            output = temp / f"{example.name}.html"
            run([
                sys.executable, str(FORGE / "scripts" / "build.py"),
                "--example", example.name,
                "--data", str(example / "data.json"),
                "--out", str(output),
            ])
            text = output.read_text(encoding="utf-8")
            if "/*__" in text:
                raise SystemExit(f"{example.name}: unreplaced template token")
            manifest_path = output.with_suffix(".manifest.json")
            if not manifest_path.is_file():
                raise SystemExit(f"{example.name}: missing build manifest")
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            if "data_csv" in manifest:
                csv_path = manifest_path.parent / manifest["data_csv"]
                if not csv_path.is_file():
                    raise SystemExit(f"{example.name}: manifest references missing {csv_path.name}")
            outputs.append(output)

        hostile_data = temp / "hostile.data.json"
        hostile_data.write_text(json.dumps({
            "data": {"rows": [{"label": "</script><script>window.BAD=1</script>", "value": 1}]},
            "meta": {"title": "injection acceptance", "source": "audit fixture"},
        }), encoding="utf-8")
        hostile_output = temp / "injection-acceptance.html"
        run([
            sys.executable, str(FORGE / "scripts" / "build.py"),
            "--example", "categorical-bars", "--data", str(hostile_data),
            "--out", str(hostile_output),
        ])
        if "</script><script>window.BAD" in hostile_output.read_text(encoding="utf-8"):
            raise SystemExit("HTML data injection was not escaped")
        outputs.append(hostile_output)

        export_dir = temp / "export"
        if args.browser:
            node = shutil.which("node")
            if not node:
                raise SystemExit("--browser requires node on PATH")
            command = [node, str(FORGE / "scripts" / "export.mjs")]
            for output in outputs:
                command.extend(("--in", str(output)))
            command.extend(("--out", str(export_dir)))
            run(command)
            pngs = list(export_dir.glob("*.png"))
            if len(pngs) != len(outputs):
                raise SystemExit(f"browser audit exported {len(pngs)}/{len(outputs)} PNG files")
            undersized = [path.name for path in pngs if path.stat().st_size < 1024]
            if undersized:
                raise SystemExit(f"undersized PNG export(s): {', '.join(undersized)}")
            for stem in ("feature-map", "structure-ptm"):
                if (export_dir / f"{stem}.svg").exists():
                    raise SystemExit(f"{stem}: canvas/hybrid figure exported an incomplete SVG")
            for stem in ("categorical-bars", "fragment-annotation"):
                if not (export_dir / f"{stem}.svg").is_file():
                    raise SystemExit(f"{stem}: pure-SVG figure did not export SVG")

        if args.keep:
            target = args.keep.resolve()
            target.mkdir(parents=True, exist_ok=True)
            for output in outputs:
                shutil.copy2(output, target / output.name)
                for sidecar in output.parent.glob(f"{output.stem}.*"):
                    if sidecar != output:
                        shutil.copy2(sidecar, target / sidecar.name)
            if export_dir.is_dir():
                shutil.copytree(export_dir, target / "export", dirs_exist_ok=True)

        mode = "build+browser" if args.browser else "build"
        print(f"OK chart-forge {mode}: {len(examples)} templates + injection acceptance")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
