"""Quiet, strict chart-forge compiler invocation for domain adapters."""
from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

FORGE = Path(__file__).resolve().parents[1] / "chart-forge"


class ChartForgeError(RuntimeError):
    pass


def compile_html(
    example: str,
    data_path: str | Path,
    out_path: str | Path,
    *,
    no_manifest: bool = False,
) -> Path:
    """Compile one example without leaking build diagnostics to adapter stdout."""
    out = Path(out_path)
    command = [
        sys.executable,
        str(FORGE / "scripts" / "build.py"),
        "--example", example,
        "--data", str(data_path),
        "--out", str(out),
    ]
    if no_manifest:
        command.append("--no-manifest")
    completed = subprocess.run(command, capture_output=True, text=True)
    if completed.returncode:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise ChartForgeError(f"chart-forge build failed for {example}: {detail}")
    return out


def emit_html(example: str, payload: dict, plot_path: str | Path, *,
              source_default: str | None = None) -> tuple[Path, Path]:
    """The adapter exit pipe: meta provenance defaults -> ``<stem>.data.json``
    -> one self-contained HTML. Returns (html_path, contract_path).

    Every figure adapter funnels through here, so the contract sidecar naming
    and the meta provenance keys have one implementation."""
    plot = Path(plot_path)
    plot.parent.mkdir(parents=True, exist_ok=True)
    stem = plot.with_suffix("") if plot.suffix else plot
    contract_path = Path(str(stem) + ".data.json")
    payload.setdefault("meta", {})
    if source_default:
        payload["meta"].setdefault("source", source_default)
    payload["meta"].setdefault("generated",
                               datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"))
    payload["meta"]["tool"] = "chart-forge"
    contract_path.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    out_html = compile_html(example, contract_path, Path(str(stem) + ".html"), no_manifest=True)
    return out_html, contract_path
