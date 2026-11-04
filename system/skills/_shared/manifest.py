"""Figure provenance sidecars: sha256 helpers and optional drawn-data CSV.

One implementation for every skill that writes ``<stem>.manifest.json`` +
``<stem>.data.csv`` beside a figure (tdp-opensearch/figures, tdp-copilot
envelope figures, isoforma). No rendering dependencies.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

VERSION = "0.1.0"


def sha256_file(path) -> str:
    """SHA-256 hex digest of a file (strict; used for provenance of inputs)."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_or_none(path: Path) -> str | None:
    try:
        return sha256_file(path)
    except OSError:
        return None


def write_manifest(
    figure_path,
    *,
    figure: str,
    script: str,
    version: str = VERSION,
    params: dict | None = None,
    inputs: list[str] | None = None,
    caption: dict | None = None,
    plot_data: list[dict] | None = None,
    data_columns: list[str] | None = None,
) -> Path:
    """Write ``<stem>.manifest.json`` and, when rows exist, a data CSV.

    The manifest records the figure's provenance (input SHA-256, script
    version, parameters, caption facts) and the exact data drawn, so a figure
    can be reproduced or audited without re-running the analysis.
    """
    figure_path = Path(figure_path)
    stem = figure_path.with_suffix("")
    inputs = inputs or []
    manifest = {
        "figure": figure,
        "script": script,
        "version": version,
        "params": params or {},
        "inputs": [
            {"path": str(path), "sha256": _sha256_or_none(Path(path))}
            for path in inputs
        ],
        "caption": caption or {},
    }
    if plot_data:
        manifest["data_csv"] = f"{stem.name}.data.csv"
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path = stem.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    if plot_data:
        columns = data_columns or list(plot_data[0])
        with stem.with_suffix(".data.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(plot_data)
    return manifest_path
