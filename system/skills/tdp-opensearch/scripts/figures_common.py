"""
Shared rendering adapter for the OpenSearch figure modules: contract sidecar
+ one chart-forge build call, and the standard stdout envelope.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from _shared.chart_forge import emit_html


def _load_json(source):
    """A panel input that is either a JSON path or the already-parsed bundle
    object; returns (payload, source display name)."""
    if isinstance(source, (str, Path)):
        return json.loads(Path(source).read_text(encoding="utf-8")), Path(source).name
    return source, "opensearch.json"


def _render(example: str, payload: dict, plot_path: str) -> list[str]:
    """Write the contract beside the plot, then one build call. Returns the
    written paths (the HTML is the deliverable; PNG/SVG live in the page's
    Save buttons)."""
    out_html, contract_path = emit_html(example, payload, plot_path)
    return [str(out_html), str(contract_path)]


def _emit(mode: str, written: list[str]) -> None:
    """Standard adapter stdout: one JSON envelope whose ``plot`` the runtime
    parses; extras beyond figure/plot/written are mode-specific."""
    print(json.dumps({"figure": mode, "plot": written[0], "written": written},
                     ensure_ascii=False))

