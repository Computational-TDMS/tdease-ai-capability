"""Single source of truth for locating TDEase read-only System assets.

Every skill script that reads a bundled reference asset (unimod.xml,
knowledge.db, the shared plotting helper) resolves it through this module —
never through a writable System-root variable or a caller-specific path hack.

Resolution order (all functions, nothing at import time):
  1. `TDEASE_SYSTEM_ASSETS` — the **read-only** System assets root
     (dev: `tdease-ai-capability/system`; packaged:
     `resources/app/assets`). This is the product's real asset root; it is
     always preferred.
  2. Repo fallback — the `system/` directory of the toolkit checkout, derived
     from this file's location (`parents[2]`, i.e.
     `<toolkit>/system/reference/scripts` -> `<toolkit>/system`): works when
     scripts run from the toolkit repo or its submodule mount.

This module is the same shared layer as `_plot.py`: it lives in
`reference/scripts/`, scripts locate the directory once via their bootstrap
`sys.path` insert, then `import system_paths`. Import has zero side effects —
every lookup happens inside a function, so a script can import this module in
an env with no `TDEASE_*` variables set.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

# error codes (stable contract, mirrored by each consumer's ToolError)
SYSTEM_ASSETS_ROOT_NOT_FOUND = "SYSTEM_ASSETS_ROOT_NOT_FOUND"
REFERENCE_ROOT_NOT_FOUND = "REFERENCE_ROOT_NOT_FOUND"
REFERENCE_MISSING = "REFERENCE_MISSING"


class PathsError(Exception):
    def __init__(self, code: str, message: str, **detail):
        super().__init__(message)
        self.code = code
        self.message = message
        self.detail = detail


def _env_dir(name: str) -> Path | None:
    value = os.environ.get(name)
    if not value:
        return None
    candidate = Path(value).expanduser()
    return candidate if candidate.is_dir() else None


def system_assets_root() -> Path:
    """The read-only System assets root (skills/, prompts/, reference/)."""
    candidate = _env_dir("TDEASE_SYSTEM_ASSETS")
    if candidate is not None:
        return candidate
    # Repo fallback: this file lives at <toolkit>/system/reference/scripts/
    fallback = Path(__file__).resolve().parents[2]
    if (fallback / "skills").is_dir():
        return fallback
    raise PathsError(
        SYSTEM_ASSETS_ROOT_NOT_FOUND,
        "set TDEASE_SYSTEM_ASSETS to the System assets root",
        env_tried=("TDEASE_SYSTEM_ASSETS",),
    )


def reference_dir() -> Path:
    """The reference/ directory under the assets root."""
    candidate = _env_dir("TDEASE_SYSTEM_ASSETS")
    if candidate is not None:
        ref = candidate / "reference"
        if ref.is_dir():
            return ref
    fallback = Path(__file__).resolve().parents[2] / "reference"
    if fallback.is_dir():
        return fallback
    raise PathsError(
        REFERENCE_ROOT_NOT_FOUND,
        "cannot locate the System reference directory",
        env_tried=("TDEASE_SYSTEM_ASSETS",),
    )


def unimod_xml() -> Path:
    """The bundled unimod.xml under the reference root."""
    path = reference_dir() / "unimod.xml"
    if not path.is_file():
        raise PathsError(REFERENCE_MISSING, f"{path} does not exist", path=str(path))
    return path


def add_reference_scripts() -> None:
    """Idempotently put reference/scripts on sys.path (for importing _plot)."""
    scripts = reference_dir() / "scripts"
    key = str(scripts)
    if key not in sys.path:
        sys.path.insert(0, key)
