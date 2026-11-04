"""Calling the tdeasenexus-toolbox kernel: one naming rule, one subprocess pipe.

There is no gateway and no capability probing: the kernel is a plain CLI and
fails loudly on its own when a command is missing. ``binary`` is the whole
naming contract — an explicit ``--ms-bin`` override wins, then
``TDEASENEXUS_TOOLBOX_BIN``, then the kernel that ships beside the skills in the
assets tree, then PATH (the harness-managed bin dir). ``run_json`` is only
shared subprocess hygiene: UTF-8 text regardless of locale (the zh-CN GBK
default crashed hand-rolled call sites), JSON-envelope parsing, and one
error surface.

Shared by every skill (PEP 723 scripts bootstrap this package via
``sys.path`` — the whole ``system/`` tree ships together).
"""
from __future__ import annotations

import json
import os
import subprocess
from functools import lru_cache
from pathlib import Path


def load_json_arg(flag: str, text: str):
    """Parse agent-supplied JSON exactly as the kernel does: an existing file
    path wins, otherwise the text itself is the JSON (mirror of the toolbox's
    ``parse_json_arg`` — the two CLIs must not diverge on this contract, or a
    path that one accepts becomes a cryptic bare ``JSONDecodeError`` in the
    other's generic error surface)."""
    path = Path(text)
    from_file = path.is_file()
    source = path.read_text(encoding="utf-8") if from_file else text
    try:
        return json.loads(source)
    except json.JSONDecodeError as exc:
        where = f"file {text}" if from_file else "inline JSON"
        raise RuntimeError(
            f"{flag} must be inline JSON or a path to a JSON file ({where} failed to parse: {exc})"
        ) from exc


def binary(requested: str | None = None) -> str:
    """--ms-bin override > TDEASENEXUS_TOOLBOX_BIN > assets tree beside the skills > PATH."""
    if requested:
        return requested
    env = os.environ.get("TDEASENEXUS_TOOLBOX_BIN")
    if env:
        return env
    name = "tdeasenexus-toolbox.exe" if os.name == "nt" else "tdeasenexus-toolbox"
    candidate = Path(__file__).resolve().parents[2] / "tools" / name
    if candidate.is_file():
        return str(candidate)
    return name


@lru_cache(maxsize=4)
def _mass_constants(kernel: str) -> dict:
    """The kernel's numerical authority, once per selected binary/process."""
    doc = run_json([kernel, "schema", "theory"], label="kernel mass constants")
    return doc["result"]["x-mass-constants"]


def mass_constants(requested: str | None = None) -> dict:
    return _mass_constants(binary(requested))


def validate_envelope(doc, *, source):
    if not isinstance(doc, dict) or doc.get("schemaVersion") != 1:
        version = doc.get("schemaVersion") if isinstance(doc, dict) else None
        raise RuntimeError(f"SCHEMA_MISMATCH: {source} returned schemaVersion={version!r}; expected 1")
    return doc


def run_json(cmd, label=None, check=True, timeout=None):
    """Run one kernel command and parse its JSON stdout envelope.

    ``cmd[0]`` is the binary path (see :func:`binary`). ``check=False``
    returns None instead of raising, for call sites that tolerate a failed
    measurement.
    """
    name = label or " ".join(cmd[:2])
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          encoding="utf-8", errors="replace", timeout=timeout)
    if proc.returncode:
        if not check:
            return None
        detail = (proc.stderr or proc.stdout).strip()[-2000:]
        raise RuntimeError(f"{name} failed (exit {proc.returncode}): {detail}")
    try:
        doc = json.loads(proc.stdout)
    except json.JSONDecodeError as exc:
        if not check:
            return None
        raise RuntimeError(f"{name} printed invalid JSON: {exc}") from exc
    try:
        validate_envelope(doc, source=name)
    except RuntimeError:
        if not check:
            return None
        raise
    return doc
