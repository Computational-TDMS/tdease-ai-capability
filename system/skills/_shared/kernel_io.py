"""Reading kernel artifacts from disk: one unwrap, one error vocabulary.

The kernel's CLI envelope is ``{command, params, source, elapsedMs, warnings,
result}`` (rust ``output.rs``). Adapters used to hand-unwrap the durable JSON
artifacts with three divergent strictness levels — a typed error, a silent
``or {}``, a bare ``KeyError``. This module is the single reader for those
artifacts and the python twin of the kernel's own ``read_anchored`` shape
dispatch (bare anchored JSON / ``{envelope, anchored}`` / full CLI envelope).

Scope line: artifacts on DISK go through here. In-memory envelopes returned by
``native.run_json`` are exit-0-gated kernel contract and stay at their call
sites, as do batch protocols that iterate per-item results.
"""
from __future__ import annotations

import json
from pathlib import Path

from .native import validate_envelope


class KernelArtifactError(ValueError):
    """A kernel artifact is missing, unreadable, malformed, or shape-wrong."""


def load_doc(path: str | Path) -> dict:
    """One kernel CLI JSON artifact from disk, IO errors included."""
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except OSError as exc:
        raise KernelArtifactError(f"{path}: cannot read ({exc.strerror or exc})") from exc
    except json.JSONDecodeError as exc:
        raise KernelArtifactError(f"{path}: malformed JSON ({exc})") from exc


def result_of(doc: dict, *, source: str) -> dict:
    """The ``result`` object of one in-memory kernel envelope."""
    try:
        validate_envelope(doc, source=source)
    except RuntimeError as exc:
        raise KernelArtifactError(str(exc)) from exc
    result = doc.get("result") if isinstance(doc, dict) else None
    if not isinstance(result, dict):
        raise KernelArtifactError(f"{source}: no 'result' object in kernel envelope")
    return result


def load_result(path: str | Path) -> dict:
    """``result`` of one kernel CLI JSON artifact on disk."""
    return result_of(load_doc(path), source=str(path))


def _theory_object(doc: dict, *, source: str) -> dict:
    """Unwrap one theory artifact to its payload object: a full CLI envelope
    resolves to its (schema-validated) ``result``; a sliced theory object or a
    bare anchored JSON passes through untouched. Batch callers legitimately
    write slices of one ``theory --models`` envelope, so every theory reader
    below accepts the same three shapes the kernel's ``read_anchored`` does."""
    if isinstance(doc, dict) and isinstance(doc.get("result"), dict):
        return result_of(doc, source=source)
    return doc


def load_anchored(path: str | Path) -> dict:
    """The ``anchored`` envelope of a theory artifact: full CLI envelope,
    sliced theory object, or bare anchored JSON."""
    doc = _theory_object(load_doc(path), source=str(path))
    if isinstance(doc, dict) and isinstance(doc.get("anchored"), dict):
        doc = doc["anchored"]
    if not isinstance(doc, dict) or not isinstance(doc.get("peaks"), list):
        raise KernelArtifactError(f"{path}: no anchored envelope in theory artifact")
    return doc


def load_theory_envelope(path: str | Path) -> dict:
    """The display envelope of a theory artifact: ``envelope`` when present,
    else ``anchored`` — at ``result.`` depth for full envelopes, top level for
    slices, or the bare envelope itself (kernel ``read_anchored`` shapes)."""
    doc = _theory_object(load_doc(path), source=str(path))
    for key in ("envelope", "anchored"):
        if isinstance(doc, dict) and isinstance(doc.get(key), dict):
            return doc[key]
    if isinstance(doc, dict) and isinstance(doc.get("peaks"), list):
        return doc
    raise KernelArtifactError(
        f"{path}: theory artifact has neither 'envelope' nor 'anchored'")
