# /// script
# requires-python = ">=3.11"
# ///
"""Build a ProForma model set from a placement table, for spectral measurement.

One transformation: placements in, ProForma strings out. A model is one way of
distributing a fixed set of modifications over a sequence; identical
composition with different placement is exactly what a spectral comparison
needs, and writing these strings by hand is where transcription errors enter.

Usage:
  models_build.py --sequence ACDGSK --model id:pos:Mod[,pos:Mod ...] [--model ...] [--out models.json]

  --sequence STR   the bare chain the positions index into (1-based)
  --model SPEC     id:position:Modification[,position:Modification ...]
                   positions are 1-based; position 0 is the N-terminus and
                   position len+1 the C-terminus (terminal modifications emit
                   the ProForma terminal form [Mod]-SEQ / SEQ-[Mod], which is
                   a different site semantics from a residue-side mass);
                   several modifications on one residue are allowed
                   (R[Mod1][Mod2], the null-control form); modification names
                   must resolve in the kernel ontology when measured
Stdout carries the JSON model array [{"id", "label", "proforma"}] — pipe it
straight into spectral.py measure --models.
"""
from __future__ import annotations

import argparse
import json
import re


def parse_model(spec: str) -> tuple[str, list[tuple[int, str]]]:
    model_id, _, body = spec.partition(":")
    if not model_id or not body:
        raise ValueError(f"model spec must be id:pos:Mod[,pos:Mod ...] — got {spec!r}")
    placements = []
    for part in body.split(","):
        pos_s, _, mod = part.partition(":")
        if not re.fullmatch(r"\d+", pos_s) or not mod:
            raise ValueError(f"placement must be pos:Mod — got {part!r} in {spec!r}")
        placements.append((int(pos_s), mod))
    return model_id, placements


def build(sequence: str, placements: list[tuple[int, str]]) -> str:
    by_pos: dict[int, list[str]] = {}
    for pos, mod in placements:
        if not 0 <= pos <= len(sequence) + 1:
            raise ValueError(f"position {pos} outside sequence of length {len(sequence)} "
                             "(0 = N-terminus, len+1 = C-terminus)")
        by_pos.setdefault(pos, []).append(mod)
    n_term = "".join(f"[{m}]" for m in sorted(by_pos.get(0, [])))
    c_term = "".join(f"[{m}]" for m in sorted(by_pos.get(len(sequence) + 1, [])))
    out = []
    for i, residue in enumerate(sequence, start=1):
        mods = by_pos.get(i)
        out.append(residue + "".join(f"[{m}]" for m in sorted(mods)) if mods else residue)
    core = "".join(out)
    prefix = f"{n_term}-" if n_term else ""
    suffix = f"-{c_term}" if c_term else ""
    return f"{prefix}{core}{suffix}"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sequence", required=True)
    p.add_argument("--model", action="append", required=True)
    p.add_argument("--out")
    a = p.parse_args()

    models = []
    for spec in a.model:
        model_id, placements = parse_model(spec)
        models.append({
            "id": model_id,
            "label": model_id.replace("_", " "),
            "proforma": build(a.sequence, placements),
        })
    seen: dict[str, str] = {}
    for m in models:
        if m["proforma"] in seen:
            raise SystemExit(f"models {seen[m['proforma']]} and {m['id']} build identical ProForma strings")
        seen[m["proforma"]] = m["id"]

    text = json.dumps(models, indent=1)
    print(text)
    if a.out:
        with open(a.out, "w", encoding="utf-8") as fh:
            fh.write(text + "\n")


if __name__ == "__main__":
    main()
