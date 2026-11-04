#!/usr/bin/env python3
"""Generate proteoform-coverage/data.json — per-proteoform sequence-coverage
bands from matched-fragment CSVs (tdease-ms fragment-match output:
ion, number, charge, mz, observed_mz, error_ppm, intensity, label, sequence).

Top-down coverage lives at BACKBONE CLEAVAGE SITES, not residue spans: a c_n
(or b_n) match reads the bond after residue n; a z_m (or y_m) match reads the
bond before residue L-m+1. Each bond cell is 0/1/2/3 = unread / c-side /
z-side / both — so full-ladder spectra still show their holes, and competing
hypotheses diverge at the bonds only one of them explains. Sequence and mods
are read back from the longest fragment's bracket notation, making the CSV
the single source of truth. The renderer only draws.

Each --row "LABEL:TRACK:CSV" is one band row (TRACK names its first spectrum;
a row may gain more tracks later — mergedSites is their union). Default demo:
a fictional dual-hypothesis demo: two fictional models differing at one site."""
import argparse
import csv
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[6]
DEMO = [
    ("model A (ModA@K47)", "ETD 2201", REPO / "examples/fictional-inputs/fictional-annotation-input"),
    ("model B (ModB@K47)", "ETD 2201", REPO / "examples/fictional-inputs/fragment-matches-phospho-match.matches.csv"),
]
MOD_NAMES = {  # delta mass -> short name (Unimod monoisotopic)
    28.0313: "me2", 16.00: "ox", 77.9269: "Br", 79.9663: "phos",
}

NONE, C_SIDE, Z_SIDE, BOTH = 0, 1, 2, 3


def parse_longest(rows: list[dict]) -> tuple[str, list[dict]]:
    """Longest fragment string -> (bare sequence, mods[{pos,name,delta}])."""
    s = max((r["sequence"] for r in rows), key=len)
    seq, mods, pos = [], [], 1
    for tok in re.findall(r"\[[^\]]*\]|.", s):
        if tok.startswith("["):
            if tok == "[Acetyl]":
                mods.append({"pos": 0, "name": "ac", "delta": 42.0106})
            else:
                delta = float(tok[1:-1].lstrip("+"))
                name = next((v for k, v in MOD_NAMES.items() if abs(k - delta) < 0.02), f"+{delta:.2f}")
                mods.append({"pos": pos, "name": name, "delta": delta})
            continue
        seq.append(tok)
        pos += 1
    return "".join(seq), mods


def site_cells(rows: list[dict], L: int) -> list[int]:
    """Bond after residue n (1-based, n = 1..L-1): c/b n → C_SIDE;
    z/y m → bond before residue L-m+1, i.e. bond L-m → Z_SIDE."""
    cells = [NONE] * (L - 1)
    for r in rows:
        n, ion = int(r["number"]), r["ion"]
        bond = n if ion in ("c", "b") else L - n
        if 1 <= bond <= L - 1:
            cells[bond - 1] |= C_SIDE if ion in ("c", "b") else Z_SIDE
    return cells


def union_cells(tracks: list[list[int]]) -> list[int]:
    out = [NONE] * len(tracks[0])
    for t in tracks:
        out = [a | b for a, b in zip(out, t)]
    return out


def main() -> None:
    here = Path(__file__).parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--row", action="append", metavar="LABEL:TRACK:CSV",
                    help="repeatable; defaults to the example dual-hypothesis demo")
    ap.add_argument("--start", type=int, default=1,
                    help="biological numbering of the first residue (default 1 — the "
                         "label convention, e.g. ModA@K47)")
    ap.add_argument("--out", type=Path, default=here / "data.json")
    args = ap.parse_args()

    rows_spec = []
    for spec in (args.row or [f"{l}:{t}:{p}" for l, t, p in DEMO]):
        label, _, rest = spec.partition(":")
        track, _, path = rest.partition(":")
        if not (label and track and path):
            raise SystemExit(f"--row needs LABEL:TRACK:CSV, got {spec!r}")
        rows_spec.append((label, track, Path(path)))

    seq_ref, payload_rows = None, []
    for label, track, path in rows_spec:
        with path.open(encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        seq, mods = parse_longest(rows)
        if seq_ref is None:
            seq_ref = seq
        elif seq != seq_ref:
            raise SystemExit(f"{label}: sequence differs from the first row — rows must share one sequence")
        L = len(seq)
        cells = site_cells(rows, L)
        merged = union_cells([cells])
        payload_rows.append({
            "label": label,
            "mods": mods,
            "mergedSites": merged,
            "coveragePct": round(100 * sum(1 for c in merged if c) / (L - 1), 1),
            "tracks": [{"name": track, "sites": cells}],
            "nFragments": len(rows),
        })

    payload = {
        "data": {
            "sequence": seq_ref,
            "start": args.start,
            "rows": payload_rows,
        },
        "meta": {
            "title": " · ".join(r["label"] for r in payload_rows),
            "subtitle": (f"{len(seq_ref)} residues numbered from {args.start} · thick band = merged backbone "
                         "cleavage sites · thin band = per-spectrum track · ticks = modifications · "
                         "cell colour: c-side / z-side / both"),
            "source": "fictional fragment-match tables (example, ETD c/z)",
            "transform": "bond after residue n read by c/b n; bond L-m read by z/y m; cells 0/1/2/3 = none/c/z/both; mods from fragment bracket notation",
            "n": f"{len(payload_rows)} proteoform rows · {payload_rows[0]['tracks'][0]['name']}",
        },
    }
    args.out.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    for r in payload_rows:
        holes = [i + 1 for i, c in enumerate(r["mergedSites"]) if not c]
        print(f"{r['label']}: {r['nFragments']} fragments, {r['coveragePct']}% bonds read, "
              f"mods at {[m['pos'] for m in r['mods']]}, unread bonds: {holes[:20]}{'...' if len(holes) > 20 else ''}")
    print(f"-> {args.out}")


if __name__ == "__main__":
    main()
