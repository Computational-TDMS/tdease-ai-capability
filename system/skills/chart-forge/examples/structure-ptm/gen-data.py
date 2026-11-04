#!/usr/bin/env python3
"""Generate structure-ptm/data.json — a protein structure with modification
sites to annotate in 3D.

--pdb takes a PDB file (AlphaFold model or RCSB entry; use the
structure-evidence skill to fetch+cache one: `uv run structure_figure.py
--accession EXAMP1 ...` caches under its reference/structures/). Sites are
--site RESI:LABEL (repeatable); resi MUST be the structure's own numbering —
for AlphaFold that is UniProt 1:1, so chain numbering (e.g. protein A
S1 = UniProt 2) needs the +1 offset. Colours are assigned from the categorical
palette in label order. PDB text is embedded ATOM/TER-only: the figure stays
self-contained and offline."""
import argparse
import json
from pathlib import Path

DEMO_PDB = Path(r"C:\Users\ray\Desktop\TDEase-zero\tdease-ai-capability\system\skills\structure-evidence\scripts\reference\structures\AF-EXAMP1.pdb")
DEMO_SITES = [
    "89:ModA@K47 — adjudicated (+77.9269)",
    "21:K20me2 (+28.0313)",
    "85:M84ox (+15.9949)",
    "2:N-term acetyl (mature S1)",
]
# AlphaFold confidence bins (standard colours) — pLDDT rides in the PDB
# B-factor column; gen-data bins CA atoms so the viewer never guesses.
PLDDT_BINS = [
    (90, 999, "#0053D6", "≥90 very high"),
    (70, 90, "#65CBF3", "70–90 confident"),
    (50, 70, "#FFDB13", "50–70 low"),
    (0, 50, "#FF7D45", "<50 very low"),
]


def plddt_bins(atoms: list[str]) -> list[dict]:
    resi_to_plddt = {}
    for line in atoms:
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        try:
            resi_to_plddt[int(line[22:26])] = float(line[60:66])
        except ValueError:
            continue
    out = []
    for lo, hi, color, label in PLDDT_BINS:
        resis = sorted(r for r, p in resi_to_plddt.items() if lo <= p < hi)
        if resis:
            out.append({"lo": lo, "hi": hi, "color": color, "label": label, "resis": resis})
    return out


def main() -> None:
    here = Path(__file__).parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb", type=Path, default=DEMO_PDB)
    ap.add_argument("--site", action="append", default=[], metavar="RESI:LABEL")
    ap.add_argument("--out", type=Path, default=here / "data.json")
    args = ap.parse_args()
    site_specs = args.site or list(DEMO_SITES)

    atoms = [l for l in args.pdb.read_text(encoding="utf-8").splitlines()
             if l.startswith(("ATOM", "TER"))]
    sites = []
    for spec in site_specs:
        resi_s, _, label = spec.partition(":")
        if not (resi_s.strip().lstrip("-").isdigit() and label):
            raise SystemExit(f"--site needs RESI:LABEL, got {spec!r}")
        sites.append({"resi": int(resi_s), "label": label})

    payload = {
        "data": {"pdb": "\n".join(atoms)},
        "meta": {
            "title": f"Protein A (fictional) — {len(sites)} modification sites in 3D",
            "subtitle": ("AlphaFold AF-EXAMP1 · grey cartoon · site spheres in family colours, "
                         "numbered = K47 adjudicated · click a residue to inspect, "
                         "legend chips zoom to sites"),
            "source": (f"{args.pdb.name} embedded verbatim (AlphaFold v6, UniProt EXAMP1, "
                       "numbering 1:1 — chain numbering as shipped)"),
            "transform": "ATOM/TER records only; sites resi-mapped, no coordinate changes; pLDDT binned from CA B-factors",
            "sites": sites,
            "plddt": plddt_bins(atoms),
            "n": sum(1 for l in atoms if l.startswith("ATOM")),
        },
    }
    args.out.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    print(f"structure data: {payload['meta']['n']} atoms, sites at "
          f"{[s['resi'] for s in sites]} -> {args.out} ({args.out.stat().st_size/1024:.0f} KB)")


if __name__ == "__main__":
    main()
