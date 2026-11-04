"""Family-10 fragment-annotation composite — TopPIC PRSM title + annotated
spectrum + ppm error panel. Fictional fragment table (396 c/z, 15 ppm);
full peptide reconstructed and asserted (same pipeline as spectrum-viewer)."""
import csv
import json
import pathlib
import re

HERE = pathlib.Path(__file__).resolve().parent
SRC = HERE.parent / "fictional-inputs" / "fictional-annotation-input"

strip = lambda s: re.sub(r"\[[^\]]*\]", "", s)
parse_residues = lambda s: [  # annotated sequence -> [{char, mod}]; mod tag attaches to the residue FOLLOWING it
    {"char": ch, "mod": mods[i] if i < len(mods) else None}
    for i, ch in enumerate(strip(s))
] if (mods := re.findall(r"\[([^\]]*)\]", s)) is not None else []

peaks = []
for r in csv.DictReader(open(SRC, encoding="utf-8")):
    peaks.append({
        "ion": r["ion"],
        "number": int(r["number"]),
        "charge": int(r["charge"]),
        "mz": float(r["observed_mz"]),
        "theo": float(r["mz"]),
        "ppm": float(r["error_ppm"]),
        "intensity": float(r["intensity"]),
        "label": r["label"],
        "sequence": r["sequence"],
    })

c_max = max((p for p in peaks if p["ion"] == "c"), key=lambda p: p["number"])
z_max = max((p for p in peaks if p["ion"] == "z"), key=lambda p: p["number"])
c_pref, z_suf = strip(c_max["sequence"]), strip(z_max["sequence"])
assert all(c_pref[:p["number"]] == strip(p["sequence"]) for p in peaks if p["ion"] == "c")
assert all(z_suf[-p["number"]:] == strip(p["sequence"]) for p in peaks if p["ion"] == "z")

full_annot = c_max["sequence"] + z_suf[-1]           # c covers 1..101; last residue from z
residues = []
pos = 0
pending_mod = None
for tok in re.findall(r"\[[^\]]*\]|[A-Z]", full_annot):
    if tok.startswith("["):
        pending_mod = tok[1:-1]
        continue
    pos += 1
    residues.append({"pos": pos, "char": tok, "mod": pending_mod})
    pending_mod = None
L = pos
assert L == c_max["number"] + 1, (L, c_max["number"])

mods = [r for r in residues if r["mod"]]

# annotations fixture: base peak + top z flag, c-ladder core band (dynamic —
# survives data refreshes; the z flag picks the strongest ion OUTSIDE the band
# so the two top labels never collide)
base = max(peaks, key=lambda p: p["intensity"])
c_mzs = sorted(p["mz"] for p in peaks if p["ion"] == "c")
q = len(c_mzs) // 4
band_lo, band_hi = c_mzs[q], c_mzs[3 * q]
z_top = max((p for p in peaks if p["ion"] == "z"
             and not band_lo <= p["mz"] <= band_hi),
            key=lambda p: p["intensity"])

out = {
    "data": {
        "peaks": peaks,
        "residues": residues,
        "tolerance_ppm": 15,
        "scan": "subunit alpha · scan 2201 · ModB case",
        "annotations": [
            {"kind": "callout", "mz": base["mz"], "label": "base peak"},
            {"kind": "callout", "mz": z_top["mz"], "label": "top z-ion"},
            {"kind": "band", "from": c_mzs[q], "to": c_mzs[3 * q], "label": "c-ion ladder"},
        ],
    },
    "meta": {
        "source": "fictional placeholder example; full peptide reconstructed from c/z rules (asserted in gen-data.py)",
        "thresholds": f"fragment tolerance 15 ppm; L={L} residues; {len(mods)} modified sites; peaks plotted at observed m/z",
        "n": str(len(peaks)),
    },
}
(HERE / "data.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print(f"peaks: {len(peaks)} | L: {L} | mods: {[(r['pos'], r['char'], r['mod']) for r in mods]}")
