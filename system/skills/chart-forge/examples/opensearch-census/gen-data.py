"""tdp-opensearch figures.py census reproduction — the mass-shift census view.

Real Chlamydomonas PNNL-batch census (3 TopPIC PrSM tables) with the evidence
marks that adjudicated its clusters. The matplotlib/plotly renderer's framing
contract is replicated exactly: evidence marks project onto the nearest
cluster within 1.5 Da (marks are overlays, not observations); the axis keeps
every shift observed >= floor = max(5, ceil(0.5% of top cluster count)) times
and force-keeps marked masses ±15 Da; the tail below the bar stays in data,
not on the axis."""
import json
import math
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
REAL = HERE / "real"

import gzip
read_json = lambda p: json.loads(gzip.decompress((REAL / p).read_bytes()))
census = read_json("sample-census-full.json.gz")
marks_doc = read_json("sample-census-marks.json.gz")
clusters = census["clusters"]
assert clusters, "census has no clusters"

shifts = [float(c["center_mass"]) for c in clusters]
marks_by_index = {}
for mark in marks_doc["marks"]:
    index = min(range(len(shifts)), key=lambda i: abs(shifts[i] - mark["mass"]))
    if abs(shifts[index] - mark["mass"]) > 1.5:
        continue
    marks_by_index.setdefault(index, []).append(mark)

floor = max(5, math.ceil(0.005 * max(int(c["count"]) for c in clusters)))
kept = [float(c["center_mass"]) for c in clusters if int(c["count"]) >= floor]
marked = [shifts[i] for i in marks_by_index]
lo = min(min(kept), min(marked) - 15.0) if marked else min(kept)
hi = max(max(kept), max(marked) + 15.0) if marked else max(kept)
pad = max(10.0, (hi - lo) * 0.03)
x_range = [lo - pad, hi + pad]

out_clusters = []
for i, c in enumerate(clusters):
    attached = marks_by_index.get(i, [])
    out_clusters.append({
        "mass": c["center_mass"], "count": c["count"], "samples": c.get("sample_count"),
        "localized": c.get("localized_count"), "proteins": c.get("distinct_proteins"),
        "dominant": c.get("dominant_protein"), "residues": c.get("candidate_residues"),
        "best_e": c.get("best_e_value"), "context": c.get("sequence_context"),
        "mark_labels": [m["label"] for m in attached],
        "mark_colors": [m.get("color") for m in attached if m.get("color")],
        "mark_verdicts": [m.get("verdict") for m in attached if m.get("verdict")],
        "inview": x_range[0] <= c["center_mass"] <= x_range[1],
    })

out_marks = []
for i, ms in sorted(marks_by_index.items()):
    for m in ms:
        out_marks.append({
            "mass": shifts[i], "label": m["label"], "role": m.get("role"),
            "verdict": m.get("verdict"), "color": m.get("color") or "grey",
            "basis": m.get("basis", ""), "observations": m.get("cluster_observations"),
        })

out = {
    "data": {
        "clusters": out_clusters,
        "marks": out_marks,
        "xlim": x_range,
        "floor": floor,
    },
    "meta": {
        "source": "example: census tables (fictional batch, 3 PrSM tables); projection + auto-xlim replicated from tdp-opensearch figures.py",
        "thresholds": f"{census['cluster_count']} clusters total, {sum(c['inview'] for c in out_clusters)} on axis (floor: >= {floor} observations); {len(out_marks)} marks projected at <= 1.5 Da; marks are overlays, not observations",
        "n": str(census["cluster_count"]),
    },
}
(HERE / "data.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print(f"clusters {len(out_clusters)} | inview {sum(c['inview'] for c in out_clusters)} | floor {floor} | marks {len(out_marks)} | xlim [{x_range[0]:.1f}, {x_range[1]:.1f}]")
