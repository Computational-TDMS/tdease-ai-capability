#!/usr/bin/env python3
"""Regenerates data.json for the candidate-discrimination example.

The committed example is FICTIONAL placeholder data: a synthetic spectrum,
two synthetic models (ModA vs ModB at a fictional site) and synthetic
matches, shaped exactly like the adapter output so the renderer develops
against the real schema. It intentionally carries no real acquisition,
sequence or chemistry. Run: python gen-data.py
"""
import json
import random
from pathlib import Path

random.seed(20261016)
HERE = Path(__file__).resolve().parent
BASE = "AQKELGDMSTVPNFQRWHYCI"
seq = (BASE * 6)[:102]
seq = seq[:46] + "K" + seq[47:]
SITE, MOD_A, MOD_B = 47, 186.10, 79.97
RES = {"A": 71.03, "Q": 128.06, "K": 128.10, "E": 129.04, "L": 113.08,
       "G": 57.02, "D": 115.03, "M": 131.04, "S": 87.03, "T": 101.05,
       "V": 99.07, "P": 97.05, "N": 114.04, "F": 147.07, "R": 156.10,
       "W": 186.08, "H": 137.06, "Y": 163.06, "C": 103.01, "I": 113.08}

peaks, mz = [], 149.0
while mz < 2000.0 and len(peaks) < 7309:
    mz += random.uniform(0.004, 0.42)
    i = random.choice([random.uniform(80, 900), random.uniform(900, 9000),
                       random.uniform(9000, 60000)])
    peaks.append({"mz": round(mz, 8), "i": round(i, 3)})


def ion_series(ion, model_delta):
    out, mass = [], 3.0 if ion == "c" else 0.0
    lo, hi = (3, 71) if ion == "c" else (9, 98)
    for n in range(1, 103):
        mass += RES[seq[n - 1]]
        if n < lo or n > hi or random.random() < 0.42:
            continue
        covered = (ion == "c" and n >= SITE) or (ion == "z" and n >= 102 - SITE + 1)
        extra = (86.04 + (model_delta if covered else 0.0)) if ion == "c" else (
            model_delta if (102 - n) < SITE else 0.0)
        for q in range(1, 15):
            if random.random() < (0.75 if q <= 3 else 0.28):
                theo = (mass + extra + (17.03 if ion == "c" else 1.0) + q * 1.007) / q
                if 140 < theo < 2000:
                    out.append({"ion": ion, "n": n, "q": q, "cleav": n,
                                "mz": round(theo * (1 + random.uniform(-4e-6, 4e-6)), 8),
                                "ppm": round(random.uniform(-2.4, 2.9), 3),
                                "i": round(random.uniform(600, 45000), 2)})
    return out


def candidate(cid, label, delta, color, gap):
    matches = ion_series("c", delta) + ion_series("z", delta)
    missing = ([{"start": gap, "end": gap + 9, "length": 10},
                {"start": 62, "end": 66, "length": 5}] if gap else [])
    counts = {i: len({m["n"] for m in matches if m["ion"] == i})
              for i in ("c", "z")}
    return {"id": cid, "label": label, "matches": matches, "counts": counts,
            "total": len(matches), "seriesMean": {},
            "seriesDivergence": {"flagged": False,
                                 "nMinusCPpm": round(random.uniform(-1.5, 1.5), 3),
                                 "thresholdPpm": 2.0},
            "missingRuns": {"c": missing, "z": []},
            "worstGap": {"c": missing[0] if missing else None, "z": None},
            "stateCounts": {},
            "coverage": round(len({m["n"] for m in matches}) / 101, 4),
            "mods": [{"pos": 1, "tag": "ac", "nTerm": True, "diff": False},
                     {"pos": 17, "tag": "me2", "nTerm": False, "diff": False},
                     {"pos": 71, "tag": "ox", "nTerm": False, "diff": False},
                     {"pos": SITE, "tag": cid, "nTerm": False, "diff": True}],
            "color": color}


doc = {"data": {
    "peaks": peaks,
    "candidates": [candidate("m1", f"model A (ModA@K{SITE})", MOD_A, "#2a9d8f", 47),
                   candidate("m2", f"model B (ModB@K{SITE})", MOD_B, "#e76f51", None)],
    "sequence": seq, "focusSite": SITE, "scan": "2201", "activation": "ETD",
    "tol": 10.0, "charges": list(range(1, 18)), "minRelIntensity": 0.01,
    "precursorCharge": 17},
    "meta": {"title": f"K{SITE} ModA-vs-ModB discrimination (fictional example data)",
             "source": "synthetic placeholder — no real acquisition",
             "thresholds": "10 ppm; Pearson >= 0.7; display >= 0.01 base peak",
             "n": "321", "generated": "2026-10-16T00:00:00Z", "tool": "chart-forge"}}
(HERE / "data.json").write_text(
    json.dumps(doc, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
print("regenerated fictional data.json")
