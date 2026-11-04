"""massql plot.py reproduction — MassQL SELECT traces on a shared RT axis.

Real pipeline output: three `tdease-ms select` statements over the real
a centroided sample mzML (fictional example run). TIC gives the run
context; the two XICs are the acetylation ladder charge states from
fictional-acetylation-ladder.json (0Ac z13, 1Ac z14). The matplotlib
renderer stacks all statements on one raw-intensity axis; TIC (~1e8) vs
XIC (~1e4) spans four decades there, so this D3 reproduction keeps the
shared-RT contract but splits magnitude domains into two stacked panels."""
import json
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
REAL = HERE / "real"

TRACES = [
    ("TIC · MS1", "SELECT spectra FROM mzml WHERE ms_level = 1 LIMIT 10000", "tic.json"),
    ("protein-A 0Ac · z13 (877.266)", "SELECT spectra FROM mzml WHERE ms_level = 1 AND product_mz WITHIN 877.266 PPM 10 LIMIT 10000", "xic-z13-unmod.json"),
    ("protein-A 1Ac · z14 (817.677)", "SELECT spectra FROM mzml WHERE ms_level = 1 AND product_mz WITHIN 817.677 PPM 10 LIMIT 10000", "xic-z14-1ac.json"),
]

out_traces = []
for label, statement, fname in TRACES:
    doc = json.loads((REAL / fname).read_text(encoding="utf-8"))["result"]
    rows = sorted(doc["spectra"], key=lambda r: r["retentionTime"])
    has_matched = any(r.get("matchedProducts") for r in rows)
    kind = "xic" if has_matched else "tic"
    if has_matched:
        ys = [float(sum(p["intensity"] for p in r["matchedProducts"])) for r in rows]
    else:
        ys = [float(r["totalIonCurrent"]) for r in rows]
    xs = [float(r["retentionTime"]) for r in rows]
    assert xs == sorted(xs) and len(xs) > 0, f"{fname}: rt must be ascending and non-empty"
    assert not doc["truncated"], f"{fname}: truncated at LIMIT — not a complete trace"
    out_traces.append({
        "label": label,
        "statement": statement,
        "kind": kind,
        "rt": [round(v, 3) for v in xs],
        "y": [round(v, 2) for v in ys],
    })

tic = out_traces[0]
assert tic["kind"] == "tic" and out_traces[1]["kind"] == "xic" and out_traces[2]["kind"] == "xic"

out = {
    "data": {
        "traces": out_traces,
        "mzml": "centroided.mzML (fictional-run)",
    },
    "meta": {
        "source": "example: tdease-ms select over sample-run/centroided.mzML; raw JSON kept in real/",
        "thresholds": f"XIC windows ±10 ppm; XIC m/z anchors from the acetylation ladder (0Ac 11391.369 Da z13, 1Ac 11433.379 Da z14); RT {out_traces[0]['rt'][0]:.1f}–{out_traces[-1] if False else out_traces[0]['rt'][-1]:.1f} min",
        "n": " · ".join(f"{t['label'].split(' ·')[0]}: {len(t['rt'])} pts" for t in out_traces),
    },
}
(HERE / "data.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print(" | ".join(f"{t['label']}: {len(t['rt'])}" for t in out_traces))
