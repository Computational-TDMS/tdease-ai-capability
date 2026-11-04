#!/usr/bin/env python3
"""Generate feature-map/data.json from a PNNL .feature table (TSV).

All MS1 features pass through compacted to {mz, rt0, rt1, apex, logI, q} for
the canvas hybrid renderer — no filtering, no thresholding: the cloud IS the
figure. Default --src is the real sample 0Percoll run."""
import argparse
import csv
import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[6]   # .../TDEase-zero
DEFAULT_SRC = REPO / "data/workspace/runs/sample-pnnl-batch/20201107_0Percoll_2_5uL_PNNL_column_notrap_ms1.feature"


def main() -> None:
    here = Path(__file__).parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--src", type=Path, default=DEFAULT_SRC,
                    help="PNNL .feature table (TSV with Rep_average_mz, Min/Max_time, Apex_*)")
    ap.add_argument("--out", type=Path, default=here / "data.json")
    args = ap.parse_args()

    with args.src.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    feats = [{
        "mz": round(float(r["Rep_average_mz"]), 2),
        "rt0": round(float(r["Min_time"]), 2),
        "rt1": round(float(r["Max_time"]), 2),
        "apex": round(float(r["Apex_time"]), 2),
        "logI": round(math.log10(float(r["Apex_intensity"])), 2),
        "q": int(float(r["Rep_charge"])),
    } for r in rows]

    payload = {
        "data": feats,
        "meta": {
            "title": "Chlamy MS1 feature cloud: the 60–110 min gradient carries the strong signal",
            "subtitle": (f"{len(feats):,} MS1 features (all shown, canvas) · each bar spans its elution window "
                         "[Min_time, Max_time] at its average m/z · colour = apex intensity "
                         "(batlow reversed, darker = stronger) · charge + apex on hover"),
            "source": str(args.src.relative_to(REPO)) if str(args.src).startswith(str(REPO)) else str(args.src),
            "transform": "Min_time/Max_time (elution window), Rep_average_mz, log10(Apex_intensity), Rep_charge; unmodified pass-through",
            "colour_scale": "crameri-batlow reversed (sequential, intensity)",
            "n": len(feats),
        },
    }
    args.out.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    print(f"feature map data: {len(feats)} features, {args.out.stat().st_size / 1024:.0f} KB -> {args.out}")


if __name__ == "__main__":
    main()
