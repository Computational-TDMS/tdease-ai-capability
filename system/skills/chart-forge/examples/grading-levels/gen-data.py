# Regenerates data.json for the grading-levels example through the REAL payload
# builder (figures_report.smith_histogram_payload) — the example is the
# adapter's own output, never hand-written.
#
# The synthetic 100-assignment record exercises the grading contract: 15
# mass-only assignments (level 5, all four ambiguities) move one level toward
# fully defined once the pipeline's direct attribution resolves the
# PTM-identity ambiguity; the other cohorts stay as searched.
#
#   python gen-data.py
#
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "tdp-opensearch" / "scripts"))

from figures_report import smith_histogram_payload


def observation(index: str, best: int, ambiguities: list[str]) -> dict:
    return {"assignmentId": index,
            "smith": {"best": best, "ambiguities": ambiguities},
            "topPic": {"carriers": []}}


observations = (
    [observation(f"A{i:03d}", 1, []) for i in range(30)]
    + [observation(f"A{i:03d}", 2, ["gene"]) for i in range(30, 55)]
    + [observation(f"A{i:03d}", 3, ["gene", "sequence"]) for i in range(55, 75)]
    + [observation(f"A{i:03d}", 4, ["gene", "sequence", "ptm_localization"]) for i in range(75, 85)]
    + [observation(f"A{i:03d}", 5, ["gene", "sequence", "ptm_identity", "ptm_localization"])
       for i in range(85, 100)]
)
annotations = [{"basis": "automatic approved UniMod mass attribution",
                "assignmentIds": [f"A{i:03d}" for i in range(85, 100)]}]

payload = smith_histogram_payload(observations, annotations)
payload["meta"]["generated"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
payload["meta"]["tool"] = "chart-forge"
(HERE / "data.json").write_text(
    json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print("data.json regenerated:", json.dumps(payload["data"]["levels"]))
