# /// script
# requires-python = ">=3.11"
# dependencies = ["pandas", "numpy"]
# ///
"""Recompute per-model discriminating-ion support from full kernel evidence dumps.

One transformation: evidence dumps in, support tables out. A compare card's
verdict lists at most a couple dozen ions; the -evidence.json dumps beside it
carry every fragment status of every model. This tool walks those full tables
and reports which ions actually separate the models — counting every accepted
match, not the card's truncated window — so a verdict can be checked against
its own complete evidence.

Usage:
  discrimination_table.py --evidence "spectra/cmp_*-evidence.json" [options]

  --evidence GLOB        kernel fragment-compare evidence dump (repeatable)
  --sample NAME          keep only these sample labels (repeatable; the label
                         is the file stem with -evidence/.json removed)
  --strip-prefix TEXT    remove this prefix from file stems before using them
                         as sample labels (e.g. --strip-prefix cmp_), default none
  --dmz DA               ion counts as differing between models only when its
                         predicted m/z spread across models exceeds this
                         (default 0.01)
  --out PATH             write the per-ion match CSV here
Stdout: per-model aggregate (matches, scans, samples, median ppm, median r),
the model x scan crosstab, and the most-supported ions per model.
"""
from __future__ import annotations

import argparse
import glob
import json
import os

import numpy as np
import pandas as pd


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--evidence", action="append", required=True)
    p.add_argument("--sample", action="append", default=[])
    p.add_argument("--strip-prefix", default="")
    p.add_argument("--dmz", type=float, default=0.01)
    p.add_argument("--out")
    a = p.parse_args()

    paths = []
    for pattern in a.evidence:
        paths.extend(sorted(glob.glob(pattern)))
    if not paths:
        raise SystemExit(f"no evidence files matched: {a.evidence}")

    rows: list[dict] = []
    scans = 0
    truncated = set()
    for path in paths:
        sample = os.path.basename(path).removesuffix("-evidence.json").removesuffix(".json")
        if a.strip_prefix and sample.startswith(a.strip_prefix):
            sample = sample.removeprefix(a.strip_prefix)
        if a.sample and sample not in a.sample:
            continue
        result = json.load(open(path, encoding="utf-8")).get("result", {})
        for comp in result.get("comparisons", []):
            scans += 1
            scan = str(comp.get("scan"))
            cands = comp["measurement"]["candidates"]
            status: dict[tuple[str, int, int], dict[str, dict]] = {}
            observed: dict[tuple[str, str, int, int], tuple] = {}
            for cand in cands:
                mid = cand["candidateId"]
                if cand.get("matchesTruncated"):
                    truncated.add(f"{sample}/{scan}/{mid}")
                for st in cand.get("fragmentStatuses", []):
                    status.setdefault((st["ion"], st["number"], st["charge"]), {})[mid] = st
                for m in cand.get("matches", []):
                    observed[(mid, m["ion"], m["number"], m["charge"])] = (
                        m.get("errorPpm"), m.get("envelopeCorrelation"), m.get("intensity"),
                    )
            for key, models in status.items():
                mzs = {mid: st["theoreticalMz"] for mid, st in models.items()}
                if max(mzs.values()) - min(mzs.values()) < a.dmz:
                    continue  # identical prediction in every model — separates nothing
                for mid, st in models.items():
                    if st["state"] != "matched":
                        continue
                    ppm, r, intensity = observed.get((mid,) + key, (None, None, None))
                    rows.append({
                        "sample": sample, "scan": scan,
                        "ion": f"{key[0]}{key[1]}:{key[2]}", "model": mid,
                        "ppm": ppm, "r": r, "intensity": intensity,
                        "n_models_matched": sum(1 for s in models.values() if s["state"] == "matched"),
                    })

    d = pd.DataFrame(rows)
    if d.empty:
        raise SystemExit(f"no discriminating matches across {scans} scans (dmz={a.dmz})")
    d["ppm"] = pd.to_numeric(d.ppm, errors="coerce").round(2)
    d["r"] = pd.to_numeric(d.r, errors="coerce").round(3)
    print(f"scans analysed: {scans} | discriminating accepted matches: {len(d)}")
    if truncated:
        print(f"WARNING matchesTruncated (counts are LOWER bounds) for: {sorted(truncated)}")
    print("\n=== accepted discriminating-ion matches per model ===")
    print(d.groupby("model").agg(
        matches=("ion", "size"), scans=("scan", "nunique"), samples=("sample", "nunique"),
        median_ppm=("ppm", "median"), median_r=("r", "median"),
    ).to_string())
    print("\n=== model x scan (count) ===")
    print(pd.crosstab([d["sample"], d["scan"]], d["model"]).to_string())
    print("\n=== most-supported ions per model ===")
    print(d.groupby(["model", "ion"]).size().sort_values(ascending=False).head(20).to_string())
    if a.out:
        d.to_csv(a.out, index=False)
        print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
