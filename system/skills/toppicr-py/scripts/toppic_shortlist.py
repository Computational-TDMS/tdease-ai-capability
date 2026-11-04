# /// script
# requires-python = ">=3.11"
# dependencies = ["pandas", "numpy"]
# ///
"""Tier a comparison result into a detection-weighted shortlist.

One transformation: the comparison JSON plus an occupancy table in, tiered
shortlist and gene concordance out. A significant fold change resting on two
detections in the weaker group is a different kind of finding from one seen in
every sample, and a flat volcano plot hides the difference. This tool splits
the significant set by per-group detection counts, orders the strongest tier
by effect size, and reports genes whose significant features move concordantly
— the inputs a "worth deeper investigation" judgement actually needs.

Usage:
  toppic_shortlist.py --comparison comparison_full.json --occupancy features_with_occupancy.csv
                      --groups WT,M [options]

  --comparison PATH   toppicr-py comparison output (candidates list)
  --occupancy PATH    CSV indexed by featureId with per-group detection
                      counts in columns named n_<group>
  --groups A,B        group names, matching the n_<group> columns and the
                      log2FoldChange direction convention of the comparison
  --fdr Q             significance cutoff, default 0.05
  --abs-log2fc F      effect cutoff, default 1
  --tier-a-det N      detections required in BOTH groups for tier A,
                      default 4 (tier B = 3, tier C = the rest)
  --out-shortlist PATH / --out-genes PATH   write the full tables here
Stdout: tier counts, the tier-A head table, concordant genes.
"""
from __future__ import annotations

import argparse
import json

import pandas as pd


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--comparison", required=True)
    p.add_argument("--occupancy", required=True)
    p.add_argument("--groups", required=True, help="A,B e.g. WT,M")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--abs-log2fc", type=float, default=1.0)
    p.add_argument("--tier-a-det", type=int, default=4)
    p.add_argument("--out-shortlist")
    p.add_argument("--out-genes")
    a = p.parse_args()

    ga, gb = (s.strip() for s in a.groups.split(",", 1))
    n_a, n_b = f"n_{ga}", f"n_{gb}"
    full = pd.DataFrame(json.load(open(a.comparison, encoding="utf-8"))["candidates"])
    occ = pd.read_csv(a.occupancy, index_col=0)
    occ = occ[~occ.index.duplicated()]
    need = [c for c in (n_a, n_b) if c not in occ.columns]
    if need:
        raise SystemExit(f"occupancy table lacks {need}; pass --groups matching its n_<group> columns")
    f = full.set_index("featureId").join(occ[[n_a, n_b]])
    f["min_det"] = f[[n_a, n_b]].min(axis=1)
    sig = (f["bhFdr"] < a.fdr) & (f["log2FoldChange"].abs() > a.abs_log2fc)
    f["significant"] = sig
    tier_a = sig & (f["min_det"] >= a.tier_a_det)
    tier_b = sig & (f["min_det"] == a.tier_a_det - 1)
    tier_c = sig & (f["min_det"] < a.tier_a_det - 1)
    print(f"tested {len(f)} | significant {int(sig.sum())} "
          f"(higher in {gb}: {int((sig & (f.log2FoldChange > 0)).sum())}, "
          f"higher in {ga}: {int((sig & (f.log2FoldChange < 0)).sum())})")
    print(f"tier A (>={a.tier_a_det} det both): {int(tier_a.sum())} | "
          f"tier B ({a.tier_a_det - 1}): {int(tier_b.sum())} | "
          f"tier C (<={a.tier_a_det - 2}): {int(tier_c.sum())}")

    cols = [c for c in ("gene", "log2FoldChange", "bhFdr", n_a, n_b,
                        "spectralCount", "proteoform_level", "unexpectedModificationCount")
            if c in f.columns]
    head = f[tier_a].reindex(f[tier_a]["log2FoldChange"].abs().sort_values(ascending=False).index)
    print("\n=== tier A, |log2FC| sorted (head) ===")
    print(head[cols].head(30).to_string())

    g = (
        f[f["significant"]].groupby("gene")
        .agg(n=("log2FoldChange", "size"),
             up=("log2FoldChange", lambda x: int((x > 0).sum())),
             dn=("log2FoldChange", lambda x: int((x < 0).sum())),
             medFC=("log2FoldChange", "median"))
        if "gene" in f.columns else None
    )
    if g is not None:
        g["min_side"] = g[["up", "dn"]].min(axis=1)
        g = g.sort_values(["n", "min_side"], ascending=False)
        print("\n=== genes with >=2 significant features of concordant direction ===")
        print(g[g["min_side"] >= 2].head(20).to_string())
        if a.out_genes:
            g.to_csv(a.out_genes)

    if a.out_shortlist:
        f.to_csv(a.out_shortlist)
        print(f"\nwrote {a.out_shortlist}")
    if a.out_genes and g is not None and not a.out_shortlist:
        pass


if __name__ == "__main__":
    main()
