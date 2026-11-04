# /// script
# requires-python = ">=3.11"
# dependencies = ["pandas", "numpy"]
# ///
"""Aggregate exported PrSM rows into a per-species mass-level inventory.

One transformation, no chemistry: rows in, grouped levels out. Feed it the
PrSM rows of one protein family (toppic-suite prsm_read.py exports them, or
any CSV superset of the TopPIC column names) and it reports, per accession,
the distinct adjusted precursor-mass levels with their support — row counts,
samples, distinct unexpected-shift masses, best q-value — plus the mass gaps
between neighbouring levels. The gaps are the hint: a family of levels spaced
by known increments is one species series, not many findings.

Usage:
  prsm_inventory.py --rows family_prsm_rows.csv [options]

  --rows PATH              CSV of PrSM rows (TopPIC column names by default)
  --accession TEXT         filter: substring or regex on the accession column
                           (repeatable; default: all rows)
  --accession-col NAME     default "Protein accession"
  --mass-col NAME          default "Adjusted precursor mass"
  --shift-col NAME         default "Unexpected modifications" ("<mass>:[range]")
  --q-col NAME             default "Spectrum-level Q-value"
  --sample-col NAME        default "sample" (add it when exporting)
  --bin DA                 mass-level bin width: levels are multiples of this
                           (default 0.02)
  --out PATH               write the full inventory CSV here
Stdout carries the inventory tables themselves; nothing is file-only.
"""
from __future__ import annotations

import argparse
import re

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rows", required=True)
    p.add_argument("--accession", action="append", default=[])
    p.add_argument("--accession-col", default="Protein accession")
    p.add_argument("--mass-col", default="Adjusted precursor mass")
    p.add_argument("--shift-col", default="Unexpected modifications")
    p.add_argument("--q-col", default="Spectrum-level Q-value")
    p.add_argument("--sample-col", default="sample")
    p.add_argument("--bin", type=float, default=0.02)
    p.add_argument("--out")
    return p.parse_args()


def main() -> None:
    a = parse_args()
    d = pd.read_csv(a.rows, dtype=str)
    missing = [c for c in (a.accession_col, a.mass_col) if c not in d.columns]
    if missing:
        raise SystemExit(f"missing columns: {missing}; override with --accession-col/--mass-col")

    acc = d[a.accession_col].fillna("").astype(str)
    # UniProt-style "tr|ACC|name" keeps the accession; bare accessions pass through.
    d = d.assign(acc=acc.str.extract(r"(?:tr\|(\w+)\|)|(^\w+)", expand=True).bfill(axis=1)[0])
    if a.accession:
        mask = np.zeros(len(d), dtype=bool)
        for pat in a.accession:
            mask |= d["acc"].str.contains(pat, regex=True, na=False)
        d = d[mask]
    if d.empty:
        raise SystemExit("no rows matched")

    d["mass"] = pd.to_numeric(d[a.mass_col], errors="coerce")
    d = d.dropna(subset=["mass"])
    d["shift"] = pd.to_numeric(
        d.get(a.shift_col, pd.Series(dtype=str)).astype(str).str.extract(r"^([+-]?[\d.]+):")[0], errors="coerce"
    )
    if a.q_col in d.columns:
        d["q"] = pd.to_numeric(d[a.q_col], errors="coerce")
    else:
        d["q"] = np.nan
    d["level"] = (d["mass"] / a.bin).round() * a.bin

    g = (
        d.groupby(["acc", "level"])
        .agg(
            n=("mass", "size"),
            samples=(a.sample_col, lambda x: ",".join(sorted(set(x.dropna().astype(str))))),
            shifts=("shift", lambda x: sorted(set(np.round(x.dropna(), 2)))),
            mass_min=("mass", "min"),
            mass_max=("mass", "max"),
            q_min=("q", "min"),
        )
        .reset_index()
        .sort_values(["acc", "level"])
    )
    print(f"rows {len(d)} | accessions {d.acc.nunique()} | mass levels {g.groupby('acc').size().to_dict()}")
    print(g.to_string(index=False, max_colwidth=48))
    print("\n--- gaps between neighbouring levels (Da) ---")
    for acc_name, s in g.groupby("acc"):
        gaps = np.round(np.diff(s.level.values), 2)
        print(acc_name, "levels:", list(s.level.values), "gaps:", list(gaps))

    if a.out:
        g.to_csv(a.out, index=False)
        print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
