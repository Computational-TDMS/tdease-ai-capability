---
name: feature-table-leading-spaces
description: "TopFD feature/msalign TSV exports can carry leading spaces in column names ('  Feature_ID') — normalize headers before dict access or the read dies in KeyError."
version: 1
signals: ["Feature_ID", "KeyError", "msalign", ".feature", "feature table"]
scope: tdp-opensearch
source: session:harness_326653bd-2b36-4b10-ad75-04d367f0d1b9
created: 2026-09-04
---

⚠ A `csv.DictReader` over a TopFD `.feature` table died
`KeyError: ' Feature_ID'` while a sibling export used clean names
(`Feature_ID`) — the leading spaces are per-file, so a script that worked on
one sample fails on the next.

Fix: normalize headers (`{k: k.strip() for k in reader.fieldnames}`) or access
via stripped aliases before any row lookup.

Bounds: TopFD/Promex tabular exports in this workspace; verified on the
chlamy-pnnl-batch run.
