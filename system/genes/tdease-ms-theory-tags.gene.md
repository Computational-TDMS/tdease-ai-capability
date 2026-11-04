---
name: tdease-ms-theory-tags
description: "tdease-ms theory rejects numeric mass-only ProForma tags — declare a named modification tag, or a --delta family when the chemistry is deliberately unclaimed."
version: 1
signals: ["tdease-ms theory", "numeric tag", "--delta", "proforma", "theory"]
scope: tdease-ms
source: session:harness_326653bd-2b36-4b10-ad75-04d367f0d1b9
created: 2026-09-04
---

⚠ `tdease-ms theory` asserts an isotope shape, which a numeric mass-only tag
cannot define — the call fails instead of guessing.

Fix: use a named modification tag when the chemistry is claimed; reserve
`--delta <Da>` for a mass-delta family whose exact chemistry is deliberately
unclaimed. They are different hypotheses — do not interchange them to make a
call pass.

Bounds: the tdease-ms theory engine; quantifier-side tag handling may differ.
