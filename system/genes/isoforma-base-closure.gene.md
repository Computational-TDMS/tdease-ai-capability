---
name: isoforma-base-closure
description: "BASE_NOT_CLOSED on an isoforma run means an undeclared sister series rides inside the family window — declare its complete ProForma base and rerun; never tune the quantifier around the gate."
version: 1
signals: ["BASE_NOT_CLOSED", "proformaBaseClosure", "sister series"]
scope: isoforma
source: session:s_b7772756-fbb5-4765-8ca5-d4668e846f52
created: 2026-09-04
---

⚠ A family window whose declared base leaves feature masses unexplained dies
`BASE_NOT_CLOSED` before fitting. The gate is a diagnosis, not an obstacle:
beyond a few registration neutrons the offset is undeclared chemistry (e.g. a
methyl/oxidation sister series) riding inside the window.

Fix: read `proformaBaseClosure.unexplainedClusters` in the selection manifest,
declare each offset as its own complete ProForma base, and rerun the
selection. Sister series are reported as separate series — never merged into
one fit.

Bounds: isoforma-family runs only. A true isotope-registration offset is
closed differently, and a small PTM delta can masquerade as deep registration.
