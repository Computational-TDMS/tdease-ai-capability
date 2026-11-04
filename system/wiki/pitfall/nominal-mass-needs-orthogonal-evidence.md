---
title: "Nominal mass shifts require orthogonal evidence"
kind: pitfall
status: accepted
updated: 2026-08-02
tags: [unknown-ptm, top-down, mass-shift, open-search, unknown-modification, top-down]
type: pitfall
subjects: [identification]
sources: [wilson-zhou-2022.md]
generated_by: human
generated_at: 2026-08-02
verified: [by: human:zhou]
---

## Literal Symptom

In an open modification search, proteoforms can carry a nominal mass shift for
which several chemically distinct candidates or combinations are plausible.
The search engine alone cannot decide identity or site.

## Likely Causes

- Isobaric or near-isobaric database candidates.
- Combinations of smaller modifications with the same net shift.

## Discriminating Checks

- Calibrate the shift before narrowing tolerance.
- Compare exact composition and site feasibility for all candidates.
- Require fragment ions that separate candidate sites or chemistries.
- Check neutral losses only when the acquisition mode can expose them.

## False Friends

- Treating a nominal-mass hit as an identification assigns the wrong chemistry
  and can mislocalize the site.
- Wide search tolerances conflate distinct candidates into one hit.

## Corrective Direction

Disambiguate by exact mass first, then confirm with fragment ions and any
diagnostic neutral loss. Do not mint a `supported` claim from the open-search
hit alone — that requires the separate adjudication activity.

## Bounds

- Applies to top-down histone work with open modification search.
- Requires high-resolution MS (exact mass separates the candidates) and
  fragment-level confirmation.
- A low-resolution instrument or a wide search tolerance invalidates the
  discriminating checks.
