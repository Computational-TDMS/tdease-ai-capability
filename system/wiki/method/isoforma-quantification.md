---
title: "IsoForma positional isomer quantification workflow"
kind: project
status: proposed
updated: 2026-08-02
tags: [isoforma, positional-isomer, quantification, top-down, ms2]
type: method
subjects: [data-analysis, quantification]
sources: [degnan-2024-isoforma.md]
generated_by: agent
generated_at: 2026-08-02
---

## Objective

Estimate the relative abundance of protein positional isomers from top-down
LC-MS/MS fragment data, replacing manual annotation of proteoform-specific ions.

## Method Structure

1. Identify the set of positional isomers (same total mass, same modification
   set, different site combinations).
2. Match proteoform-specific fragment ions to the MS2 spectra.
3. Derive per-isomer relative abundance from the matched fragments.

## Quality Criteria

- Fragment ions must be proteoform-specific (discriminate the isomers).
- Sites with no discriminating fragment are not localizable by this method.

## Bounds

- Requires top-down MS2 with fragment coverage at the modification sites.
- Applies to relative quantification, not absolute abundance.
- Sparse or ambiguous fragment evidence degrades the ratio estimate.

## Alternatives

- Manual annotation in LcMsSpectator (Informed-Proteomics) for targeted
  confirmation; slower but human-verifiable.
