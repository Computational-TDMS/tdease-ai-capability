---
title: IsoForma
kind: general
status: accepted
updated: 2026-08-02
tags: [isoforma, positional-isomer, top-down, r-package, quantification]
type: entity
subjects: [data-analysis, quantification]
sources: [degnan-2024-isoforma.md]
generated_by: human
generated_at: 2026-08-02
verified: [by: human:zhou]
---

## Identity

IsoForma is an open-source R package for the relative quantification of
positional isomers (PI) in top-down LC-MS/MS data.

## Aliases

- IsoForma

## Properties

- Input: proteoform-specific fragment ions matched to top-down MS2.
- Output: relative abundance of each positional isomer.
- Method: matches proteoform-specific fragments to localize and quantify
  modifications, replacing manual annotation.
- Benchmarked against two existing workflows: comparable results, faster.

## Relationships

- Positional isomers: proteoforms with identical total mass and modification
  set but different PTM site combinations.
- Top-down MS/MS: fragment evidence localizes the modification site.

## Bounds

- Quantifies relative abundance, not absolute.
- Depends on fragment coverage localizing the sites; sparse fragmentation
  weakens the ratio estimate.
