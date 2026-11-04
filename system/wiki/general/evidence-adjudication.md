---
title: "Evidence adjudication requires a separate activity"
kind: general
status: accepted
updated: 2026-08-02
tags: [evidence, provenance, adjudication, scientific-truth]
type: rule
subjects: [quality-control]
sources: [SCIENTIFIC-LAYER.md]
generated_by: human
generated_at: 2026-08-02
verified: [by: human:zhou]
---

## Rule

`proposed`, `challenged` and `inconclusive` are descriptions the Agent may
assert freely; `supported` and `refuted` are adjudications that may NOT be
minted in the same call that states the claim.

## Conditions

- Applies to any Evidence claim in the Workspace.
- The claim is being recorded or evaluated by an AI tool.

## Expected Outcome

`supported`/`refuted` require a separate adjudicating activity — a person or a
deterministic sufficiency check. The claim call itself stays at the descriptive
status.

## Exceptions

- A deterministic check that fully determines the claim may adjudicate.
- A human explicitly adjudicating may mint the status.

## Bounds

- This is a provenance truth: a claim and its adjudication are different acts
  and must not share a minting moment.
- No background, automatic adjudication.
