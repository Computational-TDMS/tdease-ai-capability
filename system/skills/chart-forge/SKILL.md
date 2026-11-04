---
name: chart-forge
description: Build source-backed scientific HTML figures; deliver the HTML/data pair in replies — PNG/SVG only on explicit request.
---

# chart-forge

The deliverable is a reproducible scientific figure.

## Scientific figure

Prefer the producing domain skill's adapter; it supplies scientific derivation,
source provenance and a matching renderer. Otherwise use INDEX.md to select a
shipped template and read its data contract.

    uv run --script "$TDEASE_SYSTEM_ASSETS/skills/chart-forge/scripts/build.py" \
      --example <template> --data <contract.json> --out <figure.html>

Input is { "data": ..., "meta": {...} }; CSV is accepted for flat rows.
Keep exact plotted data and source identity with the HTML and one manifest.
Include units, relevant transforms, thresholds and missingness. Derive chemistry
and statistics before rendering; the renderer displays the supplied values.

Use offline shipped assets: D3, or 3Dmol for molecular structures. Open the
result with an available browser/preview and check the main interaction and
readable labels. Report rendering failures or pending visual verification.
Export PNG/SVG when requested; canvas/WebGL use export.mjs, while pure SVG
figures have in-page controls.

If no shipped template fits, explain the missing view and continue the analysis
with inspectable data and interpretation. A new reusable figure is template
maintenance, not a prerequisite for answering the current question.

## Read when needed

- INDEX.md: shipped routes and matching contracts.
- references/engineering-contract.md: compiler, failure, export and template
  maintenance, including the audit.py --browser gate.
- references/scientific-integrity.md: domain-specific figure safeguards.
- paradigms/data-shaping.md: adapting an unfamiliar source.
- Other paradigms/: template design craft for maintenance, not routine analysis.

Product assets live in this installed skill. Developer galleries are not an
implicit runtime source.
