---
name: structure-evidence
description: Map assigned or hypothetical sites to structures and measure solvent accessibility with coordinate provenance.
---

# Structural context

Describe the environment of measured, assigned or hypothetical sites, labeling
their origin. Structural context can guide a hypothesis without prior spectral
localization; it does not establish modification identity.

Choose a structure matching the question: --accession fetches AlphaFold;
--pdb selects an experimental RCSB structure. Check accession/isoform, sequence
coverage and numbering before mapping processed-chain sites. Experimental
author numbering, gaps and local offsets may differ from UniProt coordinates.

    uv run --script "$TDEASE_SYSTEM_ASSETS/skills/structure-evidence/scripts/structure_figure.py" \
      --accession <accession> --sites <start[-end][:label]> \
      --out <analysis-dir>/structure.html

Repeat --sites as needed, or replace --accession with --pdb <id>.
Structures cache under reference/structures/ with provenance; --refresh
downloads again.

Stdout carries source mapping, per-residue SASA and model-confidence context.
The script uses Tien theoretical maxima and descriptive relative-SASA bands.
Report measured values and method, not the band alone. Exposure does not prove
modification, burial does not rule it out, and low-confidence geometry limits
the structural inference.

Keep the result JSON with the analysis and present the figure when useful:
![](workspace-relative/path/structure.html).
