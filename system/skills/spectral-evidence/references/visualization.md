# Available spectral views

Every figure is reached through `spectral.py`; there is no second renderer to
find.

| Observation | Route | View |
|---|---|---|
| Raw MS2 audit / comparison | `measure --mzml --scan --models` | Spectrum, candidate error maps, sequence strips, ppm histograms |
| MS2 fragment isotope windows | `show <card.json>` | Per-model observed/theoretical mirror panels |
| MS1 composition (feature family) | `measure --ms1 --feature ... --models` | Member envelope comparison with per-model fit statistics |
| MS1 mass refinement | `show --ms1 --refine <refine-shift.json>` | Source envelope, isotope-offset evaluations, fit metrics |
| MS1 composition (existing artifacts) | `show --ms1 --observed --theory` | Shared observation with model envelopes |
| Deconvoluted MS2 annotation | `annotate --msalign` | Neutral-mass spectrum and ppm panel |

Choose by the actual observation domain. MS1 refinement does not provide MS2
fragment errors; deconvoluted masses do not retain raw isotope clusters.

Raw MS2 views include all standard panels. `--focus-site` controls the local
summary; the default focus is a display choice, not a site assignment.
Series-divergence warnings and measurement metrics remain useful alongside the
plot; a colored map is not a substitute for their values.

Each mirror panel's subtitle carries the kernel's verdict *and* the window's own
slot count. They answer different questions — a populated window can still be
rejected on the full-slot Pearson — so quote both when the figure is the
evidence. Correlation verdicts are never inferred from the window.

The msalign view accepts repeatable `--annotate 'mz=LABEL'` or
`'mz1..mz2=LABEL'`. These are attention pointers, not measured data.

Adapters use the shared offline rendering templates with data contracts and
manifests.
Reuse recorded data to redraw: every figure leaves its `.data.json` contract,
its extraction sidecars and a manifest beside the HTML. Select the views that
make the conclusion inspectable; if a view fails or is unavailable, identify
that limit and provide the underlying measurement.
