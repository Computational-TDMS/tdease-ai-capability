# Report interpretation

Read the recorded data, method and section status. Colors and projections
display those values; they do not diagnose biology.

| Section | What it shows | What to assess |
|---|---|---|
| completeness | Features detected per sample | Detection differences alongside input amount, acquisition and biological context |
| heatmap | Selected variable features on the displayed transform | Selection bias, missing cells and sample metadata |
| correlation | Pairwise expression correlations | Replicate/context consistency and outliers without a universal cutoff |
| featuremap | Raw TopFD feature RT intervals and m/z | Acquisition structure, recorded charge/intensity and shifts between files |
| levels | Proteoform ambiguity levels, weighted by spectral_count | Identity ambiguity, separately from intact/truncated sequence bounds |
| pca | Samples projected from transformed/imputed expression | Grouping relative to batch, missingness and preprocessing |
| rt-alignment | Applied aligned-minus-raw RT corrections | Correction size, not independently measured alignment accuracy |
| mass-error | Residuals on saved calibration rows | Calibration fit, not held-out mass accuracy |
| volcano | Group B/A log2FC and BH-adjusted significance | Experimental design, normalization, replicate support and detection filters |

Read thresholds and test definitions from the comparison artifact. Presence-only
features are outside the standard volcano test; their absence from it does not
mean no difference. PCA, correlation and a group-colored heatmap are related
summaries of the same matrix, not independent validation of group effects.

Use exact recorded values rather than estimating from color. A visible pattern
can motivate an explanation; distinguish technical and biological alternatives
with the relevant design or additional evidence. Spectral measurements address
identity/localization, not the validity of a group-level statistical test.

Saved calibration may be unavailable in older runs; report the skipped view.
Missingness must remain distinguishable from measured low abundance.
