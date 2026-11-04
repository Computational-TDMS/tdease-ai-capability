# TDP Copilot action index

Choose the capability that addresses the current unknown. Investigation order
and completion belong to ../SKILL.md.

| Capability | Question | Result |
|---|---|---|
| tdp-opensearch inspect | What did the source record actually report? | IDs, backbone, mass basis, quality, proposal inputs |
| refine-shift | Is the source mass/isotope frame reliable enough to constrain candidates? | Provisional MS1 frame and QC |
| ptm-database | What annotations, motifs or homologs suggest explanations? | Source-attributed priors |
| UniMod search/solve | Which named modifications or combinations fit this mass frame? | Composition proposals |
| resolve-shift | Where does the shift sit, and can components be separated? | Intervals, component proposals, residual |
| Family analysis | Do related source-linked features constrain the explanation? | Increments, recurrence, member differences |
| spectral-evidence `measure --ms1` | Do complete compositions predict different envelopes? | Shared-observation composition comparison |
| spectral-evidence `test` | Which fragments separate the complete models or sites? | Evidence card: separating ions, contested peaks, kernel state per model |
| spectral-evidence `show` | What do the separating windows look like? | Per-model isotope-window panels beside the card's verdict |
| Characteristic-loss/diagnostic-ion check | Does a candidate predict a distinctive loss or ion — the discriminator when candidates are mass-degenerate at precursor mass? | Targeted observation (spectral-evidence tool-recipes) |
| parse-proforma | What sequence, mass/composition and coordinates does this representation encode? | Parsed measurement input, not proof of chemical feasibility or sample identity |

Use the producing skill's reader for its native source format. When additional
postprocessing is needed, use named columns and explicit joins, recording the
source and transformation. Read tool-recipes.md for unfamiliar command syntax,
or the relevant method reference for a scientific uncertainty; there is no
required tour of all references.
