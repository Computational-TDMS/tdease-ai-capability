# Scores and identity ambiguity

| Output | Meaning | Limit |
|---|---|---|
| E-value | Expected chance matches at least this good under the search model | Not a probability that this PTM/site is correct |
| Spectrum Q-value | Estimated set-level false discovery rate at the match threshold | Not an individual posterior probability |
| Proteoform Q-value | Error-rate estimate at the proteoform reporting level | Does not establish localization |
| MIScore | Localization confidence conditional on the proposed mass shift | Does not establish its chemistry |

A small E-value supports the match under its search model, not every chemical
interpretation attached to it. A missing MIScore is unavailable evidence,
not perfect localization. Read thresholds and error-control scope from the
source/version.

Compare adjusted precursor mass with the complete reported model mass in the
same frame. Near-isotope-spacing residuals motivate checking registration;
near-known chemical deltas motivate alternatives. Neither pattern identifies
the cause without further evidence.

## Definition levels

The five-level scheme counts ambiguity in gene, sequence, PTM identity and
PTM localization: level = 1 + ambiguous category count. Definition completeness
and statistical confidence are distinct.

TopPIC rows can expose multiple protein hits, bare-mass annotations and
localization spans. Sequence ambiguity is not fully derivable from one row;
the shared classifier accepts an explicit sequence-ambiguity declaration.
State that assumption rather than treating an inferred level as independent
evidence.

Sources: [five-level scheme](https://pmc.ncbi.nlm.nih.gov/articles/PMC6857706/)
and [PrSM statistical significance](https://pmc.ncbi.nlm.nih.gov/articles/PMC6484843/).
