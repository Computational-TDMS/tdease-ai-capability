# Result semantics

Resolve output paths and formats from run_manifest.json. The manifest records
inputs, parameters, stage_metrics, warnings and saved calibration; the report
does not recreate missing measurements.

| Output | Meaning | Interpretation |
|---|---|---|
| expression | Feature-by-sample intensity matrix | Zero is not detected/recovered; it is not proof of absence |
| features | Representative identity and cluster metadata | Feature identity is a pipeline estimate, with any source ambiguity retained |
| samples | Dataset identifiers and supplied metadata | Verify biological groups/replicates instead of inferring them from names |
| calibration | Saved alignment/recalibration observations | Corrections and fit residuals describe the applied fit |

Feature IDs such as gene_pcGroup identify clustered proteoform features, not
individual PrSMs or fully resolved molecular identities. The matrix covers
features identified at least once; it is not a census of the entire proteome.

Use the actual table schema. Common feature fields include gene, pc_group,
uniprot_acc, first_aa, last_aa, prot_length, mass, rt, proteoform, collision,
spectral_count and n_unexpected_mods. Their different roles matter:

- Sequence bounds describe the representative mapped interval, not fragment
  coverage. Check processing and coordinate conventions before calling it intact.
- Collision flags indicate ambiguity under mass/RT clustering.
- Spectral counts describe repeated identifications, not an identity confidence
  probability or independent biological replication.
- Modification labels and unexpected shifts retain the search/annotation
  uncertainty; counting them does not establish chemical burden.

Explain stage changes from filters, inference and merging actually recorded.
A row-count decrease may be expected processing, not improved quality by itself.
A completed build with warnings or skipped report sections needs those limits
stated. No paper's example counts define successful contraction.

## Prioritization semantics

Two selection axes are first-class and must survive into any recommendation:

- **Effect size.** The extremes of the fold-change distribution — the largest
  between-group separations in the whole comparison — are candidates on
  their own numbers, whatever their tier. Detection counts are a caveat to
  state next to the feature (n of N per group), never a reason to drop it
  from the answer: on/off-switch biology lives exactly at low-occupancy
  extremes. Surface the top |log2FC| features explicitly, with their
  detection support, before or alongside any tier-filtered view.
- **Unnamed shift families.** A numeric-bracket shift in a proteoform string
  (e.g. `(AAAA)[+123.4567]AG…`) is a chemistry-unresolved observation, not a
  failed identification and not a named modification. Group recurring shift
  values across different proteins into families and surface the family: a
  shift that recurs on unrelated proteins is either a shared modification or
  a systematic search artifact — both are investigable questions, and neither
  is settled by which search was or was not run. Do not dismiss such features
  for lacking a specialized search, and do not name the chemistry from file
  names or search configuration.

Features that combine all three signals — an unexplained shift, an extreme
between-group separation, and recurrence across features or proteins — form
the top-priority class of the comparison. Name them as a class in the
recommendation rather than distributing them across unrelated findings: the
conjunction is the finding.

