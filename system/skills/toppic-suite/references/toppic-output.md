# TopPIC output

## Files a TopPIC run produces

| File | Pattern | Content |
|---|---|---|
| PrSM table | `*_ms2_toppic_prsm.tsv` | one row per proteoform-spectrum match; the primary candidate source |
| PrSM table, filtered | `*_ms2_toppic_prsm_single.tsv` | best match per spectrum after filtering |
| Proteoform table | `*_ms2_toppic_proteoform.tsv` | matches grouped to one row per proteoform |
| MS1 features | `*_ms1.feature` | deconvoluted MS1 features: mass, charge range, RT range, intensity |
| MS2 features | `*_ms2.feature` | precursor feature per MS2 scan |
| Deconvoluted spectra | `*_ms2.msalign`, `*_ms1.msalign` | TopFD output consumed by TopPIC |
| HTML report | `*_html/` | per-PrSM visual report |

Read the actual filenames from the Run manifest. Do not construct them from a
sample name; TopFD and TopPIC prefixes differ between pipeline versions.

The PrSM table answers "which spectra matched"; the proteoform table answers
"which distinct proteoforms were identified". Counting proteoforms from the PrSM
table without grouping overstates the result — the same proteoform typically
matches several spectra across the elution peak.

## Preamble

Both tables begin with a parameter block delimited by `*` lines:

```text
********************************** Parameters **********************************
Protein database file:                      proteins.fasta
Fixed modifications:                        None
Maximum number of unexpected modifications: 2
Error tolerance:                            10 ppm
Spectrum-level FDR cutoff:                  0.01
********************************************************************************
Data file name	Prsm ID	Spectrum ID	...
```

This block is data worth reading: the error tolerance and the maximum number of
unexpected modifications bound every conclusion drawn from the table below it.
A search that allowed at most two unexpected shifts cannot report a third.

## Columns

| Meaning | Column |
|---|---|
| match identity | `Prsm ID`, `Spectrum ID`, `Proteoform ID`, `Feature ID` |
| acquisition | `Scan(s)`, `Retention time`, `Fragmentation`, `Charge`, `#peaks` |
| precursor | `Precursor mass`, `Adjusted precursor mass`, `Feature intensity`, `Feature apex time` |
| protein | `Protein accession`, `Protein description`, `#Protein hits` |
| chain | `First residue`, `Last residue`, `Protein N-terminal form`, `Database protein sequence` |
| proteoform | `Proteoform`, `Proteoform mass`, `Special amino acids` |
| modifications | `Fixed PTMs`, `#variable PTMs`, `variable PTMs`, `#unexpected modifications`, `unexpected modifications` |
| scores | `MIScore`, `#matched peaks`, `#matched fragment ions`, `E-value`, `Spectrum-level Q-value`, `Proteoform-level Q-value` |

Notes that matter:

- `Precursor mass` is the observed monoisotopic neutral mass; `Adjusted
  precursor mass` is TopPIC's corrected value after isotope-peak reassignment.
  A difference between them near 1.003 Da is a corrected monoisotopic call.
- `Proteoform mass` is theoretical for the reported proteoform. Its difference
  from the adjusted precursor mass is the residual the search did not explain.
- `Protein N-terminal form` records `NONE`, `NME` (initiator methionine removed),
  `NME_ACETYLATION` or `M_ACETYLATION`. It changes what "position 1" means.
- `Unexpected modifications` is a **residual**: TopPIC folds the assigned
  N-terminal form and the named fixed/variable PTMs in *before* reporting it.
  It is not the proteoform's total delta against the bare database sequence.
  Worked check on a real row: an `NME_ACETYLATION` proteoform of 11391.3531 Da
  over residues 2–103 (bare stretch 11319.3476) reports unexpected `+29.9947` —
  `11319.3476 + 42.0106 + 29.9947 = 11391.3529`, while its total delta against
  the unexcised bare sequence is −59.035 (−131.0405 NME + 42.0106 acetyl +
  29.9947 unexpected). In one representative run 56% of rows carried a
  non-NONE terminal form, so reading the shift as "PTM vs bare sequence" is
  wrong more often than not. Rebuild the complete frame (terminal form + named
  PTMs + unexpected) before composing ProForma candidates.
- `First residue` / `Last residue` are 1-based positions in the database
  sequence. A proteoform position is only meaningful together with them.

## Proteoform string grammar

```text
A.(GGKAK)[Acetyl]SPGKGGKAKSPGKG.V
│  └────┬────┘└──┬──┘                └ residue after the C-terminal boundary
│       │        └ modification applied to the preceding span
│       └ parenthesised span: the modification is somewhere inside it
└ residue before the N-terminal boundary
```

- The first and last dot-separated fields are the flanking residues in the
  database sequence, not part of the proteoform.
- `[Name]` is an identified modification. `[42.0106]` is a bare mass — the search
  found a shift it could not name. The distinction is exactly PTM-identity
  ambiguity.
- `(span)[mod]` means localization to a span. A single residue followed by
  `[mod]` means localization to that residue.
- A `[mod]` before any residue is N-terminal.
- Multiple modifications appear in sequence order and may mix named and bare-mass
  forms in one string.

Parse the string; do not regex the whole row. The scripts in this skill share one
parser that returns residues and modifications with explicit spans, and it is
short enough to copy if you need a variant.
