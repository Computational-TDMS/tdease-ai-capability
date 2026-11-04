# Suite And Files

## Stage map

```text
vendor RAW
  -> MSConvert
  -> sample.mzML
  -> TopFD
  -> sample_ms1.feature
  -> sample_ms2.feature
  -> sample_ms1.msalign
  -> sample_ms2.msalign
  -> TopPIC + database.fasta
  -> sample_ms2_toppic_prsm.tsv
  -> sample_ms2_toppic_proteoform.tsv
```

TopIndex prepares a reusable native FASTA index. It is not an unconditional stage:
use it only when the selected TopPIC version/workflow requires or benefits from the
native index and the Backend storage budget can support it.

TopMG is a different identification strategy based on mass graphs. Do not substitute
it for TopPIC merely because both consume deconvoluted spectra.

## Flat Run layout

```text
runs/<run>/
  Snakefile
  config.yaml
  sample.raw
  database.fasta
  sample.mzML
  sample_ms1.feature
  sample_ms2.feature
  sample_ms1.msalign
  sample_ms2.msalign
  sample_ms2_toppic_prsm.tsv
  sample_ms2_toppic_proteoform.tsv
```

TopFD and TopPIC write related files next to their inputs. Keep native names unless a
specific Node rule declares otherwise. Do not add staging/result subdirectories only
to make the workflow look organized.

## Selection questions

- Is the input vendor RAW or already valid centroided mzML?
- Does TopFD have MS1 data and the matching feature-file expectations?
- Is the FASTA target-only or target-decoy?
- Is the search closed, variable-PTM, or unrestricted mass-shift/open search?
- Are results required per sample only, or will downstream match-between-runs be used?
- Does the Run declare the actual TSV outputs the user expects to retrieve?
