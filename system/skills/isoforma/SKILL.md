---
name: isoforma
description: Estimate positional-isomer fractions from centroided mzML, declared configurations and fragment responses.
---

# IsoForma

Estimate fractions conditional on a processed sequence, modification model,
site set and acquisition window. Inputs include centroided mzML, fixed tags,
occupancies, RT bounds and activation. Use source-assigned or human-declared
modifications for sample quantification; a candidate fit is an exploratory
result until its chemistry is established independently.

## Run

    uv run --script "$TDEASE_SYSTEM_ASSETS/skills/isoforma/scripts/isoforma_run.py" \
      --mzml-dir <run> --sequence <processed-sequence> \
      --sites <residue,positions> --mod-tag <modification> \
      --fixed-tags <fixed-tags> --activation <activation> \
      --occupancies <occupancies> --rt-min <min> --rt-max <max> \
      > <analysis-dir>/isoforma.json

Create the output directory and replace placeholders from the source/model.
The runner owns configuration enumeration, scan selection, isotope gates,
centroid merging, fragment matching and fitting. Read --help for options;
inspect recorded parameters and selection warnings rather than repeating
these operations manually.

--feature-scan-map supplies source/occupancy-linked MS2 scans. Its default
isotope gate is diagnostic; without the map it is strict. --isotope-gate strict
requests strict rejection on either route. Choose from acquisition evidence,
not to recover an expected fraction.

For declared reporting categories use --isomer-model <csv>, with columns
occupancy,sites,category. Repeated category labels request one total. Categories
must be identifiable from the response matrix; merging is a model decision,
not evidence that the underlying configurations are identical.

## Interpret

Stdout is one JSON result. The kernel's full fragment-match evidence is kept
beside the mzML sources (or `--evidence-dir`) — one JSON per source and
occupancy; result.evidence names the file behind each fit, so a proportion
can be traced to its supporting ions without re-measuring. Read
acceptedScanCount, scanCountQc, modelSpace, selection QC, measured (the ion
series and activation actually used), calibration and estimator status.

- no_qualified_spectra: this route has no accepted observations.
- suggested_aggregates: candidate inseparable groups, not reported fractions.
- provisional: least-squares diagnostic result, not the formal beta estimator.
- calibration.status=unassessed: no source-specific calibration claim.

A solved fit remains conditional on its inputs and fragment-response model.
Assess identifiability, interference and calibration before quantitative claims.
Use symmetric response assumptions across configurations; do not change scans,
tolerances or the model to reproduce literature values.
