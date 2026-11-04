---
name: unidec
description: Deconvolve MS1 spectra, centroid peaks or predict isotope distributions with UniDec, IsoDec and IsoGen.
---

# UniDec tools

Use the shipped wrappers for deconvolution and isotope prediction. Cite UniDec
(Marty et al., Anal. Chem. 2015) in publications using its results. These tools
report masses and charge states, not chemical identities.

## Choose the script

| Script | Input | Output |
|---|---|---|
| `unidec_deconvolve.py` | Profile MS1 text or mzML with MS1 summed | Mass, height, area, dscore |
| `findpeaks.py` | A supported 1D spectrum | Centroid m/z and intensity |
| `isodec_deconvolve.py` | Centroided two-column text or MGF | Charge, monoisotopic, peak and average mass |
| `isogen_predict.py` | Neutral mass and optional molecule type | Isotope distribution |

Run the selected script from the Workspace:

```bash
uv run --script "$TDEASE_SYSTEM_ASSETS/skills/unidec/scripts/isodec_deconvolve.py" \
  <analysis-dir>/centroids.mgf > <analysis-dir>/masses.tsv
```

Use `--help` for its arguments. `findpeaks.py` accepts an optional threshold;
`isogen_predict.py` accepts `<neutral_mass> [type]`. Runtime depends on input
size. Use a background job when needed.

## Inputs and saved state

- Text contains m/z and intensity, separated by tabs or spaces, with an
  optional header.
- MGF input combines all blocks and bins peaks at 20 ppm. Export only the scans
  needed for the question. The MGF format alone does not prove centroiding;
  use centroided peaks for IsoDec and centroid profile input with `findpeaks.py`.
- The mzML route sums all MS1 scans. Pre-extract when scan or RT selection matters.
- The 1D engine writes `<stem>_unidecfiles/` beside its input and may reuse a
  saved `_conf.dat`. Stage inputs under the analysis directory so these files
  do not modify the source Run. For a fresh comparison, use a fresh working
  copy and preserve the previous parameter record.

The resolved `_conf.dat` is the parameter authority for a run. Check mass and
charge windows, mass grid spacing, adduct mass, smoothing and peak thresholds
before comparing outputs. The shipped scripts expose a limited CLI; do not
invent flags or edit installed scripts to tune `eng.config`. If additional
configuration is needed, use an explicitly scoped Workspace script and record
its parameters.

## Interpret the return

Successful scripts print TSV followed by `#` report lines containing counts,
warnings and timing. Failures return an `ok: false` JSON error envelope and a
non-zero exit. Read warnings for empty output and noise-dominated charge states.
The 1D peaks are height-sorted; IsoDec peaks are matched-intensity-sorted.

UniScore bands are upstream heuristics, not identification confidence.
Investigate systematic mass offsets from acquisition and calibration evidence;
do not tune a fit to reproduce an expected answer.

Only the listed wrappers are supported by this interface. IM-MS, MetaUniDec
batch and other unwrapped engines require separate implementation work.
Keep the shipped `unidec==8.2.1` dependency pin unless explicitly maintaining
this Skill. Return the peak table, parameter record where produced, warnings
and the limits of the selected source aggregation.
