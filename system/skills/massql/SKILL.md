---
name: massql
description: Query mzML spectra with SELECT predicates and plot TIC, XIC or MS2 surveys.
---

# MassQL — query and view mzML data

Use `tdeasenexus-toolbox select --statement '<SQL>' FILE` for spectrum
selection. Use `scripts/plot.py` for chromatograms from those queries; other
measurement questions belong to the relevant scientific Skill.

## Grammar

```
SELECT spectra FROM mzml WHERE <predicate> [AND <predicate>...] [LIMIT n];
```

Predicates (case-insensitive; `rt` = `retention_time`):

| field | form | meaning |
|---|---|---|
| `ms_level` | `= 2` | MS level |
| `rt` | `BETWEEN 10 AND 20` | retention time (minutes) |
| `precursor_mz` | `WITHIN 502.29 PPM 10` / `BETWEEN a AND b` | precursor ion m/z (MS2+) |
| `precursor_neutral_mass` | `BETWEEN a AND b` | precursor neutral mass = selected m/z × z − z·proton from the scan's own charge header — one interval across charges (the top-down family window; charge-less scans never match) |
| `charge` | `= 3` | precursor charge |
| `activation` | `= 'HCD'` | fragmentation method |
| `product_mz` | `WITHIN 189.12 PPM 10` / `BETWEEN a AND b` | a peak exists at this m/z — **the XIC predicate**, works on any ms_level |
| `neutral_loss` | `BETWEEN 40 AND 41` | precursor m/z − peak m/z in range |

`WITHIN <center> DA <tol>` and `WITHIN <center> PPM <tol>` both expand to a
range. `LIMIT` defaults to 100 (max 10000). For chromatograms, request
`LIMIT 10000` and check `truncated`; split the RT window if the result is still
truncated. Do not present a capped query as a complete chromatogram.

## What comes back

```json
{ "plan": {...}, "scannedSpectrumCount": N, "matchedSpectrumCount": N,
  "truncated": false, "spectra": [{
    "scan": "7235", "msLevel": 1, "retentionTime": 12.4,
    "totalIonCurrent": 7.2e7, "peakCount": 178, "centroided": true,
    "precursor": {...}, "matchedProducts": [{"mz": 189.1234, "intensity": 9000}]
  }] }
```

`totalIonCurrent` is always present (header value, or summed intensities when
the header omits it). `peakCount` is `null` unless the statement forced peak
decoding (`product_mz`/`neutral_loss`) — that null is the fast path, not an
error. `matchedProducts` lists the peaks that satisfied the m/z predicates.

## Canonical statements

```sql
-- TIC (total ion chromatogram)
SELECT spectra FROM mzml WHERE ms_level = 1 LIMIT 10000;

-- XIC at MS1 (extracted ion chromatogram)
SELECT spectra FROM mzml WHERE ms_level = 1 AND product_mz WITHIN 502.29 PPM 10 LIMIT 10000;

-- XIC of a fragment across MS2 scans of one precursor
SELECT spectra FROM mzml WHERE ms_level = 2 AND precursor_mz WITHIN 502.29 PPM 10 AND product_mz WITHIN 189.12 PPM 10 LIMIT 10000;

-- MS2 survey in an RT window
SELECT spectra FROM mzml WHERE ms_level = 2 AND rt BETWEEN 10 AND 20 LIMIT 200;

-- Spectra showing a neutral loss (e.g. precursor − 40.016 ± 0.5)
SELECT spectra FROM mzml WHERE neutral_loss BETWEEN 39.5 AND 40.5 LIMIT 200;
```

## Plotting

`scripts/plot.py` runs one or more statements and draws them on one shared
RT axis (traces are auto-classified: rows with `matchedProducts` plot summed
matched intensity = XIC; otherwise `totalIonCurrent` = TIC):

```bash
uv run --script "$TDEASE_SYSTEM_ASSETS/skills/massql/scripts/plot.py" \
  --mzml runs/sample/data.mzML \
  --statement 'SELECT spectra FROM mzml WHERE ms_level = 1 LIMIT 10000;' --label TIC \
  --statement 'SELECT spectra FROM mzml WHERE ms_level = 1 AND product_mz WITHIN 502.29 PPM 10 LIMIT 10000;' --label 'XIC 502.29' \
  --out-dir runs/sample/analysis/massql/
```

Outputs into `--out-dir` (self-contained HTML figures): `<stem>.html`
(self-contained figure on one shared RT axis), `<stem>.data.json` (the
plotting contract), `<stem>.manifest.json` (statements, binary, counts —
cite it as provenance). A PNG is a client-side export button inside the
figure, never a script product. Embed the HTML in the reply with a
relative-path markdown image (`![](path)`); it renders inline in the
conversation.

For tabular inspection use `select` directly; `plot.py` adds the figure.
