# Engineering contract

## Ownership

```text
source artifact
  -> domain adapter / gen-data.py (derivation and assertions)
  -> {data, meta} contract
  -> build.py + example/build.json (compilation only)
  -> self-contained HTML + manifest
  -> browser validation
  -> optional PNG/SVG
```

The source artifact is the truth. The JSON contract is the exact plotted data.
`chart.js` is presentation. `build.py` compiles these inputs and does not repair
or reinterpret them.

## Example Package

Every direct child of `examples/` contains:

- `chart.js`: reads `window.DATA` and renders into `#chart`.
- `data.json`: a small representative contract used by the acceptance audit.
- `build.json`: dependency declaration with only these keys:

```json
{
  "bundle": "d3",
  "scripts": ["../../components/example.js"],
  "journal": "nature"
}
```

`bundle` is `d3`, `3dmol`, or `none`; `journal` is `nature`, `cell`,
`analytical-chemistry`, or `null`. Script paths must resolve inside the shipped
chart-forge directory. Unknown fields and missing providers fail the build.

## Data Contract

JSON must be an object with `data` and optional object-valued `meta`.

- `data`: exact chart-ready values; missing values are `null`, never fabricated 0.
- `meta.source`: source artifact identity, not a generic label.
- `meta.transform`: derivation summary when values are transformed.
- `meta.thresholds`: cutoffs and adjustment methods that affect interpretation.
- `meta.colour_scale`: named scale when color carries quantitative meaning.
- `meta.n`: the scientifically meaningful count. The compiler infers it only
  for a flat list and never from the number of keys in an object.

Inline JSON is escaped as data, so source text cannot terminate the HTML script
element. CSV input is parsed as a flat list with conservative numeric coercion.

## Outputs

- `<stem>.html`: primary, offline, interactive asset.
- `<stem>.manifest.json`: provenance, input digests, build parameters, caption.
- `<stem>.data.csv`: written when the plotted contract is a flat row list.
- PNG/SVG: optional exports, never substitutes for the HTML/data provenance pair.
  SVG is emitted only for pure-SVG figures; canvas/hybrid/WebGL output is PNG.

A domain adapter that writes a richer scientific manifest passes
`--no-manifest` to `build.py` and writes the one authoritative sidecar itself.

## Failure Policy

Fail before or during build for invalid JSON shape, unsafe example path, missing
build spec, unknown spec field, missing dependency, or missing template token.
Fail browser validation for script errors, readiness timeout, absent/zero-sized
graphic, or export failure. Never downgrade an explicitly requested export to a
warning.

## Adding a Template

1. Start from the closest shipped example; do not copy from an unshipped runtime
   path.
2. Add the three required files and an adapter only when the source format needs it.
3. Keep derivation out of `chart.js`; expose thresholds and missingness in `meta`.
4. Add the route once to [`../INDEX.md`](../INDEX.md).
5. Run `python scripts/audit.py --browser` and inspect the representative PNG.

This complete gate is for template/compiler maintenance. Routine figure use
checks the produced figure, without rebuilding every shipped example.
