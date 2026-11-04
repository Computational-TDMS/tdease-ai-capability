# Vendored D3 assets

This directory contains the offline runtime assets used by ChartForge:

- `vendor/d3.min.js` — primary D3 runtime;
- `vendor/d3-sankey.min.js` — optional Sankey layout;
- `vendor/topojson-client.min.js` and `vendor/countries-110m.json` — optional
  map support.

Templates load these shipped assets through `build.json`. They do not use a
CDN or a second chart runtime. The license files stay beside the vendored
assets.

See [`../../INDEX.md`](../../INDEX.md) for routes and
[`../../paradigms/html-pipeline.md`](../../paradigms/html-pipeline.md) for the
offline build path.
