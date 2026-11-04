# Changelog

## 0.1.0 - 2026-07-29

- Package the TopPICR Python rewrite with a modern `src/` layout.
- Add the shared `TopPICPipeline` Python API.
- Add the `toppic run` command-line interface.
- Route the optional Streamlit UI through the shared service.
- Add structured inputs, configuration, progress, results, errors, manifests,
  and guarded CSV/Parquet export.
- Separate fast CI regression tests from optional external JPR validation.
- Document platform integration, formats, architecture, scientific validation,
  and known differences.
- Adopt GPL-3.0-or-later with upstream attribution.
