# Repository Instructions

- Use `TopPICPipeline` for end-to-end integrations.
- Keep scientific algorithms in their existing focused modules.
- Do not duplicate orchestration in CLI or Streamlit.
- Run `python -m pytest -m "not jpr"` before committing.
- The `jpr` suite requires `TOPPIC_JPR_DATA_DIR`.
- Never add `TopPICR/`, `TopPICR_reproducible_code_for_JPR/`,
  `_pipeline_data.pkl`, generated results, local settings, or credentials.
- Update README and relevant docs when changing the public API, CLI, input,
  output, progress, exception, or manifest contract.
- Preserve GPL-3.0-or-later licensing and upstream TopPICR attribution.
- Treat uploaded data as caller-owned local files; the package must not delete
  inputs.
- Prefer small synthetic fixtures for CI. Do not download the JPR dataset in
  GitHub Actions.
