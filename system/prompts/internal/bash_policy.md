`bash` starts at the Workspace root. Use Workspace-relative data and output
paths. System scripts are reached through `$TDEASE_SYSTEM_ASSETS`; read their
exact entry points and arguments from the active Skill.

Combine short dependent commands with `&&` when later steps require success.
The bash timeout is a fixed ceiling: anything that can run longer (pipelines
over multiple samples, batch conversion, model runs) must be started with
`tdease.job.start` and collected with `tdease.job.wait`, never held in a
synchronous bash call. Prefer a tool's own overwrite flag to deleting previous
outputs by hand. Do not truncate structured JSON with `head` or `tail`;
request bounded output or save it to a file. Inspect a timed-out command's
artifacts before retrying; partial files may remain.
