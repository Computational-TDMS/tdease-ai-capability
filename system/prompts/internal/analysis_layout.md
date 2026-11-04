Keep one directory per analysis question — `runs/<run>/analysis/<slug>/`
inside a Run, `analysis/<slug>/` outside one — with a short kebab-case slug.
Reuse the directory when continuing the same question, including from a new
session; open a new one for a different question. Keep scripts, derived data
and figures together; pass bulk data between tools by file path, using
Workspace-relative paths and tool-returned paths verbatim.

For Backend inputs, call `tdease.runs.pullFile` with the result `revision`
and the `analysisSlug`. The runtime records the source in that directory's
`sources.json`; do not edit this generated record, and use the returned
`written` path — a conflicting version may be isolated under the analysis
rather than placed at the Run root, and a `conflict` result supplies no
usable new download. Recheck recorded sources before continuing an earlier
analysis; changed content must not retain the old provenance.

Invoke Skill scripts exactly as the loaded Skill specifies, through the
`$TDEASE_SYSTEM_ASSETS` or `.tdease/skills` paths. Asset paths hold installed
code; mutable outputs stay in the Workspace.
