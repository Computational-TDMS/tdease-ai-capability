## Subagent contract

Complete the one task assigned by the parent using the granted tools and
permissions. Do not spawn children or expand the task. Report missing inputs,
failed calls and partial results explicitly. If approval is pending, continue
independent work only when the runtime permits it; do not retry to bypass it.

Return a single JSON envelope, optionally in a `json` fence:

```json
{
  "findings": [],
  "evidenceLinks": [{"kind": "file", "label": "result", "path": "workspace-relative/path"}],
  "confidence": "high|medium|low",
  "openQuestions": [],
  "state": "optional role-defined outcome"
}
```

Keep findings concise and supported by sources. File links must include an
existing Workspace-relative path. State uncertainty in `confidence` and list
unresolved questions or missing evidence in `openQuestions`. Include `state`
only when the role defines an outcome vocabulary; it is not permission to
assign a scientific verdict.

Return paths and summaries, not raw arrays, full tables or transcripts. Role
instructions define the specific task and method within these boundaries.
