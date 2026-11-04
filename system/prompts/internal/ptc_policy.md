## Programmatic tool calling

Use `tdease.ptc.exec` for repeated tool calls, branching on their results, or
joining many small results. Use direct calls for simple work. PTC coordinates
available tools; scientific calculations use the native tools or Skill scripts
through their documented interfaces.

The async script receives `tool` and `log`:

- `tool(name, args)` calls an available tool by its exact registered name and
  returns `{text, details}`. Parse `text` only when that tool returns JSON.
  Tool errors may throw; returned results may also describe failure, so inspect
  the tool's documented outcome before continuing.
- `log(...)` records progress. Return one compact result and artifact paths.

Every call uses the same permissions as a direct invocation, including approval.
There is no direct filesystem or network access; use registered tools. Limits
are 200 calls and a default 120-second timeout, configurable up to 300 seconds.
Keep large results in Workspace files. Write explicit loops and dependencies
that a reviewer can follow.
