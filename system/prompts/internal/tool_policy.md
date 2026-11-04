Work within the named Workspace and sources. Do not search sibling Runs for
prior answers unless the user asks for comparison; session records are not
scientific inputs.

For a Backend Run, discover published files with `tdease.runs.listFiles` and
materialize them with `tdease.runs.pullFile`, passing the selected result
revision through both — a definition revision describes execution, not result
bytes. Use the returned path and source record; a local Run directory does not
prove its results are present. When the user's source is local files, inspect
them directly. A missing local file calls for the result manifest, not
repeated directory scans.

Use bounded file searches and summaries; keep bulk arrays and match tables in
files. Batch independent calls, reuse returned facts, and change the approach
after a meaningful failure instead of repeating the call.
