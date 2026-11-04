# Node binding

How a Backend Node becomes executable inside a Run.

## What a Node is

A Node is a self-contained, versioned directory: executable/resources plus a `rules.smk`
that wraps them for Snakemake. The Backend owns the source package; the Client and Agent
see only the rule surface exposed by `tdease.nodes.get`. A Node's
`rules.smk` carries:

- the exact binary invocation (flags calibrated against the binary's own `-?`/`--help`),
- SMK parameters with defaults,
- resource declarations (`threads`, `resources.mem_mb`, …),
- any serialization invariants.

These override any example anywhere else.

## One config-injected include path

A Run binds a Node by id. At dry-run/execution the Backend resolves the binding and
injects the package's absolute catalog path into the engine config as
`config["nodes"][<id>]`; nothing node-related is ever written under the Run directory.
A root Snakefile pulls a bound Node through that mapping:

```python
_TOPPIC = (config.get("nodes") or {}).get("toppic")
if not _TOPPIC:
    raise FileNotFoundError("bound toppic node is missing from this Run")
include: _TOPPIC + "/rules.smk"
```

The native `module` form also reads that mapping and adds namespacing
(`use rule * from toppic as toppic_*`) — but a Snakemake module gets its OWN
empty config by default, so the passthrough line is mandatory. Without it the
module's rules silently resolve empty paths (no binding, no bin_dir) and can
fall through to a host binary:

```python
module toppic:
    snakefile:
        config["nodes"]["toppic"] + "/rules.smk"
    config:
        config          # MANDATORY: modules otherwise see config == {}

use rule * from toppic as toppic_*
```

Rules, binaries, DLLs and resources replay together from the catalog's own bytes — the
catalog IS deployment state, and an execution uses whatever bytes it currently holds.
There is no Backend-root fallback and no environment-variable path. Before submission,
the Agent reads the catalog rule text; only the Backend execution path is executable.

## Executables are addressed by path, never by name

Every node invocation is an explicit path — the `shell:` line each capability
carries in its digest shows it: the absolute catalog path resolved from
`config["nodes"]` for bundled nodes, or a quoted absolute path resolved from
`bin_dir`/`DEFAULT_BIN_DIR` for host-installed ones. When authoring or repairing a
rule, copy that line; never write a bare command name and never rely on PATH. A bare
name resolves against whatever PATH the engine happens to inherit — wrong version,
wrong architecture, or "command not found" — while a path names the exact bytes the
Run bound. `where`/`which` in a Run shell is a diagnostic, not a resolution mechanism.

## Parameter passthrough

`config.yaml` holds one block per Node; the Node's `rules.smk` reads it via
`config.get("<nodeId>", {...})`. The Run-level Snakefile does **not** redeclare Node
parameters — it only supplies `samples:` and any Run-level glue. Changing a Node
parameter means editing `config.yaml`, which makes a new revision. The Backend-reserved
key `nodes` (the injected binding map) is not a parameter block: a `nodes:` block in
`config.yaml` is ignored at execution.

## Digests are provenance, not enforcement

`prepare` records the Node digest in the Run revision — the catalog content the Run
was authored against. The catalog is deployment state: editing or redeploying a
package does not invalidate existing Runs, and resume does not require the digest to
be unchanged. The digest documents what was there; it never gates execution.

## Adding or adapting a Node package

The catalog root is user-manageable: a directory there containing a
`rules.smk` IS a node package (directory name = node id, content digest =
identity), so a package can be added or edited in place and the next Run binds
whatever the catalog now holds. Names starting with `_` are templates/drafts
and are skipped by the catalog scanner. Where the catalog lives is a Backend
deployment fact, and the Agent never assumes it: list what exists through
`tdease.nodes.get`, and ask the user (or the tdeasenexus-docs doctor) when a
package must be installed by hand — a desktop install keeps it beside the
Backend's data root, a remote Backend keeps it on the server, and the
`config["nodes"]` binding contract is identical either way.

`_template/` is the copyable starting point for a tool that is installed on the
machine and must not be bundled (vendor installer, site licence, size). Copy it
to `<id>/`, set `NODE_ID`, `EXE_NAME`, `DEFAULT_BIN_DIR` and the rule contract,
and the node calls the installed binary from that directory. The Run's
`config.yaml` overrides the directory per Run:

```yaml
<node-id>:
  bin_dir: 'D:\Other Install'
```

`bin_dir` is the **one** binary-location override, and it works for every node,
bundled or not: the node owns its executable names, the Run names the install
directory. Bundled nodes use their own catalog package by default
(`msconvert.bin_dir`, `topfd.bin_dir`, `toppic.bin_dir`, `topindex.bin_dir`,
`informed-proteomics.bin_dir`). Single-quote Windows paths in YAML — a
double-quoted `"C:\Users\..."` is an invalid escape (`\U`) and the Run's config
fails to load. Prefer a bundled node when one covers the task; a host-install
node is the fallback for tools the product does not ship, and the executed
command line lands in the Run log.

## What lives where

| thing | Client (Run folder) | Backend |
|-------|---------------------|---------|
| Snakefile, config.yaml | yes | — |
| scientific inputs/outputs | yes | — |
| Node packages | never | yes, in the Backend's node catalog (deployment-owned; on a desktop install the `nodes/` directory beside the Backend data root) — reached at execution through `config["nodes"]` absolute paths |
