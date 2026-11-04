---
name: packages
description: Discover, verify and install Node or Workspace Skill packages; check updates while preserving local edits.
---

# Packages

Use the package index for extension discovery and installation:
https://raw.githubusercontent.com/Tisen-Ray/TDEaseNexus-Packages/main/index.json

Refresh it for current availability. Entries expose kind (node/skill), id,
version, description, curation, url, size and sha256 under schemaVersion 1.
Preserve the selected entry/source. Curation and a matching digest do not
authenticate an arbitrary publisher.

## Install

A discovery request lists relevant packages; a request for a specific install
authorizes proceeding through the runtime permission flow. Stage under
.tdease/packages/<id>-<version>/:

    curl -L -o pkg.zip "<url>"
    sha256sum pkg.zip

Compare with the index sha256 before extracting with tar -xf pkg.zip. A mismatch
aborts installation; report expected and actual values.

- Node: the extracted folder contains rules.smk. Pass it to tdease.nodes.install.
  This installs into the local Backend catalog; a remote Backend requires
  installation on that host.
- Skill: extract into .tdease/skills/<id>/ in the Workspace. It is discoverable
  next turn; installed System Skills remain product assets.

Report package/version, source, verified digest and installation path.
On failure, identify the failed boundary rather than disguising it with a
renamed package or unverified substitute.

## Updates

Installed files may contain user edits. An update check reports availability;
replacement needs an upgrade request. Use nodes.list/get and any saved index
entry for installed identity. Archive SHA-256 does not describe later edits
or equal a Node tree digest; local modification status can be unknown.

Preserve edits when upgrading. One Node id holds one version; earlier Run
provenance is historical. Stage Skill upgrades beside the old version for
comparison unless direct replacement was requested.
