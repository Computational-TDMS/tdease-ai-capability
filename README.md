<p align="center">
  <img src="assets/logo.png" width="96" alt="TDEaseNexus logo">
</p>

# tdease-ai-capability

Agent Skills and the Rust mass-spectrometry kernel behind
[TDEaseNexus](https://github.com/Computational-TDMS/TDEaseNexus). Two halves meant to live
together:

| Path | What it is |
|------|------------|
| `system/` | Agent Skills (one `SKILL.md` per skill), internal prompts, UniMod/PSI-MOD reference data, Snakemake templates |
| `rust/` | `tdeasenexus-toolbox` CLI — isotopic envelopes, fragmentation theory, UniMod resolution |

`system/` works with any agent that speaks
[Agent Skills](https://docs.claude.com/en/docs/agents-and-tools/agent-skills)
(Claude, Codex, …): point it at `system/skills/` and let it read the `SKILL.md`
files. Several skills shell out to the CLI for the heavy numerics, so keep the
built binary on `PATH`.

## Build the CLI

Requires Rust ≥1.97.

```bash
cd rust
cargo build --release   # → rust/target/release/tdeasenexus-toolbox
cargo test              # unit + integration tests, incl. the CLI contract suite
```

```bash
./rust/target/release/tdeasenexus-toolbox fragment-theoretical \
  --proforma 'M.(S)[Acetyl]ESSSK.K' --ion-types C --max-charge 1
```

## Skills

`system/skills/` covers mass-shift exploration, isoforma quantification, spectral
and structure evidence, TopPIC pipelines, MassQL traces, scientific figure
compilation, and more. Each skill is self-documenting: its `SKILL.md` states when
to use it and which scripts it runs.

## License

GPL-3.0-or-later — see [LICENSE](LICENSE).
