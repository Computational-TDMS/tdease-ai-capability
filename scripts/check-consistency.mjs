#!/usr/bin/env node
/**
 * TDEase capability-package consistency lint.
 *
 * The `system/` tree carries four kinds of truth — prompts, skills, reference
 * data, and vendored code — that reference each other by hand-written strings
 * (script names, directory names, tool names, result vocabulary). This script
 * is the guardrail: it cross-checks those strings against the filesystem and
 * against the consuming Electron runtime, so a drift is caught at the editor,
 * not the moment the agent tries to run a deleted script.
 *
 * Zero dependencies (Node 18+). Run from anywhere:
 *
 *   node scripts/check-consistency.mjs [--frontend <path>]
 *
 * The frontend root defaults to the sibling `../TDEase-FrontEnd`; pass
 * `--frontend` to point elsewhere, or omit it to skip the cross-repo checks.
 *
 * Exit code: 1 when any ERROR finding is present (dangling reference, dead
 * asset, contract break), 0 otherwise. WARN/INFO do not affect the exit code.
 */

import { readFile, readdir, stat } from 'node:fs/promises'
import { execFileSync } from 'node:child_process'
import path from 'node:path'
import { fileURLToPath } from 'node:url'

// ── Paths ────────────────────────────────────────────────────────────────────

const here = path.dirname(fileURLToPath(import.meta.url))
const capRoot = path.resolve(here, '..')

const args = process.argv.slice(2)
function argValue(flag) {
  const i = args.indexOf(flag)
  return i >= 0 ? args[i + 1] : undefined
}
const frontendArg = argValue('--frontend')
const frontendRoot = frontendArg ? path.resolve(frontendArg) : path.resolve(capRoot, '..', 'TDEase-FrontEnd')

// ── Findings ─────────────────────────────────────────────────────────────────

/** @type {{severity:'ERROR'|'WARN'|'INFO', check:string, file?:string, detail:string}[]} */
const findings = []
function report(severity, check, detail, file) {
  findings.push({ severity, check, detail, file })
}
const err = (check, detail, file) => report('ERROR', check, detail, file)
const warn = (check, detail, file) => report('WARN', check, detail, file)
const info = (check, detail, file) => report('INFO', check, detail, file)

// ── Filesystem helpers ───────────────────────────────────────────────────────

async function walk(dir, out = []) {
  let entries
  try { entries = await readdir(dir, { withFileTypes: true }) } catch { return out }
  for (const e of entries) {
    if (e.name === 'node_modules' || e.name === '.git' || e.name === '.venv' || e.name === '__pycache__') continue
    const full = path.join(dir, e.name)
    if (e.isDirectory()) await walk(full, out)
    else out.push(full)
  }
  return out
}

function rel(file) { return path.relative(capRoot, file).split(path.sep).join('/') }

async function exists(p) { try { await stat(p); return true } catch { return false } }

async function readIf(p) { try { return await readFile(p, 'utf8') } catch { return null } }

function norm(p) { return p.split(path.sep).join('/') }

/** Normalize CRLF so regexes written against `\n` match on Windows checkouts. */
function nl(s) { return s.replace(/\r\n/g, '\n') }

// ── Data collection ──────────────────────────────────────────────────────────

const systemFiles = (await walk(path.join(capRoot, 'system'))).map(rel)

const docFiles = systemFiles.filter((f) => f.endsWith('.md'))
const pyFiles = systemFiles.filter((f) => f.endsWith('.py'))
// Skills may own .mjs tooling too (e.g. chart-forge export.mjs) — index both.
const scriptFiles = systemFiles.filter((f) => f.endsWith('.py') || f.endsWith('.mjs'))
const skillDirs = systemFiles
  .filter((f) => f.startsWith('system/skills/') && f.includes('/SKILL.md'))
  .map((f) => f.slice('system/skills/'.length, f.indexOf('/SKILL.md')))
  .filter((name) => !name.includes('/subagents/'))
const skillTopDirs = [...new Set(skillDirs.map((n) => n.split('/')[0]))]
const filesByBasename = new Map()
for (const f of scriptFiles) {
  const base = f.split('/').pop()
  if (!filesByBasename.has(base)) filesByBasename.set(base, [])
  filesByBasename.get(base).push(f)
}
const systemSet = new Set(systemFiles)

let gitTracked = null
try {
  gitTracked = execFileSync('git', ['-C', capRoot, 'ls-files'], { encoding: 'utf8' })
    .split('\n').filter(Boolean).map(norm)
} catch { /* not a git repo or git unavailable */ }

// Build tooling (e.g. build-unimod-database.mjs) lives in the frontend, not in
// system/. Collect its script basenames so `.mjs` references resolve correctly
// instead of being reported as dangling capability scripts.
const frontendBasenames = new Set()
for (const sub of ['scripts', 'electron']) {
  for (const f of await walk(path.join(frontendRoot, sub))) {
    frontendBasenames.add(f.split(path.sep).pop())
  }
}

// ── 1. Prompt assets the Electron runtime loads must exist ───────────────────

async function checkPromptAssets() {
  const typesPath = path.join(frontendRoot, 'src', 'agent', 'prompt', 'types.ts')
  const types = await readIf(typesPath)
  if (!types) { warn('prompt-assets', `frontend not found at ${rel(frontendRoot)}; skipping cross-repo prompt check`); return }
  const assetRefs = [...types.matchAll(/:\s*'prompts\/internal\/([\w.]+\.md)'/g)].map((m) => `prompts/internal/${m[1]}`)
  const consumed = new Set(['prompts/internal/subagent_base.md', ...assetRefs])
  for (const a of consumed) {
    if (!systemSet.has(`system/${a}`)) err('prompt-assets', `runtime consumes ${a} but it does not exist under system/`)
  }
  // Dead prompt files: present but never referenced by the runtime or subagent registry.
  const promptFiles = systemFiles.filter((f) => f.startsWith('system/prompts/internal/') && f.endsWith('.md'))
  for (const f of promptFiles) {
    const name = f.slice('system/prompts/internal/'.length)
    if (!consumed.has(`prompts/internal/${name}`)) warn('prompt-assets', `prompt file ${f} is not referenced by PROMPT_ASSETS or subagent_registry (dead?)`, f)
  }
}

// ── 2. Script references in docs must resolve ────────────────────────────────

// Scripts documented as removed; keep them out of the dangling-report so the
// migration tables don't fire false positives, but note they were seen.
const DOCUMENTED_REMOVED = new Set(['unimod_lookup.py', 'uniprot_fetch.py', 'psimod_lookup.py'])

async function checkScriptRefs() {
  for (const doc of docFiles) {
    // Only treat prompts, skill bodies, references and subagent roles as
    // authoritative instruction surface — not vendored upstream docs
    // (toppic library internals under scripts/toppic/, and anything under a
    // skill's assets/, which ships verbatim for reference).
    if (doc.includes('scripts/toppic/') || /\/assets\//.test(doc)) continue
    const text = nl(await readIf(path.join(capRoot, doc)) ?? '')
    if (!text) continue
    const ownerSkill = doc.startsWith('system/skills/')
      ? doc.slice('system/skills/'.length).split('/')[0]
      : null
    for (const match of text.matchAll(/\b([A-Za-z0-9_][A-Za-z0-9_./-]*\.(?:py|mjs))\b/g)) {
      const token = match[1].replace(/\\/g, '/')
      const base = token.split('/').pop()
      if (DOCUMENTED_REMOVED.has(base)) continue
      const isMjs = base.endsWith('.mjs')
      // A doc-scoped 'scripts/NAME.py' path must live in the owning skill's
      // scripts/ directory (skills are self-contained). Build-tooling .mjs
      // references resolve against the frontend, so skip that branch.
      if (ownerSkill && !isMjs && token.startsWith('scripts/')) {
        const localPath = `system/skills/${ownerSkill}/${token}`
        if (!systemSet.has(localPath)) {
          const elsewhere = filesByBasename.get(base)?.[0] ?? null
          if (elsewhere) warn('script-ref', `references ${token} but it is not in this skill; found at ${elsewhere}`, doc)
          else err('script-ref', `references ${token} (${base}) which does not exist anywhere under system/`, doc)
        }
        continue
      }
      // Bare or path-qualified token: resolve by basename anywhere under
      // system/, plus the frontend build scripts for .mjs.
      const found = filesByBasename.has(base) || (isMjs && frontendBasenames.has(base))
      if (!found) {
        err('script-ref', `references ${token} which does not exist anywhere under system/${isMjs ? ' or the frontend build scripts' : ''}`, doc)
      }
    }
  }
}

// ── 3. README must list every skill and every system/ subdirectory ──────────

async function checkReadme() {
  const readme = nl(await readIf(path.join(capRoot, 'README.md')) ?? '')
  if (!readme) { err('readme', 'README.md missing'); return }
  for (const skill of skillTopDirs) {
    if (!readme.includes(skill)) warn('readme', `skill dir system/skills/${skill}/ is not mentioned in README.md`)
  }
  // system/ subdirectories declared in the Layout block must exist. Parse the
  // block indentation-aware: only entries indented one level under `system/`
  // (4 spaces) count; top-level `rust/`/`system/` and skill names under
  // `skills/` do not.
  const fence = readme.match(/## Layout\s*\n```\n([\s\S]*?)\n```/)?.[1] ?? ''
  let inSystem = false
  for (const line of fence.split('\n')) {
    // A directory entry: leading spaces, a name, `/`, then an optional
    // right-aligned description (which `trim().endsWith('/')` would miss).
    const m = line.match(/^(\s*)([\w-]+)\//)
    if (!m) continue
    const indent = m[1].length
    const name = m[2]
    if (indent === 2 && name === 'system') { inSystem = true; continue }
    if (inSystem && indent <= 2) { inSystem = false; continue }
    if (!inSystem || indent !== 4) continue
    if (!systemFiles.some((f) => f.startsWith(`system/${name}/`))) {
      warn('readme', `README Layout lists system/${name}/ but no such directory exists`, 'README.md')
    }
  }
}

// ── 4. allowed_tool_refs must use bare command ids (no tdease. prefix) ───────

async function checkToolRefs() {
  const roleSk = systemFiles.filter((f) => /subagents\/[^/]+\/SKILL\.md$/.test(f) || f.endsWith('/SKILL.md'))
  for (const f of roleSk) {
    const text = nl(await readIf(path.join(capRoot, f)) ?? '')
    if (!text || !text.includes('allowed_tool_refs')) continue
    const block = text.match(/---\n([\s\S]*?)\n---/)?.[1]
    if (!block) continue
    for (const m of block.matchAll(/^\s*-\s+([^\s#]+)\s*$/gm)) {
      const ref = m[1]
      if (ref.startsWith('tdease.') || ref.includes(':')) {
        err('tool-ref', `allowed_tool_refs value "${ref}" uses a prefixed/scoped name; the convention is the bare command id (e.g. "spectra.query")`, f)
      }
    }
  }
}

// ── 5. Hygiene: committed binaries and stray cache dirs ─────────────────────
//
// NOTE: duplicate script basenames across skills (e.g. `_paper.py` in
// isoforma and mass-shift-analysis) are deliberate — each skill is
// self-contained and must not reach into another skill's scripts/. This lint
// therefore does not flag cross-skill script duplication.

async function checkHygiene() {
  if (gitTracked) {
    const binaries = gitTracked.filter((f) => /\.(db|sqlite)$/.test(f))
    for (const f of binaries) warn('hygiene', `committed binary reference ${f} — prefer committing the source (XML/OBO) and generating at build`, f)
  }
  for (const f of systemFiles) {
    if (f.includes('.ruff_cache/')) { info('hygiene', `stray .ruff_cache ${f} (not in .gitignore)`, f); break }
  }
}

// ── 6. Vendored copies of shared modules must stay in sync ───────────────────
//
// Some modules are duplicated across skill packages on purpose (each skill is
// self-contained and must not import across packages). That only works while
// the copies cannot drift silently: the first listed path is the canonical
// owner; every other copy must match it byte for byte (CRLF-normalized). Edit
// the canonical copy first, then propagate to all locations.

const VENDORED_MODULES = {
  // (currently empty: isoforma's _paper.py lost its toppicr-py twin when the
  // dead matplotlib pathway was pruned from that skill; re-register copies here
  // when a module is vendored into multiple skills again)
}

async function checkVendoredSync() {
  for (const [base, copies] of Object.entries(VENDORED_MODULES)) {
    const contents = []
    for (const f of copies) {
      const text = await readIf(path.join(capRoot, f))
      if (text == null) { err('vendored-sync', `${base}: declared copy ${f} does not exist`, f); contents.push(null) }
      else contents.push(nl(text))
    }
    const canonical = contents[0]
    if (canonical == null) continue
    for (let i = 1; i < copies.length; i++) {
      if (contents[i] != null && contents[i] !== canonical) {
        err('vendored-sync', `${base}: copy ${copies[i]} drifted from canonical owner ${copies[0]} — edit the owner, then propagate to all ${copies.length} locations`, copies[i])
      }
    }
  }
}

// ── Run ──────────────────────────────────────────────────────────────────────

await checkPromptAssets()
await checkScriptRefs()
await checkReadme()
await checkToolRefs()
await checkHygiene()
await checkVendoredSync()

// ── Report ───────────────────────────────────────────────────────────────────

const order = { ERROR: 0, WARN: 1, INFO: 2 }
findings.sort((a, b) => order[a.severity] - order[b.severity] || (a.file ?? '').localeCompare(b.file ?? ''))

const counts = { ERROR: 0, WARN: 0, INFO: 0 }
for (const f of findings) counts[f.severity]++

console.log(`TDEase capability consistency lint`)
console.log(`  capability: ${capRoot}`)
console.log(`  frontend:   ${frontendRoot}`)
console.log(`  files: ${systemFiles.length} · skills: ${skillTopDirs.length} · python: ${pyFiles.length}`)
console.log('')

for (const f of findings) {
  const loc = f.file ? `${f.file}` : '-'
  console.log(`[${f.severity}] ${f.check}`)
  console.log(`    ${loc}: ${f.detail}`)
}

console.log('')
console.log(`Summary: ${counts.ERROR} error · ${counts.WARN} warn · ${counts.INFO} info`)

process.exit(counts.ERROR > 0 ? 1 : 0)
