# /// script
# requires-python = ">=3.11"
# dependencies = ["biopython>=1.79"]
# ///
"""Fetch a protein structure, score site solvent accessibility, render one HTML figure.

One pass, three outputs:

  fetch   AlphaFold model for a UniProt accession (numbering = UniProt 1:1),
          or an RCSB entry by id. Cached under reference/structures/ with a
          provenance sidecar; --refresh re-downloads.
  score   per-residue SASA (Shrake-Rupley) mapped onto requested site spans,
          with relative exposure verdicts (Tien 2013 theoretical maxima).
  render  a self-contained 3Dmol HTML — cartoon, pLDDT coloring for AlphaFold
          models, site spheres + labels, SASA table. Inline the viewer library
          from a local cache so the figure renders offline; CDN fallback.

Prints one JSON object on stdout: structure provenance, per-site residue
verdicts, and the figure path when --out is given.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import sys
import urllib.request
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

UA = "TDEase-structure-evidence/1.0"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction"
RCSB_DOWNLOAD = "https://files.rcsb.org/download"
MOL_CACHE_CANDIDATES = [
    "https://3Dmol.org/build/3Dmol-min.js",
    "https://cdn.jsdelivr.net/npm/3dmol@2.4.0/build/3Dmol-min.js",
]

# Theoretical maximal ASA (A^2), Tien et al. 2013 (Gly-X-Gly, tripeptide).
MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

EXPOSED_CUTOFF = 0.30
BURIED_CUTOFF = 0.10


def http_get(url: str, timeout: int = 60) -> bytes:
    request = urllib.request.Request(url, headers={"user-agent": UA})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def fail(message: str) -> None:
    print(json.dumps({"ok": False, "error": message}), file=sys.stderr)
    raise SystemExit(1)


# ── fetch ────────────────────────────────────────────────────────────────────

def resolve_alphafold(accession: str) -> tuple[str, str, dict]:
    """(download_url, cache_key, provenance) for the first AlphaFold model."""
    try:
        entries = json.loads(http_get(f"{ALPHAFOLD_API}/{accession}"))
    except Exception as error:  # noqa: BLE001 - report, don't traceback
        fail(f"AlphaFold API unreachable for {accession}: {error}")
    if not entries:
        fail(f"AlphaFold has no model for {accession}")
    entry = entries[0]
    return entry["pdbUrl"], f"AF-{accession}", {
        "source": "alphafold",
        "accession": accession,
        "url": entry["pdbUrl"],
        "uniprot_start": entry.get("uniprotStart"),
        "uniprot_end": entry.get("uniprotEnd"),
        "date_modified": entry.get("latest", {}).get("dateModified"),
    }


def resolve_rcsb(pdb_id: str) -> tuple[str, str, dict]:
    url = f"{RCSB_DOWNLOAD}/{pdb_id}.pdb"
    return url, pdb_id.upper(), {"source": "rcsb", "pdb_id": pdb_id.upper(), "url": url}


def cached_fetch(url: str, cache_key: str, cache_dir: Path, refresh: bool) -> tuple[str, dict]:
    """Download to cache unless a fresh-enough copy exists; return (path, meta)."""
    target = cache_dir / f"{cache_key}.pdb"
    sidecar = cache_dir / f"{cache_key}.meta.json"
    if target.exists() and not refresh:
        meta = json.loads(sidecar.read_text()) if sidecar.exists() else {}
        return str(target), meta
    try:
        data = http_get(url)
    except Exception as error:  # noqa: BLE001
        fail(f"structure download failed ({url}): {error}")
    text = data.decode("ascii", errors="replace")
    cache_dir.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding="ascii", errors="replace")
    return str(target), {}


def write_sidecar(cache_key: str, cache_dir: Path, provenance: dict, residue_count: int, digest: str) -> None:
    sidecar = cache_dir / f"{cache_key}.meta.json"
    sidecar.write_text(json.dumps({
        **provenance,
        "fetched_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "residues": residue_count,
        "md5": digest,
    }, indent=2))


# ── structure model ──────────────────────────────────────────────────────────

def load_residues(pdb_path: str, chain_id: str | None) -> tuple[list[dict], dict]:
    """Per-residue records from one polymer chain, SASA attached.

    The chain is part of the measurement's identity, not an implementation
    detail: an explicit --chain is honored verbatim, the default picks the
    longest chain but records that choice in the output, and a tie for
    longest is refused — a successful run must mean one known object was
    analyzed. Insertion codes are kept so author-numbered experimental
    entries stay unambiguous.
    """
    structure = PDBParser(QUIET=True).get_structure("s", pdb_path)
    ShrakeRupley().compute(structure, level="R")
    chains = {chain.id: chain for chain in structure.get_chains()}
    lengths = {
        cid: sum(1 for residue in chain if residue.id[0] == " ")
        for cid, chain in chains.items()
    }
    if not any(lengths.values()):
        fail(f"{pdb_path}: no standard polymer residues found")
    menu = ", ".join(f"{cid}({count})" for cid, count in sorted(lengths.items()))
    if chain_id is not None:
        if not lengths.get(chain_id):
            fail(f"chain {chain_id!r} not found in {pdb_path}; available: {menu}")
        selected, selection = chain_id, "explicit"
    else:
        longest = max(lengths.values())
        tied = sorted(cid for cid, count in lengths.items() if count == longest)
        if len(tied) > 1:
            fail(f"{pdb_path} has chains tied for longest ({', '.join(tied)}); "
                 f"pass --chain to declare which one to analyze; available: {menu}")
        selected, selection = tied[0], "longest"
    residues = []
    for residue in chains[selected]:
        if residue.id[0] != " ":  # het/water
            continue
        record = {
            "resno": residue.id[1],
            "aa": residue.get_resname(),
            "sasa": round(float(residue.sasa), 2),
        }
        if residue.id[2] and residue.id[2] != " ":
            record["icode"] = residue.id[2]
        residues.append(record)
    return residues, {"id": selected, "selection": selection, "chainLengths": lengths}


def plddt_by_resno(pdb_path: str) -> dict[int, float]:
    """Mean B-factor per residue; for AlphaFold files that is pLDDT."""
    by_resno: dict[int, list[float]] = {}
    for line in Path(pdb_path).read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 60:
            try:
                by_resno.setdefault(int(line[22:26]), []).append(float(line[60:66]))
            except ValueError:
                continue
    return {resno: round(sum(v) / len(v), 1) for resno, v in by_resno.items()}


# ── sites ────────────────────────────────────────────────────────────────────

SITE_RE = re.compile(r"^(\d+)(?:-(\d+))?(?::(.+))?$")


def parse_site(spec: str) -> dict:
    match = SITE_RE.match(spec.strip())
    if not match:
        fail(f"bad --sites spec {spec!r}; expected START[-END][:LABEL]")
    start = int(match.group(1))
    end = int(match.group(2) or start)
    if end < start:
        fail(f"site span {spec!r} ends before it starts")
    return {"start": start, "end": end, "label": (match.group(3) or "").strip()}


def score_site(site: dict, residues: list[dict]) -> dict:
    # Insertion-code variants share a residue number; the first is scored and
    # each keeps its icode in the returned row, so the caller can see both.
    by_resno: dict[int, dict] = {}
    for residue in residues:
        by_resno.setdefault(residue["resno"], residue)
    rows = []
    for resno in range(site["start"], site["end"] + 1):
        residue = by_resno.get(resno)
        if residue is None:
            rows.append({"resno": resno, "aa": None, "in_model": False})
            continue
        rel = residue["sasa"] / MAX_ASA.get(residue["aa"], 200.0)
        verdict = "exposed" if rel >= EXPOSED_CUTOFF else "buried" if rel < BURIED_CUTOFF else "intermediate"
        rows.append({**residue, "rel_sasa": round(rel, 3), "verdict": verdict})
    modeled = [row for row in rows if row.get("in_model", True)]
    tally = {v: sum(1 for row in modeled if row["verdict"] == v) for v in ("exposed", "intermediate", "buried")}
    return {
        **site,
        "residues": rows,
        "in_model": len(modeled),
        **tally,
        "note": None if modeled else "no residue of this span is present in the model",
    }


# ── render ───────────────────────────────────────────────────────────────────

def mol_library(cache_dir: Path) -> str:
    """Inline 3Dmol from a local cache (offline figures); CDN tag as fallback."""
    cached = cache_dir / ".assets" / "3Dmol-min.js"
    if cached.exists():
        return f"<script>\n{cached.read_text(encoding='utf-8', errors='replace')}\n</script>"
    for url in MOL_CACHE_CANDIDATES:
        try:
            data = http_get(url)
        except Exception:  # noqa: BLE001 - try the next mirror
            continue
        cached.parent.mkdir(parents=True, exist_ok=True)
        cached.write_bytes(data)
        return f"<script>\n{data.decode('utf-8', errors='replace')}\n</script>"
    return '<script src="https://3Dmol.org/build/3Dmol-min.js"></script>'


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
  html, body {{ margin: 0; height: 100%; background: #ffffff; color: #1f2933;
    font: 13px/1.45 system-ui, -apple-system, sans-serif; }}
  #viewer {{ width: 100%; height: calc(100% - 168px); position: relative; }}
  header {{ padding: 8px 14px; border-bottom: 1px solid #e5e7eb; }}
  header h1 {{ font-size: 14px; margin: 0 0 2px; }}
  header p {{ margin: 0; color: #6b7280; font-size: 12px; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 12px; }}
  th, td {{ text-align: left; padding: 2px 14px; border-top: 1px solid #f0f1f3; }}
  th {{ color: #6b7280; font-weight: 500; }}
  .exposed {{ color: #b45309; font-weight: 600; }}
  .buried {{ color: #1d4ed8; font-weight: 600; }}
  .intermediate {{ color: #6b7280; }}
  footer {{ padding: 4px 14px; color: #9ca3af; font-size: 11px; }}
</style>
{mol_lib}
</head>
<body>
<header>
  <h1>{title}</h1>
  <p>{subtitle}</p>
</header>
<div id="viewer"></div>
<script type="text/x-pdb" id="pdb-data">{pdb_text}</script>
<script>
const PDB = document.getElementById('pdb-data').textContent;
const SITES = {sites_json};
const PLDDT = {plddt_json};
const IS_AF = {is_af};
// Residue numbers embedded from the Python-side PDB parse — no model
// introspection here (the atoms accessor is not stable across 3Dmol builds).
const RESNOS = {resnos_json};

let element = document.getElementById('viewer');
let viewer = $3Dmol.createViewer(element, {{ backgroundColor: 'white' }});
viewer.addModel(PDB, 'pdb');

function plddtColor(value) {{
  if (value >= 90) return '#0053d6';
  if (value >= 70) return '#65cbfe';
  if (value >= 50) return '#ffe24d';
  return '#ff7d45';
}}

// Cartoon: pLDDT palette for AlphaFold models, plain slate otherwise.
const resnos = RESNOS;
for (const resi of resnos) {{
  const color = (IS_AF && PLDDT[resi] !== undefined) ? plddtColor(PLDDT[resi]) : '#607d8b';
  viewer.setStyle({{ resi }}, {{ cartoon: {{ color }} }});
}}

for (const site of SITES) {{
  const span = [];
  for (let r = site.start; r <= site.end; r++) span.push(r);
  const present = span.filter(r => resnos.includes(r));
  if (present.length === 0) continue;
  viewer.addStyle({{ resi: present }}, {{
    sphere: {{ radius: 1.9, color: '#e11d48', opacity: 0.95 }},
    stick: {{ radius: 0.25, color: '#e11d48' }},
  }});
  viewer.addResLabels({{ resi: present }}, {{
    font: '11px sans-serif', color: '#e11d48', showBackground: false,
  }});
}}
if (SITES.length > 0) {{
  viewer.addSurface($3Dmol.SurfaceType.VDW, {{ opacity: 0.16, color: '#94a3b8' }});
}}
viewer.zoomTo();
viewer.render();
</script>
{table_html}
<footer>{footer}</footer>
</body>
</html>
"""


def site_table(sites: list[dict]) -> str:
    rows = []
    for site in sites:
        span = f"{site['start']}" if site["start"] == site["end"] else f"{site['start']}\u2013{site['end']}"
        for row in site["residues"]:
            if not row.get("in_model", True):
                cell = '<td colspan="4" style="color:#9ca3af">not in model</td>'
            else:
                cell = (
                    f'<td>{row["aa"]}</td>'
                    f'<td>{row["resno"]}</td>'
                    f'<td>{row["sasa"]:.0f} \u00c5\u00b2 / {row["rel_sasa"] * 100:.0f}%</td>'
                    f'<td class="{row["verdict"]}">{row["verdict"]}</td>'
                )
            rows.append(f'<tr><td>{span}</td><td>{site["label"] or "\u2014"}</td>{cell}</tr>')
    if not rows:
        return ""
    return (
        '<table><tr><th>span</th><th>label</th><th>aa</th><th>resno</th>'
        '<th>SASA (rel)</th><th>verdict</th></tr>' + "".join(rows) + "</table>"
    )


def render_html(out_path: Path, pdb_path: str, provenance: dict, residues: list[dict],
                sites: list[dict], plddt: dict[int, float], cache_dir: Path,
                chain_info: dict) -> None:
    is_af = provenance.get("source") == "alphafold"
    title = (f"{provenance.get('accession') or provenance.get('pdb_id')} "
             f"chain {chain_info['id']} \u2014 "
             f"{len(residues)} residues ({provenance.get('source')})")
    subtitle_bits = []
    if sites:
        subtitle_bits.append("rose = localized span; SASA verdicts below")
    if is_af:
        subtitle_bits.append("cartoon colored by pLDDT (blue \u2265 90, yellow \u2265 50, orange < 50)")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(HTML_TEMPLATE.format(
        title=title,
        subtitle="; ".join(subtitle_bits),
        mol_lib=mol_library(cache_dir),
        pdb_text=Path(pdb_path).read_text(),
        sites_json=json.dumps(sites and [dict(start=s["start"], end=s["end"], label=s["label"]) for s in sites] or []),
        plddt_json=json.dumps(plddt if is_af else {}),
        is_af="true" if is_af else "false",
        resnos_json=json.dumps([r["resno"] for r in residues]),
        table_html=site_table(sites),
        footer=f"verdicts: rel SASA \u2265 {EXPOSED_CUTOFF:.0%} exposed, < {BURIED_CUTOFF:.0%} buried",
    ), encoding="utf-8")


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--accession", help="UniProt accession; fetches the AlphaFold model (numbering = UniProt)")
    source.add_argument("--pdb", help="RCSB entry id; experimental structure (author numbering!)")
    parser.add_argument("--chain",
                        help="explicit polymer chain id; when omitted the longest chain "
                             "is analyzed and that choice is recorded — a tie for "
                             "longest is refused instead of silently resolved")
    parser.add_argument("--sites", action="append", default=[], metavar="START[-END][:LABEL]",
                        help="localized span to map and score, e.g. 24-28:+15.9949; repeatable")
    parser.add_argument("--out", type=Path, help="write the interactive HTML figure here (workspace-relative)")
    parser.add_argument("--cache-dir", type=Path, default=Path("reference/structures"),
                        help="structure cache dir, cwd-relative (default: reference/structures)")
    parser.add_argument("--refresh", action="store_true", help="re-download even if cached")
    args = parser.parse_args()

    if args.accession:
        url, cache_key, provenance = resolve_alphafold(args.accession.strip().upper())
    else:
        url, cache_key, provenance = resolve_rcsb(args.pdb.strip())

    pdb_path, meta = cached_fetch(url, cache_key, args.cache_dir, args.refresh)
    residues, chain_info = load_residues(pdb_path, args.chain.strip() if args.chain else None)
    digest = hashlib.md5(Path(pdb_path).read_bytes()).hexdigest()
    provenance["fetched_at"] = meta.get("fetched_at", dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"))
    if not meta:
        write_sidecar(cache_key, args.cache_dir, provenance, len(residues), digest)

    warnings = []
    start = provenance.get("uniprot_start")
    if isinstance(start, int) and start != 1:
        warnings.append(f"model numbering starts at {start}, not 1 — offset your spans before citing positions")

    sites = [score_site(parse_site(spec), residues) for spec in args.sites]
    result = {
        "ok": True,
        "structure": {
            "path": pdb_path,
            "chain": chain_info,
            "residues": len(residues),
            "md5": digest,
            **provenance,
        },
        "sites": sites,
        "warnings": warnings,
    }
    if args.out:
        render_html(args.out, pdb_path, provenance, residues, sites,
                    plddt_by_resno(pdb_path), args.cache_dir, chain_info)
        result["figure"] = str(args.out)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
