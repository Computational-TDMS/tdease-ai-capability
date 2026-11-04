# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Fetch the small specialized PTM/motif databases into one local SQLite.

UniProt is deliberately NOT mirrored (too large); it stays a query source.
Mirrored here: dbPTM (per-type experimental files), iPTMnet (3 txt files),
ELM (class + instance tables). Everything is normalized into one schema so
the atomic database queries (`query.py`) never adapt per-database
schemas at query time:

    site          (accession, entry, organism, position, residue, ptm_type,
                   source, pmids, seq_window)
    enzyme_score  (accession, position, enzyme, ptm_type, score)
    protein       (accession, name, organism)
    motif_class   (class_id, pattern, description)
    motif_instance(class_id, accession, start, end)

dbPTM's seq_window column (20 aa around the site) is the asset that makes
motif-cohort scans possible — no database's web interface exposes that query.

Usage:
  python fetch.py [--cache DIR] [--types LIST|default|all]
                         [--skip-elm] [--force]

Default type set excludes dbPTM Phosphorylation/Dephosphorylation (~165 MB
zipped); pass --types phosphorylation to add them, or --types all.
"""
from __future__ import annotations

import argparse
import io
import json
import re
import sqlite3
import sys
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

DBPTM = "https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment"
IPTMNET = "https://research.bioinformatics.udel.edu/iptmnet_data/files/current"
ELM_INDEX = "http://elm.eu.org/elms/elms_index.tsv"
ELM_INSTANCES = "http://elm.eu.org/instances.tsv"

# Canonical hyphen tokens that stay uppercase after normalization ("GPI-anchor").
_ACRONYM_TOKENS = {"gpi"}


def canonical_type(raw: str) -> str:
    """One spelling per chemistry across sources.

    Canonical form follows dbPTM's convention: first letter of the name
    uppercase; single-letter locus/stereo tokens after a hyphen uppercase
    (S-, N-, O-, D-); everything else lowercase ("S-palmitoylation",
    "D-glucuronoylation", "GPI-anchor"). iPTMnet's "S-PALMITOYLATION" and a
    naive .title() ("S-Palmitoylation") both normalize to the same string.
    """
    tokens = raw.strip().lower().split("-")
    first = tokens[0]
    out = [first.upper() if first in _ACRONYM_TOKENS
           else first[:1].upper() + first[1:]]
    for token in tokens[1:]:
        if len(token) == 1 or token in _ACRONYM_TOKENS:
            out.append(token.upper())
        else:
            out.append(token)
    return "-".join(out)

# Compact default: every small dbPTM type relevant to top-down mass shifts.
# The 165 MB Phosphorylation/Dephosphorylation pair is opt-in.
DEFAULT_TYPES = [
    "Acetylation", "Butyrylation", "Crotonylation", "Formylation",
    "Malonylation", "Methylation", "Propionylation", "Sumoylation",
    "Ubiquitination",
    "Myristoylation", "S-palmitoylation", "N-palmitoylation",
    "O-palmitoylation", "O-palmitoleoylation",
    "Farnesylation", "Geranylgeranylation", "GPI-anchor", "S-diacylglycerol",
    "S-nitrosylation", "Sulfation", "Sulfoxidation", "Oxidation",
    "Hydroxylation", "Deamidation", "Citrullination", "Amidation",
    "Biotinylation", "ADP-ribosylation", "AMPylation", "Carbamidation",
    "N-carbamoylation", "S-carbamoylation", "D-glucuronoylation",
]
OPTIN_TYPES = ["Phosphorylation", "Dephosphorylation"]
DBPTM_DOWNLOAD_PAGE = "https://biomics.lab.nycu.edu.tw/dbPTM/download.php"


def discover_dbptm_types() -> list[str]:
    """All experiment type names from dbPTM's download page (~60; stays
    current without a hardcoded list). Used by `--types full`."""
    import re as _re
    html = fetch(DBPTM_DOWNLOAD_PAGE, timeout=60).decode("utf-8", "replace")
    names = sorted(set(_re.findall(r"download/experiment/([A-Za-z0-9_+-]+)\.zip", html)))
    return [n for n in names if n not in ("Dephosphorylation", "Phosphorylation")] \
        + OPTIN_TYPES


def fetch(url: str, timeout: int = 180) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "tdease-ptm-database/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def init_db(db_path: Path) -> sqlite3.Connection:
    new = not db_path.exists()
    conn = sqlite3.connect(db_path)
    if new:
        conn.executescript(
            """
            CREATE TABLE site (accession TEXT, entry TEXT, organism TEXT,
              position INTEGER, residue TEXT, ptm_type TEXT, source TEXT,
              pmids TEXT, seq_window TEXT);
            CREATE INDEX idx_site_acc ON site(accession);
            CREATE INDEX idx_site_type ON site(ptm_type);
            CREATE TABLE enzyme_score (accession TEXT, position INTEGER,
              enzyme TEXT, ptm_type TEXT, score REAL);
            CREATE INDEX idx_es_acc ON enzyme_score(accession);
            CREATE TABLE protein (accession TEXT PRIMARY KEY, name TEXT, organism TEXT);
            CREATE TABLE motif_class (class_id TEXT PRIMARY KEY, pattern TEXT,
              description TEXT, implies_type TEXT, probability REAL);
            CREATE TABLE motif_instance (class_id TEXT, accession TEXT,
              start INTEGER, end INTEGER);
            CREATE INDEX idx_mi_acc ON motif_instance(accession);
            CREATE TABLE vocab (alias TEXT, source TEXT, canonical TEXT,
              PRIMARY KEY (alias, source));
            CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT);
            """
        )
    return conn


def record_alias(conn: sqlite3.Connection, raw: str, source: str) -> str:
    canonical = canonical_type(raw)
    conn.execute("INSERT OR REPLACE INTO vocab VALUES (?,?,?)",
                 (raw.strip(), source, canonical))
    return canonical


def load_dbptm_type(conn: sqlite3.Connection, ptm_type: str) -> int:
    """One dbPTM zip: TSV rows entry|accession|position|type|pmids|window.

    dbPTM's window is 21 characters: the modified residue sits at index 10,
    flanked by up to 10 residues on each side (dash-padded at the termini).
    Sites dedupe to one row per (accession, position, type) with pmids
    merged, so query-side counts are sites, not evidence lines.
    """
    raw = fetch(f"{DBPTM}/{ptm_type}.zip")
    canonical = record_alias(conn, ptm_type, "dbPTM")
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        name = archive.namelist()[0]
        text = archive.read(name).decode("utf-8", "replace")
    merged: dict[tuple[str, int], list] = {}
    for line in text.strip().splitlines():
        parts = [p.strip() for p in line.split("\t")]
        if len(parts) < 4 or parts[1] == "Accession":
            continue
        try:
            position = int(parts[2])
        except ValueError:
            continue
        window = parts[5] if len(parts) > 5 and parts[5] else None
        pmids = parts[4] if len(parts) > 4 and parts[4] else None
        key = (parts[1], position)
        if key not in merged:
            residue = window[10] if window and len(window) == 21 \
                and window[10] != "-" else None
            merged[key] = [parts[1], parts[0], None, position, residue,
                           canonical, "dbPTM", pmids, window]
        elif pmids and pmids not in (merged[key][7] or ""):
            merged[key][7] = f"{merged[key][7]};{pmids}" if merged[key][7] else pmids
    rows = [entry for entry in merged.values()]
    # delete-then-insert keeps reruns (--types <subset>) idempotent
    conn.execute("DELETE FROM site WHERE source='dbPTM' AND ptm_type=?", (canonical,))
    conn.executemany("INSERT INTO site VALUES (?,?,?,?,?,?,?,?,?)", rows)
    return len(rows)


def load_iptmnet(conn: sqlite3.Connection) -> dict[str, int]:
    """Three txt files: ptm rows + enzyme-score triplets + protein names.

    ptm.txt carries one line per (site, enzyme, pmid) — the same site returns
    many times. Sites dedupe to one row per (accession, position, type) with
    pmids merged; the per-enzyme mapping lives in enzyme_score.
    """
    counts: dict[str, int] = {}
    ptm_text = fetch(f"{IPTMNET}/ptm.txt").decode("utf-8", "replace")
    merged: dict[tuple[str, int, str], list] = {}
    seen_proteins = {}
    for line in ptm_text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        ptm_type, _source, accession, gene, organism, site = parts[:6]
        pmid = parts[9] if len(parts) > 9 and parts[9] else None
        if not re.fullmatch(r"[A-Z]\d+", site):
            continue
        canonical = record_alias(conn, ptm_type, "iPTMnet")
        position = int(site[1:])
        key = (accession, position, canonical)
        if key not in merged:
            merged[key] = [accession, gene, organism, position, site[0],
                           canonical, "iPTMnet", pmid, None]
        elif pmid and pmid not in (merged[key][7] or ""):
            merged[key][7] = f"{merged[key][7]};{pmid}" if merged[key][7] else pmid
        seen_proteins[accession] = (gene, organism)
    site_rows = list(merged.values())
    counts["iptmnet_site"] = len(site_rows)

    score_text = fetch(f"{IPTMNET}/score.txt").decode("utf-8", "replace")
    score_rows = []
    score_aliases: set[str] = set()
    for line in score_text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        accession, site, enzyme, ptm_type, score = parts[:5]
        match = re.fullmatch(r"[A-Z](\d+)", site)
        if not match:
            continue
        try:
            value = float(score)
        except ValueError:
            continue
        if ptm_type not in score_aliases:
            score_aliases.add(ptm_type)
            record_alias(conn, ptm_type, "iPTMnet")
        score_rows.append((accession, int(match.group(1)), enzyme,
                           canonical_type(ptm_type), value))
    counts["iptmnet_score"] = len(score_rows)

    protein_text = fetch(f"{IPTMNET}/protein.txt").decode("utf-8", "replace")
    for line in protein_text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) >= 6 and parts[0]:
            seen_proteins.setdefault(parts[0], (parts[3] if len(parts) > 3 else None,
                                                parts[4] if len(parts) > 4 else None))

    conn.execute("DELETE FROM site WHERE source='iPTMnet'")
    conn.execute("DELETE FROM enzyme_score")
    conn.execute("DELETE FROM protein")
    conn.executemany("INSERT INTO site VALUES (?,?,?,?,?,?,?,?,?)", site_rows)
    conn.executemany("INSERT INTO enzyme_score VALUES (?,?,?,?,?)", score_rows)
    conn.executemany(
        "INSERT OR REPLACE INTO protein VALUES (?,?,?)",
        [(acc, gene, organism) for acc, (gene, organism) in seen_proteins.items()])
    return counts


def _read_quoted_tsv(text: str) -> list[list[str]]:
    """ELM ships quoted, tab-separated tables with # comment headers."""
    import csv
    lines = [line for line in text.splitlines()
             if line and not line.startswith("#")]
    return [row for row in csv.reader(lines, delimiter="\t") if row]


ELM_INSTANCE_CATEGORIES = ["MOD_", "LIG_", "DOC_", "DEG_", "CLV_", "TRG_"]


def _elm_instances_for(query: str) -> list[tuple]:
    """instances.tsv?q=<class-prefix> returns the FULL instance set for that
    category (the bare endpoint serves only a 100-instance sample)."""
    rows = _read_quoted_tsv(
        fetch(f"{ELM_INSTANCES}?q={query}", timeout=120).decode("utf-8", "replace"))
    header_i = next((i for i, row in enumerate(rows)
                     if row and row[0] == "Accession"), None)
    if header_i is None:
        return []
    col = {name: i for i, name in enumerate(rows[header_i])}
    out = []
    for row in rows[header_i + 1:]:
        try:
            out.append((row[col["ELMIdentifier"]], row[col["Primary_Acc"]],
                        int(row[col["Start"]]), int(row[col["End"]])))
        except (KeyError, IndexError, ValueError):
            continue
    return out


def load_elm(conn: sqlite3.Connection) -> dict[str, int]:
    """ELM class index (353 regex classes + probability) + FULL instances.

    Instances download per category prefix (MOD_/LIG_/DOC_/DEG_/CLV_/TRG_,
    ~4.3 k rows total); a flaky category is skipped and counted, not fatal —
    elm.eu.org has uptime windows.
    """
    index_rows = _read_quoted_tsv(fetch(ELM_INDEX, timeout=60).decode("utf-8", "replace"))
    header = next((i for i, row in enumerate(index_rows)
                   if row and row[0] == "Accession"), None)
    if header is None:
        raise ValueError("unexpected ELM index format (no Accession header)")
    classes = []
    for row in index_rows[header + 1:]:
        if len(row) >= 5 and row[1]:
            # ELM classes are motif CONTEXTS, not measured sites: link a type
            # only where the class name/description is unambiguous.
            implies = None
            if row[1].startswith("MOD_NMyristoyl"):
                implies = "Myristoylation"
            elif row[1].startswith("MOD_") and "phosphorylation" in row[3].lower():
                implies = "Phosphorylation"
            probability = float(row[5]) if len(row) > 5 else None
            classes.append((row[1], row[4], row[3] if len(row) > 3 else None,
                            implies, probability))

    instances: list[tuple] = []
    category_counts: dict[str, int] = {}
    for category in ELM_INSTANCE_CATEGORIES:
        try:
            rows = _elm_instances_for(category)
        except Exception:  # noqa: BLE001 - one flaky category is not fatal
            category_counts[f"elm_inst_{category.strip('_').lower()}_failed"] = 1
            continue
        instances.extend(rows)
        category_counts[f"elm_inst_{category.strip('_').lower()}"] = len(rows)

    conn.execute("DELETE FROM motif_class")
    conn.execute("DELETE FROM motif_instance")
    conn.execute("DELETE FROM vocab WHERE source='ELM'")
    for class_id, _, _, implies, _probability in classes:
        if implies:
            conn.execute("INSERT OR REPLACE INTO vocab VALUES (?,?,?)",
                         (class_id, "ELM", implies))
    conn.executemany("INSERT OR REPLACE INTO motif_class VALUES (?,?,?,?,?)", classes)
    conn.executemany("INSERT INTO motif_instance VALUES (?,?,?,?)", instances)
    return {"elm_class": len(classes), "elm_instance": len(instances),
            **category_counts}


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cache", type=Path,
                        default=Path(__file__).resolve().parent.parent / "cache",
                        help="cache directory (default: skill-local cache/, gitignored)")
    parser.add_argument("--types", default="default",
                        help="dbPTM types: 'default', 'all', or comma list")
    parser.add_argument("--skip-elm", action="store_true",
                        help="skip ELM (elm.eu.org has flaky uptime windows)")
    parser.add_argument("--force", action="store_true", help="rebuild the database")
    args = parser.parse_args(argv)

    if args.types == "default":
        types = DEFAULT_TYPES
    elif args.types == "all":
        types = DEFAULT_TYPES + OPTIN_TYPES
    elif args.types == "full":
        types = discover_dbptm_types()
    else:
        types = [t.strip() for t in args.types.split(",") if t.strip()]

    args.cache.mkdir(parents=True, exist_ok=True)
    db_path = args.cache / "ptm_database.db"
    if args.force and db_path.exists():
        db_path.unlink()
    conn = init_db(db_path)
    counts: dict[str, int] = {}
    errors: list[str] = []
    for ptm_type in types:
        try:
            counts[f"dbptm_{ptm_type}"] = load_dbptm_type(conn, ptm_type)
            print(f"dbPTM {ptm_type}: {counts[f'dbptm_{ptm_type}']} sites", file=sys.stderr)
        except Exception as exc:  # noqa: BLE001 - one bad type must not kill the fetch
            errors.append(f"dbPTM {ptm_type}: {exc}")
            print(f"SKIP dbPTM {ptm_type}: {exc}", file=sys.stderr)
    try:
        counts.update(load_iptmnet(conn))
        print(f"iPTMnet: {counts}", file=sys.stderr)
    except Exception as exc:  # noqa: BLE001
        errors.append(f"iPTMnet: {exc}")
        print(f"SKIP iPTMnet: {exc}", file=sys.stderr)
    if not args.skip_elm:
        try:
            counts.update(load_elm(conn))
        except Exception as exc:  # noqa: BLE001
            errors.append(f"ELM: {exc}")
            print(f"SKIP ELM: {exc} (rerun with --skip-elm to silence)", file=sys.stderr)

    conn.executemany("INSERT OR REPLACE INTO meta VALUES (?,?)", [
        ("built", datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")),
        ("sources", json.dumps({"dbPTM": DBPTM, "iPTMnet": IPTMNET,
                                "ELM": ELM_INDEX})),
        ("errors", json.dumps(errors)),
    ])
    conn.commit()
    total = conn.execute("SELECT COUNT(*) FROM site").fetchone()[0]
    conn.close()
    print(json.dumps({"database": str(db_path), "site_rows": total,
                      "counts": counts, "errors": errors}, ensure_ascii=False))
    return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
