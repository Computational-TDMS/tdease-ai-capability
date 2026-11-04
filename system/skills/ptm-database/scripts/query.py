# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""ptm-database query CLI — one atomic interface per database, plus a
one-step flow that chains them.

Atomic commands (each prints the matched records themselves — site rows,
motif classes with their matches, exemplars, top hits — and archives the
full JSON artifact: the CLI return carries the evidence, not just counts):

  local (offline, mirror):            network (live REST):
  query.py dbptm  --accession <ACC>   query.py uniprot --accession <ACC>
            [--type <TYPE>]           query.py uniref  --accession <ACC>
  query.py iptmnet --accession <ACC>  query.py blast   --sequence <SEQ>
  query.py elm    --sequence <SEQ>               [--min-identity 20]
            [--class-prefix <PREFIX>]             [--email you@host]
  query.py cohort --regex <GROUNDED_PATTERN> [--type <TYPE>]

Flow (one step, network + local chained):

  query.py flow --accession <ACC> [--cohort-regex <GROUNDED_PATTERN>]
                [--skip-blast]

  flow queries independent local sources before slow homology searches.
  UniProt supplies sequence when absent. Source failures remain scoped;
  partial reports preserve completed results while remote queries run.

Facts only — the caller reports, the human judges. A database hit is
candidate knowledge, never identity.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sqlite3
import sys
import time
import uuid
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fetch import canonical_type  # noqa: E402 - same skill directory

UNIPROT = "https://rest.uniprot.org"
EBI_BLAST = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
DB_PATH = Path(__file__).resolve().parent.parent / "cache" / "ptm_database.db"

def http_get(url: str, timeout: int = 60) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "tdease-ptm-database/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def open_mirror(path: Path | None) -> sqlite3.Connection:
    db = path or DB_PATH
    if not db.exists():
        raise SystemExit(f"local mirror not found at {db} — "
                         "build it first: uv run --script scripts/fetch.py")
    return sqlite3.connect(f"file:{db}?mode=ro", uri=True)


def parse_sequence(sequence: str | None, fasta: Path | None) -> str:
    sequence = (sequence or "").upper().replace(" ", "").replace("\n", "")
    if not sequence and fasta:
        text = fasta.read_text(encoding="utf-8")
        if sum(line.startswith(">") for line in text.splitlines()) > 1:
            raise SystemExit("--fasta needs one sequence; query multiple entries separately")
        sequence = "".join(line.strip() for line in text
                           if not line.startswith(">")).upper()
    bad = set(sequence) - set("ACDEFGHIKLMNPQRSTVWY")
    if bad:
        raise SystemExit(f"invalid residues in sequence: {sorted(bad)}")
    return sequence


# ── mirror lookups (offline) ──
def resolve_type(conn: sqlite3.Connection, user: str) -> str:
    row = conn.execute("SELECT canonical FROM vocab WHERE alias=? LIMIT 1",
                       (user.strip(),)).fetchone()
    return row[0] if row else canonical_type(user)


def type_clause(conn: sqlite3.Connection, types: list[str] | None):
    if not types:
        return "", []
    canonical = list(dict.fromkeys(resolve_type(conn, t) for t in types))
    placeholders = ",".join("?" for _ in canonical)
    return f" AND ptm_type IN ({placeholders})", canonical


def mirror_sites(conn: sqlite3.Connection, accession: str, source: str,
                 types: list[str] | None) -> dict:
    where, params = type_clause(conn, types)
    rows = conn.execute(
        "SELECT ptm_type, position, residue, pmids, seq_window FROM site "
        f"WHERE accession=? AND source=?{where} ORDER BY ptm_type, position",
        [accession, source, *params]).fetchall()
    by_type: dict[str, int] = {}
    for row in rows:
        by_type[row[0]] = by_type.get(row[0], 0) + 1
    return {"source": source, "accession": accession, "site_count": len(rows),
            "by_type": by_type,
            "sites": [list(r) for r in rows]}


def iptmnet_scores(conn: sqlite3.Connection, accession: str,
                   types: list[str] | None) -> list[list]:
    where, params = type_clause(conn, types)
    return [list(r) for r in conn.execute(
        "SELECT ptm_type, position, enzyme, score FROM enzyme_score "
        f"WHERE accession=?{where} ORDER BY score DESC",
        [accession, *params]).fetchall()]


def cohort_scan(conn: sqlite3.Connection, regexes: list[str],
                types: list[str] | None, max_exemplars: int) -> list[dict]:
    where, params = type_clause(conn, types)
    out = []
    for pattern in regexes:
        try:
            regex = re.compile(pattern)
        except re.error as exc:
            raise SystemExit(f"--cohort-regex {pattern!r} is not a valid regex: {exc}")
        total = 0
        by_type: dict[str, int] = {}
        exemplars: list[list] = []
        for accession, entry, position, ptm_type, window in conn.execute(
                "SELECT accession, entry, position, ptm_type, seq_window FROM site "
                f"WHERE seq_window IS NOT NULL{where}", params):
            if regex.search(window):
                total += 1
                by_type[ptm_type] = by_type.get(ptm_type, 0) + 1
                if len(exemplars) < max_exemplars:
                    exemplars.append([accession, entry, position, ptm_type, window])
        out.append({"regex": pattern, "total": total, "by_type": by_type,
                    "exemplars": exemplars})
    return out


def elm_class_scan(conn: sqlite3.Connection, sequence: str,
                   class_prefix: str | None) -> dict:
    """ELM web-search semantics: every class regex is tried at EVERY start
    position (overlapping matches), all matches kept — reproduces the
    website's sequence-search counts (its species/SMART/structural filters
    are a separate annotation layer we do not mirror)."""
    classes_out = []
    match_count = 0
    skipped = 0
    for class_id, pattern, description, implies, probability in conn.execute(
            "SELECT class_id, pattern, description, implies_type, probability "
            "FROM motif_class"):
        if class_prefix and not class_id.startswith(class_prefix):
            continue
        try:
            regex = re.compile(pattern)
        except re.error:
            skipped += 1
            continue
        matches = []
        for start in range(len(sequence)):
            found = regex.match(sequence, start)
            if found:
                matches.append({"start": start + 1, "end": start + len(found.group(0)),
                                "match": found.group(0)})
        if matches:
            match_count += len(matches)
            classes_out.append({"class": class_id, "description": description,
                                "implies_type": implies, "probability": probability,
                                "matches": matches})
    classes_out.sort(key=lambda c: c["matches"][0]["start"])
    return {"class_count": len(classes_out), "match_count": match_count,
            "skipped_uncompilable": skipped, "classes": classes_out}


def elm_instances(conn: sqlite3.Connection, accession: str) -> list[list]:
    return [list(r) for r in conn.execute(
        "SELECT class_id, start, end FROM motif_instance WHERE accession=?",
        (accession,))]


# ── network lookups (UniProt / EBI) ──
def uniprot_entry(accession: str) -> dict:
    entry = json.loads(http_get(f"{UNIPROT}/uniprotkb/{accession}.json"))
    sequence = entry.get("sequence", {}).get("value", "").upper()
    mod_res = [{"position": f["location"]["start"]["value"],
                "description": f["description"]}
               for f in entry.get("features", [])
               if f["type"] == "Modified residue"]
    cc_ptm = [c.get("texts", []) for c in entry.get("comments", [])
              if c.get("commentType") == "PTM"]
    desc = entry.get("proteinDescription", {})
    return {"accession": entry.get("primaryAccession", accession),
            "name": (desc.get("recommendedName", {}).get("fullName", {}).get("value")
                     or (desc.get("submissionNames") or [{}])[0]
                     .get("fullName", {}).get("value")),
            "organism": entry.get("organism", {}).get("scientificName"),
            "sequence": sequence, "length": len(sequence),
            "mod_res": mod_res, "cc_ptm": [t for texts in cc_ptm for t in texts],
            "features": entry.get("features", []), "comments": entry.get("comments", [])}


def uniref50_members(accession: str) -> list[str]:
    # The entry endpoint /uniref/UniRef50_<acc>.json only exists when the
    # accession IS the cluster representative; a plain member (the common
    # case — e.g. an unreviewed entry inside a reviewed cluster) 404s there.
    # Search by accession first, pick the UniRef50 hit, then read its entry.
    search = json.loads(http_get(
        f"{UNIPROT}/uniref/search?query={accession}&format=json&size=10"
    ))
    cluster_id = next(
        (r["id"] for r in search.get("results", []) if r.get("entryType") == "UniRef50"),
        f"UniRef50_{accession}",  # representative form: fall back to the entry
    )
    data = json.loads(http_get(f"{UNIPROT}/uniref/{cluster_id}.json?size=500"))
    result = data.get("results", [data])[0]
    members = []
    for member in result.get("members", []):
        # The API labels UniProt entries "UniProtKB ID" (mnemonic form,
        # e.g. Q2KMJ4_SOYBN); UniParc members carry no accession to cite.
        if member.get("memberIdType") not in ("UniProtKB", "UniProtKB ID"):
            continue
        member_id = member.get("memberId", "")
        member_id = member_id.split("_")[0]
        if member_id and member_id != accession:
            members.append(member_id)
    return members


def run_blast(sequence: str, email: str, timeout_s: int, database: str,
              min_identity: float) -> list[dict]:
    params = urllib.parse.urlencode({
        "email": email, "program": "blastp", "stype": "protein",
        "database": database, "matrix": "BLOSUM62", "exp": "1e-3",
        "alignments": "100", "scores": "100", "gapalign": "true",
        "sequence": sequence,
    }).encode()
    request = urllib.request.Request(f"{EBI_BLAST}/run", data=params,
                                     headers={"User-Agent": "tdease-ptm-database/1.0"})
    job = urllib.request.urlopen(request, timeout=60).read().decode().strip()
    deadline = time.monotonic() + timeout_s
    status = "RUNNING"
    while time.monotonic() < deadline:
        status = http_get(f"{EBI_BLAST}/status/{job}", timeout=30).decode().strip()
        if status in ("FINISHED", "ERROR", "FAILURE", "NOT_FOUND"):
            break
        time.sleep(5)
    if status != "FINISHED":
        raise RuntimeError(f"BLAST job {job} ended as {status}")
    report = http_get(f"{EBI_BLAST}/result/{job}/out", timeout=120).decode("utf-8", "replace")
    hits: list[dict] = []
    for block in re.split(r"^>", report, flags=re.MULTILINE)[1:]:
        header = block.splitlines()[0]
        match = re.match(r"(?:SP|TR):([A-Z0-9]+)\s+(\S+)", header)
        if not match:
            continue
        ident = re.search(r"Identities\s*=\s*\d+/\d+\s*\((\d+)%\)", block)
        evalue = re.search(r"Expect\s*=\s*([0-9.eE+-]+)", block)
        identity = float(ident.group(1)) if ident else 0.0
        if identity < min_identity:
            continue
        hits.append({"accession": match.group(1), "entry": match.group(2),
                     "identity_pct": identity,
                     "evalue": evalue.group(1) if evalue else None})
    hits.sort(key=lambda h: -h["identity_pct"])
    return hits


def lookup_brief(conn: sqlite3.Connection, accession: str) -> dict:
    rows = conn.execute("SELECT ptm_type, position FROM site WHERE accession=?",
                        (accession,)).fetchall()
    by_type: dict[str, int] = {}
    for row in rows:
        by_type[row[0]] = by_type.get(row[0], 0) + 1
    return {"accession": accession, "site_count": len(rows), "by_type": by_type}


# ── output contract: evidence records on stdout, full report on disk ──
def emit(command: str, name: str, payload: dict, summary: dict,
         out: Path | None, full_stdout: bool) -> int:
    query = summary.pop("_query", {})
    report = {"command": command, "query": query, **payload}
    out_path = out or Path(f"ptmdb_{command}_{name}_{uuid.uuid4().hex[:12]}.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, ensure_ascii=False), encoding="utf-8")
    summary = {"command": command, "report": str(out_path.resolve()),
               "query": query, **summary}
    print(json.dumps(report if full_stdout else summary, ensure_ascii=False))
    return 0


def mirror_meta(conn: sqlite3.Connection) -> dict:
    built = conn.execute("SELECT value FROM meta WHERE key='built'").fetchone()
    total = conn.execute("SELECT COUNT(*) FROM site").fetchone()[0]
    return {"site_rows_total": total, "built": built[0] if built else None}


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--out", type=Path, help="artifact path override")
    common.add_argument("--full-stdout", action="store_true",
                        help="print the full JSON instead of the compact summary")
    sub = parser.add_subparsers(dest="command", required=True)

    def add_local(name: str, help_text: str, *, sequence_input=False):
        p = sub.add_parser(name, help=help_text, parents=[common])
        p.add_argument("--db", type=Path, help="mirror path override")
        p.add_argument("--accession")
        if sequence_input:
            source = p.add_mutually_exclusive_group()
            source.add_argument("--sequence")
            source.add_argument("--fasta", type=Path)
        else:
            p.add_argument("--type", action="append", default=[],
                           help="PTM type filter; accepts source vocabulary aliases")
        return p

    add_local("dbptm", "dbPTM rows for one accession (mirror, offline)")
    add_local("iptmnet", "iPTMnet sites + enzyme scores for one accession (mirror)")
    p_elm = add_local("elm", "ELM class regexes on a sequence (mirror); "
                             "--accession adds known instances", sequence_input=True)
    p_elm.add_argument("--class-prefix", default=None, help="e.g. MOD_ or LIG_")
    p_cohort = sub.add_parser("cohort", help="regex scan over every mirrored "
                                             "seq_window (the query no database "
                                             "website exposes)", parents=[common])
    p_cohort.add_argument("--regex", action="append", required=True)
    p_cohort.add_argument("--db", type=Path, help="mirror path override")
    p_cohort.add_argument("--type", action="append", default=[])
    p_cohort.add_argument("--max-exemplars", type=int, default=20)

    p_uni = sub.add_parser("uniprot", help="UniProt entry: sequence + own PTM "
                                           "annotations (network)", parents=[common])
    p_uni.add_argument("--accession", required=True)
    p_ref = sub.add_parser("uniref", help="UniRef50 cluster members (network)", parents=[common])
    p_ref.add_argument("--accession", required=True)

    p_blast = sub.add_parser("blast", help="BLAST vs reviewed Swiss-Prot via "
                                           "EBI (network)", parents=[common])
    p_blast.add_argument("--sequence", required=True)
    p_blast.add_argument("--blast-db", default="uniprotkb_swissprot")
    p_blast.add_argument("--email", default=os.environ.get(
        "TDEASE_BLAST_EMAIL", "tdease.database@gmail.com"),
        help="EBI requires a contact email; set TDEASE_BLAST_EMAIL ideally")
    p_blast.add_argument("--blast-timeout", type=int, default=300)
    p_blast.add_argument("--min-identity", type=float, default=20.0)

    p_flow = sub.add_parser("flow", help="one step: uniprot + uniref + blast + "
                                         "mirror three levels + cohort", parents=[common])
    p_flow.add_argument("--accession", required=True)
    p_flow.add_argument("--db", type=Path, help="mirror path override")
    p_flow.add_argument("--sequence", help="override the UniProt sequence")
    p_flow.add_argument("--cohort-regex", action="append", default=[])
    p_flow.add_argument("--type", action="append", default=[])
    p_flow.add_argument("--skip-blast", action="store_true",
                        help="skip EBI BLAST; UniProt and UniRef still use the network")
    p_flow.add_argument("--email", default=os.environ.get(
        "TDEASE_BLAST_EMAIL", "tdease.database@gmail.com"))
    p_flow.add_argument("--blast-timeout", type=int, default=300)
    p_flow.add_argument("--min-identity", type=float, default=20.0)
    p_flow.add_argument("--max-exemplars", type=int, default=20)

    args = parser.parse_args(argv)
    conn_holder: sqlite3.Connection | None = None

    def mirror() -> sqlite3.Connection:
        nonlocal conn_holder
        if conn_holder is None:
            conn_holder = open_mirror(args.db)
        return conn_holder

    if args.command in ("dbptm", "iptmnet"):
        if not args.accession:
            raise SystemExit(f"{args.command} needs --accession")
        payload = mirror_sites(mirror(), args.accession,
                               "dbPTM" if args.command == "dbptm" else "iPTMnet",
                               args.type)
        if args.command == "iptmnet":
            payload["enzyme_scores"] = iptmnet_scores(mirror(), args.accession, args.type)
        payload["meta"] = mirror_meta(mirror())
        summary = {"_query": {"accession": args.accession, "types": args.type,
                              "mirror": str((args.db or DB_PATH).resolve())},
                   **{k: payload[k] for k in ("site_count", "by_type")},
                   "sites": payload["sites"][:20],
                   "sites_truncated": payload["site_count"] > 20}
        if args.command == "iptmnet":
            summary["enzyme_scores"] = payload["enzyme_scores"][:10]
        return emit(args.command, args.accession, payload, summary,
                    args.out, args.full_stdout)

    if args.command == "elm":
        sequence = parse_sequence(args.sequence, args.fasta)
        if not sequence and not args.accession:
            raise SystemExit("elm needs --sequence/--fasta and/or --accession")
        payload: dict = {"meta": mirror_meta(mirror())}
        if sequence:
            payload["class_scan"] = elm_class_scan(mirror(), sequence, args.class_prefix)
        if args.accession:
            payload["instances"] = elm_instances(mirror(), args.accession)
        scan = payload.get("class_scan", {})
        summary = {"_query": {"accession": args.accession,
                              "class_prefix": args.class_prefix, "sequence": sequence,
                              "sequence_source": str(args.fasta.resolve()) if args.fasta else "inline",
                              "coordinate_basis": "input_sequence_1_based",
                              "mirror": str((args.db or DB_PATH).resolve())},
                   "class_count": scan.get("class_count"),
                   "match_count": scan.get("match_count"),
                   "classes": scan.get("classes", [])[:30],
                   "classes_truncated": scan.get("class_count", 0) > 30,
                   "skipped_uncompilable": scan.get("skipped_uncompilable", 0),
                   "instances": len(payload.get("instances", [])),
                   "instance_records": payload.get("instances", [])[:20],
                   "instances_truncated": len(payload.get("instances", [])) > 20}
        return emit("elm", args.accession or sequence[:8], payload, summary,
                    args.out, args.full_stdout)

    if args.command == "cohort":
        payload = {"scans": cohort_scan(mirror(), args.regex, args.type,
                                        args.max_exemplars),
                   "meta": mirror_meta(mirror())}
        summary = {"_query": {"regexes": args.regex, "types": args.type,
                              "max_exemplars": args.max_exemplars,
                              "mirror": str((args.db or DB_PATH).resolve())},
                   "scans": [{"regex": s["regex"], "total": s["total"],
                              "by_type": s["by_type"],
                              "exemplars": s["exemplars"][:5]}
                             for s in payload["scans"]]}
        return emit("cohort", "scan", payload, summary, args.out, args.full_stdout)

    if args.command == "uniprot":
        entry = uniprot_entry(args.accession)
        return emit("uniprot", args.accession,
                    {"entry": entry | {"sequence_length": entry["length"],
                        "n_terminal_20": entry["sequence"][:20]}},
                    {"_query": {"accession": args.accession},
                     "name": entry["name"], "organism": entry["organism"],
                     "length": entry["length"],
                     "mod_res_count": len(entry["mod_res"]),
                     "mod_res": entry["mod_res"][:20],
                     "feature_count": len(entry["features"]),
                     "features": entry["features"][:20],
                     "features_truncated": len(entry["features"]) > 20,
                     "cc_ptm": entry["cc_ptm"][:3]},
                    args.out, args.full_stdout)

    if args.command == "uniref":
        members = uniref50_members(args.accession)
        return emit("uniref", args.accession, {"members": members},
                    {"_query": {"accession": args.accession},
                     "member_count": len(members),
                     "members": members[:20]},
                    args.out, args.full_stdout)

    if args.command == "blast":
        sequence = parse_sequence(args.sequence, None)
        hits = run_blast(sequence, args.email, args.blast_timeout,
                         args.blast_db, args.min_identity)
        return emit("blast", sequence[:8], {"hits": hits},
                    {"_query": {"database": args.blast_db,
                                "min_identity": args.min_identity, "sequence": sequence,
                                "timeout_seconds": args.blast_timeout},
                     "hit_count": len(hits),
                     "top": [{k: h[k] for k in ("accession", "identity_pct", "evalue")}
                             for h in hits[:20]]},
                    args.out, args.full_stdout)

    if args.command == "flow":
        sequence = parse_sequence(args.sequence, None)
        query = {"accession": args.accession, "sequence": sequence or None,
                 "cohort_regexes": args.cohort_regex, "types": args.type,
                 "mirror": str((args.db or DB_PATH).resolve()),
                 "skip_blast": args.skip_blast, "min_identity": args.min_identity,
                 "max_exemplars": args.max_exemplars, "timeout_seconds": args.blast_timeout}
        if args.out is None:
            args.out = Path(f"ptmdb_flow_{args.accession}_{uuid.uuid4().hex[:12]}.json")
        for pattern in args.cohort_regex:
            try:
                re.compile(pattern)
            except re.error as exc:
                raise SystemExit(f"--cohort-regex {pattern!r} is not a valid regex: {exc}")
        source_status = {name: "pending" for name in (
            "uniprot", "dbptm", "iptmnet", "iptmnet_scores", "elm", "cohort",
            "uniref50", "blast")}
        if args.skip_blast:
            source_status["blast"] = "skipped"
        lane_errors: dict[str, str] = {}
        try:
            conn = mirror()
        except (SystemExit, OSError, sqlite3.Error) as exc:
            conn = None
            lane_errors["mirror"] = str(exc)

        def local():
            if conn is None:
                raise RuntimeError(lane_errors["mirror"])
            return conn

        def read_source(name, operation, fallback, has_hits=bool):
            try:
                value = operation()
                source_status[name] = "hit" if has_hits(value) else "empty"
                return value
            except Exception as exc:
                source_status[name] = "unavailable"
                lane_errors[name] = f"{type(exc).__name__}: {exc}"
                return fallback

        empty_sites = {"site_count": 0, "by_type": {}, "sites": []}
        empty_elm = {"class_count": 0, "match_count": 0,
                     "skipped_uncompilable": 0, "classes": []}
        entry = {"name": None, "organism": None, "sequence": "",
                 "mod_res": [], "cc_ptm": []}
        levels: dict[str, object] = {
            "uniprot": {"mod_res": [], "cc_ptm": []},
            "direct": read_source("dbptm", lambda: mirror_sites(
                local(), args.accession, "dbPTM", args.type), dict(empty_sites),
                lambda value: value["site_count"] > 0) | {
                    "iptmnet": read_source("iptmnet", lambda: mirror_sites(
                        local(), args.accession, "iPTMnet", args.type), dict(empty_sites),
                        lambda value: value["site_count"] > 0),
                    "enzyme_scores": read_source("iptmnet_scores", lambda: iptmnet_scores(
                        local(), args.accession, args.type), []),
                },
            "homolog": {"uniref50": [], "blast": None},
            "cohort": {"cohort_regexes": [], "elm_classes": dict(empty_elm)},
        }
        if args.cohort_regex:
            levels["cohort"]["cohort_regexes"] = read_source(
                "cohort", lambda: cohort_scan(local(), args.cohort_regex, args.type,
                                              args.max_exemplars), [],
                lambda value: any(row["total"] for row in value))
        else:
            source_status["cohort"] = "skipped"
        if sequence:
            levels["cohort"]["elm_classes"] = read_source(
                "elm", lambda: elm_class_scan(local(), sequence, None), dict(empty_elm),
                lambda value: value["match_count"] > 0)
        meta = {}
        if conn is not None:
            try:
                meta = mirror_meta(conn)
            except sqlite3.Error as exc:
                lane_errors["mirror_meta"] = str(exc)

        def group_status(*names):
            states = {source_status[name] for name in names} - {"skipped"}
            if "pending" in states:
                return "pending"
            if "unavailable" in states:
                return "unavailable" if states == {"unavailable"} else "partial"
            return "hit" if "hit" in states else ("empty" if states else "skipped")

        def progress(complete=False):
            lane_status = {
                "uniprot": source_status["uniprot"],
                "direct": group_status("dbptm", "iptmnet", "iptmnet_scores"),
                "motif": group_status("elm", "cohort"),
                "uniref50": source_status["uniref50"], "blast": source_status["blast"],
            }
            payload = {"complete": complete, "query": query,
                       "input": {"accession": args.accession, "name": entry["name"],
                                 "organism": entry["organism"], "sequence": sequence,
                                 "sequence_source": "supplied" if args.sequence else (
                                     "uniprot" if sequence else None),
                                 "sequence_length": len(sequence), "n_terminal_20": sequence[:20]},
                       "levels": levels, "meta": meta, "lane_status": lane_status,
                       "source_status": source_status, "lane_errors": lane_errors}
            # Same report as the final return, so a slow/failed remote leg
            # cannot erase completed local evidence or its unfinished scope.
            out_path = args.out or Path(f"ptmdb_flow_{args.accession}.json")
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text(json.dumps({"command": "flow", **payload}, ensure_ascii=False),
                                encoding="utf-8")
            return payload

        progress()
        entry = read_source("uniprot", lambda: uniprot_entry(args.accession), entry)
        levels["uniprot"] = {"mod_res": entry["mod_res"], "cc_ptm": entry["cc_ptm"],
                             "features": entry.get("features", []),
                             "comments": entry.get("comments", [])}
        sequence = sequence or parse_sequence(entry["sequence"], None)
        if source_status["elm"] == "pending":
            if sequence:
                levels["cohort"]["elm_classes"] = read_source(
                    "elm", lambda: elm_class_scan(local(), sequence, None), dict(empty_elm),
                    lambda value: value["match_count"] > 0)
            else:
                source_status["elm"] = "unavailable"
                lane_errors["elm"] = "No supplied or retrieved sequence for motif search"
        progress()
        try:
            members = uniref50_members(args.accession)
            levels["homolog"]["uniref50_members"] = members
            source_status["uniref50"] = "hit" if members else "empty"
            levels["homolog"]["uniref50_annotation_status"] = "unavailable"
            if conn is not None:
                for member in members:
                    brief = lookup_brief(conn, member)
                    if brief["site_count"]:
                        levels["homolog"]["uniref50"].append(brief)
                levels["homolog"]["uniref50_annotation_status"] = "complete"
        except Exception as exc:
            if source_status["uniref50"] == "pending":
                source_status["uniref50"] = "unavailable"
            lane_errors["uniref50"] = f"{type(exc).__name__}: {exc}"
        progress()
        if args.skip_blast:
            source_status["blast"] = "skipped"
        elif not sequence:
            source_status["blast"] = "unavailable"
            lane_errors["blast"] = "No supplied or retrieved sequence for BLAST"
        else:
            try:
                blast_level = {"hits_with_mirror_rows": [], "hits_without_rows": 0,
                               "annotation_status": "unavailable"}
                hits = run_blast(sequence, args.email, args.blast_timeout,
                                 "uniprotkb_swissprot", args.min_identity)
                blast_level["hits"] = hits
                levels["homolog"]["blast"] = blast_level
                source_status["blast"] = "hit" if hits else "empty"
                if conn is not None:
                    for hit in hits:
                        brief = lookup_brief(conn, hit["accession"])
                        if brief["site_count"]:
                            blast_level["hits_with_mirror_rows"].append({**hit, **brief})
                        else:
                            blast_level["hits_without_rows"] += 1
                    blast_level["annotation_status"] = "complete"
                else:
                    blast_level["hits_without_rows"] = None
            except Exception as exc:
                if source_status["blast"] == "pending":
                    source_status["blast"] = "unavailable"
                lane_errors["blast"] = f"{type(exc).__name__}: {exc}"
        payload = progress(complete=True)
        lane_status = payload["lane_status"]
        blast_info = levels["homolog"]["blast"]  # type: ignore[assignment]
        direct: dict = levels["direct"]  # type: ignore[assignment]
        summary = {"_query": {"accession": args.accession,
                              "cohort_regex": args.cohort_regex,
                              "skip_blast": args.skip_blast, "types": args.type},
                   "uniprot_mod_res": len(entry["mod_res"]),
                   "uniprot_mod_res_rows": entry["mod_res"][:20],
                   "uniprot_feature_count": len(entry.get("features", [])),
                   "uniprot_features": entry.get("features", [])[:20],
                   "uniprot_features_truncated": len(entry.get("features", [])) > 20,
                   "direct_sites": direct["site_count"],
                   "direct_dbptm_sites": direct["sites"][:20],
                   "direct_iptmnet_sites": direct["iptmnet"]["sites"][:20],
                   "direct_iptmnet_site_count": direct["iptmnet"]["site_count"],
                   "iptmnet_enzyme_scores": direct["enzyme_scores"][:10],
                   "uniref50_with_rows": len(levels["homolog"]["uniref50"]),  # type: ignore[index]
                   "uniref50_row_briefs": levels["homolog"]["uniref50"][:10],  # type: ignore[index]
                   "blast_hits_with_rows": (len(blast_info["hits_with_mirror_rows"])
                                            if blast_info else None),
                   "blast_top_with_rows": (blast_info["hits_with_mirror_rows"][:10]
                                           if blast_info else None),
                   "cohort_regexes": [{"regex": c["regex"], "total": c["total"]}
                                      for c in levels["cohort"]["cohort_regexes"]],  # type: ignore[index]
                   "cohort_exemplars": [{"regex": c["regex"],
                                         "exemplars": c["exemplars"][:5]}
                                        for c in levels["cohort"]["cohort_regexes"]],  # type: ignore[index]
                   "elm_class_hits": levels["cohort"]["elm_classes"]["class_count"],  # type: ignore[index]
                   "elm_match_total": levels["cohort"]["elm_classes"]["match_count"],  # type: ignore[index]
                   "elm_class_records": levels["cohort"]["elm_classes"]["classes"][:30],  # type: ignore[index]
                   "elm_classes_truncated": levels["cohort"]["elm_classes"]["class_count"] > 30,  # type: ignore[index]
                   "elm_skipped_uncompilable": levels["cohort"]["elm_classes"]["skipped_uncompilable"],  # type: ignore[index]
                   "lane_status": lane_status,
                   "source_status": source_status,
                   "complete": True,
                   "lane_errors": lane_errors,
                   **meta}
        return emit("flow", args.accession, payload, summary, args.out,
                    args.full_stdout)

    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
