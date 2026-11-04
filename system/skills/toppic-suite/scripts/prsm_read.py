# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Read a TopPIC PrSM or proteoform table without loading it into a generic table layer.

TopPIC writes a parameter preamble above the real header, so a naive TSV reader
either fails or silently treats the preamble as data. This script finds the
header line, then answers three domain questions:

  schema   which columns exist, how many data rows, which optional fields are present
  summary  identification counts, score distributions, modification prevalence
  rows     filtered rows, using domain filters only -- no query language

Prints one JSON object on stdout. Read-only unless --out is given.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _shared.toppic import parse_mods

# Columns whose presence identifies a TopPIC table. Names have drifted across
# TopPIC versions, so identification is by any-of, not by an exact set.
SIGNATURE_COLUMNS = ("Prsm ID", "Spectrum ID", "Protein accession", "Proteoform", "Proteoform ID")

CANONICAL = {
    "prsm_id": ("Prsm ID", "PrSM ID"),
    "spectrum_id": ("Spectrum ID",),
    "scan": ("Scan(s)", "Scans", "Scan"),
    "retention_time": ("Retention time", "RT"),
    "charge": ("Charge",),
    "precursor_mass": ("Precursor mass",),
    "adjusted_precursor_mass": ("Adjusted precursor mass",),
    "proteoform_id": ("Proteoform ID",),
    "feature_id": ("Feature ID",),
    "feature_intensity": ("Feature intensity",),
    "protein_hits": ("#Protein hits", "Protein hits"),
    "accession": ("Protein accession",),
    "description": ("Protein description",),
    "first_residue": ("First residue",),
    "last_residue": ("Last residue",),
    "proteoform": ("Proteoform",),
    "proteoform_mass": ("Proteoform mass",),
    "n_terminal_form": ("Protein N-terminal form",),
    "fixed_ptms": ("Fixed PTMs",),
    "unexpected_count": ("#unexpected modifications", "#Unexpected modifications"),
    "unexpected_mods": ("unexpected modifications", "Unexpected modifications"),
    "variable_count": ("#variable PTMs",),
    "variable_ptms": ("variable PTMs", "Variable PTMs"),
    "mi_score": ("MIScore", "MIScore(s)"),
    "matched_peaks": ("#matched peaks", "#Matched peaks"),
    "matched_fragments": ("#matched fragment ions", "#Matched fragment ions"),
    "e_value": ("E-value",),
    "q_value": ("Spectrum-level Q-value", "Q-value"),
    "proteoform_q_value": ("Proteoform-level Q-value",),
}


class ToolError(Exception):
    def __init__(self, code: str, message: str, **detail):
        super().__init__(message)
        self.code = code
        self.message = message
        self.detail = detail


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_table(path: Path) -> tuple[list[str], list[dict], int]:
    """Return (columns, rows, header_line_number) skipping the TopPIC preamble."""
    if not path.is_file():
        raise ToolError("FILE_NOT_FOUND", f"{path} does not exist")
    try:
        text = path.read_text(encoding="utf-8-sig", errors="replace")
    except OSError as exc:
        raise ToolError("FILE_UNREADABLE", str(exc)) from exc

    lines = text.splitlines()
    header_index = None
    for index, line in enumerate(lines):
        if not line.strip() or line.lstrip().startswith(("*", "#!")):
            continue
        if "\t" not in line:
            continue
        cells = [cell.strip() for cell in line.split("\t")]
        if any(name in cells for name in SIGNATURE_COLUMNS):
            header_index = index
            break
    if header_index is None:
        raise ToolError(
            "NOT_A_TOPPIC_TABLE",
            "no TopPIC header row found",
            path=str(path),
            hint=f"expected a tab-separated line containing one of {list(SIGNATURE_COLUMNS)}",
        )

    reader = csv.DictReader(lines[header_index:], delimiter="\t")
    columns = [name.strip() for name in (reader.fieldnames or [])]
    rows = []
    for raw in reader:
        row = {}
        for key, value in raw.items():
            if key is None:
                continue
            row[key.strip()] = (value or "").strip() if isinstance(value, str) else value
        # A trailing blank line yields a row of empty strings; drop it.
        if any(row.values()):
            rows.append(row)
    return columns, rows, header_index + 1


def resolve_columns(columns: list[str]) -> dict:
    present = {}
    lowered = {name.lower(): name for name in columns}
    for key, candidates in CANONICAL.items():
        for candidate in candidates:
            actual = lowered.get(candidate.lower())
            if actual:
                present[key] = actual
                break
    return present


def cell(row: dict, mapping: dict, key: str) -> str:
    column = mapping.get(key)
    return row.get(column, "") if column else ""


def as_float(value: str) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(result) else result


def as_int(value: str) -> int | None:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


MASS_SHIFT_PATTERN = re.compile(r"\[([+-]?\d+(?:\.\d+)?)\]")


def unexpected_shifts(row: dict, mapping: dict) -> list:
    """Numeric unexpected mass shifts for one row.

    TopPIC reports them in a dedicated column in newer versions and only inside
    the proteoform string in older ones, so read the column when it exists and
    fall back to bracketed numbers in the proteoform string.
    """
    text = cell(row, mapping, "unexpected_mods")
    count = as_int(cell(row, mapping, "unexpected_count"))
    if "unexpected_mods" in mapping:
        try:
            shifts = [item[0] for item in parse_mods(text, strict=True)]
            if count is not None and count != len(shifts):
                raise ValueError(f"declared count {count} differs from parsed count {len(shifts)}")
            return shifts
        except ValueError as exc:
            raise ToolError("UNPARSEABLE_MODIFICATIONS", str(exc),
                            prsm_id=cell(row, mapping, "prsm_id"), raw=text) from exc
    if count == 0:
        return []
    variable_count = as_int(cell(row, mapping, "variable_count"))
    if variable_count:
        # Variable PTMs and unexpected shifts render identically as bracketed
        # masses in the proteoform string. With both present the regex cannot
        # attribute values: a matching total can still be the wrong set.
        raise ToolError(
            "UNPARSEABLE_MODIFICATIONS",
            "legacy table has no unexpected-modifications column and the proteoform "
            f"string mixes {variable_count} variable PTM(s) with unexpected shifts; "
            "attribution is ambiguous — supply a table with the dedicated column",
            prsm_id=cell(row, mapping, "prsm_id"))
    shifts = [float(match) for match in MASS_SHIFT_PATTERN.findall(
        cell(row, mapping, "proteoform"))]
    if count is not None and count != len(shifts):
        raise ToolError("UNPARSEABLE_MODIFICATIONS", "numeric shifts do not match declared count",
                        prsm_id=cell(row, mapping, "prsm_id"))
    return shifts


def quantiles(values: list[float]) -> dict:
    if not values:
        return {}
    ordered = sorted(values)

    def pick(fraction: float) -> float:
        index = min(len(ordered) - 1, max(0, round(fraction * (len(ordered) - 1))))
        return ordered[index]

    return {
        "min": ordered[0],
        "p25": pick(0.25),
        "median": pick(0.5),
        "p75": pick(0.75),
        "max": ordered[-1],
    }


def cmd_schema(args, columns, rows, header_line, mapping) -> dict:
    return {
        "query": {"kind": "schema", "path": args.path},
        "header_line": header_line,
        "row_count": len(rows),
        "columns": columns,
        "recognized": mapping,
        "unrecognized": [name for name in columns if name not in set(mapping.values())],
        "missing_expected": [key for key in ("accession", "proteoform", "e_value") if key not in mapping],
    }


def cmd_summary(args, columns, rows, header_line, mapping) -> dict:
    accessions = set()
    proteoforms = set()
    e_values = []
    q_values = []
    masses = []
    # One column basis for the whole distribution: mixing proteoform and
    # precursor mass value-by-value (and `or` swallows a legitimate 0.0)
    # would quietly blend two different measurements. Older tables without
    # the proteoform-mass column fall back to precursor mass, declared.
    mass_column = "proteoform_mass" if "proteoform_mass" in mapping else "precursor_mass"
    charges = []
    with_unexpected = 0
    shift_total = 0
    multi_hit = 0

    for row in rows:
        accession = cell(row, mapping, "accession")
        if accession:
            accessions.add(accession)
        proteoform_id = cell(row, mapping, "proteoform_id") or cell(row, mapping, "proteoform")
        if proteoform_id:
            proteoforms.add(proteoform_id)
        value = as_float(cell(row, mapping, "e_value"))
        if value is not None:
            e_values.append(value)
        value = as_float(cell(row, mapping, "q_value"))
        if value is not None:
            q_values.append(value)
        value = as_float(cell(row, mapping, mass_column))
        if value is not None:
            masses.append(value)
        value = as_int(cell(row, mapping, "charge"))
        if value is not None:
            charges.append(float(value))
        hits = as_int(cell(row, mapping, "protein_hits"))
        if hits is not None and hits > 1:
            multi_hit += 1
        shifts = unexpected_shifts(row, mapping)
        if shifts:
            with_unexpected += 1
            shift_total += len(shifts)

    return {
        "query": {"kind": "summary", "path": args.path},
        "row_count": len(rows),
        "unique_proteins": len(accessions),
        "unique_proteoforms": len(proteoforms),
        "rows_with_unexpected_shift": with_unexpected,
        "unexpected_shift_total": shift_total,
        "rows_with_multiple_protein_hits": multi_hit,
        "e_value": quantiles(e_values),
        "q_value": quantiles(q_values),
        "mass": {**quantiles(masses), "column": mass_column},
        "charge": quantiles(charges),
        "note": "E-value and Q-value describe spectrum-match reliability only. They say nothing about modification identity, localization or biological reality.",
    }


def cmd_rows(args, columns, rows, header_line, mapping) -> dict:
    selected = []
    for row in rows:
        if args.feature_id and cell(row, mapping, "feature_id") != args.feature_id:
            continue
        if args.proteoform_id and cell(row, mapping, "proteoform_id") != args.proteoform_id:
            continue
        if args.prsm_id and cell(row, mapping, "prsm_id") != args.prsm_id:
            continue
        if args.max_evalue is not None:
            value = as_float(cell(row, mapping, "e_value"))
            if value is None or value > args.max_evalue:
                continue
        if args.max_qvalue is not None:
            value = as_float(cell(row, mapping, "q_value"))
            if value is None or value > args.max_qvalue:
                continue
        if args.protein:
            haystack = f"{cell(row, mapping, 'accession')} {cell(row, mapping, 'description')}".lower()
            if args.protein.lower() not in haystack:
                continue
        if args.scan and args.scan not in re.split(r"[,;\s]+", cell(row, mapping, "scan")):
            continue
        shifts = unexpected_shifts(row, mapping)
        if args.with_unexpected and not shifts:
            continue
        if args.min_abs_shift is not None and not any(abs(item) >= args.min_abs_shift for item in shifts):
            continue
        payload = dict(row) if args.all_columns else {
            key: cell(row, mapping, key) for key in mapping if cell(row, mapping, key) != ""
        }
        payload["_unexpected_shifts"] = shifts
        selected.append(payload)

    return {
        "query": {
            "kind": "rows",
            "path": args.path,
            "max_evalue": args.max_evalue,
            "max_qvalue": args.max_qvalue,
            "protein": args.protein,
            "scan": args.scan,
            "feature_id": args.feature_id,
            "proteoform_id": args.proteoform_id,
            "prsm_id": args.prsm_id,
            "with_unexpected": args.with_unexpected,
            "min_abs_shift": args.min_abs_shift,
        },
        "match_count": len(selected),
        "truncated": len(selected) > args.limit,
        "rows": selected[: args.limit],
    }


def add_shared(target) -> None:
    """Options accepted on either side of the subcommand.

    argparse resets a parent option to its default while parsing a subcommand,
    so shared options are declared with SUPPRESS and their base values come
    from set_defaults.
    """
    target.add_argument("--out", default=argparse.SUPPRESS, help="also write the JSON result to this path")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_shared(parser)
    parser.set_defaults(out=None)
    sub = parser.add_subparsers(dest="command", required=True)

    for name, help_text in (
        ("schema", "columns, row count and recognized fields"),
        ("summary", "identification counts and score distributions"),
        ("rows", "filtered rows"),
    ):
        target = sub.add_parser(name, help=help_text)
        add_shared(target)
        target.add_argument("--path", required=True, help="TopPIC *_prsm.tsv or *_proteoform.tsv")
        if name == "rows":
            target.add_argument("--max-evalue", type=float)
            target.add_argument("--max-qvalue", type=float)
            target.add_argument("--protein", help="substring of accession or description")
            target.add_argument("--scan", help="exact scan number")
            target.add_argument("--feature-id", help="exact native TopPIC Feature ID")
            target.add_argument("--proteoform-id", help="exact native TopPIC Proteoform ID")
            target.add_argument("--prsm-id", help="exact native TopPIC Prsm ID")
            target.add_argument("--with-unexpected", action="store_true", help="only rows carrying an unexpected shift")
            target.add_argument("--min-abs-shift", type=float, help="only rows with |shift| >= this value")
            target.add_argument("--all-columns", action="store_true", help="emit every column, not the recognized subset")
            target.add_argument("--limit", type=int, default=50)
    return parser


def main(argv: list[str]) -> int:
    args = build_parser().parse_args(argv)
    try:
        path = Path(args.path)
        columns, rows, header_line = read_table(path)
        mapping = resolve_columns(columns)
        handlers = {"schema": cmd_schema, "summary": cmd_summary, "rows": cmd_rows}
        result = handlers[args.command](args, columns, rows, header_line, mapping)
        result["source"] = {"path": str(path), "sha256": digest(path), "header_line": header_line}
        text = json.dumps(result, ensure_ascii=False, indent=2)
        if args.out:
            out = Path(args.out)
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(f"{text}\n", encoding="utf-8")
            result["written_to"] = str(out)
            text = json.dumps(result, ensure_ascii=False, indent=2)
        sys.stdout.write(f"{text}\n")
        return 0
    except ToolError as exc:
        json.dump({"error": {"code": exc.code, "message": exc.message, **exc.detail}}, sys.stdout, indent=2)
        sys.stdout.write("\n")
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
