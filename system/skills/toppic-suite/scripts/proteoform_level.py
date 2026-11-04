# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Assign a proteoform identification level from a TopPIC PrSM table.

Implements the counting rule of the five-level classification (Smith et al.,
Nat Methods 2019): four independent ambiguity categories, and

    level = 1 + number of ambiguous categories        (1 = fully defined, 5 = mass only)

The four categories are gene/protein identity, amino-acid sequence, PTM identity
and PTM localization, all derived from the row: a missing accession or several
protein hits, a missing database sequence or proteoform string, a bare-mass
modification, and a localization span or region wider than one residue. Sequence
ambiguity that the row cannot show -- an isobaric alternative sequence, an
unreported variant -- is yours to declare with --sequence-ambiguous.

Each row also carries a subtype letter when exactly one category is ambiguous,
and an E-value triage zone, both matching the earlier TDEase implementation so
that older result tables remain comparable.

Every row reports *why* each category was called ambiguous, so a level can be
argued with instead of trusted.

Prints one JSON object on stdout. Read-only unless --out is given.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path

SIGNATURE_COLUMNS = ("Prsm ID", "Spectrum ID", "Protein accession", "Proteoform", "Proteoform ID")
MASS_TEXT = re.compile(r"^[+-]?\d+(?:\.\d+)?$")
CATEGORIES = ("gene", "sequence", "ptm_identity", "ptm_localization")

# Subclass letters follow the TDEase convention, carried over from the previous
# implementation so historical result tables stay comparable.
SUBTYPE_BY_CATEGORY = {
    "ptm_localization": "2A",
    "ptm_identity": "2B",
    "sequence": "2C",
    "gene": "2D",
}

# A numeric region written into the unexpected-modification column, e.g. [12-18].
REGION = re.compile(r"\[(\d+)(?:-(\d+))?\]")

UNKNOWN_ACCESSIONS = {"", "unknown", "unassigned", "na", "n/a", "-"}

# E-value zones. A zone is a triage band for the spectrum match, nothing more:
# it never speaks to modification identity, localization or biology.
EVALUE_ZONES = (("Gold", 1e-7), ("Silver", 1e-3))
ZONE_ORDER = {"Gold": 0, "Silver": 1, "Trash": 2, "unknown": 3}


# The shared plotting helper and the asset-path resolver live in the System
# reference root (data/system/reference/scripts), not beside this script.
_plot_dir = str(Path(__file__).resolve().parents[3] / 'reference' / 'scripts')
if _plot_dir not in sys.path:
    sys.path.insert(0, _plot_dir)
import _plot  # noqa: E402
import system_paths  # noqa: E402


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


def read_table(path: Path) -> tuple[list[dict], int]:
    if not path.is_file():
        raise ToolError("FILE_NOT_FOUND", f"{path} does not exist")
    lines = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    header_index = None
    for index, line in enumerate(lines):
        if not line.strip() or line.lstrip().startswith(("*", "#!")) or "\t" not in line:
            continue
        if any(name in [cell.strip() for cell in line.split("\t")] for name in SIGNATURE_COLUMNS):
            header_index = index
            break
    if header_index is None:
        raise ToolError("NOT_A_TOPPIC_TABLE", "no TopPIC header row found", path=str(path))
    reader = csv.DictReader(lines[header_index:], delimiter="\t")
    rows = []
    for raw in reader:
        row = {(key or "").strip(): (value or "").strip() for key, value in raw.items() if key is not None}
        if any(row.values()):
            rows.append(row)
    return rows, header_index + 1


def column(row: dict, *names: str) -> str:
    for name in names:
        if name in row and row[name] != "":
            return row[name]
    return ""


def parse_proteoform(text: str) -> tuple[str, list]:
    body = text or ""
    parts = body.split(".")
    if len(parts) >= 3:
        body = ".".join(parts[1:-1])
    elif len(parts) == 2:
        body = parts[1]

    residues: list[str] = []
    mods: list[dict] = []
    span: tuple[int, int] | None = None
    index = 0
    while index < len(body):
        char = body[index]
        if char == "(":
            close = body.find(")", index)
            if close == -1:
                break
            start = len(residues) + 1
            residues.extend(item for item in body[index + 1 : close] if item.isalpha())
            span = (start, len(residues))
            index = close + 1
        elif char == "[":
            close = body.find("]", index)
            if close == -1:
                break
            content = body[index + 1 : close].strip()
            if span is not None:
                start, end = span
            elif residues:
                start = end = len(residues)
            else:
                start = end = 0
            is_mass = bool(MASS_TEXT.match(content))
            mods.append(
                {
                    "text": content,
                    "mass": float(content) if is_mass else None,
                    "named": not is_mass,
                    "span_start": start,
                    "span_end": end,
                    "span_length": max(1, end - start + 1) if start else 1,
                    "span_residues": "".join(residues[start - 1 : end]) if start else "",
                }
            )
            span = None
            index = close + 1
        elif char.isalpha():
            residues.append(char)
            span = None
            index += 1
        else:
            index += 1
    return "".join(residues), mods


def as_float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def evalue_zone(value: float | None) -> str:
    if value is None:
        return "unknown"
    for label, threshold in EVALUE_ZONES:
        if value < threshold:
            return label
    return "Trash"


def classify(row: dict, args) -> dict:
    _, mods = parse_proteoform(column(row, "Proteoform"))
    reasons: dict[str, str] = {}
    derived: dict[str, bool] = {}

    accession = column(row, "Protein accession").strip()
    sequence = column(row, "Database protein sequence").strip()
    proteoform = column(row, "Proteoform").strip()

    hits = as_float(column(row, "#Protein hits", "Protein hits"))
    unknown_accession = accession.lower() in UNKNOWN_ACCESSIONS
    gene_ambiguous = bool(hits and hits > 1) or unknown_accession
    if unknown_accession:
        reasons["gene"] = "no usable protein accession on this row"
    elif gene_ambiguous:
        reasons["gene"] = f"{int(hits)} protein hits reported for this spectrum"
    derived["gene"] = True

    # Derivable after all: a row with no database sequence or no proteoform
    # string does not define the molecule, whatever else it reports.
    missing_sequence = not sequence or not proteoform
    sequence_ambiguous = args.sequence_ambiguous or missing_sequence
    if args.sequence_ambiguous:
        reasons["sequence"] = "declared ambiguous by the caller"
    elif missing_sequence:
        reasons["sequence"] = "row carries no database sequence or no proteoform string"
    derived["sequence"] = missing_sequence or not args.sequence_ambiguous

    unnamed = [mod for mod in mods if mod["mass"] is not None]
    ptm_identity_ambiguous = bool(unnamed)
    if ptm_identity_ambiguous:
        reasons["ptm_identity"] = (
            f"{len(unnamed)} modification(s) reported as a bare mass: "
            + ", ".join(mod["text"] for mod in unnamed[:4])
        )
    derived["ptm_identity"] = True

    # TopPIC states the unexpected count in its own column. If the proteoform
    # string does not carry that many mass shifts, the string was not parsed as
    # written -- call the identity ambiguous and say so rather than returning a
    # confident Level 1 built on a failed parse.
    warnings = []
    declared = as_float(column(row, "#unexpected modifications", "#Unexpected modifications"))
    if declared is not None and int(declared) > len(unnamed):
        warnings.append(
            f"row declares {int(declared)} unexpected modification(s) but the proteoform string yielded {len(unnamed)}"
        )
        if not ptm_identity_ambiguous:
            ptm_identity_ambiguous = True
            reasons["ptm_identity"] = warnings[-1]

    wide = [mod for mod in mods if mod["span_length"] > 1]
    # Some TopPIC versions write the window into the modification column as a
    # numeric region instead of a parenthesised span in the proteoform string.
    regions = REGION.findall(column(row, "unexpected modifications", "Unexpected modifications"))
    wide_regions = [(start, end) for start, end in regions if end and end != start]
    mi_score = as_float(column(row, "MIScore", "MIScore(s)"))
    ptm_localization_ambiguous = bool(wide or wide_regions)
    if wide:
        reasons["ptm_localization"] = (
            "modification localized only to a span of "
            + ", ".join(f"{mod['span_length']} residues ({mod['span_residues']})" for mod in wide[:3])
        )
    elif wide_regions:
        reasons["ptm_localization"] = (
            "modification reported over region(s) " + ", ".join(f"{start}-{end}" for start, end in wide_regions[:3])
        )
    elif mods and mi_score is not None and mi_score < args.mi_score_threshold:
        ptm_localization_ambiguous = True
        reasons["ptm_localization"] = f"MIScore {mi_score} below the {args.mi_score_threshold} threshold"
    derived["ptm_localization"] = True

    flags = {
        "gene": gene_ambiguous,
        "sequence": sequence_ambiguous,
        "ptm_identity": ptm_identity_ambiguous,
        "ptm_localization": ptm_localization_ambiguous,
    }
    ambiguities = [name for name in CATEGORIES if flags[name]]
    # An observed mass with no protein, no sequence and no proteoform is the
    # bottom of the scale by definition, not by counting.
    observed_mass_only = not accession and not sequence and not proteoform
    level = 5 if observed_mass_only else min(5, 1 + len(ambiguities))
    e_value = as_float(column(row, "E-value"))
    return {
        "prsm_id": column(row, "Prsm ID", "PrSM ID"),
        "scan": column(row, "Scan(s)", "Scans", "Scan"),
        "accession": accession,
        "proteoform": proteoform,
        "e_value": e_value,
        "evalue_zone": evalue_zone(e_value),
        "mi_score": mi_score,
        "level": level,
        "subtype": SUBTYPE_BY_CATEGORY.get(ambiguities[0]) if len(ambiguities) == 1 else None,
        "ambiguities": ambiguities,
        "reasons": reasons,
        "derived_from_data": derived,
        "modification_count": len(mods),
        "observed_mass_only": observed_mass_only,
        "suggest_verification": flags["ptm_identity"] or flags["ptm_localization"],
        "warnings": warnings,
    }


def level_figure(distribution: dict, zones: dict, category_counts: dict) -> list:
    if not distribution:
        return []
    plot = _plot.Plot(
        title="Proteoform identification level",
        subtitle="level = 1 + ambiguous categories; higher means less completely defined",
        x_label="level",
        y_label="rows",
        height=320,
        width=560,
    )
    plot.bar(
        [(int(level), count) for level, count in sorted(distribution.items())],
        annotations=[(int(level), count, str(count)) for level, count in sorted(distribution.items())],
    )
    figures = [("levels", "level-distribution", "Proteoform identification level", plot, plot.series_data())]

    if category_counts:
        order = ["gene", "sequence", "ptm_identity", "ptm_localization"]
        present = [(index, category_counts.get(name, 0)) for index, name in enumerate(order)]
        ambiguity = _plot.Plot(
            title="Which category is ambiguous",
            subtitle="counted across every classified row",
            y_label="rows",
            height=300,
            width=560,
        )
        ambiguity.bar(present, annotations=[(index, value, order[index]) for index, value in present if value])
        figures.append((
            "ambiguity", "ambiguity-breakdown", "Which category is ambiguous",
            ambiguity, ambiguity.series_data(),
        ))
    return figures


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--path", required=True, help="TopPIC *_prsm.tsv")
    parser.add_argument("--max-evalue", type=float, help="classify only rows at or below this E-value")
    parser.add_argument(
        "--mi-score-threshold",
        type=float,
        default=0.5,
        help="MIScore below which a single-residue localization is still called ambiguous",
    )
    parser.add_argument(
        "--sequence-ambiguous",
        action="store_true",
        help="force sequence ambiguity on every row, on top of the rows where it is derivable",
    )
    parser.add_argument("--limit", type=int, default=200)
    parser.add_argument("--out", help="also write the JSON result to this path")
    _plot.add_plot_arguments(parser)
    args = parser.parse_args(argv)

    try:
        path = Path(args.path)
        rows, header_line = read_table(path)
        classified = []
        for row in rows:
            if args.max_evalue is not None:
                value = as_float(column(row, "E-value"))
                if value is None or value > args.max_evalue:
                    continue
            classified.append(classify(row, args))

        classified.sort(key=lambda item: (
            ZONE_ORDER.get(item["evalue_zone"], 3),
            item["e_value"] if item["e_value"] is not None else float("inf"),
        ))
        distribution = Counter(item["level"] for item in classified)
        subtypes = Counter(item["subtype"] for item in classified if item["subtype"])
        zones = Counter(item["evalue_zone"] for item in classified)
        category_counts = Counter(name for item in classified for name in item["ambiguities"])
        result = {
            "query": {
                "kind": "proteoform-level",
                "path": str(path),
                "max_evalue": args.max_evalue,
                "mi_score_threshold": args.mi_score_threshold,
                "sequence_ambiguous": args.sequence_ambiguous,
            },
            "row_count": len(classified),
            "level_distribution": {str(level): distribution[level] for level in sorted(distribution)},
            "subtype_distribution": dict(sorted(subtypes.items())),
            "evalue_zones": dict(sorted(zones.items(), key=lambda pair: ZONE_ORDER.get(pair[0], 3))),
            "ambiguity_counts": dict(category_counts),
            "rows_suggesting_verification": sum(1 for item in classified if item["suggest_verification"]),
            "truncated": len(classified) > args.limit,
            "rows": classified[: args.limit],
            "source": {"path": str(path), "sha256": digest(path), "header_line": header_line},
            "method": "level = 1 + ambiguous categories among gene, sequence, ptm_identity, ptm_localization (Smith et al. 2019); subtype letters follow the TDEase convention 2A localization, 2B identity, 2C sequence, 2D gene",
            "caution": "A level is a completeness statement about the identification, not a confidence score, and the E-value zone is a triage band for the spectrum match only. Rows are sorted by zone then E-value so the strongest matches read first.",
        }
        _plot.emit_figures(
            result, args,
            level_figure(result["level_distribution"], result["evalue_zones"], result["ambiguity_counts"]),
        )
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
