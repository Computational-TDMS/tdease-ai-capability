"""
Smith proteoform-definition completeness scale, from Smith et al. (2013).
classify_row is the importable kernel used by ingest; per-row results land in
the durable projection. The scale itself is documented in
tdp-copilot/references/evidence-rules.md (Smith proteoform level).
"""

from __future__ import annotations
import re

MASS_TEXT = re.compile(r"^[+-]?\d+(?:\.\d+)?$")

REGION = re.compile(r"\[(\d+)(?:-(\d+))?\]")

CATEGORIES = ("gene", "sequence", "ptm_identity", "ptm_localization")

SUBTYPE_BY_CATEGORY = {
    "ptm_localization": "2A",
    "ptm_identity": "2B",
    "sequence": "2C",
    "gene": "2D",
}

UNKNOWN_ACCESSIONS = {"", "unknown", "unassigned", "na", "n/a", "-"}

def column(row, *names):
    return next((row[name] for name in names if row.get(name)), "")

def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        return None

def parse_proteoform(text):
    body = text or ""
    residues = []
    modifications = []
    span = None
    index = 0
    while index < len(body):
        char = body[index]
        if char == "(":
            close = body.find(")", index)
            if close == -1:
                break
            start = len(residues) + 1
            residues.extend(item for item in body[index + 1:close] if item.isalpha())
            span = (start, len(residues))
            index = close + 1
        elif char == "[":
            close = body.find("]", index)
            if close == -1:
                break
            content = body[index + 1:close].strip()
            start, end = span or ((len(residues), len(residues)) if residues else (0, 0))
            modifications.append({
                "mass": float(content) if MASS_TEXT.match(content) else None,
                "span_length": max(1, end - start + 1) if start else 1,
            })
            span = None
            index = close + 1
        elif char.isalpha():
            residues.append(char)
            span = None
            index += 1
        else:
            index += 1
    return modifications

def classify_row(row, mi_score_threshold=0.75, sequence_ambiguous=False):
    """Return Smith 1-5 completeness and its ambiguity categories for one row."""
    modifications = parse_proteoform(column(row, "Proteoform"))
    accession = column(row, "Protein accession").strip()
    sequence = column(row, "Database protein sequence").strip()
    proteoform = column(row, "Proteoform").strip()
    hits = as_float(column(row, "# Protein hits", "#Protein hits", "Protein hits"))
    gene_ambiguous = bool(hits and hits > 1) or accession.lower() in UNKNOWN_ACCESSIONS
    sequence_is_ambiguous = sequence_ambiguous or not sequence or not proteoform

    unnamed = [item for item in modifications if item["mass"] is not None]
    identity_ambiguous = bool(unnamed)
    declared = as_float(column(row, "# Unexpected modifications", "#unexpected modifications"))
    if declared is not None and int(declared) > len(unnamed):
        identity_ambiguous = True

    wide = any(item["span_length"] > 1 for item in modifications)
    regions = REGION.findall(column(row, "Unexpected modifications", "unexpected modifications"))
    wide_region = any(end and end != start for start, end in regions)
    mi_score = as_float(column(row, "MIScore", "MIScore(s)"))
    localization_ambiguous = wide or wide_region
    if modifications and mi_score is not None and mi_score < mi_score_threshold:
        localization_ambiguous = True

    flags = {
        "gene": gene_ambiguous,
        "sequence": sequence_is_ambiguous,
        "ptm_identity": identity_ambiguous,
        "ptm_localization": localization_ambiguous,
    }
    ambiguities = [name for name in CATEGORIES if flags[name]]
    observed_mass_only = not accession and not sequence and not proteoform
    level = 5 if observed_mass_only else min(5, 1 + len(ambiguities))
    return {
        "level": level,
        "subtype": SUBTYPE_BY_CATEGORY.get(ambiguities[0]) if len(ambiguities) == 1 else "",
        "ambiguities": ambiguities,
    }
