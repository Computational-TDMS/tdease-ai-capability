"""Small, lossless readers and string helpers for TopPIC and TopFD tables."""

import csv
import re
from functools import lru_cache
from pathlib import Path

from _shared.toppic import parse_mods

@lru_cache(maxsize=1)
def c13_gap():
    """Kernel C13-C12 spacing for coarse grouping, not isotope-shape theory."""
    from _shared.native import mass_constants
    return mass_constants()["carbon13Spacing"]


def ref_accession(ref: str) -> str:
    """Extract accession from ``feature|proteoform|accession``."""
    parts = str(ref or "").split("|", 2)
    return parts[2] if len(parts) == 3 else ""


def number(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def prsm_passes(row, qvalue_max=0.01, evalue_fallback=1e-5):
    qvalues = [value for value in (
        number(row.get(field), None) for field in
        ("Spectrum-level Q-value", "Proteoform-level Q-value")
    ) if value is not None]
    if qvalues:
        return all(value <= qvalue_max for value in qvalues)
    return number(row.get("E-value"), 1.0) <= evalue_fallback


def read_toppic(path):
    """Read a TopPIC TSV after its parameter preamble."""
    with open(path, encoding="utf-8", newline="") as fh:
        lines = fh.readlines()
    header = next(
        (index for index, line in enumerate(lines)
         if line.startswith("Data file name\tPrsm ID\t")), None)
    if header is None:
        raise ValueError(f"TopPIC table header not found: {path}")
    return list(csv.DictReader(lines[header:], delimiter="\t"))


def read_toppic_parameters(path):
    parameters = {}
    for line in Path(path).read_text(encoding="utf-8-sig").splitlines():
        if line.startswith("Data file name\tPrsm ID\t"):
            break
        if "\t" not in line or line.startswith("*"):
            continue
        name, value = line.split("\t", 1)
        name, value = name.rstrip(": "), value.strip()
        if name and value:
            parameters[name] = value
    return parameters


def read_ms1(path):
    with open(path, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    out = {}
    for row in rows:
        feature_id = row.get("Feature_ID") or row.get("Feature ID") or ""
        mass = number(row.get("Mass"), None)
        if not feature_id or mass is None:
            continue
        out[feature_id] = {
            "source_feature_id": feature_id,
            "mass": mass,
            "intensity": number(row.get("Intensity")),
            "rt_start": number(row.get("Min_time")),
            "rt_end": number(row.get("Max_time")),
            "rt_apex": number(row.get("Apex_time")),
            "min_charge": row.get("Min_charge") or "",
            "max_charge": row.get("Max_charge") or "",
            "ec_score": number(row.get("EC_score")),
        }
    return out


def read_feature_rows(path):
    if not path.is_file():
        raise ValueError(f"feature table not found: {path}")
    rows = []
    with path.open(encoding="utf-8-sig", errors="replace") as fh:
        for raw in csv.DictReader(fh, delimiter="\t"):
            try:
                min_charge = int(float(raw.get("Min_charge") or 0))
                rep_charge = int(float(raw.get("Rep_charge") or 0))
                row = {
                    "feature_id": int(float(raw.get("Feature_ID") or 0)),
                    "mass": float(raw.get("Mass") or 0),
                    "intensity": float(raw.get("Intensity") or 0),
                    "rt_apex": float(raw.get("Apex_time") or 0),
                    "rt_min": float(raw.get("Min_time") or raw.get("Apex_time") or 0),
                    "rt_max": float(raw.get("Max_time") or raw.get("Apex_time") or 0),
                    "scan_min": int(float(raw.get("Min_scan") or 0)),
                    "scan_max": int(float(raw.get("Max_scan") or 0)),
                    # representative charge anchors the envelope refine (mono mz)
                    "rep_charge": rep_charge or max(1, min_charge),
                    "ec_score": float(raw.get("EC_score") or 0),
                }
            except (TypeError, ValueError):
                continue
            if row["feature_id"] >= 0 and row["mass"] > 0:
                rows.append(row)
    if not rows:
        raise ValueError(f"no parseable MS1 features: {path}")
    return sorted(rows, key=lambda row: (row["mass"], row["feature_id"]))


_FIXED_NTERM = re.compile(r"^\[([A-Za-z][^\]-]*)\]-")
_FIXED_CTERM = re.compile(r"-\[([A-Za-z][^\]]*)\]$")
_UNKNOWN_MOD = re.compile(r"\(([A-Z]+)\)\[([+-]\d+\.?\d*)\]")


def parse_proforma(proforma, full_sequence):
    """Parse a TopPIC ProForma string into structured facts.

    Returns (fixed, unknown):
      fixed:  [{"position": "nTerm"|"cTerm", "name": str}]
      unknown: [{"mass": float, "residues": str, "interval": [start, end]}
                in FULL-CHAIN 1-based coordinates (start/end inclusive),
                located by matching the residue span in full_sequence.
                None interval when the span is not found or occurs more than
                once (never guessed).
    """
    fixed = []
    m = _FIXED_NTERM.match(proforma or "")
    if m:
        fixed.append({"position": "nTerm", "name": m.group(1)})
    m = _FIXED_CTERM.search(proforma or "")
    if m:
        fixed.append({"position": "cTerm", "name": m.group(1)})
    unknown = []
    for m in _UNKNOWN_MOD.finditer(proforma or ""):
        residues = m.group(1)
        mass = float(m.group(2))
        interval = None
        # Locate by content only when the span is unambiguous: a repeated
        # residue stretch must not yield a confident-but-wrong interval.
        positions = []
        search_from = 0
        while (pos := full_sequence.find(residues, search_from)) != -1:
            positions.append(pos)
            search_from = pos + 1
        if len(positions) == 1:
            interval = [positions[0] + 1, positions[0] + len(residues)]
        unknown.append({"mass": mass, "residues": residues, "interval": interval})
    return fixed, unknown
