"""
Read TopPIC output files into pandas DataFrames.

Handles both MS2 (PrSM) and MS1 (feature) file formats across multiple
TopPIC versions (1.5.4, 1.7.0, 1.8.1+).  Version is auto-detected from
the column headers of the first file.

The main entry point is :func:`read_toppic`.
"""

import re
from pathlib import Path
from typing import Optional

import pandas as pd

from toppic.utils.fasta import extract_gene, extract_uniprot_acc

# ═══════════════════════════════════════════════════════════════════════
# Version detection
# ═══════════════════════════════════════════════════════════════════════

# v1.8.1+ MS2 columns that are NOT in the pipeline keep list
# but are kept in the output DataFrame for user reference.
_MS2_EXTRA_COLS = [
    "Proteoform intensity",
    "Feature ID",
    "Database protein sequence",
    "Previous amino acid",
    "Next amino acid",
    "Fixed PTMs",
    "Unexpected modifications",
    "Variable PTMs",
    "Envelope_num",
]


def _detect_ms2_version(header_line: str) -> str:
    """Detect TopPIC MS2 version from the tab-separated header line."""
    cols = header_line.strip().split("\t")
    if "# Masses" in cols:
        return "1.8"
    elif "#peaks" in cols:
        return "1.5"
    elif "Feature ID" in cols and "Database protein sequence" in cols:
        return "1.7"
    return "unknown"


def _detect_ms1_version(header_line: str) -> str:
    """Detect TopPIC MS1 version from the tab-separated header line."""
    cols = header_line.strip().split("\t")
    if "Apex_time" in cols:
        return "1.8"
    elif "Time_apex" in cols:
        return "1.5"
    elif "Monoisotopic_mass" in cols:
        return "1.7"
    return "unknown"


# ═══════════════════════════════════════════════════════════════════════
# Column maps: TopPIC raw names → canonical snake_case
# One map per supported version
# ═══════════════════════════════════════════════════════════════════════

# ─── MS2 v1.5.4 (JPR data) ────────────────────────────────────────────

_MS2_COLUMN_MAP_154 = {
    "Data file name": "data_file_name",
    "Prsm ID": "prsm_id",
    "Spectrum ID": "spectrum_id",
    "Fragmentation": "fragmentation",
    "Scan(s)": "scans",
    "Retention time": "retention_time",
    "#peaks": "n_peaks",
    "Charge": "charge",
    "Precursor mass": "precursor_mass",
    "Adjusted precursor mass": "adjusted_precursor_mass",
    "Proteoform ID": "proteoform_id",
    "Feature intensity": "feature_intensity",
    "Feature score": "feature_score",
    "Feature apex time": "feature_apex_time",
    "#Protein hits": "n_protein_hits",
    "Protein accession": "protein_accession",
    "Protein description": "protein_description",
    "First residue": "first_residue",
    "Last residue": "last_residue",
    "Special amino acids": "special_amino_acids",
    "Proteoform": "proteoform",
    "Proteoform mass": "proteoform_mass",
    "Protein N-terminal form": "protein_nterm_form",
    "#unexpected modifications": "n_unexpected_mods",
    "#variable PTMs": "n_variable_ptms",
    "MIScore": "mi_score",
    "#matched peaks": "n_matched_peaks",
    "#matched fragment ions": "n_matched_fragment_ions",
    "E-value": "e_value",
    "Spectrum-level Q-value": "spectrum_q_value",
    "Proteoform-level Q-value": "proteoform_q_value",
}

# ─── MS2 v1.8.1 (current TopPIC version) ──────────────────────────────

_MS2_COLUMN_MAP_181 = {
    "Data file name": "data_file_name",
    "Prsm ID": "prsm_id",
    "Spectrum ID": "spectrum_id",
    "Fragmentation": "fragmentation",
    "Scan(s)": "scans",
    "Retention time": "retention_time",
    "# Masses": "n_masses",
    "Charge": "charge",
    "Precursor mass": "precursor_mass",
    "Adjusted precursor mass": "adjusted_precursor_mass",
    "Proteoform ID": "proteoform_id",
    "Proteoform intensity": "_proteoform_intensity",
    "Feature ID": "_feature_id",
    "Feature intensity": "feature_intensity",
    "Feature score": "feature_score",
    "Feature apex time": "feature_apex_time",
    "# Protein hits": "n_protein_hits",
    "Protein accession": "protein_accession",
    "Protein description": "protein_description",
    "First residue": "first_residue",
    "Last residue": "last_residue",
    "Special amino acids": "special_amino_acids",
    "Database protein sequence": "_db_protein_sequence",
    "Previous amino acid": "_prev_aa",
    "Proteoform": "proteoform",
    "Next amino acid": "_next_aa",
    "Proteoform mass": "proteoform_mass",
    "Protein N-terminal form": "protein_nterm_form",
    "Fixed PTMs": "_fixed_ptms",
    "# Unexpected modifications": "n_unexpected_mods",
    "Unexpected modifications": "_unexpected_mods",
    "# Variable PTMs": "n_variable_ptms",
    "Variable PTMs": "_variable_ptms",
    "MIScore": "mi_score",
    "# Matched masses": "n_matched_masses",
    "# Matched fragments": "n_matched_fragments",
    "E-value": "e_value",
    "Spectrum-level Q-value": "spectrum_q_value",
    "Proteoform-level Q-value": "proteoform_q_value",
}

# Map version string → MS2 column map
_MS2_VERSION_MAP = {
    "1.5": _MS2_COLUMN_MAP_154,
    "1.7": _MS2_COLUMN_MAP_154,  # v1.7 mostly compatible with 1.5 for our cols
    "1.8": _MS2_COLUMN_MAP_181,
}

# ─── MS1 v1.5.4 ───────────────────────────────────────────────────────

_MS1_COLUMN_MAP_154 = {
    "Mass": "mass",
    "Intensity": "intensity",
    "Time_apex": "time_apex",
}

# ─── MS1 v1.8.1 ───────────────────────────────────────────────────────

_MS1_COLUMN_MAP_181 = {
    "File_name": "_file_name",
    "Fraction_ID": "_fraction_id",
    "Feature_ID": "_feature_id",
    "Mass": "mass",
    "Intensity": "intensity",
    "Min_time": "_min_time",
    "Max_time": "_max_time",
    "Min_scan": "_min_scan",
    "Max_scan": "_max_scan",
    "Min_charge": "_min_charge",
    "Max_charge": "_max_charge",
    "Apex_time": "time_apex",
    "Apex_scan": "_apex_scan",
    "Apex_intensity": "_apex_intensity",
    "Rep_charge": "_rep_charge",
    "Rep_average_mz": "_rep_average_mz",
    "Envelope_num": "_envelope_num",
    "EC_score": "_ec_score",
}

_MS1_VERSION_MAP = {
    "1.5": _MS1_COLUMN_MAP_154,
    "1.7": _MS1_COLUMN_MAP_154,
    "1.8": _MS1_COLUMN_MAP_181,
}

# ═══════════════════════════════════════════════════════════════════════
# Pipeline keep-list columns
# ═══════════════════════════════════════════════════════════════════════

# Columns the downstream pipeline requires (always present in output)
_MS2_KEEP = [
    "dataset",
    "scans",
    "retention_time",
    "feature_apex",
    "charge",
    "mz",
    "precursor_mass",
    "adjusted_precursor_mass",
    "feature_intensity",
    "e_value",
    "n_unexpected_mods",
    "acc_map",
    "mi_score",
    "proteoform",
    "special_amino_acids",
    "uniprot_acc",
    "gene",
    "is_decoy",
]

# Columns the pipeline always keeps from MS1
_MS1_KEEP = ["dataset", "mass", "intensity", "time_apex"]


# ═══════════════════════════════════════════════════════════════════════
# Public API
# ═══════════════════════════════════════════════════════════════════════


def read_toppic(
    file_path: str | Path,
    file_name: list[str],
    faims: bool = False,
    **kwargs,
) -> pd.DataFrame:
    """Read TopPIC output files into a single DataFrame.

    Auto-detects whether files are MS2 (PrSM) or MS1 (feature) output
    based on the first file's name, and auto-detects the TopPIC version
    (1.5.4 / 1.7.0 / 1.8.1+) from column headers.  All files of the
    same type are concatenated.

    Parameters
    ----------
    file_path : str | Path
        Path to the folder containing TopPIC output files.
    file_name : list[str]
        Names of the files to read.  All files must be of the **same**
        type (MS2 PrSM or MS1 feature).
    faims : bool, default=False
        If True, extract FAIMS CV from the ``dataset`` column.
    **kwargs
        Extra arguments forwarded to ``pandas.read_csv``.

    Returns
    -------
    pd.DataFrame
        Cleaned and annotated DataFrame.  Pipeline-essential columns
        are always present (see ``_MS2_KEEP`` / ``_MS1_KEEP``).
        Version-specific extra columns (v1.8.1+) are retained with
        their ``_`` -prefixed snake_case names.
    """
    file_path = Path(file_path)
    if not file_name:
        raise ValueError("file_name must not be empty")

    first = file_name[0]
    if "ms1.feature" in first:
        return _read_ms1(file_path, file_name)
    else:
        return _read_ms2(file_path, file_name, faims, **kwargs)


# ═══════════════════════════════════════════════════════════════════════
# MS2 (PrSM) reader
# ═══════════════════════════════════════════════════════════════════════


def _read_ms2(
    file_path: Path,
    file_name: list[str],
    faims: bool,
    **kwargs,
) -> pd.DataFrame:
    # Detect version from the first file
    first_path = file_path / file_name[0]
    _validate_file(first_path)
    n_skip = _find_header_line(first_path)

    with open(first_path) as fh:
        for _ in range(n_skip):
            next(fh)
        header_line = next(fh)

    version = _detect_ms2_version(header_line)
    if version == "unknown":
        # An unrecognized header read through the wrong column map yields
        # confident garbage downstream (renames miss, derived columns
        # KeyError or misalign). Refuse instead of silently guessing 1.5.4.
        raise ValueError(
            f"Unrecognized TopPIC MS2 table version in {first_path}: no known "
            "version marker (# Masses for 1.8, #peaks for 1.5, Feature ID + "
            "Database protein sequence for 1.7) found in the header. "
            "Convert the table to a supported version's column schema."
        )
    col_map = _MS2_VERSION_MAP[version]

    tables = []
    for fname in file_name:
        full_path = file_path / fname
        _validate_file(full_path)

        n_skip = _find_header_line(full_path)
        raw = pd.read_csv(
            full_path,
            sep="\t",
            skiprows=n_skip,
            dtype={"Special amino acids": str},
            **kwargs,
        )

        # Rename columns using version-specific map
        raw = raw.rename(columns=col_map, errors="ignore")

        # Add derived columns (version-independent)
        raw["dataset"] = raw["data_file_name"].str.replace(
            "_ms2.msalign", "", regex=False
        )
        raw["dataset"] = raw["dataset"].str.split("/").str[-1]

        raw["mz"] = (raw["precursor_mass"] + raw["charge"] * 1.007276466621) / raw[
            "charge"
        ]

        raw["gene"] = raw["protein_description"].fillna("").apply(extract_gene)
        raw["is_decoy"] = raw["protein_accession"].fillna("").str.startswith("DECOY")
        raw["uniprot_acc"] = raw["protein_accession"].apply(_safe_uniprot_acc)
        raw["acc_map"] = raw["n_protein_hits"]

        # Rename feature_apex_time → feature_apex (both versions use this)
        raw["feature_apex"] = raw["feature_apex_time"]

        # Filter rows without a gene in the description
        raw = raw[raw["protein_description"].fillna("").str.contains("GN=", na=False)]

        # Fill missing special amino acids
        raw["special_amino_acids"] = raw["special_amino_acids"].fillna("-")

        # TopPIC continuation records follow their complete PrSM record. R
        # fills down within each file in that original order. Sorting first
        # changes continuation assignments and downstream cluster tie order.
        _fill_cols = [
            "scans", "retention_time", "feature_apex", "charge", "mz",
            "precursor_mass", "adjusted_precursor_mass", "feature_intensity",
            "e_value", "n_unexpected_mods", "acc_map", "mi_score", "proteoform",
        ]
        raw[_fill_cols] = raw[_fill_cols].ffill()

        # Select pipeline keep-columns + any version-specific extras
        keep_cols = list(_MS2_KEEP)
        for c in raw.columns:
            if c not in keep_cols and c.startswith("_"):
                keep_cols.append(c)
        raw = raw[[c for c in keep_cols if c in raw.columns]]

        tables.append(raw)

    combined = pd.concat(tables, ignore_index=True)

    if faims:
        combined["cv"] = combined["dataset"].str.split("_").str[-1]
        combined["dataset"] = combined["dataset"].str.replace(r"_[0-9]*$", "", regex=True)
    else:
        combined["cv"] = "noFAIMS"

    return combined


# ═══════════════════════════════════════════════════════════════════════
# MS1 (feature) reader
# ═══════════════════════════════════════════════════════════════════════


def _read_ms1(file_path: Path, file_name: list[str]) -> pd.DataFrame:
    # Detect version from the first file
    first_path = file_path / file_name[0]
    _validate_file(first_path)

    with open(first_path) as fh:
        header_line = fh.readline()

    version = _detect_ms1_version(header_line)
    if version == "unknown":
        # Same refusal as the MS2 reader: an unrecognized header must not be
        # silently interpreted through the 1.5.4 column map.
        raise ValueError(
            f"Unrecognized TopPIC MS1 table version in {first_path}: no known "
            "version marker (Apex_time for 1.8, Time_apex for 1.5, "
            "Monoisotopic_mass for 1.7) found in the header."
        )
    col_map = _MS1_VERSION_MAP[version]

    tables = []
    for fname in file_name:
        full_path = file_path / fname
        _validate_file(full_path)

        raw = pd.read_csv(full_path, sep="\t")
        raw = raw.rename(columns=col_map, errors="ignore")
        raw["dataset"] = fname.replace("_ms1.feature", "")

        # Keep pipeline columns + any extra columns from newer versions
        keep_cols = list(_MS1_KEEP)
        for c in raw.columns:
            if c not in keep_cols and c.startswith("_"):
                keep_cols.append(c)
        raw = raw[[c for c in keep_cols if c in raw.columns]]

        tables.append(raw)

    combined = pd.concat(tables, ignore_index=True)
    combined["cv"] = "noFAIMS"
    return combined


# ═══════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════

_HEADER_MARKER = "Data file name"


def _find_header_line(path: Path) -> int:
    """Return the 0-based index of the line starting with ``Data file name``."""
    with open(path) as fh:
        for i, line in enumerate(fh):
            if line.startswith(_HEADER_MARKER):
                return i
    raise ValueError(f"Could not find header line ('{_HEADER_MARKER}') in {path}")


def _validate_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"TopPIC output file not found: {path}")


def _safe_uniprot_acc(acc: str) -> Optional[str]:
    if pd.isna(acc):
        return None
    return extract_uniprot_acc(acc)
