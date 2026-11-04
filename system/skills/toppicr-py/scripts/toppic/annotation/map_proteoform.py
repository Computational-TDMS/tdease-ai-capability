"""
Map proteoform sequences to protein entries in a FASTA database.

Extracts a "clean" amino-acid sequence from each proteoform string,
then locates that sequence within the matching protein entry from the
FASTA database to compute positional metadata.

Uses a hybrid strategy:
1. Primary: inner-join on ``uniprot_acc`` (fast, covers ~95%)
2. Fallback: for unmatched clean sequences, search the full FASTA
   database by sequence content (like R's ``augment_annotation``).
"""

import re
from math import floor, log10
from pathlib import Path

import numpy as np
import pandas as pd

from toppic.utils.fasta import (
    classify_annotation,
    extract_uniprot_acc,
    read_fasta,
)


def map_proteoform(
    df: pd.DataFrame,
    fasta_path: str | Path,
    fasta_name: str,
) -> pd.DataFrame:
    """Map proteoforms to their parent proteins in a FASTA database.

    Adds columns: ``clean_seq``, ``ann_type``, ``first_aa``, ``last_aa``,
    ``prot_length``, ``percent_coverage``.

    Normal and decoy entries are handled separately.  Unmatched rows
    are dropped only if there is genuinely no match in the database.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain ``proteoform``, ``uniprot_acc``, and ``is_decoy``.
    fasta_path : str | Path
        Path to the directory containing the FASTA file.
    fasta_name : str
        Name of the FASTA file (e.g. ``"human_mouse.fasta.gz"``).

    Returns
    -------
    pd.DataFrame
    """
    fasta_path = Path(fasta_path) / fasta_name
    records = read_fasta(fasta_path)
    db = _fasta_to_df(records)

    db_norm = db[~db["accession"].apply(lambda x: x.startswith("DECOY"))].copy()
    db_decoy = db[db["accession"].apply(lambda x: x.startswith("DECOY"))].copy()

    df = df.copy()
    df["clean_seq"] = df["proteoform"].apply(_clean_sequence)

    norm = df[~df["is_decoy"]].copy()
    decoy = df[df["is_decoy"]].copy()

    if not norm.empty:
        norm = _augment(norm, db_norm)
    if not decoy.empty:
        decoy = _augment(decoy, db_decoy)

    result = pd.concat([norm, decoy], ignore_index=True)
    return result


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _fasta_to_df(records: dict[str, str]) -> pd.DataFrame:
    rows = []
    for acc, seq in records.items():
        rows.append(
            {
                "accession": acc,
                "sequence": seq,
                "prot_length": len(seq),
                "uniprot_acc": extract_uniprot_acc(acc),
                "ann_type": classify_annotation(acc),
            }
        )
    return pd.DataFrame(rows)


def _clean_sequence(proteoform: str) -> str:
    """Extract the clean amino-acid sequence from a proteoform string."""
    if not isinstance(proteoform, str) or not proteoform:
        return ""
    s = re.sub(r"\[.+?\]", "", proteoform)
    s = re.sub(r"[()]", "", s)
    s = re.sub(r"^[A-Z]?\.(.+)\.[A-Z]?$", r"\1", s)
    return s


def _signif(x: float, digits: int = 2) -> float:
    """Significant-figure rounding matching R's ``signif(x, digits)``.

    R reports percent coverage as ``signif(..., 2)`` (2 significant figures);
    the historical Python used ``round(..., 2)`` (2 decimal places), which
    differs whenever the value is not in ``[1, 100)`` at 2 significant digits
    (e.g. 12.345 -> R 12, Python 12.35).
    """
    if not np.isfinite(x) or x == 0:
        return 0.0

    exponent = floor(log10(abs(x)))
    factor = 10.0 ** (digits - 1 - exponent)
    return np.round(x * factor) / factor


def _augment(df: pd.DataFrame, db: pd.DataFrame) -> pd.DataFrame:
    """Hybrid mapping: accession join then sequence-content fallback."""
    # Strategy 1: inner-join on uniprot_acc, locate clean_seq
    merged = df.merge(db, on="uniprot_acc", how="inner", suffixes=("", "_db"))
    merged[["first_aa", "last_aa"]] = merged.apply(
        lambda r: pd.Series(_str_locate(r["sequence"], r["clean_seq"])),
        axis=1,
    )
    # Compute coverage only for rows with valid positions. R uses
    # `signif(100*(lastAA-firstAA+1)/protLength, 2)` — two significant figures,
    # not two decimal places.
    merged["percent_coverage"] = (
        (merged["last_aa"].fillna(0) - merged["first_aa"].fillna(0) + 1)
        / merged["prot_length"] * 100
    ).map(lambda v: _signif(v))
    merged.loc[merged["first_aa"].isna(), "percent_coverage"] = 0.0

    # Rows that matched both accession AND position
    ok = merged[merged["first_aa"].notna()].copy()
    # Rows that matched accession but clean_seq NOT found in the protein
    no_pos = merged[merged["first_aa"].isna()].copy()

    # Strategy 2: for rows with no position match, try sequence-content search
    if not no_pos.empty:
        seq_fallback = _augment_by_sequence(no_pos, db)
        if seq_fallback is not None and not seq_fallback.empty:
            ok = pd.concat([ok, seq_fallback], ignore_index=True)
            # Rows the sequence fallback did NOT find are kept with NaN
            # positions, matching R's behaviour of retaining accession-matched
            # rows whose clean sequence cannot be located (firstAA/lastAA NA).
            # Dropping them (the historical behaviour) silently removes PrSMs
            # R keeps, which is the known PrSM-count divergence.
            found_keys = set(seq_fallback.index)
            no_pos = no_pos.loc[~no_pos.index.isin(found_keys)]

    if not no_pos.empty:
        # Keep accession-matched rows with unresolvable positions: NA first_aa/
        # last_aa, 0% coverage (R's map_proteoform keeps them as NA rows).
        no_pos = no_pos.copy()
        no_pos["percent_coverage"] = 0.0
        ok = pd.concat([ok, no_pos], ignore_index=True)

    ok = ok.drop(columns=["accession", "sequence", "sequence_db"], errors="ignore")
    return ok


def _augment_by_sequence(
    no_pos_df: pd.DataFrame, db: pd.DataFrame
) -> pd.DataFrame | None:
    """Search ALL protein sequences for matching clean sequences (fallback).

    This mirrors R's ``augment_annotation`` behaviour which uses
    ``grepl(peptide, sequence)`` across the entire FASTA.
    """
    rows = []
    for _, row in no_pos_df.iterrows():
        cs = row["clean_seq"]
        if not cs:
            continue
        # Search all db sequences for this clean_seq
        hits = db[db["sequence"].str.contains(cs, na=False, regex=False)]
        if hits.empty:
            continue
        # Pick the first match
        hit = hits.iloc[0]
        pos = _str_locate(hit["sequence"], cs)
        if pos is None:
            continue
        rows.append(
            {
                **row.drop(["sequence"], errors="ignore").to_dict(),
                "uniprot_acc": hit["uniprot_acc"],
                "ann_type": hit["ann_type"],
                "prot_length": hit["prot_length"],
                "first_aa": pos[0],
                "last_aa": pos[1],
                "percent_coverage": _signif(
                    (pos[1] - pos[0] + 1) / hit["prot_length"] * 100
                ),
            }
        )
    if not rows:
        return None
    return pd.DataFrame(rows)


def _str_locate(full_seq: str, sub_seq: str):
    """Return 1-based (start, end) or (NaN, NaN)."""
    if pd.isna(full_seq) or pd.isna(sub_seq) or not sub_seq:
        return pd.NA, pd.NA
    idx = full_seq.find(sub_seq)
    if idx == -1:
        return np.nan, np.nan
    return idx + 1, idx + len(sub_seq)
