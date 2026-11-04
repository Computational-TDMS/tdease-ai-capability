"""
Proteoform-level classification based on the Smith et al. 2019 system.

Assigns each proteoform a level (1, 2A, 2B, …, 5) based on confidence
metrics and modification evidence.
"""

import re

import numpy as np
import pandas as pd


def set_pf_level(df: pd.DataFrame) -> pd.DataFrame:
    """Assign proteoform levels according to the Smith 2019 classification.

    The level is derived from ``mi_score``, ``n_unexpected_mods``,
    ``acc_map``, ``special_amino_acids``, and a PTM score computed from
    the proteoform string.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain columns ``proteoform``, ``mi_score``,
        ``n_unexpected_mods``, ``acc_map``, ``special_amino_acids``.

    Returns
    -------
    pd.DataFrame
        Input with an added ``proteoform_level`` column.
    """
    x = df.copy()

    # Convert special_amino_acids: "-" -> NaN
    x["special_amino_acids"] = x["special_amino_acids"].replace("-", np.nan)

    # Compute PTMScore (extracts the "score" from the proteoform string)
    x["_ptm_score"] = x["proteoform"].apply(_compute_ptm_score)

    # Apply the classification rules
    levels = _classify(x)
    x["proteoform_level"] = levels
    x = x.drop(columns=["_ptm_score"])

    return x


# ---------------------------------------------------------------------------
# PTMScore computation
# ---------------------------------------------------------------------------


def _compute_ptm_score(pform: str) -> str:
    """Extract the PTM score string from a proteoform.

    Strips modification annotations, carbamidomethyl groups, digits,
    and punctuation, leaving only the parenthesised amino-acid content
    that indicates modified residues.
    """
    s = pform

    # Remove N-terminal single-residue acetylation annotations
    # e.g. .(K)[Acetyl] → .
    s = re.sub(r"\.\(([A-J,L-Z])\)\[Acetyl\]", "", s)

    # Collapse double carbamidomethyl: ((C)[CAM](C)[CAM]) -> (C)
    s = re.sub(
        r"\(\(C\)\[Carbamidomethyl\]\(C\)\[Carbamidomethyl\]\)",
        "(C)",
        s,
    )
    # Single double: ((C)[CAM]) -> (C)
    s = re.sub(r"\(\(C\)\[Carbamidomethyl\]\)", "(C)", s)
    # Regular carbamidomethyl: (C)[CAM] -> (C)  (we want to keep the (C) for scoring)
    # Actually the R code: (C)[Carbamidomethyl] → ""
    # Wait, re-reading R code: it removes `(C)[Carbamidomethyl]` entirely
    s = re.sub(r"\(C\)\[Carbamidomethyl\]", "", s)

    # Remove all bracketed annotations [...]
    s = re.sub(r"\[.*?\]", "", s)

    # Clean up double-parenthesised M
    s = re.sub(r"\(\(M", "(M", s)
    s = re.sub(r"M\)\)", "M)", s)

    # Remove digits
    s = re.sub(r"\d+", "", s)

    # Remove remaining brackets and punctuation-like characters
    s = s.replace("[", "").replace(".", "")
    s = s.replace("[", "").replace("-", "").replace(".", "")

    # Remove content between ) and ( i.e. everything between a closing
    # and the next opening paren
    s = re.sub(r"(?<=\)).*(?=\()", "", s)

    # Collapse )(
    s = s.replace(")(", "")

    # Extract inner content of outermost parentheses
    m = re.search(r"\((.+)\)", s)
    if m:
        s = m.group(1)
    else:
        s = ""

    # Remove any remaining punctuation
    s = re.sub(r"[^\w]", "", s)
    s = re.sub(r"_", "", s)

    return s


# ---------------------------------------------------------------------------
# Classification rules (Smith et al. 2019)
# ---------------------------------------------------------------------------


def _classify(df: pd.DataFrame) -> pd.Series:
    """Apply the 30+ case_when rules to assign proteoform levels."""

    mi = df["mi_score"].astype(str)
    unexp = df["n_unexpected_mods"]
    acc = df["acc_map"]
    ptm_len = df["_ptm_score"].fillna("").apply(len)
    sa = df["special_amino_acids"]

    has_mi = mi != "-"
    has_mi_str = mi.str.contains(r"\[\]", na=False)
    sa_missing = sa.isna()

    # Build conditions in order (R uses case_when = first match wins)
    level = pd.Series(index=df.index, dtype=str)

    # Level 1: MIScore = "-" & 0 unexpected & AccMap=1 & no PTMScore & no SAA
    _assign(level, (~has_mi) & (unexp == 0) & (acc == 1) & (ptm_len == 0) & sa_missing, "1")

    # Level 5: MIScore contains "[]"
    _assign(level, has_mi_str, "5")

    # Level 1: has MIScore, 0 unexpected, AccMap=1, PTMScore len=1, no SAA
    _assign(level, has_mi & (unexp == 0) & (acc == 1) & (ptm_len == 1) & sa_missing, "1")

    # Level 1: has MIScore, 0 unexpected, AccMap=1, PTMScore len=2, no SAA
    _assign(level, has_mi & (unexp == 0) & (acc == 1) & (ptm_len == 2) & sa_missing, "1")

    # Level 2A: has MISCORE, 0 unexpected, AccMap=1, PTMScore len>1, no SAA
    _assign(level, has_mi & (unexp == 0) & (acc == 1) & (ptm_len > 1) & sa_missing, "2A")

    # Level 2A: has MIScore, 1 unexpected, AccMap=1, PTMScore len>2, no SAA
    _assign(level, has_mi & (unexp == 1) & (acc == 1) & (ptm_len > 2) & sa_missing, "2A")

    # Level 2B: no MIScore, 1 unexpected, AccMap=1, PTMScore len=1, no SAA
    _assign(level, (~has_mi) & (unexp == 1) & (acc == 1) & (ptm_len == 1) & sa_missing, "2B")

    # Level 2B: no MIScore, 2 unexpected, AccMap=1, PTMScore len=2, no SAA
    _assign(level, (~has_mi) & (unexp == 2) & (acc == 1) & (ptm_len == 2) & sa_missing, "2B")

    # Level 2C: has MIScore, 0 unexpected, AccMap=1, has SAA
    _assign(level, has_mi & (unexp == 0) & (acc == 1) & (~sa_missing), "2C")

    # Level 2C: has MIScore, 1 unexpected, AccMap=1, PTMScore len=2, has SAA
    _assign(level, has_mi & (unexp == 1) & (acc == 1) & (ptm_len == 2) & (~sa_missing), "2C")

    # Level 2C: no MIScore, 0 unexpected, AccMap=1, no PTMScore, has SAA
    _assign(level, (~has_mi) & (unexp == 0) & (acc == 1) & (ptm_len == 0) & (~sa_missing), "2C")

    # Level 2D: has MIScore, 0 unexpected, AccMap>1, PTMScore len=1, no SAA
    _assign(level, has_mi & (unexp == 0) & (acc > 1) & (ptm_len == 1) & sa_missing, "2D")

    # Level 2D: has MIScore, 0 unexpected, AccMap>1, no PTMScore, no SAA
    _assign(level, has_mi & (unexp == 0) & (acc > 1) & (ptm_len == 0) & sa_missing, "2D")

    # Level 2D: has MIScore, 1 unexpected, AccMap>1, PTMScore len=2, no SAA
    _assign(level, has_mi & (unexp == 1) & (acc > 1) & (ptm_len == 2) & sa_missing, "2D")

    # Level 2D: no MIScore, 0 unexpected, AccMap>1, no PTMScore, no SAA
    _assign(level, (~has_mi) & (unexp == 0) & (acc > 1) & (ptm_len == 0) & sa_missing, "2D")

    # Level 3A: no MIScore, 1 unexpected, AccMap=1, PTMScore len>1, no SAA
    _assign(level, (~has_mi) & (unexp == 1) & (acc == 1) & (ptm_len > 1) & sa_missing, "3A")

    # Level 3B: has MIScore, 0 unexpected, AccMap=1, PTMScore len>1, has SAA
    _assign(level, has_mi & (unexp == 0) & (acc == 1) & (ptm_len > 1) & (~sa_missing), "3B")

    # Level 3B: has MIScore, 1 unexpected, AccMap=1, PTMScore len>2, has SAA
    _assign(level, has_mi & (unexp == 1) & (acc == 1) & (ptm_len > 2) & (~sa_missing), "3B")

    # Level 3C: no MIScore, 1 unexpected, AccMap=1, PTMScore len=1, has SAA
    _assign(level, (~has_mi) & (unexp == 1) & (acc == 1) & (ptm_len == 1) & (~sa_missing), "3C")

    # Level 3C: no MIScore, 2 unexpected, AccMap=1, PTMScore len=2, has SAA
    _assign(level, (~has_mi) & (unexp == 2) & (acc == 1) & (ptm_len == 2) & (~sa_missing), "3C")

    # Level 3D: has MIScore, 0 unexpected, AccMap>1, PTMScore len>1, no SAA
    _assign(level, has_mi & (unexp == 0) & (acc > 1) & (ptm_len > 1) & sa_missing, "3D")

    # Level 3D: has MIScore, 1 unexpected, AccMap>1, PTMScore len>2, no SAA
    _assign(level, has_mi & (unexp == 1) & (acc > 1) & (ptm_len > 2) & sa_missing, "3D")

    # Level 3E: no MIScore, 1 unexpected, AccMap>1, PTMScore len=1, no SAA
    _assign(level, (~has_mi) & (unexp == 1) & (acc > 1) & (ptm_len == 1) & sa_missing, "3E")

    # Level 3E: no MIScore, 2 unexpected, AccMap>1, PTMScore len=2, no SAA
    _assign(level, (~has_mi) & (unexp == 2) & (acc > 1) & (ptm_len == 2) & sa_missing, "3E")

    # Level 3F: has MIScore, 0 unexpected, AccMap>1, PTMScore len=1, has SAA
    _assign(level, has_mi & (unexp == 0) & (acc > 1) & (ptm_len == 1) & (~sa_missing), "3F")

    # Level 3F: has MIScore, 1 unexpected, AccMap>1, PTMScore len=2, has SAA
    _assign(level, has_mi & (unexp == 1) & (acc > 1) & (ptm_len == 2) & (~sa_missing), "3F")

    # Level 3F: no MIScore, 0 unexpected, AccMap>1, no PTMScore, has SAA
    _assign(level, (~has_mi) & (unexp == 0) & (acc > 1) & (ptm_len == 0) & (~sa_missing), "3F")

    # Level 4A: no MIScore, 1 unexpected, AccMap=1, PTMScore len>1, has SAA
    _assign(level, (~has_mi) & (unexp == 1) & (acc == 1) & (ptm_len > 1) & (~sa_missing), "4A")

    # Level 4A: no MIScore, 2 unexpected, AccMap=1, PTMScore len>2, has SAA
    _assign(level, (~has_mi) & (unexp == 2) & (acc == 1) & (ptm_len > 2) & (~sa_missing), "4A")

    # Level 4B: no MIScore, 1 unexpected, AccMap>1, PTMScore len>1, no SAA
    _assign(level, (~has_mi) & (unexp == 1) & (acc > 1) & (ptm_len > 1) & sa_missing, "4B")

    # Level 4B: no MIScore, 2 unexpected, AccMap>1, PTMScore len>2, no SAA
    _assign(level, (~has_mi) & (unexp == 2) & (acc > 1) & (ptm_len > 2) & sa_missing, "4B")

    # Level 4C: has MIScore, 0 unexpected, AccMap>1, PTMScore len>1, has SAA
    _assign(level, has_mi & (unexp == 0) & (acc > 1) & (ptm_len > 1) & (~sa_missing), "4C")

    # Level 4C: has MIScore, 1 unexpected, AccMap>1, PTMScore len>2, has SAA
    _assign(level, has_mi & (unexp == 1) & (acc > 1) & (ptm_len > 2) & (~sa_missing), "4C")

    # Level 4D: no MIScore, 1 unexpected, AccMap>1, PTMScore len=1, has SAA
    _assign(level, (~has_mi) & (unexp == 1) & (acc > 1) & (ptm_len == 1) & (~sa_missing), "4D")

    # Level 4D: no MIScore, 2 unexpected, AccMap>1, PTMScore len=2, has SAA
    _assign(level, (~has_mi) & (unexp == 2) & (acc > 1) & (ptm_len == 2) & (~sa_missing), "4D")

    # Level 5: no MIScore, 1 unexpected, AccMap>1, PTMScore len>1, has SAA
    _assign(level, (~has_mi) & (unexp == 1) & (acc > 1) & (ptm_len > 1) & (~sa_missing), "5")

    # Fallback for any unmatched rows
    level = level.fillna("unknown")

    return level


def _assign(series: pd.Series, condition, label: str) -> None:
    """Assign *label* where *condition* is True and no label is set yet."""
    series[(series.isna()) & (condition)] = label
