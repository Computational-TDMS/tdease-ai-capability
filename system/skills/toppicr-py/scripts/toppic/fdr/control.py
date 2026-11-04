"""
FDR control functions for TopPIC data.

Implements E-value cutoff search, multi-level FDR computation,
and decoy-protein removal.
"""

import warnings

import numpy as np
import pandas as pd


def find_evalue_cutoff(df: pd.DataFrame, fdr_threshold: float = 0.01) -> np.ndarray:
    """Find the E-value cutoff for each annotation type that achieves the
    desired FDR.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain ``ann_type``, ``e_value``, ``gene``, and ``is_decoy``.
    fdr_threshold : float, default=0.01
        Desired FDR (e.g. 0.01 = 1 %).

    Returns
    -------
    np.ndarray
        E-value cutoffs, one per unique annotation type, in the order
        returned by ``df['ann_type'].unique()``.
    """
    anns = df["ann_type"].unique()
    cutoffs = [_find_cutoff_for_ann(df, a, fdr_threshold) for a in anns]
    return np.array(cutoffs)


def _find_cutoff_for_ann(
    df: pd.DataFrame, ann_type: str, fdr_threshold: float
) -> float:
    """Binary-search style search for the best E-value cutoff for one
    annotation type.
    """
    sub = df[df["ann_type"] == ann_type][["e_value", "gene", "is_decoy"]].copy()
    if sub.empty:
        return np.nan

    min_ev = sub["e_value"].min()
    max_ev = sub["e_value"].max()

    min_fdr = _calc_fdr(sub, min_ev)
    max_fdr = _calc_fdr(sub, max_ev)

    if min_fdr > fdr_threshold:
        warnings.warn(
            f"FDR for {ann_type} cannot be below {fdr_threshold} "
            f"(minimum attainable: {min_fdr:.3f})"
        )
        return min_ev

    if max_fdr < fdr_threshold:
        return max_ev

    # Binary search using quantile-based narrowing (mirrors R logic)
    ordered = sub["e_value"].sort_values().values
    n = len(ordered)

    # Quantile-based index vector (0th, 25th, 50th, 100th)
    inds = np.floor(np.quantile(np.arange(n), [0, 0.25, 0.5, 0.75, 1.0])).astype(int)
    inds = np.clip(inds, 0, n - 1)
    inds = np.unique(inds)  # handle edge cases where quantiles overlap

    # Safe default so ``maximin`` is always defined (R's own loop leaves it
    # unassigned → NA when the group is too small or no quantile point is
    # strictly below the threshold; we must never crash on those inputs).
    maximin = _scan_max_below_fdr(sub, fdr_threshold, default=min_ev)

    while len(inds) >= 2 and (inds[-1] - inds[0]) >= 2:
        summary_evals = ordered[inds]
        fdrs = np.array([_calc_fdr(sub, ev) for ev in summary_evals])

        below = np.where(fdrs < fdr_threshold)[0]
        if len(below) == 0:
            break

        max_below = below.max()
        maximin = ordered[inds[max_below]]
        if max_below >= len(inds) - 1:
            # Hit the top — can't narrow further
            break

        # Narrow search range
        upper = max_below + 1
        new_inds = np.floor(
            np.quantile(np.arange(inds[max_below], inds[upper]), [0, 0.25, 0.5, 0.75, 1.0])
        ).astype(int)
        inds = np.clip(new_inds, 0, n - 1)
        inds = np.unique(inds)

    return maximin


def _scan_max_below_fdr(
    sub: pd.DataFrame, fdr_threshold: float, default: float
) -> float:
    """Largest observed E-value whose gene-level FDR is below the threshold.

    Exact fallback for the cases where the R-mirroring quantile search cannot
    run (too few rows, or no sampled point is strictly below the threshold):
    iterate the discrete E-value boundaries and take the largest one that
    still satisfies the FDR constraint.
    """
    evals = np.sort(sub["e_value"].unique())
    ok = [ev for ev in evals if _calc_fdr(sub, ev) < fdr_threshold]
    return ok[-1] if ok else default


def _calc_fdr(sub: pd.DataFrame, eval_cutoff: float) -> float:
    """Calculate the gene-level FDR at a given E-value cutoff."""
    filtered = sub[sub["e_value"] <= eval_cutoff]
    if filtered.empty:
        return 0.0
    gene_level = filtered[["gene", "is_decoy"]].drop_duplicates()
    return gene_level["is_decoy"].mean()


def apply_evalue_cutoff(df: pd.DataFrame, e_vals: np.ndarray) -> pd.DataFrame:
    """Filter rows to keep only those with E-value below the type-specific
    cutoff.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain ``ann_type`` and ``e_value``.
    e_vals : np.ndarray
        One cutoff per unique annotation type (in the same order as
        ``df['ann_type'].unique()``).

    Returns
    -------
    pd.DataFrame
    """
    anns = df["ann_type"].unique()
    masks = []
    for i, ann in enumerate(anns):
        masks.append((df["ann_type"] == ann) & (df["e_value"] <= e_vals[i]))
    combined = np.logical_or.reduce(masks)
    return df[combined].copy()


def compute_fdr(df: pd.DataFrame) -> dict:
    """Compute FDR at the PrSM, proteoform, sequence, accession, and gene
    levels.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain ``dataset``, ``scans``, ``is_decoy``, ``proteoform``,
        ``clean_seq``, ``uniprot_acc``, and ``gene``.

    Returns
    -------
    dict
        Keys are level names; each value is ``(fdr_pct, count)``.
    """
    levels = {}

    for name, cols in [
        ("PrSM", ["dataset", "scans", "is_decoy"]),
        ("Proteoform", ["proteoform", "is_decoy"]),
        ("Sequence", ["clean_seq", "is_decoy"]),
        ("UniProtAcc", ["uniprot_acc", "is_decoy"]),
        ("Gene", ["gene", "is_decoy"]),
    ]:
        tmp = df[cols].drop_duplicates()
        fdr_pct = round(tmp["is_decoy"].mean() * 100, 3)
        count = int((~tmp["is_decoy"]).sum())
        levels[name] = (fdr_pct, count)
        print(
            f"{name} level FDR: {fdr_pct} %\n"
            f"Number of unique {name}s: {count}\n"
        )

    return levels


def remove_decoys(df: pd.DataFrame) -> pd.DataFrame:
    """Remove decoy rows and drop the ``is_decoy`` column."""
    return df[~df["is_decoy"]].drop(columns=["is_decoy"]).copy()
