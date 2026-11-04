"""
Protein inference using Occam's-razor principle.

Finds the minimal set of UniProt accessions that explains the largest
number of observed amino-acid sequences.
"""

import pandas as pd


def infer_proteins(df: pd.DataFrame) -> pd.DataFrame:
    """Run Occam's-razor protein inference.

    Iteratively selects the accession that matches the most remaining
    sequences, then removes all sequences it explains, until no sequences
    are left.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain ``clean_seq``, ``uniprot_acc``, and ``ann_type``.
        Typically this is post-FDR data.

    Returns
    -------
    pd.DataFrame
        Subset of the input containing only rows belonging to the
        selected protein set.
    """
    # Build reduced table: unique clean_seq + uniprot_acc combos
    y = df[["clean_seq", "uniprot_acc", "ann_type"]].drop_duplicates()
    y = y.copy()
    y["ann_type"] = pd.Categorical(
        y["ann_type"],
        categories=["SwissProt", "VarSplic", "TrEMBL"],
        ordered=True,
    )
    y = y.sort_values(["ann_type", "uniprot_acc", "clean_seq"])

    selected = _occam_inference(y)
    # Keep all original rows that match the selected set
    keys = ["clean_seq", "uniprot_acc", "ann_type"]
    # R semi_join is a membership filter: preserve observation order and do
    # not multiply records when the selected table contains duplicate keys.
    keep = pd.MultiIndex.from_frame(df[keys]).isin(
        pd.MultiIndex.from_frame(selected[keys])
    )
    return df.loc[keep].copy()


def _occam_inference(y: pd.DataFrame) -> pd.DataFrame:
    """Core Occam's-razor loop."""
    res = []
    remaining = y.copy()
    while len(remaining) > 0:
        # Accession matching the most remaining sequences. R picks the winner
        # with ``x[, .N, by = UniProtAcc][which.max(N)]``; data.table's ``by=``
        # returns groups in *first-appearance* order, and ``infer_proteins``
        # sorted by ``[ann_type, uniprot_acc, clean_seq]`` — so on a count tie R
        # keeps the SwissProt > VarSplic > TrEMBL priority (then accession order
        # within a type). pandas ``groupby`` sorts by key alphabetically, so the
        # tie-break must be applied explicitly to reproduce R.
        counts = remaining.groupby("uniprot_acc").size().rename("_count")
        # remaining is sorted by ann_type, so group first() = highest priority.
        first_ann = (
            remaining.groupby("uniprot_acc")["ann_type"].first().rename("_ann")
        )
        cand = pd.concat([counts, first_ann], axis=1).reset_index()
        cand = cand.sort_values(
            ["_count", "_ann", "uniprot_acc"],
            ascending=[False, True, True],  # most seqs, then SwissProt first
        )
        top_acc = cand.iloc[0]["uniprot_acc"]

        top_rows = remaining[remaining["uniprot_acc"] == top_acc]
        res.append(top_rows)

        # Remove all sequences explained by this accession
        explained = set(top_rows["clean_seq"])
        remaining = remaining[~remaining["clean_seq"].isin(explained)]

    return pd.concat(res, ignore_index=True)
