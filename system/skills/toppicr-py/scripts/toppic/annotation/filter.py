"""
Filtering utilities for TopPIC data.

Includes filters based on feature count and gene assignment.
"""

import pandas as pd


def filter_by_count(
    df: pd.DataFrame,
    count_within: list[str],
    count_column: str,
    threshold: int,
) -> pd.DataFrame:
    """Filter rows where *count_column* appears fewer than *threshold* times.

    Parameters
    ----------
    df : pd.DataFrame
        Input data.
    count_within : list[str]
        Columns to group by when counting.
    count_column : str
        The column whose distinct values are counted.
    threshold : int
        Minimum count required to retain a row.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame.
    """
    # Deduplicate so count_column isn't listed twice in the selection
    cols = list(dict.fromkeys(count_within + [count_column]))
    counts = (
        df[cols]
        .drop_duplicates()
        .groupby(count_column)
        .size()
    )
    keep = counts[counts >= threshold].index
    return df[df[count_column].isin(keep)].copy()


def rm_false_gene(df: pd.DataFrame) -> pd.DataFrame:
    """For each feature, keep only the gene(s) with the smallest E-value.

    When a feature matches multiple genes, this function resolves the
    ambiguity by selecting the gene(s) with the best (lowest) E-value.
    If multiple genes have the same E-value, all are kept (matching R's
    ``dplyr::slice_min`` behaviour — no arbitrary tie-breaking).

    Parameters
    ----------
    df : pd.DataFrame
        Must contain ``dataset``, ``cv``, ``feature_apex``,
        ``feature_intensity``, ``gene``, and ``e_value`` columns.

    Returns
    -------
    pd.DataFrame
        Rows with non-optimal gene assignments removed.
    """
    group_cols = ["dataset", "cv", "feature_apex", "feature_intensity"]
    sub = df[group_cols + ["gene", "e_value"]].drop_duplicates()
    # Per-feature minimal e-value, then keep every gene tied at that minimum
    # (matches R's `slice_min` with_ties=TRUE). Avoid `groupby(...).apply`,
    # which in pandas 2.x drops the grouping keys from the result.
    min_ev = sub.groupby(group_cols, as_index=False)["e_value"].transform("min")
    best = sub[sub["e_value"] == min_ev][group_cols + ["gene"]]
    keys = group_cols + ["gene"]
    keep = pd.MultiIndex.from_frame(df[keys]).isin(pd.MultiIndex.from_frame(best))
    return df.loc[keep].copy()
