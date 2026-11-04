"""
Retention time alignment across datasets.

Aligns retention times from multiple datasets to a common reference
using LOESS regression on shared proteoforms.
"""

from typing import Any

import numpy as np
import pandas as pd

from toppic.clustering._loess import loess_fit, loess_predict


def find_ref_ds(df: pd.DataFrame) -> str:
    """Find the dataset with the highest number of unique proteoforms.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain ``dataset`` and ``proteoform``.

    Returns
    -------
    str
        Dataset name to use as the reference.
    """
    counts = df.groupby("dataset")["proteoform"].nunique()
    return counts.idxmax()


def form_model(
    df: pd.DataFrame,
    ref_ds: str,
    **kwargs,
) -> dict[str, Any]:
    """Build a LOESS model between each dataset and a reference dataset.

    Only proteoforms found in both the reference and the working dataset
    are used to fit the model.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain ``dataset``, ``feature_intensity``,
        ``feature_apex``, and ``proteoform``.
    ref_ds : str
        Reference dataset name.
    **kwargs
        Extra arguments passed to the LOESS fitting (see
        :func:`~toppic.clustering._loess.loess_fit`).

    Returns
    -------
    dict[str, Any]
        Keys are dataset names; values are LOESS model objects.
    """
    # Reference data: one row per proteoform (max intensity, min apex)
    x_ref = _pick_representative(df[df["dataset"] == ref_ds])

    models = {}
    for ds in df["dataset"].unique():
        x_cur = _pick_representative(df[df["dataset"] == ds])
        # Find shared proteoforms
        common = np.intersect1d(x_ref["proteoform"], x_cur["proteoform"])
        if len(common) < 2:
            models[ds] = None
            continue

        idx_ref = x_ref["proteoform"].isin(common)
        idx_cur = x_cur["proteoform"].isin(common)

        y = x_ref.loc[idx_ref, "feature_apex"].values
        x = x_cur.loc[idx_cur, "feature_apex"].values

        models[ds] = loess_fit(x, y, **kwargs)

    return models


def align_rt(
    df: pd.DataFrame,
    model: dict[str, Any],
    var_name: str = "feature_apex",
) -> pd.DataFrame:
    """Align retention times to a reference using precomputed LOESS models.

    Adds an ``rt_align`` column.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain ``dataset`` and the column specified by *var_name*.
    model : dict
        LOESS models keyed by dataset name (output of :func:`form_model`).
    var_name : str, default="feature_apex"
        Name of the retention-time column to align.

    Returns
    -------
    pd.DataFrame
    """
    # Only require at least one dataset to align; skip missing ones
    available = set(df["dataset"].unique()) & set(model.keys())
    if not available:
        raise ValueError(
            "No common datasets between model and df. "
            f"Model has {set(model.keys())}, df has {set(df['dataset'].unique())}."
        )

    result = df.copy()
    result["rt_align"] = np.nan
    for ds in result["dataset"].unique():
        mask = result["dataset"] == ds
        aligned = _apply_model(result.loc[mask], model, var_name)
        result.loc[mask, "rt_align"] = aligned

    # _apply_model handles absent models. Missing fitted predictions must
    # remain missing, matching R's no-extrapolation policy.
    result["rt_align"] = result["rt_align"].round(5)
    return result


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _pick_representative(grp: pd.DataFrame) -> pd.DataFrame:
    """Pick one representative retention time per proteoform.

    Within each proteoform, selects the row with the maximum feature
    intensity; ties are broken by the minimum feature apex. Mirrors R's
    ``slice_max(Feature intensity) %>% slice_min(Feature apex)``.
    """
    return (
        grp[["feature_intensity", "feature_apex", "proteoform"]]
        .drop_duplicates()
        .sort_values(
            ["feature_intensity", "feature_apex"], ascending=[False, True]
        )
        .groupby("proteoform", as_index=False)
        .first()
    )


def _apply_model(
    grp: pd.DataFrame,
    model: dict,
    var_name: str,
) -> np.ndarray:
    ds = grp["dataset"].iloc[0]
    mdl = model.get(ds)
    if mdl is None:
        return grp[var_name].values
    x_new = grp[var_name].values
    return loess_predict(mdl, x_new)
