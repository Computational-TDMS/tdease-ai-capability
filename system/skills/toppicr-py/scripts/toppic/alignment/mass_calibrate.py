"""
Mass recalibration using per-dataset ppm errors.

Computes the mass measurement error (MME) in ppm for each dataset,
then recalibrates precursor masses by the median offset.
"""

import numpy as np
import pandas as pd


def _mad(x) -> float:
    """Median absolute deviation, scaled to match R's ``stats::mad``.

    ``mad(x) = median(|x - median(x)|) * 1.4826``. R is the reference for the
    ppm/RT error widths (they drive clustering normalisation, pcGroup merging
    and MS1 matching windows), so this must match R exactly rather than the
    previous mean-absolute-deviation estimate.
    """
    x = np.asarray(x, dtype=float)
    return float(np.median(np.abs(x - np.median(x))) * 1.4826)


def mass_reference_mask(df: pd.DataFrame) -> pd.Series:
    """Rows used by the TopPICR mass-error estimator."""
    return ((df["n_unexpected_mods"] == 0)
            & ((df["precursor_mass"] - df["adjusted_precursor_mass"]).abs() < 0.5))


def calc_error(df: pd.DataFrame, ref_ds: str) -> dict:
    """Compute ppm and retention-time error statistics per dataset.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain ``dataset``, ``precursor_mass``,
        ``adjusted_precursor_mass``, ``n_unexpected_mods``,
        ``rt_align``, ``proteoform``, ``feature_intensity``,
        and ``feature_apex``.
    ref_ds : str
        Name of the reference dataset.

    Returns
    -------
    dict
        ``ppm_error`` — DataFrame of per-dataset ppm stats.
        ``ppm_sd`` — median of per-dataset ppm standard deviations.
        ``rt_sd`` — median absolute deviation of RT error.
    """
    # ---- ppm error ----
    ppm_sub = df.loc[mass_reference_mask(df)].copy()
    ppm_sub["ppm_error"] = (
        1e6
        * (ppm_sub["precursor_mass"] - ppm_sub["adjusted_precursor_mass"])
        / ppm_sub["adjusted_precursor_mass"]
    )

    ppm_stats = (
        ppm_sub.groupby("dataset")["ppm_error"]
        .agg(ppm_sd=_mad, ppm_median="median")
        .reset_index()
    )

    # R: global ppm_sd = median of per-dataset mad() (`stats::median`, default
    # na.rm = FALSE). NaN must propagate — R has no hard fallback, and a silent
    # 1.0 would silently widen the downstream matching/merging windows.
    ppm_sd_global = ppm_stats["ppm_sd"].median(skipna=False)

    print(f"Global ppm SD: {ppm_sd_global:.4f}")

    # ---- RT error ----
    x_ref = df[df["dataset"] == ref_ds][
        ["rt_align", "proteoform", "feature_intensity", "feature_apex"]
    ].copy()
    x_ref = _pick_rep(x_ref)

    rt_sds = []
    for ds in df["dataset"].unique():
        if ds == ref_ds:
            continue
        x_cur = df[df["dataset"] == ds][
            ["rt_align", "proteoform", "feature_intensity", "feature_apex"]
        ].copy()
        x_cur = _pick_rep(x_cur)

        merged = x_ref.merge(x_cur, on="proteoform", suffixes=("_ref", "_cur"))
        errors = (
            _mad(merged["rt_align_cur"] - merged["rt_align_ref"])
            if len(merged) > 0
            else np.nan
        )
        rt_sds.append(errors)

    # R: global rt_sd = median of per-dataset mad(). NaN propagates naturally.
    rt_sd_global = np.median(rt_sds) if rt_sds else np.nan

    return {
        "ppm_error": ppm_stats,
        "ppm_sd": ppm_sd_global,
        "rt_sd": rt_sd_global,
    }


def recalibrate_mass(
    df: pd.DataFrame,
    errors: dict,
    var_name: str,
) -> pd.DataFrame:
    """Recalibrate a mass column using per-dataset median ppm errors.

    Adds a ``recal_mass`` column.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain ``dataset`` and the column named by *var_name*.
    errors : dict
        Output of :func:`calc_error`.

    Returns
    -------
    pd.DataFrame
    """
    # TopPICR uses an inner join after asserting that every dataset has a
    # calibration row.  A permissive left join would silently create NaN
    # masses and change clustering/matching semantics.
    expected = set(df["dataset"].dropna().unique())
    observed = set(errors["ppm_error"]["dataset"].dropna().unique())
    if expected != observed:
        missing = sorted(expected - observed)
        extra = sorted(observed - expected)
        raise ValueError(
            "The data sets present in x and errors['ppm_error'] do not match "
            f"(missing={missing}, extra={extra})."
        )
    result = df.merge(errors["ppm_error"], on="dataset", how="inner", validate="many_to_one")
    result["recal_mass"] = result[var_name] * (1 - result["ppm_median"] / 1e6)
    result = result.drop(columns=["ppm_sd", "ppm_median"])
    return result


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _pick_rep(g: pd.DataFrame) -> pd.DataFrame:
    """One representative row per proteoform.

    Mirrors R's ``slice_max(Feature intensity) %>% slice_min(Feature apex)``:
    among equal max intensities pick the minimum feature apex. Sorting on both
    keys and taking ``first`` gives the same representative.
    """
    return (
        g.sort_values(
            ["feature_intensity", "feature_apex"], ascending=[False, True]
        )
        .groupby("proteoform", as_index=False)
        .first()[["rt_align", "proteoform"]]
    )
