"""Explicit sample median normalization for downstream group comparisons."""
import numpy as np
import pandas as pd


def normalize_log2_medians(
    expression: pd.DataFrame,
    strata: pd.Series,
    eligible: pd.Series | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Equalize sample log2 medians independently within each feature stratum.

    Zero, negative and nonfinite intensities are absent observations. ``strata``
    and optional ``eligible`` must describe every expression row by feature ID.
    Eligibility defines the normalization population, independently of the
    comparison cohort. No filtering, imputation or significance test is implied.
    The returned factor table records the full normalization population.
    """
    if not expression.index.is_unique or not expression.columns.is_unique:
        raise ValueError("Expression feature and sample IDs must be unique")
    if not strata.index.is_unique:
        raise ValueError("Normalization stratum IDs must be unique")
    groups = strata.reindex(expression.index)
    if groups.isna().any():
        raise ValueError("Every expression feature needs a normalization stratum")
    if eligible is None:
        population = pd.Series(True, index=expression.index)
    else:
        if not eligible.index.is_unique:
            raise ValueError("Normalization eligibility IDs must be unique")
        population = eligible.reindex(expression.index)
        if population.isna().any() or not population.isin([True, False]).all():
            raise ValueError("Normalization eligibility must be boolean for every feature")
        population = population.astype(bool)
    log = np.log2(expression.where(np.isfinite(expression) & (expression > 0)))
    normalized = log.copy()
    records = []
    for group in pd.unique(groups):
        rows = groups.eq(group)
        selected = log.loc[rows & population]
        medians = selected.median(axis=0)
        if not np.isfinite(medians).all():
            raise ValueError(f"Normalization stratum {group!r} lacks positive reference observations in a sample")
        target = float(medians.median())
        normalized.loc[rows] = log.loc[rows] - medians + target
        for sample in expression.columns:
            records.append({"stratum": group, "sample": sample,
                "population_features": len(selected),
                "positive_reference_observations": int(selected[sample].notna().sum()),
                "sample_median_log2": float(medians[sample]),
                "target_median_log2": target,
                "log2_shift": target-float(medians[sample])})
    return normalized, pd.DataFrame(records)
