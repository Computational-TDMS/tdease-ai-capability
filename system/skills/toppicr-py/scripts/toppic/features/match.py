"""
Match MS1 (unidentified) features to MS2 (identified) clusters.

For each identified cluster, finds MS1 features that fall within the
specified mass and retention-time tolerances.
"""

import concurrent.futures

import numpy as np
import pandas as pd


def match_features(
    ms2: pd.DataFrame,
    ms1: pd.DataFrame,
    errors: dict,
    n_mme_sd: float = 3.0,
    n_rt_sd: float = 3.0,
    summary_fn: str = "max",
    n_jobs: int = -1,
) -> pd.DataFrame:
    """Match unidentified MS1 features to identified MS2 clusters.

    Parameters
    ----------
    ms2 : pd.DataFrame
        Identified-feature data (output of ``create_pcg``). Must contain
        ``dataset``, ``cv``, ``gene``, ``pc_group``, ``cluster``,
        ``recal_mass``, ``rt_align``.
    ms1 : pd.DataFrame
        Unidentified-feature data (output of :func:`~toppic.alignment.
        mass_calibrate.recalibrate_mass` applied to MS1). Must contain
        ``dataset``, ``cv``, ``intensity``, ``recal_mass``, ``rt_align``.
    errors : dict
        Error dict from :func:`~toppic.alignment.mass_calibrate.calc_error`.
    n_mme_sd : float, default=3.0
        Number of mass-error SDs for the ppm tolerance.
    n_rt_sd : float, default=3.0
        Number of RT-error SDs for the RT tolerance.
    summary_fn : str, default="max"
        Function name to summarise multiple intensities: ``"max"``,
        ``"sum"``, or ``"median"``.
    n_jobs : int, default=-1
        Number of parallel jobs (-1 = all cores).

    Returns
    -------
    pd.DataFrame
        Combined MS1 + MS2 intensity table with columns
        ``dataset``, ``cv``, ``gene``, ``pc_group``, ``intensity``.
    """
    ppm_cutoff = n_mme_sd * errors["ppm_sd"]
    rt_cutoff = n_rt_sd * errors["rt_sd"]

    summary = _get_summary_fn(summary_fn)

    # ---- Centroids from MS2 ----
    centroids_mass = (
        ms2[ms2["cluster"] != 0]
        .groupby(["gene", "cluster", "pc_group"])["recal_mass"]
        .median()
        .reset_index(name="medi_mass")
    )
    centroids_rt = (
        ms2[ms2["cluster"] != 0]
        .groupby(["gene", "pc_group"])["rt_align"]
        .median()
        .reset_index(name="medi_rt")
    )
    centroids = centroids_mass.merge(centroids_rt, on=["gene", "pc_group"])

    if centroids.empty:
        return pd.DataFrame(
            columns=["dataset", "cv", "gene", "pc_group", "intensity"]
        )

    # ---- Parallel matching ----
    tasks = list(centroids.iterrows())

    if n_jobs == 1 or len(tasks) <= 1:
        results = [
            _match_one(idx, row, ms1, ppm_cutoff, rt_cutoff)
            for idx, row in tasks
        ]
    else:
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=n_jobs if n_jobs > 1 else None
        ) as pool:
            results = list(
                pool.map(
                    lambda idx_row: _match_one(
                        idx_row[0], idx_row[1], ms1, ppm_cutoff, rt_cutoff
                    ),
                    tasks,
                )
            )
    results = [r for r in results if r is not None]

    if not results:
        return pd.DataFrame(
            columns=["dataset", "cv", "gene", "pc_group", "intensity"]
        )

    output = pd.concat(results, ignore_index=True)
    output = (
        output.groupby(["dataset", "cv", "gene", "cluster", "pc_group"])["intensity"]
        .apply(summary)
        .reset_index()
    )

    # ---- Combine MS2 intensities ----
    ms2_intensity = (
        ms2.groupby(["dataset", "cv", "gene", "cluster", "pc_group"])[
            "feature_intensity"
        ]
        .max()
        .reset_index()
    )

    fused = output.merge(
        ms2_intensity, on=["dataset", "cv", "gene", "cluster", "pc_group"], how="outer"
    )
    fused["intensity"] = fused[["intensity", "feature_intensity"]].max(
        axis=1, skipna=True
    )
    fused = fused.drop(columns=["feature_intensity", "cluster"])

    # Final: one intensity per dataset/cv/gene/pc_group
    fused = (
        fused.groupby(["dataset", "cv", "gene", "pc_group"])["intensity"]
        .apply(summary)
        .reset_index()
    )

    return fused


def _match_one(
    idx: int,
    row: pd.Series,
    ms1: pd.DataFrame,
    ppm_cutoff: float,
    rt_cutoff: float,
) -> pd.DataFrame | None:
    """Match one centroid against all ms1 features."""

    ppm_mask = (
        1e6
        * (ms1["recal_mass"] - row["medi_mass"]).abs()
        / row["medi_mass"]
    ) < ppm_cutoff
    rt_mask = (ms1["rt_align"] - row["medi_rt"]).abs() < rt_cutoff
    matched = ms1[ppm_mask & rt_mask]

    if matched.empty:
        return None

    return matched[["dataset", "cv", "intensity"]].assign(
        gene=row["gene"],
        cluster=row["cluster"],
        pc_group=row["pc_group"],
    )





def _get_summary_fn(name: str):
    # R's summary functions are called with `na.rm = TRUE`; use the
    # NaN-skipping numpy variants so a single missing intensity does not poison
    # an otherwise-valid cluster.
    if name == "max":
        return np.nanmax
    elif name == "sum":
        return np.nansum
    elif name == "median":
        return np.nanmedian
    raise ValueError(f"Unknown summary function: {name}")
