"""
Merge clusters into proteoform cluster groups (pcGroups).

Implements the Metroid algorithm: starting from the largest cluster,
merge nearby clusters within a mass/RT envelope.
"""

import numpy as np
import pandas as pd

# Carbon-13 / Carbon-12 mass difference in Da
_C12C13 = 1.0033548378


def create_pcg(
    df: pd.DataFrame,
    errors: dict,
    n_mme_sd: float = 3.0,
    n_da: float = 4.0,
    n_rt_sd: float = 3.0,
) -> pd.DataFrame:
    """Merge clusters within each gene into proteoform cluster groups.

    Adds a ``pc_group`` column.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain ``gene``, ``cluster``, ``recal_mass``, and
        ``rt_align``.  Cluster 0 (noise) is skipped.
    errors : dict
        Error dict from :func:`~toppic.alignment.mass_calibrate.calc_error`.
    n_mme_sd : float, default=3.0
        Number of mass-error standard deviations for the ppm cutoff.
    n_da : float, default=4.0
        Number of Daltons for the isotopic-envelope envelope.
    n_rt_sd : float, default=3.0
        Number of RT-error standard deviations for the RT envelope.

    Returns
    -------
    pd.DataFrame
        Data with ``pc_group`` column added.
    """
    ppm_cutoff = n_mme_sd * errors["ppm_sd"]
    rt_sd = errors["rt_sd"]

    grouped = df[df["cluster"] != 0].copy()
    if grouped.empty:
        grouped["pc_group"] = grouped["cluster"]
        return grouped

    result = []
    for gene, grp in grouped.groupby("gene"):
        grp["pc_group"] = _groupate(grp, ppm_cutoff, n_da, rt_sd, n_rt_sd)
        result.append(grp)

    # Noise must also be excluded from downstream MS2 intensity fallback and
    # representative metadata, as in the R create_pcg reference function.
    return pd.concat(result, ignore_index=False)


def _groupate(
    x: pd.DataFrame,
    ppm_cutoff: float,
    n_da: float,
    rt_sd: float,
    n_rt_sd: float,
) -> np.ndarray:
    """Metroid-style cluster merging within one gene.

    Parameters
    ----------
    x : pd.DataFrame
        Data for a single gene (must have ``cluster``, ``recal_mass``,
        ``rt_align``, ``gene``).
    """
    # Build centroid table
    ridley = (
        x.groupby("cluster")
        .agg(
            n_members=("recal_mass", "size"),
            cntr_mass=("recal_mass", "median"),
            cntr_rt=("rt_align", "median"),
        )
        .reset_index()
        # R arrange is stable: equal-sized clusters retain ascending cluster
        # order from group_by. Seed order changes this greedy merge's members.
        .sort_values("n_members", ascending=False, kind="stable")
    )

    if len(ridley) == 1:
        return np.full(len(x), ridley["cluster"].values[0])

    samus = []

    while True:
        if len(ridley) == 0:
            break

        metroid = ridley.iloc[[0]]
        ridley = ridley.iloc[1:]

        if len(ridley) == 0:
            samus.append(metroid.assign(pcg=metroid["cluster"].values[0]))
            break

        # Find clusters within the RT envelope
        ridley_candidates = ridley.copy()
        ridley_candidates["diff_rt"] = (
            metroid["cntr_rt"].values[0] - ridley_candidates["cntr_rt"]
        ).abs()

        ridley_candidates = ridley_candidates[
            ridley_candidates["diff_rt"] < n_rt_sd * rt_sd
        ]

        if len(ridley_candidates) > 0:
            ridley_candidates["diff_mass"] = (
                metroid["cntr_mass"].values[0] - ridley_candidates["cntr_mass"]
            ).abs()
            ridley_candidates = ridley_candidates[
                ridley_candidates["diff_mass"] < n_da * _C12C13
            ]

        if len(ridley_candidates) > 0:
            ridley_candidates["diff_ppm"] = (
                (
                    metroid["cntr_mass"].values[0]
                    - ridley_candidates["cntr_mass"]
                    - _C12C13
                    * np.round(
                        metroid["cntr_mass"].values[0]
                        - ridley_candidates["cntr_mass"]
                    )
                )
                / metroid["cntr_mass"].values[0]
                * 1e6
            ).abs()
            ridley_candidates = ridley_candidates[
                ridley_candidates["diff_ppm"] < ppm_cutoff
            ]

        if len(ridley_candidates) == 0:
            samus.append(metroid.assign(pcg=metroid["cluster"].values[0]))
        else:
            to_remove = ridley_candidates["cluster"].values
            samus.append(
                pd.concat(
                    [
                        metroid.assign(pcg=metroid["cluster"].values[0]),
                        ridley_candidates.assign(pcg=metroid["cluster"].values[0]),
                    ]
                )
            )
            ridley = ridley[~ridley["cluster"].isin(to_remove)]

            # If ridley is now empty, we're done
            if len(ridley) == 0:
                break

    result = pd.concat(samus, ignore_index=True)
    # Merge back
    merged = x[["cluster"]].merge(
        result[["cluster", "pcg"]], on="cluster", how="left"
    )
    return merged["pcg"].fillna(x["cluster"]).values
