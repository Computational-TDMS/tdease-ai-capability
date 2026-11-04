"""
Create feature metadata (proteoform metadata / ``create_mdata``).

Detects colliding clusters that share mass/RT centroids across different
genes and computes spectral counts.
"""

import re

import numpy as np
import pandas as pd


def create_mdata(
    df: pd.DataFrame,
    errors: dict,
    n_mme_sd: float = 3.0,
    n_rt_sd: float = 3.0,
) -> pd.DataFrame:
    """Create feature-level metadata for the expression matrix.

    Detects "collisions" — clusters from different genes that share the
    same mass/RT centroid.  Also computes spectral counts and median
    mass/RT per cluster group.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain ``gene``, ``pc_group``, ``cluster``, ``uniprot_acc``,
        ``recal_mass``, ``rt_align``, ``proteoform``, ``first_aa``,
        ``last_aa``, ``prot_length``, ``ann_type``, ``e_value``,
        ``n_unexpected_mods``, ``feature_intensity``.
    errors : dict
        Output of :func:`~toppic.alignment.mass_calibrate.calc_error`.
    n_mme_sd : float, default=3.0
        Number of mass-error SDs for collision detection.
    n_rt_sd : float, default=3.0
        Number of RT-error SDs for collision detection.

    Returns
    -------
    pd.DataFrame
        One representative row per ``(gene, pc_group)`` (ties kept, matching
        R's ``slice_max(spectralCount) %>% slice_min(E-value)``).
    """
    ppm_cutoff = n_mme_sd * errors["ppm_sd"]
    rt_cutoff = n_rt_sd * errors["rt_sd"]

    # ---- Collision detection (per gene + cluster, R semantics) ----
    collisions = _find_collisions(df, ppm_cutoff, rt_cutoff)
    df = df.merge(
        collisions, on=["gene", "cluster", "pc_group"], how="left"
    )
    df["collision"] = df["collision"].fillna("-")

    df["_mod_mass"] = df["proteoform"].apply(_extract_mod_mass)

    # ---- Per (gene, pc_group, clean_seq, mod_mass) group stats ----
    # R: group_by(Gene, pcGroup, cleanSeq, ModMass) %>% mutate(
    #       spectralCount = n(), rt = median(RTalign, na.rm=TRUE),
    #       mass = median(RecalMass, na.rm=TRUE))
    grp = df.groupby(
        ["gene", "pc_group", "clean_seq", "_mod_mass"], dropna=False
    )
    df["spectral_count"] = grp["feature_intensity"].transform("count")
    df["rt"] = grp["rt_align"].transform("median")
    df["mass"] = grp["recal_mass"].transform("median")

    # ---- Representative: slice_max(spectralCount) %>% slice_min(E-value) ----
    # Both keep ties (R's slice_* default with_ties = TRUE).
    max_sc = df.groupby(["gene", "pc_group"])["spectral_count"].transform("max")
    best = df[df["spectral_count"] == max_sc]
    min_ev = best.groupby(["gene", "pc_group"])["e_value"].transform("min")
    best = best[best["e_value"] == min_ev].copy()

    best = best.drop(columns=["_mod_mass"], errors="ignore")

    # NOTE: spectral_count must be exported — the fragment-coverage and PTM
    # figures colour by it (R: plot_accession sums spectralCount per fragment).
    # It is the size of the (gene, pc_group, clean_seq, mod_mass) group that
    # this representative row was selected from.
    return best[
        [
            "gene",
            "pc_group",
            "collision",
            "uniprot_acc",
            "mass",
            "rt",
            "first_aa",
            "last_aa",
            "prot_length",
            "proteoform",
            "ann_type",
            "n_unexpected_mods",
            "spectral_count",
            "proteoform_level",
        ]
    ]


# ---------------------------------------------------------------------------
# Collision detection
# ---------------------------------------------------------------------------


def _find_collisions(
    df: pd.DataFrame, ppm_cutoff: float, rt_cutoff: float
) -> pd.DataFrame:
    """Find clusters (gene + cluster) whose centroids collide across genes.

    Mirrors R's ``find_collider``: the centroid is computed per
    ``(Gene, cluster)`` (excluding cluster 0), and a collision is recorded
    (as ``Gene_pcGroup``) when two clusters from *different* genes share a
    mass/RT centroid within the tolerance.

    Returns
    -------
    pd.DataFrame
        Columns ``gene``, ``cluster``, ``pc_group``, ``collision``.
    """
    centroids = (
        df[df["cluster"] != 0]
        .groupby(["gene", "cluster"])
        .agg(
            cntr_mass=("recal_mass", "median"),
            cntr_rt=("rt_align", "median"),
            pc_group=("pc_group", "first"),
        )
        .reset_index()
        .assign(gcc=lambda x: x["gene"] + "_" + x["cluster"].astype(str))
    )

    centroids["collision"] = "-"
    n = len(centroids)

    for i in range(n):
        for j in range(i + 1, n):
            rt_ok = (
                abs(centroids["cntr_rt"].iloc[i] - centroids["cntr_rt"].iloc[j])
                < rt_cutoff
            )
            if not rt_ok:
                continue

            diff_ppm = (
                abs(
                    centroids["cntr_mass"].iloc[i]
                    - centroids["cntr_mass"].iloc[j]
                )
                / centroids["cntr_mass"].iloc[i]
                * 1e6
            )

            if diff_ppm < ppm_cutoff and centroids["gene"].iloc[i] != centroids["gene"].iloc[j]:
                centroids.loc[centroids.index[i], "collision"] = (
                    f"{centroids['gene'].iloc[j]}_{centroids['pc_group'].iloc[j]}"
                )
                centroids.loc[centroids.index[j], "collision"] = (
                    f"{centroids['gene'].iloc[i]}_{centroids['pc_group'].iloc[i]}"
                )

    return centroids[["gene", "cluster", "pc_group", "collision"]]


# ---------------------------------------------------------------------------
# ModMass extraction
# ---------------------------------------------------------------------------


def _extract_mod_mass(proteoform: str) -> float | None:
    """Extract numeric modification mass from a proteoform string.

    Removes named modifications (e.g. ``[Acetyl]``) and isolates the
    remaining numeric mass.
    """
    s = proteoform
    # Remove named mods like [Acetyl], [Plus1Oxy], etc.
    s = re.sub(r"\[[a-zA-Z]+[0-9]*.*?[a-zA-Z]*\]", "", s)
    # Remove amino-acid letters
    s = re.sub(r"[a-zA-Z]", "", s)
    # Remove leading/trailing periods
    s = s.strip(".")
    # Remove parens and brackets
    s = s.replace("(", "").replace(")", "").replace("[", "").replace("]", "")
    if not s.strip():
        return None
    try:
        return round(float(s), 0)
    except ValueError:
        return None
