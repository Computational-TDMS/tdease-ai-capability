"""
Cluster proteoforms by mass and retention time using hierarchical
clustering.
"""

import numpy as np
import pandas as pd


def cluster(
    df: pd.DataFrame,
    errors: dict,
    method: str = "ward",
    height: float = 3.0,
    min_size: int = 2,
) -> pd.DataFrame:
    """Run hierarchical clustering on normalized mass and RT.

    Adds a ``cluster`` column (integer) per row.  Cluster 0 denotes
    "noise" — clusters with fewer than *min_size* members.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain ``gene``, ``recal_mass``, and ``rt_align``.
    errors : dict
        Output of :func:`toppic.alignment.mass_calibrate.calc_error`;
        must have keys ``ppm_sd`` and ``rt_sd``.
    method : str, default="ward"
        Agglomeration method for ``scipy.cluster.hierarchy.linkage``.
    height : float, default=3.0
        Height at which to cut the dendrogram (in standard-deviation
        units, since data are normalised).
    min_size : int, default=2
        Minimum cluster size; smaller clusters are relabelled as noise.

    Returns
    -------
    pd.DataFrame
        Clustered data with cluster assignment.  Genes with only one
        observation are dropped.
    """
    # scipy is only needed here. Importing it at module top forces every
    # `import toppic` (pipeline, CLI, metadata) to drag in scipy even when no
    # clustering runs, which broke import-level tests in environments without
    # scipy. Localise the heavy dep to the one call site that uses it.
    from scipy.cluster.hierarchy import cut_tree, linkage

    ppm_sd = errors["ppm_sd"]
    rt_sd = errors["rt_sd"]

    # Remove genes with a single observation (hclust needs >=2)
    gene_counts = df.groupby("gene")["recal_mass"].transform("size")
    result = df[gene_counts > 1].copy()

    if result.empty:
        result["cluster"] = 0
        return result

    # Normalise (guard against division by zero)
    ppm_denom = np.log10(1 + ppm_sd / 1e6)
    if ppm_denom == 0:
        ppm_denom = 1.0
    rt_denom = rt_sd if rt_sd > 0 else 1.0

    result["_norm_mass"] = np.log10(result["recal_mass"]) / ppm_denom
    result["_norm_rt"] = result["rt_align"] / rt_denom

    # Cluster per gene. R's stats::dist allows missing coordinates: each pair
    # is computed from the coordinates available for both observations, and
    # Euclidean distances are scaled by p / n_valid. Preserve that behavior
    # instead of turning an entire gene into noise when one RT is missing.
    clusters = []
    for gene, grp in result.groupby("gene"):
        X = grp[["_norm_rt", "_norm_mass"]].values
        if len(X) < 2:
            cl = np.zeros(len(X), dtype=int)
        else:
            if np.all(np.isfinite(X)):
                Z = linkage(X, method=method)
            else:
                # R's dist(..., method="euclidean") excludes missing
                # coordinate pairs and scales the squared sum by p / n.
                # Passing the condensed distances preserves the same input
                # semantics for scipy's hierarchical clustering.
                Z = linkage(_r_euclidean_condensed(X), method=method)
            cl = cut_tree(Z, height=height).flatten()
            # Convert from 0-based to 1-based numbering (like R's cutree)
            cl = _renumber_clusters(cl)

        noise = _find_noise(cl, min_size)
        cl[np.isin(cl, noise)] = 0
        clusters.append(
            pd.DataFrame({"gene": gene, "_idx": grp.index, "cluster": cl})
        )

    clust_df = pd.concat(clusters, ignore_index=True).set_index("_idx")
    result["cluster"] = clust_df["cluster"]
    result = result.drop(columns=["_norm_mass", "_norm_rt"])

    return result


def _renumber_clusters(cl: np.ndarray) -> np.ndarray:
    """Renumber clusters from 0-based to 1-based, preserving order."""
    uniq = np.unique(cl)
    mapping = {old: i + 1 for i, old in enumerate(uniq)}
    return np.array([mapping[c] for c in cl])


def _find_noise(cl: np.ndarray, min_size: int) -> np.ndarray:
    """Find cluster IDs with fewer than *min_size* members."""
    ids, counts = np.unique(cl, return_counts=True)
    # np.isin treats a set as one object, rather than a collection of IDs.
    return ids[counts < min_size]


def _r_euclidean_condensed(X: np.ndarray) -> np.ndarray:
    """Return R ``dist(..., method='euclidean')`` distances for missing data.

    R excludes coordinates missing in either member of a pair and scales the
    squared sum by ``p / n_valid``. Pairs with no jointly finite coordinate
    remain non-finite and are allowed to fail in scipy as they would in the
    native R clustering call.
    """
    n, p = X.shape
    out = np.empty(n * (n - 1) // 2, dtype=float)
    k = 0
    for i in range(n - 1):
        valid = np.isfinite(X[i + 1 :]) & np.isfinite(X[i])
        n_valid = valid.sum(axis=1)
        delta = np.where(valid, X[i + 1 :] - X[i], 0.0)
        sums = np.sum(delta * delta, axis=1)
        values = np.full(n - i - 1, np.nan, dtype=float)
        ok = n_valid > 0
        values[ok] = np.sqrt(sums[ok] * (p / n_valid[ok]))
        out[k : k + len(values)] = values
        k += len(values)
    return out
