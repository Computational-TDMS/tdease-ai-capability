"""Exercise noise exclusion through clustering, merging and quantification."""
import numpy as np
import pandas as pd

from toppic.clustering.cluster import cluster
from toppic.clustering.merge import create_pcg
from toppic.features.match import match_features


def test_subthreshold_cluster_cannot_contribute_quantitative_features():
    errors = {"ppm_sd": 1.0, "rt_sd": 0.1}
    observations = pd.DataFrame({
        "gene": ["g"] * 3,
        "dataset": ["sample"] * 3,
        "cv": [0] * 3,
        "recal_mass": [1000.0, 1000.0, 1050.0],
        "rt_align": [10.0] * 3,
        "feature_intensity": [10.0, 20.0, 999.0],
    })
    clustered = cluster(observations, errors, min_size=2)
    assert clustered.loc[2, "cluster"] == 0
    merged = create_pcg(clustered, errors)
    assert merged.index.tolist() == [0, 1]
    ms1 = observations.rename(columns={"feature_intensity": "intensity"})
    intensities = match_features(merged, ms1, errors, n_jobs=1)
    assert len(intensities) == 1
    assert intensities.intensity.tolist() == [20.0]
    empty = create_pcg(cluster(observations, errors, min_size=4), errors)
    assert empty.empty and "pc_group" in empty.columns
    assert match_features(empty, ms1, errors, n_jobs=1).empty


def test_equal_size_merge_seeds_follow_native_r_order():
    # Native TopPICR groupate fixture: 21 clusters along a chain. The largest
    # seed consumes its neighbor first; equal-sized seeds then follow their
    # ascending IDs. An unstable sort changes which overlapping pairs merge.
    sizes = [2] * 20 + [3]
    observations = pd.DataFrame({
        "gene": "g",
        "cluster": np.repeat(np.arange(1, 22), sizes),
        "recal_mass": 1000.,
        "rt_align": np.repeat(10 + np.arange(21) * .2, sizes),
    })
    expected = [1, 1, 3, 3, 5, 5, 7, 7, 9, 9, 11, 11, 13, 13,
                15, 15, 17, 17, 19, 21, 21]
    merged = create_pcg(observations, {"ppm_sd": 1., "rt_sd": .1})
    np.testing.assert_array_equal(merged.pc_group, np.repeat(expected, sizes))


def test_missing_rt_uses_available_mass_like_native_r_dist():
    """A missing RT must not turn a whole Gene into cluster-0 noise."""
    errors = {"ppm_sd": 1.0, "rt_sd": 0.1}
    observations = pd.DataFrame({
        "gene": ["g"] * 3,
        "recal_mass": [1000.0, 1000.0, 1000.0],
        "rt_align": [10.0, np.nan, 10.0],
    })

    clustered = cluster(observations, errors, min_size=2)

    # R's stats::dist uses the finite mass coordinate for every pair and
    # therefore forms one three-member cluster.
    np.testing.assert_array_equal(clustered["cluster"].to_numpy(), [1, 1, 1])
