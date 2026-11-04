"""Missing RT/mass coordinates cluster by R's pairwise-available distance."""
import numpy as np
import pandas as pd
import pytest

from toppic.clustering.cluster import _r_euclidean_condensed, cluster

# Generated with R 4.6.1: stats::dist(cbind(NormRTalign=rt, NormRecalMass=mass))
# with mass[2] and mass[6] NA, then hclust(method="ward.D2") + cutree().
RT = np.array([10.0, 10.2, 10.4, 9.9, 50.0, 50.3, 50.2, 49.8])
MASS = np.array([3.0, 3.1, np.nan, 3.05, 3.2, 3.25, np.nan, 3.3])
R_DISTANCES = np.array([
    0.22360679774997838, 0.56568542494923857, 0.11180339887498908,
    40.00049999687503544, 40.30077542678304070, 56.85138520739842249,
    39.80113063720677502, 0.28284271247462051, 0.30413812651490996,
    39.80012562794242825, 40.10028054764703853, 56.56854249492380404,
    39.60050504728442888, 0.70710678118654757, 56.00285706997456714,
    56.42712113868649482, 56.28569978244919270, 55.72001435749994158,
    40.10028054764704564, 40.40049504647188883, 56.99280656363573883,
    39.90078320033330073, 0.30413812651490818, 0.28284271247462306,
    0.22360679774998138, 0.14142135623730148, 0.50249378105604448,
    0.56568542494924612,
])


def test_condensed_distances_match_native_r_pairwise_missing_scaling():
    np.testing.assert_allclose(
        _r_euclidean_condensed(np.column_stack([RT, MASS])),
        R_DISTANCES, rtol=0, atol=1e-12,
    )


def _observations():
    errors = {"ppm_sd": 1.0, "rt_sd": 1.0}
    denom = np.log10(1 + 1e-6)
    return errors, pd.DataFrame({
        "gene": "g",
        "recal_mass": 10.0 ** (MASS * denom),
        "rt_align": RT,
    })


def test_missing_mass_rows_cluster_with_their_rt_neighbours_like_r():
    errors, observations = _observations()
    clustered = cluster(observations, errors, method="ward", height=3.0, min_size=2)
    assert clustered.cluster.tolist() == [1, 1, 1, 1, 2, 2, 2, 2]


def test_fine_cut_marks_singletons_noise_exactly_like_r():
    errors, observations = _observations()
    clustered = cluster(observations, errors, method="ward", height=0.25, min_size=2)
    assert clustered.cluster.tolist() == [1, 0, 0, 1, 4, 5, 5, 4]


def test_row_with_no_finite_coordinate_fails_like_r_hclust():
    errors, observations = _observations()
    observations.loc[0, ["recal_mass", "rt_align"]] = [np.nan, np.nan]
    with pytest.raises(ValueError):
        cluster(observations, errors, method="ward", height=3.0, min_size=2)
