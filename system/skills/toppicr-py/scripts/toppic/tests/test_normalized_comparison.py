"""Normalization must precede estimation and respect organism strata."""
import numpy as np
import pandas as pd
from toppic.utils.normalization import normalize_log2_medians


def test_sample_loading_offsets_do_not_become_fold_changes():
    log = pd.DataFrame({"M": [9., 10., 11., 23., 24., 25.],
                        "WT": [5., 6., 7., 25., 26., 27.]},
                       index=list("abcdef"))
    groups = pd.Series(["plant"]*3+["bacterium"]*3,index=log.index)
    normalized, factors = normalize_log2_medians(2.**log,groups)
    np.testing.assert_allclose(normalized.M-normalized.WT,0,atol=1e-12)
    assert len(factors)==4
    assert set(factors.population_features)=={3}
    # A loading correction must preserve a feature's biological deviation;
    # missing observations cannot silently turn into positive measurements.
    expression=2.**log
    expression.loc["a","M"] *= 2
    expression.loc["b","WT"] = 0
    population=pd.Series(True,index=log.index)
    population.loc["a"] = False
    observed, factors = normalize_log2_medians(expression,groups,population)
    assert np.isnan(observed.loc["b","WT"])
    assert factors.loc[factors.stratum.eq("plant"),"population_features"].eq(2).all()
    # Removing a reference row is an explicit population change, not a fitted
    # offset chosen from the desired fold-change outcomes.
    np.testing.assert_allclose(observed.loc["a","M"]-observed.loc["a","WT"],1.5)
