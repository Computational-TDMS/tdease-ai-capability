"""A native stats::loess fixture exercises duplicates, outliers and bounds."""
import numpy as np
import pandas as pd

from toppic.alignment.rt_align import align_rt
from toppic.clustering._loess import loess_fit, loess_predict


def test_default_alignment_matches_native_r_and_preserves_missing_predictions():
    # Generated with R 4.6.1: loess(y ~ x), then predict(model, query).
    # Repeated x=3 with conflicting y detects duplicate removal and robust
    # reweighting; off-grid queries also detect direct-vs-interpolated fitting.
    x = np.r_[np.arange(20.), 3., 3., 7., 14.]
    y = np.sin(x / 3) + x / 2
    y[[3, 20, 21]] = [8., -2., 3.]
    query = np.array([-1., 0., .5, 3., 7.3, 10., 14., 18.5, 19., 20., np.nan])
    expected = np.array([
        np.nan, .12660350879314258, .59945791700309070,
        2.55004524352534645, 4.27611076201787377, 4.82148224640566525,
        6.08389890439296011, 9.09347998719272255, 9.52668243026136885,
        np.nan, np.nan,
    ])
    model = loess_fit(x, y)
    np.testing.assert_allclose(loess_predict(model, query), expected, rtol=0, atol=1e-12)
    observations = pd.DataFrame({"dataset": "sample", "feature_apex": query})
    aligned = align_rt(observations, {"sample": model})
    np.testing.assert_allclose(aligned.rt_align, np.round(expected, 5), rtol=0, atol=0)
