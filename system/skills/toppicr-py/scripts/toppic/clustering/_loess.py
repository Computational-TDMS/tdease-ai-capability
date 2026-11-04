"""One-dimensional Netlib LOESS with R's fitting and prediction defaults."""

import numpy as np


def loess_fit(
    x: np.ndarray,
    y: np.ndarray,
    frac: float = 0.75,
    degree: int = 2,
    iterations: int = 4,
    family: str = "gaussian",
    surface: str = "interpolate",
) -> dict:
    """Fit LOESS without removing repeated predictors or changing weights.

    scikit-misc exposes the Netlib LOESS routines also underlying R's
    ``stats::loess``. The Gaussian default uses least squares; robust
    iterations apply only when ``family='symmetric'`` is explicitly selected.
    The default surface uses the fitted interpolation tree, rather than a
    separate local polynomial calculation for every prediction.
    """
    from skmisc.loess import loess

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or y.shape != x.shape or len(x) < 2:
        raise ValueError("LOESS requires paired one-dimensional observations.")
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("LOESS training coordinates must be finite.")

    fitted = loess(
        x, y, span=frac, degree=degree, iterations=iterations,
        family=family, surface=surface,
    )
    fitted.fit()
    return {"fitted": fitted, "x_min": x.min(), "x_max": x.max(),
            "surface": surface}


def loess_predict(model: dict, x_new: np.ndarray) -> np.ndarray:
    """Predict with R's missing-value policy, including no extrapolation.

    ``predict.loess`` returns NA outside the observed predictor range for an
    interpolated surface. Keep those positions missing so clustering and
    feature matching do not silently use an unaligned retention time.
    """
    query = np.asarray(x_new, dtype=float)
    valid = np.isfinite(query)
    if model["surface"] == "interpolate":
        valid &= (query >= model["x_min"]) & (query <= model["x_max"])
    prediction = np.full(query.shape, np.nan)
    if valid.any():
        prediction[valid] = model["fitted"].predict(query[valid]).values
    return prediction
