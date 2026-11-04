# /// script
# requires-python = ">=3.11"
# dependencies = ["pandas", "numpy", "scipy", "statsmodels"]
# ///
"""TopPICR run report (skill entry point).

One pipeline command, ONE self-contained D3 report: the fixed default view set
for a pipeline output dir, packed into a single offline HTML (vendored d3,
chart-forge tokens; per-panel SVG/PNG export lives in the report page's Save
buttons). No PNG pathway, no plotly, no per-figure sidecar files, no CDN.

The former per-figure PNG subcommands are gone by design: the report covers
the run overview deterministically, and a focused single figure or a special
plot is chart-forge's job (system/skills/chart-forge/examples, e.g.
proteoform-map for a proteoform sequence view) or an in-place edit of the report's D3.

Subcommands:
  report    -- one self-contained HTML report from a pipeline output dir:
               PCA, volcano when --group-a/--group-b are given, completeness,
               sample correlation, expression heatmap, feature map, proteoform
               levels (Smith 2019), RT alignment, and mass recalibration.
               The sections are presented as tabs in one self-contained page.
               Per-section status for missing inputs; stdout is one JSON
               envelope ({report, figures, sections}).
  compare   -- compact ranked group-comparison JSON (bounded agent view).

Usage (system-assets prefix, quoted):
  "$TDEASE_SYSTEM_ASSETS/skills/toppicr-py/scripts/toppic_plot.py" \
    report --project runs/<run>/analysis/<slug> \
    [--group-a <regex> --group-b <regex>] [--out <report.html>]
  "$TDEASE_SYSTEM_ASSETS/skills/toppicr-py/scripts/toppic_plot.py" \
    compare --expression .../expression.csv --features .../features.csv \
    --group-a <regex> --group-b <regex> [--top 12]
"""

import argparse
import json
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import special

_HERE = Path(__file__).resolve().parent
_SKILLS = _HERE.parents[1]
for _path in (str(_SKILLS), str(_HERE)):
    if _path not in sys.path:
        sys.path.insert(0, _path)


def _sample_token(name) -> str:
    """Compact display label for a sample column (house convention)."""
    match = re.search(r"(?:^|_)([A-Za-z]+\d+)(?:_|$)", str(name))
    return match.group(1) if match else str(name)[:24]


def _group_of(name) -> str | None:
    """Letter prefix of a sample's group token ('WT' from '_WT02_'), else None."""
    match = re.search(r"(?:^|_)([A-Za-z]+)\d+(?:_|$)", str(name))
    return match.group(1) if match else None


def _order_samples(columns) -> list:
    """Biological-group blocks (first-appearance order), within-group order kept."""
    cols = list(columns)
    groups = [_group_of(c) for c in cols]
    if all(g is None for g in groups):
        return cols
    seen = list(dict.fromkeys(g for g in groups if g is not None))
    ordered = []
    for g in seen + [None]:
        ordered.extend(c for c, gg in zip(cols, groups) if gg == g)
    return ordered


def _group_blocks(ordered_cols) -> list[dict]:
    """Group-block boundaries for the matrix panel's design separators."""
    blocks = []
    for c in ordered_cols:
        g = _group_of(c) or "other"
        if blocks and blocks[-1]["label"] == g:
            blocks[-1]["n"] += 1
        else:
            blocks.append({"label": g, "n": 1})
    return blocks


def _clustered_group_order(corr: pd.DataFrame, ordered_cols: list) -> list:
    """Within each biological-group block, order samples by average linkage so
    clustering never masks the design; no groups -> full clustering."""
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import squareform

    def leaf_order(frame: pd.DataFrame) -> list:
        dist = ((1 - frame).fillna(1.0) + (1 - frame).fillna(1.0).T) / 2
        values = dist.to_numpy(copy=True)
        np.fill_diagonal(values, 0.0)
        values = (values + values.T) / 2
        return list(leaves_list(linkage(squareform(values), method="average")))

    groups = [_group_of(c) for c in ordered_cols]
    if all(g is None for g in groups):
        return [ordered_cols[i] for i in leaf_order(corr.loc[ordered_cols, ordered_cols])]
    result = []
    for g in dict.fromkeys(g for g in groups if g is not None):
        block = [c for c in ordered_cols if _group_of(c) == g]
        if len(block) > 1:
            sub = corr.loc[block, block]
            result.extend(block[i] for i in leaf_order(sub))
        else:
            result.extend(block)
    return result


# ---------------------------------------------------------------------------
# calibration data pipeline (pipeline stages, not plotting)
# ---------------------------------------------------------------------------


def _result_table(project: Path, manifest: dict, name: str):
    filename = manifest.get("outputs", {}).get(name)
    path = project / (filename or f"{name}.csv")
    if not filename and not path.is_file():
        return None
    return (pd.read_parquet(path) if path.suffix == ".parquet"
            else pd.read_csv(path, index_col=0 if name == "expression" else None))


# ---------------------------------------------------------------------------
# report panels: each returns the panel contract for _shared.report
# ---------------------------------------------------------------------------


def _panel_completeness(expression: pd.DataFrame) -> dict:
    present = (expression > 0).sum()
    return {
        "rows": [{"label": _sample_token(c), "value": int(v),
                  "hover": f"<b>{c}</b>: {v} / {expression.shape[0]} features detected"}
                 for c, v in present.items()],
        "ylabel": "Features detected", "valueLabels": True, "height": 300,
    }


def _panel_correlation(expression: pd.DataFrame, raw: bool) -> dict:
    data = expression if raw else np.log2(expression.replace(0, np.nan))
    corr = data.corr(method="pearson")
    ordered = _clustered_group_order(corr, _order_samples(corr.columns))
    corr = corr.loc[ordered, ordered]
    return {
        "rowLabels": [_sample_token(c) for c in corr.index],
        "colLabels": [_sample_token(c) for c in corr.columns],
        "values": [[None if pd.isna(v) else round(float(v), 4) for v in row]
                   for row in corr.values],
        "diverging": True, "lo": -1.0, "hi": 1.0,
        "colGroups": _group_blocks(ordered),
        "legendTitle": "Pearson r", "fmt": ".2f", "height": 420,
    }


def _panel_heatmap(expression: pd.DataFrame, top_n: int, raw: bool) -> dict:
    data = expression if raw else np.log2(expression.replace(0, np.nan))
    var = data.var(axis=1, skipna=True)
    top = data.assign(_var=var).sort_values("_var", ascending=False).head(top_n).drop(columns="_var")
    ordered = _order_samples(top.columns)
    top = top[ordered]
    values = top.astype(float)
    return {
        "rowLabels": [str(i)[:18] for i in top.index],
        "colLabels": [_sample_token(c) for c in top.columns],
        "values": [[None if pd.isna(v) else round(float(v), 3) for v in row]
                   for row in values.values],
        "diverging": False,
        "lo": float(np.nanmin(values.values)), "hi": float(np.nanmax(values.values)),
        "colGroups": _group_blocks(ordered),
        "legendTitle": "log2 intensity" if not raw else "intensity",
        "fmt": ".2f", "height": 620, "ylabel": "top-N features",
    }


def _panel_pca(expression: pd.DataFrame) -> dict:
    X = np.log2(expression.replace(0, np.nan))
    X = X.subtract(X.median(axis=0, skipna=True), axis=1)
    X = X.T.fillna(X.median(axis=1, skipna=True)).fillna(0).values.astype(float)
    Xc = X - X.mean(axis=0)
    U, S, _ = np.linalg.svd(Xc, full_matrices=False)
    scores = U[:, :2] * S[:2]
    var = (S[:2] ** 2) / (S ** 2).sum() * 100

    groups: dict[str, list[int]] = {}
    for i, sample in enumerate(expression.columns):
        groups.setdefault(_group_of(sample) or "other", []).append(i)
    detected = (expression > 0).sum(axis=0)
    points = []
    for i, sample in enumerate(expression.columns):
        group = _group_of(sample) or "other"
        points.append({
            "sample": str(sample), "group": group,
            "x": round(float(scores[i, 0]), 4), "y": round(float(scores[i, 1]), 4),
            "sizeValue": int(detected.iloc[i]),
            "hover": (f"<b>{sample}</b><br>{group}<br>"
                      f"PC1 {scores[i, 0]:+.3f} · PC2 {scores[i, 1]:+.3f}<br>"
                      f"detected features {int(detected.iloc[i]):,}"),
        })
    return {
        "points": points, "vlines": [{"v": 0}], "hlines": [{"v": 0}],
        "pc1Pct": round(float(var[0]), 1), "pc2Pct": round(float(var[1]), 1),
        "xlabel": f"PC1 ({var[0]:.1f}%)", "ylabel": f"PC2 ({var[1]:.1f}%)",
        "sizeLabel": "detected features", "height": 520,
    }


def _featuremap_sources(project: Path, manifest: dict) -> list[Path]:
    """Resolve the original TopFD MS1 feature files recorded by the run."""
    sources = []
    for entry in manifest.get("inputs", {}).get("ms1", []):
        raw = entry.get("path") if isinstance(entry, dict) else entry
        if not raw:
            continue
        candidate = Path(raw)
        if not candidate.is_file():
            candidate = project / "project" / "topfd" / candidate.name
        if candidate.is_file() and candidate.suffix == ".feature":
            sources.append(candidate)
    return list(dict.fromkeys(sources))


def _panel_featuremap(sources: list[Path]) -> dict:
    """Read every TopFD interval; do not replace an interval with a centroid."""
    keep = ["Feature_ID", "Mass", "Min_time", "Max_time", "Apex_time",
            "Apex_intensity", "Rep_charge", "Rep_average_mz"]
    rows = []
    used = []
    for source in sources:
        frame = pd.read_csv(source, sep="\t", usecols=keep)
        frame = frame.dropna(subset=["Rep_average_mz", "Min_time", "Max_time", "Apex_time"])
        for row in frame.itertuples(index=False):
            apex_intensity = float(row.Apex_intensity) if pd.notna(row.Apex_intensity) else 0.0
            rows.append({
                "mz": round(float(row.Rep_average_mz), 6),
                "rt0": round(float(row.Min_time), 4),
                "rt1": round(float(row.Max_time), 4),
                "apex": round(float(row.Apex_time), 4),
                "logI": round(float(np.log10(max(apex_intensity, 1.0))), 5),
                "q": int(row.Rep_charge) if pd.notna(row.Rep_charge) else None,
                "mass": round(float(row.Mass), 5) if pd.notna(row.Mass) else None,
                "featureId": str(row.Feature_ID),
                "dataset": _sample_token(source.stem.replace("_ms1", "")),
            })
        used.append(str(source))
    if not rows:
        raise ValueError("TopFD MS1 feature files contain no usable intervals")
    return {"features": rows, "mode": "interval", "n": len(rows), "sourceCount": len(used),
            "sourceFiles": used}


LEVEL_ORDER = ["1", "2A", "2B", "2C", "2D", "3A", "3B", "3C", "3D", "3E",
               "3F", "4A", "4B", "4C", "4D", "5"]


def _panel_levels(features: pd.DataFrame) -> dict:
    required = {"proteoform_level", "spectral_count"}
    missing = required - set(features.columns)
    if missing:
        raise ValueError(f"features table missing {sorted(missing)}; rerun the pipeline")
    counts = features.groupby(features["proteoform_level"].astype(str))["spectral_count"].sum()
    order = LEVEL_ORDER
    total = float(counts.sum())
    if total <= 0:
        raise ValueError("features table has no positive spectral_count")
    pct = counts.reindex(order, fill_value=0) / total * 100
    return {
        "rows": [{"label": lvl, "value": round(float(pct[lvl]), 2),
                  "hover": f"level {lvl}: {int(counts.get(lvl, 0)):,} PrSMs ({pct[lvl]:.1f}%)"}
                 for lvl in order],
        "xlabel": "Proteoform level (Smith 2019)", "ylabel": "% of PrSMs",
        "valueLabels": True, "height": 300,
    }


def _panel_rt_alignment(calib) -> dict:
    series = []
    # Bound SVG size only; the calibration artifact retains every observation.
    budget = max(1, 2000 // max(1, calib["dataset"].nunique()))
    for dataset, rows in calib.groupby("dataset"):
        total = len(rows)
        rows = rows.sort_values("feature_apex")
        if total > budget:
            rows = rows.iloc[np.linspace(0, total - 1, budget, dtype=int)]
        series.append({"label": f"{dataset} ({len(rows)}/{total} rows)", "points": [
            {"x": float(raw), "y": float(aligned - raw)}
            for raw, aligned in zip(rows["feature_apex"], rows["rt_align"])
            if np.isfinite(raw) and np.isfinite(aligned)]})
    return {
        "series": series, "hlines": [{"v": 0}],
        "xlabel": "Raw retention time (min)", "ylabel": "Applied RT correction (min)",
        "legend": True, "pointR": 2.6, "height": 420,
    }


def _panel_mass_error(calib) -> dict:
    ref = calib.loc[calib["mass_reference"] == True]
    adjusted = ref["adjusted_precursor_mass"].astype(float)
    ppm_before = 1e6 * (ref["precursor_mass"].astype(float) - adjusted) / adjusted
    ppm_after = 1e6 * (ref["recal_mass"].astype(float) - adjusted) / adjusted
    order = sorted(ref["dataset"].unique())
    rows = []
    for ds in order:
        before = ppm_before[ref["dataset"] == ds].abs().median()
        after = ppm_after[ref["dataset"] == ds].abs().median()
        rows.append({"label": _sample_token(ds), "series": "before",
                     "value": round(float(before), 3),
                     "hover": f"<b>{ds}</b> before: median |ppm| {before:.3f}"})
        rows.append({"label": _sample_token(ds), "series": "after",
                     "value": round(float(after), 3),
                     "hover": f"<b>{ds}</b> after: median |ppm| {after:.3f}"})
    return {
        "rows": rows, "series": ["before", "after"],
        "xlabel": "Sample", "ylabel": "median |ppm error|", "valueLabels": False,
        "height": 320,
    }


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values, clamped to 1.0; NaN stays NaN."""
    ok = ~np.isnan(pvals)
    q = np.full(len(pvals), np.nan)
    p = pvals[ok]
    order = np.argsort(p)
    adjusted = p[order] * len(p) / np.arange(1, len(p) + 1)
    q[np.where(ok)[0][order]] = np.minimum(np.minimum.accumulate(adjusted[::-1])[::-1], 1.0)
    return q


# ---------------------------------------------------------------------------
# Official downstream chain — Python port of the TopPICR paper's statistics
# (MSnSet.utils log2_zero_center -> normalizeByGlob medpolish ->
#  normalize_by_feature_property lowess; limma lmFit+eBayes+topTable BH).
#  Ported from MSnSet.utils and limma sources; see cmd_compare definitions.
# ---------------------------------------------------------------------------

def _masked_medians(e: np.ndarray, axis: int) -> np.ndarray:
    """NaN-ignoring medians along axis; all-NaN slices contribute 0 (no shift)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        m = np.nanmedian(e, axis=axis)
    return np.where(np.isfinite(m), m, 0.0)


def _median_polish(frame: pd.DataFrame, maxiter: int = 100, tol: float = 1e-8) -> pd.DataFrame:
    """Tukey median polish (na-aware); returns the additive model's residuals —
    the R ``medpolish``-based global normalization (normalizeByGlob)."""
    e = frame.to_numpy(dtype=float)
    ok = np.isfinite(e)
    e = np.where(ok, e - _masked_medians(e, None), np.nan)
    for _ in range(maxiter):
        delta = _masked_medians(e, axis=1)  # row effects first (R medpolish order)
        e = np.subtract(e, delta[:, None])
        delta = _masked_medians(e, axis=0)
        e = e - delta
        if np.nanmax(np.abs(delta)) < tol:
            break
    return pd.DataFrame(np.where(ok, e, np.nan), index=frame.index, columns=frame.columns)


def _remove_rt_trend(frame: pd.DataFrame, rt: pd.Series, frac: float = 0.5) -> pd.DataFrame:
    """Per-sample lowess trend vs retention time, subtracted in place — the
    ``lowess`` branch of normalize_by_feature_property. Fitted on each sample's
    observed features; the trend is interpolated to every observed feature."""
    from statsmodels.nonparametric.smoothers_lowess import lowess

    x_all = rt.to_numpy(dtype=float)
    out = frame.copy()
    for sample in frame.columns:
        y = frame[sample].to_numpy(dtype=float)
        sel = np.isfinite(y) & np.isfinite(x_all)
        if sel.sum() < 8:
            continue  # too sparse to define a trend
        fitted = lowess(y[sel], x_all[sel], frac=frac)
        trend = np.interp(x_all, fitted[:, 0], fitted[:, 1])
        out[sample] = np.where(np.isfinite(y), y - trend, np.nan)
    return out


def _logmdigamma(x: float) -> float:
    return float(np.log(x) - special.digamma(x))


def _trigamma_inverse(y: float) -> float:
    """Solve trigamma(x) = y; trigamma is decreasing on (0, inf)."""
    lo, hi = 1e-8, 1e12
    for _ in range(200):
        mid = float(np.sqrt(lo * hi))
        if special.polygamma(1, mid) > y:
            lo = mid
        else:
            hi = mid
    return float(np.sqrt(lo * hi))


def _fit_fdist(x: np.ndarray, df: np.ndarray) -> tuple[float, float]:
    """Moment estimation of a scaled F prior for residual variances — the limma
    fitFDist port (no trend, not robust). Returns (df.prior, var.prior)."""
    x = np.asarray(x, dtype=float)
    df = np.asarray(df, dtype=float)
    ok = np.isfinite(x) & (x > -1e-15) & np.isfinite(df) & (df > 1e-15)
    x, df = x[ok], df[ok]
    if len(x) == 0:
        return np.inf, 1.0
    if len(x) == 1:
        return 0.0, float(x[0])
    median = float(np.median(x))
    if median <= 0:
        median = 1.0
    x = np.maximum(x, 1e-5 * median)
    e = np.log(x) + np.array([_logmdigamma(v) for v in df / 2])
    emean = float(e.mean())
    evar = float(((e - emean) ** 2).sum() / (len(e) - 1))
    evar -= float(np.mean([special.polygamma(1, v) for v in df / 2]))
    if evar > 0:
        df_prior = 2 * _trigamma_inverse(evar)
        var_prior = float(np.exp(emean - _logmdigamma(df_prior / 2)))
    else:
        df_prior = np.inf
        var_prior = float(np.mean(x))  # R uses the arithmetic mean here (2017+)
    return df_prior, var_prior


def _moderated_two_group(beta: np.ndarray, unscaled_se: np.ndarray,
                         s2: np.ndarray, df: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float]:
    """limma eBayes (trend=FALSE, robust=FALSE) over per-feature two-group OLS
    fits. Returns (p-values, moderated t, df.prior, var.prior)."""
    from scipy import stats

    df_prior, var_prior = _fit_fdist(s2, df)
    if np.isinf(df_prior):
        s2_post = np.full(len(s2), var_prior)
    else:
        s2_post = (df * s2 + df_prior * var_prior) / (df + df_prior)
    t = beta / (unscaled_se * np.sqrt(s2_post))
    df_total = np.minimum(df + df_prior, float(df.sum())) if not np.isinf(df_prior) else np.full(len(df), float(df.sum()))
    p = 2 * stats.t.sf(np.abs(t), df_total)
    return p, t, df_prior, var_prior


def _independent_groups(columns, group_a: str, group_b: str):
    a_cols = [c for c in columns if re.search(group_a, c)]
    b_cols = [c for c in columns if re.search(group_b, c)]
    if not a_cols or not b_cols:
        raise ValueError(f"cannot split samples: {group_a!r} vs {group_b!r}")
    overlap = sorted(set(a_cols) & set(b_cols))
    if overlap:
        raise ValueError(f"independent comparison groups overlap: {', '.join(overlap)}")
    return a_cols, b_cols


def _panel_volcano(expression: pd.DataFrame, features: pd.DataFrame | None,
                   group_a: str, group_b: str,
                   fdr: float, group_a_name: str | None, group_b_name: str | None) -> dict:
    a_cols, b_cols = _independent_groups(expression.columns, group_a, group_b)

    def display(cols, fallback):
        counts: dict[str, int] = {}
        for c in cols:
            g = _group_of(c)
            if g:
                counts[g] = counts.get(g, 0) + 1
        return max(counts, key=counts.get) if counts else fallback

    a_name, b_name = group_a_name or display(a_cols, "Group A"), group_b_name or display(b_cols, "Group B")

    # The SAME official chain as `compare` (default): one statistics authority
    # for the report and the agent-facing comparison — no second Welch track.
    log_matrix = np.log2(expression.where(np.isfinite(expression) & (expression > 0)))
    analyzed = _median_polish(log_matrix.sub(log_matrix.median(axis=1), axis=0))
    rt_note = "RT trend not removed (no usable RT on features.csv)"
    if features is not None and "rt" in features.columns:
        # cmd_report loads features with a RangeIndex (feature_id is a column);
        # compare loads it indexed. Normalize before reindexing onto the
        # expression rows, or every RT value silently becomes NaN.
        indexed = (features.drop_duplicates("feature_id").set_index("feature_id")
                   if "feature_id" in features.columns else features)
        rt = indexed["rt"].reindex(expression.index)
        if rt is not None and np.isfinite(rt.to_numpy(dtype=float)).sum() >= 8:
            analyzed = _remove_rt_trend(analyzed, rt.astype(float))
            rt_note = "per-sample lowess vs RT, frac 0.5, subtracted"

    A = analyzed[a_cols]
    B = analyzed[b_cols]
    Av, Bv = A.to_numpy(float), B.to_numpy(float)
    log2fc = (B.mean(axis=1, skipna=True) - A.mean(axis=1, skipna=True)).astype(float)

    beta, unscaled, s2, df, positions = [], [], [], [], []
    for pos in range(len(expression)):
        a = Av[pos][~np.isnan(Av[pos])]
        b = Bv[pos][~np.isnan(Bv[pos])]
        if len(a) >= 2 and len(b) >= 2 and np.var(a) + np.var(b) > 0:
            beta.append(b.mean() - a.mean())
            unscaled.append(np.sqrt(1 / len(a) + 1 / len(b)))
            s2.append(((len(a) - 1) * np.var(a, ddof=1) + (len(b) - 1) * np.var(b, ddof=1))
                      / (len(a) + len(b) - 2))
            df.append(len(a) + len(b) - 2)
            positions.append(pos)
    pvals = np.full(len(expression), np.nan)
    if positions:
        pvals[np.array(positions)] = _moderated_two_group(
            np.array(beta), np.array(unscaled), np.array(s2), np.array(df))[0]
    q = _bh_fdr(pvals)

    det_a = (expression[a_cols] > 0).any(axis=1)
    det_b = (expression[b_cols] > 0).any(axis=1)
    co = det_a & det_b
    mean_abundance = expression[a_cols + b_cols].replace(0, np.nan).mean(axis=1)

    CAP = 12.0
    fc_cutoff = 1.0
    points = []
    significant = higher_in_b = higher_in_a = 0
    for pos, idx in enumerate(expression.index):
        if not co.iloc[pos] or np.isnan(q[pos]):
            continue
        y = min(-np.log10(q[pos]) if q[pos] > 0 else CAP, CAP)
        p_sig = q[pos] < fdr
        fc_sig = abs(log2fc.iloc[pos]) > fc_cutoff
        if p_sig and fc_sig:
            key = "significant"
            significant += 1
            if log2fc.iloc[pos] > 0:
                higher_in_b += 1
            else:
                higher_in_a += 1
        elif p_sig:
            key = "p-only"
        elif fc_sig:
            key = "fc-only"
        else:
            key = "ns"
        points.append({
            "feature": str(idx), "x": round(float(log2fc.iloc[pos]), 4), "y": round(float(y), 4),
            "fdr": float(q[pos]),
            "sizeValue": None if pd.isna(mean_abundance.iloc[pos]) else round(float(mean_abundance.iloc[pos]), 4),
            "state": key, "capped": bool(q[pos] == 0 or -np.log10(q[pos]) > CAP),
        })
    return {
        "points": points, "fdr": fdr, "fcCutoff": fc_cutoff,
        "higherIn": b_name, "lowerIn": a_name, "yCap": CAP,
        "xlabel": f"log2 fold change ({b_name} / {a_name})", "ylabel": "-log10 BH FDR",
        "sizeLabel": "mean abundance", "height": 540,
        "method": "official chain: median-polish normalization + RT-trend removal + "
                  "moderated t-test (limma eBayes port) + BH — same as `compare`",
        "rtTrend": rt_note,
        "testedCount": len(points),
        "significantCount": significant,
        "higherInGroupB": higher_in_b,
        "higherInGroupA": higher_in_a,
    }


# ---------------------------------------------------------------------------
# report command
# ---------------------------------------------------------------------------


def cmd_report(args) -> int:
    """One self-contained D3 report; per-section status, one JSON envelope."""
    from _shared.report import compile_report, envelope, panel_script

    project = Path(args.project)
    if not project.is_dir():
        print(f"error: project dir not found: {project}", file=sys.stderr)
        return 2

    manifest_path = project / "run_manifest.json"
    manifest = (json.loads(manifest_path.read_text(encoding="utf-8"))
                if manifest_path.is_file() else {})
    expression = _result_table(project, manifest, "expression")
    features = _result_table(project, manifest, "features")
    calib = calib_error = None
    try:
        calib = _result_table(project, manifest, "calibration")
    except (OSError, ValueError, ImportError) as exc:
        calib_error = exc

    def calibration():
        if calib_error is not None:
            raise calib_error
        return calib

    ms1_sources = _featuremap_sources(project, manifest)
    no_ms1 = None if ms1_sources else "no TopFD .feature inputs recorded in run manifest"

    def build(pid, title, caption, panel_name, skip_reason, build_data):
        if skip_reason:
            return {"id": pid, "title": title, "caption": caption,
                    "status": "skipped", "note": skip_reason}
        try:
            return {"id": pid, "title": title, "caption": caption,
                    "script": panel_script(panel_name),
                    "data": build_data(), "meta": {"source": str(project)}}
        except Exception as exc:  # noqa: BLE001 - isolate failures to their report section
            return {"id": pid, "title": title, "caption": caption,
                    "status": "error", "note": str(exc)}

    group_pair = bool(args.group_a and args.group_b)
    no_expr = None if expression is not None else "no expression.csv in project"
    no_features = None if features is not None else "no features.csv in project"
    no_calib = (None if calib is not None or calib_error is not None
                else "run has no saved calibration observations")

    panels = [
        build("pca", "Sample PCA",
              "PCA of log2 expression (median centering, feature-median "
              "imputation); replicates should cluster by biological group.",
              "pca", no_expr, lambda: _panel_pca(expression)),
        build("volcano", "Differential abundance",
              "log2 fold change vs BH FDR on the paper's official chain "
              "(median-polish normalization, RT-trend removal, moderated "
              "t-test) — the same statistics as `toppic_plot.py compare`; "
              f"points encode mean abundance by area, while opacity separates the "
              f"non-significant cloud from thresholded points (BH FDR < {args.fdr}, "
              "|log2FC| > 1). Presence-only features stay outside this "
              "statistical view.",
              "volcano",
              None if (expression is not None and group_pair) else "needs --group-a and --group-b",
              lambda: _panel_volcano(expression, features, args.group_a, args.group_b,
                                     args.fdr, args.group_a_name, args.group_b_name)),
        build("completeness", "Data completeness (MBR)",
              "Features detected per sample against the feature-space total; "
              "gaps are recovery (MBR), not proven biological absence.",
              "bars", no_expr, lambda: _panel_completeness(expression)),
        build("correlation", "Sample correlation",
              "Pearson correlation of log2 intensities; within-group blocks are "
              "clustered, so off-diagonal structure exposes group separation or "
              "swapped injections.",
              "matrix", no_expr, lambda: _panel_correlation(expression, args.raw)),
        build("heatmap", "Expression heatmap",
              f"Top-{args.top_n} most-variable features (log2); rows carry the "
              "dominant abundance structure of the run.",
              "matrix", no_expr, lambda: _panel_heatmap(expression, args.top_n, args.raw)),
        build("featuremap", "Feature map (LC-MS space)",
              "All TopFD MS1 features are shown as horizontal elution-window "
              "bars at average m/z; bar width is the observed RT interval, "
              "colour is apex intensity, and charge is available on hover.",
              "featuremap", no_ms1, lambda: _panel_featuremap(ms1_sources)),
        build("levels", "Proteoform definition levels (Smith 2019)",
              "Percentage of PrSMs by proteoform definition level; low levels "
              "quantify the definition ambiguity the identification carries.",
              "bars", no_features, lambda: _panel_levels(features)),
        build("rt-alignment", "Applied RT alignment",
              "Saved aligned minus raw RT for each pipeline observation; "
              "this is the applied correction, not an accuracy estimate.",
              "scatter", no_calib, lambda: _panel_rt_alignment(calibration())),
        build("mass-error", "Mass recalibration (median |ppm|)",
              "Median absolute ppm residual on the pipeline's mass-reference "
              "rows before and after calibration; not held-out accuracy.",
              "bars", no_calib, lambda: _panel_mass_error(calibration())),
    ]

    out = compile_report(
        panels, args.out if args.out else project / "qc_report.html",
        title="TopPICR run report",
        meta_rows=[("Project", str(project)),
                   ("Samples", expression.shape[1] if expression is not None else "—"),
                   ("Expression features", expression.shape[0] if expression is not None else "—"),
                   ("Annotated features", len(features) if features is not None else "—"),
                   ("Calibration reference", manifest.get("calibration", {}).get("reference_dataset", "not recorded")),
                   *([("Groups", f"{args.group_a} vs {args.group_b}")] if group_pair else [])])
    summary = envelope(out, panels)
    # Headline numbers the prose will cite, so no consumer has to mine the
    # 20 MB HTML for them (the P1 case run regexed PC1 variance out of the
    # report four times). Generic report facts, computed by the panels anyway.
    by_id = {p["id"]: p for p in panels}
    highlights: dict[str, dict] = {}
    pca_data = by_id.get("pca", {}).get("data") or {}
    if "pc1Pct" in pca_data:
        highlights["pca"] = {"pc1Pct": pca_data["pc1Pct"], "pc2Pct": pca_data["pc2Pct"]}
    comp_rows = (by_id.get("completeness", {}).get("data") or {}).get("rows")
    if comp_rows and expression is not None:
        values = sorted(int(row["value"]) for row in comp_rows)
        highlights["detectedPerSample"] = {
            "min": values[0],
            "median": values[len(values) // 2],
            "max": values[-1],
            "featureTotal": int(expression.shape[0]),
        }
    vol_data = by_id.get("volcano", {}).get("data") or {}
    if "testedCount" in vol_data:
        highlights["volcano"] = {
            "tested": vol_data["testedCount"],
            "significant": vol_data["significantCount"],
            "higherInGroupB": vol_data["higherInGroupB"],
            "higherInGroupA": vol_data["higherInGroupA"],
            "fdr": vol_data["fdr"], "absLog2FcCutoff": vol_data["fcCutoff"],
        }
    if highlights:
        summary["highlights"] = highlights
    summary["outputDir"] = str(project.resolve())
    if manifest_path.is_file():
        summary["manifest"] = str(manifest_path.resolve())
    print(json.dumps(summary, ensure_ascii=False))
    return 0


# ---------------------------------------------------------------------------
# compare command (bounded JSON, no figures; default = the paper's official
# downstream chain, see cmd_compare's method block)
# ---------------------------------------------------------------------------


def cmd_compare(args) -> int:
    from scipy import stats

    expression = pd.read_csv(args.expression, index_col=0)
    features = pd.read_csv(args.features).drop_duplicates("feature_id").set_index("feature_id")
    a_cols, b_cols = _independent_groups(expression.columns, args.group_a, args.group_b)
    # A missing n_unexpected_mods column must not read as "no feature carries
    # one": --require-unexpected over a schema without the column would
    # return an empty candidate list and present it as a finding.
    has_unexpected_column = "n_unexpected_mods" in features.columns
    if args.require_unexpected and not has_unexpected_column:
        raise SystemExit(
            "features.csv has no n_unexpected_mods column; --require-unexpected "
            "cannot be applied to this export (regenerate with a pipeline that "
            "records it, or drop the flag)")

    # Detection counts and eligibility come from the RAW matrix; the analyzed
    # values come from the requested normalization (official chain by default).
    log_matrix = np.log2(expression.where(np.isfinite(expression) & (expression > 0)))
    method = {"normalization": ["log2"], "test": None, "rtTrend": None}
    if args.normalize == "official":
        centered = log_matrix.sub(log_matrix.median(axis=1), axis=0)
        polished = _median_polish(centered)
        method["normalization"] += ["per-feature median centering (log2_zero_center)",
                                    "median polish residuals (normalizeByGlob medpolish)"]
        rt = features["rt"].reindex(expression.index) if "rt" in features.columns else None
        if rt is not None and np.isfinite(rt.to_numpy(dtype=float)).sum() >= 8:
            analyzed = _remove_rt_trend(polished, rt.astype(float))
            method["rtTrend"] = "per-sample lowess vs RT, frac 0.5, subtracted (normalize_by_feature_property lowess branch)"
        else:
            analyzed = polished
            method["rtTrend"] = "skipped: fewer than 8 features carry a finite RT"
    else:
        analyzed = log_matrix
        method["normalization"].append("none (raw log2 detected values)")

    fits = []
    rows = []
    for feature_id, values in analyzed.iterrows():
        raw = expression.loc[feature_id]
        n_a = int((raw[a_cols] > 0).sum())
        n_b = int((raw[b_cols] > 0).sum())
        if n_a < args.min_detected or n_b < args.min_detected:
            continue
        a = values[a_cols].dropna()
        b = values[b_cols].dropna()
        meta = features.loc[feature_id] if feature_id in features.index else None
        unexpected = (int(meta.get("n_unexpected_mods", 0))
                      if meta is not None and has_unexpected_column else 0)
        if args.require_unexpected and unexpected == 0:
            continue
        row = {
            "featureId": str(feature_id),
            "log2FoldChange": None,
            "log2MedianRatio": None,
            "groupADetected": n_a,
            "groupBDetected": n_b,
            "groupAMedian": None,
            "groupBMedian": None,
            "pValue": None,
            "gene": None if meta is None else str(meta.get("gene", "")),
            "uniprotAccession": None if meta is None else str(meta.get("uniprot_acc", "")),
            "proteoform": None if meta is None else str(meta.get("proteoform", "")),
            "unexpectedModificationCount": unexpected,
            "spectralCount": None if meta is None else int(meta.get("spectral_count", 0)),
        }
        if len(a) >= 2 and len(b) >= 2:
            row["log2FoldChange"] = float(b.mean() - a.mean())
            row["log2MedianRatio"] = float(b.median() - a.median())
            row["groupAMedian"] = float(a.median())
            row["groupBMedian"] = float(b.median())
        rows.append(row)
        if len(a) >= 2 and len(b) >= 2 and np.var(a) + np.var(b) > 0:
            fits.append((len(rows) - 1, a, b))

    if args.stats == "moderated":
        method["test"] = ("two-group OLS per feature on available samples (pooled variance, df = nA+nB-2), "
                          "then limma eBayes shrinkage (fitFDist moments, trend=FALSE, robust=FALSE; Python port) "
                          "+ BH — the paper's lmFit/eBayes/topTable chain")
        beta = np.array([b.mean() - a.mean() for _, a, b in fits])
        unscaled = np.array([np.sqrt(1 / len(a) + 1 / len(b)) for _, a, b in fits])
        s2 = np.array([((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2)
                       for _, a, b in fits])
        df = np.array([len(a) + len(b) - 2 for _, a, b in fits])
        p_values, t_values, df_prior, var_prior = _moderated_two_group(beta, unscaled, s2, df)
        method["priorDf"] = float(df_prior) if np.isfinite(df_prior) else None
        method["priorVar"] = float(var_prior)
        for (index, _, _), p in zip(fits, p_values):
            rows[index]["pValue"] = float(p)
    else:
        method["test"] = ("Welch t-test on the analyzed (normalized) log2 values per feature; "
                          "requires >= 2 detections per group and nonzero variance")
        for index, a, b in fits:
            rows[index]["pValue"] = float(stats.ttest_ind(a, b, equal_var=False).pvalue)

    p_indices = [i for i, row in enumerate(rows) if row["pValue"] is not None]
    if p_indices:
        q = _bh_fdr(np.array([rows[i]["pValue"] for i in p_indices]))
        for pos, i in enumerate(p_indices):
            rows[i]["bhFdr"] = float(q[pos])
    for row in rows:
        row.setdefault("bhFdr", None)

    rows.sort(key=lambda row: (
        -(row["log2FoldChange"] if args.prefer_positive and row["log2FoldChange"] is not None else
          (abs(row["log2FoldChange"]) if row["log2FoldChange"] is not None else 0.0)),
        row["bhFdr"] is None,
        row["bhFdr"] if row["bhFdr"] is not None else 1.0,
        -row["spectralCount"],
    ))
    selected = rows[:args.top]
    significant = [
        row for row in rows
        if row["bhFdr"] is not None and row["bhFdr"] < 0.05
        and row["log2FoldChange"] is not None and abs(row["log2FoldChange"]) > 1
    ]
    payload = {
        "groupA": {"pattern": args.group_a, "samples": a_cols},
        "groupB": {"pattern": args.group_b, "samples": b_cols},
        "eligibleFeatureCount": len(rows),
        "testedCount": sum(1 for row in rows if row["pValue"] is not None),
        "bhSignificantCount": len(significant),
        "enrichedInGroupA": sum(1 for row in significant if row["log2FoldChange"] < 0),
        "enrichedInGroupB": sum(1 for row in significant if row["log2FoldChange"] > 0),
        "method": method,
        "definitions": {
            "eligible": f">= {args.min_detected} positive detections in BOTH groups (raw counts)"
                        + (" AND >= 1 unexpected modification" if args.require_unexpected else ""),
            "log2FoldChange": "mean(analyzed groupB) - mean(analyzed groupA); positive = higher in groupB",
            "log2MedianRatio": "median(analyzed groupB) - median(analyzed groupA); the analyzed scale is log2, so this is a log-ratio; robust companion to log2FoldChange",
            "groupAMedian/groupBMedian": "medians of the analyzed (normalized) values, not raw intensities",
            "pValue": "see method.test",
            "bhFdr": "Benjamini-Hochberg over features with a pValue; null when untestable",
            "bhSignificantCount": "features with bhFdr < 0.05 and |log2FoldChange| > 1 (the counts below that key use this rule)",
            "provenance": "default normalization/test reproduce the TopPICR paper's downstream chain as a Python port "
                          "(MSnSet.utils log2_zero_center, normalizeByGlob medpolish, normalize_by_feature_property "
                          "lowess; limma eBayes BH) — state the port, not R limma identity, in publications",
        },
        "candidates": selected,
    }
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, separators=(",", ":")))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="toppic_plot",
        description="TopPICR run report (one self-contained D3 HTML) + bounded group-comparison JSON.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    rep = sub.add_parser("report",
                         help="one self-contained D3 report: tabbed default view "
                              "set, per-section status")
    rep.add_argument("--project", required=True,
                     help="pipeline output dir (expression.csv, features.csv, toppic/ tables)")
    rep.add_argument("--out", default=None,
                     help="report HTML path (default: <project>/qc_report.html)")
    rep.add_argument("--top-n", type=int, default=40, help="heatmap top-N most-variable features")
    rep.add_argument("--raw", action="store_true", help="no log2 transform (default is log2)")
    rep.add_argument("--group-a", default=None,
                     help="sample regex for group A; presence enables the volcano section")
    rep.add_argument("--group-b", default=None, help="sample regex for group B")
    rep.add_argument("--group-a-name", default=None, help="display name for group A")
    rep.add_argument("--group-b-name", default=None, help="display name for group B")
    rep.add_argument("--fdr", type=float, default=0.05, help="volcano FDR threshold")
    rep.set_defaults(func=cmd_report)

    comp = sub.add_parser("compare", help="compact ranked group comparison for agent consumption")
    comp.add_argument("--expression", required=True, help="pipeline expression.csv")
    comp.add_argument("--features", required=True, help="pipeline features.csv")
    comp.add_argument("--group-a", required=True, help="sample regex for baseline group")
    comp.add_argument("--group-b", required=True, help="sample regex for comparison group")
    comp.add_argument("--normalize", choices=["official", "none"], default="official",
                      help="official: the paper's chain (log2, feature centering, median polish, "
                           "RT lowess trend removal); none: raw log2 detected values")
    comp.add_argument("--stats", choices=["moderated", "welch"], default="moderated",
                      help="moderated: per-feature two-group OLS + limma eBayes shrinkage (paper chain); "
                           "welch: per-feature Welch t-test")
    comp.add_argument("--min-detected", type=int, default=2,
                      help="minimum positive observations required in each group")
    comp.add_argument("--require-unexpected", action="store_true",
                      help="retain only features carrying unexpected modifications")
    comp.add_argument("--top", type=int, default=12, help="maximum candidates returned")
    comp.add_argument("--prefer-positive", action="store_true",
                      help="rank features enriched in group B before depleted features")
    comp.add_argument("--out", help="optional JSON output path")
    comp.set_defaults(func=cmd_compare)

    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
