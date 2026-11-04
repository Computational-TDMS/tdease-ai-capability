"""
Build an expression matrix (MSnSet equivalent) from matched features and
metadata.

The output is a dataclass with three tables: the expression matrix,
feature annotations, and sample annotations.
"""

import numpy as np
import pandas as pd

from toppic.pipeline.result import PipelineResult


class TopPICResult(PipelineResult):
    """Container analogous to Bioconductor's ``MSnSet``.

    Attributes
    ----------
    expression : pd.DataFrame
        Rows = features (``gene_pc_group``), columns = sample names,
        values = intensity.
    features : pd.DataFrame
        Feature metadata (index = feature ID).
    samples : pd.DataFrame
        Sample metadata (index = sample name).
    This compatibility class shares its fields with
    :class:`toppic.pipeline.PipelineResult`.
    """


def to_msnset(
    combined: pd.DataFrame,
    meta: pd.DataFrame | None = None,
    dataset_feature_data: pd.DataFrame | None = None,
    dataset_sample_data: pd.DataFrame | None = None,
    proteoform_def: list[str] | None = None,
    partitioning: list[str] | None = None,
) -> TopPICResult:
    """Convert matched features into an expression matrix.

    Parameters
    ----------
    combined : pd.DataFrame
        Output of :func:`~toppic.features.match.match_features`.
        Must contain ``dataset``, ``gene``, ``pc_group``, ``intensity``,
        and optionally ``cv``.
    meta : pd.DataFrame or None
        Output of :func:`~toppic.features.metadata.create_mdata`.
    dataset_feature_data : pd.DataFrame or None
        Maps ``dataset`` to additional feature-level columns
        (e.g. fraction).
    dataset_sample_data : pd.DataFrame or None
        Maps ``dataset`` to sample-level columns (e.g. ``sample_name``,
        ``group``).  If not provided, each dataset is its own sample.
    proteoform_def : list[str], optional
        Columns defining a unique proteoform (default ``["gene",
        "pc_group"]``).
    partitioning : list[str], optional
        Additional columns that split a proteoform into multiple
        features (e.g. ``cv`` for FAIMS, ``fraction``).

    Returns
    -------
    TopPICResult
    """
    if proteoform_def is None:
        proteoform_def = ["gene", "pc_group"]

    x = combined.copy()

    # Merge dataset feature data
    if dataset_feature_data is not None:
        x = x.merge(dataset_feature_data, on="dataset", how="left")

    # Sample data
    if dataset_sample_data is not None:
        x = x.merge(dataset_sample_data, on="dataset", how="left")
    else:
        x["sample_name"] = x["dataset"]
        dataset_sample_data = x[["dataset", "sample_name"]].drop_duplicates()

    # Build feature and proteoform IDs
    feature_def = proteoform_def + (partitioning or [])
    # Explicit column-by-column conversion to avoid pandas type issues
    feature_ids = x[feature_def[0]].astype(str)
    for col in feature_def[1:]:
        feature_ids = feature_ids + "_" + x[col].astype(str)
    x["feature_id"] = feature_ids

    proteoform_ids = x[proteoform_def[0]].astype(str)
    for col in proteoform_def[1:]:
        proteoform_ids = proteoform_ids + "_" + x[col].astype(str)
    x["proteoform_id"] = proteoform_ids

    # Check uniqueness
    dupes = x[["feature_id", "sample_name"]].duplicated()
    if dupes.any():
        raise ValueError(
            "Feature/sample combinations are not unique.  "
            "Check that fractions / CV values are included in the "
            "feature definition."
        )

    # Pivot to expression matrix
    expr = x.pivot_table(
        index="feature_id",
        columns="sample_name",
        values="intensity",
        aggfunc="max",
    ).fillna(0)

    expr = expr.astype(float)

    # Feature annotations
    feat = x.groupby(["feature_id", "proteoform_id"] + feature_def, as_index=False).agg(
        median_intensity=("intensity", "median"),
        count=("intensity", "count"),
    )

    if meta is not None:
        feat = feat.merge(meta, on=proteoform_def, how="left")

    feat = feat.set_index("feature_id")

    # Sample annotations
    samp = (
        dataset_sample_data.drop(columns=["dataset"])
        .drop_duplicates()
        .set_index("sample_name")
    )

    return TopPICResult(expression=expr, features=feat, samples=samp)
