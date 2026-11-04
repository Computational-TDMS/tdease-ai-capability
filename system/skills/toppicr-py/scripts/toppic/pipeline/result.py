"""Public result model returned by the shared pipeline service."""

from dataclasses import dataclass, field

import pandas as pd


@dataclass
class PipelineResult:
    """Expression matrix, annotations, and execution provenance."""

    expression: pd.DataFrame
    features: pd.DataFrame
    samples: pd.DataFrame
    stage_metrics: dict[str, dict[str, int | float]] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    manifest: dict[str, object] = field(default_factory=dict)
    calibration: pd.DataFrame = field(default_factory=pd.DataFrame)
