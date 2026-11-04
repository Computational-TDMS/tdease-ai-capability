"""TopPICR-compatible quantitative analysis for top-down proteomics."""

from .pipeline import (
    PipelineConfig,
    PipelineInputs,
    PipelineResult,
    TopPICPipeline,
)

__version__ = "0.1.0"

__all__ = [
    "PipelineConfig",
    "PipelineInputs",
    "PipelineResult",
    "TopPICPipeline",
    "__version__",
]
