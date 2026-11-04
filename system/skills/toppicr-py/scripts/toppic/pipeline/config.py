"""Validated, framework-independent pipeline configuration."""

from dataclasses import asdict, dataclass

from .errors import PipelineConfigError


@dataclass(frozen=True)
class PipelineConfig:
    """Parameters shared by the Python API, CLI, and web interface."""

    fdr_threshold: float = 0.01
    occurrence_threshold: int = 2
    cluster_height: float = 3.0
    cluster_min_size: int = 2
    cluster_method: str = "ward"
    mass_error_sd: float = 3.0
    rt_error_sd: float = 3.0
    n_jobs: int = 1

    # scipy.cluster.hierarchy.linkage agglomeration methods (R hclust method
    # names are a subset of these).
    _LINKAGE_METHODS = frozenset(
        {"single", "complete", "average", "weighted", "centroid", "median", "ward"}
    )

    def validate(self) -> None:
        """Raise :class:`PipelineConfigError` for the first invalid value."""
        checks = {
            "fdr_threshold": 0 < self.fdr_threshold <= 1,
            "occurrence_threshold": self.occurrence_threshold >= 1,
            "cluster_height": self.cluster_height > 0,
            "cluster_min_size": self.cluster_min_size >= 1,
            "cluster_method": self.cluster_method in self._LINKAGE_METHODS,
            "mass_error_sd": self.mass_error_sd > 0,
            "rt_error_sd": self.rt_error_sd > 0,
            "n_jobs": self.n_jobs == -1 or self.n_jobs >= 1,
        }
        for field, valid in checks.items():
            if not valid:
                raise PipelineConfigError(
                    f"{field} has invalid value {getattr(self, field)!r}"
                )

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-compatible parameter mapping."""
        return asdict(self)
