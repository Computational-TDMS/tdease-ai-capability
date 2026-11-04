"""Public exception types for pipeline integrations."""


class TopPICError(Exception):
    """Base class for public TopPICPy failures."""


class InputValidationError(TopPICError):
    """Raised when input files do not satisfy the pipeline contract."""


class UnsupportedTopPICVersionError(InputValidationError):
    """Raised when a TopPIC or TopFD input format cannot be recognized."""


class PipelineConfigError(TopPICError):
    """Raised when a pipeline configuration value is invalid."""


class PipelineStageError(TopPICError):
    """Raised when a named scientific processing stage fails."""

    def __init__(self, stage: str, cause: Exception):
        self.stage = stage
        self.cause = cause
        super().__init__(f"Pipeline stage '{stage}' failed: {cause}")


class ExportError(TopPICError):
    """Raised when a completed result cannot be persisted."""
