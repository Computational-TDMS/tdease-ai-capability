"""Progress events emitted by the shared pipeline service."""

from dataclasses import dataclass


@dataclass(frozen=True)
class ProgressEvent:
    """A framework-neutral pipeline progress notification."""

    stage: str
    index: int
    total: int
    message: str
