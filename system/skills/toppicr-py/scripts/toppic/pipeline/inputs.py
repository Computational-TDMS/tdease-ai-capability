"""Canonical local-file input contract and directory discovery."""

from dataclasses import dataclass
from pathlib import Path

from .errors import InputValidationError


@dataclass(frozen=True)
class PipelineInputs:
    """Resolved FASTA, MS2 PrSM, and optional MS1 feature files."""

    fasta_files: tuple[Path, ...]
    ms2_files: tuple[Path, ...]
    ms1_files: tuple[Path, ...] = ()

    def __post_init__(self) -> None:
        for name in ("fasta_files", "ms2_files", "ms1_files"):
            object.__setattr__(
                self,
                name,
                tuple(Path(path).resolve() for path in getattr(self, name)),
            )

    @classmethod
    def from_project_directory(
        cls,
        path: str | Path,
        use_fdr: bool = False,
    ) -> "PipelineInputs":
        """Discover conventional TopPIC project subdirectories."""
        root = Path(path).resolve()
        ms2_dir = root / ("toppic_fdr" if use_fdr else "toppic")
        inputs = cls(
            fasta_files=tuple(sorted((root / "fasta").glob("*.fasta*"))),
            ms2_files=tuple(sorted(ms2_dir.glob("*ms2_toppic_prsm.tsv"))),
            ms1_files=tuple(sorted((root / "topfd").glob("*ms1.feature"))),
        )
        inputs.validate()
        return inputs

    def validate(self) -> None:
        """Validate required file groups and file readability."""
        if not self.fasta_files:
            raise InputValidationError("At least one FASTA file is required")
        if len(self.fasta_files) != 1:
            raise InputValidationError(
                "PipelineInputs requires exactly one FASTA database"
            )
        if not self.ms2_files:
            raise InputValidationError("At least one MS2 PrSM file is required")
        for name, paths in (
            ("MS2", self.ms2_files),
            ("MS1", self.ms1_files),
        ):
            if len({path.parent for path in paths}) > 1:
                raise InputValidationError(
                    f"All {name} files must be staged in the same directory"
                )
        for path in self.fasta_files + self.ms2_files + self.ms1_files:
            if not path.is_file():
                raise InputValidationError(f"Input file does not exist: {path}")
            if not path.stat().st_size:
                raise InputValidationError(f"Input file is empty: {path}")
