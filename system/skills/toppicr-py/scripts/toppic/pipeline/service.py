"""Shared end-to-end TopPIC scientific workflow."""

from collections.abc import Callable
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal, TypeVar

import pandas as pd

import toppic
from toppic.alignment.mass_calibrate import (
    calc_error,
    mass_reference_mask,
    recalibrate_mass,
)
from toppic.alignment.rt_align import align_rt, find_ref_ds, form_model
from toppic.annotation.filter import filter_by_count, rm_false_gene
from toppic.annotation.map_proteoform import map_proteoform
from toppic.annotation.proteoform_level import set_pf_level
from toppic.clustering.cluster import cluster
from toppic.clustering.merge import create_pcg
from toppic.fdr.control import apply_evalue_cutoff, find_evalue_cutoff
from toppic.features.match import match_features
from toppic.features.metadata import create_mdata
from toppic.inference.protein import infer_proteins
from toppic.io.msnset import to_msnset
from toppic.io.read_toppic import read_toppic

from .config import PipelineConfig
from .errors import PipelineStageError, TopPICError
from .events import ProgressEvent
from .export import export_result
from .inputs import PipelineInputs
from .result import PipelineResult

T = TypeVar("T")
ProgressCallback = Callable[[ProgressEvent], None]


class TopPICPipeline:
    """Coordinate the validated algorithms behind every public interface."""

    def __init__(self, config: PipelineConfig | None = None):
        self.config = config or PipelineConfig()

    def run(
        self,
        inputs: PipelineInputs,
        output_dir: str | Path | None = None,
        *,
        output_format: Literal["csv", "parquet"] = "csv",
        overwrite: bool = False,
        on_progress: ProgressCallback | None = None,
    ) -> PipelineResult:
        """Execute the full workflow and optionally persist its result."""
        self.config.validate()
        inputs.validate()
        started_at = datetime.now(timezone.utc)
        stage_metrics: dict[str, dict[str, int | float]] = {}
        warnings: list[str] = []
        total = 21 if inputs.ms1_files else 17
        index = 0

        def stage(
            name: str,
            message: str,
            operation: Callable[[], T],
        ) -> T:
            nonlocal index
            index += 1
            if on_progress is not None:
                on_progress(
                    ProgressEvent(
                        stage=name,
                        index=index,
                        total=total,
                        message=message,
                    )
                )
            try:
                value = operation()
            except TopPICError:
                raise
            except Exception as exc:
                raise PipelineStageError(name, exc) from exc
            shape = getattr(value, "shape", None)
            if shape is not None and len(shape):
                stage_metrics[name] = {"rows": int(shape[0])}
            return value

        ms2_parent, ms2_names = _parent_and_names(inputs.ms2_files)
        ms2 = stage(
            "read_ms2",
            "Reading TopPIC PrSM files",
            lambda: read_toppic(ms2_parent, ms2_names),
        )
        ms2 = stage(
            "gene_resolution",
            "Resolving ambiguous genes",
            lambda: rm_false_gene(ms2),
        )
        fasta = inputs.fasta_files[0]
        ms2 = stage(
            "fasta_mapping",
            "Mapping proteoforms to FASTA",
            lambda: map_proteoform(ms2, str(fasta.parent), fasta.name),
        )
        ms2 = stage(
            "occurrence_filter",
            "Filtering low-occurrence sequences",
            lambda: filter_by_count(
                ms2,
                ["dataset", "scans"],
                "clean_seq",
                self.config.occurrence_threshold,
            ),
        )
        cutoffs = stage(
            "fdr_cutoffs",
            "Finding E-value cutoffs",
            lambda: find_evalue_cutoff(ms2, self.config.fdr_threshold),
        )
        ms2 = stage(
            "fdr_filter",
            "Applying FDR control",
            lambda: apply_evalue_cutoff(ms2, cutoffs),
        )
        ms2 = stage(
            "protein_inference",
            "Inferring proteins",
            lambda: infer_proteins(ms2),
        )
        ms2 = stage(
            "proteoform_level",
            "Classifying proteoform levels",
            lambda: set_pf_level(ms2),
        )
        reference = stage(
            "reference_selection",
            "Selecting retention-time reference",
            lambda: find_ref_ds(ms2),
        )
        model = stage(
            "rt_model",
            "Building retention-time alignment model",
            lambda: form_model(ms2, reference),
        )
        ms2 = stage(
            "align_ms2",
            "Aligning MS2 retention times",
            lambda: align_rt(ms2, model, "feature_apex"),
        )
        error = stage(
            "error_model",
            "Estimating mass and retention-time error",
            lambda: calc_error(ms2, reference),
        )
        ms2 = stage(
            "recalibrate_ms2",
            "Recalibrating MS2 precursor masses",
            lambda: recalibrate_mass(ms2, error, "precursor_mass"),
        )
        calibration = ms2[[
            "dataset", "cv", "scans", "proteoform", "feature_apex", "rt_align",
            "precursor_mass", "adjusted_precursor_mass", "recal_mass",
        ]].copy()
        calibration["mass_reference"] = mass_reference_mask(ms2)
        ms2 = stage(
            "clustering",
            "Clustering proteoform observations",
            lambda: cluster(
                ms2,
                error,
                method=self.config.cluster_method,
                height=self.config.cluster_height,
                min_size=self.config.cluster_min_size,
            ),
        )
        ms2 = stage(
            "pcgroup_merge",
            "Merging isotopologue clusters",
            lambda: create_pcg(ms2, error),
        )

        if inputs.ms1_files:
            ms1_parent, ms1_names = _parent_and_names(inputs.ms1_files)
            ms1 = stage(
                "read_ms1",
                "Reading TopFD MS1 features",
                lambda: read_toppic(ms1_parent, ms1_names),
            )
            ms1 = stage(
                "align_ms1",
                "Aligning MS1 retention times",
                lambda: align_rt(ms1, model, "time_apex"),
            )
            ms1 = stage(
                "recalibrate_ms1",
                "Recalibrating MS1 masses",
                lambda: recalibrate_mass(ms1, error, "mass"),
            )
            matched = stage(
                "feature_matching",
                "Matching MS1 and MS2 features",
                lambda: match_features(
                    ms2,
                    ms1,
                    error,
                    n_mme_sd=self.config.mass_error_sd,
                    n_rt_sd=self.config.rt_error_sd,
                    n_jobs=self.config.n_jobs,
                ),
            )
        else:
            warnings.append(
                "No MS1 feature files were supplied; expression values use "
                "aggregated MS2 feature intensity."
            )
            matched = stage(
                "ms2_intensity_fallback",
                "Aggregating MS2 feature intensities",
                lambda: _aggregate_ms2_intensity(ms2),
            )

        metadata = stage(
            "feature_metadata",
            "Building feature metadata",
            lambda: create_mdata(ms2, error),
        )
        raw_result = stage(
            "expression_matrix",
            "Building expression matrix",
            lambda: to_msnset(matched, meta=metadata),
        )
        completed_at = datetime.now(timezone.utc)
        manifest = {
            "status": "completed",
            "version": toppic.__version__,
            "python_started_at": started_at.isoformat(),
            "python_completed_at": completed_at.isoformat(),
            "inputs": _input_manifest(inputs),
            "config": self.config.to_dict(),
            "calibration": {"reference_dataset": reference,
                            "stage": "recalibrate_ms2"},
            "matrix": {
                "features": int(raw_result.expression.shape[0]),
                "samples": int(raw_result.expression.shape[1]),
            },
        }
        result = PipelineResult(
            expression=raw_result.expression,
            features=raw_result.features,
            samples=raw_result.samples,
            stage_metrics=stage_metrics,
            warnings=warnings,
            manifest=manifest,
            calibration=calibration,
        )
        if output_dir is not None:
            export_result(
                result,
                output_dir,
                format=output_format,
                overwrite=overwrite,
            )
        return result


def _parent_and_names(paths: tuple[Path, ...]) -> tuple[Path, list[str]]:
    parents = {path.parent for path in paths}
    if len(parents) != 1:
        raise ValueError("Each input file group must share one parent directory")
    return next(iter(parents)), [path.name for path in paths]


def _aggregate_ms2_intensity(ms2: pd.DataFrame) -> pd.DataFrame:
    keys = ["dataset"]
    if "cv" in ms2.columns:
        keys.append("cv")
    keys.extend(["gene", "pc_group"])
    return (
        ms2.groupby(keys, dropna=False)["feature_intensity"]
        .max()
        .reset_index()
        .rename(columns={"feature_intensity": "intensity"})
    )


def _input_manifest(inputs: PipelineInputs) -> dict[str, list[dict[str, object]]]:
    def records(paths: tuple[Path, ...]) -> list[dict[str, object]]:
        return [
            {"path": str(path), "size_bytes": int(path.stat().st_size)}
            for path in paths
        ]

    return {
        "fasta": records(inputs.fasta_files),
        "ms2": records(inputs.ms2_files),
        "ms1": records(inputs.ms1_files),
    }
