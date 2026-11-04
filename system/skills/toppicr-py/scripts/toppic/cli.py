"""Command-line adapter for the shared TopPIC pipeline service."""

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from toppic import __version__
from toppic.pipeline import (
    PipelineConfig,
    PipelineInputs,
    ProgressEvent,
    TopPICError,
    TopPICPipeline,
)


def build_parser() -> argparse.ArgumentParser:
    """Build the public command-line parser."""
    parser = argparse.ArgumentParser(
        prog="toppic",
        description="TopPICR-compatible quantitative analysis pipeline",
    )
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    run = subparsers.add_parser("run", help="run the complete pipeline",
                               formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    run.add_argument("--project-dir", type=Path, help="project with fasta/, toppic/, optional topfd/")
    run.add_argument("--use-fdr", action="store_true", help="discover toppic_fdr/ instead of toppic/; does not disable pipeline FDR")
    run.add_argument("--fasta", type=Path, action="append", default=[])
    run.add_argument("--ms2", type=Path, action="append", default=[])
    run.add_argument("--ms1", type=Path, action="append", default=[])
    run.add_argument("--output-dir", type=Path, required=True)
    run.add_argument("--format", choices=("csv", "parquet"), default="csv")
    run.add_argument("--overwrite", action="store_true")
    run.add_argument("--fdr-threshold", type=float, default=0.01, help="target gene-level decoy fraction for E-value cutoff selection")
    run.add_argument("--occurrence-threshold", type=int, default=2, help="minimum distinct dataset/scan observations per clean sequence")
    run.add_argument("--cluster-height", type=float, default=3.0)
    run.add_argument("--cluster-min-size", type=int, default=2)
    run.add_argument(
        "--cluster-method",
        default="ward",
        choices=["single", "complete", "average", "weighted", "centroid", "median", "ward"],
        help="hierarchical linkage method (paper preset: single)",
    )
    run.add_argument(
        "--paper-preset",
        action="store_true",
        help="use Martin et al. settings: single linkage, height 4, min cluster size 2",
    )
    run.add_argument("--mass-error-sd", type=float, default=3.0, help="MS1 matching window multiplier of estimated mass-error width (not ppm)")
    run.add_argument("--rt-error-sd", type=float, default=3.0, help="MS1 matching window multiplier of estimated RT-error width (not minutes)")
    run.add_argument("--n-jobs", type=int, default=1)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CLI and return a process-compatible exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(args, parser)


def run(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    """Execute already-parsed CLI options through the shared pipeline."""
    if args.command != "run":
        parser.error(f"Unsupported command: {args.command}")

    explicit = bool(args.fasta or args.ms2 or args.ms1)
    if args.project_dir is not None and explicit:
        parser.error("--project-dir cannot be combined with explicit input files")
    if args.project_dir is None and not (args.fasta and args.ms2):
        parser.error(
            "provide --project-dir or both --fasta and --ms2 input files"
        )

    cluster_method = "single" if args.paper_preset else args.cluster_method
    cluster_height = 4.0 if args.paper_preset else args.cluster_height
    cluster_min_size = 2 if args.paper_preset else args.cluster_min_size
    config = PipelineConfig(
        fdr_threshold=args.fdr_threshold,
        occurrence_threshold=args.occurrence_threshold,
        cluster_height=cluster_height,
        cluster_min_size=cluster_min_size,
        cluster_method=cluster_method,
        mass_error_sd=args.mass_error_sd,
        rt_error_sd=args.rt_error_sd,
        n_jobs=args.n_jobs,
    )

    try:
        if args.project_dir is not None:
            inputs = PipelineInputs.from_project_directory(
                args.project_dir,
                use_fdr=args.use_fdr,
            )
        else:
            inputs = PipelineInputs(
                fasta_files=tuple(args.fasta),
                ms2_files=tuple(args.ms2),
                ms1_files=tuple(args.ms1),
            )

        def show_progress(event: ProgressEvent) -> None:
            print(
                f"[{event.index}/{event.total}] {event.message}",
                file=sys.stderr,
            )

        result = TopPICPipeline(config).run(
            inputs,
            output_dir=args.output_dir,
            output_format=args.format,
            overwrite=args.overwrite,
            on_progress=show_progress,
        )
    except TopPICError as exc:
        print(f"toppic: {exc}", file=sys.stderr)
        return 2

    rows, columns = result.expression.shape
    print(
        f"Completed: {rows} features x {columns} samples; "
        f"{len(result.warnings)} warning(s); output={args.output_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
