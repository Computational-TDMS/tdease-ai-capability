"""Persistence for :class:`~toppic.pipeline.PipelineResult`."""

import importlib.util
import json
import shutil
import tempfile
from pathlib import Path
from typing import Literal

from .errors import ExportError
from .result import PipelineResult


def _has_pyarrow() -> bool:
    return importlib.util.find_spec("pyarrow") is not None


def export_result(
    result: PipelineResult,
    output_dir: str | Path,
    format: Literal["csv", "parquet"] = "csv",
    overwrite: bool = False,
) -> Path:
    """Write result tables and a JSON run manifest."""
    if format not in {"csv", "parquet"}:
        raise ExportError(f"Unsupported export format: {format}")
    if format == "parquet" and not _has_pyarrow():
        raise ExportError(
            "Parquet export requires the 'parquet' extra: "
            "install toppic-py[parquet]"
        )

    output = Path(output_dir).resolve()
    if output.exists() and any(output.iterdir()) and not overwrite:
        raise ExportError(f"Output directory is not empty: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent)
    )

    try:
        extension = "csv" if format == "csv" else "parquet"
        tables = {
            "expression": result.expression,
            "features": result.features,
            "samples": result.samples,
        }
        outputs: dict[str, str] = {}
        if not result.calibration.empty:
            tables["calibration"] = result.calibration
        for name, table in tables.items():
            filename = f"{name}.{extension}"
            destination = temporary / filename
            if format == "csv":
                table.to_csv(destination, index=True)
            else:
                table.to_parquet(destination, index=True)
            outputs[name] = filename

        manifest = dict(result.manifest)
        manifest.update(
            {
                "status": manifest.get("status", "completed"),
                "outputs": outputs,
                "warnings": result.warnings,
                "stage_metrics": result.stage_metrics,
            }
        )
        (temporary / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True, default=str) + "\n",
            encoding="utf-8",
        )

        if output.exists():
            if output.is_dir():
                shutil.rmtree(output)
            else:
                output.unlink()
        temporary.rename(output)
        return output
    except ExportError:
        raise
    except Exception as exc:
        raise ExportError(f"Could not export result to {output}: {exc}") from exc
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)
