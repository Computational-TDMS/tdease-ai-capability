# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
#!/usr/bin/env python3
"""MassQL plotting: run SELECT statements against an mzML and draw the traces.

Each statement becomes one trace on a shared RT axis. Rows carrying
matchedProducts plot their summed matched intensity (XIC); rows without them
plot totalIonCurrent (TIC). Emits one self-contained chart-forge HTML
(TIC/XIC split into stacked magnitude domains) + data.json + manifest.json.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from _shared.chart_forge import emit_html
from _shared.kernel_io import result_of
from _shared.manifest import write_manifest
from _shared.native import binary as ms_binary, run_json


def trace_from_result(result: dict) -> tuple[list[float], list[float], str]:
    """Rows with matched peaks are an XIC (summed matched intensity); else TIC."""
    rows = sorted(result["spectra"], key=lambda row: row["retentionTime"])
    has_matched = any(row["matchedProducts"] for row in rows)
    kind = "xic" if has_matched else "tic"
    xs = [float(row["retentionTime"]) for row in rows]
    if has_matched:
        ys = [float(sum(p["intensity"] for p in row["matchedProducts"])) for row in rows]
    else:
        ys = [float(row["totalIonCurrent"]) for row in rows]
    return xs, ys, kind


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mzml", required=True, type=Path, help="mzML file (read-only)")
    parser.add_argument(
        "--statement",
        required=True,
        action="append",
        help="SELECT statement; repeat for multiple traces",
    )
    parser.add_argument(
        "--label",
        action="append",
        help="trace label per --statement (default: statement text)",
    )
    parser.add_argument("--out-dir", required=True, type=Path, help="output directory")
    parser.add_argument("--title", default=None, help="figure title (default: mzML stem)")
    parser.add_argument("--binary", default=None, help="tdeasenexus-toolbox path override")
    args = parser.parse_args()

    if args.label and len(args.label) != len(args.statement):
        sys.exit("massql: --label count must match --statement count")
    labels = args.label or args.statement

    binary = ms_binary(args.binary)
    mzml = args.mzml.resolve()
    if not mzml.is_file():
        sys.exit(f"massql: mzML not found: {mzml}")

    traces = []
    manifest_rows = []
    for label, statement in zip(labels, args.statement):
        envelope = run_json([binary, "select", "--statement", statement, str(mzml)])
        result = result_of(envelope, source="massql select")
        xs, ys, kind = trace_from_result(result)
        if result.get("truncated"):
            print(
                f"massql: warning — trace {label!r} truncated at LIMIT; "
                "rerun with a higher LIMIT (max 10000) or an rt BETWEEN filter",
                file=sys.stderr,
            )
        traces.append((label, kind, xs, ys))
        manifest_rows.append(
            {
                "label": label,
                "statement": statement,
                "kind": kind,
                "points": len(xs),
                "matchedSpectrumCount": result["matchedSpectrumCount"],
                "truncated": result["truncated"],
            }
        )

    stem = args.out_dir / f"massql-{mzml.stem}"
    args.out_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "data": {
            "traces": [{"kind": kind, "label": f"{label} ({kind.upper()})", "rt": xs, "y": ys}
                       for label, kind, xs, ys in traces],
        },
        "meta": {
            "title": args.title or mzml.stem,
            "source": mzml.name,
            "transform": "MassQL SELECT rows -> TIC (totalIonCurrent) or XIC (summed matched intensity); stacked magnitude domains per kind",
            "n": str(sum(len(xs) for _, _, xs, _ in traces)),
        },
    }
    out_html, contract_path = emit_html("massql-traces", payload, stem)
    write_manifest(
        out_html, figure="massql-traces", script=Path(__file__).name,
        params={"statements": args.statement},
        inputs=[str(mzml), str(contract_path)],
        caption={"traces": manifest_rows},
        plot_data=manifest_rows,
    )

    print(json.dumps({"figure": str(out_html), "traces": manifest_rows}))


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, ValueError) as exc:
        sys.exit(f"massql: {exc}")
