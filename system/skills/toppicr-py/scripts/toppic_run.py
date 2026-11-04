# /// script
# requires-python = ">=3.11"
# dependencies = ["pandas", "numpy", "scipy", "scikit-misc>=0.5.3,<0.6", "scikit-learn", "statsmodels", "biopython", "joblib"]
# ///
"""TopPICR-Py quantitative replica (skill entry point).

Runs the full TopPICR-parity pipeline (PrSM + FASTA + optional MS1 feature
files) to a feature x sample expression matrix. This is a thin wrapper over the
vendored ``toppic`` package's ``cli.main`` — the analysis logic lives in
``toppic/`` next to this script, so the behaviour matches the standalone
``toppic run`` CLI exactly.

Usage (via the system-assets prefix, quoted):
  "$TDEASE_SYSTEM_ASSETS/skills/toppicr-py/scripts/toppic_run.py" \
    --project-dir <project> --output-dir <out> \
    --fdr-threshold 0.01 --occurrence-threshold 2
or with explicit files:
  .../toppic_run.py --fasta db.fasta.gz --ms2 a_prsm.tsv --ms2 b_prsm.tsv \
    --ms1 a.feature --output-dir results
"""

import contextlib
import sys
from pathlib import Path

# Make the vendored ``toppic`` package importable (it is the directory this
# script sits in). Must precede any ``from toppic...`` import.
_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

from toppic.cli import build_parser, run
from toppic_plot import main as report_main

if __name__ == "__main__":
    # The vendored CLI has a `run` subcommand; make it implicit so a caller can
    # pass flags directly (toppic_run.py --project-dir ... -> run --project-dir).
    argv = sys.argv[1:]
    if not argv or argv[0] not in ("run", "--version"):
        argv = ["run", *argv]
    parser = build_parser()
    args = parser.parse_args(argv)
    output_dir = str(args.output_dir)
    with contextlib.redirect_stdout(sys.stderr):
        code = run(args, parser)
    if code:
        raise SystemExit(code)
    raise SystemExit(report_main(["report", "--project", output_dir]))
