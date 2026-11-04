# /// script
# requires-python = ">=3.11"
# dependencies = ["unidec==8.2.1"]
# ///
"""Headless findpeaks: centroid/peak-detect one spectrum via the UniDec
engine's peak detector.

Consumes whatever the 1D engine opens (two-column .txt, .mzML[.gz] with MS1
summed, .mgf via the shared MGF adapter). This is the profile→centroid step
that feeds IsoDec.

Usage: findpeaks.py <spectrum> [threshold, default 0.05, fraction of max]

stdout: TSV (mz, intensity) + a `#` footer with the run report.
"""

import contextlib
import os
import sys

import _speak

# CLI-documentation contract: the skill-runtime offline assertion probes every
# dependency-bearing script with --help; the engine arg passthrough below
# would treat "--help" as an input path, so declare it head-of-file.
if "--help" in sys.argv or "-h" in sys.argv:
    print(__doc__.strip())
    sys.exit(0)

with contextlib.redirect_stdout(sys.stderr):
    try:
        from unidec.engine import UniDec

        from _formats import load

        threshold = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
        input_path, _ = load(sys.argv[1])
        eng = UniDec()
        eng.open_file(input_path, os.path.dirname(input_path) or None, silent=True)
        peaks = eng.get_spectrum_peaks(threshold=threshold)
    except Exception as exc:  # noqa: BLE001 - one error surface
        _speak.fail("ENGINE_FAILURE", f"{type(exc).__name__}: {exc}")

for mz, intensity in peaks:
    print(f"{mz:.5f}\t{intensity:.2f}")
if len(peaks) == 0:
    _speak.warn("no peaks above threshold — raise or lower the threshold and retry")
_speak.note(f"threshold={threshold} of max intensity")
_speak.finish(len(peaks))
