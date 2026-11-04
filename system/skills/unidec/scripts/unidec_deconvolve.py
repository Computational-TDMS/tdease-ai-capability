# /// script
# requires-python = ">=3.11"
# dependencies = ["unidec==8.2.1"]
# ///
"""Headless UniDec: deconvolve one MS1 spectrum (two-column m/z-intensity
.txt, .mzML[.gz] with all MS1 summed, or .mgf with all blocks summed —
the system's peak-list standard).

stdout: TSV peak table (mass, height, area, dscore) + a `#` footer with the
run report. Engine artifacts land in <stem>_unidecfiles/ next to the input;
delete that folder for a pristine re-run. Tuning = set fields on eng.config
before process_data() (masslb/massub, startz/endz, massbins, adductmass,
mzsig, peakthresh...).
"""

import contextlib
import os
import sys

import _speak

# CLI-documentation contract: the skill-runtime offline assertion probes every
# dependency-bearing script with --help; the argv passthrough below would treat
# "--help" as an input path, so declare it head-of-file.
if "--help" in sys.argv or "-h" in sys.argv:
    print(__doc__.strip())
    sys.exit(0)

with contextlib.redirect_stdout(sys.stderr):
    try:
        from unidec.engine import UniDec

        from _formats import load

        input_path, _ = load(sys.argv[1])
        eng = UniDec()
        eng.open_file(input_path, os.path.dirname(input_path) or None, silent=True)
        eng.process_data(silent=True)
        eng.run_unidec(silent=True)
        eng.pick_peaks()
        pks = sorted(eng.pks.peaks, key=lambda p: -p.height)
    except Exception as exc:  # noqa: BLE001 - one error surface
        _speak.fail("ENGINE_FAILURE", f"{type(exc).__name__}: {exc}")

for pk in pks:
    print(f"{pk.mass:.2f}\t{pk.height:.1f}\t{pk.area:.1f}\t{pk.dscore:.2f}")
if not pks:
    _speak.warn("no peaks — empty input, or mass/charge window excluded everything")
_speak.note(f"mass window {eng.config.masslb}-{eng.config.massub} Da, z {eng.config.startz}-{eng.config.endz}")
_speak.note(f"artifacts in {eng.config.udir}")
_speak.finish(len(pks))
