# /// script
# requires-python = ">=3.11"
# dependencies = ["unidec==8.2.1"]
# ///
"""Headless IsoDec: isotopic deconvolution of centroided MS data via the
bundled C engine (isodeclib.dll — no PyTorch needed).

Consumes centroided input: two-column m/z-intensity .txt or .mgf (all blocks
summed and 20 ppm-binned). Predicts charge states with the embedded model
and deconvolves each matched isotope cluster.

stdout: TSV (mz, charge, monoisotopic mass, peak mass, average mass) + a `#`
footer with the run report. Export kernel MGFs with --min-intensity to keep
the noise floor out.
"""

import contextlib
import sys
from pathlib import Path

import _speak

# CLI-documentation contract: the skill-runtime offline assertion probes every
# dependency-bearing script with --help; the argv passthrough below would treat
# "--help" as an input path, so declare it head-of-file.
if "--help" in sys.argv or "-h" in sys.argv:
    print(__doc__.strip())
    sys.exit(0)

with contextlib.redirect_stdout(sys.stderr):
    try:
        import numpy as np

        from _formats import _mgf_sum

        from unidec.IsoDec.c_interface import IsoDecWrapper

        path = Path(sys.argv[1]).resolve()
        cents = _mgf_sum(str(path)) if path.suffix.lower() == ".mgf" else np.loadtxt(path)
        pks = sorted(IsoDecWrapper().process_spectrum(cents).peaks, key=lambda p: -p.peakint)
    except Exception as exc:  # noqa: BLE001 - one error surface
        _speak.fail("ENGINE_FAILURE", f"{type(exc).__name__}: {exc}")

for pk in pks:
    print(f"{pk.mz:.4f}\t{pk.z}\t{pk.monoiso:.3f}\t{pk.peakmass:.3f}\t{pk.avgmass:.3f}")

if not pks:
    _speak.warn("no clusters matched — input may be profile (IsoDec needs centroided data)")
else:
    charges = {}
    for pk in pks:
        charges[pk.z] = charges.get(pk.z, 0) + 1
    _speak.note("charge: " + " ".join(f"z{k}={v}" for k, v in sorted(charges.items())))
    if len(charges) == 1 and 1 in charges and len(pks) > 10:
        _speak.warn("all clusters z=1 — noise-dominated input? filter low-intensity peaks and retry")
_speak.finish(len(pks))
