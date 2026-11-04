# /// script
# requires-python = ">=3.11"
# dependencies = ["unidec==8.2.1"]
# ///
"""Headless IsoGen: predict the isotope distribution of a neutral mass via
the bundled C neural net / FFT engine (isogen.dll — no PyTorch needed).

Usage: isogen_predict.py <neutral_mass> [molecule type, default PEPTIDE]

stdout: TSV isotope distribution (isotope_index = 13C offset, relative
abundance, rows >= 1% of max) + a `#` footer with the run report.
"""

import contextlib
import sys

import _speak

# CLI-documentation contract: the skill-runtime offline assertion probes every
# dependency-bearing script with --help; the argv passthrough below would treat
# "--help" as an input mass, so declare it head-of-file.
if "--help" in sys.argv or "-h" in sys.argv:
    print(__doc__.strip())
    sys.exit(0)

with contextlib.redirect_stdout(sys.stderr):
    try:
        from unidec.IsoDec.isogenwrapper import nn_gen_isodist

        mass = float(sys.argv[1])
        mol_type = sys.argv[2] if len(sys.argv) > 2 else "PEPTIDE"
        dist = nn_gen_isodist(mass, type=mol_type, isolen=64)
    except Exception as exc:  # noqa: BLE001 - one error surface
        _speak.fail("ENGINE_FAILURE", f"{type(exc).__name__}: {exc}")

top = max(dist)
rows = [(i, a) for i, a in enumerate(dist) if a >= 0.01 * top]
for i, abundance in rows:
    print(f"{i}\t{abundance:.4f}")
_speak.note(f"type={mol_type}, peak isotope index {max(rows, key=lambda r: r[1])[0]} at mass {mass}")
_speak.finish(len(rows))
