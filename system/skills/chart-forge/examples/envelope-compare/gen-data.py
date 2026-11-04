# Regenerates data.json for the envelope-compare example through the REAL
# spectral adapter (the `spectral.py show --ms1` composition overlay) on the
# committed fictional-shift-example artifacts: the observed F12 feature envelope against
# the model A and model B complete-model theories, with the
# kernel's fit metrics passing through to the stats block. The example is
# adapter output, never hand-written.
#
#   python gen-data.py
#
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[5]  # envelope-compare -> examples -> chart-forge -> skills -> system -> tdease-ai-capability -> workspace root
GEN = HERE.parents[2] / "spectral-evidence" / "scripts" / "spectral.py"
MS1 = (WORKSPACE / "data" / "workspace" / "runs" / "sample-a"
       / "analysis" / "fictional-shift-example" / "ms1")

cmd = [
    sys.executable, str(GEN), "show", "--ms1",
    "--observed", str(MS1 / "observed_F12.json"),
    "--theory", str(MS1 / "theory_oxibr.json"),
    "--fit", str(MS1 / "fit_oxibr.json"),
    "--theory", str(MS1 / "theory_meposho.json"),
    "--fit", str(MS1 / "fit_meposho.json"),
    "--theory-label", "dime+oxi+Br",
    "--theory-label", "model B",
    "--title", "F12 precursor envelope: model A vs model B candidates",
    "--out", str(HERE / "envelope.html"),
]
subprocess.run(cmd, check=True)
(HERE / "envelope.html").unlink()
(HERE / "envelope.manifest.json").unlink(missing_ok=True)
(HERE / "envelope.data.json").replace(HERE / "data.json")
print("data.json regenerated")
