"""Shared input handling for the unidec skill scripts.

The engine reads two-column m/z-intensity text and mzML natively. MGF — the
system's peak-list standard (tdeasenexus-toolbox `spectrum --out-format mgf`,
isoforma feature MGFs) — is parsed here: all TITLE blocks are summed into one
spectrum, matching the accumulate convention.

`load(path)` returns an (input_path, Nx2-array-or-None) pair: the array is
None when the engine can read the file itself, and a summed spectrum when the
input was MGF (written to a temp .txt that the engine reads natively).
"""

import os
import tempfile

import numpy as np

import _speak


def _mgf_sum(path, bin_ppm=20.0):
    blocks = []
    current = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or not line[0].isdigit():
                if current:
                    blocks.append(current)
                    current = []
                continue
            current.append([float(x) for x in line.split()[:2]])
    if current:
        blocks.append(current)
    if not blocks:
        raise ValueError(f"no peak blocks found in MGF: {path}")
    data = np.vstack([np.array(b, dtype=float) for b in blocks])
    raw = len(data)
    binned = _bin_peaks(data, bin_ppm)
    _speak.note(f"mgf: {len(blocks)} blocks, {raw} peaks -> {len(binned)} after {bin_ppm:g} ppm bin")
    return binned


def _bin_peaks(data, bin_ppm):
    """Merge peaks within bin_ppm of each other (intensity-weighted m/z),
    the accumulate contract shared with quantify.py's _sum_peak_lists."""
    data = data[np.argsort(data[:, 0])]
    tol = data[:, 0] * bin_ppm * 1e-6
    mz, inten, group = [], [], []
    for row, t in zip(data, tol):
        if group and row[0] - group[-1][0] > max(t, group[0][2]):
            mz.append(sum(p[0] * p[1] for p in group) / sum(p[1] for p in group))
            inten.append(sum(p[1] for p in group))
            group = []
        group.append((row[0], row[1], t))
    if group:
        mz.append(sum(p[0] * p[1] for p in group) / sum(p[1] for p in group))
        inten.append(sum(p[1] for p in group))
    return np.column_stack([mz, inten])


def load(path):
    path = os.path.abspath(path)
    ext = os.path.splitext(path)[1].lower()
    if ext == ".mgf":
        data = _mgf_sum(path)
        txt = os.path.join(tempfile.gettempdir(), os.path.basename(path) + ".unidec.txt")
        np.savetxt(txt, data)
        return txt, data
    if ext in (".mzml", ".gz", ".mzxml"):
        return path, None
    return path, None  # plain two-column text: engine-native
