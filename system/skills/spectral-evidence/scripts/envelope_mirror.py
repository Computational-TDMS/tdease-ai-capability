"""Mirror-view panels for fragment isotope windows — reached through
``spectral.py show``, never directly.

One panel per (model, ion): the observed window up, the theoretical comb down,
in-slot peaks ppm-labelled, so several models' predictions of the SAME window
land in one figure and a populated window in one model next to an empty one in
the other is visible rather than asserted.

Two functions: :func:`_panels` turns one ``envelope-regions`` artifact into
chart-ready panels, :func:`_quick_extractions` runs the kernel per entry and
leaves each extraction JSON beside the figure — the same auditable artifact the
artifact path would start from.

Correlation verdicts come from the `fragment` run that produced the card, never
from here: a window can look well populated and still be rejected on the
full-slot Pearson, so the caller states the kernel's verdict in the subtitle
alongside this function's own slot count.
"""

from __future__ import annotations

import json
from pathlib import Path

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from _shared.kernel_io import result_of
from _shared.native import binary, run_json

# Comparison pair order: the first source is the model under test, the second
# its challenger — teal leads, displaced/challenge reads warm.
MIRROR_PALETTE = ("#0F766E", "#C2410C", "#1D4ED8", "#7C3AED", "#B45309", "#0E7490")


class FigureError(Exception):
    pass


def _panels(doc: dict, path: str, ppm_tol: float) -> list[dict]:
    """Chart-ready panels from one envelope-regions artifact."""
    regions = result_of(doc, source=path).get("regions") or []
    panels = []
    for r in regions:
        theo = r["theoretical"]
        theo_hi = max((t["relativeIntensity"] for t in theo), default=0.0) or 1.0
        obs = r.get("observed") or []
        obs_hi = max((o["intensity"] for o in obs), default=0.0)
        # errPpm missing means the schema drifted, not that the peak missed
        # its slot: count both separately so "0/N slots" stays interpretable.
        missing_errppm = sum(1 for o in obs if o.get("errPpm") is None)
        in_slot = [o for o in obs
                   if o.get("errPpm") is not None
                   and abs(float(o["errPpm"])) <= ppm_tol]
        # in-slot apex only: the raw max of the window is routinely a
        # co-eluting off-slot peak far taller than the envelope itself
        apex = max(in_slot, key=lambda o: o["intensity"])["nearestNeutronOffset"] if in_slot else None
        name = f"{r['ion']}{r['number']} z={r['charge']}"
        sub = (f"{name} · no observed peaks in window" if not obs else
               f"{name} · {len({o['nearestNeutronOffset'] for o in in_slot})}/{len(theo)}"
               f" slots ±{ppm_tol:g} ppm"
               + (f" · in-slot apex M+{apex}" if apex is not None else "")
               + (f" · {missing_errppm} peaks lack errPpm" if missing_errppm else ""))
        scan = (doc.get("params") or {}).get("scan")
        stem = Path(path).stem
        panels.append({
            "title": f"{stem} · scan {scan}" if scan is not None else stem,
            "subtitle": sub,
            "observed": {
                "x": [o["mz"] for o in obs],
                "y": [100.0 * o["intensity"] / obs_hi if obs_hi else 0.0 for o in obs],
                "ppm": [o.get("errPpm") for o in obs],
                "matched": [o in in_slot for o in obs],
            },
            "theoretical": {
                "x": [t["mz"] for t in theo],
                "y": [100.0 * t["relativeIntensity"] / theo_hi for t in theo],
                "offsets": [t["neutronOffset"] for t in theo],
            },
            "_facts": {"ion": r["ion"], "number": r["number"], "charge": r["charge"],
                       "slots": len({o["nearestNeutronOffset"] for o in in_slot}),
                       "slotTotal": len(theo), "apex": apex},
        })
    return panels


def _quick_extractions(spec: list, out_plot: str, ppm_tol: float,
                       ms_bin: str | None = None) -> tuple[list[dict], list[str]]:
    """Run the kernel per entry; return (panels, extraction paths).

    Each extraction JSON lands beside the figure as ``<stem>.regions<i>.json``,
    so the figure always leaves the same auditable artifact behind it.
    ``ms_bin`` names the toolbox binary; omitted, it is resolved from the
    environment and the packaged copy.

    Styling (label/sub/color) needs exactly one region per entry — a bare ion
    expands to every charge state and the panel would stop being one window.
    """
    stem = Path(out_plot).with_suffix("")
    panels: list[dict] = []
    written: list[str] = []
    for i, entry in enumerate(spec):
        missing = [k for k in ("mzml", "scan", "proforma", "ions") if not entry.get(k)]
        if missing:
            raise FigureError(f"EXTRACT_ENTRY_{i}: missing {', '.join(missing)}")
        cmd = [ms_bin or binary(), "envelope-regions", str(entry["mzml"]),
               "--scan", str(entry["scan"]),
               "--proforma", str(entry["proforma"]),
               "--ions", str(entry["ions"])]
        if entry.get("pad") is not None:
            cmd += ["--pad", str(entry["pad"])]
        doc = run_json(cmd, label=f"envelope-regions #{i}")
        path = f"{stem}.regions{i}.json"
        Path(path).write_text(json.dumps(doc, ensure_ascii=False), encoding="utf-8")
        written.append(path)
        got = _panels(doc, path, ppm_tol)
        if not got:
            raise FigureError(f"EXTRACT_ENTRY_{i}: no region for --ions {entry['ions']} "
                              "(no theoretical fragment at this charge, or an empty "
                              "isotope envelope)")
        styled = [k for k in ("label", "sub", "color") if entry.get(k)]
        if styled and len(got) > 1:
            raise FigureError(
                f"EXTRACT_ENTRY_{i}: label/sub/color need exactly one region, "
                f"--ions {entry['ions']} produced {len(got)}; give each ion its "
                "charge to keep one window per panel")
        for panel in got:
            if "label" in styled:
                panel["title"] = str(entry["label"])
            if "sub" in styled:
                panel["subtitle"] = str(entry["sub"])
            if "color" in styled:
                panel["color"] = str(entry["color"])
        panels.extend(got)
    return panels, written
