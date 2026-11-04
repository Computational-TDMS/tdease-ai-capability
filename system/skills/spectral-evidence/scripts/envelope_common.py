"""
Shared machinery for the MS1 envelope figure modes: refine/theory loaders,
shape statistics, palette, and the chart-forge render adapter.
"""

import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from _shared.chart_forge import emit_html
from _shared.kernel_io import load_anchored, load_doc, result_of
from _shared.native import mass_constants


class FigureError(Exception):
    pass


@dataclass
class Envelope:
    mz: list[float]
    intensity: list[float]
    label: str


# Wong sequence, contrast orange leading: candidate 0 (best fit) leads in the
# Stable candidate colors shared by all envelope comparisons.
CANDIDATE_PALETTE = ("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")


def _stat_text(st: dict) -> str:
    """One-line fit summary for the figure's stats block (and tooltips)."""
    parts = []
    if st.get("pearson") is not None:
        parts.append(f"r={st['pearson']:.3f}")
    if st.get("delta_apex_mTh") is not None:
        parts.append(f"apex {st['delta_apex_mTh']:+.1f} mTh")
    return " · ".join(parts)


def _ranked_evaluations(result: dict) -> list[dict]:
    """Refine-shift evaluations, best-pearson first (unknown last)."""
    return sorted(result.get("evaluations") or [],
                  key=lambda e: (e.get("pearson") is None, -(e.get("pearson") or 0)))


def _evaluation_candidates(result: dict, theo: Envelope, base_mass: float | None,
                           charge: int | None) -> list[Envelope]:
    """One Envelope per refine-shift evaluation, registered at its own mass.

    Every isotope-offset candidate shares the theoretical isotope comb; a
    candidate's comb sits at ``(refinedMass - base refinedMass) / charge`` m/z
    from the base registration, so shifting reproduces each candidate exactly
    without re-running chemistry."""
    if charge in (None, 0) or base_mass is None:
        return []
    out = []
    for ev in _ranked_evaluations(result):
        mass = ev.get("refinedMass")
        shift_mz = (mass - base_mass) / charge if mass is not None else 0.0
        label = f"offset {ev.get('offset', 0):+d}" if ev.get("offset") is not None else "candidate"
        out.append(Envelope([mz + shift_mz for mz in theo.mz], list(theo.intensity), label))
    return out


def _load_refine(path: str) -> dict:
    doc = load_doc(path)
    result = result_of(doc, source=path)
    if "theoretical" not in result or "observed" not in result:
        raise FigureError(f"{path}: refine-shift JSON missing 'theoretical'/'observed' envelopes")
    result["_params"] = doc.get("params") or {}
    return result


def _load_theory_distribution(path: str) -> list[float]:
    anchored = load_anchored(path)
    peaks = anchored.get("peaks") or []
    center_masses = anchored.get("centerMasses") or []
    anchor_mass = anchored.get("anchorMass")
    if anchor_mass is None or not peaks or len(peaks) != len(center_masses):
        raise FigureError(f"{path}: theory JSON has incomplete anchored envelope")
    indexed: dict[int, float] = {}
    spacing = mass_constants()["carbon13Spacing"]
    for peak, center_mass in zip(peaks, center_masses):
        isotope_index = round((float(center_mass) - float(anchor_mass))
                              / spacing)
        # Peaks below the nominal anchor (M-1, M-2...) are part of the
        # distribution's shape, not noise to truncate: re-anchor the whole
        # distribution at its lowest isotope index instead of dropping them.
        indexed[isotope_index] = indexed.get(isotope_index, 0.0) + float(peak["intensity"])
    if not indexed:
        raise FigureError(f"{path}: theory JSON has no isotope peaks")
    values = [indexed.get(index, 0.0)
              for index in range(min(indexed), max(indexed) + 1)]
    total = sum(values)
    if total <= 0:
        raise FigureError(f"{path}: theory envelope intensity is zero")
    return [value / total for value in values]


def _compare_distributions(left: list[float], right: list[float], label: str) -> dict:
    size = max(len(left), len(right))
    left = left + [0.0] * (size - len(left))
    right = right + [0.0] * (size - len(right))
    return {
        "label": label,
        "left": left,
        "right": right,
        "tv_distance": 0.5 * sum(abs(a - b) for a, b in zip(left, right)),
        "left_apex_index": max(range(size), key=lambda index: left[index]),
        "right_apex_index": max(range(size), key=lambda index: right[index]),
    }


def _load_theory_profile(path: str) -> tuple[float, list[float]]:
    """(anchor neutral mass, true per-peak center masses, peak intensities
    normalized to max 1) from a `tdeasenexus-toolbox theory` JSON's anchored envelope."""
    anchored = load_anchored(path)
    peaks = anchored.get("peaks") or []
    anchor = anchored.get("anchorMass")
    if anchor is None or not peaks:
        raise FigureError(f"{path}: theory JSON has incomplete anchored envelope")
    center_masses = anchored.get("centerMasses") or []
    masses = [float(m) for m in center_masses]
    if len(masses) != len(peaks):
        # The true per-peak center masses are required: reconstructing a regular
        # anchor + i*NEUTRON grid drops heavy-isotope mass defects and silently
        # misplaces every peak. Refuse instead of drawing a wrong axis.
        raise FigureError(
            f"{path}: theory JSON missing per-peak centerMasses "
            f"({len(peaks)} peaks, {len(masses)} center masses)"
        )
    intensities = [float(peak["intensity"]) for peak in peaks]
    apex = max(intensities) or 1.0
    return float(anchor), masses, [value / apex for value in intensities]


def _render(plot_path: str, payload: dict) -> list[str]:
    """Emit the contract and one self-contained chart-forge HTML figure.

    The renderer never sees chemistry; the contract contains chart-ready numbers.
    """
    out_html, _ = emit_html("envelope-compare", payload, plot_path,
                            source_default="tdeasenexus-toolbox envelope artifacts")
    return [str(out_html)]


def _norm(peaks: list[dict]) -> tuple[list[float], list[float]]:
    """(mz list, intensity list) — observed normalized to its own max so both
    envelopes share one y-axis (0..1) in the mirror layout."""
    if not peaks:
        return [], []
    mz = [p["mz"] for p in peaks]
    raw = [p["intensity"] for p in peaks]
    hi = max(raw) or 1.0
    return mz, [v / hi for v in raw]


def _centroid(mz: list[float], intensity: list[float]) -> float | None:
    """Intensity-weighted mean m/z — where the shape's intensity is centred."""
    total = sum(intensity)
    if not mz or total <= 0:
        return None
    return sum(m * v for m, v in zip(mz, intensity)) / total


def _apex(mz: list[float], intensity: list[float]) -> float | None:
    """The m/z of the most intense peak, a descriptive coordinate only."""
    if not mz:
        return None
    return mz[max(range(len(mz)), key=lambda i: intensity[i])]


def _stats(mz: list[float], intensity: list[float]) -> dict:
    return {
        "apex": _apex(mz, intensity),
        "centroid": _centroid(mz, intensity),
    }


def _metrics_text(metrics: dict) -> str:
    parts = []
    for key, fmt in [
        ("pearson", lambda v: f"pearson {v:.3f}"),
        ("refinedMass", lambda v: f"refined {v:.4f} Da"),
        ("matchedPeaks", lambda v: f"{v} peaks"),
    ]:
        if key in metrics and metrics[key] is not None:
            parts.append(fmt(metrics[key]))
    return "  |  ".join(parts)


def _compute_pair_stats(
    obs: Envelope,
    candidates: list[Envelope],
    metrics_list: list[dict],
) -> dict:
    """Project source metrics and descriptive coordinates; never infer QC states."""

    obs_stats = _stats(obs.mz, obs.intensity)
    cand_stats = []
    for cand, metrics in zip(candidates, metrics_list):
        st = _stats(cand.mz, cand.intensity)
        st["name"] = cand.label
        st["pearson"] = metrics.get("pearson")
        st["refinedMass"] = metrics.get("refinedMass")
        st["matchedPeaks"] = metrics.get("matchedPeaks")
        st["delta_apex_mTh"] = (
            (obs_stats["apex"] - st["apex"]) * 1000.0
            if obs_stats["apex"] is not None and st["apex"] is not None
            else None
        )
        st["delta_centroid_mTh"] = (
            (obs_stats["centroid"] - st["centroid"]) * 1000.0
            if obs_stats["centroid"] is not None and st["centroid"] is not None
            else None
        )
        cand_stats.append(st)

    return {
        "observed": obs_stats,
        "candidates": cand_stats,
    }


def _comparison_payload(obs: Envelope, cands: list[Envelope], stats: dict,
                        title: str) -> dict:
    """The one-panel envelope-compare contract: observed sticks, palette-colored
    candidate series with per-candidate stat lines, apex registration."""
    return {
        "data": {
            "panels": [{
                "observed": {"x": obs.mz, "y": obs.intensity},
                "series": [{"label": c.label,
                            "color": CANDIDATE_PALETTE[i % len(CANDIDATE_PALETTE)],
                            "x": c.mz, "y": c.intensity,
                            "stat": _stat_text(stats["candidates"][i])}
                           for i, c in enumerate(cands)],
                "apex": {"observed": stats["observed"]["apex"],
                         "series": [st["apex"] for st in stats["candidates"]]},
            }],
            "xlabel": "m/z", "ylabel": "Relative intensity",
        },
        "meta": {"title": title, "n": len(obs.mz)},
    }


def build_figure(result: dict) -> dict:
    """Build normalized data for one refine-shift result."""
    obs_mz, obs_i = _norm(result["observed"])
    theo_mz, theo_i = _norm(result["theoretical"])
    metrics = {
        "command": result.get("command") or "refine-shift",
        "pearson": result.get("pearson"),
        "refinedMass": result.get("refinedMass"),
        "matchedPeaks": result.get("matchedPeaks"),
    }
    return {
        "figure": "envelope-mirror",
        "observed": [{"mz": m, "intensity": v} for m, v in zip(obs_mz, obs_i)],
        "theoretical": [{"mz": m, "intensity": v} for m, v in zip(theo_mz, theo_i)],
        "centroid": {"observed": _centroid(obs_mz, obs_i), "theoretical": _centroid(theo_mz, theo_i)},
        "metrics": metrics,
    }


def build_figures(results: list[dict]) -> dict:
    """Multi-refine figure data: observed (first result's, or shared) + every
    candidate's theoretical envelope and metrics."""
    observed = results[0]["observed"]
    if any(result["observed"] != observed for result in results[1:]):
        raise FigureError("OBSERVATION_MISMATCH: paired candidates must use the same observed envelope")
    candidates = []
    metrics_list = []
    for result in results:
        _, theo_i = _norm(result["theoretical"])
        theo_mz, _ = _norm(result["theoretical"])
        candidates.append({
            "label": result.get("label") or result.get("command") or f"candidate {len(candidates) + 1}",
            "theoretical": [{"mz": m, "intensity": v} for m, v in zip(theo_mz, theo_i)],
        })
        metrics_list.append({
            "command": result.get("command") or "refine-shift",
            "pearson": result.get("pearson"),
            "refinedMass": result.get("refinedMass"),
            "matchedPeaks": result.get("matchedPeaks"),
        })
    obs_mz, obs_i = _norm(observed)
    return {
        "figure": "envelope-paired",
        "observed": [{"mz": m, "intensity": v} for m, v in zip(obs_mz, obs_i)],
        "candidates": candidates,
        "metrics": metrics_list,
        "centroid": {"observed": _centroid(obs_mz, obs_i)},
    }

