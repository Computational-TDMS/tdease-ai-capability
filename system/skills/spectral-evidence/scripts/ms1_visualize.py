"""MS1 envelope figures — reached through ``spectral.py``, never directly.

Two routes, both drawing one shared overlay: the observed MS1 envelope as gray
centroid sticks, with any number of candidate theoretical envelopes as colored
lines on the same m/z axis. Detailed fit metrics remain in the manifest.

- :func:`run_theory_overlay` — one observed feature envelope against N complete
  models' theoretical envelopes (the composition flow).
- :func:`render_refine` — refine-shift's fitted precursor envelope against its
  own isotope-offset evaluations, or several candidates against one observation.

The feature-ladder route lives in ``envelope_ladder`` and is what
``spectral.py measure --ms1`` calls.

Apex and centroid differences are descriptive coordinates, not QC or chemical
decisions. Fit metrics come from the source artifact; these renderers apply no
calibration, isotope-registration or spacing acceptance threshold.

The figures never scan integer isotope offsets: candidate models remain fixed
against the one observed physical envelope.

Rendered as self-contained D3 HTML through chart-forge. PNG/SVG export remains
an explicit action from the figure page.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from _shared.manifest import write_manifest
from _shared.kernel_io import load_doc, load_result, load_theory_envelope, result_of
from envelope_common import (
    Envelope,
    FigureError,
    _comparison_payload,
    _compute_pair_stats,
    _evaluation_candidates,
    _load_refine,
    _ranked_evaluations,
    _render,
    build_figure,
    build_figures,
)


def _bound_fit_metrics(doc: dict, fit_path: str, evaluation: dict, observed_path: str,
                       theory_path: str, theory_charge) -> dict:
    """Bind one fit-envelope artifact to the observed/theory pair it claims.

    Positional order states intent; identity is verified from the artifact's
    own recorded provenance (``params.observed``/``params.theoretical`` paths
    and ``result.charge``). A definite mismatch — the fit names a different
    readable input, or a different charge — is rejected rather than presented
    as this candidate's score. Missing provenance degrades to an explicit
    ``unverified`` marker instead of borrowed numbers passing silently.
    """
    result = doc.get("result") or {}
    params = doc.get("params") if isinstance(doc.get("params"), dict) else {}
    fit_charge = result.get("charge")
    if theory_charge and fit_charge and fit_charge != theory_charge:
        raise FigureError(
            f"FIT_CHARGE_MISMATCH: {fit_path} was measured at z{fit_charge} but "
            f"{theory_path} is z{theory_charge}")
    reasons = []

    def check_input(role: str, recorded, actual: str) -> None:
        if not recorded or not isinstance(recorded, str):
            reasons.append(f"fit params record no {role} provenance")
            return
        if Path(recorded).resolve() == Path(actual).resolve():
            return
        if Path(recorded).is_file():
            raise FigureError(
                f"FIT_{role.upper()}_MISMATCH: {fit_path} was measured from "
                f"{recorded}, not from {actual}")
        reasons.append(
            f"fit {role} reference {recorded} is not readable; "
            f"cannot confirm it is {actual}")

    check_input("observed", params.get("observed"), observed_path)
    check_input("theoretical", params.get("theoretical"), theory_path)
    metrics = {
        "pearson": evaluation.get("pearson"),
        "matchedPeaks": evaluation.get("matchedPeaks"),
        "coverage": evaluation.get("coverage"),
        "verification": "verified" if not reasons else "unverified",
    }
    if reasons:
        metrics["verificationReasons"] = reasons
    return metrics


def run_theory_overlay(args) -> dict:
    """Composition flow: one observed feature envelope vs N complete-model
    theory envelopes (optionally with fit-envelope metrics), rendered as the
    SAME interactive envelope-compare figure every MS1 view uses.

    Returns its stdout payload; ``spectral.py show --ms1`` owns the printing."""
    if not args.theory:
        raise FigureError("--observed needs at least one --theory JSON "
                          "(one per candidate, same order as --fit)")
    observed = load_result(args.observed)
    if not observed.get("peaks"):
        raise FigureError(f"{args.observed}: expected an envelope-extract result "
                          "with a peaks array")
    obs = Envelope([p["mz"] for p in observed["peaks"]],
                   [p["intensity"] for p in observed["peaks"]], "observed")
    if len(args.fit) > len(args.theory):
        raise FigureError("at most one --fit per --theory, repeated in the same order")

    cands, metrics_list, theory_charges = [], [], []
    for index, theory_path in enumerate(args.theory):
        envelope = load_theory_envelope(theory_path)
        peaks = envelope.get("peaks")
        if not peaks:
            raise FigureError(f"{theory_path}: no theory envelope peaks")
        charge = envelope.get("charge")
        theory_charges.append(charge)
        if observed.get("charge") and charge and charge != observed["charge"]:
            raise FigureError(
                f"CHARGE_MISMATCH: {theory_path} is z{charge} but the observed "
                f"feature is z{observed['charge']}")
        label = (args.theory_label[index] if index < len(args.theory_label)
                 else Path(theory_path).stem.removeprefix("theory_"))
        cands.append(Envelope([p["mz"] for p in peaks],
                              [p["intensity"] for p in peaks], label))
        metrics_list.append({})
    for index, fit_path in enumerate(args.fit):
        doc = load_doc(fit_path)
        evaluation = result_of(doc, source=str(fit_path)).get("evaluation") or {}
        metrics_list[index] = _bound_fit_metrics(
            doc, fit_path, evaluation, args.observed,
            args.theory[index], theory_charges[index])

    stats = _compute_pair_stats(obs, cands, metrics_list)
    plot_path = args.plot or (Path(args.observed).stem + "-envelope.png")
    _render(plot_path, _comparison_payload(
        obs, cands, stats, args.label or "MS1 envelope comparison"))
    out = {
        "figure": str(Path(plot_path).with_suffix(".html").resolve()),
        "mediaType": "text/html",
        "renderer": "chart-forge d3",
        "observedPeaks": len(obs.mz),
        "candidateCount": len(cands),
        "candidates": stats["candidates"],
        "fits": metrics_list,
        "shapeMetric": "Pearson intensity-shape correlation; not an identity score",
    }
    write_manifest(
        plot_path,
        figure="envelope-theory-overlay",
        script=Path(__file__).name,
        params={"observed": args.observed,
                "theories": list(args.theory), "fits": list(args.fit)},
        inputs=[args.observed, *args.theory, *args.fit],
        caption={"apex_delta_mTh": stats["candidates"][0].get("delta_apex_mTh")
                 if stats["candidates"] else None},
        plot_data=stats["candidates"],
    )
    return out


def render_refine(args) -> dict:
    """Render refine-shift results: the observed precursor envelope against one
    candidate's theoretical envelope (with its isotope-offset evaluations), or
    several candidates against the same observation.

    Returns its stdout payload; ``spectral.py show --ms1`` owns the printing."""
    results = [_load_refine(r) for r in args.refine]
    charges = {result.get("_params", {}).get("charge") for result in results}
    if len(charges) != 1:
        raise FigureError("CHARGE_MISMATCH: paired candidates must use the same precursor charge")
    charge = next(iter(charges))
    for result, refine_path in zip(results, args.refine):
        result["_path"] = refine_path
        if args.label and len(results) == 1:
            result["label"] = args.label
        elif len(results) > 1:
            # Derive a short per-candidate label from the file stem.
            result["label"] = Path(refine_path).stem.replace("refine-", "").replace("_refine", "")

    if len(results) == 1:
        result = results[0]
        fig = build_figure(result)
        obs = Envelope([p["mz"] for p in fig["observed"]], [p["intensity"] for p in fig["observed"]],
                       "observed")
        theo = Envelope([p["mz"] for p in fig["theoretical"]], [p["intensity"] for p in fig["theoretical"]],
                        fig.get("_label", results[0].get("label", "theoretical (candidate)")))
        # All isotope-offset candidates on one axis (the paper's dual-hypothesis
        # read): best-pearson candidate leads in orange, runners-up follow the
        # palette. Falls back to the single refined envelope for refine JSONs
        # without a multi-offset evaluation list.
        cands = _evaluation_candidates(result, theo, result.get("refinedMass"), charge)
        metrics_list = [{"pearson": ev.get("pearson"), "refinedMass": ev.get("refinedMass"),
                         "matchedPeaks": ev.get("matchedPeaks")}
                        for ev in _ranked_evaluations(result)]
        if not cands:
            cands = [theo]
            metrics_list = [fig["metrics"]]
        plot_path = args.plot or (Path(args.refine[0]).stem + "-envelope.png")
        stats = _compute_pair_stats(obs, cands, metrics_list)
        _render(plot_path, _comparison_payload(
            obs, cands, stats, args.label or "MS1 envelope comparison"))
        out = {
            "figure": fig["figure"],
            "plot": str(Path(plot_path).with_suffix(".html").resolve()),
            "mediaType": "text/html",
            "renderer": "chart-forge d3",
            "observedPeaks": len(fig["observed"]),
            "candidateCount": len(cands),
            "theoreticalPeaks": len(fig["theoretical"]),
            "centroid": fig["centroid"],
            "metrics": fig["metrics"],
            "apex": stats["observed"],
            "shapeMetric": "Pearson intensity-shape correlation; not an identity score",
        }
        if args.plot_data:
            out["data"] = {"observed": fig["observed"], "theoretical": fig["theoretical"]}
        write_manifest(
            plot_path,
            figure="envelope-mirror",
            script=Path(__file__).name,
            params={"mode": args.mode, "label": args.label,
                    "candidateCount": len(cands)},
            inputs=[results[0].get("_path", "")],
            caption={"apex_delta_mTh": stats["candidates"][0].get("delta_apex_mTh") if stats["candidates"] else None},
            plot_data=stats["candidates"],
        )
    else:
        fig = build_figures(results)
        obs = Envelope([p["mz"] for p in fig["observed"]], [p["intensity"] for p in fig["observed"]],
                       "observed")
        cands = []
        for cand in fig["candidates"]:
            cands.append(Envelope([p["mz"] for p in cand["theoretical"]],
                                  [p["intensity"] for p in cand["theoretical"]],
                                  cand["label"]))
        plot_path = args.plot or (Path(args.refine[0]).stem + "-envelope.png")
        stats = _compute_pair_stats(obs, cands, fig["metrics"])
        _render(plot_path, _comparison_payload(
            obs, cands, stats, args.label or "MS1 envelope comparison"))
        out = {
            "figure": fig["figure"],
            "plot": str(Path(plot_path).with_suffix(".html").resolve()),
            "mediaType": "text/html",
            "renderer": "chart-forge d3",
            "observedPeaks": len(fig["observed"]),
            "candidateCount": len(cands),
            "candidatePeaks": [len(c["theoretical"]) for c in fig["candidates"]],
            "centroid": fig["centroid"],
            "metrics": fig["metrics"],
            "candidates": stats["candidates"],
            "shapeMetric": "Pearson intensity-shape correlation; not an identity score",
        }
        if args.plot_data:
            out["data"] = {
                "observed": fig["observed"],
                "candidates": [{"label": c["label"], "theoretical": c["theoretical"]} for c in fig["candidates"]],
            }
        write_manifest(
            plot_path,
            figure="envelope-paired",
            script=Path(__file__).name,
            params={"mode": args.mode, "label": args.label},
            inputs=[r.get("_path", "") for r in results],
            caption={},
            plot_data=stats["candidates"],
        )
    return out
