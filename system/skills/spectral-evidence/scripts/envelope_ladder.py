"""Feature-ladder mode: RT-window-aggregated observed feature envelopes
against candidate theoretical series (the only kernel-touching mode)."""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from _shared.manifest import write_manifest
from _shared.kernel_io import result_of
from envelope_common import (
    CANDIDATE_PALETTE,
    FigureError,
    _load_theory_profile,
    _render,
)


def _feature_ladder_contract(
    mzml: str,
    features: list[list[str]],
    *,
    cli: str | None,
    theories: list[list[str]],
    theory_labels: list[str],
    tolerance_ppm: float,
    isotopes: int,
    artifact_stem: Path,
) -> tuple[dict, list[dict], list[str], str | None]:
    """Observed-envelope ladder contract with any number of theoretical
    hypotheses. Single column, one panel per feature (1x/2x/3x ...). Each panel
    overlays the observed MS1 envelope (gray sticks) with each candidate's
    theoretical isotope envelope as a colored line. The shared axis is the mass
    offset from the feature's neutral mass (Da): observed m/z bins and neutral
    theory masses only meet when both are mapped there — an isotope-index axis
    would silently misalign coarse observed bins against theory peaks."""
    from _shared.native import binary as ms_binary, mass_constants, run_json

    n_features = len(features)
    n_candidates = len(theories)
    if n_candidates and len(theory_labels) != n_candidates:
        raise FigureError("THEORY_LABEL_MISMATCH: provide one --theory-label per candidate")
    for index, candidate_paths in enumerate(theories):
        if len(candidate_paths) != n_features:
            raise FigureError(
                "THEORY_COUNT_MISMATCH: candidate "
                f"{index} has {len(candidate_paths)} JSONs; expected {n_features} (one per feature)"
            )
    panels = []
    rows = []
    # one mzML traversal for the whole ladder: the kernel batches every window
    # (per-request matching is exactly the single-window command's)
    command = [
        ms_binary(cli), "envelope-windows", mzml,
        "--windows", json.dumps([
            {"label": label, "mass": float(mass), "charge": int(charge),
             "rtMin": float(rt_min), "rtMax": float(rt_max)}
            for label, mass, charge, rt_min, rt_max in features
        ]),
        "--tolerance-ppm", str(tolerance_ppm), "--isotopes", str(isotopes),
    ]
    try:
        payload = run_json(command, label="tdeasenexus-toolbox envelope-windows")
    except RuntimeError as error:
        raise FigureError(str(error)) from error
    windows = result_of(payload, source="tdeasenexus-toolbox envelope-windows")["windows"]
    artifact_stem.parent.mkdir(parents=True, exist_ok=True)
    observations = []
    pairs = []
    for i, window in enumerate(windows):
        if window.get("error"):
            raise FigureError(str(window["error"]))
        if not window.get("peaks"):
            raise FigureError(f"NO_OBSERVED_ENVELOPE: {features[i][0]}")
        observed = artifact_stem.with_name(f"{artifact_stem.name}-observed-{i}.json")
        observed.write_text(json.dumps({**window, "sourceFile": mzml}), encoding="utf-8")
        observations.append(str(observed))
        for j, paths in enumerate(theories):
            pairs.append({"id": f"f{i}-m{j}", "observed": str(observed),
                          "theoretical": paths[i], "matchWindowPpm": tolerance_ppm,
                          "minimumWindowMz": 0})
    fit_path = artifact_stem.with_name(artifact_stem.name + "-fits.json")
    fits = {}
    if pairs:
        fitted = run_json([ms_binary(cli), "fit-envelope", "--pairs", json.dumps(pairs)], label="batch precursor-envelope fits")
        fit_path.write_text(json.dumps(fitted, ensure_ascii=False), encoding="utf-8")
        fits = {row["id"]: row["evaluation"] for row in fitted["result"]["fits"]}
    for member_index, (label, mass, charge, rt_min, rt_max) in enumerate(features):
        window = windows[member_index]
        # The batch contract is positional; the kernel echoes each request's
        # label, so verify it — a reordered response must never silently
        # attach one feature's aggregated envelope to another feature's panel.
        if window.get("label") != label:
            raise FigureError(
                f"WINDOW_ORDER_MISMATCH: envelope-windows returned label "
                f"{window.get('label')!r} at position {member_index}; expected {label!r}")
        if window.get("error"):
            raise FigureError(str(window["error"]))
        peaks = window.get("peaks") or []
        if not peaks:
            raise FigureError(f"NO_OBSERVED_ENVELOPE: {label}")
        intensities = [float(peak["intensity"]) for peak in peaks]
        apex = max(intensities) or 1.0
        z = int(charge)
        proton = mass_constants(cli)["protonMass"]
        mass_offset = lambda mz: mz * z - z * proton - float(mass)  # noqa: E731
        row = {
            "label": label,
            "mass": float(mass),
            "charge": z,
            "rt_min": float(rt_min),
            "rt_max": float(rt_max),
            "scan_count": int(window.get("scanCount", 0)),
            "observed_offset": [mass_offset(float(peak["mz"])) for peak in peaks],
            "observed": [value / apex for value in intensities],
            "theories": [],
        }
        series = []
        for candidate_index, candidate_paths in enumerate(theories):
            anchor_mass, center_masses, probabilities = _load_theory_profile(candidate_paths[member_index])
            # The theory's anchor records which model mass it was computed
            # for; surfacing the offset to the declared feature mass makes a
            # mis-ordered --theory list visible instead of silently rendered
            # as a plausible-looking comparison.
            anchor_offset = float(anchor_mass) - float(mass)
            theory_apex = max(probabilities) or 1.0
            row["theories"].append({
                "path": candidate_paths[member_index],
                "label": theory_labels[candidate_index],
                "anchor_mass": float(anchor_mass),
                "anchor_offset_da": anchor_offset,
                "x": [
                    # True per-peak center mass, not a regular anchor + i*NEUTRON
                    # grid: heavy-isotope mass defects (e.g. Br, ~0.998 Da per
                    # atom) would otherwise shift every peak off by ~1 isotope
                    # index and compound with multiplicity.
                    center_mass - float(mass)
                    for center_mass in center_masses
                ],
                "normalized_intensity": [value / theory_apex for value in probabilities],
            })
            tx = row["theories"][-1]["x"]
            ty = row["theories"][-1]["normalized_intensity"]
            evaluation = fits[f"f{member_index}-m{candidate_index}"]
            r = evaluation.get("pearson")
            row["theories"][-1]["matched_peaks"] = evaluation["matchedPeaks"]
            apex_delta = None
            if ty and row["observed"]:
                apex_delta = (tx[ty.index(max(ty))]
                              - row["observed_offset"][row["observed"].index(max(row["observed"]))])
            row["theories"][-1]["pearson_r"] = r
            row["theories"][-1]["apex_delta_da"] = apex_delta
            stat = f"r={r:.3f}" if r is not None else "r n/a"
            if apex_delta is not None:
                stat += f" · apex {apex_delta:+.2f} Da"
            series.append({
                "label": theory_labels[candidate_index],
                "color": CANDIDATE_PALETTE[candidate_index % len(CANDIDATE_PALETTE)],
                "x": tx,
                "y": ty,
                "stat": stat,
            })
        rows.append(row)

        def apex_index(values: list[float]) -> int | None:
            return values.index(max(values)) if values else None

        panels.append({
            "title": f"{label} · {float(mass):.4f} Da · z={int(charge)}",
            "observed": {"x": row["observed_offset"], "y": row["observed"]},
            "series": series,
            "apex": {"observed": apex_index(row["observed"]),
                     "series": [apex_index(s["y"]) for s in series]},
        })

    payload = {
        "data": {
            "panels": panels,
            "xlabel": "Mass offset from feature mass (Da)",
            "ylabel": "Normalized intensity",
        },
        "meta": {
            "title": "RT-window-aggregated precursor envelopes",
            "source": mzml,
            "n": f"{len(rows)} feature panels",
        },
    }
    return payload, rows, observations, str(fit_path) if pairs else None



def ladder(args) -> dict:
    """Build and render the feature ladder; return the stdout payload.

    The caller emits it — one JSON object per process, whether the ladder was
    reached through ``spectral.py measure --ms1`` or through the retired
    ms1_visualize entry.
    """
    if not args.feature:
        raise FigureError("FEATURE_REQUIRED: provide at least one --feature")
    n_features = len(args.feature)
    if args.theory:
        if len(args.theory) % n_features != 0:
            raise FigureError(
                "THEORY_COUNT_MISMATCH: expected --theory entries in multiples of "
                f"the number of features ({n_features}); got {len(args.theory)}"
            )
        n_candidates = len(args.theory) // n_features
        if n_candidates and len(args.theory_label) != n_candidates:
            raise FigureError("THEORY_LABEL_MISMATCH: provide one --theory-label per candidate")
        theories = [
            args.theory[candidate * n_features:(candidate + 1) * n_features]
            for candidate in range(n_candidates)
        ]
    else:
        theories = []
    plot_path = args.plot or "observed-feature-envelopes.png"
    payload, rows, observations, fits = _feature_ladder_contract(
        args.feature_ladder, args.feature,
        cli=args.cli, theories=theories, theory_labels=args.theory_label,
        tolerance_ppm=args.tolerance_ppm, isotopes=args.isotopes,
        artifact_stem=Path(plot_path).with_suffix(""),
    )
    if args.label:
        payload["meta"]["title"] = args.label
    render = not getattr(args, "no_figures", False)
    written = (_render(plot_path, payload) if render else []) + observations + ([fits] if fits else [])
    caption = [{
        "label": row["label"], "mass": row["mass"], "charge": row["charge"],
        "rt_min": row["rt_min"], "rt_max": row["rt_max"], "scan_count": row["scan_count"],
        # per-candidate identity binding and statistics, flat in stdout: the
        # theory's own anchor mass vs this feature's declared mass, plus the
        # comparison numbers — auditable without opening the figure
        "theories": [{
            "label": theory["label"],
            "anchor_mass": theory["anchor_mass"],
            "anchor_offset_da": theory["anchor_offset_da"],
            "pearson_r": theory["pearson_r"],
            "matched_peaks": theory["matched_peaks"],
            "apex_delta_da": theory["apex_delta_da"],
        } for theory in row["theories"]],
    } for row in rows]
    if render:
        write_manifest(
            plot_path, figure="feature-envelope-ladder", script=Path(__file__).name,
            params={"features": args.feature, "theoryLabels": args.theory_label,
                    "tolerancePpm": args.tolerance_ppm, "isotopes": args.isotopes},
            inputs=[args.feature_ladder, *args.theory, *observations, *([fits] if fits else [])], caption={"members": caption},
            plot_data=rows,
        )
    return {
        "figure": "feature-envelope-ladder",
        "plot": str(Path(plot_path).with_suffix(".html").resolve()) if render else None,
        "fits": fits, "observations": observations,
        "renderer": "chart-forge d3", "members": caption, "written": written,
    }


def run(args) -> int:
    print(json.dumps(ladder(args), ensure_ascii=False))
    return 0
