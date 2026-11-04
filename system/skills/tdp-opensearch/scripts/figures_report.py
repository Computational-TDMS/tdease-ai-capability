"""One-page OpenSearch report with the statistical discovery view set."""

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from _shared.manifest import write_manifest
from figures_census import census_payload, shift_spectrum_payload
from figures_featuremap import featuremap_payload
from figures_proteoform import proteoform_payload


def smith_histogram_payload(observations):
    """Preserve source-record ambiguity grades independently of reference overlays."""
    counts = {level: 0 for level in range(1, 6)}
    total = 0
    for obs in observations:
        smith = obs.get("smith") or {}
        try:
            level = int(smith.get("best"))
        except (TypeError, ValueError):
            continue
        if level in counts:
            total += 1
            counts[level] += 1
    if not total:
        raise ValueError("bundle carries no graded assignments")
    levels = [{"level": level, "searched": counts[level]} for level in sorted(counts)]
    return {"data": {"levels": levels,
                     "series": [{"key": "searched", "label": "source-record grade"}]},
            "meta": {"title": "Smith proteoform grading",
                     "source": "opensearch.json",
                     "n": f"{total} graded assignments"}}


def report_main(argv=None):
    """Build one self-contained discovery report.

    Statistical by intent: mass shifts, MS1 families and the Smith grading
    histogram render by default from the ingest bundle; detailed proteoform
    views stay explicit opt-ins.
    """
    from _shared.report import compile_report, envelope, example_script

    ap = argparse.ArgumentParser(
        description="One self-contained OpenSearch report: mass shifts, "
                    "families and Smith grading statistics.")
    ap.add_argument("--out", "--plot", dest="out", required=True, help="report HTML path")
    ap.add_argument("--bundle", help="opensearch.json produced by ingest build")
    ap.add_argument("--census", help="clusters JSON -> mass-shift panel (no-bundle legacy path)")
    ap.add_argument("--marks", help="evidence marks JSON for --census")
    ap.add_argument("--annotations", help="stable modification annotations JSON; "
                                          "overlays are opt-in, never bundled by default")
    ap.add_argument("--mark-role", choices=("all", "routine", "deep"), default="all")
    ap.add_argument("--routine-limit", type=int)
    ap.add_argument("--xlim", nargs=2, type=float, metavar=("LO", "HI"))
    ap.add_argument("--checklist", help="optional shift checklist CSV")
    ap.add_argument("--include-shift-spectrum", action="store_true",
                    help="add the optional rounded 1-Da view")
    ap.add_argument("--families", help="MS1 families JSON -> featuremap panel (no-bundle legacy path)")
    ap.add_argument("--features", help="TopPICR features.csv -> proteoform panel")
    ap.add_argument("--accession", help="accession for the proteoform panel")
    args = ap.parse_args(argv)

    raw_inputs = [args.census, args.marks, args.checklist, args.annotations,
                  args.families]
    if args.bundle and any(raw_inputs):
        ap.error("--bundle cannot be combined with raw census/families/checklist "
                 "inputs; the bundle already carries those views")
    if bool(args.features) != bool(args.accession):
        ap.error("--features and --accession are required together")

    bundle = None
    if args.bundle:
        bundle = json.loads(Path(args.bundle).read_text(encoding="utf-8"))
        if bundle.get("schema") != "tdease.opensearch.v1":
            ap.error("--bundle is not a tdease.opensearch.v1 file")
    census_source = args.census or (bundle.get("census") if bundle else None)
    families_source = args.families or (
        {"families": bundle.get("families", [])} if bundle else None
    )
    # Automatic UniMod attributions color the attributed sticks (label in the
    # tooltip); scattering them as on-axis dots/labels stays an opt-in.
    annotations_source = args.annotations or (
        {"annotations": bundle.get("annotations", [])} if bundle else None
    )

    def build(pid, title, caption, example, available, build_payload):
        if not available:
            return {"id": pid, "title": title, "caption": caption,
                    "status": "skipped", "note": "input not provided"}
        try:
            payload = build_payload()
        except Exception as exc:
            return {"id": pid, "title": title, "caption": caption,
                    "status": "error", "note": str(exc)}
        return {"id": pid, "title": title, "caption": caption,
                "script": example_script(example),
                "data": payload["data"], "meta": payload.get("meta", {})}

    proteoform_pid = "proteoform"
    if args.accession:
        proteoform_pid = "proteoform-" + re.sub(r"[^a-z0-9-]+", "-",
                                                args.accession.lower()).strip("-")
    panels = [
        build("massshifts", "Mass shifts",
              "Unexpected shift clusters observed after fixed forms; evidence "
              "marks are overlays, not observations.",
              "opensearch-census", bool(census_source),
              lambda: census_payload(census_source, args.marks, annotations_source,
                                     args.mark_role, args.routine_limit, args.xlim)[0]),
        build("featuremap", "MS1 feature families",
              "MS1 isotope families located by RT window and isotope spacing; area is member count and colour indicates MS2 support.",
              "ms1-family-map", bool(families_source),
              lambda: featuremap_payload(families_source)),
        build("grading-levels", "Smith proteoform grading",
              "Source-record ambiguity grades. Automatic mass-reference proposals "
              "do not remove identity or localization uncertainty.",
              "grading-levels", bool(bundle and bundle.get("observations")),
              lambda: smith_histogram_payload(bundle.get("observations"))),
    ]
    if args.include_shift_spectrum and args.checklist:
        panels.append(build("shift-spectrum", "Mass shifts (1 Da bins)",
                            "Optional rounded view for paper comparison; the precise "
                            "mass-shift panel above remains authoritative.",
                            "opensearch-census", True,
                            lambda: shift_spectrum_payload(args.checklist)))
    if args.features:
        panels.append(
            build(proteoform_pid, f"Proteoform map ({args.accession})",
                  "Sequence intervals with modification boxes for the accession; "
                  "the focused per-protein view.",
                  "proteoform-map", True,
                  lambda: proteoform_payload(args.features, args.accession, annotations_source)),
        )

    out = compile_report(
        panels, args.out, title="TDP OpenSearch report",
        meta_rows=[("Generated by", Path(__file__).name),
                   ("Inputs", ", ".join(filter(None, [
                       args.bundle, args.census, args.checklist,
                       args.families, args.features]))) or "none"])
    write_manifest(
        out, figure="opensearch_report", script=Path(__file__).name,
        params={"mark_role": args.mark_role},
        inputs=[p for p in (args.bundle, args.census, args.marks, args.annotations,
                            args.checklist, args.families, args.features) if p],
        caption={"meaning": "OpenSearch discovery report: mass shifts, family "
                            "context and Smith grading statistics."},
        plot_data=[],
    )
    print(json.dumps(envelope(out, panels), ensure_ascii=False))
    return 0
