"""Census panels: cluster sticks + evidence marks, and the 1-Da-binned
shift spectrum (both render the opensearch-census template)."""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from _shared.manifest import write_manifest
from figures_common import _load_json, _render, _emit


def auto_xlim(clusters: list[dict], *, base: int = 5,
              min_fraction: float = 0.005) -> tuple[list[float], int]:
    """Data-driven census axis frame: every shift observed at least
    ``max(base, min_fraction * top-cluster count)`` times; evidence marks are
    force-kept by the caller. Everything below the bar stays in the companion
    data CSV: the axis is a view, the CSV is the record."""
    if not clusters:
        return [-50.0, 50.0], base
    floor = max(base, math.ceil(min_fraction * max(int(c["count"]) for c in clusters)))
    kept = [float(c["center_mass"]) for c in clusters if int(c["count"]) >= floor]
    if not kept:
        return [-50.0, 50.0], floor
    lo, hi = min(kept), max(kept)
    pad = max(10.0, (hi - lo) * 0.03)
    return [lo - pad, hi + pad], floor


def census_payload(census_path, marks_path=None,
                   annotations_path: str | None = None, mark_role: str = "all",
                   routine_limit: int | None = None,
                   xlim: list[float] | None = None):
    """Build the census contract; returns (payload, plot_data) — plot_data is
    the manifest CSV record (census + attached evidence per cluster)."""
    census, source_name = _load_json(census_path)
    marks, _ = _load_json(marks_path) if marks_path else ({"marks": []}, "")
    annotations, _ = (_load_json(annotations_path)
                      if annotations_path else ({"annotations": []}, ""))
    clusters = census["clusters"]
    if not clusters:
        raise SystemExit("census has no clusters")

    shifts = [float(c["center_mass"]) for c in clusters]
    selected_marks = [
        mark for mark in marks["marks"]
        if mark_role == "all" or mark.get("role") == mark_role
    ]
    if routine_limit is not None:
        routine = [mark for mark in selected_marks if mark.get("role") == "routine"]
        selected_marks = [mark for mark in selected_marks
                          if mark.get("role") != "routine"] + routine[:routine_limit]
    selected_marks.extend(
        mark for mark in annotations.get("annotations", [])
        if mark_role == "all" or mark.get("role") == mark_role
    )
    marks_by_index = defaultdict(list)
    for mark in selected_marks:
        index = min(range(len(shifts)), key=lambda item: abs(shifts[item] - mark["mass"]))
        if abs(shifts[index] - mark["mass"]) > 1.5:
            continue
        marks_by_index[index].append(mark)

    def flat(mark: dict, key: str) -> str:
        """Marks files carry these fields as str or list; the CSV row wants one str."""
        value = mark.get(key, "")
        if isinstance(value, list):
            return " > ".join(str(part) for part in value)
        return str(value)

    plot_data = []
    contract_clusters = []
    for index, cluster in enumerate(clusters):
        attached = marks_by_index.get(index, [])
        automatic_labels = [mark.get("label", "") for mark in attached
                            if mark.get("role") == "automatic" and mark.get("label")]
        plot_data.append({
            "mass_shift_da": cluster["center_mass"],
            "prsm_count": cluster["count"],
            "evidence_role": "; ".join(sorted({mark.get("role", "") for mark in attached} - {""})),
            "evidence_label": "; ".join(mark.get("label", "") for mark in attached if mark.get("label")),
            "verdict": "; ".join(sorted({mark.get("verdict", "") for mark in attached} - {""})),
            "basis": " | ".join(b for mark in attached if (b := flat(mark, "basis"))),
            "sources": " | ".join(s for mark in attached if (s := flat(mark, "sources"))),
        })
        contract_clusters.append({
            **{k: v for k, v in cluster.items() if k != "center_mass"},
            "mass": cluster["center_mass"],
            "count": cluster["count"],
            "annotation": " / ".join(automatic_labels),
            "mark_labels": [mark.get("label", "") for mark in attached if mark.get("label")],
            "mark_colors": [mark.get("color", "") for mark in attached if mark.get("color")],
            "mark_verdicts": [mark.get("verdict", "") for mark in attached if mark.get("verdict")],
        })
    auto_range, xlim_floor = auto_xlim(clusters)
    if xlim:
        x_range = [float(xlim[0]), float(xlim[1])]
        xlim_mode = "manual"
    else:
        x_range = auto_range
        marked = [shifts[index] for index in marks_by_index]
        if marked:
            x_range[0] = min(x_range[0], min(marked) - 15.0)
            x_range[1] = max(x_range[1], max(marked) + 15.0)
        xlim_mode = f"auto: shifts with >= {xlim_floor} observations"
    for cluster in contract_clusters:
        cluster["inview"] = x_range[0] <= cluster["mass"] <= x_range[1]

    return {"data": {"clusters": contract_clusters,
                     "marks": [m for m in selected_marks
                               if any(abs(m["mass"] - c["mass"]) < 0.01
                                      for c in contract_clusters if c["inview"])],
                     "xlim": x_range, "floor": xlim_floor},
            "meta": {"title": "Mass shifts",
                     "source": source_name,
                     "thresholds": xlim_mode,
                     "n": f"{sum(c['inview'] for c in contract_clusters)} clusters on axis"}}, plot_data


def census_main(argv=None):
    ap = argparse.ArgumentParser(description="Observed mass-shift census with evidence marks.")
    ap.add_argument("--census", required=True, help="clusters JSON from mass_shift_census.py")
    ap.add_argument("--marks", help="optional evidence marks JSON for this census")
    ap.add_argument("--annotations",
                    help="optional stable modification annotations JSON; labels are overlays, not observations")
    ap.add_argument("--mark-role", choices=("all", "routine", "deep"), default="all",
                    help="which evidence marks to project (default: all)")
    ap.add_argument("--routine-limit", type=int,
                    help="project only the first N routine marks; deep marks are kept")
    ap.add_argument("--out", "--plot", dest="out", required=True, help="output HTML path")
    ap.add_argument("--xlim", nargs=2, type=float, metavar=("LO", "HI"),
                    help="restrict the mass axis")
    args = ap.parse_args(argv)

    payload, plot_data = census_payload(args.census, args.marks, args.annotations,
                                        args.mark_role, args.routine_limit, args.xlim)
    written = _render("opensearch-census", payload, args.out)
    write_manifest(
        args.out,
        figure="mass_shift_census",
        script=Path(__file__).name,
        params={"mark_role": args.mark_role, "routine_limit": args.routine_limit,
                "xlim": payload["data"]["xlim"]},
        inputs=[path for path in (args.census, args.marks, args.annotations) if path],
        caption={
            "cluster_count": len(payload["data"]["clusters"]),
            "meaning": "Observed mass-shift census; evidence marks are overlays, not source observations.",
        },
        plot_data=plot_data,
    )
    _emit("census", written + [str(Path(args.out).with_suffix(".manifest.json"))])


SPECTRUM_MARK_CYCLE = ("orange", "green")   # adjacent marks alternate for legibility


def load_rows(path):
    with open(path, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["mass_shift_id"]].append(row)
    output = []
    for mass_shift_id, members in grouped.items():
        output.append({
            "mass_shift_id": mass_shift_id,
            "shift_mass": float(members[0]["shift_mass"]),
            "prsm_count": sum(int(row["prsm_count"]) for row in members),
        })
    return sorted(output, key=lambda row: row["shift_mass"])


def spectrum_bins(rows):
    """Paper convention: sum PrSM counts after rounding shifts to 1 Da."""
    bins = defaultdict(int)
    for row in rows:
        bins[round(row["shift_mass"])] += row["prsm_count"]
    return [
        {"shift_mass_rounded": mass, "prsm_count": bins[mass]}
        for mass in sorted(bins)
    ]


def annotation_rows(rows, annotate_ids):
    selected = [row for row in rows if row["mass_shift_id"] in annotate_ids]
    return list({row["mass_shift_id"]: row for row in selected}.values())


def shift_spectrum_payload(checklist_path: str, annotate_ids=None) -> dict:
    """Build the 1-Da-binned shift-census contract (figure + report)."""
    rows = load_rows(checklist_path)
    bins = spectrum_bins(rows)
    by_mass = {row["shift_mass_rounded"]: row for row in bins}
    selected = annotation_rows(rows, set(annotate_ids or []))
    mark_rows = []
    for index, row in enumerate(sorted(selected, key=lambda item: item["shift_mass"])):
        mass = round(row["shift_mass"])
        if mass not in by_mass:
            continue
        mark_rows.append({"mass": float(mass),
                          "label": f"{row['shift_mass']:+.2f}",
                          "role": "routine",
                          "color": SPECTRUM_MARK_CYCLE[index % 2]})

    x_range, floor = auto_xlim([{"center_mass": b["shift_mass_rounded"], "count": b["prsm_count"]}
                                for b in bins])
    return {"data": {
        "clusters": [{"mass": float(b["shift_mass_rounded"]), "count": b["prsm_count"],
                      "annotation": "",
                      "mark_labels": [], "mark_colors": [], "mark_verdicts": [],
                      "inview": x_range[0] <= b["shift_mass_rounded"] <= x_range[1]}
                     for b in bins],
        "marks": mark_rows, "xlim": x_range, "floor": floor},
        "meta": {"title": "Mass shifts (1 Da bins)",
                 "source": Path(checklist_path).name,
                 "n": f"{len(bins)} bins"}}


def shift_spectrum_main(argv=None):
    parser = argparse.ArgumentParser(description="Shift census with evidence marks (rounded 1 Da bins).")
    parser.add_argument("--checklist", required=True)
    parser.add_argument("--out", "--plot", dest="out", required=True, help="output HTML path")
    parser.add_argument("--annotate-mass-shift-id", action="append", default=[])
    args = parser.parse_args(argv)

    payload = shift_spectrum_payload(args.checklist, args.annotate_mass_shift_id)
    written = _render("opensearch-census", payload, args.out)
    _emit("shift-spectrum", written + [str(Path(args.out).with_suffix(".manifest.json"))])

