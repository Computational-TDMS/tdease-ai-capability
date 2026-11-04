"""MS1 feature-family map panel (ms1-family-map template)."""

import argparse
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from _shared.manifest import write_manifest
from figures_common import _load_json, _render, _emit


def featuremap_payload(families_path) -> dict:
    """Build the flat family-row contract consumed by ms1-family-map."""
    payload_in, source_name = _load_json(families_path)
    families = payload_in.get("families", [])
    rows = []
    for family in families:
        members = family.get("members") or family.get("member_records") or []
        feature_ids = family.get("feature_ids") or [item.get("feature_id")
                                                      for item in members
                                                      if item.get("feature_id") is not None]
        masses = family.get("masses") or [item.get("mass") for item in members
                                           if item.get("mass") is not None]
        rt_apexes = family.get("rt_apexes") or family.get("rtApexes") or [
            item.get("rt_apex") for item in members if item.get("rt_apex") is not None
        ]
        start = family.get("common_rt_start")
        end = family.get("common_rt_end")
        if start is None:
            start = min(rt_apexes) if rt_apexes else None
        if end is None:
            end = max(rt_apexes) if rt_apexes else None
        support = family.get("ms2_support")
        if support is None:
            support = "yes" if (family.get("prsm_count") or 0) > 0 else "none"
        rows.append({
            "family_id": str(family.get("family_id", "")),
            "spacing_center_da": family.get("spacing_center_da"),
            "spacing_spread_da": family.get("spacing_spread_da"),
            "member_count": family.get("member_count", len(masses)),
            "feature_ids": "; ".join(str(value) for value in feature_ids),
            "masses": "; ".join(str(value) for value in masses),
            "mass_gaps": "; ".join(str(value) for value in (
                family.get("mass_gaps") or [])),
            "gap_error_ppm_intact_max": family.get("gap_error_ppm_intact_max"),
            "common_rt_start": start,
            "common_rt_end": end,
            "intensity_total": family.get("intensity_total"),
            "ms2_support": support,
        })
    return {"data": rows,
            "meta": {"title": "MS1 feature families",
                     "source": source_name,
                     "transform": "one row per discovered MS1 isotope family; RT bounds are retained",
                     "n": len(rows)}}


def featuremap_main(argv=None):
    ap = argparse.ArgumentParser(description="Plot MS1 feature-family ladders from OpenSearch JSON.")
    ap.add_argument("--families", required=True)
    ap.add_argument("--out", "--plot", dest="out", required=True, help="output HTML path")
    args = ap.parse_args(argv)
    payload = featuremap_payload(args.families)
    written = _render("ms1-family-map", payload, args.out)
    write_manifest(
        args.out, figure="ms1_feature_families", script=Path(__file__).name,
        params={}, inputs=[args.families],
        caption={"meaning": "One row per MS1 isotope family; position is RT and isotope spacing, area is member count, colour indicates MS2 support."},
        plot_data=[],
    )
    _emit("featuremap", written + [str(Path(args.out).with_suffix(".manifest.json"))])
    return 0

