"""Proteoform map panel (proteoform-map template, packed lanes + boxes)."""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from _shared.manifest import write_manifest
from figures_common import _load_json, _render, _emit


def _pack_proteoform_rows(rows):
    """Pack non-overlapping sequence intervals onto shared lanes (greedy
    first-fit; preserves TopPICR's first-asc/last-desc ordering)."""
    layer_ends = []
    packed = []
    for row in rows:
        lane = next(
            (index for index, end in enumerate(layer_ends)
             if row["first_aa"] > end + 1),
            None,
        )
        if lane is None:
            layer_ends.append(row["last_aa"])
            lane = len(layer_ends) - 1
        else:
            layer_ends[lane] = max(layer_ends[lane], row["last_aa"])
        packed.append(lane)
    return packed


def _proteoform_mods(proteoform: str, labels_by_mass=None):
    """Return (label, position, explicit, raw mass shift) for ``[...]`` mods."""
    if not proteoform:
        return []
    core = re.sub(r"^.??\.", "", proteoform)
    core = re.sub(r"\..??$", "", core)
    out = []
    for match in re.finditer(r"\[(.+?)\]", core):
        clean_before = re.sub(r"[()]", "", core[:match.start()])
        position = max(1, len(clean_before))
        label = match.group(1)
        raw_mass = None
        try:
            raw_mass = float(label)
            hit = min(labels_by_mass.items(), key=lambda item: abs(item[0] - raw_mass)) if labels_by_mass else None
            if hit and abs(hit[0] - raw_mass) <= 0.1:
                label, explicit = hit[1], True
            else:
                label, explicit = f"{raw_mass:+.3f} Da", False
        except ValueError:
            if label == "Acetyl" and position == 1:
                label = "N-Acetyl"
            explicit = True
        out.append((label, position, explicit, raw_mass))
    return out


def proteoform_payload(features_path: str, accession: str,
                       annotations_path: str | None = None) -> dict:
    """Build the proteoform-map contract (figure + report)."""
    rows = []
    with open(features_path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("uniprot_acc") == accession:
                rows.append(row)
    if not rows:
        raise SystemExit(f"accession {accession} not found")
    rows.sort(key=lambda row: (float(row.get("first_aa") or 0), -float(row.get("last_aa") or 0)))
    protein_length = max(float(row.get("prot_length") or 0) for row in rows)

    def first_of(r):
        return float(r["first_aa"]) if r.get("first_aa") else 1.0

    def last_of(r):
        # TopPICR leaves first/last empty on full-length forms — the whole
        # sequence, not a 2-aa stub.
        return float(r["last_aa"]) if r.get("last_aa") else protein_length

    lane_of = _pack_proteoform_rows(
        [{**r, "first_aa": first_of(r), "last_aa": last_of(r)} for r in rows])
    n_lanes = max(lane_of) + 1

    labels_by_mass = {}
    if annotations_path:
        payload, _ = _load_json(annotations_path)
        for item in payload.get("annotations", []):
            try:
                role = item.get("role")
                basis = str(item.get("basis") or "")
                explicit = role in {"manual", "deep"} or (
                    role == "automatic"
                    and basis.startswith("automatic approved UniMod mass attribution")
                )
                if explicit:
                    labels_by_mass[float(item["mass"])] = str(item["label"])
            except (KeyError, TypeError, ValueError):
                continue

    # Explicit families: one numbered box per distinct label, masses matched
    # within 0.1 Da (the plotly version's explicit-attribution contract).
    families = [{"n": i + 1, "masses": [mass], "label": label}
                for i, (mass, label) in enumerate(sorted(labels_by_mass.items()))]

    payload_rows, fam_bars = [], {f["n"]: 0 for f in families}
    for i, r in enumerate(rows):
        mods = []
        for label, pos, explicit, raw_mass in _proteoform_mods(r.get("proteoform", ""), labels_by_mass):
            entry = {"pos": pos}
            if raw_mass is not None:
                entry["mass"] = raw_mass
                fam = next((f for f in families
                            if any(abs(fm - raw_mass) <= 0.1 for fm in f["masses"])), None)
                if fam:
                    entry["box"] = fam["n"]
                    fam_bars[fam["n"]] += 1
            if not explicit and raw_mass is None:
                entry["name"] = label
            mods.append(entry)
        payload_rows.append({
            "first": first_of(r),
            "last": last_of(r),
            "count": float(r.get("spectral_count") or 0),
            "lane": lane_of[i],
            "proforma": r.get("proteoform", ""),
            "mods": mods,
        })

    box_notes = [{**f, "nBars": fam_bars[f["n"]]} for f in families if fam_bars[f["n"]]]
    return {"data": {"accession": accession, "proteinLength": protein_length,
                     "nLanes": n_lanes, "rows": payload_rows, "boxes": box_notes},
            "meta": {"title": f"{accession} proteoform map — {len(rows)} forms, {n_lanes} lanes",
                     "subtitle": ("bars = sequence intervals (first–last aa), colour = PrSM count (batlow, darker = more) · "
                                  "grey dot = unresolved shift · boxed number = identified family (side notes right) · "
                                  "non-overlapping bars share a lane"),
                     "source": f"{Path(features_path).name} · accession {accession}",
                     "colour_scale": "crameri-batlow (sequential, sqrt PrSM count)",
                     "n": f"{len(rows)} proteoforms"}}


def proteoform_main(argv=None):
    ap = argparse.ArgumentParser(description="Proteoform sequence intervals + PTM boxes (chart-forge).")
    ap.add_argument("--features", required=True, help="TopPICR features.csv")
    ap.add_argument("--accession", required=True)
    ap.add_argument("--out", "--plot", dest="out", required=True, help="output HTML path")
    ap.add_argument("--annotations", help="OpenSearch automatic/manual annotation JSON")
    args = ap.parse_args(argv)

    payload = proteoform_payload(args.features, args.accession, args.annotations)
    written = _render("proteoform-map", payload, args.out)
    write_manifest(
        args.out, figure="proteoform_ptm_map", script=Path(__file__).name,
        params={"accession": args.accession},
        inputs=[args.features] + ([args.annotations] if args.annotations else []),
        caption={"meaning": "Sequence intervals and modification overlays; separate from mass-shift census."},
        plot_data=[],
    )
    _emit("proteoform", written + [str(Path(args.out).with_suffix(".manifest.json"))])
    return 0

