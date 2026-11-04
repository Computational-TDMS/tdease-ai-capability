#!/usr/bin/env python3
"""Generate proteoform-map/data.json — all proteoforms of one accession as
packed horizontal sequence bars (TopPICR features.csv), with modifications
marked on the bars.

Lane packing ports tdp-opensearch `_pack_proteoform_rows` verbatim: rows keep
TopPICR's (first_aa asc, last_aa desc) order and non-overlapping intervals
share a lane, so short terminal forms tile onto the same row. Mod parsing
ports `_proteoform_mods`: `(SEG)[+delta]` places the mark after the segment
(counting the paren-stripped prefix); named mods (Acetyl…) are explicit.

--box "MASS|MASS|…:LABEL" (repeatable) declares the identified/adjudicated
modifications: masses matched within 0.15 Da become boxed family marks
carrying one shared number + side annotation; everything else stays a grey
ambiguous dot. Boxes are analyst input only — they assert an adjudication
this script never invents; with no --box every modification stays a grey
ambiguous dot."""
import argparse
import csv
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[6]


def parse_mods(proteoform: str) -> list[dict]:
    """Port of figures.py _proteoform_mods: (pos, mass) per [...] block."""
    core = re.sub(r"^.??\.", "", proteoform or "")
    core = re.sub(r"\..??$", "", core)
    out = []
    for m in re.finditer(r"\[(.+?)\]", core):
        pos = max(1, len(re.sub(r"[()]", "", core[:m.start()])))
        try:
            out.append({"pos": pos, "mass": float(m.group(1))})
        except ValueError:
            out.append({"pos": pos, "mass": None, "name": m.group(1)})
    return out


def pack_lanes(rows: list[dict]) -> list[int]:
    """Greedy first-fit over (first asc, last desc) — non-overlapping share a lane."""
    lane_ends: list[float] = []
    lane_of = [0] * len(rows)
    for i in sorted(range(len(rows)), key=lambda i: (rows[i]["first"], -rows[i]["last"])):
        lane = next((l for l, end in enumerate(lane_ends) if rows[i]["first"] > end + 1), None)
        if lane is None:
            lane_ends.append(0.0)
            lane = len(lane_ends) - 1
        lane_ends[lane] = rows[i]["last"]
        lane_of[i] = lane
    return lane_of


def main() -> None:
    here = Path(__file__).parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--features", type=Path, required=True,
                    help="TopPICR features.csv to read")
    ap.add_argument("--accession", required=True)
    ap.add_argument("--box", action="append", default=[], metavar="MASS|MASS:LABEL",
                    help="identified mods: boxed + side-annotated (repeatable)")
    ap.add_argument("--out", type=Path, default=here / "data.json")
    args = ap.parse_args()
    boxes = args.box
    try:
        source_display = str(args.features.resolve().relative_to(REPO))
    except ValueError:
        source_display = str(args.features)

    families = []
    for spec in boxes:
        masses_s, _, label = spec.partition(":")
        families.append({
            "n": len(families) + 1,
            "masses": [float(x) for x in masses_s.split("|")],
            "label": label,
        })

    rows = []
    with args.features.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row.get("uniprot_acc") == args.accession:
                rows.append(row)
    if not rows:
        raise SystemExit(f"accession {args.accession} not in {args.features}")

    prot_length = max(float(r.get("prot_length") or 0) for r in rows)

    def first_of(r):
        return float(r["first_aa"]) if r.get("first_aa") else 1.0

    def last_of(r):
        # TopPICR leaves first/last empty on full-length forms — that means the
        # whole sequence, not a 2-aa stub.
        return float(r["last_aa"]) if r.get("last_aa") else prot_length

    lane_of = pack_lanes([{**r, "first": first_of(r), "last": last_of(r)} for r in rows])
    n_lanes = max(lane_of) + 1

    payload_rows, fam_bars = [], {f["n"]: 0 for f in families}
    for i, r in enumerate(rows):
        mods = []
        for m in parse_mods(r.get("proteoform", "")):
            if m["mass"] is not None:
                fam = next((f for f in families
                            if any(abs(fM - m["mass"]) <= 0.15 for fM in f["masses"])), None)
            else:
                fam = None
            entry = {"pos": m["pos"], "mass": m.get("mass")}
            if m.get("name"):
                entry["name"] = m["name"]
            if fam:
                entry["box"] = fam["n"]
                fam_bars[fam["n"]] += 1
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
    payload = {
        "data": {
            "accession": args.accession,
            "proteinLength": prot_length,
            "nLanes": n_lanes,
            "rows": payload_rows,
            "boxes": box_notes,
        },
        "meta": {
            "title": f"{args.accession} proteoform map — {len(rows)} forms, {n_lanes} lanes",
            "subtitle": ("bars = sequence intervals (first–last aa), colour = PrSM count (batlow, darker = more) · "
                         "grey dot = unresolved shift · boxed number = identified family (side notes right) · "
                         "non-overlapping bars share a lane"),
            "source": f"{source_display} · accession {args.accession}",
            "transform": "lane packing = _pack_proteoform_rows port; mods = _proteoform_mods port; box families matched within 0.15 Da",
            "colour_scale": "crameri-batlow (sequential, sqrt PrSM count)",
            "n": f"{len(rows)} proteoforms · {sum(fam_bars.values())} boxed mods",
        },
    }
    args.out.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    print(f"{args.accession}: {len(rows)} forms -> {n_lanes} lanes; boxes: "
          + "; ".join(f"#{b['n']} x{b['nBars']}" for b in box_notes))
    print(f"-> {args.out}")


if __name__ == "__main__":
    main()
