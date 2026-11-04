"""MS1 feature-family discovery used by OpenSearch."""

import argparse
import bisect
import json
import math
import statistics
import sys
import tempfile
from pathlib import Path

from toppic_io import read_feature_rows

REFINE_MIN_PEARSON = 0.3
REFINE_MIN_MATCHED_PEAKS = 3
REFINE_MIN_SCAN_COUNT = 2


def refine_feature_masses(features, mzml, binary):
    """Optional RT-stacked mass calibration for every feature (salvaged from
    the parked two-pass machine — removed, in git history — because it is
    generally useful).

    The kernel's refine-shift averagine path (no sequence needed) runs per
    feature with the extraction window set to the feature's OWN RT span --
    isotope-envelope evidence stacked along the whole elution, minimum
    half-width 0.5 min (single-scan spans cannot build a coherent isotope
    chain). The fit corrects the monoisotopic mass and the isotope index;
    shape-QC failures keep their raw mass (no evidence, no change).
    """
    from _shared import native  # deferred: the no-mzml path stays kernel-free
    usable = [index for index, row in enumerate(features)
              if row["rep_charge"] > 0 and row["rt_max"] >= row["rt_min"]]
    if not usable:
        return {"skipped": True}
    proton = native.mass_constants(binary)["protonMass"]

    def window(row):
        half = max(0.5, (row["rt_max"] - row["rt_min"]) / 2)
        center = (row["rt_min"] + row["rt_max"]) / 2
        return center - half, center + half

    requests = [{
        "id": f"Q{k}",
        "sequenceMass": round(features[index]["mass"], 6),
        "modificationMass": 0.0,
        "mz": round((features[index]["mass"]
                     + features[index]["rep_charge"] * proton)
                    / features[index]["rep_charge"], 6),
        "charge": features[index]["rep_charge"],
        "rtMin": round(window(features[index])[0], 4),
        "rtMax": round(window(features[index])[1], 4),
    } for k, index in enumerate(usable)]
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8",
                                     suffix=".json", delete=False) as fh:
        json.dump(requests, fh, separators=(",", ":"))
        query_path = Path(fh.name)
    try:
        kernel = native.binary(binary)
        batch = native.run_json([kernel, "refine-shift-batch", str(Path(mzml)),
                                 "--input", str(query_path)],
                                label="tdeasenexus-toolbox refine-shift-batch")
    finally:
        query_path.unlink(missing_ok=True)
    results = {item["id"]: item
               for item in batch.get("result", {}).get("results", [])}
    refined = 0
    isotope_corrected = 0
    abs_shifts = []
    for k, index in enumerate(usable):
        item = results.get(f"Q{k}")
        if not item or not item.get("ok"):
            continue
        result = item["result"]
        pearson = result.get("pearson")
        if not (isinstance(pearson, (int, float))
                and pearson >= REFINE_MIN_PEARSON
                and int(result.get("matchedPeaks") or 0) >= REFINE_MIN_MATCHED_PEAKS
                and int(result.get("scanCount") or 0) >= REFINE_MIN_SCAN_COUNT):
            continue
        abs_shifts.append(abs(result["refinedMass"] - features[index]["mass"]))
        if int(result.get("offset") or 0) != 0:
            isotope_corrected += 1
        features[index]["mass"] = float(result["refinedMass"])
        refined += 1
    abs_shifts.sort()
    return {
        "features": len(usable),
        "refined": refined,
        "isotope_corrected": isotope_corrected,
        "median_abs_shift_mda": round(statistics.median(abs_shifts) * 1000, 2)
        if abs_shifts else None,
    }


def _ppm(mass, ppm):
    return mass * ppm / 1_000_000.0


def _extend(features, masses, first, second, spacing_min, spacing_max, ppm, maximum):
    members = [first, second]
    gaps = [features[second]["mass"] - features[first]["mass"]]
    common_start = max(features[first]["rt_min"], features[second]["rt_min"])
    common_end = min(features[first]["rt_max"], features[second]["rt_max"])
    while len(members) < maximum:
        spacing = sum(gaps) / len(gaps)
        expected = features[members[-1]]["mass"] + spacing
        lo = max(expected - _ppm(expected, ppm), features[members[-1]]["mass"] + spacing_min)
        hi = min(expected + _ppm(expected, ppm), features[members[-1]]["mass"] + spacing_max)
        candidates = []
        for index in range(bisect.bisect_left(masses, lo, members[-1] + 1), bisect.bisect_right(masses, hi)):
            item = features[index]
            start = max(common_start, item["rt_min"])
            end = min(common_end, item["rt_max"])
            if start < end:
                candidates.append((abs(item["mass"] - expected), -item["intensity"], index, start, end))
        if not candidates:
            break
        _, _, chosen, common_start, common_end = min(candidates)
        gaps.append(features[chosen]["mass"] - features[members[-1]]["mass"])
        members.append(chosen)
    return tuple(members)


def discover(features, spacing_min, spacing_max, ppm, min_members, max_members, overlap_fraction):
    """MS1-only ladder discovery over a COMMON RT intersection.

    NOT the same kernel as ingest.find_identified_serial_families despite
    the shared ladder vocabulary: MS1-only features must coelute to be one
    family, while identified serial families allow bounded per-step RT
    drift. Do not merge them.
    """
    masses = [row["mass"] for row in features]
    paths = set()
    for first, left in enumerate(features):
        start = bisect.bisect_left(masses, left["mass"] + spacing_min, first + 1)
        end = bisect.bisect_right(masses, left["mass"] + spacing_max, start)
        for second in range(start, end):
            right = features[second]
            if max(left["rt_min"], right["rt_min"]) >= min(left["rt_max"], right["rt_max"]):
                continue
            path = _extend(features, masses, first, second, spacing_min, spacing_max, ppm, max_members)
            if len(path) >= min_members:
                paths.add(path)
    maximal = []
    for path in sorted(paths, key=lambda item: (-len(item), item)):
        if not any(set(path) < set(existing) for existing in maximal):
            maximal.append(path)
    output = []
    for path in maximal:
        members = [features[index] for index in path]
        gaps = [b["mass"] - a["mass"] for a, b in zip(members, members[1:])]
        center = sum(gaps) / len(gaps)
        output.append({
            "member_count": len(members),
            "feature_ids": [item["feature_id"] for item in members],
            "masses": [round(item["mass"], 6) for item in members],
            "mass_gaps": [round(gap, 6) for gap in gaps],
            "spacing_center_da": round(center, 6),
            "spacing_spread_da": round(max(gaps) - min(gaps), 6),
            "gap_error_ppm_intact": [round((gap - center) / item["mass"] * 1e6, 3) for gap, item in zip(gaps, members[1:])],
            "common_rt": [round(max(item["rt_min"] for item in members), 6), round(min(item["rt_max"] for item in members), 6)],
            "rt_apexes": [round(item["rt_apex"], 6) for item in members],
            "intensities": [item["intensity"] for item in members],
            "ec_scores": [round(item["ec_score"], 6) for item in members],
        })
    output.sort(key=lambda item: (-item["member_count"], item["spacing_spread_da"], item["masses"][0]))
    result = []
    for family in output:
        ids = set(family["feature_ids"])
        if any(overlap_fraction > 0 and len(ids & set(old["feature_ids"])) / min(len(ids), len(old["feature_ids"])) >= overlap_fraction for old in result):
            continue
        family["family_id"] = f"FAM{len(result) + 1:04d}"
        result.append(family)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature", required=True)
    parser.add_argument("--spacing-min", type=float, required=True,
                        help="lower searched neutral-mass spacing in Da")
    parser.add_argument("--spacing-max", type=float, required=True,
                        help="upper searched neutral-mass spacing in Da")
    parser.add_argument("--ppm", type=float, default=20.0)
    parser.add_argument("--min-members", type=int, default=3)
    parser.add_argument("--max-members", type=int, default=20)
    parser.add_argument("--min-intensity", type=float, default=0.0)
    parser.add_argument("--min-ec-score", type=float, default=0.0)
    parser.add_argument("--overlap-fraction", type=float, default=0.75)
    parser.add_argument("--limit", type=int, default=20)
    parser.add_argument("--mzml",
                        help="raw mzML for optional RT-stacked mass/isotope refinement; "
                             "auto-discover the same-stem .mzML when omitted")
    parser.add_argument("--no-refine", action="store_true",
                        help="skip refinement even if mzML is available, retaining source masses")
    parser.add_argument("--tdeasenexus-toolbox", default=None)
    parser.add_argument("--out")
    args = parser.parse_args(argv)
    if not (math.isfinite(args.spacing_min) and math.isfinite(args.spacing_max)
            and 0 < args.spacing_min <= args.spacing_max):
        parser.error("spacing bounds must be finite and satisfy 0 < min <= max")
    features = [row for row in read_feature_rows(Path(args.feature)) if row["intensity"] >= args.min_intensity and row["ec_score"] >= args.min_ec_score]
    layer_refine = {"skipped": True}
    mzml = args.mzml
    if mzml is None and not args.no_refine:
        stem = Path(args.feature).with_suffix("")          # drop .feature
        if stem.name.endswith("_ms1"):
            stem = stem.with_name(stem.name[:-len("_ms1")])  # drop _ms1
        candidate = stem.with_suffix(".mzML")
        if candidate.is_file():
            mzml = str(candidate)
    if mzml and not args.no_refine:
        print(f"# RT-stacked refine: {mzml}", file=sys.stderr)
        layer_refine = refine_feature_masses(features, mzml, args.tdeasenexus_toolbox)
        features.sort(key=lambda row: (row["mass"], row["feature_id"]))
        print(f"# refined {layer_refine.get('refined', 0)}/{layer_refine.get('features', 0)} "
              f"features (median |shift| {layer_refine.get('median_abs_shift_mda')} mDa)",
              file=sys.stderr)
    elif not args.no_refine:
        print("# no sibling mzML found; running on the uncalibrated feature table "
              "(pass --mzml to calibrate)", file=sys.stderr)
    families = discover(features, args.spacing_min, args.spacing_max, args.ppm, args.min_members, args.max_members, args.overlap_fraction)
    payload = {"feature_file": args.feature, "mzml": mzml or None,
               "layer_refine": layer_refine,
               "features_scanned": len(features), "total_family_count": len(families), "family_count": len(families if args.limit == 0 else families[:args.limit]), "families": families if args.limit == 0 else families[:args.limit]}
    text = json.dumps(payload, ensure_ascii=False, indent=2) + "\n"
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(text, encoding="utf-8")
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
