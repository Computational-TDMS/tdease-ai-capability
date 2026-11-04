"""
Ingest family discovery: the numeric/isotope helpers and the two ladder
kernels. Grouping uses kernel-provided constants; no spectrum or artifact IO.

Two deliberately different ladder kernels live here (do not merge them):
- find_identified_serial_families — same-protein, PrSM-backed, bounded
  per-step RT drift;
- find_feature_families — MS1-only, common-RT intersection, delegates the
  path search to family_discovery.discover.
"""

import bisect
import statistics
from collections import Counter, defaultdict

from toppic_io import c13_gap, prsm_passes
from family_discovery import discover as discover_feature_families


def ppm_tolerance(mass: float, ppm: float) -> float:
    return mass * ppm / 1_000_000.0

def c13_equivalent(left, right, ppm_tol=5.0, n_max=1):
    delta = left - right
    isotope = round(delta / c13_gap())
    if abs(isotope) > n_max:
        return False
    residual = delta - isotope * c13_gap()
    return abs(residual) / max(left, right, 1.0) * 1e6 <= ppm_tol

def rt_gap(left, right):
    if left["rt_end"] < right["rt_start"]:
        return right["rt_start"] - left["rt_end"]
    if right["rt_end"] < left["rt_start"]:
        return left["rt_start"] - right["rt_end"]
    return 0.0

def normalized_mass(masses):
    reference = statistics.median(masses)
    corrected = [mass - round((mass - reference) / c13_gap()) * c13_gap() for mass in masses]
    return statistics.median(corrected)

def shift_equivalent(left, right, total_mass, ppm_tol=10.0, n_max_c13=1):
    delta = left - right
    tolerance = total_mass * ppm_tol / 1e6
    return any(abs(delta - isotope * c13_gap()) <= tolerance
               for isotope in range(-n_max_c13, n_max_c13 + 1))

def format_counts(values):
    counts = Counter(values)
    return "; ".join(f"{key}:{counts[key]}" for key in sorted(counts))


SERIAL_MIN_MEMBERS = 3
SERIAL_MAX_MEMBERS = 12
SERIAL_SPACING_MIN_DA = 10.0
SERIAL_SPACING_MAX_DA = 600.0
SERIAL_RT_BACKTRACK = 0.5
SERIAL_MAX_RT_DRIFT_DA = 12.0


def find_identified_serial_families(
    ms1, prsms, qvalue_max=0.01, evalue_fallback=1e-5,
    ppm_tol=10.0, n_max_c13=1,
    max_families_per_protein=3,
):
    """Same-protein, PrSM-backed serial families with bounded RT drift.

    NOT the same kernel as family_discovery.discover despite the shared
    ladder vocabulary: that one requires a common RT intersection (MS1-only
    features must coelute to be one family); this one allows per-step RT
    drift because identified proteoforms can shift apex with each
    modification. Do not merge them.

    This is a separate evidence channel from the MS1-only family search.
    Identified proteoforms need not share one RT intersection: successive
    modifications can move the chromatographic apex. Qualifying PrSMs and a
    common protein replace that MS1-only constraint. The result records a
    measured family and never assigns its chemistry.
    """
    support = defaultdict(lambda: defaultdict(list))
    for row in prsms:
        feature_id = row.get("Feature ID") or ""
        protein = row.get("Protein accession") or ""
        if feature_id not in ms1 or not protein:
            continue
        if prsm_passes(row, qvalue_max, evalue_fallback):
            support[protein][feature_id].append(row)

    raw_paths = set()
    for protein, by_feature in support.items():
        features = [
            dict(ms1[feature_id], feature_id=feature_id)
            for feature_id in by_feature
        ]
        features.sort(key=lambda row: (row["mass"], str(row["feature_id"])))
        if len(features) < SERIAL_MIN_MEMBERS:
            continue
        masses = [row["mass"] for row in features]

        def normalize_gap(raw_gap, center):
            isotope = round((center - raw_gap) / c13_gap())
            if abs(isotope) > n_max_c13:
                return None
            corrected = raw_gap + isotope * c13_gap()
            return corrected, isotope

        buckets = defaultdict(set)
        for first, left in enumerate(features):
            start = bisect.bisect_left(
                masses, left["mass"] + SERIAL_SPACING_MIN_DA, first + 1)
            end = bisect.bisect_right(
                masses, left["mass"] + SERIAL_SPACING_MAX_DA, start)
            for second in range(start, end):
                right = features[second]
                if abs(right["rt_apex"] - left["rt_apex"]) > SERIAL_MAX_RT_DRIFT_DA:
                    continue
                raw_gap = right["mass"] - left["mass"]
                for isotope in range(-n_max_c13, n_max_c13 + 1):
                    corrected = raw_gap + isotope * c13_gap()
                    if SERIAL_SPACING_MIN_DA <= corrected <= SERIAL_SPACING_MAX_DA:
                        buckets[round(corrected * 2) / 2].add(first)

        # This is a discovery kernel. Do not discard a chemically relevant
        # window because a noisy pair-count ranking pushed it past a display
        # cap; cap the returned families after paths have been built.
        windows = sorted(
            ((len(firsts), bucket) for bucket, firsts in buckets.items()
             if len(firsts) >= SERIAL_MIN_MEMBERS - 1),
            key=lambda item: (-item[0], item[1]),
        )
        for _support_count, bucket in windows:
            for first, left in enumerate(features):
                search_width = n_max_c13 * c13_gap() + 0.3
                start = bisect.bisect_left(
                    masses, left["mass"] + bucket - search_width, first + 1)
                end = bisect.bisect_right(
                    masses, left["mass"] + bucket + search_width, start)
                for second in range(start, end):
                    right = features[second]
                    normalized = normalize_gap(right["mass"] - left["mass"], bucket)
                    if normalized is None or abs(normalized[0] - bucket) > 0.3:
                        continue
                    direction = 1 if right["rt_apex"] >= left["rt_apex"] else -1
                    members = [first, second]
                    gaps = [normalized[0]]
                    while len(members) < SERIAL_MAX_MEMBERS:
                        spacing = sum(gaps) / len(gaps)
                        previous = features[members[-1]]
                        tolerance = max(ppm_tolerance(previous["mass"] + spacing, ppm_tol), 0.1)
                        lo = previous["mass"] + spacing - n_max_c13 * c13_gap() - tolerance
                        hi = previous["mass"] + spacing + n_max_c13 * c13_gap() + tolerance
                        candidate_start = bisect.bisect_left(
                            masses, lo, members[-1] + 1)
                        candidate_end = bisect.bisect_right(
                            masses, hi, candidate_start)
                        candidates = []
                        for index in range(candidate_start, candidate_end):
                            candidate = features[index]
                            rt_step = direction * (
                                candidate["rt_apex"] - previous["rt_apex"])
                            if rt_step < -SERIAL_RT_BACKTRACK:
                                continue
                            if abs(candidate["rt_apex"] - left["rt_apex"]) > SERIAL_MAX_RT_DRIFT_DA:
                                continue
                            normalized = normalize_gap(
                                candidate["mass"] - previous["mass"], spacing)
                            if normalized is None or abs(normalized[0] - spacing) > tolerance:
                                continue
                            candidates.append((
                                abs(normalized[0] - spacing),
                                -candidate["intensity"], index, normalized[0],
                            ))
                        if not candidates:
                            break
                        chosen_candidate = min(candidates)
                        chosen = chosen_candidate[2]
                        gaps.append(chosen_candidate[3])
                        members.append(chosen)
                    if len(members) >= SERIAL_MIN_MEMBERS:
                        raw_paths.add((
                            protein,
                            tuple(features[index]["feature_id"] for index in members),
                        ))

    maximal_paths = []
    for path in sorted(raw_paths, key=lambda item: (-len(item[1]), item)):
        protein, feature_ids = path
        ids = set(feature_ids)
        if any(protein == existing_protein and ids < set(existing_ids)
               for existing_protein, existing_ids in maximal_paths):
            continue
        maximal_paths.append(path)

    output = []
    for protein, feature_ids in maximal_paths:
        members = [ms1[feature_id] for feature_id in feature_ids]
        raw_gaps = [right["mass"] - left["mass"]
                    for left, right in zip(members, members[1:])]
        centers = {
            gap + isotope * c13_gap()
            for gap in raw_gaps
            for isotope in range(-n_max_c13, n_max_c13 + 1)
            if SERIAL_SPACING_MIN_DA <= gap + isotope * c13_gap() <= SERIAL_SPACING_MAX_DA
        }
        scored = []
        for candidate_center in centers:
            candidate = [normalize_gap(gap, candidate_center) for gap in raw_gaps]
            if any(item is None for item in candidate):
                continue
            corrected = [item[0] for item in candidate]
            fitted = sum(corrected) / len(corrected)
            scored.append((
                max(corrected) - min(corrected),
                sum(abs(gap - fitted) for gap in corrected),
                abs(fitted),
                fitted,
                candidate,
            ))
        if not scored:
            continue
        _spread, _residual, _magnitude, center, normalized = min(scored)
        corrected_gaps = [item[0] for item in normalized]
        isotope_offsets = [item[1] for item in normalized]
        member_records = []
        for feature_id, member in zip(feature_ids, members):
            rows = support[protein][feature_id]
            member_records.append({
                "featureId": str(feature_id),
                "mass": round(member["mass"], 6),
                "rtApex": round(member["rt_apex"], 6),
                "prsmCount": len({row.get("Prsm ID") for row in rows
                                  if row.get("Prsm ID")}),
                "prsmIds": sorted({row.get("Prsm ID") for row in rows
                                   if row.get("Prsm ID")}),
                "proteoformIds": sorted({row.get("Proteoform ID") for row in rows
                                         if row.get("Proteoform ID")}),
            })
        output.append({
            "evidence_type": "identified_same_protein_serial",
            "chemistry_status": "unassigned",
            "member_count": len(members),
            "feature_ids": list(feature_ids),
            "protein": protein,
            "masses": [round(member["mass"], 6) for member in members],
            "raw_mass_gaps": [round(gap, 6) for gap in raw_gaps],
            "mass_gaps": [round(gap, 6) for gap in corrected_gaps],
            "isotope_offsets": isotope_offsets,
            "spacing_center_da": round(center, 6),
            "spacing_spread_da": round(max(corrected_gaps) - min(corrected_gaps), 6),
            "gap_error_ppm_intact": [
                round((gap - center) / member["mass"] * 1_000_000.0, 3)
                for gap, member in zip(corrected_gaps, members[1:])
            ],
            "rt_apexes": [round(member["rt_apex"], 6) for member in members],
            "rt_drift": round(members[-1]["rt_apex"] - members[0]["rt_apex"], 6),
            "intensity_total": sum(member["intensity"] for member in members),
            "prsm_count": sum(record["prsmCount"] for record in member_records),
            "member_records": member_records,
        })
    output.sort(key=lambda family: (
        -family["member_count"], family["spacing_spread_da"],
        -family["prsm_count"], -family["intensity_total"], family["masses"][0],
    ))
    returned = []
    protein_counts = Counter()
    for family in output:
        protein = family["protein"]
        if protein_counts[protein] >= max_families_per_protein:
            continue
        returned.append(family)
        protein_counts[protein] += 1
    for family_index, family in enumerate(returned, 1):
        family["family_id"] = f"SERIAL{family_index:04d}"
    return returned


def find_feature_families(ms1, identified, ppm_tol=10.0, min_members=3,
                          max_members=12, spacing_min=10.0, spacing_max=600.0,
                          max_families=24):
    """Repeated-spacing families among MS1-only features (discovery scaffold).

    Spacing windows come from the data itself -- a histogram of pairwise gaps
    between RT-intersecting MS1-only features -- so any ladder chemistry
    (multi-acetyl, multiphospho, glycan series) surfaces without telling the
    screen what to look for. Families are Agent reference only: MS1-envelope
    evidence with no assignment weight, features never merge. Uses the
    feature-family kernel's path search verbatim.
    """
    features = [
        {"feature_id": fid, "mass": node["mass"], "intensity": node["intensity"],
         "rt_min": node["rt_start"], "rt_max": node["rt_end"],
         "rt_apex": node["rt_apex"], "ec_score": node.get("ec_score", 0.0)}
        for fid, node in ms1.items()
        if fid not in identified and node["mass"] > 0
    ]
    features.sort(key=lambda row: (row["mass"], str(row["feature_id"])))
    if len(features) < min_members:
        return [], {
            "spacingWindowsInHistogram": 0,
            "spacingWindowsSearched": 0,
            "familiesFoundBeforeCap": 0,
            "familiesKept": 0,
            "familyCap": max_families,
            "note": f"fewer than min_members={min_members} MS1-only features",
        }
    masses = [row["mass"] for row in features]

    def rt_overlap(left, right):
        return max(left["rt_min"], right["rt_min"]) < min(left["rt_max"], right["rt_max"])

    buckets = defaultdict(set)
    for i, left in enumerate(features):
        start = bisect.bisect_left(masses, left["mass"] + spacing_min, i + 1)
        end = bisect.bisect_right(masses, left["mass"] + spacing_max, start)
        for k in range(start, end):
            if rt_overlap(left, features[k]):
                buckets[round((features[k]["mass"] - left["mass"]) * 2) / 2].add(i)
    ranked_windows = sorted(((len(firsts), bucket) for bucket, firsts in buckets.items()
                             if len(firsts) >= 2),
                            key=lambda item: (-item[0], item[1]))
    # The histogram shortlist and the family cap are search budgets, not the
    # whole chemistry space: return what was searched so the parts never
    # examined stay visible in the bundle instead of silently absent.
    search_coverage = {
        "spacingWindowsInHistogram": len(ranked_windows),
        "spacingWindowsSearched": min(len(ranked_windows), 12),
        "familiesFoundBeforeCap": None,
        "familiesKept": None,
        "familyCap": max_families,
    }
    windows = ranked_windows[:12]
    by_id = {str(row["feature_id"]): row for row in features}
    families = []
    seen = set()
    for _support, bucket in windows:
        for fam in discover_feature_families(features, bucket - 0.3, bucket + 0.3,
                                     ppm_tol, min_members, max_members, 0.0):
            key = tuple(fam["feature_ids"])
            if key in seen:
                continue
            seen.add(key)
            fam["intensity_total"] = sum(
                by_id[str(fid)]["intensity"] for fid in fam["feature_ids"])
            families.append(fam)
    families.sort(key=lambda fam: (-fam["member_count"], -fam["intensity_total"]))
    search_coverage["familiesFoundBeforeCap"] = len(families)
    kept = families[:max_families]
    search_coverage["familiesKept"] = len(kept)
    return kept, search_coverage
