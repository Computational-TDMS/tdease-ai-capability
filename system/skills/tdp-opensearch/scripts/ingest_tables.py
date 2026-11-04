"""
Ingest table construction: MS1 feature merging, per-sample typed rows, and
the checklist projection used by the bundle. Reads TopPIC tables; callers
pass open paths. No kernel subprocess here.
"""

import hashlib
import json
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from levels import classify_row
from toppic_io import (
    c13_gap,
    number,
    parse_mods,
    parse_proforma,
    prsm_passes,
    read_ms1,
    read_toppic,
)
from ingest_families import (
    c13_equivalent,
    find_feature_families,
    find_identified_serial_families,
    format_counts,
    normalized_mass,
    rt_gap,
    shift_equivalent,
)


def merge_source_features(nodes, ppm_tol=5.0, rt_gap_limit=0.2, n_max_c13=1):
    """Merge only identification-compatible, isotope-equivalent RT segments."""
    parent = list(range(len(nodes)))
    minimum = [node["mass"] for node in nodes]
    maximum = [node["mass"] for node in nodes]

    def find(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    pairs = []
    for left in range(len(nodes)):
        for right in range(left + 1, len(nodes)):
            if not nodes[left]["proteins"].intersection(nodes[right]["proteins"]):
                continue
            if not c13_equivalent(nodes[left]["mass"], nodes[right]["mass"], ppm_tol, n_max_c13):
                continue
            gap = rt_gap(nodes[left], nodes[right])
            if gap <= rt_gap_limit:
                pairs.append((gap, abs(nodes[left]["mass"] - nodes[right]["mass"]), left, right))
    pairs.sort()

    for _gap, _mass_gap, left, right in pairs:
        root_left, root_right = find(left), find(right)
        if root_left == root_right:
            continue
        new_min = min(minimum[root_left], minimum[root_right])
        new_max = max(maximum[root_left], maximum[root_right])
        if new_max - new_min > c13_gap() + 0.1:
            continue
        parent[root_right] = root_left
        minimum[root_left], maximum[root_left] = new_min, new_max

    groups = defaultdict(list)
    for index, node in enumerate(nodes):
        groups[find(index)].append(node)
    return list(groups.values())

def cluster_modification_masses(records, ppm_tol=10.0):
    """Assign lookup-only mass clusters without single-linkage chaining.

    The anchor is the majority mass (median), never the lowest member. This
    prevents one low-quality observation from seeding a cluster and shifting
    the dominant monoisotopic family by one C13 interval. Members are still
    isotope-equivalence adjusted around that observed anchor."""
    remaining = set(range(len(records)))
    clusters = []
    while remaining:
        pool = sorted(remaining)
        median_mass = statistics.median(
            [records[index]["shift_mass"] for index in pool])
        # The anchor must be an actual record: the median can fall between
        # two distant masses and match nothing.
        anchor = min(pool, key=lambda index: abs(
            records[index]["shift_mass"] - median_mass))
        center = records[anchor]["shift_mass"]
        anchor_mass = records[anchor]["proteoform_mass"]
        members = sorted(
            index for index in remaining
            if shift_equivalent(
                records[index]["shift_mass"], center,
                max(records[index]["proteoform_mass"], anchor_mass), ppm_tol,
            )
        )
        adjusted = [
            records[index]["shift_mass"]
            - round((records[index]["shift_mass"] - center) / c13_gap()) * c13_gap()
            for index in members
        ]
        center = statistics.median(adjusted)
        refined = sorted(
            index for index in remaining
            if shift_equivalent(
                records[index]["shift_mass"], center,
                records[index]["proteoform_mass"], ppm_tol,
            )
        )
        # A refined center between two distant masses can match nothing;
        # the first-pass membership is still a valid cluster.
        members = refined or members
        identity = "|".join(sorted(records[index]["observation_id"] for index in members))
        shift_id = "M" + hashlib.sha1(identity.encode("utf-8")).hexdigest()[:12]
        for index in members:
            records[index]["mass_shift_id"] = shift_id
            records[index]["shift_mass_center"] = f"{center:.4f}"
        clusters.append((shift_id, center, members))
        remaining.difference_update(members)
    return clusters


def representative_proteoform_rows(prsm_rows):
    """Derive the representative-proteoform view from a full PrSM table.

    TopPIC writes the one-row-per-proteoform table as its own file; runs that
    do not ship it still carry every needed column in the PrSM table, so the
    representative is the best-E-value row per Proteoform ID (ties broken by
    Prsm ID for determinism). Ingest never reads third-party identification
    tables (e.g. Informed Proteomics IcTarget): the columns are unrelated.
    """
    best = {}
    for row in prsm_rows:
        proteoform_id = row.get("Proteoform ID") or ""
        if not proteoform_id:
            continue
        rank = (number(row.get("E-value"), 1.0), row.get("Prsm ID") or "")
        current = best.get(proteoform_id)
        if current is None or rank < (number(current.get("E-value"), 1.0),
                                      current.get("Prsm ID") or ""):
            best[proteoform_id] = row
    return [best[proteoform_id] for proteoform_id in sorted(best)]


def build_single_sample(ms1_path, proteoform_path, prsm_path,
                        qvalue_max=0.01, evalue_fallback=1e-5,
                        ppm_tol=10.0, rt_gap_limit=0.2, n_max_c13=1):
    ms1 = read_ms1(ms1_path)
    prsms = read_toppic(prsm_path)
    if proteoform_path is None:
        proteoforms = representative_proteoform_rows(prsms)
    else:
        proteoforms = read_toppic(proteoform_path)
    # Join losses are part of "what was not examined": proteoform rows whose
    # Feature ID has no MS1 feature would otherwise vanish without a trace
    # and silently deflate every downstream count.
    join_summary = {
        "proteoform_rows_total": len(proteoforms),
        "proteoform_rows_dropped_no_ms1_feature": sum(
            1 for row in proteoforms if row.get("Feature ID") not in ms1),
    }
    proteoforms = [row for row in proteoforms if row.get("Feature ID") in ms1]

    proteoforms_by_feature = defaultdict(list)
    for row in proteoforms:
        proteoforms_by_feature[row["Feature ID"]].append(row)
    prsms_by_proteoform = defaultdict(list)
    for row in prsms:
        key = (
            row.get("Feature ID") or "",
            row.get("Proteoform ID") or "",
            row.get("Protein accession") or "",
        )
        prsms_by_proteoform[key].append(row)

    def qualifying_prsms(proteoform_key):
        """PrSM evidence gates which features are identified."""
        return [row for row in prsms_by_proteoform.get(proteoform_key, [])
                if prsm_passes(row, qvalue_max, evalue_fallback)]

    proteoforms = [
        row for row in proteoforms
        if qualifying_prsms((
            row.get("Feature ID") or "", row.get("Proteoform ID") or "",
            row.get("Protein accession") or "",
        ))
    ]

    # The representative-proteoform table names ONE feature per proteoform
    # (the best-E-value one), but the full PrSM table can carry qualifying
    # PrSMs on OTHER features of the same proteoform. Those orphan features
    # are real identified carriers (a 2x-stronger MS1 feature with Gold-zone
    # PrSMs can vanish this way); rescue them by cloning the proteoform row
    # relabeled to the PrSM-side feature. Feature -> Proteoform -> Prsm
    # links stay the source of truth throughout.
    pf_by_id = {}
    for row in proteoforms:
        if row.get("Proteoform ID"):
            pf_by_id.setdefault(row["Proteoform ID"], row)
    for row in prsms:
        fid = row.get("Feature ID") or ""
        if (fid in ms1 and fid not in proteoforms_by_feature
                and prsm_passes(row, qvalue_max, evalue_fallback)):
            template = pf_by_id.get(row.get("Proteoform ID") or "")
            if template is None:
                continue
            clone = dict(template)
            clone["Feature ID"] = fid
            proteoforms_by_feature[fid].append(clone)
            proteoforms.append(clone)

    # Identified carriers: the PrSM-qualified feature -> protein set. Ladder
    # detection anchors MS1-only rungs to these so each rung's shift is an
    # absolute Nx-delta_m composition rather than an unanchored spacing.
    identified = defaultdict(set)
    for row in proteoforms:
        identified[row.get("Feature ID") or ""].add(row.get("Protein accession") or "")

    nodes = []
    for source_id, rows in proteoforms_by_feature.items():
        source = dict(ms1[source_id])
        source["proteoform_rows"] = rows
        source["proteins"] = {row.get("Protein accession") or "" for row in rows}
        nodes.append(source)

    groups = merge_source_features(nodes, ppm_tol, rt_gap_limit, n_max_c13)
    feature_rows = []
    modification_rows = []
    for members in groups:
        source_ids = sorted({member["source_feature_id"] for member in members}, key=lambda value: int(value))
        feature_id = source_ids[0]
        pf_rows = [row for member in members for row in member["proteoform_rows"]]
        smith_by_row = {id(row): classify_row(row) for row in pf_rows}
        smith_levels = [smith_by_row[id(row)]["level"] for row in pf_rows]
        proteoform_keys = sorted({
            (row.get("Feature ID") or "", row.get("Proteoform ID") or "",
             row.get("Protein accession") or "")
            for row in pf_rows if row.get("Proteoform ID")
        })
        proteoform_refs = ["|".join(key) for key in proteoform_keys]
        linked_prsms = {
            row.get("Prsm ID") or "": row
            for key in proteoform_keys
            for row in prsms_by_proteoform.get(key, [])
            if row.get("Prsm ID")
        }
        if not linked_prsms:
            # Accession formatting can differ between the proteoform and PrSM
            # tables; drop ONLY the accession from the key before giving up.
            # Never fall back to the bare Feature ID: that counts PrSMs from
            # unrelated proteoform interpretations as support (the composite
            # key exists precisely to prevent that inflation).
            relaxed = defaultdict(list)
            for key, grouped in prsms_by_proteoform.items():
                relaxed[(key[0], key[1])].extend(grouped)
            linked_prsms = {
                row.get("Prsm ID") or "": row
                for key in proteoform_keys
                for row in relaxed.get((key[0], key[1]), [])
                if row.get("Prsm ID")
            }
        masses = [member["mass"] for member in members]
        feature_mass = normalized_mass(masses)
        segments = sorted((member["rt_start"], member["rt_end"]) for member in members)
        feature_rows.append({
            "feature_id": feature_id,
            "source_feature_ids": "; ".join(source_ids),
            "mass": f"{feature_mass:.4f}",
            "rt_start": f"{min(start for start, _ in segments):.4f}",
            "rt_apex": f"{statistics.median(member['rt_apex'] for member in members):.4f}",
            "rt_end": f"{max(end for _, end in segments):.4f}",
            "rt_segments": "; ".join(f"{start:.4f}-{end:.4f}" for start, end in segments),
            "intensity": f"{sum(member['intensity'] for member in members):.6f}",
            "proteoform_refs": "; ".join(proteoform_refs),
            "prsm_count": len(linked_prsms),
            "protein_accessions": "; ".join(sorted({row.get("Protein accession") or "" for row in pf_rows})),
            "smith_level_counts": format_counts(smith_levels),
            "smith_best_level": min(smith_levels),
            "smith_worst_level": max(smith_levels),
            "best_evalue": min(number(row.get("E-value"), 1.0) for row in pf_rows),
            "source_ms1": str(Path(ms1_path).resolve()),
            "source_proteoforms": str(Path(proteoform_path or prsm_path).resolve()),
            "source_prsms": str(Path(prsm_path).resolve()),
        })

        for pf_row in pf_rows:
            proteoform_id = pf_row.get("Proteoform ID") or ""
            smith = smith_by_row[id(pf_row)]
            proteoform_key = (
                pf_row.get("Feature ID") or "", proteoform_id,
                pf_row.get("Protein accession") or "",
            )
            proteoform_ref = "|".join(proteoform_key)
            all_support = qualifying_prsms(proteoform_key)
            sequence = pf_row.get("Database protein sequence") or ""
            fixed_forms, _ = parse_proforma(pf_row.get("Proteoform") or "", sequence)
            # TopPIC writes an internal sequence gap as "(RESIDUES)[+/-mass]"
            # in the proforma. A negative gap means the observed proteoform is
            # LIGHTER than any database segment can explain -- the classic
            # signature of a MISSING database entry (the true carrier is not
            # in the FASTA, so TopPIC forced a homolog with a hole).
            gap = re.search(r"\(([A-Z]+)\)\[([+-][\d.]+)\]",
                            pf_row.get("Proteoform") or "")
            seq_gap = (f"gap{gap.group(2)}[{gap.group(1)}]"
                       if gap and gap.group(2).startswith("-") else "")
            for index, (mass, start, end, is_nterm) in enumerate(
                parse_mods(pf_row.get("Unexpected modifications")), 1
            ):
                support = [
                    row for row in all_support
                    if any(shift_equivalent(observed_mass, mass, feature_mass, ppm_tol)
                           for observed_mass, *_ in parse_mods(row.get("Unexpected modifications")))
                ]
                if not support:
                    continue
                identity = f"{feature_id}|{proteoform_id}|{index}|{mass:.6f}|{start}-{end}"
                modification_rows.append({
                    "observation_id": "S" + hashlib.sha1(identity.encode("utf-8")).hexdigest()[:12],
                    "mass_shift_id": "",
                    "shift_mass_center": "",
                    "feature_id": feature_id,
                    "protein": pf_row.get("Protein accession") or "",
                    "proteoform_ref": proteoform_ref,
                    "smith_level": smith["level"],
                    "smith_subtype": smith["subtype"],
                    "smith_ambiguities": "; ".join(smith["ambiguities"]),
                    "shift_mass": mass,
                    "raw_shift_mass": mass,
                    "fixed_forms_excluded": fixed_forms,
                    "localizations": f"{start}-{end}",
                    "localization_residues": sequence[max(0, start - 1):end],
                    "sequence": sequence,
                    "is_nterm": is_nterm,
                    "proteoform_mass": feature_mass,
                    "rt_start": min(start for start, _ in segments),
                    "rt_apex": statistics.median(member["rt_apex"] for member in members),
                    "rt_end": max(end for _, end in segments),
                    "ppm_tolerance": ppm_tol,
                    "supporting_prsms": len({row.get("Prsm ID") for row in support if row.get("Prsm ID")}),
                    "supporting_prsm_ids": "; ".join(sorted({row.get("Prsm ID") for row in support if row.get("Prsm ID")})),
                    "supporting_scans": "; ".join(sorted({row.get("Scan(s)") for row in support if row.get("Scan(s)")})),
                    "prsm_facts": [{
                        "prsmId": row.get("Prsm ID") or "",
                        "featureId": row.get("Feature ID") or "",
                        "proteoformRef": "|".join((
                            row.get("Feature ID") or "",
                            row.get("Proteoform ID") or "",
                            row.get("Protein accession") or "",
                        )),
                        "sourceProforma": row.get("Proteoform") or "",
                        "fixedPtms": row.get("Fixed PTMs") or "",
                        "variablePtms": row.get("Variable PTMs") or "",
                        "unexpectedModificationsRaw": row.get("Unexpected modifications") or "",
                        "scan": row.get("Scan(s)") or "",
                        "firstResidue": int(number(row.get("First residue"), 1)),
                        "lastResidue": int(number(row.get("Last residue"), 0)),
                        "proteinNTerminalForm": row.get("Protein N-terminal form") or "NONE",
                        "proteoformMass": number(row.get("Proteoform mass")),
                        "precursorMass": number(row.get("Precursor mass")),
                        "adjustedPrecursorMass": number(row.get("Adjusted precursor mass")),
                        "fragmentation": row.get("Fragmentation") or "",
                        "retentionTime": number(row.get("Retention time")),
                        "unexpectedModifications": [
                            {"value": observed_mass,
                             "frame": "unexpected_shift_after_fixed_forms",
                             "interval": [observed_start, observed_end]}
                            for observed_mass, observed_start, observed_end, _
                            in parse_mods(row.get("Unexpected modifications"))
                        ],
                        "evalue": number(row.get("E-value"), 1.0),
                        "spectrumQValue": number(
                            row.get("Spectrum-level Q-value"), None),
                        "proteoformQValue": number(
                            row.get("Proteoform-level Q-value"), None),
                        "matchedMasses": int(number(row.get("# Matched masses"))),
                        "matchedFragments": int(number(row.get("# Matched fragments"))),
                        "sequenceCoverage": None,
                    } for row in all_support],
                    "fragmentations": "; ".join(sorted({row.get("Fragmentation") for row in support if row.get("Fragmentation")})),
                    "best_evalue": min((number(row.get("E-value"), 1.0) for row in support), default=number(pf_row.get("E-value"), 1.0)),
                    "best_fragments": max((int(number(row.get("# Matched fragments"))) for row in support), default=int(number(pf_row.get("# Matched fragments")))),
                    "ms2_support": "yes",
                    "seq_gap": seq_gap,
                })
    # ms1 masses here are the C13-normalized products of merge_source_features,
    # so the feature-family ladder search needs no isotope slack of its own.
    ms1_families, family_search_coverage = find_feature_families(
        ms1, identified, ppm_tol=ppm_tol)
    join_summary["ms1FamilySearch"] = family_search_coverage
    families = list(ms1_families)
    families.extend(find_identified_serial_families(
        ms1, prsms, qvalue_max=qvalue_max,
        evalue_fallback=evalue_fallback, ppm_tol=ppm_tol,
        n_max_c13=n_max_c13,
    ))
    cluster_modification_masses(modification_rows, ppm_tol)
    feature_rows.sort(key=lambda row: (float(row["mass"]), float(row["rt_start"])))
    return feature_rows, modification_rows, families, join_summary


def split_values(value):
    return {item for item in (value or "").split("; ") if item}

def _proforma_with_fixed_forms(sequence, fixed_forms):
    """Materialize the record's named terminal fixed forms into the sequence.

    Name tags only, never hand-written masses: "[Acetyl]-SEQ" for an nTerm
    form, "SEQ-[Amidated]" for cTerm. The result is the frame-correct base
    for refine-shift — TopPIC's unexpected shift is measured after these
    forms, so the calibration anchor must sit on sequence+forms.
    """
    prefix = sorted({form["name"] for form in (fixed_forms or [])
                     if form.get("position") == "nTerm" and form.get("name")})
    suffix = sorted({form["name"] for form in (fixed_forms or [])
                     if form.get("position") == "cTerm" and form.get("name")})
    out = sequence or ""
    for name in prefix:
        out = f"[{name}]-{out}"
    for name in suffix:
        out = f"{out}-[{name}]"
    return out

def ordered_values(value):
    return [item for item in (value or "").split("; ") if item]

def aggregate_checklist(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["mass_shift_id"]].append(row)
    output = []
    for mass_shift_id, members in grouped.items():
        # A shared measured shift does not establish a shared chemical identity.
        # The identified feature is the physical carrier boundary; candidate
        # chemistry cannot split, merge, retain, or prioritize it.
        by_feature = defaultdict(list)
        for row in members:
            by_feature[row["feature_id"]].append(row)
        for feature_id, feature_rows in by_feature.items():
            assignment_id = f"{mass_shift_id}:F{feature_id}"
            output.append(_aggregate_partition(mass_shift_id, assignment_id, feature_rows))
    output.sort(key=lambda row: (-row["prsm_count"], float(row["best_evalue"])))
    return output


def _aggregate_partition(mass_shift_id, assignment_id, members):
        feature_ids = sorted({row["feature_id"] for row in members}, key=lambda value: int(value))
        proteoform_refs = sorted({row["proteoform_ref"] for row in members if row["proteoform_ref"]})
        smith_by_ref = {row["proteoform_ref"]: int(row["smith_level"])
                        for row in members if row["proteoform_ref"]}
        smith_counts = Counter(smith_by_ref.values())
        prsm_ids = sorted({item for row in members for item in split_values(row["supporting_prsm_ids"])})
        localizations = sorted({row["localizations"] for row in members if row["localizations"]})
        residues = sorted({row["localization_residues"] for row in members if row["localization_residues"]})
        carrier_facts = []
        for feature_id in feature_ids:
            feature_rows = [row for row in members if row["feature_id"] == feature_id]
            fixed_forms = {
                (form.get("position"), form.get("name"))
                for row in feature_rows
                for form in (row.get("fixed_forms_excluded") or [])
                if form.get("position") and form.get("name")
            }
            # Frame-correct calibration base: the same named fixed forms the
            # record declares, materialized into the sequence (name tags only,
            # never hand-written masses). TopPIC's unexpected shift is measured
            # after these forms, so the refine anchor must sit on
            # sequence+forms, not the bare slice.
            sequences_with_forms = sorted({
                _proforma_with_fixed_forms(row["sequence"], row.get("fixed_forms_excluded"))
                for row in feature_rows if row.get("sequence")
            })
            carrier_facts.append({
                "featureId": feature_id,
                "proteins": sorted({row["protein"] for row in feature_rows if row["protein"]}),
                "proteoformRefs": sorted({
                    row["proteoform_ref"] for row in feature_rows if row["proteoform_ref"]
                }),
                "sequences": sorted({row["sequence"] for row in feature_rows if row.get("sequence")}),
                "sequencesWithForms": sequences_with_forms,
                "rawObservedShifts": sorted({
                    round(float(row.get("raw_shift_mass", row["shift_mass"])), 6)
                    for row in feature_rows
                }),
                "fixedFormsExcluded": [
                    {"position": position, "name": name}
                    for position, name in sorted(fixed_forms)
                ],
                "featureMass": round(float(feature_rows[0]["proteoform_mass"]), 6),
                "retentionTime": {
                    "start": round(min(float(row["rt_start"]) for row in feature_rows), 6),
                    "apex": round(statistics.median(
                        float(row["rt_apex"]) for row in feature_rows
                    ), 6),
                    "end": round(max(float(row["rt_end"]) for row in feature_rows), 6),
                    "unit": "min",
                },
                "localizations": sorted({
                    row["localizations"] for row in feature_rows if row["localizations"]
                }),
                "residues": sorted({
                    row["localization_residues"] for row in feature_rows
                    if row["localization_residues"]
                }),
                "prsmIds": sorted({
                    item for row in feature_rows
                    for item in split_values(row.get("supporting_prsm_ids", ""))
                }),
                "prsmFacts": sorted(
                    {
                        (fact["prsmId"], fact["proteoformRef"]): fact
                        for row in feature_rows
                        for fact in (row.get("prsm_facts") or [])
                        if fact.get("prsmId")
                    }.values(),
                    key=lambda fact: (fact["proteoformMass"], fact["prsmId"]),
                ),
                "scans": sorted({
                    item for row in feature_rows
                    for item in split_values(row.get("supporting_scans", ""))
                }),
                "fragmentations": sorted({
                    item for row in feature_rows
                    for item in split_values(row.get("fragmentations", ""))
                }),
                "bestEvalue": min(float(row["best_evalue"]) for row in feature_rows),
                "bestFragments": max(int(row["best_fragments"]) for row in feature_rows),
            })
        prsm_flags = {bool(split_values(row.get("supporting_prsm_ids", ""))) for row in members}
        ms2_support = "yes" if prsm_flags == {True} else ("none" if prsm_flags == {False} else "partial")
        # Ingest has no chemistry candidate. PTM identity therefore remains
        # open; localization presence only distinguishes level 4 from 5.
        status = "4" if localizations else "5"
        carrier_facts = json.dumps(carrier_facts, separators=(",", ":"))
        # Clustering merges C13-equivalent raw shifts into one mass_shift_id
        # (mono-isotopic mispick tolerance) and candidates are searched at the
        # cluster center only — declare the merge and the member spread so a
        # minority raw value relabeled to the center stays visible.
        raw_shifts = sorted({
            round(float(row.get("raw_shift_mass", row["shift_mass"])), 4)
            for row in members
        })
        isotope_merged = len(raw_shifts) > 1 and (
            raw_shifts[-1] - raw_shifts[0] >= 0.5)
        return {
            "assignment_id": assignment_id,
            "mass_shift_id": mass_shift_id,
            "shift_mass": members[0]["shift_mass_center"],
            "cluster_isotope_merged": isotope_merged,
            "cluster_raw_shifts": "; ".join(f"{value:.4f}" for value in raw_shifts)
            if isotope_merged else "",
            "feature_count": len(feature_ids),
            "proteoform_count": len(proteoform_refs),
            "prsm_count": len(prsm_ids),
            "protein_count": len({row["protein"] for row in members if row["protein"]}),
            "feature_ids": "; ".join(feature_ids),
            "proteoform_refs": "; ".join(proteoform_refs),
            "localizations": "; ".join(localizations),
            "localization_residues": "; ".join(residues),
            "smith_level_counts": "; ".join(
                f"{level}:{smith_counts[level]}" for level in sorted(smith_counts)
            ),
            "smith_best_level": min(smith_by_ref.values()) if smith_by_ref else "",
            "smith_worst_level": max(smith_by_ref.values()) if smith_by_ref else "",
            "smith_level_1_2_count": sum(count for level, count in smith_counts.items() if level <= 2),
            "smith_level_3_count": smith_counts[3],
            "smith_level_4_5_count": sum(count for level, count in smith_counts.items() if level >= 4),
            "smith_ambiguities": "; ".join(sorted({
                item for row in members for item in split_values(row["smith_ambiguities"])
            })),
            "best_evalue": min(float(row["best_evalue"]) for row in members),
            "best_fragments": max(int(row["best_fragments"]) for row in members),
            "ms2_support": ms2_support,
            "status": status,
            "carrier_facts": carrier_facts,
            "seq_gap": "; ".join(sorted({
                row["seq_gap"] for row in members if row.get("seq_gap")
            })),
            "calibration_status": "unavailable",
            "calibrated_shift": "",
            "isotope_offset": "",
            "envelope_fit": "",
            "evidence_refs": "",
            "notes": "",
        }
