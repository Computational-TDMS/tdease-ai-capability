"""Read-only compact discovery summary for one or more TopPIC PrSM tables."""

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

# Standalone-script invocation (uv run --script) puts only scripts/ on sys.path;
# _shared lives in the skills parent (same bootstrap as figures.py), so the
# path fix must precede the local imports.
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from toppic_io import parse_mods, prsm_passes, read_ms1, read_toppic
from _shared import native


def sibling_path(source, old_suffix, new_suffix):
    path = Path(source)
    name = path.name
    return path.with_name(name[:-len(old_suffix)] + new_suffix) if name.endswith(old_suffix) else None


def sibling_feature_path(source):
    """Resolve the TopFD feature sibling for either TopPIC PrSM filename."""
    for suffix in ("_ms2_toppic_prsm_single.tsv", "_ms2_toppic_prsm.tsv"):
        candidate = sibling_path(source, suffix, "_ms1.feature")
        if candidate is not None:
            return candidate
    return None


def resolve_prsm_sources(prsms, run_dir):
    """Resolve one primary PrSM table per sample without row-level parsing."""
    if run_dir:
        root = Path(run_dir)
        if not root.is_dir():
            raise ValueError(f"run directory does not exist: {root}")
        single = sorted(root.glob("*_ms2_toppic_prsm_single.tsv"))
        combined = sorted(root.glob("*_ms2_toppic_prsm.tsv"))
        if single:
            return [str(path) for path in single]
        if combined:
            return [str(path) for path in combined]
        raise ValueError(f"no TopPIC primary PrSM tables in run directory: {root}")
    if not prsms:
        raise ValueError("provide --prsms or --run-dir")
    return list(prsms)


def window_intersection(windows):
    """Return the shared TopPIC localization corridor, when one exists."""
    parsed = []
    for value in windows:
        try:
            start, end = (int(item) for item in value.split("-", 1))
        except (AttributeError, ValueError):
            continue
        parsed.append((start, end))
    if not parsed:
        return None
    start = max(item[0] for item in parsed)
    end = min(item[1] for item in parsed)
    return f"{start}-{end}" if start <= end else None


def representative_record(row):
    """Project one PrSM row into the representative record (top families only)."""
    return {
        "dataFile": row.get("Data file name"),
        "sourcePrsm": row.get("_sourcePrsm"),
        "prsmId": row.get("Prsm ID") or "",
        "scan": row.get("Scan(s)"),
        "fragmentation": row.get("Fragmentation"),
        "retentionTime": row.get("Retention time"),
        "charge": row.get("Charge"),
        "proteoformId": row.get("Proteoform ID"),
        "featureId": row.get("Feature ID"),
        "proteinAccession": row.get("Protein accession"),
        "proteoformMass": row.get("Proteoform mass"),
        # The feature-group fallbacks read these: without the projection the
        # fallback silently yields nulls even when the source table has them.
        "precursorMass": row.get("Precursor mass"),
        "featureApexTime": row.get("Feature apex time"),
        "firstResidue": row.get("First residue"),
        "lastResidue": row.get("Last residue"),
        "databaseSequence": row.get("Database protein sequence"),
        "sourceProforma": row.get("Proteoform"),
        "fixedPtms": row.get("Fixed PTMs"),
        "variablePtms": row.get("Variable PTMs"),
        "terminalForm": row.get("Protein N-terminal form"),
        "unexpectedModification": row.get("Unexpected modifications"),
        "matchedMasses": row.get("# Matched masses"),
        "matchedFragments": row.get("# Matched fragments"),
        "evalue": row.get("E-value"),
        "spectrumQvalue": row.get("Spectrum-level Q-value"),
        "proteoformQvalue": row.get("Proteoform-level Q-value"),
    }


def proposal_input(row):
    """Bind source facts to native argv; never reconstruct chemistry or a shell string."""
    source = row.get("sourcePrsm")
    msalign = None
    if source:
        for suffix in ("_ms2_toppic_prsm_single.tsv", "_ms2_toppic_prsm.tsv"):
            candidate = sibling_path(source, suffix, "_ms2.msalign")
            if candidate is not None and candidate.is_file():
                msalign = str(candidate.resolve())
                break
    proforma = row.get("sourceProforma")
    scan = row.get("scan")
    missing = [name for name, value in (
        ("msalign", msalign), ("sourceProforma", proforma), ("scan", scan),
    ) if not value]
    limitation = ("Select one unresolved shift before using resolve-shift; "
                  "this source does not declare exactly one unexpected modification."
                  if len(parse_mods(row.get("unexpectedModification"))) != 1 else None)
    return {
        "argv": (["tdeasenexus-toolbox", "resolve-shift", msalign,
                  "--scan", scan, "--proforma", proforma, "--library", "builtin"]
                 if not missing and not limitation else None),
        "missingInputs": missing,
        "limitation": limitation,
    }


def carrier_summaries(representatives):
    """Keep protein-local localization facts separate from family recurrence."""
    grouped = defaultdict(list)
    for row in representatives:
        grouped[row.get("proteinAccession") or "unknown"].append(row)
    output = []
    for accession, rows in grouped.items():
        windows = [
            f"{start}-{end}"
            for row in rows
            for _mass, start, end, _raw in parse_mods(row.get("unexpectedModification"))
        ]
        values = [
            mass
            for row in rows
            for mass, _start, _end, _raw in parse_mods(row.get("unexpectedModification"))
        ]
        best = min(rows, key=lambda row: float(row.get("evalue") or "inf"))
        output.append({
            "proteinAccession": accession,
            "prsmCount": len({
                (row.get("sourcePrsm"), row.get("prsmId")) for row in rows
            }),
            "featureCount": len({
                (row.get("sourcePrsm"), row.get("featureId")) for row in rows
            }),
            "meanShiftDa": sum(values) / len(values),
            "localizationIntersection": window_intersection(windows),
            "terminalSequence": (best.get("databaseSequence") or "")[:24],
            "terminalForms": sorted({row.get("terminalForm") or "NONE" for row in rows}),
            "bestSourcePrsm": best.get("sourcePrsm"),
            "bestPrsmId": best.get("prsmId"),
            "bestEvalue": best.get("evalue"),
        })
    return sorted(output, key=lambda row: (-row["prsmCount"], float(row["bestEvalue"] or "inf")))


def attention_score(family):
    """Rank compact, recurring, carrier-coherent families for first-pass review.

    Log-product formula over PrSM-family facts. It ranks discovery families;
    it is not a chemistry confidence score.
    """
    count = len(family["prsmIds"])
    proteins = len(family["proteins"])
    proteoforms = len(family["proteoformIds"])
    windows = family["windows"]
    concentration = max(windows.values(), default=0) / max(count, 1)
    carrier_counts = list(family["proteins"].values())
    smallest_carrier = min(carrier_counts, default=0)
    carrier_balance = smallest_carrier / max(count / max(proteins, 1), 1)
    # Broad routine shifts spread across many proteins; a repeated, localized
    # family deserves attention even when its raw count is lower.
    return (
        math.log1p(count)
        * math.log1p(proteins)
        * math.log1p(proteoforms)
        * concentration
        * min(carrier_balance, 1.0)
        / math.sqrt(max(proteins, 1))
    )


def shift_series(families, tolerance=0.15, min_step=5.0, max_level=6, limit=3):
    """Find recurring arithmetic relations between observed shift families.

    This is deliberately a mass-space relation, not a modification classifier.
    Candidate step sizes come from observed family differences and are retained
    only when two or more levels share a carrier protein.  The output points
    back to the existing family bins; it does not create another data store.
    """
    nodes = []
    for shift, family in families.items():
        mean = sum(family["values"]) / len(family["values"])
        if abs(mean) < min_step:
            continue
        nodes.append({
            "shiftBinDa": shift,
            "meanShiftDa": mean,
            "magnitudeDa": abs(mean),
            "direction": 1 if mean > 0 else -1,
            "prsmIds": family["prsmIds"],
            "proteins": set(family["proteins"]),
            "sourceRefs": {
                row.get("sourcePrsm") for row in family["representatives"]
                if row.get("sourcePrsm")
            },
        })
    if len(nodes) < 2:
        return []

    # Derive possible steps from observations that share a carrier.  The old
    # implementation formed every pair of global shift bins and only later
    # asked whether the resulting relation had a common protein.  That made a
    # broad Run quadratic in all unrelated residuals.  Carrier-local pairs
    # preserve the same relation while keeping the search proportional to the
    # actual evidence that can support a serial family.
    steps = set()
    by_protein = defaultdict(list)
    for node in nodes:
        for protein in node["proteins"]:
            by_protein[protein].append(node)
    for protein_nodes in by_protein.values():
        for direction in (-1, 1):
            directed = sorted(
                (node for node in protein_nodes if node["direction"] == direction),
                key=lambda node: node["magnitudeDa"],
            )
            for left_index, left in enumerate(directed):
                for right in directed[left_index + 1:]:
                    difference = right["magnitudeDa"] - left["magnitudeDa"]
                    if difference <= tolerance:
                        continue
                    for divisor in range(1, max_level):
                        steps.add((direction, round(difference / divisor, 6)))

    candidates = []
    for direction, step in steps:
        if step <= tolerance:
            continue
        matched = []
        for node in nodes:
            if node["direction"] != direction:
                continue
            level = max(1, round(node["magnitudeDa"] / step))
            if level > max_level:
                continue
            residual = abs(node["magnitudeDa"] - level * step)
            if residual <= tolerance:
                matched.append((level, residual, node))
        by_level = {}
        for level, residual, node in matched:
            current = by_level.get(level)
            if current is None or residual < current[0]:
                by_level[level] = (residual, node)
        if len(by_level) < 2:
            continue

        levels = sorted(by_level)
        shared = set.intersection(*(by_level[level][1]["proteins"] for level in levels))
        if not shared:
            continue
        consecutive = sum(
            1 for left, right in zip(levels, levels[1:]) if right == left + 1
        )
        if not consecutive:
            continue
        weights = {
            level: len(by_level[level][1]["prsmIds"])
            for level in levels
        }
        fitted_step = sum(
            by_level[level][1]["magnitudeDa"] / level * weights[level]
            for level in levels
        ) / sum(weights.values())
        residual = sum(
            abs(by_level[level][1]["magnitudeDa"] - level * fitted_step)
            * weights[level]
            for level in levels
        ) / sum(weights.values())
        prsm_ids = set().union(*(by_level[level][1]["prsmIds"] for level in levels))
        source_refs = set().union(*(
            by_level[level][1]["sourceRefs"] for level in levels
        ))
        baseline_observed = 1 in levels
        # Prefer a relation anchored at its first observed level, then use
        # recurrence and adjacency.  This suppresses d/2 aliases while keeping
        # incomplete ladders available behind complete-from-level-one ladders.
        score = (
            (1000 if baseline_observed else 0)
            + len(prsm_ids) * 5
            + len(levels) * 20
            + consecutive * 25
            + len(shared) * 10
            - residual * 10
        )
        candidates.append({
            "stepDa": direction * fitted_step,
            "levels": [{
                "level": level,
                "shiftBinDa": by_level[level][1]["shiftBinDa"],
                "meanShiftDa": by_level[level][1]["meanShiftDa"],
                "prsmCount": len(by_level[level][1]["prsmIds"]),
            } for level in levels],
            "sharedProteins": sorted(shared),
            "levelCount": len(levels),
            "baselineObserved": baseline_observed,
            "consecutiveLevelLinks": consecutive,
            "prsmCount": len(prsm_ids),
            "sourceRefs": sorted(source_refs),
            "residualDa": residual,
            "attentionScore": score,
            "measurementRole": "shift_series_prior",
        })

    # Multiple divisors can describe the same observed levels.  Keep the best
    # relation for a level set and carrier set before returning a compact view.
    candidates.sort(key=lambda item: (-item["attentionScore"], item["stepDa"]))
    selected = []
    seen = set()
    for item in candidates:
        key = (
            tuple(level["shiftBinDa"] for level in item["levels"]),
            tuple(item["sharedProteins"]),
        )
        if key in seen:
            continue
        seen.add(key)
        item.pop("attentionScore", None)
        selected.append(item)
        if len(selected) >= max(0, limit):
            break
    return selected


def row_mods(row):
    """Parsed 'Unexpected modifications' for one row, cached on the row dict:
    the same rows flow through family grouping, background candidates and
    anchors, and the regex parse is the per-row hot cost."""
    cached = row.get("_mods")
    if cached is None:
        cached = parse_mods(row.get("Unexpected modifications"))
        row["_mods"] = cached
    return cached


def unimod_priors(mass, tolerance, limit, _binary=None, _cache={}):
    """Mass-only UniMod candidates (display overlays, not observations)."""
    key = round(mass, 6)
    if key not in _cache:
        doc = native.run_json([
            native.binary(None), "unimod-search", "--mass", str(mass),
            "--tolerance", str(tolerance), "--tolerance-unit", "da",
            "--keep-classification", "Post-translational",
            "--limit", str(limit),
        ], label="tdeasenexus-toolbox unimod-search")
        _cache[key] = [{
            "title": match["title"],
            "theoreticalDeltaMass": match["theoreticalDeltaMass"],
            "absoluteErrorDa": match["absoluteErrorDa"],
        } for match in doc.get("result", {}).get("massMatches", [])
          if match.get("absoluteErrorDa", float("inf")) <= tolerance]
    return _cache[key]


def annotate_serial_spacings(families, tolerance=0.02, limit=3):
    """Attach mass-only UniMod priors to measured proteoform ladders."""
    for family in families:
        family["unimodPriors"] = unimod_priors(family["spacing_center_da"],
                                               tolerance, limit)


def annotate_shift_series(series, tolerance=0.02, limit=3):
    """Attach the same mass-only UniMod priors to recurring shift steps."""
    projected = [
        {"spacing_center_da": item["stepDa"]}
        for item in series
    ]
    annotate_serial_spacings(projected, tolerance=tolerance, limit=limit)
    for item, annotation in zip(series, projected):
        item["unimodPriors"] = annotation["unimodPriors"]


def background_candidates(rows, limit=6, tolerance=0.02, unimod_limit=3):
    """Project recurring carrier-scoped shifts as background candidates.

    A recurrence may be shared background, the declared variable family, or a
    separate physical population. This projection surfaces those alternatives
    (with terminal forms) before the model is declared; it does not choose one.
    """
    terminal_forms = Counter(
        row.get("Protein N-terminal form") or "unknown" for row in rows
    )
    total_prsms = len({
        (row.get("_sourcePrsm"), row.get("Prsm ID") or "") for row in rows
    })
    clusters = defaultdict(lambda: {
        "prsmIds": set(), "proteoformIds": set(), "values": [],
        "intervals": [], "sequence": None, "representative": None,
    })
    for row in rows:
        prsm_ref = (row.get("_sourcePrsm"), row.get("Prsm ID") or "")
        for value, start, end, _ in row_mods(row):
            cluster = clusters[round(value, 1)]
            cluster["prsmIds"].add(prsm_ref)
            if row.get("Proteoform ID"):
                cluster["proteoformIds"].add((row.get("_sourcePrsm"), row["Proteoform ID"]))
            cluster["values"].append(value)
            cluster["intervals"].append((start, end))
            if cluster["sequence"] is None:
                cluster["sequence"] = row.get("Database protein sequence") or ""
                cluster["representative"] = representative_record(row)
    candidates = []
    for shift, cluster in sorted(
        clusters.items(), key=lambda item: -len(item[1]["prsmIds"])
    ):
        if len(cluster["prsmIds"]) < 2 or len(candidates) >= limit:
            continue
        mass = sum(cluster["values"]) / len(cluster["values"])
        start = max(low for low, _ in cluster["intervals"])
        end = min(high for _, high in cluster["intervals"])
        converged = start <= end
        priors = unimod_priors(mass, tolerance, unimod_limit)
        windows = Counter(f"{low}-{high}" for low, high in cluster["intervals"])
        sequence = cluster["sequence"] or ""
        coverage = (len(cluster["prsmIds"]) / total_prsms) if total_prsms else 0.0
        candidates.append({
            "shiftBinDa": shift,
            "meanShiftDa": mass,
            "prsmCount": len(cluster["prsmIds"]),
            "proteoformCount": len(cluster["proteoformIds"]),
            "representativePrsm": cluster["representative"],
            "populationCoverage": round(coverage, 3),
            "windowConverged": converged,
            "windowIntersection": f"{start}-{end}" if converged else None,
            "windowResidues": sequence[start - 1:end] if converged and sequence else None,
            "candidateResidues": (
                [{"residue": residue, "position": start + offset}
                 for offset, residue in enumerate(sequence[start - 1:end])]
                if converged and sequence else None
            ),
            "topWindows": windows.most_common(3),
            "unimodSinglePriors": priors,
            "needsCombination": not priors,
        })
        representative = candidates[-1]["representativePrsm"]
        representative["proposalInput"] = proposal_input(representative)
    return {
        "terminalForms": terminal_forms.most_common(4),
        "qualifiedPrsmsInScope": total_prsms,
        "coverageBasis": "explicit shift annotations among qualified carrier PrSMs; "
                         "not modification prevalence or physical population membership",
        "recurringShifts": candidates,
    }


def serial_projection(family):
    """Public schema for one identified proteoform ladder."""
    return {
        "evidenceType": family["evidence_type"],
        "chemistryStatus": family["chemistry_status"],
        "familyId": family.get("family_id"),
        "proteinAccession": family["protein"],
        "memberCount": family["member_count"],
        "prsmCount": family["prsm_count"],
        "featureIds": family["feature_ids"],
        "featureRefs": [{
            "sourcePrsm": family["sourcePrsm"],
            "featureId": str(feature_id),
        } for feature_id in family["feature_ids"]],
        "masses": family["masses"],
        "rawMassGaps": family["raw_mass_gaps"],
        "massGaps": family["mass_gaps"],
        "isotopeOffsets": family["isotope_offsets"],
        "spacingDa": family["spacing_center_da"],
        "spacingSpreadDa": family["spacing_spread_da"],
        "gapErrorPpmIntact": family["gap_error_ppm_intact"],
        "rtApexes": family["rt_apexes"],
        "rtDrift": family["rt_drift"],
        "intensityTotal": family["intensity_total"],
        "memberRecords": family["member_records"],
        "unimodPriors": family["unimodPriors"],
        "sourcePrsm": family["sourcePrsm"],
        "massSource": family["massSource"],
    }


def unified_attention_candidates(shift_families, serials, shift_series_items,
                                 tolerance=0.02, limit=8):
    """One compact ranking across residual and identified-series evidence."""
    candidates = []
    for family in serials:
        priors = family["unimodPriors"]
        best_error = min((item["absoluteErrorDa"] for item in priors), default=tolerance)
        prior_quality = max(0.0, 1.0 - best_error / tolerance) if priors else 0.0
        score = (
            5.0 * prior_quality
            + math.log1p(family["prsmCount"])
            + math.log1p(family["memberCount"])
            - family["spacingSpreadDa"] / tolerance
        )
        candidates.append({
            "kind": "identified_proteoform_series",
            "id": family["familyId"],
            "proteinAccession": family["proteinAccession"],
            "massPattern": {"stepDa": family["spacingDa"], "levels": family["memberCount"]},
            "support": {
                "prsmCount": family["prsmCount"],
                "featureRefs": family["featureRefs"],
            },
            "unimodPriors": priors,
            "attentionScore": round(score, 6),
            "sourceRefs": [family["sourcePrsm"]],
            "massSource": family["massSource"],
            "notIdentity": "A repeated mass step and UniMod prior do not establish PTM identity or localization.",
        })
    for family in shift_families:
        score = family["attentionScore"]
        candidates.append({
            "kind": "unexpected_shift_family",
            "id": f"shift:{family['shiftBinDa']}",
            "proteinAccessions": [item[0] for item in family["topProteins"]],
            "massPattern": {"meanShiftDa": family["meanShiftDa"]},
            "support": {
                "prsmCount": family["prsmCount"],
                "proteoformCount": family["proteoformCount"],
                "featureRefs": [{
                    "sourcePrsm": group.get("sourcePrsm"),
                    "featureId": group.get("featureId"),
                } for group in family["featureGroups"]],
            },
            "attentionScore": round(score, 6),
            "sourceRefs": sorted({
                group["sourcePrsm"] for group in family["featureGroups"]
                if group.get("sourcePrsm")
            }),
            "notIdentity": "A recurring unexpected shift is an attention candidate, not a chemical identity.",
        })
    for series in shift_series_items:
        priors = series.get("unimodPriors", [])
        best_error = min((item["absoluteErrorDa"] for item in priors), default=tolerance)
        prior_quality = max(0.0, 1.0 - best_error / tolerance) if priors else 0.0
        score = (
            5.0 * prior_quality
            + math.log1p(series["prsmCount"])
            + math.log1p(series["levelCount"])
            + series["consecutiveLevelLinks"]
            - series["residualDa"] / max(tolerance, 1e-9)
        )
        candidates.append({
            "kind": "repeating_shift_series",
            "id": f"series:{series['stepDa']:.4f}",
            "proteinAccessions": series["sharedProteins"],
            "massPattern": {
                "stepDa": series["stepDa"],
                "levels": series["levelCount"],
            },
            "support": {
                "prsmCount": series["prsmCount"],
                "levelSummary": series["levels"],
            },
            "unimodPriors": priors,
            "attentionScore": round(score, 6),
            "sourceRefs": series.get("sourceRefs", []),
            "notIdentity": "A repeated shift step proposes a modification series; it does not identify its chemistry or sites.",
        })
    candidates.sort(key=lambda item: (-item["attentionScore"], item["kind"], item["id"]))
    return candidates[:max(0, limit)]


def terminal_frame_equivalences(rows, anchors=None, mass_tolerance=0.02):
    """Report alternate TopPIC frames of one measured proteoform mass.

    Read-side inspect view with its own schema and anchor handling;
    ingest.terminal_frame_evidence is the durable-projection writer for the
    same physics — two contracts, kept separate on purpose.
    """
    anchor_keys = None
    if anchors is not None:
        anchor_keys = [
            (row.get("_sourcePrsm"),
             row.get("Protein accession") or "unknown",
             float(row.get("Proteoform mass") or 0))
            for row in anchors if row.get("Proteoform mass")
        ]
    groups = []
    ordered = sorted(
        (row for row in rows if row.get("Proteoform mass")),
        key=lambda row: (
            row.get("_sourcePrsm") or "",
            row.get("Protein accession") or "",
            float(row.get("Proteoform mass") or 0),
        ),
    )
    for row in ordered:
        accession = row.get("Protein accession") or "unknown"
        source_prsm = row.get("_sourcePrsm")
        mass = float(row.get("Proteoform mass") or 0)
        group = next((item for item in reversed(groups)
                      if item["protein"] == accession
                      and item["sourcePrsm"] == source_prsm
                      and abs(item["proteoformMass"] - mass) <= mass_tolerance), None)
        if group is None:
            group = {
                "protein": accession,
                "sourcePrsm": source_prsm,
                "proteoformMass": mass,
                "rows": [],
            }
            groups.append(group)
        group["rows"].append(row)

    output = []
    for group in groups:
        if anchor_keys is not None and not any(
            source_prsm == group["sourcePrsm"]
            and accession == group["protein"]
            and abs(mass - group["proteoformMass"]) <= mass_tolerance
            for source_prsm, accession, mass in anchor_keys
        ):
            continue
        representations = defaultdict(list)
        for row in group["rows"]:
            mods = row_mods(row)
            if not mods:
                continue
            value = mods[0][0]
            key = (
                int(row.get("First residue") or 1),
                row.get("Protein N-terminal form") or "NONE",
            )
            representations[key].append((value, row))
        terminal_keys = {(key[0], key[1]) for key in representations}
        if len(terminal_keys) < 2:
            continue
        frame_values = {
            key: [value for value, _row in members]
            for key, members in representations.items()
        }
        frame_means = {
            key: sum(values) / len(values)
            for key, values in frame_values.items()
        }
        output.append({
            "protein": group["protein"],
            "sourcePrsm": group["sourcePrsm"],
            "proteoformMass": round(group["proteoformMass"], 6),
            "representations": [{
                "firstResidue": key[0],
                "terminalForm": key[1],
                "meanObservedShiftDa": round(frame_means[key], 6),
                "shiftRangeDa": [round(min(frame_values[key]), 6), round(max(frame_values[key]), 6)],
                "prsmCount": len(members),
                "representativePrsmId": min(
                    (row for _value, row in members),
                    key=lambda row: float(row.get("E-value") or "inf"),
                ).get("Prsm ID"),
                "representativeScan": min(
                    (row for _value, row in members),
                    key=lambda row: float(row.get("E-value") or "inf"),
                ).get("Scan(s)"),
            } for key, members in sorted(representations.items())],
            "observedFrameDifferences": [{
                "from": {"firstResidue": left[0], "terminalForm": left[1]},
                "to": {"firstResidue": right[0], "terminalForm": right[1]},
                "shiftDifferenceDa": round(frame_means[right] - frame_means[left], 6),
            } for index, left in enumerate(sorted(frame_means))
              for right in sorted(frame_means)[index + 1:]],
            "measurementRole": "terminal_frame_equivalence",
            "notIdentity": "Equivalent TopPIC mass frames do not identify terminal processing or PTM chemistry.",
        })
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument("--prsms", nargs="+",
                              help="explicit primary TopPIC PrSM tables")
    source_group.add_argument("--run-dir",
                              help="Run directory; choose one primary PrSM table per sample")
    parser.add_argument("--qvalue-max", type=float, default=0.01)
    parser.add_argument("--evalue-fallback", type=float, default=1e-5)
    parser.add_argument("--top", type=int, default=3)
    parser.add_argument("--protein", help="accession or accession substring")
    parser.add_argument("--prsm-id", help="read one exact PrSM from one --prsms table; "
                        "return source facts and proposal inputs without a census or quality exclusion")
    parser.add_argument("--shift", type=float, help="center of one observed shift family")
    parser.add_argument("--shift-tolerance", type=float, default=0.1)
    parser.add_argument("--around-shift", type=float,
                        help="rounded user-reported shift; retains all observed clusters in a broad window")
    parser.add_argument("--around-window", type=float, default=5.0,
                        help="Da half-window for --around-shift (default: 5)")
    parser.add_argument("--representatives", type=int, default=1)
    parser.add_argument("--feature-groups", type=int, default=2,
                        help="representative feature groups per shift cluster (default: 2)")
    parser.add_argument("--terminal-frame-limit", type=int, default=3,
                        help="terminal frame-equivalence records to return (default: 3)")
    parser.add_argument("--series-tolerance", type=float, default=0.15,
                        help="mass tolerance for recurring shift relations (default: 0.15 Da)")
    parser.add_argument("--series-top", type=int, default=3,
                        help="number of recurring shift relations to return (default: 3)")
    parser.add_argument("--proteoform-series-top", type=int, default=12,
                        help="number of identified proteoform ladders to return (default: 12)")
    parser.add_argument("--out", help="optional path for the complete JSON report")
    parser.add_argument("--summary-out", help="optional path for the bounded Agent summary")
    parser.add_argument("--background-out",
                        help="optional path for the carrier shared-background projection "
                             "(requires --protein; source-linked proposals, NOT validated fixed tags)")
    args = parser.parse_args(argv)
    if args.prsm_id is not None:
        if not args.prsms or len(args.prsms) != 1:
            parser.error("PRSM_SOURCE_REQUIRED: --prsm-id requires exactly one --prsms table; IDs are sample-local")
        if args.protein or args.shift is not None or args.around_shift is not None or args.background_out:
            parser.error("PRSM_SCOPE: --prsm-id cannot be combined with census filters or --background-out")
        source = args.prsms[0]
        try:
            matches = [row for row in read_toppic(source) if row.get("Prsm ID") == args.prsm_id]
        except (OSError, ValueError) as error:
            parser.error(f"PRSM_SOURCE_INVALID: {error}")
        if not matches:
            parser.error(f"PRSM_NOT_FOUND: {args.prsm_id} in {source}")
        if len(matches) != 1:
            parser.error(f"PRSM_AMBIGUOUS: {args.prsm_id} has {len(matches)} rows in {source}; provide a primary single table")
        row = matches[0]
        record = representative_record({**row, "_sourcePrsm": source})
        record["proposalInput"] = proposal_input(record)
        report = {
            "prsm": record,
            "quality": {"passes": prsm_passes(row, args.qvalue_max, args.evalue_fallback),
                        "qvalueMax": args.qvalue_max, "evalueFallback": args.evalue_fallback},
        }
        for destination in (args.out, args.summary_out):
            if destination:
                path = Path(destination)
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(report, separators=(",", ":")))
        return 0
    if args.shift is not None and args.around_shift is not None:
        parser.error("--shift and --around-shift are mutually exclusive")
    if args.around_window <= 0:
        parser.error("--around-window must be positive")
    if args.background_out and not args.protein:
        parser.error("--background-out requires --protein")
    explicit_target = args.shift is not None or args.around_shift is not None

    try:
        sources = resolve_prsm_sources(args.prsms, args.run_dir)
    except ValueError as error:
        parser.error(str(error))

    rows = []
    feature_maps = {}
    serial_families = []
    for source in sources:
        source_rows = [
            row for row in read_toppic(source)
            if prsm_passes(row, args.qvalue_max, args.evalue_fallback)
        ]
        for row in source_rows:
            row["_sourcePrsm"] = source
        rows.extend(source_rows)
        feature_path = sibling_feature_path(source)
        feature_map = read_ms1(feature_path) if feature_path and feature_path.is_file() else {}
        feature_maps[source] = feature_map
        if feature_map and not explicit_target:
            # Deferred and sourced from the families module directly: inspect
            # never pulls the ingest command's import graph.
            from ingest_families import find_identified_serial_families
            for family in find_identified_serial_families(
                feature_map,
                source_rows,
                qvalue_max=args.qvalue_max,
                evalue_fallback=args.evalue_fallback,
                max_families_per_protein=20,
            ):
                family["sourcePrsm"] = source
                family["massSource"] = "topfd_ms1_feature"
                for member in family["member_records"]:
                    member["sourcePrsm"] = source
                serial_families.append(family)
    if args.protein:
        needle = args.protein.casefold()
        rows = [row for row in rows if needle in (row.get("Protein accession") or "").casefold()]
        serial_families = [
            family for family in serial_families
            if needle in (family.get("protein") or "").casefold()
        ]
    source_scope_count = len(rows)
    proteins = Counter(row.get("Protein accession") or "unknown" for row in rows)
    terminal_forms = Counter(row.get("Protein N-terminal form") or "unknown" for row in rows)
    proteoforms = {
        (row.get("_sourcePrsm"), row.get("Proteoform ID"))
        for row in rows if row.get("Proteoform ID")
    }

    families = defaultdict(lambda: {
        "prsmIds": set(), "proteoformIds": set(), "proteins": Counter(),
        "windows": Counter(), "values": [], "minEvalue": None, "representatives": [],
    })
    modified_prsms = set()
    for row in rows:
        prsm_id = row.get("Prsm ID") or ""
        source_prsm = row.get("_sourcePrsm")
        prsm_ref = (source_prsm, prsm_id)
        for value, start, end, _ in row_mods(row):
            modified_prsms.add(prsm_ref)
            if args.shift is not None:
                if abs(value - args.shift) > args.shift_tolerance:
                    continue
                key = args.shift
            elif args.around_shift is not None:
                if abs(value - args.around_shift) > args.around_window:
                    continue
                key = round(value, 1)
            else:
                key = round(value, 1)
            family = families[key]
            family["prsmIds"].add(prsm_ref)
            if row.get("Proteoform ID"):
                family["proteoformIds"].add((source_prsm, row["Proteoform ID"]))
            family["proteins"][row.get("Protein accession") or "unknown"] += 1
            family["windows"][f"{start}-{end}"] += 1
            family["values"].append(value)
            try:
                evalue = float(row.get("E-value") or "nan")
                if family["minEvalue"] is None or evalue < family["minEvalue"]:
                    family["minEvalue"] = evalue
            except ValueError:
                pass
            family["representatives"].append(row)

    ranked = sorted(
        families.items(),
        key=(
            (lambda item: (-len(item[1]["prsmIds"]), item[0]))
            if explicit_target
            else (lambda item: (-attention_score(item[1]), item[0]))
        ),
    )
    output_families = []
    for shift, family in ranked[:max(1, args.top)]:
        values = family["values"]
        representatives = sorted(
            (representative_record(row) for row in family["representatives"]),
            key=lambda row: float(row["evalue"] or "inf"),
        )
        feature_groups = []
        by_feature = defaultdict(list)
        for row in representatives:
            by_feature[(row.get("sourcePrsm"), row.get("featureId") or "unknown")].append(row)
        for (source_prsm, feature_id), linked in by_feature.items():
            feature = feature_maps.get(source_prsm, {}).get(feature_id)
            linked_mods = parse_mods(linked[0].get("unexpectedModification"))
            # A carrier row can hold several unexpected modifications; the
            # family is ranked by ONE of them.  Attribute the feature group
            # to the mod this family is about (nearest to the family shift),
            # never to whichever mod happens to be first in the string.
            shift_value, shift_start, shift_end, _raw = min(
                linked_mods, key=lambda mod: abs(mod[0] - shift))
            feature_groups.append({
                "featureId": feature_id,
                "featureIdMissing": feature_id == "unknown",
                "sourcePrsm": source_prsm,
                "proteinAccession": linked[0].get("proteinAccession"),
                "shiftDa": shift_value,
                "allShiftsDa": [mod[0] for mod in linked_mods],
                "localizationWindow": f"{shift_start}-{shift_end}",
                "observedMass": feature["mass"] if feature else linked[0].get("precursorMass"),
                "rtStart": feature["rt_start"] if feature else None,
                "rtEnd": feature["rt_end"] if feature else None,
                "rtApex": feature["rt_apex"] if feature else linked[0].get("featureApexTime"),
                "scans": [{"scan": row["scan"], "fragmentation": row["fragmentation"], "charge": row["charge"]} for row in linked],
                "terminalSequence": (linked[0].get("databaseSequence") or "")[:24],
                "msalignSource": linked[0].get("dataFile"),
            })
        output_families.append({
            "shiftBinDa": shift,
            "meanShiftDa": sum(values) / len(values),
            "minShiftDa": min(values),
            "maxShiftDa": max(values),
            "prsmCount": len(family["prsmIds"]),
            "proteoformCount": len(family["proteoformIds"]),
            "proteinCount": len(family["proteins"]),
            "topProteins": family["proteins"].most_common(5),
            "topWindows": family["windows"].most_common(5),
            "carriers": carrier_summaries(representatives),
            "minEvalue": family["minEvalue"],
            "attentionScore": round(attention_score(family), 6),
            "featureGroups": feature_groups[:max(0, args.feature_groups)],
            "representativePrsms": [
                {**row, "proposalInput": proposal_input(row)}
                for row in representatives[:max(0, args.representatives)]
            ],
        })

    series = [] if explicit_target else shift_series(
        families,
        tolerance=args.series_tolerance,
        limit=args.series_top,
    )
    annotate_shift_series(series)

    # UniMod is a prior annotation, not required to discover recurrence.  A
    # complete Run can contain hundreds of serial ladders, so annotate only a
    # bounded recurrence shortlist.  Focused shift queries do not need these
    # unrelated ladder priors at all.
    for family in serial_families:
        family.setdefault("unimodPriors", [])
    if not explicit_target:
        serial_families.sort(key=lambda family: (
            -(
                math.log1p(family["prsm_count"])
                + math.log1p(family["member_count"])
                - family["spacing_spread_da"] / 0.02
            ),
            family["protein"],
        ))
        annotation_limit = max(
            args.proteoform_series_top * 2,
            args.top * 4,
            12,
        )
        annotate_serial_spacings(serial_families[:annotation_limit])
    serial_families.sort(key=lambda family: (
        0 if family["unimodPriors"] else 1,
        min((prior["absoluteErrorDa"] for prior in family["unimodPriors"]),
            default=float("inf")),
        -family["member_count"],
        -family["prsm_count"],
        family["spacing_spread_da"],
        family["protein"],
    ))

    public_serials = [
        serial_projection(family)
        for family in serial_families[:max(0, args.proteoform_series_top)]
    ]
    attention_candidates = unified_attention_candidates(
        output_families, public_serials, series, limit=min(8, max(1, args.top)),
    )
    background = background_candidates(rows) if args.protein else None
    report = {
        "source": sources,
        "proteinFilter": args.protein,
        "shiftFilterDa": args.shift,
        "shiftToleranceDa": args.shift_tolerance if args.shift is not None else None,
        "aroundShiftDa": args.around_shift,
        "aroundWindowDa": args.around_window if args.around_shift is not None else None,
        "qualifiedPrsmsInProteinScope": source_scope_count,
        "matchingShiftPrsms": sum(len(item[1]["prsmIds"]) for item in ranked[:max(1, args.top)]),
        "matchingShiftPrsmsScope": f"top {max(1, args.top)} families only; "
                                   "totalFamilyPrsms is the global figure",
        "proteoforms": len(proteoforms),
        "proteins": len(proteins),
        "prsmsWithUnexpectedModification": len(modified_prsms),
        "totalShiftFamilies": len(ranked),
        "totalFamilyPrsms": sum(len(item[1]["prsmIds"]) for item in ranked),
        "terminalForms": dict(terminal_forms.most_common()),
        "backgroundCandidates": background,
        "topProteins": proteins.most_common(args.top),
        "attentionCandidates": attention_candidates,
        "shiftFamilies": output_families,
        "shiftSeries": series,
        "identifiedProteoformSeries": public_serials,
        "rankingMethod": "raw_recurrence" if explicit_target else "cross_evidence_attention",
        "terminalFrameEquivalences": terminal_frame_equivalences(
            rows,
            anchors=[
                row for row in rows
                if (args.shift is not None or args.around_shift is not None) and any(
                    abs(value - (args.shift if args.shift is not None else args.around_shift))
                    <= (args.shift_tolerance if args.shift is not None else args.around_window)
                    for value, _start, _end, _raw in row_mods(row)
                )
            ] if args.shift is not None or args.around_shift is not None else None,
        )[:max(0, args.terminal_frame_limit)],
        "calibrationApplied": False,
        "measurementRole": "discovery_summary",
        "notIdentity": "No C13 calibration or chemical assignment was performed; rounded recurring mass shifts are attention candidates only.",
    }
    if args.background_out:
        background_path = Path(args.background_out)
        background_path.parent.mkdir(parents=True, exist_ok=True)
        background_path.write_text(
            json.dumps(background, indent=2) + "\n", encoding="utf-8")
    if args.out or args.summary_out:
        # The compact summary is a projection of the report, never a
        # replacement for it: persist the complete report first, then derive
        # the summary so every derived fact stays traceable on disk.
        if args.out:
            report_path = Path(args.out)
        else:
            requested = Path(args.summary_out)
            report_path = requested.with_name(
                f"{requested.stem}.report{requested.suffix}")
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        summary_path = Path(args.summary_out) if args.summary_out else None
        compact_attention = [{
            key: candidate[key]
            for key in (
                "kind", "id", "proteinAccession", "proteinAccessions",
                "massPattern", "support", "unimodPriors", "attentionScore",
                "sourceRefs",
            )
            if key in candidate
        } for candidate in report["attentionCandidates"][:3]]
        compact_shift_families = []
        if args.protein:
            compact_shift_families = [{
                "id": f"shift:{family['shiftBinDa']}",
                "meanShiftDa": family["meanShiftDa"],
                "prsmCount": family["prsmCount"],
                "proteoformCount": family["proteoformCount"],
                "topWindows": family["topWindows"][:3],
            } for family in report["shiftFamilies"][:min(10, max(1, args.top))]]

        def coverage_entry(total, returned):
            return {"total": total, "returned": returned, "truncated": returned < total}

        compact = {
            "source": report["source"],
            "qualifiedPrsmsInProteinScope": report["qualifiedPrsmsInProteinScope"],
            "totalShiftFamilies": report["totalShiftFamilies"],
            "topProteins": report["topProteins"][:3],
            "attentionCandidates": compact_attention,
            "terminalFrameEquivalences": report["terminalFrameEquivalences"],
            "report": str(report_path),
            "coverage": {
                "shiftFamilies": coverage_entry(
                    report["totalShiftFamilies"], len(compact_shift_families)),
                "attentionCandidates": coverage_entry(
                    len(report["attentionCandidates"]), len(compact_attention)),
                "terminalFrameEquivalences": coverage_entry(
                    len(report["terminalFrameEquivalences"]),
                    len(report["terminalFrameEquivalences"])),
                "shiftSeries": coverage_entry(len(report["shiftSeries"]), 0),
                "identifiedProteoformSeries": coverage_entry(
                    len(report["identifiedProteoformSeries"]), 0),
                "topProteins": coverage_entry(
                    len(report["topProteins"]), min(3, len(report["topProteins"]))),
            },
        }
        if args.protein:
            compact["backgroundCandidates"] = background
            compact["shiftFamilies"] = compact_shift_families
        if summary_path:
            summary_path.parent.mkdir(parents=True, exist_ok=True)
            summary_path.write_text(json.dumps(compact, indent=2) + "\n", encoding="utf-8")
            print(json.dumps({"summary": str(summary_path), **compact}, separators=(",", ":")))
        else:
            print(json.dumps(compact, separators=(",", ":")))
    else:
        print(json.dumps(report, separators=(",", ":")))


if __name__ == "__main__":
    raise SystemExit(main())
