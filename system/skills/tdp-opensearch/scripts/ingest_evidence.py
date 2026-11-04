"""
Ingest evidence and calibration: the residual-queue handoff artifact,
evidence records, UniMod candidate generation, and the embedded MS1
refine-shift calibration pass. The only ingest module that talks to the
tdeasenexus-toolbox binary.
"""

import hashlib
import json
import statistics
import tempfile
from collections import defaultdict
from pathlib import Path

from family_discovery import (
    REFINE_MIN_MATCHED_PEAKS,
    REFINE_MIN_PEARSON,
    REFINE_MIN_SCAN_COUNT,
)
from toppic_io import number, ref_accession
from _shared import native
from ingest_tables import ordered_values


CALIBRATION_RT_HALF_WIDTH = 0.5

# Same physical envelope-refine QC gates as family_discovery.REFINE_* — one
# constant home, aliased to keep the calibration code's vocabulary.
CALIBRATION_MIN_PEARSON = REFINE_MIN_PEARSON
CALIBRATION_MIN_MATCHED_PEAKS = REFINE_MIN_MATCHED_PEAKS
CALIBRATION_MIN_SCAN_COUNT = REFINE_MIN_SCAN_COUNT


def residual_queue_records(rows, families=None, candidates=None):
    """Build the single preprocessing handoff stream in memory.

    The records feed the compact bundle directly; keeping construction here
    prevents a second serialization/readback pass.
    """
    seen = set()
    direct_families = defaultdict(set)
    protein_families = defaultdict(set)
    for family in families or []:
        family_id = family.get("family_id")
        if not family_id:
            continue
        for feature_id in family.get("feature_ids", []):
            direct_families[str(feature_id)].add(family_id)
        protein = family.get("protein")
        if protein:
            protein_families[protein].add(family_id)
    records = []
    # Family records first: discovery facts and routing rows share one stream.
    for family in families or []:
        family_id = family.get("family_id")
        if not family_id:
            continue
        metrics = {key: family[key] for key in (
            "member_count", "masses", "raw_mass_gaps", "mass_gaps",
            "isotope_offsets", "spacing_center_da",
            "spacing_spread_da", "common_rt", "intensity_total",
        ) if key in family}
        if family.get("protein"):
            metrics["protein"] = family["protein"]
        records.append({
                "familyId": str(family_id),
                "familyType": family.get("evidence_type") or "ms1_rt_coherent_ladder",
                "definition": family.get("definition") or "",
                "featureIds": [str(item) for item in family.get("feature_ids", [])],
                "metrics": metrics,
            })
    # Candidate sets ride the same artifact: priors only, never identity.
    for record in candidates or []:
        metrics = record.get("metrics") or {}
        singles = metrics.get("singleModCandidates") or []
        combos = metrics.get("combinationCandidates") or []
        records.append({
                "candidateSet": {
                    "evidenceId": record.get("evidenceId", ""),
                    "kind": (record.get("evidenceId", "").split(":")[1]
                             if record.get("evidenceId", "").count(":") >= 1 else "assignment"),
                    "measurementRole": "candidate_generation",
                    "notIdentity": True,
                    "question": (record.get("assignmentIds") or [None])[0],
                    "familyId": (record.get("familyRefs") or [None])[0],
                    "featureIds": [str(item) for item in record.get("featureIds", [])],
                    "proteoformRefs": list(record.get("proteoformRefs", [])),
                    "metrics": {
                        "target": metrics.get("target"),
                        "toleranceDa": metrics.get("toleranceDa"),
                        "singleModCount": len(singles),
                        "combinationCount": len(combos),
                        "singleModCandidates": singles,
                        "combinationCandidates": combos,
                    },
                },
            })
    for row in rows:
        try:
            carriers = json.loads(row.get("carrier_facts") or "[]")
        except (TypeError, ValueError):
            carriers = []
        for carrier in carriers:
            refs = carrier.get("proteoformRefs") or []
            for proteoform_ref in refs:
                protein = ref_accession(proteoform_ref)
                prsm_facts = [
                    fact for fact in (carrier.get("prsmFacts") or [])
                    if fact.get("proteoformRef") == proteoform_ref
                ]
                prsm_count = len(prsm_facts)
                calibrated = row.get("calibration_status") == "calibrated"
                observed_shift_value = float(row.get("shift_mass", 0.0))
                if observed_shift_value < -5.0 or row.get("seq_gap"):
                    state = "sequence_event_review"
                    reason = "negative mass difference or sequence gap; review truncation before modification search"
                elif prsm_facts:
                    state = "prsm_refine"
                    reason = "PrSM-backed unknown shift; construct and compare complete hypotheses"
                else:
                    state = "composition"
                    reason = "composition/candidate review before spectra"
                if state == "prsm_refine":
                    tier = 1
                elif state == "composition":
                    tier = 3
                elif state == "sequence_event_review":
                    tier = 4
                else:
                    tier = 5
                key = (row["assignment_id"], str(carrier.get("featureId", "")),
                       protein, proteoform_ref)
                if key in seen:
                    continue
                seen.add(key)
                records.append({
                    "assignmentId": row["assignment_id"],
                    "massShiftId": row["mass_shift_id"],
                    "featureId": carrier.get("featureId"),
                    "proteoformRef": proteoform_ref,
                    "protein": protein,
                    "observedShift": {
                        "value": observed_shift_value,
                        "valueBasis": "raw",
                        "frame": "unexpected_shift_after_fixed_forms",
                        "source": "TopPIC Unexpected modifications",
                        "rawValues": carrier.get("rawObservedShifts", []),
                        "calibratedValues": [
                            carrier["calibration"]["calibratedShift"]
                        ] if carrier.get("calibration", {}).get("status") == "calibrated" else [],
                        "fixedFormsExcluded": carrier.get("fixedFormsExcluded", []),
                        "terminalForms": sorted({
                            fact.get("proteinNTerminalForm") or "NONE"
                            for fact in prsm_facts
                        }),
                    },
                    "carrier": {
                        "featureMass": carrier.get("featureMass"),
                        "retentionTime": carrier.get("retentionTime"),
                        "prsmFacts": prsm_facts,
                    },
                    "prsmIds": sorted({fact.get("prsmId") for fact in prsm_facts
                                       if fact.get("prsmId")}),
                    "prsmCount": prsm_count,
                    "calibrationStatus": row.get("calibration_status", "unavailable"),
                    "measurementAdequate": calibrated,
                    "supportClass": "strong" if prsm_count >= 3 else "weak",
                    "scans": sorted({fact.get("scan") for fact in prsm_facts
                                     if fact.get("scan")}),
                    "smithBestLevel": row.get("smith_best_level")
                                      if row.get("smith_best_level") != "" else None,
                    "smithWorstLevel": row.get("smith_worst_level")
                                       if row.get("smith_worst_level") != "" else None,
                    "smithLevelCounts": row.get("smith_level_counts") or "",
                    "smithAmbiguities": row.get("smith_ambiguities") or "",
                    "topPicIntervals": carrier.get("localizations", []),
                    "routingState": state,
                    "attentionTier": tier,
                    "knowledgeState": "not_checked",
                    "knowledgeRefs": [],
                    "evidenceRefs": [],
                    "familyRefs": sorted(direct_families[str(carrier.get("featureId", ""))]),
                    "relatedFamilyRefs": sorted(
                        protein_families[protein]
                        - direct_families[str(carrier.get("featureId", ""))]
                    ),
                    "routingReason": reason,
                    "source": "TopPIC carrier_facts + same-protein recurring-shift family",
                })
    return records

def terminal_frame_evidence(rows, ppm_tol=10.0):
    """Group alternative TopPIC terminal frames by measured physical mass.

    Durable-projection schema (rows keyed for tdp-check evidence_refs);
    census_inspect.terminal_frame_equivalences is the read-side inspect
    view with its own schema and anchor handling — same physics, two
    contracts; keep them separate.

    This records representation relationships only. It deliberately does not
    name a removed residue or assign chemistry.
    """
    facts = []
    for row in rows:
        for carrier in json.loads(row.get("carrier_facts") or "[]"):
            for fact in carrier.get("prsmFacts") or []:
                facts.append((row, fact))
    grouped = []
    for row, fact in sorted(facts, key=lambda item: float(item[1]["proteoformMass"])):
        protein = ref_accession(fact.get("proteoformRef", ""))
        mass = float(fact.get("proteoformMass") or 0.0)
        match = next((group for group in grouped
                      if group["protein"] == protein
                      and abs(group["mass"] - mass) <= max(group["mass"], mass, 1.0) * ppm_tol / 1e6), None)
        if match is None:
            match = {"protein": protein, "mass": mass, "items": []}
            grouped.append(match)
        match["items"].append((row, fact))

    evidence = []
    for group in grouped:
        by_frame = defaultdict(list)
        for row, fact in group["items"]:
            mods = fact.get("unexpectedModifications") or []
            if not mods:
                continue
            key = (int(fact.get("firstResidue") or 1),
                   fact.get("proteinNTerminalForm") or "NONE")
            by_frame[key].append((row, fact))
        if len(by_frame) < 2:
            continue
        representations = []
        assignments, feature_ids, proteoform_refs, prsm_ids = set(), set(), set(), set()
        for (first, form), members in sorted(by_frame.items()):
            member_prsms = sorted({fact["prsmId"] for _, fact in members})
            shifts = [float(mod["value"])
                      for _, fact in members
                      for mod in (fact.get("unexpectedModifications") or [])]
            # A frame median mixes every modification carried by its member
            # PrSMs; when the frame's values span more than measurement noise
            # the number is a mixed population, and the frame comparison must
            # say so instead of presenting it as one chemistry.
            mixed = bool(shifts) and (max(shifts) - min(shifts)) > 0.05
            representations.append({
                "firstResidue": first,
                "nTerminalForm": form,
                "observedShift": {"value": round(statistics.median(shifts), 6),
                                  "frame": "unexpected_shift_after_fixed_forms",
                                  **({"mixedPopulation": True,
                                      "shiftRangeDa": [round(min(shifts), 6),
                                                       round(max(shifts), 6)],
                                      "note": "frame median spans >0.05 Da of "
                                              "distinct modification values"}
                                     if mixed else {})},
                "prsmIds": member_prsms,
            })
            assignments.update(row["assignment_id"] for row, _ in members)
            feature_ids.update(str(fact["featureId"]) for _, fact in members)
            proteoform_refs.update(fact["proteoformRef"] for _, fact in members)
            prsm_ids.update(member_prsms)
        identity = json.dumps([group["protein"], round(group["mass"], 4), representations],
                              sort_keys=True, separators=(",", ":"))
        evidence_id = "terminal-frame:" + hashlib.sha1(identity.encode("utf-8")).hexdigest()[:16]
        evidence.append({
            "evidenceId": evidence_id,
            "type": "terminal_frame_equivalence",
            "family": {"familyId": evidence_id.replace("terminal-frame:", "TFAM"),
                       "familyType": "terminal_frame_family"},
            "source": "TopPIC PrSM",
            "sourceRefs": sorted(prsm_ids),
            "assignmentIds": sorted(assignments),
            "featureIds": sorted(feature_ids),
            "proteoformRefs": sorted(proteoform_refs),
            "measurementRole": "representation_relationship",
            "notIdentity": True,
            "metrics": {"protein": group["protein"],
                        "proteoformMass": round(group["mass"], 6),
                        "representations": representations,
                        "shiftDifferencesDa": sorted({
                            round(abs(left["observedShift"]["value"]
                                      - right["observedShift"]["value"]), 6)
                            for index, left in enumerate(representations)
                            for right in representations[index + 1:]
                        })},
        })
    return evidence


def build_evidence_records(rows, source_parameters=None):
    """Pipeline evidence records, in memory — the preprocessing's own
    intermediates feeding candidate generation. The source-linked observation
    handoff is written by ingest.py; durable adjudication evidence lives in the
    TS PROV store."""
    records = {}
    for row in rows:
        if not row.get("prsm_count"):
            continue
        records[f"{row['assignment_id']}:toppic_prsm"] = {
            "evidenceId": f"{row['assignment_id']}:toppic_prsm",
            "assignmentId": row["assignment_id"],
            "proteoformRefs": ordered_values(row.get("proteoform_refs", "")),
            "massShiftId": row["mass_shift_id"],
            "featureIds": ordered_values(row.get("feature_ids", "")),
            "type": "toppic_prsm", "source": "TopPIC/pipeline",
            "sourceRefs": [],
            "measurementRole": "observation",
            "notIdentity": True,
            "metrics": {"prsmCount": int(row.get("prsm_count") or 0),
                        "smith": row.get("smith_level_counts", "")},
        }
    for record in terminal_frame_evidence(rows):
        records[record["evidenceId"]] = record
    if source_parameters:
        source = source_parameters.get("source", "TopPIC")
        digest = hashlib.sha1(json.dumps(source_parameters, sort_keys=True).encode(
            "utf-8")).hexdigest()[:16]
        records[f"toppic-search:{digest}"] = {
            "evidenceId": f"toppic-search:{digest}",
            "type": "toppic_search_parameters",
            "source": source,
            "sourceRefs": source_parameters.get("sourceRefs", []),
            "measurementRole": "source_search_contract",
            "notIdentity": True,
            "metrics": source_parameters.get("parameters", {}),
        }
    return list(records.values())


def write_unimod_candidate_evidence(checklist, families, prior_records, binary,
                                     tolerance_da=0.02, max_components=2):
    """Generate first-pass UniMod candidates; never assign identity or
    chemistry state. Takes the pipeline's in-memory evidence records and
    RETURNS the candidate sets — they ride the queue artifact into the TS
    PROV store; no ledger file is written."""
    prior_records = list(prior_records or [])
    targets = []
    for row in checklist:
        calibrated = [number(value, None) for value in
                      ordered_values(row.get("calibrated_shift", ""))]
        calibrated = [value for value in calibrated if value is not None]
        raw_mass = float(row["shift_mass"])
        targets.append({
            "kind": "assignment", "ref": row["assignment_id"],
            "mass": statistics.median(calibrated) if calibrated else raw_mass,
            "rawMass": raw_mass,
            "calibrationStatus": row.get("calibration_status") or "unavailable",
            "assignmentIds": [row["assignment_id"]],
            "featureIds": ordered_values(row.get("feature_ids", "")),
            "proteoformRefs": ordered_values(row.get("proteoform_refs", "")),
        })
    for family in families:
        if family.get("spacing_center_da") is None or not family.get("family_id"):
            continue
        targets.append({
            "kind": "family_spacing", "ref": family["family_id"],
            "mass": float(family["spacing_center_da"]),
            "familyRefs": [family["family_id"]],
            "featureIds": [str(value) for value in family.get("feature_ids", [])],
            "spread": float(family.get("spacing_spread_da") or 0.0),
        })
    for evidence in prior_records:
        if evidence.get("type") != "terminal_frame_equivalence":
            continue
        offsets = []
        for assignment_id in evidence.get("assignmentIds", []):
            row = next((item for item in checklist
                        if item["assignment_id"] == assignment_id), None)
            if row is None:
                continue
            for carrier in json.loads(row.get("carrier_facts") or "[]"):
                calibration = carrier.get("calibration") or {}
                if calibration.get("status") == "calibrated":
                    offsets.append(float(calibration.get("axisShiftDa") or 0.0))
        calibration_offset = statistics.median(offsets) if offsets else 0.0
        for index, representation in enumerate(
                (evidence.get("metrics") or {}).get("representations") or []):
            observed = representation.get("observedShift") or {}
            if observed.get("value") is None:
                continue
            raw_mass = float(observed["value"])
            targets.append({
                "kind": "terminal_frame",
                "ref": f"{evidence['evidenceId']}:{index}",
                "mass": raw_mass + calibration_offset,
                "rawMass": raw_mass,
                "calibrationStatus": "calibrated" if offsets else "unavailable",
                "assignmentIds": evidence.get("assignmentIds", []),
                "featureIds": evidence.get("featureIds", []),
                "proteoformRefs": evidence.get("proteoformRefs", []),
                "sourceRefs": [evidence["evidenceId"]],
                "representation": representation,
            })
    if not targets:
        return []

    queries = []
    query_ids = {}
    for target in targets:
        tolerance = max(tolerance_da, target.get("spread", 0.0) / 2.0)
        key = (round(target["mass"], 6), round(tolerance, 6))
        if key not in query_ids:
            query_ids[key] = f"Q{len(query_ids) + 1}"
            queries.append({"id": query_ids[key], "mass": target["mass"],
                            "tolerance": tolerance})
        target["queryId"] = query_ids[key]
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8",
                                     suffix=".json", delete=False) as fh:
        json.dump(queries, fh, separators=(",", ":"))
        query_path = Path(fh.name)
    try:
        search_binary = native.binary(binary)
        solve_binary = native.binary(binary)
        search_doc = native.run_json(
            [search_binary, "unimod-batch-search", "--input", str(query_path),
             "--limit", "20"], label="tdeasenexus-toolbox unimod-batch-search")
        solve_doc = native.run_json(
            [solve_binary, "unimod-batch-solve", "--input", str(query_path),
             "--max-modifications", str(max_components), "--limit", "20"],
            label="tdeasenexus-toolbox unimod-batch-solve")
    finally:
        query_path.unlink(missing_ok=True)

    searches = {item["id"]: item for item in
                (search_doc.get("result") or {}).get("queries", [])}
    solutions = {item["id"]: item for item in
                 (solve_doc.get("result") or {}).get("queries", [])}
    records = {
        record["evidenceId"]: record for record in prior_records
        if (record.get("evidenceId")
            and record.get("type") != "unimod_candidate_set")
    }
    tolerance_by_id = {item["id"]: item["tolerance"] for item in queries}
    for target in targets:
        query_id = target["queryId"]
        evidence_id = f"unimod-candidates:{target['kind']}:{target['ref']}"
        records[evidence_id] = {
            "evidenceId": evidence_id,
            "type": "unimod_candidate_set",
            "source": "tdeasenexus-toolbox UniMod kernel",
            "sourceRefs": target.get("sourceRefs", [target["ref"]]),
            "assignmentIds": target.get("assignmentIds", []),
            "familyRefs": target.get("familyRefs", []),
            "featureIds": target.get("featureIds", []),
            "proteoformRefs": target.get("proteoformRefs", []),
            "measurementRole": "candidate_generation",
            "notIdentity": True,
            "metrics": {
                "target": {
                    "value": target["mass"],
                    "rawValue": target.get("rawMass", target["mass"]),
                    "calibrationStatus": target.get("calibrationStatus", "not_applicable"),
                    "frame": ("unexpected_shift_after_fixed_forms"
                              if target["kind"] in {"assignment", "terminal_frame"}
                              else "measured_family_spacing"),
                },
                "terminalRepresentation": target.get("representation"),
                "toleranceDa": tolerance_by_id[query_id],
                "singleModCandidates": searches.get(query_id, {}).get("massMatches", []),
                "combinationCandidates": (solutions.get(query_id, {}).get("solution") or {}).get(
                    "candidates", []),
            },
        }
    return [records[evidence_id] for evidence_id in sorted(records)
            if records[evidence_id].get("type") == "unimod_candidate_set"]


def calibrate_pending(checklist, mzml, binary, ppm_tol):
    spectrum_binary = native.binary(binary)
    # A carrier-level correction is propagated to its source rows. It is only
    # defined when those rows share one baseline and one residual component.
    # Never cross-product independently sorted sequences and residuals.
    bases_by_feature = defaultdict(set)
    multi_component_features = set()
    for row in checklist:
        for carrier in json.loads(row["carrier_facts"]):
            feature_id = str(carrier["featureId"])
            bases_by_feature[feature_id].update(
                carrier.get("sequencesWithForms") or carrier.get("sequences", []))
            if any(len(fact.get("unexpectedModifications", [])) > 1
                   for fact in carrier.get("prsmFacts", [])):
                multi_component_features.add(feature_id)
    scans = sorted({
        scan.removeprefix("scan=")
        for row in checklist
        for carrier in json.loads(row["carrier_facts"])
        for scan in carrier.get("scans", [])
    }, key=lambda value: int(value))
    if not scans:
        return
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", suffix=".json", delete=False) as fh:
        json.dump(scans, fh, separators=(",", ":"))
        scans_path = fh.name
    try:
        metadata = native.run_json([
            spectrum_binary, "spectrum", mzml, "--mode", "scans", "--scans-input", scans_path,
            "--max-peaks", "1",
        ], "tdeasenexus-toolbox spectrum")
    finally:
        Path(scans_path).unlink(missing_ok=True)
    spectra = metadata.get("result", {}).get("spectra", [])
    by_scan = {item["scan"]: item["spectrum"] for item in spectra}
    requests = []
    carriers_by_id = {}
    for row in checklist:
        carriers = json.loads(row["carrier_facts"])
        for carrier in carriers:
            carrier["calibration"] = {"status": "unavailable"}
            feature_id = str(carrier["featureId"])
            if (len(bases_by_feature[feature_id]) != 1
                    or feature_id in multi_component_features):
                carrier["calibration"]["reason"] = (
                    "Carrier-level residual calibration requires one declared baseline "
                    "and one residual component; refine the source models separately."
                )
                continue
            scan = next((value.removeprefix("scan=") for value in carrier.get("scans", [])
                         if by_scan.get(value.removeprefix("scan="), {}).get("precursor", {}).get("charge")), None)
            spectrum = by_scan.get(scan) if scan else None
            precursor = spectrum.get("precursor") if spectrum else None
            shifts = carrier.get("rawObservedShifts", [])
            sequences = carrier.get("sequencesWithForms") or carrier.get("sequences", [])
            if not precursor or not precursor.get("charge") or not shifts or not sequences:
                continue
            request_id = f'{row["assignment_id"]}:{carrier["featureId"]}'
            shift = float(shifts[0])
            ms2_rt = float(spectrum["retentionTime"])
            rt_apex = carrier.get("retentionTime", {}).get("apex")
            # Center the MS1 window on the feature apex nearest the MS2 scan,
            # not on the single MS2 scan's RT — a fixed ±3 s window catches
            # only one MS1 scan and produces single-scan calibration noise.
            rt = float(rt_apex) if rt_apex is not None else ms2_rt
            request = {
                "id": request_id,
                "sequence": sequences[0],
                "mz": float(precursor["selectedMz"]),
                "charge": int(precursor["charge"]),
                "rtMin": rt - CALIBRATION_RT_HALF_WIDTH,
                "rtMax": rt + CALIBRATION_RT_HALF_WIDTH,
                "tolerancePpm": ppm_tol,
                "priorShift": shift,
                "matchWindowPpm": 30.0,
                "modificationMass": shift,
            }
            requests.append(request)
            carriers_by_id[request_id] = carrier
        row["carrier_facts"] = json.dumps(carriers, separators=(",", ":"))


    if not requests:
        return
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", suffix=".json", delete=False) as fh:
        json.dump(requests, fh, separators=(",", ":"))
        request_path = fh.name
    try:
        batch = native.run_json([spectrum_binary, "refine-shift-batch", mzml, "--input", request_path],
                         "tdeasenexus-toolbox refine-shift-batch")
    finally:
        Path(request_path).unlink(missing_ok=True)
    for item in batch.get("result", {}).get("results", []):
        carrier = carriers_by_id.get(item.get("id"))
        if not carrier:
            continue
        if not item.get("ok"):
            carrier["calibration"] = {"status": "failed", "error": item.get("error", {})}
            continue
        result = item["result"]
        ambiguity = result.get("ambiguity")
        pearson = result.get("pearson")
        matched = int(result["matchedPeaks"])
        scans = int(result["scanCount"])
        quality_ok = (
            isinstance(pearson, (int, float))
            and pearson >= CALIBRATION_MIN_PEARSON
            and matched >= CALIBRATION_MIN_MATCHED_PEAKS
            and scans >= CALIBRATION_MIN_SCAN_COUNT
        )
        status = "ambiguous" if ambiguity else ("calibrated" if quality_ok else "weak")
        selected_evaluation = next(
            (evaluation for evaluation in result.get("evaluations", [])
             if int(evaluation.get("offset", 0)) == int(result["offset"])),
            {},
        )
        carrier["calibration"] = {
            "status": status,
            # Post-declared-forms frame (TopPIC's unexpected-shift frame):
            # refinedShiftAfterDeclared is measured after the record's named
            # fixed forms; refinedShift (bare frame) would silently absorb
            # the terminal background. declaredDeltaMass records what the
            # forms contributed, tool-derived.
            "calibratedShift": round(float(
                result.get("refinedShiftAfterDeclared", result["refinedShift"])), 6),
            "declaredDeltaMass": result.get("declaredDeltaMass"),
            "isotopeOffset": int(result["offset"]),
            "axisShiftDa": round(float(selected_evaluation.get("fittedDa", 0.0)), 6),
            "matchedPeaks": matched,
            "pearson": pearson,
            "scanCount": scans,
            "ambiguity": ambiguity,
            "quality": None if status == "calibrated" else (
                f"{status}: pearson={pearson if isinstance(pearson, (int, float)) else 'n/a'} "
                f"peaks={matched} scans={scans}"
            ),
        }
    # Copy the batch calibration back into the pending rows' carrier facts:
    # recalibrate_observation_clusters reads the carrier-level fit.
    for row in checklist:
        carriers = json.loads(row["carrier_facts"])
        for carrier in carriers:
            request_id = f'{row["assignment_id"]}:{carrier["featureId"]}'
            if request_id in carriers_by_id:
                carrier["calibration"] = carriers_by_id[request_id]["calibration"]
        row["carrier_facts"] = json.dumps(carriers, separators=(",", ":"))

def recalibrate_observation_clusters(observations, preliminary):
    accepted, calibration = {}, {}
    for row in preliminary:
        for carrier in json.loads(row.get("carrier_facts") or "[]"):
            feature_id = str(carrier.get("featureId", ""))
            fit = carrier.get("calibration") or {}
            shifts = carrier.get("rawObservedShifts") or []
            if not shifts:
                continue
            # One MS1 feature can carry several single-shift proteoform
            # interpretations (e.g. +42 and +79 readings of one envelope).
            # A fit was measured on ONE of them: keying by feature alone lets
            # the last interpretation's value silently overwrite every other
            # interpretation's observation rows.
            key = (feature_id, round(float(shifts[0]), 6))
            calibration[key] = fit
            if fit.get("status") == "calibrated":
                accepted[key] = float(fit["calibratedShift"])

    def accepted_for(feature_id, raw_shift):
        exact = accepted.get((feature_id, round(float(raw_shift), 6)))
        if exact is not None:
            return exact
        near = [(abs(key[1] - float(raw_shift)), key) for key in accepted
                if key[0] == feature_id and abs(key[1] - float(raw_shift)) <= 0.05]
        return accepted[min(near)[1]] if near else None

    original_values = defaultdict(list)
    for row in observations:
        value = accepted_for(str(row.get("feature_id", "")),
                             row.get("raw_shift_mass") or row["shift_mass"])
        if value is not None:
            original_values[row["mass_shift_id"]].append(value)
    original_centers = {key: statistics.median(values) for key, values in original_values.items()}
    shifted_values = []
    for row in observations:
        value = accepted_for(str(row.get("feature_id", "")),
                             row.get("raw_shift_mass") or row["shift_mass"])
        if value is not None and abs(value - float(row["shift_mass_center"])) > 0.05:
            shifted_values.append(value)
    groups = []
    for value in sorted(set(shifted_values)):
        if not groups or value - groups[-1][-1] > 0.05:
            groups.append([value])
        else:
            groups[-1].append(value)
    center_for = {value: statistics.median(group) for group in groups for value in group}
    for row in observations:
        value = accepted_for(str(row.get("feature_id", "")),
                             row.get("raw_shift_mass") or row["shift_mass"])
        if value is not None:
            original = float(row["shift_mass_center"])
            if abs(value - original) <= 0.05:
                center = original_centers[row["mass_shift_id"]]
                cluster_id = row["mass_shift_id"]
            else:
                center = center_for[value]
                cluster_id = f"MCAL{center:+.6f}"
            row.update(shift_mass=f"{value:.6f}", shift_mass_center=f"{center:.6f}",
                       mass_shift_id=cluster_id)
    return calibration

def attach_calibration(checklist, calibration):
    for row in checklist:
        carriers = json.loads(row.get("carrier_facts") or "[]")
        for carrier in carriers:
            shifts = carrier.get("rawObservedShifts") or []
            key = ((str(carrier.get("featureId", "")), round(float(shifts[0]), 6))
                   if shifts else None)
            carrier["calibration"] = calibration.get(
                key, {"status": "unavailable"})
        statuses = {c.get("calibration", {}).get("status", "unavailable") for c in carriers}
        settled = [c["calibration"] for c in carriers if c.get("calibration", {}).get("status") == "calibrated"]
        row["calibration_status"] = "calibrated" if statuses == {"calibrated"} else "; ".join(sorted(statuses))
        row["calibrated_shift"] = "; ".join(str(x["calibratedShift"]) for x in settled)
        row["isotope_offset"] = "; ".join(str(x["isotopeOffset"]) for x in settled)
        row["carrier_facts"] = json.dumps(carriers, separators=(",", ":"))
