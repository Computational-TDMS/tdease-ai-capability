"""
Ingest command: TopPIC tables -> typed feature/proteoform/PrSM facts,
relationships, and bounded durable projections. No chemistry assignment
occurs here. Orchestration only — the machinery lives in ingest_families
(ladders), ingest_tables (rows/checklist), and ingest_evidence (queue/
evidence/calibration).
"""

import argparse
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

from toppic_io import read_toppic_parameters
from ingest_tables import (
    aggregate_checklist,
    build_single_sample,
    split_values,
)

from ingest_evidence import (
    attach_calibration,
    build_evidence_records,
    calibrate_pending,
    recalibrate_observation_clusters,
    residual_queue_records,
    write_unimod_candidate_evidence,
)

BUNDLE_FILE = "opensearch.json"


def build_annotation_projection(candidate_records):
    """Project automatic UniMod attribution into the shared display contract.

    The labels are generated automatically from the measured mass and are
    therefore an explicit first-pass attribution. They remain separate from
    source observations and carry ``notIdentity`` until the AI/deep-review
    stage supplies a final adjudication or manual overlay.
    """
    annotations = []
    for record in candidate_records or []:
        metrics = record.get("metrics") or {}
        target = metrics.get("target") or {}
        mass = target.get("value")
        if mass is None:
            continue
        candidates = metrics.get("singleModCandidates") or []
        labels = []
        for candidate in candidates:
            if not isinstance(candidate, dict):
                continue
            label = (
                candidate.get("name")
                or candidate.get("title")
                or candidate.get("label")
                or candidate.get("fullName")
            )
            if label and label not in labels:
                labels.append(str(label))
        if not labels:
            continue
        if len(labels) == 1:
            label = labels[0]
            basis = "automatic approved UniMod mass attribution within tolerance"
        else:
            label = " / ".join(labels[:4])
            if len(labels) > 4:
                label += " / ..."
            basis = f"automatic UniMod attribution: {len(labels)} candidates remain within tolerance"
        annotations.append({
            "mass": float(mass),
            "label": label,
            "role": "automatic",
            "stage": "automatic_unimod",
            "basis": basis,
            "sources": [record.get("evidenceId", "")],
            "assignmentIds": sorted(record.get("assignmentIds") or []),
            "familyRefs": target.get("familyRefs", []),
            "featureIds": target.get("featureIds", []),
            "proteoformRefs": target.get("proteoformRefs", []),
            "frame": target.get("frame", "unexpected_shift_after_fixed_forms"),
            "verdict": "automatic",
            "notIdentity": True,
        })
    return annotations


def observation_records(rows):
    """Build compact source-linked assignment records for census/Copilot/HTML.

    The raw TopPIC tables remain authoritative.  This projection only joins
    their IDs and preserves the mass frame, TopPIC background and nested PrSM
    facts so downstream code never has to reconstruct the chain.
    """
    records = []
    for row in rows:
        try:
            carriers = json.loads(row.get("carrier_facts") or "[]")
        except (TypeError, ValueError):
            carriers = []
        raw_values = sorted({
            value
            for carrier in carriers
            for value in (carrier.get("rawObservedShifts") or [])
        })
        records.append({
            "schema": "tdease.observation.v1",
            "assignmentId": row.get("assignment_id", ""),
            "massShift": {
                "value": row.get("shift_mass", ""),
                "frame": "unexpected_shift_after_fixed_forms",
                "rawValues": raw_values or [row.get("shift_mass", "")],
            },
            "source": {
                "ms1": row.get("source_ms1", ""),
                "proteoforms": row.get("source_proteoforms", ""),
                "prsms": row.get("source_prsms", ""),
            },
            "ids": {
                "feature": sorted(split_values(row.get("feature_ids", ""))),
                "proteoform": sorted(split_values(row.get("proteoform_refs", ""))),
                "prsm": sorted({p for c in carriers for p in c.get("prsmIds", [])}),
            },
            "smith": {
                "best": row.get("smith_best_level", ""),
                "worst": row.get("smith_worst_level", ""),
                "ambiguities": sorted(split_values(row.get("smith_ambiguities", ""))),
            },
            "topPic": {
                "status": row.get("status", ""),
                "ms2Support": row.get("ms2_support", ""),
                "families": sorted(split_values(row.get("family_refs", ""))),
                "unimodProposals": sorted(split_values(row.get("unimod_refs", ""))),
                "carriers": carriers,
            },
        })
    return records


def _compact_census(rows):
    grouped = defaultdict(lambda: {"count": 0, "assignmentIds": []})
    for row in rows:
        item = grouped[row["mass_shift_id"]]
        item["count"] += int(row.get("prsm_count") or 0)
        item["assignmentIds"].append(row["assignment_id"])
        item["mass"] = float(row["shift_mass"])
    return {
        "clusters": [
            {"id": mass_shift_id, "center_mass": item["mass"],
             "count": item["count"],
             "assignmentIds": sorted(item["assignmentIds"])}
            for mass_shift_id, item in sorted(
                grouped.items(), key=lambda pair: (pair[1]["mass"], pair[0]))
        ]
    }


def _compact_candidate_set(record):
    """Keep enough UniMod detail for routine assignment without a dump."""
    metrics = record.get("metrics") or {}
    candidates = []
    for key in ("singleModCandidates", "combinationCandidates"):
        for item in metrics.get(key) or []:
            if not isinstance(item, dict):
                continue
            summary = {
                name: item[name]
                for name in ("name", "title", "fullName", "deltaMass",
                             "residues", "location", "positions")
                if item.get(name) is not None
            }
            if summary and summary not in candidates:
                candidates.append(summary)
            if len(candidates) >= 3:
                break
        if len(candidates) >= 3:
            break
    return {
        "evidenceId": record.get("evidenceId", ""),
        "kind": record.get("kind") or record.get("type") or "assignment",
        "measurementRole": "candidate_generation",
        "notIdentity": True,
        "question": (record.get("assignmentIds") or [record.get("question")])[0],
        "familyRefs": record.get("familyRefs") or (
            [record["familyId"]] if record.get("familyId") else []
        ),
        "featureIds": record.get("featureIds", []),
        "proteoformRefs": record.get("proteoformRefs", []),
        "target": metrics.get("target"),
        "singleModCount": len(metrics.get("singleModCandidates") or []),
        "combinationCount": len(metrics.get("combinationCandidates") or []),
        "topCandidates": candidates,
    }


def write_compact_bundle(rows, families, queue_records, annotations, library,
                         args, path, join_summary=None,
                         calibration_outcome=None, candidate_outcome=None):
    """Write the single compact batch ledger used by the default workflow.

    Assignment observations carry their routing record inline. Family and
    candidate records remain typed collections in the same file, so consumers
    do not need to join three duplicate assignment projections.
    """
    observations = observation_records(rows)
    routes = defaultdict(list)
    for record in queue_records:
        if record.get("assignmentId"):
            routes[record["assignmentId"]].append(record)
    for observation in observations:
        route_records = routes.get(observation["assignmentId"], [])
        if route_records:
            # The observation already owns the full mass/carrier evidence.
            # Keep only routing fields here; copying the queue's carrier facts
            # would recreate the large queue inside the bundle.
            observation["routing"] = [{
                key: route.get(key)
                for key in (
                    "routingState", "attentionTier", "routingReason",
                    "massShiftId", "featureId", "proteoformRef", "protein",
                    "familyRefs", "relatedFamilyRefs", "prsmIds", "scans",
                    "knowledgeState", "knowledgeRefs", "evidenceRefs",
                )
                if route.get(key) not in (None, "", [], {})
            } for route in route_records]
    candidate_sets = [record["candidateSet"] for record in queue_records
                      if record.get("candidateSet")]
    if not args.include_unimod_candidates:
        candidate_sets = [_compact_candidate_set(record) for record in candidate_sets]
    payload = {
        "schema": "tdease.opensearch.v1",
        "artifactPolicy": "compact",
        "sources": {
            "ms1": str(Path(args.ms1).resolve()),
            "proteoforms": str(Path(args.proteoforms or args.prsms).resolve()),
            "prsms": str(Path(args.prsms).resolve()),
            "mzml": str(Path(args.mzml).resolve()),
        },
        "parameters": {
            "qvalueMax": args.qvalue_max,
            "evalueFallback": args.evalue_fallback,
            "ppmTolerance": args.ppm_tol,
            "rtGap": args.rt_gap,
            "nMaxC13": args.n_max_c13,
            "calibration": calibration_outcome
            or ("skipped" if args.skip_calibration else "no_carrier_calibrated"),
            **({"proteoformTable": "derived-from-prsms"}
               if args.proteoforms is None else {}),
        },
        "summary": {
            "observations": len(observations),
            "massShifts": len({row["mass_shift_id"] for row in rows}),
            "families": len(families),
            "candidateSets": len(candidate_sets),
            **({"sourceJoin": join_summary} if join_summary else {}),
        },
        "cost": {
            "finalArtifacts": 1,
            "candidatePass": candidate_outcome
            or ("skipped" if args.skip_unimod_candidates else "searched_no_hits"),
            "routing": "embedded",
            "duplicateProjections": 0,
        },
        "census": _compact_census(rows),
        "observations": observations,
        "families": families,
        "candidateSets": candidate_sets,
        "annotations": annotations,
        "modificationLibrary": library,
    }
    Path(path).write_text(
        json.dumps(payload, ensure_ascii=False, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )


def build_modification_library(candidate_records):
    """Project first-pass single-mod UniMod hits into resolve-shift library input.

    This is a search library, not a verdict: only records returned by the
    Rust UniMod search are copied, and positional constraints are retained
    when the source result provides them.
    """
    entries = {}
    for record in candidate_records or []:
        metrics = record.get("metrics") or {}
        for item in metrics.get("singleModCandidates") or []:
            if not isinstance(item, dict):
                continue
            name = item.get("name") or item.get("title") or item.get("fullName")
            if not name:
                continue
            key = str(item.get("unimodId") or name)
            entry = {"name": str(name)}
            if item.get("deltaMass") is not None:
                entry["deltaMass"] = float(item["deltaMass"])
            if item.get("residues"):
                entry["residues"] = list(item["residues"])
            if item.get("location"):
                entry["location"] = item["location"]
            if item.get("positions"):
                entry["positions"] = list(item["positions"])
            entries[key] = entry
    return [entries[key] for key in sorted(entries)]

def build(args):
    output_dir = Path(args.out_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    features, modifications, families, join_summary = build_single_sample(
        args.ms1, args.proteoforms, args.prsms,
        qvalue_max=args.qvalue_max,
        evalue_fallback=args.evalue_fallback,
        ppm_tol=args.ppm_tol,
        rt_gap_limit=args.rt_gap,
        n_max_c13=args.n_max_c13,
    )

    if modifications:
        binary = args.tdeasenexus_toolbox
        preliminary = aggregate_checklist(modifications)
        calibration = {}
        if not args.skip_calibration:
            calibrate_pending(preliminary, args.mzml, binary, args.ppm_tol)
            calibration = recalibrate_observation_clusters(modifications, preliminary)
        checklist = aggregate_checklist(modifications)
        if not args.skip_calibration:
            attach_calibration(checklist, calibration)
    else:
        checklist = []
    # Pipeline evidence stays in memory: it feeds candidate generation, and
    # durable evidence lives in the TS PROV store — no ledger file.
    evidence_records = build_evidence_records(checklist, {
        "source": "TopPIC PrSM table preamble",
        "sourceRefs": [str(Path(args.prsms).resolve())],
        "parameters": read_toppic_parameters(args.prsms),
    })
    candidates = []
    if not args.skip_unimod_candidates:
        candidates = write_unimod_candidate_evidence(
            checklist, families, evidence_records, args.tdeasenexus_toolbox,
            args.candidate_tolerance_da, args.candidate_max_components)
    candidate_by_assignment = defaultdict(list)
    for candidate in candidates:
        for assignment_id in candidate.get("assignmentIds", []):
            candidate_by_assignment[assignment_id].append(candidate.get("evidenceId", ""))
    for row in checklist:
        refs = sorted(set(candidate_by_assignment.get(row["assignment_id"], [])))
        row["unimod_status"] = "available" if refs else "none"
        row["unimod_refs"] = "; ".join(refs)
        row_features = set(split_values(row.get("feature_ids", "")))
        row["family_refs"] = "; ".join(sorted({
            family.get("family_id", "") for family in families
            if row_features.intersection(str(value) for value in family.get("feature_ids", []))
        }))
        row["source_ms1"] = str(Path(args.ms1).resolve())
        row["source_proteoforms"] = str(Path(args.proteoforms or args.prsms).resolve())
        row["source_prsms"] = str(Path(args.prsms).resolve())
    annotations = build_annotation_projection(candidates)
    library = build_modification_library(candidates)
    bundle_path = output_dir / BUNDLE_FILE
    queue_records = residual_queue_records(checklist, families, candidates)
    # Report outcomes, not intents: "calibration: true" when every carrier
    # failed would claim a correction that never happened; a candidate search
    # that ran and found nothing is not "skipped".
    calibration_outcome = (
        "skipped" if args.skip_calibration else
        "applied" if checklist and any(
            row.get("calibration_status") == "calibrated" for row in checklist)
        else "no_carrier_calibrated")
    candidate_outcome = (
        "skipped" if args.skip_unimod_candidates else
        "included" if candidates else "searched_no_hits")
    write_compact_bundle(checklist, families, queue_records, annotations,
                         library, args, bundle_path, join_summary=join_summary,
                         calibration_outcome=calibration_outcome,
                         candidate_outcome=candidate_outcome)
    queue_states = Counter()
    queue_tiers = Counter()
    for record in queue_records:
        if "familyId" in record or "candidateSet" in record:
            continue  # family and candidate records are not assignment rows
        queue_states[record.get("routingState", "unknown")] += 1
        queue_tiers[str(record.get("attentionTier", "unknown"))] += 1
    queue_sha = hashlib.sha256(bundle_path.read_bytes()).hexdigest()
    print(json.dumps({
        "artifactPolicy": "compact",
        **({"proteoformTable": "derived-from-prsms"}
           if args.proteoforms is None else {}),
        "cost": {
            "finalArtifacts": 1,
            "candidatePass": candidate_outcome,
            "routing": "embedded",
        },
        "compactFile": str(bundle_path),
        "families": len(families),
        "familyTypes": dict(sorted(Counter(
            family.get("evidence_type", "ms1_rt_coherent")
            for family in families
        ).items())),
        "sourceFeatures": len({item for row in features for item in split_values(row["source_feature_ids"])}),
        "features": len(features),
        "sourceJoin": join_summary,
        "massShifts": len({row["mass_shift_id"] for row in checklist}),
        "assignmentRows": len(checklist),
        "observationRecords": len(checklist),
        "modificationLibraryEntries": len(library),
        "annotations": len(annotations),
        "residualQueueRows": sum(queue_states.values()),
        "residualQueueStates": dict(sorted(queue_states.items())),
        "residualQueueTiers": dict(sorted(queue_tiers.items())),
        "unimodCandidateSets": len(candidates),
        # The single-line handoff block the TS controller ingests by
        # reference: batch payloads ride the artifact file, never stdout.
        "censusQueue": {
            "path": str(bundle_path),
            "sha256": queue_sha,
            "rows": sum(queue_states.values()),
        },
    }, separators=(",", ":")))

def parser():
    root = argparse.ArgumentParser(description=__doc__)
    root.add_argument("--ms1", required=True,
                      help="TopPIC ms1 feature table (e.g. <run>_ms1.feature)")
    root.add_argument("--proteoforms", default=None,
                      help="TopPIC proteoform table (e.g. <run>_ms2_toppic_proteoform_"
                           "single.tsv); derived from the PrSM table when omitted")
    root.add_argument("--prsms", required=True,
                      help="TopPIC prsm table; the full <run>_ms2_toppic_prsm.tsv is "
                           "accepted as-is")
    root.add_argument("--mzml", required=True,
                      help="raw spectrum file (e.g. <run>.mzML)")
    root.add_argument("--out-dir", required=True,
                      help="chain working directory (share it with `run "
                           "--dir`; e.g. <runs>/<sample>/analysis/mass-shift-chain)")
    root.add_argument("--qvalue-max", type=float, default=0.01,
                      help="TopPIC spectrum/proteoform FDR gate")
    root.add_argument("--evalue-fallback", type=float, default=1e-5,
                      help="used only when the source table has no Q-values")
    root.add_argument("--ppm-tol", type=float, default=10.0,
                      help="ppm on full proteoform mass after integer-C13 correction")
    root.add_argument("--rt-gap", type=float, default=0.2,
                      help="maximum gap between RT intervals in minutes")
    root.add_argument("--n-max-c13", type=int, default=1)
    root.add_argument("--skip-calibration", action="store_true",
                      help="skip the embedded full-ledger MS1 calibration pass; "
                           "focused Interpret builds finish inside one bash timeout "
                           "and deep-dive tools measure envelopes on demand")
    root.add_argument("--tdeasenexus-toolbox", default=None)
    root.add_argument("--skip-unimod-candidates", action="store_true")
    root.add_argument("--candidate-tolerance-da", type=float, default=0.02)
    root.add_argument("--candidate-max-components", type=int, default=2)
    root.add_argument(
        "--include-unimod-candidates", action="store_true",
        help=("embed full UniMod candidate arrays instead of the default "
              "top-three summaries"),
    )
    return root

def build_main(argv=None):
    args = parser().parse_args(argv)
    try:
        build(args)
    except (OSError, ValueError, RuntimeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0




def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in {"-h", "--help"}:
        print("usage: main.py ingest build [args...]\n")
        print("Batch ingestion: TopPIC tables -> opensearch.json.")
        return 0 if argv else 2
    mode, rest = argv[0], argv[1:]
    if mode != "build":
        print(f"unknown ingest mode {mode!r}", file=sys.stderr)
        return 2
    return build_main(rest)
