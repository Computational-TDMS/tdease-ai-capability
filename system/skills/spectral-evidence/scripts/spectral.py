#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""One entry for spectral evidence: measure, read the card, look at the windows.

Three verbs, and every one of them answers a question you can ask in a
sentence:

  annotate  --msalign F --scan S --proforma M
        Where are the leads in the deconvoluted spectrum? One complete model
        annotated on a TopFD msalign scan: matched fragments labelled with
        their sequence, ppm error panel, sequence coverage.

  measure      --mzml F --scan S[,S..] --models M.json
        Do these models actually differ, and on which ions? Runs the kernel's
        envelope matcher once per scan for EVERY model on the same observed
        peaks, and returns an evidence card: the ions whose predicted m/z
        differs between the models, which model matched each one, and — for
        the ones nobody matched — the kernel's reason. The default output is
        the full figure set — overview figure plus the separating isotope
        windows — so nothing downstream needs rendering by hand;
        --no-figures keeps the measurement only. Add --ms1 to ask the
        same question of precursor envelopes instead of fragments.

  show      CARD.json [--ions <ion:charge,...>]
        Look at the windows themselves: each requested ion's theoretical
        isotope comb against the raw peaks, one panel per model, so a
        populated window in one model and an empty one in the other is
        visible rather than asserted.

WHAT THIS TOOL WILL NOT DO. It does not judge existence. Every match, state
and correlation in the card is emitted by the Rust kernel's isotope-envelope
matcher; this entry reports them. Writing your own monoisotopic-peak window
over a peak list to decide whether an ion is "really" there produces false
negatives (an M+2 hit is a hit) and contradicts the criterion the rest of the
analysis uses. If a match surprises you, read its state in the card and look
at the window with `show` — do not re-match by hand.

READ THE CARD THIS WAY. `verdict.compositionDeltaDa` is how far apart the
models' complete chemistry is, from the kernel's own parse of each full
model; ions whose span covers a differing site are the ones that can separate
them, and `verdict.ions` lists those the kernel measured, most separating
first (match count, then correlation and error — never the isotope slot a
match landed in, which is a fact to report, not a rank), each carrying the
kernel's own state per model. `verdict.differences` separates predicted
difference from measured difference. `uninformative` counts what the card
did not list.

The card is the whole stdout: what you see is what is written to --out, so it
can be handed on as the measurement it is.

Usage:
  spectral.py annotate --msalign <sample_ms2.msalign> --scan <scan> \\
      --proforma "<complete-model>" --out annotation.html
  spectral.py measure --mzml <source.mzML> --scan <scan[,scan...]> --models <models.json> \\
      --out <card.json>
  spectral.py measure --mzml <source.mzML> --ms1 \\
      --feature <label> <mass> <charge> <rtMin> <rtMax> --models <models.json> \\
      --out <precursor-card.json>
  spectral.py show <card.json> --ions <ion:charge[,ion:charge...]>
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import types
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _shared.chart_forge import emit_html
from _shared.kernel_io import load_doc
from _shared.manifest import write_manifest
from _shared.native import load_json_arg, run_json

import envelope_ladder
import envelope_mirror
import ms2_visualize
from _ms2render import (
    AgentArgumentParser, candidate_rows, ion_types_label, measure_candidates,
    ms_binary,
)

CODE = "SPECTRAL_FAILED"
# Bookkeeping split between "numerically identical prediction" and
# "predicted differently" — float-identity only, NOT a separability claim.
# Whether a real difference is resolvable depends on the measurement
# conditions (tolerance, charge ceiling, interference), which the card
# reports; identical m/z can still hide different envelope predictions.
MIN_DELTA_DA = 1e-8
# The card stays the size of a verdict, not of a match table — and what it
# drops is counted, never silently missing.
MAX_VERDICT_IONS = 24
MAX_SITE_ROWS = 16
MAX_CONTESTED = 8
DEFAULT_SHOW_IONS = 6
MAX_DEGENERATE_SHOWN = 3
SHOW_COLORS = envelope_mirror.MIRROR_PALETTE


# --------------------------------------------------------------- card helpers

def _key(row: dict) -> tuple:
    return (row["ion"], int(row["number"]), int(row["charge"]))


def _label(key: tuple) -> str:
    return f"{key[0]}{key[1]}:{key[2]}"


def _by_scan(values: dict, single: str, multiple: str) -> dict:
    """One scan's value inline, several keyed by scan.

    ``observation.scans`` always lists what was measured, and the two key
    names are distinct, so the reader never has to guess which shape it got."""
    if len(values) == 1:
        return {single: next(iter(values.values()))}
    return {multiple: values}


def _diff_sites(parsed: list[dict]) -> tuple[list, bool]:
    """Different modification claims on a common sequence, including terminals."""
    if any(p["sequence"] != parsed[0]["sequence"] for p in parsed[1:]):
        return [], False
    tags = [ms2_visualize.mod_tags(model) for model in parsed]
    diff = sorted({site for own in tags for site in own
                   if len({other.get(site) for other in tags}) > 1}, key=str)
    return diff, "nTerm" in diff


def _reported_sites(ion: str, number: int, length: int) -> dict:
    """The residues whose modification state moves this fragment's m/z.

    An N-terminal a/b/c ion of length n holds residues 1..n; a C-terminal
    x/y/z ion of length n holds the last n. A modification outside that span
    cannot shift the fragment — a difference moves only the ions whose span
    includes it, so an N-terminal difference reaches every a/b/c ion but only
    the x/y/z ions long enough to reach it."""
    if ion in {"a", "b", "c"}:
        return {"from": 1, "to": min(number, length)}
    return {"from": max(1, length - number + 1), "to": length}


def _cell(status: dict | None, match: dict | None) -> dict:
    """One model's answer for one ion: the kernel's own state.

    A non-match still carries the kernel's measurements — the envelope
    correlation it measured before rejecting, and which isotope slot the apex
    was expected at — because "correlationRejected" without the number invites
    the reader to re-measure it by hand."""
    cell: dict = {}
    if status is not None:
        cell["state"] = status.get("state")
        if status.get("theoreticalMz") is not None:
            cell["theoreticalMz"] = status["theoreticalMz"]
        if status.get("apexOffset") is not None:
            cell["apexOffset"] = status["apexOffset"]
        if status.get("envelopeCorrelation") is not None:
            cell["correlation"] = status["envelopeCorrelation"]
    if match is not None:
        cell.update({
            "state": "matched",
            "observedMz": match.get("observedMz"),
            "errorPpm": match.get("errorPpm"),
            "isotopeOffset": match.get("observedIsotopeOffset"),
            "correlation": match.get("envelopeCorrelation"),
            "intensity": match.get("intensity"),
        })
    return cell


def _measurement(row: dict) -> dict:
    """The per-model numbers, as the kernel emitted them.

    A one-model run nests them under ``summary``; a comparison puts them on the
    candidate itself. Same numbers either way — read both places rather than
    reporting an empty block for one of the two modes."""
    summary = row.get("summary") or {}
    return {key: (row.get(key, summary.get(key)))
            for key in ("matchedIons", "coverage", "discriminatingMatched",
                        "uniqueCleavages", "medianAbsErrorPpm", "p95AbsErrorPpm")
            if row.get(key, summary.get(key)) is not None}


def _composition_delta(parsed: list[dict], model_ids: list[str]) -> list[dict]:
    """How far apart the models' complete chemistry is, from the kernel's own
    parse of each full model (complete-molecule monoisotopic mass). Not
    derived from fragment differences, which sample only the spans covering
    the differing sites; per-fragment deltas stay readable in each verdict
    row's per-model theoreticalMz.

    This is the context the ion list needs — a ~2 Da difference between two
    near-isobar compositions means isotope borrowing is the thing to watch,
    and an ~80 Da difference means the differing region reaches most of the
    spectrum."""
    base = parsed[0].get("monoisotopicMass")
    out = []
    for model, cid in zip(parsed[1:], model_ids[1:]):
        mass = model.get("monoisotopicMass")
        if base is not None and mass is not None:
            out.append({"model": cid, "vs": model_ids[0],
                        "deltaDa": mass - base})
    return out


def _contested_peaks(cells: dict, tolerance_ppm: float) -> list[dict]:
    """Observed peaks that two different (model, ion) pairs both matched.

    Two models whose compositions differ by a near-integer number of neutrons
    have offset isotope combs: one model's monoisotopic prediction sits where
    the other's M+2 sits, at every charge. The same observed peak is then
    claimed twice — by different ions, at different isotope offsets. The
    kernel's ``isotopeBorrowPairs`` flags one ion matched at two peaks; this
    flags the converse, which no single-model view can see. One peak explained
    by two candidates is a normal state — competing explanations, or
    unresolved contributions within one peak. It means the matches are not
    independent support for either model and the ambiguity deserves a look;
    it does not by itself eliminate any candidate.

    Detected from the kernel's emitted ``observedMz`` only — nothing is
    re-matched here.
    """
    claims: dict[str, list[dict]] = {}
    for key, per_scan in cells.items():
        for scan, model_cells in per_scan.items():
            for model_id, cell in model_cells.items():
                if cell.get("observedMz") is None:
                    continue
                claims.setdefault(scan, []).append({
                    "mz": cell["observedMz"], "model": model_id,
                    "label": _label(key), "isotopeOffset": cell.get("isotopeOffset"),
                    "errorPpm": cell.get("errorPpm"),
                    "intensity": cell.get("intensity"),
                })
    contested = []
    for scan, entries in claims.items():
        entries.sort(key=lambda entry: entry["mz"])
        group: list[dict] = []
        for entry in entries:
            if group and entry["mz"] - group[-1]["mz"] > entry["mz"] * tolerance_ppm * 1e-6:
                if _is_contested(group):
                    contested.append({"scan": scan, "observedMz": group[0]["mz"],
                                      "claims": group})
                group = []
            group.append(entry)
        if _is_contested(group):
            contested.append({"scan": scan, "observedMz": group[0]["mz"],
                              "claims": group})
    contested.sort(key=lambda entry: -max(claim["intensity"] or 0.0
                                          for claim in entry["claims"]))
    return contested[:MAX_CONTESTED]


def _is_contested(group: list[dict]) -> bool:
    """One peak claimed by more than one (model, ion): the same peak cannot be
    two different fragments' evidence at once."""
    return len({(claim["model"], claim["label"]) for claim in group}) > 1


def _verdict_ions(measure: dict, length: int, diff_sites: list,
                  model_ids: list[str], tolerance_ppm: float) -> dict:
    """Every ion whose prediction differs between the models, with the kernel's
    state per model per scan — the discriminating verdict, in one pass."""
    cells: dict[tuple, dict[str, dict[str, dict]]] = {}
    discriminating: set[tuple] = set()
    differing: set[tuple] = set()
    degenerate: set[tuple] = set()
    identical = 0
    for entry in measure["comparisons"]:
        scan = entry["scan"]
        rows = {row["candidateId"]: row for row in candidate_rows(entry)}
        discriminating |= {_key(ion)
                           for ion in entry["measurement"].get("discriminatingIons") or []}
        grids = {
            cid: ({_key(s): s for s in row.get("fragmentStatuses") or []},
                  {_key(m): m for m in row.get("matches") or []})
            for cid, row in rows.items()
        }
        for key in {k for statuses, _m in grids.values() for k in statuses}:
            predicted = [statuses[key]["theoreticalMz"]
                         for statuses, _m in grids.values()
                         if statuses.get(key, {}).get("theoreticalMz") is not None]
            if len(predicted) < 2:
                continue
            if max(predicted) - min(predicted) <= MIN_DELTA_DA:
                identical += 1
                degenerate.add(key)
            else:
                differing.add(key)
            cells.setdefault(key, {})[scan] = {
                cid: _cell(statuses.get(key), matches.get(key))
                for cid, (statuses, matches) in grids.items()
            }

    verdict = []
    unmatched = 0
    for key in sorted(differing, key=lambda key: (key not in discriminating,
                                                  key[0], key[1], key[2])):
        per_scan = cells[key]
        matched_by: dict[str, list[str]] = {}
        for scan, model_cells in per_scan.items():
            for cid, cell in model_cells.items():
                if cell.get("state") == "matched":
                    matched_by.setdefault(cid, []).append(scan)
        if not matched_by:
            unmatched += 1
        reports = _reported_sites(key[0], key[1], length)
        verdict.append({
            "ion": key[0], "number": key[1], "charge": key[2], "label": _label(key),
            "role": ("discriminator" if key in discriminating else
                     "sharedMatch" if matched_by else "noAcceptedMatch"),
            "reportsSites": reports,
            "coversDiffSites": [site for site in diff_sites
                                if (isinstance(site, int) and reports["from"] <= site <= reports["to"])
                                or (site == "nTerm" and key[0] in {"a", "b", "c"})
                                or (site == "cTerm" and key[0] in {"x", "y", "z"})],
            "matchedBy": matched_by,
            **_by_scan({scan: value for scan, value in per_scan.items()},
                       "cells", "cellsByScan"),
        })
    # Most separating first — an ion one model matched and the others did not
    # carries the comparison. Among those: match quality (correlation, then
    # error), then the strongest peak. The isotope slot a match landed in is a
    # fact each row reports, never a rank criterion — envelope matching
    # accepts any slot by design. Unmatched predictions (noAcceptedMatch)
    # last.
    def order(row: dict) -> tuple:
        if row["role"] != "discriminator":
            return (99, 0.0, 0.0, 0.0, row["ion"], row["number"], row["charge"])
        cells = [cell for group in _cell_groups(row) for cell in group.values()
                 if cell.get("state") == "matched"]
        correlations = [cell["correlation"] for cell in cells
                        if cell.get("correlation") is not None]
        errors = [abs(cell["errorPpm"]) for cell in cells
                  if cell.get("errorPpm") is not None]
        peak = max((cell.get("intensity") or 0.0 for cell in cells), default=0.0)
        mean_correlation = sum(correlations) / len(correlations) if correlations else 0.0
        worst_error = max(errors) if errors else 0.0
        return (len(row["matchedBy"]), -mean_correlation, worst_error, -peak,
                row["ion"], row["number"], row["charge"])

    verdict.sort(key=order)
    # Degenerate-but-matched ions stay visible: every model predicts them at
    # the same m/z, so they separate nothing — but they ARE measured, and a
    # span statement ("the whole differing block lies inside this window")
    # legitimately quotes them. Surfaced with the models they cannot split,
    # so the reader never re-derives what an ion excludes by hand.
    degenerate_rows = []
    for key in sorted(degenerate, key=lambda key: (key[0], key[1], key[2])):
        per_scan = cells[key]
        matched_by: dict[str, list[str]] = {}
        for scan, model_cells in per_scan.items():
            for cid, cell in model_cells.items():
                if cell.get("state") == "matched":
                    matched_by.setdefault(cid, []).append(scan)
        if not matched_by:
            continue
        reports = _reported_sites(key[0], key[1], length)
        degenerate_rows.append({
            "ion": key[0], "number": key[1], "charge": key[2], "label": _label(key),
            "role": "degenerateMatch",
            "reportsSites": reports,
            "coversDiffSites": [site for site in diff_sites
                                if (isinstance(site, int) and reports["from"] <= site <= reports["to"])
                                or (site == "nTerm" and key[0] in {"a", "b", "c"})
                                or (site == "cTerm" and key[0] in {"x", "y", "z"})],
            "degenerateAmong": sorted({cid for group in per_scan.values() for cid in group}),
            "matchedBy": matched_by,
            **_by_scan({scan: value for scan, value in per_scan.items()},
                       "cells", "cellsByScan"),
        })
    return {"ions": verdict, "shown": verdict[:MAX_VERDICT_IONS],
            "degenerateMatched": degenerate_rows[:MAX_DEGENERATE_SHOWN],
            "degenerateMatchedTotal": len(degenerate_rows),
            "differingIonKeys": len(differing), "identicalPrediction": identical,
            "noModelMatched": unmatched, "discriminatingIonKeys": len(discriminating),
            "contestedPeaks": _contested_peaks(cells, tolerance_ppm),
            "differences": {
                "predictedDiffer": len(differing),
                "predictedIdentical": identical,
                "measuredDiffer": len(discriminating),
                "note": ("predictedDiffer counts ions whose theoretical m/z "
                         "differs between models; predictedIdentical counts "
                         "ions every model predicts at the same m/z (numerically "
                         "identical) — those "
                         "share m/z predictions; isotope compositions may still differ. "
                         "measuredDiffer counts ions where the kernel's "
                         "accepted matches differ between models. Whether a "
                         "predicted difference is resolvable under this "
                         "measurement depends on tolerance, charge ceiling "
                         "and interference — reported as conditions, not "
                         "pre-judged here."),
                "ambiguity": ("identity, sequence, modification identity and "
                              "site localization are different questions: "
                              "matchedBy on an ion covering a differing site "
                              "speaks to placement; compositionDeltaDa speaks "
                              "to composition; neither alone establishes "
                              "identity."),
            }}


def _cell_groups(row: dict) -> list[dict]:
    """The per-model cell maps of a verdict row, whatever shape it has."""
    return [value for key in ("cells", "cellsByScan")
            if isinstance(row.get(key), dict)
            for value in (row[key].values() if key == "cellsByScan" else [row[key]])]


def _site_windows(measure: dict, diff_sites: list) -> dict:
    """The kernel's cleavage-window counts at the sites the models disagree on."""
    rows = []
    for entry in measure["comparisons"]:
        for row in candidate_rows(entry):
            for window in row.get("siteWindow") or []:
                if window.get("position") in diff_sites:
                    rows.append({"scan": entry["scan"],
                                 "model": row.get("candidateId"), **window})
    return {
        "diffSites": diff_sites,
        "rows": rows[:MAX_SITE_ROWS],
        "truncated": max(0, len(rows) - MAX_SITE_ROWS),
        "note": ("leftCleavage/rightCleavage count the model's own cleavages around "
                 "the site; matchedLeft/matchedRight count how many were observed"),
    }


def _ms2_next(result: dict, borrow: list) -> list[dict]:
    """Capabilities the evidence leaves open — reported, not ordered.

    The card does not decide the next experiment. It states what this
    observation could and could not do, and which capabilities remain
    usable with which inputs; choosing among them belongs to the
    investigation, not to the measurement."""
    capabilities = []
    separated = any(row["role"] == "discriminator" and row["matchedBy"]
                    for row in result["ions"])
    if not separated:
        capabilities.append({
            "status": ("no ion whose prediction differs obtained an accepted "
                       "match separation on the measured scans — this "
                       "observation does not separate these models"),
            "stillUsable": [
                {"capability": "other MS2 scans of the same species",
                 "inputs": "spectral.py measure --mzml <file> --scan <scan[,scan...]> --models <models.json>"},
                {"capability": "precursor-envelope comparison, when the models differ in composition",
                 "inputs": "spectral.py measure --ms1 --feature <label> <mass> <charge> <rtMin> <rtMax> --models <models.json>"},
            ],
        })
    if borrow or result.get("contestedPeaks"):
        capabilities.append({
            "status": ("some separating matches share observed peaks "
                       "(contested or isotope-borrowed) — the ambiguity is "
                       "visible, not resolved"),
            "stillUsable": [
                {"capability": "window inspection of the shared/borrowed ions",
                 "inputs": "spectral.py show <card> --ions <ion:charge,...>"},
            ],
        })
    return capabilities


def build_ms2_card(measure: dict, *, mzml: str, tolerance_ppm: float,
                   title: str | None, evidence: dict, figure: str,
                   warnings: list[str]) -> dict:
    """The MS2 evidence card: what was compared, how far apart the models are,
    which ions separate them, and what the kernel says about the rest."""
    models = measure["candidates"]
    parsed = measure["parsed"]
    model_ids = [model["id"] for model in models]
    shared_backbone = all(p["sequence"] == parsed[0]["sequence"] for p in parsed)
    ion_comparison = shared_backbone and all(
        e["measurement"].get("ionComparisonAvailable", True) for e in measure["comparisons"])
    length = len(parsed[0]["sequence"])
    diff_sites, nterm_differs = _diff_sites(parsed)
    result = _verdict_ions(measure if ion_comparison else {"comparisons": []},
                           length, diff_sites, model_ids, tolerance_ppm)
    borrow = [pair for entry in measure["comparisons"]
              for pair in entry["measurement"].get("isotopeBorrowPairs") or []]
    contested = result["contestedPeaks"]

    notes = list(warnings)
    if len(models) == 1:
        notes.append(
            "one model: this card audits compatibility with the spectrum "
            "(coverage, matched ions, per-model error) and can separate nothing "
            "— pass the competing models together to get a comparison")
    if contested:
        notes.append(
            f"{len(contested)} observed peaks are claimed by more than one "
            "model/ion pair: the models' isotope combs are offset, so the same "
            "peak can be one model's monoisotopic fragment and another's M+2. "
            "Those matches are not independent support for either model — "
            "competing explanations of one peak are a normal state that marks "
            "the evidence as shared or ambiguous and asks for a look, not an "
            "elimination of either candidate. See verdict.contestedPeaks and "
            "the windows via `show`")
    if nterm_differs:
        notes.append(
            "the models differ at the N-terminus: N-terminal (a/b/c) ions "
            "contain that terminal group; C-terminal (x/y/z) fragments do not. "
            "Other modifications may cancel a shift — check coversDiffSites before "
            "reading any ion as separating; a whole-spectrum mass shift is "
            "not a site localization")
    if result["differingIonKeys"] > len(result["shown"]):
        notes.append(f"verdict list truncated to {len(result['shown'])} of "
                     f"{result['differingIonKeys']} differing ions (the most "
                     f"separating first) — full tables: {evidence}")

    if not ion_comparison:
        notes.append("positional ion comparison unavailable: each model was measured on its own backbone; no cross-model site mapping was assumed")

    card = {
        "card": "measure/ms2",
        "title": title or f"{Path(mzml).name} · scan {', '.join(measure['scans'])}",
        "observation": {
            "mzml": mzml,
            "scans": measure["scans"],
            "ionTypes": ion_types_label(
                measure["params"].get("ionTypes")
                or (measure["comparisons"][0].get("ionTypes") if len(measure["comparisons"]) == 1 else None)),
            "activation": measure["params"].get("resolvedActivation"),
            "conditionsByScan": [{key: entry.get(key) for key in
                ("scan", "activation", "ionTypes", "maxCharge", "maxChargeSource",
                 "pooledScans", "binPpm", "minimumAbundancePct")}
                for entry in measure["comparisons"]],
            "maxCharge": measure["params"].get("maxCharge"),
            "maxChargeSource": measure["params"].get("maxChargeSource"),
            "tolerancePpm": tolerance_ppm,
            "minCorrelation": measure["params"].get("minCorrelation"),
            "minIntensity": measure["params"].get("minIntensity"),
            "criterion": ("exact isotope-envelope match on raw centroided MS2: the "
                          "observed cluster is anchored on its most abundant isotope "
                          "and walked with a monotone cursor, missing slots count as "
                          "zero, and the full-slot Pearson must clear the kernel's "
                          "correlation floor. A hit on M+2 is a hit; the isotope "
                          "offset is reported, not normalized away"),
        },
        "sequence": parsed[0]["sequence"] if shared_backbone else None,
        "length": length if shared_backbone else None,
        "sharedBackbone": shared_backbone,
        "ionComparisonAvailable": ion_comparison,
        "models": [{
            "id": model["id"],
            "label": model.get("label", model["id"]),
            "proforma": model["proforma"],
            "bareSequence": parsed[index]["sequence"],
            "modifications": parsed[index]["modifications"],
            **_by_scan({entry["scan"]: _measurement(
                next(row for row in candidate_rows(entry)
                     if row["candidateId"] == model["id"]))
                for entry in measure["comparisons"]},
                "measurement", "measurementByScan"),
        } for index, model in enumerate(models)],
        "verdict": {
            "diffSites": diff_sites,
            "compositionDeltaDa": _composition_delta(measure["parsed"], model_ids),
            "differences": result["differences"],
            "differingIonKeys": result["differingIonKeys"],
            "discriminatingIonKeys": result["discriminatingIonKeys"],
            "ions": result["shown"],
            "degenerateMatched": result["degenerateMatched"],
            "degenerateMatchedTotal": result["degenerateMatchedTotal"],
            "roleGlossary": {
                "discriminator": ("models predict different m/z for this ion; a "
                                  "match here separates them"),
                "sharedMatch": ("matched by some model but not marked "
                                "discriminating by the kernel"),
                "noAcceptedMatch": ("prediction differs but no model obtained an "
                                    "accepted envelope match — absence of "
                                    "acceptance, not evidence of absence"),
                "degenerateMatch": ("every model predicts the same m/z and the "
                                    "ion is matched — measured, but it separates "
                                    "nothing among degenerateAmong; it grounds "
                                    "span/window statements only"),
            },
            "contestedPeaks": contested,
            "contestedPeakNote": ("observed peaks matched by more than one "
                                  "model/ion pair — read them before counting "
                                  "either model's match as independent support"),
            "byScan": [{
                "scan": entry["scan"],
                "discriminatingIonCount": entry["measurement"].get("discriminatingIonCount"),
                "perModel": {row["candidateId"]: _measurement(row)
                             for row in candidate_rows(entry)},
            } for entry in measure["comparisons"]],
        },
        "siteWindows": _site_windows(measure, diff_sites),
        "uninformative": {
            "identicalPrediction": result["identicalPrediction"],
            "noModelMatched": result["noModelMatched"],
            "note": ("identicalPrediction: ions every model predicts at the same "
                     "m/z (numerically identical — not a separability claim) — "
                     "they can "
                     "share m/z, not necessarily isotope composition. noModelMatched: ions whose "
                     "prediction differs but which no model obtained an "
                     "accepted envelope match on any measured scan — an "
                     "unaccepted match is absence of acceptance, not evidence "
                     "of absence, unless observability is argued. The matched "
                     "subset of identicalPrediction is not dropped silently — "
                     "verdict.degenerateMatched lists those ions with the "
                     "models they cannot separate and their full measured "
                     "values. This block "
                     "counts the verdict's remainder; it is not a second "
                     "verdict"),
        },
        "isotopeBorrowPairs": borrow,
        "warnings": notes,
        "audit": {"evidence": evidence, "figure": figure},
    }
    card["next"] = _ms2_next(result, borrow)
    return card


# ------------------------------------------------------------------- verb: measure

def cmd_measure(args, ms_bin: str) -> int:
    if args.ms1:
        return _test_ms1(args, ms_bin)
    if not args.mzml:
        raise ValueError("measure requires --mzml")
    models = _models(args.models)
    scans = ([str(scan) for scan in load_json_arg("--scans-input", args.scans_input)]
             if args.scans_input else
             [scan.strip() for scan in (args.scan or "").split(",") if scan.strip()])
    out = Path(args.out or f"{Path(args.mzml).stem}-card.json")
    _prepare_output(out, args.mzml, args.models)
    figure_path = out.with_name(out.stem + "-figure.html").with_suffix(".html")
    evidence_path = out.with_name(out.stem + "-evidence.json")

    measure = measure_candidates(
        ms_bin, args.mzml, scans, models, str(evidence_path),
        args.ion_types, args.max_charge, args.tolerance_ppm, args.activation,
        selection={"--rt-min": args.rt_min, "--rt-max": args.rt_max,
                   "--scan-min": args.scan_min, "--scan-max": args.scan_max,
                   "--sum-scans": args.sum_scans, "--bin-ppm": args.bin_ppm,
                   "--min-correlation": args.min_correlation, "--min-intensity": args.min_intensity},
    )
    # the figure draws ONE scan; the card covers every measured scan and says
    # which scan the figure shows, rather than leaving them to be assumed equal
    scans = measure["scans"]
    figure_scan = measure["comparisons"][0]["scan"]
    shared_backbone = all(p["sequence"] == measure["parsed"][0]["sequence"] for p in measure["parsed"])
    figure = None
    if not args.no_figures and shared_backbone:
        emitted = ms2_visualize.render_discrimination(
            _namespace(
                mzml=args.mzml, scan=figure_scan, out=str(figure_path),
                candidates=json.dumps(models), ion_types=args.ion_types,
                activation=args.activation, max_charge=args.max_charge,
                tolerance_ppm=args.tolerance_ppm, focus_site=args.focus_site,
                min_rel_intensity=args.min_rel_intensity, title=args.title,
            ),
            ms_bin, measure=measure,
        )
        figure = emitted["figure"]
    card = build_ms2_card(
        measure, mzml=args.mzml, tolerance_ppm=args.tolerance_ppm,
        title=args.title, evidence=_evidence_ref(measure["evidence"]),
        figure=figure, warnings=measure["warnings"],
    )
    card["observation"]["figureScan"] = figure_scan
    if len(scans) > 1:
        card["observation"]["figureScanNote"] = (
            f"the figure draws scan {figure_scan}; the card covers "
            f"{', '.join(scans)}")
    card["audit"]["written"] = [str(out.resolve())]
    if figure is not None:
        card["audit"]["written"].append(figure)
        card["audit"]["figureManifest"] = str(
            Path(figure).with_suffix(".manifest.json"))
    else:
        card["audit"]["figuresSkipped"] = ("--no-figures: measurement only" if args.no_figures
                                             else "overview needs a shared sequence; independent model measurements are retained")
    if not args.no_figures:
        # Default window coverage: the same selection `show` uses — top
        # separating ions plus unmatched windows — so one `measure` call yields
        # card, overview figure and windows, and nothing is left to render
        # by hand. Named ions still go through `show`.
        ions = _show_ions(card, None)
        if ions:
            windows = _render_windows(
                card, ions, out.with_name(out.stem + "-windows.html"),
                args.tolerance_ppm, ms_bin, args.title,
                card_source=str(out.resolve()))
            card["audit"]["windowsFigure"] = windows["figure"]
            card["audit"]["windowsManifest"] = str(
                Path(windows["figure"]).with_suffix(".manifest.json"))
            card["audit"]["written"] += windows["written"]
    _emit_card(card, out)
    return 0


def _namespace(**kwargs) -> types.SimpleNamespace:
    return types.SimpleNamespace(**kwargs)


def _models(value: str) -> list[dict]:
    models = load_json_arg("--models", value)
    if not isinstance(models, list) or not models:
        raise ValueError("--models must be a non-empty JSON array of complete models")
    ids = set()
    for index, model in enumerate(models):
        if not isinstance(model, dict) or not isinstance(model.get("proforma"), str) or not model["proforma"].strip():
            raise ValueError(f"--models[{index}] needs a complete proforma")
        unknown = model.keys() - {"id", "label", "proforma"}
        if unknown:
            raise ValueError(f"--models[{index}] unknown fields: {', '.join(sorted(unknown))}")
        model.setdefault("id", f"model{index + 1}")
        if not isinstance(model["id"], str) or not model["id"].strip() or model["id"] in ids:
            raise ValueError("model ids must be unique nonempty strings")
        if "label" in model and not isinstance(model["label"], str):
            raise ValueError("model labels must be strings")
        ids.add(model["id"])
    return models


def _prepare_output(out: Path, *inputs: str) -> None:
    for value in inputs:
        try:
            if Path(value).is_file() and out.resolve() == Path(value).resolve():
                raise ValueError("output must not overwrite an input")
        except OSError:
            pass  # inline JSON is not a path on every platform
    out.parent.mkdir(parents=True, exist_ok=True)


def _evidence_ref(path: str) -> dict:
    """The artifact's own provenance: where it is and how big it is."""
    return {"path": str(Path(path).resolve()), "size": Path(path).stat().st_size}


def _emit_card(card: dict, out: Path) -> None:
    """The card is the stdout — one document, the same bytes on disk."""
    text = json.dumps(card, ensure_ascii=False, indent=1)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    print(text)


# ------------------------------------------------------------------- verb: MS1

def _test_ms1(args, ms_bin: str) -> int:
    if not args.mzml:
        raise ValueError("measure --ms1 requires --mzml")
    if not args.feature:
        raise ValueError(
            "measure --ms1 requires --feature LABEL MASS CHARGE RT_MIN RT_MAX "
            "(one per observed precursor feature)")
    models = _models(args.models)
    out = Path(args.out or f"{Path(args.mzml).stem}-ms1-card.json")
    figure_path = out.with_name(out.stem + "-figure.html").with_suffix(".html")
    _prepare_output(out, args.mzml, args.models)
    labels = set()
    for feature in args.feature:
        label, mass, charge, rt_min, rt_max = feature
        if not label or label in labels:
            raise ValueError("feature labels must be unique and nonempty")
        labels.add(label)
        values = [float(mass), float(rt_min), float(rt_max)]
        if not all(math.isfinite(v) for v in values) or values[0] <= 0 or int(charge) <= 0 or values[1] > values[2]:
            raise ValueError("feature needs positive mass/charge and ordered finite RT bounds")
    charges = list(dict.fromkeys(int(f[2]) for f in args.feature))
    requests = [{"id": f"m{i}-z{charge}", "proforma": model["proforma"], "charge": charge}
                for i, model in enumerate(models) for charge in charges]
    batch = run_json([ms_bin, "theory", "--models", json.dumps(requests)], label="batch complete-model theory")
    theory_by_id = {row["id"]: row["measurement"] for row in batch["result"]["models"]}
    parsed = [ms2_visualize.parse_proforma(ms_bin, model["proforma"]) for model in models]
    theories = []
    written = []
    for i, model in enumerate(models):
        paths = {}
        for charge in charges:
            path = out.with_name(f"{out.stem}-theory-m{i}-z{charge}.json")
            path.write_text(json.dumps(theory_by_id[f"m{i}-z{charge}"], ensure_ascii=False), encoding="utf-8")
            paths[charge] = str(path)
            written.append(str(path))
        theories.extend(paths[int(feature[2])] for feature in args.feature)
    emitted = envelope_ladder.ladder(_namespace(
        feature_ladder=args.mzml, feature=args.feature, theory=theories,
        theory_label=[model["id"] for model in models],
        cli=ms_bin, tolerance_ppm=args.tolerance_ppm, isotopes=args.isotopes,
        label=args.title or "", plot=str(figure_path), no_figures=args.no_figures,
    ))
    card = _ms1_card(args, models, parsed, emitted, figure=emitted["plot"])
    card["audit"]["written"] = [str(out.resolve()), *written,
                                *emitted.get("written", [])]
    _emit_card(card, out)
    return 0


def _ms1_card(args, models: list[dict], parsed: list[dict], emitted: dict,
              *, figure: str) -> dict:
    """The precursor card.

    MS1 separates compositions, not sites: the numbers that matter are each
    model's composition mass (its theoretical envelope's anchor) against the
    observed feature mass, and the shape agreement of the two isotope
    distributions. No acceptance threshold is invented here — the offsets and
    correlations are reported, and the calibration of the observed mass is the
    run's, not this tool's.
    """
    by_model: dict[str, dict] = {model["id"]: {} for model in models}
    features = []
    for member in emitted["members"]:
        ranking = []
        for theory in member["theories"]:
            key = theory["label"]
            row = {"model": key, "anchorMass": theory["anchor_mass"],
                   "anchorOffsetDa": theory["anchor_offset_da"],
                   "pearsonR": theory["pearson_r"],
                   "matchedPeaks": theory["matched_peaks"],
                   "fitMethod": "kernel fit-envelope; matched slots; minimumWindowMz=0",
                   "apexDeltaDa": theory["apex_delta_da"]}
            by_model.setdefault(key, {})[member["label"]] = row
            ranking.append(row)
        ranking.sort(key=lambda row: (row["anchorOffsetDa"] is None,
                                      abs(row["anchorOffsetDa"] or 0.0)))
        offsets = [abs(row["anchorOffsetDa"]) for row in ranking
                   if row["anchorOffsetDa"] is not None]
        correlations = [row["pearsonR"] for row in ranking
                        if row["pearsonR"] is not None]
        features.append({
            "label": member["label"], "observedMass": member["mass"],
            "charge": member["charge"], "rtMin": member["rt_min"],
            "rtMax": member["rt_max"], "scanCount": member["scan_count"],
            "ranking": ranking,
            # descriptive spreads, not acceptance thresholds: the card says how
            # far apart the models' predictions are, not which one is true
            "anchorOffsetSpreadDa": (max(offsets) - min(offsets)
                                     if len(offsets) > 1 else None),
            "pearsonSpread": (max(correlations) - min(correlations)
                              if len(correlations) > 1 else None),
        })
    return {
        "card": "measure/ms1",
        "title": args.title or f"{Path(args.mzml).name} · precursor envelopes",
        "observation": {
            "mzml": args.mzml,
            "features": [feature["label"] for feature in features],
            "tolerancePpm": args.tolerance_ppm,
            "isotopes": args.isotopes,
            "criterion": ("the RT-window-aggregated observed precursor envelope "
                          "against each model's theoretical isotope distribution, "
                          "on one shared axis of mass offset from the feature's "
                          "neutral mass. Pearson r is shape agreement, not an "
                          "identity score; the anchor offset is where the model's "
                          "composition sits relative to the measured feature mass"),
        },
        "sequence": parsed[0]["sequence"] if all(p["sequence"] == parsed[0]["sequence"] for p in parsed) else None,
        "length": len(parsed[0]["sequence"]) if all(p["sequence"] == parsed[0]["sequence"] for p in parsed) else None,
        "models": [{
            "id": model["id"], "label": model.get("label", model["id"]),
            "proforma": model["proforma"],
            "bareSequence": parsed[index]["sequence"],
            "modifications": parsed[index]["modifications"],
            "perFeature": by_model.get(model["id"], {}),
        } for index, model in enumerate(models)],
        "verdict": {"features": features},
        "siteWindows": None,
        "uninformative": {
            "note": ("precursor envelopes constrain composition and mass, never "
                     "site: two models that differ only in WHERE a modification "
                     "sits have the same precursor envelope, and this card cannot "
                     "separate them"),
        },
        "warnings": [],
        "next": [{
            "status": ("MS1 ranks compositions; it cannot place a modification — "
                       "models differing only in placement are indistinguishable "
                       "here by construction"),
            "stillUsable": [
                {"capability": "MS2 fragment comparison for site placement",
                 "inputs": ("spectral.py measure --mzml <file> "
                            "--scan <scan[,scan...]> --models <models.json>")},
            ],
        }],
        "audit": {"figure": figure, "fits": emitted["fits"], "observations": emitted["observations"]},
    }


# ----------------------------------------------------------------- verb: annotate

def cmd_annotate(args, ms_bin: str) -> int:
    if not args.msalign:
        raise ValueError("annotate requires --msalign")
    if not args.proforma:
        raise ValueError("annotate requires --proforma (one complete model)")
    out = Path(args.out or f"{Path(args.msalign).stem}-scan{args.scan}-annotation.html")
    evidence_path = out.with_suffix(".fragment-match.json")
    emitted = ms2_visualize.render_msalign_annotation(
        _namespace(
            msalign=args.msalign, scan=str(args.scan), proforma=args.proforma,
            out=str(out), ion_types=args.ion_types, activation=args.activation,
            max_charge=args.max_charge, tolerance_ppm=args.tolerance_ppm,
            annotate=args.annotate, evidence_out=str(evidence_path),
        ),
        ms_bin,
    )
    card = {
        "card": "annotate/msalign",
        "title": f"{Path(args.msalign).name} · scan {args.scan}",
        "observation": {
            "msalign": args.msalign, "scan": str(args.scan),
            "tolerancePpm": args.tolerance_ppm,
            "criterion": ("charge-aware neutral-mass matching in the deconvoluted "
                          "domain: TopFD has already deconvoluted the spectra, so "
                          "a miss here carries no existence verdict — collapse, "
                          "deconvolution failure, filtering, interference and "
                          "insufficient signal all look the same. Raw "
                          "isotope-envelope matching does not apply here"),
        },
        "model": {"proforma": args.proforma},
        "verdict": {"matchedFragments": emitted["matchedFragments"],
                    "coverage": emitted["coverage"]},
        "next": [{
            "status": ("the annotation shows where the leads are, not whether "
                       "the model beats an alternative"),
            "stillUsable": [
                {"capability": "raw-mzML envelope matching, the criterion",
                 "inputs": ("spectral.py measure --mzml <source mzML> "
                            "--scan <linked scan> --models <models.json>")},
            ],
        }],
        "audit": {
            "evidence": _evidence_ref(str(evidence_path)),
            "figure": emitted["figure"],
            "written": [emitted["figure"], str(evidence_path)],
        },
    }
    _emit_card(card, out.with_name(out.stem + "-card.json"))
    return 0


# --------------------------------------------------------------------- verb: show

def _parse_ion_spec(spec: str) -> tuple:
    """``c2:1`` -> ``("c", 2, 1)``.

    The charge is required: a bare ion expands to every charge state, and the
    window would stop being a window."""
    head, sep, charge = spec.partition(":")
    if (not sep or not head[:1].isalpha() or not head[1:].isdigit()
            or not charge.isdigit() or head[0].lower() not in "abcxyz"
            or int(head[1:]) <= 0 or int(charge) <= 0):
        raise ValueError(f"--ions expects ION:CHARGE such as c2:1 or y3:2, got {spec!r}")
    return (head[0].lower(), int(head[1:]), int(charge))


def _show_ions(card: dict, requested: list[str] | None) -> list[tuple[str, str]]:
    """``(ion label, the scan to extract it from)`` for the windows to draw.

    Default: the card's most separating verdict ions, then enough controls to
    show that the window is not simply empty everywhere."""
    if requested:
        # one --ions flag per call or a comma list, as the kernel's --ions takes
        specs = [spec.strip() for value in requested
                 for spec in str(value).split(",") if spec.strip()]
        for spec in specs:
            _parse_ion_spec(spec)
        return [(spec, "") for spec in specs]
    rows = (card.get("verdict") or {}).get("ions") or []
    picked = [row for row in rows
              if row["role"] == "discriminator" and row.get("matchedBy")][:DEFAULT_SHOW_IONS]
    # windows with no accepted match make the absence visible; whether they
    # act as controls depends on an observability argument, not on this pick
    picked += [row for row in rows if row["role"] == "noAcceptedMatch"][:2]
    if not picked:
        picked = rows[:DEFAULT_SHOW_IONS]
    out = []
    for row in picked:
        scans = [scan for scans in row.get("matchedBy", {}).values() for scan in scans]
        scan = scans[0] if scans else (card["observation"]["scans"] or [""])[0]
        out.append((row["label"], str(scan)))
    return out


def _show_sub(row: dict | None, model_id: str, ion_label: str, scan: str) -> str:
    """The subtitle: what the kernel measured here, not what to conclude."""
    if row is None:
        return ion_label
    group = row.get("cells") or (row.get("cellsByScan") or {}).get(scan, {})
    cell = group.get(model_id)
    if not cell:
        return f"{ion_label} · not a differing prediction for this model"
    bits = [f"kernel: {cell['state']}"]
    if cell.get("errorPpm") is not None:
        bits.append(f"{cell['errorPpm']:+.1f} ppm")
    if cell.get("isotopeOffset") is not None:
        bits.append(f"M+{cell['isotopeOffset']}")
    if cell.get("correlation") is not None:
        bits.append(f"r={cell['correlation']:.2f}")
    elif cell["state"] == "missingApex" and cell.get("apexOffset") is not None:
        bits.append(f"apex expected at M+{cell['apexOffset']}")
    elif cell["state"] == "correlationRejected":
        # CorrelationRejected with no number is the degenerate-variance gate,
        # not a low score: the kernel never formed a Pearson here
        bits.append("no correlation could be formed")
    return " · ".join(bits)


def _show_ms1(args) -> int:
    """Re-expose the MS1 figure routes: an already-fitted precursor envelope
    (refine-shift) or an observed feature against complete-model theories.

    Rendering a measurement that already exists on disk is the same verb as
    rendering MS2 windows — look at what was measured — so it lives here
    rather than behind a second entry point the agent has to discover."""
    from ms1_visualize import render_refine, run_theory_overlay

    if bool(args.refine) == bool(args.observed):
        raise ValueError(
            "show --ms1 takes exactly one source: --refine <refine-shift JSON> "
            "(repeat to compare candidates) or --observed <feature JSON> with "
            "one --theory per model")
    plot = args.out or (
        f"{Path(args.refine[0]).stem}-envelope.html" if args.refine
        else f"{Path(args.observed).stem}-envelope.html")
    ns = _namespace(
        refine=args.refine or [], observed=args.observed,
        theory=args.theory or [], fit=args.fit or [],
        theory_label=args.theory_label or [], label=args.title or "",
        mode="evidence", plot=plot, plot_data=True,
    )
    payload = run_theory_overlay(ns) if args.observed else render_refine(ns)
    payload["written"] = [payload.get("figure") or payload.get("plot")]
    print(json.dumps(payload, ensure_ascii=False))
    return 0


def _render_windows(card: dict, ions: list[tuple[str, str]], out: Path, ppm: float,
                    ms_bin: str, title: str | None, card_source: str) -> dict:
    """Render the isotope-window mirror figure for `(ion, scan)` pairs.

    The one renderer behind both verbs: `show` for ions the caller names,
    `measure` for its default selection, so there is exactly one way these
    windows are drawn and no reason to hand-roll a variant."""
    observation = card["observation"]
    mzml = observation["mzml"]
    # the proformas travel in the card, so rendering needs nothing but the card
    proformas = {model["id"]: model["proforma"] for model in card["models"]}
    rows = {(row["ion"], int(row["number"]), int(row["charge"])): row
            for row in (card.get("verdict") or {}).get("ions") or []}
    model_ids = [model["id"] for model in card["models"]]
    entries = [(ion_label, scan or observation["scans"][0], model_id)
               for ion_label, scan in ions for model_id in model_ids]
    # Theory and observation are independent atoms: generate theory in one
    # batch, crop the exact saved observation for display, and keep the kernel's
    # existing match verdict. This also handles pooled and saved spectra.
    artifact = ms2_visualize.load_evidence(card["audit"]["evidence"]["path"])
    measured = {str(entry["scan"]): (measurement, entry)
                for measurement, entry in ms2_visualize.comparison_entries(artifact)}
    requested = [_parse_ion_spec(label) for label, _ in ions]
    theory = run_json([ms_bin, "fragment-theoretical", "--models",
        json.dumps([{"id": m["id"], "proforma": m["proforma"]} for m in card["models"]]),
        "--ion-types", ",".join(sorted({key[0].upper() for key in requested})),
        "--max-charge", str(max(key[2] for key in requested))], label="batch window theory")
    theory_path = out.with_name(out.stem + "-theory.json")
    theory_path.parent.mkdir(parents=True, exist_ok=True)
    theory_path.write_text(json.dumps(theory, ensure_ascii=False), encoding="utf-8")
    grids = {model["id"]: {_key(f): f for f in model["measurement"]["fragments"]}
             for model in theory["result"]["models"]}
    states = {}
    for scan, (measurement, entry) in measured.items():
        for candidate in measurement["candidates"]:
            matches = {_key(m): m for m in candidate.get("matches") or []}
            for status in candidate.get("fragmentStatuses") or []:
                states[(scan, candidate["candidateId"], _key(status))] = _cell(status, matches.get(_key(status)))
    panels, extractions = [], []
    for i, (label, scan, model_id) in enumerate(entries):
        key = _parse_ion_spec(label)
        if scan not in measured:
            raise ValueError(f"scan {scan} is not in this measurement")
        _, entry = measured[scan]
        peaks = entry.get("pooledPeaks") or entry.get("observedPeaks")
        if peaks is None:
            raise ValueError("window inspection needs observed peaks in the full evidence; remeasure the older artifact")
        sticks = (grids[model_id].get(key) or {}).get("isotopeEnvelope") or []
        observed = []
        if sticks:
            lo, hi = min(p["mz"] for p in sticks) - 0.5, max(p["mz"] for p in sticks) + 0.5
            for peak in peaks:
                if lo <= peak["mz"] <= hi:
                    nearest = min(sticks, key=lambda t: abs(t["mz"] - peak["mz"]))
                    observed.append({**peak, "nearestNeutronOffset": nearest["neutronOffset"],
                                     "errPpm": (peak["mz"] - nearest["mz"]) / nearest["mz"] * 1e6})
        region = {"ion": key[0], "number": key[1], "charge": key[2],
                  "theoretical": sticks, "observed": observed}
        path = out.with_name(f"{out.stem}.regions{i}.json")
        projection = {"schemaVersion": 1, "params": {"scan": scan},
                      "source": {"evidence": card["audit"]["evidence"]["path"], "theory": str(theory_path)},
                      "result": {"regions": [region]}}
        path.write_text(json.dumps(projection, ensure_ascii=False), encoding="utf-8")
        panel = envelope_mirror._panels(projection, str(path), ppm)[0]
        panel["title"] = f"{model_id} · {label} · scan {scan}"
        if not sticks:
            panel["subtitle"] = f"{label} · no defined fragment for this model"
        panels.append(panel)
        extractions.append(str(path))
    for panel, (ion_label, scan, model_id) in zip(panels, entries):
        # the window's own count AND the kernel's verdict on it, on one line:
        # "10/14 slots" and "kernel: correlationRejected" are different
        # statistics, and reading either alone misleads
        panel["subtitle"] = (
            f"{_show_sub({'cells': {model_id: states.get((scan, model_id, _parse_ion_spec(ion_label)), {})}}, model_id, ion_label, scan)}"
            f" · {panel['subtitle']}")
    for index, panel in enumerate(panels):
        panel["color"] = SHOW_COLORS[index % len(SHOW_COLORS)]

    payload = {
        "data": {
            "panels": [{key: value for key, value in panel.items()
                        if not key.startswith("_")} for panel in panels],
            "xlabel": "m/z",
            "ylabel": "rel. intensity (%)",
        },
        "meta": {"title": title or f"MS2 isotope windows · {Path(mzml).name}",
                 "source": mzml,
                 "n": f"{len(panels)} panels · {len(ions)} ions"},
    }
    html, _ = emit_html("envelope-mirror", payload, out,
                        source_default="tdeasenexus-toolbox envelope-regions artifacts")
    facts = [{**panel["_facts"], "label": panel["title"],
              "subtitle": panel["subtitle"]} for panel in panels]
    write_manifest(out, figure="envelope-mirror", script=Path(__file__).name,
                   params={"card": card_source, "ppm": ppm,
                           "ions": [ion for ion, _ in ions]},
                   inputs=[card_source, card["audit"]["evidence"]["path"], str(theory_path), *extractions], caption={"panels": facts},
                   plot_data=facts)
    return {
        "figure": str(Path(html).resolve()),
        "renderer": "chart-forge d3",
        "card": card.get("title"),
        "ions": [{"ion": ion, "scan": scan,
                  "role": (rows.get(_parse_ion_spec(ion)) or {}).get("role")}
                 for ion, scan in ions],
        "panels": facts,
        "extractions": extractions,
        "written": [str(html), str(theory_path), *extractions],
        "note": ("each panel is one model's theoretical isotope comb against the "
                 "raw peaks in the same window; the verdict is in the card, this "
                 "figure is where you check it"),
    }


def cmd_show(args, ms_bin: str) -> int:
    if args.ms1:
        return _show_ms1(args)
    if not args.card:
        raise ValueError("show requires a card.json written by `measure --mzml --scan` "
                         "(or --ms1 to render an existing MS1 measurement)")
    card = load_doc(args.card)
    if card.get("card") != "measure/ms2":
        raise ValueError(
            f"{args.card} is a {card.get('card')!r} card; show renders MS2 "
            "isotope windows from a `measure --mzml --scan` card (MS1's figure is "
            "already the observed-vs-theory overlay, and the msalign annotation "
            "is drawn by `annotate`)")
    ions = _show_ions(card, args.ions)
    if not ions:
        raise ValueError("no default separating windows; name the desired ions with --ions ION:CHARGE")
    out = Path(args.out or Path(args.card).with_name(Path(args.card).stem + "-windows.html"))
    result = _render_windows(card, ions, out, args.ppm, ms_bin, args.title,
                             card_source=args.card)
    print(json.dumps(result, ensure_ascii=False))
    return 0


# -------------------------------------------------------------------- dispatch

def main(argv: list[str] | None = None) -> int:
    # stdout/stderr carry the card: they are UTF-8 regardless of the console's
    # code page, so a non-ASCII sequence, title or OS error message survives
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8")

    # an invocation mistake belongs on the same JSON surface as a tool failure:
    # argparse's own exit-2 usage text would be the one output shape the agent
    # has no reader for
    parser = AgentArgumentParser(
        prog="spectral.py",
        description="One spectral-evidence entry: annotate (msalign), measure "
                    "(the criterion, as a card), show (the windows).")
    verbs = parser.add_subparsers(dest="verb", required=True)

    def common(sub):
        sub.add_argument("--ms-bin", help="tdeasenexus-toolbox binary")
        sub.add_argument("--ion-types",
                         help="fragment series; inferred from the activation when omitted")
        sub.add_argument("--activation",
                         help="activation label override; otherwise read from source metadata")
        sub.add_argument("--max-charge", type=int,
                         help="theoretical fragment charge ceiling; omitted, it follows "
                              "the precursor's declared charge (the kernel's physical ceiling)")
        sub.add_argument("--tolerance-ppm", type=float, default=10.0)
        sub.add_argument("--out")
        sub.add_argument("--title")

    annotate = verbs.add_parser("annotate",
                                help="annotate one model on a deconvoluted msalign scan")
    annotate.add_argument("--msalign", required=True)
    annotate.add_argument("--scan", required=True)
    annotate.add_argument("--proforma", required=True)
    annotate.add_argument("--annotate", action="append", default=[], metavar="SPEC",
                          help="attention pointer on the figure ('mz=LABEL' or "
                               "'mz1..mz2=LABEL'); not measured data")
    common(annotate)

    measure = verbs.add_parser("measure",
                            help="the criterion: measure models against one observation")
    measure.add_argument("--mzml", "--input", dest="mzml", required=True,
                      help="centroided mzML or saved spectrum sum/window JSON")
    measure.add_argument("--models", required=True,
                      help='JSON array [{"id","label","proforma"}, ...] — inline or a '
                           "path; every model complete (sequence, terminal forms, all "
                           "modifications); each model owns its backbone")
    measure.add_argument("--scan", help="scan number, or a comma list to measure several scans")
    measure.add_argument("--scans-input", help="inline JSON or a file of scan IDs")
    measure.add_argument("--rt-min", type=float)
    measure.add_argument("--rt-max", type=float)
    measure.add_argument("--scan-min", type=int)
    measure.add_argument("--scan-max", type=int)
    measure.add_argument("--sum-scans", action="store_true", help="explicitly pool compatible scans")
    measure.add_argument("--bin-ppm", type=float, default=10.0, help="pooling bin tolerance, separate from matching")
    measure.add_argument("--ms1", action="store_true",
                      help="ask about precursor envelopes instead of fragments")
    measure.add_argument("--feature", action="append", nargs=5,
                      metavar=("LABEL", "MASS", "CHARGE", "RT_MIN", "RT_MAX"),
                      help="--ms1: observed precursor feature; repeat for a family")
    measure.add_argument("--isotopes", type=int, default=0,
                      help="--ms1: isotope count; 0 lets the kernel derive it")
    measure.add_argument("--focus-site", type=int,
                      help="the site the figure's local summaries centre on; defaults "
                           "to the first site where the models' chemistry differs")
    measure.add_argument("--min-correlation", type=float, default=0.7)
    measure.add_argument("--min-intensity", type=float, help="absolute matching intensity floor, separate from display")
    measure.add_argument("--min-rel-intensity", type=float, default=0.01,
                      help="spectrum display floor, as a fraction of the base peak")
    measure.add_argument("--no-figures", action="store_true",
                      help="measurement only: write the card and evidence, skip "
                           "the overview figure and the default isotope windows")
    common(measure)

    show = verbs.add_parser("show", help="draw the windows of a measurement")
    show.add_argument("card", nargs="?",
                      help="the card.json written by `measure --mzml --scan`")
    show.add_argument("--ms1", action="store_true",
                      help="render an existing MS1 measurement instead: a "
                           "refine-shift result, or an observed feature against "
                           "complete-model theories")
    show.add_argument("--ions", action="append", metavar="ION:CHARGE",
                      help="ion to draw, e.g. c2:1; repeat or comma-separate. "
                           "Default: the card's separating ions plus two "
                           "unmatched windows")
    show.add_argument("--ppm", type=float, default=10.0,
                      help="in-slot tolerance for the panel annotations (default 10)")
    show.add_argument("--refine", action="append", metavar="JSON",
                      help="--ms1: refine-shift JSON (observed + theoretical "
                           "envelopes); repeat to compare candidates")
    show.add_argument("--observed", metavar="JSON",
                      help="--ms1: one envelope feature JSON; overlays every "
                           "--theory against it")
    show.add_argument("--theory", action="append", metavar="JSON",
                      help="--ms1: theoretical envelope JSON, one per model")
    show.add_argument("--fit", action="append", metavar="JSON",
                      help="--ms1: fit-envelope JSON, one per --theory; metrics "
                           "are bound to the pair by the fit's own recorded inputs")
    show.add_argument("--theory-label", action="append", metavar="LABEL",
                      help="--ms1: one label per --theory")
    show.add_argument("--ms-bin", help="tdeasenexus-toolbox binary")
    show.add_argument("--out", metavar="PLOT", help="figure stem; writes <stem>.html")
    show.add_argument("--title")

    try:
        args = parser.parse_args(argv)
        if args.verb in {"measure", "annotate"}:
            if not math.isfinite(args.tolerance_ppm) or args.tolerance_ppm <= 0:
                raise ValueError("--tolerance-ppm must be positive and finite")
            if args.max_charge is not None and args.max_charge <= 0:
                raise ValueError("--max-charge must be positive")
        if args.verb == "measure":
            if not math.isfinite(args.min_correlation) or not -1 <= args.min_correlation <= 1:
                raise ValueError("--min-correlation must be finite and between -1 and 1")
            if args.min_intensity is not None and (not math.isfinite(args.min_intensity) or args.min_intensity < 0):
                raise ValueError("--min-intensity must be finite and nonnegative")
            if args.scans_input:
                scans = load_json_arg("--scans-input", args.scans_input)
                if not isinstance(scans, list) or not scans or any(
                    isinstance(v, bool) or not isinstance(v, (str, int)) or not str(v).strip() for v in scans):
                    raise ValueError("--scans-input needs a nonempty array of integer or native scan IDs")
            for lo, hi in ((args.rt_min, args.rt_max), (args.scan_min, args.scan_max)):
                if any(v is not None and (not math.isfinite(v) or v < 0) for v in (lo, hi)) or (lo is not None and hi is not None and lo > hi):
                    raise ValueError("selection bounds must be finite, nonnegative and ordered")
            if args.ms1 and (args.scan or args.scans_input or args.sum_scans or any(v is not None for v in (args.rt_min, args.rt_max, args.scan_min, args.scan_max))):
                raise ValueError("--ms1 selects RT windows with --feature; fragment scan selectors do not apply")
            if args.scan and args.scans_input:
                raise ValueError("choose --scan or --scans-input")
            if not math.isfinite(args.min_rel_intensity) or not 0 <= args.min_rel_intensity <= 1:
                raise ValueError("--min-rel-intensity must be between 0 and 1")
        if args.verb == "show" and (not math.isfinite(args.ppm) or args.ppm <= 0):
            raise ValueError("--ppm must be positive and finite")
        ms_bin = ms_binary(args.ms_bin)
        if args.verb == "measure":
            return cmd_measure(args, ms_bin)
        if args.verb == "annotate":
            return cmd_annotate(args, ms_bin)
        return cmd_show(args, ms_bin)
    except SystemExit:
        raise
    except Exception as exc:  # noqa: BLE001 — one error contract for the agent
        print(json.dumps({"error": {"code": CODE, "message": str(exc)}},
                         ensure_ascii=False))
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
