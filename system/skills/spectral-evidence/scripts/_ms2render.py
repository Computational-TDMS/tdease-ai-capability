# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Shared machinery for the spectral-evidence MS2 entry.

``spectral.py`` is the one agent-facing entry; every MS2 route runs the same
kernel-backed pipeline through :func:`measure_candidates`: parse the ProForma,
run one batched `fragment --mode compare` invocation for all models and scans, read the evidence artifact back, summarize
per-ion-series continuity, then write data.json -> chart-forge HTML -> manifest
-> one stdout JSON. This module is that pipeline; the renderers are import
targets, never CLIs themselves.

The evidence artifact and the stdout envelope are ONE schema (the same
envelope, stdout a bounded projection of it), so a single reader serves both.
:func:`load_evidence` stays tolerant of the earlier two-shape artifacts
(schemaVersion 2 artifacts, and comparison entries that carried their payload
under ``result`` rather than ``measurement``) so an analysis directory written
before the unify keeps working.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from _shared.chart_forge import emit_html
from _shared.kernel_io import load_doc
from _shared.manifest import write_manifest
from _shared.native import binary as ms_binary, run_json  # noqa: F401 (re-export)


class AgentArgumentParser(argparse.ArgumentParser):
    """Keep invocation mistakes on the same JSON error surface as tool failures."""

    def error(self, message: str) -> None:
        raise ValueError(f"invalid arguments: {message}; use --help for the command contract")


def run(ms_bin: str, args: list[str]) -> dict:
    return run_json([ms_bin, *args], label=f"tdeasenexus-toolbox {args[0]}")


def parse_proforma(ms_bin: str, proforma: str) -> dict:
    return run(ms_bin, ["parse-proforma", "--proforma", proforma])["result"]


def mzml_scan(
    ms_bin: str, mzml: str, scan, *, max_peaks: int | None = None
) -> tuple[list[dict], str | None, int | None]:
    argv = [
        "spectrum", mzml, "--mode", "scans", "--scans", str(scan),
        "--peak-order", "mzAscending",
    ]
    if max_peaks is not None:
        argv += ["--max-peaks", str(max_peaks)]
    result = run(ms_bin, argv)["result"]
    spectra = result.get("spectra") or []
    if len(spectra) != 1:
        raise ValueError(f"scan {scan} was not found exactly once in {mzml}")
    spectrum = spectra[0]["spectrum"]
    precursor = spectrum.get("precursor") or {}
    peaks = [
        {"mz": float(peak["mz"]), "i": float(peak["intensity"])}
        for peak in spectrum.get("peaks") or []
    ]
    return peaks, precursor.get("activation"), precursor.get("charge")


def fragment_match(ms_bin: str, spectrum_file: str, scan, proforma: str,
                   evidence_out: str, ion_types: str | None,
                   max_charge: int | None, tolerance_ppm: float,
                   activation: str | None = None) -> dict:
    """One kernel match; returns the parsed evidence-file ``result``.

    mzML input always matches at the exact isotope-envelope level; msalign
    input is the deconvoluted domain (charge-aware neutral-mass one-to-one)
    — the kernel picks the route from the file, there is no toggle.
    ``max_charge=None`` omits the kernel flag so the fragment charge ceiling
    follows the precursor's declared charge (the physical ceiling); the
    effective ceiling and its provenance are attached as ``_maxCharge`` /
    ``_maxChargeSource``.
    """
    argv = [
        "fragment", spectrum_file, "--mode", "match", "--scan", str(scan),
        "--proforma", proforma, "--evidence-out", evidence_out,
    ]
    if ion_types:
        argv += ["--ion-types", ion_types]
    if max_charge is not None:
        argv += ["--max-charge", str(max_charge)]
    argv += ["--tolerance-ppm", str(tolerance_ppm)]
    if activation:
        argv += ["--activation", activation]
    response = run(ms_bin, argv)
    result = match_measurement(load_evidence(response["result"]["evidenceRef"]["path"]))
    result["_resolvedIonTypes"] = response.get("params", {}).get("ionTypes")
    result["_resolvedActivation"] = response.get("params", {}).get("resolvedActivation")
    result["_maxCharge"] = response.get("params", {}).get("maxCharge")
    result["_maxChargeSource"] = response.get("params", {}).get("maxChargeSource")
    return result


def fragment_compare(ms_bin: str, spectrum_file: str, scan, candidate_proformas: list,
                     evidence_out: str, ion_types: str | None,
                     max_charge: int | None,
                     tolerance_ppm: float,
                     activation: str | None = None) -> tuple[list[dict], dict]:
    """One kernel comparison for every candidate — one shared peak extraction.

    Returns (per-candidate results in ``--candidate-proformas`` order, the
    response params). Full match tables are file-backed: the comparison
    evidence sidecar (schemaVersion 2, no ``result`` wrapper) holds one
    ``comparisons[]`` entry per scan with the complete
    ``result.candidates[]``; the stdout envelope carries only the compact
    measurement view. mzML comparison is always isotope-envelope level.
    ``max_charge=None`` follows the precursor's declared
    charge; the effective ceiling arrives back in ``params["maxCharge"]``.
    """
    argv = [
        "fragment", spectrum_file, "--mode", "compare", "--scan", str(scan),
        "--candidate-proformas", json.dumps(candidate_proformas),
        "--evidence-out", evidence_out,
    ]
    if ion_types:
        argv += ["--ion-types", ion_types]
    if max_charge is not None:
        argv += ["--max-charge", str(max_charge)]
    argv += ["--tolerance-ppm", str(tolerance_ppm)]
    if activation:
        argv += ["--activation", activation]
    response = run(ms_bin, argv)
    evidence = load_evidence(response["result"]["evidenceRef"]["path"])
    candidates = [candidate for measurement, _entry in comparison_entries(evidence)
                  for candidate in measurement["candidates"]]
    return candidates, response.get("params", {})


# ------------------------------------------------------------- evidence readers

def load_evidence(path: str) -> dict:
    """One fragment evidence artifact, or a captured stdout envelope.

    They are the same document (the artifact is the envelope with the
    unbounded tables kept in place); the tolerant ``measurement``/``result``
    seam below keeps artifacts written before that unify readable.
    """
    doc = load_doc(path)
    if isinstance(doc.get("result"), dict):
        return doc
    # An artifact written before the envelope unify: schemaVersion 2, the
    # comparison list at the top level rather than inside a `result`. Reading
    # an older analysis directory must not require re-measuring it.
    if isinstance(doc.get("comparisons"), list):
        return {"result": {"comparisons": doc["comparisons"],
                           "aggregation": doc.get("aggregation")}}
    raise ValueError(
        f"{path}: not a fragment artifact — no 'result' object "
        "(expected a `fragment` --evidence-out artifact or a captured "
        "stdout envelope)")


def comparison_entries(doc: dict) -> list[tuple[dict, dict]]:
    """``(measurement, entry)`` per scan from a compare artifact."""
    entries = doc["result"].get("comparisons")
    if not isinstance(entries, list) or not entries:
        raise ValueError("not a fragment comparison artifact — no comparisons[]")
    out = []
    for index, entry in enumerate(entries):
        payload = entry.get("measurement") or entry.get("result")
        if not isinstance(payload, dict):
            raise ValueError(f"comparisons[{index}]: no measurement object")
        out.append((payload, entry))
    return out


def match_measurement(doc: dict) -> dict:
    """The single-model measurement of a match artifact."""
    result = doc["result"]
    measurement = result.get("measurement") or result
    if not isinstance(measurement.get("matches"), list):
        raise ValueError("not a fragment match artifact — no matches[]")
    return measurement


def ion_types_label(value) -> str:
    """The fragment series as one display string.

    The kernel reports them lowercase inside a comparison entry and uppercase
    in ``params``; the figure and the card record the uppercase form whichever
    side the answer came from, so the same measurement does not read two ways.
    """
    if isinstance(value, str):
        return value
    return ",".join(str(part).upper() for part in value or [])


def candidate_rows(entry: dict) -> list[dict]:
    """The per-candidate summaries of one comparison entry.

    A one-model run reports its single summary inline instead of a list;
    normalizing here keeps the card and the figure indifferent to which
    kernel mode produced them."""
    measurement = entry["measurement"]
    candidates = measurement.get("candidates")
    if isinstance(candidates, list) and candidates:
        return candidates
    return [measurement]


# --------------------------------------------------------- one measurement pass

def measure_candidates(ms_bin: str, spectrum_file: str, scans: list[str],
                       cands: list[dict], evidence_out: str,
                       ion_types: str | None = None,
                       max_charge: int | None = None,
                       tolerance_ppm: float = 10.0,
                       activation: str | None = None, *,
                       selection: dict | None = None) -> dict:
    """One batch kernel invocation for all models and selected observations.

    Each model owns its sequence. Positional comparisons are available only
    when the kernel reports a shared backbone; independent fits remain useful.
    Full observed peaks travel with the evidence, so rendering never reopens
    the source to reconstruct an observation.
    """
    if not cands:
        raise ValueError("no models: pass at least one complete model")
    parsed = [parse_proforma(ms_bin, cand["proforma"]) for cand in cands]
    argv = ["fragment", spectrum_file, "--mode", "compare",
            "--candidate-proformas", json.dumps(cands),
            "--evidence-out", evidence_out,
            "--tolerance-ppm", str(tolerance_ppm)]
    if scans:
        argv += ["--scans-input", json.dumps(scans)]
    for flag, value in (selection or {}).items():
        if value is True:
            argv.append(flag)
        elif value is not None and value is not False:
            argv += [flag, str(value)]
    if ion_types:
        argv += ["--ion-types", ion_types]
    if max_charge is not None:
        argv += ["--max-charge", str(max_charge)]
    if activation:
        argv += ["--activation", activation]
    response = run(ms_bin, argv)
    artifact = load_evidence(response["result"]["evidenceRef"]["path"])
    entries = []
    for measurement, entry in comparison_entries(artifact):
        entries.append({**entry, "scan": str(entry["scan"]),
                        "measurement": measurement})
    for entry in entries:
        found = {c.get("candidateId") for c in candidate_rows(entry)}
        if found != {cand["id"] for cand in cands}:
            raise ValueError("CANDIDATE_ID_MISMATCH: kernel results do not cover requested model ids")
    return {
        "scans": [entry["scan"] for entry in entries],
        "candidates": cands, "parsed": parsed, "comparisons": entries,
        "params": response.get("params", {}),
        "warnings": list(response.get("warnings") or []),
        "evidence": evidence_out,
    }


def for_scan(measure: dict, scan: str | None = None) -> dict:
    """One comparison entry of a measurement, by scan (default: the first)."""
    entries = measure["comparisons"]
    if scan is None:
        return entries[0]
    for entry in entries:
        if entry["scan"] == str(scan):
            return entry
    raise ValueError(f"scan {scan} is not part of this measurement")


def missing_runs(values: set[int], start: int, end: int) -> list[list[int]]:
    runs: list[list[int]] = []
    current: list[int] = []
    for value in range(start, end + 1):
        if value not in values:
            current.append(value)
        elif current:
            runs.append(current)
            current = []
    if current:
        runs.append(current)
    return runs


def series_start(ion: str, focus: int | None, length: int) -> int:
    """First fragment number the focus site makes observable.

    C-terminal series (z/x/y) count from the end, so the site is observable
    from ``length - focus + 1`` onward; N-terminal series (a/b/c) from
    ``focus`` onward.
    """
    if focus is None:
        return 1
    if ion.startswith("z") or ion in {"x", "y"}:
        return max(1, length - focus + 1)
    if ion in {"a", "b", "c"}:
        return focus
    return 1


def worst_interior_gap(all_runs: list[list[int]], last_match: int) -> dict | None:
    """The worst INTERIOR run (a match exists above it), if any.

    Tail runs beyond a series' last match are the ladder's natural reach
    limit — the same for every candidate — so they are never discrimination
    evidence.
    """
    interior = [r for r in all_runs if r[-1] < last_match]
    if not interior:
        return None
    worst = max(interior, key=len)
    return {"start": worst[0], "end": worst[-1], "length": len(worst)}


def state_counts(statuses: list[dict]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for status in statuses:
        state = status["state"]
        counts[state] = counts.get(state, 0) + 1
    return counts


def generated_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def parse_annotations(items: list[str]) -> list[dict]:
    """Parse --annotate flags into the standard `annotations` contract field.

    'mz=LABEL' -> one dashed callout at that m/z; 'mz1..mz2=LABEL' -> one
    tinted band over the region. Spectrum templates render both through the
    forgeMs2 annBands/annCallouts atoms.
    """
    out: list[dict] = []
    for raw in items or []:
        spec, sep, label = raw.partition("=")
        label = label.strip()
        if not sep or not label:
            raise ValueError(
                "--annotate expects 'mz=LABEL' or 'mz1..mz2=LABEL', "
                f"got {raw!r}")
        if ".." in spec:
            lo, _, hi = spec.partition("..")
            out.append({"kind": "band", "from": float(lo), "to": float(hi),
                        "label": label})
        else:
            out.append({"kind": "callout", "mz": float(spec), "label": label})
    return out


def finish_figure(out, payload: dict, template: str, *, figure: str, script: str,
                  params: dict, inputs: list, caption: dict, stdout: dict) -> dict:
    """data.json -> chart-forge HTML -> manifest; returns the stdout payload.

    The caller emits it: ``spectral.py`` composes one document per call
    (the card) instead of letting each figure print its own, so there is
    exactly one JSON object on stdout.
    """
    out = Path(out)
    emit_html(template, payload, out)
    write_manifest(out, figure=figure, script=script, params=params,
                   inputs=inputs, caption=caption)
    return stdout


def fail(code: str, exc: Exception) -> int:
    print(json.dumps({"error": {"code": code, "message": str(exc)}}))
    return 1
