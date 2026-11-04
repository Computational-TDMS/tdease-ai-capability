"""MS2 figures — reached through ``spectral.py``, never directly.

Two render routes over the same kernel measurement:

- :func:`render_discrimination` — the raw-mzML one-flow discrimination figure:
  annotated MS2 spectrum (matched stems coloured by ion series, labelled by ion
  and number; hover audits every peak across candidates), per-candidate error
  maps sandwiching the sequence, per-hypothesis ppm-error histograms. One model
  and several models share the figure; only the second carries a relative
  comparison. ``spectral.py measure`` passes the measurement it already took.
- :func:`render_msalign_annotation` — the deconvoluted TopFD MS2 annotation:
  the spectrum with the ppm error panel (fragment-annotation template).

BOTH domains annotate spectra; pick the route by the observation domain you
hold (raw mzML vs TopFD msalign), never by "annotation vs error map".

Matching authority is the Rust kernel in every branch: exact isotope-envelope
matching on raw mzML, charge-aware neutral-mass one-to-one on msalign (TopFD
already collapsed the envelopes; envelope machinery does not apply there).
10 ppm by default; fragment series follow the activation (ETD -> C,Z;
HCD/CID -> B,Y). Neither raw-mzML branch sends its full match table through
stdout.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from _ms2render import (
    candidate_rows, finish_figure, for_scan, fragment_match, generated_stamp,
    ion_types_label, measure_candidates, missing_runs, load_evidence, comparison_entries,
    parse_annotations, parse_proforma, series_start, state_counts,
    worst_interior_gap,
)
from _shared.native import load_json_arg

SHORT_TAG = {
    "Acetyl": "ac", "Dimethyl": "me2", "Oxidation": "ox",
    "Bromo": "Br", "Phospho": "ph", "Methyl": "me", "Trimethyl": "me3",
    "Carbamidomethyl": "cam", "Sulfo": "su",
}
CANDIDATE_COLORS = ["#2a9d8f", "#e8871e", "#7a5195", "#4c78a8"]


# ------------------------------------------------------------ discrimination

def parse_position(value) -> int | None:
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        match = re.fullmatch(r"residue:(\d+)", value)
        if match:
            return int(match.group(1))
    return None


def mod_site(mod: dict) -> int | None:
    return parse_position(mod.get("position"))


def mod_key(mod: dict) -> int | str:
    position = mod.get("position")
    if position in {"nTerm", "cTerm"}:
        return position
    site = mod_site(mod)
    if site is None:
        raise ValueError(f"unknown modification position: {position!r}")
    return site


def mod_tags(parsed: dict) -> dict:
    tags = {}
    for mod in parsed["modifications"]:
        tags.setdefault(mod_key(mod), []).append(mod["tag"])
    return {site: tuple(sorted(values)) for site, values in tags.items()}


def build_slice(cand: dict, parsed: dict, result: dict, focus_site: int | None) -> dict:
    """Figure slice from one kernel result: the match table becomes chart
    data plus per-ion-series continuity summaries."""
    matches = [
        {
            "ion": m["ion"].lower(),
            "n": int(m["number"]),
            "q": int(m["charge"]),
            "cleav": int(m["cleavageAfter"]),
            "mz": float(m["observedMz"]),
            "ppm": float(m["errorPpm"]),
            "i": float(m["intensity"]),
        }
        for m in result["matches"]
    ]
    matches.sort(key=lambda m: (m["ion"], m["q"], m["cleav"]))

    length = len(parsed["sequence"])
    by_ion: dict[str, set[int]] = {}
    for m in matches:
        by_ion.setdefault(m["ion"], set()).add(m["n"])
    focus = focus_site
    runs: dict[str, list[dict]] = {}
    worst_gap: dict[str, list[int] | None] = {}
    for ion, numbers in by_ion.items():
        all_runs = missing_runs(numbers, series_start(ion, focus, length), length - 1)
        runs[ion] = [
            {"start": r[0], "end": r[-1], "length": len(r)}
            for r in all_runs if len(r) >= 2
        ]
        last_match = max(numbers) if numbers else 0
        worst_gap[ion] = worst_interior_gap(all_runs, last_match)

    series_mean = {
        ion: series.get("meanPpm")
        for ion, series in (result.get("perSeriesError") or {}).items()
        if isinstance(series, dict) and series.get("meanPpm") is not None
    }

    counts = {ion: len(numbers) for ion, numbers in by_ion.items()}
    summary = result.get("summary") or {}
    coverage = summary.get("coverage")
    if coverage is None:
        coverage = result.get("coverage")
    measurement = {
        key: result.get(key, summary.get(key))
        for key in ("coverage", "matchedIons", "discriminatingMatched",
                    "matchedIntensity", "uniqueCleavages", "medianAbsErrorPpm",
                    "p95AbsErrorPpm")
        if result.get(key, summary.get(key)) is not None
    }
    return {
        "id": cand["id"],
        "label": cand.get("label", cand["id"]),
        "matches": matches,
        "counts": counts,
        "total": len(matches),
        "seriesMean": series_mean,
        "seriesDivergence": result.get("seriesDivergence"),
        "missingRuns": runs,
        "worstGap": worst_gap,
        "stateCounts": state_counts(result.get("fragmentStatuses", [])),
        "coverage": coverage,
        "measurement": measurement,
        "siteWindow": result.get("siteWindow"),
    }


def render_discrimination(args, ms_bin: str, measure: dict | None = None) -> int:
    """The one-flow discrimination figure for ONE scan.

    ``measure`` supplies an already-taken measurement (``spectral.py measure``
    passes the one it built the card from, so the figure and the card are the
    same observation); without it the figure takes its own.
    """
    if not args.out:
        args.out = f"{Path(args.mzml).stem}-scan{args.scan}-comparison.html"
    if measure is None:
        out = Path(args.out)
        cands = load_json_arg("--candidates", args.candidates)
        if not cands:
            raise ValueError("--candidates needs at least one model")
        measure = measure_candidates(
            ms_bin, args.mzml, [str(args.scan)], cands,
            getattr(args, "evidence_out", None) or str(out.with_name(
                out.stem + ("_compare.json" if len(cands) > 1
                            else f"_match_{cands[0]['id']}.json"))),
            args.ion_types, args.max_charge, args.tolerance_ppm, args.activation,
        )
    out = Path(args.out)
    # pairs the requested model with its own kernel summary by identity, never
    # by position: a reordered response must not attach one model's matches to
    # another model's figure
    entry = for_scan(measure, str(args.scan))
    cands = measure["candidates"]
    parsed_models = measure["parsed"]
    result_by_id = {result["candidateId"]: result for result in candidate_rows(entry)}
    slices = [(build_slice(cand, parsed, result_by_id[cand["id"]], args.focus_site),
               parsed)
              for cand, parsed in zip(cands, parsed_models)]
    observed = entry.get("pooledPeaks") or entry.get("observedPeaks")
    if observed is None:
        raise ValueError("figure needs observed peaks in the fragment evidence; remeasure an older artifact")
    peaks = [{"mz": p["mz"], "i": p["intensity"]} for p in observed]
    declared = [(m.get("precursor") or {}).get("charge") for m in entry.get("spectraMetadata") or []]
    prec_charge = max((z for z in declared if z is not None), default=None)
    resolved = measure["params"].get("ionTypes") or entry.get("ionTypes")
    resolved_ions = (ion_types_label(resolved) if resolved
                     else (args.ion_types or "auto"))
    activation = entry.get("activation") or measure["params"].get("resolvedActivation") \
        or args.activation
    effective_max_charge = entry.get("maxCharge") or measure["params"].get("maxCharge")
    evidence_files = [measure["evidence"]]
    sequence = slices[0][1]["sequence"]

    # per-hypothesis modification claims: every mod of every candidate goes to
    # its own strip; sites where the candidates differ are flagged (diff)
    per_cand_tags = [mod_tags(parsed) for _, parsed in slices]
    diff_sites = {site for tags in per_cand_tags for site in tags
                  if len({t.get(site) for t in per_cand_tags}) > 1}
    focus = args.focus_site or next(
        (s for s in sorted(diff_sites, key=str)
         if isinstance(s, int)), None)

    max_i = max((p["i"] for p in peaks), default=1.0)
    display = [p for p in peaks if p["i"] >= args.min_rel_intensity * max_i]

    # Boundary ions must follow the kernel's actual series (e.g. c/z for ETD).
    series = (resolved_ions.split(",") if isinstance(resolved_ions, str)
              else resolved_ions)
    core_targets = []
    for ion in series:
        ion = ion.strip().lower()
        if focus is not None:
            if ion in {"a", "b", "c"}:
                core_targets.append((ion, focus))
                if focus > 1:
                    core_targets.append((ion, focus - 1))
            elif ion in {"x", "y", "z"}:
                core_targets.append((ion, len(sequence) - focus + 1))
                if focus < len(sequence):
                    core_targets.append((ion, len(sequence) - focus))
        elif "nTerm" in diff_sites and ion in {"a", "b", "c"}:
            core_targets.extend([(ion, 1), (ion, 2)])

    candidates = []
    for index, (slice_, parsed) in enumerate(slices):
        mods = []
        for site, tags in mod_tags(parsed).items():
            mods.append({
                "pos": site if isinstance(site, int) else (1 if site == "nTerm" else len(sequence)),
                "tag": "/".join(SHORT_TAG.get(tag, tag) for tag in tags),
                "nTerm": site == "nTerm", "cTerm": site == "cTerm",
                "diff": site in diff_sites,
            })
        mods.sort(key=lambda m: m["pos"])
        slice_["mods"] = mods
        slice_["color"] = CANDIDATE_COLORS[index % len(CANDIDATE_COLORS)]
        slice_["siteIons"] = [
            {"ion": ion, "number": number,
             "states": [{"q": m["q"], "ppm": m["ppm"], "mz": m["mz"]}
                        for m in slice_["matches"]
                        if m["ion"] == ion and m["n"] == number]}
            for ion, number in core_targets
        ]
        candidates.append(slice_)

    payload = {
        "data": {
            "peaks": display,
            "candidates": candidates,
            "sequence": sequence,
            "focusSite": focus,
            "scan": args.scan,
            "activation": activation,
            "tol": args.tolerance_ppm,
            # the kernel's effective fragment charge ceiling (follows the
            # precursor's declared charge unless --max-charge overrides it)
            "charges": list(range(1, int(effective_max_charge or prec_charge or 0) + 1)),
            "minRelIntensity": args.min_rel_intensity,
            "precursorCharge": prec_charge,
        },
        "meta": {
            "title": args.title or f"One-flow discrimination - scan {args.scan}",
            "source": (f"raw mzML MS2 scan {args.scan} - per-candidate exact "
                       "isotope-envelope matching"),
            "thresholds": (f"{args.tolerance_ppm:g} ppm; Pearson >= {measure['params'].get('minCorrelation', 0.7):g}; "
                           f"display >= {args.min_rel_intensity:g} base peak"),
            "n": str(sum(c["total"] for c in candidates)),
            "generated": generated_stamp(),
        },
    }
    return finish_figure(
        out, payload, "candidate-discrimination",
        figure="candidate-discrimination-html",
        script=Path(__file__).name,
        params={
            "scan": args.scan,
            "candidates": [{"id": c["id"], "label": c.get("label"), "proforma": c["proforma"]}
                           for c in cands],
            "ion_types": resolved_ions,
            "activation": activation,
            "max_charge": args.max_charge,
            "tolerance_ppm": args.tolerance_ppm,
            "focus_site": focus,
            "min_rel_intensity": args.min_rel_intensity,
        },
        inputs=[
            args.mzml,
            *evidence_files,
            str(out.with_suffix(".data.json")),
        ],
        caption={
            "candidates": [
                {"id": c["id"], "total": c["total"], "counts": c["counts"],
                 "missingRuns": c["missingRuns"],
                 "siteIons": c["siteIons"],
                 **({"seriesMean": c["seriesMean"]} if c["seriesMean"] else {}),
                 **({"seriesDivergence": c["seriesDivergence"]}
                    if c["seriesDivergence"] else {})}
                for c in candidates
            ],
            "lineage": "raw MS2 exact-envelope one-flow discrimination",
        },
        stdout={"figure": str(out.resolve()), "scan": args.scan,
                "activation": activation, "ionTypes": resolved_ions,
                "maxCharge": effective_max_charge,
                "tolerancePpm": args.tolerance_ppm,
                "candidates": [{"id": c["id"], "total": c["total"],
                                "counts": c["counts"],
                                "measurement": c["measurement"],
                                "siteWindow": c["siteWindow"],
                                "seriesMean": c["seriesMean"],
                                "seriesDivergence": c["seriesDivergence"],
                                "missingRuns": c["missingRuns"],
                                "siteIons": {f"{s['ion']}{s['number']}": s["states"]
                                             for s in c["siteIons"]}}
                               for c in candidates]},
    )


# ------------------------------------------------------- msalign annotation


def render_msalign_annotation(args, ms_bin: str) -> int:
    out = Path(args.out or f"{Path(args.msalign).stem}-scan{args.scan}-annotation.html")
    evidence_out = str(out.with_suffix(".fragment-match.json"))
    activation = args.activation
    ion_types = args.ion_types or "auto"

    parsed = parse_proforma(ms_bin, args.proforma)
    sequence = parsed["sequence"]
    mods = parsed.get("modifications", [])
    tags = mod_tags(parsed)
    residues = [{"pos": i, "char": ch, "mod": "/".join(tags.get(i, ())) or None}
                for i, ch in enumerate(sequence, start=1)]
    tokens = [r["char"] + (f"[{r['mod']}]" if r["mod"] else "") for r in residues]
    if "nTerm" in tags:
        tokens[0] = "[" + "/".join(tags["nTerm"]) + "]-" + tokens[0]
    if "cTerm" in tags:
        tokens[-1] += "-[" + "/".join(tags["cTerm"]) + "]"

    match_result = fragment_match(
        ms_bin, args.msalign, args.scan, args.proforma, evidence_out,
        args.ion_types, args.max_charge, args.tolerance_ppm,
        activation=args.activation,
    )
    resolved_ions = match_result.get("_resolvedIonTypes") or ion_types
    activation = match_result.get("_resolvedActivation") or activation

    peaks = []
    for m in match_result["matches"]:
        ion = m["ion"].lower()
        number = int(m["number"])
        label = f"{ion}{number}" + (f"({m['neutralLoss']})" if m.get("neutralLoss") else "")
        fragment_sequence = (
            "".join(tokens[:number]) if ion in ("a", "b", "c")
            else "".join(tokens[-number:]) if number <= len(tokens) else "".join(tokens)
        )
        peaks.append({
            "ion": ion,
            "number": number,
            "charge": int(m["charge"]),
            "mz": float(m["observedMz"]),
            "theo": float(m["theoreticalMz"]),
            "ppm": float(m["errorPpm"]),
            "intensity": float(m["intensity"]),
            "label": label,
            "sequence": fragment_sequence,
        })

    activation_label = f" - {activation}" if activation else ""
    contract = {
        "peaks": peaks,
        "residues": residues,
        "tolerance_ppm": args.tolerance_ppm,
        "scan": f"{args.scan} (deconvoluted msalign)",
        "coverage": match_result.get("summary", {}).get("coverage"),
    }
    annotations = parse_annotations(args.annotate)
    if annotations:
        contract["annotations"] = annotations
    meta = {
        "title": f"Annotated msalign MS/MS - scan {args.scan}{activation_label}",
        "subtitle": (f"{Path(args.msalign).name} - {len(peaks)} matched deconvoluted "
                     f"fragments ({resolved_ions}, {args.tolerance_ppm:g} ppm) - "
                     "MATCHED FRAGMENTS ONLY: unmatched deconvoluted peaks are not "
                     "drawn, so the figure is not spectrum coverage - hover "
                     "cleavage marks, peaks or error points - drag to zoom, "
                     "double-click to reset"),
        "source": (f"tdeasenexus-toolbox fragment --mode match (msalign domain) - "
                   f"{Path(args.msalign).name} scan {args.scan}"),
        "thresholds": (f"fragment tolerance {args.tolerance_ppm:g} ppm; "
                       f"L={len(residues)} residues; ion types {resolved_ions}"),
        "n": str(len(peaks)),
        "generated": generated_stamp(),
    }
    return finish_figure(
        out, {"data": contract, "meta": meta}, "fragment-annotation",
        figure="msalign-fragment-annotation-html",
        script=Path(__file__).name,
        params={"scan": args.scan, "proforma": args.proforma,
                "ion_types": resolved_ions, "max_charge": args.max_charge,
                "tolerance_ppm": args.tolerance_ppm,
                "activation": activation,
                "annotate": args.annotate},
        inputs=[args.msalign, evidence_out, str(out.with_suffix(".data.json"))],
        caption={"peaks": len(peaks),
                 "residues": len(residues),
                 "lineage": "chart-forge fragment-annotation template (msalign domain)"},
        stdout={
            "figure": str(out.resolve()),
            "scan": args.scan,
            "peaks": len(peaks),
            "matchedFragments": len(peaks),
            "coverage": contract["coverage"],
        },
    )
