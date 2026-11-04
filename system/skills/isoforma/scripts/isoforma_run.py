# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "scipy"]
# ///
"""IsoForma: original mzML scan gates, summed MS2, Rust envelope matching."""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import minimize
from scipy.special import betaln

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _shared.native import binary as native_binary, mass_constants


def die(code: str, message: str) -> None:
    print(json.dumps({"error": {"code": code, "message": message}}), file=sys.stderr)
    raise SystemExit(1)


def make_kernel(explicit: str | None) -> Path:
    name = native_binary(explicit)
    path = Path(shutil.which(name) or name).resolve()
    if not path.is_file():
        die("KERNEL", f"Rust toolbox not found: {path}")
    return path


def run_json(argv: list[str | Path]) -> dict:
    completed = subprocess.run([str(value) for value in argv], capture_output=True,
                               text=True, encoding="utf-8", check=False)
    if completed.returncode:
        detail = completed.stderr.strip() or completed.stdout.strip()
        die("KERNEL", detail.splitlines()[-1] if detail else "kernel failed")
    return json.loads(completed.stdout)


def parse_sites(spec: str, sequence: str) -> list[int]:
    parts = [part.strip() for part in spec.split(",") if part.strip()]
    if len(parts) < 2:
        die("SITES", "declare residue and positions, for example K,3,9,14")
    residue = parts[0].upper()
    sites = sorted({int(value) for value in parts[1:]})
    for site in sites:
        if site < 1 or site > len(sequence) or sequence[site - 1] != residue:
            actual = sequence[site - 1] if 1 <= site <= len(sequence) else "outside sequence"
            die("SITES", f"site {site} is {actual}, not {residue}")
    return sites


def parse_fixed(spec: str) -> dict[int, str]:
    fixed = {}
    for item in filter(None, (part.strip() for part in spec.split(","))):
        key, tag = item.split(":", 1)
        fixed[0 if key.lower() in {"n-term", "nterm"} else int(key)] = tag.strip()
    return fixed


def load_isomer_model(path: str, occupancy: int, allowed_sites: list[int]) -> tuple[list[tuple[int, ...]], list[dict]]:
    by_category: dict[str, list[tuple[int, ...]]] = {}
    seen: set[tuple[int, ...]] = set()
    with Path(path).open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or not {"occupancy", "sites"}.issubset(reader.fieldnames):
            die("ISOMER_MODEL", "model CSV requires occupancy and sites columns")
        for row_number, row in enumerate(reader, start=2):
            try:
                row_occupancy = int((row.get("occupancy") or "").strip())
            except ValueError:
                die("ISOMER_MODEL", f"invalid occupancy at row {row_number}")
            if row_occupancy != occupancy:
                continue
            try:
                config = tuple(sorted(int(value) for value in
                                      (row.get("sites") or "").split("+") if value.strip()))
            except ValueError:
                die("ISOMER_MODEL", f"invalid sites at row {row_number}")
            if len(config) != occupancy or len(set(config)) != occupancy:
                die("ISOMER_MODEL", f"row {row_number} does not declare {occupancy} unique sites")
            if not set(config).issubset(allowed_sites):
                die("ISOMER_MODEL", f"row {row_number} uses undeclared sites")
            if config in seen:
                die("ISOMER_MODEL", f"duplicate configuration at row {row_number}")
            seen.add(config)
            label = (row.get("category") or "").strip() or "+".join(f"K{site}" for site in config)
            by_category.setdefault(label, []).append(config)
    if not seen:
        die("ISOMER_MODEL", f"model has no rows for occupancy {occupancy}")
    configs = [config for members in by_category.values() for config in members]
    categories = [{"id": label, "configs": members}
                  for label, members in by_category.items()]
    return configs, categories


def load_feature_scan_union(path: str) -> dict[tuple[str, int], set[str]]:
    """Collect the union of feature-linked MS2 scans per source and occupancy."""
    required = {"source_file", "occupancy", "feature_group", "scan"}
    selected: dict[tuple[str, int], set[str]] = defaultdict(set)
    with Path(path).open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            die("FEATURE_MAP", "feature map requires source_file, occupancy, feature_group, scan columns")
        for row_number, row in enumerate(reader, start=2):
            source = (row.get("source_file") or "").strip()
            scan = (row.get("scan") or "").strip()
            try:
                occupancy = int((row.get("occupancy") or "").strip())
            except ValueError:
                die("FEATURE_MAP", f"invalid occupancy at row {row_number}")
            if not source or not scan or occupancy < 0:
                die("FEATURE_MAP", f"invalid source, occupancy, or scan at row {row_number}")
            selected[(source, occupancy)].add(scan)
    if not selected:
        die("FEATURE_MAP", "feature map has no scan rows")
    return dict(selected)


def proforma(sequence: str, variable_sites: tuple[int, ...], variable_tag: str,
             fixed: dict[int, str]) -> str:
    mods = dict(fixed)
    for site in variable_sites:
        if site in mods:
            die("DECLARATION", f"fixed and variable modifications overlap at {site}")
        mods[site] = variable_tag
    nterm = mods.pop(0, None)
    output = [f"[{nterm}]-"] if nterm else []
    for site, residue in enumerate(sequence, 1):
        output.append(residue)
        if site in mods:
            output.append(f"[{mods[site]}]")
    return "".join(output)


def activation_code(value: str | None) -> str:
    text = (value or "").strip().lower()
    if "electron transfer" in text or text == "etd":
        return "etd"
    if "electron capture" in text or text == "ecd":
        return "ecd"
    if "higher-energy" in text or text in {"hcd", "beam-type collision-induced dissociation"}:
        return "hcd"
    if "collision-induced" in text or text == "cid":
        return "cid"
    return text


def list_ms2(kernel: Path, source: Path) -> list[dict]:
    result = run_json([kernel, "select", source, "--statement",
                       "SELECT spectra FROM mzml WHERE ms_level = 2 LIMIT 10000"])["result"]
    if result.get("truncated"):
        die("SOURCE", f"MS2 listing truncated for {source.name}")
    return result.get("spectra", [])


def gate_scans(kernel: Path, source: Path, spectra: list[dict], model: str,
               mass: float, rt_min: float, rt_max: float, activation: str | None,
               mass_window: float, precursor_ppm: float, min_isotopes: int,
               requested: set[str] | None, isotope_gate: str,
               proton_mass: float) -> tuple[list[str], dict]:
    candidates = []
    initial_rejections = Counter()
    for spectrum in spectra:
        scan = str(spectrum.get("scan"))
        precursor = spectrum.get("precursor") or {}
        charge = precursor.get("charge")
        selected_mz = precursor.get("selectedMz")
        rt = float(spectrum.get("retentionTime") or 0.0)
        if requested is not None and scan not in requested:
            initial_rejections["scan"] += 1
        elif not rt_min <= rt <= rt_max:
            initial_rejections["rt"] += 1
        elif activation is not None and activation_code(precursor.get("activation")) != activation:
            initial_rejections["activation"] += 1
        elif selected_mz is None or charge is None:
            # A missing charge cannot support a neutral-mass test: an invented
            # z=1 would reject real candidates as "mass" mismatches.
            initial_rejections["precursor"] += 1
        elif abs(float(selected_mz) * int(charge)
                 - proton_mass * int(charge) - mass) > mass_window:
            initial_rejections["mass"] += 1
        else:
            candidates.append(scan)

    accepted = []
    gate_rejections = Counter()
    isotope_warnings = []
    validated = [] if not candidates else run_json([
            kernel, "validate-precursor", source, "--scan", ",".join(candidates),
            "--proforma", model,
            "--tolerance-ppm", str(precursor_ppm),
            "--min-matched-peaks", str(min_isotopes),
        ])["result"]
    if isinstance(validated, dict):
        validated = [validated]
    for result in validated:
        scan = result["ms2Scan"]
        if result["isValid"]:
            accepted.append(scan)
        else:
            warning = {
                "scan": scan,
                "matchedPeaks": result["matchedPeaks"],
                "matchedPeaksInIsolationWindow": result["matchedPeaksInIsolationWindow"],
                "minimumProfileMatches": result["minimumProfileMatches"],
            }
            if isotope_gate == "diagnostic":
                accepted.append(scan)
                isotope_warnings.append(warning)
            else:
                if result["matchedPeaksInIsolationWindow"] < min_isotopes:
                    gate_rejections["isolation_isotopes"] += 1
                if result["matchedPeaks"] < result["minimumProfileMatches"]:
                    gate_rejections["profile_coverage"] += 1
    accepted.sort(key=int)
    return accepted, {
        "initialCandidates": len(candidates),
        "acceptedScans": accepted,
        "acceptedScanCount": len(accepted),
        "initialRejections": dict(initial_rejections),
        "gateRejections": dict(gate_rejections),
        "isotopeGate": isotope_gate,
        "isotopeWarnings": isotope_warnings,
    }


def entry_candidates(entry: dict) -> list[dict]:
    return ((entry.get("measurement") or {}).get("candidates")
            or (entry.get("result") or {}).get("candidates") or [])


def measure_family(kernel: Path, source: Path, configs: list[dict], scans: list[str],
                   ions: list[str] | None, activation: str | None, ppm: float,
                   max_charge: int, correlation: float, bin_ppm: float,
                   minimum_abundance_pct: float, evidence_path: Path) -> tuple[dict, dict]:
    """Run the kernel's summed-fragment comparison and KEEP its evidence file.

    The kernel only emits full match tables through ``--evidence-out`` (they
    are file-backed by contract), so the measurement is only auditable if the
    file survives; it is written next to the analysis outputs, never into a
    scratch directory. Activation/series selection stays with the kernel:
    ``--activation`` derives the series there, and with neither option the
    kernel reads the source's own acquisition metadata.
    """
    evidence_path.parent.mkdir(parents=True, exist_ok=True)
    command = [
        kernel, "fragment", source, "--mode", "compare", "--sum-scans",
        "--scan", ",".join(scans), "--max-charge", str(max_charge),
        "--tolerance-ppm", str(ppm), "--min-correlation", str(correlation),
        "--bin-ppm", str(bin_ppm),
        "--minimum-abundance-pct", str(minimum_abundance_pct),
        "--evidence-out", evidence_path, "--candidate-proformas",
        json.dumps([{"id": f"cfg{index}", "proforma": config["proforma"]}
                    for index, config in enumerate(configs)]),
    ]
    if ions:
        command += ["--ion-types", ",".join(ions)]
    if activation:
        command += ["--activation", activation]
    doc = run_json(command)
    params = doc.get("params") if isinstance(doc.get("params"), dict) else {}
    return json.loads(evidence_path.read_text(encoding="utf-8")), params


def family_equations(evidence: dict) -> list[tuple[dict[int, float], float]]:
    rows = []
    # The evidence artifact is a full CLI envelope: comparisons live under
    # result, each entry carrying its complete measurement.
    for comparison in evidence.get("result", {}).get("comparisons") or []:
        by_bond = defaultdict(dict)
        for candidate in entry_candidates(comparison):
            index = int(candidate["candidateId"].removeprefix("cfg"))
            matches = candidate.get("matches", [])
            for state in candidate.get("fragmentStatuses") or matches:
                bond = (state["ion"], state["cleavageAfter"], state["charge"])
                mass = round(float(state["theoreticalMz"]), 6)
                by_bond[bond].setdefault(
                    mass, {"configs": set(), "intensity": None, "peaks": set()}
                )["configs"].add(index)
            for match in matches:
                quant = match.get("quantification") or {}
                observed = float(quant.get("observedIntensity") or 0.0)
                fraction = float(quant.get("theoreticalFraction") or 0.0)
                if observed <= 0 or not 0 < fraction <= 1.0000001:
                    continue
                bond = (match["ion"], match["cleavageAfter"], match["charge"])
                group = by_bond[bond][round(float(match["theoreticalMz"]), 6)]
                group["intensity"] = observed / fraction
                group["peaks"].update(quant.get("observedPeakIndices", []))
        owners = defaultdict(set)
        for bond, groups in by_bond.items():
            for mass, group in groups.items():
                for peak in group["peaks"]:
                    owners[peak].add((*bond, mass))
        for groups in by_bond.values():
            for group in groups.values():
                if any(len(owners[peak]) > 1 for peak in group["peaks"]):
                    group["intensity"] = None
            present = [group for group in groups.values() if group["intensity"] is not None]
            pairs = [(left, right) for left, right in itertools.combinations(present, 2)
                     if left["peaks"].isdisjoint(right["peaks"])]
            for left, right in pairs:
                total = left["intensity"] + right["intensity"]
                row = {index: right["intensity"] / total for index in left["configs"]}
                row.update({index: -left["intensity"] / total for index in right["configs"]})
                rows.append((row, 1.0 / math.sqrt(len(pairs))))
    return rows


def ambiguity_components(null: np.ndarray, count: int) -> list[list[int]]:
    projector = null @ null.T
    pending, components = set(range(count)), []
    while pending:
        component, frontier = set(), {min(pending)}
        while frontier:
            index = frontier.pop()
            if index in component:
                continue
            component.add(index)
            frontier.update(other for other in pending
                            if abs(projector[index, other]) > 1e-8)
        pending.difference_update(component)
        components.append(sorted(component))
    return [component for component in components if len(component) > 1]


def solve_family(rows: list[tuple[dict[int, float], float]], count: int,
                 estimator: str = "beta", reporting_groups: list[dict] | None = None) -> dict:
    if count == 1:
        return {"status": "single_configuration", "classes": [[0]], "proportions": [1.0],
                "rank": 1, "nRatios": 0, "reportability": "reportable",
                "fitMethod": "none"}
    if not rows:
        return {"status": "no_discriminating_evidence", "classes": [list(range(count))],
                "rank": 1, "nRatios": 0, "reportability": "diagnostic_only",
                "fitMethod": "none"}
    raw = np.array([[row.get(index, 0.0) for index in range(count)] for row, _ in rows])
    left, right = (raw > 0).astype(float), (raw < 0).astype(float)
    ratios = -raw.min(axis=1) / (raw.max(axis=1) - raw.min(axis=1))
    weights = np.array([weight ** 2 for _, weight in rows])
    size = count

    def predictions(values):
        first, second = left @ values, right @ values
        return first / np.maximum(first + second, 1e-12)

    def jacobian(values):
        first, second = left @ values, right @ values
        return (left * second[:, None] - right * first[:, None]) / np.maximum(
            (first + second)[:, None] ** 2, 1e-24)

    constraint = {"type": "eq", "fun": lambda values: values.sum() - 1.0,
                  "jac": lambda values: np.ones(size)}
    support = {"type": "ineq", "fun": lambda values: (left + right) @ values - 1e-9}
    fit = minimize(lambda values: 0.5 * float(np.dot(
        weights, (predictions(values) - ratios) ** 2)), np.full(size, 1.0 / size),
        jac=lambda values: jacobian(values).T @ (weights * (predictions(values) - ratios)),
        bounds=[(0.0, 1.0)] * size, constraints=[constraint, support], method="SLSQP",
        options={"ftol": 1e-12, "maxiter": 1000})
    if not fit.success:
        return {"status": "solver_failed", "message": str(fit.message)}
    values, used, precision = fit.x, "least-squares", None
    beta_eligible = len(rows) >= 3 and np.all((ratios > 0) & (ratios < 1))
    beta_fit_succeeded = False
    if estimator == "beta" and beta_eligible:
        def objective(parameters):
            means = np.clip(predictions(parameters[:size]), 1e-9, 1 - 1e-9)
            phi = np.exp(parameters[-1])
            alpha, beta = means * phi, (1 - means) * phi
            return float(np.dot(weights, betaln(alpha, beta)
                - (alpha - 1) * np.log(ratios) - (beta - 1) * np.log1p(-ratios)))
        beta_fit = minimize(objective, np.r_[values, math.log(20.0)],
            bounds=[(0.0, 1.0)] * size + [(math.log(0.1), math.log(1e5))],
            constraints=[{"type": "eq", "fun": lambda p: p[:size].sum() - 1.0},
                         {"type": "ineq", "fun": lambda p: (left + right) @ p[:size] - 1e-9}],
            method="SLSQP", options={"ftol": 1e-10, "maxiter": 1000})
        if beta_fit.success:
            values, precision, used = beta_fit.x[:size], float(np.exp(beta_fit.x[-1])), "beta"
            beta_fit_succeeded = True
    predictions_at_fit = predictions(values)
    ratio_rmse = float(np.sqrt(np.average(
        (predictions_at_fit - ratios) ** 2, weights=weights)))
    fit_quality = {
        "ratioRows": len(rows),
        "interiorRatioRows": int(np.sum((ratios > 0) & (ratios < 1))),
        "betaEligible": bool(beta_eligible),
        "betaUsed": beta_fit_succeeded,
        "ratioRmse": ratio_rmse,
    }
    design = np.vstack([jacobian(values), np.ones(size)])
    null = null_space(design, rcond=1e-9)
    rank = size - null.shape[1]
    suggestions = ambiguity_components(null, size) if null.size else []
    provisional = used != "beta"
    if reporting_groups:
        group_vectors = []
        for group in reporting_groups:
            vector = np.zeros(size)
            vector[group["members"]] = 1.0
            group_vectors.append(vector)
        identifiable = [not null.size or float(np.max(np.abs(vector @ null))) <= 1e-8
                        for vector in group_vectors]
        if not all(identifiable):
            return {"status": "declared_model_unidentifiable",
                    "suggested_classes": suggestions,
                    "categoryIdentifiable": dict(zip(
                        (group["id"] for group in reporting_groups), identifiable)),
                    "rank": rank, "nRatios": len(rows), "estimator": used,
                    "fitQuality": fit_quality, "reportability": "diagnostic_only"}
        return {"status": "provisional" if provisional else "solved",
                "classes": [group["members"] for group in reporting_groups],
                "labels": [group["id"] for group in reporting_groups],
                "proportions": [float(vector @ values) for vector in group_vectors],
                "rank": rank, "nRatios": len(rows), "estimator": used,
                "betaPrecision": precision, "fitQuality": fit_quality,
                "reportability": "diagnostic_only" if provisional else "reportable"}
    if rank < size:
        return {"status": "suggested_aggregates", "suggested_classes": suggestions,
                "rank": rank, "nRatios": len(rows), "estimator": used,
                "mergeReason": "coupled abundance ambiguity", "fitQuality": fit_quality,
                "reportability": "diagnostic_only"}
    return {"status": "provisional" if provisional else "solved",
            "classes": [[index] for index in range(count)],
            "proportions": [float(value) for value in values], "rank": rank,
            "nRatios": len(rows), "estimator": used, "betaPrecision": precision,
            "fitQuality": fit_quality,
            "reportability": "diagnostic_only" if provisional else "reportable"}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=(
        "Select mzML MS2 scans with IsoForma gates, sum them, match fragment isotope "
        "envelopes in Rust, and estimate positional-isomer proportions."))
    parser.add_argument("--mzml-dir", required=True)
    parser.add_argument("--source", action="append")
    parser.add_argument("--sequence", required=True)
    parser.add_argument("--sites", required=True, help="residue,pos1,pos2,...")
    parser.add_argument("--isomer-model",
                        help="CSV with occupancy,sites,category for declared reporting groups")
    parser.add_argument("--mod-tag", required=True)
    parser.add_argument("--fixed-tags", default="")
    parser.add_argument("--occupancies", default="auto")
    parser.add_argument("--rt-min", required=True, type=float)
    parser.add_argument("--rt-max", required=True, type=float)
    parser.add_argument("--activation",
                        help="fragmentation method (ETD/ECD/HCD/CID); when omitted, "
                             "scan gating skips the activation filter and the kernel "
                             "reads each source's own acquisition metadata")
    parser.add_argument("--scans", help="optional comma-separated scan numbers")
    parser.add_argument(
        "--feature-scan-map",
        help="CSV of source_file,occupancy,feature_group,scan; union feature-linked scans per source/occupancy",
    )
    parser.add_argument("--mass-window-da", type=float, default=10.0)
    parser.add_argument("--precursor-ppm", type=float, default=15.0)
    parser.add_argument("--min-ms1-isotopes", type=int, default=3)
    parser.add_argument("--isotope-gate", choices=("strict", "diagnostic"),
                        help="strictly reject failed MS1 isotope checks, or retain them as warnings")
    parser.add_argument("--ion-types",
                        help="override the derived fragment series (comma-separated); "
                             "by default the kernel derives them from --activation "
                             "or the source acquisition metadata")
    parser.add_argument("--fragment-ppm", type=float, default=10.0)
    parser.add_argument("--max-fragment-charge", type=int, default=4)
    parser.add_argument("--min-correlation", type=float, default=0.7)
    parser.add_argument("--bin-ppm", type=float, default=5.0)
    parser.add_argument("--minimum-abundance-pct", type=float, default=0.01)
    parser.add_argument("--evidence-dir",
                        help="directory for the kernel's full fragment-match evidence "
                             "files, one per source and occupancy (default: "
                             "<mzml-dir>/isoforma-evidence)")
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--ms-bin")
    args = parser.parse_args(argv)

    if args.rt_min > args.rt_max or args.min_ms1_isotopes < 1:
        die("PARAMETERS", "invalid RT window or minimum isotope count")
    kernel = make_kernel(args.ms_bin)
    sequence = args.sequence.strip().upper()
    sites = parse_sites(args.sites, sequence)
    fixed = parse_fixed(args.fixed_tags)
    occupancies = (list(range(len(sites) + 1)) if args.occupancies.lower() == "auto"
                   else sorted({int(value) for value in args.occupancies.split(",")}))
    if any(occupancy < 0 or occupancy > len(sites) for occupancy in occupancies):
        die("OCCUPANCIES", f"occupancies must be within 0..{len(sites)}")
    activation = activation_code(args.activation) if args.activation else None
    if args.activation and activation not in {"etd", "ecd", "hcd", "cid"}:
        die("ACTIVATION", f"unknown activation {args.activation!r}; "
            "expected ETD/ECD/HCD/CID, or omit to read source metadata")
    ions = ([value.strip() for value in args.ion_types.split(",") if value.strip()]
            if args.ion_types else None)
    if args.scans and args.feature_scan_map:
        die("PARAMETERS", "use --scans or --feature-scan-map, not both")
    requested = set(args.scans.split(",")) if args.scans else None
    feature_scan_union = (load_feature_scan_union(args.feature_scan_map)
                          if args.feature_scan_map else None)
    isotope_gate = args.isotope_gate or ("diagnostic" if feature_scan_union else "strict")

    directory = Path(args.mzml_dir)
    sources = sorted(directory.glob("*.mzML"))
    if args.source:
        names = set(args.source)
        sources = [source for source in sources if source.name in names]
        if names != {source.name for source in sources}:
            die("SOURCE", "a requested source was not found")
    if not sources:
        die("SOURCE", "no mzML files found")

    proton_mass = float(mass_constants(args.ms_bin)["protonMass"])
    evidence_dir = (Path(args.evidence_dir) if args.evidence_dir
                    else directory / "isoforma-evidence")

    configurations = {}
    reporting_models = {}
    for occupancy in occupancies:
        if args.isomer_model:
            combinations, categories = load_isomer_model(
                args.isomer_model, occupancy, sites)
        else:
            combinations = list(itertools.combinations(sites, occupancy))
            categories = None
        configs = [{"sites": list(combination),
                    "proforma": proforma(sequence, combination, args.mod_tag, fixed)}
                   for combination in combinations]
        mass = run_json([kernel, "parse-proforma", "--proforma",
                         configs[0]["proforma"]])["result"]["monoisotopicMass"]
        for config in configs:
            config["mass"] = mass
        configurations[occupancy] = configs
        if categories:
            index = {tuple(config["sites"]): position
                     for position, config in enumerate(configs)}
            reporting_models[occupancy] = [
                {"id": category["id"],
                 "members": [index[tuple(config)] for config in category["configs"]]}
                for category in categories]

    results = []
    for source in sources:
        spectra = list_ms2(kernel, source)
        for occupancy in occupancies:
            configs = configurations[occupancy]
            if feature_scan_union is not None:
                requested_for_occupancy = feature_scan_union.get(
                    (source.name, occupancy), set())
                selection_basis = "feature_scan_union"
            else:
                requested_for_occupancy = requested
                selection_basis = "explicit_scan_set" if requested is not None else "rt_window"
            scans, selection = gate_scans(
                kernel, source, spectra, configs[0]["proforma"], configs[0]["mass"],
                args.rt_min, args.rt_max, activation, args.mass_window_da,
                args.precursor_ppm, args.min_ms1_isotopes, requested_for_occupancy,
                isotope_gate, proton_mass)
            result = {"source": source.name, "occupancy": occupancy, **selection,
                      "selectionBasis": selection_basis,
                      "requestedScanCount": (len(requested_for_occupancy)
                                              if requested_for_occupancy is not None else None),
                      "modelSpace": {"kind": ("declared_categories" if args.isomer_model
                                                else "complete_combination_lattice"),
                                     "candidateCount": len(configs),
                                     "sites": sites,
                                     **({"categories": [group["id"] for group in
                                                        reporting_models[occupancy]]}
                                        if args.isomer_model else {})},
                      "scanCountQc": "at_least_one" if scans else "none"}
            if not scans:
                result["status"] = "no_qualified_spectra"
                result["reportability"] = "diagnostic_only"
                result["reportabilityReason"] = "no_qualified_ms2_scans"
            elif args.preflight:
                result["status"] = "ready"
            else:
                evidence_path = evidence_dir / f"{source.stem}.occ{occupancy}.json"
                evidence, measured = measure_family(
                    kernel, source, configs, scans, ions, activation,
                    args.fragment_ppm, args.max_fragment_charge, args.min_correlation,
                    args.bin_ppm, args.minimum_abundance_pct, evidence_path)
                result["evidence"] = str(evidence_path)
                result["measured"] = {
                    "scans": scans,
                    "ionTypes": measured.get("ionTypes"),
                    "activation": measured.get("resolvedActivation") or activation,
                    "activationBasis": ("requested" if activation
                                        else ("source_metadata"
                                              if measured.get("resolvedActivation")
                                              else "unresolved")),
                }
                solved = solve_family(
                    family_equations(evidence), len(configs), "beta",
                    reporting_models.get(occupancy))
                if solved.get("suggested_classes"):
                    solved["suggestedGroups"] = [
                        [configs[index]["sites"] for index in group]
                        for group in solved.pop("suggested_classes")]
                elif solved.get("proportions"):
                    solved["proportions"] = [
                        {**({"category": solved["labels"][position]}
                            if solved.get("labels") else {}),
                         "sites": [configs[index]["sites"] for index in group],
                         "fraction": round(fraction, 6)}
                        for position, (group, fraction) in enumerate(
                            zip(solved["classes"], solved["proportions"]))]
                solved.pop("classes", None)
                solved.pop("labels", None)
                result.update(solved)
            results.append(result)

    print(json.dumps({
        "tool": "isoforma_run",
        "route": "IsoForma gates -> summed mzML -> Rust informed envelopes",
        "calibration": {"status": "unassessed",
                        "reason": "no source reference profile was supplied"},
        "parameters": {"sequence": sequence, "modTag": args.mod_tag,
                       "fixedTags": args.fixed_tags, "sites": sites,
                       "activation": args.activation, "occupancies": occupancies,
                       "rt": [args.rt_min, args.rt_max], "massWindowDa": args.mass_window_da,
                       "precursorPpm": args.precursor_ppm,
                       "minMs1Isotopes": args.min_ms1_isotopes,
                       "fragmentPpm": args.fragment_ppm,
                       "maxFragmentCharge": args.max_fragment_charge,
                       "minCorrelation": args.min_correlation,
                       "binPpm": args.bin_ppm,
                       "minimumAbundancePct": args.minimum_abundance_pct,
                       "isotopeGate": isotope_gate,
                       "featureScanMap": args.feature_scan_map,
                       "evidenceDir": str(evidence_dir),
                       "protonMass": proton_mass,
                       "protonMassSource": "tdeasenexus-toolbox schema theory"},
        "results": results,
    }, ensure_ascii=False, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
