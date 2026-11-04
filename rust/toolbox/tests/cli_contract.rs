//! Divergence tests: the `tdeasenexus-toolbox` CLI must produce **byte-identical**
//! `result` JSON to a direct core call for the same inputs.
//!
//! This is the executable proof that the CLI is an adapter, not a second
//! implementation (design D2 / D11). The comparison is on `result` only —
//! the envelope's `elapsedMs` is scheduling noise, and `params`/`source` are
//! provenance, not the answer.

use std::process::Command;

use tdeasenexus_toolbox::{
    aggregate_window_inner, anchored_envelope, calibrate_mass_shift,
    calibrate_mass_shift_from_observed, chemical_from_counts, fit_envelope, parse_proforma,
    query_candidate_comparison_mzml, query_fragment_match_mzml, query_index_mzml,
    query_select_mzml, theoretical_envelope, theoretical_fragment_list,
    theoretical_fragment_list_proforma_variants, CandidateComparisonQuery, CandidateModel,
    EnvelopeMode, FitEnvelopeQuery, FragmentMatchQuery, IonType, IsotopeEnvelopeQuery,
    MassModification, RefineMassShiftQuery, SpectrumSelector, PROTON,
};

/// Locate the `tdeasenexus-toolbox` binary. Cargo exposes `CARGO_BIN_EXE_<name>` for
/// integration tests; try both the hyphenated and mangled spellings, then fall
/// back to the build dir.
fn bin_path() -> std::path::PathBuf {
    for var in [
        "CARGO_BIN_EXE_tdeasenexus_toolbox",
        "CARGO_BIN_EXE_tdeasenexus-toolbox",
    ] {
        if let Ok(path) = std::env::var(var) {
            return path.into();
        }
    }
    let mut path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("target");
    let profile = if std::env::var("PROFILE").as_deref() == Ok("release") {
        "release"
    } else {
        "debug"
    };
    path.push(profile);
    path.push(if cfg!(windows) {
        "tdeasenexus-toolbox.exe"
    } else {
        "tdeasenexus-toolbox"
    });
    path
}

/// Run the CLI; returns `(success, parsed stdout)`.
fn run(args: &[&str]) -> (bool, serde_json::Value) {
    let out = Command::new(bin_path())
        .args(args)
        .output()
        .expect("tdeasenexus-toolbox binary runs");
    let json: serde_json::Value =
        serde_json::from_slice(&out.stdout).expect("CLI stdout is a JSON envelope");
    (out.status.success(), json)
}

/// Structural JSON comparison that tolerates ulp-level float drift.
///
/// The CLI and the test link *separately compiled* copies of the core. BRAIN's
/// polynomial center-of-mass accumulation is float-order-sensitive, so two
/// compilations of the same source can agree to within ~1 ulp on one value.
/// This comparator keeps the divergence test's real job — proving the CLI is
/// the same algorithm, not a second implementation — without being hostage to
/// cross-build float reproducibility: strings, integers, field shapes, and
/// array/object structure must match exactly; floats must match to 1e-12
/// relative (or 1e-9 absolute), far tighter than any real semantic difference
/// (a wrong PROTON or spacing constant is ~1e-3 Da).
fn assert_result_matches(cli: &serde_json::Value, core: &serde_json::Value) {
    match (cli, core) {
        (serde_json::Value::Number(a), serde_json::Value::Number(b)) => {
            if let (Some(a), Some(b)) = (a.as_f64(), b.as_f64()) {
                let scale = a.abs().max(b.abs()).max(1.0);
                let tolerance = (1e-12 * scale).max(1e-9);
                assert!(
                    (a - b).abs() <= tolerance,
                    "float mismatch: {a} vs {b} (delta {})",
                    (a - b).abs()
                );
            } else {
                assert_eq!(a, b, "non-float number mismatch");
            }
        }
        (serde_json::Value::Array(a), serde_json::Value::Array(b)) => {
            assert_eq!(a.len(), b.len(), "array length");
            for (left, right) in a.iter().zip(b.iter()) {
                assert_result_matches(left, right);
            }
        }
        (serde_json::Value::Object(a), serde_json::Value::Object(b)) => {
            assert_eq!(a.len(), b.len(), "object key count");
            for (key, left) in a.iter() {
                let right = b.get(key).unwrap_or_else(|| panic!("missing key {key}"));
                assert_result_matches(left, right);
            }
        }
        _ => assert_eq!(cli, core, "non-float value mismatch"),
    }
}

#[test]
fn cli_parse_proforma_result_matches_core() {
    let (ok, json) = run(&["parse-proforma", "--proforma", "[Bromo]-PEPTIDE"]);
    assert!(ok, "exit code 0");
    let core = parse_proforma("[Bromo]-PEPTIDE", &[]).unwrap();
    assert_result_matches(&json["result"], &serde_json::to_value(core).unwrap());
}

#[test]
fn cli_parse_proforma_toppic_flanking_matches_core() {
    // TopPIC proteoform format: `M.(S)[Acetyl]...K.K` with flanking dots and a
    // `(S)` modRange.
    let proforma = "M.(S)[Acetyl]ESSSKSSQPLASKQEKDGTEKRGRGRPRKQPPK.K";
    let (ok, json) = run(&["parse-proforma", "--proforma", proforma]);
    assert!(ok, "exit code 0");
    let core = parse_proforma(proforma, &[]).unwrap();
    assert_result_matches(&json["result"], &serde_json::to_value(core).unwrap());
    // Flanking residues are stripped; the (S) range residue is kept.
    assert_eq!(
        json["result"]["sequence"],
        "SESSSKSSQPLASKQEKDGTEKRGRGRPRKQPPK"
    );
    assert_eq!(json["result"]["modifications"][0]["tag"], "Acetyl");
    assert_eq!(json["result"]["modifications"][0]["position"], "site:S");
}

#[test]
fn cli_parse_proforma_unimod_resolves_without_db() {
    // [UNIMOD:n] now resolves through the embedded ontology even without --db.
    let (ok, json) = run(&["parse-proforma", "--proforma", "AC[UNIMOD:7]KDEF"]);
    assert!(ok, "exit code 0");
    assert_eq!(json["result"]["modifications"][0]["tag"], "UNIMOD:7");
    let delta = json["result"]["modifications"][0]["deltaMass"]
        .as_f64()
        .unwrap();
    assert!((delta - 0.984).abs() < 1e-2, "delta {delta}");
}

#[test]
fn cli_parse_proforma_with_external_mods_matches_core() {
    let (ok, json) = run(&[
        "parse-proforma",
        "--proforma",
        "[UniModStyle]-PEPTIDE",
        "--mods",
        "tests/fixtures/mods-unimod-style.json",
    ]);
    assert!(ok, "exit code 0");
    let core = parse_proforma(
        "[UniModStyle]-PEPTIDE",
        &[tdeasenexus_toolbox::ExternalModification {
            name: "UniModStyle".to_owned(),
            formula: "S -1 Se 1".to_owned(),
        }],
    )
    .unwrap();
    assert_result_matches(&json["result"], &serde_json::to_value(core).unwrap());
    assert_eq!(
        json["source"][0]["path"],
        "tests/fixtures/mods-unimod-style.json"
    );
}

#[test]
fn cli_theory_result_matches_core() {
    let (ok, json) = run(&[
        "theory",
        "--proforma",
        "[Bromo]-PEPTIDE",
        "--charge",
        "14",
        "--anchor",
        "11411.315",
    ]);
    assert!(ok, "exit code 0");
    let parsed = parse_proforma("[Bromo]-PEPTIDE", &[]).unwrap();
    let comp = chemical_from_counts(&parsed.composition);
    let envelope = theoretical_envelope(&comp, 14, 0.99, 0.001);
    let anchored = anchored_envelope(&envelope, 11411.315);
    assert_result_matches(
        &json["result"]["envelope"],
        &serde_json::to_value(envelope).unwrap(),
    );
    assert_result_matches(
        &json["result"]["anchored"],
        &serde_json::to_value(anchored).unwrap(),
    );
    assert_eq!(json["result"]["proforma"], "[Bromo]-PEPTIDE");
    assert_eq!(json["result"]["sequence"], parsed.sequence);
    assert_result_matches(
        &json["result"]["composition"],
        &serde_json::to_value(parsed.composition).unwrap(),
    );
    assert_result_matches(
        &json["result"]["modifications"],
        &serde_json::to_value(parsed.modifications).unwrap(),
    );
    assert_eq!(
        json["result"]["monoisotopicMass"].as_f64().unwrap(),
        parsed.monoisotopic_mass
    );
}

#[test]
fn cli_theory_rejects_mass_only_proforma() {
    let out = Command::new(bin_path())
        .args([
            "theory",
            "--proforma",
            "[+42.0106]-PEPTIDE",
            "--charge",
            "2",
        ])
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(1));
    assert!(out.stdout.is_empty());
    let error: serde_json::Value = serde_json::from_slice(&out.stderr).unwrap();
    assert_eq!(error["error"]["code"], "MASS_SPEC_INVALID_QUERY");
    assert!(error["error"]["message"]
        .as_str()
        .unwrap()
        .contains("mass-only"));
}

#[test]
fn cli_theory_isotope_doublet_discriminates_bromo_from_phospho() {
    // The Br-vs-Phospho isotope doublet: 79Br/81Br (~50/50) makes a
    // Bromo-carrying species' M+2/M0 near 1, while a Phospho (C/H/N/O/P-only)
    // species' M+2 is 13C-driven (~0.11 on this peptide). This is the
    // discriminator the mass-shift workflow relies on, and it requires a
    // *formula-bearing* tag — a bare `[+mass]` delta carries no element
    // information and cannot change the isotope shape.
    let ratios: Vec<f64> = ["PEPTIDE", "[Bromo]-PEPTIDE", "[Phospho]-PEPTIDE"]
        .iter()
        .map(|proforma| {
            let (ok, json) = run(&["theory", "--proforma", proforma, "--charge", "2"]);
            assert!(ok, "exit code 0 for {proforma}");
            let peaks = &json["result"]["envelope"]["peaks"];
            let m0 = peaks[0]["intensity"].as_f64().unwrap();
            let m2 = peaks[2]["intensity"].as_f64().unwrap();
            m2 / m0
        })
        .collect();
    assert!(
        ratios[1] > 0.7,
        "Bromo M+2/M0 {:.3} must be near-equal (79/81 Br)",
        ratios[1]
    );
    assert!(
        ratios[2] < 0.3,
        "Phospho M+2/M0 {:.3} must be 13C-only",
        ratios[2]
    );
    assert!(
        ratios[1] > ratios[2] * 3.0,
        "Bromo must dominate Phospho on M+2/M0: {ratios:?}"
    );
}

#[test]
fn cli_theory_named_dimethyl_keeps_nitrogen_in_the_carrier_only() {
    let (ok, json) = run(&[
        "theory",
        "--proforma",
        "[Dimethyl]-PEPTIDE",
        "--charge",
        "1",
        "--truncate-after",
        "1",
        "--ignore-below",
        "0",
    ]);
    assert!(ok, "exit code 0");
    assert_eq!(json["result"]["modifications"][0]["formula"], "C2H4");
    assert_eq!(json["result"]["chemicalFormula"], "C36H57N7O15");
    assert_eq!(
        json["result"]["modifications"][0]["deltaMass"],
        28.031300127592
    );
}

#[test]
fn cli_fragment_theoretical_series_are_method_decided() {
    // Theoretical fragment generation is decided by the fragmentation
    // method's chemistry: ETD -> plain c/z (the paper-tools convention);
    // HCD/CID -> b/y. The z•/z' hydrogen variants are no longer generated
    // and the ZDot/ZPrime tokens are refused with an instructive error.
    let (ok, json) = run(&[
        "fragment-theoretical",
        "--proforma",
        "AK",
        "--ion-types",
        "Z",
        "--max-charge",
        "1",
    ]);
    assert!(ok, "exit code 0");
    let rows = json["result"]["fragments"].as_array().unwrap();
    let z = rows
        .iter()
        .find(|r| r["ion"] == "z" && r["number"] == 1)
        .unwrap_or_else(|| panic!("no z1 row"));
    // K residue 128.094963 + 0.984016 (plain z).
    assert!((z["neutralMass"].as_f64().unwrap() - 129.078_979).abs() < 1e-4);
    let labels: Vec<&str> = rows.iter().filter_map(|r| r["ion"].as_str()).collect();
    assert!(
        !labels.contains(&"z*") && !labels.contains(&"z'"),
        "no variant rows are generated"
    );

    // Kernel contract: errors are a JSON object on STDERR, stdout stays empty.
    let out = Command::new(bin_path())
        .args([
            "fragment-theoretical",
            "--proforma",
            "AK",
            "--ion-types",
            "C,ZDOT",
            "--max-charge",
            "1",
        ])
        .output()
        .expect("tdeasenexus-toolbox binary runs");
    assert!(!out.status.success(), "ZDot is refused");
    let err: serde_json::Value =
        serde_json::from_slice(&out.stderr).expect("CLI errors are a JSON envelope on stderr");
    assert!(
        err["error"]["message"]
            .as_str()
            .unwrap_or("")
            .contains("decided"),
        "the refusal explains the method-decided rule: {err}"
    );
}

#[test]
fn cli_fragment_defaults_tolerance_and_derives_method_decided_ions() {
    let fixture = "fixtures/fragment-match-gav.mzML";
    let evidence = std::env::temp_dir().join(format!(
        "tdeasenexus-fragment-default-tolerance-{}.json",
        std::process::id()
    ));
    let evidence_text = evidence.to_string_lossy().to_string();
    let (ok, json) = run(&[
        "fragment",
        fixture,
        "--mode",
        "match",
        "--scan",
        "2",
        "--proforma",
        "GAV",
        "--evidence-out",
        evidence_text.as_str(),
        "--ion-types",
        "B,Y",
    ]);
    assert!(ok, "10 ppm is the recoverable default: {json}");
    assert_eq!(json["params"]["tolerancePpm"], 10.0);
    let _ = std::fs::remove_file(evidence);

    let out = Command::new(bin_path())
        .args([
            "fragment",
            fixture,
            "--mode",
            "match",
            "--scan",
            "2",
            "--proforma",
            "GAV",
            "--evidence-out",
            evidence_text.as_str(),
        ])
        .output()
        .expect("tdeasenexus-toolbox binary runs");
    assert!(
        out.status.success(),
        "activation metadata selects the series"
    );
    let derived: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(derived["params"]["ionTypes"], "B,Y");
}

#[test]
fn cli_fragment_compare_accepts_msalign_in_mass_domain() {
    let evidence =
        std::env::temp_dir().join(format!("msalign-comparison-{}.json", std::process::id()));
    let (ok, output) = run(&[
        "fragment",
        "fixtures/histone48-scan1164.msalign",
        "--scan",
        "1164",
        "--candidate-proformas",
        r#"[{"id":"a","proforma":"GAV"},{"id":"b","proforma":"G[Oxidation]AV"}]"#,
        "--evidence-out",
        evidence.to_str().unwrap(),
        "--ion-types",
        "C,Z",
    ]);
    assert!(ok);
    assert_eq!(
        output["params"]["observationDomain"],
        "deconvolutedNeutralMass"
    );
    assert_eq!(
        output["result"]["comparisons"][0]["measurement"]["matchingStrategy"],
        "deconvoluted_neutral_mass_matching"
    );
    let full: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&evidence).unwrap()).unwrap();
    for candidate in full["result"]["comparisons"][0]["measurement"]["candidates"]
        .as_array()
        .unwrap()
    {
        for matched in candidate["matches"].as_array().unwrap() {
            assert!(matched["envelopeCorrelation"].is_null());
        }
    }
    std::fs::remove_file(evidence).unwrap();
}

#[test]
fn cli_inspect_matches_core_on_real_fixture() {
    let fixture = "fixtures/histone48-scan1157-z14.mzML";
    assert!(
        std::path::Path::new(fixture).exists(),
        "{fixture} must exist"
    );
    let bytes = std::fs::read(fixture).unwrap();
    let core = query_index_mzml(&bytes, None).unwrap();

    let (ok, json) = run(&["inspect", fixture]);
    assert!(ok, "exit code 0");
    assert_result_matches(&json["result"], &serde_json::to_value(core).unwrap());
    assert_eq!(json["source"][0]["path"], fixture);
    assert!(json["source"][0]["digest"]["value"].as_str().unwrap().len() == 64);
}

#[test]
fn cli_parse_proforma_error_uses_stable_code_and_empty_stdout() {
    let out = Command::new(bin_path())
        .args(["parse-proforma", "--proforma", "[Unknown]-PEPTIDE"])
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(1), "runtime error exits 1");
    assert!(out.stdout.is_empty(), "stdout stays empty on error");
    let err: serde_json::Value = serde_json::from_slice(&out.stderr).unwrap();
    assert_eq!(err["error"]["code"], "PROFORMA_PARSE");
    assert!(err["error"]["message"]
        .as_str()
        .unwrap()
        .contains("unknown ProForma modification"));
}

#[test]
fn cli_schema_lists_every_command_with_its_envelope() {
    let help = Command::new(bin_path()).arg("--help").output().unwrap();
    assert!(help.status.success(), "--help exits 0");
    let help = String::from_utf8(help.stdout).unwrap();
    let commands = help
        .lines()
        .skip_while(|line| *line != "Commands:")
        .skip(1)
        .take_while(|line| *line != "Options:")
        .filter_map(|line| {
            let trimmed = line.trim();
            (!trimmed.is_empty()).then(|| trimmed.split_whitespace().next().unwrap())
        })
        .filter(|command| *command != "help")
        .collect::<Vec<_>>();

    assert!(!commands.is_empty(), "--help exposes the command list");
    for command in commands {
        let (ok, json) = run(&["schema", command]);
        assert!(ok, "public command {command} has a schema");
        let schema = &json["result"];
        assert_eq!(schema["title"], format!("tdeasenexus-toolbox {command}"));
        assert_eq!(schema["properties"]["command"]["const"], command);
        let required = schema["required"].as_array().unwrap();
        assert!(required.len() >= 5, "envelope required fields present");
    }
}

#[test]
fn cli_schema_unknown_command_is_a_clean_error() {
    let out = Command::new(bin_path())
        .args(["schema", "nope"])
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(1));
    let err: serde_json::Value = serde_json::from_slice(&out.stderr).unwrap();
    assert_eq!(err["error"]["code"], "MASS_SPEC_INVALID_QUERY");
}

#[test]
fn cli_usage_error_is_json_exit_2() {
    // parse-proforma without the required --proforma argument.
    let out = Command::new(bin_path())
        .arg("parse-proforma")
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(2), "usage error exits 2");
    let err: serde_json::Value = serde_json::from_slice(&out.stderr).unwrap();
    assert_eq!(err["error"]["code"], "USAGE");
}

// ---- stage 2: envelope / fit-envelope / refine-shift ----------------------

const BR_FIXTURE: &str = "fixtures/histone48-scan1157-z14.mzML";
const CHLAMY_FIXTURE: &str = "fixtures/chlamy-h4-scan4113-z11.mzML";

fn br_mz() -> f64 {
    (11393.314 + 14.0 * PROTON) / 14.0
}

fn br_prior() -> f64 {
    11393.314 - 11315.3778
}

#[test]
fn cli_envelope_result_matches_core() {
    let bytes = std::fs::read(BR_FIXTURE).unwrap();
    let query = IsotopeEnvelopeQuery {
        mass: Some(11393.314),
        proforma: None,
        external_mods: Vec::new(),
        charge: 14,
        rt_min: -1.0,
        rt_max: 1.0,
        tolerance_ppm: 100.0,
        isotopes: 10,
        max_offset: 0,
        theory_locked: false,
    };
    let (peaks, scan_count) = aggregate_window_inner(&bytes, &query, 11393.314).unwrap();

    let (ok, json) = run(&[
        "envelope",
        BR_FIXTURE,
        "--mass",
        "11393.314",
        "--charge",
        "14",
        "--tolerance-ppm",
        "100",
        "--isotopes",
        "10",
    ]);
    assert!(ok, "exit code 0");
    let expected = serde_json::json!({
        "peaks": peaks, "scanCount": scan_count, "charge": 14, "mass": 11393.314
    });
    assert_result_matches(&json["result"], &expected);
}

#[test]
fn cli_envelope_windows_match_single_window_runs() {
    // One mzML traversal must reproduce per-window results identical to the
    // single-window command (the ladder's batch contract).
    let windows = r#"[{"label":"a","mass":11393.314,"charge":14,"rtMin":-1.0,"rtMax":1.0},
                      {"label":"b","mass":11393.314,"charge":14,"rtMin":-1.0,"rtMax":1.0}]"#;
    let (ok, batch) = run(&[
        "envelope-windows",
        BR_FIXTURE,
        "--windows",
        windows,
        "--tolerance-ppm",
        "100",
        "--isotopes",
        "10",
    ]);
    assert!(ok, "exit code 0");
    let (ok, single) = run(&[
        "envelope",
        BR_FIXTURE,
        "--mass",
        "11393.314",
        "--charge",
        "14",
        "--tolerance-ppm",
        "100",
        "--isotopes",
        "10",
    ]);
    assert!(ok, "exit code 0");
    let windows_out = batch["result"]["windows"]
        .as_array()
        .expect("windows array");
    assert_eq!(windows_out.len(), 2);
    assert_eq!(windows_out[0]["label"], "a");
    assert_eq!(windows_out[1]["label"], "b");
    for window in windows_out {
        assert_result_matches(
            &serde_json::json!({
                "peaks": window["peaks"], "scanCount": window["scanCount"],
                "charge": window["charge"], "mass": window["mass"]
            }),
            &single["result"],
        );
    }
}

#[test]
fn cli_refine_shift_result_matches_core() {
    let bytes = std::fs::read(BR_FIXTURE).unwrap();
    let query = RefineMassShiftQuery {
        proforma: None,
        external_mods: Vec::new(),
        sequence_mass: Some(11315.3778),
        modification_formula: Some("Br-H".to_owned()),
        modification_mass: Some(77.9105),
        mz: br_mz(),
        charge: 14,
        rt_min: -1.0,
        rt_max: 1.0,
        tolerance_ppm: 100.0,
        prior_shift: br_prior(),
        match_window_ppm: 30.0,
    };
    let core = calibrate_mass_shift(&bytes, &query).unwrap();

    let (ok, json) = run(&[
        "refine-shift",
        BR_FIXTURE,
        "--sequence-mass",
        "11315.3778",
        "--mod-formula",
        "Br-H",
        "--mod-mass",
        "77.9105",
        "--mz",
        &br_mz().to_string(),
        "--charge",
        "14",
        "--tolerance-ppm",
        "100",
        "--prior-shift",
        &br_prior().to_string(),
    ]);
    assert!(ok, "exit code 0");
    assert_result_matches(&json["result"], &serde_json::to_value(core).unwrap());
}

#[test]
fn cli_refine_shift_observed_reuse_matches_core() {
    let bytes = std::fs::read(BR_FIXTURE).unwrap();
    let query = IsotopeEnvelopeQuery {
        mass: Some(11393.314),
        proforma: None,
        external_mods: Vec::new(),
        charge: 14,
        rt_min: -1.0,
        rt_max: 1.0,
        tolerance_ppm: 100.0,
        isotopes: 10,
        max_offset: 0,
        theory_locked: false,
    };
    let (peaks, scan_count) = aggregate_window_inner(&bytes, &query, 11393.314).unwrap();

    // The observed envelope written as the CLI envelope-result shape.
    let dir = std::env::temp_dir();
    let obs_path = dir.join(format!("tdease-refine-obs-{}.json", std::process::id()));
    let obs_json = serde_json::json!({ "peaks": peaks, "scanCount": scan_count });
    std::fs::write(&obs_path, serde_json::to_vec(&obs_json).unwrap()).unwrap();

    let refine = RefineMassShiftQuery {
        proforma: None,
        external_mods: Vec::new(),
        sequence_mass: Some(11315.3778),
        modification_formula: Some("Br-H".to_owned()),
        modification_mass: Some(77.9105),
        mz: br_mz(),
        charge: 14,
        rt_min: -1.0,
        rt_max: 1.0,
        tolerance_ppm: 100.0,
        prior_shift: br_prior(),
        match_window_ppm: 30.0,
    };
    let core = calibrate_mass_shift_from_observed(&peaks, scan_count, &refine).unwrap();

    let (ok, json) = run(&[
        "refine-shift",
        BR_FIXTURE,
        "--observed",
        obs_path.to_str().unwrap(),
        "--sequence-mass",
        "11315.3778",
        "--mod-formula",
        "Br-H",
        "--mod-mass",
        "77.9105",
        "--mz",
        &br_mz().to_string(),
        "--charge",
        "14",
        "--tolerance-ppm",
        "100",
        "--prior-shift",
        &br_prior().to_string(),
    ]);
    let _ = std::fs::remove_file(&obs_path);
    assert!(ok, "exit code 0");
    assert_result_matches(&json["result"], &serde_json::to_value(core).unwrap());
}

#[test]
fn cli_fit_envelope_matches_core() {
    let bytes = std::fs::read(BR_FIXTURE).unwrap();
    let q = IsotopeEnvelopeQuery {
        mass: Some(11393.314),
        proforma: None,
        external_mods: Vec::new(),
        charge: 14,
        rt_min: -1.0,
        rt_max: 1.0,
        tolerance_ppm: 100.0,
        isotopes: 10,
        max_offset: 0,
        theory_locked: false,
    };
    let (peaks, _) = aggregate_window_inner(&bytes, &q, 11393.314).unwrap();
    let parsed = parse_proforma("[Bromo]-PEPTIDE", &[]).unwrap();
    let comp = chemical_from_counts(&parsed.composition);
    let env = theoretical_envelope(&comp, 14, 0.99, 0.001);
    let anchored = anchored_envelope(&env, 11393.314);

    let dir = std::env::temp_dir();
    let obs_path = dir.join(format!("tdease-fit-obs-{}.json", std::process::id()));
    let theo_path = dir.join(format!("tdease-fit-theo-{}.json", std::process::id()));
    std::fs::write(&obs_path, serde_json::to_vec(&peaks).unwrap()).unwrap();
    std::fs::write(&theo_path, serde_json::to_vec(&anchored).unwrap()).unwrap();

    let fq = FitEnvelopeQuery {
        minimum_window_mz: tdeasenexus_toolbox::DEFAULT_ENVELOPE_WINDOW_MZ,
        match_window_ppm: 30.0,
    };
    let core_evaluation = fit_envelope(&peaks, &anchored, &fq).unwrap();

    let (ok, json) = run(&[
        "fit-envelope",
        "--observed",
        obs_path.to_str().unwrap(),
        "--theoretical",
        theo_path.to_str().unwrap(),
    ]);
    let _ = std::fs::remove_file(&obs_path);
    let _ = std::fs::remove_file(&theo_path);
    assert!(ok, "exit code 0");
    assert_result_matches(
        &json["result"]["evaluation"],
        &serde_json::to_value(core_evaluation).unwrap(),
    );
}

// ---- fixture regressions through the CLI ----

#[test]
fn cli_refine_shift_converges_on_bromine_fixture() {
    let (ok, json) = run(&[
        "refine-shift",
        BR_FIXTURE,
        "--sequence-mass",
        "11315.3778",
        "--mod-formula",
        "Br-H",
        "--mod-mass",
        "77.9105",
        "--mz",
        &br_mz().to_string(),
        "--charge",
        "14",
        "--prior-shift",
        &br_prior().to_string(),
    ]);
    assert!(ok, "exit code 0");
    let r = &json["result"];
    let refined = r["refinedMass"].as_f64().unwrap();
    assert!(
        (refined - 11393.314).abs() < 0.05,
        "Br refined mass {refined} must land on feature mass 11393.314"
    );
    assert_eq!(
        r["offset"], 0,
        "0.06 Da prior error must not be a whole-neutron re-assignment"
    );
    assert!(
        r["matchedPeaks"].as_u64().unwrap() >= 10,
        "enough peaks matched"
    );
    // On clean data the two leading offsets are separated by a whole neutron.
    assert!(
        r["ambiguity"].is_null(),
        "no false ambiguity: {:?}",
        r["ambiguity"]
    );
}

#[test]
fn cli_refine_shift_is_stable_across_processes() {
    let args = [
        "refine-shift",
        BR_FIXTURE,
        "--sequence-mass",
        "11315.3778",
        "--mod-formula",
        "Br-H",
        "--mod-mass",
        "77.9105",
        "--mz",
        "814.772354181069",
        "--charge",
        "14",
        "--prior-shift",
        "77.9362",
    ];
    let signatures = (0..64)
        .map(|_| {
            let (ok, json) = run(&args);
            assert!(ok, "exit code 0");
            let result = &json["result"];
            serde_json::json!({
                "observed": result["observed"].as_array().map(Vec::len),
                "theoretical": result["theoretical"].as_array().map(Vec::len),
                "matchedPeaks": result["matchedPeaks"],
                "refinedMass": result["refinedMass"],
            })
        })
        .collect::<Vec<_>>();
    assert!(
        signatures
            .iter()
            .all(|signature| signature == &signatures[0]),
        "same input changed across CLI processes: {signatures:?}"
    );
}

#[test]
fn cli_refine_shift_converges_on_chlamy_fixture() {
    let true_shift = 11411.315 - 11319.3476;
    let mz = (11411.315 + 11.0 * PROTON) / 11.0;
    let (ok, json) = run(&[
        "refine-shift",
        CHLAMY_FIXTURE,
        "--sequence-mass",
        "11319.3476",
        "--mod-mass",
        &true_shift.to_string(),
        "--mz",
        &mz.to_string(),
        "--charge",
        "11",
        "--prior-shift",
        &(true_shift + 0.04).to_string(),
    ]);
    assert!(ok, "exit code 0");
    let r = &json["result"];
    let refined = r["refinedMass"].as_f64().unwrap();
    assert!(
        (refined - 11411.315).abs() < 0.05,
        "chlamy refined mass {refined} must converge to TopPIC feature 11411.315"
    );
    assert!(r["matchedPeaks"].as_u64().unwrap() >= 3);
}

// ---- UniMod reference query (D12-R): CLI result == lib kernel result -------
//
// The divergence proof for the unimod commands: the CLI reads the real bundled
// database (rusqlite) and feeds the lib kernel, so `unimod-search` / `unimod-solve`
// must be byte-identical (mod floats) to a direct kernel call. The `source`
// provenance key is injected by the CLI only, so it is stripped before the
// comparison.

fn unimod_db_path() -> std::path::PathBuf {
    let mut path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    // rust/toolbox -> toolkit root is two levels up (toolbox -> rust -> root),
    // and system/ is a sibling of rust/ within the toolkit.
    path.push("../../system/reference/unimod.sqlite");
    path
}

/// Load the same in-memory records the CLI builds, via rusqlite.
fn load_unimod_records(
    db: &rusqlite::Connection,
) -> (
    Vec<tdeasenexus_toolbox::UniModRecord>,
    Vec<tdeasenexus_toolbox::UniModNameEntry>,
    Vec<tdeasenexus_toolbox::UniModSpecificity>,
    Vec<tdeasenexus_toolbox::UniModNeutralLoss>,
) {
    use tdeasenexus_toolbox::{
        NameKind, UniModNameEntry, UniModNeutralLoss, UniModRecord, UniModSpecificity,
    };
    let mut stmt = db.prepare("SELECT record_id, title, full_name, approved, delta_mono, delta_average, composition FROM modifications ORDER BY record_id").unwrap();
    let records: Vec<UniModRecord> = stmt
        .query_map([], |row| {
            Ok(UniModRecord {
                record_id: row.get(0)?,
                title: row.get(1)?,
                full_name: row.get(2)?,
                approved: row.get::<_, i64>(3)? != 0,
                delta_mono: row.get(4)?,
                delta_average: row.get(5)?,
                composition: row.get(6)?,
                source: "unimod".to_owned(),
            })
        })
        .unwrap()
        .collect::<Result<_, _>>()
        .unwrap();
    let mut stmt = db
        .prepare("SELECT modification_id, name, kind FROM modification_names")
        .unwrap();
    let names: Vec<UniModNameEntry> = stmt
        .query_map([], |row| {
            let kind: String = row.get(2)?;
            Ok(UniModNameEntry {
                modification_id: row.get(0)?,
                name: row.get(1)?,
                kind: match kind.as_str() {
                    "title" => NameKind::Title,
                    "full_name" => NameKind::FullName,
                    _ => NameKind::AltName,
                },
            })
        })
        .unwrap()
        .collect::<Result<_, _>>()
        .unwrap();
    let mut stmt = db.prepare("SELECT id, modification_id, site, position, classification, spec_group_id, hidden FROM specificities").unwrap();
    let specificities: Vec<UniModSpecificity> = stmt
        .query_map([], |row| {
            Ok(UniModSpecificity {
                id: row.get(0)?,
                modification_id: row.get(1)?,
                site: row.get(2)?,
                position: row.get(3)?,
                classification: row.get(4)?,
                spec_group_id: row.get(5)?,
                hidden: row.get::<_, i64>(6)? != 0,
            })
        })
        .unwrap()
        .collect::<Result<_, _>>()
        .unwrap();
    let mut stmt = db
        .prepare("SELECT specificity_id, mono_mass, average_mass, composition FROM neutral_losses")
        .unwrap();
    let neutral_losses: Vec<UniModNeutralLoss> = stmt
        .query_map([], |row| {
            Ok(UniModNeutralLoss {
                specificity_id: row.get(0)?,
                mono_mass: row.get(1)?,
                average_mass: row.get(2)?,
                composition: row.get(3)?,
            })
        })
        .unwrap()
        .collect::<Result<_, _>>()
        .unwrap();
    (records, names, specificities, neutral_losses)
}

/// Drop the `source` provenance key the CLI injects into `result`.
fn strip_source(mut value: serde_json::Value) -> serde_json::Value {
    if let serde_json::Value::Object(map) = &mut value {
        map.remove("source");
    }
    value
}

#[test]
fn cli_unimod_solve_matches_kernel_on_real_database() {
    let db_path = unimod_db_path();
    assert!(db_path.exists(), "{db_path:?} must exist");
    let conn =
        rusqlite::Connection::open_with_flags(&db_path, rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY)
            .unwrap();
    let (records, _names, specificities, neutral_losses) = load_unimod_records(&conn);
    let query = tdeasenexus_toolbox::SolveQuery {
        target_mass: 159.932662, // 2 × Phospho
        tolerance: tdeasenexus_toolbox::Tol {
            unit: tdeasenexus_toolbox::TolUnit::Ppm,
            value: 20.0,
        },
        max_modifications: 3,
        max_distinct_modifications: None,
        allow_repeated_modification: true,
        approved: None,
        allowed_sites: Vec::new(),
        fixed_modifications: Vec::new(),
        limit: 50,
    };
    let core =
        tdeasenexus_toolbox::solve_combinations(&records, &specificities, &neutral_losses, &query)
            .expect("kernel solve succeeds");

    let (ok, json) = run(&[
        "unimod-solve",
        "--db",
        db_path.to_str().unwrap(),
        "--target-mass",
        "159.932662",
        "--tolerance",
        "20",
        "--tolerance-unit",
        "ppm",
        "--max-modifications",
        "3",
        // The kernel query has no element whitelist and no classification
        // exclusions; disable the CLI's default gates so both sides solve
        // the identical record pool.
        "--element-whitelist",
        "",
        "--exclude-classification",
        "",
        "--limit",
        "50",
    ]);
    assert!(ok, "exit code 0");
    assert_result_matches(
        &strip_source(json["result"].clone()),
        &serde_json::to_value(core).unwrap(),
    );
}

#[test]
fn cli_unimod_solve_uses_sequence_sites_to_find_the_feasible_combination() {
    let db_path = unimod_db_path();
    let (ok, json) = run(&[
        "unimod-solve",
        "--db",
        db_path.to_str().unwrap(),
        "--target-mass",
        "317.3935",
        "--tolerance",
        "0.02",
        "--tolerance-unit",
        "da",
        "--max-modifications",
        "3",
        "--max-distinct",
        "3",
        "--allowed-site",
        "M:1",
        "--allowed-site",
        "G:2",
        "--allowed-site",
        "C:4",
        "--limit",
        "20",
    ]);
    assert!(ok, "exit code 0");
    let candidates = json["result"]["candidates"].as_array().unwrap();
    assert_eq!(
        candidates.len(),
        1,
        "sequence sites should remove mass-only combination noise"
    );
    assert_eq!(candidates[0]["siteFeasible"], true);
    let titles = candidates[0]["modifications"]
        .as_array()
        .unwrap()
        .iter()
        .map(|item| item["title"].as_str().unwrap())
        .collect::<std::collections::HashSet<_>>();
    assert_eq!(
        titles,
        std::collections::HashSet::from(["Met-loss", "Myristoyl", "Palmitoyl"])
    );
}

#[test]
fn cli_unimod_solve_pins_fixed_modifications_to_search_the_residual() {
    let db_path = unimod_db_path();
    let (ok, json) = run(&[
        "unimod-solve",
        "--db",
        db_path.to_str().unwrap(),
        // Carbamidomethyl (57.021464) + Phospho (79.966331): the pin leaves
        // only Phospho for the variable search.
        "--target-mass",
        "136.987795",
        "--tolerance",
        "0.02",
        "--tolerance-unit",
        "da",
        "--max-modifications",
        "2",
        "--fixed-mod",
        "Carbamidomethyl",
        "--limit",
        "20",
    ]);
    assert!(ok, "exit code 0");
    let candidates = json["result"]["candidates"].as_array().unwrap();
    assert!(
        !candidates.is_empty(),
        "the pinned CAM leaves Phospho (and near-neighbours) for the residual"
    );
    // Every candidate carries the pinned member; the best is CAM + Phospho.
    assert_eq!(candidates[0]["rank"], 1);
    assert_eq!(candidates[0]["absoluteErrorDa"], 0.0);
    let best_mods = candidates[0]["modifications"].as_array().unwrap();
    let cam = best_mods
        .iter()
        .find(|m| m["title"] == "Carbamidomethyl")
        .expect("pinned member reported");
    assert_eq!(cam["fixed"], true);
    assert_eq!(cam["count"], 1);
    let phospho = best_mods
        .iter()
        .find(|m| m["title"] == "Phospho")
        .expect("variable member reported");
    assert_eq!(phospho["fixed"], false);
    for candidate in candidates {
        let mods = candidate["modifications"].as_array().unwrap();
        assert!(
            mods.iter()
                .any(|m| m["title"] == "Carbamidomethyl" && m["fixed"] == true),
            "every candidate includes the pinned CAM: {candidate:?}"
        );
    }
    assert_eq!(json["params"]["fixedMods"][0], "Carbamidomethyl");

    // An unknown pin is a query error, not an empty answer.
    let out = Command::new(bin_path())
        .args([
            "unimod-solve",
            "--db",
            db_path.to_str().unwrap(),
            "--target-mass",
            "136.987795",
            "--tolerance",
            "0.02",
            "--fixed-mod",
            "NoSuchModification",
        ])
        .output()
        .expect("tdeasenexus-toolbox binary runs");
    assert!(!out.status.success());
    let err: serde_json::Value =
        serde_json::from_slice(&out.stderr).expect("CLI domain errors are JSON on stderr");
    let message = err["error"]["message"].as_str().unwrap_or("");
    assert!(
        message.contains("NoSuchModification") && message.contains("unknown modification"),
        "error names the unknown pin: {message}"
    );
}

#[test]
fn cli_unimod_solve_window_keeps_each_parsimony_levels_best() {
    // Regression (2026-10-17 case evidence): solving 448.4311 ± 0.03 Da with
    // max-modifications 3, the ppm-first window of 8 was entirely 3-component
    // sets (labels, detergents, substitutions, all < 0.3 mDa error) while the
    // two-component Myristoyl + Palmitoyl answer (3.07 mDa) fell outside the
    // returned window — the reader concluded "no 1-2 component hit". Deeper
    // combinations have more free parameters and fit any mass more tightly,
    // so the window guarantees the ppm-best of every distinct-modification
    // level a leading slot; rank stays the pure ppm order.
    let db_path = unimod_db_path();
    let (ok, json) = run(&[
        "unimod-solve",
        "--db",
        db_path.to_str().unwrap(),
        "--target-mass",
        "448.4311",
        "--tolerance",
        "0.03",
        "--tolerance-unit",
        "da",
        "--max-modifications",
        "3",
        "--limit",
        "8",
    ]);
    assert!(ok, "exit code 0");
    let candidates = json["result"]["candidates"].as_array().unwrap();
    assert!(
        !candidates.is_empty(),
        "the target has matches at several parsimony levels"
    );
    let titles = |c: &serde_json::Value| -> Vec<String> {
        c["modifications"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|m| m["fixed"].as_bool().unwrap_or(false) == false)
            .map(|m| m["title"].as_str().unwrap_or("").to_string())
            .collect()
    };
    let pair = candidates
        .iter()
        .find(|c| {
            let t = titles(c);
            t.len() == 2
                && t.iter().any(|t| t == "Myristoyl")
                && t.iter().any(|t| t == "Palmitoyl")
        })
        .expect("the two-component Myristoyl+Palmitoyl answer stays inside the window");
    assert!(
        pair["distinctModificationCount"].as_u64().unwrap() <= 2,
        "the parsimonious answer carries its own level, not a deeper one"
    );
    // The window leads with each present level's ppm-best, fewest components
    // first; rank remains the ppm order (no rank-1 tie with the leader).
    let leader = &candidates[0];
    assert!(
        leader["distinctModificationCount"].as_u64().unwrap()
            <= pair["distinctModificationCount"].as_u64().unwrap(),
        "fewest-component level's best leads the window"
    );
    let ranks: Vec<u64> = candidates
        .iter()
        .filter_map(|c| c["rank"].as_u64())
        .collect();
    let mut sorted = ranks.clone();
    sorted.sort_unstable();
    assert_eq!(ranks.len(), sorted.len(), "ranks are unique");
}

#[test]
fn cli_unimod_commands_accept_negative_mass_values() {
    let db_path = unimod_db_path();
    let (search_ok, _) = run(&[
        "unimod-search",
        "--db",
        db_path.to_str().unwrap(),
        "--mass",
        "-17.0265",
        "--tolerance",
        "0.02",
        "--tolerance-unit",
        "da",
    ]);
    assert!(search_ok, "negative search mass must parse");

    let (solve_ok, _) = run(&[
        "unimod-solve",
        "--db",
        db_path.to_str().unwrap(),
        "--target-mass",
        "-17.0265",
        "--tolerance",
        "0.02",
        "--tolerance-unit",
        "da",
        "--max-modifications",
        "1",
        "--allowed-site",
        "Q:1",
    ]);
    assert!(solve_ok, "negative solve mass must parse");
}

#[test]
fn cli_unimod_search_matches_kernel_on_real_database() {
    let db_path = unimod_db_path();
    assert!(db_path.exists(), "{db_path:?} must exist");
    let conn =
        rusqlite::Connection::open_with_flags(&db_path, rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY)
            .unwrap();
    let (records, names, specificities, _neutral_losses) = load_unimod_records(&conn);

    let query = tdeasenexus_toolbox::SearchQuery {
        text: Some("phosph".to_owned()),
        mass: Some(77.9105),
        tolerance: Some(tdeasenexus_toolbox::Tol {
            unit: tdeasenexus_toolbox::TolUnit::Da,
            value: 0.05,
        }),
        filters: tdeasenexus_toolbox::UniModFilters {
            residues: Vec::new(),
            approved: None,
            contains: vec!["Br".to_owned()],
            excludes: Vec::new(),
            keep_classification: Vec::new(),
        },
        limit: 100,
    };
    let name_matches =
        tdeasenexus_toolbox::search_by_name(&records, &names, &specificities, &query);
    let mass_matches = tdeasenexus_toolbox::search_by_mass(&records, &specificities, &query);
    let core = serde_json::json!({ "matches": name_matches, "massMatches": mass_matches });

    let (ok, json) = run(&[
        "unimod-search",
        "--db",
        db_path.to_str().unwrap(),
        "--text",
        "phosph",
        "--mass",
        "77.9105",
        "--tolerance",
        "0.05",
        "--tolerance-unit",
        "da",
        "--contains",
        "Br",
        // Same query as the kernel side (limit 100), not the interactive
        // default of 10.
        "--limit",
        "100",
    ]);
    assert!(ok, "exit code 0");
    assert_result_matches(&strip_source(json["result"].clone()), &core);
}

#[test]
fn cli_unimod_batch_search_matches_individual_searches() {
    let db_path = unimod_db_path();
    let input = std::env::temp_dir().join(format!(
        "tdease-unimod-batch-{}-{}.json",
        std::process::id(),
        std::thread::current().name().unwrap_or("test")
    ));
    std::fs::write(
        &input,
        r#"[{"id":"positive","mass":15.9949,"tolerance":0.05},{"id":"negative","mass":-17.0265,"tolerance":0.02}]"#,
    )
    .unwrap();

    let (ok, batch) = run(&[
        "unimod-batch-search",
        "--input",
        input.to_str().unwrap(),
        "--db",
        db_path.to_str().unwrap(),
        "--tolerance-unit",
        "da",
    ]);
    assert!(ok, "batch search exit 0");
    let queries = batch["result"]["queries"].as_array().unwrap();
    assert_eq!(queries.len(), 2);

    for (id, mass, tolerance) in [
        ("positive", "15.9949", "0.05"),
        ("negative", "-17.0265", "0.02"),
    ] {
        let (single_ok, single) = run(&[
            "unimod-search",
            "--db",
            db_path.to_str().unwrap(),
            "--mass",
            mass,
            "--tolerance",
            tolerance,
            "--tolerance-unit",
            "da",
        ]);
        assert!(single_ok);
        let batch_matches = &queries.iter().find(|query| query["id"] == id).unwrap()["massMatches"];
        assert_result_matches(batch_matches, &single["result"]["massMatches"]);
    }
    let _ = std::fs::remove_file(input);
}

#[test]
fn cli_unimod_batch_solve_preserves_single_query_constraints() {
    let db_path = unimod_db_path();
    let input = std::env::temp_dir().join(format!(
        "tdease-unimod-batch-solve-{}.json",
        std::process::id()
    ));
    std::fs::write(
        &input,
        r#"[{"id":"soybean","mass":317.3935,"tolerance":0.02}]"#,
    )
    .unwrap();
    let common = [
        "--db",
        db_path.to_str().unwrap(),
        "--tolerance-unit",
        "da",
        "--max-modifications",
        "3",
        "--max-distinct",
        "3",
        "--allowed-site",
        "M:1",
        "--allowed-site",
        "G:2",
        "--allowed-site",
        "C:4",
    ];
    let mut batch_args = vec!["unimod-batch-solve", "--input", input.to_str().unwrap()];
    batch_args.extend(common);
    let (batch_ok, batch) = run(&batch_args);
    assert!(batch_ok);

    let mut single_args = vec![
        "unimod-solve",
        "--target-mass",
        "317.3935",
        "--tolerance",
        "0.02",
    ];
    single_args.extend(common);
    let (single_ok, single) = run(&single_args);
    assert!(single_ok);
    assert_result_matches(
        &batch["result"]["queries"][0]["solution"],
        &strip_source(single["result"].clone()),
    );
    let _ = std::fs::remove_file(input);
}

#[test]
fn cli_rejects_unknown_tolerance_unit() {
    let output = Command::new(bin_path())
        .args([
            "unimod-search",
            "--mass",
            "15.9949",
            "--tolerance",
            "0.01",
            "--tolerance-unit",
            "pmm",
        ])
        .output()
        .expect("tdeasenexus-toolbox binary runs");
    assert!(!output.status.success(), "an unknown unit must fail");
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("invalid value 'pmm'"),
        "Clap must report the rejected unit"
    );
}

#[test]
fn cli_parse_proforma_with_db_matches_resolver_on_real_database() {
    let db_path = unimod_db_path();
    assert!(db_path.exists(), "{db_path:?} must exist");
    let conn =
        rusqlite::Connection::open_with_flags(&db_path, rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY)
            .unwrap();
    let (records, names, _specificities, _neutral_losses) = load_unimod_records(&conn);
    let amino_acids = {
        let mut stmt = conn
            .prepare(
                "SELECT title, three_letter, full_name, mono_mass, average_mass FROM amino_acids",
            )
            .unwrap();
        stmt.query_map([], |row| {
            Ok(tdeasenexus_toolbox::UniModAminoAcid {
                title: row.get(0)?,
                three_letter: row.get(1)?,
                full_name: row.get(2)?,
                mono_mass: row.get(3)?,
                average_mass: row.get(4)?,
            })
        })
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap()
    };
    let resolver = tdeasenexus_toolbox::UniModProFormaResolver::new(&records, &names, &amino_acids);
    let core = resolver.resolve("[Acetyl]-AC[UNIMOD:35]").unwrap();

    let (ok, json) = run(&[
        "parse-proforma",
        "--db",
        db_path.to_str().unwrap(),
        "--proforma",
        "[Acetyl]-AC[UNIMOD:35]",
    ]);
    assert!(ok, "exit code 0");
    assert_result_matches(
        &strip_source(json["result"].clone()),
        &serde_json::to_value(core).unwrap(),
    );
}

#[test]
fn cli_dual_source_user_mod_is_searchable_and_tagged() {
    let db_path = unimod_db_path();
    // A throwaway user store under temp.
    let dir = std::env::temp_dir().join(format!("tdease-usermods-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let user_db = dir.join("user-mods.sqlite");
    let _ = std::fs::remove_file(&user_db);

    // Add a user modification.
    let (ok, json) = run(&[
        "unimod-add",
        "--user-db",
        user_db.to_str().unwrap(),
        "--title",
        "MyCustomMod",
        "--formula",
        "Br H -1",
        "--mono-mass",
        "77.9105",
        "--sites",
        "K,S",
    ]);
    assert!(ok, "unimod-add exit 0");
    assert_eq!(json["result"]["added"], true);
    assert_eq!(json["result"]["unimodId"], "USER:1");

    // Dual-source mass search: the user mod and the official Bromo both carry
    // Br and match the same delta; the user record is tagged source/custom.
    let (ok, json) = run(&[
        "unimod-search",
        "--db",
        db_path.to_str().unwrap(),
        "--user-db",
        user_db.to_str().unwrap(),
        "--mass",
        "77.9105",
        "--tolerance",
        "0.01",
        "--tolerance-unit",
        "da",
        "--contains",
        "Br",
    ]);
    assert!(ok, "unimod-search exit 0");
    let mass_matches = json["result"]["massMatches"].as_array().unwrap();
    let my = mass_matches
        .iter()
        .find(|m| m["title"] == "MyCustomMod")
        .expect("user mod found");
    assert_eq!(my["source"], "user");
    assert_eq!(my["custom"], true);
    let bromo = mass_matches
        .iter()
        .find(|m| m["title"] == "Bromo")
        .expect("official Bromo found");
    assert_eq!(bromo["source"], "unimod");
    assert_eq!(bromo["custom"], false);

    // parse-proforma with a colliding name must refuse (D17 floor).
    let (ok_add, _) = run(&[
        "unimod-add",
        "--user-db",
        user_db.to_str().unwrap(),
        "--title",
        "Acetyl",
        "--formula",
        "C2 H3 O2",
        "--mono-mass",
        "43.018",
    ]);
    assert!(ok_add, "second unimod-add exit 0");
    let out = std::process::Command::new(bin_path())
        .args([
            "parse-proforma",
            "--db",
            db_path.to_str().unwrap(),
            "--user-db",
            user_db.to_str().unwrap(),
            "--proforma",
            "AC[Acetyl]",
        ])
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(1), "ambiguity exits 1");
    let err: serde_json::Value = serde_json::from_slice(&out.stderr).unwrap();
    let msg = err["error"]["message"].as_str().unwrap();
    assert!(
        msg.contains("PROFORMA_MODIFICATION_AMBIGUOUS"),
        "message: {msg}"
    );
    assert!(msg.contains("USER:2"), "mentions the user id: {msg}");

    // An explicit USER:<n> resolves without ambiguity.
    let (ok, json) = run(&[
        "parse-proforma",
        "--db",
        db_path.to_str().unwrap(),
        "--user-db",
        user_db.to_str().unwrap(),
        "--proforma",
        "AC[USER:1]",
    ]);
    assert!(ok, "USER:1 parse exit 0");
    let mods = json["result"]["modifications"].as_array().unwrap();
    assert_eq!(mods[0]["unimodId"], "USER:1");
    assert_eq!(mods[0]["source"], "user");
    assert_eq!(mods[0]["custom"], true);

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn cli_unimod_get_returns_a_record_with_specificities() {
    let db_path = unimod_db_path();
    let (ok, json) = run(&[
        "unimod-get",
        "--db",
        db_path.to_str().unwrap(),
        "--title",
        "Acetyl",
    ]);
    assert!(ok, "exit code 0");
    let r = &json["result"];
    assert_eq!(r["source"]["database"], "unimod");
    assert_eq!(r["record"]["title"], "Acetyl");
    assert_eq!(r["record"]["unimodId"], "UNIMOD:1");
    assert_eq!(r["record"]["deltaMono"].as_f64().unwrap(), 42.010565);
    assert!(!r["record"]["specificities"].as_array().unwrap().is_empty());
}

// ---- fragment matching (D-stage 6): CLI == core, and real matching ---------

const GAV_FIXTURE: &str = "fixtures/fragment-match-gav.mzML";

#[test]
fn cli_named_chemistry_agrees_with_ontology_ids() {
    for (name, id) in [
        ("Deamidation", 7),
        ("Acetyl", 1),
        ("Phospho", 21),
        ("Oxidation", 35),
    ] {
        let (ok, named) = run(&["parse-proforma", "--proforma", &format!("AQ[{name}]K")]);
        let (id_ok, identified) =
            run(&["parse-proforma", "--proforma", &format!("AQ[UNIMOD:{id}]K")]);
        assert!(ok && id_ok);
        assert_eq!(
            named["result"]["composition"],
            identified["result"]["composition"]
        );
        assert_eq!(
            named["result"]["monoisotopicMass"],
            identified["result"]["monoisotopicMass"]
        );
    }
}

#[test]
fn cli_batch_comparison_uses_each_scans_activation() {
    let xml = std::fs::read_to_string(GAV_FIXTURE).unwrap();
    let start = xml.find("<spectrum id=").unwrap();
    let end = xml.find("</spectrum>").unwrap() + "</spectrum>".len();
    let second = xml[start..end]
        .replace("scan=2", "scan=3")
        .replace("index=\"0\"", "index=\"1\"")
        .replace("MS:1000133", "MS:1000598")
        .replace(
            "collision-induced dissociation",
            "electron transfer dissociation",
        );
    let xml = xml
        .replace("<spectrumList count=\"1\"", "<spectrumList count=\"2\"")
        .replace("</spectrumList>", &format!("{second}</spectrumList>"));
    let input = std::env::temp_dir().join(format!("mixed-activation-{}.mzML", std::process::id()));
    let evidence = input.with_extension("evidence.json");
    std::fs::write(&input, xml).unwrap();
    let (ok, output) = run(&[
        "fragment",
        input.to_str().unwrap(),
        "--mode",
        "compare",
        "--scan",
        "2,3",
        "--candidate-proformas",
        r#"[{"id":"a","proforma":"GAV"},{"id":"b","proforma":"G[Acetyl]AV"}]"#,
        "--evidence-out",
        evidence.to_str().unwrap(),
    ]);
    assert!(ok, "{output}");
    let comparisons = &output["result"]["comparisons"];
    assert!(output["params"]["resolvedActivation"].is_null());
    assert!(output["params"]["ionTypes"].is_null());
    assert_eq!(comparisons[0]["ionTypes"], serde_json::json!(["b", "y"]));
    assert_eq!(comparisons[1]["ionTypes"], serde_json::json!(["c", "z"]));
    for comparison in comparisons.as_array().unwrap() {
        assert!(comparison["measurement"].get("adjudication").is_none());
    }
    std::fs::remove_file(input).unwrap();
    std::fs::remove_file(evidence).unwrap();
}

#[test]
fn cli_fragment_match_result_matches_core() {
    let fixture = GAV_FIXTURE;
    assert!(
        std::path::Path::new(fixture).exists(),
        "{fixture} must exist"
    );
    let bytes = std::fs::read(fixture).unwrap();
    let query = FragmentMatchQuery {
        selector: SpectrumSelector::Id("scan=2".to_owned()),
        sequence: "GAV".to_owned(),
        modifications: vec![MassModification {
            position: "1".to_owned(),
            delta_mass: 14.015_650_063_796,
            formula: None,
        }],
        ion_types: vec![IonType::B, IonType::Y],
        max_charge: 1,
        tolerance_ppm: 20.0,
        min_intensity: None,
        max_matches: None,
        neutral_losses: Vec::new(),
        engine_neutral_losses: false,
        min_correlation: 0.7,
    };
    let core = query_fragment_match_mzml(&bytes, &query).unwrap();
    let evidence =
        std::env::temp_dir().join(format!("tdease-fragment-match-{}.json", std::process::id()));
    let evidence_text = evidence.to_string_lossy().into_owned();

    let (ok, json) = run(&[
        "fragment",
        fixture,
        "--mode",
        "match",
        "--scan",
        "2",
        "--proforma",
        "G[Formula:CH2]AV",
        "--ion-types",
        "B,Y",
        "--max-charge",
        "1",
        "--tolerance-ppm",
        "20",
        "--evidence-out",
        &evidence_text,
    ]);
    assert!(ok, "exit code 0");
    assert!(json["result"]["measurement"]["matches"].is_null());
    assert_eq!(
        json["result"]["measurement"]["matchTableCount"],
        core.matches.len()
    );
    let sidecar: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&evidence).unwrap()).unwrap();
    assert_result_matches(
        &sidecar["result"]["measurement"],
        &serde_json::to_value(core).unwrap(),
    );
    let _ = std::fs::remove_file(evidence);
}

#[test]
fn cli_fragment_pool_conserves_intensity_and_rejects_duplicate_scans() {
    let xml = std::fs::read_to_string(GAV_FIXTURE).unwrap();
    let start = xml.find("<spectrum id=").unwrap();
    let end = xml.find("</spectrum>").unwrap() + "</spectrum>".len();
    let second = xml[start..end]
        .replace("scan=2", "scan=3")
        .replace("index=\"0\"", "index=\"1\"");
    let xml = xml
        .replace("<spectrumList count=\"1\"", "<spectrumList count=\"2\"")
        .replace("</spectrumList>", &format!("{second}</spectrumList>"));
    let input = std::env::temp_dir().join(format!("fragment-pool-{}.mzML", std::process::id()));
    let evidence = input.with_extension("json");
    std::fs::write(&input, xml).unwrap();
    let candidates = r#"[{"id":"a","proforma":"GAV"},{"id":"b","proforma":"G[Acetyl]AV"}]"#;
    let mut totals = Vec::new();
    for scans in ["2", "2,3"] {
        let (ok, output) = run(&[
            "fragment",
            input.to_str().unwrap(),
            "--mode",
            "compare",
            "--scan",
            scans,
            "--sum-scans",
            "--candidate-proformas",
            candidates,
            "--evidence-out",
            evidence.to_str().unwrap(),
        ]);
        assert!(ok, "{output}");
        let payload: serde_json::Value =
            serde_json::from_slice(&std::fs::read(&evidence).unwrap()).unwrap();
        let comparisons = payload["result"]["comparisons"].as_array().unwrap();
        let peaks = comparisons[0]["pooledPeaks"].as_array().unwrap();
        totals.push(
            peaks
                .iter()
                .map(|p| p["intensity"].as_f64().unwrap())
                .sum::<f64>(),
        );
        assert_eq!(comparisons.len(), 1);
    }
    assert!((totals[1] / totals[0] - 2.0).abs() < 1e-6);
    let result = Command::new(bin_path())
        .args([
            "fragment",
            input.to_str().unwrap(),
            "--mode",
            "compare",
            "--scan",
            "2,2",
            "--sum-scans",
            "--candidate-proformas",
            candidates,
            "--evidence-out",
            evidence.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("duplicate scan"));
    std::fs::remove_file(input).unwrap();
    std::fs::remove_file(evidence).unwrap();
}

#[test]
fn cli_fragment_compare_result_matches_core() {
    let fixture = GAV_FIXTURE;
    assert!(
        std::path::Path::new(fixture).exists(),
        "{fixture} must exist"
    );
    let bytes = std::fs::read(fixture).unwrap();
    let query = CandidateComparisonQuery {
        max_matches: None,
        neutral_losses: Vec::new(),
        engine_neutral_losses: false,
        selector: SpectrumSelector::Id("scan=2".to_owned()),
        sequence: "GAV".to_owned(),
        candidate_models: vec![
            CandidateModel {
                id: "m1".to_owned(),
                sequence: None,
                modifications: vec![MassModification {
                    position: "1".to_owned(),
                    delta_mass: 14.015_650_063_796,
                    formula: Some("CH2".to_owned()),
                }],
            },
            CandidateModel {
                id: "m2".to_owned(),
                sequence: None,
                modifications: vec![MassModification {
                    position: "2".to_owned(),
                    delta_mass: 14.015_650_063_796,
                    formula: Some("CH2".to_owned()),
                }],
            },
        ],
        ion_types: vec![IonType::B, IonType::Y],
        max_charge: 1,
        tolerance_ppm: 20.0,
        min_intensity: None,
        min_correlation: 0.7,
    };
    let core = query_candidate_comparison_mzml(&bytes, &query).unwrap();

    let evidence = std::env::temp_dir().join(format!(
        "tdease-fragment-compare-{}.json",
        std::process::id()
    ));
    let evidence_text = evidence.to_string_lossy().into_owned();
    let (ok, _json) = run(&[
        "fragment",
        fixture,
        "--mode",
        "compare",
        "--scan",
        "2",
        "--candidate-proformas",
        r#"[{"id":"m1","proforma":"G[Formula:CH2]AV"},{"id":"m2","proforma":"GA[Formula:CH2]V"}]"#,
        "--ion-types",
        "B,Y",
        "--max-charge",
        "1",
        "--tolerance-ppm",
        "20",
        "--evidence-out",
        &evidence_text,
    ]);
    assert!(ok, "exit code 0");
    let sidecar: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&evidence).unwrap()).unwrap();
    assert_result_matches(
        &sidecar["result"]["comparisons"][0]["measurement"],
        &serde_json::to_value(core).unwrap(),
    );
    let _ = std::fs::remove_file(evidence);
}

/// The adjudication of a fixed model SET must not depend on candidate input
/// order. Run `fragment --mode compare` under every ordering and require the
/// full adjudication payload to be identical.
#[test]
fn cli_fragment_compare_measurements_are_permutation_invariant() {
    let fixture = GAV_FIXTURE;
    assert!(
        std::path::Path::new(fixture).exists(),
        "{fixture} must exist"
    );
    let models: Vec<(&str, &str)> = vec![
        ("m1", "G[Formula:CH2]AV"),
        ("m2", "GA[Formula:CH2]V"),
        ("m3", "GAV[Formula:CH2]"),
    ];
    let compare = |ids: &[&str], tag: &str| -> serde_json::Value {
        let candidates = serde_json::to_string(
            &ids.iter()
                .map(|id| {
                    let (_, proforma) = models.iter().find(|(model, _)| model == id).unwrap();
                    serde_json::json!({"id": id, "proforma": proforma})
                })
                .collect::<Vec<_>>(),
        )
        .unwrap();
        let evidence = std::env::temp_dir().join(format!(
            "tdease-fragment-compare-perm-{tag}-{}.json",
            std::process::id()
        ));
        let evidence_text = evidence.to_string_lossy().into_owned();
        let (ok, stdout) = run(&[
            "fragment",
            fixture,
            "--mode",
            "compare",
            "--scan",
            "2",
            "--candidate-proformas",
            &candidates,
            "--ion-types",
            "B,Y",
            "--max-charge",
            "1",
            "--tolerance-ppm",
            "20",
            "--evidence-out",
            &evidence_text,
        ]);
        assert!(ok, "exit code 0 for {ids:?}");
        let measurement = &stdout["result"]["comparisons"][0]["measurement"];
        assert!(measurement.get("adjudication").is_none());
        let sidecar: serde_json::Value =
            serde_json::from_slice(&std::fs::read(&evidence).unwrap()).unwrap();
        assert!(sidecar["comparisons"][0]["result"]
            .get("adjudication")
            .is_none());
        let metrics: std::collections::BTreeMap<_, _> = measurement["candidates"]
            .as_array()
            .unwrap()
            .iter()
            .map(|candidate| {
                (
                    candidate["candidateId"].as_str().unwrap().to_owned(),
                    candidate.clone(),
                )
            })
            .collect();
        let _ = std::fs::remove_file(evidence);
        serde_json::to_value(metrics).unwrap()
    };

    // Two-candidate: the transposition.
    let reference = compare(&["m1", "m2"], "2a");
    assert_eq!(compare(&["m2", "m1"], "2b"), reference);

    // Three-candidate: all six orderings of the same model set.
    let reference = compare(&["m1", "m2", "m3"], "3a");
    for (ids, tag) in [
        (["m1", "m3", "m2"], "3b"),
        (["m2", "m1", "m3"], "3c"),
        (["m2", "m3", "m1"], "3d"),
        (["m3", "m1", "m2"], "3e"),
        (["m3", "m2", "m1"], "3f"),
    ] {
        assert_eq!(
            compare(&ids, tag),
            reference,
            "candidate input order changed the adjudication"
        );
    }
}

#[test]
fn cli_fragment_compare_zero_evidence_abstains() {
    let fixture = GAV_FIXTURE;
    assert!(
        std::path::Path::new(fixture).exists(),
        "{fixture} must exist"
    );
    let evidence = std::env::temp_dir().join(format!(
        "tdease-fragment-compare-zero-{}.json",
        std::process::id()
    ));
    let evidence_text = evidence.to_string_lossy().into_owned();
    // Tolerance below every error in the fixture: no candidate matches
    // anything, so the spectra carry zero evidence about the models.
    let (ok, stdout) = run(&[
        "fragment",
        fixture,
        "--mode",
        "compare",
        "--scan",
        "2",
        "--candidate-proformas",
        r#"[{"id":"m1","proforma":"G[Formula:CH2]AV"},{"id":"m2","proforma":"GA[Formula:CH2]V"}]"#,
        "--ion-types",
        "B,Y",
        "--max-charge",
        "1",
        "--tolerance-ppm",
        "0.000001",
        "--evidence-out",
        &evidence_text,
    ]);
    assert!(ok, "exit code 0");
    let measurement = &stdout["result"]["comparisons"][0]["measurement"];
    assert!(measurement.get("adjudication").is_none());
    assert_eq!(measurement["discriminatingIonCount"], 0);
    // Zero matched intensity serializes as one zero shape, never -0.0.
    for candidate in measurement["candidates"].as_array().unwrap() {
        let intensity = candidate["matchedIntensity"].as_f64().unwrap();
        assert_eq!(intensity, 0.0);
        assert!(!intensity.is_sign_negative(), "-0.0 reached the output");
    }
    let _ = std::fs::remove_file(evidence);
}

// ---- fragment-theoretical (pure theory: no spectrum needed) ---------------

#[test]
fn cli_fragment_theoretical_matches_core() {
    // The CLI is a thin adapter over `theoretical_fragment_list` — same
    // fragment set (m/z, neutral mass, ion, number) for the same ProFormA.
    let sequence = "SGRGKGGKGL";
    let modifications = vec![MassModification {
        position: "nTerm".to_owned(),
        delta_mass: 42.010_565,
        formula: None,
    }];
    let core = theoretical_fragment_list(
        sequence,
        &modifications,
        &[IonType::C, IonType::Z],
        1,
        &[],
        false,
    )
    .unwrap();

    let (ok, json) = run(&[
        "fragment-theoretical",
        "--proforma",
        "[Acetyl]-SGRGKGGKGL",
        "--ion-types",
        "C,Z",
        "--max-charge",
        "1",
    ]);
    assert!(ok, "exit code 0");
    let cli_frags = json["result"]["fragments"]
        .as_array()
        .expect("fragments array");
    assert_eq!(cli_frags.len(), core.len(), "same fragment count");
    // The CLI resolves `[Acetyl]` through the built-in named table (a slightly
    // more precise delta than the 42.010565 approximation here), so compare
    // within sub-mDa tolerance rather than bit-identical.
    let tolerance = 1e-4;
    for (cli, core_f) in cli_frags.iter().zip(&core) {
        assert!(
            (cli["mz"].as_f64().unwrap() - core_f.mz).abs() < tolerance,
            "m/z {} != {} (tolerance {tolerance})",
            cli["mz"].as_f64().unwrap(),
            core_f.mz
        );
        assert_eq!(cli["number"].as_u64().unwrap(), core_f.number as u64);
        assert_eq!(
            cli["ion"].as_str().unwrap().to_ascii_uppercase(),
            format!("{:?}", core_f.ion)
        );
    }
    assert!(
        (json["result"]["monoisotopicMass"].as_f64().unwrap() - 957.535_577_409).abs() < 1e-6,
        "parsed monoisotopic mass"
    );
}

#[test]
fn cli_resolve_shift_default_score_uses_deconvoluted_fragment_intensity() {
    let delta = 14.015_650_064;
    let fragment = theoretical_fragment_list(
        "GGGKG",
        &[MassModification {
            position: "4".to_owned(),
            delta_mass: delta,
            formula: None,
        }],
        &[IonType::C],
        1,
        &[],
        false,
    )
    .unwrap()
    .into_iter()
    .find(|fragment| fragment.number == 4)
    .expect("modified c4 fragment");

    let path = std::env::temp_dir().join(format!(
        "tdeasenexus-resolve-shift-intensity-{}.msalign",
        std::process::id()
    ));
    let library = r#"[{"name":"Methyl","deltaMass":14.015650064,"residues":["K"]}]"#;
    let score = |intensity: f64| {
        let noise = (0..8)
            .map(|index| format!("{}\t100\t1\t1.0", 10_000 + index * 100))
            .collect::<Vec<_>>()
            .join("\n");
        std::fs::write(
            &path,
            format!(
                "BEGIN IONS\nSCANS=1\nLEVEL=2\nACTIVATION=ETD\n\
                 PRECURSOR_MASS=400.0\n{}\t{}\t1\t1.0\n{}\nEND IONS\n",
                fragment.neutral_mass, intensity, noise
            ),
        )
        .unwrap();
        let (ok, json) = run(&[
            "resolve-shift",
            path.to_str().unwrap(),
            "--scan",
            "1",
            "--sequence",
            "GGGKG",
            "--shift",
            "14.015650064",
            "--library",
            library,
            "--max-components",
            "1",
            "--no-free-shift",
            "--mods",
            "[]",
            "--ion-types",
            "c",
            "--details",
        ]);
        assert!(ok, "resolve-shift succeeds: {json}");
        assert_eq!(json["params"]["scorer"], "weighted");
        json["result"]["candidates"][0]["score"]
            .as_f64()
            .expect("candidate score")
    };

    let weak = score(10.0);
    let strong = score(90.0);
    let _ = std::fs::remove_file(path);
    assert!(
        strong > weak,
        "stronger matched fragment must score higher: weak={weak}, strong={strong}"
    );
}

#[test]
fn fragment_theoretical_golden_mz_and_sequence_annotation() {
    // c1 = Acetyl + S + NH3 = 42.010565 + 87.032028 + 17.026549 = 146.069142 m/z;
    // the sequence annotation places the N-terminal Acetyl on the prefix.
    let (ok, json) = run(&[
        "fragment-theoretical",
        "--proforma",
        "[Acetyl]-SGRGKGGKGL",
        "--ion-types",
        "C,Z",
    ]);
    assert!(ok, "exit code 0");
    let frags = json["result"]["fragments"].as_array().unwrap();
    let c1 = frags
        .iter()
        .find(|f| f["ion"] == "c" && f["number"] == 1)
        .unwrap();
    assert!((c1["mz"].as_f64().unwrap() - 147.076_418_2).abs() < 1e-6);
    assert_eq!(c1["sequence"].as_str().unwrap(), "[Acetyl]S");
    let c2 = frags
        .iter()
        .find(|f| f["ion"] == "c" && f["number"] == 2)
        .unwrap();
    assert!((c2["mz"].as_f64().unwrap() - 204.097_881_9).abs() < 1e-6);
    // z9 carries the C-terminal suffix without the N-term Acetyl.
    let z9 = frags
        .iter()
        .find(|f| f["ion"] == "z" && f["number"] == 9)
        .unwrap();
    assert_eq!(z9["sequence"].as_str().unwrap(), "GRGKGGKGL");
}

#[test]
fn fragment_theoretical_preserves_exact_toppic_range_modification() {
    let proforma = "M.(S)[Acetyl]ESSSKSSQPLASKQEKDGTEKRGRGRPRKQPPK.K";
    let (ok, json) = run(&[
        "fragment-theoretical",
        "--proforma",
        proforma,
        "--ion-types",
        "C,Z",
        "--max-charge",
        "1",
    ]);
    assert!(ok, "exit code 0");
    let modification = &json["result"]["modifications"][0];
    assert_eq!(modification["localizationStart"], 1);
    assert_eq!(modification["localizationEnd"], 1);
    let c1 = json["result"]["fragments"]
        .as_array()
        .unwrap()
        .iter()
        .find(|fragment| fragment["ion"] == "c" && fragment["number"] == 1)
        .unwrap();
    let mz = c1["mz"].as_f64().unwrap();
    assert!((mz - 147.076_419).abs() < 1e-6, "acetylated c1 m/z {mz}");
    assert!(json["result"]["warnings"].as_array().unwrap().is_empty());
}

#[test]
fn fragment_theoretical_keeps_toppic_range_localization_ambiguous() {
    let (_, plain) = run(&[
        "fragment-theoretical",
        "--proforma",
        "PEPTIDE",
        "--ion-types",
        "C",
        "--max-charge",
        "1",
    ]);
    let out = Command::new(bin_path())
        .args([
            "fragment-theoretical",
            "--proforma",
            "PEP(TID)[Formula:CH2]E",
            "--ion-types",
            "C",
            "--max-charge",
            "1",
        ])
        .output()
        .unwrap();
    assert!(out.status.success());
    assert!(
        out.stderr.is_empty(),
        "successful TopPIC range parsing stays quiet"
    );
    let ranged: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    let fragments = ranged["result"]["fragments"].as_array().unwrap();
    assert!(!fragments
        .iter()
        .any(|f| f["ion"] == "c" && (f["number"] == 4 || f["number"] == 5)));
    let plain_c6 = plain["result"]["fragments"]
        .as_array()
        .unwrap()
        .iter()
        .find(|f| f["ion"] == "c" && f["number"] == 6)
        .unwrap();
    let ranged_c6 = fragments
        .iter()
        .find(|f| f["ion"] == "c" && f["number"] == 6)
        .unwrap();
    let delta = ranged_c6["mz"].as_f64().unwrap() - plain_c6["mz"].as_f64().unwrap();
    assert!(
        (delta - 14.015650064).abs() < 1e-9,
        "the range modification is applied exactly once"
    );
}

#[test]
fn fragment_theoretical_mzpaf_rejects_unrepresentable_custom_loss() {
    let out = Command::new(bin_path())
        .args([
            "fragment-theoretical",
            "--proforma",
            "PEPTIDE",
            "--ion-types",
            "B,Y",
            "--neutral-losses",
            r#"[{"label":"custom","deltaMass":10.0,"position":1}]"#,
            "--out-format",
            "mzpaf",
        ])
        .output()
        .unwrap();
    assert!(!out.status.success());
    let error: serde_json::Value = serde_json::from_slice(&out.stderr).unwrap();
    assert!(error["error"]["message"]
        .as_str()
        .unwrap()
        .contains("not representable by the mzPAF writer"));
}

#[test]
fn fragment_sequence_annotates_only_the_modified_residue() {
    // Regression: a single site modification must be tagged on exactly that
    // residue, not on every residue of a fragment that merely contains it.
    let (ok, json) = run(&[
        "fragment-theoretical",
        "--proforma", "[Bromo@Y88]-SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG",
        "--ion-types", "Y",
        "--max-charge", "1",
    ]);
    assert!(ok, "exit code 0");
    let frags = json["result"]["fragments"].as_array().unwrap();
    // y15 = C-terminal 15 residues = positions 88..102 (sequence length 102),
    // so it contains Y88 exactly once.
    let y15 = frags
        .iter()
        .find(|f| f["ion"] == "y" && f["number"] == 15)
        .unwrap();
    let seq = y15["sequence"].as_str().unwrap();
    let n_br = seq.matches("[Bromo]").count();
    assert_eq!(n_br, 1, "exactly one Bromo tag on y15, got {n_br}: {seq}");
    assert!(
        seq.starts_with("Y[Bromo]"),
        "Y88 (the N-terminal residue of y15) is the tagged site: {seq}"
    );
    // y14 = positions 89..102, does NOT contain Y88 — no Bromo tag at all.
    let y14 = frags
        .iter()
        .find(|f| f["ion"] == "y" && f["number"] == 14)
        .unwrap();
    assert_eq!(
        y14["sequence"].as_str().unwrap().matches("[Bromo]").count(),
        0,
        "y14 does not contain Y88"
    );
    // c-ion ladder carries the modification only at its site too.
    let (ok2, json2) = run(&[
        "fragment-theoretical",
        "--proforma", "[Bromo@Y88]-SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG",
        "--ion-types", "C",
        "--max-charge", "1",
    ]);
    assert!(ok2, "exit code 0");
    let cfrags = json2["result"]["fragments"].as_array().unwrap();
    let c88 = cfrags
        .iter()
        .find(|f| f["ion"] == "c" && f["number"] == 88)
        .unwrap();
    let cseq = c88["sequence"].as_str().unwrap();
    assert_eq!(
        cseq.matches("[Bromo]").count(),
        1,
        "exactly one Bromo on c88: {cseq}"
    );
}

#[test]
fn fragment_theoretical_json_includes_isotope_envelope_by_default() {
    // The json dialect carries each fragment's BRAIN isotope stick. Invariants:
    //   - every fragment has an isotopeEnvelope array (never absent);
    //   - envelope M0 (neutronOffset 0) sits at the fragment's monoisotopic m/z;
    //   - the peak at apexOffset is apex-normalized to relativeIntensity 1.0;
    //   - neutron offsets are strictly increasing.
    let (ok, json) = run(&[
        "fragment-theoretical",
        "--proforma",
        "[Acetyl]-SGRGKGGKGLGKGGAKRHRK",
        "--ion-types",
        "C,Z",
        "--max-charge",
        "14",
    ]);
    assert!(ok, "exit code 0");
    let frags = json["result"]["fragments"].as_array().unwrap();
    assert!(!frags.is_empty());
    for f in frags.iter().take(50) {
        let envelope = f["isotopeEnvelope"]
            .as_array()
            .unwrap_or_else(|| panic!("fragment {} has no isotopeEnvelope", f["ion"]));
        // Small early fragments keep only M0..M+1 (M+2 truncates away); large
        // fragments keep many more (covered by the core unit test).
        assert!(
            envelope.len() >= 2,
            "envelope keeps M0..M+1: {:?}",
            envelope
        );
        let m0 = &envelope[0];
        assert_eq!(m0["neutronOffset"].as_u64().unwrap(), 0);
        let mz = f["mz"].as_f64().unwrap();
        let env_m0 = m0["mz"].as_f64().unwrap();
        assert!(
            (env_m0 - mz).abs() < 1e-6,
            "envelope M0 {env_m0:.6} vs fragment m/z {mz:.6}"
        );
        // neutron offsets strictly increasing
        for window in envelope.windows(2) {
            let a = window[0]["neutronOffset"].as_u64().unwrap();
            let b = window[1]["neutronOffset"].as_u64().unwrap();
            assert_eq!(b, a + 1, "offsets increasing by one");
        }
        // apex normalization
        let apex_offset = f["apexOffset"].as_u64().unwrap() as usize;
        let apex = &envelope[apex_offset];
        let rel = apex["relativeIntensity"].as_f64().unwrap();
        assert!(
            (rel - 1.0).abs() < 1e-9,
            "apex-normalized intensity at apexOffset must be 1.0, got {rel}"
        );
        // apexMz agrees with the envelope apex peak
        let apex_mz = f["apexMz"].as_f64().unwrap();
        let env_apex_mz = apex["mz"].as_f64().unwrap();
        assert!(
            (apex_mz - env_apex_mz).abs() < 1e-6,
            "apexMz {apex_mz:.6} vs envelope apex {env_apex_mz:.6}"
        );
    }
}

#[test]
fn fragment_theoretical_no_envelope_omits_the_isotope_stick_keeps_apex() {
    // `--no-envelope` drops isotopeEnvelope (empty array) for bulk batch
    // output, but must keep the apex anchor the matching pipeline relies on.
    let (ok, json) = run(&[
        "fragment-theoretical",
        "--proforma",
        "[Acetyl]-SGRGKGGKGLGKGGAKRHRK",
        "--ion-types",
        "C,Z",
        "--max-charge",
        "14",
        "--no-envelope",
    ]);
    assert!(ok, "exit code 0");
    let frags = json["result"]["fragments"].as_array().unwrap();
    assert!(!frags.is_empty());
    for f in frags.iter().take(20) {
        assert!(
            f["isotopeEnvelope"].as_array().unwrap().is_empty(),
            "--no-envelope must emit an empty isotopeEnvelope for {}",
            f["ion"]
        );
        assert!(f["apexMz"].is_number(), "apexMz still present");
        assert!(f["apexOffset"].is_number(), "apexOffset still present");
    }
}

#[test]
fn fragment_theoretical_mztab_output_is_standard_and_parseable() {
    // mzTab 2.0-M: SMH header + SML row, `mzTab_id` first column, and the
    // site-resolved modification trace (`Dimethyl(residue:20)`).
    let out = Command::new(bin_path())
        .args([
            "fragment-theoretical",
            "--proforma",
            "[Dimethyl@K20]-SGRGKGGKGLGKGGAKRHRK",
            "--ion-types",
            "C,Z",
            "--out-format",
            "mztab",
        ])
        .output()
        .expect("fragment-theoretical --out-format mztab runs");
    assert!(out.status.success(), "exit code 0");
    let text = String::from_utf8(out.stdout).expect("UTF-8 mzTab output");
    assert!(
        text.starts_with("MTD\tmzTab-version\t2.0.0-M\n"),
        "mzTab 2.0-M header"
    );
    assert!(
        text.contains("SMH\tmzTab_id\tSMILES\tchemical_formula"),
        "SMH header"
    );
    assert!(text.contains("\nSML\t1\t\tC"), "SML row with mzTab_id 1");
    assert!(
        text.contains("Dimethyl(residue:20)"),
        "site-resolved modification"
    );
}

#[test]
fn cli_fragment_match_finds_all_b_y_ions_on_synthetic_fixture() {
    let evidence = std::env::temp_dir().join(format!(
        "tdease-fragment-all-ions-{}.json",
        std::process::id()
    ));
    let evidence_text = evidence.to_string_lossy().into_owned();
    let (ok, json) = run(&[
        "fragment",
        GAV_FIXTURE,
        "--mode",
        "match",
        "--scan",
        "2",
        "--proforma",
        "GAV",
        "--ion-types",
        "B,Y",
        "--max-charge",
        "1",
        "--tolerance-ppm",
        "20",
        "--evidence-out",
        &evidence_text,
    ]);
    assert!(ok, "exit code 0");
    let s = &json["result"]["measurement"]["summary"];
    assert_eq!(s["matchedIons"], 4, "all four b/y ions of GAV match");
    assert_eq!(s["coverage"], 1.0);
    assert_eq!(s["perIonType"]["b"]["longestConsecutiveRun"], 2);
    assert_eq!(s["perIonType"]["y"]["longestConsecutiveRun"], 2);
    assert_eq!(
        json["result"]["measurement"]["unmatchedPeakCount"], 0,
        "every placed peak was explained"
    );
    let _ = std::fs::remove_file(evidence);
}

#[test]
fn cli_fragment_bare_scan_resolves_controller_prefixed_native_id() {
    // Thermo data uses `controllerType=0 controllerNumber=1 scan=N`; a bare
    // `--scan 2` must resolve to that full nativeID, not fail SpectrumNotFound.
    let fixture = "fixtures/fragment-match-gav-controller.mzML";
    assert!(
        std::path::Path::new(fixture).exists(),
        "{fixture} must exist"
    );
    let evidence = std::env::temp_dir().join(format!(
        "tdease-fragment-controller-{}.json",
        std::process::id()
    ));
    let evidence_text = evidence.to_string_lossy().into_owned();
    let (ok, json) = run(&[
        "fragment",
        fixture,
        "--mode",
        "match",
        "--scan",
        "2",
        "--proforma",
        "GAV",
        "--ion-types",
        "B,Y",
        "--max-charge",
        "1",
        "--tolerance-ppm",
        "20",
        "--evidence-out",
        &evidence_text,
    ]);
    assert!(ok, "exit code 0");
    assert_eq!(
        json["result"]["measurement"]["summary"]["matchedIons"], 4,
        "controller-prefixed scan resolves"
    );
    let _ = std::fs::remove_file(evidence);
}

#[test]
fn cli_fragment_compare_reports_discriminating_ions() {
    let evidence = std::env::temp_dir().join(format!(
        "tdease-fragment-discriminating-{}.json",
        std::process::id()
    ));
    let evidence_text = evidence.to_string_lossy().into_owned();
    let (ok, json) = run(&[
        "fragment",
        GAV_FIXTURE,
        "--mode",
        "compare",
        "--scan",
        "2",
        "--candidate-proformas",
        r#"[{"id":"m1","proforma":"G[Formula:CH2]AV"},{"id":"m2","proforma":"GA[Formula:CH2]V"}]"#,
        "--ion-types",
        "B,Y",
        "--max-charge",
        "1",
        "--tolerance-ppm",
        "20",
        "--evidence-out",
        &evidence_text,
    ]);
    assert!(ok, "exit code 0");
    assert_eq!(
        json["result"]["measurementRole"],
        "candidate_model_evaluation"
    );
    // The verdict travels in the bounded stdout itself, not only the artifact.
    let bounded = &json["result"]["comparisons"][0]["measurement"];
    assert_eq!(bounded["discriminatingIonCount"], 2);
    let ions = bounded["discriminatingIons"].as_array().unwrap();
    // b1 discriminates position 2, y2 discriminates position 1.
    assert!(ions.iter().any(|i| i["ion"] == "b"
        && i["number"] == 1
        && i["matchedCandidates"] == serde_json::json!(["m2"])));
    assert!(ions.iter().any(|i| i["ion"] == "y"
        && i["number"] == 2
        && i["matchedCandidates"] == serde_json::json!(["m1"])));
    let sidecar: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&evidence).unwrap()).unwrap();
    let r = &sidecar["result"]["comparisons"][0]["measurement"];
    assert_eq!(r["discriminatingIons"], bounded["discriminatingIons"]);
    let _ = std::fs::remove_file(evidence);
}

#[test]
fn cli_fragment_compare_without_candidates_is_a_clean_usage_error() {
    let out = Command::new(bin_path())
        .args([
            "fragment",
            GAV_FIXTURE,
            "--mode",
            "compare",
            "--scan",
            "2",
            "--ion-types",
            "B,Y",
            "--tolerance-ppm",
            "20",
        ])
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(1), "runtime error exits 1");
    assert!(out.stdout.is_empty(), "stdout stays empty on error");
    let err: serde_json::Value = serde_json::from_slice(&out.stderr).unwrap();
    assert_eq!(err["error"]["code"], "USAGE");
    assert!(err["error"]["message"]
        .as_str()
        .unwrap()
        .contains("--candidate-proformas"));
}

#[test]
fn cli_fragment_match_accepts_mass_delta_and_labels_the_rung() {
    // Evidence lattice: fragment m/z depends only on summed deltas, so a
    // mass-only proforma is a legal match input — the result carries the
    // massDelta rung so downstream adjudication knows no composition or
    // identity was tested.
    let evidence = std::env::temp_dir().join(format!(
        "tdease-fragment-massdelta-{}.json",
        std::process::id()
    ));
    let evidence_text = evidence.to_string_lossy().into_owned();
    let (ok, json) = run(&[
        "fragment",
        GAV_FIXTURE,
        "--mode",
        "match",
        "--scan",
        "2",
        "--proforma",
        "G[+14.01565]AV",
        "--ion-types",
        "B,Y",
        "--tolerance-ppm",
        "20",
        "--evidence-out",
        &evidence_text,
    ]);
    assert!(ok, "exit code 0");
    assert_eq!(json["result"]["hypothesisType"], "massDelta");
    assert_eq!(
        json["result"]["measurementRole"],
        "candidate_model_evaluation"
    );
    let sidecar: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&evidence).unwrap()).unwrap();
    assert_eq!(
        sidecar["result"]["hypothesisType"], "massDelta",
        "sidecar labeled too"
    );
    assert!(
        !sidecar["result"]["measurement"]["matches"]
            .as_array()
            .unwrap()
            .is_empty(),
        "same delta mass matches the same fragments as a Formula tag"
    );
    let _ = std::fs::remove_file(evidence);
}

#[test]
fn cli_fragment_theoretical_mass_delta_degrades_to_averagine() {
    // The massDelta rung keeps exact fragment m/z but substitutes the
    // averagine model for the isotope shape: neutron positions are exact,
    // intensities approximate, and the payload says so.
    let proforma = "G[+14.01565]AV".to_owned();
    let core = theoretical_fragment_list_proforma_variants(
        std::slice::from_ref(&proforma),
        &[],
        &[IonType::B, IonType::Y],
        1,
        &[],
        false,
        EnvelopeMode::Averagine,
    )
    .unwrap();
    assert!(
        core.iter()
            .all(|fragment| !fragment.isotope_envelope.is_empty()),
        "averagine keeps an envelope for matching positions"
    );
    let (ok, json) = run(&[
        "fragment-theoretical",
        "--proforma",
        &proforma,
        "--ion-types",
        "B,Y",
        "--max-charge",
        "1",
    ]);
    assert!(ok, "exit code 0");
    assert_eq!(json["result"]["hypothesisType"], "massDelta");
    let warnings = json["result"]["warnings"].as_array().unwrap();
    assert!(
        warnings
            .iter()
            .any(|warning| warning.as_str().unwrap().contains("averagine")),
        "the averagine degradation is declared, not silent: {warnings:?}"
    );
    let cli_frags = json["result"]["fragments"].as_array().unwrap();
    assert_eq!(cli_frags.len(), core.len(), "same fragment count");
    for (cli, core_fragment) in cli_frags.iter().zip(&core) {
        assert!(
            (cli["mz"].as_f64().unwrap() - core_fragment.mz).abs() < 1e-9,
            "massDelta ladder is exact"
        );
        assert_eq!(
            cli["isotopeEnvelope"].as_array().unwrap().len(),
            core_fragment.isotope_envelope.len(),
            "averagine envelope emitted per fragment"
        );
    }
}

#[test]
fn spectrum_scans_honors_the_ms_level_filter() {
    // GAV scan 2 is an MS2 spectrum: the default (0 = any) returns it, an
    // explicit --ms-level 1 sends it to missingScans, --ms-level 2 returns it.
    let (ok, json) = run(&["spectrum", GAV_FIXTURE, "--mode", "scans", "--scans", "2"]);
    assert!(ok);
    assert_eq!(
        json["result"]["spectra"].as_array().unwrap().len(),
        1,
        "any-level default"
    );

    // An explicit --ms-level 1 makes the MS2 scan "not exist" at that level:
    // the same typed NOT_FOUND error as a nonexistent scan (empty stdout,
    // error envelope on stderr).
    let out = Command::new(bin_path())
        .args([
            "spectrum",
            GAV_FIXTURE,
            "--mode",
            "scans",
            "--scans",
            "2",
            "--ms-level",
            "1",
        ])
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert!(out.stdout.is_empty());
    let err: serde_json::Value = serde_json::from_slice(&out.stderr).unwrap();
    assert_eq!(err["error"]["code"], "MASS_SPEC_SPECTRUM_NOT_FOUND");

    let (ok, json) = run(&[
        "spectrum",
        GAV_FIXTURE,
        "--mode",
        "scans",
        "--scans",
        "2",
        "--ms-level",
        "2",
    ]);
    assert!(ok);
    assert_eq!(
        json["result"]["spectra"].as_array().unwrap().len(),
        1,
        "matching level"
    );
}

#[test]
fn mgf_export_prefers_the_explicit_charge_override() {
    // PEPMASS stays the raw acquisition fact; an explicit --charge writes the
    // CHARGE header (acquisition headers are known to carry wrong charges).
    let out = Command::new(bin_path())
        .args([
            "spectrum",
            GAV_FIXTURE,
            "--mode",
            "scans",
            "--scans",
            "2",
            "--charge",
            "17",
            "--out-format",
            "mgf",
        ])
        .output()
        .unwrap();
    assert!(out.status.success());
    let text = String::from_utf8(out.stdout).unwrap();
    assert!(text.contains("CHARGE=17+"), "override wins: {text}");
}

#[test]
fn cli_validate_precursor_rejects_mass_only_proforma() {
    // The precursor validation is an isotope-envelope test: a mass-only delta
    // cannot define the shape, so it is rejected instead of silently matching
    // a base-composition envelope.
    let out = Command::new(bin_path())
        .args([
            "validate-precursor",
            CHLAMY_FIXTURE,
            "--scan",
            "4113",
            "--charge",
            "11",
            "--proforma",
            "SGRGK[+42.010565]GGAKRHRK",
        ])
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(1));
    assert!(out.stdout.is_empty());
    let err: serde_json::Value = serde_json::from_slice(&out.stderr).unwrap();
    assert_eq!(err["error"]["code"], "MASS_SPEC_INVALID_QUERY");
    assert!(err["error"]["message"]
        .as_str()
        .unwrap()
        .contains("isotope-envelope test"));
}

#[test]
fn cli_complete_fragment_compare_is_bounded_and_file_backed() {
    let evidence = std::env::temp_dir().join(format!(
        "tdease-fragment-complete-{}-{}.json",
        std::process::id(),
        std::thread::current().name().unwrap_or("test")
    ));
    let evidence_text = evidence.to_string_lossy().into_owned();
    let (ok, json) = run(&[
        "fragment",
        GAV_FIXTURE,
        "--mode",
        "compare",
        "--scan",
        "2",
        "--candidate-proformas",
        r#"[{"id":"m1","proforma":"G[Formula:CH2]AV"},{"id":"m2","proforma":"GA[Formula:CH2]V"}]"#,
        "--ion-types",
        "B,Y",
        "--max-charge",
        "1",
        "--tolerance-ppm",
        "20",
        "--evidence-out",
        &evidence_text,
    ]);
    assert!(ok, "complete candidate comparison succeeds");
    let result = &json["result"];
    assert_eq!(result["measurementRole"], "candidate_model_evaluation");
    assert_eq!(result["notIdentity"], true);
    assert_eq!(result["scans"].as_array().unwrap().len(), 1);
    let comparisons = result["comparisons"].as_array().unwrap();
    assert_eq!(comparisons.len(), 1);
    for entry in comparisons {
        let measurement = &entry["measurement"];
        assert!(measurement["candidates"][0].get("matches").is_none());
        assert!(measurement["candidates"][0]["matchTableCount"].is_number());
    }
    assert!(result["evidenceRef"]["sha256"].is_string());

    let sidecar: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&evidence).unwrap()).unwrap();
    // The artifact is the same CLI envelope as stdout (one schema); full
    // tables live at the same paths the bounded view compacts.
    assert_eq!(sidecar["schemaVersion"], 1);
    assert_eq!(sidecar["command"], "fragment");
    assert_eq!(sidecar["result"]["notIdentity"], true);
    let entries = sidecar["result"]["comparisons"].as_array().unwrap();
    assert_eq!(entries.len(), 1);
    assert!(entries[0]["measurement"]["candidates"][0]["matches"].is_array());
    let _ = std::fs::remove_file(evidence);
}

#[test]
fn cli_fragment_match_rejects_duplicate_scans() {
    let out = Command::new(bin_path())
        .args([
            "fragment",
            GAV_FIXTURE,
            "--mode",
            "match",
            "--scan",
            "2,2",
            "--proforma",
            "GAV",
            "--ion-types",
            "B,Y",
            "--tolerance-ppm",
            "20",
            "--evidence-out",
            "unused.json",
        ])
        .output()
        .unwrap();
    assert!(!out.status.success());
    let json: serde_json::Value =
        serde_json::from_slice(&out.stderr).expect("typed error envelope on stderr");
    assert_eq!(json["error"]["code"], "MASS_SPEC_INVALID_QUERY");
    assert!(json["error"]["message"]
        .as_str()
        .unwrap()
        .contains("duplicate scan"));
}

#[test]
fn cli_spectrum_scans_exports_mgf() {
    // The isoforma skill consumes `tdeasenexus-toolbox spectrum --out-format mgf`
    // directly (raw-m/z peak exchange, community standard). The MGF must carry
    // one BEGIN IONS block per scan with PEPMASS/CHARGE/RT and the peak pairs.
    let out = Command::new(bin_path())
        .args([
            "spectrum",
            GAV_FIXTURE,
            "--mode",
            "scans",
            "--scans",
            "2",
            "--out-format",
            "mgf",
        ])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "exit code 0: {:?}",
        String::from_utf8_lossy(&out.stderr)
    );
    let text = String::from_utf8_lossy(&out.stdout);
    assert!(
        text.contains("BEGIN IONS"),
        "has at least one spectrum block"
    );
    assert!(text.contains("END IONS"), "blocks are terminated");
    assert!(text.contains("TITLE=scan=2"), "scan number in TITLE");
    assert!(text.contains("PEPMASS="), "precursor m/z present");
    assert!(text.contains("CHARGE="), "charge present");
    assert!(text.contains("RTINSECONDS="), "retention time present");
    // the first peak line must be "mz intensity" (two whitespace-separated floats)
    let first_peak = text
        .lines()
        .find(|line| {
            !line.starts_with("BEGIN")
                && !line.starts_with("END")
                && !line.contains("=")
                && !line.contains("MGF from")
                && line.split_whitespace().count() == 2
        })
        .expect("a peak line of the form 'mz intensity'");
    let parts: Vec<&str> = first_peak.split_whitespace().collect();
    let _: f64 = parts[0].parse().expect("m/z is numeric");
    let _: f64 = parts[1].parse().expect("intensity is numeric");
}

#[test]
fn cli_select_matches_core_on_real_fixture() {
    let fixture = "fixtures/fragment-match-gav.mzML";
    assert!(
        std::path::Path::new(fixture).exists(),
        "{fixture} must exist"
    );
    let statement = "SELECT spectra FROM mzml WHERE product_mz WITHIN 189.1234 PPM 10 LIMIT 5;";
    let bytes = std::fs::read(fixture).unwrap();
    let core = query_select_mzml(&bytes, statement).unwrap();

    let (ok, json) = run(&["select", "--statement", statement, fixture]);
    assert!(ok, "exit code 0");
    assert_result_matches(&json["result"], &serde_json::to_value(core).unwrap());
    assert_eq!(json["result"]["matchedSpectrumCount"], 1);
    assert_eq!(
        json["result"]["spectra"][0]["matchedProducts"][0]["intensity"],
        9000.0
    );
}

#[test]
fn cli_select_tic_uses_header_without_peak_decode() {
    // This fixture declares MS:1000285 in the spectrum header: a metadata-only
    // statement must answer from the header alone (peakCount null, no
    // matchedProducts) instead of decoding the binary arrays.
    let statement = "SELECT spectra FROM mzml WHERE ms_level = 1 LIMIT 10;";
    let (ok, json) = run(&[
        "select",
        "--statement",
        statement,
        "fixtures/select-tic-header.mzML",
    ]);
    assert!(ok, "exit code 0");
    let spectra = &json["result"]["spectra"];
    assert_eq!(spectra.as_array().unwrap().len(), 1);
    assert_eq!(spectra[0]["totalIonCurrent"], 12345.6);
    assert!(spectra[0]["peakCount"].is_null(), "peaks were not decoded");
    assert_eq!(spectra[0]["matchedProducts"].as_array().unwrap().len(), 0);
}

#[test]
fn cli_select_reports_native_fragment_series() {
    let (ok, json) = run(&[
        "select",
        "fixtures/histone48-scan1164.mzML",
        "--statement",
        "SELECT spectra FROM mzml WHERE ms_level = 2 LIMIT 1",
    ]);
    assert!(ok);
    let precursor = &json["result"]["spectra"][0]["precursor"];
    assert_eq!(precursor["activation"], "electron transfer dissociation");
    assert_eq!(precursor["fragmentIonTypes"], serde_json::json!(["c", "z"]));
}

#[test]
fn cli_select_precursor_neutral_mass_windows_before_peak_decode() {
    // The family-window predicate: neutral mass = selected m/z × z − z·proton
    // from the scan's own charge header, one interval across charges. It is
    // evaluated at the metadata stage (the original IsoForma
    // pull_scan_numbers filter) and excludes scans without a charge header.
    // fixture: scan 2 is MS2, selected m/z 229.14 at z=1 → neutral ≈ 228.13.
    let fixture = "fixtures/fragment-match-gav.mzML";
    let inside = "SELECT spectra FROM mzml WHERE ms_level = 2 AND precursor_neutral_mass BETWEEN 227.0 AND 229.5;";
    let bytes = std::fs::read(fixture).unwrap();
    let core = query_select_mzml(&bytes, inside).unwrap();
    let (ok, json) = run(&["select", "--statement", inside, fixture]);
    assert!(ok, "exit code 0");
    assert_result_matches(&json["result"], &serde_json::to_value(core).unwrap());
    assert_eq!(json["result"]["matchedSpectrumCount"], 1);
    assert_eq!(json["result"]["spectra"][0]["scan"], "2");
    // No metadata-only assertion here: this fixture omits the TIC header, so
    // the matched spectrum decodes peaks for the TIC fallback (existing
    // contract, covered by cli_select_tic_uses_header_without_peak_decode).
    // The window's cost contract is that NON-matching spectra never decode.

    let outside = "SELECT spectra FROM mzml WHERE ms_level = 2 AND precursor_neutral_mass BETWEEN 300.0 AND 400.0;";
    let (ok, json) = run(&["select", "--statement", outside, fixture]);
    assert!(ok, "exit code 0");
    assert_eq!(json["result"]["matchedSpectrumCount"], 0);
}

// ---- neutral losses: site gating is the scientific-correctness contract ----
//
// A modification-defined loss (H3PO4 from a phospho) exists only on fragments
// whose residue span contains the modified site. Emitting it elsewhere
// manufactures chemically impossible theoretical peaks that match noise inside
// the tolerance window — the false-positive channel evidence-bounds.md warns
// about. `position` is therefore required and gated.

#[test]
fn cli_neutral_losses_are_gated_to_the_stated_site() {
    // L = 8, phospho-style +79.9663 on residue 5: c_k spans 1..=k (needs
    // k >= 5), z_k spans 9-k..=8 (needs k >= 4). c4/z3 must have no -H3PO4
    // variant; c5/z4 must.
    let proforma = "SGRGK[Phospho]GAK";
    let loss = r#"[{"label":"H3PO4","deltaMass":97.9769,"position":5}]"#;
    let (ok, json) = run(&[
        "fragment-theoretical",
        "--proforma",
        proforma,
        "--ion-types",
        "C,Z",
        "--max-charge",
        "2",
        "--neutral-losses",
        loss,
    ]);
    assert!(ok, "exit code 0");
    let rows = json["result"]["fragments"].as_array().unwrap();
    let has_loss = |ion: &str, number: u64, loss: &str| {
        rows.iter().any(|r| {
            r["ion"] == ion && r["number"].as_u64() == Some(number) && r["neutralLoss"] == loss
        })
    };
    assert!(!has_loss("c", 4, "H3PO4"), "c4 does not contain residue 5");
    assert!(!has_loss("z", 3, "H3PO4"), "z3 does not contain residue 5");
    assert!(has_loss("c", 5, "H3PO4"), "c5 brackets residue 5");
    assert!(has_loss("z", 4, "H3PO4"), "z4 brackets residue 5");
    // No loss variant may sit on a fragment that does not contain the site.
    for r in rows.iter().filter(|r| r["neutralLoss"] == "H3PO4") {
        let number = r["number"].as_u64().unwrap();
        match r["ion"].as_str().unwrap() {
            "c" => assert!(number >= 5, "c{number} cannot lose H3PO4"),
            "z" => assert!(number >= 4, "z{number} cannot lose H3PO4"),
            other => panic!("unexpected ion {other}"),
        }
    }
    // The position accepts the string form too, with identical output.
    let loss_str = r#"[{"label":"H3PO4","deltaMass":97.9769,"position":"5"}]"#;
    let (ok2, json2) = run(&[
        "fragment-theoretical",
        "--proforma",
        proforma,
        "--ion-types",
        "C,Z",
        "--max-charge",
        "2",
        "--neutral-losses",
        loss_str,
    ]);
    assert!(ok2);
    assert_eq!(
        json["result"]["fragmentCount"], json2["result"]["fragmentCount"],
        "integer and string position forms agree"
    );
}

#[test]
fn cli_neutral_losses_require_a_position() {
    let run_err = |losses: &str| {
        let out = Command::new(bin_path())
            .args([
                "fragment-theoretical",
                "--proforma",
                "SGRGKGAK",
                "--ion-types",
                "C,Z",
                "--neutral-losses",
                losses,
            ])
            .output()
            .unwrap();
        let error: serde_json::Value = serde_json::from_slice(&out.stderr).unwrap();
        (!out.status.success(), error)
    };
    let (failed, error) = run_err(r#"[{"label":"H3PO4","deltaMass":97.9769}]"#);
    assert!(failed, "position-less loss must fail");
    assert!(
        error["error"]["message"]
            .as_str()
            .unwrap_or_default()
            .contains("position"),
        "error names the missing position"
    );
    let (failed, _) = run_err(r#"[{"label":"H3PO4","deltaMass":97.9769,"position":0}]"#);
    assert!(failed, "position 0 must fail (1-based)");
    let (failed, _) = run_err(r#"[{"label":"H3PO4","deltaMass":97.9769,"position":99}]"#);
    assert!(failed, "position beyond the sequence must fail");
}

#[test]
fn cli_engine_neutral_losses_are_opt_in() {
    // Off (the default): zero loss rows — recorded no-loss runs stay
    // byte-for-byte reproducible. On (`engine`): the fragmentation engine's
    // own chemically-gated variants appear, labelled in Hill notation.
    let proforma = "SGRGKGAK";
    let (ok, json) = run(&[
        "fragment-theoretical",
        "--proforma",
        proforma,
        "--ion-types",
        "C,Z",
        "--max-charge",
        "1",
    ]);
    assert!(ok);
    let rows = json["result"]["fragments"].as_array().unwrap();
    assert!(
        rows.iter().all(|r| r["neutralLoss"].is_null()),
        "default output carries no loss variants"
    );

    let (ok, json) = run(&[
        "fragment-theoretical",
        "--proforma",
        proforma,
        "--ion-types",
        "C,Z",
        "--max-charge",
        "1",
        "--neutral-losses",
        "engine",
    ]);
    assert!(ok);
    let rows = json["result"]["fragments"].as_array().unwrap();
    let labelled: Vec<&serde_json::Value> = rows
        .iter()
        .filter(|r| !r["neutralLoss"].is_null())
        .collect();
    assert!(
        !labelled.is_empty(),
        "engine mode emits the engine's loss variants"
    );
    assert!(
        labelled
            .iter()
            .all(|r| r["neutralLoss"].as_str().unwrap().contains("H2O")),
        "engine labels are Hill notation (e.g. -H2O1)"
    );
    // Symmetric: both c and z series carry the engine loss.
    for ion in ["c", "z"] {
        assert!(
            labelled.iter().any(|r| r["ion"] == ion),
            "{ion} series carries engine loss variants"
        );
    }
}
