//! Real CLI paths: selection -> pooling -> saved observation -> measurement.
use serde_json::{json, Value};
use std::{
    path::{Path, PathBuf},
    process::Command,
};

struct Workspace(PathBuf);
impl Workspace {
    fn new() -> Self {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path =
            std::env::temp_dir().join(format!("tdp-observation-{}-{nonce}", std::process::id()));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn file(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
    fn mzml(&self, name: &str, level: u8, count: usize) -> PathBuf {
        let xml = std::fs::read_to_string(
            Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/fragment-match-gav.mzML"),
        )
        .unwrap();
        let start = xml.find("<spectrum id=").unwrap();
        let end = xml.find("</spectrum>").unwrap() + "</spectrum>".len();
        let spectra: String = (0..count)
            .map(|i| {
                xml[start..end]
                    .replace("scan=2", &format!("scan={}", 2 + i * 2))
                    .replace("index=\"0\"", &format!("index=\"{i}\""))
                    .replace(
                        "name=\"ms level\" value=\"2\"",
                        &format!("name=\"ms level\" value=\"{level}\""),
                    )
                    .replace(
                        "name=\"scan start time\" value=\"0\"",
                        &format!("name=\"scan start time\" value=\"{i}\""),
                    )
            })
            .collect();
        let content = format!("{}{}{}", &xml[..start], spectra, &xml[end..]).replace(
            "<spectrumList count=\"1\"",
            &format!("<spectrumList count=\"{count}\""),
        );
        let path = self.file(name);
        std::fs::write(&path, content).unwrap();
        path
    }
    fn save(&self, name: &str, value: &Value) -> PathBuf {
        let path = self.file(name);
        std::fs::write(&path, serde_json::to_vec(value).unwrap()).unwrap();
        path
    }
}
impl Drop for Workspace {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn invoke(args: &[&str]) -> (bool, Value) {
    let output = Command::new(env!("CARGO_BIN_EXE_tdeasenexus-toolbox"))
        .args(args)
        .output()
        .unwrap();
    let data = if output.status.success() {
        &output.stdout
    } else {
        &output.stderr
    };
    (
        output.status.success(),
        serde_json::from_slice(data)
            .unwrap_or_else(|_| panic!("{}", String::from_utf8_lossy(data))),
    )
}
fn ok(args: &[&str]) -> Value {
    let (success, value) = invoke(args);
    assert!(success, "{value}");
    value
}
fn read(path: &Path) -> Value {
    serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap()
}

#[test]
fn selection_and_pooling_are_consistent_for_ms1_and_ms2() {
    let workspace = Workspace::new();
    for level in [1, 2] {
        let input = workspace.mzml(&format!("ms{level}.mzML"), level, 2);
        let level = level.to_string();
        let mut reference = None;
        for selection in [
            vec!["--scans", "2,4"],
            vec!["--scan-min", "2", "--scan-max", "4"],
            vec!["--rt-min", "0", "--rt-max", "1"],
        ] {
            let mut args = vec![
                "spectrum",
                input.to_str().unwrap(),
                "--mode",
                "sum",
                "--ms-level",
                &level,
            ];
            args.extend(selection);
            let output = ok(&args);
            assert_eq!(
                output["result"]["spectraMetadata"]
                    .as_array()
                    .unwrap()
                    .len(),
                2
            );
            assert_eq!(output["result"]["truncated"], false);
            if let Some(expected) = &reference {
                assert_eq!(&output["result"]["peaks"], expected);
            } else {
                reference = Some(output["result"]["peaks"].clone());
            }
        }
        let window = ok(&[
            "spectrum",
            input.to_str().unwrap(),
            "--mode",
            "window",
            "--ms-level",
            &level,
            "--rt-min",
            "0",
            "--rt-max",
            "1",
        ]);
        assert_eq!(window["result"]["peaks"], reference.unwrap());
        // Predicates intersect rather than ignoring the explicit scan list.
        let intersection = ok(&[
            "spectrum",
            input.to_str().unwrap(),
            "--mode",
            "sum",
            "--scans",
            "2,4",
            "--rt-min",
            "1",
            "--ms-level",
            &level,
        ]);
        assert_eq!(
            intersection["result"]["spectraMetadata"]
                .as_array()
                .unwrap()
                .len(),
            1
        );
        assert_eq!(intersection["result"]["spectraMetadata"][0]["id"], "scan=4");
    }
    let input = workspace.mzml("long.mzML", 1, 105);
    let output = ok(&[
        "spectrum",
        input.to_str().unwrap(),
        "--mode",
        "window",
        "--rt-min",
        "0",
        "--rt-max",
        "104",
    ]);
    assert_eq!(output["result"]["scansSummed"], 105);
    assert_eq!(output["result"]["truncated"], false);
    let limited = ok(&[
        "spectrum",
        input.to_str().unwrap(),
        "--mode",
        "window",
        "--rt-min",
        "0",
        "--rt-max",
        "104",
        "--max-scans",
        "5",
    ]);
    assert_eq!(limited["result"]["scansSummed"], 5);
    assert_eq!(limited["result"]["truncated"], true);
}

#[test]
fn one_model_batches_and_pooling_use_the_single_measurement_engine() {
    let workspace = Workspace::new();
    let input = workspace.mzml("batch.mzML", 2, 2);
    let evidence = workspace.file("evidence.json");
    let batch = ok(&[
        "fragment",
        input.to_str().unwrap(),
        "--scan",
        "2,4",
        "--proforma",
        "GAV",
        "--evidence-out",
        evidence.to_str().unwrap(),
    ]);
    assert_eq!(batch["result"]["comparisons"].as_array().unwrap().len(), 2);
    assert_eq!(batch["params"]["preparedTheorySets"], 1);
    let full = read(&evidence);
    for (i, scan) in ["2", "4"].iter().enumerate() {
        ok(&[
            "fragment",
            input.to_str().unwrap(),
            "--scan",
            scan,
            "--proforma",
            "GAV",
            "--evidence-out",
            evidence.to_str().unwrap(),
        ]);
        let single = read(&evidence);
        let candidate = &full["result"]["comparisons"][i]["measurement"]["candidates"][0];
        assert_eq!(
            candidate["matches"],
            single["result"]["measurement"]["matches"]
        );
        assert_eq!(
            candidate["fragmentStatuses"],
            single["result"]["measurement"]["fragmentStatuses"]
        );
    }
    ok(&[
        "fragment",
        input.to_str().unwrap(),
        "--scan-min",
        "2",
        "--scan-max",
        "4",
        "--sum-scans",
        "--proforma",
        "GAV",
        "--evidence-out",
        evidence.to_str().unwrap(),
    ]);
    let pooled = read(&evidence);
    assert_eq!(
        pooled["result"]["observation"]["pooledScans"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    assert!(!pooled["result"]["measurement"]["matches"]
        .as_array()
        .unwrap()
        .is_empty());
}

#[test]
fn saved_observations_and_raw_pools_agree_for_independent_backbones() {
    let workspace = Workspace::new();
    let input = workspace.mzml("source.mzML", 2, 2);
    let models = r#"[{"id":"full","proforma":"GAV"},{"id":"short","proforma":"AV"},{"id":"terminal","proforma":"GAV-[Amidated]"}]"#;
    let evidence = workspace.file("evidence.json");
    ok(&[
        "fragment",
        input.to_str().unwrap(),
        "--scan",
        "2,4",
        "--sum-scans",
        "--candidate-proformas",
        models,
        "--evidence-out",
        evidence.to_str().unwrap(),
    ]);
    let integrated = read(&evidence);
    let output = ok(&[
        "spectrum",
        input.to_str().unwrap(),
        "--mode",
        "sum",
        "--scans",
        "2,4",
        "--bin-ppm",
        "10",
    ]);
    let saved = workspace.save("observation.json", &output);
    let measured = ok(&[
        "fragment",
        saved.to_str().unwrap(),
        "--candidate-proformas",
        models,
        "--evidence-out",
        evidence.to_str().unwrap(),
    ]);
    let reused = read(&evidence);
    assert_eq!(
        reused["result"]["comparisons"][0]["measurement"],
        integrated["result"]["comparisons"][0]["measurement"]
    );
    let comparison = &measured["result"]["comparisons"][0]["measurement"];
    assert_eq!(comparison["sharedBackbone"], false);
    assert_eq!(comparison["candidates"][1]["sequence"], "AV");
    assert_eq!(comparison["discriminatingIonCount"], 0);
    assert!(comparison["candidates"][0]["theoreticalNeutralMass"].is_number());
    assert!(!reused["params"]["observationSource"]
        .as_array()
        .unwrap()
        .is_empty());
}

#[test]
fn comparison_forwards_neutral_loss_settings() {
    let workspace = Workspace::new();
    let input = workspace.mzml("loss.mzML", 2, 1);
    let evidence = workspace.file("evidence.json");
    ok(&[
        "fragment",
        input.to_str().unwrap(),
        "--scan",
        "2",
        "--proforma",
        "GAV",
        "--neutral-losses",
        "engine",
        "--evidence-out",
        evidence.to_str().unwrap(),
    ]);
    let single = read(&evidence);
    ok(&[
        "fragment",
        input.to_str().unwrap(),
        "--scan",
        "2",
        "--candidate-proformas",
        r#"[{"id":"a","proforma":"GAV"}]"#,
        "--neutral-losses",
        "engine",
        "--evidence-out",
        evidence.to_str().unwrap(),
    ]);
    let comparison = read(&evidence);
    assert_eq!(
        comparison["result"]["comparisons"][0]["measurement"]["candidates"][0]["matches"],
        single["result"]["measurement"]["matches"]
    );
    assert_eq!(
        comparison["result"]["comparisons"][0]["measurement"]["candidates"][0]["fragmentStatuses"],
        single["result"]["measurement"]["fragmentStatuses"]
    );
}

#[test]
fn invalid_or_incomplete_observations_do_not_become_measurements() {
    let workspace = Workspace::new();
    let input = workspace.mzml("source.mzML", 2, 2);
    let evidence = workspace.file("evidence.json");
    let output = ok(&[
        "spectrum",
        input.to_str().unwrap(),
        "--mode",
        "sum",
        "--scans",
        "2,4",
    ]);
    for (label, modified) in [
        ("truncated", {
            let mut d = output.clone();
            d["result"]["truncated"] = json!(true);
            d
        }),
        ("version", {
            let mut d = output.clone();
            d["schemaVersion"] = json!(2);
            d
        }),
        ("ms1", {
            let mut d = output.clone();
            d["result"]["spectraMetadata"][0]["msLevel"] = json!(1);
            d
        }),
        ("profile", {
            let mut d = output.clone();
            d["result"]["spectraMetadata"][0]["centroided"] = json!(false);
            d
        }),
        ("duplicate", {
            let mut d = output.clone();
            d["result"]["spectraMetadata"][1] = d["result"]["spectraMetadata"][0].clone();
            d
        }),
    ] {
        let path = workspace.save(&format!("{label}.json"), &modified);
        let (success, error) = invoke(&[
            "fragment",
            path.to_str().unwrap(),
            "--proforma",
            "GAV",
            "--evidence-out",
            evidence.to_str().unwrap(),
        ]);
        assert!(!success, "{label}: {error}");
        assert!(!evidence.exists());
    }
    let (success, _) = invoke(&[
        "spectrum",
        input.to_str().unwrap(),
        "--mode",
        "sum",
        "--scans",
        "2,999",
    ]);
    assert!(!success);
    let (success, _) = invoke(&[
        "spectrum",
        input.to_str().unwrap(),
        "--mode",
        "sum",
        "--scans",
        "2,scan=2",
    ]);
    assert!(!success);
    let mut mixed = output.clone();
    mixed["result"]["spectraMetadata"][1]["precursor"]["activation"] = json!("ETD");
    let mixed = workspace.save("mixed.json", &mixed);
    let (success, _) = invoke(&[
        "fragment",
        mixed.to_str().unwrap(),
        "--activation",
        "CID",
        "--proforma",
        "GAV",
        "--evidence-out",
        evidence.to_str().unwrap(),
    ]);
    assert!(!success);
    let (success, _) = invoke(&[
        "fragment",
        input.to_str().unwrap(),
        "--scan",
        "2",
        "--proforma",
        "GAV",
        "--evidence-out",
        input.to_str().unwrap(),
    ]);
    assert!(!success);
    let mut unknown = output;
    for metadata in unknown["result"]["spectraMetadata"].as_array_mut().unwrap() {
        metadata["precursor"]["activation"] = Value::Null;
    }
    let path = workspace.save("unknown.json", &unknown);
    ok(&[
        "fragment",
        path.to_str().unwrap(),
        "--activation",
        "CID",
        "--proforma",
        "GAV",
        "--evidence-out",
        evidence.to_str().unwrap(),
    ]);
}

#[test]
fn per_observation_charge_and_export_limits_do_not_change_measurements() {
    let workspace = Workspace::new();
    let input = workspace.mzml("charges.mzML", 2, 2);
    let xml = std::fs::read_to_string(&input).unwrap();
    let boundary = xml.find("<spectrum id=\"scan=4\"").unwrap();
    let second = xml[boundary..].replace(
        "name=\"charge state\" value=\"1\"",
        "name=\"charge state\" value=\"2\"",
    );
    std::fs::write(&input, format!("{}{}", &xml[..boundary], second)).unwrap();
    let evidence = workspace.file("evidence.json");
    let models = r#"[{"id":"a","proforma":"GAV"}]"#;
    let output = ok(&[
        "fragment",
        input.to_str().unwrap(),
        "--scan",
        "2,4",
        "--candidate-proformas",
        models,
        "--evidence-out",
        evidence.to_str().unwrap(),
    ]);
    assert_eq!(output["result"]["comparisons"][0]["maxCharge"], 1);
    assert_eq!(output["result"]["comparisons"][1]["maxCharge"], 2);
    assert!(output["params"]["maxCharge"].is_null());
    assert_eq!(output["params"]["preparedTheorySets"], 2);
    let baseline = read(&evidence);
    ok(&[
        "fragment",
        input.to_str().unwrap(),
        "--scan",
        "2,4",
        "--candidate-proformas",
        models,
        "--max-matches",
        "1",
        "--evidence-out",
        evidence.to_str().unwrap(),
    ]);
    let limited = read(&evidence);
    for i in 0..2 {
        let all = &baseline["result"]["comparisons"][i]["measurement"];
        let part = &limited["result"]["comparisons"][i]["measurement"];
        assert_eq!(
            part["candidates"][0]["matches"].as_array().unwrap().len(),
            1
        );
        assert_eq!(part["candidates"][0]["matchesTruncated"], true);
        assert_eq!(
            part["candidates"][0]["matchedIons"],
            all["candidates"][0]["matchedIons"]
        );
        assert_eq!(
            part["candidates"][0]["medianAbsErrorPpm"],
            all["candidates"][0]["medianAbsErrorPpm"]
        );
        assert_eq!(part["discriminatingIons"], all["discriminatingIons"]);
    }
}

#[test]
fn file_backed_batch_inputs_and_early_parameter_errors() {
    let workspace = Workspace::new();
    let input = workspace.mzml("batch.mzML", 2, 2);
    let scans = workspace.save("scans.json", &json!([2, "scan=4"]));
    let models = workspace.save(
        "models.json",
        &json!([{"id":"a", "proforma":"GAV"}, {"id":"b", "proforma":"AV"}]),
    );
    let evidence = workspace.file("evidence.json");
    let result = ok(&[
        "fragment",
        input.to_str().unwrap(),
        "--scans-input",
        scans.to_str().unwrap(),
        "--candidate-proformas",
        models.to_str().unwrap(),
        "--evidence-out",
        evidence.to_str().unwrap(),
    ]);
    assert_eq!(result["result"]["comparisons"].as_array().unwrap().len(), 2);
    assert_eq!(result["params"]["preparedTheorySets"], 1);
    let pooled = ok(&[
        "spectrum",
        input.to_str().unwrap(),
        "--mode",
        "sum",
        "--scans-input",
        scans.to_str().unwrap(),
    ]);
    assert_eq!(
        pooled["result"]["spectraMetadata"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    for config in [
        vec!["--tolerance-ppm", "0"],
        vec!["--min-correlation", "2"],
        vec!["--max-charge", "0"],
    ] {
        let mut args = vec![
            "fragment",
            "file-does-not-exist.mzML",
            "--proforma",
            "GAV",
            "--scan",
            "2",
            "--evidence-out",
            evidence.to_str().unwrap(),
        ];
        args.extend(config);
        let (success, error) = invoke(&args);
        assert!(!success);
        assert_eq!(error["error"]["code"], "MASS_SPEC_INVALID_QUERY");
    }
    // Unknown model fields are rejected before reading the large source file.
    let (success, error) = invoke(&[
        "fragment",
        "file-does-not-exist.mzML",
        "--scan",
        "2",
        "--candidate-proformas",
        r#"[{"id":"a","proforma":"GAV","maxChrage":2}]"#,
        "--evidence-out",
        evidence.to_str().unwrap(),
    ]);
    assert!(!success);
    assert!(error["error"]["message"]
        .as_str()
        .unwrap()
        .contains("unknown field"));
    let db = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../system/reference/unimod.sqlite");
    let input =
        r#"[{"id":"q1","mass":20,"tolerance":0.01},{"id":"q2","mass":21,"tolerance":0.01}]"#;
    let searched = ok(&[
        "unimod-batch-search",
        "--input",
        input,
        "--db",
        db.to_str().unwrap(),
    ]);
    assert_eq!(searched["result"]["queries"].as_array().unwrap().len(), 2);
    let (success, error) = invoke(&[
        "refine-shift-batch",
        "file-does-not-exist.mzML",
        "--input",
        "[]",
    ]);
    assert!(!success);
    assert_eq!(error["error"]["code"], "MASS_SPEC_INVALID_QUERY");
}

#[test]
fn theory_fragment_theory_and_fits_batch_without_changing_single_results() {
    let workspace = Workspace::new();
    let models = workspace.save(
        "models.json",
        &json!([
            {"id":"a", "proforma":"GAV", "charge":2}, {"id":"b", "proforma":"AV", "charge":3}
        ]),
    );
    let batch = ok(&["theory", "--models", models.to_str().unwrap()]);
    for (index, (proforma, charge)) in [("GAV", "2"), ("AV", "3")].iter().enumerate() {
        let single = ok(&["theory", "--proforma", proforma, "--charge", charge]);
        assert_eq!(
            single["result"],
            batch["result"]["models"][index]["measurement"]
        );
    }
    let fragment_models = workspace.save(
        "fragment-models.json",
        &json!([
            {"id":"a", "proforma":"GAV"}, {"id":"b", "proforma":"AV"}
        ]),
    );
    let fragments = ok(&[
        "fragment-theoretical",
        "--models",
        fragment_models.to_str().unwrap(),
        "--ion-types",
        "B,Y",
    ]);
    for (index, proforma) in ["GAV", "AV"].iter().enumerate() {
        let single = ok(&[
            "fragment-theoretical",
            "--proforma",
            proforma,
            "--ion-types",
            "B,Y",
        ]);
        assert_eq!(
            single["result"],
            fragments["result"]["models"][index]["measurement"]
        );
    }
    let theory = workspace.save("theory.json", &batch["result"]["models"][0]["measurement"]);
    let observed: Vec<_> = batch["result"]["models"][0]["measurement"]["anchored"]["peaks"]
        .as_array()
        .unwrap()
        .iter()
        .map(|p| json!({"mz":p["mz"].as_f64().unwrap()+0.01,"intensity":p["intensity"]}))
        .collect();
    let observed = workspace.save("observed.json", &json!(observed));
    let pairs = workspace.save(
        "pairs.json",
        &json!([
            {"id":"legacy","observed":observed,"theoretical":theory},
            {"id":"strict","observed":observed,"theoretical":theory,"minimumWindowMz":0}
        ]),
    );
    let fitted = ok(&["fit-envelope", "--pairs", pairs.to_str().unwrap()]);
    let single = ok(&[
        "fit-envelope",
        "--observed",
        observed.to_str().unwrap(),
        "--theoretical",
        theory.to_str().unwrap(),
    ]);
    assert_eq!(
        fitted["result"]["fits"][0]["evaluation"],
        single["result"]["evaluation"]
    );
    assert!(
        fitted["result"]["fits"][0]["evaluation"]["matchedPeaks"]
            .as_u64()
            .unwrap()
            > 0
    );
    assert_eq!(fitted["result"]["fits"][1]["evaluation"]["matchedPeaks"], 0);
    assert_eq!(fitted["source"].as_array().unwrap().len(), 2);
    // A typo in a later row is a configuration error, never a defaulted value.
    let (success, error) = invoke(&[
        "theory",
        "--models",
        r#"[{"proforma":"GAV"},{"proforma":"AV","chagre":3}]"#,
    ]);
    assert!(!success);
    assert!(error["error"]["message"]
        .as_str()
        .unwrap()
        .contains("unknown field"));
}
