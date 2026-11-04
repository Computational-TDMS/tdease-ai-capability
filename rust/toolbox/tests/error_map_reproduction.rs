use std::collections::BTreeSet;
use std::path::Path;
use std::process::Command;

const BEFORE_Y88: &str =
    "[Acetyl]-SGRGKGGKGLGKGGAKRHRK[Dimethyl]VLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAM[Oxidation]DVVY";
const AFTER_Y88: &str = "ALKRQGRTLYGFGG";

fn run_candidate(chemistry: &str) -> serde_json::Value {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures")
        .join("histone48-scan1164.mzML");
    let evidence = std::env::temp_dir().join(format!(
        "tdease-error-map-{}-{}.json",
        chemistry.to_ascii_lowercase(),
        std::process::id()
    ));
    let proforma = format!("{BEFORE_Y88}[{chemistry}]{AFTER_Y88}");
    let output = Command::new(env!("CARGO_BIN_EXE_tdeasenexus-toolbox"))
        .args([
            "fragment",
            fixture.to_str().unwrap(),
            "--mode",
            "match",
            "--scan",
            "1164",
            "--proforma",
            &proforma,
            "--evidence-out",
            evidence.to_str().unwrap(),
            "--ion-types",
            "C,Z",
            "--max-charge",
            "14",
            "--tolerance-ppm",
            "10",
        ])
        .output()
        .expect("tdeasenexus-toolbox runs");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let result = serde_json::from_slice::<serde_json::Value>(
        &std::fs::read(&evidence).expect("evidence sidecar exists"),
    )
    .expect("evidence sidecar is JSON")["result"]["measurement"]
        .clone();
    std::fs::remove_file(evidence).expect("temporary evidence is removable");
    result
}

fn matched_z_numbers(result: &serde_json::Value) -> BTreeSet<u64> {
    result["matches"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|matched| matched["ion"] == "z")
        .map(|matched| matched["number"].as_u64().unwrap())
        .collect()
}

#[test]
fn wilson_y88_phospho_has_the_missing_z_segment_that_bromo_retains() {
    let bromo = run_candidate("Bromo");
    let phospho = run_candidate("Phospho");

    assert_eq!(bromo["matches"].as_array().unwrap().len(), 151);
    assert_eq!(phospho["matches"].as_array().unwrap().len(), 170);
    assert_eq!(bromo["fragmentStatuses"].as_array().unwrap().len(), 2828);
    assert_eq!(phospho["fragmentStatuses"].as_array().unwrap().len(), 2828);

    let bromo_z = matched_z_numbers(&bromo);
    let phospho_z = matched_z_numbers(&phospho);
    assert!((17..=29).all(|number| !phospho_z.contains(&number)));
    assert!([19, 23, 24, 25, 27, 29]
        .into_iter()
        .all(|number| bromo_z.contains(&number)));
}
