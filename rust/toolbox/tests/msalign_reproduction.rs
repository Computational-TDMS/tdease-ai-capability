//! E2E reproduction of the deconvoluted-msalign fragment annotation on the
//! real scan: TopFD `_ms2.msalign` peaks are isotope-collapsed observations,
//! so matching is a charge-aware one-to-one assignment in neutral-mass space
//! (a DIFFERENT observation domain from the raw-mzML error map). This test
//! pins that domain contract end to end: the fixture is the exact scan-1164
//! block of MZ20160222DS_histone48_ms2.msalign, the models are the complete
//! Y88 candidates, and the counts below are the recorded Wilson-case values.

use std::collections::BTreeSet;
use std::path::Path;
use std::process::Command;

const BEFORE_Y88: &str =
    "[Acetyl]-SGRGKGGKGLGKGGAKRHRK[Dimethyl]VLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAM[Oxidation]DVVY";
const AFTER_Y88: &str = "ALKRQGRTLYGFGG";
const PEAKS_IN_SCAN: usize = 189;

fn run_candidate(chemistry: &str) -> serde_json::Value {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures")
        .join("histone48-scan1164.msalign");
    let evidence = std::env::temp_dir().join(format!(
        "tdease-msalign-{}-{}.json",
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

fn positions(result: &serde_json::Value, ion: &str) -> BTreeSet<u64> {
    result["matches"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|m| m["ion"].as_str().map(|s| s.to_lowercase()) == Some(ion.to_owned()))
        .map(|m| m["number"].as_u64().unwrap())
        .collect()
}

#[test]
fn msalign_domain_matches_charge_aware_neutral_masses() {
    let bromo = run_candidate("Bromo");
    let bytes = std::fs::read(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/histone48-scan1164.msalign"),
    )
    .unwrap();
    let source = tdeasenexus_toolbox::parse_msalign_spectrum(&bytes, "1164").unwrap();
    let mut assigned = BTreeSet::new();
    for matched in bromo["matches"].as_array().unwrap() {
        let index = matched["observedPeakIndex"].as_u64().unwrap() as usize;
        assert!(
            assigned.insert(index),
            "one deconvoluted observation is used once"
        );
        let peak = &source.peaks[index];
        assert_eq!(matched["charge"].as_u64(), peak.charge.map(u64::from));
        assert!((matched["observedNeutralMass"].as_f64().unwrap() - peak.mass).abs() < 1e-8);
        assert!((matched["intensity"].as_f64().unwrap() - f64::from(peak.intensity)).abs() < 1e-3);
        assert!(
            matched.get("quantification").is_none(),
            "no second isotope integration"
        );
    }

    // Scan-1164 carries 189 deconvoluted peaks; the recorded Bromo match set.
    assert_eq!(bromo["matches"].as_array().unwrap().len(), 60);
    let c_positions = positions(&bromo, "c");
    assert_eq!(c_positions.len(), 42);
    // TopFD deconvolutes the ETD products as they are; the plain-z series
    // records no match on this scan in this domain.
    assert!(positions(&bromo, "z").is_empty());
    // Every non-matching theoretical carries the deconvoluted-domain state.
    let states: BTreeSet<String> = bromo["fragmentStatuses"]
        .as_array()
        .unwrap()
        .iter()
        .map(|s| s["state"].as_str().unwrap().to_owned())
        .collect();
    assert_eq!(states, BTreeSet::from(["missingPeak".to_owned()]));
}

#[test]
fn msalign_domain_rejects_variant_requests() {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures")
        .join("histone48-scan1164.msalign");
    let proforma = format!("{BEFORE_Y88}[Bromo]{AFTER_Y88}");

    // Series are decided by the fragmentation method; ZDot is not generated.
    let out = Command::new(env!("CARGO_BIN_EXE_tdeasenexus-toolbox"))
        .args([
            "fragment",
            fixture.to_str().unwrap(),
            "--mode",
            "match",
            "--scan",
            "1164",
            "--proforma",
            &proforma,
            "--ion-types",
            "C,ZDOT",
            "--tolerance-ppm",
            "10",
        ])
        .output()
        .expect("tdeasenexus-toolbox binary runs");
    assert!(!out.status.success(), "ZDot is refused");
    let err: serde_json::Value =
        serde_json::from_slice(&out.stderr).expect("errors are a JSON envelope on stderr");
    assert!(
        err["error"]["message"]
            .as_str()
            .unwrap_or("")
            .contains("decided"),
        "refusal explains the method-decided rule: {err}"
    );
}

#[test]
fn msalign_scan1164_fixture_holds_the_real_scan() {
    // The fixture is the extracted scan block: guard its identity so a
    // re-extraction cannot silently change what the counts above pin.
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures")
        .join("histone48-scan1164.msalign");
    let text = std::fs::read_to_string(fixture).unwrap();
    assert!(text.contains("SCANS=1164"));
    assert!(text.contains("ACTIVATION=ETD"));
    assert_eq!(
        text.lines().filter(|l| !l.trim().is_empty()).count()
            - text.lines().filter(|l| l.contains('=')).count()
            - 2,
        PEAKS_IN_SCAN,
        "peak line count of the real scan"
    );
}
