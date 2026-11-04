//! Golden regression: `fragment-theoretical --out-format mzpaf` emits standard
//! HUPO-PSI mzPAF fragment annotations, filtered to the requested ion series.
use tdeasenexus_toolbox::{theoretical_fragments_mzpaf, IonType, MassModification};

#[test]
fn mzpaf_output_is_standard_and_filtered() {
    let mzpaf = theoretical_fragments_mzpaf(
        "SGRGKGGKGL",
        &[MassModification {
            position: "nTerm".to_owned(),
            delta_mass: 42.010_565,
            formula: None,
        }],
        &[IonType::C, IonType::Z],
        1,
    )
    .unwrap();
    let lines: Vec<&str> = mzpaf.lines().collect();
    assert!(!lines.is_empty(), "mzPAF is not empty");
    // First line is c1, annotated with the N-term acetyl (the peptide is
    // acetylated) — c1 = Acetyl + Ser.
    assert_eq!(lines[0], "1@c1", "first annotation is c1");
    // Only c and z ions appear (no a/b/x/y/v/w).
    for line in &lines {
        let kind = line
            .split('@')
            .nth(1)
            .unwrap_or("")
            .chars()
            .next()
            .unwrap_or(' ');
        assert!(
            kind == 'c' || kind == 'z',
            "unexpected ion kind {kind}: {line}"
        );
    }
    // A z ion carries the -1H variant (z· radical → z-H).
    assert!(
        lines.iter().any(|l| l.contains("z") && l.contains("-1H")),
        "z ions carry -1H"
    );
    // Same fragment set as the JSON path (c1..c9 + z1..z9 = 18).
    let c_count = lines.iter().filter(|l| l.contains("@c")).count();
    let z_count = lines.iter().filter(|l| l.contains("@z")).count();
    assert_eq!(c_count, 9, "9 c ions");
    assert_eq!(z_count, 9, "9 z ions");
}
