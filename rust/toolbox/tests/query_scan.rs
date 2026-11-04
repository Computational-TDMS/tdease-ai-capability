use tdeasenexus_toolbox::{query_mzml, PeakOrder, ScanQuery, SpectrumSelector};

const FIXTURE: &[u8] = include_bytes!("../fixtures/single-spectrum.mzML");

#[test]
fn scan_query_filters_and_caps_peaks_through_the_public_operator() {
    let result = query_mzml(
        FIXTURE,
        &ScanQuery {
            selector: SpectrumSelector::Index(0),
            mz_min: Some(150.0),
            mz_max: Some(600.0),
            min_intensity: Some(2_500.0),
            max_peaks: Some(2),
            peak_order: PeakOrder::IntensityDescending,
            ms_level: None,
        },
    )
    .expect("fixture query should succeed");

    assert_eq!(result.id, "scan=1");
    assert_eq!(result.index, 0);
    assert_eq!(result.ms_level, 1);
    assert_eq!(result.total_peak_count, 5);
    assert_eq!(result.matched_peak_count, 3);
    assert!(result.truncated);
    assert_eq!(
        result.peaks,
        vec![(500.2, 20_000.0).into(), (200.0, 5_000.0).into()]
    );
}
