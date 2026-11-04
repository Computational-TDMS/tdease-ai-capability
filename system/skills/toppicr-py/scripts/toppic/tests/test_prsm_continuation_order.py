"""Read real TopPIC continuation layout without changing observation order."""
import numpy as np
import pandas as pd

from toppic.io.read_toppic import read_toppic


def test_continuation_fill_precedes_any_sequence_sort(tmp_path):
    # R read_toppic uses fill(Scan(s):Proteoform, .direction='down') in file
    # order. A continuation lacks Proteoform itself, so grouping by that
    # missing value cannot recover the parent observation.
    values = {
        "Data file name": ["sample_ms2.msalign"] * 3,
        "#peaks": [10] * 3,
        "Protein accession": ["sp|P1|entry", "sp|P2|entry", "sp|P1|entry"],
        "Protein description": ["GN=g1", "GN=g2", "GN=g1"],
        "Proteoform": ["K.ZZ.A", np.nan, "K.AA.A"],
        "Scan(s)": [200., np.nan, 100.],
        "Feature intensity": [20., np.nan, 10.],
        "Special amino acids": [np.nan] * 3,
    }
    for col in ["Retention time", "Feature apex time", "Charge",
                "Precursor mass", "Adjusted precursor mass", "E-value",
                "#unexpected modifications", "#Protein hits", "MIScore"]:
        values[col] = [1., np.nan, 2.]
    pd.DataFrame(values).to_csv(tmp_path/'sample_ms2_toppic_prsm.tsv',sep='\t',index=False)
    actual = read_toppic(tmp_path, ['sample_ms2_toppic_prsm.tsv'])
    assert actual.scans.tolist() == [200., 200., 100.]
    assert actual.proteoform.tolist() == ["K.ZZ.A", "K.ZZ.A", "K.AA.A"]
    assert actual.feature_intensity.tolist() == [20., 20., 10.]
    assert actual.gene.tolist() == ["g1", "g2", "g1"]
