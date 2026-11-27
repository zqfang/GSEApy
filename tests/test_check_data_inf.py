import numpy as np
import pandas as pd

import gseapy as gp
from gseapy.base import GSEAbase


def test_check_data_replaces_inf_per_column():
    df = pd.DataFrame(
        {
            "gene": ["G1", "G2", "G3", "G4", "G5"],
            "s1": [np.inf, 2.0, 3.0, 5.0, -np.inf],
            "s2": [np.inf, 10.0, 20.0, 30.0, -np.inf],
            "s3": [1, 2, 3, 4, 5],
            "s4": [np.inf, -10.0, -5.0, -2.0, -np.inf],
        }
    )
    res = GSEAbase(outdir=None)._check_data(df)

    expected = pd.DataFrame(
        {
            "s1": [5.0, 2.0, 3.0, 5.0, 2.0],
            "s2": [30.0, 10.0, 20.0, 30.0, 10.0],
            "s3": [1, 2, 3, 4, 5],
            "s4": [-2.0, -10.0, -5.0, -2.0, -10.0],
        },
        index=pd.Index(["G1", "G2", "G3", "G4", "G5"], name="gene"),
    )
    pd.testing.assert_frame_equal(res, expected)


def test_check_data_replaces_inf_single_column():
    df = pd.DataFrame({"gene": ["G1", "G2", "G3", "G4"], "s1": [np.inf, 2.0, 5.0, -np.inf]})
    res = GSEAbase(outdir=None)._check_data(df)
    assert res["s1"].tolist() == [5.0, 2.0, 5.0, 2.0]


def test_ssgsea_inf_preserves_sample_ranking():
    genes = [f"G{i}" for i in range(1, 21)]
    df_inf = pd.DataFrame(
        {
            "s1": [np.inf] + list(np.linspace(5.0, 2.0, 18)) + [-np.inf],
            "s2": [np.inf] + list(np.linspace(30.0, 10.0, 18)) + [-np.inf],
        },
        index=genes,
    )
    df_finite = pd.DataFrame(
        {
            "s1": [5.0] + list(np.linspace(5.0, 2.0, 18)) + [2.0],
            "s2": [30.0] + list(np.linspace(30.0, 10.0, 18)) + [10.0],
        },
        index=genes,
    )
    gs = {"TOP": genes[:5], "BOT": genes[-5:]}

    res_inf = gp.ssgsea(data=df_inf, gene_sets=gs, sample_norm_method=None, min_size=3, outdir=None, no_plot=True)
    res_fin = gp.ssgsea(data=df_finite, gene_sets=gs, sample_norm_method=None, min_size=3, outdir=None, no_plot=True)

    pd.testing.assert_frame_equal(res_inf.ranking, res_fin.ranking)
    pd.testing.assert_frame_equal(res_inf.res2d, res_fin.res2d)
