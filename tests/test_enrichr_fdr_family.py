import math
from fractions import Fraction

import numpy as np
import pytest
from scipy.stats import false_discovery_control

import gseapy as gp


def run(query, sets, background):
    return (
        gp.enrichr(gene_list=query, gene_sets=sets, background=background, outdir=None, no_plot=True)
        .results.set_index("Term")
        .sort_index()
    )


def exact_tail(N, m, k, x):
    return float(
        sum(
            (
                Fraction(math.comb(m, i) * math.comb(N - m, k - i), math.comb(N, k))
                for i in range(x, min(m, k) + 1)
                if 0 <= k - i <= N - m
            ),
            Fraction(),
        )
    )


@pytest.mark.parametrize("background_kind", ["list", "integer", "implicit"])
def test_fixed_family_singleton_adjustment(background_kind):
    bg = [f"G{i}" for i in range(100)]
    sets = {f"S{i}": [g] for (i, g) in enumerate(bg)}
    background = {"list": bg, "integer": 100, "implicit": None}[background_kind]
    result = run(["G0"], sets, background)
    assert result.loc["S0", "P-value"] == pytest.approx(0.01)
    assert result.loc["S0", "Adjusted P-value"] == pytest.approx(1.0)
    assert len(result) == 1


def test_exhaustive_all_null_single_gene_queries():
    bg = [f"G{i}" for i in range(100)]
    sets = {f"S{i}": [g] for (i, g) in enumerate(bg)}
    rejected = sum((bool((run([gene], sets, bg)["Adjusted P-value"] <= 0.05).any()) for gene in bg))
    assert rejected == 0, f"{rejected}/100 equiprobable all-null queries reject"


@pytest.mark.parametrize("seed", range(8))
def test_exact_hypergeometric_and_full_family_bh(seed):
    rng = np.random.default_rng(seed)
    bg = [f"G{i}" for i in range(60)]
    query = list(rng.choice(bg, size=5, replace=False))
    sets = {f"T{i}": list(rng.choice(bg, size=int(rng.integers(1, 9)), replace=False)) for i in range(40)}
    sets["guaranteed_hit"] = query[:2]
    p = {
        term: exact_tail(60, len(set(genes)), len(query), len(set(genes) & set(query)))
        for (term, genes) in sets.items()
    }
    terms = sorted(p)
    expected = dict(zip(terms, false_discovery_control([p[t] for t in terms], method="bh")))
    result = run(query, sets, bg)
    for term, row in result.iterrows():
        assert row["P-value"] == pytest.approx(p[term])
        assert row["Adjusted P-value"] == pytest.approx(expected[term])


def test_all_overlapping_family_control():
    sets = {"A": ["G0"], "B": ["G0", "G1"], "C": ["G0", "G2", "G3"]}
    result = run(["G0"], sets, 100)
    np.testing.assert_allclose(result["Adjusted P-value"], false_discovery_control([0.01, 0.02, 0.03]))


def test_empty_background_term_does_not_change_family():
    result = run(["G0"], {"H": ["G0"], "outside": ["X"]}, ["G0", "G1"])
    assert result.loc["H", "Adjusted P-value"] == pytest.approx(0.5)


def test_no_overlaps_preserves_empty_public_result():
    result = gp.enrichr(
        gene_list=["G0"], gene_sets={"other": ["G1"]}, background=["G0", "G1"], outdir=None, no_plot=True
    )
    assert result.results == []


def test_gmt_file_keeps_zero_hit_terms_in_family(tmp_path):
    gmt = tmp_path / "fixed.gmt"
    gmt.write_text("".join(f"S{i}\tna\tG{i}\n" for i in range(100)))
    result = run(["G0"], str(gmt), 100)
    assert result.loc["S0", "Adjusted P-value"] == pytest.approx(1.0)


def test_each_local_library_retains_its_own_family():
    libraries = [{"S0": ["G0"]}, {f"S{i}": [f"G{i}"] for i in range(100)}]
    result = gp.enrichr(gene_list=["G0"], gene_sets=libraries, background=100, outdir=None, no_plot=True).results
    np.testing.assert_allclose(result["P-value"], [0.01, 0.01])
    np.testing.assert_allclose(result["Adjusted P-value"], [0.01, 1.0])
