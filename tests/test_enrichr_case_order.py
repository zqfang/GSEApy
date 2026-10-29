import pytest

import gseapy as gp


def run(query, sets, background):
    return (
        gp.enrichr(gene_list=query, gene_sets=sets, background=background, outdir=None, no_plot=True)
        .results.set_index("Term")
        .sort_index()
    )


@pytest.mark.parametrize("position", range(11))
def test_mixed_case_library_order_does_not_change_exact_match(position):
    upper = [(f"U{i}", ["ACTB"]) for i in range(10)]
    upper.insert(position, ("M", ["Actb"]))
    result = run(["Actb"], dict(upper), ["ACTB", "Actb", "OTHER"])
    assert result.index.tolist() == ["M"]
    assert result.loc["M", "P-value"] == pytest.approx(1 / 3)


def test_uniform_uppercase_conversion_control():
    result = run(["Actb"], {"U": ["ACTB"]}, ["ACTB", "OTHER"])
    assert result.loc["U", "P-value"] == pytest.approx(0.5)
    assert result.loc["U", "Genes"] == "ACTB"


def test_uniform_mixed_case_exact_matching_control():
    result = run(["Actb"], {"M": ["Actb"]}, ["Actb", "Other"])
    assert result.loc["M", "P-value"] == pytest.approx(0.5)


@pytest.mark.parametrize("mixed_first", [False, True])
def test_gmt_file_order_preserves_exact_matching(tmp_path, mixed_first):
    entries = [(f"U{i}", ["ACTB"]) for i in range(10)]
    entries.insert(0 if mixed_first else 10, ("M", ["Actb"]))
    gmt = tmp_path / "case.gmt"
    gmt.write_text("".join(f"{term}\tna\t{genes[0]}\n" for term, genes in entries))
    result = run(["Actb"], str(gmt), ["ACTB", "Actb", "OTHER"])
    assert result.index.tolist() == ["M"]
    assert result.loc["M", "P-value"] == pytest.approx(1 / 3)


def test_empty_term_does_not_prevent_uppercase_conversion():
    result = run(["Actb"], {"empty": [], "U": ["ACTB"]}, ["ACTB", "OTHER"])
    assert result.index.tolist() == ["U"]
    assert result.loc["U", "P-value"] == pytest.approx(0.5)
