import sys

import numpy as np

from gseapy import enrich

# gseapy/__init__.py defines a top-level `enrichr` function that shadows the
# submodule attribute, so `import gseapy.enrichr as m` binds the function.
enrichr_mod = sys.modules["gseapy.enrichr"]

SMALLEST_SUBNORMAL = np.nextafter(0, 1)  # 5e-324


def _fake_calc_pvalues(**kwargs):
    """Force one p-value to exactly 0.0 and one to the smallest subnormal.

    hypergeom.sf() underflows to 0.0 for extreme overlaps, passing through the
    subnormal range on the way. Faking it keeps the test fast and exact.
    Same output shape as calc_pvalues: a zip of columns.
    """
    vals = [
        ("TERM_ZERO", 0.0, 5.0, 3, 5, {"A", "B", "C"}),
        ("TERM_SUBNORMAL", SMALLEST_SUBNORMAL, 4.0, 3, 5, {"A", "B", "C"}),
        ("TERM_NORMAL", 0.2, 1.5, 2, 5, {"D", "E"}),
    ]
    return zip(*vals)


def _run(monkeypatch):
    monkeypatch.setattr(enrichr_mod, "calc_pvalues", _fake_calc_pvalues)
    return enrich(
        gene_list=["A", "B", "C", "D", "E"],
        gene_sets={"dummy": ["A", "B", "C", "D", "E"]},
        background=None,
        outdir=None,
        cutoff=1.0,
        no_plot=True,
    )


def test_combined_score_finite_when_pvalue_underflows_to_zero(monkeypatch):
    """A p-value of 0.0 must not make Combined Score inf."""
    with np.errstate(divide="raise"):
        result = _run(monkeypatch)

    combined = result.res2d["Combined Score"]
    assert np.isfinite(combined).all()

    row = result.res2d[result.res2d["Term"] == "TERM_ZERO"]
    assert np.isfinite(row["Combined Score"].iloc[0])
    assert row["P-value"].iloc[0] == 0.0


def test_subnormal_pvalue_is_not_floored(monkeypatch):
    """A real subnormal p-value must reach the log untouched.

    Flooring at np.finfo(float).tiny (2.2e-308, the smallest normal) instead
    of the smallest subnormal would round these up and shift the score.
    """
    result = _run(monkeypatch)

    row = result.res2d[result.res2d["Term"] == "TERM_SUBNORMAL"]
    assert row["Combined Score"].iloc[0] == -1 * np.log(SMALLEST_SUBNORMAL) * 4.0
    assert row["P-value"].iloc[0] == SMALLEST_SUBNORMAL
