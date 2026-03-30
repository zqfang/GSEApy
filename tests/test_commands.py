from tempfile import TemporaryDirectory, mkdtemp

import numpy as np
import pandas as pd
import pytest

from gseapy.__init__ import enrich, enrichr, gsea, gsva, prerank, replot, ssgsea


@pytest.fixture
def edbDIR():
    return "tests/data"


@pytest.fixture
def genelist():
    return "tests/data/gene_list.txt"


@pytest.fixture
def gseaGCT():
    return "tests/extdata/Leukemia_hgu95av2.trim.txt"


@pytest.fixture
def gseaCLS():
    return "tests/extdata/Leukemia.cls"


@pytest.fixture
def prernk():
    return "tests/data/edb/gsea_data.gsea_data.rnk"


@pytest.fixture
def geneGMT():
    return "tests/data/genes.gmt"


@pytest.fixture
def ssGMT():
    return "tests/data/temp.gmt"


@pytest.fixture
def ssGCT():
    return "tests/extdata/Leukemia_hgu95av2.trim.txt"


@pytest.fixture
def gsvaExpr():
    return "tests/data/expr.gsva.csv"


@pytest.fixture
def gsvaGMT():
    return "tests/data/geneset.gsva.gmt"


def test_gsea(gseaGCT, gseaCLS, geneGMT):
    # Only tests of the command runs successfully,
    # doesnt't check the image
    tmpdir = TemporaryDirectory(dir="tests")
    gs_res = gsea(
        data=gseaGCT,
        gene_sets=geneGMT,
        cls=gseaCLS,
        method="t_test",
        outdir=tmpdir.name,
        permutation_type="phenotype",
        permutation_num=47,
    )
    tmpdir.cleanup()


def test_fdr_gsea(gseaGCT, gseaCLS, geneGMT):
    # Runs and verifies a reasonable result for pval
    # when method="t_test" and permutation_type="gene_set"
    tmpdir = TemporaryDirectory(dir="tests")
    gs_res = gsea(
        data=gseaGCT,
        gene_sets=geneGMT,
        cls=gseaCLS,
        method="t_test",
        outdir=tmpdir.name,
        permutation_type="gene_set",
        permutation_num=47,
        seed=7,
    )
    assert gs_res.res2d["NOM p-val"].max() > 0.0
    tmpdir.cleanup()


def test_prerank(prernk, geneGMT):
    # Only tests of the command runs successfully,
    # doesnt't check the image
    tmpdir = TemporaryDirectory(dir="tests")
    prerank(prernk, geneGMT, outdir=tmpdir.name, permutation_num=10)
    tmpdir.cleanup()
    prerank(
        prernk, ["KEGG_2016", "GO_Biological_Process_2021"], outdir=None, permutation_num=20
    )


def test_prerank_fdr_per_library(prernk):
    """FDR values for a term should be the same whether the gene set library
    is provided alone or combined with other libraries (issue: FDR per library)."""
    # Build two gene set dicts using genes from the rnk file so they overlap
    rnk = pd.read_csv(prernk, header=None, names=["gene", "score"], sep="\t")
    genes = rnk["gene"].tolist()
    # Create two non-overlapping gene sets from the ranked list
    lib1 = {"SetA": genes[0:20], "SetB": genes[20:40]}
    lib2 = {"SetC": genes[200:220], "SetD": genes[220:240]}

    # Run with lib1 only (as list of one dict, so terms get "0__" prefix)
    res_single = prerank(prernk, [lib1], outdir=None, permutation_num=100, seed=42)

    # Run with lib1 + lib2 combined (terms from lib1 get "0__" prefix)
    res_multi = prerank(
        prernk, [lib1, lib2], outdir=None, permutation_num=100, seed=42
    )

    single_df = res_single.res2d.copy()
    multi_df = res_multi.res2d.copy()

    # Terms from lib1 are prefixed "0__" in both single and multi runs
    single_terms = set(single_df["Term"])
    multi_lib1 = multi_df[multi_df["Term"].str.startswith("0__")].copy()
    multi_terms = set(multi_lib1["Term"])

    common_terms = single_terms & multi_terms
    assert len(common_terms) > 0, "No common terms found between single and multi runs"

    single_fdr = single_df.set_index("Term")["FDR q-val"]
    multi_fdr = multi_lib1.set_index("Term")["FDR q-val"]

    for term in common_terms:
        assert np.isclose(single_fdr[term], multi_fdr[term], atol=1e-6), (
            f"FDR mismatch for term '{term}': "
            f"single={single_fdr[term]}, multi={multi_fdr[term]}"
        )


def test_ssgsea1(ssGCT, geneGMT):
    # Only tests of the command runs successfully,
    # doesnt't check the image
    tmpdir = TemporaryDirectory(dir="tests")
    ssgsea(ssGCT, geneGMT, outdir=tmpdir.name, permutation_num=100)
    tmpdir.cleanup()


def test_ssgsea2(ssGCT, geneGMT):
    # Only tests of the command runs successfully,
    # doesnt't check the image
    tmpdir = TemporaryDirectory(dir="tests")
    ssgsea(ssGCT, geneGMT, outdir=tmpdir.name, permutation_num=0)
    tmpdir.cleanup()
    ssgsea(ssGCT, geneGMT, outdir=None, permutation_num=0, correl_norm_type="symrank")
    ssgsea(ssGCT, geneGMT, outdir=None, permutation_num=0, correl_norm_type="zscore")


def test_gsva(gsvaExpr, gsvaGMT):
    # Only tests of the command runs successfully,
    # doesnt't check the image
    tmpdir = TemporaryDirectory(dir="tests")
    gsva(gsvaExpr, gsvaGMT, outdir=tmpdir.name, kcdf="Gaussian")
    gsva(gsvaExpr, gsvaGMT, outdir=tmpdir.name, kcdf="Poisson")
    gsva(gsvaExpr, gsvaGMT, outdir=tmpdir.name, kcdf=None)
    tmpdir.cleanup()


def test_gsva_uppercase_conversion():
    """Test GSVA when gene sets are uppercase but expression genes are lowercase.

    Regression test for: 'numpy.float64' object has no attribute 'upper'
    This ensures str() is called before .upper() on gene names.
    """
    np.random.seed(42)
    n_genes = 100
    n_samples = 4
    gene_names = ["gene%d" % i for i in range(n_genes)]
    expr = pd.DataFrame(
        np.random.randn(n_genes, n_samples),
        index=gene_names,
        columns=["s%d" % i for i in range(n_samples)],
    )
    # gene sets with UPPERCASE gene names trigger the _gene_toupper path
    gene_sets = {
        "SET_A": [g.upper() for g in gene_names[:30]],
        "SET_B": [g.upper() for g in gene_names[20:60]],
    }
    tmpdir = TemporaryDirectory(dir="tests")
    result = gsva(data=expr, gene_sets=gene_sets, outdir=tmpdir.name)
    assert result.res2d is not None
    tmpdir.cleanup()


def test_enrichr(genelist, geneGMT):
    # Only tests of the command runs successfully,
    # doesnt't check the image
    tmpdir = TemporaryDirectory(dir="tests")
    enrichr(
        genelist,
        gene_sets=["KEGG_2016", geneGMT],
        background="hsapiens_gene_ensembl",
        outdir=None,
        cutoff=0.5,
    )
    tmpdir.cleanup()
    tmpdir = TemporaryDirectory(dir="tests")
    enrich(genelist, gene_sets=geneGMT, background=None, outdir=tmpdir.name, cutoff=0.5)
    tmpdir.cleanup()


# ---------------------------------------------------------------------------
# Enrichr / enrich — additional test cases
# ---------------------------------------------------------------------------


@pytest.fixture
def small_gene_sets():
    """Small in-memory gene sets dict for local enrichment tests."""
    return {
        "PATHWAY_A": ["IGKC", "CD55", "ABHD4", "PCSK6", "PGD", "ARHGDIB"],
        "PATHWAY_B": ["ITGB2", "CARD6", "MNDA", "ATE1", "LGALS8", "HOMER3"],
        "PATHWAY_C": ["VDAC3", "MAP4K4", "HS2ST1", "FAKE1", "FAKE2", "FAKE3"],
    }


@pytest.fixture
def small_gene_list():
    """A small gene list that partially overlaps with small_gene_sets."""
    return ["IGKC", "CD55", "ABHD4", "ITGB2", "CARD6", "VDAC3", "MAP4K4", "NOTFOUND"]


class TestEnrichLocal:
    """Tests for offline / local enrichment (enrich function with dict gene sets)."""

    def test_enrich_dict_gene_sets(self, small_gene_list, small_gene_sets):
        """enrich() with dict gene sets should return results with expected columns."""
        tmpdir = TemporaryDirectory(dir="tests")
        result = enrich(
            gene_list=small_gene_list,
            gene_sets=small_gene_sets,
            background=None,
            outdir=tmpdir.name,
            cutoff=1.0,
            no_plot=True,
        )
        assert result.res2d is not None
        assert not result.res2d.empty
        expected_cols = {
            "Gene_set",
            "Term",
            "Overlap",
            "P-value",
            "Adjusted P-value",
            "Odds Ratio",
            "Combined Score",
            "Genes",
        }
        assert expected_cols.issubset(set(result.res2d.columns))
        tmpdir.cleanup()

    def test_enrich_returns_all_terms(self, small_gene_list, small_gene_sets):
        """All terms with at least one hit should appear in results."""
        result = enrich(
            gene_list=small_gene_list,
            gene_sets=small_gene_sets,
            background=None,
            outdir=None,
            cutoff=1.0,
            no_plot=True,
        )
        returned_terms = set(result.res2d["Term"])
        # Each pathway has at least some overlap with small_gene_list
        for term in small_gene_sets:
            assert term in returned_terms

    def test_enrich_overlap_format(self, small_gene_list, small_gene_sets):
        """Overlap column should be in 'k/K' format."""
        result = enrich(
            gene_list=small_gene_list,
            gene_sets=small_gene_sets,
            background=None,
            outdir=None,
            cutoff=1.0,
            no_plot=True,
        )
        for overlap in result.res2d["Overlap"]:
            parts = overlap.split("/")
            assert len(parts) == 2
            assert int(parts[0]) > 0
            assert int(parts[1]) > 0
            assert int(parts[0]) <= int(parts[1])

    def test_enrich_pvalue_range(self, small_gene_list, small_gene_sets):
        """P-values and adjusted p-values should be in [0, 1]."""
        result = enrich(
            gene_list=small_gene_list,
            gene_sets=small_gene_sets,
            background=None,
            outdir=None,
            cutoff=1.0,
            no_plot=True,
        )
        assert (result.res2d["P-value"] >= 0).all()
        assert (result.res2d["P-value"] <= 1).all()
        assert (result.res2d["Adjusted P-value"] >= 0).all()
        assert (result.res2d["Adjusted P-value"] <= 1).all()

    def test_enrich_with_integer_background(self, small_gene_list, small_gene_sets):
        """enrich() should work when background is an integer."""
        result = enrich(
            gene_list=small_gene_list,
            gene_sets=small_gene_sets,
            background=20000,
            outdir=None,
            cutoff=1.0,
            no_plot=True,
        )
        assert result.res2d is not None
        assert not result.res2d.empty

    def test_enrich_with_list_background(self, small_gene_list, small_gene_sets):
        """enrich() should work when background is a list of genes."""
        bg = list(
            {g for genes in small_gene_sets.values() for g in genes}
            | set(small_gene_list)
            | {"EXTRA1", "EXTRA2", "EXTRA3", "EXTRA4", "EXTRA5"}
        )
        result = enrich(
            gene_list=small_gene_list,
            gene_sets=small_gene_sets,
            background=bg,
            outdir=None,
            cutoff=1.0,
            no_plot=True,
        )
        assert result.res2d is not None
        assert not result.res2d.empty

    def test_enrich_with_background_file(self, small_gene_list, small_gene_sets):
        """enrich() should work when background is a file path."""
        result = enrich(
            gene_list=small_gene_list,
            gene_sets=small_gene_sets,
            background="tests/data/background.txt",
            outdir=None,
            cutoff=1.0,
            no_plot=True,
        )
        assert result.res2d is not None

    def test_enrich_genes_column_contains_hits(self, small_gene_list, small_gene_sets):
        """Genes column should only contain genes from the input list."""
        result = enrich(
            gene_list=small_gene_list,
            gene_sets=small_gene_sets,
            background=None,
            outdir=None,
            cutoff=1.0,
            no_plot=True,
        )
        input_set = set(small_gene_list)
        for genes_str in result.res2d["Genes"]:
            genes = set(genes_str.split(";"))
            assert genes.issubset(input_set)

    def test_enrich_gmt_file(self, genelist, geneGMT):
        """enrich() should work with a GMT file path."""
        result = enrich(
            gene_list=genelist,
            gene_sets=geneGMT,
            background=None,
            outdir=None,
            cutoff=1.0,
            no_plot=True,
        )
        assert result.res2d is not None

    def test_enrich_no_outdir(self, small_gene_list, small_gene_sets):
        """enrich() should work when outdir=None (no file output)."""
        result = enrich(
            gene_list=small_gene_list,
            gene_sets=small_gene_sets,
            background=None,
            outdir=None,
            cutoff=1.0,
            no_plot=True,
        )
        assert result.res2d is not None

    def test_enrich_cutoff_filters_plot_not_data(
        self, small_gene_list, small_gene_sets
    ):
        """Cutoff should not remove rows from res2d (only affects plotting)."""
        result_strict = enrich(
            gene_list=small_gene_list,
            gene_sets=small_gene_sets,
            background=None,
            outdir=None,
            cutoff=0.0001,
            no_plot=True,
        )
        result_loose = enrich(
            gene_list=small_gene_list,
            gene_sets=small_gene_sets,
            background=None,
            outdir=None,
            cutoff=1.0,
            no_plot=True,
        )
        # Both should contain the same number of terms
        assert len(result_strict.res2d) == len(result_loose.res2d)


class TestEnrichInputFormats:
    """Tests for various input formats accepted by enrich/enrichr."""

    def test_enrich_list_input(self, small_gene_sets):
        """enrich() should accept a plain list of gene names."""
        result = enrich(
            gene_list=["IGKC", "CD55", "ABHD4"],
            gene_sets=small_gene_sets,
            outdir=None,
            cutoff=1.0,
            no_plot=True,
        )
        assert result.res2d is not None

    def test_enrich_series_input(self, small_gene_sets):
        """enrich() should accept a pd.Series of gene names."""
        gene_series = pd.Series(["IGKC", "CD55", "ABHD4", "ITGB2"])
        result = enrich(
            gene_list=gene_series,
            gene_sets=small_gene_sets,
            outdir=None,
            cutoff=1.0,
            no_plot=True,
        )
        assert result.res2d is not None

    def test_enrich_dataframe_input(self, small_gene_sets):
        """enrich() should accept a single-column pd.DataFrame of gene names."""
        gene_df = pd.DataFrame({"genes": ["IGKC", "CD55", "ABHD4", "ITGB2"]})
        result = enrich(
            gene_list=gene_df,
            gene_sets=small_gene_sets,
            outdir=None,
            cutoff=1.0,
            no_plot=True,
        )
        assert result.res2d is not None

    def test_enrich_file_input(self, small_gene_sets):
        """enrich() should accept a file path to a gene list."""
        result = enrich(
            gene_list="tests/data/gene_list.txt",
            gene_sets=small_gene_sets,
            outdir=None,
            cutoff=1.0,
            no_plot=True,
        )
        assert result.res2d is not None

    def test_enrich_comma_separated_gene_sets(self, genelist):
        """enrich() should accept comma-separated GMT file paths."""
        gmt = "tests/data/genes.gmt"
        result = enrich(
            gene_list=genelist,
            gene_sets=gmt,
            outdir=None,
            cutoff=1.0,
            no_plot=True,
        )
        assert result.res2d is not None


class TestEnrichUpperCase:
    """Tests for gene name case handling in local enrichment."""

    def test_enrich_uppercase_gene_sets_lowercase_input(self):
        """When gene sets are UPPERCASE and input is lowercase, should still match."""
        gene_list = ["gene1", "gene2", "gene3", "gene4"]
        gene_sets = {
            "SET_A": ["GENE1", "GENE2", "GENE5", "GENE6"],
            "SET_B": ["GENE3", "GENE4", "GENE7", "GENE8"],
        }
        result = enrich(
            gene_list=gene_list,
            gene_sets=gene_sets,
            outdir=None,
            cutoff=1.0,
            no_plot=True,
        )
        assert result.res2d is not None
        assert not result.res2d.empty

    def test_enrich_mixed_case_query(self):
        """Mixed-case gene lists should still produce results."""
        gene_list = ["Gene1", "GENE2", "gene3"]
        gene_sets = {
            "SET_A": ["GENE1", "GENE2", "GENE3", "GENE4"],
        }
        result = enrich(
            gene_list=gene_list,
            gene_sets=gene_sets,
            outdir=None,
            cutoff=1.0,
            no_plot=True,
        )
        assert result.res2d is not None


class TestEnrichValidation:
    """Tests for input validation and error handling."""

    def test_enrich_empty_gene_list_raises(self, small_gene_sets):
        """enrich() should raise on empty gene list."""
        with pytest.raises(Exception):
            enrich(
                gene_list=[],
                gene_sets=small_gene_sets,
                outdir=None,
                no_plot=True,
            )

    def test_enrich_invalid_gene_sets_type_raises(self):
        """enrich() should raise when gene_sets is an invalid type."""
        with pytest.raises(Exception):
            enrich(
                gene_list=["GENE1", "GENE2"],
                gene_sets=12345,
                outdir=None,
                no_plot=True,
            )

class TestEnrichrAPI:
    """Tests for the EnrichrAPI class directly."""

    def test_enrichr_api_invalid_organism_raises(self):
        from gseapy.enrichr import EnrichrAPI

        with pytest.raises(ValueError, match="Invalid organism"):
            EnrichrAPI(organism="alien")

    def test_enrichr_api_ensure_list_string(self):
        from gseapy.enrichr import EnrichrAPI

        api = EnrichrAPI(organism="human")
        assert api._ensure_list("KEGG_2016") == ["KEGG_2016"]

    def test_enrichr_api_ensure_list_already_list(self):
        from gseapy.enrichr import EnrichrAPI

        api = EnrichrAPI(organism="human")
        assert api._ensure_list(["A", "B"]) == ["A", "B"]

    def test_enrichr_api_organism_aliases(self):
        """Different aliases for the same organism should resolve to the same URL."""
        from gseapy.enrichr import EnrichrAPI

        human_aliases = ["human", "hsapiens", "homo sapiens", "hs"]
        urls = {EnrichrAPI(organism=a).base_url for a in human_aliases}
        assert len(urls) == 1  # all should map to the same URL


class TestEnrichrClass:
    """Tests for the Enrichr class instantiation and helpers."""

    def test_enrichr_class_parse_gmt(self, geneGMT):
        from gseapy.enrichr import Enrichr

        enr = Enrichr(
            gene_list=["A", "B"],
            gene_sets=geneGMT,
            outdir=None,
            no_plot=True,
        )
        gmt = enr._parse_gmt(geneGMT)
        assert isinstance(gmt, dict)
        assert len(gmt) > 0
        for term, genes in gmt.items():
            assert isinstance(term, str)
            assert isinstance(genes, list)
            assert len(genes) > 0
        enr.close()

    def test_enrichr_class_check_uppercase(self):
        from gseapy.enrichr import Enrichr

        enr = Enrichr(
            gene_list=["A", "B"],
            gene_sets={"X": ["A"]},
            outdir=None,
            no_plot=True,
        )
        assert enr.check_uppercase(["BRCA1", "TP53", "EGFR"]) is True
        assert enr.check_uppercase(["brca1", "tp53", "egfr"]) is False
        assert enr.check_uppercase(["123", "456"]) is False  # entrez IDs
        enr.close()

    def test_enrichr_class_is_entrez_id(self):
        from gseapy.enrichr import Enrichr

        enr = Enrichr(
            gene_list=["A"],
            gene_sets={"X": ["A"]},
            outdir=None,
            no_plot=True,
        )
        assert enr._is_entrez_id("12345") is True
        assert enr._is_entrez_id(12345) is True
        assert enr._is_entrez_id("BRCA1") is False
        assert enr._is_entrez_id("") is False
        enr.close()

    def test_enrichr_class_filter_gmt(self):
        from gseapy.enrichr import Enrichr

        enr = Enrichr(
            gene_list=["A"],
            gene_sets={"X": ["A"]},
            outdir=None,
            no_plot=True,
        )
        gmt = {"TERM1": ["A", "B", "C"], "TERM2": ["D", "E"]}
        bg = {"A", "B", "D"}
        filtered = enr.filter_gmt(gmt, bg)
        assert filtered["TERM1"] == ["A", "B"]
        assert filtered["TERM2"] == ["D"]
        enr.close()

    def test_enrichr_class_filter_gmt_removes_empty(self):
        from gseapy.enrichr import Enrichr

        enr = Enrichr(
            gene_list=["A"],
            gene_sets={"X": ["A"]},
            outdir=None,
            no_plot=True,
        )
        gmt = {"TERM1": ["X", "Y"], "TERM2": ["A"]}
        bg = {"A"}
        filtered = enr.filter_gmt(gmt, bg)
        assert "TERM1" not in filtered
        assert filtered["TERM2"] == ["A"]
        enr.close()

    def test_enrichr_class_context_manager(self):
        from gseapy.enrichr import Enrichr

        with Enrichr(
            gene_list=["A", "B"],
            gene_sets={"X": ["A"]},
            outdir=None,
            no_plot=True,
        ) as enr:
            assert enr is not None

    def test_enrichr_class_parse_genelists_list(self):
        from gseapy.enrichr import Enrichr

        enr = Enrichr(
            gene_list=["BRCA1", "TP53", "EGFR"],
            gene_sets={"X": ["BRCA1"]},
            outdir=None,
            no_plot=True,
        )
        result = enr.parse_genelists()
        assert "BRCA1" in result
        assert "TP53" in result
        assert enr._isezid is False
        enr.close()

    def test_enrichr_class_parse_genelists_entrez(self):
        from gseapy.enrichr import Enrichr

        enr = Enrichr(
            gene_list=["672", "7157", "1956"],
            gene_sets={"X": ["672"]},
            outdir=None,
            no_plot=True,
        )
        result = enr.parse_genelists()
        assert enr._isezid is True
        assert isinstance(enr._gls, set)
        enr.close()

    def test_enrichr_class_parse_background_none(self):
        from gseapy.enrichr import Enrichr

        enr = Enrichr(
            gene_list=["A"],
            gene_sets={"X": ["A"]},
            outdir=None,
            background=None,
            no_plot=True,
        )
        gmt = {"TERM1": ["A", "B", "C"]}
        bg = enr.parse_background(gmt)
        assert isinstance(bg, set)
        assert bg == {"A", "B", "C"}
        enr.close()

    def test_enrichr_class_parse_background_int(self):
        from gseapy.enrichr import Enrichr

        enr = Enrichr(
            gene_list=["A"],
            gene_sets={"X": ["A"]},
            outdir=None,
            background=20000,
            no_plot=True,
        )
        bg = enr.parse_background()
        assert bg == 20000
        enr.close()

    def test_enrichr_class_parse_background_set(self):
        from gseapy.enrichr import Enrichr

        bg_genes = ["A", "B", "C", "D"]
        enr = Enrichr(
            gene_list=["A"],
            gene_sets={"X": ["A"]},
            outdir=None,
            background=bg_genes,
            no_plot=True,
        )
        bg = enr.parse_background()
        assert isinstance(bg, set)
        assert bg == {"A", "B", "C", "D"}
        enr.close()

    def test_extract_go_ids_with_go_id(self):
        """_extract_go_ids should return the GO ID when present in the term string."""
        from gseapy.enrichr import Enrichr

        enr = Enrichr(
            gene_list=["A"],
            gene_sets={"X": ["A"]},
            outdir=None,
            no_plot=True,
        )
        import pandas as pd

        terms = pd.Series(
            [
                "response to stimulus (GO:0050896)",
                "DNA repair (GO:0006281)",
                "cell cycle",
            ]
        )
        result = enr._extract_go_ids(terms)
        assert result[0] == "GO:0050896"
        assert result[1] == "GO:0006281"
        assert result[2] is None
        enr.close()

    def test_go_filter_no_go_ids_returns_unchanged(self):
        """go_filter should return the input unchanged when no GO IDs are found."""
        from gseapy.enrichr import Enrichr

        enr = Enrichr(
            gene_list=["A"],
            gene_sets={"X": ["A"]},
            outdir=None,
            no_plot=True,
        )
        import pandas as pd

        df = pd.DataFrame(
            {
                "Gene_set": ["KEGG"],
                "Term": ["Glycolysis"],
                "P-value": [0.01],
                "Adjusted P-value": [0.05],
                "Genes": ["GENE1;GENE2"],
            }
        )
        result = enr.go_filter(df=df, min_level=2, max_level=8)
        # No GO IDs → unchanged
        assert len(result) == len(df)
        enr.close()

    def test_go_filter_empty_df_returns_unchanged(self):
        """go_filter should return an empty DataFrame when given one."""
        from gseapy.enrichr import Enrichr

        enr = Enrichr(
            gene_list=["A"],
            gene_sets={"X": ["A"]},
            outdir=None,
            no_plot=True,
        )
        import pandas as pd

        empty_df = pd.DataFrame(columns=["Gene_set", "Term", "P-value"])
        result = enr.go_filter(df=empty_df)
        assert result.empty
        enr.close()

    def test_go_filter_with_mocked_levels(self):
        """go_filter should correctly filter by min/max level using mocked API."""
        from unittest.mock import patch

        import pandas as pd

        from gseapy.enrichr import Enrichr

        enr = Enrichr(
            gene_list=["A"],
            gene_sets={"X": ["A"]},
            outdir=None,
            no_plot=True,
        )
        df = pd.DataFrame(
            {
                "Gene_set": ["GO_BP"] * 3,
                "Term": [
                    # Levels below are mocked - not actual GO hierarchy depths
                    "biological_process (GO:0008150)",   # mocked level 0 (root)
                    "cellular process (GO:0009987)",     # mocked level 1
                    "DNA repair (GO:0006281)",            # mocked level 5
                ],
                "P-value": [0.001, 0.002, 0.003],
                "Adjusted P-value": [0.01, 0.02, 0.03],
                "Genes": ["G1", "G2", "G3"],
            }
        )
        mocked_levels = {
            "GO:0008150": 0,
            "GO:0009987": 1,
            "GO:0006281": 5,
        }
        with patch.object(enr, "_get_go_levels", return_value=mocked_levels):
            # Keep only levels 1-5 → exclude the root (level 0)
            result = enr.go_filter(df=df, min_level=1, max_level=5)
            assert len(result) == 2
            assert "GO:0008150" not in result["Term"].values
            assert any("GO:0009987" in t for t in result["Term"].values)
            assert any("GO:0006281" in t for t in result["Term"].values)

        enr.close()

    def test_gofilter_function(self):
        """gseapy.gofilter convenience function should delegate to Enrichr.go_filter."""
        from unittest.mock import patch

        import pandas as pd

        import gseapy

        df = pd.DataFrame(
            {
                "Gene_set": ["GO_BP"] * 2,
                "Term": [
                    "biological_process (GO:0008150)",
                    "DNA repair (GO:0006281)",
                ],
                "P-value": [0.001, 0.003],
                "Adjusted P-value": [0.01, 0.03],
                "Genes": ["G1", "G2"],
            }
        )
        mocked_levels = {"GO:0008150": 0, "GO:0006281": 5}

        with patch(
            "gseapy.enrichr.Enrichr._get_go_levels", return_value=mocked_levels
        ):
            result = gseapy.gofilter(df, min_level=1, max_level=20)
            # Root term (level 0) should be excluded
            assert len(result) == 1
            assert "GO:0006281" in result["Term"].values[0]


class TestEnrichOnline:
    """Tests for Enrichr online API calls (these require network access)."""

    @pytest.mark.network
    def test_enrichr_kegg(self, genelist):
        """enrichr() with KEGG_2016 library should return results."""
        result = enrichr(
            gene_list=genelist,
            gene_sets="KEGG_2016",
            organism="human",
            outdir=None,
            cutoff=0.5,
            no_plot=True,
        )
        assert result.res2d is not None
        assert not result.res2d.empty

    @pytest.mark.network
    def test_enrichr_multiple_libraries(self, genelist):
        """enrichr() should handle multiple Enrichr library names."""
        result = enrichr(
            gene_list=genelist,
            gene_sets=["KEGG_2016", "GO_Biological_Process_2021"],
            organism="human",
            outdir=None,
            cutoff=0.5,
            no_plot=True,
        )
        assert result.res2d is not None
        assert "Gene_set" in result.res2d.columns

    @pytest.mark.network
    def test_go_filter_online(self, genelist):
        """go_filter should return a non-empty subset of GO results."""
        result = enrichr(
            gene_list=genelist,
            gene_sets="GO_Biological_Process_2021",
            organism="human",
            outdir=None,
            cutoff=1.0,  # keep all terms
            no_plot=True,
        )
        assert result.res2d is not None
        all_results = result.res2d
        filtered = result.go_filter(min_level=3, max_level=10)
        # The filtered set should be a subset of the original
        assert len(filtered) <= len(all_results)
        # Filtered result should only contain terms with GO IDs at the right levels
        assert isinstance(filtered, pd.DataFrame)


def test_replot(edbDIR):
    # Only tests of the command runs successfully,
    # doesnt't check the image
    tmpdir = TemporaryDirectory(dir="tests")
    replot(edbDIR, tmpdir.name)
    tmpdir.cleanup()


# ---------------------------------------------------------------------------
# CLI argument parsing tests
# ---------------------------------------------------------------------------

from gseapy.__main__ import prepare_argparser


class TestCLIArgParsing:
    """Tests that CLI argument parsing produces correct argument values."""

    def test_gsea_required_args(self):
        parser = prepare_argparser()
        args = parser.parse_args([
            "gsea",
            "-d", "expr.txt",
            "-c", "test.cls",
            "-g", "gene_sets.gmt",
        ])
        assert args.subcommand_name == "gsea"
        assert args.data == "expr.txt"
        assert args.cls == "test.cls"
        assert args.gmt == "gene_sets.gmt"

    def test_gsea_defaults(self):
        parser = prepare_argparser()
        args = parser.parse_args([
            "gsea", "-d", "expr.txt", "-c", "test.cls", "-g", "sets.gmt",
        ])
        assert args.organism == "human"
        assert args.type == "phenotype"
        assert args.n == 1000
        assert args.mins == 15
        assert args.maxs == 500
        assert args.weight == 1.0
        assert args.method == "signal_to_noise"
        assert args.ascending is False
        assert args.seed == 123
        assert args.threads == 4
        assert args.noplot is False
        assert args.verbose is False
        assert args.format == "pdf"

    def test_gsea_custom_args(self):
        parser = prepare_argparser()
        args = parser.parse_args([
            "gsea",
            "-d", "expr.txt",
            "-c", "test.cls",
            "-g", "sets.gmt",
            "--org", "mouse",
            "-t", "gene_set",
            "-n", "500",
            "--min-size", "5",
            "--max-size", "1000",
            "-w", "2.0",
            "-m", "t_test",
            "-a",
            "-s", "42",
            "-p", "8",
            "--no-plot",
            "-f", "png",
        ])
        assert args.organism == "mouse"
        assert args.type == "gene_set"
        assert args.n == 500
        assert args.mins == 5
        assert args.maxs == 1000
        assert args.weight == 2.0
        assert args.method == "t_test"
        assert args.ascending is True
        assert args.seed == 42
        assert args.threads == 8
        assert args.noplot is True
        assert args.format == "png"

    def test_gsea_missing_required(self):
        parser = prepare_argparser()
        with pytest.raises(SystemExit):
            parser.parse_args(["gsea", "-d", "expr.txt"])

    def test_prerank_required_args(self):
        parser = prepare_argparser()
        args = parser.parse_args([
            "prerank", "-r", "ranked.rnk", "-g", "sets.gmt",
        ])
        assert args.subcommand_name == "prerank"
        assert args.rnk == "ranked.rnk"
        assert args.gmt == "sets.gmt"

    def test_prerank_defaults(self):
        parser = prepare_argparser()
        args = parser.parse_args([
            "prerank", "-r", "ranked.rnk", "-g", "sets.gmt",
        ])
        assert args.organism == "human"
        assert args.label == ("Pos", "Neg")
        assert args.n == 1000
        assert args.mins == 15
        assert args.maxs == 500
        assert args.weight == 1.0
        assert args.ascending is False
        assert args.seed == 123

    def test_prerank_custom_label(self):
        parser = prepare_argparser()
        args = parser.parse_args([
            "prerank", "-r", "ranked.rnk", "-g", "sets.gmt",
            "-l", "Up", "Down",
        ])
        assert args.label == ["Up", "Down"]

    def test_ssgsea_required_args(self):
        parser = prepare_argparser()
        args = parser.parse_args([
            "ssgsea", "-d", "expr.txt", "-g", "sets.gmt",
        ])
        assert args.subcommand_name == "ssgsea"
        assert args.data == "expr.txt"
        assert args.gmt == "sets.gmt"

    def test_ssgsea_defaults(self):
        parser = prepare_argparser()
        args = parser.parse_args([
            "ssgsea", "-d", "expr.txt", "-g", "sets.gmt",
        ])
        assert args.organism == "human"
        assert args.norm == "rank"
        assert args.correl == "rank"
        assert args.n == 0
        assert args.weight == 0.25
        assert args.mins == 15
        assert args.maxs == 2000

    def test_ssgsea_custom_args(self):
        parser = prepare_argparser()
        args = parser.parse_args([
            "ssgsea", "-d", "expr.txt", "-g", "sets.gmt",
            "--sn", "log", "-c", "zscore", "-n", "100",
        ])
        assert args.norm == "log"
        assert args.correl == "zscore"
        assert args.n == 100

    def test_gsva_required_args(self):
        parser = prepare_argparser()
        args = parser.parse_args([
            "gsva", "-d", "expr.txt", "-g", "sets.gmt",
        ])
        assert args.subcommand_name == "gsva"
        assert args.data == "expr.txt"
        assert args.gmt == "sets.gmt"

    def test_gsva_defaults(self):
        parser = prepare_argparser()
        args = parser.parse_args([
            "gsva", "-d", "expr.txt", "-g", "sets.gmt",
        ])
        assert args.organism == "human"
        assert args.kcdf == "Gaussian"
        assert args.weight == 1.0
        assert args.mx_diff is True
        assert args.abs_rnk is False
        assert args.mins == 15
        assert args.maxs == 2000

    def test_gsva_custom_args(self):
        parser = prepare_argparser()
        args = parser.parse_args([
            "gsva", "-d", "expr.txt", "-g", "sets.gmt",
            "-k", "Poisson", "-m", "-a", "-w", "0.5",
        ])
        assert args.kcdf == "Poisson"
        assert args.mx_diff is False
        assert args.abs_rnk is True
        assert args.weight == 0.5

    def test_enrichr_required_args(self):
        parser = prepare_argparser()
        args = parser.parse_args([
            "enrichr", "-i", "genes.txt", "-g", "KEGG_2016",
        ])
        assert args.subcommand_name == "enrichr"
        assert args.gene_list == "genes.txt"
        assert args.library == "KEGG_2016"

    def test_replot_required_args(self):
        parser = prepare_argparser()
        args = parser.parse_args([
            "replot", "-i", "gsea_results_dir",
        ])
        assert args.subcommand_name == "replot"
        assert args.indir == "gsea_results_dir"

    def test_no_subcommand(self):
        parser = prepare_argparser()
        args = parser.parse_args([])
        assert args.subcommand_name is None
