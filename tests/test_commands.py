import pytest
from tempfile import NamedTemporaryFile, TemporaryDirectory, mkdtemp
from gseapy.gsea import gsea, prerank, ssgsea, replot
from gseapy.enrichr import enrichr

@pytest.fixture
def edbDIR():
    return "tests/data"

@pytest.fixture
def genelist():
    return "tests/data/gene_list.txt"

@pytest.fixture
def gseaGCT():
    return "tests/data/P53_resampling_data.txt"

@pytest.fixture
def gseaCLS():
    return "tests/data/P53.cls"

@pytest.fixture
def prernk():
    return "tests/data/edb/gsea_data.gsea_data.rnk"

@pytest.fixture
def geneGMT():
    return "tests/data/genes.gmt"

@pytest.fixture
def ssGCT():
    return "tests/data/ss/ss_test.txt"

def test_gsea(gseaGCT, gseaCLS, geneGMT):
    # Only tests of the command runs successfully,
    # doesnt't check the image
    tmpdir= TemporaryDirectory(dir="tests")
    gsea(gseaGCT, geneGMT, gseaCLS, tmpdir.name)
    tmpdir.cleanup()


def test_prerank(prernk, geneGMT):
    # Only tests of the command runs successfully,
    # doesnt't check the image
    tmpdir= TemporaryDirectory(dir="tests")
    prerank(prernk, geneGMT, tmpdir.name)
    tmpdir.cleanup()

def test_ssgsea(ssGCT, geneGMT):
    # Only tests of the command runs successfully,
    # doesnt't check the image
    tmpdir= TemporaryDirectory(dir="tests")
    ssgsea(ssGCT, geneGMT, tmpdir.name)
    tmpdir.cleanup()

def test_enrichr(genelist):
    # Only tests of the command runs successfully,
    # doesnt't check the image
    tmpdir= TemporaryDirectory(dir="tests")
    enrichr(genelist, gene_sets='KEGG_2016', outdir=tmpdir.name)
    tmpdir.cleanup()


def test_replot(edbDIR):
    # Only tests of the command runs successfully,
    # doesnt't check the image
    tmpdir= TemporaryDirectory(dir="tests")
    replot(edbDIR, tmpdir.name)
    tmpdir.cleanup()
