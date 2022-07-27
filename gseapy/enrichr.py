#!/usr/bin/env python
# -*- coding: utf-8 -*-
# see: http://amp.pharm.mssm.edu/Enrichr/help#api for API docs

import json
import logging
import os
from collections import OrderedDict
from io import StringIO
from tempfile import TemporaryDirectory
from time import sleep
from typing import Any, AnyStr, Dict, Iterable, List, Optional, Set, Tuple, Union

import pandas as pd
import requests
from numpy import isscalar

from gseapy.biomart import Biomart
from gseapy.plot import barplot
from gseapy.stats import calc_pvalues, multiple_testing_correction
from gseapy.utils import *


class Enrichr(object):
    """Enrichr API"""

    def __init__(
        self,
        gene_list: Iterable[str],
        gene_sets: Union[List[str], str, Dict[str, str]],
        organism: str = "human",
        outdir: Optional[str] = "Enrichr",
        background: Union[List[str], int, str] = "hsapiens_gene_ensembl",
        cutoff: float = 0.05,
        format: str = "pdf",
        figsize: Tuple[float, float] = (6.5, 6),
        top_term: int = 10,
        no_plot: bool = False,
        verbose: bool = False,
    ):

        self.gene_list = gene_list
        self.gene_sets = gene_sets
        self.descriptions = ""
        self.outdir = outdir
        self.cutoff = cutoff
        self.format = format
        self.figsize = figsize
        self.__top_term = int(top_term)
        self.__no_plot = no_plot
        self.verbose = bool(verbose)
        self.module = "enrichr"
        self.res2d = None
        self.background = background
        self._bg = None
        self.organism = organism
        self._organism = None
        self.ENRICHR_URL = "http://maayanlab.cloud"
        # init logger
        logfile = self.prepare_outdir()
        self._logger = log_init(
            outlog=logfile, log_level=logging.INFO if self.verbose else logging.WARNING
        )

    def prepare_outdir(self):
        """create temp directory."""
        self._outdir = self.outdir
        if self._outdir is None:
            self._tmpdir = TemporaryDirectory()
            self.outdir = self._tmpdir.name
        elif isinstance(self.outdir, str):
            mkdirs(self.outdir)
        else:
            raise Exception("Error parsing outdir: %s" % type(self.outdir))

        # handle gene_sets
        logfile = os.path.join(
            self.outdir, "gseapy.%s.%s.log" % (self.module, self.descriptions)
        )
        return logfile

    def __parse_gmt(self, g: str):
        with open(g) as genesets:
            g_dict = {
                line.strip().split("\t")[0]: line.strip().split("\t")[2:]
                for line in genesets.readlines()
            }
        return g_dict

    def __gmt2dict(self, gene_sets: List[str]) -> List[Dict[str, List[str]]]:
        """helper function, only convert gmt to dict and keep strings"""
        gss = []
        for g in gene_sets:
            # only convert gmt to dict. local mode
            if isinstance(g, str) and g.lower().endswith(".gmt"):
                if os.path.exists(g):
                    self._logger.info("User Defined gene sets is given: %s" % g)
                    gss.append(self.__parse_gmt(g))
                else:
                    self._logger.warning("User Defined gene sets is not found: %s" % g)
            else:
                gss.append(g)

        return gss

    def parse_genesets(self):
        """parse gene_sets input file type"""

        gss = []
        if isinstance(self.gene_sets, list):
            gss = self.__gmt2dict(self.gene_sets)

        elif isinstance(self.gene_sets, str):
            gss = [g.strip() for g in self.gene_sets.strip().split(",")]
            gss = self.__gmt2dict(gss)

        elif isinstance(self.gene_sets, dict):
            gss = [self.gene_sets]
        else:
            raise Exception(
                "Error parsing enrichr libraries, please provided corrected one"
            )

        # now, gss[LIST] contains dict or strings.
        if len(gss) < 1:
            raise Exception("No GeneSets are valid !!! Check your gene_sets input.")
        gss_exist = []
        enrichr_library = []
        # if all local gmts (local mode), skip connect to enrichr server
        if not all([isinstance(g, dict) for g in gss]):
            enrichr_library = self.get_libraries()

        # check enrichr libraries are valid

        for g in gss:
            if isinstance(g, dict):
                gss_exist.append(g)
                continue
            if isinstance(g, str):
                if g in enrichr_library:
                    gss_exist.append(g)
                else:
                    self._logger.warning("Enrichr library not found: %s" % g)

        return gss_exist

    def parse_genelists(self) -> str:
        """parse gene list"""
        if isinstance(self.gene_list, list):
            genes = self.gene_list
        elif isinstance(self.gene_list, pd.DataFrame):
            # input type is bed file
            if self.gene_list.shape[1] >= 3:
                genes = (
                    self.gene_list.iloc[:, :3]
                    .apply(lambda x: "\t".join([str(i) for i in x]), axis=1)
                    .tolist()
                )
            # input type with weight values
            elif self.gene_list.shape[1] == 2:
                genes = self.gene_list.apply(
                    lambda x: ",".join([str(i) for i in x]), axis=1
                ).tolist()
            else:
                genes = self.gene_list.squeeze().tolist()
        elif isinstance(self.gene_list, pd.Series):
            genes = self.gene_list.squeeze().tolist()
        else:
            # get gene lists or bed file, or gene list with weighted values.
            genes = []
            with open(self.gene_list) as f:
                for gene in f:
                    genes.append(gene.strip())

        self._isezid = all(map(self._is_entrez_id, genes))
        if self._isezid:
            self._gls = set(map(int, genes))
        else:
            self._gls = genes

        return "\n".join(genes)

    def send_genes(self, gene_list, url) -> AnyStr:
        """send gene list to enrichr server"""
        payload = {"list": (None, gene_list), "description": (None, self.descriptions)}
        # response
        response = requests.post(url, files=payload, verify=True)
        if not response.ok:
            raise Exception("Error analyzing gene list")
        sleep(1)
        job_id = json.loads(response.text)

        return job_id

    def check_genes(self, gene_list: List[str], usr_list_id: str):
        """
        Compare the genes sent and received to get successfully recognized genes
        """
        response = requests.get(
            "%s/%s/view?userListId=%s"
            % (self.ENRICHR_URL, self._organism, usr_list_id),
            verify=True,
        )
        if not response.ok:
            raise Exception("Error getting gene list back")
        returnedL = json.loads(response.text)["genes"]
        returnedN = sum([1 for gene in gene_list if gene in returnedL])
        self._logger.info(
            "{} genes successfully recognized by Enrichr".format(returnedN)
        )

    def get_results(self, gene_list: List[str]) -> Tuple[AnyStr, pd.DataFrame]:
        """Enrichr API"""
        ADDLIST_URL = "%s/%s/addList" % (self.ENRICHR_URL, self._organism)
        job_id = self.send_genes(gene_list, ADDLIST_URL)
        user_list_id = job_id["userListId"]

        RESULTS_URL = "%s/%s/export" % (self.ENRICHR_URL, self._organism)
        query_string = "?userListId=%s&filename=%s&backgroundType=%s"
        # set max retries num =5
        s = retry(num=5)
        filename = "%s.%s.reports" % (self._gs, self.descriptions)
        url = RESULTS_URL + query_string % (user_list_id, filename, self._gs)
        response = s.get(url, stream=True)
        response.encoding = "utf-8"
        if not response.ok:
            self._logger.error("Error fetching enrichment results: %s" % self._gs)

        try:
            res = pd.read_csv(StringIO(response.text), sep="\t")
        except pd.errors.ParserError as e:
            RESULTS_URL = "%s/Enrichr/enrich" % self.ENRICHR_URL
            query_string = "?userListId=%s&backgroundType=%s"
            url = RESULTS_URL + query_string % (user_list_id, self._gs)
            response = s.get(url)
            if not response.ok:
                self._logger.error("Error fetching enrichment results: %s" % self._gs)
            data = json.loads(response.text)
            colnames = [
                "Rank",
                "Term",
                "P-value",
                "Z-score",
                "Combined Score",
                "Genes",
                "Adjusted P-value",
                "Old P-value",
                "Old adjusted P-value",
            ]
            res = pd.DataFrame(data[self._gs], columns=colnames)
            res["Genes"] = res["Genes"].apply(";".join)

        return [job_id["shortId"], res]

    def _is_entrez_id(self, idx: Union[int, str]) -> bool:
        try:
            int(idx)
            return True
        except:
            return False

    def get_libraries(self) -> List[str]:
        """return active enrichr library name. Official API"""

        lib_url = "%s/%s/datasetStatistics" % (self.ENRICHR_URL, self._organism)
        response = requests.get(lib_url, verify=True)
        if not response.ok:
            raise Exception("Error getting the Enrichr libraries")
        libs_json = json.loads(response.text)
        libs = [lib["libraryName"] for lib in libs_json["statistics"]]

        return sorted(libs)

    def get_background(self) -> Set[str]:
        """get background gene"""

        # input is a file
        if os.path.isfile(self.background):
            with open(self.background) as b:
                bg2 = b.readlines()
            bg = [g.strip() for g in bg2]
            return set(bg)

        # package included data
        DB_FILE = os.path.join(
            DEFAULT_CACHE_PATH, "{}.background.genes.txt".format(self.background)
        )
        if os.path.exists(DB_FILE):
            df = pd.read_csv(DB_FILE, sep="\t")
        else:
            # background is a biomart database name
            self._logger.warning(
                "Downloading %s for the first time. It might take a couple of miniutes."
                % self.background
            )
            bm = Biomart()
            df = bm.query(
                dataset=self.background,
                attributes=["ensembl_gene_id", "external_gene_name", "entrezgene_id"],
                filename=os.path.join(
                    DEFAULT_CACHE_PATH,
                    "{}.background.genes.txt".format(self.background),
                ),
            )
        self._logger.info(
            "Using all annotated genes with GO_ID as background: %s" % self.background
        )
        df.dropna(subset=["entrezgene_id"], inplace=True)
        # input id type: entrez or gene_name
        if self._isezid:
            bg = df["entrezgene_id"].astype(int)
        else:
            bg = df["external_gene_name"]

        return set(bg)

    def set_organism(self):
        """Select Enrichr organism from below:

        Human & Mouse, H. sapiens & M. musculus
        Fly, D. melanogaster
        Yeast, S. cerevisiae
        Worm, C. elegans
        Fish, D. rerio

        """

        default = [
            "human",
            "mouse",
            "hs",
            "mm",
            "homo sapiens",
            "mus musculus",
            "h. sapiens",
            "m. musculus",
        ]

        if self.organism.lower() in default:
            self._organism = "Enrichr"
            return

        organism = {
            "Fly": ["fly", "d. melanogaster", "drosophila melanogaster"],
            "Yeast": ["yeast", "s. cerevisiae", "saccharomyces cerevisiae"],
            "Worm": ["worm", "c. elegans", "caenorhabditis elegans", "nematode"],
            "Fish": ["fish", "d. rerio", "danio rerio", "zebrafish"],
        }

        for k, v in organism.items():
            if self.organism.lower() in v:
                self._organism = k + "Enrichr"
                return

        if self._organism is None:
            raise Exception("No supported organism found !!!")

        ENRICHR_SERVER = "%s/%s" % (self.ENRICHR_URL, self._organism)

        if requests.get(ENRICHR_SERVER, verify=True).ok:
            return

        # self.ENRICHR_URL = 'http://amp.pharm.mssm.edu'
        ENRICHR_SERVER = "%s/%s" % (self.ENRICHR_URL, self._organism)

        if requests.get(ENRICHR_SERVER, verify=True).ok:
            return
        else:
            raise Exception("Please check Enrichr URL is OK: %s" % self.ENRICHR_URL)

        return

    def filter_gmt(self, gmt, background):
        """the gmt values should be filtered only for genes that exist in background
        this substantially affect the significance of the test, the hypergeometric distribution.

        :param gmt: a dict of gene sets.
        :param background: list, set, or tuple. A list of custom backgound genes.
        """
        gmt2 = {}
        for term, genes in gmt.items():
            # If value satisfies the condition, then store it in new_dict
            newgenes = [g for g in genes if g in background]
            if len(newgenes) > 0:
                gmt2[term] = newgenes
        return gmt2

    def enrich(self, gmt: Dict[str, List[str]]):
        """use local mode

        p = p-value computed using the Fisher exact test (Hypergeometric test)

        Not implemented here:

            combine score = log(p)·z

        see here: http://amp.pharm.mssm.edu/Enrichr/help#background&q=4

        columns contain:

            Term Overlap P-value Adjusted_P-value Genes

        """
        if isscalar(self.background):
            if isinstance(self.background, int) or self.background.isdigit():
                self._bg = int(self.background)
            elif isinstance(self.background, str):
                # self.background = set(reduce(lambda x,y: x+y, gmt.values(),[]))
                self._bg = self.get_background()
                self._logger.info("Background: found %s genes" % (len(self._bg)))
            else:
                raise Exception("Unsupported background data type")
        else:
            # handle array object: nd.array, list, tuple, set, Series
            try:
                it = iter(self.background)
                self._bg = set(self.background)
            except TypeError:
                self._logger.error("Unsupported background data type")
        # statistical testing
        hgtest = list(calc_pvalues(query=self._gls, gene_sets=gmt, background=self._bg))
        if len(hgtest) > 0:
            terms, pvals, oddr, olsz, gsetsz, genes = hgtest
            fdrs, rej = multiple_testing_correction(
                ps=pvals, alpha=self.cutoff, method="benjamini-hochberg"
            )
            # save to a dataframe
            odict = OrderedDict()
            odict["Term"] = terms
            odict["Overlap"] = list(map(lambda h, g: "%s/%s" % (h, g), olsz, gsetsz))
            odict["P-value"] = pvals
            odict["Adjusted P-value"] = fdrs
            odict["Odds Ratio"] = oddr
            # odict['Reject (FDR< %s)'%self.cutoff ] = rej
            odict["Genes"] = [";".join(g) for g in genes]
            res = pd.DataFrame(odict)
            return res
        return

    def run(self):
        """run enrichr for one sample gene list but multi-libraries"""

        # read input file
        genes_list = self.parse_genelists()

        self._logger.info("Connecting to Enrichr Server to get latest library names")
        gss = self.parse_genesets()
        if len(gss) < 1:
            self._logger.error(
                "None of your input gene set matched ! %s" % self.gene_sets
            )
            self._logger.error(
                "Hint: Current organism = %s, is this correct?\n" % self.organism
                + "Hint: use get_library_name() to view full list of supported names."
            )
            raise LookupError(
                "Not validated Enrichr library ! Please provide correct organism and library name!"
            )
        self.results = []

        for g in gss:
            if isinstance(g, dict):
                ## local mode
                res = self.enrich(g)
                shortID, self._gs = str(id(g)), "CUSTOM%s" % id(g)
                if res is None:
                    self._logger.info(
                        "No hits return, for gene set: Custom%s" % shortID
                    )
                    continue
            else:
                ## online mode
                self._gs = str(g)
                self._logger.debug("Start Enrichr using library: %s" % (self._gs))
                self._logger.info(
                    "Analysis name: %s, Enrichr Library: %s"
                    % (self.descriptions, self._gs)
                )
                shortID, res = self.get_results(genes_list)
                # Remember gene set library used
            res.insert(0, "Gene_set", self._gs)
            # Append to master dataframe
            self.results.append(res)
            self.res2d = res
            if self._outdir is None:
                continue
            self._logger.info("Save file of enrichment results: Job Id:" + str(shortID))
            outfile = "%s/%s.%s.%s.reports.txt" % (
                self.outdir,
                self._gs,
                self.organism,
                self.module,
            )
            self.res2d.to_csv(
                outfile, index=False, encoding="utf-8", float_format="%.6e", sep="\t"
            )
            # plotting
            if not self.__no_plot:
                msg = barplot(
                    df=res,
                    cutoff=self.cutoff,
                    figsize=self.figsize,
                    top_term=self.__top_term,
                    color="salmon",
                    title=self._gs,
                    ofname=outfile.replace("txt", self.format),
                )
                if msg is not None:
                    self._logger.warning(msg)
            self._logger.info("Done.\n")
        self.results = pd.concat(self.results, ignore_index=True)
        # clean up tmpdir
        if self._outdir is None:
            self._tmpdir.cleanup()

        return
