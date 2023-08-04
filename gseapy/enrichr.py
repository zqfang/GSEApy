#!/usr/bin/env python
# -*- coding: utf-8 -*-
# see: http://amp.pharm.mssm.edu/Enrichr/help#api for API docs

import json
import logging
import os
from collections import OrderedDict
from io import StringIO
from typing import Any, AnyStr, Dict, Iterable, List, Optional, Set, Tuple, Union

import pandas as pd
import requests
from numpy import isscalar, log

from gseapy.biomart import Biomart
from gseapy.plot import barplot
from gseapy.stats import calc_pvalues, multiple_testing_correction
from gseapy.utils import DEFAULT_CACHE_PATH, log_init, mkdirs, retry


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
        self.ENRICHR_URL_SPEED = "https://maayanlab.cloud/speedrichr"
        # init logger
        self.prepare_outdir()

    def __del__(self):
        handlers = self._logger.handlers[:]
        for handler in handlers:
            handler.close()  # close file
            self._logger.removeHandler(handler)

    def prepare_outdir(self):
        """create temp directory."""
        self._outdir = self.outdir

        logfile = None
        if isinstance(self.outdir, str):
            mkdirs(self.outdir)
            logfile = os.path.join(
                self.outdir, "gseapy.%s.%s.log" % (self.module, id(self))
            )
        self._logger = log_init(
            name=str(self.module) + str(id(self)),
            log_level=logging.INFO if self.verbose else logging.WARNING,
            filename=logfile,
        )

    def __parse_gmt(self, g: str):
        with open(g) as genesets:
            g_dict = {
                line.strip().split("\t")[0]: line.strip().split("\t")[2:]
                for line in genesets.readlines()
            }
        return g_dict

    def __gs2dict(self, gene_sets: List[str]) -> List[Dict[str, List[str]]]:
        """helper function, only convert gmt to dict and keep strings"""
        gss = []
        self._gs_name = []
        for i, g in enumerate(gene_sets):
            # only convert gmt to dict. local mode
            if isinstance(g, str) and g.lower().endswith(".gmt"):
                if os.path.exists(g):
                    self._logger.info("User defined gene sets is given: %s" % g)
                    gss.append(self.__parse_gmt(g))
                    self._gs_name.append(os.path.basename(g))
                else:
                    self._logger.warning(
                        "User defined gene sets is not found: %s, skip." % g
                    )
            else:
                gss.append(g)
                _name = g
                if isinstance(g, dict):
                    _name = "gs_ind_" + str(i)
                    self._logger.info("Input dict object named with %s" % _name)
                self._gs_name.append(_name)
        return gss

    def parse_genesets(self, gene_sets=None):
        """parse gene_sets input file type"""
        if gene_sets is None:
            gene_sets = self.gene_sets

        gss = []
        if isinstance(gene_sets, list):
            gss = self.__gs2dict(gene_sets)
        elif isinstance(self.gene_sets, str):
            gss = [g.strip() for g in gene_sets.strip().split(",")]
            gss = self.__gs2dict(gss)

        elif isinstance(gene_sets, dict):
            gss = self.__gs2dict([gene_sets.copy()])
        else:
            raise Exception(
                "Error parsing enrichr libraries, please provided corrected one"
            )

        # now, gss[LIST] contains dict or strings.
        if len(gss) < 1:
            raise Exception("No GeneSets are valid !!! Check your gene_sets input.")
        gss_exist = []
        gss_name = []
        enrichr_library = []
        # if all local gmts (local mode), skip connect to enrichr server
        if not all([isinstance(g, dict) for g in gss]):
            enrichr_library = self.get_libraries()

        # check enrichr libraries are valid
        for n, g in zip(self._gs_name, gss):
            if isinstance(g, dict):
                gss_exist.append(g)
                gss_name.append(n)
                continue
            if isinstance(g, str):
                if g in enrichr_library:
                    gss_exist.append(g)
                    gss_name.append(n)
                else:
                    self._logger.warning("Input library not found: %s. Skip" % g)

        self._gs_name = gss_name  # update names
        if len(gss_exist) < 1:
            raise Exception("No GeneSets are valid !!! Check your gene_sets input.")
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

    def send_genes(self, payload, url) -> Dict:
        """send gene list to enrichr server"""
        # payload = {"list": (None, gene_list), "description": (None, self.descriptions)}
        # response
        s = retry(num=5)
        response = s.post(url, files=payload, verify=True)
        if not response.ok:
            self._logger.debug(url)
            self._logger.debug(payload)
            raise Exception("Error sending gene list, try again later")
        job_id = json.loads(response.text)
        # response.text
        # {
        # 	"userListId": 667152768,
        # 	"shortId": "27c3f180"
        # }
        return job_id

    def send_background(self, payload, url) -> Dict:
        s = retry(num=5)
        res = s.post(url, data=payload)
        if not res.ok:
            self._logger.debug(url)
            self._logger.debug(payload)
            raise Exception("Error sending background list, try again later")
        background_response = res.json()
        # {
        #     "backgroundid": "3ff7ef9d"
        # }
        return background_response

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

    def get_results_with_background(
        self, gene_list: List[str], background: List[str]
    ) -> Tuple[AnyStr, pd.DataFrame]:
        ## add gene list
        ADDLIST_URL = "%s/api/addList" % (self.ENRICHR_URL_SPEED)
        payload = {"list": (None, gene_list), "description": (None, self.descriptions)}
        job_id = self.send_genes(payload, ADDLIST_URL)
        ## add background list
        ADDBG_URL = "%s/api/addbackground" % (self.ENRICHR_URL_SPEED)
        payload = dict(background="\n".join(background))
        bg_id = self.send_background(payload, ADDBG_URL)

        # now get background enrich result
        BGENR_URL = "%s/api/backgroundenrich" % (self.ENRICHR_URL_SPEED)
        payload = dict(
            userListId=job_id["userListId"],
            backgroundid=bg_id["backgroundid"],
            backgroundType=self._gs,
        )
        s = retry(num=5)
        response = s.post(BGENR_URL, data=payload)
        if not response.ok:
            self._logger.error("Error fetching enrichment results: %s" % self._gs)

        data = response.json()
        # Note: missig Overlap column
        colnames = [
            "Rank",
            "Term",
            "P-value",
            "Odds Ratio",  # Z-Score
            "Combined Score",
            "Genes",
            "Adjusted P-value",
            "Old P-value",
            "Old adjusted P-value",
        ]
        res = pd.DataFrame(data[self._gs], columns=colnames)
        # res.drop(columns=["Rank"], inplace=True)
        res["Genes"] = res["Genes"].apply(";".join)
        colord = [
            "Term",
            "P-value",
            "Adjusted P-value",
            "Old P-value",
            "Old adjusted P-value",
            "Odds Ratio",  # Z-Score
            "Combined Score",
            "Genes",
        ]
        res = res.loc[:, colord]
        return (job_id["shortId"], res)

    def get_results(self, gene_list: List[str]) -> Tuple[AnyStr, pd.DataFrame]:
        """Enrichr API"""
        ADDLIST_URL = "%s/%s/addList" % (self.ENRICHR_URL, self._organism)
        payload = {"list": (None, gene_list), "description": (None, self.descriptions)}
        job_id = self.send_genes(payload, ADDLIST_URL)
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
                "Odds Ratio",  # 'oddsratio'
                "Combined Score",
                "Genes",
                "Adjusted P-value",
                "Old P-value",
                "Old adjusted P-value",
            ]
            res = pd.DataFrame(data[self._gs], columns=colnames)
            # res.drop(columns=["Rank"], inplace=True)
            res["Genes"] = res["Genes"].apply(";".join)
            colord = [
                "Term",
                "P-value",
                "Adjusted P-value",
                "Old P-value",
                "Old adjusted P-value",
                "Odds Ratio",  # Z-Score
                "Combined Score",
                "Genes",
            ]
            res = res.loc[:, colord]

        return (job_id["shortId"], res)

    def _is_entrez_id(self, idx: Union[int, str]) -> bool:
        try:
            int(idx)
            return True
        except:
            return False

    def get_libraries(self) -> List[str]:
        """return active enrichr library name. Official API"""

        lib_url = "%s/%s/datasetStatistics" % (self.ENRICHR_URL, self._organism)
        s = retry(num=5)
        response = s.get(lib_url, verify=True)
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
        # input id type: entrez or gene_name
        if self._isezid:
            df.dropna(subset=["entrezgene_id"], inplace=True)
            bg = df["entrezgene_id"].astype(int)
        else:
            df.dropna(subset=["external_gene_name"], inplace=True)
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

    def parse_background(self, gmt: Dict[str, List[str]] = None):
        """
        set background genes
        """
        if hasattr(self, "_bg") and self._bg:
            return self._bg

        self._bg = set()
        if self.background is None:
            # use all genes in the dict input as background if background is None
            if gmt:
                bg = set()
                for term, genes in gmt.items():
                    bg = bg.union(set(genes))
                self._logger.info(
                    "Background is not set! Use all %s genes in %s."
                    % (len(bg), self._gs)
                )
                self._bg = bg
        elif isscalar(self.background):
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

        return self._bg

    def enrich(self, gmt: Dict[str, List[str]]):
        """use local mode

        p = p-value computed using the Fisher exact test (Hypergeometric test)
        z = z-score (Odds Ratio)
        combine score = - log(p)·z

        see here: http://amp.pharm.mssm.edu/Enrichr/help#background&q=4

        columns contain:

            Term Overlap P-value Odds Ratio Combinde Score Adjusted_P-value Genes

        """
        bg = self.parse_background(gmt)
        # statistical testing
        hgtest = list(calc_pvalues(query=self._gls, gene_sets=gmt, background=bg))
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
            odict["Combined Score"] = -1 * log(pvals) * oddr
            # odict['Reject (FDR< %s)'%self.cutoff ] = rej
            odict["Genes"] = [";".join(map(str, g)) for g in genes]
            res = pd.DataFrame(odict)
            return res
        return

    def run(self):
        """run enrichr for one sample gene list but multi-libraries"""

        # read input file
        genes_list = self.parse_genelists()

        # self._logger.info("Connecting to Enrichr Server to get latest library names")
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

        for name, g in zip(self._gs_name, gss):
            self._logger.info("Run: %s " % name)
            if isinstance(g, dict):
                ## local mode
                shortID, self._gs = str(id(g)), name
                self._logger.debug(
                    "Off-line enrichment analysis with library: %s" % (self._gs)
                )
                if self._isezid:
                    g = {k: list(map(int, v)) for k, v in g.items()}
                res = self.enrich(g)
                if res is None:
                    self._logger.info(
                        "No hits return, for gene set: Custom%s" % shortID
                    )
                    continue
            else:
                ## online mode
                self._gs = name
                self._logger.debug("Enrichr service using library: %s" % (name))
                # self._logger.info("Enrichr Library: %s"% self._gs)
                bg = self.parse_background()
                # whether user input background
                if isinstance(bg, set) and len(bg) > 0:
                    shortID, res = self.get_results_with_background(
                        genes_list, self._bg
                    )
                else:
                    shortID, res = self.get_results(genes_list)

            # Remember gene set library used
            res.insert(0, "Gene_set", name)
            # Append to master dataframe
            self.results.append(res)
            self.res2d = res
            if self._outdir is None:
                continue
            outfile = "%s/%s.%s.%s.reports.txt" % (
                self.outdir,
                self._gs,
                self.organism,
                self.module,
            )
            self.res2d.to_csv(
                outfile, index=False, encoding="utf-8", float_format="%.6e", sep="\t"
            )
            self._logger.info("Save enrichment results for %s " % name)
            # plotting
            if not self.__no_plot:
                ax = barplot(
                    df=res,
                    cutoff=self.cutoff,
                    figsize=self.figsize,
                    top_term=self.__top_term,
                    color="salmon",
                    title=self._gs,
                    ofname=outfile.replace("txt", self.format),
                )
                self._logger.debug("Generate figures")
        self.results = pd.concat(self.results, ignore_index=True)
        self._logger.info("Done.")

        return
