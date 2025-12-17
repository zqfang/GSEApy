#!/usr/bin/env python
# -*- coding: utf-8 -*-
# see: http://amp.pharm.mssm.edu/Enrichr/help#api for API docs

import json
import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
from io import StringIO
from typing import Dict, Iterable, List, Optional, Set, Tuple, Union

import numpy as np
import pandas as pd
import requests

from gseapy.biomart import Biomart
from gseapy.plot import barplot
from gseapy.stats import calc_pvalues, multiple_testing_correction
from gseapy.utils import DEFAULT_CACHE_PATH, log_init, mkdirs, retry

# Module-level constants
ENRICHR_URL = "http://maayanlab.cloud"
ENRICHR_URL_SPEED = "https://maayanlab.cloud/speedrichr"

DEFAULT_ORGANISMS = [
    "human",
    "mouse",
    "hs",
    "mm",
    "homo sapiens",
    "mus musculus",
    "h. sapiens",
    "m. musculus",
]

ORGANISM_MAPPINGS = {
    "Fly": ["fly", "d. melanogaster", "drosophila melanogaster"],
    "Yeast": ["yeast", "s. cerevisiae", "saccharomyces cerevisiae"],
    "Worm": ["worm", "c. elegans", "caenorhabditis elegans", "nematode"],
    "Fish": ["fish", "d. rerio", "danio rerio", "zebrafish"],
}


class EnrichrError(Exception):
    """Base exception for Enrichr-related errors."""

    pass


class EnrichrValidationError(EnrichrError):
    """Raised when input validation fails."""

    pass


class EnrichrAPIError(EnrichrError):
    """Raised when Enrichr API requests fail."""

    pass


class EnrichrNetworkError(EnrichrError):
    """Raised when network connectivity issues occur."""

    pass


class EnrichrParseError(EnrichrError):
    """Raised when parsing responses or files fails."""

    pass


class Enrichr(object):
    """Enrichr API"""

    def __init__(
        self,
        gene_list: Iterable[str],
        gene_sets: Union[List[str], str, Dict[str, str]],
        organism: str = "human",
        outdir: Optional[str] = "Enrichr",
        background: Union[List[str], int, str] = None,
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
        self.res2d: Optional[pd.DataFrame] = None
        self.background = background
        self._bg: Union[Set[str], int, None] = None
        self.organism = organism
        self._organism: Optional[str] = None
        self._gene_isupper = True
        self._gene_toupper = False
        self._gls: Union[Set[int], List[str]] = []
        self._isezid = False
        self._gs: str = ""
        self._gs_name: List[str] = []
        # init logger
        self.prepare_outdir()

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()
        return False

    def close(self):
        """Explicitly close logger handlers."""
        if hasattr(self, "_logger"):
            handlers = self._logger.handlers[:]
            for handler in handlers:
                handler.close()
                self._logger.removeHandler(handler)

    def prepare_outdir(self):
        """create temp directory."""
        self._outdir = self.outdir

        logfile = None
        if isinstance(self.outdir, str):
            mkdirs(self.outdir)
            logfile = os.path.join(self.outdir, f"gseapy.{self.module}.{id(self)}.log")
        self._logger = log_init(
            name=f"{self.module}{id(self)}",
            log_level=logging.INFO if self.verbose else logging.WARNING,
            filename=logfile,
        )

    def _parse_gmt(self, g: str) -> Dict[str, List[str]]:
        """Parse GMT file efficiently with single split per line."""
        g_dict = {}
        with open(g) as genesets:
            for line in genesets:
                parts = line.strip().split("\t")
                if len(parts) > 2:
                    g_dict[parts[0]] = parts[2:]
        return g_dict

    def _gs2dict(self, gene_sets: List[str]) -> List[Dict[str, List[str]]]:
        """helper function, only convert gmt to dict and keep strings"""
        gss = []
        self._gs_name = []
        for i, g in enumerate(gene_sets):
            # only convert gmt to dict. local mode
            if isinstance(g, str) and g.lower().endswith(".gmt"):
                if os.path.exists(g):
                    self._logger.info(f"User defined gene sets is given: {g}")
                    gss.append(self._parse_gmt(g))
                    self._gs_name.append(os.path.basename(g))
                else:
                    self._logger.warning(
                        f"User defined gene sets is not found: {g}, skip."
                    )
            else:
                gss.append(g)
                _name = g
                if isinstance(g, dict):
                    _name = f"gs_ind_{i}"
                    self._logger.info(f"Input dict object named with {_name}")
                self._gs_name.append(_name)
        return gss

    def parse_genesets(self, gene_sets=None) -> List[Union[Dict[str, List[str]], str]]:
        """parse gene_sets input file type"""
        if gene_sets is None:
            gene_sets = self.gene_sets

        gss = []
        if isinstance(gene_sets, list):
            gss = self._gs2dict(gene_sets)
        elif isinstance(self.gene_sets, str):
            gss = [g.strip() for g in gene_sets.strip().split(",")]
            gss = self._gs2dict(gss)
        elif isinstance(gene_sets, dict):
            gss = self._gs2dict([gene_sets.copy()])
        else:
            raise EnrichrValidationError(
                "Error parsing enrichr libraries, please provided corrected one"
            )

        # now, gss[LIST] contains dict or strings.
        if len(gss) < 1:
            raise EnrichrValidationError(
                "No GeneSets are valid !!! Check your gene_sets input."
            )
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
                    self._logger.warning(f"Input library not found: {g}. Skip")

        self._gs_name = gss_name  # update names
        if len(gss_exist) < 1:
            raise EnrichrValidationError(
                "No GeneSets are valid !!! Check your gene_sets input."
            )
        return gss_exist

    def parse_genelists(self) -> str:
        """Parse gene list with single-pass processing."""
        if isinstance(self.gene_list, list):
            genes = self.gene_list
        elif isinstance(self.gene_list, pd.DataFrame):
            # input type is bed file
            if self.gene_list.shape[1] >= 3:
                genes = (
                    self.gene_list.iloc[:, :3]
                    .apply(lambda x: "\t".join([str(i) for i in x]), axis=1)
                    .to_list()
                )
            # input type with weight values
            elif self.gene_list.shape[1] == 2:
                genes = self.gene_list.apply(
                    lambda x: ",".join([str(i) for i in x]), axis=1
                ).to_list()
            else:
                genes = self.gene_list.squeeze().to_list()
        elif isinstance(self.gene_list, pd.Series):
            genes = self.gene_list.squeeze().to_list()
        else:
            # get gene lists or bed file, or gene list with weighted values.
            genes = []
            with open(self.gene_list) as f:
                for gene in f:
                    genes.append(gene.strip())

        if not genes:
            raise EnrichrValidationError("Gene list cannot be empty")

        # Single-pass processing: strip whitespace and check if entrez IDs
        cleaned_genes = []
        all_entrez = True
        for g in genes:
            cleaned = g.strip()
            cleaned_genes.append(cleaned)
            if all_entrez and not self._is_entrez_id(cleaned):
                all_entrez = False

        self._isezid = all_entrez
        if self._isezid:
            self._gls = set(map(int, cleaned_genes))
        else:
            self._gls = cleaned_genes

        return "\n".join(cleaned_genes)

    def send_genes(self, payload, url) -> Dict[str, Union[int, str]]:
        """Send gene list to enrichr server."""
        s = retry(num=5)
        try:
            response = s.post(url, files=payload, verify=True)
            if not response.ok:
                self._logger.debug(f"URL: {url}")
                self._logger.debug(f"Payload: {payload}")
                raise EnrichrAPIError(
                    f"Error sending gene list, status code: {response.status_code}"
                )
            job_id = json.loads(response.text)
            # response.text format:
            # {"userListId": 667152768, "shortId": "27c3f180"}
            return job_id
        except requests.exceptions.RequestException as e:
            raise EnrichrNetworkError(f"Network error sending gene list: {e}")

    def send_background(self, payload, url) -> Dict[str, str]:
        """Send background gene list to enrichr server."""
        s = retry(num=5)
        try:
            res = s.post(url, data=payload)
            if not res.ok:
                self._logger.debug(f"URL: {url}")
                self._logger.debug(f"Payload: {payload}")
                raise EnrichrAPIError(
                    f"Error sending background list, status code: {res.status_code}"
                )
            background_response = res.json()
            # response format: {"backgroundid": "3ff7ef9d"}
            return background_response
        except requests.exceptions.RequestException as e:
            raise EnrichrNetworkError(f"Network error sending background list: {e}")

    def check_genes(self, gene_list: List[str], usr_list_id: str) -> None:
        """Compare the genes sent and received to get successfully recognized genes."""
        s = retry(num=5)
        url = f"{ENRICHR_URL}/{self._organism}/view?userListId={usr_list_id}"
        try:
            response = s.get(url, verify=True)
            if not response.ok:
                raise EnrichrAPIError(
                    f"Error getting gene list back, status code: {response.status_code}"
                )
            returnedL = json.loads(response.text)["genes"]
            returnedN = sum(1 for gene in gene_list if gene in returnedL)
            self._logger.info(f"{returnedN} genes successfully recognized by Enrichr")
        except requests.exceptions.RequestException as e:
            raise EnrichrNetworkError(f"Network error checking genes: {e}")

    def check_uppercase(self, gene_list: List[str]) -> bool:
        """
        Check whether a list of gene names are mostly in uppercase.

        Parameters
        ----------
        gene_list : list
            A list of gene names

        Returns
        -------
        bool
            Whether the list of gene names are mostly in uppercase
        """
        if all(self._is_entrez_id(s) for s in gene_list):
            return False
        is_upper = [str(s).isupper() for s in gene_list]
        return sum(is_upper) / len(is_upper) >= 0.9

    def _parse_enrichr_response(self, data: Dict, gene_set_key: str) -> pd.DataFrame:
        """Parse Enrichr JSON response into DataFrame (extracted to avoid duplication)."""
        colnames = [
            "Rank",
            "Term",
            "P-value",
            "Odds Ratio",
            "Combined Score",
            "Genes",
            "Adjusted P-value",
            "Old P-value",
            "Old adjusted P-value",
        ]
        res = pd.DataFrame(data[gene_set_key], columns=colnames)
        res["Genes"] = res["Genes"].apply(";".join)
        colord = [
            "Term",
            "P-value",
            "Adjusted P-value",
            "Old P-value",
            "Old adjusted P-value",
            "Odds Ratio",
            "Combined Score",
            "Genes",
        ]
        return res.loc[:, colord]

    def get_results_with_background(
        self, gene_list: List[str], background: List[str]
    ) -> Tuple[str, pd.DataFrame]:
        """Get enrichment results with custom background."""
        # Add gene list
        addlist_url = f"{ENRICHR_URL_SPEED}/api/addList"
        payload = {"list": (None, gene_list), "description": (None, self.descriptions)}
        job_id = self.send_genes(payload, addlist_url)

        # Add background list
        addbg_url = f"{ENRICHR_URL_SPEED}/api/addbackground"
        payload = dict(background="\n".join(background))
        bg_id = self.send_background(payload, addbg_url)

        # Get background enrich result
        bgenr_url = f"{ENRICHR_URL_SPEED}/api/backgroundenrich"
        payload = dict(
            userListId=job_id["userListId"],
            backgroundid=bg_id["backgroundid"],
            backgroundType=self._gs,
        )
        s = retry(num=5)
        response = s.post(bgenr_url, data=payload)
        if not response.ok:
            self._logger.error(f"Error fetching enrichment results: {self._gs}")
            raise EnrichrAPIError(
                f"Error fetching enrichment results for {self._gs}, status: {response.status_code}"
            )
        data = json.loads(response.content)
        res = self._parse_enrichr_response(data, self._gs)
        return (job_id["shortId"], res)

    def get_results(self, gene_list: List[str]) -> Tuple[str, pd.DataFrame]:
        """Get enrichment results from Enrichr API."""
        addlist_url = f"{ENRICHR_URL}/{self._organism}/addList"
        payload = {"list": (None, gene_list), "description": (None, self.descriptions)}
        job_id = self.send_genes(payload, addlist_url)
        user_list_id = job_id["userListId"]

        results_url = f"{ENRICHR_URL}/{self._organism}/export"
        filename = f"{self._gs}.{self.descriptions}.reports"
        url = f"{results_url}?userListId={user_list_id}&filename={filename}&backgroundType={self._gs}"

        s = retry(num=5)
        response = s.get(url, stream=True)
        response.encoding = "utf-8"
        if not response.ok:
            self._logger.error(f"Error fetching enrichment results: {self._gs}")

        try:
            res = pd.read_csv(StringIO(response.text), sep="\t")
        except pd.errors.ParserError as parse_err:
            # Fallback to JSON endpoint
            self._logger.warning(
                f"CSV parsing failed: {parse_err}. Trying JSON endpoint."
            )
            fallback_url = f"{ENRICHR_URL}/Enrichr/enrich?userListId={user_list_id}&backgroundType={self._gs}"
            response = s.get(fallback_url)
            if not response.ok:
                self._logger.error(f"Error fetching enrichment results: {self._gs}")
                raise EnrichrAPIError(
                    f"Error fetching enrichment results for {self._gs}, status: {response.status_code}"
                )
            data = json.loads(response.text)
            res = self._parse_enrichr_response(data, self._gs)

        return (job_id["shortId"], res)

    def _is_entrez_id(self, idx: Union[int, str]) -> bool:
        """Check if the input is a valid Entrez ID (integer)."""
        try:
            int(idx)
            return True
        except (ValueError, TypeError):
            return False

    @lru_cache(maxsize=128)
    def _get_libraries_cached(self, organism: str) -> Tuple[str, ...]:
        """Cached helper to fetch libraries (returns tuple for hashability)."""
        lib_url = f"{ENRICHR_URL}/{organism}/datasetStatistics"
        s = retry(num=5)
        try:
            response = s.get(lib_url, verify=True)
            if not response.ok:
                raise EnrichrAPIError(
                    f"Error getting the Enrichr libraries, status: {response.status_code}"
                )
            libs_json = json.loads(response.text)
            libs = [lib["libraryName"] for lib in libs_json["statistics"]]
            return tuple(sorted(libs))
        except requests.exceptions.RequestException as e:
            raise EnrichrNetworkError(f"Network error getting libraries: {e}")

    def get_libraries(self) -> List[str]:
        """Return active enrichr library names (Official API with caching)."""
        return list(self._get_libraries_cached(self._organism))

    def get_background(self) -> Set[str]:
        """Get background genes from file or BioMart."""
        # input is a file
        if os.path.isfile(self.background):
            with open(self.background) as b:
                bg = [g.strip() for g in b]
            return set(bg)

        # package included data
        db_file = os.path.join(
            DEFAULT_CACHE_PATH, f"{self.background}.background.genes.txt"
        )
        if os.path.exists(db_file):
            df = pd.read_csv(db_file, sep="\t")
        else:
            # background is a biomart database name
            self._logger.warning(
                f"Downloading {self.background} for the first time. It might take a couple of minutes."
            )
            bm = Biomart()
            df = bm.query(
                dataset=self.background,
                attributes=["ensembl_gene_id", "external_gene_name", "entrezgene_id"],
                filename=db_file,
            )
        self._logger.info(
            f"Using all annotated genes with GO_ID as background: {self.background}"
        )
        # input id type: entrez or gene_name
        if self._isezid:
            df.dropna(subset=["entrezgene_id"], inplace=True)
            bg = df["entrezgene_id"].astype(int)
        else:
            df.dropna(subset=["external_gene_name"], inplace=True)
            bg = df["external_gene_name"]

        return set(bg)

    def set_organism(self) -> None:
        """Select Enrichr organism from supported list.

        Supported organisms:
        - Human & Mouse, H. sapiens & M. musculus
        - Fly, D. melanogaster
        - Yeast, S. cerevisiae
        - Worm, C. elegans
        - Fish, D. rerio
        """
        if self.organism.lower() in DEFAULT_ORGANISMS:
            self._organism = "Enrichr"
            return

        for k, v in ORGANISM_MAPPINGS.items():
            if self.organism.lower() in v:
                self._organism = f"{k}Enrichr"
                return

        if self._organism is None:
            raise EnrichrValidationError(
                f"No supported organism found for: {self.organism}. "
                f"Supported organisms: {', '.join(DEFAULT_ORGANISMS + [x for v in ORGANISM_MAPPINGS.values() for x in v])}"
            )

        # Verify organism endpoint is accessible
        enrichr_server = f"{ENRICHR_URL}/{self._organism}"
        s = retry(num=3)
        try:
            if s.get(enrichr_server, verify=True).ok:
                return
        except requests.exceptions.RequestException:
            pass

        raise EnrichrNetworkError(
            f"Cannot connect to Enrichr server for organism: {self._organism}"
        )

    def filter_gmt(
        self, gmt: Dict[str, List[str]], background: Set[str]
    ) -> Dict[str, List[str]]:
        """Filter GMT to only include genes that exist in background.

        This substantially affects the significance of the hypergeometric test.

        Parameters
        ----------
        gmt : dict
            A dict of gene sets
        background : set
            A set of custom background genes

        Returns
        -------
        dict
            Filtered gene sets
        """
        gmt_filtered = {}
        for term, genes in gmt.items():
            filtered_genes = [g for g in genes if g in background]
            if filtered_genes:
                gmt_filtered[term] = filtered_genes
        return gmt_filtered

    def parse_background(
        self, gmt: Optional[Dict[str, List[str]]] = None
    ) -> Union[Set[str], int]:
        """Parse and set background genes."""
        self._bg = set()
        if self.background is None:
            # use all genes in the dict input as background if background is None
            if gmt:
                bg = set()
                for term, genes in gmt.items():
                    bg = bg.union(set(genes))
                self._logger.info(
                    f"  Background is not set! Use all {len(bg)} genes in {self._gs}."
                )
                self._bg = bg
        elif np.isscalar(self.background):
            if isinstance(self.background, int) or self.background.isdigit():
                self._bg = int(self.background)
                self._logger.info(f"  Background input is a number: {self._bg}")
            elif isinstance(self.background, str):
                self._bg = self.get_background()
                self._logger.info(f"  Background: found {len(self._bg)} genes")
            else:
                raise EnrichrValidationError("Unsupported background data type")
        else:
            # handle array object: nd.array, list, tuple, set, Series
            try:
                _ = iter(self.background)
                self._bg = set(self.background)
            except TypeError:
                self._logger.error("  Unsupported background data type")
                raise EnrichrValidationError("Unsupported background data type")

        return self._bg

    def enrich(self, gmt: Dict[str, List[str]]) -> Optional[pd.DataFrame]:
        """Perform local enrichment analysis using hypergeometric test.

        p-value: computed using the Fisher exact test (Hypergeometric test)
        z-score: Odds Ratio
        combined score: -log(p) * z

        See: http://amp.pharm.mssm.edu/Enrichr/help#background&q=4

        Columns: Term, Overlap, P-value, Odds Ratio, Combined Score, Adjusted_P-value, Genes
        """
        # Check gene case consistency
        top10 = min(len(gmt), 10)
        top10_keys = list(gmt.keys())[:top10]
        _gls = self._gls
        ups = [self.check_uppercase(gmt[key]) for key in top10_keys]
        _gene_toupper = False

        if all(ups) and not self._gene_isupper:
            _gls = [str(x).upper() for x in self._gls]
            _gene_toupper = True
            self._logger.info(
                "  Genes in GMT file are all in upper case, convert query to upper case."
            )

        bg = self.parse_background(gmt)
        if isinstance(bg, set) and not self._gene_isupper and _gene_toupper:
            bg = {str(s).upper() for s in bg}

        # Statistical testing
        hgtest = list(calc_pvalues(query=_gls, gene_sets=gmt, background=bg))
        if len(hgtest) > 0:
            terms, pvals, oddr, olsz, gsetsz, genes = hgtest
            fdrs, _rej = multiple_testing_correction(
                ps=pvals, alpha=self.cutoff, method="benjamini-hochberg"
            )
            # Build result DataFrame (dict maintains insertion order in Python 3.7+)
            res = pd.DataFrame(
                {
                    "Term": terms,
                    "Overlap": [f"{h}/{g}" for h, g in zip(olsz, gsetsz)],
                    "P-value": pvals,
                    "Adjusted P-value": fdrs,
                    "Odds Ratio": oddr,
                    "Combined Score": -1 * np.log(pvals) * oddr,
                    "Genes": [";".join(map(str, g)) for g in genes],
                }
            )
            return res
        return None

    def _process_single_geneset(
        self, name: str, geneset: Union[Dict[str, List[str]], str], genes_list: str
    ) -> Optional[pd.DataFrame]:
        """Process a single gene set enrichment (extracted for parallelization)."""
        self._logger.info(f"Run: {name}")
        self._gs = name

        if isinstance(geneset, dict):
            # Local mode - offline enrichment
            self._logger.debug(f"Off-line enrichment analysis with library: {name}")
            if self._isezid:
                geneset = {k: list(map(int, v)) for k, v in geneset.items()}
            res = self.enrich(geneset)
            if res is None:
                self._logger.info(f"No hits returned for library: {name}")
                return None
        else:
            # Online mode - API enrichment
            self._logger.debug(f"Enrichr service using library: {name}")
            bg = self.parse_background()
            # Whether user input background
            if isinstance(bg, set) and len(bg) > 0:
                _shortID, res = self.get_results_with_background(genes_list, list(bg))
            else:
                _shortID, res = self.get_results(genes_list)

        # Add gene set library name to results
        res.insert(0, "Gene_set", name)
        return res

    def _save_results(self, res: pd.DataFrame, name: str) -> None:
        """Save enrichment results to file and generate plots."""
        if self._outdir is None:
            return

        outfile = f"{self.outdir}/{name}.{self.organism}.{self.module}.reports.txt"
        res.to_csv(
            outfile, index=False, encoding="utf-8", float_format="%.6e", sep="\t"
        )
        self._logger.info(f"Save enrichment results for {name}")

        # Generate plots
        if not self.__no_plot:
            _ax = barplot(
                df=res,
                cutoff=self.cutoff,
                figsize=self.figsize,
                top_term=self.__top_term,
                color="salmon",
                title=name,
                ofname=outfile.replace("txt", self.format),
            )
            self._logger.debug("Generate figures")

    def _validate_inputs(self) -> Tuple[str, List[Union[Dict[str, List[str]], str]]]:
        """Validate and parse inputs."""
        # Parse gene list
        genes_list = self.parse_genelists()

        if not self._isezid:
            self._gene_isupper = self.check_uppercase(self._gls)

        # Parse gene sets
        gss = self.parse_genesets()
        if len(gss) < 1:
            self._logger.error(f"None of your input gene set matched! {self.gene_sets}")
            self._logger.error(
                f"Hint: Current organism = {self.organism}, is this correct?\n"
                "Hint: use get_library_name() to view full list of supported names."
            )
            raise LookupError(
                "Not validated Enrichr library! Please provide correct organism and library name!"
            )

        return genes_list, gss

    def run(self) -> None:
        """Run enrichr for one sample gene list against multiple libraries.

        This method now supports parallel processing of multiple gene set libraries.
        """
        # Validate inputs
        genes_list, gss = self._validate_inputs()

        self.results = []
        all_online = all(isinstance(g, str) for g in gss)

        # Process gene sets (with parallelization for online mode)
        if all_online and len(gss) > 1:
            # Parallel processing for online API calls
            self._logger.info(f"Processing {len(gss)} gene sets in parallel")
            with ThreadPoolExecutor(max_workers=min(len(gss), 10)) as executor:
                future_to_name = {
                    executor.submit(
                        self._process_single_geneset, name, g, genes_list
                    ): name
                    for name, g in zip(self._gs_name, gss)
                }
                for future in as_completed(future_to_name):
                    name = future_to_name[future]
                    try:
                        res = future.result()
                        if res is not None:
                            self.results.append(res)
                            self.res2d = res
                            self._save_results(res, name)
                    except Exception as exc:
                        self._logger.error(
                            f"Gene set {name} generated an exception: {exc}"
                        )
        else:
            # Sequential processing for local mode or single gene set
            for name, g in zip(self._gs_name, gss):
                res = self._process_single_geneset(name, g, genes_list)
                if res is not None:
                    self.results.append(res)
                    self.res2d = res
                    self._save_results(res, name)

        # Combine all results
        if len(self.results) == 0:
            self._logger.error("No hits returned for all input gene sets!")
            return

        self.results = pd.concat(self.results, ignore_index=True)
        self._logger.info("Done.")
