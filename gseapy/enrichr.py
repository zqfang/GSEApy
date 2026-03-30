#!/usr/bin/env python
# -*- coding: utf-8 -*-
# see: http://amp.pharm.mssm.edu/Enrichr/help#api for API docs
import io
import json
import logging
import os
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple, Union

import numpy as np
import pandas as pd
import requests

from gseapy.biomart import Biomart
from gseapy.plot import barplot
from gseapy.stats import calc_pvalues, multiple_testing_correction
from gseapy.utils import DEFAULT_CACHE_PATH, log_init, mkdirs, retry


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


class EnrichrAPI:
    """A Python client for the modEnrichr suite and Speedrichr REST APIs."""

    def __init__(self, organism: str = "human"):
        """
        Initializes the API client for a specific organism.

        Args:
            organism (str): The target organism. Options include 'human', 'mouse',
                            'fly', 'yeast', 'worm', 'fish'. Defaults to 'human'.
        """
        organism = organism.lower().strip()

        # Map common names to the correct modEnrichr instance
        self.db_map = {
            "human": "Enrichr",
            "hsapiens": "Enrichr",
            "homo sapiens": "Enrichr",
            "hs": "Enrichr",
            "h. sapiens": "Enrichr",
            "mus musculus": "Enrichr",
            "m. musculus": "Enrichr",
            "mouse": "Enrichr",
            "mm": "Enrichr",
            "enrichr": "Enrichr",
            "fly": "FlyEnrichr",
            "drosophila": "FlyEnrichr",
            "drosophila melanogaster": "FlyEnrichr",
            "d. melanogaster": "FlyEnrichr",
            "yeast": "YeastEnrichr",
            "saccharomyces": "YeastEnrichr",
            "s. cerevisiae": "YeastEnrichr",
            "saccharomyces cerevisiae": "YeastEnrichr",
            "worm": "WormEnrichr",
            "celegans": "WormEnrichr",
            "c. elegans": "WormEnrichr",
            "caenorhabditis elegans": "WormEnrichr",
            "nematode": "WormEnrichr",
            "fish": "FishEnrichr",
            "zebrafish": "FishEnrichr",
            "danio rerio": "FishEnrichr",
            "d. rerio": "FishEnrichr",
        }
        if organism not in self.db_map:
            valid_opts = ", ".join(sorted(set(self.db_map.keys())))
            raise ValueError(
                f"Invalid organism '{organism}'. Valid options are: {valid_opts}"
            )

        instance = self.db_map[organism]

        # Set the dynamic base URL
        self.base_url = f"https://maayanlab.cloud/{instance}"

        # Speedrichr URL (primarily supports standard Enrichr)
        self.speedrichr_url = "https://maayanlab.cloud/speedrichr/api"

        # Shared session with retry/pooling for all HTTP calls
        self._session = retry()

    def _ensure_list(self, item: Union[str, List[str]]) -> List[str]:
        """Utility to convert a single string to a list for uniform processing."""
        return [item] if isinstance(item, str) else item

    # ------------------------------------------------------------------------
    # Info & Utility Endpoints
    # ------------------------------------------------------------------------

    def get_libraries(self) -> List[str]:
        """Fetches a list of all available gene set library names for the current organism."""
        url = f"{self.base_url}/datasetStatistics"
        response = self._session.get(url)
        response.raise_for_status()
        data = response.json()

        # The API returns a dictionary with a 'statistics' key containing a list of library info
        return [library["libraryName"] for library in data.get("statistics", [])]

    def download_libraries(self, libname: str) -> Dict[str, List[str]]:
        """Download enrichr libraries"""
        # query string
        url = f"{self.base_url}/geneSetLibrary?mode=text&libraryName={libname}"
        # get
        response = self._session.get(url, timeout=None, stream=True)
        if not response.ok:
            raise Exception(
                "Error fetching gene set library. Please verify that the library name and organism are correct."
            )
        # reformat to dict and save to disk
        mkdirs(DEFAULT_CACHE_PATH)
        genesets_dict = {}
        outname = "Enrichr.%s.gmt" % libname  # pattern: database.library.gmt
        output_path = os.path.join(DEFAULT_CACHE_PATH, outname)
        temp_path = output_path + ".tmp"
        with open(temp_path, "w") as gmtout:
            for line in response.iter_lines(chunk_size=1024, decode_unicode="utf-8"):
                line = line.strip().split("\t")
                k = line[0]
                v = map(lambda x: x.split(",")[0], line[2:])
                v = list(filter(lambda x: True if len(x) else False, v))
                genesets_dict[k] = v
                outline = "%s\t%s\t%s\n" % (k, line[1], "\t".join(v))
                gmtout.write(outline)
        os.replace(temp_path, output_path)

        return genesets_dict

    def find_terms_by_gene(
        self, gene: str, include_json: bool = True, include_setup: bool = True
    ) -> Dict[str, Any]:
        """Finds terms that contain a given gene across Enrichr libraries."""
        url = f"{self.base_url}/genemap"
        params = {
            "gene": gene,
            "json": str(include_json).lower(),
            "setup": str(include_setup).lower(),
        }
        response = self._session.get(url, params=params)
        response.raise_for_status()
        return response.json()

    # ------------------------------------------------------------------------
    # Standard Enrichment API
    # ------------------------------------------------------------------------

    def add_list(
        self, genes: List[str], description: str = "Gene list"
    ) -> Dict[str, Any]:
        """Analyzes a gene set and returns a userListId."""
        url = f"{self.base_url}/addList"
        payload = {
            "list": (None, "\n".join(genes)),
            "description": (None, description),
        }
        response = self._session.post(url, files=payload)
        response.raise_for_status()
        return response.json()

    def view_list(self, user_list_id: int) -> Dict[str, Any]:
        """Views an added gene set using its userListId."""
        url = f"{self.base_url}/view"
        params = {"userListId": user_list_id}
        response = self._session.get(url, params=params)
        response.raise_for_status()
        return response.json()

    def get_enrichment(
        self, user_list_id: int, gene_set_library: Union[str, List[str]]
    ) -> Dict[str, Any]:
        """Gets enrichment results for one or more gene set libraries."""
        url = f"{self.base_url}/enrich"
        libraries = self._ensure_list(gene_set_library)
        all_results = {}

        for lib in libraries:
            params = {"userListId": user_list_id, "backgroundType": lib}
            response = self._session.get(url, params=params)
            response.raise_for_status()
            all_results[lib] = response.json()

        return all_results

    def get_results_dataframe(
        self, user_list_id: int, gene_set_library: Union[str, List[str]]
    ) -> pd.DataFrame:
        """Fetches enrichment analysis results and returns a combined pandas DataFrame."""
        url = f"{self.base_url}/export"
        libraries = self._ensure_list(gene_set_library)

        all_dfs = []

        for lib in libraries:
            params = {
                "userListId": user_list_id,
                "filename": "temp",
                "backgroundType": lib,
            }
            response = self._session.get(url, params=params)
            response.raise_for_status()

            df = pd.read_csv(io.StringIO(response.text), sep="\t")
            cols = df.columns.to_list()
            df["Gene_set"] = lib
            cols = ["Gene_set"] + cols  # reorder columns to put Gene_set first
            df = df[cols]
            all_dfs.append(df)
        if all_dfs:
            return pd.concat(all_dfs, ignore_index=True)
        else:
            return pd.DataFrame()

    # ------------------------------------------------------------------------
    # Background Enrichment API (Speedrichr)
    # ------------------------------------------------------------------------

    def add_list_speedrichr(
        self, genes: List[str], description: str = "Gene list with background"
    ) -> Dict[str, Any]:
        """Uploads a gene set to Speedrichr for background analysis."""
        url = f"{self.speedrichr_url}/addList"
        payload = {
            "list": (None, "\n".join(genes)),
            "description": (None, description),
        }
        response = self._session.post(url, files=payload)
        response.raise_for_status()
        return response.json()

    def add_background(self, background_genes: List[str]) -> Dict[str, Any]:
        """Uploads a background gene set to Speedrichr."""
        url = f"{self.speedrichr_url}/addbackground"
        payload = {"background": "\n".join(background_genes)}
        response = self._session.post(url, data=payload)
        response.raise_for_status()
        return response.json()

    def get_background_enrichment(
        self,
        user_list_id: int,
        background_id: str,
        gene_set_library: Union[str, List[str]],
    ) -> Dict[str, Any]:
        """Gets enrichment results calculated against a custom background for one or more libraries."""
        url = f"{self.speedrichr_url}/backgroundenrich"
        libraries = self._ensure_list(gene_set_library)
        all_results = []

        for lib in libraries:
            payload = {
                "userListId": user_list_id,
                "backgroundid": background_id,
                "backgroundType": lib,
            }
            response = self._session.post(url, data=payload)
            response.raise_for_status()
            df = self._parse_json_response(response.json(), lib)
            cols = df.columns.to_list()
            df["Gene_set"] = lib
            cols = ["Gene_set"] + cols  # reorder columns to put Gene_set first
            df = df[cols]
            all_results.append(df)
        all_results = (
            pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()
        )
        return all_results

    def _parse_json_response(self, data: Dict, gene_set_key: str) -> pd.DataFrame:
        """Parse Enrichr JSON response into DataFrame."""
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


class Enrichr(EnrichrAPI):
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
        self._gene_isupper = True
        self._gene_toupper = False
        self._gls: Union[Set[int], List[str]] = []
        self._isezid = False
        self._gs: str = ""
        self._gs_name: List[str] = []
        super().__init__(organism)
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

    def set_organism(self) -> None:
        """Initialize EnrichrAPI base with the selected organism, setting base_url."""
        if self.organism not in self.db_map:
            valid_opts = ", ".join(sorted(set(self.db_map.keys())))
            raise ValueError(
                f"Invalid organism '{self.organism}'. Valid options are: {valid_opts}"
            )
        instance = self.db_map[self.organism]
        # Set the dynamic base URL
        self.base_url = f"https://maayanlab.cloud/{instance}"

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

    def _is_entrez_id(self, idx: Union[int, str]) -> bool:
        """Check if the input is a valid Entrez ID (integer)."""
        try:
            int(idx)
            return True
        except (ValueError, TypeError):
            return False

    def send_genes(self, payload, url) -> Dict[str, Union[int, str]]:
        """Send gene list to enrichr server."""
        try:
            response = self._session.post(url, files=payload, verify=True)
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

    def check_genes(self, gene_list: List[str], usr_list_id: str) -> None:
        """Compare the genes sent and received to get successfully recognized genes."""
        try:
            returned = self.view_list(usr_list_id)
            returnedL = returned["genes"]
            returnedN = sum(1 for gene in gene_list if gene in returnedL)
            self._logger.info(f"{returnedN} genes successfully recognized by Enrichr")
        except requests.exceptions.RequestException as e:
            raise EnrichrNetworkError(f"Network error checking genes: {e}")

    def get_results_with_background(
        self, gene_list: str, background: List[str], gene_set_libraries: List[str]
    ) -> Tuple[str, pd.DataFrame]:
        """Get enrichment results with custom background."""
        payload = {"list": (None, gene_list), "description": (None, self.descriptions)}
        job_id = self.send_genes(payload, f"{self.speedrichr_url}/addList")

        bg_id = self.add_background(background)

        res = self.get_background_enrichment(
            job_id["userListId"], bg_id["backgroundid"], gene_set_libraries
        )
        return (job_id["shortId"], res)

    def get_results(
        self, gene_list: str, gene_set_libraries: List[str]
    ) -> Tuple[str, pd.DataFrame]:
        """Get enrichment results from Enrichr API."""
        payload = {"list": (None, gene_list), "description": (None, self.descriptions)}
        job_id = self.send_genes(payload, f"{self.base_url}/addList")

        res = self.get_results_dataframe(job_id["userListId"], gene_set_libraries)
        return (job_id["shortId"], res)

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
            if df is None:
                raise ConnectionError(
                    f"Failed to download background gene set '{self.background}' from BioMart. "
                    "Check your internet connection or try again later."
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
                for _, genes in gmt.items():
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

    def enrich_local(self, gmt: Dict[str, List[str]]) -> Optional[pd.DataFrame]:
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
            fdrs, _ = multiple_testing_correction(
                ps=pvals, alpha=self.cutoff, method="benjamini-hochberg"
            )
            # Build result DataFrame (dict maintains insertion order in Python 3.7+)
            res = pd.DataFrame(
                {
                    "Gene_set": self._gs,
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

    def enrich_online(
        self, genes_list: str, geneset_libraries: List[str]
    ) -> Optional[pd.DataFrame]:
        """Perform online enrichment analysis using Enrichr API."""

        # Online mode - API enrichment
        self._logger.debug(f"Enrichr service using library: {geneset_libraries}")
        bg = self.parse_background()
        # Whether user input background
        if isinstance(bg, set) and len(bg) > 0:
            _shortID, res = self.get_results_with_background(
                genes_list, list(bg), geneset_libraries
            )
        else:
            _shortID, res = self.get_results(genes_list, geneset_libraries)
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
        """Run enrichr for one sample gene list against multiple libraries."""
        # Validate inputs
        genes_list, gss = self._validate_inputs()

        self.results = []

        online_genesets = [g for g in gss if isinstance(g, str)]
        # run online enrichment for all string libraries together to minimize API calls
        if online_genesets:
            self._logger.info(
                f"Online enrichment analysis with libraries: {', '.join(online_genesets)}"
            )
            res_online = self.enrich_online(genes_list, online_genesets)
            if res_online is not None:
                self.results.append(res_online)

        ## Process local gene sets (dicts) separately to avoid unnecessary API calls
        for name, geneset in zip(self._gs_name, gss):
            if isinstance(geneset, dict):
                # Local mode - offline enrichment
                self._logger.info(f"Off-line enrichment analysis with library: {name}")
                if self._isezid:
                    geneset = {k: list(map(int, v)) for k, v in geneset.items()}
                self._gs = (
                    name  # assign the name of local gene set for background parsing
                )
                res_local = self.enrich_local(geneset)
                if res_local is not None:
                    self.results.append(res_local)

        if len(self.results) == 0:
            self._logger.error("No hits returned for all input gene sets!")
            return

        self.results = pd.concat(self.results, ignore_index=True)
        self.res2d = self.results
        self._save_results(self.results, name="Enrichr")
        self._logger.info("Done.")
