import errno
import logging
import os
import re
import sys
from typing import Dict, List, Optional

import pandas as pd
import requests
from requests.adapters import HTTPAdapter

try:
    # Prefer direct urllib3 import
    from urllib3.util.retry import Retry
except Exception:  # pragma: no cover
    # Fallback for environments relying on requests' vendored urllib3
    from requests.packages.urllib3.util.retry import Retry  # type: ignore

DEFAULT_CACHE_PATH = os.path.join(os.path.expanduser("~"), ".cache/gseapy")


def unique(seq):
    """Remove duplicates from a list in Python while preserving order.

    :param seq: a python list object.
    :return: a list without duplicates while preserving order.

    """

    seen = set()
    seen_add = seen.add
    """
    The fastest way to sovle this problem is here
    Python is a dynamic language, and resolving seen.add each iteration
    is more costly than resolving a local variable. seen.add could have
    changed between iterations, and the runtime isn't smart enough to rule
    that out. To play it safe, it has to check the object each time.
    """

    return [x for x in seq if x not in seen and not seen_add(x)]


def mkdirs(outdir):
    """create new directory"""
    try:
        os.makedirs(outdir)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise exc
        pass


def log_init(name, log_level=logging.INFO, filename=None):
    """logging

    :param name: logger name
    :log_level: refer to logging module
    :filename: if given a filename, write log to a file (only works in commandline)

    """
    # inside python console environment. only need one root logger
    if hasattr(sys, "ps1"):
        name = "gseapy"
    logger = logging.getLogger(name)
    # clear old root logger handlers
    if logger.hasHandlers():
        log_close(logger)
        # logger.handlers = []
    logger.setLevel(logging.DEBUG)
    logger.propagate = False  # don't find root logger
    # init a root logger
    # logging.basicConfig(
    #     level=logging.DEBUG,
    #     format="%(asctime)s %(name)s::[%(levelname)-8s] %(message)s",
    #     # filename=filename,
    #     # filemode="w",
    # )
    # # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(log_level)
    # set a format which is simpler for console use
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add handlers
    logger.addHandler(console)
    # only write log file when in command line
    if (not hasattr(sys, "ps1")) and filename:
        fhandler = logging.FileHandler(filename=filename, mode="w")
        fhandler.setLevel(logging.DEBUG)
        formatter = logging.Formatter("%(asctime)s %(name)s::[%(levelname)-8s] %(message)s")
        fhandler.setFormatter(formatter)
        logger.addHandler(fhandler)
    # logger.handlers.clear()
    # logger.removeHandler(fh)
    return logger


def log_close(logger):
    handlers = logger.handlers[:]
    for handler in handlers:
        logger.removeHandler(handler)
        handler.close()


def retry(num=5, pool_maxsize: int = 50):
    """retry connection with keep-alive and pooling.

    - Retries on common 5xx errors
    - Enables connection pooling for both HTTP and HTTPS
    - Sets a larger pool size to support concurrent requests
    - Accepts gzip/deflate to speed up transfers
    """
    s = requests.Session()
    retries = Retry(total=num, backoff_factor=1, status_forcelist=[500, 502, 503, 504])
    adapter = HTTPAdapter(max_retries=retries, pool_connections=pool_maxsize, pool_maxsize=pool_maxsize)
    s.mount("http://", adapter)
    s.mount("https://", adapter)
    # Encourage compressed responses
    s.headers.update({"Accept-Encoding": "gzip, deflate"})
    return s


class GOFilter:
    """Filter GO enrichment results by GO term level (depth in the hierarchy).

    The GO level of a term is defined as the number of ``is_a`` ancestors in
    the Gene Ontology hierarchy (root = level 0, direct children of root =
    level 1, etc.).  Ancestor counts are fetched in batches from the
    `EBI QuickGO <https://www.ebi.ac.uk/QuickGO/>`_ REST API.

    This is analogous to the ``gofilter`` function in R's clusterProfiler
    (https://rdrr.io/bioc/clusterProfiler/man/gofilter.html).

    Examples
    --------
    >>> import gseapy
    >>> enr = gseapy.enrichr(gene_list=my_genes,
    ...                      gene_sets="GO_Biological_Process_2021",
    ...                      organism="human", no_plot=True, outdir=None)
    >>> gf = gseapy.GOFilter()
    >>> filtered = gf.filter(enr.results, min_level=3, max_level=8)
    """

    _QUICKGO_URL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{ids}/ancestors"
    _BATCH_SIZE = 100  # QuickGO accepts comma-separated IDs in one request

    def __init__(self):
        self._logger = log_init(__name__)

    def _extract_go_ids(self, terms: "pd.Series") -> List[Optional[str]]:
        """Extract GO IDs from a Series of term strings.

        GO IDs are expected to appear in the format ``GO:XXXXXXX`` anywhere
        within the term string, e.g. ``"response to stimulus (GO:0050896)"``.

        Parameters
        ----------
        terms : pd.Series
            Series of GO term strings from enrichment results.

        Returns
        -------
        list of str or None
            A list with a GO ID string for each term, or ``None`` when no
            GO ID is found in that term.
        """
        go_ids = []
        for term in terms:
            match = re.search(r"(GO:\d+)", str(term))
            go_ids.append(match.group(1) if match else None)
        return go_ids

    def _get_go_levels(self, go_ids: List[str]) -> Dict[str, int]:
        """Fetch GO term levels from the QuickGO REST API.

        The level of a GO term is the number of ``is_a`` ancestors it has in
        the Gene Ontology hierarchy (root = level 0, direct children of root =
        level 1, etc.).

        Parameters
        ----------
        go_ids : list of str
            Unique GO IDs to look up (e.g. ``["GO:0008150", "GO:0009987"]``).

        Returns
        -------
        dict
            Mapping of GO ID -> level (int).  IDs that could not be resolved
            are omitted from the result.
        """
        unique_ids = list(dict.fromkeys(go_ids))  # deduplicate, preserve order
        levels: Dict[str, int] = {}

        for i in range(0, len(unique_ids), self._BATCH_SIZE):
            batch = unique_ids[i : i + self._BATCH_SIZE]
            ids_str = ",".join(batch)
            url = self._QUICKGO_URL.format(ids=ids_str)
            try:
                resp = requests.get(
                    url,
                    params={"relations": "is_a"},
                    headers={"Accept": "application/json"},
                    timeout=30,
                )
                if not resp.ok:
                    self._logger.warning(
                        "QuickGO API returned status %s for batch starting at index %d.",
                        resp.status_code,
                        i,
                    )
                    continue
                data = resp.json()
                for result in data.get("results", []):
                    go_id = result.get("id")
                    ancestors = result.get("ancestors", [])
                    if go_id:
                        levels[go_id] = len(ancestors)
            except Exception as exc:
                self._logger.warning(
                    "Failed to retrieve GO levels for batch %s...: %s",
                    batch[:3],
                    exc,
                )

        return levels

    def filter(
        self,
        df: "pd.DataFrame",
        min_level: int = 1,
        max_level: int = 20,
    ) -> "pd.DataFrame":
        """Filter enrichment results by GO term level.

        Only terms whose GO level falls within ``[min_level, max_level]``
        (inclusive) are retained.  Terms without a recognisable GO ID in the
        ``Term`` column are kept unchanged.

        .. note::
           Requires internet access to query the QuickGO API.
           Has no effect when the ``Term`` column contains no ``GO:XXXXXXX``
           identifiers (e.g. non-GO enrichment libraries).

        Parameters
        ----------
        df : pd.DataFrame
            Enrichment result DataFrame (must have a ``Term`` column), e.g.
            ``enr.res2d`` or ``enr.results``.
        min_level : int
            Minimum GO level to keep (inclusive).  Raise this value to exclude
            very general (high-level) terms.  Default: ``1``.
        max_level : int
            Maximum GO level to keep (inclusive).  Lower this value to exclude
            very specific (low-level) terms.  Default: ``20``.

        Returns
        -------
        pd.DataFrame
            A filtered copy of *df* with index reset.
        """
        if df is None or df.empty:
            self._logger.warning("GOFilter.filter: No results to filter.")
            return df

        go_ids = self._extract_go_ids(df["Term"])
        valid_ids = [gid for gid in go_ids if gid is not None]

        if not valid_ids:
            self._logger.warning(
                "GOFilter.filter: No GO IDs found in the Term column. "
                "This filter only applies to GO enrichment results."
            )
            return df

        levels = self._get_go_levels(valid_ids)
        if not levels:
            self._logger.warning(
                "GOFilter.filter: Could not retrieve GO levels from QuickGO API. Returning unfiltered results."
            )
            return df

        keep_mask = []
        for go_id in go_ids:
            if go_id is None:
                keep_mask.append(True)
                continue
            level = levels.get(go_id)
            if level is None:
                # Level unknown - keep the term to avoid silent data loss
                keep_mask.append(True)
                continue
            keep_mask.append(min_level <= level <= max_level)

        return df.loc[keep_mask].reset_index(drop=True)


# CONSTANT
DEFAULT_LIBRARY = [
    "GO_Biological_Process_2013",
    "GO_Biological_Process_2015",
    "GO_Cellular_Component_2013",
    "GO_Cellular_Component_2015",
    "GO_Molecular_Function_2013",
    "GO_Molecular_Function_2015",
    "GeneSigDB",
    "HumanCyc_2015",
    "Human_Gene_Atlas",
    "Human_Phenotype_Ontology",
    "Humancyc_2016",
    "KEGG_2013",
    "KEGG_2015",
    "KEGG_2016",
    "MGI_Mammalian_Phenotype_2013",
    "MGI_Mammalian_Phenotype_Level_3",
    "MGI_Mammalian_Phenotype_Level_4",
    "MSigDB_Computational",
    "MSigDB_Oncogenic_Signatures",
    "Mouse_Gene_Atlas",
    "Panther_2015",
    "Panther_2016",
    "Reactome_2013",
    "Reactome_2015",
    "Reactome_2016",
    "WikiPathways_2013",
    "WikiPathways_2015",
    "WikiPathways_2016",
]
