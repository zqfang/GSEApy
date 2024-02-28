#! python
# -*- coding: utf-8 -*-

import json
import logging
import os
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from gseapy.plot import GSEAPlot, TracePlot, gseaplot, heatmap
from gseapy.utils import DEFAULT_CACHE_PATH, log_init, mkdirs, retry


class GMT(dict):
    def __init__(
        self,
        mapping: Optional[Dict[str, str]] = None,
        description: Optional[str] = None,
    ):
        """
        wrapper of dict. this helps merge multiple dict into one
        the original key will changed to new key with suffix '__{description}'
        """
        if description is None:
            description = ""
        self.description = description
        _mapping = {}
        if mapping is not None:
            for key, value in mapping.items():
                k = key + "__" + self.description
                _mapping[k] = value
        super().__init__(_mapping)

    def apply(self, func):
        """apply function in place"""
        for key, value in self.items():
            self[key] = func(value)

    def is_empty(self):
        return len(self) == 0

    def write(self, ofname: str):
        """
        write gmt file to disk
        """
        with open(ofname, "w") as out:
            for key, value in self.items():
                collections = key.split("__")
                collections += list(value)
                out.write("\t".join(collections) + "\n")

    def _read(path):
        mapping = {}
        with open(path, "r") as inp:
            for line in inp:
                items = line.strip().split("\t")
                key = items[0]
                if items[1] != "":
                    key += "__" + items[1]
                mapping[key] = items[2:]
        return mapping

    @classmethod
    def read(cls, paths):
        paths = paths.strip().split(",")
        # mapping
        mapping = {}
        for path in paths:
            mapping.update(cls._read(path))
        return cls(mapping)


class GSEAbase(object):
    """base class of GSEA."""

    def __init__(
        self,
        outdir: Optional[str] = None,
        gene_sets: Union[List[str], str, Dict[str, str]] = "KEGG_2016",
        module: str = "base",
        threads: int = 1,
        enrichr_url: str = "http://maayanlab.cloud",
        verbose: bool = False,
    ):
        self.outdir = outdir
        self.gene_sets = gene_sets
        self.fdr = 0.05
        self.module = module
        self.res2d = None
        self.ranking = None
        self.ascending = False
        self.verbose = verbose
        self._threads = threads
        self.ENRICHR_URL = enrichr_url
        self.pheno_pos = ""
        self.pheno_neg = ""
        self.permutation_num = 0
        self._LIBRARY_LIST_URL = "https://maayanlab.cloud/speedrichr/api/listlibs"

        self._set_cores()
        # init logger
        self.prepare_outdir()

    def __del__(self):
        if hasattr(self, "_logger"):
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
        self._logfile = logfile
        self._logger = log_init(
            name=str(self.module) + str(id(self)),
            log_level=logging.INFO if self.verbose else logging.WARNING,
            filename=logfile,
        )

    def _set_cores(self):
        """set cpu numbers to be used"""

        cpu_num = os.cpu_count() - 1
        if self._threads > cpu_num:
            cores = cpu_num
        elif self._threads < 1:
            cores = 1
        else:
            cores = self._threads
        # have to be int if user input is float
        self._threads = int(cores)

    def _read_file(self, path: str) -> pd.DataFrame:
        """
        read file, and return dataframe (first column are gene IDs)
        """
        # just txt file like input
        header, sep = "infer", "\t"
        # GCT input format?
        if path.endswith(".gct"):
            rank_metric = pd.read_csv(
                path, skiprows=1, comment="#", index_col=0, sep=sep
            )
        else:
            if path.endswith(".csv"):
                sep = ","
            if path.endswith(".rnk"):
                header = None
            rank_metric = pd.read_csv(
                path, comment="#", index_col=0, sep=sep, header=header
            )
            if rank_metric.shape[1] == 1:
                # rnk file like input
                rank_metric.columns = rank_metric.columns.astype(str)

        return rank_metric.select_dtypes(include=[np.number]).reset_index()

    def _load_data(self, exprs: Union[str, pd.Series, pd.DataFrame]) -> pd.DataFrame:
        """
        helper function to read data
        """
        # load data
        if isinstance(exprs, pd.DataFrame):
            rank_metric = exprs.copy()
            # handle dataframe with gene_name as index.
            self._logger.debug("Input data is a DataFrame with gene names")
            # handle index is already gene_names
            if rank_metric.index.dtype == "O":
                rank_metric = rank_metric.reset_index()
                # rank_metric.set_index(keys=rank_metric.columns[0], inplace=True)
            # if rank_metric.columns.dtype != "O":
            rank_metric.columns = rank_metric.columns.astype(str)

        elif isinstance(exprs, pd.Series):
            # change to DataFrame
            self._logger.debug("Input data is a Series with gene names")
            if not isinstance(exprs.name, str):
                if exprs.name is None:
                    # rename col if name attr is none
                    exprs.name = "sample1"
                elif hasattr(exprs.name, "dtype"):
                    if exprs.name.dtype != "O":
                        exprs.name = exprs.name.astype(str)
                else:
                    exprs.name = str(exprs.name)
            rank_metric = exprs.reset_index()
        elif os.path.isfile(exprs):
            rank_metric = self._read_file(exprs)

        else:
            raise Exception("Error parsing expression values!")
        # select numbers
        # rank_metric = rank_metric.select_dtypes(include=[np.number])
        return rank_metric

    def _check_data(self, exprs: pd.DataFrame) -> pd.DataFrame:
        """
        check NAs, duplicates.
        exprs: dataframe, the frist column must be gene identifiers

        return: dataframe, index is gene ids
        """
        ## if gene names contain NA, drop them
        if exprs.iloc[:, 0].isnull().any():
            exprs.dropna(subset=[exprs.columns[0]])
        ## then fill na for numeric columns
        if exprs.isnull().any().sum() > 0:
            self._logger.warning("Input data contains NA, filled NA with 0")
            exprs.dropna(how="all", inplace=True)  # drop rows with all NAs
            exprs = exprs.fillna(0)
        ## check duplicated IDs
        # set gene name as index
        exprs.set_index(keys=exprs.columns[0], inplace=True)
        # select numberic columns
        df = exprs.select_dtypes(include=[np.number])
        # microarray data may contained multiple probs of same gene, average them
        if df.index.duplicated().sum() > 0:
            self._logger.warning(
                "Found duplicated gene names, values averaged by gene names!"
            )
            df = df.groupby(level=0).mean()
        # check whether contains infinity values
        if np.isinf(df).values.sum() > 0:
            self._logger.warning("Input gene rankings contains inf values!")
            col_min_max = {
                np.inf: df[np.isfinite(df)].max(),  # column-wise max
                -np.inf: df[np.isfinite(df)].min(),  # column-wise min
            }
            df = df.replace({col: col_min_max for col in df.columns})
        return df

    def make_unique(self, rank_metric: pd.DataFrame, col_idx: int) -> pd.DataFrame:
        """
        make gene id column unique by adding a digit, similar to R's make.unique
        """
        id_col = rank_metric.columns[col_idx]
        if rank_metric.duplicated(subset=id_col).sum() > 0:
            self._logger.info("Input gene rankings contains duplicated IDs")
            mask = rank_metric.duplicated(subset=id_col, keep=False)
            dups = (
                rank_metric.loc[mask, id_col]
                .to_frame()
                .groupby(id_col)
                .cumcount()
                .map(lambda c: "_" + str(c) if c else "")
            )
            rank_metric.loc[mask, id_col] = rank_metric.loc[mask, id_col] + dups
        return rank_metric

    def load_gmt_only(
        self, gmt: Union[List[str], str, Dict[str, str]]
    ) -> Dict[str, List[str]]:
        """parse gene_sets.
        gmt: List, Dict, Strings

        However,this function will merge different gene sets into one big dict to
        save computation time for later.
        """
        genesets_dict = dict()

        if isinstance(gmt, dict):
            genesets_dict = gmt.copy()
        elif isinstance(gmt, str):
            gmts = gmt.split(",")
            if len(gmts) > 1:
                for gm in gmts:
                    tdt = self.parse_gmt(gm)
                    for k, v in tdt.items():
                        new_k = os.path.split(gm)[-1] + "__" + k
                        genesets_dict[new_k] = v
            else:
                genesets_dict = self.parse_gmt(gmt)
        elif isinstance(gmt, list):
            for i, gm in enumerate(gmt):
                prefix = str(i)
                if isinstance(gm, dict):
                    tdt = gm.copy()
                elif isinstance(gm, str):
                    tdt = self.parse_gmt(gm)
                    prefix = os.path.split(gm)[-1]
                else:
                    continue
                for k, v in tdt.items():
                    new_k = prefix + "__" + k
                    genesets_dict[new_k] = v
        else:
            raise Exception("Error parsing gmt parameter for gene sets")

        if len(genesets_dict) == 0:
            raise Exception("Error parsing gmt parameter for gene sets")
        return genesets_dict

    def load_gmt(
        self, gene_list: Iterable[str], gmt: Union[List[str], str, Dict[str, str]]
    ) -> Dict[str, List[str]]:
        """load gene set dict"""

        genesets_dict = self.load_gmt_only(gmt)

        subsets = list(genesets_dict.keys())
        entry1st = genesets_dict[subsets[0]]
        gene_dict = {g: i for i, g in enumerate(gene_list)}
        for subset in subsets:
            subset_list = set(genesets_dict.get(subset))  # remove duplicates
            # drop genes not found in the gene_dict
            gene_overlap = [g for g in subset_list if g in gene_dict]
            genesets_dict[subset] = gene_overlap
            tag_len = len(gene_overlap)
            if (self.min_size <= tag_len <= self.max_size) and tag_len < len(gene_list):
                # tag_len should < gene_list
                continue
            del genesets_dict[subset]

        filsets_num = len(subsets) - len(genesets_dict)
        self._logger.info(
            "%04d gene_sets have been filtered out when max_size=%s and min_size=%s"
            % (filsets_num, self.max_size, self.min_size)
        )

        if filsets_num == len(subsets):
            msg = (
                "No gene sets passed through filtering condition !!! \n"
                + "Hint 1: Try to lower min_size or increase max_size !\n"
                + "Hint 2: Check gene symbols are identifiable to your gmt input.\n"
                + "Hint 3: Gene symbols curated in Enrichr web services are all upcases.\n"
            )
            self._logger.error(msg)
            dict_head = "{ %s: [%s]}" % (subsets[0], ", ".join(entry1st))
            self._logger.error(
                "The first entry of your gene_sets (gmt) look like this : %s"
                % dict_head
            )
            self._logger.error(
                "The first 5 genes look like this : [ %s ]"
                % (", ".join(list(gene_list)[:5]))
            )
            raise LookupError(msg)

        # self._gmtdct = genesets_dict
        return genesets_dict

    def parse_gmt(self, gmt: str) -> Dict[str, List[str]]:
        """gmt parser when input is a string"""

        if gmt.lower().endswith(".gmt"):
            genesets_dict = {}
            with open(gmt) as genesets:
                for line in genesets:
                    entries = line.strip().split("\t")
                    key = entries[0]
                    genes = [g.split(",")[0] for g in entries[2:]]
                    genesets_dict[key] = genes
            return genesets_dict

        tmpname = "Enrichr." + gmt + ".gmt"
        tempath = os.path.join(DEFAULT_CACHE_PATH, tmpname)
        # if file already download
        if os.path.isfile(tempath):
            self._logger.info(
                "Enrichr library gene sets already downloaded in: %s, use local file"
                % DEFAULT_CACHE_PATH
            )
            return self.parse_gmt(tempath)

        elif gmt in self.get_libraries():
            return self._download_libraries(gmt)
        else:
            self._logger.error("No supported gene_sets: %s" % gmt)
        return dict()

    def get_libraries(self) -> List[str]:
        """return active enrichr library name.Offical API"""

        lib_url = self.ENRICHR_URL + "/Enrichr/datasetStatistics"
        # "LIBRARY_LIST_URL": "https://maayanlab.cloud/speedrichr/api/listlibs",
        s = retry(num=5)
        response = s.get(lib_url, verify=True)
        if not response.ok:
            raise Exception("Error getting the Enrichr libraries")
        libs_json = json.loads(response.text)
        libs = [lib["libraryName"] for lib in libs_json["statistics"]]

        return sorted(libs)

    def _download_libraries(self, libname: str) -> Dict[str, List[str]]:
        """Download enrichr libraries. Only Support Enrichr libraries now"""
        self._logger.info("Downloading and generating Enrichr library gene sets......")
        s = retry(5)
        # queery string
        ENRICHR_URL = self.ENRICHR_URL + "/Enrichr/geneSetLibrary"
        query_string = "?mode=text&libraryName=%s"
        # get
        response = s.get(
            ENRICHR_URL + query_string % libname, timeout=None, stream=True
        )
        if not response.ok:
            raise Exception(
                "Error fetching gene set library, input name is correct for the organism you've set?."
            )
        # reformat to dict and save to disk
        mkdirs(DEFAULT_CACHE_PATH)
        genesets_dict = {}
        outname = "Enrichr.%s.gmt" % libname  # pattern: database.library.gmt
        gmtout = open(os.path.join(DEFAULT_CACHE_PATH, outname), "w")
        for line in response.iter_lines(chunk_size=1024, decode_unicode="utf-8"):
            line = line.strip().split("\t")
            k = line[0]
            v = map(lambda x: x.split(",")[0], line[2:])
            v = list(filter(lambda x: True if len(x) else False, v))
            genesets_dict[k] = v
            outline = "%s\t%s\t%s\n" % (k, line[1], "\t".join(v))
            gmtout.write(outline)
        gmtout.close()

        return genesets_dict

    def _heatmat(self, df: pd.DataFrame, classes: List[str]):
        """only use for gsea heatmap"""

        cls_booA = list(map(lambda x: True if x == self.pheno_pos else False, classes))
        cls_booB = list(map(lambda x: True if x == self.pheno_neg else False, classes))
        datA = df.loc[:, cls_booA]
        datB = df.loc[:, cls_booB]
        datAB = pd.concat([datA, datB], axis=1)
        self.heatmat = datAB
        return

    def _plotting(self, metric: Dict[str, Union[pd.Series, pd.DataFrame]]):
        """Plotting API.
        :param metric: sorted pd.Series with rankings values.
        """

        # no values need to be returned
        if self._outdir is None:
            return
        # indices = self.res2d["NES"].abs().sort_values(ascending=False).index
        indices = self.res2d.index
        # Plotting
        for i, idx in enumerate(indices):
            record = self.res2d.iloc[idx]
            if self.module != "ssgsea" and record["FDR q-val"] > 0.25:
                continue
            if i >= self.graph_num:  # already sorted by abs(NES) in descending order
                break
            # if self.res2d["Name"].nunique() > 1 and hasattr(
            #     self, "_metric_dict"
            # ):  # self.module != "ssgsea":
            #     key = record["Name"]
            #     rank_metric = metric[key]
            #     hit = self._results[key][record["Term"]]["hits"]
            #     RES = self._results[key][record["Term"]]["RES"]
            # else:
            #     rank_metric = metric[self.module]
            #     hit = self._results[record["Term"]]["hits"]
            #     RES = self._results[record["Term"]]["RES"]
            key = record["Name"]
            rank_metric = metric[key]
            hit = self._results[key][record["Term"]]["hits"]
            RES = self._results[key][record["Term"]]["RES"]

            outdir = os.path.join(self.outdir, record["Name"])
            mkdirs(outdir)
            term = record["Term"].replace("/", "-").replace(":", "_")
            outfile = os.path.join(outdir, "{0}.{1}".format(term, self.format))
            if self.module == "gsea":
                outfile2 = "{0}/{1}.heatmap.{2}".format(outdir, term, self.format)
                heatmat = self.heatmat.iloc[hit, :]
                width = np.clip(heatmat.shape[1], 4, 20)
                height = np.clip(heatmat.shape[0], 4, 20)
                heatmap(
                    df=heatmat,
                    title=record["Term"].split("__")[-1],
                    ofname=outfile2,
                    z_score=0,
                    figsize=(width, height),
                    xticklabels=True,
                    yticklabels=True,
                )

            if self.permutation_num > 0:
                # skip plotting when nperm=0
                gseaplot(
                    term=record["Term"].split("__")[-1],
                    hits=hit,
                    nes=record["NES"],
                    pval=record["NOM p-val"],
                    fdr=record["FDR q-val"],
                    RES=RES,
                    rank_metric=rank_metric,
                    pheno_pos=self.pheno_pos,
                    pheno_neg=self.pheno_neg,
                    figsize=self.figsize,
                    ofname=outfile,
                )

    def _to_df(
        self,
        gsea_summary: List[Dict],
        gmt: Dict[str, List[str]],
        metric: Dict[str, pd.Series],
    ) -> pd.DataFrame:
        """Convernt GSEASummary to DataFrame

        rank_metric: Must be sorted in descending order already
        """

        res_df = pd.DataFrame(
            index=range(len(gsea_summary)),
            columns=[
                "name",
                "term",
                "es",
                "nes",
                "pval",
                "fdr",
                "fwerp",
                "tag %",
                "gene %",
                "lead_genes",
                "matched_genes",
                "hits",
                "RES",
            ],
        )
        # res = OrderedDict()

        for i, gs in enumerate(gsea_summary):
            # reformat gene list.
            name = (
                self._metric_dict[str(gs.index)]
                if (gs.index is not None)
                else self.module
            )
            _genes = metric[name].index.values[gs.hits]
            genes = ";".join([str(g).strip() for g in _genes])
            RES = np.array(gs.run_es)
            lead_genes = ""
            tag_frac = ""
            gene_frac = ""
            if len(RES) > 1:
                # extract leading edge genes
                if float(gs.es) >= 0:
                    # RES -> ndarray, ind -> list
                    es_i = RES.argmax()
                    ldg_pos = list(filter(lambda x: x <= es_i, gs.hits))
                    gene_frac = (es_i + 1) / len(metric[name])
                else:
                    es_i = RES.argmin()
                    ldg_pos = list(filter(lambda x: x >= es_i, gs.hits))
                    ldg_pos.reverse()
                    gene_frac = (len(metric[name]) - es_i) / len(metric[name])

                # tag_frac = len(ldg_pos) / len(gmt[gs.term])
                gene_frac = "{0:.2%}".format(gene_frac)
                lead_genes = ";".join(list(map(str, metric[name].iloc[ldg_pos].index)))
                tag_frac = "%s/%s" % (len(ldg_pos), len(gmt[gs.term]))

            e = pd.Series(
                [
                    name,
                    gs.term,
                    gs.es,
                    gs.nes,
                    gs.pval,
                    gs.fdr,
                    gs.fwerp,
                    tag_frac,
                    gene_frac,
                    lead_genes,
                    genes,
                    gs.hits,
                    gs.run_es,
                ],
                index=res_df.columns,
            )
            res_df.iloc[i, :] = e
        return res_df

    def to_df(
        self,
        gsea_summary: List[Dict],
        gmt: Dict[str, List[str]],
        rank_metric: Union[pd.Series, pd.DataFrame],
        indices: Optional[List] = None,
    ):
        """Convernt GSEASummary to DataFrame

        rank_metric: if a Series, then it must be sorted in descending order already
                     if a DataFrame, indices must not None.
        indices: Only works for DataFrame input. Stores the indices of sorted array
        """
        if isinstance(rank_metric, pd.DataFrame) and (indices is not None):
            self._metric_dict = {str(c): n for c, n in enumerate(rank_metric.columns)}
            metric = {
                n: rank_metric.iloc[indices[i], i]  # .sort_values(ascending=False)
                for i, n in enumerate(rank_metric.columns)  # indices is a 2d list
            }
        else:
            metric = {self.module: rank_metric}
            self._metric_dict = {self.module: self.module}

        res_df = self._to_df(gsea_summary, gmt, metric)
        self._results = {}
        # save dict
        # if res_df["name"].nunique() >= 2:
        #     for name, dd in res_df.groupby(["name"]):
        #         self._results[name] = dd.set_index("term").to_dict(orient="index")
        # else:
        #     self._results = res_df.set_index("term").to_dict(orient="index")
        for name, dd in res_df.groupby("name"):
            self._results[name] = dd.set_index("term").to_dict(orient="index")
        # trim
        res_df.rename(
            columns={
                "name": "Name",
                "term": "Term",
                "es": "ES",
                "nes": "NES",
                "pval": "NOM p-val",
                "fdr": "FDR q-val",
                "fwerp": "FWER p-val",
                "tag %": "Tag %",
                "gene %": "Gene %",
                "lead_genes": "Lead_genes",
            },
            inplace=True,
        )
        # res_df["Gene %"] = res_df["Gene %"].map(lambda x: "{0:.2%}".format(x) if x !="" else "")
        # trim
        dc = ["RES", "hits", "matched_genes"]
        if self.permutation_num == 0:
            dc += [
                "NOM p-val",
                "FWER p-val",
                "FDR q-val",
                "Tag %",
                "Gene %",
                "Lead_genes",
            ]
        if self.module == "gsva":
            dc += ["NES"]
        # re-order by NES
        # for pandas > 1.1, use df.sort_values(by='B', key=abs) will sort by abs value
        self.res2d = res_df.reindex(
            res_df["NES"].abs().sort_values(ascending=False).index
        ).reset_index(drop=True)
        self.res2d.drop(dc, axis=1, inplace=True)

        if self._outdir is not None:
            out = os.path.join(
                self.outdir,
                "gseapy.{b}.{c}.report.csv".format(
                    b=self.permutation_type, c=self.module
                ),
            )
            self.res2d.to_csv(out, index=False, float_format="%.6e")
            with open(os.path.join(self.outdir, "gene_sets.gmt"), "w") as gout:
                for term, genes in gmt.items():
                    collection = ""
                    if term.find("__") > -1:
                        collections = term.split("__")
                        collection = collections[0]
                        term = collections[1]
                    gg = "\t".join(genes)
                    gout.write(f"{term}\t{collection}\t{gg}\n")
        # generate gseaplots
        if not self._noplot:
            self._plotting(metric)
        return

    @property
    def results(self):
        """
        compatible to old style
        """
        keys = list(self._results.keys())
        if len(keys) == 1:
            return self._results[keys[0]]
        return self._results

    def enrichment_score(
        self,
        gene_list: Iterable[str],
        correl_vector: Iterable[float],
        gene_set: Dict[str, List[str]],
        weight: float = 1.0,
        nperm: int = 1000,
        seed: int = 123,
        single: bool = False,
        scale: bool = False,
    ):
        """This is the most important function of GSEApy. It has the same algorithm with GSEA and ssGSEA.

        :param gene_list:       The ordered gene list gene_name_list, rank_metric.index.values
        :param gene_set:        gene_sets in gmt file, please use gmt_parser to get gene_set.
        :param weight:  It's the same with gsea's weighted_score method. Weighting by the correlation
                                is a very reasonable choice that allows significant gene sets with less than perfect coherence.
                                options: 0(classic),1,1.5,2. default:1. if one is interested in penalizing sets for lack of
                                coherence or to discover sets with any type of nonrandom distribution of tags, a value p < 1
                                might be appropriate. On the other hand, if one uses sets with large number of genes and only
                                a small subset of those is expected to be coherent, then one could consider using p > 1.
                                Our recommendation is to use p = 1 and use other settings only if you are very experienced
                                with the method and its behavior.

        :param correl_vector:   A vector with the correlations (e.g. signal to noise scores) corresponding to the genes in
                                the gene list. Or rankings, rank_metric.values
        :param nperm:           Only use this parameter when computing esnull for statistical testing. Set the esnull value
                                equal to the permutation number.
        :param seed:            Random state for initializing gene list shuffling. Default: seed=None

        :return:

        ES: Enrichment score (real number between -1 and +1)

        ESNULL: Enrichment score calculated from random permutations.

        Hits_Indices: Index of a gene in gene_list, if gene is included in gene_set.

        RES: Numerical vector containing the running enrichment score for all locations in the gene list .

        """
        N = len(gene_list)
        # Test whether each element of a 1-D array is also present in a second array
        # It's more intuitive here than original enrichment_score source code.
        # use .astype to covert bool to integer
        tag_indicator = np.in1d(gene_list, gene_set, assume_unique=True).astype(
            int
        )  # notice that the sign is 0 (no tag) or 1 (tag)

        if weight == 0:
            correl_vector = np.repeat(1, N)
        else:
            correl_vector = np.abs(correl_vector) ** weight

        # get indices of tag_indicator
        hit_ind = np.flatnonzero(tag_indicator).tolist()
        # if used for compute esnull, set esnull equal to permutation number, e.g. 1000
        # else just compute enrichment scores
        # set axis to 1, because we have 2D array
        axis = 1
        tag_indicator = np.tile(tag_indicator, (nperm + 1, 1))
        correl_vector = np.tile(correl_vector, (nperm + 1, 1))
        # gene list permutation
        rs = np.random.RandomState(seed)
        for i in range(nperm):
            rs.shuffle(tag_indicator[i])
        # np.apply_along_axis(rs.shuffle, 1, tag_indicator)

        Nhint = tag_indicator.sum(axis=axis, keepdims=True)
        sum_correl_tag = np.sum(correl_vector * tag_indicator, axis=axis, keepdims=True)
        # compute ES score, the code below is identical to gsea enrichment_score method.
        no_tag_indicator = 1 - tag_indicator
        Nmiss = N - Nhint
        norm_tag = 1.0 / sum_correl_tag
        norm_no_tag = 1.0 / Nmiss

        RES = np.cumsum(
            tag_indicator * correl_vector * norm_tag - no_tag_indicator * norm_no_tag,
            axis=axis,
        )

        if scale:
            RES = RES / N
        if single:
            es_vec = RES.sum(axis=axis)
        else:
            max_ES, min_ES = RES.max(axis=axis), RES.min(axis=axis)
            es_vec = np.where(np.abs(max_ES) > np.abs(min_ES), max_ES, min_ES)
        # extract values
        es, esnull, RES = es_vec[-1], es_vec[:-1], RES[-1, :]

        return es, esnull, hit_ind, RES

    def plot(
        self,
        terms: Union[str, List[str]],
        colors: Optional[Union[str, List[str]]] = None,
        legend_kws: Optional[Dict[str, Any]] = None,
        figsize: Tuple[float, float] = (4, 5),
        show_ranking: bool = True,
        ofname: Optional[str] = None,
    ):
        """
        terms: str, list.  terms/pathways to show
        colors: str, list. list of colors for each term/pathway
        legend_kws: kwargs to pass to ax.legend. e.g. `loc`, `bbox_to_achor`.
        ofname: savefig
        """
        # if hasattr(self, "results"):
        if self.module in ["ssgsea", "gsva"]:
            raise NotImplementedError("not for ssgsea")
        keys = list(self._results.keys())
        if len(keys) > 1:
            raise NotImplementedError("Multiple Dataset input No supported yet!")

        ranking = self.ranking if show_ranking else None
        if isinstance(terms, str):
            gsdict = self.results[terms]
            g = GSEAPlot(
                term=terms,
                tag=gsdict["hits"],
                rank_metric=ranking,
                runes=gsdict["RES"],
                nes=gsdict["nes"],
                pval=gsdict["pval"],
                fdr=gsdict["fdr"],
                ofname=ofname,
                pheno_pos=self.pheno_pos,
                pheno_neg=self.pheno_neg,
                color=colors,
                figsize=figsize,
            )
            g.add_axes()
            g.savefig()
            return g.fig

        elif hasattr(terms, "__len__"):  # means iterable
            terms = list(terms)
            tags = [self.results[t]["hits"] for t in terms]
            runes = [self.results[t]["RES"] for t in terms]
            t = TracePlot(
                terms=terms,
                tags=tags,
                runes=runes,
                rank_metric=ranking,
                colors=colors,
                legend_kws=legend_kws,
                ofname=ofname,
            )
            t.add_axes()
            t.savefig(ofname)
            return t.fig
        else:
            print("not supported input: terms")
