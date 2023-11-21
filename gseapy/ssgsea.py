#! python
# -*- coding: utf-8 -*-
import os
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from gseapy.base import GSEAbase
from gseapy.gse import CorrelType, ssgsea_rs
from gseapy.utils import mkdirs


class SingleSampleGSEA(GSEAbase):
    """GSEA extension: single sample GSEA"""

    def __init__(
        self,
        data: Union[pd.DataFrame, pd.Series, str],
        gene_sets: Union[List[str], str, Dict[str, str]],
        outdir: Optional[str] = None,
        sample_norm_method: Optional[str] = "rank",
        correl_norm_type: Optional[str] = None,
        min_size: int = 15,
        max_size: int = 500,
        permutation_num: Optional[int] = None,
        weight: float = 0.25,
        ascending: bool = False,
        threads: int = 1,
        figsize: Tuple[float, float] = (6.5, 6),
        format: str = "pdf",
        graph_num: int = 20,
        no_plot: bool = True,
        seed: int = 123,
        verbose: bool = False,
        **kwargs,
    ):
        super(SingleSampleGSEA, self).__init__(
            outdir=outdir,
            gene_sets=gene_sets,
            module="ssgsea",
            threads=threads,
            verbose=verbose,
        )
        self.data = data
        self.sample_norm_method = sample_norm_method
        self.correl_type = self.norm_correl(correl_norm_type)
        self.weight = weight
        self.min_size = min_size
        self.max_size = max_size
        self.permutation_num = int(permutation_num) if permutation_num else None
        self.ascending = ascending
        self.figsize = figsize
        self.format = format
        self.graph_num = int(graph_num)
        self.seed = seed
        self.ranking = None
        self._noplot = no_plot
        self.permutation_type = "gene_set"

    def corplot(self):
        """NES Correlation plot
        TODO
        """

    def setplot(self):
        """ranked genes' location plot
        TODO
        """

    def load_data(self) -> pd.DataFrame:
        # load data
        exprs = self._load_data(self.data)
        return self._check_data(exprs)

    def norm_samples(self, dat: pd.DataFrame) -> pd.DataFrame:
        """normalization samples
        see here: https://github.com/broadinstitute/ssGSEA2.0/blob/f682082f62ae34185421545f25041bae2c78c89b/src/ssGSEA2.0.R#L237
        """
        if self.sample_norm_method is None:
            self._logger.info("Use input as rank metric for ssGSEA")
            return dat

        if self.sample_norm_method == "rank":
            data = dat.rank(axis=0, method="average", na_option="bottom")
            data = 10000 * data / data.shape[0]
        elif self.sample_norm_method == "log_rank":
            data = dat.rank(axis=0, method="average", na_option="bottom")
            data = np.log(10000 * data / data.shape[0] + np.exp(1))
        elif self.sample_norm_method == "log":
            dat[dat < 1] = 1
            data = np.log(dat + np.exp(1))
        elif self.sample_norm_method == "custom":
            self._logger.info("Use custom rank metric for ssGSEA")
            data = dat
        else:
            raise Exception("No supported method: %s" % self.sample_norm_method)

        return data

    def norm_correl(self, cortype):
        """
        After norm_samples, input value will be further nomizalized before using them to calculate scores.
        see source code:
        https://github.com/broadinstitute/ssGSEA2.0/blob/f682082f62ae34185421545f25041bae2c78c89b/src/ssGSEA2.0.R#L396
        """
        if cortype is None:
            return CorrelType.Rank

        if cortype == "zscore":
            return CorrelType.ZScore
        elif cortype == "rank":
            return CorrelType.Rank
        elif cortype == "symrank":
            return CorrelType.SymRank
        else:
            raise Exception("unsupported correl type input")

    def run(self):
        """run entry"""
        if self.permutation_num is None:
            self.permutation_num = 0
        self._logger.info("Parsing data files for ssGSEA...........................")
        # load data
        data = self.load_data()
        # normalized samples, and rank
        normdat = self.norm_samples(data)
        self.ranking = normdat
        # filtering out gene sets and build gene sets dictionary
        gmt = self.load_gmt(gene_list=normdat.index.values, gmt=self.gene_sets)
        self.gmt = gmt
        self._logger.info(
            "%04d gene_sets used for further statistical testing....." % len(gmt)
        )
        # start analysis
        self._logger.info("Start to run ssGSEA...Might take a while................")
        if self.permutation_num > 0:
            # run permutation procedure and calculate pvals, fdrs
            self._logger.warning(
                "Run ssGSEA with permutation procedure, don't use the pval, fdr results for publication."
            )
        self.runSamplesPermu(df=normdat, gmt=gmt)

    def runSamplesPermu(
        self, df: pd.DataFrame, gmt: Optional[Dict[str, List[str]]] = None
    ):
        """Single Sample GSEA workflow with permutation procedure"""

        assert self.min_size <= self.max_size
        if self._outdir:
            mkdirs(self.outdir)
        gsum = ssgsea_rs(
            df.index.values.tolist(),
            df.values.tolist(),
            gmt,
            self.weight,
            self.min_size,
            self.max_size,
            self.permutation_num,  # permutate just like runing prerank analysis
            self.correl_type,
            self._threads,
            self.seed,
        )
        self.to_df(gsum.summaries, gmt, df, gsum.indices)
        return
