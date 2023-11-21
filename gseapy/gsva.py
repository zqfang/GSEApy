#! python
# -*- coding: utf-8 -*-

import os
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd

from gseapy.base import GSEAbase
from gseapy.gse import gsva_rs
from gseapy.utils import mkdirs


class GSVA(GSEAbase):
    """GSVA"""

    def __init__(
        self,
        data: Union[pd.DataFrame, pd.Series, str],
        gene_sets: Union[List[str], str, Dict[str, str]],
        outdir: Optional[str] = None,
        kcdf: Optional[str] = "Gaussian",
        weight: float = 1.0,
        mx_diff: bool = True,
        abs_rnk: bool = False,
        min_size: int = 15,
        max_size: int = 500,
        threads: int = 1,
        seed: int = 123,
        verbose: bool = False,
        **kwargs,
    ):
        super(GSVA, self).__init__(
            outdir=outdir,
            gene_sets=gene_sets,
            module="gsva",
            threads=threads,
            verbose=verbose,
        )
        self.data = data
        self.tau = weight
        self.min_size = min_size
        self.max_size = max_size
        self.seed = seed
        self.mx_diff = mx_diff
        self.abs_rnk = abs_rnk
        self.ranking = None
        self.permutation_num = 0
        self._noplot = True
        if kcdf in ["Gaussian", "gaussian"]:
            self.kernel = True
            self.rnaseq = False
        elif kcdf in ["Poisson", "poisson"]:
            self.kernel = True
            self.rnaseq = True
        else:
            self.kernel = False
            self.rnaseq = False

        # self.figsize = figsize
        # self.format = format
        # self.graph_num = int(graph_num)
        # self.seed = seed
        self.ranking = None
        self.permutation_type = "gene_set"

    def load_data(self) -> pd.DataFrame:
        # load data
        data = self._load_data(self.data)
        data = self._check_data(data)
        return data

    def run(self):
        """run entry"""
        assert self.min_size <= self.max_size
        if self._outdir:
            mkdirs(self.outdir)

        self._logger.info("Parsing data files for GSVA.............................")
        # load data
        df = self.load_data()
        # kernel
        if self.kernel:
            if self.rnaseq:
                self._logger.info(
                    "Estimating ECDFs with Poisson kernels. Clip negative values to 0 !"
                )
                df = df.clip(lower=0)
            else:
                self._logger.info("Estimating ECDFs with Gaussian kernels.")
        else:
            self._logger.info("Estimating ECDFs with directly.")
        # save data
        self.data = df
        # normalized samples, and rank
        # filtering out gene sets and build gene sets dictionary
        gmt = self.load_gmt(gene_list=df.index.values, gmt=self.gene_sets)
        self.gmt = gmt
        self._logger.info(
            "%04d gene_sets used for further statistical testing....." % len(gmt)
        )
        # start analysis
        self._logger.info("Start to run GSVA...Might take a while................")
        # run
        gsum = gsva_rs(
            df.index.values.tolist(),
            df.values.tolist(),
            gmt,
            self.kernel,
            self.rnaseq,
            self.mx_diff,
            self.abs_rnk,
            self.tau,
            self.min_size,
            self.max_size,
            self._threads,
        )
        self.to_df(gsum.summaries, gmt, df, gsum.indices)
        self.ranking = gsum.rankings
        self._logger.info("Done")
