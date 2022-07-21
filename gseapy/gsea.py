#! python
# -*- coding: utf-8 -*-

import os
import logging
import glob
import numpy as np
import pandas as pd
from gseapy.base import GSEAbase
from gseapy.parser import gsea_edb_parser, gsea_cls_parser
from gseapy.plot import gseaplot
from gseapy.utils import mkdirs, log_init
from gseapy.gse import prerank_rs, gsea_rs, ssgsea_rs, Metric  # import gseapy rust lib
from typing import AnyStr, Tuple, Union, List, Dict, Iterable, Optional


class GSEA(GSEAbase):
    """GSEA main tool"""

    def __init__(self, data: Union[pd.DataFrame, str],
                 gene_sets: Union[List[str], str, Dict[str, str]],
                 classes: Union[List[str], str],
                 outdir: Optional[str] = None,
                 min_size: int = 15,
                 max_size: int = 500,
                 permutation_num: int = 1000,
                 weight: float = 1.0,
                 permutation_type: str = 'gene_set',
                 method: str = 'log2_ratio_of_classes',
                 ascending: bool = False,
                 threads: int = 1,
                 figsize: Tuple[float, float] = (6.5, 6),
                 format: str = 'pdf',
                 graph_num: int = 20,
                 no_plot: bool = False,
                 seed: int = 123,
                 verbose: bool = False):
        super(GSEA, self).__init__(outdir=outdir,
                                   gene_sets=gene_sets,
                                   module='gsea',
                                   threads=threads,
                                   verbose=verbose)
        self.data = data
        self.classes = classes
        self.permutation_type = permutation_type
        self.method = method
        self.min_size = min_size
        self.max_size = max_size
        self.permutation_num = int(permutation_num) if int(
            permutation_num) > 0 else 0
        self.weight = weight
        self.ascending = ascending
        self.figsize = figsize
        self.format = format
        self.graph_num = int(graph_num)
        self.seed = seed
        self.ranking = None
        self._noplot = no_plot
        self.pheno_pos = "pos"
        self.pheno_neg = "neg"

    def load_data(self, cls_vec: List) -> Tuple[pd.DataFrame, Dict]:
        """pre-processed the data frame.new filtering methods will be implement here.
        """
        # read data in
        if isinstance(self.data, pd.DataFrame):
            exprs = self.data.copy()
            # handle index is gene_names
            if exprs.index.dtype == 'O':
                exprs = exprs.reset_index()
        elif os.path.isfile(self.data):
            # GCT input format?
            if self.data.endswith("gct"):
                exprs = pd.read_csv(self.data, skiprows=2, sep="\t")
            else:
                exprs = pd.read_csv(self.data, comment='#', sep="\t")
        else:
            raise Exception('Error parsing gene expression DataFrame!')

        # drop duplicated gene names
        if exprs.iloc[:, 0].duplicated().sum() > 0:
            self._logger.warning(
                "Warning: dropping duplicated gene names, only keep the first values")
            # drop duplicate gene_names.
            exprs.drop_duplicates(subset=exprs.columns[0], inplace=True)
        if exprs.isnull().any().sum() > 0:
            self._logger.warning(
                "Warning: Input data contains NA, filled NA with 0")
            exprs.dropna(how='all', inplace=True)  # drop rows with all NAs
            exprs = exprs.fillna(0)
        # set gene name as index
        exprs.set_index(keys=exprs.columns[0], inplace=True)
        # select numberic columns
        df = exprs.select_dtypes(include=[np.number])

        # in case the description column is numeric
        if len(cls_vec) == (df.shape[1] - 1):
            df = df.iloc[:, 1:]
        # drop any genes which std ==0
        cls_dict = {k: v for k, v in zip(df.columns, cls_vec)}
        df_std = df.groupby(by=cls_dict, axis=1).std()
        df = df[~df_std.isin([0]).any(axis=1)]
        df = df + 1e-08  # we don't like zeros!!!

        return df, cls_dict

    def calculate_metric(self,
                         df: pd.DataFrame,
                         method: str,
                         pos: str,
                         neg: str,
                         classes: Dict,
                         ascending: bool) -> pd.Series:
        """The main function to rank an expression table. works for 2d array.

        :param df:      gene_expression DataFrame.
        :param method:  The method used to calculate a correlation or ranking. Default: 'log2_ratio_of_classes'.
                        Others methods are:

                        1. 'signal_to_noise' (s2n) or 'abs_signal_to_noise' (abs_s2n)

                            You must have at least three samples for each phenotype.
                            The more distinct the gene expression is in each phenotype,
                            the more the gene acts as a “class marker”.

                        2. 't_test'

                            Uses the difference of means scaled by the standard deviation and number of samples.
                            Note: You must have at least three samples for each phenotype to use this metric.
                            The larger the t-test ratio, the more distinct the gene expression is in each phenotype
                            and the more the gene acts as a “class marker.”

                        3. 'ratio_of_classes' (also referred to as fold change).

                            Uses the ratio of class means to calculate fold change for natural scale data.

                        4. 'diff_of_classes'

                            Uses the difference of class means to calculate fold change for natural scale data

                        5. 'log2_ratio_of_classes'

                            Uses the log2 ratio of class means to calculate fold change for natural scale data.
                            This is the recommended statistic for calculating fold change for log scale data.

        :param str pos: one of labels of phenotype's names.
        :param str neg: one of labels of phenotype's names.
        :param dict classes: column id to group mapping.
        :param bool ascending:  bool or list of bool. Sort ascending vs. descending.
        :return: returns a pd.Series of correlation to class of each variable. Gene_name is index, and value is rankings.

        visit here for more docs: http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html
        """

        # exclude any zero stds.
        df_mean = df.groupby(by=classes, axis=1).mean()
        df_std = df.groupby(by=classes, axis=1).std()
        n_pos = np.sum(classes == pos)
        n_neg = np.sum(classes == neg)
        if method in ['signal_to_noise', 's2n']:
            ser = (df_mean[pos] - df_mean[neg])/(df_std[pos] + df_std[neg])
        elif method in ['abs_signal_to_noise', 'abs_s2n']:
            ser = ((df_mean[pos] - df_mean[neg]) /
                   (df_std[pos] + df_std[neg])).abs()
        elif method == 't_test':
            ser = (df_mean[pos] - df_mean[neg]) / \
                np.sqrt(df_std[pos]**2/n_pos+df_std[neg]**2/n_neg)
        elif method == 'ratio_of_classes':
            ser = df_mean[pos] / df_mean[neg]
        elif method == 'diff_of_classes':
            ser = df_mean[pos] - df_mean[neg]
        elif method == 'log2_ratio_of_classes':
            ser = np.log2(df_mean[pos] / df_mean[neg])
        else:
            logging.error("Please provide correct method name!!!")
            raise LookupError("Input method: %s is not supported" % method)
        ser = ser.sort_values(ascending=ascending)

        return ser

    def run(self):
        """GSEA main procedure"""
        m = self.method.lower()
        if m in ['signal_to_noise', 's2n']:
            method = Metric.Signal2Noise
        elif m in ['s2n', 'abs_signal_to_noise', 'abs_s2n']:
            method = Metric.AbsSignal2Noise
        elif m == 't_test':
            method = Metric.Ttest
        elif m == 'ratio_of_classes':
            method = Metric.RatioOfClasses
        elif m == 'diff_of_classes':
            method = Metric.DiffOfClasses
        elif m == 'log2_ratio_of_classes':
            method = Metric.Log2RatioOfClasses
        else:
            raise Exception("Sorry, input method %s is not supported" % m)

        assert self.permutation_type in ["phenotype", "gene_set"]
        assert self.min_size <= self.max_size

        # Start Analysis
        self._logger.info(
            "Parsing data files for GSEA.............................")
        # phenotype labels parsing
        phenoPos, phenoNeg, cls_vector = gsea_cls_parser(self.classes)
        # select correct expression genes and values.
        dat, cls_dict = self.load_data(cls_vector)
        self.cls_dict = cls_dict
        self.pheno_pos = phenoPos
        self.pheno_neg = phenoNeg
        # data frame must have length > 1
        assert len(dat) > 1
        # ranking metrics calculation.
        dat2 = self.calculate_metric(df=dat, method=self.method, pos=phenoPos, neg=phenoNeg,
                                     classes=cls_dict, ascending=self.ascending)
        self.ranking = dat2
        # filtering out gene sets and build gene sets dictionary
        gmt = self.load_gmt_only(gmt=self.gene_sets)

        self._logger.info(
            "%04d gene_sets used for further statistical testing....." % len(gmt))
        self._logger.info(
            "Start to run GSEA...Might take a while..................")
        # cpu numbers
        # compute ES, NES, pval, FDR, RES
        if self.permutation_type == "gene_set":
            gsum = prerank_rs(dat2.index.to_list(),  # gene list
                              dat2.squeeze().to_list(),  # ranking values
                              gmt,  # must be a dict object
                              self.weight,
                              self.min_size,
                              self.max_size,
                              self.permutation_num,
                              self._threads,
                              self.seed,
                              )
        else:  # phenotype permutation
            group = list(map(lambda x: True if x ==
                         phenoPos else False, cls_vector))
            gsum = gsea_rs(dat.index.to_list(),
                           dat.values.tolist(),  # each row is gene values across samples
                           gmt,
                           group,
                           method,
                           self.weight,
                           self.min_size,
                           self.max_size,
                           self.permutation_num,
                           self._threads,
                           self.seed)

        if self._outdir is not None:
            self._logger.info(
                "Start to generate GSEApy reports and figures............")

        self._to_df(gsum, gmt, dat2)
        # reorder datarame for heatmap
        self._heatmat(df=dat.loc[dat2.index], classes=cls_vector)
        # Plotting
        if not self._noplot:
            self._plotting(rank_metric=dat2)

        self._logger.info(
            "Congratulations. GSEApy ran successfully.................\n")
        if self._outdir is None:
            self._tmpdir.cleanup()

        return





class Prerank(GSEAbase):
    """GSEA prerank tool"""

    def __init__(self, rnk: Union[pd.DataFrame, pd.Series, str], 
                 gene_sets: Union[List[str], str, Dict[str, str]],
                 outdir: Optional[str] = None,
                 pheno_pos='Pos', 
                 pheno_neg='Neg', 
                 min_size: int = 15,
                 max_size: int = 500,
                 permutation_num: int = 1000,
                 weight: float = 1.0,
                 ascending: bool = False,
                 threads: int = 1,
                 figsize: Tuple[float, float] = (6.5, 6),
                 format: str = 'pdf',
                 graph_num: int = 20,
                 no_plot: bool = False,
                 seed: int = 123,
                 verbose: bool = False):
        super(Prerank, self).__init__(outdir=outdir,
                                      gene_sets=gene_sets,
                                      module='prerank',
                                      threads=threads,
                                      verbose=verbose)
        self.rnk = rnk
        self.pheno_pos = pheno_pos
        self.pheno_neg = pheno_neg
        self.min_size = min_size
        self.max_size = max_size
        self.permutation_num = int(permutation_num) if int(
            permutation_num) > 0 else 0
        self.weight = weight
        self.ascending = ascending
        self.figsize = figsize
        self.format = format
        self.graph_num = int(graph_num)
        self.seed = seed
        self.ranking = None
        self._noplot = no_plot
        self.permutation_type = 'gene_set'

    def run(self):
        """GSEA prerank workflow"""

        assert self.min_size <= self.max_size

        # parsing rankings
        dat2 = self._load_ranking(self.rnk)
        assert len(dat2) > 1

        # Start Analysis
        self._logger.info(
            "Parsing data files for GSEA.............................")
        # filtering out gene sets and build gene sets dictionary
        gmt = self.load_gmt_only(gmt=self.gene_sets)

        self._logger.info(
            "%04d gene_sets used for further statistical testing....." % len(gmt))
        self._logger.info(
            "Start to run GSEA...Might take a while..................")
        # compute ES, NES, pval, FDR, RES
        gsum = prerank_rs(dat2.index.to_list(),  # gene list
                          dat2.squeeze().to_list(),  # ranking values
                          gmt,  # must be a dict object
                          self.weight,
                          self.min_size,
                          self.max_size,
                          self.permutation_num,
                          self._threads,
                          self.seed,
                          )

        if self._outdir is not None:
            self._logger.info(
                "Start to generate gseapy reports, and produce figures...")
        self._to_df(gsum, gmt, rank_metric=dat2)

        # Plotting
        if not self._noplot:
            self._plotting(rank_metric=dat2)

        self._logger.info(
            "Congratulations. GSEApy runs successfully................\n")
        if self._outdir is None:
            self._tmpdir.cleanup()

        return


class SingleSampleGSEA(GSEAbase):
    """GSEA extension: single sample GSEA"""

    def __init__(self, data, 
                 gene_sets: Union[List[str], str, Dict[str, str]],
                 outdir: Optional[str] = None,
                 sample_norm_method: str='rank', 
                 min_size: int = 15,
                 max_size: int = 500,
                 permutation_num: Optional[int]=None, 
                 weight: float =0.25,
                 ascending: bool = False,
                 threads: int = 1,
                 figsize: Tuple[float, float] = (6.5, 6),
                 format: str = 'pdf',
                 graph_num: int = 20,
                 no_plot: bool = True,
                 seed: int = 123,
                 verbose: bool = False):
        super(SingleSampleGSEA, self).__init__(outdir=outdir,
                                               gene_sets=gene_sets,
                                               module='ssgsea',
                                               threads=threads,
                                               verbose=verbose)
        self.data = data
        self.sample_norm_method = sample_norm_method
        self.weight = weight
        self.min_size = min_size
        self.max_size = max_size
        self.permutation_num = permutation_num if permutation_num is None else int(permutation_num)
        self.ascending = ascending
        self.figsize = figsize
        self.format = format
        self.graph_num = int(graph_num)
        self.seed = seed
        self.ranking = None
        self._noplot = no_plot
        self.permutation_type = 'gene_set'

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
        exprs = self.data
        if isinstance(exprs, pd.DataFrame):
            rank_metric = exprs.copy()
            # handle dataframe with gene_name as index.
            self._logger.debug("Input data is a DataFrame with gene names")
            # handle index is not gene_names
            if rank_metric.index.dtype != 'O':
                rank_metric.set_index(
                    keys=rank_metric.columns[0], inplace=True)

            rank_metric = rank_metric.select_dtypes(include=[np.number])
        elif isinstance(exprs, pd.Series):
            # change to DataFrame
            self._logger.debug("Input data is a Series with gene names")
            rank_metric = pd.DataFrame(exprs)
            # rename col if name attr is none
            rank_metric.columns = ["sample1"]
        elif os.path.isfile(exprs):
            # GCT input format?
            if exprs.endswith("gct"):
                rank_metric = pd.read_csv(
                    exprs, skiprows=1, comment='#', index_col=0, sep="\t")
            else:
                # just txt file like input
                rank_metric = pd.read_csv(
                    exprs, comment='#', index_col=0, sep="\t")
                if rank_metric.shape[1] == 1:
                    # rnk file like input
                    rank_metric = pd.read_csv(exprs, header=None, comment='#',
                                              names=['sample1'], index_col=0, sep="\t")
            # select numbers
            rank_metric = rank_metric.select_dtypes(include=[np.number])
        else:
            raise Exception('Error parsing gene ranking values!')

        if rank_metric.index.duplicated().sum() > 0:
            self._logger.warning(
                "Warning: dropping duplicated gene names, only keep the first values")
            rank_metric = rank_metric.loc[rank_metric.index.drop_duplicates(
                keep='first')]
            rank_metric = rank_metric.loc[rank_metric.index.dropna()]
        if rank_metric.isnull().any().sum() > 0:
            self._logger.warning(
                "Warning: Input data contains NA, filled NA with 0")
            rank_metric = rank_metric.fillna(0)

        return rank_metric

    def norm_samples(self, dat: pd.DataFrame) -> pd.DataFrame:
        """normalization samples
           see here: http://rowley.mit.edu/caw_web/ssGSEAProjection/ssGSEAProjection.Library.R
        """

        if self.sample_norm_method == 'rank':
            data = dat.rank(axis=0, method='average', na_option='bottom')
            data = 10000*data / data.shape[0]
        elif self.sample_norm_method == 'log_rank':
            data = dat.rank(axis=0, method='average', na_option='bottom')
            data = np.log(10000*data / data.shape[0] + np.exp(1))
        elif self.sample_norm_method == 'log':
            dat[dat < 1] = 1
            data = np.log(dat + np.exp(1))
        elif self.sample_norm_method == 'custom':
            self._logger.info("Use custom rank metric for ssGSEA")
            data = dat
        else:
            raise Exception("No supported method: %s" %
                            self.sample_norm_method)

        return data

    def run(self):
        """run entry"""
        self._logger.info(
            "Parsing data files for ssGSEA...........................")
        # load data
        data = self.load_data()
        # normalized samples, and rank
        normdat = self.norm_samples(data)
        # filtering out gene sets and build gene sets dictionary
        gmt = self.load_gmt(gene_list=normdat.index.values, gmt=self.gene_sets)
        self._logger.info(
            "%04d gene_sets used for further statistical testing....." % len(gmt))
        # start analysis
        self._logger.info(
            "Start to run ssGSEA...Might take a while................")
        if (self.permutation_num is not None) and self.permutation_num > 0:
            # run permutation procedure and calculate pvals, fdrs
            self._logger.warning(
                "run ssGSEA with permutation procedure, don't use the pval, fdr results for publication.")
        self.runSamplesPermu(df=normdat, gmt=gmt)
        # clean up all outputs if _outdir is None
        if self._outdir is None:
            self._tmpdir.cleanup()

    def runSamplesPermu(self, 
                        df: pd.DataFrame, 
                        gmt: Optional[Dict[str, List[str]]]= None):
        """Single Sample GSEA workflow with permutation procedure"""

        assert self.min_size <= self.max_size
        mkdirs(self.outdir)
        gsum = ssgsea_rs(df.index.to_list(),
                         df.values.tolist(),
                         gmt,
                         df.columns.astype(str).to_list(),  # sample name
                         self.weight,
                         self.min_size,
                         self.max_size,
                         self.permutation_num,  # permutate just like runing prerank analysis
                         self._threads,
                         self.seed
                         )

        outcol = ['name', 'term', 'es', 'nes', 'pval', 'fdr', 'fwer p-val', 'tag %', 'gene %',
                  'lead_genes', 'matched_genes','hits', 'RES']
        out = []
        for i, gs in enumerate(gsum):
            e = [gs.name, gs.term, gs.es, gs.nes, gs.pval,
                 gs.fdr, gs.fwerp, '', '', '', '', gs.hits, gs.run_es]
            out.append(e)
        out = pd.DataFrame(out, columns=outcol)
        # plotting
        idx_o = outcol.index('tag %')
        idx_g = outcol.index('gene %')
        idx_m = outcol.index('matched_genes')
        idx_l = outcol.index('lead_genes')
        for name, dd in out.groupby(['name']):
            rnk = df[name].sort_values(ascending=False, ignore_index=True)
            mkdirs(os.path.join(self.outdir, name))
            dd2 = dd.sort_values(by=['nes'])
            for i, gs in dd2.iterrows():
                RES = np.array(gs.RES)
                # extract leading edge genes
                if float(gs.es) >= 0:
                    es_i = RES.argmax()
                    # RES -> ndarray, ind -> list
                    ldg_pos = list(filter(lambda x: x <= es_i, gs.hits))
                    gene_frac = (es_i + 1) / len(rnk)
                else:
                    es_i = RES.argmin()
                    ldg_pos = list(filter(lambda x: x >= es_i, gs.hits))
                    ldg_pos.reverse()
                    gene_frac = (len(rnk) - es_i)/ len(rnk)

                lead_genes = ';'.join(
                    list(map(str, rnk.iloc[ldg_pos].index)))
                tag_frac = "%s/%s" % (len(ldg_pos), len(gmt[gs.term]))
                # reformat gene list.
                _genes = rnk.index.values[gs.hits]
                genes = ";".join([str(g).strip() for g in _genes])
                out.iloc[i, idx_o] = tag_frac
                out.iloc[i, idx_g] = gene_frac
                out.iloc[i, idx_m] = genes
                out.iloc[i, idx_l] = lead_genes
                if (not self._noplot) and i < self.graph_num:
                    term = gs.term.replace('/', '_').replace(":", "_")
                    outfile = "%s/%s/%s.%s.%s" % (self.outdir,
                                                  name, term, self.module, self.format)
                    gseaplot(rank_metric=rnk, term=term,
                             hit_indices=gs.hits, nes=gs.nes,
                             pval=gs.pval, fdr=gs.fdr,
                             RES=gs.RES, pheno_pos='', pheno_neg='',
                             figsize=self.figsize, ofname=outfile)

        out.drop(['RES', 'hits', 'matched_genes'], axis=1, inplace=True)
        if self.permutation_num < 1:
            out.drop(['pval', 'fdr'], axis=1, inplace=True)
        self.res2d = out
        if self._outdir is None:
            return
        outcsv = "gseapy.{a}ssgsea.report.csv".format(
            a="permut_genes." if self.permutation_num > 0 else '')
        out.to_csv(os.path.join(self.outdir, outcsv),
                   index=False, float_format='%.6e')
        return


class Replot(GSEAbase):
    """To reproduce GSEA desktop output results."""

    def __init__(self, 
                 indir: str, 
                 outdir: str='GSEApy_Replot', 
                 weight: float=1.0,
                 min_size: int=3, 
                 max_size: int=1000, 
                 figsize: Tuple[float, float]=(6.5, 6), 
                 format: str='pdf', 
                 verbose: bool=False):
        self.indir = indir
        self.outdir = outdir
        self.weight = weight
        self.min_size = min_size
        self.max_size = max_size
        self.figsize = figsize
        self.format = format
        self.verbose = bool(verbose)
        self.module = 'replot'
        self.gene_sets = None
        self.ascending = False
        # init logger
        mkdirs(self.outdir)
        outlog = os.path.join(self.outdir, "gseapy.%s.%s.log" %
                              (self.module, "run"))
        self._logger = log_init(outlog=outlog,
                                log_level=logging.INFO if self.verbose else logging.WARNING)

    def run(self):
        """main replot function"""
        assert self.min_size <= self.max_size

        # parsing files.......
        try:
            results_path = glob.glob(self.indir+'*/edb/results.edb')[0]
            rank_path = glob.glob(self.indir+'*/edb/*.rnk')[0]
            gene_set_path = glob.glob(self.indir+'*/edb/gene_sets.gmt')[0]
        except IndexError as e:
            raise Exception(
                "Could not locate GSEA files in the given directory!")
        # extract sample names from .cls file
        cls_path = glob.glob(self.indir+'*/edb/*.cls')
        if cls_path:
            pos, neg, classes = gsea_cls_parser(cls_path[0])
        else:
            # logic for prerank results
            pos, neg = '', ''
        # start reploting
        self.gene_sets = gene_set_path
        # obtain gene sets
        gene_set_dict = self.parse_gmt(gmt=gene_set_path)
        # obtain rank_metrics
        rank_metric = self._load_ranking(rank_path)
        correl_vector = rank_metric.values
        gene_list = rank_metric.index.values
        # extract each enriment term in the results.edb files and plot.

        database = gsea_edb_parser(results_path)
        for enrich_term, data in database.items():
            # extract statistical resutls from results.edb file
            hit_ind, nes, pval, fdr = data
            gene_set = gene_set_dict.get(enrich_term)
            # calculate enrichment score
            RES = self.enrichment_score(gene_list=gene_list,
                                        correl_vector=correl_vector,
                                        gene_set=gene_set,
                                        weight=self.weight,
                                        nperm=0)[-1]
            # plotting
            term = enrich_term.replace('/', '_').replace(":", "_")
            outfile = '{0}/{1}.{2}.{3}'.format(self.outdir,
                                               term, self.module, self.format)
            gseaplot(rank_metric=rank_metric, term=enrich_term,
                     hits=hit_ind, nes=nes, pval=pval, fdr=fdr,
                     RES=RES, pheno_pos=pos, pheno_neg=neg,
                     figsize=self.figsize, ofname=outfile)

        self._logger.info(
            "Congratulations! Your plots have been reproduced successfully!\n")
