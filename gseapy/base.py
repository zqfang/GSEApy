#! python
# -*- coding: utf-8 -*-

import os
import logging
import json
from collections import OrderedDict
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
import requests

from gseapy.plot import gseaplot, heatmap
from gseapy.utils import mkdirs, log_init, retry, DEFAULT_LIBRARY, DEFAULT_CACHE_PATH
from typing import AnyStr, Tuple, Union, List, Dict, Iterable, Optional


class GSEAbase(object):
    """base class of GSEA."""

    def __init__(self, outdir: Optional[str] = None,
                 gene_sets: Union[List[str], str, Dict[str, str]] = 'KEGG_2016',
                 module: str = 'base',
                 threads: int = 1,
                 enrichr_url: str = "http://maayanlab.cloud",
                 verbose: bool = False):
        self.outdir = outdir
        self.gene_sets = gene_sets
        self.fdr = 0.05
        self.module = module
        self.results = None
        self.res2d = None
        self.ranking = None
        self.ascending = False
        self.verbose = verbose
        self._threads = threads
        self.ENRICHR_URL = enrichr_url
        self.pheno_pos = None
        self.pheno_neg = None

        self._set_cores()
        # init logger
        logfile = self.prepare_outdir()
        self._logger = log_init(outlog=logfile,
                                log_level=logging.INFO if self.verbose else logging.WARNING)

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

        # handle gmt type
        if isinstance(self.gene_sets, str):
            _gset = os.path.split(self.gene_sets)[-1].lower().rstrip(".gmt")
        elif isinstance(self.gene_sets, dict):
            _gset = "blank_name"
        else:
            raise Exception("Error parsing gene_sets parameter for gene sets")

        logfile = os.path.join(
            self.outdir, "gseapy.%s.%s.log" % (self.module, _gset))
        return logfile

    def _set_cores(self):
        """set cpu numbers to be used"""

        cpu_num = os.cpu_count()-1
        if self._threads > cpu_num:
            cores = cpu_num
        elif self._threads < 1:
            cores = 1
        else:
            cores = self._threads
        # have to be int if user input is float
        self._threads = int(cores)

    def _load_ranking(self, rnk: Union[pd.DataFrame, pd.Series, str]) -> pd.Series:
        """Parse ranking file. This file contains ranking correlation vector( or expression values)
           and gene names or ids.

            :param rnk: the .rnk file of GSEA input or a Pandas DataFrame, Series instance.
            :return: a Pandas Series with gene name indexed rankings

        """
        # load data
        if isinstance(rnk, pd.DataFrame):
            rank_metric = rnk.copy()
            # handle dataframe with gene_name as index.
            if rnk.shape[1] == 1:
                rank_metric = rnk.reset_index()
        elif isinstance(rnk, pd.Series):
            rank_metric = rnk.reset_index()
        elif os.path.isfile(rnk):
            rank_metric = pd.read_csv(rnk, header=None, sep="\t")
        else:
            raise Exception('Error parsing gene ranking values!')
        # sort ranking values from high to low
        rank_metric.sort_values(
            by=rank_metric.columns[1], ascending=self.ascending, inplace=True)
        # drop na values
        if rank_metric.isnull().any(axis=1).sum() > 0:
            self._logger.warning(
                "Input gene rankings contains NA values(gene name and ranking value), drop them all!")
            # print out NAs
            NAs = rank_metric[rank_metric.isnull().any(axis=1)]
            self._logger.debug('NAs list:\n'+NAs.to_string())
            rank_metric.dropna(how='any', inplace=True)
        # drop duplicate IDs, keep the first
        if rank_metric.duplicated(subset=rank_metric.columns[0]).sum() > 0:
            self._logger.warning(
                "Input gene rankings contains duplicated IDs, Only use the duplicated ID with highest value!")
            # print out duplicated IDs.
            dups = rank_metric[rank_metric.duplicated(
                subset=rank_metric.columns[0])]
            self._logger.debug('Dups list:\n'+dups.to_string())
            rank_metric.drop_duplicates(
                subset=rank_metric.columns[0], inplace=True, keep='first')
        # reset ranking index, because you have sort values and drop duplicates.
        rank_metric.reset_index(drop=True, inplace=True)
        rank_metric.columns = ['gene_name', 'rank']
        rankser = rank_metric.set_index('gene_name')['rank']
        self.ranking = rankser
        # return series
        return rankser

    def load_gmt_only(self, gmt: Union[List[str], str, Dict[str, str]]) -> Dict[str, List[str]]:
        if isinstance(gmt, dict):
            genesets_dict = gmt
        elif isinstance(gmt, str):
            genesets_dict = self.parse_gmt(gmt)
        else:
            raise Exception("Error parsing gmt parameter for gene sets")
        # self._gmtdct = genesets_dict
        return genesets_dict

    def load_gmt(self,
                 gene_list: Iterable[str],
                 gmt: Union[List[str], str, Dict[str, str]]) -> Dict[str, List[str]]:
        """load gene set dict"""

        if isinstance(gmt, dict):
            genesets_dict = gmt
        elif isinstance(gmt, str):
            genesets_dict = self.parse_gmt(gmt)
        else:
            raise Exception("Error parsing gmt parameter for gene sets")

        subsets = list(genesets_dict.keys())
        self.n_genesets = len(subsets)
        for subset in subsets:
            subset_list = genesets_dict.get(subset)
            if isinstance(subset_list, set):
                subset_list = list(subset_list)
                genesets_dict[subset] = subset_list
            tag_indicator = np.in1d(gene_list, subset_list, assume_unique=True)
            tag_len = tag_indicator.sum()
            if (self.min_size <= tag_len <= self.max_size) and tag_len < len(gene_list):
                # tag_len should not < gene_list
                continue
            del genesets_dict[subset]

        filsets_num = len(subsets) - len(genesets_dict)
        self._logger.info("%04d gene_sets have been filtered out when max_size=%s and min_size=%s" % (
            filsets_num, self.max_size, self.min_size))

        if filsets_num == len(subsets):
            self._logger.error("No gene sets passed through filtering condition!!! " +
                               "Try to set min_size or max_size parameters again!\n" +
                               "Note: check gene name, gmt format, or filtering size.")
            raise Exception("No gene sets passed through filtering condition")

        # self._gmtdct = genesets_dict
        return genesets_dict

    def parse_gmt(self, gmt: Union[List[str], str, Dict[str, str]]) -> Dict[str, List[str]]:
        """gmt parser"""

        if gmt.lower().endswith(".gmt"):
            with open(gmt) as genesets:
                genesets_dict = {line.strip().split("\t")[0]: line.strip().split("\t")[2:]
                                 for line in genesets.readlines()}
            return genesets_dict

        elif gmt in DEFAULT_LIBRARY:
            pass
        elif gmt in self.get_libraries():
            pass
        else:
            self._logger.error("No supported gene_sets: %s" % gmt)
            raise Exception("No supported gene_sets: %s" % gmt)

        tmpname = "enrichr." + gmt + ".gmt"
        tempath = os.path.join(DEFAULT_CACHE_PATH, tmpname)
        # if file already download
        if os.path.isfile(tempath):
            self._logger.info(
                "Enrichr library gene sets already downloaded in: %s, use local file" % DEFAULT_CACHE_PATH)
            return self.parse_gmt(tempath)
        else:
            return self._download_libraries(gmt)

    def get_libraries(self) -> List[str]:
        """return active enrichr library name.Offical API """

        lib_url = self.ENRICHR_URL+'/Enrichr/datasetStatistics'
        response = requests.get(lib_url, verify=True)
        if not response.ok:
            raise Exception("Error getting the Enrichr libraries")
        libs_json = json.loads(response.text)
        libs = [lib['libraryName'] for lib in libs_json['statistics']]

        return sorted(libs)

    def _download_libraries(self, libname: str) -> Dict[str, List[str]]:
        """ download enrichr libraries."""
        self._logger.info(
            "Downloading and generating Enrichr library gene sets......")
        s = retry(5)
        # queery string
        ENRICHR_URL = self.ENRICHR_URL+'/Enrichr/geneSetLibrary'
        query_string = '?mode=text&libraryName=%s'
        # get
        response = s.get(ENRICHR_URL + query_string % libname, timeout=None)
        if not response.ok:
            raise Exception(
                'Error fetching enrichment results, check internet connection first.')
        # reformat to dict and save to disk
        mkdirs(DEFAULT_CACHE_PATH)
        genesets_dict = {}
        outname = "enrichr.%s.gmt" % libname
        gmtout = open(os.path.join(DEFAULT_CACHE_PATH, outname), "w")
        for line in response.iter_lines(chunk_size=1024, decode_unicode='utf-8'):
            line = line.strip()
            k = line.split("\t")[0]
            v = list(map(lambda x: x.split(",")[0], line.split("\t")[2:]))
            genesets_dict.update({k: v})
            outline = "%s\t\t%s\n" % (k, "\t".join(v))
            gmtout.write(outline)
        gmtout.close()

        return genesets_dict

    def _heatmat(self, df: pd.DataFrame, classes: List[str]):
        """only use for gsea heatmap"""

        cls_booA = list(map(lambda x: True if x ==
                        self.pheno_pos else False, classes))
        cls_booB = list(map(lambda x: True if x ==
                        self.pheno_neg else False, classes))
        datA = df.loc[:, cls_booA]
        datB = df.loc[:, cls_booB]
        datAB = pd.concat([datA, datB], axis=1)
        self.heatmat = datAB
        return

    def _plotting(self, rank_metric: pd.Series):
        """ Plotting API.
            :param rank_metric: sorted pd.Series with rankings values.
        """

        # no values need to be returned
        if self._outdir is None:
            return
        # Plotting
        for i, record in self.res2d.iterrows():
            if self.module != 'ssgsea' and record['fdr'] > 0.05:
                continue
            if i >= self.graph_num:
                break
            hit = record['hits']
            NES = 'nes' if self.module != 'ssgsea' else 'es'
            term = record['term'].replace('/', '_').replace(":", "_")
            outfile = '{0}/{1}.{2}.{3}'.format(self.outdir,
                                               term, self.module, self.format)
            gseaplot(rank_metric=rank_metric, term=term, hit_indices=hit,
                     nes=record[NES], pval=record['pval'],
                     fdr=record['fdr'], RES=record['RES'],
                     pheno_pos=self.pheno_pos,
                     pheno_neg=self.pheno_neg,
                     figsize=self.figsize,
                     ofname=outfile)

            if self.module == 'gsea':
                outfile2 = "{0}/{1}.heatmap.{2}".format(
                    self.outdir, term, self.format)
                heatmat = self.heatmat.iloc[hit, :]
                width = np.clip(heatmat.shape[1], 4, 20)
                height = np.clip(heatmat.shape[0], 4, 20)
                heatmap(df=heatmat, title=term, ofname=outfile2,
                        z_score=0, figsize=(width, height),
                        xticklabels=True, yticklabels=True)

    def _to_df(self, gsea_summary,
               gmt: Dict[str, List[str]],
               rank_metric: pd.Series):
        """Convernt GSEASummary to DataFrame"""

        outcol = ['term', 'es', 'nes', 'pval', 'fdr', 'overlap',
                  'genes', 'lead_genes', 'hits', 'RES']
        res_df = []
        # res = OrderedDict()
        for gs in gsea_summary:
            # sample = '' if gs.name is None else gs.name
            # reformat gene list.
            _genes = rank_metric.index.values[gs.hits]
            genes = ";".join([str(g).strip() for g in _genes])

            RES = np.array(gs.run_es)
            # extract leading edge genes
            if float(gs.es) > 0:
                # RES -> ndarray, ind -> list
                ldg_pos = list(filter(lambda x: x <= RES.argmax(), gs.hits))
            elif float(gs.es) < 0:
                ldg_pos = list(filter(lambda x: x >= RES.argmin(), gs.hits))
            else:
                ldg_pos = gs.hits  # es == 0 ?
            lead_genes = ';'.join(
                list(map(str, rank_metric.iloc[ldg_pos].index)))
            overlap = "%s/%s" % (len(gs.hits), len(gmt[gs.term]))
            e = [gs.term, gs.es, gs.nes, gs.pval, gs.fdr,
                 overlap, genes, lead_genes, gs.hits, gs.run_es]
            res_df.append(e)

        # save to dataframe
        res_df = pd.DataFrame(res_df, columns=outcol)
        self.results = res_df.set_index('term').to_dict(orient='index')
        # save
        # res_df.set_index('term', inplace=True)
        res_df.drop(['RES', 'hits'], axis=1, inplace=True)
        res_df.sort_values(by=['fdr', 'nes'], inplace=True)
        self.res2d = res_df
        # self.results = res
        if self._outdir is None:
            return
        out = os.path.join(self.outdir, 'gseapy.{b}.{c}.report.csv'.format(
            b=self.permutation_type, c=self.module))
        res_df.to_csv(out, index=False, float_format='%.6e')

        return

    def enrichment_score(self,
                         gene_list: Iterable[str],
                         correl_vector: Iterable[float],
                         gene_set: Dict[str, List[str]],
                         weight: float = 1.0,
                         nperm: int = 1000,
                         seed: int = 123,
                         single: bool = False,
                         scale: bool = False):
        """This is the most important function of GSEApy. It has the same algorithm with GSEA and ssGSEA.

        :param gene_list:       The ordered gene list gene_name_list, rank_metric.index.values
        :param gene_set:        gene_sets in gmt file, please use gsea_gmt_parser to get gene_set.
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
            int)  # notice that the sign is 0 (no tag) or 1 (tag)

        if weight == 0:
            correl_vector = np.repeat(1, N)
        else:
            correl_vector = np.abs(correl_vector)**weight

        # get indices of tag_indicator
        hit_ind = np.flatnonzero(tag_indicator).tolist()
        # if used for compute esnull, set esnull equal to permutation number, e.g. 1000
        # else just compute enrichment scores
        # set axis to 1, because we have 2D array
        axis = 1
        tag_indicator = np.tile(tag_indicator, (nperm+1, 1))
        correl_vector = np.tile(correl_vector, (nperm+1, 1))
        # gene list permutation
        rs = np.random.RandomState(seed)
        for i in range(nperm):
            rs.shuffle(tag_indicator[i])
        # np.apply_along_axis(rs.shuffle, 1, tag_indicator)

        Nhint = tag_indicator.sum(axis=axis, keepdims=True)
        sum_correl_tag = np.sum(
            correl_vector*tag_indicator, axis=axis, keepdims=True)
        # compute ES score, the code below is identical to gsea enrichment_score method.
        no_tag_indicator = 1 - tag_indicator
        Nmiss = N - Nhint
        norm_tag = 1.0/sum_correl_tag
        norm_no_tag = 1.0/Nmiss

        RES = np.cumsum(tag_indicator * correl_vector *
                        norm_tag - no_tag_indicator * norm_no_tag, axis=axis)

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
