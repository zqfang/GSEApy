#! python
# -*- coding: utf-8 -*-

import os,sys,logging,json
import requests
import numpy as np
import pandas as pd
from collections import OrderedDict
from multiprocessing import Pool, cpu_count
from numpy import log, exp
from gseapy.parser import *
from gseapy.algorithm import enrichment_score, gsea_compute, gsea_compute_ss, ranking_metric
from gseapy.plot import gsea_plot, heatmap
from gseapy.utils import mkdirs, DEFAULT_LIBRARY


class GSEAbase(object):
    """base class of GSEA."""
    def __init__(self):
        self.gene_sets='KEGG_2016'
        self.module='base'
        self.results=None
        self.res2d=None
        self.ranking=None
        self.ascending=False
        self.verbose=False
        self._processes=1
        self._logger=None

    def _log_init(self, module='GSEA', log_level=logging.INFO):
        """logging start"""

        # clear old root logger handlers
        logging.getLogger("").handlers = []
        # log file naming rules
        gene_set =os.path.split(self.gene_sets)[-1].lower().rstrip(".gmt")
        # init a root logger
        logging.basicConfig(level    = logging.DEBUG,
                            format   = 'LINE %(lineno)-4d: %(asctime)s [%(levelname)-8s] %(message)s',
                            filename = "%s/gseapy.%s.%s.log"%(self.outdir, module, gene_set),
                            filemode = 'w')

        # define a Handler which writes INFO messages or higher to the sys.stderr
        console = logging.StreamHandler()
        console.setLevel(log_level)
        # set a format which is simpler for console use
        formatter = logging.Formatter('%(asctime)s %(message)s')
        # tell the handler to use this format
        console.setFormatter(formatter)
        # add handlers
        logging.getLogger("").addHandler(console)
        logger = logging.getLogger("")
        #logger.setLevel(log_level)
        self._logger=logger
        return logger

    def _log_stop(self):
        """log stop"""

        handlers = self._logger.handlers[:]
        for handler in handlers:
            handler.close()
            self._logger.removeHandler(handler)

    def _set_cores(self):
        """set cpu numbers to be used"""

        cpu_num = cpu_count()-1
        if self._processes > cpu_num:
            cores = cpu_num
        elif self._processes < 1:
            cores = 1
        else:
            cores = self._processes
        # have to be int if user input is float
        self._processes = int(cores)

    def _load_ranking(self, rnk):
        """Parse ranking file. This file contains ranking correlation vector( or expression values)
           and gene names or ids.

        :param rnk: the .rnk file of GSEA input or a pandas DataFrame, Series instance.
        :return: a pandas DataFrame with 3 columns names are: 'gene_name','rank',rank2'

        """
        #load data
        if isinstance(rnk, pd.DataFrame):
            rank_metric = rnk.copy()
            # handle dataframe with gene_name as index.
            if rnk.shape[1] == 1:
                rank_metric = rnk.reset_index()
        elif isinstance(rnk, pd.Series):
            rank_metric = rnk.reset_index()
        elif os.path.isfile(rnk):
            rank_metric = pd.read_table(rnk, header=None, comment='#')
        else:
            raise Exception('Error parsing gene ranking values!')

        #sort ranking values from high to low
        if rank_metric.shape[1] == 2: # if cols > 2, means this is a expression matrixs of gct format
            rank_metric.sort_values(by=rank_metric.columns[1], ascending=self.ascending, inplace=True)
        #drop na values
        if rank_metric.isnull().any(axis=1).sum() >0:
            self._logger.warning("Input gene rankings contains NA values(gene name and ranking value), drop them all!")
            #print out NAs
            NAs = rank_metric[rank_metric.isnull().any(axis=1)]
            self._logger.debug('NAs list:\n'+NAs.to_string())
            rank_metric.dropna(how='any', inplace=True)
        #drop duplicate IDs, keep the first
        if rank_metric.duplicated(subset=rank_metric.columns[0]).sum() >0:
            self._logger.warning("Input gene rankings contains duplicated IDs, Only use the duplicated ID with highest value!")
            #print out duplicated IDs.
            dups = rank_metric[rank_metric.duplicated(subset=rank_metric.columns[0])]
            self._logger.debug('Dups list:\n'+dups.to_string())
            rank_metric.drop_duplicates(subset=rank_metric.columns[0], inplace=True, keep='first')


        if rank_metric.shape[1] == 2:
            #reset ranking index, because you have sort values and drop duplicates.
            rank_metric.reset_index(drop=True, inplace=True)
            rank_metric.columns = ['gene_name','rank']
            self.ranking = rank_metric['rank']
            #use for plotting, need 2d array
            rank_metric['rank2'] = rank_metric['rank']
        return rank_metric

    def parse_gmt(self, min_size = 1, max_size = 1000):
        """gmt parser"""

        gmt = self.gene_set
        if gmt.lower().endswith(".gmt"):
            logging.info("User Defined gene sets is given.......continue..........")
            with open(gmt) as genesets:
                 genesets_dict = { line.strip("\n").split("\t")[0]: line.strip("\n").split("\t")[2:]
                                  for line in genesets.readlines()}
        else:
            logging.info("Downloading and generating Enrichr library gene sets...")
            if gmt in DEFAULT_LIBRARY:
                names = DEFAULT_LIBRARY
            else:
                names = self.get_libraries()

            if gmt in names:
                """
                define max tries num
                if the backoff_factor is 0.1, then sleep() will sleep for
                [0.1s, 0.2s, 0.4s, ...] between retries.
                It will also force a retry if the status code returned is 500, 502, 503 or 504.
                """
                s = requests.Session()
                retries = Retry(total=5, backoff_factor=0.1,
                                status_forcelist=[ 500, 502, 503, 504 ])
                s.mount('http://', HTTPAdapter(max_retries=retries))
                #queery string
                ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary'
                query_string = '?mode=text&libraryName=%s'
                #get
                response = s.get( ENRICHR_URL + query_string % gmt, timeout=None)
            else:
                raise Exception("gene_set files(.gmt) not found")

            if not response.ok:
                raise Exception('Error fetching enrichment results, check internet connection first.')
            ### reformat to dict and wirte to disk
            genesets_dict = {}
            with open(os.path.join(self.outdir,"gene_sets.gmt"),"w") as g:
                for line in response.iter_lines(chunk_size=1024, decode_unicode='utf-8'):
                    g.write(line)
                    genesets_dict.update({ line.split("\t")[0]:
                                           list(map(lambda x: x.split(",")[0], line.split("\t")[2:-1]))})


        return genesets_dict

    def tag_indi_mat(self, gsets_dict, cor_keys, min_size=1, max_size=1000):
        """build geneset to gene_keys matrix """

        keys = sorted(gsets_dict.keys())
        tag_indicator = np.vstack([np.in1d(cor_keys, gsets_dict[key], assume_unique=True) for key in keys])
        #filtering gene sets
        Nhits = tag_indicator.sum(axis=1)
        # python keyword ``and`` is not work for ndarray
        mask = (min_size <= Nhits) & (Nhits <= max_size)
        tag_mat = tag_indicator[mask]
        filsets_num = tag_indicator.shape[0] - tag_mat.shape[0]
        logging.info("%04d gene_sets have been filtered out when max_size=%s and min_size=%s"%(filsets_num, max_size, min_size))

        if filsets_num == len(keys):
            logging.error("No gene sets passed throught filtering condition!!!, try new paramters again!\n" +\
                             "Note: Gene names for GSEApy is case sensitive." )
            sys.exit(1)

        return tag_mat

    def _plotting(self, rank_metric, results, res2d,
                 graph_num, outdir, format, figsize, module=None, data=None,
                 classes=None, phenoPos='', phenoNeg=''):
        """
        :param rank_metric: dat2.
        :param results: self.results
        :param data: preprocessed expression table

        """
        #Plotting
        top_term = res2d.head(graph_num).index

        #multi-threading
        pool = Pool(processes=self._processes)

        for gs in top_term:
            hit = results.get(gs)['hit_index']
            NES = 'nes' if module != 'ssgsea' else 'es'
            """
            gsea_plot(rank_metric=rank_metric, enrich_term=gs, hit_ind=hit,
                      nes=results.get(gs)['nes'], pval=results.get(gs)['pval'], fdr=results.get(gs)['fdr'],
                      RES=results.get(gs)['rank_ES'], phenoPos=phenoPos, phenoNeg=phenoNeg, figsize=figsize,
                      format=format, outdir=outdir, module=module)
            """
            pool.apply_async(gsea_plot, args=(rank_metric, gs, hit, results.get(gs)[NES],
                                              results.get(gs)['pval'],results.get(gs)['fdr'],
                                              results.get(gs)['rank_ES'],
                                              phenoPos, phenoNeg, figsize, self.format,
                                              self.outdir,self.module))
        pool.close()
        pool.join()

        if module == 'gsea':
            width = len(classes) if len(classes) >= 6 else  5
            cls_booA =list(map(lambda x: True if x == phenoPos else False, classes))
            cls_booB =list(map(lambda x: True if x == phenoNeg else False, classes))
            datA = data.loc[:, cls_booA]
            datB = data.loc[:, cls_booB]
            datAB=pd.concat([datA,datB], axis=1)
            pool_heat = Pool(self._processes)

            #no values need to be returned
            for gs in top_term:
                hit = results.get(gs)['hit_index']
                pool_heat.apply_async(heatmap, args=(datAB.iloc[hit], gs, outdir, 0,
                                                    (width, len(hit)/2), format))
                #heatmap(datAB.iloc[hit], gs, outdir, 0, (width, len(hit)/2), format)

            pool_heat.close()
            pool_heat.join()

    def _save_results(self, zipdata, outdir, module, gmt, rank_metric, permutation_type):
        """reformat gsea results, and save to txt"""

        res = OrderedDict()
        for gs,gseale,ind,RES in zipdata:
            rdict = OrderedDict()
            rdict['es'] = gseale[0]
            rdict['nes'] = gseale[1]
            rdict['pval'] = gseale[2]
            rdict['fdr'] = gseale[3]
            rdict['gene_set_size'] = len(gmt[gs])
            rdict['matched_size'] = len(ind)
            #reformat gene list.
            _genes = rank_metric.iloc[ind, rank_metric.columns.get_loc('gene_name')]
            _genes = _genes.to_string(header=False, index=False).split("\n")
            rdict['genes'] = ",".join([ g.strip() for g in _genes ])

            rdict['rank_ES'] = RES
            rdict['hit_index'] = ind
            res[gs] = rdict

        res_df = pd.DataFrame.from_dict(res, orient='index')
        res_df.index.name = 'Term'
        res_df.sort_values(by='fdr', inplace=True)
        res_df.drop(['rank_ES','hit_index'], axis=1, inplace=True)

        out = '{a}/gseapy.{b}.{c}.report.csv'.format(a=outdir, b=module, c=permutation_type)
        if self.module == 'ssgsea':
            with open(out, 'a') as f:
                f.write('# normalize enrichment scores by random permutation procddure\n')
                f.write("# Same method to the orignial GSEA method, and it's not proper to use these values in your publication\n")
                res_df.to_csv(f, sep='\t')
        else:
            res_df.to_csv(out)

        self.res2d = res_df
        self.results  = res
        return

    def get_libraries(self):
        """return enrichr active enrichr library name.Offical API """
        lib_url='http://amp.pharm.mssm.edu/Enrichr/datasetStatistics'
        libs_json = json.loads(requests.get(lib_url).text)
        libs = [lib['libraryName'] for lib in libs_json['statistics']]
        return sorted(libs)


class GSEA(GSEAbase):
    """GSEA main tool"""
    def __init__(self, data, gene_sets, classes, outdir='GSEA_ouput',
                 min_size=15, max_size=500, permutation_num=1000,
                 weighted_score_type=1, permutation_type='gene_set',
                 method='log2_ratio_of_classes', ascending=False,
                 processes=1, figsize=[6.5,6], format='pdf', graph_num=20,
                 seed=None, verbose=False):

        self.data = data
        self.gene_sets=gene_sets
        self.classes=classes
        self.outdir=outdir
        self.permutation_type=permutation_type
        self.method=method
        self.min_size=min_size
        self.max_size=max_size
        self.permutation_num=permutation_num
        self.weighted_score_type=weighted_score_type
        self.ascending=ascending
        self._processes=processes
        self.figsize=figsize
        self.format=format
        self.graph_num=graph_num
        self.seed=seed
        self.verbose=verbose
        self.module='gsea'
        self.ranking=None

    def load_data(self, cls_vec):
        """pre-processed the data frame.new filtering methods will be implement here.
        """
        # read data in
        if isinstance(self.data, pd.DataFrame) :
            exprs = self.data.copy()
        elif os.path.isfile(self.data) :
            # GCT input format?
            if self.data.endswith("gct"):
                exprs = pd.read_table(self.data, skiprows=1, comment='#')
            else:
                exprs = pd.read_table(self.data, comment='#')
        else:
            raise Exception('Error parsing gene expression dataframe!')
            sys.exit(1)

        #drop duplicated gene names
        if exprs.iloc[:,0].duplicated().sum() > 0:
            self._logger.warning("Warning: dropping duplicated gene names, only keep the first values")
            exprs.drop_duplicates(subset=exprs.columns[0], inplace=True) #drop duplicate gene_names.
        if exprs.isnull().any().sum() > 0:
            self._logger.warning("Warning: Input data contains NA, filled NA with 0")
            exprs.dropna(how='all', inplace=True) #drop rows with all NAs
            exprs = exprs.fillna(0)
        # set gene name as index
        exprs.set_index(keys=exprs.columns[0], inplace=True)
        # select numberic columns
        df = exprs.select_dtypes(include=[np.number])
        # drop any genes which std ==0
        df_std =  df.groupby(by=cls_vec, axis=1).std()
        df =  df[~df_std.isin([0]).any(axis=1)]
        df = df + 0.00001 # we don't like zeros!!!
        return df

    def run(self):
        """GSEA main procedure"""

        assert self.permutation_type in ["phenotype", "gene_set"]
        assert self.min_size <= self.max_size
        # creat output dirs
        mkdirs(self.outdir)
        logger = self._log_init(module=self.module,
                               log_level=logging.INFO if self.verbose else logging.WARNING)
        #Start Analysis
        logger.info("Parsing data files for GSEA.............................")
        # phenotype labels parsing
        phenoPos, phenoNeg, cls_vector = gsea_cls_parser(self.classes)
        #select correct expression genes and values.
        dat = self.load_data(cls_vector)
        #data frame must have lenght > 1
        assert len(dat) > 1
        #ranking metrics calculation.
        dat2 = ranking_metric(df=dat, method=self.method, phenoPos=phenoPos, phenoNeg=phenoNeg,
                              classes=cls_vector, ascending=self.ascending)
        #filtering out gene sets and build gene sets dictionary
        gmt = gsea_gmt_parser(self.gene_sets, min_size=self.min_size, max_size=self.max_size,
                              gene_list=dat2['gene_name'].values)

        logger.info("%04d gene_sets used for further statistical testing....."% len(gmt))
        logger.info("Start to run GSEA...Might take a while..................")
        #cpu numbers
        self._set_cores()
        #compute ES, NES, pval, FDR, RES
        dataset = dat if self.permutation_type =='phenotype' else dat2
        gsea_results,hit_ind,rank_ES, subsets = gsea_compute(data=dataset, n=self.permutation_num, gmt=gmt,
                                                             weighted_score_type=self.weighted_score_type,
                                                             permutation_type=self.permutation_type,
                                                             method=self.method,
                                                             phenoPos=phenoPos, phenoNeg=phenoNeg,
                                                             classes=cls_vector, ascending=self.ascending,
                                                             seed=self.seed, processes=self._processes)

        logger.info("Start to generate gseapy reports, and produce figures...")
        res_zip = zip(subsets, list(gsea_results), hit_ind, rank_ES)
        self._save_results(zipdata=res_zip, outdir=self.outdir, module=self.module,
                                   gmt=gmt, rank_metric=dat2, permutation_type="gene_sets")

        #Plotting
        heat_dat = dat.loc[dat2.gene_name]
        self._plotting(rank_metric=dat2, results=self.results, res2d=self.res2d,
                       graph_num=self.graph_num, outdir=self.outdir,
                       figsize=self.figsize, format=self.format, module=self.module,
                       data=heat_dat, classes=cls_vector, phenoPos=phenoPos, phenoNeg=phenoNeg)

        logger.info("Congratulations. GSEApy run successfully................\n")

        return


class Prerank(GSEAbase):
    """GSEA prerank tool"""
    def __init__(self, rnk, gene_sets, outdir='GSEA_prerank',
                 pheno_pos='Pos', pheno_neg='Neg', min_size=15, max_size=500,
                 permutation_num=1000, weighted_score_type=1,
                 ascending=False, processes=1, figsize=[6.5,6], format='pdf',
                 graph_num=20, seed=None, verbose=False):

        self.rnk =rnk
        self.gene_sets=gene_sets
        self.outdir=outdir
        self.pheno_pos=pheno_pos
        self.pheno_neg=pheno_neg
        self.min_size=min_size
        self.max_size=max_size
        self.permutation_num=permutation_num
        self.weighted_score_type=weighted_score_type
        self.ascending=ascending
        self.figsize=figsize
        self.format=format
        self.graph_num=graph_num
        self.seed=seed
        self.verbose=verbose
        self.ranking=None
        self.module='prerank'
        self._processes=processes


    def run(self):
        """GSEA prerank workflow"""

        assert self.min_size <= self.max_size
        mkdirs(self.outdir)
        logger = self._log_init(module=self.module,
                               log_level=logging.INFO if self.verbose else logging.WARNING)
        #parsing rankings
        dat = self._load_ranking(self.rnk)
        dat2 = dat.set_index("gene_name")
        assert len(dat2) > 1

        #cpu numbers
        self._set_cores()
        #Start Analysis
        logger.info("Parsing data files for GSEA.............................")
        #filtering out gene sets and build gene sets dictionary
        gmt = gsea_gmt_parser(self.gene_sets, min_size=self.min_size, max_size=self.max_size,
                              gene_list=dat['gene_name'].values)
        logger.info("%04d gene_sets used for further statistical testing....."% len(gmt))
        logger.info("Start to run GSEA...Might take a while..................")
        #compute ES, NES, pval, FDR, RES
        gsea_results, hit_ind,rank_ES, subsets = gsea_compute(data=dat2, n=self.permutation_num, gmt=gmt,
                                                              weighted_score_type=self.weighted_score_type,
                                                              permutation_type='gene_set', method=None,
                                                              phenoPos=self.pheno_pos, phenoNeg=self.pheno_neg,
                                                              classes=None, ascending=self.ascending, seed=self.seed,
                                                              processes=self._processes)

        logger.info("Start to generate gseapy reports, and produce figures...")
        res_zip = zip(subsets, list(gsea_results), hit_ind, rank_ES)
        self._save_results(zipdata=res_zip, outdir=self.outdir, module=self.module,
                                   gmt=gmt, rank_metric=dat, permutation_type="gene_sets")

        #Plotting
        self._plotting(rank_metric=dat, results=self.results, res2d=self.res2d,
                       graph_num=self.graph_num, outdir=self.outdir,
                       figsize=self.figsize, format=self.format,
                       module=self.module, phenoPos=self.pheno_pos, phenoNeg=self.pheno_neg)

        logger.info("Congratulations. GSEApy run successfully................\n")

        return


class SingleSampleGSEA(GSEAbase):
    """GSEA extention: single sample GSEA"""
    def __init__(self, data, gene_sets, outdir="GSEA_SingleSample", sample_norm_method='rank',
                 min_size=15, max_size=2000, permutation_num=1000, weighted_score_type=0.25,
                 scale=True, ascending=False, processes=1, figsize=[7,6], format='pdf',
                 graph_num=20, seed=None, verbose=False):
        self.data=data
        self.gene_sets=gene_sets
        self.outdir=outdir
        self.sample_norm_method=sample_norm_method
        self.weighted_score_type=weighted_score_type
        self.scale = scale
        self.min_size=min_size
        self.max_size=max_size
        self.permutation_num=permutation_num
        self.ascending=ascending
        self.figsize=figsize
        self.format=format
        self.graph_num=graph_num
        self.seed=seed
        self.verbose=verbose
        self.ranking=None
        self.module='ssgsea'
        self._processes=processes
        self._imat=None

    def corplot(self):
        """NES Correlation plot
        """
    def setplot(self):
        """ranked genes' location plot
        """

    def load_data(self):
        #load data
        exprs = self.data
        if isinstance(exprs, pd.DataFrame):
            rank_metric = exprs.copy()
            # handle dataframe with gene_name as index.
            self._logger.debug("Input data is a DataFrame with gene names")
            #handle index is not gene_names
            if rank_metric.index.dtype != 'O':
                rank_metric.set_index(keys=rank_metric.columns[0], inplace=True)

            rank_metric = rank_metric.select_dtypes(include=[np.number])
        elif isinstance(exprs, pd.Series):
            #change to DataFrame
            self._logger.debug("Input data is a Series with gene names")
            rank_metric = pd.DataFrame(exprs)
        elif os.path.isfile(exprs):
            # GCT input format?
            if exprs.endswith("gct"):
                rank_metric = pd.read_table(exprs, skiprows=1, comment='#', index_col=0)
            else:
                #just rnk file input
                rank_metric = pd.read_table(exprs, header=None, comment='#', index_col=0)
                if rank_metric.shape[1] >=2:
                    # txt file input
                    rank_metric = pd.read_table(exprs, comment='#', index_col=0)
            #select numbers
            rank_metric = rank_metric.select_dtypes(include=[np.number])
        else:
            raise Exception('Error parsing gene ranking values!')

        if rank_metric.index.duplicated().sum() > 0:
            self._logger.warning("Warning: dropping duplicated gene names, only keep the first values")
            rank_metric = rank_metric.loc[rank_metric.index.drop_duplicates(keep='first')]
            rank_metric = rank_metric.loc[rank_metric.index.dropna()]
        if rank_metric.isnull().any().sum() > 0:
            self._logger.warning("Warning: Input data contains NA, filled NA with 0")
            rank_metric = rank_metric.fillna(0)
        # if single sample input, set ranking is not None for temp.
        if rank_metric.shape[1] <= 2: self.ranking=1

        return rank_metric

    def norm_samples(self, dat):
        """normalizatin samples
           see here: http://rowley.mit.edu/caw_web/ssGSEAProjection/ssGSEAProjection.Library.R
        """

        #set index of gene_names
        #data = self.data.set_index(keys=data.columns[0], inplace=True)

        if self.sample_norm_method == 'rank':
            data = dat.rank(axis=0, method='average', na_option='bottom')
            data = 10000*data / data.shape[0]
        elif self.sample_norm_method == 'log_rank':
            data = dat.rank(axis=0, method='average', na_option='bottom')
            data = log(10000*data / data.shape[0] + exp(1))
        elif self.sample_norm_method == 'log':
            dat[dat < 1] = 1
            data = log(dat + exp(1))
        elif self.sample_norm_method == 'custom':
            self._logger.info("Set user defined rank metric for ssGSEA")
        else:
            sys.stderr.write("No supported method: %s"%self.sample_norm_type)
            sys.exit(0)

        return data

    def run(self):
        """
        """
        mkdirs(self.outdir)
        logger = self._log_init(module=self.module,
                                log_level=logging.INFO if self.verbose else logging.WARNING)
        #load data
        data = self.load_data()

        # normalized samples, and rank
        normdat = self.norm_samples(data)

        # logic to process gct expression matrix
        if self.ranking is None:
            #gct expression matrix support for ssGSEA
            self.runOnSamples(df=normdat)
        else:
            #only for one sample
            self.runSample(df=normdat)

    def runSample(self, df, gmt=None, multisamples=False):
        """Single Sample GSEA workflow"""

        assert self.min_size <= self.max_size

        mkdirs(self.outdir)

        #Start Analysis
        self._logger.info("Parsing data files for GSEA.............................")
        #select correct expression genes and values.
        if isinstance(df, pd.DataFrame):
            if df.shape[1] == 1:
                df = df.reset_index()
        elif isinstance(df, pd.Series):
            df = df.reset_index()
        else:
            raise Exception('Error parsing gene ranking values!')

        #sort ranking values from high to low or reverse
        df.sort_values(by=df.columns[1], ascending=self.ascending, inplace=True)
        df.columns = ['gene_name','rank']
        df['rank2'] = df['rank']
        #reset interger index, or caused unwanted problems
        df.reset_index(drop=True, inplace=True)
        # revmove rank2
        dat2 = df.set_index('gene_name')
        dat2 = dat2['rank']
        #cpu numbers
        self._set_cores()
        #filtering out gene sets and build gene sets dictionary
        if gmt is None:
            gmt = gsea_gmt_parser(self.gene_sets,
                                  min_size=self.min_size,
                                  max_size=self.max_size,
                                  gene_list=dat2.index.values)
        self._logger.info("%04d gene_sets used for further statistical testing....."% len(gmt))
        self._logger.info("Start to run GSEA...Might take a while..................")
        #compute ES, NES, pval, FDR, RES
        gsea_results, hit_ind, rank_ES, subsets = gsea_compute_ss(data=dat2, n=self.permutation_num,
                                                                  gmt=gmt,
                                                                  weighted_score_type=self.weighted_score_type,
                                                                  scale=self.scale,
                                                                  seed=self.seed,
                                                                  processes=self._processes)

        self._logger.info("Start to generate gseapy reports, and produce figures...")
        res_zip = zip(subsets, list(gsea_results), hit_ind, rank_ES)

        self._save_results(zipdata=res_zip, outdir=self.outdir, module=self.module,
                                   gmt=gmt, rank_metric=df, permutation_type="gene_sets")

        #Plotting
        self._imat = df
        if not multisamples:
            self._plotting(rank_metric=df, results=self.results, res2d=self.res2d,
                           graph_num=self.graph_num, outdir=self.outdir,
                           figsize=self.figsize, format=self.format, module=self.module)

        self._logger.info("Congratulations. GSEApy run successfully................\n")

        return

    def runOnSamples(self, df):
        """ssGSEA for gct expression matrix
        """

        # df.index.values are gene_names
        #filtering out gene sets and build gene sets dictionary
        gmt = gsea_gmt_parser(self.gene_sets,
                              min_size=self.min_size,
                              max_size=self.max_size,
                              gene_list=df.index.values)

        #Save each sample results to ordereddict
        self.resultsOnSamples = {}
        outdir = self.outdir
        #run ssgsea for gct expression matrixs
        for name, ser in df.iteritems():
            self.outdir= os.path.join(outdir, str(name))
            self._logger.info("Run Sample: %s "%name)
            self.runSample(df=ser, gmt=gmt)
            self.resultsOnSamples[name] = self.res2d.es
        #save raw ES to one csv file
        samplesRawES = pd.DataFrame(self.resultsOnSamples)
        samplesRawES.index.name = 'Term'
        # write es
        outESfile=os.path.join(outdir, "gseapy.samples.raw.es.txt")
        with open(outESfile, 'a') as f:
            if self.scale :
                f.write('# scale the enrichment scores by number of genes in the gene sets\n')
                f.write('# this normalization has not effects on the finall NES ' +\
                        'as indicated by Barbie et al., 2009, online methods, pg. 2\n')
            else:
                f.write('# raw enrichment scores of all data\n')
                f.write('# no scale es by numbers of genes in the gene sets\n')
            samplesRawES.to_csv(f, sep='\t')
        ## normalize enrichment scores by using the entire data set, as indicated
        ## by Barbie et al., 2009, online methods, pg. 2
        samplesNES = samplesRawES / (samplesRawES.values.max() - samplesRawES.values.min())
        outNESfile = os.path.join(outdir, "gseapy.samples.normalized.es.txt")
        with open(outNESfile, 'a') as f:
            f.write('# normalize enrichment scores by using the entire data set\n')
            f.write('# as indicated by Barbie et al., 2009, online methods, pg. 2\n')
            samplesNES.to_csv(f, sep='\t')
        return


class Replot(GSEAbase):
    """To Reproduce GSEA desktop output results."""
    def __init__(self, indir, outdir='GSEApy_Replot', weighted_score_type=1,
                  min_size=3, max_size=1000, figsize=[6.5,6], format='pdf', verbose=False):
        self.indir=indir
        self.outdir=outdir
        self.weighted_score_type=weighted_score_type
        self.min_size=min_size
        self.max_size=max_size
        self.figsize=figsize
        self.format=format
        self.verbose=verbose
        self.module='replot'
        self.gene_sets=None
        self.ascending=False

    def run(self):
        """main replot function"""
        assert self.min_size <= self.max_size

        import glob
        from bs4 import BeautifulSoup

        #parsing files.......
        try:
            results_path = glob.glob(self.indir+'*/edb/results.edb')[0]
            rank_path =  glob.glob(self.indir+'*/edb/*.rnk')[0]
            gene_set_path =  glob.glob(self.indir+'*/edb/gene_sets.gmt')[0]
        except IndexError as e:
            #logger.debug(e)
            sys.stderr.write("Could not locate GSEA files in the given directory!")
            sys.exit(1)
        #extract sample names from .cls file
        cls_path = glob.glob(self.indir+'*/edb/*.cls')
        if cls_path:
            phenoPos, phenoNeg, classes = gsea_cls_parser(cls_path[0])
        else:
            # logic for prerank results
            phenoPos, phenoNeg = '',''
        #start reploting
        self.gene_sets=gene_set_path
        mkdirs(self.outdir)
        logger = self._log_init(module=self.module,
                                log_level=logging.INFO if self.verbose else logging.WARNING)

        #obtain gene sets
        gene_set_dict = gsea_gmt_parser(gene_set_path, min_size=self.min_size, max_size=self.max_size)
        #obtain rank_metrics
        rank_metric = self._load_ranking(rank_path)
        correl_vector =  rank_metric['rank'].values
        gene_list = rank_metric['gene_name']
        #extract each enriment term in the results.edb files and plot.
        database = BeautifulSoup(open(results_path), features='xml')
        length = len(database.findAll('DTG'))

        for idx in range(length):
            #extract statistical resutls from results.edb file
            enrich_term, hit_ind, nes, pval, fdr= gsea_edb_parser(results_path, index=idx)
            gene_set = gene_set_dict.get(enrich_term)
            #calculate enrichment score
            RES = enrichment_score(gene_list=gene_list, gene_set=gene_set,
                                   weighted_score_type=self.weighted_score_type,
                                   correl_vector=correl_vector)[2]
            #plotting
            gsea_plot(rank_metric, enrich_term, hit_ind, nes, pval,
                            fdr, RES, phenoPos, phenoNeg, self.figsize,
                            self.format, self.outdir, self.module)

        logger.info("Congratulations! Your plots have been reproduced successfully!\n")


def call(data, gene_sets, cls, outdir='GSEA_', min_size=15, max_size=500, permutation_num=1000,
          weighted_score_type=1,permutation_type='gene_set', method='log2_ratio_of_classes',
      ascending=False, figsize=[6.5,6], format='pdf', graph_num=20, seed=None, verbose=False):

    sys.stderr.write("DeprecationWarning: "+\
                     "call function has been deprecated, plesea use gseapy.gsea() instead.")

    return


def gsea(data, gene_sets, cls, outdir='GSEA_', min_size=15, max_size=500, permutation_num=1000,
          weighted_score_type=1,permutation_type='gene_set', method='log2_ratio_of_classes',
	  ascending=False, processes=1, figsize=[6.5,6], format='pdf', graph_num=20, seed=None, verbose=False):
    """ Run Gene Set Enrichment Analysis.

    :param data: Gene expression data table, pandas DataFrame, gct file.
    :param gene_sets: Enrichr Library name or .gmt gene sets file. Same input with GSEA.
    :param cls: a list or a .cls file format required for GSEA.
    :param str outdir: Results output directory.
    :param int permutation_num: Number of permutations for significance computation. Default: 1000.
    :param str permutation_type: Permutation type, "phenotype" for phenotypes, "gene_set" for genes.
    :param int min_size: Minimum allowed number of genes from gene set also the data set. Defaut: 15.
    :param int max_size: Maximum allowed number of genes from gene set also the data set. Defaults: 500.
    :param float weighted_score_type: Refer to :func:`algorithm.enrichment_socre`. Default:1.
    :param method:  The method used to calculate a correlation or ranking. Default: 'log2_ratio_of_classes'.
                   Others methods are:

                   1. 'signal_to_noise'

                      You must have at least three samples for each phenotype to use this metric.
                      The larger the signal-to-noise ratio, the larger the differences of the means (scaled by the standard deviations);
                      that is, the more distinct the gene expression is in each phenotype and the more the gene acts as a “class marker.”

                   2. 't_test'

                      Uses the difference of means scaled by the standard deviation and number of samples.
                      Note: You must have at least three samples for each phenotype to use this metric.
                      The larger the tTest ratio, the more distinct the gene expression is in each phenotype
                      and the more the gene acts as a “class marker.”

                   3. 'ratio_of_classes' (also referred to as fold change).

                      Uses the ratio of class means to calculate fold change for natural scale data.

                   4. 'diff_of_classes'


                      Uses the difference of class means to calculate fold change for nature scale data


                   5. 'log2_ratio_of_classes'

                      Uses the log2 ratio of class means to calculate fold change for natural scale data.
                      This is the recommended statistic for calculating fold change for log scale data.


    :param bool ascending: Sorting order of rankings. Default: False.
    :param int processes: Number of Processes you are going to use. Default: 1.
    :param list figsize: Matplotlib figsize, accept a tuple or list, e.g. [width,height]. Default: [6.5,6].
    :param str format: Matplotlib figure format. Default: 'pdf'.
    :param int graph_num: Plot graphs for top sets of each phenotype
    :param seed: Random seed. expect an interger. Defalut:None.
    :param bool verbose: Bool, increase output verbosity, print out progress of your job, Default: False.

    :return: Return a GSEA obj. All results store to a dictionary, obj.results,
             where contains::

                 | {es: enrichment score,
                 |  nes: normalized enrichment score,
                 |  p: P-value,
                 |  fdr: FDR,
                 |  size: gene set size,
                 |  matched_size: genes matched to the data,
                 |  genes: gene names from the data set}

    """
    gs = GSEA(data, gene_sets, cls, outdir, min_size, max_size, permutation_num,
              weighted_score_type, permutation_type, method, ascending, processes,
               figsize, format, graph_num, seed, verbose)
    gs.run()

    return gs


def ssgsea(data, gene_sets, outdir="GSEA_SingleSample", sample_norm_method='rank', min_size=15, max_size=2000,
           permutation_num=1000, weighted_score_type=0.25, scale=True, ascending=False, processes=1,
           figsize=[7,6], format='pdf', graph_num=20, seed=None, verbose=False):
    """Run Gene Set Enrichment Analysis with single sample GSEA tool

    :param data: expression table, pd.Series, pd.DataFrame, GCT file, or .rnk file formate.
    :param gene_sets: Enrichr Library name or .gmt gene sets file. Same input with GSEA.
    :param outdir: results output directory.
    :param str sample_norm_method: "Sample normalization method. Choose from {'rank', 'log', 'log_rank'}. Default: rank.

               1. 'rank': Rank your expression data, and transformed by 10000*rank_dat/gene_numbers
               2. 'log' : Do not rank, but transformed data by log(data + exp(1)), while  data = data[data<1] =1.
               3. 'log_rank': Rank your expression data, and transformed by log(10000*rank_dat/gene_numbers+ exp(1))
               4. 'custom': Do nothing, and use your own rank value to calulate enrichment score.
               see here: https://github.com/GSEA-MSigDB/ssGSEAProjection-gpmodule/blob/master/src/ssGSEAProjection.Library.R, line 86

    :param int min_size: Minimum allowed number of genes from gene set also the data set. Defaut: 15.
    :param int max_size: Maximum allowed number of genes from gene set also the data set. Defaults: 2000.
    :param int permutation_num: Number of permutations for significance computation. Default: 1000.
    :param str weighted_score_type: Refer to :func:`algorithm.enrichment_socre`. Default:0.25.
    :param bool scale: If True, normalize the scores by number of genes in the gene sets.
    :param bool ascending: Sorting order of rankings. Default: False.
    :param int processes: Number of Processes you are going to use. Default: 1.
    :param list figsize: Matplotlib figsize, accept a tuple or list, e.g. [width,height]. Default: [7,6].
    :param str format: Matplotlib figure format. Default: 'pdf'.
    :param int graph_num: Plot graphs for top sets of each phenotype
    :param seed: Random seed. expect an interger. Defalut:None.
    :param bool verbose: Bool, increase output verbosity, print out progress of your job, Default: False.

    :return: Return a ssGSEA obj. All results store to  a dictionary, obj.results,
             where contains::

                 | {es: enrichment score,
                 |  nes: normalized enrichment score,
                 |  p: P-value,
                 |  fdr: FDR,
                 |  size: gene set size,
                 |  matched_size: genes matched to the data,
                 |  genes: gene names from the data set}

    """

    ss = SingleSampleGSEA(data, gene_sets, outdir, sample_norm_method, min_size, max_size,
                          permutation_num, weighted_score_type, scale, ascending,
                          processes, figsize, format, graph_num, seed, verbose)
    ss.run()
    return ss


def prerank(rnk, gene_sets, outdir='GSEA_Prerank', pheno_pos='Pos', pheno_neg='Neg',
            min_size=15, max_size=500, permutation_num=1000, weighted_score_type=1,
            ascending=False, processes=1, figsize=[6.5,6], format='pdf',
            graph_num=20, seed=None, verbose=False):
    """ Run Gene Set Enrichment Analysis with pre-ranked correlation defined by user.

    :param rnk: pre-ranked correlation table or pandas DataFrame. Same input with ``GSEA`` .rnk file.
    :param gene_sets: Enrichr Library name or .gmt gene sets file. Same input with GSEA.
    :param outdir: results output directory.
    :param int permutation_num: Number of permutations for significance computation. Default: 1000.
    :param int min_size: Minimum allowed number of genes from gene set also the data set. Defaut: 15.
    :param int max_size: Maximum allowed number of genes from gene set also the data set. Defaults: 500.
    :param str weighted_score_type: Refer to :func:`algorithm.enrichment_socre`. Default:1.
    :param bool ascending: Sorting order of rankings. Default: False.
    :param int processes: Number of Processes you are going to use. Default: 1.
    :param list figsize: Matplotlib figsize, accept a tuple or list, e.g. [width,height]. Default: [6.5,6].
    :param str format: Matplotlib figure format. Default: 'pdf'.
    :param int graph_num: Plot graphs for top sets of each phenotype
    :param seed: Random seed. expect an interger. Defalut:None.
    :param bool verbose: Bool, increase output verbosity, print out progress of your job, Default: False.

    :return: Return a Prerank obj. All results store to  a dictionary, obj.results,
             where contains::

                 | {es: enrichment score,
                 |  nes: normalized enrichment score,
                 |  p: P-value,
                 |  fdr: FDR,
                 |  size: gene set size,
                 |  matched_size: genes matched to the data,
                 |  genes: gene names from the data set}

    """
    pre = Prerank(rnk, gene_sets, outdir, pheno_pos, pheno_neg,
                  min_size, max_size, permutation_num, weighted_score_type,
                  ascending, processes, figsize, format, graph_num, seed, verbose)
    pre.run()
    return pre



def replot(indir, outdir='GSEA_Replot', weighted_score_type=1,
           min_size=3, max_size=1000, figsize=[6.5,6], format='pdf', verbose=False):
    """The main fuction to reproduce GSEA desktop outputs.

    :param indir: GSEA desktop results directory. In the sub folder, you must contain edb file foder.
    :param outdir: Output directory.
    :param float weighted_score_type: weighted score type. choose from {0,1,1.5,2}. Default: 1.
    :param list figsize: matplotlib output figure figsize. Defult: [6.5,6].
    :param str format: matplotlib output figure format. Default: 'pdf'.
    :param int min_size: min size of input genes presented in Gene Sets. Default: 3.
    :param int max_size: max size of input genes presented in Gene Sets. Default: 5000.
                     you will not encourage to use min_size, or max_size argment in :func:`replot` function.
                     Because gmt file has already been filter.
    :param verbose: Bool, increase output verbosity, print out progress of your job, Default: False.

    :return: Generate new figures with seleted figure format. Default: 'pdf'.

    """
    rep = Replot(indir, outdir, weighted_score_type,
                 min_size, max_size, figsize, format, verbose)
    rep.run()

    return
