#! python
# -*- coding: utf-8 -*-

import os,sys, logging
import pandas as pd
from collections import OrderedDict
from numpy import number
from gseapy.parser import *
from gseapy.algorithm import enrichment_score, gsea_compute, gsea_compute_ss, ranking_metric
from gseapy.plot import gsea_plot, heatmap
from gseapy.utils import mkdirs


class GSEAbase:
    """base class of GSEA."""
    def __init__(self):
        self.gene_sets='KEGG_2016'
        self.module='base'
        self.results=None
        self.res2d=None
        self.ranking=None
        self.verbose=False
        self._logger=None
       

    def _log_init(self, module='GSEA', log_level=logging.INFO):
        """logging start"""
        gene_set =os.path.split(self.gene_sets)[-1].lower().rstrip(".gmt")
        logging.basicConfig(
                            level    = logging.DEBUG,
                            format   = 'LINE %(lineno)-4d: %(asctime)s [%(levelname)-8s] %(message)s',
                            filename = "%s/gseapy.%s.%s.log"%(self.outdir, module, gene_set),
                            filemode = 'w')

        logger = logging.getLogger(__name__)
        #logger.setLevel(logging.DEBUG)
        # define a Handler which writes INFO messages or higher to the sys.stderr
        console = logging.StreamHandler()
        console.setLevel(log_level)
        # set a format which is simpler for console use
        formatter = logging.Formatter('%(asctime)s %(message)s')
        # tell the handler to use this format
        console.setFormatter(formatter)
        logger.addHandler(console)
        #if you want information print to the console, uisng logger.info....
        self._logger=logger

        return self._logger
    def _log_stop(self):
        """log stop"""

        handlers = self._logger.handlers[:]
        for handler in handlers:
            handler.close()
            self._logger.removeHandler(handler)     
    def _rank_metric(self, rnk):
        """Parse ranking file. This file contains ranking correlation vector and gene names or ids. 
        
        :param rnk: the .rnk file of GSEA input or a ranked pandas DataFrame instance.
        :return: a pandas DataFrame with 3 columns names are: 'gene_name','rank',rank2'
                     
        """
        if isinstance(rnk, pd.DataFrame) :
            rank_metric = rnk.copy()
        elif isinstance(rnk, str) :
            rank_metric = pd.read_table(rnk, header=None)
        else:
            raise Exception('Error parsing gene ranking values!')
            
        rank_metric.columns = ['gene_name','rank']
        rank_metric['rank2'] = rank_metric['rank']
               
        self.ranking = rank_metric
        return self.ranking

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
        if module == 'gsea':
            width = len(classes) if len(classes) >= 6 else  5
            datA = data.loc[:, cls_booA]
            datB = data.loc[:, cls_booB]
            datAB=pd.concat([datA,datB], axis=1)

        cls_booA =list(map(lambda x: True if x == phenoPos else False, classes))
        cls_booB =list(map(lambda x: True if x == phenoNeg else False, classes))
        
        for gs in top_term:
            hit = results.get(gs)['hit_index']
            fig = gsea_plot(rank_metric=rank_metric, enrich_term=gs, hit_ind=hit,
                            nes=results.get(gs)['nes'], pval=results.get(gs)['pval'], fdr=results.get(gs)['fdr'], 
                            RES=results.get(gs)['rank_ES'], phenoPos=phenoPos, phenoNeg=phenoNeg, figsize=figsize)        
            gs = gs.replace('/','_').replace(":","_")
            fig.savefig('{a}/{b}.{c}.{d}'.format(a=outdir, b=gs, c=module, d=format), bbox_inches='tight', dpi=300,)
             
            if module == 'gsea':
                heatmap(df=datAB.iloc[hit], term=gs, outdir=outdir, 
                        figsize=(width, len(hit)/2), format=format)
      


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
            rdict['genes'] = _genes.to_string(header=False, index=False).replace("\n", ",")

            rdict['rank_ES'] = RES
            rdict['hit_index'] = ind
            res[gs] = rdict           

        res_df = DataFrame.from_dict(res, orient='index')
        res_df.index.name = 'Term'
        res_df.sort_values(by='fdr', inplace=True)
        res_df.drop(['rank_ES','hit_index'], axis=1, inplace=True)
        res_df.to_csv('{a}/gseapy.{b}.{c}.report.csv'.format(a=outdir, b=module, c=permutation_type),
                      float_format ='%.7f')    

        self.res2d = res_df
        self.results  = res
        return 

    
class GSEA(GSEAbase):
    """GSEA main tool"""
    def __init__(self, data, gene_sets, classes, outdir='GSEA_ouput', 
                 min_size=15, max_size=500, permutation_num=1000, 
                 weighted_score_type=1,permutation_type='gene_set', 
                 method='log2_ratio_of_classes', ascending=False, 
                 figsize=[6.5,6], format='pdf', graph_num=20, 
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
        self.figsize=figsize
        self.format=format
        self.graph_num=graph_num
        self.seed=seed
        self.verbose=verbose
        self.module='gsea'
    def __preprocess(self, df, cls_vector):
        """pre-processed the data frame.new filtering methods will be implement here.
        """    
        
        df.drop_duplicates(subset=df.columns[0], inplace=True) #drop duplicate gene_names.    
        df.set_index(keys=df.columns[0], inplace=True)
        df.dropna(how='all', inplace=True)                     #drop rows with all NAs
        df2 = df.select_dtypes(include=[number])  + 0.00001 #select numbers in DataFrame
        
        #drop any genes which std ==0
        df_std =  df2.groupby(by=cls_vector, axis=1).std()    
        df2 =  df2[~df_std.isin([0]).any(axis=1)]     
        
        return df2

    def run(self):
        """GSEA main procedure"""

        assert self.permutation_type in ["phenotype", "gene_set"]
        assert self.min_size <= self.max_size

        if isinstance(self.data, pd.DataFrame) :
            df = self.data.copy()

        elif isinstance(self.data, str) :
            df = pd.read_table(self.data)
        else:
            raise Exception('Error parsing gene expression dataframe!')
            sys.exit(1)
   
        assert len(df) > 1

        mkdirs(self.outdir)
        logger = self._log_init(module=self.module, 
                               log_level=logging.INFO if self.verbose else logging.WARNING)

        #Start Analysis
        logger.info("Parsing data files for GSEA.............................")     

        # phenotype labels parsing
        phenoPos, phenoNeg, cls_vector = gsea_cls_parser(self.classes)
        #select correct expression genes and values.
        dat = self.__preprocess(df, cls_vector)
        
        
        #ranking metrics calculation.    
        dat2 = ranking_metric(df=dat, method=self.method, phenoPos=phenoPos, phenoNeg=phenoNeg, 
                              classes=cls_vector, ascending=self.ascending)

        #filtering out gene sets and build gene sets dictionary
        gmt = gsea_gmt_parser(self.gene_sets, min_size=self.min_size, max_size=self.max_size, 
                              gene_list=dat2['gene_name'].values)

        logger.info("%04d gene_sets used for further statistical testing....."% len(gmt))

        logger.info("Start to run GSEA...Might take a while..................")   
        #compute ES, NES, pval, FDR, RES
        gsea_results,hit_ind,rank_ES, subsets = gsea_compute(data=dat, n=self.permutation_num, gmt=gmt, 
                                                             weighted_score_type=self.weighted_score_type,
                                                             permutation_type=self.permutation_type, 
                                                             method=self.method,
                                                             phenoPos=phenoPos, phenoNeg=phenoNeg, 
                                                             classes=cls_vector, 
                                                             ascending=self.ascending, seed=self.seed)                                                             
        
        logger.info("Start to generate gseapy reports, and produce figures...")
        res_zip = zip(subsets, list(gsea_results), hit_ind, rank_ES)

        self._save_results(zipdata=res_zip, outdir=self.outdir, module=self.module, 
                                   gmt=gmt, rank_metric=dat2, permutation_type="gene_sets")
        
        #Plotting
        self._plotting(rank_metric=dat2, results=self.results, res2d=self.res2d, 
                       graph_num=self.graph_num, outdir=self.outdir,
                       figsize=self.figsize, format=self.format, module=self.module,
                       data=dat, classes=cls_vector, phenoPos=phenoPos, phenoNeg=phenoNeg)
        
        logger.info("Congratulations. GSEApy run successfully................")
        self._log_stop()

        return         

class Prerank(GSEAbase):
    """GSEA prerank tool"""
    def __init__(self, rnk, gene_sets, outdir='GSEA_prerank', 
                 pheno_pos='Pos', pheno_neg='Neg', min_size=15, max_size=500, 
                 permutation_num=1000, weighted_score_type=1,
                 ascending=False, figsize=[6.5,6], format='pdf', 
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
        self.module='prerank'

    def run(self):
        """GSEA prerank workflow"""

        assert self.min_size <= self.max_size
        mkdirs(self.outdir)
        logger = self._log_init(module=self.module, 
                               log_level=logging.INFO if self.verbose else logging.WARNING)


        dat2 = self._rank_metric(self.rnk)  
        #drop duplicates in ranking metrics. 
        dat2.drop_duplicates(subset='gene_name',inplace=True, keep='first')     
        assert len(dat2) > 1
        
        #Start Analysis
        logger.info("Parsing data files for GSEA.............................")     

        #filtering out gene sets and build gene sets dictionary
        gmt = gsea_gmt_parser(self.gene_sets, min_size=self.min_size, max_size=self.max_size,
                              gene_list=dat2['gene_name'].values)
        logger.info("%04d gene_sets used for further statistical testing....."% len(gmt))

        
        logger.info("Start to run GSEA...Might take a while..................")   
        #compute ES, NES, pval, FDR, RES
        gsea_results, hit_ind,rank_ES, subsets = gsea_compute(data=dat2, n=self.permutation_num, gmt=gmt, 
                                                              weighted_score_type=self.weighted_score_type,
                                                              permutation_type='gene_set', method=None, 
                                                              phenoPos=self.pheno_pos, phenoNeg=self.pheno_neg,
                                                              classes=None, ascending=self.ascending, seed=self.seed, 
                                                              prerank=True)

        logger.info("Start to generate gseapy reports, and produce figures...")
        res_zip = zip(subsets, list(gsea_results), hit_ind, rank_ES)

        self._save_results(zipdata=res_zip, outdir=self.outdir, module=self.module, 
                                   gmt=gmt, rank_metric=dat2, permutation_type="gene_sets")
        
        #Plotting
        self._plotting(rank_metric=dat2, results=self.results, res2d=self.res2d, 
                       graph_num=self.graph_num, outdir=self.outdir,
                       figsize=self.figsize, format=self.format, module=self.module)
        
        logger.info("Congratulations. GSEApy run successfully................")
        self._log_stop()

        return 
    

class SingleSampleGSEA(GSEAbase):
    """GSEA extention: single sample GSEA"""
    def __init__(self, data, gene_sets, outdir="GSEA_SingleSample",
                 min_size=15, max_size=500, permutation_num=1000, weighted_score_type=0.25,
                 ascending=False, figsize=[6.5,6], format='pdf',
                 graph_num=20, seed=None, verbose=False):
        self.data = data
        self.gene_sets=gene_sets
        self.outdir=outdir
        self.weighted_score_type=weighted_score_type
        self.min_size=min_size
        self.max_size=max_size
        self.permutation_num=permutation_num
        self.ascending=ascending
        self.figsize=figsize
        self.format=format
        self.graph_num=graph_num
        self.seed=seed
        self.verbose=verbose
        self.module='SingleSample'

    def run(self):
        """Single Sample GSEA workflow"""

        assert self.min_size <= self.max_size

        mkdirs(self.outdir)
        logger = self._log_init(module=self.module, 
                               log_level=logging.INFO if self.verbose else logging.WARNING)


        dat = self._rank_metric(self.data)  
   
        assert len(dat) > 1
        
        #Start Analysis
        logger.info("Parsing data files for GSEA.............................")     
        #select correct expression genes and values.

        dat2 = dat.set_index('gene_name')
        del dat2['rank2']
        #filtering out gene sets and build gene sets dictionary
        gmt = gsea_gmt_parser(self.gene_sets, min_size=self.min_size, max_size=self.max_size, gene_list=dat2.index.values)
        logger.info("%04d gene_sets used for further statistical testing....."% len(gmt))

        logger.info("Start to run GSEA...Might take a while..................") 
        #compute ES, NES, pval, FDR, RES
        gsea_results, hit_ind, rank_ES, subsets = gsea_compute_ss(data=dat2, n=self.permutation_num, gmt=gmt,
                                                                  weighted_score_type=self.weighted_score_type,
                                                                  seed=self.seed)
          
        logger.info("Start to generate gseapy reports, and produce figures...")
        res_zip = zip(subsets, list(gsea_results), hit_ind, rank_ES)

        self._save_results(zipdata=res_zip, outdir=self.outdir, module=self.module, 
                                   gmt=gmt, rank_metric=dat, permutation_type="gene_sets")
        
        #Plotting
        self._plotting(rank_metric=dat, results=self.results, res2d=self.res2d, 
                       graph_num=self.graph_num, outdir=self.outdir,
                       figsize=self.figsize, format=self.format, module=self.module)
        
        logger.info("Congratulations. GSEApy run successfully................")
        self._log_stop()

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

    def run(self):
        """main replot function"""
        assert self.min_size <= self.max_size

        mkdirs(self.outdir)
        logger = self._log_init(module=self.module,  
                                log_level=logging.INFO if self.verbose else logging.WARNING)
        import glob
        from bs4 import BeautifulSoup
        
        #parsing files.......    
        try:
            results_path = glob.glob(self.indir+'*/edb/results.edb')[0]
            rank_path =  glob.glob(self.indir+'*/edb/*.rnk')[0]
            gene_set_path =  glob.glob(self.indir+'*/edb/gene_sets.gmt')[0]
        except IndexError as e: 
            logger.debug(e) 
            logger.error("Could not locate GSEA files in the given directory!")
            sys.exit(1)    
        #extract sample names from .cls file
        cls_path = glob.glob(self.indir+'*/edb/*.cls')
        if cls_path:
            phenoPos, phenoNeg, classes = gsea_cls_parser(cls_path[0])
        else:
            # logic for prerank results
            phenoPos, phenoNeg = '',''  
        #obtain gene sets
        gene_set_dict = gsea_gmt_parser(gene_set_path, min_size=self.min_size, max_size=self.max_size)
        #obtain rank_metrics
        rank_metric = self._rank_metric(rank_path)
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
            fig = gsea_plot(rank_metric, enrich_term, hit_ind, nes, pval,
                            fdr, RES, phenoPos, phenoNeg, figsize=self.figsize)    
            fig.savefig('{a}/{b}.gsea.{c}.{d}'.format(a=self.outdir, b=enrich_term,
                                                      c=self.module, d=self.format),
                        bbox_inches='tight', dpi=300,)

          
        logger.info("Congratulations! Your plots have been reproduced successfully!")
        self._log_stop()


def call(data, gene_sets, cls, outdir='GSEA_', min_size=15, max_size=500, permutation_num=1000, 
          weighted_score_type=1,permutation_type='gene_set', method='log2_ratio_of_classes',
      ascending=False, figsize=[6.5,6], format='pdf', graph_num=20, seed=None, verbose=False):

    sys.stderr.write("DeprecationWarning: "+\
                     "call function has been deprecated, plesea use gseapy.gsea() instead.")

    return

        
def gsea(data, gene_sets, cls, outdir='GSEA_', min_size=15, max_size=500, permutation_num=1000, 
          weighted_score_type=1,permutation_type='gene_set', method='log2_ratio_of_classes',
	  ascending=False, figsize=[6.5,6], format='pdf', graph_num=20, seed=None, verbose=False):
    """ Run Gene Set Enrichment Analysis.

    :param data: Gene expression data table.  
    :param gene_sets: Gene sets file. e.g. gmt files. Same input with GSEA.
    :param permutation_n: Number of permutations for significance computation. Default: 1000.
    :param permutation_type: Permutation type, "phenotype" for phenotypes, "gene_set" for genes.
    :param int min_size: Minimum allowed number of genes from gene set also the data set. Defaut: 15.
    :param int max_size: Maximum allowed number of genes from gene set also the data set. Defaults: 500.
    :param weighted_score_type: Refer to :func:`algorithm.enrichment_socre`. Default:1.
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
      
                      Uses the difference of class means to calculate fold change for log scale data
    
                   5. 'log2_ratio_of_classes' 
      
                      Uses the log2 ratio of class means to calculate fold change for natural scale data.
                      This is the recommended statistic for calculating fold change for natural scale data.
   
      	
    :param ascending: Sorting order of rankings. Default: False.
    :param outdir: Results output directory.
    :param figsize: Matplotlib figsize, accept a tuple or list, e.g. [width,height]. Default: [6.5,6].
    :param format: Matplotlib figure format. Default: 'pdf'.
    :param graph_num: Plot graphs for top sets of each phenotype
    :param seed: Random seed. expect an interger. Defalut:None.
    :param verbose: Bool, increase output verbosity, print out progress of your job, Default: False.

    :return: Return a DataFrame when inside python console.
             Generate ``GSEA`` plots and store a dictionary into csv file,
             where  each column represents::

                 | {es: enrichment score, 
                 |  nes: normalized enrichment score, 
                 |  p: P-value, 
                 |  fdr: FDR, 
                 |  size: gene set size,
                 |  matched_size: genes matched to the data, 
                 |  genes: gene names from the data set}
    
    """
    gs = GSEA(data, gene_sets, cls, outdir, min_size, max_size, permutation_num, 
              weighted_score_type,permutation_type, method, ascending,
               figsize, format, graph_num, seed, verbose)
    gs.run()

    return gs


def ssgsea(data, gene_sets, outdir="GSEA_SingleSample", min_size=15, max_size=500,
           permutation_num=1000, weighted_score_type=0.25, ascending=False, 
           figsize=[6.5,6], format='pdf', graph_num=20, seed=None, verbose=False):
    """single sample GSEA tool"""
    ss = SingleSampleGSEA(data, gene_sets, outdir, min_size, max_size, 
                          permutation_num, weighted_score_type,
                          ascending, figsize, format, graph_num, seed, verbose)
    ss.run()
    return ss
    

def prerank(rnk, gene_sets, outdir='GSEA_Prerank', pheno_pos='Pos', pheno_neg='Neg',
            min_size=15, max_size=500, permutation_num=1000, weighted_score_type=1,
            ascending=False, figsize=[6.5,6], format='pdf', graph_num=20, seed=None, verbose=False):
    """ Run Gene Set Enrichment Analysis with pre-ranked correlation defined by user.

    :param rnk: pre-ranked correlation table, Same input with ``GSEA`` .rnk file.  
    :param gene_sets: Gene sets file. e.g. gmt files. Same input with GSEA.
    :param outdir: results output directory.
    :param permutation_n: Number of permutations for significance computation. Default: 1000.
    :param int min_size: Minimum allowed number of genes from gene set also the data set. Defaut: 15.
    :param int max_size: Maximum allowed number of genes from gene set also the data set. Defaults: 500.
    :param weighted_score_type: Refer to :func:`algorithm.enrichment_socre`. Default:1.
    :param ascending: Sorting order of rankings. Default: False.
    :param figsize: Matplotlib figsize, accept a tuple or list, e.g. [width,height]. Default: [6.5,6].
    :param format: Matplotlib figure format. Default: 'pdf'.
    :param graph_num: Plot graphs for top sets of each phenotype
    :param seed: Random seed. expect an interger. Defalut:None. 
    :param verbose: Bool, increase output verbosity, print out progress of your job, Default: False.
       
    :return: Return a DataFrame when inside python console.
             Generate ``GSEA`` plots and store a dictionary into csv file,
             where each column represents::

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
                  ascending, figsize, format, graph_num, seed, verbose)
    pre.run()
    return pre



def replot(indir, outdir='GSEA_Replot', weighted_score_type=1, 
           min_size=3, max_size=1000, figsize=[6.5,6], format='pdf', verbose=False):
    """The main fuction to reproduce GSEA desktop outputs.
          
    :param indir: GSEA desktop results directory. In the sub folder, you must contain edb file foder.    
    :param outdir: Output directory.
    :param weight: weighted score type. choose from {0,1,1.5,2}. Default: 1.
    :param figsize: matplotlib output figure figsize. Defult: [6.5,6].
    :param format: matplotlib output figure format. Default: 'pdf'.
    :param min_size: min size of input genes presented in Gene Sets. Default: 3.
    :param max_size: max size of input genes presented in Gene Sets. Default: 5000.
                     you will not encourage to use min_size, or max_size argment in :func:`replot` function.
                     Because gmt file has already been filter.
    :param verbose: Bool, increase output verbosity, print out progress of your job, Default: False.

    :return: Generate new figures with seleted figure format. Default: 'pdf'.   
    """
    rep = Replot(indir, outdir, weighted_score_type, 
                 min_size, max_size, figsize, format, verbose)
    rep.run()

    return 