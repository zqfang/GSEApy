#! python
# -*- coding: utf-8 -*-
from __future__ import  division


import os,sys, logging
from .parser import *
from .algorithm import enrichment_score, gsea_compute, gsea_compute_ss, preprocess, ranking_metric
from .plot import gsea_plot, heatmap
from .utils import log_init, log_remove, mkdirs, save_results
import pandas as pd



class GSEAbase:
    def __init__(self):
        self.results=None
        self.verbose=False
        self.module=None
        self.logger=None
       
    def savefig(self, fig):
        #fig.savefig()
        return
    def to_csv(self, res):
        self.results = res
        return self.results

    def log_start(self):
        verobse = logging.INFO if self.verbose else logging.WARNING
        self.logger = log_init(self.outdir, module=self.module, log_level=logging.INFO)
        return self.logger
    def log_stop(self):
         log_remove(self.logger)
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
        self.weight_score_type=weight_score_type
        self.ascending=ascending
        self.figsize=figsize
        self.format=format
        self.graph_num=graph_num
        self.seed=seed
        self.verbose=verbose

    def run(self):

        self.results=call(self.data, self.gene_sets, self.classes, self.outdir, 
                          self.min_size, self.max_size, self.permutation_num, 
                          self.weighted_score_type, self.permutation_type,self.method,
                          self.ascending, self.figsize, self.format, self.graph_num, 
                          self.seed, self.verbose)

        return  self.results


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
        self.weight_score_type=weight_score_type
        self.ascending=ascending
        self.figsize=figsize
        self.format=format
        self.graph_num=graph_num
        self.seed=seed
        self.verbose=verbose
    def run(self):
        self.results = prerank(self.rnk, self.gene_sets, self.outdir, 
                               self.pheno_pos, self.pheno_neg, self.min_size, self.max_size, 
                               self.permutation_num, self.weighted_score_type, 
                               self.ascending, self.figsize, self.format, 
                               self.graph_num, self.seed, self.verbose)
        return self.results
    

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
        self.results=None

    def run(self):
        
        mkdirs(self.outdir)
        logger = self.log_start()

        if isinstance(self.data, pd.DataFrame) :
            df = self.data.copy()
        elif isinstance(self.data, str):
            df = pd.read_table(self.data)
        else:
            raise Exception('Error parsing gene expression dataframe!')
            sys.exit(1)
       
        assert len(df) > 1
        #write command to log file

        #Start Analysis
        logger.info("Parsing data files for GSEA.............................")     
        #select correct expression genes and values.
        dat = gsea_rank_metric(df)
        dat2 = dat.set_index('gene_name')
        del dat2['rank2']
        #filtering out gene sets and build gene sets dictionary
        gmt = gsea_gmt_parser(self.gene_sets, min_size=self.min_size, max_size=self.max_size, gene_list=dat2.index.values)
        logger.info("%04d gene_sets used for further statistical testing....."% len(gmt))

        
        logger.info("Start to run GSEA...Might take a while..................")   
        #compute ES, NES, pval, FDR, RES
        results,hit_ind,rank_ES, subsets = gsea_compute_ss(data=dat2, n=self.permutation_num, gmt=gmt,
                                                           weighted_score_type=self.weighted_score_type,
                                                           seed=self.seed)
        logger.info("Start to generate gseapy reports, and produce figures...")
        res_zip = zip(subsets, list(results), hit_ind, rank_ES)
        res, res_df = save_results(zipdata=res_zip, outdir=self.outdir, module='SingleSample', 
                                   gmt=gmt, data=dat, permutation_type="gene_sets")

        #Plotting
        top_term = res_df.head(self.graph_num).index
        
        for gs in top_term:
            hit = res.get(gs)['hit_index']
            gene_symbol = res.get(gs)['genes']
            fig = gsea_plot(rank_metric=dat, enrich_term=gs, hit_ind=hit,
                            nes=res.get(gs)['nes'], pval=res.get(gs)['pval'], fdr=res.get(gs)['fdr'], 
                            RES=res.get(gs)['rank_ES'], phenoPos="", phenoNeg="", figsize=self.figsize)        
            gs = gs.replace('/','_').replace(":","_")
            fig.savefig('{a}/{b}.gsea.{c}'.format(a=outdir, b=gs, c=format), bbox_inches='tight', dpi=300,)

        
        logger.info("Congratulations. GSEApy run successfully................")
        self.log_stop()
        self.results = res
        return self.results 
        



        
def call(data, gene_sets, cls, outdir='gseapy_out', min_size=15, max_size=500, permutation_n=1000, 
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
    argument = locals()
    assert permutation_type in ["phenotype", "gene_set"]
    assert min_size <= max_size

    mkdirs(outdir)
    logger = log_init(outdir, module='call', log_level = logging.INFO if verbose else logging.WARNING)

    if isinstance(data, pd.DataFrame) :
        df = data.copy()
        argument['data'] = 'DataFrame'
    elif isinstance(data, str) :
        df = pd.read_table(data)
    else:
        raise Exception('Error parsing gene expression dataframe!')
        sys.exit(1)
   
    assert len(df) > 1
    #write command to log file
    argument = OrderedDict(sorted(argument.items(), key=lambda t:t[0]))
    logger.debug("Command: call, "+str(argument))
    #Start Analysis
    logger.info("Parsing data files for GSEA.............................")     
    #select correct expression genes and values.
    dat = preprocess(df)
    
    # phenotype labels parsing
    phenoPos, phenoNeg, classes = gsea_cls_parser(cls)
    
    #ranking metrics calculation.    
    dat2 = ranking_metric(df=dat, method=method, phenoPos=phenoPos, phenoNeg=phenoNeg, classes=classes, ascending=ascending)
    
    #filtering out gene sets and build gene sets dictionary
    gmt = gsea_gmt_parser(gene_sets, min_size=min_size, max_size=max_size,gene_list=dat2['gene_name'].values)
    logger.info("%04d gene_sets used for further statistical testing....."% len(gmt))

    logger.info("Start to run GSEA...Might take a while..................")   
    #compute ES, NES, pval, FDR, RES
    results,hit_ind,rank_ES, subsets = gsea_compute(data=dat, n=permutation_n,gmt=gmt, weighted_score_type=weighted_score_type,
                                                    permutation_type=permutation_type, method=method,
                                                    phenoPos=phenoPos, phenoNeg=phenoNeg, classes=classes, ascending=ascending,
                                                    seed=seed)
    logger.info("Start to generate gseapy reports, and produce figures...")
    res_zip = zip(subsets, list(results), hit_ind, rank_ES)
    res, res_df = save_results(zipdata=res_zip, outdir=outdir, 
                               module='GSEA', gmt=gmt,data=dat2, permutation_type=permutation_type)

    #Plotting
    top_term = res_df.head(graph_num).index
    width = len(classes) if len(classes) >= 6 else  5
    for gs in top_term:
        hit = res.get(gs)['hit_index']
        gene_symbol = res.get(gs)['genes']
        fig = gsea_plot(rank_metric=dat2, enrich_term=gs, hit_ind=hit,
                        nes=res.get(gs)['nes'], pval=res.get(gs)['pval'], fdr=res.get(gs)['fdr'], 
                        RES=res.get(gs)['rank_ES'], phenoPos=phenoPos, phenoNeg=phenoNeg, figsize=figsize)        
        gs = gs.replace('/','_').replace(":","_")
        fig.savefig('{a}/{b}.gsea.{c}'.format(a=outdir, b=gs, c=format), bbox_inches='tight', dpi=300,)

        heatmap(df=dat.loc[gene_symbol], term=gs, outdir=outdir, 
                figsize=(width, len(gene_symbol)/2), format=format)
      
    logger.info("Congratulations. GSEApy run successfully................")
    
	# return dataframe if run gsea inside python console
    #if isinstance(data, pd.DataFrame) or isinstance(cls, list):

    log_remove(logger)

    return res


def prerank(rnk, gene_sets, outdir='gseapy_out', pheno_pos='Pos', pheno_neg='Neg',
            min_size=15, max_size=500, permutation_n=1000, weighted_score_type=1,
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
    argument = locals()
    assert min_size <= max_size
    
    mkdirs(outdir)
    logger = log_init(outdir, module='prerank', log_level= logging.INFO if verbose else logging.WARNING)
    if isinstance(rnk, pd.DataFrame) :       
        argument['rnk'] = 'DataFrame'
    #write command to log file
    argument = OrderedDict(sorted(argument.items(), key=lambda t:t[0]))
    logger.debug("Command: prerank, "+str(argument))
    #Start Analysis
    logger.info("Parsing data files for GSEA.............................") 
    dat2 = gsea_rank_metric(rnk)
    assert len(dat2) > 1
    #drop duplicates in ranking metrics. 
    dat2.drop_duplicates(subset='gene_name',inplace=True, keep='first')   
    #filtering out gene sets and build gene sets dictionary
    gmt = gsea_gmt_parser(gene_sets, min_size=min_size, max_size=max_size, gene_list=dat2['gene_name'].values)
    logger.info("%04d gene_sets used for further statistical testing....."% len(gmt))   
    logger.info("Start to run GSEA...Might take a while..................") 
    #compute ES, NES, pval, FDR, RES
    results,hit_ind,rank_ES, subsets = gsea_compute(data=dat2, n=permutation_n, gmt=gmt, weighted_score_type=weighted_score_type,
                                                    permutation_type='gene_set', method=None, phenoPos=pheno_pos, phenoNeg=pheno_neg,
                                                    classes=None, ascending=ascending, seed=seed, prerank=True)
   
    logger.info("Start to generate gseapy reports, and produce figures...")
    res_zip = zip(subsets, list(results), hit_ind, rank_ES)
    res, res_df = save_results(zipdata=res_zip, outdir=outdir, 
                               module='prerank', gmt=gmt, data=dat2, permutation_type=permutation_type)
    

    #Plotting
    top_term = res_df.head(graph_num).index
    for gs in top_term:
        fig = gsea_plot(rank_metric=dat2, enrich_term=gs, hit_ind=res.get(gs)['hit_index'],
                        nes=res.get(gs)['nes'], pval=res.get(gs)['pval'], fdr=res.get(gs)['fdr'], 
                        RES=res.get(gs)['rank_ES'], phenoPos=pheno_pos, phenoNeg=pheno_neg, figsize=figsize)        
        gs = gs.replace('/','_').replace(":","_")
        fig.savefig('{a}/{b}.gsea.{c}'.format(a=outdir, b=gs, c=format), bbox_inches='tight', dpi=300,)

   
    logger.info("Congratulations...GSEApy run successfully...............")
    
    
    #return dataframe if run gsea inside python console
    #if isinstance(rnk, pd.DataFrame):
    log_remove(logger)
    if hasattr(sys, 'ps1'):
        return res

def replot(indir, outdir='gseapy_replot', weight=1, figsize=[6.5,6], format='pdf', min_size=3, max_size=5000, verbose=False):
    """The main fuction to run inside python.
          
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
    argument = locals()

    mkdirs(outdir)   
    logger = log_init(outdir, module='replot', log_level= logging.INFO if verbose else logging.WARNING)
    #write command to log file
    argument = OrderedDict(sorted(argument.items(), key=lambda t:t[0]))
    logger.debug("Command: replot, "+str(argument))
    import glob
    from bs4 import BeautifulSoup
    
    #parsing files.......    
    try:
        results_path = glob.glob(indir+'*/edb/results.edb')[0]
        rank_path =  glob.glob(indir+'*/edb/*.rnk')[0]
        gene_set_path =  glob.glob(indir+'*/edb/gene_sets.gmt')[0]
    except IndexError as e: 
        logger.debug(e) 
        logger.error("Could not locate GSEA files in the given directory!")
        sys.exit(1)    
    #extract sample names from .cls file
    cls_path = glob.glob(indir+'*/edb/*.cls')
    if cls_path:
        phenoPos, phenoNeg, classes = gsea_cls_parser(cls_path[0])
    else:
        # logic for prerank results
        phenoPos, phenoNeg = '',''  
    #obtain gene sets
    gene_set_dict = gsea_gmt_parser(gene_set_path, min_size=min_size, max_size=max_size)
    #obtain rank_metrics
    rank_metric = gsea_rank_metric(rank_path)
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
        RES = enrichment_score(gene_list=gene_list, gene_set=gene_set, weighted_score_type=weight, 
                               correl_vector=correl_vector)[2]
        #plotting
        fig = gsea_plot(rank_metric, enrich_term,hit_ind, nes, pval,
                        fdr, RES, phenoPos, phenoNeg, figsize=figsize)    
        fig.savefig('{a}/{b}.gsea.replot.{c}'.format(a=outdir, b=enrich_term, c=format),
                    bbox_inches='tight', dpi=300,)

      
    logger.info("Congratulations! Your plots have been reproduced successfully!")
    log_remove(logger)

    return 