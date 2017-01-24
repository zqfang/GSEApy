#! python
# -*- coding: utf-8 -*-
from __future__ import  absolute_import, division


import os,sys, logging
from .parser import *
from .algorithm import enrichment_score, gsea_compute, preprocess, ranking_metric
from .plot import gsea_plot, heatmap
from collections import OrderedDict
from .utils import log_init, log_remove, mkdirs
import pandas as pd


def replot(indir, outdir='gseapy_out', weight=1, figsize=[6.5,6], format='pdf', min_size=3, max_size=5000):
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
    
    :return: Generate new figures with seleted figure format. Default: 'pdf'.   
    """
    argument = locals()

    mkdirs(outdir)   
    logger = log_init(outdir, module='replot')
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
        fig.savefig('{a}/.gsea.replot.{b}.{c}'.format(a=outdir, b=enrich_term, c=format),
                    bbox_inches='tight', dpi=300,)

      
    logger.info("Congratulations! Your plots have been reproduced successfully!")
    log_remove(logger)
    return 
        
def call(data, gene_sets, cls, outdir='gseapy_out', min_size=15, max_size=500, permutation_n=1000, 
          weighted_score_type=1,permutation_type='gene_set', method='log2_ratio_of_classes',
	  ascending=False, figsize=[6.5,6], format='pdf', graph_num=20, seed=None):
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
    logger = log_init(outdir, module='call')

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
    res = OrderedDict()
    for gs,gseale,ind,RES in zip(subsets, list(results), hit_ind, rank_ES):        
        rdict = OrderedDict()      
        rdict['es'] = gseale[0]
        rdict['nes'] = gseale[1]
        rdict['pval'] = gseale[2]
        rdict['fdr'] = gseale[3]
        rdict['gene_set_size'] = len(gmt[gs])
        rdict['matched_size'] = len(ind)
        rdict['rank_ES'] = RES
        rdict['genes'] = dat2.ix[ind,'gene_name'].tolist()
        rdict['hit_index'] = ind
        res[gs] = rdict           
    
    res_df = pd.DataFrame.from_dict(res,orient='index')
    res_df.index.name = 'Term'
    res_df.sort_values(by='fdr', inplace=True)
    
    res_df.drop(['rank_ES','hit_index'], axis=1, inplace=True)
    res_df.to_csv('{a}/{b}.{c}.gsea.reports.csv'.format(a=outdir, b='gseapy', c=permutation_type), float_format ='%.7f')
    


    #Plotting
    top_term = res_df.head(graph_num).index
    width = len(classes) if len(classes) >= 6 else  5
    for gs in top_term:
        hit = res.get(gs)['hit_index']
        gene_symbol = res.get(gs)['genes']
        fig = gsea_plot(rank_metric=dat2, enrich_term=gs, hit_ind=hit,
                        nes=res.get(gs)['nes'], pval=res.get(gs)['pval'], fdr=res.get(gs)['fdr'], 
                        RES=res.get(gs)['rank_ES'], phenoPos=phenoPos, phenoNeg=phenoNeg, figsize=figsize)        
        gs = gs.replace('/','_')
        fig.savefig('{a}/{b}.gsea.{c}'.format(a=outdir, b=gs, c=format), bbox_inches='tight', dpi=300,)

        heatmap(df=dat.loc[gene_symbol], term=gs, outdir=outdir, 
                figsize=(width, len(gene_symbol)/2), format=format)
      
    logger.info("Congratulations. GSEAPY run successfully................")
    
	# return dataframe if run gsea inside python console
    #if isinstance(data, pd.DataFrame) or isinstance(cls, list):
    log_remove(logger)
    if hasattr(sys, 'ps1'):
        return res_df 

def prerank(rnk, gene_sets, outdir='gseapy_out', pheno_pos='Pos', pheno_neg='Neg',
            min_size=15, max_size=500, permutation_n=1000, weighted_score_type=1,
            ascending=False, figsize=[6.5,6], format='pdf', graph_num=20, seed=None):
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
    logger = log_init(outdir, module='prerank')
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
    
    res = OrderedDict()
    for gs,gseale,ind,RES in zip(subsets, list(results), hit_ind, rank_ES):        
        rdict = OrderedDict()       
        rdict['es'] = gseale[0]
        rdict['nes'] = gseale[1]
        rdict['pval'] = gseale[2]
        rdict['fdr'] = gseale[3]
        rdict['gene_set_size'] = len(gmt[gs])
        rdict['matched_size'] = len(ind)
        rdict['rank_ES'] = RES
        rdict['genes'] = dat2.ix[ind,'gene_name'].tolist()
        rdict['hit_index'] = ind
        res[gs] = rdict           


    res_df = pd.DataFrame.from_dict(res, orient='index')
    res_df.index.name = 'Term'
    res_df.sort_values(by='fdr', inplace=True)
    
    res_df.drop(['rank_ES','hit_index'], axis=1, inplace=True)
    res_df.to_csv('{a}/{b}.prerank.reports.csv'.format(a=outdir, b='gseapy'), float_format ='%.7f')
    

    #Plotting
    top_term = res_df.head(graph_num).index
    for gs in top_term:
        fig = gsea_plot(rank_metric=dat2, enrich_term=gs, hit_ind=res.get(gs)['hit_index'],
                        nes=res.get(gs)['nes'], pval=res.get(gs)['pval'], fdr=res.get(gs)['fdr'], 
                        RES=res.get(gs)['rank_ES'], phenoPos=pheno_pos, phenoNeg=pheno_neg, figsize=figsize)        
        gs = gs.replace('/','_')
        fig.savefig('{a}/{b}.gsea.{c}'.format(a=outdir, b=gs, c=format), bbox_inches='tight', dpi=300,)

   
    logger.info("Congratulations...GSEAPY run successfully...............")
    
    
    #return dataframe if run gsea inside python console
    #if isinstance(rnk, pd.DataFrame):
    log_remove(logger)
    if hasattr(sys, 'ps1'):
        return res_df 
