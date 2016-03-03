#! python
# -*- coding: utf-8 -*-
import os
import sys
import time


from .parser import gsea_edb_parser,gsea_rank_metric,gsea_gmt_parser,gsea_cls_parser
from .algorithm import enrichment_score,gsea_compute,preprocess,ranking_metric
from .gsea_plot import gsea_plot
from collections import OrderedDict

import pandas as pd

def replot(indir,outdir='gseapy_out', weight=1,figsize=[6.5,6], format='pdf',min_size=3, max_size=5000):
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
    import glob
    from bs4 import BeautifulSoup   
    #parsing files.......    
    results_path = glob.glob(indir+'*/edb/results.edb')[0]
    rank_path =  glob.glob(indir+'*/edb/*.rnk')[0]
    gene_set_path =  glob.glob(indir+'*/edb/gene_sets.gmt')[0]
    cls_path = glob.glob(indir+'*/edb/*.cls')[0]
    file_list = [results_path, rank_path, gene_set_path, cls_path]      
    for file in file_list: 
        if not os.path.isfile(file):
            print("Incorrect Input %s !" %file)
            sys.exit(1)    
    #extract sample names from .cls file
    phenoPos, phenoNeg, classes = gsea_cls_parser(cls_path)  
    
    #extract each enriment term in the results.edb files and plot.
    database = BeautifulSoup(open(results_path),features='xml')
    length = len(database.findAll('DTG'))
    os.system("mkdir "+ outdir)
    for idx in range(length):
        #extract statistical resutls from results.edb file
        enrich_term, hit_ind, nes, pval, fdr= gsea_edb_parser(results_path, index=idx)
        
        #obtain rank_metrics
        rank_metric = gsea_rank_metric(rank_path)
        correl_vector =  rank_metric['rank'].values

        #obtain gene sets
        gene_set_dict = gsea_gmt_parser(gene_set_path, min_size=min_size, max_size=max_size)
        gene_set = gene_set_dict.get(enrich_term)
        gene_list = rank_metric['gene_name']

        #calculate enrichment score    
        RES = enrichment_score(gene_list=gene_list, gene_set=gene_set, weighted_score_type=weight, 
                               correl_vector=correl_vector)[2]

        #plotting
        fig = gsea_plot(rank_metric, enrich_term,hit_ind, nes, pval,
                        fdr, RES, phenoPos, phenoNeg, figsize=figsize)    
        fig.savefig('{a}/{b}.{c}'.format(a=outdir, b=enrich_term, c=format), dpi=300,)
        
    print("Congratulations! Your plots have been reproduced successfully!")

def call(data, gene_sets, cls, outdir='gseapy_out', min_size=15, max_size=1000, permutation_n=1000, weighted_score_type=1,
        permutation_type='gene_set', method='log2_ratio_of_classes', ascending=False, figsize=[6.5,6], format='pdf', 
        graph_num=20, seed=None):
    """ Run Gene Set Enrichment Analysis.

    :param data: Gene expression data table.  
    :param gene_sets: Gene sets file. e.g. gmt files. Same input with GSEA.
    :param permutation_n: Number of permutations for significance computation. Default: 1000.
    :param permutation_type: Permutation type, "phenotype" (default) for phenotypes, "gene_set" for genes.
    :param int min_size: Minimum allowed number of genes from gene set also the data set. Defaut: 15.
    :param int max_size: Maximum allowed number of genes from gene set also the data set. Defaults: 15 and 1000.
    :param weighted_score_type: refer to :func:`algorithm.enrichment_socre`. Default:1.
    :param method: ranking metric method, refer to :func:`algorithm.ranking_metric`.
    :param ascending: sorting order of rankings. Default: False.
    :param outdir: results output directory.
    :param figsize: matplotlib figsize, accept a tuple or list, e.g. [width,height]. Default: [6.5,6].
    :param format: matplotlib figure format. Default: pdf.
    :param graph_num: plot graphs for top sets of each phenotype
    :param seed: random seed. expect an interger. Defalut:None.
    :return: Generate ``GSEA`` plots and store a dictionary into csv file,
             where dictionary key is a gene set and values are::

                 | {es: enrichment score, 
                 |  nes: normalized enrichment score, 
                 |  p: P-value, 
                 |  fdr: FDR, 
                 |  size: gene set size,
                 |  matched_size: genes matched to the data, 
                 |  genes: gene names from the data set}
    
    """
    assert permutation_type in ["phenotype", "gene_set"]
    df = pd.read_table(data)    
    assert len(df) > 1   
    
    #select correct expression genes and values.
    dat = preprocess(df)
    
    # phenotype labels parsing
    phenoPos, phenoNeg, classes = gsea_cls_parser(cls)
    
    #ranking metrics calculation.    
    dat2 = ranking_metric(df=dat, method=method, phenoPos=phenoPos, phenoNeg=phenoNeg, classes=classes, ascending=ascending)
    
    #filtering out gene sets and build gene sets dictionary
    gmt = gsea_gmt_parser(gene_sets, min_size=min_size, max_size=max_size,gene_list=dat2['gene_name'].values)
    
    #compute ES, NES, pval, FDR, RES
    results,hit_ind,rank_ES, subsets = gsea_compute(data=dat, n=permutation_n,gmt=gmt, weighted_score_type=weighted_score_type,
                    permutation_type=permutation_type, method=method, phenoPos=phenoPos, phenoNeg=phenoNeg, classes=classes, ascending=ascending, seed=seed)
   
    print("Start to generate gseapy reports, and produce figures...",time.ctime())
    os.system("mkdir "+ outdir)
    res = OrderedDict()
    for gs, gseale,ind,RES in zip(subsets, list(results), hit_ind, rank_ES):        
        rdict = OrderedDict()      
        rdict['es'] = gseale[0]
        rdict['nes'] = gseale[1]
        rdict['pval'] = gseale[2]
        rdict['fdr'] = gseale[3]
        rdict['gene_set_size'] = len(gmt[gs])
        rdict['matched_size'] = len(ind)
        rdict['rank_ES'] = RES
        rdict['genes'] = dat.iloc[ind].index.tolist()
        res[gs] = rdict           

        #plotting
        fig = gsea_plot(rank_metric=dat2, enrich_term=gs, hit_ind=ind, nes=gseale[1], pval=gseale[2],
                        fdr=gseale[3], RES=RES, phenoPos=phenoPos, phenoNeg=phenoNeg, figsize=figsize)        
        fig.savefig('{a}/{b}.{c}'.format(a=outdir, b=gs, c=format), dpi=300,)

    res_df =pd.DataFrame.from_dict(res,orient='index')
    res_df.index.name = 'Enrich_terms'
    #res_df = res_df[['es','nes','pval','fdr','gene_set_size','matched_size','rank_ES','genes']]
    res_df.sort_values(by='fdr', inplace=True)
    res_final = res_df.head(graph_num)
    res_final.to_csv('{a}/{b}.csv'.format(a=outdir, b='gseapy_reports'), float_format ='%.7f')
    
    #print(res_df.head(10))
    print("...Congratulations. GSEAPY run successfully!!!\n...The Job is done...........................Goodbye!")
    
    return res_df

def prerank(rnk, gene_sets, outdir='gseapy_out', pheno_pos='Postive', pheno_neg='Negative',
            min_size=15, max_size=1000, permutation_n=1000, weighted_score_type=1,
            ascending=False, figsize=[6.5,6], format='pdf',graph_num=20, seed=None):
    """ Run Gene Set Enrichment Analysis with pre-ranked correlation defined by user.

    :param rnk: pre-ranked correlation table, Same input with ``GSEA`` .rnk file.  
    :param gene_sets: Gene sets file. e.g. gmt files. Same input with GSEA.
    :param outdir: results output directory.
    :param permutation_n: Number of permutations for significance computation. Default: 1000.
    :param int min_size: Minimum allowed number of genes from gene set also the data set. Defaut: 15.
    :param int max_size: Maximum allowed number of genes from gene set also the data set. Defaults: 15 and 1000.
    :param weighted_score_type: refer to :func:`algorithm.enrichment_socre`. Default:1.
    :param ascending: sorting order of rankings. Default: False.
    :param figsize: matplotlib figsize, accept a tuple or list, e.g. [width,height]. Default: [6.5,6].
    :param format: matplotlib figure format. Default: pdf.
    :param graph_num: plot graphs for top sets of each phenotype
    :param seed: random seed. expect an interger. Defalut:None.    
    :return: Generate ``GSEA`` plots and store a dictionary into csv file,
             where dictionary key is a gene set and values are::

                 | {es: enrichment score, 
                 |  nes: normalized enrichment score, 
                 |  p: P-value, 
                 |  fdr: FDR, 
                 |  size: gene set size,
                 |  matched_size: genes matched to the data, 
                 |  genes: gene names from the data set}
    
    """
    
    dat2 = gsea_rank_metric(rnk)    
    assert len(dat2) > 1   
       
    #ranking metrics calculation.    
    
    #filtering out gene sets and build gene sets dictionary
    gmt = gsea_gmt_parser(gene_sets, min_size=min_size, max_size=max_size, gene_list=dat2['gene_name'].values)
    
    #compute ES, NES, pval, FDR, RES
    results,hit_ind,rank_ES, subsets = gsea_compute(data=dat2, n=permutation_n, gmt=gmt, weighted_score_type=weighted_score_type,
                                                    permutation_type='gene_set',method=None, phenoPos=pheno_pos, phenoNeg=pheno_neg,
                                                    classes=None, ascending=ascending,seed=seed, prerank=True)
   
    print("Start to generate gseapy reports, and produce figures...",time.ctime())
    os.system("mkdir "+ outdir)
    res = OrderedDict()
    for gs, gseale,ind,RES in zip(subsets, list(results), hit_ind, rank_ES):        
        rdict = OrderedDict()       
        rdict['es'] = gseale[0]
        rdict['nes'] = gseale[1]
        rdict['pval'] = gseale[2]
        rdict['fdr'] = gseale[3]
        rdict['gene_set_size'] = len(gmt[gs])
        rdict['matched_size'] = len(ind)
        rdict['rank_ES'] = RES
        rdict['genes'] = dat2.ix[ind,'gene_name'].tolist()
        res[gs] = rdict           

        #plotting
        fig = gsea_plot(rank_metric=dat2, enrich_term=gs,hit_ind=ind,nes=gseale[1],pval= gseale[2],
                        fdr=gseale[3], RES=RES, phenoPos=pheno_pos, phenoNeg=pheno_neg, figsize=figsize)        
        fig.savefig('{a}/{b}.{c}'.format(a=outdir, b=gs, c=format), dpi=300,)

    res_df = pd.DataFrame.from_dict(res, orient='index')
    res_df.index.name = 'Enrich_terms'
    res_df.sort_values(by='fdr', inplace=True)
    #res_df = res_df[['es','nes','pval','fdr','gene_set_size','matched_size','rank_ES','genes']]
    res_final = res_df.head(graph_num)
    res_final.to_csv('{a}/{b}.csv'.format(a=outdir, b='gseapy_reports'), float_format ='%.7f')
    #print(res_df.head(10))
    print("...Congratulations. GSEAPY run successfully!!!\n...The Job is done...........................Goodbye!")
    
    return res_df
