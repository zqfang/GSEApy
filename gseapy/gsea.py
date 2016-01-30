#! python
# -*- coding: utf-8 -*-
import os
import sys



from bs4 import BeautifulSoup
from .parser import gsea_edb_parser,gsea_rank_metric,gsea_gmt_parser,gsea_cls_parser
from .algorithm import enrichment_score,gsea_compute,preprocess,ranking_metric
from .gsea_plot import gsea_plot

import glob
import pandas as pd

def replot(indir,outdir,weight=1,figsize=[6.5,6],format='pdf',):
    """The main fuction to run inside python."""
        
    #parsing files.......
    
    results_path = glob.glob(indir+'*/edb/results.edb')[0]
    rank_path =  glob.glob(indir+'*/edb/*.rnk')[0]
    gene_set_path =  glob.glob(indir+'*/edb/gene_sets.gmt')[0]
    cls_path = glob.glob(indir+'*/edb/*.cls')[0]
    file_list = [results_path ,rank_path,gene_set_path,cls_path]  
    
    for file in file_list: 
        if not os.path.isfile(file):
            print("Incorrect Input %s !" %file)
            sys.exit(1)


   
    
    
    #extract sample names from .cls file
    phenoPos,phenoNeg,classes = gsea_cls_parser(cls_path)  
    
    #extract each enriment term in the results.edb files and plot.
    database = BeautifulSoup(open(results_path),features='xml')
    length = len(database.findAll('DTG'))
    os.system("mkdir "+ outdir)
    for idx in range(length):
        #extract statistical resutls from results.edb file
        enrich_term,hit_ind, nes,pval,fdr,rank_es = gsea_edb_parser( results_path,index=idx)
        
        #obtain rank_metrics
        rank_metric = gsea_rank_metric(rank_path)
        correl_vector =  rank_metric['rank'].values

        #obtain gene sets
        gene_set_dict = gsea_gmt_parser(gene_set_path)
        gene_set = gene_set_dict.get(enrich_term)
        gene_list = rank_metric['gene_name']

        #calculate enrichment score    
        RES = enrichment_score(gene_list = gene_list, gene_set = gene_set, weighted_score_type = weight, 
                               correl_vector = correl_vector)[2]



        #plotting
        fig = gsea_plot(rank_metric, enrich_term,hit_ind,nes,pval,
                        fdr, RES, phenoPos,phenoNeg,figsize=figsize)
    
        fig.savefig('{a}/{b}.{c}'.format(a=outdir,b=enrich_term,c=format),dpi=300,)


def run(data, gene_sets,cls, min_size, max_size, permutation_n, weighted_score_type,
        permutation_type, method,ascending, outdir,figsize,format):
    """ Run Gene Set Enrichment Analysis.

    :param data: Gene expression data table.  
    :paramg gene_sets: Gene sets file. e.g. gmt files. Same input with GSEA.
    :param permutation_n: Number of permutations for significance computation. Default: 1000.
    :param permutation_type: Permutation type, "phenotype" (default) for phenotypes, "gene_set" for genes.
    :param int min_size:
    :param int max_size: Minimum and maximum allowed number of genes from gene set also the data set. 
                         Defaults: 15 and 1000.
    :param weighted_score_type: default:1
    :param ascending: sorting order of rankings. Default: False.
    :param outdir: results output directory.
    :param figsize: matplotlib figsize, accept a tuple or list, e.g. [width,height]. Default: [6.5,6].
    
    :return: | a dictionary where key is a gene set and values are:
        | { es: enrichment score, 
        | nes: normalized enrichment score, 
        | p: P-value, 
        | fdr: FDR, 
        | size: gene set size,
        | matched_size: genes matched to the data, 
        | genes: gene names from the data set }

    """
    assert permutation_type in ["phenotype", "gene_set"]
    df = pd.read_table(data)
    
    assert len(df) > 1   
    assert permutation_type in ["phenotype", "gene_set"]
    #select correct expression genes and values.
    dat = preprocess(df)
    # phenotype labels parsing
    phenoPos, phenoNeg, classes = gsea_cls_parser(cls)
    #ranking metrics calculation.    
    dat2 = ranking_metric(df = dat,method= method,phenoPos=phenoPos,phenoNeg=phenoNeg,classes = classes ,ascending=ascending)
    #filtering out gene sets and build gene sets dictionary
    gmt = gsea_gmt_parser(gene_sets, min_size = min_size, max_size = max_size,gene_list=dat2['gene_name'].values)
    
    #compute ES, NES, pval, FDR, RES
    results,hit_ind,rank_ES, subsets = gsea_compute(data = dat, n=permutation_n,gmt = gmt, weighted_score_type=weighted_score_type,
                    permutation_type=permutation_type,method=method,phenoPos=phenoPos,phenoNeg=phenoNeg,classes = classes,ascending=ascending)
   
   
    res = {}
    for gs, gseale,ind,RES in zip(subsets,list(results),hit_ind,rank_ES):        
        rdict ={}       
        rdict['es'] = gseale[0]
        rdict['nes'] = gseale[1]
        rdict['pval'] = gseale[2]
        rdict['fdr'] = gseale[3]
        rdict['size'] = len(gmt[gs])
        rdict['matched_size'] = len(ind)
        rdict['rank_ES'] = RES
        rdict['genes'] = dat.iloc[ind].index.tolist()
        res[gs] = rdict           
        #plotting
        
    
        fig = gsea_plot(rank_metric = dat2, enrich_term = gs,hit_ind = ind,nes = gseale[1],pval= gseale[2],
                        fdr = gseale[3], RES = RES, phenoPos =phenoPos ,phenoNeg = phenoNeg,figsize=figsize)
        
        fig.savefig('{a}/{b}.{c}'.format(a= outdir,b= gs,c= format),dpi=300,)

    res_df =pd.DataFrame.from_dict(res,orient='index')
    res_df.index.name = 'Enrich_terms'
    res_df.to_csv('{a}/{b}.csv'.format(a= outdir,b='gseapy_reports' ))
    
    
    return res

