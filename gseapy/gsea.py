#! python
# -*- coding: utf-8 -*-
import os
import sys



from bs4 import BeautifulSoup
from .parser import gsea_edb_parser,gsea_rank_metric,gsea_gmt_parser,gsea_cls_parser
from .algorithm import enrichment_score,gsea_compute,preprocess
from .gsea_plot import gsea_plot

import glob
import pandas as pd

def re_plot(indir,outdir,weight=1,figsize=[6.5,6],format='pdf',):
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
    phenoPos,phenoNeg = gsea_cls_parser(cls_path)  
    
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
                               correl_vector = correl_vector)



        #plotting
        fig = gsea_plot(rank_metric, enrich_term,hit_ind,nes,pval,
                        fdr, RES, phenoPos,phenoNeg,figsize=figsize)
    
        fig.savefig('{a}/{b}.{c}'.format(a=outdir,b=enrich_term,c=format),dpi=300,)


def run(data, gene_sets,cls, min_size, max_size, permutation_n, weighted_score_type,
        permutation_type, method,ascending, out,figsize=[6.5,6]):
    """ Run Gene Set Enrichment Analysis.

    :param data.Table data: Gene expression data.  
    :paramg GSEA gene_sets: Gene sets. e.g. gmt files  
    
.
    :param permutation_n: Number of permutations for significance computation. Default: 100.
    :param str permutation_type: Permutation type, "phenotype" (default) for 
        phenotypes, "gene_set" for genes.
    :param int min_size:
    :param int max_size: Minimum and maximum allowed number of genes from
        gene set also the data set. Defaults: 15 and 1000.

    :param float min_part: Minimum fraction of genes from the gene set
        also in the data set. Default: 0.1.

    :param permutation_type: 
    
    :param ranking_metric:

    :return: | a dictionary where key is a gene set and values are:
        | { es: enrichment score, 
        | nes: normalized enrichment score, 
        | p: P-value, 
        | fdr: FDR, 
        | size: gene set size,
        | matched_size: genes matched to the data, 
        | genes: gene names from the data set }

    """
    assert len(data) > 1
    assert permutation_type in ["phenotype", "gene_set"]
    
    data = pd.read_table(data,index_col=0)
    phenoA, phenoB, classes = gsea_cls_parser(cls)
    gmt = gsea_gmt_parser(gene_sets, min_size = min_size, max_size = max_size)
    #gmt.sort()

    dat = preprocess(data)
    #dat2 = ranking_metric(dat,method= method,classes = classes ,ascending=ascending)
    #compute ES, NES, pval, FDR, RES
    results,hit_ind,rank_ES, subsets = gsea_compute(data = dat, n=permutation_n,gmt = gmt, weighted_score_type=weighted_score_type,
                    permutation_type=permutation_type,method=method,classes = classes,ascending=ascending)
   
   
    res = {}
    for gs, gseale,ind,rank in zip(subsets,list(results),hit_ind,rank_ES):        
        rdict ={}       
        rdict['es'] = gseale[0]
        rdict['nes'] = gseale[1]
        rdict['pval'] = gseale[2]
        rdict['fdr'] = gseale[3]
        rdict['size'] = len(gmt[gs])
        rdict['matched_size'] = len(ind)
        rdict['rank_ES'] = rank
        rdict['genes'] = dat.iloc[ind].index.tolist()
        res[gs] = rdict
    
    
    
        #plotting

        #fig = gsea_plot(rank_metric = dat, enrich_term = gs,hit_ind = ind,nes = gseale[1],pval= gseale[2],
        #                fdr = gseale[3], RES = rnk, phenoPos =phenoA ,phenoNeg = phenoB,figsize=figsize)
        
        #fig.savefig('{a}/{b}.{c}'.format(a= out,b= gs,c= 'pdf'),dpi=300,)
    
    
    
    return res
if __name__ == '__main__':
    run()

