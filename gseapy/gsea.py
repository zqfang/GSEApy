#! python
# -*- coding: utf-8 -*-
import os
import sys



from bs4 import BeautifulSoup
from .parser import gsea_edb_parser,gsea_rank_metric,gsea_gmt_parser,gsea_cls_parser

from .algorithm import enrichment_score
from .gsea_plot import gsea_plot

import glob



def gsea(indir,outdir,weight=1,figsize=[6.5,6],format='pdf',):
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
        enrich_term,es_profile,hit_ind, nes,pval,fdr,rank_es = gsea_edb_parser( results_path,index=idx)

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
        fig = gsea_plot(rank_metric, enrich_term,es_profile,hit_ind,nes,pval,
                        fdr, RES, phenoPos,phenoNeg,figsize=figsize)
    
        fig.savefig(outdir+'/'+enrich_term+'./'+format,dpi=300,)

if __name__ == "__main__":
    gsea()
