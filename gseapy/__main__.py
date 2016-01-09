#! python
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import matplotlib as mpl







from bs4 import BeautifulSoup
from .parser import gsea_edb_parser,gsea_rank_metric,gsea_gmt_parser,gsea_cls_parser

from .algorithm import enrichment_score
from .gsea_plot import gsea_plot

import glob


__version__ = '0.2.2'
__author__ = 'Zhuoqing Fang'

def main():
    """The main routine."""
    

    # parse command line args
    parser = argparse.ArgumentParser(description="Python wrapper of Gene Set Enrichment Analysis tool")
    parser.add_argument("-i","--InDir", action="store", dest="file",
                        help="the GSEA desktop output dir that you want to reproduce the figure ")
    parser.add_argument("-o","--outDir",action="store",default="foo",dest="out",\
                     help="the output directory")
    parser.add_argument("--version",action="version",version="%(prog)s {}".format(__version__))
    
    args = parser.parse_args()

    print("Input_directroy        =", args.file)
    print("Output_directory      =", args.out)
    
    

    file_name = args.file
    
    # checking flies and parameters.
    if not os.path.exists(args.file) :
        print("Input_Directory doesn't exist, please check your file path!")
        sys.exit(1)    
    
   
    
    print("parsing files.......")
    
    results_path = glob.glob(file_name+'*/edb/results.edb')[0]
    rank_path =  glob.glob(file_name+'*/edb/*.rnk')[0]
    gene_set_path =  glob.glob(file_name+'*/edb/gene_sets.gmt')[0]
    cls_path = glob.glob(file_name+'*/edb/*.cls')[0]
    file_list = [results_path ,rank_path,gene_set_path,cls_path]  
    
    for file in file_list: 
        if not os.path.isfile(file):
            print("Incorrect Input %s !" %file)
            sys.exit(1)

    os.system('mkdir '+ args.out)
   

    
    
    #extract sample names from .cls file
    phenoPos,phenoNeg = gsea_cls_parser(cls_path)  
    
    #extract each enriment term in the results.edb files and plot.
    database = BeautifulSoup(open(results_path),features='xml')
    length = len(database.findAll('DTG'))
    
    
    print("Generate Plots......Please wait.....")
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
        RES = enrichment_score(gene_list = gene_list, gene_set = gene_set, weighted_score_type = 1, 
                               correl_vector = correl_vector)

        #plotting
        fig = gsea_plot(rank_metric, enrich_term,es_profile,hit_ind,nes,pval,fdr,
                        RES, phenoPos,phenoNeg,figsize=(6.5,6))
        fig.savefig(args.out+'/'+enrich_term+'.png',format='png',dpi=300,)
    
    print("Congratulations! The job is done!")

if __name__ == "__main__":
    #do not show the figure
    mpl.use('Agg')
    main()
