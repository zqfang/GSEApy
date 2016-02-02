# -*- coding: utf-8 -*-

from __future__ import absolute_import, print_function
from bs4 import BeautifulSoup
from numpy import in1d
from pandas import read_table

import sys


def gsea_cls_parser(cls):
    """Extact class(phenotype) name from .cls file.
    
    :param cls: the .cls file where located inside edb folder.
    :return: phenotype name and a list of class vector. 
    """
    
    with open(cls) as cls:
        file = cls.readlines()
    sample_name = file[1].strip('\n').split(" ")        
    classes = file[2].strip('\n').split(" ")
    phenoPos = sample_name[1]
    phenoNeg = sample_name[2]
    
    return phenoPos,phenoNeg,classes



def gsea_edb_parser(results_path,index = 0,):
    """Parse results.edb file stored under **edb** file folder.            

    :param results_path: the .results file where lcoated inside edb folder.
    :param index: gene_set index of gmt database, used for iterating items.
   
    :return: enrichment_term, hit_index,nes, pval, fdr.
    """

    
    soup = BeautifulSoup(open(results_path),features='xml')
    tag = soup.findAll('DTG')
   
    term = dict(tag[index].attrs)
    # dict_keys(['RANKED_LIST', 'GENESET', 'FWER', 'ES_PROFILE', 
    # 'HIT_INDICES', 'ES', 'NES', 'TEMPLATE', 'RND_ES', 'RANK_SCORE_AT_ES',
    #'NP', 'RANK_AT_ES', 'FDR'])
      
    
    enrich_term = term.get('GENESET').split("#")[1]
    es_profile = term.get('ES_PROFILE').split(" ")
    #rank_es = term.get('RND_ES').split(" ")
    hit_ind =term.get('HIT_INDICES').split(" ")

    es_profile = [float(i) for i in es_profile ]
    hit_ind = [float(i) for i in hit_ind ]
    #rank_es = [float(i) for i in rank_es ]
    nes = term.get('NES')
    pval = term.get('NP')
    fdr =  term.get('FDR')
    #fwer = term.get('FWER')
   
   
    
    #index_range = len(tag)-1
    print("Enriched Gene set is: ", enrich_term)
    return enrich_term,hit_ind, nes,pval,fdr
    

def gsea_rank_metric(rnk):
    """Parse .rnk file. This file contains ranking correlation vector and gene names or ids. 
    
    :param rnk: the .rnk file where located inside the edb folder.
    :return: a pandas DataFrame with 3 columns names are::
             
                 'gene_name','rank',rank2'
                 
    """
    
    
    rank_metric = read_table(rnk,header=None)
    rank_metric.columns = ['gene_name','rank']
    rank_metric['rank2'] = rank_metric['rank']
     
    return rank_metric
    
def gsea_gmt_parser(gmt, min_size = 3, max_size = 5000, gene_list=None):
    """Parse gene_sets.gmt(gene set database) file. 
    
    :param gmt: the gene_sets.gmt file where loacated inside edb folder.
    :param min_size: Minimum allowed number of genes from gene set also the data set. Default: 3. 
    :param max_size: Maximum allowed number of genes from gene set also the data set. Default: 5000.
    :param gene_list: Used for filtering gene set. Only used this argument for :func:`run` method.
    :return: Return a new filtered gene set database dictionary. 

    **DO NOT** filter gene sets, when use :func:`replot`. Because ``GSEA`` Desktop have already
    do this for you.
            
    """
    
    
    with open(gmt) as genesets:
        genesets_dict = { line.rstrip("\n").split("\t")[0]:  
                          line.rstrip("\n").split("\t")[2:] 
                          for line in genesets.readlines()}
    
    
    #filtering dict
    if sys.version_info[0] == 3 :
        genesets_filter =  {k: v for k, v in genesets_dict.items() if len(v) >= min_size and len(v) <= max_size}
    elif sys.version_info[0] == 2:
        genesets_filter =  {k: v for k, v in genesets_dict.iteritems() if len(v) >= min_size and len(v) <= max_size}
    else:
        print("system failure. Please Provide correct input files")
        sys.exit(1)
    
    
    if gene_list is not None:
        subsets = sorted(genesets_filter.keys())             
        for subset in subsets:            
            tag_indicator = in1d(gene_list,genesets_filter.get(subset),assume_unique=True)
            tag_len = sum(tag_indicator)      
            if tag_len <= min_size and tag_len >= max_size:                    
                del genesets_filter[subset]
     #some_dict = {key: value for key, value in some_dict.items() if value != value_to_remove}
     #use np.intersect1d() may be faster???    
    filsets_num = len(genesets_dict) - len(genesets_filter)
    print("{a} gene_sets have been filtered out for max_size = {b} and min_size = {c}".format(a=filsets_num,b=max_size,c=min_size))
          
    return genesets_filter
    
