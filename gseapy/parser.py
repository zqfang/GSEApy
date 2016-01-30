# -*- coding: utf-8 -*-

from __future__ import absolute_import, print_function
from bs4 import BeautifulSoup
import pandas as pd
import numpy as np
import sys


def gsea_cls_parser(cls_path):
    '''
    extact class(sample) name
    :param cls_path: location of GSEA-P .cls file
    '''
    
    with open(cls_path) as cls:
        file = cls.readlines()
    sample_name = file[1].strip('\n').split(" ")        
    classes = file[2].strip('\n').split(" ")
    phenoPos = sample_name[1]
    phenoNeg = sample_name[2]
    
    return phenoPos,phenoNeg,classes



def gsea_edb_parser(results_path,index = 0,):
    '''
    parse results.edb files.
            

    
    :param results_path: location of GSEA-P results.gmt file
    :param index: gene_sets_length.
   
    '''

    
    soup = BeautifulSoup(open(results_path),features='xml')
    tag = soup.findAll('DTG')
   
    term = dict(tag[index].attrs)
    # dict_keys(['RANKED_LIST', 'GENESET', 'FWER', 'ES_PROFILE', 
    # 'HIT_INDICES', 'ES', 'NES', 'TEMPLATE', 'RND_ES', 'RANK_SCORE_AT_ES',
    #'NP', 'RANK_AT_ES', 'FDR'])
      
    
    enrich_term = term.get('GENESET').split("#")[1]
    es_profile = term.get('ES_PROFILE').split(" ")
    rank_es = term.get('RND_ES').split(" ")
    hit_ind =term.get('HIT_INDICES').split(" ")

    es_profile = [float(i) for i in es_profile ]
    hit_ind = [float(i) for i in hit_ind ]
    rank_es = [float(i) for i in rank_es ]
    nes = term.get('NES')
    pval = term.get('NP')
    fdr =  term.get('FDR')
    #fwer = term.get('FWER')
   
   
    
    #index_range = len(tag)-1
    print("Enriched Gene set is: ", enrich_term)
    return enrich_term,hit_ind, nes,pval,fdr,rank_es
    

def gsea_rank_metric(rank_path):
    '''
    parse rank_metric file
    
    :param rank_path: location of GSEA-P .rnk file
    '''
    
    
    rank_metric = pd.read_table(rank_path,header=None)
    rank_metric.columns = ['gene_name','rank']
    rank_metric['rank2'] = rank_metric['rank']
     
    return rank_metric
    
def gsea_gmt_parser(gmt_path, min_size = 3, max_size = 1000, gene_list=None):
    '''
    parser gene sets file
    
    :param gmt_path: location of GSEA-P .gmt file
    '''
    
    
    with open(gmt_path) as gmt:
        genesets_dict = { line.rstrip("\n").split("\t")[0]:  
                          line.rstrip("\n").split("\t")[2:] 
                          for line in gmt.readlines()}
    
    
    #filtering dict
    if sys.version_info[0] == 3 :
        genesets_filter =  {k: v for k, v in genesets_dict.items() if len(v) >= min_size and len(v) <= max_size}
    elif sys.version_info[0] == 2:
        genesets_filter =  {k: v for k, v in genesets_dict.iteritems() if len(v) >= min_size and len(v) <= max_size}
    else:
        print("system failure. Please Provide correct input files")
        sys.exit(1)
    
    
    if gene_list != None:
        subsets = sorted(genesets_filter.keys())  
        keys_new = []    
        for subset in subsets:            
            tag_indicator = np.in1d(gene_list,genesets_filter.get(subset),assume_unique=True)
            tag_len = np.sum(tag_indicator)      
            if tag_len <= min_size and tag_len >= max_size:    
                keys_new.append(subset)
                del genesets_filter[subset]
     #some_dict = {key: value for key, value in some_dict.items() if value != value_to_remove}
    print("{a} gene_sets have been filtered out for max_size = {b} and min_size = {c}".format(a=len(keys_new),
          b=max_size,c=min_size))
          
    return genesets_filter
    
