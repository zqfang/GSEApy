# -*- coding: utf-8 -*-


from bs4 import BeautifulSoup
import pandas as pd


def gsea_cls_parser(cls_path):
    '''
    extact sample name
    '''
    
    with open(cls_path) as cls:
        sample_name = cls.readlines()[1].strip('\n').split(" ")
    phenoPos = sample_name[1]
    phenoNeg = sample_name[2]
    return phenoPos,phenoNeg



def gsea_edb_parser(results_path,index = 0,):
    '''
    paramter
    -----------
    results_path:
    index:
   
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
    return enrich_term,es_profile,hit_ind, nes,pval,fdr,rank_es
    

def gsea_rank_metric(rank_path):
    '''
    parser rank_metric file
    '''
    
    
    rank_metric = pd.read_table(rank_path,header=None)
    rank_metric.columns = ['gene_name','rank']
    rank_metric['rank2'] = rank_metric['rank']
     
    return rank_metric
    
def gsea_gmt_parser(gene_set_path):
    '''
    parser gene set file
    '''
    
    
    with open(gene_set_path) as gene_sets:
        gene_set_dict = { line.rstrip("\n").split("\t")[0]:  
                          line.rstrip("\n").split("\t")[2:] 
                          for line in gene_sets.readlines()}
    return gene_set_dict
    
