# -*- coding: utf-8 -*-
import numpy as np

def enrichment_score(gene_list, gene_set, weighted_score_type = 1, correl_vector = None):
    '''
    this is the most important function of GSEApy. It has the same algorithm with GSEA.
    
    paramter
    -----------------
    correl_vector:                  rank_metric['rank'].values
    gene_list:                      rank_metric['gene_name']
    gene_set:                       gene_sets in gmt file, please used gsea_gmt_parser to get gene_set.
    weighted_sore_type:             0, 1, 1.5, 2. It's indentical to gsea's weighted_sore_method. defalut: 1.
    
    '''
    
    #Test whether each element of a 1-D array is also present in a second array
    #use astype covert bool to intergers  
    tag_indicator = np.in1d(gene_list,gene_set).astype(int)  # notice that the sign is 0 (no tag) or 1 (tag)
    no_tag_indicator = 1 - tag_indicator
     
    
    #compute ES score, the code below is identical to gsea enrichment_score method.
    N = len(gene_list) 
    Nh = len(gene_set) 
    Nm =  N - Nh 
    if (weighted_score_type == 0 ): 
        correl_vector = np.repeat(1, N)
  
    alpha = weighted_score_type
    correl_vector = np.abs(correl_vector**alpha)
    sum_correl_tag = np.sum(correl_vector[tag_indicator.astype(bool)])
    norm_tag =  1.0/sum_correl_tag
    norm_no_tag = 1.0/Nm
    RES = np.cumsum(tag_indicator * correl_vector * norm_tag - no_tag_indicator * norm_no_tag)      
    max_ES = max(RES)
    min_ES = min(RES)
    
    print("Max ES is ",max_ES)
    print("Min ES is ",min_ES)
    print("length ES is", len(RES))
    
    return RES

