# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 12:07:06 2016

@author: Bioninja
"""

from .parser import gsea_rank_metric,gsea_gmt_parser,gsea_cls_parser
from .algorithm import gsea_compute,ranking_metric

import pandas as pd


def run(data, gene_sets,cls, min_size= 15, max_size=1000, permutation_n=1000, weighted_score_type=1,
        permutation_type="gene_set", method='log2_ratio_of_classes',ascending=False,rank_metric=None):
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
    
    data = pd.read_table(data)
    classes = gsea_cls_parser(cls)[2]
    gmt = gsea_gmt_parser(gene_sets)
    gmt.sort()
    #Ecompute ES, NES, pval, FDR, RES
    if rank_metric is None:
        dat = ranking_metric(data,method= method,classes = classes ,ascending=ascending)
        results,hit_ind,RES = gsea_compute(data = dat, gene_list = None,rankings = None,
                    n=permutation_n,gmt = gmt, weighted_score_type=weighted_score_type,
                    permutation_type=permutation_type)
    else:
        dat = pd.read_table(rank_metric)
        results,hit_ind,RES = gsea_compute(data = None, gene_list = rank_metric['gene_name'],rankings = rank_metric['rank'].values,
                                           n=permutation_n,gmt = gmt, weighted_score_type=weighted_score_type,
                                           permutation_type=permutation_type)
    
    res = {}

    for gs, gseale in zip(gmt.keys(), list(results)):
        rdict = {}
        rdict['es'] = gseale[0]
        rdict['nes'] = gseale[1]
        rdict['pval'] = gseale[2]
        rdict['fdr'] = gseale[3]
        rdict['size'] = len(gmt[gs])
        #rdict['matched_size'] = len(gseale[5])
        #rdict['genes'] = rankings.ix[gseale[5],'gene_name']
        res[gs] = rdict

    return res, hit_ind, RES