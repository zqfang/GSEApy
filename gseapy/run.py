# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 12:07:06 2016

@author: Bioninja
"""

from .parser import gsea_rank_metric,gsea_gmt_parser
from .algorithm import gsea_compute


def run(data, gene_sets, min_size= 15, max_size=1000, min_part=0.1,
     permutation_n=1000, permutation_type="gene_set", ranking_metric=None):
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
    
    if ranking_metric is not None:
        rankings = ranking_metric
        
    rankings = gsea_rank_metric("./gseapy/data/edb/gsea_data.gsea_data.rnk")
    gmt = gsea_gmt_parser("./gseapy/data/edb/gene_sets.gmt")
    
    
    #ES, RES = enrichment_score()
    results = gsea_compute(gene_list = rankings['gene_name'],rankings = rankings['rank'].values,
                    n=100,gmt = gmt, weighted_score_type=1,permutation_type='gene_set',expression_data=None)
    
    res = {}

    for gs, gseale in zip(gmt.keys(), list(results)):
        rdict = {}
        rdict['es'] = gseale[0]
        rdict['nes'] = gseale[1]
        rdict['p'] = gseale[2]
        rdict['fdr'] = gseale[3]
        rdict['size'] = len(gmt[gs])
        #rdict['matched_size'] = len(gs)
        #rdict['genes'] = subset
        res[gs] = rdict

        return res