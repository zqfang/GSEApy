# -*- coding: utf-8 -*-


from __future__ import absolute_import, print_function,division
from functools import reduce

import time
import numpy as np
import pandas as pd
import random





def enrichment_score(gene_list, gene_set, weighted_score_type = 1, correl_vector = None):
    '''
    this is the most important function of GSEApy. It has the same algorithm with GSEA.
    
    :param gene_list:       The ordered gene list 
                            gene_name_list, rank_metric['gene_name']

    :param gene_set:        gene_sets in gmt file, please used gsea_gmt_parser to get gene_set.
                             
    :param weighted_score_type: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list.
                                 It's indentical to gsea's weighted_sore_method. options: 0(classic),1,1.5,2. default:1.
    :param correl_vector:   A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list.
                            or rankings, rank_metric['rank'].values
    
    
    :return: 
    |ES: Enrichment score (real number between -1 and +1) 
    |RES: Numerical vector containing the running enrichment score for all locations in the gene list .
             
    '''
    
    #Test whether each element of a 1-D array is also present in a second array
    #It's more intuitived here than orginal enrichment_score source code.
    #use .astype to covert bool to intergers  
    tag_indicator = np.in1d(gene_list,gene_set).astype(int)  # notice that the sign is 0 (no tag) or 1 (tag)
    no_tag_indicator = 1 - tag_indicator
     
    
    #compute ES score, the code below is identical to gsea enrichment_score method.
    N = len(gene_list) 
    Nhint = np.sum(tag_indicator)
    Nmiss =  N - Nhint 
    if (weighted_score_type == 0 ): 
        correl_vector = np.repeat(1, N)
  
    alpha = weighted_score_type
    correl_vector = np.abs(correl_vector**alpha)
    sum_correl_tag = np.sum(correl_vector[tag_indicator.astype(bool)])
    norm_tag =  1.0/sum_correl_tag
    norm_no_tag = 1.0/Nmiss
    RES = np.cumsum(tag_indicator * correl_vector * norm_tag - no_tag_indicator * norm_no_tag)      
    max_ES = max(RES)
    min_ES = min(RES)
    
    #print("The length of ES is", len(RES))
    
    return (max_ES if abs(max_ES) > abs(min_ES) else min_ES, RES)
 

def shuffle_list(gene_list, rand=random.Random(0)):
    """
    Returns a copy of a shuffled input gene_list.
    
    :gene_list: rank_metric['gene_name'].values
    """
    
    l2 = gene_list.values.copy()
    rand.shuffle(l2)

    return l2
        
    
def shuffle_class(df,n):
    """
    a function to shuffle rows and columns
    
    :param df:
    :return: dataframe
    """
    
    df2 = df.T.copy()
    for _ in range(n):    
        df2.apply(np.random.shuffle,axis=0)
    return df2.T       
    
    
    
def ranking_metric(df, method='log2_ratio_of_classes'):
    """
    :param df: gene_expression DataFrame.
    
    
    
    :param method:
     1. signal_to_noise 
         You must have at least three samples for each phenotype to use this metric.
         The larger the signal-to-noise ratio, the larger the differences of the means (scaled by the standard deviations);
         that is, the more distinct the gene expression is in each phenotype and the more the gene acts as a “class marker.” 
     2. t-test
         uses the difference of means scaled by the standard deviation and number of samples. 
         Note: You must have at least three samples for each phenotype to use this metric.
         The larger the tTest ratio, the more distinct the gene expression is in each phenotype 
         and the more the gene acts as a “class marker.”
     3. ratio_of_classes (also referred to as fold change) 
         uses the ratio of class means to calculate fold change for natural scale data.
     4. Diff_of_classes 
         uses the difference of class means to calculate fold change for log scale data
     5. log2_ratio_of_classes 
         uses the log2 ratio of class means to calculate fold change for natural scale data.
         This is the recommended statistic for calculating fold change for natural scale data.
          
     
    :return: returns correlation to class of each variable.
             same format with .rnk file. gene_name in first coloum, correlation
             in second column.
    """ 
    
    #To be complete
    
def gsea_pval(es, esnull):
    """
    From article (PNAS):
    estimate nominal p-value for S from esnull by using the positive
    or negative portion of the distribution corresponding to the sign 
    of the observed ES(S).
    """
    
    try:
        if es < 0:
            return float(len([ a for a in esnull if a <= es ]))/ \
                len([ a for a in esnull if a < 0])    
        else: 
            return float(len([ a for a in esnull if a >= es ]))/ \
                len([ a for a in esnull if a >= 0])
    except:
        return 1.0




def ordered_pointers_corr(correl_vector ):
    """
    Return a list of integers: indexes in original
    lcor. Elements in the list are ordered by
    their lcor[i] value. Higher correlations first.
    """
    ordered = [ (i,a) for i,a in enumerate(correl_vector) ] #original pos + correlation
    ordered.sort(key=lambda x: -x[1]) #sort by correlation, descending
    index = [ i[0] for i in ordered] #contains positions in the original list
    
    return index
    
    
def gsea_compute(gene_list, rankings, gmt, n, weighted_score_type=1,permutation_type='gene_set',expression_data=None):
    """
    compute enrichment scores and enrichment nulls. 
    
    :param gene_list: rank_metric['gene_name']
    :param rankings: correl_vector, ranking_metric['rank']
    :param subsets: all gene sets in .gmt file. need gmt_parser() results 
    :param n: permutation number. default: 1000
    :param weighted_score_type: default:1
    """
    enrichment_scores = []
    w = weighted_score_type
    subsets = gmt.keys()
    for subset in subsets:
        es = enrichment_score(gene_list = gene_list, gene_set=gmt.get(subset), 
                              weighted_score_type = w, correl_vector = rankings)[0]
        enrichment_scores.append(es)
    
    

    enrichment_nulls = [ [] for a in range(len(subsets)) ]
    
    index = ordered_pointers_corr(rankings)
    for i in range(n):
        if permutation_type == "phenotype":
            d2 = shuffle_class(expression_data,n=2000+i) #fixed permutation
            r2 = ranking_metric(d2)
            ranking2=r2['rank']
            gene_list2=r2['gene_name']
        else:
           
            index2 = shuffle_list(index, random.Random(2000+i))        
            gene_list2 = gene_list[index2]
            ranking2 = rankings[index2]
       
        
        for si,subset in enumerate(subsets):
            esn = enrichment_score(gene_list = gene_list2, gene_set=gmt.get(subset), 
                              weighted_score_type = w, correl_vector = ranking2)[0]            
            
            enrichment_nulls[si].append(esn)

        

    return gsea_significance(enrichment_scores, enrichment_nulls)


def gsea_significance(enrichment_scores, enrichment_nulls):
    """
    Computing p-vals, normalized ES, FDR
    """
    

    

    tb1 = time.time()
    print("Start to compute enrichment socres..........................",tb1)
    enrichmentPVals = []
    nEnrichmentScores = []
    nEnrichmentNulls = []

    for i in range(len(enrichment_scores)):
        es = enrichment_scores[i]
        enrNull = enrichment_nulls[i]
        #print es, enrNull

        enrichmentPVals.append(gsea_pval(es, enrNull))

        #normalize the ES(S,pi) and the observed ES(S), separetely rescaling
        #the positive and negative scores by divident by the mean of the 
        #ES(S,pi)

        #print es, enrNull

        def normalize(s):
            try:
                if s == 0:
                    return 0.0
                if s >= 0:
                    meanPos = np.mean([a for a in enrNull if a >= 0])
                    #print s, meanPos
                    return s/meanPos
                else:
                    meanNeg = np.mean([a for a in enrNull if a < 0])
                    #print s, meanNeg
                    return -s/meanNeg
            except:
                return 0.0 #return if according mean value is uncalculable


        nes = normalize(es)
        nEnrichmentScores.append(nes)
        
        nenrNull = [ normalize(s) for s in enrNull ]
        nEnrichmentNulls.append(nenrNull)
 

    print("Enrichment Score computing finished............. ", time.time() - tb1)

    #FDR computation
    #create a histogram of all NES(S,pi) over all S and pi
    vals = reduce(lambda x,y: x+y, nEnrichmentNulls, [])

    """
    def shorten(l, p=10000):
        
        #Take each len(l)/p element, if len(l)/p >= 2.
        
        e = len(l)/p
        if e <= 1:
            return l
        else:
            return [ l[i] for i in range(0, len(l), e) ]

    #vals = shorten(vals) -> this can speed up second part. is it relevant TODO?
    """
   
   
   
   
   
   
    """
    Use this null distribution to compute an FDR q value, for a given NES(S) =
    NES* >= 0. The FDR is the ratio of the percantage of all (S,pi) with
    NES(S,pi) >= 0, whose NES(S,pi) >= NES*, divided by the percentage of
    observed S wih NES(S) >= 0, whose NES(S) >= NES*, and similarly if NES(S)
    = NES* <= 0.
    """

    nvals = np.array(sorted(vals))
    nnes = np.array(sorted(nEnrichmentScores))

    #print("LEN VALS", len(vals), len(nEnrichmentScores))

    fdrs = []


    for i in range(len(enrichment_scores)):

        nes = nEnrichmentScores[i]


        #this could be speed up twice with the same accuracy! 
        if nes >= 0:
            allPos = int(len(vals) - np.searchsorted(nvals, 0, side="left"))
            allHigherAndPos = int(len(vals) - np.searchsorted(nvals, nes, side="left"))
            nesPos = len(nnes) - int(np.searchsorted(nnes, 0, side="left"))
            nesHigherAndPos = len(nnes) - int(np.searchsorted(nnes, nes, side="left"))
        else:
            allPos = int(np.searchsorted(nvals, 0, side="left"))
            allHigherAndPos = int(np.searchsorted(nvals, nes, side="right"))
            nesPos = int(np.searchsorted(nnes, 0, side="left"))
            nesHigherAndPos = int(np.searchsorted(nnes, nes, side="right"))
           


        try:
            top = allHigherAndPos/float(allPos) #p value
            down = nesHigherAndPos/float(nesPos)

            fdrs.append(top/down)
        except:
            fdrs.append(1000000000.0)
    
    print("Statistial testing finished.", time.time() - tb1)

    return zip(enrichment_scores, nEnrichmentScores, enrichmentPVals, fdrs)

