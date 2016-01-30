# -*- coding: utf-8 -*-


from __future__ import absolute_import, print_function,division
from functools import reduce

import time
import numpy as np
import pandas as pd
import sys
import random


def preprocess(data):
    """pre-processed the data frame.
       new filtering methods will be implement here.
    """
    df = pd.read_table(data,index_col=0)
    
    assert len(df) > 1
    
    ids = df.index.drop_duplicates(keep='first') #drop duplicate gene_names.    
    df = df.loc[ids]
    df.dropna(how='all',inplace=True)
    df2 = df.select_dtypes(include=['float64'])  + 0.0001 #select numbers in DataFrame      
    
    return df2


def enrichment_score(gene_list, gene_set, weighted_score_type = 1, correl_vector = None):
    '''
    This is the most important function of GSEAPY. It has the same algorithm with GSEA.
    
    :param gene_list:       The ordered gene list 
                            gene_name_list, rank_metric['gene_name']

    :param gene_set:        gene_sets in gmt file, please used gsea_gmt_parser to get gene_set.
                             
    :param weighted_score_type: A vector with the coorelations (e.g. signal to noise scores) 
                                corresponding to the genes in the gene list.
                                It's indentical to gsea's weighted_sore_method. options: 0(classic),1,1.5,2. default:1.
    :param correl_vector:   A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in 
                            the gene list. Or rankings, rank_metric['rank'].values
    
    
    :return: 
    |ES: Enrichment score (real number between -1 and +1) 
    |hit_index: index of a gene in gene_list, if gene included in gene_set.
    |RES: Numerical vector containing the running enrichment score for all locations in the gene list .
             
    '''
    
    #Test whether each element of a 1-D array is also present in a second array
    #It's more intuitived here than orginal enrichment_score source code.
    #use .astype to covert bool to intergers  
    tag_indicator = np.in1d(gene_list,gene_set,assume_unique=True).astype(int)  # notice that the sign is 0 (no tag) or 1 (tag)
    no_tag_indicator = 1 - tag_indicator
    hit_ind = np.flatnonzero(tag_indicator).tolist()
     
    
    #compute ES score, the code below is identical to gsea enrichment_score method.
    N = len(gene_list) 
    Nhint = np.sum(tag_indicator)
    Nmiss =  N - Nhint 
    if (weighted_score_type == 0 ): 
        correl_vector = np.repeat(1, N)
    else:
        correl_vector = np.abs(correl_vector**weighted_score_type)
    sum_correl_tag = np.sum(correl_vector[tag_indicator.astype(bool)])
    norm_tag =  1.0/sum_correl_tag
    norm_no_tag = 1.0/Nmiss
    RES = np.cumsum(tag_indicator * correl_vector * norm_tag - no_tag_indicator * norm_no_tag)      
    max_ES = max(RES)
    min_ES = min(RES)
    
    #print("The length of ES is", len(RES))
    
    return max_ES if abs(max_ES) > abs(min_ES) else min_ES, hit_ind, RES.tolist()
 

def shuffle_list(gene_list, rand=random.Random(0)):
    """
    Returns a copy of a shuffled input gene_list.
    
    :gene_list: rank_metric['gene_name'].values
    """
    
    l2 = gene_list.copy()
    rand.shuffle(l2)

    return l2        
         
    
    
    
def ranking_metric(df, method,phenoPos,phenoNeg,classes,ascending):
    """
     The main function to rank a expression table.
    
    :param df: gene_expression DataFrame.    
    :param method:
     1. signal_to_noise 
         You must have at least three samples for each phenotype to use this metric.
         The larger the signal-to-noise ratio, the larger the differences of the means (scaled by the standard deviations);
         that is, the more distinct the gene expression is in each phenotype and the more the gene acts as a “class marker.” 
     2. t_test
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
    :param phenoPos: one of lables of phenotype's names 
    :param phenoNeg: one of lable of phenotype's names     
    :param classes: a list of phenotype labels, to specify which column of dataframe belongs to what catogry of phenotype.
    :param ascending:  bool or list of bool. Sort ascending vs. descending.
    :param cls_path: only use when classes is not None.
 
   :return: returns correlation to class of each variable.
             same format with .rnk file. gene_name in first coloum, correlation
             in second column.
    """ 
    
    
    A = phenoPos
    B = phenoNeg
    df2 = df.T   
    df2['class'] = classes
    df_mean= df2.groupby('class').mean().T
    df_std = df2.groupby('class').std().T
    
    if method == 'signal_to_noise':
        sr = (df_mean[A] - df_mean[B])/(df_std[A] + df_std[B])
    elif method == 't_test':
        sr = (df_mean[A] - df_mean[B])/ np.sqrt(df_std[A]**2/len(df_std)+df_std[B]**2/len(df_std) )
    elif method == 'ratio_of_classes':
        sr = df_mean[A] / df_mean[B]
    elif method == 'diff_of_classes':
        sr  = df_mean[A] - df_mean[B]
    elif method == 'log2_ratio_of_classes':
        sr  =  np.log2(df_mean[A] / df_mean[B])
    else:
        print("Please provide correct method name!!!")        
        sys.exit()
    sr.sort_values(ascending=ascending,inplace=True)
    df3 = sr.to_frame().reset_index()
    df3.columns = ['gene_name','rank']
    df3['rank2'] = df3['rank']

    return df3
        
    
def gsea_compute(data, gmt, n, weighted_score_type,permutation_type,method,phenoPos,phenoNeg,classes,ascending):
    """
    compute enrichment scores and enrichment nulls. 
    
    :param data: prepreocessed expression dataframe.
    :param gmt: all gene sets in .gmt file. need to call gsea_gmt_parser() to get results. 
    :param n: permutation number. default: 1000
    :param method: ranking_metric method. see above.
    :param phenoPos: one of lables of phenotype's names 
    :param phenoNeg: one of lable of phenotype's names     
    :param classes: a list of phenotype labels, to specify which column of dataframe belongs to what catogry of phenotype.
    :param weighted_score_type: default:1
    :param ascending: sorting order of rankings. Default: False.
    """
    enrichment_scores = []

    w = weighted_score_type
    subsets = sorted(gmt.keys())
    
    dat = data.copy()
    r2 = ranking_metric(df=dat,method = method,phenoPos=phenoPos,phenoNeg=phenoNeg,classes=classes, ascending= ascending)
    ranking=r2['rank'].values
    gene_list=r2['gene_name']
        

    print("Start to compute enrichment socres............",time.ctime())

    rank_ES = []
    hit_ind = []
    
    for subset in subsets:
        es,ind,RES = enrichment_score(gene_list = gene_list, gene_set=gmt.get(subset), 
                              weighted_score_type = w, correl_vector = ranking)
        enrichment_scores.append(es)
        rank_ES.append(RES)
        hit_ind.append(ind)
    
    
    print("Star to compute esnulls.......................",time.ctime())

    enrichment_nulls = [ [] for a in range(len(subsets)) ]
    

    for i in range(n):
        if permutation_type == "phenotype":
            dat2 = dat.T   
            dat2.apply(np.random.shuffle,axis=0) #permutation classes
            r2 = ranking_metric(df=dat2.T,method = method,phenoPos=phenoPos,phenoNeg=phenoNeg,classes=classes,ascending= ascending )
            ranking2=r2['rank']
            gene_list2=r2['gene_name'].values
        else:
            
            #gene_list.apply(np.random.shuffle,axis=0) #permutation genes
            #r2 = ranking_metric(df=dat,method = method, classes=classes,ascending= ascending)
            gene_list2 = shuffle_list(gene_list,random.Random(2000+i))
            ranking2=ranking
            #gene_list2=gene_list
            
        for si,subset in enumerate(subsets):
            esn = enrichment_score(gene_list = gene_list2, gene_set= gmt.get(subset), 
                              weighted_score_type = w, correl_vector = ranking2)[0]                        
            enrichment_nulls[si].append(esn)

        

    return gsea_significance(enrichment_scores, enrichment_nulls),hit_ind,rank_ES, subsets


def gsea_pval(es, esnull):
    """
    From article (PNAS):
    estimate nominal p-value for S from esnull by using the positive
    or negative portion of the distribution corresponding to the sign 
    of the observed ES(S).
    """
    
    try:
        if es < 0:
            return float(len([ a for a in esnull if a <= es ]))/len([ a for a in esnull if a < 0])    
        else: 
            return float(len([ a for a in esnull if a >= es ]))/len([ a for a in esnull if a >= 0])
    except:
        return 1.0

def gsea_significance(enrichment_scores, enrichment_nulls):
    """
    Computing p-vals, normalized ES, FDR
    """
    
    print("Start to compute pvals........................",time.ctime())
    
    enrichmentPVals = []
    nEnrichmentScores = []
    nEnrichmentNulls = []

    for i in range(len(enrichment_scores)):
        es = enrichment_scores[i]
        enrNull = enrichment_nulls[i]
       

        enrichmentPVals.append(gsea_pval(es, enrNull))

        #normalize the ES(S,pi) and the observed ES(S), separetely rescaling
        #the positive and negative scores by divident by the mean of the 
        #ES(S,pi)

        def normalize(s):
            try:
                if s == 0:
                    return 0.0
                if s >= 0:
                    meanPos = np.mean([a for a in enrNull if a >= 0])                   
                    return s/meanPos
                else:
                    meanNeg = np.mean([a for a in enrNull if a < 0])                    
                    return -s/meanNeg
            except:
                return 0.0 #return if according mean value is uncalculable


        nes = normalize(es)
        nEnrichmentScores.append(nes)
        
        nenrNull = [ normalize(s) for s in enrNull ]
        nEnrichmentNulls.append(nenrNull)
 

    print("start to compute fdrs.........................", time.ctime())

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

    print("Statistial testing finished...................", time.ctime())

    return zip(enrichment_scores, nEnrichmentScores, enrichmentPVals, fdrs)

