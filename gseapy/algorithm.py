# -*- coding: utf-8 -*-

from __future__ import  division

import sys, logging
import numpy as np
from functools import reduce
from multiprocessing import Pool
from math import ceil

def enrichment_score(gene_list, correl_vector, gene_set, weighted_score_type=1, 
                     nperm=1000, rs=np.random.RandomState(), single=False, scale=False):
    """This is the most important function of GSEApy. It has the same algorithm with GSEA and ssGSEA.

    :param gene_list:       The ordered gene list gene_name_list, rank_metric.index.values
    :param gene_set:        gene_sets in gmt file, please used gsea_gmt_parser to get gene_set.
    :param weighted_score_type:  It's same with gsea's weighted_score method. weighting by the correlation
                            is a very reasonable choice that allows significant gene sets with less than perfect coherence.
                            options: 0(classic),1,1.5,2. default:1. if one is interested in penalizing sets for lack of
                            coherence or to discover sets with any type of nonrandom distribution of tags, a value p < 1
                            might be appropriate. On the other hand, if one uses sets with large number of genes and only
                            a small subset of those is expected to be coherent, then one could consider using p > 1.
                            Our recommendation is to use p = 1 and use other settings only if you are very experienced
                            with the method and its behavior.

    :param correl_vector:   A vector with the correlations (e.g. signal to noise scores) corresponding to the genes in
                            the gene list. Or rankings, rank_metric.values
    :param nperm:           Only used this parameter when computing esnull for statistical testing. set the esnull value
                            equal to the permutation number.
    :param rs:              Random state for initialize gene list shuffling. Default: np.random.RandomState(seed=None)

    :return:

     ES: Enrichment score (real number between -1 and +1)

     ESNULL: Enrichment score calculated from random permutation.

     Hits_Indices: index of a gene in gene_list, if gene included in gene_set.

     RES: Numerical vector containing the running enrichment score for all locations in the gene list .

    """

    N = len(gene_list)
    # Test whether each element of a 1-D array is also present in a second array
    # It's more intuitived here than orginal enrichment_score source code.
    # use .astype to covert bool to intergers
    tag_indicator = np.in1d(gene_list, gene_set, assume_unique=True).astype(int)  # notice that the sign is 0 (no tag) or 1 (tag)

    if weighted_score_type == 0 :
        correl_vector = np.repeat(1, N)
    else:
        correl_vector = np.abs(correl_vector)**weighted_score_type

    # get indices of tag_indicator
    hit_ind = np.flatnonzero(tag_indicator).tolist()
    # if used for compute esnull, set esnull equal to permutation number, e.g. 1000
    # else just compute enrichment scores
    # set axis to 1, because we have 2 dimentional array
    axis = 1
    tag_indicator = np.tile(tag_indicator, (nperm+1,1))
    correl_vector = np.tile(correl_vector,(nperm+1,1))
    # gene list permutation
    for i in range(nperm): rs.shuffle(tag_indicator[i])
    # np.apply_along_axis(rs.shuffle, 1, tag_indicator)

    Nhint = tag_indicator.sum(axis=axis, keepdims=True)
    sum_correl_tag = np.sum(correl_vector*tag_indicator, axis=axis, keepdims=True)
    # compute ES score, the code below is identical to gsea enrichment_score method.
    no_tag_indicator = 1 - tag_indicator
    Nmiss =  N - Nhint
    norm_tag =  1.0/sum_correl_tag
    norm_no_tag = 1.0/Nmiss

    RES = np.cumsum(tag_indicator * correl_vector * norm_tag - no_tag_indicator * norm_no_tag, axis=axis)

    if scale: RES = RES / N
    if single:
        es_vec = RES.sum(axis=axis)
    else:
        max_ES, min_ES =  RES.max(axis=axis), RES.min(axis=axis)
        es_vec = np.where(np.abs(max_ES) > np.abs(min_ES), max_ES, min_ES)
    # extract values
    es, esnull, RES = es_vec[-1], es_vec[:-1], RES[-1,:]

    return es, esnull, hit_ind, RES



def enrichment_score_tensor(gene_mat, cor_mat, gene_sets, weighted_score_type, nperm=1000,
                            rs=np.random.RandomState(), single=False, scale=False):
    """Next generation algorithm of GSEA and ssGSEA.

        :param gene_mat:        the ordered gene list(vector) with or without gene indices matrix.
        :param cor_mat:         correlation vector or matrix  (e.g. signal to noise scores)
                                corresponding to the genes in the gene list or matrix.
        :param dict gene_sets:  gmt file dict.
        :param float weighted_score_type:     weighting by the correlation.
                                options: 0(classic), 1, 1.5, 2. default:1 for GSEA and 0.25 for ssGSEA.
        :param int nperm:       permutation times.
        :param bool scale:      If True, normalize the scores by number of genes_mat.
        :param bool single:     If True, use ssGSEA algorithm, otherwise use GSEA.
        :param rs:              Random state for initialize gene list shuffling.
                                Default: np.random.RandomState(seed=None)
        :return: a tuple contains::

                 | ES: Enrichment score (real number between -1 and +1), for ssGSEA, set scale eq to True.
                 | ESNULL: Enrichment score calculated from random permutation.
                 | Hits_Indices: Indices of genes if genes are included in gene_set.
                 | RES: The running enrichment score for all locations in the gene list.

    """
    # gene_mat -> 1d: prerank, ssSSEA or 2d: GSEA
    keys = sorted(gene_sets.keys())

    if weighted_score_type == 0:
        # don't bother doing calcuation, just set to 1
        cor_mat = np.ones(cor_mat.shape)
    elif weighted_score_type > 0:
        pass
    else:
        logging.error("Using negative values of weighted_score_type, not allowed")
        sys.exit(0)

    cor_mat = np.abs(cor_mat)
    if cor_mat.ndim ==1:
        # ssGSEA or Prerank
        # genestes->M, genes->N, perm-> axis=2
        N, M = len(gene_mat), len(keys)
        # generate gene hits matrix
        # for 1d ndarray of gene_mat, set assume_unique=True,
        # means the input arrays are both assumed to be unique,
        # which can speed up the calculation.
        tag_indicator = np.vstack([np.in1d(gene_mat, gene_sets[key], assume_unique=True) for key in keys])
        tag_indicator = tag_indicator.astype(int)
        # index of hits
        hit_ind = [ np.flatnonzero(tag).tolist() for tag in tag_indicator ]
        # generate permutated hits matrix
        perm_tag_tensor = np.repeat(tag_indicator, nperm+1).reshape((M,N,nperm+1))
        # shuffle matrix, last matrix is not shuffled when nperm > 0
        if nperm: np.apply_along_axis(lambda x: np.apply_along_axis(rs.shuffle,0,x),1, perm_tag_tensor[:,:,:-1])
        # missing hits
        no_tag_tensor = 1 - perm_tag_tensor
        # calculate numerator, denominator of each gene hits
        rank_alpha = (perm_tag_tensor*cor_mat[np.newaxis,:,np.newaxis])** weighted_score_type

    elif cor_mat.ndim == 2:
        # GSEA
        # 2d ndarray, gene_mat and cor_mat are shuffled already
        # reshape matrix
        cor_mat = cor_mat.T
        # gene_mat is a tuple contains (gene_name, permuate_gene_name_indices)
        genes, genes_ind = gene_mat
        # genestes->M, genes->N, perm-> axis=2
        # don't use assume_unique=True in 2d array when use np.isin().
        # elements in gene_mat are not unique, or will cause unwanted results
        tag_indicator = np.vstack([np.in1d(genes, gene_sets[key], assume_unique=True) for key in keys])
        tag_indicator = tag_indicator.astype(int)
        perm_tag_tensor = np.stack([tag.take(genes_ind).T for tag in tag_indicator], axis=0)
        #index of hits
        hit_ind = [ np.flatnonzero(tag).tolist() for tag in perm_tag_tensor[:,:,-1] ]
        # nohits
        no_tag_tensor = 1 - perm_tag_tensor
        # calculate numerator, denominator of each gene hits
        rank_alpha = (perm_tag_tensor*cor_mat[np.newaxis,:,:])** weighted_score_type
    else:
        logging.error("Program die because of unsupported input")
        sys.exit(0)

    # Nhint = tag_indicator.sum(1)
    # Nmiss =  N - Nhint
    axis=1
    P_GW_denominator = np.sum(rank_alpha, axis=axis, keepdims=True)
    P_NG_denominator = np.sum(no_tag_tensor, axis=axis, keepdims=True)
    REStensor = np.cumsum(rank_alpha / P_GW_denominator - no_tag_tensor / P_NG_denominator, axis=axis)
    # ssGSEA: scale es by gene numbers ?
    # https://gist.github.com/gaoce/39e0907146c752c127728ad74e123b33
    if scale: REStensor = REStensor / len(gene_mat)
    if single:
        #ssGSEA
        esmatrix = REStensor.sum(axis=axis)
    else:
        #GSEA
        esmax, esmin = REStensor.max(axis=axis), REStensor.min(axis=axis)
        esmatrix = np.where(np.abs(esmax)>np.abs(esmin), esmax, esmin)

    es, esnull, RES = esmatrix[:,-1], esmatrix[:,:-1], REStensor[:,:,-1]

    return es, esnull, hit_ind, RES


def ranking_metric_tensor(exprs, method, permutation_num, pos, neg, classes,
                          ascending, rs=np.random.RandomState()):
    """Build shuffled ranking matrix when permutation_type eq to phenotype.

       :param exprs:   gene_expression DataFrame, gene_name indexed.
       :param str method:  calculate correlation or ranking. methods including:
                           1. 'signal_to_noise'.
                           2. 't_test'.
                           3. 'ratio_of_classes' (also referred to as fold change).
                           4. 'diff_of_classes'.
                           5. 'log2_ratio_of_classes'.
       :param int permuation_num: how many times of classes is being shuffled
       :param str pos: one of labels of phenotype's names.
       :param str neg: one of labels of phenotype's names.
       :param list classes:  a list of phenotype labels, to specify which column of
                             dataframe belongs to what class of phenotype.
       :param bool ascending:  bool. Sort ascending vs. descending.

       :return:
                returns two 2d ndarray with shape (nperm, gene_num).

                | cor_mat_indices: the indices of sorted and permutate (exclude last row) ranking matrix.
                | cor_mat: sorted and permutate (exclude last row) ranking matrix.

    """
    # S: samples, G: gene number
    G, S = exprs.shape
    # genes = exprs.index.values
    expr_mat = exprs.values.T
    perm_cor_tensor = np.tile(expr_mat, (permutation_num+1,1,1))
    # random shuffle on the first dim, last matrix is not shuffled
    for arr in perm_cor_tensor[:-1]: rs.shuffle(arr)
    classes = np.array(classes)
    pos = classes == pos
    neg = classes == neg
    pos_cor_mean = perm_cor_tensor[:,pos,:].mean(axis=1)
    neg_cor_mean = perm_cor_tensor[:,neg,:].mean(axis=1)
    pos_cor_std = perm_cor_tensor[:,pos,:].std(axis=1, ddof=1)
    neg_cor_std = perm_cor_tensor[:,neg,:].std(axis=1, ddof=1)

    if method == 'signal_to_noise':
        cor_mat = (pos_cor_mean - neg_cor_mean)/(pos_cor_std + neg_cor_std)
    elif method == 't_test':
        denom = 1.0/G
        cor_mat = (pos_cor_mean - neg_cor_mean)/ np.sqrt(denom*pos_cor_std**2 + denom*neg_cor_std**2)
    elif method == 'ratio_of_classes':
        cor_mat = pos_cor_mean / neg_cor_mean
    elif method == 'diff_of_classes':
        cor_mat  = pos_cor_mean - neg_cor_mean
    elif method == 'log2_ratio_of_classes':
        cor_mat  =  np.log2(pos_cor_mean / neg_cor_mean)
    else:
        logging.error("Please provide correct method name!!!")
        sys.exit(0)
    # return matix[nperm+1, perm_cors]
    cor_mat_ind = cor_mat.argsort()
    # ndarray: sort in place
    cor_mat.sort()
    # genes_mat = genes.take(cor_mat_ind)
    if ascending: return cor_mat_ind, cor_mat
    # descending order of ranking and genes
    # return genes_mat[:,::-1], cor_mat[:,::-1]
    return cor_mat_ind[:, ::-1], cor_mat[:, ::-1]

def ranking_metric(df, method, pos, neg, classes, ascending):
    """The main function to rank an expression table.

       :param df:      gene_expression DataFrame.
       :param method:  The method used to calculate a correlation or ranking. Default: 'log2_ratio_of_classes'.
                       Others methods are:

                       1. 'signal_to_noise'

                          You must have at least three samples for each phenotype to use this metric.
                          The larger the signal-to-noise ratio, the larger the differences of the means (scaled by the standard deviations);
                          that is, the more distinct the gene expression is in each phenotype and the more the gene acts as a “class marker.”

                       2. 't_test'

                          Uses the difference of means scaled by the standard deviation and number of samples.
                          Note: You must have at least three samples for each phenotype to use this metric.
                          The larger the tTest ratio, the more distinct the gene expression is in each phenotype
                          and the more the gene acts as a “class marker.”

                       3. 'ratio_of_classes' (also referred to as fold change).

                          Uses the ratio of class means to calculate fold change for natural scale data.

                       4. 'diff_of_classes'

                          Uses the difference of class means to calculate fold change for natural scale data

                       5. 'log2_ratio_of_classes'

                          Uses the log2 ratio of class means to calculate fold change for natural scale data.
                          This is the recommended statistic for calculating fold change for log scale data.


       :param str pos: one of labels of phenotype's names.
       :param str neg: one of labels of phenotype's names.
       :param list classes:  a list of phenotype labels, to specify which column of dataframe belongs to what category of phenotype.
       :param bool ascending:  bool or list of bool. Sort ascending vs. descending.
       :return:

            returns a pd.Series of correlation to class of each variable. Gene_name is index, and value is rankings.

            visit here for more docs: http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html
    """

    # exclude any zero stds.
    df_mean = df.groupby(by=classes, axis=1).mean()
    df_std =  df.groupby(by=classes, axis=1).std()


    if method == 'signal_to_noise':
        ser = (df_mean[pos] - df_mean[neg])/(df_std[pos] + df_std[neg])
    elif method == 't_test':
        ser = (df_mean[pos] - df_mean[neg])/ np.sqrt(df_std[pos]**2/len(df_std)+df_std[neg]**2/len(df_std) )
    elif method == 'ratio_of_classes':
        ser = df_mean[pos] / df_mean[neg]
    elif method == 'diff_of_classes':
        ser  = df_mean[pos] - df_mean[neg]
    elif method == 'log2_ratio_of_classes':
        ser  =  np.log2(df_mean[pos] / df_mean[neg])
    else:
        logging.error("Please provide correct method name!!!")
        sys.exit(0)
    ser = ser.sort_values(ascending=ascending)

    return ser


def gsea_compute_tensor(data, gmt, n, weighted_score_type, permutation_type,
                 method, pheno_pos, pheno_neg, classes, ascending,
                 processes=1, seed=None, single=False, scale=False):
    """compute enrichment scores and enrichment nulls.

        :param data: preprocessed expression dataframe or a pre-ranked file if prerank=True.
        :param dict gmt: all gene sets in .gmt file. need to call load_gmt() to get results.
        :param int n: permutation number. default: 1000.
        :param str method: ranking_metric method. see above.
        :param str pheno_pos: one of labels of phenotype's names.
        :param str pheno_neg: one of labels of phenotype's names.
        :param list classes: a list of phenotype labels, to specify which column of dataframe belongs to what category of phenotype.
        :param float weighted_score_type: default:1
        :param bool ascending: sorting order of rankings. Default: False.
        :param seed: random seed. Default: np.random.RandomState()
        :param bool scale: if true, scale es by gene number.

        :return: a tuple contains::

                | zipped results of es, nes, pval, fdr.
                | nested list of hit indices of input gene_list.
                | nested list of ranked enrichment score of each input gene_sets.
                | list of enriched terms

    """
    w = weighted_score_type
    subsets = sorted(gmt.keys())
    rs = np.random.RandomState(seed)
    genes_mat, cor_mat = data.index.values, data.values
    base = 5 if data.shape[0] >= 5000 else 10
    block = ceil(len(subsets) / base)

    if permutation_type == "phenotype":
        # shuffling classes and generate random correlation rankings
        logging.debug("Start to permutate classes..............................")
        genes_ind = []
        cor_mat = []
        temp_rnk = []
        pool_rnk = Pool(processes=processes)
        # split large array into smaller blocks to avoid memory overflow
        i=1
        while i <= block:
            rs = np.random.RandomState(seed)
            temp_rnk.append(pool_rnk.apply_async(ranking_metric_tensor,
                                         args=(data, method, base, pheno_pos, pheno_neg, classes,
                                                ascending, rs)))
            i +=1
        pool_rnk.close()
        pool_rnk.join()

        for k, temp in enumerate(temp_rnk):
            gi, cor = temp.get()
            if k+1 == block:
               genes_ind.append(gi)
               cor_mat.append(cor)
            else:
                genes_ind.append(gi[:-1])
                cor_mat.append(cor[:-1])

        genes_ind, cor_mat = np.vstack(genes_ind), np.vstack(cor_mat)
        # convert to tuple
        genes_mat = (data.index.values, genes_ind)

    logging.debug("Start to compute es and esnulls........................")
    # Prerank, ssGSEA, GSEA
    es = []
    RES = []
    hit_ind = []
    esnull = []
    temp_esnu = []
    pool_esnu = Pool(processes=processes)
    # split large array into smaller blocks to avoid memory overflow
    i, m = 1, 0
    while i <= block:
        # you have to reseed, or all your processes are sharing the same seed value
        rs = np.random.RandomState(seed)
        gmtrim = {k: gmt.get(k) for k in subsets[m:base * i]}
        temp_esnu.append(pool_esnu.apply_async(enrichment_score_tensor,
                                               args=(genes_mat, cor_mat,
                                                     gmtrim, w, n, rs,
                                                     single, scale)))
        m = base * i
        i += 1
    pool_esnu.close()
    pool_esnu.join()
    # esn is a list, don't need to use append method.
    for si, temp in enumerate(temp_esnu):
        e, enu, hit, rune = temp.get()
        esnull.append(enu)
        es.append(e)
        RES.append(rune)
        hit_ind += hit
    # concate results
    es, esnull, RES = np.hstack(es), np.vstack(esnull), np.vstack(RES)

    return gsea_significance(es, esnull), hit_ind, RES, subsets



def gsea_compute(data, gmt, n, weighted_score_type, permutation_type,
                 method, pheno_pos, pheno_neg, classes, ascending,
                 processes=1, seed=None, single=False, scale=False):
    """compute enrichment scores and enrichment nulls.

        :param data: preprocessed expression dataframe or a pre-ranked file if prerank=True.
        :param dict gmt: all gene sets in .gmt file. need to call load_gmt() to get results.
        :param int n: permutation number. default: 1000.
        :param str method: ranking_metric method. see above.
        :param str pheno_pos: one of labels of phenotype's names.
        :param str pheno_neg: one of labels of phenotype's names.
        :param list classes: a list of phenotype labels, to specify which column of dataframe belongs to what category of phenotype.
        :param float weighted_score_type: default:1
        :param bool ascending: sorting order of rankings. Default: False.
        :param seed: random seed. Default: np.random.RandomState()
        :param bool scale: if true, scale es by gene number.

        :return: a tuple contains::

                | zipped results of es, nes, pval, fdr.
                | nested list of hit indices of input gene_list.
                | nested list of ranked enrichment score of each input gene_sets.
                | list of enriched terms

    """
    
    w = weighted_score_type
    subsets = sorted(gmt.keys())
    es = []
    RES=[]
    hit_ind=[]
    esnull = [ [] for a in range(len(subsets)) ]

    logging.debug("Start to compute enrichment scores......................")

    if permutation_type == "phenotype":
        logging.debug("Start to permutate classes..............................")
        # shuffling classes and generate random correlation rankings
        rs = np.random.RandomState(seed)
        genes_mat, cor_mat = ranking_metric_tensor(exprs=data, method=method,
                                                   permutation_num=n,
                                                   pos=pheno_pos, neg=pheno_neg,
                                                   classes=classes,
                                                   ascending=ascending, rs=rs)

        # compute es, esnulls. hits, RES
        logging.debug("Start to compute enrichment nulls.......................")
        es, esnull, hit_ind, RES = enrichment_score_tensor(gene_mat=genes_mat,
                                                           cor_mat=cor_mat,
                                                           gene_sets=gmt,
                                                           weighted_score_type=w,
                                                           nperm=n, rs=rs,
                                                           single=False, scale=False,)

    else:
        # Prerank, ssGSEA, GSEA with gene_set permutation
        gl, cor_vec = data.index.values, data.values
        logging.debug("Start to compute es and esnulls........................")
        # es, esnull, hit_ind, RES = enrichment_score_tensor(gene_mat=gl,
        #                                                    cor_mat=cor_mat,
        #                                                    gene_sets=gmt,
        #                                                    weighted_score_type=w,
        #                                                    nperm=n, rs=rs
        #                                                    single=single, scale=scale)

        # split large array into smaller blocks to avoid memory overflow
        temp_esnu=[]
        pool_esnu = Pool(processes=processes)
        for subset in subsets:
            # you have to reseed, or all your processes are sharing the same seed value
            rs = np.random.RandomState(seed)
            temp_esnu.append(pool_esnu.apply_async(enrichment_score,
                                                   args=(gl, cor_vec, gmt.get(subset), w,
                                                         n, rs, single, scale)))

        pool_esnu.close()
        pool_esnu.join()
        # esn is a list, don't need to use append method.
        for si, temp in enumerate(temp_esnu):
            e, enu, hit, rune = temp.get()
            esnull[si] = enu
            es.append(e)
            RES.append(rune)
            hit_ind.append(hit)

    return gsea_significance(es, esnull), hit_ind, RES, subsets

def gsea_pval(es, esnull):
    """Compute nominal p-value.

    From article (PNAS):
    estimate nominal p-value for S from esnull by using the positive
    or negative portion of the distribution corresponding to the sign
    of the observed ES(S).
    """

    # to speed up, using numpy function to compute pval in parallel.
    es = np.array(es)
    esnull = np.array(esnull)
    # try:
    condlist = [ es < 0, es >=0]
    choicelist = [np.sum(esnull < es.reshape(len(es),1), axis=1)/ np.sum(esnull < 0, axis=1),
                  np.sum(esnull >= es.reshape(len(es),1), axis=1)/ np.sum(esnull >= 0, axis=1)]
    pval = np.select(condlist, choicelist)

    return pval
    # except:
    #    return np.repeat(1.0 ,len(es))

def normalize(es, esnull):
    """normalize the ES(S,pi) and the observed ES(S), separately rescaling
       the positive and negative scores by dividing the mean of the ES(S,pi).
    """

    try:
        if es == 0:
            return 0.0
        if es >= 0:
            meanPos = np.mean([a for a in esnull if a >= 0])
            #print es, meanPos
            return es/meanPos
        else:
            meanNeg = np.mean([a for a in esnull if a < 0])
            #print es, meanNeg
            return -es/meanNeg
    except:
        return 0.0 #return if according mean value is uncalculable


def gsea_significance(enrichment_scores, enrichment_nulls):
    """Compute nominal pvals, normalized ES, and FDR q value.

        For a given NES(S) = NES* >= 0. The FDR is the ratio of the percentage of all (S,pi) with
        NES(S,pi) >= 0, whose NES(S,pi) >= NES*, divided by the percentage of
        observed S wih NES(S) >= 0, whose NES(S) >= NES*, and similarly if NES(S) = NES* <= 0.
    """
    # For a zero by zero division (undetermined, results in a NaN),
    # np.seterr(divide='ignore', invalid='ignore')
    import warnings
    warnings.simplefilter("ignore")
    
    logging.debug("Start to compute pvals..................................")
    # compute pvals.
    enrichmentPVals = gsea_pval(enrichment_scores, enrichment_nulls).tolist()

    # new normalize enrichment score calculating method. this could speed up significantly.
    esnull_meanPos = []
    esnull_meanNeg = []

    es = np.array(enrichment_scores)
    esnull = np.array(enrichment_nulls)

    for i in range(len(enrichment_scores)):
        enrNull = esnull[i]
        meanPos = enrNull[enrNull >= 0].mean()
        esnull_meanPos.append(meanPos)
        meanNeg = enrNull[enrNull < 0 ].mean()
        esnull_meanNeg.append(meanNeg)

    pos = np.array(esnull_meanPos).reshape(len(es), 1)
    neg = np.array(esnull_meanNeg).reshape(len(es), 1)

    # compute normalized enrichment score and normalized esnull
    logging.debug("Compute normalized enrichment score and normalized esnull")

    try:
        condlist1 = [ es >= 0, es < 0]
        choicelist1 = [ es/esnull_meanPos, -es/esnull_meanNeg ]
        nEnrichmentScores = np.select(condlist1, choicelist1).tolist()

        condlist2 = [ esnull >= 0, esnull < 0]
        choicelist2 = [ esnull/pos, -esnull/neg ]
        nEnrichmentNulls = np.select(condlist2, choicelist2)

    except:  #return if according nes, nesnull is uncalculable
        nEnrichmentScores = np.repeat(0.0, es.size).tolist()
        nEnrichmentNulls = np.repeat(0.0 , es.size).reshape(esnull.shape)


    logging.debug("start to compute fdrs..................................")

    # FDR null distribution histogram
    # create a histogram of all NES(S,pi) over all S and pi
    # Use this null distribution to compute an FDR q value,
    # vals = reduce(lambda x,y: x+y, nEnrichmentNulls, [])
    # nvals = np.array(sorted(vals))
    # or
    nvals = np.sort(nEnrichmentNulls.flatten())
    nnes = np.array(sorted(nEnrichmentScores))
    fdrs = []
    # FDR computation
    for i in range(len(enrichment_scores)):
        nes = nEnrichmentScores[i]

        if nes >= 0:
            allPos = int(len(nvals) - np.searchsorted(nvals, 0, side="left"))
            allHigherAndPos = int(len(nvals) - np.searchsorted(nvals, nes, side="left"))
            nesPos = len(nnes) - int(np.searchsorted(nnes, 0, side="left"))
            nesHigherAndPos = len(nnes) - int(np.searchsorted(nnes, nes, side="left"))
        else:
            allPos = int(np.searchsorted(nvals, 0, side="left"))
            allHigherAndPos = int(np.searchsorted(nvals, nes, side="right"))
            nesPos = int(np.searchsorted(nnes, 0, side="left"))
            nesHigherAndPos = int(np.searchsorted(nnes, nes, side="right"))
        try:
            pi_norm = allHigherAndPos/float(allPos) #p value
            pi_obs = nesHigherAndPos/float(nesPos)

            fdr = pi_norm/pi_obs if pi_norm/pi_obs < 1.0  else 1.0
            fdrs.append(fdr)
        except:
            fdrs.append(1000000000.0)

    logging.debug("Statistical testing finished.............................")

    return zip(enrichment_scores, nEnrichmentScores, enrichmentPVals, fdrs)
