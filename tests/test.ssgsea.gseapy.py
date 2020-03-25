
# coding: utf-8

# In[2]:
import numpy as np
import pandas as pd
import gseapy as gp
import logging, sys


# In[3]:
np.seterr(divide='ignore')
print("GSEApy version: %s"%gp.__version__)

# In[15]:

# identical function with gseapy.algorithm.enrichment_score_tensor
def enrichment_score_tensor(gene_mat, cor_mat, gene_sets, weighted_score_type, nperm=1000,
                            scale=False, single=False, rs=np.random.RandomState()):
    """Next generation algorithm of GSEA and ssGSEA.

    :param gene_mat:        the ordered gene list(vector) or gene matrix.
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
    :return:
             ES: Enrichment score (real number between -1 and +1), it's true for ssGSEA, only scaled

             ESNULL: Enrichment score calculated from random permutation

             Hits_Indices: Indices of genes if genes are included in gene_set.

             RES: Numerical vector containing the running enrichment score for
                  all locations in the gene list .

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
        #genestes->M, genes->N, perm-> axis=2
        N, M = len(gene_mat), len(keys)
        # generate gene hits matrix
        # for 1d ndarray of gene_mat, set assume_unique=True,
        # means the input arrays are both assumed to be unique,
        # which can speed up the calculation.
        tag_indicator = np.vstack([np.in1d(gene_mat, gene_sets[key], assume_unique=True) for key in keys])
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
        cor_mat, gene_mat = cor_mat.T, gene_mat.T
        # genestes->M, genes->N, perm-> axis=2
        # don't use assume_unique=True in 2d array when use np.isin().
        # elements in gene_mat are not unique, or will cause unwanted results
        perm_tag_tensor = np.stack([np.isin(gene_mat, gene_sets[key]) for key in keys], axis=0)
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
        esmatrix = np.sum(REStensor, axis=axis)
    else:
        #GSEA
        esmax, esmin = REStensor.max(axis=axis), REStensor.min(axis=axis)
        esmatrix = np.where(np.abs(esmax)>np.abs(esmin), esmax, esmin)

    es, esnull, RES = esmatrix[:,-1], esmatrix[:,:-1], REStensor[:,:,-1]

    return es, esnull, hit_ind, RES


# In[16]:
gex = pd.read_table("./data/testSet_rand1200.gct", comment='#', index_col=0)
# In[18]:
gmt = gp.parser.gsea_gmt_parser("./data/randomSets.gmt")
# In[24]:



df = pd.DataFrame(index=sorted(gmt.keys()))
rs = np.random.RandomState(0)
gexrnk = gex.rank(axis=0, method='average', na_option='bottom')
for name, ser in gexrnk.iteritems():
    ser = ser.sort_values(ascending=False)
    est = enrichment_score_tensor(gene_mat=ser.index.values, cor_mat=ser.values, gene_sets=gmt,
                                  weighted_score_type=0.25, nperm=0,
                                  scale=True, single=True, rs=rs)[0]
    df[name]=est



# In[25]:
print("\n\n\n")
print("Scaled Enrichment Scores (ES):")
print(df)

## no scale es
df2 = pd.DataFrame(index=sorted(gmt.keys()))
rs = np.random.RandomState(0)
gexrnk = gex.rank(axis=0, method='average', na_option='bottom')
for name, ser in gexrnk.iteritems():
    ser = ser.sort_values(ascending=False)
    est = enrichment_score_tensor(gene_mat=ser.index.values, cor_mat=ser.values, gene_sets=gmt,
                                  weighted_score_type=0.25, nperm=0,
                                  scale=True, single=True, rs=rs)[0]
    df2[name]=est

# In[27]:
## no scale es values and norm by all samples
nes = df2/(df2.values.max() -df2.values.min())
print("\n\n\n")
print("Normalized Enrichment Scores (NES):")
# for n, e  in zip(names, nes):
#     print(n, ": ", e)
print(df2)
