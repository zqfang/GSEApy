
# -*- coding: utf-8 -*-
from __future__ import  division

import sys, logging
import numpy as np
from scipy.stats import hypergeom


def calc_pvalues(query, gene_sets, background=20000, **kwargs):
    """ calculate pvalues for all categories in the graph

    :param set query: set of identifiers for which the p value is calculated
    :param dict gene_sets: gmt file dict after background was set
    :param set background: total number of genes in your annoated database.
    :returns: pvalues
              x: overlaped gene number
              n:legth of gene_set belongs to each terms
              hits: overlaped gene names.


    For 2*2 contingency table: 
    =============================================================================
                         |   in  query  |  not in query |    row total
    =>      in gene_set  |        a     |       b       |       a+b  
    =>  not in gene_set  |        c     |       d       |       c+d  
           column total                                 | a+b+c+d = anno database
    =============================================================================
    backgroud genes number = a + b + c + d.

    Then, in R
        x=a     the number of white balls drawn without replacement 
                from an urn which contains both black and white balls.
        m=a+b   the number of white balls in the urn    
        n=c+d   the number of black balls in the urn    
        k=a+c   the number of balls drawn from the urn  
   
    In Scipy:
    for args in scipy.hypergeom.sf(k, M, n, N, loc=0):
        M: the total number of objects, 
        n: the total number of Type I objects. 
        k: the random variate represents the number of Type I objects in N drawn 
           without replacement from the total population.
    
    Therefore, these two fuction are same when using parameters from 2*2 table:
    R:     >   phyper(x-1, m, n, k, lower.tail=FALSE)
    Scipy: >>> hypergeom.sf(x-1, m+n, m, k)
     
    """

    # number of genes in your query data
    k = len(query) 
    query = set(query)
    vals = []
    # background should be all genes in annotated database
    # such as go, kegg et.al.
    if isinstance(background, set): 
        bg = len(background) # total number in your annotated database 
        # filter genes that not found in annotated database
        query = query.intersection(background)
    elif isinstance(background, int):
        bg = background
    else:
        raise ValueError("backgroud should be set or int object")
    # pval
    subsets = sorted(gene_sets.keys())
    for s in subsets:
        category = gene_sets.get(s)
        m = len(category)
        hits = query.intersection(set(category))
        x = len(hits)
        if x < 1 : continue
        # pVal = hypergeom.sf(hitCount-1,popTotal,bgHits,queryTotal) 
        # p(X >= hitCounts)
        vals.append((s, hypergeom.sf(x-1, bg, m, k), x, m, hits))

    return zip(*vals)


def _ecdf(x):
    nobs = len(x)
    return np.arange(1,nobs+1)/float(nobs)


def fdrcorrection(pvals, alpha=0.05):
    """ benjamini hocheberg fdr correction. inspired by statsmodels 
    """
    # Implement copy from GOATools.
    pvals = np.asarray(pvals)
    pvals_sortind = np.argsort(pvals)
    pvals_sorted = np.take(pvals, pvals_sortind)

    ecdffactor = _ecdf(pvals_sorted)
    reject = pvals_sorted <= ecdffactor*alpha
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0])
        reject[:rejectmax] = True
    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    del pvals_corrected_raw
    pvals_corrected[pvals_corrected>1] = 1
    pvals_corrected_ = np.empty_like(pvals_corrected)
    pvals_corrected_[pvals_sortind] = pvals_corrected
    del pvals_corrected
    reject_ = np.empty_like(reject)
    reject_[pvals_sortind] = reject
    return reject_, pvals_corrected_


def multiple_testing_correction(ps, alpha=0.05, method='benjamini-hochberg', **kwargs):
    """ correct pvalues for multiple testing and add corrected `q` value

    :param ps: list of pvalues
    :param alpha: significance level default : 0.05
    :param method: multiple testing correction method [bonferroni|benjamini-hochberg]
    :returns (q, rej): two lists of q-values and rejected nodes
    """
    # Implement copy from GOATools.
    _p = np.array(ps)
    q = _p.copy()
    rej = _p.copy()
    mask = ~np.isnan(_p)
    p = _p[mask]
    if method == 'bonferroni':
        q[mask] = p * len(p)
        rej[mask] = q[mask] < alpha
    elif method == 'benjamini-hochberg':
        _rej, _q = fdrcorrection(p, alpha)
        rej[mask] = _rej
        q[mask] = _q
    else:
        raise ValueError(method)
    return q, rej