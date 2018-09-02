
import os, errno, logging
# import numpy as np
# import pandas as pd
import requests
from requests.packages.urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
# from scipy.stats import hypergeom


# CONSTANT
DEFAULT_LIBRARY=['GO_Biological_Process_2013',
				 'GO_Biological_Process_2015',
				 'GO_Cellular_Component_2013',
				 'GO_Cellular_Component_2015',
				 'GO_Molecular_Function_2013',
				 'GO_Molecular_Function_2015',
				 'GeneSigDB',
				 'HumanCyc_2015',
				 'Human_Gene_Atlas',
				 'Human_Phenotype_Ontology',
				 'Humancyc_2016',
				 'KEGG_2013',
				 'KEGG_2015',
				 'KEGG_2016',
				 'MGI_Mammalian_Phenotype_2013',
				 'MGI_Mammalian_Phenotype_Level_3',
				 'MGI_Mammalian_Phenotype_Level_4',
				 'MSigDB_Computational',
				 'MSigDB_Oncogenic_Signatures',
				 'Mouse_Gene_Atlas',
				 'Panther_2015',
				 'Panther_2016',
				 'Reactome_2013',
				 'Reactome_2015',
				 'Reactome_2016',
				 'WikiPathways_2013',
				 'WikiPathways_2015',
				 'WikiPathways_2016']

def unique(seq):
    """Remove duplicates from a list in Python while preserving order.

    :param seq: a python list object.
    :return: a list without duplicates while preserving order.

    """

    seen = set()
    seen_add = seen.add
    """
    The fastest way to sovle this problem is here
    Python is a dynamic language, and resolving seen.add each iteration
    is more costly than resolving a local variable. seen.add could have
    changed between iterations, and the runtime isn't smart enough to rule
    that out. To play it safe, it has to check the object each time.
    """

    return [x for x in seq if x not in seen and not seen_add(x)]

def mkdirs(outdir):

    try:
        os.makedirs(outdir)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise exc
        pass

def log_init(outlog, log_level=logging.INFO):
    """logging start"""

    # clear old root logger handlers
    logging.getLogger("gseapy").handlers = []
    # init a root logger
    logging.basicConfig(level    = logging.DEBUG,
                        format   = 'LINE %(lineno)-4d: %(asctime)s [%(levelname)-8s] %(message)s',
                        filename = outlog,
                        filemode = 'w')

    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(log_level)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(asctime)s %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add handlers
    logging.getLogger("gseapy").addHandler(console)
    logger = logging.getLogger("gseapy")
    # logger.setLevel(log_level)
    return logger

def log_stop(logger):
    """log stop"""

    handlers = logger.handlers[:]
    for handler in handlers:
        handler.close()
        logger.removeHandler(handler)


def retry(num=5):
    """"retry connection.
    
        define max tries num
        if the backoff_factor is 0.1, then sleep() will sleep for
        [0.1s, 0.2s, 0.4s, ...] between retries.
        It will also force a retry if the status code returned is 500, 502, 503 or 504.    
    
    """
    s = requests.Session()
    retries = Retry(total=num, backoff_factor=0.1,
                    status_forcelist=[500, 502, 503, 504])
    s.mount('http://', HTTPAdapter(max_retries=retries))

    return s


# def _ecdf(x):
#     nobs = len(x)
#     return np.arange(1,nobs+1)/float(nobs)


# def fdrcorrection(pvals, alpha=0.05):
#     """ benjamini hocheberg fdr correction. inspired by statsmodels """
#     pvals = np.asarray(pvals)
#     pvals_sortind = np.argsort(pvals)
#     pvals_sorted = np.take(pvals, pvals_sortind)
#
#     ecdffactor = _ecdf(pvals_sorted)
#     reject = pvals_sorted <= ecdffactor*alpha
#     if reject.any():
#         rejectmax = max(np.nonzero(reject)[0])
#         reject[:rejectmax] = True
#     pvals_corrected_raw = pvals_sorted / ecdffactor
#     pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
#     del pvals_corrected_raw
#     pvals_corrected[pvals_corrected>1] = 1
#     pvals_corrected_ = np.empty_like(pvals_corrected)
#     pvals_corrected_[pvals_sortind] = pvals_corrected
#     del pvals_corrected
#     reject_ = np.empty_like(reject)
#     reject_[pvals_sortind] = reject
#     return reject_, pvals_corrected_
#
#
# def calc_pvalues(nodes, query, background_attribute, M,
#                       min_category_size=3, max_category_size=500,
#         max_category_depth=5, **kwargs):
#     """ calculate pvalues for all categories in the graph
#     :param nodes: nodes dictionary from the ontology graph after background was set
#     :param query: set of identifiers for which the p value is calculated
#     :param background_attribute: node attribute assoc. with the background set
#     :param M: background size, total number of genes in the data
#     :param min_category_size: categories smaller than this number are ignored
#     :param max_category_size: categories larger than this number are ignored
#     :param max_category_depth: categories lower in the hierarchy (more specific) will be ignored
#     :returns: pvalues, x, n
#     """
#     N = len(query)
#     vals = []
#     for node in nodes:
#         category = node[background_attribute]
#         n = len(category)
#         hits = query.intersection(category)
#         x = len(hits)
#         if ((node.get('depth', 0) > max_category_depth)
#             or (n <= min_category_size)
#             or (n > max_category_size)):
#             vals.append((float('NaN'), x, n))
#         else:
#             vals.append((hypergeom.sf(x-1, M, n, N), x, n))
#     return zip(*vals)
#
#
# def multiple_testing_correction(ps, alpha=0.05,
#         method='benjamini-hochberg', **kwargs):
#     """ correct pvalues for multiple testing and add corrected `q` value
#     :param ps: list of pvalues
#     :param alpha: significance level default : 0.05
#     :param method: multiple testing correction method [bonferroni|benjamini-hochberg]
#     :returns (q, rej): two lists of q-values and rejected nodes
#     """
#     _p = np.array(ps)
#     q = _p.copy()
#     rej = _p.copy()
#     mask = ~np.isnan(_p)
#     p = _p[mask]
#     if method == 'bonferroni':
#         q[mask] = p * len(p)
#         rej[mask] = q[mask] < alpha
#     elif method == 'benjamini-hochberg':
#         _rej, _q = fdrcorrection(p, alpha)
#         rej[mask] = _rej
#         q[mask] = _q
#     else:
#         raise ValueError(method)
#     return q, rej