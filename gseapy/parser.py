# -*- coding: utf-8 -*-

import sys, logging, json, os
import requests
from io import StringIO

from numpy import in1d
from pandas import read_table, DataFrame
from requests.packages.urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
from bioservices import BioMart, BioServicesError
from gseapy.utils import unique, DEFAULT_LIBRARY, DEFAULT_CACHE_PATH, mkdirs

def gsea_cls_parser(cls):
    """Extract class(phenotype) name from .cls file.

    :param cls: the a class list instance or .cls file which is identical to GSEA input .
    :return: phenotype name and a list of class vector.
    """

    if isinstance(cls, list) :
        classes = cls
        sample_name= unique(classes)
    elif isinstance(cls, str) :
        with open(cls) as c:
            file = c.readlines()
        classes = file[2].strip('\n').split(" ")
        sample_name = file[1].lstrip("# ").strip('\n').split(" ")
    else:
        raise Exception('Error parsing sample name!')

    return sample_name[0], sample_name[1], classes

def gsea_edb_parser(results_path, index=0):
    """Parse results.edb file stored under **edb** file folder.

    :param results_path: the .results file where lcoated inside edb folder.
    :param index: gene_set index of gmt database, used for iterating items.
    :return: enrichment_term, hit_index,nes, pval, fdr.
    """
    from bs4 import BeautifulSoup

    soup = BeautifulSoup(open(results_path), features='xml')
    tag = soup.findAll('DTG')
    term = dict(tag[index].attrs)
    # dict_keys(['RANKED_LIST', 'GENESET', 'FWER', 'ES_PROFILE',
    # 'HIT_INDICES', 'ES', 'NES', 'TEMPLATE', 'RND_ES', 'RANK_SCORE_AT_ES',
    # 'NP', 'RANK_AT_ES', 'FDR'])
    enrich_term = term.get('GENESET').split("#")[1]
    es_profile = term.get('ES_PROFILE').split(" ")
    # rank_es = term.get('RND_ES').split(" ")
    hit_ind =term.get('HIT_INDICES').split(" ")
    es_profile = [float(i) for i in es_profile ]
    hit_ind = [float(i) for i in hit_ind ]
    #r ank_es = [float(i) for i in rank_es ]
    nes = term.get('NES')
    pval = term.get('NP')
    fdr =  term.get('FDR')
    # fwer = term.get('FWER')
    # index_range = len(tag)-1
    logging.debug("Enriched Gene set is: "+ enrich_term)

    return enrich_term, hit_ind, nes, pval, fdr


def gsea_gmt_parser(gmt, min_size = 3, max_size = 1000, gene_list=None):
    """Parse gene_sets.gmt(gene set database) file or download from enrichr server.

    :param gmt: the gene_sets.gmt file of GSEA input or an enrichr library name.
                checkout full enrichr library name here: http://amp.pharm.mssm.edu/Enrichr/#stats

    :param min_size: Minimum allowed number of genes from gene set also the data set. Default: 3.
    :param max_size: Maximum allowed number of genes from gene set also the data set. Default: 5000.
    :param gene_list: Used for filtering gene set. Only used this argument for :func:`call` method.
    :return: Return a new filtered gene set database dictionary.

    **DO NOT** filter gene sets, when use :func:`replot`. Because ``GSEA`` Desktop have already
    do this for you.

    """

    if gmt.lower().endswith(".gmt"):
        logging.info("User Defined gene sets is given.......continue..........")
        with open(gmt) as genesets:
             genesets_dict = { line.strip().split("\t")[0]: line.strip().split("\t")[2:]
                              for line in genesets.readlines()}
    else:
        logging.info("Downloading and generating Enrichr library gene sets...")
        if gmt in DEFAULT_LIBRARY:
            names = DEFAULT_LIBRARY
        else:
            names = get_library_name()
        if gmt in names:
            """
            define max tries num
            if the backoff_factor is 0.1, then sleep() will sleep for
            [0.1s, 0.2s, 0.4s, ...] between retries.
            It will also force a retry if the status code returned is 500, 502, 503 or 504.
            """
            s = requests.Session()
            retries = Retry(total=5, backoff_factor=0.1,
                            status_forcelist=[ 500, 502, 503, 504 ])
            s.mount('http://', HTTPAdapter(max_retries=retries))
            # queery string
            ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary'
            query_string = '?mode=text&libraryName=%s'
            # get
            response = s.get( ENRICHR_URL + query_string % gmt, timeout=None)
        else:
            raise Exception("gene_set files(.gmt) not found")
        if not response.ok:
            raise Exception('Error fetching enrichment results, check internet connection first.')

        genesets_dict = { line.strip().split("\t")[0]:
                          list(map(lambda x: x.split(",")[0], line.strip().split("\t")[2:]))
                          for line in response.iter_lines(chunk_size=1024, decode_unicode='utf-8')}



    # filtering dict
    if sys.version_info[0] >= 3 :
        genesets_filter =  {k: v for k, v in genesets_dict.items() if len(v) >= min_size and len(v) <= max_size}
    elif sys.version_info[0] == 2:
        genesets_filter =  {k: v for k, v in genesets_dict.iteritems() if len(v) >= min_size and len(v) <= max_size}
    else:
        logging.error("System failure. Please Provide correct input files")
        sys.exit(1)
    if gene_list is not None:
        subsets = sorted(genesets_filter.keys())
        for subset in subsets:
            tag_indicator = in1d(gene_list, genesets_filter.get(subset), assume_unique=True)
            tag_len = sum(tag_indicator)
            if tag_len <= min_size or tag_len >= max_size:
                del genesets_filter[subset]
            else:
                continue
    # some_dict = {key: value for key, value in some_dict.items() if value != value_to_remove}
    # use np.intersect1d() may be faster???
    filsets_num = len(genesets_dict) - len(genesets_filter)
    logging.info("%04d gene_sets have been filtered out when max_size=%s and min_size=%s"%(filsets_num, max_size, min_size))

    if filsets_num == len(genesets_dict):
        logging.error("No gene sets passed throught filtering condition!!!, try new paramters again!\n" +\
                         "Note: Gene names for gseapy is case sensitive." )
        sys.exit(1)
    else:
        return genesets_filter

def get_library_name():
    """return enrichr active enrichr library name. """

    # make a get request to get the gmt names and meta data from Enrichr
    """old code
    response = requests.get('http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=meta')
    gmt_data = response.json()
    # generate list of lib names
    libs = []
    # get library names
    for inst_gmt in gmt_data['libraries']:
        # only include active gmts
        if inst_gmt['isActive'] == True:
            libs.append(inst_gmt['libraryName'])
    """
    lib_url='http://amp.pharm.mssm.edu/Enrichr/datasetStatistics'
    libs_json = json.loads(requests.get(lib_url).text)
    libs = [lib['libraryName'] for lib in libs_json['statistics']]

    return sorted(libs)

def get_mart(dataset='hsapiens_gene_ensembl', attributes=[], 
             filters=[], filename=None, host="www.ensemble.org"):
    """mapping ids using BioMart. A wrapper of BioMart() from bioseverices 

    :param dataset: str, default: 'hsapiens_gene_ensembl'
    :param attributes: str, list, tuple
    :param filters: str, list, tuple
    :param host: www.ensemble.org, asia.ensembl.org, useast.ensembl.org
    :return: a dataframe contains all attributes you selected.

    **Note**: it will take a couple of minutes to get the results.

    How to see validated dataset, attributes, filters:
    example:

        >>> from bioservices import BioMart 
        >>> bm = BioMart(verbose=False, host="asia.ensembl.org")
        >>> b.valid_attributes ## view validated datasets, select one from dict values
        >>> bm.add_dataset_to_xml('hsapiens_gene_ensembl') # set dataset
        >>> bm.attributes('hsapiens_gene_ensembl') # view validated attrs after dataset was set
        >>> bm.filters('hsapiens_gene_ensembl') # view validated filters after dataset was set
    
    
    """
    # hosts=["www.ensemble.org", "asia.ensembl.org", "useast.ensembl.org"]
    if not attributes: 
        attributes = ['ensembl_gene_id', 'external_gene_name', 'entrezgene', 'go_id']
    
    # init
    bm = BioMart(verbose=False, host=host)
    bm.new_query()
    # 'mmusculus_gene_ensembl'
    bm.add_dataset_to_xml(dataset)
    for at in attributes:
        bm.add_attribute_to_xml(at)
    #bm.add_filter_to_xml('with_entrezgene',[])
    #bm.add_filter_to_xml('with_go',[]) # or 'go'
    #bm.add_filter_to_xml('with_hgnc',[]) # or 'hgnc_symbol'
    xml_query = bm.get_xml()
    results = bm.query(xml_query)
    df = read_table(StringIO(results), header=None, names=['gene_name', 'entrez_id','go_id'])
    # save genes with entrez and go id
    # df.dropna(inplace=True)
    # save file to cache path.
    df.dropna(subset=['entrez_id'], inplace=True)
    df['entrez'] = df.entrez.astype(int)
    df.index.name = 'gene_id'
    if filename is None: 
        mkdirs(DEFAULT_CACHE_PATH)
        filename = os.path.join(DEFAULT_CACHE_PATH, "{}.background.genes.txt".format(dataset))
    df.to_csv(filename, sep="\t")
    return df

