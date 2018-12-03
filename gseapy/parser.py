# -*- coding: utf-8 -*-

import sys, logging, json, os
import requests
import pandas as pd
from io import StringIO
from numpy import in1d
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


class Biomart(BioMart):
    """query from BioMart"""
    def __init__(self, host="www.ensembl.org", verbose=False):
        """A wrapper of BioMart() from bioseverices.

        How to query validated dataset, attributes, filters:
        example:
        >>> from gseapy.parser import Biomart 
        >>> bm = Biomart(verbose=False, host="asia.ensembl.org")
        >>> ## view validated marts
        >>> marts = bm.get_marts()
        >>> ## view validated dataset
        >>> datasets = bm.get_datasets(mart='ENSEMBL_MART_ENSEMBL')
        >>> ## view validated attributes
        >>> attrs = bm.get_attributes(dataset='hsapiens_gene_ensembl') 
        >>> ## view validated filters
        >>> filters = bm.get_filters(dataset='hsapiens_gene_ensembl')
        >>> ## query results
        >>> results = bm.query(dataset='hsapiens_gene_ensembl', 
                               attributes=['entrezgene', ‘go_id'],
                               filters={'ensemble_gene_id': [your input list]}
                              )         
        """
        super(Biomart, self).__init__(host=host, verbose=verbose)
        hosts=["www.ensemble.org", "asia.ensembl.org", "useast.ensembl.org"]
        # if host not work, select next
        i=0
        while (self.host is None) and (i < 3):
            self.host = hosts[i]
            i +=1 
     
    def get_marts(self):
        """Get available marts and their names."""

        mart_names = pd.Series(self.names, name="Name")
        mart_descriptions = pd.Series(self.displayNames, name="Description")

        return pd.concat([mart_names, mart_descriptions], axis=1)

    def get_datasets(self, mart='ENSEMBL_MART_ENSEMBL'):
        """Get available datasets from mart you've selected"""
        datasets = self.datasets(mart, raw=True)
        return pd.read_table(StringIO(datasets), header=None, usecols=[1, 2],
                            names = ["Name", "Description"])

    def get_attributes(self, dataset):
        """Get available attritbutes from dataset you've selected"""
        attributes = self.attributes(dataset)
        attr_ = [ (k, v[0]) for k, v in attributes.items()]
        return pd.DataFrame(attr_, columns=["Attribute","Description"])

    def get_filters(self, dataset):
        """Get available filters from dataset you've selected"""
        filters = self.filters(dataset)
        filt_ = [ (k, v[0]) for k, v in filters.items()]
        return pd.DataFrame(filt_, columns=["Filter", "Description"])
    
    def query(self, dataset='hsapiens_gene_ensembl', attributes=[], 
              filters={}, filename=None):
        """mapping ids using BioMart.  

        :param dataset: str, default: 'hsapiens_gene_ensembl'
        :param attributes: str, list, tuple
        :param filters: dict, {'filter name': list(filter value)}
        :param host: www.ensemble.org, asia.ensembl.org, useast.ensembl.org
        :return: a dataframe contains all attributes you selected.

        **Note**: it will take a couple of minutes to get the results.
        A xml template for querying biomart. (see https://gist.github.com/keithshep/7776579)
        
        exampleTaxonomy = "mmusculus_gene_ensembl"
        exampleGene = "ENSMUSG00000086981,ENSMUSG00000086982,ENSMUSG00000086983"
        urlTemplate = \
        '''http://ensembl.org/biomart/martservice?query=''' \
        '''<?xml version="1.0" encoding="UTF-8"?>''' \
        '''<!DOCTYPE Query>''' \
        '''<Query virtualSchemaName="default" formatter="CSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.6">''' \
        '''<Dataset name="%s" interface="default"><Filter name="ensembl_gene_id" value="%s"/>''' \
        '''<Attribute name="ensembl_gene_id"/><Attribute name="ensembl_transcript_id"/>''' \
        '''<Attribute name="transcript_start"/><Attribute name="transcript_end"/>''' \
        '''<Attribute name="exon_chrom_start"/><Attribute name="exon_chrom_end"/>''' \
        '''</Dataset>''' \
        '''</Query>''' 
        
        exampleURL = urlTemplate % (exampleTaxonomy, exampleGene)
        req = requests.get(exampleURL, stream=True)
                   
        """
        if not attributes: 
            attributes = ['ensembl_gene_id', 'external_gene_name', 'entrezgene', 'go_id'] 
        # i=0
        # while (self.host is None) and (i < 3):
        #     self.host = self.ghosts[i]
        #     i +=1 
        self.new_query()
        # 'mmusculus_gene_ensembl'
        self.add_dataset_to_xml(dataset)
        for at in attributes:
            self.add_attribute_to_xml(at)
        # add filters
        if filters:
            for k, v in filters.items(): 
                if isinstance(v, list): v = ",".join(v)
                self.add_filter_to_xml(k, v)

        xml_query = self.get_xml()
        results = super(Biomart, self).query(xml_query)
        df = pd.read_table(StringIO(results), header=None, names=attributes, index_col=None)
        # save file to cache path.
        if filename is None: 
            mkdirs(DEFAULT_CACHE_PATH)
            filename = os.path.join(DEFAULT_CACHE_PATH, "{}.background.genes.txt".format(dataset))
        df.to_csv(filename, sep="\t", index=False)
      
        return df
