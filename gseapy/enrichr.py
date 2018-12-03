#!/usr/bin/env python
# -*- coding: utf-8 -*-
# see: http://amp.pharm.mssm.edu/Enrichr/help#api for API docs

import sys, json, os, logging
import requests
import pandas as pd
from io import StringIO
from collections import OrderedDict
from functools import reduce
from pkg_resources import resource_filename
from time import sleep
from tempfile import TemporaryDirectory
from gseapy.plot import barplot
from gseapy.parser import get_library_name, Biomart
from gseapy.utils import *
from gseapy.stats import calc_pvalues, multiple_testing_correction


class Enrichr(object):
    """Enrichr API"""
    def __init__(self, gene_list, gene_sets, descriptions='', outdir='Enrichr', 
                 cutoff=0.05, background='hsapiens_gene_ensembl', 
                 format='pdf', figsize=(6.5,6), top_term=10, no_plot=False, 
                 verbose=False):

        self.gene_list=gene_list
        self.gene_sets=gene_sets
        self.descriptions=str(descriptions)
        self.outdir=outdir
        self.cutoff=cutoff
        self.format=format
        self.figsize=figsize
        self.__top_term=int(top_term)
        self.__no_plot=no_plot
        self.verbose=bool(verbose)
        self.module="enrichr"
        self.res2d=None
        self._processes=1
        self.background=background
        self._bg=None
        # init logger
        logfile = self.prepare_outdir()
        self._logger = log_init(outlog=logfile,
                                log_level=logging.INFO if self.verbose else logging.WARNING)

    def prepare_outdir(self):
        """create temp directory."""
        self._outdir = self.outdir
        if self._outdir is None:
            self._tmpdir = TemporaryDirectory()
            self.outdir = self._tmpdir.name
        elif isinstance(self.outdir, str):
            mkdirs(self.outdir)
        else:
            raise Exception("Error parsing outdir: %s"%type(self.outdir))

        # handle gene_sets
        logfile = os.path.join(self.outdir, "gseapy.%s.%s.log" % (self.module, self.descriptions))
        return logfile

    def parse_genesets(self):
        """parse gene_sets input file type"""

        enrichr_library = get_library_name()
        if isinstance(self.gene_sets, list):
            gss = self.gene_sets
        elif isinstance(self.gene_sets, str):
            gss = [ g.strip() for g in self.gene_sets.strip().split(",") ]
        elif isinstance(self.gene_sets, dict):
            gss = [self.gene_sets]
        else:
            raise Exception("Error parsing enrichr libraries, please provided corrected one")
        
        # gss: a list contain .gmt, dict, enrichr_liraries.
        # now, convert .gmt to dict
        gss_exist = [] 
        for g in gss:
            if isinstance(g, dict): 
                gss_exist.append(g)
                continue

            if isinstance(g, str): 
                if g in enrichr_library: 
                    gss_exist.append(g)
                    continue
                if g.lower().endswith(".gmt") and os.path.exists(g):
                    self._logger.info("User Defined gene sets is given: %s"%g)
                    with open(g) as genesets:
                        g_dict = { line.strip().split("\t")[0]: line.strip().split("\t")[2:]
                                        for line in genesets.readlines() }
                    gss_exist.append(g_dict)
        return gss_exist

    def parse_genelists(self):
        """parse gene list"""
        if isinstance(self.gene_list, list):
            genes = self.gene_list
        elif isinstance(self.gene_list, pd.DataFrame):
            # input type is bed file
            if self.gene_list.shape[1] >=3:
                genes= self.gene_list.iloc[:,:3].apply(lambda x: "\t".join([str(i) for i in x]), axis=1).tolist()
            # input type with weight values
            elif self.gene_list.shape[1] == 2:
               genes= self.gene_list.apply(lambda x: ",".join([str(i) for i in x]), axis=1).tolist()
            else:
               genes = self.gene_list.squeeze().tolist()
        elif isinstance(self.gene_list, pd.Series):
            genes = self.gene_list.squeeze().tolist()
        else:
            # get gene lists or bed file, or gene list with weighted values.
            genes=[]
            with open(self.gene_list) as f:
                for gene in f:
                    genes.append(gene.strip())

        self._isezid = all(map(self._is_entrez_id, genes))
        if self._isezid: 
            self._gls = set(map(int, self._gls))
        else:
            self._gls = genes

        return '\n'.join(genes)

    def send_genes(self, gene_list, url):
        """ send gene list to enrichr server"""
        payload = {
          'list': (None, gene_list),
          'description': (None, self.descriptions)
           }
        # response
        response = requests.post(url, files=payload)
        if not response.ok:
            raise Exception('Error analyzing gene list')
        sleep(1)
        job_id = json.loads(response.text)

        return job_id

    def check_genes(self, gene_list, usr_list_id):
        '''
        Compare the genes send and received to get succesfully recognized genes
        '''
        response = requests.get('http://amp.pharm.mssm.edu/Enrichr/view?userListId=%s' % usr_list_id)
        if not response.ok:
            raise Exception('Error getting gene list back')
        returnedL = json.loads(response.text)["genes"]
        returnedN = sum([1 for gene in gene_list if gene in returnedL])
        self._logger.info('{} genes successfully recognized by Enrichr'.format(returnedN))

    def get_results(self, gene_list):
        """Enrichr API"""
        ADDLIST_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
        # RESULTS_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
        # query_string = '?userListId=%s&backgroundType=%s'
        job_id = self.send_genes(gene_list, ADDLIST_URL)
        user_list_id = job_id['userListId']

        RESULTS_URL = 'http://amp.pharm.mssm.edu/Enrichr/export'
        query_string = '?userListId=%s&filename=%s&backgroundType=%s'
        # set max retries num =5
        s = retry(num=5)
        filename = "%s.%s.reports" % (self._gs, self.descriptions)
        url = RESULTS_URL + query_string % (user_list_id, filename, self._gs)
        response = s.get(url, stream=True, timeout=None)
        # response = requests.get(RESULTS_URL + query_string % (user_list_id, gene_set))
        sleep(1)
        res = pd.read_table(StringIO(response.content.decode('utf-8')))
        return [job_id['shortId'], res]

    def _is_entrez_id(self, idx):
        try:
            int(idx)
            return True
        except:
            return False   

    def get_background(self):
        """get background gene"""
        DB_FILE = resource_filename("gseapy", "data/{}.background.genes.txt".format(self.background))
        filename = os.path.join(DEFAULT_CACHE_PATH, "{}.background.genes.txt".format(self.background))
        if os.path.exists(filename):
            df = pd.read_table(filename)
        elif os.path.exists(DB_FILE):
            df = pd.read_table(DB_FILE)
        else:
            self._logger.warning("Downloading %s for the first time. It might take a couple of miniutes."%self.background)
            bm = Biomart()
            df = bm.query(dataset=self.background)
            df.dropna(subset=['go_id'], inplace=True)
        self._logger.info("using all annotated genes with GO_ID as background genes")
        df.dropna(subset=['entrezgene'], inplace=True)     

        return df

    def enrich(self, gmt):
        """use local mode
         
        p = p-value computed using the Fisher exact test (Hypergeometric test)  

        Not implemented here:

            combine score = log(p)·z

        see here: http://amp.pharm.mssm.edu/Enrichr/help#background&q=4
        
        columns contain:
            
            Term Overlap P-value Adjusted_P-value Genes

        """
        if isinstance(self.background, str): 
            # self.background = set(reduce(lambda x,y: x+y, gmt.values(),[]))
            df = self.get_background()
            # input id type: entrez or gene_name
            if self._isezid:
                bg = df['entrezgene'].astype(int)
            else:
                bg = df['external_gene_name']

            self._bg = set(bg.unique())
            self._logger.warning("Background: %s %s genes with GO_IDs. "%(self._bg, self.background))
            self._logger.warning("If this is not you wanted, please give a number to background argument")
        hgtest = list(calc_pvalues(query=self._gls, gene_sets=gmt, 
                                   background=self._bg))
        if len(hgtest) > 0:
            terms, pvals, olsz, gsetsz, genes = hgtest
            fdrs, rej = multiple_testing_correction(ps = pvals, 
                                                    alpha=self.cutoff,
                                                    method='benjamini-hochberg')
            # save to a dataframe
            odict = OrderedDict()
            odict['Term'] = terms
            odict['Overlap'] = list(map(lambda h,g: "%s/%s"%(h, g), olsz, gsetsz))
            odict['P-value'] = pvals
            odict['Adjusted P-value'] = fdrs
            # odict['Reject (FDR< %s)'%self.cutoff ] = rej
            odict['Genes'] = [";".join(g) for g in genes]
            res = pd.DataFrame(odict)
            return res
        return 

    def run(self):
        """run enrichr for one sample gene list but multi-libraries"""

        # read input file
        genes_list = self.parse_genelists()
        gss = self.parse_genesets()
        # if gmt
        self._logger.info("Connecting to Enrichr Server to get latest library names")
        if len(gss) < 1:
            sys.stderr.write("Not validated Enrichr library name provided\n")
            sys.stdout.write("Hint: use get_library_name() to view full list of supported names")
            sys.exit(1)
        self.results = pd.DataFrame()

        for g in gss: 
            if isinstance(g, dict): 
                ## local mode
                res = self.enrich(g)
                shortID, self._gs = str(id(g)), "CUSTOM%s"%id(g)
                if res is None: 
                    self._logger.info("No hits return, for gene set: Custom%s"%shortID)
                    continue
            else:
                ## online mode
                self._gs = str(g)
                self._logger.debug("Start Enrichr using library: %s" % (self._gs))
                self._logger.info('Analysis name: %s, Enrichr Library: %s' % (self.descriptions, self._gs))
                shortID, res = self.get_results(genes_list)
                # Remember gene set library used
            res.insert(0, "Gene_set", self._gs)
            # Append to master dataframe
            self.results = self.results.append(res, ignore_index=True)
            self.res2d = res
            if self._outdir is None: continue
            self._logger.info('Save file of enrichment results: Job Id:' + str(shortID))
            outfile = "%s/%s.%s.%s.reports.txt" % (self.outdir, self._gs, self.descriptions, self.module)
            self.res2d.to_csv(outfile, index=False, encoding='utf-8', sep="\t")
            # plotting
            if not self.__no_plot:
                msg = barplot(df=res, cutoff=self.cutoff, figsize=self.figsize,
                              top_term=self.__top_term, color='salmon',
                              title=self._gs,
                              ofname=outfile.replace("txt", self.format))
                if msg is not None : self._logger.warning(msg)
            self._logger.info('Done.\n')
        # clean up tmpdir
        if self._outdir is None: self._tmpdir.cleanup()

        return


def enrichr(gene_list, gene_sets, description='', outdir='Enrichr', cutoff=0.05, 
            background='hsapiens_gene_ensembl', format='pdf', 
            figsize=(8,6), top_term=10, no_plot=False, verbose=False):
    """Enrichr API.

    :param gene_list: Flat file with list of genes, one gene id per row, or a python list object
    :param gene_sets: Enrichr Library to query. Required enrichr library name(s). Separate each name by comma.
    :param description: name of analysis. optional.
    :param outdir: Output file directory
    :param float cutoff: Adjust P-value (benjamini-hochberg correction)cutoff. Default: 0.05
    :param int background: BioMart dataset name for retrieving background gene information.
                           This argument only works when gene_sets input is a gmt file or python dict.
                           You could also specify a number by yourself, e.g. total expressed genes number.
                           In this case, you will skip retrieving background infos from biomart.
    
    use the code below to see validated background dataset name from BioMart.
    Here are examle code:
    >>> from gseapy.parser import Biomart 
    >>> bm = Biomart(verbose=False, host="asia.ensembl.org")
    >>> ## view validated marts
    >>> marts = bm.get_marts()
    >>> ## view validated dataset
    >>> datasets = bm.get_datasets(mart='ENSEMBL_MART_ENSEMBL')

    :param str format: Output figure format supported by matplotlib,('pdf','png','eps'...). Default: 'pdf'.
    :param list figsize: Matplotlib figsize, accept a tuple or list, e.g. (width,height). Default: (6.5,6).
    :param bool no_plot: if equal to True, no figure will be draw. Default: False.
    :param bool verbose: Increase output verbosity, print out progress of your job, Default: False.

    :return: An Enrichr object, which obj.res2d stores your last query, obj.results stores your all queries.
    
    """
    enr = Enrichr(gene_list, gene_sets, description, outdir,
                  cutoff, background, format, figsize, top_term, no_plot, verbose)
    enr.run()

    return enr

