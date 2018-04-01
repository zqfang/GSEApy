#!/usr/bin/env python
# -*- coding: utf-8 -*-
# see: http://amp.pharm.mssm.edu/Enrichr/help#api for API docs

import sys, json, os, logging
import requests
from time import sleep
from pandas import read_table, DataFrame, Series
from gseapy.plot import barplot
from gseapy.parser import get_library_name
from gseapy.utils import *


class Enrichr(object):
    """Enrichr API"""
    def __init__(self, gene_list, gene_sets, descriptions='foo', 
                 outdir='Enrichr', cutoff=0.05, format='pdf', 
                 figsize=(6.5,6), top_term=10, no_plot=False, verbose=False):

        self.gene_list=gene_list
        self.gene_sets=gene_sets
        self.descriptions=descriptions
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
        # init logger
        mkdirs(self.outdir)
        _gset =os.path.split(self.gene_sets)[-1].lower().rstrip(".gmt")
        outlog = os.path.join(self.outdir,"gseapy.%s.%s.log"%(self.module, _gset))
        self._logger = log_init(outlog=outlog,
                                log_level=logging.INFO if self.verbose else logging.WARNING)

    def parse_input(self):
        if isinstance(self.gene_list, list):
            genes = [str(gene) for gene in self.gene_list]
        elif isinstance(self.gene_list, DataFrame):
            # input type is bed file
            if self.gene_list.shape[1] >=3:
                genes= self.gene_list.iloc[:,:3].apply(lambda x: "\t".join([str(i) for i in x]), axis=1).tolist()
            # input type with weight values
            elif self.gene_list.shape[1] == 2:
               genes= self.gene_list.apply(lambda x: ",".join([str(i) for i in x]), axis=1).tolist()
            else:
               genes = self.gene_list.squeeze().tolist()
        elif isinstance(self.gene_list, Series):
            genes = self.gene_list.squeeze().tolist()
        else:
            # get gene lists or bed file, or gene list with weighted values.
            genes=[]
            with open(self.gene_list) as f:
                for gene in f:
                    genes.append(gene.strip())

        genes_str = '\n'.join(genes)
        return genes_str

    def run(self):
        """run enrichr"""

        # read input file
        genes_str=self.parse_input()
        
        # name of analysis or list
        description = str(self.descriptions)
        gene_set = str(self.gene_sets)

        self._logger.info("Connecting to Enrichr Server to get latest library names")
        if gene_set in DEFAULT_LIBRARY:
            enrichr_library = DEFAULT_LIBRARY
        else:
            enrichr_library = get_library_name()
            if gene_set not in enrichr_library:
                sys.stderr.write("%s is not a enrichr library name\n"%gene_set)
                sys.stdout.write("Hint: use get_library_name() to veiw full list of supported names")
                sys.exit(1)

        self._logger.info('Analysis name: %s, Enrichr Library: %s'%(description, gene_set))

        # enrichr url
        ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
        # payload
        payload = {
          'list': (None, genes_str),
          'description': (None, description)
           }
        # response
        response = requests.post(ENRICHR_URL, files=payload)
        if not response.ok:
            raise Exception('Error analyzing gene list')

        job_id = json.loads(response.text)

        self._logger.debug('Job ID:'+ str(job_id))
        ENRICHR_URL_A = 'http://amp.pharm.mssm.edu/Enrichr/view?userListId=%s'
        user_list_id = job_id['userListId']
        response_gene_list = requests.get(ENRICHR_URL_A % str(user_list_id), timeout=None)
        # wait for 1s
        sleep(1)
        if not response_gene_list.ok:
            raise Exception('Error getting gene list')

        self._logger.info('Submitted gene list:' + str(job_id))
        # Get enrichment results
        ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
        query_string = '?userListId=%s&backgroundType=%s'
        # get id data
        user_list_id = job_id['userListId']
        response = requests.get(ENRICHR_URL + query_string % (str(user_list_id), gene_set))
        if not response.ok:
            raise Exception('Error fetching enrichment results')

        self._logger.debug('Get enrichment results: Job Id:'+ str(job_id))
        # Download file of enrichment results
        ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/export'
        query_string = '?userListId=%s&filename=%s&backgroundType=%s'
        user_list_id = str(job_id['userListId'])
        filename = "%s.%s.%s.reports"%(gene_set, description, self.module)
        url = ENRICHR_URL + query_string % (user_list_id, filename, gene_set)

        # set max retries num =5
        s = retry(num=5)
        response = s.get(url, stream=True, timeout=None)

        self._logger.info('Downloading file of enrichment results: Job Id:'+ str(job_id))
        outfile="%s/%s.%s.%s.reports.txt"%(self.outdir, gene_set, description, self.module)

        with open(outfile, 'wb') as f:
            for chunk in response.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)

        self._logger.debug('Results written to: ' + outfile)
        # save results
        df =  read_table(outfile)
        self.res2d = df

        # plotting
        if not self.__no_plot:
            fig = barplot(df=df, cutoff=self.cutoff,
                          figsize=self.figsize, 
                          top_term=self.__top_term,
                          color='salmon',
                          title='')
            if fig is None:
                self._logger.warning("Warning: No enrich terms using library %s when cuttoff = %s"%(gene_set, self.cutoff))
            else:
                fig.savefig(outfile.replace("txt", self.format),
                            bbox_inches='tight', dpi=300)
        self._logger.info('Done.\n')
        return


def enrichr(gene_list, gene_sets, description='foo', outdir='Enrichr',
            cutoff=0.05, format='pdf', figsize=(8,6), top_term=10, no_plot=False, verbose=False):
    """Enrichr API.

    :param gene_list: Flat file with list of genes, one gene id per row, or a python list object
    :param gene_sets: Enrichr Library to query. Required enrichr library name
    :param description: name of analysis. optional.
    :param outdir: Output file directory
    :param float cutoff: Adjust P-value cutoff, for plotting. Default: 0.05
    :param str format: Output figure format supported by matplotlib,('pdf','png','eps'...). Default: 'pdf'.
    :param list figsize: Matplotlib figsize, accept a tuple or list, e.g. (width,height). Default: (6.5,6).
    :param bool no_plot: if equal to True, no figure will be draw. Default: False.
    :param bool verbose: Increase output verbosity, print out progress of your job, Default: False.

    :return: An Enrichr object, which obj.res2d contains your enrichr query.
    """
    enr = Enrichr(gene_list, gene_sets, description, outdir,
                  cutoff, format, figsize, top_term, no_plot, verbose)
    enr.run()

    return enr
