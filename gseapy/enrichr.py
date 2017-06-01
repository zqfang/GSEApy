#!/usr/bin/env python
# python 
# see: http://amp.pharm.mssm.edu/Enrichr/help#api for API docs

import sys, json, requests, logging
from pandas import read_table
from gseapy.gsea import GSEAbase
from gseapy.plot import barplot
from gseapy.utils import *

class Enrichr(GSEAbase):
    """Enrichr API"""
    def __init__(self, gene_list, gene_sets, descriptions='foo', outdir='Enrichr',
            cutoff=0.05, format='pdf', figsize=(12,6), top_term=10, no_plot=False, verbose=False):

        self.gene_list=gene_list
        self.gene_sets=gene_sets
        self.descriptions=descriptions
        self.outdir=outdir
        self.cutoff=cutoff
        self.format=format
        self.figsize=figsize
        self.__top_term=top_term
        self.__no_plot=no_plot
        self.verbose=verbose
        self.module="enrichr"
        self.res2d=None
        self._processes=1
    def run(self):
        """run enrichr"""

        mkdirs(self.outdir)

        if isinstance(self.gene_list, list):
            genes = [str(gene) for gene in self.gene_list]
            genes_str = '\n'.join(genes)
        else:
            # get gene lists
            with open(self.gene_list) as f:
                genes = f.read()
        
            genes_str = str(genes)

        # name of analysis or list
        description = str(self.descriptions)
        
        #library validaty confirmationi
        gene_set = str(self.gene_sets) 
        #logging start
        logger = self._log_init(module=self.module,  
                                log_level=logging.INFO if self.verbose else logging.WARNING)
    
        logger.info("Connecting to Enrichr Server to get latest library names")
        if gene_set in DEFAULT_LIBRARY:
            enrichr_library = DEFAULT_LIBRARY
        else:
            enrichr_library = get_library_name()
            if gene_set not in enrichr_library:
                sys.stderr.write("%s is not a enrichr library name\n"%gene_set)
                sys.stdout.write("Hint: use get_library_name() to veiw full list of supported names")
                sys.exit(1)
            
        logger.info('Analysis name: %s, Enrichr Library: %s'%(description, gene_set))

        ## enrichr url
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

        logger.debug('Job ID:'+ str(job_id))   
        ENRICHR_URL_A = 'http://amp.pharm.mssm.edu/Enrichr/view?userListId=%s'
        user_list_id = job_id['userListId']
        response_gene_list = requests.get(ENRICHR_URL_A % str(user_list_id))

        if not response_gene_list.ok:
            raise Exception('Error getting gene list')

        logger.info('Submitted gene list:' + str(job_id))
        # Get enrichment results
        ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
        query_string = '?userListId=%s&backgroundType=%s'
        ## get id data
        user_list_id = job_id['userListId']
        response = requests.get(
            ENRICHR_URL + query_string % (str(user_list_id), gene_set)
              )
        if not response.ok:
            raise Exception('Error fetching enrichment results')

        logger.debug('Get enrichment results: Job Id:'+ str(job_id))
        ## Download file of enrichment results
        ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/export'
        query_string = '?userListId=%s&filename=%s&backgroundType=%s'
        user_list_id = str(job_id['userListId'])
        filename = "%s.%s.%s.reports"%(gene_set, description, self.module)
        
        url = ENRICHR_URL + query_string % (user_list_id, filename, gene_set)
        response = requests.get(url, stream=True)

        logger.info('Downloading file of enrichment results: Job Id:'+ str(job_id)) 
        outfile="%s/%s.%s.%s.reports.txt"%(self.outdir, gene_set, description, self.module)

        with open(outfile, 'wb') as f:
            for chunk in response.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)

        logger.debug('Results written to: ' + outfile)
        #save results
        df =  read_table(outfile)
        self.res2d = df

        #plotting
        if not self.__no_plot:
            fig = barplot(df=df, cutoff=self.cutoff, 
                        figsize=self.figsize, top_term=self.__top_term,)
            if fig is None:
                logger.warning("Warning: No enrich terms using library %s when cuttoff = %s"%(gene_set, self.cutoff))
            else:
                fig.savefig(outfile.replace("txt", self.format),
                            bbox_inches='tight', dpi=300)
        
        self._log_stop()

        return 
def enrichr(gene_list, gene_sets, description='foo', outdir='Enrichr',
            cutoff=0.05, format='pdf', figsize=(12,6), top_term=10, no_plot=False, verbose=False):
    """Enrichr API.

    :param gene_list: Flat file with list of genes, one gene id per row.
                      or a python list object
    :param gene_sets: Enrichr Library to query.
    :param description: name of analysis
    :param outdir: Output file directory
    :param cutoff: Adjust P-value cutoff, for plotting. Default: 0.05
    :param format: Output figure format supported by matplotlib,('pdf','png','eps'...). Default: 'pdf'.
    :param no_plot: Bool, if equal to True, no figure will be draw. This is useful only if data are interested. Default: False.
    :param verbose: Bool, increase output verbosity, print out progress of your job, Default: False.
    :return: A DataFrame of enrchment results, only if call ``enrichr`` inside python console.
    """
    enr = Enrichr(gene_list, gene_sets, description, outdir,
                  cutoff, format, figsize, top_term, no_plot, verbose)
    enr.run()

    return enr

    
