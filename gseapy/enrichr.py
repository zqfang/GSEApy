#!/usr/bin/env python
# python 
# see: http://amp.pharm.mssm.edu/Enrichr/help#api for API docs

import sys, json, requests, logging
from pandas import read_table
from .plot import barplot
from .utils import *

    
def enrichr(gene_list, gene_sets, description='foo', outdir='Enrichr',
            cutoff=0.05, format='pdf', figsize=(3,6), top_term=10, scale=0.8, no_plot=False):
    """Enrichr API.

    :param gene_list: Flat file with list of genes, one gene id per row.
                      or a python list object
    :param gene_sets: Enrichr Library to query.
    :param description: name of analysis
    :param outdir: Output file directory
    :param cutoff: Adjust P-value cutoff, for plotting. Default: 0.05
    :param format: Output figure format supported by matplotlib,('pdf','png','eps'...). Default: 'pdf'.
    :param no_plot: Bool, if equal to True, no figure will be draw. This is useful only if data are interested. Default: False.
    :return: A DataFrame of enrchment results, only if call ``enrichr`` inside python console.
    """
    mkdirs(outdir)

    logger = log_init(outdir, module='enrichr')
    if isinstance(gene_list, list):
        genes = [str(gene) for gene in gene_list]
        genes_str = '\n'.join(genes)
    else:
        # get gene lists
        with open(gene_list) as f:
            genes = f.read()
	
        genes_str = str(genes)
       
    # name of analysis or list
    description = str(description)
    
    #library validaty confirmationi
    gene_set = str(gene_sets)  
    
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
    outfile='Enrichr.reports.'
    url = ENRICHR_URL + query_string % (user_list_id, outfile, gene_set)
    response = requests.get(url, stream=True)

    logger.info('Downloading file of enrichment results: Job Id:'+ str(job_id))	
    with open(outdir+'/'+ outfile + description + '.txt', 'wb') as f:
        for chunk in response.iter_content(chunk_size=1024):
            if chunk:
                f.write(chunk)

    logger.debug('Results written to: ' + outfile + description + ".txt")

    df =  read_table(outdir+'/'+ outfile + description + '.txt')
    if not no_plot:
        fig = barplot(df, top_term=top_term,)
        if fig is not None:
            fig.savefig(outdir+'/'+"enrichr.reports.%s.%s"%(description, format),
                        bbox_inches='tight', dpi=300)

    if hasattr(sys, 'ps1'):
        logger.info("Enrichr: You are inside python console, a dataframe is returned.")
        logger.info("Enrichr: Job Done!")
        log_remove(logger)
        return df        
    else:
        logger.info("Enrichr: Job Done!")
    log_remove(logger)
        
    

    
