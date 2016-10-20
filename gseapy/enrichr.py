#!/usr/bin/env python
# python 
# see: http://amp.pharm.mssm.edu/Enrichr/help#api for API docs

import sys, os, errno, json, requests, logging
from pandas import read_table
from .plot import dotplot
from .__main__ import log_init

default_gene_set_libraries = [
    'GO_Biological_Process_2015',
    "ChEA_2015",
    "KEGG_2016",
    "ESCAPE",
    "Epigenomics_Roadmap_HM_ChIP-seq",
    "ENCODE_TF_ChIP-seq_2015",
    "ENCODE_Histone_Modifications_2015",
    "OMIM_Expanded",
    "TF-LOF_Expression_from_GEO",
    "Single_Gene_Perturbations_from_GEO_down",
    "Single_Gene_Perturbations_from_GEO_up",
    "Disease_Perturbations_from_GEO_down",
    "Disease_Perturbations_from_GEO_up",
    "Drug_Perturbations_from_GEO_down",
    "Drug_Perturbations_from_GEO_up",
    "WikiPathways_2016",
    "Reactome_2016",
    "BioCarta_2016",
    "NCI-Nature_2016"]


def get_library_name():
    """return enrichr active enrichr library name. """

	# make a get request to get the gmt names and meta data from Enrichr
	#python 2
    if sys.version_info[0] == 2 :
        import urllib2
        x = urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=meta')
        response = x.read()
        gmt_data = json.loads(response)

	# python 3
    elif sys.version_info[0] == 3:
        import urllib
        x = urllib.request.urlopen('http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=meta')
        response = x.read()
        gmt_data = json.loads(response.decode('utf-8'))
    else:
        sys.stderr.write("System failure. Please Provide correct input files")
        sys.exit(1) 
	# generate list of gmts 
    gmt_names = []

	# get library names 
    for inst_gmt in gmt_data['libraries']:

		# only include active gmts 
        if inst_gmt['isActive'] == True:
            gmt_names.append(inst_gmt['libraryName'])
    
    return sorted(gmt_names)
    
def enrichr(gene_list, gene_sets, description='foo', outdir='gseapy_out', cutoff=0.05, format='png', figsize=(3,6)):
    """Enrichr API.

    :param gene_list: flat file with list of genes, one gene id per row.
                      or a python list object, which make it easy to use 
                      inside python console.
    :param gene_sets: Enrichr Library to query.
    :param description: name of analysis
    :param outdir: out put file directory
	:param cutoff: Adjust P-value cutoff, for plotting. Default: 0.05
    :return: A DataFrame of enrchment results, only if call ``enrichr`` inside python console.
    """
    try:
        os.makedirs(outdir)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise exc
        pass

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
    enrichr_library = get_library_name()
    
    while gene_set not in enrichr_library:
        logger.error("%s is a invalid library name, Please enter a new one here!!!\n"%gene_set)
        logger.info("Use get_library_name() to veiw full list of supported names")
        gene_set = str(input())
    ## logging.debug options
    #logging.debug('Enrichr API : Input file is:', genelist)
    logger.info('Analysis name: %s, Enrichr Library: %s'%(description, gene_set))
    #logging.('Enrichr API : Enrichr Results File: ', enrichr_results)


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

    logging.debug('Job ID:'+ str(job_id))
    
    ENRICHR_URL_A = 'http://amp.pharm.mssm.edu/Enrichr/view?userListId=%s'

    user_list_id = job_id['userListId']
    response_gene_list = requests.get(ENRICHR_URL_A % str(user_list_id))

    if not response_gene_list.ok:
        raise Exception('Error getting gene list')

    logging.info('Submitted gene list:' + str(job_id))


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

    logging.debug('Get enrichment results: Job Id:'+ str(job_id))


    ## Download file of enrichment results
    #
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/export'
    query_string = '?userListId=%s&filename=%s&backgroundType=%s'
    user_list_id = str(job_id['userListId'])
    outfile='enrichr.reports.'
    url = ENRICHR_URL + query_string % (user_list_id, outfile, gene_set)
    response = requests.get(url, stream=True)

    logger.info('Downloading file of enrichment results: Job Id:'+ str(job_id))

	
    with open(outdir+'/'+ outfile + description + '.txt', 'wb') as f:
        for chunk in response.iter_content(chunk_size=1024):
            if chunk:
                f.write(chunk)

    logging.debug('Results written to: ' + outfile + description + ".txt")

    df =  read_table(outdir+'/'+ outfile + description + '.txt')
    fig = dotplot(df, cutoff=cutoff, figsize=figsize)
    if fig is not None:
        fig.savefig(outdir+'/'+"enrichr.reports.%s"%format, bbox_inches='tight', dpi=300)
		
    # convinient for viewing results inside python console. 
    #if isinstance(gene_list, list):
    if hasattr(sys, 'ps1'):
        logger.info("Enrichr: You are inside python console, a dataframe is returned.")
        logger.info("Enrichr: Job Done!")

        handlers = logger.handlers[:]
        for handler in handlers:
            handler.close()
            logger.removeHandler(handler)

        return df
        
    logger.info("Enrichr: Job Done!")

    
