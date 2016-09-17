#!/usr/bin/env python
# python 
# see: http://amp.pharm.mssm.edu/Enrichr/help#api for API docs


from __future__ import print_function


import json
import requests
import sys, os
from pandas import read_table
from .plot import dotplot

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


def get_libary_name():
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
        print("System failure. Please Provide correct input files")
        sys.exit(1) 
	# generate list of gmts 
    gmt_names = []

	# get library names 
    for inst_gmt in gmt_data['libraries']:

		# only include active gmts 
        if inst_gmt['isActive'] == True:
            gmt_names.append(inst_gmt['libraryName'])


    return sorted(gmt_names)
    
def enrichr(gene_list, description, gene_sets, outdir, cutoff=0.05, format='png', figsize=(3,6)):
    """Enrichr API.

    :param gene_list: flat file with list of genes, one gene id per row.
                      or a python list object, which make it easy to use 
                      inside python console.
    :param description: name of analysis
    :param gene_set: Enrichr Library to query.
    :param outdir: out put file directory
	:param cutoff: Adjust P-value cutoff, for plotting. Default: 0.05
    :return: A DataFrame of enrchment results, only if call ``enrichr`` inside python console.
    """
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
    
    enrichr_library = get_libary_name()
    
    while gene_set not in enrichr_library:
        print("You have choose a invalidity library, Please enter a worked one here!!!\n" )
        print(enrichr_library,"\n")
        gene_set = str(input())
    ## Print options
    #print('Enrichr API : Input file is:', genelist)
    print('Enrichr API : Analysis name: ', description)
    print('Enrichr API : Enrichr Library: ', gene_set)
    #print('Enrichr API : Enrichr Results File: ', enrichr_results)


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

    print('Enrichr API : Job ID:', job_id)
    
    ENRICHR_URL_A = 'http://amp.pharm.mssm.edu/Enrichr/view?userListId=%s'

    user_list_id = job_id['userListId']
    response_gene_list = requests.get(ENRICHR_URL_A % str(user_list_id))

    if not response_gene_list.ok:
        raise Exception('Error getting gene list')

    print('Enrichr API : Submitted gene list:', job_id)


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

    print('Enrichr API : Get enrichment results: Job Id:', job_id)


    ## Download file of enrichment results
    #
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/export'
    query_string = '?userListId=%s&filename=%s&backgroundType=%s'
    user_list_id = str(job_id['userListId'])
    outfile='enrichr.reports'
    url = ENRICHR_URL + query_string % (user_list_id, outfile, gene_set)
    response = requests.get(url, stream=True)

    print('Enrichr API : Downloading file of enrichment results: Job Id:', job_id)
    os.system("mkdir "+ outdir)
    with open(outdir+'/'+ outfile + description + '.txt', 'wb') as f:
        for chunk in response.iter_content(chunk_size=1024):
            if chunk:
                f.write(chunk)

    print('Enrichr API : Results written to:', outfile + description + ".txt")

    df =  read_table(outdir+'/'+outfile + '.txt')
    fig = dotplot(df, cutoff=cutoff, figsize=figsize)
    if fig is not None:
        fig.savefig(outdir+'/'+"enrichr.reports.%s"%format, bbox_inches='tight', dpi=300)
		
    # convinient for viewing results inside python console. 
    if isinstance(gene_list, list):
        print("Enrichr API : You are seeing this message, because you are inside python console.\n"+\
              "Enrichr API : It will return a pandas dataframe for veiwing results."  )
        print("Enrichr API : Done")

        return df
        
    print("Enrichr API : Done")

    
