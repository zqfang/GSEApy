#!/usr/bin/env python
# python 
# see: http://amp.pharm.mssm.edu/Enrichr/help#api for API docs


from __future__ import print_function
import json
import requests
import sys


def enrichr(gene_list, description, enrichr_library, outdir):
    """Enrichr API.

    :param gene_list: flat file with list of genes, one gene id per row
    :param description: name of analysis
    :param enrich_library: Enrichr Library to query.
    :param outdir: out put file prefix
    
    """
    
    genelist = gene_list
    list_desrciption = description
    enrichr_library = enrichr_library
    enrichr_results = outdir


    ## Print options
    #print('Enrichr API : Input file is:', genelist)
    print('Enrichr API : Analysis name: ', list_desrciption)
    print('Enrichr API : Enrichr Library: ', enrichr_library)
    #print('Enrichr API : Enrichr Results File: ', enrichr_results)

    # get gene lits
    #with open(genelist) as f:
    #    genes = f.read()
 
    ## enrichr url
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'

    # stick gene list here
    #genes_str = str(genes)
    genes_str = genelist
    # genes_str = '\n'.join(genelist)

    # name of analysis or list
    description = str(list_desrciption)

    # payload
    payload = {
      'list': (None, genes_str),
      'description': (None, description)
       }   

    # response
    #print("Enrichr API : requests.post")
    response = requests.post(ENRICHR_URL, files=payload)

    if not response.ok:
        raise Exception('Error analyzing gene list')

    job_id = json.loads(response.text)

    print('Enrichr API : Job ID:', job_id)

    ################################################################################
    # View added gene list
    #
    ENRICHR_URL_A = 'http://amp.pharm.mssm.edu/Enrichr/view?userListId=%s'

    user_list_id = job_id['userListId']
    #print(user_list_id)

    response_gene_list = requests.get(ENRICHR_URL_A % str(user_list_id))

    if not response_gene_list.ok:
        raise Exception('Error getting gene list')

    print('Enrichr API : View added gene list:', job_id)
    added_gene_list = json.loads(response_gene_list.text)
    #print(added_gene_list)

    ################################################################################
    # Get enrichment results
    #
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'

    ## get id data
    user_list_id = job_id['userListId']

    ## Libraray
    gene_set_library = str(enrichr_library)

    response = requests.get(
        ENRICHR_URL + query_string % (str(user_list_id), gene_set_library)
          )
    if not response.ok:
        raise Exception('Error fetching enrichment results')

    print('Enrichr API : Get enrichment results: Job Id:', job_id)
    data = json.loads(response.text)
    #print(data)

    ################################################################################
    ## Download file of enrichment results
    #
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/export'

    query_string = '?userListId=%s&filename=%s&backgroundType=%s'

    user_list_id = str(job_id['userListId'])

    filename = enrichr_results

    gene_set_library = str(enrichr_library)

    url = ENRICHR_URL + query_string % (user_list_id, filename, gene_set_library)

    response = requests.get(url, stream=True)

    print('Enrichr API : Downloading file of enrichment results: Job Id:', job_id)
    with open(filename + '.txt', 'wb') as f:
        for chunk in response.iter_content(chunk_size=1024):
            if chunk:
                f.write(chunk)
    ################################################
    print('Enrichr API : Results written to:', enrichr_results + ".txt")
    print("Enrichr API : Done")

