
import os, errno

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
