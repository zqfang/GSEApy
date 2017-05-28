
import os,errno,logging
from collections import OrderedDict
from pandas import DataFrame

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

def log_init(outdir, module='foo', log_level=logging.INFO):
    logging.basicConfig(
                level    = logging.DEBUG,
                format   = 'LINE %(lineno)-4d: %(asctime)s [%(levelname)-8s] %(message)s',
                filename = "%s/gseapy.%s.log"%(outdir, module),
                filemode = 'w')
    logger = logging.getLogger(__name__)
    #logger.setLevel(logging.DEBUG)

    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(log_level)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(asctime)s %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    logger.addHandler(console)
    #if you want information print to the console, uisng logger.info....
    return logger

def log_remove(logger):
    handlers = logger.handlers[:]
    for handler in handlers:
    	   handler.close()
    	   logger.removeHandler(handler)

def mkdirs(outdir):
    
    try:
        os.makedirs(outdir)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise exc
        pass   
 
def save_results(zipdata, outdir, module, gmt, data, permutation_type):
    """Save GSEA results to a csv file
    """
    res = OrderedDict()
    for gs,gseale,ind,RES in zipdata:        
        rdict = OrderedDict()      
        rdict['es'] = gseale[0]
        rdict['nes'] = gseale[1]
        rdict['pval'] = gseale[2]
        rdict['fdr'] = gseale[3]
        rdict['gene_set_size'] = len(gmt[gs])
        rdict['matched_size'] = len(ind)
        rdict['rank_ES'] = RES
        rdict['genes'] = data.iloc[ind, data.columns.get_loc('gene_name')].tolist()
        rdict['hit_index'] = ind
        res[gs] = rdict           

    res_df = DataFrame.from_dict(res, orient='index')
    res_df.index.name = 'Term'
    res_df.sort_values(by='fdr', inplace=True)

    res_df.drop(['rank_ES','hit_index'], axis=1, inplace=True)
    res_df.to_csv('{a}/gseapy.{c}.gsea.reports.csv'.format(a=outdir, b=module, c=permutation_type), float_format ='%.7f')

    return res, res_df

DEFAULT_LIBRARY =   ['Achilles_fitness_decrease', 
					 'Achilles_fitness_increase',
					 'Aging_Perturbations_from_GEO_down',
					 'Aging_Perturbations_from_GEO_up',
 					 'Allen_Brain_Atlas_down',
 					 'Allen_Brain_Atlas_up',
 					 'BioCarta_2013',
 					 'BioCarta_2015',
 					 'BioCarta_2016',
 					 'CORUM',
					 'Cancer_Cell_Line_Encyclopedia',
					 'ChEA_2013',
					 'ChEA_2015',
					 'ChEA_2016',
					 'Chromosome_Location',
					 'Disease_Perturbations_from_GEO_down',
					 'Disease_Perturbations_from_GEO_up',
					 'Disease_Signatures_from_GEO_down_2014',
					 'Disease_Signatures_from_GEO_up_2014',
					 'DrugMatrix',
					 'Drug_Perturbations_from_GEO_2014',
					 'Drug_Perturbations_from_GEO_down',
					 'Drug_Perturbations_from_GEO_up',
					 'ENCODE_Histone_Modifications_2013',
					 'ENCODE_Histone_Modifications_2015',
					 'ENCODE_TF_ChIP-seq_2014',
					 'ENCODE_TF_ChIP-seq_2015',
					 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',
					 'ESCAPE',
					 'Epigenomics_Roadmap_HM_ChIP-seq',
					 'GO_Biological_Process_2013',
					 'GO_Biological_Process_2015',
					 'GO_Cellular_Component_2013',
					 'GO_Cellular_Component_2015',
					 'GO_Molecular_Function_2013',
					 'GO_Molecular_Function_2015',
					 'GTEx_Tissue_Sample_Gene_Expression_Profiles_down',
					 'GTEx_Tissue_Sample_Gene_Expression_Profiles_up',
					 'GeneSigDB',
					 'Genes_Associated_with_NIH_Grants',
					 'Genome_Browser_PWMs',
					 'HMDB_Metabolites',
					 'HomoloGene',
					 'HumanCyc_2015',
					 'Human_Gene_Atlas',
					 'Human_Phenotype_Ontology',
					 'Humancyc_2016',
					 'KEA_2013',
					 'KEA_2015',
					 'KEGG_2013',
					 'KEGG_2015',
					 'KEGG_2016',
					 'Kinase_Perturbations_from_GEO_down',
					 'Kinase_Perturbations_from_GEO_up',
					 'LINCS_L1000_Chem_Pert_down',
					 'LINCS_L1000_Chem_Pert_up',
					 'LINCS_L1000_Kinase_Perturbations_down',
					 'LINCS_L1000_Kinase_Perturbations_up',
					 'LINCS_L1000_Ligand_Perturbations_down',
					 'LINCS_L1000_Ligand_Perturbations_up',
					 'Ligand_Perturbations_from_GEO_down',
					 'Ligand_Perturbations_from_GEO_up',
					 'MCF7_Perturbations_from_GEO_down',
					 'MCF7_Perturbations_from_GEO_up',
					 'MGI_Mammalian_Phenotype_2013',
					 'MGI_Mammalian_Phenotype_Level_3',
					 'MGI_Mammalian_Phenotype_Level_4',
					 'MSigDB_Computational',
					 'MSigDB_Oncogenic_Signatures',
					 'Microbe_Perturbations_from_GEO_down',
					 'Microbe_Perturbations_from_GEO_up',
					 'Mouse_Gene_Atlas',
					 'NCI-60_Cancer_Cell_Lines',
					 'NCI-Nature_2015',
					 'NCI-Nature_2016',
					 'NURSA_Human_Endogenous_Complexome',
					 'OMIM_Disease',
					 'OMIM_Expanded',
					 'Old_CMAP_down',
					 'Old_CMAP_up',
					 'PPI_Hub_Proteins',
					 'Panther_2015',
					 'Panther_2016',
					 'Pfam_InterPro_Domains',
					 'Phosphatase_Substrates_from_DEPOD',
					 'Reactome_2013',
					 'Reactome_2015',
					 'Reactome_2016',
					 'SILAC_Phosphoproteomics',
					 'Single_Gene_Perturbations_from_GEO_down',
					 'Single_Gene_Perturbations_from_GEO_up',
					 'TF-LOF_Expression_from_GEO',
					 'TRANSFAC_and_JASPAR_PWMs',
					 'TargetScan_microRNA',
					 'Tissue_Protein_Expression_from_Human_Proteome_Map',
					 'Tissue_Protein_Expression_from_ProteomicsDB',
					 'Transcription_Factor_PPIs',
					 'VirusMINT',
					 'Virus_Perturbations_from_GEO_down',
					 'Virus_Perturbations_from_GEO_up',
					 'WikiPathways_2013',
					 'WikiPathways_2015',
					 'WikiPathways_2016',
					 'dbGaP'] 