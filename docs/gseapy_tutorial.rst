.. _tutorial:

======================================
A Protocol to Prepare files for GSEApy
======================================

As a biological researcher, I like protocols.

Here is a short tutorial for you to walk you through gseapy.

For file format explanation, please see `here <http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html.>`_

In order to run gseapy successfully, install gseapy use pip.

.. code:: bash

    pip install gseapy

    # if you have conda
    conda install -c bioconda gseapy



Use ``gsea`` command, or :func:`gsea`
================================================

Follow the steps blow.

One thing you should know is that the gseapy input files are the same as
``GSEA`` desktop required. You can use these files below to run ``GSEA`` desktop, too.


1. Prepare an tabular text file of gene expression like this:
------------------------------------------------------------------

**RNA-seq,ChIP-seq, Microarry data** are all supported.

Here is to see what the structure of expression table looks like

.. code:: python

    import pandas as pd
    df = pd.read_table('./test/gsea_data.txt')
    df.head()

    #or assign dataframe to the parameter 'data'


.. raw:: html

    <div>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>Gene_Symbol</th>
          <th>ACD2</th>
          <th>BCD2</th>
          <th>CCD2</th>
          <th>APYD2</th>
          <th>BPYD2</th>
          <th>CPYD2</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>SCARA3</td>
          <td>6.287</td>
          <td>6.821</td>
          <td>6.005</td>
          <td>2.525</td>
          <td>1.911</td>
          <td>1.500</td>
        </tr>
        <tr>
          <th>1</th>
          <td>POU5F1</td>
          <td>11.168</td>
          <td>11.983</td>
          <td>10.469</td>
          <td>7.795</td>
          <td>7.911</td>
          <td>6.174</td>
        </tr>
        <tr>
          <th>2</th>
          <td>CTLA2B</td>
          <td>4.362</td>
          <td>5.708</td>
          <td>4.633</td>
          <td>1.493</td>
          <td>0.000</td>
          <td>1.369</td>
        </tr>
        <tr>
          <th>3</th>
          <td>CRYAB</td>
          <td>11.339</td>
          <td>11.662</td>
          <td>11.714</td>
          <td>7.698</td>
          <td>7.928</td>
          <td>7.779</td>
        </tr>
        <tr>
          <th>4</th>
          <td>PMP22</td>
          <td>7.259</td>
          <td>7.548</td>
          <td>6.803</td>
          <td>4.418</td>
          <td>2.239</td>
          <td>3.071</td>
        </tr>
      </tbody>
    </table>
    </div>





2. An cls file is also expected.
-----------------------------------------------

This file is used to specify column attributes in step 1, just like ``GSEA`` asked.

An example of cls file looks like below.

.. code:: python

    with open('gsea/edb/C1OE.cls') as cls:
        print(cls.read())

    # or assign a list object to parameter 'cls' like this
    # cls=['C1OE', 'C1OE', 'C1OE', 'Vector', 'Vector', 'Vector']

.. parsed-literal::

    6 2 1
    # C1OE Vector
    C1OE C1OE C1OE Vector Vector Vector


| The first line specify the total samples and phenotype numbers. Leave number 1 always be 1.
| The second line specify the phenotype class(name).
| The third line specify column attributes in step 1.





3. Gene_sets file in gmt format.
-----------------------------------------------------

All you need to do is to download gene set database file from ``GSEA`` website.

Or you could use enrichr library. In this case, just provide library name to parameter 'gene_sets'

If you would like to use you own gene_sets.gmts files, build such a file use excel:


An example of gmt file looks like below:


.. code:: python

    with open('gsea/edb/gene_sets.gmt') as gmt:
        print(gmt.read())


.. parsed-literal::

    ES-SPECIFIC	Arid3a_used	ACTA1	CALML4	CORO1A	DHX58	DPYS	EGR1	ESRRB	GLI2	GPX2	HCK	INHBB
    HDAC-UNIQUE     Arid3a_used	1700017B05RIK	8430427H17RIK	ABCA3	ANKRD44	ARL4A	BNC2	CLDN3
    XEN-SPECIFIC	Arid3a_used	1110036O03RIK	A130022J15RIK	B2M	B3GALNT1	CBX4	CITED1	CLU	CTSH	CYP26A1
    GATA-SPECIFIC	Arid3a_used	1200009I06RIK	5430407P10RIK	BAIAP2L1	BMP8B	CITED1	CLDN3	COBLL1	CORO1A	CRYAB	CTDSPL	DKKL1
    TS-SPECIFIC	Arid3a_used	5430407P10RIK	AFAP1L1	AHNAK	ANXA2	ANXA3	ANXA5	B2M	BIK	BMP8B	CAMK1D	CBX4	CLDN3	CSRP1	DKKL1	DSC2


Use ``enrichr`` command, or :func:`enrichr`
===============================================================

The only thing you need to prepare is a gene list file.

**Note**: Enrichr uses a list of Entrez gene symbols as input.


For ``enrichr`` , you could assign a list object

.. code:: python

    # assign a list object to enrichr
    l = ['SCARA3', 'LOC100044683', 'CMBL', 'CLIC6', 'IL13RA1', 'TACSTD2', 'DKKL1', 'CSF1',
         'SYNPO2L', 'TINAGL1', 'PTX3', 'BGN', 'HERC1', 'EFNA1', 'CIB2', 'PMP22', 'TMEM173']

    gseapy.enrichr(gene_list=l, description='pathway', gene_sets='KEGG_2016', outfile='test')




or a gene list file in txt format(one gene id per row)

.. code:: python

   gseapy.enrichr(gene_list='gene_list.txt', description='pathway', gene_sets='KEGG_2016', outfile='test')


Let's see what the txt file looks like.

.. code:: python

    with open('data/gene_list.txt') as genes:
        print(genes.read())

.. code:: python

    CTLA2B
    SCARA3
    LOC100044683
    CMBL
    CLIC6
    IL13RA1
    TACSTD2
    DKKL1
    CSF1
    CITED1
    SYNPO2L
    TINAGL1
    PTX3


Select the library you want to do enrichment analysis. To get a list of all available libraries, run

.. code:: python

   #s get_library_name(), it will print out all library names.
   import gseapy
   names = gseapy.get_library_name()
   print(names)


.. code:: python

   ['Genome_Browser_PWMs',
  'TRANSFAC_and_JASPAR_PWMs',
  'ChEA_2013',
  'Drug_Perturbations_from_GEO_2014',
  'ENCODE_TF_ChIP-seq_2014',
  'BioCarta_2013',
  'Reactome_2013',
  'WikiPathways_2013',
  'Disease_Signatures_from_GEO_up_2014',
  'KEGG_2013',
  'TF-LOF_Expression_from_GEO',
  'TargetScan_microRNA',
  'PPI_Hub_Proteins',
  'GO_Molecular_Function_2015',
  'GeneSigDB',
  'Chromosome_Location',
  'Human_Gene_Atlas',
  'Mouse_Gene_Atlas',
  'GO_Cellular_Component_2015',
  'GO_Biological_Process_2015',
  'Human_Phenotype_Ontology',
  'Epigenomics_Roadmap_HM_ChIP-seq',
  'KEA_2013',
  'NURSA_Human_Endogenous_Complexome',
  'CORUM',
  'SILAC_Phosphoproteomics',
  'MGI_Mammalian_Phenotype_Level_3',
  'MGI_Mammalian_Phenotype_Level_4',
  'Old_CMAP_up',
  'Old_CMAP_down',
  'OMIM_Disease',
  'OMIM_Expanded',
  'VirusMINT',
  'MSigDB_Computational',
  'MSigDB_Oncogenic_Signatures',
  'Disease_Signatures_from_GEO_down_2014',
  'Virus_Perturbations_from_GEO_up',
  'Virus_Perturbations_from_GEO_down',
  'Cancer_Cell_Line_Encyclopedia',
  'NCI-60_Cancer_Cell_Lines',
  'Tissue_Protein_Expression_from_ProteomicsDB',
  'Tissue_Protein_Expression_from_Human_Proteome_Map',
  'HMDB_Metabolites',
  'Pfam_InterPro_Domains',
  'GO_Biological_Process_2013',
  'GO_Cellular_Component_2013',
  'GO_Molecular_Function_2013',
  'Allen_Brain_Atlas_up',
  'ENCODE_TF_ChIP-seq_2015',
  'ENCODE_Histone_Modifications_2015',
  'Phosphatase_Substrates_from_DEPOD',
  'Allen_Brain_Atlas_down',
  'ENCODE_Histone_Modifications_2013',
  'Achilles_fitness_increase',
  'Achilles_fitness_decrease',
  'MGI_Mammalian_Phenotype_2013',
  'BioCarta_2015',
  'HumanCyc_2015',
  'KEGG_2015',
  'NCI-Nature_2015',
  'Panther_2015',
  'WikiPathways_2015',
  'Reactome_2015',
  'ESCAPE',
  'HomoloGene',
  'Disease_Perturbations_from_GEO_down',
  'Disease_Perturbations_from_GEO_up',
  'Drug_Perturbations_from_GEO_down',
  'Genes_Associated_with_NIH_Grants',
  'Drug_Perturbations_from_GEO_up',
  'KEA_2015',
  'Single_Gene_Perturbations_from_GEO_up',
  'Single_Gene_Perturbations_from_GEO_down',
  'ChEA_2015',
  'dbGaP',
  'LINCS_L1000_Chem_Pert_up',
  'LINCS_L1000_Chem_Pert_down',
  'GTEx_Tissue_Sample_Gene_Expression_Profiles_down',
  'GTEx_Tissue_Sample_Gene_Expression_Profiles_up',
  'Ligand_Perturbations_from_GEO_down',
  'Aging_Perturbations_from_GEO_down',
  'Aging_Perturbations_from_GEO_up',
  'Ligand_Perturbations_from_GEO_up',
  'MCF7_Perturbations_from_GEO_down',
  'MCF7_Perturbations_from_GEO_up',
  'Microbe_Perturbations_from_GEO_down',
  'Microbe_Perturbations_from_GEO_up',
  'LINCS_L1000_Ligand_Perturbations_down',
  'LINCS_L1000_Ligand_Perturbations_up',
  'LINCS_L1000_Kinase_Perturbations_down',
  'LINCS_L1000_Kinase_Perturbations_up',
  'Reactome_2016',
  'KEGG_2016',
  'WikiPathways_2016',
  'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',
  'Kinase_Perturbations_from_GEO_down',
  'Kinase_Perturbations_from_GEO_up',
  'BioCarta_2016',
  'Humancyc_2016',
  'NCI-Nature_2016',
  'Panther_2016']


For more details, please track the official links: http://amp.pharm.mssm.edu/Enrichr/


Use ``replot`` Command, or :func:`replot`
===============================================================

You may also want to use :func:`replot()` to reproduce ``GSEA`` desktop plots.

The only input of :func:`replot` is the directory of ``GSEA`` desktop output.

The input directory(e.g. gsea), must contained **edb** folder, gseapy need 4 data files
inside edb folder.The gsea document tree looks like this::

    gsea
    └─edb
        └─test.cls
        └─gene_sets.gmt
        └─gsea_data.rnk
        └─results.edb

After this, you can start to run gseapy.

.. code:: python

    import gseapy
    gseapy.replot(indir ='gsea', outdir = 'gseapy_out')


If you prefer to run in command line, it's more simple.

.. code:: bash

   gseapy replot -i gsea -o gseapy_out


| For advanced usage of library, see the :ref:`run`.
