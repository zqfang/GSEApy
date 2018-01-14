
GSEAPY
========

GSEAPY: Gene Set Enrichment Analysis in Python.
------------------------------------------------

This is a fork of `GSEApy <https://github.com/BioNinja/GSEApy>`_ (`original documentation here <http://gseapy.rtfd.io/>`_). We have added a new tool ``GSEA_PEN`` which runs the GSEAPreranked algorithm but uses a background set of ranked lists to calculate an empirical null distribution for informing the permutation procedure.

See below for examples on running the GSEA-PEN algorithm.

For examples using the original GSEApy library, `visit this page <http://gseapy.readthedocs.io/en/master/gseapy_example.html>`_.


GSEAPY and GSEA-PEN
--------------------------------------------------------------------------------------------

GSEAPY can be used for RNA-seq, ChIP-seq, and microarry data. It's used for convenient GO enrichment and produces publishable quality figures in Python.


Original GSEAPY offers five sub-commands: ``gsea``, ``prerank``, ``ssgsea``, ``replot`` ``enrichr``.
We have added one additional sub-command: ``gsea_pen``.


:prerank: The ``prerank`` module produces **GSEAPreranked** results.  Required inputs: a pre-ranked gene list dataset with correlation values in .rnk format, and gene_sets file in gmt format.  The ``prerank`` module is an API to the `GSEA` pre-rank tools.

:gsea_pen:    The ``gsea_pen`` module **uses an informed permutation strategy** when running the GSEAPreranked algorithm.  Required inputs: a pre-ranked gene list dataset with correlation values in .rnk format, a gene_sets file in gmt format, and a path to a directory containing .rnk lists to be used as the background distribution.

:replot: The ``replot`` module reproduces GSEA desktop version results.  The only input for GSEApy is the location to ``GSEA`` Desktop output results.

:enrichr: The ``enrichr`` module enable you perform gene set enrichment analysis using ``Enrichr`` API. Enrichr is open source and freely available online at: http://amp.pharm.mssm.edu/Enrichr . It runs very fast.

Please use 'gseapy COMMAND -h' to see the detail description for each option of each module. For descriptions of the other commands, see the original `GSEApy <https://github.com/BioNinja/GSEApy>`_ repository.


The full ``GSEA`` is too extensive to describe here; see
`GSEA  <http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Main_Page>`_ documentation for more information. All files' formats for GSEApy are identical to ``GSEA`` desktop version.


About GSEA-PEN
-----------------------------------------------------

TODO add descriptions and pseudo-code


Installation
------------

| To install the gseapy package with the GSEA-PEN extended functionality:


.. code:: shell

   $ pip install git+git://github.com/CostelloLab/gseapy.git#egg=gseapy


Dependencies & Requirements
--------------
* Python 2.7 or 3.4+
* Numpy >= 1.13.0
* Pandas
* Matplotlib
* Beautifulsoup4
* Requests (for enrichr API)

You may also need to install lxml and html5lib to parse xml files.


Running GSEAPY and GSEA-PEN
-----------------

Before you start:
~~~~~~~~~~~~~~~~~~~~~~

Convert all gene symbol names to uppercase. The ranked lists input to ``prerank`` or ``gsea_pen`` can be supplied as file paths (.rnk) or a two-column Pandas DataFrame (columns ``gene_name`` and ``fold_change``). The background ranked lists input to ``gsea_pen`` can be supplied as a path to a directory containing the desired .rnk files, or as a Pandas DataFrame with each column containing the rank-ordered list of genes for a given background experiment.


For command line usage:
~~~~~~~~~~~~~~~~~~~~~~~

| Test installation of GSEAPY:

.. code:: bash


  # Test the prerank module
  $ gseapy prerank -r data/example_expt1.rnk -g gene_sets.gmt -o test

  # Test the gsea_pen module 
  $ gseapy prerank -r data/example_expt1.rnk -g gene_sets.gmt -o test

| Run GSEAPreranked or GSEA-PEN

.. code:: bash


  # Run GSEAPreranked with the prerank module
  $ gseapy prerank -r data/example_expt1.rnk -g gene_sets.gmt

  # Run GSEA-PEN with the gsea_pen module 
  $ gseapy gsea_pen -r data/example_expt1.rnk -g gene_sets.gmt -b data/example_background_rnks


Run GSEAPY inside python console:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Running GSEA, GSEAPreranked or GSEA-PEN in Python using file paths as input

.. code:: python

    import gseapy

    # Run GSEA
    gseapy.gsea(data='expression.txt', gene_sets='gene_sets.gmt', cls='test.cls', outdir='test')

    # Run GSEA Prerank
    gseapy.prerank(rnk='gsea_data.rnk', gene_sets='gene_sets.gmt', outdir='test')

    # Run GSEA-PEN
    gseapy.gsea_pen(rnk='gsea_data.rnk', gene_sets='gene_sets.gmt', backround_rnks = 'background_dir', outdir='test')

    # Plot figures using replot module.
    gseapy.replot(indir='./Gsea.reports', outdir='test')


| Running GSEA, GSEAPreranked or GSEA-PEN in Python using DataFrames as input

.. code:: python


    # TODO finish small reproducible examples
    
    gene_ranked_dataframe = pd.DataFrame()

    # Run GSEA (using local gene set file)
    expression_dataframe = pd.DataFrame()
    sample_names = ['Expt1','Expt1','Expt1','Ctrl','Ctrl','Ctrl'] # must contain exactly two groups
    gseapy.gsea(data=expression_dataframe, gene_sets='gene_set.gmt', cls=sample_names, outdir='test')

    # Run GSEAPreranked (using Enrichr to find gene sets)
    gene_ranked_dataframe = pd.DataFrame()
    gseapy.prerank(rnk=gene_ranked_dataframe, gene_sets='KEGG_2016', outdir='test')

    # Run GSEA-PEN (using directory path for background rnk lists)
    gene_ranked_dataframe = pd.DataFrame()
    gseapy.gsea_pen(rnk=gene_ranked_dataframe, gene_sets='KEGG_2016', background_rnks = 'background_dir' outdir='test')

    # Run GSEA-PEN (using DataFrame for background rnk lists)
    gene_ranked_dataframe = pd.DataFrame()
    bg_ranked_dataframe = pd.DataFrame()
    gseapy.gsea_pen(rnk=gene_ranked_dataframe, gene_sets='KEGG_2016', background_rnks = bg_ranked_dataframe, outdir='test')


GSEAPY / Enrichr supported gene set libaries :
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To see the full list of gseapy supported gene set libraries, please click here: `Library <http://amp.pharm.mssm.edu/Enrichr/#stats>`_

Or use ``get_library_name`` function inside python console.

.. code:: python

    # See full list of latest enrichr library names
    names = gseapy.get_library_name()

    # Preview top 20 entries
    print(names[:20])


   ['Genome_Browser_PWMs',
   'TRANSFAC_and_JASPAR_PWMs',
   'ChEA_2013',
   'Drug_Perturbations_from_GEO_2014',
   'ENCODE_TF_ChIP-seq_2014',
   'BioCarta_2013',
   'Reactome_2013',
   'WikiPathways_2013',
   'Disease_Signatures_from_GEO_up_2014',
   'KEGG_2016',
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
   'Human_Phenotype_Ontology',]


Bug Reports
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you would like to report any bugs when you running the ``gsea_pen`` module, please create an issue on GitHub here. For issues relating to other modules, you may wish to visit the `original GSEAPY repo <https://github.com/BioNinja/GSEApy>`_.
