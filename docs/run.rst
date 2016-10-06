.. _run:

How to Use GSEAPY
=====================================

.. module:: gseapy


For command line usage:
--------------------------------

The ``enrichr`` module will call enrichr-api. You need a gene_list file for input. That's all.

.. code:: bash

  # An example to use enrichr api
  $ gseapy enrichr -i gene_list.txt -g KEGG_2016 -d pathway_enrichment -o test



The ``replot`` module reproduce GSEA desktop version results. The only input for replot module is the location to GSEA results.

.. code:: bash
  
    # An example to reproduce figures using replot module.
    $ gseapy replot -i ./gsea -o gseapy_out



The ``call`` module produce GSEAPY results. 

The input requries a txt file( FPKM, Expected Counts, TPM, et.al), a cls file,
and gene_sets file in gmt format.


.. code:: bash
    
    # an example to compute using gseapy call module
    $ gseapy call -d gsea_data.txt -c test.cls -g gene_sets.gmt



The ``prerank`` module input expects a pre-ranked gene list dataset with correlation values, which in .rnk format,
and gene_sets file in gmt format.  ``prerank`` module is an API to `GSEA` pre-rank tools.


.. code:: bash

    $ gseapy prerank -r gsea_data.rnk -g gene_sets.gmt -o test


The input files' formats are identical to ``GSEA`` desktop version. 
See :ref:`example` for details, or, see `GSEA  <http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Main_Page>`_ documentation for more information.



Run gseapy inside python:
-----------------------------------

.. code:: python

   import gseapy

This will import the following:

* The :func:`enrichr` function to perform gene set enrichment analysis by calling enrichr online api.

.. code:: python

    # call enrichr
    gseapy.enrichr(gene_list='gene_list.txt', description='pathway_analysis', gene_set='KEGG_2016', outdir='test')


* The :func:`replot` function to reproduce GSEA desktop results

.. code:: python

   # An example to reproduce figures using replot module.
   gseapy.replot('gsea','gseapy_out')



* The :func:`call` function to computing es,nes,pval,fdr,and generate plots *de novo*.

.. code:: python

   # An example to calculate es, nes, pval,fdrs, and produce figures using gseapy.
   gseapy.call(data='gsea_dat.txt', gene_sets='gene_sets.gmt', cls='test.cls', outdir='gseapy_out', 
             min_size=15, max_size=1000, permutation_n = 1000, weighted_score_type=1,
             permutation_type = 'gene_set', method='log2_ratio_of_classes', ascending=False, 
             figsize=(6.5,6), format='png')


* The :func:`prerank` function to computing es,nes,pval,fdr,and generate plots using a pre-ranked gene list.

.. code:: python

   # An example to calculate es, nes, pval,fdrs, and produce figures using gseapy.
   gseapy.prerank(rnk='gsea_data.rnk', gene_sets='ene_sets.gmt', outdir='gseapy_out', min_size=15,
                  max_size=1000, permutation_n=1000, weighted_score_type=1, ascending=False, 
                  figsize=(6.5,6), format='png')

			 
To See help information of GSEAPY
--------------------------------------

1. gseapy subcommands
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code:: bash
   
   $ gseapy --help 
   
    usage: gseapy [-h] [--version] {call,prerank,replot,enrichr} ...

    gseapy -- Gene Set Enrichment Analysis in Python

    positional arguments:
      {call,replot}
        call       Main GSEAPY Function: run GSEAPY instead of GSEA.
        prerank    Using pre-ranked tool to run GSEAPY.
        replot     Reproduce GSEA desktop figures.
        enrichr    Peform GSEA using enrichr API.
    optional arguments:
      -h, --help   show this help message and exit
      --version    show program's version number and exit



For command line options of each command, type: gseapy COMMAND -h


2. The subcommands help
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

   $ gseapy replot -h
   # or
   $ gseapy call -h
   # or
   $ gseapy prerank -h
   # or 
   $ gseapy enrichr -h


							


							
Module APIs 
---------------------


.. autofunction:: replot()


.. autofunction:: call()


.. autofunction:: prerank()


.. autofunction:: enrichr()


Algorithm 
-------------------------


.. automodule:: gseapy.algorithm
   :members:


Enrichr
--------------------------

.. automodule:: gseapy.enrichr
   :members:


Parser 
--------------------------

.. automodule:: gseapy.parser
   :members:


Graph 
--------------------------

.. automodule:: gseapy.plot
   :members:   

   