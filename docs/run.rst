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
    gseapy.enrichr(gene_list='gene_list.txt', description='pathway_analysis, gene_set='KEGG_2016', outfile='test')


* The :func:`replot` function to reproduce GSEA desktop results

.. code:: python

   # An example to reproduce figures using replot module.
   gseapy.replot('gsea','gseapy_out')



* The :func:`call` function to computing es,nes,pval,fdr,and generate plots *de novo*.

.. code:: python

   # An example to calculate es, nes, pval,fdrs, and produce figures using gseapy.
   gseapy.call(data=gsea_dat.txt, gene_sets=gene_sets.gmt, cls=test.cls, outdir='gseapy_out', 
             min_size=15, max_size=1000, permutation_n = 1000, weighted_score_type=1,
             permutation_type = 'gene_set', method='log2_ratio_of_classes',ascending=False, 
             figsize=(6.5,6), format='png')


* The :func:`prerank` function to computing es,nes,pval,fdr,and generate plots using a pre-ranked gene list.

.. code:: python

   # An example to calculate es, nes, pval,fdrs, and produce figures using gseapy.
   gseapy.prerank(rnk=gsea_data.rnk, gene_sets=gene_sets.gmt, outdir='gseapy_out', min_size = 15,
                  max_size = 1000, permutation_n=1000, weighted_score_type=1, ascending = False, 
                  figsize = (6.5,6), format = 'png')

			 
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


2. The ``replot`` Command
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

   $ gseapy replot -h

   usage: gseapy replot [-h] -i [-o] [-f] [--figsize width height] [-w]

    optional arguments:
      -h, --help            show this help message and exit

    Positional arguments:
      -i , --indir          The GSEA desktop results directroy that you want to
                            reproduce the figure
      -o , --outdir         The gseapy output directory. Default: the current
                            working directory
      -f , --format         Format of output figures, choose from {'pdf', 'png',
                            'jpeg', 'eps'}. Default: 'pdf'.
      --figsize width height
                            The figsize keyword argument need two parameter to
                            define. Default: [6.5, 6]
      -w , --weight         Weighted_score of rank_metric. Please use the same 
                            value in GSEA. Choose from (0, 1, 1.5, 2),default: 1



3. The ``call`` Command
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

   $ gseapy call -h

    usage: gseapy call [-h] -i DATA -c CLS -g GMT [-p] [-o] [-f]
                        [--figsize width height] [--min-size] [--max-size] [-n]
                        [-w] [-m] [-a] [-s]
    
    optional arguments:
      -h, --help            show this help message and exit
        
    Input files arguments:
      -d DATA, --data DATA
                            Input gene expression Affymetrix dataset file in txt
                            format.Same with GSEA.
      -c CLS, --cls CLS     Input class vector (phenotype) file in CLS format.
                            Same with GSEA.
      -g GMT, --gmt GMT     Gene set database in GMT format. Same with GSEA.
      -p , --permu-type     Permutation type. Same with GSEA, choose from
                            {'gene_set', 'phenotype'}
    
    Output arguments:
      -o , --outdir         The GSEAPY output directory. Default: the current
                            working directory
      -f , --format         Format of output figures, choose from {'pdf', 'png',
                            'jpeg', 'eps'}. Default: 'pdf'.
      --figsize width height
                            The figsize keyword argument need two parameter to
                            define. Default: [6.5,6]
    
    GSEA advanced arguments:
      --min-size            Min size of input genes presented in Gene Sets.
                            Default: 15
      --max-size            Max size of input genes presented in Gene Sets.
                            Default: 1000
      -n , --permu-num      Number of random permutations. For calculating
                            esnulls. Default: 1000
      -w , --weight         Weighted_score of rank_metrics.For weighting input
                            genes. Choose from {0, 1, 1.5, 2},default: 1
      -m , --method         Methods to calculate correlations of ranking metrics.
                            Choose from {'signal_to_noise', 't_test',
                            'ratio_of_classes',
                            'diff_of_classes','log2_ratio_of_classes'}. Default:
                            'log2_ratio_of_classes'
      -a , --ascending      Rank metric sorting order. If the -a flag was chosen,
                            then ascending equals to True. Default: False.
      -t , --top-graph      Plot graphs for top sets of each phenotype. Default:
                            20
      -s , --seed           Number of random seed. Default: None

4. The ``prerank`` Command
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash
   
  $ gseapy prerank -h
    
    usage: gseapy prerank       [-h] -r RNK -g GMT [-l pos neg] [-o] [-f]
                                [--figsize width height] [--min-size]
                                [--max-size] [-n] [-w] [-a] [-s]

    optional arguments:
      -h, --help            show this help message and exit
    
    Input files arguments:
      -r RNK, --rnk RNK     ranking dataset file in .rnk format.Same with GSEA.
      -g GMT, --gmt GMT     Gene set database in GMT format. Same with GSEA.
      -l pos neg, --label pos neg
                            The phenotype label argument need two parameter to
                            define. Default: ['Postive','Negative']
    
    Output arguments:
      -o , --outdir         The GSEAPY output directory. Default: the current
                            working directory
      -f , --format         Format of output figures, choose from {'pdf', 'png',
                            'jpeg', 'eps'}. Default: 'pdf'.
      --figsize width height
                            The figsize keyword argument need two parameter to
                            define. Default: [6.5,6]
    
    GSEA advanced arguments:
      --min-size            Min size of input genes presented in Gene Sets.
                            Default: 15
      --max-size            Max size of input genes presented in Gene Sets.
                            Default: 1000
      -n , --permu-num      Number of random permutations. For calculating
                            esnulls. Default: 1000
      -w , --weight         Weighted_score of rank_metrics.For weighting input
                            genes. Choose from {0, 1, 1.5, 2},default: 1
      -a , --ascending      Rank metric sorting order. If the -a flag was chosen,
                            then ascending equals to True. Default: False.
      -t , --top-graph      Plot graphs for top sets of each phenotype. Default:
                            20
      -s , --seed           Number of random seed. Default: None


4. The ``prerank`` Command
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash
   
   $ gseapy enrichr -h

   usage: gseapy enrichr [-h] -i  [-d] -g  (-o OFILE | --o-prefix OFILE)

    optional arguments:
      -h, --help            show this help message and exit

    Positional arguments:
      -i , --input-list     Enrichr uses a list of Entrez gene symbols as input.
      -d , --description    It is recommended to enter a description for your list
                            so that multiple lists can be differentiated from each
                            other if you choose to save or share your list
      -g , --gene-sets      Enrichr uses a list of Entrez gene symbols as input.
      -o OFILE, --ofile OFILE
                            Output file name. Mutually exclusive with --o-prefix.
      --o-prefix OFILE      Output file prefix. Mutually exclusive with
                            -o/--ofile.  
							


							
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

.. automodule:: gseapy.gsea_plot
   :members:   

   