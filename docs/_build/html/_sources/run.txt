.. _run:


How to Use GSEAPY
=====================================

.. module:: gseapy


For command line usage:
--------------------------------

The ``replot`` module reproduce GSEA desktop version results. The only input for GSEAPY is the location to GSEA results.

The ``call`` module produce GSEAPY results. The input requries a txt file( FPKM, Expected Counts, TPM, et.al), a cls file,
and gene_sets file in gmt format. The input files' formats are identical to ``GSEA`` desktop version. 
See `GSEA  <http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Main_Page>`_ documentation for more information.



.. code:: bash
  
  $ gseapy replot -i path/to/GSEA_resutls_folder -o gesapy_out

  # An example to reproduce figures using replot module.
  $ gseapy replot -i ./Gsea.reports -o test
  
  
  # an example to compute using gseapy call module
  $ gseapy call -d exptable.txt -c test.cls -g test.gmt



Run gseapy inside python:
-----------------------------------

.. code:: python

   import gseapy

This will import the following:

* The :func:`replot` function to reproduce GSEA desktop results

.. code:: python

   # An example to reproduce figures using replot module.
   gseapy.replot('./Gsea.reports','test')



* The :func:`run` function to computing es,nes,pval,fdr,and generate plots *de novo*.

.. code:: python

   # An example to calculate es, nes, pval,fdrs, and produce figures using gseapy.
   gseapy.run(data = expression.txt, gene_sets = gene_sets.gmt,cls = test.cls, 
             min_size = 15, max_size = 1000, permutation_n = 1000, weighted_score_type = 1,
             permutation_type = 'gene_set', method = 'log2_ratio_of_classes',ascending = False, 
             outdir='test',figsize = (6.5,6),format = 'png')



			 
To See help information of GSEAPY
--------------------------------------

1. gseapy subcommands
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code:: bash
   
   $ gseapy --help 
   
    usage: gseapy [-h] [--version] {call,replot} ...

    gseapy -- Gene Set Enrichment Analysis in Python

    positional arguments:
      {call,replot}
        call       Main GSEAPY Function: run GSEAPY instead of GSEA.
        replot     Reproduce GSEA desktop figures.

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
                        [-w] [-m] [-a]
    
    optional arguments:
      -h, --help            show this help message and exit
        
    Input files arguments:
      -d DATA, --datab DATA
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
      -a, --ascending       Rank metric sorting order. If the -a flag was chosen,
                            then ascending equals to True. Default: False.



							
							
Module APIs 
---------------------


.. autofunction:: replot()



.. autofunction:: run()






Algorithms
-------------------------


.. automodule:: gseapy.algorithm
   :members:


Parsers
--------------------------

.. automodule:: gseapy.parser
   :members:

