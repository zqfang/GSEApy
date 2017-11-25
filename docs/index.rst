.. GSEAPY documentation master file, created by
   sphinx-quickstart on Mon Feb  1 13:52:57 2016.
   You can adapt this file completely to your linking, but it should at least
   contain the root `toctree` directive.

Welcome to GSEAPY's documentation!
=====================================================

GSEAPY: Gene Set Enrichment Analysis in Python.
------------------------------------------------

.. image:: https://badge.fury.io/py/gseapy.svg
    :target: https://badge.fury.io/py/gseapy

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square
    :target: http://bioconda.github.io

.. image:: https://travis-ci.org/BioNinja/GSEApy.svg?branch=master
    :target: https://travis-ci.org/BioNinja/GSEApy

.. image:: http://readthedocs.org/projects/gseapy/badge/?version=latest
    :target: http://gseapy.readthedocs.org/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target:  https://img.shields.io/badge/license-MIT-blue.svg
.. image:: https://img.shields.io/badge/python-3.5-blue.svg
    :target:   https://img.shields.io/badge/python-3.5-blue.svg
.. image:: https://img.shields.io/badge/python-2.7-blue.svg
    :target:  https://img.shields.io/badge/python-2.7-blue.svg



**Release notes** : https://github.com/BioNinja/gseapy/releases 

GSEAPY is a python wrapper for **GSEA** and **Enrichr**. 
--------------------------------------------------------------------------------------------

GSEAPY has five subcommands: ``gsea``, ``prerank``, ``ssgsea``, ``replot`` ``enrichr``.

1. The ``gsea`` module produce **GSEA** results.    
The input requries a txt file(FPKM, Expected Counts, TPM, et.al), a cls file, and gene_sets file in gmt format. 

2. The ``prerank`` module produce **Prerank tool** results.  
The input expects a pre-ranked gene list dataset with correlation values, which in .rnk format, and gene_sets file in gmt format.  ``prerank`` module is an API to `GSEA` pre-rank tools.

3. The ``ssgsea`` module perform **single sample GSEA(ssGSEA)** analysis.  
The input expects a gene list with expression values(same with ``.rnk`` file, and gene_sets file in gmt format. ssGSEA enrichment score for the gene set as described by `D. Barbie et al 2009 <http://www.nature.com/nature/journal/v462/n7269/abs/nature08460.html>`_.

4. The ``replot`` module reproduce GSEA desktop version results.  
The only input for GSEAPY is the location to GSEA Desktop output results.

5. The ``enrichr`` module enable you perform gene set enrichment analysis using ``Enrichr`` API.
Enrichr is open source and freely available online at: http://amp.pharm.mssm.edu/Enrichr . It runs very fast and generates results in txt format.

GSEAPY could be used for **RNA-seq, ChIP-seq, Microarry** data. It's used for convenient GO enrichments and produce **publishable quality figures** in python. 


The full ``GSEA`` is far too extensive to describe here; see
`GSEA  <http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Main_Page>`_ documentation for more information. All files' formats for GSEApy are identical to ``GSEA`` desktop version. 


**If you use gseapy, you should cite the original ``GSEA`` and ``Enrichr`` paper.**

Why GSEAPY
-----------------------------------------------------

I would like to use Pandas to explore my data, but I did not find a  convenient tool to
do gene set enrichment analysis in python. So, here is my reason: 

* **Running inside python interactive console without switch to R!!!**
* User friendly for both wet and dry lab usrers.
* Produce and reproduce pubilishable figures.
* Perform batch jobs easy(using for loops).
* Easy to use in bash shell or your  data analysis workflow, e.g. snakemake.  



.. toctree::
    :maxdepth: 3
   
    introduction.rst
    run.rst
    gseapy_tutorial
    gseapy_example.ipynb
	



   
   


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


