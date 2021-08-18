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

.. image:: https://img.shields.io/conda/vn/bioconda/GSEApy.svg?style=plastic
    :target: http://bioconda.github.io

.. image:: https://github.com/zqfang/GSEApy/workflows/GSEApy/badge.svg?branch=master
    :target: https://github.com/zqfang/GSEApy/actions
    :alt: Action Status

.. image:: http://readthedocs.org/projects/gseapy/badge/?version=master
    :target: http://gseapy.readthedocs.io/en/master/?badge=master
    :alt: Documentation Status

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target:  https://img.shields.io/badge/license-MIT-blue.svg

.. image:: https://img.shields.io/pypi/pyversions/gseapy.svg
    :alt: PyPI - Python Version

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3748085.svg
   :target: https://doi.org/10.5281/zenodo.3748085


**Release notes** : https://github.com/zqfang/GSEApy/releases 

GSEApy is a python wrapper for **GSEA** and **Enrichr**. 
--------------------------------------------------------------------------------------------

GSEApy has six subcommands: ``gsea``, ``prerank``, ``ssgsea``, ``replot`` ``enrichr``, ``biomart``.

1. The ``gsea`` module produces **GSEA** results.    
The input requries a txt file(FPKM, Expected Counts, TPM, et.al), a cls file, and gene_sets file in gmt format. 

2. The ``prerank`` module produces **Prerank tool** results.  
The input expects a pre-ranked gene list dataset with correlation values, which in .rnk format, and gene_sets file in gmt format.  ``prerank`` module is an API to `GSEA` pre-rank tools.

3. The ``ssgsea`` module performs **single sample GSEA(ssGSEA)** analysis.  
The input expects a gene list with expression values(same with ``.rnk`` file, and gene_sets file in gmt format. ssGSEA enrichment score for the gene set as described by `D. Barbie et al 2009 <http://www.nature.com/nature/journal/v462/n7269/abs/nature08460.html>`_.

4. The ``replot`` module reproduces GSEA desktop version results.  
The only input for GSEAPY is the location to GSEA Desktop output results.

5. The ``enrichr`` module enables you to perform gene set enrichment analysis using ``Enrichr`` API.
Enrichr is open source and freely available online at: http://amp.pharm.mssm.edu/Enrichr . It runs very fast and generates results in txt format.

6. The ``biomart`` module helps you convert gene ids using BioMart API.


GSEApy could be used for **RNA-seq, ChIP-seq, Microarry** data. It's used for convenient GO enrichments and produce **publishable quality figures** in python. 


The full ``GSEA`` is far too extensive to describe here; see
`GSEA  <http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Main_Page>`_ documentation for more information. All files' formats for GSEApy are identical to ``GSEA`` desktop version. 


**If you use gseapy, you should cite the original ``GSEA`` and ``Enrichr`` paper.**

Why GSEAPY
-----------------------------------------------------

I would like to use Pandas to explore my data, but I did not find a convenient tool to
do gene set enrichment analysis in python. So, here are my reasons:

* **Ability to run inside python interactive console without having to switch to R!!!**
* User friendly for both wet and dry lab users.
* Produce or reproduce publishable figures.
* Perform batch jobs easy.
* Easy to use in bash shell or your data analysis workflow, e.g. snakemake.


.. toctree::
    :maxdepth: 3
   
    introduction.rst
    gseapy_tutorial
    gseapy_example.ipynb
    run.rst
    faq.rst
	



   
   


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


