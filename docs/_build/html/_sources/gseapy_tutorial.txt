.. _example:

======================================
A Protocol to walk through GSEAPY
======================================

As a biological reseacher, I like protocols, so as other reseachers, too.

Here is an short tutorial to walk you through gseapy.

In order to run gseapy successfully, install gseapy use pip.

.. code:: bash

    pip install gseapy



Use ``call`` command, or :func:`call`
================================================

Follow the steps blow.

One thing you should know is that the gseapy input files are totally the same as
``GSEA`` desktop requried. You can use these files below to run ``GSEA`` desktop, too.


1. Prepare an tabular text file of gene expression like this:
------------------------------------------------------------------

**RNA-seq,ChIP-seq, Microarry data** are all supported.

Here is to see what the structure of expression table looks like, you don't have to run
commands below:

.. code:: python

    import pandas as pd
    data = pd.read_table('./test/gsea_data.txt')
    data.head()




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


.. parsed-literal::

    6 2 1
    # C1OE Vector
    C1OE C1OE C1OE Vector Vector Vector
    
    
| The first line specify the total samples and phenotype numbers. Leave number 1 alway be 1.
| The second line specify the phenotype class(name).
| The third line specify column attributes in setp 1.     





3. Gene_sets file in gmt format. 
-----------------------------------------------------

All you need to do is to download gene set database file from ``GSEA`` website.

If you would like to use you own gene_sets.gmts files, build such a file use excel,
and then rename to gene_sets.gmt.

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
    
    

4. Run gseapy inside python
-------------------------------------------------------

At least 3 files are required to run gseapy.

.. code:: python

    import gseapy
    gseapy.call(data=gsea_data.txt, cls=gsea.cls, gmt=gene_sets.gmt, outdir='gseapy_out')


5. Command line 
---------------------------------------------------------

.. code:: bash

    gseapy call -d gsea_data.txt -c test.cls -g gene_sets.gmt -o gseapy_out


Use ``prerank``, or :func:`prerank`
===============================================================

If you would like to use a pre-ranked gene list to run GSEAPY, ``prerank`` module expects
a pre-ranked gene list dataset with correlation values, which in .rnk format,
and gene_sets file in gmt format.  ``prerank`` module the same API to `GSEA` pre-rank tools.

After this, you can start to run gseapy.

.. code:: bash
 
    gseapy prerank -r gsea_data.rnk -g gene_sets.gmt -o test

Or run inside python

..code:: python

    import gseapy
    # using prerank tool
    gseapy.prerank(rnk=gsea_data.rnk, gene_sets=gene_sets.gmt, outdir='test')







Use ``replot``, or :func:`replot`
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


| For advanced usage of library,see the :ref:`run`. 