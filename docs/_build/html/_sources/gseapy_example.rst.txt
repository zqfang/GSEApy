
GSEAPY Example
==============

Examples to walk through ``GSEApy``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Load essential packages
--------------------------

.. code:: ipython2

    %matplotlib inline
    import pandas as pd
    import gseapy as gp
    import matplotlib.pyplot as plt

\*\* Check gseapy version \*\*

.. code:: ipython2

    gp.__version__




.. parsed-literal::

    '0.8.7'



See all gseapy supported enrichr library names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Enrichr library could be used for ``gsea``, ``ssgsea``, and ``prerank``,
too

.. code:: ipython2

    names = gp.get_library_name()
    names[:10]




.. parsed-literal::

    [u'ARCHS4_Cell-lines',
     u'ARCHS4_IDG_Coexp',
     u'ARCHS4_Kinases_Coexp',
     u'ARCHS4_TFs_Coexp',
     u'ARCHS4_Tissues',
     u'Achilles_fitness_decrease',
     u'Achilles_fitness_increase',
     u'Aging_Perturbations_from_GEO_down',
     u'Aging_Perturbations_from_GEO_up',
     u'Allen_Brain_Atlas_down']



2. Enrichr Example
------------------

1) Assign enrichr with ``pd.Series``, ``pd.DataFrame``, or ``list`` object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython2

    gene_list = pd.read_table("./gene_list.txt",header=None)
    gene_list.head()




.. raw:: html

    <div>
    <style>
        .dataframe thead tr:only-child th {
            text-align: right;
        }
    
        .dataframe thead th {
            text-align: left;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>0</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>CTLA2B</td>
        </tr>
        <tr>
          <th>1</th>
          <td>SCARA3</td>
        </tr>
        <tr>
          <th>2</th>
          <td>LOC100044683</td>
        </tr>
        <tr>
          <th>3</th>
          <td>CMBL</td>
        </tr>
        <tr>
          <th>4</th>
          <td>CLIC6</td>
        </tr>
      </tbody>
    </table>
    </div>



.. code:: ipython2

    type(gene_list)




.. parsed-literal::

    pandas.core.frame.DataFrame



.. code:: ipython2

    # convert dataframe or series to list
    glist = gene_list.squeeze().tolist()
    print(glist[:10])


.. parsed-literal::

    ['CTLA2B', 'SCARA3', 'LOC100044683', 'CMBL', 'CLIC6', 'IL13RA1', 'TACSTD2', 'DKKL1', 'CSF1', 'CITED1']


.. code:: ipython2

    # run enrichr
    # if you are only intrested in dataframe that enrichr returned, please set no_plot=True
    
    # list, dataframe, series inputs are supported
    enr = gp.enrichr(gene_list="./gene_list.txt", 
                     # or gene_list='./gene_list.txt', or gene_list=glist
                     description='test_name', 
                     gene_sets='KEGG_2016', 
                     outdir='enrichr_kegg', 
                     cutoff=0.5 # test dataset, use lower value of range(0,1)
                    )


.. code:: ipython2

    enr.res2d.head()




.. raw:: html

    <div>
    <style>
        .dataframe thead tr:only-child th {
            text-align: right;
        }
    
        .dataframe thead th {
            text-align: left;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>Term</th>
          <th>Overlap</th>
          <th>P-value</th>
          <th>Adjusted P-value</th>
          <th>Old P-value</th>
          <th>Old Adjusted P-value</th>
          <th>Z-score</th>
          <th>Combined Score</th>
          <th>Genes</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>Rap1 signaling pathway_Homo sapiens_hsa04015</td>
          <td>19/211</td>
          <td>0.000148</td>
          <td>0.035223</td>
          <td>0.000436</td>
          <td>0.103734</td>
          <td>-1.961363</td>
          <td>17.295956</td>
          <td>PDGFRB;CSF1;FLT4;VEGFC;ARAP3;LPAR4;ADCY7;ADCY6...</td>
        </tr>
        <tr>
          <th>1</th>
          <td>Pathways in cancer_Homo sapiens_hsa05200</td>
          <td>27/397</td>
          <td>0.000729</td>
          <td>0.066282</td>
          <td>0.001816</td>
          <td>0.152127</td>
          <td>-2.083086</td>
          <td>15.046848</td>
          <td>RET;LEF1;TGFA;LPAR4;ADCY7;ETS1;ADCY6;GLI2;FGF4...</td>
        </tr>
        <tr>
          <th>2</th>
          <td>Ras signaling pathway_Homo sapiens_hsa04014</td>
          <td>18/227</td>
          <td>0.000999</td>
          <td>0.066282</td>
          <td>0.002351</td>
          <td>0.152127</td>
          <td>-1.956845</td>
          <td>13.519663</td>
          <td>PDGFRB;CSF1;FLT4;VEGFC;ETS1;GNG13;FGF4;PLD2;EF...</td>
        </tr>
        <tr>
          <th>3</th>
          <td>Dilated cardiomyopathy_Homo sapiens_hsa05414</td>
          <td>10/90</td>
          <td>0.001114</td>
          <td>0.066282</td>
          <td>0.002557</td>
          <td>0.152127</td>
          <td>-1.805957</td>
          <td>12.280169</td>
          <td>DES;SGCB;TPM2;TNNC1;LMNA;TPM1;ITGAV;ADCY7;ADCY...</td>
        </tr>
        <tr>
          <th>4</th>
          <td>HTLV-I infection_Homo sapiens_hsa05166</td>
          <td>19/258</td>
          <td>0.001747</td>
          <td>0.083151</td>
          <td>0.003877</td>
          <td>0.184562</td>
          <td>-1.843079</td>
          <td>11.703417</td>
          <td>PDGFRB;STAT5B;EGR1;JUN;CD40;FZD2;CRTC3;NFATC1;...</td>
        </tr>
      </tbody>
    </table>
    </div>



2) Command line usage
~~~~~~~~~~~~~~~~~~~~~

You may also want to use enrichr in command line

the option **-v** will print out the progress of your job

.. code:: ipython2

    !gseapy enrichr -i ./gene_list.txt \
                   --description BP2017 \
                   -g GO_Biological_Process_2017 \
                   -v -o enrichr_BP


.. parsed-literal::

    2017-11-24 13:11:55,413 Connecting to Enrichr Server to get latest library names
    2017-11-24 13:11:56,232 Analysis name: BP2017, Enrichr Library: GO_Biological_Process_2017
    2017-11-24 13:11:58,805 Submitted gene list:{'shortId': '350iz', 'userListId': 6127777}
    2017-11-24 13:12:04,922 Downloading file of enrichment results: Job Id:{'shortId': '350iz', 'userListId': 6127777}
    2017-11-24 13:12:08,329 Warning: No enrich terms using library GO_Biological_Process_2017 when cuttoff = 0.05
    2017-11-24 13:12:08,329 Done.
    


3. Prerank example
------------------

1) Assign prerank() with a pd.DataFrame, pd.Series , or a txt file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Do not include header !
| GSEApy will skip any comment lines startswith “#”.
| Only contains two columns, or one cloumn with gene_name indexed when
  assign a ``DataFrame`` to prerank

.. code:: ipython2

    rank = pd.read_table("./edb/gsea_data.gsea_data.rnk", header=None)
    rank.head()




.. raw:: html

    <div>
    <style>
        .dataframe thead tr:only-child th {
            text-align: right;
        }
    
        .dataframe thead th {
            text-align: left;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>0</th>
          <th>1</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>CTLA2B</td>
          <td>2.502482</td>
        </tr>
        <tr>
          <th>1</th>
          <td>SCARA3</td>
          <td>2.095578</td>
        </tr>
        <tr>
          <th>2</th>
          <td>LOC100044683</td>
          <td>1.116398</td>
        </tr>
        <tr>
          <th>3</th>
          <td>CMBL</td>
          <td>0.877640</td>
        </tr>
        <tr>
          <th>4</th>
          <td>CLIC6</td>
          <td>0.822181</td>
        </tr>
      </tbody>
    </table>
    </div>



.. code:: ipython2

    # run prerank
    # enrichr libraries are supported by prerank module. Just provide the name
    pre=[]
    for s, n in zip(['./genes.gmt', 'KEGG_2016'],['bp','kegg']):
        #use 4 process to acceralate the permutation speed
        pre_res = gp.prerank(rnk=rank, 
                             gene_sets=s, 
                             processes=4,
                             permutation_num=100, # reduce number to speed up test
                             outdir='prerank_report_'+n,format='png')
        pre.append(pre_res)

.. code:: ipython2

    #access results through res2d attribute
    pre[0].res2d.head()




.. raw:: html

    <div>
    <style>
        .dataframe thead tr:only-child th {
            text-align: right;
        }
    
        .dataframe thead th {
            text-align: left;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>es</th>
          <th>nes</th>
          <th>pval</th>
          <th>fdr</th>
          <th>gene_set_size</th>
          <th>matched_size</th>
          <th>genes</th>
        </tr>
        <tr>
          <th>Term</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>DvA_UpIN_A</th>
          <td>0.405880</td>
          <td>1.590404</td>
          <td>0.015873</td>
          <td>0.075188</td>
          <td>284</td>
          <td>19</td>
          <td>ABHD14B,VNN1,NELF,MARVELD2,LAMB3,TMPRSS2,TM6SF...</td>
        </tr>
        <tr>
          <th>DvA_UpIN_D</th>
          <td>0.166924</td>
          <td>0.626715</td>
          <td>0.842857</td>
          <td>0.872180</td>
          <td>236</td>
          <td>21</td>
          <td>PMP22,STBD1,DUSP14,RET,GPX8,CHRNB1,PRKD1,COL7A...</td>
        </tr>
      </tbody>
    </table>
    </div>



.. code:: ipython2

    pre_res = pre[0]
    prerank_results = pre_res.res2d
    prerank_results = prerank_results.reset_index()
    prerank_results.head(5).plot.barh(y='fdr',x='Term',fontsize=16)




.. parsed-literal::

    <matplotlib.axes._subplots.AxesSubplot at 0x10dbb69d0>




.. image:: output_20_1.png


2) Command line usage
~~~~~~~~~~~~~~~~~~~~~

You may also want to use prerank in command line

.. code:: ipython2

    # ! gseapy prerank -r temp.rnk -g temp.gmt -o prerank_report_temp

4. GSEA Example
---------------

1) Assign gsea() with a pandas DataFrame, .gct format file, or a text file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

and cls with a list object or just .cls format file

.. code:: ipython2

    phenoA, phenoB, class_vector =  gp.parser.gsea_cls_parser("./P53.cls")

.. code:: ipython2

    #class_vector used to indicate group attributes for each sample
    print(class_vector)


.. parsed-literal::

    ['MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'MUT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT']


.. code:: ipython2

    gene_exp = pd.read_table("./P53_resampling_data.txt")
    gene_exp.head()




.. raw:: html

    <div>
    <style>
        .dataframe thead tr:only-child th {
            text-align: right;
        }
    
        .dataframe thead th {
            text-align: left;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>NAME</th>
          <th>786-0</th>
          <th>BT-549</th>
          <th>CCRF-CEM</th>
          <th>COLO 205</th>
          <th>EKVX</th>
          <th>HCC-2998</th>
          <th>HCT-15</th>
          <th>HOP-62</th>
          <th>HOP-92</th>
          <th>...</th>
          <th>MCF7</th>
          <th>MOLT-4</th>
          <th>NCI-H460</th>
          <th>OVCAR-4</th>
          <th>SF-539</th>
          <th>SK-MEL-5</th>
          <th>SR</th>
          <th>UACC-257</th>
          <th>UACC-62</th>
          <th>UO-31</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>CTLA2B</td>
          <td>111.19</td>
          <td>86.22</td>
          <td>121.85</td>
          <td>75.19</td>
          <td>208.62</td>
          <td>130.59</td>
          <td>124.72</td>
          <td>324.09</td>
          <td>242.71</td>
          <td>...</td>
          <td>163.76</td>
          <td>59.50</td>
          <td>134.12</td>
          <td>152.09</td>
          <td>197.46</td>
          <td>137.79</td>
          <td>81.53</td>
          <td>123.37</td>
          <td>81.41</td>
          <td>180.78</td>
        </tr>
        <tr>
          <th>1</th>
          <td>SCARA3</td>
          <td>460.30</td>
          <td>558.34</td>
          <td>183.55</td>
          <td>37.29</td>
          <td>158.00</td>
          <td>43.61</td>
          <td>80.83</td>
          <td>300.08</td>
          <td>1250.25</td>
          <td>...</td>
          <td>109.91</td>
          <td>120.42</td>
          <td>73.06</td>
          <td>115.03</td>
          <td>95.12</td>
          <td>37.56</td>
          <td>76.16</td>
          <td>41.10</td>
          <td>77.51</td>
          <td>519.17</td>
        </tr>
        <tr>
          <th>2</th>
          <td>LOC100044683</td>
          <td>97.25</td>
          <td>118.94</td>
          <td>81.17</td>
          <td>119.51</td>
          <td>119.88</td>
          <td>107.73</td>
          <td>165.57</td>
          <td>203.97</td>
          <td>135.43</td>
          <td>...</td>
          <td>222.84</td>
          <td>124.98</td>
          <td>114.75</td>
          <td>141.66</td>
          <td>170.19</td>
          <td>147.70</td>
          <td>157.48</td>
          <td>152.18</td>
          <td>98.89</td>
          <td>118.06</td>
        </tr>
        <tr>
          <th>3</th>
          <td>CMBL</td>
          <td>33.45</td>
          <td>55.10</td>
          <td>221.67</td>
          <td>50.30</td>
          <td>35.12</td>
          <td>75.70</td>
          <td>84.01</td>
          <td>44.12</td>
          <td>79.96</td>
          <td>...</td>
          <td>51.32</td>
          <td>117.11</td>
          <td>59.46</td>
          <td>78.46</td>
          <td>45.55</td>
          <td>49.07</td>
          <td>96.69</td>
          <td>33.09</td>
          <td>10.38</td>
          <td>52.89</td>
        </tr>
        <tr>
          <th>4</th>
          <td>CLIC6</td>
          <td>35.75</td>
          <td>41.26</td>
          <td>63.04</td>
          <td>219.86</td>
          <td>42.53</td>
          <td>54.19</td>
          <td>86.98</td>
          <td>71.20</td>
          <td>53.89</td>
          <td>...</td>
          <td>154.05</td>
          <td>31.62</td>
          <td>37.66</td>
          <td>32.64</td>
          <td>63.35</td>
          <td>27.95</td>
          <td>70.99</td>
          <td>36.25</td>
          <td>17.50</td>
          <td>49.41</td>
        </tr>
      </tbody>
    </table>
    <p>5 rows × 51 columns</p>
    </div>



.. code:: ipython2

    print("positively correlated: ", phenoA)


.. parsed-literal::

    ('positively correlated: ', 'MUT')


.. code:: ipython2

    print("negtively correlated: ", phenoB)


.. parsed-literal::

    ('negtively correlated: ', 'WT')


.. code:: ipython2

    # run gsea
    # enrichr libraries are supported by gsea module. Just provide the name
    
    gs_res = gp.gsea(data=gene_exp, # or data='./P53_resampling_data.txt'
                     gene_sets='KEGG_2016', # enrichr library names
                     cls=class_vector, # or cls= './P53.cls'
                     #set permutation_type to phenotype if samples >=15
                     permutation_type='phenotype', 
                     permutation_num=100, # reduce number to speed up test
                     outdir='gsea_reprot', 
                     method='signal_to_noise', 
                     format='png')

.. code:: ipython2

    #access the dataframe results throught res2d attribute
    gs_res.res2d.head()




.. raw:: html

    <div>
    <style>
        .dataframe thead tr:only-child th {
            text-align: right;
        }
    
        .dataframe thead th {
            text-align: left;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>es</th>
          <th>nes</th>
          <th>pval</th>
          <th>fdr</th>
          <th>gene_set_size</th>
          <th>matched_size</th>
          <th>genes</th>
        </tr>
        <tr>
          <th>Term</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>MAPK signaling pathway_Homo sapiens_hsa04010</th>
          <td>-0.392928</td>
          <td>-1.270760</td>
          <td>0.166667</td>
          <td>0.474684</td>
          <td>255</td>
          <td>18</td>
          <td>GADD45B,RRAS,SOS2,FGF17,PPP3CC,TNFRSF1A,PDGFRB...</td>
        </tr>
        <tr>
          <th>HTLV-I infection_Homo sapiens_hsa05166</th>
          <td>-0.249752</td>
          <td>-0.790485</td>
          <td>0.818182</td>
          <td>0.743671</td>
          <td>258</td>
          <td>19</td>
          <td>FZD2,ETS1,STAT5B,RRAS,LTBR,PPP3CC,TNFRSF1A,EGR...</td>
        </tr>
        <tr>
          <th>Rap1 signaling pathway_Homo sapiens_hsa04015</th>
          <td>-0.285975</td>
          <td>-0.914519</td>
          <td>0.609756</td>
          <td>0.873418</td>
          <td>211</td>
          <td>19</td>
          <td>RRAS,VEGFC,CSF1,FGF17,PDGFRB,FGF4,PDGFC,SIPA1L...</td>
        </tr>
        <tr>
          <th>PI3K-Akt signaling pathway_Homo sapiens_hsa04151</th>
          <td>0.182245</td>
          <td>0.590397</td>
          <td>0.978723</td>
          <td>0.968750</td>
          <td>341</td>
          <td>22</td>
          <td>GNG13,VEGFC,GNB4,CSF1,SOS2,FGF17,THBS4,PDGFRB,...</td>
        </tr>
        <tr>
          <th>Cytokine-cytokine receptor interaction_Homo sapiens_hsa04060</th>
          <td>0.229069</td>
          <td>0.670014</td>
          <td>0.884615</td>
          <td>1.000000</td>
          <td>265</td>
          <td>18</td>
          <td>IL10RB,VEGFC,CSF1,TNFSF12,LTBR,CXCL10,TNFRSF1A...</td>
        </tr>
      </tbody>
    </table>
    </div>



.. code:: ipython2

    gsea_results= gs_res.res2d
    with plt.style.context('ggplot'):
        gsea_results = gsea_results.reset_index()
        gsea_results.head(5).plot.barh(y='fdr',x='Term',fontsize=16)



.. image:: output_32_0.png


2) Show the gsea plots
~~~~~~~~~~~~~~~~~~~~~~

The **gsea** module will generate heatmap for genes in each gene sets in
the backgroud.

.. code:: ipython2

    from IPython.display import Image
    
    #erich plot
    Image("./gsea_reprot/MAPK signaling pathway_Homo sapiens_hsa04010.gsea.png",width=650, height=600)




.. image:: output_34_0.png
   :width: 650px
   :height: 600px



.. code:: ipython2

    #corresponding heatmap
    Image("./gsea_reprot/MAPK signaling pathway_Homo sapiens_hsa04010.heatmap.png")




.. image:: output_35_0.png



3) Command line usage
~~~~~~~~~~~~~~~~~~~~~

You may also want to use gsea in command line

.. code:: ipython2

    # !gseapy gsea -d ./P53_resampling_data.txt -g KEGG_2016 -c ./P53.cls -o gsea_reprot_2 -v -t phenotype

5. Single Sample GSEA example
-----------------------------

**Note: When you run ssGSEA, all genes names in your gene_sets file
should be found in your expression table**

1) Assign ssgsea() with a txt file, dataframe, or Seires(gene name as index).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython2

    # txt file input
    ss = gp.ssgsea(data="./testSet_rand1200.gct",
                   gene_sets="./randomSets.gmt", 
                   outdir='ssgsea_report', 
                   sample_norm_method='rank', # choose 'custom' for your own rank list 
                   permutation_num=100, # reduce number to speed up test
                   processes=4, format='png')

.. code:: ipython2

    # or assign a dataframe, or Series to ssgsea()
    ssdf = pd.read_table("./temp.txt",header=None)
    ssdf.head()




.. raw:: html

    <div>
    <style>
        .dataframe thead tr:only-child th {
            text-align: right;
        }
    
        .dataframe thead th {
            text-align: left;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>0</th>
          <th>1</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>ATXN1</td>
          <td>16.456753</td>
        </tr>
        <tr>
          <th>1</th>
          <td>UBQLN4</td>
          <td>13.989493</td>
        </tr>
        <tr>
          <th>2</th>
          <td>CALM1</td>
          <td>13.745533</td>
        </tr>
        <tr>
          <th>3</th>
          <td>DLG4</td>
          <td>12.796588</td>
        </tr>
        <tr>
          <th>4</th>
          <td>MRE11A</td>
          <td>12.787631</td>
        </tr>
      </tbody>
    </table>
    </div>



.. code:: ipython2

    # dataframe with one column is also supported by ssGSEA or Prerank
    # But you have to set gene_names as index
    ssdf2 = ssdf.set_index(0)
    ssdf2.head()




.. raw:: html

    <div>
    <style>
        .dataframe thead tr:only-child th {
            text-align: right;
        }
    
        .dataframe thead th {
            text-align: left;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>1</th>
        </tr>
        <tr>
          <th>0</th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>ATXN1</th>
          <td>16.456753</td>
        </tr>
        <tr>
          <th>UBQLN4</th>
          <td>13.989493</td>
        </tr>
        <tr>
          <th>CALM1</th>
          <td>13.745533</td>
        </tr>
        <tr>
          <th>DLG4</th>
          <td>12.796588</td>
        </tr>
        <tr>
          <th>MRE11A</th>
          <td>12.787631</td>
        </tr>
      </tbody>
    </table>
    </div>



.. code:: ipython2

    type(ssdf2)




.. parsed-literal::

    pandas.core.frame.DataFrame



.. code:: ipython2

    ssSeries = ssdf2.squeeze()
    type(ssSeries)




.. parsed-literal::

    pandas.core.series.Series



.. code:: ipython2

    #Series Example
    # supports dataframe and series
    for dat in [ssdf, ssdf2, ssSeries]:
        ss = gp.ssgsea(data=ssdf, 
                       gene_sets="./temp.gmt", 
                       outdir='ssgsea_report_series',
                       permutation_num=100, # reduce number to speed up test
                       processes=4, format='png')

.. code:: ipython2

    ss.res2d.head(5)




.. raw:: html

    <div>
    <style>
        .dataframe thead tr:only-child th {
            text-align: right;
        }
    
        .dataframe thead th {
            text-align: left;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>es</th>
          <th>nes</th>
          <th>pval</th>
          <th>fdr</th>
          <th>gene_set_size</th>
          <th>matched_size</th>
          <th>genes</th>
        </tr>
        <tr>
          <th>Term</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>hsa05205</th>
          <td>0.341007</td>
          <td>13.452099</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>203</td>
          <td>201</td>
          <td>CTNNB1,PRKACA,GRB2,EGFR,RAC1,PRKCA,KRAS,CD44,M...</td>
        </tr>
        <tr>
          <th>hsa05412</th>
          <td>0.290588</td>
          <td>9.710331</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>74</td>
          <td>74</td>
          <td>CTNNB1,ACTB,ITGB1,CACNG3,RYR2,CTNNA1,CACNA2D3,...</td>
        </tr>
        <tr>
          <th>hsa05410</th>
          <td>0.270626</td>
          <td>10.461864</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>83</td>
          <td>83</td>
          <td>ACTB,ITGB1,TPM3,CACNG3,RYR2,CACNA2D3,ITGAV,ITG...</td>
        </tr>
        <tr>
          <th>hsa05323</th>
          <td>0.166596</td>
          <td>5.752460</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>89</td>
          <td>89</td>
          <td>JUN,ITGB2,ATP6V1B2,ATP6V1E1,IL1A,TGFB1,TEK,ATP...</td>
        </tr>
        <tr>
          <th>hsa05322</th>
          <td>0.176818</td>
          <td>7.742502</td>
          <td>0.0</td>
          <td>0.0</td>
          <td>134</td>
          <td>134</td>
          <td>GRIN2B,H2AFX,ACTN1,HIST4H4,SNRPD1,C3,GRIN2A,SS...</td>
        </tr>
      </tbody>
    </table>
    </div>



2) ``ssgsea`` supports gene expression matix in gct format.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| if gene expression matrix is provided, ssgsea works like pandas
  apply(),
| which means it will compute NES,FDR … for every sample pairwise.
| finally, you can assces the reuslts through **resultsOnSamples**
  attribute.

Take previous gene_exp dataframe for example

.. code:: ipython2

    df = pd.read_table("./P53_resampling_data.txt")
    df.head()




.. raw:: html

    <div>
    <style>
        .dataframe thead tr:only-child th {
            text-align: right;
        }
    
        .dataframe thead th {
            text-align: left;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>NAME</th>
          <th>786-0</th>
          <th>BT-549</th>
          <th>CCRF-CEM</th>
          <th>COLO 205</th>
          <th>EKVX</th>
          <th>HCC-2998</th>
          <th>HCT-15</th>
          <th>HOP-62</th>
          <th>HOP-92</th>
          <th>...</th>
          <th>MCF7</th>
          <th>MOLT-4</th>
          <th>NCI-H460</th>
          <th>OVCAR-4</th>
          <th>SF-539</th>
          <th>SK-MEL-5</th>
          <th>SR</th>
          <th>UACC-257</th>
          <th>UACC-62</th>
          <th>UO-31</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>CTLA2B</td>
          <td>111.19</td>
          <td>86.22</td>
          <td>121.85</td>
          <td>75.19</td>
          <td>208.62</td>
          <td>130.59</td>
          <td>124.72</td>
          <td>324.09</td>
          <td>242.71</td>
          <td>...</td>
          <td>163.76</td>
          <td>59.50</td>
          <td>134.12</td>
          <td>152.09</td>
          <td>197.46</td>
          <td>137.79</td>
          <td>81.53</td>
          <td>123.37</td>
          <td>81.41</td>
          <td>180.78</td>
        </tr>
        <tr>
          <th>1</th>
          <td>SCARA3</td>
          <td>460.30</td>
          <td>558.34</td>
          <td>183.55</td>
          <td>37.29</td>
          <td>158.00</td>
          <td>43.61</td>
          <td>80.83</td>
          <td>300.08</td>
          <td>1250.25</td>
          <td>...</td>
          <td>109.91</td>
          <td>120.42</td>
          <td>73.06</td>
          <td>115.03</td>
          <td>95.12</td>
          <td>37.56</td>
          <td>76.16</td>
          <td>41.10</td>
          <td>77.51</td>
          <td>519.17</td>
        </tr>
        <tr>
          <th>2</th>
          <td>LOC100044683</td>
          <td>97.25</td>
          <td>118.94</td>
          <td>81.17</td>
          <td>119.51</td>
          <td>119.88</td>
          <td>107.73</td>
          <td>165.57</td>
          <td>203.97</td>
          <td>135.43</td>
          <td>...</td>
          <td>222.84</td>
          <td>124.98</td>
          <td>114.75</td>
          <td>141.66</td>
          <td>170.19</td>
          <td>147.70</td>
          <td>157.48</td>
          <td>152.18</td>
          <td>98.89</td>
          <td>118.06</td>
        </tr>
        <tr>
          <th>3</th>
          <td>CMBL</td>
          <td>33.45</td>
          <td>55.10</td>
          <td>221.67</td>
          <td>50.30</td>
          <td>35.12</td>
          <td>75.70</td>
          <td>84.01</td>
          <td>44.12</td>
          <td>79.96</td>
          <td>...</td>
          <td>51.32</td>
          <td>117.11</td>
          <td>59.46</td>
          <td>78.46</td>
          <td>45.55</td>
          <td>49.07</td>
          <td>96.69</td>
          <td>33.09</td>
          <td>10.38</td>
          <td>52.89</td>
        </tr>
        <tr>
          <th>4</th>
          <td>CLIC6</td>
          <td>35.75</td>
          <td>41.26</td>
          <td>63.04</td>
          <td>219.86</td>
          <td>42.53</td>
          <td>54.19</td>
          <td>86.98</td>
          <td>71.20</td>
          <td>53.89</td>
          <td>...</td>
          <td>154.05</td>
          <td>31.62</td>
          <td>37.66</td>
          <td>32.64</td>
          <td>63.35</td>
          <td>27.95</td>
          <td>70.99</td>
          <td>36.25</td>
          <td>17.50</td>
          <td>49.41</td>
        </tr>
      </tbody>
    </table>
    <p>5 rows × 51 columns</p>
    </div>



.. code:: ipython2

    # dataframe support for multisamples 
    ss = gp.ssgsea(data=df, 
                   gene_sets="edb/gene_sets.gmt", 
                   outdir='ssgsea_df_test', 
                   permutation_num=100, # reduce number to speed up test
                   processes=4, format='png')

| Results for all samples are saves to a dataframe,
| you can assces the reuslts through resultsOnSamples attribute.

.. code:: ipython2

    # es results for all samples are saves to dict.
    # convert to dataframe 
    ss2 = pd.DataFrame(ss.resultsOnSamples)
    ss2.head()




.. raw:: html

    <div>
    <style>
        .dataframe thead tr:only-child th {
            text-align: right;
        }
    
        .dataframe thead th {
            text-align: left;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>786-0</th>
          <th>A498</th>
          <th>A549/ATCC</th>
          <th>ACHN</th>
          <th>BT-549</th>
          <th>CAKI-1</th>
          <th>CCRF-CEM</th>
          <th>COLO 205</th>
          <th>EKVX</th>
          <th>HCC-2998</th>
          <th>...</th>
          <th>SN12C</th>
          <th>SNB-19</th>
          <th>SNB-75</th>
          <th>SR</th>
          <th>SW-620</th>
          <th>T-47D</th>
          <th>U251</th>
          <th>UACC-257</th>
          <th>UACC-62</th>
          <th>UO-31</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>CDX2 OE-SPECIFIC</th>
          <td>0.029866</td>
          <td>0.009337</td>
          <td>0.009991</td>
          <td>-0.006401</td>
          <td>0.021200</td>
          <td>0.004605</td>
          <td>0.040562</td>
          <td>-0.000645</td>
          <td>-0.038442</td>
          <td>0.034446</td>
          <td>...</td>
          <td>0.004149</td>
          <td>0.070998</td>
          <td>0.013445</td>
          <td>0.024796</td>
          <td>0.040438</td>
          <td>0.020149</td>
          <td>0.046420</td>
          <td>0.068213</td>
          <td>0.074007</td>
          <td>0.004254</td>
        </tr>
        <tr>
          <th>ES-SPECIFIC</th>
          <td>-0.034765</td>
          <td>-0.064445</td>
          <td>-0.093547</td>
          <td>-0.079831</td>
          <td>-0.087688</td>
          <td>-0.029811</td>
          <td>-0.134692</td>
          <td>-0.072833</td>
          <td>-0.094930</td>
          <td>-0.041915</td>
          <td>...</td>
          <td>-0.105478</td>
          <td>-0.070089</td>
          <td>-0.052564</td>
          <td>-0.147355</td>
          <td>-0.084841</td>
          <td>-0.109212</td>
          <td>-0.068695</td>
          <td>-0.063946</td>
          <td>-0.109282</td>
          <td>-0.031117</td>
        </tr>
        <tr>
          <th>GATA3 OE-SPECIFIC</th>
          <td>0.006017</td>
          <td>0.020498</td>
          <td>-0.002182</td>
          <td>-0.031979</td>
          <td>0.019852</td>
          <td>-0.013829</td>
          <td>-0.000989</td>
          <td>0.031977</td>
          <td>0.016550</td>
          <td>0.040318</td>
          <td>...</td>
          <td>0.006542</td>
          <td>0.007366</td>
          <td>-0.011697</td>
          <td>0.000496</td>
          <td>0.017734</td>
          <td>0.018743</td>
          <td>-0.013939</td>
          <td>0.004336</td>
          <td>0.010753</td>
          <td>-0.029823</td>
        </tr>
        <tr>
          <th>HDAC1 UNIQUE TARGETS</th>
          <td>-0.031270</td>
          <td>-0.003659</td>
          <td>-0.012545</td>
          <td>-0.042056</td>
          <td>-0.025924</td>
          <td>-0.037497</td>
          <td>0.010344</td>
          <td>-0.030242</td>
          <td>-0.023006</td>
          <td>-0.030458</td>
          <td>...</td>
          <td>-0.054006</td>
          <td>-0.042326</td>
          <td>-0.039184</td>
          <td>-0.000971</td>
          <td>-0.052314</td>
          <td>-0.022150</td>
          <td>-0.013882</td>
          <td>0.028718</td>
          <td>-0.010532</td>
          <td>-0.036023</td>
        </tr>
        <tr>
          <th>OCT4 KD-SPECIFIC</th>
          <td>-0.032486</td>
          <td>-0.041851</td>
          <td>-0.032199</td>
          <td>-0.041748</td>
          <td>-0.007263</td>
          <td>-0.033633</td>
          <td>0.009335</td>
          <td>-0.008391</td>
          <td>-0.027491</td>
          <td>-0.006586</td>
          <td>...</td>
          <td>-0.039470</td>
          <td>-0.016508</td>
          <td>-0.035164</td>
          <td>0.025732</td>
          <td>-0.008062</td>
          <td>-0.011654</td>
          <td>-0.022169</td>
          <td>-0.009253</td>
          <td>-0.010288</td>
          <td>-0.047723</td>
        </tr>
      </tbody>
    </table>
    <p>5 rows × 50 columns</p>
    </div>



3) command line usage of single sample gsea
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython2

    # !gseapy ssgsea -d ./testSet_rand1200.gct -g temp.gmt -o ssgsea_report2  -p 4
