import gseapy as gp

gene_list = ['IGKV4-1', 'CD55', 'IGKC', 'PPFIBP1',
             'ABHD4', 'PCSK6', 'PGD', 'ARHGDIB', 'ITGB2', 'CARD6']
enr = gp.enrichr(gene_list=gene_list,
                 gene_sets=['KEGG_2016', 'KEGG_2013'],
                 organism='Human',  # don't forget to set organism to the one you desired! e.g. Yeast
                 description='test_name',
                 outdir='test/enrichr_kegg',
                 no_plot=True,
                 cutoff=0.5  # test dataset, use lower value from range(0,1)
                 )
