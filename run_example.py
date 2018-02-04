import gseapy
import sys, logging

rnk = "GSE4773_DEG_Expt1_Control_vs_Group1_gene.rnk"

# Run GSEA preranked - Kegg
prerank_results = gseapy.prerank(
	rnk='data/rnk_lists/' + rnk,
	gene_sets='data/gene_sets/c2.cp.kegg.v6.0.symbols.gmt', 
	outdir='out/Prerank_KEGG/' + rnk[:-4], 
	permutation_num=100, 
	no_plot=True,
	processes=4
)

# Run GSEA-InContext - Kegg
gseapen_results = gseapy.incontext(
	rnk='data/rnk_lists/' + rnk,
	gene_sets='data/gene_sets/c2.cp.kegg.v6.0.symbols.gmt', 
	background_rnks='data/bg_rnk_lists/all_442_lists_permuted_x100.csv', 
	outdir='out/InContext_KEGG/' + rnk[:-4], 
	permutation_num=100,
	no_plot=True,
	processes=4
)

# Run GSEA preranked - Hallmarks
prerank_results = gseapy.prerank(
	rnk='data/rnk_lists/' + rnk,
	gene_sets='data/gene_sets/hallmarks.gmt', 
	outdir='out/Prerank_Hallmarks/' + rnk[:-4], 
	permutation_num=100, 
	no_plot=True,
	processes=4
)

# Run GSEA-InContext - Hallmarks
gseapen_results = gseapy.incontext(
	rnk='data/rnk_lists/' + rnk,
	gene_sets='data/gene_sets/hallmarks.gmt', 
	background_rnks='data/bg_rnk_lists/all_442_lists_permuted_x100.csv', 
	outdir='out/InContext_Hallmarks/' + rnk[:-4], 
	permutation_num=100,
	no_plot=True,
	processes=4
)

print('Done!')