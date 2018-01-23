import gseapy
import os
import numpy as np
from numpy.random import choice
from gseapy.algorithm import enrichment_score_pen, make_background_dist
from collections import defaultdict
import csv
import sys, logging

#bg_lst = [['A', 'B', 'C', 'D', 'E', 'F'], ['F', 'B', 'C', 'A', 'E', 'D'], ['A', 'F', 'B', 'D', 'E', 'C']]
#permed = make_background_dist(bg_lst, 3)
#print(permed)

"""
rnk='data/rnk_lists/GSE11791_DEG_Expt1_Control_vs_Group1_gene.rnk'
gene_sets='data/gene_sets/hallmarks.gmt'
background_rnks='data/bg_rnk_lists/all_443_background.txt'
outdir='out/all_permuted_data'
permutation_num=100
processes=4

# Defaults
pheno_pos='Pos'
pheno_neg='Neg'
min_size=15
max_size=500
weighted_score_type=1
ascending=False
figsize=(6.5,6)
mpl_format='pdf'
graph_num=20
no_plot=False
seed=None
verbose=False

pen = gseapy.GSEA_PEN(rnk, gene_sets, background_rnks, outdir, pheno_pos, pheno_neg,
	min_size, max_size, permutation_num, weighted_score_type,
	ascending, processes, figsize, mpl_format, graph_num, no_plot, seed, verbose)
dat2 = pen._load_ranking(rnk)
bg_lists = pen._load_background_ranking(background_rnks)
gmt = pen.load_gmt(gene_list=dat2.index.values, gmt=gene_sets)

permed = make_background_dist(bg_lists, permutation_num)

with open("out/all_443_lists_permuted_x100.csv", "w") as f:
	writer = csv.writer(f)
	writer.writerows(permed)

"""

rnks = os.listdir('data/rnk_lists')
rnks = [r for r in rnks if r.startswith('GSE') and r.endswith('.rnk')]

for rnk in rnks:
	# Run GSEA preranked on a certain rnk list
	"""
	prerank_results = gseapy.prerank(
		rnk='data/rnk_lists/' + rnk,
		gene_sets='data/gene_sets/c2.cp.kegg.v6.0.symbols.gmt', 
		outdir='out/Prerank_KEGG/' + rnk[:-4], 
		permutation_num=100, 
		graph_num=50,
		processes=4
	)
	"""

	# Run GSEA preranked on a certain rnk list
	gseapen_results = gseapy.gsea_pen(
		rnk='data/rnk_lists/' + rnk,
		gene_sets='data/gene_sets/c2.cp.kegg.v6.0.symbols.gmt', 
		background_rnks='out/all_443_lists_permuted_x100.csv', 
		outdir='out/GSEAPEN_KEGG/' + rnk[:-4], 
		permutation_num=100,
		graph_num=50,
		processes=4
	)

print('Done!')