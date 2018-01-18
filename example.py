import gseapy
import numpy as np
from numpy.random import choice
from gseapy.algorithm import enrichment_score_pen
from collections import defaultdict
import csv
import sys, logging

"""
def make_background_dist(background_rnk_lists, nperm):

	# First gather all existing indices for each gene in the background lists
	gene_dists = defaultdict(list)
	for lst in background_rnk_lists:
		for i, gene in enumerate(lst):
			gene_dists[gene].append(i)

	# Then make n permutations of that
	n_genes = len(background_rnk_lists[0])
	tag_indicator = []
	logging.debug("permutation #")
	
	for i in range(nperm):
		logging.debug(i)
		unseen_idx = list(range(n_genes))
		genes = [None] * n_genes
		for gene in gene_dists.keys():
			gene_dist = gene_dists[gene]
			gene_idx = int(choice(gene_dist, 1))
			if gene_idx not in unseen_idx:
				diffs = [abs(gene_idx - us) for us in unseen_idx]
				gene_idx = unseen_idx[diffs.index(min(diffs))]
			unseen_idx.remove(gene_idx)
			genes[gene_idx] = gene
		tag_indicator.append(genes)
		assert len(unseen_idx) == 0

	return tag_indicator

#bg_lst = [['A', 'B', 'C', 'D', 'E', 'F'], ['F', 'B', 'C', 'A', 'E', 'D'], ['A', 'F', 'B', 'D', 'E', 'C']]
#permed = make_background_dist(bg_lst, 3)
#print(permed)

rnk='data/rnk_lists/GSE11791_DEG_Expt1_Control_vs_Group1_gene.rnk'
gene_sets='data/gene_sets/hallmarks_small.gmt'
background_rnks='data/bg_rnk_lists/mcf7_estrogen.txt'
outdir='out/MCF7_permuted_data'
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

with open("out/MCF7_22_background_lists_permuted_x100.csv", "w") as f:
	writer = csv.writer(f)
	writer.writerows(permed)

# Run GSEA preranked on a certain rnk list

"""
prerank_results = gseapy.prerank(
	rnk='data/rnk_lists/GSE11791_DEG_Expt1_Control_vs_Group1_gene.rnk',  # MCF7 estradiol
	gene_sets='data/gene_sets/hallmarks.gmt', 
	outdir='out/MCF7_prerank', 
	permutation_num=100, 
	graph_num=50,
	processes=4
)

"""
# Run GSEA preranked on a certain rnk list
gseapen_results = gseapy.gsea_pen(
	rnk='data/rnk_lists/GSE11791_DEG_Expt1_Control_vs_Group1_gene.rnk',  # MCF7 estradiol
	gene_sets='data/gene_sets/hallmarks.gmt', 
	background_rnks='out/MCF7_22_background_lists_permuted_x100.csv', 
	outdir='out/MCF7_gseapen', 
	permutation_num=100,
	graph_num=50,
	processes=4
)


# Edit these
rnk='data/rnk_lists/GSE11791_DEG_Expt1_Control_vs_Group1_gene.rnk'
gene_sets='data/gene_sets/hallmarks_small.gmt'
background_rnks='data/bg_rnk_lists/mcf7_estrogen_small.txt'
outdir='out/MCF7_gseapen'
permutation_num=1000
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

#pen.run()
dat2 = pen._load_ranking(rnk)
print(len(dat2))

# parse background gene lists
bg_lists = pen._load_background_ranking(background_rnks)
print(len(bg_lists[0]))

gmt = pen.load_gmt(gene_list=dat2.index.values, gmt=gene_sets)

# gsea_compute()
w = weighted_score_type
subsets = sorted(gmt.keys())
es = []
RES=[]
hit_ind=[]
esnull = [ [] for a in range(len(subsets)) ]

data=dat2
n=permutation_num
gmt=gmt
weighted_score_type=weighted_score_type
bg_lists=bg_lists
permutation_type='gsea_pen'

gl, cor_vec = data.index.values, data.values

# Try on just one gene set
subset = subsets[0]
rs = np.random.RandomState(seed)
es, esnull, hit_ind, RES = enrichment_score_pen(
	rnk_list=gl, 
	correl_vector=cor_vec, 
	gene_set=gmt.get(subset), 
	background_rnk_lists=bg_lists, 
	dist='original', 
	weighted_score_type=w, 
	nperm=n, rs=rs, single=False, scale=False)
"""
