from .gsea import Replot, Prerank, GSEA, SingleSampleGSEA
from .enrichr import Enrichr
from .parser import get_library_name
from .plot import dotplot, barplot, heatmap, gseaplot
from .__main__ import __version__



def gsea(data, gene_sets, cls, outdir='GSEA_', min_size=15, max_size=500, permutation_num=1000,
         weighted_score_type=1, permutation_type='gene_set', method='log2_ratio_of_classes',
         ascending=False, processes=1, figsize=(6.5, 6), format='pdf',
         graph_num=20, no_plot=False, seed=123, verbose=False):
    """ Run Gene Set Enrichment Analysis.

    :param data: Gene expression data table, Pandas DataFrame, gct file.
    :param gene_sets: Enrichr Library name or .gmt gene sets file or dict of gene sets. Same input with GSEA.
    :param cls: A list or a .cls file format required for GSEA.
    :param str outdir: Results output directory.
    :param int permutation_num: Number of permutations for significance computation. Default: 1000.
    :param str permutation_type: Permutation type, "phenotype" for phenotypes, "gene_set" for genes.
    :param int min_size: Minimum allowed number of genes from gene set also the data set. Default: 15.
    :param int max_size: Maximum allowed number of genes from gene set also the data set. Default: 500.
    :param float weighted_score_type: Refer to :func:`algorithm.enrichment_score`. Default:1.
    :param method:  The method used to calculate a correlation or ranking. Default: 'log2_ratio_of_classes'.
                   Others methods are:

                   1. 'signal_to_noise'

                      You must have at least three samples for each phenotype to use this metric.
                      The larger the signal-to-noise ratio, the larger the differences of the means (scaled by the standard deviations);
                      that is, the more distinct the gene expression is in each phenotype and the more the gene acts as a “class marker.”

                   2. 't_test'

                      Uses the difference of means scaled by the standard deviation and number of samples.
                      Note: You must have at least three samples for each phenotype to use this metric.
                      The larger the tTest ratio, the more distinct the gene expression is in each phenotype
                      and the more the gene acts as a “class marker.”

                   3. 'ratio_of_classes' (also referred to as fold change).

                      Uses the ratio of class means to calculate fold change for natural scale data.

                   4. 'diff_of_classes'


                      Uses the difference of class means to calculate fold change for nature scale data


                   5. 'log2_ratio_of_classes'

                      Uses the log2 ratio of class means to calculate fold change for natural scale data.
                      This is the recommended statistic for calculating fold change for log scale data.


    :param bool ascending: Sorting order of rankings. Default: False.
    :param int processes: Number of Processes you are going to use. Default: 1.
    :param list figsize: Matplotlib figsize, accept a tuple or list, e.g. [width,height]. Default: [6.5,6].
    :param str format: Matplotlib figure format. Default: 'pdf'.
    :param int graph_num: Plot graphs for top sets of each phenotype.
    :param bool no_plot: If equals to True, no figure will be drawn. Default: False.
    :param seed: Random seed. expect an integer. Default:None.
    :param bool verbose: Bool, increase output verbosity, print out progress of your job, Default: False.

    :return: Return a GSEA obj. All results store to a dictionary, obj.results,
             where contains::

                 | {es: enrichment score,
                 |  nes: normalized enrichment score,
                 |  p: P-value,
                 |  fdr: FDR,
                 |  size: gene set size,
                 |  matched_size: genes matched to the data,
                 |  genes: gene names from the data set
                 |  ledge_genes: leading edge genes}


    """
    gs = GSEA(data, gene_sets, cls, outdir, min_size, max_size, permutation_num,
              weighted_score_type, permutation_type, method, ascending, processes,
              figsize, format, graph_num, no_plot, seed, verbose)
    gs.run()

    return gs


def ssgsea(data, gene_sets, outdir="ssGSEA_", sample_norm_method='rank', min_size=15, max_size=2000,
           permutation_num=0, weighted_score_type=0.25, scale=True, ascending=False, processes=1,
           figsize=(7, 6), format='pdf', graph_num=20, no_plot=True, seed=123, verbose=False):
    """Run Gene Set Enrichment Analysis with single sample GSEA tool

    :param data: Expression table, pd.Series, pd.DataFrame, GCT file, or .rnk file format.
    :param gene_sets: Enrichr Library name or .gmt gene sets file or dict of gene sets. Same input with GSEA.
    :param outdir: Results output directory.
    :param str sample_norm_method: "Sample normalization method. Choose from {'rank', 'log', 'log_rank'}. Default: rank.

               1. 'rank': Rank your expression data, and transform by 10000*rank_dat/gene_numbers
               2. 'log' : Do not rank, but transform data by log(data + exp(1)), while data = data[data<1] =1.
               3. 'log_rank': Rank your expression data, and transform by log(10000*rank_dat/gene_numbers+ exp(1))
               4. 'custom': Do nothing, and use your own rank value to calculate enrichment score.

    see here: https://github.com/GSEA-MSigDB/ssGSEAProjection-gpmodule/blob/master/src/ssGSEAProjection.Library.R, line 86

    :param int min_size: Minimum allowed number of genes from gene set also the data set. Default: 15.
    :param int max_size: Maximum allowed number of genes from gene set also the data set. Default: 2000.
    :param int permutation_num: Number of permutations for significance computation. Default: 0.
    :param str weighted_score_type: Refer to :func:`algorithm.enrichment_score`. Default:0.25.
    :param bool scale: If True, normalize the scores by number of genes in the gene sets.
    :param bool ascending: Sorting order of rankings. Default: False.
    :param int processes: Number of Processes you are going to use. Default: 1.
    :param list figsize: Matplotlib figsize, accept a tuple or list, e.g. [width,height]. Default: [7,6].
    :param str format: Matplotlib figure format. Default: 'pdf'.
    :param int graph_num: Plot graphs for top sets of each phenotype.
    :param bool no_plot: If equals to True, no figure will be drawn. Default: False.
    :param seed: Random seed. expect an integer. Default:None.
    :param bool verbose: Bool, increase output verbosity, print out progress of your job, Default: False.

    :return: Return a ssGSEA obj. 
             All results store to  a dictionary, access enrichment score by obj.resultsOnSamples,
             and normalized enrichment score by obj.res2d.
             if permutation_num > 0, additional results contain::

                 | {es: enrichment score,
                 |  nes: normalized enrichment score,
                 |  p: P-value,
                 |  fdr: FDR,
                 |  size: gene set size,
                 |  matched_size: genes matched to the data,
                 |  genes: gene names from the data set
                 |  ledge_genes: leading edge genes, if permutation_num >0}


    """
    ss = SingleSampleGSEA(data, gene_sets, outdir, sample_norm_method, min_size, max_size,
                          permutation_num, weighted_score_type, scale, ascending,
                          processes, figsize, format, graph_num, no_plot, seed, verbose)
    ss.run()
    return ss


def prerank(rnk, gene_sets, outdir='GSEA_Prerank', pheno_pos='Pos', pheno_neg='Neg',
            min_size=15, max_size=500, permutation_num=1000, weighted_score_type=1,
            ascending=False, processes=1, figsize=(6.5, 6), format='pdf',
            graph_num=20, no_plot=False, seed=123, verbose=False):
    """ Run Gene Set Enrichment Analysis with pre-ranked correlation defined by user.

    :param rnk: pre-ranked correlation table or pandas DataFrame. Same input with ``GSEA`` .rnk file.
    :param gene_sets: Enrichr Library name or .gmt gene sets file or dict of gene sets. Same input with GSEA.
    :param outdir: results output directory.
    :param int permutation_num: Number of permutations for significance computation. Default: 1000.
    :param int min_size: Minimum allowed number of genes from gene set also the data set. Default: 15.
    :param int max_size: Maximum allowed number of genes from gene set also the data set. Defaults: 500.
    :param str weighted_score_type: Refer to :func:`algorithm.enrichment_score`. Default:1.
    :param bool ascending: Sorting order of rankings. Default: False.
    :param int processes: Number of Processes you are going to use. Default: 1.
    :param list figsize: Matplotlib figsize, accept a tuple or list, e.g. [width,height]. Default: [6.5,6].
    :param str format: Matplotlib figure format. Default: 'pdf'.
    :param int graph_num: Plot graphs for top sets of each phenotype.
    :param bool no_plot: If equals to True, no figure will be drawn. Default: False.
    :param seed: Random seed. expect an integer. Default:None.
    :param bool verbose: Bool, increase output verbosity, print out progress of your job, Default: False.

    :return: Return a Prerank obj. All results store to  a dictionary, obj.results,
             where contains::

                 | {es: enrichment score,
                 |  nes: normalized enrichment score,
                 |  p: P-value,
                 |  fdr: FDR,
                 |  size: gene set size,
                 |  matched_size: genes matched to the data,
                 |  genes: gene names from the data set
                 |  ledge_genes: leading edge genes}


    """
    pre = Prerank(rnk, gene_sets, outdir, pheno_pos, pheno_neg,
                  min_size, max_size, permutation_num, weighted_score_type,
                  ascending, processes, figsize, format, graph_num, no_plot, seed, verbose)
    pre.run()
    return pre


def replot(indir, outdir='GSEA_Replot', weighted_score_type=1,
           min_size=3, max_size=1000, figsize=(6.5, 6), format='pdf', verbose=False):
    """The main function to reproduce GSEA desktop outputs.

    :param indir: GSEA desktop results directory. In the sub folder, you must contain edb file folder.
    :param outdir: Output directory.
    :param float weighted_score_type: weighted score type. choose from {0,1,1.5,2}. Default: 1.
    :param list figsize: Matplotlib output figure figsize. Default: [6.5,6].
    :param str format: Matplotlib output figure format. Default: 'pdf'.
    :param int min_size: Min size of input genes presented in Gene Sets. Default: 3.
    :param int max_size: Max size of input genes presented in Gene Sets. Default: 5000.
                     You are not encouraged to use min_size, or max_size argument in :func:`replot` function.
                     Because gmt file has already been filtered.
    :param verbose: Bool, increase output verbosity, print out progress of your job, Default: False.

    :return: Generate new figures with selected figure format. Default: 'pdf'.

    """
    rep = Replot(indir, outdir, weighted_score_type,
                 min_size, max_size, figsize, format, verbose)
    rep.run()

    return


def enrichr(gene_list, gene_sets, organism='human', description='',
            outdir='Enrichr', background='hsapiens_gene_ensembl', cutoff=0.05,
            format='pdf', figsize=(8, 6), top_term=10, no_plot=False, verbose=False):
    """Enrichr API.

    :param gene_list: str, list, tuple, series, dataframe. Also support input txt file with one gene id per row. 
                      The input `identifier` should be the same type to `gene_sets`.

    :param gene_sets: str, list, tuple of Enrichr Library name(s). 
                      or custom defined gene_sets (dict, or gmt file). 

                      Examples: 

                      Input Enrichr Libraries (https://maayanlab.cloud/Enrichr/#stats):
                        str: 'KEGG_2016'
                        list: ['KEGG_2016','KEGG_2013']
                        Use comma to separate each other, e.g. "KEGG_2016,huMAP,GO_Biological_Process_2018"

                      Input custom files:
                        dict: gene_sets={'A':['gene1', 'gene2',...],
                                        'B':['gene2', 'gene4',...], ...}
                        gmt: "genes.gmt"

                      see also the online docs: 
                      https://gseapy.readthedocs.io/en/latest/gseapy_example.html#2.-Enrichr-Example


    :param organism: Enrichr supported organism. Select from (human, mouse, yeast, fly, fish, worm).
                     This argument only affects the Enrichr library names you've chosen.
                     No any affects to gmt or dict input of `gene_sets`.

                     see here for more details: https://maayanlab.cloud/modEnrichr/.

    :param description: optional. name of the job.
    :param outdir:   Output file directory

    :param background: int, list, str. Please ignore this argument if your input are just Enrichr library names.

                       However, this argument is not straightforward when `gene_sets` is given a custom input (a gmt file or dict).
                       There are 3 ways to set this argument:

                       (1) (Recommended) Input a list of background genes. 
                           The background gene list is defined by your experment. e.g. the expressed genes in your RNA-seq.
                           The gene identifer in gmt/dict should be the same type to the backgound genes.  

                       (2) Specify a number, e.g. the number of total expressed genes.
                           This works, but not recommend. It assumes that all your genes could be found in background.
                           If genes exist in gmt but not included in background, 
                           they will affect the significance of the statistical test.  

                       (3) (Default) Set a Biomart dataset name.
                           The background will be all annotated genes from the `BioMart datasets` you've choosen. 
                           The program will try to retrieve the background information automatically.

                           Please Use the example code below to choose the correct dataset name:
                            >>> from gseapy.parser import Biomart 
                            >>> bm = Biomart()
                            >>> datasets = bm.get_datasets(mart='ENSEMBL_MART_ENSEMBL')

    :param cutoff:   Show enriched terms which Adjusted P-value < cutoff. 
                     Only affects the output figure, not the final output file. Default: 0.05
    :param format:  Output figure format supported by matplotlib,('pdf','png','eps'...). Default: 'pdf'.
    :param figsize: Matplotlib figsize, accept a tuple or list, e.g. (width,height). Default: (6.5,6).
    :param bool no_plot: If equals to True, no figure will be drawn. Default: False.
    :param bool verbose: Increase output verbosity, print out progress of your job, Default: False.

    :return: An Enrichr object, which obj.res2d stores your last query, obj.results stores your all queries.

    """
    enr = Enrichr(gene_list, gene_sets, organism, description, outdir,
                  cutoff, background, format, figsize, top_term, no_plot, verbose)
    # set organism
    enr.set_organism()
    enr.run()

    return enr


__all__ = ['dotplot', 'barplot', 'heatmap', 'gseaplot',
           'replot', 'prerank', 'gsea', 'ssgsea',
           'enrichr', 'get_library_name']
