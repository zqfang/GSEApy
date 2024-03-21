import warnings
from typing import AnyStr, Dict, Iterable, List, Optional, Tuple, Union

import pandas as pd

from .__main__ import __version__
from .biomart import Biomart
from .enrichr import Enrichr
from .gsea import GSEA, Prerank, Replot
from .gsva import GSVA
from .msigdb import Msigdb
from .parser import get_library, get_library_name, read_gmt
from .plot import barplot, dotplot, enrichment_map, gseaplot, gseaplot2, heatmap
from .ssgsea import SingleSampleGSEA


def gsea(
    data: Union[pd.DataFrame, str],
    gene_sets: Union[List[str], str, Dict[str, str]],
    cls: Union[List[str], str],
    outdir: Optional[str] = None,
    min_size: int = 15,
    max_size: int = 500,
    permutation_num: int = 1000,
    weight: float = 1.0,
    permutation_type: str = "phenotype",
    method: str = "signal_to_noise",
    ascending: bool = False,
    threads: int = 4,
    figsize: Tuple[float, float] = (6.5, 6),
    format: str = "pdf",
    graph_num: int = 20,
    no_plot: bool = False,
    seed: int = 123,
    verbose: bool = False,
    *args,
    **kwargs,
) -> GSEA:
    """Run Gene Set Enrichment Analysis.

    :param data: Gene expression data table, Pandas DataFrame, gct file.

    :param gene_sets: Enrichr Library name or .gmt gene sets file or dict of gene sets. Same input with GSEA.

    :param cls: A list or a .cls file format required for GSEA.

    :param str outdir: Results output directory. If None, nothing will write to disk.

    :param int permutation_num: Number of permutations. Default: 1000.
                                Minimial possible nominal p-value is about 1/nperm.

    :param str permutation_type: Type of permutation reshuffling,
                                 choose from {"phenotype": 'sample.labels' , "gene_set" : gene.labels}.

    :param int min_size: Minimum allowed number of genes from gene set also the data set. Default: 15.

    :param int max_size: Maximum allowed number of genes from gene set also the data set. Default: 500.

    :param float weight: Refer to :func:`algorithm.enrichment_score`. Default:1.

    :param method:  The method used to calculate a correlation or ranking. Default: 'signal_to_noise'.
                   Others methods are:

                   1. 'signal_to_noise'

                      You must have at least three samples for each phenotype to use this metric.
                      The larger the signal-to-noise ratio, the larger the differences of the means
                      (scaled by the standard deviations); that is, the more distinct
                      the gene expression is in each phenotype and the more the gene acts as a “class marker.”

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

    :param int threads: Number of threads you are going to use. Default: 4.

    :param list figsize: Matplotlib figsize, accept a tuple or list, e.g. [width,height]. Default: [6.5,6].

    :param str format: Matplotlib figure format. Default: 'pdf'.

    :param int graph_num: Plot graphs for top sets of each phenotype.

    :param bool no_plot: If equals to True, no figure will be drawn. Default: False.

    :param seed: Random seed. expect an integer. Default:None.

    :param bool verbose: Bool, increase output verbosity, print out progress of your job, Default: False.

    :return: Return a GSEA obj. All results store to a dictionary, obj.results,
             where contains::

                 | {
                 |  term: gene set name,
                 |  es: enrichment score,
                 |  nes: normalized enrichment score,
                 |  pval:  Nominal p-value (from the null distribution of the gene set,
                 |  fdr: FDR qvalue (adjusted False Discory Rate),
                 |  fwerp: Family wise error rate p-values,
                 |  tag %: Percent of gene set before running enrichment peak (ES),
                 |  gene %: Percent of gene list before running enrichment peak (ES),
                 |  lead_genes: leading edge genes (gene hits before running enrichment peak),
                 |  matched genes: genes matched to the data,
                 | }


    """
    if "processes" in kwargs:
        warnings.warn("processes is deprecated; use threads", DeprecationWarning, 2)
        threads = kwargs["processes"]
    if "weighted_score_type" in kwargs:
        warnings.warn(
            "weighted_score_type is deprecated; use weight", DeprecationWarning, 2
        )
        weight = kwargs["weighted_score_type"]

    gs = GSEA(
        data,
        gene_sets,
        cls,
        outdir,
        min_size,
        max_size,
        permutation_num,
        weight,
        permutation_type,
        method,
        ascending,
        threads,
        figsize,
        format,
        graph_num,
        no_plot,
        seed,
        verbose,
    )
    gs.run()

    return gs


def ssgsea(
    data: Union[pd.Series, pd.DataFrame, str],
    gene_sets: Union[List[str], str, Dict[str, str]],
    outdir: Optional[str] = None,
    sample_norm_method: Optional[str] = "rank",
    correl_norm_type: Optional[str] = "rank",
    min_size: int = 15,
    max_size: int = 500,
    permutation_num: Optional[int] = None,
    weight: float = 0.25,
    ascending: bool = False,
    threads: int = 4,
    figsize: Tuple[float, float] = (6.5, 6),
    format: str = "pdf",
    graph_num: int = 20,
    no_plot: bool = True,
    seed: int = 123,
    verbose: bool = False,
    *args,
    **kwargs,
) -> SingleSampleGSEA:
    """Run Gene Set Enrichment Analysis with single sample GSEA tool

    :param data: Expression table, pd.Series, pd.DataFrame, GCT file, or .rnk file format.

    :param gene_sets: Enrichr Library name or .gmt gene sets file or dict of gene sets. Same input with GSEA.

    :param outdir: Results output directory. If None, nothing will write to disk.

    :param str sample_norm_method: Sample normalization method. Choose from {'rank', 'log', 'log_rank', None}. Default: rank.
               this argument will be used for ordering genes.

               1. 'rank': Rank your expression data, and transform by 10000*rank_dat/gene_numbers
               2. 'log' : Do not rank, but transform data by log(data + exp(1)), while data = data[data<1] =1.
               3. 'log_rank': Rank your expression data, and transform by log(10000*rank_dat/gene_numbers+ exp(1))
               4. None or 'custom': Do nothing, and use your own rank value to calculate enrichment score.

    see here: https://github.com/GSEA-MSigDB/ssGSEAProjection-gpmodule/blob/master/src/ssGSEAProjection.Library.R, line 86

    :param str correl_norm_type: correlation normalization type. Choose from {'rank', 'symrank', 'zscore', None}. Default: rank.
            After ordering genes by sample_norm_method, further data transformed could be applied to get enrichment score.

            when weight == 0, sample_norm_method and correl_norm_type do not matter;
            when weight > 0, the combination of sample_norm_method and correl_norm_type
            dictate how the gene expression values in input data are transformed
            to obtain the score -- use this setting with care (the transformations
            can skew scores towards +ve or -ve values)

            sample_norm_method will first transformed and rank original data. the data is named correl_vector for each sample.
            then correl_vector is transformed again by

            1. correl_norm_type is None or 'rank' :  do nothing, genes are weighted by actual correl_vector.
            2. correl_norm_type =='symrank': symmetric ranking.
            3. correl_norm_type =='zscore':  standardizes the correl_vector before using them to calculate scores.


    :param int min_size: Minimum allowed number of genes from gene set also the data set. Default: 15.

    :param int max_size: Maximum allowed number of genes from gene set also the data set. Default: 2000.

    :param int permutation_num: For ssGSEA, default is 0.
                                However, if you try to use ssgsea method to get pval and fdr, set to an interger.

    :param str weight: Refer to :func:`algorithm.enrichment_score`. Default:0.25.

    :param bool ascending: Sorting order of rankings. Default: False.

    :param int threads: Number of threads you are going to use. Default: 4.

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

                 | {
                 |  term: gene set name,
                 |  es: enrichment score,
                 |  nes: normalized enrichment score,
                 |  pval:  Nominal p-value (from the null distribution of the gene set (if permutation_num > 0),
                 |  fdr: FDR qvalue (adjusted FDR) (if permutation_num > 0),
                 |  fwerp: Family wise error rate p-values (if permutation_num > 0),
                 |  tag %: Percent of gene set before running enrichment peak (ES),
                 |  gene %: Percent of gene list before running enrichment peak (ES),
                 |  lead_genes: leading edge genes (gene hits before running enrichment peak),
                 |  matched genes: genes matched to the data,
                 | }


    """
    if "processes" in kwargs:
        warnings.warn("processes is deprecated; use threads", DeprecationWarning, 2)
        threads = kwargs["processes"]
    if "weighted_score_type" in kwargs:
        warnings.warn(
            "weighted_score_type is deprecated; use weight", DeprecationWarning, 2
        )
        weight = kwargs["weighted_score_type"]

    ss = SingleSampleGSEA(
        data=data,
        gene_sets=gene_sets,
        outdir=outdir,
        sample_norm_method=sample_norm_method,
        correl_norm_type=correl_norm_type,
        min_size=min_size,
        max_size=max_size,
        permutation_num=permutation_num,
        weight=weight,
        ascending=ascending,
        threads=threads,
        figsize=figsize,
        format=format,
        graph_num=graph_num,
        no_plot=no_plot,
        seed=seed,
        verbose=verbose,
    )
    ss.run()
    return ss


def prerank(
    rnk: Union[pd.DataFrame, pd.Series, str],
    gene_sets: Union[List[str], str, Dict[str, str]],
    outdir: Optional[str] = None,
    pheno_pos: str = "Pos",
    pheno_neg: str = "Neg",
    min_size: int = 15,
    max_size: int = 500,
    permutation_num: int = 1000,
    weight: float = 1.0,
    ascending: bool = False,
    threads: int = 4,
    figsize: Tuple[float, float] = (6.5, 6),
    format: str = "pdf",
    graph_num: int = 20,
    no_plot: bool = False,
    seed: int = 123,
    verbose: bool = False,
    *arg,
    **kwargs,
) -> Prerank:
    """Run Gene Set Enrichment Analysis with pre-ranked correlation defined by user.

    :param rnk: pre-ranked correlation table or pandas DataFrame. Same input with ``GSEA`` .rnk file.

    :param gene_sets: Enrichr Library name or .gmt gene sets file or dict of gene sets. Same input with GSEA.

    :param outdir: results output directory. If None, nothing will write to disk.

    :param int permutation_num: Number of permutations. Default: 1000.
                                Minimial possible nominal p-value is about 1/nperm.

    :param int min_size: Minimum allowed number of genes from gene set also the data set. Default: 15.

    :param int max_size: Maximum allowed number of genes from gene set also the data set. Defaults: 500.

    :param str weight: Refer to :func:`algorithm.enrichment_score`. Default:1.

    :param bool ascending: Sorting order of rankings. Default: False.

    :param int threads: Number of threads you are going to use. Default: 4.

    :param list figsize: Matplotlib figsize, accept a tuple or list, e.g. [width,height]. Default: [6.5,6].

    :param str format: Matplotlib figure format. Default: 'pdf'.

    :param int graph_num: Plot graphs for top sets of each phenotype.

    :param bool no_plot: If equals to True, no figure will be drawn. Default: False.

    :param seed: Random seed. expect an integer. Default:None.

    :param bool verbose: Bool, increase output verbosity, print out progress of your job, Default: False.

    :return: Return a Prerank obj. All results store to  a dictionary, obj.results,
             where contains::

                 | {
                 |  term: gene set name,
                 |  es: enrichment score,
                 |  nes: normalized enrichment score,
                 |  pval:  Nominal p-value (from the null distribution of the gene set,
                 |  fdr: FDR qvalue (adjusted False Discory Rate),
                 |  fwerp: Family wise error rate p-values,
                 |  tag %: Percent of gene set before running enrichment peak (ES),
                 |  gene %: Percent of gene list before running enrichment peak (ES),
                 |  lead_genes: leading edge genes (gene hits before running enrichment peak),
                 |  matched genes: genes matched to the data,
                 | }


    """
    if "processes" in kwargs:
        warnings.warn("processes is deprecated; use threads", DeprecationWarning, 2)
        threads = kwargs["processes"]

    if "weighted_score_type" in kwargs:
        warnings.warn(
            "weighted_score_type is deprecated; use weight", DeprecationWarning, 2
        )
        weight = kwargs["weighted_score_type"]

    pre = Prerank(
        rnk,
        gene_sets,
        outdir,
        pheno_pos,
        pheno_neg,
        min_size,
        max_size,
        permutation_num,
        weight,
        ascending,
        threads,
        figsize,
        format,
        graph_num,
        no_plot,
        seed,
        verbose,
    )
    pre.run()
    return pre


def replot(
    indir,
    outdir="GSEA_Replot",
    weight=1,
    min_size=3,
    max_size=1000,
    figsize=(6.5, 6),
    format="pdf",
    verbose=False,
    *args,
    **kwargs,
):
    """The main function to reproduce GSEA desktop outputs.

    :param indir: GSEA desktop results directory. In the sub folder, you must contain edb file folder.

    :param outdir: Output directory.

    :param float weight: weighted score type. choose from {0,1,1.5,2}. Default: 1.

    :param list figsize: Matplotlib output figure figsize. Default: [6.5,6].

    :param str format: Matplotlib output figure format. Default: 'pdf'.

    :param int min_size: Min size of input genes presented in Gene Sets. Default: 3.

    :param int max_size: Max size of input genes presented in Gene Sets. Default: 5000.
                     You are not encouraged to use min_size, or max_size argument in :func:`replot` function.
                     Because gmt file has already been filtered.

    :param verbose: Bool, increase output verbosity, print out progress of your job, Default: False.

    :return: Generate new figures with selected figure format. Default: 'pdf'.

    """
    if "weighted_score_type" in kwargs:
        warnings.warn(
            "weighted_score_type is deprecated; use weight", DeprecationWarning, 2
        )
        weight = kwargs["weighted_score_type"]

    rep = Replot(indir, outdir, weight, min_size, max_size, figsize, format, verbose)
    rep.run()

    return


def enrichr(
    gene_list: Iterable[str],
    gene_sets: Union[List[str], str, Dict[str, str]],
    organism: str = "human",
    outdir: Optional[str] = None,
    background: Union[List[str], int, str] = None,
    cutoff: float = 0.05,
    format: str = "pdf",
    figsize: Tuple[float, float] = (6.5, 6),
    top_term: int = 10,
    no_plot: bool = False,
    verbose: bool = False,
) -> Enrichr:
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

    :param outdir:   Output file directory

    :param background: int, list, str.
                       Background genes. This argument works only if `gene_sets` has a type Dict or gmt file.
                       If your input are just Enrichr library names, this argument will be ignored.


                       However, this argument is not straightforward when `gene_sets` is given a custom input (a gmt file or dict).

                       By default, all genes listed in the `gene_sets` input will be used as background.

                       There are 3 ways to tune this argument:

                       (1) (Recommended) Input a list of background genes: ['gene1', 'gene2',...]
                           The background gene list is defined by your experment. e.g. the expressed genes in your RNA-seq.
                           The gene identifer in gmt/dict should be the same type to the backgound genes.

                       (2) Specify a number: e.g. 20000. (the number of total expressed genes).
                           This works, but not recommend. It assumes that all your genes could be found in background.
                           If genes exist in gmt but not included in background provided,
                           they will affect the significance of the statistical test.

                       (3) Set a Biomart dataset name: e.g. "hsapiens_gene_ensembl"
                           The background will be all annotated genes from the `BioMart datasets` you've choosen.
                           The program will try to retrieve the background information automatically.

                           Enrichr module use the code below to get the background genes:
                            >>> from gseapy.parser import Biomart
                            >>> bm = Biomart()
                            >>> df = bm.query(dataset=background, #  e.g. 'hsapiens_gene_ensembl'
                                         attributes=['ensembl_gene_id', 'external_gene_name', 'entrezgene_id'],
                                         filename=f'~/.cache/gseapy/{background}.background.genes.txt')
                            >>> df.dropna(subset=["entrezgene_id"], inplace=True)

                           So only genes with entrezid above will be the background genes if not input specify by user.

    :param cutoff:   Show enriched terms which Adjusted P-value < cutoff.
                     Only affects the output figure, not the final output file. Default: 0.05
    :param format:  Output figure format supported by matplotlib,('pdf','png','eps'...). Default: 'pdf'.

    :param figsize: Matplotlib figsize, accept a tuple or list, e.g. (width,height). Default: (6.5,6).

    :param bool no_plot: If equals to True, no figure will be drawn. Default: False.

    :param bool verbose: Increase output verbosity, print out progress of your job, Default: False.

    :return: An Enrichr object, which obj.res2d stores your last query, obj.results stores your all queries.

    """
    enr = Enrichr(
        gene_list,
        gene_sets,
        organism,
        outdir,
        background,
        cutoff,
        format,
        figsize,
        top_term,
        no_plot,
        verbose,
    )
    # set organism
    enr.set_organism()
    enr.run()

    return enr


def enrich(
    gene_list: Iterable[str],
    gene_sets: Union[List[str], str, Dict[str, str]],
    background: Union[List[str], int, str] = None,
    outdir: Optional[str] = None,
    cutoff: float = 0.05,
    format: str = "pdf",
    figsize: Tuple[float, float] = (6.5, 6),
    top_term: int = 10,
    no_plot: bool = False,
    verbose: bool = False,
) -> Enrichr:
    """Perform over-representation analysis (hypergeometric test).

    :param gene_list: str, list, tuple, series, dataframe. Also support input txt file with one gene id per row.
                      The input `identifier` should be the same type to `gene_sets`.

    :param gene_sets: str, list, tuple of Enrichr Library name(s).
                      or custom defined gene_sets (dict, or gmt file).

                      Examples:

                        dict: gene_sets={'A':['gene1', 'gene2',...],
                                        'B':['gene2', 'gene4',...], ...}
                        gmt: "genes.gmt"

    :param outdir:   Output file directory

    :param background: None | int | list | str.
                       Background genes. This argument works only if `gene_sets` has a type Dict or gmt file.

                       However, this argument is not straightforward when `gene_sets` is given a custom input (a gmt file or dict).

                       By default, all genes listed in the `gene_sets` input will be used as background.

                       There are 3 ways to tune this argument:

                       (1) (Recommended) Input a list of background genes: ['gene1', 'gene2',...]
                           The background gene list is defined by your experment. e.g. the expressed genes in your RNA-seq.
                           The gene identifer in gmt/dict should be the same type to the backgound genes.

                       (2) Specify a number: e.g. 20000. (the number of total expressed genes).
                           This works, but not recommend. It assumes that all your genes could be found in background.
                           If genes exist in gmt but not included in background provided,
                           they will affect the significance of the statistical test.

                       (3) Set a Biomart dataset name: e.g. "hsapiens_gene_ensembl"
                           The background will be all annotated genes from the `BioMart datasets` you've choosen.
                           The program will try to retrieve the background information automatically.

                           Enrichr module use the code below to get the background genes:
                            >>> from gseapy.parser import Biomart
                            >>> bm = Biomart()
                            >>> df = bm.query(dataset=background, #  e.g. 'hsapiens_gene_ensembl'
                                         attributes=['ensembl_gene_id', 'external_gene_name', 'entrezgene_id'],
                                         filename=f'~/.cache/gseapy/{background}.background.genes.txt')
                            >>> df.dropna(subset=["entrezgene_id"], inplace=True)

                           So only genes with entrezid above will be the background genes if not input specify by user.

    :param cutoff:   Show enriched terms which Adjusted P-value < cutoff.
                     Only affects the output figure, not the final output file. Default: 0.05
    :param format:  Output figure format supported by matplotlib,('pdf','png','eps'...). Default: 'pdf'.

    :param figsize: Matplotlib figsize, accept a tuple or list, e.g. (width,height). Default: (6.5,6).

    :param bool no_plot: If equals to True, no figure will be drawn. Default: False.

    :param bool verbose: Increase output verbosity, print out progress of your job, Default: False.

    :return: An Enrichr object, which obj.res2d stores your last query, obj.results stores your all queries.

    """
    organism = "human"  # has not any effects here
    _gene_sets = gene_sets
    if isinstance(gene_sets, str):
        _gene_sets = gene_sets.split(",")

    enr = Enrichr(
        gene_list,
        _gene_sets,
        organism,
        outdir,
        background,
        cutoff,
        format,
        figsize,
        top_term,
        no_plot,
        verbose,
    )
    # set organism
    enr.set_organism()
    enr.run()

    return enr


def gsva(
    data: Union[pd.DataFrame, pd.Series, str],
    gene_sets: Union[List[str], str, Dict[str, str]],
    outdir: Optional[str] = None,
    kcdf: Optional[str] = "Gaussian",
    weight: float = 1.0,
    mx_diff: bool = True,
    abs_rnk: bool = False,
    min_size: int = 15,
    max_size: int = 1000,
    threads: int = 4,
    seed: int = 123,
    verbose: bool = False,
    **kwargs,
) -> GSVA:
    """Run Gene Set Enrichment Analysis with single sample GSEA tool

    :param data: Expression table, pd.Series, pd.DataFrame, GCT file

    :param gene_sets: Enrichr Library name or .gmt gene sets file or dict of gene sets. Same input with GSEA.

    :param outdir: Results output directory. If None, nothing will write to disk.

    :param str kcdf: "Gaussian", "Possion" or None.
                     The non-parametric estimation of the cumulative distribution function

    :param str weight: tau of gsva. Default: 1.

    :param bool mx_diff: Offers two approaches to calculate the enrichment statistic (ES) from the KS random walk statistic.
                         If True (default), ES is calculated as the difference
                         between the maximum positive and negative random walk deviations.
                         If False, ES is calculated as the maximum positive to 0.

    :param bool abs_rnk: Flag used only when mx_diff=True.
                         If False (default), the enrichment statistic (ES) is calculated taking the magnitude difference
                         between the largest positive and negative random walk deviations.
                         If True, feature sets with features enriched on either extreme (high or low)
                         will be regarded as 'highly' activated.

    see GSVA doc: https://rdrr.io/bioc/GSVA/man/gsva.html

    :param int min_size: Minimum allowed number of genes from gene set also the data set. Default: 15.

    :param int max_size: Maximum allowed number of genes from gene set also the data set. Default: 1000.


    :param int threads: Number of threads you are going to use. Default: 4.

    :param seed: Random seed. expect an integer. Default:None.

    :param bool verbose: Bool, increase output verbosity, print out progress of your job, Default: False.

    :return: Return a GSVA obj.
             All results can be accessed by obj.results or obj.res2d.
    """
    gv = GSVA(
        data=data,
        gene_sets=gene_sets,
        outdir=outdir,
        kcdf=kcdf,
        weight=weight,
        mx_diff=mx_diff,
        abs_rnk=abs_rnk,
        min_size=min_size,
        max_size=max_size,
        threads=threads,
        seed=seed,
        verbose=verbose,
    )
    gv.run()
    return gv


__all__ = [
    "dotplot",
    "barplot",
    "enrichment_map",
    "heatmap",
    "gseaplot",
    "gseaplot2",
    "replot",
    "prerank",
    "gsea",
    "gsva",
    "ssgsea",
    "enrichr",
    "enrich",
    "Replot",
    "Prerank",
    "GSEA",
    "GSVA",
    "SingleSampleGSEA",
    "Enrichr",
    "Biomart",
    "Msigdb",
    "get_library",
    "get_library_name",
    "read_gmt",
]
