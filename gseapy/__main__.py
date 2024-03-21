import argparse as ap
import logging
import os
import sys

# ------------------------------------
# Main function
# ------------------------------------

# there is a bug in add_argument(required=True), for hacking, don't set metavar='' when required=True,
# or args = argparser.parse_args() will throw bugs!!!


__version__ = "1.1.3"


def main():
    """The Main function/pipeline for GSEApy."""

    # Parse options...
    argparser = prepare_argparser()
    args = argparser.parse_args()
    subcommand = args.subcommand_name

    if subcommand == "replot":
        # reproduce plots using GSEAPY
        from .gsea import Replot

        Replot(
            indir=args.indir,
            outdir=args.outdir,
            weight=args.weight,
            figsize=args.figsize,
            format=args.format,
            verbose=args.verbose,
        ).run()

    elif subcommand == "gsea":
        # compute using GSEAPY
        from .gsea import GSEA

        gs = GSEA(
            args.data,
            args.gmt,
            args.cls,
            args.outdir,
            args.mins,
            args.maxs,
            args.n,
            args.weight,
            args.type,
            args.method,
            args.ascending,
            args.threads,
            args.figsize,
            args.format,
            args.graph,
            args.noplot,
            args.seed,
            args.verbose,
        )
        gs.run()
    elif subcommand == "prerank":
        from .gsea import Prerank

        pre = Prerank(
            args.rnk,
            args.gmt,
            args.outdir,
            args.label[0],
            args.label[1],
            args.mins,
            args.maxs,
            args.n,
            args.weight,
            args.ascending,
            args.threads,
            args.figsize,
            args.format,
            args.graph,
            args.noplot,
            args.seed,
            args.verbose,
        )
        pre.run()

    elif subcommand == "ssgsea":
        from .ssgsea import SingleSampleGSEA

        ss = SingleSampleGSEA(
            data=args.data,
            gene_sets=args.gmt,
            outdir=args.outdir,
            sample_norm_method=args.norm,
            correl_norm_type=args.correl,
            min_size=args.mins,
            max_size=args.maxs,
            permutation_num=args.n,
            weight=args.weight,
            ascending=args.ascending,
            threads=args.threads,
            figsize=args.figsize,
            format=args.format,
            graph_num=args.graph,
            no_plot=args.noplot,
            seed=args.seed,
            verbose=args.verbose,
        )
        ss.run()

    elif subcommand == "gsva":
        from .gsva import GSVA

        gv = GSVA(
            data=args.data,
            gene_sets=args.gmt,
            outdir=args.outdir,
            kcdf=args.kcdf,
            weight=args.weight,
            mx_diff=args.mx_diff,
            abs_rnk=args.abs_rnk,
            min_size=args.mins,
            max_size=args.maxs,
            threads=args.threads,
            seed=args.seed,
            verbose=args.verbose,
        )
        gv.run()
    elif subcommand == "enrichr":
        # calling enrichr API
        from .enrichr import Enrichr

        enr = Enrichr(
            gene_list=args.gene_list,
            gene_sets=args.library,
            organism=args.organism,
            outdir=args.outdir,
            format=args.format,
            cutoff=args.thresh,
            background=args.bg,
            figsize=args.figsize,
            top_term=args.term,
            no_plot=args.noplot,
            verbose=args.verbose,
        )
        # set organism
        enr.set_organism()
        enr.run()
    elif subcommand == "biomart":
        from .biomart import Biomart

        # read input file or a argument
        name, value = args.filter
        if os.path.isfile(value):
            with open(value, "r") as val:
                lines = val.readlines()
            value = [l.strip() for l in lines]
        # run query
        bm = Biomart(host=args.host, verbose=args.verbose)
        bm.query(
            dataset=args.bg,
            attributes=args.attrs.split(","),
            filters={name: value},
            filename=args.ofile,
        )
    else:
        argparser.print_help()
        sys.exit(0)


def prepare_argparser():
    """Prepare argparser object. New options will be added in this function first."""
    description = "%(prog)s -- Gene Set Enrichment Analysis in Python"
    epilog = "For command line options of each command, type: %(prog)s COMMAND -h"

    # top-level parser
    argparser = ap.ArgumentParser(description=description, epilog=epilog)
    argparser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )
    subparsers = argparser.add_subparsers(
        dest="subcommand_name"
    )  # help="sub-command help")

    # command for 'gsea'
    add_gsea_parser(subparsers)
    # command for 'prerank'
    add_prerank_parser(subparsers)
    # command for 'ssgsea'
    add_singlesample_parser(subparsers)
    # command for 'gsva'
    add_gsva_parser(subparsers)
    # command for 'plot'
    add_plot_parser(subparsers)
    # command for 'enrichr'
    add_enrichr_parser(subparsers)
    # command for 'biomart'
    add_biomart_parser(subparsers)

    return argparser


def add_output_option(parser):
    """output option"""

    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        type=str,
        default="GSEApy_reports",
        metavar="",
        action="store",
        help="The GSEApy output directory. Default: the current working directory",
    )
    parser.add_argument(
        "-f",
        "--format",
        dest="format",
        type=str,
        metavar="",
        action="store",
        choices=("pdf", "png", "jpeg", "eps", "svg"),
        default="pdf",
        help="File extensions supported by Matplotlib active backend,\
                              choose from {'pdf', 'png', 'jpeg','ps', 'eps','svg'}. Default: 'pdf'.",
    )
    parser.add_argument(
        "--fs",
        "--figsize",
        action="store",
        nargs=2,
        dest="figsize",
        metavar=("width", "height"),
        type=float,
        default=(6.5, 6),
        help="The figsize keyword argument need two parameters to define. Default: (6.5, 6)",
    )
    parser.add_argument(
        "--graph",
        dest="graph",
        action="store",
        type=int,
        default=20,
        metavar="int",
        help="Numbers of top graphs produced. Default: 20",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        dest="noplot",
        default=False,
        help="Speed up computing by suppressing the plot output."
        + "This is useful only if data are interested. Default: False.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        dest="verbose",
        help="Increase output verbosity, print out progress of your job",
    )


def add_output_group(parser, required=True):
    """output group"""

    # output_group = parser.add_mutually_exclusive_group(required=required)
    # output_group.add_argument(
    #     "-o",
    #     "--ofile",
    #     dest="ofile",
    #     type=str,
    #     default="GSEApy_reports",
    #     help="Output file name. Mutually exclusive with --o-prefix.",
    # )
    # output_group.add_argument(
    #     "--o-prefix",
    #     dest="ofile",
    #     type=str,
    #     default="GSEApy_reports",
    #     help="Output file prefix. Mutually exclusive with -o/--ofile.",
    # )
    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        type=str,
        default="GSEApy_reports",
        metavar="",
        action="store",
        help="The GSEApy output directory. Default: the current working directory",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        dest="verbose",
        help="Increase output verbosity, print out progress of your job",
    )


def add_gsea_parser(subparsers):
    """Add main function 'gsea' argument parsers."""

    argparser_gsea = subparsers.add_parser(
        "gsea", help="Main GSEApy Function: run GSEApy instead of GSEA."
    )

    # group for input files
    group_input = argparser_gsea.add_argument_group("Input files arguments")
    group_input.add_argument(
        "-d",
        "--data",
        dest="data",
        action="store",
        type=str,
        required=True,
        help="Input gene expression dataset file in txt format.Same with GSEA.",
    )
    group_input.add_argument(
        "-c",
        "--cls",
        dest="cls",
        action="store",
        type=str,
        required=True,
        help="Input class vector (phenotype) file in CLS format. Same with GSEA.",
    )
    group_input.add_argument(
        "-g",
        "--gmt",
        dest="gmt",
        action="store",
        type=str,
        required=True,
        help="Gene set database in GMT format. Same with GSEA.",
    )
    group_input.add_argument(
        "-t",
        "--permu-type",
        action="store",
        dest="type",
        type=str,
        metavar="perType",
        choices=("gene_set", "phenotype"),
        default="phenotype",
        help="Type of permutation reshuffling, "
        + "Choose from {'phenotype': 'sample.labels' , 'gene_set' : gene.labels}. "
        + "Default: phenotype",
    )

    # group for output files
    group_output = argparser_gsea.add_argument_group("Output arguments")
    add_output_option(group_output)

    # group for General options.
    group_opt = argparser_gsea.add_argument_group("GSEA advanced arguments")
    group_opt.add_argument(
        "-n",
        "--permu-num",
        dest="n",
        action="store",
        type=int,
        default=1000,
        metavar="nperm",
        help="Number of random permutations. For calculating esnulls. Default: 1000",
    )
    group_opt.add_argument(
        "--min-size",
        dest="mins",
        action="store",
        type=int,
        default=15,
        metavar="int",
        help="Min size of input genes presented in Gene Sets. Default: 15",
    )
    group_opt.add_argument(
        "--max-size",
        dest="maxs",
        action="store",
        type=int,
        default=500,
        metavar="int",
        help="Max size of input genes presented in Gene Sets. Default: 500",
    )
    group_opt.add_argument(
        "-w",
        "--weight",
        action="store",
        dest="weight",
        default=1.0,
        type=float,
        metavar="float",
        help="Weighted_score of rank_metrics. For weighting input genes. Choose from {0, 1, 1.5, 2}. Default: 1",
    )
    group_opt.add_argument(
        "-m",
        "--method",
        action="store",
        dest="method",
        type=str,
        metavar="",
        choices=(
            "signal_to_noise",
            "s2n",
            "abs_signal_to_noise",
            "abs_s2n",
            "t_test",
            "ratio_of_classes",
            "diff_of_classes",
            "log2_ratio_of_classes",
        ),
        default="signal_to_noise",
        help="Methods to calculate correlations of ranking metrics. \
                           Choose from {'signal_to_noise', 'abs_signal_to_noise', 't_test', 'ratio_of_classes', 'diff_of_classes','log2_ratio_of_classes'}.\
                           Default: 'signal_to_noise'",
    )
    group_opt.add_argument(
        "-a",
        "--ascending",
        action="store_true",
        dest="ascending",
        default=False,
        help="Rank metric sorting order. If the -a flag was chosen, then ascending equals to True. Default: False.",
    )
    group_opt.add_argument(
        "-s",
        "--seed",
        dest="seed",
        action="store",
        type=int,
        default=123,
        metavar="",
        help="Number of random seed. Default: 123",
    )
    group_opt.add_argument(
        "-p",
        "--threads",
        dest="threads",
        action="store",
        type=int,
        default=4,
        metavar="procs",
        help="Number of threads you are going to use. Default: 4",
    )

    return


def add_prerank_parser(subparsers):
    """Add function 'prerank' argument parsers."""

    argparser_prerank = subparsers.add_parser(
        "prerank", help="Run GSEApy Prerank tool on preranked gene list."
    )

    # group for input files
    prerank_input = argparser_prerank.add_argument_group("Input files arguments")
    prerank_input.add_argument(
        "-r",
        "--rnk",
        dest="rnk",
        action="store",
        type=str,
        required=True,
        help="Ranking metric file in .rnk format. Same with GSEA.",
    )
    prerank_input.add_argument(
        "-g",
        "--gmt",
        dest="gmt",
        action="store",
        type=str,
        required=True,
        help="Gene set database in GMT format. Same with GSEA.",
    )
    prerank_input.add_argument(
        "-l",
        "--label",
        action="store",
        nargs=2,
        dest="label",
        metavar=("pos", "neg"),
        type=str,
        default=("Pos", "Neg"),
        help="The phenotype label argument need two parameters to define. Default: ('Pos','Neg')",
    )

    # group for output files
    prerank_output = argparser_prerank.add_argument_group("Output arguments")
    add_output_option(prerank_output)

    # group for General options.
    prerank_opt = argparser_prerank.add_argument_group("GSEA advanced arguments")
    prerank_opt.add_argument(
        "-n",
        "--permu-num",
        dest="n",
        action="store",
        type=int,
        default=1000,
        metavar="nperm",
        help="Number of random permutations. For calculating esnulls. Default: 1000",
    )
    prerank_opt.add_argument(
        "--min-size",
        dest="mins",
        action="store",
        type=int,
        default=15,
        metavar="int",
        help="Min size of input genes presented in Gene Sets. Default: 15",
    )
    prerank_opt.add_argument(
        "--max-size",
        dest="maxs",
        action="store",
        type=int,
        default=500,
        metavar="int",
        help="Max size of input genes presented in Gene Sets. Default: 500",
    )
    prerank_opt.add_argument(
        "-w",
        "--weight",
        action="store",
        dest="weight",
        default=1.0,
        type=float,
        metavar="float",
        help="Weighted_score of rank_metrics. For weighting input genes. Choose from {0, 1, 1.5, 2}. Default: 1",
    )
    prerank_opt.add_argument(
        "-a",
        "--ascending",
        action="store_true",
        dest="ascending",
        default=False,
        help="Rank metric sorting order. If the -a flag was chosen, then ascending equals to True. Default: False.",
    )
    prerank_opt.add_argument(
        "-s",
        "--seed",
        dest="seed",
        action="store",
        type=int,
        default=123,
        metavar="",
        help="Number of random seed. Default: 123",
    )
    prerank_opt.add_argument(
        "-p",
        "--threads",
        dest="threads",
        action="store",
        type=int,
        default=4,
        metavar="procs",
        help="Number of threads you are going to use. Default: 4",
    )

    return


def add_singlesample_parser(subparsers):
    """Add function 'singlesample' argument parsers."""

    argparser_gsea = subparsers.add_parser("ssgsea", help="Run Single Sample GSEA.")

    # group for input files
    group_input = argparser_gsea.add_argument_group("Input files arguments")
    group_input.add_argument(
        "-d",
        "--data",
        dest="data",
        action="store",
        type=str,
        required=True,
        help="Input gene expression dataset file in txt format. Same with GSEA.",
    )
    group_input.add_argument(
        "-g",
        "--gmt",
        dest="gmt",
        action="store",
        type=str,
        required=True,
        help="Gene set database in GMT format. Same with GSEA.",
    )
    # group for output files
    group_output = argparser_gsea.add_argument_group("Output arguments")
    add_output_option(group_output)

    # group for General options.
    group_opt = argparser_gsea.add_argument_group(
        "Single Sample GSEA advanced arguments"
    )
    group_opt.add_argument(
        "--sn",
        "--sample-norm",
        dest="norm",
        action="store",
        type=str,
        default="rank",
        metavar="normalize",
        choices=("rank", "log", "log_rank", "custom"),
        help="Sample normalization method. Choose from {'rank', 'log', 'log_rank','custom'}. Default: rank",
    )

    group_opt.add_argument(
        "-c",
        "--correl-type",
        dest="correl",
        action="store",
        type=str,
        default="rank",
        metavar="transform",
        choices=("rank", "symrank", "zscore"),
        help="Input data transformation after sample normalization. Choose from {'rank','symrank', 'zscore'}. Default: rank",
    )
    group_opt.add_argument(
        "--ns",
        "--no-scale",
        action="store_false",
        dest="scale",
        default=True,
        help="If the flag was set, don't normalize the enrichment scores by number of genes.",
    )
    group_opt.add_argument(
        "-n",
        "--permu-num",
        dest="n",
        action="store",
        type=int,
        default=0,
        metavar="nperm",
        help="Number of random permutations. For calculating esnulls. Default: 0",
    )
    group_opt.add_argument(
        "--min-size",
        dest="mins",
        action="store",
        type=int,
        default=15,
        metavar="int",
        help="Min size of input genes presented in Gene Sets. Default: 15",
    )
    group_opt.add_argument(
        "--max-size",
        dest="maxs",
        action="store",
        type=int,
        default=2000,
        metavar="int",
        help="Max size of input genes presented in Gene Sets. Default: 2000",
    )
    group_opt.add_argument(
        "-w",
        "--weight",
        action="store",
        dest="weight",
        default=0.25,
        type=float,
        metavar="float",
        help="Weighted_score of rank_metrics. For weighting input genes. Default: 0.25",
    )
    group_opt.add_argument(
        "-a",
        "--ascending",
        action="store_true",
        dest="ascending",
        default=False,
        help="Rank metric sorting order. If the -a flag was chosen, then ascending equals to True. Default: False.",
    )
    group_opt.add_argument(
        "-s",
        "--seed",
        dest="seed",
        action="store",
        type=int,
        default=123,
        metavar="",
        help="Number of random seed. Default: 123",
    )
    group_opt.add_argument(
        "-p",
        "--threads",
        dest="threads",
        action="store",
        type=int,
        default=4,
        metavar="procs",
        help="Number of Processes you are going to use. Default: 4",
    )

    return


def add_gsva_parser(subparsers):
    """Add function 'GSVA' argument parsers."""

    argparser_gsva = subparsers.add_parser("gsva", help="Run GSVA.")

    # group for input files
    group_input = argparser_gsva.add_argument_group("Input files arguments")
    group_input.add_argument(
        "-d",
        "--data",
        dest="data",
        action="store",
        type=str,
        required=True,
        help="Input gene expression dataset file in txt format. Same with GSEA.",
    )
    group_input.add_argument(
        "-g",
        "--gmt",
        dest="gmt",
        action="store",
        type=str,
        required=True,
        help="Gene set database in GMT format. Same with GSEA.",
    )
    # group for output files
    group_output = argparser_gsva.add_argument_group("Output arguments")
    # add_output_option(group_output)
    add_output_group(group_output)
    # group for General options.
    group_opt = argparser_gsva.add_argument_group("GSVA advanced arguments")
    group_opt.add_argument(
        "-m",
        "--mx-diff",
        dest="mx_diff",
        action="store_false",
        default=True,
        help="When set, ES is calculated as the maximum distance of the random walk from 0. Default: False",
    )

    group_opt.add_argument(
        "-k",
        "--kernel-cdf",
        dest="kcdf",
        action="store",
        type=str,
        default="Gaussian",
        metavar="",
        choices=("Gaussian", "Poisson", "None"),
        help="Gaussian is suitable when input expression values are continuous. "
        + "If input integer counts, then this argument should be set to 'Poisson'",
    )

    group_opt.add_argument(
        "-a",
        "--abs-ranking",
        action="store_true",
        dest="abs_rnk",
        default=False,
        help="Flag used only when --mx-diff is not set. When set, the original Kuiper statistic is used",
    )
    group_opt.add_argument(
        "--min-size",
        dest="mins",
        action="store",
        type=int,
        default=15,
        metavar="int",
        help="Min size of input genes presented in Gene Sets. Default: 15",
    )
    group_opt.add_argument(
        "--max-size",
        dest="maxs",
        action="store",
        type=int,
        default=2000,
        metavar="int",
        help="Max size of input genes presented in Gene Sets. Default: 2000",
    )
    group_opt.add_argument(
        "-w",
        "--weight",
        action="store",
        dest="weight",
        default=1,
        type=float,
        metavar="float",
        help="tau in the random walk performed by the gsva. Default: 1",
    )
    group_opt.add_argument(
        "-s",
        "--seed",
        dest="seed",
        action="store",
        type=int,
        default=123,
        metavar="",
        help="Number of random seed. Default: 123",
    )
    group_opt.add_argument(
        "-p",
        "--threads",
        dest="threads",
        action="store",
        type=int,
        default=4,
        metavar="int",
        help="Number of Processes you are going to use. Default: 4",
    )

    return


def add_plot_parser(subparsers):
    """Add function 'plot' argument parsers."""

    argparser_replot = subparsers.add_parser(
        "replot", help="Reproduce GSEA desktop output figures."
    )

    group_replot = argparser_replot.add_argument_group("Input arguments")

    group_replot.add_argument(
        "-i",
        "--indir",
        action="store",
        dest="indir",
        required=True,
        metavar="GSEA_dir",
        help="The GSEA desktop results directroy that you want to reproduce the figure ",
    )
    add_output_option(group_replot)
    # add_output_group( argparser_plot )
    group_replot.add_argument(
        "-w",
        "--weight",
        action="store",
        dest="weight",
        default=1.0,
        type=float,
        metavar="float",
        help="Weighted_score of rank_metrics. Please Use the same value in GSEA. Choose from (0, 1, 1.5, 2),default: 1",
    )

    return


def add_enrichr_parser(subparsers):
    """Add function 'enrichr' argument parsers."""

    argparser_enrichr = subparsers.add_parser(
        "enrichr", help="Using Enrichr API to perform GO analysis."
    )

    # group for required options.
    enrichr_opt = argparser_enrichr.add_argument_group("Input arguments")
    enrichr_opt.add_argument(
        "-i",
        "--input-list",
        action="store",
        dest="gene_list",
        type=str,
        required=True,
        metavar="IDs",
        help="Enrichr uses a list of gene names as input.",
    )
    enrichr_opt.add_argument(
        "-g",
        "--gene-sets",
        action="store",
        dest="library",
        type=str,
        required=True,
        metavar="GMT",
        help="Enrichr library name(s) required. Separate each name by comma.",
    )
    enrichr_opt.add_argument(
        "--org",
        "--organism",
        action="store",
        dest="organism",
        type=str,
        default="human",
        help="Enrichr supported organism name. Default: human. See here: https://amp.pharm.mssm.edu/modEnrichr.",
    )
    enrichr_opt.add_argument(
        "--ds",
        "--description",
        action="store",
        dest="descrip",
        type=str,
        default="enrichr",
        metavar="STRING",
        help="It is recommended to enter a short description for your list so that multiple lists \
                              can be differentiated from each other if you choose to save or share your list.",
    )
    enrichr_opt.add_argument(
        "--cut",
        "--cut-off",
        action="store",
        dest="thresh",
        metavar="float",
        type=float,
        default=0.05,
        help="Adjust-Pval cutoff, used for generating plots. Default: 0.05.",
    )
    enrichr_opt.add_argument(
        "--bg",
        "--background",
        action="store",
        dest="bg",
        default=None,
        metavar="BG",
        help="Choose background from one of the following. \n"
        + "(1) A BioMart Dataset name, e.g. 'hsapiens_gene_ensembl' . \n"
        + "(2) A total gene number, e.g. 20000. Only works for GMT file input. \n"
        + "(3) A text file contains the background gene list (one gene per row). Gene identifier should be the same to your input (-i). \n"
        + "(4) Default: None. It means all genes in the (-g) input as the background.",
    )
    enrichr_opt.add_argument(
        "-t",
        "--top-term",
        dest="term",
        action="store",
        type=int,
        default=10,
        metavar="int",
        help="Numbers of top terms shown in the plot. Default: 10",
    )
    # enrichr_opt.add_argument("--scale", dest = "scale", action="store", type=float, default=0.5, metavar='float',
    #                          help="scatter dot scale in the dotplot. Default: 0.5")
    # enrichr_opt.add_argument("--no-plot", action='store_true', dest='no_plot', default=False,
    #                           help="Suppress the plot output.This is useful only if data are interested. Default: False.")

    enrichr_output = argparser_enrichr.add_argument_group("Output figure arguments")
    add_output_option(enrichr_output)
    return


def add_biomart_parser(subparsers):
    """Add function 'biomart' argument parsers."""

    argparser_biomart = subparsers.add_parser(
        "biomart", help="Using BioMart API to convert gene ids."
    )

    # group for required options.
    biomart_opt = argparser_biomart.add_argument_group("Input arguments")
    biomart_opt.add_argument(
        "-f",
        "--filter",
        action="store",
        nargs=2,
        dest="filter",
        required=True,
        metavar=("NAME", "VALUE"),
        help="""Which filter to use. Input filter name, and value.
                                     If multi-value required, separate each value by comma.
                                     If value is a txt file, then one ID per row, exclude header.""",
    )
    biomart_opt.add_argument(
        "-a",
        "--attributes",
        action="store",
        dest="attrs",
        type=str,
        required=True,
        metavar="ATTR",
        help="Which attribute(s) to retrieve. Separate each attr by comma.",
    )
    biomart_opt.add_argument(
        "-o", "--ofile", dest="ofile", type=str, required=True, help="Output file name"
    )
    biomart_opt.add_argument(
        "-d",
        "--dataset",
        action="store",
        dest="bg",
        type=str,
        default="hsapiens_gene_ensembl",
        metavar="DATA",
        help="Which dataset to use. Default: hsapiens_gene_ensembl",
    )
    biomart_opt.add_argument(
        "--host",
        action="store",
        dest="host",
        type=str,
        default="www.ensembl.org",
        metavar="HOST",
        help="Which host to use. Select from {'www.ensembl.org', 'asia.ensembl.org', 'useast.ensembl.org'}.",
    )
    biomart_opt.add_argument(
        "-m",
        "--mart",
        action="store",
        dest="mart",
        type=str,
        metavar="MART",
        default="ENSEMBL_MART_ENSEMBL",
        help="Which mart to use. Default: ENSEMBL_MART_ENSEMBL.",
    )
    biomart_opt.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        dest="verbose",
        help="Increase output verbosity, print out progress of your job",
    )


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted me! ;-) Bye!\n")
        sys.exit(0)
