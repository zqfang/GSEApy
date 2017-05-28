
import sys, logging
import argparse as ap


# ------------------------------------
# Main function
# ------------------------------------


# there is a bug in add_argument(required=True), for hacking, don't set metavar='' when required=True,
# or args = argparser.parse_args() will throw bugs!!!


__version__ = '0.8.0'

def main():
    """The Main function/pipeline for GSEAPY."""

    # Parse options...
    argparser = prepare_argparser()
    args = argparser.parse_args()
    subcommand = args.subcommand_name

    if subcommand == "replot":
        # reproduce plots using GSEAPY
        from .gsea import replot

        replot(indir=args.indir, outdir=args.outdir, weight=args.weight,
               figsize=args.figsize, format=args.format, verbose=args.verbose)

    elif subcommand == "gsea":
        # compute using GSEAPY
        from .gsea import GSEA

        gs = GSEA(args.data, args.gmt, args.cls, args.outdir, args.mins, args.maxs, args.n, args.weight,
                  args.type, args.method, args.ascending, args.figsize, args.format, args.graph, args.seed, args.verbose)
        gs.run()
    elif subcommand == "prerank":
        from .gsea import Prerank
        
        pre = Prerank(args.rnk, args.gmt, args.outdir, args.label[0], args.label[1], args.mins, args.maxs, args.n, args.weight,
                args.ascending, args.figsize, args.format, args.graph, args.seed, args.verbose)
        pre.run()

    elif subcommand == "single":
        from .gsea import SingleSampleGSEA
        ss = SingleSampleGSEA(data=args.data, gene_sets=args.gmt, outdir=args.outdir,
                              min_size=args.mins, max_size=args.max, permutation_num=args.n, 
                              weighted_score_type=args.weight, ascending=args.ascending, 
                              figsize=args.figsize, format='pdf', graph_num=args.graph, 
                              seed=args.seed, verbose=args.verbose)
        ss.run()

    elif subcommand == "enrichr":
        # calling enrichr API
        from .enrichr import enrichr
        enrichr(gene_list= args.gene_list, description=args.descrip, gene_sets=args.library,
                outdir=args.outdir, format=args.format, cutoff=args.thresh, figsize=args.figsize,
                top_term=args.term, scale=args.scale, no_plot=args.no_plot, verbose=args.verbose)    
    else:
        argparser.print_help()
        sys.exit(0)


def prepare_argparser():
    """Prepare argparser object. New options will be added in this function first."""
    description = "%(prog)s -- Gene Set Enrichment Analysis in Python"
    epilog = "For command line options of each command, type: %(prog)s COMMAND -h"

    # top-level parser
    argparser = ap.ArgumentParser(description=description, epilog=epilog)
    argparser.add_argument("--version", action="version", version="%(prog)s "+ __version__)
    subparsers = argparser.add_subparsers(dest='subcommand_name') #help="sub-command help")

    # command for 'gsea'
    add_gsea_parser(subparsers)
    # command for 'prerank'
    add_prerank_parser(subparsers)
    # command for 'single'
    add_singlesample_parser(subparsers)
    # command for 'plot'
    add_plot_parser(subparsers)
    # command for 'enrichr'
    add_enrichr_parser(subparsers)

    return argparser

def add_output_option(parser):
    """output option"""

    parser.add_argument("-o", "--outdir", dest="outdir", type=str, default='GSEApy_reports',
                        metavar='', action="store", help="The GSEApy output directory. Default: the current working directory")
    parser.add_argument("-f", "--format", dest="format", type=str, metavar='', action="store",
                        choices=("pdf", "png", "jpeg", "eps","svg"), default="pdf",
                        help="File extensions supported by Matplotlib active backend,\
                              choose from {'pdf', 'png', 'jpeg','ps', 'eps','svg'}. Default: 'pdf'.")
    parser.add_argument("--figsize", action='store', nargs=2, dest='figsize',
                        metavar=('width', 'height'),type=float, default=(6.5, 6),
                        help="The figsize keyword argument need two parameter to define. Default: (6.5, 6)")
    parser.add_argument("-v", "--verbose", action="store_true", default=False, dest='verbose',
                        help="increase output verbosity, print out progress of your job", )
def add_output_group(parser, required=True):
    """output group"""

    output_group = parser.add_mutually_exclusive_group(required=required)
    output_group.add_argument("-o", "--ofile", dest="ofile", type=str, default='GSEApy_reports',
                              help="Output file name. Mutually exclusive with --o-prefix.")
    output_group.add_argument("--o-prefix", dest="ofile", type=str, default='GSEApy_reports',
                              help="Output file prefix. Mutually exclusive with -o/--ofile.")





def add_gsea_parser(subparsers):
    """Add main function 'gsea' argument parsers."""

    argparser_gsea = subparsers.add_parser("gsea", help="Main GSEApy Function: run GSEApy instead of GSEA.")

    # group for input files
    group_input = argparser_gsea.add_argument_group("Input files arguments")
    group_input.add_argument("-d", "--data", dest="data", action="store", type=str, required=True, 
                             help="Input gene expression dataset file in txt format.Same with GSEA.")
    group_input.add_argument("-c", "--cls", dest="cls", action="store", type=str, required=True,
                             help="Input class vector (phenotype) file in CLS format. Same with GSEA.")
    group_input.add_argument("-g", "--gmt", dest="gmt", action="store", type=str, required=True,
                             help="Gene set database in GMT format. Same with GSEA.")
    group_input.add_argument("-p", "--permu-type", action="store", dest="type", type=str, metavar='',
                             choices=("gene_set", "phenotype"), default="gene_set",
                             help="Permutation type. Same with GSEA, choose from {'gene_set', 'phenotype'}")

    # group for output files
    group_output = argparser_gsea.add_argument_group("Output arguments")
    add_output_option(group_output)

     # group for General options.
    group_opt = argparser_gsea.add_argument_group("GSEA advanced arguments")
    group_opt.add_argument("--min-size",  dest="mins", action="store", type=int, default=15, metavar='int',
                           help="Min size of input genes presented in Gene Sets. Default: 15")
    group_opt.add_argument("--max-size", dest = "maxs", action="store", type=int, default=500, metavar='int',
                           help="Max size of input genes presented in Gene Sets. Default: 500")
    group_opt.add_argument("-n", "--permu-num", dest = "n", action="store", type=int, default=1000, metavar='PerMut',
                           help="Number of random permutations. For calculating esnulls. Default: 1000")
    group_opt.add_argument("-w", "--weight", action='store', dest='weight', default=1.0, type=float, metavar='',
                           help='Weighted_score of rank_metrics.For weighting input genes. Choose from {0, 1, 1.5, 2},default: 1',)
    group_opt.add_argument("-m", "--method", action="store", dest="method", type=str, metavar='',
                           choices=("signal_to_noise", "t_test", "ratio_of_classes", "diff_of_classes", "log2_ratio_of_classes"),
                           default="log2_ratio_of_classes",
                           help="Methods to calculate correlations of ranking metrics. \
                           Choose from {'signal_to_noise', 't_test', 'ratio_of_classes', 'diff_of_classes','log2_ratio_of_classes'}.\
                           Default: 'log2_ratio_of_classes'")
    group_opt.add_argument("-a", "--ascending", action='store_true', dest='ascending', default=False,
                           help='Rank metric sorting order. If the -a flag was chosen, then ascending equals to True. Default: False.')
    group_opt.add_argument("-t", "--top-graphNum", dest = "graph", action="store", type=int, default=20, metavar='int',
                           help="Numbers of top graphs of each phenotype. Default: 20")
    group_opt.add_argument("-s", "--seed", dest = "seed", action="store", type=int, default=None, metavar='',
                           help="Number of random seed. Default: None")

    return
    
def add_prerank_parser(subparsers):
    """Add function 'prerank' argument parsers."""

    argparser_prerank = subparsers.add_parser("prerank", help="Run GSEApy Prerank tool on preranked gene list.")

    # group for input files
    prerank_input = argparser_prerank.add_argument_group("Input files arguments")
    prerank_input.add_argument("-r", "--rnk", dest="rnk", action="store", type=str, required=True,  
                             help="ranking dataset file in .rnk format.Same with GSEA.")
    prerank_input.add_argument("-g", "--gmt", dest="gmt", action="store", type=str, required=True, 
                             help="Gene set database in GMT format. Same with GSEA.")
    prerank_input.add_argument("-l", "--label", action='store', nargs=2, dest='label',
                             metavar=('pos', 'neg'), type=str, default=('Pos','Neg'),
                             help="The phenotype label argument need two parameter to define. Default: ('Pos','Neg')")

    # group for output files
    prerank_output = argparser_prerank.add_argument_group("Output arguments")
    add_output_option(prerank_output)

     # group for General options.
    prerank_opt = argparser_prerank.add_argument_group("GSEA advanced arguments")
    prerank_opt.add_argument("--min-size",  dest="mins", action="store", type=int, default=15, metavar='int',
                             help="Min size of input genes presented in Gene Sets. Default: 15")
    prerank_opt.add_argument("--max-size", dest = "maxs", action="store", type=int, default=500, metavar='int',
                             help="Max size of input genes presented in Gene Sets. Default: 500")
    prerank_opt.add_argument("-n", "--permu-num", dest = "n", action="store", type=int, default=1000, metavar='',
                             help="Number of random permutations. For calculating esnulls. Default: 1000")
    prerank_opt.add_argument("-w", "--weight", action='store', dest='weight', default=1.0, type=float, metavar='',
                             help='Weighted_score of rank_metrics.For weighting input genes. Choose from {0, 1, 1.5, 2},default: 1',)
    prerank_opt.add_argument("-a", "--ascending", action='store_true', dest='ascending', default=False,
                             help='Rank metric sorting order. If the -a flag was chosen, then ascending equals to True. Default: False.')
    prerank_opt.add_argument("-t", "--top-graphNum", dest = "graph", action="store", type=int, default=20, metavar='int',
                             help="Numbers of top graphs of each phenotype. Default: 20")
    prerank_opt.add_argument("-s", "--seed", dest = "seed", action="store", type=int, default=None, metavar='',
                             help="Number of random seed. Default: None")
    
    return

def add_singlesample_parser(subparsers):
    """Add function 'singlesample' argument parsers."""

    argparser_gsea = subparsers.add_parser("single", help="Run Single Sample GSEA.")

    # group for input files
    group_input = argparser_gsea.add_argument_group("Input files arguments")
    group_input.add_argument("-d", "--data", dest="data", action="store", type=str, required=True, 
                             help="Input gene expression dataset file in txt format. Same with GSEA.")
    group_input.add_argument("-g", "--gmt", dest="gmt", action="store", type=str, required=True,
                             help="Gene set database in GMT format. Same with GSEA.")
    # group for output files
    group_output = argparser_gsea.add_argument_group("Output arguments")
    add_output_option(group_output)

    # group for General options.
    group_opt = argparser_gsea.add_argument_group("GSEA advanced arguments")
    group_opt.add_argument("--min-size", dest="mins", action="store", type=int, default=15, metavar='int',
                           help="Min size of input genes presented in Gene Sets. Default: 15")
    group_opt.add_argument("--max-size", dest = "maxs", action="store", type=int, default=500,metavar='int',
                           help="Max size of input genes presented in Gene Sets. Default: 500")
    group_opt.add_argument("-n", "--permu-num", dest = "n", action="store", type=int, default=1000, metavar='PerMut',
                           help="Number of random permutations. For calculating esnulls. Default: 1000")
    group_opt.add_argument("-w", "--weight", action='store', dest='weight', default=0.25, type=float, metavar='weight',
                           help='Weighted_score of rank_metrics.For weighting input genes.default: 0.25',)
    group_opt.add_argument("-a", "--ascending", action='store_true', dest='ascending', default=False,
                           help='Rank metric sorting order. If the -a flag was chosen, then ascending equals to True. Default: False.')
    group_opt.add_argument("-t", "--top-graphNum", dest = "graph", action="store", type=int, default=20, metavar='int',
                           help="Numbers of top graphs of each phenotype. Default: 20")
    group_opt.add_argument("-s", "--seed", dest = "seed", action="store", type=int, default=None, metavar='',
                           help="Number of random seed. Default: None")

    return


def add_plot_parser(subparsers):
    """Add function 'plot' argument parsers."""

    argparser_replot = subparsers.add_parser("replot", help="Reproduce GSEA desktop output figures.")

    group_replot = argparser_replot.add_argument_group("Input files arguments")

    group_replot.add_argument("-i", "--indir", action="store", dest="indir", required=True, metavar='',
                              help="The GSEA desktop results directroy that you want to reproduce the figure ")
    add_output_option(group_replot)
    #add_output_group( argparser_plot )
    group_replot.add_argument("-w", "--weight", action='store', dest='weight', default=1.0, type=float, metavar='float',
                              help='Weighted_score of rank_metrics. Please Use the same value in GSEA. Choose from (0, 1, 1.5, 2),default: 1',)

    return
def add_enrichr_parser(subparsers):
    """Add function 'enrichr' argument parsers."""
    
    argparser_enrichr = subparsers.add_parser("enrichr", help="Using enrichr API to perform GO analysis.")

    # group for required options.
    enrichr_opt = argparser_enrichr.add_argument_group("Input arguments")
    enrichr_opt.add_argument("-i", "--input-list", action="store", dest="gene_list", type=str, required=True, metavar='file',
                              help="Enrichr uses a list of Entrez gene symbols as input.")
    enrichr_opt.add_argument("-g", "--gene-sets", action="store", dest="library", type=str, required=True, metavar='gmt',
                              help="Enrichr library name required. see online tool for libarry names.")
    enrichr_opt.add_argument("-d", "--description", action="store", dest="descrip", type=str, default='foo', metavar='',
                              help="It is recommended to enter a description for your list so that multiple lists \
                              can be differentiated from each other if you choose to save or share your list.") 
    enrichr_opt.add_argument("--cut-off", action="store", dest="thresh", metavar='float', type=float, default=0.05,
                              help="Adjust-Pval cutoff, used for generating plots. Default: 0.05.")
    enrichr_opt.add_argument("-t", "--top-term", dest="term", action="store", type=int, default=10, metavar='int',
                              help="Numbers of top terms showed in the plot. Default: 10")
    #enrichr_opt.add_argument("--scale", dest = "scale", action="store", type=float, default=0.5, metavar='float',
    #                          help="scatter dot scale in the dotplot. Default: 0.5")    
    enrichr_opt.add_argument("--no-plot", action='store_true', dest='no_plot', default=False, 
                              help="Suppress the plot output.This is useful only if data are intrested. Default: False.")


    enrichr_output = argparser_enrichr.add_argument_group("Output figure arguments")
    add_output_option(enrichr_output)


    return

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted me! ;-) Bye!\n")
        sys.exit(0)
