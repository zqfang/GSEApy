
import sys
import argparse as ap


# ------------------------------------
# Main function
# ------------------------------------

__version__ = '0.6.1'

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
               figsize=args.figsize, format=args.format,)

    elif subcommand == "call":
        # compute using GSEAPY
        from .gsea import call

        call(args.data, args.gmt, args.cls, args.outdir, args.mins, args.maxs, args.n, args.weight,
            args.type, args.method, args.ascending, args.figsize, args.format, args.graph, args.seed)
    
    elif subcommand == "prerank":
        from .gsea import prerank
        
        prerank(args.rnk, args.gmt, args.outdir, args.label[0], args.label[1], args.mins, args.maxs, args.n, args.weight,
                args.ascending, args.figsize, args.format, args.graph, args.seed)
    
    elif subcommand == 'enrichr':
        # calling enrichr API
        from .enrichr import enrichr
        enrichr(gene_list= args.gene_list, description=args.description, gene_sets=args.library, outfile=args.ofile)    
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

    # command for 'call'
    add_call_parser(subparsers)
    # command for 'prerank'
    add_prerank_parser(subparsers)
    # command for 'plot'
    add_plot_parser(subparsers)
    # command for 'enrichr'
    add_enrichr_parser(subparsers)

    return argparser

def add_output_option(parser):
    """output option"""

    parser.add_argument("-o", "--outdir", dest="outdir", type=str, default='gseapy_out',
                        metavar='', action="store", help="The GSEAPY output directory. Default: the current working directory")
    parser.add_argument("-f", "--format", dest="format", type=str, metavar='', action="store",
                        choices=("pdf", "png", "jpeg", "eps"), default="png",
                        help="Format of output figures, choose from {'pdf', 'png', 'jpeg', 'eps'}. Default: 'pdf'.")
    parser.add_argument("--figsize", action='store', nargs=2, dest='figsize',
                        metavar=('width', 'height'), type=float, default=[6.5, 6],
                        help="The figsize keyword argument need two parameter to define. Default: [6.5,6]")

def add_output_group(parser, required=True):
    """output group"""

    output_group = parser.add_mutually_exclusive_group(required=required)
    output_group.add_argument("-o", "--ofile", dest="ofile", type=str, default='enrich_report',
                              help="Output file name. Mutually exclusive with --o-prefix.")
    output_group.add_argument("--o-prefix", dest="ofile", type=str, default='enrich_report',
                              help="Output file prefix. Mutually exclusive with -o/--ofile.")



def add_call_parser(subparsers):
    """Add main function 'call' argument parsers."""

    argparser_call = subparsers.add_parser("call", help="Main GSEAPY Function: run GSEAPY instead of GSEA.")

    # group for input files
    group_input = argparser_call.add_argument_group("Input files arguments")
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
    group_output = argparser_call.add_argument_group("Output arguments")
    add_output_option(group_output)

     # group for General options.
    group_opt = argparser_call.add_argument_group("GSEA advanced arguments")
    group_opt.add_argument("--min-size",  dest="mins", action="store", type=int, default=15, metavar='',
                           help="Min size of input genes presented in Gene Sets. Default: 15")
    group_opt.add_argument("--max-size", dest = "maxs", action="store", type=int, default =1000, metavar='',
                           help="Max size of input genes presented in Gene Sets. Default: 1000")
    group_opt.add_argument("-n", "--permu-num", dest = "n", action="store", type=int, default=1000, metavar='',
                           help="Number of random permutations. For calculating esnulls. Default: 1000")
    group_opt.add_argument("-w", "--weight", action='store', dest='weight', default=1, type=float, metavar='',
                           help='Weighted_score of rank_metrics.For weighting input genes. Choose from {0, 1, 1.5, 2},default: 1',)
    group_opt.add_argument("-m", "--method", action="store", dest="method", type=str, metavar='',
                           choices=("signal_to_noise", "t_test", "ratio_of_classes", "diff_of_classes", "log2_ratio_of_classes"),
                           default="log2_ratio_of_classes",
                           help="Methods to calculate correlations of ranking metrics. \
                           Choose from {'signal_to_noise', 't_test', 'ratio_of_classes', 'diff_of_classes','log2_ratio_of_classes'}.\
                           Default: 'log2_ratio_of_classes'")
    group_opt.add_argument("-a", "--ascending", action='store_true', dest='ascending', default=False,
                           help='Rank metric sorting order. If the -a flag was chosen, then ascending equals to True. Default: False.')
    group_opt.add_argument("-t", "--top-graphNum", dest = "graph", action="store", type=int, default=20, metavar='',
                           help="Plot graphs for top sets of each phenotype. Default: 20")
    group_opt.add_argument("-s", "--seed", dest = "seed", action="store", type=int, default=None, metavar='',
                           help="Number of random seed. Default: 2000")

    return
    
def add_prerank_parser(subparsers):
    """Add function 'prerank' argument parsers."""

    argparser_prerank = subparsers.add_parser("prerank", help="Run GSEAPY on preranked gene list.")

    # group for input files
    prerank_input = argparser_prerank.add_argument_group("Input files arguments")
    prerank_input.add_argument("-r", "--rnk", dest="rnk", action="store", type=str, required=True,  
                             help="ranking dataset file in .rnk format.Same with GSEA.")
    prerank_input.add_argument("-g", "--gmt", dest="gmt", action="store", type=str, required=True, 
                             help="Gene set database in GMT format. Same with GSEA.")
    prerank_input.add_argument("-l", "--label", action='store', nargs=2, dest='label',
                             metavar=('pos', 'neg'), type=float, default=['Postive','Negative'],
                             help="The phenotype label argument need two parameter to define. Default: ['Postive','Negative']")

    # group for output files
    prerank_output = argparser_prerank.add_argument_group("Output arguments")
    add_output_option(prerank_output)

     # group for General options.
    prerank_opt = argparser_prerank.add_argument_group("GSEA advanced arguments")
    prerank_opt.add_argument("--min-size",  dest="mins", action="store", type=int, default=15, metavar='',
                             help="Min size of input genes presented in Gene Sets. Default: 15")
    prerank_opt.add_argument("--max-size", dest = "maxs", action="store", type=int, default =1000, metavar='',
                             help="Max size of input genes presented in Gene Sets. Default: 1000")
    prerank_opt.add_argument("-n", "--permu-num", dest = "n", action="store", type=int, default=1000, metavar='',
                             help="Number of random permutations. For calculating esnulls. Default: 1000")
    prerank_opt.add_argument("-w", "--weight", action='store', dest='weight', default=1, type=float, metavar='',
                             help='Weighted_score of rank_metrics.For weighting input genes. Choose from {0, 1, 1.5, 2},default: 1',)
    prerank_opt.add_argument("-a", "--ascending", action='store_true', dest='ascending', default=False,
                             help='Rank metric sorting order. If the -a flag was chosen, then ascending equals to True. Default: False.')
    prerank_opt.add_argument("-t", "--top-graphNum", dest = "graph", action="store", type=int, default=20, metavar='',
                             help="Plot graphs for top sets of each phenotype. Default: 20")
    prerank_opt.add_argument("-s", "--seed", dest = "seed", action="store", type=int, default=None, metavar='',
                             help="Number of random seed. Default: 2000")
    
    return
    
def add_plot_parser(subparsers):
    """Add function 'plot' argument parsers."""

    argparser_replot = subparsers.add_parser("replot", help="Reproduce GSEA desktop figures.")

    group_replot = argparser_replot.add_argument_group("Input files arguments")

    group_replot.add_argument("-i", "--indir", action="store", dest="indir", required=True, metavar='',
                              help="The GSEA desktop results directroy that you want to reproduce the figure ")
    add_output_option(group_replot)
    #add_output_group( argparser_plot )
    group_replot.add_argument("-w", "--weight", action='store', dest='weight', default=1, type=float, metavar='',
                              help='Weighted_score of rank_metrics. Please Use the same value in GSEA. Choose from (0, 1, 1.5, 2),default: 1',)

    return
def add_enrichr_parser(subparsers):
    """Add function 'enrichr' argument parsers."""
    
    argparser_enrichr = subparsers.add_parser("enrichr", help="Peform GSEA using enrichr API.")

    group_enrichr = argparser_enrichr.add_argument_group("Input files arguments")

    group_enrichr.add_argument("-i", "--input-list", action="store", dest="gene_list", required=True, metavar='',
                              help="Enrichr uses a list of Entrez gene symbols as input.  ")
    group_enrichr.add_argument("-d", "--description", action='store', dest='description', default='foo',metavar='',
                              help="It is recommended to enter a description for your list so that" +
                                   " multiple lists can be differentiated from each other if you choose to save or share your list")
    group_enrichr.add_argument("-g", "--gene-sets", action="store", dest="library", required=True, metavar='',
                              help="Enrichr library name requires. see online tool for names.  ")

    add_output_group(group_enrichr)

    return

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted me! ;-) Bye!\n")
        sys.exit(0)
