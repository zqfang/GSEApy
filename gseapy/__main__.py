

import os
import sys
import argparse as ap


# ------------------------------------
# Main function
# ------------------------------------

__version__ = '0.3.0'

def main():
    """The Main function/pipeline for GSEAPY.
    
    """
    # Parse options...
    argparser = prepare_argparser()
    args = argparser.parse_args()

    if args.outdir:
        # use a output directory to store GSEAPY output
        if not os.path.exists( args.outdir ):
            try:
                os.makedirs( args.outdir )
            except:
                sys.exit( "Output directory (%s) could not be created. Terminating program." % args.outdir )

    subcommand  = args.subcommand_name

    if subcommand == "replot":
        # reproduce plots using GSEAPY
        from .gsea import replot               

        replot(indir=args.indir,outdir=args.outdir,weight=args.weight,figsize=args.figsize,format=args.format,) 
    
    elif subcommand == "call":
        # compute using GSEAPY
        from .gsea import run
                
        run(args.data, args.gmt,args.cls, args.mins, args.maxs, args.n, args.weight,
            args.type, args.method,args.ascending, args.outdir,args.figsize,args.format,)


def prepare_argparser ():
    """Prepare argparser object. New options will be added in this
    function first.
    
    """
    description = "%(prog)s -- GSEAPY for Gene Set Enrichment Analysis"
    epilog = "For command line options of each command, type: %(prog)s COMMAND -h"

    # top-level parser
    argparser = ap.ArgumentParser( description = description, epilog = epilog )
    argparser.add_argument("--version", action="version", version="%(prog)s "+ __version__)
    subparsers = argparser.add_subparsers( dest = 'subcommand_name' ) #help="sub-command help")
    
    # command for 'call'
    add_call_parser( subparsers )

    
    # command for 'plot'
    add_plot_parser( subparsers )
    
    
    return argparser

def add_output_option ( parser ):
    parser.add_argument("-o", "--outdir", dest = "outdir", type = str, default = 'gseapy_out',
                        metavar='',action="store",
                        help = "The gseapy output directory. Default: the current working directory")

    parser.add_argument( "-f", "--format", dest = "format", type = str, metavar='',action="store",
                              choices = ("pdf", "png", "jpeg", "eps"),default = "pdf",
                              help = "Format of output figures, choose from {'pdf', 'png', 'jpeg', 'eps'}. Default: 'pdf'." ) 
    parser.add_argument("--figsize",action='store',nargs=2,dest='figsize',
                        metavar=('width', 'height'),type=float,default=[6.5,6],
                        help="The figsize keyword argument need two parameter to define. Default: [6.5,6]") 

def add_output_group ( parser, required = True ):
    output_group = parser.add_mutually_exclusive_group( required = required )
    output_group.add_argument( "--ofile", dest = "ofile", type = str,
                               help = "Output file name. Mutually exclusive with --o-prefix." )
    output_group.add_argument( "--o-prefix", dest = "oprefix", type = str,
                               help = "Output file prefix. Mutually exclusive with -o/--ofile." )



def add_call_parser( subparsers ):
    """Add main function 'call' argument parsers.
    """
    argparser_call = subparsers.add_parser("call", help="Main GSEAPY Function: run GSEAPY instead of GSEA.")
    
    # group for input files
    group_input = argparser_call.add_argument_group( "Input files arguments" )
    group_input.add_argument( "-d", "--datab", dest = "data",  action="store",type = str, required = True, 
                              help = "Expression table of phenotypes. Expected a txt file.Same with GSEA." )
    group_input.add_argument( "-c", "--cls", dest = "cls",  action="store", type = str, required = True,
                                    help = ".cls files. Same with GSEA.")
    group_input.add_argument( "-g", "--gmt", dest = "gmt",  action="store", type = str, required = True,
                              help = "Gene Sets in gmt format. Same with GSEA." )
    group_input.add_argument( "-p", "--permu-type",  action="store",dest = "type", type = str,metavar='',
                              choices = ("gene_set", "phnotype"),default = "gene_set",
                              help = "Gene Sets in gmt format. Same with GSEA, choose from {'gene_set', 'phenotype'}")

    # group for output files
    group_output = argparser_call.add_argument_group( "Output arguments" )
    add_output_option( group_output )    
       
     # group for General options.
    group_opt = argparser_call.add_argument_group( "GSEA advance arguments" )
    group_opt.add_argument( "--min-size",  dest = "mins",  action="store",type = int, default =15,metavar='',
                            help = "Min size of gene sets. Default: 15")
    group_opt.add_argument( "--max-size", dest = "maxs",  action="store",type = int, default = 1000,metavar='',
                            help = "Max size of gene sets. Default: 1000")
    group_opt.add_argument( "-n", "--permu-num", dest = "n",  action="store",type = int, default = 1000, metavar='',
                            help = "Permutation number. Default: 1000" )
    group_opt.add_argument("-w","--weight",action='store',dest='weight',default= 1, type= float,metavar='',
                            help='Weighted_score type of rank_metrics.Choose from {0, 1, 1.5, 2},default: 1',)

    group_opt.add_argument( "-m", "--method",  action="store",dest = "method", type = str, metavar='',
                            choices = ("signal_to_noise", "t_test", "ratio_of_classes", "diff_of_classes","log2_ratio_of_classes"),
                            default = "log2_ratio_of_classes",
                            help = "Methods to calculate correlations of ranking metrics. \
                            Choose from {'signal_to_noise', 't_test', 'ratio_of_classes', 'diff_of_classes','log2_ratio_of_classes'}.\
                            Default: 'log2_ratio_of_classes'" )   
    group_opt.add_argument("-a","--ascending",action='store_true',dest='ascending',default= False ,
                            help='Rank metrice sorting order. If the -a flag was chosen, then ascending equals to True. Default: False.')

    
     
    return


def add_plot_parser( subparsers ):
    """Add function 'plot' argument parsers.
    """    
    argparser_replot = subparsers.add_parser( "replot",help = "Reproduce GSEA desktop figures." )
    
    group_replot = argparser_replot.add_argument_group( "Positional arguments" )

    group_replot.add_argument("-i","--indir", action="store", dest="indir", required=True, metavar='',
                        help="The GSEA desktop results directroy that you want to reproduce the figure ")
    add_output_option( group_replot)
    #add_output_group( argparser_plot )
    group_replot.add_argument("-w","--weight",action='store',dest='weight',default= 1, type= float,metavar='',
                        help='Weighted_score type of rank_metrics.Choose from (0, 1, 1.5, 2),default: 1',)
       
    return




if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted me! ;-) Bye!\n")
        sys.exit(0)