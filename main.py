from glasso_pca import GLASSO_PCA
import interactive_flow as flow
import os.path
import argparse
'''
Wednesday March 15 2017
TODO:
finish up with interactive console
 - especially: running pca with types, overlaying PCA data on glasso
investigate 'without types' functionality
clean up all functions in use
remove non-used functions

'''


def main():
    args = parseArgs()
    
    gl_alpha = args.glasso_alpha if args.glasso_alpha else 0.58 #TODO: -1
    pca_comp_var_threshold = args.pca_threshold if args.pca_threshold else 0.1
    
    gp = GLASSO_PCA(inputfile = args.filename,
                   verbose = not args.quiet,
                   gl_alpha = gl_alpha,
                   pca_comp_var_threshold = pca_comp_var_threshold)
    gp.load()
    if not args.no_glasso_plots:
        gp.glasso_plots()
    
    if not args.no_pca_plots:
        gp.pca_plots()
    
    if not args.no_boxplots:
        gp.boxplots()
    
    p = flow.Prompt(gp)
    p.prompt = '> '
    p.cmdloop('Starting prompt:')

def parseArgs():
    '''Parses input arguments. Run 'python main.py -h' for descriptions
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="the input file to be analyzed")
    parser.add_argument("--quiet", 
                        help="Toggle stdout printouts through analysis, calculation. Default = verbose",
                        action="store_true")
    parser.add_argument("--glasso_alpha", 
                        help="Alpha parameter for GLASSO. Must be >= 0. Default: selected by cross validation",
                        type=float)
    parser.add_argument("--pca_threshold", 
                        help="Threshold of explained variance for PCA components. Must be [0,1]. Only components above threshold will be analyzed. Default = 0.1",
                        type=float)
    parser.add_argument("--no_boxplots", 
                        help="Do not generate boxplots on load. Creating boxplots will be an option from interactive console.",
                        action="store_true")
    parser.add_argument("--no_pca_plots", 
                        help="Do not generate pca plots on load. Creating pca plots will be an option from interactive console.",
                        action="store_true")
    parser.add_argument("--no_glasso_plots", 
                        help="Do not generate glasso plots on load. Creating glasso plots will be an option from interactive console.",
                        action="store_true")
    args = parser.parse_args()
    
    if not os.path.isfile(args.filename):
        print 'invalid filename... see \'-h\' for help'
        raise SystemExit
    
    if args.glasso_alpha and args.glasso_alpha < 0:
        print 'glasso_alpha must be >= 0'
        raise SystemExit
    
    if args.pca_threshold and not(0 <= args.pca_threshold <= 1):
        print 'pca_threshold must be within interval [0,1]'
        raise SystemExit
    return args

    
if __name__ == '__main__':
    main()
