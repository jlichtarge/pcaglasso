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

boxplots
5/7/17
Still TODO:
Clean / finish interactive component
Go through file by file and clean up code
investigate 'without types' functionality
finish up with interactive console
 - especially: running pca with types, overlaying PCA data on glasso
'''


def main():
    args = parseArgs()
    
    gl_alpha = args.glalpha if args.glalpha else -1
    pca_comp_var_threshold = args.pca_threshold if args.pca_threshold else 0.1
    
    gp = GLASSO_PCA(inputfile=args.filename,
                   verbose=not(args.q or args.quiet),
                   gl_alpha=gl_alpha,
                   pca_comp_var_threshold=pca_comp_var_threshold)
    gp.load()
    if not args.no_plots:
        if not args.no_glasso_plots: gp.glasso_plots()
        if not args.no_pca_plots: gp.pca_plots()
        if not args.no_boxplots: gp.boxplots()
            
    if args.quick:
        gp.quick()
    
    if args.i or args.interact:
        p = flow.Prompt(gp)
        p.prompt = '> '
        p.cmdloop('Starting prompt:')

def parseArgs():
    '''Parses input arguments. Run 'python main.py -h' for descriptions
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='the input file to be analyzed')
    parser.add_argument('-quiet',
                        help='Suppress stdout printouts through analysis, calculation. Default = verbose',
                        action='store_true')
    parser.add_argument('-i','-interact',
                        help='Following load, enter interactive console to generate graphs with more finely-tuned parameters.',
                        action='store_true')
    parser.add_argument('-quick',
                        help='Run PCA over sequential pairs of input types, generating PCA-glasso plots for each.\
                        \nFor 5 input types, pairs will be (0,1) (1,2) (2,3) (3,4)',
                        action='store_true')
    parser.add_argument('-glalpha',
                        help='Set alpha parameter for GLASSO. Must be >= 0. Default: selected by cross validation',
                        type=float)
    parser.add_argument('-pca_threshold',
                        help='Threshold of explained variance for PCA components. Must be [0,1]. Only components above threshold will be analyzed. Default: 0.1',
                        type=float)
    parser.add_argument('-no_plots',
                        help='Do not generate any plots on load. Creating plots will be an option from interactive console.',
                        action='store_true')
    parser.add_argument('-no_boxplots',
                        help='Do not generate boxplots on load. Creating boxplots will be an option from interactive console.',
                        action='store_true')
    parser.add_argument('-no_pca_plots',
                        help='Do not generate pca plots on load. Creating pca plots will be an option from interactive console.',
                        action='store_true')
    parser.add_argument('-no_glasso_plots',
                        help='Do not generate glasso plots on load. Creating glasso plots will be an option from interactive console.',
                        action='store_true')
    args = parser.parse_args()
    
    if not os.path.isfile(args.filename):
        print 'invalid filename... see \'-h\' for help'
        raise SystemExit
    
    if args.glalpha and args.glalpha < 0:
        print 'glalpha must be >= 0\nDefaulting to assigning alpha by cross validation'
        args.glalpha = -1
    
    if args.pca_threshold and not(0 <= args.pca_threshold <= 1):
        print 'pca_threshold must be within interval [0,1]\nDefaulting to 0.1'
        args.pca_threshold = 0.1
    return args

    
if __name__ == '__main__':
    main()
