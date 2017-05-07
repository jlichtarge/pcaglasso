from cmd import Cmd

class Prompt(Cmd):
       
    def help_network(self):
        print 'Generates glasso network'
        print 'usage: \'network\' [options]'
        print 'options:\
        \n\t-gl_only\
        \n\t\tchange network draw style to focus on graph structure:\
        \n\t\tgreen lines are positive GL correlation, red are negative\
        \n\t\tdefault: solid lines are positive GL correlation, dashed are negative\
        \n\t-no_labels\
        \n\t\tturn off node labeling\
        \n\t-only_connected_nodes\
        \n\t\tonly show major connected component of graph'
#         \n\t-nc_pca\
#         \n\t\tset node color: PCA score for intensity, PCA-associated type for color\
#         \n\t-ns_delta_conc\
#         \n\t\tset node size to be change in concentration over last PCA pairing\
#         \n\t\tdefault: constant'
     
    def do_network(self, args):
        """Generates glasso network"""
        args = args.split()
        print 'Parsing args:', args
        glasso_only = True if '-gl_only' in args else False        
        no_labels = True if '-no_labels' in args else False   
        only_connected_nodes = True if '-only_connected_nodes' in args else False
        ncolor = 'delta_conc_scaled' if '-ns_delta_conc' in args else 'default' 
        nsize = 'pca_type' if '-nc_pca' in args else 'default' 

        print 'Generating network...'
        self.gp.gl.network_plot(glasso_only = glasso_only,
                                node_color_selector    = ncolor,
                                node_size_selector     = nsize,
                                draw_labels            = not no_labels,
                                show_disconnected_nodes= not only_connected_nodes,
                                tag                    = '')
        print '...done'
    
    def help_pca_pair(self):
        print "Runs PCA on selected type pair"
        print "Generates PCA type-specific graphs, and PCA-overlaid glasso network"
        print 'selected types must be space-separated ints (ex: \'pca_pair 0 1\')'
        print 'selection options',range(len(self.gp.pca.types))
        print 'usage: pca_pair [0 1 2 3 4] --> choose 2'
        print 'example: \'pca_pair 0 1\''

    def do_pca_pair(self, args):
        """Runs PCA on selected type pair"""
        print len(self.gp.pca.types)
        if len(args) == 0:
            print "wrong number of args, try again"
            print "example: \'run_pca 0 1\'"
            return
        else:
            args = args.split()
            try:
                select_types = [int(x) for x in args]
                select_types = list(set(select_types))
                select_types.sort()
            except ValueError:
                print 'selected types must be space-separated ints (ex: \'run_pca 0 1\')'
                return
            if (select_types[0] < 0) or (select_types[-1] > len(self.gp.pca.types)):
                print 'selected types out of range...\
                        \n\tmust be selected from available types: ',range(len(self.gp.pca.types))
                return
            print 'Parsing selected types:', args,'->',select_types
        
        print 'Running PCA...'
        self.gp.pca.run_pca(select_types = select_types)
        print '...done'
        
        print 'Generating Graphs...'
        pca_folder = 'pca_plots/'+self.gp.pca.types[select_types[0]]['name']+'_'+self.gp.pca.types[select_types[1]]['name']+'/'
        self.gp.pca.comps_plot_1d(folder=pca_folder)
        self.gp.pca.comps_explained_var_plot(folder=pca_folder, tag='')
        self.gp.gl.load_pca_data(self.gp.pca.feature_scores_dict, self.gp.pca.run_types_by_num, pca_component = 0)
        self.gp.gl.network_plot(glasso_only = False,
                                node_color_selector    = 'pca_type',
                                node_size_selector     = 'delta_conc_scaled',
                                draw_labels            = True,
                                show_disconnected_nodes= True,
                                folder                 = pca_folder,
                                tag                    = '')
        print '...done'

    def do_q(self, args):
        """Quits the program."""
        print "Quitting."
        raise SystemExit
    
    def do_quit(self, args):
        """Quits the program."""
        print "Quitting."
        raise SystemExit
    
    def preloop(self):
        '''Explains usage of interactive console
        '''
        print 'Input data has been successfully loaded.'
        print 'Select from the below options which graphs to generate:'
        self.do_help('')
        
    def __init__(self, glasso_pca_obj, completekey='tab', stdin=None, stdout=None):
        """Instantiate a line-oriented interpreter framework.

        The optional argument 'completekey' is the readline name of a
        completion key; it defaults to the Tab key. If completekey is
        not None and the readline module is available, command completion
        is done automatically. The optional arguments stdin and stdout
        specify alternate input and output file objects; if not specified,
        sys.stdin and sys.stdout are used.

        """
        import sys
        if stdin is not None:
            self.stdin = stdin
        else:
            self.stdin = sys.stdin
        if stdout is not None:
            self.stdout = stdout
        else:
            self.stdout = sys.stdout
        self.cmdqueue = []
        self.completekey = completekey
        self.gp = glasso_pca_obj