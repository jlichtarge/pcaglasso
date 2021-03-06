from cmd import Cmd


class Prompt(Cmd):
    current_type_combo = []
    
    network_help_summary = 'Generates glasso network, saves image as .pdf and network as .gml'
    network_help_usage = 'network [-no_lab] [-cc] [-gl_only | [[-color_pca] [-size_pca]] ]'
    network_help_body = \
'Network Display Options:\n\
  Network-specific options:\n\
    -no_lab\n\
        Turn off node labeling.\n\
    -cc\n\
        Only show connected components of graph; nodes without edges are not displayed.\n\
    -gl_only\n\
        Change network draw style to highlight correlative structure:\n\
            Green lines are positive GL correlation, red are negative.\n\
        Default:\n\
            Solid lines are positive GL correlation, dashed are negative.\n\
        Note: if -gl_only is used, -color_pca and -size_pca will be ignored.\n\
  PCA-pair options:\n\
    Note:\n\
      These options visualize PCA data generated by running PCA on pairwise combinations\n\
      of sample types. To visualize a particular pair of types for PCA analysis, the PCA data \n\
      must first be initialized by running the \'pca_pair X Y\' command. When the \'pca_pair\'\n\
      command is run, the pairwise PCA data visualized by the \'network\' command is overwritten\n\
      with the data of the new run\'s type pair. If pca_pair has not yet been run, the data is\n\
      initialized to a PCA analysis of the first two types: [0,1].\n\
    -color_pca\n\
        Set node color as PCA-associated type (Ex: all green nodes same type).\n\
        Set color intensity as PCA score magnitude (Ex: higher PCA score => more vibrant color).\n\
    -size_pca\n\
        Set node size to be change in concentration over last PCA pairing.\n\
            A large node size indicates a larger magnitude change in value/concentration for\n\
            that feature between types X and Y (X and Y being the two types most recently run\n\
            via the \'pca_pair\' command). A small node size indicates a smaller magnitude\n\
            change in concentration.\n\
            Nodes with no change will have a marginal node size for visualization purposes.\n\
        Set node icon to be arrow showing direction of change in concentration.\n\
            A node represented by an upward-facing arrow has a higher concentration in the type\n\
            represented by the larger of the type-designating numbers. A node with a\n\
            downward-facing arrow has a higher concentration in the type with the lower number.\n\
            Ex: If the PCA data being visualized is generated by comparison of types X and Y,\n\
                where Y > X, then a downward arrow representing a node indicates a higher\n\
                concentration for that feature in type X, and an upward arrow indicates a higher\n\
                concentration in type Y. The numbers associated with each type are assigned in\n\
                order of the each type\'s precedence in the input data.'
    network_args = ['-no_lab', '-cc', '-gl_only', '-size_pca', '-color_pca']

    def help_network(self):
        print 'Summary:', self.network_help_summary
        print 'Usage:', self.network_help_usage
        print self.network_help_body
        print 'Current PCA type pair:', self.gp.gl.current_type_combo
     
    def do_network(self, args):
        """Generates glasso network"""
        args = args.split()
        print 'Parsing args:', args
        
        args = self._filter_network_args(args, self.network_args)
        
        no_labels = bool('-no_lab' in args)
        only_conn_comp = bool('-cc' in args)
        glasso_only = bool('-gl_only' in args)  
        nsize = 'delta_conc_scaled' if '-size_pca' in args and '-gl_only' not in args \
            else 'default'
        ncolor = 'pca_type' if '-color_pca' in args and '-gl_only' not in args \
            else 'default'
        
        if nsize != 'default' or ncolor != 'default':
            print 'Current PCA type pair:', self.gp.gl.current_type_combo
            print 'To use a different pair of types [X,Y], run \'pca_pair X Y\'.'
                        
        print 'Generating network...'
        self.gp.gl.network_plot(glasso_only=glasso_only,
                                node_color_selector=ncolor,
                                node_size_selector=nsize,
                                draw_labels=not no_labels,
                                only_conn_comp=only_conn_comp,
                                tag='')
        print '...done'
    
    def help_pca_pair(self):
        print 'Summary: Runs PCA on selected pair of types'
        print '\tGenerates pairwise PCA analysis graphs, and PCA-overlaid glasso network'
        print '\tInitializes PCA data for further graph generation via \'network\' command.'
        print 'Usage: pca_pair {choose 2 of (0 1 ... number_of_sample_types)}'
        print '\tEx: \'pca_pair 0 1\' runs PCA on pair of types 0 and 1'
        print '\tEx: \'pca_pair 4 2\' runs PCA on pair of types 2 and 4'

    def do_pca_pair(self, args):
        """Runs PCA on selected type pair"""
        print len(self.gp.pca.types)
        if len(args) != 3:
            print 'Wrong number of args:'
            print '\tWanted two types for pairwise PCA analysis, got ',len(args)-1,' args.'
            self.help_pca_pair()
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
                        \n\tmust be selected from available types: ', range(len(self.gp.pca.types))
                return
            print 'Parsing selected types:', args, '->', select_types
        
        print 'Running PCA...'
        self.gp.pca.run_pca(select_types=select_types)
        print '...done'
        
        print 'Generating Graphs...'
        pca_folder = 'pca_plots/' + \
            self.gp.pca.types[select_types[0]]['name'] + \
            '_' + self.gp.pca.types[select_types[1]]['name'] + '/'
        self.gp.pca.comps_plot_1d(folder=pca_folder)
        self.gp.pca.comps_explained_var_plot(folder=pca_folder, tag='')
        self.gp.gl.load_pca_data(self.gp.pca.feature_scores_dict,
                                 self.gp.pca.run_types_by_num,
                                 pca_component=0)
        self.gp.gl.network_plot(glasso_only=False,
                                node_color_selector='pca_type',
                                node_size_selector='delta_conc_scaled',
                                draw_labels=True,
                                only_conn_comp=False,
                                folder=pca_folder,
                                tag='')
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
        
#         Initialize PCA pairwise type data with first two types. 
        self.gp.pca.run_pca(select_types=[0,1])
        self.gp.gl.load_pca_data(self.gp.pca.feature_scores_dict,
                                 self.gp.pca.run_types_by_num,
                                 pca_component=0)
        print 'Select from the below options which graphs to generate:'
        self.do_help('')
        
    def _filter_network_args(self, got, want):
        bad = [b for b in got if b not in want]
        good = list(set([g for g in got if g in want]))
        if len(bad) != 0:
            print 'Warning: unrecognized args \'',bad,'\'will be ignored'
            print 'Using args:',good
        return good

    def __init__(self,
                 glasso_pca_obj,
                 completekey='tab',
                 stdin=None,
                 stdout=None):
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
