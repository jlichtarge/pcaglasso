'''
Created on Jun 11, 2016
TODO:
network graph:
    add option to remove unconnected notes from network graph
boxplots:
    lots of stuff for aesthetics

consider:
    function that combines boxplots/network graph
    
    
12/16/16:
TODO:
graph structure from glasso on all data G = {V,E} (edges same across all k)
Gk = HESC, EB, Rosette, NPC, Neuron (V varies with k)
Vik = Value of node i in cell k
size of Vik = pca coefficient
color of Vik = magnitude of concentration


Generate graphs G1/2/3/4/5 with above constraints

Figure 3:
Two cell types (1/2):
graph of PCA separation
boxplots of major metabolites
G1 and G2
difference between g1 and g2
--clearly identify subnetworks

Figure 4:
Composite of all transitions, subnetworks clearly identified

Writeup:
generate graphs
assemble into figures
subtext explanation for each figure

1/16 changes:
Vik:
size  = magnitude of change in concentration
color = pca coefficient
shape = arrow pointing in direction of change of concentration over transition between cell types
implemented glasso cross validation

TODO:
sparse pca?
show non-connected nodes alongside network graph
refactor graph functions
boxplots
make different set of graphs, one for each cell type, in which size = absolute concentration, and color = sum of PCA scores?

to discuss:
how to show difference between g(1/2) and g(2/3)?
is sparse PCA a priority?
what is the most important information for biological interpretation?
-pca score
-change in concentration

-absolute concentration
-change in pca score?

what functionalities need to be exposed to user in software package?


meeting notes
switch solid/dashed lines
command line based python script
allows pca/graphical lasso to generate network plots
output network structure into common format (hashlist, GML)
output graph in pdf
output dump in txt file

mock the figures

explore glasso - pca relationshipn w positive/negative correlations


2/12/17 - Pre-meeting notes
TODO:
refactor box plots
-line showing 0 value
-choose cutoff, specify number, or all
-toggle pca color/type color
-tile boxplots or separate files? consider new folder if separate files
-toggle raw data, scaled data
-standardize y scale (at least within a batch)

save network structure for later use


node color from 
@author: Jared
'''

from input_output_functions import dataio
from pca_functions import pca_bundle
from glasso_functions import glasso_bundle
import sys

class GLASSO_PCA(object):
    '''orchestrates all analytical functionality of the package
    '''
    
    # these vars can be set for programmatic use (rather than cmd-line interactive use) 
    INPUTFILE = 'all_clean.txt'
    PCA_COMP_VAR_THRESHOLD = 0.1
    GLASSO_ALPHA = -1
    
    
    def load(self):
        print 'starting glasso_pca load ....'
        self.io = dataio(self.inputfile, use_types=True, verbose=self.verbose)
    
        self.pca = pca_bundle(self.io.savedir,
                            self.io.data,
                            self.io.rawdata,
                            self.io.feature_names,
                            self.io.sample_names,
                            self.io.types,
                            comp_var_threshold=self.pca_comp_var_threshold,
                            verbose=self.verbose)
        
        self.gl = glasso_bundle(self.io.savedir,
                                self.io.data,
                                self.io.rawdata,
                                self.io.types,
                                self.io.feature_names,
                                alpha=self.gl_alpha,
                                verbose=self.verbose)
#         gl.network_plot(draw_labels=True, glasso_only = True, node_size_selector ='default')
#     #     gl.network_plot(draw_labels=True, glasso_only = True, node_size_selector ='network_edge_weight')
#     #     gl.network_plot(draw_labels=True, glasso_only = True, node_size_selector ='delta_conc_scaled')
#     #     gl.network_plot(draw_labels=True, glasso_only = True, node_size_selector ='delta_conc_raw')
#     #     gl.network_node_attributes_plot()
#     
#         pca.run_pca(select_types = [0,1,2,3,4])
#     
#     #     pca.boxplot(0, number_boxplots = 4, tag = '_all')
#         pca.comps_plot_1d(tag = '_all')
#         pca.comps_plot_3d(tag = '_all', projections = False, show_legend = False)
#         pca.comps_explained_var_plot(tag = '')
#     #     pca.gen_boxplots()
#         print '~done~'
#         sys.exit()
#     
#         feature_data_list = []
#         for typ in range(len(io.types)-1):
#             type_combo = [typ,typ+1]
#             tag = ''.join([str(x) for x in type_combo])+'_'
#             print '\n\n+++++++++++++++++++++++++++'
#             #=======================================================================
#             #=======TODO:================================================================
#             # # http://stackoverflow.com/questions/16834861/create-own-colormap-using-matplotlib-and-plot-color-scale
#             #=======================================================================
#             #=======================================================================
#             pca.run_pca(select_types = type_combo)
#             
#     #         pca.boxplot(0, number_boxplots = 4, tag = tag)
#     #         pca.boxplot_single(0, number_boxplots = 3, tag = tag)
#     #         pca.comps_explained_var_plot(tag = tag)
#     #         pca.comps_plot_1d(components = [0], tag = tag)
#     #         feature_data_dict, fd, cm = pca.loadings_plot_1d(component_num = 0, tag = '', no_figure = False, cutoff_percentile = 0)
#     #         feature_data_list.append(fd)
#     #         _, nodelist = gl.network_plot1(component = 0, colordict = feature_data_dict, pos = pos, tag = tag)
#             gl.load_pca_data(pca.feature_scores_dict, pca.run_types_by_num, pca_component = 0)
#             gl.network_plot(node_color_selector = 'pca_type',
#                             node_size_selector = 'delta_conc_scaled',
#                             draw_labels = False)
#             gl.network_plot(node_color_selector = 'pca_raw',
#                             node_size_selector = 'delta_conc_raw',
#                             draw_labels = False)
    
    def __init__(self,
                 inputfile=INPUTFILE,
                 pca_comp_var_threshold=PCA_COMP_VAR_THRESHOLD,
                 gl_alpha=GLASSO_ALPHA,
                 verbose=True):
        
        self.inputfile = inputfile
        self.pca_comp_var_threshold = pca_comp_var_threshold
        self.gl_alpha = gl_alpha
        self.verbose = verbose
    
        if self.verbose:
            print 'init glasso_pca:'
            print '\tinput file:', inputfile
            print '\tgl alpha:', gl_alpha
            print '\tpca threshold:', pca_comp_var_threshold
    
    def quick(self):
        print "Running PCA on sequential type pairs:"
        for x in range(len(self.pca.types) - 1):
            select_types = [x, x + 1]
            print 'Running PCA...'
            self.pca.run_pca(select_types=select_types)
            print '...done'
            
            print 'Generating Graphs...'
            pca_folder = 'pca_plots/' + self.pca.types[select_types[0]]['name'] + '_' + self.pca.types[select_types[1]]['name'] + '/'
            self.pca.comps_plot_1d(folder=pca_folder)
            self.pca.comps_explained_var_plot(folder=pca_folder, tag='')
            self.gl.load_pca_data(self.pca.feature_scores_dict, self.pca.run_types_by_num, pca_component=0)
            self.gl.network_plot(glasso_only=False,
                                    node_color_selector='pca_type',
                                    node_size_selector='delta_conc_scaled',
                                    draw_labels=True,
                                    show_disconnected_nodes=True,
                                    folder=pca_folder,
                                    tag='')
            print '...done'
            print"-----------------------"

        
    def glasso_plots(self):
        glasso_folder = 'glasso_plots/'
        self.gl.network_plot(draw_labels=True,
                             glasso_only=True,
                             node_size_selector='default',
                             folder=glasso_folder)
        self.gl.network_node_attributes_plot(folder=glasso_folder)
        self.gl.network_node_attributes_file(folder=glasso_folder)
        self.gl._save_network_graph(folder=glasso_folder)

    def pca_plots(self, select_types=None):
        pca_folder = 'pca_plots/'
        if select_types is None:
            select_types = range(len(self.pca.types))
            
        self.pca.run_pca(select_types=select_types)
        self.pca.comps_plot_1d(folder=pca_folder)
        self.pca.comps_plot_3d(folder=pca_folder,
                               tag='_all',
                               projections=False)
        self.pca.comps_explained_var_plot(folder=pca_folder, tag='')
        
    def boxplots(self):
        self.pca.gen_boxplots()

    
#     def old_code(self):
        # sorted(zip(features, loadvals, color_range), key = lambda x: x[1])
    #     for j, fd in enumerate(feature_data_list):
    #         if j %2 ==0:
    #             feats = [x[0] for x in fd]
    #             loads = [x[1] for x in fd]
    #             cm = plt.get_cmap('RdBu')
    #             cols = [cm((-x)+0.5)[:-1] for x in loads]# '+0.5' centers the colormap around white
    #             fd = sorted(zip(feats, loads, cols), key = lambda x: x[1])
    #         fdd = {}
    #         for i in range(len(fd)):
    #             fdd[fd[i][0]] = fd[i][1:]
    #                 
    #         _, nodelist = gl.network_plot(component = 0, colordict = fdd, pos = pos, tag ='_test'+str(j)+str(j+1))
    
    #     pca.run_pca(max_comps=PCA_MAX_COMPS, comp_var_threshold = PCA_COMP_VAR_THRESHOLD)
    #     pca.boxplot(0, number_boxplots = 6, tag = tag)
    
    
    #         gl.network_node_attributes_plot(nodelist, colordict = feature_data_dict, tag = 'color'+tag)
    #         gl.network_node_attributes_plot(nodelist, tag = 'none'+tag)
    #         gl.network_node_attributes_file(nodelist, tag = tag)
    
    #     pca.loadings
    #     pca.scores
    #     pca.good_comps
    #     pca.comp_var
    
    
    #=====OLD PCA CODE==========================================================================
    #     for i in range(len(pca.good_comps)):
    #         pca.boxplot(i, number_boxplots = 3)
    # 
    #     #===========================================================================
    #     pca.comps_explained_var_plot()
    #     pca.comps_plot_1d()
    #     pca.pca_comps_plot_2d(io.sample_names)
    #     pca.comps_plot_3d(io.sample_names)
    #     #===========================================================================
    #     feature_data_dict = pca.loadings_plot_1d(io.feature_names, component_num = 0)
    #     feature_data_dict1 = pca.loadings_plot_1d(io.feature_names, component_num = 1)
    #     feature_data_dict2 = pca.loadings_plot_1d(io.feature_names, component_num = 2)
    #     feature_data_dict3 = pca.loadings_plot_1d(io.feature_names, component_num = 3)
    #===============================================================================
    #========OLD GL CODE=======================================================================
    #     pos1, nodelist = gl.network_plot(component = 0, colordict = feature_data_dict, pos = pos1)
    # 
    #     
    #     gl = glasso_bundle(io.savedir, io.data, io.feature_names, GLASSO_ALPHA, GLASSO_MAX_ITERATIONS)
    #     pos1, nodelist = gl.network_plot(component = 0, colordict = feature_data_dict)
    #     pos1, nodelist = gl.network_plot(component = 0, colordict = feature_data_dict, pos = pos1, drawlbls = False, tag = 'nolbl')
    # #     pos1, nodelist = gl.network_plot(component = 1, colordict = feature_data_dict)
    #     pos1, nodelist = gl.network_plot(component = 1, colordict = feature_data_dict1, pos = pos1, drawlbls = False, tag = 'nolbl')
    #     pos1, nodelist = gl.network_plot(component = 2, colordict = feature_data_dict2, pos = pos1, drawlbls = False, tag = 'nolbl')
    #     pos1, nodelist = gl.network_plot(component = 3, colordict = feature_data_dict3, pos = pos1, drawlbls = False, tag = 'nolbl')
    # 
    # 
    # #     pos, nodelist = gl.network_plot(pos = pos1, component = 1, colordict = feature_data_dict1)
    # #     pos, nodelist = gl.network_plot(pos = pos1, component = 2, colordict = feature_data_dict2)
    # #     pos, nodelist = gl.network_plot(pos = pos1, component = 3, colordict = feature_data_dict3)
    # #     pos, nodelist = gl.network_plot(io.feature_names, colordict = feature_data_dict)
    # 
    #     #===========================================================================
    #     gl.network_node_attributes_plot(nodelist, colordict = feature_data_dict, tag = 'color')
    #     gl.network_node_attributes_plot(nodelist, tag = 'none')
    #     gl.network_node_attributes_file(nodelist)
    #     #===========================================================================
    #===============================================================================
        return
