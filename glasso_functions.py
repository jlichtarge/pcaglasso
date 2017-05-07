'''
Created on Jun 11, 2016

@author: Jared
'''
import numpy as np
from sklearn.covariance import GraphLasso, GraphLassoCV
import matplotlib.pyplot as plt
import networkx as nx
from math import sqrt, exp, sin, pi
from math import e as math_e
from os import path, mkdir

class glasso_bundle(object):
    '''
    TODO: classdocs
    '''
    def _count_nz_offd(self,data):
        """Count off-diagonal non-zero elements of a matrix
        """
        data = data - np.diag(np.diag(data))
        count = 0
        for n in data.flatten():
            if n != 0:
                count += 1
        return count 

    def distance_calc(self,matrix):
        """    calculates distance value of n-dimensional coordinate system
            example:
                distance_calc(1.5) = 1.5
                distance_calc([1,2] = sqrt(1^2+2^2) = sqrt(3)
                distance_calc([1,2,5] = sqrt(1^2+2^2+5^2) = sqrt(30)
        """
        dist = 0
        npmatrix = np.array(matrix)
        if npmatrix.size >1:
            for i in range(npmatrix.size):
                dist += npmatrix[i]**2
            return sqrt(dist)
        else: 
            return np.absolute(matrix) 

#     def _generate_network_node_attributes_dict(self):
#         D = nx.Graph(np.absolute(self.prec))
#         
#         node_attributes_dict = {}
#         for i, node in enumerate(self.feature_names):
#             node_attributes_dict[node] = {}
#             node_attributes_dict[node]['weight'] = []
#             node_attributes_dict[node]['samp_indices'].append(i)
#         
#         nodedict = zip([self.feature_names[n] for n in D.nodes()],D.degree(D.nodes(), 'weight').values(), D.degree(D.nodes()).values())
#         nodelist = np.array(sorted(nodelist,key=lambda x: x[1], reverse = True))

    def network_node_attributes_file(self,folder = ''):
        D = self.D
        nodelist = zip([self.feature_names[n] for n in D.nodes()],D.degree(D.nodes(), 'weight').values(), D.degree(D.nodes()).values())
        nodelist = np.array(sorted(nodelist,key=lambda x: x[1], reverse = True))
        
        if not path.isdir(self.savedir + folder):
            mkdir(self.savedir + folder)
        nodelist_output = open(self.savedir+folder+'ntwrk_node_list.txt','w')
        nodelist_output.write('GLASSO NETWORK:\nFeature\tsum_of_abs_edge_weight\tdegree\n')
        for node in nodelist:
            nodelist_output.write('\t'.join(node))
            nodelist_output.write('\n')
    
    def network_node_attributes_plot(self, 
                                     colordict = None, 
                                     component = 1, 
                                     folder = '',
                                     tag = ''):
        
        D = nx.Graph(np.absolute(self.prec))
        nodelist = zip([self.feature_names[n] for n in D.nodes()],D.degree(D.nodes(), 'weight').values(), D.degree(D.nodes()).values())
        nodelist = np.array(sorted(nodelist,key=lambda x: x[1], reverse = True))
        
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111)
        for i in range(len(nodelist)):
            if float(nodelist[i][2] > 0):
                ax.plot(float(nodelist[i][1]),float(nodelist[i][2]), 's',color=colordict[nodelist[i][0]][1] if colordict else 'k')
                ax.annotate(s = '  '+nodelist[i][0], 
                                xy = (float(nodelist[i][1]),float(nodelist[i][2])), 
                                verticalalignment = 'left',
                                rotation = 45,
                                fontsize = '8')
        plt.title('features/loadings network: node attributes')
        ax.set_ylabel('degree')
        ax.set_xlabel('sum of abs edge weights')
        saveinfo = str(component)
        if tag is not '':
            saveinfo += '_'+tag
        if not path.isdir(self.savedir + folder):
            mkdir(self.savedir + folder)
        plt.savefig(self._avoid_overwrite(self.savedir+folder+'network_node_attrb_'+saveinfo+'.pdf'), bbox_inches="tight") #comp_'+str(component_num+1)
        plt.clf()
    
    def _set_edge_attributes(self, D):
        '''adds edge display attributes to D
        
        Parameters
        ----------
        D : networkx.classes.graph.Graph
            graph to be modified
    
        Notes
        ----------
        Adds edge display attributes to each edge in D:
            edge_style (dict)
                bold: all edges solid for equal visual representation
                subtle: positive correlations are solid, negative are dashed
            edge_color (dict)
                bold: positive correlations are colored GREEN ('g'), negative are RED ('r')
                subtle: all edges are colored GRAY ('k')
                
        The dictionary allows either display style (bold/subtle) to be selected when the 
        graph is actually drawn. Both styles carry the same information, but when PCA data
        is being overlayed on the network, it is easier to interpret if the edges are more
        subtly rendered.
        
        '''
        #make a new graph to see sign of edges
        Dprime = nx.Graph(self.prec)
        for a,b in Dprime.edges_iter(): 
            D[a][b]['edge_style'] = {'bold': 'solid', 'subtle': 'solid' if Dprime[a][b]['weight'] < 0.0 else 'dashed'}
            D[a][b]['edge_color'] = {'bold': 'g' if Dprime[a][b]['weight'] < 0.0 else 'r', 'subtle': 'k'}

    def _set_node_positions(self, D):
        '''assigns node positions for visualization
        
        Parameters
        ----------
        D : the network
                        
        Returns
        ----------
        pos : dict
            maps connected nodes (by number) to 2d positions in visualized graph
            keys : int
                number associated with each node in graph
            vals : numpy.ndarray [xpos, ypos]
                2d position between (0,0) and (1,1)
                
        Notes
        ----------
        unless specified by flag show_disconnected_nodes, unconnected nodes
            will not be displayed. Though they are not captured in the GLASSO network,
            it may be valuable to display the disconnected nodes as well 
            in order to visualize their PCA scores, and relative concentrations.
        
        connected nodes map to a graph entirely within the square from (0,0) to (1,1)
        unconnected nodes are roughly evenly spaced in a rounded-off square with side 
            length slightly longer than 1, so as to encircle the connected component
        '''
        num_iterations = 200
        spring_strength = 0.5
        layout = 'spring'
        
        #remove non-edge nodes so they don't clutter the layout
        nodes_with_edges = list(set(sum(D.edges(), ())))
        nodes_without_edges = [i for i in D.nodes() if i not in nodes_with_edges]

        D1 = D.copy()
        D1.remove_nodes_from(nodes_without_edges)
        
        if self.verbose: print '\t\tstarting layout...',layout,
        if layout == "circular":
            pos = nx.circular_layout(D1)
        else:
            pos = nx.fruchterman_reingold_layout(D1,
                                                 weight = 'weight',
                                                 iterations = num_iterations,
                                                 k = (spring_strength)*(1/sqrt(len(D.nodes()))))
        if self.verbose: print '...done' 
        
        D2 = D.copy()
        D2.remove_nodes_from(nodes_with_edges)
        pos2 = nx.circular_layout(D2)
        for k,v in pos2.items():
            pos2[k] = map(lambda w: w+0.5, map(lambda y: sin(pi*y)/sqrt(2), map(lambda x: x-0.5, v)))
        pos2.update(pos) #there is no overlap, so this just adds mappings in pos to pos2
        
        return pos2

    def _update_node_attributes(self, D,  node_size_scalar = 600):
        '''updates node size and shape attributes
        
        necessary to recalculate if there is a change in current_type_combo
        '''
        type_combo = self.current_type_combo
        change_in_conc_scaled = np.mean(self.data[self.types_dict[type_combo[0]]['samp_indices'],], axis=0) -\
              np.mean(self.data[self.types_dict[type_combo[1]]['samp_indices'],], axis=0)
        change_in_conc_raw = np.mean(self.rawdata[self.types_dict[type_combo[0]]['samp_indices'],], axis=0) -\
              np.mean(self.rawdata[self.types_dict[type_combo[1]]['samp_indices'],], axis=0)
        
        for n in D.nodes():
            D.node[n]['size']['delta_conc_scaled'] = node_size_scalar*change_in_conc_scaled[n]
            D.node[n]['size']['delta_conc_raw'] = node_size_scalar/5*change_in_conc_raw[n]
        
        for n in D.nodes():
            shape_scaled = 'v' if D.node[n]['size']['delta_conc_scaled'] > 0 else '^'
            if D.node[n]['size']['delta_conc_scaled'] == 0: 
                shape_scaled = 'o'
            
            shape_raw = 'v' if D.node[n]['size']['delta_conc_raw'] > 0 else '^'
            if D.node[n]['size']['delta_conc_raw'] == 0: 
                shape_raw = 'o'
                
            D.node[n]['shape']['arrow_raw'] = shape_raw
            D.node[n]['shape']['arrow_scaled'] = shape_scaled

    def _set_node_attributes(self, D, node_size_scalar = 600):
        '''TODO: document
        '''
        max_node_weight = max(D.degree(D.nodes(), 'weight').values())

        for n in D.nodes():
            D.node[n]['size'] = {}
            D.node[n]['size']['default'] = 150
            D.node[n]['size']['network_edge_weight'] = node_size_scalar*exp(D.degree(n, 'weight')/max_node_weight)

            D.node[n]['shape'] = {}
            D.node[n]['shape']['default'] = 'o'
            
            D.node[n]['color'] = {}
            D.node[n]['color']['default'] = (0.9, 0.9, 0.9, 0.25)
            
            feat = self.feature_names[n]
            D.node[n]['label'] = feat
        
        self._update_node_attributes(D, node_size_scalar = node_size_scalar)

    def _set_pca_node_colors(self, D, pca_data, component = 0):
        '''load colors from pca data into node attributes of graph for later access
        '''
        for n in D.nodes():
            feat = self.feature_names[n]
            D.node[n]['color']['pca_raw'] = pca_data[feat]['pca_raw_color'][component]
            D.node[n]['color']['pca_type'] = pca_data[feat]['pca_type_color'][component]
                
    def _draw_network(self, 
                      ax, 
                      nodes_to_show, 
                      draw_labels,
                      glasso_only,
                      node_color_selector = 'default',
                      node_size_selector = 'default',
                      label_vert_offset = 0.03):
        '''TODO: document
        ''' 
        D = self.D
        node_shape_selector = 'default'
        if node_size_selector == 'delta_conc_raw':
            node_shape_selector = 'arrow_raw'
        elif node_size_selector == 'delta_conc_scaled':
            node_shape_selector = 'arrow_scaled'

            
        edge_emphasis = 'subtle'
        max_edge_thickness = 8
        edge_alpha = 0.25
        if glasso_only:
            edge_emphasis = 'bold'
            max_edge_thickness += 2
            edge_alpha = 0.75
        
        edge_size_scalar = max_edge_thickness/max([D[x][y]['weight'] for x, y in D.edges()])
        for edge in D.edges():
            nx.draw_networkx_edges(D, self.pos,
                                   edgelist = [edge],
                                   ax=ax,
                                   edge_color = D.get_edge_data(*edge)['edge_color'][edge_emphasis],
                                   width = edge_size_scalar*D.get_edge_data(*edge)['weight'],
                                   alpha = edge_alpha,
                                   style = D.get_edge_data(*edge)['edge_style'][edge_emphasis])
        
        for nod in nodes_to_show:
            nx.draw_networkx_nodes(D, self.pos,
                                   nodelist = [nod],
                                   ax=ax,
                                   node_size = abs(D.node[nod]['size'][node_size_selector]),
                                   node_color = D.node[nod]['color'][node_color_selector],
                                   node_shape = D.node[nod]['shape'][node_shape_selector])
            if draw_labels:
                x,y = self.pos[nod]
                labelpos = (x,y+label_vert_offset)
                nx.draw_networkx_labels(D,
                                        {nod: labelpos},
                                        labels={nod: D.node[nod]['label']} ,
                                        ax=ax)
    
    def load_pca_data(self, 
                      pca_data,
                      pca_types_by_num, 
                      pca_component):
        '''TODO: document
        '''
        self.pca_data = pca_data
        self.current_type_combo = pca_types_by_num
        self.pca_component = pca_component
        
        D = self.D
        self._set_pca_node_colors(D, pca_data)
        self._update_node_attributes(D)

    def network_plot(self,
                     glasso_only            = False,
                     node_color_selector    = 'default',
                     node_size_selector     = 'default',
                     draw_labels            = True,
                     show_disconnected_nodes= True,
                     folder                 = '',
                     tag                    = ''):
        '''TODO: document, modularize further, add legend w colors/type for PCA ones
        '''
        D = self.D
        
        nodes_with_edges = list(set(sum(D.edges(), ())))
        nodes_to_show = D.nodes() if show_disconnected_nodes else nodes_with_edges

        if self.verbose:
            print '\tprecision matrix:', self.prec.shape 
            print '\tgraph:\t',len(nodes_with_edges),'connected nodes, ',self.prec.shape[0]-len(nodes_with_edges),'unconnected'
            print '\t\t',D.number_of_edges(), 'edges'
                    
        fig = plt.figure(figsize=(15,15))
        plt.axis("off")
        ax = fig.add_subplot(111)
        ax.set_xlim([-0.05,1.05])
        ax.set_ylim([-0.05,1.05])
        if show_disconnected_nodes:
            ax.set_xlim([-0.25,1.25])
            ax.set_ylim([-0.25,1.25])
        
        self._draw_network(ax, 
                           nodes_to_show, 
                           draw_labels, 
                           glasso_only = glasso_only, 
                           node_color_selector = node_color_selector, 
                           node_size_selector = node_size_selector)

        saveinfo = 'nwk'
        if node_size_selector is not 'default':
            saveinfo +='_'+node_size_selector
        if node_color_selector is not 'default':
            saveinfo +='_'+node_color_selector
        if show_disconnected_nodes:
            saveinfo +='_all'
        if not glasso_only:
            saveinfo +='_comp'+str(self.component+1)+'_types'+reduce(lambda w,y: w+y, map(lambda x: str(x), self.current_type_combo))
        
        if tag is not '':
            saveinfo += '_'+tag

        plt.title('glasso network graph: GL alpha '+str(self.alpha)+'\n'+saveinfo, fontsize = 20)
        if not path.isdir(self.savedir + folder):
            mkdir(self.savedir + folder)
        plt.savefig(self._avoid_overwrite(self.savedir+folder+saveinfo+'.pdf'), bbox_inches="tight")
        plt.clf()

    def _avoid_overwrite(self, save_location):
        inc = 1
        while path.isfile(save_location):
            save_location = save_location.rstrip('.pdf')
            if inc != 1:
                save_location = save_location.rsplit('_',1)[0]
            save_location = save_location + '_'+str(inc)+'.pdf'
            inc += 1
        return save_location
    
    def _initialize_network(self):
        if self.verbose:
            print "\tinitializing network ..."
        D = nx.Graph(np.absolute(self.prec))
        
        if self.pos == None:
            self.pos = self._set_node_positions(D)
            self._save_network_graph(D)
            
        self._set_edge_attributes(D)
        self._set_node_attributes(D)
        
        if self.verbose:
            print "\t... done"
        return D

    def _save_network_graph(self):
        nx.write_gml(self.D,self._avoid_overwrite(self.savedir+"network.gml"))
        return

    def _run_glasso_cv(self):
        if self.verbose:
            print "\tfinding GL alpha with cross-validation..."
            print "\t\tworking.... (this can take a while)"
        glcv = GraphLassoCV()
        glcv.fit(self.data)
        if self.verbose: 
            print "\t\t...done"
        return glcv.alpha_

    def __init__(self, 
                 savedir, 
                 data, 
                 rawdata, 
                 types, 
                 feature_names, 
                 alpha = 0, 
                 pos = None,
                 verbose = True):
        '''
        Constructor
        '''

            
        self.savedir = savedir
        self.data = data
        self.rawdata = rawdata
        self.types_dict = types
        self.feature_names = feature_names
        self.verbose = verbose
        self.pos = pos

        self.current_type_combo = [0,1]    
        self.component = 0
        
        if self.verbose: 
            print "GLASSO:"
        
        if alpha <= 0:
            alpha = self._run_glasso_cv()


        gl = GraphLasso(alpha = alpha)
        gl.fit(data)

        glcov = gl.covariance_
        glprec = gl.precision_
        #remove diagonal values
        glprec = glprec - np.diag(np.diag(glprec))
        
        self.prec = glprec
        self.cov = glcov
        self.alpha = alpha
        
        if verbose:
            print '\tglasso alpha:', alpha
            print "\tprec:",glprec.shape, ", nonzero vals:",self._count_nz_offd(glprec)
            print "\tcov:",glcov.shape, ", nonzero vals:",self._count_nz_offd(glcov)

        self.D = self._initialize_network()
