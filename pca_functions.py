'''
Created on Jun 11, 2016

@author: Jared
'''
from sklearn.decomposition import PCA as skl_pca
from sklearn.decomposition import SparsePCA as skl_sparse_pca
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
from os import path, mkdir

# import mpl_toolkits.mplot3d


class pca_bundle(object):
    '''bundle PCA functions for readability
    '''
    def ordered_set(self, seq):
        '''collapse list into set, preserving order
        from: http://www.peterbe.com/plog/uniqifiers-benchmark and comments there
        '''
        seen = {}
        result = []
        for item in seq:
            if not seen.__contains__(item):
                seen[item] = 1
                result.append(item)
        return result
    
    def type_to_color(self, type_label):
        return [x['color'] for x in self.types if x['name'] == type_label][0]
    
#     def init_type_to_color(self, sample_types_set, COLORMAP = 'spectral'):
#         #COLORMAP options: spectral, viridis, plamsa, magma
#     ##strip sample names after first '_' to isolate group labels 
#         sample_labels_list = map(lambda x: x.split("_")[0], self.sample_names)
#         sample_labels_set = self.ordered_set(sample_labels_list)        
#         print "\t", len(sample_labels_set),"labels found in samples:",sample_labels_set
#         sample_labels_set = sample_types_set
#         ##GENERATE LIST OF COLORS to be assigned to different types: range (0,1)
#         label_colors_list = [(float(i) +.5) / len(sample_labels_set) for i in range(len(sample_labels_set))]
# #         print "\tcolorrange:",label_colors_list
#         
#         #MAP TO COLORMAP (RGB space)
#         #'spectral' chosen because different points correspond to different categories
#         cm = plt.get_cmap(COLORMAP)
#         label_colors_list = [cm(x) for x in label_colors_list]
# #         print "\tlabel_colors_list mapped:",label_colors_list
#          
#         return dict(zip(sample_labels_set,label_colors_list))
    
    def comps_explained_var_plot(self, folder='', tag=''):
        types = 'all' if self.run_types_by_num == range(len(self.types)) else str(self.run_types_by_num)
        
        fig = plt.figure(figsize=(4, 4))
        ax = fig.add_subplot(111) 
        ax.plot(range(1, len(self.comp_var) + 1), self.comp_var)
        ax.set_xlabel('component #')
        ax.set_title('PCA Explained Var: types '+types)
        ax.set_ylabel('explained variance')
        ax.set_xbound((0.5,len(self.comp_var) + 0.5))
        ax.set_ybound((0,max(self.comp_var) + 0.1))
        ax.set_xticks(range(1,len(self.comp_var)+1))
        ax.axhline(y=self.comp_var_threshold, xmin = 0, xmax = len(self.comp_var) + 0.5, linestyle='dotted')        
        if not path.isdir(self.savedir + folder):
            mkdir(self.savedir + folder)
        plt.savefig(self._avoid_overwrite(self.savedir + folder + 'pca_explained_var' + tag + '.pdf'), bbox_inches="tight")
        plt.clf()
        if self.verbose:
            print 'PCA Explained Variance Graph\n'

    def comps_plot_1d(self,
                      folder='',
                      show_legend = True,
                      tag=''):
        '''plot pca samples by component score (in 1 dimension)
        
            INPUTS:
                sample_scores: 2d matrix of samples sample_scores (num samples by num components)
                samples: list of sample names
                components: list of components to be plotted (can be list of 1 element)
            
        '''
        scores = self.sample_scores
        comps = self.sample_scores.shape[1]
        
        SAMPLE_ALPHA = 1.0
        print 'PCA Component Scores - 1d: \n\tcomponents', comps
        
        fig = plt.figure(figsize=(10, 1.5 + comps))
        ax = fig.add_subplot(111)  # , projection = '3d')
        seen = set()

        for compnum in range(comps):  # for each component
            for n in range(scores.shape[0]):  # iterate through sample score positions
                type_label = self.sample_name_to_type(self.run_sample_names[n])
                colorRBGA = self.type_to_color(type_label)
                

        #         print type_label, " : ", colorRBGA
    #                ax.plot(sample_scores[0:5,n],np.ones(5)*n,'d', markersize=30, color='blue',\
    #                      alpha=0.5, type_label='hESC')
                ax.plot([scores[n, compnum]], [compnum], 'd', markersize=15, color=colorRBGA, alpha=SAMPLE_ALPHA, label=type_label if type_label not in seen else "")
                seen.add(type_label)
    #     ax.set_xlabel('comp: ' + str(comps[0]+1))
    #     ax.set_ylabel('comp: '+ str(comps[1]+1))
    #     ax.set_zlabel('comp: '+ str(comps[2]+1))
    #     plt.subplots_adjust(left=.05, right=0.95, top=.95, bottom=0.05)
    #     plt.legend(loc = 'center left',bbox_to_anchor = (.85,.4,1,1))
         
    #     ##PLOT Y AXIS
        yb = ax.get_ybound()
        xbtop, xbbot = ax.get_xbound()
        xbm = max(abs(xbtop), abs(xbbot))
        ax.set_xlim(-7, 7)
        yAxisLine = ((0, 0), (min(yb[0], 0) - 1, max(yb[1], 0) + 1))
        ax.plot(yAxisLine[0], yAxisLine[1], 'k', alpha=0.3)
        
        plt.title("Score Distributions for First " + str(comps) + " Components")
        ax.set_xlabel('Score Distribution')
        ax.set_ylabel('Component #')
        ax.set_ybound(lower=-1, upper=comps)
        
        ax.set_yticks(range(comps))
        ylabels = [x + 1 for x in range(comps)]
        ax.set_yticklabels(ylabels)
    #     plt.tight_layout()
        if show_legend:
            plt.legend(loc='center left', bbox_to_anchor=(.85, .4, 1, 1))
    #     if tag is not '' : tag = tag + " " # so no awkward space at start of file
        if not path.isdir(self.savedir + folder):
            mkdir(self.savedir + folder)
        plt.savefig(self._avoid_overwrite(self.savedir + folder + 'comps_plot_1d' + tag + '.pdf'), bbox_inches="tight")
        plt.clf()
    
    def pca_comps_plot_2d(self, samples, components='12', tag=None):
        scores = self.sample_scores
        
        comps = map(lambda x: int(x) - 1, components)  # MINUS 1 because indexing starts at 1
        print "\nscores plot 2d: \n\tanalyzing components", comps
        figname = '2d sample_scores plot for comps:' + components
        
        
#         sample_types = [self.sample_name_to_type(x) for x in samples]    
#         sample_set = self.ordered_set(sample_types)
#         print "\t", len(sample_set),"labels found in samples:",sample_set
#         
#         #generate list of colors to be assigned to different types: range (0,1)
#         color_range = [(float(i) / len(sample_set)) +.5 / len(sample_set) for i in range(len(sample_set))]
#     #     print "\ncolorrange:",color_range
#         #map to colormap, 'spectral' chosen because different points correspond to different categories
#         cm = plt.get_cmap('spectral')
#         color_range = [cm(x) for x in color_range]
#     #     print "\ncolorrange mapped:",color_range
#         
#         color_lookup = dict(zip(sample_set,color_range))
#     #     print '\ncolor_lookup:'
#     #     for i in color_lookup.items():
#     #         print i
#         
#         type_loc = {pos:val for pos,val in enumerate(sample_types)}
    #     print "\ntype loc:"
    #     for i in type_loc.iteritems():
    #         print i
    #     print "\n\n"
        
        fig = plt.figure(figsize=(10, 10))
        cm = plt.get_cmap('spectral')
        ax = fig.add_subplot(111)
        seen = set()
    
        for n in range(scores.shape[0]):
            type_label = self.sample_name_to_type(self.run_sample_names[n])
            colorRBGA = self.type_to_color(type_label)
    #         print type_label, " : ", colorRBGA
            ax.plot([scores[n, comps[0]]], [scores[n, comps[1]]], 'o', markersize=8, \
                color=colorRBGA, alpha=0.5, label=type_label if type_label not in seen else "")
            seen.add(type_label)
        ax.set_xlabel('comp: ' + str(comps[0] + 1))
        ax.set_ylabel('comp: ' + str(comps[1] + 1))
        plt.title(figname)
        plt.subplots_adjust(left=.05, right=0.95, top=.95, bottom=0.05)
        plt.legend(loc='best')  # ,bbox_to_anchor = (.85,.4,1,1))
        
        # #PLOT AXES
        xb = ax.get_xbound()
        yb = ax.get_ybound()
        xAxisLine = ((min(xb[0], 0), max(xb[1], 0)), (0, 0), (0, 0))
        yAxisLine = ((0, 0), (min(yb[0], 0), max(yb[1], 0)), (0, 0))
        ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'k', alpha=0.3)
        ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'k', alpha=0.3)
        
        plt.savefig(self._avoid_overwrite(self.savedir + tag + 'pca_comps_plot_2d_' + components + '.pdf'), bbox_inches="tight")
        plt.clf()

    def comps_plot_3d(self,
                      components='123',
                      projections=True,
                      show_legend=True,
                      folder='',
                      tag=None):
        from mpl_toolkits.mplot3d import Axes3D as ax3d
#         samples = self.sample_names,
        scores = self.sample_scores
        
        comps = [int(components[0]) - 1, int(components[1]) - 1, int(components[2]) - 1]  # MINUS 1 because indexing starts at 1
        print "\nscores plot 3d: \n\tanalyzing components", comps
        figname = '3d sample_scores plot for comps:' + components
    #     fig = plt.figure(figsize=(10,10))
    #     print "sample_scores:\n",sample_scores
    #     print '\nsamples',samples
        
#         sample_types = [self.sample_name_to_type(x) for x in samples]    
#         sample_set = self.ordered_set(sample_types)
#         print "\t", len(sample_set),"labels found in samples:",sample_set
#         
#         #generate list of colors to be assigned to different types: range (0,1)
#         color_range = [(float(i) / len(sample_set)) +.5 / len(sample_set) for i in range(len(sample_set))]
#     #     print "\ncolorrange:",color_range
#         #map to colormap, 'spectral' chosen because different points correspond to different categories
#         cm = plt.get_cmap('spectral')
#         #TODO: choose best colormap
#         color_range = [cm(x) for x in color_range]
#     #     print "\ncolorrange mapped:",color_range
#         
#         color_lookup = dict(zip(sample_set,color_range))
#     #     print '\ncolor_lookup:'
#     #     for i in color_lookup.items():
#     #         print i
#         
#         type_loc = {pos:val for pos,val in enumerate(sample_types)}
    #     print "\ntype loc:"
    #     for i in type_loc.iteritems():
    #         print i
    #     print "\n\n"
        
        fig = plt.figure(figsize=(10, 10))
        cm = plt.get_cmap('spectral')
        ax = fig.add_subplot(111, projection='3d')
        seen = set()
    
        for n in range(scores.shape[0]):
            type_label = self.sample_name_to_type(self.run_sample_names[n])
            colorRBGA = self.type_to_color(type_label)
    #         print type_label, " : ", colorRBGA
            ax.plot([scores[n, comps[0]]], [scores[n, comps[1]]], [scores[n, comps[2]]], 'o', markersize=8, \
                color=colorRBGA, alpha=1, label=type_label if type_label not in seen else "", zorder=2)
            seen.add(type_label)
            
        xb = ax.get_xbound()
        yb = ax.get_ybound()
        zb = ax.get_zbound()
        if projections == True:
            for n in range(scores.shape[0]):
                type_label = self.sample_name_to_type(self.run_sample_names[n])
                colorRBGA = self.type_to_color(type_label)
                ax.plot([scores[n, comps[0]]], [scores[n, comps[2]]], '<', color=colorRBGA, zdir='y', \
                        zs=yb[1], alpha=.3, markersize=4, zorder=1, mec=colorRBGA)
                ax.plot([scores[n, comps[1]]], [scores[n, comps[2]]], '>', color=colorRBGA, zdir='x', \
                        zs=xb[0], alpha=.3, markersize=4, zorder=1, mec=colorRBGA)
                ax.plot([scores[n, comps[0]]], [scores[n, comps[1]]], '^', color=colorRBGA, zdir='z', \
                        zs=zb[0], alpha=.3, markersize=4, zorder=1, mec=colorRBGA)
        ax.set_xlabel('comp: ' + str(comps[0] + 1))
        ax.set_ylabel('comp: ' + str(comps[1] + 1))
        ax.set_zlabel('comp: ' + str(comps[2] + 1))
    
    #     plt.subplots_adjust(left=.05, right=0.95, top=.95, bottom=0.05)
        if show_legend:
            plt.legend(loc='center left', bbox_to_anchor=(.85, .4, 1, 1))
        
        # #PLOT AXES
    
        xAxisLine = ((min(xb[0], 0), max(xb[1], 0)), (0, 0), (0, 0))
        yAxisLine = ((0, 0), (min(yb[0], 0), max(yb[1], 0)), (0, 0))
        zAxisLine = ((0, 0), (0, 0), (min(zb[0], 0), max(zb[1], 0)))
        ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'k', alpha=0.3, zorder=0)
        ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'k', alpha=0.3, zorder=0)
        ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'k', alpha=0.3, zorder=0)
        ax.set_xbound(lower=xb[0], upper=xb[1])
        ax.set_ybound(lower=yb[0], upper=yb[1])
        ax.set_zbound(lower=zb[0], upper=zb[1])
        ax.set_autoscale_on(False)
        if show_legend:
            plt.title(figname)
        if not path.isdir(self.savedir + folder):
            mkdir(self.savedir + folder)
        plt.savefig(self._avoid_overwrite(self.savedir + folder + 'pca_comps_plot_3d_' + components + tag + '.pdf'), bbox_inches="tight")
    #     plt.show()
        plt.clf()
    
    def loadings_plot_1d(self, component_num=0, cutoff_percentile=66, tag=None, no_figure=False):
        features = self.feature_names
        feature_scores = self.feature_scores
        label_nodes = False
        
        print "loadings_stem_plot: \n\tanalyzing component", component_num
#         print "feature_scores", feature_scores
    #     print "feat:",features
    #     OLD LOADVALS: loadvals = [x for x in feature_scores[:,component_num] if abs(x) > cutoff] #ONLY TAKES > CUTOFF VALUES
        loadvals = [x for x in feature_scores[:, component_num]]  # select component column only
        abs_loadvals = [abs(x) for x in loadvals]
        cm = plt.get_cmap('RdBu')
        pca_color_vals = [cm(x + 0.5)[:-1] for x in loadvals]  # '+0.5' centers the colormap around white
        type_colors = [self.type_to_color(x) for x in self.current_run_types]
        print 'type to color assignment: ', self.current_run_types, [self.type_to_color(x) for x in self.current_run_types]
        cm = self._diverge_map(high=type_colors[1][:-1], low=type_colors[0][:-1])

        type_color_vals = [cm(x + 0.5)[:-1] for x in loadvals]  # '+0.5' centers the colormap around white
        print 'loadvals', len(loadvals), loadvals
#         print 'features',len(features),features
#         print 'colors',len(pca_color_vals),pca_color_vals
        feature_data = sorted(zip(features, loadvals, pca_color_vals, type_color_vals), key=lambda x: x[1])

        feature_data_dict = {}
        for i in range(len(feature_data)):
            feature_data_dict[feature_data[i][0]] = feature_data[i][1:]

        if no_figure:
            return feature_data_dict, feature_data, cm
        # #generate plot
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111)
        cutoff = np.percentile([abs(x) for x in loadvals], cutoff_percentile)
    
        print '\tcutoff percentile:', cutoff_percentile
        print '\tabsolute cutoff:', cutoff
    
        for i in range(len(feature_data)):
            if abs(feature_data[i][1]) > cutoff:  # #if this point beats cutoff
                ax.plot(feature_data[i][1], i, 's', color=feature_data[i][3])
    #             if feature_data[i][1] > 0: 
    #                 h_align = 'left' #vertical
    #             else:
    #                 h_align = 'right'
    #             ax.text(feature_data[i][1],i, ' '+feature_data[i][0], rotation = 0, fontsize = '8', rotation_mode = 'anchor')
#                 if abs(feature_data[i][1]) > 0.15:
                if label_nodes:
                    ax.annotate(s='  ' + feature_data[i][0] + '  ',
                                xy=(feature_data[i][1], i),
                                horizontalalignment='left' if feature_data[i][1] > 0 else 'right',
                                verticalalignment='center',
                                fontsize='10')  # if feature_data[i][0] not in ['Acetate', '2-Phosphoglycerate','4-Aminobutyrate','Ethanolamine'] else '6')
            else:
                ax.plot(feature_data[i][1], i, 'o', color=feature_data[i][2])
            mn = min(0, feature_data[i][1])
            mx = max(0, feature_data[i][1])
            ax.hlines(i, mn, mx, linestyles='dotted')
        
        ax.vlines(0, -1, len(features), linestyles='solid')
        ax.set_yticks(range(len(features)))
        ax.set_yticklabels([x[0] for x in feature_data], fontsize='6')
        
        cuttitle = ''
        if cutoff != 0:
            ax.vlines(cutoff, -1, len(features), linestyles='--')
            ax.vlines(-cutoff, -1, len(features), linestyles='--')
            cuttitle = " cut " + str(cutoff_percentile) + '%'
    
    
        # TODO: add colorbar to stemplot
    
        plt.title('PCA Component ' + str(component_num + 1) + " Loadings: " + cuttitle)
        ax.set_ylabel(str(len(features)) + ' nonzero features')
        ax.set_ybound(lower=-1, upper=len(features))
        ax.set_xbound(lower=-max(abs_loadvals) - .1, upper=max(abs_loadvals) + .1)
        ax.set_xlabel('Component ' + str(component_num + 1) + " Scores")
        #     plt.subplots_adjust(wspace = 0.7, hspace = 0.7) 
        #     plt.subplots_adjust(left=.1, right=0.9, top=.9, bottom=0.2)
    #     plt.tight_layout()
    #     cut = ''
    #     if cutoff is not 0: 
    #     cut = ' cut_'+str(cutoff)
        plt.savefig(self._avoid_overwrite(self.savedir + tag + 'pca_loadings_1d_comp' + str(component_num) + '.pdf'), bbox_inches="tight")  # comp_'+str(component_num+1)
        plt.clf()
        #=======multiplot experimentation===========================================================
    #     print type(fig)
    #     print type(ax)
    #     fig1 = plt.figure(figsize=(10,10))
    # #     ax1 = fig1.add_subplot(figure = fig)
    # #     ax2 = fig1.add_subplot(222)
    # #     ax3 = fig1.add_subplot(212)
    #     plt.savefig(savedir+'ZZZpca_loadings_1d_comp'+str(component_num)+'.pdf', bbox_inches="tight") #comp_'+str(component_num+1)
    
        #===========================================================================
    #=====horizontal version========================================================
    #     for i in range(len(feature_data)):
    #         if abs(feature_data[i][1]) > cutoff: ##if this point beats cutoff
    #             ax.plot(i,feature_data[i][1],'s',color=feature_data[i][2])
    #             if feature_data[i][1] > 0: 
    #                 rot = 90 #vertical
    #             else:
    #                 rot = 300
    #             ax.text(i,feature_data[i][1], ' '+feature_data[i][0], rotation = rot, fontsize = '8', rotation_mode = 'anchor')
    #         else:
    #             ax.plot(i,feature_data[i][1],'o',color=feature_data[i][2])
    #         mn = min(0,feature_data[i][1])
    #         mx = max(0,feature_data[i][1])
    #         ax.vlines(i,mn,mx,linestyles = 'dotted')
    #     
    #     ax.hlines(0,-1,len(features),linestyles = 'solid')
    #     ax.set_xticks(range(len(features)))
    #     ax.set_xticklabels([x[0] for x in feature_data], rotation=270, fontsize = '6')
    #     
    #     cuttitle = ''
    #     if cutoff != 0:
    #         ax.hlines(cutoff,-1,len(features),linestyles = '--')
    #         ax.hlines(-cutoff,-1,len(features),linestyles = '--')
    #         cuttitle = " cutoff: "+str(cutoff)
    # 
    # 
    # 
    #     plt.title('PCA Component '+str(component_num+1)+" "+tag + cuttitle)
    #     ax.set_xlabel(str(len(features))+' nonzero features')
    #     ax.set_xbound(lower = -1, upper = len(features))
    #     ax.set_ybound(lower = -max(abs_loadvals)-.1, upper = max(abs_loadvals)+.1)
    #     ax.set_ylabel('Component '+str(component_num+1)+" Scores")
    #     #     plt.subplots_adjust(wspace = 0.7, hspace = 0.7) 
    #     #     plt.subplots_adjust(left=.1, right=0.9, top=.9, bottom=0.2)
    # #     plt.tight_layout()
    # #     cut = ''
    # #     if cutoff is not 0: 
    #     cut = ' cut_'+str(cutoff)
    #     if tag is not '' : tag = tag + " " # so no awkward space at start of file
    #     plt.savefig(savedir+'loadings_plot_1d.pdf', bbox_inches="tight") #comp_'+str(component_num+1)
    #     plt.clf()
    #===============================================================================

        return feature_data_dict, feature_data, cm
    
    def _boxplot_colors(self, bp):
        '''manages custom coloring of boxplots
        '''
        # set all defaults to black
        for k in bp.keys():
            for x in bp[k]:
                x.set(color='k')
                
        # set boxes to type color, with light alpha
        for i, box in enumerate(bp['boxes']):
            col = tuple(list(self.types[i]['color'][0:3]) + [0.25])
            box.set(color=col)
            
        # set medians to type color
        for i, med in enumerate(bp['medians']):
            med.set(color=self.types[i]['color'], linewidth=2)

    def single_boxplot(self,
                       feature_name,
                       boxplot_folder='boxplots/',
                       show_raw=False):
        '''creates a boxplot of given feature, and saves it in boxplot_folder
        '''
        if not path.isdir(self.savedir + boxplot_folder):
            mkdir(self.savedir + boxplot_folder)
        
        fig = plt.figure(figsize=(4, 4) if show_raw else (4, 2))
        ax = fig.add_subplot(212 if show_raw else 111)
        
        boxplot_data = []
        boxplot_rawdata = []
        feat_index = self.feature_names.index(feature_name)
        for t in self.types:
            boxplot_data.append(list(self.data[t['samp_indices'], [feat_index]]))
            boxplot_rawdata.append(list(self.rawdata[t['samp_indices'], [feat_index]]))

        ax.set_title('normalized', loc='right', fontsize=11)
        bp = ax.boxplot(boxplot_data, patch_artist=True)
        self._boxplot_colors(bp) 
        
        ax.set_xticklabels([x['name'] for x in self.types])
        ax.get_xaxis().set_ticks_position('none')  # remove y ticks
        ax.get_yaxis().set_ticks_position('none')  # remove x ticks
        ax.set_ybound(lower=self.min_data * 1.1, upper=self.max_data * 1.1)
        
        if show_raw:
            ax1 = fig.add_subplot(211)
            ax1.set_xticklabels(['' for x in self.types])
            ax1.set_title('raw', loc='right', fontsize=11)
            ax1.get_xaxis().set_ticks_position('none')  # remove y ticks
            ax1.get_yaxis().set_ticks_position('none')  # remove x ticks
            bp1 = ax1.boxplot(boxplot_rawdata, patch_artist=True)
            self._boxplot_colors(bp1) 

#             ax1.set_ybound(lower = self.min_rawdata*1.1, upper = self.max_rawdata*1.1)
        fig.suptitle(feature_name, y=1 if show_raw else 1.05)
        plt.savefig(self._avoid_overwrite(self.savedir + boxplot_folder + feature_name + '.pdf'), bbox_inches="tight")
        plt.close(fig)
        
    def gen_boxplots(self, show_raw=False):
        if self.verbose: print 'generating boxplots...(can take some time)',
        for feat in self.feature_names:
            self.single_boxplot(feat, show_raw=show_raw)
        if self.verbose: print '...done'

    def boxplot(self, component=0, number_boxplots=3, tag=None):
        SAMPLE_ALPHA = 1.0
        print 'BOXPLOT:'
        print '\tcomponent:', component
        np_feat_names = np.array([self.feature_names]).T
        comp_score = np.array([self.feature_scores[:, component]]).T
#         print comp_score.shape
#         print np_feat_names.shape
        feature_names_comp = np.concatenate((np_feat_names, comp_score), axis=1)
        print 'feature names:'
        print self.feature_names[0:3]
        print self.feature_scores[0:3]
        print comp_score[0:3]
        a = feature_names_comp.tolist()
        print 'feature names comp'
        for i in a[0:3]:
            print i
        a = sorted(a, key=lambda x: abs(float(x[1])), reverse=True)
        print 'feature names comp sorted'
        for i in a[0:3]:
            print i
        feature_names_to_plot = [x[0] for x in a[0:number_boxplots]]
        print feature_names_to_plot
        
        fig, axarr = plt.subplots(number_boxplots, 2, sharex=True, sharey=False)
        fig.suptitle('boxplot comp' + str(component) + ' significant features')

        for j, feat in enumerate(feature_names_to_plot):
            boxplot_data = []
            feat_index = self.feature_names.index(feat)
            type_data_dict = {}
            orderseen = []
            for i, sample in enumerate(self.run_sample_names):
                sample_type = self.sample_name_to_type(sample)
                if sample_type in type_data_dict.keys():
                    type_data_dict[sample_type] += [self.data[i][feat_index]]
                else:
                    type_data_dict[sample_type] = [self.data[i][feat_index]]
                    orderseen += [sample_type]
#             print feat,
#             print feat_index
#             print type_data_dict.items()
#             print orderseen
#             print '++++\n'
            for typesamp in orderseen:
                boxplot_data += [type_data_dict[typesamp]]
            ax = axarr[j][0]
            ax.boxplot(boxplot_data)
            ax.set_xticklabels(orderseen)
            ax.get_xaxis().tick_bottom()  # remove top ticks
            ax.get_yaxis().tick_left()  # remove right ticks
            ax.set_title(feat)
            axarr[j][1].set_axis_off()

        # ADD COMPONENT SAMPLE PCA SCORES, COLORBAR
        #=======================================================================
        ax_samp = fig.add_subplot(1, 2, 2)

        scores = self.sample_scores
        seen = set()
#         print "SAMPLE SCORES:",scores.shape
        for n in range(scores.shape[0]):  # iterate through sample score positions
            type_label = self.sample_name_to_type(self.run_sample_names[n])
            colorRBGA = self.type_to_color(type_label)
    #         print type_label, " : ", colorRBGA
#                ax.plot(sample_scores[0:5,n],np.ones(5)*n,'d', markersize=30, color='blue',\
#                      alpha=0.5, type_label='hESC')
#             print [0,scores[n,component]]
            ax_samp.plot(0, scores[n, component], 'D', markersize=15, color=colorRBGA, alpha=SAMPLE_ALPHA, label=type_label if type_label not in seen else "")
            seen.add(type_label)  #=======================================================================
        ax_samp.set_xlabel('Sample Scores\tFeature Scores')
        ax_samp.get_xaxis().set_ticks([])
        
        feat_scores = self.feature_scores
#         print "FEAT SCORES:",feat_scores.shape
        
        loadvals = [x for x in feat_scores[:, component]]
#         abs_loadvals = [abs(x) for x in loadvals]
        cm = plt.get_cmap('RdBu')
        color_range = [cm(x + .5)[:-1] for x in loadvals]
    #     print 'loadvals',len(loadvals),loadvals
    #     print 'features',len(features),features
    #     print 'colors',len(color_range),color_range
        feature_data = sorted(zip(self.feature_names, loadvals, color_range), key=lambda x: x[1])
        #=======================================================================
#         for i in range(len(feature_data)):
#             if abs(feature_data[i][1]) > cutoff: ##if this point beats cutoff
#                 ax.plot(feature_data[i][1],i,'s',color=feature_data[i][2])
#     #             if feature_data[i][1] > 0: 
#     #                 h_align = 'left' #vertical
#     #             else:
#     #                 h_align = 'right'
#     #             ax.text(feature_data[i][1],i, ' '+feature_data[i][0], rotation = 0, fontsize = '8', rotation_mode = 'anchor')
#                 ax.annotate(s = '  '+feature_data[i][0]+'  ', 
#                             xy = (feature_data[i][1],i), 
#                             horizontalalignment = 'left' if feature_data[i][1] > 0 else 'right', 
#                             verticalalignment = 'center',
#                             fontsize = '8')
#             else:
#                 ax.plot(feature_data[i][1],i,'o',color=feature_data[i][2])
#             mn = min(0,feature_data[i][1])
#             mx = max(0,feature_data[i][1])
#             ax.hlines(i,mn,mx,linestyles = 'dotted')
#         
#         ax.vlines(0,-1,len(features),linestyles = 'solid')
#         ax.set_yticks(range(len(features)))
#         ax.set_yticklabels([x[0] for x in feature_data], fontsize = '6')
        #=======================================================================
        SCALAR = 5
        for n in range(len(feature_data)):  # iterate through feature score positions
    #         print type_label, " : ", colorRBGA
#                ax.plot(sample_scores[0:5,n],np.ones(5)*n,'d', markersize=30, color='blue',\
#                      alpha=0.5, type_label='hESC')
#             print [0,scores[n,component]]
            this_feature = feature_data[n][0]
#             if this_feature in feature_names_to_plot:
#                 print this_feature,'DFADFDS', feature_data[n][1]
            ax_samp.plot(1, feature_data[n][1] * SCALAR, \
                         's' if this_feature in feature_names_to_plot else '.', \
                         color=feature_data[n][2], \
                         label=this_feature if this_feature in feature_names_to_plot else None)
            if this_feature in feature_names_to_plot:
                ax_samp.annotate(s='  ' + feature_data[n][0] + '  ',
                            xy=(1, feature_data[n][1] * SCALAR),
#                             horizontalalignment = 'left' if feature_data[i][1] > 0 else 'right', 
#                             verticalalignment = 'center',
                            fontsize='8')
#             ax_feat.plot([0,scores[n,component]],'D',markersize=15, color=colorRBGA, alpha=0.5, label = type_label if type_label not in seen else "")
#             seen.add(type_label)
        ax_samp.set_xbound((-0.5, 1.5))
        plt.savefig(self._avoid_overwrite(self.savedir + str(tag) + 'boxplot_comp' + str(component) + '.pdf'), bbox_inches="tight")
        plt.legend(loc='best', bbox_to_anchor=(.85, .4, 1, 1))

        plt.clf()
#             print self.data[:,self.feature_names.index(feat)]
#         print self.data.shape        

    def boxplot_old(self, component=0, number_boxplots=3, tag=None):
        SAMPLE_ALPHA = 1.0
        print 'BOXPLOT:'
        print '\tcomponent:', component
        np_feat_names = np.array([self.feature_names]).T
        comp_score = np.array([self.feature_scores[:, component]]).T
#         print comp_score.shape
#         print np_feat_names.shape
        feature_names_comp = np.concatenate((np_feat_names, comp_score), axis=1)
#         print self.feature_names[0:3]
#         print self.feature_scores[0:3]
#         print comp_score[0:3]
        a = feature_names_comp.tolist()
#         for i in a:
#             print i
        a = sorted(a, key=lambda x: abs(float(x[1])), reverse=True)
#         for i in a:
#             print i
        feature_names_to_plot = [x[0] for x in a[0:number_boxplots]]
#         print feature_names_to_plot
        
        fig, axarr = plt.subplots(number_boxplots, 2, sharex=True, sharey=False)
        fig.suptitle('boxplot comp' + str(component) + ' significant features')

        for j, feat in enumerate(feature_names_to_plot):
            boxplot_data = []
            feat_index = self.feature_names.index(feat)
            type_data_dict = {}
            orderseen = []
            for i, sample in enumerate(self.run_sample_names):
                sample_type = self.sample_name_to_type(sample)
                if sample_type in type_data_dict.keys():
                    type_data_dict[sample_type] += [self.data[i][feat_index]]
                else:
                    type_data_dict[sample_type] = [self.data[i][feat_index]]
                    orderseen += [sample_type]
#             print feat,
#             print feat_index
#             print type_data_dict.items()
#             print orderseen
#             print '++++\n'
            for typesamp in orderseen:
                boxplot_data += [type_data_dict[typesamp]]
            ax = axarr[j][0]
            ax.boxplot(boxplot_data)
            ax.set_xticklabels(orderseen)
            ax.get_xaxis().tick_bottom()  # remove top ticks
            ax.get_yaxis().tick_left()  # remove right ticks
            ax.set_title(feat)
            axarr[j][1].set_axis_off()

        # ADD COMPONENT SAMPLE PCA SCORES, COLORBAR
        #=======================================================================
        ax_samp = fig.add_subplot(1, 2, 2)

        scores = self.sample_scores
        seen = set()
#         print "SAMPLE SCORES:",scores.shape
        for n in range(scores.shape[0]):  # iterate through sample score positions
            type_label = self.sample_name_to_type(self.run_sample_names[n])
            colorRBGA = self.type_to_color(type_label)
    #         print type_label, " : ", colorRBGA
#                ax.plot(sample_scores[0:5,n],np.ones(5)*n,'d', markersize=30, color='blue',\
#                      alpha=0.5, type_label='hESC')
#             print [0,scores[n,component]]
            ax_samp.plot(0, scores[n, component], 'D', markersize=15, color=colorRBGA, alpha=SAMPLE_ALPHA, label=type_label if type_label not in seen else "")
            seen.add(type_label)  #=======================================================================
        ax_samp.set_xlabel('Sample Scores\tFeature Scores')
        ax_samp.get_xaxis().set_ticks([])
        
        feat_scores = self.feature_scores
#         print "FEAT SCORES:",feat_scores.shape
        
        loadvals = [x for x in feat_scores[:, component]]
#         abs_loadvals = [abs(x) for x in loadvals]
        cm = plt.get_cmap('RdBu')
        color_range = [cm(x + .5)[:-1] for x in loadvals]
    #     print 'loadvals',len(loadvals),loadvals
    #     print 'features',len(features),features
    #     print 'colors',len(color_range),color_range
        feature_data = sorted(zip(self.feature_names, loadvals, color_range), key=lambda x: x[1])
        #=======================================================================
#         for i in range(len(feature_data)):
#             if abs(feature_data[i][1]) > cutoff: ##if this point beats cutoff
#                 ax.plot(feature_data[i][1],i,'s',color=feature_data[i][2])
#     #             if feature_data[i][1] > 0: 
#     #                 h_align = 'left' #vertical
#     #             else:
#     #                 h_align = 'right'
#     #             ax.text(feature_data[i][1],i, ' '+feature_data[i][0], rotation = 0, fontsize = '8', rotation_mode = 'anchor')
#                 ax.annotate(s = '  '+feature_data[i][0]+'  ', 
#                             xy = (feature_data[i][1],i), 
#                             horizontalalignment = 'left' if feature_data[i][1] > 0 else 'right', 
#                             verticalalignment = 'center',
#                             fontsize = '8')
#             else:
#                 ax.plot(feature_data[i][1],i,'o',color=feature_data[i][2])
#             mn = min(0,feature_data[i][1])
#             mx = max(0,feature_data[i][1])
#             ax.hlines(i,mn,mx,linestyles = 'dotted')
#         
#         ax.vlines(0,-1,len(features),linestyles = 'solid')
#         ax.set_yticks(range(len(features)))
#         ax.set_yticklabels([x[0] for x in feature_data], fontsize = '6')
        #=======================================================================
        SCALAR = 5
        for n in range(len(feature_data)):  # iterate through feature score positions
    #         print type_label, " : ", colorRBGA
#                ax.plot(sample_scores[0:5,n],np.ones(5)*n,'d', markersize=30, color='blue',\
#                      alpha=0.5, type_label='hESC')
#             print [0,scores[n,component]]
            this_feature = feature_data[n][0]
#             if this_feature in feature_names_to_plot:
#                 print this_feature,'DFADFDS', feature_data[n][1]
            ax_samp.plot(1, feature_data[n][1] * SCALAR, \
                         's' if this_feature in feature_names_to_plot else '.', \
                         color=feature_data[n][2], \
                         label=this_feature if this_feature in feature_names_to_plot else None)
            if this_feature in feature_names_to_plot:
                ax_samp.annotate(s='  ' + feature_data[n][0] + '  ',
                            xy=(1, feature_data[n][1] * SCALAR),
#                             horizontalalignment = 'left' if feature_data[i][1] > 0 else 'right', 
#                             verticalalignment = 'center',
                            fontsize='8')
#             ax_feat.plot([0,scores[n,component]],'D',markersize=15, color=colorRBGA, alpha=0.5, label = type_label if type_label not in seen else "")
#             seen.add(type_label)
        ax_samp.set_xbound((-0.5, 1.5))
        plt.savefig(self._avoid_overwrite(self.savedir + str(tag if tag else '') + 'boxplot_comp' + str(component) + '.pdf'), bbox_inches="tight")
        plt.legend(loc='best', bbox_to_anchor=(.85, .4, 1, 1))

        plt.clf()
#             print self.data[:,self.feature_names.index(feat)]
#         print self.data.shape
        
    def sample_name_to_type(self, sample_name):
        return sample_name.split("_")[0]

    def run_pca(self, select_types=None):
        data = self.data
        self.current_run_types = [self.types[x]['name'] for x in select_types]
        if select_types != None:
            select_indices = [self.types[x]['samp_indices'] for x in select_types]
            select_indices = [x for sublist in select_indices for x in sublist]
            data = self.data[select_indices, :]
            run_sample_names = np.array(self.sample_names)[select_indices]
            if self.verbose:
                print '\tselected types:'
                for x in select_types:
                    print '\t\t[', x, ']', self.types[x]['name'], ':', len(self.types[x]['samp_indices']), \
                    'samples at indices', self.types[x]['samp_indices']

        print '\tsample labels:', len(run_sample_names)  # ,sample_names
        print '\tdata', np.shape(data)
#         print data[:,[0,-1]]

        pca = skl_pca()
        pca = pca.fit(data)
        num_good_comps = len([x for x in pca.explained_variance_ratio_ if x > self.comp_var_threshold])
        
        # run again, only calculating/preserving components > threshold
        pca = skl_pca(n_components=num_good_comps+1)
        pca = pca.fit(data)
        
        feature_scores = pca.components_.T  # transform to get n_features rows by n_components columns
        sample_scores = pca.transform(data)
        comp_variance = pca.explained_variance_ratio_ 
        if self.verbose is True:
            print "\tsample_scores:", sample_scores.shape
            print "\tfeature_scores: ", feature_scores.shape
            print '\tnum components: ', num_good_comps, ' (threshold: ', self.comp_var_threshold, ')', comp_variance

        self.feature_scores = feature_scores
        self.sample_scores = sample_scores
        self.comp_var = comp_variance
        self.run_sample_names = run_sample_names
        self.run_types_by_num = select_types

        self.sample_scores_dict = {}
        for i, sample in enumerate(run_sample_names):
            self.sample_scores_dict[sample] = list(sample_scores[i])

        self.feature_scores_dict = {}
        for i, feat in enumerate(self.feature_names):
            self.feature_scores_dict[feat] = {}
            self.feature_scores_dict[feat]['scores'] = list(feature_scores[i])
            cm = plt.get_cmap('RdBu')
            pca_color_vals = [cm(x + 0.5)[:-1] for x in feature_scores[i]]  # '+0.5' centers the colormap around white
            if len(self.current_run_types) == 2:
                type_colors = [self.type_to_color(x) for x in self.current_run_types]
                cm = self._diverge_map(high=type_colors[1][:-1], low=type_colors[0][:-1])
                type_color_vals = [cm(x + 0.5)[:-1] for x in feature_scores[i]]
            else:
                # if not binary PCA, default to light gray
                type_color_vals = [(0.9, 0.9, 0.9, 0.25)] * len(pca.components_)
            self.feature_scores_dict[feat]['pca_raw_color'] = pca_color_vals
            self.feature_scores_dict[feat]['pca_type_color'] = type_color_vals

    def _make_colormap(self, seq):
        """Return a LinearSegmentedColormap
        seq: a sequence of floats and RGB-tuples. The floats should be increasing
        and in the interval (0,1).
        """
        seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
        cdict = {'red': [], 'green': [], 'blue': []}
        for i, item in enumerate(seq):
            if isinstance(item, float):
                r1, g1, b1 = seq[i - 1]
                r2, g2, b2 = seq[i + 1]
                cdict['red'].append([item, r1, r2])
                cdict['green'].append([item, g1, g2])
                cdict['blue'].append([item, b1, b2])
        return mcolors.LinearSegmentedColormap('CustomMap', cdict)

    def _diverge_map(self, high=(0.565, 0.392, 0.173), low=(0.094, 0.310, 0.635)):
        '''
        low and high are colors that will be used for the two
        ends of the spectrum. they can be either color strings
        or rgb color tuples
        '''
        c = mcolors.ColorConverter().to_rgb
        if isinstance(low, basestring): low = c(low)
        if isinstance(high, basestring): high = c(high)
        return self._make_colormap([low, c('white'), 0.5, c('white'), high])

    def boxplot_single1(self, component=0, number_boxplots=3, tag=None):
        number_boxplots = 1
        SAMPLE_ALPHA = 1.0
        print 'BOXPLOT:'
        print '\tcomponent:', component
        np_feat_names = np.array([self.feature_names]).T
        comp_score = np.array([self.feature_scores[:, component]]).T
#         print comp_score.shape
#         print np_feat_names.shape
        feature_names_comp = np.concatenate((np_feat_names, comp_score), axis=1)
#         print self.feature_names[0:3]
#         print self.feature_scores[0:3]
#         print comp_score[0:3]
        a = feature_names_comp.tolist()
#         for i in a:
#             print i
        a = sorted(a, key=lambda x: abs(float(x[1])), reverse=True)
#         for i in a:
#             print i
        feature_names_to_plot = [x[0] for x in a[0:number_boxplots]]
#         print feature_names_to_plot
        
        fig, axarr = plt.subplots(number_boxplots, 2, sharex=True, sharey=False)
        fig.suptitle('boxplot comp' + str(component) + ' significant features')

        for j, feat in enumerate(feature_names_to_plot):
            boxplot_data = []
            feat_index = self.feature_names.index(feat)
            type_data_dict = {}
            orderseen = []
            for i, sample in enumerate(self.run_sample_names):
                sample_type = self.sample_name_to_type(sample)
                if sample_type in type_data_dict.keys():
                    type_data_dict[sample_type] += [self.data[i][feat_index]]
                else:
                    type_data_dict[sample_type] = [self.data[i][feat_index]]
                    orderseen += [sample_type]
            print feat,
            print feat_index
            print type_data_dict.items()
            print orderseen
            print '++++\n'
            for typesamp in orderseen:
                boxplot_data += [type_data_dict[typesamp]]
            ax = axarr[j][0]
            ax.boxplot(boxplot_data)
            ax.set_xticklabels(orderseen)
            ax.get_xaxis().tick_bottom()  # remove top ticks
            ax.get_yaxis().tick_left()  # remove right ticks
            ax.set_title(feat)
            axarr[j][1].set_axis_off()

        # ADD COMPONENT SAMPLE PCA SCORES, COLORBAR
        #=======================================================================
        ax_samp = fig.add_subplot(1, 2, 2)

        scores = self.sample_scores
        seen = set()
        print "SAMPLE SCORES:", scores.shape
        for n in range(scores.shape[0]):  # iterate through sample score positions
            type_label = self.sample_name_to_type(self.run_sample_names[n])
            colorRBGA = self.type_to_color(type_label)
    #         print type_label, " : ", colorRBGA
#                ax.plot(sample_scores[0:5,n],np.ones(5)*n,'d', markersize=30, color='blue',\
#                      alpha=0.5, type_label='hESC')
#             print [0,scores[n,component]]
            ax_samp.plot(0, scores[n, component], 'D', markersize=15, color=colorRBGA, alpha=SAMPLE_ALPHA, label=type_label + '000' if type_label not in seen else "")
            seen.add(type_label)  #=======================================================================
        ax_samp.set_xlabel('Sample Scores\tFeature Scores')
        ax_samp.get_xaxis().set_ticks([])
        
        feat_scores = self.feature_scores
#         print "FEAT SCORES:",feat_scores.shape
        
        loadvals = [x for x in feat_scores[:, component]]
#         abs_loadvals = [abs(x) for x in loadvals]
        cm = plt.get_cmap('RdBu')
        color_range = [cm(x + .5)[:-1] for x in loadvals]
    #     print 'loadvals',len(loadvals),loadvals
    #     print 'features',len(features),features
    #     print 'colors',len(color_range),color_range
        feature_data = sorted(zip(self.feature_names, loadvals, color_range), key=lambda x: x[1])
        SCALAR = 5
        for n in range(len(feature_data)):  # iterate through feature score positions
    #         print type_label, " : ", colorRBGA
#                ax.plot(sample_scores[0:5,n],np.ones(5)*n,'d', markersize=30, color='blue',\
#                      alpha=0.5, type_label='hESC')
#             print [0,scores[n,component]]
            this_feature = feature_data[n][0]
#             if this_feature in feature_names_to_plot:
#                 print this_feature,'DFADFDS', feature_data[n][1]
            ax_samp.plot(1, feature_data[n][1] * SCALAR, \
                         's' if this_feature in feature_names_to_plot else '.', \
                         color=feature_data[n][2], \
                         label=this_feature if this_feature in feature_names_to_plot else None)
            if this_feature in feature_names_to_plot:
                ax_samp.annotate(s='  ' + feature_data[n][0] + '  ',
                            xy=(1, feature_data[n][1] * SCALAR),
#                             horizontalalignment = 'left' if feature_data[i][1] > 0 else 'right', 
#                             verticalalignment = 'center',
                            fontsize='8')
#             ax_feat.plot([0,scores[n,component]],'D',markersize=15, color=colorRBGA, alpha=0.5, label = type_label if type_label not in seen else "")
#             seen.add(type_label)
        ax_samp.set_xbound((-0.5, 1.5))
        plt.savefig(self._avoid_overwrite(self.savedir + tag + 'boxplot_comp' + str(component) + '.pdf'), bbox_inches="tight")
        plt.legend(loc='best', bbox_to_anchor=(.85, .4, 1, 1))

        plt.clf()

    def _avoid_overwrite(self, save_location):
        inc = 1
        while path.isfile(save_location):
            save_location = save_location.rstrip('.pdf')
            if inc != 1:
                save_location = save_location.rsplit('_', 1)[0]
            save_location = save_location + '_' + str(inc) + '.pdf'
            inc += 1
        return save_location
    
    def _set_global_max_min(self, data):
        return min(data.flatten()), max(data.flatten())

    def __init__(self, savedir, data, rawdata, feature_names, sample_names, types, \
                 select_types=None, comp_var_threshold=0.1, verbose=True):
        '''
        '''
        self.savedir = savedir
        self.feature_names = feature_names
        self.sample_names = sample_names
        self.data = data
        self.rawdata = rawdata
        self.types = types
        self.min_data, self.max_data = self._set_global_max_min(self.data)
        self.min_rawdata, self.max_rawdata = self._set_global_max_min(self.rawdata)
        
        self.verbose = verbose
        self.comp_var_threshold = comp_var_threshold
#         self.run_pca(select_types, max_comps, comp_var_threshold, verbose)
        if verbose:
            print 'PCA: \n\tinitialized'
        
# TODO: switch all plotting methods from using self.feature_scores/sample_scores to using\
        # TODO:  self.feature_scores_dict/sample_scores_dict
        # TODO: ^
        # TODO:^
            

#         print feature_scores.shape
#         feature_names = np.array([feature_names]).T
#         self.np_features_name_scores = np.concatenate((feature_names, feature_scores), axis = 1)
#         print feature_names.shape
#         print feature_names[:3]
#         print self.np_features_name_scores[:3,:]
#         
#         
#         sample_names = np.array([sample_names]).T
#         self.np_samples_name_scores = np.concatenate((sample_names, sample_scores), axis = 1)
#         print self.np_samples_name_scores[:3,:]

        
        
        
        
        
        
        
