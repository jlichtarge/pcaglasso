'''
Created on Jun 11, 2016

@author: Jared
'''
import time
import sys
import os
import numpy as np
from sklearn.preprocessing import scale
from matplotlib.pyplot import get_cmap 


class dataio(object):
    '''bundles I/O functions for readability
        instance vars:
            savedir: string
                output directory path
            data: numpy.ndarray
                (nsamples x nfeatures) of standardized input data
            feature_names: list
                list of feature names (nfeatures long)
                features_names[k] corresponds to k-th column of data: data[:][k]
            sample_names: list
                list of sample names (nsamples long)
                samples_names[j] corresponds to j-th row of data: data[j][:]
    '''
    
    #potential cm options: spectral, viridis, plamsa, magma
    type_colormap = 'spectral'

    def generate_savedir(self):
        '''generates directory for run output files
        
        output/ folder is created in same dir as main.py
        further nesting folders specify day, specific time of run
        
        Example:
            perform run on 5:35:03 PM, January 22, 2017 -->
                files saved to /output/Jan_22_17/17_34_03/
        '''
        timestamp_day = time.strftime('%b_%d_%y')
        timestamp_run = time.strftime('%H_%M_%S')
        output_dir = sys.path[0]+'/output/'
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        output_dir = sys.path[0]+'/output/'+timestamp_day+'/'
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        output_dir += timestamp_run + '/'
        os.mkdir(output_dir)
        return output_dir
    
    def _open_file(self):
        """helper for data_from_file:
        open self.file and return its contents, split by line 
        """
        f = open(self.file, 'r')
        file_lines = f.read().splitlines()
        f.close()
        return file_lines
    
    def data_from_file(self):
        """read in data from self.file
        
            Input File Specs
            ----------
                split file by newline
                remove comments (lines starting with #)
                parse data
                    first row is a whitespace-separated list of feature labels (labeling the columns beneath them)
                    the remaining rows correspond to each sample:
                        first element is string representing sample name for current row
                            each sample label is a concatenation of the TYPE of the sample, 
                            and a unique identifier, separated by an underscore (eg. type1_sampleName)
                        remaining elements in row are data points
            Returns
            ----------
                sample_lbls : list
                    labels for each sample [row]
                
                feature_lbls : list
                    labels corresponding to each feature [column]
          
                rawdata : np.ndarray
                    raw data matrix
                
            Notes
            ----------
                sample_lbls and feature_lbls retain order such that the index of an element
                    in each corresponds to that row (for samples) or column (for features)
                    of the output data matrix
           
            Example
            ----------
                input file:
                         ""feat1    feat2    feat3
                    t1_s1  5        1        8
                    t1_s2  6        2        7                    
                    t2_s3  4        1        6
                    t2_s4  5        1        8""
                    
                outputs:
                                        
                    sample_lbls:
                        ['t1_s1', 't1_s2', 't2_s3', 't2_s4']
                    
                    feature_lbls:
                        ['feat1', 'feat2', 'feat3']

                    rawdata:
                        [[ 5.  1.  8.]
                         [ 6.  2.  7.]
                         [ 4.  1.  6.]
                         [ 5.  1.  8.]]
                    
        """
        #open file, split into lines
        file_lines = self._open_file()
        
        #remove comments
        file_lines = [x for x in file_lines if x and x.strip()[0]!='#']
        #pop first line to retrieve sample names, split by space
        feature_lbls = file_lines.pop(0).strip().split()

        #parse body of data: feature_name \t data.1 \t data.2 ... data.n_samples
        sample_lbls, rawdata = [], []
        for line in [l.split() for l in file_lines]:
            sample_lbls.append(line.pop(0))
            #check correct number of data points
            assert (len(line) == len(feature_lbls)),"Data Format Error: incorrect number of samples for "+\
                sample_lbls.pop()+":"+str(len(line))+"!="+str(len(feature_lbls))
            rawdata.append(map(np.float64,line))
        rawdata = np.array(rawdata)   
        
        if self.verbose:
            print '\tdata:',rawdata.shape
            print '\tsamples:',len(sample_lbls),'\t', sample_lbls[0],'...',sample_lbls[-1]
            print '\tfeatures:',len(feature_lbls),'\t',feature_lbls[0],'...',feature_lbls[-1]  
        
        return sample_lbls, feature_lbls, rawdata

    def types_to_indices_dictionary(self):
        '''extracts types from sample_names in order of appearance
        returns a list of type dictionaries preserving that order
            
            Returns
            ----------
            types : list of types, in order of appearance in sample_names
                each element in the list is a dictionary corresponding to that type
                    dictionary keys : 
                        'name' : string
                            the label of the type
                        'samp_indices' : list
                            list of sample_names indices whose type is that of the key,
                            also corresponds to rows in self.data of the type specified
                            by the same key
                        'color' : tuple
                            RGBA color tuple (RED, GREEN, BLUE, ALPHA)
                            generated by _init_colors_for_types
                    
            Example
            ----------
            input file :
                ['fea1', 'fea2', 'fea3']
                A_s1 [ 1.  2.  3.]        <-- first type A  -->     types[0]['name'] == A
                B_s2 [ 4.  5.  6.]        <-- second type B -->     types[1]['name'] == B
                B_s3 [ 7.  8.  9.]
                A_s4 [ 10.  11.  12.]
                C_s5 [ 13.  14.  15.]     <-- third type C  -->     types[2]['name'] == C
            
            temp_type_dict :
                dict{
                    A {'indices': [0, 3]}
                    B {'indices': [1, 2]} 
                    C {'indices': [4]}
                }
            
            types :
                list(
                    dict{
                        'name' : 'A'
                        'samp_indices' : [0, 3]
                        'color' : (1,0,0,1)
                    },
                    dict{
                        'name' : 'B'
                        'samp_indices' : [1, 2]
                        'color' : (0,1,0,1)
                    },
                    dict{
                        'name' : 'C'
                        'samp_indices' : [4]
                        'color' : (0,,1,1)
                    }
                )
            
            
            Notes
            ----------
            the indices specified by temp_type_dict map to both the locations within 
                self.sample_names and the corresponding data points within self.data
            
            each type in the list is specified by a dictionary so that multiple attributes 
                may be associated with each type (currently just indices and color)
                
            it is necessary to return a list of dictionaries (types) rather than simply
                the intermediary dictionary temp_type_dict, in order to preserve the 
                order of appearance of the types, which will be used to specify
                type comparisons in PCA analysis
        '''
        temp_type_dict = {}
        for i, sample in enumerate(self.sample_names):
            typ = sample.split("_")[0]
            if not temp_type_dict.__contains__(typ):
                temp_type_dict[typ] = {}
                temp_type_dict[typ]['samp_indices'] = []
            temp_type_dict[typ]['samp_indices'].append(i)

        list_of_types = self._remove_duplicates_inorder([x.split('_')[0] for x in self.sample_names])
        types = []
        for i,t in enumerate(list_of_types):
            types.append({})
            types[i]['name'] = t
            types[i]['samp_indices'] = temp_type_dict[t]['samp_indices']

        types = self._init_colors_for_types(types)
        if self.verbose:
            print '\tsplitting samples into',len(types),'types:'
            for i, t in enumerate(types):
                print '\t  [',i,']', t['name'],':'
                print '\t\t  samp_indices:',len(t['samp_indices']),'samples, indices',t['samp_indices']
                print '\t\t  color: (%1.2f, %1.2f, %1.2f)' % (t['color'][0], t['color'][1],t['color'][2])
                
        return types

    def _init_colors_for_types(self, types):
        '''associates each type in types with a color value for visualization
        '''
        cm = get_cmap(self.type_colormap)

        label_colors_list = [(float(i) +.5) / len(types) for i in range(len(types))]
        label_colors_list = [cm(x) for x in label_colors_list]
        
        for i, t in enumerate(types):
            t['color'] = label_colors_list[i]
        
        return types

    def _remove_duplicates_inorder(self,seq):
        '''remove duplicates from list, preserving order
        from: http://www.peterbe.com/plog/uniqifiers-benchmark and comments there
        '''
        seen = {}
        result = []
        for item in seq:
            if not seen.__contains__(item):
                seen[item] = 1
                result.append(item)
        return result
    
    def types_to_avg_feature_value_dictionary(self):
        '''TODO: decide to store these values in the dict or generate them later
        '''
        for i, tdict in enumerate(self.types):
            tdict['feat_avg'] = np.mean(self.data[self.types[i]['samp_indices'],], axis=0)
            tdict['feat_avg_raw'] = np.mean(self.rawdata[self.types[i]['samp_indices'],], axis=0)
    
    def __init__(self, inputfile, use_types = True, verbose = True):
        '''sets vars: savedir & (from inputfile:) data, feature_names, sample_names
        '''
        np.set_printoptions(precision=2, threshold = 50)
        
        if verbose: 
            print "\nINPUT: ", inputfile
        
        self.verbose = verbose
        self.savedir = self.generate_savedir()
        self.file = inputfile
        self.sample_names, self.feature_names, self.rawdata = self.data_from_file()
        self.data = scale(self.rawdata, axis=1, with_mean=True, with_std=True)
        self.types = self.types_to_indices_dictionary()
        self.types_to_avg_feature_value_dictionary()
        
        
        
        
        