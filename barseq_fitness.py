#%%
from operator import index
import sys
import pandas as pd
import numpy as np
import argparse as args
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import axes3d
import seaborn as sns

#%%
class LoadData:
    def __init__(self, file):
        self.file = file
        self.get_header()
        self.read_in()
        
    def read_in(self):
        with open(self.file) as f:
            header = f.readlines()[0].strip('\n').split('\t')
        self.df = pd.read_csv(self.file,sep='\t',skiprows=1,names=header)

    def get_header(self):
        with open(self.file) as f:
            header = f.readlines()[0].strip('\n').split('\t')
        self.header = header

class HandleGeneCounts(LoadData):
    def __init__(self,gene_counts_file):
        super().__init__(gene_counts_file)
        self.streamline_sample_names()
    
    def streamline_sample_names(self):
        # get only the set#INDEX eg. set1IT002
        self.header = [s.split()[0] for s in self.header]
        self.df.columns=self.header
    
class HandleLogRations(HandleGeneCounts):
    # if handling logratios_unnormalized file, it should parse the same way as gene_counts.tab
    def __init__(self, gene_counts_file):
        super().__init__(gene_counts_file)

class HandleExps(LoadData):
    def __init__(self, exps_file):
        super().__init__(exps_file)
        self.define_set()
        self.clean_conditional_nans()
        self.block_by_condition()
    
    def define_set(self):
        # this is to only get the set# to be fused with index
        # this makes parsing gene_counts.tab easier
        set_num = [s.split('_')[-1] for s in self.df['SetName']]
        self.df['SetNumber'] = set_num 
    
    def clean_conditional_nans(self):
        self.df['Condition_1'].fillna('.', inplace=True)
        self.df['Concentration_1'].fillna('.', inplace=True)
        self.df['Units_1'].fillna('.', inplace=True)
        self.df['Condition_2'].fillna('.', inplace=True)
        self.df['Concentration_2'].fillna('.', inplace=True)
        self.df['Units_2'].fillna('.', inplace=True)

    def block_by_condition(self):
        # subset the dataframe by conditions
        # hard coding for now, but may be sufficient to keep this way
        # need a way to detect if some fields are empty unless experiments are kept constant and we always know there are only 2 fields
        # return set/index names based on condition
        meta_conditions = {}
        conditions = self.df[['Condition_1','Condition_2']].drop_duplicates().to_numpy()
        for c in conditions:

            c1 = c[0]
            c2 = c[1]
            where_condition1 = np.where(c1 == self.df['Condition_1'].to_numpy())[0]
            where_condition2 = np.where(c2 == self.df['Condition_2'].to_numpy())[0]
            loci_list = np.intersect1d(where_condition1,where_condition2)
            set_numbers = self.df['SetNumber'].iloc[loci_list].to_numpy()
            index_list = self.df['Index'].iloc[loci_list].to_numpy()
            concentration1 = self.df['Concentration_1'].iloc[loci_list].to_numpy()
            concentration2 = self.df['Concentration_2'].iloc[loci_list].to_numpy()
            units1_list = self.df['Units_1'].iloc[loci_list].to_numpy()
            units2_list = self.df['Units_2'].iloc[loci_list].to_numpy()
            if len(np.unique(units1_list))>1:
                print('Warning: Condition 1 has multiple units: {}'.format(np.unique(concentration1)))
            if len(np.unique(units2_list))>1:
                print('Warning: Condition 2 has multiple units: {}'.format(np.unique(concentration2)))


            meta_conditions[tuple(c)] = {'exps_loci':loci_list,
                                         'condition1':c1,
                                         'concentration1':concentration1,
                                         'units1':units1_list,
                                         'condition2':c2,
                                         'concentration2':concentration2,
                                         'units2':units2_list,
                                         'set':set_numbers,
                                         'index':index_list,
                                         'set_index':[s+index_list[i] for i,s in enumerate(set_numbers)]
                                         }
        self.meta_conditions = meta_conditions
        # add function to convert mM to M and other conc stuff
        # in general input should all be in the same units anyways
        
        
class Plotting():
    def __init__(self,df,conditions,meta_data):
        # meta_data is a cut of the meta_conditions dictionary with conditions of interest only
        # this dataframe should be any datafile from rbseq pipeline
        # gene_counts.tab, or fit_logrations_unnormalized.tab
        # df should be dataframe parsed for pre parsed for plotting
        #   this also means data should be transformed or manipulated before being plotted
        # hand plot functions proper headers to plot the df
        self.df = df 
        self.conditions = conditions
        self.meta = meta_data
        conc1 = self.meta['concentration1']
        conc2 = self.meta['concentration2']
        sample = self.meta['set_index']
        conc_pairs = np.array(list(zip(conc1,conc2,sample)),dtype=object)
        self.to_plot_df = self.df[conc_pairs[:,2]].transpose()
        self.to_plot_df.columns = self.df['sysName'] # choose a better gene name
           
        
    def condition_gene(self,cond_number,quant):
        # cond = 1 or 2
        # quant is gene_count or fitness, ... [depends on file input]
        # 2d line plot of condition concentration vs gene abundance
        # need to fix if condition 1 because [0,0,0,0...,1,1,1,1...,2,2,2...]
        conc1 = self.meta['concentration1']
        conc2 = self.meta['concentration2']
        sample = self.meta['set_index']
        conc_pairs = np.array(list(zip(conc1,conc2,sample)),dtype=object)
    
        if cond_number == 1:
            conc_pairs=conc_pairs[np.lexsort((conc_pairs[:,0], conc_pairs[:,1]))]
            # track changes in column 2 of conc_pairs to make blocks
            blocks_mark = set(conc_pairs[:,1])
            block_index = [np.where(b==conc_pairs[:,1])[0] for b in blocks_mark]
            plot_conc = np.sort(np.unique(conc_pairs[:,0]))


        if cond_number == 2:
            conc_pairs=conc_pairs[np.lexsort((conc_pairs[:,1], conc_pairs[:,0]))]
            # track changes in column 1 of conc_pairs to make blocks
            blocks_mark = set(conc_pairs[:,0])
            block_index = [np.where(b==conc_pairs[:,0])[0] for b in blocks_mark]
            plot_conc = np.sort(np.unique(conc_pairs[:,1]))
            #print(block_index)

        condition_name = self.conditions[cond_number-1]
  
        # find columns to plot from meta
        units = self.meta['units{}'.format(cond_number)]
        xlab = '{} {}'.format(condition_name, units[0])
        plt.figure(figsize=[10,10])
        plt.title('{} vs {} Single Condition Visualization'.format(xlab,quant))
        plt.xlabel(xlab)
        plt.ylabel(quant)
        for block in block_index:
            columns = conc_pairs[block,2]
            self.to_plot_df = self.df[columns].transpose()
            self.to_plot_df.columns = self.df['sysName'] # choose a better gene name
            plot_matrix = self.to_plot_df.to_numpy() # numpy plots faster
            # plot matrix dims:
                # rows are the index
                # columns are the gene
            for i,col in enumerate(self.to_plot_df.columns):
                # at some point add a legend, label=col (gene name)
                plt.plot(plot_conc, plot_matrix[:,i],'o-',alpha=.02)
                pass

        #plt.legend()
        
        plt.show()
        plt.close()
    
    def gene_isocline(self,gene):
        plt.plot(self.to_plot_df[gene])
        plt.show()
        plt.close()


    def condition1_condition2_gene(self,gene):
        # 3d plot of condition1 vs condition2 vs gene abundance
        conc1 = self.meta['concentration1']
        conc2 = self.meta['concentration2']
        sample = self.meta['set_index']
        conc_pairs = np.array(list(zip(conc1,conc2,sample)),dtype=object)

        fig = plt.figure(figsize=[10,10])
        ax = fig.add_subplot(121, projection='3d')
        gene_df_index = np.where(gene==self.df['sysName'])[0]
        x = conc_pairs[:,0]
        y = conc_pairs[:,1]
        top = self.df[conc_pairs[:,2]].iloc[gene_df_index].to_numpy()[0]
        bottom = np.zeros((len(x)))
        width=.2
        depth=.2

        ax.bar3d(x, y, bottom, width, depth, top, shade=True,alpha=.2)
        ax.set_xlabel('{} {}'.format(self.meta['condition1'], self.meta['units1'][0]))
        ax.set_ylabel('{} {}'.format(self.meta['condition2'], self.meta['units2'][0]))
        plt.show()

    def calculate_zmap(self,gene):
        # change this to calculate (8,12) shape
        conc1 = self.meta['concentration1']
        conc2 = self.meta['concentration2']
        sample = self.meta['set_index']
        conc_pairs = np.array(list(zip(conc1,conc2,sample)),dtype=object)
        a = np.unique(np.sort(self.meta['concentration1'])) # len 8
        b = np.unique(np.sort(self.meta['concentration2'])) # len 12
        x,y = np.meshgrid(a,b) # x is shape (12,8)
        print(x.shape,y.shape)
        # x--> conc1:
            # [0 .013 .026 .039 ...]
            # [0 .013 .026 .039 ...]
        # y--> conc2
            # [0  0  0  0  0]
            # [.0138 .0138 .0138 .0138]
        z = np.zeros((x.shape[0],x.shape[1]))
    
        for i in range(x.shape[0]):
            for j in range(x.shape[1]):
                find_z = np.intersect1d(np.where(x[i,j]==conc_pairs[:,0])[0],np.where(y[i,j]==conc_pairs[:,1])[0])[0]
                sample = conc_pairs[find_z,2]
                z[i,j] = self.df[sample].iloc[np.where(gene==self.df['sysName'])[0]]
        #z = np.flip(z)
        return x,y,z

    def wireframe(self,gene):
        x,y,z = self.calculate_zmap(gene)
        print(x.shape)
        print(y.shape)
        # x--> conc1:
            # [0 .013 .026 .039 ...]
            # [0 .013 .026 .039 ...]
        # y--> conc2
            # [0  0  0  0  0]
            # [.0138 .0138 .0138 .0138]
        #z = np.zeros((x.shape[0],x.shape[1]))
    

        fig = plt.figure(figsize=[15,15])
        wf = fig.add_subplot(projection ='3d')
        wf.plot_wireframe(x, y, z, color ='blue')

        gene_df_index = np.where(gene==self.df['sysName'])[0]
        #x_ = conc_pairs[:,0]
        #y_ = conc_pairs[:,1]
        #top = self.df[conc_pairs[:,2]].iloc[gene_df_index].to_numpy()[0]
 
        #bottom = np.zeros((len(x_)))
        #width=.5
        #depth=.5
        wf.set_title('3D wireframe barplot {} abundance \n dummy data'.format(gene))
        #wf.bar3d(x_, y_, bottom, width, depth, top, shade=True, alpha=.15,color='cyan')
        

        wf.set_xlabel(self.meta['condition1'])
        wf.set_ylabel(self.meta['condition2'])
        wf.set_xlabel('{} {}'.format(self.meta['condition1'], self.meta['units1'][0]),fontsize=10)
        wf.set_ylabel('{} {}'.format(self.meta['condition2'], self.meta['units2'][0]),fontsize=10)
        wf.set_zlabel('Gene {} Abundance'.format(gene))
        plt.show()
        plt.close()
    
    def heatmap_condition1_condition2_gene(self,gene):
        x,y,z = self.calculate_zmap(gene)
        fig,ax = plt.subplots(figsize=(8,7))
        ax = sns.heatmap(z,
                        cmap=plt.get_cmap('cividis'),
                        cbar_kws={'label': 'power'})
                        #cividis is a colorblind friendly palette
        title = '{}'.format(gene)
        xlab = '{} ({})'.format(self.meta['condition1'], self.meta['units1'][0])
        ylab = '{} ({})'.format(self.meta['condition2'], self.meta['units2'][0])
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(title)

        plt.show()
        #plt.savefig(title+'.png')
        plt.close()

#%%
gene_counts_file = 'gene_counts.tab'
exps_file = 'exps'
GenesCounts = HandleGeneCounts(gene_counts_file)
Experiments = HandleExps(exps_file) 

genes_counts_df = GenesCounts.df 
conditions = ('Phosphate','Nitrate')
section = Experiments.meta_conditions[conditions]
plot_df = GenesCounts.df 
#plot_df.iloc[:,4:] = np.log(plot_df.iloc[:,4:])
P = Plotting(plot_df, conditions, section)
condition1_mat, condition2_mat, z = P.calculate_zmap('b0001')
#P.condition_gene(cond_number=2,quant='gene counts')
#P.condition_gene(cond_number=1,quant='gene counts')
P.condition1_condition2_gene('b0001')
P.wireframe('b0001')
P.heatmap_condition1_condition2_gene('b0001')
P.gene_isocline('b0001')
'''
Useful headers for exps file:
SetName and Index for sample names
Condition_#, Concentration_#, Units_#
Temperature
Description
StartOD, EndOD
Total Generations

Note: When making the exps file, for SetName use the naming scheme:
        [Experiment_name]_[set_#]
^ the name will be fused with index [IT###]
--> in gene_counts.tab the sample names become: [set{}IT### condition]
# Assumption that you are running this program with a rbseq file of only a single 
  96 well sequencing run
'''


#%%
# exps file has information on the set names and experimental conditions


# %%
