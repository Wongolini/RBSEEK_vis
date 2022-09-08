#%%
import sys
import pandas as pd
import numpy as np
import argparse as args
import matplotlib.pyplot as plt 


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
                                         'concentration1':concentration1,
                                         'units1':units1_list,
                                         'concentration2':concentration2,
                                         'units2':units2_list,
                                         'set':set_numbers,
                                         'index':index_list,
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
        
    def condition_gene(self,cond_number,quant):
        # cond = 1 or 2
        # quant is gene_count or fitness, ... [depends on file input]
        # 2d line plot of condition concentration vs gene abundance
        
        condition_name = self.conditions[cond_number-1]
        concentration = 'concentration{}'.format(cond_number)
        print(condition_name)
        # find columns to plot from meta
        units = self.meta['units{}'.format(cond_number)]
        xlab = '{} {}'.format(condition_name, units[0])
        columns = []
        plt.figure(figsize=[10,10])
        plt.xlabel(xlab)
        plt.ylabel(quant)
    
        for set_num,index in zip(self.meta['set'],self.meta['index']):
            sample = set_num+index
            columns.append(sample)
        self.to_plot_df = self.df[columns].transpose()
        self.to_plot_df.columns = self.df['sysName']
        self.to_plot_df['concentration'] = self.meta[concentration]
        plot_matrix = self.to_plot_df.to_numpy()

        # concentrations will repeat if plotting on one axis
        single_conc = np.unique(self.meta[concentration])
        repeats = len(single_conc)
        repeats_list = [int(r*repeats) for r in range(len(single_conc)+1)]

        for i,col in enumerate(self.to_plot_df.columns[:-1]):
            for j in range(1,len(repeats_list)):
                plt.plot(self.to_plot_df['concentration'][repeats_list[j-1]:repeats_list[j]], 
                                        plot_matrix[repeats_list[j-1]:repeats_list[j],i],'o-',alpha=.5,linewidth=.2,label=col)
        plt.show()
        plt.close()
       

    def condition1_condition2_gene():
        # 3d plot of condition1 vs condition2 vs gene abundance
        pass 


gene_counts_file = 'gene_counts.tab'
exps_file = 'exps'
GenesCounts = HandleGeneCounts(gene_counts_file)
Experiments = HandleExps(exps_file) 

genes_counts_df = GenesCounts.df 
section = Experiments.meta_conditions[('Phosphate','Nitrate')]
to_plot_df = genes_counts_df.iloc[section['exps_loci']]
to_plot_df = to_plot_df.iloc[[x for x in range(len(to_plot_df)) if np.all(to_plot_df.iloc[:,4:].to_numpy()[x]!=0)]]
P = Plotting(GenesCounts.df, ('Phosphate','Nitrate'), section)
P.condition_gene(cond_number=2,quant='counts')

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
'''


#%%
# exps file has information on the set names and experimental conditions



# %%
