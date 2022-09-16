#!/usr/bin/env python
import pandas as pd 
import numpy as np
import string
import argparse
import os, sys
'''
read in a config file for the experiment
config file explains condition 1 and 2
eg.
## SetName =  
## Date_pool_expt_started =
## Person =
## Mutant Library =
## Description = 
## gDNA plate = 
## gDNA well = 
## Sequenced At = 
## Media = 
## Growth Method = 
## Group = 
## Temperature =
## pH =
## Liquid v. solid =
## Aerobic_v_Anaerobic =  
## Shaking = 
## Condition 1 = 
## Concentration 1 low = 
## Concentration 1 high = 
## Concentration 1 n_increments = 
## Units 1 =
## Concentration 2 low =
## Concentration 2 high =
## Concentration 2 n_ncrements =
## Units 2 =
## Timecourse = 
## Timecourse Sample = 
## Growth Plate ID = 
## Growth Plate wells = 

config file will allow this program to write Indexes and what concnetration describes what well in sequencing
output will be a 96 well table with columns and rows detailing condition1 and condition2 info useful for dna extraction and library prep
condition 1 must be max 8 increments (rows)
condition 2 must be max 12 increments (columns)
# '''
class ReadConfig():
    def __init__(self,config_file):
        self.config_file = config_file
        self.read_config()

    def read_config(self):
        config_dict = {}
        with open(self.config_file) as f:
            for line in f.readlines():
                if line.startswith('##'):
                    try:
                        line=line.lstrip('## ').rstrip('\n')
                        line = line.split(' = ')
                        variable=line[0]
                        value=line[1]
                        config_dict[variable] = value 
                    except:
                        continue
                else:
                    pass 
                
        self.config_dict = config_dict
    
class DefineExperimentMeta():
    def __init__(self,config_dict):
        self.config_dict = config_dict
        self.create_conditions_matrix()
        self.header = ['SetName', 'Date_pool_expt_started', 'Drop', 
                        'Person', 'Mutant Library', 'Description', 
                        'gDNA plate', 'gDNA well', 'Index', 'Sequenced At', 
                        'Media', 'Growth Method', 'Group', 'Temperature', 'pH', 
                        'Liquid v. solid', 'Aerobic_v_Anaerobic', 'Shaking', 
                        'Condition_1', 'Concentration_1', 'Units_1', 
                        'Condition_2', 'Concentration_2', 'Units_2', 
                        'Timecourse', 'Timecourse Sample', 'Growth Plate ID', 
                        'Growth Plate wells', 'StartOD', 'EndOD', 
                        'Total Generations']
        self.df = pd.DataFrame(columns=self.header)
        self.write_table()
    
    def create_conditions_matrix(self):
        self.conditions_dict = {}
        condition1 = self.config_dict['Condition_1']
        concentration1_low = float(self.config_dict['Concentration 1 low'])
        concentration1_high = float(self.config_dict['Concentration 1 high'])
        concentration1_incr = int(self.config_dict['Concentration 1 n_increments'])
        condition1_units = self.config_dict['Units_1']
        condition2 = self.config_dict['Condition_2']
        concentration2_low = float(self.config_dict['Concentration 2 low'])
        concentration2_high = float(self.config_dict['Concentration 2 high'])
        condition2_units = self.config_dict['Units_2']
        concentration2_incr = int(self.config_dict['Concentration 2 n_increments'])

        # define the microplate
        row_letters = list(string.ascii_uppercase[:8])
        col_numbers = list(range(1,13))

        # number of rows == concentration1 n_increments
        # number of cols == concentration2 n_increments
        condition1_mat = np.zeros((concentration1_incr,concentration2_incr))
        condition2_mat = np.zeros((concentration1_incr,concentration2_incr))

        condition1_vector = [x*(concentration1_high-concentration1_low)/(concentration1_incr-1) for x in range(concentration1_incr)]
        # condition1_vector size 8
        condition2_vector = [x*(concentration2_high-concentration2_low)/(concentration2_incr-1) for x in range(concentration2_incr)]
        # condition2_vector size 12
        # x--> conc1 (8,12)
            # [0  0  0  0  0]
            # [.0138 .0138 .0138 .0138]

        # y--> conc2: (8,12)
            # [0 .013 .026 .039 ...]
            # [0 .013 .026 .039 ...]

        for col in range(concentration2_incr):
            condition1_mat[:,col] = condition1_vector
        for row in range(concentration1_incr):
            condition2_mat[row] = condition2_vector

        for col in range(concentration2_incr):
            for row in range(concentration1_incr):
                well_key = '{}{}'.format(row_letters[row],col_numbers[col])
                self.conditions_dict[well_key] = { 
                                                  condition1:condition1_mat[row,col],
                                                  'Units_1':condition1_units,
                                                  condition2:condition2_mat[row,col],
                                                  'Units_2':condition2_units
                                                 }
    
    def write_table(self):
        well_list = []
        condition1_list = []
        condition2_list = []
        index_list = []
        for i,well in enumerate(self.conditions_dict.keys()):
            well_list.append(well)
            condition1_list.append(self.conditions_dict[well][self.config_dict['Condition_1']])
            condition2_list.append(self.conditions_dict[well][self.config_dict['Condition_2']])
            index_str = str(i+1)
            if len(index_str)==1:
                number='00'+index_str
            elif len(index_str)==2:
                number='0'+index_str
            elif len(index_str)==3:
                number=index_str
            index = 'IT'+number
            index_list.append(index)
        self.df['Index'] = index_list

        for h in self.header:
            try:
                self.df[h] = [self.config_dict[h]]*len(self.conditions_dict)
            except:
                if h not in ['Concentration_1','Condition_1','Concentration_2',
                             'Condition_2','Index','gDNA well','StartOD','EndOD']:
                    print('Missing header: {}'.format(h))
                pass 
        self.df['Concentration_1'] = condition1_list 
        self.df['Concentration_2'] = condition2_list
        self.df['Growth Plate wells'] = well_list 
        self.df['gDNA well'] = well_list 
        
    def write_OD_data(self,OD_file):
        # optional function if you have OD data file to write in
        # OD_file should be a tsv table read from a microplate reader
        # It will assume it is in a 8x12 microplate format and will align to 'Growth Plate wells' in DEM.df
        pass
    
    def compile(self,outdir):
        outfile_path = '{}/exps'.format(outdir)
        self.df.to_csv(outfile_path,index=False,sep='\t')
        self.df.to_csv('{}.tsv'.format(outfile_path),index=False,sep='\t')

class DefineExperimentalSetup():
    def __init__(self,setup_df,config_dict):
        self.setup_df = setup_df # experimental set up data
        self.config_dict = config_dict
        self.parse_units()
        
    def parse_units(self):
        # conc_df = slice of DEM.df[['Concentration_1','Units_1','Concentration_2','Units_2']]
        conversion_dict_c = {'M': 1,
                             'mM':0.001,
                             'uM':1e-06} # converts to M
        conversion_dict_v = {'mL':1,
                             'uL':.001} # converts to mL
        # 1) convert working conc. to M
        # 2) convert stock conc. M
        exp_volume = float(self.config_dict['Working Volume'])
        exp_volume_units = self.config_dict['Working Volume Units']
        exp_volume_mL = conversion_dict_v[exp_volume_units]*exp_volume

        cond1_stock_conc = float(self.config_dict['Condition_1 Stock Conc.'])
        cond1_stock_unit = self.config_dict['Condition_1 Stock Units']
        cond1_stock_conc_M = cond1_stock_conc*conversion_dict_c[cond1_stock_unit]
        cond2_stock_conc = float(self.config_dict['Condition_2 Stock Conc.'])
        cond2_stock_unit = self.config_dict['Condition_2 Stock Units']
        cond2_stock_conc_M = cond2_stock_conc*conversion_dict_c[cond2_stock_unit]

        well_name = []
        volumes_cond1 = []
        volumes_cond2 = [] 
        volume_media = []
        for row in self.setup_df.to_numpy():
            exp_cond1_conc_M = row[2]*conversion_dict_c[row[3]] 
            exp_cond2_conc_M = row[5]*conversion_dict_c[row[6]]

            volume_needed_cond1 = exp_cond1_conc_M*exp_volume_mL/cond1_stock_conc_M
            volume_needed_cond2 = exp_cond2_conc_M*exp_volume_mL/cond2_stock_conc_M

            well_name.append(row[0])
            volumes_cond1.append(volume_needed_cond1)
            volumes_cond2.append(volume_needed_cond2)
            volume_media.append(exp_volume_mL-volume_needed_cond1-volume_needed_cond2)
        
        volumes_setup_df = pd.DataFrame()
        volumes_setup_df['well_name'] = well_name 
        volumes_setup_df['{} {} stock volume needed mL'.format(self.config_dict['Condition_1'], self.config_dict['Condition_1 Stock Conc.']+self.config_dict['Condition_1 Stock Units'])] = volumes_cond1 
        volumes_setup_df['{} {} stock volume needed mL'.format(self.config_dict['Condition_2'], self.config_dict['Condition_2 Stock Conc.']+self.config_dict['Condition_2 Stock Units'])] = volumes_cond2 
        volumes_setup_df['media_vol_needed mL'] = volume_media
        self.volumes_setup_df = volumes_setup_df
  
def check_essential_args(args):
    for arg in vars(args):
        f = getattr(args,arg)
        if f is not None and os.path.exists(f)==False:
            print('{} does not exist'.format(f))
            sys.exit()
        else:
            pass
    
def main(args):
    config_file = args.config
    outdir = args.outdir
    if args.od_start is not None:
        od_start = args.od_start
    if args.od_end is not None:
        od_end = args.od_end
    R = ReadConfig(config_file)
    config_dict = R.config_dict
    DEM = DefineExperimentMeta(config_dict)
    DEM.compile(outdir)
    setup_df = DEM.df[['gDNA well','Condition_1','Concentration_1','Units_1',
                       'Condition_2','Concentration_2','Units_2']]
    DES = DefineExperimentalSetup(setup_df, config_dict)
    DES.volumes_setup_df.to_csv('{}/experiment_setup.tsv'.format(outdir), sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Creates exps file for FEBA rbseq and experimental setup file')
    parser.add_argument('-c','--config',help='exps_setup.config',required=True)
    parser.add_argument('-o','--outdir',help='output directory',required=True)
    parser.add_argument('-d1','--od_start',help='OD data inoculation',required=False)
    parser.add_argument('-d2','--od_end',help='OD data endpoint',required=False)

    args = parser.parse_args()
    check_essential_args(args)
    main(args)


