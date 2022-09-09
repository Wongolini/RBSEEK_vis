#%%
import pandas as pd 
import numpy as np 
exps = '/Users/cmdb/Documents/rbseq_ecoli/exps'
gene_counts = '/Users/cmdb/Documents/rbseq_ecoli/gene_counts.tab'
gene_counts_df = pd.read_csv(gene_counts,sep='\t')
exps_df = pd.read_csv(exps,sep='\t')
# %%

exps_dummy = pd.DataFrame(np.empty((96,len(exps_df.columns))),columns=exps_df.columns)
Condition1 = ['Phosphate']*96
Condition2 = ['Nitrate']*96
phosphate = []
nitrate = []
for row in range(8):
    for col in range(12):
        phosphate.append(row/76)
        nitrate.append(col/72)

units1 = ['mM']*96
units2 = ['mM']*96 

# %%
exps_dummy['Condition_1'] = Condition1
exps_dummy['Condition_2'] = Condition2
exps_dummy['Concentration_1'] = phosphate
exps_dummy['Concentration_2'] = nitrate
exps_dummy['SetName'] = ['set1']*96
exps_dummy['Person'] = ['Nick']*96
exps_dummy['Mutant Library'] = ['dummy']*96 
exps_dummy['Units_1'] = units1
exps_dummy['Units_2'] = units2 
exps_dummy['Index'] = np.unique(exps_df['Index'])
exps_dummy['StartOD'] = np.zeros((96))
exps_dummy['EndOD'] = np.random.uniform(size=96,low=.1,high=.5)

# %%
header = ['locusId','sysName','desc','comb']
for setname,index in zip(exps_dummy['SetName'],exps_dummy['Index']):
    header.append(setname+index)
gene_counts_dummy = pd.DataFrame((np.empty((gene_counts_df.shape[0],100))),columns=header)
gene_counts_dummy['locusId'] = gene_counts_df['locusId']
gene_counts_dummy['sysName'] = gene_counts_df['sysName']
gene_counts_dummy['desc'] = gene_counts_df['desc']
gene_counts_dummy['comb'] = gene_counts_df['comb']

# %%
# create a relationship between some genes and phsophates
# choose 10,000 random genes
def stochastic_log_decay_nitrate(c1):
    part1 = (np.random.normal(100000,1000))/(1+np.exp(-1/(c1+.01))) # decay from phosphate
    return part1

def decreasing_relationship(c2):
    p = (np.random.normal(60000,1000))/(1+np.exp(-1/(c2+.01))) 
    if p<0:
        p=np.random.normal(10,1)
    return p 

def increasing_relationship(c1,c2):
    part1 = np.random.normal(10,1)*10*np.exp(c1)
    part2 = np.random.normal(20,2)*c2
    return part1+part2+np.random.normal(10000,1000)


random_gene_i = np.random.choice(range(len(gene_counts_dummy)), 2500, replace=False)
random_gene_j = np.random.choice(range(len(gene_counts_dummy)), 1500, replace=False)
random_gene_k = np.random.choice(range(len(gene_counts_dummy)), 100, replace=False)
random_gene_l = np.random.choice(range(len(gene_counts_dummy)), 100, replace=False)

dummy_mat = np.ones((len(gene_counts_df),len(exps_dummy['Index'])))*np.random.normal(20000,100)

for r in random_gene_i:
    f_list = []
    for index,c1,c2 in zip(exps_dummy['Index'],exps_dummy['Concentration_1'],exps_dummy['Concentration_2']):
        
        f = stochastic_log_decay_nitrate(c1)
        f_list.append(f)
        #print(c1,c2,index,f)
    dummy_mat[r] = f_list


    

for r in random_gene_j:
    f_list = []
    for index,c1,c2 in zip(exps_dummy['Index'],exps_dummy['Concentration_1'],exps_dummy['Concentration_2']):
        f = decreasing_relationship(c2)
        f_list.append(f)
        #print(c2,f,index)
    dummy_mat[r] = f_list
    

for r in random_gene_k:
    f_list = []
    for index,c1,c2 in zip(exps_dummy['Index'],exps_dummy['Concentration_1'],exps_dummy['Concentration_2']):
        f = increasing_relationship(c1,c2)
        f_list.append(f)
    dummy_mat[r] = f_list
 

for r in random_gene_l:
    f_list = np.random.normal()
    for index,c1,c2 in zip(exps_dummy['Index'],exps_dummy['Concentration_1'],exps_dummy['Concentration_2']):
        f_list = np.random.normal(10000,1000,len(exps_dummy))
    dummy_mat[r] = f_list

#%%
gene_counts_dummy.iloc[:,4:] = dummy_mat
gene_counts_dummy.to_csv('gene_counts.tab',sep='\t',index=False)
exps_dummy.to_csv('exps',sep='\t',index=False)

# for the dummy data we should see two deterministic curves:
    # ne exponentially increasing
    # one exponentially decaying
# %%
