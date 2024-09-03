import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('pdf') # do this because environment does not have GUI backend
import matplotlib.pyplot as plt
import seaborn as sns

# ACTB - chr7
# KCNQ1 - chr11
#coordinates of elements to be quanitified 
# EDIT THE CONTENTS OF THIS DICTIONARY HERE
elemDict = {'AAVS1':['chr19', 55600286, 55631005],
			'ACTB':['chr7', 5564779, 5572232],
			'KCNQ1':['chr11', 2464238, 2872335]}

#save input directory and initialize dictionary of elements
inputDir = sys.argv[1]
initDict = {key:0 for key in elemDict}

#initialize dictionary of sample conditions (antibody, dox)
antibodies = ['H3K9me3', 'H3K36me2', 'H3K36me3', 'H3K27me3', 'H3K4me3', 'HA-tag']
celltypes = ['K36M', 'PARENT']
conditions = ['nodox', 'dox']

sigDict_keynames = []
for a in antibodies:
	for c in celltypes:
		for d in conditions:
			sigDict_keynames.append('_'.join([a, c, d]))

sigDict = {k:initDict.copy() for k in sigDict_keynames}

AbDict = {'H3K9me3':'H3K9me3',
		  'H3K36me2':'H3K36me2',
		  'H3K36me3':'H3K36me3',
		  'H3K27me3':'H3K27me3',
		  'H3K4me3':'H3K4me3',
		  'HA-tag':'HA-tag',}

doxDict = {'ND':'nodox',
		   'DOX':'dox'}

#read in bedgraph files and sum signal at different elements in various sample conditions 
for file in os.listdir(inputDir):
	if file[-8:] != 'bedgraph':
		continue
	filepath = os.path.join(inputDir, file)
	with open(filepath, 'r') as bedgraph:
		print(file)
		file_contents = file.split('_')
		celltype = file_contents[1]
		condition = doxDict[file_contents[2]]
		antibody = AbDict[file_contents[3].split('.')[0]]
		sigDictkey = '_'.join([antibody, celltype, condition])
		for line in bedgraph:
			chrom = line.rstrip('\n').split('\t')[0]
			start = int(line.rstrip('\n').split('\t')[1])
			end = int(line.rstrip('\n').split('\t')[2])
			signal = float(line.rstrip('\n').split('\t')[3])
			for key in elemDict:
				if (chrom == elemDict[key][0]) & (start >= elemDict[key][1]) & ((end - 1) <= (elemDict[key][2] + 1)):
					sigDict[sigDictkey][key] += signal
	bedgraph.close()

print(sigDict)
#format signal into lists to make dataframe
AbList = []
cellList = []
conditionList= []
elemList = []
sigList = []
for sample in sigDict:
	Abs, celltype, condition = sample.split('_')
	for element in sigDict[sample]:
		AbList.append(Abs)
		cellList.append(celltype)
		conditionList.append(condition)
		elemList.append(element)
		sigList.append(sigDict[sample][element])
#construct dataframe
df = pd.DataFrame({'Antibody':AbList, 'Cell Type':cellList, 'Condition':conditionList,
				   'Element':elemList, 'Signal':sigList})

#define order of elements and conditions
elemOrder = ['AAVS1', 'ACTB', 'KCNQ1']
condOrder = ['nodox', 'dox']
#reorder and sort dataframe
df.Element = pd.Categorical(df.Element, categories=elemOrder, ordered=True)
df.Condition = pd.Categorical(df.Condition, categories=condOrder, ordered=True)
df = df.sort_values(by=['Condition', 'Cell Type', 'Antibody', 'Element'])
print(df)
#calculate length of element and signal/kb
elemStart = {}
elemEnd = {}
for key in elemDict:
	elemStart[key] = elemDict[key][1]
	elemEnd[key] = elemDict[key][2]
df['Start'] = df['Element'].map(elemStart)
df['End'] = df['Element'].map(elemEnd)
df['Length'] = df['End'].astype(int) - df['Start'].astype(int)
df['Signal per kb'] = 1000*df['Signal']/df['Length']

#save dataframe for plotting signal 
df.to_csv('bedgraphAUC.csv', index=False)
print(df)

#plot signal for all conditions 
pal = ['#999999', '#333333']
# g = sns.FacetGrid(data=df, col='Antibody', hue='Condition', aspect=3, height=5)
g = sns.catplot(data=df, x='Element', y='Signal per kb', hue='Condition',
				col='Antibody', row='Cell Type', kind='bar', height=4, aspect=1, palette=pal)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.tight_layout()
plt.savefig('bedgraphAUC_plot_col-antibody_sigperkb.png')
plt.close()