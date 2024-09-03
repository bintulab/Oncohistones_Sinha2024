import sys
import os
import csv
import numpy as np
import pandas as pd

mod_list = ['Unmodified', 'Unmodified', 'H3K4me1', 'H3K4me1', 'H3K4me2', 'H3K4me2', 'H3K4me3', 'H3K4me3', 'H3K9me1', 'H3K9me1',
			'H3K9me2', 'H3K9me2','H3K9me3','H3K9me3','H3K27me1','H3K27me1','H3K27me2','H3K27me2','H3K27me3','H3K27me3',
			'H3K36me1','H3K36me1','H3K36me2','H3K36me2','H3K36me3','H3K36me3','H4K20me1','H4K20me1','H4K20me2','H4K20me2','H4K20me3','H4K20me3']

qc_matrix = sys.argv[1]
analysis_summary = sys.argv[2]
antibody_field = int(sys.argv[3])

qc_matrix = pd.read_csv(qc_matrix)
analysis_summary = pd.read_csv(analysis_summary)
# simplify qc_matrix by summing across barcodes and reads
qc_matrix = qc_matrix.groupby(by='Mod').sum() # sum across barcodes
current_cols = list(qc_matrix.columns) # get current column names (files with R1 or R2)
new_cols = ['_'.join(c.split('_')[:-1]) for c in current_cols] # make new column names (files without R1 and R2 extensions)
qc_matrix = qc_matrix.rename(columns=dict(zip(current_cols, new_cols))).T # rename qc_matrix columns and transpose
qc_matrix = qc_matrix.groupby(by=qc_matrix.index).sum().T # sum across reads and transpose
print(qc_matrix)


#print(analysis_summary['Sample'].tolist())
analysis_summary['Sample'] = analysis_summary['Sample'].str.split('_').str[:-1].str.join('_') # NOTE: THIS IS TO CORRECT FOR HARD-CODED ERROR THAT PRESERVED '_R2' - CHANGE IF NEEDED
analysis_summary['Antibody'] = analysis_summary['Sample'].str.split('_').str[antibody_field] # pull the antibody info from this field of the name (0-indexed)
analysis_summary['KmetStat'] = 0
# print(analysis_summary)

samples = list(analysis_summary['Sample'])
for s, ab in zip(list(analysis_summary['Sample']), list(analysis_summary['Antibody'])):
	if ab not in mod_list:
		continue
	if s in qc_matrix.columns:
		analysis_summary.loc[analysis_summary['Sample']==s, 'KmetStat'] = sum(list(qc_matrix[s]))
print(analysis_summary)

analysis_summary.to_csv('analysis_summary_kmet-stat.csv', index=False)

