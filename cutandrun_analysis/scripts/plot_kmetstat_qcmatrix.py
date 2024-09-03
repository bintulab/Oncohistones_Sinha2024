import subprocess
import sys
import os
import csv
import pandas as pd
import matplotlib
matplotlib.use('pdf') # do this because environment does not have GUI backend
import matplotlib.pyplot as plt
import seaborn as sns

mod_list = ['Unmodified', 'Unmodified', 'H3K4me1', 'H3K4me1', 'H3K4me2', 'H3K4me2', 'H3K4me3', 'H3K4me3', 'H3K9me1', 'H3K9me1',
			'H3K9me2', 'H3K9me2','H3K9me3','H3K9me3','H3K27me1','H3K27me1','H3K27me2','H3K27me2','H3K27me3','H3K27me3',
			'H3K36me1','H3K36me1','H3K36me2','H3K36me2','H3K36me3','H3K36me3','H4K20me1','H4K20me1','H4K20me2','H4K20me2','H4K20me3','H4K20me3']

barcode_list = ['TTCGCGCGTAACGACGTACCGT', 'CGCGATACGACCGCGTTACGCG', 'CGACGTTAACGCGTTTCGTACG', 'CGCGACTATCGCGCGTAACGCG', 'CCGTACGTCGTGTCGAACGACG', 'CGATACGCGTTGGTACGCGTAA',
				'TAGTTCGCGACACCGTTCGTCG', 'TCGACGCGTAAACGGTACGTCG', 'TTATCGCGTCGCGACGGACGTA', 'CGATCGTACGATAGCGTACCGA', 'CGCATATCGCGTCGTACGACCG', 'ACGTTCGACCGCGGTCGTACGA',
				'ACGATTCGACGATCGTCGACGA', 'CGATAGTCGCGTCGCACGATCG', 'CGCCGATTACGTGTCGCGCGTA', 'ATCGTACCGCGCGTATCGGTCG', 'CGTTCGAACGTTCGTCGACGAT', 'TCGCGATTACGATGTCGCGCGA',
				'ACGCGAATCGTCGACGCGTATA', 'CGCGATATCACTCGACGCGATA', 'CGCGAAATTCGTATACGCGTCG', 'CGCGATCGGTATCGGTACGCGC', 'GTGATATCGCGTTAACGTCGCG', 'TATCGCGCGAAACGACCGTTCG',
				'CCGCGCGTAATGCGCGACGTTA', 'CCGCGATACGACTCGTTCGTCG', 'GTCGCGAACTATCGTCGATTCG', 'CCGCGCGTATAGTCCGAGCGTA', 'CGATACGCCGATCGATCGTCGG', 'CCGCGCGATAAGACGCGTAACG',
				'CGATTCGACGGTCGCGACCGTA', 'TTTCGACGCGTCGATTCGGCGA']


qc_matrix = pd.read_csv(sys.argv[1])
analysis_summary = pd.read_csv(sys.argv[2])
antibody_field = int(sys.argv[3])


# simplify qc_matrix
qc_matrix = qc_matrix.groupby(by='Mod').sum() # sum across barcodes
current_cols = list(qc_matrix.columns) # get current column names (files with R1 or R2)
new_cols = ['_'.join(c.split('_')[:-1]) for c in current_cols] # make new column names (files without R1 and R2 extensions)
qc_matrix = qc_matrix.rename(columns=dict(zip(current_cols, new_cols))).T # rename qc_matrix columns and transpose
qc_matrix = qc_matrix.groupby(by=qc_matrix.index).sum().T # sum across reads and transpose

qc_matrix2 = qc_matrix.T.div(qc_matrix.T.max(axis=1), axis=0)*100 # normalize row by max value in each row (assumed)
qc_matrix2['Order'] = qc_matrix2.index.str.split('_').str[-1].str[1:].astype(int)
qc_matrix2 = qc_matrix2.sort_values(by='Order', ascending=True)
qc_matrix2 = qc_matrix2.drop(columns=['Order'])

# plot matrix
plt.figure(figsize=(8,12))
ax = sns.heatmap(qc_matrix2, cmap='magma', annot=True, fmt='.0f', cbar_kws={'shrink': 0.6})
# ax.hlines(range(len(qc_matrix2)), *ax.get_xlim(), color='#FFFFFF')
ax.vlines([3, 6, 9, 12, 15], *ax.get_ylim(), color='#FFFFFF')
plt.yticks(rotation=0)
plt.xticks(rotation=90)
plt.title('Antibody Specificity (% Target)')
plt.tight_layout()
plt.savefig('KmetStat_qc_matrix.png', dpi=300)

# end of code
