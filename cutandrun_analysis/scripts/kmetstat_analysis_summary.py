# script to count kmetstat barcodes in sequencing data and update analysis summary for normalization

import sys
import os
import csv
import numpy as np
import pandas as pd
 
# csv_data = pd.read_csv('MetStat_Reads.csv')
mod_list = ['Unmodified', 'Unmodified', 'H3K4me1', 'H3K4me1', 'H3K4me2', 'H3K4me2', 'H3K4me3', 'H3K4me3', 'H3K9me1', 'H3K9me1',
			'H3K9me2', 'H3K9me2','H3K9me3','H3K9me3','H3K27me1','H3K27me1','H3K27me2','H3K27me2','H3K27me3','H3K27me3',
			'H3K36me1','H3K36me1','H3K36me2','H3K36me2','H3K36me3','H3K36me3','H4K20me1','H4K20me1','H4K20me2','H4K20me2','H4K20me3','H4K20me3']

barcode_list = ['TTCGCGCGTAACGACGTACCGT', 'CGCGATACGACCGCGTTACGCG', 'CGACGTTAACGCGTTTCGTACG', 'CGCGACTATCGCGCGTAACGCG', 'CCGTACGTCGTGTCGAACGACG', 'CGATACGCGTTGGTACGCGTAA',
				'TAGTTCGCGACACCGTTCGTCG', 'TCGACGCGTAAACGGTACGTCG', 'TTATCGCGTCGCGACGGACGTA', 'CGATCGTACGATAGCGTACCGA', 'CGCATATCGCGTCGTACGACCG', 'ACGTTCGACCGCGGTCGTACGA',
				'ACGATTCGACGATCGTCGACGA', 'CGATAGTCGCGTCGCACGATCG', 'CGCCGATTACGTGTCGCGCGTA', 'ATCGTACCGCGCGTATCGGTCG', 'CGTTCGAACGTTCGTCGACGAT', 'TCGCGATTACGATGTCGCGCGA',
				'ACGCGAATCGTCGACGCGTATA', 'CGCGATATCACTCGACGCGATA', 'CGCGAAATTCGTATACGCGTCG', 'CGCGATCGGTATCGGTACGCGC', 'GTGATATCGCGTTAACGTCGCG', 'TATCGCGCGAAACGACCGTTCG',
				'CCGCGCGTAATGCGCGACGTTA', 'CCGCGATACGACTCGTTCGTCG', 'GTCGCGAACTATCGTCGATTCG', 'CCGCGCGTATAGTCCGAGCGTA', 'CGATACGCCGATCGATCGTCGG', 'CCGCGCGATAAGACGCGTAACG',
				'CGATTCGACGGTCGCGACCGTA', 'TTTCGACGCGTCGATTCGGCGA']

# mod_list = ['H3K36me2','H3K36me2','H3K36me3','H3K36me3']
# barcode_list = ['GTGATATCGCGTTAACGTCGCG', 'TATCGCGCGAAACGACCGTTCG', 'CCGCGCGTAATGCGCGACGTTA', 'CCGCGATACGACTCGTTCGTCG']

inputDir_fastqzip = sys.argv[1] # directory containing subdirectories with zipped fastq files
analysis_summary = sys.argv[2]
analysis_summary = pd.read_csv(analysis_summary)# csv output from part 1 of the analysis pipeline
antibody_field = int(sys.argv[3]) # where antibody information is in the sample string (0-indexed: e.g. fourth field = 3)
inputDir = os.path.dirname(os.path.dirname(inputDir_fastqzip))
kmetstat_path = os.path.join(inputDir, 'kmetstat_counts')

Ab_mat = []
sample_Ab = []
fileID_list = []

qc_matrix = pd.DataFrame({'Mod':mod_list, 'Barcode':barcode_list})

# make kmetstat_counts folder if one doesn't exist
if os.path.isdir(kmetstat_path) == False:
	print('Making directory called kmetstat_counts to store deduplicated and normalized bedgraph files')
	# print('Making directory called dedup_scaled to store deduplicated and normalized bedgraph files', file=records)
	os.mkdir(kmetstat_path)

# iterate through folders with fastq.gz files
for folder in os.listdir(inputDir_fastqzip):
	# Define temporary path for subdirectory with fastq.gz files
	tempDir_path = os.path.join(inputDir_fastqzip, folder)
	Ab_spec = []
	# iterate through each file in the subdirectory
	for file in os.listdir(tempDir_path):
		print('~ Searching for KmetStat barcodes in %s' % file)
		fastqfile_path = os.path.join(tempDir_path, file)
		# get the file components
		fileID_components = file.split('_')[:-1]
		fileID = '_'.join(fileID_components)
		fileID_list.append(fileID)
		kmetfile_path = os.path.join(kmetstat_path, fileID + '.txt')

		# # iterate over the modifications in the mod_list
		# for mod in mod_list:
		# 	# if the file contains the modification in its name, then get the index of the modification from the mod_list above
		# 	if mod in fileID:
		# 		Ab_index = mod_list.index(mod)
		# 	# if not, then report none
		# 	else:
		# 		Ab_index = None

		# search for each barcode in each file and use the index determined above to store the appropriate read counts
		for barcode in barcode_list:
			cmd = 'zgrep -c %s %s >> %s' % (barcode, fastqfile_path, kmetfile_path)
			os.system(cmd)

		# read in saved counts and merge with qc_matrix data frame above
		temp_counts = pd.read_csv(kmetfile_path, header=None).rename(columns={0:fileID})
		temp_counts['Mod'] = mod_list
		temp_counts['Barcode'] = barcode_list
		qc_matrix = pd.merge(qc_matrix, temp_counts, how='inner', on=['Mod', 'Barcode'])

# save qc_matrix
qc_matrix.to_csv('kmetstat_counts_matrix.csv', index=False)

# simplify qc_matrix by summing across barcodes and reads
qc_matrix = qc_matrix.groupby(by='Mod').sum() # sum across barcodes
current_cols = list(qc_matrix.columns) # get current column names (files with R1 or R2)
new_cols = ['_'.join(c.split('_')[:-1]) for c in current_cols] # make new column names (files without R1 and R2 extensions)
qc_matrix = qc_matrix.rename(columns=dict(zip(current_cols, new_cols))).T # rename qc_matrix columns and transpose
qc_matrix = qc_matrix.groupby(by=qc_matrix.index).sum().T # sum across reads and transpose
print(qc_matrix)

# print(analysis_summary)
print(analysis_summary['Sample'].tolist())
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

# end of code
