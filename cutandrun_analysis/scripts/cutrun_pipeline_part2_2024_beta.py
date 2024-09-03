# Script: cutrun_pipeline_part2_20211019.py
# Author: Connor Ludwig
# Organization: Bintu Lab, Stanford University
# Package versions: samtools 1.7 ; bamCoverage 3.3.1 
# Date: 10/19/2021

import sys
import os
import numpy as np
import pandas as pd
from datetime import datetime

# NOTE: run this script after running the first part of the CUT&RUN pipeline
# NOTE: the scripts folder should be parallel with the other folders

# usage: if no argument is specified, this will use CPM
# 		 if an argument is specified, this will assume that is a CSV file with E coli stats

# find the master directory with the structure described above
outputDir = sys.argv[1]

# user specifies normalization method:
#    'cpm' = normalize to total counts (sequencing depth normalization)
#    'ecoli' = normalize to E. coli spike-in DNA for all samples
#    'kmet' = for samples with non-zero values in the KmetStat column, normalize to these counts; otherwise, normalize to E. coli spike-in DNA

norm_method = sys.argv[2] # mandatory
genome_analysis_summary = sys.argv[3] # needed for 'ecoli' and 'kmet' normalization methods
ecoli_analysis_summary = sys.argv[4] # needed for 'ecoli' and 'kmet' normalization methods
kmet_analysis_summary = sys.argv[5] # only needed for 'kmet' normalization method

picardDir_path = os.path.join(outputDir, 'picard_outputs')
dedupscaled_path = os.path.join(outputDir, 'dedup_scaled')

if (norm_method == 'ecoli') | (norm_method == 'kmet'):
	ecoli_stats = pd.read_csv(ecoli_analysis_summary)
	ecoli_stats['Ecoli'] = ecoli_stats['Aligned'] - ecoli_stats['Duplicates']
	genome_stats = pd.read_csv(genome_analysis_summary)
	if list(ecoli_stats['Sample'])[0][-3:] == '_R2':
		ecoli_stats['Sample'] = ecoli_stats['Sample'].str.split('_').str[:-1].str.join('_')
		genome_stats['Sample'] = genome_stats['Sample'].str.split('_').str[:-1].str.join('_')

	if 'Removed' in genome_stats.columns:
		genome_stats['Uniquely Aligned'] = genome_stats['Aligned'] - genome_stats['Duplicates'] - genome_stats['Removed']
	else:
		genome_stats['Uniquely Aligned'] = genome_stats['Aligned'] - genome_stats['Duplicates']
	genome_stats = pd.merge(genome_stats[['Sample', 'Uniquely Aligned']], ecoli_stats[['Sample', 'Ecoli']], on='Sample', how='left')
	genome_stats['Spike-in Percent'] = 100*genome_stats['Ecoli']/genome_stats['Uniquely Aligned']
	genome_stats['Ecoli Scale Factor'] = 1/genome_stats['Spike-in Percent']
	scale_factors = dict(zip(list(genome_stats['Sample']), list(genome_stats['Ecoli Scale Factor'])))

if norm_method == 'kmet':
	kmet_stats = pd.read_csv(kmet_analysis_summary)
	kmet_stats['Ecoli Scale Factor'] = list(genome_stats['Ecoli Scale Factor'])
	kmet_stats['Uniquely Aligned'] = list(genome_stats['Uniquely Aligned'])
	kmet_stats['Spike-in Percent'] = 100*kmet_stats['KmetStat']/kmet_stats['Uniquely Aligned']
	
	#calculate minimal subset
	grouped = kmet_stats.groupby('Antibody')
	
	for antiname, group in grouped:
		min_val = group['Spike-in Percent'].min()
		for index, row in group.iterrows():
			# print(f"antiname={antiname}, min_val={min_val}, index={index}, row={row}, result={row['Spike-in Percent'] / min_val}")
			kmet_stats.at[index, 'Scale_Factor'] = row['Spike-in Percent'] / min_val
			# print("scale factor computed for one item", kmet_stats['Scale_Factor'])
	
	#kmet_stats['Kmet Scale Factor'] = 1/kmet_stats['Spike-in Percent']
	#updated testing for manual scale factor 
	#kmet_stats['Kmet Scale Factor'] = list(kmet_stats['Scale_Factor']) #newly added testing this:
	
	kmet_stats['Scale Factor'] = np.where(kmet_stats['Spike-in Percent']!=0, kmet_stats['Scale_Factor'], kmet_stats['Ecoli Scale Factor'])
	scale_factors = dict(zip(list(kmet_stats['Sample']), list(kmet_stats['Scale Factor'])))
	print("scale_factors", scale_factors)

# print(scale_factors)
# save outputs to a records file
records = open(os.path.join(outputDir, 'normalization_records.txt'), 'a')

# make dedup_scaled folder if one doesn't exist
if os.path.isdir(dedupscaled_path) == False:
	print('Making directory called dedup_scaled to store deduplicated and normalized bedgraph files')
	print('Making directory called dedup_scaled to store deduplicated and normalized bedgraph files', file=records)
	os.mkdir(dedupscaled_path)

# for bam file deduplicated by Picard, index with with samtools and make CPM-normalized bedgraphs with bamCoverage
for file in os.listdir(picardDir_path):
	if file[-7:] == 'pic.bam':
		fileID_components = file.split('_')[:-1]
		fileID = '_'.join(fileID_components)
		file_path = os.path.join(picardDir_path, file)
		dedupfile = fileID + '_dedup.scaled.bedgraph'
		dedupfile_path = os.path.join(dedupscaled_path, dedupfile)

		index_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
		print(' ~ Indexing deduplicated BAM file %s @ %s ~' % (file, index_time))
		print(' ~ Indexing deduplicated BAM file %s @ %s ~' % (file, index_time), file=records)
		index_cmd = 'samtools index ' + file_path
		print(index_cmd)
		print(index_cmd, file=records)
		os.system(index_cmd)
		
		norm_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
		if norm_method == 'cpm':
			print(' ~ Normalizing deduplicated BAM file using CPM %s @ %s ~' % (file, norm_time))
			print(' ~ Normalizing deduplicated BAM file using CPM %s @ %s ~' % (file, norm_time), file=records)
			norm_cmd = 'bamCoverage --bam ' + file_path + ' -o ' + dedupfile_path + ' --outFileFormat bedgraph --extendReads --centerReads -p 12 --binSize 10 --normalizeUsing CPM'
		elif norm_method == 'ecoli':
			print(' ~ Normalizing deduplicated BAM file based on E. coli spike-in %s @ %s ~' % (file, norm_time))
			print(' ~ Normalizing deduplicated BAM file based on E. coli spike-in %s @ %s ~' % (file, norm_time), file=records)
			file_root = '_'.join(fileID_components[:-2]) # CHANGE HERE DEPENDING ON HOW FILES WERE GENERATED IN FIRST PART
			norm_cmd = 'bamCoverage --bam ' + file_path + ' -o ' + dedupfile_path + ' --outFileFormat bedgraph --extendReads --centerReads -p 12 --binSize 10 --scaleFactor ' + str(scale_factors[file_root])
		elif norm_method == 'kmet':
			print(' ~ Normalizing deduplicated BAM file based on KmetStat spike-in %s @ %s ~' % (file, norm_time))
			print(' ~ Normalizing deduplicated BAM file based on KmetStat spike-in %s @ %s ~' % (file, norm_time), file=records)
			file_root = '_'.join(fileID_components[:-2]) # CHANGE HERE DEPENDING ON HOW FILES WERE GENERATED IN FIRST PART
			norm_cmd = 'bamCoverage --bam ' + file_path + ' -o ' + dedupfile_path + ' --outFileFormat bedgraph --extendReads --centerReads -p 12 --binSize 10 --scaleFactor ' + str(scale_factors[file_root])
		print(norm_cmd)
		print(norm_cmd, file=records)
		os.system(norm_cmd)

finish_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
print(' ~ Batch normalization of CUT&RUN files completed @ %s ~' % finish_time)
print(' ~ Batch normalization of CUT&RUN files completed @ %s ~' % finish_time, file=records)
records.close()

# New line

