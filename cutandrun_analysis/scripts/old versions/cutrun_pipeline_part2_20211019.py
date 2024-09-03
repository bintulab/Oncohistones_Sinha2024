# Script: cutrun_pipeline_part2_20211019.py
# Author: Connor Ludwig
# Organization: Bintu Lab, Stanford University
# Package versions: samtools 1.7 ; bamCoverage 3.3.1 
# Date: 10/19/2021

import sys
import os
from datetime import datetime

# NOTE: run this script after running the first part of the CUT&RUN pipeline
# NOTE: the scripts folder should be parallel with the other folders

# find the master directory with the structure described above
inputDir = os.path.dirname(os.path.dirname(sys.argv[0]))
picardDir_path = os.path.join(inputDir, 'picard_outputs')
dedupscaled_path = os.path.join(inputDir, 'dedup_scaled')

# save outputs to a records file
records = open(os.path.join(inputDir, 'normalization_records.txt'), 'a')

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
		print(' ~ Normalizing deduplicated BAM file %s @ %s ~' % (file, norm_time))
		print(' ~ Normalizing deduplicated BAM file %s @ %s ~' % (file, norm_time), file=records)
		norm_cmd = 'bamCoverage --bam ' + file_path + ' -o ' + dedupfile_path + ' --outFileFormat bedgraph --extendReads --centerReads -p max --binSize 10 --normalizeUsing CPM'
		print(norm_cmd)
		print(norm_cmd, file=records)
		os.system(norm_cmd)

finish_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
print(' ~ Batch normalization of CUT&RUN files completed @ %s ~' % finish_time)
print(' ~ Batch normalization of CUT&RUN files completed @ %s ~' % finish_time, file=records)
records.close()

# New line
