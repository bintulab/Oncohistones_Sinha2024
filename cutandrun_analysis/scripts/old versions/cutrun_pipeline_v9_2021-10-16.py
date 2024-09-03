# Script: cutrun_pipeline_v9_2021-10-16.py
# Author: Connor Ludwig
# Organization: Bintu Lab, Stanford University
# Date: 10/16/2021
# Package versions:Bowtie 2 2.3.4.1 ; samtools 1.7 ;  Picard 11.0.7
# NOTE: before starting, need to organize paired read files into subfolders within a parent folder that will be specified as the first argument

import os
import sys
import pandas as pd
from datetime import datetime

# Check if the correct number of arguments were passed to the script
argv_len = len(sys.argv)
if (argv_len < 4) | (argv_len > 5):
	print('ERROR: INCORRECT NUMBER OF ARGUMENTS DETECTED')
	print('Please specify the following in order:')
	print('   1) Path to directory containing sub-directories with pairs of zipped fastq files')
	print('   2) Path to reference genome (common root - exclude .bt2 extensions)')
	print('   3) Path to the Picard tool for quantifying and marking duplicates')
	print('   4) Optional: Path to script for editing sam files')
	sys.exit()

# Assign inputs to variable names
inputDir_fastqzip = sys.argv[1]
inputRoot_refgenome = sys.argv[2]
picardTool = sys.argv[3]

inputDir_parent = os.path.dirname(os.path.dirname(inputDir_fastqzip))
records = open(os.path.join(inputDir_parent, 'analysis_records.txt'), 'a')

# Initialize
init_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
print(' ~ Initializing for bowtie2 alignment @ %s ~' % init_time)
print(' ~ Initializing for bowtie2 alignment @ %s ~' % init_time, file=records)

# Make directories and get paths; store commands in variables
bowtie2Dir = 'bowtie2-aligned'
bowtie2Dir_path = os.path.join(inputDir_parent, bowtie2Dir)
os.mkdir(bowtie2Dir_path)
bowtie2_cmd = 'bowtie2 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x '

alnstatsDir = 'alignment_stats'
alnstatsDir_path = os.path.join(inputDir_parent, alnstatsDir)
os.mkdir(alnstatsDir_path)

##########################################################################################

# Process each zipped fastq file in input directory
for folder in os.listdir(inputDir_fastqzip):
	# Define temporary path for subdirectory with alignment files
	tempDir_path = os.path.join(inputDir_fastqzip, folder)
	files = []
	for file in os.listdir(tempDir_path):
		files.append(os.path.join(tempDir_path, file))

	fileID_components = file.split('_')[:-1]
	fileID = '_'.join(fileID_components)
	
	# Create file to store alignment stats
	alnstats = fileID + '_alnstats.txt'
	alnstats_path = os.path.join(alnstatsDir_path, alnstats)

	# Use bowtie2 to align reads
	aligned_name = fileID + '_aligned.sam'
	aligned_path = os.path.join(bowtie2Dir_path, aligned_name)
	bowtie2_cmd_full = bowtie2_cmd + inputRoot_refgenome + ' -p 16 -1 ' + files[0] + ' -2 ' + files[1] + ' -S ' + aligned_path + ' 2> ' + alnstats_path
	
	bowtie2_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
	print(' ~ Performing alignment with bowtie2 for sample %s @ %s ~' % (fileID, bowtie2_time))
	print(' ~ Performing alignment with bowtie2 for sample %s @ %s ~' % (fileID, bowtie2_time), file=records)
	print(bowtie2_cmd_full)
	print(bowtie2_cmd_full, file=records)
	os.system(bowtie2_cmd_full)

##########################################################################################

# If a script for editing sam files is specified, call this
if argv_len == 5:
	sam_edit_script = sys.argv[4]
	sam_edit_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
	sam_edit_cmd = 'python ' + sam_edit_script + ' ' + bowtie2Dir_path + ' ' + os.path.join(inputDir_parent, 'analysis_records.txt')
	print(' ~ Preparing to edit sam files with %s' % sam_edit_script)
	print(' ~ Preparing to edit sam files with %s' % sam_edit_script, file=records)
	records.close()
	os.system(sam_edit_cmd)
	records = open(os.path.join(inputDir_parent, 'analysis_records.txt'), 'a')

# If not, proceed with normal conversion to bam file and indexing
elif argv_len == 4:
	sortedbamDir = 'sortedbam'
	sortedbamDir_path = os.path.join(inputDir_parent, sortedbamDir)
	os.mkdir(sortedbamDir_path)

	index_cmd = 'samtools index '

	for file in os.listdir(bowtie2Dir_path):
		file_path = os.path.join(bowtie2Dir_path, file)
		# Convert sam to bam file
		sortedbam_name = file[:-3] + 'sorted.bam'
		sortedbam_path = os.path.join(sortedbamDir_path, sortedbam_name)
		sortedbam_cmd = 'samtools sort ' + file_path + ' -o ' + sortedbam_path

		# Sort sam file and convert to bam file
		sortedbam_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
		print(' ~ Sorting and converting SAM file to BAM file @ %s ~' % sortedbam_time)
		print(' ~ Sorting and converting SAM file to BAM file @ %s ~' % sortedbam_time, file=records)
		print(sortedbam_cmd)
		print(sortedbam_cmd, file=records)
		os.system(sortedbam_cmd)
	
		# Index bam file
		index_cmd_full = index_cmd + sortedbam_path
		index_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
		print(' ~ Indexing BAM file @ %s ~' % index_time)
		print(' ~ Indexing BAM file @ %s ~' % index_time, file=records)
		print(index_cmd_full)
		print(index_cmd_full, file=records)
		os.system(index_cmd_full)

##########################################################################################

# Quantify and mark duplicates with Picard
if argv_len == 5:
	sortedbamDir_path = os.path.join(inputDir_parent, 'sortedbam_edited')

# Initialize
init_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
print(' ~ Initializing for marking duplicates with Picard @ %s ~' % init_time)
print(' ~ Initializing for marking duplicates with Picard @ %s ~' % init_time, file=records)

# Make directories and get paths; store commands in variables
picardDir = 'picard_outputs'
picardDir_path = os.path.join(inputDir_parent, picardDir)
os.mkdir(picardDir_path)

# Iterate over sorted.bam files
for file in os.listdir(sortedbamDir_path):
	if file[:-3] == 'bai':
		continue
	fileID_components = file.split('_')[:-1]
	fileID = '_'.join(fileID_components)
	file_path = os.path.join(sortedbamDir_path, file)
	picbam_name = fileID + '_pic.bam'
	picbam_path = os.path.join(picardDir, picbam_name)
	metric_name = fileID + '_metrics.txt'
	metric_path = os.path.join(picardDir, metric_name)
	picard_cmd = 'java -jar %s MarkDuplicates -I %s -O %s -M %s --REMOVE_DUPLICATES true' % (picardTool, file_path, picbam_path, metric_path)

	picard_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
	print(' ~ Marking duplicates with Picard for sample %s @ %s ~' % (fileID, picard_time))
	print(' ~ Marking duplicates with Picard for sample %s @ %s ~' % (fileID, picard_time), file=records)
	print(picard_cmd)
	print(picard_cmd, file=records)
	os.system(picard_cmd)

# Compile quantification of duplicates into a single text file
dupfile = 'marked_duplicates.txt'
dupfile_path = os.path.join(inputDir_parent, dupfile)
with open(dupfile_path, 'w') as outFile:
	for file in os.listdir(picardDir_path):
		if file[-11:] == 'metrics.txt':
			fileID_components = file.split('_')[:-1]
			fileID = '_'.join(fileID_components)
			file_path = os.path.join(picardDir_path, file)
		
			with open(file_path, 'r') as inFile:
				parse = False
				for line in inFile:
					if parse:
						readDups = line.rstrip('\n').split('\t')[6]
						outLine = fileID + '\t' + readDups + '\n'
						outFile.write(outLine)
						parse = False
					if len(line) < 7:
						continue
					elif line[:7] == 'LIBRARY':
						parse = True
					else:
						continue

##########################################################################################

# Extract read information from the alignment stats files
sampleList1 = []
unaligned = []
single = []
multi = []

for file in os.listdir(alnstatsDir_path):
	fileID_components = file.split('_')[:-1]
	fileID = '_'.join(fileID_components)
	file_path = os.path.join(alnstatsDir_path, file)
	sampleList1.append(fileID)

	with open(file_path, 'r') as inFile:
		for line in inFile:
			if '0 times' in line:
				unaligned.append(int(line.lstrip().split(' (')[0]))
			elif 'exactly 1 time' in line:
				single.append(int(line.lstrip().split(' (')[0]))
			elif '>1 times' in line:
				multi.append(int(line.lstrip().split(' (')[0]))
			else:
				continue

df = pd.DataFrame({'Sample':sampleList1,
				   'Unaligned':unaligned,
				   'Single':single,
				   'Multi':multi})
df['Aligned'] = df['Single'] + df['Multi']

##########################################################################################

# Subtract out the detected duplicate reads
sampleList2 = []
duplicates = []

with open(dupfile_path, 'r') as inFile:
	for line in inFile:
		line_elements = line.rstrip('\n').split('\t')
		if argv_len == 5:
			sample = '_'.join(line_elements[0].split('_')[:-1])
		else:
			sample = line_elements[0]
		
		sampleList2.append(sample)
		duplicates.append(int(line_elements[1]))
df2 = pd.DataFrame({'Sample':sampleList2,
					'Duplicates':duplicates})

df = pd.merge(df, df2, on='Sample', how='outer')

##########################################################################################

# Remove reads from sam editing if applicable
if argv_len == 5:
	# initialize lists to store file contents
	sampleList3 = []
	removedList = []
	sam_edit_output = os.path.join(inputDir_parent, 'sam_editing.txt')
	# iterate through file and extract information
	with open(sam_edit_output, 'r') as inFile:
		for line in inFile:
			if line[:6] == 'Number':
				sample = line.rstrip('\n').split(' ')[5].split('_aligned.sam:')[0]
				sampleList3.append(sample)
				removed = int(line.rstrip('\n').split(' ')[-1])
				removedList.append(removed)
			else:
				continue
	inFile.close()

	df3 = pd.DataFrame({'Sample':sampleList3,
						'Removed':removedList})

	df = pd.merge(df, df3, on='Sample', how='outer')

	# Expand sample name to file name
	df['File'] = df['Sample'] + '_edited_aligned.sorted.bam'

##########################################################################################

# Compute reads per million scaling factor and generate text file for scaling
summary_path = os.path.join(inputDir_parent, 'analysis_summary.csv')
df.to_csv(summary_path, index=False)

# Tell user that processing is complete
remove_cmd = 'rm marked_duplicates.txt sam_editing.txt'
os.system(remove_cmd)

finish_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
print(' ~ Processing of CUT&RUN files completed @ %s ~' % finish_time)
print(' ~ Processing of CUT&RUN files completed @ %s ~' % finish_time, file=records)

records.close()

# New line
