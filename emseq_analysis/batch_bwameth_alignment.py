# Script Name: batch_bwameth_alignment.py
# Author: Connor Ludwig, with modifications by Abby Thurm
# Organization: Bintu Lab, Stanford University
# Date: 08/30/2022


# NOTE: need to organize paired read files into a subfolder in a parent 'fastazip' folder

# Pipeline to process and analyze EM-seq data, specifically mapping reads to amplicon
# Start with fastq files of paired-end reads
# Usage: python {path to batch_bwameth_alignment.py} {path to directory with zipped fasta files} {path to reference amplicon}

import os
import sys
from datetime import datetime

# Assign inputs to variable names
inputDir = sys.argv[1]
inputRoot_amplicon = sys.argv[2]

inputDir_parent = os.path.dirname(os.path.dirname(inputDir))

# Initialize
init_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
print(' ~ Initializing for bwameth alignment @ %s ~' % init_time)

# Make directories and get paths; store commands in variables
samDir = 'sam'
samDir_path = os.path.join(inputDir_parent, samDir)
os.mkdir(samDir_path)

# updated 10/3/2021 to output unaligned files in case transgene transcripts are mapped there
# also allows identification of contamination
unalDir = 'unaligned'
unalDir_path = os.path.join(inputDir_parent, unalDir)
os.mkdir(unalDir_path)

bamDir = 'bam'
bamDir_path = os.path.join(inputDir_parent, bamDir)
os.mkdir(bamDir_path)

for folder in os.listdir(inputDir):
	# Define temporary path for subdirectory with alignment files
	tempDir_path = os.path.join(inputDir, folder)
	files = []
	for file in os.listdir(tempDir_path):
		files.append(os.path.join(tempDir_path, file))

	root = '_'.join(file.split('_')[:-4])
	sam_name = root + '.sam'
	sam_path = os.path.join(samDir_path, sam_name)
	unal_name = root + '_unaligned.fastq'
	unal_path = os.path.join(unalDir_path, unal_name)

	print(inputRoot_amplicon, files[0], files[1], sam_path)
	bwameth_cmd = 'python ~/win_f/Abby/EM-Seq/bwa-meth-master/bwameth.py --reference %s %s %s > %s' % (inputRoot_amplicon, files[0], files[1], sam_path)
	bwameth_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
	print(' ~ Aligning file @ %s ~' % bwameth_time)
	print(bwameth_cmd)
	os.system(bwameth_cmd)

	# Convert sam to bam file
	bam_name = root + '.sorted.bam'
	bam_path = os.path.join(bamDir_path, bam_name)
	# bam_cmd = 'samtools view -S -b ' + sam_path + ' > ' + bam_path
	bam_cmd = 'samtools sort ' + sam_path + ' -o ' + bam_path # updated 10/3/2021 to reduce later sort steps and extra files
	bam_index_cmd = 'samtools index ' + bam_path

	bam_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
	print(' ~ Converting SAM file to BAM file @ %s ~' % bam_time)
	print(bam_cmd)
	print(bam_index_cmd)
	os.system(bam_cmd)
	os.system(bam_index_cmd)

# Tell user that processing is complete
finish_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
print(' ~ Processing of EM-seq files completed @ %s ~' % finish_time)
