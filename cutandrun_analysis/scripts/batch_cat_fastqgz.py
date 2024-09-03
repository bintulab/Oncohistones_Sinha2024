# Script: batch_cat_fastqgz.py
# Author: Connor Ludwig
# Organization: Bintu Lab, Stanford University
# Date: 11/9/2021

import os
import sys

# input should be folder containing split fastq.gz files
inputDir = sys.argv[1]
parentDir = os.path.dirname(os.path.dirname(inputDir))

# make new directory to store combined fastq.gz files
cat_fastqDir = 'cat_fastq'
if os.path.isdir(cat_fastqDir) == False:
	print('Making directory called cat_fastq to store concatenated fastq.gz files')
	os.mkdir(cat_fastqDir)

# extract roots (part of file name that is shared)
rootList = []
for file in os.listdir(inputDir):
	file_root = file.split('_L00')[0]
	if file_root not in rootList:
		rootList.append(file_root)

# for each root, find corresponding files, define path for new output file, and execute cat command
for root in rootList:
	catList = []
	for file in os.listdir(inputDir):
		if root in file:
			catList.append(os.path.join(inputDir, file))
	catList = ' '.join(catList)
	cat_path = os.path.join(cat_fastqDir, root + '.fastq.gz')
	cat_cmd = 'cat ' + catList + ' > ' + cat_path
	print(cat_cmd)
	os.system(cat_cmd)

# new line
