import os
import sys
from datetime import datetime

# accessory script to be used with first part of the CUT&RUN pipeline to remove ambiguous reads that fall entirely within duplicated regions

### EDIT THIS ###
# coordinates for pEF mag-sep reporter
readLen = 150
reporter = 'pJT039'
repElmts = {'polyA1':[756, 980], 'pEF':[1396, 2550], 'polyA2':[4474, 4599]}

# specify directory with sam files and records file to continue documentation
samDir = sys.argv[1]
records = open(sys.argv[2], 'a')
parentDir = os.path.dirname(os.path.dirname(samDir))
sam_edit_output = open(os.path.join(parentDir, 'sam_editing.txt'), 'a')

# make directory for edited sam files if one does not exist
editedDir = os.path.join(parentDir, 'sam_edited')
if os.path.isdir(editedDir) == False:
	print('Making directory for edited .sam files (no pEF-, H2B-, or polyA-internal reads)')
	print('Making directory for edited .sam files (no pEF-, H2B-, or polyA-internal reads)', file=records)
	os.mkdir(editedDir)

# make new sam file without reads that fall within pEF, H2B, or the polyA in their entirety
for file in os.listdir(samDir):
	removeList = []
	if file[-11:] == 'aligned.sam':
		samedit_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
		print(' ~ Editing SAM file %s @ %s ~' % (file, samedit_time))
		print(' ~ Editing SAM file %s @ %s ~' % (file, samedit_time), file=records)
		
		# find reads to remove
		new_file = file[:-11] + 'edited_aligned.sam'
		temp_pathIn = os.path.join(samDir, file)
		temp_pathOut = os.path.join(editedDir, new_file)
		with open(temp_pathOut, 'w') as outFile:
			with open(temp_pathIn, 'r') as inFile:
				for line in inFile:
					if line[0] == '@':
						outFile.write(line)
					else:
						line_elements = line.rstrip('\n').split('\t')
						chrom = line_elements[2]
						if chrom not in reporter:
							outFile.write(line)
						else:
							start = int(line_elements[3])
							size = int(line_elements[8])
							if size < 0:
								readLen = -readLen
							end = start + size - readLen
							internal = False
							for elmt in repElmts:
								if (start > repElmts[elmt][0]) & (end < repElmts[elmt][1]):
									internal = True
							if internal == False:
								outFile.write(line)
							else:
								removeList.append(line)

		print('Number of reads removed for %s: %d' % (file, len(removeList)))
		print('Number of reads removed for %s: %d' % (file, len(removeList)), file=sam_edit_output)

sam_edit_output.close()

# make directory for sorted and edited bam files if one does not exist
sortedbamDir = os.path.join(parentDir, 'sortedbam_edited')
if os.path.isdir(sortedbamDir) == False:
	print('Making directory for edited .bam files')
	print('Making directory for edited .bam files', file=records)
	os.mkdir(sortedbamDir)

for file in os.listdir(editedDir):
	edited_path = os.path.join(editedDir, file)
	
	# Sort bam file
	sortedbam_name = file[:-3] + 'sorted.bam'
	sortedbam_path = os.path.join(sortedbamDir, sortedbam_name)
	sortedbam_cmd = 'samtools sort ' + edited_path + ' -o ' + sortedbam_path

	sortedbam_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
	print(' ~ Sorting and converting SAM file to BAM file @ %s ~' % sortedbam_time)
	print(' ~ Sorting and converting SAM file to BAM file @ %s ~' % sortedbam_time, file=records)
	print(sortedbam_cmd)
	print(sortedbam_cmd, file=records)
	os.system(sortedbam_cmd)
	
	# Index bam file
	index_cmd = 'samtools index ' + sortedbam_path
	index_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
	print(' ~ Indexing BAM file @ %s ~' % index_time)
	print(' ~ Indexing BAM file @ %s ~' % index_time, file=records)
	print(index_cmd)
	print(index_cmd, file=records)
	os.system(index_cmd)

records.close()

# New line
