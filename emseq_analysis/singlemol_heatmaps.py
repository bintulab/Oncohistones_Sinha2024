#second part of EM-Seq pipeline: make single-molecule methylation plots once you have (indexed) .bam files
#usage: pass in folder that holds your bam/sam folders, then your amplicon
#if you haven't made a bed file for your amplicon yet, pass in 'bed' as third argument

import sys
import os
from datetime import datetime

# find the master directory with the structure described above
inputDir = os.path.dirname(sys.argv[1])
heatmapDir_path = os.path.join(inputDir, 'heatmaps')
amplicon = sys.argv[2]
amplicon_root = amplicon.split('.')[0]
bam_path = os.path.join(inputDir, 'bam')
os.mkdir(heatmapDir_path)
records = open(os.path.join(heatmapDir_path, 'heatmap_records.txt'), 'a')

if len(sys.argv) >3:
	#make bed file to pass in to clustering script if not made
	bed_cmd = 'bioawk -c fastx \'{print $name \"\t0\t\"length($seq)}\''+amplicon + '>' + amplicon_root + '.bed'
	print('Making BED file from amplicon FASTA')
	print(bed_cmd)
	os.system(bed_cmd)
	bed_file = amplicon_root + '.bed'
else:
	bed_file = amplicon_root + '.bed'

os.chdir(heatmapDir_path)
amplicon_path = os.path.join(inputDir, amplicon)
bed_path = os.path.join(inputDir, bed_file)
for file in os.listdir(bam_path):
	if file.endswith('.bam'):
		file_path = os.path.join(bam_path, file)
		file_root = file.split('.')[0]

		heatmap_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')[:-7]
		print(' ~ Making single-molecule heatmap for sample %s @ %s ~' % (file, heatmap_time))
		print(' ~ Making single-molecule heatmap for sample %s @ %s ~' % (file, heatmap_time), file=records)
		heatmap_cmd = 'python2.7 ~/win_f/Abby/EM-Seq/dSMF-footprints_optional_clustering.py %s %s CG %s 0 1 2 3 %s %s -unstranded -subset 1000 -cluster -heatmap ~/win_f/Abby/EM-Seq/heatmap.py 10 3 binary 10,100 -minCov 0.7' %(file_path, amplicon_path, bed_path, file_root, amplicon_root)
		print(heatmap_cmd)
		os.system(heatmap_cmd)

#MethylDackel
os.chdir(inputDir)
for file in os.listdir(bam_path):
	if file.endswith('.bam'):
		file_path = os.path.join(bam_path, file)

		print(' ~Creating MethylDackel BedGraphs for sample %s' % file_path)
		methyldackel_cmd = 'MethylDackel extract --CHG --CHH %s %s' % (amplicon, file_path)
		print(methyldackel_cmd)
		os.system(methyldackel_cmd)