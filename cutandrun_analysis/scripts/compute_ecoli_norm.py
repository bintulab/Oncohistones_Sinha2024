# Script: cutrun_pipeline_part2_20211019.py
# Author: Connor Ludwig
# Organization: Bintu Lab, Stanford University
# Package versions: samtools 1.7 ; bamCoverage 3.3.1 
# Date: 10/19/2021

import sys
import os
import matplotlib
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
from datetime import datetime

# NOTE: run this script after running the first part of the CUT&RUN pipeline
# NOTE: the scripts folder should be parallel with the other folders

# specify analysis summary CSVs for 1) all samples aligned to human genome, and 2) all samples aligned to E. coli genome
analysis_summary = pd.read_csv(sys.argv[1])
ecoli_stats = pd.read_csv(sys.argv[2])

# supply Admera sample sheet for converting Admera names to sample names
samplekey = pd.read_csv(sys.argv[3])
file_root = list(samplekey['AdmeraID'])
sample = list(samplekey['SampleID'])
sampleDict = dict(zip(file_root, sample))

ecoli_rename = {i:('Ecoli '+i) for i in ecoli_stats.columns if i != 'Sample'}
ecoli_stats = ecoli_stats.rename(columns=ecoli_rename)
human_rename = {i:('Human '+i) for i in analysis_summary.columns if i != 'Sample'}
analysis_summary = analysis_summary.rename(columns=human_rename)
analysis_summary = pd.merge(analysis_summary, ecoli_stats, on='Sample', how='left')

analysis_summary['Ecoli Uniquely Aligned'] = analysis_summary['Ecoli Aligned'] - analysis_summary['Ecoli Duplicates']
analysis_summary['Human Uniquely Aligned'] = analysis_summary['Human Aligned'] - analysis_summary['Human Duplicates'] - analysis_summary['Human Removed']

analysis_summary['Spike-in Percent'] = round(100*analysis_summary['Ecoli Uniquely Aligned']/analysis_summary['Human Uniquely Aligned'], 4)
analysis_summary['Scale Factor'] = round(1/analysis_summary['Spike-in Percent'], 4)
analysis_summary['Human Normalized Reads'] = round(analysis_summary['Human Uniquely Aligned']*analysis_summary['Scale Factor'], 0)
analysis_summary['AdmeraID'] = analysis_summary['Sample'].str.split('_').str[0]

analysis_summary['Sample'] = analysis_summary['AdmeraID'].map(sampleDict)
analysis_summary['Antibody'] = analysis_summary['Sample'].str.split('_').str[4]
analysis_summary = analysis_summary.sort_values(by=['Antibody', 'AdmeraID'])
print(analysis_summary)

save_path = os.path.join(os.path.dirname(os.path.dirname(sys.argv[2])), 'ecoli-normalized_read_stats.csv')
analysis_summary.to_csv(save_path, index=False)

sns.barplot(data=analysis_summary, x='Sample', y='Human Normalized Reads')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(save_path[:-3] + 'png')

# New line
