#!bin/bash
##########################################
# Written by Dr. Bryan Venters, EpiCypher Inc. 
# Updated 29 OCT 2021
#
# Purpose of script: Use "grep -c" to count exact match to CUTANA spike-in nucleosome barcodes from unzipped paired-end (R1 & R2) fastqs. Both R1 and R2 are searched because the barcoded side of nucleosome may be ligated in either direction relative to the P5 or P7 adapters.
# 
# Caution: avoid spaces in file paths and filenames because it is a frequent source of syntax errors.
##########################################

## Instructions ##
# In Finder (Mac):
# 1. Duplicate this template shell script to create an experiment-specific copy. Save to a desired folder. 
# 2. Unzip (extract) all fastq.gz to fastq on Mac or Linux by double-clicking on *.gz files. Save these files to the same folder containing shell script.

# In TextEdit:
# 3. Copy/Paste the lines between "# template loop begin ##" and "# template loop end ##" (below; from the first "echo" to the last "done") in this shell script as many times as needed (1 template loop per paired-end R1 & R2 data set). This is the script needed to align samples to the spike-in barcodes.
# 4. Add sample filenames (eg: H3K4me3_dNuc_R1_100k.fastq) in place of "sample1_R1.fastq" in copied/pasted loops. Each loop should contain R1 and R2 files for the matched reaction.

# In Terminal:
# 5. In Terminal, cd (change directory) to directory (folder) containing unzipped_fastq files (eg: cd path_to_fastq). This can also be done in a Mac by dragging your file onto Terminal.
# 6. To execute shell script in the Terminal application type "sh" followed by the file path to your saved shell script, or by dragging shell script file into Terminal: sh path_to_shell.sh
#7: Press enter. The shell script will output the barcode counts in the order listed under "# Barcode identities" (below) and datasets will be annotated based on the filenames. For each annotated loop (R1 and R2 sample set), it will generate all R1 counts and then all R2 counts.

# In Excel:
# Note: The Excel template provides space to copy in read count data for the IgG negative control, the H3K4me3 positive control, and 6 additional samples. Copy and make additional analyses as needed.
# 8. Copy and paste R1 and R2 barcode counts for control reactions from Terminal to yellow highlighted cells.  
# 9. For reactions using a K-methyl antibody represented in the panel, select the on-target name from the dropdown list in column B. Then copy the R1 and R2 barcode read counts generated from running the script and paste into the appropriate highlighted yellow cells. Make sure the R1 and R2 files are matched (i.e. from the same reaction).
# 10. The Excel file will automatically generate a heatmap visualization of results, where read counts are normalized relative to the on-target PTM specified in column B. Condensed results, summarizing antibody binding data for each PTM in panel, is provided to right.
# 11. Enter the total "Unique align reads" from your sequencing reaction at the bottom of each table. The total barcode reads will auto-fill. The target for "% total barcode reads" is 1-10%. Lower than 1% and there may be insufficient reads to determine specificity. Greater than 10% indicates spike-ins should be diluted further in future experiments.
# 12. Go to the "Output Table" sheet in the Excel Workbook for a combined heat map showing antibody specificity data for all 8 samples.

# Notes:
# EpiCypher considers an antibody with <20% binding to all off-target PTMs specific and suitable for downstream data analysis. 
# For IgG, data is normalized to the sum of total barcode reads.



# Barcode identities
# Unmodified (A & B)
# TTCGCGCGTAACGACGTACCGT 
# CGCGATACGACCGCGTTACGCG 

# H3K4me1 (A & B)
# CGACGTTAACGCGTTTCGTACG 
# CGCGACTATCGCGCGTAACGCG 

# H3K4me2 (A & B)
# CCGTACGTCGTGTCGAACGACG 
# CGATACGCGTTGGTACGCGTAA 

# H3K4me3 (A & B)
# TAGTTCGCGACACCGTTCGTCG 
# TCGACGCGTAAACGGTACGTCG

# H3K9me1 (A & B)
# TTATCGCGTCGCGACGGACGTA
# CGATCGTACGATAGCGTACCGA

# H3K9me2 (A & B)
# CGCATATCGCGTCGTACGACCG
# ACGTTCGACCGCGGTCGTACGA

# H3K9me3 (A & B)
# ACGATTCGACGATCGTCGACGA
# CGATAGTCGCGTCGCACGATCG

# H3K27me1 (A & B)
# CGCCGATTACGTGTCGCGCGTA
# ATCGTACCGCGCGTATCGGTCG

# H3K27me2 (A & B)
# CGTTCGAACGTTCGTCGACGAT
# TCGCGATTACGATGTCGCGCGA

# H3K27me3 (A & B)
# ACGCGAATCGTCGACGCGTATA
# CGCGATATCACTCGACGCGATA

# H3K36me1 (A & B)
# CGCGAAATTCGTATACGCGTCG
# CGCGATCGGTATCGGTACGCGC

# H3K36me2 (A & B)
# GTGATATCGCGTTAACGTCGCG
# TATCGCGCGAAACGACCGTTCG

# H3K36me3 (A & B)
# CCGCGCGTAATGCGCGACGTTA
# CCGCGATACGACTCGTTCGTCG

# H4K20me1 (A & B)
# GTCGCGAACTATCGTCGATTCG
# CCGCGCGTATAGTCCGAGCGTA

# H4K20me2 (A & B)
# CGATACGCCGATCGATCGTCGG
# CCGCGCGATAAGACGCGTAACG

# H4K20me3 (A & B)
# CGATTCGACGGTCGCGACCGTA
# TTTCGACGCGTCGATTCGGCGA


# template loop begin ##

#updated shellscript to be able to work with .gz files which are defaulted in cut&run pipeline 


for i in *.gz; do
	echo "$i finding match against 250bp 22nt barcodes = 8 counts" #>> MetStat_Reads.csv

	for barcode in TTCGCGCGTAACGACGTACCGT CGCGATACGACCGCGTTACGCG CGACGTTAACGCGTTTCGTACG CGCGACTATCGCGCGTAACGCG CCGTACGTCGTGTCGAACGACG CGATACGCGTTGGTACGCGTAA TAGTTCGCGACACCGTTCGTCG TCGACGCGTAAACGGTACGTCG TTATCGCGTCGCGACGGACGTA CGATCGTACGATAGCGTACCGA CGCATATCGCGTCGTACGACCG ACGTTCGACCGCGGTCGTACGA ACGATTCGACGATCGTCGACGA CGATAGTCGCGTCGCACGATCG CGCCGATTACGTGTCGCGCGTA ATCGTACCGCGCGTATCGGTCG CGTTCGAACGTTCGTCGACGAT TCGCGATTACGATGTCGCGCGA ACGCGAATCGTCGACGCGTATA CGCGATATCACTCGACGCGATA CGCGAAATTCGTATACGCGTCG CGCGATCGGTATCGGTACGCGC GTGATATCGCGTTAACGTCGCG TATCGCGCGAAACGACCGTTCG CCGCGCGTAATGCGCGACGTTA CCGCGATACGACTCGTTCGTCG GTCGCGAACTATCGTCGATTCG CCGCGCGTATAGTCCGAGCGTA CGATACGCCGATCGATCGTCGG CCGCGCGATAAGACGCGTAACG CGATTCGACGGTCGCGACCGTA TTTCGACGCGTCGATTCGGCGA ;
	do
		zgrep -c $barcode "$i" >> MetStat_Reads.csv
done

done



# template loop end ##




