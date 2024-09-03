if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install('GenomicFeatures')
BiocManager::install('ChIPseeker')
BiocManager::install('clusterProfiler')
# BiocManager::install('annotables')

BiocManager::install("stephenturner/annotables")

# purpose of script: annotate peaks with nearest gene
# code based on snippets from https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionV/lessons/12_annotation_functional_analysis.html
library(GenomicFeatures)
library(ChIPseeker) # this does the peak annotation! see annotatePeak function below
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
library(annotables)
library(org.Hs.eg.db)

### CHANGE ONLY THE THREE PATHS BELOW FOR YOUR EXPERIMENT ###
infile = file.path('F:/Joydeb/20230124_CNR_6_0/PEAKS_K36M_DOX_HA-tag_intersected.bed') # change path here; this should be final output of process_peaks.py
annofile = file.path('F:/Joydeb/20230124_CNR_6_0/PEAKS_K36M_DOX_HA-tag_annotated.txt') # change path here; this is the annotated peaks file produced by this script
GOcsv = file.path('F:/Joydeb/20230124_CNR_6_0/PEAKS_K36M_DOX_HA-tag_GOterms.csv') # change path here; this is the GO terms file produced by this script
### NO NEED TO TOUCH ANYTHING BELOW THIS LINE ###

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene # transcript database for pulling annotations
peakAnnoList <- lapply(infile, annotatePeak, TxDb=txdb, # how you annotate peaks! Simple function
                       tssRegion=c(-1000, 1000), verbose=FALSE)
peakAnnoList # this just prints preview of annotated peaks

# plots
plotAnnoPie(peakAnnoList[[1]]) # pie chart showing what fraction of peaks are in various regulatory or coding regions
plotDistToTSS(peakAnnoList, title="Distribution of incorporated oncohistone relative to TSS") # x-axis is wrong, should be % of oncohistone peaks!

# writing to file and add gene symbols
onco_annot <- as.data.frame(peakAnnoList[[1]]@anno)
entrezids <- unique(onco_annot$geneId)
entrez2gene <- grch37 %>% 
  filter(entrez %in% entrezids) %>% 
  dplyr::select(entrez, symbol)

m <- match(onco_annot$geneId, entrez2gene$entrez)
onco_annot <- cbind(onco_annot[,1:14], geneSymbol=entrez2gene$symbol[m], onco_annot[,15:16]) # add in gene symbol column to main dataframe

write.table(onco_annot, annofile, sep="\t", quote=F, row.names=F) # saves dataframe to file name specified above
onco_annot

# gene enrichment - e.g. do oncohistones incorporate preferentially into genes associated with certain biological processes?
entrezids <- onco_annot$geneId %>% # get your entrez IDs
  as.character() %>% 
  unique()

ego <- enrichGO(gene = entrezids,     # GO enrichment function
                keyType = "ENTREZID", 
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
write.csv(cluster_summary, GOcsv)

dotplot(ego, showCategory=10) #, cex.axis=2)

# KEGG pathway enrichment analysis <- not working for some reason even though species is correct
# failing due to package error
# ekegg <- enrichKEGG(gene = entrezids,
#                     organism = 'hsa',
#                     keyType = 'kegg',
#                     pvalueCutoff = 0.05)
# dotplot(ekegg)

