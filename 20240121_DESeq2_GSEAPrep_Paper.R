## DIFFERENTIAL EXPRESSION ANALYSIS & PREP FOR GSEA

## SETUP -----------------------------------------------------------------------------------------

library("DESeq2")
library("ggplot2")
#prevent errors due to max overlaps globally for session
options(ggrepel.max.overlaps = Inf)
library("pheatmap")
library("reshape")
library(ggrepel)
library(RColorBrewer)
library(dplyr)
library(readxl)


# SETUP & RUN DESEQ2 ------------------------------------------------------------------------------

#Set seed to allow reproducibility in L2FC
set.seed(13)

#Create DESeq Object from counts matrix
#Input Counts Matrix
cts <- as.matrix(read.csv("Manuscript/count_matrix.csv", row.names="featureID"))
coldata <- read.csv("Manuscript/coldata.csv", row.names = 1)
coldata$Mucus <- factor(coldata$Mucus)
coldata$MucusType <- factor(coldata$MucusType)

#Ensure count matrix and column metadata have same sample order
head(cts, 5)
coldata
#Good to go!

#Round values in cts to the nearest integer value (required for DESeq2 input)
round_cts <- round(cts, digits=0)
head(round_cts) #rounded correctly

dds <- DESeqDataSetFromMatrix(countData = round_cts, colData = coldata, design = ~ MucusType)

#Specify the reference level (which sample counts as the control group). Default is alphabetical order.
#Alternatively, this can be specified using contrast() when building results table.
dds$MucusType <- relevel(dds$MucusType, ref = "None")

#Run the statistical tests
dds <- DESeq(dds, test = "Wald", fitType = "parametric") 

#Results for IEC Mucus vs. CDMM test only
res <- results(dds, name="MucusType_IEC_Mucus_vs_None")
head(res)

#Determine total number of genes meeting LFC and p-value cutoff
sum(res$log2FoldChange >= 1 | res$log2FoldChange <= -1 & res$padj < 0.05, na.rm = TRUE) #612 genes

#Perform LFC Shrinkage, as recommended by Michael Love for bulk RNASeq
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="MucusType_IEC_Mucus_vs_None")
#check columns of results table
head(resLFC)

#Sort results by Log2FC
resLFC_sorted <- resLFC[order(resLFC$log2FoldChange),]
head(resLFC_sorted)
tail(resLFC_sorted)

#How many pass LFC threshold?
sum(resLFC_sorted$log2FoldChange >= 1 | resLFC_sorted$log2FoldChange <= -1 & resLFC_sorted$padj < 0.05, na.rm=TRUE) #567 genes

#Confirm statistical tests used
mcols(res)$description
mcols(resLFC)$description

#Write results to CSV files
write.csv(as.data.frame(res), file="Manuscript/IEC_vs_CDMM_DESeq2_NoShrinkage.csv")
write.csv(as.data.frame(resLFC_sorted), file="Manuscript/IEC_vs_CDMM_DESeq2.csv")

#Subset list of differentially expressed genes and export
resLFC_sig = subset(resLFC_sorted, resLFC_sorted$log2FoldChange >= 1 | resLFC_sorted$log2FoldChange <= -1 & resLFC_sorted$padj < 0.05)
write.csv(as.data.frame(resLFC_sig), file="Manuscript/IEC_vs_CDMM_Sig.csv")


## DATA VISUALIZATION OPTIONS (HEATMAPS, VOLCANO PLOTS) ---------------------------------------------------------

#Make heatmaps and Volcano plots. Includes scripts to narrow plotted results or label by specified padj or L2FC cutoffs.
#Using Data Visualization tutorial from HBC Training Website:

### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 1
lfc.cutoff.dn = -1

#create logical vector of values to determine whether a gene meets significance threshold or not
threshold_c <- ifelse(res$padj < padj.cutoff & abs(res$log2FoldChange) > lfc.cutoff, ifelse(res$log2FoldChange > lfc.cutoff, "U", "D"), "N")
threshold_bin = res$padj < padj.cutoff & abs(res$log2FoldChange) > lfc.cutoff

threshold_up = res$padj < padj.cutoff & res$log2FoldChange > lfc.cutoff
threshold_dn = res$padj < padj.cutoff & res$log2FoldChange < lfc.cutoff.dn

#Create vector and append to results table
res$threshold_c <- threshold_c

#subset just the significant results
sig <- data.frame(subset(res, threshold_bin==TRUE))

#remove rRNA genes (clutter), these are first 5 rows
sig = sig[-c(1:5),]

#Create heatmap using pheatmap
#extract normalized counts (not VST)
normalized_counts <- counts(dds, normalized=T)
norm_sig <- normalized_counts[rownames(sig),]
head(norm_sig, n=20)

### Annotate our heatmap (optional)
annotation <- data.frame(MucusType=coldata[,'MucusType'],
                         row.names=rownames(coldata))

### Set a color palette
heat.colors <- rev(brewer.pal(9, "RdBu"))

### Run pheatmap
#note: scale = "row" option breaks the colors in the scale by Z-score, which helps with visualization. Z-score is determined after clustering, so samples are still clustered by their counts.
heatmap_lfc2 <- pheatmap(mat = norm_sig, color = heat.colors, clustering_method="complete", cluster_rows = T, show_rownames=T, annotation= annotation, border_color="white", fontsize = 10, scale = "row", fontsize_row = 10, height=20)
heatmap_lfc1 <- pheatmap(mat = norm_sig, color = heat.colors, clustering_method="complete", cluster_rows = T, show_rownames=F, annotation= annotation, fontsize = 10, scale = "row", fontsize_row = 10, height=20)
#Save the heatmap as .svg

##Making Volcano plots
#Add the threshold: 
#create logical vector of values to determine whether a gene meets significance threshold or not
threshold_iec <- res$padj < padj.cutoff & abs(res$log2FoldChange) > lfc.cutoff
threshold_iec_c = ifelse(res$padj < padj.cutoff & abs(res$log2FoldChange) > lfc.cutoff, ifelse(res$log2FoldChange > lfc.cutoff, "U", "D"), "N")
threshold_iec_bin=res$padj < padj.cutoff & abs(res$log2FoldChange) > lfc.cutoff

#Determine how many genes pass threshold
length(which(threshold_iec))

#Create vector and append to results table
res$threshold <- threshold_iec_c

#Start with the data frame of our results
res_df <- data.frame(res)
View(res_df)
#Confirm that there is a "threshold" column before plotting (created above)

##Volcano plot labeling largest fold change
res_df_ordered_lfc <- res_df[order(res_df$log2FoldChange), ] 
#Remove last row, which is a gene with no counts
res_df_ordered_lfc <- res_df_ordered_lfc %>% filter(row_number() <= n()-1)

res_df_ordered_lfc$genelabels <- rownames(res_df_ordered_lfc) %in% rownames(res_df_ordered_lfc[c(1:5, 3530:3534),])

res_df_ordered_lfc$labelnames = unlist(lapply(strsplit(rownames(res_df_ordered_lfc),"_"),'[',2))

View(res_df_ordered_lfc)

iec_vs_cdmm_volc_lfc <- ggplot(res_df_ordered_lfc) + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold, alpha = 0.8, size = ifelse(genelabels==T, 5, 3))) +
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = ifelse(genelabels == T, labelnames,""), size=6)) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "gray100"), 
        panel.border = element_rect(colour = "gray50", fill=NA, size=0.5)) +
  scale_color_manual(values=c("#0077BB","#BBBBBB","#CC3311")) +
  scale_x_continuous(limits=c(-9,5), breaks=c(-8,-6,-4,-2,0,2,4)) +
  geom_hline(yintercept=1.3, linetype="dashed") +
  geom_vline(xintercept=c(-1, 1), linetype="dashed")
iec_vs_cdmm_volc_lfc
#Can export as .svg direcly from plot viewer in RStudio, or use ggsave below

##GGSave
ggsave(file="20240115_Volcano_Plot.svg", plot=iec_vs_cdmm_volc_lfc, width=7, height=4.5, dpi=300)


## FILE PREPARATION FOR GSEA -------------------------------------------------------------------------------------

#Create tables compatible with the Desktop version of GSEA. More information can be found here:
#https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Using_RNA-seq_Datasets_with_GSEA

#First, make .gct file with dataset to be analyzed.

#Making a gct file for input to GSEA
#Need normalized counts for each sample, to start
norm_counts <- counts(dds, normalized = T)
head(norm_counts)
tail(norm_counts)
norm_counts <- as.data.frame(norm_counts)
head(norm_counts)
#Create description and name columns, as required for .gct format. Description can be as simple as the gene names
norm_counts$description <- rownames(norm_counts)
#norm_counts$NAME <- norm_counts$description
tail(norm_counts)
#Move the description column to be the second in the dataset.
require(dplyr)
norm_counts <- norm_counts %>% relocate(description, .before = CDMM_A)
#norm_counts <- norm_counts %>% relocate(NAME, .before = description)
head(norm_counts)
#Save file as a .txt
write.table(norm_counts, "20231227_full_norm_dataset.txt", sep = "\t", quote = F, row.names = F)
#Upon export, remove "gene-" prefix from the beginning of each locus tag and save as .txt, which can then be imported to GSEA.

#Then, obtain information from KEGG servers to match C. difficile pathways/gene sets with locus tags.

#Obtain pathway information from KEGG for gene set curation, compatible with GSEA
#End goal is to create a .gmt or .gmx file compatible with GSEA
#information obtained from: https://www.researchgate.net/post/How_i_can_get_a_list_of_KEGG_pathways_and_its_list_of_genes
#Note: I had to update R to 4.1.2 before this would work; BiocManager needs to be compatible with R

#step 1: install related packages and open them in R:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("KEGGREST")
BiocManager::install("EnrichmentBrowser")
library("KEGGREST")
library("EnrichmentBrowser")
#step2: check and obtain a list of entry identifiers (in this case: cdl for C. difficile R20291) and associated definition for a given database or a given set of database entries.
R20291 <- keggList("cdl")
#step 3: download the pathways of that organism:
cdlpathway <- downloadPathways("cdl")
#step 4: retrieve gene sets for an organism from databases such as GO and KEGG:
cdl <- getGenesets(org = "cdl", db = "kegg", cache = TRUE, return.type="list")
#step5: Parse and write the gene sets to a flat text file in GMT format for other pathway enrichment analysis programs (e.g., GSEA):
writeGMT(cdl, gmt.file = "20231227_kegg_cdl.gmt")

#Used exported .txt and .gmt files to run GSEA on desktop


## PREPARING RANK FILE FOR GSEAPRERANKED WITH 567 MOST DE GENES ONLY ------------------------------------------------------------------

#To run GSEA with a focused, pre-ranked list of signficant genes only, created a .rnk file with the subsetted, significant genes.
#Making a pre-ranked .rnk file for preranked GSEA analysis.
#Using the IEC vs. CDMM dataset, ranking based on shrunken LFC values >1 or < -1.
#Note that the values in the .rnk file need not actually be in order.

#Create list of all genes passing p-adj < 0.05 threshold
resLFC_p = subset(resLFC_sorted, resLFC_sorted$padj < 0.05)
sig_df = as.data.frame(resLFC_p)

#Create a gene IDs column
sig_df$geneIDs <- rownames(sig_df)
head(sig_df)

#Create a table with gene IDs as the first column, Log2FC as the second column.
#shuffle column orders
require(dplyr)
sig_df <- sig_df %>% relocate(geneIDs, .before = log2FoldChange)
head(sig_df)
#Select out just the geneIDs and log2FoldChange columns
sig_rnk <- sig_df %>% select(geneIDs, log2FoldChange)
head(sig_rnk)
#export table as tab delimited txt file, then open in text editor to remove column names and add # to top row. Then same as .rnk.
write.table(sig_rnk, "20231227_iec_sig_ranked.txt", sep = "\t", quote = F, row.names = F)

#Narrow further to only include genes meeting p-adj and LFC cutoff
sig2_df = as.data.frame(resLFC_sig)
#Create a gene IDs column
sig2_df$geneIDs <- rownames(sig2_df)
head(sig2_df)
#Create a table with gene IDs as the first column, Log2FC as the second column.
#shuffle column orders
require(dplyr)
sig2_df <- sig2_df %>% relocate(geneIDs, .before = log2FoldChange)
head(sig2_df)
#Select out just the geneIDs and log2FoldChange columns
sig2_rnk <- sig2_df %>% select(geneIDs, log2FoldChange)
head(sig2_rnk)
#export table as tab delimited txt file, then open in text editor to remove column names and add # to top row. Then same as .rnk.
write.table(sig2_rnk, "20231227_iec_sig2_ranked.txt", sep = "\t", quote = F, row.names = F)


## HEATMAPS FROM FOCUSED GSEA ANALYSIS ------------------------------------------------------------------------------

#Link gene set IDs to genes in normalized counts
#ENRICHED GENE SET FROM SELF-RANKED GENE LIST (BASED ON LFC) WITH DE GENES ONLY
#Input genes from cdl00051 (only DE genes)
cdl00051_DE = read.delim("~/2020_RNASeq_RT_Muc/CDL00051_DE_Core_2.txt", sep='\t', header=TRUE)
cdl00051_DE = merge(normalized_counts, cdl00051_DE, by.x=0, by.y="SYMBOL", all=FALSE)
cdl00051_DE_names = cdl00051_DE[,-1]
rownames(cdl00051_DE_names)=cdl00051_DE[,1]
#Make a new heatmap
heatmap_cdl00051_DE <- pheatmap(mat = cdl00051_DE_names, color = heat.colors, clustering_method="complete", cluster_rows = T, show_rownames=T, annotation= annotation, border_color="white", fontsize = 10, scale = "row", fontsize_row = 10, height=20)

#Save
ggsave(file="20240115_heatmap.svg", plot=heatmap_cdl00051_DE, width=7, height=5, dpi=300)


## CREATE FILE WITH EXPRESSION DATA FOR EVERY GENE AND ANNOTATIONS FROM GENEIOUS --------------------------------------

#Merge list of significant genes from IEC vs. CDMM comparison with their annotations from Geneious
geneious = read_excel("~/2020_RNASeq_RT_Muc/R20291_FN545816_Annotations.xlsx", col_names=TRUE)
geneious = subset(geneious, geneious$Type=="CDS")
#Remove "gene-" preceding each locus tag.
resLFC_sorted = as.data.frame(resLFC_sorted)
rownames(resLFC_sorted) = gsub("gene-", "", rownames(resLFC_sorted))
#Merge the two datasets together
annotated = merge(resLFC_sorted, geneious, by.x=0, by.y="locus_tag", all=FALSE)
annotated_names = annotated[,-1]
rownames(annotated_names)=annotated[,1]
#Add TPM values for individual IEC & CDMM samples
tpm = read.delim("~/2020_RNASeq_RT_Muc/tpm_matrix.txt", sep='\t', header=TRUE)
annotated_tpm = merge(annotated_names, tpm, by.x=0, by.y="featureID", all=FALSE)
annotated_tpm_names = annotated_tpm[,-1]
rownames(annotated_tpm_names) = annotated_tpm[,1]
#Export for further fine tuning, sorting, etc.
write.csv(annotated_tpm_names, "./IEC_vs_CDMM_Annotations_TPM_All.csv")


## END OF ANALYSIS
