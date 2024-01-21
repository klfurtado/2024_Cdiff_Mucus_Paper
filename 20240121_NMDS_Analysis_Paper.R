## NMDS, TESTING, & ANALYSIS

## Load Libraries---------------------------------------------
#vegan may have some dependency issues with DESeq2; if running at same time as differential expression analysis, may need to restart R.
library(vegan)
library(randomForest)
library(ggplot2)

## Complete NMDS and statistical testing for complete media conditions +/- mucus---------------------------------
#Set seed for reproducibility
set.seed(13)

# Load flux sampling files
cdmm_fs <- '~/RIPTIDE/riptide_cdmm_max_complete/flux_samples.tsv'
iec_fs <- '~/RIPTIDE/riptide_iec_max_complete/flux_samples.tsv'

# Read in data
cdmm_fs <- read.delim(cdmm_fs, sep='\t', header=TRUE) #323 rxns
cdmm_fs$X <- NULL
iec_fs <- read.delim(iec_fs, sep='\t', header=TRUE) #320 rxns
iec_fs$X <- NULL

# Subsample data (subsample to 250 random samplings without replacement)
sample_size <- min(c(nrow(cdmm_fs), nrow(iec_fs), 250))
sub_sample <- sample(1:min(c(nrow(cdmm_fs), nrow(iec_fs))), sample_size, replace=FALSE)
cdmm_fs <- cdmm_fs[sub_sample,]
iec_fs <- iec_fs[sub_sample,]
rm(sample_size, sub_sample)

# Limit to overlapping reactions
overlap <- intersect(colnames(cdmm_fs), colnames(iec_fs))
cdmm_fs <- cdmm_fs[, overlap]
iec_fs <- iec_fs[, overlap]
rm(overlap)

# Format row names
cdmm_names <- paste('cdmm_', 1:nrow(cdmm_fs), sep='')
rownames(cdmm_fs) <- cdmm_names
iec_names <- paste('iec_', 1:nrow(iec_fs), sep='')
rownames(iec_fs) <- iec_names

# Create metadata
cdmm_metadata <- cbind(cdmm_names, rep('NoMucus', length(cdmm_names)))
iec_metadata <- cbind(iec_names, rep('Mucus', length(iec_names)))
metadata <- rbind(cdmm_metadata, iec_metadata)
colnames(metadata) <- c('label', 'mucus')
metadata <- as.data.frame(metadata)
rm(cdmm_metadata, iec_metadata)

# Merge data and prep for unsupervised learning
all_samples <- rbind(cdmm_fs, iec_fs)
flux_groups <- as.factor(c(rep('CDMM', nrow(cdmm_fs)), rep('IEC', nrow(iec_fs))))
all_samples_adj <- all_samples + abs(min(all_samples)) #adds 1000 to flux value for every sample (makes all values positive)
#Run randomForest on true flux values for all samples
rf_obj <- randomForest(flux_groups ~ ., data=all_samples, importance=TRUE, err.rate=TRUE, ntree=1500, mtry=15)
rf_obj <- importance(rf_obj, type=1, scale=TRUE) #specifies importance by mean decrease accuracy
rf_mda <- as.data.frame(subset(rf_obj, rf_obj > (abs(min(rf_obj))))) #only include reactions with MDA > 0
all_samples_adj <- all_samples_adj[,rownames(rf_mda)] #pares down all_samples_adj to reactions with MDA > 0

# Calculate dissimilarity (Bray-Curtis)
# may need to restart R before vegan will load; had some depency issues thanks to DESeq2
flux_dist <- vegdist(all_samples_adj, method='bray') # Bray-Curtis

# Mean within-group dissimilarity
meandist(flux_dist, grouping=flux_groups)

# Unsupervised learning (NMDS)
flux_nmds <- as.data.frame(metaMDS(flux_dist, k=2, trymax=50)$points)

# Center points
flux_x <- (abs(max(flux_nmds$MDS1)) - abs(min(flux_nmds$MDS1))) / 2
flux_y <- (abs(max(flux_nmds$MDS2)) - abs(min(flux_nmds$MDS2))) / 2
flux_nmds$MDS1 <- flux_nmds$MDS1 - flux_x
flux_nmds$MDS2 <- flux_nmds$MDS2 - flux_y
flux_xlim <- max(abs(max(flux_nmds$MDS1)), abs(min(flux_nmds$MDS1))) + 0.01
flux_ylim <- max(abs(max(flux_nmds$MDS2)), abs(min(flux_nmds$MDS2))) + 0.01
write.table(flux_nmds, file='~/RIPTIDE/CDMMvsIEC_nmds_complete_max.tsv', 
            quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE) #Creates output file of MDS1 and MDS2 values

# Statistical testing (permANOVA)
test <- merge(x=metadata, y=all_samples_adj, by.x='label', by.y='row.names')
rownames(test) <- test$label
test$label <- NULL
media_pval <- adonis(flux_dist ~ mucus, data=test, perm=999, method='bray')
media_pval <- media_pval$aov.tab[[6]][1] #p-val = 0.001, significant!
#double check p-value with adonis2 (adonis deprecated)
adonis2(flux_dist ~ mucus, data=test, perm=999, method='bray') #p-value = 0.001


## Use NMDS data to make plots in ggplot2---------------------------------------------------

#Add mucus column to flux_nmds for plotting
flux_nmds_metadata <- merge(flux_nmds, metadata, by.x=0, by.y="label")
rownames(flux_nmds_metadata) <- flux_nmds_metadata$Row.names 
flux_nmds_metadata <- flux_nmds_metadata[,-1]

#Generate the plot
NMDSplot <- ggplot(flux_nmds_metadata,aes(x=MDS1,y=MDS2)) + 
  geom_point(size=3, shape=21, color="black") +
  aes(fill=mucus) + #label by presence of mucus
  scale_fill_manual(values=c("NoMucus" = "#004488", "Mucus" = "#5AAE61")) + #specify color of dots
  theme_bw() +
  labs(fill=NULL) + #remove unnecessary label on top of legend
  theme(legend.position = c(.18,.88), axis.text.x = element_text(size = 14), axis.text.y = element_text(size=14), axis.title = element_text(size=16), legend.text=element_text(size=16)) +
  annotate(geom="text", x=0.03, y=-0.03, label="p = 0.001***", size=5)
plot(NMDSplot)

#Save plot
ggsave(file="NMDS_Complete.svg", plot=NMDSplot, width=6, height=4, dpi=300)
