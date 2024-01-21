#RANDOM FORESTS ANALYSIS

##Setup---------------------------------------
#Load libraries for session
library(randomForest)
library(ggplot2)
library(ggbreak)
library(reshape2)
library(svglite)
#Set working directory
setwd("~/RIPTIDE")


## Run Random Forests Analysis (generates MDA values) for complete media conditions +/- mucus------------------------
#Set seed for reproducibility
set.seed(13)

# Flux sampling files
# data generated from riptide.maxfit_contextualize()
cdmm_fluxes <- '~/RIPTIDE/riptide_cdmm_max_complete/flux_samples.tsv'
iec_fluxes <- '~/RIPTIDE/riptide_iec_max_complete/flux_samples.tsv'

# Read in data
cdmm_fluxes <- read.delim(cdmm_fluxes, sep='\t', header=TRUE)
cdmm_fluxes$X <- NULL
iec_fluxes <- read.delim(iec_fluxes, sep='\t', header=TRUE)
iec_fluxes$X <- NULL

# Subsample
sample_size <- min(c(nrow(cdmm_fluxes), nrow(iec_fluxes), 250))
sub_sample <- sample(1:min(c(nrow(cdmm_fluxes), nrow(iec_fluxes))), sample_size, replace=FALSE)
cdmm_fluxes <- cdmm_fluxes[sub_sample,]
iec_fluxes <- iec_fluxes[sub_sample,]
rm(sample_size, sub_sample)

# Format data for random forest
shared_rxns <- intersect(colnames(cdmm_fluxes), colnames(iec_fluxes))
cdmm_fluxes <- cdmm_fluxes[,shared_rxns]
iec_fluxes <- iec_fluxes[,shared_rxns]
#remove unnecessary columns
cdmm_fluxes[,c('biomass','dna_rxn','rna_rxn','peptidoglycan_rxn','protein_rxn','lipid_rxn','cofactor_rxn','teichoicacid_rxn','cellwall_rxn')] <- NULL
iec_fluxes[,c('biomass','dna_rxn','rna_rxn','peptidoglycan_rxn','protein_rxn','lipid_rxn','cofactor_rxn','teichoicacid_rxn','cellwall_rxn')] <- NULL
rm(shared_rxns)

# Add column and merge
cdmm_fluxes$media <- rep('cdmm', nrow(cdmm_fluxes))
iec_fluxes$media <- rep('iec', nrow(iec_fluxes))
media_fluxes <- as.data.frame(rbind(cdmm_fluxes, iec_fluxes))
rownames(media_fluxes) <- c(paste0('cdmm_', c(1:nrow(cdmm_fluxes))), paste0('iec_', c(1:nrow(iec_fluxes))))

# Find most informative metabolites with RF
condition <- as.factor(media_fluxes$media)
media_fluxes$media <- NULL
media_fluxes <- droplevels(media_fluxes)
rf_obj <- randomForest(condition ~ ., data=media_fluxes, importance=TRUE, err.rate=TRUE, ntree=1500, mtry=15)
rf_obj <- importance(rf_obj, type=1, scale=TRUE)
rm(condition)

# Subset to most informative features
rf_mda <- as.data.frame(subset(rf_obj, rf_obj > (abs(min(rf_obj)))))
rf_mda$reaction <- rownames(rf_mda)
rf_mda <- rf_mda[order(-rf_mda$MeanDecreaseAccuracy),]
rf_mda <- rf_mda[c(1:20),]
rm(rf_obj)

#To this list of most informative features, add columns for average fluxes in each condition:
rxns <- as.vector(rf_mda$reaction)
cdmm_rxns <- c()
iec_rxns <- c()
for (r in rxns) {
  cdmm_rxns[r]<- median(cdmm_fluxes[,r])
  iec_rxns[r]<- median(iec_fluxes[,r])
}
mda_fluxes <- data.frame(rxn = rxns, cdmm = cdmm_rxns, iec = iec_rxns)
mda_fluxes <- merge(rf_mda, mda_fluxes, by.x="reaction", by.y="rxn")
mda_fluxes <- melt(mda_fluxes, id=c("reaction", "MeanDecreaseAccuracy"), value.name = "flux")

#Generate file
write.table(mda_fluxes, file='~/RIPTIDE/mda_complete.tsv', 
            quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

#Take generated .tsv file and do the following before proceeding:
#Add a 'name' column. Copy names for corresponding reaction codes from model into column.
#Save spreadsheet with '_name' at end of filename prior to re-importing.

## Plots for complete media--------------------------------------------------
#Reimport final file for Plotting: 
MDA_flux <- read.delim('~/RIPTIDE/mda_complete_name.txt', sep='\t', header=TRUE)

#Reorder based on increasing MDA! (Annoying, but this resolves ggplot2's automatic factor ordering which makes no sense)
MDA_flux$name <- factor(MDA_flux$name, levels=c("N-Acetyl-D-glucosamine exchange",
                                                "L-Glutamate exchange",
                                                "L-Leucine hydrogen symport",
                                                "H+ exchange",
                                                "trehalose transport via PEP:Pyr PTS",
                                                "2,3,4,5-Tetrahydrodipicolinate:NAD+ oxidoreductase",
                                                "Amino acid transporter (ala-L)",
                                                "Na+/glutamate symport",
                                                "L-Leucine exchange",
                                                "UDP-N-acetyl-D-glucosamine pyrophosphohydrolase (periplasm)",
                                                "ATP:L-aspartate 4-phosphotransferase",
                                                "Alanine transaminase",
                                                "L-Alanine exchange",
                                                "L-Aspartate-4-semialdehyde hydro-lyase (adding pyruvate and cyclizing)",
                                                "N-Acetyl-D-glucosamine transport via PEP:Pyr PTS",
                                                "L-glutamate:NAD+ oxidoreductase (deaminating)",
                                                "N-Acetyl-D-glucosamine 1-phosphate 1,6-phosphomutase",
                                                "NADPH:oxidized-thioredoxin oxidoreductase",
                                                "glycine synthase",
                                                "5,10-Methylenetetrahydrofolate:glycine hydroxymethyltransferase"))

MDA_fluxes <- ggplot(MDA_flux, aes(y=name, x=flux, fill=variable)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5)+
  geom_point(size=4, shape=21, color="black")+
  scale_fill_manual(values=c("cdmm" = "#332288", "iec" = "#5AAE61"))+
  theme_bw() +
  labs(fill=NULL) + #remove unnecessary label on top of legend
  theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust=0.5)) +
  scale_x_break(c(-500, -200), scales=4) +scale_x_break(c(200, 500), scales=1)
MDA_fluxes

#Repeat for only top 10 MDA fluxes:
#Reimport final file for Plotting: 
MDA10_flux <- read.delim('~/RIPTIDE/mda10_complete_name.txt', sep='\t', header=TRUE)

#Reorder based on increasing MDA! (Annoying, but this resolves ggplot2's automatic factor ordering which makes no sense)
MDA10_flux$name <- factor(MDA10_flux$name, levels=c("ATP:L-aspartate 4-phosphotransferase",
                                                    "Alanine transaminase",
                                                    "L-Alanine exchange",
                                                    "L-Aspartate-4-semialdehyde hydro-lyase (adding pyruvate and cyclizing)",
                                                    "N-Acetyl-D-glucosamine transport via PEP:Pyr PTS",
                                                    "L-glutamate:NAD+ oxidoreductase (deaminating)",
                                                    "N-Acetyl-D-glucosamine 1-phosphate 1,6-phosphomutase",
                                                    "NADPH:oxidized-thioredoxin oxidoreductase",
                                                    "glycine synthase",
                                                    "5,10-Methylenetetrahydrofolate:glycine hydroxymethyltransferase"))

MDA10_fluxes <- ggplot(MDA10_flux, aes(y=name, x=flux, fill=variable)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5)+
  geom_point(size=4, shape=21, color="black")+
  scale_fill_manual(values=c("cdmm" = "#332288", "iec" = "#5AAE61"))+
  theme_bw() +
  labs(fill=NULL) + #remove unnecessary label on top of legend
  theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust=0.5)) +
  scale_x_break(c(-100, -50), scales=3) +scale_x_break(c(50, 100), scales=1)
MDA10_fluxes

#Save Plots
ggsave(file="MDA_fluxes_complete.svg", plot=MDA_fluxes, width=12, height=5)
ggsave(file="MDA10_fluxes_complete.svg", plot=MDA10_fluxes, width=12, height=4)
