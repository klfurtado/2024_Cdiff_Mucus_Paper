## BIOMASS & METABOLITE FLUXES


## SETUP ------------------------------------------------------------------------------------
#Load necesary libraries for session
#Session run under R 4.3.2
library(reshape2)
library(ggplot2)


## LOAD DATA --------------------------------------------------------------------------------
#Set seed for reproducibility
set.seed(13)

#Read in data from riptide.maxfit_contextualize
cdmm_fluxes <- read.delim('~/RIPTIDE/riptide_cdmm_max_complete/flux_samples.tsv', sep='\t', header=TRUE, row.names=1)
iec_fluxes <- read.delim('~/RIPTIDE/riptide_iec_max_complete/flux_samples.tsv', sep='\t', header=TRUE, row.names=1)

#Subsample to 250 samplings, or minimum number for any condition.
sample_size <- min(c(nrow(cdmm_fluxes), nrow(iec_fluxes), 250))
sub_sample <- sample(1:min(c(nrow(cdmm_fluxes), nrow(iec_fluxes))), sample_size, replace=FALSE)
cdmm_fluxes <- cdmm_fluxes[sub_sample,]
iec_fluxes <- iec_fluxes[sub_sample,]
rm(sample_size, sub_sample)

## DETERMINE BIOMASS FLUX --------------------------------------------------------------------
cdmm_biomass <- as.vector(cdmm_fluxes[,'biomass']) 
iec_biomass <- as.vector(iec_fluxes[,'biomass']) 
mean(cdmm_biomass)
mean(iec_biomass)
sd(cdmm_biomass)
sd(iec_biomass)

#Test normality of distributions
shapiro.test(cdmm_biomass)$p.value #not normal
shapiro.test(iec_biomass)$p.value #not normal

#Test differences in biomass
pvals <- c()
for (x in c(1:1000)) {
  test_1 <- sample(cdmm_biomass, size=12)
  test_2 <- sample(iec_biomass, size=12)
  pvals[x] <- wilcox.test(test_1, test_2, exact=FALSE)$p.value}
biomass_pval <- round(median(pvals), 4) #p=0.000004; significant

#Biomass Vplot
biomass <- data.frame(cdmm_biomass, iec_biomass)
colnames(biomass) <- c("No Mucus", "IEC Mucus")
biomass <- melt(biomass)
colnames(biomass) <- c("Mucus", "Biomass")
biomass_vplot <- ggplot(biomass, aes(x=Mucus, y=Biomass, fill=Mucus)) +
  geom_violin(trim=FALSE) +
  theme_bw() +
  theme(legend.position='none') +
  stat_summary(fun = median, 
               geom = "point",
               color="white") +
  scale_fill_manual(values=c("No Mucus" = "#332288", "IEC Mucus" = "#5AAE61")) +
  annotate(geom="text", x= 1.5, y=max(biomass$Biomass), label="****", hjust="center", vjust="top", size=8)
biomass_vplot

ggsave(file="Biomass_Violin.svg", plot=biomass_vplot, width=4, height=4, dpi=300)

## DETERMINE UPTAKE/EFFLUX FOR METABOLITES OF INTEREST---------------------------------------------------
#Glycine-------------------------------------------------------
#see increased gly synthesis from serine and 5,10-methylene-THF with mucus
#cdmm_gly_uptake <- as.vector(cdmm_fluxes[,'EX_cpd00033_e']) * -1.0 #pruned; no uptake or efflux
#iec_gly_uptake <- as.vector(iec_fluxes[,'EX_cpd00033_e']) * -1.0 #pruned; no uptake or efflux

#Serine--------------------------------------------------------
#See increased conversion of serine to glycine in mucus
#cdmm_ser_uptake <- as.vector(cdmm_fluxes[,'EX_cpd00054_e']) * -1.0 #pruned; no uptake or efflux
iec_ser_uptake <- as.vector(iec_fluxes[,'EX_cpd00054_e']) * -1.0
mean(iec_ser_uptake) #7.415264
sd(iec_ser_uptake)
#Assess distribution
shapiro.test(iec_ser_uptake)$p.value #not normal

#Serine Vplot
cdmm_ser_null <- rep(0, 250)
serine <- data.frame(cdmm_ser_null, iec_ser_uptake)
colnames(serine) <- c("None", "IEC Mucus")
serine <- melt(serine)
colnames(serine) <- c("Mucus", "Uptake")
serine <- serine[-c(1:249),] #removes all but one sample; allows condition without mucus to be "plotted" without disrupting appearance of IEC Mucus violin distribution.
serine_vplot <- ggplot(serine, aes(x=Mucus, y=Uptake, fill=Mucus)) +
  geom_violin(trim=FALSE) +
  theme_bw() +
  theme(legend.position='none') +
  stat_summary(fun = median, 
               geom = "point",
               color="white") +
  scale_fill_manual(values=c("None" = "#332288", "IEC Mucus" = "#5AAE61")) +
  annotate(geom="text", x= 1, y=1.5, label="No Uptake", hjust="center", vjust="top", size=5)
serine_vplot

ggsave(file="Serine_Violin.svg", plot=serine_vplot, width=3, height=4, dpi=300)

#Aspartate---------------------------------------------------------
#See production of dihydrodipicolinate from L-Asp-4-semialdehyde; precursor to lysine biosynthesis
cdmm_asp_uptake <- as.vector(cdmm_fluxes[,'EX_cpd00041_e']) * -1.0
iec_asp_uptake <- as.vector(iec_fluxes[,'EX_cpd00041_e']) * -1.0
mean(cdmm_asp_uptake) #369.9667
mean(iec_asp_uptake) #207.3622
#Assess distributions
shapiro.test(cdmm_asp_uptake)$p.value #not normal
shapiro.test(iec_asp_uptake)$p.value #not normal
#Differences
pvals <- c()
for (x in c(1:1000)) {
  test_1 <- sample(cdmm_asp_uptake, size=12)
  test_2 <- sample(iec_asp_uptake, size=12)
  pvals[x] <- wilcox.test(test_1, test_2, exact=FALSE)$p.value}
asp_pval <- round(median(pvals), 4) #p=0.0351; significant!

#Aspartate Vplot
asp <- data.frame(cdmm_asp_uptake, iec_asp_uptake)
colnames(asp) <- c("No Mucus", "Mucus")
asp <- melt(asp)
colnames(asp) <- c("Mucus", "Uptake")
asp_vplot <- ggplot(asp, aes(x=Mucus, y=Uptake, fill=Mucus)) +
  geom_violin(trim=FALSE) +
  theme_bw() +
  theme(legend.position='none') +
  stat_summary(fun = median, 
               geom = "point",
               color="white") +
  scale_fill_manual(values=c("No Mucus" = "#332288", "Mucus" = "#5AAE61")) +
  annotate(geom="text", x= 1.5, y=max(asp$Uptake), label="*", hjust="center", vjust="top", size=8)
asp_vplot
ggsave(file="Asp_Violin.svg", plot=asp_vplot, width=3, height=4, dpi=300)

#L-glutamate------------------------------------------------------------- 
#see conversion from 2-oxoglutarate to NAD+ and glutamate with mucus
cdmm_glu_uptake <- as.vector(cdmm_fluxes[,'EX_cpd00023_e']) * -1.0
iec_glu_uptake <- as.vector(iec_fluxes[,'EX_cpd00023_e']) * -1.0
mean(cdmm_glu_uptake) #-37.99017; more L-glutamate efflux
mean(iec_glu_uptake) #28.56968; more L-glutamate uptake
#Assess distributions
shapiro.test(cdmm_glu_uptake)$p.value #not normal
shapiro.test(iec_glu_uptake)$p.value #not normal
#Differences
pvals <- c()
for (x in c(1:1000)) {
  test_1 <- sample(cdmm_glu_uptake, size=12)
  test_2 <- sample(iec_glu_uptake, size=12)
  pvals[x] <- wilcox.test(test_1, test_2, exact=FALSE)$p.value}
glu_pval <- round(median(pvals), 4) #p=0.0001; significant!

#Glutamate Vplot
glu <- data.frame(cdmm_glu_uptake, iec_glu_uptake)
colnames(glu) <- c("No Mucus", "Mucus")
glu <- melt(glu)
colnames(glu) <- c("Mucus", "Uptake")
glu_vplot <- ggplot(glu, aes(x=Mucus, y=Uptake, fill=Mucus)) +
  geom_violin(trim=FALSE) +
  theme_bw() +
  theme(legend.position='none') +
  stat_summary(fun = median, 
               geom = "point",
               color="white") +
  scale_fill_manual(values=c("No Mucus" = "#332288", "Mucus" = "#5AAE61")) +
  annotate(geom="text", x= 1.5, y=max(glu$Uptake), label="***", hjust="center", vjust="top", size=8)
glu_vplot
ggsave(file="Glu_Violin.svg", plot=glu_vplot, width=3, height=4, dpi=300)

#L-alanine-------------------------------------------------------
#see increased uptake and conversion to pyruvate with mucus
cdmm_ala_uptake <- as.vector(cdmm_fluxes[,'EX_cpd00035_e']) * -1.0
iec_ala_uptake <- as.vector(iec_fluxes[,'EX_cpd00035_e']) * -1.0
mean(cdmm_ala_uptake) #174.7038
mean(iec_ala_uptake) #674.6835
#Assess distributions
shapiro.test(cdmm_ala_uptake)$p.value #not normal
shapiro.test(iec_ala_uptake)$p.value #not normal
#Differences
pvals <- c()
for (x in c(1:1000)) {
  test_1 <- sample(cdmm_ala_uptake, size=12)
  test_2 <- sample(iec_ala_uptake, size=12)
  pvals[x] <- wilcox.test(test_1, test_2, exact=FALSE)$p.value}
ala_pval <- round(median(pvals), 4) #p=0.0001; significant!

#Alanine Vplot
ala <- data.frame(cdmm_ala_uptake, iec_ala_uptake)
colnames(ala) <- c("No Mucus", "Mucus")
ala <- melt(ala)
colnames(ala) <- c("Mucus", "Uptake")
ala_vplot <- ggplot(ala, aes(x=Mucus, y=Uptake, fill=Mucus)) +
  geom_violin(trim=FALSE) +
  theme_bw() +
  theme(legend.position='none') +
  stat_summary(fun = median, 
               geom = "point",
               color="white") +
  scale_fill_manual(values=c("No Mucus" = "#332288", "Mucus" = "#5AAE61")) +
  annotate(geom="text", x= 1.5, y=max(ala$Uptake), label="***", hjust="center", vjust="top", size=8)
ala_vplot
ggsave(file="Ala_Violin.svg", plot=ala_vplot, width=3, height=4, dpi=300)

#N-acetyl-D-glucosamine (GlcNAc)------------------------------------------------
#uptake only allowed for mucus condition, but see flux in pathways containing this metabolite both with and without mucus
cdmm_glcNAc_efflux <- as.vector(cdmm_fluxes[,'EX_cpd00122_e']) 
iec_glcNAc_efflux <- as.vector(iec_fluxes[,'EX_cpd00122_e']) 
mean(cdmm_glcNAc_efflux) #137.79
mean(iec_glcNAc_efflux) #711.7133
sd(cdmm_glcNAc_efflux)
sd(iec_glcNAc_efflux)
#Assess distributions
shapiro.test(cdmm_glcNAc_efflux)$p.value #not normal
shapiro.test(iec_glcNAc_efflux)$p.value #not normal
#Differences
pvals <- c()
for (x in c(1:1000)) {
  test_1 <- sample(cdmm_glcNAc_efflux, size=12)
  test_2 <- sample(iec_glcNAc_efflux, size=12)
  pvals[x] <- wilcox.test(test_1, test_2, exact=FALSE)$p.value}
glcNAc_pval <- round(median(pvals), 4) #p=0.0003; significant!

#GlcNAc Vplot
glcNAc <- data.frame(cdmm_glcNAc_efflux, iec_glcNAc_efflux)
colnames(glcNAc) <- c("No Mucus", "Mucus")
glcNAc <- melt(glcNAc)
colnames(glcNAc) <- c("Mucus", "Efflux")
glcNAc_vplot <- ggplot(glcNAc, aes(x=Mucus, y=Efflux, fill=Mucus)) +
  geom_violin(trim=FALSE) +
  theme_bw() +
  theme(legend.position='none') +
  stat_summary(fun = median, 
               geom = "point",
               color="white") +
  scale_fill_manual(values=c("No Mucus" = "#332288", "Mucus" = "#5AAE61")) +
  annotate(geom="text", x= 1.5, y=max(glcNAc$Efflux), label="***", hjust="center", vjust="top", size=8)
glcNAc_vplot
ggsave(file="GlcNAc_Violin.svg", plot=glcNAc_vplot, width=3, height=4, dpi=300)

#Proline: important stickland metabolite, also in mucin protein seq -------------------------- 
cdmm_pro_uptake <- as.vector(cdmm_fluxes[,'EX_cpd00129_e']) * -1.0
iec_pro_uptake <- as.vector(iec_fluxes[,'EX_cpd00129_e']) * -1.0
mean(cdmm_pro_uptake) #57.84925 #fits observations of increased proline stickland fermentation without mucus?
mean(iec_pro_uptake) #1.967352
sd(cdmm_pro_uptake)
sd(iec_pro_uptake)
#Assess distributions
shapiro.test(cdmm_pro_uptake)$p.value #not normal
shapiro.test(iec_pro_uptake)$p.value #not normal
#Differences
pvals <- c()
for (x in c(1:1000)) {
  test_1 <- sample(cdmm_pro_uptake, size=12)
  test_2 <- sample(iec_pro_uptake, size=12)
  pvals[x] <- wilcox.test(test_1, test_2, exact=FALSE)$p.value}
pro_pval <- round(median(pvals), 4) #p=0.014; significant!

#Pro Vplot
pro <- data.frame(cdmm_pro_uptake, iec_pro_uptake)
colnames(pro) <- c("No Mucus", "Mucus")
pro <- melt(pro)
colnames(pro) <- c("Mucus", "Uptake")
pro_vplot <- ggplot(pro, aes(x=Mucus, y=Uptake, fill=Mucus)) +
  geom_violin(trim=FALSE) +
  theme_bw() +
  theme(legend.position='none') +
  scale_y_continuous(trans="log10")+
  stat_summary(fun = median, 
               geom = "point",
               color="white") +
  scale_fill_manual(values=c("No Mucus" = "#332288", "Mucus" = "#5AAE61")) +
  annotate(geom="text", x= 1.5, y=max(pro$Uptake), label="*", hjust="center", vjust="top", size=8)
pro_vplot
ggsave(file="Pro_Violin.svg", plot=pro_vplot, width=3, height=4, dpi=300)

#Threonine is also prevalent in mucin protein seq --------------------------------------
cdmm_thr_uptake <- as.vector(cdmm_fluxes[,'EX_cpd00161_e']) * -1.0
iec_thr_uptake <- as.vector(iec_fluxes[,'EX_cpd00161_e']) * -1.0
mean(cdmm_thr_uptake) #2.316037
mean(iec_thr_uptake) #3.561337
sd(cdmm_thr_uptake)
sd(iec_thr_uptake)
#Assess distributions
shapiro.test(cdmm_thr_uptake)$p.value #not normal
shapiro.test(iec_thr_uptake)$p.value #not normal
#Differences
pvals <- c()
for (x in c(1:1000)) {
  test_1 <- sample(cdmm_thr_uptake, size=12)
  test_2 <- sample(iec_thr_uptake, size=12)
  pvals[x] <- wilcox.test(test_1, test_2, exact=FALSE)$p.value}
thr_pval <- round(median(pvals), 4) #p=0.0005; significant!
#Thr Vplot
thr <- data.frame(cdmm_thr_uptake, iec_thr_uptake)
colnames(thr) <- c("No Mucus", "Mucus")
thr <- melt(thr)
colnames(thr) <- c("Mucus", "Uptake")
thr_vplot <- ggplot(thr, aes(x=Mucus, y=Uptake, fill=Mucus)) +
  geom_violin(trim=FALSE) +
  theme_bw() +
  theme(legend.position='none') +
  stat_summary(fun = median, 
               geom = "point",
               color="white") +
  scale_fill_manual(values=c("No Mucus" = "#332288", "Mucus" = "#5AAE61")) +
  annotate(geom="text", x= 1.5, y=max(thr$Uptake), label="***", hjust="center", vjust="top", size=8)
thr_vplot
ggsave(file="thr_Violin.svg", plot=thr_vplot, width=3, height=4, dpi=300)

#Doubling Times--------------------------------------------------
#How to determine doubling times. Model has 3600 seconds defined as growth period, hence * 3600 used here.
cdmm_doubling <- (1/cdmm_biomass) * 3600.0
iec_doubling <- (1/iec_biomass) * 3600.0
mean(cdmm_doubling) #169.7202
mean(iec_doubling) #104.0111 C. difficile doubles faster with mucus.

#Doubling Time Vplot
dbl <- data.frame(cdmm_doubling, iec_doubling)
colnames(dbl) <- c("No Mucus", "Mucus")
dbl <- melt(dbl)
colnames(dbl) <- c("Mucus", "DoublingTime")
dbl_vplot <- ggplot(dbl, aes(x=Mucus, y=DoublingTime, fill=Mucus)) +
  geom_violin(trim=FALSE) +
  theme_bw() +
  stat_summary(fun = median, 
               geom = "point",
               color="white") +
  scale_fill_manual(values=c("No Mucus" = "#332288", "Mucus" = "#5AAE61")) +
  annotate(geom="text", x= 1.5, y=max(dbl$DoublingTime), label="***", hjust="center", vjust="top", size=8)
dbl_vplot

## CONFIRM BIOMASS DIFFERENCES ARE TRUE -----------------------------------------------------------------------------------
#Confirm differences are not due to inclusion of mucus monosaccharide exchanges in mucus condition only
set.seed(13)

# Read in data
cdmm_mono_fluxes <- read.delim('~/RIPTIDE/riptide_cdmm_monosacc_test_max/flux_samples.tsv', sep='\t', header=TRUE, row.names=1)

#Subsample IEC mucus condition to same number of samples as CDMM condition
sample_size <- min(c(nrow(cdmm_mono_fluxes), nrow(iec_fluxes), 250))
sub_sample <- sample(1:min(c(nrow(cdmm_mono_fluxes), nrow(iec_fluxes))), sample_size, replace=FALSE)
cdmm_mono_fluxes <- cdmm_mono_fluxes[sub_sample,]
iec_check_fluxes <- iec_fluxes[sub_sample,]
rm(sample_size, sub_sample)

#New Biomass
cdmm_mono_biomass <- as.vector(cdmm_mono_fluxes[,'biomass'])
iec_check_biomass <- as.vector(iec_check_fluxes[,'biomass'])
mean(cdmm_mono_biomass) #23.88783
mean(iec_check_biomass) #35.90058
sd(cdmm_mono_biomass)

#Confirm p-value:
pvals <- c()
for (x in c(1:1000)) {
  test_1 <- sample(cdmm_mono_biomass, size=12)
  test_2 <- sample(iec_check_biomass, size=12)
  pvals[x] <- wilcox.test(test_1, test_2, exact=FALSE)$p.value}
biomass_mono_pval <- round(median(pvals), 4) #p-0.0007, still significant!

#New violin plot
biomass_mono <- data.frame(cdmm_mono_biomass, iec_check_biomass)
colnames(biomass_mono) <- c("No Mucus", "IEC Mucus")
biomass_mono <- melt(biomass_mono)
colnames(biomass_mono) <- c("Mucus", "Biomass")

biomass_mono_vplot <- ggplot(biomass_mono, aes(x=Mucus, y=Biomass, fill=Mucus)) +
  geom_violin(trim=FALSE) +
  stat_summary(fun = median, 
               geom = "point",
               color="white") +
  theme_bw() +
  scale_fill_manual(values=c("No Mucus" = "#332288", "IEC Mucus" = "#5AAE61")) +
  annotate(geom="text", x= 1.5, y=max(biomass_mono$Biomass), label="***", hjust="center", vjust="top", size=8)
biomass_mono_vplot

