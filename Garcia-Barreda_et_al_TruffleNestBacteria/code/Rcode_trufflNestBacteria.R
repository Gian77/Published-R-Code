# ************ DATA ANALYSIS ********************************************** -----
# Project name: bacteria in truffle nests 
# Manuscript:   Bacterial communities show distinctive spatial diversity patterns
#               in productive truffle orchards amended with peat-based substrate
# Authors:      P Marco, S Sánchez, S Garcia-Barreda, J Parladé, M Rondolini,
#               V González, GMN Benucci, G Bonito
# Affiliation:  CITA, IRTA, Università di Perugia, ICA-CSIC, MSU
# Journal:      submitted to Environmental Microbiome
# Date:         June 16, 2020
# Code developer: Gian M. N. Benucci, Sergi Garcia-Barreda
# ************************************************************************* -----

# WORKING ENVIRONMENT SETUP -----------------
options(scipen = 9999) #to use decimals
options(max.print=100000000) # to print more lines on the display
setwd("path/to/workingdirectory")

# loading required packages -----------------
library(phyloseq); packageVersion("phyloseq")
library(Biostrings)
library(ggplot2)
library(ape)

library(dplyr)
library(tidyr)

# A) UPARSE ITS PAIRED -----------------------------
otus_ITS_uparse <- read.delim("ITS/otu_table_ITS_UPARSE_R1.txt", row.names=1) 
otus_phy_ITS_uparse <-otu_table(otus_ITS_uparse, taxa_are_rows = TRUE)

metadata_ITS_uparse <-read.delim("ITS/map_ITS_EXP2.txt", row.names=1, header=TRUE, sep="\t")
metadata_phy_ITS_uparse <-sample_data(metadata_ITS_uparse)

taxonomy_ITS_uparse_cons <-read.delim("ITS/consensus_taxonomy.txt", header=TRUE, row.names=1)
taxonomy_phy_ITS_uparse_cons <- tax_table(as.matrix(taxonomy_ITS_uparse_cons))

otus_seq_ITS_uparse <- readDNAStringSet("ITS/otus_R1.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

physeq_ITS_uparse <- phyloseq(otus_phy_ITS_uparse,
                              metadata_phy_ITS_uparse,
                              taxonomy_phy_ITS_uparse_cons,
                              otus_seq_ITS_uparse) 

physeq_ITS_uparse
tax_table(physeq_ITS_uparse)[tax_table(physeq_ITS_uparse)==""]<- NA
head(tax_table(physeq_ITS_uparse))
sample_data(physeq_ITS_uparse)

# Importing taxonomies at 0.6 confidence --------------------------------------
taxonomy_rdp <-read.delim("ITS/otu_taxonomy_rdp_final.txt", header=TRUE, row.names=1)
tail(taxonomy_rdp)
str(taxonomy_rdp)

taxonomy_blast <-read.delim("ITS/otu_taxonomy_blast_final.txt", header=TRUE, row.names=1)
tail(taxonomy_blast)

taxonomy_sintax <-read.delim("ITS/otu_taxonomy_sintax_final.txt", header=TRUE, row.names=1)
tail(taxonomy_sintax)

identical(rownames(taxonomy_rdp), rownames(taxonomy_blast))
taxonomy_all <- taxonomy_rdp
taxonomy_all <- taxonomy_all[,2:5]
colnames(taxonomy_all) <- c("Kingdom_rdp", "K_score_rdp", "Phylum_rdp", "P_score_rdp")
head(taxonomy_all)

taxonomy_all$Kingdom_blast <- taxonomy_blast$Kingdom
taxonomy_all$K_score_blast <- taxonomy_blast$K_score
taxonomy_all$Phylum_blast <- taxonomy_blast$Phylum
taxonomy_all$P_score_blast <- taxonomy_blast$P_score

taxonomy_all$Kingdom_sintax <- taxonomy_sintax$Kingdom
taxonomy_all$K_score_sintax <- taxonomy_sintax$K_score
taxonomy_all$Phylum_sintax <- taxonomy_sintax$Phylum
taxonomy_all$P_score_sintax <- taxonomy_sintax$P_score

head(taxonomy_all)

# Generating a new phyloseq object for easier manipulation --------
physeq_ITS <- phyloseq(otus_phy_ITS_uparse,
                       metadata_phy_ITS_uparse,
                       tax_table(as.matrix(taxonomy_all)),
                       otus_seq_ITS_uparse) 

physeq_ITS
head(tax_table(physeq_ITS))
tail(tax_table(physeq_ITS))

# Uing CONSTAX voting rules in here! 
# All OTUs classified as fungi (2 out 3 classifiers) ----------
physeq_fungi <- subset_taxa(physeq_ITS, Kingdom_rdp=="Fungi" & Kingdom_blast=="Fungi" |
                              Kingdom_rdp=="Fungi" & Kingdom_sintax=="Fungi" |
                              Kingdom_sintax=="Fungi" & Kingdom_blast=="Fungi")

# removing other taxa according RDP only
physeq_fungi <- subset_taxa(physeq_fungi, !(Kingdom_rdp=="Choanoflagellozoa" |
                                              Kingdom_rdp=="Euglenozoa" |
                                              Kingdom_rdp=="Ichthyosporia" |
                                              Kingdom_rdp=="Metazoa" |
                                              Kingdom_rdp=="Protista" |
                                              Kingdom_rdp=="Rhizaria"))

physeq_fungi
head(tax_table(physeq_fungi))
tail(tax_table(physeq_fungi))
any(tax_table(physeq_fungi)=="Euglenozoa")

sort(unique(as.data.frame(tax_table(physeq_fungi))$Kingdom_rdp))

# FILTERED FUNGI - PHYLOSEQ OBJECT --------------------------------------------
subset_taxa(physeq_ITS_uparse, 
            rownames(tax_table(physeq_ITS_uparse)) %in% 
              taxa_names(physeq_fungi)) -> physeq_fungi_uparse

physeq_fungi_uparse
head(tax_table(physeq_fungi_uparse))

# ADDITIONAL INFO --------------------------------------
# non-Fungi just according RDP

sort(unique(as.data.frame(tax_table(physeq_ITS))$Kingdom_rdp))

physeq_non_fungi <- subset_taxa(physeq_ITS, 
                                Kingdom_rdp=="Choanoflagellozoa" |
                                  Kingdom_rdp=="Euglenozoa" | 
                                  Kingdom_rdp=="Ichthyosporia" | 
                                  Kingdom_rdp=="Metazoa" |
                                  Kingdom_rdp=="Protista" |
                                  Kingdom_rdp=="Rhizaria")

physeq_non_fungi
head(tax_table(physeq_non_fungi))

# Fungi not classified at Phylum level -------------------------------------
# (2 out of 3 classifiers)
physeq_just_fungi <- subset_taxa(physeq_fungi, Phylum_rdp=="" & Phylum_blast=="" |
                                   Phylum_rdp=="" & Phylum_sintax=="" |
                                   Phylum_sintax=="" & Phylum_blast=="")
physeq_just_fungi
head(tax_table(physeq_just_fungi))

# Unclassified OTUs --------------------------------------------------------
# (2 out of 3 classifiers)
physeq_NA <- subset_taxa(physeq_ITS,  is.na(K_score_rdp) &  is.na(K_score_blast) |
                           is.na(K_score_rdp) &  is.na(K_score_sintax) |
                           is.na(K_score_blast) &  is.na(K_score_sintax))
physeq_NA
head(tax_table(physeq_NA))

# B) UPARSE 16S PAIRED -----------------------------------------
otus_16s_uparse <- read.delim("16S/otu_table_16s_UPARSE.txt",row.names=1) 
otus_phy_16s_uparse <- otu_table(otus_16s_uparse, taxa_are_rows = TRUE)

metadata_16s_uparse <- read.delim("16S/map_16s.txt", row.names=1, header=TRUE, sep="\t")
metadata_phy_16s_uparse <- sample_data(metadata_16s_uparse)

# importing RDP taxonomy ---------------------------------------
taxonomy_16s_uparse_RDP <- read.delim("16S/otus_taxonomy_stand_RDP.txt", header = TRUE, row.names = 1)
head(taxonomy_16s_uparse_RDP)

ifelse(taxonomy_16s_uparse_RDP$D_Score>=0.7, paste(taxonomy_16s_uparse_RDP$Domain), NA) -> Kingdom
ifelse(taxonomy_16s_uparse_RDP$P_Score>=0.7, paste(taxonomy_16s_uparse_RDP$Phylum), NA) -> Phylum
ifelse(taxonomy_16s_uparse_RDP$C_Score>=0.7, paste(taxonomy_16s_uparse_RDP$Class), NA) -> Class
ifelse(taxonomy_16s_uparse_RDP$O_Score>=0.7, paste(taxonomy_16s_uparse_RDP$Order), NA) -> Order
ifelse(taxonomy_16s_uparse_RDP$F_Score>=0.7, paste(taxonomy_16s_uparse_RDP$Family), NA) -> Family
ifelse(taxonomy_16s_uparse_RDP$G_Score>=0.7, paste(taxonomy_16s_uparse_RDP$Genus), NA) -> Genus
taxonomy <- cbind(Kingdom, Phylum, Class, Order, Family, Genus)
rownames(taxonomy) <- rownames(taxonomy_16s_uparse_RDP)
taxonomy_phy_16s_uparse_RDP <- tax_table(as.matrix(taxonomy))
head(taxonomy_phy_16s_uparse_RDP)

otus_seq_16s_uparse <- readDNAStringSet("16S/otus.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

physeq_obj_16s_uparse <- phyloseq(otus_phy_16s_uparse, 
                                  metadata_phy_16s_uparse,
                                  taxonomy_phy_16s_uparse_RDP,
                                  otus_seq_16s_uparse) 

physeq_obj_16s_uparse
head(tax_table(physeq_obj_16s_uparse))
tax_table(physeq_obj_16s_uparse)[tax_table(physeq_obj_16s_uparse)==""]<- NA
sample_data(physeq_obj_16s_uparse)

# looking for chloroplast in the RDP taxonomy 
any(tax_table(physeq_obj_16s_uparse) == "Chloroplast")
sort(unique(as.data.frame(tax_table(physeq_obj_16s_uparse))$Phylum))
physeq_obj_16s_uparse <- subset_taxa(physeq_obj_16s_uparse, Phylum != "Cyanobacteria/Chloroplast")

# looking for Chloroplast and Mitochondria in SILVA -------------------
library("dplyr")
library("tidyr")
library("stringr")

connect_16s_uparse_SILVA <-readLines("16s/otu_taxonomy_SILVA.txt")
connect_16s_uparse_SILVA <- gsub("\tk:.*\t\\+\tk:", ",", connect_16s_uparse_SILVA) 
head(connect_16s_uparse_SILVA)
str(connect_16s_uparse_SILVA)

max(count.fields(textConnection(connect_16s_uparse_SILVA), sep = ","))

taxonomy_16s_uparse_SILVA <-read.table(textConnection(connect_16s_uparse_SILVA),
                                       header = FALSE, 
                                       sep = ",",
                                       col.names = paste0("V",seq_len(8)),
                                       fill = TRUE)

rownames(taxonomy_16s_uparse_SILVA) <- taxonomy_16s_uparse_SILVA$V1
colnames(taxonomy_16s_uparse_SILVA) <- c("ZOTU_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
str(taxonomy_16s_uparse_SILVA)
head(taxonomy_16s_uparse_SILVA)

otu_chlo_mito <- rownames(taxonomy_16s_uparse_SILVA[taxonomy_16s_uparse_SILVA$Class=="c:Chloroplast" |
                                                      taxonomy_16s_uparse_SILVA$Family=="f:Mitochondria", ])
otu_chlo_mito

# REMOVING UNWANTED taxa ---------------------------------------------------
remove_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

physeq_prok_uparse <- remove_taxa(physeq_obj_16s_uparse, otu_chlo_mito)
physeq_prok_uparse

# * REMOVING CONTAMINANTS  -------------------------------------------------------
library(decontam)

# Fungi UPARSE -------------------------------------------------------------------------
# Inspect Library Sizes
df_fungi_uparse <- as.data.frame(as.matrix(sample_data(physeq_fungi_uparse)))
df_fungi_uparse$LibrarySize <- sample_sums(physeq_fungi_uparse)
df_fungi_uparse <- df_fungi_uparse[order(df_fungi_uparse$LibrarySize),]
df_fungi_uparse$Index <- seq(nrow(df_fungi_uparse)) # sample numbering
df_fungi_uparse

ggplot(data=df_fungi_uparse, aes(x=Index, y=LibrarySize, color=NULL)) + 
  geom_point() + 
  theme_classic()

ggplot(df_fungi_uparse, aes(x = LibrarySize)) + # Histogram of sample read counts
  theme_classic() +
  geom_histogram(color = "black", fill = "indianred", binwidth = 1000) +
  #facet_grid(~Treatment, scales = "free_x", space="free_x") +
  labs(title="Distribution of Sample Libraries", x="Read number", y="Sample number") + 
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  scale_color_manual(values="green")

# detecting contaminants 
sample_data(physeq_fungi_uparse)$is.neg <- sample_data(physeq_fungi_uparse)$Description == "control"
contam_prev_fungi <- isContaminant(physeq_fungi_uparse, method="prevalence", neg="is.neg", threshold=0.5)
table(contam_prev_fungi$contaminant)
head(which(contam_prev_fungi$contaminant))

# removing contaminants form the phyloseq object ----------------------
contaminants_fungi = rownames(subset(contam_prev_fungi, contaminant%in%c("TRUE")))
physeq_fungi_uparse_clean <- remove_taxa(physeq_fungi_uparse, contaminants_fungi)
physeq_fungi_uparse_clean
sample_data(physeq_fungi_uparse_clean)

# Now removing controls samples ----------------------------------------
physeq_fungi_uparse_clean <- subset_samples(physeq_fungi_uparse_clean, !Description%in%c("control"))
otu_table(physeq_fungi_uparse_clean) <- otu_table(physeq_fungi_uparse_clean)[which(rowSums(otu_table(physeq_fungi_uparse_clean)) > 0),] 
physeq_fungi_uparse_clean
sample_data(physeq_fungi_uparse_clean)

any(taxa_sums(physeq_fungi_uparse_clean) == 0)
sort(taxa_sums(physeq_fungi_uparse_clean))
any(sample_sums(physeq_fungi_uparse_clean) == 0)
sort(sample_sums(physeq_fungi_uparse_clean))

# Bacteria UPARSE -------------------------------------------------------
# Inspect Library Sizes
df_prok_uparse <- as.data.frame(as.matrix(sample_data(physeq_prok_uparse)))
df_prok_uparse$LibrarySize <- sample_sums(physeq_prok_uparse)
df_prok_uparse <- df_prok_uparse[order(df_prok_uparse$LibrarySize),]
df_prok_uparse$Index <- seq(nrow(df_prok_uparse)) # sample numbering
df_prok_uparse

ggplot(data=df_prok_uparse, aes(x=Index, y=LibrarySize, color=NULL)) + geom_point() + theme_classic()

ggplot(df_prok_uparse, aes(x = LibrarySize)) + # Histogram of sample read counts
  geom_histogram(color = "black", fill = "indianred", binwidth = 2000) +
  #facet_grid(~Treatment, scales = "free_x", space="free_x") +
  labs(title="Distribution of Sample Libraries", x="Read number", y="Sample number") + 
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  scale_color_manual(values="green")

# detecting contaminants 
sample_data(physeq_prok_uparse)$is.neg <- sample_data(physeq_prok_uparse)$Description == "control"
contam_prev_prokaryote <- isContaminant(physeq_prok_uparse, method="prevalence", neg="is.neg", threshold=0.5)
table(contam_prev_prokaryote$contaminant)
head(which(contam_prev_prokaryote$contaminant))
sample_data(physeq_prok_uparse)

# removing contaminants form the phyloseq object ----------------------
contaminants_prokaryote = rownames(subset(contam_prev_prokaryote, contaminant%in%c("TRUE")))
physeq_prok_uparse_clean <- remove_taxa(physeq_prok_uparse, contaminants_prokaryote)
physeq_prok_uparse_clean

# Now removing control samples
physeq_prok_uparse_clean <- subset_samples(physeq_prok_uparse_clean, !Description%in%c("control"))
otu_table(physeq_prok_uparse_clean) <- otu_table(physeq_prok_uparse_clean)[which(rowSums(otu_table(physeq_prok_uparse_clean)) > 0),] 
physeq_prok_uparse_clean
sample_data(physeq_prok_uparse_clean)

any(taxa_sums(physeq_prok_uparse_clean) == 0)
sort(taxa_sums(physeq_prok_uparse_clean))
any(sample_sums(physeq_prok_uparse_clean) == 0)
sort(sample_sums(physeq_prok_uparse_clean))

# FILTERED DATASTES following quality  -----------------------------
physeq_prok_uparse_clean
physeq_fungi_uparse_clean

sample_sums(physeq_fungi_uparse_clean)

physeq_fungi_uparse_qc <- physeq_fungi_uparse_clean
otu_table(physeq_fungi_uparse_qc) <- otu_table(physeq_fungi_uparse_qc)[which(rowSums(otu_table(physeq_fungi_uparse_qc)) > 0),]
physeq_fungi_uparse_qc
sample_sums(physeq_fungi_uparse_qc) 
any(taxa_sums(physeq_fungi_uparse_qc)==0)
any(sample_sums(physeq_fungi_uparse_qc)==0)

physeq_prok_uparse_qc <- physeq_prok_uparse_clean
otu_table(physeq_prok_uparse_qc) <- otu_table(physeq_prok_uparse_qc)[which(rowSums(otu_table(physeq_prok_uparse_qc)) > 0),] 
physeq_prok_uparse_qc
any(sample_sums(physeq_prok_uparse_qc)==0) 
any(taxa_sums(physeq_prok_uparse_qc)==0)

# removing samples with few reads -------------------
subset_samples(physeq_fungi_uparse_qc, sample_sums(physeq_fungi_uparse_qc) > 157) -> physeq_fungi_uparse_qc
otu_table(physeq_fungi_uparse_qc) <- otu_table(physeq_fungi_uparse_qc)[which(rowSums(otu_table(physeq_fungi_uparse_qc)) > 0),]

subset_samples(physeq_prok_uparse_qc, sample_sums(physeq_prok_uparse_qc) > 184) -> physeq_prok_uparse_qc
otu_table(physeq_prok_uparse_qc) <- otu_table(physeq_prok_uparse_qc)[which(rowSums(otu_table(physeq_prok_uparse_qc)) > 0),] 

# GOOD DATASETS to WORK WITH!! --------------------------------
physeq_fungi_uparse_qc
physeq_prok_uparse_qc

# DATA NORMALIZATION BEFORE ANALYSIS ----------------------

# TESTING CSS TRANSFORMATION ------------------------------
library(metagenomeSeq)

CSSNorm <-function(dataframe){
  require(metagenomeSeq)
  dataframe %>% 
    phyloseq_to_metagenomeSeq() -> physeq_CSS
  p_biom <-cumNormStatFast(physeq_CSS)
  biom_quant <-cumNorm(physeq_CSS, p=p_biom)
  physeq_CSS <- MRcounts(biom_quant, norm=T)
  physeq_mSeq <- dataframe
  otu_table(physeq_mSeq) <- otu_table(physeq_CSS, taxa_are_rows=TRUE)
  return(physeq_mSeq)
}

# MetagenomeSeq normalization - Gaussian model -----------------------------------------------------
CSSNorm(physeq_fungi_uparse_qc) -> physeq_fungi_uparse_mSeq
head(otu_table(physeq_fungi_uparse_mSeq))

CSSNorm(physeq_prok_uparse_qc) -> physeq_prok_uparse_mSeq

# PREPARATION FOR STATISTICAL ANALYSES
metadatos_bact <- read.delim("metadatos_bact.txt", row.names=1, sep = "\t", header=TRUE)
sample_data(physeq_prok_uparse_mSeq) <- metadatos_bact
head(sample_data(physeq_prok_uparse_mSeq))

# ALPHA-DIVERSITY
# only samples of the experiment
physeq_prok_uparse_mSeq_PCoA <- subset_samples(physeq_prok_uparse_mSeq, PCoA%in%c(1))
# for statistical analyses, without samples from year 0 and year 3
physeq_prok_uparse_mSeq_perman <- subset_samples(physeq_prok_uparse_mSeq_PCoA, anova%in%c(1))

# alpha-diversity
physeq_prok_uparse_mSeq_PCoA2 <- transform_sample_counts(physeq_prok_uparse_mSeq_PCoA, function(x) round(x, digits=0))
sample_data(physeq_prok_uparse_mSeq_PCoA2)$etiq <- factor(sample_data(physeq_prok_uparse_mSeq_PCoA2)$etiq,
                                                          levels = c("P0M", "P0T", "P1M", "P1T", "P2M", "P2T", "P3*",
                                                                     "S2M", "S2T"))
alpha_diversity <- estimate_richness(physeq_prok_uparse_mSeq_PCoA2, measure = c("Observed", "Shannon"))
df <- data.frame(alpha_diversity, sample_data(physeq_prok_uparse_mSeq_PCoA2))
df$Evenness <- df$Shannon/log(df$Observed)
df2 <- tidyr::gather(df, key = "Measure", value = "Value", Observed, Shannon, Evenness)

df2$Measure <- factor(df2$Measure)
levels(df2$Measure) <- list(Evenness="Evenness", Richness="Observed", Shannon="Shannon")
df2$Measure <- factor(df2$Measure, levels=c("Richness", "Evenness", "Shannon"))
df2$origin <- factor(df2$origin, levels=c("peat", "interphase", "soil"))

df3 <- subset(df2, etiq%in%c("P0M", "P0T", "P1M", "P1T", "P2M", "P2T", "S2M", "S2T"))
df3$Measure <- factor(df3$Measure)
levels(df3$Measure) <- list(Evenness="Evenness", Richness="Observed", Shannon="Shannon")
df3$Measure <- factor(df3$Measure, levels=c("Richness", "Evenness", "Shannon"))
df3$origin <- factor(df3$origin, levels=c("peat", "interphase", "soil"))

ggplot(data = df3, aes(x = etiq, y = Value)) +
  geom_jitter(position="identity", aes(color=origin)) +
  facet_wrap(~Measure, scale = "free") +
  geom_boxplot(fill=NA, outlier.color=NA) +
  ylab("Alpha Diversity Measure") +
  scale_color_manual(values=c("grey50", "red")) +
  scale_x_discrete(labels=c("N0M", "N0T", "N1M", "N1T", "N2M", "N2T", "S2M", "S2T")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none", axis.title.x=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black", linewidth=0.5, linetype = "solid"))

divers.index <- estimate_richness(physeq_prok_uparse_mSeq_PCoA2, measures = c("Observed", "Shannon"))
divers.index.2 <- cbind(sample_data(physeq_prok_uparse_mSeq_PCoA2), divers.index)
divers.perman2 <- subset(divers.index.2, anova%in%c(1))
physeq_prok_uparse_mSeq_perman2 <- subset_samples(physeq_prok_uparse_mSeq_PCoA2, anova%in%c(1))
richness.anova <- aov(Observed ~ location*year, divers.perman2)
# non significant predictors are removed from the model
richness.anova <- aov(Observed ~ location, divers.perman2)
summary(richness.anova)
# Table S3
par(mfcol=c(2,2))
plot(richness.anova)
shapiro.test(residuals(richness.anova))

# mixed model for richness
library(nlme)
lme.rich0 <- gls(Observed ~ location*year,
                 data = divers.perman2)
lme.rich1 <- lme(Observed ~ location*year,
                 random= ~1|tree_id, data=divers.perman2)
anova(lme.rich0, lme.rich1)
# AIC min in lme.rich0

# mixed models for Shannon index
lme.Shannon0 <- gls(Shannon ~ location*year,
                    data = divers.perman2)
lme.Shannon1 <- lme(Shannon ~ location*year,
                    random= ~1|tree_id, data=divers.perman2)
anova(lme.Shannon0, lme.Shannon1)

shannon.anova <- aov(Shannon ~ year*location, divers.perman2)
summary(shannon.anova)
par(mfcol=c(2,2))
plot(shannon.anova)
shapiro.test(residuals(shannon.anova))

# evenness
divers.perman2$evenness <- divers.perman2$Shannon/log(divers.perman2$Observed)
lme.even0 <- gls(evenness ~ location*year,
                 data = divers.perman2)
lme.even1 <- lme(evenness ~ location*year,
                 random= ~1|tree_id, data=divers.perman2)
anova(lme.even0, lme.even1)

even.anova <- aov(evenness ~ year*location, divers.perman2)
summary(even.anova)
par(mfcol=c(2,2))
plot(even.anova)
shapiro.test(residuals(even.anova))

# BETA-DIVERSITY
library(ggrepel)
physeq_prok_uparse_mSeq_PCoA %>% ordinate(method ="PCoA", distance="bray") -> pcoa_prok_mSeq_PCoA
pcoa_prok_mSeq_PCoA$values$Eigenvalues
# negative eigenvalues are not admited

physeq_prok_uparse_mSeq_PCoAsinP3_noround <- subset_samples(physeq_prok_uparse_mSeq_PCoA, origin%in%c("peat", "soil"))

physeq_prok_uparse_mSeq_PCoAsinP3_noround %>% ordinate(method ="PCoA", distance="bray") -> pcoa_prok_mSeq_PCoA
pcoa_prok_mSeq_PCoA$values$Eigenvalues
# negative eigenvalues are not admited
bray_dist <- phyloseq::distance(physeq_prok_uparse_mSeq_PCoAsinP3_noround, "bray")
ordinate(physeq_prok_uparse_mSeq_PCoAsinP3_noround, method ="PCoA", distance=sqrt(bray_dist)) -> pcoa_prok_mSeq_PCoA
plot_scree(pcoa_prok_mSeq_PCoA)

# Figure 3
library(glue)
pretty_pe <- format(round(100*pcoa_prok_mSeq_PCoA$values$Relative_eig[1:2], digits=1), nsmall=1, trim=T)
labs <- c(glue("PCoA 1 ({pretty_pe[1]}%)"), glue("PCoA 2 ({pretty_pe[2]}%)"))

PlotOrdin <-function(dataframe, ord){
  ord <- plot_ordination(dataframe, ord) + 
    geom_point(size=2, aes(color=origin, shape=year, alpha=location)) +
    theme_classic() +
    scale_color_manual("Soil position" ,values=c("blue3", "red"), labels=c("Truffle nests", "Bulk soil")) +
    scale_shape_manual("Year", values=c(1, 15, 17)) +
    scale_alpha_manual("Plot", values=c(1, 0.3)) +
    xlab(labs[1]) + ylab(labs[2]) +
    guides(shape = guide_legend(order = 2), color = guide_legend(order = 1))
  return(ord)
}

p <- PlotOrdin(physeq_prok_uparse_mSeq_PCoAsinP3_noround, pcoa_prok_mSeq_PCoA)
p$layers <- p$layers[-1]
p + ggpubr::grids(linetype = "dashed")

# permanova
otu_prok <- as.data.frame(otu_table(physeq_prok_uparse_mSeq_perman))
metadata_prok <- as.data.frame(as.matrix(sample_data(physeq_prok_uparse_mSeq_perman)[,2:16]))
identical(colnames(otu_prok), rownames(metadata_prok))
metadata_prok$mycelium <- as.numeric(metadata_prok$mycelium)

library(vegan)
prok_jc <- as.matrix(otu_table(physeq_prok_uparse_mSeq_perman))
vegan::vegdist(t(prok_jc), method="bray") -> dist_prok_jc
adonis2(dist_prok_jc ~ location * year, data=metadata_prok, permutations=9999, sqrt.dist=T, by="terms")

# COMMUNITY COMPOSITION AT PHYLUM/GENUS LEVEL (Fig. 1a)
new.etiq <- c("P0", "P0", "P1-2", "P1-2", "P1-2", "P1-2",
              "P3*", "S", "S")
sample_data(physeq_prok_uparse_mSeq_PCoA2)$etiq2 <- factor(new.etiq[sample_data(physeq_prok_uparse_mSeq_PCoA2)$etiq])
physeq_prok_uparse_mSeq_PCoAsinP3 <- subset_samples(physeq_prok_uparse_mSeq_PCoA2, etiq2%in%c("P0", "P1-2", "S"))

library(microbiome)
ps1.rel <- microbiome::transform(physeq_prok_uparse_mSeq_PCoAsinP3, "compositional")
ps1.phy.rel <-aggregate_rare(ps1.rel, level = "Phylum", detection = 0.005, prevalence = 0.5)
# igual resultat si canviem a detection=0.01, prevalence=0.3

dat.trans = transform_sample_counts(ps1.phy.rel, function(x) 100*x/sum(x))
dat.phyl = psmelt(dat.trans)
dat.agr = aggregate(Abundance~etiq+Phylum, data=dat.phyl, FUN=mean)
dat.agr2 <- subset(dat.agr, etiq%in%c("P0M", "P0T", "P1M", "P1T", "P2M", "P2T", "S2M", "S2T"))
dat.agr2$Phylum <- factor(dat.agr2$Phylum, levels = c("Acidobacteria", "Actinobacteria", "Bacteroidetes", "Chloroflexi", "Gemmatimonadetes", "Planctomycetes", "Proteobacteria", "Thaumarchaeota", "Verrucomicrobia", "Other"))
aggregate(dat.agr2$Abundance, list(dat.agr2$Phylum), FUN=mean)
getPalette2 = brewer.pal(10, "Set3")
ggplot(dat.agr2, aes(x=etiq, y=Abundance, fill=Phylum)) +
  geom_bar(stat="identity") +
  ggtitle("(a)") +
  scale_y_continuous(name="Relative abundance (%)", limits=c(0, 100), expand = c(0.01, 0)) +
  scale_x_discrete(labels=c("N0M", "N0T", "N1M", "N1T", "N2M", "N2T", "S2M", "S2T")) +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black", linewidth = 1, linetype = "solid"), plot.title=element_text(face="bold"))+
  guides(fill=guide_legend(ncol =1)) +
  scale_fill_manual(values=getPalette2)

aggregate(dat.agr2$Abundance, list(dat.agr2$Phylum, dat.agr2$etiq), FUN=mean)

# genus level (Fig. 1b)
ps1.gen.rel <-aggregate_rare(ps1.rel, level = "Genus", detection = 0.01, prevalence = 0.10)
dat.trans = transform_sample_counts(ps1.gen.rel, function(x) 100*x/sum(x))
dat.dataframe = psmelt(dat.trans)
dat.agr = aggregate(Abundance~etiq+Genus, data=dat.dataframe, FUN=mean)
dat.agr2 <- subset(dat.agr, etiq%in%c("P0M", "P0T", "P1M", "P1T", "P2M", "P2T", "S2M", "S2T"))
dat.agr3 <- dat.agr2[!(dat.agr2$Genus == "Other"),]
dat.agr3 <- dat.agr3[!(dat.agr3$Genus == "Unknown"),]
dat.agr3$Genus[dat.agr3$Genus=="Gp16"]<-"Acidobacteria subdiv. 16"
dat.agr3$Genus[dat.agr3$Genus=="Gp17"]<-"Acidobacteria subdiv. 17"
dat.agr3$Genus[dat.agr3$Genus=="Gp4"]<-"Acidobacteria subdiv. 4"
dat.agr3$Genus[dat.agr3$Genus=="Gp6"]<-"Acidobacteria subdiv. 6"
aggregate(dat.agr3$Abundance, list(dat.agr3$Genus), FUN=mean)

getPalette = colorRampPalette(brewer.pal(25, "Paired"))
ggplot(dat.agr3, aes(x=etiq, y=Abundance, fill=Genus)) +
  geom_bar(stat="identity") +
  ggtitle("(b)") +
  scale_y_continuous(name="Relative abundance (%)", limits=c(0, 50), expand = c(0.01, 0)) +
  scale_x_discrete(labels=c("N0M", "N0T", "N1M", "N1T", "N2M", "N2T", "S2M", "S2T")) +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black", linewidth = 1, linetype = "solid"), plot.title=element_text(face="bold"))+
  guides(fill=guide_legend(ncol =1)) +
  scale_fill_manual(values=getPalette(25))

# Rel. freq at order level
tax_table(physeq_prok_uparse_mSeq_perman)[1:5, 1:4]
actino <- subset_taxa(physeq_prok_uparse_mSeq_perman, Phylum=="Actinobacteria")
proteo <- subset_taxa(physeq_prok_uparse_mSeq_perman, Phylum=="Proteobacteria")
acido <- subset_taxa(physeq_prok_uparse_mSeq_perman, Phylum=="Acidobacteria")
bacteroi <- subset_taxa(physeq_prok_uparse_mSeq_perman, Phylum=="Bacteroidetes")
sum(sample_sums(actino))/sum(sample_sums(physeq_prok_uparse_mSeq_perman))
sum(sample_sums(proteo))/sum(sample_sums(physeq_prok_uparse_mSeq_perman))
sum(sample_sums(acido))/sum(sample_sums(physeq_prok_uparse_mSeq_perman))
sum(sample_sums(bacteroi))/sum(sample_sums(physeq_prok_uparse_mSeq_perman))

#Venn diagrams
# P0 vs P1-P2 vs Soil
new.etiq <- c("P0", "P0", "P1-2", "P1-2", "P1-2", "P1-2",
              NA, "S", "S")
sample_data(physeq_prok_uparse_mSeq_PCoA2)$etiq2 <- factor(new.etiq[sample_data(physeq_prok_uparse_mSeq_PCoA2)$etiq])

physeq_pcoa_etiq <- merge_samples(physeq_prok_uparse_mSeq_PCoA2, "etiq2")
otu_table_etiq <- as.data.frame(t(otu_table(subset_samples(physeq_pcoa_etiq))))
otu_table_etiq <- otu_table_etiq[(otu_table_etiq[,1] + otu_table_etiq[,2] + otu_table_etiq[,3]) > 0, ]
#colnames(otu_table_etiq) <- c("Peat substrate (n = 2)", "Nests (n = 45)", "Bulk soil (n = 4)")
venn_counts_etiq <- vennCounts(otu_table_etiq, include="both")
venn_counts_etiq
vennDiagram(venn_counts_etiq, cex=c(0.8,0.8,0.8))
# interes: nest mes en comu amb soil o substrat...

colnames(otu_table_etiq) <- c("Peat", "Nests", "Soil")
otu_table_etiq$OTU <- row.names(otu_table_etiq)
attach(otu_table_etiq)
otu_table_etiq$Peat[Peat > 0] <- otu_table_etiq$OTU[Peat > 0]
otu_table_etiq$Peat[Peat == 0] <- NA
otu_table_etiq$Nests[Nests > 0] <- otu_table_etiq$OTU[Nests > 0]
otu_table_etiq$Nests[Nests == 0] <- NA
otu_table_etiq$Soil[Soil > 0] <- otu_table_etiq$OTU[Soil > 0]
otu_table_etiq$Soil[Soil == 0] <- NA
detach(otu_table_etiq)
otu_table_etiq <- otu_table_etiq[,-4]
mylist <- list(Peat=otu_table_etiq$Peat[!is.na(otu_table_etiq$Peat)],
               Nest=otu_table_etiq$Nests[!is.na(otu_table_etiq$Nests)],
               Soil=otu_table_etiq$Soil[!is.na(otu_table_etiq$Soil)])
library(ggVennDiagram)
venn1 <- ggVennDiagram(mylist, label_alpha = 0, category.names = c("Peat\n(n=2)", "Nests\n(n=45)", "Bulk soil (n=4)"), label = c("count"), set_size=5) +
  ggplot2::scale_fill_gradient(low="white",high = "red") +
  theme(legend.position = "none")
venn1

# P0 vs P1 vs P2
new.etiq <- c("P0", "P0", "P1", "P1", "P2", "P2", NA, NA, NA)
sample_data(physeq_prok_uparse_mSeq_PCoA2)$etiq2 <- factor(new.etiq[sample_data(physeq_prok_uparse_mSeq_PCoA2)$etiq])

physeq_pcoa_etiq <- merge_samples(physeq_prok_uparse_mSeq_PCoA2, "etiq2")
otu_table_etiq <- as.data.frame(t(otu_table(subset_samples(physeq_pcoa_etiq))))
otu_table_etiq <- otu_table_etiq[(otu_table_etiq[,1] + otu_table_etiq[,2] + otu_table_etiq[,3]) > 0, ]
# colnames(otu_table_etiq) <- c("Peat substrate (n = 2)", "Nests year 1 (n = 22)", "Nests year 2 (n = 23)")
venn_counts_etiq <- vennCounts(otu_table_etiq, include="both")
venn_counts_etiq
vennDiagram(venn_counts_etiq, cex=c(0.8,0.8))

colnames(otu_table_etiq) <- c("Peat", "Nests1", "Nests2")
otu_table_etiq$OTU <- row.names(otu_table_etiq)
attach(otu_table_etiq)
otu_table_etiq$Peat[Peat > 0] <- otu_table_etiq$OTU[Peat > 0]
otu_table_etiq$Peat[Peat == 0] <- NA
otu_table_etiq$Nests1[Nests1 > 0] <- otu_table_etiq$OTU[Nests1 > 0]
otu_table_etiq$Nests1[Nests1 == 0] <- NA
otu_table_etiq$Nests2[Nests2 > 0] <- otu_table_etiq$OTU[Nests2 > 0]
otu_table_etiq$Nests2[Nests2 == 0] <- NA
detach(otu_table_etiq)
otu_table_etiq <- otu_table_etiq[,-4]
mylist <- list(Peat=otu_table_etiq$Peat[!is.na(otu_table_etiq$Peat)],
               Nests1=otu_table_etiq$Nests1[!is.na(otu_table_etiq$Nests1)],
               Nests2=otu_table_etiq$Nests2[!is.na(otu_table_etiq$Nests2)])

venn2 <- ggVennDiagram(mylist, label_alpha = 0, category.names = c("Peat\n(n=2)", "Nests year 1\n(n=22)", "Nests year 2 (n=23)"), label = c("count"), set_size=5) +
  ggplot2::scale_fill_gradient(low="white",high = "red") +
  theme(legend.position = "none")
venn2

# P0 vs PM vs PT
new.etiq <- c("P0", "P0", "PM", "PT", "PM", "PT", NA, NA, NA)
sample_data(physeq_prok_uparse_mSeq_PCoA2)$etiq2 <- factor(new.etiq[sample_data(physeq_prok_uparse_mSeq_PCoA2)$etiq])

physeq_pcoa_etiq <- merge_samples(physeq_prok_uparse_mSeq_PCoA2, "etiq2")
otu_table_etiq <- as.data.frame(t(otu_table(subset_samples(physeq_pcoa_etiq))))
otu_table_etiq <- otu_table_etiq[(otu_table_etiq[,1] + otu_table_etiq[,2] + otu_table_etiq[,3]) > 0, ]
#colnames(otu_table_etiq) <- c("Peat substrate (n = 2)", "Nests in Mora plantation (n = 20)", "Nests in Teruel plantation (n = 25)")
venn_counts_etiq <- vennCounts(otu_table_etiq, include="both")
venn_counts_etiq
vennDiagram(venn_counts_etiq, cex=c(0.8,0.8))

colnames(otu_table_etiq) <- c("Peat", "NestsM", "NestsT")
otu_table_etiq$OTU <- row.names(otu_table_etiq)
attach(otu_table_etiq)
otu_table_etiq$Peat[Peat > 0] <- otu_table_etiq$OTU[Peat > 0]
otu_table_etiq$Peat[Peat == 0] <- NA
otu_table_etiq$NestsM[NestsM > 0] <- otu_table_etiq$OTU[NestsM > 0]
otu_table_etiq$NestsM[NestsM == 0] <- NA
otu_table_etiq$NestsT[NestsT > 0] <- otu_table_etiq$OTU[NestsT > 0]
otu_table_etiq$NestsT[NestsT == 0] <- NA
detach(otu_table_etiq)
otu_table_etiq <- otu_table_etiq[,-4]
mylist <- list(Peat=otu_table_etiq$Peat[!is.na(otu_table_etiq$Peat)],
               NestsM=otu_table_etiq$NestsM[!is.na(otu_table_etiq$NestsM)],
               NestsT=otu_table_etiq$NestsT[!is.na(otu_table_etiq$NestsT)])

venn3 <- ggVennDiagram(mylist, label_alpha = 0, category.names = c("Peat\n(n=2)", "Nests Mora\n(n=20)", "Nests Teruel (n=25)"), label = c("count"), set_size=5) +
  ggplot2::scale_fill_gradient(low="white",high = "red") +
  theme(legend.position = "none")
venn3

# S vs P1 vs P2
new.etiq <- c(NA, NA, "P1", "P1", "P2", "P2", NA, "S", "S")
sample_data(physeq_prok_uparse_mSeq_PCoA2)$etiq2 <- factor(new.etiq[sample_data(physeq_prok_uparse_mSeq_PCoA2)$etiq])

physeq_pcoa_etiq <- merge_samples(physeq_prok_uparse_mSeq_PCoA2, "etiq2")
otu_table_etiq <- as.data.frame(t(otu_table(subset_samples(physeq_pcoa_etiq))))
otu_table_etiq <- otu_table_etiq[(otu_table_etiq[,1] + otu_table_etiq[,2] + otu_table_etiq[,3]) > 0, ]
#colnames(otu_table_etiq) <- c("Nests year 1 (n = 22)", "Nests year 2 (n = 23)", "Bulk soil (n = 4)")
venn_counts_etiq <- vennCounts(otu_table_etiq, include="both")
venn_counts_etiq
vennDiagram(venn_counts_etiq, cex=c(0.8,0.8))

colnames(otu_table_etiq) <- c("Nests1", "Nests2", "Soil")
otu_table_etiq$OTU <- row.names(otu_table_etiq)
attach(otu_table_etiq)
otu_table_etiq$Soil[Soil > 0] <- otu_table_etiq$OTU[Soil > 0]
otu_table_etiq$Soil[Soil == 0] <- NA
otu_table_etiq$Nests1[Nests1 > 0] <- otu_table_etiq$OTU[Nests1 > 0]
otu_table_etiq$Nests1[Nests1 == 0] <- NA
otu_table_etiq$Nests2[Nests2 > 0] <- otu_table_etiq$OTU[Nests2 > 0]
otu_table_etiq$Nests2[Nests2 == 0] <- NA
detach(otu_table_etiq)
otu_table_etiq <- otu_table_etiq[,-4]
mylist <- list(Soil=otu_table_etiq$Soil[!is.na(otu_table_etiq$Soil)],
               Nests1=otu_table_etiq$Nests1[!is.na(otu_table_etiq$Nests1)],
               Nests2=otu_table_etiq$Nests2[!is.na(otu_table_etiq$Nests2)])

venn4 <- ggVennDiagram(mylist, label_alpha = 0, category.names = c("Bulk soil\n(n=4)", "Nests year 1\n(n=22)", "Nests year 2 (n=23)"), label = c("count"), set_size=5) +
  ggplot2::scale_fill_gradient(low="white",high = "red") +
  theme(legend.position = "none")
venn4

# S vs PM vs PT
new.etiq <- c(NA, NA, "PM", "PT", "PM", "PT", NA, "S", "S")
sample_data(physeq_prok_uparse_mSeq_PCoA2)$etiq2 <- factor(new.etiq[sample_data(physeq_prok_uparse_mSeq_PCoA2)$etiq])

physeq_pcoa_etiq <- merge_samples(physeq_prok_uparse_mSeq_PCoA2, "etiq2")
otu_table_etiq <- as.data.frame(t(otu_table(subset_samples(physeq_pcoa_etiq))))
otu_table_etiq <- otu_table_etiq[(otu_table_etiq[,1] + otu_table_etiq[,2] + otu_table_etiq[,3]) > 0, ]
#colnames(otu_table_etiq) <- c("Nests in Mora plantation (n = 20)", "Nests in Teruel plantation (n = 25)", "Bulk soil (n = 4)")
venn_counts_etiq <- vennCounts(otu_table_etiq, include="both")
venn_counts_etiq
vennDiagram(venn_counts_etiq, cex=c(0.8,0.8))

colnames(otu_table_etiq) <- c("NestsM", "NestsT", "Soil")
otu_table_etiq$OTU <- row.names(otu_table_etiq)
attach(otu_table_etiq)
otu_table_etiq$Soil[Soil > 0] <- otu_table_etiq$OTU[Soil > 0]
otu_table_etiq$Soil[Soil == 0] <- NA
otu_table_etiq$NestsM[NestsM > 0] <- otu_table_etiq$OTU[NestsM > 0]
otu_table_etiq$NestsM[NestsM == 0] <- NA
otu_table_etiq$NestsT[NestsT > 0] <- otu_table_etiq$OTU[NestsT > 0]
otu_table_etiq$NestsT[NestsT == 0] <- NA
detach(otu_table_etiq)
otu_table_etiq <- otu_table_etiq[,-4]
mylist <- list(Soil=otu_table_etiq$Soil[!is.na(otu_table_etiq$Soil)],
               NestsM=otu_table_etiq$NestsM[!is.na(otu_table_etiq$NestsM)],
               NestsT=otu_table_etiq$NestsT[!is.na(otu_table_etiq$NestsT)])

venn5 <- ggVennDiagram(mylist, label_alpha = 0, category.names = c("Bulk soil\n(n=4)", "Nests Mora\n(n=20)", "Nests Teruel (n=25)"), label = c("count"), set_size=5) +
  ggplot2::scale_fill_gradient(low="white",high = "red") +
  theme(legend.position = "none")
venn5

# CORRELATIONS WITH TRUFFLE MYCELIUM
physeq_prok_uparse_mSeq_PCoA_percent <- transform_sample_counts(physeq_prok_uparse_mSeq_PCoA, function(x) x / sum(x) )
physeq_prok_uparse_mSeq_perman_percent <- subset_samples(physeq_prok_uparse_mSeq_PCoA_percent, anova%in%c(1))
dat.fam.brady <- subset_taxa(physeq_prok_uparse_mSeq_perman_percent,
                             Family=="Bradyrhizobiaceae")
tax_table(dat.fam.brady)
dat.fam.brady.merged <- merge_taxa(dat.fam.brady, taxa_names(dat.fam.brady))
mel.myc <- data.frame(sample_data(dat.fam.brady.merged)[,15])
mel.myc$n_sample <- rownames(mel.myc)

divers.index$n_sample <- rownames(divers.index)
divers.mycel <- merge(x = divers.index, y = mel.myc, by = "n_sample", all = TRUE)
colnames(divers.mycel)[2:3] <- c("prok_rich", "prok_Shannon")

lm.prok.rich <- lm(prok_rich ~ log(mycelium+0.001), divers.mycel)
summary(lm.prok.rich)
par(mfcol=c(2,2))
plot(lm.prok.rich)
shapiro.test(residuals(lm.prok.rich))

cor.test(divers.mycel$prok_rich, log(divers.mycel$mycelium+0.001))

lm.prok.shannon <- lm(prok_Shannon ~ log(mycelium+0.001), divers.mycel)
summary(lm.prok.shannon)
par(mfcol=c(2,2))
plot(lm.prok.shannon)
shapiro.test(residuals(lm.prok.shannon))

cor.test(divers.mycel$prok_Shannon, log(divers.mycel$mycelium+0.001))


actino2 <- transform_sample_counts(actino, function(x) round(x, digits=0))
divers.actino <- estimate_richness(actino2, measure = c("Observed", "Shannon"))
divers.actino$n_sample <- rownames(divers.actino)
actino.mycel <- merge(x = divers.actino, y = mel.myc, by = "n_sample", all = TRUE)
colnames(actino.mycel)[2:3] <- c("prok_rich", "prok_Shannon")

lm.actino.rich <- lm(sqrt(prok_rich+250) ~ log(mycelium+0.001), actino.mycel)
summary(lm.actino.rich)
par(mfcol=c(2,2))
plot(lm.actino.rich)
shapiro.test(residuals(lm.actino.rich))

cor.test(sqrt(actino.mycel$prok_rich+250), log(actino.mycel$mycelium+0.001))

lm.actino.shannon <- lm(prok_Shannon ~ log(mycelium+0.001), actino.mycel)
summary(lm.actino.shannon)
par(mfcol=c(2,2))
plot(lm.actino.shannon)
shapiro.test(residuals(lm.actino.shannon))

cor.test(actino.mycel$prok_Shannon, log(actino.mycel$mycelium+0.001))


proteo2 <- transform_sample_counts(proteo, function(x) round(x, digits=0))
divers.proteo <- estimate_richness(proteo2, measure = c("Observed", "Shannon"))
divers.proteo$n_sample <- rownames(divers.proteo)
proteo.mycel <- merge(x = divers.proteo, y = mel.myc, by = "n_sample", all = TRUE)
colnames(proteo.mycel)[2:3] <- c("prok_rich", "prok_Shannon")

lm.proteo.rich <- lm(sqrt(prok_rich+350) ~ log(mycelium+0.001), proteo.mycel)
summary(lm.proteo.rich)
par(mfcol=c(2,2))
plot(lm.proteo.rich)
shapiro.test(residuals(lm.proteo.rich))

cor.test(sqrt(proteo.mycel$prok_rich+350), log(proteo.mycel$mycelium+0.001))

lm.proteo.shannon <- lm(prok_Shannon ~ log(mycelium+0.001), proteo.mycel)
summary(lm.proteo.shannon)
par(mfcol=c(2,2))
plot(lm.proteo.shannon)
shapiro.test(residuals(lm.proteo.shannon))

cor.test(proteo.mycel$prok_Shannon, log(proteo.mycel$mycelium+0.001))


acido2 <- transform_sample_counts(acido, function(x) round(x, digits=0))
divers.acido <- estimate_richness(acido2, measure = c("Observed", "Shannon"))
divers.acido$n_sample <- rownames(divers.acido)
acido.mycel <- merge(x = divers.acido, y = mel.myc, by = "n_sample", all = TRUE)
colnames(acido.mycel)[2:3] <- c("prok_rich", "prok_Shannon")

lm.acido.rich <- lm(prok_rich ~ log(mycelium+0.001), acido.mycel)
# lm.acido.rich <- lm(sqrt(prok_rich+150) ~ log(mycelium+0.001), acido.mycel)
summary(lm.acido.rich)
par(mfcol=c(2,2))
plot(lm.acido.rich)
shapiro.test(residuals(lm.acido.rich))

cor.test(acido.mycel$prok_rich, log(acido.mycel$mycelium+0.001))

lm.acido.shannon <- lm(prok_Shannon ~ log(mycelium+0.001), acido.mycel)
summary(lm.acido.shannon)
par(mfcol=c(2,2))
plot(lm.acido.shannon)
shapiro.test(residuals(lm.acido.shannon))

cor.test(acido.mycel$prok_Shannon, log(acido.mycel$mycelium+0.001))


bacteroi2 <- transform_sample_counts(bacteroi, function(x) round(x, digits=0))
divers.bacteroi <- estimate_richness(actino2, measure = c("Observed", "Shannon"))
divers.bacteroi$n_sample <- rownames(divers.bacteroi)
bacteroi.mycel <- merge(x = divers.bacteroi, y = mel.myc, by = "n_sample", all = TRUE)
colnames(bacteroi.mycel)[2:3] <- c("prok_rich", "prok_Shannon")

lm.bacteroi.rich <- lm(prok_rich ~ log(mycelium+0.001), bacteroi.mycel)
# lm.bacteroi.rich <- lm(sqrt(prok_rich+250) ~ log(mycelium+0.001), bacteroi.mycel)
summary(lm.bacteroi.rich)
par(mfcol=c(2,2))
plot(lm.bacteroi.rich)
shapiro.test(residuals(lm.bacteroi.rich))

cor.test(bacteroi.mycel$prok_rich, log(bacteroi.mycel$mycelium+0.001))

lm.bacteroi.shannon <- lm(prok_Shannon ~ log(mycelium+0.001), bacteroi.mycel)
summary(lm.bacteroi.shannon)
par(mfcol=c(2,2))
plot(lm.bacteroi.shannon)
shapiro.test(residuals(lm.bacteroi.shannon))

cor.test(bacteroi.mycel$prok_Shannon, log(bacteroi.mycel$mycelium+0.001))

# figure
mel.divers <- data.frame(taxon = c(rep("Prokaryotes", 2), rep("Actinobacteria", 2), rep("Proteobacteria", 2), rep("Acidobacteria", 2), rep("Bacteroidetes", 2)),
                         index=rep(c("richness", "Shannon"), 5),
                         mean=rep("NA", 10), ci.lower=rep("NA", 10), ci.upper=rep("NA", 10))

mel.divers$mean[1] <- cor.test(divers.mycel$prok_rich, log(divers.mycel$mycelium+0.001))$estimate
mel.divers$ci.lower[1] <- cor.test(divers.mycel$prok_rich, log(divers.mycel$mycelium+0.001))$conf.int[1]
mel.divers$ci.upper[1] <- cor.test(divers.mycel$prok_rich, log(divers.mycel$mycelium+0.001))$conf.int[2]

mel.divers$mean[2] <- cor.test(divers.mycel$prok_Shannon, log(divers.mycel$mycelium+0.001))$estimate
mel.divers$ci.lower[2] <- cor.test(divers.mycel$prok_Shannon, log(divers.mycel$mycelium+0.001))$conf.int[1]
mel.divers$ci.upper[2] <- cor.test(divers.mycel$prok_Shannon, log(divers.mycel$mycelium+0.001))$conf.int[2]

mel.divers$mean[3] <- cor.test(sqrt(actino.mycel$prok_rich+250), log(actino.mycel$mycelium+0.001))$estimate
mel.divers$ci.lower[3] <- cor.test(sqrt(actino.mycel$prok_rich+250), log(actino.mycel$mycelium+0.001))$conf.int[1]
mel.divers$ci.upper[3] <- cor.test(sqrt(actino.mycel$prok_rich+250), log(actino.mycel$mycelium+0.001))$conf.int[2]

mel.divers$mean[4] <- cor.test(actino.mycel$prok_Shannon, log(actino.mycel$mycelium+0.001))$estimate
mel.divers$ci.lower[4] <- cor.test(actino.mycel$prok_Shannon, log(actino.mycel$mycelium+0.001))$conf.int[1]
mel.divers$ci.upper[4] <- cor.test(actino.mycel$prok_Shannon, log(actino.mycel$mycelium+0.001))$conf.int[2]

mel.divers$mean[5] <- cor.test(sqrt(proteo.mycel$prok_rich+350), log(proteo.mycel$mycelium+0.001))$estimate
mel.divers$ci.lower[5] <- cor.test(sqrt(proteo.mycel$prok_rich+350), log(proteo.mycel$mycelium+0.001))$conf.int[1]
mel.divers$ci.upper[5] <- cor.test(sqrt(proteo.mycel$prok_rich+350), log(proteo.mycel$mycelium+0.001))$conf.int[2]

mel.divers$mean[6] <- cor.test(proteo.mycel$prok_Shannon, log(proteo.mycel$mycelium+0.001))$estimate
mel.divers$ci.lower[6] <- cor.test(proteo.mycel$prok_Shannon, log(proteo.mycel$mycelium+0.001))$conf.int[1]
mel.divers$ci.upper[6] <- cor.test(proteo.mycel$prok_Shannon, log(proteo.mycel$mycelium+0.001))$conf.int[2]

mel.divers$mean[7] <- cor.test(acido.mycel$prok_rich, log(acido.mycel$mycelium+0.001))$estimate
mel.divers$ci.lower[7] <- cor.test(acido.mycel$prok_rich, log(acido.mycel$mycelium+0.001))$conf.int[1]
mel.divers$ci.upper[7] <- cor.test(acido.mycel$prok_rich, log(acido.mycel$mycelium+0.001))$conf.int[2]

mel.divers$mean[8] <- cor.test(acido.mycel$prok_Shannon, log(acido.mycel$mycelium+0.001))$estimate
mel.divers$ci.lower[8] <- cor.test(acido.mycel$prok_Shannon, log(acido.mycel$mycelium+0.001))$conf.int[1]
mel.divers$ci.upper[8] <- cor.test(acido.mycel$prok_Shannon, log(acido.mycel$mycelium+0.001))$conf.int[2]

mel.divers$mean[9] <- cor.test(bacteroi.mycel$prok_rich, log(bacteroi.mycel$mycelium+0.001))$estimate
mel.divers$ci.lower[9] <- cor.test(bacteroi.mycel$prok_rich, log(bacteroi.mycel$mycelium+0.001))$conf.int[1]
mel.divers$ci.upper[9] <- cor.test(bacteroi.mycel$prok_rich, log(bacteroi.mycel$mycelium+0.001))$conf.int[2]

mel.divers$mean[10] <- cor.test(bacteroi.mycel$prok_Shannon, log(bacteroi.mycel$mycelium+0.001))$estimate
mel.divers$ci.lower[10] <- cor.test(bacteroi.mycel$prok_Shannon, log(bacteroi.mycel$mycelium+0.001))$conf.int[1]
mel.divers$ci.upper[10] <- cor.test(bacteroi.mycel$prok_Shannon, log(bacteroi.mycel$mycelium+0.001))$conf.int[2]

mel.divers$mean <- as.numeric(mel.divers$mean)
mel.divers$ci.lower <- as.numeric(mel.divers$ci.lower)
mel.divers$ci.upper <- as.numeric(mel.divers$ci.upper)
mel.divers$taxon <- factor(mel.divers$taxon, levels=c("Prokaryotes", "Actinobacteria", "Proteobacteria", "Acidobacteria", "Bacteroidetes"))
mel.divers$index <- factor(mel.divers$index, levels=c("richness", "Shannon"))

ggplot(mel.divers, aes(x=taxon, y=mean, fill=index)) +
  geom_bar(position="dodge", stat="identity", color="black") +
  geom_errorbar(aes(x=taxon, ymin=ci.lower, ymax=ci.upper), width=.3, position=position_dodge(.9)) +
  labs(x="", y="Correl. coeff. with truffle mycelium abundance") +
  scale_fill_discrete(labels=c("Richness", "Shannon index")) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line.y = element_line(colour = "black",linetype = "solid"),
        axis.ticks.x = element_blank()) +
  geom_hline(yintercept=0) +
  ggtitle("(b)")

# other phylla

chlor <- subset_taxa(physeq_prok_uparse_mSeq_perman, Phylum=="Chloroflexi")
gemmati <- subset_taxa(physeq_prok_uparse_mSeq_perman, Phylum=="Gemmatimonadetes")
plancton <- subset_taxa(physeq_prok_uparse_mSeq_perman, Phylum=="Planctomycetes")
thauma <- subset_taxa(physeq_prok_uparse_mSeq_perman, Phylum=="Thaumarchaeota")
verru <- subset_taxa(physeq_prok_uparse_mSeq_perman, Phylum=="Verrucomicrobia")

# Chloroflexi

chlor2 <- transform_sample_counts(chlor, function(x) round(x, digits=0))
divers.chlor <- estimate_richness(chlor2, measure = c("Observed", "Shannon"))
divers.chlor$n_sample <- rownames(divers.chlor)
chlor.mycel <- merge(x = divers.chlor, y = mel.myc, by = "n_sample", all = TRUE)
colnames(chlor.mycel)[2:3] <- c("prok_rich", "prok_Shannon")

lm.chlor.rich <- lm(prok_rich ~ log(mycelium+0.001), chlor.mycel)
summary(lm.chlor.rich)
par(mfcol=c(2,2))
plot(lm.chlor.rich)
shapiro.test(residuals(lm.chlor.rich))

cor.test(chlor.mycel$prok_rich, log(chlor.mycel$mycelium+0.001))

lm.chlor.shannon <- lm(prok_Shannon ~ log(mycelium+0.001), chlor.mycel)
summary(lm.chlor.shannon)
par(mfcol=c(2,2))
plot(lm.chlor.shannon)
shapiro.test(residuals(lm.chlor.shannon))
# no-norm, pero sqrt i log empitjoren

cor.test(chlor.mycel$prok_Shannon, log(chlor.mycel$mycelium+0.001))

# Gemmatimonadetes

gemmati2 <- transform_sample_counts(gemmati, function(x) round(x, digits=0))
divers.gemmati <- estimate_richness(gemmati2, measure = c("Observed", "Shannon"))
divers.gemmati$n_sample <- rownames(divers.gemmati)
gemmati.mycel <- merge(x = divers.gemmati, y = mel.myc, by = "n_sample", all = TRUE)
colnames(gemmati.mycel)[2:3] <- c("prok_rich", "prok_Shannon")

lm.gemmati.rich <- lm(prok_rich ~ log(mycelium+0.001), gemmati.mycel)
summary(lm.gemmati.rich)
par(mfcol=c(2,2))
plot(lm.gemmati.rich)
shapiro.test(residuals(lm.gemmati.rich))
# no-normal, pero sqrt i log empitjoren

cor.test(gemmati.mycel$prok_rich, log(gemmati.mycel$mycelium+0.001))

lm.gemmati.shannon <- lm(prok_Shannon ~ log(mycelium+0.001), gemmati.mycel)
summary(lm.gemmati.shannon)
par(mfcol=c(2,2))
plot(lm.gemmati.shannon)
shapiro.test(residuals(lm.gemmati.shannon))

cor.test(gemmati.mycel$prok_Shannon, log(gemmati.mycel$mycelium+0.001))

# Planctomycetes

plancton2 <- transform_sample_counts(plancton, function(x) round(x, digits=0))
divers.plancton <- estimate_richness(plancton2, measure = c("Observed", "Shannon"))
divers.plancton$n_sample <- rownames(divers.plancton)
plancton.mycel <- merge(x = divers.plancton, y = mel.myc, by = "n_sample", all = TRUE)
colnames(plancton.mycel)[2:3] <- c("prok_rich", "prok_Shannon")

lm.plancton.rich <- lm(sqrt(prok_rich+250) ~ log(mycelium+0.001), plancton.mycel)
summary(lm.plancton.rich)
par(mfcol=c(2,2))
plot(lm.plancton.rich)
shapiro.test(residuals(lm.plancton.rich))

cor.test(plancton.mycel$prok_rich, log(plancton.mycel$mycelium+0.001))

lm.plancton.shannon <- lm(prok_Shannon ~ log(mycelium+0.001), plancton.mycel)
summary(lm.plancton.shannon)
par(mfcol=c(2,2))
plot(lm.plancton.shannon)
shapiro.test(residuals(lm.plancton.shannon))

cor.test(plancton.mycel$prok_Shannon, log(plancton.mycel$mycelium+0.001))

# Thaumarchaeota

thauma2 <- transform_sample_counts(thauma, function(x) round(x, digits=0))
divers.thauma <- estimate_richness(thauma2, measure = c("Observed", "Shannon"))
divers.thauma$n_sample <- rownames(divers.thauma)
thauma.mycel <- merge(x = divers.thauma, y = mel.myc, by = "n_sample", all = TRUE)
colnames(thauma.mycel)[2:3] <- c("prok_rich", "prok_Shannon")

lm.thauma.rich <- lm(prok_rich ~ log(mycelium+0.001), thauma.mycel)
summary(lm.thauma.rich)
par(mfcol=c(2,2))
plot(lm.thauma.rich)
shapiro.test(residuals(lm.thauma.rich))
#no-normal, pero sqrt i log empitjoren

cor.test(thauma.mycel$prok_rich, log(thauma.mycel$mycelium+0.001))

lm.thauma.shannon <- lm(prok_Shannon ~ log(mycelium+0.001), thauma.mycel)
summary(lm.thauma.shannon)
par(mfcol=c(2,2))
plot(lm.thauma.shannon)
shapiro.test(residuals(lm.thauma.shannon))

cor.test(thauma.mycel$prok_Shannon, log(thauma.mycel$mycelium+0.001))

# Verrucomicrobia

verru2 <- transform_sample_counts(verru, function(x) round(x, digits=0))
divers.verru <- estimate_richness(verru2, measure = c("Observed", "Shannon"))
divers.verru$n_sample <- rownames(divers.verru)
verru.mycel <- merge(x = divers.verru, y = mel.myc, by = "n_sample", all = TRUE)
colnames(verru.mycel)[2:3] <- c("prok_rich", "prok_Shannon")

lm.verru.rich <- lm(prok_rich ~ log(mycelium+0.001), verru.mycel)
summary(lm.verru.rich)
par(mfcol=c(2,2))
plot(lm.verru.rich)
shapiro.test(residuals(lm.verru.rich))

cor.test(verru.mycel$prok_rich, log(verru.mycel$mycelium+0.001))

lm.verru.shannon <- lm(prok_Shannon ~ log(mycelium+0.001), verru.mycel)
summary(lm.verru.shannon)
par(mfcol=c(2,2))
plot(lm.verru.shannon)
shapiro.test(residuals(lm.verru.shannon))

cor.test(verru.mycel$prok_Shannon, log(verru.mycel$mycelium+0.001))

# figure Verrucomicrobia
verru3 <- subset_taxa(physeq_prok_uparse_mSeq_PCoA, Phylum=="Verrucomicrobia")

ps1.verru.rel <- microbiome::transform(verru3, "compositional")

ps1.verru.gen.rel <-aggregate_rare(ps1.verru.rel, level = "Genus", detection = 0.01, prevalence = 0.01)
dat.trans = transform_sample_counts(ps1.verru.gen.rel, function(x) 100*x/sum(x))
dat.dataframe = psmelt(dat.trans)
dat.agr = aggregate(Abundance~etiq+Genus, data=dat.dataframe, FUN=mean)
dat.agr2 <- subset(dat.agr, etiq%in%c("P0M", "P0T", "P1M", "P1T", "P2M", "P2T", "S2M", "S2T"))
dat.agr3 <- dat.agr2[!(dat.agr2$Genus == "Other"),]
dat.agr3 <- dat.agr3[!(dat.agr3$Genus == "Unknown"),]
aggregate(dat.agr3$Abundance, list(dat.agr3$Genus), FUN=mean)

ggplot(dat.agr3, aes(x=etiq, y=Abundance, fill=Genus)) +
  geom_bar(stat="identity") +
  ylab("Relative abundance within Verrucomicrobia (%)") +
  scale_x_discrete(labels=c("N0M", "N0T", "N1M", "N1T", "N2M", "N2T", "S2M", "S2T")) +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black", size = 1, linetype = "solid"))+
  guides(fill=guide_legend(ncol =1)) +
  scale_fill_brewer(palette = "Paired")

# proteobacteria
proteo2 <- subset_taxa(physeq_prok_uparse_mSeq_PCoA, Phylum=="Proteobacteria")
ps1.proteo.rel <- microbiome::transform(proteo2, "compositional")
ps1.proteo.rel.class <-aggregate_rare(ps1.proteo.rel, level = "Class", detection = 0.005, prevalence = 0.1)
p <- plot_composition(ps1.proteo.rel.class,
                      average_by = "etiq") + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = "", 
       y = "Relative abundance",
       title = "", 
       subtitle = "",
       caption = "") 
print(p + scale_fill_brewer("Class", palette = "Paired") + theme_bw())

# analysis per class
physeq_prok_uparse_mSeq_perman_percent_clas <- tax_glom(physeq_prok_uparse_mSeq_perman_percent, taxrank="Class")
df.class.merged <- data.frame(t(otu_table(physeq_prok_uparse_mSeq_perman_percent_clas)))
df.class.merged$n_sample <- rownames(df.class.merged)
df.class.merged <- merge(x = df.class.merged, y = mel.myc, by = "n_sample", all = TRUE)
# tax_table(physeq_prok_uparse_mSeq_perman_percent_clas)
lm.clas.betaproteo <- lm(log(BOTU_4) ~ log(mycelium+0.001), df.class.merged)
summary(lm.clas.betaproteo)
par(mfcol=c(2,2))
plot(lm.clas.betaproteo)
shapiro.test(residuals(lm.clas.betaproteo))

cor.test(log(df.class.merged$BOTU_4), log(df.class.merged$mycelium+0.001))

lm.clas.alphaproteo <- lm(BOTU_6 ~ log(mycelium+0.001), df.class.merged)
summary(lm.clas.alphaproteo)
par(mfcol=c(2,2))
plot(lm.clas.alphaproteo)
shapiro.test(residuals(lm.clas.alphaproteo))

cor.test(df.class.merged$BOTU_6, log(df.class.merged$mycelium+0.001))

lm.clas.gammaproteo <- lm(log(BOTU_7) ~ log(mycelium+0.001), df.class.merged)
summary(lm.clas.gammaproteo)
par(mfcol=c(2,2))
plot(lm.clas.gammaproteo)
shapiro.test(residuals(lm.clas.gammaproteo))

cor.test(log(df.class.merged$BOTU_7), log(df.class.merged$mycelium+0.001))
# no normalitzable, pero P=0.51

lm.clas.deltaproteo <- lm(sqrt(BOTU_31+0.005) ~ log(mycelium+0.001), df.class.merged)
summary(lm.clas.deltaproteo)
par(mfcol=c(2,2))
plot(lm.clas.deltaproteo)
shapiro.test(residuals(lm.clas.deltaproteo))

cor.test(sqrt(df.class.merged$BOTU_31+0.005), log(df.class.merged$mycelium+0.001))

# same thing with orders...
ps1.proteo.rel.order <-aggregate_rare(ps1.proteo.rel, level = "Order", detection = 0.01, prevalence = 0.3)
p <- plot_composition(ps1.proteo.rel.order,
                      average_by = "etiq") + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = "", 
       y = "Relative abundance",
       title = "", 
       subtitle = "",
       caption = "") 
print(p + scale_fill_brewer("Class", palette = "Paired") + theme_bw())

dat.trans = transform_sample_counts(ps1.proteo.rel.order, function(x) 100*x/sum(x))
dat.proteo.ord = psmelt(dat.trans)
dat.agr = aggregate(Abundance~etiq+Order, data=dat.proteo.ord, FUN=mean)
dat.agr2 <- subset(dat.agr, etiq%in%c("P0M", "P0T", "P1M", "P1T", "P2M", "P2T", "S2M", "S2T"))
aggregate(dat.agr2$Abundance, list(dat.agr2$Order), FUN=mean)

# analitzar principals ordres: Rhizobiales BOTU_6, Sphingomonadales BOTU_76, Xanthomonadales BOTU_7, Rhodospirillales BOTU_224, Burkholderiales BOTU_92, Myxococcales BOTU_31, Caulobacterales BOTU_56, Pseudomonadales BOTU_48
# Rhizobiales: en abundancia relativa NO hay dif sign, en richness SI
physeq_prok_uparse_mSeq_perman_percent_order <- tax_glom(physeq_prok_uparse_mSeq_perman_percent, taxrank="Order")
df.order.merged <- data.frame(t(otu_table(physeq_prok_uparse_mSeq_perman_percent_order)))
df.order.merged$n_sample <- rownames(df.order.merged)
df.order.merged <- merge(x = df.order.merged, y = mel.myc, by = "n_sample", all = TRUE)
# tax_table(physeq_prok_uparse_mSeq_perman_percent_order)
lm.ord.rhiz <- lm(BOTU_6 ~ log(mycelium+0.001), df.order.merged)
summary(lm.ord.rhiz)
par(mfcol=c(2,2))
plot(lm.ord.rhiz)
shapiro.test(residuals(lm.ord.rhiz))

cor.test(df.order.merged$BOTU_6, log(df.order.merged$mycelium+0.001))

rhiz <- subset_taxa(physeq_prok_uparse_mSeq_perman, Order=="Rhizobiales")

rhiz2 <- transform_sample_counts(rhiz, function(x) round(x, digits=0))
divers.rhiz <- estimate_richness(rhiz2, measure = c("Observed", "Shannon"))
divers.rhiz$n_sample <- rownames(divers.rhiz)
rhiz.mycel <- merge(x = divers.rhiz, y = mel.myc, by = "n_sample", all = TRUE)
colnames(rhiz.mycel)[2:3] <- c("prok_rich", "prok_Shannon")

lm.rhiz.rich <- lm(prok_rich ~ log(mycelium+0.001), rhiz.mycel)
summary(lm.rhiz.rich)
par(mfcol=c(2,2))
plot(lm.rhiz.rich)
shapiro.test(residuals(lm.rhiz.rich))
# tranf empitjoren normalitat, deixe var. original

cor.test(rhiz.mycel$prok_rich, log(rhiz.mycel$mycelium+0.001))

lm.rhiz.shannon <- lm(prok_Shannon ~ log(mycelium+0.001), rhiz.mycel)
summary(lm.rhiz.shannon)
par(mfcol=c(2,2))
plot(lm.rhiz.shannon)
shapiro.test(residuals(lm.rhiz.shannon))

cor.test(rhiz.mycel$prok_Shannon, log(rhiz.mycel$mycelium+0.001))

# Sphingomonadales
sphin <- subset_taxa(physeq_prok_uparse_mSeq_perman, Order=="Sphingomonadales")

sphin2 <- transform_sample_counts(sphin, function(x) round(x, digits=0))
divers.sphin <- estimate_richness(sphin2, measure = c("Observed", "Shannon"))
divers.sphin$n_sample <- rownames(divers.sphin)
sphin.mycel <- merge(x = divers.sphin, y = mel.myc, by = "n_sample", all = TRUE)
colnames(sphin.mycel)[2:3] <- c("prok_rich", "prok_Shannon")

lm.sphin.rich <- lm(prok_rich ~ log(mycelium+0.001), sphin.mycel)
summary(lm.sphin.rich)
par(mfcol=c(2,2))
plot(lm.sphin.rich)
shapiro.test(residuals(lm.sphin.rich))

cor.test(sphin.mycel$prok_rich, log(sphin.mycel$mycelium+0.001))

lm.sphin.shannon <- lm(prok_Shannon ~ log(mycelium+0.001), sphin.mycel)
summary(lm.sphin.shannon)
par(mfcol=c(2,2))
plot(lm.sphin.shannon)
shapiro.test(residuals(lm.sphin.shannon))
# tranf empitjoren normalitat, deixe var. original

cor.test(sphin.mycel$prok_Shannon, log(sphin.mycel$mycelium+0.001))

# Xanthomonadales
xantho <- subset_taxa(physeq_prok_uparse_mSeq_perman, Order=="Xanthomonadales")

xantho2 <- transform_sample_counts(xantho, function(x) round(x, digits=0))
divers.xantho <- estimate_richness(xantho2, measure = c("Observed", "Shannon"))
divers.xantho$n_sample <- rownames(divers.xantho)
xantho.mycel <- merge(x = divers.xantho, y = mel.myc, by = "n_sample", all = TRUE)
colnames(xantho.mycel)[2:3] <- c("prok_rich", "prok_Shannon")

lm.xantho.rich <- lm(prok_rich ~ log(mycelium+0.001), xantho.mycel)
summary(lm.xantho.rich)
par(mfcol=c(2,2))
plot(lm.xantho.rich)
shapiro.test(residuals(lm.xantho.rich))

cor.test(xantho.mycel$prok_rich, log(xantho.mycel$mycelium+0.001))

lm.xantho.shannon <- lm(prok_Shannon ~ log(mycelium+0.001), xantho.mycel)
summary(lm.xantho.shannon)
par(mfcol=c(2,2))
plot(lm.xantho.shannon)
shapiro.test(residuals(lm.xantho.shannon))
# tranf empitjoren normalitat, deixe var. original

cor.test(xantho.mycel$prok_Shannon, log(xantho.mycel$mycelium+0.001))

# Rhodospirillales
rhodo <- subset_taxa(physeq_prok_uparse_mSeq_perman, Order=="Rhodospirillales")

rhodo2 <- transform_sample_counts(rhodo, function(x) round(x, digits=0))
divers.rhodo <- estimate_richness(rhodo2, measure = c("Observed", "Shannon"))
divers.rhodo$n_sample <- rownames(divers.rhodo)
rhodo.mycel <- merge(x = divers.rhodo, y = mel.myc, by = "n_sample", all = TRUE)
colnames(rhodo.mycel)[2:3] <- c("prok_rich", "prok_Shannon")

lm.rhodo.rich <- lm(prok_rich ~ log(mycelium+0.001), rhodo.mycel)
summary(lm.rhodo.rich)
par(mfcol=c(2,2))
plot(lm.rhodo.rich)
shapiro.test(residuals(lm.rhodo.rich))

cor.test(rhodo.mycel$prok_rich, log(rhodo.mycel$mycelium+0.001))

lm.rhodo.shannon <- lm(prok_Shannon ~ log(mycelium+0.001), rhodo.mycel)
summary(lm.rhodo.shannon)
par(mfcol=c(2,2))
plot(lm.rhodo.shannon)
shapiro.test(residuals(lm.rhodo.shannon))
# tranf empitjoren normalitat, deixe var. original

cor.test(rhodo.mycel$prok_Shannon, log(rhodo.mycel$mycelium+0.001))

# Burkholderiales
burk <- subset_taxa(physeq_prok_uparse_mSeq_perman, Order=="Burkholderiales")

burk2 <- transform_sample_counts(burk, function(x) round(x, digits=0))
divers.burk <- estimate_richness(burk2, measure = c("Observed", "Shannon"))
divers.burk$n_sample <- rownames(divers.burk)
burk.mycel <- merge(x = divers.burk, y = mel.myc, by = "n_sample", all = TRUE)
colnames(burk.mycel)[2:3] <- c("prok_rich", "prok_Shannon")

lm.burk.rich <- lm(prok_rich ~ log(mycelium+0.001), burk.mycel)
summary(lm.burk.rich)
par(mfcol=c(2,2))
plot(lm.burk.rich)
shapiro.test(residuals(lm.burk.rich))
# transf empitjora, deixe vaar. original

cor.test(burk.mycel$prok_rich, log(burk.mycel$mycelium+0.001))

lm.burk.shannon <- lm(prok_Shannon ~ log(mycelium+0.001), burk.mycel)
summary(lm.burk.shannon)
par(mfcol=c(2,2))
plot(lm.burk.shannon)
shapiro.test(residuals(lm.burk.shannon))
# tranf empitjoren normalitat, deixe var. original

cor.test(burk.mycel$prok_Shannon, log(burk.mycel$mycelium+0.001))

# Myxococcales
myxo <- subset_taxa(physeq_prok_uparse_mSeq_perman, Order=="Myxococcales")

myxo2 <- transform_sample_counts(myxo, function(x) round(x, digits=0))
divers.myxo <- estimate_richness(myxo2, measure = c("Observed", "Shannon"))
divers.myxo$n_sample <- rownames(divers.myxo)
myxo.mycel <- merge(x = divers.myxo, y = mel.myc, by = "n_sample", all = TRUE)
colnames(myxo.mycel)[2:3] <- c("prok_rich", "prok_Shannon")

lm.myxo.rich <- lm(prok_rich ~ log(mycelium+0.001), myxo.mycel)
summary(lm.myxo.rich)
par(mfcol=c(2,2))
plot(lm.myxo.rich)
shapiro.test(residuals(lm.myxo.rich))

cor.test(myxo.mycel$prok_rich, log(myxo.mycel$mycelium+0.001))

lm.myxo.shannon <- lm(prok_Shannon ~ log(mycelium+0.001), myxo.mycel)
summary(lm.myxo.shannon)
par(mfcol=c(2,2))
plot(lm.myxo.shannon)
shapiro.test(residuals(lm.myxo.shannon))

cor.test(myxo.mycel$prok_Shannon, log(myxo.mycel$mycelium+0.001))

# Caulobacterales
caulo <- subset_taxa(physeq_prok_uparse_mSeq_perman, Order=="Caulobacterales")

caulo2 <- transform_sample_counts(caulo, function(x) round(x, digits=0))
divers.caulo <- estimate_richness(caulo2, measure = c("Observed", "Shannon"))
divers.caulo$n_sample <- rownames(divers.caulo)
caulo.mycel <- merge(x = divers.caulo, y = mel.myc, by = "n_sample", all = TRUE)
colnames(caulo.mycel)[2:3] <- c("prok_rich", "prok_Shannon")

lm.caulo.rich <- lm(prok_rich ~ log(mycelium+0.001), caulo.mycel)
summary(lm.caulo.rich)
par(mfcol=c(2,2))
plot(lm.caulo.rich)
shapiro.test(residuals(lm.caulo.rich))
# transf empitjora, deixe vaar. original

cor.test(caulo.mycel$prok_rich, log(caulo.mycel$mycelium+0.001))

lm.caulo.shannon <- lm(prok_Shannon ~ log(mycelium+0.001), caulo.mycel)
summary(lm.caulo.shannon)
par(mfcol=c(2,2))
plot(lm.caulo.shannon)
shapiro.test(residuals(lm.caulo.shannon))
# tranf empitjoren normalitat, deixe var. original

cor.test(caulo.mycel$prok_Shannon, log(caulo.mycel$mycelium+0.001))

# Pseudomonadales
pseudo <- subset_taxa(physeq_prok_uparse_mSeq_perman, Order=="Pseudomonadales")

pseudo2 <- transform_sample_counts(pseudo, function(x) round(x, digits=0))
divers.pseudo <- estimate_richness(pseudo2, measure = c("Observed", "Shannon"))
divers.pseudo$n_sample <- rownames(divers.pseudo)
pseudo.mycel <- merge(x = divers.pseudo, y = mel.myc, by = "n_sample", all = TRUE)
colnames(pseudo.mycel)[2:3] <- c("prok_rich", "prok_Shannon")

lm.pseudo.rich <- lm(log(prok_rich) ~ log(mycelium+0.001), pseudo.mycel)
summary(lm.pseudo.rich)
par(mfcol=c(2,2))
plot(lm.pseudo.rich)
shapiro.test(residuals(lm.pseudo.rich))

cor.test(log(pseudo.mycel$prok_rich), log(pseudo.mycel$mycelium+0.001))

lm.pseudo.shannon <- lm(log(prok_Shannon+0.1) ~ log(mycelium+0.001), pseudo.mycel)
summary(lm.pseudo.shannon)
par(mfcol=c(2,2))
plot(lm.pseudo.shannon)
shapiro.test(residuals(lm.pseudo.shannon))
# tranf empitjoren normalitat, deixe var. original

cor.test(log(pseudo.mycel$prok_Shannon + 0.1), log(pseudo.mycel$mycelium+0.001))

# figures per proteobacteria orders
rhiz3 <- subset_taxa(physeq_prok_uparse_mSeq_PCoA, Order=="Rhizobiales")

ps1.rhiz.rel <- microbiome::transform(rhiz3, "compositional")

ps1.rhiz.gen.rel <-aggregate_rare(ps1.rhiz.rel, level = "Genus", detection = 0.01, prevalence = 0.01)
dat.trans = transform_sample_counts(ps1.rhiz.gen.rel, function(x) 100*x/sum(x))
dat.dataframe = psmelt(dat.trans)
dat.agr = aggregate(Abundance~etiq+Genus, data=dat.dataframe, FUN=mean)
dat.agr2 <- subset(dat.agr, etiq%in%c("P0M", "P0T", "P1M", "P1T", "P2M", "P2T", "S2M", "S2T"))
dat.agr3 <- dat.agr2[!(dat.agr2$Genus == "Other"),]
dat.agr3 <- dat.agr3[!(dat.agr3$Genus == "Unknown"),]
aggregate(dat.agr3$Abundance, list(dat.agr3$Genus), FUN=mean)

getPalette = colorRampPalette(brewer.pal(18, "Paired"))
ggplot(dat.agr3, aes(x=etiq, y=Abundance, fill=Genus)) +
  geom_bar(stat="identity") +
  ylab("Relative abundance within Rhizobiales (%)") +
  scale_x_discrete(labels=c("N0M", "N0T", "N1M", "N1T", "N2M", "N2T", "S2M", "S2T")) +
  ggtitle("(a)") +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black", size = 1, linetype = "solid"))+
  guides(fill=guide_legend(ncol =1)) +
  scale_fill_manual(values=getPalette(18))

xantho3 <- subset_taxa(physeq_prok_uparse_mSeq_PCoA, Order=="Xanthomonadales")
ps1.xantho.rel <- microbiome::transform(xantho3, "compositional")
ps1.xantho.gen.rel <-aggregate_rare(ps1.xantho.rel, level = "Genus", detection = 0.01, prevalence = 0.01)
dat.trans = transform_sample_counts(ps1.xantho.gen.rel, function(x) 100*x/sum(x))
dat.dataframe = psmelt(dat.trans)
dat.agr = aggregate(Abundance~etiq+Genus, data=dat.dataframe, FUN=mean)
dat.agr2 <- subset(dat.agr, etiq%in%c("P0M", "P0T", "P1M", "P1T", "P2M", "P2T", "S2M", "S2T"))
dat.agr3 <- dat.agr2[!(dat.agr2$Genus == "Other"),]
dat.agr3 <- dat.agr3[!(dat.agr3$Genus == "Unknown"),]
aggregate(dat.agr3$Abundance, list(dat.agr3$Genus), FUN=mean)

getPalette = colorRampPalette(brewer.pal(14, "Paired"))
ggplot(dat.agr3, aes(x=etiq, y=Abundance, fill=Genus)) +
  geom_bar(stat="identity") +
  ylab("Relative abundance within Xanthomonadales (%)") +
  scale_x_discrete(labels=c("N0M", "N0T", "N1M", "N1T", "N2M", "N2T", "S2M", "S2T")) +
  ggtitle("(b)") +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black", size = 1, linetype = "solid"))+
  guides(fill=guide_legend(ncol =1)) +
  scale_fill_manual(values=getPalette(14))

rhodo3 <- subset_taxa(physeq_prok_uparse_mSeq_PCoA, Order=="Rhodospirillales")
ps1.rhodo.rel <- microbiome::transform(rhodo3, "compositional")
ps1.rhodo.gen.rel <-aggregate_rare(ps1.rhodo.rel, level = "Genus", detection = 0.01, prevalence = 0.01)
dat.trans = transform_sample_counts(ps1.rhodo.gen.rel, function(x) 100*x/sum(x))
dat.dataframe = psmelt(dat.trans)
dat.agr = aggregate(Abundance~etiq+Genus, data=dat.dataframe, FUN=mean)
dat.agr2 <- subset(dat.agr, etiq%in%c("P0M", "P0T", "P1M", "P1T", "P2M", "P2T", "S2M", "S2T"))
dat.agr3 <- dat.agr2[!(dat.agr2$Genus == "Other"),]
dat.agr3 <- dat.agr3[!(dat.agr3$Genus == "Unknown"),]
aggregate(dat.agr3$Abundance, list(dat.agr3$Genus), FUN=mean)

ggplot(dat.agr3, aes(x=etiq, y=Abundance, fill=Genus)) +
  geom_bar(stat="identity") +
  ylab("Relative abundance within Rhodospirillales (%)") +
  scale_x_discrete(labels=c("N0M", "N0T", "N1M", "N1T", "N2M", "N2T", "S2M", "S2T")) +
  ggtitle("(c)") +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black", size = 1, linetype = "solid"))+
  guides(fill=guide_legend(ncol =1)) +
  scale_fill_brewer(palette = "Paired")

burk3 <- subset_taxa(physeq_prok_uparse_mSeq_PCoA, Order=="Burkholderiales")
ps1.burk.rel <- microbiome::transform(burk3, "compositional")
ps1.burk.gen.rel <-aggregate_rare(ps1.burk.rel, level = "Genus", detection = 0.01, prevalence = 0.01)
dat.trans = transform_sample_counts(ps1.burk.gen.rel, function(x) 100*x/sum(x))
dat.dataframe = psmelt(dat.trans)
dat.agr = aggregate(Abundance~etiq+Genus, data=dat.dataframe, FUN=mean)
dat.agr2 <- subset(dat.agr, etiq%in%c("P0M", "P0T", "P1M", "P1T", "P2M", "P2T", "S2M", "S2T"))
dat.agr3 <- dat.agr2[!(dat.agr2$Genus == "Other"),]
dat.agr3 <- dat.agr3[!(dat.agr3$Genus == "Unknown"),]
aggregate(dat.agr3$Abundance, list(dat.agr3$Genus), FUN=mean)

getPalette = colorRampPalette(brewer.pal(14, "Paired"))
ggplot(dat.agr3, aes(x=etiq, y=Abundance, fill=Genus)) +
  geom_bar(stat="identity") +
  ylab("Relative abundance within Burkholderiales (%)") +
  scale_x_discrete(labels=c("N0M", "N0T", "N1M", "N1T", "N2M", "N2T", "S2M", "S2T")) +
  ggtitle("(d)") +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black", size = 1, linetype = "solid"))+
  guides(fill=guide_legend(ncol =1)) +
  scale_fill_manual(values=getPalette(14))

myxo3 <- subset_taxa(physeq_prok_uparse_mSeq_PCoA, Order=="Myxococcales")
ps1.myxo.rel <- microbiome::transform(myxo3, "compositional")
ps1.myxo.gen.rel <-aggregate_rare(ps1.myxo.rel, level = "Genus", detection = 0.01, prevalence = 0.01)
dat.trans = transform_sample_counts(ps1.myxo.gen.rel, function(x) 100*x/sum(x))
dat.dataframe = psmelt(dat.trans)
dat.agr = aggregate(Abundance~etiq+Genus, data=dat.dataframe, FUN=mean)
dat.agr2 <- subset(dat.agr, etiq%in%c("P0M", "P0T", "P1M", "P1T", "P2M", "P2T", "S2M", "S2T"))
dat.agr3 <- dat.agr2[!(dat.agr2$Genus == "Other"),]
dat.agr3 <- dat.agr3[!(dat.agr3$Genus == "Unknown"),]
aggregate(dat.agr3$Abundance, list(dat.agr3$Genus), FUN=mean)

getPalette = colorRampPalette(brewer.pal(17, "Paired"))
ggplot(dat.agr3, aes(x=etiq, y=Abundance, fill=Genus)) +
  geom_bar(stat="identity") +
  ylab("Relative abundance within Myxococcales (%)") +
  scale_x_discrete(labels=c("N0M", "N0T", "N1M", "N1T", "N2M", "N2T", "S2M", "S2T")) +
  ggtitle("(e)") +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black", size = 1, linetype = "solid"))+
  guides(fill=guide_legend(ncol =1)) +
  scale_fill_manual(values=getPalette(17))

# CANONICAL CORRESPONDENCE ANALYSIS
otu_prok_p <- otu_prok
for (i in 1:ncol(otu_prok)){
  otu_prok_p[,i] <- otu_prok[,i]/sum(otu_prok[,i])
}
otu_prok_t <- as.data.frame(t(otu_prok_p))
otu_prok_t <- otu_prok_t[, colSums(otu_prok_t) != 0]

set.seed(3087)
cca_veg <- cca(otu_prok_t ~ log(mycelium+0.001), metadata_prok, scale = T)
cca_veg

anova.cca(cca_veg, by="terms", permutations=2000)
anova.cca(cca_veg, by="axis")

library(ggvegan)
# OTUs linked to truffle mycelium
summar <- summary(cca_veg)
summar.sp <- as.data.frame(summar$species[summar$species[,1] < -0.878 , ]) # en valor absoluto, superior al valor maximo positivo de una OTU en el eje CCA1 (sin valor estadistico)
summar.sp <- summar.sp[, -c(3:6)]
temptax <- as.data.frame(tax_table(physeq_prok_uparse_mSeq))
temptax <- subset(temptax, row.names(temptax) %in% colnames(otu_prok_t))
temptax$OTU <- rownames(temptax)
summar.sp$OTU <- rownames(summar.sp)
summar.sp_sign <- merge(summar.sp, temptax, by="OTU")
otu_prok_sign <- subset(otu_prok, row.names(otu_prok) %in% summar.sp_sign$OTU)
otu_prok_sign$OTU <- rownames(otu_prok_sign)
summar.sp_sign <- merge(summar.sp_sign, otu_prok_sign, by="OTU")
summar.sp_sign$freq <- rowMeans(summar.sp_sign[,10:54]) # atencion, no es freq rel, es promedio de No reads
summar.sp_sign_domin <- summar.sp_sign[rowSums(summar.sp_sign[,10:54] > 0) > 23, ]
# with presence > 75%, only 3 OTUs, with 50%, 21 OTUs

cca.sites <-data.frame(summar$sites[,1:2])
round(cor(otu_prok_t, cca.sites[,1:2]), 3)

num <- 45
limit.sign_r <- tanh(3.29/ sqrt(num - 3))
limit.sign_r
# alpha=0.001

otu_prok_t_inv <- as.data.frame(t(otu_prok_t))
otu_prok_t_inv$rel.abu <- rowSums(otu_prok_t_inv)/45
otu_prok_t_inv$rel.freq <- rowSums(otu_prok_t_inv[1:45]>0)

FSr  = transform_sample_counts(physeq_prok_uparse_mSeq, function(x) x / sum(x) )
FSfr = filter_taxa(FSr, function(x) sum(x) > .005, TRUE)

# customize autoplot
temptax2 <- temptax
temptax2$Phylum[temptax2$Phylum == "candidate division WPS-1"] <- "Other"
temptax2$Phylum[temptax2$Phylum == "Firmicutes"] <- "Other"
temptax2$Phylum[temptax2$Phylum == "Microgenomates"] <- "Other"
temptax2$Phylum[temptax2$Phylum == "Candidatus Saccharibacteria"] <- "Other"
temptax2$Phylum[temptax2$Phylum == "Nitrospirae"] <- "Other"
temptax2$Phylum[temptax2$Phylum == "Armatimonadetes"] <- "Other"
temptax2$Phylum[temptax2$Phylum == "Chlamydiae"] <- "Other"
temptax2$Phylum[temptax2$Phylum == "Hydrogenedentes"] <- "Other"
temptax2$Phylum[temptax2$Phylum == "Parcubacteria"] <- "Other"
temptax2$Phylum[temptax2$Phylum == "Ignavibacteriae"] <- "Other"
temptax2$Phylum[temptax2$Phylum == "BRC1"] <- "Other"
temptax2$Phylum[temptax2$Phylum == "Euryarchaeota"] <- "Other"
temptax2$Phylum[temptax2$Phylum == "Latescibacteria"] <- "Other"
temptax2$Phylum[temptax2$Phylum == "Aminicenantes"] <- "Other"
temptax2$Phylum[temptax2$Phylum == "candidate division WPS-2"] <- "Other"
temptax2$Phylum[temptax2$Phylum == "Crenarchaeota"] <- "Other"
temptax2$Phylum[temptax2$Phylum == "Deinococcus-Thermus"] <- "Other"
temptax2$Phylum[temptax2$Phylum == "Pacearchaeota"] <- "Other"
temptax2$Phylum[temptax2$Phylum == "Spirochaetes"] <- "Other"
temptax2$Phylum[temptax2$Phylum == "SR1"] <- "Other"
temptax2$Phylum[temptax2$Phylum == "Woesearchaeota"] <- "Other"
temptax2$Phylum <- factor(temptax2$Phylum, levels = c("Acidobacteria", "Actinobacteria", "Bacteroidetes", "Chloroflexi", "Gemmatimonadetes", "Planctomycetes", "Proteobacteria", "Thaumarchaeota", "Verrucomicrobia", "Other"))

fmod <- fortify(cca_veg)
size <- 1.8

cca.species <- data.frame(summar$species)
cca.species2 <- data.frame(CCA1=cca.species$CCA1,CA1=cca.species$CA1)
rownames(cca.species2) <- rownames(cca.species)
fmod[8642, 2] <- "T. melanosporum"

# Fig 5b
ggplot(fmod, aes(x = CCA1, y = CA1)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_point(data = subset(fmod, score == "species"), aes(colour = temptax2$Phylum), size = size) +
  scale_colour_brewer("score", palette = "Set3") +
  coord_fixed() +
  xlim(-4, 2) +
  theme_classic() +
  theme(legend.position = "bottom", legend.title=element_blank()) +
  geom_text(data = subset(fmod, score == "biplot"), aes(x=CCA1,y=CA1,label=label), size=3,
            hjust = 2.5, vjust = -1.5) +
  xlab("CCA1 (4.0%)") + ylab("CA1 (23.0%)") +
  geom_segment(data = subset(fmod, score == "biplot"),
               aes(x = 0, y = 0, xend = CCA1*3.3, yend = CA1), arrow = arrow(length = unit(1/2, 'picas')), colour = "red", linewidth=1.2) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  ggtitle("(c)")

datamyc <- sample_data(physeq_prok_uparse_mSeq_PCoAsinP3_noround)
0.0001 -> datamyc[3,15]
# Fig 5a
ggplot(data = datamyc, aes(x = etiq, y = mycelium)) +
  geom_jitter(position="identity", aes(color=origin)) +
  geom_boxplot(fill=NA, outlier.color=NA) +
  ylab(expression(italic(T.)~italic(melanosporum)~"mycelium abundance (mg g"^-1~"sample)")) +
  scale_color_manual(values=c("grey50", "red")) +
  scale_x_discrete(labels=c("N0M", "N0T", "N1M", "N1T", "N2M", "N2T", "S2M", "S2T")) +
  scale_y_continuous(trans="log10", breaks=c(0.001, 0.01, 0.1, 1, 10)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none", axis.title.x=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black", linewidth=0.5, linetype = "solid")) +
  ggtitle("(a)")

lm.myc <- aov(log(mycelium+0.001) ~ location*year, divers.perman2)

# NETWORK ANALYSIS
otu_prok_count_t <- as.data.frame(t(otu_prok))
otu_prok_count_t <- otu_prok_count_t[, colSums(otu_prok_count_t) != 0]
library(NetCoMi)
netcomi_net <- netConstruct(otu_prok_count_t, dataType = "counts", # otus in columns, samples in rows
                            filtTax = c("numbSamp", "relFreq"), # only taxa present in 75% of samples with a rel freq > 0.05%
                            filtTaxPar = list(numbSamp = 34, relFreq = 0.001),
                            # normMethod = "clr",
                            measure = "spieceasi",
                            measurePar = list(method = "mb", sel.criterion = "stars",
                                              pulsar.params = list(rep.num = 20)), # opciones por defecto, solo explicitadas
                            sparsMethod = "none", # sparsification and transformation are done internally by the function
                            dissFunc = "signed", # makes more biological sense for microbial commun
                            seed = 13075)
# netConstruct returns an object of the class microNet
edgelist <- netcomi_net$edgelist1[order(netcomi_net$edgelist1$adja, 
                                        decreasing = TRUE), ]
head(edgelist)
# convert the adjacency matrix into an igraph object
netcomi_graph <- SpiecEasi::adj2igraph(abs(netcomi_net$adjaMat1))

for (i in 1:ncol(otu_prok_count_t)){
  otu_prok_count_t[,i] <- log(otu_prok_count_t[,i]+1)
}
vsize <- (colMeans(otu_prok_count_t))*2 # node sizes

# Fruchterman-Reingold layout from igraph package
library(igraph)
set.seed(13075)
lay_fr <- layout_with_fr(netcomi_graph)
plot(netcomi_graph, layout = lay_fr, vertex.size = vsize, 
     vertex.label = NA, main = "NetCoMi network\n(with Spiec.Easi associations)")

# centrality measures
get_centr <- function(graph_obj) {
  # We access igraph directly with "::" because there are more packages loaded in 
  # this chapter that contain a degree() function.
  df <- data.frame(Degree = igraph::degree(graph_obj))
  df$Betweenness <- betweenness(graph_obj)
  df$Closeness <- closeness(graph_obj, normalized = TRUE)
  df$Eigenvector <- eigen_centrality(graph_obj)$vector
  return(df)
}

centr_df <- get_centr(netcomi_graph)
rownames(centr_df) <- rownames(netcomi_net$assoMat1)
head(centr_df, 15)

# visualize centrality measures
get_vsizes <- function(centr_df) {
  df <- as.matrix(centr_df)
  df[, "Degree"] <- df[, "Degree"]
  df[, "Betweenness"] <- log(df[, "Betweenness"])
  df[, "Closeness"] <- df[, "Closeness"] * 5
  df[, "Eigenvector"] <- df[, "Eigenvector"] * 10
  df[is.infinite(df) | is.na(df)] <- 0
  return(df)
}

vsize_df <- get_vsizes(centr_df )
head(vsize_df)

par(mfrow = c(2,2))
for (i in seq_along(centr_df)) {
  plot(netcomi_graph, layout = lay_fr, vertex.size = vsize_df[, i], 
       vertex.label = NA, main = colnames(centr_df)[i])
}

# Degree distribution
ddist<- igraph::degree_distribution(netcomi_graph)
df <- data.frame(Degree = as.factor((seq_along(ddist)) - 1),
                 Fraction = ddist)
ggplot(data = df, aes(x = Degree, y = Fraction, group = 1)) +
  geom_line() +
  geom_point() +
  theme_bw()
# The network has no singletons and just a few sparsely connected nodes, and most nodes with a degree of 4-10

# Global network measures: density, transitivity, and average path length
edge_density(netcomi_graph)
transitivity(netcomi_graph)
mean_distance(netcomi_graph)

par(mfrow = c(1,1))
netcomi_netprops <- netAnalyze(netcomi_net, 
                               clustMethod = "cluster_fast_greedy",
                               hubPar = "eigenvector",
                               normDeg = FALSE)
summary(netcomi_netprops, numbNodes = 5)

#plot with NetCoMi (not igraph), OTUs coloured according to cluster
renamed <- read.table(file="temptax2b.txt", header=T, sep="\t")
new.names <- renamed$Renamed
names(new.names) <- renamed$Original

plot(netcomi_netprops,
     repulsion = 0.98,
     rmSingles = "all",
     #shortenLabels = "intelligent", here with names OTUxxx this is not needed
     labels=new.names, labelScale = F,
     cexLabels = 0, cexHubLabels = 0.8, hubBorderWidth=1,
     nodeSize = "eigenvector",
     nodeSizeSpread = 3,
     nodeColor = "cluster", 
     hubBorderCol = "gray40",
     cexNodes = 1.2,
     edgeTranspHigh = 20,
     #title1 = "(a)", 
     #showTitle = TRUE,
     #cexTitle = 2,
     mar = c(1, 1, 4.5, 4))

legend(0.7, 1.1, cex = 0.9, title = "Estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 2, col = c("#009900","red"), 
       bty = "n", horiz = F)
title("(a)", adj=0.025, line=0, cex.main=1.2)

# OTUs coloured according to tax rank
library(RColorBrewer)
# Generate vector with phylum names for node coloring
otu.netprop <- netcomi_netprops$lccNames1
otu.netprop <- as.data.frame(otu.netprop)
colnames(otu.netprop)[1] <- "OTU"
temptax2 <- as.data.frame(tax_table(physeq_prok_uparse_mSeq))
temptax2 <- subset(temptax2, row.names(temptax2) %in% otu.netprop$OTU)

temptax2$Phylum[temptax2$Phylum=="Firmicutes"]<-"Other"
phyla <- as.factor(temptax2$Phylum)
phyla <- factor(phyla, levels = c("Acidobacteria", "Actinobacteria", "Bacteroidetes", "Chloroflexi", "Gemmatimonadetes", "Planctomycetes", "Proteobacteria", "Thaumarchaeota", "Verrucomicrobia", "Other"))
names(phyla) <- rownames(temptax2)

dominant.names <- renamed$Dominant
names(dominant.names) <- renamed$Original

# Create color vector
colvec <- RColorBrewer::brewer.pal(length(levels(phyla)), "Set3")
plot(netcomi_netprops,
     repulsion = 0.98,
     rmSingles = "all",
     #shortenLabels = "intelligent",
     labels=new.names, labelScale = F,
     cexLabels = 0, cexHubLabels = 0.8, hubBorderWidth=0.8,
     #labelScale = FALSE, labels=F, (linia per deixar graf sense labels)
     #labels = dominant.names, labelFont = 2, cexLabel = 0.9, (linia per posar labels de spp dominants, enlloc de hubs)
     nodeSize = "mclr",
     nodeColor = "feature", 
     featVecCol = phyla, 
     colorVec =  colvec,
     nodeTransp = 20,
     highlightHubs = T,
     cexNodes = 1.2,
     edgeTranspHigh = 20,
     #title1 = "Data features highlighted", 
     #showTitle = TRUE,
     #cexTitle = 2,
     mar = c(1, 1, 4.5, 4))

col_transp <- colToTransp(colvec, 20)

legend(1.0, 0, cex = 0.75, pt.cex = 2, title = "Phylum:", y.intersp = 0.8, title.adj = 0.25,
       legend=levels(phyla), col = col_transp, bty = "n", pch = 16, xpd = T)
title("(b)", adj=0.025, line=0, cex.main=1.2)

# null model
set.seed(197702)
nullmodel <- sample_gnm(n = 185, m = 620, directed = FALSE, loops = FALSE)
plot(nullmodel)
# Degree distribution
degree_distribution(nullmodel)
ddist.null<- igraph::degree_distribution(nullmodel)
df.null <- data.frame(Degree = as.factor((seq_along(ddist.null)) - 1),
                      Fraction = ddist.null)
ggplot(data = df.null, aes(x = Degree, y = Fraction, group = 1)) +
  geom_line() +
  geom_point() +
  theme_bw()

testnull <- matrix(NA, nrow=1000, ncol=4)
colnames(testnull) <- c("degree.centr", "betw.centr", "path.length", "transitivity") 
for(i in 1:1000) {
  set.seed(197702+i)
  nullmodel <- sample_gnm(n = 185, m = 620, directed = FALSE, loops = FALSE)
  testnull[i,1] <- centr_degree(nullmodel, mode="in", normalized=T)$centralization
  testnull[i,2] <- centr_betw(nullmodel, directed=F, normalized=T)$centralization
  testnull[i,3] <- mean_distance(nullmodel)
  testnull[i,4] <- transitivity(nullmodel)
}
testnull <- as.data.frame(testnull)
quantile(testnull$degree.centr, probs = c(0.025, 0.5, 0.975))
centr_degree(netcomi_graph, mode="in", normalized=T)$centralization

quantile(testnull$betw.centr, probs = c(0.025, 0.5, 0.975))
centr_betw(netcomi_graph, directed=F, normalized=T)$centralization

quantile(testnull$path.length, probs = c(0.025, 0.5, 0.975))
mean_distance(netcomi_graph)

quantile(testnull$transitivity, probs = c(0.025, 0.5, 0.975))
transitivity(netcomi_graph)

# CO-CORRESPONDENCE ANALYSIS
metadatos_parc2 <- read.delim("metadatos_parc2.txt", row.names=1, sep = "\t", header=TRUE)
sample_data(physeq_fungi_uparse_mSeq) <- metadatos_parc2
head(sample_data(physeq_fungi_uparse_mSeq))
physeq_fungi_uparse_mSeq_peat <- subset_samples(physeq_fungi_uparse_mSeq, origin%in%c("peat"))
otu_fung_peat <- as.data.frame(otu_table(physeq_fungi_uparse_mSeq_peat))
otu_fung_peat_p <- otu_fung_peat
for (i in 1:ncol(otu_fung_peat)){
  otu_fung_peat_p[,i] <- otu_fung_peat[,i]/sum(otu_fung_peat[,i])
}
otu_fung_peat_t <- as.data.frame(t(otu_fung_peat_p))
otu_fung_peat_t <- otu_fung_peat_t[, colSums(otu_fung_peat_t) != 0]
otu_fung_peat_t <- otu_fung_peat_t[-c(31, 41), ]
otu_prok_coca_t <- otu_prok_t[-c(31, 32, 40), ]
otu_fung_peat_t <- otu_fung_peat_t[ order(row.names(otu_fung_peat_t)), ]
otu_prok_coca_t <- otu_prok_coca_t[ order(row.names(otu_prok_coca_t)), ]
otu_fung_peat_t <- otu_fung_peat_t[, colSums(otu_fung_peat_t) != 0]
otu_prok_coca_t <- otu_prok_coca_t[, colSums(otu_prok_coca_t) != 0]

otu_fung_peat_t2 <- otu_fung_peat_t[, colSums(otu_fung_peat_t > 0) > 10 ] #appear.freq>25%, 144 fungi
otu_prok_coca_t2 <- otu_prok_coca_t[, colSums(otu_prok_coca_t > 0) > 31 ] #appear.freq>75%, 186 prok
otu_prok_coca_t3 <- otu_prok_coca_t2[, colMeans(otu_prok_coca_t2) > 0.001]
colnames(otu_prok_coca_t3) <-gsub("BOTU_","",as.character(colnames(otu_prok_coca_t3)))
colnames(otu_fung_peat_t2) <-gsub("OTU_","",as.character(colnames(otu_fung_peat_t2)))

library(cocorresp)
bact.fung <- cocorresp::coca(y = otu_prok_coca_t3, x = otu_fung_peat_t2, method = "symmetric")
summary(bact.fung)
par(mfcol=c(1,1))
screeplot(bact.fung)
cocorresp::corAxis(bact.fung)
cocorresp::crossval(otu_prok_coca_t3, otu_fung_peat_t2)

#percent variance explained = eigenvalue/sumatori(eigenvalues)
bact.fung$lambda[1]/sum(bact.fung$lambda)
bact.fung$lambda[2]/sum(bact.fung$lambda)
#percent variance explained for each community with 2 components
bact.fung2 <- cocorresp::coca(y = otu_prok_coca_t3, x = otu_fung_peat_t2, method = "symmetric", n.axes=2)
summary(bact.fung2)
# 44.3% explained for prok, 20.8% for fungi

# correlation sites scores between bact-fung in firs two axis
cor.test(bact.fung$scores$site$Y[,1], bact.fung$scores$site$X[,1])
cor.test(bact.fung$scores$site$Y[,2], bact.fung$scores$site$X[,2])

coca.bact.sign <- data.frame(bact.fung.scores$species$Y[,1:2], colnames(otu_prok_coca_t3))
colnames(coca.bact.sign)[3] <- "OTU"
coca.fung.sign <- data.frame(bact.fung.scores$species$X[,1:2], colnames(otu_fung_peat_t2))
colnames(coca.fung.sign)[3] <- "OTU"
attach(coca.bact.sign)
coca.bact.sign$sign1[abs(COCA.1) >= 0.5] <- 1
coca.bact.sign$sign1[abs(COCA.1) < 0.5] <- 0
detach(coca.bact.sign)

attach(coca.bact.sign)
coca.bact.sign$sign2[abs(COCA.2) >= 0.75] <- 1
coca.bact.sign$sign2[abs(COCA.2) < 0.75] <- 0
detach(coca.bact.sign)
coca.bact.sign$sign <- coca.bact.sign$sign1 + coca.bact.sign$sign2

attach(coca.fung.sign)
coca.fung.sign$sign1[abs(COCA.1) >= 0.7] <- 1
coca.fung.sign$sign1[abs(COCA.1) < 0.7] <- 0
detach(coca.fung.sign)

attach(coca.fung.sign)
coca.fung.sign$sign2[abs(COCA.2) >= 1] <- 1
coca.fung.sign$sign2[abs(COCA.2) < 1] <- 0
detach(coca.fung.sign)
coca.fung.sign$sign <- coca.fung.sign$sign1 + coca.fung.sign$sign2

site.scores <- data.frame(bact.fung$scores$site$Y[,1:2])
site.scores.fung <- data.frame(bact.fung$scores$site$X[,1:2])

# Figure 6
library(ggrepel)
library(ggnewscale)
site.scores2 <- merge(x = site.scores, y = metadatos_parc2, by = "row.names", all.x=TRUE)
site.scores2$year <- as.factor(site.scores2$year)
site.scores2$location <- as.factor(site.scores2$location)
colnames(site.scores2)[c(11, 12)] <- c("Year", "Plot")

site.scores.fung2 <- merge(x = site.scores.fung, y = metadatos_parc2, by = "row.names", all.x=TRUE)
site.scores.fung2$year <- as.factor(site.scores.fung2$year)
site.scores.fung2$location <- as.factor(site.scores.fung2$location)
colnames(site.scores.fung2)[c(11, 12)] <- c("Year", "Plot")

ggplot(coca.bact.sign, aes(x=COCA.1, y=COCA.2)) +
  theme_classic() +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_point(aes(), color="brown1", size=2, shape=20) +
  geom_text_repel(data = subset(coca.bact.sign, sign > 0),
                  aes(label = OTU), color="black", size=4, max.overlaps = 12) +
  theme(legend.position = "none") +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  xlab("CoCA 1 (53.6%)") + ylab("CoCA 2 (20.4%)") + new_scale("shape") +
  geom_point(data=site.scores2, aes(x=COCA.1, y=COCA.2, shape=Plot, fill=Year),
             color="grey60", size=2.5) +
  scale_shape_manual(values=c(21, 24)) +
  scale_fill_manual(values=c("white", "grey80")) +
  ggtitle("(a) Prokaryotes")

ggplot(coca.fung.sign, aes(x=COCA.1, y=COCA.2)) +
  theme_classic() +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_point(aes(), color="brown1",size=2, shape=20) +
  geom_text_repel(data = subset(coca.fung.sign, sign > 0), 
                  aes(label = OTU), color="black", size=4, max.overlaps = 12) +
  theme(legend.position = "none") +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  xlab("CoCA 1 (53.6%)") + ylab("CoCA 2 (20.4%)") + new_scale("shape") +
  geom_point(data=site.scores.fung2, aes(x=COCA.1, y=COCA.2, shape=Plot, fill=Year),
             color="grey60", size=2.5) +
  scale_shape_manual(values=c(21, 24)) +
  scale_fill_manual(values=c("white", "grey80")) +
  ggtitle("(b) Fungi")

temptax[rownames(temptax) %in% c("BOTU_33", "BOTU_10", "BOTU_7", "BOTU_62", "BOTU_369", "BOTU_66", "BOTU_108", "BOTU_234"), ]
temptax.f <- as.data.frame(tax_table(physeq_fungi_uparse_mSeq))
temptax.f[rownames(temptax.f) %in% c("OTU_64", "OTU_19", "OTU_1183", "OTU_241", "OTU_112", "OTU_172"), ]
temptax[rownames(temptax) %in% c("BOTU_3", "BOTU_14", "BOTU_71", "BOTU_147", "BOTU_75"), ]
temptax.f[rownames(temptax.f) %in% c("OTU_80", "OTU_102", "OTU_332", "OTU_180", "OTU_883"), ]
temptax[rownames(temptax) %in% c("BOTU_82", "BOTU_120", "BOTU_166", "BOTU_86", "BOTU_54", "BOTU_711", "BOTU_48", "BOTU_7962"), ]
temptax.f[rownames(temptax.f) %in% c("OTU_32", "OTU_70", "OTU_1774", "OTU_15", "OTU_3285", "OTU_2052"), ]

# Correlations between richness/Shannon index in the bact-fung communities
metadatos_parc2 <- read.delim("metadatos_parc2.txt", row.names=1, sep = "\t", header=TRUE)
sample_data(physeq_fungi_uparse_mSeq) <- metadatos_parc2
head(sample_data(physeq_fungi_uparse_mSeq))
physeq_fungi_uparse_mSeq_PCoA <- subset_samples(physeq_fungi_uparse_mSeq, PCoA%in%c(1))
physeq_fungi_uparse_mSeq_PCoA2 <- transform_sample_counts(physeq_fungi_uparse_mSeq_PCoA, function(x) round(x, digits=0))
sample_data(physeq_fungi_uparse_mSeq_PCoA2)$etiq <- factor(sample_data(physeq_fungi_uparse_mSeq_PCoA2)$etiq,
                                                           levels = c("P0M", "P0T", "P1M", "P1T", "P2M", "P2T", "P3", "P3*",
                                                                      "I1M", "I1T", "I2M", "I2T", "I3", "S1M", "S1T", "S2M", "S2T"))
divers.index.fung <- estimate_richness(physeq_fungi_uparse_mSeq_PCoA2, measures = c("Observed", "Shannon"))
divers.index.fung$n_sample <- rownames(divers.index.fung)
parlade140 <- read.table(file="parlade_nidos_140m.txt",
                         header=T, sep="\t")
parlade140$soiltype <- factor(parlade140$soiltype)
parlade140$age <- factor(parlade140$age)
parlade140$etiq <- factor(parlade140$etiq)
parlade140$plot <- factor(parlade140$plot)
parlade140$logmyc <- log(parlade140$mycelium + 0.001)

divers.index.fung2 <- merge(x = divers.index.fung, y = parlade140, by = "n_sample", all = TRUE)
divers.index.fung2$Observed[c(9, 10, 13)] <- 0
divers.index.fung2$Shannon[c(9, 10, 13)] <- 0
colnames(divers.index.fung2)[2:3] <- c("fung_rich", "fung_Shannon")

divers.bact.fung <- merge(x = divers.index, y = divers.index.fung2, by = "n_sample", all.x = TRUE)
colnames(divers.bact.fung)[2:3] <- c("prok_rich", "prok_Shannon")
divers.bact.fung <- divers.bact.fung[-c(43:44, 48:52), ] #elimino pozo8 pero tambien las muestras de sustrato

lm.prok.fung.rich <- lm(prok_rich ~ sqrt(fung_rich+0.5), divers.bact.fung)
summary(lm.prok.fung.rich)
par(mfcol=c(2,2))
plot(lm.prok.fung.rich)
shapiro.test(residuals(lm.prok.fung.rich))
cor.test(divers.bact.fung$prok_rich, sqrt(divers.bact.fung$fung_rich+0.5))

lm.prok.fung.shannon <- lm(prok_Shannon ~ fung_Shannon, divers.bact.fung)
summary(lm.prok.fung.shannon)
par(mfcol=c(2,2))
plot(lm.prok.fung.shannon)
shapiro.test(residuals(lm.prok.fung.shannon))
cor.test(divers.bact.fung$prok_Shannon, divers.bact.fung$fung_Shannon)
