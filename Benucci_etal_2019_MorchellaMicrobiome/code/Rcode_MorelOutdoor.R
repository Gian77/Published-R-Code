# ******* DATA ANALYSIS ********************************************** ----------------------------
# Manuscript:   Microbial Communities Associated with Cultivated Morel OUTDOOR
# Authors:      Benucci GMN, Longley R, ... Bonito G.
# Affiliation:  Michigan State University
# Journal:      PeerJ
# Date:         June 6, 2019
# ******************************************************************** ----------------------------

# ___________WORKING ENVIRONMENT SETUP _____________ ----------------------------------------------

options(scipen = 999) #to use decimals
options(max.print=100000000) # to print more lines on the display
options(verbose=TRUE)
#rm(list = ls()) # use with caution

# >>> IMPORTING DATASETS -------------------------------------------------------------------------
library(phyloseq); packageVersion("phyloseq")
library(ggplot2)
library(Biostrings)
library(ape)

# Import funbgal data ----------------------------------------------------------------------------
biom_ITS_uparse = import_biom("otu_table_tax_R1_ITS_json.biom")
map_ITS = import_qiime_sample_data("mapping_ITS_new.txt")
sample_data(biom_ITS_uparse) <- map_ITS
colnames(tax_table(biom_ITS_uparse)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
otus_rep_ITS_uparse <- readDNAStringSet("otus_ITS.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
biom_ITS_uparse <- merge_phyloseq(biom_ITS_uparse, otus_rep_ITS_uparse)
refseq(biom_ITS_uparse)
biom_ITS_uparse

tax_table(biom_ITS_uparse)[tax_table(biom_ITS_uparse)==""]<- NA
tax_table(biom_ITS_uparse)[is.na(tax_table(biom_ITS_uparse))]<-"Unclassified"

# filtering taxonomies ---------------------------------------------------------------------------
unique(as.data.frame(tax_table(biom_ITS_uparse))$Kingdom)
unique(as.data.frame(tax_table(biom_ITS_uparse))$Phylum) # found Cercozoa

biom_ITS_uparse <- subset_taxa(biom_ITS_uparse, Phylum!="Cercozoa")

# import prokaryotic data ------------------------------------------------------------------------
biom_16s_uparse = import_biom("otu_table_tax_R1_16S_json.biom")
map_16s = import_qiime_sample_data("mapping_16s_new.txt")
sample_data(biom_16s_uparse) <- map_16s
colnames(tax_table(biom_16s_uparse)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
otus_rep_16s_uparse <- readDNAStringSet("otus_16s.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
biom_16s_uparse <- merge_phyloseq(biom_16s_uparse, otus_rep_16s_uparse)
biom_16s_uparse

# relabeling taxonomies --------------------------------------------------------------------------
tax_table(biom_16s_uparse)[, "Kingdom"] <- gsub("k__", "", tax_table(biom_16s_uparse)[, "Kingdom"])
tax_table(biom_16s_uparse)[, "Phylum"] <- gsub("p__", "", tax_table(biom_16s_uparse)[, "Phylum"])
tax_table(biom_16s_uparse)[, "Class"] <- gsub("c__", "", tax_table(biom_16s_uparse)[, "Class"])
tax_table(biom_16s_uparse)[, "Order"] <- gsub("o__", "", tax_table(biom_16s_uparse)[, "Order"])
tax_table(biom_16s_uparse)[, "Family"] <- gsub("f__", "", tax_table(biom_16s_uparse)[, "Family"])
tax_table(biom_16s_uparse)[, "Genus"] <- gsub("g__", "", tax_table(biom_16s_uparse)[, "Genus"])
tax_table(biom_16s_uparse)[, "Species"] <- gsub("s__", "", tax_table(biom_16s_uparse)[, "Species"])

tax_table(biom_16s_uparse)[tax_table(biom_16s_uparse)==""]<- NA
tax_table(biom_16s_uparse)[is.na(tax_table(biom_16s_uparse))]<-"Unclassified"

# filtering taxonomies ---------------------------------------------------------------------------
sort(unique(as.data.frame(tax_table(biom_16s_uparse))$Kingdom))
sort(unique(as.data.frame(tax_table(biom_16s_uparse))$Phylum))
sort(unique(as.data.frame(tax_table(biom_16s_uparse))$Class)) # found Chloroplast
sort(unique(as.data.frame(tax_table(biom_16s_uparse))$Class))
sort(unique(as.data.frame(tax_table(biom_16s_uparse))$Family)) # found mitochondria

biom_16s_uparse <- subset_taxa(biom_16s_uparse, Class!="Chloroplast")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Family!="mitochondria")

# Filtering out OTUs < than 10 reads ------------------------------------------------------------

biom_ITS_uparse -> biom_ITS_uparse_filt
otu_table(biom_ITS_uparse_filt) <- otu_table(biom_ITS_uparse_filt)[which(rowSums(otu_table(biom_ITS_uparse_filt)) >= 10),] ### PCR Errors 
biom_ITS_uparse_filt

biom_16s_uparse -> biom_16s_uparse_filt
otu_table(biom_16s_uparse_filt) <- otu_table(biom_16s_uparse_filt)[which(rowSums(otu_table(biom_16s_uparse_filt)) >= 10),] ### PCR Errors 
biom_16s_uparse_filt

# function to remove taxa present in the negative controls ---------------------------------------
remove_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}


badTaxa_ITS =c("OTU_17", "OTU_51")
biom_ITS_uparse_filt = remove_taxa(biom_ITS_uparse_filt, badTaxa_ITS)
biom_ITS_uparse_filt

badTaxa_16s =c("OTU_260","OTU_1933", "OTU_245","OTU_1155","OTU_2789","OTU_2693","OTU_9938","OTU_2645",
               "OTU_348","OTU_1790","OTU_414","OTU_2021","OTU_11","OTU_3102","OTU_5590","OTU_1530",
               "OTU_111","OTU_811","OTU_78","OTU_9041","OTU_104","OTU_8262","OTU_2701","OTU_5","OTU_150",
               "OTU_112","OTU_68","OTU_460","OTU_161","OTU_786","OTU_883","OTU_532","OTU_1149","OTU_889",
               "OTU_1244","OTU_482","OTU_1001","OTU_1172","OTU_692","OTU_957","OTU_1127","OTU_4226",
               "OTU_601","OTU_977")
biom_16s_uparse_filt = remove_taxa(biom_16s_uparse_filt, badTaxa_16s)
biom_16s_uparse_filt

# ****************************************************************--------------------------------
# ___ Extractring samples of the OUTDOOR Study ____-----------------------------------------------
library("dplyr")

biom_ITS_uparse_out -> physeq_fungi
rownames_fungi_out <- c("samplefun17","samplefun19","samplefun15","samplefun16","samplefun25",
"samplefun24","samplefun22","samplefun13","samplefun20","samplefun23")
otu_table(physeq_fungi) <- subset(otu_table(physeq_fungi), select = rownames_fungi_out)
otu_table(physeq_fungi) <- otu_table(physeq_fungi)[which(rowSums(otu_table(physeq_fungi)) >= 1),] 
physeq_fungi
sample_data(physeq_fungi)

sample_data(biom_16s_uparse_filt)
biom_16s_uparse_out -> physeq_prokaryote
rownames_prokaryote_out<- c("samplebac26","samplebac28","samplebac295","samplebac297","samplebac30",
                           "samplebac37","samplebac39","samplebac43","samplebac45","samplebac49",
                           "samplebac64","samplebac65","samplebac66","samplebac67","samplebac69",
                           "samplebac70","samplebac72","samplebac73","samplebac74","samplebac75",
                           "samplebac27","samplebac29","samplebac31","samplebac34","samplebac38",
                           "samplebac40","samplebac44","samplebac46","samplebac48","samplebac50") 
otu_table(physeq_prokaryote) <- subset(otu_table(physeq_prokaryote), select = rownames_prokaryote_out)
otu_table(physeq_prokaryote) <- otu_table(physeq_prokaryote)[which(rowSums(otu_table(physeq_prokaryote)) >= 1),] 
physeq_prokaryote
sample_data(physeq_prokaryote)

# create vegan objects 
otu_fungi_out <- as.data.frame(otu_table(physeq_fungi))
taxa_fungi_out <- as.data.frame(as.matrix(tax_table(physeq_fungi)))
metadata_fungi_out <- as.data.frame(as.matrix(sample_data(physeq_fungi)))
dim(otu_fungi_out)

otu_prokaryote_out <- as.data.frame(otu_table(physeq_prokaryote))
taxa_prokaryote_out <- as.data.frame(as.matrix(tax_table(physeq_prokaryote)))
metadata_prokaryote_out <- as.data.frame(as.matrix(sample_data(physeq_prokaryote)))
dim(otu_prokaryote_out)

# >>> RAREFACTION CURVES -----------------------------------------------------------------------------
library("vegan")

# *** FIGURE S2 rarecurve prokaryotes  -----------------------------------------------------------
rarecurve(t(otu_prokaryote_out), col = metadata_prokaryote_out$Origin, label = FALSE, 
          sample=min(colSums(otu_prokaryote_out)), step = 50,
          main="Prokaryotes", ylab = "Number of OTUs", xlab = "Number of DNA reads") -> rare_fungi_out2
legend("bottomright", legend=c("Cap", "Stem", "Soil"),
       col=c("black", "green", "red"), lty=1, cex=0.8, box.lty=1) #box.lty=0 remove the legend border

# *** FIGURE S3 rarecurve prokaryotes  -----------------------------------------------------------
rarecurve(t(otu_fungi_out), col = metadata_fungi_out$Stage, label = FALSE, sample=min(colSums(otu_fungi_out)), step = 50,
          main="Fungi", ylab = "Number of OTUs", xlab = "Number of DNA reads", cex=0.6) -> rare_fungi
legend("bottomright", legend=c("Young", "Mature"),
       col=c("red", "black"), lty=1, cex=0.8, box.lty=1) #box.lty=0 remove the legend border


# MetagenomeSeq normalization --------------------------------------------------------------------
library("metagenomeSeq")

# fitting into a Gaussian Model using metagenomeSeq Fungi 
set.seed(2019)
physeq_fungi_CSS = phyloseq_to_metagenomeSeq(physeq_fungi)
p_physeq_fungi<-round(cumNormStatFast(physeq_fungi_CSS,pFlag=TRUE),digits=2)
physeq_quant<-cumNorm(physeq_fungi_CSS, p=p_physeq_fungi)
physeq_quant
normFactors(physeq_quant)
physeq_fungi_CSS<-MRcounts(physeq_quant, norm=T)
head(physeq_fungi_CSS)

physeq_fungi_mSeq <- physeq_fungi
otu_table(physeq_fungi_mSeq) <- otu_table(physeq_fungi_CSS, taxa_are_rows=TRUE)
physeq_fungi_mSeq

otu_fungi_out_mSeq <- as.data.frame(otu_table(physeq_fungi_mSeq))
taxa_fungi_out_mSeq <- as.data.frame(as.matrix(tax_table(physeq_fungi_mSeq)))
metadata_fungi_out_mSeq <- as.data.frame(as.matrix(sample_data(physeq_fungi_mSeq)))
dim(otu_fungi_out_mSeq)

# fitting into a Gaussian Model using metagenomeSeq Prokaryotes
physeq_prokaryote_CSS = phyloseq_to_metagenomeSeq(physeq_prokaryote)
p_physeq_prokaryote<-round(cumNormStatFast(physeq_prokaryote_CSS,pFlag=TRUE),digits=2)
physeq_quant<-cumNorm(physeq_prokaryote_CSS, p=p_physeq_prokaryote)
physeq_quant
normFactors(physeq_quant)
physeq_prokaryote_CSS<-MRcounts(physeq_quant, norm=T)
head(physeq_prokaryote_CSS)

physeq_prokaryote_mSeq <- physeq_prokaryote
otu_table(physeq_prokaryote_mSeq) <- otu_table(physeq_prokaryote_CSS, taxa_are_rows=TRUE)
physeq_prokaryote_mSeq

otu_prokaryote_out_mSeq <- as.data.frame(otu_table(physeq_prokaryote_mSeq))
taxa_prokaryote_out_mSeq <- as.data.frame(as.matrix(tax_table(physeq_prokaryote_mSeq)))
metadata_prokaryote_out_mSeq <- as.data.frame(as.matrix(sample_data(physeq_prokaryote_mSeq)))
dim(otu_prokaryote_out_mSeq)

#>>> BETA DIVERSITY ------------------------------------------------------------------------------
library("ggrepel")

colnames(sample_data(physeq_fungi_mSeq))[6] <- "Origin"
sample_data(physeq_fungi_mSeq) <- sample_data(physeq_fungi_mSeq)[,c(1,6,8:9)]
sample_data(physeq_fungi_mSeq)

pcoa_fungi_out = phyloseq::ordinate(physeq_fungi_mSeq, method ="PCoA", distance="bray")

p_pcoa_ITS_out = plot_ordination(physeq_fungi_mSeq, pcoa_fungi_out, color="Stage", shape = "Origin") + 
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2, alpha=1) + # ,aes(shape=Age))
  scale_colour_manual(values=c("red","grey")) +
  scale_shape_manual(values=c(2, 2)) +
  #stat_ellipse(aes(group=Stage), type="norm", alpha=0.8, linetype = 2, show.legend = FALSE) +
  #geom_text_repel(aes(label=sample_data(physeq_fungi_mSeq)$Description), size = 2) + 
  theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  theme(legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8)) +
  #theme(legend.title=element_blank())
  theme(legend.position="bottom") 

p_pcoa_ITS_out

sample_data(physeq_prokaryote_mSeq)
sample_data(physeq_prokaryote_mSeq) <- sample_data(physeq_prokaryote_mSeq)[,c(1,4,6:7)]
sample_data(physeq_prokaryote_mSeq)$Origin <- factor(sample_data(physeq_prokaryote_mSeq)$Origin,levels=c("Cap","Stem","Soil"))

pcoa_prokaryote_out = phyloseq::ordinate(physeq_prokaryote_mSeq, method ="PCoA", distance="bray")

p_pcoa_16s_out = plot_ordination(physeq_prokaryote_mSeq, pcoa_prokaryote_out, color="Stage", shape = "Origin") + 
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2, alpha=1) + # ,aes(shape=Age))
  scale_colour_manual(values=c("red", "grey")) +
  scale_shape_manual(values=c(0,1,2)) +
  #geom_text_repel(aes(label=sample_data(physeq_prokaryote_mSeq)$Description), size = 2) + 
  #stat_ellipse(aes(group=Origin), type="norm", alpha=0.8, linetype = 2, show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8)) +
  #theme(legend.title=element_blank())
  theme(legend.position="bottom") 

p_pcoa_16s_out

# *** FIGURE 2 - ordinations ---------------------------------------------------------------------
library("ggpubr")

ggarrange(p_pcoa_16s_out,
          p_pcoa_ITS_out,
          labels = c("A", "B"),
          widths = c(1,1),
          align = "none", 
          ncol = 2, nrow = 1,
          common.legend = TRUE,
          legend = c("right"))


# >>> PERMANOVA ----------------------------------------------------------------------------------

# fungal communities 
vegan::vegdist(t(otu_fungi_out_mSeq), method="bray") -> dist_fungi

model.matrix(~ Stage, data=metadata_fungi_out_mSeq)
adonis(dist_fungi ~ Stage, data=metadata_fungi_out_mSeq, permutations=999) -> adonis_fungi
adonis_fungi

densityplot(permustats(adonis_fungi), main=list(label="Fungi var. Stage", cex=1)) -> density_fungi
density_fungi

betadisper(dist_fungi, metadata_fungi_out_mSeq$Stage) -> betadisper_fungi
anova(betadisper_fungi, permutations = 9999)
permutest(betadisper_fungi, permutations = 9999, pairwise = T)
plot(betadisper_fungi)
plot(TukeyHSD(betadisper_fungi), las=0)
boxplot(betadisper_fungi)

# prokaryotic communities 
vegan::vegdist(t(otu_prokaryote_out_mSeq), method="bray") -> dist_prokaryote

model.matrix(~ Stage * Origin, data=metadata_prokaryote_out_mSeq)
adonis(dist_prokaryote ~ Stage * Origin, data=metadata_prokaryote_out_mSeq, permutations=999) -> adonis_prokaryote
adonis_prokaryote
densityplot(permustats(adonis_prokaryote), main=list(label="Prokaryote var. Stage * Origin", cex=1)) -> density_prokaryote
density_prokaryote

betadisper(dist_prokaryote, metadata_prokaryote_out_mSeq$Stage) -> betadisper_prokaryote_S
anova(betadisper_prokaryote_S, permutations = 9999)
permutest(betadisper_prokaryote_S, permutations = 9999, pairwise = T)
plot(betadisper_prokaryote_S)
plot(TukeyHSD(betadisper_prokaryote_S), las=0)
boxplot(betadisper_prokaryote_S)

betadisper(dist_prokaryote, metadata_prokaryote_out_mSeq$Origin) -> betadisper_prokaryote_O
anova(betadisper_prokaryote_O, permutations = 9999)
permutest(betadisper_prokaryote_O, permutations = 9999, pairwise = T)
plot(betadisper_prokaryote_O)
plot(TukeyHSD(betadisper_prokaryote_O), las=0)
boxplot(betadisper_prokaryote_O)

# Figure S4 - betadisper -------------------------------------------------------------------------
par(mfrow=c(3,2))
plot(betadisper_prokaryote_S, main="(A) Prokaryotes\n\nPCoA (Stage)", las=1) #, cex.main=1
boxplot(betadisper_prokaryote_S, main="\nboxplot (Stage)", xlab="Stage")
#plot(TukeyHSD(betadisper_prokaryote_S), las=0)
plot(betadisper_prokaryote_O, main="\nPCoA (Origin)", las=1) #, cex.main=1
boxplot(betadisper_prokaryote_O, main="\nboxplot (Origin)", xlab="Stage")
#plot(TukeyHSD(betadisper_prokaryote_O), las=0)
plot(betadisper_fungi, main="(B) Fungi\n\nPCoA (Stage)", las=1) #, cex.main=1
boxplot(betadisper_fungi, main="\nboxplot (Stage)", xlab="Stage")
#plot(TukeyHSD(betadisper_fungi), las=0)
dev.off()


# >>> ALPHA DIVERSITY ----------------------------------------------------------------------------
library("BiodiversityR")

dim(otu_fungi_out)
dim(metadata_fungi_out)
count(metadata_fungi_out, vars = Stage)

alpha_div_fungi <- metadata_fungi_out[,c(1,8)]
alpha_div_fungi$readNO <- sample_sums(physeq_fungi)
alpha_div_fungi$Observed <- specnumber(otu_fungi_out, MARGIN = 2)
alpha_div_fungi$Rarefied <- rarefy(otu_fungi_out,sample=min(alpha_div_fungi$readNO), MARGIN = 2)
alpha_div_fungi$Shannon <- diversity(otu_fungi_out, index="shannon", MARGIN = 2)
jevenness_fungi <- diversityresult(t(otu_fungi_out), method = "each site", index = "Jevenness")
alpha_div_fungi$Jevenness <- jevenness_fungi$Jevenness
alpha_div_fungi <- alpha_div_fungi[order(sums_fungi$ReadNO), ]
alpha_div_fungi

# get descriptive stats
library("psych")
describeBy(alpha_div_fungi, alpha_div_fungi$Stage)

# check homogenity of variances 
# check homogenity of variances 
fligner.test(Observed ~ Stage, data=alpha_div_fungi)
fligner.test(Shannon ~ Stage, data=alpha_div_fungi)
fligner.test(Jevenness ~ Stage, data=alpha_div_fungi)

# get significant differences
aov_fungi_rich <- aov(Observed ~ Stage, data=alpha_div_fungi)
summary(aov_fungi_rich)
aov_fungi_shan <- aov(Shannon ~ Stage, data=alpha_div_fungi)
summary(aov_fungi_shan)
aov_fungi_Jeven <- aov(Jevenness ~ Stage, data=alpha_div_fungi)
summary(aov_fungi_Jeven)

dim(otu_prokaryote_out)
dim(metadata_prokaryote_out)
count(metadata_prokaryote_out, vars = Stage)
count(metadata_prokaryote_out, vars = Origin)

alpha_div_prokaryote <- metadata_prokaryote_out[,c(1,4,6:7)]
alpha_div_prokaryote$readNO <- sample_sums(physeq_prokaryote)
alpha_div_prokaryote$Observed <- specnumber(otu_prokaryote_out, MARGIN = 2)
alpha_div_prokaryote$Rarefied <- rarefy(otu_prokaryote_out,sample=min(alpha_div_prokaryote$readNO), MARGIN = 2)
alpha_div_prokaryote$Shannon <- diversity(otu_prokaryote_out, index="shannon", MARGIN = 2)
jevenness_prokaryote <- diversityresult(t(otu_prokaryote_out), method = "each site", index = "Jevenness")
alpha_div_prokaryote$Jevenness <- jevenness_prokaryote$Jevenness
alpha_div_prokaryote <- alpha_div_prokaryote[order(sums_prokaryote$ReadNO), ]
alpha_div_prokaryote

describeBy(alpha_div_prokaryote, alpha_div_prokaryote$Stage)
describeBy(alpha_div_prokaryote, alpha_div_prokaryote$Origin)

# check homogenity of variances 
fligner.test(Observed ~ Origin, data=alpha_div_prokaryote)
fligner.test(Shannon ~ Origin, data=alpha_div_prokaryote)
fligner.test(Jevenness ~ Origin, data=alpha_div_prokaryote)

# get significant differences
library("agricolae")

aov_prokaryote_rich <- aov(Observed ~ Stage, data=alpha_div_prokaryote)
summary(aov_prokaryote_rich)
aov_prokaryote_shan <- aov(Shannon ~ Stage, data=alpha_div_prokaryote)
summary(aov_prokaryote_shan)
aov_prokaryote_Jeven <- aov(Jevenness ~ Stage, data=alpha_div_prokaryote)
summary(aov_prokaryote_Jeven)

aov_prokaryote_rich <- aov(Observed ~ Origin, data=alpha_div_prokaryote)
summary(aov_prokaryote_rich)
HSD.test(aov_prokaryote_rich, "Origin") -> tukeyHSD_prokaryote_rich
tukeyHSD_prokaryote_rich

aov_prokaryote_shan <- aov(Shannon ~ Origin, data=alpha_div_prokaryote)
summary(aov_prokaryote_shan)
HSD.test(aov_prokaryote_shan, "Origin") -> tukeyHSD_prokaryote_shan
tukeyHSD_prokaryote_shan

aov_prokaryote_Jeven <- aov(Jevenness ~ Origin, data=alpha_div_prokaryote)
summary(aov_prokaryote_Jeven)
HSD.test(aov_prokaryote_Jeven, "Origin") -> tukeyHSD_prokaryote_Jeven
tukeyHSD_prokaryote_Jeven


# >>> INDICATOR SPECIES ANALYSIS (ISA) -----------------------------------------------------------
library("indicspecies")

isa_fungi <- multipatt(as.data.frame(t(otu_fungi_out)), metadata_fungi_out$Stage, control=how(nperm=9999))
summary(isa_fungi, indvalcomp=TRUE)
isa_fungi -> isa_fungi_fdr
isa_fungi_fdr$sign$p.value<-p.adjust(isa_fungi$sign$p.value, "fdr")
summary(isa_fungi_fdr)

isa_prokaryote_ST <- multipatt(as.data.frame(t(otu_prokaryote_out)), metadata_prokaryote_out$Stage, control=how(nperm=9999))
summary(isa_prokaryote_ST, indvalcomp=TRUE)
isa_prokaryote_ST -> isa_prokaryote_fdr_ST
isa_prokaryote_fdr_ST$sign$p.value<-p.adjust(isa_prokaryote_ST$sign$p.value, "fdr")
summary(isa_prokaryote_fdr_ST)

isa_prokaryote_OR <- multipatt(as.data.frame(t(otu_prokaryote_out)), metadata_prokaryote_out$Origin, control=how(nperm=9999))
summary(isa_prokaryote_OR, indvalcomp=TRUE)
isa_prokaryote_OR -> isa_prokaryote_fdr_OR
isa_prokaryote_fdr_OR$sign$p.value<-p.adjust(isa_prokaryote_OR$sign$p.value, "fdr")
summary(isa_prokaryote_fdr_OR)

sink(file="isa_prokaryote_fdr_Origin.csv") 
summary(isa_prokaryote_fdr_OR)
sink()

# extracting ISA OTUs ----------------------------------------------------------------------------
result_isa_prokaryote_fdr_OR <- isa_prokaryote_fdr_OR$sign[which(isa_prokaryote_fdr_OR$sign$p.value <= 0.05), ]
head(result_isa_prokaryote_fdr_OR)
dim(result_isa_prokaryote_fdr_OR)

result_isa_prokaryote_fdr_OR[result_isa_prokaryote_fdr_OR$s.Cap==1 &
                               result_isa_prokaryote_fdr_OR$s.Stem==0 &
                               result_isa_prokaryote_fdr_OR$s.Soil==0 ,] -> Cap

result_isa_prokaryote_fdr_OR[result_isa_prokaryote_fdr_OR$s.Cap==0 &
                               result_isa_prokaryote_fdr_OR$s.Stem==1 &
                               result_isa_prokaryote_fdr_OR$s.Soil==0 ,] -> Stem

result_isa_prokaryote_fdr_OR[result_isa_prokaryote_fdr_OR$s.Cap==1 &
                               result_isa_prokaryote_fdr_OR$s.Stem==1 &
                               result_isa_prokaryote_fdr_OR$s.Soil==0 ,] -> Cap_Stem

result_isa_prokaryote_fdr_OR[result_isa_prokaryote_fdr_OR$s.Cap==1 &
                               result_isa_prokaryote_fdr_OR$s.Stem==0 &
                               result_isa_prokaryote_fdr_OR$s.Soil==1 ,] -> Cap_Soil

result_isa_prokaryote_fdr_OR[result_isa_prokaryote_fdr_OR$s.Cap==0 &
                               result_isa_prokaryote_fdr_OR$s.Stem==1 &
                               result_isa_prokaryote_fdr_OR$s.Soil==1 ,] -> Soil_Stem

result_isa_prokaryote_fdr_OR[result_isa_prokaryote_fdr_OR$s.Cap==0 &
                               result_isa_prokaryote_fdr_OR$s.Stem==0 &
                               result_isa_prokaryote_fdr_OR$s.Soil==1 ,] -> Soil


isa_df <- rbind(Cap, Stem, Cap_Stem, Cap_Soil, Soil_Stem)
dim(isa_df)
isa_df

# phyloseq objets of ISA OTUs --------------------------------------------------------------------
physeq_prokaryote -> physeq_prokaryote_isa
physeq_prokaryote_isa = transform_sample_counts(physeq_prokaryote_isa, function(x) 100*x/sum(x)) # transform to relativa abundances 
otu_table(physeq_prokaryote_isa) <-otu_table(physeq_prokaryote_isa_ab)[rownames(isa_df), ]
physeq_prokaryote_isa
sample_data(physeq_prokaryote_isa)

# improving taxonomy for indicator OTUs by BLAST
#write.dna(refseq(physeq_prokaryote_isa), format="fasta", colsep="", file="physeq_prokaryote_isa.fasta")
#write.csv(tax_table(physeq_prokaryote_isa), file = "tax_table_physeq_prokaryote_isa.csv")
tax_physeq_prokaryote_isa <- read.csv("tax_table_physeq_prokaryote_isa.csv", header=T, row.names =1)
tax_table(physeq_prokaryote_isa) <- tax_table(as.matrix(tax_physeq_prokaryote_isa))
#write.csv(sample_data(physeq_prokaryote_isa) , file = "sample_data_physeq_prokaryote_isa.csv")
#sample_data_physeq_prokaryote_isa <- read.csv("sample_data_physeq_prokaryote_isa.csv", header=T, row.names =1)
sample_data(physeq_prokaryote_isa) <- sample_data(sample_data_physeq_prokaryote_isa)
sample_data(physeq_prokaryote_isa)$Isa <- factor(sample_data(physeq_prokaryote_isa)$Isa,
                                                 levels=c("CapS1","CapS2","CapS3","CapS4","CapS5",
                                                          "CapS6","CapS7","CapS8","CapS9","CapS10",
                                                          "StemS1","StemS2","StemS3","StemS4","StemS5",
                                                          "StemS6","StemS7","StemS8","StemS9","StemS10", 
                                                          "SoilS1","SoilS2","SoilS3","SoilS4","SoilS5",
                                                          "SoilS6","SoilS7","SoilS8","SoilS9","SoilS10"))
levels(sample_data(physeq_prokaryote_isa)$Isa)
sample_data(physeq_prokaryote_isa)
otu_table(physeq_prokaryote_isa)

isa_obj <- as.data.frame(otu_table(physeq_prokaryote_isa))
dim(isa_obj)
isa_obj

# Creating a data.frame for plotting the heatmap 
identical(colnames(isa_obj), rownames(sample_data(physeq_prokaryote_isa)))
sample_data(physeq_prokaryote_isa)$Isa
colnames(isa_obj) <- sample_data(physeq_prokaryote_isa)$Isa
colnames(isa_df) <- c("Cap", "Soil", "Stem", "Index", "Stat", "p.value")
identical(rownames(isa_obj), rownames(isa_df))

isa_obj <- cbind(isa_obj, isa_df)
isa_obj$readNo <- rowSums(otu_prokaryote_out[rownames(isa_df),])
isa_obj$relAb <- (isa_obj$readNo/sum(colSums(otu_prokaryote_out))) * 100
isa_obj$logAb <- log(isa_obj$readNo)
isa_obj$sqrtAb <- sqrt(isa_obj$readNo)
identical(rownames(tax_table(physeq_prokaryote_isa)), rownames(isa_obj))
isa_obj$taxa <- paste(rownames(isa_obj), as.data.frame(as.matrix(tax_table(physeq_prokaryote_isa)))$GenBank, sep = " ")
isa_obj$GenBank <- as.data.frame(as.matrix(tax_table(physeq_prokaryote_isa)))$GenBank
isa_obj

# >>> HEATMAP -----------------------------------------------------------------------
library("ComplexHeatmap")
library("circlize")

ht1 = Heatmap(as.matrix(sqrt(isa_obj[,1:30]*10)), col = colorRamp2(c(0, 5), c("white","red")), 
              cluster_rows = FALSE, cluster_columns = TRUE, name = "Abundance",
              row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
              show_heatmap_legend = FALSE, row_labels = isa_obj$taxa)


ht2 = Heatmap(as.matrix(isa_obj[,31:33]), col = structure(c("red","white"), names = c("1","0")),
              cluster_rows = FALSE, cluster_columns = FALSE, name = "Group",
              row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
              show_heatmap_legend = FALSE, row_labels = isa_obj$taxa)

ha_bar = HeatmapAnnotation("Abundance (%)" = anno_barplot(isa_obj$sqrtAb/30, axis = TRUE, width = unit(2, "cm")), 
                           which = "row", show_annotation_name = TRUE, annotation_name_gp = gpar(fontsize = 8),
                           annotation_name_rot = c(0)) # annotation_width = unit(20, "cm"), annotation_height = unit(10, "cm"),

# *** FIGURE 3 - HEATMAP -------------------------------------------------------------------------
ht1 + ha_bar + ht2 

# >>> COMMUNITY COMPOSITION **---------------------------------------------------------------------

# fungal composition -----------------------------------------------------------------------------
physeq_fungi_phylum <- tax_glom(physeq_fungi, "Phylum")
otu_fungi_abund_phy <- taxa_sums(physeq_fungi_phylum)/sum(taxa_sums(physeq_fungi_phylum))*100
tax_abund_fungi_phy <- as(tax_table(physeq_fungi_phylum), "matrix")
tax_abund_fungi_phy <- as.data.frame(tax_abund_fungi_phy)
tax_abund_fungi_phy <- tax_abund_fungi_phy[c(2)]
tax_abund_fungi_phy$abundance <- as.vector(otu_fungi_abund_phy)
tax_abund_fungi_phy <- tax_abund_fungi_phy[order(tax_abund_fungi_phy$abundance, decreasing = TRUE),] 
tax_abund_fungi_phy

ggplot(tax_abund_fungi_phy, aes(x=reorder(Phylum, -abundance), y=abundance)) + 
  geom_bar(stat="identity") +
  labs(title="Soil", x="Taxa", y="Relative\nabundance %") +
  theme(axis.text.x = element_text(size=10, angle = -90, hjust = 0, vjust = 0.5))


physeq_fungi_order = tax_glom(physeq_fungi, "Order")
otu_fungi_abund_ord = taxa_sums(physeq_fungi_order)/sum(taxa_sums(physeq_fungi_order))*100
tax_abund_fungi_ord <- as(tax_table(physeq_fungi_order), "matrix")
tax_abund_fungi_ord <- as.data.frame(tax_abund_fungi_ord)
tax_abund_fungi_ord <- tax_abund_fungi_ord[c(1:4)]
tax_abund_fungi_ord$abundance <- as.vector(otu_fungi_abund_ord)
tax_abund_fungi_ord <- tax_abund_fungi_ord[order(tax_abund_fungi_ord$abundance, decreasing = TRUE),] 
tax_abund_fungi_ord

physeq_fungi_genus = tax_glom(physeq_fungi, "Genus")
otu_fungi_abund_gen = taxa_sums(physeq_fungi_genus)/sum(taxa_sums(physeq_fungi_genus))*100
tax_abund_fungi_gen <- as(tax_table(physeq_fungi_genus), "matrix")
tax_abund_fungi_gen <- as.data.frame(tax_abund_fungi_gen)
tax_abund_fungi_gen <- tax_abund_fungi_gen[c(1:6)]
tax_abund_fungi_gen$abundance <- as.vector(otu_fungi_abund_gen)
tax_abund_fungi_gen <- tax_abund_fungi_gen[order(tax_abund_fungi_gen$abundance, decreasing = TRUE),] 
tax_abund_fungi_gen

# extracting abundances for different stages
physeq_fungi_Young <- subset_samples(physeq_fungi, Stage%in%c("Young"))
otu_table(physeq_fungi_Young) <- otu_table(physeq_fungi_Young)[which(rowSums(otu_table(physeq_fungi_Young)) >= 1),] 
physeq_fungi_Young

physeq_fungi_Young_phylum = tax_glom(physeq_fungi_Young, "Phylum")
otu_fungi_abund_Young_phy = taxa_sums(physeq_fungi_Young_phylum)/sum(taxa_sums(physeq_fungi_Young_phylum))*100
tax_abund_fungi_Young_phy <- as(tax_table(physeq_fungi_Young_phylum), "matrix")
tax_abund_fungi_Young_phy <- as.data.frame(tax_abund_fungi_Young_phy)
tax_abund_fungi_Young_phy <- tax_abund_fungi_Young_phy[c(2)]
tax_abund_fungi_Young_phy$abundance <- as.vector(otu_fungi_abund_Young_phy)
tax_abund_fungi_Young_phy <- tax_abund_fungi_Young_phy[order(tax_abund_fungi_Young_phy$abundance, decreasing = TRUE),] 
tax_abund_fungi_Young_phy

physeq_fungi_Young_genus = tax_glom(physeq_fungi_Young, "Genus")
otu_fungi_abund_Young_gen = taxa_sums(physeq_fungi_Young_genus)/sum(taxa_sums(physeq_fungi_Young_genus))*100
tax_abund_fungi_Young_gen <- as(tax_table(physeq_fungi_Young_genus), "matrix")
tax_abund_fungi_Young_gen <- as.data.frame(tax_abund_fungi_Young_gen)
tax_abund_fungi_Young_gen <- tax_abund_fungi_Young_gen[c(1:6)]
tax_abund_fungi_Young_gen$abundance <- as.vector(otu_fungi_abund_Young_gen)
tax_abund_fungi_Young_gen <- tax_abund_fungi_Young_gen[order(tax_abund_fungi_Young_gen$abundance, decreasing = TRUE),] 
tax_abund_fungi_Young_gen[1:50, ]

physeq_fungi_Mature <- subset_samples(physeq_fungi, Stage%in%c("Mature"))
otu_table(physeq_fungi_Mature) <- otu_table(physeq_fungi_Mature)[which(rowSums(otu_table(physeq_fungi_Mature)) >= 1),] 
physeq_fungi_Mature

physeq_fungi_Mature_phylum = tax_glom(physeq_fungi_Mature, "Phylum")
otu_fungi_abund_Mature_phy = taxa_sums(physeq_fungi_Mature_phylum)/sum(taxa_sums(physeq_fungi_Mature_phylum))*100
tax_abund_fungi_Mature_phy <- as(tax_table(physeq_fungi_Mature_phylum), "matrix")
tax_abund_fungi_Mature_phy <- as.data.frame(tax_abund_fungi_Mature_phy)
tax_abund_fungi_Mature_phy <- tax_abund_fungi_Mature_phy[c(2)]
tax_abund_fungi_Mature_phy$abundance <- as.vector(otu_fungi_abund_Mature_phy)
tax_abund_fungi_Mature_phy <- tax_abund_fungi_Mature_phy[order(tax_abund_fungi_Mature_phy$abundance, decreasing = TRUE),] 
tax_abund_fungi_Mature_phy

physeq_fungi_Mature_genus = tax_glom(physeq_fungi_Mature, "Genus")
otu_fungi_abund_Mature_gen = taxa_sums(physeq_fungi_Mature_genus)/sum(taxa_sums(physeq_fungi_Mature_genus))*100
tax_abund_fungi_Mature_gen <- as(tax_table(physeq_fungi_Mature_genus), "matrix")
tax_abund_fungi_Mature_gen <- as.data.frame(tax_abund_fungi_Mature_gen)
tax_abund_fungi_Mature_gen <- tax_abund_fungi_Mature_gen[c(1:6)]
tax_abund_fungi_Mature_gen$abundance <- as.vector(otu_fungi_abund_Mature_gen)
tax_abund_fungi_Mature_gen <- tax_abund_fungi_Mature_gen[order(tax_abund_fungi_Mature_gen$abundance, decreasing = TRUE),] 
tax_abund_fungi_Mature_gen[1:50, ]

# prokaryotic composition ------------------------------------------------------------------------
physeq_prokaryote_phylum = tax_glom(physeq_prokaryote, "Phylum")
otu_prokaryote_abund_phy = taxa_sums(physeq_prokaryote_phylum)/sum(taxa_sums(physeq_prokaryote_phylum))*100
tax_abund_prokaryote_phy <- as(tax_table(physeq_prokaryote_phylum), "matrix")
tax_abund_prokaryote_phy <- as.data.frame(tax_abund_prokaryote_phy)
tax_abund_prokaryote_phy <- tax_abund_prokaryote_phy[c(2)]
tax_abund_prokaryote_phy$abundance <- as.vector(otu_prokaryote_abund_phy)
tax_abund_prokaryote_phy <- tax_abund_prokaryote_phy[order(tax_abund_prokaryote_phy$abundance, decreasing = TRUE),] 
tax_abund_prokaryote_phy

physeq_prokaryote_class = tax_glom(physeq_prokaryote, "Class")
otu_prokaryote_abund_class = taxa_sums(physeq_prokaryote_class)/sum(taxa_sums(physeq_prokaryote_class))*100
tax_prokaryote_abund_class <- as(tax_table(physeq_prokaryote_class), "matrix")
tax_prokaryote_abund_class <- as.data.frame(tax_prokaryote_abund_class)
tax_prokaryote_abund_class <- tax_prokaryote_abund_class[c(1:3)]
tax_prokaryote_abund_class$abundance <- as.vector(otu_prokaryote_abund_class)
tax_prokaryote_abund_class <- tax_prokaryote_abund_class[order(tax_prokaryote_abund_class$abundance, decreasing = TRUE),] 
tax_prokaryote_abund_class[1:50, ]

physeq_prokaryote_order = tax_glom(physeq_prokaryote, "Order")
otu_prokaryote_abund_ord = taxa_sums(physeq_prokaryote_order)/sum(taxa_sums(physeq_prokaryote_order))*100
tax_prokaryote_abund_ord <- as(tax_table(physeq_prokaryote_order), "matrix")
tax_prokaryote_abund_ord <- as.data.frame(tax_prokaryote_abund_ord)
tax_prokaryote_abund_ord <- tax_prokaryote_abund_ord[c(1:4)]
tax_prokaryote_abund_ord$abundance <- as.vector(otu_prokaryote_abund_ord)
tax_prokaryote_abund_ord <- tax_prokaryote_abund_ord[order(tax_prokaryote_abund_ord$abundance, decreasing = TRUE),] 
tax_prokaryote_abund_ord[1:50, ]

physeq_prokaryote_genus = tax_glom(physeq_prokaryote, "Genus")
otu_prokaryote_abund_gen = taxa_sums(physeq_prokaryote_genus)/sum(taxa_sums(physeq_prokaryote_genus))*100
tax_abund_prokaryote_gen <- as(tax_table(physeq_prokaryote_genus), "matrix")
tax_abund_prokaryote_gen <- as.data.frame(tax_abund_prokaryote_gen)
tax_abund_prokaryote_gen <- tax_abund_prokaryote_gen[c(1:6)]
tax_abund_prokaryote_gen$abundance <- as.vector(otu_prokaryote_abund_gen)
tax_abund_prokaryote_gen <- tax_abund_prokaryote_gen[order(tax_abund_prokaryote_gen$abundance, decreasing = TRUE),] 
tax_abund_prokaryote_gen[1:50, ]

# extracting abundances for different origins
physeq_prokaryote_Cap <- subset_samples(physeq_prokaryote, Origin%in%c("Cap"))
otu_table(physeq_prokaryote_Cap) <- otu_table(physeq_prokaryote_Cap)[which(rowSums(otu_table(physeq_prokaryote_Cap)) >= 1),] 
physeq_prokaryote_Cap

physeq_prokaryote_Cap_phylum = tax_glom(physeq_prokaryote_Cap, "Phylum")
otu_prokaryote_abund_Cap_phy = taxa_sums(physeq_prokaryote_Cap_phylum)/sum(taxa_sums(physeq_prokaryote_Cap_phylum))*100
tax_abund_prokaryote_Cap_phy <- as(tax_table(physeq_prokaryote_Cap_phylum), "matrix")
tax_abund_prokaryote_Cap_phy <- as.data.frame(tax_abund_prokaryote_Cap_phy)
tax_abund_prokaryote_Cap_phy <- tax_abund_prokaryote_Cap_phy[c(2)]
tax_abund_prokaryote_Cap_phy$abundance <- as.vector(otu_prokaryote_abund_Cap_phy)
tax_abund_prokaryote_Cap_phy <- tax_abund_prokaryote_Cap_phy[order(tax_abund_prokaryote_Cap_phy$abundance, decreasing = TRUE),] 
tax_abund_prokaryote_Cap_phy

physeq_prokaryote_Cap_genus = tax_glom(physeq_prokaryote_Cap, "Genus")
otu_prokaryote_abund_Cap_gen = taxa_sums(physeq_prokaryote_Cap_genus)/sum(taxa_sums(physeq_prokaryote_Cap_genus))*100
tax_abund_prokaryote_Cap_gen <- as(tax_table(physeq_prokaryote_Cap_genus), "matrix")
tax_abund_prokaryote_Cap_gen <- as.data.frame(tax_abund_prokaryote_Cap_gen)
tax_abund_prokaryote_Cap_gen <- tax_abund_prokaryote_Cap_gen[c(1:6)]
tax_abund_prokaryote_Cap_gen$abundance <- as.vector(otu_prokaryote_abund_Cap_gen)
tax_abund_prokaryote_Cap_gen <- tax_abund_prokaryote_Cap_gen[order(tax_abund_prokaryote_Cap_gen$abundance, decreasing = TRUE),] 
tax_abund_prokaryote_Cap_gen[1:50, ]

physeq_prokaryote_Stem <- subset_samples(physeq_prokaryote, Origin%in%c("Stem"))
otu_table(physeq_prokaryote_Stem) <- otu_table(physeq_prokaryote_Stem)[which(rowSums(otu_table(physeq_prokaryote_Stem)) >= 1),] 
physeq_prokaryote_Stem

physeq_prokaryote_Stem_phylum = tax_glom(physeq_prokaryote_Stem, "Phylum")
otu_prokaryote_abund_Stem_phy = taxa_sums(physeq_prokaryote_Stem_phylum)/sum(taxa_sums(physeq_prokaryote_Stem_phylum))*100
tax_abund_prokaryote_Stem_phy <- as(tax_table(physeq_prokaryote_Stem_phylum), "matrix")
tax_abund_prokaryote_Stem_phy <- as.data.frame(tax_abund_prokaryote_Stem_phy)
tax_abund_prokaryote_Stem_phy <- tax_abund_prokaryote_Stem_phy[c(2)]
tax_abund_prokaryote_Stem_phy$abundance <- as.vector(otu_prokaryote_abund_Stem_phy)
tax_abund_prokaryote_Stem_phy <- tax_abund_prokaryote_Stem_phy[order(tax_abund_prokaryote_Stem_phy$abundance, decreasing = TRUE),] 
tax_abund_prokaryote_Stem_phy

physeq_prokaryote_Stem_genus = tax_glom(physeq_prokaryote_Stem, "Genus")
otu_prokaryote_abund_Stem_gen = taxa_sums(physeq_prokaryote_Stem_genus)/sum(taxa_sums(physeq_prokaryote_Stem_genus))*100
tax_abund_prokaryote_Stem_gen <- as(tax_table(physeq_prokaryote_Stem_genus), "matrix")
tax_abund_prokaryote_Stem_gen <- as.data.frame(tax_abund_prokaryote_Stem_gen)
tax_abund_prokaryote_Stem_gen <- tax_abund_prokaryote_Stem_gen[c(1:6)]
tax_abund_prokaryote_Stem_gen$abundance <- as.vector(otu_prokaryote_abund_Stem_gen)
tax_abund_prokaryote_Stem_gen <- tax_abund_prokaryote_Stem_gen[order(tax_abund_prokaryote_Stem_gen$abundance, decreasing = TRUE),] 
tax_abund_prokaryote_Stem_gen[1:50, ]

physeq_prokaryote_Soil <- subset_samples(physeq_prokaryote, Origin%in%c("Soil"))
otu_table(physeq_prokaryote_Soil) <- otu_table(physeq_prokaryote_Soil)[which(rowSums(otu_table(physeq_prokaryote_Soil)) >= 1),] 
physeq_prokaryote_Soil

physeq_prokaryote_Soil_phylum = tax_glom(physeq_prokaryote_Soil, "Phylum")
otu_prokaryote_abund_Soil_phy = taxa_sums(physeq_prokaryote_Soil_phylum)/sum(taxa_sums(physeq_prokaryote_Soil_phylum))*100
tax_abund_prokaryote_Soil_phy <- as(tax_table(physeq_prokaryote_Soil_phylum), "matrix")
tax_abund_prokaryote_Soil_phy <- as.data.frame(tax_abund_prokaryote_Soil_phy)
tax_abund_prokaryote_Soil_phy <- tax_abund_prokaryote_Soil_phy[c(2)]
tax_abund_prokaryote_Soil_phy$abundance <- as.vector(otu_prokaryote_abund_Soil_phy)
tax_abund_prokaryote_Soil_phy <- tax_abund_prokaryote_Soil_phy[order(tax_abund_prokaryote_Soil_phy$abundance, decreasing = TRUE),] 
tax_abund_prokaryote_Soil_phy

physeq_prokaryote_Soil_genus = tax_glom(physeq_prokaryote_Soil, "Genus")
otu_prokaryote_abund_Soil_gen = taxa_sums(physeq_prokaryote_Soil_genus)/sum(taxa_sums(physeq_prokaryote_Soil_genus))*100
tax_abund_prokaryote_Soil_gen <- as(tax_table(physeq_prokaryote_Soil_genus), "matrix")
tax_abund_prokaryote_Soil_gen <- as.data.frame(tax_abund_prokaryote_Soil_gen)
tax_abund_prokaryote_Soil_gen <- tax_abund_prokaryote_Soil_gen[c(1:6)]
tax_abund_prokaryote_Soil_gen$abundance <- as.vector(otu_prokaryote_abund_Soil_gen)
tax_abund_prokaryote_Soil_gen <- tax_abund_prokaryote_Soil_gen[order(tax_abund_prokaryote_Soil_gen$abundance, decreasing = TRUE),] 
tax_abund_prokaryote_Soil_gen[1:50, ]

# looking for specific genera
tax_abund_prokaryote_Soil_gen[tax_abund_prokaryote_Soil_gen$Genus=="Pedobacter",]

# plottign bars ----------------------------------------------------------------------------------
library("scales")
library("grid")
library("reshape2")

# melt to long format (for ggploting) ------------------------------------------------------------
# and  prune out phyla below 2% in each sample
otu_fungi_phylum <- physeq_fungi %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

otu_fungi_phylum
head(otu_fungi_phylum)
unique(otu_fungi_phylum$Family)

# Set colors for plotting
phylum_colors <- c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
                   "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5",
                   "#5E738F","#D1A33D", "#8A7C64", "#599861")

# *** FIGURE 1A barplot FUNGI 
barplot_fungi = ggplot(otu_fungi_phylum, aes(x = Sample, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity") +
  labs(title="", x="Samples", y = "Relative Abundance") +
  scale_fill_manual(values = phylum_colors) +
  facet_grid(~Stage, scales = "free_x", space="free_x") +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
  theme(strip.text.x = element_text(size = 8, face = "bold")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")

barplot_fungi

# composition Prokaryotic communities ------------------------------------------------------------
sample_data(physeq_prokaryote)$Origin <- factor(sample_data(physeq_prokaryote)$Origin,
                                                levels=c("Cap","Stem","Soil"))

otu_prokaryote_phylum <- physeq_prokaryote %>%
  tax_glom(taxrank = "Phylum") %>%                     
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.01) %>%                         
  arrange(Phylum)                                     

otu_prokaryote_phylum
head(otu_prokaryote_phylum)
unique(otu_prokaryote_phylum$Phylum)

# *** FIGURE 1B barplot prokaryote 
barplot_prokaryote = ggplot(otu_prokaryote_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  labs(title="", x="Samples", y = "Relative Abundance") +
  scale_fill_manual(values = phylum_colors) +
  facet_grid(~Origin, scales = "free_x", space="free_x") +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
  theme(strip.text.x = element_text(size = 8, face = "bold")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")

barplot_prokaryote

# *** FIGURE 1 - BARPLOTS ------------------------------------------------------------------------

ggarrange(barplot_fungi,
          barplot_prokaryote,
          labels = c("A", "B"),
          widths = c(1.05,2),
          align = "none", 
          ncol = 2, nrow = 1,
          common.legend = FALSE,
          legend = c("right"))

# >>> VENN DIAGRAM -------------------------------------------------------------------------------
library("limma")
source("../my_venn_diag.R")

physeq_fungi_St = merge_samples(physeq_fungi, "Stage")
otu_fungi_St <- as.data.frame(t(otu_table(physeq_fungi_St)))
venn_counts_otu_fungi_St <- vennCounts(otu_fungi_St, include="both")
venn_counts_otu_fungi_St

physeq_prokaryote_St = merge_samples(physeq_prokaryote, "Stage")
otu_prokaryote_St <- as.data.frame(t(otu_table(physeq_prokaryote_St)))
venn_counts_otu_prokaryote_St <- vennCounts(otu_prokaryote_St, include="both")
venn_counts_otu_prokaryote_St

physeq_prokaryote_Ts = merge_samples(physeq_prokaryote, "Origin")
otu_prokaryote_Ts <- as.data.frame(t(otu_table(physeq_prokaryote_Ts)))
venn_counts_otu_prokaryote_Ts <- vennCounts(otu_prokaryote_Ts, include="both")
venn_counts_otu_prokaryote_Ts

# *** FIGURE 4 - venn diagrams -------------------------------------------------------------------

layout(matrix(1:3, ncol=3))
venn_Gian(venn_counts_otu_prokaryote_Ts,
          cex=c(1),
          circle.col =c("#bdbdbd", "#525252", "#000000"),
          mar = c(1,1,1,1),
          lwd = 2, main="A")

venn_Gian(venn_counts_otu_prokaryote_St,
          cex=c(1),
          circle.col =c("red", "grey"),
          mar = c(1,1,1,1),
          lwd = 2, main="B")

venn_Gian(venn_counts_otu_fungi_St,
          cex=c(1),
          circle.col =c("red", "grey"),
          mar = c(1,1,1,1),
          lwd = 2, main="C")
dev.off()


# >>> MICROBIAL NETWORK - prokaryotic communitie -----------------------------------------------------

# Note: This is just an example fo the code used to plot the 
# networks included in the manuscript. For this reason, the network
# figure (the plot) will not be the exact same of the one included 
# in the manuscript.

# extarcting taxa present in >= 17 samples. This value is choosen 
# according to: 1) group of samples present in the variable with the most 
# sample per treatment, 2) by visually estimating the best level of 
# sparsity that will make clear to distinguish  important network
# properties, and 3) by an optimal network stability.

physeq_prokaryote -> physeq_prokaryote_filt
cntNonZero <- apply(as.data.frame(otu_table(physeq_prokaryote)), 1, function(x) sum(x > 0))
otu_table(physeq_prokaryote_filt) <- otu_table(physeq_prokaryote_filt)[which(cntNonZero >=17),]
rowSums(otu_table(physeq_prokaryote_filt))
physeq_prokaryote_filt

# Network analysis -------------------------------------------------------------------------------
library("SpiecEasi"); packageVersion("SpiecEasi")
library("igraph"); packageVersion("igraph")
library("huge"); packageVersion("huge")
library("qgraph"); packageVersion("qgraph")
library("MCL")

set.seed(0619)

spiec_prok <- spiec.easi(physeq_prokaryote_filt, 
                         method="mb",
                         lambda.min.ratio=1e-2, 
                         nlambda=50, 
                         sel.criterion ="stars", 
                         pulsar.select = TRUE,
                         pulsar.params=list(rep.num=99))

spiec_prok
getStability(spiec_prok)

# creating the network object 
plot_spiec_prok <- adj2igraph(getRefit(spiec_prok),
                              vertex.attr=list(name=taxa_names(physeq_prokaryote_filt)))

plot_spiec_prok


# PLOT function ----------------------------------------------------------------------------------
# modified from:
# Beiko, R. G., Hsiao, W., and Parkinson, J.  
# Microbiome Analysis: Methods and Protocols.
# Humana Press, 2018.

plot.net.cls <- function(net, scores, cls, art_point, physeq_obj) {
  # Get size of clusters to find isolated nodes.
  cls_sizes <- sapply(groups(cls), length)
  # Randomly choosing node colors. Users can provide their own vector of colors.
  colors <- sample(colours(), length(cls))
  # Nodes in clusters will be color coded. Isolated nodes will be white.
  V(net)$color <- sapply(membership(cls),
                         function(x) {ifelse(cls_sizes[x]>1,
                                             colors[x], "white")})
  # Convert node label from names to numerical IDs.
  node.names <- V(net)$name
  #col_ids <- seq(1, length(node.names))
  #V(net)$name <- col_ids
  # name nodes with taxonomic names 
  #taxa_physeq_obj <- as.data.frame(as.matrix(tax_table(physeq_obj)))
  #V(net)$name <- taxa_physeq_obj$Phylum
  # To draw a halo around articulation points.
  art_point <- lapply(names(art_point), function(x) x)
  marks <- lapply(1:length(art_point), function(x) which(node.names == art_point[[x]]))
  # vertex size 
  sum_list <- taxa_sums(physeq_obj)
  vsize = log(sum_list)/2
  # set size of vertex proportional to clr-mean
  #vsize <- rowMeans(clr(t(otu_table(physeq_obj)), 1)) + 7
  # Customized layout to avoid nodes overlapping.
  #e <- get.edgelist(net)
  #class(e) <- "numeric"
  #l <- qgraph.layout.fruchtermanreingold(e, vcount=vcount(net),
  #                                         area=8*(vcount(net)^2),
  #                                        repulse.rad=(vcount(net)^3.1))
  # Main plot function.
  plot.igraph(net, 
              vertex.size = vsize, 
              vertex.label.cex=0.4,
              vertex.label.color = "black",
              mark.border="black",
              mark.groups = marks,
              mark.col = "white",
              mark.expand = 10,
              mark.shape = 1,
              layout = NULL)
}

# *** FIGURE 5 - the network ---------------------------------------------------------------------
plot.net.cls(plot_spiec_prok, 
             hub_score(plot_spiec_prok)$vector, 
             walktrap.community(plot_spiec_prok), 
             articulation.points(plot_spiec_prok), 
             physeq_prokaryote_filt)

# >>> NETOWRK TOPOLOGY INDEXES ------------------------------------------------------------------
cent_res <- igraph::centrality(plot_spiec_prok, all.shortest.paths = TRUE)
cent_res$OutDegree
cent_res$InDegree
cent_res$Closeness
cent_res$Betweenness

# Degree 
igraph::degree(plot_spiec_prok, mode="all")

# Degree distribution 
deg.dist1 <- degree_distribution(plot_spiec_prok, mode = "all")
deg.dist1
plot(deg.dist1, type='b', ylab="Frequency", xlab="Degree", main="Degree Distribution")

# Transitivity or clustering coefficient
clustering_coeff <- transitivity(plot_spiec_prok, type = "global")
clustering_coeff

                  
# _______________ IMPORTANT NOTE ___________________-----------------

 # To relabel the Origin factor with the correct scientific terminology 
 # we used the code below and then we re-sun the R scripts 
sample_data(physeq_prokaryote_mSeq)$Origin <- c(rep("Pileus", 10),rep("Soil",10), rep("Stipe",10)) 
sample_data(physeq_prokaryote_mSeq)$Origin <- factor(sample_data(physeq_prokaryote_mSeq)$Origin,levels=c("Pileus","Stipe","Soil"))
                  
                  
                  
                  
                  
