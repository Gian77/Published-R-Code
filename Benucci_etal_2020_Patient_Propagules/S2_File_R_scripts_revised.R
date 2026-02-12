# Patient propagules: do long-term soil archives preserve legacy fungal and bacterial communities?
# Benucci GMN, Rennick B, Bonito G.
# Michigan State University
# PloSONE journal

# IF YOU USE THIS CODE, PLEASE CITE THE RELATED PAPER. YOU CAN FIND MORE
# DETAILS ON GITHUB AT THIS LINK https://github.com/Gian77/Scientific-Papers-R-Code

# Benucci et al. - manuscripot R scripts -----------------------------------------------------------
# WORKING ENVIRONMENT SETUP ------------------------------------------------------------------------
options(scipen = 999) 
options(max.print=100000000)

#rm(list = ls(all=TRUE)) # removes all variables in the global environment so you start fresh
Sys.time() # prints out the time and date you ran the code
set.seed(1977) #to make reproduceble results
detach(package:phyloseq, unload=TRUE) #to reorder packages
search() #to search into the environment
#rm(list= ls()[!(ls() %in% c('keepThis','andThis'))]) #remove all but ...
session.info() #to get session information

# NECESSARY PACKAGES -------------------------------------------------------------------------------

library(phyloseq)
library(Biostrings)
library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(ggpubr)

# color palettes -----------------------------------------------------------------------------------
palette_CB6 <-c("#332288","#88CCEE","#117733","#DDCC77","#FF899d","#AA4499")
paletteCB3 = c("#332288","#88CCEE","#D55E00")


paletteCB3_fungi = c("#332288","#88CCEE","#117733")
paletteCB3_prok = c("#E69F00","#FF899d","#663F05")

paletteCB3 = c("#332288","#88CCEE","#D55E00")


# importing datasets -------------------------------------------------------------------------------
biom_ITS_uparse = import_biom("Suppl_file_S4_otutab_UPARSE_ITS_tax_json.biom")
map_ITS = import_qiime_sample_data("Suppl_file_S6_metadata_its.txt")
sample_data(biom_ITS_uparse) <- map_ITS
colnames(tax_table(biom_ITS_uparse)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
otus_rep_ITS_uparse <- readDNAStringSet("otus_ITS.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
biom_ITS_uparse <- merge_phyloseq(biom_ITS_uparse, otus_rep_ITS_uparse)
head(tax_table(biom_ITS_uparse))
refseq(biom_ITS_uparse)
biom_ITS_uparse

biom_16s_uparse = import_biom("Suppl_file_S5_otutab_UPARSE_16s_tax_json.biom")
map_16s = import_qiime_sample_data("Suppl_file_S7_metadata_16s.txt")
sample_data(biom_16s_uparse) <- map_16s
colnames(tax_table(biom_16s_uparse)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
otus_rep_16s_uparse <- readDNAStringSet("otus_16s.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
biom_16s_uparse <- merge_phyloseq(biom_16s_uparse, otus_rep_16s_uparse)
sample_data(biom_16s_uparse)
head(tax_table(biom_16s_uparse))
biom_16s_uparse

# cleaning and filtering ---------------------------------------------------------------------------
tax_table(biom_ITS_uparse)[, "Kingdom"] <- gsub("PMI", "", tax_table(biom_ITS_uparse)[, "Kingdom"]) # based on the reference taxonomy
tax_table(biom_ITS_uparse)[, "Kingdom"] <- gsub("NVP", "", tax_table(biom_ITS_uparse)[, "Kingdom"])

tax_table(biom_ITS_uparse)[, "Kingdom"] <- gsub("k__", "", tax_table(biom_ITS_uparse)[, "Kingdom"])
tax_table(biom_ITS_uparse)[, "Phylum"] <- gsub("p__", "", tax_table(biom_ITS_uparse)[, "Phylum"])
tax_table(biom_ITS_uparse)[, "Class"] <- gsub("c__", "", tax_table(biom_ITS_uparse)[, "Class"])
tax_table(biom_ITS_uparse)[, "Order"] <- gsub("o__", "", tax_table(biom_ITS_uparse)[, "Order"])
tax_table(biom_ITS_uparse)[, "Family"] <- gsub("f__", "", tax_table(biom_ITS_uparse)[, "Family"])
tax_table(biom_ITS_uparse)[, "Genus"] <- gsub("g__", "", tax_table(biom_ITS_uparse)[, "Genus"])
tax_table(biom_ITS_uparse)[, "Species"] <- gsub("s__", "", tax_table(biom_ITS_uparse)[, "Species"])

biom_ITS_uparse <- subset_taxa(biom_ITS_uparse, Kingdom == "Fungi")
tax_table(biom_ITS_uparse)[tax_table(biom_ITS_uparse)=="unidentified"]<- NA
tax_table(biom_ITS_uparse)[is.na(tax_table(biom_ITS_uparse))]<-"Unclassified"

head(otu_table(biom_ITS_uparse))
head(tax_table(biom_ITS_uparse))
head(sample_data(biom_ITS_uparse))
biom_ITS_uparse

tax_table(biom_16s_uparse)[, "Kingdom"] <- gsub("k__", "", tax_table(biom_16s_uparse)[, "Kingdom"])
tax_table(biom_16s_uparse)[, "Phylum"] <- gsub("p__", "", tax_table(biom_16s_uparse)[, "Phylum"])
tax_table(biom_16s_uparse)[, "Class"] <- gsub("c__", "", tax_table(biom_16s_uparse)[, "Class"])
tax_table(biom_16s_uparse)[, "Order"] <- gsub("o__", "", tax_table(biom_16s_uparse)[, "Order"])
tax_table(biom_16s_uparse)[, "Family"] <- gsub("f__", "", tax_table(biom_16s_uparse)[, "Family"])
tax_table(biom_16s_uparse)[, "Genus"] <- gsub("g__", "", tax_table(biom_16s_uparse)[, "Genus"])
tax_table(biom_16s_uparse)[, "Species"] <- gsub("s__", "", tax_table(biom_16s_uparse)[, "Species"])

biom_16s_uparse <- subset_taxa(biom_16s_uparse, Phylum!="Chloroplast")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Class!="Chloroplast")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Order!="Chloroplast")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Family!="Chloroplast")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Genus!="Chloroplast")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Phylum!="chloroplast")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Class!="chloroplast")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Order!="chloroplast")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Family!="chloroplast")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Genus!="chloroplast")

biom_16s_uparse <- subset_taxa(biom_16s_uparse, Phylum!="Mitochondria")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Class!="Mitochondria")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Order!="Mitochondria")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Family!="Mitochondria")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Genus!="Mitochondria")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Phylum!="mitochondria")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Class!="mitochondria")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Order!="mitochondria")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Family!="mitochondria")
biom_16s_uparse <- subset_taxa(biom_16s_uparse, Genus!="mitochondria")

tax_table(biom_16s_uparse)[tax_table(biom_16s_uparse)==""]<- NA
tax_table(biom_16s_uparse)[is.na(tax_table(biom_16s_uparse))]<-"Unclassified"
head(otu_table(biom_16s_uparse))
head(tax_table(biom_16s_uparse))
head(sample_data(biom_16s_uparse))
biom_16s_uparse

# function for filtering contaminants --------------------------------------------------------------
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

# taxa to filter out 
badTaxa_uparse_ITS =c("OTU_8","OTU_39","OTU_269")

biom_ITS_uparse_qf_prun = pop_taxa(biom_ITS_uparse, badTaxa_uparse_ITS)
biom_ITS_uparse_qf_prun

badTaxa_uparse_16s =c("OTU_4","OTU_9","OTU_8176","OTU_5073","OTU_37","OTU_47","OTU_109","OTU_135","OTU_133",
                      "OTU_405","OTU_603","OTU_397","OTU_865","OTU_1303","OTU_970","OTU_2480",
                      "OTU_2015","OTU_2488","OTU_2642","OTU_2815","OTU_2390","OTU_12514",
                      "OTU_9394","OTU_7279","OTU_8176")

biom_16s_uparse_qf_prun = pop_taxa(biom_16s_uparse, badTaxa_uparse_16s)
biom_16s_uparse_qf_prun

# Filtering out OTUs -------------------------------------------------------------------------------
# Oliver et al. 2015, PCR errors and tag switching
# Lindhal et al. 2013, tag switching - that's a good  one!

biom_ITS_uparse_qf_prun -> biom_ITS_uparse_filt
otu_table(biom_ITS_uparse_filt) <- otu_table(biom_ITS_uparse_filt)[which(rowSums(otu_table(biom_ITS_uparse_filt)) >= 10),]
biom_ITS_uparse_filt

biom_16s_uparse_qf_prun -> biom_16s_uparse_filt
otu_table(biom_16s_uparse_filt) <- otu_table(biom_16s_uparse_filt)[which(rowSums(otu_table(biom_16s_uparse_filt)) >= 10),] 
biom_16s_uparse_filt

# filtering out PCR sequencing control samples -----------------------------------------------------

sample_data(biom_ITS_uparse_filt)
otu_table(biom_ITS_uparse_filt) <- subset(otu_table(biom_ITS_uparse_filt), select = -c(DF33bis, control2))
otu_table(biom_ITS_uparse_filt) <- otu_table(biom_ITS_uparse_filt)[which(rowSums(otu_table(biom_ITS_uparse_filt)) >= 1),] 
biom_ITS_uparse_filt

sum(rowSums(otu_table(biom_ITS_uparse_filt)))
any(taxa_sums(biom_ITS_uparse_filt) == 0)

sample_data(biom_16s_uparse_filt)
otu_table(biom_16s_uparse_filt) <- subset(otu_table(biom_16s_uparse_filt),
                          select = -c(DF33bis, FJ34B, P615A, C4FC, p589, empty1, empty3, control3, control2, control1))
otu_table(biom_16s_uparse_filt) <- otu_table(biom_16s_uparse_filt)[which(rowSums(otu_table(biom_16s_uparse_filt)) >= 1),] 
biom_16s_uparse_filt

sum(rowSums(otu_table(biom_16s_uparse_filt)))
any(taxa_sums(biom_16s_uparse_filt) == 0)

# RAREFACTION CURVES -------------------------------------------------------------------------------
library(ggpubr)

get_palette(palette = "jco", 3)
curve_colors_fungi <- rep("grey", ncol(as.data.frame(otu_table(biom_ITS_uparse_filt))))
curve_colors_fungi[as.data.frame(as.matrix(sample_data(biom_ITS_uparse_filt)))$exp1=="exp1"] <- "#0073C2FF"
curve_colors_fungi[as.data.frame(as.matrix(sample_data(biom_ITS_uparse_filt)))$exp2=="exp2"] <- "#EFC000FF"
curve_colors_fungi[as.data.frame(as.matrix(sample_data(biom_ITS_uparse_filt)))$exp3=="exp3"] <- "#868686FF"
curve_colors_fungi

curve_colors_prok <- rep("grey", ncol(as.data.frame(otu_table(biom_16s_uparse_filt))))
curve_colors_prok[as.data.frame(as.matrix(sample_data(biom_16s_uparse_filt)))$exp1=="exp1"] <- "#0073C2FF"
curve_colors_prok[as.data.frame(as.matrix(sample_data(biom_16s_uparse_filt)))$exp2=="exp2"] <- "#EFC000FF"
curve_colors_prok[as.data.frame(as.matrix(sample_data(biom_16s_uparse_filt)))$exp3=="exp3"] <- "#868686FF"
curve_colors_prok

# *** FIGURE S4 ******------------------------------------------------------------------------------
par(mfrow=c(1,2)) # Change the panel layout to 2 x 1
rarecurve(t(as.data.frame(otu_table(biom_ITS_uparse_filt))), 
          col = curve_colors_fungi, 
          label = FALSE, step = 500, 
          main="Fungi", ylab = "Number of OTUs", xlab = "Number of DNA reads") -> rare_fungi
legend("bottomright", legend=c("Exp.1","Exp.2","Exp.3"),
       col=c("#0073C2FF", "#EFC000FF", "#868686FF"), lty=1, cex=0.8, bty = "n") 
abline(v = min(sample_sums(biom_ITS_uparse_filt)), col="red", lwd=3, lty=2)

rarecurve(t(as.data.frame(otu_table(biom_16s_uparse_filt))), 
          col = curve_colors_prok, 
          label = FALSE, step = 500, 
          main="Prokaryotes", ylab = "Number of OTUs", xlab = "Number of DNA reads") -> rare_prok
legend("bottomright", legend=c("Exp.1","Exp.2","Exp.3"),
       col=c("#0073C2FF", "#EFC000FF", "#868686FF"), lty=1, cex=0.8,  bty = "n") 
abline(v = min(sample_sums(biom_16s_uparse_filt)), col="red", lwd=3, lty=2)
dev.off()

# INSEPCTIN LIBRARY SIZES --------------------------------------------------------------------------
sample_data(biom_ITS_uparse_filt)$LibrarySize <- sample_sums(biom_ITS_uparse_filt)
sample_data(biom_ITS_uparse_filt)$richness <- specnumber(as.data.frame(otu_table(biom_ITS_uparse_filt)) , MARGIN = 2)
sample_data(biom_ITS_uparse_filt)$rarefied <- rarefy(as.data.frame(otu_table(biom_ITS_uparse_filt)), 
                                                     min(sample_sums(biom_ITS_uparse_filt)), se = FALSE, MARGIN = 2)
sample_data(biom_ITS_uparse_filt)$shannon <- estimate_richness(biom_ITS_uparse_filt, split = TRUE, measures = "Shannon")$Shannon
sample_data(biom_ITS_uparse_filt)
str(sample_data(biom_ITS_uparse_filt))

# extracting sample names
unique(sample_data(biom_ITS_uparse_filt)$X.SampleID)

sample_data(biom_16s_uparse_filt)$LibrarySize <- sample_sums(biom_16s_uparse_filt)
sample_data(biom_16s_uparse_filt)$richness <- specnumber(as.data.frame(otu_table(biom_16s_uparse_filt)) , MARGIN = 2)
sample_data(biom_16s_uparse_filt)$rarefied <- rarefy(as.data.frame(otu_table(biom_16s_uparse_filt)), 
                                                     min(sample_sums(biom_16s_uparse_filt)), se = FALSE, MARGIN = 2)
sample_data(biom_16s_uparse_filt)$shannon <- estimate_richness(biom_16s_uparse_filt, split = TRUE, measures = "Shannon")$Shannon
sample_data(biom_16s_uparse_filt)

PlotDistrib <-function(dataframe){
  distrplot <- ggplot(sample_data(dataframe), aes(x = LibrarySize)) + 
  geom_histogram(color = "grey", fill = "grey", binwidth = 500) +
  labs(x="Read number", y="Sample number") + 
  geom_vline(aes(xintercept=mean(LibrarySize, na.rm=T)), color="black", linetype="dashed", size=0.8) + 
    theme_classic() +
    xlim(0, NA)+
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) 
return(distrplot)
}

# Figure SX Library distribution -------------------------------------------------------------------
ggarrange(PlotDistrib(biom_ITS_uparse_filt) + labs(title = "Fungi Soil\nDistribution of Sample Libraries"), 
          PlotDistrib(biom_16s_uparse_filt) + labs(title = "Prokaryotes Soil\nDistribution of Sample Libraries"),
          labels = c("A", "B"),
          widths = c(1, 1),
          align = "h", 
          ncol = 2, nrow = 1) -> figure_S1

figure_S1

ggarrange(PlotDistrib(biom_ITS_uparse_filt) + labs(title = "Fungi Soil\nDistribution of Sample Libraries")+
            facet_grid(~year~habitat), 
          PlotDistrib(biom_16s_uparse_filt) + labs(title = "Prokaryotes Soil\nDistribution of Sample Libraries")+
            facet_grid(~year~habitat),
          labels = c("A", "B"),
          widths = c(1, 1),
          align = "h", 
          ncol = 2, nrow = 1) -> figure_S2

figure_S2


# Extracting data of different experiments ---------------------------------------------------------
# 1) this is EXPERIMENT 1 ITS dataset --------------------------------------------------------------
library(metagenomeSeq)

# MetagenomeSeq function ---------------------------------------------------------------------------
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

# Extracting datasets 
biom_ITS_exp1 <- subset_samples(biom_ITS_uparse_filt, exp1%in%c("exp1"))
otu_table(biom_ITS_exp1) <- otu_table(biom_ITS_exp1)[which(rowSums(otu_table(biom_ITS_exp1)) >= 1),] 
biom_ITS_exp1
sample_data(biom_ITS_exp1)

# rarefaction at minimum library size --------------------------------------------------------------
set.seed(2020)
biom_ITS_exp1_ev = rarefy_even_depth(biom_ITS_exp1,
                   sample.size = min(sample_sums(biom_ITS_exp1)),
                   rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
biom_ITS_exp1_ev
sample_sums(biom_ITS_exp1_ev)

# MetagenomeSeq normalization - Gaussian model -----------------------------------------------------
CSSNorm(biom_ITS_exp1) -> biom_ITS_exp1_mSeq
head(otu_table(biom_ITS_exp1_mSeq))

# 2) this is EXPERIMENT 2 ITS dataset --------------------------------------------------------------
biom_ITS_exp2 <- subset_samples(biom_ITS_uparse_filt, exp2%in%c("exp2"))
otu_table(biom_ITS_exp2) <- otu_table(biom_ITS_exp2)[which(rowSums(otu_table(biom_ITS_exp2)) >=1),] 
biom_ITS_exp2
sample_data(biom_ITS_exp2)

set.seed(2022)
biom_ITS_exp2_ev = rarefy_even_depth(biom_ITS_exp2,
                   sample.size = min(sample_sums(biom_ITS_exp2)),
                   rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
biom_ITS_exp2_ev
sample_sums(biom_ITS_exp2_ev)

CSSNorm(biom_ITS_exp2) -> biom_ITS_exp2_mSeq
head(otu_table(biom_ITS_exp2_mSeq))

# 3) this is EXPERIMENT 3 ITS dataset --------------------------------------------------------------
biom_ITS_exp3 <- subset_samples(biom_ITS_uparse_filt, exp3%in%c("exp3"))
otu_table(biom_ITS_exp3) <- otu_table(biom_ITS_exp3)[which(rowSums(otu_table(biom_ITS_exp3)) >=1),] 
biom_ITS_exp3
sample_data(biom_ITS_exp3)

set.seed(2023)
biom_ITS_exp3_ev = rarefy_even_depth(biom_ITS_exp3,
                  sample.size = min(sample_sums(biom_ITS_exp3)),
                  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
biom_ITS_exp3_ev
sample_sums(biom_ITS_exp3_ev)

CSSNorm(biom_ITS_exp3) -> biom_ITS_exp3_mSeq
head(otu_table(biom_ITS_exp3_mSeq))

# 1) this is EXPERIMENT 1 16s dataset --------------------------------------------------------------
biom_16s_exp1 <- subset_samples(biom_16s_uparse_filt, exp1%in%c("exp1"))
otu_table(biom_16s_exp1) <- otu_table(biom_16s_exp1)[which(rowSums(otu_table(biom_16s_exp1)) >= 1),] 
biom_16s_exp1
sample_data(biom_16s_exp1)

set.seed(2024)
biom_16s_exp1_ev = rarefy_even_depth(biom_16s_exp1,
                   sample.size = min(sample_sums(biom_16s_exp1)),
                   rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
biom_16s_exp1_ev
sample_sums(biom_16s_exp1_ev)

CSSNorm(biom_16s_exp1) -> biom_16s_exp1_mSeq
head(otu_table(biom_16s_exp1_mSeq))

# 2) this is EXPERIMENT 2 16s dataset --------------------------------------------------------------
biom_16s_exp2 <- subset_samples(biom_16s_uparse_filt, exp2%in%c("exp2"))
otu_table(biom_16s_exp2) <- otu_table(biom_16s_exp2)[which(rowSums(otu_table(biom_16s_exp2)) >=1),] 
biom_16s_exp2
sample_data(biom_16s_exp2)

set.seed(2025)
biom_16s_exp2_ev = rarefy_even_depth(biom_16s_exp2,
                   sample.size = min(sample_sums(biom_16s_exp2)),
                   rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
biom_16s_exp2_ev
sample_sums(biom_16s_exp2_ev)

CSSNorm(biom_16s_exp2) -> biom_16s_exp2_mSeq
head(otu_table(biom_16s_exp2_mSeq))

# 3) this is EXPERIMENT 3 16s dataset --------------------------------------------------------------
biom_16s_exp3 <- subset_samples(biom_16s_uparse_filt, exp3%in%c("exp3"))
otu_table(biom_16s_exp3) <- otu_table(biom_16s_exp3)[which(rowSums(otu_table(biom_16s_exp3)) >=1),] 
biom_16s_exp3
sample_data(biom_16s_exp3)

set.seed(2026)
biom_16s_exp3_ev = rarefy_even_depth(biom_16s_exp3,
                   sample.size = min(sample_sums(biom_16s_exp3)),
                   rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
biom_16s_exp3_ev
sample_sums(biom_16s_exp3_ev)

CSSNorm(biom_16s_exp3) -> biom_16s_exp3_mSeq
head(otu_table(biom_16s_exp3_mSeq))

# >>> ALPHA DIVERSITY METRICS ----------------------------------------------------------------------

# coearcing into factors for plotting 
sample_data(biom_ITS_exp1)$year <- as.factor(sample_data(biom_ITS_exp1)$year)
sample_data(biom_ITS_exp1)$year <- factor(sample_data(biom_ITS_exp1)$year, levels=c("2015","2014","2010","2005","2000","1995"))
sample_data(biom_ITS_exp2)$year <- as.factor(sample_data(biom_ITS_exp2)$year)
sample_data(biom_ITS_exp2)$year <- factor(sample_data(biom_ITS_exp2)$year, levels=c("2015","2014","2010","2005","2000"))
sample_data(biom_16s_exp1)$year <- as.factor(sample_data(biom_16s_exp1)$year)
sample_data(biom_16s_exp1)$year <- factor(sample_data(biom_16s_exp1)$year, levels=c("2015","2014","2010","2005","2000","1995"))
sample_data(biom_16s_exp2)$year <- as.factor(sample_data(biom_16s_exp2)$year)
sample_data(biom_16s_exp2)$year <- factor(sample_data(biom_16s_exp2)$year, levels=c("2015","2014","2010","2005","2000"))

PlotAlpha <-function(dataframe, x, y, facet, maxyAxis, colbox, fillbox){
  sample_data(dataframe)$time <- sample_data(dataframe)$year %>% 
    recode_factor('2015' = "0", '2014' = "1", '2010' = "5", '2005' = "10", '2000' ="15", '1995'="20")
    as.data.frame(sample_data(dataframe)) %>%
    rename(facet=!!facet) %>% # you have to rename the facet becasue of the facet_grid
    ggplot(aes(x=get(x), y=get(y))) + # using get() to pass aes to ggplot 
        geom_boxplot(lwd = 0.7, color=colbox, fill=fillbox, alpha=0.6) +
          stat_summary(fun=mean, geom="point", shape=18, size=1.9, color="red", fill="red") +
          facet_grid(~facet, scales = "free_x", space="free_x") + 
    labs(x="Year") +
    ylim(0, maxyAxis) +
    theme_classic() +
    theme(strip.text.x = element_text(size = 9, face = "bold")) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 0, face = "bold", size = 8,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 9, face = "bold")) +
    grids(linetype = "dashed") +
    theme(legend.position = "none") -> alphaplot
  return(alphaplot)
}


# *** FIGURE 1 ------------------------------------------------------------------------------------
library(ggpubr)
library(cowplot)

# ggarrange(PlotAlpha(biom_ITS_exp1, x = "year",y = "richness","habitat", 650,"#332288","#332288") + labs(title = "Experiment 1",y="Observed richness"),
#           PlotAlpha(biom_ITS_exp2, x = "year",y = "richness","site", 650,"#332288","#332288") + labs(title = "Experiment 2", y=NULL),
#           PlotAlpha(biom_ITS_exp3, x = "site",y = "richness","habitat", 650,"#332288","#332288") + labs(title = "Experiment 3", y=NULL),
#           PlotAlpha(biom_16s_exp1, x = "year",y = "richness","habitat", 4500,"#D55E00","#D55E00") + labs(title = "Experiment 1",y="Observed richness"),
#           PlotAlpha(biom_16s_exp2, x = "year",y = "richness","site", 4500,"#D55E00","#D55E00") + labs(title = "Experiment 2", y=NULL),
#           PlotAlpha(biom_16s_exp3, x = "site",y = "richness","habitat", 4500,"#D55E00","#D55E00") + labs(title = "Experiment 3", y=NULL),
#           labels = c("A",NULL,NULL,"B",NULL,NULL),
#           #widths = c(1.1, 0.45, 0.75, 1.1, 0.45, 0.75),
#           widths = c(1.15, 0.505, 0.71, 1.15, 0.505, 0.71),
#           align = "none", ncol = 3, nrow = 2) -> Fig1_alpha_box
# 
# Fig1_alpha_box

plot_grid(PlotAlpha(biom_ITS_exp1, x = "time",y = "richness","habitat", 650,"#332288","#332288") + labs(title = "Experiment 1",y="Observed richness"),
          PlotAlpha(biom_ITS_exp2, x = "time",y = "richness","site", 650,"#332288","#332288") + labs(title = "Experiment 2", y=NULL),
          PlotAlpha(biom_ITS_exp3, x = "site",y = "richness","habitat", 650,"#332288","#332288") + labs(title = "Experiment 3", y=NULL, x="Site"),
          PlotAlpha(biom_16s_exp1, x = "time",y = "richness","habitat", 4500,"#D55E00","#D55E00") + labs(title = "Experiment 1",y="Observed richness"),
          PlotAlpha(biom_16s_exp2, x = "time",y = "richness","site", 4500,"#D55E00","#D55E00") + labs(title = "Experiment 2", y=NULL),
          PlotAlpha(biom_16s_exp3, x = "site",y = "richness","habitat", 4500,"#D55E00","#D55E00") + labs(title = "Experiment 3", y=NULL, x="Site"),
          labels = c("A","", "","B","",""),
          rel_widths = c(1.1, 0.555, 0.76, 1.1, 0.555, 0.76),
          ncol = 3, nrow = 2, 
          align = 'hv', 
          axis = 'l') -> Fig1_richness

Fig1_richness

plot_grid(PlotAlpha(biom_ITS_exp1, x = "time",y = "rarefied","habitat", 550,"#332288","#332288") + labs(title = "Experiment 1",y="Rarefied richness"),
          PlotAlpha(biom_ITS_exp2, x = "time",y = "rarefied","site", 550,"#332288","#332288") + labs(title = "Experiment 2", y=NULL),
          PlotAlpha(biom_ITS_exp3, x = "site",y = "rarefied","habitat", 550,"#332288","#332288") + labs(title = "Experiment 3", y=NULL, x="Site"),
          PlotAlpha(biom_16s_exp1, x = "time",y = "rarefied","habitat", 1750,"#D55E00","#D55E00") + labs(title = "Experiment 1",y="Rarefied richness"),
          PlotAlpha(biom_16s_exp2, x = "time",y = "rarefied","site", 1750,"#D55E00","#D55E00") + labs(title = "Experiment 2", y=NULL),
          PlotAlpha(biom_16s_exp3, x = "site",y = "rarefied","habitat", 1750,"#D55E00","#D55E00") + labs(title = "Experiment 3", y=NULL, x="Site"),
          labels = c("A","", "","B","",""),
          rel_widths = c(1.1, 0.555, 0.76, 1.1, 0.555, 0.76),
          ncol = 3, nrow = 2, 
          align = 'hv', 
          axis = 'l') -> Fig1_rarefied_rich

Fig1_rarefied_rich

plot_grid(PlotAlpha(biom_ITS_exp1, x = "time",y = "shannon","habitat", 5,"#332288","#332288") + labs(title = "Experiment 1",y="Rarefied richness"),
          PlotAlpha(biom_ITS_exp2, x = "time",y = "shannon","site", 5,"#332288","#332288") + labs(title = "Experiment 2", y=NULL),
          PlotAlpha(biom_ITS_exp3, x = "site",y = "shannon","habitat", 5,"#332288","#332288") + labs(title = "Experiment 3", y=NULL, x="Site"),
          PlotAlpha(biom_16s_exp1, x = "time",y = "shannon","habitat", 7.5,"#D55E00","#D55E00") + labs(title = "Experiment 1",y="Rarefied richness"),
          PlotAlpha(biom_16s_exp2, x = "time",y = "shannon","site", 7.5,"#D55E00","#D55E00") + labs(title = "Experiment 2", y=NULL),
          PlotAlpha(biom_16s_exp3, x = "site",y = "shannon","habitat", 7.5,"#D55E00","#D55E00") + labs(title = "Experiment 3", y=NULL, x="Site"),
          labels = c("A","", "","B","",""),
          rel_widths = c(1.1, 0.555, 0.76, 1.1, 0.555, 0.76),
          ncol = 3, nrow = 2, 
          align = 'hv', 
          axis = 'l') -> Fig1_shannon

Fig1_shannon

# Modeling Richness -------------------------------------------------------------------------------
library(nlcor)
library(mgcv)
library(MASS)
library(dplyr)

# create a new numeric variable called time
sample_data(biom_ITS_exp1)$time <- sample_data(biom_ITS_exp1)$year %>% 
  recode_factor('2015' = "0", '2014' = "1", '2010' = "5", '2005' = "10", '2000' ="15", '1995'="20")
sample_data(biom_16s_exp1)$time <- sample_data(biom_16s_exp1)$year %>% 
  recode_factor('2015' = "0", '2014' = "1", '2010' = "5", '2005' = "10", '2000' ="15", '1995'="20")

sample_data(biom_ITS_exp2)$time <- sample_data(biom_ITS_exp2)$year %>% 
  recode_factor('2015' = "0", '2014' = "1", '2010' = "5", '2005' = "10", '2000' ="15", '1995'="20")
sample_data(biom_16s_exp2)$time <- sample_data(biom_16s_exp2)$year %>% 
  recode_factor('2015' = "0", '2014' = "1", '2010' = "5", '2005' = "10", '2000' ="15", '1995'="20")


#Extracting datasets by habitat
ExtractData <-function(dataframe, variable){
  if (is.null(variable)){
      as.data.frame(as.matrix(sample_data(dataframe))) -> output
      output$richness <- as.numeric(as.character(output$richness))
      output$rarefied <- as.numeric(as.character(output$rarefied))
      output$roundRare <- as.numeric(as.character(round(output$rarefied, 0)))
      output$shannon <- as.numeric(as.character(output$shannon))
      output$time <- as.numeric(as.character(output$time))
        return(output)
    } else {
    dataframe %>% 
      subset_samples(habitat%in%c("PS")) -> output_PS
      as.data.frame(as.matrix(sample_data(output_PS))) -> output_PS
      output_PS$richness <- as.numeric(as.character(output_PS$richness))
      output_PS$rarefied <- as.numeric(as.character(output_PS$rarefied))
      output_PS$roundRare <- as.numeric(as.character(round(output_PS$rarefied, 0)))
      output_PS$shannon <- as.numeric(as.character(output_PS$shannon))
      output_PS$time <- as.numeric(as.character(output_PS$time))
    dataframe %>% 
      subset_samples(habitat%in%c("DF")) -> output_DF
      as.data.frame(as.matrix(sample_data(output_DF))) -> output_DF
      output_DF$richness <- as.numeric(as.character(output_DF$richness))
      output_DF$rarefied <- as.numeric(as.character(output_DF$rarefied))
      output_DF$roundRare <- as.numeric(as.character(round(output_DF$rarefied, 0)))
      output_DF$shannon <- as.numeric(as.character(output_DF$shannon))
      output_DF$time <- as.numeric(as.character(output_DF$time))
      if (variable == "PS"){
          return(output_PS)
      } else {
      return(output_DF)
    }
  }
}


ExtractData(biom_ITS_exp1, "PS") ->df_ITS_exp1_PS
ExtractData(biom_ITS_exp1, "DF") ->df_ITS_exp1_DF
ExtractData(biom_16s_exp1, "PS") ->df_16s_exp1_PS
ExtractData(biom_16s_exp1, "DF") ->df_16s_exp1_DF

df_ITS_exp1_PS
df_ITS_exp1_DF
df_16s_exp1_PS
df_16s_exp1_DF

ExtractData(biom_ITS_exp2, NULL) -> df_ITS_exp2
ExtractData(biom_16s_exp2, NULL) -> df_16s_exp2

df_ITS_exp2
df_16s_exp2

# >>> MODELING TIME --------------------------------------------------------------------------------

# # A more realistic example: split a data frame into pieces, fit a
# # model to each piece, summarise and extract R^2
# mtcars %>%
#   split(.$cyl) %>%
#   map(~ lm(mpg ~ wt, data = .x)) %>%
#   map(summary) %>%
#   map_dbl("r.squared")


# Extracting coefficients form the different models ------------------------------------------------
ExtractCoeff <-function(dataframe, y){
  fit_linear <- lm(get(y) ~ time, data = dataframe)
  fit_quad <- lm(get(y) ~ poly(time, 2, raw=TRUE), data = dataframe)
  fit_cubic <- lm(get(y) ~ poly(time, 3, raw=TRUE), data = dataframe)
  fit_4 <- lm(get(y) ~ poly(time, 4, raw=TRUE), data = dataframe)
  fit_log <- lm(get(y) ~ log(time +1), data = dataframe)
  rbind(summary(fit_linear)$coefficients,
        summary(fit_quad)$coefficients,
        summary(fit_cubic)$coefficients,
        summary(fit_4)$coefficients,
        summary(fit_log)$coefficients)-> df
  df <- as.data.frame(df)
  df$model <- c(rep("linear", 2),rep("quadratic", 3),
                rep("cubic", 4),rep("IV degree", 5),
                rep("log", 2))
  result <- c()
  for (x in 1:length(df$`Pr(>|t|)`)){
    if (df$`Pr(>|t|)`[x] >0.05 &
        df$`Pr(>|t|)`[x] <= 0.1){
      result[x] <- "."
    }else if (df$`Pr(>|t|)`[x]>0.01 &
              df$`Pr(>|t|)`[x] <= 0.05){
      result[x] <- "*"
    } else if (df$`Pr(>|t|)`[x] > 0.001 &
               df$`Pr(>|t|)`[x] <= 0.01) {
      result[x] <- "**"
    } else if (df$`Pr(>|t|)`[x] <= 0.001) {
      result[x] <- "***"
    } else {
      result[x] <- "ns"
    }
  }
  df$signif <- as.character(result)
  df2 <- data.frame(matrix(ncol=1,nrow=5, 
              dimnames=list(c("linear", "quandratic",
              "cubic","IV degree", "log"), c("V1"))))
        df2[1,1] <-anova(fit_linear)$'F value'[1]
        df2[2,1] <-anova(fit_quad)$'F value'[1]
        df2[3,1] <-anova(fit_cubic)$'F value'[1]
        df2[4,1] <-anova(fit_4)$'F value'[1]
        df2[5,1] <-anova(fit_log)$'F value'[1]
        df2[1,2] <-anova(fit_linear)$'Pr(>F)'[1]
        df2[2,2] <-anova(fit_quad)$'Pr(>F)'[1]
        df2[3,2] <-anova(fit_cubic)$'Pr(>F)'[1]
        df2[4,2] <-anova(fit_4)$'Pr(>F)'[1]
        df2[5,2] <-anova(fit_log)$'Pr(>F)'[1]
        colnames(df2) <- c("F.value", "p.value")
  require(rcompanion)
  compareLM(fit_linear, fit_quad, 
            fit_cubic, fit_4, fit_log) -> result2
  cbind(result2$Models, result2$Fit.criteria) -> result2
  return(list(coefficients =df, anova=df2, Fit = result2))
}


ExtractCoeff(df_ITS_exp1_PS, "roundRare") -> fit_info_ITS_exp1_PS_rare 
fit_info_ITS_exp1_PS_rare
ExtractCoeff(df_ITS_exp1_DF,"roundRare") -> fit_info_ITS_exp1_DF_rare 
fit_info_ITS_exp1_DF_rare 
ExtractCoeff(df_16s_exp1_PS,"roundRare") -> fit_info_16s_exp1_PS_rare  
fit_info_16s_exp1_PS_rare 
ExtractCoeff(df_16s_exp1_DF, "roundRare") -> fit_info_16s_exp1_DF_rare  
fit_info_16s_exp1_DF_rare 

ExtractCoeff(df_ITS_exp1_PS, "richness") -> fit_info_ITS_exp1_PS_obs
fit_info_ITS_exp1_PS_obs
ExtractCoeff(df_ITS_exp1_DF,"richness") -> fit_info_ITS_exp1_DF_obs
fit_info_ITS_exp1_DF_obs
ExtractCoeff(df_16s_exp1_PS,"richness") -> fit_info_16s_exp1_PS_obs 
fit_info_16s_exp1_PS_obs
ExtractCoeff(df_16s_exp1_DF, "richness") -> fit_info_16s_exp1_DF_obs 
fit_info_16s_exp1_DF_obs

ExtractCoeff(df_ITS_exp1_PS, "shannon") -> fit_info_ITS_exp1_PS_shan 
fit_info_ITS_exp1_PS_shan
ExtractCoeff(df_ITS_exp1_DF,"shannon") -> fit_info_ITS_exp1_DF_shan 
fit_info_ITS_exp1_DF_shan 
ExtractCoeff(df_16s_exp1_PS,"shannon") -> fit_info_16s_exp1_PS_shan 
fit_info_16s_exp1_PS_shan 
ExtractCoeff(df_16s_exp1_DF, "shannon") -> fit_info_16s_exp1_DF_shan  
fit_info_16s_exp1_DF_shan 


# calculating R2 and AICc form glm and glm.nb models  
CalculateAICR2 <- function(list){
  require(rsq)
  df_p <- data.frame(matrix(ncol=4,nrow=6, 
          dimnames=list(c("ITS_PS", "ITS_DF","16s_PS","16s_DF", "ITS_exp2", "16s_exp2"))))
  colnames(df_p) <- c("AIC", "AICc","BIC","Adj.R2")
    df_p$AIC <- unlist(lapply(list,
           function(x) AIC(glm(round(x$rarefied, 0) ~ time, data = x, family="poisson"))))
    df_p$AICc <- unlist(lapply(list,
           function(x) AICc(glm(round(x$rarefied, 0) ~ time, data = x, family="poisson"))))
    df_p$BIC <- unlist(lapply(list,
           function(x) BIC(glm(round(x$rarefied, 0) ~ time, data = x, family="poisson"))))
    df_p$Adj.R2 <- unlist(lapply(list,
           function(x) rsq(glm(round(x$rarefied, 0) ~ time, data = x, family="poisson"), 
                           adj=TRUE)))
    df_nb <- data.frame(matrix(ncol=4,nrow=6, 
            dimnames=list(c("ITS_PS", "ITS_DF","16s_PS","16s_DF", "ITS_exp2", "16s_exp2"))))
    colnames(df_nb) <- c("AIC", "AICc","BIC","Adj.R2")
    df_nb$AIC <- unlist(lapply(list,
           function(x) AIC(glm.nb(round(x$rarefied, 0) ~ time , data = x))))
    df_nb$AICc <- unlist(lapply(list,
            function(x) AICc(glm.nb(round(x$rarefied, 0) ~ time , data = x))))
    df_nb$BIC <- unlist(lapply(list,
            function(x) BIC(glm.nb(round(x$rarefied, 0) ~ time , data = x))))
    df_nb$Adj.R2 <- unlist(lapply(list,
            function(x) rsq(glm.nb(round(x$rarefied, 0) ~ time , data = x), 
                         adj=TRUE)))
  return(list(poisson =df_p, negaive.binomial = df_nb))
}

# Calculate values
CalculateAICR2(list(df_ITS_exp1_PS,
                    df_ITS_exp1_DF,
                    df_16s_exp1_PS,
                    df_16s_exp1_DF,
                    df_ITS_exp2,
                    df_16s_exp2)) -> result_pois_nb
result_pois_nb

# MODELS PLOTS -------------------------------------------------------------------------------------
# function for potting regression models -----------------------------------------------------------

PlotRegression <-function(dataframe, x, y, method, family, formula, color){
  regplot <- ggplot(dataframe, aes(x=get(x), y=get(y))) +
    geom_point(alpha = 0.8, shape=16, size=2, color=color) +
    expand_limits(y = 0) +
    scale_x_continuous(breaks=c(0,1,5,10,15,20)) +
    #labs(x="Time", y="Shannon index") +
    theme_classic() +
    theme(strip.text.x = element_text(size = 9, face = "bold")) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 9, face = "bold")) +
    grids(linetype = "dashed") -> regplot
  if (is.null(family)) 
    regplot + geom_smooth(method = method, formula = formula, color="red", size=0.7) -> regplot1
  if (is.null(formula))
    regplot + geom_smooth(method = method, method.args = list(family = family), color="red", size=0.7) -> regplot1
  if (is.null(family) & is.null(formula))
    regplot + geom_smooth(method = method, color="red", size=0.7) -> regplot1
    return(regplot1)
}


# testing for different models ---------------------------------------------------------------------
PlotRegression(df_ITS_exp1_PS, x = "time", y = "richness" ,"lm", NULL, formula = y ~ x, "#332288") + 
  labs(title="y ~ x") + xlab("Time") + ylab("Observed Richness") # linear

PlotRegression(df_ITS_exp1_PS, x = "time", y = "richness", "glm", "poisson", NULL, "#332288") 

PlotRegression(df_ITS_exp1_PS, x = "time", y = "richness" , "glm", "quasipoisson", NULL, "#332288") + 
  labs(title="quasipoisson") + xlab("Time") + ylab("Observed Richness") + 
  annotate("text", x=xAxis,y=yAaxis,label=paste("italic(R)[Adj]^2 ==", round(rsq(fit_quasip, adj=TRUE), 2),
                                                "~AICc==", round(AICc(fit_quasip), 2)), size=2.5, parse=TRUE)# quasipoisson

PlotRegression(df_ITS_exp1_PS, x = "time", y = "richness" , "glm.nb", NULL, NULL, "#332288") + 
  labs(title="Negative Binomial") + xlab("Time") + ylab("Observed Richness") # negative binomial


# Function to multiplot using ggarrange 
PlotModel <-function(dataframe, x, y, method, color, xAxis, yAaxis, xLabel){
  noquote(y) -> y2
  require(wiqid)
  require(rsq)
      fit_linear <- lm(get(y) ~ get(x), data = dataframe)
      fit_quad <- lm(get(y) ~ poly(get(x), 2, raw=TRUE), data = dataframe)
      fit_cubic <- lm(get(y) ~ poly(get(x), 3, raw=TRUE), data = dataframe)
      fit_4 <- lm(get(y) ~ poly(get(x), 4, raw=TRUE), data = dataframe)
      fit_log <- lm(get(y) ~ log(get(x) +1), data = dataframe)
      #fit_pois <- glm(y2 ~ get(x), data=dataframe, family="poisson") # here is the problem
      fit_nb <- glm.nb(get(y) ~ get(x), data = dataframe)
  MultiPlots <- ggpubr::ggarrange(
    PlotRegression(dataframe, x = "time", y = y, "lm", NULL, y ~ x, color) + labs(title="y ~ x") + 
       xlab("Time") + ylab(xLabel) +
      annotate("text", x=xAxis,y=yAaxis,label=paste("italic(R)[Adj]^2 ==", round(summary(fit_linear)$adj.r.squared, 2),
                                                    "~AICc==", round(AICc(fit_linear), 2)), size=2.5, parse=TRUE),
    PlotRegression(dataframe, x = "time", y = y, "lm", NULL, y ~ poly(x, 2), color) + labs(title="y ~ poly(x, 2)") +
      xlab("Time") + ylab(xLabel) +
      annotate("text", x=xAxis,y=yAaxis,label=paste("italic(R)[Adj]^2 ==", round(summary(fit_quad)$adj.r.squared, 2),
                                                    "~AICc==", round(AICc(fit_quad), 2)), size=2.5, parse=TRUE),
    PlotRegression(dataframe, x = "time", y = y, "lm", NULL, y ~ poly(x, 3), color) + labs(title="y ~ poly(x, 3)") +
      xlab("Time") + ylab(xLabel) +
      annotate("text", x=xAxis,y=yAaxis,label=paste("italic(R)[Adj]^2 ==", round(summary(fit_cubic)$adj.r.squared, 2),
                                                    "~AICc==", round(AICc(fit_cubic), 2)), size=2.5, parse=TRUE),
    PlotRegression(dataframe, x = "time", y = y, "lm", NULL, y ~ poly(x, 4) ,  color) + labs(title="y ~ poly(x, 4)") +
      xlab("Time") + ylab(xLabel) +
      annotate("text", x=xAxis,y=yAaxis,label=paste("italic(R)[Adj]^2 ==", round(summary(fit_4)$adj.r.squared, 2),
                                                    "~AICc==", round(AICc(fit_4), 2)), size=2.5, parse=TRUE),
    PlotRegression(dataframe, x = "time", y = y, "lm", NULL, y ~ log(x+1), color) + labs(title="y ~ log(x+1)") +
      xlab("Time") + ylab(xLabel) +
      annotate("text", x=xAxis,y=yAaxis,label=paste("italic(R)[Adj]^2 ==", round(summary(fit_log)$adj.r.squared, 2),
                                                    "~AICc==", round(AICc(fit_log), 2)), size=2.5, parse=TRUE),
#    PlotRegression(dataframe, x = "time", y = y, "glm", "poisson", NULL, color) + labs(title="poisson") +
#      xlab("Time") + ylab(xLabel) +
#      annotate("text", x=xAxis,y=yAaxis,label=paste("italic(R)[Adj]^2 ==", round(rsq(fit_pois, adj=TRUE), 2),
#                                                    "~AICc==", round(AICc(fit_pois), 2)), size=2.5, parse=TRUE),
    PlotRegression(dataframe, x = "time", y = y, "glm.nb", NULL, NULL, color) + labs(title="negative binomial") +
      xlab("Time") + ylab(xLabel) +
      annotate("text", x=xAxis,y=yAaxis,label=paste("italic(R)[Adj]^2 ==", round(rsq(fit_nb, adj=TRUE), 2),
                                                    "~AICc==", round(AICc(fit_nb), 2)), size=2.5, parse=TRUE),
    widths = c(1,1,1,1,1,1),
    align = "hv", ncol = 3, nrow = 2)
  return(MultiPlots)
}

# testing 
PlotModel(dataframe = df_ITS_exp1_PS, x = "time", y = "roundRare", method = "lm",
                  color = "#332288", xAxis = 10, yAaxis = 700, xLabel = "Observed")

# FIGURE S5 to S8 Rarefied Richness ---------------------------------------------------------------- 
PlotModel(df_ITS_exp1_PS,x = "time", y = "roundRare", "lm","#332288", 10, 500, "Rarefied richness") -> models_PS_fungi_rare
models_PS_fungi_rare
PlotModel(df_ITS_exp1_DF, x = "time", y = "roundRare", "lm","#332288", 7.5, 500, "Rarefied richness") -> models_DF_fungi_rare
models_DF_fungi_rare
PlotModel(df_16s_exp1_PS, x = "time", y = "roundRare", "lm", "#D55E00", 10, 100, "Rarefied richness") -> models_PS_prok_rare
models_PS_prok_rare
PlotModel(df_16s_exp1_DF, x = "time", y = "roundRare", "lm", "#D55E00",7.5, 100, "Rarefied richness") -> models_DF_prok_rare
models_DF_prok_rare

# FIGURE XXX Observed Richness ---------------------------------------------------------------------
PlotModel(df_ITS_exp1_PS,x = "time", y = "richness", "lm","#332288", 10, 600,"Observed richness") -> models_PS_fungi
models_PS_fungi
PlotModel(df_ITS_exp1_DF,x = "time", y = "richness", "lm","#332288", 7.5, 600, "Observed richness") -> models_DF_fungi
models_DF_fungi
PlotModel(df_16s_exp1_PS, x = "time", y = "richness", "lm", "#D55E00", 10, 100, "Observed richness") -> models_PS_prok
models_PS_prok
PlotModel(df_16s_exp1_DF, x = "time", y = "richness", "lm", "#D55E00",7.5, 100, "Observed richness") -> models_DF_prok
models_DF_prok


# plotting models for Exp2
PlotModel(df_ITS_exp2,x = "time", y = "roundRare", "lm","#332288", 7.5, 600,"Rarefied richness") -> models_exp2_fungi
models_exp2_fungi

PlotModel(df_16s_exp2,x = "time", y = "roundRare", "lm","#D55E00", 7.5, 100,"Rarefied richness") -> models_exp2_prok
models_exp2_prok


# # FIGURE S5 to S8 Shannon Index ------------------------------------------------------------------
# PlotModel(df_ITS_exp1_PS,x = "time", y = "shannon", "lm","#332288", 10, 0.2, "Shannon index") -> models_PS_fungi_Shan
# models_PS_fungi_Shan
# PlotModel(df_ITS_exp1_DF,x = "time", y = "shannon", "lm","#332288", 7.5, 0.2, "Shannon index") -> models_DF_fungi_Shan
# models_DF_fungi_Shan
# PlotModel(df_16s_exp1_PS, x = "time", y = "shannon", "lm", "#D55E00", 10, 0.2, "Shannon index") -> models_PS_prok_Shan
# models_PS_prok_Shan
# PlotModel(df_16s_exp1_DF, x = "time", y = "shannon", "lm", "#D55E00",7.5, 0.2, "Shannon index") -> models_DF_prok_Shan
# models_DF_prok_Shan


# compare models 
CompareModel <-function(dataframe, x, y){
  require(rcompanion)
  fit_linear <- lm(get(y) ~ get(x), data = dataframe)
  fit_quad <- lm(get(y) ~ poly(get(x), 2, raw=TRUE), data = dataframe)
  fit_cubic <- lm(get(y) ~ poly(get(x), 3, raw=TRUE), data = dataframe)
  fit_4 <- lm(get(y) ~ poly(get(x), 4, raw=TRUE), data = dataframe)
  fit_log <- lm(get(y) ~ log(get(x) +1), data = dataframe)
    compareLM(fit_linear, fit_quad, 
              fit_cubic, fit_4, fit_log) -> result
    cbind(result$Models, result$Fit.criteria) -> result
return(result)
}


list(exp1_ITS_PS = CompareModel(df_ITS_exp1_PS, "time", "roundRare"), 
     exp1_ITS_DF = CompareModel(df_ITS_exp1_DF, "time", "roundRare"),
     exp1_16s_PS = CompareModel(df_16s_exp1_PS, "time", "roundRare"),
     exp1_16s_DF = CompareModel(df_16s_exp1_DF, "time", "roundRare"), 
     exp2_ITS_DF1 = CompareModel(df_ITS_exp2, "time", "roundRare"),
     exp2_16s_DF1 = CompareModel(df_16s_exp2, "time", "roundRare")) -> exp1_CompareModel

exp1_CompareModel
result_pois_nb

CompareModel(df_ITS_exp1_PS, "time", "richness")
CompareModel(df_ITS_exp1_DF, "time", "richness")
CompareModel(df_16s_exp1_PS, "time", "richness")
CompareModel(df_16s_exp1_DF, "time", "richness")

# >>> IDENTIFYING BEST MODELS ----------------------------------------------------------------------
# for the fungi ------------------------------------------------------------------------------------

library(car)
fit_nb_ITS_DF <- glm.nb(roundRare ~ time, data = df_ITS_exp1_DF) # negative binomial fungi_DF
summary(fit_nb_ITS_DF)
car::Anova(fit_nb_ITS_DF)
anova(fit_nb_ITS_DF)
rsq::rsq(fit_nb_ITS_DF, adj=FALSE)
rsq::rsq(fit_nb_ITS_DF, adj=TRUE)

fit_nb_ITS_PS <- glm.nb(roundRare ~ time, data = df_ITS_exp1_PS) # negative binomial fungi_PS
summary(fit_nb_ITS_PS)
anova(fit_nb_ITS_PS)
car::Anova(fit_nb_ITS_PS)
rsq::rsq(fit_nb_ITS_PS , adj=FALSE)
rsq::rsq(fit_nb_ITS_PS , adj=TRUE)

# diagnostic plots 
par(mfrow=c(2,2)); plot(fit_nb_ITS_DF)
par(mfrow=c(2,2)); plot(fit_nb_ITS_PS)

# Calculating scaled residuals
library(DHARMa)
simulres_ITS_DF <- simulateResiduals(fittedModel = fit_nb_ITS_DF, n = 1000)
plot(simulres_ITS_DF)

simulres_ITS_PS <- simulateResiduals(fittedModel = fit_nb_ITS_PS, n = 1000)
plot(simulres_ITS_PS)

# Negative binomial equations 
estimate_nb_ITS_DF <- cbind(Estimate = coef(fit_nb_ITS_DF), confint(fit_nb_ITS_DF))
estimate_nb_ITS_DF
summary(fit_nb_ITS_DF)

library(effects)
summary(effects::allEffects(fit_nb_ITS_DF))
# same as doing the follwoing 
# fit_nb_ITS_DF_mean <- exp(predict(fit_nb_ITS_DF, newdata = data.frame(time = as.numeric(c("0","1", "5","10","15")))))

# to interpet the results with the log or logit links, we need to 
# exponentiate the coefficients to get the incidence rate ratios
library(jtools)
summ(fit_nb_ITS_DF, exp = T)
summ(fit_nb_ITS_PS, exp = T)

# figure diagnostic plots fungi 
par(mfrow=c(2,3))
testUniformity(simulres_ITS_DF)
plotResiduals(simulres_ITS_DF$fittedPredictedResponse, simulres_ITS_DF$scaledResiduals,
              main = "Residual vs. predicted\nlines should match",
              xlab="Predicted values (rank transformed)", 
              ylab= "Standardized residual")
testDispersion(simulres_ITS_DF)

testUniformity(simulres_ITS_PS)
plotResiduals(simulres_ITS_PS$fittedPredictedResponse, simulres_ITS_PS$scaledResiduals,
              main = "Residual vs. predicted\nlines should match",
              xlab= "Predicted values (rank transformed)", 
              ylab= "Standardized residual")
testDispersion(simulres_ITS_PS)

# for 16s ------------------------------------------------------------------------------------------
fit_quad_16s_DF <- lm(roundRare ~ poly(time, 2, raw=TRUE), data = df_16s_exp1_DF) #quadratic for bacteria_DF
coef(lm(roundRare ~ poly(time, 2, raw=TRUE), data = df_16s_exp1_DF))
coef(lm(roundRare ~ time + I(time^2), data = df_16s_exp1_DF))

summary(fit_quad_16s_DF)
anova(fit_quad_16s_DF)
confint(fit_quad_16s_DF, level=0.95)

# # calculate the new x values using predict.poly()
# x_poly <- stats:::predict.poly(object = poly(df_16s_exp1_DF$time, 2), newdata = 3)
# coef(fit_quad_16s_DF)[1] + coef(fit_quad_16s_DF)[2] * x_poly[1] + coef(fit_quad_16s_DF)[3] * x_poly[2]

par(mfrow=c(2,2)); plot(fit_quad_16s_DF)
 
fit_log_16s_PS <- lm(roundRare ~ log(time +1), data = df_16s_exp1_PS) # log for bacteira_PS
coef(lm(roundRare ~ log(time +1), data = df_16s_exp1_PS))

summary(fit_log_16s_PS)
anova(fit_log_16s_PS)
par(mfrow=c(2,2)); plot(fit_log_16s_PS)

# *** FIGURE 2 fungal models -----------------------------------------------------------------------
PlotBestModel <-function(x, y){
ggarrange(PlotRegression(df_ITS_exp1_DF, x = x, y = y, "glm.nb", NULL, NULL, "#332288") + 
            labs(title="Fungi DF", x="Time", y="Rarefied richness") +
            annotate("text", x=7.5, y=500,
                     label=paste("italic(R)[Adj]^2 ==",round(rsq(fit_nb_ITS_DF, adj=TRUE), 2),
                                 "~italic(P)<=", round(0.0001, 4)), #0.0001182
                     size=2.5, parse=TRUE),
            PlotRegression(df_ITS_exp1_PS, x = x, y = y, "glm.nb", NULL, NULL, "#332288") + 
            labs(title="Fungi PS", x="Time", y="Rarefied richness") +
            annotate("text", x=10,y=500, 
                     label=paste("italic(R)[Adj]^2 ==", round(rsq(fit_nb_ITS_PS, adj=TRUE), 2),
                                 "~italic(P)<=", round(0.0001, 4)), #0.000000004632
                     size=2.5, parse=TRUE),
          PlotRegression(df_16s_exp1_DF, x = x, y = y, "lm", NULL, y ~ poly(x, 2), "#D55E00") + 
            labs(title="Prokaryotes DF", x="Time", y="Rarefied richness") +
            annotate("text", x=7.5,y=100,
                     label=paste("italic(R)[Adj]^2 ==", round(summary(fit_quad_16s_DF)$adj.r.squared, 2),
                                 "~italic(P)<=",round(0.0020, 4)), #0.001936
                     size=2.5, parse=TRUE),
          PlotRegression(df_16s_exp1_PS, x = x, y = y, "lm", NULL, y ~ log(x + 1), "#D55E00") + 
            labs(title="Prokaryotes PS", x="Time", y="Rarefied richness") +
            annotate("text", x=10, y=100,
                     label=paste("italic(R)[Adj]^2 ==", round(summary(fit_log_16s_PS)$adj.r.squared, 2),
                                 "~italic(P)<=", round(0.0003, 4)), #0.0002506
                     size=2.5, parse=TRUE),
          labels = c("A", "B", "C", "D"),
          widths = c(1,1,1,1),
          align = "hv", 
          ncol = 2,
          nrow = 2) -> Fig_models
      return(Fig_models)
}
  
PlotBestModel(x = "time", y = "roundRare") -> Fig_2_models
Fig_2_models


# INDICATOR SPECIES ANALYSIS -----------------------------------------------------------------------
subset_samples(biom_ITS_exp1, habitat%in%c("PS")) -> biom_ITS_exp1_PS
otu_table(biom_ITS_exp1_PS) <- otu_table(biom_ITS_exp1_PS)[which(rowSums(otu_table(biom_ITS_exp1_PS)) > 0),] 

subset_samples(biom_ITS_exp1, habitat%in%c("DF")) -> biom_ITS_exp1_DF
otu_table(biom_ITS_exp1_DF) <- otu_table(biom_ITS_exp1_DF)[which(rowSums(otu_table(biom_ITS_exp1_DF)) > 0),] 

subset_samples(biom_16s_exp1, habitat%in%c("PS")) -> biom_16s_exp1_PS
otu_table(biom_16s_exp1_PS) <- otu_table(biom_16s_exp1_PS)[which(rowSums(otu_table(biom_16s_exp1_PS)) > 0),] 

subset_samples(biom_16s_exp1, habitat%in%c("DF")) -> biom_16s_exp1_DF
otu_table(biom_16s_exp1_DF) <- otu_table(biom_16s_exp1_DF)[which(rowSums(otu_table(biom_16s_exp1_DF)) > 0),] 

#indicator species analysis (isa) Fungi by year
library(indicspecies)

GetIndicators <-function(dataframe){
  require(phyloseq); require(indicspecies)
  otu <- as.data.frame(otu_table(dataframe))
  metadata = as(sample_data(dataframe), "data.frame")
  #taxa <- as.data.frame(as.matrix(tax_table(dataframe)))
      multipatt <- multipatt(t(otu), metadata$year, func = "r.g",
                       control=how(nperm=999))
  multipatt -> multipatt_fdr
  multipatt_fdr$sign$p.value <- p.adjust(multipatt_fdr$sign$p.value, "fdr")
  multipatt_fdr$sign[which(
    multipatt_fdr$sign$p.value <= 0.05), ] -> indicator_taxa
  return(indicator_taxa)
}

GetIndicators(biom_ITS_exp1_PS) -> ind_ITS_PS
GetIndicators(biom_ITS_exp1_DF)-> ind_ITS_DF
GetIndicators(biom_16s_exp1_PS)-> ind_16s_PS
GetIndicators(biom_16s_exp1_DF)-> ind_16s_DF
# no indicators found for variable time after pvalue correction

# BETA DIVERSITY -----------------------------------------------------------------------------------
# add another variable with the site
colnames(sample_data(biom_ITS_exp1_ev))[5] <- "block"
colnames(sample_data(biom_16s_exp1_ev))[5] <- "block"
# sample_data(biom_ITS_exp1_ev)$block <- dplyr::recode(sample_data(
#   biom_ITS_exp1_ev)$site, PS1 = "1", PS2 = "2", PS3 = "3",DF1 ="1", DF2="2", DF3="3")
# sample_data(biom_16s_exp1_ev)$block <- dplyr::recode(sample_data(
#   biom_16s_exp1_ev)$site, PS1 = "1", PS2 = "2", PS3 = "3",DF1 ="1", DF2="2", DF3="3")

# CAP model fungi ----------------------------------------------------------------------------------
otu_ITS_exp1_ev <- as.data.frame(otu_table(biom_ITS_exp1_ev))
metadata_ITS_exp1_ev = as(sample_data(biom_ITS_exp1_ev), "data.frame")

p.adjust(anova(fullModel_ITS)$`Pr(>F)`, "bonferroni")

# calculating full model with random factor
fullModel_ITS <- vegan::capscale(t(otu_ITS_exp1_ev) ~ year * habitat + Condition(block), metadata_ITS_exp1_ev, distance = "bray")
anova(fullModel_ITS, by="terms")
anova(fullModel_ITS)
round(p.adjust(anova(fullModel_ITS)$`Pr(>F)`, "bonferroni"), 4)
RsquareAdj(fullModel_ITS)

# backward selection
parsimoniousModel_ITS <- ordistep(fullModel_ITS, direction = "backward", Pout = 0.05, permutations = how(nperm = 999))
parsimoniousModel_ITS

# extract the parsimonious model formula for later use
FormulaParsMod_ITS <- formula(parsimoniousModel_ITS)[1:3]
FormulaParsMod_ITS

# the full model = the most parsimonious 
fullModel_ITS_1 <- vegan::capscale(t(otu_ITS_exp1_ev) ~ year + habitat + Condition(block), metadata_ITS_exp1_ev, distance = "bray")
fullModel_ITS_1
anova(fullModel_ITS_1, permutations=999)
round(p.adjust(anova(fullModel_ITS_1)$`Pr(>F)`, "bonferroni"), 4)
anova(fullModel_ITS_1, permutations=999, by = "terms")
RsquareAdj(fullModel_ITS_1)

# to test terms order in the  model ...
fullModel_ITS_2 <- vegan::capscale(t(otu_ITS_exp1_ev) ~ habitat + Condition(block + year), metadata_ITS_exp1_ev, distance = "bray")
anova(fullModel_ITS_2, permutations=999)
round(p.adjust(anova(fullModel_ITS_2)$`Pr(>F)`, "bonferroni"), 4)
anova(fullModel_ITS_2, permutations=999, by = "terms")
RsquareAdj(fullModel_ITS_2)

fullModel_ITS_3 <- vegan::capscale(t(otu_ITS_exp1_ev) ~ year + Condition(block + habitat), metadata_ITS_exp1_ev, distance = "bray")
anova(fullModel_ITS_3, permutations=999)
round(p.adjust(anova(fullModel_ITS_3)$`Pr(>F)`, "bonferroni"), 4)
anova(fullModel_ITS_3, permutations=999, by = "terms")
RsquareAdj(fullModel_ITS_3)

# fullModel_ITS_4 <- vegan::capscale(t(otu_ITS_exp1_ev) ~ year:habitat + Condition(site + year + habitat), metadata_ITS_exp1_ev, distance = "bray")
# anova(fullModel_ITS_4, permutations=999)
# anova(fullModel_ITS_4, permutations=999, by = "terms")
# RsquareAdj(fullModel_ITS_4)

# CAP model prokaryotes ----------------------------------------------------------------------------
otu_16s_exp1_ev <- as.data.frame(otu_table(biom_16s_exp1_ev))
metadata_16s_exp1_ev = as(sample_data(biom_16s_exp1_ev), "data.frame")

# calculating full model with random factor
fullModel_16s <- vegan::capscale(t(otu_16s_exp1_ev) ~ year * habitat + Condition(block), metadata_16s_exp1_ev, distance = "bray")
anova(fullModel_16s, by="terms")
anova(fullModel_16s)
round(p.adjust(anova(fullModel_16s)$`Pr(>F)`, "bonferroni"), 4)
RsquareAdj(fullModel_16s)

# backward selection
parsimoniousModel_16s <- ordistep(fullModel_16s, direction = "backward", Pout = 0.05, permutations = how(nperm = 999))
parsimoniousModel_16s

# extract the parsimonious model formula for later use
FormulaParsMod_16s <- formula(parsimoniousModel_16s)[1:3]
FormulaParsMod_16s

# the full model = the most parsimonious 
fullModel_16s_1 <- vegan::capscale(t(otu_16s_exp1_ev) ~ year + habitat + year:habitat, metadata_16s_exp1_ev, distance = "bray")
fullModel_16s_1
anova(fullModel_16s_1, permutations=999)
round(p.adjust(anova(fullModel_16s_1)$`Pr(>F)`, "bonferroni"), 4)
anova(fullModel_16s_1, permutations=999, by = "terms")
RsquareAdj(fullModel_16s_1)

# to test terms order in the  model ...
fullModel_16s_2 <- vegan::capscale(t(otu_16s_exp1_ev) ~ habitat + Condition(year), metadata_16s_exp1_ev, distance = "bray")
anova(fullModel_16s_2, permutations=999)
round(p.adjust(anova(fullModel_16s_2)$`Pr(>F)`, "bonferroni"), 4)
anova(fullModel_16s_2, permutations=999, by = "terms")
RsquareAdj(fullModel_16s_2)

fullModel_16s_3 <- vegan::capscale(t(otu_16s_exp1_ev) ~ year + Condition(habitat), metadata_16s_exp1_ev, distance = "bray")
anova(fullModel_16s_3, permutations=999)
round(p.adjust(anova(fullModel_16s_3)$`Pr(>F)`, "bonferroni"), 4)
anova(fullModel_16s_3, permutations=999, by = "terms")
RsquareAdj(fullModel_16s_3)

fullModel_16s_4 <- vegan::capscale(t(otu_16s_exp1_ev) ~ year:habitat + Condition(year + habitat), metadata_16s_exp1_ev, distance = "bray")
anova(fullModel_16s_4, permutations=999)
round(p.adjust(anova(fullModel_16s_4)$`Pr(>F)`, "bonferroni"), 4)
anova(fullModel_16s_4, permutations=999, by = "terms")
RsquareAdj(fullModel_16s_4)

# position of the annotation layer -----------------------------------------------------------------
# https://github.com/tidyverse/ggplot2/issues/1244

# Annotation position and hjust, vjust use to modify it
# xpos ypos        annotateText        hjust    vjust
# 1 -Inf -Inf Bottom Left (h0,v0)        0        0
# 2 -Inf  Inf    Top Left (h0,v1)        0        1
# 3  Inf -Inf  Bottom Right h1,v0        1        0
# 4  Inf  Inf     Top Right h1,v1        1        1

paletteCB6 = c("#CC743C", "#A49C94", "#FCC41C", "#009E73", "#040404", "#0464A4") #"#F4F4FC"
paletteCB6 = c("#2D3142","#058ED9","#848FA2","#E1DAAE","#FF934F", "#CC2D35")

sample_data(biom_ITS_exp1_ev)$year <- factor(sample_data(biom_ITS_exp1_ev)$year, levels=c("2015","2014","2010","2005","2000","1995"))
sample_data(biom_16s_exp1_ev)$year <- factor(sample_data(biom_ITS_exp1_ev)$year, levels=c("2015","2014","2010","2005","2000","1995"))

set.seed(2020)
nmds_ITS_exp1_ev <- ordinate(biom_ITS_exp1_ev, method ="NMDS", distance="bray", autotransform = FALSE)
nmds_ITS_exp1_ev$stress

nmds_16s_exp1_ev <- ordinate(biom_16s_exp1_ev, method ="NMDS", distance="bray", autotransform = FALSE)
nmds_16s_exp1_ev$stress

cap_ITS_exp1_ev <- ordinate(biom_ITS_exp1_ev, ~ year + habitat + Condition(block), method ="CAP", distance="bray")
RsquareAdj(cap_ITS_exp1_ev)
RsquareAdj(fullModel_ITS)
anova(cap_ITS_exp1_ev)

cap_16s_exp1_ev <- ordinate(biom_16s_exp1_ev, ~ year + habitat + year:habitat, method ="CAP", distance="bray")
RsquareAdj(cap_16s_exp1_ev)
RsquareAdj(fullModel_16s_1)
anova(cap_16s_exp1_ev)

PlotOrdin <-function(dataframe, ord){
ord <- plot_ordination(dataframe, ord, color="year", shape="habitat") + 
  geom_point(size=1.7, alpha=0.9, aes(shape=habitat)) +
  scale_shape_manual(values = c(17, 16), labels=c("DF","PS")) +
  scale_colour_manual(values=paletteCB6) +
    theme_classic() +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
    theme(strip.text.x = element_text(size = 9, face = "bold")) +
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 0, face = "bold", size = 8,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 9, face = "bold")) +
  grids(linetype = "dashed") +
  theme(legend.position="right")
return(ord)
}

# *** FIGURE 3 - Ordination ------------------------------------------------------------------------
Fi3_ordin <- ggarrange(PlotOrdin(biom_ITS_exp1_ev, nmds_ITS_exp1_ev) + labs(title = "Fungi") +
             annotate("text", -Inf, -Inf, label = "stress== 0.143", parse=TRUE, size = 3, hjust = -0.1, vjust = -0.8), 
             PlotOrdin(biom_16s_exp1_ev, nmds_16s_exp1_ev) + labs(title = "Prokaryotes") +
               annotate("text", Inf, Inf, label = "stress== 0.082", parse=TRUE, size = 3, hjust = 1.1, vjust = 1.2), 
             PlotOrdin(biom_ITS_exp1_ev, cap_ITS_exp1_ev) + labs(title = "Fungi") +
               annotate("text", -Inf, Inf, label="italic(R)[Adj]^2 == 0.266", parse=TRUE, size = 3, hjust = -0.1, vjust = 1.2), 
             PlotOrdin(biom_16s_exp1_ev, cap_16s_exp1_ev) + labs(title = "Prokaryotes") +
               annotate("text", Inf, -Inf, label ="italic(R)[Adj]^2 == 0.448", parse=TRUE, size = 3, hjust = 1.1, vjust = -0.2), 
             labels = c("A", "B", "C", "D"),
             widths = c(1, 1),
             common.legend = TRUE, legend = "bottom",
             align = "hv", 
             ncol = 2,
             nrow = 2)

Fi3_ordin

# PlotOrdin(biom_ITS_exp1_ev, nmds_ITS_exp1_ev) +
#   annotate("text", -Inf, -Inf, label = "stress== 0.143", parse=TRUE, size = 3, hjust = -0.1, vjust = -0.8) 
# PlotOrdin(biom_16s_exp1_ev, nmds_16s_exp1_ev) +
#   annotate("text", -Inf, -Inf, label = "stress== 0.082", parse=TRUE, size = 3, hjust = -0.1, vjust = -0.8) 
# PlotOrdin(biom_ITS_exp1_ev, cap_ITS_exp1_ev) +
#   annotate("text", Inf, -Inf, label="italic(R)[Adj]^2 == 0.266", parse=TRUE, size = 3, hjust = 1.1, vjust = -0.2) 
# PlotOrdin(biom_16s_exp1_ev, cap_16s_exp1_ev) +
#   annotate("text", Inf, -Inf, label ="italic(R)[Adj]^2 == 0.448", parse=TRUE, size = 3, hjust = 1.1, vjust = -0.2) 


# FIGURE SX - Stress plots -------------------------------------------------------------------------
par(mfrow=c(1,2)) 
stressplot(nmds_ITS_exp1_ev, main="NMDS Fungi")
stressplot(nmds_16s_exp1_ev, main="NMDS Prokaryotes")
dev.off()


# Ordination Experiemnt 2 and 3 -------------------------------------------------------------------- 
# renaming for plotting 
colnames(sample_data(biom_ITS_exp2_ev))[5] <- "block"
colnames(sample_data(biom_16s_exp2_ev))[5] <- "block"

colnames(sample_data(biom_ITS_exp3_ev))[5] <- "block"
head(sample_data(biom_ITS_exp3_ev))
colnames(sample_data(biom_16s_exp3_ev))[5] <- "block"
head(sample_data(biom_16s_exp3_ev))

otu_ITS_exp2_ev <- as.data.frame(otu_table(biom_ITS_exp2_ev))
metadata_ITS_exp2_ev = as(sample_data(biom_ITS_exp2_ev), "data.frame")

otu_ITS_exp3_ev <- as.data.frame(otu_table(biom_ITS_exp3_ev))
metadata_ITS_exp3_ev = as(sample_data(biom_ITS_exp3_ev), "data.frame")

otu_16s_exp2_ev <- as.data.frame(otu_table(biom_16s_exp2_ev))
metadata_16s_exp2_ev = as(sample_data(biom_16s_exp2_ev), "data.frame")

otu_16s_exp3_ev <- as.data.frame(otu_table(biom_16s_exp3_ev))
metadata_16s_exp3_ev = as(sample_data(biom_16s_exp3_ev), "data.frame")

# Experiemnt 2 - this is easy, one factor only
cap_ITS_exp2_ev <- ordinate(biom_ITS_exp2_ev, ~ year, method ="CAP", distance="bray")
RsquareAdj(cap_ITS_exp2_ev)
anova(cap_ITS_exp2_ev)

cap_16s_exp2_ev <- ordinate(biom_16s_exp2_ev, ~ year, method ="CAP", distance="bray")
RsquareAdj(cap_16s_exp2_ev)
anova(cap_16s_exp2_ev)

# Experiment 3 - this have two, site and habitat, no interactions.
# Calculating full model with random factor. I want to evaluate the
# effect of sites so I put habitat as a preudo random factor
cap_ITS_exp3_ev <- ordinate(biom_ITS_exp3_ev, ~ site, method ="CAP", distance="bray")
RsquareAdj(cap_ITS_exp3_ev)
anova(cap_ITS_exp3_ev)
anova(cap_ITS_exp3_ev, by="terms")

cap_16s_exp3_ev <- ordinate(biom_16s_exp3_ev, ~ site, method ="CAP", distance="bray")
RsquareAdj(cap_16s_exp3_ev)
anova(cap_16s_exp3_ev)

PlotOrdin2 <-function(dataframe, ord){
  ord <- plot_ordination(dataframe, ord, color="site", shape="habitat") + 
    geom_point(size=1.7, alpha=0.9, aes(shape=habitat)) +
    scale_shape_manual(values = c(17, 16, 15), labels=c("DF","PS")) +
    scale_colour_manual(values=paletteCB6) +
    theme_classic() +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
    theme(strip.text.x = element_text(size = 9, face = "bold")) +
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 0, face = "bold", size = 8,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 9, face = "bold")) +
    grids(linetype = "dashed") +
    theme(legend.position="right")
  return(ord)
}


sample_data(biom_ITS_exp2_ev)$year <- factor(sample_data(biom_ITS_exp2_ev)$year, levels=c("2015","2010","2005","2000"))
sample_data(biom_16s_exp2_ev)$year <- factor(sample_data(biom_16s_exp2_ev)$year, levels=c("2015","2010","2005","2000"))

sample_data(biom_16s_exp3)
sample_data(biom_16s_exp3_ev)
write.csv(sample_data(biom_16s_exp3), "sample_data_biom_16s_exp3.cvs")
read.csv("sample_data_biom_16s_exp3.cvs", row.names=1, header=T) -> new_map
sample_data(biom_16s_exp3_ev) <- sample_data(new_map)

# *** FIGURE SX - Ordinaiton exp2 and exp3 ---------------------------------------------------------
Fig_S10_ordin <- ggarrange(
                 ggarrange(PlotOrdin(biom_ITS_exp2_ev, cap_ITS_exp2_ev) + labs(title = "Fungi") + 
                             #geom_text_repel(aes(label= X.SampleID), size=2, force = 5) +
                           annotate("text", Inf, Inf, label = "italic(R)[Adj]^2 ==0.616", parse=TRUE, size = 3, hjust = 1.1, vjust = 1.2), 
                           PlotOrdin(biom_16s_exp2_ev, cap_16s_exp2_ev) + labs(title = "Prokaryotes") +
                             #geom_text_repel(aes(label= X.SampleID), size=2, force = 5) +
                           annotate("text", Inf, Inf, label = "italic(R)[Adj]^2 == 0.721", parse=TRUE, size = 3, hjust = 1.1, vjust = 1.2), 
                         labels = c("A", "B"),
                         widths = c(1, 1),
                         common.legend = TRUE, 
                         legend = "bottom",
                         align = "hv", 
                         ncol = 2, 
                         nrow = 1),
                 ggarrange(PlotOrdin2(biom_ITS_exp3_ev, cap_ITS_exp3_ev) + labs(title = "Fungi") +
                             #geom_text_repel(aes(label= X.SampleID), size=2, force = 5) +
                           annotate("text", -Inf, Inf, label="italic(R)[Adj]^2 == 0.901", parse=TRUE, size = 3, hjust = -0.1, vjust = 1.2), 
                           PlotOrdin2(biom_16s_exp3_ev, cap_16s_exp3_ev) + labs(title = "Prokaryotes") +
                             #geom_text_repel(aes(label= X.SampleID), size=2, force = 5) +
                           annotate("text", -Inf, Inf, label ="italic(R)[Adj]^2 == 0.697", parse=TRUE, size = 3, hjust = -0.1, vjust = 1.2), 
                       labels = c("C", "D"),
                       widths = c(1, 1),
                       common.legend = TRUE, 
                       legend = "bottom",
                       align = "hv", 
                       ncol = 2,
                       nrow = 1),
                 ncol = 1,
                 nrow = 2, 
                 widths = c(1, 1),
                 common.legend = FALSE)

Fig_S10_ordin


PlotOrdin2(biom_16s_exp3_ev, cap_16s_exp3_ev) +
  geom_text_repel(aes(label= X.SampleID), size=2, force = 10, )

PlotOrdin2(biom_ITS_exp3_ev, cap_ITS_exp3_ev) +
  geom_text_repel(aes(label= X.SampleID), size=2, force = 10, )



#  PERMANOVA  --------------------------------------------------------------------------------------
library(vegan)
adonis_ITS_exp1_ev <- adonis(t(otu_ITS_exp1_ev) ~ block + year + habitat + year:habitat, data=metadata_ITS_exp1_ev, permutations=999)
adonis_ITS_exp1_ev
round(p.adjust(adonis_ITS_exp1_ev$aov.tab$`Pr(>F)`, "bonferroni"), 4)
round(p.adjust(adonis_ITS_exp1_ev$aov.tab$`Pr(>F)`, "fdr"), 4)

adonis_ITS_exp2_ev <- adonis(t(otu_ITS_exp2_ev) ~ year, data=metadata_ITS_exp2_ev, permutations=999)
adonis_ITS_exp2_ev

adonis_ITS_exp3_ev <- adonis(t(otu_ITS_exp3_ev) ~  habitat + site, data=metadata_ITS_exp3_ev, permutations=999)
adonis_ITS_exp3_ev

adonis_16s_exp1_ev <- adonis(t(otu_16s_exp1_ev) ~ block + year + habitat + year:habitat, data=metadata_16s_exp1_ev, permutations=999)
adonis_16s_exp1_ev
round(p.adjust(adonis_16s_exp1_ev$aov.tab$`Pr(>F)`, "bonferroni"), 4)
round(p.adjust(adonis_16s_exp1_ev$aov.tab$`Pr(>F)`, "fdr"), 4)

adonis_16s_exp2_ev <- adonis(t(otu_16s_exp2_ev) ~ year, data=metadata_16s_exp2_ev, permutations=999)
adonis_16s_exp2_ev

adonis_16s_exp3_ev <- adonis(t(otu_16s_exp3_ev) ~  habitat + site, data=metadata_16s_exp3_ev, permutations=999)
adonis_16s_exp3_ev


# Test the homogenity of group variances -----------------------------------------------------------
permdisp_ITS_exp1_ev_year <- betadisper(vegdist(t(otu_ITS_exp1_ev), method="bray"), metadata_ITS_exp1_ev$year) 
anova(permdisp_ITS_exp1_ev_year, permutations = 999)
permutest(permdisp_ITS_exp1_ev_year, permutations = 999, pairwise = T)
round(p.adjust(anova(permdisp_ITS_exp1_ev_year, permutations = 999)$`Pr(>F)`, "bonferroni"), 4)

permdisp_ITS_exp1_ev_habitat <- betadisper(vegdist(t(otu_ITS_exp1_ev), method="bray"), metadata_ITS_exp1_ev$habitat) 
anova(permdisp_ITS_exp1_ev_habitat , permutations = 999)
permutest(permdisp_ITS_exp1_ev_habitat , permutations = 999, pairwise = T)

permdisp_16s_exp1_ev_year <- betadisper(vegdist(t(otu_16s_exp1_ev), method="bray"), metadata_16s_exp1_ev$year) 
anova(permdisp_16s_exp1_ev_year, permutations = 999)
permutest(permdisp_16s_exp1_ev_year, permutations = 999, pairwise = T)
round(p.adjust(anova(permdisp_16s_exp1_ev_year, permutations = 999)$`Pr(>F)`, "bonferroni"), 4)


permdisp_16s_exp1_ev_habitat <- betadisper(vegdist(t(otu_16s_exp1_ev), method="bray"), metadata_16s_exp1_ev$habitat) 
anova(permdisp_16s_exp1_ev_habitat , permutations = 999)
permutest(permdisp_16s_exp1_ev_habitat , permutations = 999, pairwise = T)


# FIGURE S10 ---------------------------------------------------------------------------------------
PlotBetadisper <- function(betadisp, value, color){
  data.frame(betadisp$group, betadisp$distances) -> df
  colnames(df) <- c("value", "distance")
  betaplot <- ggplot(df, aes(x=value, y=distance)) +
    geom_boxplot(outlier.colour="black", outlier.shape = 8,
                 alpha=0.6, lwd = 0.7, color=color, fill=color) +
    theme_classic() +
    theme(strip.text.x = element_text(size = 9, face = "bold")) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 90, face = "bold", size = 8,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 9, face = "bold")) +
    grids(linetype = "dashed") 
  return(betaplot)
}

# *** FIGURE - betadisper distance form centroids 
ggarrange(
          PlotBetadisper(permdisp_ITS_exp1_ev_year, "year", "#332288") + 
            labs(title = "Fungi\n(Year)", y ="Distance to centroids", x= "Year") +
            annotate("text", Inf, Inf, label = "italic(ANOVA)~P == 0.023", parse=TRUE, size = 2, hjust = 1.2, vjust = 0.9),
          PlotBetadisper(permdisp_ITS_exp1_ev_habitat, "habitat","#332288") + 
            labs(title = "Fungi\n(Habitat)", y ="", x= "Habitat") + 
            annotate("text", Inf, Inf, label = "italic(ANOVA)~P == 0.153", parse=TRUE, size = 2, hjust = 1.2, vjust = 0.9),
          PlotBetadisper(permdisp_16s_exp1_ev_year, "year","#D55E00") + 
            labs(title = "Prokaryotes\n(Year)", y ="", x= "Year") +
            annotate("text", Inf, Inf, label = "italic(ANOVA)~P == 0.01", parse=TRUE, size = 2, hjust = 1.2, vjust = 0.9),
          PlotBetadisper(permdisp_16s_exp1_ev_habitat, "habitat","#D55E00") + 
            labs(title = "Prokaryotes\n(Habitat)", y ="", x= "Habitat") +
            annotate("text", Inf, Inf, label = "italic(ANOVA)~P == 0.276", parse=TRUE, size = 2, hjust = 1.2, vjust = 0.9),
          labels = c("A","B","C","D"),
          widths = c(1,1,1,1),
          align = "hv", 
          ncol = 4, nrow = 1) -> betadisp_plot

betadisp_plot


# >>> HEATMAPS  ------------------------------------------------------------------------------------
# convert to long format using gather and then make it wide with pivot_wider
otu_ITS_exp1_ev_ab <- decostand(otu_ITS_exp1_ev, method="total", MARGIN=2)
head(otu_ITS_exp1_ev_ab)

otu_16s_exp1_ev_ab <- decostand(otu_16s_exp1_ev, method="total", MARGIN=2)
head(otu_16s_exp1_ev_ab)


LongForm <- function(otu, taxa, meta, rank){
            data.frame(OTU = as.factor(row.names(otu)), otu) %>% 
            tidyr::gather(X.SampleID, Abundance, -OTU) %>% 
                  dplyr::left_join(meta[, c("X.SampleID","site","year","habitat")], by = "X.SampleID") %>%
                  dplyr::left_join(taxa, by="OTU") %>%
                  dplyr::select(OTU, X.SampleID, Abundance, year, habitat, Phylum, Class, Genus) -> res_long
                  # adding presence-absence variable
                  ifelse(res_long$Abundance>0, 1, 0) -> res_long$Richness
                  if (rank == "Class"){
                     res_long %>% dplyr::group_by(year, habitat, Class) %>% 
                      dplyr::summarise_if(is.numeric, funs(sum)) %>%
                      as.data.frame() -> res_Class
                       return(res_Class)
                  }else if (rank == "Phylum"){
                      res_long %>% dplyr::group_by(year, habitat, Phylum) %>% 
                        dplyr::summarise_if(is.numeric, funs(sum)) %>%
                        as.data.frame() -> res_Class
                      return(res_Class)
                  }else {
                     res_long %>% group_by(year, habitat, Genus) %>% 
                        dplyr::summarise_if(is.numeric, funs(sum))%>%
                        as.data.frame() -> res_Genus
                    return(res_Genus)
            }
}

LongForm(otu_ITS_exp1_ev_ab, taxa_ITS_exp1_ev, metadata_ITS_exp1_ev, "Class") -> heatdf_ITS_exp1_long_class
LongForm(otu_ITS_exp1_ev_ab, taxa_ITS_exp1_ev, metadata_ITS_exp1_ev, "Genus") -> heatdf_ITS_exp1_long_genus
head(heatdf_ITS_exp1_long_class)
head(heatdf_ITS_exp1_long_genus)

LongForm(otu_16s_exp1_ev_ab, taxa_16s_exp1_ev, metadata_16s_exp1_ev, "Class") -> heatdf_16s_exp1_long_class
head(heatdf_16s_exp1_long_class)

# plot heatmap by class
heatdf_ITS_exp1_long_class$sqrt <- sqrt(heatdf_ITS_exp1_long_class$Abundance*100)
heatdf_ITS_exp1_long_class$year <- factor(heatdf_ITS_exp1_long_class$year, levels=c("2015","2014","2010","2005","2000","1995"))
heatdf_ITS_exp1_long_class$time <- heatdf_ITS_exp1_long_class$year %>% 
recode_factor("2015" = "0", "2014" = "1", "2010" = "5", "2005" = "10", "2000" ="15", "1995"="20")

heatdf_16s_exp1_long_class$sqrt <- sqrt(heatdf_16s_exp1_long_class$Abundance*100)
heatdf_16s_exp1_long_class$year <- factor(heatdf_16s_exp1_long_class$year, levels=c("2015","2014","2010","2005","2000","1995"))
heatdf_16s_exp1_long_class$time <- heatdf_16s_exp1_long_class$year %>% 
  recode_factor("2015" = "0", "2014" = "1", "2010" = "5", "2005" = "10", "2000" ="15", "1995"="20")

# reformat Class taxonomy 
heatdf_ITS_exp1_long_class$Class

heatdf_ITS_exp1_long_class$Class <- heatdf_ITS_exp1_long_class$Class %>% 
    recode_factor('Mortierellomycotina_cls_Incertae_sedis' = "Mortierellomycotina*", 
                  'Pezizomycotina_cls_Incertae_sedis' = "Pezizomycotina*", 
                  'Pucciniomycotina_cls_Incertae_sedis' = "Pucciniomycotina*", 
                  'Rozellomycota_cls_Incertae_sedis' = "Rozellomycota*")

library(stringr)
heatdf_16s_exp1_long_class$Class %>%
  str_replace_all("[[:punct:]]", "") -> heatdf_16s_exp1_long_class$Class

# ranking orders according richness 
heatdf_ITS_exp1_long_class %>%
      dplyr::select(Richness, Class) %>% 
      dplyr::group_by(Class) %>% 
      dplyr::summarise_if(is.numeric, funs(sum)) %>%
      arrange(Richness) %>%
      as.data.frame() -> class_order
rownames(class_order) <- class_order$Class

heatdf_ITS_exp1_long_class$Class <- factor(heatdf_ITS_exp1_long_class$Class, levels=rownames(class_order))
heatdf_ITS_exp1_long_class$sqrt[heatdf_ITS_exp1_long_class$sqrt==0]<- NA
heatdf_ITS_exp1_long_class$Richness[heatdf_ITS_exp1_long_class$Richness==0]<- NA

heatdf_16s_exp1_long_class %>%
  dplyr::select(Richness, Class) %>% 
  dplyr::group_by(Class) %>% 
  dplyr::summarise_if(is.numeric, funs(sum)) %>%
  arrange(Richness) %>%
  as.data.frame() -> class_order_prok
rownames(class_order_prok) <- class_order_prok$Class

heatdf_16s_exp1_long_class$Class <- factor(heatdf_16s_exp1_long_class$Class, levels=rownames(class_order_prok))
heatdf_16s_exp1_long_class$sqrt[heatdf_16s_exp1_long_class$sqrt==0]<- NA
heatdf_16s_exp1_long_class$Richness[heatdf_16s_exp1_long_class$Richness==0]<- NA

# plotting heatmaps function -----------------------------------------------------------------------
PlotHeat <- function(df, rank, type){
plot_heat <- ggplot(df, aes(x=time, y=get(rank), fill=get(type)))+
  geom_tile(aes(height = 0.92, width = 0.92))+ #colour="white",size=0.1
  scale_y_discrete(expand=c(0,0),) +
  facet_grid(~habitat, scales = "free_x", space="free_x") + 
    theme_bw() +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7, color = "black")) +
    theme(strip.text.x = element_text(size = 9, face = "bold")) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5, color = "black")) +
    theme(axis.text.y = element_text(angle = 0, face = "italic", size = 6, hjust = 1, vjust = 0.5, colour = "black")) +
    theme(axis.title = element_text(angle = 0, size = 9, face = "bold")) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.2, "lines"))
return(plot_heat)
}

max(heatdf_ITS_exp1_long_class$sqrt, na.rm = TRUE)
max(heatdf_ITS_exp1_long_class$Richness, na.rm = TRUE)

# *** FIGURE 4 - Heatmap fungi ---------------------------------------------------------------------
library(ggpubr)

brewer.pal(n = 9, name = "RdBu")
display.brewer.pal(n = 9, name = "RdBu")

ggarrange(
PlotHeat(heatdf_ITS_exp1_long_class, "Class", "sqrt") + 
    labs(x="Years",y="Class", title = "Abundance") +
    scale_fill_gradient2(low = "#B2182B", high = "#2166AC", mid = "#F7F7F7", 
                         midpoint = 6, limit = c(0, 12), space = "Lab", na.value="black",
                         name=expression(sqrt(Abundance))) +
    theme(legend.position = "bottom",
          legend.margin=margin(t = 0, unit='cm')) +
    theme(axis.text.y = element_text(face = "plain")),
PlotHeat(heatdf_ITS_exp1_long_class, "Class", "Richness") + 
    labs(x="Years",y="", title = "Richness") +
    scale_fill_gradient2(low = "#B2182B", high = "#2166AC", mid = "#F7F7F7", 
                         midpoint = 150, limit = c(0, 300), space = "Lab", na.value="black",
                         name="Richness") + 
    theme(legend.position = "bottom",
          legend.margin=margin(t = 0, unit='cm')) +
    theme(axis.text.y = element_text(face = "plain")),
          labels = c("A", "B"),
          widths = c(1, 1),
          align = "hv", ncol = 2, nrow = 1) -> Fig_4_heat

Fig_4_heat

max(heatdf_16s_exp1_long_class$sqrt, na.rm = TRUE)
max(heatdf_16s_exp1_long_class$Richness, na.rm = TRUE)

# FIGURE S15 - Heeatmap bacteria -------------------------------------------------------------------
ggarrange(
PlotHeat(heatdf_16s_exp1_long_class, "Class", "sqrt") + 
  labs(x="Year",y="Cass", title = "Abundance") +
  scale_fill_gradient2(low = "#B2182B", high = "#2166AC", mid = "#F7F7F7", 
                       midpoint = 6, limit = c(0, 12), space = "Lab", na.value="white",
                       name=expression(sqrt(Abundance))) +
  theme(legend.position = "bottom",
        legend.margin=margin(t = 0, unit='cm')) +
  theme(axis.text.y = element_text(face = "plain")),
PlotHeat(heatdf_16s_exp1_long_class, "Class", "Richness") + 
  labs(x="Year",y="", title = "Richness") +
  scale_fill_gradient2(low = "#B2182B", high = "#2166AC", mid = "#F7F7F7", 
                       midpoint = 325, limit = c(0, 650), space = "Lab", na.value="white",
                       name="Richness") + 
  theme(legend.position = "bottom",
        legend.margin=margin(t = 0, unit='cm')) +
  theme(axis.text.y = element_text(face = "plain")),
labels = c("A", "B"),
widths = c(1, 1),
align = "hv", ncol = 2, nrow = 1) -> Fig_4_heat_prok

Fig_4_heat_prok

# IMPORTING GUILDS DATA ---------------------------------------------------------------------------
# >>> Mycorrhizal ---------------------------------------------------------------------------------
fungi_symb <- read.csv("fungi_symbiotrophs.csv", header=T, row.names =1)
unique(fungi_symb$Genus)
fungi_symb

# >>> Plant Endophytes ----------------------------------------------------------------------------
fungi_endo <- read.csv("fungi_endophytes.csv", header=T, row.names =1)
unique(fungi_endo$Genus)
fungi_endo

# >>> Plant Pathogen ------------------------------------------------------------------------------
fungi_patho <- read.csv("fungi_pathogen.csv", header=T, row.names =1)
unique(fungi_patho$Genus)
fungi_patho


# Plot heatmap by Guilds --------------------------------------------------------------------------
heatdf_ITS_exp1_long_genus[heatdf_ITS_exp1_long_genus$Genus%in%fungi_symb$Genus, ] -> heatdf_ITS_exp1_long_genus_symb
heatdf_ITS_exp1_long_genus[heatdf_ITS_exp1_long_genus$Genus%in%fungi_patho$Genus, ] -> heatdf_ITS_exp1_long_genus_patho
heatdf_ITS_exp1_long_genus[heatdf_ITS_exp1_long_genus$Genus%in%fungi_endo$Genus, ] -> heatdf_ITS_exp1_long_genus_endo


ReformHeatDf <- function(heat_df){
    as.factor(heat_df$Genus %>% str_replace(" sp.", "")) -> heat_df$Genus
    heat_df$sqrt <- sqrt(heat_df$Abundance*100)
    heat_df$year <- factor(heat_df$year, levels=c("2015","2014","2010","2005","2000","1995"))
    heat_df$time <- heat_df$year %>% 
    recode_factor("2015" = "0", "2014" = "1", "2010" = "5", "2005" = "10", "2000" ="15", "1995"="20")
    # ranking orders according richness 
    heat_df %>%
          dplyr::select(Richness, Genus) %>% 
          dplyr::group_by(Genus) %>% 
          dplyr::summarise_if(is.numeric, funs(sum)) %>%
          arrange(Richness) %>%
          as.data.frame() -> genus_order
    rownames(genus_order) <- genus_order$Genus
    heat_df$Genus <- factor(heat_df$Genus, levels=rownames(genus_order))
    heat_df$sqrt[heat_df$sqrt==0]<- NA
    heat_df$Richness[heat_df$Richness==0]<- NA
return(heat_df)
}

ReformHeatDf(heatdf_ITS_exp1_long_genus_symb) -> plot_heat_ITS_exp1_symb
ReformHeatDf(heatdf_ITS_exp1_long_genus_patho) -> plot_heat_ITS_exp1_patho
ReformHeatDf(heatdf_ITS_exp1_long_genus_endo) -> plot_heat_ITS_exp1_endo

str(plot_heat_ITS_exp1_symb)
levels(plot_heat_ITS_exp1_symb$Genus)

# *** FIGURE 5A - mycorrhizal taxa -----------------------------------------------------------------
max(plot_heat_ITS_exp1_symb$sqrt, na.rm = TRUE)
max(plot_heat_ITS_exp1_symb$Richness, na.rm = TRUE)

library(ggpubr)
ggarrange(
PlotHeat(plot_heat_ITS_exp1_symb, "Genus", "sqrt") + 
   labs(x="Year",y="Genus", title = "Abundance") +
  scale_fill_gradient2(low = "#B2182B", high = "#2166AC", mid = "#F7F7F7", 
  midpoint = 3.25, limit = c(0, 7), space = "Lab", na.value="white",
  name=expression(sqrt(Abundance))) +
  theme(legend.position = "bottom",
  legend.margin=margin(t = 0, unit='cm')),
PlotHeat(plot_heat_ITS_exp1_symb, "Genus", "Richness") + 
   labs(x="Year",y="", title = "Richness") +
  scale_fill_gradient2(low = "#B2182B", high = "#2166AC", mid = "#F7F7F7", 
  midpoint = 10, limit = c(0, 20), space = "Lab", na.value="white",
  name="Richness") +
  theme(legend.position = "bottom",
  legend.margin=margin(t = 0, unit='cm')),
          labels = c("A", "B"),
          widths = c(1, 1),
          align = "h", ncol = 2, nrow = 1) -> Fig_5_heat_symb
 
Fig_5_heat_symb

max(plot_heat_ITS_exp1_patho$sqrt, na.rm = TRUE)
max(plot_heat_ITS_exp1_patho$Richness, na.rm = TRUE)

ggarrange(
PlotHeat(plot_heat_ITS_exp1_patho,"Genus", "sqrt") + 
   labs(x="Year",y="Genus", title = "Abundance") +
  scale_fill_gradient2(low = "#B2182B", high = "#2166AC", mid = "#F7F7F7", 
  midpoint = 2.5, limit = c(0, 5), space = "Lab", na.value="white",
  name=expression(sqrt(Abundance))) +
  theme(legend.position = "bottom",
  legend.margin=margin(t = 0, unit='cm')), 
PlotHeat(plot_heat_ITS_exp1_patho, "Genus", "Richness") + 
   labs(x="Year",y="", title = "Richness") +
  scale_fill_gradient2(low = "#B2182B", high = "#2166AC", mid = "#F7F7F7", 
  midpoint = 7.5, limit = c(0, 15), space = "Lab", na.value="white",
  name="Richness") +
  theme(legend.position = "bottom",
  legend.margin=margin(t = 0, unit='cm')),
          labels = c("A", "B"),
          widths = c(1, 1),
          align = "h", ncol = 2, nrow = 1) -> Fig_5_heat_patho

Fig_5_heat_patho


max(plot_heat_ITS_exp1_endo$sqrt, na.rm = TRUE)
max(plot_heat_ITS_exp1_endo$Richness, na.rm = TRUE)

ggarrange(
PlotHeat(plot_heat_ITS_exp1_endo,"Genus", "sqrt") + 
   labs(x="Year",y="Genus", title = "Abundance") +
  scale_fill_gradient2(low = "#B2182B", high = "#2166AC", mid = "#F7F7F7", 
  midpoint = 3.75, limit = c(0, 7.5), space = "Lab", na.value="white",
  name=expression(sqrt(Abundance))) +
  theme(legend.position = "bottom",
  legend.margin=margin(t = 0, unit='cm')), 

PlotHeat(plot_heat_ITS_exp1_endo, "Genus", "Richness") + 
   labs(x="Year",y="", title = "Richness") +
  scale_fill_gradient2(low = "#B2182B", high = "#2166AC", mid = "#F7F7F7", 
  midpoint = 35.5, limit = c(0, 71), space = "Lab", na.value="white",
  name="Richness") +
  theme(legend.position = "bottom",
  legend.margin=margin(t = 0, unit='cm')),
          labels = c("A", "B"),
          widths = c(1, 1),
          align = "h", ncol = 2, nrow = 1) -> Fig_5_heat_endo

Fig_5_heat_endo


# *** FIGURE 5A, B, C - fungal guilds --------------------------------------------------------------
myco_ab <- PlotHeat(plot_heat_ITS_exp1_symb, "Genus", "sqrt") + 
    labs(x="Years",y="Genus", title = "Abundance") +
    scale_fill_gradient2(low = "#B2182B", high = "#2166AC", mid = "#F7F7F7", 
                         midpoint = 3.25, limit = c(0, 7), space = "Lab", na.value="white",
                         name=expression(sqrt(Abundance))) +
  theme(axis.text.x = element_text(angle = 90, size = 6)) +
  theme(axis.text.y = element_text(size = 5.5)) +
    theme(legend.position = "bottom",
          legend.margin=margin(t = 0, unit='cm'))

myco_rich <- PlotHeat(plot_heat_ITS_exp1_symb, "Genus", "Richness") + 
    labs(x="Years",y="", title = "Richness") +
    scale_fill_gradient2(low = "#B2182B", high = "#2166AC", mid = "#F7F7F7", 
                         midpoint = 10, limit = c(0, 20), space = "Lab", na.value="white",
                         name="Richness") +
  theme(axis.text.x = element_text(angle = 90, size = 6)) +
  theme(axis.text.y = element_text(size = 5.5)) +
    theme(legend.position = "bottom",
          legend.margin=margin(t = 0, unit='cm'))

patho_ab <- PlotHeat(plot_heat_ITS_exp1_patho,"Genus", "sqrt") + 
    labs(x="Years",y="", title = "Abundance") +
    scale_fill_gradient2(low = "#B2182B", high = "#2166AC", mid = "#F7F7F7", 
                         midpoint = 2.5, limit = c(0, 5), space = "Lab", na.value="white",
                         name=expression(sqrt(Abundance))) +
  theme(axis.text.x = element_text(angle = 90, size = 6)) +
  theme(axis.text.y = element_text(size = 5.5)) +
    theme(legend.position = "bottom",
          legend.margin=margin(t = 0, unit='cm'))

patho_rich <- PlotHeat(plot_heat_ITS_exp1_patho, "Genus", "Richness") + 
    labs(x="Years",y="", title = "Richness") +
    scale_fill_gradient2(low = "#B2182B", high = "#2166AC", mid = "#F7F7F7", 
                         midpoint = 7.5, limit = c(0, 15), space = "Lab", na.value="white",
                         name="Richness") +
  theme(axis.text.x = element_text(angle = 90, size = 6)) +
  theme(axis.text.y = element_text(size = 5.5)) +
    theme(legend.position = "bottom",
          legend.margin=margin(t = 0, unit='cm'))

endo_ab <- PlotHeat(plot_heat_ITS_exp1_endo,"Genus", "sqrt") + 
    labs(x="Years",y="", title = "Abundance") +
    scale_fill_gradient2(low = "#B2182B", high = "#2166AC", mid = "#F7F7F7", 
                         midpoint = 3.75, limit = c(0, 7.5), space = "Lab", na.value="white",
                         name=expression(sqrt(Abundance))) +
  theme(axis.text.x = element_text(angle = 90, size = 6)) +
  theme(axis.text.y = element_text(size = 5.5)) +
    theme(legend.position = "right",
          legend.margin=margin(t = 0, unit='cm'))

endo_rich <- PlotHeat(plot_heat_ITS_exp1_endo, "Genus", "Richness") + 
    labs(x="Years",y="", title = "Richness") +
    scale_fill_gradient2(low = "#B2182B", high = "#2166AC", mid = "#F7F7F7", 
                         midpoint = 35.5, limit = c(0, 71), space = "Lab", na.value="white",
                         name="Richness") +
  theme(axis.text.x = element_text(angle = 90, size = 6)) +
  theme(axis.text.y = element_text(size = 5.5)) +
    theme(legend.position = "right",
          legend.margin=margin(t = 0, unit='cm'))


# ggarrange(myco_ab, myco_rich,
#           patho_ab, patho_rich,
#           endo_ab, endo_rich,
#           labels = c("A", "", "B","","C"),
#           widths = c(1,1,1,1,1,1),
#           heights =  c(0.9815, 1,  0.4444, 0.9815, 1, 0.4444),
#           ncol = 2, 
#           nrow = 3) -> Fig_6_guilds
# 
# Fig_6_guilds

# ggdraw() +
#   draw_plot(myco_ab, x = 0, y = 0.01, hjust = 0, vjust = 0, width = 0.2, height = 0.98) +
#   draw_plot(myco_rich, x = 0.18, y = 0, hjust = 0, vjust = 0, width = 0.2, height = 0.98) +
#   draw_plot(patho_ab, x = 0.36, y = 0, hjust = 0, vjust = 0, width = 0.2, height = 1) +
#   draw_plot(patho_rich, x = 0.54, y = 0, hjust = 0, vjust = 0, width = 0.2, height = 1) +
#   draw_plot(endo_ab, x = 0.72, y = 0, hjust = 0, vjust = 0, width = 0.27, height = 0.525) +
#     theme(legend.position = "right") +
#   draw_plot(endo_rich, x = 0.72, y = 0.5,  hjust = 0, vjust = 0,width = 0.27, height = 0.525) +
#     theme(legend.position = "right")
# 
#   draw_plot_label(label = c("A", "", "B","","C"),
#                   size = 15,
#                   x = c(0, 0.5, 0), 
#                   y = c(1, 1, 0.5))

# EXTRACT HEATMAPS TABLES --------------------------------------------------------------------------

# convert to long format using gather and then make it wide with pivot_wider
WideForm <- function(otu, taxa, meta){
    wide_res <- data.frame(OTU = as.factor(row.names(otu)), otu) %>%
      tidyr::gather(X.SampleID, Abundance, -OTU) %>%
      dplyr::left_join(meta[, c("X.SampleID","site","year","habitat")], by = "X.SampleID") %>%
      dplyr::left_join(taxa, by="OTU") %>%
      dplyr::select(OTU, X.SampleID, Abundance, site, year, habitat, Class) %>%
      tidyr::pivot_wider(names_from = year,
                         values_from = Abundance,
                         values_fill = list(Abundance = 0)) %>%
      dplyr::group_by(habitat, Class) %>%
      dplyr::summarise_if(is.numeric, funs(sum)) %>%
      as.data.frame() -> wide_res
    return(wide_res)
}
  
WideForm(otu_ITS_exp1_ev_ab, taxa_ITS_exp1_ev, metadata_ITS_exp1_ev) -> heatdf_ITS_exp1_ev_wide
write.csv(heatdf_ITS_exp1_ev_wide, "heatdf_ITS_exp1_ev_wide.csv")
head(heatdf_ITS_exp1_ev_wide)
dim(heatdf_ITS_exp1_ev_wide)  

# extracting OTU richness per Class, habitat, year -------------------------------------------------
WideFormCount <- function(otu, taxa, meta){
  data.frame(OTU = as.factor(row.names(otu)), otu) %>% 
    tidyr::gather(X.SampleID, Abundance, -OTU) %>% 
    dplyr::left_join(meta[, c("X.SampleID","site","year","habitat","sample")], by = "X.SampleID") %>%
    dplyr::left_join(taxa, by="OTU") %>%
    dplyr::filter(Abundance>0) %>%
    dplyr::group_by(habitat, year, Class) %>%
    dplyr::tally() %>%
    tidyr::pivot_wider(names_from = year,
                       values_from = n,
                       values_fill = list(n = 0)) %>% 
    as.data.frame() -> wide_res_count
  return(wide_res_count)
}

WideFormCount(otu_ITS_exp1_ev, taxa_ITS_exp1_ev, metadata_ITS_exp1_ev) -> heatdf_ITS_exp1_ev_wide_count
head(heatdf_ITS_exp1_ev_wide_count)
write.csv(heatdf_ITS_exp1_ev_wide_count, "heatdf_ITS_exp1_ev_wide_count.csv")


# convert to long format using gather and then make it wide with pivot_wider
WideFormSub <- function(otu, taxa, meta){
  wide_res <- data.frame(OTU = as.factor(row.names(otu)), otu) %>% 
    tidyr::gather(X.SampleID, Abundance, -OTU) %>% 
    dplyr::left_join(meta[, c("X.SampleID","site","year","habitat")], by = "X.SampleID") %>%
    dplyr::left_join(taxa, by="OTU") %>%
    dplyr::select(OTU, X.SampleID, Abundance, site, year, habitat, Genus) %>%
    tidyr::pivot_wider(names_from = year, 
                       values_from = Abundance,
                       values_fill = list(Abundance = 0)) %>%
    dplyr::group_by(habitat, Genus) %>% 
    dplyr::summarise_if(is.numeric, funs(sum)) %>%
    as.data.frame() -> wide_res
  return(wide_res)
}

WideFormSub(otu_ITS_exp1_ev_ab, taxa_ITS_exp1_ev, metadata_ITS_exp1_ev) -> heatdf_ITS_exp1_ev_wide_gen
dim(heatdf_ITS_exp1_ev_wide_gen)

heatdf_ITS_exp1_ev_wide_gen[heatdf_ITS_exp1_ev_wide_gen$Genus%in%fungi_symb$Genus, ] -> heatdf_ITS_exp1_ev_wide_symb
write.csv(heatdf_ITS_exp1_ev_wide_symb, "heatdf_ITS_exp1_ev_wide_symb.csv")

heatdf_ITS_exp1_ev_wide_gen[heatdf_ITS_exp1_ev_wide_gen$Genus%in%fungi_patho$Genus, ] -> heatdf_ITS_exp1_ev_wide_patho
write.csv(heatdf_ITS_exp1_ev_wide_patho, "heatdf_ITS_exp1_ev_wide_patho.csv")

heatdf_ITS_exp1_ev_wide_gen[heatdf_ITS_exp1_ev_wide_gen$Genus%in%fungi_endo$Genus, ] -> heatdf_ITS_exp1_ev_wide_endo
write.csv(heatdf_ITS_exp1_ev_wide_endo, "heatdf_ITS_exp1_ev_wide_endo.csv")


WideFormCountSub <- function(otu, taxa, meta){
  data.frame(OTU = as.factor(row.names(otu)), otu) %>% 
    tidyr::gather(X.SampleID, Richness, -OTU) %>% 
    dplyr::left_join(meta[, c("X.SampleID","site","year","habitat","sample")], by = "X.SampleID") %>%
    dplyr::left_join(taxa, by="OTU") %>%
    dplyr::filter(Richness>0) %>%
    dplyr::group_by(habitat, year, Genus) %>%
    dplyr::tally() %>%
    tidyr::pivot_wider(names_from = year,
                       values_from = n,
                       values_fill = list(n = 0)) %>% 
    as.data.frame() -> wide_res_count
  return(wide_res_count)
}

WideFormCountSub(otu_ITS_exp1_ev, taxa_ITS_exp1_ev, metadata_ITS_exp1_ev) -> heatdf_ITS_exp1_ev_wide_count
head(heatdf_ITS_exp1_ev_wide_count)

heatdf_ITS_exp1_ev_wide_count[heatdf_ITS_exp1_ev_wide_count$Genus%in%fungi_symb$Genus, ] -> heatdf_ITS_exp1_ev_rich_symb
write.csv(heatdf_ITS_exp1_ev_rich_symb, "heatdf_ITS_exp1_ev_rich_symb.csv")

heatdf_ITS_exp1_ev_wide_count[heatdf_ITS_exp1_ev_wide_count$Genus%in%fungi_patho$Genus, ] -> heatdf_ITS_exp1_ev_rich_patho
write.csv(heatdf_ITS_exp1_ev_rich_patho, "heatdf_ITS_exp1_ev_rich_patho.csv")

heatdf_ITS_exp1_ev_wide_count[heatdf_ITS_exp1_ev_wide_count$Genus%in%fungi_endo$Genus, ] -> heatdf_ITS_exp1_ev_rich_endo
write.csv(heatdf_ITS_exp1_ev_rich_endo, "heatdf_ITS_exp1_ev_rich_endo.csv")


# CALCULATE TURNOVER -------------------------------------------------------------------------------
library(codyn)

# reformat the data --------------------------------------------------------------------------------
GenLongTurn <- function(otu_ab, taxa, meta){
data.frame(OTU = as.factor(row.names(otu_ab)), otu_ab) %>% 
  tidyr::gather(X.SampleID, Abundance, -OTU) %>% 
  dplyr::left_join(meta[, c("X.SampleID", "site","year","habitat")], by = "X.SampleID") %>%
  dplyr::left_join(taxa, by="OTU") %>%
  group_by(habitat, year, OTU, Taxonomy) %>% 
  summarise_if(is.numeric, funs(sum)) %>%
  as.data.frame() -> long_form
  long_form$year <- as.numeric(as.character(long_form$year))
return(long_form)
}

GenLongTurn(otu_ITS_exp1_ev_ab, taxa_ITS_exp1_ev, metadata_ITS_exp1_ev) -> turndf_ITS_exp1_long
head(turndf_ITS_exp1_long)

GenLongTurn(otu_16s_exp1_ev_ab, taxa_16s_exp1_ev, metadata_16s_exp1_ev) -> turndf_16s_exp1_long
head(turndf_16s_exp1_long)

# Calculating turnover -----------------------------------------------------------------------------
CalcTurn <- function(df){
    turnover(df = df, 
         time.var = "year", 
         species.var = "Taxonomy", 
         abundance.var = "Abundance", 
         replicate.var = "habitat") -> df_tot
    # calcualting appear
    turnover(df = df, 
         time.var = "year",
         species.var = "Taxonomy",
         abundance.var = "Abundance",
         replicate.var = "habitat",
         metric = "appearance") -> df_app
    # calcualting drops 
    turnover(df = df, 
         time.var = "year",
         species.var = "Taxonomy",
         abundance.var = "Abundance",
         replicate.var = "habitat",
         metric = "disappearance") -> df_dis
    df_tot$metric<-"total"
    names(df_tot)[1]="turnover"
    df_app$metric<-"gain"
    names(df_app)[1]="turnover"
    df_dis$metric<-"drop"
    names(df_dis)[1]="turnover"
    df_all <- rbind(df_tot, df_app, df_dis)
    df_all$years <- df_all$year %>% 
      recode_factor("2015" = "2015-2014", 
                    "2014" = "2014-2010", 
                    "2010" = "2010-2005",
                    "2005" = "2005-2000",
                    "2000" = "2000-1995")
    df_all$time <- df_all$year %>% 
      recode_factor("2015" = "1", 
                    "2014" = "5", 
                    "2010" = "10", 
                    "2005" = "15", 
                    "2000" = "20")
    df_all$time <- as.numeric(as.character(df_all$time))
    return(df_all)
}

CalcTurn(turndf_ITS_exp1_long) -> turn_ITS_exp1_ev
CalcTurn(turndf_16s_exp1_long) -> turn_16s_exp1_ev
head(turn_ITS_exp1_ev)

paletteCB3_fungi = c("#332288","#88CCEE","black")
paletteCB3_prok = c("#D55E00","#DDCC77","black")

# *** FIGURE 9A Turnover fungi ---------------------------------------------------------------------
PlotTurn <- function(dataframe){
ggplot(dataframe, aes(x=time, y=turnover, color=metric)) + 
  labs(x= "", y="Turnover") +
  geom_line(size=0.8) +
  theme_classic() +
  expand_limits(y=c(0,1)) +
  scale_x_continuous(breaks = dataframe$time,
                     labels = paste0(dataframe$time)) + 
  #scale_x_continuous(trans = "reverse") +
    facet_grid(~habitat, scales = "free_x", space="free_x", ) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 8)) +
    theme(strip.text.x = element_text(size = 9, face = "bold")) +
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 9, face = "bold")) +
    grids(linetype = "dashed") +
    #guides(color=guide_legend(nrow=3)) +
    theme(legend.position="bottom") -> plot_turn
return(plot_turn)
}
  
PlotTurn(turn_ITS_exp1_ev) + scale_color_manual(values = paletteCB3_fungi) + labs(title = "Fungi") 
PlotTurn(turn_16s_exp1_ev) + scale_color_manual(values = paletteCB3_prok) + labs(title = "Prokaryotes") 

# Rank shifts --------------------------------------------------------------------------------------
df_ITS_exp1_rankshift <- rank_shift(df=turndf_ITS_exp1_long,
                                    time.var = "year",
                                    species.var = "Taxonomy",
                                    abundance.var = "Abundance", 
                                    replicate.var = "habitat")
df_ITS_exp1_rankshift$time <- df_ITS_exp1_rankshift$year %>% 
  recode_factor("2014-2015" = "1", "2010-2014" = "5","2005-2010" = "10", 
                "2000-2005" = "15", "1995-2000" = "20")
df_ITS_exp1_rankshift$time <- as.numeric(as.character(df_ITS_exp1_rankshift$time))
df_ITS_exp1_rankshift

# separate year_pair for plotting -----------------------------------------------------------------
# df_ITS_exp1_rankshift %>% select(year_pair) %>% separate(year_pair, c(NA, "time")) 

df_16s_exp1_rankshift <- rank_shift(df=turndf_16s_exp1_long,
                                    time.var = "year",
                                    species.var = "Taxonomy",
                                    abundance.var = "Abundance", 
                                    replicate.var = "habitat")
df_16s_exp1_rankshift$time <- df_16s_exp1_rankshift$year %>% 
  recode_factor("2014-2015" = "1", "2010-2014" = "5","2005-2010" = "10", 
                "2000-2005" = "15", "1995-2000" = "20")
df_16s_exp1_rankshift$time <- as.numeric(as.character(df_16s_exp1_rankshift$time))
df_16s_exp1_rankshift$

PlotRank <- function(dataframe){
  ggplot(dataframe, aes(x=time, y=MRS)) + 
    labs(x= "", y="Mean rank shift") +
    geom_line(size=0.8) +
    theme_classic() +
    expand_limits(y=c(0,1)) +
    scale_x_continuous(breaks = dataframe$time,
                       labels = paste0(dataframe$time)) + 
    #scale_x_continuous(trans = "reverse") +
    facet_grid(~habitat, scales = "free_x", space="free_x", ) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
    theme(strip.text.x = element_text(size = 9, face = "bold")) +
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 9, face = "bold")) +
    grids(linetype = "dashed") +
    theme(legend.position="bottom") -> plot_turn
  return(plot_turn)
}

PlotRank(df_ITS_exp1_rankshift) + scale_color_manual(values = "black")
PlotRank(df_16s_exp1_rankshift) + scale_color_manual(values = "black")

# Rate of community change -------------------------------------------------------------------------
df_ITS_exp1_comm_change <- rate_change_interval(turndf_ITS_exp1_long, 
                                                time.var = "year",
                                                species.var = "Taxonomy",
                                                abundance.var = "Abundance",
                                                replicate.var = "habitat")
df_ITS_exp1_comm_change

df_16s_exp1_comm_change <- rate_change_interval(turndf_16s_exp1_long, 
                                                time.var = "year",
                                                species.var = "Taxonomy",
                                                abundance.var = "Abundance",
                                                replicate.var = "habitat")
df_16s_exp1_comm_change

PlotCommChange <- function(dataframe, color){
  dataframe$time <- dataframe$interval %>% 
      recode_factor("5" = "20", "4" = "15", "3" = "10", "2" = "5", "1" ="1")
  dataframe$time <- as.numeric(as.character(dataframe$time))
  ggplot(dataframe, aes(x=time, y=distance)) + 
    labs(x= "Years", y="Euclidean distance") +
    geom_point(alpha = 0.8, shape=16, size=2, color=color) +
    theme_classic() +
    stat_smooth(method = "lm", se = FALSE, size = 1, formula = y ~ x, color ="red") +
    scale_x_continuous(breaks = dataframe$time) +
    #                   labels = paste0(dataframe$time)) + 
    #scale_x_continuous(trans = "reverse") +
    facet_grid(~habitat, scales = "free_x", space="free_x", ) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
    theme(strip.text.x = element_text(size = 9, face = "bold")) +
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 9, face = "bold")) +
    grids(linetype = "dashed") +
    theme(legend.position="bottom") -> plot_turn
  return(plot_turn)
}

PlotCommChange(df_ITS_exp1_comm_change, "#332288") + scale_color_manual(values = "black")
PlotCommChange(df_16s_exp1_comm_change, "#D55E00") + scale_color_manual(values = "black")

# > *** FIGURE 9 Turnover --------------------------------------------------------------------------
library(ggpubr)
library(cowplot)

ggarrange(PlotTurn(turn_ITS_exp1_ev) + scale_color_manual(values = paletteCB3_fungi) + labs(title = "Fungi"), 
          PlotTurn(turn_16s_exp1_ev) + scale_color_manual(values = paletteCB3_prok) + labs(title = "Prokaryotes"),
          PlotRank(df_ITS_exp1_rankshift) + scale_color_manual(values = "black"),
          PlotRank(df_16s_exp1_rankshift) + scale_color_manual(values = "black"),
          PlotCommChange(df_ITS_exp1_comm_change, "#332288") + scale_color_manual(values = "black"),
          PlotCommChange(df_16s_exp1_comm_change, "#D55E00") + scale_color_manual(values = "black"),
          labels = c("A","B","C","D","E","F"),
          heights = c(1.35,0.95,0.95,1.35,0.95,0.95),
          ncol = 2, nrow = 3, 
          align = 'v', 
          common.legend = FALSE,
          legend = "bottom") -> Fig9_turnover

Fig9_turnover

# CULTURABLE FUNGI ---------------------------------------------------------------------------------

read.csv("fungi_cultures.csv", header = TRUE, row.names=1, sep = ",") -> fungi_count
fungi_count

fungi_count$year <- as.factor(fungi_count$year)
fungi_count$years <- fungi_count$year %>% 
  recode_factor("2015" = "0", 
                "2014" = "1", 
                "2010" = "5",
                "2005" = "10",
                "2000" = "15",
                "1995" = "20")

# *** FIGURE S16 ----------------------------------------------------------------------------------
ggplot(fungi_count, aes(x=years, y=counts, color= habitat)) + 
  labs(title = "Fungal isolates", x="Years", y="Number of isolates") +
  geom_boxplot(lwd = 0.7, alpha=0.6) +
  theme_classic() +
  #scale_x_discrete(breaks=c(0,1,5,10,15,20)) +
  theme_classic() +
  scale_color_manual(values = c("#332288","#88CCEE")) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
  theme(strip.text.x = element_text(size = 9, face = "bold")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, face = "bold", size = 8,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 9, face = "bold")) +
  grids(linetype = "dashed") +
  theme(legend.position = c(0.9, 0.85)) -> fig_fungal_cult
  
fig_fungal_cult


# Extracting Wallemiomycetes OTUs rep seq ----------------------------------------------------------

tax_table(subset_taxa(biom_ITS_exp1_ev, Class=="Wallemiomycetes"))


WriteRepSeq <- function(physeq, keepTaxa, filename){
  require(ape)
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[(allTaxa %in% keepTaxa)]
  prune_taxa(myTaxa, physeq) -> prunedTaxa
  write.dna(refseq(prunedTaxa), 
            format="fasta", 
            colsep="", 
            file=filename)
}

WriteRepSeq(biom_ITS_exp1_ev, 
            taxa_names(tax_table(subset_taxa(biom_ITS_exp1_ev, Class=="Wallemiomycetes"))), 
            "Wallemiomycetes_exp1.fasta")

WriteRepSeq(biom_ITS_exp2_ev, 
            taxa_names(tax_table(subset_taxa(biom_ITS_exp1_ev, Class=="Wallemiomycetes"))), 
            "Wallemiomycetes_exp2.fasta")


# ***************************************************-----------------------------------------------
# MINOR REVISIONS ----------------------------------------------------------------------------------
library(scales)

# Reformat taxonomy on the non-rarefied data  
# This is an additional step - not necessary for plotting
source("../R_functions/ReformatTaxonomy.R")

ReformatTaxonomy(biom_ITS_exp1, "ITS") -> biom_ITS_exp1
head(tax_table(biom_ITS_exp1))

tax_table(biom_ITS_exp1)[is.na(tax_table(biom_ITS_exp1))]<-"Unclassified"

otu_ITS_exp1 <- as.data.frame(otu_table(biom_ITS_exp1))
metadata_ITS_exp1 <- as(sample_data(biom_ITS_exp1), "data.frame")
taxa_ITS_exp1 <- as.data.frame(as.matrix(tax_table(biom_ITS_exp1)))
head(taxa_ITS_exp1)

# Rescaling OTU counts to 0-1 ----------------------------------------------------------------------
apply(otu_ITS_exp1, 1, function(x) rescale(x, to = c(0, 1))) -> otu_ITS_exp1_scale
t(otu_ITS_exp1_scale) -> otu_ITS_exp1_scale
head(otu_ITS_exp1_scale)
rowSums(otu_ITS_exp1_scale)
colSums(otu_ITS_exp1_scale)

identical(rownames(otu_ITS_exp1_scale), rownames(taxa_ITS_exp1))
identical(colnames(otu_ITS_exp1_scale), rownames(metadata_ITS_exp1))

# Getting long format
LongForm(otu_ITS_exp1_scale, taxa_ITS_exp1, metadata_ITS_exp1, "Class") -> heatdf_ITS_exp1_scale_class
LongForm(otu_ITS_exp1_scale, taxa_ITS_exp1, metadata_ITS_exp1, "Genus") -> heatdf_ITS_exp1_scale_genus
head(heatdf_ITS_exp1_scale_class)
head(heatdf_ITS_exp1_scale_genus)

heatdf_ITS_exp1_scale_genus$Genus <- gsub(" sp.", "", heatdf_ITS_exp1_scale_genus$Genus)

# preparing to plot the new heatmap ----------------------------------------------------------------

# plot heatmap by class
heatdf_ITS_exp1_scale_class$sqrt <- sqrt(heatdf_ITS_exp1_scale_class$Abundance)
heatdf_ITS_exp1_scale_class$year <- factor(heatdf_ITS_exp1_scale_class$year, levels=c("2015","2014","2010","2005","2000","1995"))
heatdf_ITS_exp1_scale_class$time <- heatdf_ITS_exp1_scale_class$year %>% 
  recode_factor("2015" = "0", "2014" = "1", "2010" = "5", "2005" = "10", "2000" ="15", "1995"="20")

# reformat Class taxonomy 
# heatdf_ITS_exp1_scale_class$Class
# heatdf_ITS_exp1_scale_class$Class <-as.character(heatdf_ITS_exp1_scale_class$Class)
# heatdf_ITS_exp1_scale_class[is.na(heatdf_ITS_exp1_scale_class)] <-"Unclassified"
# str(heatdf_ITS_exp1_scale_class)

heatdf_ITS_exp1_scale_class$Class <- heatdf_ITS_exp1_scale_class$Class %>% 
  recode_factor('Mortierellomycotina_cls_Incertae_sedis' = "Mortierellomycotina*", 
                'Pezizomycotina_cls_Incertae_sedis' = "Pezizomycotina*", 
                'Pucciniomycotina_cls_Incertae_sedis' = "Pucciniomycotina*", 
                'Rozellomycota_cls_Incertae_sedis' = "Rozellomycota*")

# ranking orders according richness 
heatdf_ITS_exp1_scale_class %>%
            dplyr::select(Richness, Class) %>% 
            dplyr::group_by(Class) %>% 
            dplyr::summarise_if(is.numeric, funs(sum)) %>%
                    arrange(Richness) %>%
                          as.data.frame() -> class_order_scaled

rownames(class_order_scaled) <- class_order_scaled$Class
class_order_scaled

#  check variabel range for plotting 
range(heatdf_ITS_exp1_scale_class$Abundance)
range(heatdf_ITS_exp1_scale_class$Richness)

heatdf_ITS_exp1_scale_class$Class <- factor(heatdf_ITS_exp1_scale_class$Class, levels=rownames(class_order_scaled))
heatdf_ITS_exp1_scale_class$sqrt[heatdf_ITS_exp1_scale_class$sqrt==0]<- NA
heatdf_ITS_exp1_scale_class$Richness[heatdf_ITS_exp1_scale_class$Richness==0]<- NA
heatdf_ITS_exp1_scale_class$Abundance[heatdf_ITS_exp1_scale_class$Abundance==0]<- NA
head(heatdf_ITS_exp1_scale_class)

# *** FIGURE 4 new - Heatmap fungi -----------------------------------------------------------------
library(ggpubr)

brewer.pal(n = 9, name = "RdBu")
display.brewer.pal(n = 9, name = "RdBu")

# checking color-gradint scale reange
range(heatdf_ITS_exp1_scale_class$sqrt[!is.na(heatdf_ITS_exp1_scale_class$sqrt)])
range(heatdf_ITS_exp1_scale_class$Richness[!is.na(heatdf_ITS_exp1_scale_class$Richness)])

ggarrange(
  PlotHeat(heatdf_ITS_exp1_scale_class, "Class", "sqrt") + 
    labs(x="Years",y="Class", title = "Abundance") +
    scale_fill_gradient(low = "white", high = "#2166AC",
                        space = "Lab", na.value="black",
                        name=expression(sqrt(Abundance))) + 
    theme(legend.position = "bottom",
          legend.margin=margin(t = 0, unit='cm')) +
    theme(axis.text.y = element_text(face = "plain")) +
    geom_text(aes(x = time, y = Class, label = round(sqrt, 2)), 
              color = "black", size = 1.5),
  
  PlotHeat(heatdf_ITS_exp1_scale_class, "Class", "Richness") + 
    labs(x="Years",y="", title = "Richness") +
    scale_fill_gradient(low = "white", high = "#2166AC", 
                        breaks=c(1, 200, 400),
                         space = "Lab", na.value="black",
                         name="Richness") + 
    theme(legend.position = "bottom",
          legend.margin=margin(t = 0, unit='cm')) +
    theme(axis.text.y = element_text(face = "plain")) +
    geom_text(aes(x = time, y = Class, label = Richness), 
              color = "black", size = 1.5),
  labels = c("A", "B"),
  widths = c(1, 1),
  align = "hv", ncol = 2, nrow = 1) -> Fig_4_heat_new

Fig_4_heat_new

# re-plotting heatmap for 16S data -----------------------------------------------------------------
ReformatTaxonomy(biom_16s_exp1, "16s") -> biom_16s_exp1
tax_table(biom_16s_exp1)[is.na(tax_table(biom_16s_exp1))]<-"Unclassified"

otu_16s_exp1 <- as.data.frame(otu_table(biom_16s_exp1))
metadata_16s_exp1 <- as(sample_data(biom_16s_exp1), "data.frame")
taxa_16s_exp1 <- as.data.frame(as.matrix(tax_table(biom_16s_exp1)))

# Rescaling OTU counts to 0-1 ----------------------------------------------------------------------
apply(otu_16s_exp1, 1, function(x) rescale(x, to = c(0, 1))) -> otu_16s_exp1_scale
t(otu_16s_exp1_scale) -> otu_16s_exp1_scale
head(otu_16s_exp1_scale)
rowSums(otu_16s_exp1_scale)
colSums(otu_16s_exp1_scale)

identical(rownames(otu_16s_exp1_scale), rownames(taxa_16s_exp1))
identical(colnames(otu_16s_exp1_scale), rownames(metadata_16s_exp1))

# get long format 
LongForm(otu_16s_exp1_scale, taxa_16s_exp1, metadata_16s_exp1, "Class") -> heatdf_16s_exp1_scale_class
head(heatdf_16s_exp1_scale_class)

# plot heatmap by class
heatdf_16s_exp1_scale_class$sqrt <- sqrt(heatdf_16s_exp1_scale_class$Abundance)
heatdf_16s_exp1_scale_class$year <- factor(heatdf_16s_exp1_scale_class$year, levels=c("2015","2014","2010","2005","2000","1995"))
heatdf_16s_exp1_scale_class$time <- heatdf_16s_exp1_scale_class$year %>% 
  recode_factor("2015" = "0", "2014" = "1", "2010" = "5", "2005" = "10", "2000" ="15", "1995"="20")

# reformat Class taxonomy 
library(stringr)
heatdf_16s_exp1_scale_class$Class %>%
  str_replace_all("[[:punct:]]", "") -> heatdf_16s_exp1_scale_class$Class

# ranking orders according richness 
heatdf_16s_exp1_scale_class %>%
  dplyr::select(Richness, Class) %>% 
  dplyr::group_by(Class) %>% 
  dplyr::summarise_if(is.numeric, funs(sum)) %>%
  arrange(Richness) %>%
  as.data.frame() -> class_order_prok

rownames(class_order_prok) <- class_order_prok$Class

heatdf_16s_exp1_scale_class$Class <- factor(heatdf_16s_exp1_scale_class$Class, levels=rownames(class_order_prok))
heatdf_16s_exp1_scale_class$sqrt[heatdf_16s_exp1_scale_class$sqrt==0]<- NA
heatdf_16s_exp1_scale_class$Richness[heatdf_16s_exp1_scale_class$Richness==0]<- NA

# checking color-gradint scale reange
range(heatdf_16s_exp1_scale_class$sqrt[!is.na(heatdf_16s_exp1_scale_class$sqrt)])
range(heatdf_16s_exp1_scale_class$Richness[!is.na(heatdf_16s_exp1_scale_class$Richness)])

# *** FIGURE S15 - Heeatmap bacteria ---------------------------------------------------------------
ggarrange(
  PlotHeat(heatdf_16s_exp1_scale_class, "Class", "sqrt") + 
    labs(x="Years",y="Class", title = "Abundance") +
    scale_fill_gradient(low = "white", high = "#2166AC",
                        space = "Lab", na.value="black",
                        name=expression(sqrt(Abundance))) + 
    theme(legend.position = "right",
          legend.margin=margin(t = 0, unit='cm')) +
    theme(axis.text.y = element_text(face = "plain")) +
    geom_text(aes(x = time, y = Class, label = round(sqrt, 2)), 
              color = "black", size = 1.5) +
    theme(axis.text.y = element_text(angle = 0, size = 5)),
  
  PlotHeat(heatdf_16s_exp1_scale_class, "Class", "Richness") + 
    labs(x="Years",y="", title = "Richness") +
    scale_fill_gradient(low = "white", high = "#2166AC", 
                        breaks=c(1, 600, 1200),
                        space = "Lab", na.value="black",
                        name="Richness") + 
    theme(legend.position = "right",
          legend.margin=margin(t = 0, unit='cm')) +
    theme(axis.text.y = element_text(face = "plain")) +
    geom_text(aes(x = time, y = Class, label = Richness), 
              color = "black", size = 1.5) +
    theme(axis.text.y = element_text(angle = 0, size = 5)),
  labels = c("A", "B"),
  widths = c(1, 1),
  align = "hv", ncol = 2, nrow = 1) -> Fig_S15_heat_prok

Fig_S15_heat_prok


# REPLOTTING GUILDS DATA ---------------------------------------------------------------------------
# Plot heatmap by Guilds ---------------------------------------------------------------------------
fungi_symb$Genus <- gsub(" sp.", "", fungi_symb$Genus)
fungi_patho$Genus <- gsub(" sp.", "", fungi_patho$Genus)
fungi_endo$Genus <- gsub(" sp.", "", fungi_endo$Genus)

heatdf_ITS_exp1_scale_genus[heatdf_ITS_exp1_scale_genus$Genus%in%fungi_symb$Genus, ] -> heatdf_ITS_exp1_scale_genus_symb
heatdf_ITS_exp1_scale_genus[heatdf_ITS_exp1_scale_genus$Genus%in%fungi_patho$Genus, ] -> heatdf_ITS_exp1_scale_genus_patho
heatdf_ITS_exp1_scale_genus[heatdf_ITS_exp1_scale_genus$Genus%in%fungi_endo$Genus, ] -> heatdf_ITS_exp1_scale_genus_endo

ReformHeatDf(heatdf_ITS_exp1_scale_genus_symb) -> plot_heat_ITS_exp1_scaled_symb
ReformHeatDf(heatdf_ITS_exp1_scale_genus_patho) -> plot_heat_ITS_exp1_scaled_patho
ReformHeatDf(heatdf_ITS_exp1_scale_genus_endo) -> plot_heat_ITS_exp1_scaled_endo

# *** FIGURE 5A - mycorrhizal taxa -----------------------------------------------------------------
require(ggpubr)

max(plot_heat_ITS_exp1_scaled_symb$sqrt, na.rm = TRUE)
max(plot_heat_ITS_exp1_scaled_symb$Richness, na.rm = TRUE)

ggarrange(
PlotHeat(plot_heat_ITS_exp1_scaled_symb, "Genus", "sqrt") + 
  labs(x="Years",y="Genus", title = "Abundance") +
  scale_fill_gradient(low = "white", high = "#2166AC",
                      space = "Lab", na.value="black",
                      name=expression(sqrt(Abundance))) + 
  theme(legend.position = "bottom",
        legend.margin=margin(t = 0, unit='cm')) +
  theme(axis.text.y = element_text(face = "plain")) +
  geom_text(aes(x = time, y = Genus, label = round(sqrt, 2)), 
            color = "black", size = 1.5) +
  theme(axis.text.y = element_text(face = "italic", size = 6)),

PlotHeat(plot_heat_ITS_exp1_scaled_symb, "Genus", "Richness") + 
  labs(x="Years",y="", title = "Richness") +
  scale_fill_gradient(low = "white", high = "#2166AC", 
                      breaks=c(1, 10, 20),
                      space = "Lab", na.value="black",
                      name="Richness") + 
  theme(legend.position = "bottom",
        legend.margin=margin(t = 0, unit='cm')) +
  theme(axis.text.y = element_text(face = "plain")) +
  geom_text(aes(x = time, y = Genus, label = Richness), 
            color = "black", size = 1.5) +
  theme(axis.text.y = element_text(face = "italic", size = 6)),
labels = c("A", "B"),
widths = c(1, 1),
align = "hv", ncol = 2, nrow = 1)-> Fig_5_heat_scaled_symb

Fig_5_heat_scaled_symb

max(plot_heat_ITS_exp1_scaled_patho$sqrt, na.rm = TRUE)
max(plot_heat_ITS_exp1_scaled_patho$Richness, na.rm = TRUE)

ggarrange(
  PlotHeat(plot_heat_ITS_exp1_scaled_patho, "Genus", "sqrt") + 
    labs(x="Years",y="Genus", title = "Abundance") +
    scale_fill_gradient(low = "white", high = "#2166AC",
                        space = "Lab", na.value="black",
                        name=expression(sqrt(Abundance))) + 
    theme(legend.position = "bottom",
          legend.margin=margin(t = 0, unit='cm')) +
    theme(axis.text.y = element_text(face = "plain")) +
    geom_text(aes(x = time, y = Genus, label = round(sqrt, 2)), 
              color = "black", size = 1.5) +
    theme(axis.text.y = element_text(face = "italic", size = 6)),
  
  PlotHeat(plot_heat_ITS_exp1_scaled_patho, "Genus", "Richness") + 
    labs(x="Years",y="", title = "Richness") +
    scale_fill_gradient(low = "white", high = "#2166AC", 
                        breaks=c(1, 5, 10, 15),
                        space = "Lab", na.value="black",
                        name="Richness") + 
    theme(legend.position = "bottom",
          legend.margin=margin(t = 0, unit='cm')) +
    theme(axis.text.y = element_text(face = "plain")) +
    geom_text(aes(x = time, y = Genus, label = Richness), 
              color = "black", size = 1.5) +
    theme(axis.text.y = element_text(face = "italic", size = 6)),
  labels = c("A", "B"),
  widths = c(1, 1),
  align = "hv", ncol = 2, nrow = 1)-> Fig_5_heat_scaled_patho

Fig_5_heat_scaled_patho

max(plot_heat_ITS_exp1_scaled_endo$sqrt, na.rm = TRUE)
max(plot_heat_ITS_exp1_scaled_endo$Richness, na.rm = TRUE)

ggarrange(
  PlotHeat(plot_heat_ITS_exp1_scaled_endo, "Genus", "sqrt") + 
    labs(x="Years",y="Genus", title = "Abundance") +
    scale_fill_gradient(low = "white", high = "#2166AC",
                        space = "Lab", na.value="black",
                        name=expression(sqrt(Abundance))) + 
    theme(legend.position = "right",
          legend.margin=margin(t = 0, unit='cm')) +
    theme(axis.text.y = element_text(face = "plain")) +
    geom_text(aes(x = time, y = Genus, label = round(sqrt, 2)), 
              color = "black", size = 1.5) +
    theme(axis.text.y = element_text(face = "italic", size = 6)),
  
  PlotHeat(plot_heat_ITS_exp1_scaled_endo, "Genus", "Richness") + 
    labs(x="Years",y="", title = "Richness") +
    scale_fill_gradient(low = "white", high = "#2166AC", 
                        space = "Lab", na.value="black",
                        name="Richness") + 
    theme(legend.position = "right",
          legend.margin=margin(t = 0, unit='cm')) +
    theme(axis.text.y = element_text(face = "plain")) +
    geom_text(aes(x = time, y = Genus, label = Richness), 
              color = "black", size = 1.5) +
    theme(axis.text.y = element_text(face = "italic", size = 6)),
  labels = c("A", "B"),
  widths = c(1, 1),
  align = "hv", 
  ncol = 1, 
  nrow = 2)-> Fig_5_heat_scaled_endo

Fig_5_heat_scaled_endo


# SUPPORTING FIGURE FOR OTU RESCALING --------------------------------------------------------------
# Comparing read number (depth) and rescaled data (0-1) --------------------------------------------
tax_table(subset_taxa(biom_ITS_exp1, Class=="Wallemiomycetes"))
otu_table(subset_taxa(biom_ITS_exp1, Class=="Wallemiomycetes"))
otu_table(subset_taxa(biom_ITS_exp1, Class=="Wallemiomycetes"))
otu_ITS_exp1_scale[taxa_names(subset_taxa(biom_ITS_exp1, Class=="Wallemiomycetes")), ]

plot()

otu_table(subset_taxa(biom_ITS_exp1, Class=="Wallemiomycetes"))[1,] -> pre
otu_ITS_exp1_scale[taxa_names(subset_taxa(biom_ITS_exp1, Class=="Wallemiomycetes")), ][1,] -> post
as.data.frame(t(pre)) -> plot_wall
as.data.frame(post) -> post

identical(rownames(plot_wall), rownames(post))
cbind(plot_wall, post) -> plot_wall
sum(plot_wall)


ggplot(plot_wall, aes(OTU_2, post)) + 
  geom_point() +
  geom_segment(aes(yend=post), xend=0, colour = "red", size = 0.5, linetype="dashed") +
  geom_segment(aes(xend=OTU_2), yend=0, colour = "grey", size = 0.5, linetype="dashed") +
  theme_classic() +
  labs(title = "OTU 2 - Geminibasidium sp." , x="Read number", y= "Rescaled 0-1") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, face = "bold", size = 8,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 9, face = "bold")) -> rescaled_gemini

rescaled_gemini


