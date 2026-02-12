# **** DATA ANALYSIS UPARSE *****-------------------------------------------------------------------
# Project name: Tuber magantum project Giorgio
# Manuscript:   
# Authors:      
# Affiliation:  University of Perugia - Michigan State University
# Journal:     
# Date:         May 24, 2019
# ************************************--------------------------------------------------------------

# WORKING ENVIRONMENT SETUP ------------------------------------------------------------------------
library(styler) # hilight the code and press CTRL+SHIFT+A to style the code.
options(scipen = 9999) #to use decimals
options(max.print=100000000) # to print more lines on the display
# rm(list= ls()) # remove all objects

# main required packages ---------------------------------------------------------------------------
library(phyloseq); packageVersion("phyloseq")
library(Biostrings)
library(ape)
library(ggplot2)
library(vegan); packageVersion("vegan")
library(dplyr)
library(indicspecies)
library(multcompView)
library(ggpubr)

# Color Palettes -----------------------------------------------------------------------------------
palette_site = c("#332288", "#88CCEE", "#117733")
pie(
  rep(1, length(palette_site)),
  labels = sprintf("%d (%s)",
                   seq_along(palette_site), palette_site),
  col = palette_site
)

palette_year = c("#FF899d", "#663F05") #,  "#E69F00",
pie(
  rep(1, length(palette_year)),
  labels = sprintf("%d (%s)",
                   seq_along(palette_year), palette_year),
  col = palette_year
)


# ************************************--------------------------------------------------------------
# IMPORTING DATASETS -------------------------------------------------------------------------------
# A) uparse ITS R1 ---------------------------------------------------------------------------------
# Importing taxonomies at 0.6 confidence -----------------------------------------------------------
taxonomy_ITS06 <-
  read.delim(
    "taxonomies/consensus_taxonomy_ITS06.txt",
    header = TRUE,
    row.names = 1,
    sep = "\t")

head(taxonomy_ITS06)

# Importing taxonomies at 0.8 confidence -----------------------------------------------------------
taxonomy_ITS08 <-
  read.delim(
    "taxonomies/consensus_taxonomy_ITS08.txt",
    header = TRUE,
    row.names = 1, 
    sep = "\t")

head(taxonomy_ITS08)

# check for identical ordering
identical(rownames(taxonomy_ITS06), rownames(taxonomy_ITS08))
# sample_order <- match(rownames(taxonomy_ITS08), rownames(taxonomy_ITS06))
# taxonomy_ITS06 <- taxonomy_ITS06[,sample_order]

taxonomy_ITS08$Kingdom_06 <- taxonomy_ITS06$Kingdom
head(taxonomy_ITS08)
dim(taxonomy_ITS08)

levels(taxonomy_ITS08$Kingdom)

nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom!="Fungi",])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom=="Choanoflagellozoa",])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom=="Metazoa",])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom=="Stramenopila",])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom=="Rhizaria",])

# removing non-fungal OTUs
subset(taxonomy_ITS08, taxonomy_ITS08$Kingdom_06 == "Fungi") -> taxonomy_ITS08_filt
dim(taxonomy_ITS08_filt)
str(taxonomy_ITS08_filt)

# checking filtering results
nrow(taxonomy_ITS08_filt[taxonomy_ITS08_filt$Kingdom != "Fungi", ])
nrow(taxonomy_ITS08_filt[taxonomy_ITS08_filt$Kingdom_06 != "Fungi", ])

# classify all OTUs at 60% identity to Fungi
taxonomy_ITS08_filt$Kingdom <- taxonomy_ITS08_filt$Kingdom_06

# formatting ranks
taxonomy_ITS08_filt <- taxonomy_ITS08_filt[, 1:9]
taxonomy_ITS08_filt$OTU_ID <- rownames(taxonomy_ITS08_filt)

taxonomy_ITS08_filt <- taxonomy_ITS08_filt[, c(
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species",
  "OTU_ID",
  "Isolate",
  "Isolate_percent_id"
)]

head(taxonomy_ITS08_filt)

# importing other tables 
otus_ITS_uparse_R1 <-
  read.delim(
    "otu_tables/otu_table_ITS_UPARSE_R1.txt",
    header = TRUE,
    row.names = 1,)

metadata_ITS_uparse_R1 <-
  read.delim(
    "metadata/map_new_ITS.txt",
    row.names = 1,
    header = TRUE,
    sep = "\t")

otus_seq_ITS_uparse_R1 <-
  readDNAStringSet(
    "fasta/otus_ITS_R1.fasta",
    format = "fasta",
    seek.first.rec = TRUE,
    use.names = TRUE)

physeq_ITS_uparse <- phyloseq(otu_table(otus_ITS_uparse_R1, taxa_are_rows = TRUE),
                              sample_data(metadata_ITS_uparse_R1),
                              tax_table(as.matrix(taxonomy_ITS08_filt)),
                              otus_seq_ITS_uparse_R1) 

physeq_ITS_uparse
str(physeq_ITS_uparse)
head(sample_data(physeq_ITS_uparse))
tax_table(physeq_ITS_uparse)[tax_table(physeq_ITS_uparse)==""]<- NA
head(tax_table(physeq_ITS_uparse))

# checking the phyloseq object
sort(unique(as.data.frame(tax_table(physeq_ITS_uparse))$Kingdom)) # everything is Fungi
nrow(as.data.frame(tax_table(physeq_ITS_uparse))[as.data.frame(tax_table(physeq_ITS_uparse))$Kingdom!="Fungi",])

# ************************************--------------------------------------------------------------
# B) uparse 16S PAIRED -----------------------------------------------------------------------------
# Importing CONSTAX SILVA taxonomy -----------------------------------------------------------------
silva_taxonomy_16s <-
  read.delim("taxonomies/consensus_taxonomy_16s.txt",
             header = TRUE,
             row.names = 1)

head(silva_taxonomy_16s)

# cleaning taxonomy labels -------------------------------------------------------------------------
colnames(silva_taxonomy_16s) <-
  c("Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
    "OTU_ID",
    "Isolate",
    "Isolate_percent_id")

silva_taxonomy_16s$OTU_ID <- rownames(silva_taxonomy_16s)

silva_taxonomy_16s[, "Kingdom"] <- as.factor(gsub("_1", "", silva_taxonomy_16s[, "Kingdom"]))
silva_taxonomy_16s[, "Phylum"] <- as.factor(gsub("_1", "", silva_taxonomy_16s[, "Phylum"]))
silva_taxonomy_16s[, "Class"] <- as.factor(gsub("_1", "", silva_taxonomy_16s[, "Class"]))
silva_taxonomy_16s[, "Order"] <- as.factor(gsub("_1", "",silva_taxonomy_16s[, "Order"]))
silva_taxonomy_16s[, "Family"] <- as.factor(gsub("_1", "", silva_taxonomy_16s[, "Family"]))
silva_taxonomy_16s[, "Genus"] <- as.factor(gsub("_1", "", silva_taxonomy_16s[, "Genus"]))
silva_taxonomy_16s[, "Species"] <- as.factor(gsub("_1", "", silva_taxonomy_16s[, "Species"]))

head(silva_taxonomy_16s)
# silva_taxonomy_16s[1:50,]
# endsWith("aceae", as.character(silva_taxonomy_16s$Family))

any(silva_taxonomy_16s$Kingdom == "Chloroplast" | silva_taxonomy_16s$Kingdom == "Mitochondria")
any(silva_taxonomy_16s$Phylum == "Chloroplast" | silva_taxonomy_16s$Phylum == "Mitochondria")
any(silva_taxonomy_16s$Class == "Chloroplast" | silva_taxonomy_16s$Class == "Mitochondria")
any(silva_taxonomy_16s$Family == "Chloroplast" | silva_taxonomy_16s$Family == "Mitochondria")
any(silva_taxonomy_16s$Genus == "Chloroplast" | silva_taxonomy_16s$Genus == "Mitochondria")

silva_taxonomy_16s %>% dplyr::filter(silva_taxonomy_16s$Family == "Chloroplast")
silva_taxonomy_16s %>% dplyr::filter(silva_taxonomy_16s$Family == "Mitochondria")
silva_taxonomy_16s <- subset(silva_taxonomy_16s, Family != "Mitochondria")

# Check for unclassified OTUs
any(silva_taxonomy_16s$Kingdom == "")

silva_taxonomy_16s %>% dplyr::filter(silva_taxonomy_16s$Kingdom == "")
silva_taxonomy_16s <- subset(silva_taxonomy_16s, Kingdom != "")

# importing other tables
otus_16s_uparse <-
  read.delim("otu_tables/otu_table_16s_UPARSE.txt",
             row.names = 1,
             header = TRUE)

metadata_16s_uparse <-
  read.delim(
    "metadata/map_new_16s.txt",
    row.names = 1,
    header = TRUE,
    sep = "\t"
  )

otus_seq_16s_uparse <-
  readDNAStringSet(
    "fasta/otus_16s.fasta",
    format = "fasta",
    seek.first.rec = TRUE,
    use.names = TRUE
  )


physeq_16s_uparse <- phyloseq(otu_table(otus_16s_uparse, taxa_are_rows = TRUE), 
                              sample_data(metadata_16s_uparse),
                              tax_table(as.matrix(silva_taxonomy_16s)),
                              otus_seq_16s_uparse) 

physeq_16s_uparse
head(sample_data(physeq_16s_uparse))
tax_table(physeq_16s_uparse)[tax_table(physeq_16s_uparse)==""]<- NA
head(tax_table(physeq_16s_uparse))

# ************************************--------------------------------------------------------------
# CREATIN NEW PHYLOSEQ OBJETCS ---------------------------------------------------------------------
physeq_ITS_uparse -> physeq_fungi_uparse
physeq_16s_uparse -> physeq_prok_uparse

# REMOVING SAMPLES non included in this study ------------------------------------------------------
sample_data(physeq_fungi_uparse)
physeq_fungi_uparse
sample_data(physeq_prok_uparse)
physeq_prok_uparse

saveRDS(physeq_fungi_uparse, "physeq_fungi_uparse_prefiltering.rds")
saveRDS(physeq_prok_uparse, "physeq_prok_uparse_prefiltering.rds")


physeq_fungi_uparse <-
  subset_samples(
    physeq_fungi_uparse,
    Description != "2_100_2015" &
      Description != "2_T_2015" &
      Description != "3_100_2015" &
      Description != "3_T_2015" &
      Description != "4_100_2015" &
      Description != "4_T_2015" &
      Description != "5_100_2015" &
      Description != "5_T_2015" &
      Description != "6_100_2015" &
      Description != "6_T_2015"  &
      Description != "7_100_2015" &
      Description != "7_T_2015" &
      Description != "8_100_2015" &
      Description != "8_T_2015"
  )
otu_table(physeq_fungi_uparse) <-
  otu_table(physeq_fungi_uparse)[which(rowSums(otu_table(physeq_fungi_uparse)) > 0), ]

physeq_fungi_uparse

physeq_prok_uparse <-
  subset_samples(
    physeq_prok_uparse,
    Description != "2_100_2015" &
      Description != "2_T_2015" &
      Description != "3_100_2015" &
      Description != "3_T_2015" &
      Description != "4_100_2015" &
      Description != "4_T_2015" &
      Description != "5_100_2015" &
      Description != "5_T_2015" &
      Description != "6_100_2015" &
      Description != "6_T_2015"  &
      Description != "7_100_2015" &
      Description != "7_T_2015" &
      Description != "8_100_2015" &
      Description != "8_T_2015"
  )
otu_table(physeq_prok_uparse) <-
  otu_table(physeq_prok_uparse)[which(rowSums(otu_table(physeq_prok_uparse)) > 0), ]

physeq_prok_uparse

# ************************************--------------------------------------------------------------
# FILTERING OUT CONTAMINANTS -----------------------------------------------------------------------
library(decontam)

# detecting contaminants Fungi
sample_data(physeq_fungi_uparse)$is.neg <-
  sample_data(physeq_fungi_uparse)$Site == "control"

contam_prev_fungi <-
  isContaminant(
    physeq_fungi_uparse,
    method = "prevalence",
    neg = "is.neg",
    threshold = 0.5)

table(contam_prev_fungi$contaminant)

# detecting contaminants Prokaryotes
sample_data(physeq_prok_uparse)$is.neg <-
  sample_data(physeq_prok_uparse)$Site == "control"
contam_prev_prokaryote <-
  isContaminant(
    physeq_prok_uparse,
    method = "prevalence",
    neg = "is.neg",
    threshold = 0.5)

table(contam_prev_prokaryote$contaminant)

# INSPECTING LIBRARY SIZES -------------------------------------------------------------------------
# Before removing controls
df_fungi_uparse <-
  as.data.frame(as.matrix(sample_data(physeq_fungi_uparse)))
df_fungi_uparse$LibrarySize <- sample_sums(physeq_fungi_uparse)
df_fungi_uparse <-
  df_fungi_uparse[order(df_fungi_uparse$LibrarySize), ]
df_fungi_uparse$Index <-
  seq(nrow(df_fungi_uparse)) # sample numbering
df_fungi_uparse

# reorder the factor
df_fungi_uparse$Site <- factor(
  df_fungi_uparse$Site,
  levels = c(
    "Citta_di_Castello",
    "Ripabianca",
    "San_Giovanni_dAsso",
    "control"))

df_prok_uparse <-
  as.data.frame(as.matrix(sample_data(physeq_prok_uparse)))
df_prok_uparse$LibrarySize <- sample_sums(physeq_prok_uparse)
df_prok_uparse <- df_prok_uparse[order(df_prok_uparse$LibrarySize), ]
df_prok_uparse$Index <- seq(nrow(df_prok_uparse)) # sample numbering
df_prok_uparse

df_prok_uparse$Site <- factor(
  df_prok_uparse$Site,
  levels = c(
    "Citta_di_Castello",
    "Ripabianca",
    "San_Giovanni_dAsso",
    "control"))

# *** FIGURES S2 - Distribution of Sample Libraries ------------------------------------------------
ggarrange(
  ggplot(data=df_fungi_uparse, aes(x=Index, y=LibrarySize, color=is.neg)) +
    geom_point(alpha =0.7) +
    labs(title="Fungi", x="Sample number", y="Read number") + 
    theme_classic() +
    scale_colour_manual("Negative control", values = c("grey", "red")) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)),
  ggplot(data=df_prok_uparse, aes(x=Index, y=LibrarySize, color=is.neg)) +
    geom_point(alpha =0.7) +
    labs(title="Prokaryotes", x="Sample number", y="Read number") + 
    theme_classic() +
    scale_colour_manual("Negative control", values = c("grey", "red")) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)),
  ggplot(df_fungi_uparse, aes(x = LibrarySize)) + # Histogram of sample read counts
    geom_histogram(color = "indianred", fill = "indianred", binwidth = 1000) +
    theme_classic() +
    #facet_grid(~Treatment, scales = "free_x", space="free_x") +
    labs(title=NULL, x="Read number", y="Sample number") +
    theme(axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)),
  ggplot(df_prok_uparse, aes(x = LibrarySize)) + # Histogram of sample read counts
    geom_histogram(color = "indianred", fill = "indianred", binwidth = 1000) +
    theme_classic() +
    #facet_grid(~Treatment, scales = "free_x", space="free_x") +
    labs(title=NULL, x="Read number", y="Sample number") + 
    theme(axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)),
  labels = c("A", "B", "C", "D"),
  widths = c(1,1,1,1),
  align = "v", ncol = 2, nrow = 2,
  common.legend = TRUE,
  legend = "bottom") -> lib_size_all

lib_size_all


# removing contaminants form the phyloseq object - Fungi -------------------------------------------
remove_taxa = function(physeq, badTaxa) {
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

contaminants_fungi = rownames(subset(contam_prev_fungi, contaminant %in% c("TRUE")))
physeq_fungi_uparse_clean <-
  remove_taxa(physeq_fungi_uparse, contaminants_fungi)
physeq_fungi_uparse_clean
sample_data(physeq_fungi_uparse_clean)

# Now removing controls samples
physeq_fungi_uparse_clean <-
  subset_samples(physeq_fungi_uparse_clean, Type %in% c("truffle", "soil"))
otu_table(physeq_fungi_uparse_clean) <-
  otu_table(physeq_fungi_uparse_clean)[which(rowSums(otu_table(physeq_fungi_uparse_clean)) > 0), ]
physeq_fungi_uparse_clean
levels(sample_data(physeq_fungi_uparse_clean)$Site)

any(taxa_sums(physeq_fungi_uparse_clean) == 0)
any(sample_sums(physeq_fungi_uparse_clean) == 0)
sort(sample_sums(physeq_fungi_uparse_clean))

# removing contaminants form the phyloseq object - Prokaryotes -------------------------------------
contaminants_prokaryote = rownames(subset(contam_prev_prokaryote, contaminant %in%
                                            c("TRUE")))
length(contaminants_prokaryote)

physeq_prok_uparse_clean <-
  remove_taxa(physeq_prok_uparse, contaminants_prokaryote)
physeq_prok_uparse_clean

# Now removing controls samples
physeq_prok_uparse_clean <-
  subset_samples(physeq_prok_uparse_clean, Type %in% c("truffle", "soil"))
otu_table(physeq_prok_uparse_clean) <-
  otu_table(physeq_prok_uparse_clean)[which(rowSums(otu_table(physeq_prok_uparse_clean)) > 0), ]
physeq_prok_uparse_clean
levels(sample_data(physeq_prok_uparse_clean)$Site)

any(taxa_sums(physeq_prok_uparse_clean) == 0)
any(sample_sums(physeq_prok_uparse_clean) == 0)
sort(sample_sums(physeq_prok_uparse_clean))

# FILTERING DATASTES following quality  ------------------------------------------------------------

# Oliver et al. 2015, PCR errors and tag switching
# Barberan et al. 2012, removing OTUs that appear in less than x samples
# Lindhal et al. 2013, tag switching - that's a good  one!

physeq_fungi_uparse_qc <- physeq_fungi_uparse_clean
otu_table(physeq_fungi_uparse_qc) <-
  otu_table(physeq_fungi_uparse_qc)[which(rowSums(otu_table(physeq_fungi_uparse_qc)) >= 10), ] ### PCR errors
otu_table(physeq_fungi_uparse_qc)[otu_table(physeq_fungi_uparse_qc) < 4] <- 0 ### tag switching
otu_table(physeq_fungi_uparse_qc) <-
  otu_table(physeq_fungi_uparse_qc)[which(rowSums(otu_table(physeq_fungi_uparse_qc)) > 0), ]

physeq_fungi_uparse_qc

physeq_prok_uparse_qc <- physeq_prok_uparse_clean
otu_table(physeq_prok_uparse_qc) <-
  otu_table(physeq_prok_uparse_qc)[which(rowSums(otu_table(physeq_prok_uparse_qc)) >= 10), ]
otu_table(physeq_prok_uparse_qc)[otu_table(physeq_prok_uparse_qc) < 4] <- 0 ### tag switching
otu_table(physeq_prok_uparse_qc) <-
  otu_table(physeq_prok_uparse_qc)[which(rowSums(otu_table(physeq_prok_uparse_qc)) > 0), ]

physeq_prok_uparse_qc

# Change OTU names ---------------------------------------------------------------------------------
ChangeName <- function(physeq, org){
  otu_names<- rownames(otu_table(physeq))
  if (startsWith(otu_names, org)[1] == TRUE) {
    stop("OTUs are already renamed")
  } else {
    otu_names<- paste(org, otu_names, sep="")
    taxa_names(physeq) <- paste(otu_names)
    return(physeq)
  }
}

ChangeName(physeq_fungi_uparse_qc, "F") -> physeq_fungi_uparse_qc
ChangeName(physeq_prok_uparse_qc, "B") -> physeq_prok_uparse_qc

# EXTRACT LAS TAXONOMIC LEVEL ----------------------------------------------------------------------
# thanks to https://rdrr.io/github/jerryzhujian9/ezR/src/R/basic.R

source("../R_functions/ReformatTaxonomy.R")

head(tax_table(physeq_fungi_uparse_qc))
ReformatTaxonomy(physeq_fungi_uparse_qc) -> physeq_fungi_uparse_filt
head(tax_table(physeq_fungi_uparse_filt))

head(tax_table(physeq_prok_uparse_qc))
ReformatTaxonomy(physeq_prok_uparse_qc) -> physeq_prok_uparse_filt
head(tax_table(physeq_prok_uparse_filt))


# RAREFACTION CURVES -------------------------------------------------------------------------------
# Sample40 has zero reads that's why the warning!
# converting to table objetcs ----------------------------------------------------------------------
sample_sums(physeq_fungi_uparse_filt) >0 
subset_samples(physeq_fungi_uparse_filt, sample_sums(physeq_fungi_uparse_filt) >0) -> physeq_fungi_uparse_filt

otu_fungi_up <- as.data.frame(otu_table(physeq_fungi_uparse_filt))
taxa_fungi_up <- as.data.frame(as.matrix(tax_table(physeq_fungi_uparse_filt)))
metadata_fungi_up <- as.data.frame(as.matrix(sample_data(physeq_fungi_uparse_filt)))
dim(otu_fungi_up)
metadata_fungi_up

sample_sums(physeq_prok_uparse_filt) >0 
subset_samples(physeq_prok_uparse_filt, sample_sums(physeq_prok_uparse_filt) >0) -> physeq_prok_uparse_filt

otu_prok_up<- as.data.frame(otu_table(physeq_prok_uparse_filt))
taxa_prok_up <- as.data.frame(as.matrix(tax_table(physeq_prok_uparse_filt)))
metadata_prok_up <- as.data.frame(as.matrix(sample_data(physeq_prok_uparse_filt)))
dim(metadata_prok_up)
metadata_prok_up

# *** FIGURE S3 - Rarefaction curves ---------------------------------------------------------------

# creating and empty character 
curve_colors_fungi <- vector(mode="character", length=ncol(as.data.frame(otu_table(physeq_fungi_uparse_filt))))
curve_colors_fungi[as.data.frame(as.matrix(sample_data(physeq_fungi_uparse_filt)))$Type=="truffle"] <- "gold"
curve_colors_fungi[as.data.frame(as.matrix(sample_data(physeq_fungi_uparse_filt)))$Type=="soil"] <- "brown"

curve_colors_prok <- vector(mode="character", length=ncol(as.data.frame(otu_table(physeq_prok_uparse_filt))))
curve_colors_prok[as.data.frame(as.matrix(sample_data(physeq_prok_uparse_filt)))$Type=="truffle"] <- "gold"
curve_colors_prok[as.data.frame(as.matrix(sample_data(physeq_prok_uparse_filt)))$Type=="soil"] <- "brown"


par(mar = c(3.5, 3.5, 2, 0) + 0.1,mgp = c(2.4, 0.8, 0),las = 0,mfrow = c(1, 2))
# Change the panel layout to 2 x 1
rarecurve(t(otu_fungi_up), col = curve_colors_fungi, label = FALSE, step = 1000,
          main="Fungi", ylab = "Number of OTUs", xlab = "Number of DNA reads")
legend("topright", legend=c("Truffle", "Soil"),
       col=c("gold", "brown"), lty=2, cex=0.8, box.lty=0,lwd=3,bty = "n") 
# if you rarefy you can add at what depth you rarefy here 
#abline(v = min(sample_sums(physeq_fungi_uparse_qc)), col="red", lwd=3, lty=2)
rarecurve(t(otu_prok_up), col = curve_colors_prok, label = FALSE,  step = 1000,
          main="Prokaryotes", ylab = "Number of OTUs", xlab = "Number of DNA reads")
legend("topright", legend=c("Truffle", "Soil"),
       col=c("gold", "brown"), lty=2, cex=0.8, box.lty=0,lwd=3,bty = "n") 
#abline(v = min(sample_sums(physeq_prok_uparse_qc)), col="red", lwd=3, lty=2)
dev.off()

# Save filtered datasets ---------------------------------------------------------------------------
saveRDS(physeq_fungi_uparse_filt, "physeq_fungi_uparse.rds")
saveRDS(physeq_prok_uparse_filt, "physeq_prok_uparse.rds")

# ************************************--------------------------------------------------------------
# CSS TRANSFORMATION -------------------------------------------------------------------------------
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

CSSNorm(physeq_fungi_uparse_filt) -> physeq_fungi_uparse_mSeq
head(otu_table(physeq_fungi_uparse_mSeq))

CSSNorm(physeq_prok_uparse_filt) -> physeq_prok_uparse_mSeq
head(otu_table(physeq_prok_uparse))


# >>> ALPHA DIVERSITY ------------------------------------------------------------------------------
library(multcompView)
library(ggpubr)

# creating an additional variable for filtering
ExtrDataFrame <- function(physeq, var){
  sample_data(physeq)$LocYear <- paste(
    sample_data(physeq)$Site,
    sample_data(physeq)$Year,
    sep = "_")
  if (var == "Site"){
    physeq %>% subset_samples(LocYear%in%c("Citta_di_Castello_2014",
                                           "San_Giovanni_dAsso_2015",
                                           "Ripabianca_2014")) -> physeq_filt
  }else{
    physeq %>% subset_samples(LocYear%in%c("Citta_di_Castello_2014",
                                           "Citta_di_Castello_2015")) -> physeq_filt
  }
  otu_table(physeq_filt) <- otu_table(physeq_filt)[which(rowSums(otu_table(physeq_filt)) > 0),] 
  # adding aplha metrics to metadata file
  df_alpha <- as.data.frame(as.matrix(sample_data(physeq_filt))) 
  df_alpha$Observed <- specnumber(as.data.frame(otu_table(physeq_filt)), MARGIN = 2)
  df_alpha$Shannon <- diversity(as.data.frame(otu_table(physeq_filt)), index="shannon", MARGIN = 2)
  # re-order factor levels
  df_alpha$Distance <- factor(df_alpha$Distance,levels=c("Inside", "0", "40", "100", "200"))
  return(df_alpha)
}


ExtrDataFrame(physeq_fungi_uparse_filt, "Site") -> df_alpha_fungi
df_alpha_fungi
ExtrDataFrame(physeq_prok_uparse_filt, "Site") -> df_alpha_prok
df_alpha_prok

ExtrDataFrame(physeq_fungi_uparse_filt, "Year") -> df_alpha_fungi_yr
df_alpha_fungi_yr
ExtrDataFrame(physeq_prok_uparse_filt, "Year") -> df_alpha_prok_yr
df_alpha_prok_yr


# count variables levels
df_alpha_fungi %>% group_by(Site) %>% count(vars = Distance)
df_alpha_fungi %>% group_by(Distance) %>% count(vars = Site)
df_alpha_prok %>% group_by(Site) %>% count(vars = Distance)
df_alpha_prok %>% group_by(Distance) %>% count(vars = Site)


# comparing Distance across sites using Wilcox.test
# Since we have done many independent tests (one for each taxon), 
# we should correct for multiple comparisions. However if we use p.adj,
# even even Inside form 0 is not significant 
CompSampl <- function(dataframe, formula, level, var){
  require(ggpubr)
  if (var == "Site"){
    compare_means(formula, data = dataframe[dataframe$Site==level,],
                  p.adjust.method = "bonferroni") -> test_CC 
    test_CC <- as.data.frame(test_CC)[,c(2,3,4)] # for p.adj 4 to 5
    test_CC2 <- data.frame(test_CC[,2], test_CC[,1], test_CC[,3])
    colnames(test_CC2) <- c("group1", "group2", "p") # change p to p.adj
    rbind(test_CC, test_CC2) -> test_all
    as.dist(xtabs(test_all[, 3] ~ (test_all[, 2] + test_all[, 1])), diag = TRUE) -> dist_CC
    data.frame(multcompLetters(dist_CC)['Letters']) -> res_CC
    res_CC$sample <- rownames(res_CC)
    res_CC %>% slice(match(c("Inside", "0","40","100","200"), sample)) -> res_CC
    return(res_CC)
  }else{
    compare_means(formula, data = dataframe[dataframe$Distance==level,],
                  p.adjust.method = "bonferroni") -> test_CC 
    test_CC <- as.data.frame(test_CC)[,c(2,3,4)] 
    test_CC2 <- data.frame(test_CC[,2], test_CC[,1], test_CC[,3])
    colnames(test_CC2) <- c("group1", "group2", "p") 
    rbind(test_CC, test_CC2) -> test_all
    as.dist(xtabs(test_all[, 3] ~ (test_all[, 2] + test_all[, 1])), diag = TRUE) -> dist_CC
    data.frame(multcompLetters(dist_CC)['Letters']) -> res_CC
    res_CC$sample <- rownames(res_CC)
    res_CC %>% slice(match(c("2014", "2015"), sample)) -> res_CC
    return(res_CC)
  }
}


# fungi --------------------------------------------------------------------------------------------
CompSampl(df_alpha_fungi, formula(Observed ~ Distance), "Citta_di_Castello", "Site") -> Wtest_fungi_Cast_rich
CompSampl(df_alpha_fungi, formula(Observed ~ Distance), "Ripabianca", "Site") -> Wtest_fungi_Rip_rich
CompSampl(df_alpha_fungi, formula(Observed ~ Distance), "San_Giovanni_dAsso", "Site") -> Wtest_fungi_Gio_rich

my_lables_fungi_rich <-c(as.character(Wtest_fungi_Cast_rich$Letters), 
                         as.character(Wtest_fungi_Rip_rich$Letters),
                         as.character(Wtest_fungi_Gio_rich$Letters))
my_lables_fungi_rich

CompSampl(df_alpha_fungi, formula(Shannon ~ Distance), "Citta_di_Castello", "Site") -> Wtest_fungi_Cast_shan
CompSampl(df_alpha_fungi, formula(Shannon ~ Distance), "Ripabianca", "Site") -> Wtest_fungi_Rip_shan
CompSampl(df_alpha_fungi, formula(Shannon ~ Distance), "San_Giovanni_dAsso", "Site") -> Wtest_fungi_Gio_shan

my_lables_fungi_shan <-c(as.character(Wtest_fungi_Cast_shan$Letters), 
                         as.character(Wtest_fungi_Rip_shan$Letters),
                         as.character(Wtest_fungi_Gio_shan$Letters))
my_lables_fungi_shan

# prok ---------------------------------------------------------------------------------------------
CompSampl(df_alpha_prok, formula(Observed ~ Distance), "Citta_di_Castello", "Site") -> Wtest_prok_Cast_rich
CompSampl(df_alpha_prok, formula(Observed ~ Distance), "Ripabianca", "Site") -> Wtest_prok_Rip_rich
CompSampl(df_alpha_prok, formula(Observed ~ Distance), "San_Giovanni_dAsso", "Site") -> Wtest_prok_Gio_rich

my_lables_prok_rich <-c(as.character(Wtest_prok_Cast_rich$Letters), 
                        as.character(Wtest_prok_Rip_rich$Letters),
                        as.character(Wtest_prok_Gio_rich$Letters))
my_lables_prok_rich


CompSampl(df_alpha_prok, formula(Shannon ~ Distance), "Citta_di_Castello", "Site") -> Wtest_prok_Cast_shan
CompSampl(df_alpha_prok, formula(Shannon ~ Distance), "Ripabianca", "Site") -> Wtest_prok_Rip_shan
CompSampl(df_alpha_prok, formula(Shannon ~ Distance), "San_Giovanni_dAsso", "Site") -> Wtest_prok_Gio_shan

my_lables_prok_shan <-c(as.character(Wtest_prok_Cast_shan$Letters), 
                        as.character(Wtest_prok_Rip_shan$Letters),
                        as.character(Wtest_prok_Gio_shan$Letters))
my_lables_prok_shan


# Plot Richness ------------------------------------------------------------------------------------
PlotRich <- function(dataframe, var, palette, my_labels){
  #calculating where to put the letters
  max(dataframe[,var] + 0.1* max(dataframe[,var])) -> labels_y
  ggplot(dataframe, aes(x = Distance, y = get(var), color = Site)) +
    geom_boxplot(outlier.colour="black", outlier.size=1.2, outlier.shape = 8, lwd=0.5, width = 0.7) +
    stat_summary(fun=mean, geom="point", shape=18, size=1.5, color="red", fill="red") +
    stat_summary(geom = 'text', label = my_labels, fun= max, aes(y = labels_y), size=2.2, color="black") +
    facet_grid(~Site) +
    #         facet_grid(~Site, scales = "free_x", 
    #          labeller = as_labeller(c(Citta_di_Castello = "Site 1",
    #                                   Ripabianca = "Site 2",
    #                                   San_Giovanni_dAsso = "Site 3"))) +
    theme_classic() +
    #scale_colour_manual(values = palette, 
    #                    labels = c("Site 1", "Site 2", "Site 3")) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) +
    theme(strip.text.x = element_text(size = 8)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 8)) +
    theme(axis.text.x = element_text(angle = 45, size = 7, hjust = 0.9, vjust = 1)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) + 
    grids(linetype = "dashed") +
    theme(legend.position="bottom") -> rich_plot 
  return(rich_plot)
}

PlotRich(df_alpha_fungi, "Observed", palette_site, my_lables_fungi_rich) + 
  labs(x=NULL, y="Richness", title = "Fungi") -> plot_rich_fungi_rich
plot_rich_fungi_rich
PlotRich(df_alpha_prok, "Observed", palette_site, my_lables_prok_rich) + 
  labs(x=NULL, y=NULL, title = "Prokaryotes")  -> plot_rich_prok_rich
plot_rich_prok_rich


PlotRich(df_alpha_fungi, "Shannon", palette_site, my_lables_fungi_rich) + 
  labs(x="Distance", y="Shannon", title = NULL) -> plot_rich_fungi_shan
plot_rich_fungi_shan
PlotRich(df_alpha_prok,"Shannon", palette_site, my_lables_prok_rich) + 
  labs(x="Distance", y=NULL, title = NULL)  -> plot_rich_prok_shan
plot_rich_prok_shan


ggarrange(
  df_alpha_fungi %>%
    subset(Site%in%c("Citta_di_Castello")) %>%
    ggplot(aes(Distance, Observed)) +
    geom_boxplot(),
  
  df_alpha_fungi %>%
    subset(Site%in%c("Ripabianca")) %>%
    ggplot(aes(Distance, Observed)) +
    geom_boxplot(),
  
  df_alpha_fungi %>%
    subset(Site%in%c("San_Giovanni_dAsso")) %>%
    ggplot(aes(Distance, Observed)) +
    geom_boxplot(),
  ncol = 3)


# *** FIGURE 1 - Richness --------------------------------------------------------------------------
Fig_1A <- ggarrange(plot_rich_fungi_rich, 
                    plot_rich_prok_rich,
                    plot_rich_fungi_shan, 
                    plot_rich_prok_shan,
                    labels = c("A", "B", "C", "D"),
                    widths = c(1,1,1,1),
                    align = "hv", 
                    ncol = 2, 
                    nrow = 2,
                    common.legend = TRUE,
                    legend = "none")

Fig_1A

title1=text_grob("Site",size=12, face=2)
grid.arrange(Fig_1A, top = title1) -> Fig_1A_alpha_site
Fig_1A_alpha_site


# now the plotting between years
Calclabel <- function(df, formula){
  CompSampl(df, formula(Observed ~ Year), "Inside", "Year") -> pval_Inside
  CompSampl(df, formula(Observed ~ Year), "0", "Year") -> pval_0
  CompSampl(df, formula(Observed ~ Year), "40", "Year") -> pval_40
  CompSampl(df, formula(Observed ~ Year), "100", "Year") -> pval_100
  CompSampl(df, formula(Observed ~ Year), "200", "Year") -> pval_200
  my_lables <-c(as.character(pval_Inside$Letters), 
                as.character(pval_0$Letters),
                as.character(pval_40$Letters),
                as.character(pval_100$Letters),
                as.character(pval_200$Letters))
  return(my_lables)
}


Calclabel(df_alpha_fungi_yr, formula(Observed ~ Year)) -> Wtest_fungi_year_rich
Calclabel(df_alpha_fungi_yr, formula(Shannon ~ Year)) -> Wtest_fungi_year_shan

Calclabel(df_alpha_prok_yr, formula(Observed ~ Year)) -> Wtest_prok_year_rich
Calclabel(df_alpha_prok_yr, formula(Shannon ~ Year)) -> Wtest_prok_year_shan



#Plotting
PlotRichY <- function(dataframe, var, palette, my_labels){
  #calculating where to put the letters
  max(dataframe[,var] + 0.1* max(dataframe[,var])) -> labels_y
  ggplot(dataframe, aes(x = Year, y = get(var), color = Year)) +
    geom_boxplot(outlier.colour="black", outlier.size=1.2, outlier.shape = 8, lwd=0.5, width = 0.7) +
    stat_summary(fun=mean, geom="point", shape=18, size=1.5, color="red", fill="red") +
    stat_summary(geom = 'text', label = my_labels, fun= max, aes(y = labels_y), size=2.2, color="black") +
    facet_grid(~factor(dataframe$Distance, levels = c("Inside", "0", "40", "100", "200"))) +
    theme_classic() +
    scale_colour_manual(values = palette, 
                        labels = c("Inside", "0", "40", "100", "200")) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) +
    theme(strip.text.x = element_text(size = 8)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 8)) +
    theme(axis.text.x = element_text(angle = 45, size = 7, hjust = 0.9, vjust = 1)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) + 
    grids(linetype = "dashed") +
    theme(legend.position="none") -> rich_plot 
  return(rich_plot)
}

PlotRichY(df_alpha_fungi_yr, "Observed", palette_year, Wtest_fungi_year_rich) + 
  labs(x=NULL, y="Richness", title = "Fungi") 

# **** FIGURE 2 richness year --------------------------------------------------------------------- 
Fig_2 <- ggarrange(PlotRichY(df_alpha_fungi_yr, "Observed", palette_year, Wtest_fungi_year_rich) + 
                      labs(x=NULL, y="Richness", title = "Fungi"),
                    PlotRichY(df_alpha_prok_yr, "Observed", palette_year, Wtest_prok_year_rich) + 
                      labs(x=NULL, y=NULL, title = "Prokaryotes"),
                    PlotRichY(df_alpha_fungi_yr, "Shannon", palette_year, Wtest_fungi_year_shan) + 
                      labs(x="Distance", y="Shannon", title = NULL),
                    PlotRichY(df_alpha_prok_yr,"Shannon", palette_year, Wtest_prok_year_shan) + 
                      labs(x="Distance", y=NULL, title = NULL),
                    labels = c("A", "B", "C", "D"),
                    widths = c(1,1,1,1),
                    align = "hv", 
                    ncol = 2, 
                    nrow = 2,
                    common.legend = TRUE,
                    legend = "none")

Fig_2

title2=text_grob("Year",size=12, face=2)
grid.arrange(Fig_2, top = title2) -> Fig_2_year
Fig_2_year


# ************************************--------------------------------------------------------------
# >>> BETA DIVERSITY -------------------------------------------------------------------------------
library(ggpubr)

# >>> PERMANOVA ------------------------------------------------------------------------------------
count(as.data.frame(as.matrix(
  sample_data(
    physeq_fungi_uparse_mSeq %>%
      subset_samples(
        Type == "soil" & LocYear %in% c(
          "Citta_di_Castello_2014",
          "San_Giovanni_dAsso_2015",
          "Ripabianca_2014"))))), Site)

count(as.data.frame(as.matrix(
  sample_data(
    physeq_fungi_uparse_mSeq %>%
      subset_samples(
        Type == "soil" & LocYear %in% c(
          "Citta_di_Castello_2014",
          "San_Giovanni_dAsso_2015",
          "Ripabianca_2014"))))), Distance)

ExtractOTU <- function(physeq){
  physeq %>%
    subset_samples(Type=="soil" & LocYear%in%c("Citta_di_Castello_2014",
                                               "San_Giovanni_dAsso_2015",
                                               "Ripabianca_2014")) -> Physeq_filt
  otu <- as.data.frame(otu_table(Physeq_filt))
  return(otu)
}

ExtractMETA <- function(physeq){
  physeq %>%
    subset_samples(Type=="soil" & LocYear%in%c("Citta_di_Castello_2014",
                                               "San_Giovanni_dAsso_2015",
                                               "Ripabianca_2014")) -> Physeq_filt
  meta <- as(sample_data(Physeq_filt), "data.frame")
  return(meta)
}

ExtractOTU(physeq_fungi_uparse_mSeq) -> otu_fungi_uparse_mSeq
ExtractMETA(physeq_fungi_uparse_mSeq) -> metadata_fungi_uparse_mSeq

ExtractOTU(physeq_prok_uparse_mSeq) -> otu_prok_uparse_mSeq
ExtractMETA(physeq_prok_uparse_mSeq) -> metadata_prok_uparse_mSeq

# test the effect Site 
adonis_fungi <- adonis(t(otu_fungi_uparse_mSeq) ~ Site,
                       metadata_fungi_uparse_mSeq, 
                       strata=metadata_fungi_uparse_mSeq$Distance, 
                       permutations=9999)
adonis_fungi

# test for the effect of Distance 
adonis_fungi1 <- adonis(t(otu_fungi_uparse_mSeq) ~ Distance,  
                        metadata_fungi_uparse_mSeq, 
                        strata=metadata_fungi_uparse_mSeq$Site, 
                        permutations=9999)
adonis_fungi1

# test for interactions - in these case I do no include the strata
# since Distance which is at the beginning of the model is treated as a random
# effect 
adonis_fungi2 <- adonis(t(otu_fungi_uparse_mSeq) ~ Distance * Site,  
                        metadata_fungi_uparse_mSeq, 
                        permutations=9999)
adonis_fungi2

adonis_fungi3 <- adonis(t(otu_fungi_uparse_mSeq) ~ Site * Distance,
                        metadata_fungi_uparse_mSeq, 
                        permutations=9999)
adonis_fungi3
round(p.adjust(adonis_fungi3$aov.tab$`Pr(>F)`, "bonferroni"), 4)
# [1] 0.0027 1.0000 1.0000     NA     NA


# *** TABLE 1 - PERMANOVA --------------------------------------------------------------------------

# Fungi
#                   Df  SumsOfSqs  MeanSqs  F.Model      R2       Pr(>F)    
#Site           2    2.2735    1.13677  3.02143   0.12183     0.0001 ***
#Distance           3    1.1452    0.38173  1.01462   0.06137     0.3989    
#Site:Distance  6    2.0745    0.34574  0.91896   0.11116     0.9213    
#Residuals         35    13.1682   0.37623            0.70564           
#Total             46    18.6614                      1.00000 


# test the effect Site 
adonis_prok <- adonis(t(otu_prok_uparse_mSeq) ~ Site,
                      metadata_prok_uparse_mSeq, 
                      strata=metadata_prok_uparse_mSeq$Distance, 
                      permutations=9999)
adonis_prok

# test for the effect of Distance 
adonis_prok1 <- adonis(t(otu_prok_uparse_mSeq) ~ Distance,  
                       metadata_prok_uparse_mSeq, 
                       strata=metadata_prok_uparse_mSeq$Site, 
                       permutations=9999)
adonis_prok1

# test for interactions - in these case I do no include the strata
# since Distance which is at the beginning of the model is treated as a random
# effect 
adonis_prok2 <- adonis(t(otu_prok_uparse_mSeq) ~ Distance * Site,  
                       metadata_prok_uparse_mSeq, 
                       permutations=9999)
adonis_prok2

adonis_prok3 <- adonis(t(otu_prok_uparse_mSeq) ~ Site * Distance,
                       metadata_prok_uparse_mSeq, 
                       permutations=9999)
adonis_prok3
round(p.adjust(adonis_prok3$aov.tab$`Pr(>F)`, "bonferroni"), 4)
# [1] 0.0003 0.0891 1.0000     NA     NA

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Site           2    1.5820 0.79099  7.7525 0.24844 0.0001 ***
#   Distance       3    0.5178 0.17260  1.6916 0.08132 0.0077 ** 
#   Site:Distance  6    0.5949 0.09914  0.9717 0.09342 0.5378    
# Residuals     36    3.6731 0.10203         0.57683           
# Total         47    6.3677                 1.00000      
# 


# Now testing the years ----------------------------------------------------------------------------
count(as.data.frame(as.matrix(
  sample_data(
    physeq_fungi_uparse_mSeq %>%
      subset_samples(
        Type == "soil" & LocYear %in% c(
          "Citta_di_Castello_2014",
          "Citta_di_Castello_2015"))))), Year)

ExtractOTU_y <- function(physeq){
  physeq %>%
    subset_samples(Type=="soil" & LocYear%in%c("Citta_di_Castello_2014",
                                               "Citta_di_Castello_2015")) -> Physeq_filt
  otu <- as.data.frame(otu_table(Physeq_filt))
  return(otu)
}

ExtractMETA_y <- function(physeq){
  physeq %>%
    subset_samples(Type=="soil" & LocYear%in%c("Citta_di_Castello_2014",
                                               "Citta_di_Castello_2015")) -> Physeq_filt
  meta <- as(sample_data(Physeq_filt), "data.frame")
  return(meta)
}


ExtractOTU_y(physeq_fungi_uparse_mSeq) -> otu_fungi_uparse_mSeq_y
ExtractMETA_y(physeq_fungi_uparse_mSeq) -> metadata_fungi_uparse_mSeq_y

ExtractOTU_y(physeq_prok_uparse_mSeq) -> otu_prok_uparse_mSeq_y
ExtractMETA_y(physeq_prok_uparse_mSeq) -> metadata_prok_uparse_mSeq_y

# Fungi --------------------------------------------------------------------------------------------
# test the effect Site 
adonis_fungi <- adonis(t(otu_fungi_uparse_mSeq_y) ~ Year,
                       metadata_fungi_uparse_mSeq_y, 
                       strata=metadata_fungi_uparse_mSeq_y$Distance, 
                       permutations=9999)
adonis_fungi

# test for the effect of Distance 
adonis_fungi1 <- adonis(t(otu_fungi_uparse_mSeq_y) ~ Distance,  
                        metadata_fungi_uparse_mSeq_y, 
                        strata=metadata_fungi_uparse_mSeq_y$Site, 
                        permutations=9999)
adonis_fungi1

# test for interactions - in these case I do no include the strata
# since Distance which is at the beginning of the model is treated as a random
# effect 
adonis_fungi2 <- adonis(t(otu_fungi_uparse_mSeq_y) ~ Distance * Year,  
                        metadata_fungi_uparse_mSeq_y, 
                        permutations=9999)
adonis_fungi2

adonis_fungi3_y <- adonis(t(otu_fungi_uparse_mSeq_y) ~ Year * Distance,
                          metadata_fungi_uparse_mSeq_y, 
                          permutations=9999)
adonis_fungi3_y
round(p.adjust(adonis_fungi3_y$aov.tab$`Pr(>F)`, "bonferroni"), 4)
# [1] 0.0033 1.0000 1.0000     NA     NA

#                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Year           1    0.7042 0.70416 1.76779 0.05688 0.0009 ***
#   Distance       3    1.0710 0.35699 0.89622 0.08650 0.8611    
# Year:Distance  3    1.0455 0.34850 0.87491 0.08445 0.9075    
# Residuals     24    9.5598 0.39833         0.77217           
# Total         31   12.3804                 1.00000  


# PRokaryotes --------------------------------------------------------------------------------------
# test the effect Site 
adonis_prok <- adonis(t(otu_prok_uparse_mSeq_y) ~ Year,
                      metadata_prok_uparse_mSeq_y, 
                      strata=metadata_prok_uparse_mSeq_y$Distance, 
                      permutations=9999)
adonis_prok

# test for the effect of Distance 
adonis_prok1 <- adonis(t(otu_prok_uparse_mSeq_y) ~ Distance,  
                       metadata_prok_uparse_mSeq_y, 
                       strata=metadata_prok_uparse_mSeq_y$Site, 
                       permutations=9999)
adonis_prok1

# test for interactions - in these case I do no include the strata
# since Distance which is at the beginning of the model is treated as a random
# effect 
adonis_prok2 <- adonis(t(otu_prok_uparse_mSeq_y) ~ Distance * Year,  
                       metadata_prok_uparse_mSeq_y, 
                       permutations=9999)
adonis_prok2

adonis_prok3_y <- adonis(t(otu_prok_uparse_mSeq_y) ~ Year * Distance,
                         metadata_prok_uparse_mSeq_y, 
                         permutations=9999)
adonis_prok3_y
round(p.adjust(adonis_prok3_y$aov.tab$`Pr(>F)`, "bonferroni"), 4)
#[1] 0.0012 0.0912 1.0000     NA     NA

#                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Year           1    0.4289 0.42890  3.2038 0.09383 0.0004 ***
#   Distance       3    0.5744 0.19147  1.4303 0.12567 0.0304 *  
#   Year:Distance  3    0.3547 0.11824  0.8832 0.07760 0.7023    
# Residuals     24    3.2130 0.13387         0.70290           
# Total         31    4.5710                 1.00000  


# BETADISPER - analysis of betadispersion ----------------------------------------------------------

# Gio', I think you have the code for this, right?
# It needs to be run on both fungi and bacteria on Site
# and Distance. You could also plot the distance to centroids
# as boxplots, it is nice. 

# Fungi --------------------------------------------------------------------------------------------
sample_data(physeq_fungi_uparse_mSeq)$LocYear <- paste(sample_data(physeq_fungi_uparse_mSeq)$Site,
                                                       sample_data(physeq_fungi_uparse_mSeq)$Year,
                                                       sep = "_")

sample_data(physeq_fungi_uparse_mSeq)$Distance <- factor(sample_data(physeq_fungi_uparse_mSeq)$Distance, 
                                                         levels=c("Inside", "0", "40", "100", "200"))
physeq_fungi_uparse_mSeq %>%
  subset_samples(Type=="soil" & LocYear%in%c("Citta_di_Castello_2014","San_Giovanni_dAsso_2015","Ripabianca_2014")) %>%
  ordinate(method ="PCoA", distance="bray") -> pcoa_fungi_mSeq_site

physeq_fungi_uparse_mSeq %>%
  subset_samples(Type=="soil" & LocYear%in%c("Citta_di_Castello_2014","Citta_di_Castello_2015")) %>%
  ordinate(method ="PCoA", distance="bray") -> pcoa_fungi_mSeq_year


# Prokaryotes --------------------------------------------------------------------------------------
sample_data(physeq_prok_uparse_mSeq)$LocYear <- paste(sample_data(physeq_prok_uparse_mSeq)$Site,
                                                      sample_data(physeq_prok_uparse_mSeq)$Year,
                                                      sep = "_")

sample_data(physeq_prok_uparse_mSeq)$Distance <- factor(sample_data(physeq_prok_uparse_mSeq)$Distance, 
                                                        levels=c("Inside", "0", "40", "100", "200"))
physeq_prok_uparse_mSeq %>%
  subset_samples(Type=="soil" & LocYear%in%c("Citta_di_Castello_2014","San_Giovanni_dAsso_2015","Ripabianca_2014")) %>%
  ordinate(method ="PCoA", distance="bray") -> pcoa_prok_mSeq_site

physeq_prok_uparse_mSeq %>%
  subset_samples(Type=="soil" & LocYear%in%c("Citta_di_Castello_2014","Citta_di_Castello_2015")) %>%
  ordinate(method ="PCoA", distance="bray") -> pcoa_prok_mSeq_year

physeq_prok_uparse_mSeq %>%
  subset_samples(Type=="truffle" & LocYear%in%c("Citta_di_Castello_2014","San_Giovanni_dAsso_2015","Ripabianca_2014")) %>%
  ordinate(method ="PCoA", distance="bray") -> pcoa_prok_mSeq_truf


PlotOrdin <-function(dataframe, ord, var){
  ord <- plot_ordination(dataframe, ord, color=var, shape="Distance") + 
    geom_point(size=1.2, alpha=0.9, aes(shape=Distance)) +
    scale_shape_manual(values = c(15, 16, 17, 18)) + #  labels=c("","", ...)
    theme_classic() +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8)) +
    theme(strip.text.x = element_text(size = 8, face = "bold")) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) +
    grids(linetype = "dashed") +
    theme(legend.position="bottom") +
    guides(color=guide_legend(nrow=2, order = 2), shape=guide_legend(nrow=2, order = 1) )
  return(ord)
}

PlotOrdin(physeq_fungi_uparse_mSeq, pcoa_fungi_mSeq_site, "Site") + 
  scale_colour_manual("Site", values=palette_site, 
                      labels = c("Site 1", "Site 2", "Site 3")) + 
  labs(title = "Fungi") +
  annotate("text", Inf, -Inf, 
           label = expression(paste("Site ", italic(R) ^ 2,"= 7.2%***"), parse=TRUE), 
           size = 3, hjust = 1.2, vjust = -0.5) +
  annotate("text", -Inf, -Inf, 
           label = expression(paste("Depth ", italic(R) ^ 2,"= 7.2%***"), parse=TRUE),
           size = 3, hjust = -0.2, vjust = -0.5)


PlotOrdin(physeq_fungi_uparse_mSeq, pcoa_fungi_mSeq_year, "Year") + 
  scale_colour_manual("Year", values=palette_year, 
                      labels = c("2014", "2015")) + 
  labs(title = "Fungi")

PlotOrdin(physeq_prok_uparse_mSeq, pcoa_prok_mSeq_truf, "Site") + 
  scale_colour_manual("Site", values=paletteCB3_prok, 
                      labels = c("Citta di Castello'", "San Giovanni", "Ripabianca")) + 
  labs(title = "Prokaryotes")


# *** FIGURE 3 - Ordination ------------------------------------------------------------------------
library(grid)
library(gridExtra)

plots1 <- ggarrange(PlotOrdin(physeq_fungi_uparse_mSeq, pcoa_fungi_mSeq_site, "Site") + 
                      scale_colour_manual("Site", values=palette_site, 
                                          labels = c("Site 1", "Site 2", "Site 3")) + 
                      labs(title = "Fungi") +
                      annotate("text", -Inf, Inf, 
                               label = expression(paste("Site ", italic(R) ^ 2,"= 12.2% **"), parse=TRUE),
                               size = 2.7, hjust = -0.05, vjust =  1),
                    PlotOrdin(physeq_prok_uparse_mSeq, pcoa_prok_mSeq_site, "Site") + 
                      scale_colour_manual("Site", values=palette_site, 
                                          labels = c("Site 1", "Site 2", "Site 3")) + 
                      labs(title = "Prokaryotes") +
                      annotate("text", -Inf, Inf, 
                               label = expression(paste("Site ", italic(R) ^ 2,"= 24.8% ***"), parse=TRUE),
                               size = 2.7, hjust = -0.05, vjust = 1) +
                      annotate("text", -Inf, Inf, 
                               label = expression(paste("Distance ", italic(R) ^ 2,"= 8.1% **"), parse=TRUE),
                               size = 2.7, hjust = -0.04, vjust = 2.5),
                    labels = c("A", "B"),
                    widths = c(1,1),
                    align = "hv", 
                    ncol = 2, 
                    nrow = 1,
                    common.legend = TRUE,
                    legend = "bottom")

title1=text_grob("Site",size=12, face=2)
grid.arrange(plots1, top = title1)

plots2 <- ggarrange(PlotOrdin(physeq_fungi_uparse_mSeq, pcoa_fungi_mSeq_year, "Year") + 
                      scale_colour_manual("Site", values=palette_year, 
                                          labels = c("2014", "2015")) + 
                      labs(title = "Fungi") +
                      annotate("text", Inf, -Inf, 
                               label = expression(paste("Year ", italic(R) ^ 2,"= 5.7% **"), parse=TRUE), 
                               size = 2.7, hjust = 1.2, vjust = -0.5) , 
                    PlotOrdin(physeq_prok_uparse_mSeq, pcoa_prok_mSeq_year, "Year") + 
                      scale_colour_manual("Site", values=palette_year, 
                                          labels = c("2014", "2015")) + 
                      labs(title = "Prokaryotes") +
                      annotate("text", Inf, -Inf, 
                               label = expression(paste("Year ", italic(R) ^ 2,"= 9.4% **"), parse=TRUE), 
                               size = 2.7, hjust = 1.2, vjust = -0.5) ,
                    labels = c("C", "D"),
                    widths = c(1,1),
                    align = "hv", 
                    ncol = 2, 
                    nrow = 1,
                    common.legend = TRUE,
                    legend = "bottom")

title2=text_grob("Year",size=12, face=2)
grid.arrange(plots2, top = title2)


Fig_3 = ggarrange(grid.arrange(plots1, top = title1),
                  grid.arrange(plots2, top = title2),
                  widths = c(1,1),
                  align = "hv", 
                  ncol = 1, 
                  nrow = 2)

Fig_3

# ************************************--------------------------------------------------------------
# ANALYZING TRUFFLES MICROBIOME --------------------------------------------------------------------

count(as.data.frame(as.matrix(
  sample_data(
    physeq_fungi_uparse_mSeq %>%
      subset_samples(
        Type == "truffle" & LocYear %in% c(
          "Citta_di_Castello_2014",
          "San_Giovanni_dAsso_2015",
          "Ripabianca_2014"))))), Site)

count(as.data.frame(as.matrix(
  sample_data(
    physeq_fungi_uparse_mSeq %>%
      subset_samples(
        Type == "truffle" & LocYear %in% c(
          "Citta_di_Castello_2014",
          "Citta_di_Castello_2015"))))), Year)


ExtractOTU_t <- function(physeq){
  physeq %>%
    subset_samples(Type=="truffle" & LocYear%in%c("Citta_di_Castello_2014",
                                                  "San_Giovanni_dAsso_2015",
                                                  "Ripabianca_2014")) -> Physeq_filt
  otu <- as.data.frame(otu_table(Physeq_filt))
  return(otu)
}

ExtractMETA_t <- function(physeq){
  physeq %>%
    subset_samples(Type=="truffle" & LocYear%in%c("Citta_di_Castello_2014",
                                                  "San_Giovanni_dAsso_2015",
                                                  "Ripabianca_2014")) -> Physeq_filt
  meta <- as(sample_data(Physeq_filt), "data.frame")
  return(meta)
}


ExtractOTU_t(physeq_fungi_uparse_mSeq) -> otu_fungi_truf
ExtractMETA_t(physeq_fungi_uparse_mSeq) -> metadata_fungi_truf

ExtractOTU_t(physeq_prok_uparse_mSeq) -> otu_prok_truf
ExtractMETA_t(physeq_prok_uparse_mSeq) -> metadata_prok_truf

# test the effect Site 
adonis_fungi_t <- adonis(t(otu_fungi_truf) ~ Site,
                         metadata_fungi_truf, 
                         permutations=9999)
adonis_fungi_t

adonis_prok_t <- adonis(t(otu_prok_truf) ~ Site,
                        metadata_prok_truf, 
                        permutations=9999)
adonis_prok_t


ExtractOTU_t_y <- function(physeq){
  physeq %>%
    subset_samples(Type=="truffle" & LocYear%in%c("Citta_di_Castello_2014",
                                                  "Citta_di_Castello_2015")) -> Physeq_filt
  otu <- as.data.frame(otu_table(Physeq_filt))
  return(otu)
}

ExtractMETA_t_y <- function(physeq){
  physeq %>%
    subset_samples(Type=="truffle" & LocYear%in%c("Citta_di_Castello_2014",
                                                  "Citta_di_Castello_2015")) -> Physeq_filt
  meta <- as(sample_data(Physeq_filt), "data.frame")
  return(meta)
}


ExtractOTU_t_y(physeq_fungi_uparse_mSeq) -> otu_fungi_truf_year
ExtractMETA_t_y(physeq_fungi_uparse_mSeq) -> metadata_fungi_truf_year

ExtractOTU_t_y(physeq_prok_uparse_mSeq) -> otu_prok_truf_year
ExtractMETA_t_y(physeq_prok_uparse_mSeq) -> metadata_prok_truf_year


# test the effect Year 
adonis_fungi_t_year <- adonis(t(otu_fungi_truf_year) ~ Year,
                              metadata_fungi_truf_year, 
                              permutations=9999)
adonis_fungi_t_year

adonis_prok_t_year <- adonis(t(otu_prok_truf_year) ~ Year,
                             metadata_prok_truf_year, 
                             permutations=9999)
adonis_prok_t_year

# Test the homogenity of group variances
permdisp_fungi_truf_year <- 
  betadisper(vegdist(t(otu_fungi_truf_year), method="bray"), metadata_fungi_truf_year$Year) 
anova(permdisp_fungi_truf_year, permutations = 9999)


permdisp_prok_truf_year <- 
  betadisper(vegdist(t(otu_prok_truf_year), method="bray"), metadata_prok_truf_year$Year) 
anova(permdisp_prok_truf_year, permutations = 9999)


# Plotting
PlotOrdinTruf <-function(dataframe, ord, var){
  ord <- plot_ordination(dataframe, ord, color=var) + 
    geom_point(size=1.2, alpha=0.9) +
    scale_shape_manual(values = c(15, 16, 17, 18)) + #  labels=c("","", ...)
    theme_classic() +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8)) +
    theme(strip.text.x = element_text(size = 8, face = "bold")) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) +
    grids(linetype = "dashed") +
    theme(legend.position="bottom") +
    guides(color=guide_legend(nrow=2, order = 2), shape=guide_legend(nrow=2, order = 1) )
  return(ord)
}



physeq_fungi_uparse_mSeq %>%
  subset_samples(Type=="truffle" & LocYear%in%c("Citta_di_Castello_2014",
                                                "San_Giovanni_dAsso_2015",
                                                "Ripabianca_2014")) %>%
  ordinate(method ="PCoA", distance="bray") -> pcoa_fungi_truf_site


physeq_prok_uparse_mSeq %>%
  subset_samples(Type=="truffle" & LocYear%in%c("Citta_di_Castello_2014",
                                                "San_Giovanni_dAsso_2015",
                                                "Ripabianca_2014")) %>%
  ordinate(method ="PCoA", distance="bray") -> pcoa_prok_truf_site



physeq_fungi_uparse_mSeq %>%
  subset_samples(Type=="truffle" & LocYear%in%c("Citta_di_Castello_2014",
                                                "Citta_di_Castello_2015")) %>%
  ordinate(method ="PCoA", distance="bray") -> pcoa_fungi_truf_year


physeq_prok_uparse_mSeq %>%
  subset_samples(Type=="truffle" & LocYear%in%c("Citta_di_Castello_2014",
                                                "Citta_di_Castello_2015")) %>%
  ordinate(method ="PCoA", distance="bray") -> pcoa_prok_truf_year



PlotOrdinTruf(physeq_fungi_uparse_mSeq, pcoa_fungi_truf_year, "Year") +
  scale_colour_manual("Year", values=palette_year, labels = c("2014", "2015")) + 
  labs(title = "Fungi")

PlotOrdinTruf(physeq_prok_uparse_mSeq, pcoa_prok_truf_year, "Year") +
  scale_colour_manual("Year", values=palette_year, labels = c("2014", "2015")) + 
  labs(title = "Fungi")



# FIGURE S4 - truffle geography and year -----------------------------------------------------------

plots_truf_site <- ggarrange(
  PlotOrdinTruf(physeq_fungi_uparse_mSeq, pcoa_fungi_truf_site, "Site") +
    scale_colour_manual("Site", values=palette_site, 
                        labels = c("Site 1", "Site 2", "Site 3")) +
    labs(title = "Fungi") +
    annotate("text", -Inf, Inf, 
             label = expression(paste("Site ", italic(R) ^ 2,"= 6.8% ns"), parse=TRUE),
             size = 2.3, hjust = -0.05, vjust =  1),
  PlotOrdinTruf(physeq_prok_uparse_mSeq, pcoa_prok_truf_site, "Site") +
    scale_colour_manual("Site", values=palette_site, 
                        labels = c("Site 1", "Site 2", "Site 3")) +
    labs(title = "Prokaryotes") +
    annotate("text", -Inf, Inf, 
             label = expression(paste("Site ", italic(R) ^ 2,"= 19.9% ns"), parse=TRUE),
             size = 2.3, hjust = -0.05, vjust = 1),
  labels = c("A", "B"),
  widths = c(1,1),
  align = "hv", 
  ncol = 2, 
  nrow = 1,
  common.legend = TRUE,
  legend = "bottom")

plots_truf_site


plots_truf_year <- 
  ggarrange(
    PlotOrdinTruf(physeq_fungi_uparse_mSeq, pcoa_fungi_truf_year, "Year") +
      scale_colour_manual("Year", values=palette_year, labels = c("2014", "2015")) +
      labs(title = "Fungi") +
      annotate("text", -Inf, Inf, 
               label = expression(paste("Year ", italic(R) ^ 2,"= 30.0% *"), parse=TRUE),
               size = 2.3, hjust = -0.05, vjust =  1) +
      annotate("text", -Inf, Inf, 
               label = expression(paste("Dispersion ", italic(p),"= 0.038"), parse=TRUE),
               size = 2.3, hjust = -0.05, vjust =  2.5),
    PlotOrdinTruf(physeq_prok_uparse_mSeq, pcoa_prok_truf_year, "Year") +
      scale_colour_manual("Year", values=palette_year, labels = c("2014", "2015")) + 
      labs(title = "Prokaryotes") +
      annotate("text", Inf, -Inf, 
               label = expression(paste("Year ", italic(R) ^ 2,"= 27.9% *"), parse=TRUE),
               size = 2.3, hjust = 1, vjust = -1.5) +
      annotate("text", Inf, -Inf, 
               label = expression(paste("Dispersion ", italic(p),"= 0.036"), parse=TRUE),
               size = 2.3, hjust = 1, vjust = -0.5),
    labels = c("C", "D"),
    widths = c(1,1),
    align = "hv", 
    ncol = 2, 
    nrow = 1,
    common.legend = TRUE,
    legend = "bottom")


plots_truf_year


Fig_3_suppl = ggarrange(
  grid.arrange(plots_truf_site, top = title1),
  grid.arrange(plots_truf_year, top = title2),
  widths = c(1,1),
  align = "hv", 
  ncol = 1, 
  nrow = 2)

Fig_3_suppl


# ************************************--------------------------------------------------------------
# >>> HEATMAP TREES --------------------------------------------------------------------------------
library(metacoder)
library(tibble)

# Extracting Dataset at for Tuffle and all soils ---------------------------------------------------
physeq_prok_uparse_filt_fam <- physeq_prok_uparse_filt
physeq_prok_uparse_filt_fam@tax_table <- physeq_prok_uparse_filt_fam@tax_table[, 2:7]

physeq_prok_uparse_filt_fam %>% 
  tax_glom(taxrank="Family") -> physeq_prok_fam

otu_table(physeq_prok_fam) <- otu_table(physeq_prok_fam)[which(
  rowSums(otu_table(physeq_prok_fam)) > 0),]

physeq_prok_fam
head(tax_table(physeq_prok_fam))
head(sample_data(physeq_prok_fam))

# HEAT TREE Bacteria by Family ---------------------------------------------------------------------
bacteria_family_gio <- physeq_prok_fam

# otu table
otu_bacteria_fam <- as.data.frame(otu_table(bacteria_family_gio))
head(otu_bacteria_fam)

metadata_bacteria_fam <- as.data.frame(as.matrix(sample_data(bacteria_family_gio)))
metadata_bacteria_fam$sample_id <- rownames(metadata_bacteria_fam)
metadata_bacteria_fam

# remove NA
tax_table(bacteria_family_gio)[is.na(tax_table(bacteria_family_gio))] <- ""

taxa_bacteria_fam <- as.matrix(tax_table(bacteria_family_gio))
taxa_bacteria_fam[,c(2:7)] -> taxa_bacteria_fam
as.data.frame(taxa_bacteria_fam) -> taxa_bacteria_fam
head(taxa_bacteria_fam)
taxa_bacteria_fam[, 6]

# Reformatting taxonomy

cols <- c("Domain", "Kingdom","Phylum","Class","Order","Family")
taxa_bacteria_fam$Taxonomy <- do.call(paste, c(taxa_bacteria_fam[cols], sep=";"))
head(taxa_bacteria_fam)

identical(rownames(taxa_bacteria_fam), rownames(otu_bacteria_fam))
otu_bacteria_fam$Taxonomy <- taxa_bacteria_fam$Taxonomy
head(otu_bacteria_fam)
dim(otu_bacteria_fam)

# transform a data.frame into a tibble
tibbl_bacteria_fam <- tibble::as_tibble(otu_bacteria_fam)
tibbl_bacteria_fam

tibbl_metadata_bacteria_fam <- tibble::as_tibble(metadata_bacteria_fam)
tibbl_metadata_bacteria_fam

identical(colnames(tibbl_bacteria_fam[-61]), tibbl_metadata_bacteria_fam$sample_id)
#sample_order <- match(tibbl_metadata_fungi$sample_id, colnames(tibbl_fungi))
#tibbl_fungi <- tibbl_fungi[,sample_order]

# create metacoder object
heat_tree_16s_fam <- parse_tax_data(tibbl_bacteria_fam,
                                    class_cols = "Taxonomy",
                                    class_sep = ";")
print(heat_tree_16s_fam)
heat_tree_16s_fam$data$tax_data

# Removing low-abundance counts
# We need to thing about this very well. Maybe better usinf the
# non-normalized data!
heat_tree_16s_fam$data$tax_data <- zero_low_counts(heat_tree_16s_fam, "tax_data", min_count = 50)
no_reads_fam <- rowSums(heat_tree_16s_fam$data$tax_data[, tibbl_metadata_bacteria_fam$sample_id]) == 0
sum(no_reads_fam)
# removing OTUs with no reads left 
heat_tree_16s_fam <- 
  filter_obs(heat_tree_16s_fam, data  = "tax_data", ! no_reads_fam, drop_taxa = TRUE)
print(heat_tree_16s_fam)
# calculate taxa proprostion
heat_tree_16s_fam$data$tax_data <- 
  calc_obs_props(heat_tree_16s_fam, "tax_data")
# abundance per-taxon
heat_tree_16s_fam$data$tax_abund <- 
  calc_taxon_abund(heat_tree_16s_fam, "tax_data", cols = tibbl_metadata_bacteria_fam$sample_id)
# number of samples have reads for each taxon
heat_tree_16s_fam$data$tax_occ <- calc_n_samples(heat_tree_16s_fam, "tax_abund", 
                                                 cols = tibbl_metadata_bacteria_fam$sample_id,
                                                 groups = tibbl_metadata_bacteria_fam$Type)


print(heat_tree_16s_fam$data$tax_occ)
# I think the column order match but we should verify! 
# Maybe adding the cols parameter would be better.
# adding taxon abundance
heat_tree_16s_fam$data$type_abund <- calc_group_mean(heat_tree_16s_fam, "tax_abund",
                                                     cols = tibbl_metadata_bacteria_fam$sample_id,
                                                     groups = tibbl_metadata_bacteria_fam$Type)
print(heat_tree_16s_fam$data$type_abund)

# Comparing groups
tibbl_metadata_bacteria_fam$Type <- factor(tibbl_metadata_bacteria_fam$Type, levels = c("soil", "truffle"))
levels(tibbl_metadata_bacteria_fam$Type)

heat_tree_16s_fam$data$diff_table <- compare_groups(heat_tree_16s_fam, data = "tax_abund",
                                                    cols = tibbl_metadata_bacteria_fam$sample_id, # What columns of sample data to use
                                                    groups = tibbl_metadata_bacteria_fam$Type) # What category each sample is assigned to

# compare results
print(heat_tree_16s_fam$data$diff_table)
print(heat_tree_16s$data$diff_table)

# Invert the significance to plot to the same direction and color
heat_tree_16s$data$diff_table$log2_median_ratio <-
  heat_tree_16s$data$diff_table$log2_median_ratio*-1

# Correcting P-values
heat_tree_16s_fam$data$diff_table$wilcox_p_value <- p.adjust(heat_tree_16s_fam$data$diff_table$wilcox_p_value, method = "fdr")

# *** FIGURE 4 A -----------------------------------------------------------------------------------
set.seed(12)

Fig_4A <-
  heat_tree(heat_tree_16s_fam,
            node_size = heat_tree_16s_fam$n_obs(), # n_obs is a function that calculates, in this case, the number of OTUs per taxon
            node_label = heat_tree_16s_fam$taxon_names(),
            node_color = heat_tree_16s_fam$data$diff_table$log2_median_ratio, # A column from `obj$data$diff_table`
            node_color_range = c("gold", "grey","brown"), 
            node_color_trans = "linear", # The default is scaled by circle area
            node_label_size_range = c(0.0115, 0.026), #range of nodes label size
            node_label_size_trans = "linear",
            edge_size_range = c(0.003, 0.003),
            node_size_range = c(0.01, 0.05),
            tree_label_size = as.numeric(0.5),
            #node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
            #edge_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
            node_size_axis_label = "Number of OTUs",
            node_color_axis_label = "Log2 ratio median proportions",
            layout = "davidson-harel", # The primary layout algorithm
            initial_layout = "reingold-tilford",
            output_file = "Fig_4_A_Fam_red_blue.pdf") # The layout algorithm that initializes node Sites

Fig_4A


# Extracting Dataset at for Tuffle iand 0 cm -------------------------------------------------------
physeq_prok_uparse_filt_gen <- physeq_prok_uparse_filt
physeq_prok_uparse_filt_gen@tax_table <- physeq_prok_uparse_filt_gen@tax_table[, 2:7]

physeq_prok_uparse_filt_gen %>% 
  subset_samples(Distance%in%c("Inside","0")) %>%
tax_glom(taxrank="Genus") -> physeq_prok_gen

otu_table(physeq_prok_gen) <- otu_table(physeq_prok_gen)[which(
  rowSums(otu_table(physeq_prok_gen)) > 0),] 

physeq_prok_gen
head(tax_table(physeq_prok_gen))
head(sample_data(physeq_prok_gen))

# HEAT TREE Bacteria by Genus ----------------------------------------------------------------------
bacteria_genera_gio <- physeq_prok_gen

# otu table
otu_bacteria_gen <- as.data.frame(otu_table(bacteria_genera_gio))
head(otu_bacteria_gen)
dim(otu_bacteria_gen)

metadata_bacteria_gen <- as.data.frame(as.matrix(sample_data(bacteria_genera_gio)))
metadata_bacteria_gen$sample_id <- rownames(metadata_bacteria_gen)
metadata_bacteria_gen
dim(metadata_bacteria_gen)

# remove NA
tax_table(bacteria_genera_gio)[is.na(tax_table(bacteria_genera_gio))] <- ""

taxa_bacteria_gen <- as.matrix(tax_table(bacteria_genera_gio))
taxa_bacteria_gen[,c(2:8)] -> taxa_bacteria_gen
as.data.frame(taxa_bacteria_gen) -> taxa_bacteria_gen
head(taxa_bacteria_gen)
str(taxa_bacteria_gen)

# Reformatting and correcting taxonomy

taxa_bacteria_gen[, colnames(taxa_bacteria_gen)][,7]
taxa_bacteria_gen$Genus <- as.character(taxa_bacteria_gen$Genus)

taxa_bacteria_gen[taxa_bacteria_gen == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Rhizobium"
taxa_bacteria_gen[taxa_bacteria_gen == "Allorhizobium-neorhizobium-pararhizobium-rhizobium"] <- "Rhizobium"
taxa_bacteria_gen[taxa_bacteria_gen == "Burkholderia-caballeronia-paraburkholderia"] <- "Burkholderia"
taxa_bacteria_gen[taxa_bacteria_gen == "Burkholderia"] 

cols <- c("Domain", "Kingdom","Phylum","Class","Order","Family", "Genus")
taxa_bacteria_gen$Taxonomy <- do.call(paste, c(taxa_bacteria_gen[cols], sep=";"))
head(taxa_bacteria_gen)
taxa_bacteria_gen[, 1:7]
taxa_bacteria_gen[, 7:8]


identical(rownames(taxa_bacteria_gen), rownames(otu_bacteria_gen))
otu_bacteria_gen$Taxonomy <- taxa_bacteria_gen$Taxonomy
head(otu_bacteria_gen)

# transform a data.frame into a tibble
tibbl_bacteria <- tibble::as_tibble(otu_bacteria_gen)
tibbl_bacteria

tibbl_metadata_bacteria <- tibble::as_tibble(metadata_bacteria_gen)
tibbl_metadata_bacteria

colnames(tibbl_bacteria)[25] # is the taxonomy columns

identical(colnames(tibbl_bacteria)[-25], tibbl_metadata_bacteria$sample_id)
sample_order <- match(tibbl_metadata_fungi$sample_id, colnames(tibbl_fungi))
tibbl_fungi <- tibbl_fungi[,sample_order]

# create metacoder object
heat_tree_16s <- parse_tax_data(tibbl_bacteria,
                                class_cols = "Taxonomy",
                                class_sep = ";")
print(heat_tree_16s)
heat_tree_16s$data$tax_data

# Removing low-abundance counts
# We need to thing about this very well. Maybe better usinf the
# non-normalized data!
heat_tree_16s$data$tax_data <- zero_low_counts(heat_tree_16s, "tax_data", min_count = 50)
no_reads <- rowSums(heat_tree_16s$data$tax_data[, tibbl_metadata_bacteria$sample_id]) == 0
sum(no_reads)
# removing OTUs with no reads left 
heat_tree_16s <- filter_obs(heat_tree_16s, data  = "tax_data", ! no_reads, drop_taxa = TRUE)
print(heat_tree_16s)
# calculate taxa proprostion
heat_tree_16s$data$tax_data <- calc_obs_props(heat_tree_16s, "tax_data")
# abundance per-taxon
heat_tree_16s$data$tax_abund <- calc_taxon_abund(heat_tree_16s, "tax_data", cols = tibbl_metadata_bacteria$sample_id)
# number of samples have reads for each taxon
heat_tree_16s$data$tax_occ <- calc_n_samples(heat_tree_16s, "tax_abund", 
                                             cols = tibbl_metadata_bacteria$sample_id,
                                             groups = tibbl_metadata_bacteria$Type)

print(heat_tree_16s$data$tax_occ)
# I think the column order match but we should verify! 
# Maybe adding the cols parameter would be better.
# adding taxon abundance
heat_tree_16s$data$type_abund <- calc_group_mean(heat_tree_16s, "tax_abund",
                                                 cols = tibbl_metadata_bacteria$sample_id,
                                                 groups = tibbl_metadata_bacteria$Type)
print(heat_tree_16s$data$type_abund)

# Comparing groups
tibbl_metadata_bacteria$Type <- factor(tibbl_metadata_bacteria$Type, levels = c("soil", "truffle"))
levels(tibbl_metadata_bacteria$Type)

heat_tree_16s$data$diff_table <- compare_groups(heat_tree_16s, data = "tax_abund",
                                                cols = tibbl_metadata_bacteria$sample_id, # What columns of sample data to use
                                                groups = tibbl_metadata_bacteria$Type) # What category each sample is assigned to

print(heat_tree_16s$data$diff_table)

# Adjusting p-values
heat_tree_16s$data$diff_table$wilcox_p_value <- p.adjust(heat_tree_16s$data$diff_table$wilcox_p_value, method = "fdr")

# *** FIGURE 4 B -----------------------------------------------------------------------------------
set.seed(13)

Fig_4B <-
  heat_tree(heat_tree_16s,
            node_size = heat_tree_16s$n_obs(), # n_obs is a function that calculates, in this case, the number of OTUs per taxon
            node_label = heat_tree_16s$taxon_names(),
            node_color = heat_tree_16s$data$diff_table$log2_median_ratio, # A column from `obj$data$diff_table`
            node_color_range = c("gold", "grey","brown"), # The built-in palette for diverging data
            node_color_trans = "linear", # The default is scaled by circle area
            node_label_size_range = c(0.0115, 0.026), #range of nodes label size c(0.0125, 0.03) c(0.0095, 0.025)
            node_label_size_trans = "linear",
            edge_size_range = c(0.003, 0.003),
            node_size_range = c(0.01, 0.05),
            tree_label_size = as.numeric(0.5),
            #node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
            #edge_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
            node_size_axis_label = "Number of OTUs",
            node_color_axis_label = "Log2 ratio median proportions",
            layout = "davidson-harel", # The primary layout algorithm
            initial_layout = "reingold-tilford",
            output_file = "Figure_heattree_Sites.pdf") # The layout algorithm that initializes node Sites

Fig_4B

# *** Figure 4 - complete ---------------------------------------------------------------------------
library(patchwork)

Fig_4A + ggtitle("A") +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  Fig_4B + ggtitle("B") + 
  theme(plot.title = element_text(size = 20, face = "bold"))


# ************************************--------------------------------------------------------------
# Supplementary HEAT TREE --------------------------------------------------------------------------
#Soil+truffle - Years Dataset -Prokaryotes "Family"
Phy_bacteria_CdC

bacteria_family_year %>% 
  subset_samples(Distance%in%c("Inside","0")) %>%
  tax_glom(taxrank="Family") -> bacteria_family_year

otu_table(bacteria_family_year) <- otu_table(bacteria_family_year)[which(
  rowSums(otu_table(bacteria_family_year)) > 0),]

bacteria_family_year
head(tax_table(bacteria_family_year))
head(sample_data(bacteria_family_year))

# otu table
otu_bacteria_fam_y <- as.data.frame(otu_table(bacteria_family_year))
head(otu_bacteria_fam_y)

metadata_bacteria_fam_y <- as.data.frame(as.matrix(sample_data(bacteria_family_year)))
metadata_bacteria_fam_y$sample_id <- rownames(metadata_bacteria_fam_y)
metadata_bacteria_fam_y

# remove NA
tax_table(bacteria_family_year)[is.na(tax_table(bacteria_family_year))] <- ""

taxa_bacteria_fam_y <- as.matrix(tax_table(bacteria_family_year))
taxa_bacteria_fam_y[,c(2:7)] -> taxa_bacteria_fam_y
as.data.frame(taxa_bacteria_fam_y) -> taxa_bacteria_fam_y
head(taxa_bacteria_fam_y)

# Reformatting taxonomy
cols_y <- c("Domain", "Kingdom","Phylum","Class","Order","Family")
taxa_bacteria_fam_y$Taxonomy <- do.call(paste, c(taxa_bacteria_fam_y[cols_y], sep=";"))
head(taxa_bacteria_fam_y)

identical(rownames(taxa_bacteria_fam_y), rownames(otu_bacteria_fam_y))
otu_bacteria_fam_y$Taxonomy <- taxa_bacteria_fam_y$Taxonomy
head(otu_bacteria_fam_y)

# transform a data.frame into a tibble
tibbl_bacteria_y <- tibble::as_tibble(otu_bacteria_fam_y)
tibbl_bacteria_y

tibbl_metadata_bacteria_y <- tibble::as_tibble(metadata_bacteria_fam_y)
tibbl_metadata_bacteria_y

identical(colnames(tibbl_bacteria_y[-41]), metadata_bacteria_fam_y$sample_id)
#sample_order <- match(tibbl_metadata_fungi$sample_id, colnames(tibbl_fungi))
#tibbl_fungi <- tibbl_fungi[,sample_order]

heat_tree_16s_y <- parse_tax_data(tibbl_bacteria_y,
                                  class_cols = "Taxonomy",
                                  class_sep = ";")
print(heat_tree_16s_y)
heat_tree_16s_y$data$tax_data

heat_tree_16s_y$data$tax_data <- zero_low_counts(heat_tree_16s_y, "tax_data", min_count = 50)
no_reads_y <- rowSums(heat_tree_16s_y$data$tax_data[, tibbl_metadata_bacteria_y$sample_id]) == 0
sum(no_reads_y)

heat_tree_16s_y <- filter_obs(heat_tree_16s_y, data  = "tax_data", ! no_reads_y, drop_taxa = TRUE)
print(heat_tree_16s_y)

heat_tree_16s_y$data$tax_data <- calc_obs_props(heat_tree_16s_y, "tax_data")

heat_tree_16s_y$data$tax_abund <- calc_taxon_abund(heat_tree_16s_y, "tax_data", cols = tibbl_metadata_bacteria_y$sample_id)

heat_tree_16s_y$data$tax_occ <- calc_n_samples(heat_tree_16s_y, "tax_abund", 
                                               cols = tibbl_metadata_bacteria_y$sample_id,
                                               groups = tibbl_metadata_bacteria_y$Year)
print(heat_tree_16s_y$data$tax_occ)

heat_tree_16s_y$data$type_abund <- calc_group_mean(heat_tree_16s_y, "tax_abund",
                                                   cols = tibbl_metadata_bacteria_y$sample_id,
                                                   groups = tibbl_metadata_bacteria_y$Year)
print(heat_tree_16s_y$data$type_abund)


# Comparing groups
heat_tree_16s_y$data$diff_table <- compare_groups(heat_tree_16s_y, data = "tax_abund",
                                                  cols = tibbl_metadata_bacteria_y$sample_id, # What columns of sample data to use
                                                  groups = tibbl_metadata_bacteria_y$Year) # What category each sample is assigned to


print(heat_tree_16s_y$data$diff_table)

#Adjusting p-values 
heat_tree_16s_y$data$diff_table$wilcox_p_value <- p.adjust(heat_tree_16s_y$data$diff_table$wilcox_p_value, method = "fdr")

# *** FIGURE S5 - Heat tree by year -----------------------------------------------------------------
set.seed(999)

Fig_S5 <-
  heat_tree(heat_tree_16s_y,
            node_size = heat_tree_16s_y$n_obs(), # n_obs is a function that calculates, in this case, the number of OTUs per taxon
            node_label = heat_tree_16s_y$taxon_names(),
            node_color = heat_tree_16s_y$data$diff_table$log2_median_ratio, # A column from `obj$data$diff_table`
            node_color_range = c("gold", "grey","brown"), # The built-in palette for diverging data
            node_color_trans = "linear", # The default is scaled by circle area
            node_label_size_range = c(0.0115, 0.026), #range of nodes label size
            node_label_size_trans = "linear",
            edge_size_range = c(0.003, 0.003),
            node_size_range = c(0.01, 0.04),
            tree_label_size = as.numeric(0.5),
            #node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
            #edge_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
            node_size_axis_label = "Number of OTUs",
            node_color_axis_label = "Log2 ratio median proportions",
            layout = "davidson-harel", # The primary layout algorithm
            initial_layout = "reingold-tilford",
            #title = "Heat Tree Comparing years 2014 and 2015 (prokaryotes)",
            #title_size = 0.017,
            output_file = "Fig_S5_2014Vs2015.pdf") # The layout algorithm that initializes node Sites

Fig_S5

GetIndicators <-function(dataframe, var){
  require(phyloseq); require(indicspecies); require(dplyr)
  otu <- as.data.frame(otu_table(dataframe))
  metadata = as(sample_data(dataframe), "data.frame")
  taxa <- as.data.frame(as.matrix(tax_table(dataframe)))
  multipatt <- multipatt(t(otu), metadata[,var], func = "r.g",
                         control=how(nperm=9999), duleg=TRUE)
  multipatt -> multipatt_fdr
  multipatt_fdr$sign$p.value <- p.adjust(multipatt_fdr$sign$p.value, "fdr")
  multipatt_fdr$sign[which(
    multipatt_fdr$sign$p.value <= 0.05), ] -> indicator_taxa
  taxa$OTU <- rownames(taxa)
  data.frame(OTU = as.factor(row.names(indicator_taxa)), indicator_taxa) %>%
    dplyr::left_join(taxa, by="OTU") -> indicator_taxa
  rownames(indicator_taxa) <- indicator_taxa$OTU
  indicator_taxa <- arrange(indicator_taxa, desc(stat))
  return(indicator_taxa)
}

# indicator value >0.5 and p-value <0.05 after fdr correction
GetIndicators(Phy_bacteria_CdC_mSeq, "Year") -> ind_16s_Year

ind_ITS_Site
ind_16s_Site
ind_16s_Distance
ind_ITS_Distance
ind_ITS_Year
ind_16s_Year
write.csv(ind_16s_Year, "Indicspecies_tax_16s_Year.csv")

# ************************************---------------------------------------------------------------
# VENN DIAGRAMS -------------------------------------------------------------------------------------
library(limma)

#modified form limma package
VennGraph <- 
  function (object, include = "both", names = NULL, mar = rep(1,
                                                              4), cex = c(1.5, 1, 0.7), lwd = 1, circle.col = NULL, counts.col = NULL,
            show.include = NULL, ...)
  {
    include <- as.character(include)
    LenInc <- min(length(include), 2)
    if (is(object, "VennCounts")) {
      include <- include[1]
      LenInc <- 1
    }
    else {
      if (LenInc > 1)
        z2 <- vennCounts(object, include = include[2])[,
                                                       "Counts"]
      object <- vennCounts(object, include = include[1])
    }
    z <- object[, "Counts"]
    nsets <- ncol(object) - 1
    if (nsets > 5)
      stop("Can't plot Venn diagram for more than 5 sets")
    VennZone <- object[, 1:nsets, drop = FALSE]
    VennZone <- apply(VennZone, 1, function(x) paste(x, sep = "",
                                                     collapse = ""))
    names(z) <- VennZone
    if (length(include) == 2)
      names(z2) <- VennZone
    if (is.null(names))
      names <- colnames(object)[1:nsets]
    FILL.COL <- TRUE
    if (is.null(circle.col)) {
      circle.col <- par("col")
      FILL.COL <- FALSE
    }
    if (length(circle.col) < nsets)
      circle.col <- rep(circle.col, length.out = nsets)
    if (is.null(counts.col))
      counts.col <- par("col")
    if (length(counts.col) < LenInc)
      counts.col <- rep(counts.col, length.out = LenInc)
    if (is.null(show.include))
      show.include <- as.logical(LenInc - 1)
    old.par <- par()$mar
    on.exit(par(mar = old.par))
    par(mar = mar)
    if (nsets <= 3) {
      plot(x = 0, y = 0, type = "n", xlim = c(-4, 4), ylim = c(-4,
                                                               4), xlab = "", ylab = "", axes = FALSE, ...)
      theta <- 2 * pi * (0:360)/360
      xcentres <- switch(nsets, 0, c(-1, 1), c(-1, 1, 0))
      ycentres <- switch(nsets, 0, c(0, 0), c(1, 1, -2)/sqrt(3))
      r <- 1.5
      xtext <- switch(nsets, -1.2, c(-1.2, 1.2), c(-1.2, 1.2,
                                                   0))
      ytext <- switch(nsets, 1.8, c(1.8, 1.8), c(2.4, 2.4,
                                                 -3))
      for (circle in 1:nsets) {
        if (!FILL.COL)
          lines(xcentres[circle] + r * cos(theta), ycentres[circle] +
                  r * sin(theta), lwd = lwd, col = circle.col[circle])
        if (FILL.COL) {
          RGB <- col2rgb(circle.col[circle])/255
          ALPHA <- 0.06
          RGB.ALP <- rgb(RGB[1, 1], RGB[2, 1], RGB[3, 1],
                         alpha = ALPHA)
          polygon(xcentres[circle] + r * cos(theta), ycentres[circle] +
                    r * sin(theta), border = circle.col[circle],
                  lwd = lwd, col = RGB.ALP)
        }
        text(xtext[circle], ytext[circle], names[circle],
             cex = cex)
      }
      #switch(nsets, rect(-3, -2.5, 3, 2.5), rect(-3, -2.5,
      #    3, 2.5), rect(-3, -3.5, 3, 3.3))
      showCounts <- switch(nsets, function(counts, cex, adj,
                                           col, leg) {
        #     text(2.3, -2.1, sprintf("Not in any =  %i", counts[1]), cex = cex, col = col,
        #          adj = adj)
        text(0, 0, counts[2], cex = cex, col = col, adj = adj)
        if (show.include) text(-2.3, -2.1, leg, cex = cex,
                               col = col, adj = adj)
      }, function(counts, cex, adj, col, leg) {
        #    text(2.3, -2.1, sprintf("Not in any = %i", counts[1]), cex = cex, col = col,
        #          adj = adj)
        text(1.5, 0.1, counts[2], cex = cex, col = col, adj = adj)
        text(-1.5, 0.1, counts[3], cex = cex, col = col,
             adj = adj)
        text(0, 0.1, counts[4], cex = cex, col = col, adj = adj)
        if (show.include) text(-2.3, -2.1, leg, cex = cex,
                               col = col, adj = adj)
      }, function(counts, cex, adj, col, leg) {
        #     text(2.5, -3, sprintf("Not in any = %i", counts[1]), cex = cex, col = col, adj = adj)
        text(0, -1.7, counts[2], cex = cex, col = col, adj = adj)
        text(1.5, 1, counts[3], cex = cex, col = col, adj = adj)
        text(0.75, -0.35, counts[4], cex = cex, col = col,
             adj = adj)
        text(-1.5, 1, counts[5], cex = cex, col = col, adj = adj)
        text(-0.75, -0.35, counts[6], cex = cex, col = col,
             adj = adj)
        text(0, 0.9, counts[7], cex = cex, col = col, adj = adj)
        text(0, 0, counts[8], cex = cex, col = col, adj = adj)
        if (show.include) text(-2.5, -3, leg, cex = cex,
                               col = col, adj = adj)
      })
      if (LenInc == 1)
        adj <- c(0.5, 0.5)
      else adj <- c(0.5, 0)
      print(z)
      showCounts(counts = z, cex = cex[1], adj = adj, col = counts.col[1],
                 leg = include[1])
      if (LenInc == 2)
        showCounts(counts = z2, cex = cex[1], adj = c(0.5,
                                                      1), col = counts.col[2], leg = include[2])
      return(invisible())
    }
    plot(c(-20, 420), c(-20, 420), type = "n", axes = FALSE,
         ylab = "", xlab = "", ...)
    relocate_elp <- function(e, alpha, x, y) {
      phi <- (alpha/180) * pi
      xr <- e[, 1] * cos(phi) + e[, 2] * sin(phi)
      yr <- -e[, 1] * sin(phi) + e[, 2] * cos(phi)
      xr <- x + xr
      yr <- y + yr
      cbind(xr, yr)
    }
    if (4 == nsets) {
      #rect(-20, -20, 420, 400)
      elps <- cbind(162 * cos(seq(0, 2 * pi, len = 1000)),
                    108 * sin(seq(0, 2 * pi, len = 1000)))
      if (!FILL.COL) {
        polygon(relocate_elp(elps, 45, 130, 170), border = circle.col[1],
                lwd = lwd)
        polygon(relocate_elp(elps, 45, 200, 200), border = circle.col[2],
                lwd = lwd)
        polygon(relocate_elp(elps, 135, 200, 200), border = circle.col[3],
                lwd = lwd)
        polygon(relocate_elp(elps, 135, 270, 170), border = circle.col[4],
                lwd = lwd)
      }
      if (FILL.COL) {
        RGB <- col2rgb(circle.col)/255
        ALPHA <- 0.06
        RGB.ALP1 <- rgb(RGB[1, 1], RGB[2, 1], RGB[3, 1],
                        alpha = ALPHA)
        RGB.ALP2 <- rgb(RGB[1, 2], RGB[2, 2], RGB[3, 2],
                        alpha = ALPHA)
        RGB.ALP3 <- rgb(RGB[1, 3], RGB[2, 3], RGB[3, 3],
                        alpha = ALPHA)
        RGB.ALP4 <- rgb(RGB[1, 4], RGB[2, 4], RGB[3, 4],
                        alpha = ALPHA)
        polygon(relocate_elp(elps, 45, 130, 170), border = circle.col[1],
                lwd = lwd, col = RGB.ALP1)
        polygon(relocate_elp(elps, 45, 200, 200), border = circle.col[2],
                lwd = lwd, col = RGB.ALP2)
        polygon(relocate_elp(elps, 135, 200, 200), border = circle.col[3],
                lwd = lwd, col = RGB.ALP3)
        polygon(relocate_elp(elps, 135, 270, 170), border = circle.col[4],
                lwd = lwd, col = RGB.ALP4)
      }
      text(35, 315, names[1], cex = cex[1])
      text(138, 350, names[2], cex = cex[1])
      text(262, 347, names[3], cex = cex[1])
      text(365, 315, names[4], cex = cex[1])
      text(35, 250, z["1000"], cex = cex[2], col = counts.col[1])
      text(140, 315, z["0100"], cex = cex[2], col = counts.col[1])
      text(260, 315, z["0010"], cex = cex[2], col = counts.col[1])
      text(365, 250, z["0001"], cex = cex[2], col = counts.col[1])
      text(90, 282, z["1100"], cex = cex[3], col = counts.col[1])
      text(95, 110, z["1010"], cex = cex[2], col = counts.col[1])
      text(200, 52, z["1001"], cex = cex[3], col = counts.col[1])
      text(200, 292, z["0110"], cex = cex[2], col = counts.col[1])
      text(300, 110, z["0101"], cex = cex[2], col = counts.col[1])
      text(310, 282, z["0011"], cex = cex[3], col = counts.col[1])
      text(130, 230, z["1110"], cex = cex[2], col = counts.col[1])
      text(245, 81, z["1101"], cex = cex[3], col = counts.col[1])
      text(155, 81, z["1011"], cex = cex[3], col = counts.col[1])
      text(270, 230, z["0111"], cex = cex[2], col = counts.col[1])
      text(200, 152, z["1111"], cex = cex[2], col = counts.col[1])
      # text(400, 15, sprintf("Not in any = %i", z["0000"]), cex = cex[1], col = counts.col[1])
      if (length(include) == 2) {
        text(35, 238, z2["1000"], cex = cex[2], col = counts.col[2])
        text(140, 304, z2["0100"], cex = cex[2], col = counts.col[2])
        text(260, 304, z2["0010"], cex = cex[2], col = counts.col[2])
        text(365, 238, z2["0001"], cex = cex[2], col = counts.col[2])
        text(90, 274, z2["1100"], cex = cex[3], col = counts.col[2])
        text(95, 100, z2["1010"], cex = cex[2], col = counts.col[2])
        text(200, 43, z2["1001"], cex = cex[3], col = counts.col[2])
        text(200, 280, z2["0110"], cex = cex[2], col = counts.col[2])
        text(300, 100, z2["0101"], cex = cex[2], col = counts.col[2])
        text(310, 274, z2["0011"], cex = cex[3], col = counts.col[2])
        text(130, 219, z2["1110"], cex = cex[2], col = counts.col[2])
        text(245, 71, z2["1101"], cex = cex[3], col = counts.col[2])
        text(155, 72, z2["1011"], cex = cex[3], col = counts.col[2])
        text(270, 219, z2["0111"], cex = cex[2], col = counts.col[2])
        text(200, 140, z2["1111"], cex = cex[2], col = counts.col[2])
        #  text(400, -2, sprintf("Not in any = %i", z2["0000"]), cex = cex[1], col = counts.col[2])
        if (show.include) {
          text(10, 15, include[1], cex = cex[1], col = counts.col[1])
          text(10, -2, include[2], cex = cex[1], col = counts.col[2])
        }
      }
      return(invisible())
    }
    #rect(-20, -30, 430, 430)
    elps <- cbind(150 * cos(seq(0, 2 * pi, len = 1000)), 60 *
                    sin(seq(0, 2 * pi, len = 1000)))
    if (!FILL.COL) {
      polygon(relocate_elp(elps, 90, 200, 250), border = circle.col[1],
              lwd = lwd)
      polygon(relocate_elp(elps, 162, 250, 220), border = circle.col[2],
              lwd = lwd)
      polygon(relocate_elp(elps, 234, 250, 150), border = circle.col[3],
              lwd = lwd)
      polygon(relocate_elp(elps, 306, 180, 125), border = circle.col[4],
              lwd = lwd)
      polygon(relocate_elp(elps, 378, 145, 200), border = circle.col[5],
              lwd = lwd)
    }
    if (FILL.COL) {
      RGB <- col2rgb(circle.col)/255
      ALPHA <- 0.06
      RGB.ALP1 <- rgb(RGB[1, 1], RGB[2, 1], RGB[3, 1], alpha = ALPHA)
      RGB.ALP2 <- rgb(RGB[1, 2], RGB[2, 2], RGB[3, 2], alpha = ALPHA)
      RGB.ALP3 <- rgb(RGB[1, 3], RGB[2, 3], RGB[3, 3], alpha = ALPHA)
      RGB.ALP4 <- rgb(RGB[1, 4], RGB[2, 4], RGB[3, 4], alpha = ALPHA)
      RGB.ALP5 <- rgb(RGB[1, 5], RGB[2, 5], RGB[3, 5], alpha = ALPHA)
      polygon(relocate_elp(elps, 90, 200, 250), border = circle.col[1],
              lwd = lwd, col = RGB.ALP1)
      polygon(relocate_elp(elps, 162, 250, 220), border = circle.col[2],
              lwd = lwd, col = RGB.ALP2)
      polygon(relocate_elp(elps, 234, 250, 150), border = circle.col[3],
              lwd = lwd, col = RGB.ALP3)
      polygon(relocate_elp(elps, 306, 180, 125), border = circle.col[4],
              lwd = lwd, col = RGB.ALP4)
      polygon(relocate_elp(elps, 378, 145, 200), border = circle.col[5],
              lwd = lwd, col = RGB.ALP5)
    }
    text(50, 285, names[1], cex = cex[1])
    text(200, 415, names[2], cex = cex[1])
    text(350, 305, names[3], cex = cex[1])
    text(350, 20, names[4], cex = cex[1])
    text(100, -10, names[5], cex = cex[1])
    text(61, 231, z["10000"], cex = cex[2], col = counts.col[1])
    text(200, 332, z["01000"], cex = cex[2], col = counts.col[1])
    text(321, 248, z["00100"], cex = cex[2], col = counts.col[1])
    text(290, 84, z["00010"], cex = cex[2], col = counts.col[1])
    text(132, 72, z["00001"], cex = cex[2], col = counts.col[1])
    text(146, 253, z["11000"], cex = cex[3], col = counts.col[1])
    text(123, 191, z["10100"], cex = cex[3], col = counts.col[1])
    text(275, 155, z["10010"], cex = cex[3], col = counts.col[1])
    text(137, 149, z["10001"], cex = cex[3], col = counts.col[1])
    text(243, 271, z["01100"], cex = cex[3], col = counts.col[1])
    text(175, 270, z["01010"], cex = cex[3], col = counts.col[1])
    text(187, 120, z["01001"], cex = cex[3], col = counts.col[1])
    text(286, 193, z["00110"], cex = cex[3], col = counts.col[1])
    text(267, 238, z["00101"], cex = cex[3], col = counts.col[1])
    text(228, 108, z["00011"], cex = cex[3], col = counts.col[1])
    text(148, 213, z["11100"], cex = cex[3], col = counts.col[1])
    text(159, 255, z["11010"], cex = cex[3], col = counts.col[1])
    text(171, 144, z["11001"], cex = cex[3], col = counts.col[1])
    text(281, 178, z["10110"], cex = cex[3], col = counts.col[1])
    text(143, 166, z["10101"], cex = cex[3], col = counts.col[1])
    text(252, 148, z["10011"], cex = cex[3], col = counts.col[1])
    text(205, 258, z["01110"], cex = cex[3], col = counts.col[1])
    text(254, 248, z["01101"], cex = cex[3], col = counts.col[1])
    text(211, 121, z["01011"], cex = cex[3], col = counts.col[1])
    text(267, 214, z["00111"], cex = cex[3], col = counts.col[1])
    text(170, 234, z["11110"], cex = cex[3], col = counts.col[1])
    text(158, 172, z["11101"], cex = cex[3], col = counts.col[1])
    text(212, 142, z["11011"], cex = cex[3], col = counts.col[1])
    text(263, 183, z["10111"], cex = cex[3], col = counts.col[1])
    text(239, 235, z["01111"], cex = cex[3], col = counts.col[1])
    text(204, 193, z["11111"], cex = cex[2], col = counts.col[1])
    # text(400, 7, sprintf("Not in any = %i", z["00000"]), cex = cex[1], col = counts.col[1])
    if (length(include) == 2) {
      text(61, 220, z2["10000"], cex = cex[2], col = counts.col[2])
      text(200, 321, z2["01000"], cex = cex[2], col = counts.col[2])
      text(321, 237, z2["00100"], cex = cex[2], col = counts.col[2])
      text(290, 73, z2["00010"], cex = cex[2], col = counts.col[2])
      text(132, 61, z2["00001"], cex = cex[2], col = counts.col[2])
      text(146, 244, z2["11000"], cex = cex[3], col = counts.col[2])
      text(123, 180, z2["10100"], cex = cex[3], col = counts.col[2])
      text(275, 144, z2["10010"], cex = cex[3], col = counts.col[2])
      text(137, 143, z2["10001"], cex = cex[3], col = counts.col[2])
      text(243, 260, z2["01100"], cex = cex[3], col = counts.col[2])
      text(175, 259, z2["01010"], cex = cex[3], col = counts.col[2])
      text(187, 110, z2["01001"], cex = cex[3], col = counts.col[2])
      text(286, 186, z2["00110"], cex = cex[3], col = counts.col[2])
      text(267, 230, z2["00101"], cex = cex[3], col = counts.col[2])
      text(228, 97, z2["00011"], cex = cex[3], col = counts.col[2])
      text(148, 203, z2["11100"], cex = cex[3], col = counts.col[2])
      text(159, 249, z2["11010"], cex = cex[3], col = counts.col[2])
      text(171, 137, z2["11001"], cex = cex[3], col = counts.col[2])
      text(281, 171, z2["10110"], cex = cex[3], col = counts.col[2])
      text(143, 155, z2["10101"], cex = cex[3], col = counts.col[2])
      text(252, 137, z2["10011"], cex = cex[3], col = counts.col[2])
      text(205, 247, z2["01110"], cex = cex[3], col = counts.col[2])
      text(254, 242, z2["01101"], cex = cex[3], col = counts.col[2])
      text(211, 112, z2["01011"], cex = cex[3], col = counts.col[2])
      text(267, 207, z2["00111"], cex = cex[3], col = counts.col[2])
      text(170, 223, z2["11110"], cex = cex[3], col = counts.col[2])
      text(158, 162, z2["11101"], cex = cex[3], col = counts.col[2])
      text(212, 133, z2["11011"], cex = cex[3], col = counts.col[2])
      text(263, 172, z2["10111"], cex = cex[3], col = counts.col[2])
      text(239, 228, z2["01111"], cex = cex[3], col = counts.col[2])
      text(204, 182, z2["11111"], cex = cex[2], col = counts.col[2])
      text(400, -10, sprintf("Not in any = %i", z2["00000"]), cex = cex[1], col = counts.col[2])
      if (show.include) {
        text(10, 7, include[1], cex = cex[1], col = counts.col[1])
        text(10, -10, include[2], cex = cex[1], col = counts.col[2])
      }
    }
    invisible()
  }


taxa_fungi_up[taxa_fungi_up$Genus%in%"Tuber", ]
taxa_fungi_up[!(rownames(taxa_fungi_up) %in% c("FOTU_1","FOTU_3","FOTU_809","FOTU_3308",
                                               "FOTU_1844","FOTU_1987","FOTU_2741","FOTU_2998")), ] 
physeq_fungi_uparse_filt -> physeq_fungi_uparse_ven

tax_table(physeq_fungi_uparse_ven) <- tax_table(as.matrix(
  taxa_fungi_up[!(rownames(taxa_fungi_up) %in% c("FOTU_1","FOTU_3","FOTU_809","FOTU_3308",
                                                 "FOTU_1844","FOTU_1987","FOTU_2741","FOTU_2998")), ]))

tax_table(physeq_fungi_uparse_ven) <-
  tax_table(physeq_fungi_uparse_ven)[,!(colnames(tax_table(physeq_fungi_uparse_ven)) %in% c("OTU_ID"))]

# Fungi
physeq_fungi_uparse_ven %>% 
  tax_glom(taxrank="Genus") %>% #how many shared genera
  subset_samples(Distance%in%c("Inside","0")) %>%
  merge_samples("Type") -> physeq_fungi_0inside

otu_table(physeq_fungi_0inside) <- otu_table(physeq_fungi_0inside)[,which(
  colSums(otu_table(physeq_fungi_0inside)) > 0)] 

otu_fungi_0inside <- as.data.frame(t(otu_table(physeq_fungi_0inside)))
counts_fungi_0inside <- vennCounts(otu_fungi_0inside, include="both")
counts_fungi_0inside

VennGraph(counts_fungi_0inside,
          cex=c(1),
          circle.col =c("brown", "gold"),
          mar = c(1,1,1,1),
          lwd = 2, main="Fungi")

# Prokaryotes
physeq_prok_uparse_filt -> physeq_prok_uparse_ven
tax_table(physeq_prok_uparse_ven) <-
  tax_table(physeq_prok_uparse_ven)[,!(colnames(tax_table(physeq_prok_uparse_ven)) %in% c("OTU_ID"))]

physeq_prok_uparse_ven %>% 
  tax_glom(taxrank="Genus") %>%
  subset_samples(Distance%in%c("Inside","0")) %>%
  merge_samples("Type") -> physeq_prok_0inside

otu_table(physeq_prok_0inside) <- otu_table(physeq_prok_0inside)[, which(
  colSums(otu_table(physeq_prok_0inside)) > 0)] 

otu_prok_0inside <- as.data.frame(t(otu_table(physeq_prok_0inside)))
counts_prok_0inside <- vennCounts(otu_prok_0inside, include="both")
counts_prok_0inside

VennGraph(counts_prok_0inside,
          cex=c(1),
          circle.col =c("brown", "gold"),
          mar = c(1,1,1,1),
          lwd = 2, main="Prokaryotes")

# Extract OTU names for shared and truffle only compartment
library(reshape2)

head(otu_fungi_0inside)
subset(otu_fungi_0inside, truffle>0) -> fungi_core
fungi_core$niche <- ifelse(fungi_core$soil>0, "shared", "truffle")
left_join(rownames_to_column(fungi_core), 
          rownames_to_column(taxa_fungi_up), by = ("rowname" = "rowname")) -> fungi_core
fungi_core$rel_soil <- fungi_core$soil/sum(fungi_core$soil)*100
fungi_core$rel_truffle <- fungi_core$truffle/sum(fungi_core$truffle)*100
head(fungi_core)
melt(fungi_core[, c(11,17:18)]) -> melt_fungi_core


head(otu_prok_0inside)
subset(otu_prok_0inside, truffle>0) -> prok_core
prok_core$niche <- ifelse(prok_core$soil>0, "shared", "truffle")
left_join(rownames_to_column(prok_core), 
          rownames_to_column(taxa_prok_up), by = ("rowname" = "rowname")) -> prok_core
prok_core$rel_soil <- prok_core$soil/sum(prok_core$soil)*100
prok_core$rel_truffle <- prok_core$truffle/sum(prok_core$truffle)*100

# selecting top 30 bacteria
prok_core$sum <- prok_core$rel_soil + prok_core$rel_truffle
prok_core[order(-prok_core$sum),][1:30,] -> prok_core

# correctin taxonomy 
refseq(physeq_prok_uparse_filt)[prok_core$rowname,] 

write.dna(refseq(refseq(physeq_prok_uparse_filt)[prok_core$rowname,]), 
          format="fasta", 
          colsep="",
          file="bacteria_core30.fasta")

#write.csv(prok_core, "prok_core.csv")
read.csv("prok_core_new.csv", sep = ",", header = TRUE, row.names = 1) -> prok_core
head(prok_core)
melt(prok_core[, c(11,17:18)]) -> melt_prok_core

melt_fungi_core
range(melt_fungi_core$value)
melt_fungi_core$value <- ifelse(melt_fungi_core$value==0, NA, paste(melt_fungi_core$value))
melt_fungi_core$value <- as.numeric(as.character(melt_fungi_core$value))

melt_prok_core
range(melt_prok_core$value)
melt_prok_core$value <- ifelse(melt_prok_core$value==0, NA, paste(melt_prok_core$value))
melt_prok_core$value <- as.numeric(as.character(melt_prok_core$value))


# Plot Baloons 
PlotBaloon <- function(df){
  plotb <- ggplot(df, aes(x = variable, y = reorder(Genus, desc(Genus)))) +
    geom_point(aes(size=value, fill=value, color=value), shape=16)+
    #scale_size_area("Abundance %", max_size = 20) +
    #scale_size("Abundance %", range = c(1,20)) +
    scale_x_discrete(labels= c(rel_soil = "Soil", rel_truffle="Truffle")) +
    scale_colour_gradientn(colours = c("#d73027","darkgrey","#4575b4")) +
    theme_classic() +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8)) +
    theme(strip.text.x = element_text(size = 8, face = "bold")) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 0, face = "bold", size = 9, hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, face ="italic", size = 7, hjust = 1, vjust = 0.5)) +
    theme(axis.title.y = element_text(angle = 90, size = 9, face = "bold", hjust = 0.5, vjust = 0.5)) +
    theme(axis.title.x = element_blank()) +
    theme(legend.position="bottom") +
    guides(color = FALSE, fill = FALSE) +
    grids(linetype = "dashed") 
  return(plotb)
}

PlotBaloon(melt_fungi_core) +
  labs(title = "Fungi", x=NULL, y="Genus")+
  scale_size_continuous("Abundance %",range = c(0,15), breaks = c(1, 10, 20, 30, 40))

PlotBaloon(melt_prok_core) + 
  labs(title = "Prokaryotes",  x=NULL, y="Genus") +
  scale_size_continuous("Abundance %",range = c(0,15), breaks = c(1, 10, 20, 30, 40))

get_legend(PlotBaloon(melt_fungi_core) +
             labs(title = "Fungi", x=NULL, y="Genus")+
             scale_size_continuous("Abundance %",range = c(0,15), breaks = c(1, 10, 20, 30, 40))) -> legend_baloon
as_ggplot(legend_baloon)

write.csv(melt_fungi_core, "melt_fungi_core.csv")
write.csv(melt_prok_core, "melt_prok_core.csv")

# plotting ballons
plots_baloon <- ggarrange(PlotBaloon(melt_fungi_core) +
                            labs(title = "Fungi", x=NULL, y="Genus")+
                            scale_size_continuous("Abundance %",range = c(0,15), breaks = c(1, 10, 20, 30, 40)),
                          PlotBaloon(melt_prok_core) + 
                            labs(title = "Prokaryotes",  x=NULL, y="Genus") +
                            scale_size_continuous("Abundance %",range = c(0,15), breaks = c(1, 10, 20, 30, 40)),
                          labels = c("A", "B"),
                          widths = c(1,1),
                          #align = "hv", 
                          ncol = 2, 
                          nrow = 1,
                          common.legend = TRUE,
                          legend = "bottom")

plots_baloon

# *** FIGURE 5 -Balloons ---------------------------------------------------------------------------
# Puttig two ggplot2 graphs with diffeent size together using cowplot ------------------------------
library(cowplot)

ggdraw() +
  draw_plot(PlotBaloon(melt_fungi_core) +
              labs(title = "Fungi", x=NULL, y="Genus")+
              scale_size_continuous("Abundance %",range = c(0,15), breaks = c(1, 10, 20, 30, 40)) +
              theme(legend.position="none"), 
            x = 0, y = 0.24, width = 0.45, height = 0.6) +
  draw_plot(PlotBaloon(melt_prok_core) + 
              labs(title = "Prokaryotes",  x=NULL, y="Genus") +
              scale_size_continuous("Abundance %",range = c(0,15), breaks = c(1, 10, 20, 30, 40)) +
              theme(legend.position="none"), 
            x = 0.45, y = 0.1, width = 0.5, height = 0.74) +
  draw_plot(as_ggplot(legend_baloon), 
            x = 0, y = 0, width = 1, height = 0.1) +
  draw_plot_label(c("A", "B"), x = c(0, 0.45), y = c(0.84, 0.84), size = 12)

