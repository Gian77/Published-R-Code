# Manuscript:   Contrasting effects of bioenergy crops on biodiversity
# Authors:  Nathan L. Haan1,2,4, Gian N.M. Benucci2,3, Cynthia M. Fiser1, Gregory Bonito2,3, and Douglas A. Landis1,2
# Affiliation:  1 Department of Entomology, Michigan State University, East Lansing MI 48824
#               2 Great Lakes Bioenergy Research Center, Michigan State University, East Lansing MI 48824
#               3 Department of Plant Soil and Microbial Sciences, Michigan State University, East Lansing MI 48824
# Journal:      PNAS 
# Date:         Dec 19, 2022

# WORKING ENVIRONMENT SETUP ----------------------------------------------------
options(scipen = 999)
options(max.print = 100000000)

# R PACKAGES -------------------------------------------------------------------
library(styler) 
library(evoPalette)
library(phyloseq)
library(Biostrings)
library(ape)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(decontam)

library(CNVRG) # library normalization
library(mikropml)
library(mbImpute)

devtools::install_github("Gian77/Misfunk")
library(Misfunk)

# *******************************************************************-----------
# IMPORT DATASETS --------------------------------------------------------------
mapping_ITS <-
  read.delim("otu_tables/mapping_ITS.txt", row.names=1, header=TRUE, sep="\t")
mapping_ITS
table_ITS <- 
  read.delim("otu_tables/otu_table_ITS_UPARSE.txt",row.names=1, header=TRUE, sep="\t") 
head(table_ITS)
constax_ITS <-
  read.delim("otu_tables/constax_taxonomy_ITS.txt",header=TRUE, row.names=1, sep="\t")
head(constax_ITS)
otu_ITS <- 
  readDNAStringSet("otu_tables/otus_ITS.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
otu_ITS

mapping_16s <-
  read.delim("otu_tables/mapping_16s.txt", row.names=1, header=TRUE, sep="\t")
head(mapping_16s)
table_16s <- 
  read.delim("otu_tables/otu_table_16S_UPARSE.txt",row.names=1, header=TRUE, sep="\t") 
head(table_16s)
constax_16s <-
  read.delim("otu_tables/constax_taxonomy_16S.txt",header=TRUE, row.names=1, sep="\t")
head(constax_16s)
otu_16s <- 
  readDNAStringSet("otu_tables/otus_16s.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
otu_16s

constax_16s <- constax_16s[, 1:7]
colnames(constax_16s) <- 
  c("Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species")


mapping_18s <-
  read.delim("otu_tables/mapping_18s.txt", row.names=1, header=TRUE, sep="\t")
head(mapping_18s)
table_18s <- 
  read.delim("otu_tables/otu_table_18S_UPARSE.txt",row.names=1, header=TRUE, sep="\t") 
head(table_18s)
constax_18s <-
  read.delim("otu_tables/constax_taxonomy_18S.txt",header=TRUE, row.names=1, sep="\t")
head(constax_18s)
otu_18s <- 
  readDNAStringSet("otu_tables/otus_18s.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
otu_18s


# Filter taxonomy UPARSE -------------------------------------------------------
CleanTax <- function(tax){
  tax[, "Kingdom"] <- as.factor(gsub("_1", "", tax[, "Kingdom"]))
  tax[, "Phylum"] <- as.factor(gsub("_1", "", tax[, "Phylum"]))
  tax[, "Class"] <- as.factor(gsub("_1", "", tax[, "Class"]))
  tax[, "Order"] <- as.factor(gsub("_1", "",tax[, "Order"]))
  tax[, "Family"] <- as.factor(gsub("_1", "", tax[, "Family"]))
  tax[, "Genus"] <- as.factor(gsub("_1", "", tax[, "Genus"]))
  tax[, "Species"] <- as.factor(gsub("_1", "", tax[, "Species"]))
  print("High Level taxonomic groups in CONSTAX")
  print(table(tax$High_level_taxonomy))
  print("Taxonomic groups at Kingdom level")
  print(table(tax$Kingdom))
  return(tax)
}

constax_ITS_clean <- 
  CleanTax(constax_ITS)
constax_16s_clean <- 
  CleanTax(constax_16s)
constax_18s_clean <- 
  CleanTax(constax_18s)


HighLevel <- function(tax, group){
high_level <-
  unique(
    c(
      names(table(tax$High_level_taxonomy)),
      names(table(tax$Kingdom)))
    )
high_level <- 
  high_level[!high_level %in% c(group)]
return(high_level)
}

high_level_ITS <-
  HighLevel(constax_ITS_clean, "Fungi") 
high_level_ITS

high_level_18s <-
  HighLevel(constax_18s_clean, "") 
high_level_18s



# ITS --------------------------------------------------------------------------
GetBadTaxa <- function(tax, hl){
# Remove non-target taxa and non-classified taxa that did not have any hit
# when blasted against UNITE at 60% conevrage and identity. ALso, remove 
# fungal OTUs that were classified as Fungi bu had low hit coverage < 60%
# and identity < 60% 
bad_taxa <- 
  subset(
    tax,
    High_level_taxonomy%in%hl |
      Kingdom%in%hl |
      High_level_taxonomy%in%c("") &
      Kingdom%in%c("") |
      Kingdom%in%c("Fungi") &
      HL_hit_query_cover<60 |
      HL_hit_percent_id<60)
print(str(bad_taxa))
return(bad_taxa)
}

bad_taxa_ITS <-
  GetBadTaxa(constax_ITS_clean, high_level_ITS)
write.csv(bad_taxa_ITS, "bad_taxa_ITS.csv")
str(bad_taxa_ITS)

table(bad_taxa_ITS$High_level_taxonomy)

constax_ITS_clean <-
  constax_ITS_clean[!(rownames(constax_ITS_clean) %in% rownames(bad_taxa_ITS)), ]
dim(constax_ITS_clean)
head(constax_ITS_clean)

table(constax_ITS$High_level_taxonomy)
table(constax_ITS$Kingdom)
table(constax_ITS_clean$High_level_taxonomy)
table(constax_ITS_clean$Kingdom)

# 18S --------------------------------------------------------------------------
constax_18s_clean <- 
  subset(
    constax_18s_clean,
    High_level_taxonomy%in%high_level_18s |
      Kingdom%in%high_level_18s & 
      HL_hit_query_cover>60 |
      HL_hit_percent_id>60)

str(constax_18s_clean)
head(constax_18s_clean)

table(constax_18s$High_level_taxonomy)
table(constax_18s$Kingdom)
table(constax_18s_clean$High_level_taxonomy)
table(constax_18s_clean$Kingdom)
constax_18s[constax_18s$Kingdom%in%"Planomonada_1",]
nrow(constax_18s_clean[constax_18s_clean$High_level_taxonomy%in%c(""),])

# 16S --------------------------------------------------------------------------
any(constax_16s_clean$Kingdom == "Chloroplast")
any(constax_16s_clean$Kingdom == "Mitochondria")
any(constax_16s_clean$Phylum == "Chloroplast") 
any(constax_16s_clean$Phylum == "Mitochondria") #TRUE
any(constax_16s_clean$Class == "Chloroplast")
any(constax_16s_clean$Class == "Mitochondria")
any(constax_16s_clean$Order == "Chloroplast") #TRUE
any(constax_16s_clean$Order == "Mitochondria")
any(constax_16s_clean$Family == "Chloroplast") 
any(constax_16s_clean$Family == "Mitochondria")#TRUE
any(constax_16s_clean$Genus == "Chloroplast")
any(constax_16s_clean$Genus == "Mitochondria")

unwanted_amp <- c("Mitochondria", "Chloroplast", "mitochondria", "chloroplast")
apply(constax_16s_clean, 2, function(x) which(x %in% unwanted_amp))

constax_16s_clean %>% dplyr::filter(constax_16s_clean$Phylum%in%unwanted_amp)
constax_16s_clean %>% dplyr::filter(constax_16s_clean$Order%in%unwanted_amp)
constax_16s_clean %>% dplyr::filter(constax_16s_clean$Family%in%unwanted_amp)
str(constax_16s_clean)

bad_taxa_16s <- 
  subset(constax_16s_clean, 
         Phylum %in% "Mitochondria" | 
         Order %in% "Chloroplast" |
         Family %in% "Mitochondria")
str(bad_taxa_16s)
bad_taxa_16s

constax_16s_clean <-
  constax_16s_clean[!(rownames(constax_16s_clean) %in% rownames(bad_taxa_16s)), ]
dim(constax_16s_clean)
head(constax_16s_clean)

# Reformat taxonomy ------------------------------------------------------------
source("../R_functions/ReformatTaxonomy.R")

FinalizeTaxonomy <- function(constax){
  constax <- data.frame(OTU_ID = rownames(constax), constax)
  #constax[constax == "Incertae_sedis" ] <- "" # remove Incertae sedis
  constax[, "Species"] <- gsub(" sp", "", constax[, "Species"]) # remove sp at Species level
  constax[] = lapply(constax, blank2na, na.strings=c('','NA','na','N/A','n/a','NaN','nan'))
  lastValue <- function(x) tail(x[!is.na(x)], 1)
  last_taxons<- apply(constax[,1:8], 1, lastValue)
  constax$BestMatch <- last_taxons
  constax[, "BestMatch"] <-
    gsub("_", " ", constax[, "BestMatch"])
  constax$Taxonomy <-
    paste(constax$OTU_ID, constax$BestMatch, sep = "-")
  constax[, "Genus"] <- gsub(" sp.", "", constax[, "Genus"])
  return(constax)
}

constax_16s_cor <- 
  FinalizeTaxonomy(constax_16s_clean)
constax_16s_cor[1:20,]

constax_ITS_cor <- 
  FinalizeTaxonomy(constax_ITS_clean)
constax_ITS_cor[1:20,]

constax_18s_cor <- 
  FinalizeTaxonomy(constax_18s_clean)
constax_18s_cor[1:20,]


# **************************************************************************----
# *********NEW MODIFIED CODE STARTS ****************----------------------------
# **************************************************************************----

# Generate phyloseq objects ----------------------------------------------------
physeq_ITS <-
  phyloseq(
    otu_table(table_ITS, taxa_are_rows = TRUE),
    sample_data(mapping_ITS),
    tax_table(as.matrix(constax_ITS_cor)),
    otu_ITS)

physeq_ITS
head(physeq_ITS@sam_data)
physeq_ITS@sam_data$Niche

identical(
  rownames(data.frame(LibSize = colSums(physeq_ITS@otu_table))),
  rownames(as(physeq_ITS@sam_data, "data.frame")))
physeq_ITS@sam_data$LibSize <- 
  colSums(physeq_ITS@otu_table)


physeq_16s <-
  phyloseq(
    otu_table(table_16s, taxa_are_rows = TRUE),
    sample_data(mapping_16s),
    tax_table(as.matrix(constax_16s_cor)),
    otu_16s)

physeq_16s
head(physeq_16s@sam_data)
physeq_16s@sam_data$Niche
identical(
  rownames(data.frame(LibSize = colSums(physeq_16s@otu_table))),
  rownames(as(physeq_16s@sam_data, "data.frame")))
physeq_16s@sam_data$LibSize <- 
  colSums(physeq_16s@otu_table)

# remove samples with low reads that have been reamplified
subset_samples(physeq_16s, 
        Description%in%c("LFG1r4","LFG6r5","LFG7r2",
        "LFG1r4_rec","LFG6r5_rec","LFG7r2_rec"))@sam_data

physeq_16s <-
  subset_samples(physeq_16s, 
        !Description%in%c("LFG1r4","LFG6r5","LFG7r2_rec"))

physeq_18s <-
  phyloseq(
    otu_table(table_18s, taxa_are_rows = TRUE),
    sample_data(mapping_18s),
    tax_table(as.matrix(constax_18s_cor)),
    otu_18s)

physeq_18s
head(physeq_18s@sam_data)
physeq_18s@sam_data$Niche
identical(
  rownames(data.frame(LibSize = colSums(physeq_18s@otu_table))),
  rownames(as(physeq_18s@sam_data, "data.frame")))
physeq_18s@sam_data$LibSize <- 
  colSums(physeq_18s@otu_table)



# DECONTAMINATION --------------------------------------------------------------
physeq_ITS@sam_data

sample_data(physeq_ITS)$is.neg <-
  sample_data(physeq_ITS)$Niche == "Control"

sample_data(physeq_16s)$is.neg <-
  sample_data(physeq_16s)$Niche == "Control"

sample_data(physeq_18s)$is.neg <-
  sample_data(physeq_18s)$Niche == "Control"

contam_ITS <-
  decontam::isContaminant(
    physeq_ITS,
    method = "prevalence",
    batch = "Plate",
    neg = "is.neg",
    threshold = 0.5)

table(contam_ITS$contaminant) 


contam_16s <-
  decontam::isContaminant(
    physeq_16s,
    method = "prevalence",
    batch = "Plate",
    neg = "is.neg",
    threshold = 0.5)

table(contam_16s$contaminant) 


contam_18s <-
  decontam::isContaminant(
    physeq_18s,
    method = "prevalence",
    batch = "Plate",
    neg = "is.neg",
    threshold = 0.5)

table(contam_18s$contaminant) 

# plotting contaminant OTUs ----------------------------------------------------
PlotContam <- function(df, contam){
  # Make phyloseq object of presence-absence in negative controls and true samples
  physeq_pa <- transform_sample_counts(df, function(abund) 1*(abund>0))
  physeq_pa_neg <- subset_samples(physeq_pa, is.neg%in%c("TRUE"))
  physeq_pa_pos <- subset_samples(physeq_pa, is.neg%in%c("FALSE"))
  # Make data.frame of prevalence in positive and negative samples
  df_contam <- data.frame(pa.pos=taxa_sums(physeq_pa_pos), 
                          pa.neg=taxa_sums(physeq_pa_neg),
                          contaminant=contam$contaminant, 
                          Pvalue=contam$p)
  head(df_contam) %T>% print()
  # plotting 
  ggplot(data=df_contam, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
    geom_point(size=2, alpha=0.7) +
    labs(x="Prevalence in negative controls", y="Prevalence in true samples") +
    theme_classic() +
    scale_colour_manual("Contaminant OTUs", values = c("grey", "red")) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.x = element_text(angle =0, size = 8, hjust = 0.5, vjust = 1), 
          axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8)) -> plot_cont
  return(plot_cont)
}

PlotContam(physeq_ITS, contam_ITS)
PlotContam(physeq_16s, contam_16s) 
PlotContam(physeq_18s, contam_18s) 

sample_data(physeq_ITS)$Index <- as.numeric(as.character(seq(nrow(sample_data(physeq_ITS)))))
sample_data(physeq_16s)$Index <- as.numeric(as.character(seq(nrow(sample_data(physeq_16s)))))
sample_data(physeq_18s)$Index <- as.numeric(as.character(seq(nrow(sample_data(physeq_18s)))))

# Function to plot sample depth ------------------------------------------------
PlotDepth <- function(physeq){
  df <- as(sample_data(physeq), "matrix")
  df <- as.data.frame(df)
  # reconvert to numeric
  df$LibSize <- as.numeric(as.character(df$LibSize))
  df$Index <- as.numeric(as.character(df$Index))
  # order
  df <- df[order(df$LibSize), ]
  df$Index <- seq(nrow(df))
  # inspect
  str(df) %T>% print()
  head(df) %T>% print()
  ggplot(data=df, aes(x=Index, y=LibSize, color=is.neg)) +
    geom_point(alpha =0.7, size=2) +
    theme_classic() +
    scale_colour_manual("Negative control", values = c("grey", "red")) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.x = element_text(angle =0, size = 8, hjust = 0.5, vjust = 1), 
          axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8)) -> plot_dist
  return(plot_dist)  
}

PlotDepth(physeq_ITS)



# *** FIGURE S1 - dcontam ------------------------------------------------------
ggarrange(
  ggarrange(PlotDepth(physeq_ITS) +
              labs(title="ITS", 
                   subtitle = "Samples read depth", 
                   x="Sample index", 
                   y="Read number"),
            PlotDepth(physeq_16s) +
              labs(title="16S", 
                   subtitle = "Samples read depth",
                   x="Sample index",
                   y="Read number"),
            PlotDepth(physeq_18s) +
              labs(title="18S", 
                   subtitle = "Samples read depth",
                   x="Sample index",
                   y="Read number"),
            labels = c("A","B","C"),
            widths = c(1,1,1),
            align = "hv" ,
            ncol = 3, 
            nrow = 1, 
            common.legend = TRUE, 
            legend = c("bottom")),
  ggarrange(
    PlotContam(physeq_ITS, contam_ITS) +
      labs(subtitle="Contaminants n. 13"),
    PlotContam(physeq_16s, contam_16s) +
      labs(subtitle="Contaminants n. 38"),
    PlotContam(physeq_18s, contam_18s) +
      labs(subtitle="Contaminants n. 40"),
    widths = c(1,1,1),
    labels = c("C","D","E"),
    align = "hv" ,
    ncol = 3, 
    nrow = 1,
    common.legend = TRUE,
    legend = c("bottom")),
  widths =  c(1, 1.2),
  ncol = 1, 
  nrow = 2) -> Fig_S1

Fig_S1


# ***************************************************---------------------------
# REMOVING CONTAMINANTS --------------------------------------------------------
# function to remove bad taxa
remove_taxa = function(physeq, badTaxa) {
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

physeq_ITS_clean <-
  remove_taxa(physeq_ITS, rownames(subset(
    contam_ITS, contaminant %in% c("TRUE")
  )))

physeq_16s_clean <-
  remove_taxa(physeq_16s, rownames(subset(
    contam_16s, contaminant %in% c("TRUE")
  )))

physeq_18s_clean <-
  remove_taxa(physeq_18s, rownames(subset(
    contam_18s, contaminant %in% c("TRUE")
  )))


# REMOVE CONTROLS and non-POPLAR LEAF SAMPLES ----------------------------------
physeq_ITS_filt <-
  subset_samples(physeq_ITS_clean, is.neg %in% c("FALSE"))
otu_table(physeq_ITS_filt) <-
  otu_table(physeq_ITS_filt)[which(rowSums(otu_table(physeq_ITS_filt)) > 0), ]
otu_table(physeq_ITS_filt) <-
  otu_table(physeq_ITS_filt)[which(colSums(otu_table(physeq_ITS_filt)) > 0), ]

physeq_ITS_filt@sam_data

physeq_ITS_filt@sam_data$Niche

physeq_16s_filt <-
  subset_samples(physeq_16s_clean, is.neg %in% c("FALSE"))
otu_table(physeq_16s_filt) <-
  otu_table(physeq_16s_filt)[which(rowSums(otu_table(physeq_16s_filt)) > 0), ]
physeq_16s_filt@sam_data

physeq_18s_filt <-
  subset_samples(physeq_18s_clean, is.neg %in% c("FALSE"))
otu_table(physeq_18s_filt) <-
  otu_table(physeq_18s_filt)[which(rowSums(otu_table(physeq_18s_filt)) > 0), ]
physeq_18s_filt@sam_data

# RAREFACTION CURVES -----------------------------------------------------------
# RAREFACTION CURVES -----------------------------------------------------------
# https://stat.ethz.ch/pipermail/r-sig-ecology/2018-December/005867.html
# https://community.rstudio.com/t/problem-with-names-function/66198
# https://github.com/riffomonas/distances/blob/7945a48b9f404242f3ae3ca28624df31ddf7a422/code/vegan_rarefy.R

library(vegan)

rarecurve(t(as.data.frame(otu_table(physeq_ITS_filt))), step = 1000, label = FALSE) -> rare_its
rarecurve(t(as.data.frame(otu_table(physeq_16s_filt))), step = 1000, label = FALSE) -> rare_16s
rarecurve(t(as.data.frame(otu_table(physeq_18s_filt))), step = 1000, label = FALSE) -> rare_18s

ExtrRare <- function(rare_obj, physeq){
  require(tidyverse)
  # The dataframe has to have a column with SampleID
  if ("SampleID" %in% colnames(sample_data(physeq))) {
    cat("\nSampleID variable is present. Good to go!!\n")
  } else {
    cat("\nSampleID variable is absent. Creating it from rownames.\n")
    sample_data(physeq)$SampleID <-
      rownames(sample_data(physeq))
  }
  rare_joint <-
    map_dfr(rare_obj, bind_rows) %>%
    bind_cols(SampleID = rownames(sample_data(physeq)), .) %>%
    pivot_longer(-SampleID) %>%
    drop_na() %>%
    mutate(n_seqs = as.numeric(str_replace(name, "N", ""))) %>%
    left_join(x = .,
              y = as.data.frame(sample_data(physeq)),
              by = "SampleID")
  return(rare_joint)
}


ExtrRare(rare_its, physeq_ITS_filt)

# Then plotting rarecurve in ggplot2
PlotRareCurve <- function(rare_joint){
  ggplot(rare_joint, aes(x=n_seqs, y=value, group=SampleID, color=Niche)) +  
    geom_line() +
    scale_colour_manual(values = c("darkgreen", "green", "brown")) +
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 9, face = "bold")) + 
    theme(axis.text.x = element_text(angle = 0, size = 9,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 9, hjust = 0.5, vjust = 0.5)) +
    theme(legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.7, "cm")) +
    theme(legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 10)) +
    grids(linetype = "dashed") -> plot_rare
  return(plot_rare)
}


PlotRareCurve(ExtrRare(rare_its, physeq_ITS_filt)) +
  labs(title="ITS", x= "Number of DNA reads", y= "Number of OTUs")

# *** FIGURE S2 - Rarefaction curves **********************************---------
ggarrange(PlotRareCurve(ExtrRare(rare_its, physeq_ITS_filt)) +
            labs(title="ITS", x= "Number of DNA sequences", y= "Number of OTUs"),
          PlotRareCurve(ExtrRare(rare_16s, physeq_16s_filt)) +
            labs(title="16S", x= "Number of DNA sequences", y= "Number of OTUs"),
          PlotRareCurve(ExtrRare(rare_18s, physeq_18s_filt)) +
            labs(title="ITS", x= "Number of DNA sequences", y= "Number of OTUs"),
          heights = c(1, 1, 1),
          align = "hv",
          ncol = 3, 
          nrow = 1,
          common.legend = TRUE, 
          legend = "bottom") -> Fig_S2_rare

Fig_S2_rare


# NORMALIZATION ----------------------------------------------------------------
set.seed(20221312)

# ITS region
sort(colSums(otu_table(physeq_ITS_filt)), decreasing = TRUE)

physeq_ITS_ev <- 
  rarefy_even_depth(physeq_ITS_filt, 
                    sample.size = 12774, 
                    rngseed = FALSE, 
                    replace = TRUE, 
                    trimOTUs = TRUE, 
                    verbose = TRUE) 
otu_table(physeq_ITS_ev) <- 
  prune_taxa(taxa_sums(physeq_ITS_ev) > 0, physeq_ITS_ev) 

colSums(otu_table(physeq_ITS_ev))
any(taxa_sums(physeq_ITS_ev) == 0)
physeq_ITS_ev@sam_data$Niche

# 16S region
sort(colSums(otu_table(physeq_16s_filt)), decreasing = TRUE)

physeq_16s_ev <- 
  rarefy_even_depth(physeq_16s_filt, 
                    sample.size = 1000, 
                    rngseed = FALSE, 
                    replace = TRUE, 
                    trimOTUs = TRUE, 
                    verbose = TRUE) 
otu_table(physeq_16s_ev) <- 
  prune_taxa(taxa_sums(physeq_16s_ev) > 0, physeq_16s_ev) 

colSums(otu_table(physeq_16s_ev))
any(taxa_sums(physeq_16s_ev) == 0)


# 18S region
sort(colSums(otu_table(physeq_18s_filt)), decreasing = TRUE)

physeq_18s_ev <- 
  rarefy_even_depth(physeq_18s_filt, 
                    sample.size = 100, 
                    rngseed = FALSE, 
                    replace = TRUE, 
                    trimOTUs = TRUE, 
                    verbose = TRUE) 
otu_table(physeq_18s_ev) <- 
  prune_taxa(taxa_sums(physeq_18s_ev) > 0, physeq_18s_ev) 

colSums(otu_table(physeq_18s_ev))
any(taxa_sums(physeq_18s_ev) == 0)


# Writing tables
write.csv(physeq_ITS_mSeq@otu_table, "fungi_mSeq.csv")
write.csv(physeq_ITS_mSeq@tax_table, "fungi_tax_mSeq.csv")
write.csv(physeq_ITS_mSeq@sam_data, "fungi_meta_mSeq.csv")

write.csv(physeq_16s_mSeq@otu_table, "bacteria_mSeq.csv")
write.csv(physeq_16s_mSeq@tax_table, "bacteria_tax_mSeq.csv")
write.csv(physeq_16s_mSeq@sam_data, "bacteria_meta_mSeq.csv")



# TESITNG CSS TRANSFORMATION ---------------------------------------------------
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

physeq_ITS_mSeq <-
  CSSNorm(physeq_ITS_filt) 
head(otu_table(physeq_ITS_mSeq))

physeq_16s_mSeq <-
  CSSNorm(physeq_16s_filt) 
head(otu_table(physeq_16s_mSeq))

physeq_18s_mSeq <-
  CSSNorm(physeq_18s_filt) 
head(otu_table(physeq_18s_mSeq))


# TEST quick ordination --------------------------------------------------------
pcoa_ITS <- ordinate(physeq_ITS, method ="PCoA", distance="bray") 
pcoa_ITS_ev <- ordinate(physeq_ITS_ev, method ="PCoA", distance="bray") 
pcoa_ITS_mSeq <- ordinate(physeq_ITS_mSeq, method ="PCoA", distance="bray") 

p_its <-
ggarrange(
  plot_ordination(physeq_ITS, pcoa_ITS, type="samples", 
                  color="Niche", shape = "Niche", title="Raw") +
    scale_colour_manual(values = c("black", "darkgreen", "green", "brown")) + 
    scale_shape_manual(values = c(15, 16, 17, 18)) +
    theme_classic() + theme(legend.position =  "right"),
  ggarrange(
      plot_ordination(physeq_ITS_ev, pcoa_ITS_ev, type="samples", 
                      color="Niche",shape = "Niche", title="Rarefaction")+
        scale_colour_manual(values = c("darkgreen", "green", "brown")) + 
        theme_classic() + 
        scale_shape_manual(values = c(16, 17, 18)),
      plot_ordination(physeq_ITS_mSeq, pcoa_ITS_mSeq, type="samples", 
                      color="Niche", shape = "Niche", title="metagenomeSeq")+
        scale_colour_manual(values = c("darkgreen", "green", "brown"))+ 
        theme_classic() + 
        scale_shape_manual(values = c(16, 17, 18)) ,
      ncol = 2, 
      nrow=1, 
      common.legend = TRUE, 
      legend = "right"),
  ncol = 2,
  nrow = 1,
  widths = c(0.55, 0.9),
  common.legend = FALSE)

p_its

pcoa_16s <- ordinate(physeq_16s, method ="PCoA", distance="bray") 
pcoa_16s_ev <- ordinate(physeq_16s_ev, method ="PCoA", distance="bray") 
pcoa_16s_mSeq <- ordinate(physeq_16s_mSeq, method ="PCoA", distance="bray") 

p_16s <-
  ggarrange(
    plot_ordination(physeq_16s, pcoa_16s, type="samples", 
                    color="Niche", shape = "Niche", title="Raw") +
      scale_colour_manual(values = c("black", "darkgreen", "green", "brown")) + 
      scale_shape_manual(values = c(15, 16, 17, 18)) +
      theme_classic() + theme(legend.position =  "right"),
    ggarrange(
      plot_ordination(physeq_16s_ev, pcoa_16s_ev, type="samples", 
                      color="Niche",shape = "Niche", title="Rarefaction")+
        scale_colour_manual(values = c("darkgreen", "green", "brown")) + 
        theme_classic() + 
        scale_shape_manual(values = c(16, 17, 18)),
      plot_ordination(physeq_16s_mSeq, pcoa_16s_mSeq, type="samples", 
                      color="Niche", shape = "Niche", title="metagenomeSeq")+
        scale_colour_manual(values = c("darkgreen", "green", "brown"))+ 
        theme_classic() + 
        scale_shape_manual(values = c(16, 17, 18)) ,
      ncol = 2, 
      nrow=1, 
      common.legend = TRUE, 
      legend = "right"),
    ncol = 2,
    nrow = 1,
    widths = c(0.55, 0.9),
    common.legend = FALSE)

p_16s


pcoa_18s <- ordinate(physeq_18s, method ="PCoA", distance="bray") 
pcoa_18s_ev <- ordinate(physeq_18s_ev, method ="PCoA", distance="bray") 
pcoa_18s_mSeq <- ordinate(physeq_18s_mSeq, method ="PCoA", distance="bray") 

p_18s <-
  ggarrange(
    plot_ordination(physeq_18s, pcoa_18s, type="samples", 
                    color="Niche", shape = "Niche", title="Raw") +
      scale_colour_manual(values = c("black", "darkgreen", "green", "brown")) + 
      scale_shape_manual(values = c(15, 16, 17, 18)) +
      theme_classic() + theme(legend.position =  "right"),
    ggarrange(
      plot_ordination(physeq_18s_ev, pcoa_18s_ev, type="samples", 
                      color="Niche",shape = "Niche", title="Rarefaction")+
        scale_colour_manual(values = c("darkgreen", "green", "brown")) + 
        theme_classic() + 
        scale_shape_manual(values = c(16, 17, 18)),
      plot_ordination(physeq_18s_mSeq, pcoa_18s_mSeq, type="samples", 
                      color="Niche", shape = "Niche", title="metagenomeSeq")+
        scale_colour_manual(values = c("darkgreen", "green", "brown"))+ 
        theme_classic() + 
        scale_shape_manual(values = c(16, 17, 18)) ,
      ncol = 2, 
      nrow=1, 
      common.legend = TRUE, 
      legend = "right"),
    ncol = 2,
    nrow = 1,
    widths = c(0.55, 0.9),
    common.legend = FALSE)

p_18s


library(gridExtra)
library(grid)

title1=text_grob("ITS",size=12, face=2)
title2=text_grob("16S",size=12, face=2)
title3=text_grob("18S",size=12, face=2)

ggarrange(
    grid.arrange(p_its, top = title1),
    grid.arrange(p_16s, top = title2),
    grid.arrange(p_18s, top = title3),
    ncol = 1,
    nrow = 3)

grid.arrange(
  grid.arrange(p_its, top = title1),
  grid.arrange(p_16s, top = title2),
  grid.arrange(p_18s, top = title3))



# ***************************************************---------------------------

# Subsetting to only microinvertebrates. Animal groups will be left out for 
# another project.
# Separating Kingdoms and checking assignments through NCBI BLAST --------------
table(data.frame(
  as(physeq_18s_filt@tax_table, "matrix"))$High_level_taxonomy)

table(data.frame(
  as(physeq_18s_filt@tax_table, "matrix"))$Kingdom)

# Select taxa form phyloseq
selectTaxa = function(physeq, goodTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[(allTaxa %in% goodTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

# Filter taxa form phyloseq
filterTaxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}


# Alveolata --------------------------------------------------------------------
physeq_18s_Alve <-
  physeq_18s_mSeq %>%
  subset_taxa(High_level_taxonomy%in%"Alveolata" | Kingdom%in%"Alveolata")

physeq_18s_Alve
physeq_18s_Alve@tax_table

# Write sequences for BLAST 
write.dna(
  refseq(physeq_18s_Alve),
  format="fasta", 
  colsep="", 
  file="otus_Alveolata.fasta")

# filtering amoebozoa,unidentified and dubious
bad_otus_Alve <- 
  c("OTU_3448","OTU_8298","OTU_441", "OTU_4675", "OTU_1021",
    "OTU_5678", "OTU_5010","OTU_830", "OTU_8013", "OTU_7833",
    "OTU_1138", "OTU_7934", "OTU_10649", "OTU_6513", "OTU_10228",
    "OTU_4816","OTU_2988" ,"OTU_9808", "OTU_7776", "OTU_1059",
    "OTU_6496", "OTU_1325", "OTU_4136", "OTU_9199", "OTU_6212",
    "OTU_10198", "OTU_6706","OTU_4415", "OTU_5912","OTU_2470",
    "OTU_5712",  "OTU_7383","OTU_6297"
  ) 

rhizaria_otus <-
  c("OTU_1403","OTU_111","OTU_3815", "OTU_566", "OTU_5351",
    "OTU_776","OTU_340","OTU_1792","OTU_705", "OTU_2970",
    "OTU_2247","OTU_1622", "OTU_2977", "OTU_4151", "OTU_2712", 
    "OTU_3297","OTU_10558", "OTU_3815","OTU_6112", "OTU_6319", 
    "OTU_6380", "OTU_6999"
  )

# Filtering taxa
physeq_18s_Alve_filt <- 
  filterTaxa(physeq_18s_Alve, c(bad_otus_Alve, rhizaria_otus))
physeq_18s_Alve_filt

# Select taxa
physeq_18s_Alve_filt_2 <- 
  selectTaxa(physeq_18s_Metazoa, alveolata_otus_1)
physeq_18s_Alve_filt_2

# Merge phyloseq objects
physeq_18s_Alve_filt <- 
  merge_phyloseq(physeq_18s_Alve_filt, physeq_18s_Alve_filt_2)
physeq_18s_Alve_filt


# Protista ---------------------------------------------------------------------
physeq_18s_Prot <-
  physeq_18s_mSeq %>%
  subset_taxa(High_level_taxonomy%in%"Protista" | Kingdom%in%"Protista")
physeq_18s_Prot

write.dna(
  refseq(physeq_18s_Prot),
  format="fasta", 
  colsep="", 
  file="otus_Protista.fasta")


bad_protista_otus <-
  c("OTU_1915", "OTU_2935","OTU_275", "OTU_2990", "OTU_3835",
    "OTU_3678", "OTU_5356", "OTU_4526","OTU_621", "OTU_8139",
    "OTU_9883", "OTU_6372", "OTU_5390", "OTU_10739", "OTU_6633",
    "OTU_7294", "OTU_8651"
  )

stramenopiles_otus_2 <- c("OTU_6403", "OTU_10212","OTU_3792", "OTU_4583", "OTU_10951",
                          "OTU_9888", "OTU_3177", "OTU_8418","OTU_10620", "OTU_9253"
)

# Filtering taxa
physeq_18s_Prot_filt <- 
  filterTaxa(physeq_18s_Prot, c(bad_protista_otus, stramenopiles_otus_2))

# Stramenopila -----------------------------------------------------------------
physeq_18s_Stra <-
  physeq_18s_mSeq %>%
  subset_taxa(High_level_taxonomy%in%"Stramenopila" | Kingdom%in%"Stramenopila")
physeq_18s_Stra

write.dna(
  refseq(physeq_18s_Stra),
  format="fasta", 
  colsep="", 
  file="otus_Stramenopila.fasta")


bad_stramenopiles_otus <- 
  c("OTU_1414", "OTU_1437", "OTU_9102", "OTU_5112", "OTU_6203",
    "OTU_5689", "OTU_942","OTU_2344", "OTU_2808", "OTU_2374", "
    OTU_6741","OTU_3847","OTU_2734", "OTU_6376", "OTU_1818",
    "OTU_7224", "OTU_2979", "OTU_5886", "OTU_7887", "OTU_6488",
    "OTU_4139", "OTU_6968", "OTU_7294","OTU_6718","OTU_8007",
    "OTU_10097","OTU_712","OTU_2037", "OTU_1723")


stramenopiles_otus_2

# Filtering taxa
physeq_18s_Stra_filt <- 
  filterTaxa(physeq_18s_Stra, bad_stramenopiles_otus)
physeq_18s_Stra_filt

physeq_18s_Stra_filt_2 <- 
  selectTaxa(physeq_18s_Prot, stramenopiles_otus_2)
physeq_18s_Stra_filt_2

# Merge phyloseq objects
physeq_18s_Stra_filt <- 
  merge_phyloseq(physeq_18s_Stra_filt, physeq_18s_Stra_filt_2)
physeq_18s_Stra_filt


# Rhizaria ---------------------------------------------------------------------
physeq_18s_Rhizaria <-
  physeq_18s_mSeq %>%
  subset_taxa(High_level_taxonomy%in%"Rhizaria" | Kingdom%in%"Rhizaria")
physeq_18s_Rhizaria

write.dna(
  refseq(physeq_18s_Rhizaria),
  format="fasta", 
  colsep="", 
  file="otus_Rhizaria.fasta")

bad_rhizaria_otus <- 
  c("OTU_702", "OTU_2367", "OTU_10588","OTU_2838", "OTU_10405",
    "OTU_2907", "OTU_4145", "OTU_10181", "OTU_7191"
  )


physeq_18s_Rhizaria_filt <- 
  filterTaxa(physeq_18s_Rhizaria, bad_rhizaria_otus)
physeq_18s_Rhizaria_filt

physeq_18s_Rhizaria_filt_2 <- 
  selectTaxa(physeq_18s_Alve, rhizaria_otus)
physeq_18s_Rhizaria_filt_2

physeq_18s_Rhizaria_filt_3 <- 
  selectTaxa(physeq_18s_Ichthyosporia, rhizaria_otus_2)
physeq_18s_Rhizaria_filt_3

physeq_18s_Rhizaria_filt_4 <- 
  selectTaxa(physeq_18s_Metazoa, rhizaria_otus_3)
physeq_18s_Rhizaria_filt_4


physeq_18s_Rhizaria_filt <- 
  merge_phyloseq(physeq_18s_Rhizaria_filt, 
                 physeq_18s_Rhizaria_filt_2,
                 physeq_18s_Rhizaria_filt_3,
                 physeq_18s_Rhizaria_filt_4)
physeq_18s_Rhizaria_filt



# Cryptista --------------------------------------------------------------------
physeq_18s_Cryptista <-
  physeq_18s_mSeq %>%
  subset_taxa(High_level_taxonomy%in%"Cryptista" | Kingdom%in%"Cryptista")
physeq_18s_Cryptista

write.dna(
  refseq(physeq_18s_Cryptista),
  format="fasta", 
  colsep="", 
  file="otus_Cryptista.fasta")

bad_cryptista_otys <-
  c("OTU_8012", "OTU_2916", "OTU_8029", "OTU_2916", "OTU_2152",
    "OTU_4639", "OTU_8747", "OTU_2870", "OTU_4399", "OTU_10731",
    "OTU_5191", "OTU_4326", "OTU_6622", "OTU_4953", "OTU_5417", 
    "OTU_4744", "OTU_6593", "OTU_8635", "OTU_10630"
  )

physeq_18s_Cryptista_filt <- 
  filterTaxa(physeq_18s_Cryptista, bad_cryptista_otys)
physeq_18s_Cryptista_filt


# Euglenozoa -------------------------------------------------------------------
physeq_18s_Euglenozoa <-
  physeq_18s_mSeq %>%
  subset_taxa(High_level_taxonomy%in%"Euglenozoa" | Kingdom%in%"Euglenozoa")
physeq_18s_Euglenozoa

write.dna(
  refseq(physeq_18s_Euglenozoa),
  format="fasta", 
  colsep="", 
  file="otus_Euglenozoa.fasta")


bad_euglenozoa_otus <-
  c("OTU_8611")

physeq_18s_Euglenozoa_filt <- 
  filterTaxa(physeq_18s_Euglenozoa, bad_euglenozoa_otus)
physeq_18s_Euglenozoa_filt

# Ichthyosporia ----------------------------------------------------------------

physeq_18s_Ichthyosporia <-
  physeq_18s_mSeq %>%
  subset_taxa(High_level_taxonomy%in%"Ichthyosporia" | Kingdom%in%"Ichthyosporia")
physeq_18s_Ichthyosporia

write.dna(
  refseq(physeq_18s_Ichthyosporia),
  format="fasta", 
  colsep="", 
  file="otus_Ichthyosporia.fasta")

bad_Ichthyosporia_otus <-
  c("OTU_3399", "OTU_6439", "OTU_943", "OTU_9172", "OTU_1039",
    "OTU_6929",  "OTU_10758", "OTU_6815", "OTU_2750",  "OTU_8908", 
    "OTU_3773", "OTU_7653",  "OTU_10827", "OTU_8985", "OTU_5150",
    "OTU_9324"
  )


metazoa_otus_2 <-
  c("OTU_162", "OTU_8212","OTU_11004", "OTU_5266", "OTU_2947",
    "OTU_2527","OTU_3803","OTU_8128","OTU_10356", "OTU_10337"
  )

rhizaria_otus_2 <- 
  c("OTU_270", "OTU_594", "OTU_3419", "OTU_1229", 
    "OTU_2741", "OTU_8548", "OTU_8180", "OTU_5387", "OTU_609",
    "404700", "OTU_7413", "OTU_4481", "OTU_1170", "OTU_6304", 
    "OTU_11059", "OTU_5109", "OTU_2751", "OTU_7254", "OTU_3842", 
    "OTU_3518", "OTU_9212"
  )


physeq_18s_Ichthyosporia_filt <- 
  filterTaxa(physeq_18s_Ichthyosporia, c(
    bad_Ichthyosporia_otus,
    metazoa_otus_2,
    rhizaria_otus_2))
physeq_18s_Ichthyosporia_filt

physeq_18s_Ichthyosporia_filt_2 <- 
  selectTaxa(physeq_18s_Choanoflagellozoa, Ichthyosporea_otus)
physeq_18s_Ichthyosporia_filt_2

physeq_18s_Ichthyosporia_filt <- 
  merge_phyloseq(physeq_18s_Ichthyosporia_filt, 
                 physeq_18s_Ichthyosporia_filt_2)
physeq_18s_Ichthyosporia_filt


# Metazoa ----------------------------------------------------------------------

physeq_18s_Metazoa <-
  physeq_18s_mSeq %>%
  subset_taxa(High_level_taxonomy%in%"Metazoa" | Kingdom%in%"Metazoa")
physeq_18s_Metazoa

write.dna(
  refseq(physeq_18s_Metazoa),
  format="fasta", 
  colsep="", 
  file="otus_Metazoa.fasta")


bad_metazoa_otus <-
  c("OTU_4096", "OTU_4439", "OTU_52", "OTU_99", "OTU_458", 
    "OTU_7592","OTU_3483", "OTU_3483", "OTU_102", "OTU_158", 
    "OTU_6076", "OTU_5809", "OTU_8548", "OTU_3979", "OTU_1136",
    "OTU_2367", "OTU_4216","OTU_8597", "OTU_1139", "OTU_6927", 
    "OTU_11004", "OTU_3863", "OTU_10653", "OTU_8871","OTU_9172", 
    "OTU_4408", "OTU_10706", "OTU_3418", "OTU_6370", "OTU_2779",
    "OTU_2645", "OTU_10491", "OTU_1033", "OTU_8175", "OTU_9549", 
    "OTU_3108", "OTU_7793", "OTU_3083", "OTU_10541", "OTU_3702",
    "OTU_6752", "OTU_5928", "OTU_5675", "OTU_3025", "OTU_10356",
    "OTU_6584", "OTU_5834", "OTU_9964", "OTU_7280", "OTU_10789",
    "OTU_8739", "OTU_7191", "OTU_8318"
  )


rhizaria_otus_3 <-
  c("OTU_575", "OTU_1192", "OTU_270", "OTU_482", "OTU_330", 
    "OTU_1003", "OTU_250", "OTU_4624", "OTU_186",  "OTU_2133",
    "OTU_146", "OTU_5072", "OTU_151", "OTU_635", "OTU_1229", 
    "OTU_562", "OTU_1729", "OTU_2164", "OTU_410", "OTU_574", 
    "OTU_2192", "OTU_859", "OTU_1940", "OTU_297", "OTU_1291",
    "OTU_2345", "OTU_998", "OTU_729", "OTU_7843", "OTU_5387", 
    "OTU_1732", "OTU_973", "OTU_1045", "OTU_1491", "OTU_8328", 
    "OTU_1883", "OTU_4481", "OTU_2938", "OTU_1170", "OTU_2142",
    "OTU_5542", "OTU_2331", "OTU_2959", "OTU_10718", "OTU_11059", 
    "OTU_2727", "OTU_1410", "OTU_3035", "OTU_2683", "OTU_1689", 
    "OTU_9588", "OTU_6932", "OTU_3794", "OTU_1381", "OTU_5210", 
    "OTU_8692", "OTU_7254", "OTU_2394", "OTU_2306", "OTU_1441", 
    "OTU_2601", "OTU_5486", "OTU_9742", "OTU_7742", "OTU_3025", 
    "OTU_5133", "OTU_8065", "OTU_6809", "OTU_5565", "OTU_5840", 
    "OTU_6425", "OTU_5468", "OTU_8638", "OTU_9041", "OTU_9212"
  )

alveolata_otus_1 <-
  c("OTU_2224", "OTU_4938", "OTU_4494", "OTU_6503", "OTU_9261",
    "OTU_8023", "OTU_2275", "OTU_7084"
  )

physeq_18s_Metazoa_filt <- 
  filterTaxa(physeq_18s_Metazoa, c(
    bad_metazoa_otus,
    rhizaria_otus_3,
    alveolata_otus_1))
physeq_18s_Metazoa_filt


# Choanoflagellozoa ------------------------------------------------------------
physeq_18s_Choanoflagellozoa <-
  physeq_18s_mSeq %>%
  subset_taxa(High_level_taxonomy%in%"Choanoflagellozoa" | Kingdom%in%"Choanoflagellozoa")
physeq_18s_Choanoflagellozoa

write.dna(
  refseq(physeq_18s_Choanoflagellozoa),
  format="fasta", 
  colsep="", 
  file="otus_Choanoflagellozoa.fasta")

bad_Choanoflagellozoa_otus <-
  c("OTU_346","OTU_1025", "OTU_5459", "OTU_1561","OTU_1910",
    "OTU_278"
  )

Ichthyosporea_otus <-
  c("OTU_346", "OTU_1487", "OTU_5301","OTU_9007", "OTU_4998",
    "OTU_3714", "OTU_8449", "OTU_5347", "OTU_10642", "OTU_8511"
  )


physeq_18s_Choanoflagellozoa_filt <- 
  filterTaxa(physeq_18s_Choanoflagellozoa, c(
    Ichthyosporea_otus,
    bad_Choanoflagellozoa_otus))
physeq_18s_Choanoflagellozoa_filt

physeq_18s_Choanoflagellozoa_filt_2 <- 
  selectTaxa(physeq_18s_Filasteriae, Choanoflagellata_otus)
physeq_18s_Choanoflagellozoa_filt_2

physeq_18s_Choanoflagellozoa_filt <- 
  merge_phyloseq(physeq_18s_Choanoflagellozoa_filt, 
                 physeq_18s_Choanoflagellozoa_filt_2)
physeq_18s_Choanoflagellozoa_filt


# Filasteriae ------------------------------------------------------------------
physeq_18s_Filasteriae <-
  physeq_18s_mSeq %>%
  subset_taxa(High_level_taxonomy%in%"Filasteriae" | Kingdom%in%"Filasteriae")
physeq_18s_Filasteriae

write.dna(
  refseq(physeq_18s_Filasteriae),
  format="fasta", 
  colsep="", 
  file="otus_Filasteriae.fasta")

# I can porbably put this together with the Choanoflagellata and Ichthyosporea

bad_filasteriae_otus <-
  c("OTU_7491", "OTU_7144", "OTU_5834", "OTU_7610"
  )

Choanoflagellata_otus <-
  c("OTU_4439", "OTU_2976", "OTU_6373", "OTU_2669", "OTU_5621"
  )


physeq_18s_Filasteriae_filt <- 
  filterTaxa(physeq_18s_Filasteriae, c(
    bad_filasteriae_otus,
    Choanoflagellata_otus))
physeq_18s_Filasteriae_filt

# Filtered datasets by organism group ------------------------------------------

physeq_18s_Alve_filt
physeq_18s_Prot_filt
physeq_18s_Stra_filt
physeq_18s_Rhizaria_filt
physeq_18s_Cryptista_filt
physeq_18s_Euglenozoa_filt
physeq_18s_Ichthyosporia_filt
physeq_18s_Metazoa_filt
physeq_18s_Choanoflagellozoa_filt
physeq_18s_Filasteriae_filt

# Writing Alveolata - Rhizaria tablee
write.csv(physeq_18s_Alve_filt@otu_table, "table_Alveolata.csv")
write.csv(physeq_18s_Alve_filt@tax_table, "tax_Alveolata.csv")
write.csv(physeq_18s_Alve_filt@sam_data, "meta_Alveolata.csv")


addGroup <- function(physeq, Var){
  if (length(colnames(physeq@tax_table))> 14){
    stop("Group is already present in sam_data") 
  }else{
    physeq@tax_table <-
      tax_table(
        as(
          data.frame(
            as.data.frame(
              as(physeq@tax_table, "matrix")), 
            Group = rep(Var, 
                        length(taxa_names(physeq))
            )), "matrix")
      )
    return(physeq)
  }
}

phy_alveolata <- addGroup(physeq_18s_Alve_filt, "Alveolata")
phy_protista <- addGroup(physeq_18s_Prot_filt, "Protista")
phy_stramenopiles <- addGroup(physeq_18s_Stra_filt, "Stramenopiles")
phy_rhizaria <- addGroup(physeq_18s_Rhizaria_filt, "Rhizaria")
phy_cryptista <- addGroup(physeq_18s_Cryptista_filt, "Cryptista")
phy_euglenozoa <- addGroup(physeq_18s_Euglenozoa_filt, "Euglenozoa")
phy_ichthyosporia <- addGroup(physeq_18s_Ichthyosporia_filt, "Ichthyosporia")
phy_metazoa <- addGroup(physeq_18s_Metazoa_filt, "Metazoa")
phy_choanoflagellozoa <- addGroup(physeq_18s_Choanoflagellozoa_filt, "Choanoflagellozoa")
phy_filasteriae <- addGroup(physeq_18s_Filasteriae_filt, "Filasteriae")

physeq_18s_all_groups <-
  merge_phyloseq(phy_alveolata, 
                 phy_protista,
                 phy_stramenopiles,
                 phy_rhizaria,
                 phy_cryptista,
                 phy_euglenozoa,
                 phy_ichthyosporia,
                 phy_metazoa,
                 phy_choanoflagellozoa,
                 phy_filasteriae
  )


physeq_18s_all_groups
physeq_18s_all_groups@otu_table

nrow(phy_protista@otu_table)+
  nrow(phy_alveolata@otu_table)+
  nrow(phy_protista@otu_table)+
  nrow(phy_stramenopiles@otu_table)+
  nrow(phy_rhizaria@otu_table)+
  nrow(phy_cryptista@otu_table)+
  nrow(phy_euglenozoa@otu_table)+
  nrow(phy_ichthyosporia@otu_table)+
  nrow(phy_metazoa@otu_table)+
  nrow(phy_choanoflagellozoa@otu_table)+
  nrow(phy_filasteriae@otu_table)


# Writing Alveolata - Rhizaria tablee
write.csv(physeq_18s_all_groups@otu_table, "table_18s_by_groups.csv")
write.csv(physeq_18s_all_groups@tax_table, "tax_18s_by_groups.csv")
write.csv(physeq_18s_all_groups@sam_data, "meta_18s_by_groups.csv")

physeq_18s_all_groups@sam_data

as.data.frame(as.matrix(tax_table(physeq_18s_all_groups)))$High_level_taxonomy

# removing samples with less than 2 reads
physeq_18s_all_groups_filt <-
  physeq_18s_all_groups

otu_table(physeq_18s_all_groups_filt) <- 
  prune_samples(sample_sums(physeq_18s_all_groups) > 5, physeq_18s_all_groups) 

# **************************************************************************----
# *********NEW MODIFIED CODE ENDS ****************------------------------------
# **************************************************************************----
