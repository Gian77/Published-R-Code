# ****** DEEP CORE ALPHA AND BETA ******************************************************----
# Manuscript: Stochastic and deterministic processes shape bioenergy crop microbiomes 
#             along a vertical soil niche
# Authors:      Benucci GMN, ... Bonito G.
# Affiliation:  Michigan State University - GLBRC
# Journal:      environmental Microbiolgy Journal
# Date:         Oct 23, 2022
# *************************************************************************************-----

# WORKING ENVIRONMENT SETUP ----------------------------------------------------------------
options(scipen = 9999)
options(max.print=100000000)

# loading required packages ----------------------------------------------------------------
library(styler)
library(evoPalette)

library(phyloseq)
library(Biostrings)
library(ape)

library(tidyr)
library(tidyverse)
library(reshape2)
library(tidygraph)
library(ggpubr)

library(vegan)
library(decontam)
library(mctoolsr)
library(RVAideMemoire)
library(multcompView)
library(grid)
library(gridExtra)

library(magrittr)
library(SpiecEasi)
library(igraph)
library(ggraph)

library(Boruta)
library(car)
library(broom)
library(lme4)
library(MASS)

# Generating color palettes for plotting ---------------------------------------------------
paletteCB4 = c("#E69F00", "#D55E00", "#CBD588", "#0072B2")

pie(rep(1, length(paletteCB4)), labels = sprintf("%d (%s)",
    seq_along(paletteCB4),paletteCB4), col = paletteCB4)

paletteCB3 = c("#D0CED1", "#A27F5E", "#000000")
paletteCB3 = c("#57B4AE", "#CB54D6", "#21191A")
paletteCB3 = c("#BBC9E5","#5FBB68","#EB1E2C")
paletteCB3 = c("#CB54D6","gray60", "#09104d")

pie(rep(1, length(paletteCB3)), labels = sprintf("%d (%s)",
    seq_along(paletteCB3),paletteCB3), col = paletteCB3)


# IMPORTING DATASETS -----------------------------------------------------------------------
# A) uparse ITS R1 -------------------------------------------------------------------------
# Importing taxonomies at 0.6 confidence ---------------------------------------------------
taxonomy_ITS06 <-
  read.delim(
    "taxonomy_assignments_euk_06/consensus_taxonomy.txt",
    header = TRUE,
    row.names = 1,
    sep = "\t")

head(taxonomy_ITS06)
taxonomy_ITS06[1:100, ]

# Importing taxonomies at 0.8 confidence ---------------------------------------------------
taxonomy_ITS08 <-
  read.delim(
    "taxonomy_assignments_euk_08/consensus_taxonomy.txt",
    header = TRUE,
    row.names = 1, 
    sep = "\t")

head(taxonomy_ITS08)
taxonomy_ITS08[1:100, ]

# check for identical ordering
identical(rownames(taxonomy_ITS06), rownames(taxonomy_ITS08))
# sample_order <- match(rownames(taxonomy_ITS08), rownames(taxonomy_ITS06))
# taxonomy_ITS06 <- taxonomy_ITS06[,sample_order]

taxonomy_ITS08$Kingdom_06 <- taxonomy_ITS06$Kingdom
head(taxonomy_ITS08)
taxonomy_ITS08[1:100, ]
dim(taxonomy_ITS08)

levels(taxonomy_ITS08$Kingdom)

# how many unclassified OTUs in the two taxonomies?
nrow(taxonomy_ITS06[taxonomy_ITS06$Kingdom!="Fungi",])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom!="Fungi",])

# Non-target taxa
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom=="Ichthyosporia",])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom=="Metazoa",])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom=="Protista",])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom=="Rhizaria",])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom=="Stramenopila",])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom=="Viridiplantae",])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom=="",])

nrow(taxonomy_ITS08[taxonomy_ITS06$Kingdom=="Fungi",])

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
dim(taxonomy_ITS08_filt)
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

# Reduce isolate name length, to escape the special charater use \\
taxonomy_ITS08_filt$Isolate <- gsub("\\|.*","",taxonomy_ITS08_filt$Isolate)
# cut before the last |
#gsub("(.*)\\|.*","\\1",taxonomy_ITS08_filt$Isolate) 
head(taxonomy_ITS08_filt)

# Removing Lab contaminants in the ITS data ------------------------------------------------
str(taxonomy_ITS08_filt)

taxonomy_ITS08_filt[taxonomy_ITS08_filt$Order=="Mortierellales",]

subset(taxonomy_ITS08_filt, taxonomy_ITS08_filt$Isolate_percent_id<=99) -> taxonomy_ITS08_clean
dim(taxonomy_ITS08_clean)

taxonomy_ITS08_clean[taxonomy_ITS08_clean$Family=="Mortierellaceae",]

# importing other tables -------------------------------------------------------------------
otus_ITS_uparse_R1 <-
  read.delim(
    "otu_table_ITS_UPARSE_R1_new.txt",
    header = TRUE,
    row.names = 1,)

metadata_ITS_uparse_R1 <-
  read.delim(
    "map_ITS.txt",
    row.names = 1,
    header = TRUE,
    sep = "\t")

otus_seq_ITS_uparse_R1 <-
  readDNAStringSet(
    "otus_ITS.fasta",
    format = "fasta",
    seek.first.rec = TRUE,
    use.names = TRUE)

physeq_ITS_uparse <- phyloseq(otu_table(otus_ITS_uparse_R1, taxa_are_rows = TRUE),
                              sample_data(metadata_ITS_uparse_R1),
                              tax_table(as.matrix(taxonomy_ITS08_clean)),
                              otus_seq_ITS_uparse_R1) 

physeq_ITS_uparse
str(physeq_ITS_uparse)
head(sample_data(physeq_ITS_uparse))
tax_table(physeq_ITS_uparse)[tax_table(physeq_ITS_uparse)==""]<- NA
head(tax_table(physeq_ITS_uparse))

saveRDS(object = physeq_ITS_uparse, file = "physeq_ITS_uparse.RDS")

# checking the phyloseq object
sort(unique(as.data.frame(tax_table(physeq_ITS_uparse))$Kingdom)) # everything is Fungi
nrow(as.data.frame(tax_table(physeq_ITS_uparse))[as.data.frame(tax_table(physeq_ITS_uparse))$Kingdom!="Fungi",])


# A) uparse 16S paired ---------------------------------------------------------------------
# Importing CONSTAX SILVA taxonomy ---------------------------------------------------------
silva_taxonomy_16s <-
  read.delim("taxonomy_assignments_prok_08/consensus_taxonomy.txt",
             header = TRUE,
             row.names = 1)

head(silva_taxonomy_16s)

# cleaning taxonomy labels -----------------------------------------------------------------
colnames(silva_taxonomy_16s) <-
  c("Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
    "OTU_ID")

silva_taxonomy_16s$OTU_ID <- rownames(silva_taxonomy_16s)

silva_taxonomy_16s[, "Kingdom"] <- as.factor(gsub("_1", "", silva_taxonomy_16s[, "Kingdom"]))
silva_taxonomy_16s[, "Phylum"] <- as.factor(gsub("_1", "", silva_taxonomy_16s[, "Phylum"]))
silva_taxonomy_16s[, "Class"] <- as.factor(gsub("_1", "", silva_taxonomy_16s[, "Class"]))
silva_taxonomy_16s[, "Order"] <- as.factor(gsub("_1", "",silva_taxonomy_16s[, "Order"]))
silva_taxonomy_16s[, "Family"] <- as.factor(gsub("_1", "", silva_taxonomy_16s[, "Family"]))
silva_taxonomy_16s[, "Genus"] <- as.factor(gsub("_1", "", silva_taxonomy_16s[, "Genus"]))
silva_taxonomy_16s[, "Species"] <- as.factor(gsub("_1", "", silva_taxonomy_16s[, "Species"]))

head(silva_taxonomy_16s)
str(silva_taxonomy_16s)

# silva_taxonomy_16s[1:50,]
# endsWith("aceae", as.character(silva_taxonomy_16s$Family))

any(silva_taxonomy_16s$Kingdom == "Chloroplast")
any(silva_taxonomy_16s$Kingdom == "Mitochondria")
any(silva_taxonomy_16s$Phylum == "Chloroplast") 
any(silva_taxonomy_16s$Phylum == "Mitochondria")
any(silva_taxonomy_16s$Class == "Chloroplast")
any(silva_taxonomy_16s$Class == "Mitochondria")
any(silva_taxonomy_16s$Order == "Chloroplast") #TRUE
any(silva_taxonomy_16s$Order == "Mitochondria")
any(silva_taxonomy_16s$Family == "Chloroplast") 
any(silva_taxonomy_16s$Family == "Mitochondria")#TRUE
any(silva_taxonomy_16s$Genus == "Chloroplast")
any(silva_taxonomy_16s$Genus == "Mitochondria")

silva_taxonomy_16s %>% dplyr::filter(silva_taxonomy_16s$Order == "Chloroplast")
silva_taxonomy_16s %>% dplyr::filter(silva_taxonomy_16s$Family == "Mitochondria")
silva_taxonomy_16s <- subset(silva_taxonomy_16s, Family != "Mitochondria" | Order!= "Chloroplast")
dim(silva_taxonomy_16s)

# Check for unclassified OTUs and remove them
any(silva_taxonomy_16s$Kingdom == "")
nrow(silva_taxonomy_16s[silva_taxonomy_16s$Kingdom == "", ])

silva_taxonomy_16s %>% dplyr::filter(silva_taxonomy_16s$Kingdom == "")
silva_taxonomy_16s <- subset(silva_taxonomy_16s, Kingdom != "")

# importing other tables
otus_16s_uparse <-
  read.delim("otu_table_16S_UPARSE_R1_new.txt",
             row.names = 1,
             header = TRUE)

metadata_16s_uparse <-
  read.delim(
    "map_16S.txt",
    row.names = 1,
    header = TRUE,
    sep = "\t"
  )

otus_seq_16s_uparse <-
  readDNAStringSet(
    "otus_16s.fasta",
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

saveRDS(object = physeq_16s_uparse, file = "physeq_16s_uparse.RDS")

# ********************************************----------------------------------------------
# CREATING NEW PHYLOSEQ OBJETCS ------------------------------------------------------------
physeq_ITS_uparse -> physeq_fungi_uparse
physeq_fungi_uparse

physeq_16s_uparse -> physeq_prok_uparse
physeq_prok_uparse

# Change OTU names -------------------------------------------------------------------------
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

ChangeName(physeq_ITS_uparse, "F") -> physeq_fungi_uparse
ChangeName(physeq_16s_uparse, "P") -> physeq_prok_uparse

# REFORMAT TAXONOMY ------------------------------------------------------------------------
blank2na = function(x, na.strings=c('','.','NA','na','N/A','n/a','NaN','nan')) {
  if (is.factor(x)) {
    lab = attr(x, 'label', exact = T)
    labs1 <- attr(x, 'labels', exact = T)
    labs2 <- attr(x, 'value.labels', exact = T)
    # trimws will convert factor to character
    x = trimws(x,'both')
    if (! is.null(lab)) lab = trimws(lab,'both')
    if (! is.null(labs1)) labs1 = trimws(labs1,'both')
    if (! is.null(labs2)) labs2 = trimws(labs2,'both')
    if (!is.null(na.strings)) {
      # convert to NA
      x[x %in% na.strings] = NA
      # also remember to remove na.strings from value labels 
      labs1 = labs1[! labs1 %in% na.strings]
      labs2 = labs2[! labs2 %in% na.strings]
    }
    # the levels will be reset here
    x = factor(x)
    if (! is.null(lab)) attr(x, 'label') <- lab
    if (! is.null(labs1)) attr(x, 'labels') <- labs1
    if (! is.null(labs2)) attr(x, 'value.labels') <- labs2
  } else if (is.character(x)) {
    lab = attr(x, 'label', exact = T)
    labs1 <- attr(x, 'labels', exact = T)
    labs2 <- attr(x, 'value.labels', exact = T)
    # trimws will convert factor to character
    x = trimws(x,'both')
    if (! is.null(lab)) lab = trimws(lab,'both')
    if (! is.null(labs1)) labs1 = trimws(labs1,'both')
    if (! is.null(labs2)) labs2 = trimws(labs2,'both')
    if (!is.null(na.strings)) {
      # convert to NA
      x[x %in% na.strings] = NA
      # also remember to remove na.strings from value labels 
      labs1 = labs1[! labs1 %in% na.strings]
      labs2 = labs2[! labs2 %in% na.strings]
    }
    if (! is.null(lab)) attr(x, 'label') <- lab
    if (! is.null(labs1)) attr(x, 'labels') <- labs1
    if (! is.null(labs2)) attr(x, 'value.labels') <- labs2
  } else {
    x = x
  }
  return(x)
}

# In the tax_table add a column naming the highest resolution taxonomy 
# achieved for each OTU, remove _ and add sp. to genera
ReformatTaxonomy <- function(dataframe){
  taxa_table <- as.data.frame(as.matrix(tax_table(dataframe)))
  # remember to do run this function only once on your dataframe
  taxa_table$Genus <- as.character(taxa_table$Genus)
  taxa_table[taxa_table=="Unclassified"]<- NA
  taxa_table[taxa_table=="Unidentified"]<- NA
  taxa_table[taxa_table==""]<- NA
  taxa_table[which(is.na(taxa_table$Genus) == FALSE), ]$Genus <- paste(
    taxa_table$Genus[is.na(taxa_table$Genus)==FALSE], "sp.", sep = " ")
  taxa_table$OTU_ID <- rownames(taxa_table)
  taxa_table <- taxa_table[c(8,1,2,3,4,5,6,7)] 
  taxa_table[] = lapply(taxa_table, blank2na, na.strings=c('','NA','na','N/A','n/a','NaN','nan'))
  lastValue <- function(x) tail(x[!is.na(x)], 1)
  last_taxons<- apply(taxa_table[,1:8], 1, lastValue)
  taxa_table$BestMatch <- last_taxons
  taxa_table[, "BestMatch"] <- gsub("_", " ", taxa_table[, "BestMatch"])
  taxa_table$Taxonomy <- paste(taxa_table$OTU_ID, taxa_table$BestMatch, sep="-")
  taxa_table[, "Genus"] <- gsub(" sp.", "", taxa_table[, "Genus"])
  tax_table(dataframe) <- tax_table(as.matrix(taxa_table))
  return(dataframe)
}


ReformatTaxonomy(physeq_fungi_uparse) -> physeq_fungi_uparse
head(tax_table(physeq_fungi_uparse))

ReformatTaxonomy(physeq_prok_uparse) -> physeq_prok_uparse
head(tax_table(physeq_prok_uparse))

# ********************************************----------------------------------------------
# INSPECTING LIBRARY SIZES -----------------------------------------------------------------
# Before removing controls
df_fungi_uparse <-
  as.data.frame(as.matrix(sample_data(physeq_fungi)))
df_fungi_uparse$LibrarySize <- sample_sums(physeq_fungi)
df_fungi_uparse <-
  df_fungi_uparse[order(df_fungi_uparse$LibrarySize), ]
df_fungi_uparse$Index <-
  seq(nrow(df_fungi_uparse)) # sample numbering
head(df_fungi_uparse)

df_prok_uparse <-
  as.data.frame(as.matrix(sample_data(physeq_prok)))
df_prok_uparse$LibrarySize <- sample_sums(physeq_prok)
df_prok_uparse <-
  df_prok_uparse[order(df_prok_uparse$LibrarySize),]
df_prok_uparse$Index <- seq(nrow(df_prok_uparse)) # sample numbering
head(df_prok_uparse)

# Plotting sample read depth --------------------------------------------------------------
PlotDist <- function(df, width){
  ggplot(df, aes(x=LibrarySize)) +
    labs(title="Distribution of Sample Libraries", x="Read number", y="Sample number") + 
    facet_grid(Depth_Core ~ .) +
    geom_histogram(binwidth=width, colour="grey80", fill="grey80") +
    geom_vline(aes(xintercept=mean(LibrarySize, na.rm=T)), color="red", linetype="dashed", size=0.8) +   # Ignore NA values for mean
    #geom_vline(aes(xintercept=median(ReadNO, na.rm=T)), color="blue", linetype="dashed", size=0.8) + 
    theme_classic() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) + 
    theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
    theme(legend.position="bottom") -> plot_dist
  return(plot_dist)  
}

# Calculating binwidths to make them all the same size
max(df_fungi_uparse[df_fungi_uparse$Origin=="Soil",]$LibrarySize) #27500
1000*27500/145492 #189.0138
max(df_fungi_uparse[df_fungi_uparse$Origin=="Root",]$LibrarySize) #44281
1000*44281/145492 #304.3535

max(df_prok_uparse[df_prok_uparse$Origin=="Soil",]$LibrarySize) #145492
1000*145492/145492 #1000
max(df_prok_uparse[df_prok_uparse$Origin=="Root",]$LibrarySize) #108502
1000*108502/145492 #745.7592

PlotDist(df_fungi_uparse[df_fungi_uparse$Origin=="Soil",], 189.0138) +
  labs(title="Fungi - soil\nDistribution of sample libraries")

# FILTERING OUT CONTAMINANTS ---------------------------------------------------------------
head(sample_data(physeq_fungi_uparse))
subset_samples(physeq_fungi_uparse, sample_sums(physeq_fungi_uparse) >0) -> physeq_fungi

# detecting contaminants Fungi
sample_data(physeq_fungi)$is.neg <-
  sample_data(physeq_fungi)$Origin == "Control"

contam_prev_fungi <-
  isContaminant(
    physeq_fungi,
    method = "prevalence",
    neg = "is.neg",
    threshold = 0.5)

table(contam_prev_fungi$contaminant)

head(sample_data(physeq_prok_uparse))
subset_samples(physeq_prok_uparse, sample_sums(physeq_prok_uparse) >0) -> physeq_prok

# detecting contaminants Prokaryotes
sample_data(physeq_prok)$is.neg <-
  sample_data(physeq_prok)$Origin == "Control"

contam_prev_prokaryote <-
  isContaminant(
    physeq_prok,
    method = "prevalence",
    neg = "is.neg",
    threshold = 0.5)

table(contam_prev_prokaryote$contaminant)

# plotting contaminant OTUs ----------------------------------------------------------------
PlotContam <- function(df, contam){
  # Make phyloseq object of presence-absence in negative controls and true samples
  physeq_pa <- transform_sample_counts(df, function(abund) 1*(abund>0))
  physeq_pa_neg <- subset_samples(physeq_pa, Origin%in%c("Control"))
  physeq_pa_pos <- subset_samples(physeq_pa, Origin%in%c("Soil", "Root"))
  # Make data.frame of prevalence in positive and negative samples
  df_contam <- data.frame(pa.pos=taxa_sums(physeq_pa_pos), 
                          pa.neg=taxa_sums(physeq_pa_neg),
                          contaminant=contam$contaminant, 
                          Pvalue=contam$p)
  # plotting 
  ggplot(data=df_contam, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
    geom_point(size=0.8, alpha=0.7) +
    labs(x="Prevalence in negative controls", y="Prevalence in true samples") +
    scale_colour_manual("Contaminant OTUs", values = c("grey", "red")) +
    theme_classic() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) + 
    theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
    theme(legend.position="bottom") -> plot_cont
  return(plot_cont)
}

PlotContam(physeq_fungi, contam_prev_fungi)
PlotContam(physeq_prok, contam_prev_prokaryote) 

# *** FIGURE S1 - dcontam ------------------------------------------------------------------
ggarrange(
  ggarrange(PlotDist(df_fungi_uparse[df_fungi_uparse$Origin=="Soil",], 189.0138) +
              labs(title="Fungi - soil\nDistribution of sample libraries"),
            PlotDist(df_prok_uparse[df_prok_uparse$Origin=="Soil",], 1000) +
              labs(title="Prokaryote - soil\nDistribution of sample libraries"),
            PlotDist(df_fungi_uparse[df_fungi_uparse$Origin=="Root",], 304.3535) +
              labs(title="Fungi - root\nDistribution of sample libraries"),
            PlotDist(df_prok_uparse[df_prok_uparse$Origin=="Root",], 745.7592) +
              labs(title="Prokaryote - root\nDistribution of sample libraries"),
            labels = c("A","B","",""),
            widths = c(1,1,1,1),
            align = "hv" ,
            ncol = 2, 
            nrow = 2, 
            legend = c("none")),
  ggarrange(PlotContam(physeq_fungi, contam_prev_fungi) + labs(title="Fungi\nContaminants"),
            PlotContam(physeq_prok, contam_prev_prokaryote) + labs(title="Prokaryotes\nContaminants"),
            widths = c(1,1),
            labels = c("C",""),
            align = "hv" ,
            ncol = 1, 
            nrow = 2,
            common.legend = TRUE,
            legend = c("bottom")),
  widths =  c(1,0.5),
  ncol = 2, 
  nrow = 1) -> Fig_S1

Fig_S1

# REMOVING CONTAMINANTS --------------------------------------------------------------------
# function to remove dab taxa
remove_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

physeq_fungi_clean <-
  remove_taxa(physeq_fungi, rownames(subset(
    contam_prev_fungi, contaminant %in% c("TRUE"))))

physeq_prok_clean <-
  remove_taxa(physeq_prok, rownames(subset(
    contam_prev_prokaryote, contaminant %in% c("TRUE"))))


write.csv(physeq_fungi_clean@sam_data, "samples_used_ITS.csv")
write.csv(physeq_prok_clean@sam_data, "samples_used_16S.csv")

intersect(rownames(physeq_fungi_clean@sam_data),rownames(physeq_prok_clean@sam_data)) -> shared_samples

physeq_fungi_clean@sam_data$rownames <- rownames(physeq_fungi_clean@sam_data)
physeq_prok_clean@sam_data$rownames <- rownames(physeq_prok_clean@sam_data)

physeq_fungi_used <-
  subset_samples(physeq_fungi_clean, rownames%in%shared_samples)
physeq_fungi_used@otu_table <- 
  physeq_fungi_used@otu_table[which(rowSums(physeq_fungi_used@otu_table) > 0),] 
physeq_fungi_used

physeq_prok_used <-
  subset_samples(physeq_prok_clean, rownames%in%shared_samples)
physeq_prok_used@otu_table <- 
  physeq_prok_used@otu_table[which(rowSums(physeq_prok_used@otu_table) > 0),] 
physeq_prok_used


write.csv(physeq_fungi_used@sam_data, "samples_used_ITS_new.csv")
write.csv(physeq_prok_used@sam_data, "samples_used_16S_new.csv")

# REMOVING ALL CONTROLS ------------------------------------------------------------------- 
# fungi
physeq_fungi_clean <-
  subset_samples(physeq_fungi_clean, Origin %in% c("Root", "Soil"))

# Checking for residual Mock OTUs
head(taxonomy_ITS08_clean)
c(rownames(taxonomy_ITS08_clean[taxonomy_ITS08_clean$Isolate == "mock_1",]),
  rownames(taxonomy_ITS08_clean[taxonomy_ITS08_clean$Isolate == "mock_5",])) -> mock_otu

physeq_fungi_clean <-remove_taxa(physeq_fungi_clean, mock_otu)

otu_table(physeq_fungi_clean) <-
  otu_table(physeq_fungi_clean)[which(rowSums(otu_table(physeq_fungi_clean)) > 0), ]
physeq_fungi_clean

#looking for tag-switching using mock OTUs
otu_table(physeq_ITS_uparse)[mock_otu,]
sample_data(physeq_ITS_uparse)[rownames(sample_data(physeq_ITS_uparse))=="SAM1",]
sample_data(physeq_ITS_uparse)[rownames(sample_data(physeq_ITS_uparse))=="SAM28",]

# Analysing mock samples
sample_data(physeq_ITS_uparse)

#prokaryotes
physeq_prok_clean <-
  subset_samples(physeq_prok_clean, Origin %in% c("Root", "Soil"))

otu_table(physeq_prok_clean) <-
  otu_table(physeq_prok_clean)[which(rowSums(otu_table(physeq_prok_clean)) > 0),]
physeq_prok_clean


# FILTERING DATASTES following quality  ----------------------------------------------------
# Oliver et al. 2015, PCR errors and tag switching
# Barberan et al. 2012, removing OTUs that appear in less than x samples
# Lindhal et al. 2013, tag switching - that's a good  one!

otu_table(physeq_fungi_clean) <-
  otu_table(physeq_fungi_clean)[which(rowSums(otu_table(physeq_fungi_clean)) >= 10), ] ### PCR errors
physeq_fungi_clean

sum(taxa_sums(physeq_fungi_clean))

otu_table(physeq_prok_clean) <-
  otu_table(physeq_prok_clean)[which(rowSums(otu_table(physeq_prok_clean)) >= 10), ]
physeq_prok_clean

sum(taxa_sums(physeq_prok_clean))

# RAREFACTION CURVES -----------------------------------------------------------------------
rarecurve(t(as.data.frame(otu_table(physeq_fungi_clean))), step = 1000, label = FALSE) -> rare_fungi
rarecurve(t(as.data.frame(otu_table(physeq_prok_clean))), step = 1000, label = FALSE) -> rare_prok

# convert rarecurve() output to data.frame, bind_rows doesn't work because items are 
# different lengths, also need to extract sample sizes from attribute. Allocate result
# into a dataframe.
ExtrRare <- function(x){
  require(purrr)
  require(vegan)
  nsamples <- map_int(x, length)
  total_samples <- sum(nsamples)
  if(!is.null(names(x))){
    sites <- names(x)
  } else {
    sites <- as.character(1:length(nsamples))
  }
  result <- tibble(Site = rep("", total_samples),
                   Sample_size = rep(0, total_samples),
                   Species = rep(0, total_samples))
  start <- 1
  for (i in 1:length(nsamples)){
    result[start:(start + nsamples[i]-1), "Site"] <- sites[i]
    result[start:(start + nsamples[i]-1), "Sample_size"] <- attr(x[[i]], "Subsample")
    result[start:(start + nsamples[i]-1), "Species"] <- x[[i]]
    start <- start + nsamples[i]
  }
  return(result)
}


ExtrRare(rare_fungi)
ExtrRare(rare_prok)

# Then plotting Rarecurve in ggplot2
PlotRareCurve <- function(rare_obj, physeq){
  names(rare_obj) <- rownames(sample_data(physeq))
  ExtrRare(rare_obj) -> rare_df
  sample_data(physeq)$Site <- rownames(sample_data(physeq))
  left_join(rare_df, as.data.frame(sample_data(physeq)), by = "Site") -> rare_joint
  rare_joint$Depth_Core <- paste(rare_joint$Depth_Core, "cm", sep = "")
  rare_joint$Depth_Core <- factor(rare_joint$Depth_Core,levels=c("0-10cm","10-25cm", "25-50cm","50-100cm"))
  ggplot(rare_joint, aes(x = Sample_size, y = Species.x, color = Depth_Core, group = Site)) +
    geom_line() +
    scale_colour_manual("Depth", values = paletteCB4) +
    theme_classic() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) + 
    theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
    grids(linetype = "dashed") -> plot_rare
  return(plot_rare)
}


PlotRareCurve(rare_fungi, physeq_fungi_clean) +
  labs(title="Fungi", x= "Number of DNA reads", y= "Number of OTUs")

# *** FIGURE S2 - Rarefaction curves -------------------------------------------------------
ggarrange(PlotRareCurve(rare_fungi, physeq_fungi_clean) +
            labs(title="Fungi", x= "Number of DNA reads", y= "Number of OTUs") +
            theme(legend.position="none"), 
          PlotRareCurve(rare_prok, physeq_prok_clean) +
            labs(title="Prokaryotes", x= "Number of DNA reads", y= "Number of OTUs") +
            theme(legend.position = c(0.9, 0.2)),
          heights = c(1, 1),
          align = "hv",
          ncol = 2, 
          nrow = 1,
          common.legend = FALSE) -> Fig_S2_rare

Fig_S2_rare

# ********************************************----------------------------------------------
# BETA DIVERSITY ---------------------------------------------------------------------------
sort(sample_sums(physeq_fungi_clean))
sort(sample_sums(physeq_prok_clean))

# removed ~3% of samples with the lower library 1000/358 =2.793296
subset_samples(physeq_fungi_clean, sample_sums(physeq_fungi_clean) > 1781) -> physeq_fungi_filt
otu_table(physeq_fungi_filt) <-
  otu_table(physeq_fungi_filt)[which(rowSums(otu_table(physeq_fungi_filt)) > 0), ]
physeq_fungi_filt

subset_samples(physeq_prok_clean, sample_sums(physeq_prok_clean) > 6585 ) -> physeq_prok_filt
otu_table(physeq_prok_filt) <-
  otu_table(physeq_prok_filt)[which(rowSums(otu_table(physeq_prok_filt)) > 0), ]
physeq_prok_filt

# MetagenomeSeq normalization - Gaussian model ---------------------------------------------
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

CSSNorm(physeq_fungi_filt) -> physeq_fungi_mSeq
head(otu_table(physeq_fungi_mSeq))
head(sample_data(physeq_fungi_mSeq))

CSSNorm(physeq_prok_filt) -> physeq_prok_mSeq
head(otu_table(physeq_prok_mSeq))
head(sample_data(physeq_fungi_mSeq))

# Get a sense of the distribution of values in the dataset
count(as.data.frame(as.matrix(sample_data(physeq_fungi_mSeq))), Crop, Origin, Depth)
count(as.data.frame(as.matrix(sample_data(physeq_prok_mSeq))), Crop, Origin, Depth)

# recode mapping files
names(sample_data(physeq_fungi_mSeq))[names(sample_data(physeq_fungi_mSeq)) == 'Depth_Core'] <- 'Section'
sample_data(physeq_fungi_mSeq)$Section <- paste(sample_data(physeq_fungi_mSeq)$Section, "cm", sep = "")
head(sample_data(physeq_fungi_mSeq))

sample_data(physeq_fungi_mSeq)$Section <-
  factor(
    sample_data(physeq_fungi_mSeq)$Section,
    levels = c("0-10cm", "10-25cm", "25-50cm", "50-100cm"))

sample_data(physeq_fungi_mSeq)$Depth <-
  factor(
    sample_data(physeq_fungi_mSeq)$Depth,
    levels = c("10", "25", "50", "100"))

sample_data(physeq_fungi_mSeq)$Crop <-
  factor(
    sample_data(physeq_fungi_mSeq)$Crop,
    levels = c("Switchgrass", "Prairie","Poplar"))

names(sample_data(physeq_prok_mSeq))[names(sample_data(physeq_prok_mSeq)) == 'Depth_Core'] <- 'Section'
sample_data(physeq_prok_mSeq)$Section <- paste(sample_data(physeq_prok_mSeq)$Section, "cm", sep = "")
head(sample_data(physeq_prok_mSeq))

sample_data(physeq_prok_mSeq)$Section <-
  factor(
    sample_data(physeq_prok_mSeq)$Section,
    levels = c("0-10cm", "10-25cm", "25-50cm", "50-100cm"))

sample_data(physeq_prok_mSeq)$Depth <-
  factor(
    sample_data(physeq_prok_mSeq)$Depth,
    levels = c("10", "25", "50", "100"))

sample_data(physeq_prok_mSeq)$Crop <-
  factor(
    sample_data(physeq_prok_mSeq)$Crop,
    levels = c("Switchgrass", "Prairie","Poplar"))

# Extract dataset --------------------------------------------------------------------------
otu_fungi_mSeq <- t(as.data.frame(otu_table(physeq_fungi_mSeq)))
metadata_fungi_mSeq <- as.data.frame(as.matrix(sample_data(physeq_fungi_mSeq)))

otu_prok_mSeq <- t(as.data.frame(otu_table(physeq_prok_mSeq)))
metadata_prok_mSeq <- as.data.frame(as.matrix(sample_data(physeq_prok_mSeq)))

# Calculating 3D PCoA in vegan
cmdscale(vegdist(otu_fungi_mSeq, method = "bray"), k = 3, eig=TRUE) -> pcoa_k3_fungi
as.data.frame(pcoa_k3_fungi$points) -> pcoa_k3_fungi
colnames(pcoa_k3_fungi) <- c("Axis.1", "Axis.2", "Axis.3")
identical(rownames(pcoa_k3_fungi), rownames(metadata_fungi_mSeq))
str(pcoa_k3_fungi)
metadata_fungi_mSeq$Axis.1 <- pcoa_k3_fungi$Axis.1
metadata_fungi_mSeq$Axis.2 <- pcoa_k3_fungi$Axis.2
metadata_fungi_mSeq$Axis.3 <- pcoa_k3_fungi$Axis.3

library(plotly)
library(gganimate)

plot_ly(
  x = metadata_fungi_mSeq$Axis.1,
  y = metadata_fungi_mSeq$Axis.2,
  z = metadata_fungi_mSeq$Axis.3,
  type = "scatter3d",
  mode = "markers",
  color = metadata_fungi_mSeq$Origin) -> plot_fungi_gif


plot_fungi_gif %>%
  layout(scene = list(x = list(title = 'Axis.1'),
                      y = list(title = 'Axis.2'),
                      z = list(title = 'Axis.3')))

gganimate(plot_fungi_gif, "animatedgraph.gif")


cmdscale(vegdist(otu_prok_mSeq, method = "bray"), k = 3, eig=TRUE) -> pcoa_k3_prok
as.data.frame(pcoa_k3_prok$points) -> pcoa_k3_prok
colnames(pcoa_k3_prok) <- c("Axis.1", "Axis.2", "Axis.3")
identical(rownames(pcoa_k3_prok), rownames(metadata_prok_mSeq))
str(pcoa_k3_prok)
metadata_prok_mSeq$Axis.1 <- pcoa_k3_prok$Axis.1
metadata_prok_mSeq$Axis.2 <- pcoa_k3_prok$Axis.2
metadata_prok_mSeq$Axis.3 <- pcoa_k3_prok$Axis.3

plot_ly(
  x = metadata_prok_mSeq$Axis.1,
  y = metadata_prok_mSeq$Axis.2,
  z = metadata_prok_mSeq$Axis.3,
  type = "scatter3d",
  mode = "markers",
  color = metadata_prok_mSeq$Origin)



# pcoa using the ape package
require(ape)
pcoa(vegdist(otu_fungi_mSeq, method = "bray"), correction="none", rn=NULL) -> pcoa_fungi
identical(rownames(pcoa_fungi$vectors[,1:2]), rownames(metadata_fungi_mSeq))
metadata_fungi_mSeq <- cbind(metadata_fungi_mSeq, pcoa_fungi$vectors[,1:2])
head(pcoa_fungi$values)

pcoa(vegdist(otu_prok_mSeq, method = "bray"), correction="none", rn=NULL) -> pcoa_prok
identical(rownames(pcoa_prok$vectors[,1:2]), rownames(metadata_prok_mSeq))
metadata_prok_mSeq <- cbind(metadata_prok_mSeq, pcoa_prok$vectors[,1:2])
head(pcoa_prok$values)


# plot ordination 
PlotOrdin <-function(dataframe){
  ggplot(dataframe,
         aes(x=Axis.1,
             y=Axis.2,
             color = Depth,
             shape = Crop,
             fill = factor(ifelse(Origin == "Soil",Depth, NA)))) +
    geom_point(aes(size = Origin)) +
    scale_colour_manual("Depth", values = paletteCB4) +
    scale_shape_manual(name = "Crop", values = c(21,22,24)) +
    scale_fill_manual(name = "Niche", values = c(paletteCB4, "white")) +
    scale_size_manual(name = "Niche", values = c(1.5, 1.51)) +
    theme_classic() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
    theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
    theme(legend.position="bottom")+
    guides(
      fill = "none",
      size = guide_legend(override.aes = list(shape = c(1, 16)), nrow = 2),
      color = guide_legend(nrow = 2),
      shape = guide_legend(nrow = 2)) +
    grids(linetype = "dashed") -> plot_ord
  return(plot_ord)
}


PlotOrdin(metadata_fungi_mSeq) + labs(title = "Fungi")
PlotOrdin(metadata_prok_mSeq) + labs(title = "Prokaryote")



# PERMANOVA --------------------------------------------------------------------------------
GettingPerm <- function(dataframe, meta){
  adonis1 <- adonis(dataframe ~ Origin, 
                    meta, method = "bray",
                    permutations=9999)
  adonis2 <- adonis(dataframe ~ Crop, 
                    meta, method = "bray",
                    permutations=9999)
  adonis3 <- adonis(dataframe ~ Section, 
                    meta, method = "bray",
                    permutations=9999)
  data.frame(1:3, row.names = c("Origin","Crop","Depth")) -> df
  colnames(df) <- "R2"
  df$R2 <- c(
    round(adonis1$aov.tab$R2[1]*100, 1) ,
    round(adonis2$aov.tab$R2[1]*100, 1),
    round(adonis3$aov.tab$R2[1]*100, 1))
  df$Adj.p <-
    c(round(p.adjust(adonis1$aov.tab$`Pr(>F)`, "bonferroni"), 4)[1],
      round(p.adjust(adonis2$aov.tab$`Pr(>F)`, "bonferroni"), 4)[1],
      round(p.adjust(adonis3$aov.tab$`Pr(>F)`, "bonferroni"), 4)[1])
  return(df)
}


GettingPerm(otu_fungi_mSeq, metadata_fungi_mSeq) -> adonis_fungi
adonis_fungi

GettingPerm(otu_prok_mSeq, metadata_prok_mSeq) -> adonis_prok
adonis_prok

# PAIRWISE PERMANOVA -----------------------------------------------------------------------
library(mctoolsr) # post hoc test for permanova and other tools
library(RVAideMemoire)

mctoolsr::calc_pairwise_permanovas(vegdist(otu_fungi_mSeq, method="bray"), metadata_fungi_mSeq, "Origin")
mctoolsr::calc_pairwise_permanovas(vegdist(otu_fungi_mSeq, method="bray"), metadata_fungi_mSeq, "Crop")
mctoolsr::calc_pairwise_permanovas(vegdist(otu_fungi_mSeq, method="bray"), metadata_fungi_mSeq, "Depth")

# I like the RVAideMemoise output better
pairwise.perm.manova(vegdist(otu_fungi_mSeq, method="bray"), metadata_fungi_mSeq$Origin, p.method = "bonferroni")
pairwise.perm.manova(vegdist(otu_fungi_mSeq, method="bray"), metadata_fungi_mSeq$Crop, p.method = "bonferroni")
pairwise.perm.manova(vegdist(otu_fungi_mSeq, method="bray"), metadata_fungi_mSeq$Depth, p.method = "bonferroni")

pairwise.perm.manova(vegdist(otu_prok_mSeq, method="bray"), metadata_prok_mSeq$Origin, p.method = "bonferroni")
pairwise.perm.manova(vegdist(otu_prok_mSeq, method="bray"), metadata_prok_mSeq$Crop, p.method = "bonferroni")
pairwise.perm.manova(vegdist(otu_prok_mSeq, method="bray"), metadata_prok_mSeq$Depth, p.method = "bonferroni")

# *** FIGURE 2A - PCOA ---------------------------------------------------------------------
ggarrange(PlotOrdin(metadata_fungi_mSeq) + 
            labs(title = "Fungi") +
            annotate("text", -Inf, Inf, 
                     label = expression(paste("Origin ", italic(R) ^ 2,"= 11.0%***"),
                                        parse=TRUE), size = 2.15, hjust = -0.05, vjust = 1) +
            annotate("text", -Inf, Inf, 
                     label = expression(paste("Crop", italic(R) ^ 2,"= 3.7%***"),
                                        parse=TRUE), size = 2.15, hjust = -0.055, vjust = 1.8) +
            annotate("text", -Inf, Inf, 
                     label = expression(paste("Depth", italic(R) ^ 2,"= 7.2%***"),
                                        parse=TRUE), size = 2.15, hjust = -0.05, vjust = 2.6),
          PlotOrdin(metadata_prok_mSeq) + 
            labs(title = "Prokaryote") +
            annotate("text", -Inf, Inf, 
                     label = expression(paste("Origin ", italic(R) ^ 2,"= 25.9%***"),
                                        parse=TRUE), size = 2.15, hjust = -0.05, vjust = 1) +
            annotate("text", -Inf, Inf, 
                     label = expression(paste("Crop", italic(R) ^ 2,"= 2.3%***"),
                                        parse=TRUE), size = 2.15, hjust = -0.055, vjust = 1.8) +
            annotate("text", -Inf, Inf, 
                     label = expression(paste("Depth", italic(R) ^ 2,"= 10.1%***"),
                                        parse=TRUE), size = 2.15, hjust = -0.05, vjust = 2.6),
          heights = c(1, 1),
          labels = c("A", "B"),
          align = "hv",
          ncol = 1, 
          nrow = 2, 
          common.legend = TRUE, 
          legend = c("none")) -> plot_beta_div

plot_beta_div


# MULTIVARIATE DISPERSION ------------------------------------------------------------------
# Test the homogenity of group variances
library(multcompView)

pcoa(vegdist(otu_fungi_mSeq, method = "bray"), correction="none", rn=NULL) -> pcoa_fungi

permdisp_fungi_origin <- betadisper(vegdist(otu_fungi_mSeq, method="bray"), metadata_fungi_mSeq$Origin)

par(mfrow=c(1,2))
plot(permdisp_fungi_origin, main="PCoA")
boxplot(permdisp_fungi_origin, main="Distance to centroids")

boxplot(permdisp_fungi_origin$distances, permdisp_fungi_origin$group)
boxplot(permdisp_fungi_origin$distances ~ permdisp_fungi_origin$group)

permdisp_fungi_origin$group


par(mfrow=c(1,2))
plot(permdisp_prok_origin,main="PCoA")
boxplot(permdisp_prok_origin, main="Distance to centroids")
boxplot(permdisp$distances ~ permdisp$group)

anova(permdisp_fungi_origin, permutations = 9999)
permutest(permdisp_fungi_origin, permutations = 9999, pairwise = T) -> dist_fungi_origin
p.adjust(anova(permdisp_fungi_origin, permutations = 9999)$`Pr(>F)`, "bonferroni") # p= 0.0000000000000000000000000002661006 
data.frame(multcompLetters(p.adjust(dist_fungi_origin$pairwise$observed, 
                                    method="bonferroni"))['Letters']) -> pair_fungi_origin

permdisp_fungi_crop <- betadisper(vegdist(otu_fungi_mSeq, method = "bray"), metadata_fungi_mSeq$Crop) 
boxplot(permdisp_fungi_crop)
anova(permdisp_fungi_crop, permutations = 9999)
permutest(permdisp_fungi_crop, permutations = 9999, pairwise = T) -> dist_fungi_crop
p.adjust(anova(permdisp_fungi_crop, permutations = 9999)$`Pr(>F)`, "bonferroni") # p=0.002029742
data.frame(multcompLetters(p.adjust(dist_fungi_crop$pairwise$observed, 
                                    method="bonferroni"))['Letters']) -> pair_fungi_crop

permdisp_fungi_depth <- betadisper(vegdist(otu_fungi_mSeq, method="bray"), metadata_fungi_mSeq$Depth) 
boxplot(permdisp_fungi_depth)
anova(permdisp_fungi_depth, permutations = 9999)
permutest(permdisp_fungi_depth, permutations = 9999, pairwise = T) -> dist_fungi_depth
p.adjust(anova(permdisp_fungi_depth, permutations = 9999)$`Pr(>F)`, "bonferroni") # p=0.01688317
data.frame(multcompLetters(p.adjust(dist_fungi_depth$pairwise$observed, 
                                    method="bonferroni"))['Letters']) -> pair_fungi_depth

par(mfrow=c(1,2))
plot(permdisp_fungi_depth)
boxplot(permdisp_fungi_depth, main="Distance to centroids")

par(mfrow=c(1,2))
plot(permdisp_prok_depth)
boxplot(permdisp_prok_depth, main="Distance to centroids")


permdisp_prok_origin <- betadisper(vegdist(otu_prok_mSeq, method="bray"), metadata_prok_mSeq$Origin) 
boxplot(permdisp_prok_origin)
anova(permdisp_prok_origin, permutations = 9999)
permutest(permdisp_prok_origin, permutations = 9999, pairwise = T) -> dist_prok_origin
p.adjust(anova(permdisp_prok_origin, permutations = 9999)$`Pr(>F)`, "bonferroni") # p= 0.000000001663054
data.frame(multcompLetters(p.adjust(dist_prok_origin$pairwise$observed, 
                                    method="bonferroni"))['Letters']) -> pair_prok_origin

permdisp_prok_crop <- betadisper(vegdist(otu_prok_mSeq, method = "bray"), metadata_prok_mSeq$Crop) 
boxplot(permdisp_prok_crop)
anova(permdisp_prok_crop, permutations = 9999)
permutest(permdisp_prok_crop, permutations = 9999, pairwise = T) -> dist_prok_crop
p.adjust(anova(permdisp_prok_crop, permutations = 9999)$`Pr(>F)`, "bonferroni") # p=0.03036966
data.frame(multcompLetters(p.adjust(dist_prok_crop$pairwise$observed, 
                                    method="bonferroni"))['Letters']) -> pair_prok_crop

permdisp_prok_depth <- betadisper(vegdist(otu_prok_mSeq, method="bray"), metadata_prok_mSeq$Depth) 
boxplot(permdisp_prok_depth)
anova(permdisp_prok_depth, permutations = 9999)
permutest(permdisp_prok_depth, permutations = 9999, pairwise = T) -> dist_prok_depth
p.adjust(anova(permdisp_prok_depth, permutations = 9999)$`Pr(>F)`, "bonferroni") # p=0.000000018291
data.frame(multcompLetters(p.adjust(dist_prok_depth$pairwise$observed, 
                                    method="bonferroni"))['Letters']) -> pair_prok_depth


# max(dataframe$Observed + 0.1* max(dataframe$Observed)) -> labels_y
# stat_summary(geom = 'text', label = my_labels, fun= max, aes(y = labels_y), size=3, color="black") +

# plot multivariate dispersion 
PlotBetadisper <- function(betadisp, value, my_labels, metadata){
  # creating a label and dataframe
  max(betadisp$distances + 0.1* max(betadisp$distances)) -> labels_y
  #labels_y = 0.9
  data.frame(betadisp$group, betadisp$distances) -> df
  colnames(df) <- c("value", "distance")
  if (identical(rownames(df), rownames(metadata))==TRUE){
    df$Crop <- metadata$Crop
    df$Origin <- metadata$Origin
    df$Depth <- metadata$Depth
    # plotting
    ggplot(df, aes(x=value, y=distance)) +
      geom_jitter(aes(
        fill = factor(ifelse(Origin == "Soil", Depth, NA)),
        size = Origin,
        shape = Crop,
        colour = Depth),
        alpha = 0.8) +
      geom_boxplot(outlier.colour="black", outlier.shape = 8, alpha=0.6, lwd = 0.3) +
      stat_summary(fun=mean, geom="point", shape=18, size=1.9, color="red", fill="red") +
      stat_summary(geom = 'text', label = my_labels, fun= max, aes(y = labels_y), size=3, color="black") +
      scale_fill_manual(name = "Niche", values = c(paletteCB4, "white")) +
      scale_size_manual(name = "Niche", values = c(0.9, 0.91)) +
      scale_shape_manual(name = "Crop", values = c(21,22,24)) +
      scale_colour_manual("Depth", values = paletteCB4) +
      theme_classic() +
      theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
      theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
      theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1, vjust = 1)) +
      theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
      theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
      theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
      theme(legend.position="bottom")+
      guides(
        size = guide_legend(title = "Niche", override.aes = list(shape = c(1, 16)), nrow = 2, order = 1),
        fill = "none",
        shape = guide_legend(title = "Crop", nrow = 2, order = 2), 
        color = guide_legend(title = "Depth", nrow = 2, order = 3)) +
      grids(linetype = "dashed") -> betaplot
    return(betaplot)
  }else{
    stop("Error: dataframe and metadata rownames are not matching!")
  }
}


# permdisp_fungi_depth$group <- factor(permdisp_fungi_depth$group, levels = c("10", "25", "50", "100"))
# permdisp_prok_depth$group <-factor(permdisp_prok_depth$group, levels = c("10", "25", "50", "100"))
names(metadata_fungi_mSeq)[names(metadata_fungi_mSeq) == "Depth"] <- "Bottom"
names(metadata_prok_mSeq)[names(metadata_prok_mSeq) == "Depth"] <- "Bottom"
names(metadata_fungi_mSeq)[names(metadata_fungi_mSeq) == "Section"] <- "Depth"
names(metadata_prok_mSeq)[names(metadata_prok_mSeq) == "Section"] <- "Depth"


PlotBetadisper(permdisp_fungi_origin, "Origin", 
               as.character(pair_fungi_origin$Letters), metadata_fungi_mSeq) + 
  labs(title = "Niche", y ="Distance to centroids", x= "")


# extracting the legend for plotting 
get_legend(
  PlotBetadisper(permdisp_fungi_origin,"Origin",
                 as.character(pair_fungi_origin$Letters),
                 metadata_fungi_mSeq)) -> legend_beta_div

as_ggplot(legend_beta_div)

# *** FIGURE 2B - Beta Dispersion ----------------------------------------------------------
ggarrange(PlotBetadisper(permdisp_fungi_origin, "Origin", as.character(pair_fungi_origin$Letters), metadata_fungi_mSeq) + 
            labs(title = "Niche", y ="Distance to centroid", x= NULL),
          PlotBetadisper(permdisp_fungi_crop, "Crop", as.character(pair_fungi_crop$Letters), metadata_fungi_mSeq) + 
            labs(title = "Crop", y ="Distance to centroid", x= NULL),
          PlotBetadisper(permdisp_fungi_depth, "Depth", as.character(pair_fungi_depth[c(1,3,4,2),]), metadata_fungi_mSeq) + 
            labs(title = "Depth", y ="Distance to centroid", x= NULL) +
            scale_x_discrete(breaks=c("10", "25", "50", "100"), labels= c("0-10cm", "10-25cm", "25-50cm", "50-100cm")),
          PlotBetadisper(permdisp_prok_origin, "Origin", as.character(pair_prok_origin$Letters), metadata_prok_mSeq) + 
            labs(title = "Niche", y ="Distance to centroid", x= NULL),
          PlotBetadisper(permdisp_prok_crop, "Crop", as.character(pair_prok_crop$Letters), metadata_prok_mSeq) + 
            labs(title = "Crop", y ="Distance to centroid", x= NULL),
          PlotBetadisper(permdisp_prok_depth, "Depth", as.character(pair_prok_depth[c(1,3,4,2),]), metadata_prok_mSeq) + 
            labs(title = "Depth", y ="Distance to centroid", x= NULL) +
            scale_x_discrete(breaks=c("10", "25", "50", "100"), labels= c("0-10cm", "10-25cm", "25-50cm", "50-100cm")),
          widths = c(1,1,1,1,1,1),
          labels = c("C","","","D","",""),
          align = "hv",
          ncol = 3, 
          nrow = 2,
          common.legend = TRUE, 
          legend = c("none")) -> betadisp_plot

betadisp_plot

# *** FIGURE 2 - complete ------------------------------------------------------------------

Fig2_beta_div <-
ggarrange(
  ggarrange(
    plot_beta_div,
    betadisp_plot,
    widths = c(0.6, 1),
    ncol = 2),
  as_ggplot(legend_beta_div),
  nrow = 2, 
  heights = c(1, 0.06))

Fig2_beta_div


# Importing soil data --------------------------------------------------------------------
head(sample_data(physeq_fungi_mSeq))
sample_data(physeq_fungi_mSeq)[sample_data(physeq_fungi_mSeq)$Crop=="Switchgrass" & # see the error here?
                                 sample_data(physeq_fungi_mSeq)$Plot=="R1" & 
                                 sample_data(physeq_fungi_mSeq)$Section=="10-25cm",]

sample_data(physeq_fungi_mSeq)[rownames(sample_data(physeq_fungi_mSeq))=="SAM34",]$ph <- "7.3"
sample_data(physeq_fungi_mSeq)[rownames(sample_data(physeq_fungi_mSeq))=="SAM34",]$P <- "20"
sample_data(physeq_fungi_mSeq)[rownames(sample_data(physeq_fungi_mSeq))=="SAM34",]$K <- "59"
sample_data(physeq_fungi_mSeq)[rownames(sample_data(physeq_fungi_mSeq))=="SAM34",]$Ca <- "944"
sample_data(physeq_fungi_mSeq)[rownames(sample_data(physeq_fungi_mSeq))=="SAM34",]$Mg <- "119"
sample_data(physeq_fungi_mSeq)[rownames(sample_data(physeq_fungi_mSeq))=="SAM34",]$CEC <- "5.86295"

# add a new variable 
sample_data(physeq_fungi_mSeq)$Filt <- paste(sample_data(physeq_fungi_mSeq)$Crop, 
                                             sample_data(physeq_fungi_mSeq)$Plot,
                                             sample_data(physeq_fungi_mSeq)$Section,
                                             sep = "_") # add a composite variable
head(sample_data(physeq_fungi_mSeq))

head(sample_data(physeq_prok_mSeq))
sample_data(physeq_prok_mSeq)[sample_data(physeq_prok_mSeq)$Crop=="Switchgrass" & # see the error here?
                                sample_data(physeq_prok_mSeq)$Plot=="R1" & 
                                sample_data(physeq_prok_mSeq)$Section=="10-25cm",]

sample_data(physeq_prok_mSeq)[rownames(sample_data(physeq_prok_mSeq))=="SAM30",]$ph <- "7.3"
sample_data(physeq_prok_mSeq)[rownames(sample_data(physeq_prok_mSeq))=="SAM30",]$P <- "20"
sample_data(physeq_prok_mSeq)[rownames(sample_data(physeq_prok_mSeq))=="SAM30",]$K <- "59"
sample_data(physeq_prok_mSeq)[rownames(sample_data(physeq_prok_mSeq))=="SAM30",]$Ca <- "944"
sample_data(physeq_prok_mSeq)[rownames(sample_data(physeq_prok_mSeq))=="SAM30",]$Mg <- "119"
sample_data(physeq_prok_mSeq)[rownames(sample_data(physeq_prok_mSeq))=="SAM30",]$CEC <- "5.86295"

# add a new variable 
sample_data(physeq_prok_mSeq)$Filt <- paste(sample_data(physeq_prok_mSeq)$Crop, 
                                            sample_data(physeq_prok_mSeq)$Plot,
                                            sample_data(physeq_prok_mSeq)$Section,
                                            sep = "_") # add a composite variable
head(sample_data(physeq_prok_mSeq))

# ********************************************----------------------------------------------
# CAP MODELS -------------------------------------------------------------------------------
# Mergin pseudorep to match chemistry data -----------------------------------------------
# Divide root from soil samples 
MergeData <- function(physeq, Var){
  sub_set <- subset(sample_data(physeq), Origin == Var)
  physeq_filt <- merge_phyloseq(otu_table(physeq),
                                tax_table(physeq),
                                refseq(physeq),
                                sub_set)
  #sample_data(subset_samples(physeq_ITS_scaled, Type=="in")
  otu_table(physeq_filt) <-
    otu_table(physeq_filt)[which(rowSums(otu_table(physeq_filt)) > 0),]
  physeq_filt %T>% print()
  #merge data
  physeq1 <- merge_samples(physeq_filt, "Filt", fun = "sum") # mean or sum depends on the goal
  sample_data(physeq1) <- sample_data(physeq1)[, 13:18]
  sample_data(physeq1)$Filt <- rownames(sample_data(physeq1))
  #head(otu_table(physeq1)) %T>% print()
  # re add variables
  sample_data(physeq1)$Crop <-
    gsub("\\_.*", "", sample_data(physeq1)$Filt)
  sample_data(physeq1)$Depth <-
    gsub(".*_", "", sample_data(physeq1)$Filt)
  sample_data(physeq1)
  return(physeq1)
}


MergeData(physeq_fungi_mSeq, "Root") ->  physeq_fungi_Root_filt
physeq_fungi_mSeq
sample_data(physeq_fungi_mSeq)
sample_data(physeq_fungi_Root_filt)
otu_table(physeq_fungi_Root_filt)

MergeData(physeq_fungi_mSeq, "Soil") ->  physeq_fungi_Soil_filt
physeq_fungi_Root_filt

MergeData(physeq_prok_mSeq, "Root") ->  physeq_prok_Root_filt
physeq_prok_Root_filt

MergeData(physeq_prok_mSeq, "Soil") ->  physeq_prok_Soil_filt
physeq_prok_Soil_filt


# 1) Fungi Root ----------------------------------------------------------------------------
otu_fungi_Root <- as.data.frame(otu_table(physeq_fungi_Root_filt))
metadata_fungi_Root <- as.data.frame(as.matrix(sample_data(physeq_fungi_Root_filt)))

cap_fungi_Root_full <- vegan::capscale(otu_fungi_Root ~ Crop + Depth + Crop:Depth, metadata_fungi_Root,distance = "bray")
cap_fungi_Root_full
anova(cap_fungi_Root_full, permutations = 9999, by = "terms")
#p.adjust(anova(cap_fungi_Root_full, permutations = 9999)$`Pr(>F)`, "bonferroni")

pars_model_fungi_Root <- vegan::ordistep(cap_fungi_Root_full, direction = "backward", Pout = 0.05, permutations = how(nperm = 9999))
pars_model_fungi_Root <- formula(pars_model_fungi_Root)[1:3]
pars_model_fungi_Root

# Now extracting variances for each term
cap_fungi_Root_Crop <-capscale(otu_fungi_Root ~ Crop, metadata_fungi_Root, distance = "bray")
RsquareAdj(cap_fungi_Root_Crop)
anova(cap_fungi_Root_Crop, permutations = 9999)
cap_fungi_Root_Depth <-capscale(otu_fungi_Root ~ Depth, metadata_fungi_Root, distance = "bray")
RsquareAdj(cap_fungi_Root_Depth)
anova(cap_fungi_Root_Depth, permutations = 9999)
cap_fungi_Root_int <-capscale(otu_fungi_Root ~ Crop:Depth + Condition(Crop) + Condition(Depth), 
                              metadata_fungi_Root, distance = "bray")
RsquareAdj(cap_fungi_Root_int)
anova(cap_fungi_Root_int, permutations = 9999)

RsquareAdj(cap_fungi_Root_full)

sum(unlist(RsquareAdj(cap_fungi_Root_Crop)[1]),unlist(RsquareAdj(cap_fungi_Root_Depth)[1]),
    unlist(RsquareAdj(cap_fungi_Root_int)[1]))
sum(unlist(RsquareAdj(cap_fungi_Root_Crop)[2]),unlist(RsquareAdj(cap_fungi_Root_Depth)[2]),
    unlist(RsquareAdj(cap_fungi_Root_int)[2]))

# prepare for plotting
identical(rownames(scores(cap_fungi_Root_full)$sites), rownames(metadata_fungi_Root))
metadata_fungi_Root <- cbind(metadata_fungi_Root, scores(cap_fungi_Root_full)$sites)
head(metadata_fungi_Root)

# RDA models using Hellinger transformed data ----------------------------------------------
# Bray-Curtis metric is insensitive to double absences (i.e., an
# OTU has a 0 count in any pair of samples). A metric that accounts for
# double absences will group samples sharing 0 for an OTU closer together
# than Bray-Curtis will to. If some of your OTUs are absent in a certain
# soil type or do not like a certain plant species, then it would make
# sense to take double absences into account.

# add a new variable to the raw dataframe
sample_data(physeq_fungi_filt)$Filt <- paste(sample_data(physeq_fungi_filt)$Crop, 
                                             sample_data(physeq_fungi_filt)$Plot,
                                             sample_data(physeq_fungi_filt)$Depth_Core,
                                             sep = "_") # add a composite variable
head(sample_data(physeq_fungi_filt))

MergeData(physeq_fungi_filt, "Root") -> physeq_fungi_raw_Root
otu_table(physeq_fungi_raw_Root)[1:10, 1:10]

otu_fungi_raw_Root = as(otu_table(physeq_fungi_raw_Root), "matrix")
otu_fungi_raw_Root = as.data.frame(otu_fungi_raw_Root)
metadata_fungi_raw_Root <- as.data.frame(as.matrix(sample_data(physeq_fungi_raw_Root)))
otu_fungi_raw_Root[1:10, 1:10]

otu_fungi_raw_Root_hel <- decostand(otu_fungi_raw_Root, method = "hellinger")
otu_fungi_raw_Root_hel[1:10, 1:10]

rda_fungi_Root_full <- vegan::rda(otu_fungi_raw_Root_hel ~ Crop + Depth + Crop:Depth, metadata_fungi_raw_Root)
rda_fungi_Root_full
anova(rda_fungi_Root_full, permutations = 9999, by = "terms")

pars_model_fungi_Root_hel <- vegan::ordistep(rda_fungi_Root_full, direction = "backward", Pout = 0.05, permutations = how(nperm = 9999))
pars_model_fungi_Root_hel <- formula(pars_model_fungi_Root_hel)[1:3]
pars_model_fungi_Root_hel

# Now extracting variances for each term
rda_fungi_Root_Crop <-rda(otu_fungi_raw_Root_hel ~ Crop, metadata_fungi_raw_Root)
RsquareAdj(rda_fungi_Root_Crop)
rda_fungi_Root_Depth <-rda(otu_fungi_raw_Root_hel ~ Depth, metadata_fungi_raw_Root)
RsquareAdj(rda_fungi_Root_Depth)
rda_fungi_Root_int <-rda(otu_fungi_raw_Root_hel ~ Crop:Depth + Condition(Crop) + Condition(Depth), metadata_fungi_raw_Root)
RsquareAdj(rda_fungi_Root_int)

RsquareAdj(rda_fungi_Root_full)

sum(unlist(RsquareAdj(rda_fungi_Root_Crop)[1]),unlist(RsquareAdj(rda_fungi_Root_Depth)[1]),
    unlist(RsquareAdj(rda_fungi_Root_int)[1]))

sum(unlist(RsquareAdj(rda_fungi_Root_Crop)[2]),unlist(RsquareAdj(rda_fungi_Root_Depth)[2]),
    unlist(RsquareAdj(rda_fungi_Root_int)[2]))

# Testing for dispersion
anova(betadisper(vegdist(otu_fungi_raw_Root_hel, method = "euclidean"), metadata_fungi_raw_Root$Depth), permutations = 9999)
anova(betadisper(vegdist(otu_fungi_raw_Root_hel, method = "euclidean"), metadata_fungi_raw_Root$Crop), permutations = 9999)

# prepare for plotting
identical(rownames(scores(rda_fungi_Root_full)$sites), rownames(metadata_fungi_raw_Root))
metadata_fungi_raw_Root <- cbind(metadata_fungi_raw_Root, scores(rda_fungi_Root_full)$sites)
head(metadata_fungi_raw_Root)

# protest ----------------------------------------------------------------------------------
protest(
  cap_fungi_Root_full,
  rda_fungi_Root_full,
  scores = "sites",
  permutations = how(nperm = 9999))


# 2) Fungi Soil ----------------------------------------------------------------------------
otu_fungi_Soil <- as.data.frame(otu_table(physeq_fungi_Soil_filt))
metadata_fungi_Soil <- as.data.frame(as.matrix(sample_data(physeq_fungi_Soil_filt)))

cap_fungi_Soil_full <- vegan::capscale(otu_fungi_Soil ~ Crop + Depth + Crop:Depth, metadata_fungi_Soil, distance = "bray")
anova(cap_fungi_Soil_full, permutations = 9999, by = "terms")

pars_model_fungi_Soil <- vegan::ordistep(cap_fungi_Soil_full, direction = "backward", Pout = 0.05, permutations = how(nperm = 9999))
pars_model_fungi_Soil <- formula(pars_model_fungi_Soil)[1:3]
pars_model_fungi_Soil

# re-run with the most parsimonious model
cap_fungi_Soil_pars <- vegan::capscale(pars_model_fungi_Soil, metadata_fungi_Soil, distance = "bray")
anova(cap_fungi_Soil_pars, permutations = 9999, by = "terms")

# Now extracting variances for each term
cap_fungi_Soil_Crop <-capscale(otu_fungi_Soil ~ Crop, metadata_fungi_Soil, distance = "bray")
RsquareAdj(cap_fungi_Soil_Crop)
anova(cap_fungi_Soil_Crop, permutations = 9999)
cap_fungi_Soil_Depth <-capscale(otu_fungi_Soil ~ Depth, metadata_fungi_Soil, distance = "bray")
RsquareAdj(cap_fungi_Soil_Depth)
anova(cap_fungi_Soil_Depth, permutations = 9999)

RsquareAdj(cap_fungi_Soil_pars)

sum(unlist(RsquareAdj(cap_fungi_Soil_Crop)[1]),unlist(RsquareAdj(cap_fungi_Soil_Depth)[1]))
sum(unlist(RsquareAdj(cap_fungi_Soil_Crop)[2]),unlist(RsquareAdj(cap_fungi_Soil_Depth)[2]))

# Testing for dispersion
anova(betadisper(vegdist(otu_fungi_Soil, method="bray"), metadata_fungi_Soil$Depth), permutations = 9999)
anova(betadisper(vegdist(otu_fungi_Soil, method="bray"), metadata_fungi_Soil$Crop), permutations = 9999)

identical(rownames(scores(cap_fungi_Soil_pars)$sites), rownames(metadata_fungi_Soil))
metadata_fungi_Soil <- cbind(metadata_fungi_Soil, scores(cap_fungi_Soil_pars)$sites)
head(metadata_fungi_Soil)

# RDA models using Hellinger transformed data ----------------------------------------------
MergeData(physeq_fungi_filt, "Soil") -> physeq_fungi_raw_Soil
otu_table(physeq_fungi_raw_Soil)[1:10, 1:10]

otu_fungi_raw_Soil = as(otu_table(physeq_fungi_raw_Soil), "matrix")
otu_fungi_raw_Soil = as.data.frame(otu_fungi_raw_Soil)
metadata_fungi_raw_Soil <- as.data.frame(as.matrix(sample_data(physeq_fungi_raw_Soil)))
otu_fungi_raw_Soil[1:10, 1:10]

otu_fungi_raw_Soil_hel <- decostand(otu_fungi_raw_Soil, method = "hellinger")
otu_fungi_raw_Soil_hel[1:10, 1:10]

rda_fungi_Soil_full <- vegan::rda(otu_fungi_raw_Soil_hel ~ Crop + Depth + Crop:Depth, metadata_fungi_raw_Soil)
rda_fungi_Soil_full
anova(rda_fungi_Soil_full, permutations = 9999, by = "terms")

pars_model_fungi_Soil_hel <- vegan::ordistep(rda_fungi_Soil_full, direction = "backward", Pout = 0.05, permutations = how(nperm = 9999))
pars_model_fungi_Soil_hel <- formula(pars_model_fungi_Soil_hel)[1:3]
pars_model_fungi_Soil_hel

# re-run with the most parsimonious model
rda_fungi_Soil_pars <- vegan::rda(pars_model_fungi_Soil_hel, metadata_fungi_raw_Soil)
rda_fungi_Soil_pars
anova(rda_fungi_Soil_pars, permutations = 9999, by = "terms")

# Now extracting variances for each term
rda_fungi_Soil_Crop <-rda(otu_fungi_raw_Soil_hel ~ Crop, metadata_fungi_raw_Soil)
RsquareAdj(rda_fungi_Soil_Crop)
rda_fungi_Soil_Depth <-rda(otu_fungi_raw_Soil_hel ~ Depth, metadata_fungi_raw_Soil)
RsquareAdj(rda_fungi_Soil_Depth)

RsquareAdj(rda_fungi_Soil_full)
RsquareAdj(rda_fungi_Soil_pars)

sum(unlist(RsquareAdj(rda_fungi_Soil_Crop)[1]),unlist(RsquareAdj(rda_fungi_Soil_Depth)[1]))
sum(unlist(RsquareAdj(rda_fungi_Soil_Crop)[2]),unlist(RsquareAdj(rda_fungi_Soil_Depth)[2]))

# Testing for dispersion
anova(betadisper(vegdist(otu_fungi_raw_Soil_hel, method = "euclidean"), metadata_fungi_raw_Soil$Depth), permutations = 9999)
anova(betadisper(vegdist(otu_fungi_raw_Soil_hel, method = "euclidean"), metadata_fungi_raw_Soil$Crop), permutations = 9999)

# prepare for plotting
identical(rownames(scores(rda_fungi_Soil_pars)$sites), rownames(metadata_fungi_raw_Soil))
metadata_fungi_raw_Soil <- cbind(metadata_fungi_raw_Soil, scores(rda_fungi_Soil_pars)$sites)
head(metadata_fungi_raw_Soil)

# protest ----------------------------------------------------------------------------------
protest(
  cap_fungi_Soil_pars,
  rda_fungi_Soil_pars,
  scores = "sites",
  permutations = how(nperm = 9999))


# 3) Prokaryote Root -----------------------------------------------------------------------
otu_prok_Root = as(otu_table(physeq_prok_Root_filt), "matrix")
otu_prok_Root = as.data.frame(otu_prok_Root)
otu_prok_Root[1:10, 1:10]

metadata_prok_Root <- as.data.frame(as.matrix(sample_data(physeq_prok_Root_filt)))

cap_prok_Root_full <- vegan::capscale(otu_prok_Root ~ Crop + Depth + Crop:Depth, metadata_prok_Root,distance = "bray")
cap_prok_Root_full
anova(cap_prok_Root_full, permutations = 9999, by = "terms")

pars_model_prok_Root <- vegan::ordistep(cap_prok_Root_full, direction = "backward", Pout = 0.05, permutations = how(nperm = 9999))
pars_model_prok_Root <- formula(pars_model_prok_Root)[1:3]
pars_model_prok_Root

# Now extracting variances for each term
cap_prok_Root_Crop <-capscale(otu_prok_Root ~ Crop, metadata_prok_Root, distance = "bray")
RsquareAdj(cap_prok_Root_Crop)
anova(cap_prok_Root_Crop, permutations = 9999)
cap_prok_Root_Depth <-capscale(otu_prok_Root ~ Depth, metadata_prok_Root, distance = "bray")
RsquareAdj(cap_prok_Root_Depth)
anova(cap_prok_Root_Depth, permutations = 9999)
cap_prok_Root_int <-capscale(otu_prok_Root ~ Crop:Depth + Condition(Crop) + Condition(Depth), 
                             metadata_prok_Root, distance = "bray")
RsquareAdj(cap_prok_Root_int)
anova(cap_prok_Root_int, permutations = 9999)

RsquareAdj(cap_prok_Root_full)

sum(unlist(RsquareAdj(cap_prok_Root_Crop)[1]),unlist(RsquareAdj(cap_prok_Root_Depth)[1]),
    unlist(RsquareAdj(cap_prok_Root_int)[1]))

sum(unlist(RsquareAdj(cap_prok_Root_Crop)[2]),unlist(RsquareAdj(cap_prok_Root_Depth)[2]),
    unlist(RsquareAdj(cap_prok_Root_int)[2]))

# Testing for dispersion
anova(betadisper(vegdist(otu_prok_Root, method="bray"), metadata_prok_Root$Depth), permutations = 9999)
anova(betadisper(vegdist(otu_prok_Root, method="bray"), metadata_prok_Root$Crop), permutations = 9999)

# prepare for plotting
identical(rownames(scores(cap_prok_Root_full)$sites), rownames(metadata_prok_Root))
metadata_prok_Root <- cbind(metadata_prok_Root, scores(cap_prok_Root_full)$sites)
head(metadata_prok_Root)

# RDA models using Hellinger transformed data ----------------------------------------------
# add a new variable to the raw dataframe
sample_data(physeq_prok_filt)$Filt <- paste(sample_data(physeq_prok_filt)$Crop, 
                                            sample_data(physeq_prok_filt)$Plot,
                                            sample_data(physeq_prok_filt)$Depth_Core,
                                            sep = "_") # add a composite variable
head(sample_data(physeq_prok_filt))

MergeData(physeq_prok_filt, "Root") -> physeq_prok_raw_Root
otu_table(physeq_prok_raw_Root)[1:10, 1:10]

otu_prok_raw_Root = as(otu_table(physeq_prok_raw_Root), "matrix")
otu_prok_raw_Root = as.data.frame(otu_prok_raw_Root)
metadata_prok_raw_Root <- as.data.frame(as.matrix(sample_data(physeq_prok_raw_Root)))
otu_prok_raw_Root[1:10, 1:10]

otu_prok_raw_Root_hel <- decostand(otu_prok_raw_Root, method = "hellinger")
otu_prok_raw_Root_hel[1:10, 1:10]

rda_prok_Root_full <- vegan::rda(otu_prok_raw_Root_hel ~ Crop + Depth + Crop:Depth, metadata_prok_raw_Root)
rda_prok_Root_full
anova(rda_prok_Root_full, permutations = 9999, by = "terms")

pars_model_prok_Root_hel <- vegan::ordistep(rda_prok_Root_full, direction = "backward", Pout = 0.05, permutations = how(nperm = 9999))
pars_model_prok_Root_hel <- formula(pars_model_prok_Root_hel)[1:3]
pars_model_prok_Root_hel

# Now extracting variances for each term
rda_prok_Root_Crop <-rda(otu_prok_raw_Root_hel ~ Crop, metadata_prok_raw_Root)
RsquareAdj(rda_prok_Root_Crop)
rda_prok_Root_Depth <-rda(otu_prok_raw_Root_hel ~ Depth, metadata_prok_raw_Root)
RsquareAdj(rda_prok_Root_Depth)
rda_prok_Root_int <-rda(otu_prok_raw_Root_hel ~ Crop:Depth + Condition(Crop) + Condition(Depth), metadata_prok_raw_Root)
RsquareAdj(rda_prok_Root_int)

RsquareAdj(rda_prok_Root_full)

sum(unlist(RsquareAdj(rda_prok_Root_Crop)[1]),unlist(RsquareAdj(rda_prok_Root_Depth)[1]),
    unlist(RsquareAdj(rda_prok_Root_int)[1]))

sum(unlist(RsquareAdj(rda_prok_Root_Crop)[2]),unlist(RsquareAdj(rda_prok_Root_Depth)[2]),
    unlist(RsquareAdj(rda_prok_Root_int)[2]))

# Testing for dispersion
anova(betadisper(vegdist(otu_prok_raw_Root_hel, method = "euclidean"), metadata_prok_raw_Root$Depth), permutations = 9999)
anova(betadisper(vegdist(otu_prok_raw_Root_hel, method = "euclidean"), metadata_prok_raw_Root$Crop), permutations = 9999)

# prepare for plotting
identical(rownames(scores(rda_prok_Root_full)$sites), rownames(metadata_prok_raw_Root))
metadata_prok_raw_Root <- cbind(metadata_prok_raw_Root, scores(rda_prok_Root_full)$sites)
head(metadata_prok_raw_Root)

# protest ----------------------------------------------------------------------------------
protest(
  cap_prok_Root_full,
  rda_prok_Root_full,
  scores = "sites",
  permutations = how(nperm = 9999))


# 4) Prokaryote Soil -----------------------------------------------------------------------
otu_prok_Soil = as(otu_table(physeq_prok_Soil_filt), "matrix")
otu_prok_Soil = as.data.frame(otu_prok_Soil)
otu_prok_Soil[1:10, 1:10]

metadata_prok_Soil <- as.data.frame(as.matrix(sample_data(physeq_prok_Soil_filt)))
head(metadata_prok_Soil)

cap_prok_Soil_full <- vegan::capscale(otu_prok_Soil ~ Crop + Depth + Crop:Depth, metadata_prok_Soil, distance = "bray")
cap_prok_Soil_full
anova(cap_prok_Soil_full, permutations = 9999, by = "terms")

pars_model_prok_Soil <- vegan::ordistep(cap_prok_Soil_full, direction = "backward", Pout = 0.05, permutations = how(nperm = 9999))
pars_model_prok_Soil <- formula(pars_model_prok_Soil)[1:3]
pars_model_prok_Soil

# Now extracting variances for each term
cap_prok_Soil_Crop <-capscale(otu_prok_Soil ~ Crop, metadata_prok_Soil, distance = "bray")
RsquareAdj(cap_prok_Soil_Crop)
anova(cap_prok_Soil_Crop, permutations = 9999)
cap_prok_Soil_Depth <-capscale(otu_prok_Soil ~ Depth, metadata_prok_Soil, distance = "bray")
RsquareAdj(cap_prok_Soil_Depth)
anova(cap_prok_Soil_Depth, permutations = 9999)
cap_prok_Soil_int <-capscale(otu_prok_Soil ~ Crop:Depth + Condition(Crop) + Condition(Depth), 
                             metadata_prok_Soil, distance = "bray")
RsquareAdj(cap_prok_Soil_int)
anova(cap_prok_Soil_int, permutations = 9999)

RsquareAdj(cap_prok_Soil_full)

sum(unlist(RsquareAdj(cap_prok_Soil_Crop)[1]),unlist(RsquareAdj(cap_prok_Soil_Depth)[1]),
    unlist(RsquareAdj(cap_prok_Soil_int)[1]))
sum(unlist(RsquareAdj(cap_prok_Soil_Crop)[2]),unlist(RsquareAdj(cap_prok_Soil_Depth)[2]),
    unlist(RsquareAdj(cap_prok_Soil_int)[1]))

# Testing for dispersion
anova(betadisper(vegdist(otu_prok_Soil, method="bray"), metadata_prok_Soil$Depth), permutations = 9999)
anova(betadisper(vegdist(otu_prok_Soil, method="bray"), metadata_prok_Soil$Crop), permutations = 9999)

identical(rownames(scores(cap_prok_Soil_full)$sites), rownames(metadata_prok_Soil))
metadata_prok_Soil <- cbind(metadata_prok_Soil, scores(cap_prok_Soil_full)$sites)
head(metadata_prok_Soil)

#************************************************************************************-------
# ENVFIT -----------------------------------------------------------------------------------
# For continuous environmental variables it actually uses ordination to predict the environmental variable.
# Under the hood (bonnet), it fits a first degree trend surface (plane in 2D) for environmental variable 
# over the ordination scores, and the R2 is the proportion that surface explains of the variable. 
# The arrow shown is the gradient of this fitted trend surface and shows the direction to which the
# variable changes most rapidly in a first degree linear model.

# Calculating fit of environemtnal variables -----------------------------------------------
GenerateEnvFit <- function(cap, metadata){
  metadata$ph <- as.numeric(as.character(metadata$pH))
  #names(metadata)[names(metadata) == "ph"] <- "pH"
  metadata$P <- as.numeric(as.character(metadata$P))
  metadata$K <- as.numeric(as.character(metadata$K))
  metadata$Ca <- as.numeric(as.character(metadata$Ca))
  metadata$Mg <- as.numeric(as.character(metadata$Mg))
  metadata$CEC <- as.numeric(as.character(metadata$CEC))
  metadata$C <- as.numeric(as.character(metadata$C))
  metadata$N <- as.numeric(as.character(metadata$N))
  metadata$Sand <- as.numeric(as.character(metadata$Sand))
  metadata$Clay <- as.numeric(as.character(metadata$Clay))
  metadata$Silt <- as.numeric(as.character(metadata$Silt))
  envfit(cap,
         metadata[, 1:11],
         perm = 9999,
         display = "lc") -> res_envfit
  df_envfit <- as.data.frame(scores(res_envfit, display = "vectors"))
  df_envfit <- cbind(df_envfit, Var = rownames(df_envfit))
  df_envfit$pvals <-res_envfit$vectors$pvals
  df_envfit <- subset(df_envfit, df_envfit$pvals<=0.05)
  df_envfit
  return(df_envfit)
}

data_chem_new <-
  data_chem %>% 
  mutate(Merge = paste(Crop, Plot, sep = "_"),
         Filt = paste(Merge, Depth, sep = "_")) %>% 
  select(-Merge, -Crop, -Depth, -Plot)
data_chem_new

dim(metadata_fungi_Root)
dim(data_chem_new)
rownames(data_chem_new) <- data_chem_new$Filt
data_chem_new <-
  data_chem_new %>%  select(-Filt)

metadata_fungi_Root_new <-
  left_join(
    select(metadata_fungi_Root, CAP1, CAP2, Filt, Crop, Depth, Origin), 
    rownames_to_column(data_chem_new, var="Filt"), by="Filt")
metadata_fungi_Root_new

metadata_fungi_Soil_new <-
left_join(
  select(metadata_fungi_Soil, CAP1, CAP2, Filt, Crop, Depth, Origin),
  rownames_to_column(data_chem_new, var="Filt"), by="Filt")
metadata_fungi_Soil_new

metadata_prok_Root_new <-
left_join(
  select(metadata_prok_Root, CAP1, CAP2, Filt, Crop, Depth, Origin),
  rownames_to_column(data_chem_new, var="Filt"), by="Filt")
metadata_prok_Root_new

metadata_prok_Soil_new <-
left_join(
  select(metadata_prok_Soil, CAP1, CAP2, Filt, Crop, Depth, Origin),
  rownames_to_column(data_chem_new, var="Filt"), by="Filt")
metadata_prok_Soil_new


set.seed(20)
GenerateEnvFit(cap_fungi_Root_full, data_chem_new) -> res_envfit_fungi_Root
res_envfit_fungi_Root
GenerateEnvFit(cap_fungi_Soil_pars, data_chem_new) -> res_envfit_fungi_Soil
res_envfit_fungi_Soil
GenerateEnvFit(cap_prok_Root_full, data_chem_new) -> res_envfit_prok_Root
res_envfit_prok_Root
GenerateEnvFit(cap_prok_Soil_full, data_chem_new) -> res_envfit_prok_Soil
res_envfit_prok_Soil


# plot ordination 
PlotCap <-function(dataframe, res_envfit, niche){
  ggplot(dataframe,
         aes(dataframe[,10],dataframe[,11], color = Depth, shape = Crop)) +
    geom_point(size = 1.5) +
    scale_colour_manual("Depth", values = paletteCB4) +
    theme_classic() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
    theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
    theme(legend.position="bottom")+
    guides(fill = "none",
           color = guide_legend(nrow = 2),
           shape = guide_legend(nrow = 2)) +
            grids(linetype = "dashed") -> plot_ord
  # to make sure the arrows are plotted above the point, set the points layer here.
  if (niche == "Root"){
    plot_ord + scale_shape_manual(name = "Crop", values = c(1,0,2)) -> plot_ord1
  } else {
    plot_ord +
      scale_shape_manual(name = "Crop", values = c(16,15,17)) +
      geom_point(size = 1.8) -> plot_ord1
  } 
  plot_ord1 + geom_segment(data = res_envfit, inherit.aes = FALSE,
                           mapping = aes(x = 0, xend = res_envfit[,1], y = 0, yend = res_envfit[,2]), 
                           color = "darkgrey",
                           arrow = arrow(length = unit(0.02, "npc"))) +
    geom_text(data = res_envfit,inherit.aes = FALSE,
              mapping = aes(x = 1.2 * res_envfit[,1], y = 1.2 * res_envfit[,2], label = Var),
              size = 3) -> plot_ord2
  return(plot_ord2)
}


PlotCap(as.data.frame(metadata_fungi_Root_new), res_envfit_fungi_Root, "Root") + 
  labs(title = "Fungi - Root") 

# Annotation position and hjust, vjust use to modify it
# xpos ypos        annotateText        hjust    vjust
# 1 -Inf -Inf Bottom Left (h0,v0)        0        0
# 2 -Inf  Inf    Top Left (h0,v1)        0        1
# 3  Inf -Inf  Bottom Right h1,v0        1        0
# 4  Inf  Inf     Top Right h1,v1        1        1


# *** FIGURE 3A - CAP models ---------------------------------------------------------------
ggarrange(PlotCap(as.data.frame(metadata_fungi_Root), res_envfit_fungi_Root, "Root") + 
            labs(title = "Fungi - Root", x="CAP1", y="CAP2"),
          PlotCap(as.data.frame(metadata_fungi_Soil), res_envfit_fungi_Soil, "Soil") +
            labs(title = "Fungi - Soil", x="CAP1", y="CAP2"),
          PlotCap(as.data.frame(metadata_prok_Root), res_envfit_prok_Root, "Root") + 
            labs(title = "Prokaryotes - Root", x="CAP1", y="CAP2"),
          PlotCap(as.data.frame(metadata_prok_Soil), res_envfit_prok_Soil, "Soil") +
            labs(title = "Prokaryotes - Soil", x="CAP1", y="CAP2"),
          labels = c("A", "B", "C", "D"),
          align = "hv",
          ncol = 4, 
          nrow = 1, 
          common.legend = TRUE, 
          legend = c("none")) -> plot_cap

plot_cap

# ********************************************--------------------------------------------
# ANALYSIS of SOIL DATA ------------------------------------------------------------------
head(metadata_fungi_Root)
str(metadata_fungi_Root)

metadata_fungi_Root$ph <- as.numeric(as.character(metadata_fungi_Root$ph))
metadata_fungi_Root$P <- as.numeric(as.character(metadata_fungi_Root$P))
metadata_fungi_Root$K <- as.numeric(as.character(metadata_fungi_Root$K))
metadata_fungi_Root$Mg <- as.numeric(as.character(metadata_fungi_Root$Mg))
metadata_fungi_Root$Ca <- as.numeric(as.character(metadata_fungi_Root$Ca))
metadata_fungi_Root$CEC <- as.numeric(as.character(metadata_fungi_Root$CEC))

# Add soil structure to the dataset 
soil_carb <-
  left_join(
    read.table("chemistry/soil_carbon.csv", header = TRUE, sep=",") %>% 
      mutate(Merge = paste(Crop, Plot, sep = "_"),
             Filt= paste(Merge, Depth, sep = "_")) %>% 
      select(Filt, Carbon) %>% 
      nest(data= -Filt) %>% 
      mutate(C = map(.x=data, ~sum(x=.), data = .x)) %>% 
      unnest(C) %>% 
      select(Filt, C), 
    read.table("chemistry/soil_carbon.csv", header = TRUE, sep=",") %>% 
      drop_na() %>% 
      mutate(Merge = paste(Crop, Plot, sep = "_"),
             Filt= paste(Merge, Depth, sep = "_")) %>% 
      select(Filt, Nitrogen) %>% 
      nest(data= -Filt) %>% 
      mutate(N = map(.x=data, ~sum(x=.), data = .x)) %>% 
      unnest(N) %>% 
      select(Filt, N), 
    by="Filt")

soil_carb 

soil_carb_new <-
  soil_carb %>% 
  separate(Filt, c("Crop", "Plot", "Depth"), sep = "_") 
soil_carb_new
write.csv(soil_carb_new, "chemistry/soil_carb.csv")


soil_str <-
  read.table("chemistry/soil_structure.csv", header = TRUE, sep=",") %>% 
  mutate(Merge = paste(Crop, Plot, sep = "_"),
         Filt = paste(Merge, Depth, sep = "_")) %>% 
  select(Filt, Sand, Silt, Clay) %>%  as_tibble()

soil_str

data_chem <-
  metadata_fungi_Root %>% 
  select(ph, P, K, Mg, Ca, CEC, Filt, Crop, Depth) %>% 
  left_join(
    left_join(soil_carb, 
              soil_str, by = "Filt"), 
    by = "Filt") %>% 
  as_tibble()

data_chem
as.data.frame(data_chem)

data_chem$Mg
metadata_fungi_mSeq$Mg
data_chem$C
range(data_chem$C)
data_chem$N
range(data_chem$N)

new_chem_data <-
  read.table("chemistry/new_chem_data.csv", header = TRUE, sep=",") 
new_chem_data

rm(data_chem)
data_chem <- new_chem_data

ExtractLetters <- function(dataframe, formula, Var){
  require(multcompView)
  compare_means(formula, data = subset(dataframe, Crop==Var), 
                p.adjust.method = "BH") -> test_CC 
  test_CC <- as.data.frame(test_CC)[,c(2,3,5)]
  colnames(test_CC) <- c("group1", "group2", "p.adj")
  test_CC
  test_CC2 <- data.frame(test_CC[,2], test_CC[,1], test_CC[,3])
  colnames(test_CC2) <- c("group1", "group2", "p.adj") # change p to p.adj
  rbind(test_CC, test_CC2) -> test_all
  as.dist(xtabs(test_all[, 3] ~ (test_all[, 2] + test_all[, 1])), diag = TRUE) -> dist_CC
  data.frame(multcompLetters(dist_CC)['Letters']) -> res_CC
  res_CC$sample <- rownames(res_CC)
  res_CC$Crop <- rep(Var, times=4)
  return(res_CC)
}


ExtractLetters(data_chem, formula(pH ~ Depth), "Poplar")

CalcType <- function(df, formula) {
  ExtractLetters(df, formula, "Poplar") -> pval_pop
  ExtractLetters(df, formula, "Prairie") -> pval_pra
  ExtractLetters(df, formula, "Switchgrass") -> pval_sw
  my_lables <- c(
    as.character(pval_pop$Letters),
    as.character(pval_pra$Letters),
    as.character(pval_sw$Letters))
  return(my_lables)
}


CalcType(data_chem, formula(pH ~ Depth)) -> ph_label
CalcType(data_chem, formula(P ~ Depth)) -> P_label
CalcType(data_chem, formula(K ~ Depth)) -> K_label
CalcType(data_chem, formula(Ca ~ Depth)) -> Ca_label
CalcType(data_chem, formula(Mg ~ Depth)) -> Mg_label
CalcType(data_chem, formula(CEC ~ Depth)) -> CEC_label

CalcType(data_chem, formula(Sand ~ Depth)) -> Sand_label
CalcType(data_chem, formula(Silt ~ Depth)) -> Silt_label
CalcType(data_chem, formula(Clay ~ Depth)) -> Clay_label
CalcType(data_chem, formula(C ~ Depth)) -> C_label
CalcType(data_chem, formula(N ~ Depth)) -> N_label

# plotting richness x Depth ----------------------------------------------------------------
PlotRich <- function(dataframe, my_labels, Var, y_max){
  #max(dataframe[,Var] + 0.001* max(dataframe[,Var])) -> y_max
  # plot
  ggline(dataframe, "Depth", Var, shape = NA, color = "Crop", add = "mean_se",) +
    stat_summary(geom = 'text', label = my_labels, fun= max, aes(y = y_max), size=2.5, color="black") +
    facet_grid(Crop ~ ., switch = "y") +
    theme_classic() +
    #expand_limits(y = 0) +
    #scale_shape_manual(name = "Crop", values = c(21,22,24)) +
    theme(strip.text.x = element_text(size = 7.5), strip.text.y = element_text(size = 7.5)) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 8)) +
    theme(axis.text.x = element_text(angle = 45, size = 7, hjust = 1.1, vjust = 1.15)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) + 
    grids(linetype = "dashed") +
    theme(legend.position="none") -> rich_plot 
  return(rich_plot)
}

# Removed 12 rows containing missing values (geom_point). 
# This is because I removed the points.
PlotRich(data_chem, ph_label, "pH", 7.2) + 
  labs(title = "pH", x=NULL, y="Richness") +
  scale_color_manual(values=c("black","black","black"))


# **** FIGURE S3A - Chemistry Wilcox tests--------------------------------------------------
ggarrange(
  PlotRich(data_chem, ph_label, "pH", 7.5) + 
    labs(title = "pH", x=NULL, y="pH") +
    scale_color_manual(values=c("black","black","black")),
 
   PlotRich(data_chem, P_label, "P", 55) + 
    labs(title = "P", x=NULL, y="P(ppm)") +
    scale_color_manual(values=c("black","black","black")),
 
   PlotRich(data_chem, K_label, "K", 225) + 
    labs(title = "K", x=NULL, y="K(ppm)") +
    scale_color_manual(values=c("black","black","black")),
 
   PlotRich(data_chem, Ca_label, "Ca", 1300) + 
    labs(title = "Ca", x=NULL, y="Ca(ppm)") +
    scale_color_manual(values=c("black","black","black")),
 
   PlotRich(data_chem, Mg_label, "Mg", 200) + 
    labs(title = "Mg", x=NULL, y="Mg(ppm)") +
    scale_color_manual(values=c("black","black","black")),
 
   PlotRich(data_chem, CEC_label, "CEC", 8.5) + 
    labs(title = "CEC", x=NULL, y="CEC (meq/L)") + 
    scale_color_manual(values=c("black","black","black")),
  
  PlotRich(data_chem, C_label, "C", 6) + 
    labs(title = "C", x=NULL, y="C (%)") + 
    scale_color_manual(values=c("black","black","black")),
 
   PlotRich(data_chem, N_label, "N", 0.6) + 
    labs(title = "N", x=NULL, y="N (%)") + 
    scale_color_manual(values=c("black","black","black")),
  
  PlotRich(data_chem, Sand_label, "Sand", 99) + 
    labs(title = "Sand", x=NULL, y="Sand (%)") + 
    scale_color_manual(values=c("black","black","black")),
 
   PlotRich(data_chem, Silt_label, "Silt", 50) + 
    labs(title = "Silt", x=NULL, y="Silt (%)") + 
    scale_color_manual(values=c("black","black","black")),
 
   PlotRich(data_chem, Clay_label, "Clay", 50) + 
    labs(title = "Clay", x=NULL, y="Clay (%)") + 
    scale_color_manual(values=c("black","black","black")),
  labels = c("A","","","","","","","","","",""),
  heights =  c(1,1,1,1,1,1,1,1,1,1,1),
  ncol = 4, 
  nrow = 3) -> Fig_SX_chem

Fig_SX_chem

title_chem = text_grob("Soil chemistry", size = 12, face = 2)
Fig_SX_chem <-
  grid.arrange(Fig_SX_chem, top = title_chem)


# Corrrelation between Env variables -------------------------------------------------------
library(corrplot)

names(data_chem)[names(data_chem) == "ph"] <- "pH"

cor_chem <-
  cor(data_chem[,4:14], method = "spearman")

# *** FIGURE S3B - correlation plot --------------------------------------------------------
corrplot.mixed(cor_chem, upper = "ellipse",
               title="Correlation plot", mar=c(0,0,1,0)) -> corr_plot
corr_plot

# # RDA models using Hellinger transformed data ----------------------------------------------
# MergeData(physeq_prok_filt, "Soil") -> physeq_prok_raw_Soil
# otu_table(physeq_prok_raw_Soil)[1:10, 1:10]
# 
# otu_prok_raw_Soil = as(otu_table(physeq_prok_raw_Soil), "matrix")
# otu_prok_raw_Soil = as.data.frame(otu_prok_raw_Soil)
# metadata_prok_raw_Soil <- as.data.frame(as.matrix(sample_data(physeq_prok_raw_Soil)))
# otu_prok_raw_Soil[1:10, 1:10]
# 
# otu_prok_raw_Soil_hel <- decostand(otu_prok_raw_Soil, method = "hellinger")
# otu_prok_raw_Soil_hel[1:10, 1:10]
# 
# rda_prok_Soil_full <- vegan::rda(otu_prok_raw_Soil_hel ~ Crop + Depth + Crop:Depth, metadata_prok_raw_Soil)
# rda_prok_Soil_full
# anova(rda_prok_Soil_full, permutations = 9999, by = "terms")
# 
# pars_model_prok_Soil_hel <- vegan::ordistep(rda_prok_Soil_full, direction = "backward", Pout = 0.05, permutations = how(nperm = 9999))
# pars_model_prok_Soil_hel <- formula(pars_model_prok_Soil_hel)[1:3]
# pars_model_prok_Soil_hel
# 
# # re-run with the most parsimonious model
# rda_prok_Soil_pars <- vegan::rda(pars_model_prok_Soil_hel, metadata_prok_raw_Soil)
# rda_prok_Soil_pars
# anova(rda_prok_Soil_pars, permutations = 9999, by = "terms")
# 
# # Now extracting variances for each term
# rda_prok_Soil_Crop <-rda(otu_prok_raw_Soil_hel ~ Crop, metadata_prok_raw_Soil)
# RsquareAdj(rda_prok_Soil_Crop)
# rda_prok_Soil_Depth <-rda(otu_prok_raw_Soil_hel ~ Depth, metadata_prok_raw_Soil)
# RsquareAdj(rda_prok_Soil_Depth)
# 
# RsquareAdj(rda_prok_Soil_full)
# RsquareAdj(rda_prok_Soil_pars)
# 
# sum(unlist(RsquareAdj(rda_prok_Soil_Crop)[1]),unlist(RsquareAdj(rda_prok_Soil_Depth)[1]))
# sum(unlist(RsquareAdj(rda_prok_Soil_Crop)[2]),unlist(RsquareAdj(rda_prok_Soil_Depth)[2]))
# 
# # Testing for dispersion
# anova(betadisper(vegdist(otu_prok_raw_Soil_hel, method = "euclidean"), metadata_prok_raw_Soil$Depth), permutations = 9999)
# anova(betadisper(vegdist(otu_prok_raw_Soil_hel, method = "euclidean"), metadata_prok_raw_Soil$Crop), permutations = 9999)
# 
# # prepare for plotting
# identical(rownames(scores(rda_prok_Soil_pars)$sites), rownames(metadata_prok_raw_Soil))
# metadata_prok_raw_Soil <- cbind(metadata_prok_raw_Soil, scores(rda_prok_Soil_pars)$sites)
# head(metadata_prok_raw_Soil)
# 
# # protest ----------------------------------------------------------------------------------
# protest(
#   cap_prok_Soil_full,
#   rda_prok_Soil_pars,
#   scores = "sites",
#   permutations = how(nperm = 9999))

# ********************************************--------------------------------------------
# MULTIVARIATE DISPERSION ------------------------------------------------------------------
# Test the homogenity of group variances
p.adjust(anova(betadisper(vegdist(otu_fungi_Root, method="bray"), 
                          metadata_fungi_Root$Crop), permutations = 9999)$`Pr(>F)`, "bonferroni") #0.5192738
boxplot(betadisper(vegdist(otu_fungi_Root, method="bray"), metadata_fungi_Root$Crop))
p.adjust(anova(betadisper(vegdist(otu_fungi_Root, method="bray"),
                          metadata_fungi_Root$Depth), permutations = 9999)$`Pr(>F)`, "bonferroni") #0.01685156 **
boxplot(betadisper(vegdist(otu_fungi_Root, method="bray"), metadata_fungi_Root$Depth))

p.adjust(anova(betadisper(vegdist(otu_fungi_Soil, method="bray"), 
                          metadata_fungi_Soil$Crop), permutations = 9999)$`Pr(>F)`, "bonferroni") #0.09345761 .
boxplot(betadisper(vegdist(otu_fungi_Soil, method="bray"), metadata_fungi_Soil$Crop))
p.adjust(anova(betadisper(vegdist(otu_fungi_Soil, method="bray"),
                          metadata_fungi_Soil$Depth), permutations = 9999)$`Pr(>F)`, "bonferroni") #0.0000025 ***
boxplot(betadisper(vegdist(otu_fungi_Soil, method="bray"), metadata_fungi_Soil$Depth))

p.adjust(anova(betadisper(vegdist(otu_prok_Root, method="bray"), 
                          metadata_prok_Root$Crop), permutations = 9999)$`Pr(>F)`, "bonferroni") #0.7285613
boxplot(betadisper(vegdist(otu_prok_Root, method="bray"), metadata_prok_Root$Crop))
p.adjust(anova(betadisper(vegdist(otu_prok_Root, method="bray"),
                          metadata_prok_Root$Depth), permutations = 9999)$`Pr(>F)`, "bonferroni") #0.2515306
boxplot(betadisper(vegdist(otu_prok_Root, method="bray"), metadata_prok_Root$Depth))

p.adjust(anova(betadisper(vegdist(otu_prok_Soil, method="bray"), 
                          metadata_prok_Soil$Crop), permutations = 9999)$`Pr(>F)`, "bonferroni") #0.2428832 
boxplot(betadisper(vegdist(otu_prok_Soil, method="bray"), metadata_prok_Soil$Crop))
p.adjust(anova(betadisper(vegdist(otu_prok_Soil, method="bray"),
                          metadata_prok_Soil$Depth), permutations = 9999)$`Pr(>F)`, "bonferroni") #0.0000000027 ***
boxplot(betadisper(vegdist(otu_prok_Soil, method="bray"), metadata_prok_Soil$Depth))

metadata_fungi_Root$Origin <- rep("Root", nrow(metadata_fungi_Root))
metadata_prok_Root$Origin <- rep("Root", nrow(metadata_prok_Root))
metadata_fungi_Soil$Origin <- rep("Soil", nrow(metadata_fungi_Soil))
metadata_prok_Soil$Origin <- rep("Soil", nrow(metadata_prok_Soil))


ExtractLabel <- function(dataframe, metadata, Var) {
  require(multcompView)
  metadata$Depth <- gsub("-", "_", metadata_fungi_Root$Depth)
  permdisp <-
    betadisper(vegdist(dataframe, method = "bray"), metadata[, Var])
  dist <-
    permutest(permdisp, permutations = 9999, pairwise = T)
  label_disp <-
    data.frame(multcompLetters(p.adjust(dist$pairwise$observed,
                                        method = "bonferroni"))['Letters'])
  return(label_disp)
}


ExtractLabel(otu_fungi_Root, metadata_fungi_Root, "Crop")

# Plotting Betadisper for CAP - soil and root separated 
PlotBetadisperCap <- function(betadisp,  my_labels, metadata, var, niche){
  # creating a label and dataframe
  max(betadisp$distances + 0.1* max(betadisp$distances)) -> labels_y
  data.frame(betadisp$distances) -> df
  colnames(df) <- c("distance")
  if (identical(rownames(df), rownames(metadata))==TRUE){
    df$Crop <- metadata$Crop
    df$Origin <- metadata$Origin
    df$Depth <- metadata$Depth
    # plotting
    ggplot(df, aes(x=get(var), y=distance)) +
      geom_jitter(aes(x=get(var), y=distance, shape = Crop, color = Depth),size = 0.9, alpha = 0.8) +
      geom_boxplot(outlier.colour="black", outlier.shape = 8, outlier.size = 0.8, alpha=0.6, lwd = 0.3) +
      stat_summary(fun=mean, geom="point", shape=18, size=1.9, color="red", fill="red") +
      stat_summary(geom = 'text', label = my_labels, fun= max, aes(y = labels_y), size=3, color="black") +
      scale_colour_manual("Depth", values = paletteCB4) +
      theme_classic() +
      theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) +
      theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) +
      theme(axis.text.x = element_text(angle = 45, size = 7,hjust = 1, vjust = 1.1)) +
      theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
      theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
      theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
      theme(legend.position="bottom")+
      guides(fill = "none",
             color = guide_legend(nrow = 2),
             shape = guide_legend(nrow = 2)) +
      grids(linetype = "dashed") -> plot_ord
    if (niche == "Root"){
      plot_ord + scale_shape_manual(name = "Crop", values = c(1,0,2)) -> plot_ord1
    } else {
      plot_ord +
        scale_shape_manual(name = "Crop", values = c(16,15,17)) -> plot_ord1
    } 
    return(plot_ord1)
  }else{
    stop("Error: dataframe and metadata rownames are not matching!")
  }
}


PlotBetadisperCap(betadisp = betadisper(vegdist(otu_fungi_Root, method="bray"), metadata_fungi_Root$Crop),
                  my_labels = as.character(ExtractLabel(otu_fungi_Root, metadata_fungi_Root, "Crop")$Letters),
                  metadata =metadata_fungi_Root, 
                  var = "Crop",
                  niche="Root") 


# *** FIGURE 3B - betasiper for CAP --------------------------------------------------------
ggarrange(PlotBetadisperCap(betadisp = betadisper(vegdist(otu_fungi_Root, method="bray"), metadata_fungi_Root$Crop),
                            my_labels = "",
                            metadata =metadata_fungi_Root, var = "Crop",niche="Root") +
            labs(title = "Crop", y ="Distance to centroid", x= NULL),
          PlotBetadisperCap(betadisp = betadisper(vegdist(otu_fungi_Root, method="bray"), metadata_fungi_Root$Depth),
                            my_labels = as.character(ExtractLabel(otu_fungi_Root, metadata_fungi_Root, "Depth")$Letters),
                            metadata =metadata_fungi_Root, var = "Depth",niche="Root") +
            labs(title = "Depth", y ="Distance to centroid", x= NULL),
          PlotBetadisperCap(betadisp = betadisper(vegdist(otu_fungi_Soil, method="bray"), metadata_fungi_Soil$Crop),
                            my_labels = "",
                            metadata =metadata_fungi_Soil, var = "Crop",niche="Soil") +
            labs(title = "Crop", y ="Distance to centroid", x= NULL),
          PlotBetadisperCap(betadisp = betadisper(vegdist(otu_fungi_Soil, method="bray"), metadata_fungi_Soil$Depth),
                            my_labels = as.character(ExtractLabel(otu_fungi_Soil, metadata_fungi_Soil, "Depth")$Letters),
                            metadata =metadata_fungi_Soil, var = "Depth",niche="Soil") +
            labs(title = "Depth", y ="Distance to centroid", x= NULL),
          PlotBetadisperCap(betadisp = betadisper(vegdist(otu_prok_Root, method="bray"), metadata_prok_Root$Crop),
                            my_labels = "",
                            metadata =metadata_prok_Root, var = "Crop",niche="Root") +
            labs(title = "Crop", y ="Distance to centroid", x= NULL),
          PlotBetadisperCap(betadisp = betadisper(vegdist(otu_prok_Root, method="bray"), metadata_prok_Root$Depth),
                            my_labels = "",
                            metadata =metadata_prok_Root, var = "Depth",niche="Root") +
            labs(title = "Depth", y ="Distance to centroid", x= NULL),
          PlotBetadisperCap(betadisp = betadisper(vegdist(otu_prok_Soil, method="bray"), metadata_prok_Soil$Crop),
                            my_labels = "",
                            metadata =metadata_prok_Soil, var = "Crop",niche="Soil") +
            labs(title = "Crop", y ="Distance to centroid", x= NULL),
          PlotBetadisperCap(betadisp = betadisper(vegdist(otu_prok_Soil, method="bray"), metadata_prok_Soil$Depth),
                            my_labels = as.character(ExtractLabel(otu_prok_Soil, metadata_prok_Soil, "Depth")$Letters),
                            metadata =metadata_prok_Soil, var = "Depth",niche="Soil") +
            labs(title = "Depth", y ="Distance to centroid", x= NULL),
          widths = c(1,1,1,1,1,1,1,1),
          labels = c("E","","F","","G","","H",""),
          align = "hv",
          ncol = 8, 
          nrow = 1,
          common.legend = TRUE, 
          legend = c("none")) -> betadisp_cap_plot

betadisp_cap_plot

# *** FIGURE 3 - complete ------------------------------------------------------------------
ggarrange(
  ggarrange(
    plot_cap,
    betadisp_cap_plot,
    ncol = 1,
    nrow = 2),
  as_ggplot(legend_beta_div),
  nrow = 2, 
  heights = c(1, 0.06)) -> Fig_3_cap

Fig_3_cap


# # *** FIGURE S3 - RDA models 
# ggarrange(PlotCap(metadata_fungi_raw_Root, res_envfit_rda_fungi_Root, "Root") + 
#             labs(title = "Fungi - Root", x="RDA1", y="RDA2"),
#           PlotCap(metadata_fungi_raw_Soil, res_envfit_rda_fungi_Soil, "Soil") +
#             labs(title = "Fungi - Soil", x="RDA1", y="RDA2"),
#           PlotCap(metadata_prok_raw_Root, res_envfit_rda_prok_Root, "Root") + 
#             labs(title = "Prokaryotes - Root", x="RDA1", y="RDA2"),
#           PlotCap(metadata_prok_raw_Soil, res_envfit_prok_Soil, "Soil") +
#             labs(title = "Prokaryotes - Soil", x="RDA1", y="RDA2"),
#           labels = c("A", "B", "C", "D"),
#           align = "hv",
#           ncol = 2, 
#           nrow = 2, 
#           common.legend = TRUE, 
#           legend = c("none")) -> plot_rda
# 
# plot_rda

# ********************************************----------------------------------------------
# BETA-DIVERSITY group - by - group --------------------------------------------------------
BetaGroup <- function(physeq, Var){
  if (Var == "Poplar" | Var == "Prairie" | Var == "Switchgrass") {
    sub_set <- subset(sample_data(physeq), Crop == Var)
    physeq_filt <- merge_phyloseq(otu_table(physeq),
                                  tax_table(physeq),
                                  refseq(physeq),
                                  sub_set)
    otu_table(physeq_filt) <-
      otu_table(physeq_filt)[which(rowSums(otu_table(physeq_filt)) > 0),]
    otu <- as(t(otu_table(physeq_filt)), "matrix")
    otu <- as.data.frame(otu)
    metadata <- as.data.frame(as.matrix(sample_data(physeq_filt)))
    table(metadata$Section) %T>% print()
    model <- vegan::capscale(otu ~ Section + Condition(Plot), metadata ,distance = "bray")
    pval <- anova(model, permutations = 9999)
    adjr <- RsquareAdj(model)
    anosim <- anosim(vegdist(otu, method="bray"), metadata$Section, permutations = 9999, distance = "bray")
    betadisp <- anova(betadisper(vegdist(otu, method="bray"), metadata$Section), permutations = 9999)
  }else{
    sub_set <- subset(sample_data(physeq), Section == Var)
    physeq_filt <- merge_phyloseq(otu_table(physeq),
                                  tax_table(physeq),
                                  refseq(physeq),
                                  sub_set)
    otu_table(physeq_filt) <-
      otu_table(physeq_filt)[which(rowSums(otu_table(physeq_filt)) > 0),]
    otu <- as(t(otu_table(physeq_filt)), "matrix")
    otu <- as.data.frame(otu)
    metadata <- as.data.frame(as.matrix(sample_data(physeq_filt)))
    table(metadata$Crop) %T>% print()
    model <- vegan::capscale(otu ~ Crop + Condition(Plot), metadata ,distance = "bray")
    pval <- anova(model, permutations = 9999)
    adjr <- RsquareAdj(model)
    anosim <- anosim(vegdist(otu, method="bray"), metadata$Crop, permutations = 9999, distance = "bray")
    betadisp <- anova(betadisper(vegdist(otu, method="bray"), metadata$Crop), permutations = 9999)
  }
  return(list(model, pval, adjr, anosim, betadisp))
}

BetaGroup(physeq_fungi_mSeq_Root, "0-10cm") 
BetaGroup(physeq_fungi_mSeq_Root, "10-25cm")
BetaGroup(physeq_fungi_mSeq_Root, "25-50cm")
BetaGroup(physeq_fungi_mSeq_Root, "50-100cm")

# Check Marti Anderson 2017 _ PERMANOVA
# In ANOSIM. An R value close to "1.0" suggests dissimilarity between groups while 
#an R value close to "0" suggests an even distribution of high and low ranks within
#and between groups. R values below "0" suggest that dissimilarities are greater 
#within groups than between groups.

test_df[[2]]$`Pr(>F)`[1]
test_df[[3]]$adj.r.squared
test_df[[4]]$`Pr(>F)`[1]

anosim(dist_fungi, metadata_fungi$Origin, permutations = 999, distance = "bray")

# dataframe for Crop effects ---------------------------------------------------------------
data.frame(Crop=c("Poplar", "Prairie", "Switchgrass"),
           Origin=rep(c("Root", "Soil"), each=3),
           Organism=rep(c("Fungi", "Prokaryotes"), each=6),
           #Color=rep(c("#9A7145","#57C4AD","#4B2D0B","#006164"), each=3),
           Color=rep(c("black","grey","black","grey"), each=3),
           Adj.R2=c(BetaGroup(physeq_fungi_mSeq_Root, "Poplar")[[3]]$adj.r.squared,
                    BetaGroup(physeq_fungi_mSeq_Root, "Prairie")[[3]]$adj.r.squared,
                    BetaGroup(physeq_fungi_mSeq_Root, "Switchgrass")[[3]]$adj.r.squared,
                    BetaGroup(physeq_fungi_mSeq_Soil, "Poplar")[[3]]$adj.r.squared,
                    BetaGroup(physeq_fungi_mSeq_Soil, "Prairie")[[3]]$adj.r.squared,
                    BetaGroup(physeq_fungi_mSeq_Soil, "Switchgrass")[[3]]$adj.r.squared,
                    BetaGroup(physeq_prok_mSeq_Root, "Poplar")[[3]]$adj.r.squared,
                    BetaGroup(physeq_prok_mSeq_Root, "Prairie")[[3]]$adj.r.squared,
                    BetaGroup(physeq_prok_mSeq_Root, "Switchgrass")[[3]]$adj.r.squared,
                    BetaGroup(physeq_prok_mSeq_Soil, "Poplar")[[3]]$adj.r.squared,
                    BetaGroup(physeq_prok_mSeq_Soil, "Prairie")[[3]]$adj.r.squared,
                    BetaGroup(physeq_prok_mSeq_Soil, "Switchgrass")[[3]]$adj.r.squared)) -> AdjR2_depth
AdjR2_depth$Line <- paste(AdjR2_depth$Origin, AdjR2_depth$Organism)
AdjR2_depth$Crop <- factor(AdjR2_depth$Crop,levels=c("Switchgrass", "Prairie", "Poplar"))
AdjR2_depth$Line <- factor(AdjR2_depth$Line,levels=c("Soil Fungi","Root Fungi", "Soil Prokaryotes", "Root Prokaryotes"))
#AdjR2_depth$Color <- factor(AdjR2_depth$Color,levels=c("#9A7145","#57C4AD","#4B2D0B","#006164"))
AdjR2_depth$Color <- factor(AdjR2_depth$Color,levels=c("black","grey"))
AdjR2_depth

# dataframe for Depth effects --------------------------------------------------------------
data.frame(Depth=c("0-10cm","10-25cm", "25-50cm","50-100cm"),
           Origin=rep(c("Root", "Soil"), each=4),
           Organism=rep(c("Fungi", "Prokaryotes"), each=8),
           #Color=rep(c("#9A7145","#57C4AD","#4B2D0B","#006164"), each=4),
           Color=rep(c("black","grey","black","grey"), each=4),
           Adj.R2=c(BetaGroup(physeq_fungi_mSeq_Root, "0-10cm")[[3]]$adj.r.squared,
                    BetaGroup(physeq_fungi_mSeq_Root, "10-25cm")[[3]]$adj.r.squared,
                    BetaGroup(physeq_fungi_mSeq_Root, "25-50cm")[[3]]$adj.r.squared,
                    BetaGroup(physeq_fungi_mSeq_Root, "50-100cm")[[3]]$adj.r.squared,
                    BetaGroup(physeq_fungi_mSeq_Soil, "0-10cm")[[3]]$adj.r.squared,
                    BetaGroup(physeq_fungi_mSeq_Soil, "10-25cm")[[3]]$adj.r.squared,
                    BetaGroup(physeq_fungi_mSeq_Soil, "25-50cm")[[3]]$adj.r.squared,
                    BetaGroup(physeq_fungi_mSeq_Soil, "50-100cm")[[3]]$adj.r.squared,
                    BetaGroup(physeq_prok_mSeq_Root, "0-10cm")[[3]]$adj.r.squared,
                    BetaGroup(physeq_prok_mSeq_Root, "10-25cm")[[3]]$adj.r.squared,
                    BetaGroup(physeq_prok_mSeq_Root, "25-50cm")[[3]]$adj.r.squared,
                    BetaGroup(physeq_prok_mSeq_Root, "50-100cm")[[3]]$adj.r.squared,
                    BetaGroup(physeq_prok_mSeq_Soil, "0-10cm")[[3]]$adj.r.squared,
                    BetaGroup(physeq_prok_mSeq_Soil, "10-25cm")[[3]]$adj.r.squared,
                    BetaGroup(physeq_prok_mSeq_Soil, "25-50cm")[[3]]$adj.r.squared,
                    BetaGroup(physeq_prok_mSeq_Soil, "50-100cm")[[3]]$adj.r.squared)) -> AdjR2_crop
AdjR2_crop$Line <- paste(AdjR2_crop$Origin, AdjR2_crop$Organism)
AdjR2_crop$Depth <- factor(AdjR2_crop$Depth,levels=c("50-100cm","25-50cm","10-25cm","0-10cm"))
AdjR2_crop$Line <- factor(AdjR2_crop$Line,levels=c("Soil Fungi","Root Fungi", "Soil Prokaryotes", "Root Prokaryotes"))
#AdjR2_crop$Color <- factor(AdjR2_crop$Color,levels=c("#9A7145","#57C4AD","#4B2D0B","#006164"))
AdjR2_crop$Color <- factor(AdjR2_crop$Color,levels=c("black","grey"))
AdjR2_crop

grid.table(AdjR2_crop)

# plotting AdjR2 lines -------------------------------------------------------------------- 
PlotAdjrLines <- function(dataframe, Var){
  Fig_R2_depth <-
    ggplot(dataframe, aes(x=get(Var), y=Adj.R2, group=Line)) +
    #labs(title=expression(paste(italic(Adj.R^{2}),"of Depth")), y=expression(italic(Adj.R^{2})), x="") +
    geom_line(aes(linetype=Organism, color=Line), size=0.7, show.legend = FALSE) +
    geom_point(aes(shape=Line, color=Line, fill=Line), size=2.5, stroke = 1) + 
    scale_color_manual(values=c("black","black", "grey", "grey")) +
    #scale_color_manual(values=c("black","grey", "black", "grey")) +
    theme_classic() +
    coord_flip() +
    ylim(0, 0.5) +
    #theme(panel.grid.minor = element_blank()) +
    #scale_y_continuous(breaks = axis_scale(0.1)) +
    theme(plot.title = element_text(size = 9, face="bold", hjust = 0.5)) +
    theme(legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 7)) +
    theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 1, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
    grids(linetype = "dashed", color = ) +
    # theme(legend.position="bottom") +
    guides(color=guide_legend(ncol=1))
  return(Fig_R2_depth)
}

PlotAdjrLines(AdjR2_depth, "Crop") + 
  scale_shape_manual(values = c(18,1,18,1)) +
  labs(title="Effect of depth across crops", y=expression(italic(Adj.R^{2})), x=NULL)

PlotAdjrLines(AdjR2_crop, "Depth") + 
  scale_shape_manual(values = c(18,1,18,1)) +
  labs(title="Effect of crop across depth", y=expression(italic(Adj.R^{2})), x=NULL) +
  theme(legend.position=c(0.7, 0.8))

# *** FIGURE 4 - Adj R2 lines --------------------------------------------------------------
line_bars_AdjR2 = ggarrange(PlotAdjrLines(AdjR2_depth, "Crop") + 
                              scale_shape_manual(values = c(18,1,18,1)) +
                              labs(title="Effect of depth across crops", y=expression(italic(Adj.R^{2})), x=NULL) +
                              theme(legend.position="none"),
                            PlotAdjrLines(AdjR2_crop, "Depth") + 
                              scale_shape_manual(values = c(18,1,18,1)) +
                              labs(title="Effect of crop across depths", y=expression(italic(Adj.R^{2})), x=NULL) +
                              theme(legend.position=c(0.75, 0.8)),
                            labels = c("A", "B"), 
                            heights =  c(1,1),
                            align = "hv", 
                            ncol = 2,
                            nrow = 1,
                            common.legend = FALSE)

line_bars_AdjR2

ggarrange(Fig_3_cap,
          line_bars_AdjR2,
          ncol = 2,
          nrow = 1, 
          heights = c(1, 1),
          widths = c(1.3, 0.4)) -> Fig_3_add

Fig_3_add

# ********************************************----------------------------------------------
# ALPHA DIVERSITY --------------------------------------------------------------------------

head(metadata_fungi_mSeq)
names(metadata_fungi_mSeq)[names(metadata_fungi_mSeq) == "Depth"] <- "Bottom"
names(metadata_fungi_mSeq)[names(metadata_fungi_mSeq) == "Section"] <- "Depth"

metadata_fungi_mSeq$richness <- specnumber(otu_fungi_mSeq, MARGIN = 1)
metadata_fungi_mSeq$shannon <- diversity(otu_fungi_mSeq,index = "shannon", MARGIN=1) 
metadata_fungi_mSeq$Eshannon <- 1 - metadata_fungi_mSeq$shannon/log(metadata_fungi_mSeq$richness)
metadata_fungi_mSeq$Bottom <- as.numeric(as.character(metadata_fungi_mSeq$Bottom))
metadata_fungi_mSeq$ReadNo <- sample_sums(physeq_fungi_filt)
head(metadata_fungi_mSeq)
str(metadata_fungi_mSeq)

head(metadata_prok_mSeq)
metadata_prok_mSeq$richness <- specnumber(otu_prok_mSeq, MARGIN = 1)
metadata_prok_mSeq$shannon <- diversity(otu_prok_mSeq,index = "shannon", MARGIN=1) 
metadata_prok_mSeq$Eshannon <- 1 - metadata_prok_mSeq$shannon/log(metadata_prok_mSeq$richness)
metadata_prok_mSeq$Bottom <- as.numeric(as.character(metadata_prok_mSeq$Bottom))
head(metadata_prok_mSeq)
str(metadata_prok_mSeq)

# Wilcox test ----------------------------------------------------------------------------------
# Calculating differences in alpha diversity - pairwise wilcox tests

CompSampl <- function(df, formula){
  require(multcompView)
  require(lazyeval)
  compare_means(formula, data = df, 
                method = "wilcox.test", p.adjust.method = "BH") -> test_CC
  print(test_CC)
  test_CC <- as.data.frame(test_CC)[,c(2,3,5)] # to change form p to p.adj do 4 to 5
  test_CC2 <- data.frame(test_CC[,2], test_CC[,1], test_CC[,3])
  colnames(test_CC2) <- c("group1", "group2", "p.adj") # change p to p.adj
  rbind(test_CC, test_CC2) -> test_all
  as.dist(xtabs(test_all[, 3] ~ (test_all[, 2] + test_all[, 1])), diag = TRUE) -> dist_CC
  #print(dist_CC)
  data.frame(multcompLetters(dist_CC)['Letters']) -> res_CC
  res_CC$Var <- rownames(res_CC)
  return(res_CC)
}


CompSampl(df = metadata_fungi_mSeq, formula(richness ~ Depth))


df_letters <-
metadata_fungi_mSeq %>% 
  select(Crop, Origin, Depth, richness) %>% 
  mutate(Group=paste(Origin, Crop, sep="_")) %>% 
  nest(-Group) %>% 
  mutate(Sig = map(.x=data, ~CompSampl(df = ., formula(richness ~ Depth)), data = .x)) %>% 
  unnest(Sig) %>% 
  select(Group, Letters, Var) %>%
  separate(Group, c("Origin", "Crop"), sep = "_") %>% 
  rename(Depth=Var) %>% 
  select(Crop, Origin, Depth, Letters) %>% 
  arrange(Crop, Depth, Origin)

df_letters

df_top <-
metadata_fungi_mSeq %>% 
  select(Crop, Origin, Depth, richness) %>% 
  mutate(Split=paste(Depth, Origin, sep="_"),
         Group=paste(Split, Crop, sep = "_")) %>%
  select(richness, Group) %>% 
  nest(-Group) %>%  
  mutate(Max = map(.x=data, ~max(x=.), data = .x)) %>% 
  unnest(Max) %>% 
  separate(Group, c("Depth", "Origin", "Crop"), sep = "_") %>% 
  select(Crop, Origin, Depth,  Max) %>%
  arrange(Crop, Depth, Origin)

df_top

df_join <-
  df_letters %>% 
  left_join(df_top) 

metadata_fungi_mSeq %>% 
  ggplot(aes(x = Depth, y=richness, color=Origin)) +
  geom_boxplot(outlier.colour="black", outlier.shape = 8, outlier.size = 0.8, alpha=0.6, lwd = 0.3) +
  geom_text(data = df_join, aes(x=Depth, y = 5 + Max, label = Letters), size=3,
            vjust=-0.3, position = position_dodge(width = 1.2), show.legend = F)+
  theme_classic() +
  scale_color_manual(values = c("red", "black")) +
  coord_cartesian(ylim = c(0, 650))+
  facet_wrap(~Crop) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.title = element_text(angle = 0, size = 7, face = "bold"),
        axis.text.x = element_text(angle = 45, size = 7,hjust = 1, vjust = 1.1),
        axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
        legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"),
        legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7),
        legend.position="right")

# Function to calculate Significant differences across groups - wilcox.test and p.adjust BH
getAlphaStats <- function(metadata, formula, metric){
  metric <- enquo(metric) 
letters <-
  metadata %>% 
  select(Crop, Origin, Depth, !!metric) %>% 
  mutate(Group=paste(Origin, Crop, sep="_")) %>% 
  nest(-Group) %>% 
  mutate(Sig = map(.x=data, ~CompSampl(df = ., formula), data = .x)) %>% 
  unnest(Sig) %>% 
  select(Group, Letters, Var) %>%
  separate(Group, c("Origin", "Crop"), sep = "_") %>% 
  rename(Depth=Var) %>% 
  select(Crop, Origin, Depth, Letters) %>% 
  arrange(Crop, Depth, Origin)

top <-
  metadata %>% 
  select(Crop, Origin, Depth, !!metric) %>% 
  mutate(Split=paste(Depth, Origin, sep="_"),
         Group=paste(Split, Crop, sep = "_")) %>%
  select(!!metric, Group) %>% 
  nest(-Group) %>%  
  mutate(Max = map(.x=data, ~max(x=.), data = .x)) %>% 
  unnest(Max) %>% 
  separate(Group, c("Depth", "Origin", "Crop"), sep = "_") %>% 
  select(Crop, Origin, Depth,  Max) %>%
  arrange(Crop, Depth, Origin)
new_df <-
  letters %>% 
  left_join(top) 
return(new_df)
}

stat_fungi_rich <-
  getAlphaStats(metadata_fungi_mSeq, formula(richness ~ Depth), "richness")

stat_fungi_Eshan <-
  getAlphaStats(metadata_fungi_mSeq, formula(Eshannon ~ Depth), "Eshannon")

stat_prok_rich <-
  getAlphaStats(metadata_prok_mSeq, formula(richness ~ Depth), "richness")

stat_prok_Eshan <-
  getAlphaStats(metadata_prok_mSeq, formula(Eshannon ~ Depth), "Eshannon")


# metadata- dataframe with richness data, metric= richness or shannon, stat_df= the dataframe with the significant levels,
# y_max= the y axis where the letters will be plotted, const= a constant added to y_max to plot the letters a little higher
PlotAlphaDiv <-function(metadata, metric, stat_df, y_max, const){
  metadata %>% 
  ggplot(aes(x = Depth, y=get(metric), color=Origin)) +
  geom_boxplot(outlier.colour="black", outlier.shape = 8, outlier.size = 0.8, 
               alpha=0.6, lwd = 0.3, show.legend = T) +
  geom_text(data = stat_df, aes(x=Depth, y = Max + const, label = Letters), size=3,
            vjust=-0.3, position = position_dodge(width = 1.2), show.legend = F) +
    guides(color = guide_legend("Niche", override.aes = list(shape = 15))) +
    theme_classic() +
    scale_color_manual(values = c("red", "black")) +
    coord_cartesian(ylim = c(0, y_max))+
    facet_wrap(~Crop) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.title = element_text(angle = 0, size = 7, face = "bold"),
        axis.text.x = element_text(angle = 45, size = 7,hjust = 1, vjust = 1.1),
        axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
        legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"),
        legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7),
        legend.position="right") -> plot_alpha_div
return(plot_alpha_div)
}

PlotAlphaDiv(metadata_fungi_mSeq, "richness", stat_fungi_rich, 650, 5) +
  labs(title = "Fungi", x="Depth", y="OTU richness")

PlotAlphaDiv(metadata_fungi_mSeq, "Eshannon", stat_fungi_Eshan, 1, 0.05)+
  labs(title = "Fungi", x="Depth", y="OTU richness")

PlotAlphaDiv(metadata_prok_mSeq, "richness", stat_prok_rich, 6000, 20)+
  labs(title = "Fungi", x="Depth", y="OTU richness")

PlotAlphaDiv(metadata_prok_mSeq, "Eshannon", stat_prok_Eshan, 1, 0.05)+
  labs(title = "Fungi", x="Depth", y="OTU richness")


ggarrange(PlotAlphaDiv(metadata_fungi_mSeq, "richness", stat_fungi_rich, 650, 5) +
            labs(title = "Fungi", x=NULL, y="OTU richness"),
          PlotAlphaDiv(metadata_prok_mSeq, "richness", stat_prok_rich, 6000, 20)+
            labs(title = "Prokaryotes", x=NULL, y=NULL),
          PlotAlphaDiv(metadata_fungi_mSeq, "Eshannon", stat_fungi_Eshan, 1, 0.05)+
            labs(x="Depth", y="OTU Shannon"),
          PlotAlphaDiv(metadata_prok_mSeq, "Eshannon", stat_prok_Eshan, 1, 0.05)+
            labs(x="Depth", y=NULL),
            labels = c("A", "B", "C","D"),
          align = "hv",
          widths = c(1,1,1,1), 
          ncol = 2, 
          nrow = 2,
          common.legend = TRUE, 
          legend = c("bottom")) -> Fig1_rich

Fig1_rich

# Significant differences in Alpha div -----------------------------------------------------
library(car)
library(tidyverse)
library(broom)
library(lme4)
library(MASS)

head(metadata_fungi_mSeq)
count(metadata_fungi_mSeq, Crop)
dim(subset(metadata_fungi_mSeq, Crop=="Poplar"))

lm(subset(metadata_fungi_mSeq, Crop=="Poplar") ~ sqrt(ReadNo) + richness)

# Richness fungi
kruskal.test(richness ~ Origin, subset(metadata_fungi_mSeq, Crop=="Poplar"))
kruskal.test(richness ~ Depth, subset(metadata_fungi_mSeq, Crop=="Poplar"))
summary(aov(richness ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Poplar")))
par(mfcol=c(2,2))
plot(aov(richness ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Poplar")))
shapiro.test(residuals(aov(richness ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Poplar"))))
leveneTest(richness ~ Origin * Depth, data=subset(metadata_fungi_mSeq, Crop=="Poplar"))

kruskal.test(richness ~ Origin, subset(metadata_fungi_mSeq, Crop=="Prairie"))
kruskal.test(richness ~ Depth, subset(metadata_fungi_mSeq, Crop=="Prairie"))
summary(aov(richness ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Prairie")))
par(mfcol=c(2,2))
plot(aov(richness ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Prairie")))
shapiro.test(residuals(aov(richness ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Prairie"))))
leveneTest(richness ~ Origin * Depth, data=subset(metadata_fungi_mSeq, Crop=="Prairie"))

kruskal.test(richness ~ Origin, subset(metadata_fungi_mSeq, Crop=="Switchgrass"))
kruskal.test(richness ~ Depth, subset(metadata_fungi_mSeq, Crop=="Switchgrass"))
aov(richness ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Switchgrass"))
summary(aov(richness ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Switchgrass")))
par(mfcol=c(2,2))
plot(aov(richness ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Switchgrass")))
shapiro.test(residuals(aov(richness ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Switchgrass"))))
leveneTest(richness ~ Origin * Depth, data=subset(metadata_fungi_mSeq, Crop=="Prairie"))

# Richness Prokaryotes
kruskal.test(richness ~ Origin, subset(metadata_prok_mSeq, Crop=="Poplar"))
kruskal.test(richness ~ Depth, subset(metadata_prok_mSeq, Crop=="Poplar"))
aov(richness ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Poplar"))
summary(aov(richness ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Poplar")))
par(mfcol=c(2,2))
plot(aov(richness ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Poplar")))
shapiro.test(residuals(aov(richness ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Poplar"))))
leveneTest(richness ~ Origin * Depth, data=subset(metadata_prok_mSeq, Crop=="Poplar"))

kruskal.test(richness ~ Origin, subset(metadata_prok_mSeq, Crop=="Prairie"))
kruskal.test(richness ~ Depth, subset(metadata_prok_mSeq, Crop=="Prairie"))
aov(richness ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Prairie"))
summary(aov(richness ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Prairie")))
par(mfcol=c(2,2))
plot(aov(richness ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Prairie")))
shapiro.test(residuals(aov(richness ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Prairie"))))
leveneTest(richness ~ Origin * Depth, data=subset(metadata_prok_mSeq, Crop=="Prairie"))

kruskal.test(richness ~ Origin, subset(metadata_prok_mSeq, Crop=="Switchgrass"))
kruskal.test(richness ~ Depth, subset(metadata_prok_mSeq, Crop=="Switchgrass"))
aov(richness ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Switchgrass"))
summary(aov(richness ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Switchgrass")))
par(mfcol=c(2,2))
plot(aov(richness ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Switchgrass")))
shapiro.test(residuals(aov(richness ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Switchgrass"))))
leveneTest(richness ~ Origin * Depth, data=subset(metadata_prok_mSeq, Crop=="Prairie"))

# Shannon equitabiliy index fungi
kruskal.test(Eshannon ~ Origin, subset(metadata_fungi_mSeq, Crop=="Poplar"))
kruskal.test(Eshannon ~ Depth, subset(metadata_fungi_mSeq, Crop=="Poplar"))
summary(aov(Eshannon ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Poplar")))
par(mfcol=c(2,2))
plot(aov(Eshannon ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Poplar")))
shapiro.test(residuals(aov(Eshannon ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Poplar"))))
leveneTest(Eshannon ~ Origin * Depth, data=subset(metadata_fungi_mSeq, Crop=="Poplar"))

kruskal.test(Eshannon ~ Origin, subset(metadata_fungi_mSeq, Crop=="Prairie"))
kruskal.test(Eshannon ~ Depth, subset(metadata_fungi_mSeq, Crop=="Prairie"))
aov(Eshannon ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Prairie"))
summary(aov(Eshannon ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Prairie")))
par(mfcol=c(2,2))
plot(aov(Eshannon ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Prairie")))
shapiro.test(residuals(aov(Eshannon ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Prairie"))))
leveneTest(Eshannon ~ Origin * Depth, data=subset(metadata_fungi_mSeq, Crop=="Prairie"))

kruskal.test(Eshannon ~ Origin, subset(metadata_fungi_mSeq, Crop=="Switchgrass"))
kruskal.test(Eshannon ~ Depth, subset(metadata_fungi_mSeq, Crop=="Switchgrass"))
aov(Eshannon ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Switchgrass"))
summary(aov(Eshannon ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Switchgrass")))
par(mfcol=c(2,2))
plot(aov(Eshannon ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Switchgrass")))
shapiro.test(residuals(aov(Eshannon ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Switchgrass"))))
leveneTest(Eshannon ~ Origin * Depth, data=subset(metadata_fungi_mSeq, Crop=="Prairie"))

# Shannon equitabiliy index Prokaryotes
kruskal.test(Eshannon ~ Origin, subset(metadata_prok_mSeq, Crop=="Poplar"))
kruskal.test(Eshannon ~ Depth, subset(metadata_prok_mSeq, Crop=="Poplar"))
aov(Eshannon ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Poplar"))
summary(aov(Eshannon ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Poplar")))
par(mfcol=c(2,2))
plot(aov(Eshannon ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Poplar")))
shapiro.test(residuals(aov(Eshannon ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Poplar"))))
leveneTest(Eshannon ~ Origin * Depth, data=subset(metadata_prok_mSeq, Crop=="Poplar"))

kruskal.test(Eshannon ~ Origin, subset(metadata_prok_mSeq, Crop=="Prairie"))
kruskal.test(Eshannon ~ Depth, subset(metadata_prok_mSeq, Crop=="Prairie"))
aov(Eshannon ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Prairie"))
summary(aov(Eshannon ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Prairie")))
par(mfcol=c(2,2))
plot(aov(Eshannon ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Prairie")))
shapiro.test(residuals(aov(Eshannon ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Prairie"))))
leveneTest(Eshannon ~ Origin * Depth, data=subset(metadata_prok_mSeq, Crop=="Prairie"))

kruskal.test(Eshannon ~ Origin, subset(metadata_prok_mSeq, Crop=="Switchgrass"))
kruskal.test(Eshannon ~ Depth, subset(metadata_prok_mSeq, Crop=="Switchgrass"))
aov(Eshannon ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Switchgrass"))
summary(aov(Eshannon ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Switchgrass")))
par(mfcol=c(2,2))
plot(aov(Eshannon ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Switchgrass")))
shapiro.test(residuals(aov(Eshannon ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Switchgrass"))))
leveneTest(Eshannon ~ Origin * Depth, data=subset(metadata_prok_mSeq, Crop=="Prairie"))

# Diagnostic plots in ggplot2 
# In a loess fit, the alpha parameter determines the width of the sliding window. More 
# specifically, alpha gives the proportion of observations that is to be used in each 
# local regression. Accordingly, this parameter is specified as a value between 0 and 1.
# If the alpha value used for a loess curve is 0.65, each of the local
# regressions used to produce that curve incorporates 65% of the total data points.

# Plotting Diagnostic plots in ggplot2 -----------------------------------------------------
diagPlot<-function(model){
  require(grid)
  require(gridExtra)
  broom::augment(model) -> df
  ggplot(df, aes(.fitted, .resid)) + 
    geom_point(shape=1, size= 1,  na.rm=TRUE) +
    stat_smooth(method="loess",col="red", linetype="solid", span = 1) + 
    geom_hline(yintercept=0, col="black", linetype="dashed") +
    theme_classic() +
    labs(x="Fitted values", y="Residuals", title = "Residual vs Fitted") +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) + 
    theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) -> p1
  ggplot(df, aes(x=.resid)) + 
    geom_density(color="black") +
    geom_vline(aes(xintercept=mean(.resid)),color="red", linetype="dashed", size=1) +
    theme_classic() +
    labs(title  = "Density of Residuals", x = "Residuals", y="Density") +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) + 
    theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7))  -> p2
  ggplot(df, aes(.fitted, sqrt(abs(.std.resid)))) + 
    geom_point(shape=1, size=1, na.rm=TRUE) + 
    stat_smooth(method="loess", col="red", linetype="solid", na.rm = TRUE, span = 1) + 
    theme_classic() +
    labs(x="Fitted Value", y= expression(sqrt("|Standardized residuals|")), title="Scale-Location") +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) + 
    theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7))  -> p3
  
  ggplot(df, aes(sample = .std.resid)) +
    geom_qq(color = "black", shape=1, size=1) + 
    stat_qq_line(color="red") +
    theme_classic() +
    labs(x="Theoretical Quantiles", y="Standardized Residuals", title="Normal Q-Q") +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) + 
    theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7))  -> p4
  ggplot(df, aes(seq_along(cooks.distance(model)), cooks.distance(model))) +
    geom_bar(stat="identity", position="identity", color="grey", fill="grey") +
    theme_classic() +
    labs(x="Obs. Number", y="Cook's distance", title="Cook's distance") +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) + 
    theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7))  -> p5
  ggplot(df, aes(.hat, .std.resid)) +
    geom_point(size=1, shape=1, na.rm = TRUE) +
    expand_limits(x = 0) + 
    geom_smooth(method="lm",col="red", linetype="solid", span = 1, na.rm = TRUE) + 
    geom_hline(yintercept=0, linetype="dashed",color="grey") +
    geom_vline(xintercept=0, linetype="dashed",color="grey") +
    theme_classic() + 
    labs(x="Leverage", y= "Standardized Residuals", title="Residual vs Leverage") +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) + 
    theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7))  -> p6
  return(list(rvfPlot=p1, qqPlot=p2, sclLocPlot=p3, cdPlot=p4, rvlevPlot=p5, cvlPlot=p6))
}

title=text_grob("Diagnostic Plots",size=12,face=2)

diagPlts <- diagPlot(aov(richness ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Poplar")))
grid.arrange(grobs = diagPlts, top = title, ncol=3)
diagPlts <- diagPlot(aov(richness ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Prairie")))
grid.arrange(grobs = diagPlts, top = title, ncol=3)
diagPlts <- diagPlot(aov(richness ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Switchgrass")))
grid.arrange(grobs = diagPlts, top = title, ncol=3)

diagPlts <- diagPlot(aov(richness ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Poplar")))
grid.arrange(grobs = diagPlts, top = title, ncol=3)
diagPlts <- diagPlot(aov(richness ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Prairie")))
grid.arrange(grobs = diagPlts, top = title, ncol=3)
diagPlts <- diagPlot(aov(richness ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Switchgrass")))
grid.arrange(grobs = diagPlts, top = title, ncol=3)


diagPlts <- diagPlot(aov(Eshannon ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Poplar")))
grid.arrange(grobs = diagPlts, top = title, ncol=3)
diagPlts <- diagPlot(aov(Eshannon ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Prairie")))
grid.arrange(grobs = diagPlts, top = title, ncol=3)
diagPlts <- diagPlot(aov(Eshannon ~ Origin * Depth, subset(metadata_fungi_mSeq, Crop=="Switchgrass")))
grid.arrange(grobs = diagPlts, top = title, ncol=3)

diagPlts <- diagPlot(aov(Eshannon ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Poplar")))
grid.arrange(grobs = diagPlts, top = title, ncol=3)
diagPlts <- diagPlot(aov(Eshannon ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Prairie")))
grid.arrange(grobs = diagPlts, top = title, ncol=3)
diagPlts <- diagPlot(aov(Eshannon ~ Origin * Depth, subset(metadata_prok_mSeq, Crop=="Switchgrass")))
grid.arrange(grobs = diagPlts, top = title, ncol=3)

# ********************************************----------------------------------------------
# INDICATOR SPECIES ANALYSIS ---------------------------------------------------------------
library(indicspecies)

GetIndicators <-function(dataframe, var){
  require(phyloseq) 
  require(indicspecies)
  require(dplyr)
  otu <- as(otu_table(dataframe), "matrix")
  otu <- as.data.frame(otu)
  metadata = as(sample_data(dataframe), "data.frame")
  taxa <- as(tax_table(dataframe), "matrix")
  taxa <- as.data.frame(taxa)
  multipatt <- multipatt(t(otu), metadata[,var], func = "IndVal.g",
                         control=how(nperm=999), duleg=TRUE)
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


# indicator p.value <=0.05 after fdr correction
GetIndicators(physeq_fungi_mSeq, "Origin") -> ind_fungi_niche
dim(ind_fungi_niche)
subset(ind_fungi_niche, ind_fungi_niche$stat>=0.75)-> ind_fungi_niche

GetIndicators(physeq_fungi_mSeq, "Crop") -> ind_fungi_crop
dim(ind_fungi_crop)
subset(ind_fungi_crop, ind_fungi_crop$stat>=0.75) -> ind_fungi_crop

GetIndicators(physeq_fungi_mSeq, "Depth") -> ind_fungi_depth
dim(ind_fungi_depth)
subset(ind_fungi_depth, ind_fungi_depth$stat>=0.75) -> ind_fungi_depth


GetIndicators(physeq_prok_mSeq, "Origin") -> ind_prok_niche
dim(ind_prok_niche)
subset(ind_prok_niche, ind_prok_niche$stat>=0.75) -> ind_prok_niche

GetIndicators(physeq_prok_mSeq, "Crop") -> ind_prok_crop
dim(ind_prok_crop)
subset(ind_prok_crop, ind_prok_crop$stat>=0.75) -> ind_prok_crop

GetIndicators(physeq_prok_mSeq, "Depth") -> ind_prok_depth
dim(ind_prok_depth)
subset(ind_prok_depth, ind_prok_depth$stat>=0.75) -> ind_prok_depth


# ********************************************----------------------------------------------
# EXPORT OBJECTS ---------------------------------------------------------------------------
save(physeq_fungi_filt, physeq_prok_filt,
     physeq_fungi_mSeq, physeq_prok_mSeq,
     ind_fungi_crop, ind_fungi_depth, ind_fungi_niche,
     ind_prok_crop, ind_prok_depth, ind_prok_niche,
     file = "filtered_data_for_networks.RData")


# ********************************************----------------------------------------------
# EXTRACTING CORE OTUs ---------------------------------------------------------------------
# Modified from Stopnisek and Shade 2019 - Current Opinion in Microbiology
source("../R_functions/ExtarctCore.R")

# Fungi
physeq_fungi_Pop <-
  physeq_fungi_filt %>%
  subset_samples(Crop %in% c("Poplar"))
otu_table(physeq_fungi_Pop) <-
  otu_table(physeq_fungi_Pop)[which(taxa_sums(physeq_fungi_Pop) > 0),]

dplyr::count(as.data.frame(as.matrix(sample_data(physeq_fungi_Pop))), vars = Depth, Origin)

ExtractCore(physeq_fungi_Pop, Var = Origin, method = "increase", Group="Depth", Level="10") -> otu_core_fungi_Pop_10
otu_core_fungi_Pop_10[1]
ExtractCore(physeq_fungi_Pop, Var = Origin, method = "increase", Group="Depth", Level="25") -> otu_core_fungi_Pop_25
otu_core_fungi_Pop_25[1]
ExtractCore(physeq_fungi_Pop, Var = Origin, method = "increase", Group="Depth", Level="50") -> otu_core_fungi_Pop_50
otu_core_fungi_Pop_50[1]
ExtractCore(physeq_fungi_Pop, Var = Origin, method = "increase", Group="Depth", Level="100") -> otu_core_fungi_Pop_100
otu_core_fungi_Pop_100[1]


physeq_fungi_Pra <-
  physeq_fungi_filt %>%
  subset_samples(Crop %in% c("Prairie"))
otu_table(physeq_fungi_Pra) <-
  otu_table(physeq_fungi_Pra)[which(taxa_sums(physeq_fungi_Pra) > 0),]

dplyr::count(as.data.frame(as.matrix(sample_data(physeq_fungi_Pra))), vars = Depth, Origin)

ExtractCore(physeq_fungi_Pra, Var = Origin, method = "increase",Group="Depth", Level="10") -> otu_core_fungi_Pra_10
otu_core_fungi_Pra_10[1]
ExtractCore(physeq_fungi_Pra, Var = Origin, method = "increase",Group="Depth", Level="25") -> otu_core_fungi_Pra_25
otu_core_fungi_Pra_25[1]
ExtractCore(physeq_fungi_Pra, Var = Origin, method = "increase",Group="Depth", Level="50") ->  otu_core_fungi_Pra_50
otu_core_fungi_Pra_50[1]
ExtractCore(physeq_fungi_Pra, Var = Origin, method = "increase",Group="Depth", Level="100") -> otu_core_fungi_Pra_100
otu_core_fungi_Pra_100[1]

physeq_fungi_Sw <-
  physeq_fungi_filt %>%
  subset_samples(Crop %in% c("Switchgrass"))
otu_table(physeq_fungi_Sw) <-
  otu_table(physeq_fungi_Sw)[which(taxa_sums(physeq_fungi_Sw) > 0),]

dplyr::count(as.data.frame(as.matrix(sample_data(physeq_fungi_Sw))), vars = Depth, Origin)

ExtractCore(physeq_fungi_Sw, Var = Origin, method = "increase", Group="Depth", Level="10") -> otu_core_fungi_Sw_10
otu_core_fungi_Sw_10[1]
ExtractCore(physeq_fungi_Sw, Var = Origin, method = "increase", Group="Depth", Level="25") -> otu_core_fungi_Sw_25
otu_core_fungi_Sw_25[1]
ExtractCore(physeq_fungi_Sw, Var = Origin, method = "increase", Group="Depth", Level="50") ->  otu_core_fungi_Sw_50
otu_core_fungi_Sw_50[1]
ExtractCore(physeq_fungi_Sw, Var = Origin, method = "increase", Group="Depth", Level="100")-> otu_core_fungi_Sw_100 
otu_core_fungi_Sw_100[1]

# Prokaryotes 
physeq_prok_Pop <-
  physeq_prok_filt %>%
  subset_samples(Crop %in% c("Poplar"))
otu_table(physeq_prok_Pop) <-
  otu_table(physeq_prok_Pop)[which(taxa_sums(physeq_prok_Pop) > 0),]

dplyr::count(as.data.frame(as.matrix(sample_data(physeq_prok_Pop))), vars = Depth, Origin)

ExtractCore(physeq_prok_Pop, Var = Origin, method = "increase", Group="Depth", Level="10") -> otu_core_prok_Pop_10
otu_core_prok_Pop_10[1]
ExtractCore(physeq_prok_Pop, Var = Origin, method = "increase", Group="Depth", Level="25") -> otu_core_prok_Pop_25
otu_core_prok_Pop_25[1]
ExtractCore(physeq_prok_Pop, Var = Origin, method = "increase", Group="Depth", Level="50") ->  otu_core_prok_Pop_50
otu_core_prok_Pop_50[1]
ExtractCore(physeq_prok_Pop, Var = Origin, method = "increase", Group="Depth", Level="100") -> otu_core_prok_Pop_100
otu_core_prok_Pop_100[1]


physeq_prok_Pra <-
  physeq_prok_filt %>%
  subset_samples(Crop %in% c("Prairie"))
otu_table(physeq_prok_Pra) <-
  otu_table(physeq_prok_Pra)[which(taxa_sums(physeq_prok_Pra) > 0),]

dplyr::count(as.data.frame(as.matrix(sample_data(physeq_prok_Pra))), vars = Depth, Origin)

ExtractCore(physeq_prok_Pra, Var = Origin, method = "increase", Group="Depth", Level="10") -> otu_core_prok_Pra_10
otu_core_prok_Pra_10[1]
ExtractCore(physeq_prok_Pra, Var = Origin, method = "increase", Group="Depth", Level="25") -> otu_core_prok_Pra_25
otu_core_prok_Pra_25[1]
ExtractCore(physeq_prok_Pra, Var = Origin, method = "increase", Group="Depth", Level="50") ->  otu_core_prok_Pra_50
otu_core_prok_Pra_50[1]
ExtractCore(physeq_prok_Pra, Var = Origin, method = "increase", Group="Depth", Level="100") -> otu_core_prok_Pra_100
otu_core_prok_Pra_100[1]


physeq_prok_Sw <-
  physeq_prok_filt %>%
  subset_samples(Crop %in% c("Switchgrass"))
otu_table(physeq_prok_Sw) <-
  otu_table(physeq_prok_Sw)[which(taxa_sums(physeq_prok_Sw) > 0),]

dplyr::count(as.data.frame(as.matrix(sample_data(physeq_prok_Sw))), vars = Depth, Origin)

ExtractCore(physeq_prok_Sw, Var = Origin, method = "increase", Group="Depth", Level="10") -> otu_core_prok_Sw_10
otu_core_prok_Sw_10[1]
ExtractCore(physeq_prok_Sw, Var = Origin, method = "increase", Group="Depth", Level="25") -> otu_core_prok_Sw_25
otu_core_prok_Sw_25[1]
ExtractCore(physeq_prok_Sw, Var = Origin, method = "increase", Group="Depth", Level="50") ->  otu_core_prok_Sw_50
otu_core_prok_Sw_50[1]
ExtractCore(physeq_prok_Sw, Var = Origin, method = "increase", Group="Depth", Level="100") -> otu_core_prok_Sw_100


# Plotting Bray-Curtis Increase ------------------------------------------------------------
# Threshold set to 1.03 equal to 3% increase
ggarrange(PlotBCincrease(otu_core_fungi_Pop_10, 400) + labs(title = "0-10cm"),
          PlotBCincrease(otu_core_fungi_Pop_25, 400) + labs(title = "10-25cm"),
          PlotBCincrease(otu_core_fungi_Pop_50, 400) + labs(title = "25-50cm"),
          PlotBCincrease(otu_core_fungi_Pop_100, 400) + labs(title = "50-100cm"),
          labels = c("A","","",""),
          ncol = 1, 
          nrow = 4)

ggarrange(PlotBCincrease(otu_core_prok_Pop_10, 400) + labs(title = "0-10cm"),
          PlotBCincrease(otu_core_prok_Pop_25, 400) + labs(title = "10-25cm"),
          PlotBCincrease(otu_core_prok_Pop_50, 400) + labs(title = "25-50cm"),
          PlotBCincrease(otu_core_prok_Pop_100, 400) + labs(title = "50-100cm"),
          labels = c("A","","",""),
          ncol = 1, 
          nrow = 4)

ggarrange(PlotBCincrease(otu_core_fungi_Pra_10, 400) + labs(title = "0-10cm"),
          PlotBCincrease(otu_core_fungi_Pra_25, 400) + labs(title = "10-25cm"),
          PlotBCincrease(otu_core_fungi_Pra_50, 400) + labs(title = "25-50cm"),
          PlotBCincrease(otu_core_fungi_Pra_100, 400) + labs(title = "50-100cm"),
          labels = c("A","","",""),
          ncol = 1, 
          nrow = 4)

ggarrange(PlotBCincrease(otu_core_prok_Pra_10, 400) + labs(title = "0-10cm"),
          PlotBCincrease(otu_core_prok_Pra_25, 400) + labs(title = "10-25cm"),
          PlotBCincrease(otu_core_prok_Pra_50, 400) + labs(title = "25-50cm"),
          PlotBCincrease(otu_core_prok_Pra_100, 400) + labs(title = "50-100cm"),
          labels = c("A","","",""),
          ncol = 1, 
          nrow = 4)

ggarrange(PlotBCincrease(otu_core_fungi_Sw_10, 400) + labs(title = "0-10cm"),
          PlotBCincrease(otu_core_fungi_Sw_25, 400) + labs(title = "10-25cm"),
          PlotBCincrease(otu_core_fungi_Sw_50, 400) + labs(title = "25-50cm"),
          PlotBCincrease(otu_core_fungi_Sw_100, 400) + labs(title = "50-100cm"),
          labels = c("A","","",""),
          ncol = 1, 
          nrow = 4)

ggarrange(PlotBCincrease(otu_core_prok_Sw_10, 400) + labs(title = "0-10cm"),
          PlotBCincrease(otu_core_prok_Sw_25, 400) + labs(title = "10-25cm"),
          PlotBCincrease(otu_core_prok_Sw_50, 400) + labs(title = "25-50cm"),
          PlotBCincrease(otu_core_prok_Sw_100, 400) + labs(title = "50-100cm"),
          labels = c("A","","",""),
          ncol = 1, 
          nrow = 4)


# NEUTRAL MODELS ---------------------------------------------------------------------------
#Fitting neutral model (Burns et al., 2016 (ISME J) - functions are in the sncm.fit.R)
source("../R_functions/sncm.fit.R")
source("../R_functions/ExtractCore.R")

# Poplar
FitNeutral(otu_core_fungi_Pop_10) -> nfit_fungi_Pop_10
FitNeutral(otu_core_fungi_Pop_25) -> nfit_fungi_Pop_25
FitNeutral(otu_core_fungi_Pop_50) -> nfit_fungi_Pop_50
FitNeutral(otu_core_fungi_Pop_100) -> nfit_fungi_Pop_100

FitNeutral(otu_core_prok_Pop_10) -> nfit_prok_Pop_10
FitNeutral(otu_core_prok_Pop_25) -> nfit_prok_Pop_25
FitNeutral(otu_core_prok_Pop_50) -> nfit_prok_Pop_50
FitNeutral(otu_core_prok_Pop_100) -> nfit_prok_Pop_100

# Prairie
FitNeutral(otu_core_fungi_Pra_10) -> nfit_fungi_Pra_10
FitNeutral(otu_core_fungi_Pra_25) -> nfit_fungi_Pra_25
FitNeutral(otu_core_fungi_Pra_50) -> nfit_fungi_Pra_50
FitNeutral(otu_core_fungi_Pra_100) -> nfit_fungi_Pra_100

FitNeutral(otu_core_prok_Pra_10) -> nfit_prok_Pra_10
FitNeutral(otu_core_prok_Pra_25) -> nfit_prok_Pra_25
FitNeutral(otu_core_prok_Pra_50) -> nfit_prok_Pra_50
FitNeutral(otu_core_prok_Pra_100) -> nfit_prok_Pra_100

# Switchgrass
FitNeutral(otu_core_fungi_Sw_10) -> nfit_fungi_Sw_10
FitNeutral(otu_core_fungi_Sw_25) -> nfit_fungi_Sw_25
FitNeutral(otu_core_fungi_Sw_50) -> nfit_fungi_Sw_50
FitNeutral(otu_core_fungi_Sw_100) -> nfit_fungi_Sw_100 

FitNeutral(otu_core_prok_Sw_10) -> nfit_prok_Sw_10
FitNeutral(otu_core_prok_Sw_25) -> nfit_prok_Sw_25
FitNeutral(otu_core_prok_Sw_50) -> nfit_prok_Sw_50
FitNeutral(otu_core_prok_Sw_100) -> nfit_prok_Sw_100


# plotting neutral models Poplar
ggarrange(PlotNeutral(nfit_fungi_Pop_10) + labs(title = "0-10cm"),
          PlotNeutral(nfit_fungi_Pop_25) + labs(title = "10-25cm"),
          PlotNeutral(nfit_fungi_Pop_50) + labs(title = "25-50cm"),
          PlotNeutral(nfit_fungi_Pop_100) + labs(title = "50-100cm"),
          labels = c("A","","",""),
          ncol = 4, 
          nrow = 1)

ggarrange(PlotNeutral(nfit_prok_Pop_10) + labs(title = "0-10cm"),
          PlotNeutral(nfit_prok_Pop_25) + labs(title = "10-25cm"),
          PlotNeutral(nfit_prok_Pop_50) + labs(title = "25-50cm"),
          PlotNeutral(nfit_prok_Pop_100) + labs(title = "50-100cm"),
          labels = c("A","","",""),
          ncol = 4, 
          nrow = 1)

# plotting neutral models Prairie
ggarrange(PlotNeutral(nfit_fungi_Pra_10) + labs(title = "0-10cm"),
          PlotNeutral(nfit_fungi_Pra_25) + labs(title = "10-25cm"),
          PlotNeutral(nfit_fungi_Pra_50) + labs(title = "25-50cm"),
          PlotNeutral(nfit_fungi_Pra_100) + labs(title = "50-100cm"),
          labels = c("A","","",""),
          ncol = 4, 
          nrow = 1)

ggarrange(PlotNeutral(nfit_prok_Pra_10) + labs(title = "0-10cm"),
          PlotNeutral(nfit_prok_Pra_25) + labs(title = "10-25cm"),
          PlotNeutral(nfit_prok_Pra_50) + labs(title = "25-50cm"),
          PlotNeutral(nfit_prok_Pra_100) + labs(title = "50-100cm"),
          labels = c("A","","",""),
          ncol = 4, 
          nrow = 1)

# plotting neutral models Switchgrass
ggarrange(PlotNeutral(nfit_fungi_Sw_10) + labs(title = "0-10cm"),
          PlotNeutral(nfit_fungi_Sw_25) + labs(title = "10-25cm"),
          PlotNeutral(nfit_fungi_Sw_50) + labs(title = "25-50cm"),
          PlotNeutral(nfit_fungi_Sw_100) + labs(title = "50-100cm"),
          labels = c("A","","",""),
          ncol = 4, 
          nrow = 1)

ggarrange(PlotNeutral(nfit_prok_Sw_10) + labs(title = "0-10cm"),
          PlotNeutral(nfit_prok_Sw_25) + labs(title = "10-25cm"),
          PlotNeutral(nfit_prok_Sw_50) + labs(title = "25-50cm"),
          PlotNeutral(nfit_prok_Sw_100) + labs(title = "50-100cm"),
          labels = c("A","","",""),
          ncol = 4, 
          nrow = 1)


title1=text_grob("Poplar Fungi",size=10, face=2)
title2=text_grob("Poplar Prokaryotes",size=10, face=2)

title3=text_grob("Prairie Fungi",size=10, face=2)
title4=text_grob("Prairie Prokaryotes",size=10, face=2)


# FIGURE S4 - Neutral models ---------------------------------------------------------------
Neutral_all <-
  ggarrange(
    grid.arrange(ggarrange(PlotNeutral(nfit_fungi_Pop_10) + labs(title = "0-10cm"),
                           PlotNeutral(nfit_fungi_Pop_25) + labs(title = "10-25cm"),
                           PlotNeutral(nfit_fungi_Pop_50) + labs(title = "25-50cm"),
                           PlotNeutral(nfit_fungi_Pop_100) + labs(title = "50-100cm"),
                           labels = "A",
                           ncol = 4,
                           nrow = 1),
                 top = title1),
    grid.arrange(ggarrange(PlotNeutral(nfit_prok_Pop_10) + labs(title = "0-10cm"),
                           PlotNeutral(nfit_prok_Pop_25) + labs(title = "10-25cm"),
                           PlotNeutral(nfit_prok_Pop_50) + labs(title = "25-50cm"),
                           PlotNeutral(nfit_prok_Pop_100) + labs(title = "50-100cm"),
                           labels = "B",
                           ncol = 4,
                           nrow = 1),
                 top = title2),
    grid.arrange(ggarrange(PlotNeutral(nfit_fungi_Pra_10) + labs(title = "0-10cm"),
                           PlotNeutral(nfit_fungi_Pra_25) + labs(title = "10-25cm"),
                           PlotNeutral(nfit_fungi_Pra_50) + labs(title = "25-50cm"),
                           PlotNeutral(nfit_fungi_Pra_100) + labs(title = "50-100cm"),
                           labels = "C",
                           ncol = 4,
                           nrow = 1),
                 top = title3),
    grid.arrange(ggarrange(PlotNeutral(nfit_prok_Pra_10) + labs(title = "0-10cm"),
                           PlotNeutral(nfit_prok_Pra_25) + labs(title = "10-25cm"),
                           PlotNeutral(nfit_prok_Pra_50) + labs(title = "25-50cm"),
                           PlotNeutral(nfit_prok_Pra_100) + labs(title = "50-100cm"),
                           labels = "D",
                           ncol = 4,
                           nrow = 1),
                 top = title4),
    ncol = 1, 
    nrow = 4)

Neutral_all

# ********************************************----------------------------------------------
# MICROBIAL NETWORKS -----------------------------------------------------------------------
source("../R_functions/MicrobNet.R")

# Extracting phyloseq objects
subset_taxa(physeq_fungi_Pop, OTU_ID %in% unlist(otu_core_fungi_Pop_10[1])) -> physeq_fungi_core_Pop_10
subset_taxa(physeq_fungi_Pop, OTU_ID %in% unlist(otu_core_fungi_Pop_25[1])) -> physeq_fungi_core_Pop_25
subset_taxa(physeq_fungi_Pop, OTU_ID %in% unlist(otu_core_fungi_Pop_50[1])) -> physeq_fungi_core_Pop_50
subset_taxa(physeq_fungi_Pop, OTU_ID %in% unlist(otu_core_fungi_Pop_100[1])) -> physeq_fungi_core_Pop_100
subset_taxa(physeq_prok_Pop, OTU_ID %in% unlist(otu_core_prok_Pop_10[1])) -> physeq_prok_core_Pop_10
subset_taxa(physeq_prok_Pop, OTU_ID %in% unlist(otu_core_prok_Pop_25[1])) -> physeq_prok_core_Pop_25
subset_taxa(physeq_prok_Pop, OTU_ID %in% unlist(otu_core_prok_Pop_50[1])) -> physeq_prok_core_Pop_50
subset_taxa(physeq_prok_Pop, OTU_ID %in% unlist(otu_core_prok_Pop_100[1])) -> physeq_prok_core_Pop_100

subset_taxa(physeq_fungi_Pra, OTU_ID %in% unlist(otu_core_fungi_Pra_10[1])) -> physeq_fungi_core_Pra_10
subset_taxa(physeq_fungi_Pra, OTU_ID %in% unlist(otu_core_fungi_Pra_25[1])) -> physeq_fungi_core_Pra_25
subset_taxa(physeq_fungi_Pra, OTU_ID %in% unlist(otu_core_fungi_Pra_50[1])) -> physeq_fungi_core_Pra_50
subset_taxa(physeq_fungi_Pra, OTU_ID %in% unlist(otu_core_fungi_Pra_100[1])) -> physeq_fungi_core_Pra_100
subset_taxa(physeq_prok_Pra, OTU_ID %in% unlist(otu_core_prok_Pra_10[1])) -> physeq_prok_core_Pra_10
subset_taxa(physeq_prok_Pra, OTU_ID %in% unlist(otu_core_prok_Pra_25[1])) -> physeq_prok_core_Pra_25
subset_taxa(physeq_prok_Pra, OTU_ID %in% unlist(otu_core_prok_Pra_50[1])) -> physeq_prok_core_Pra_50
subset_taxa(physeq_prok_Pra, OTU_ID %in% unlist(otu_core_prok_Pra_100[1])) -> physeq_prok_core_Pra_100

subset_taxa(physeq_fungi_Sw, OTU_ID %in% unlist(otu_core_fungi_Sw_10[1])) -> physeq_fungi_core_Sw_10
subset_taxa(physeq_fungi_Sw, OTU_ID %in% unlist(otu_core_fungi_Sw_25[1])) -> physeq_fungi_core_Sw_25
subset_taxa(physeq_fungi_Sw, OTU_ID %in% unlist(otu_core_fungi_Sw_50[1])) -> physeq_fungi_core_Sw_50
subset_taxa(physeq_fungi_Sw, OTU_ID %in% unlist(otu_core_fungi_Sw_100[1])) -> physeq_fungi_core_Sw_100
subset_taxa(physeq_prok_Sw, OTU_ID %in% unlist(otu_core_prok_Sw_10[1])) -> physeq_prok_core_Sw_10
subset_taxa(physeq_prok_Sw, OTU_ID %in% unlist(otu_core_prok_Sw_25[1])) -> physeq_prok_core_Sw_25
subset_taxa(physeq_prok_Sw, OTU_ID %in% unlist(otu_core_prok_Sw_50[1])) -> physeq_prok_core_Sw_50
subset_taxa(physeq_prok_Sw, OTU_ID %in% unlist(otu_core_prok_Sw_100[1])) -> physeq_prok_core_Sw_100


# SpiecEasi NETWORKS -----------------------------------------------------------------------
source("../R_functions/MicrobNet.R")

MakeSENet(physeq_fungi_core_Pop_10, physeq_prok_core_Pop_10) -> spiec_net_Pop_10 #100, 1e-1, thresh=0.03
getStability(spiec_net_Pop_10)
plot(adj2igraph(getRefit(spiec_net_Pop_10)), 
     vertex.color=c(rep(1,ntaxa(physeq_fungi_core_Pop_10)), 
                    rep(2,ntaxa(physeq_prok_core_Pop_10))) +1,vertex.size=6)

MakeSENet(physeq_fungi_core_Pop_25, physeq_prok_core_Pop_25) -> spiec_net_Pop_25
getStability(spiec_net_Pop_25)
plot(adj2igraph(getRefit(spiec_net_Pop_10)), 
     vertex.color=c(rep(1,ntaxa(physeq_fungi_core_Pop_25)), 
                    rep(2,ntaxa(physeq_prok_core_Pop_25))) +1,vertex.size=6)

MakeSENet(physeq_fungi_core_Pop_50, physeq_prok_core_Pop_50) -> spiec_net_Pop_50 
getStability(spiec_net_Pop_50)
plot(adj2igraph(getRefit(spiec_net_Pop_50)), 
     vertex.color=c(rep(1,ntaxa(physeq_fungi_core_Pop_50)), 
                    rep(2,ntaxa(physeq_prok_core_Pop_50))) +1,vertex.size=6)

MakeSENet(physeq_fungi_core_Pop_100, physeq_prok_core_Pop_100) -> spiec_net_Pop_100
getStability(spiec_net_Pop_100)
plot(adj2igraph(getRefit(spiec_net_Pop_100)), 
     vertex.color=c(rep(1,ntaxa(physeq_fungi_core_Pop_100)), 
                    rep(2,ntaxa(physeq_prok_core_Pop_100))) +1,vertex.size=6)

save(spiec_net_Pop_10,
     spiec_net_Pop_25,
     spiec_net_Pop_50,
     spiec_net_Pop_100,
     file = "spiec_net_pop.RData")


MakeSENet(physeq_fungi_core_Pra_10, physeq_prok_core_Pra_10) -> spiec_net_Pra_10 #100, 1e-1
getStability(spiec_net_Pra_10)
plot(adj2igraph(getRefit(spiec_net_Pra_10)), 
     vertex.color=c(rep(1,ntaxa(physeq_fungi_core_Pra_10)), 
                    rep(2,ntaxa(physeq_prok_core_Pra_10))) +1,vertex.size=6)


MakeSENet(physeq_fungi_core_Pra_25, physeq_prok_core_Pra_25) -> spiec_net_Pra_25 #100, 1e-1
getStability(spiec_net_Pra_25)
plot(adj2igraph(getRefit(spiec_net_Pra_25)), 
     vertex.color=c(rep(1,ntaxa(physeq_fungi_core_Pra_25)), 
                    rep(2,ntaxa(physeq_prok_core_Pra_25))) +1,vertex.size=6)


MakeSENet(physeq_fungi_core_Pra_50, physeq_prok_core_Pra_50) -> spiec_net_Pra_50 #100, 1e-1
getStability(spiec_net_Pra_50)
plot(adj2igraph(getRefit(spiec_net_Pra_50)), 
     vertex.color=c(rep(1,ntaxa(physeq_fungi_core_Pra_50)), 
                    rep(2,ntaxa(physeq_prok_core_Pra_50))) +1,vertex.size=6)


MakeSENet(physeq_fungi_core_Pra_100, physeq_prok_core_Pra_100) -> spiec_net_Pra_100 #100, 1e-1
getStability(spiec_net_Pra_100)
plot(adj2igraph(getRefit(spiec_net_Pra_100)), 
     vertex.color=c(rep(1,ntaxa(physeq_fungi_core_Pra_100)), 
                    rep(2,ntaxa(physeq_prok_core_Pra_100))) +1,vertex.size=6)

save(spiec_net_Pra_10,
     spiec_net_Pra_25,
     spiec_net_Pra_50,
     spiec_net_Pra_100,
     file = "spiec_net_pra.RData")

MakeSENet(physeq_fungi_core_Sw_10, physeq_prok_core_Sw_10) -> spiec_net_Sw_10 #100, 1e-1
getStability(spiec_net_Sw_10)
plot(adj2igraph(getRefit(spiec_net_Sw_10)), 
     vertex.color=c(rep(1,ntaxa(physeq_fungi_core_Sw_10)), 
                    rep(2,ntaxa(physeq_prok_core_Sw_10))) +1,vertex.size=6)

MakeSENet(physeq_fungi_core_Sw_25, physeq_prok_core_Sw_25) -> spiec_net_Sw_25 #100, 1e-1
getStability(spiec_net_Sw_25)
plot(adj2igraph(getRefit(spiec_net_Sw_25)), 
     vertex.color=c(rep(1,ntaxa(physeq_fungi_core_Sw_25)), 
                    rep(2,ntaxa(physeq_prok_core_Sw_25))) +1,vertex.size=6)

MakeSENet(physeq_fungi_core_Sw_50, physeq_prok_core_Sw_50) -> spiec_net_Sw_50 #100, 1e-1
getStability(spiec_net_Sw_50)
plot(adj2igraph(getRefit(spiec_net_Sw_50)), 
     vertex.color=c(rep(1,ntaxa(physeq_fungi_core_Sw_50)), 
                    rep(2,ntaxa(physeq_prok_core_Sw_50))) +1,vertex.size=6)

MakeSENet(physeq_fungi_core_Sw_100, physeq_prok_core_Sw_100) -> spiec_net_Sw_100 #100, 1e-1
getStability(spiec_net_Sw_100)
plot(adj2igraph(getRefit(spiec_net_Sw_100)), 
     vertex.color=c(rep(1,ntaxa(physeq_fungi_core_Sw_100)), 
                    rep(2,ntaxa(physeq_prok_core_Sw_100))) +1,vertex.size=6)

save(spiec_net_Sw_10,
     spiec_net_Sw_25,
     spiec_net_Sw_50,
     spiec_net_Sw_100,
     file = "spiec_net_sw.RData")

# Generating the networks
GetNetwork(physeq_fungi_core_Pop_10, physeq_prok_core_Pop_10, spiec_net_Pop_10) -> network_Pop_10
GetNetwork(physeq_fungi_core_Pop_25, physeq_prok_core_Pop_25, spiec_net_Pop_25) -> network_Pop_25
GetNetwork(physeq_fungi_core_Pop_50, physeq_prok_core_Pop_50, spiec_net_Pop_50) -> network_Pop_50
GetNetwork(physeq_fungi_core_Pop_100, physeq_prok_core_Pop_100, spiec_net_Pop_100) -> network_Pop_100

GetNetwork(physeq_fungi_core_Pra_10, physeq_prok_core_Pra_10, spiec_net_Pra_10) -> network_Pra_10
GetNetwork(physeq_fungi_core_Pra_25, physeq_prok_core_Pra_25, spiec_net_Pra_25) -> network_Pra_25
GetNetwork(physeq_fungi_core_Pra_50, physeq_prok_core_Pra_50, spiec_net_Pra_50) -> network_Pra_50
GetNetwork(physeq_fungi_core_Pra_100, physeq_prok_core_Pra_100, spiec_net_Pra_100) -> network_Pra_100

GetNetwork(physeq_fungi_core_Sw_10, physeq_prok_core_Sw_10, spiec_net_Sw_10) -> network_Sw_10
GetNetwork(physeq_fungi_core_Sw_25, physeq_prok_core_Sw_25, spiec_net_Sw_25) -> network_Sw_25
GetNetwork(physeq_fungi_core_Sw_50, physeq_prok_core_Sw_50, spiec_net_Sw_50) -> network_Sw_50
GetNetwork(physeq_fungi_core_Sw_100, physeq_prok_core_Sw_100, spiec_net_Sw_100) -> network_Sw_100

# COMPARE TO RANDOM NETWORK ----------------------------------------------------------------
TesToNULL(network_Pop_10)
TesToNULL(network_Pop_25)
TesToNULL(network_Pop_50)
TesToNULL(network_Pop_100)

TesToNULL(network_Pra_10)
TesToNULL(network_Pra_25)
TesToNULL(network_Pra_50)
TesToNULL(network_Pra_100)

TesToNULL(network_Sw_10)
TesToNULL(network_Sw_25)
TesToNULL(network_Sw_50)
TesToNULL(network_Sw_100)

# NETWORK ATTRIBUTES -----------------------------------------------------------------------
source("../R_functions/MicrobNet.R")

# Extract Taxonomy -------------------------------------------------------------------------
ExtarctTaxa(physeq_fungi_filt, physeq_prok_filt) -> taxa_all
head(taxa_all)

# Generate compostie phyloseq objects ----------------------------------------------------- 
GeneratePhyseqNet(physeq_fungi_mSeq, physeq_prok_mSeq, "Crop", "Poplar") -> physeq_fb_Pop
GeneratePhyseqNet(physeq_fungi_mSeq, physeq_prok_mSeq, "Crop", "Prairie") -> physeq_fb_Pra
GeneratePhyseqNet(physeq_fungi_mSeq, physeq_prok_mSeq, "Crop", "Switchgrass") -> physeq_fb_Sw

# Calculate Abundances ---------------------------------------------------------------------
ExtractAbund(physeq_fb_Pop, "Poplar", "0-10cm") -> abund_Pop_10
ExtractAbund(physeq_fb_Pop, "Poplar", "10-25cm") -> abund_Pop_25
ExtractAbund(physeq_fb_Pop, "Poplar", "25-50cm") -> abund_Pop_50
ExtractAbund(physeq_fb_Pop, "Poplar", "50-100cm") -> abund_Pop_100

ExtractAbund(physeq_fb_Pra, "Prairie", "0-10cm") -> abund_Pra_10
ExtractAbund(physeq_fb_Pra, "Prairie", "10-25cm") -> abund_Pra_25
ExtractAbund(physeq_fb_Pra, "Prairie", "25-50cm") -> abund_Pra_50
ExtractAbund(physeq_fb_Pra, "Prairie", "50-100cm") -> abund_Pra_100

ExtractAbund(physeq_fb_Sw, "Switchgrass", "0-10cm") -> abund_Sw_10
ExtractAbund(physeq_fb_Sw, "Switchgrass", "10-25cm") -> abund_Sw_25
ExtractAbund(physeq_fb_Sw, "Switchgrass", "25-50cm") -> abund_Sw_50
ExtractAbund(physeq_fb_Sw, "Switchgrass", "50-100cm") -> abund_Sw_100


# NODE attribtues --------------------------------------------------------------------------
head(CalcNetAttr(network_Pop_10, "10")) 

# Pi-Zi dataframes -------------------------------------------------------------------------
df_zipi_10
df_zipi_10[is.na(df_zipi_10$P), ]
df_zipi_10[is.nan(df_zipi_10$Z), ]


# NODE TABLE -------------------------------------------------------------------------------
CalcNodes(abund_Pop_10, network_Pop_10, nfit_fungi_Pop_10, nfit_prok_Pop_10, "0-10cm", "Poplar") -> nodes_Pop_10
CalcNodes(abund_Pop_25, network_Pop_25, nfit_fungi_Pop_25, nfit_prok_Pop_25, "10-25cm","Poplar") -> nodes_Pop_25
CalcNodes(abund_Pop_50, network_Pop_50, nfit_fungi_Pop_50, nfit_prok_Pop_50, "25-50cm", "Poplar") -> nodes_Pop_50
CalcNodes(abund_Pop_100, network_Pop_100, nfit_fungi_Pop_100, nfit_prok_Pop_100, "50-100cm", "Poplar") -> nodes_Pop_100

CalcNodes(abund_Pra_10, network_Pra_10, nfit_fungi_Pra_10, nfit_prok_Pra_10, "0-10cm", "Prairie") -> nodes_Pra_10
CalcNodes(abund_Pra_25, network_Pra_25, nfit_fungi_Pra_25, nfit_prok_Pra_25, "10-25cm","Prairie") -> nodes_Pra_25
CalcNodes(abund_Pra_50, network_Pra_50, nfit_fungi_Pra_50, nfit_prok_Pra_50, "25-50cm","Prairie") -> nodes_Pra_50
CalcNodes(abund_Pra_100, network_Pra_100, nfit_fungi_Pra_100, nfit_prok_Pra_100, "50-100cm","Prairie") -> nodes_Pra_100

CalcNodes(abund_Sw_10, network_Sw_10, nfit_fungi_Sw_10, nfit_prok_Sw_10, "0-10cm", "Switchgrass") -> nodes_Sw_10
CalcNodes(abund_Sw_25, network_Sw_25, nfit_fungi_Sw_25, nfit_prok_Sw_25, "10-25cm", "Switchgrass") -> nodes_Sw_25
CalcNodes(abund_Sw_50, network_Sw_50, nfit_fungi_Sw_50, nfit_prok_Sw_50, "25-50cm", "Switchgrass") -> nodes_Sw_50
CalcNodes(abund_Sw_100, network_Sw_100, nfit_fungi_Sw_100, nfit_prok_Sw_100, "50-100cm", "Switchgrass") -> nodes_Sw_100


# EDGE TABLE -------------------------------------------------------------------------------
CalcEdges(network_Pop_10, nodes_Pop_10,"0-10cm") -> edges_Pop_10
CalcEdges(network_Pop_25, nodes_Pop_25,"10-25cm") -> edges_Pop_25
CalcEdges(network_Pop_50, nodes_Pop_50,"25-50cm") -> edges_Pop_50
CalcEdges(network_Pop_100, nodes_Pop_100,"50-100cm") -> edges_Pop_100

CalcEdges(network_Pra_10, nodes_Pra_10,"0-10cm") -> edges_Pra_10
CalcEdges(network_Pra_25, nodes_Pra_25,"10-25cm") -> edges_Pra_25
CalcEdges(network_Pra_50, nodes_Pra_50,"25-50cm") -> edges_Pra_50
CalcEdges(network_Pra_100, nodes_Pra_100,"50-100cm") -> edges_Pra_100

CalcEdges(network_Sw_10, nodes_Sw_10,"0-10cm") -> edges_Sw_10
CalcEdges(network_Sw_25, nodes_Sw_25,"10-25cm") -> edges_Sw_25
CalcEdges(network_Sw_50, nodes_Sw_50,"25-50cm") -> edges_Sw_50
CalcEdges(network_Sw_100, nodes_Sw_100,"50-100cm") -> edges_Sw_100

# # Expoering tables 
# write.csv(edges_Pop, "edges_Pop.csv")
# write.csv(nodes_Pop, "nodes_Pop.csv")

# ********************************************----------------------------------------------
# NETWORK STATS SUMMARY --------------------------------------------------------------------
cbind(NetStats(nodes_Pop_10, edges_Pop_10, network_Pop_10, 
               nfit_fungi_Pop_10,nfit_prok_Pop_10,"0-10cm"),
      NetStats(nodes_Pop_25, edges_Pop_25, network_Pop_25, 
               nfit_fungi_Pop_25,nfit_prok_Pop_25,"10-25cm"),
      NetStats(nodes_Pop_50, edges_Pop_50, network_Pop_50, 
               nfit_fungi_Pop_50,nfit_prok_Pop_50,"25-50cm"),
      NetStats(nodes_Pop_100, edges_Pop_100, network_Pop_100,
               nfit_fungi_Pop_100,nfit_prok_Pop_100,"50-100cm")) -> net_stats_pop


cbind(NetStatsCore(nodes_Pop_10, edges_Pop_10, network_Pop_10, 
                   nfit_fungi_Pop_10,nfit_prok_Pop_10,"0-10cm","Poplar"),
      NetStatsCore(nodes_Pop_25, edges_Pop_25, network_Pop_25, 
                   nfit_fungi_Pop_25,nfit_prok_Pop_25,"10-25cm","Poplar"),
      NetStatsCore(nodes_Pop_50, edges_Pop_50, network_Pop_50, 
                   nfit_fungi_Pop_50,nfit_prok_Pop_50,"25-50cm","Poplar"),
      NetStatsCore(nodes_Pop_100, edges_Pop_100, network_Pop_100,
                   nfit_fungi_Pop_100,nfit_prok_Pop_100,"50-100cm","Poplar")) -> net_stats_pop_core

net_stats_pop_core

net_stats_pop
grid.table(net_stats_pop)
grid.draw(tableGrob(net_stats_pop,
                    theme=ttheme_minimal(base_size = 8, 
                                         base_colour = "black", base_family = "",
                                         padding = unit(c(2, 2), "mm"))))
grid.arrange(tableGrob(net_stats_pop))


cbind(NetStats(nodes_Pra_10, edges_Pra_10, network_Pra_10,
               nfit_fungi_Pra_10, nfit_prok_Pra_10, "0-10cm"),
      NetStats(nodes_Pra_25, edges_Pra_25, network_Pra_25, 
               nfit_fungi_Pra_25,nfit_prok_Pra_25,"10-25cm"),
      NetStats(nodes_Pra_50, edges_Pra_50, network_Pra_50,
               nfit_fungi_Pra_50,nfit_prok_Pra_50,"25-50cm"),
      NetStats(nodes_Pra_100, edges_Pra_100, network_Pra_100,
               nfit_fungi_Pra_100,nfit_prok_Pra_100,"50-100cm")) -> net_stats_pra

net_stats_pra

cbind(NetStatsCore(nodes_Pra_10, edges_Pra_10, network_Pra_10,
                   nfit_fungi_Pra_10, nfit_prok_Pra_10, "0-10cm", "Prairie"),
      NetStatsCore(nodes_Pra_25, edges_Pra_25, network_Pra_25, 
                   nfit_fungi_Pra_25,nfit_prok_Pra_25,"10-25cm", "Prairie"),
      NetStatsCore(nodes_Pra_50, edges_Pra_50, network_Pra_50,
                   nfit_fungi_Pra_50,nfit_prok_Pra_50,"25-50cm", "Prairie"),
      NetStatsCore(nodes_Pra_100, edges_Pra_100, network_Pra_100,
                   nfit_fungi_Pra_100,nfit_prok_Pra_100,"50-100cm", "Prairie")) -> net_stats_pra_core

net_stats_pra_core


cbind(NetStats(nodes_Sw_10, edges_Sw_10, network_Sw_10,
               nfit_fungi_Sw_10, nfit_prok_Sw_10, "0-10cm"),
      NetStats(nodes_Sw_25, edges_Sw_25, network_Sw_25, 
               nfit_fungi_Sw_25,nfit_prok_Sw_25,"10-25cm"),
      NetStats(nodes_Sw_50, edges_Sw_50, network_Sw_50,
               nfit_fungi_Sw_50,nfit_prok_Sw_50,"25-50cm"),
      NetStats(nodes_Sw_100, edges_Sw_100, network_Sw_100,
               nfit_fungi_Sw_100,nfit_prok_Sw_100,"50-100cm")) -> net_stats_sw

net_stats_sw

cbind(NetStatsCore(nodes_Sw_10, edges_Sw_10, network_Sw_10,
                   nfit_fungi_Sw_10, nfit_prok_Sw_10, "0-10cm", "Switchgrass"),
      NetStatsCore(nodes_Sw_25, edges_Sw_25, network_Sw_25, 
                   nfit_fungi_Sw_25,nfit_prok_Sw_25,"10-25cm", "Switchgrass"),
      NetStatsCore(nodes_Sw_50, edges_Sw_50, network_Sw_50,
                   nfit_fungi_Sw_50,nfit_prok_Sw_50,"25-50cm", "Switchgrass"),
      NetStatsCore(nodes_Sw_100, edges_Sw_100, network_Sw_100,
                   nfit_fungi_Sw_100,nfit_prok_Sw_100,"50-100cm", "Switchgrass")) -> net_stats_sw_core

net_stats_sw_core

write.csv(net_stats_pop,"net_stats_pop.csv")
write.csv(net_stats_pra,"net_stats_pra.csv")
write.csv(net_stats_sw,"net_stats_sw.csv")

write.csv(net_stats_pop_core,"net_stats_pop_core.csv")
write.csv(net_stats_pra_core,"net_stats_pra_core.csv")
write.csv(net_stats_sw_core,"net_stats_sw_core.csv")

# ********************************************----------------------------------------------
# PLOTTING THE NETWORK ---------------------------------------------------------------------

# Generate the table graph -----------------------------------------------------------------
MakeTabGraph(nodes_Pop_10, edges_Pop_10) -> ggraph_Pop_10
MakeTabGraph(nodes_Pop_25, edges_Pop_25) -> ggraph_Pop_25
MakeTabGraph(nodes_Pop_50, edges_Pop_50) -> ggraph_Pop_50
MakeTabGraph(nodes_Pop_100, edges_Pop_100) -> ggraph_Pop_100

MakeTabGraph(nodes_Pra_10, edges_Pra_10) -> ggraph_Pra_10
MakeTabGraph(nodes_Pra_25, edges_Pra_25) -> ggraph_Pra_25
MakeTabGraph(nodes_Pra_50, edges_Pra_50) -> ggraph_Pra_50
MakeTabGraph(nodes_Pra_100, edges_Pra_100) -> ggraph_Pra_100

MakeTabGraph(nodes_Sw_10, edges_Sw_10) -> ggraph_Sw_10
MakeTabGraph(nodes_Sw_25, edges_Sw_25) -> ggraph_Sw_25
MakeTabGraph(nodes_Sw_50, edges_Sw_50) -> ggraph_Sw_50
MakeTabGraph(nodes_Sw_100, edges_Sw_100) -> ggraph_Sw_100

# palette15 <- c("#F99FDF","#C55498","#BF227F","#911F78","#6E5A34","#2A5446",
#                "#117E8D","#2345BC","#2D31A1","#2787DE","#8A808F","#BCB670",
#                "#E4A36D","#E97E65","#20D454", "black") 
# 
# messy_microbial <-c("#E8C0D2","#A8C7A3","#E3C344","#8ED52D",
#                     "#AF6F5D","#C9691F","#121360D2","#0ABECD29")
# 
# palette <- c("#CBD588", "#599861", "#508578","#FFD700","#FFA500",
#              "#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD", 
#              "#D14285", "#C84248", "#8569D5","#673770","#1163cf",
#              "#A9CDFC","#5E738F","#000000")

# network components
components(network_Pop_10)
components(network_Pop_25)
components(network_Pop_50)
components(network_Pop_100)


# Extracting a common legent for plotting --------------------------------------------------
# Generate a Color Palette -----------------------------------------------------------------

palette_phylum <- c(
  Ascomycota = "red",
  Basidiomycota = "orange",
  Calcarisporiellomycota = "gold",  #"cornflowerblue"
  Chytridiomycota = "yellow",
  Glomeromycota = "papayawhip",
  Mortierellomycota = "peachpuff", #"pink"
  Monoblepharomycota = "burlywood",
  Mucoromycota = "salmon4", #"darkgreen"
  Olpidiomycota = "burlywood4",
  Rozellomycota = "brown", #"peachpuff",
  Zoopagomycota =  "gray24",#"burlywood",
  Unclassified = "navy",
  Acidobacteriota = "blue",
  Actinobacteriota = "dodgerblue2",
  Bacteroidota = "cyan", #"wheat",
  Chloroflexi = "lightcyan", #"plum1",
  Crenarchaeota = "plum1", #"#A20E42",
  Desulfobacterota = "magenta",
  Firmicutes = "purple", #"limegreen",
  Gemmatimonadota = "slateblue4", #"olivedrab",
  Methylomirabilota = "darkslategray", #"cyan",
  Myxococcota = "darkseagreen3", #"palevioletred1",
  `Nb1-j` = "forestgreen",
  Nitrospirota = "darkgreen",
  Planctomycetota = "springgreen",
  Proteobacteria = "greenyellow", #"royalblue",
  Verrucomicrobiota  = "grey") #"purple"

factor_levels = c(
  "Ascomycota",
  "Basidiomycota",
  "Calcarisporiellomycota",
  "Chytridiomycota",
  "Glomeromycota",
  "Mortierellomycota",
  "Monoblepharomycota",
  "Mucoromycota",
  "Olpidiomycota",
  "Rozellomycota",
  "Zoopagomycota" ,
  "Unclassified",
  "Acidobacteriota",
  "Actinobacteriota",
  "Bacteroidota",
  "Chloroflexi",
  "Crenarchaeota",
  "Desulfobacterota",
  "Firmicutes",
  "Gemmatimonadota",
  "Methylomirabilota",
  "Myxococcota",
  "Nb1-j",
  "Nitrospirota" ,
  "Planctomycetota",
  "Proteobacteria",
  "Verrucomicrobiota")

# Get common Taxa Legend -------------------------------------------------------------------
PlotLegend <- function(){
  data.frame(Phylum =
               unique(c(nodes_Pop_10$Phylum,
                        nodes_Pop_25$Phylum,
                        nodes_Pop_50$Phylum,
                        nodes_Pop_100$Phylum,
                        nodes_Sw_10$Phylum,
                        nodes_Sw_25$Phylum,
                        nodes_Sw_50$Phylum,
                        nodes_Sw_100$Phylum,
                        nodes_Pra_10$Phylum,
                        nodes_Pra_25$Phylum,
                        nodes_Pra_50$Phylum,
                        nodes_Pra_100$Phylum)),
             Abund = runif(27, 1, 27)) -> df
  print("Unique Phyla")
  df$Phylum %T>% print()
  df$Phylum <- 
    factor(df$Phylum, 
           levels = factor_levels)
  df %>%
    ggplot(aes(x=Abund, y=Phylum, fill=Phylum, color =Phylum)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = palette_phylum) + 
    scale_color_manual(values = palette_phylum) +
    theme_classic() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), 
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title = element_text(angle = 0, size = 7, face = "bold"),
          axis.text.y = element_text(size = 7), 
          legend.key.height = unit(0.25, "cm"), legend.key.width = unit(0.25, "cm"), 
          legend.title = element_text(size = 8, face = "bold", hjust = 0, vjust = 0.5), 
          legend.text = element_text(size = 7), legend.position = "right") +
    guides(fill = guide_legend(ncol = 1, title.position="top"),
           color = guide_legend(ncol = 1, title.position="top")) -> plot
  legend <-
    as_ggplot(get_legend(plot))
  return(legend)
}

PlotLegend()

# Get network legend 
as_ggplot(get_legend(PlotNetwTax(ggraph_Pop_25, "AbundDepth", "10-25cm")))


# **** Poplar Network **** -----------------------------------------------------------------
# Module hubs and indicators
nodes_Pop_10[nodes_Pop_10$Key=="Module hubs",]
nodes_Pop_10[nodes_Pop_10$Key=="Network hubs",]
nodes_Pop_50[nodes_Pop_50$Key=="Module hubs",]
nodes_Pop_100[nodes_Pop_100$Key=="Module hubs",]

nodes_Pop_10[nodes_Pop_10$Depth_ind%in%"0-10cm",]
nodes_Pop_25[nodes_Pop_25$Depth_ind%in%"10-25cm",]
nodes_Pop_50[nodes_Pop_50$Depth_ind%in%"25-50cm",]
nodes_Pop_100[nodes_Pop_100$Depth_ind%in%"50-100cm",]

nodes_Pop_10[nodes_Pop_10$Niche_ind%in%"Root",]
nodes_Pop_25[nodes_Pop_25$Niche_ind%in%"Root",]
nodes_Pop_50[nodes_Pop_50$Depth_ind%in%"25-50cm",]
nodes_Pop_100[nodes_Pop_100$Depth_ind%in%"50-100cm",]

sort(unique(c(nodes_Pop_10$Phylum,
              nodes_Pop_25$Phylum,
              nodes_Pop_50$Phylum,
              nodes_Pop_100$Phylum)))

# >>> Plot Network  - OTU taxonomy ---------------------------------------------------------
PlotNetwTax <-function(tab_raph, Var, title){
  set.seed(022621)
  plot <-
    tab_raph %>%
    ggraph(layout = "fr") +
    geom_edge_link(aes(edge_colour = Direction), edge_width = 0.05) +
    geom_node_point(aes(color = Phylum, fill = Phylum, shape = Kingdom, size = get(Var))) + # log(abundance)
    geom_node_text(aes(filter = Key%in%"Module hubs" | Key%in%"Network hubs",
                       label = OTU_ID), repel=TRUE, color="black",  size =2) +
    #facet_edges(~InterKing+Direction) +
    scale_edge_color_manual(values = c("red", "black")) + 
    scale_edge_linetype_manual(values = c("solid", "dashed")) +
    scale_edge_width(range = c(0.005, 0.5))+
    scale_size(range = c(0.0005, 1)) +
    #scale_size(range = c(0.001, 1.5)) +
    scale_shape_manual(values = c(Archaea = 18, 
                                  Bacteria = 16, 
                                  Fungi = 17)) +
    # Note! If the color is set to taxonomy should specify each color 
    # for each taxa or they are going to be different acorss networks.
    scale_color_manual(values = palette_phylum) +
    theme_graph() +
    theme(plot.title = element_text(size = 10, family ="", face = "bold", hjust = 0.5),
          plot.margin = unit(c(0.05,0,0.05,0.05), "cm"), 
          axis.title.y = element_text(size = 10, family ="", face = "bold"),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 8, family ="", face = "bold", hjust = 0, vjust = 0.5), 
          legend.text = element_text(size = 7, family = ""), legend.position = "bottom",
          legend.title.align=0.5) +
    labs(title = title) +
    guides(size = FALSE, 
           shape = guide_legend(nrow = 1, order = 3, title.position="top"),
           fill = FALSE, #guide_legend(ncol = 4),
           colour = FALSE, #guide_legend(ncol=4),
           edge_colour = guide_legend(nrow = 1, order = 1, title.position="top"), 
           edge_width = FALSE)
  return(plot)
}


PlotNetwTax(ggraph_Pop_10, "AbundDepth", "0-10cm") + labs(y="Poplar")

as_ggplot(
  get_legend(PlotNetwTax(ggraph_Pop_25, "AbundDepth", "10-25cm")))


net_graph_Pop <-
  ggarrange(
    ggarrange(PlotNetwTax(ggraph_Pop_10, "AbundDepth", "0-10cm") + 
                geom_node_point(aes(filter = Key%in%"Module hubs"),shape = 1, size = 3, show.legend = F) +
                geom_node_point(aes(filter = Key%in%"Network hubs"), shape = 0, size = 3, show.legend = F), 
              PlotNetwTax(ggraph_Pop_25, "AbundDepth", "10-25cm"),
              PlotNetwTax(ggraph_Pop_50, "AbundDepth", "25-50cm") +
                geom_node_point(aes(filter = Key%in%"Module hubs"),shape = 1, size = 3, show.legend = F),
              PlotNetwTax(ggraph_Pop_100, "AbundDepth", "50-100cm") +
                geom_node_point(aes(filter = Key%in%"Module hubs"),shape = 1, size = 3, show.legend = F),
              ncol = 4,
              nrow = 1,
              #labels = c("A", "B", "C", "D"),
              widths = c(1, 1, 1, 1),
              common.legend = TRUE,
              legend = "none"),
    as_ggplot(get_legend(PlotNetwTax(ggraph_Pop_25, "AbundDepth", "10-25cm"))),
    ncol = 2,
    nrow = 1,
    widths  = c(1, 0.2))

net_graph_Pop


# >>> Plot Network  - OTU abundance --------------------------------------------------------
PlotNetwAbund <- function(tab_raph, Var, title){
  set.seed(022621)
  plot <-
    tab_raph %>%
    ggraph(layout = "fr") +
    geom_edge_link(aes(edge_colour = Direction), edge_width = 0.05) +
    geom_node_point(aes(color = RootAbund, fill = RootAbund, 
                        shape = fit_class, size = get(Var))) + # log(abundance)
    geom_node_text(aes(filter = Depth_ind%in%title,
                       label = OTU_ID), repel=TRUE, color="black",  size =2) +
    # #facet_edges(~InterKing+Direction) +
    scale_edge_color_manual(values = c("red", "black")) + 
    scale_edge_linetype_manual(values = c("solid", "dashed")) +
    scale_edge_width(range = c(0.005, 0.5))+
    scale_size(range = c(0.0005, 1)) +
    #scale_size(range = c(0.001, 1.5)) +
    scale_shape_manual(values = c(`Above prediction` = 17, 
                                  `As predicted` = 16, 
                                  `Below prediction` = 15)) +
    # scale_color_manual(values = c(Connectors = "#F8766D", 
    #                               `Module hubs` = "#7CAE00", 
    #                               `Network hubs` = "#00BFC4", 
    #                               Peripherals = "#C77CFF")) +
    # Note! If the color is set to taxonomy should specify each color 
    # for each taxa or they are going to be different acorss networks.
    theme_graph() +
    theme(plot.title = element_text(size = 10, family ="", face = "bold", hjust = 0.5),
          plot.margin = unit(c(0.4,0.01,0.4,0.01), "cm"), 
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 8, family ="", face = "bold", hjust = 0, vjust = 0.5), 
          legend.text = element_text(size = 7, family = ""), legend.position = "right") +
    labs(title = title) +
    guides(size = FALSE, 
           shape = guide_legend(ncol = 1, order = 3, title = "Neutral fit"),
           edge_colour = guide_legend(ncol = 1, order = 1), 
           edge_width = FALSE)
  return(plot)
}


PlotNetwAbund(ggraph_Pop_10, "AbundDepth", "0-10cm")


net_graph_Pop_Abund <-
  ggarrange(
    ggarrange(
      PlotNetwAbund(ggraph_Pop_10, "AbundDepth", "0-10cm") +
        geom_node_point(aes(filter = Depth_ind%in%"0-10cm"), 
                        shape = 1, size = 3, show.legend = F),
      PlotNetwAbund(ggraph_Pop_25, "AbundDepth", "10-25cm") +
        geom_node_point(aes(filter = Depth_ind%in%"10-25cm"), 
                        shape = 1, size = 3, show.legend = F),
      PlotNetwAbund(ggraph_Pop_50, "AbundDepth", "25-50cm") +
        geom_node_point(aes(filter = Depth_ind%in%"25-50cm"), 
                        shape = 1, size = 3, show.legend = F),
      PlotNetwAbund(ggraph_Pop_100, "AbundDepth", "50-100cm") +
        geom_node_point(aes(filter =Depth_ind%in%"50-100cm"), 
                        shape = 1, size = 3, show.legend = F),
      ncol = 4,
      nrow = 1,
      #labels = c("A", "B", "C", "D"),
      widths = c(1, 1, 1, 1),
      legend = "none"),
    as_ggplot(get_legend(PlotNetwAbund(ggraph_Pop_10, "AbundDepth", "0-10cm"))),
    ncol = 2,
    nrow = 1,
    widths = c(1, 0.2))

net_graph_Pop_Abund

# >>> Plot Network  - Neutral model fit ----------------------------------------------------
PlotNetwNeutral <- function(tab_raph, Var, title){
  set.seed(022621)
  plot <-
    tab_raph %>%
    ggraph(layout = "fr") +
    geom_edge_link(aes(edge_colour = Direction), edge_width = 0.05) +
    geom_node_point(aes(color = fit_class, fill = fit_class, 
                        shape = fit_class, size = get(Var))) + # log(abundance)
    geom_node_text(aes(filter = Crop_ind%in%"Poplar",
                       label = OTU_ID), repel=TRUE, color="black",  size =2) +
    # #facet_edges(~InterKing+Direction) +
    scale_edge_color_manual(values = c("red", "black")) + 
    scale_edge_linetype_manual(values = c("solid", "dashed")) +
    scale_edge_width(range = c(0.005, 0.5))+
    scale_size(range = c(0.005, 2)) +
    scale_shape_manual(values = c(`Above prediction` = 17, 
                                  `As predicted` = 16, 
                                  `Below prediction` = 15)) +
    scale_color_manual(values = c(`Above prediction` = "#F8766D", 
                                  `As predicted` = "#7CAE00", 
                                  `Below prediction` = "#00BFC4")) +
    # Note! If the color is set to taxonomy should specify each color 
    # for each taxa or they are going to be different acorss networks.
    theme_graph() +
    theme(plot.title = element_text(size = 10, family ="", face = "bold", hjust = 0.5),
          plot.margin = unit(c(0.4,0.01,0.4,0.01), "cm"), 
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 8, family ="", face = "bold", hjust = 0, vjust = 0.5), 
          legend.text = element_text(size = 7, family = ""), legend.position = "right") +
    labs(title = title) +
    guides(size = FALSE, 
           shape = guide_legend(ncol = 1, order = 3, title = "Neutral fit"),
           edge_colour = guide_legend(ncol = 1, order = 1), 
           edge_width = FALSE)
  return(plot)
}


PlotNetwNeutral(ggraph_Pop_10, "AbundCrop", "0-10cm")


net_graph_Pop_Neutral <-
  ggarrange(
    ggarrange(
      PlotNetwNeutral(ggraph_Pop_10, "AbundCrop", "0-10cm") +
        geom_node_point(aes(filter = Crop_ind%in%"Poplar"), 
                        shape = 1, size = 3, show.legend = F),
      PlotNetwNeutral(ggraph_Pop_25, "AbundCrop", "10-25cm") +
        geom_node_point(aes(filter = Crop_ind%in%"Poplar"), 
                        shape = 1, size = 3, show.legend = F),
      PlotNetwNeutral(ggraph_Pop_50, "AbundCrop", "25-50cm") +
        geom_node_point(aes(filter = Crop_ind%in%"Poplar"), 
                        shape = 1, size = 3, show.legend = F),
      PlotNetwNeutral(ggraph_Pop_100, "AbundCrop", "50-100cm") +
        geom_node_point(aes(filter =Crop_ind%in%"Poplar"), 
                        shape = 1, size = 3, show.legend = F),
      ncol = 4,
      nrow = 1,
      #labels = c("A", "B", "C", "D"),
      widths = c(1, 1, 1, 1),
      legend = "none"),
    as_ggplot(get_legend(PlotNetwNeutral(ggraph_Pop_10, "AbundCrop", "0-10cm"))),
    ncol = 2,
    nrow = 1,
    widths = c(1, 0.2))

net_graph_Pop_Neutral


# >>> Plot Network  - Interkigdom links -----------------------------------------------------
PlotNetwEdge <- function(tab_raph, Var, title){
  set.seed(022621)
  plot <-
    tab_raph %>%
    ggraph(layout = "fr") +
    geom_edge_link0(aes(edge_colour = RootV1V2_Sum), edge_width = 0.05) +
    geom_node_point(aes(shape = Kingdom, size = get(Var)), color = "grey", fill = "grey", ) + # log(abundance)
    # geom_node_text(aes(filter = Crop_ind%in%"Poplar",
    #                    label = OTU_ID), repel=TRUE, color="black",  size =2) +
    # #facet_edges(~InterKing+Direction) +
    #scale_edge_color_manual(values = c("green", "yellow")) + 
    scale_edge_colour_gradient(low = "red", high = "blue") +
    scale_edge_linetype_manual(values = c("solid", "dashed")) +
    scale_edge_width(range = c(1, 3))+
    scale_size(range = c(0.001, 1.5)) +
    scale_shape_manual(values = c(Archaea = 18,
                                  Bacteria = 16,
                                  Fungi = 17)) +
    # scale_color_manual(values = c(Connectors = "#F8766D",
    #                               `Module hubs` = "#7CAE00",
    #                               `Network hubs` = "#00BFC4",
    #                               Peripherals = "#C77CFF")) +
    # Note! If the color is set to taxonomy should specify each color 
    # for each taxa or they are going to be different acorss networks.
    theme_graph() +
    theme(plot.title = element_text(size = 10, family ="", face = "bold", hjust = 0.5),
          plot.margin = unit(c(0.4,0.01,0.4,0.01), "cm"), 
          legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.5, "cm"), 
          legend.title = element_text(size = 8, family ="", face = "bold", hjust = 0, vjust = 0.5), 
          legend.text = element_text(size = 7, family = ""), legend.position = "right") +
    labs(title = title) +
    guides(size = FALSE, 
           shape = guide_legend(ncol = 1, order = 3, title = "Neutral fit"),
           edge_colour = guide_legend(ncol = 1, order = 1), 
           edge_width = FALSE)
  return(plot)
}


PlotNetwEdge(ggraph_Pop_10, "AbundDepth", "0-10cm")


net_graph_Pop_Edge <-
  ggarrange(
    ggarrange(
      PlotNetwEdge(ggraph_Pop_10, "AbundDepth", "0-10cm"),
      # geom_node_point(aes(filter = Crop_ind%in%"Poplar"), 
      #                 shape = 1, size = 3, show.legend = F),
      PlotNetwEdge(ggraph_Pop_25, "AbundDepth", "10-25cm"),
      PlotNetwEdge(ggraph_Pop_50, "AbundDepth", "25-50cm"),
      PlotNetwEdge(ggraph_Pop_100, "AbundDepth", "50-100cm"),
      ncol = 4,
      nrow = 1,
      #labels = c("A", "B", "C", "D"),
      widths = c(1, 1, 1, 1),
      legend = "none"),
    as_ggplot(get_legend(PlotNetwEdge(ggraph_Pop_10, "AbundDepth", "0-10cm"))),
    ncol = 2,
    nrow = 1,
    widths = c(1, 0.2))

net_graph_Pop_Edge

# **** Prairie Network **** ----------------------------------------------------------------
# Extract Taxonomic Names ------------------------------------------------------------------
sort(unique(c(nodes_Pra_10$Phylum,
              nodes_Pra_25$Phylum,
              nodes_Pra_50$Phylum,
              nodes_Pra_100$Phylum)))

PlotNetwTax(ggraph_Pra_10, "AbundDepth", "0-10cm")

as_ggplot(
  get_legend(PlotNetwTax(ggraph_Pra_50, "AbundDepth", "25-50cm")))


net_graph_Pra <-
  ggarrange(
    ggarrange(PlotNetwTax(ggraph_Pra_10, "AbundDepth", "0-10cm") + 
                geom_node_point(aes(filter = Key%in%"Module hubs"),shape = 1, size = 3, show.legend = F) +
                geom_node_point(aes(filter = Key%in%"Network hubs"), shape = 0, size = 3, show.legend = F), 
              PlotNetwTax(ggraph_Pra_25, "AbundDepth", "10-25cm"),
              PlotNetwTax(ggraph_Pra_50, "AbundDepth", "25-50cm") +
                geom_node_point(aes(filter = Key%in%"Module hubs"),shape = 1, size = 3, show.legend = F),
              PlotNetwTax(ggraph_Pra_100, "AbundDepth", "50-100cm") +
                geom_node_point(aes(filter = Key%in%"Module hubs"),shape = 1, size = 3, show.legend = F),
              ncol = 4,
              nrow = 1,
              #labels = c("A", "B", "C", "D"),
              widths = c(1, 1, 1, 1),
              common.legend = TRUE,
              legend = "none"),
    as_ggplot(get_legend(PlotNetwTax(ggraph_Pra_50, "AbundDepth", "25-50cm"))),
    ncol = 2,
    nrow = 1,
    widths  = c(1, 0.2))

net_graph_Pra


# **** Switchgrass Network **** ------------------------------------------------------------
sort(unique(c(nodes_Sw_10$Phylum,
              nodes_Sw_25$Phylum,
              nodes_Sw_50$Phylum,
              nodes_Sw_100$Phylum)))

as_ggplot(
  get_legend(PlotNetwTax(ggraph_Sw_100, "AbundDepth", "50-100cm")))

net_graph_Sw <-
  ggarrange(
    ggarrange(PlotNetwTax(ggraph_Sw_10, "AbundDepth", "0-10cm") + 
                geom_node_point(aes(filter = Key%in%"Module hubs"),shape = 1, size = 3, show.legend = F) +
                geom_node_point(aes(filter = Key%in%"Network hubs"), shape = 0, size = 3, show.legend = F), 
              PlotNetwTax(ggraph_Sw_25, "AbundDepth", "10-25cm"),
              PlotNetwTax(ggraph_Sw_50, "AbundDepth", "25-50cm") +
                geom_node_point(aes(filter = Key%in%"Module hubs"),shape = 1, size = 3, show.legend = F),
              PlotNetwTax(ggraph_Sw_100, "AbundDepth", "50-100cm") +
                geom_node_point(aes(filter = Key%in%"Module hubs"),shape = 1, size = 3, show.legend = F),
              ncol = 4,
              nrow = 1,
              #labels = c("A", "B", "C", "D"),
              widths = c(1, 1, 1, 1),
              common.legend = TRUE,
              legend = "none"),
    as_ggplot(get_legend(PlotNetwTax(ggraph_Sw_25, "AbundDepth", "25-50cm"))),
    ncol = 2,
    nrow = 1,
    widths  = c(1, 0.2))

net_graph_Sw

# **** FIGURE 6A - Network plot ------------------------------------------------------------
net_graph_all <-
  ggarrange(
    ggarrange(
      PlotNetwTax(ggraph_Pop_10, "AbundDepth", "0-10cm") + 
        geom_node_point(aes(filter = Key%in%"Module hubs"),shape = 1, size = 3, show.legend = F) +
        geom_node_point(aes(filter = Key%in%"Network hubs"), shape = 0, size = 3, show.legend = F) +
        labs(y="Poplar"), 
      PlotNetwTax(ggraph_Pop_25, "AbundDepth", "10-25cm") + theme(axis.title.y = element_blank()),
      PlotNetwTax(ggraph_Pop_50, "AbundDepth", "25-50cm") +
        geom_node_point(aes(filter = Key%in%"Module hubs"),shape = 1, size = 3, show.legend = F) +
        theme(axis.title.y = element_blank()),
      PlotNetwTax(ggraph_Pop_100, "AbundDepth", "50-100cm") +
        geom_node_point(aes(filter = Key%in%"Module hubs"),shape = 1, size = 3, show.legend = F) +
        theme(axis.title.y = element_blank()),
      
      PlotNetwTax(ggraph_Pra_10, "AbundDepth", "0-10cm") + 
        geom_node_point(aes(filter = Key%in%"Module hubs"),shape = 1, size = 3, show.legend = F) +
        geom_node_point(aes(filter = Key%in%"Network hubs"), shape = 0, size = 3, show.legend = F) +
        labs(y="Prairie") + theme(plot.title = element_blank()), 
      PlotNetwTax(ggraph_Pra_25, "AbundDepth", "10-25cm")+ 
        theme(plot.title = element_blank(), axis.title.y = element_blank()),
      PlotNetwTax(ggraph_Pra_50, "AbundDepth", "25-50cm") +
        geom_node_point(aes(filter = Key%in%"Module hubs"),shape = 1, size = 3, show.legend = F)+ 
        theme(plot.title = element_blank(), axis.title.y = element_blank()),
      PlotNetwTax(ggraph_Pra_100, "AbundDepth", "50-100cm") +
        geom_node_point(aes(filter = Key%in%"Module hubs"),shape = 1, size = 3, show.legend = F) + 
        theme(plot.title = element_blank(), axis.title.y = element_blank()),
      
      PlotNetwTax(ggraph_Sw_10, "AbundDepth", "0-10cm") + 
        geom_node_point(aes(filter = Key%in%"Module hubs"),shape = 1, size = 3, show.legend = F) +
        geom_node_point(aes(filter = Key%in%"Network hubs"), shape = 0, size = 3, show.legend = F) +
        labs(y= "Switchgrass")+ theme(plot.title = element_blank()), 
      PlotNetwTax(ggraph_Sw_25, "AbundDepth", "10-25cm")+ 
        theme(plot.title = element_blank(), axis.title.y = element_blank()),
      PlotNetwTax(ggraph_Sw_50, "AbundDepth", "25-50cm") +
        geom_node_point(aes(filter = Key%in%"Module hubs"),shape = 1, size = 3, show.legend = F)+ 
        theme(plot.title = element_blank(), axis.title.y = element_blank()),
      PlotNetwTax(ggraph_Sw_100, "AbundDepth", "50-100cm") +
        geom_node_point(aes(filter = Key%in%"Module hubs"),shape = 1, size = 3, show.legend = F)+
        theme(plot.title = element_blank(), axis.title.y = element_blank()),
      ncol = 4,
      nrow = 3,
      labels = c("A","","","", "B","","","", "C","","",""),
      widths = c(1, 0.9, 0.9, 0.9),
      heights = c(1,0.85,0.85),
      common.legend = TRUE,
      legend = "none"),
    as_ggplot(get_legend(PlotNetwTax(ggraph_Pop_25, "AbundDepth", "10-25cm"))),
    ncol = 1,
    nrow = 2,
    heights = c(1, 0.1))

net_graph_all

Fig_5_net <- 
  ggarrange(net_graph_all,
            PlotLegend(),
            ncol = 2,
            nrow = 1,
            widths = c(1, 0.2))

title_net=text_grob("Microbial networks", size=10, face=2)

ggarrange(
  grid.arrange(net_graph_all, top=title_net),
  BaloonPlot(df_bar_all_rich) +
    geom_text(aes(label=Freq, x=variable), # position = position_dodge(-0.9),
              alpha=0.7, size=2, color="black") +
    labs(title = "OTU richness of network phyla", y=NULL, x=NULL),
  ncol = 2,
  nrow = 1,
  widths = c(1, 0.8),
  labels = c("", "D"))



# ********************************************----------------------------------------------
# BARPLOT OF NETWORK TAXA (not included) ---------------------------------------------------
# general plot abundance at Phylum level
source("../R_functions/ReformatTaxonomy.R")

ExtractBarDF <- function(nodes1, nodes2, nodes3, nodes4, Var){
  df1 <- nodes1[, c("OTU_ID", "AbundDepth", "Phylum")]
  df1$RelAb <-
    df1$AbundDepth/sum(df1$AbundDepth)*100
  df2 <- nodes2[, c("OTU_ID", "AbundDepth", "Phylum")]
  df2$RelAb <-
    df2$AbundDepth/sum(df2$AbundDepth)*100
  df3 <- nodes3[, c("OTU_ID", "AbundDepth", "Phylum")]
  df3$RelAb <-
    df3$AbundDepth/sum(df3$AbundDepth)*100
  df4 <- nodes4[, c("OTU_ID", "AbundDepth", "Phylum")]
  df4$RelAb <-
    df4$AbundDepth/sum(df4$AbundDepth)*100
  df <-
    full_join(
      full_join(
        full_join(df1, df2, by="OTU_ID"),
        df3, by="OTU_ID"), 
      df4, by="OTU_ID")
  print("head dataframe")
  head(df) %T>% print()
  lastValue <- function(x) tail(x[!is.na(x)], 1) # All Phyla
  last_taxon <- apply(df[, c(3,6,9,12)], 1, lastValue)
  df$BestMatch <- last_taxon
  df <-
    df[, c("OTU_ID", "BestMatch",
           "RelAb.x", "RelAb.y",
           "RelAb.x.x", "RelAb.y.y")]
  colnames(df) <-
    c("OTU_ID","Phylum", "0-10cm", "10-25cm",
      "25-50cm", "50-100cm")
  # New levels 
  df$Phylum <- 
    factor(df$Phylum, 
           levels = factor_levels)
  #df[is.na(df)] <- 0
  df <-
    melt(df)
  df$Crop <-
    rep(Var, nrow(df))
  return(df)
}

# The abudndance is standardized to the tot abund at each depth
df_bar_all <-
  rbind(ExtractBarDF(nodes_Pop_10, nodes_Pop_25, nodes_Pop_50,nodes_Pop_100, "Poplar"),
        ExtractBarDF(nodes_Pra_10, nodes_Pra_25, nodes_Pra_50,nodes_Pra_100, "Prairie"),
        ExtractBarDF(nodes_Sw_10, nodes_Sw_25, nodes_Sw_50,nodes_Sw_100, "Switchgrass"))
head(df_bar_all)
dim(df_bar_all)
unique(df_bar_all$Phylum)

PlotBar <- function(df){
  plot_bars <- 
    ggplot(df, aes(x = Phylum, y = value, color = Phylum, fill=Phylum)) + 
    geom_bar(stat = "identity", width = 0.6) +
    #ylim(0, 300) +
    theme_classic() +
    coord_flip() +
    scale_x_discrete(limits = rev(levels(df$Phylum))) +
    facet_grid(Crop~variable) +
    scale_fill_manual(values = palette_phylum) +
    scale_color_manual(values = palette_phylum) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), 
          axis.text.x = element_text(size=7),
          axis.title = element_text(angle = 0, size = 7, face = "bold"),
          axis.text.y = element_text(size = 7), 
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.25, "cm"), 
          legend.title = element_text(size = 8, face = "bold", hjust = 0, vjust = 0.5), 
          legend.text = element_text(size = 7), legend.position = "right") +
    #theme(axis.text.x = element_text(angle = 90, size = 7, vjust = 0.5, hjust = 1)) +
    guides(color=guide_legend(ncol=1),
           fill=guide_legend(ncol=1)) 
  return(plot_bars)
}

PlotBar(df_bar_all) +
  labs(title = "Relative abundance of network phyla", y=NULL, x="Phylum")

# BARPLOT OF NETWORK RICHNESS --------------------------------------------------------------
ExtractBarDF_Rich <- function(nodes1, nodes2, nodes3, nodes4, Var){
  # generate phylum df
  df1 <- as.data.frame(table(nodes1$Phylum))
  df2 <- as.data.frame(table(nodes2$Phylum))
  df3 <- as.data.frame(table(nodes3$Phylum))
  df4 <- as.data.frame(table(nodes4$Phylum))
  df <-
    full_join(
      full_join(
        full_join(df1, df2, by="Var1"),
        df3, by="Var1"), 
      df4, by="Var1")
  colnames(df) <-
    c("Phylum", "0-10cm", "10-25cm",
      "25-50cm", "50-100cm")
  df$Phylum <- 
    factor(df$Phylum, 
           levels = factor_levels)
  df <- melt(df)
  df$Crop <-
    rep(Var, nrow(df))
  # rename variables 
  names(nodes1)[6] <- "Depth"
  names(nodes2)[6] <- "Depth"
  names(nodes3)[6] <- "Depth"
  names(nodes4)[6] <- "Depth"
  # extract hubs values
  df_hub <-
    rbind(nodes1[nodes1$Key == "Network hubs" |
                   nodes1$Key == "Module hubs", ][, c("OTU_ID","Phylum", "Depth", "Key")],
          nodes2[nodes2$Key == "Network hubs" |
                   nodes2$Key == "Module hubs", ][, c("OTU_ID","Phylum", "Depth", "Key")],
          nodes3[nodes3$Key == "Network hubs" |
                   nodes3$Key == "Module hubs", ][, c("OTU_ID","Phylum", "Depth", "Key")],
          nodes4[nodes4$Key == "Network hubs" |
                   nodes4$Key == "Module hubs", ][, c("OTU_ID","Phylum", "Depth", "Key")])
  df_hub$Crop <-
    rep(Var, nrow(df_hub))
  df_hub$Hub <-
    paste(df_hub$Phylum, df_hub$Depth, df_hub$Crop, sep = "_")
  df$Hub <- paste(df$Phylum, df$variable, df$Crop, sep = "_")
  df_num <- as.data.frame(table(df_hub$Hub))
  colnames(df_num) <- c("Hub", "Freq")
  df_all <- 
    df_hub[!duplicated(df_hub$Hub), ]
  df_all <- 
    left_join(df_all, df_num, by="Hub")
  head(df_all) %T>% print()
  df_all <-
    full_join(df, df_all[,c(1,6:7)], by="Hub")
  return(df_all)
}

ExtractBarDF_Rich(nodes_Pop_10, nodes_Pop_25, nodes_Pop_50,nodes_Pop_100, "Poplar") 

# Richness at each Phylum
df_bar_all_rich <-
  rbind(ExtractBarDF_Rich(nodes_Pop_10, nodes_Pop_25, nodes_Pop_50,nodes_Pop_100, "Poplar"),
        ExtractBarDF_Rich(nodes_Pra_10, nodes_Pra_25, nodes_Pra_50,nodes_Pra_100, "Prairie"),
        ExtractBarDF_Rich(nodes_Sw_10, nodes_Sw_25, nodes_Sw_50,nodes_Sw_100, "Switchgrass"))
head(df_bar_all_rich)
dim(df_bar_all_rich)
unique(df_bar_all_rich$Phylum)

df_bar_all_rich$variable <- 
  factor(df_bar_all_rich$variable, 
         levels = c("0-10cm", "10-25cm",
                    "25-50cm", "50-100cm"))

PlotBar(df_bar_all_rich) +
  labs(title = "OTU richness of network phyla", y=NULL, x="Phylum")

# ********************************************----------------------------------------------
# BALOON PLOTS -----------------------------------------------------------------------------
# Plot Relative Abundance or Richness with baloons
BaloonPlot <- function(df, Var){
  plot_bal <-
    ggplot(df, aes(x=fct_rev(variable), y=Phylum, color=Phylum, size=value)) + #alpha=value
    geom_point() +
    #scale_alpha_continuous(range=c(0.6, 0.9)) +
    scale_x_discrete(limits = rev) +
    scale_y_discrete(limits = rev) +
    facet_grid(~Crop) +
    scale_color_manual(values = palette_phylum) +
    scale_size_continuous(name = Var,
                          breaks = c(1, 25, 50, 75, 100),
                          range = c(1, 6)) +
    theme_classic() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), 
          axis.text.x = element_text(angle=45, size=7, hjust = 1, vjust = 1.1),
          axis.title = element_text(angle = 0, size = 7, face = "bold"),
          axis.text.y = element_text(size = 7), 
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.25, "cm"), 
          legend.title = element_text(size = 8, face = "bold", hjust = 0, vjust = 0.5), 
          legend.text = element_text(size = 7), legend.position = "right") +
    guides(color=guide_legend(ncol=1, override.aes = list(size=2)),
           size=guide_legend(ncol=1, order = 1)) 
  return(plot_bal)
}


# **** FIGURE 6D - Network Baloon plots ----------------------------------------------------
BaloonPlot(df_bar_all_rich, "Richness") +  
  geom_text(aes(label=Freq, x=variable), # position = position_dodge(-0.9),
            alpha=1.0, size=2, color="black") +
  labs(title = "OTU richness of network phyla", y=NULL, x="Phylum")

# I have to remove NA and summarize Phyla before plotting abundance
df_bar_all %>%
  drop_na() %>%
  group_by(Phylum, variable, Crop) %>%
  summarise(value = round(sum(value), 1)) %>%
  BaloonPlot(df = ., "Richness") +
  geom_text(aes(label=value, x=variable), # position = position_dodge(-0.9),
            alpha=1.0, size=2, color="black") +
  labs(title = "Relative abundance of network phyla", y=NULL, x="Phylum")


# ********************************************----------------------------------------------
# EDGES HEATMAP ----------------------------------------------------------------------------
# Phylum-to-Phylum edges 
head(nodes_Pop_10)
head(edges_Pop_10)

# Generate DF for heatmap at Phylum level --------------------------------------------------
HeatEdge <- function(df_edges, Var, Depth){
  df_edges <-
    df_edges %>%
    mutate(ModuleConnectivity = ifelse(V1_Phylum==V2_Phylum, "in-taxon", "out-taxon"))
  df_count <-
    dplyr::count(df_edges, 
                 vars= Direction, V1_Phylum, V2_Phylum, ModuleConnectivity) 
  df_count_filt <-
    subset(df_count, vars%in%Var)
  df_count_filt <-
    df_count_filt[, c(2,3,5)]
  print("Subsetted dataframe dimensions")
  dim(df_count_filt) %T>% print()
  # extract Class data
  df_count_filt$V1_Phylum <- as.character(df_count_filt$V1_Phylum)
  df_count_filt$V2_Phylum <- as.character(df_count_filt$V2_Phylum)
  faclist <- vector("list", nrow(df_count_filt))
  for (i in 1:nrow(df_count_filt)) {
    faclist[[i]] <-
      sort(c(df_count_filt$V1_Phylum[i], df_count_filt$V2_Phylum[i]))
    faclist[[i]] <- paste(faclist[[i]][1], faclist[[i]][2], sep = "_")
  }
  df_count_filt$v1v2 <- unlist(faclist)
  print(df_count_filt)
  df_count_filt %>% 
    dplyr::group_by(v1v2) %>% 
    dplyr::summarize(Sum=sum(n)) %>%
    as.data.frame() -> df_new
  print("Final dataframe dimensions")
  dim(df_new)  %T>% print()
  df_final <-
    left_join(df_new, df_count_filt[,c(1:2,4)], "v1v2")
  print("Final dataframe")
  head(df_final)  %T>% print()
  df_final$Depth <- rep(Depth, times=nrow(df_final))
  print("Phyla involved")
  df_final <-
    df_final[!duplicated(df_final$v1v2), ]
  unique(c(as.character(df_final$V1_Phylum), 
           as.character(df_final$V2_Phylum))) %T>%
    print()
  return(df_final)
}


HeatEdge(edges_Pop_10, "positive", "0-10cm")
HeatEdge(edges_Pop_10, "negative", "0-10cm")

df_Pop_positive <-
  rbind(
    HeatEdge(edges_Pop_10, "positive", "0-10cm"),
    HeatEdge(edges_Pop_25, "positive", "10-25cm"),
    HeatEdge(edges_Pop_50, "positive", "25-50cm"),
    HeatEdge(edges_Pop_100, "positive", "50-100cm")) 

df_Pop_negative <-
  rbind(
    HeatEdge(edges_Pop_10, "negative", "0-10cm"),
    HeatEdge(edges_Pop_25, "negative", "10-25cm"),
    HeatEdge(edges_Pop_50, "negative", "25-50cm"),
    HeatEdge(edges_Pop_100, "negative", "50-100cm")) 

df_Pop_all <-
  rbind(data.frame(Direction = rep("positive", nrow(df_Pop_positive)),
                   df_Pop_positive),
        data.frame(Direction = rep("negative", nrow(df_Pop_negative)), 
                   df_Pop_negative))


# Plotting Edge Heatmap --------------------------------------------------------------------
PlotHeat <- function(dataframe){
  heat_plot <-
    ggplot(dataframe, aes(x=factor(V1_Phylum, levels = levels), 
                          y=factor(V2_Phylum, levels = levels), fill=Sum)) +
    geom_point()+ # I have to trick R with geom point
    geom_rect(aes(xmin=0, xmax=9.5, ymin=0, ymax=10.5),fill="wheat", alpha=0.2) + # for all phyla heatmap
    annotate("rect", xmin=0, xmax=9.5, ymin=0, ymax=10.5, alpha=0.4, fill="wheat") +
    annotate("rect", xmin=9.5, xmax=25.5, ymin=10.5, ymax=25.5, alpha=0.4, fill="wheat") +
    geom_tile(aes(height = 0.92, width = 0.92)) + 
    geom_text(aes(label = Sum), size=2, color="white") +
    theme_classic() +
    facet_grid(~Direction~Depth) +
    labs(title = "Within- and between-Class connections") +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), 
          axis.text.x =  element_text(size = 6, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 6), 
          axis.title = element_blank(),
          legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"), 
          legend.title = element_text("Edges", size = 9, face = "bold",  
                                      hjust = 0, vjust = 0.5), 
          legend.text = element_text(size = 9), legend.position = "right") +
    #guides(fill = guide_legend(title = "Edges")) +
    scale_fill_gradient2("Edges", low="blue", mid="grey", high="red", #colors in the scale
                         midpoint=max(dataframe$Sum)/2,   #same midpoint for plots (mean of the range)
                         limits=c(0, max(dataframe$Sum))) # 70 for positive
  return(heat_plot)
}


PlotHeat(df_Pop_all) 
#geom_hline(yintercept=10.5, linetype="dashed", color = "black", size=0.4) +
#geom_vline(xintercept=9.5, linetype="dashed", color = "black", size=0.4) 


df_Pra_positive <-
  rbind(
    HeatEdge(edges_Pra_10, "positive", "0-10cm"),
    HeatEdge(edges_Pra_25, "positive", "10-25cm"),
    HeatEdge(edges_Pra_50, "positive", "25-50cm"),
    HeatEdge(edges_Pra_100, "positive", "50-100cm")) 

df_Pra_negative <-
  rbind(
    HeatEdge(edges_Pra_10, "negative", "0-10cm"),
    HeatEdge(edges_Pra_25, "negative", "10-25cm"),
    HeatEdge(edges_Pra_50, "negative", "25-50cm"),
    HeatEdge(edges_Pra_100, "negative", "50-100cm")) 

df_Pra_all <-
  rbind(data.frame(Direction = rep("positive", nrow(df_Pra_positive)),
                   df_Pra_positive),
        data.frame(Direction = rep("negative", nrow(df_Pra_negative)), 
                   df_Pra_negative))

PlotHeat(df_Pra_all) 


df_Sw_positive <-
  rbind(
    HeatEdge(edges_Sw_10, "positive", "0-10cm"),
    HeatEdge(edges_Sw_25, "positive", "10-25cm"),
    HeatEdge(edges_Sw_50, "positive", "25-50cm"),
    HeatEdge(edges_Sw_100, "positive", "50-100cm")) 

df_Sw_negative <-
  rbind(
    HeatEdge(edges_Sw_10, "negative", "0-10cm"),
    HeatEdge(edges_Sw_25, "negative", "10-25cm"),
    HeatEdge(edges_Sw_50, "negative", "25-50cm"),
    HeatEdge(edges_Sw_100, "negative", "50-100cm")) 

df_Sw_all <-
  rbind(data.frame(Direction = rep("positive", nrow(df_Sw_positive)),
                   df_Sw_positive),
        data.frame(Direction = rep("negative", nrow(df_Sw_negative)), 
                   df_Sw_negative))

PlotHeat(df_Sw_all) 

# Rank most connected families -------------------------------------------------------------
df_Pop_all %>%
  group_by(v1v2) %>%
  summarize(Tot = sum(Sum)) %>%
  arrange(desc(Tot)) %>%
  as.data.frame()

df_Pra_all %>%
  group_by(v1v2) %>%
  summarize(Tot = sum(Sum)) %>%
  arrange(desc(Tot)) %>%
  as.data.frame()

df_Sw_all %>%
  group_by(v1v2) %>%
  summarize(Tot = sum(Sum)) %>%
  arrange(desc(Tot)) %>%
  as.data.frame()


# Summarize to most important Phyla --------------------------------------------------------
ExtractPhyla <- function(dataframe){
  dataframe_filt <-
    subset(dataframe, 
           V1_Phylum%in%c("Proteobacteria", "Actinobacteriota",
                          "Acidobacteriota", "Chloroflexi", 
                          "Ascomycota", "Basidiomycota") &
             V2_Phylum%in%c("Proteobacteria", "Actinobacteriota",
                            "Acidobacteriota", "Chloroflexi", 
                            "Ascomycota", "Basidiomycota"))
  dataframe_filt$V1_Phylum <-
    factor(dataframe_filt$V1_Phylum, 
           levels=c("Ascomycota", "Basidiomycota",
                    "Actinobacteriota","Chloroflexi",
                    "Acidobacteriota","Proteobacteria"))
  dataframe_filt$V2_Phylum <-
    factor(dataframe_filt$V2_Phylum, 
           levels=c("Ascomycota", "Basidiomycota",
                    "Actinobacteriota","Chloroflexi",
                    "Acidobacteriota","Proteobacteria"))
  melt_df <-
    dataframe %>%
    group_by(Depth) %>%
    summarise(All=sum(Sum)) %>%
    melt()
  head(melt_df) %T>% print()
  dataframe_filt <-
    inner_join(dataframe_filt, melt_df[-2], by="Depth")
  dataframe_filt$Mean <- 
    round(dataframe_filt[,"Sum"]/dataframe_filt$value*100, digits = 1)
  return(dataframe_filt)
}

ExtractPhyla(df_Pop_all) -> df_Pop_all_filt
df_Pop_all_filt
df_Pop_all[order(df_Pop_all$Sum),]

PlotHeatPhy <- function(dataframe){
  heat_plot <-
    ggplot(dataframe, aes(x=V1_Phylum, y=V2_Phylum, fill=Mean)) +
    geom_tile(aes(height = 0.92, width = 0.92)) + 
    theme_classic() +
    facet_grid(~Direction~Depth) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
          plot.margin = unit(c(0.1,0.3,0.1,0.1), "cm"), 
          axis.text.x =  element_text(size = 7, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 7), 
          axis.title = element_blank(),
          legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text("Edges", size = 7, face = "bold",  
                                      hjust = 0, vjust = 0.5), 
          legend.text = element_text(size = 8), legend.position = "right") 
  #guides(fill = guide_legend(title = "Edges")) +
  return(heat_plot)
}

# The % of Edges is calculated on the total edges of a network/depth
PlotHeatPhy(df_Pop_all_filt) +
  geom_text(aes(label = Mean), size=2, color="white") +
  scale_fill_gradient2("Edges %", low="blue", mid="grey", high="red", 
                       midpoint=max(df_Pop_all_filt$Mean)/2,  
                       limits=c(0, max(df_Pop_all_filt$Mean))) 


# **** FIGURE S8 - heatmap with to connected Phyla  ----------------------------------------
heat_Phyl_all <-
  ggarrange(PlotHeatPhy(ExtractPhyla(df_Pop_all)) +
              geom_text(aes(label = Mean), size=2, color="white") +
              scale_fill_gradient2("Links %", low="blue", mid="grey", high="red", 
                                   midpoint=max(df_Pop_all_filt$Mean)/2,  
                                   limits=c(0, max(df_Pop_all_filt$Mean))) +
              theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
              labs(title = "Poplar"), 
            PlotHeatPhy(ExtractPhyla(df_Pra_all)) +
              geom_text(aes(label = Mean), size=2, color="white") +
              scale_fill_gradient2("Links %", low="blue", mid="grey", high="red", 
                                   midpoint=max(ExtractPhyla(df_Pra_all)$Mean)/2,  
                                   limits=c(0, max(ExtractPhyla(df_Pra_all)$Mean))) +
              theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
              labs(title = "Prairie"), 
            PlotHeatPhy(ExtractPhyla(df_Sw_all)) +
              geom_text(aes(label = Mean), size=2, color="white") +
              scale_fill_gradient2("Links %", low="blue", mid="grey", high="red", 
                                   midpoint=max(ExtractPhyla(df_Sw_all)$Mean)/2,  
                                   limits=c(0, max(ExtractPhyla(df_Sw_all)$Mean))) +
              labs(title = "Switchgrass"), 
            ncol = 1,
            nrow = 3,
            labels = c("A", "B", "C"),
            heights = c(0.74, 0.74, 0.96),
            legend = "right")

heat_Phyl_all

title_heat=text_grob("Intra- and inter-phylum edges", size=12, face=2)
grid.arrange(heat_Phyl_all, top=title_heat)

# Root-to-Soil Phyla connections -----------------------------------------------------------
HeatEdgeAbund <- function(df_edges, Var, Depth){
  df_count <-
    df_edges %>% 
    dplyr::group_by(Direction, V1_Phylum, V2_Phylum) %>% 
    dplyr::summarize(n=round(sum(RootV1V2_Sum), 1)) %>%
    as.data.frame() 
  print(df_count)
  df_count_filt <-
    subset(df_count, Direction%in%Var)
  # df_count_filt <-
  #   subset(df_count, vars%in%Var)
  # df_count_filt <-
  #   df_count_filt[, c(2,3,5)]
  print("Subsetted dataframe dimensions")
  head(df_count_filt) %T>% print()
  # extract Class data
  df_count_filt$V1_Phylum <- as.character(df_count_filt$V1_Phylum)
  df_count_filt$V2_Phylum <- as.character(df_count_filt$V2_Phylum)
  faclist <- vector("list", nrow(df_count_filt))
  for (i in 1:nrow(df_count_filt)) {
    faclist[[i]] <-
      sort(c(df_count_filt$V1_Phylum[i], df_count_filt$V2_Phylum[i]))
    faclist[[i]] <- paste(faclist[[i]][1], faclist[[i]][2], sep = "_")
  }
  df_count_filt$v1v2 <- unlist(faclist)
  df_count_filt %>% 
    dplyr::group_by(v1v2) %>% 
    dplyr::summarize(Sum=sum(n)) %>%
    as.data.frame() -> df_new
  print("Final dataframe dimensions")
  head(df_new)  %T>% print()
  df_final <-
    left_join(df_new, df_count_filt, "v1v2")
  print("Final dataframe")
  head(df_final)  %T>% print()
  df_final <-
    df_final[!duplicated(df_final$v1v2), ]
  df_final <-
    df_final[, c(3:5,2)]
  df_final$Depth <- rep(Depth, times=nrow(df_final))
  print("Phyla involved")
  unique(c(as.character(df_final$V1_Phylum), 
           as.character(df_final$V2_Phylum))) %T>%
    print()
  return(df_final)
}


HeatEdgeAbund(edges_Pop_10, "positive", "0-10cm")

df_Pop_positive_abund <-
  rbind(
    HeatEdgeAbund(edges_Pop_10, "positive", "0-10cm"),
    HeatEdgeAbund(edges_Pop_25, "positive", "10-25cm"),
    HeatEdgeAbund(edges_Pop_50, "positive", "25-50cm"),
    HeatEdgeAbund(edges_Pop_100, "positive", "50-100cm")) 

df_Pop_negative_abund <-
  rbind(
    HeatEdgeAbund(edges_Pop_10, "negative", "0-10cm"),
    HeatEdgeAbund(edges_Pop_25, "negative", "10-25cm"),
    HeatEdgeAbund(edges_Pop_50, "negative", "25-50cm"),
    HeatEdgeAbund(edges_Pop_100, "negative", "50-100cm")) 

df_Pop_all_abund <-
  rbind(df_Pop_positive_abund, df_Pop_negative_abund)

PlotHeat(df_Pop_all_abund) 

df_Pop_all_abund <-rbind(
  HeatEdgeAbund(edges_Pop_10, "positive", "0-10cm"),
  HeatEdgeAbund(edges_Pop_25, "positive", "10-25cm"),
  HeatEdgeAbund(edges_Pop_50, "positive", "25-50cm"),
  HeatEdgeAbund(edges_Pop_100, "positive", "50-100cm"),
  HeatEdgeAbund(edges_Pop_10, "negative", "0-10cm"),
  HeatEdgeAbund(edges_Pop_25, "negative", "10-25cm"),
  HeatEdgeAbund(edges_Pop_50, "negative", "25-50cm"),
  HeatEdgeAbund(edges_Pop_100, "negative", "50-100cm")) 


df_Pra_all_abund <-rbind(
  HeatEdgeAbund(edges_Pra_10, "positive", "0-10cm"),
  HeatEdgeAbund(edges_Pra_25, "positive", "10-25cm"),
  HeatEdgeAbund(edges_Pra_50, "positive", "25-50cm"),
  HeatEdgeAbund(edges_Pra_100, "positive", "50-100cm"),
  HeatEdgeAbund(edges_Pra_10, "negative", "0-10cm"),
  HeatEdgeAbund(edges_Pra_25, "negative", "10-25cm"),
  HeatEdgeAbund(edges_Pra_50, "negative", "25-50cm"),
  HeatEdgeAbund(edges_Pra_100, "negative", "50-100cm")) 


df_sw_all_abund <-rbind(
  HeatEdgeAbund(edges_Sw_10, "positive", "0-10cm"),
  HeatEdgeAbund(edges_Sw_25, "positive", "10-25cm"),
  HeatEdgeAbund(edges_Sw_50, "positive", "25-50cm"),
  HeatEdgeAbund(edges_Sw_100, "positive", "50-100cm"),
  HeatEdgeAbund(edges_Sw_10, "negative", "0-10cm"),
  HeatEdgeAbund(edges_Sw_25, "negative", "10-25cm"),
  HeatEdgeAbund(edges_Sw_50, "negative", "25-50cm"),
  HeatEdgeAbund(edges_Sw_100, "negative", "50-100cm")) 


# Extarct most interesting Phyla
ExtractPhyla(df_Pop_all_abund) -> df_Pop_all_filt_abund
head(df_Pop_all_filt_abund)

df_Pop_all_filt_abund

# **** FIGURE S7 - root-to-root Phyla  ------------------------------------------------------
PlotHeatPhy(df_Pop_all_filt_abund) +
  geom_text(aes(label = Mean), size=2.5, color="black") +
  scale_fill_gradient2("Abundance %", low="blue", mid="grey", high="yellow", 
                       midpoint=max(df_Pop_all_filt_abund$Mean)/2,  
                       limits=c(0, max(df_Pop_all_filt_abund$Mean))) +
  labs(title = "Intra- and Inter-Phyla Root-to-Soil edges")

# **** FIGURE S9 - root-to-root Phyla  -----------------------------------------------------
# The idea is that some root important function performed by
# root-associated taxa may be mediated by soil taxa. The mean
# is calculated as the sum of the root abundances in a class
# divided by the total root abundance for that specific network/depth

df_Pop_all_abund$Direction <- factor(df_Pop_all_abund$Direction, levels = c("positive", "negative"))
df_Pra_all_abund$Direction <- factor(df_Pra_all_abund$Direction, levels = c("positive", "negative"))
df_sw_all_abund$Direction <- factor(df_sw_all_abund$Direction, levels = c("positive", "negative"))

heat_Phyl_all_root <-
  ggarrange(PlotHeatPhy(ExtractPhyla(df_Pop_all_abund)) +
              geom_text(aes(label = Mean), size=2, color="black") +
              scale_fill_gradient2("Abund. %", low="blue", mid="grey", high="yellow", 
                                   midpoint=max(df_Pop_all_filt_abund$Mean)/2,  
                                   limits=c(0, max(df_Pop_all_filt_abund$Mean))) +
              theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
              labs(title = "Poplar"), 
            PlotHeatPhy(ExtractPhyla(df_Pra_all_abund)) +
              geom_text(aes(label = Mean), size=2, color="black") +
              scale_fill_gradient2("Abund. %", low="blue", mid="grey", high="yellow", 
                                   midpoint=max(ExtractPhyla(df_Pra_all_abund)$Mean)/2,  
                                   limits=c(0, max(ExtractPhyla(df_Pra_all_abund)$Mean))) +
              theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
              labs(title = "Prairie"), 
            PlotHeatPhy(ExtractPhyla(df_sw_all_abund)) +
              geom_text(aes(label = Mean), size=2, color="black") +
              scale_fill_gradient2("Abund. %", low="blue", mid="grey", high="yellow", 
                                   midpoint=max(ExtractPhyla(df_sw_all_abund)$Mean)/2,  
                                   limits=c(0, max(ExtractPhyla(df_sw_all_abund)$Mean))) +
              labs(title = "Switchgrass"), 
            ncol = 1,
            nrow = 3,
            labels = c("A", "B", "C"),
            heights = c(0.74, 0.74, 0.96),
            legend = "right")

heat_Phyl_all_root

title_heat_abund=text_grob("Intra- and Inter-Kingdom root-to-root links", size=12, face=2)
grid.arrange(heat_Phyl_all_root, top=title_heat_abund)


# Summarize by just Kingdom ----------------------------------------------------------------
SummInterKing <- function(edges, Depth){
  df <-
    edges %>%
    dplyr::count(vars= Direction, V1_taxa, V2_taxa) 
  df$Mean <- 
    round(df$n/nrow(edges)*100, digits = 1)
  df$Depth <- 
    rep(Depth, times=nrow(df))
  colnames(df)[1] <- "Direction"
  return(df)
}

SummInterKing(edges_Pop_10, "0-10cm")

df_Pop_king <-
  rbind(
    SummInterKing(edges_Pop_10, "0-10cm"),
    SummInterKing(edges_Pop_25, "10-25cm"),
    SummInterKing(edges_Pop_50, "25-50cm"),
    SummInterKing(edges_Pop_100, "50-100cm")) 

df_Pop_king$Crop <-
  rep("Poplar", nrow(df_Pop_king))

df_Pra_king <-
  rbind(
    SummInterKing(edges_Pra_10, "0-10cm"),
    SummInterKing(edges_Pra_25, "10-25cm"),
    SummInterKing(edges_Pra_50, "25-50cm"),
    SummInterKing(edges_Pra_100, "50-100cm")) 
df_Pra_king$Crop <-
  rep("Prairie", nrow(df_Pra_king))

df_Sw_king <-
  rbind(
    SummInterKing(edges_Sw_10, "0-10cm"),
    SummInterKing(edges_Sw_25, "10-25cm"),
    SummInterKing(edges_Sw_50, "25-50cm"),
    SummInterKing(edges_Sw_100, "50-100cm")) 
df_Sw_king$Crop <-
  rep("Switchgrass", nrow(df_Sw_king))

df_king_all <-
  rbind(df_Pop_king, df_Pra_king, df_Sw_king)


PlotHeatKing <- function(dataframe){
  heat_plot <-
    ggplot(dataframe, aes(x=V1_taxa, y=V2_taxa, fill=Mean)) +
    geom_tile(aes(height = 0.92, width = 0.92)) + 
    theme_classic() +
    facet_grid(Crop~Direction+Depth) +
    labs(title = "Intra- and Inter-Kingdom edges") +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.margin = unit(c(0.1,0,0.1,0.1), "cm"), 
          axis.text.x =  element_text(size = 7, angle = 45, vjust = 1.2, hjust = 1.1),
          axis.text.y = element_text(size = 7), 
          axis.title = element_blank(),
          strip.text.x = element_text(size = 7), 
          strip.text.y = element_text(size = 7),
          legend.key.height = unit(0.3, "cm"), #legend.key.width = unit(0.3, "cm"), 
          legend.margin = margin(0,0,0,0, "cm"),
          legend.title = element_text(size = 8, face = "bold",  
                                      hjust = 0, vjust = 0.5), 
          legend.text = element_text(size = 8), legend.position = "bottom") 
  return(heat_plot)
}

# **** FIGURE 7A - inta- inter-kingdom links -----------------------------------------------
PlotHeatKing(df_king_all) +
  geom_text(aes(label = Mean), size=2.5, color="white") +
  scale_fill_gradient2("Edges No.", low="blue", mid="grey", high="red", 
                       midpoint=max(df_king_all$Mean)/2,  
                       limits=c(0, max(df_king_all$Mean))) 


# ****FIGURE 6 - Baloon/Heatmap ------------------------------------------------------------
Fig_6_balon_heat_all <-
  ggarrange(
    PlotHeatKing(df_king_all) +
      geom_text(aes(label = Mean), size=2, color="white") +
      scale_fill_gradient2("Edges No.", low="blue", mid="grey", high="red", 
                           midpoint=max(df_king_all$Mean)/2,  
                           limits=c(0, max(df_king_all$Mean))) +
      labs(title = "Intra- and inter-phylum links"),
    BaloonPlot(df_bar_all_rich, "Richness") +
      geom_text(aes(label=Freq, x=variable), # position = position_dodge(-0.9),
                alpha=0.7, size=2, color="black") +
      labs(title = "OTU richness of network phyla", y=NULL, x=NULL),
    ncol = 2,
    nrow = 1,
    widths = c(1,1.5),
    labels = c("A", "B"))

Fig_6_balon_heat_all

# Root-to-Soil Kingdom connections ---------------------------------------------------------
# The higher the RooAbund the more connections happen
# between taxa abundand in Root comapre to soil.
SummRootToSoil <- function(edges, Depth){
  df <-
    edges %>% 
    dplyr::group_by(Direction, V1_taxa, V2_taxa) %>% 
    dplyr::summarize(n=round(sum(RootV1V2_Sum), 0)) %>%
    as.data.frame() 
  df$Mean <-
    round(df$n/sum(df$n) * 100,1)
  df$Depth <- 
    rep(Depth, times=nrow(df))
  return(df)
}

SummRootToSoil(edges_Pop_10, "0-10cm")

df_Pop_root_abund <-
  rbind(
    SummRootToSoil(edges_Pop_10, "0-10cm"),
    SummRootToSoil(edges_Pop_25, "10-25cm"),
    SummRootToSoil(edges_Pop_50, "25-50cm"),
    SummRootToSoil(edges_Pop_100, "50-100cm")) 
df_Pop_root_abund$Crop <-
  rep("Poplar", nrow(df_Pop_root_abund))

df_Pra_root_abund <-
  rbind(
    SummRootToSoil(edges_Pra_10, "0-10cm"),
    SummRootToSoil(edges_Pra_25, "10-25cm"),
    SummRootToSoil(edges_Pra_50, "25-50cm"),
    SummRootToSoil(edges_Pra_100, "50-100cm")) 
df_Pra_root_abund$Crop <-
  rep("Prairie", nrow(df_Pra_root_abund))

df_Sw_root_abund <-
  rbind(
    SummRootToSoil(edges_Sw_10, "0-10cm"),
    SummRootToSoil(edges_Sw_25, "10-25cm"),
    SummRootToSoil(edges_Sw_50, "25-50cm"),
    SummRootToSoil(edges_Sw_100, "50-100cm")) 
df_Sw_root_abund$Crop <-
  rep("Switchgrass", nrow(df_Sw_root_abund))

df_root_abund_all <-
  rbind(df_Pop_root_abund, df_Pra_root_abund, df_Sw_root_abund)

# **** FIGURE S7 - inta- and inter-kingdom root-to-root links ------------------------------
PlotHeatKing(df_root_abund_all) +
  geom_text(aes(label = Mean), size=2.5, color="black") +
  scale_fill_gradient2("Abund. %", low="blue", mid="grey", high="yellow", 
                       midpoint=max(df_root_abund_all$Mean)/2,  
                       limits=c(0, max(df_root_abund_all$Mean))) +
  labs(title = "Intra- and Inter-Kingdom root-to-root links")

# root-to-soil heatmap plus baloons
Fig_7_balon_heat_abund_all <-
  ggarrange(
    PlotHeatKing(df_root_abund_all) +
      geom_text(aes(label = Mean), size=2, color="black") +
      scale_fill_gradient2("Links %", low="blue", mid="grey", high="yellow", 
                           midpoint=max(df_root_abund_all$Mean)/2,  
                           limits=c(0, max(df_root_abund_all$Mean))) +
      labs(title = "Intra- and inter-phylum root-to-soil links"),
    df_bar_all %>%
      drop_na() %>%
      group_by(Phylum, variable, Crop) %>%
      summarise(value = round(sum(value), 1)) %>%
      BaloonPlot(df = ., Var = "Mean") +
      geom_text(aes(label=value, x=variable), # position = position_dodge(-0.9),
                alpha=0.7, size=2, color="black") +
      labs(title = "Relative abundance of network phyla", y=NULL, x="Phylum"),
    ncol = 2,
    nrow = 1,
    widths = c(1,1.5),
    labels = c("A", "B"))

Fig_7_balon_heat_abund_all


# **** FIGURE S6 - baloon relative abundance -----------------------------------------------
df_bar_all %>%
  drop_na() %>%
  group_by(Phylum, variable, Crop) %>%
  summarise(value = round(sum(value), 1)) %>%
  BaloonPlot("Abundance %") +
  geom_text(aes(label=value, x=variable), # position = position_dodge(-0.9),
            alpha=0.7, size=2, color="black") +
  labs(title = "Relative abundance of network phyla", y=NULL, x="Phylum")

# HUBS - INDICATORS OVERLAP ----------------------------------------------------------------
library(ggVennDiagram)

KeyOverlap <- function(nodes, Niche, Crop, Depth){
  lista <-
    list(#Hubs = rownames(nodes[nodes$Key%in%c("Connectors", 
      #                                     "Module hubs",
      #                                     "Network hubs"), ]),
      Fit = rownames(nodes[nodes$fit_class%in%"Above prediction", ]),
      #Root = rownames(nodes[nodes$Niche_ind%in%Niche, ]),
      Depth = rownames(nodes[nodes$Depth_ind%in%Depth, ]), 
      Crop =rownames(nodes[nodes$Crop_ind%in%Crop, ]))
  return(lista)
}

KeyOverlap <- function(nodes, Niche, Crop, Depth){
  lista <-
    list(Mod = rownames(nodes[nodes$Key%in%"Module hubs", ]),
         Net = rownames(nodes[nodes$Key%in%"Network hubs", ]))
  #Con = rownames(nodes[nodes$Key%in%"Connectors", ]),
  #Per = rownames(nodes[nodes$Key%in%"Peripherals", ]))
  return(lista)
}

KeyOverlap(nodes_Pop_10, "Root", "Poplar", "0-10cm")

ggVennDiagram(KeyOverlap(nodes_Pop_10, "Root", "Poplar", "0-10cm"), label_alpha = 0) 
ggVennDiagram(KeyOverlap(nodes_Pop_25, "Root", "Poplar", "10-25cm"), label_alpha = 0)
ggVennDiagram(KeyOverlap(nodes_Pop_50, "Root", "Poplar", "25-50cm"), label_alpha = 0)
ggVennDiagram(KeyOverlap(nodes_Pop_100, "Root", "Poplar", "50-100cm"), label_alpha = 0)


# ********************************************----------------------------------------------
# GENERAL STATS for ALL networks -----------------------------------------------------------
# Indicators-Hub plot ----------------------------------------------------------------------
head(nodes_Pop_10)

list(rownames(nodes_Pop_10[nodes_Pop_10$Key%in%"Module hubs", ]),
     rownames(nodes_Pop_25[nodes_Pop_25$Key%in%"Module hubs", ]),
     rownames(nodes_Pop_50[nodes_Pop_50$Key%in%"Module hubs", ]),
     rownames(nodes_Pop_100[nodes_Pop_100$Key%in%"Module hubs", ]),
     rownames(nodes_Pop_10[nodes_Pop_10$Key%in%"Network hubs", ]),
     rownames(nodes_Pop_25[nodes_Pop_25$Key%in%"Network hubs", ]),
     rownames(nodes_Pop_50[nodes_Pop_50$Key%in%"Network hubs", ]),
     rownames(nodes_Pop_100[nodes_Pop_100$Key%in%"Network hubs", ]))

list(Depth10 = rownames(nodes_Pop_10[nodes_Pop_10$Key%in%"Connectors", ]),
     Depth25 = rownames(nodes_Pop_25[nodes_Pop_25$Key%in%"Connectors", ]),
     Depth50 = rownames(nodes_Pop_50[nodes_Pop_50$Key%in%"Connectors", ]),
     Depth100 =rownames(nodes_Pop_100[nodes_Pop_100$Key%in%"Connectors", ])) %>%
  ggVennDiagram(label_alpha = 0)


ExtarctHub <-function(nodes1, nodes2, nodes3, nodes4, Var){ 
  rownames(nodes1) <- 1:nrow(nodes1)
  nodes1$Crop <- rep(Var, times=nrow(nodes1))
  #head(nodes1) %T>% print()
  rownames(nodes2) <- 1:nrow(nodes2)
  nodes2$Crop <- rep(Var, times=nrow(nodes2))
  rownames(nodes3) <- 1:nrow(nodes3)
  nodes3$Crop <- rep(Var, times=nrow(nodes3))
  rownames(nodes4) <- 1:nrow(nodes4)
  nodes4$Crop <- rep(Var, times=nrow(nodes4))
  rbind(nodes1[nodes1$Key%in%"Module hubs", ],
        nodes2[nodes2$Key%in%"Module hubs", ],
        nodes3[nodes3$Key%in%"Module hubs", ],
        nodes4[nodes4$Key%in%"Module hubs", ],
        nodes1[nodes1$Key%in%"Network hubs", ],
        nodes2[nodes2$Key%in%"Network hubs", ],
        nodes3[nodes3$Key%in%"Network hubs", ],
        nodes4[nodes4$Key%in%"Network hubs", ]) -> key
  key <- key[, c(1,4:7,9:11,13:15, 18:20)]
  key$RootAbund <- round(key$RootAbund, 2)
  key$Betweenness <- round(key$Betweenness, 2)
  return(key)
}

# Extract sequences for BLAST --------------------------------------------------------------
filterTaxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

filterTaxa(physeq_fb_Pop, unique(hubs_Pop$OTU_ID))

write.dna(refseq(filterTaxa(physeq_fb_Pop, unique(hubs_Pop$OTU_ID))),
          format="fasta", 
          colsep="", 
          file="poplar_hubs.fasta")

write.dna(refseq(filterTaxa(physeq_fb_Pra, unique(hubs_Pra$OTU_ID))),
          format="fasta", 
          colsep="", 
          file="prairie_hubs.fasta")

write.dna(refseq(filterTaxa(physeq_fb_Sw, unique(hubs_Sw$OTU_ID))),
          format="fasta", 
          colsep="", 
          file="switchgrass_hubs.fasta")

# Correcting Taxonomies --------------------------------------------------------------------
str(hubs_Pop)
hubs_Pop <-
  ExtarctHub(nodes_Pop_10,nodes_Pop_25,nodes_Pop_50,nodes_Pop_100, "Poplar")
hubs_Pop$Taxonomy <-
  as.character(hubs_Pop$Taxonomy)
hubs_Pop[hubs_Pop == "FOTU_550-Rhizophydiales"] <- "FOTU_550-Operculomyces sp."
hubs_Pop[hubs_Pop == "POTU_37-Uncultured 76 sp."] <- "POTU_37-Gemmatimonadetes"
hubs_Pop[hubs_Pop == "POTU_26-Comamonadaceae"] <- "POTU_26-Rhizobacter sp."
hubs_Pop[hubs_Pop == "POTU_6310-Uncultured 48"] <- "POTU_6310-Actinobacteria"
hubs_Pop[hubs_Pop == "POTU_1451-Uncultured 78 sp."] <- "POTU_1451-Chloroflexi"

hubs_Pop %>%
  select(Crop,Taxonomy, RootAbund,Niche_ind, Poplar, fit_class, Key, Degree, Betweenness)

hubs_Pra <-
  ExtarctHub(nodes_Pra_10,nodes_Pra_25,nodes_Pra_50,nodes_Pra_100, "Prairie")
hubs_Pra$Taxonomy <-
  as.character(hubs_Pra$Taxonomy)
hubs_Pra[hubs_Pra == "POTU_5-Kd4-96"] <- "POTU_5-Chloroflexi"
hubs_Pra[hubs_Pra == "POTU_26-Comamonadaceae"] <- "POTU_26-Rhizobacter sp."
hubs_Pra[hubs_Pra == "POTU_35-Micromonosporaceae"] <- "POTU_35-Actinoplanes sp."
hubs_Pra[hubs_Pra == "POTU_131-Uncultured 48"] <- "POTU_131-Gaiella sp."
hubs_Pra[hubs_Pra == "POTU_752-Ad3"] <- "POTU_752-Chloroflexi"

hubs_Pra %>%
  select(Crop, Taxonomy, RootAbund, Niche_ind, Prairie, fit_class, Key, Degree, Betweenness)

hubs_Sw <-
  ExtarctHub(nodes_Sw_10,nodes_Sw_25,nodes_Sw_50,nodes_Sw_100, "Switchgrass")
hubs_Sw$Taxonomy <-
  as.character(hubs_Sw$Taxonomy)
hubs_Sw[hubs_Sw == "FOTU_25-Gs20"] <- "FOTU_25-Mucoromycotina"
hubs_Sw[hubs_Sw == "FOTU_51-Fungi"] <- "FOTU_51-Cladosporium sp."
hubs_Sw[hubs_Sw == "POTU_18-Rb41 sp."] <- "POTU_18-Acidobacteriota"

hubs_Sw %>%
  select(Crop, Taxonomy, RootAbund, Niche_ind, Switchgrass, fit_class, Key, Degree, Betweenness)

# Plot Hubs table --------------------------------------------------------------------------
hubs_all_crops <-
  rbind(hubs_Pop %>% 
          select(Taxonomy, Key, Crop, Poplar, RootAbund, Niche_ind,  fit_class, Degree, Betweenness) %>%
          rename(Depth = Poplar),
        rbind(hubs_Pra %>%
                select(Taxonomy, Key, Crop, Prairie, RootAbund, Niche_ind, fit_class, Degree, Betweenness)%>%
                rename(Depth = Prairie),
              hubs_Sw %>%
                select(Taxonomy, Key, Crop, Switchgrass, RootAbund, Niche_ind,  fit_class, Degree, Betweenness) %>%
                rename(Depth = Switchgrass)))

hubs_all_crops$Key <-
  recode(hubs_all_crops$Key,
         "Module hubs" = "Module hub",
         "Network hubs" = "Network hub")
hubs_all_crops$fit_class <-
  recode(
    hubs_all_crops$fit_class,
    "Above prediction" = "Above",
    "As predicted" = "Neutral",
    "Below prediction" = "Below")

# in hubs table Depth and Crop represent wwhich network these hubs were found.
# Indicators only for niche
hubs_all_crops

grid.draw(tableGrob(hubs_all_crops,
                    theme = ttheme_minimal(
                      base_size = 7,
                      base_colour = "black",
                      base_family = "",
                      padding = unit(c(2, 2), "mm")
                    )))


write.csv(hubs_all_crops, "hubs_all_crops.csv")

# ********************************************----------------------------------------------
# >>> Neutral models and PCA plot ----------------------------------------------------------
library(ggfortify)
library(cluster)

# ALL OTUs neutral and non-neutral taxa ----------------------------------------------------
MakeNeutralDf <- function(nfit1, nfit2, nfit3, nfit4, Var){
  if (class(nfit1)=="list"){
    nfit1 <- nfit1[[2]]
    nfit2 <- nfit2[[2]]
    nfit3 <- nfit3[[2]]
    nfit4 <- nfit4[[2]]
  }else{}
  # generate the dataframe
  df <-
    cbind(nfit1 %>%
            group_by(fit_class) %>%
            summarise(Count = table(fit_class)/nrow(nfit1)*100),
          nfit2 %>%
            group_by(fit_class) %>%
            summarise(Count = table(fit_class)/nrow(nfit2)*100),
          nfit3 %>%
            group_by(fit_class) %>%
            summarise(Count = table(fit_class)/nrow(nfit3)*100),
          nfit4 %>%
            group_by(fit_class) %>%
            summarise(Count = table(fit_class)/nrow(nfit4)*100))
  df <- df[, c(2, 4, 6, 8)]
  
  colnames(df) <-
    c("0-10cm", "10-25cm", 
      "25-50cm", "50-100cm")
  df_new <- as.data.frame(t(as.matrix(df)))
  colnames(df_new) <- c("Above", "Neutral", "Below")
  df_new$Crop <- rep(Var, 4)
  df_new$Depth <- rownames(df_new)
  rownames(df_new) <- 1:nrow(df_new)
  return(df_new)
}

MakeNeutralDf(nfit_fungi_Pop_10, nfit_fungi_Pop_25, nfit_fungi_Pop_50, nfit_fungi_Pop_100, "Poplar")


df_neutral_all_fun <-
  data.frame(
    bind_rows(
      MakeNeutralDf(nfit_fungi_Pop_10, nfit_fungi_Pop_25, nfit_fungi_Pop_50, nfit_fungi_Pop_100, "Poplar"),
      MakeNeutralDf(nfit_fungi_Pra_10, nfit_fungi_Pra_25, nfit_fungi_Pra_50, nfit_fungi_Pra_100, "Prairie"),
      MakeNeutralDf(nfit_fungi_Sw_10, nfit_fungi_Sw_25, nfit_fungi_Sw_50, nfit_fungi_Sw_100, "Switchgrass")),
    Kingdom = rep("Fungi", 12))

# no significant differences across groups -------------------------------------------------
CompSampl <- function(df, level, formula){
  require(multcompView)
  df %>% 
    melt() %>%
    group_by(Depth) %>% 
    as.data.frame() -> df_new
  compare_means(formula, data = df_new[df_new$variable==level,],
                p.adjust.method = "none") -> test_CC 
  return(test_CC)
}

CompSampl(df_neutral_all_fun, "Above", formula(value ~ Depth)) 
CompSampl(df_neutral_all_fun, "Below", formula(value ~ Depth)) 
CompSampl(df_neutral_all_fun, "Neutral", formula(value ~ Depth)) 


df_neutral_all_fun %>%
  melt() %>%
    ggplot(aes(x=Depth, y=log(value), color = variable)) +
    geom_boxplot() +  
    facet_wrap(~variable ) +
    theme_classic() +
    scale_fill_manual("Prediction", values = paletteCB3) +
    scale_color_manual("Prediction", values = paletteCB3) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), 
          axis.text.x = element_text(size = 7, angle = 45, vjust = 1.1, hjust = 1),
          axis.title = element_text(angle = 0, size = 8, face = "bold"),
          axis.text.y = element_text(size = 7), 
          legend.margin = margin(0,0,0,0, unit='cm'),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.25, "cm"), 
          legend.title = element_text(size = 8, face = "bold", hjust = 0, vjust = 0.5), 
          legend.text = element_text(size = 8), legend.position = "bottom") +
    labs(title = "Fungi\nNeutral and non-Neutral core OTUs",
         x=NULL, y="Proportion") 



# Summing Fungi and Prokaryotes together ---------------------------------------------------
df_neutral_all_fun_all <-
  df_neutral_all_fun %>%
  select(Above, Neutral, Below, Depth) %>%
  group_by(Depth) %>% 
  dplyr::summarise(Above = sum(Above),
                   Neutral = sum(Neutral),
                   Below = sum(Below))

df_neutral_all_fun_all <-
  apply(df_neutral_all_fun_all[-1], 1, function(x){x/sum(x)} ) %>% t() %>% 
  melt()

df_neutral_all_fun_all$Depth <- df_neutral_all_fun$Depth

fig_S5_neut_fungi <-
  df_neutral_all_fun_all %>%
  ggplot(aes(x = Depth, y = value, color = Var2, fill=Var2)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  facet_wrap(~Var2) +
  # scale_x_discrete(limits = rev(levels(df$Depth))) +
  #coord_flip() +
  ylim(0, 1) +
  scale_fill_manual("Prediction", values = paletteCB3) +
  scale_color_manual("Prediction", values = paletteCB3) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), 
        axis.text.x = element_text(size = 7, angle = 45, vjust = 1.1, hjust = 1),
        axis.title = element_text(angle = 0, size = 8, face = "bold"),
        axis.text.y = element_text(size = 7), 
        legend.margin = margin(0,0,0,0, unit='cm'),
        legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.25, "cm"), 
        legend.title = element_text(size = 8, face = "bold", hjust = 0, vjust = 0.5), 
        legend.text = element_text(size = 8), legend.position = "bottom") +
  labs(title = "Fungi\nNeutral and non-Neutral core OTUs",
       x=NULL, y="Proportion") 

fig_S5_neut_fungi


df_neutral_all_prok <-
  data.frame(
    bind_rows(
      MakeNeutralDf(nfit_prok_Pop_10, nfit_prok_Pop_25, nfit_prok_Pop_50, nfit_prok_Pop_100, "Poplar"),
      MakeNeutralDf(nfit_prok_Pra_10, nfit_prok_Pra_25, nfit_prok_Pra_50, nfit_prok_Pra_100, "Prairie"),
      MakeNeutralDf(nfit_prok_Sw_10, nfit_prok_Sw_25, nfit_prok_Sw_50, nfit_prok_Sw_100, "Switchgrass")),
    Kingdom = rep("Prokaryotes", 12))

df_neutral_all_prok


df_neutral_all_prok %>%
  melt() %>%
  ggplot(aes(x=Depth, y=log(value), color = variable)) +
  geom_boxplot() +  
  facet_wrap(~variable ) +
  theme_classic() +
  scale_fill_manual("Prediction", values = paletteCB3) +
  scale_color_manual("Prediction", values = paletteCB3) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), 
        axis.text.x = element_text(size = 7, angle = 45, vjust = 1.1, hjust = 1),
        axis.title = element_text(angle = 0, size = 8, face = "bold"),
        axis.text.y = element_text(size = 7), 
        legend.margin = margin(0,0,0,0, unit='cm'),
        legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.25, "cm"), 
        legend.title = element_text(size = 8, face = "bold", hjust = 0, vjust = 0.5), 
        legend.text = element_text(size = 8), legend.position = "bottom") +
  labs(title = "Fungi\nNeutral and non-Neutral core OTUs",
       x=NULL, y="Proportion") 

df_neutral_all_prok_all <-
  df_neutral_all_prok %>%
  select(Above, Neutral, Below, Depth) %>%
  group_by(Depth) %>% 
  dplyr::summarise(Above = sum(Above),
                   Neutral = sum(Neutral),
                   Below = sum(Below))

df_neutral_all_prok_all <-
  apply(df_neutral_all_prok_all[-1], 1, function(x){x/sum(x)*100} ) %>% t() %>% 
  melt()

df_neutral_all_prok_all$Depth <- df_neutral_all_prok$Depth

fig_S5_neut_prok <-
  df_neutral_all_prok_all %>%
  ggplot(aes(x = Depth, y = value, color = Var2, fill=Var2)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  facet_wrap(~Var2) +
  # scale_x_discrete(limits = rev(levels(df$Depth))) +
  #coord_flip() +
  ylim(0, 1) +
  scale_fill_manual("Prediction", values = paletteCB3) +
  scale_color_manual("Prediction", values = paletteCB3) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), 
        axis.text.x = element_text(size = 7, angle = 45, vjust = 1.1, hjust = 1),
        axis.title = element_text(angle = 0, size = 8, face = "bold"),
        axis.text.y = element_text(size = 7), 
        legend.margin = margin(0,0,0,0, unit='cm'),
        legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.25, "cm"), 
        legend.title = element_text(size = 8, face = "bold", hjust = 0, vjust = 0.5), 
        legend.text = element_text(size = 8), legend.position = "bottom") +
  labs(title = "Prokaryotes\nNeutral and non-Neutral core OTUs",
       x=NULL, y="Proportion") 

fig_S5_neut_prok

ggarrange(fig_S5_neut_fungi,
          fig_S5_neut_prok,
          labels = c("A", "B"),
          ncol = 2,
          nrow = 1,
          common.legend = TRUE,
          legend = "bottom")

# PERMNOVA on ALL --------------------------------------------------------------------------
# Permanova and betadisper - Fungi
# Depth
adonis(scale(df_neutral_all_fun[,1:3]) ~ Depth, #strata =  df_neutral_all_fun[4:5]$Crop,
       df_neutral_all_fun[4:5], method = "euclidean",
       permutations=9999)

pairwise.perm.manova(vegdist(df_neutral_all_fun[,1:3], method="euclidean"), 
                     df_neutral_all_fun[4:5]$Depth, p.method = "none", nperm=9999)

pairwise.adonis(vegdist(df_neutral_all_fun[,1:3], method="euclidean"),
                df_neutral_all_fun[4:5]$Depth, 
                p.adjust.m="BH")

anova(betadisper(dist(scale(df_neutral_all_fun[,1:3]), method="euclidean"), 
                 df_neutral_all_fun[4:5]$Depth), permutations = 9999)

#Crop
adonis(scale(df_neutral_all_fun[,1:3]) ~ Crop, #strata =  df_neutral_all_fun[5:7]$Depth,
       df_neutral_all_fun[4:5], method = "euclidean",
       permutations=9999)

anova(betadisper(dist(scale(df_neutral_all_fun[,1:3]), method="euclidean"), 
                 df_neutral_all_fun[4:5]$Crop), permutations = 9999)

# Permanova and betadisper - Prokaryotes
# Depth
adonis(scale(df_neutral_all_prok[,1:3]) ~ Depth, #strata =  df_neutral_all_prok[4:5]$Crop,  
       df_neutral_all_prok[4:5], method = "euclidean",
       permutations=9999)

anova(betadisper(dist(scale(df_neutral_all_prok[,1:3]), method="euclidean"), 
                 df_neutral_all_prok[4:5]$Depth), permutations = 9999)

# Crop
adonis(scale(df_neutral_all_prok[,1:3]) ~ Crop, #strata =  df_neutral_all_prok[4:5]$Depth, 
       df_neutral_all_prok[4:5], method = "euclidean",
       permutations=9999)

anova(betadisper(dist(scale(df_neutral_all_prok[,1:3]), method="euclidean"), 
                 df_neutral_all_prok[4:5]$Crop), permutations = 9999)


# CORE neutral and non-neutral taxa ---------------------------------------------------------
net_stats_neutral_fungi <-
  rbind(
    StatsNeutralCoreFun(nfit_fungi_Pop_10, "0-10cm","Poplar"),
    StatsNeutralCoreFun(nfit_fungi_Pop_25, "10-25cm","Poplar"),
    StatsNeutralCoreFun(nfit_fungi_Pop_50, "25-50cm","Poplar"),
    StatsNeutralCoreFun(nfit_fungi_Pop_100, "50-100cm","Poplar"),
    StatsNeutralCoreFun(nfit_fungi_Pra_10, "0-10cm","Prairie"),
    StatsNeutralCoreFun(nfit_fungi_Pra_25, "10-25cm","Prairie"),
    StatsNeutralCoreFun(nfit_fungi_Pra_50, "25-50cm","Prairie"),
    StatsNeutralCoreFun(nfit_fungi_Pra_100, "50-100cm","Prairie"),
    StatsNeutralCoreFun(nfit_fungi_Sw_10, "0-10cm","Switchgrass"),
    StatsNeutralCoreFun(nfit_fungi_Sw_25, "10-25cm","Switchgrass"),
    StatsNeutralCoreFun(nfit_fungi_Sw_50, "25-50cm","Switchgrass"),
    StatsNeutralCoreFun(nfit_fungi_Sw_100, "50-100cm","Switchgrass"))

rownames(net_stats_neutral_fungi) <- rep(1:nrow(net_stats_neutral_fungi))
net_stats_neutral_fungi <-
  spread(net_stats_neutral_fungi, Var1, value)
colnames(net_stats_neutral_fungi) <-
  c("Kingdom", "Depth", "Crop", "Above", "Neutral", "Below")
net_stats_neutral_fungi

net_stats_neutral_prok <-
  rbind(
    StatsNeutralCoreProk(nfit_prok_Pop_10, "0-10cm","Poplar"),
    StatsNeutralCoreProk(nfit_prok_Pop_25, "10-25cm","Poplar"),
    StatsNeutralCoreProk(nfit_prok_Pop_50, "25-50cm","Poplar"),
    StatsNeutralCoreProk(nfit_prok_Pop_100, "50-100cm","Poplar"),
    StatsNeutralCoreProk(nfit_prok_Pra_10, "0-10cm","Prairie"),
    StatsNeutralCoreProk(nfit_prok_Pra_25, "10-25cm","Prairie"),
    StatsNeutralCoreProk(nfit_prok_Pra_50, "25-50cm","Prairie"),
    StatsNeutralCoreProk(nfit_prok_Pra_100, "50-100cm","Prairie"),
    StatsNeutralCoreProk(nfit_prok_Sw_10, "0-10cm","Switchgrass"),
    StatsNeutralCoreProk(nfit_prok_Sw_25, "10-25cm","Switchgrass"),
    StatsNeutralCoreProk(nfit_prok_Sw_50, "25-50cm","Switchgrass"),
    StatsNeutralCoreProk(nfit_prok_Sw_100, "50-100cm","Switchgrass"))

rownames(net_stats_neutral_prok) <- 
  rep(1:nrow(net_stats_neutral_prok))
net_stats_neutral_prok <-
  rbind(net_stats_neutral_prok,
        c("Below prediction",  "0", "Prokaryote","0-10cm", "Poplar"),
        c("Below prediction",  "0", "Prokaryote","10-25cm", "Poplar"),
        c("Below prediction",  "0", "Prokaryote","0-10cm", "Prairie"),
        c("Below prediction",  "0", "Prokaryote","10-25cm", "Prairie"),
        c("Below prediction",  "0", "Prokaryote","50-100cm", "Prairie"),
        c("Below prediction",  "0", "Prokaryote","0-10cm", "Switchgrass"))

net_stats_neutral_prok <-
  spread(net_stats_neutral_prok, Var1, value)
colnames(net_stats_neutral_prok) <-
  c("Kingdom", "Depth", "Crop", "Above", "Neutral", "Below")
net_stats_neutral_prok <-
  as.data.frame(net_stats_neutral_prok)

# PERMNOVA on CORE -------------------------------------------------------------------------
# Permanova and betadisper - Fungi
df_net_stats_neutral_fun <-
  as.data.frame(scale(net_stats_neutral_fungi[,4:6]))
metadata_net_stats_neutral_fun <- 
  net_stats_neutral_fungi[,1:3]

adonis(df_net_stats_neutral_fun ~ Depth, #strata = metadata_net_stats_neutral_fun$Crop,
       metadata_net_stats_neutral_fun, method = "euclidean",
       permutations=9999)

anova(betadisper(dist(df_net_stats_neutral_fun, method="euclidean"), 
                 metadata_net_stats_neutral_fun$Depth), permutations = 9999)


adonis(df_net_stats_neutral_fun ~ Crop, #strata = metadata_net_stats_neutral_fun$Depth,
       metadata_net_stats_neutral_fun, method = "euclidean",
       permutations=9999)

anova(betadisper(dist(df_net_stats_neutral_fun, method="euclidean"), 
                 metadata_net_stats_neutral_fun$Crop), permutations = 9999)

# Permanova and betadisper - Prok
net_stats_neutral_prok$Above <- as.numeric(net_stats_neutral_prok$Above)
net_stats_neutral_prok$Neutral <- as.numeric(net_stats_neutral_prok$Neutral)
net_stats_neutral_prok$Below <- as.numeric(net_stats_neutral_prok$Below)

df_net_stats_neutral_prok <-
  as.data.frame(scale(net_stats_neutral_prok[,4:6]))
metadata_net_stats_neutral_prok <- 
  net_stats_neutral_prok[,1:3]

adonis(df_net_stats_neutral_prok ~ Depth, #strata = metadata_net_stats_neutral_prok$Crop,
       metadata_net_stats_neutral_prok, method = "euclidean",
       permutations=9999)

anova(betadisper(dist(df_net_stats_neutral_prok, method="euclidean"), 
                 metadata_net_stats_neutral_prok$Depth), permutations = 9999)

adonis(df_net_stats_neutral_prok ~ Crop, #strata = metadata_net_stats_neutral_prok$Depth,
       metadata_net_stats_neutral_prok, method = "euclidean",
       permutations=9999)

anova(betadisper(dist(df_net_stats_neutral_prok, method="euclidean"), 
                 metadata_net_stats_neutral_prok$Crop), permutations = 9999)


# Plot PCA ---------------------------------------------------------------------------------
PlotNeutralPCA <- function(df, Var){
  pca_plot <-
    autoplot(
      prcomp(x = df, scale= TRUE, center=TRUE), # calculates principal compenents and pltos with ggplot2
      data = Var, label = FALSE, shape = "Crop", # add metadate, labels of objects
      loadings = TRUE, loadings.colour = "black", size=2,
      frame = TRUE, frame.colour = "Depth", loadings.label.colour = "black",
      loadings.label = TRUE, loadings.label.size = 4, loadings.label.repel = TRUE) +
    labs(title = "PCA") +
    scale_colour_manual(values = paletteCB4) +
    scale_fill_manual(values = paletteCB4) +
    scale_shape_manual(values = c(21,22,24)) +
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
    theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
          legend.spacing.y = unit(0.1, "cm"), legend.box.margin = ) +
    theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
    theme(legend.position="right") +
    guides(shape = guide_legend(keywidth=0.6, keyheight=0.6))
  # guides(color = guide_legend(ncol = 2), #title.position="top"
  #         fill = guide_legend(ncol= 2),
  #         shape = guide_legend(ncol = 1)) +
  # theme(legend.margin = margin(0,-0.5,0,0, unit="cm")) # reduce space betweem legends
  return(pca_plot)
}

PlotNeutralPCA(net_stats_neutral_fungi[4:6], net_stats_neutral_fungi[1:3])
PlotNeutralPCA(net_stats_neutral_prok[4:6], net_stats_neutral_prok[1:3])

# **** FIGURE 5E-G - PCA on All and Core proportions -------------------------------------- 
ggarrange(PlotNeutralPCA(net_stats_neutral_fungi[4:6], net_stats_neutral_fungi[1:3]) +
            labs(title = "PCA - Fungi"),
          PlotNeutralPCA(net_stats_neutral_prok[4:6], net_stats_neutral_prok[1:3]) +
            labs(title = "PCA - Prokaryotes"),
          ncol = 2,
          nrow = 1, 
          common.legend = TRUE,
          legend = "bottom")

ggarrange(PlotNeutralPCA(net_stats_neutral_fungi[4:6], net_stats_neutral_fungi[1:3]) +
            labs(title = "PCA - Fungi"),
          PlotNeutralPCA(net_stats_neutral_prok[4:6], net_stats_neutral_prok[1:3]) +
            labs(title = "PCA - Prokaryotes"),
          ncol = 2,
          nrow = 1, 
          common.legend = TRUE,
          legend = "bottom")

# **** FIGURE 5E - PCA on All and Core proportions ---------------------------------------- 
PlotNeutralPCA(df_neutral_all_fun[,1:3], df_neutral_all_fun[,4:5]) +
  annotate("text", Inf, -Inf, 
           label = expression(paste("adonis ", italic(R) ^ 2,"= 53.1% ", italic(p), "= 0.035"),
                              parse=TRUE), size = 2.2, hjust = 1, vjust = -1.5) +
  annotate("text", Inf, -Inf, 
           label = expression(paste("betadisper ", italic(p), "= 0.896"),
                              parse=TRUE), size = 2.2, hjust = 1, vjust = -0.2) 

# **** FIGURE 5G - PCA on All and Core proportions ---------------------------------------- 
PlotNeutralPCA(df_neutral_all_prok[,1:3], df_neutral_all_prok[,4:5])


# PLOT BARPLOT NEUTRAL MODEL ---------------------------------------------------------------
BarPlotNeutral <- function(df){
  Bar_Neutral <-
    ggplot(df, aes(x = Depth, y = value, color = variable, fill=variable)) + 
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    facet_wrap(~Crop) +
    # scale_x_discrete(limits = rev(levels(df$Depth))) +
    #coord_flip() +
    scale_fill_manual("Prediction", values = paletteCB3) +
    scale_color_manual("Prediction", values = paletteCB3) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), 
          axis.text.x = element_text(size = 7, angle = 45, vjust = 1.1, hjust = 1),
          axis.title = element_text(angle = 0, size = 8, face = "bold"),
          axis.text.y = element_text(size = 7), 
          legend.margin = margin(0,0,0,0, unit='cm'),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.25, "cm"), 
          legend.title = element_text(size = 8, face = "bold", hjust = 0, vjust = 0.5), 
          legend.text = element_text(size = 8), legend.position = "bottom") +
    labs(title = "Neutral and non-Neutral core OTUs",
         x=NULL, y="Proportion") 
  # guides(color=guide_legend(ncol=1,title.position="top",),
  #        fill=guide_legend(ncol=1,title.position="top",))
  return(Bar_Neutral)
}

# **** FIGURE 5C and D -neutral and non-neutral core OTUs --------------------------------- 
BarPlotNeutral(melt(net_stats_neutral_fungi))
BarPlotNeutral(melt(net_stats_neutral_prok))

# **** FIGURE 5S - neutral and non-neutral total OTUs ------------------------------------- 
BarPlotNeutral(melt(df_neutral_all_fun))
BarPlotNeutral(melt(df_neutral_all_prok))

Fig_5S_bars <-
  ggarrange(BarPlotNeutral(melt(df_neutral_all_fun)) +
            labs(title = "Neutral and non-Neutral Fungi"),
          BarPlotNeutral(melt(df_neutral_all_prok))+
            labs(title = "Neutral and non-Neutral Prokaryotes"),
          ncol = 2,
          nrow = 1, 
          common.legend = TRUE,
          legend = "bottom")

Fig_5S_bars

# Extracting R2 and m values ---------------------------------------------------------------
df_fungi_Pop_m_r2 <-
  rbind(
    nfit_fungi_Pop_10[[1]][c(1,7, 23,24)],
    nfit_fungi_Pop_25[[1]][c(1,7, 23,24)],
    nfit_fungi_Pop_50[[1]][c(1,7, 23,24)],
    nfit_fungi_Pop_100[[1]][c(1,7, 23,24)])
df_fungi_Pop_m_r2$Depth <- as.factor(c("0-10cm", "10-25cm", "25-50cm", "50-100cm"))
df_fungi_Pop_m_r2$Crop <- as.factor(rep("Poplar",4))
df_fungi_Pop_m_r2$Kingdom <- as.factor(rep("Fungi",4))
df_fungi_Pop_m_r2

df_fungi_Pra_m_r2 <-
  rbind(
    nfit_fungi_Pra_10[[1]][c(1,7, 23,24)],
    nfit_fungi_Pra_25[[1]][c(1,7, 23,24)],
    nfit_fungi_Pra_50[[1]][c(1,7, 23,24)],
    nfit_fungi_Pra_100[[1]][c(1,7, 23,24)])
df_fungi_Pra_m_r2$Depth <- as.factor(c("0-10cm", "10-25cm", "25-50cm", "50-100cm"))
df_fungi_Pra_m_r2$Crop <- as.factor(rep("Prairie",4))
df_fungi_Pra_m_r2$Kingdom <- as.factor(rep("Fungi",4))
df_fungi_Pra_m_r2

df_fungi_Sw_m_r2 <-
  rbind(
    nfit_fungi_Sw_10[[1]][c(1,7, 23,24)],
    nfit_fungi_Sw_25[[1]][c(1,7, 23,24)],
    nfit_fungi_Sw_50[[1]][c(1,7, 23,24)],
    nfit_fungi_Sw_100[[1]][c(1,7, 23,24)])
df_fungi_Sw_m_r2$Depth <- as.factor(c("0-10cm", "10-25cm", "25-50cm", "50-100cm"))
df_fungi_Sw_m_r2$Crop <- as.factor(rep("Switchgrass",4))
df_fungi_Sw_m_r2$Kingdom <- as.factor(rep("Fungi",4))
df_fungi_Sw_m_r2


df_prok_Pop_m_r2 <-
  rbind(
    nfit_prok_Pop_10[[1]][c(1,7, 23,24)],
    nfit_prok_Pop_25[[1]][c(1,7, 23,24)],
    nfit_prok_Pop_50[[1]][c(1,7, 23,24)],
    nfit_prok_Pop_100[[1]][c(1,7, 23,24)])
df_prok_Pop_m_r2$Depth <- as.factor(c("0-10cm", "10-25cm", "25-50cm", "50-100cm"))
df_prok_Pop_m_r2$Crop <- as.factor(rep("Poplar",4))
df_prok_Pop_m_r2$Kingdom <- as.factor(rep("Prokaryotes",4))
df_prok_Pop_m_r2

df_prok_Pra_m_r2 <-
  rbind(
    nfit_prok_Pra_10[[1]][c(1,7, 23,24)],
    nfit_prok_Pra_25[[1]][c(1,7, 23,24)],
    nfit_prok_Pra_50[[1]][c(1,7, 23,24)],
    nfit_prok_Pra_100[[1]][c(1,7, 23,24)])
df_prok_Pra_m_r2$Depth <- as.factor(c("0-10cm", "10-25cm", "25-50cm", "50-100cm"))
df_prok_Pra_m_r2$Crop <- as.factor(rep("Prairie",4))
df_prok_Pra_m_r2$Kingdom <- as.factor(rep("Prokaryotes",4))
df_prok_Pra_m_r2

df_prok_Sw_m_r2 <-
  rbind(
    nfit_prok_Sw_10[[1]][c(1,7, 23,24)],
    nfit_prok_Sw_25[[1]][c(1,7, 23,24)],
    nfit_prok_Sw_50[[1]][c(1,7, 23,24)],
    nfit_prok_Sw_100[[1]][c(1,7, 23,24)])
df_prok_Sw_m_r2$Depth <- as.factor(c("0-10cm", "10-25cm", "25-50cm", "50-100cm"))
df_prok_Sw_m_r2$Crop <- as.factor(rep("Switchgrass",4))
df_prok_Sw_m_r2$Kingdom <- as.factor(rep("Prokaryotes",4))
df_prok_Sw_m_r2

df_fungi_m_r2 <-
  rbind(df_fungi_Pop_m_r2,
        df_fungi_Pra_m_r2,
        df_fungi_Sw_m_r2)

df_prok_m_r2 <-
  rbind(df_prok_Pop_m_r2,
        df_prok_Pra_m_r2,
        df_prok_Sw_m_r2)

# Plot PCA m and R2 ------------------------------------------------------------------------
# Depth
adonis(scale(df_fungi_m_r2[,1:2]) ~ Depth, #strata =  df_fungi_m_r2[5:7]$Crop,
       df_fungi_m_r2[5:7], method = "euclidean",
       permutations=9999)

anova(betadisper(dist(scale(df_fungi_m_r2[,1:2]), method="euclidean"), 
                 df_fungi_m_r2[5:7]$Depth), permutations = 9999)

#Crop
adonis(scale(df_fungi_m_r2[,1:2]) ~ Crop, #strata =  df_fungi_m_r2[5:7]$Depth,
       df_fungi_m_r2[5:7], method = "euclidean",
       permutations=9999)

anova(betadisper(dist(scale(df_fungi_m_r2[,1:2]), method="euclidean"), 
                 df_fungi_m_r2[5:7]$Crop), permutations = 9999)


# Depth
adonis(scale(df_prok_m_r2[,1:2]) ~ Depth, #strata =  df_prok_m_r2[5:7]$Crop,  
       df_prok_m_r2[5:7], method = "euclidean",
       permutations=9999)

pairwise.adonis(vegdist(scale(df_prok_m_r2[,1:2]), method="euclidean"),
                df_prok_m_r2[5:7]$Depth, 
                p.adjust.m="none")


anova(betadisper(dist(scale(df_prok_m_r2[,1:2]), method="euclidean"), 
                 df_prok_m_r2[5:7]$Depth), permutations = 9999)

# Crop
adonis(scale(df_prok_m_r2[,1:2]) ~ Crop,#strata =  df_prok_m_r2[5:7]$Depth, 
       df_prok_m_r2[5:7], method = "euclidean",
       permutations=9999)

anova(betadisper(dist(scale(df_prok_m_r2[,1:2]), method="euclidean"), 
                 df_prok_m_r2[5:7]$Crop), permutations = 9999)


# **** FIGURE 5F - R squared and m fungi ---------------------------------------------------
PlotNeutralPCA(df_fungi_m_r2[,1:2], df_fungi_m_r2[5:7]) +
  annotate("text", Inf, -Inf, 
           label = expression(paste("adonis ", italic(R) ^ 2,"= 45.1% ", italic(p), "= 0.117"),
                              parse=TRUE), size = 2.2, hjust = 1, vjust = -1.5) +
  annotate("text", Inf, -Inf, 
           label = expression(paste("betadisper ", italic(p), "= 0.742"),
                              parse=TRUE), size = 2.2, hjust = 1, vjust = -0.2) 

# **** FIGURE 5H - R squared and m fungi ---------------------------------------------------
PlotNeutralPCA(df_prok_m_r2[,1:2], df_prok_m_r2[5:7]) +
  annotate("text", Inf, -Inf, 
           label = expression(paste("adonis ", italic(R) ^ 2,"= 66.8% ", italic(p) ^ 2, "= 0.03"),
                              parse=TRUE), size = 2.2, hjust = 1.7, vjust = -1.5) +
  annotate("text", Inf, -Inf, 
           label = expression(paste("betadisper ", italic(p), "= 0.574"),
                              parse=TRUE), size = 2.2, hjust = 2.1, vjust = -0.2) 

# ********************************************----------------------------------------------
# Remember, the bar-plot and the PCA are made on the core taxa

title1=text_grob("Switchgrass Fungi",size=10, face=2)
title2=text_grob("Switchgrass Prokaryotes",size=10, face=2)

# **** FIGURE 5 - neutral models plus PCA----------------------------------------------------
Fig_4_neutral <-
  ggarrange(
    ggarrange(
      grid.arrange(ggarrange(PlotNeutral(nfit_fungi_Sw_10) + labs(title = "0-10cm"),
                             PlotNeutral(nfit_fungi_Sw_25) + labs(title = "10-25cm"),
                             PlotNeutral(nfit_fungi_Sw_50) + labs(title = "25-50cm"),
                             PlotNeutral(nfit_fungi_Sw_100) + labs(title = "50-100cm"),
                             ncol = 4, 
                             nrow = 1), top = title1),
      grid.arrange(ggarrange(PlotNeutral(nfit_prok_Sw_10) + labs(title = "0-10cm"),
                             PlotNeutral(nfit_prok_Sw_25) + labs(title = "10-25cm"),
                             PlotNeutral(nfit_prok_Sw_50) + labs(title = "25-50cm"),
                             PlotNeutral(nfit_prok_Sw_100) + labs(title = "50-100cm"),
                             ncol = 4, 
                             nrow = 1), top = title2),
      labels = c("G","H"),
      ncol = 1,
      nrow =2), # First row with bar and PCA
    ggarrange(
      ggarrange(BarPlotNeutral(melt(net_stats_neutral_fungi)) +
                  labs(title = "Neutral and non-Neutral core Fungi"),
                BarPlotNeutral(melt(net_stats_neutral_prok))+
                  labs(title = "Neutral and non-Neutral core Prokaryotes"),
                ncol = 1, 
                nrow = 2,
                align = "hv",
                labels = c("A","B"),
                common.legend = TRUE,
                legend = "bottom"),
      ggarrange(PlotNeutralPCA(df_neutral_all_fun[,1:3], df_neutral_all_fun[,4:5]) +
                  annotate("text", Inf, -Inf, 
                           label = expression(paste("adonis ", italic(R) ^ 2,"= 53.1% ", italic(p), "= 0.038"),
                                              parse=TRUE), size = 2.2, hjust = 1, vjust = -1, color="red") +
                  annotate("text", Inf, -Inf, 
                           label = expression(paste("betadisper ", italic(p), "= 0.896"),
                                              parse=TRUE), size = 2.2, hjust = 1, vjust = -0.2, color="red") +
                  theme(plot.title = element_text(size = 10)) +
                  labs(title = "PCA - Fungi"), 
                PlotNeutralPCA(df_neutral_all_prok[,1:3], df_neutral_all_prok[,4:5]) +
                  theme(plot.title = element_text(size = 10)) +
                  labs(title = "PCA - Prokaryotes"),
                PlotNeutralPCA(df_fungi_m_r2[,1:2], df_fungi_m_r2[5:7]) +
                  labs(title = NULL),
                PlotNeutralPCA(df_prok_m_r2[,1:2], df_prok_m_r2[5:7]) +
                  annotate("text", Inf, -Inf, 
                           label = expression(paste("adonis ", italic(R) ^ 2,"= 66.8% ", italic(p), "= 0.029"),
                                              parse=TRUE), size = 2.2, hjust = 1, vjust = -1, color="red") +
                  annotate("text", Inf, -Inf, 
                           label = expression(paste("betadisper ", italic(p), "= 0.574"),
                                              parse=TRUE), size = 2.2, hjust = 1, vjust = -0.2, color="red") +
                  labs(title = NULL),
                ncol = 2,
                nrow = 2,
                align = "hv",
                labels = c("C", "D", "E", "F"),
                common.legend = TRUE,
                legend = "right"),
      ncol = 2,
      nrow = 1,
      widths = c(0.8, 1.2)),
    nrow = 2, 
    heights = c(1, 1.1)) 

Fig_4_neutral


# **** FIGURE S6A - PCA Network stats ------------------------------------------------------
library(Boruta)
stats_all
dim(stats_all)
str(stats_all)

# loop to store 100 boruta runs
sign_var <- vector(mode = "character")
for(i in 1:99) {
  sel_attr <- Boruta(Depth ~., stats_all[-c(6:12,16,18,20,30)], pValue = 0.05, 
                     mcAdj = TRUE, maxRuns=100, doTrace = 3)
  if (identical(getSelectedAttributes(sel_attr, withTentative = TRUE), character(0))==TRUE){
  }else{
    sign_var[i] <- getSelectedAttributes(sel_attr, withTentative = TRUE)
  }
}

sign_var
unique(sign_var)

# no attributes found for Crops
Boruta(Crop ~., stats_all[-c(6:12,31)], pValue = 0.05, 
       mcAdj = TRUE, maxRuns=100, doTrace = 3) -> test_bor

# Filter to selected variables 
stats_all[, unique(sign_var)[-5]]

Net_stats_plot <-
  autoplot(
    prcomp(stats_all[, unique(sign_var)[-5]], scale= TRUE, center=TRUE), 
    data= stats_all[, c(30,31)], label = F, shape = "Crop", # add metadate, labels of objects
    loadings = TRUE, loadings.colour = "black", size=2,
    frame = TRUE, frame.colour = "Depth", loadings.label.colour = "black",
    loadings.label = TRUE, loadings.label.size = 4, loadings.label.repel = TRUE) +
  labs(title = "Network properties") +
  scale_colour_manual(values = paletteCB4) +
  scale_fill_manual(values = paletteCB4) +
  scale_shape_manual(values = c(21,22,24)) +
  theme_classic() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
  theme(legend.position="right") +
  annotate("text", Inf, -Inf, 
           label = expression(paste("adonis ", italic(R) ^ 2,"= 51.5% ", italic(p), "= 0.028"),
                              parse=TRUE), size = 3, hjust = 1, vjust = -1.3) +
  annotate("text", Inf, -Inf, 
           label = expression(paste("betadisper ", italic(p), "= 0.178"),
                              parse=TRUE), size = 3, hjust = 1, vjust = -0.2) 

Net_stats_plot

df_net_stats <-
  as.data.frame(scale(stats_all[, unique(sign_var)[-5]]))
metadata_net_stats <- 
  stats_all[, c(30,31)]

adonis(df_net_stats ~ Depth, 
       metadata_net_stats, method = "euclidean",
       permutations=9999)

anova(betadisper(dist(df_net_stats, method="euclidean"), 
                 metadata_net_stats$Depth), permutations = 9999)

boxplot(betadisper(dist(df_net_stats, method="euclidean"), 
                   metadata_net_stats$Depth))

# Plot Figure 6 
Fig_6_pca_heat <-
  ggarrange(
    PlotHeatKing(df_king_all) +
      geom_text(aes(label = Mean), size=2, color="white") +
      scale_fill_gradient2("Edges No.", low="blue", mid="grey", high="red", 
                           midpoint=max(df_king_all$Mean)/2,  
                           limits=c(0, max(df_king_all$Mean))) +
      labs(title = "Intra- and inter-phylum links"),
    Net_stats_plot,
    ncol = 2,
    nrow = 1,
    widths = c(1.4,1),
    labels = c("A", "B"))

Fig_6_pca_heat

# plot Pi-Zi 
PlotZiPi(`df_zipi_0-10cm_Poplar`)

# ********************************************----------------------------------------------
# PLOTTING NETWORKS ATTRIBUTES -------------------------------------------------------------
library(tibble)
library(tidyr)
library(viridis)
library(ggpubr)
library(dplyr)


# Stacked bar plot for taxonomic compositon of modules -------------------------------------
PlotBarsMod <- function(abund, network, neutral_fit_f=NULL, neutral_fit_b=NULL, Var1, Var2){
  require(scales)
  igraph::delete.vertices(network, which(degree(network)==0)) -> network_new
  print(network_new)
  #CalcNodes(physeq_poplar, whole_abund, network_new, Crop) -> nodes
  CalcNodes(abund, network_new, neutral_fit_f=NULL, neutral_fit_b=NULL, Var1, Var2)-> nodes
  head(nodes) %T>% print()
  palette <- c("red","orange","gold","yellow","papayawhip",
               "peachpuff", "burlywood","salmon4", "burlywood4",
               "brown", "gray","navy","blue","dodgerblue2",
               "cyan","wheat", "black", "olivedrab", "grey20", "lightblue")
  #pie(rep(1, length(palette2)), labels = sprintf("%d (%s)",
  #   seq_along(palette2),palette2), col = palette2) %T>% print()
  #nodes$Order <- as.factor(nodes$Order)
  # generate dataframe for labels
  as.data.frame(table(nodes$Modules)) -> df
  colnames(df) <- c("Modules", "OTUs")
  df$label <- paste(df$OTUs, "OTUs")
  
  df <-
    subset(df, OTUs>1)
  
  head(df) %T>% print()
  # plot the graph
  nodes%>% 
    as_tibble() %>%
    select(AbundDepth, Phylum, BestMatch, Modules) %>%
    reshape2::melt() %>%
    mutate(Modules = as.numeric(Modules)) %>%
    assign(paste("nodes_df", Var1, Var2, sep = ""),., envir = .GlobalEnv) %>% # saving the plot 
    group_by(Modules) %>%
    nest() %>% 
    mutate(Abund = purrr::map(data, ~.x$value *100/sum(.x$value))) %>% # this generates within modules abundances
    unnest() %>%
    as.data.frame() -> nodes_new
  
  nodes_new$Modules <- as.factor(nodes_new$Modules)
  
  head(nodes_new) %T>% print()
  str(nodes_new)
  
  nodes_new[duplicated(nodes_new$Modules), ] -> nodes_new
  nodes_new %T>% print() 
  
  #nodes_new$Modules <- factor(nodes_new$Modules, 
  #                            levels = as.character(rep(1:length(unique(nodes_new$Modules)))))
  nodes_new$Phylum <- as.factor(nodes_new$Phylum)
  nodes_new$Phylum <- 
    factor(nodes_new$Phylum, 
           levels = factor_levels)
  # using cumulative abundances it is better, then use value and not Abund
  unique(nodes_new$Phylum) %T>% print()
  
  nodes_new %>% 
    # ggplot() + 
    # geom_bar(aes(x=100, y = Abund, 
    #              fill=fct_rev(Phylum), width=100),
    #              stat="identity", position="fill") +
    # geom_text(aes(x=sqrt(value), y=Inf, label=""), size=7) +
    # theme_graph() +
    # coord_polar(theta="y") +
    # scale_fill_hue(name="Modules", l=80) +
    # facet_wrap(~Modules, ncol=4) +
    ggplot(aes(x = Modules, y = value, fill = Phylum)) +
    geom_bar(stat = "identity") +
    #facet_wrap(~Modules, nrow = 1) +
    facet_grid(~Modules, scales = "free", space='free') +
    #ylim(0, 1500) +
    geom_text(data = df, 
              aes(x = Inf, y = Inf, label = OTUs), inherit.aes = FALSE,
              size = 2.5, 
              hjust = 1.5,
              vjust = 1.2) +
    scale_fill_manual(values = palette_phylum) + 
    scale_color_manual(values = palette_phylum) +
    labs(title = "Module taxonomic composition", x="Module", y="Cumulative relative abundance") +
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.key.height = unit(0.2, "cm"),legend.key.width = unit(0.3, "cm"),
          legend.title = element_text(size = 8, face = "bold"), 
          legend.position = "right",
          legend.text = element_text(size = 8),
          axis.text.x = element_text(angle = 0, size = 9, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(angle = 0,size = 9,hjust = 0.5,vjust = 0.5),
          axis.title = element_text(angle = 0, size = 9, face = "bold")) -> plot
  return(plot)
}


nodes_Pop_10[nodes_Pop_10$OTU_ID=="FOTU_550", ]
nodes_Pop_100[nodes_Pop_100$OTU_ID=="POTU_26", ]
nodes_Pop_25[nodes_Pop_25$OTU_ID=="POTU_22", ]

PlotBarsMod(abund_Pop_10, network_Pop_10, nfit_fungi_Pop_10, nfit_prok_Pop_10, "0-10cm", "Poplar") 
PlotBarsMod(abund_Pop_25, network_Pop_25, nfit_fungi_Pop_25, nfit_prok_Pop_25, "10-25cm", "Poplar") 
PlotBarsMod(abund_Pop_50, network_Pop_50, nfit_fungi_Pop_50, nfit_prok_Pop_50, "10-25cm","Poplar")
PlotBarsMod(abund_Pop_100, network_Pop_100, nfit_fungi_Pop_100, nfit_prok_Pop_100, "50-100cm","Poplar")


# Neutral fit on both bacteria and fungi together -------------------------------------------
# Netwrok stats summary ad Depths ----------------------------------------------------------
CombineStats <- function(){
  t(net_stats_pop)[,1:29] -> df_pop #[c(1:4,6:17),]
  t(net_stats_pra)[,1:29] -> df_pra
  t(net_stats_sw)[,1:29] -> df_sw
  rownames(df_pop) <- paste(rep("Pop_", 4), rownames(df_pop), sep = "")
  rownames(df_pra) <- paste(rep("Pra_", 4), rownames(df_pra), sep = "")
  rownames(df_sw) <- paste(rep("Sw_", 4), rownames(df_sw), sep = "")
  df_all <-
    rbind(df_pop, df_pra, df_sw)
  df_all <-
    as.data.frame(df_all[, colSums(df_all != 0) > 0])
  df_all$Crop <-
    factor(c(rep("Poplar",4),rep("Prairie",4), rep("Switchgrass",4)))
  df_all$Depth <-
    factor(c(rep(c("0-10cm", "10-25cm","25-50cm","50-100cm"), 3)))
  return(df_all) 
}

stats_all <-
  as.data.frame(CombineStats())
str(stats_all)

# Extract Neutral fit values ---------------------------------------------------------------
df_neutral_all_core <-
  stats_all[, c(6:8, 30:31)]
colnames(df_neutral_all_core) <-
  c("Above", "Neutral", "Below", "Crop", "Depth")

df_neutral_all_core

df_neutral_all_core$Depth <-
  factor(df_neutral_all_core$Depth, 
         levels = c("0-10cm", "10-25cm", "25-50cm", "50-100cm"))

# Plot PCA on nutral proportions ----------------------------------------------------------------------------------------
# PCA on standardized data: mean =0 and sd =1
Fit_stats_plot <-
  PlotNeutralPCA(df_neutral_all_core[,1:3], df_neutral_all_core[,4:5])

Fit_stats_plot

df_net_stats_neutral <-
  as.data.frame(scale(df_neutral_all_core[,1:3]))
metadata_net_stats_neutral <- 
  df_neutral_all_core[,4:5]

adonis(df_net_stats_neutral ~ Depth, 
       metadata_net_stats_neutral, method = "euclidean",
       permutations=9999)

anova(betadisper(dist(df_net_stats_neutral, method="euclidean"), 
                 metadata_net_stats_neutral$Depth), permutations = 9999)

boxplot(betadisper(dist(df_net_stats_neutral, method="euclidean"), 
                   metadata_net_stats_neutral$Depth))

# crop
adonis(as.data.frame(scale(df_neutral_all_core[,1:3])) ~ df_neutral_all_core[,4], 
       df_neutral_all_core[,4:5], method = "euclidean",
       permutations=9999)

anova(betadisper(dist(scale(df_neutral_all_core[, 1:3]), method="euclidean"), 
                 df_neutral_all_core[,4]), permutations = 9999)

CombineStatsNet <- function(){
  t(net_stats_pop_core) -> df_pop 
  t(net_stats_pra_core) -> df_pra
  t(net_stats_sw_core) -> df_sw
  rownames(df_pop) <- paste(rep("Pop_", 4), rownames(df_pop), sep = "")
  rownames(df_pra) <- paste(rep("Pra_", 4), rownames(df_pra), sep = "")
  rownames(df_sw) <- paste(rep("Sw_", 4), rownames(df_sw), sep = "")
  df_all <-
    rbind(df_pop, df_pra, df_sw)
  df_all %T>% print()
  df_all <-
    as.data.frame(df_all[, colSums(df_all != 0) > 0])
  df_all$Crop <-
    factor(c(rep("Poplar",4),rep("Prairie",4), rep("Switchgrass",4)))
  df_all$Depth <-
    factor(c(rep(c("0-10cm", "10-25cm","25-50cm","50-100cm"), 3)))
  return(df_all) 
}

stats_all <-
  as.data.frame(CombineStatsNet())
str(stats_all)


# Neutral core groups
df_neutral_all_core %>%
  melt() -> df_neutral

df_neutral$Depth <-
  factor(df_neutral$Depth, 
         levels = c("0-10cm", "10-25cm", "25-50cm", "50-100cm"))

df_neutral$variable <-
  factor(df_neutral$variable, 
         levels = c("Above", "Neutral", "Below"))

# ********************************************----------------------------------------------
