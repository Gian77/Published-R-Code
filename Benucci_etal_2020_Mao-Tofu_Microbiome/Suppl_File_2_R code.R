# ************ DATA ANALYSIS ***************************************** ------------------------
# Project name: hairy tofu - Mao tofu Microbiome
# Manuscript:   
# Authors:      Benucci GMN, Wang X, ..., Bonito G
# Affiliation:  Michigan State University
# Journal:      Food Microbiology
# Date:         February 6, 2021
# ******************************************************************** ------------------------

# WORKING ENVIRONMENT SETUP -------------------------------------------------------------------
library(styler) # hilight the code and press CTRL+SHIFT+A to style the code.
options(scipen = 9999) #to use decimals
options(max.print=100000000) # to print more lines on the display
#rm(list = ls())

# loading required packages -------------------------------------------------------------------
library(phyloseq); packageVersion("phyloseq")
library(Biostrings)
library(ggplot2)
library(ape)
library(dplyr)
library(ggpubr)

# PALETTE -------------------------------------------------------------------------------------
paletteCB6 = c("#2D3142","#058ED9","#848FA2","#599861","#FF934F", "#CC2D35")
pie(rep(1, length(paletteCB6)), labels = sprintf("%d (%s)",
                   seq_along(paletteCB6),paletteCB6), col = paletteCB6)

paletteCB4 = c("#E69F00", "#D55E00", "#CBD588", "#0072B2")
pie(rep(1, length(paletteCB4)), labels = sprintf("%d (%s)",
                   seq_along(paletteCB4),paletteCB4), col = paletteCB4)


paletteCB2 = c("grey17", "grey70")
pie(rep(1, length(paletteCB2)), labels = sprintf("%d (%s)",
                  seq_along(paletteCB2),paletteCB2), col = paletteCB2)


palette_heat <-c("#BBBBBB","#BBDF27FF","#85D44AFF","#54C568FF",
                 "#2FB47CFF","#1FA188FF","#228C8DFF","#2A788EFF")

palette_heat <-c("#BBBBBB","#ff3333","#dd4455","#bb5577",
                 "#996699","#7777bb","#5588dc","#3399ff")

pie(rep(1, length(palette_heat)), labels = sprintf("%d (%s)",
     seq_along(palette_heat),palette_heat), col = palette_heat)


# IMPORTING DATASETS --------------------------------------------------------------------------

# A) ITS UPARSE -------------------------------------------------------------------------------
otus_ITS_uparse_R1 <- read.delim("dataset_its/otu_table_ITS_UPARSE_R1.txt",row.names=1, header=TRUE, sep="\t") 
head(otus_ITS_uparse_R1)

metadata_ITS_uparse_R1 <-read.delim("mapping.txt", row.names=1, header=TRUE, sep="\t")
head(metadata_ITS_uparse_R1)

taxonomy_ITS_uparse_R1 <-read.delim("dataset_its/ITS_taxonomy.txt",header=TRUE, row.names=1, sep="\t")
head(taxonomy_ITS_uparse_R1)

otus_seq_ITS_uparse_R1 <- readDNAStringSet("dataset_its/otus_ITS.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

physeq_ITS_uparse_R1 <- phyloseq(otu_table(otus_ITS_uparse_R1, taxa_are_rows = TRUE),
                                     sample_data(metadata_ITS_uparse_R1),
                                     tax_table(as.matrix(taxonomy_ITS_uparse_R1)),
                                     otus_seq_ITS_uparse_R1) 

physeq_ITS_uparse_R1
tax_table(physeq_ITS_uparse_R1)
sample_data(physeq_ITS_uparse_R1)

tax_table(physeq_ITS_uparse_R1)[tax_table(physeq_ITS_uparse_R1)==""]<- NA
#tax_table(physeq_ITS_uparse_R1)[is.na(tax_table(physeq_ITS_uparse_R1))]<-"Unclassified"

physeq_ITS_uparse_R1 -> physeq_ITS
physeq_ITS

# # Filtering Taxonomy -------------------------------------------------------------------------
# BadTaxa <- function(physeq){
#   names(apply(tax_table(physeq), 2, function(x) which(x == "Glomeromycetes"))$Class) -> x
#   names(apply(tax_table(physeq), 2, function(x) which(x == "Pezizomycetes"))$Class) -> y
#   names(apply(tax_table(physeq), 2, function(x) which(x == "Agaricomycetes"))$Class) -> z 
#   names(apply(tax_table(physeq_ITS), 2, function(x) which(is.na(x)))$Kingdom) -> w
#   lista <- c(x,y,z,w)
#   bad_taxa <- tax_table(physeq)[taxa_names(physeq)%in%c(lista),]
#   return(bad_taxa)
# }
# 
# rownames(BadTaxa(physeq_ITS))

rownames(tax_table(
  subset_taxa(
    physeq_ITS,
    Class == "Glomeromycetes" |
      Class == "Pezizomycetes" |
      Class == "Agaricomycetes" |
      is.na(Kingdom)))) -> bad_taxa_ITS

length(bad_taxa_ITS)

physeq_ITS <- remove_taxa(physeq_ITS, bad_taxa_ITS)
otu_table(physeq_ITS) <- otu_table(physeq_ITS)[which(rowSums(otu_table(physeq_ITS)) > 0),] 
physeq_ITS


# FILTERING OUT CONTAMINANTS ITS --------------------------------------------------------------
# detecting contaminants 
sample_data(physeq_ITS)$is.neg <- sample_data(physeq_ITS)$Market == "control"
contam_prev_ITS <- isContaminant(physeq_ITS, method="prevalence", neg="is.neg", threshold=0.5)
table(contam_prev_ITS$contaminant)

# Make phyloseq object of presence-absence in negative controls and true samples
physeq_ITS_pa <- transform_sample_counts(physeq_ITS, function(abund) 1*(abund>0))
physeq_ITS_pa_neg <- prune_samples(sample_data(physeq_ITS_pa)$is.neg == "TRUE", physeq_ITS_pa)
physeq_ITS_pa_pos <- prune_samples(sample_data(physeq_ITS_pa)$is.neg == "FALSE", physeq_ITS_pa)

# Make data.frame of prevalence in positive and negative samples
df_ITS_pa_pa <- data.frame(pa.pos=taxa_sums(physeq_ITS_pa_pos),pa.neg=taxa_sums(physeq_ITS_pa_neg),
                             Contaminant=contam_prev_ITS$contaminant)

# removing contaminants form the phyloseq object ----------------------------------------------
# function to remove dab taxa
remove_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

contaminants_ITS = rownames(subset(contam_prev_ITS, contaminant%in%c("TRUE")))
length(contaminants_ITS)

physeq_ITS_clean <- remove_taxa(physeq_ITS, contaminants_ITS)
otu_table(physeq_ITS_clean) <- otu_table(physeq_ITS_clean)[which(rowSums(otu_table(physeq_ITS_clean)) > 0),] 
physeq_ITS_clean

count(as.data.frame(as.matrix(sample_data(physeq_ITS_uparse_R1))), vars = Market)
table(sample_data(physeq_ITS_uparse_R1)$Market)

#tax_table(subset_taxa(physeq_ITS_uparse_R1, Phylum=="Rozellomycota"))

# # Now removing controls samples
# physeq_ITS_clean <- subset_samples(physeq_ITS_clean, is.neg%in%c("FALSE"))
# otu_table(physeq_ITS_clean) <- otu_table(physeq_ITS_clean)[which(rowSums(otu_table(physeq_ITS_clean)) > 0),] 
# physeq_ITS_clean
# 
# any(taxa_sums(physeq_ITS_clean) == 0)
# sort(taxa_sums(physeq_ITS_clean))
# any(sample_sums(physeq_ITS_clean) == 0)
# sort(sample_sums(physeq_ITS_clean))
# 
# head(sample_data(physeq_ITS_clean))

# B) 16s UPARSE R1 ----------------------------------------------------------------------------
otus_16s_uparse_R1 <- read.delim("dataset_bact/otu_table_16s_UPARSE.txt",row.names=1, header=TRUE, sep="\t") 
head(otus_16s_uparse_R1)

metadata_16s_uparse_R1 <-read.delim("mapping.txt", row.names=1, header=TRUE, sep="\t")
head(metadata_16s_uparse_R1)

taxonomy_16s_uparse_R1 <-read.delim("dataset_bact/16s_taxonomy.txt", header=TRUE, row.names=1, sep="\t")
head(taxonomy_16s_uparse_R1)

otus_seq_16s_uparse_R1 <- readDNAStringSet("dataset_bact/otus_16S.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

physeq_16s_uparse_R1 <- phyloseq(otu_table(otus_16s_uparse_R1, taxa_are_rows = TRUE),
                                 sample_data(metadata_16s_uparse_R1),
                                 tax_table(as.matrix(taxonomy_16s_uparse_R1)),
                                 otus_seq_16s_uparse_R1) 

physeq_16s_uparse_R1
tax_table(physeq_16s_uparse_R1)
sample_data(physeq_16s_uparse_R1)

tax_table(physeq_16s_uparse_R1)[tax_table(physeq_16s_uparse_R1)==""]<- NA
tax_table(physeq_16s_uparse_R1)[is.na(tax_table(physeq_16s_uparse_R1))]<-"Unclassified"

physeq_16s_uparse_R1 -> physeq_16s
physeq_16s

# Filtering Taxonomy: Cloroplast, Mitochondria and Taq contaminants in the entire dataset
apply(tax_table(physeq_16s), 2, function(x) which(x == "Chloroplast"))
apply(tax_table(physeq_16s), 2, function(x) which(x == "Mitochondria"))
apply(tax_table(physeq_16s), 2, function(x) which(x == "Oceanospirillales"))
apply(tax_table(physeq_16s), 2, function(x) which(is.na(x)))

rownames(tax_table(
  subset_taxa(
    physeq_16s,
    Order == "Chloroplast" |
      Family == "Mitochondria" |
      Order == "Oceanospirillales" |
      is.na(Kingdom)))) -> bad_taxa_16s

physeq_16s <- remove_taxa(physeq_16s, bad_taxa_16s)
otu_table(physeq_16s) <- otu_table(physeq_16s)[which(rowSums(otu_table(physeq_16s)) > 0),] 
physeq_16s

# FILTERING OUT CONTAMINANTS 16S --------------------------------------------------------------
sample_data(physeq_16s)$is.neg <- sample_data(physeq_16s)$Market == "control"
contam_prev_16s <- isContaminant(physeq_16s, method="prevalence", neg="is.neg", threshold=0.5)
table(contam_prev_16s$contaminant)

# Make phyloseq object of presence-absence in negative controls and true samples
physeq_16s_pa <- transform_sample_counts(physeq_16s, function(abund) 1*(abund>0))
physeq_16s_pa_neg <- prune_samples(sample_data(physeq_16s_pa)$is.neg == "TRUE", physeq_16s_pa)
physeq_16s_pa_pos <- prune_samples(sample_data(physeq_16s_pa)$is.neg == "FALSE", physeq_16s_pa)

# Make data.frame of prevalence in positive and negative samples
df_16s_pa_pa <- data.frame(pa.pos=taxa_sums(physeq_16s_pa_pos),pa.neg=taxa_sums(physeq_16s_pa_neg),
                          Contaminant=contam_prev_16s$contaminant)

# removing contaminants form the phyloseq object ----------------------------------------------
contaminants_16s = rownames(subset(contam_prev_16s, contaminant%in%c("TRUE")))
length(contaminants_16s)

physeq_16s_clean <- remove_taxa(physeq_16s, contaminants_16s)
otu_table(physeq_16s_clean) <- otu_table(physeq_16s_clean)[which(rowSums(otu_table(physeq_16s_clean)) > 0),] 
physeq_16s_clean

count(as.data.frame(as.matrix(sample_data(physeq_16s_uparse_R1))), vars = Market)

# C) LSU UPARSE R1 ----------------------------------------------------------------------------
taxonomy_LSU_uparse_R1 <-read.delim("dataset_lsu/LSU_taxonomy.txt", header=TRUE, row.names=1, sep="\t")
head(taxonomy_LSU_uparse_R1)

# # Step to filter first then save, fix with BLAST on GenBank and reimport.
# 
# taxonomy_LSU_uparse_R1 <-read.delim("dataset_lsu/consensus_taxonomy.txt", header=TRUE, row.names=1, sep="\t")
# head(taxonomy_LSU_uparse_R1)
# 
# ifelse(taxonomy_LSU_uparse_R1$D_Score>=0.7, paste(taxonomy_LSU_uparse_R1$Domain), NA) -> Kingdom
# ifelse(taxonomy_LSU_uparse_R1$P_Score>=0.7, paste(taxonomy_LSU_uparse_R1$Phylum), NA) -> Phylum
# ifelse(taxonomy_LSU_uparse_R1$C_Score>=0.7, paste(taxonomy_LSU_uparse_R1$Class), NA) -> Class
# ifelse(taxonomy_LSU_uparse_R1$O_Score>=0.7, paste(taxonomy_LSU_uparse_R1$Order), NA) -> Order
# ifelse(taxonomy_LSU_uparse_R1$F_Score>=0.7, paste(taxonomy_LSU_uparse_R1$Family), NA) -> Family
# ifelse(taxonomy_LSU_uparse_R1$G_Score>=0.7, paste(taxonomy_LSU_uparse_R1$Genus), NA) -> Genus
# taxonomy_LSU <- cbind(Kingdom, Phylum, Class, Order, Family, Genus)
# rownames(taxonomy_LSU) <- rownames(taxonomy_LSU_uparse_R1)
# taxonomy_LSU <- as.data.frame(taxonomy_LSU)
# head(taxonomy_LSU)
# 
# write.csv(tax_table(physeq_LSU_uparse_R1), "LSU_taxonomy.csv")
# 
# # OTU_44 Cercospora sojina
# # OTU_89 and OTU_119 Cystofilobasidium macerans

otus_LSU_uparse_R1 <- read.delim("dataset_lsu/otu_table_LSU_UPARSE_R1.txt",row.names=1, header=TRUE, sep="\t") 
head(otus_LSU_uparse_R1)

metadata_LSU_uparse_R1 <-read.delim("mapping.txt", row.names=1, header=TRUE, sep="\t")
head(metadata_LSU_uparse_R1)

otus_seq_LSU_uparse_R1 <- readDNAStringSet("dataset_lsu/otus_LSU.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

physeq_LSU_uparse_R1 <- phyloseq(otu_table(otus_LSU_uparse_R1, taxa_are_rows = TRUE),
                                 sample_data(metadata_LSU_uparse_R1),
                                 tax_table(as.matrix(taxonomy_LSU_uparse_R1)),
                                 otus_seq_LSU_uparse_R1) 

physeq_LSU_uparse_R1
tax_table(physeq_LSU_uparse_R1)
sample_data(physeq_LSU_uparse_R1)

tax_table(physeq_LSU_uparse_R1)[tax_table(physeq_LSU_uparse_R1)==""]<- NA
tax_table(physeq_LSU_uparse_R1)[is.na(tax_table(physeq_LSU_uparse_R1))]<-"Unclassified"

physeq_LSU_uparse_R1 -> physeq_LSU_all
physeq_LSU_all

# Filtering Taxonomy
rownames(tax_table(
  subset_taxa(
    physeq_LSU_all,
    Class == "Glomeromycetes" |
      Class == "Pezizomycetes" |
      Class == "Agaricomycetes" |
      is.na(Kingdom)))) -> bad_taxa_LSU

physeq_LSU_all <- remove_taxa(physeq_LSU_all, bad_taxa_LSU)
otu_table(physeq_LSU_all) <- otu_table(physeq_LSU_all)[which(rowSums(otu_table(physeq_LSU_all)) > 0),] 
physeq_LSU_all

# FILTERING OUT CONTAMINANTS LSU --------------------------------------------------------------
sample_data(physeq_LSU_all)$is.neg <- sample_data(physeq_LSU_all)$Market == "control"
contam_prev_LSU <- isContaminant(physeq_LSU_all, method="prevalence", neg="is.neg", threshold=0.5)
table(contam_prev_LSU$contaminant)
which(contam_prev_LSU$contaminant)

# Make phyloseq object of presence-absence in negative controls and true samples
physeq_LSU_pa <- transform_sample_counts(physeq_LSU_all, function(abund) 1*(abund>0))
physeq_LSU_pa_neg <- prune_samples(sample_data(physeq_LSU_pa)$is.neg == "TRUE", physeq_LSU_pa)
physeq_LSU_pa_pos <- prune_samples(sample_data(physeq_LSU_pa)$is.neg == "FALSE", physeq_LSU_pa)

# Make data.frame of prevalence in positive and negative samples
df_LSU_pa_pa <- data.frame(pa.pos=taxa_sums(physeq_LSU_pa_pos),pa.neg=taxa_sums(physeq_LSU_pa_neg),
                           Contaminant=contam_prev_LSU$contaminant)

# removing contaminants form the phyloseq object ----------------------------------------------
contaminants_LSU = rownames(subset(contam_prev_LSU, contaminant%in%c("TRUE")))

physeq_LSU_clean <- remove_taxa(physeq_LSU_all, contaminants_LSU)
otu_table(physeq_LSU_clean) <- otu_table(physeq_LSU_clean)[which(rowSums(otu_table(physeq_LSU_clean)) > 0),] 
physeq_LSU_clean

count(as.data.frame(as.matrix(sample_data(physeq_LSU_clean))), vars = Market)

# ******************************************************---------------------------------------
# CONTAMINANTS and LIBRARY SIZES --------------------------------------------------------------

# plotting contaminants vs. true otus
PlotContam <- function(dataframe){
  ggplot(data=dataframe, aes(x=pa.neg, y=pa.pos, color=Contaminant)) + 
    geom_point(size=1) +
    theme_classic() +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 7)) +
    theme(axis.text.x = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) + 
    theme(legend.position="bottom") +
    labs(x="Prevalence (Negative Controls)",
         y="Prevalence (True Samples)") -> plot
  return(plot)
}

PlotContam(df_ITS_pa_pa) + labs(title = "Contaminants ITS") -> Fig_contam_ITS
Fig_contam_ITS
PlotContam(df_16s_pa_pa) + labs(title = "Contaminants 16S") -> Fig_contam_16s
Fig_contam_16s
PlotContam(df_LSU_pa_pa) + labs(title = "Contaminants LSU") -> Fig_contam_LSU
Fig_contam_LSU

# extract dataframe
LibrarySize <- function(physeq){
df <- as.data.frame(as.matrix(sample_data(physeq))) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(physeq)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df)) # sample numbering
return(df)
}

LibrarySize(physeq_ITS_clean) -> df_ITS
LibrarySize(physeq_16s_clean) -> df_16s
LibrarySize(physeq_LSU_clean) -> df_LSU

# histogram 
PlotLib <- function(dataframe, wid_bin){
ggplot(dataframe, aes(x = LibrarySize)) + # Histogram of sample read counts
  geom_histogram(color = "indianred", fill = "indianred", binwidth = wid_bin) +
  #facet_grid(~Treatment, scales = "free_x", space="free_x") +
  labs(x="Read number", y="Sample number") + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, size = 7, hjust = 1, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) -> plot_lib
return(plot_lib)
}

PlotLib(df_ITS, 500) + ylim(0, 15) + labs(title="Distribution of ITS\nSample Libraries") -> Fig_lib_ITS
Fig_lib_ITS
PlotLib(df_16s, 500) +ylim(0, 15) + labs(title="Distribution of 16S\nSample Libraries") -> Fig_lib_16s
Fig_lib_16s
PlotLib(df_LSU, 500) + ylim(0, 15) + labs(title="Distribution of LSU\nSample Libraries") -> Fig_lib_LSU
Fig_lib_LSU


# *** Figure S1 - library check ---------------------------------------------------------------
library(ggpubr)

ggarrange(ggarrange(Fig_contam_ITS, Fig_contam_16s, Fig_contam_LSU,
                    labels = c("A","B","C"),
                    widths = c(1,1,1),
                    align = "hv" ,
                    ncol = 3, 
                    nrow = 1, 
                    common.legend = TRUE,
                    legend = c("bottom")),
          ggarrange(Fig_lib_ITS, Fig_lib_16s, Fig_lib_LSU,
                    widths = c(1,1,1),
                    align = "hv" ,
                    ncol = 3, 
                    nrow = 1),
          heights =  c(1,1),
          ncol = 1, 
          nrow = 2) -> Fig_S1

Fig_S1


# DATASETS ------------------------------------------------------------------------------------

#ITS
physeq_ITS_clean
count(as.data.frame(as.matrix(sample_data(physeq_ITS))), vars = Market)
count(as.data.frame(as.matrix(sample_data(physeq_ITS))), vars = Type)

#16S
physeq_16s_clean
count(as.data.frame(as.matrix(sample_data(physeq_16s))), vars = Market)
count(as.data.frame(as.matrix(sample_data(physeq_16s))), vars = Type)

#LSU fungi only
physeq_LSU <- subset_taxa(physeq_LSU_clean, Kingdom=="Fungi")
otu_table(physeq_LSU) <- otu_table(physeq_LSU)[which(rowSums(otu_table(physeq_LSU)) > 0),] 
physeq_LSU

count(as.data.frame(as.matrix(sample_data(physeq_LSU))), vars = Market)
count(as.data.frame(as.matrix(sample_data(physeq_LSU))), vars = Type)

# LSU animals
physeq_Animal <- subset_taxa(physeq_LSU_clean, Kingdom!="Fungi")
otu_table(physeq_Animal) <- otu_table(physeq_Animal)[which(rowSums(otu_table(physeq_Animal)) > 0),] 
physeq_Animal
tax_table(physeq_Animal)

count(as.data.frame(as.matrix(sample_data(physeq_Animal))), vars = Market)
count(as.data.frame(as.matrix(sample_data(physeq_Animal))), vars = Type)

# *** TABLE S1 animal taxonomy and relative abundance --------------------------------------------------
otu_animal_abund = taxa_sums(physeq_Animal)/sum(taxa_sums(physeq_LSU_clean))*100
taxa_animal <- as.data.frame(as.matrix(tax_table(physeq_Animal)))
identical(rownames(as.data.frame(otu_animal_abund)), rownames(taxa_animal)) 
taxa_animal$Abundance <- as.vector(otu_animal_abund)
taxa_animal <- taxa_animal[order(taxa_animal$Abundance, decreasing = TRUE),] 

taxa_animal

# ******************************************************---------------------------------------
# Extracting data for 4 markets only ----------------------------------------------------------
FilterMarket <- function(physeq){
  require(dplyr)
  physeq1 <- subset_samples(physeq, Market%in%c("Ciba","Longjin", "Wanyao","New_District"))
  otu_table(physeq1) <- otu_table(physeq1)[which(
    rowSums(otu_table(physeq1)) > 0),] 
  sample_data(physeq1)$Market <- sample_data(physeq1)$Market %>% 
    recode_factor('Ciba' = "Ciba",
                  'Longjin' = "Longjin",
                  'Wanyao' = "Wanyao",
                  "New_District" = "New District")
return(physeq1)
}

FilterMarket(physeq_ITS_clean) -> physeq_ITS_filt
count(as.data.frame(as.matrix(sample_data(physeq_ITS_filt))), vars = Market)
sort(sample_sums(physeq_ITS_filt))

FilterMarket(physeq_16s_clean) -> physeq_16s_filt
count(as.data.frame(as.matrix(sample_data(physeq_16s_filt))), vars = Market)
sort(sample_sums(physeq_16s_filt))

FilterMarket(physeq_LSU) -> physeq_LSU_filt # Using just the fungi here
count(as.data.frame(as.matrix(sample_data(physeq_LSU_filt))), vars = Market)
sort(sample_sums(physeq_LSU_filt))

FilterMarket(physeq_Animal) -> physeq_Ani_filt
# removing samples sum = 0
subset_samples(physeq_Ani_filt, sample_sums(physeq_Ani_filt) >0) -> physeq_Ani_filt
count(as.data.frame(as.matrix(sample_data(physeq_Ani_filt))), vars = Market)
sort(sample_sums(physeq_Ani_filt))
physeq_Ani_filt


# otu_tables results --------------------------------------------------
sum(sample_sums(physeq_ITS_filt))
mean(sample_sums(physeq_ITS_filt))
sd(sample_sums(physeq_ITS_filt))

sum(sample_sums(physeq_LSU_filt))
mean(sample_sums(physeq_LSU_filt))
sd(sample_sums(physeq_LSU_filt))

sum(sample_sums(physeq_16s_filt))
mean(sample_sums(physeq_16s_filt))
sd(sample_sums(physeq_16s_filt))

sum(sample_sums(physeq_Ani_filt))
mean(sample_sums(physeq_Ani_filt))
sd(sample_sums(physeq_Ani_filt))


# ----------------------------------------------------------------------------------------------
# Plottong samples reads 
library(ggrepel)
library(tidyr)
library(grid)
library(gridExtra)


PlotBox <- function(dataframe){
require(ggrepel)
  sample_data(dataframe)$LibrarySize <- sample_sums(dataframe)
  meta <- as.data.frame(as.matrix(sample_data(dataframe)))
  meta$LibrarySize <- as.integer(meta$LibrarySize)
  meta$label <- paste(rownames(meta), meta$LibrarySize, sep = " ")
  
  plot_box <- ggplot(meta, aes(x=reorder(Market, LibrarySize), y = LibrarySize)) + 
          geom_point() + 
          geom_boxplot(outlier.colour="red", outlier.shape=16, outlier.size=3) +
          geom_text_repel(aes(label=label), size = 2) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
        theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
        theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
        theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) +
      labs(x="Market", y="Read number")
return(plot_box)
}

PlotBox(physeq_ITS_filt) + labs(title = "Read distribution of\nITS sample replicates")
PlotBox(physeq_16s_filt) + labs(title = "Read distribution of\n16S sample replicates")
PlotBox(physeq_LSU_filt) + labs(title = "Read distribution of\nLSU sample replicates")
PlotBox(physeq_Ani_filt) + labs(title = "Read distribution of\nAnimal sample replicates")

# I amnot goint to rmove any sample, just trying to rescale 
# read number. This is justified by the fact that most samples with
# low number of reads or that dropped were form the inside
# that is pretty clean of any organism.

PlotBar <- function(physeq){
sample_data(physeq)$Market_Niche <- paste(
                        sample_data(physeq)$Market,
                        sample_data(physeq)$Type,
                        sep = "_")
merge_samples(physeq, "Market_Niche") -> physeq_merg
sample_data(physeq_merg)$Abund <- sample_sums(physeq_merg)
sample_data(physeq_merg)$Market <- c("Ciba","Ciba", "Longjin","Longjin",
                                         "New District","New District", "Wanyao","Wanyao")
sample_data(physeq_merg)$Niche <- c("Inside", "Outside", "Inside", "Outside",
                                        "Inside", "Outside","Inside", "Outside")
meta <- as.data.frame(as.matrix(sample_data(physeq_merg)[, c(2,11,10)]))
meta$Abund <- as.numeric(as.character(meta$Abund))
ggplot(meta, aes(x=Market, y=Abund, fill=Niche)) + 
  labs(y="Cumulative DNA read number") +
  geom_bar(stat="identity", position="dodge") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 8)) +
  scale_fill_manual(values = paletteCB2) -> plot
return(plot)
}

# The warnign is about NA in the sample_data, can be ignored
PlotBar(physeq_ITS_filt) + labs(title = "ITS")
PlotBar(physeq_LSU_filt) + labs(title = "LSU")
PlotBar(physeq_16s_filt) + labs(title = "16S")


# *** FIGURE S2 top - PCoA ordination ---------------------------------------------------------
ggarrange(PlotBar(physeq_ITS_filt) + labs(title = "ITS"),
          PlotBar(physeq_LSU_filt) + labs(title = "LSU"),
          PlotBar(physeq_16s_filt) + labs(title = "16S"),
          labels = c("A","B", "C"),
          widths = c(1, 1, 1),
          align = "hv" ,
          common.legend = TRUE,
          ncol = 3, nrow = 1, 
          legend = c("bottom")) -> Fig_S2_bars

Fig_S2_bars


# RAREFACTION CURVES --------------------------------------------------------------------------
# Sample40 has zero reads that's why the warning!
# converting to table objetcs -----------------------------------------------------------------

sample_sums(physeq_ITS_filt) >0 

otu_fungi_its <- as.data.frame(otu_table(physeq_ITS_filt))
taxa_fungi_its <- as.data.frame(as.matrix(tax_table(physeq_ITS_filt)))
metadata_fungi_its <- as.data.frame(as.matrix(sample_data(physeq_ITS_filt)))
dim(otu_fungi_its)
metadata_fungi_its

sample_sums(physeq_LSU_filt) >0 

otu_fungi_lsu <- as.data.frame(otu_table(physeq_LSU_filt))
taxa_fungi_lsu <- as.data.frame(as.matrix(tax_table(physeq_LSU_filt)))
metadata_fungi_lsu <- as.data.frame(as.matrix(sample_data(physeq_LSU_filt)))
dim(otu_fungi_lsu)
metadata_fungi_lsu

sample_sums(physeq_16s_filt) >0 

otu_prok <- as.data.frame(otu_table(physeq_16s_filt))
taxa_prok <- as.data.frame(as.matrix(tax_table(physeq_16s_filt)))
metadata_prok <- as.data.frame(as.matrix(sample_data(physeq_16s_filt)))
dim(otu_prok)
metadata_prok

# *** FIGURE S3 - Rarefaction curves ******----------------------------------------------------
library(vegan)

# creating and empty character 
curve_colors_market <- vector(mode="character", length=ncol(as.data.frame(otu_table(physeq_ITS_filt))))
curve_colors_market[as.data.frame(as.matrix(sample_data(physeq_ITS_filt)))$Market=="Ciba"] <- "#E69F00"
curve_colors_market[as.data.frame(as.matrix(sample_data(physeq_ITS_filt)))$Market=="Longjin"] <- "#D55E00"
curve_colors_market[as.data.frame(as.matrix(sample_data(physeq_ITS_filt)))$Market=="Wanyao"] <- "#CBD588"
curve_colors_market[as.data.frame(as.matrix(sample_data(physeq_ITS_filt)))$Market=="New District"] <- "#0072B2"
curve_colors_market

par(mfrow=c(1,3)) # Change the panel layout to 3 x 1
rarecurve(t(otu_fungi_its), col = curve_colors_market, label = FALSE, step = 100,
          main="ITS", ylab = "Number of OTUs", xlab = "Number of DNA reads")
legend("topright", legend=c("Ciba", "Longjin", "Wanyao", "New District"),
       col= c("#E69F00", "#D55E00", "#CBD588", "#0072B2"), lty=2, cex=0.8, box.lty=0,lwd=3,bty = "n") 
# if you rarefy you can add at what depth you rarefy here 
#abline(v = min(sample_sums(physeq_fungi_uparse_qc)), col="red", lwd=3, lty=2)
rarecurve(t(otu_fungi_lsu), col = curve_colors_market, label = FALSE,  step = 100,
          main="LSU", ylab = "Number of OTUs", xlab = "Number of DNA reads")
legend("bottomright", legend=c("Ciba", "Longjin", "Wanyao", "New District"),
       col= c("#E69F00", "#D55E00", "#CBD588", "#0072B2"), lty=2, cex=0.8, box.lty=0,lwd=3,bty = "n") 
#abline(v = min(sample_sums(physeq_prok_uparse_qc)), col="red", lwd=3, lty=2)
rarecurve(t(otu_prok), col = curve_colors_market, label = FALSE,  step = 100,
          main="16S", ylab = "Number of OTUs", xlab = "Number of DNA reads")
legend("bottomright", legend=c("Ciba", "Longjin", "Wanyao", "New District"),
       col= c("#E69F00", "#D55E00", "#CBD588", "#0072B2"), lty=2, cex=0.8, box.lty=0,lwd=3,bty = "n") 
#abline(v = min(sample_sums(physeq_prok_uparse_qc)), col="red", lwd=3, lty=2)
dev.off()


# ******************************************************---------------------------------------
# ALPHA DIVERSITY -----------------------------------------------------------------------------
# The maximum value of the Shannon index is log(k), with k denoting the number of groups. 
# This value occurs when each group has the same frequency (i.e., maximum eveness).The Shannon
# equitability index is simply the Shannon diversity index divided by the maximum diversity
# EH=H/log(k). This normalizes the Shannon diversity index to a value between 0 and 1. Note
# that lower values indicate more diversity while higher values indicate less diversity. 
# Specifically, an index value of 1 means that all groups have the same frequency. Here we
# use 1 - EH so that higher values indicate higher diversity.

library(multcompView)
library(vegan)

# extracting main dataframes
ExtrDataFrame <- function(physeq){
  require(vegan)
  # adding aplha metrics to metadata file
  df_alpha <- as.data.frame(as.matrix(sample_data(physeq))) 
  df_alpha$Observed <- specnumber(as.data.frame(otu_table(physeq)), MARGIN = 2)
  df_alpha$Shannon <- diversity(as.data.frame(otu_table(physeq)), index="shannon", MARGIN = 2)
  df_alpha$EH <- 1 - df_alpha$Shannon/log(df_alpha$Observed)
  # re-order factor levels
  df_alpha$Market <- factor(df_alpha$Market,levels=c("Ciba","Longjin", "Wanyao","New District"))
  return(df_alpha)
}


ExtrDataFrame(physeq_ITS_filt) -> df_alpha_ITS
ExtrDataFrame(physeq_16s_filt) -> df_alpha_16s
ExtrDataFrame(physeq_LSU_filt) -> df_alpha_LSU
head(df_alpha_ITS)

# comparing Markets across Type using Wilcox.test
# Since we have done many independent tests (one for each taxon), 
# we should correct for multiple comparisions.
CompSampl <- function(dataframe, level, formula){
  require(multcompView)
    compare_means(formula, data = dataframe[dataframe$Type==level,],
                  p.adjust.method = "bonferroni") -> test_CC 
    test_CC <- as.data.frame(test_CC)[,c(2,3,5)] # to change form p to p.adj do 4 to 5
    test_CC2 <- data.frame(test_CC[,2], test_CC[,1], test_CC[,3])
    colnames(test_CC2) <- c("group1", "group2", "p.adj") # change p to p.adj
    rbind(test_CC, test_CC2) -> test_all
    as.dist(xtabs(test_all[, 3] ~ (test_all[, 2] + test_all[, 1])), diag = TRUE) -> dist_CC
    data.frame(multcompLetters(dist_CC)['Letters']) -> res_CC
      res_CC$sample <- rownames(res_CC)
      res_CC %>% slice(match(c("Ciba", "Longjin","Wanyao","New District"), sample)) -> res_CC
  return(res_CC)
}


CalcType <- function(df, formula){
  CompSampl(df, "in", formula) -> pval_in
  CompSampl(df, "out", formula) -> pval_out
  my_lables <-c(as.character(pval_in$Letters), 
                as.character(pval_out$Letters))
  return(my_lables)
}

CalcType(df_alpha_ITS, formula(Observed ~ Market))-> lables_ITS
CalcType(df_alpha_16s, formula(Observed ~ Market)) -> lables_16s
CalcType(df_alpha_LSU, formula(Observed ~ Market)) -> lables_LSU

CalcType(df_alpha_ITS, formula(EH ~ Market)) -> lables_ITS_shan
CalcType(df_alpha_16s, formula(EH ~ Market)) -> lables_16s_shan
CalcType(df_alpha_LSU, formula(EH ~ Market)) -> lables_LSU_shan

max(df_alpha_ITS[,"EH"] + 0.1* max(df_alpha_ITS[,"EH"])) 

# plotting richness x Markets ------------------------------------------------------------------
PlotRich <- function(dataframe, palette, my_labels, Var){
  #calculating where to put the letters
  dataframe$Type <- dataframe$Type %>% 
    recode_factor('in' = "Inside", 'out'="Outside")
  # label place
  if (Var == "EH"){
    labels_y = 1.1
  }else{
    max(dataframe[,Var] + 0.1* max(dataframe[,Var])) -> labels_y
  }
  # plot
  ggplot(dataframe, aes(x = Market, y = get(Var), color = Type)) +
    geom_boxplot(outlier.colour="black", outlier.size=1.2, outlier.shape = 8, lwd=0.5, width = 0.7) +
    stat_summary(fun=mean, geom="point", shape=18, size=1.9, color="red", fill="red") +
    stat_summary(geom = 'text', label = my_labels, fun= max, aes(y = labels_y), size=3, color="black") +
    facet_grid(~Type, scales = "free_x")+
    theme_classic() +
    expand_limits(y = 0) +
    scale_colour_manual("Origin", values = palette) +
    theme(strip.text.x = element_text(size = 8)) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 8)) +
    theme(axis.text.x = element_text(angle = -90, size = 8, hjust = 0, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
    grids(linetype = "dashed") +
    theme(legend.position="none") -> rich_plot 
  return(rich_plot)
}


PlotRich(df_alpha_ITS, paletteCB2, lables_ITS, "Observed") + labs(title = "ITS", x=NULL, y="Richness")
PlotRich(df_alpha_LSU, paletteCB2, lables_LSU, "Observed") + labs(title = "LSU", x=NULL, y="Richness")
PlotRich(df_alpha_16s, paletteCB2, lables_16s, "Observed") + labs(title = "16S", x=NULL, y="Richness")

PlotRich(df_alpha_ITS, paletteCB2, lables_ITS_shan, "EH") + 
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) + 
  labs(title = "ITS", x=NULL, y="Shannon")
PlotRich(df_alpha_LSU, paletteCB2, lables_LSU_shan, "EH") + 
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) + 
  labs(title = "LSU", x=NULL, y="Shannon")
PlotRich(df_alpha_16s, paletteCB2, lables_16s_shan, "EH") + 
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) + 
  labs(title = "16S", x=NULL, y="Shannon")

# now calculating differences in Richness in vs out
# for each market
CompSampl2 <- function(dataframe, level, formula){
  require(multcompView)
  compare_means(formula, data = dataframe[dataframe$Market==level,],
                p.adjust.method = "bonferroni") -> test_CC 
  test_CC <- as.data.frame(test_CC)[,c(2,3,5)] # to change form p to p.adj do 4 to 5
  test_CC2 <- data.frame(test_CC[,2], test_CC[,1], test_CC[,3])
  colnames(test_CC2) <- c("group1", "group2", "p.adj") # change p to p.adj
  rbind(test_CC, test_CC2) -> test_all
  as.dist(xtabs(test_all[, 3] ~ (test_all[, 2] + test_all[, 1])), diag = TRUE) -> dist_CC
  data.frame(multcompLetters(dist_CC)['Letters']) -> res_CC
  res_CC$sample <- rownames(res_CC)
  res_CC %>% slice(match(c("in", "out"), sample)) -> res_CC
  return(res_CC)
}


Calclabel2<- function(df,formula){
  CompSampl2(df, "Ciba", formula) -> pval_ciba
  CompSampl2(df, "Longjin", formula) -> pval_long
  CompSampl2(df, "Wanyao", formula) -> pval_wan
  CompSampl2(df, "New District", formula) -> pval_newdis
  my_lables <-c(as.character(pval_ciba$Letters), 
                as.character(pval_long$Letters),
                as.character(pval_wan$Letters),
                as.character(pval_newdis$Letters))
  return(my_lables)
}

Calclabel2(df_alpha_ITS, formula(Observed ~ Type)) -> lables_ITS_type
Calclabel2(df_alpha_LSU, formula(Observed ~ Type)) -> lables_LSU_type
Calclabel2(df_alpha_16s, formula(Observed ~ Type)) -> lables_16s_type

Calclabel2(df_alpha_ITS, formula(EH ~ Type)) -> lables_ITS_type_shan
Calclabel2(df_alpha_LSU, formula(EH ~ Type)) -> lables_LSU_type_shan
Calclabel2(df_alpha_16s, formula(EH ~ Type)) -> lables_16s_type_shan


# plotting richness x Markets -----------------------------------------------------------------
PlotRich2 <- function(dataframe, palette, my_labels, Var){
  #calculating where to put the letters
  dataframe$Type <- dataframe$Type %>% 
    recode_factor('in' = "Inside", 'out'="Outside")
  # label place
  if (Var == "EH"){
    labels_y = 1.1
  }else{
    max(dataframe[,Var] + 0.1* max(dataframe[,Var])) -> labels_y
  }
  # plot
  ggplot(dataframe, aes(x = Type, y = get(Var), color = Market)) +
    geom_boxplot(outlier.colour="black", outlier.size=1.2, outlier.shape = 8, lwd=0.5, width = 0.7) +
    stat_summary(fun=mean, geom="point", shape=18, size=1.9, color="red", fill="red") +
    stat_summary(geom = 'text', label = my_labels, fun= max, aes(y = labels_y), size=3, color="black") +
    facet_grid(~Market, scales = "free_x")+
    theme_classic() +
    expand_limits(y = 0) +
    scale_colour_manual("Origin", values = palette) +
    theme(strip.text.x = element_text(size = 8, angle = 90)) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 8)) +
    theme(axis.text.x = element_text(angle = -90, size = 8, hjust = 0, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
    grids(linetype = "dashed") +
    theme(legend.position="none") -> rich_plot 
  return(rich_plot)
}


PlotRich2(df_alpha_ITS, paletteCB4, lables_ITS_type, "Observed") + labs(title = "ITS", x=NULL, y="Richness")
PlotRich2(df_alpha_LSU, paletteCB4, lables_LSU_type, "Observed") + labs(title = "LSU", x=NULL, y="Richness")
PlotRich2(df_alpha_16s, paletteCB4, lables_16s_type, "Observed") + labs(title = "16S", x=NULL, y="Richness")

PlotRich2(df_alpha_ITS, paletteCB4, lables_ITS_type_shan, "EH") + 
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) + 
  labs(title = "ITS", x=NULL, y="Shannon")
PlotRich2(df_alpha_LSU, paletteCB4, lables_LSU_type_shan, "EH") +
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) + 
  labs(title = "LSU", x=NULL, y="Shannon")
PlotRich2(df_alpha_16s, paletteCB4, lables_16s_type_shan, "EH") + 
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) + 
  labs(title = "16S", x=NULL, y="Shannon")


# FIGURE 1 - Alpha diversity -------------------------------------------------------------------
ggarrange(
  ggarrange(PlotRich(df_alpha_ITS, paletteCB2, lables_ITS, "Observed") + labs(title = "ITS", x=NULL, y="Richness"),
            PlotRich(df_alpha_LSU, paletteCB2, lables_LSU, "Observed") + labs(title = "LSU", x=NULL, y=NULL),
            PlotRich(df_alpha_16s, paletteCB2, lables_16s, "Observed") + labs(title = "16S", x=NULL, y=NULL),
            labels = c("A","B","C"),
            widths = c(1,1,1),
            align = "hv" ,
            ncol = 3, 
            nrow = 1, 
            #common.legend = TRUE,
            legend = c("none")),
  ggarrange(PlotRich2(df_alpha_ITS, paletteCB4, lables_ITS_type, "Observed") + labs(title = NULL, x=NULL, y="Richness"),
            PlotRich2(df_alpha_LSU, paletteCB4, lables_LSU_type, "Observed") + labs(title = NULL, x=NULL, y=NULL),
            PlotRich2(df_alpha_16s, paletteCB4, lables_16s_type, "Observed") + labs(title = NULL, x=NULL, y=NULL),
            widths = c(1,1,1),
            align = "hv" ,
            ncol = 3, 
            nrow = 1),
  heights =  c(1,1),
  ncol = 1, 
  nrow = 2) -> Fig_1

Fig_1


ggarrange(
  ggarrange(PlotRich(df_alpha_ITS, paletteCB2, lables_ITS_shan, "EH") +
              scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) + 
              labs(title = "ITS", x=NULL, y="Shannon equitability index"),
            PlotRich(df_alpha_LSU, paletteCB2, lables_LSU_shan, "EH") +
              scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) + 
              labs(title = "LSU", x=NULL, y=NULL),
            PlotRich(df_alpha_16s, paletteCB2, lables_16s_shan, "EH") +
              scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) + 
              labs(title = "16S", x=NULL, y=NULL),
            labels = c("A","B","C"),
            widths = c(1,1,1),
            align = "hv" ,
            ncol = 3, 
            nrow = 1, 
            #common.legend = TRUE,
            legend = c("none")),
  ggarrange(PlotRich2(df_alpha_ITS, paletteCB4, lables_ITS_type_shan, "EH") +
              scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) + 
              labs(title = NULL, x=NULL, y="Shannon equitability index"),
            PlotRich2(df_alpha_LSU, paletteCB4, lables_LSU_type_shan, "EH") +
              scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) + 
              labs(title = NULL, x=NULL, y=NULL),
            PlotRich2(df_alpha_16s, paletteCB4, lables_16s_type_shan, "EH") +
              scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) + 
              labs(title = NULL, x=NULL, y=NULL),
            widths = c(1,1,1),
            align = "hv" ,
            ncol = 3, 
            nrow = 1),
  heights =  c(1,1),
  ncol = 1, 
  nrow = 2) -> Fig_1_2

Fig_1_2


# ggarrange(Fig_1, Fig_1_2,
#           widths = c(1,1),
#           ncol = 2, 
#           nrow = 1) -> Fig_1_all
# 
# Fig_1_all

# *******************************************--------------------------------------------------
# BETA DIVERSITY ------------------------------------------------------------------------------
# trying rescaling OTUs abundances
library(scales)

RescaleOTU <- function(dataframe){
  require(scales)
    otu <- as.data.frame(otu_table(dataframe))
    # Rescaling OTU counts to 0-1
    apply(otu, 1, function(x) rescale(x, to = c(0, 1))) -> otu_scaled
    t(otu_scaled) -> otu_scaled
    otu_table(dataframe) <- otu_table(otu_scaled, taxa_are_rows = TRUE)
return(dataframe)
}

RescaleOTU(physeq_ITS_filt) -> physeq_ITS_scaled
any(sample_sums(physeq_ITS_scaled)==0)
sort(sample_sums(physeq_ITS_scaled))
sort(sample_sums(physeq_ITS_filt))
sort(taxa_sums(physeq_ITS_filt))

RescaleOTU(physeq_16s_filt) -> physeq_16s_scaled
RescaleOTU(physeq_LSU_filt) -> physeq_LSU_scaled
RescaleOTU(physeq_Ani_filt) -> physeq_Ani_scaled

# sample groups -------------------------------------------------------------------------------
table(sample_data(physeq_ITS_scaled)$Market)
table(sample_data(physeq_ITS_scaled)$Type)
table(sample_data(subset_samples(physeq_ITS_scaled, Type=="in"))$Market)
table(sample_data(subset_samples(physeq_ITS_scaled, Type=="out"))$Market)

table(sample_data(physeq_LSU_scaled)$Market)
table(sample_data(physeq_LSU_scaled)$Type)
table(sample_data(subset_samples(physeq_LSU_scaled, Type=="in"))$Market)
table(sample_data(subset_samples(physeq_LSU_scaled, Type=="out"))$Market)

table(sample_data(physeq_16s_scaled)$Market)
table(sample_data(physeq_16s_scaled)$Type)
table(sample_data(subset_samples(physeq_16s_scaled, Type=="in"))$Market)
table(sample_data(subset_samples(physeq_16s_scaled, Type=="out"))$Market)

# Change OTU names ----------------------------------------------------------------------------
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

ChangeName(physeq_ITS_scaled, "F") -> physeq_ITS_scaled
ChangeName(physeq_16s_scaled, "B") -> physeq_16s_scaled
ChangeName(physeq_LSU_scaled, "E") -> physeq_LSU_scaled
ChangeName(physeq_Ani_scaled, "E") -> physeq_Ani_scaled

# calculating PCoA 
pcoa_physeq_ITS_scaled = ordinate(physeq_ITS_scaled, method ="PCoA", distance="bray")
pcoa_physeq_ITS_scaled

pcoa_physeq_16s_scaled = ordinate(physeq_16s_scaled, method ="PCoA", distance="bray")
pcoa_physeq_16s_scaled

pcoa_physeq_LSU_scaled = ordinate(physeq_LSU_scaled, method ="PCoA", distance="bray")
pcoa_physeq_LSU_scaled


# plot ordination 
PlotOrdin <-function(dataframe, ord){
  ord <- plot_ordination(dataframe, ord, color="Market", shape="Type") + 
    geom_point(size=1.7, alpha=0.9, aes(shape=Type)) +
    scale_shape_manual("Niche", values = c(17, 16), labels=c("Inside","Outside")) +
    scale_colour_manual(values=paletteCB4) +
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

PlotOrdin(physeq_ITS_scaled, pcoa_physeq_ITS_scaled) + labs(title = "ITS")
PlotOrdin(physeq_LSU_scaled, pcoa_physeq_LSU_scaled) + labs(title = "LSU")
PlotOrdin(physeq_16s_scaled, pcoa_physeq_16s_scaled) + labs(title = "16S")

# *** FIGURE 2 top - PCoA ordination -----------------------------------------------------------
ggarrange(PlotOrdin(physeq_ITS_scaled, pcoa_physeq_ITS_scaled) + labs(title = "ITS"),
          PlotOrdin(physeq_LSU_scaled, pcoa_physeq_LSU_scaled) + labs(title = "LSU"),
          PlotOrdin(physeq_16s_scaled, pcoa_physeq_16s_scaled) + labs(title = "16S"),
          labels = c("A","B", "C"),
          widths = c(1, 1, 1),
          align = "hv" ,
          common.legend = TRUE,
          ncol = 3, nrow = 1, 
          legend = c("bottom")) -> Fig_2_ord

Fig_2_ord

# ******************************************************----------------------------------------
# PERMANOVA ------------------------------------------------------------------------------------
require(vegan)

otu_ITS <- as.data.frame(otu_table(physeq_ITS_scaled))
metadata_ITS = as(sample_data(physeq_ITS_scaled), "data.frame")
taxa_ITS <- as.data.frame(as.matrix(tax_table(physeq_ITS_scaled)))

otu_16s <- as.data.frame(otu_table(physeq_16s_scaled))
metadata_16s = as(sample_data(physeq_16s_scaled), "data.frame")
taxa_16s <- as.data.frame(as.matrix(tax_table(physeq_16s_scaled)))

otu_LSU <- as.data.frame(otu_table(physeq_LSU_scaled))
metadata_LSU = as(sample_data(physeq_LSU_scaled), "data.frame")
taxa_LSU <- as.data.frame(as.matrix(tax_table(physeq_LSU_scaled)))

otu_Ani <- as.data.frame(otu_table(physeq_Ani_scaled))
metadata_Ani = as(sample_data(physeq_Ani_scaled), "data.frame")
taxa_Ani <- as.data.frame(as.matrix(tax_table(physeq_Ani_scaled)))


# permanova ITS -------------------------------------------------------------------------------
# Test the effect Market 
head(metadata_ITS)

adonis_fungi <- adonis(t(otu_ITS) ~ Market,metadata_ITS, method = "bray",
                       strata=metadata_ITS$Type, permutations=9999)
adonis_fungi
round(p.adjust(adonis_fungi$aov.tab$`Pr(>F)`, "bonferroni"), 4)

# test for the effect of Type 
adonis_fungi1 <- adonis(t(otu_ITS) ~ Type, metadata_ITS, method = "bray",
                        strata=metadata_ITS$Market, permutations=9999)
adonis_fungi1
round(p.adjust(adonis_fungi1$aov.tab$`Pr(>F)`, "bonferroni"), 4)

# test for interactions - in these case I do no include the strata
# since Distance which is at the beginning of the model is treated as a random
# effect 
adonis_fungi2 <- adonis(t(otu_ITS) ~ Market * Type, metadata_ITS, method = "bray",
                        strata=metadata_ITS$Market, permutations=9999)
adonis_fungi2

adonis_fungi2 <- adonis(t(otu_ITS) ~  Type * Market, metadata_ITS, method = "bray",
                        strata=metadata_ITS$Market, permutations=9999)
adonis_fungi2

adonis_fungi3 <- adonis(t(otu_ITS) ~ Type * Market,metadata_ITS, method = "bray",
                        strata=metadata_ITS$Type, permutations=9999)
adonis_fungi3

adonis_fungi3 <- adonis(t(otu_ITS) ~  Market * Type,metadata_ITS, method = "bray",
                        strata=metadata_ITS$Type, permutations=9999)
adonis_fungi3

#               Df  SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Type          1    1.2932   1.29322  4.7078 0.04986 0.0001 ***
# Market        3    5.0024   1.66748  6.0703 0.19286 0.0001 ***
# Type:Market   3    2.0616   0.68720  2.5017 0.07948 0.0001 ***
# Residuals     64   17.5805  0.27469         0.67780           
# Total         71   25.9377                  1.00000   

round(p.adjust(adonis_fungi2$aov.tab$`Pr(>F)`, "bonferroni"), 4)

# permanova 16S -------------------------------------------------------------------------------
# Test the effect Market 
head(metadata_16s)

adonis_bact <- adonis(t(otu_16s) ~ Market,metadata_16s, method = "bray",
                      strata=metadata_16s$Type,permutations=9999)
adonis_bact

# test for the effect of Distance 
adonis_bact1 <- adonis(t(otu_16s) ~ Type,metadata_16s, method = "bray",
                       strata=metadata_16s$Market, permutations=9999)
adonis_bact1

# test for interactions - in these case I do no include the strata
# since Distance which is at the beginning of the model is treated as a random
# effect 
adonis_bact2 <- adonis(t(otu_16s) ~ Market * Type, metadata_16s, method = "bray",
                       strata=metadata_16s$Market, permutations=9999)
adonis_bact2

adonis_bact3 <- adonis(t(otu_16s) ~ Type * Market,metadata_16s, method = "bray",
                       strata=metadata_16s$Type, permutations=9999)
adonis_bact3

#             Df SumsOfSqs MeanSqs F.Model    R2    Pr(>F)    
# Market       3    4.9351 1.64503  5.4226 0.17040 0.0001 ***
# Type         1    1.6443 1.64434  5.4203 0.05678 0.0001 ***
# Market:Type  3    2.9664 0.98881  3.2595 0.10243 0.0001 ***
# Residuals   64   19.4154 0.30337         0.67039           
# Total       71   28.9613                 1.00000         

round(p.adjust(adonis_bact2$aov.tab$`Pr(>F)`, "bonferroni"), 4)

# permanova LSU -------------------------------------------------------------------------------
# Test the effect Market 
head(metadata_LSU)

adonis_lsu <- adonis(t(otu_LSU) ~ Market,metadata_LSU, method = "bray",
                     strata=metadata_LSU$Type,permutations=9999)
adonis_lsu

# test for the effect of Distance 
adonis_lsu1 <- adonis(t(otu_LSU) ~ Type,metadata_LSU, method = "bray",
                      strata=metadata_LSU$Market, permutations=9999)
adonis_lsu1

# test for interactions - in these case I do no include the strata
# since Distance which is at the beginning of the model is treated as a random
# effect 
adonis_lsu2 <- adonis(t(otu_LSU) ~ Market * Type, metadata_LSU, method = "bray",
                      strata=metadata_LSU$Market, permutations=9999)
adonis_lsu2

adonis_lsu3 <- adonis(t(otu_LSU) ~ Type * Market,metadata_LSU, method = "bray",
                      strata=metadata_LSU$Type, permutations=9999)
adonis_lsu3

#             Df SumsOfSqs MeanSqs F.Model    R2    Pr(>F)    
# Market       3    4.7287 1.57623  4.8910 0.15825 0.0001 ***
# Type         1    1.6900 1.69004  5.2442 0.05656 0.0001 ***
# Market:Type  3    2.8369 0.94564  2.9343 0.09494 0.0001 ***
# Residuals   64   20.6253 0.32227         0.69025           
# Total       71   29.8810                 1.00000   

round(p.adjust(adonis_lsu2$aov.tab$`Pr(>F)`, "bonferroni"), 4)

# PAIRWISE PERMANOVA --------------------------------------------------------------------------
library(mctoolsr) # post hoc test for permanova and other tools

mctoolsr::calc_pairwise_permanovas(otu_ITS, metadata_ITS, "Market")
mctoolsr::calc_pairwise_permanovas(otu_ITS, metadata_ITS, "Type")

mctoolsr::calc_pairwise_permanovas(otu_LSU, metadata_LSU, "Market")
mctoolsr::calc_pairwise_permanovas(otu_LSU, metadata_LSU, "Type")

mctoolsr::calc_pairwise_permanovas(otu_16s, metadata_16s, "Market")
mctoolsr::calc_pairwise_permanovas(otu_16s, metadata_16s, "Type")


# ******************************************************---------------------------------------

# MULTIVARIATE DISPERSION ---------------------------------------------------------------------
# Test the homogenity of group variances
permdisp_ITS_market <- betadisper(vegdist(t(otu_ITS), method="bray"), metadata_ITS$Market) 
anova(permdisp_ITS_market, permutations = 9999)
permutest(permdisp_ITS_market, permutations = 9999, pairwise = T) -> dist_ITS_market
round(p.adjust(anova(permdisp_ITS_market, permutations = 9999)$`Pr(>F)`, "bonferroni"), 4) # p=0.00472

permdisp_ITS_type <- betadisper(vegdist(t(otu_ITS), method="bray"), metadata_ITS$Type) 
anova(permdisp_ITS_type, permutations = 9999)
permutest(permdisp_ITS_type, permutations = 9999, pairwise = T) -> dist_ITS_type
round(p.adjust(anova(permdisp_ITS_type, permutations = 999)$`Pr(>F)`, "bonferroni"), 4) # p= 0.0074
  
  
permdisp_LSU_market <- betadisper(vegdist(t(otu_LSU), method="bray"), metadata_LSU$Market) 
anova(permdisp_LSU_market, permutations = 9999)
permutest(permdisp_LSU_market, permutations = 9999, pairwise = T)-> dist_LSU_market
round(p.adjust(anova(permdisp_LSU_market, permutations = 999)$`Pr(>F)`, "bonferroni"), 4)
  
permdisp_LSU_type <- betadisper(vegdist(t(otu_LSU), method="bray"), metadata_LSU$Type) 
anova(permdisp_LSU_type, permutations = 9999)
permutest(permdisp_LSU_type, permutations = 9999, pairwise = T) -> dist_LSU_type
round(p.adjust(anova(permdisp_LSU_type, permutations = 9999)$`Pr(>F)`, "bonferroni"), 4)
  

permdisp_16s_market <- betadisper(vegdist(t(otu_16s), method="bray"), metadata_16s$Market) 
anova(permdisp_16s_market, permutations = 9999)
permutest(permdisp_16s_market, permutations = 9999, pairwise = T)-> dist_16s_market
round(p.adjust(anova(permdisp_16s_market, permutations = 9999)$`Pr(>F)`, "bonferroni"), 4) #p=0.00004
  
permdisp_16s_type <- betadisper(vegdist(t(otu_16s), method="bray"), metadata_16s$Type) 
anova(permdisp_16s_type, permutations = 9999)
permutest(permdisp_16s_type, permutations = 9999, pairwise = T) -> dist_16s_type
round(p.adjust(anova(permdisp_16s_type, permutations = 9999)$`Pr(>F)`, "bonferroni"), 4)
  
  
multcompLetters(p.adjust(dist_ITS_market$pairwise$observed, method="bonferroni"))['Letters']

data.frame(multcompLetters(p.adjust(dist_ITS_market$pairwise$observed, 
                                      method="bonferroni"))['Letters']) -> pair_ITS_mark
data.frame(multcompLetters(p.adjust(dist_ITS_type$pairwise$observed,
                                      method="bonferroni"))['Letters']) -> pair_ITS_type
data.frame(multcompLetters(p.adjust(dist_16s_market$pairwise$observed
                                      , method="bonferroni"))['Letters']) -> pair_16s_market
  
# max(dataframe$Observed + 0.1* max(dataframe$Observed)) -> labels_y
# stat_summary(geom = 'text', label = my_labels, fun= max, aes(y = labels_y), size=3, color="black") +

# plot multivariate dispersion 
PlotBetadisper <- function(betadisp, value, palette, my_labels, metadata){
    # creating a label and dataframe
    max(betadisp$distances + 0.1* max(betadisp$distances)) -> labels_y
        #labels_y = 0.9
    data.frame(betadisp$group, betadisp$distances) -> df
    colnames(df) <- c("value", "distance")
    if (identical(rownames(df), rownames(metadata))==TRUE){
      df$Type <- metadata$Type
      # plotting
      ggplot(df, aes(x=value, y=distance, color=value)) +
          geom_jitter(size=1.2, alpha=0.5, aes(shape=Type)) +
            geom_boxplot(outlier.colour="black", outlier.shape = 8, alpha=0.6, lwd = 0.7) +
        stat_summary(fun=mean, geom="point", shape=18, size=1.9, color="red", fill="red") +
        stat_summary(geom = 'text', label = my_labels, fun= max, aes(y = labels_y), size=3, color="black") +
        theme_classic() +
        scale_colour_manual("Origin", values = palette) +
          theme(strip.text.x = element_text(size = 9, face = "bold")) +
          theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
          theme(axis.text.x = element_text(angle = -90, size = 8, hjust = 0, vjust = 0.5)) +
          theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
          theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
        grids(linetype = "dashed") +
        theme(legend.position="none") -> betaplot
    return(betaplot)
    }else{
      stop("Error: dataframe and metadata rownames are not matching!")
    }
}
  

PlotBetadisper(permdisp_ITS_market, "Market", paletteCB4, as.character(pair_ITS_mark$Letters), metadata_ITS) + 
    labs(title = "ITS\n(Market)", y ="Distance to centroids", x= "Market")

PlotBetadisper(permdisp_ITS_type, "Type",paletteCB2, as.character(pair_ITS_type$Letters), metadata_ITS) + 
    labs(title = "ITS\n(Origin)", y ="Distance to centroids", x= "Origin") +
    scale_x_discrete(breaks=c("in", "out"), labels=c("inside", "outside"))

  
# #*** FIGURE 2 bottom - betadisper ----------------------------------------------------------
# ggarrange(PlotBetadisper(permdisp_ITS_market, "Market", paletteCB4, as.character(pair_ITS_mark$Letters)) +
#               labs(title = "ITS\n(Market)", y ="Distance to centroids", x= "Market"),
#             PlotBetadisper(permdisp_ITS_type, "Type", paletteCB2, as.character(pair_ITS_type$Letters)) +
#               labs(title = "ITS\n(Origin)", y ="Distance to centroids", x= "Origin") +
#               scale_x_discrete(breaks=c("in", "out"), labels=c("inside", "outside")),
#             PlotBetadisper(permdisp_LSU_market, "Market", paletteCB4, my_labels = "") +
#               labs(title = "LSU\n(Market)", y ="Distance to centroids", x= "Market"),
#             PlotBetadisper(permdisp_LSU_type, "Type", paletteCB2, as.character(pair_LSU_type$Letters)) +
#               labs(title = "LSU\n(Origin)", y ="Distance to centroids", x= "Origin") +
#               scale_x_discrete(breaks=c("in", "out"), labels=c("inside", "outside")),
#             PlotBetadisper(permdisp_16s_market, "Market", paletteCB4, as.character(pair_16s_market$Letters)) +
#               labs(title = "16S\n(Market)", y ="Distance to centroids", x= "Market"),
#             PlotBetadisper(permdisp_16s_type, "Type", paletteCB2, my_labels = "") +
#               labs(title = "16S\n(Origin)", y ="Distance to centroids", x= "Origin") +
#               scale_x_discrete(breaks=c("in", "out"), labels=c("inside", "outside")),
#     labels = c("A","","B","","C",""),
#     widths = c(1,1,1,1,1,1),
#     align = "hv",
#     ncol = 6,
#     nrow = 1) -> betadisp_plot
  
ggarrange(PlotBetadisper(permdisp_ITS_market, "Market", paletteCB4, as.character(pair_ITS_mark$Letters), metadata_ITS) + 
              labs(title = "Market", y ="Distance to centroids", x= NULL) +
            annotate("text", Inf, -Inf, label = "italic(p) == 0.0047", parse=TRUE, size = 2.15, hjust = 1, vjust = -0.1),
          PlotBetadisper(permdisp_ITS_type, "Type", paletteCB2, as.character(pair_ITS_type$Letters), metadata_ITS) + 
              labs(title = "Niche", y =NULL, x= NULL) +
            annotate("text", Inf, -Inf, label = "italic(p) == 0.0074", parse=TRUE, size = 2.15, hjust = 1, vjust = -0.1)+
              scale_x_discrete(breaks=c("in", "out"), labels=c("Inside", "Outside")),
          PlotBetadisper(permdisp_LSU_market, "Market", paletteCB4, my_labels = "", metadata_LSU) + 
              labs(title = "Market", y ="Distance to centroids", x= NULL),
          PlotBetadisper(permdisp_LSU_type, "Type", paletteCB2, my_labels = "", metadata_LSU) + 
              labs(title = "Niche", y =NULL, x= NULL) +
              scale_x_discrete(breaks=c("in", "out"), labels=c("Inside", "Outside")),
          PlotBetadisper(permdisp_16s_market, "Market", paletteCB4, as.character(pair_16s_market$Letters), metadata_16s) + 
              labs(title = "Market", y ="Distance to centroids", x= NULL) +
            annotate("text", Inf, -Inf, label = "italic(p) == 0.0001", parse=TRUE, size = 2.15, hjust = 1, vjust = -0.1),
          PlotBetadisper(permdisp_16s_type, "Type", paletteCB2, my_labels = "", metadata_16s) + 
              labs(title = "Niche", y =NULL, x= NULL) +
              scale_x_discrete(breaks=c("in", "out"), labels=c("Inside", "Outside")),
            widths = c(1,1,1,1,1,1),
            labels = c("D","","E","","F",""),
            align = "hv",
            ncol = 6, 
            nrow = 1) -> betadisp_plot
  
betadisp_plot
  
# *** FIGURE 2 - complete -------------------------------------------------------------------
ggarrange(Fig_2_ord, betadisp_plot,
            heights = c(1, 0.8),
            align = "hv",
            ncol = 1, 
            nrow = 2, 
            legend = c("bottom")) -> Fig2_beta_div
  
Fig2_beta_div


# COMBNING DATAFRAMES -------------------------------------------------------------------------
#http://adv-r.had.co.nz/Functions.html#lexical-scoping
# ufoMerged <- do.call("rbind", list(ufo1, ufo2, ufo3, ufo4))

BindData <- function(x, y, type) {
  if (type == "otu" | type == "taxa"){
     x.diff <- setdiff(colnames(x), colnames(y))
      y.diff <- setdiff(colnames(y), colnames(x))
     x[, c(as.character(y.diff))] <- NA
     y[, c(as.character(x.diff))] <- NA
     rbind(x, y) -> dataframe
      return(dataframe)
  }else {
    x.diff <- setdiff(rownames(x), rownames(y))
    y.diff <- setdiff(rownames(y), rownames(x))
    x[c(as.character(y.diff)),] <- NA
    y[c(as.character(x.diff)),] <- NA
    cbind(x, y) -> dataframe
    return(dataframe[,1:8])
  }
}

# BindData(otu_ITS, otu_16s, "otu") -> otu_merged
# BindData(taxa_ITS, taxa_16s, "taxa") -> taxa_merged
# BindData(metadata_ITS, metadata_16s, "meta") -> meta_merged

# I think I will just leave ITS, LSU 9clean) and 16S for the model --------------------------------------
# and will remove the non-fungal OTUs since are non-target amplifications

MergeDataframe <- function(x, y, z, type){
  BindData(x, y, type = type) -> dataframe1
  BindData(dataframe1, z, type = type) -> dataframe2
  #BindData(dataframe2, v, type = type) -> dataframe3
     return(dataframe2)
}

# going to keep Animals/Rhizaria out as they are non target PCR amplifications
MergeDataframe(otu_16s, otu_ITS, otu_LSU, "otu") -> otu_merged
dim(otu_merged)

MergeDataframe(taxa_16s, taxa_ITS, taxa_LSU[, 2:8], "taxa") -> taxa_merged
dim(taxa_merged)

MergeDataframe(metadata_16s, metadata_ITS, metadata_LSU, "meta") -> meta_merged
dim(meta_merged)

# recreating the object with filtered data ----------------------------------------------------
physeq_all <- merge_phyloseq(otu_table(otu_merged, taxa_are_rows = TRUE),
                               tax_table(as.matrix(taxa_merged)),
                               sample_data(meta_merged))

physeq_all

otu_table(physeq_all)[is.na(tax_table(physeq_all))]<-0
tax_table(physeq_all)[is.na(tax_table(physeq_all))]<-"Unclassified"

# Tesing a PCoA on merged data ----------------------------------------------------------------
pcoa_physeq_all = ordinate(physeq_all, method ="PCoA", distance="bray")
pcoa_physeq_all

PlotOrdin(physeq_all, pcoa_physeq_all) + labs(title = "Merged")

# Testing permanova on merged data ------------------------------------------------------------
# Test the effect Market 
adonis_merged <- adonis(t(otu_merged) ~ Market, meta_merged, method = "bray",
                       strata=meta_merged$Type, permutations=9999)
adonis_merged
round(p.adjust(adonis_fungi$aov.tab$`Pr(>F)`, "bonferroni"), 4)

# test for the effect of Distance 
adonis_merged1 <- adonis(t(otu_merged) ~ Type, meta_merged, method = "bray",
                        strata=meta_merged$Market, permutations=9999)
adonis_merged1
round(p.adjust(adonis_fungi1$aov.tab$`Pr(>F)`, "bonferroni"), 4)

# test for interactions - in these case I do no include the strata
# since Distance which is at the beginning of the model is treated as a random
# effect 
adonis_merged2 <- adonis(t(otu_merged) ~ Market * Type, meta_merged, method = "bray",
                        strata=meta_merged$Market, permutations=9999)
adonis_merged2

adonis_merged3 <- adonis(t(otu_merged) ~ Type * Market,meta_merged, method = "bray",
                        strata=meta_merged$Type, permutations=9999)
adonis_merged3


# REFORMAT TAXONOMY --------------------------- 
taxa_ITS <- as.data.frame(as.matrix(tax_table(physeq_ITS_scaled)))[, 1:7]

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

# **************************************-------------------------------------------------------
# PROCRUSTES ROTATION ANALYSIS ----------------------------------------------------------------
library(ade4)
library(vegan)

# Extract vectors and incorporate sample data -------------------------------------------------
pcoa_its_df <- data.frame(pcoa_physeq_ITS_scaled$vectors)
pcoa_lsu_df <- data.frame(pcoa_physeq_LSU_scaled$vectors)
pcoa_16s_df <- data.frame(pcoa_physeq_16s_scaled$vectors)

dim1_its <- dim(sample_data(physeq_ITS_scaled))[2] + 1
dim2_its <-
  as.numeric(dim(sample_data(physeq_ITS_scaled))[2] + dim(pcoa_its_df)[2])
dim1_lsu <- dim(sample_data(physeq_LSU_scaled))[2] + 1
dim2_lsu <-
  as.numeric(dim(sample_data(physeq_LSU_scaled))[2] + dim(pcoa_lsu_df)[2])
dim1_16s <- dim(sample_data(physeq_16s_scaled))[2] + 1
dim2_16s <-
  as.numeric(dim(sample_data(physeq_16s_scaled))[2] + dim(pcoa_16s_df)[2])

sample_data(physeq_ITS_scaled)[, c(dim1_its:dim2_its)] <-
  pcoa_its_df
sample_data(physeq_LSU_scaled)[, c(dim1_lsu:dim2_lsu)] <-
  pcoa_lsu_df
sample_data(physeq_16s_scaled)[, c(dim1_16s:dim2_16s)] <-
  pcoa_16s_df

ord_its <- data.frame(sample_data(physeq_ITS_scaled))
ord_lsu <- data.frame(sample_data(physeq_LSU_scaled))
ord_16s <- data.frame(sample_data(physeq_16s_scaled))

# Order both tables by sample ID ---------------------------------------------------------------
ord_its <- ord_its[order(row.names(ord_its)),]
ord_lsu <- ord_lsu[order(row.names(ord_lsu)),]
ord_16s <- ord_16s[order(row.names(ord_16s)),]

identical(rownames(ord_its), rownames(ord_lsu)) # double check
identical(rownames(ord_its), rownames(ord_16s))

# Run procrustes analysis ----------------------------------------------------------------------
pro_its_lsu <- protest(ord_its[, c("Axis.1", "Axis.2")],
                       ord_lsu[, c("Axis.1", "Axis.2")],
                       scores = "sites",
                       permutations = how(nperm = 9999))
pro_its_lsu
pro_its_lsu$ss

pro_its_16s <- protest(ord_its[, c("Axis.1", "Axis.2")],
                       ord_16s[, c("Axis.1", "Axis.2")],
                       scores = "sites",
                       permutations = how(nperm = 9999))
pro_its_16s
pro_its_16s$ss

pro_16s_lsu <- protest(ord_lsu[, c("Axis.1", "Axis.2")],
                       ord_16s[, c("Axis.1", "Axis.2")],
                       scores = "sites",
                       permutations = how(nperm = 9999))
pro_16s_lsu
pro_16s_lsu$ss


# Plotting -------------------------------------------------------------------------------------
PlotProcrust <- function(pro_df, ord){
  pro_plot <- data.frame(
    pcoa1 = pro_df$Yrot[, 1],
    pcoa2 = pro_df$Yrot[, 2],
    xpcoa1 = pro_df$X[, 1],
    xpcoa2 = pro_df$X[, 2])
  if (identical(rownames(pro_plot), rownames(ord))==TRUE){
  pro_plot$Market <- ord$Market
  pro_plot$Niche <- ord$Type
  ggplot(pro_plot) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
    theme(strip.text.x = element_text(size = 9, face = "bold")) +
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 9, face = "bold")) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.05)) + # to adjust decimals
    scale_y_continuous(labels = scales::number_format(accuracy = 0.05)) +
    grids(linetype = "dashed") +
    theme(panel.background = element_blank(),
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line()) +
  geom_point(aes(x=pcoa1, y=pcoa2, colour=Market, shape=Niche), size=1.5, alpha=0.9) +
  geom_point(aes(x=xpcoa1, y=xpcoa2, colour=Market, shape=Niche), size=1.5, alpha=0.9) +
  geom_segment(aes(x=xpcoa1,y=xpcoa2,xend=pcoa1,yend=pcoa2,colour=Market), size=0.5,
               arrow=arrow(length=unit(0.3,"cm")), show.legend=F) +
  geom_hline(yintercept = 0,linetype="dashed") +
  geom_vline(xintercept = 0,linetype="dashed") +
  geom_abline(slope = -1/pro_df$rotation[1,2]) +
  geom_abline(slope = pro_df$rotation[1,2]) +
    scale_colour_manual(values=paletteCB4, name = "Market") +
    scale_shape_manual(values = c(17, 16), 
                       labels=c("Inside","Outside"),
                       name="Niche") -> plot_proc
return(plot_proc)
}else{
stop("Error: pro_plot and ord rownames are not matching!")
  }
}


PlotProcrust(pro_its_lsu, ord_its) + 
  labs(title = "Procrustes Rotation\nITS - LSU", x="Axis.1", y="Axis.2") +
annotate("text", Inf, -Inf, label = "italic(m) ^ 2 == 0.622", parse = TRUE, size = 3, hjust = 1.1, vjust = -0.5) +
annotate("text", Inf, -Inf, label = "italic(r) ^ 2 == 0.581", parse = TRUE, size = 3, hjust = 1.1, vjust = -1.6)

PlotProcrust(pro_its_16s, ord_its) +
  labs(title = "Procrustes Rotation\nITS - 16S", x="Axis.1", y="Axis.2") +
  annotate("text", Inf, -Inf, label = "italic(m) ^ 2 == 0.843", parse = TRUE, size = 3, hjust = 1.1, vjust = -0.5)

PlotProcrust(pro_16s_lsu, ord_16s) + 
  labs(title = "Procrustes Rotation\nLSU - 16S", x="Axis.1", y="Axis.2") +
  annotate("text", Inf, -Inf, label = "italic(m) ^ 2 == 0.742", parse = TRUE, size = 3, hjust = 1.1, vjust = -0.5)

# *** FIGURE 4 Top - Procustes Rotation Ordinations ------------------------------------------------
ggarrange(PlotProcrust(pro_its_lsu, ord_its) + 
            labs(title = "Procrustes Rotation\nITS - LSU", x="Axis.1", y="Axis.2") +
            annotate("text", Inf, -Inf, label = "italic(m) ^ 2 == 0.622", parse = TRUE, size = 3, hjust = 1.1, vjust = -0.5) +
            annotate("text", Inf, -Inf, label = "italic(r) == 0.581", parse = TRUE, size = 3, hjust = 1.3, vjust = -2.2),
          PlotProcrust(pro_its_16s, ord_its) +
            labs(title = "Procrustes Rotation\nITS - 16S", x="Axis.1", y="Axis.2") +
            annotate("text", Inf, -Inf, label = "italic(m) ^ 2 == 0.843", parse = TRUE, size = 3, hjust = 1.4, vjust = -0.5) +
            annotate("text", Inf, -Inf, label = "italic(r) == 0.396", parse = TRUE, size = 3, hjust = 1.7, vjust = -2.2),
          PlotProcrust(pro_16s_lsu, ord_16s) + 
            labs(title = "Procrustes Rotation\nLSU - 16S", x="Axis.1", y="Axis.2") +
            annotate("text", Inf, -Inf, label = "italic(m) ^ 2 == 0.742", parse = TRUE, size = 3, hjust = 1.1, vjust = -0.5) +
            annotate("text", Inf, -Inf, label = "italic(r) == 0.508", parse = TRUE, size = 3, hjust = 1.3, vjust = -2.2),
          labels = c("A","B", "C"),
          widths = c(1, 1, 1),
          align = "hv" ,
          common.legend = TRUE,
          ncol = 3, nrow = 1, 
          legend = c("bottom")) -> Fig_4_ord

Fig_4_ord

# m2 statistic provides an overall measure of the concordance between two data sets 
# (i.e. how close the two data configuration match),
# Procrustes analysis allows one to determine how much variance in one matrix is 
# attributable to the variance in the other. Further, visual inspection of a 
# Procrustes plot, in which the residuals between points from each matrix are
# mapped, can allow the identification of individual objects that have (relatively)
# unusual concordance.

# getting ggplot font colors 
calc_element("axis.text.x", theme_classic())
calc_element("axis.line", theme_classic())

# Residual Plot - This allows identification of samples with the worst fit. The horizontal
# lines, from bottom to top, are the 25% (dashed), 50% (solid), and 75% (dashed) quantiles 
# of the residuals. 

# General function for residual plot exploration
ProcrustRes <- function(pro_df, ord){
as.data.frame(residuals(pro_df)) -> res_df
colnames(res_df) <- "Res"
res_df$Market <- ord$Market
res_df$Niche <- ord$Type # use this to map Niche
if (identical(rownames(res_df), rownames(ord))==TRUE){
res_df <- arrange(res_df, Market)
res_df$index <- as.numeric(1:72)
res_df$sample_ID <- rownames(res_df)
res_df[, "sample_ID"] <- gsub("sample", "", res_df[, "sample_ID"])
res_df$Market <- factor(res_df$Market,levels=c("Ciba","Longjin", "Wanyao","New District"))
res_df$sample_ID <- factor(res_df$sample_ID, levels=res_df$sample_ID)
ggplot(res_df, aes(x=reorder(sample_ID, Market), y=Res, fill = Niche)) + # change to Niche for sum by market
  geom_bar(stat="identity")+
  theme_classic() +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
  theme(strip.text.x = element_text(size = 9, face = "bold")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, size = 5.5, hjust = 1, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 9, face = "bold")) +
  ylim(0, 0.20) +
  labs(y="Procrustes residual", x="Samples") +
  geom_hline(yintercept = quantile(res_df$Res)[2], linetype="dashed", color="grey30") +
  geom_hline(yintercept = quantile(res_df$Res)[3], linetype="solid", color="grey30") +
  geom_hline(yintercept = quantile(res_df$Res)[4], linetype="dashed", color="grey30") +
  theme(legend.position = "bottom") -> re_plot
return(re_plot)
}else{
  stop("Error: pro_plot and ord rownames are not matching!")
  }
}

ProcrustRes(pro_16s_lsu, ord_its) + scale_fill_manual(values=paletteCB2) +
  facet_grid(~Market, scales = "free_x", space="free_x") 

# *** FIGURE S4 - Procustes errors residuals ---------------------------------------------------
ggarrange(ProcrustRes(pro_its_lsu, ord_its) + 
            scale_fill_manual(values=paletteCB2) +
            facet_grid(~Market, scales = "free_x", space="free_x") + 
            labs(title = "Procrustes error ITS - LSU"), 
          ProcrustRes(pro_its_16s, ord_its) + 
            scale_fill_manual(values=paletteCB2) +
            facet_grid(~Market, scales = "free_x", space="free_x") +
            labs(title = "Procrustes error ITS - 16S"),
          ProcrustRes(pro_16s_lsu, ord_its) + 
            scale_fill_manual(values=paletteCB2) +
            facet_grid(~Market, scales = "free_x", space="free_x") +
            labs(title = "Procrustes error LSU - 16S"),
          labels = c("A","B", "C"),
          widths = c(1, 1, 1),
          align = "hv" ,
          common.legend = TRUE,
          ncol = 1, nrow = 3, 
          legend = c("bottom")) -> Fig_S4_proc

Fig_S4_proc

# # if you use Market as a variable
# ProcrustRes(pro_its_lsu, ord_its) + scale_fill_manual(values=paletteCB4)
# ProcrustRes(pro_its_16s, ord_its) + scale_fill_manual(values=paletteCB4)
# ProcrustRes(pro_16s_lsu, ord_its) + scale_fill_manual(values=paletteCB4)

# ANALYSIS of PROCRUSTES ERRORS ---------------------------------------------------------------

# extracting dataset
ExtrProcDf <- function(procdf, map){
res_df <- as.data.frame(residuals(procdf))
colnames(res_df) <- "Res_Err"
res_df$Sample_ID <- rownames(res_df)
identical(rownames(map), rownames(res_df))
order_proc_its_lsu <- match(rownames(map), rownames(res_df))
res_df <- res_df[order_proc_its_lsu,]
res_df$Market <- map$Market
res_df$Type <- map$Type
#res_df$Type <- res_df$Type %>% 
#  recode_factor('in' = "Inside", 'out'="Outside")
return(res_df)
}

ExtrProcDf(pro_its_lsu, metadata_ITS) -> res_its_lsu
ExtrProcDf(pro_its_16s, metadata_ITS) -> res_its_16s
ExtrProcDf(pro_16s_lsu, metadata_ITS) -> res_16s_lsu

compare_means(Res_Err ~ Market, res_its_lsu, method = "kruskal.test", p.adjust.method = "bonferroni") # p.adj = 0.012 *
compare_means(Res_Err ~ Type, res_its_lsu, method = "kruskal.test", p.adjust.method = "bonferroni") # p.adj = 0.24 ns

compare_means(Res_Err ~ Market, res_its_16s, method = "kruskal.test", p.adjust.method = "bonferroni") # p.adj = 0.00019 ***
compare_means(Res_Err ~ Type, res_its_16s, method = "kruskal.test", p.adjust.method = "bonferroni") # p.adj = 0.00034 ***

compare_means(Res_Err ~ Market, res_16s_lsu, method = "kruskal.test", p.adjust.method = "bonferroni") # p.adj = 0.00016 ***
compare_means(Res_Err ~ Type, res_16s_lsu, method = "kruskal.test", p.adjust.method = "bonferroni") # p.adj = 0.93 ns

# calculating residuals means/median ----------------------------------------------------------

ExtProcLab <- function(dataframe, formula, method){
  require(multcompView)
  compare_means(formula, data = dataframe, method = method, p.adjust.method = "bonferroni") -> test_CC 
  test_CC <- as.data.frame(test_CC)[,c(2,3,5)] # to change form p to p.adj do 4 to 5
  test_CC2 <- data.frame(test_CC[,2], test_CC[,1], test_CC[,3])
  colnames(test_CC2) <- c("group1", "group2", "p.adj") # change p to p.adj
  rbind(test_CC, test_CC2) -> test_all
  as.dist(xtabs(test_all[, 3] ~ (test_all[, 2] + test_all[, 1])), diag = TRUE) -> dist_CC
  data.frame(multcompLetters(dist_CC)['Letters']) -> res_CC
  res_CC$sample <- rownames(res_CC)
  if (formula == formula(Res_Err ~ Market)){
    res_CC %>% slice(match(c("Ciba", "Longjin","Wanyao","New District"), sample)) -> res_CC
    as.character(res_CC$Letters) -> res_CC
  }else{
    res_CC %>% slice(match(c("in", "out"), sample)) -> res_CC
    as.character(res_CC$Letters) -> res_CC
  }
 return(res_CC)
}

ExtProcLab(res_its_lsu, formula(Res_Err ~ Market), "wilcox.test") -> labels_its_lsu_mark
ExtProcLab(res_its_lsu, formula(Res_Err ~ Type), "wilcox.test") -> labels_its_lsu_niche

ExtProcLab(res_its_16s, formula(Res_Err ~ Market), "wilcox.test") -> labels_its_16s_mark
ExtProcLab(res_its_16s, formula(Res_Err ~ Type), "wilcox.test") -> labels_its_16s_niche

ExtProcLab(res_16s_lsu, formula(Res_Err ~ Market), "wilcox.test") -> labels_16s_lsu_mark
ExtProcLab(res_16s_lsu, formula(Res_Err ~ Type), "wilcox.test") -> labels_16s_lsu_niche

# plotting functions for Residuals ------------------------------------------------------------
ProcrustResAver <- function(dataframe, palette, my_labels, Var1, Var2){
  #calculating where to put the letters
  dataframe$Type <- dataframe$Type %>% 
    recode_factor('in' = "Inside", 'out'="Outside")
  # label place
  if (Var1 == "EH"){
    labels_y = 1.1
  }else{
    max(dataframe[,Var1] + 0.1* max(dataframe[,Var1])) -> labels_y
  }
  # plot
  ggplot(dataframe, aes(x = get(Var2), y = get(Var1), color = get(Var2))) +
    geom_boxplot(outlier.colour="black", outlier.size=1.2, outlier.shape = 8, lwd=0.5, width = 0.7) +
    stat_summary(fun=mean, geom="point", shape=18, size=1.9, color="red", fill="red") +
    stat_summary(geom = 'text', label = my_labels, fun= max, aes(y = labels_y), size=3, color="black") +
        geom_jitter(size=1.2, alpha=0.5, aes(shape=Type)) +
    #facet_grid(~Type, scales = "free_x")+
    theme_classic() +
    expand_limits(y = 0) +
      scale_colour_manual("Origin", values = palette) +
    theme(strip.text.x = element_text(size = 8)) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 8)) +
    theme(axis.text.x = element_text(angle = -90, size = 8, hjust = 0, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
    grids(linetype = "dashed") +
    theme(legend.position="none") -> rich_plot 
  return(rich_plot)
}


ProcrustResAver(res_its_lsu, paletteCB4, labels_its_lsu_mark, "Res_Err", "Market") + 
  labs(title = "ITS-LSU", x=NULL, y="Procrustes residual") 

# *** FIGURE 4 Bottom - Procrustes ------------------------------------------------------------
ggarrange(ProcrustResAver(res_its_lsu, paletteCB4, labels_its_lsu_mark, "Res_Err", "Market") + 
            labs(title = "Market", x=NULL, y="Procrustes residual") +
            annotate("text", Inf, -Inf, label = "italic(p) == 0.012", parse=TRUE, size = 2.15, hjust = 1, vjust = -0.1),
          ProcrustResAver(res_its_lsu, paletteCB2, my_labels = "", "Res_Err", "Type") + 
            labs(title = "Niche", x=NULL, y=NULL),
          ProcrustResAver(res_its_16s, paletteCB4, labels_its_16s_mark, "Res_Err", "Market") + 
            annotate("text", Inf, -Inf, label = "italic(p) == 0.00019", parse=TRUE, size = 2.15, hjust = 1, vjust = -0.1) +
            labs(title = "Market", x=NULL, y="Procrustes residual"),
          ProcrustResAver(res_its_16s, paletteCB2, labels_its_16s_niche, "Res_Err", "Type") +
            annotate("text", Inf, -Inf, label = "italic(p) == 0.00034", parse=TRUE, size = 2.15, hjust = 1, vjust = -0.1) +
            labs(title = "Niche", x=NULL, y=NULL),
          ProcrustResAver(res_16s_lsu, paletteCB4, labels_16s_lsu_mark, "Res_Err", "Market") + 
            annotate("text", Inf, -Inf, label = "italic(p) == 0.00016", parse=TRUE, size = 2.15, hjust = 1, vjust = -0.1) +
            labs(title = "Market", x=NULL, y="Procrustes residual"),
          ProcrustResAver(res_16s_lsu, paletteCB2, my_labels = "", "Res_Err", "Type") + 
            labs(title = "Niche", x=NULL, y=NULL),
          widths = c(1,1,1,1,1,1),
          labels = c("D","","E","","F",""),
          align = "hv",
          ncol = 6, 
          nrow = 1) -> prcrust_err

prcrust_err

# *** FIGURE 4 - complete -------------------------------------------------------------------
ggarrange(Fig_4_ord, 
          prcrust_err,
          heights = c(1, 0.8),
          align = "hv",
          ncol = 1, 
          nrow = 2, 
          legend = c("bottom")) -> Fig4_procr

Fig4_procr


# **************************************-------------------------------------------------------
# RANDOM FOREST -------------------------------------------------------------------------------
library(randomForest)
library(scales)

ReformatTaxonomy(physeq_ITS_scaled) -> physeq_ITS_scaled_filt
tax_table(physeq_ITS_scaled_filt)

otu_ITS <- as.data.frame(otu_table(physeq_ITS_scaled_filt))
metadata_ITS = as(sample_data(physeq_ITS_scaled_filt), "data.frame")
taxa_ITS <- as.data.frame(as.matrix(tax_table(physeq_ITS_scaled_filt)))
head(taxa_ITS)

# ITS - market --------------------------------------------------------------------------------
# adding the classification variable and sample names -----------------------------------------
otu_ITS_rf <- data.frame(t(otu_ITS))
otu_ITS_rf$Sample <- rownames(otu_ITS_rf)
otu_ITS_rf$Market <- metadata_ITS[rownames(otu_ITS_rf), "Market"]
head(otu_ITS_rf)

# RUNNING THE MODEL ---------------------------------------------------------------------------
# try tuning the model first
set.seed(1)
sqrt(nrow(otu_ITS))

bestmtry_market <- tuneRF(x = otu_ITS_rf[,1:(ncol(otu_ITS_rf)-2)], 
                          y = otu_ITS_rf$Market,
                          stepFactor=1.5,
                          mtryStart = 13, 
                          improve=0.0001,
                          ntree=501, 
                          nodesize=1,
                          doBest = TRUE,
                          importance=TRUE)

bestmtry_market

set.seed(2)
RF_otu_ITS_Mark_501 <- randomForest(x=otu_ITS_rf[,1:(ncol(otu_ITS_rf)-2)],
                                    y=otu_ITS_rf$Market,
                                    ntree=501,
                                    mtry = 9,
                                    importance=TRUE,
                                    proximity=TRUE)

RF_otu_ITS_Mark_501
plot(RF_otu_ITS_Mark_501)

# removing unim-ortant variables using MeanDecreaseGini ---------------------------------------
head(importance(RF_otu_ITS_Mark_501))

##PlotMean Decrease Gini and Mean Decrease Accuracy 
varImpPlot(RF_otu_ITS_Mark_501) 

# ASSESSING MODEL FIT -------------------------------------------------------------------------
# Now check to see if the random forest is actually big enough... Up to a point, the more 
# trees in the forest, the better. You can tell when you've made enough when the OOB
# no longer improves.

oob_error_ITS_Mark <- data.frame(
  Trees=rep(1:nrow(RF_otu_ITS_Mark_501$err.rate), times=5),
  Market=rep(c("Ciba","Longjin","Wanyao","New District","OOB"), 
             each=nrow(RF_otu_ITS_Mark_501$err.rate)),
  Error=c(RF_otu_ITS_Mark_501$err.rate[,"OOB"], 
          RF_otu_ITS_Mark_501$err.rate[,"Ciba"], 
          RF_otu_ITS_Mark_501$err.rate[,"Longjin"],
          RF_otu_ITS_Mark_501$err.rate[,"Wanyao"], 
          RF_otu_ITS_Mark_501$err.rate[,"New District"]))

oob_error_ITS_Mark$Market <- factor(oob_error_ITS_Mark$Market,
                                   levels = c("Ciba","Longjin","Wanyao","New District", "OOB"))

p_oob_error_ITS_Mark = ggplot(data=oob_error_ITS_Mark, aes(x=Trees, y=Error)) +
  ggtitle("Model Errors") +
  theme_classic() +
  annotate("text", x=250, y=0.9, label="OOB estimate: 4.17%", size=3) +
  geom_line(aes(color=Market)) + 
  scale_color_manual(values = c(paletteCB4, "black")) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size =8)) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="bottom")

p_oob_error_ITS_Mark

# ASSESS MODEL PERFORMANCE WITH PERMUTATIONS --------------------------------------------------
# This is slow!!
library(rfUtilities)

set.seed(3)
perm_RF_otu_ITS_Mark_501 <- rf.significance(x=RF_otu_ITS_Mark_501, xdata=otu_ITS_rf[,1:(ncol(otu_ITS_rf)-2)],
                                            nperm=999,
                                            nmtry=9, 
                                            ntree=501)  
perm_RF_otu_ITS_Mark_501

# PLOTTING THE RESULTS ------------------------------------------------------------------------
# Now let's create an MDS-plot to show how the samples are related to each 
# other. Start by converting the proximity matrix into a distance matrix.
dist_otu_ITS_Mark <- as.dist(1-RF_otu_ITS_Mark_501$proximity)
mds_otu_ITS_Mark <- cmdscale(dist_otu_ITS_Mark, eig=TRUE, x.ret=TRUE)

## calculate the percentage of variation that each MDS axis accounts for...
var_mds_otu_ITS_Mark <- round(mds_otu_ITS_Mark$eig/sum(mds_otu_ITS_Mark$eig)*100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
values_mds_otu_ITS_Mark <- mds_otu_ITS_Mark$points
mds_data_ITS_Mark <- data.frame(Sample=rownames(values_mds_otu_ITS_Mark),
                                X=values_mds_otu_ITS_Mark[,1],
                                Y=values_mds_otu_ITS_Mark[,2],
                                Market=otu_ITS_rf$Market)

# Annotation position and hjust, vjust use to modify it
# xpos ypos        annotateText        hjust    vjust
# 1 -Inf -Inf Bottom Left (h0,v0)        0        0
# 2 -Inf  Inf    Top Left (h0,v1)        0        1
# 3  Inf -Inf  Bottom Right h1,v0        1        0
# 4  Inf  Inf     Top Right h1,v1        1        1

p_nmds_data_ITS_Mark = ggplot(data=mds_data_ITS_Mark, aes(x=X, y=Y, label=Sample, color=Market)) + 
  geom_point() +
  #geom_text_repel(label=rownames(mds_data), size = 2, min.segment.length = 0.5) +
  theme_classic() +
  xlab(paste("MDS1 - ", var_mds_otu_ITS_Mark[1], "%", sep="")) +
  ylab(paste("MDS2 - ", var_mds_otu_ITS_Mark[2], "%", sep="")) +
  ggtitle("ITS Market") +
  scale_color_manual(values = paletteCB4) +
  annotate("text", Inf, -Inf,  hjust = 1, vjust = - 0.8, 
           label="OOB = 4.17%", size=2.5) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) + # to adjust decimals
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  #guides(color=guide_legend(ncol=1)) +
  theme(legend.position="bottom")

p_nmds_data_ITS_Mark

# Identifying Important Features --------------------------------------------------------------
imp_RF_otu_ITS_Mark <- as.data.frame(RF_otu_ITS_Mark_501$importance)
imp_RF_otu_ITS_Mark$features <- rownames(imp_RF_otu_ITS_Mark)
imp_RF_otu_ITS_Mark <- arrange(imp_RF_otu_ITS_Mark, desc(MeanDecreaseAccuracy))
rownames(imp_RF_otu_ITS_Mark) <- imp_RF_otu_ITS_Mark$features
head(imp_RF_otu_ITS_Mark)

taxa_ITS[rownames(taxa_ITS)%in%rownames(imp_RF_otu_ITS_Mark), ] -> taxa_otu_ITS_Mark
head(taxa_otu_ITS_Mark)

# # adding marker gene info
# taxa_otu_ITS_Mark$Gene <- rownames(taxa_otu_ITS_Mark)
# taxa_otu_ITS_Mark$Gene <- gsub("POTU.*","16S",taxa_otu_ITS_Mark$Gene)
# taxa_otu_ITS_Mark$Gene <- gsub("FOTU.*","ITS",taxa_otu_ITS_Mark$Gene)
# taxa_otu_ITS_Mark$Gene <- gsub("LOTU.*","LSU",taxa_otu_ITS_Mark$Gene)
# taxa_otu_ITS_Mark$Gene <- gsub("AOTU.*","LSU",taxa_otu_ITS_Mark$Gene)

identical(rownames(taxa_otu_ITS_Mark), rownames(imp_RF_otu_ITS_Mark))
order_taxa_otu_ITS_Mark <- match(rownames(taxa_otu_ITS_Mark), rownames(imp_RF_otu_ITS_Mark))
imp_RF_otu_ITS_Mark <- imp_RF_otu_ITS_Mark[order_taxa_otu_ITS_Mark,]

head(imp_RF_otu_ITS_Mark)

imp_RF_otu_ITS_Mark$Taxonomy <- taxa_otu_ITS_Mark$Taxonomy
imp_RF_otu_ITS_Mark <- imp_RF_otu_ITS_Mark[order(imp_RF_otu_ITS_Mark$MeanDecreaseAccuracy, decreasing = TRUE),] 
head(imp_RF_otu_ITS_Mark)


# LSU - market --------------------------------------------------------------------------------
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
  taxa_table <- taxa_table[c(9,1,2,3,4,5,6,7,8)] 
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

ReformatTaxonomy(physeq_LSU_scaled) -> physeq_LSU_scaled_filt
head(tax_table(physeq_LSU_scaled_filt))

otu_LSU <- as.data.frame(otu_table(physeq_LSU_scaled_filt))
metadata_LSU = as(sample_data(physeq_LSU_scaled_filt), "data.frame")
taxa_LSU <- as.data.frame(as.matrix(tax_table(physeq_LSU_scaled_filt)))

ReformatTaxonomy(physeq_Ani_scaled) -> physeq_Ani_scaled_filt
head(tax_table(physeq_Ani_scaled_filt))

# adding the classification variable and sample names -----------------------------------------
otu_LSU_rf <- data.frame(t(otu_LSU))
otu_LSU_rf$Sample <- rownames(otu_LSU_rf)
otu_LSU_rf$Market <- metadata_LSU[rownames(otu_LSU_rf), "Market"]
head(otu_LSU_rf)

# RUNNING THE MODEL ---------------------------------------------------------------------------
# try tuning the model first
set.seed(8)
sqrt(nrow(otu_LSU))

bestmtry_market <- tuneRF(x = otu_LSU_rf[,1:(ncol(otu_LSU_rf)-2)], 
                          y = otu_LSU_rf$Market,
                          stepFactor=1.5,
                          mtryStart = 13, 
                          improve=0.0001,
                          ntree=501, 
                          nodesize=1,
                          doBest = TRUE,
                          importance=TRUE)

bestmtry_market

set.seed(9)
RF_otu_LSU_Mark_501 <- randomForest(x=otu_LSU_rf[,1:(ncol(otu_LSU_rf)-2)],
                                    y=otu_LSU_rf$Market,
                                    ntree=501,
                                    mtry = 13,
                                    nodesize=1,
                                    importance=TRUE,
                                    proximity=TRUE)

RF_otu_LSU_Mark_501
plot(RF_otu_LSU_Mark_501)

# removing unim-ortant variables using MeanDecreaseGini ---------------------------------------
head(importance(RF_otu_LSU_Mark_501))

##PlotMean Decrease Gini and Mean Decrease Accuracy 
varImpPlot(RF_otu_LSU_Mark_501) 

# ASSESSING MODEL FIT -------------------------------------------------------------------------
# Now check to see if the random forest is actually big enough... Up to a point, the more 
# trees in the forest, the better. You can tell when you've made enough when the OOB
# no longer improves.

oob_error_LSU_Mark <- data.frame(
  Trees=rep(1:nrow(RF_otu_LSU_Mark_501$err.rate), times=5),
  Market=rep(c("OOB","Ciba","Longjin","Wanyao","New District"), 
             each=nrow(RF_otu_LSU_Mark_501$err.rate)),
  Error=c(RF_otu_LSU_Mark_501$err.rate[,"OOB"], 
          RF_otu_LSU_Mark_501$err.rate[,"Ciba"], 
          RF_otu_LSU_Mark_501$err.rate[,"Longjin"],
          RF_otu_LSU_Mark_501$err.rate[,"Wanyao"], 
          RF_otu_LSU_Mark_501$err.rate[,"New District"]))

oob_error_LSU_Mark$Market <- factor(oob_error_LSU_Mark$Market,
                                   levels = c("Ciba","Longjin","Wanyao","New District", "OOB"))

p_oob_error_LSU_Mark = ggplot(data=oob_error_LSU_Mark, aes(x=Trees, y=Error)) +
  ggtitle("Model Errors") +
  theme_classic() +
  annotate("text", x=250, y=0.9, label="OOB estimate: 1.39%", size=3) +
  geom_line(aes(color=Market)) + 
  scale_color_manual(values = c(paletteCB4, "black")) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size =8)) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="bottom")

p_oob_error_LSU_Mark

# ASSESS MODEL PERFORMANCE WITH PERMUTATIONS --------------------------------------------------
set.seed(10)
perm_RF_otu_LSU_Mark_501 <- rf.significance(x=RF_otu_LSU_Mark_501, xdata=otu_LSU_rf[,1:(ncol(otu_LSU_rf)-2)],
                                            nperm=999,
                                            nmtry=9, 
                                            ntree=501)  
perm_RF_otu_LSU_Mark_501

# PLOTTING THE RESULTS ------------------------------------------------------------------------
# Now let's create an MDS-plot to show how the samples are related to each 
# other. Start by converting the proximity matrix into a distance matrix.
dist_otu_LSU_Mark <- as.dist(1-RF_otu_LSU_Mark_501$proximity)
mds_otu_LSU_Mark <- cmdscale(dist_otu_LSU_Mark, eig=TRUE, x.ret=TRUE)

## calculate the percentage of variation that each MDS axis accounts for...
var_mds_otu_LSU_Mark <- round(mds_otu_LSU_Mark$eig/sum(mds_otu_LSU_Mark$eig)*100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
values_mds_otu_LSU_Mark <- mds_otu_LSU_Mark$points
mds_data_LSU_Mark <- data.frame(Sample=rownames(values_mds_otu_LSU_Mark),
                                X=values_mds_otu_LSU_Mark[,1],
                                Y=values_mds_otu_LSU_Mark[,2],
                                Market=otu_LSU_rf$Market)

p_nmds_data_LSU_Mark = ggplot(data=mds_data_LSU_Mark, aes(x=X, y=Y, label=Sample, color=Market)) + 
  geom_point() +
  #geom_text_repel(label=rownames(mds_data), size = 2, min.segment.length = 0.5) +
  theme_classic() +
  xlab(paste("MDS1 - ", var_mds_otu_LSU_Mark[1], "%", sep="")) +
  ylab(paste("MDS2 - ", var_mds_otu_LSU_Mark[2], "%", sep="")) +
  ggtitle("LSU Market") +
  scale_color_manual(values = paletteCB4) +
  annotate("text", -Inf, -Inf,  hjust = - 0.1, vjust = - 0.8, 
           label="OOB = 12.5%", size=2.5) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) + # to adjust decimals
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  #guides(color=guide_legend(ncol=1)) +
  theme(legend.position="bottom")

p_nmds_data_LSU_Mark

# Identifying Important Features --------------------------------------------------------------
imp_RF_otu_LSU_Mark <- as.data.frame(RF_otu_LSU_Mark_501$importance)
imp_RF_otu_LSU_Mark$features <- rownames(imp_RF_otu_LSU_Mark)
imp_RF_otu_LSU_Mark <- arrange(imp_RF_otu_LSU_Mark, desc(MeanDecreaseAccuracy))
rownames(imp_RF_otu_LSU_Mark) <- imp_RF_otu_LSU_Mark$features
head(imp_RF_otu_LSU_Mark)

taxa_LSU[rownames(taxa_LSU)%in%rownames(imp_RF_otu_LSU_Mark), ] -> taxa_otu_LSU_Mark
head(taxa_otu_LSU_Mark)

# # adding marker gene info
# taxa_otu_LSU_Mark$Gene <- rownames(taxa_otu_LSU_Mark)
# taxa_otu_LSU_Mark$Gene <- gsub("POTU.*","16S",taxa_otu_LSU_Mark$Gene)
# taxa_otu_LSU_Mark$Gene <- gsub("FOTU.*","ITS",taxa_otu_LSU_Mark$Gene)
# taxa_otu_LSU_Mark$Gene <- gsub("LOTU.*","LSU",taxa_otu_LSU_Mark$Gene)
# taxa_otu_LSU_Mark$Gene <- gsub("AOTU.*","LSU",taxa_otu_LSU_Mark$Gene)

identical(rownames(taxa_otu_LSU_Mark), rownames(imp_RF_otu_LSU_Mark))
order_taxa_otu_LSU_Mark <- match(rownames(taxa_otu_LSU_Mark), rownames(imp_RF_otu_LSU_Mark))
imp_RF_otu_LSU_Mark <- imp_RF_otu_LSU_Mark[order_taxa_otu_LSU_Mark,]

head(imp_RF_otu_LSU_Mark)

imp_RF_otu_LSU_Mark$Taxonomy <- taxa_otu_LSU_Mark$Taxonomy
imp_RF_otu_LSU_Mark <- imp_RF_otu_LSU_Mark[order(imp_RF_otu_LSU_Mark$MeanDecreaseAccuracy, decreasing = TRUE),] 
head(imp_RF_otu_LSU_Mark)


# 16S - market --------------------------------------------------------------------------------
ReformatTaxonomy(physeq_16s_scaled) -> physeq_16s_scaled_filt
head(tax_table(physeq_16s_scaled_filt))

otu_16s <- as.data.frame(otu_table(physeq_16s_scaled_filt))
metadata_16s = as(sample_data(physeq_16s_scaled_filt), "data.frame")
taxa_16s <- as.data.frame(as.matrix(tax_table(physeq_16s_scaled_filt)))

# adding the classification variable and sample names -----------------------------------------
otu_16s_rf <- data.frame(t(otu_16s))
otu_16s_rf$Sample <- rownames(otu_16s_rf)
otu_16s_rf$Market <- metadata_16s[rownames(otu_16s_rf), "Market"]
head(otu_16s_rf)

# RUNNING THE MODEL ---------------------------------------------------------------------------
# try tuning the model first
set.seed(4)
sqrt(nrow(otu_16s))

bestmtry_market <- tuneRF(x = otu_16s_rf[,1:(ncol(otu_16s_rf)-2)], 
                          y = otu_16s_rf$Market,
                          stepFactor=1.5,
                          mtryStart = 19, 
                          improve=0.0001,
                          ntree=501, 
                          nodesize=1,
                          doBest = TRUE,
                          importance=TRUE)

bestmtry_market

set.seed(5)
RF_otu_16s_Mark_501 <- randomForest(x=otu_16s_rf[,1:(ncol(otu_16s_rf)-2)],
                                    y=otu_16s_rf$Market,
                                    ntree=501,
                                    mtry = 19,
                                    nodesize=1,
                                    importance=TRUE,
                                    proximity=TRUE)

RF_otu_16s_Mark_501
plot(RF_otu_16s_Mark_501)

# removing unim-ortant variables using MeanDecreaseGini ---------------------------------------
head(importance(RF_otu_16s_Mark_501))

##PlotMean Decrease Gini and Mean Decrease Accuracy 
varImpPlot(RF_otu_16s_Mark_501) 

# ASSESSING MODEL FIT -------------------------------------------------------------------------
# Now check to see if the random forest is actually big enough... Up to a point, the more 
# trees in the forest, the better. You can tell when you've made enough when the OOB
# no longer improves.

oob_error_16s_Mark <- data.frame(
  Trees=rep(1:nrow(RF_otu_16s_Mark_501$err.rate), times=5),
  Market=rep(c("OOB","Ciba","Longjin","Wanyao","New District"), 
             each=nrow(RF_otu_16s_Mark_501$err.rate)),
  Error=c(RF_otu_16s_Mark_501$err.rate[,"OOB"], 
          RF_otu_16s_Mark_501$err.rate[,"Ciba"], 
          RF_otu_16s_Mark_501$err.rate[,"Longjin"],
          RF_otu_16s_Mark_501$err.rate[,"Wanyao"], 
          RF_otu_16s_Mark_501$err.rate[,"New District"]))

oob_error_16s_Mark$Market <- factor(oob_error_16s_Mark$Market,
                                    levels = c("Ciba","Longjin","Wanyao","New District", "OOB"))

p_oob_error_16s_Mark = ggplot(data=oob_error_16s_Mark, aes(x=Trees, y=Error)) +
  ggtitle("Model Errors") +
  theme_classic() +
  annotate("text", x=250, y=0.9, label="OOB estimate: 1.39%", size=3) +
  geom_line(aes(color=Market)) + 
  scale_color_manual(values = c(paletteCB4, "black")) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size =8)) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="bottom")

p_oob_error_16s_Mark

# ASSESS MODEL PERFORMANCE WITH PERMUTATIONS --------------------------------------------------
set.seed(6)
perm_RF_otu_16s_Mark_501 <- rf.significance(x=RF_otu_16s_Mark_501, xdata=otu_16s_rf[,1:(ncol(otu_16s_rf)-2)],
                                            nperm=999,
                                            nmtry=9, 
                                            ntree=501)  
perm_RF_otu_16s_Mark_501

# PLOTTING THE RESULTS ------------------------------------------------------------------------
# Now let's create an MDS-plot to show how the samples are related to each 
# other. Start by converting the proximity matrix into a distance matrix.
dist_otu_16s_Mark <- as.dist(1-RF_otu_16s_Mark_501$proximity)
mds_otu_16s_Mark <- cmdscale(dist_otu_16s_Mark, eig=TRUE, x.ret=TRUE)

## calculate the percentage of variation that each MDS axis accounts for...
var_mds_otu_16s_Mark <- round(mds_otu_16s_Mark$eig/sum(mds_otu_16s_Mark$eig)*100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
values_mds_otu_16s_Mark <- mds_otu_16s_Mark$points
mds_data_16s_Mark <- data.frame(Sample=rownames(values_mds_otu_16s_Mark),
                                X=values_mds_otu_16s_Mark[,1],
                                Y=values_mds_otu_16s_Mark[,2],
                                Market=otu_16s_rf$Market)

p_nmds_data_16s_Mark = ggplot(data=mds_data_16s_Mark, aes(x=X, y=Y, label=Sample, color=Market)) + 
  geom_point() +
  #geom_text_repel(label=rownames(mds_data), size = 2, min.segment.length = 0.5) +
  theme_classic() +
  xlab(paste("MDS1 - ", var_mds_otu_16s_Mark[1], "%", sep="")) +
  ylab(paste("MDS2 - ", var_mds_otu_16s_Mark[2], "%", sep="")) +
  ggtitle("16S Market") +
  scale_color_manual(values = paletteCB4) +
  annotate("text", Inf, -Inf,  hjust = 1, vjust = - 0.8, 
           label="OOB = 1.39%", size=2.5) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) + # to adjust decimals
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  #guides(color=guide_legend(ncol=1)) +
  theme(legend.position="bottom")

p_nmds_data_16s_Mark

# Identifying Important Features --------------------------------------------------------------
imp_RF_otu_16s_Mark <- as.data.frame(RF_otu_16s_Mark_501$importance)
imp_RF_otu_16s_Mark$features <- rownames(imp_RF_otu_16s_Mark)
imp_RF_otu_16s_Mark <- arrange(imp_RF_otu_16s_Mark, desc(MeanDecreaseAccuracy))
rownames(imp_RF_otu_16s_Mark) <- imp_RF_otu_16s_Mark$features
head(imp_RF_otu_16s_Mark)

taxa_16s[rownames(taxa_16s)%in%rownames(imp_RF_otu_16s_Mark), ] -> taxa_otu_16s_Mark
head(taxa_otu_16s_Mark)

# # adding marker gene info
# taxa_otu_16s_Mark$Gene <- rownames(taxa_otu_16s_Mark)
# taxa_otu_16s_Mark$Gene <- gsub("POTU.*","16S",taxa_otu_16s_Mark$Gene)
# taxa_otu_16s_Mark$Gene <- gsub("FOTU.*","ITS",taxa_otu_16s_Mark$Gene)
# taxa_otu_16s_Mark$Gene <- gsub("LOTU.*","LSU",taxa_otu_16s_Mark$Gene)
# taxa_otu_16s_Mark$Gene <- gsub("AOTU.*","LSU",taxa_otu_16s_Mark$Gene)

identical(rownames(taxa_otu_16s_Mark), rownames(imp_RF_otu_16s_Mark))
order_taxa_otu_16s_Mark <- match(rownames(taxa_otu_16s_Mark), rownames(imp_RF_otu_16s_Mark))
imp_RF_otu_16s_Mark <- imp_RF_otu_16s_Mark[order_taxa_otu_16s_Mark,]

head(imp_RF_otu_16s_Mark)

imp_RF_otu_16s_Mark$Taxonomy <- taxa_otu_16s_Mark$Taxonomy
imp_RF_otu_16s_Mark <- imp_RF_otu_16s_Mark[order(imp_RF_otu_16s_Mark$MeanDecreaseAccuracy, decreasing = TRUE),] 
head(imp_RF_otu_16s_Mark)


# ITS - niche ---------------------------------------------------------------------------------
# adding the classification variable and sample names -----------------------------------------
otu_ITS_ty <- data.frame(t(otu_ITS))
otu_ITS_ty$Sample <- rownames(otu_ITS_ty)
otu_ITS_ty$Type <- metadata_ITS[rownames(otu_ITS_ty), "Type"] 
head(otu_ITS_ty)

# RUNNING THE MODEL ---------------------------------------------------------------------------
set.seed(10)
bestmtry_market_ty <- tuneRF(x = otu_ITS_ty[,1:(ncol(otu_ITS_ty)-2)], 
                             y = otu_ITS_ty$Type,
                             stepFactor=1.5,
                             mtryStart = 13, 
                             improve=0.0001,
                             ntree=501, 
                             nodesize=1,
                             doBest = TRUE,
                             importance=TRUE)

bestmtry_market_ty

set.seed(12)
RF_otu_ITS_Type_501 <- randomForest(x=otu_ITS_ty[,1:(ncol(otu_ITS_ty)-2)],
                                    y=otu_ITS_ty$Type,
                                    ntree=501,
                                    mtry = 19,
                                    importance=TRUE,
                                    proximity=TRUE)

RF_otu_ITS_Type_501
plot(RF_otu_ITS_Type_501)

# Now check to see if the random forest is actually big enough... Up to a point, the more 
# trees in the forest, the better. You can tell when you've made enough when the OOB
# no longer improves.

oob_error_ITS_Type <- data.frame(
  Trees=rep(1:nrow(RF_otu_ITS_Type_501$err.rate), times=3),
  Niche=rep(c("OOB","in", "out"), 
            each=nrow(RF_otu_ITS_Type_501$err.rate)),
  Error=c(RF_otu_ITS_Type_501$err.rate[,"OOB"], 
          RF_otu_ITS_Type_501$err.rate[,"in"], 
          RF_otu_ITS_Type_501$err.rate[,"out"]))

oob_error_ITS_Type$Niche <- factor(oob_error_ITS_Type$Niche,
                                   levels = c("in","out", "OOB"))

oob_error_ITS_Type$Niche <- oob_error_ITS_Type$Niche %>% 
  recode_factor('in' = "Inside", 'out'="Outside", "OOB"="OOB")

p_oob_error_ITS_Type = ggplot(data=oob_error_ITS_Type, aes(x=Trees, y=Error)) +
  ggtitle("Model Errors") +
  theme_classic() +
  annotate("text", x=250, y=0.9, label="OOB estimate: 8.33%", size=3) +
  geom_line(aes(color=Niche)) + 
  scale_color_manual(values = c(paletteCB2, "red")) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size =8)) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="bottom")

p_oob_error_ITS_Type

# ASSESS MODEL PERFORMANCE WITH PERMUTATIONS --------------------------------------------------
perm_RF_otu_ITS_Type_501 <- rf.significance(x=RF_otu_ITS_Type_501, xdata=otu_ITS_ty[,1:(ncol(otu_ITS_ty)-2)],
                                            nperm=999,
                                            mtry = 19,
                                            ntree=501)
perm_RF_otu_ITS_Type_501

# PLOTTING THE RESULTS ------------------------------------------------------------------------
# Now let's create an MDS-plot to show how the samples are related to each 
# other. Start by converting the proximity matrix into a distance matrix.
dist_otu_ITS_Type <- as.dist(1-RF_otu_ITS_Type_501$proximity)
mds_otu_ITS_Type <- cmdscale(dist_otu_ITS_Type, eig=TRUE, x.ret=TRUE)

## calculate the percentage of variation that each MDS axis accounts for...
var_mds_otu_ITS_Type <- round(mds_otu_ITS_Type$eig/sum(mds_otu_ITS_Type$eig)*100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
values_mds_otu_ITS_Type <- mds_otu_ITS_Type$points
mds_data_ITS_Type <- data.frame(Sample=rownames(values_mds_otu_ITS_Type),
                                X=values_mds_otu_ITS_Type[,1],
                                Y=values_mds_otu_ITS_Type[,2],
                                Niche=otu_ITS_ty$Type)

mds_data_ITS_Type$Niche <- mds_data_ITS_Type$Niche %>% 
        recode_factor('in' = "Inside", 'out'="Outside")

p_nmds_data_ITS_Type = ggplot(data=mds_data_ITS_Type, aes(x=X, y=Y, label=Sample, color=Niche)) + 
  geom_point() +
  #geom_text_repel(label=rownames(mds_data), size = 2, min.segment.length = 0.5) +
  theme_classic() +
  xlab(paste("MDS1 - ", var_mds_otu_ITS_Type[1], "%", sep="")) +
  ylab(paste("MDS2 - ", var_mds_otu_ITS_Type[2], "%", sep="")) +
  ggtitle("ITS Niche") +
  scale_color_manual("Niche", values = paletteCB2) +
  annotate("text", -Inf, -Inf,  hjust = - 0.1, vjust = - 0.8,
           label="OOB = 8.33%", size=2.5) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) + # to adjust decimals
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  #guides(color=guide_legend(ncol=2)) +
  theme(legend.position="right")

p_nmds_data_ITS_Type

# Identifying Important Features --------------------------------------------------------------
imp_RF_otu_ITS_Type <- as.data.frame(RF_otu_ITS_Type_501$importance)
imp_RF_otu_ITS_Type$features <- rownames(imp_RF_otu_ITS_Type)
imp_RF_otu_ITS_Type <- arrange(imp_RF_otu_ITS_Type, desc(MeanDecreaseAccuracy))
rownames(imp_RF_otu_ITS_Type) <- imp_RF_otu_ITS_Type$features
head(imp_RF_otu_ITS_Type)

taxa_ITS[rownames(taxa_ITS)%in%rownames(imp_RF_otu_ITS_Type), ] -> taxa_otu_ITS_Type
head(taxa_otu_ITS_Type)

# # adding typeer gene info
# taxa_otu_ITS_Type$Gene <- rownames(taxa_otu_ITS_Type)
# taxa_otu_ITS_Type$Gene <- gsub("POTU.*","16S",taxa_otu_ITS_Type$Gene)
# taxa_otu_ITS_Type$Gene <- gsub("FOTU.*","ITS",taxa_otu_ITS_Type$Gene)
# taxa_otu_ITS_Type$Gene <- gsub("LOTU.*","LSU",taxa_otu_ITS_Type$Gene)
# taxa_otu_ITS_Type$Gene <- gsub("AOTU.*","LSU",taxa_otu_ITS_Type$Gene)

identical(rownames(taxa_otu_ITS_Type), rownames(imp_RF_otu_ITS_Type))
order_taxa_otu_ITS_Type <- match(rownames(taxa_otu_ITS_Type), rownames(imp_RF_otu_ITS_Type))
imp_RF_otu_ITS_Type <- imp_RF_otu_ITS_Type[order_taxa_otu_ITS_Type,]

head(imp_RF_otu_ITS_Type)

imp_RF_otu_ITS_Type$Taxonomy <- taxa_otu_ITS_Type$Taxonomy
imp_RF_otu_ITS_Type <- imp_RF_otu_ITS_Type[order(imp_RF_otu_ITS_Type$MeanDecreaseAccuracy, decreasing = TRUE),] 
head(imp_RF_otu_ITS_Type)


# LSU - niche ---------------------------------------------------------------------------------
# adding the classification variable and sample names -----------------------------------------
otu_LSU_ty <- data.frame(t(otu_LSU))
otu_LSU_ty$Sample <- rownames(otu_LSU_ty)
otu_LSU_ty$Type <- metadata_LSU[rownames(otu_LSU_ty), "Type"] 
head(otu_LSU_ty)

# RUNNING THE MODEL ---------------------------------------------------------------------------
set.seed(15)
bestmtry_market_ty <- tuneRF(x = otu_LSU_ty[,1:(ncol(otu_LSU_ty)-2)], 
                             y = otu_LSU_ty$Type,
                             stepFactor=1.5,
                             mtryStart = 13, 
                             improve=0.0001,
                             ntree=501, 
                             nodesize=1,
                             doBest = TRUE,
                             importance=TRUE)

bestmtry_market_ty

set.seed(16)
RF_otu_LSU_Type_501 <- randomForest(x=otu_LSU_ty[,1:(ncol(otu_LSU_ty)-2)],
                                    y=otu_LSU_ty$Type,
                                    ntree=501,
                                    mtry = 13,
                                    importance=TRUE,
                                    proximity=TRUE)

RF_otu_LSU_Type_501
plot(RF_otu_LSU_Type_501)

# Now check to see if the random forest is actually big enough... Up to a point, the more 
# trees in the forest, the better. You can tell when you've made enough when the OOB
# no longer improves.

oob_error_LSU_Type <- data.frame(
  Trees=rep(1:nrow(RF_otu_LSU_Type_501$err.rate), times=3),
  Niche=rep(c("OOB","in", "out"), 
            each=nrow(RF_otu_LSU_Type_501$err.rate)),
  Error=c(RF_otu_LSU_Type_501$err.rate[,"OOB"], 
          RF_otu_LSU_Type_501$err.rate[,"in"], 
          RF_otu_LSU_Type_501$err.rate[,"out"]))

oob_error_LSU_Type$Niche <- factor(oob_error_LSU_Type$Niche,
                                   levels = c("in","out", "OOB"))

oob_error_LSU_Type$Niche <- oob_error_LSU_Type$Niche %>% 
  recode_factor('in' = "Inside", 'out'="Outside", "OOB"="OOB")

p_oob_error_LSU_Type = ggplot(data=oob_error_LSU_Type, aes(x=Trees, y=Error)) +
  ggtitle("Model Errors") +
  theme_classic() +
  annotate("text", x=250, y=0.9, label="OOB estimate: 12.5%", size=3) +
  geom_line(aes(color=Niche)) + 
  scale_color_manual(values = c(paletteCB2, "red")) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size =8)) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="bottom")

p_oob_error_LSU_Type

# ASSESS MODEL PERFORMANCE WITH PERMUTATIONS --------------------------------------------------
perm_RF_otu_LSU_Type_501 <- rf.significance(x=RF_otu_LSU_Type_501, xdata=otu_LSU_ty[,1:(ncol(otu_LSU_ty)-2)],
                                            nperm=999,
                                            mtry = 13,
                                            ntree=501)
perm_RF_otu_LSU_Type_501

# PLOTTING THE RESULTS ------------------------------------------------------------------------
# Now let's create an MDS-plot to show how the samples are related to each 
# other. Start by converting the proximity matrix into a distance matrix.
dist_otu_LSU_Type <- as.dist(1-RF_otu_LSU_Type_501$proximity)
mds_otu_LSU_Type <- cmdscale(dist_otu_LSU_Type, eig=TRUE, x.ret=TRUE)

## calculate the percentage of variation that each MDS axis accounts for...
var_mds_otu_LSU_Type <- round(mds_otu_LSU_Type$eig/sum(mds_otu_LSU_Type$eig)*100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
values_mds_otu_LSU_Type <- mds_otu_LSU_Type$points
mds_data_LSU_Type <- data.frame(Sample=rownames(values_mds_otu_LSU_Type),
                                X=values_mds_otu_LSU_Type[,1],
                                Y=values_mds_otu_LSU_Type[,2],
                                Niche=otu_LSU_ty$Type)

mds_data_LSU_Type$Niche <- mds_data_LSU_Type$Niche %>% 
  recode_factor('in' = "Inside", 'out'="Outside")

p_nmds_data_LSU_Type = ggplot(data=mds_data_LSU_Type, aes(x=X, y=Y, label=Sample, color=Niche)) + 
  geom_point() +
  #geom_text_repel(label=rownames(mds_data), size = 2, min.segment.length = 0.5) +
  theme_classic() +
  xlab(paste("MDS1 - ", var_mds_otu_LSU_Type[1], "%", sep="")) +
  ylab(paste("MDS2 - ", var_mds_otu_LSU_Type[2], "%", sep="")) +
  ggtitle("LSU Niche") +
  scale_color_manual("Niche", values = paletteCB2) +
  annotate("text", -Inf, -Inf,  hjust = - 0.1, vjust = - 0.8,
           label="OOB = 12.5%", size=2.5) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) + # to adjust decimals
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  #guides(color=guide_legend(ncol=2)) +
  theme(legend.position="right")

p_nmds_data_LSU_Type

# Identifying Important Features --------------------------------------------------------------
imp_RF_otu_LSU_Type <- as.data.frame(RF_otu_LSU_Type_501$importance)
imp_RF_otu_LSU_Type$features <- rownames(imp_RF_otu_LSU_Type)
imp_RF_otu_LSU_Type <- arrange(imp_RF_otu_LSU_Type, desc(MeanDecreaseAccuracy))
rownames(imp_RF_otu_LSU_Type) <- imp_RF_otu_LSU_Type$features
head(imp_RF_otu_LSU_Type)

taxa_LSU[rownames(taxa_LSU)%in%rownames(imp_RF_otu_LSU_Type), ] -> taxa_otu_LSU_Type
head(taxa_otu_LSU_Type)

# # adding typeer gene info
# taxa_otu_LSU_Type$Gene <- rownames(taxa_otu_LSU_Type)
# taxa_otu_LSU_Type$Gene <- gsub("POTU.*","16S",taxa_otu_LSU_Type$Gene)
# taxa_otu_LSU_Type$Gene <- gsub("FOTU.*","LSU",taxa_otu_LSU_Type$Gene)
# taxa_otu_LSU_Type$Gene <- gsub("LOTU.*","LSU",taxa_otu_LSU_Type$Gene)
# taxa_otu_LSU_Type$Gene <- gsub("AOTU.*","LSU",taxa_otu_LSU_Type$Gene)

identical(rownames(taxa_otu_LSU_Type), rownames(imp_RF_otu_LSU_Type))
order_taxa_otu_LSU_Type <- match(rownames(taxa_otu_LSU_Type), rownames(imp_RF_otu_LSU_Type))
imp_RF_otu_LSU_Type <- imp_RF_otu_LSU_Type[order_taxa_otu_LSU_Type,]

head(imp_RF_otu_LSU_Type)

imp_RF_otu_LSU_Type$Taxonomy <- taxa_otu_LSU_Type$Taxonomy
imp_RF_otu_LSU_Type <- imp_RF_otu_LSU_Type[order(imp_RF_otu_LSU_Type$MeanDecreaseAccuracy, decreasing = TRUE),] 
head(imp_RF_otu_LSU_Type)


# 16S - niche ---------------------------------------------------------------------------------
# adding the classification variable and sample names -----------------------------------------
otu_16s_ty <- data.frame(t(otu_16s))
otu_16s_ty$Sample <- rownames(otu_16s_ty)
otu_16s_ty$Type <- metadata_16s[rownames(otu_16s_ty), "Type"] 
head(otu_16s_ty)

# RUNNING THE MODEL ---------------------------------------------------------------------------
sqrt(nrow(otu_16s))
set.seed(21)
bestmtry_market_ty <- tuneRF(x = otu_16s_ty[,1:(ncol(otu_16s_ty)-2)], 
                             y = otu_16s_ty$Type,
                             stepFactor=1.5,
                             mtryStart = 19, 
                             improve=0.0001,
                             ntree=501, 
                             nodesize=1,
                             doBest = TRUE,
                             importance=TRUE)

bestmtry_market_ty

set.seed(22)
RF_otu_16s_Type_501 <- randomForest(x=otu_16s_ty[,1:(ncol(otu_16s_ty)-2)],
                                    y=otu_16s_ty$Type,
                                    ntree=501,
                                    mtryStart = 19,
                                    importance=TRUE,
                                    proximity=TRUE)

RF_otu_16s_Type_501
plot(RF_otu_16s_Type_501)

# Now check to see if the random forest is actually big enough... Up to a point, the more 
# trees in the forest, the better. You can tell when you've made enough when the OOB
# no longer improves.

oob_error_16s_Type <- data.frame(
  Trees=rep(1:nrow(RF_otu_16s_Type_501$err.rate), times=3),
  Niche=rep(c("OOB","in", "out"), 
            each=nrow(RF_otu_16s_Type_501$err.rate)),
  Error=c(RF_otu_16s_Type_501$err.rate[,"OOB"], 
          RF_otu_16s_Type_501$err.rate[,"in"], 
          RF_otu_16s_Type_501$err.rate[,"out"]))

oob_error_16s_Type$Niche <- factor(oob_error_16s_Type$Niche,
                                   levels = c("in","out", "OOB"))

oob_error_16s_Type$Niche <- oob_error_16s_Type$Niche %>% 
  recode_factor('in' = "Inside", 'out'="Outside", "OOB"="OOB")

p_oob_error_16s_Type = ggplot(data=oob_error_16s_Type, aes(x=Trees, y=Error)) +
  ggtitle("Model Errors") +
  theme_classic() +
  annotate("text", x=250, y=0.9, label="OOB estimate: 9.72%", size=3) +
  geom_line(aes(color=Niche)) + 
  scale_color_manual(values = c(paletteCB2, "red")) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size =8)) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="bottom")

p_oob_error_16s_Type

# ASSESS MODEL PERFORMANCE WITH PERMUTATIONS --------------------------------------------------
perm_RF_otu_16s_Type_501 <- rf.significance(x=RF_otu_16s_Type_501, xdata=otu_16s_ty[,1:(ncol(otu_16s_ty)-2)],
                                            nperm=999,
                                            mtryStart = 19,
                                            ntree=501)  
perm_RF_otu_16s_Type_501

# PLOTTING THE RESULTS ------------------------------------------------------------------------
# Now let's create an MDS-plot to show how the samples are related to each 
# other. Start by converting the proximity matrix into a distance matrix.
dist_otu_16s_Type <- as.dist(1-RF_otu_16s_Type_501$proximity)
mds_otu_16s_Type <- cmdscale(dist_otu_16s_Type, eig=TRUE, x.ret=TRUE)

## calculate the percentage of variation that each MDS axis accounts for...
var_mds_otu_16s_Type <- round(mds_otu_16s_Type$eig/sum(mds_otu_16s_Type$eig)*100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
values_mds_otu_16s_Type <- mds_otu_16s_Type$points
mds_data_16s_Type <- data.frame(Sample=rownames(values_mds_otu_16s_Type),
                                X=values_mds_otu_16s_Type[,1],
                                Y=values_mds_otu_16s_Type[,2],
                                Niche=otu_16s_ty$Type)

mds_data_16s_Type$Niche <- mds_data_16s_Type$Niche %>% 
  recode_factor('in' = "Inside", 'out'="Outside")

p_nmds_data_16s_Type = ggplot(data=mds_data_16s_Type, aes(x=X, y=Y, label=Sample, color=Niche)) + 
  geom_point() +
  #geom_text_repel(label=rownames(mds_data), size = 2, min.segment.length = 0.5) +
  theme_classic() +
  xlab(paste("MDS1 - ", var_mds_otu_16s_Type[1], "%", sep="")) +
  ylab(paste("MDS2 - ", var_mds_otu_16s_Type[2], "%", sep="")) +
  ggtitle("16S Niche") +
  scale_color_manual("Niche", values = paletteCB2) +
  annotate("text",-Inf, -Inf,  hjust = - 0.1, vjust = - 0.8,
           label="OOB = 9.72%", size=2.5) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) + # to adjust decimals
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  #guides(color=guide_legend(ncol=2)) +
  theme(legend.position="right")

p_nmds_data_16s_Type

# Identifying Important Features --------------------------------------------------------------
imp_RF_otu_16s_Type <- as.data.frame(RF_otu_16s_Type_501$importance)
imp_RF_otu_16s_Type$features <- rownames(imp_RF_otu_16s_Type)
imp_RF_otu_16s_Type <- arrange(imp_RF_otu_16s_Type, desc(MeanDecreaseAccuracy))
rownames(imp_RF_otu_16s_Type) <- imp_RF_otu_16s_Type$features
head(imp_RF_otu_16s_Type)

taxa_16s[rownames(taxa_16s)%in%rownames(imp_RF_otu_16s_Type), ] -> taxa_otu_16s_Type
head(taxa_otu_16s_Type)

identical(rownames(taxa_otu_16s_Type), rownames(imp_RF_otu_16s_Type))
order_taxa_otu_16s_Type <- match(rownames(taxa_otu_16s_Type), rownames(imp_RF_otu_16s_Type))
imp_RF_otu_16s_Type <- imp_RF_otu_16s_Type[order_taxa_otu_16s_Type,]

head(imp_RF_otu_16s_Type)

imp_RF_otu_16s_Type$Taxonomy <- taxa_otu_16s_Type$Taxonomy
imp_RF_otu_16s_Type <- imp_RF_otu_16s_Type[order(imp_RF_otu_16s_Type$MeanDecreaseAccuracy, decreasing = TRUE),] 
head(imp_RF_otu_16s_Type)


# *** FIGURE 4 - Ranfom Forest Model ----------------------------------------------------------
# ggarrange(ggarrange(ggarrange(p_nmds_data_ITS_Mark,
#                               p_nmds_data_LSU_Mark,
#                               p_nmds_data_16s_Mark,
#                               labels = c("A","B", "C"), 
#                               align = "hv" ,
#                               widths = c(1,1,1),
#                               ncol = 3, 
#                               nrow = 1, 
#                               common.legend = TRUE,
#                               legend = c("bottom")),
#                     ggarrange(p_nmds_data_ITS_Type,
#                               p_nmds_data_LSU_Type,
#                               p_nmds_data_16s_Type,
#                               #labels = c("C","D", "E"), 
#                               align = "hv" ,
#                               widths = c(1,1,1),
#                               ncol = 3, 
#                               nrow = 1, 
#                               common.legend = TRUE,
#                               legend = c("bottom")),
#                     align = "hv" ,
#                     widths = c(1,1),
#                     ncol = 1, 
#                     nrow = 2),
#           ggarrange(p_bar_otu_ITS_Mark,
#                     p_bar_otu_LSU_Mark,
#                     p_bar_otu_16s_Mark,
#                     p_bar_otu_ITS_Type,
#                     p_bar_otu_LSU_Type,
#                     p_bar_otu_16s_Type,
#                     labels = c("C","D","E"),
#                     widths = c(1,1,1,1,1,1),
#                     align = "hv" ,
#                     ncol = 3, 
#                     nrow = 2,
#                     common.legend = TRUE,
#                     legend = c("bottom")),
#           widths = c(0.76,1),
#           ncol = 2, 
#           nrow = 1) -> Fig_4_rf
# 
# Fig_4_rf

ggarrange(
  ggarrange(
    p_nmds_data_ITS_Mark,
    p_nmds_data_LSU_Mark,
    p_nmds_data_16s_Mark,
    labels = c("A", "B", "C"),
    align = "hv" ,
    widths = c(1, 1, 1),
    ncol = 3,
    nrow = 1,
    common.legend = TRUE,
    legend = c("bottom")
  ),
  ggarrange(
    p_nmds_data_ITS_Type,
    p_nmds_data_LSU_Type,
    p_nmds_data_16s_Type,
    #labels = c("C","D", "E"),
    align = "hv" ,
    widths = c(1, 1, 1),
    ncol = 3,
    nrow = 1,
    common.legend = TRUE,
    legend = c("bottom")
  ),
  #align = "" ,
  widths = c(1, 1),
  ncol = 1,
  nrow = 2) -> Fig_4

Fig_4

# ******************************************************---------------------------------------
# INDICATOR SPECIES ANALYSIS ------------------------------------------------------------------
library(indicspecies)


GetIndicators <-function(dataframe, var){
    require(phyloseq); require(indicspecies); require(dplyr)
    otu <- as.data.frame(otu_table(dataframe))
    metadata = as(sample_data(dataframe), "data.frame")
    taxa <- as.data.frame(as.matrix(tax_table(dataframe)))
    multipatt <- multipatt(t(otu), metadata[,var], func = "r.g",
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


# indicator value >0.5 and p-value <0.05 after fdr correction
GetIndicators(physeq_ITS_scaled_filt, "Market") -> ind_ITS_mark
ind_ITS_mark
GetIndicators(physeq_ITS_scaled_filt, "Type") -> ind_ITS_type
ind_ITS_type

GetIndicators(physeq_LSU_scaled_filt, "Market") -> ind_LSU_mark
ind_LSU_mark
GetIndicators(physeq_LSU_scaled_filt, "Type") -> ind_LSU_type
ind_LSU_type

GetIndicators(physeq_16s_scaled_filt, "Market") -> ind_16s_mark
ind_16s_mark
GetIndicators(physeq_16s_scaled_filt, "Type") -> ind_16s_type
ind_16s_type


# RF AND INDICATORS MATCH ---------------------------------------------------------------------
InditoRF <- function(df_imp, df_r, var) {
  #subset(ind_ITS_mark, ind_ITS_mark$stat >= 0.5) -> ind_ITS_mark_filt
  df_imp[1:20, ] -> df_imp20
  df_imp20$index <- df_r[rownames(df_imp20),]$index
  df_imp20$stat <- df_r[rownames(df_imp20),]$stat
  df_imp20
  if (var == "Market") {
    df_imp20$index <- df_imp20$index %>%
      recode_factor(
        "1" = "Ciba",
        "2" = "Longjin",
        "3" = "Wanyao",
        "4" = "New District")
  } else {
    df_imp20$index <- df_imp20$index %>%
      recode_factor("1" = "Inside", "2" = "Outside")}
  return(df_imp20)
}

# Match indicators with RF models features accuracy -------------------------------------------
InditoRF(df_imp = imp_RF_otu_ITS_Mark, df_r = ind_ITS_mark, var = "Market") -> ind_RF_ITS_mark
ind_RF_ITS_mark
InditoRF(df_imp = imp_RF_otu_ITS_Type, df_r = ind_ITS_type, var = "Niche") -> ind_RF_ITS_niche
ind_RF_ITS_niche

InditoRF(df_imp = imp_RF_otu_LSU_Mark, df_r = ind_LSU_mark, var = "Market") -> ind_RF_LSU_mark
ind_RF_LSU_mark
InditoRF(df_imp = imp_RF_otu_LSU_Type, df_r = ind_LSU_type, var = "Niche") -> ind_RF_LSU_niche
ind_RF_LSU_niche

InditoRF(df_imp = imp_RF_otu_16s_Mark, df_r = ind_16s_mark, var = "Market") -> ind_RF_16s_mark
ind_RF_16s_mark
InditoRF(df_imp = imp_RF_otu_16s_Type, df_r = ind_16s_type, var = "Niche") -> ind_RF_16s_niche
ind_RF_16s_niche


# library(purrr)
# pal_market = list(values=paletteCB4, na.value="lightblue")
# 
# # plotting bars
# p_bar_otu_ITS_Mark = ggplot(data=ind_RF_ITS_mark, aes(x= reorder(Taxonomy, -MeanDecreaseAccuracy),
#                                                       y= MeanDecreaseAccuracy, color=index, fill=index)) + 
#   geom_bar(stat="identity") +
#   labs(title = "Top ITS OTUs to\nclassify markets", 
#        x= "OTU", y= "Mean Decrease\nin Accuracy") +
#   theme_classic() +
#   coord_flip() +
#   ylim(0, 0.06) +
#   geom_text(aes(label=round(stat, digits = 2)), position = position_stack(vjust = 0.5), color="black", size=1.8) +
#   #scale_colour_manual("Market", values = c(values=paletteCB4, na.value="grey")) +
#   invoke(scale_colour_manual, pal_market) + invoke(scale_fill_manual, pal_market) +
#   #scale_fill_manual("Market", values = c(paletteCB4, "lightblue")) +
#   theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
#   theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
#   theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
#   theme(axis.text.x = element_text(angle = 0, size = 7 ,hjust = 0.5, vjust = 0.5)) +
#   theme(axis.text.y = element_text(angle = 0, size = 6, hjust = 1, vjust = 0.5)) +
#   theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
#   guides(color=guide_legend(title="r.g. value"), fill=guide_legend(title="r.g. value")) +
#   theme(legend.position="right")
# 
# p_bar_otu_ITS_Mark


# PLOT INDICATORS match to RF -----------------------------------------------------------------
PlotIndtoRF <- function(df, palette){
  require(purrr)
  pal_market = list(values=palette, na.value="lightblue")
  ggplot(data=df, aes(x= reorder(Taxonomy,-MeanDecreaseAccuracy),
                                   y= MeanDecreaseAccuracy, color=index, fill=index)) + 
  geom_bar(stat="identity") +
  theme_classic() +
  coord_flip() +
  ylim(0, 0.065) +
  geom_text(aes(label=round(stat, digits = 2)), fontface  = "bold",
            position = position_stack(vjust = 0.5), color="white", size=1.8) +
  invoke(scale_colour_manual, pal_market) + invoke(scale_fill_manual, pal_market) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7 ,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 6, hjust = 1, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(title="r.g. value",nrow=2), 
         fill=guide_legend(title="r.g. value", nrow=2)) +
  theme(legend.position="bottom") -> plot_IND
return(plot_IND)
}



PlotIndtoRF(ind_RF_ITS_mark, paletteCB4) + labs(title = "Top ITS OTUs to\nclassify markets",  x= "OTU", y= "Mean Decrease\nin Accuracy")
PlotIndtoRF(ind_RF_ITS_niche, paletteCB2) + labs(title = "Top ITS OTUs to\nclassify niches",  x= "OTU", y= "Mean Decrease\nin Accuracy")
PlotIndtoRF(ind_RF_LSU_mark, paletteCB4) + labs(title = "Top LSU OTUs to\nclassify markets",  x= "OTU", y= "Mean Decrease\nin Accuracy")
PlotIndtoRF(ind_RF_LSU_niche, paletteCB2) + labs(title = "Top LSU OTUs to\nclassify niches",  x= "OTU", y= "Mean Decrease\nin Accuracy")
PlotIndtoRF(ind_RF_16s_mark, paletteCB4) + labs(title = "Top 16s OTUs to\nclassify markets",  x= "OTU", y= "Mean Decrease\nin Accuracy")
PlotIndtoRF(ind_RF_16s_niche, paletteCB2) + labs(title = "Top 16s OTUs to\nclassify niches",  x= "OTU", y= "Mean Decrease\nin Accuracy")
  
# *** FIGURE 5 - Ranfom Forest Model ----------------------------------------------------------
require(ggpubr)

# This is a trick to make sue the legend is plotted.
get_legend(PlotIndtoRF(ind_RF_LSU_mark, paletteCB4)) -> legend_mark
as_ggplot(legend_mark)

get_legend(PlotIndtoRF(ind_RF_ITS_niche, paletteCB2)) -> legend_niche
as_ggplot(legend_niche)


ggarrange(
  ggarrange(
    PlotIndtoRF(ind_RF_ITS_mark, paletteCB4) + 
      labs(title = "Top ITS OTUs to\nclassify markets",  x = "OTU", y = "Mean Decrease\nin Accuracy") +
      theme(legend.position = "none"),
    PlotIndtoRF(ind_RF_LSU_mark, paletteCB4) + 
      labs(title = "Top LSU OTUs to\nclassify markets",  x = "OTU", y = "Mean Decrease\nin Accuracy") +
      theme(legend.position = "none"),
    PlotIndtoRF(ind_RF_16s_mark, paletteCB4) + 
      labs(title = "Top 16s OTUs to\nclassify markets",  x = "OTU", y = "Mean Decrease\nin Accuracy") +
      theme(legend.position = "none"),
    as_ggplot(legend_mark),
    labels = c("A", "B", "C", ""),
    widths = c(1, 1, 1, 1),
    heights = c(1, 1, 1, 0.2),
    align = "hv" ,
    ncol = 1,
    nrow = 4),
  ggarrange(
    PlotIndtoRF(ind_RF_ITS_niche, paletteCB2) + 
      labs(title = "Top ITS OTUs to\nclassify niches",  x = NULL, y = "Mean Decrease\nin Accuracy") +
      theme(legend.position = "none"),
    PlotIndtoRF(ind_RF_LSU_niche, paletteCB2) + 
      labs(title = "Top LSU OTUs to\nclassify niches",  x = NULL, y = "Mean Decrease\nin Accuracy") +
      theme(legend.position = "none"),
    PlotIndtoRF(ind_RF_16s_niche, paletteCB2) + 
      labs(title = "Top 16s OTUs to\nclassify niches",  x = NULL, y = "Mean Decrease\nin Accuracy") +
      theme(legend.position = "none"),
    as_ggplot(legend_niche),
    widths = c(1, 1, 1, 1),
    heights = c(1, 1, 1, 0.2),
    align = "hv" ,
    ncol = 1,
    nrow = 4),
  widths = c(1, 1),
  align = "h",
  ncol = 2,
  nrow = 1,
  legend = c("bottom")) -> Fig_5

Fig_5


  
# ******************************************************---------------------------------------
# HEATMAP TREES -------------------------------------------------------------------------------
# ITS -----------------------------------------------------------------------------------------
library(metacoder)
library(tibble)

# otu table
otu_ITS$OTU_ID <- rownames(otu_ITS)
head(otu_ITS)

# metadata
metadata_ITS$sample_id <- rownames(metadata_ITS)
head(metadata_ITS)

# taxonomy
taxa_ITS$Domain <- rep("Eukaryota", nrow(taxa_ITS))
taxa_ITS$Kingdom <- as.character(taxa_ITS$Kingdom)
taxa_ITS$Phylum <- as.character(taxa_ITS$Phylum)
taxa_ITS$Class <- as.character(taxa_ITS$Class)
taxa_ITS$Order <- as.character(taxa_ITS$Order)
taxa_ITS$Family <- as.character(taxa_ITS$Family)
taxa_ITS$Genus <- as.character(taxa_ITS$Genus)
taxa_ITS$Species <- as.character(taxa_ITS$Species)
taxa_ITS[taxa_ITS=="Unclassified"]<- ""
taxa_ITS[is.na(taxa_ITS)]<- ""
head(taxa_ITS)

# Reformatting taxonomy
# taxa_ITS[, "Genus"] <- gsub(" sp.", "", taxa_ITS[, "Genus"])
cols <- c("Domain","Kingdom","Phylum","Class","Order","Family","Genus", "Species")
taxa_ITS$Taxonomy <- do.call(paste, c(taxa_ITS[cols], sep=";"))
head(taxa_ITS)

identical(colnames(otu_ITS), metadata_ITS$sample_id)
sample_order <- match(metadata_ITS$sample_id, colnames(otu_ITS))
otu_ITS <- otu_ITS[,sample_order]

identical(rownames(taxa_ITS), rownames(otu_ITS))
otu_ITS$Taxonomy <- taxa_ITS$Taxonomy
otu_ITS$OTU_ID <- rownames(otu_ITS)
head(otu_ITS)

# transform a data.frame into a tibble --------------------------------------------------------
tib_ITS <- tibble::as_tibble(otu_ITS)
# move OTU_ID to first column
tib_ITS %>% select(OTU_ID, everything()) -> tib_ITS

tib_meta_ITS <- tibble::as_tibble(metadata_ITS)
tib_meta_ITS

# create metacoder object 
heat_tree_ITS <- parse_tax_data(tib_ITS,
                                class_cols = "Taxonomy",
                                class_sep = c(";"))

# rename tax_data 
names(heat_tree_ITS$data) <- "otu_count"
heat_tree_ITS

# calculate taxa proprotions ------------------------------------------------------------------
heat_tree_ITS$data$otu_prop <- calc_obs_props(heat_tree_ITS,
                                              data = "otu_count",
                                              cols = tib_meta_ITS$sample_id)

heat_tree_ITS

# abundance per-taxon -------------------------------------------------------------------------
heat_tree_ITS$data$tax_prop <- calc_taxon_abund(heat_tree_ITS, 
                                                 "otu_prop", 
                                                 cols = tib_meta_ITS$sample_id)

# number of samples have reads for each taxon --------------------------------------------------
heat_tree_ITS$data$tax_occ_type <- calc_n_samples(heat_tree_ITS, 
                                             "tax_prop", 
                                             cols = tib_meta_ITS$sample_id,
                                             groups = tib_meta_ITS$Type)

heat_tree_ITS$data$tax_occ_mark <- calc_n_samples(heat_tree_ITS,
                                             "tax_prop",
                                             cols = tib_meta_ITS$sample_id,
                                             groups = tib_meta_ITS$Market)

heat_tree_ITS

# *** FIGURE 5A heat tree ITS -----------------------------------------------------------------
set.seed(31)

heat_tree(heat_tree_ITS, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = heat_tree_ITS$data$tax_occ_type$`in`,
          node_label_size_range = c(0.0125, 0.05), 
          tree_label_size = as.numeric(0.5),
          node_size_range = c(0.01, 0.05),
          edge_size_range = c(0.003, 0.003),
          node_color_range = palette_heat,
          #node_color_trans = "linear",
          #node_color_interval = c(0, 36), 
          node_size_axis_label = "OTU number",
          node_color_axis_label = "Samples with counts",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford") -> plot_tree_ITS # The layout algorithm that initializes node locations

plot_tree_ITS

# # Comparing groups --------------------------------------------------------------------------
# heat_tree_ITS$data$diff_table <- compare_groups(heat_tree_ITS, dataset = "tax_abund",
#                                                 cols = tib_meta_ITS$sample_id, # What columns of sample data to use
#                                                 groups = tib_meta_ITS$Market) # What category each sample is assigned to
# print(heat_tree_ITS$data$diff_table)
# 
# heat_tree_ITS$data$diff_table$wilcox_p_value <- p.adjust(heat_tree_ITS$data$diff_table$wilcox_p_value,
#                                                          method = "fdr")
# 
# 
# heat_tree_matrix(heat_tree_ITS,
#                  data = "diff_table",
#                  node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
#                  node_label = taxon_names,
#                  node_color = log2_median_ratio, # A column from `obj$data$diff_table`
#                  node_color_range = diverging_palette(), # The built-in palette for diverging data
#                  node_color_trans = "linear", # The default is scaled by circle area
#                  node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
#                  edge_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
#                  node_size_axis_label = "Number of OTUs",
#                  node_color_axis_label = "Log2 ratio median proportions",
#                  layout = "davidson-harel", # The primary layout algorithm
#                  initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
#                  output_file = "differential_heat_tree.pdf") # Saves the plot as a pdf file
# 
# 
# heat_tree_ITS$data$diff_type <- compare_groups(heat_tree_ITS, dataset = "tax_abund",
#                                                 cols = tib_meta_ITS$sample_id, # What columns of sample data to use
#                                                 groups = tib_meta_ITS$Type) # What category each sample is assigned to
# print(heat_tree_ITS$data$diff_type)
# 
# heat_tree_ITS$data$diff_type$wilcox_p_value <- p.adjust(heat_tree_ITS$data$diff_type$wilcox_p_value,
#                                                          method = "fdr")
# 
# 
# heat_tree_matrix(heat_tree_ITS,
#                  data = "diff_type",
#                  node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
#                  node_label = taxon_names,
#                  node_color = log2_median_ratio, # A column from `obj$data$diff_type`
#                  node_color_range = diverging_palette(), # The built-in palette for diverging data
#                  node_color_trans = "linear", # The default is scaled by circle area
#                  node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
#                  edge_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
#                  node_size_axis_label = "Number of OTUs",
#                  node_color_axis_label = "Log2 ratio median proportions",
#                  layout = "davidson-harel", # The primary layout algorithm
#                  initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
#                  output_file = "differential_heat_tree.pdf") # Saves the plot as a pdf file

# LSU -----------------------------------------------------------------------------------------
merge_phyloseq(physeq_Ani_scaled_filt, physeq_LSU_scaled_filt) -> physeq_LSU_all
physeq_LSU_all

otu_LSU_all <- as.data.frame(otu_table(physeq_LSU_all))
metadata_LSU_all = as(sample_data(physeq_LSU_all), "data.frame")
taxa_LSU_all <- as.data.frame(as.matrix(tax_table(physeq_LSU_all)))
dim(taxa_LSU_all)

# otu table
otu_LSU_all$OTU_ID <- rownames(otu_LSU_all)
head(otu_LSU_all)

# metadata
metadata_LSU_all$sample_id <- rownames(metadata_LSU_all)
head(metadata_LSU_all)

#taxa_LSU_all -> taxa_old_lsu

# taxonomy
taxa_LSU_all$Domain <- rep("Eukaryota", nrow(taxa_LSU_all))
taxa_LSU_all$Kingdom <- as.character(taxa_LSU_all$Kingdom)
taxa_LSU_all$Phylum <- as.character(taxa_LSU_all$Phylum)
taxa_LSU_all$Class <- as.character(taxa_LSU_all$Class)
taxa_LSU_all$Order <- as.character(taxa_LSU_all$Order)
taxa_LSU_all$Family <- as.character(taxa_LSU_all$Family)
taxa_LSU_all$Genus <- as.character(taxa_LSU_all$Genus)
taxa_LSU_all$Species <- as.character(taxa_LSU_all$Species)
taxa_LSU_all[taxa_LSU_all=="Unclassified"]<- ""
taxa_LSU_all[is.na(taxa_LSU_all)]<- ""
head(taxa_LSU_all)

cols <- c("Domain","Kingdom","Phylum","Class","Order","Family","Genus", "Species")
taxa_LSU_all$Taxonomy <- do.call(paste, c(taxa_LSU_all[cols], sep=";"))
head(taxa_LSU_all)

identical(colnames(otu_LSU_all), metadata_LSU_all$sample_id)
sample_order <- match(metadata_LSU_all$sample_id, colnames(otu_LSU_all))
otu_LSU_all <- otu_LSU_all[,sample_order]

identical(rownames(taxa_LSU_all), rownames(otu_LSU_all))
otu_LSU_all$Taxonomy <- taxa_LSU_all$Taxonomy
head(otu_LSU_all)

# transform a data.frame into a tibble --------------------------------------------------------
# transform a data.frame into a tibble --------------------------------------------------------
tib_LSU <- tibble::as_tibble(otu_LSU_all)

# move OTU_ID to first column
tib_LSU %>% select(OTU_ID, everything()) -> tib_LSU

tib_meta_LSU <- tibble::as_tibble(metadata_LSU_all)
tib_meta_LSU

# create metacoder object 
heat_tree_LSU <- parse_tax_data(tib_LSU,
                                class_cols = "Taxonomy",
                                class_sep = c(";"))

# rename tax_data 
names(heat_tree_LSU$data) <- "otu_count"
heat_tree_LSU

# calculate taxa proprotions ------------------------------------------------------------------
heat_tree_LSU$data$otu_prop <- calc_obs_props(heat_tree_LSU,
                                              data = "otu_count",
                                              cols = tib_meta_LSU$sample_id)

heat_tree_LSU

# abundance per-taxon -------------------------------------------------------------------------
heat_tree_LSU$data$tax_prop <- calc_taxon_abund(heat_tree_LSU, 
                                                "otu_prop", 
                                                cols = tib_meta_LSU$sample_id)

# number of samples have reads for each taxon --------------------------------------------------
heat_tree_LSU$data$tax_occ_type <- calc_n_samples(heat_tree_LSU, 
                                                  "tax_prop", 
                                                  cols = tib_meta_LSU$sample_id,
                                                  groups = tib_meta_LSU$Type)

heat_tree_LSU$data$tax_occ_mark <- calc_n_samples(heat_tree_LSU,
                                                  "tax_prop",
                                                  cols = tib_meta_LSU$sample_id,
                                                  groups = tib_meta_LSU$Market)

heat_tree_LSU

# *** FIGURE 5B heat tree LSU -----------------------------------------------------------------
set.seed(34)

heat_tree(heat_tree_LSU, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = heat_tree_LSU$data$tax_occ_type$`in`,
          node_label_size_range = c(0.0125, 0.05), 
          tree_label_size = as.numeric(0.5),
          node_size_range = c(0.01, 0.05),
          edge_size_range = c(0.003, 0.003),
          node_color_range = palette_heat,
          #node_color_trans = "linear",
          node_color_interval = c(0, 36), 
          node_size_axis_label = "OTU number",
          node_color_axis_label = "Samples with counts",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford") -> plot_tree_LSU # The layout algorithm that initializes node locations

plot_tree_LSU


# 16S -----------------------------------------------------------------------------------------
# otu table
otu_16s$OTU_ID <- rownames(otu_16s)
head(otu_16s)

# metadata
metadata_16s$sample_id <- rownames(metadata_16s)
head(metadata_16s)

# taxonomy
taxa_16s$Domain <- rep("Prokaryota", nrow(taxa_16s))
taxa_16s$Kingdom <- as.character(taxa_16s$Kingdom)
taxa_16s$Phylum <- as.character(taxa_16s$Phylum)
taxa_16s$Class <- as.character(taxa_16s$Class)
taxa_16s$Order <- as.character(taxa_16s$Order)
taxa_16s$Family <- as.character(taxa_16s$Family)
taxa_16s$Genus <- as.character(taxa_16s$Genus)
taxa_16s$Species <- as.character(taxa_16s$Species)
taxa_16s[taxa_16s=="Unclassified"]<- ""
taxa_16s[is.na(taxa_16s)]<- ""
head(taxa_16s)

# Reformatting taxonomy
# taxa_16s[, "Genus"] <- gsub(" sp.", "", taxa_16s[, "Genus"])
cols <- c("Domain","Kingdom","Phylum","Class","Order","Family","Genus", "Species")
taxa_16s$Taxonomy <- do.call(paste, c(taxa_16s[cols], sep=";"))
head(taxa_16s)

identical(colnames(otu_16s), metadata_16s$sample_id)
sample_order <- match(metadata_16s$sample_id, colnames(otu_16s))
otu_16s <- otu_16s[,sample_order]

identical(rownames(taxa_16s), rownames(otu_16s))
otu_16s$Taxonomy <- taxa_16s$Taxonomy
head(otu_16s)

# transform a data.frame into a tibble --------------------------------------------------------
# transform a data.frame into a tibble --------------------------------------------------------
tib_16s <- tibble::as_tibble(otu_16s)
# move OTU_ID to first column
tib_16s %>% select(OTU_ID, everything()) -> tib_16s

tib_meta_16s <- tibble::as_tibble(metadata_16s)
tib_meta_16s

# create metacoder object 
heat_tree_16s <- parse_tax_data(tib_16s,
                                class_cols = "Taxonomy",
                                class_sep = c(";"))

# rename tax_data 
names(heat_tree_16s$data) <- "otu_count"
heat_tree_16s

# calculate taxa proprotions ------------------------------------------------------------------
heat_tree_16s$data$otu_prop <- calc_obs_props(heat_tree_16s,
                                              data = "otu_count",
                                              cols = tib_meta_16s$sample_id)

heat_tree_16s

# abundance per-taxon -------------------------------------------------------------------------
heat_tree_16s$data$tax_prop <- calc_taxon_abund(heat_tree_16s, 
                                                "otu_prop", 
                                                cols = tib_meta_16s$sample_id)

# number of samples have reads for each taxon --------------------------------------------------
heat_tree_16s$data$tax_occ_type <- calc_n_samples(heat_tree_16s, 
                                                  "tax_prop", 
                                                  cols = tib_meta_16s$sample_id,
                                                  groups = tib_meta_16s$Type)

heat_tree_16s$data$tax_occ_mark <- calc_n_samples(heat_tree_16s,
                                                  "tax_prop",
                                                  cols = tib_meta_16s$sample_id,
                                                  groups = tib_meta_16s$Market)

heat_tree_16s

# *** FIGURE 5C heat tree 16S -----------------------------------------------------------------
set.seed(33) 

heat_tree_16s %>%
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = heat_tree_16s$data$tax_occ_type$`in`,
            node_label_size_range = c(0.0125, 0.05), 
            tree_label_size = as.numeric(0.5),
            node_size_range = c(0.01, 0.05),
            edge_size_range = c(0.003, 0.003),
            node_color_range = palette_heat,
            #node_color_trans = "linear",
            node_color_interval = c(0, 36), 
            node_size_axis_label = "OTU number",
            node_color_axis_label = "Samples with counts",
            layout = "davidson-harel",
            initial_layout = "reingold-tilford") -> plot_tree_16s 

plot_tree_16s


# EXTRACT HEAT-TREE DATA ----------------------------------------------------------------------
ExractHeatData <- function(heat_tree){
  heat_tree$taxon_names() -> taxa_names
      taxa_names <- as.data.frame(taxa_names)
  heat_tree$n_obs() -> n_obs
      n_obs <- as.data.frame(n_obs)
  heat_tree$data$tax_occ_type -> tax_occ_type
      tax_occ_type <- as.data.frame(tax_occ_type)
      rownames(tax_occ_type) <- tax_occ_type$taxon_id
  heat_tree$data$tax_occ_mark -> tax_occ_mark
      tax_occ_mark <- as.data.frame(tax_occ_mark)
      rownames(tax_occ_mark) <- tax_occ_mark$taxon_id
  heat_tree$data$tax_prop-> tax_prop
      tax_prop <- as.data.frame(tax_prop)
      rownames(tax_prop) <- tax_prop$taxon_id
  tax_prop$taxa_abund <-
      rowSums(tax_prop[, 2:ncol(tax_prop)]) / rowSums(tax_prop[, 2:ncol(tax_prop)])[1] *
      100
  dplyr::bind_cols(list(taxa_names, n_obs, tax_occ_type, tax_occ_mark, tax_prop)) -> df_heat
  df_heat <- as.data.frame(df_heat[, c(1,2,4,5, 7:10,84)])
  arrange(df_heat, desc(taxa_abund)) -> df_heat
  head(df_heat)
  #(df_heat$'in'* 100)/36
return(df_heat)
}


ExractHeatData(heat_tree = heat_tree_ITS) -> heat_table_ITS
ExractHeatData(heat_tree = heat_tree_LSU) -> heat_table_LSU
ExractHeatData(heat_tree = heat_tree_16s) -> heat_table_16s

write.csv(heat_table_ITS, "heat_table_ITS.csv")
write.csv(heat_table_LSU, "heat_table_LSU.csv")
write.csv(heat_table_16s, "heat_table_16s.csv")


# EXTRACT SEQUENCES ---------------------------------------------------------------------------
tax_table(subset_taxa(physeq_ITS_uparse_R1, Order=="Mortierellales"))

write.dna(refseq(subset_taxa(physeq_ITS_uparse_R1, Order=="Mortierellales")), 
          format="fasta", 
          colsep="",
          file="mortierellales.fasta")
