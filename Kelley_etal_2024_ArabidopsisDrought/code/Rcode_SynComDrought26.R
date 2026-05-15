# *********************************************---------------------------------
# ************ USEARCH ANALYSIS ***************---------------------------------#
# Manuscript:   The impact of drought on microbiome assembly in Arabidopsis
# Authors:      ...
# Affiliation:  1 Department of Plant Soil and Microbial Sciences, Michigan State University, East Lansing MI 48824
#               2 Great Lakes Bioenergy Research Center, Michigan State University, East Lansing MI 48824
#               ...
# Journal:      ...
# Date:         April 1, 2026
# *********************************************---------------------------------

# R options --------------------------------------------------------------------
options(scipen = 9999, pillar.sigfig=6, digits=6, max.print=10000000) 
#rm(list = ls())

# R packages -------------------------------------------------------------------
suppressWarnings({
  if (!require(pacman)) install.packages("pacman")
})

pacman::p_load(
  conflicted,
  styler,
  tidyverse,
  vegan,
  ggpubr,
  ggtext,
  tidytext,
  phyloseq,
  speedyseq,
  Biostrings,
  metagenomeSeq,
  decontam,
  mvabund,
  gridExtra, 
  install=FALSE
)

# Import tables ----------------------------------------------------------------
asv_table <-
  read.csv("datasets/asvtable_UNOISE_233bp.txt", row.names = 1, sep = "\t")
head(asv_table)
colnames(asv_table)
dim(asv_table)

otu_table_uparse <-
  read.csv("datasets/otutable_UPARSE_233bp.txt", row.names = 1, sep = "\t")
head(otu_table_uparse)
colnames(otu_table_uparse)
dim(otu_table_uparse)

otu_table_closedref <-
  read.csv("datasets/otutable_233bp_closedRef.txt", row.names = 1, sep = "\t")
head(otu_table_closedref)
colnames(otu_table_closedref)
dim(otu_table_closedref)

otu_table_from_asv <-
  read.csv("datasets/otutable_asv_233bp_to97otus.txt", row.names = 1, sep = "\t")
head(otu_table_from_asv)
colnames(otu_table_from_asv)
dim(otu_table_from_asv)

otu_table_swarm <-
  read.csv("datasets/otutable_SWARM_233bp.txt", row.names = 1, sep = "\t")
head(otu_table_swarm)
colnames(otu_table_swarm)
dim(otu_table_swarm)





# Metadata ---------------------------------------------------------------------
metadata_usearch <-
  read.csv("dataset/R1/metadata.csv", sep = ",") %>% 
  dplyr::rename(SampleName=1)

metadata_usearch
head(metadata_usearch)
dim(metadata_usearch)

renamed_sampl <-
  read.delim("dataset/R1/renamed_sample.txt", header = FALSE, sep="\t") %>% 
  dplyr::rename(Fastq = V1, SampleID = V2) %>% 
  separate(Fastq, c("SampleName", "FastqID"), sep = "_", remove = FALSE) %>% 
  mutate(SampleName = gsub("-", "", SampleName))

renamed_sampl
head(renamed_sampl)
dim(renamed_sampl)
length(sort(unique(renamed_sampl$FastqID)))
renamed_sampl

metadata <-
  metadata_usearch %>% 
  full_join(., renamed_sampl, by="SampleName") %>% 
  column_to_rownames("SampleID")

head(metadata)

# Taxonomy ---------------------------------------------------------------------
constax_asv <-
  read.delim("dataset/R1/asv_constax_taxonomy.txt", row.names=1, header = TRUE,sep="\t") %>% 
  dplyr::rename(
    domain = Rank_1, phylum = Rank_2, class = Rank_3,
    order = Rank_4, family = Rank_5, genus = Rank_6,
    species = Rank_7) %>% 
  dplyr::select(-Rank_8) %>% 
  mutate(across(c("domain", "phylum", "class","order",
                  "family","genus", "species"), ~gsub("_1", "", x=.))) 

constax_asv
head(constax_asv)





