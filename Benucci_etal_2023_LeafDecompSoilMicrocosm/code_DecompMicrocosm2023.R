# *** Script for data analyses ******************************** ----------------
# Manuscript:   Soil microbiome in detritosphere hot-spots are driven by soil pores and water contents
# Authors:      ...
# Affiliation:  Michigan State University
# Journal:      ... 
# Date:         March 24, 2022
# ******************************************************************** ---------

options(scipen = 999) 
options(max.print=100000) 
library(styler)
library(targets)
library(job)

### SETTING UP WORKING ENVIRONEMNT ---------------------------------------------
library(tidyverse)
library(tidymodels)
library(phyloseq)
library(Biostrings)
library(ape)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(magrittr)
library(ggpubr)
library(patchwork)
library(vegan)
library(grid)
library(gridExtra)
library(multcompView)
library(indicspecies)

library(ggtext)

# PALETTE ----------------------------------------------------------------------
paletteCB3 = c("#599861","#FF934F","black")
pie(rep(1, length(paletteCB3)), labels = sprintf("%d (%s)",
     seq_along(paletteCB3),paletteCB3), col = paletteCB3)

paletteCB5 = c("#599861","#FF934F","#825121","#CC2D35","black")
pie(rep(1, length(paletteCB5)), labels = sprintf("%d (%s)",
      seq_along(paletteCB5),paletteCB5), col = paletteCB5)

paletteCB5_bis = c("#60ffaf", "#058ED9","#521899","#CC2D35","black")
pie(rep(1, length(paletteCB5_bis)), labels = sprintf("%d (%s)",
        seq_along(paletteCB5_bis),paletteCB5_bis), col = paletteCB5_bis)

paletteCB5_tris = c("#CCCCCC","#999999","#666666","#000000","#CC2D35")
pie(rep(1, length(paletteCB5_tris)), labels = sprintf("%d (%s)",
      seq_along(paletteCB5_tris),paletteCB5_tris), col = paletteCB5_tris)

palette_CB30 = c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77","#283dff","#058ED9","#636bb7",
                  "#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
                  "#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#5b5b19",
                  "#fcfc00","#ffff9e","#ffb7ef","#fa7efc","#ae09ea","#521899","#000000")
pie(rep(1, length(palette_CB30)), labels = sprintf("%d (%s)",
       seq_along(palette_CB30),palette_CB30), col = palette_CB30)


#"#60ffaf","#82807f", "#3f3e3d","#CC2D35","black"
#"#60ffaf", "#283dff","#636bb7","#CC2D35","black"
#"#058ED9","#848FA2","#2D3142","#CC2D35","black"

paletteCB6 = c("#2D3142","#058ED9","#848FA2","#599861","#FF934F", "#CC2D35")
pie(rep(1, length(paletteCB6)), labels = sprintf("%d (%s)",
        seq_along(paletteCB6),paletteCB6), col = paletteCB6)


paletteCB7 = c("#599861","#FF934F","#825121","#CC2D35","black","#fcfc00", "#058ED9", "#848FA2", "#ae09ea")
pie(rep(1, length(paletteCB7)), labels = sprintf("%d (%s)",
       seq_along(paletteCB7),paletteCB7), col = paletteCB7)


palette_new30 = c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77","#283dff","#636bb7",
                  "#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
                  "#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#5b5b19",
                  "#fcfc00","#ffff9e","#ffb7ef","#fa7efc","#ae09ea","#521899","#000000") #"#1e0047"
pie(rep(1, length(palette_new30)), labels = sprintf("%d (%s)",
     seq_along(palette_new30),palette_new30), col = palette_new30)

paletteCB2 <- c("#60ffaf", "#ffb7ef")
pie(rep(1, length(paletteCB2)), labels = sprintf("%d (%s)",
     seq_along(paletteCB2),paletteCB2), col = paletteCB2)

# ****************** PREFILTERING *************************---------------------
# Fungi ------------------------------------------------------------------------

otus_ITS_uparse <- read.delim("ITS/otu_table.txt",row.names=1, header=TRUE, sep="\t") 
head(otus_ITS_uparse)
dim(otus_ITS_uparse)

metadata_ITS_uparse <-read.delim("ITS/metadata_fungi_corrected_TS.txt", row.names=1, header=TRUE, sep="\t")
head(metadata_ITS_uparse)

constax_taxonomy_07 <-read.delim("ITS/constax_taxonomy.txt",header=TRUE, row.names=1, sep="\t")
head(constax_taxonomy_07)
dim(constax_taxonomy_07)

constax_taxonomy_07_filt <-
  subset(
    constax_taxonomy_07,
    High_level_taxonomy == "Fungi" &
      HL_hit_percent_id >= 60 &
      HL_hit_query_cover >= 85
  )

# check for non-target amplifications
dim(constax_taxonomy_07_filt)
table(constax_taxonomy_07_filt$High_level_taxonomy)


otus_seq_ITS_uparse <- readDNAStringSet("ITS/otus_new.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

physeq_ITS_uparse <- phyloseq(otu_table(otus_ITS_uparse, taxa_are_rows = TRUE),
                              sample_data(metadata_ITS_uparse),
                              tax_table(as.matrix(constax_taxonomy_07_filt)),
                              otus_seq_ITS_uparse) 

physeq_ITS_uparse
tax_table(physeq_ITS_uparse)
sample_data(physeq_ITS_uparse)
tax_table(physeq_ITS_uparse)[tax_table(physeq_ITS_uparse)==""]<- NA
#tax_table(physeq_ITS_uparse)[is.na(tax_table(physeq_ITS_uparse))]<-"Unclassified"

# rename OTUs
paste("FOTU", 1:length(taxa_names(physeq_ITS)), sep="_") %>% as.data.frame()

physeq_ITS_uparse -> physeq_ITS
taxa_names(physeq_ITS) <- paste("FOTU", 1:length(taxa_names(physeq_ITS)), sep="_")
otu_table(physeq_ITS)

write.csv(
as.data.frame(taxa_names(physeq_ITS_uparse)) %>% 
mutate(OTUID = paste("FOTU", 1:length(taxa_names(physeq_ITS_uparse)), sep="_")),
"otu_list_fungi.csv")


# removing samples
physeq_ITS_filt <- subset_samples(physeq_ITS, to_filter%in%c("no"))
otu_table(physeq_ITS_filt) <- otu_table(physeq_ITS_filt)[which(rowSums(otu_table(physeq_ITS_filt)) >= 2),] 
physeq_ITS_filt
head(tax_table(physeq_ITS_filt))
head(sample_data(physeq_ITS_filt))

write.csv(physeq_ITS_filt@sam_data, "finalmap_ITS.csv")
sample_data(physeq_ITS_filt) <- 
  read.csv("finalmap_ITS_new.csv", header = TRUE, row.names = 1) 

# Sequencing results and contaminants
sum(sample_data(physeq_ITS_filt)$read_number)
sample_data(physeq_ITS_filt)$SampleID <- rownames(sample_data(physeq_ITS_filt))
sample_data(physeq_ITS_filt)
table(sample_data(physeq_ITS_filt)$Treatment)

sum(sample_sums(
  subset_samples(physeq_ITS_filt, Slice1 %in% c("Leaf", "Nearsoil", "Farsoil"))))

mean(sample_sums(
  subset_samples(physeq_ITS_filt, Slice1 %in% c("Leaf", "Nearsoil", "Farsoil"))))

sd(sample_sums(
  subset_samples(physeq_ITS_filt, Slice1 %in% c("Leaf", "Nearsoil", "Farsoil"))))

tax_table(physeq_ITS_filt)[is.na(tax_table(physeq_ITS_filt))]<-"Unclassified"
unique(as.data.frame(tax_table(physeq_ITS_filt))$Kingdom)
unique(as.data.frame(tax_table(physeq_ITS_filt))$High_level_taxonomy)

# Prokaryotes ------------------------------------------------------------------
# import RDS file from Terry ---------------------------------------------------
readRDS(file = "16S/sasha_16S_min_otu_5_min_14600_seqs_rarefied_to_14800_samples_only_physeq.RDS") -> physeq_16s
physeq_16s
head(tax_table(physeq_16s))

# rename OTUs
taxa_names(physeq_16s) <- paste("POTU", 1:length(taxa_names(physeq_16s)), sep="_")
head(otu_table(physeq_16s))

# add new metadata file
map_16s <- read.csv("16S/sample_data_16s_corrected.csv", header=T, row.names =1)
sample_data(physeq_16s) <- sample_data(map_16s)
sample_data(physeq_16s)
physeq_16s

# removing samples
physeq_16s_filt <- subset_samples(physeq_16s, to_filter%in%c("no"))
otu_table(physeq_16s_filt) <- otu_table(physeq_16s_filt)[which(rowSums(otu_table(physeq_16s_filt)) >= 2),] 
physeq_16s_filt
head(sample_data(physeq_16s_filt))
head(tax_table(physeq_16s_filt))

write.csv(physeq_16s_filt@sam_data, "finalmap_16s.csv")
sample_data(physeq_16s_filt) <- 
  read.csv("finalmap_16s_new.csv", header = TRUE, row.names = 1) 


# removing non-target amplifications
apply(tax_table(physeq_16s_filt), 2, function(x) which(x == "Chloroplast"))
apply(tax_table(physeq_16s_filt), 2, function(x) which(x == "Cyanobacteria/Chloroplast"))
apply(tax_table(physeq_16s_filt), 2, function(x) which(x == "Mitochondria"))
apply(tax_table(physeq_16s_filt), 2, function(x) which(x == "Oceanospirillales"))
apply(tax_table(physeq_16s_filt), 2, function(x) which(x = is.na(x)))

rownames(tax_table(
  subset_taxa(
    physeq_16s_filt,
    class == "Chloroplast" |
      family == "Chloroplast" |
      phylum == "Cyanobacteria/Chloroplast" |
      phylum == "Mitochondria" |
      order == "Oceanospirillales" |
      is.na(domain)))) -> bad_taxa_16s

remove_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

physeq_16s_filt <- remove_taxa(physeq_16s_filt, bad_taxa_16s)
otu_table(physeq_16s_filt) <- otu_table(physeq_16s_filt)[which(rowSums(otu_table(physeq_16s_filt)) > 0),] 
physeq_16s_filt

# Sequencing results and contaminants for fungi
sample_data(physeq_16s_filt)$SampleID <- rownames(sample_data(physeq_16s_filt))
sample_data(physeq_16s_filt)

sum(sample_sums(
  subset_samples(physeq_16s_filt, Slice1 %in% c("Leaf", "Nearsoil", "Farsoil"))))

mean(sample_sums(
  subset_samples(physeq_16s_filt, Slice1 %in% c("Leaf", "Nearsoil", "Farsoil"))))

sd(sample_sums(
  subset_samples(physeq_16s_filt, Slice1 %in% c("Leaf", "Nearsoil", "Farsoil"))))

tax_table(physeq_16s_filt)[is.na(tax_table(physeq_16s_filt))]<-"Unclassified"
unique(as.data.frame(tax_table(physeq_16s_filt))$Kingdom)
unique(as.data.frame(tax_table(physeq_16s_filt))$High_level_taxonomy)


# ******************** DATA ANALYSIS ****************************---------------
# Filter to project data sets --------------------------------------------------
# > Data set WITHOUT controls --------------------------------------------------
physeq_ITS_filt@sam_data
physeq_ITS_filt_no_contr <- subset_samples(physeq_ITS_filt, Treatment1%in%c("T1", "T4") & Plant%in%c("Corn", "Soy"))
otu_table(physeq_ITS_filt_no_contr) <- otu_table(physeq_ITS_filt_no_contr)[which(rowSums(otu_table(physeq_ITS_filt_no_contr)) > 0),]
physeq_ITS_filt_no_contr
sort(sample_sums(physeq_ITS_filt_no_contr), decreasing = TRUE)

physeq_16s_filt_no_contr <- subset_samples(physeq_16s_filt, Treatment1%in%c("T1", "T4") & Plant%in%c("Corn", "Soy"))
otu_table(physeq_16s_filt_no_contr) <- otu_table(physeq_16s_filt_no_contr)[which(rowSums(otu_table(physeq_16s_filt_no_contr)) > 0),]
physeq_16s_filt_no_contr
sort(sample_sums(physeq_16s_filt_no_contr), decreasing = TRUE)

# > Data set WITH controls -----------------------------------------------------
physeq_ITS_filt@sam_data
physeq_ITS_ctr <- subset_samples(physeq_ITS_filt, Treatment1%in%c("T1","T4","Start","Control"))
otu_table(physeq_ITS_ctr) <- otu_table(physeq_ITS_ctr)[which(rowSums(otu_table(physeq_ITS_ctr)) > 0),]
physeq_ITS_ctr
sort(sample_sums(physeq_ITS_ctr), decreasing = TRUE)

physeq_16s_filt@sam_data
physeq_16s_ctr <- subset_samples(physeq_16s_filt, Treatment1%in%c("T1","T4","Start","Control"))
otu_table(physeq_16s_ctr) <- otu_table(physeq_16s_ctr)[which(rowSums(otu_table(physeq_16s_ctr)) > 0),]
physeq_16s_ctr
sort(sample_sums(physeq_16s_ctr), decreasing = TRUE)
table(physeq_16s_ctr@sam_data$Slice)

# RAREFYING EVEN DEPTH ---------------------------------------------------------
# > Data set WITHOUT controls --------------------------------------------------
#fungi
set.seed(2021)
sort(colSums(otu_table(physeq_ITS_filt_no_contr)), decreasing = TRUE)
physeq_fungi_ev = rarefy_even_depth(physeq_ITS_filt_no_contr, sample.size = 1010, 
                                    rngseed = FALSE, replace = TRUE, 
                                    trimOTUs = TRUE, verbose = TRUE) 
otu_table(physeq_fungi_ev) <- otu_table(physeq_fungi_ev)[which(rowSums(otu_table(physeq_fungi_ev)) >= 1),]
physeq_fungi_ev
colSums(otu_table(physeq_fungi_ev))
any(taxa_sums(physeq_fungi_ev) == 0)


head(sample_data(physeq_fungi_ev))
otu_fungi_ev <- as.data.frame(otu_table(physeq_fungi_ev))
taxa_fungi_ev <- as.data.frame(as.matrix(tax_table(physeq_fungi_ev)))
metadata_fungi_ev <- as.data.frame(as.matrix(sample_data(physeq_fungi_ev)))
head(metadata_fungi_ev)


# prokaryotes
sort(colSums(otu_table(physeq_16s_filt_no_contr)), decreasing = TRUE)
physeq_prok_ev = rarefy_even_depth(physeq_16s_filt_no_contr, sample.size = 13377, 
                                   rngseed = FALSE, replace = TRUE, 
                                   trimOTUs = TRUE, verbose = TRUE) 
otu_table(physeq_prok_ev) <- otu_table(physeq_prok_ev)[which(rowSums(otu_table(physeq_prok_ev)) >= 1),]
physeq_prok_ev
colSums(otu_table(physeq_prok_ev))
any(taxa_sums(physeq_prok_ev) == 0)


otu_prok_ev <- as.data.frame(otu_table(physeq_prok_ev))
taxa_prok_ev <- as(tax_table(physeq_prok_ev), "matrix")
taxa_prok_ev <- as.data.frame(taxa_prok_ev)
metadata_prok_ev <- as.data.frame(as.matrix(sample_data(physeq_prok_ev)))

# > Data set WITH controls -----------------------------------------------------
physeq_ITS_ctr@sam_data

physeq_fungi_ctr_ev = rarefy_even_depth(physeq_ITS_ctr, sample.size = 1010, 
                                    rngseed = FALSE, replace = TRUE, 
                                    trimOTUs = TRUE, verbose = TRUE) 
otu_table(physeq_fungi_ctr_ev) <- 
  otu_table(physeq_fungi_ctr_ev)[which(rowSums(otu_table(physeq_fungi_ctr_ev)) >= 1),]
physeq_fungi_ctr_ev
colSums(otu_table(physeq_fungi_ctr_ev))
any(taxa_sums(physeq_fungi_ctr_ev) == 0)

physeq_16s_ctr@sam_data

sort(colSums(otu_table(physeq_16s_ctr)), decreasing = TRUE)
physeq_prok_ctr_ev = rarefy_even_depth(physeq_16s_ctr, sample.size = 13377, 
                                        rngseed = FALSE, replace = TRUE, 
                                        trimOTUs = TRUE, verbose = TRUE) 
otu_table(physeq_prok_ctr_ev) <- 
  otu_table(physeq_prok_ctr_ev)[which(rowSums(otu_table(physeq_prok_ctr_ev)) >= 1),]
physeq_prok_ctr_ev
colSums(otu_table(physeq_prok_ctr_ev))
any(taxa_sums(physeq_prok_ctr_ev) == 0)


write.dna(refseq(physeq_fungi_ctr_ev), 
          format="fasta", 
          colsep="",
          file="fungal_otus.fasta")

# Reload dataframes ------------------------------------------------------------
#saveRDS(physeq_prok_ev, "bacteria_physeq.rds")
#saveRDS(physeq_fungi_ev, "fungi_physeq.rds")

readRDS("Random_Forest/fungi_physeq.rds") -> physeq_fungi_ev
physeq_fungi_ev
sampleSums(physeq_fungi_ev)
table(physeq_fungi_ev@sam_data$Slice)

readRDS("Random_Forest/bacteria_physeq.rds") -> physeq_prok_ev
physeq_prok_ev
sampleSums(physeq_prok_ev)

# BETA DIVERSITY----------------------------------------------------------------

# Calculating PERMANOVAs
Adjadonis <- function(physeq, formula, strata=NULL){
  require(tidyverse)
  otu <- as.data.frame(otu_table(physeq))
  metadata <- as.matrix(sample_data(physeq))
  metadata <- as.data.frame(metadata)
  metadata_new <- metadata %>% 
    mutate(across(where(is_character),as_factor))
  print(str(metadata_new))
  adonis_test <-
    adonis(formula,
           strata = NULL,
           metadata_new, method = "bray",
           permutations=9999)
  df_adonis <- as.data.frame(adonis_test$aov.tab)
  df_adonis$Padj <-
    p.adjust(adonis_test$aov.tab$`Pr(>F)`, "bonferroni")
  df_adonis <-
    subset(df_adonis, Padj<=0.05)
  ad_list <-
    list(adonis_test, df_adonis)
  #print(adonis_test)
  return(ad_list)
}


Adjadonis(physeq_fungi_ctr_start, t(otu) ~ Treatment * Pore)

# ***********************************************---------------------------------
# > DRY CONTROLS - START ---------------------------------------------------------
# PERMANOVA for Fungi ----------------------------------------------------------
physeq_fungi_ctr_ev@sam_data

physeq_fungi_ctr_start <-
  physeq_fungi_ctr_ev %>% 
  subset_samples(physeq = .,Time %in% c("0"))
otu_table(physeq_fungi_ctr_start) <- 
  otu_table(physeq_fungi_ctr_start)[which(rowSums(otu_table(physeq_fungi_ctr_start)) > 0),]
physeq_fungi_ctr_start

physeq_fungi_ctr_start@sam_data

adonis_fungi_ctr_start <-
  Adjadonis(physeq_fungi_ctr_start, t(otu) ~ Treatment * Pore)
adonis_fungi_ctr_start[[2]]
adonis_fungi_ctr_start[[1]]$aov.tab

adonis_fungi_ctr_start <-
  Adjadonis(physeq_fungi_ctr_start, t(otu) ~ Pore * Treatment)
adonis_fungi_ctr_start[[2]]
adonis_fungi_ctr_start[[1]]$aov.tab

write.csv(adonis_fungi_ctr_start[[1]]$aov.tab, "New_Manuscript_draft_2023/adonis-dry-contr-fungi.csv")

# PERMANOVA for Bacteria ---------------------------------------------------------
physeq_prok_ctr_ev@sam_data

physeq_prok_ctr_start <-
  physeq_prok_ctr_ev %>% 
  subset_samples(physeq = .,Time %in% c("0"))
otu_table(physeq_prok_ctr_start) <- 
  otu_table(physeq_prok_ctr_start)[which(rowSums(otu_table(physeq_prok_ctr_start)) > 0),]
physeq_prok_ctr_start

physeq_prok_ctr_start@sam_data

adonis_prok_ctr_start <-
  Adjadonis(physeq_prok_ctr_start, t(otu) ~ Treatment * Pore)
adonis_prok_ctr_start[[2]]
adonis_prok_ctr_start[[1]]$aov.tab

adonis_prok_ctr_start <-
  Adjadonis(physeq_prok_ctr_start, t(otu) ~ Pore * Treatment)
adonis_prok_ctr_start[[2]]
adonis_prok_ctr_start[[1]]$aov.tab

write.csv(adonis_prok_ctr_start[[1]]$aov.tab, "New_Manuscript_draft_2023/adonis-dry-contr-bact.csv")

# BETA-DISPERSION ---------------------------------------------------------------
BetadispExtr <- function(physeq, Var){
  otu <- as.data.frame(otu_table(physeq))
  metadata <- as.matrix(sample_data(physeq))
  metadata <- as.data.frame(metadata)
  disp <-
    betadisper(
      vegan::vegdist(t(otu), method="bray"),
      metadata[,Var])
  anova_d <-
    anova(disp,
          permutations = how(nperm=9999))
  p_adj <-
    round(p.adjust(anova_d$`Pr(>F)`,
                   "bonferroni"), 4)
  dist_var <-
    vegan::permutest(disp,
                     permutations = 9999,
                     pairwise = T)
  return(list(dist_var, p_adj, disp))
}

BetadispExtr(physeq_fungi_ctr_ev, "Treatment")[[1]]$tab
BetadispExtr(physeq_fungi_ctr_ev, "Pore")[[1]]$tab

BetadispExtr(physeq_prok_ctr_ev, "Treatment")[[1]]$tab
BetadispExtr(physeq_prok_ctr_ev, "Pore")[[1]]$tab

write.csv(BetadispExtr(physeq_fungi_ctr_ev, "Treatment")[[1]]$tab, "New_Manuscript_draft_2023/betadisper-dry-contr-fungi_Treat.csv")
write.csv(BetadispExtr(physeq_fungi_ctr_ev, "Pore")[[1]]$tab, "New_Manuscript_draft_2023/betadisper-dry-contr-fungi_Pore.csv")

write.csv(BetadispExtr(physeq_prok_ctr_ev, "Treatment")[[1]]$tab, "New_Manuscript_draft_2023/betadisper-dry-contr-bact_Treat.csv")
write.csv(BetadispExtr(physeq_prok_ctr_ev, "Pore")[[1]]$tab, "New_Manuscript_draft_2023/betadisper-dry-contr-bact_Pore.csv")

# ***********************************************-------------------------------
# > Data set without controls --------------------------------------------------
# PERMANOVA for Fungi ----------------------------------------------------------
adonis_fungi_inter <-
  Adjadonis(physeq_fungi_ev,
          t(otu) ~ Slice * Treatment * Plant * Pore * Moisture)

adonis_fungi_inter_2 <-
  Adjadonis(physeq_fungi_ev,
            t(otu) ~ Moisture * Pore * Plant * Treatment * Slice)

adonis_fungi_add <-
  Adjadonis(physeq_fungi_ev,
          t(otu) ~ Moisture + Pore + Plant + Treatment + Slice)


# adding time as a block/random effect in the model
adonis_fungi_inter_t <-
  Adjadonis(physeq_fungi_ev,
            t(otu) ~ Time + Slice * Treatment * Plant * Pore * Moisture)
write.csv(adonis_fungi_inter_t[[2]], "adonis_fungi_inter_t.csv")

adonis_fungi_inter_t[[1]]$aov.tab

adonis_fungi_inter_str <-
  Adjadonis(physeq_fungi_ev,
            t(otu) ~ Slice * Treatment * Plant * Pore * Moisture,
            strata = metadata$Time)


# PERMANOVA for Bacteria -------------------------------------------------------
adonis_prok_inter <-
  Adjadonis(physeq_prok_ev,
          t(otu) ~ Slice * Treatment * Plant * Pore * Moisture)

adonis_prok_inter_2 <-
  Adjadonis(physeq_prok_ev,
            t(otu) ~ Moisture * Pore * Plant * Treatment * Slice)

adonis_prok_add <-
  Adjadonis(physeq_prok_ev,
          t(otu) ~ Moisture + Pore + Plant + Treatment + Slice)

# adding time as a block/random effect in the model
adonis_prok_inter_t <-
  Adjadonis(physeq_prok_ev,
            t(otu) ~ Time + Slice * Treatment * Plant * Pore * Moisture)
write.csv(adonis_prok_inter_t[[2]], "adonis_prok_inter_t.csv")

adonis_prok_inter_t[[1]]$aov.tab

adonis_prok_inter_str <-
  Adjadonis(physeq_prok_ev,
            t(otu) ~ Slice * Treatment * Plant * Pore * Moisture,
            strata = metadata$Time)


# BETA-DISPERSION ---------------------------------------------------------------
BetadispExtr(physeq_fungi_ev, "Slice")

betadisp_fungi <-
  data.frame(
    rbind(
      BetadispExtr(physeq_fungi_ev, "Slice")[[1]]$tab,
      BetadispExtr(physeq_fungi_ev, "Treatment")[[1]]$tab,
      BetadispExtr(physeq_fungi_ev, "Plant")[[1]]$tab,
      BetadispExtr(physeq_fungi_ev, "Pore")[[1]]$tab,
      BetadispExtr(physeq_fungi_ev, "Moisture")[[1]]$tab),
    Padj = c(
      BetadispExtr(physeq_fungi_ev, "Slice")[[2]],
      BetadispExtr(physeq_fungi_ev, "Treatment")[[2]],
      BetadispExtr(physeq_fungi_ev, "Plant")[[2]],
      BetadispExtr(physeq_fungi_ev, "Pore")[[2]],
      BetadispExtr(physeq_fungi_ev, "Moisture")[[2]]),
    Variable = c(rep("Silce",2),rep("Treatment",2),rep("Plant",2),
                 rep("Pore",2), rep("Moisture", 2))
  )

betadisp_fungi
write.csv(betadisp_fungi, "betadisp_fungi.csv")

betadisp_bact <-
  data.frame(
    rbind(BetadispExtr(physeq_prok_ev, "Slice")[[1]]$tab,
        BetadispExtr(physeq_prok_ev, "Treatment")[[1]]$tab,
        BetadispExtr(physeq_prok_ev, "Plant")[[1]]$tab,
        BetadispExtr(physeq_prok_ev, "Pore")[[1]]$tab,
        BetadispExtr(physeq_prok_ev, "Moisture")[[1]]$tab,
        BetadispExtr(physeq_prok_ev, "Time")[[1]]$tab),
    Padj = c(
      BetadispExtr(physeq_prok_ev, "Slice")[[2]],
      BetadispExtr(physeq_prok_ev, "Treatment")[[2]],
      BetadispExtr(physeq_prok_ev, "Plant")[[2]],
      BetadispExtr(physeq_prok_ev, "Pore")[[2]],
      BetadispExtr(physeq_prok_ev, "Moisture")[[2]],
      BetadispExtr(physeq_prok_ev, "Time")[[2]]),
        Variable = c(rep("Silce",2),rep("Treatment",2),rep("Plant",2),
                     rep("Pore",2), rep("Moisture", 2), rep("Time", 2))
    )

write.csv(betadisp_bact, "betadisp_bact.csv")

# ***********************************************-------------------------------
# Split PERMANOVA --------------------------------------------------------------
# Data sets WITHOUT Controls ---------------------------------------------------
SplitAdonis <- function(physeq){
p1 <-
  subset_samples(physeq,
               Treatment %in% c("T1") &
                 Slice %in% c("Leaf")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  Adjadonis(physeq = .,
          t(otu) ~ Time + Plant * Pore * Moisture)
p2 <-
subset_samples(physeq,
               Treatment %in% c("T1") &
                 Slice %in% c("Nearsoil")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  Adjadonis(physeq = .,
            t(otu) ~ Time + Plant * Pore * Moisture)
p3 <-
subset_samples(physeq,
               Treatment %in% c("T1") &
                 Slice %in% c("Farsoil")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  Adjadonis(physeq = .,
            t(otu) ~ Time + Plant * Pore * Moisture)
p4 <-
subset_samples(physeq,
               Treatment %in% c("T4") &
                 Slice %in% c("Leaf")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  Adjadonis(physeq = .,
            t(otu) ~ Time + Plant * Pore * Moisture)
p5 <-
subset_samples(physeq,
               Treatment %in% c("T4") &
                 Slice %in% c("Nearsoil")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  Adjadonis(physeq = .,
            t(otu) ~ Time + Plant * Pore * Moisture)
p6 <-
subset_samples(physeq,
               Treatment %in% c("T4") &
                 Slice %in% c("Farsoil")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  Adjadonis(physeq = .,
            t(otu) ~ Time + Plant * Pore * Moisture)

p <- rbind(data.frame(p1[[2]], Model=rep("T1-Leaf", nrow(p1[[2]]))),
           data.frame(p2[[2]], Model=rep("T1-Nearsoil", nrow(p2[[2]]))),
           data.frame(p3[[2]], Model=rep("T1-Farsoil", nrow(p3[[2]]))),
           data.frame(p4[[2]], Model=rep("T4-Leaf", nrow(p4[[2]]))),
           data.frame(p5[[2]], Model=rep("T4-Nearsoil", nrow(p5[[2]]))),
           data.frame(p6[[2]], Model=rep("T4-Farsoil", nrow(p6[[2]]))))
return(p)

}

split_adonis_fungi <-
  SplitAdonis(physeq_fungi_ev)
split_adonis_fungi
write.csv(split_adonis_fungi, "split_adonis_fungi.csv")

split_adonis_bacteria <-
  SplitAdonis(physeq_prok_ev)
split_adonis_bacteria
write.csv(split_adonis_bacteria, "split_adonis_bacteria.csv")


# Data sets WITH Controls ------------------------------------------------------
SplitAdonisContr <- function(physeq){
  p1 <-
    subset_samples(physeq,
                   Treatment1 %in% c("T1", "Control") & Slice %in% c("Leaf", "Control")) %>% 
    subset_samples(physeq = ., Treatment %in% c("T1")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    Adjadonis(physeq = .,
              t(otu) ~ Time + Plant * Pore1 * Moisture1)
  p2 <-
    subset_samples(physeq,
                   Treatment1 %in% c("T1", "Control") & Slice %in% c("Nearsoil", "Control")) %>% 
    subset_samples(physeq = ., Treatment %in% c("T1")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    Adjadonis(physeq = .,
              t(otu) ~ Time + Plant * Pore1 * Moisture1)
  p3 <-
    subset_samples(physeq,
                   Treatment1 %in% c("T1", "Control") & Slice %in% c("Farsoil", "Control")) %>% 
    subset_samples(physeq = ., Treatment %in% c("T1")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    Adjadonis(physeq = .,
              t(otu) ~ Time + Plant * Pore1 * Moisture1)
  p4 <-
    subset_samples(physeq,
                   Treatment1 %in% c("Control"))  %>% 
    subset_samples(physeq = ., Treatment %in% c("T1")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    Adjadonis(physeq = .,
              t(otu) ~ Time + Pore * Moisture)
  p5 <-
    subset_samples(physeq,
                   Treatment1 %in% c("T4", "Control") & Slice %in% c("Leaf", "Control")) %>% 
    subset_samples(physeq = ., Treatment %in% c("T4")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    Adjadonis(physeq = .,
              t(otu) ~ Time + Plant * Pore1 * Moisture1)
  p6 <-
    subset_samples(physeq,
                   Treatment1 %in% c("T4", "Control") & Slice %in% c("Nearsoil", "Control")) %>% 
    subset_samples(physeq = ., Treatment %in% c("T4")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    Adjadonis(physeq = .,
              t(otu) ~ Time + Plant * Pore1 * Moisture1)
  p7 <-
    subset_samples(physeq,
                   Treatment1 %in% c("T4", "Control") & Slice %in% c("Farsoil", "Control")) %>% 
    subset_samples(physeq = ., Treatment %in% c("T4")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    Adjadonis(physeq = .,
              t(otu) ~ Time + Plant * Pore1 * Moisture1)
  p8 <-
    subset_samples(physeq,
    Treatment1 %in% c("Control"))  %>% 
    subset_samples(physeq = ., Treatment %in% c("T4")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    Adjadonis(physeq = .,
              t(otu) ~ Time + Pore * Moisture)
                   
  p <- rbind(data.frame(p1[[2]], Model=rep("T1-Leaf", nrow(p1[[2]]))),
             data.frame(p2[[2]], Model=rep("T1-Nearsoil", nrow(p2[[2]]))),
             data.frame(p3[[2]], Model=rep("T1-Farsoil", nrow(p3[[2]]))),
             data.frame(p4[[2]], Model=rep("T1-Control", nrow(p4[[2]]))),
             data.frame(p5[[2]], Model=rep("T4-Leaf", nrow(p5[[2]]))),
             data.frame(p6[[2]], Model=rep("T4-Nearsoil", nrow(p6[[2]]))),
             data.frame(p7[[2]], Model=rep("T4-Farsoil", nrow(p7[[2]]))),
             data.frame(p8[[2]], Model=rep("T4-Control", nrow(p8[[2]]))))
  return(p)
}



split_adonis_fungi_ctr <-SplitAdonisContr(physeq_fungi_ctr_ev)
write.csv(split_adonis_fungi_ctr, "New_Manuscript_draft_2023/split_adonis_fungi_ctr.csv")

split_adonis_bacteria_ctr <-SplitAdonisContr(physeq_prok_ctr_ev)
write.csv(split_adonis_bacteria_ctr, "New_Manuscript_draft_2023/split_adonis_bacteria_ctr.csv")


# Splitted BETADISPER ----------------------------------------------------------
# Data sets with CONTROLS ------------------------------------------------------
SplitBetadisper <- function(physeq, Var1, Var2=NULL, Var3=NULL, Var4=NULL){
  if (missing(Var2) & missing(Var3) & missing(Var4)){
    b1 <- BetadispExtr(physeq, Var1)
    b <-
      data.frame(b1[[1]]$tab, Var=rep(Var1, nrow(b1[[1]]$tab)))
    b <- data.frame(b, Padj=b1[[2]])
  }else if (missing(Var3) & missing(Var4)) {
    b1 <- BetadispExtr(physeq, Var1)
    b2 <- BetadispExtr(physeq, Var2)
    b <- rbind(
      data.frame(b1[[1]]$tab, Var=rep(Var1, nrow(b1[[1]]$tab))),
      data.frame(b2[[1]]$tab, Var=rep(Var2, nrow(b2[[1]]$tab))))
    b <- data.frame(b, Padj=c(b1[[2]], b2[[2]]))
  } else if (missing(Var4)){
    b1 <- BetadispExtr(physeq, Var1)
    b2 <- BetadispExtr(physeq, Var2)
    b3 <- BetadispExtr(physeq, Var3)
    b <- rbind(
    data.frame(b1[[1]]$tab, Var=rep(Var1, nrow(b1[[1]]$tab))),
    data.frame(b2[[1]]$tab, Var=rep(Var2, nrow(b2[[1]]$tab))),
    data.frame(b3[[1]]$tab, Var=rep(Var3, nrow(b3[[1]]$tab))))
    b <- data.frame(b, Padj=c(b1[[2]], b2[[2]], b3[[2]]))
  } else {
    b1 <- BetadispExtr(physeq, Var1)
    b2 <- BetadispExtr(physeq, Var2)
    b3 <- BetadispExtr(physeq, Var3)
    b4 <- BetadispExtr(physeq, Var4)
    b <- rbind(
      data.frame(b1[[1]]$tab, Var=rep(Var1, nrow(b1[[1]]$tab))),
      data.frame(b2[[1]]$tab, Var=rep(Var2, nrow(b2[[1]]$tab))),
      data.frame(b3[[1]]$tab, Var=rep(Var3, nrow(b3[[1]]$tab))),
      data.frame(b4[[1]]$tab, Var=rep(Var4, nrow(b4[[1]]$tab))))
    b <- data.frame(b, Padj=c(b1[[2]], b2[[2]], b3[[2]], b4[[2]]))
    }
    return(b)
}


subset_samples(physeq_fungi_ev,
               Treatment %in% c("T1") &
                 Slice %in% c("Leaf")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  SplitBetadisper(physeq = ., "Plant", "Pore") #T1 Leaf

subset_samples(physeq_fungi_ev,
               Treatment %in% c("T1") &
                 Slice %in% c("Nearsoil")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  SplitBetadisper(physeq = ., "Plant", "Pore", "Moisture") #T1 Nearsoil

subset_samples(physeq_fungi_ctr_ev,
               Treatment1 %in% c("Control"))  %>% 
  subset_samples(physeq = ., Treatment %in% c("T1")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  SplitBetadisper(physeq = ., "Pore", "Moisture")
  

# Data sets WITH Controls ------------------------------------------------------
ResultBetadisperSplitFungi <- function(){
  p1 <-
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment1 %in% c("T1", "Control") & Slice %in% c("Leaf", "Control")) %>% 
    subset_samples(physeq = ., Treatment %in% c("T1")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Plant", "Pore") #T1 Leaf
  p2 <-
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment1 %in% c("T1", "Control") & Slice %in% c("Nearsoil", "Control")) %>% 
    subset_samples(physeq = ., Treatment %in% c("T1")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Plant", "Pore") #T1 Nearsoil
  p3 <-
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment1 %in% c("T1", "Control") & Slice %in% c("Farsoil", "Control")) %>% 
    subset_samples(physeq = ., Treatment %in% c("T1")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Pore", "Moisture") #T1 Farsoil
  p4 <-
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment1 %in% c("Control"))  %>% 
    subset_samples(physeq = ., Treatment %in% c("T1")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Pore") #T1 Controls
  
  p5 <-
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment1 %in% c("T4", "Control") & Slice %in% c("Leaf", "Control")) %>% 
    subset_samples(physeq = ., Treatment %in% c("T4")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = .,"Time", "Plant", "Pore", "Moisture") #T4 Leaf
  p6 <-
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment1 %in% c("T4", "Control") & Slice %in% c("Nearsoil", "Control")) %>% 
    subset_samples(physeq = ., Treatment %in% c("T4")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Plant", "Pore", "Moisture") #T4 Nearsoil
  p7 <-
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment1 %in% c("T4", "Control") & Slice %in% c("Farsoil", "Control")) %>% 
    subset_samples(physeq = ., Treatment %in% c("T4")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Pore") #T4 Farsoil
  p8 <-
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment1 %in% c("Control"))  %>% 
    subset_samples(physeq = ., Treatment %in% c("T4")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Pore") #T4 Controls
  
  p <- rbind(
    data.frame(p1, Model=rep("T1-Leaf", nrow(p1))),
    data.frame(p2, Model=rep("T1-Nearsoil", nrow(p2))),
    data.frame(p3, Model=rep("T1-Farsoil", nrow(p3))),
    data.frame(p4, Model=rep("T1-Control", nrow(p4))),
             
    data.frame(p5, Model=rep("T4-Leaf", nrow(p5))),
    data.frame(p6, Model=rep("T4-Nearsoil", nrow(p6))),
    data.frame(p7, Model=rep("T4-Farsoil", nrow(p7))),
    data.frame(p8, Model=rep("T4-Control", nrow(p8)))
  )
  return(p)
}

split_betadisper_fungi_ctr <- ResultBetadisperSplitFungi()
write.csv(split_betadisper_fungi_ctr, "New_Manuscript_draft_2023/split_betadisper_fungi_ctr.csv")


ResultBetadisperSplitBact <- function(){
  p1 <-
    subset_samples(physeq_prok_ctr_ev,
                   Treatment1 %in% c("T1", "Control") & Slice %in% c("Leaf", "Control")) %>% 
    subset_samples(physeq = ., Treatment %in% c("T1")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = .,"Time", "Plant", "Pore", "Moisture") #T1 Leaf
  p2 <-
    subset_samples(physeq_prok_ctr_ev,
                   Treatment1 %in% c("T1", "Control") & Slice %in% c("Nearsoil", "Control")) %>% 
    subset_samples(physeq = ., Treatment %in% c("T1")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Time", "Plant", "Pore", "Moisture") #T1 Nearsoil
  p3 <-
    subset_samples(physeq_prok_ctr_ev,
                   Treatment1 %in% c("T1", "Control") & Slice %in% c("Farsoil", "Control")) %>% 
    subset_samples(physeq = ., Treatment %in% c("T1")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Time","Plant","Pore") #T1 Farsoil
  p4 <-
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment1 %in% c("Control"))  %>% 
    subset_samples(physeq = ., Treatment %in% c("T1")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = .,"Time", "Pore") #T1 Controls
  
  p5 <-
    subset_samples(physeq_prok_ctr_ev,
                   Treatment1 %in% c("T4", "Control") & Slice %in% c("Leaf", "Control")) %>% 
    subset_samples(physeq = ., Treatment %in% c("T4")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Time", "Plant", "Pore", "Moisture") #T4 Leaf
  p6 <-
    subset_samples(physeq_prok_ctr_ev,
                   Treatment1 %in% c("T4", "Control") & Slice %in% c("Nearsoil", "Control")) %>% 
    subset_samples(physeq = ., Treatment %in% c("T4")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Time", "Plant", "Pore", "Moisture") #T4 Nearsoil
  p7 <-
    subset_samples(physeq_prok_ctr_ev,
                   Treatment1 %in% c("T4", "Control") & Slice %in% c("Farsoil", "Control")) %>% 
    subset_samples(physeq = ., Treatment %in% c("T4")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Time", "Plant", "Pore", "Moisture") #T4 Farsoil
  p8 <-
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment1 %in% c("Control"))  %>% 
    subset_samples(physeq = ., Treatment %in% c("T4")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Time", "Pore", "Moisture") #T4 Controls
  
  px <- rbind(p1, p2, p3, p4, p5, p6, p7, p8)
  py <- rbind(data.frame(p1, Model=rep("T1-Leaf", nrow(p1))),
             data.frame(p2, Model=rep("T1-Nearsoil", nrow(p2))),
             data.frame(p3, Model=rep("T1-Farsoil", nrow(p3))),
             data.frame(p4, Model=rep("T1-Control", nrow(p4))),
             data.frame(p5, Model=rep("T4-Leaf", nrow(p5))),
             data.frame(p6, Model=rep("T4-Nearsoil", nrow(p6))),
             data.frame(p7, Model=rep("T4-Farsoil", nrow(p7))),
             data.frame(p8, Model=rep("T4-Control", nrow(p8)))
             )
  return(py)
}


split_betadisper_bact_ctr <- ResultBetadisperSplitBact()
write.csv(split_betadisper_bact_ctr, "New_Manuscript_draft_2023/split_betadisper_bact_ctr.csv")

# Data sets WITHOUT Controls ---------------------------------------------------
ResultBetadisperSplitFungi <- function(){
  p1 <-
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment %in% c("T1") &
                     Slice %in% c("Leaf")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Plant", "Pore") #T1 Leaf
  p2 <-
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment %in% c("T1") &
                     Slice %in% c("Nearsoil")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Plant", "Pore", "Moisture") #T1 Nearsoil
  p3 <-
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment %in% c("T1") &
                     Slice %in% c("Farsoil")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Pore", "Moisture") #T1 Farsoil
  p4 <-
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment1 %in% c("Control"))  %>% 
    subset_samples(physeq = ., Treatment %in% c("T1")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Pore", "Moisture") #T1 Controls
  
  p5 <-
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment %in% c("T4") &
                     Slice %in% c("Leaf")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Plant", "Pore", "Moisture") #T4 Leaf
  p6 <-
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment %in% c("T4") &
                     Slice %in% c("Nearsoil")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Plant", "Pore", "Moisture") #T4 Nearsoil
  p7 <-
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment %in% c("T4") &
                     Slice %in% c("Farsoil")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Pore") #T4 Farsoil
  p4 <-
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment1 %in% c("Control"))  %>% 
    subset_samples(physeq = ., Treatment %in% c("T4")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Pore", "Moisture") #T4 Controls
  
  p <- rbind(data.frame(p1, Model=rep("T1-Leaf", nrow(p1))),
             data.frame(p2, Model=rep("T1-Nearsoil", nrow(p2))),
             data.frame(p3, Model=rep("T1-Farsoil", nrow(p3))),
             data.frame(p4, Model=rep("T4-Control", nrow(p4))),
             data.frame(p5, Model=rep("T4-Leaf", nrow(p5))),
             data.frame(p6, Model=rep("T4-Nearsoil", nrow(p6))),
             data.frame(p7, Model=rep("T4-Farsoil", nrow(p7))),
             data.frame(p8, Model=rep("T4-Control", nrow(p8))),
  )
  return(p)
}

split_betadisper_fungi <-
  ResultBetadisperSplitFungi()
write.csv(split_betadisper_fungi, "split_betadisper_fungi.csv")


ResultBetadisperSplitBact <- function(){
  p1 <-
    subset_samples(physeq_prok_ev,
                   Treatment %in% c("T1") &
                     Slice %in% c("Leaf")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = .,"Time", "Pore", "Moisture", "Plant") #T1 Leaf
  p2 <-
    subset_samples(physeq_prok_ev,
                   Treatment %in% c("T1") &
                     Slice %in% c("Nearsoil")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Time","Pore", "Moisture", "Plant") #T1 Nearsoil
  p3 <-
    subset_samples(physeq_prok_ev,
                   Treatment %in% c("T1") &
                     Slice %in% c("Farsoil")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Time","Pore") #T1 Farsoil
  p4 <-
    subset_samples(physeq_prok_ev,
                   Treatment %in% c("T4") &
                     Slice %in% c("Leaf")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Time","Pore", "Moisture", "Plant") #T4 Leaf
  p5 <-
    subset_samples(physeq_prok_ev,
                   Treatment %in% c("T4") &
                     Slice %in% c("Nearsoil")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Time", "Pore", "Moisture", "Plant") #T4 Nearsoil
  p6 <-
    subset_samples(physeq_prok_ev,
                   Treatment %in% c("T4") &
                     Slice %in% c("Farsoil")) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    SplitBetadisper(physeq = ., "Time","Pore", "Moisture") #T4 Farsoil
  px <- rbind(p1, p2, p3, p4,p5, p6)
  py <- rbind(data.frame(p1, Model=rep("T1-Leaf", nrow(p1))),
              data.frame(p2, Model=rep("T1-Nearsoil", nrow(p2))),
              data.frame(p3, Model=rep("T1-Farsoil", nrow(p3))),
              data.frame(p4, Model=rep("T4-Leaf", nrow(p4))),
              data.frame(p5, Model=rep("T4-Nearsoil", nrow(p5))),
              data.frame(p6, Model=rep("T4-Farsoil", nrow(p6))))
  return(py)
}

subset_samples(physeq_prok_ev,
               Treatment %in% c("T4") &
                 Slice %in% c("Farsoil")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  SplitBetadisper(physeq = ., "Time","Plant", "Pore", "Moisture") #T1 Nearsoil


split_betadisper_bact <-
  ResultBetadisperSplitBact()
write.csv(split_betadisper_bact, "split_betadisper_bact.csv")

# ***********************************************---------------------------------
# PCOA ---------------------------------------------------------------------------
# The most important factors are Slice and Treatment. However Pore size, Plant,
# and Moisture are also significant.

# Reorder factors
sample_data(physeq_fungi_ev)$Slice <-
  factor(sample_data(physeq_fungi_ev)$Slice, levels = c('Leaf', 'Nearsoil', 'Farsoil'))

sample_data(physeq_prok_ev)$Slice <-
  factor(sample_data(physeq_prok_ev)$Slice, levels = c('Leaf', 'Nearsoil', 'Farsoil'))

pcoa <- ordinate(physeq_fungi_ev, method ="PCoA", distance="bray")
str(pcoa)


# Plot PCoA
PlotOrdin <-function(dataframe){
  pcoa <- ordinate(dataframe, method ="PCoA", distance="bray")
  plot_ord = plot_ordination(dataframe, pcoa,
                             type="sites",
                             shape="Treatment") +
    #geom_point(size=2, aes(colour=Slice), show.legend = FALSE) +
    geom_point(size=2, aes(fill=Slice)) +
    scale_shape_manual(values = c(21, 24)) +
    scale_size_manual(values=c(1, 2.5))+
    #scale_color_manual(values=paletteCB5_tris, 
    #                    guide = guide_legend(show = FALSE)) +
    scale_fill_manual(values=paletteCB5_tris) +
    theme_bw() +
    theme(plot.caption = element_text(size = 10, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 9, face = "bold"),
          legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.3, "cm"),
          #legend.direction = "vertical", legend.spacing.y = unit(0.0001, "cm"),
          legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8),
          legend.position = ("bottom")) +
    guides(fill = guide_legend(order = 1, nrow = 1, title = "Niche",
                                override.aes = list(shape = 22, size=3)),
           shape = guide_legend(order = 3, nrow=1, 
                                override.aes = list(color = "black")))
  return(plot_ord)
}


PlotOrdin(physeq_fungi_ev) +
  stat_ellipse(aes(group=Slice), type="t", alpha=0.8, linetype = 2, show.legend = FALSE)

PlotOrdin(physeq_prok_ev, Var1 = "Slice", Var2 = "Treatment") +
  stat_ellipse(aes(group=Slice), type="t", alpha=0.8, linetype = 2, show.legend = FALSE)


sample_data(physeq_fungi_ctr_ev)$Slice <-
  factor(sample_data(physeq_fungi_ctr_ev)$Slice, levels = c('Leaf', 'Nearsoil', 'Farsoil', "Control", "Start"))
sample_data(physeq_fungi_ctr_ev)$Treatment1 <-
  factor(sample_data(physeq_fungi_ctr_ev)$Treatment1, levels = c('T1', 'T4', 'Control', "Start"))

sample_data(physeq_prok_ctr_ev)$Slice <-
  factor(sample_data(physeq_prok_ctr_ev)$Slice, levels = c('Leaf', 'Nearsoil', 'Farsoil', "Control", "Start"))
sample_data(physeq_prok_ctr_ev)$Treatment1 <-
  factor(sample_data(physeq_prok_ctr_ev)$Treatment1, levels = c('T1', 'T4', 'Control', "Start"))


PlotOrdin(physeq_fungi_ctr_ev, Var1 = "Slice", Var2 = "Treatment") +
  stat_ellipse(aes(group=Slice), type="t", alpha=0.8, linetype = 2, show.legend = FALSE)

PlotOrdin(physeq_fungi_ctr_ev, Var1 = "Slice", Var2 = "Treatment") 

# **** FIGURE 1 - main ordination ----------------------------------------------
Fig1_pcoa <-
  ggarrange(
    PlotOrdin(physeq_fungi_ctr_ev) +
      labs(title="Fungi"),
    PlotOrdin(physeq_prok_ctr_ev) +
    labs(title="Bacteria"),
    labels = c("A", "B"),
    nrow = 1,
    ncol = 2, 
    common.legend = TRUE,
    legend = "bottom")
 
Fig1_pcoa

title1 <- text_grob("Microbiome Compostion", size = 12, face = 2)
grid.arrange(Fig1_pcoa, top=title1)


# **** FIGURE S1 - main ordination WITHOUT Controls ----------------------------
FigS1_pcoa <-
  ggarrange(
    PlotOrdin(physeq_fungi_ev) +
      labs(title="Fungi"),
    PlotOrdin(physeq_prok_ev) +
      labs(title="Bacteria"),
    labels = c("A", "B"),
    nrow = 1,
    ncol = 2, 
    common.legend = TRUE,
    legend = "bottom")

FigS1_pcoa

title1 <- text_grob("Microbiome Compostion Without Control Samples", size = 12, face = 2)
grid.arrange(FigS1_pcoa, top=title1)



# plotting main ordination
Fig_1_A_pcoa <-
  ggarrange(
    PlotOrdin(physeq_fungi_ev, Var1 = "Slice", Var2 = "Treatment") +
      stat_ellipse(aes(group=Slice), type="t", alpha=0.8, linetype = 2, show.legend = FALSE) +
      annotate("text", -Inf, Inf, label = "Slice", color = "black", fontface = 2, size = 3, hjust = -0.4, vjust = 2)  +
      labs(title="Fungi", ),
    labels = c("A"),
    nrow = 1,
    ncol = 1)

Fig_1_A_pcoa

Fig_2_A_pcoa <-
  ggarrange(
    PlotOrdin(physeq_prok_ev, Var1 = "Slice", Var2 = "Treatment") +
      stat_ellipse(aes(group=Slice), type="t", alpha=0.8, linetype = 2, show.legend = FALSE) +
      annotate("text", -Inf, Inf, label = "Slice", color = "black", fontface = 2, size = 3, hjust = -0.4, vjust = 2)  +
      labs(title="Bacteria", ),
    labels = c("A"),
    nrow = 1,
    ncol = 1)

Fig_2_A_pcoa


# position of the annotation layer ---------------------------------------------
# Annotation position and hjust, vjust use to modify it
# xpos ypos        annotateText        hjust    vjust
# 1 -Inf -Inf Bottom Left (h0,v0)        0        0
# 2 -Inf  Inf    Top Left (h0,v1)        0        1
# 3  Inf -Inf  Bottom Right h1,v0        1        0
# 4  Inf  Inf     Top Right h1,v1        1        1

PlotOrdinSplit <-function(dataframe){
  pcoa <- ordinate(dataframe, method ="PCoA", distance="bray")
  perc <-
    c(round(pcoa$values$Relative_eig[1:2] * 100, digits =  1))
  names(perc) <- c("fist", "second")
  perc <- as.data.frame(t(as.data.frame(perc)))
  print(perc)
  
  df <-
  left_join(
    pcoa$vectors %>% 
      as.data.frame() %>% 
      select(Axis.1, Axis.2) %>% 
      rownames_to_column("SampleID"),
    as.data.frame(as.matrix(
      dataframe@sam_data)) %>% 
      rownames_to_column("SampleID"),
    by="SampleID") %>% 
    mutate(Moisture = recode(Moisture, CBP = "Low", FC = "High"))
  
  df$Plant <-
    factor(df$Plant, levels = c('Corn', 'Soy', 'Control'))
  #df$Moisture <-
  #  factor(df$Moisture, levels = c('CBP', 'FC'))
  df$Moisture <-
    factor(df$Moisture, levels = c('Low', 'High'))
  df$Pore <-
    factor(df$Pore, levels = c("Small", 'Large'))
  
  plot_ord <-
    ggplot(df, aes(x=Axis.1, y=Axis.2, 
                   shape = Moisture, 
                   size = Pore)) +
    #geom_point(aes(color=Plant)) +
    geom_point(aes(fill = Plant)) +
    scale_shape_manual(values = c(25, 24)) +
    scale_size_manual(breaks = c("Small", "Large"),
                      values = c(1, 2.5)) +
    scale_fill_manual(values = paletteCB3) +
    #scale_color_manual(values = c("black", "black","black")) +
    theme_bw() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
            axis.title = element_text(angle = 0, size = 9, face = "bold"),
            legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.3, "cm"),
            legend.title = element_text(size = 9, face = "bold"), 
            legend.text = element_text(size = 8),
            legend.position = ("bottom")) +
    guides(fill = guide_legend(order = 1, nrow = 1, title = "Niche",
                  override.aes = list(shape = 22, size=3),
          shape = guide_legend(order = 3, nrow=1,
                  override.aes = list(color = "black"))))
  
  return(plot_ord)
}


subset_samples(physeq_fungi_ctr_ev,
               Treatment1 %in% c("T1", "Control") & Slice %in% c("Leaf", "Control")) %>% 
  subset_samples(physeq = ., Treatment %in% c("T1")) %>% 
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  PlotOrdinSplit(dataframe = .) +
  stat_ellipse(aes(x=Axis.1, y=Axis.2, group=Plant), inherit.aes = FALSE, type="norm", level=0.75, 
               alpha=0.8, color = "black", linetype = 2, show.legend = FALSE) +
  annotate("text", -Inf, Inf, label = expression(paste(bold("Plant "), italic(R)^2,"= 29.5%")), 
           color = "black",size = 3, hjust = -0.1, vjust = 1.2)  +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  labs(title="T1 - Leaf") 




#labs(x=expression(paste("Axis.1 [ ", perc[1], "% ]"), parse=TRUE),
#      y=expression(paste("Axis.2 [ ", perc[1], "% ]"), parse=TRUE))
#labs(x="Axis.1 [24.8%]", y= "Axis.2 [14.2%]")

# plot PCOA split ordination fungi ------------------------------------------------------------
library(ggforce)

control_T1 <-
subset_samples(physeq_fungi_ctr_ev,
               Treatment1 %in% c("Control")) %>% 
  subset_samples(physeq = ., Treatment %in% c("T1")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  PlotOrdinSplit(dataframe = .) + 
  scale_fill_manual(values = c("black", "black")) +
  stat_ellipse(aes(x=Axis.1, y=Axis.2, group=Pore), inherit.aes = FALSE,  type="norm", level=0.75, 
               alpha=0.8, color = "black", linetype = 2, show.legend = FALSE) +
  annotate("text", Inf, -Inf, label = expression(paste(bold("Pore "), italic(R)^2,"= 10.4%"),parse=TRUE), 
            color = "black",  size = 3, hjust = 1.1, vjust = -0.3) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  guides(fill = guide_legend(order = 1, nrow = 1, title = "Niche",
                             override.aes = list(shape = 15, size=2.5,color="black"),
                             shape = guide_legend(order = 3, nrow=1))) +
  labs(title="T1 - Control") +
  labs(x="Axis.1 [14.1%]", y= "Axis.2 [13.1%]")

control_T1


control_T4 <-
subset_samples(physeq_fungi_ctr_ev,
               Treatment1 %in% c("Control")) %>% 
  subset_samples(physeq = ., Treatment %in% c("T4")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  PlotOrdinSplit(dataframe = .) +
  scale_color_manual(values = "black") +
  scale_fill_manual(values = "black") +
  stat_ellipse(aes(x=Axis.1, y=Axis.2, group=Pore), inherit.aes = FALSE, type="norm", level=0.75, 
               alpha=0.8, color = "black", linetype = 2, show.legend = FALSE) +
  annotate("text", Inf, -Inf, label = expression(paste(bold("Pore "), italic(R)^2,"= 10.1%")), 
           color = "black", size = 3, hjust = 1.1, vjust = -0.3) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  guides(fill = guide_legend(order = 1, nrow = 1, title = "Niche",
                             override.aes = list(shape = 15, size=2.5,color="black"),
                             shape = guide_legend(order = 3, nrow=1))) +
  labs(title="T4 - Control") +
  labs(x="Axis.1 [20.2%]", y= "Axis.2 [18.3%]")

control_T4


pcoa_split_fungi <-
  ggarrange(
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment1 %in% c("T1", "Control") & Slice %in% c("Leaf", "Control")) %>% 
      subset_samples(physeq = ., Treatment %in% c("T1")) %>% 
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    PlotOrdinSplit(dataframe = .) +
    stat_ellipse(aes(group=Plant), type="norm", level=0.75, alpha=0.8, color = "black",linetype = 2, show.legend = FALSE) +
    annotate("text", -Inf, Inf, label = expression(paste(bold("Plant "), italic(R)^2,"= 29.5%")), 
             color = "black",size = 3, hjust = -0.05,  vjust = 1.2)  +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+
    labs(title="T1 - Leaf") + 
      labs(x="Axis.1 [24.8%]", y= "Axis.2 [14.2%]"),
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment1 %in% c("T1", "Control") & Slice %in% c("Nearsoil", "Control")) %>% 
      subset_samples(physeq = ., Treatment %in% c("T1")) %>% 
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    PlotOrdinSplit(dataframe = .) +
    stat_ellipse(aes(group=Plant), type="norm", level=0.75, alpha=0.8, color = "black",linetype = 2, show.legend = FALSE) +
    annotate("text", -Inf, Inf, label = expression(paste(bold("Plant "), italic(R)^2,"= 24.3%")),
             color = "black", size = 3, hjust = -0.05,  vjust = 1.2) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    labs(title="T1 - Nearsoil")+
      labs(x="Axis.1 [25.8%]", y= "Axis.2 [17%]"),
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment1 %in% c("T1", "Control") & Slice %in% c("Farsoil", "Control")) %>% 
      subset_samples(physeq = ., Treatment %in% c("T1")) %>% 
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    PlotOrdinSplit(dataframe = .) +
    stat_ellipse(aes(x=Axis.1, y=Axis.2, group=Pore), inherit.aes = FALSE,type="norm", level=0.75, alpha=0.8, color = "black",linetype = 2, show.legend = FALSE) +
    annotate("text", Inf, Inf, label = expression(paste(bold("Pore "), italic(R)^2,"= 7.3%")), 
             color = "black", size = 3, hjust = 1.1, vjust = 1.1) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    labs(title="T1 - Farsoil") +
      labs(x="Axis.1 [19%]", y= "Axis.2 [10.4%]"),
    control_T1,
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment1 %in% c("T4", "Control") & Slice %in% c("Leaf", "Control")) %>% 
      subset_samples(physeq = ., Treatment %in% c("T4")) %>% 
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    PlotOrdinSplit(dataframe = .) +
    stat_ellipse(aes(group=Plant), type="norm", level=0.75, alpha=0.8, color = "black",linetype = 2, show.legend = FALSE) +
    annotate("text", -Inf, Inf, label = expression(paste(bold("Plant "), italic(R)^2,"= 28.7%")), 
             color = "black", size = 3, hjust = -0.05,  vjust = 1.2) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    labs(title="T4 - Leaf") +
      labs(x="Axis.1 [24.3%]", y= "Axis.2 [15.1%]"),
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment1 %in% c("T4", "Control") & Slice %in% c("Nearsoil", "Control")) %>% 
      subset_samples(physeq = ., Treatment %in% c("T4")) %>% 
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    PlotOrdinSplit(dataframe = .) +
    stat_ellipse(aes(group=Plant), type="norm", level=0.75, alpha=0.8,color = "black", linetype = 2, show.legend = FALSE) +
    annotate("text", -Inf, Inf, label = expression(paste(bold("Plant "), italic(R)^2,"= 15.1%")),
             color = "black", size = 3, hjust = -0.05, vjust = 1.2) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    labs(title="T4 - Nearsoil")+
      labs(x="Axis.1 [15.6%]", y= "Axis.2 [14.4%]"),
    subset_samples(physeq_fungi_ctr_ev,
                   Treatment1 %in% c("T4", "Control") & Slice %in% c("Farsoil", "Control")) %>% 
      subset_samples(physeq = ., Treatment %in% c("T4")) %>% 
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    PlotOrdinSplit(dataframe = .) +
      stat_ellipse(aes(x=Axis.1, y=Axis.2, group=Pore), inherit.aes = FALSE, type="norm", level=0.75, alpha=0.8, color = "black",linetype = 2, show.legend = FALSE) +
    annotate("text", Inf, Inf,label = expression(paste(bold("Pore "), italic(R)^2,"= 4.2%")),
             color = "black", fontface = 2, size = 3, hjust = 1.1, vjust = 1.1) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    labs(title="T4 - Farsoil")+
      labs(x="Axis.1 [18.9%]", y= "Axis.2 [11.5%]"),
    control_T4,
  labels = c("A","","","","B","","",""),
  ncol = 4,
  nrow = 2,
  legend = "bottom",
  common.legend=TRUE)

pcoa_split_fungi

library(grid)
library(gridExtra)

title1=text_grob("Fungal Microbiome Composition",size=12, face=2)
grid.arrange(pcoa_split_fungi, top = title1)


# plot PCOA split ordination bacteria ------------------------------------------------------------
control_T1_prok <-
  subset_samples(physeq_prok_ctr_ev,
                 Treatment1 %in% c("Control")) %>% 
  subset_samples(physeq = ., Treatment %in% c("T1")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  PlotOrdinSplit(dataframe = .) + 
  scale_fill_manual(values = c("black", "black")) +
  stat_ellipse(aes(x=Axis.1, y=Axis.2, group=Pore), inherit.aes = FALSE,  type="norm", level=0.75, 
               alpha=0.8, color = "black", linetype = 2, show.legend = FALSE) +
  annotate("text", Inf, -Inf, label = expression(paste(bold("Pore "), italic(R)^2,"= 10.4%"),parse=TRUE), 
           color = "black",  size = 3, hjust = 1.1, vjust = -0.3) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  labs(title="T1 - Control") +
  labs(x="Axis.1 [16.9%]", y= "Axis.2 [13.9%]")

control_T1_prok


control_T4_prok <-
  subset_samples(physeq_prok_ctr_ev,
                 Treatment1 %in% c("Control")) %>% 
  subset_samples(physeq = ., Treatment %in% c("T4")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  PlotOrdinSplit(dataframe = .) +
  scale_color_manual(values = "black") +
  scale_fill_manual(values = "black") +
  stat_ellipse(aes(x=Axis.1, y=Axis.2, group=Pore), inherit.aes = FALSE, type="norm", level=0.75, 
               alpha=0.8, color = "black", linetype = 2, show.legend = FALSE) +
  annotate("text", Inf, -Inf, label = expression(paste(bold("Pore "), italic(R)^2,"= 10.1%")), 
           color = "black", size = 3, hjust = 1.1, vjust = -0.3) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  labs(title="T4 - Control") +
  labs(x="Axis.1 [10.7%]", y= "Axis.2 [9.4%]")

control_T4_prok


# main plot
pcoa_split_prok <-
  ggarrange(
    subset_samples(physeq_prok_ctr_ev,
                   Treatment1 %in% c("T1", "Control") & Slice %in% c("Leaf", "Control")) %>% 
      subset_samples(physeq = ., Treatment %in% c("T1")) %>% 
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
      PlotOrdinSplit(dataframe = .) +
      stat_ellipse(aes(x=Axis.1, y=Axis.2, group=Plant), inherit.aes = FALSE, type="norm", 
                   level=0.75, alpha=0.8, color = "black",linetype = 2, show.legend = FALSE) +
      annotate("text", -Inf, -Inf,label = expression(paste(bold("Plant "), italic(R)^2,"= 39.5%")),
               color = "black", size = 3, hjust = -0.05, vjust = -0.3)  +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
      labs(title="T1 - Leaf")+
      labs(x="Axis.1 [37.5%]", y= "Axis.2 [12.7%]"),
    subset_samples(physeq_prok_ctr_ev,
                   Treatment1 %in% c("T1", "Control") & Slice %in% c("Nearsoil", "Control")) %>% 
      subset_samples(physeq = ., Treatment %in% c("T1")) %>% 
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
      PlotOrdinSplit(dataframe = .) +
      stat_ellipse(aes(x=Axis.1, y=Axis.2, group=Plant), inherit.aes = FALSE, type="norm", 
                   level=0.75, alpha=0.8, color = "black",linetype = 2, show.legend = FALSE) +
      annotate("text", -Inf, -Inf, label = expression(paste(bold("Plant "), italic(R)^2,"= 12.5%")),
               color = "black", size = 3, hjust = -0.05, vjust = -0.3) +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
      labs(title="T1 - Nearsoil")+
      labs(x="Axis.1 [11.4%]", y= "Axis.2 [9.7%]"),
    subset_samples(physeq_prok_ctr_ev,
                   Treatment1 %in% c("T1", "Control") & Slice %in% c("Farsoil", "Control")) %>% 
      subset_samples(physeq = ., Treatment %in% c("T1")) %>% 
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
      PlotOrdinSplit(dataframe = .) +
      stat_ellipse(aes(x=Axis.1, y=Axis.2, group=Plant), inherit.aes = FALSE, type="norm", 
                   level=0.75, alpha=0.8, color = "black",linetype = 2, show.legend = FALSE) +
      annotate("text", -Inf, -Inf, label = expression(paste(bold("Plant "), italic(R)^2,"= 6.5%")),
               color = "black",  size = 3, hjust = -0.05, vjust = -0.3) +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
      labs(title="T1 - Farsoil")+
      labs(x="Axis.1 [11.9%]", y= "Axis.2 [7.2%]"),
    control_T1,
    subset_samples(physeq_prok_ctr_ev,
                   Treatment1 %in% c("T4", "Control") & Slice %in% c("Leaf", "Control")) %>% 
      subset_samples(physeq = ., Treatment %in% c("T4")) %>% 
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
      PlotOrdinSplit(dataframe = .) +
      stat_ellipse(aes(x=Axis.1, y=Axis.2, group=Plant), inherit.aes = FALSE, type="norm", 
                   level=0.75, alpha=0.8, color = "black",linetype = 2, show.legend = FALSE) +
      annotate("text", -Inf, -Inf, label = expression(paste(bold("Plant "), italic(R)^2,"= 40.9%")),
               color = "black", size = 3, hjust = -0.05, vjust = -0.3) +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
      labs(title="T4 - Leaf")+
      labs(x="Axis.1 [41.8%]", y= "Axis.2 [13.4%]"),
    subset_samples(physeq_prok_ctr_ev,
                   Treatment1 %in% c("T4", "Control") & Slice %in% c("Nearsoil", "Control")) %>% 
      subset_samples(physeq = ., Treatment %in% c("T4")) %>% 
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
      PlotOrdinSplit(dataframe = .) +
      stat_ellipse(aes(x=Axis.1, y=Axis.2, group=Plant), inherit.aes = FALSE, type="norm", 
                   level=0.75, alpha=0.8, color = "black",linetype = 2, show.legend = FALSE) +
      annotate("text", -Inf, Inf, label = expression(paste(bold("Plant "), italic(R)^2,"= 12.2%")),
               color = "black",  size = 3, hjust = -0.1, vjust = 1.2) +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
      labs(title="T4 - Nearsoil")+
      labs(x="Axis.1 [12.6%]", y= "Axis.2 [8.1%]"),
    subset_samples(physeq_prok_ctr_ev,
                   Treatment1 %in% c("T4", "Control") & Slice %in% c("Farsoil", "Control")) %>% 
      subset_samples(physeq = ., Treatment %in% c("T4")) %>% 
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
      PlotOrdinSplit(dataframe = .) +
      stat_ellipse(aes(x=Axis.1, y=Axis.2, group=Pore), inherit.aes = FALSE, type="norm", 
                   level=0.75, alpha=0.8, color = "black",linetype = 2, show.legend = FALSE) +
      annotate("text", -Inf, Inf, label = expression(paste(bold("Pore "), italic(R)^2,"= 5.2%")),
               color = "black", size = 3, hjust = -0.1, vjust = 1.2) +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
      labs(title="T4 - Farsoil")+
      labs(x="Axis.1 [8.7%]", y= "Axis.2 [6.5%]"),
    control_T4,
    labels = c("A", "","","","B","","",""),
    ncol = 4,
    nrow = 2,
    legend = "bottom",
    common.legend=TRUE)

pcoa_split_prok

title2=text_grob("Bacterial Microbiome Composition",size=12, face=2)
grid.arrange(pcoa_split_prok, top = title2)


# ***********************************************---------------------------------
# ALPHA DIVERSITY ----------------------------------------------------------------

# extracting main data frames
ExtrDataFrame <- function(physeq){
  require(vegan)
  df_alpha <- as.data.frame(as.matrix(sample_data(physeq)))
  df_alpha$Observed <- specnumber(as.data.frame(otu_table(physeq)), MARGIN = 2)
  df_alpha$Shannon <- diversity(as.data.frame(otu_table(physeq)), index="shannon", MARGIN = 2)
  df_alpha$EH <- 1 - df_alpha$Shannon/log(df_alpha$Observed)
  return(df_alpha)
}

ExtrDataFrame(physeq_fungi_ctr_ev) %>% 
  subset(Treatment1 %in% c("T1", "T4", "Control"))-> df_alpha_fungi_ctr_all
df_alpha_fungi_ctr_all
dim(df_alpha_fungi_ctr_all)

ExtrDataFrame(physeq_prok_ctr_ev) %>% 
  subset(Treatment1 %in% c("T1", "T4", "Control"))-> df_alpha_prok_ctr_all
df_alpha_prok_ctr_all
dim(df_alpha_prok_ctr_all)

# This below works correctly
# New modified function to extract significance letters across groups
CompSampl <- function(df, formula){
  require(multcompView)
  require(tidyverse)
  df1 <-
    df %>%
    compare_means(data = ., formula, 
                  method = "wilcox.test", p.adjust.method = "bonferroni") %>% 
    select(group1, group2, p.adj)
  df2 <-
    df1 %>% 
    select(group2, group1, p.adj)
  colnames(df2) <- 
    c("group1", "group2", "p.adj") 
  rbind(df1, df2) %>% 
    as.data.frame() %>% 
    xtabs(p.adj ~ group2 + group1, data=.) %>% 
    as.dist(diag = TRUE) -> dist
  res <-
    as.data.frame(multcompLetters(x=dist, reversed = FALSE)['Letters'])
  return(res)
}


df_alpha_fungi_ctr_all %>%
  mutate(Slice = recode(Slice, Leaf = "A", Nearsoil = "B", Farsoil= "C", Control="D")) %>% 
  subset(x = ., Treatment1 %in% c("T1", "Control")) %>% 
  subset(x = ., Treatment %in% c("T1")) %>% 
  select(Observed, Slice, Pore) %>%
  unite("Variable", Slice, Pore, sep = "_", remove = FALSE) %>% 
  CompSampl(df = ., formula(Observed ~ Variable))

#df_alpha_fungi_ctr_all %>%
#  subset(x = ., Treatment1 %in% c("T1", "Control")) %>% 
#  subset(x = ., Treatment %in% c("T1")) %>% 
#  select(Observed, Slice, Pore) %>%
#  unite("Variable", Slice, Pore, sep = "_", remove = FALSE) %>% 
#  CompSampl(df = ., formula(Observed ~ Variable))

letters_fungi_T1_pore <-
df_alpha_fungi_ctr_all %>%
  mutate(Slice = recode(Slice, Leaf = "A", Nearsoil = "B", Farsoil= "C", Control="D")) %>% 
  subset(x = ., Treatment1 %in% c("T1", "Control")) %>% 
  subset(x = ., Treatment %in% c("T1")) %>% 
  select(Observed, Slice, Pore) %>%
  unite("Variable", Slice, Pore, sep = "_", remove = FALSE) %>% 
CompSampl(df = ., formula(Observed ~ Variable)) %>% 
  rownames_to_column("Variable") %>% 
  separate(Variable, c("Slice", "Pore"), sep = "_") %>% 
  mutate(Slice = recode(Slice, A = "Leaf", B = "Nearsoil", C = "Farsoil", D = "Control")) %>% 
mutate(Slice = as.factor(Slice),
       Slice = fct_relevel(Slice, "Leaf", "Nearsoil", "Farsoil", "Control"),
       Pore = as.factor(Pore),
       Pore = fct_relevel(Pore, "Small", "Large"), 
       Letters = as.factor(Letters)) %>% 
  arrange(factor(Slice, levels=c("Leaf", "Nearsoil", "Farsoil", "Control")))

letters_fungi_T1_pore


# Use enquo() and !! to pass function arguments to tidyverse
# functions. You can modify this and include and ordering argument
# and call is as, for example: Leaf = ordering[1] etc...
MakeLetters <- function(dataframe, formula, Treat, Div, Var, N){
  x <- enquo(Div)
  y <- enquo(Var)
  df_letters <-
    dataframe %>%
    mutate(Slice = recode(Slice, Leaf = "A", Nearsoil = "B", Farsoil= "C", Control="D")) %>% 
    subset(x = ., Treatment1 %in% c(Treat, "Control")) %>% 
    subset(x = ., Treatment %in% c(Treat)) %>% 
    select(Slice, !!x, !!y) %>% 
    unite("Variable", Slice, !!y, sep = "_", remove = FALSE) %>% 
    CompSampl(df = ., formula) %>% 
    rownames_to_column("Variable") %>% 
    separate(Variable, c("Slice", Var), sep = "_") %>% 
    mutate(Slice = recode(Slice, A = "Leaf", B = "Nearsoil", C = "Farsoil", D = "Control")) %>% 
    mutate(Group = rep(Treat, times=N)) %>% 
    mutate(Slice = as.factor(Slice),
           Slice = fct_relevel(Slice, "Leaf", "Nearsoil", "Farsoil", "Control"),
           Letters = as.factor(Letters)) %>% 
    arrange(factor(Slice, levels=c("Leaf", "Nearsoil", "Farsoil", "Control")))
  return(df_letters)
}


# Observed fungi
letters_fungi_T1_plant <-
  MakeLetters(df_alpha_fungi_ctr_all, formula(Observed ~ Variable), "T1", "Observed", "Plant", 7)
letters_fungi_T1_pore <-
  MakeLetters(df_alpha_fungi_ctr_all, formula(Observed ~ Variable),"T1", "Observed", "Pore", 8)
letters_fungi_T1_moist <-
  MakeLetters(df_alpha_fungi_ctr_all, formula(Observed ~ Variable),"T1", "Observed", "Moisture", 8)

letters_fungi_T4_plant <-
  MakeLetters(df_alpha_fungi_ctr_all,formula(Observed ~ Variable), "T4", "Observed", "Plant", 7)
letters_fungi_T4_pore <-
  MakeLetters(df_alpha_fungi_ctr_all,formula(Observed ~ Variable), "T4", "Observed", "Pore", 8)
letters_fungi_T4_moist <-
  MakeLetters(df_alpha_fungi_ctr_all, formula(Observed ~ Variable),"T4", "Observed", "Moisture", 8)

# observed bacteria
letters_prok_T1_plant <-
  MakeLetters(df_alpha_prok_ctr_all, formula(Observed ~ Variable),"T1", "Observed", "Plant", 7)
letters_prok_T1_pore <-
  MakeLetters(df_alpha_prok_ctr_all, formula(Observed ~ Variable),"T1", "Observed", "Pore", 8)
letters_prok_T1_moist <-
  MakeLetters(df_alpha_prok_ctr_all, formula(Observed ~ Variable),"T1", "Observed", "Moisture", 8)

letters_prok_T4_plant <-
  MakeLetters(df_alpha_prok_ctr_all,formula(Observed ~ Variable), "T4", "Observed", "Plant", 7)
letters_prok_T4_pore <-
  MakeLetters(df_alpha_prok_ctr_all,formula(Observed ~ Variable), "T4", "Observed", "Pore", 8)
letters_prok_T4_moist <-
  MakeLetters(df_alpha_prok_ctr_all,formula(Observed ~ Variable), "T4", "Observed", "Moisture", 8)

# shannon fungi
letters_fungi_T1_plant_shan <-
  MakeLetters(df_alpha_fungi_ctr_all,formula(EH ~ Variable), "T1", "EH", "Plant", 7)
letters_fungi_T1_pore_shan <-
  MakeLetters(df_alpha_fungi_ctr_all, formula(EH ~ Variable),"T1", "EH", "Pore", 8)
letters_fungi_T1_moist_shan <-
  MakeLetters(df_alpha_fungi_ctr_all, formula(EH ~ Variable),"T1", "EH", "Moisture", 8)

letters_fungi_T4_plant_shan <-
  MakeLetters(df_alpha_fungi_ctr_all, formula(EH ~ Variable),"T4", "EH", "Plant", 7)
letters_fungi_T4_pore_shan <-
  MakeLetters(df_alpha_fungi_ctr_all, formula(EH ~ Variable),"T4", "EH", "Pore", 8)
letters_fungi_T4_moist_shan <-
  MakeLetters(df_alpha_fungi_ctr_all, formula(EH ~ Variable),"T4", "EH", "Moisture", 8)

# shannon bacteria
letters_prok_T1_plant_shan <-
  MakeLetters(df_alpha_prok_ctr_all, formula(EH ~ Variable),"T1", "EH", "Plant", 7)
letters_prok_T1_pore_shan <-
  MakeLetters(df_alpha_prok_ctr_all, formula(EH ~ Variable),"T1", "EH", "Pore", 8)
letters_prok_T1_moist_shan <-
  MakeLetters(df_alpha_prok_ctr_all,formula(EH ~ Variable), "T1", "EH", "Moisture", 8)

letters_prok_T4_plant_shan <-
  MakeLetters(df_alpha_prok_ctr_all,formula(EH ~ Variable), "T4", "EH", "Plant", 7)
letters_prok_T4_pore_shan <-
  MakeLetters(df_alpha_prok_ctr_all,formula(EH ~ Variable), "T4", "EH", "Pore", 8)
letters_prok_T4_moist_shan <-
  MakeLetters(df_alpha_prok_ctr_all,formula(EH ~ Variable), "T4", "EH", "Moisture", 8)

# plotting function 
PlotRich <- function(dataframe, my_labels, Var, y_limit, y_labels){
  require(ggtext)
  #calculating where to put the letters
  #labels_y <-
  #  max(dataframe[,"mean"] + max(dataframe[,"sd"]) + 
  #        0.1 * max(dataframe[,"mean"])) 
  #y_limit <- labels_y * 1.15
  # plot
  rich_plot <-
    ggplot(dataframe, aes(x = get(Var), y = mean, 
               fill = Slice), 
               colour="black") +
    geom_bar(stat = "identity", position = position_dodge(), color="black") +
    geom_errorbar(aes(ymin = mean , ymax = mean + sd), color="black",
                  width = 0.2,
                  position = position_dodge(0.9)) +
    stat_summary(geom = "text", label = my_labels, fun= max, aes(y = y_labels), size=3, color="black") +
    facet_grid(~ Slice, scales = "free_x", space="free_x") +
    theme_bw() +
    expand_limits(y = c(0, y_limit)) +
    #scale_color_manual(values="black", 
    #                   labels = c("**LF** = Leaf", "**NS** = Nearsoil", "**FS** = Farsoil", "**C** = Control")) +
    scale_fill_manual(values=paletteCB5_tris, 
                      labels = c("**LF** = Leaf", "**NS** = Nearsoil", "**FS** = Farsoil", "**C** = Control")) +
    theme(strip.text.x = element_text(size = 9, angle = 0, hjust = 0.5, vjust = 0.5),
          strip.text = element_text(size=8, face = "bold"),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          legend.key.height = unit(0.2, "cm"), 
          legend.key.width = unit(0.3, "cm"),
          legend.title = element_blank(), 
          legend.text = element_markdown(size = 9),
          axis.text.x = element_text(angle = 45, size = 8, hjust = 1, vjust = 1.05),
          axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5),
          axis.title = element_text(angle = 0, size = 8, face = "bold"),
          strip.background = element_blank(), 
          legend.position = "none") +
    grids(linetype = "dashed") 
    #guides(fill = guide_legend(order = 1, nrow = 1, title = "Niche",
    #                          override.aes = list(shape = 22, size=3)))
  return(rich_plot)
}



test_alpha <-
  df_alpha_fungi_ctr_all %>%
  subset(Treatment1 %in% c("T1", "Control")) %>% 
  subset(Treatment %in% c("T1")) %>% 
  mutate(Slice = as.factor(Slice),
         Slice = recode(Slice, Leaf = "T1\nLeaf",  Nearsoil = "T1\nNearsoil", Farsoil= "T1\nfarsoil", Control="T1\nControl"),
         Slice = fct_relevel(Slice, "T1\nLeaf", "T1\nNearsoil", "T1\nfarsoil", "T1\nControl"),
         Moisture = as.factor(Moisture), 
         Moisture = recode(Moisture, CBP = "Low", FC = "High"),
         Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
  select(Slice, Observed, Pore) %>%
  group_by(Slice, Pore) %>%
  summarise_each(funs(mean, sd)) %>%
  PlotRich(Var = "Pore", letters_fungi_T1_pore$Letters, y_limit = 150, y_labels=150) 

test_alpha

# Fungi Observed richness ------------------------------------------------------
observed_fungi <-
ggarrange(
df_alpha_fungi_ctr_all %>%
  subset(Treatment1 %in% c("T1", "Control")) %>% 
  subset(Treatment %in% c("T1")) %>% 
  mutate(Slice = as.factor(Slice),
         Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
         Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
         Moisture = as.factor(Moisture), 
         Moisture = recode(Moisture, CBP = "Low", FC = "High"),
         Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
  select(Slice, Observed, Plant) %>%
  group_by(Slice, Plant) %>%
  summarise_each(funs(mean, sd)) %>%
  PlotRich(Var = "Plant", letters_fungi_T1_plant$Letters, y_limit = 150, y_labels=150) +
  labs(title = NULL, x="Plant", y="Mean richness"), 
df_alpha_fungi_ctr_all %>%
  subset(Treatment1 %in% c("T1", "Control")) %>% 
  subset(Treatment %in% c("T1")) %>%
  mutate(Slice = as.factor(Slice),
         Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
         Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
         Moisture = as.factor(Moisture), 
         Moisture = recode(Moisture, CBP = "Low", FC = "High"),
         Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
  select(Slice, Observed, Pore) %>%
  group_by(Slice, Pore) %>%
  summarise_each(funs(mean, sd)) %>%
  PlotRich(Var = "Pore", letters_fungi_T1_pore$Letters, y_limit = 150, y_labels=150) +
  labs(title = "T1", x="Pore", y="Mean richness"),
df_alpha_fungi_ctr_all %>%
  subset(Treatment1 %in% c("T1", "Control")) %>% 
  subset(Treatment %in% c("T1")) %>%
  mutate(Slice = as.factor(Slice),
         Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
         Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
         Moisture = as.factor(Moisture), 
         Moisture = recode(Moisture, CBP = "Low", FC = "High"),
         Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
  select(Slice, Observed, Moisture) %>%
  group_by(Slice, Moisture) %>%
  summarise_each(funs(mean, sd)) %>%
  PlotRich(Var = "Moisture", letters_fungi_T1_moist$Letters, y_limit = 150, y_labels=150) +
  labs(title = NULL, x="Moisture", y="Mean richness"),
df_alpha_fungi_ctr_all %>%
  subset(Treatment1 %in% c("T4", "Control")) %>% 
  subset(Treatment %in% c("T4")) %>% 
  mutate(Slice = as.factor(Slice),
         Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
         Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
         Moisture = as.factor(Moisture), 
         Moisture = recode(Moisture, CBP = "Low", FC = "High"),
         Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
  select(Slice, Observed, Plant) %>%
  group_by(Slice, Plant) %>%
  summarise_each(funs(mean, sd)) %>%
  PlotRich(Var = "Plant", letters_fungi_T4_plant$Letters, y_limit = 150, y_labels=150) +
  labs(title = NULL, x="Plant", y="Mean richness"),
df_alpha_fungi_ctr_all %>%
  subset(Treatment1 %in% c("T4", "Control")) %>% 
  subset(Treatment %in% c("T4")) %>% 
  mutate(Slice = as.factor(Slice),
         Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
         Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
         Moisture = as.factor(Moisture), 
         Moisture = recode(Moisture, CBP = "Low", FC = "High"),
         Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
  select(Slice, Observed, Pore) %>%
  group_by(Slice, Pore) %>%
  summarise_each(funs(mean, sd)) %>%
  PlotRich(Var = "Pore", letters_fungi_T4_pore$Letters, y_limit = 150, y_labels=150) +
  labs(title = "T4", x="Pore", y="Mean richness"),
df_alpha_fungi_ctr_all %>%
  subset(Treatment1 %in% c("T4", "Control")) %>% 
  subset(Treatment %in% c("T4")) %>% 
  mutate(Slice = as.factor(Slice),
         Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
         Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
         Moisture = as.factor(Moisture), 
         Moisture = recode(Moisture, CBP = "Low", FC = "High"),
         Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
  select(Slice, Observed, Moisture) %>%
  group_by(Slice, Moisture) %>%
  summarise_each(funs(mean, sd)) %>%
  PlotRich(Var = "Moisture", letters_fungi_T4_moist$Letters, y_limit = 150, y_labels=150) +
  labs(title = NULL, x="Moisture", y="Mean richness"),
labels = c("A","","","B","",""),
align = "h",
widths = c(0.31,0.34,0.34),
ncol = 3,
nrow = 2,
legend = "bottom",
common.legend=TRUE)
 
observed_fungi

title1=text_grob("Fungal Observed Richness",size=12, face=2)
grid.arrange(observed_fungi, top = title1)

observed_fungi + 
  plot_annotation(
    title = "Fungal Observed Richness", 
    theme = theme(plot.title = element_text(
      size = 12, face="bold", hjust = 0.5))) +
  theme(legend.position = "none")


# Bacteria Observed richness ------------------------------------------------------
observed_bact <-
  ggarrange(
    df_alpha_prok_ctr_all %>%
      subset(Treatment1 %in% c("T1", "Control")) %>% 
      subset(Treatment %in% c("T1")) %>% 
      mutate(Slice = as.factor(Slice),
             Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
             Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
             Moisture = as.factor(Moisture), 
             Moisture = recode(Moisture, CBP = "Low", FC = "High"),
             Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
      select(Slice, Observed, Plant) %>%
      group_by(Slice, Plant) %>%
      summarise_each(funs(mean, sd)) %>%
      PlotRich(Var = "Plant", letters_prok_T1_plant$Letters, y_limit = 2500, y_labels=2500) +
      labs(title = "T1 - Plant", x=NULL, y="Mean richness"), 
    df_alpha_prok_ctr_all %>%
      subset(Treatment1 %in% c("T1", "Control")) %>% 
      subset(Treatment %in% c("T1")) %>% 
      mutate(Slice = as.factor(Slice),
             Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
             Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
             Moisture = as.factor(Moisture), 
             Moisture = recode(Moisture, CBP = "Low", FC = "High"),
             Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
      select(Slice, Observed, Pore) %>%
      group_by(Slice, Pore) %>%
      summarise_each(funs(mean, sd)) %>%
      PlotRich(Var = "Pore", letters_prok_T1_pore$Letters, y_limit = 2500, y_labels=2500) +
      labs(title = "T1 - Pore", x=NULL, y="Mean richness"),
    df_alpha_prok_ctr_all %>%
      subset(Treatment1 %in% c("T1", "Control")) %>% 
      subset(Treatment %in% c("T1")) %>% 
      mutate(Slice = as.factor(Slice),
             Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
             Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
             Moisture = as.factor(Moisture), 
             Moisture = recode(Moisture, CBP = "Low", FC = "High"),
             Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
      select(Slice, Observed, Moisture) %>%
      group_by(Slice, Moisture) %>%
      summarise_each(funs(mean, sd)) %>%
      PlotRich(Var = "Moisture", letters_prok_T1_moist$Letters, y_limit = 2500, y_labels=2500) +
      labs(title = "T1 - Moisture", x=NULL, y="Mean richness"),
    df_alpha_prok_ctr_all %>%
      subset(Treatment1 %in% c("T4", "Control")) %>% 
      subset(Treatment %in% c("T4")) %>% 
      mutate(Slice = as.factor(Slice),
             Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
             Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
             Moisture = as.factor(Moisture), 
             Moisture = recode(Moisture, CBP = "Low", FC = "High"),
             Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
      select(Slice, Observed, Plant) %>%
      group_by(Slice, Plant) %>%
      summarise_each(funs(mean, sd)) %>%
      PlotRich(Var = "Plant", letters_prok_T4_plant$Letters, y_limit = 2500, y_labels=2500) +
      labs(title = "T4 - Plant", x=NULL, y="Mean richness"),
    df_alpha_prok_ctr_all %>%
      subset(Treatment1 %in% c("T4", "Control")) %>% 
      subset(Treatment %in% c("T4")) %>% 
      mutate(Slice = as.factor(Slice),
             Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
             Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
             Moisture = as.factor(Moisture), 
             Moisture = recode(Moisture, CBP = "Low", FC = "High"),
             Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
      select(Slice, Observed, Pore) %>%
      group_by(Slice, Pore) %>%
      summarise_each(funs(mean, sd)) %>%
      PlotRich(Var = "Pore", letters_prok_T4_pore$Letters, y_limit = 2500, y_labels=2500) +
      labs(title = "T4 - Pore", x=NULL, y="Mean richness"),
    df_alpha_prok_ctr_all %>%
      subset(Treatment1 %in% c("T4", "Control")) %>% 
      subset(Treatment %in% c("T4")) %>% 
      mutate(Slice = as.factor(Slice),
             Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
             Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
             Moisture = as.factor(Moisture), 
             Moisture = recode(Moisture, CBP = "Low", FC = "High"),
             Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
      select(Slice, Observed, Moisture) %>%
      group_by(Slice, Moisture) %>%
      summarise_each(funs(mean, sd)) %>%
      PlotRich(Var = "Moisture", letters_prok_T4_moist$Letters, y_limit = 2500, y_labels=2500) +
      labs(title = "T4 - Moisture", x=NULL, y="Mean richness"),
    labels = c("C","","","D","",""),
    align = "h",
    widths = c(0.31,0.34,0.34),
    ncol = 3,
    nrow = 2,
    legend = "bottom",
    common.legend=TRUE)

observed_bact

title2=text_grob("Bacterial Observed Richness",size=12, face=2)
grid.arrange(observed_bact, top = title2)

observed_bact + 
  plot_annotation(
    title = "Bacterial Observed Richness", 
    theme = theme(plot.title = element_text(
      size = 12, face="bold", hjust = 0.5)))


# Fungi EH Shannon index ------------------------------------------------------
EH_fungi <-
  ggarrange(
    df_alpha_fungi_ctr_all %>%
      subset(Treatment1 %in% c("T1", "Control")) %>% 
      subset(Treatment %in% c("T1")) %>% 
      mutate(Slice = as.factor(Slice),
             Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
             Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
             Moisture = as.factor(Moisture), 
             Moisture = recode(Moisture, CBP = "Low", FC = "High"),
             Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
      select(Slice, EH, Plant) %>%
      group_by(Slice, Plant) %>%
      summarise_each(funs(mean, sd)) %>%
      PlotRich(Var = "Plant", letters_fungi_T1_plant_shan$Letters, y_limit = 0.6, y_labels=0.6) +
      labs(title = "T1 - Plant", x=NULL, y="Mean Shannon"), 
    df_alpha_fungi_ctr_all %>%
      subset(Treatment1 %in% c("T1", "Control")) %>% 
      subset(Treatment %in% c("T1")) %>% 
      mutate(Slice = as.factor(Slice),
             Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
             Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
             Moisture = as.factor(Moisture), 
             Moisture = recode(Moisture, CBP = "Low", FC = "High"),
             Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
      select(Slice, EH, Pore) %>%
      group_by(Slice, Pore) %>%
      summarise_each(funs(mean, sd)) %>%
      PlotRich(Var = "Pore", letters_fungi_T1_pore_shan$Letters, y_limit = 0.6, y_labels=0.6) +
      labs(title = "T1 - Pore", x=NULL, y="Mean Shannon"),
    df_alpha_fungi_ctr_all %>%
      subset(Treatment1 %in% c("T1", "Control")) %>% 
      subset(Treatment %in% c("T1")) %>% 
      mutate(Slice = as.factor(Slice),
             Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
             Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
             Moisture = as.factor(Moisture), 
             Moisture = recode(Moisture, CBP = "Low", FC = "High"),
             Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
      select(Slice, EH, Moisture) %>%
      group_by(Slice, Moisture) %>%
      summarise_each(funs(mean, sd)) %>%
      PlotRich(Var = "Moisture", letters_fungi_T1_moist_shan$Letters, y_limit = 0.6, y_labels=0.6) +
      labs(title = "T1 - Moisture", x=NULL, y="Mean Shannon"),
    df_alpha_fungi_ctr_all %>%
      subset(Treatment1 %in% c("T4", "Control")) %>% 
      subset(Treatment %in% c("T4")) %>% 
      mutate(Slice = as.factor(Slice),
             Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
             Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
             Moisture = as.factor(Moisture), 
             Moisture = recode(Moisture, CBP = "Low", FC = "High"),
             Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
      select(Slice, EH, Plant) %>%
      group_by(Slice, Plant) %>%
      summarise_each(funs(mean, sd)) %>%
      PlotRich(Var = "Plant", letters_fungi_T4_plant_shan$Letters, y_limit = 0.6, y_labels=0.6) +
      labs(title = "T4 - Plant", x=NULL, y="Mean Shannon"),
    df_alpha_fungi_ctr_all %>%
      subset(Treatment1 %in% c("T4", "Control")) %>% 
      subset(Treatment %in% c("T4")) %>% 
      mutate(Slice = as.factor(Slice),
             Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
             Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
             Moisture = as.factor(Moisture), 
             Moisture = recode(Moisture, CBP = "Low", FC = "High"),
             Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
      select(Slice, EH, Pore) %>%
      group_by(Slice, Pore) %>%
      summarise_each(funs(mean, sd)) %>%
      PlotRich(Var = "Pore", letters_fungi_T4_pore_shan$Letters, y_limit = 0.6, y_labels=0.6) +
      labs(title = "T4 - Pore", x=NULL, y="Mean Shannon"),
    df_alpha_fungi_ctr_all %>%
      subset(Treatment1 %in% c("T4", "Control")) %>% 
      subset(Treatment %in% c("T4")) %>% 
      mutate(Slice = as.factor(Slice),
             Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
             Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
             Moisture = as.factor(Moisture), 
             Moisture = recode(Moisture, CBP = "Low", FC = "High"),
             Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
      select(Slice, EH, Moisture) %>%
      group_by(Slice, Moisture) %>%
      summarise_each(funs(mean, sd)) %>%
      PlotRich(Var = "Moisture", letters_fungi_T4_moist_shan$Letters, y_limit = 0.6, y_labels=0.6) +
      labs(title = "T4 - Moisture", x=NULL, y="Mean Shannon"),
    labels = c("A", "B","C","D","E","F"),
    align = "h",
    widths = c(0.31,0.34,0.34),
    ncol = 3,
    nrow = 2,
    legend = "bottom",
    common.legend=TRUE)

EH_fungi

title4=text_grob("Fungal Shannon Index",size=12, face=2)
grid.arrange(EH_fungi, top = title4)


# bacteria EH Shannon index ------------------------------------------------------
EH_bacteria <-
  ggarrange(
    df_alpha_prok_ctr_all %>%
      subset(Treatment1 %in% c("T1", "Control")) %>% 
      subset(Treatment %in% c("T1")) %>% 
      mutate(Slice = as.factor(Slice),
             Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
             Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
             Moisture = as.factor(Moisture), 
             Moisture = recode(Moisture, CBP = "Low", FC = "High"),
             Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
      select(Slice, EH, Plant) %>%
      group_by(Slice, Plant) %>%
      summarise_each(funs(mean, sd)) %>%
      PlotRich(Var = "Plant", letters_prok_T1_plant_shan$Letters, y_limit = 0.6, y_labels=0.6) +
      labs(title = "T1 - Plant", x=NULL, y="Mean Shannon"), 
    df_alpha_prok_ctr_all %>%
      subset(Treatment1 %in% c("T1", "Control")) %>% 
      subset(Treatment %in% c("T1")) %>% 
      mutate(Slice = as.factor(Slice),
             Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
             Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
             Moisture = as.factor(Moisture), 
             Moisture = recode(Moisture, CBP = "Low", FC = "High"),
             Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
      select(Slice, EH, Pore) %>%
      group_by(Slice, Pore) %>%
      summarise_each(funs(mean, sd)) %>%
      PlotRich(Var = "Pore", letters_prok_T1_pore_shan$Letters, y_limit = 0.6, y_labels=0.6) +
      labs(title = "T1 - Pore", x=NULL, y="Mean Shannon"),
    df_alpha_prok_ctr_all %>%
      subset(Treatment1 %in% c("T1", "Control")) %>% 
      subset(Treatment %in% c("T1")) %>% 
      mutate(Slice = as.factor(Slice),
             Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
             Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
             Moisture = as.factor(Moisture), 
             Moisture = recode(Moisture, CBP = "Low", FC = "High"),
             Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
      select(Slice, EH, Moisture) %>%
      group_by(Slice, Moisture) %>%
      summarise_each(funs(mean, sd)) %>%
      PlotRich(Var = "Moisture", letters_prok_T1_moist_shan$Letters, y_limit = 0.6, y_labels=0.6) +
      labs(title = "T1 - Moisture", x=NULL, y="Mean Shannon"),
    df_alpha_prok_ctr_all %>%
      subset(Treatment1 %in% c("T4", "Control")) %>% 
      subset(Treatment %in% c("T4")) %>% 
      mutate(Slice = as.factor(Slice),
             Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
             Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
             Moisture = as.factor(Moisture), 
             Moisture = recode(Moisture, CBP = "Low", FC = "High"),
             Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
      select(Slice, EH, Plant) %>%
      group_by(Slice, Plant) %>%
      summarise_each(funs(mean, sd)) %>%
      PlotRich(Var = "Plant", letters_prok_T4_plant_shan$Letters, y_limit = 0.6, y_labels=0.6) +
      labs(title = "T4 - Plant", x=NULL, y="Mean Shannon"),
    df_alpha_prok_ctr_all %>%
      subset(Treatment1 %in% c("T4", "Control")) %>% 
      subset(Treatment %in% c("T4")) %>% 
      mutate(Slice = as.factor(Slice),
             Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
             Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
             Moisture = as.factor(Moisture), 
             Moisture = recode(Moisture, CBP = "Low", FC = "High"),
             Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
      select(Slice, EH, Pore) %>%
      group_by(Slice, Pore) %>%
      summarise_each(funs(mean, sd)) %>%
      PlotRich(Var = "Pore", letters_prok_T4_pore_shan$Letters, y_limit = 0.6, y_labels=0.6) +
      labs(title = "T4 - Pore", x=NULL, y="Mean Shannon"),
    df_alpha_prok_ctr_all %>%
      subset(Treatment1 %in% c("T4", "Control")) %>% 
      subset(Treatment %in% c("T4")) %>% 
      mutate(Slice = as.factor(Slice),
             Slice = recode(Slice, Leaf = "LF",  Nearsoil = "NS", Farsoil= "FS", Control="C"),
             Slice = fct_relevel(Slice, "LF", "NS", "FS", "C"),
             Moisture = as.factor(Moisture), 
             Moisture = recode(Moisture, CBP = "Low", FC = "High"),
             Moisture = fct_relevel(Moisture, "Low", "High")) %>% 
      select(Slice, EH, Moisture) %>%
      group_by(Slice, Moisture) %>%
      summarise_each(funs(mean, sd)) %>%
      PlotRich(Var = "Moisture", letters_prok_T4_moist_shan$Letters, y_limit = 0.6, y_labels=0.6) +
      labs(title = "T4 - Moisture", x=NULL, y="Mean Shannon"),
    labels = c("A", "B","C","D","E","F"),
    align = "h",
    widths = c(0.31,0.34,0.34),
    ncol = 3,
    nrow = 2,
    legend = "bottom",
    common.legend=TRUE)

EH_bacteria

title4=text_grob("Bacterial Shannon Index",size=12, face=2)
grid.arrange(EH_bacteria, top = title4)


# ***********************************************---------------------------------
# Reformat Taxonomy --------------------------------------------------------------
source("../R_functions/ReformatTaxonomy.R")

#readRDS("../project_Sasha&Terry_microcosms/Random_Forest/fungi_physeq.rds") -> physeq_fungi_ev
#readRDS("../project_Sasha&Terry_microcosms/Random_Forest/bacteria_physeq.rds") -> physeq_prok_ev

constax_fungi <-as.matrix(physeq_fungi_ev@tax_table)
constax_fungi <- as.data.frame(constax_fungi)
constax_fungi <- 
  constax_fungi %>% rownames_to_column(var = "OTU")
head(constax_fungi)

constax_prok <-as.matrix(physeq_prok_ev@tax_table)
constax_prok <- as.data.frame(constax_prok)


FinalizeTaxonomy <- function(constax){
  constax$Species <- 
    gsub(" sp ", "", constax$Species)
  constax[] = lapply(constax, blank2na, na.strings=c('','NA','na','N/A','n/a','NaN','nan'))
  lastValue <- function(x) tail(x[!is.na(x)], 1)
  constax$Genus <- as.character(constax$Genus)
  #constax[which(is.na(constax$Genus) == FALSE),]$Genus <-
  #  paste(constax$Genus[is.na(constax$Genus) == FALSE], "sp.", sep = " ")
  last_taxons<- apply(constax[,c(1:7)], 1, lastValue)
  constax$BestMatch <- last_taxons
  constax$BestMatch <-
    gsub("_", " ", constax$BestMatch)
  constax$Taxonomy <-
    paste(constax[,1], constax$BestMatch, sep = "-")
  constax$Taxonomy <- 
    gsub("_", "", constax$Taxonomy)
  constax$BestMatch <- 
    gsub(" sp.", "", constax$BestMatch)
  constax$Genus <- 
    gsub(" sp.", "", constax$Genus)
  return(constax)
}

# Fungi ------------------------------------------------------------------------

# Fixing a problem at genus level that has species - Class missing!
filter(constax_fungi, grepl(" ",Genus))

constax_fungi_new <-
  constax_fungi %>% 
  mutate(Genus = ifelse(Genus %in% "Coprinellus verrucispermus", paste("Coprinellus"), paste(Genus))) %>%
  mutate(Genus = ifelse(Genus %in% "Rhizophlyctis rosea", paste("Rhizophlyctis"), paste(Genus))) %>%
  mutate(Genus = ifelse(Genus %in% "Sphaerostilbella sp", paste("Sphaerostilbella"), paste(Genus))) %>%
  mutate(Genus = ifelse(Genus %in% "Coprinellus micaceus", paste("Coprinellus"), paste(Genus)))

constax_fungi_new$Genus

taxonomy_fungi <- 
  FinalizeTaxonomy(constax = constax_fungi_new)
taxonomy_fungi[1:50, ]

# Bacteria ---------------------------------------------------------------------
constax_prok_filt <-
  constax_prok %>%
  select(-repseq,-OTUS) %>%
  rownames_to_column(var = "OTU") %>%
  select(OTU = OTU,
         Kingdom = domain, 
         Phylum = phylum, 
         Class = class, 
         Order= order, 
         Family = family,
         Genus= genus,
         Species = genus)

  
unique(constax_prok_filt$Class)
table(constax_prok_filt$Class)

taxonomy_prok <- 
  FinalizeTaxonomy(constax = constax_prok_filt) %>% 
  mutate(Taxonomy = gsub("unclassified_Root", "Bacteria", Taxonomy),
         Taxonomy = gsub(" genera incertae sedis", "", Taxonomy),
         Taxonomy = gsub("unclassified_", "", Taxonomy),
         Taxonomy = gsub("candidate division", "cand. div.", Taxonomy),
         Taxonomy = gsub("Candidatus", "Cand.", Taxonomy),
         Class = gsub("unclassified_Root", "Bacteria", Class),
         Class = gsub("unclassified_", "", Class),
         Class = gsub("candidate division", "cand. div.", Class),
         Class = gsub("Candidatus", "Cand.", Class))

head(taxonomy_prok)
unique(taxonomy_prok$Class)
table(taxonomy_prok$Class)
dim(taxonomy_prok)

taxonomy_prok %>% 
subset(Class %in% c("Root")) %>% 
  as.data.frame()



# ***********************************************---------------------------------
# DIFFERENTIALLY ABUNDANT OTUS ----------------------------------------------------
# Identifying populations associated with a particular variable

# >>>> FUNGI <<<<<<---------------------------------------------------------------------
# Extracting otu_table -----------------------------------------------------------
fungi_otu_ctr_rf <-
  physeq_fungi_ctr_ev@otu_table %>%
  t(x = .) %>%
  as.data.frame(x = .) %>%
  mutate(Slice = as.factor(sample_data(physeq_fungi_ctr_ev)$Slice),
         Treatment = as.factor(sample_data(physeq_fungi_ctr_ev)$Treatment),
         Treatment1 = as.factor(sample_data(physeq_fungi_ctr_ev)$Treatment1),
         Pore = as.factor(sample_data(physeq_fungi_ctr_ev)$Pore),
         Moisture = as.factor(sample_data(physeq_fungi_ctr_ev)$Moisture),
         Plant = as.factor(sample_data(physeq_fungi_ctr_ev)$Plant),
         Time = as.factor(sample_data(physeq_fungi_ctr_ev)$Time)) %>% 
  subset(x = ., Treatment1 %in% c("T1", "T4", "Control")) %>% 
  droplevels(x = .)

colnames(fungi_otu_ctr_rf)
class(fungi_otu_ctr_rf)
tail(fungi_otu_ctr_rf)
dim(fungi_otu_ctr_rf)
head(fungi_otu_ctr_rf)

table(fungi_otu_ctr_rf$Treatment1)
table(droplevels(fungi_otu_ctr_rf)$Treatment1)

fungi_otu_ctr_rf$Slice

# Adding taxonomy ----------------------------------------------------------------
df_abund_fungi <-
  fungi_otu_ctr_rf %>%
  select(.data = ., starts_with("FOTU")) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Sum = rowSums(x = .)) %>%
  full_join(x = rownames_to_column(., var = "OTU"),
            y = taxonomy_fungi,
            by = "OTU") %>%
  arrange(., desc(Sum))

head(df_abund_fungi)
dim(df_abund_fungi)

# transforming to long format --------------------------------------------------
# deprecated: summarise_each(funs(sum))
df_abund_fungi_otu <-
  left_join(
    df_abund_fungi %>%
    select(.data = ., starts_with("R") | starts_with("Control"), Taxonomy) %>%
    na.omit() %>% 
    as_tibble() %>%
    group_by(Taxonomy) %>% 
    summarize(across(everything(), sum)) %>%
    ungroup(),
    df_abund_fungi %>% 
      select(Genus, BestMatch, Class, Taxonomy),
    by="Taxonomy") %>%
    pivot_longer(cols=c(-Taxonomy, -Class, -Genus, -BestMatch),
                 names_to="SampleID",
                 values_to = "Count") %>% 
  full_join(
    x = .,
    y = rownames_to_column(
      select(fungi_otu_ctr_rf,
             Slice, Treatment, Treatment1, 
             Pore, Plant, Moisture, Time), var = "SampleID"),
    by = "SampleID") %>% 
  group_by(Taxonomy)

df_abund_fungi_otu


table(
  df_abund_fungi_otu$Slice,
  df_abund_fungi_otu$Treatment1,
  df_abund_fungi_otu$Pore
)


table(
  df_abund_fungi_otu$Slice,
  df_abund_fungi_otu$Treatment1
)

# sanity check
df_abund_fungi_otu %>% 
  subset(df_abund_fungi_otu$Taxonomy=="FOTU279-Pyronemataceae") %>% 
  as.data.frame() 

df_abund_fungi_otu %>% 
  subset(df_abund_fungi_otu$BestMatch=="Pleosporales") %>% 
  as.data.frame()

df_abund_fungi_otu %>% 
  subset(df_abund_fungi_otu$Taxonomy=="FOTU1-Pleosporales") %>% 
  as.data.frame()

df_abund_fungi_otu

# Reintroducing the control samples --------------------------------------------
# filter(Class!="incertae_sedis") # Remove incertae_sedis at class level 

df_abund_fungi_otu_add <-
  rbind(
    df_abund_fungi_otu %>% 
      filter(.data = ., Slice %in% "Control") %>% 
      unite("Pore", Slice, Pore, sep = "_", remove = FALSE) %>% 
      unite("Moisture", Slice, Moisture, sep = "_", remove = FALSE),
    df_abund_fungi_otu %>% 
      filter(.data = ., !Slice %in% "Control")
  ) %>% 
  mutate(Pore = as.factor(Pore),
         Moisture = as.factor(Moisture)) 
  
df_abund_fungi_otu_add

table(
  df_abund_fungi_otu_add$Pore,
  df_abund_fungi_otu_add$Treatment1
  )

table(
  df_abund_fungi_otu_add$Treatment1,
  df_abund_fungi_otu_add$Slice
)


# The nth percentile of a data set is the value that cuts off the first n percent 
# of the data values when all of the values are sorted from least to greatest.
# We used the 50th percentile for the CV, to reduce the number of taxa.

# CV = coefficient of variation
cv <- function(x){
  y <- sd(x) / mean(x) * 100
  y[is.nan(y)] <- NA
  return(y)
  }

#otu_fungi_highcv <-
#df_abund_fungi_otu_add %>% 
#  select(Taxonomy, Count) %>% 
#  group_by(Taxonomy) %>% 
#  summarise_each(funs(mean, sd, cv)) %>% 
#  ungroup() %>% 
#  select(Taxonomy, cv) %>% 
#  na.omit() %>% 
#  filter(cv > quantile(cv, probs = c(0.50)))
#otu_fungi_highcv

# METHOD USED IN THE ANALYSIS. Different approach here. I filtered to the 50th 
# percentile for the mean and of the CV, to reduce the number of taxa. Then I 
# sorted to the most abundant taxa within each Class and pick the first 5.

otu_fungi_highcv <-
  df_abund_fungi_otu_add %>% 
  select(Taxonomy, Class, Count) %>% 
  group_by(Taxonomy, Class) %>% 
  summarise_each(funs(mean, sd, cv)) %>% 
  ungroup() %>% 
  na.omit() %>% 
  filter(cv > quantile(cv, probs = c(0.50))) %>% 
  arrange(Class, desc(mean)) %>% 
  group_by(Class) %>% 
  top_n(n = 5, wt = mean) %>% 
  ungroup()

otu_fungi_highcv$Class  


otu_fungi_highcv_top <-
  df_abund_fungi_otu_add %>% 
  select(Taxonomy, Class, Count) %>% 
  group_by(Taxonomy, Class) %>% 
  summarise_each(funs(mean, sd, cv)) %>% 
  ungroup() %>% 
  na.omit() %>% 
  filter(cv > quantile(cv, probs = c(0.50))) %>% 
  arrange(desc(mean)) %>% 
  top_n(n = 25, wt = mean)


otu_fungi_highcv_top$Class  

#otu_fungi_highcv <-
#df_abund_fungi_otu_add %>% 
#  select(Taxonomy, SampleID, Count) %>% 
#  pivot_wider(Taxonomy, names_from = SampleID, values_from = Count) %>% 
#  group_by(Taxonomy) %>% 
#  mutate(Sum = rowSums(across(where(is.numeric)))) %>% 
#  ungroup() %>% 
#  filter(Sum > quantile(Sum, probs = 0.50)) 


# >>>> BACTERIA <<<<<<---------------------------------------------------------------------
# Extracting otu_table -----------------------------------------------------------
physeq_prok_ctr_ev@sam_data

prok_otu_ctr_rf <-
  physeq_prok_ctr_ev@otu_table %>%
  t(x = .) %>%
  as.data.frame(x = .) %>%
  mutate(Slice = as.factor(sample_data(physeq_prok_ctr_ev)$Slice),
         Treatment = as.factor(sample_data(physeq_prok_ctr_ev)$Treatment),
         Treatment1 = as.factor(sample_data(physeq_prok_ctr_ev)$Treatment1),
         Pore = as.factor(sample_data(physeq_prok_ctr_ev)$Pore),
         Moisture = as.factor(sample_data(physeq_prok_ctr_ev)$Moisture),
         Plant = as.factor(sample_data(physeq_prok_ctr_ev)$Plant),
         Time = as.factor(sample_data(physeq_prok_ctr_ev)$Time), 
         SampleID = as.factor(sample_data(physeq_prok_ctr_ev)$Description)) %>% 
  subset(x = ., Treatment1 %in% c("T1", "T4", "Control")) %>% 
  rownames_to_column("OriginalID") %>% 
  column_to_rownames("SampleID") %>% 
  droplevels(x = .)

colnames(prok_otu_ctr_rf)
class(prok_otu_ctr_rf)
tail(prok_otu_ctr_rf)
dim(prok_otu_ctr_rf)
head(prok_otu_ctr_rf)
table(prok_otu_ctr_rf$Treatment1)
table(droplevels(prok_otu_ctr_rf)$Treatment1)
prok_otu_ctr_rf$Slice
prok_otu_ctr_rf[1:6, 0:20]
prok_otu_ctr_rf[1:6, 30000:length(prok_otu_ctr_rf)]

# Adding taxonomy ----------------------------------------------------------------
# Needed to fix the sample names before workin on this.
head(taxonomy_prok)

df_abund_prok <-
  prok_otu_ctr_rf %>%
  select(.data = ., starts_with("POTU")) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Sum = rowSums(x = .)) %>%
  full_join(x = rownames_to_column(., var = "OTU"),
            y = taxonomy_prok,
            by = "OTU") %>%
  arrange(., desc(Sum))

head(df_abund_prok)
dim(df_abund_prok)

df_abund_prok %>% 
  subset(df_abund_prok$Class %in% c("Root")) %>% 
  as.data.frame() 


# transforming to long format --------------------------------------------------
# deprecated: summarise_each(funs(sum))
df_abund_prok_otu <-
  left_join(
    df_abund_prok %>%
      select(.data = ., starts_with("R") | starts_with("Control"), Taxonomy) %>% 
      na.omit() %>% 
      as_tibble() %>%
      group_by(Taxonomy) %>% 
      summarize(across(everything(), sum)) %>%
      ungroup(),
    df_abund_prok %>% 
      select(Genus, Class, BestMatch, Taxonomy),
    by="Taxonomy") %>%
  pivot_longer(cols=c(-Taxonomy, -Class, -Genus, -BestMatch),
               names_to="SampleID",
               values_to = "Count") %>% 
  full_join(
    x = .,
    y = rownames_to_column(
      select(prok_otu_ctr_rf,
             Slice, Treatment, Treatment1, 
             Pore, Plant, Moisture, Time), var = "SampleID"),
    by = "SampleID") %>% 
  group_by(Taxonomy)

df_abund_prok_otu

table(
  df_abund_prok_otu$Slice,
  df_abund_prok_otu$Treatment1,
  df_abund_prok_otu$Pore
)


table(
  df_abund_prok_otu$Slice,
  df_abund_prok_otu$Treatment1
)

# Keep the controls in each variable 
df_abund_prok_otu_add <-
  rbind(
    df_abund_prok_otu %>% 
      filter(.data = ., Slice %in% "Control") %>% 
      unite("Pore", Slice, Pore, sep = "_", remove = FALSE) %>% 
      unite("Moisture", Slice, Moisture, sep = "_", remove = FALSE),
    df_abund_prok_otu %>% 
      filter(.data = ., !Slice %in% "Control")
  ) %>% 
  mutate(Pore = as.factor(Pore),
         Moisture = as.factor(Moisture))

df_abund_prok_otu_add

table(
  df_abund_prok_otu_add$Pore,
  df_abund_prok_otu_add$Treatment1
)

unique(df_abund_prok_otu_add$Class)


# Filtering dataset using Mean and CV ------------------------------------------
otu_prok_highcv <-
  df_abund_prok_otu_add %>% 
  select(Taxonomy, Class, Count) %>% 
  group_by(Taxonomy, Class) %>% 
  summarise_each(funs(mean, sd, cv)) %>% 
  ungroup() %>% 
  na.omit() %>% 
  filter(cv > quantile(cv, probs = c(0.75))) %>% 
  arrange(Class, desc(mean)) %>% 
  group_by(Class) %>% 
  top_n(n = 5, wt = mean) %>% 
  ungroup()

otu_prok_highcv$Class  

otu_prok_highcv %>% 
  subset(otu_prok_highcv$Class  %in% c("Root")) %>% 
  as.data.frame() 


otu_prok_highcv_top <-
  df_abund_prok_otu_add %>% 
  select(Taxonomy, Class, Count) %>% 
  group_by(Taxonomy, Class) %>% 
  summarise_each(funs(mean, sd, cv)) %>% 
  ungroup() %>% 
  na.omit() %>% 
  filter(cv > quantile(cv, probs = c(0.75))) %>% 
  arrange(desc(mean)) %>% 
  top_n(n = 25, wt = mean)

otu_prok_highcv_top$Class  


# >>>> PLOTTING FUNCTIONS <<<<< ------------------------------------------------
# Function to generate DF with mean and significance 
DiffAbundSet <- function(dataframe, Treat_lev, Slice_lev, Comp_var, Col_name, formula){
  require(tidyverse)
  Var <- enquo(Comp_var)
  df_mean <-
    dataframe %>% 
    subset(Treatment1 %in% c(Treat_lev, "Control") & Slice %in% c(Slice_lev, "Control")) %>% 
    subset(x=., Treatment %in% c(Treat_lev)) %>% 
    select(Taxonomy, Count, !!Var) %>% 
    group_by(Taxonomy, !!Var) %>%
    summarise_each(funs(sum, median, mean, sd)) %>% 
    select(Taxonomy, !!Var, mean) %>% 
    pivot_wider(names_from = Col_name, values_from = "mean") %>% 
    filter(!if_all(where(is.numeric), function(x) { x == 0})) %>% 
    filter(if_any(!contains("Control"), function(x) { x > 0})) %>% 
    pivot_longer(cols=-Taxonomy, names_to=Col_name, values_to = "Mean") %>% 
    arrange(.data = ., Taxonomy)
  
  # I also removed OTUs that have no statistical difference acorss groups, 
  # and those that show up only in the controls samples- see the filter 
  
  df_sig <-
    dataframe %>% 
    subset(Treatment1 %in% c(Treat_lev, "Control") & Slice %in% c(Slice_lev, "Control")) %>% 
    subset(x=., Treatment %in% c(Treat_lev)) %>% 
    select(Taxonomy, Count, !!Var) %>% 
    mutate(Taxonomy = ifelse(is.na(Taxonomy), paste("Unclassified"),paste(Taxonomy))) %>% 
    nest(data=-Taxonomy) %>% 
    mutate(Sig = map(.x=data, ~CompSampl(df = ., formula), data = .x),
           Group = map(.x=Sig, ~rownames(x = .), data=.x)) %>%
    unnest(c(Sig, Group)) %>% 
    select(Taxonomy, Group, Letters) %>% 
    pivot_wider(names_from = Group, values_from = Letters) %>% 
    filter(!if_all(!contains("Taxonomy"), function(x) { x == "a"})) %>% 
    pivot_longer(cols=-Taxonomy, names_to=Col_name, values_to = "Sig") %>% 
    arrange(.data = ., !!Var)
  
  df_final <-
    right_join(
      df_sig %>% 
        unite(df_ord, c(Col_name, "Taxonomy"), remove = FALSE),
      df_mean %>% 
        unite(df_ord, c(Col_name, "Taxonomy"), remove = FALSE),
      by="df_ord") %>% 
    select(2,3,4,7) %>% 
    mutate(Niche = rep(Col_name, times= nrow(x=.))) %>% 
    rename_with(.cols=1:2, ~ c("Taxonomy", "Group"))
  
  return(df_final)
}


# test
df_abund_fungi_otu_add %>% 
  filter(Taxonomy %in% otu_fungi_highcv$Taxonomy) %>% 
  DiffAbundSet(dataframe = ., "T1", "Leaf", Pore, "Pore", formula(Count ~ Pore)) %>% 
  mutate(Mean = ifelse(Mean == 0, NA, Mean)) %>%
    ggplot(aes(x = Group, y = Taxonomy, size=Mean)) +
    geom_point(shape = 21, colour="black", stroke=0.5) +
    geom_point(aes(fill=Group),shape = 21, alpha = 0.7) 

df_abund_prok_otu_add %>% 
  filter(Taxonomy %in% otu_fungi_highcv$Taxonomy) %>% 
  DiffAbundSet(dataframe = ., "T1", "Leaf", Pore, "Pore", formula(Count ~ Pore)) %>% 
  mutate(Mean = ifelse(Mean == 0, NA, Mean)) 




# Generate dataframe for plotting ----------------------------------------------
DftoPlot<-function(dataframe, Treat_lev, Slice_lev, df_cv, taxonomy_db){
  require(tidyverse)
  
  df_slice <- rbind(
    dataframe %>% 
    filter(Taxonomy %in% df_cv$Taxonomy) %>% 
    DiffAbundSet(dataframe = ., Treat_lev, Slice_lev, Pore, "Pore", formula(Count ~ Pore)) %>% 
    mutate(RF = ifelse(Taxonomy %in% sel_df_fungi_tax$Taxonomy, Taxonomy, NA)) %>% 
    left_join(y = taxonomy_db %>%select(Taxonomy, Class),by = "Taxonomy")%>% 
    arrange(Taxonomy) %>% 
    drop_na(Taxonomy),
    dataframe %>% 
    filter(Taxonomy %in% df_cv$Taxonomy) %>% 
    DiffAbundSet(dataframe = ., Treat_lev, Slice_lev, Plant, "Plant", formula(Count ~ Plant)) %>% 
    mutate(RF = ifelse(Taxonomy %in% sel_df_fungi_tax$Taxonomy, Taxonomy, NA)) %>% 
    left_join(y = taxonomy_db %>%select(Taxonomy, Class),by = "Taxonomy")%>% 
    arrange(Taxonomy) %>% 
    drop_na(Taxonomy),
    dataframe %>% 
    filter(Taxonomy %in% df_cv$Taxonomy) %>% 
    DiffAbundSet(dataframe = ., Treat_lev, Slice_lev, Moisture, "Moisture", formula(Count ~ Moisture)) %>% 
    mutate(RF = ifelse(Taxonomy %in% sel_df_fungi_tax$Taxonomy, Taxonomy, NA)) %>% 
    left_join(y = taxonomy_db %>%select(Taxonomy, Class),by = "Taxonomy")%>% 
    arrange(Taxonomy) %>% 
    drop_na(Taxonomy))
  
  #test4 <-
  #  dataframe %>% 
  #  filter(Taxonomy %in% df_cv$Taxonomy) %>% 
  #  DiffAbundSet(dataframe = ., Treat_lev, Slice_lev, Time, "Time", formula(Count ~ Time)) %>% 
  #  mutate(RF = ifelse(Taxonomy %in% sel_df_fungi_tax$Taxonomy, Taxonomy, NA)) %>% 
  #  left_join(y = taxonomy_db %>%select(Taxonomy, Class),by = "Taxonomy")%>% 
  #  arrange(Taxonomy) %>% 
  #  drop_na(Taxonomy)
  
  #df_slice <- rbind(test1, test2, test3) %>% 
  #  mutate(Mean = ifelse(Mean == 0, NA, Mean))
  
  # match the controls to the groups OTUs to reduce the list
  #df_all <- 
  #  rbind(df_slice, 
  #        test4 %>% 
  #          filter(Taxonomy %in% df_slice$Taxonomy)) 
  # print Range for bubble size 
  df_slice %>% 
    mutate(Mean = ifelse(Mean == 0, NA, Mean)) %>% 
    pull(Mean) %>% range(na.rm = TRUE) %>% print()
  
  df_plot <-
    df_slice %>%
    mutate(Mean = ifelse(Mean == 0, NA, Mean)) %>% 
    mutate(Class = ifelse(is.na(Class), "Unclassified", Class)) %>% 
    mutate(Niche = fct_relevel(Niche, 
          "Plant","Pore","Moisture")) %>% 
    mutate(Group = as.factor(Group),
           Group = recode(Group,
           Control = "C",
           CBP = "Low", FC = "High", 
           Control_CBP = "C-Low", Control_FC = "C-High",
           Control_Large = "C-Large", Control_Small = "C-Small"), 
           Group = fct_relevel(Group, 
           "Small", "Large", "C-Large","C-Small",
           "Low", "High","C-Low","C-High",
           "Corn", "Soy", "C"))
  return(df_plot)
}


bbl_fungi_T1_leaf <-
  df_abund_fungi_otu_add %>% 
  DftoPlot(dataframe = , "T1", "Leaf", otu_fungi_highcv_top, taxonomy_fungi)
bbl_fungi_T4_leaf <-
  df_abund_fungi_otu_add %>% 
  DftoPlot(dataframe = , "T4", "Leaf", otu_fungi_highcv_top, taxonomy_fungi)


bbl_fungi_T1_Near <-
  df_abund_fungi_otu_add %>% 
  DftoPlot(dataframe = , "T1", "Nearsoil",otu_fungi_highcv_top, taxonomy_fungi)
bbl_fungi_T4_Near <-
  df_abund_fungi_otu_add %>% 
  DftoPlot(dataframe = , "T4", "Nearsoil",otu_fungi_highcv_top, taxonomy_fungi)

bbl_fungi_T1_Far <-
  df_abund_fungi_otu_add %>% 
  DftoPlot(dataframe = , "T1", "Farsoil",otu_fungi_highcv_top, taxonomy_fungi)
bbl_fungi_T4_Far <-
  df_abund_fungi_otu_add %>% 
  DftoPlot(dataframe = , "T4", "Farsoil",otu_fungi_highcv_top, taxonomy_fungi)




bbl_prok_T1_leaf <-
  df_abund_prok_otu_add %>% 
  DftoPlot(dataframe = , "T1", "Leaf",otu_prok_highcv_top, taxonomy_prok)
bbl_prok_T4_leaf <-
  df_abund_prok_otu_add %>% 
  DftoPlot(dataframe = , "T4", "Leaf",otu_prok_highcv_top, taxonomy_prok)


bbl_prok_T1_Near <-
  df_abund_prok_otu_add %>% 
  DftoPlot(dataframe = , "T1", "Nearsoil",otu_prok_highcv_top, taxonomy_prok)
bbl_prok_T4_Near <-
  df_abund_prok_otu_add %>% 
  DftoPlot(dataframe = , "T4", "Nearsoil",otu_prok_highcv_top, taxonomy_prok)

bbl_prok_T1_Far <-
  df_abund_prok_otu_add %>% 
  DftoPlot(dataframe = , "T1", "Farsoil",otu_prok_highcv_top, taxonomy_prok)
bbl_prok_T4_Far <-
  df_abund_prok_otu_add %>% 
  DftoPlot(dataframe = , "T4", "Farsoil",otu_prok_highcv_top, taxonomy_prok)



# Plotting OTU with differences in abundance -----------------------------------
PlotBubble <- function(dataframe){
  require(forcats)
  bubble_plot <-
  dataframe %>% 
    mutate(Class = as.factor(Class)) %>% 
    mutate(Taxonomy = as.factor(Taxonomy)) %>% 
    mutate(ordering = -as.numeric(Class),
           Taxonomy = fct_reorder(Taxonomy, ordering, .desc = T)) %>% 
    ggplot(aes(x = Group, y = Taxonomy, size=Mean)) +
    geom_point(shape = 21, colour="black", stroke=0.5) +
    geom_point(aes(fill=Class), shape = 21, alpha = 0.7) +
    facet_grid(~Niche, scales = "free_x", space="free_x") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, size = 8, hjust = 1, vjust = 1, 
                                 face = "bold", color = "black"),
      axis.text.y = element_text(angle = 0, size = 7, hjust = 1, vjust = 0.5, color = "black"),
      axis.title = element_blank(),
      strip.text = element_text(face="bold"),
      legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm")) +
    guides(fill = guide_legend(title.theme = element_text(size=10, face="bold"),
                  order = 1,ncol = 1,title = "Class",
                  override.aes = list(shape = 22, size = 4),
           shape = guide_legend(title.theme = element_text(size=10, face="bold"),
                  order = 3, ncol = 1,
                  override.aes = list(color = "black"))
    ))
  return(bubble_plot)
}



# Pooling FUNGI bbl in the same graph ans use faceting -------------------------
combined_fungi_bbl <-
  rbind(
    bbl_fungi_T1_leaf %>%
      mutate(Set = rep("T1-Leaf", times= nrow(x=.))),
    bbl_fungi_T4_leaf %>%
      mutate(Set = rep("T4-Leaf", times= nrow(x=.))),
   bbl_fungi_T1_Near %>%
    mutate(Set = rep("T1-Nearsoil", times= nrow(x=.))),
   bbl_fungi_T4_Near %>%
    mutate(Set = rep("T4-Nearsoil", times= nrow(x=.))),
   bbl_fungi_T1_Far %>%
    mutate(Set = rep("T1-Farsoil", times= nrow(x=.))),
  bbl_fungi_T4_Far %>%
    mutate(Set = rep("T4-Farsoil", times= nrow(x=.)))
)

# range across all Sets
range(combined_fungi_bbl$Mean, na.rm = TRUE)

# Fix taxonomy through BLAST
leaf_otu_list <-
combined_fungi_bbl_fix %>% 
  subset(Set %in% c("T1-Leaf", "T4-Leaf")) %>% 
  separate(Taxonomy, c("OTU_ID", "Taxon"), sep = "-", remove = FALSE) %>% 
  mutate(OTU_ID = gsub("FOTU", "FOTU_", OTU_ID)) %>% 
  select(OTU_ID) %>% unique()

refseq(physeq_fungi_ctr_ev)[leaf_otu_list$OTU_ID,] 

write.dna(refseq(refseq(physeq_fungi_ctr_ev)[leaf_otu_list$OTU_ID,]), 
          format="fasta", 
          colsep="",
          file="fungi_leaf.fasta")


combined_fungi_bbl_fix <-
  combined_fungi_bbl %>% 
  mutate(Taxonomy = as.character(Taxonomy)) %>% 
  mutate(Taxonomy = ifelse(Taxonomy %in%c("FOTU58-Lasiosphaeriaceae"), paste("FOTU58-Apodus"), paste(Taxonomy)),
         Taxonomy = ifelse(Taxonomy %in%c("FOTU100-Lasiosphaeriaceae"), paste("FOTU100-Podospora"), paste(Taxonomy)),
         Taxonomy = ifelse(Taxonomy %in%c("FOTU615-Auriculariales"), paste("FOTU615-Exidiopsis"), paste(Taxonomy)),
         Taxonomy = ifelse(Taxonomy %in%c("FOTU157-Nectriaceae"), paste("FOTU157-Fusarium"), paste(Taxonomy)),
         Taxonomy = ifelse(Taxonomy %in%c("FOTU267-Sordariales"), paste("FOTU216-Trichoderma"), paste(Taxonomy)),
         Taxonomy = ifelse(Taxonomy %in%c("FOTU296-Lasiosphaeriaceae"), paste("FOTU296-Apisordaria"), paste(Taxonomy)),
         Taxonomy = ifelse(Taxonomy %in%c("FOTU527-Sordariales"), paste("FOTU527-Zopfiella"), paste(Taxonomy)),
         Taxonomy = ifelse(Taxonomy %in%c("FOTU368-Sordariales"), paste("FOTU368-Chaetomium"), paste(Taxonomy)),
         Taxonomy = ifelse(Taxonomy %in%c("FOTU366-Agaricostilbales"), paste("FOTU366-Ballistosporomyces"), paste(Taxonomy))
         )

Chaetomium
# Pooling BACTERIA bbl in the same graph ans use faceting ----------------------
combined_prok_bbl <-
  rbind(
    bbl_prok_T1_leaf %>%
      mutate(Set = rep("T1-Leaf", times= nrow(x=.))),
    bbl_prok_T4_leaf %>%
      mutate(Set = rep("T4-Leaf", times= nrow(x=.))),
    bbl_prok_T1_Near %>%
      mutate(Set = rep("T1-Nearsoil", times= nrow(x=.))),
    bbl_prok_T4_Near %>%
      mutate(Set = rep("T4-Nearsoil", times= nrow(x=.))),
    bbl_prok_T1_Far %>%
      mutate(Set = rep("T1-Farsoil", times= nrow(x=.))),
    bbl_prok_T4_Far %>%
      mutate(Set = rep("T4-Farsoil", times= nrow(x=.)))
  )

# range across all Sets
range(combined_prok_bbl$Mean, na.rm = TRUE)

# Fix taxonomy through BLAST
combined_prok_bbl_fix <-
  combined_prok_bbl %>% 
  mutate(Taxonomy = as.character(Taxonomy)) %>% 
  mutate(Taxonomy = ifelse(Taxonomy %in%c("POTU6122−unclassified Gammaproteobacteria"), paste("POTU6122-Gammaproteobacteria"), paste(Taxonomy))
  )


# ***** FINAL FIGURES ----------------------------------------------------------

# ***** LEAF plot ---------------------------------------------------------------

combined_fungi_bbl %>% 
  subset(Set %in% c("T1-Leaf", "T4-Leaf")) %>% 
  pull(Mean) %>% 
  range(na.rm = TRUE)

combined_prok_bbl %>% 
  subset(Set %in% c("T1-Leaf", "T4-Leaf")) %>% 
  pull(Mean) %>% 
  range(na.rm = TRUE)


p_leaf_fungi <-
  combined_fungi_bbl_fix %>% 
  subset(Set %in% c("T1-Leaf", "T4-Leaf")) %>% 
  PlotBubble(dataframe = .) +
  #facet_wrap(~Set+Niche, nrow = 1, scales = "free_x") +
  facet_grid(~Set+Niche, scales = "free_x", space="free_x") +
  labs(title="Fungal OTU Abundance", x="Niche", y="Average OTU abundance") +
  scale_fill_manual(values = palette_CB30) +
  theme(strip.background = element_blank(), 
        legend.title = element_text(size=10, face="bold", color = "black")) +
  scale_size(name = "Mean\nAbundance",
             range = c(0.1, 10),
             limits = c(0, NA), 
             breaks = c(0.1, 2.5, 5.0, 15, 30),
             labels = c(0.1, 2.5, 5.0, 15, 30))

p_leaf_fungi

p_leaf_bact <-
combined_prok_bbl %>% 
  subset(Set %in% c("T1-Leaf", "T4-Leaf")) %>% 
  PlotBubble(dataframe = .) +
  #facet_wrap(~Set+Niche, nrow = 1, scales = "free_x") +
  facet_grid(~Set+Niche, scales = "free_x", space="free_x") +
  labs(title="Bacterial OTU Abundance", x="Niche", y="Average OTU abundance") +
  scale_fill_manual("", values = palette_CB30) +
  theme(strip.background = element_blank(),
        legend.title = element_text(size=10, face="bold", color = "black")) +
  scale_size(name = "Mean\nAbundance",
             range = c(0.1, 10),
             limits = c(0, NA), 
             breaks = c(0.1, 2.5, 5.0, 10, 15),
             labels = c(0.1, 2.5, 5.0, 10, 15))
p_leaf_bact


Leaf_plot <- 
ggarrange(
  p_leaf_fungi, p_leaf_bact,
  ncol = 1,
  nrow = 2,
  labels = c("A", "B"), 
  heights = c(0.46, 1.2),
  align = "hv",
  common.legend = FALSE
)

Leaf_plot



# ****** LEAF plot with top OTUS -----------------------------------------------
combined_fungi_bbl %>% 
  subset(Set %in% c("T1-Leaf", "T4-Leaf")) %>% 
  pull(Mean) %>% 
  range(na.rm = TRUE)

combined_prok_bbl %>% 
  subset(Set %in% c("T1-Leaf", "T4-Leaf")) %>% 
  pull(Mean) %>% 
  range(na.rm = TRUE)


p_leaf_fungi_top <-
  combined_fungi_bbl_fix %>% 
  subset(Set %in% c("T1-Leaf", "T4-Leaf")) %>% 
  PlotBubble(dataframe = .) +
  #facet_wrap(~Set+Niche, nrow = 1, scales = "free_x") +
  facet_grid(~Set+Niche, scales = "free_x", space="free_x") +
  labs(title="Fungal OTU Abundance", x="Niche", y="Average OTU abundance") +
  scale_fill_manual(values = palette_CB30) +
  theme(strip.background = element_blank(), 
        legend.title = element_text(size=10, face="bold", color = "black")) +
  scale_size(name = "Mean\nAbundance",
             range = c(0.1, 10),
             limits = c(0, NA), 
             breaks = c(0.1, 2.5, 5.0, 15, 30),
             labels = c(0.1, 2.5, 5.0, 15, 30))

p_leaf_fungi_top

p_leaf_bact_top <-
  combined_prok_bbl %>% 
  subset(Set %in% c("T1-Leaf", "T4-Leaf")) %>% 
  PlotBubble(dataframe = .) +
  #facet_wrap(~Set+Niche, nrow = 1, scales = "free_x") +
  facet_grid(~Set+Niche, scales = "free_x", space="free_x") +
  labs(title="Bacterial OTU Abundance", x="Niche", y="Average OTU abundance") +
  scale_fill_manual("", values = palette_CB30) +
  theme(strip.background = element_blank(),
        legend.title = element_text(size=10, face="bold", color = "black")) +
  scale_size(name = "Mean\nAbundance",
             range = c(0.1, 10),
             limits = c(0, NA), 
             breaks = c(0.1, 2.5, 5.0, 10, 15),
             labels = c(0.1, 2.5, 5.0, 10, 15))
p_leaf_bact_top


Leaf_plot_top <- 
  ggarrange(
    p_leaf_fungi_top, p_leaf_bact_top,
    ncol = 1,
    nrow = 2,
    labels = c("A", "B"), 
    heights = c(0.86, 1.05),
    align = "hv",
    common.legend = FALSE
  )

Leaf_plot_top

# ***** NEARSOIL plot ----------------------------------------------------------
combined_fungi_bbl %>% 
  subset(Set %in% c("T1-Nearsoil", "T4-Nearsoil")) %>% 
  pull(Mean) %>% 
  range(na.rm = TRUE)

combined_prok_bbl %>% 
  subset(Set %in% c("T1-Nearsoil", "T4-Nearsoil")) %>% 
  pull(Mean) %>% 
  range(na.rm = TRUE)


near_otu_list <-
  combined_fungi_bbl_fix %>% 
  subset(Set %in% c("T1-Nearsoil", "T4-Nearsoil")) %>% 
  separate(Taxonomy, c("OTU_ID", "Taxon"), sep = "-", remove = FALSE) %>% 
  mutate(OTU_ID = gsub("FOTU", "FOTU_", OTU_ID)) %>% 
  select(OTU_ID) %>% unique()

refseq(physeq_fungi_ctr_ev)[near_otu_list$OTU_ID,] 

write.dna(refseq(refseq(physeq_fungi_ctr_ev)[near_otu_list$OTU_ID,]), 
          format="fasta", 
          colsep="",
          file="fungi_near.fasta")


p_near_fungi <-
  combined_fungi_bbl_fix %>% 
  subset(Set %in% c("T1-Nearsoil", "T4-Nearsoil")) %>% 
  PlotBubble(dataframe = .) +
  #facet_wrap(~Set+Niche, nrow = 1, scales = "free_x") +
  facet_grid(~Set+Niche, scales = "free_x", space="free_x") +
  labs(title="Fungal OTU Abundance", x="Niche", y="Average OTU abundance") +
  scale_fill_manual(values = palette_CB30) +
  theme(strip.background = element_blank(), 
        legend.title = element_text(size=10, face="bold", color = "black")) +
  scale_size(name = "Mean\nAbundance",
             range = c(0.5, 15),
             limits = c(0, NA), 
             breaks = c(0.1, 10, 30, 70),
             labels = c(0.1, 10, 30, 70))

p_near_fungi

p_near_bact <-
  combined_prok_bbl %>% 
  subset(Set %in% c("T1-Nearsoil", "T4-Nearsoil")) %>% 
  PlotBubble(dataframe = .) +
  #facet_wrap(~Set+Niche, nrow = 1, scales = "free_x") +
  facet_grid(~Set+Niche, scales = "free_x", space="free_x") +
  labs(title="Bacterial OTU Abundance", x="Niche", y="Average OTU abundance") +
  scale_fill_manual("", values = palette_CB30) +
  theme(strip.background = element_blank(),
        legend.title = element_text(size=10, face="bold", color = "black")) +
  scale_size(name = "Mean\nAbundance",
             range = c(0.1, 10),
             limits = c(0, NA), 
             breaks = c(0.1, 1, 2.5),
             labels = c(0.1, 1, 2.5))
p_near_bact


Nearsoil_plot <- 
  ggarrange(
    p_near_fungi, p_near_bact,
    ncol = 1,
    nrow = 2,
    labels = c("A", "B"), 
    heights = c(1.05, 1),
    align = "hv",
    common.legend = FALSE
  )

Nearsoil_plot


# ***** FARSOIL plot -----------------------------------------------------------
combined_fungi_bbl %>% 
  subset(Set %in% c("T1-Farsoil", "T4-Farsoil")) %>% 
  pull(Mean) %>% 
  range(na.rm = TRUE)

combined_prok_bbl %>% 
  subset(Set %in% c("T1-Farsoil", "T4-Farsoil")) %>% 
  pull(Mean) %>% 
  range(na.rm = TRUE)


far_otu_list <-
  combined_fungi_bbl_fix %>% 
  subset(Set %in% c("T1-Farsoil", "T4-Farsoil")) %>% 
  separate(Taxonomy, c("OTU_ID", "Taxon"), sep = "-", remove = FALSE) %>% 
  mutate(OTU_ID = gsub("FOTU", "FOTU_", OTU_ID)) %>% 
  select(OTU_ID) %>% unique()

refseq(physeq_fungi_ctr_ev)[far_otu_list$OTU_ID,] 

write.dna(refseq(refseq(physeq_fungi_ctr_ev)[far_otu_list$OTU_ID,]), 
          format="fasta", 
          colsep="",
          file="fungi_far.fasta")


p_far_fungi <-
  combined_fungi_bbl_fix %>% 
  subset(Set %in% c("T1-Farsoil", "T4-Farsoil")) %>% 
  PlotBubble(dataframe = .) +
  #facet_wrap(~Set+Niche, nrow = 1, scales = "free_x") +
  facet_grid(~Set+Niche, scales = "free_x", space="free_x") +
  labs(title="Fungal OTU Abundance", x="Niche", y="Average OTU abundance") +
  scale_fill_manual(values = palette_CB30) +
  theme(strip.background = element_blank(), 
        legend.title = element_text(size=10, face="bold", color = "black")) +
  scale_size(name = "Mean\nAbundance",
             range = c(0.5, 15),
             limits = c(0, NA), 
             breaks = c(0.1, 10, 30, 70),
             labels = c(0.1, 10, 30, 70))

p_far_fungi

p_far_bact <-
  combined_prok_bbl %>% 
  subset(Set %in% c("T1-Farsoil", "T4-Farsoil")) %>% 
  PlotBubble(dataframe = .) +
  #facet_wrap(~Set+Niche, nrow = 1, scales = "free_x") +
  facet_grid(~Set+Niche, scales = "free_x", space="free_x") +
  labs(title="Bacterial OTU Abundance", x="Niche", y="Average OTU abundance") +
  scale_fill_manual("", values = palette_CB30) +
  theme(strip.background = element_blank(),
        legend.title = element_text(size=10, face="bold", color = "black")) +
  scale_size(name = "Mean\nAbundance",
             range = c(0.1, 10),
             limits = c(0, NA), 
             breaks = c(0.1, 2.5, 5, 10),
             labels = c(0.1, 2.5, 5, 10))
p_far_bact


Farsoil_plot <- 
  ggarrange(
    p_far_fungi, p_far_bact,
    ncol = 1,
    nrow = 2,
    labels = c("A", "B"), 
    heights = c(1, 0.75),
    align = "hv",
    common.legend = FALSE
  )

Farsoil_plot
