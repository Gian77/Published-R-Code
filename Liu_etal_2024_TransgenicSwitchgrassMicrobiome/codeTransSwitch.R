# **** TRANSGENIC SWITCHGRASS ****----------------------------------------------
# Manuscript:   transgenic switchgrass
# Authors:      Shuang, Ming, Nico, Greg
# Affiliation:  Michigan State University
# Journal:      ... 
# Date:         Feb, 21 2024
# ***********************************************-------------------------------

options(
  digits = 4,
  pillar.sigfig = 7,
  scipen = 9999,
  max.print = 10000000
)

### SETTING UP WORKING ENVIRONEMNT ---------------------------------------------
library(styler)
library(tidyverse)
library(tidysq)
library(data.table)
library(vegan)
library(MASS)
library(ggpubr)
library(magrittr)
library(ggrepel)
library(ggtext)
library(scales)
library(grid)
library(gridExtra)

# COLOR PLAETTES ---------------------------------------------------------------
paletteCB2 = c("#CC2D35","#848FA2")

paletteCB6 = c("#CC2D35","#599861","#FF934F","#848FA2","#825121") #"#058ED9"

paletteCB15 = c("#560d0d","#dba4a4", "#cc1c1c","#111b77","#058ED9","#014443","#117744","#60ffaf",
                "#825121","#fcb067","#fcfc00","#82807f", "#5b5b19","#ffb7ef","#ae09ea","#000000")

paletteCB30 = c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77","#283dff","#058ED9","#636bb7",
                         "#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
                         "#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#5b5b19",
                         "#fcfc00","#ffff9e","#ffb7ef","#fa7efc","#ae09ea","#521899","#000000")

par(mfrow = c(1, 3))
pie(rep(1, length(paletteCB6)), labels = sprintf("%d (%s)",
    seq_along(paletteCB6),paletteCB6), col = paletteCB6, main="paletteCB6")
pie(rep(1, length(paletteCB2)), labels = sprintf("%d (%s)",
    seq_along(paletteCB2),paletteCB2), col = paletteCB2, main="paletteCB2")
pie(rep(1, length(paletteCB30)), labels = sprintf("%d (%s)",
    seq_along(paletteCB30),paletteCB30), col = paletteCB30, main="paletteCB30")
dev.off()

# CUSTOM TEHEMES ---------------------------------------------------------------
themePlain <- function(){ 
  require(ggtext)
  theme_bw() +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          strip.text = element_text(angle=90, size = 10, face = "bold", hjust = 0, vjust = 0.5),
          legend.title = element_blank(), legend.position = "right",
          strip.background = element_blank())
}


themeStats <- function(){ 
  require(ggtext)
  theme_bw() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 10, face = "bold"),
          axis.text.x = element_markdown(angle = 45, size = 8, hjust = 1, vjust = 1.05),
          axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
          axis.title = element_blank(),
          strip.background = element_blank())
}

themePlot <- function(){ 
  require(ggtext)
  theme_bw() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 10, face = "bold"),
          axis.text.x = element_markdown(angle = 45, size = 8, hjust = 1, vjust = 1.05),
          axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
          axis.title = element_blank(),
          strip.background = element_blank())
}

# ***********************************************-------------------------------
# Import fungal datasets -------------------------------------------------------

fungi_table <- data.frame(data.table::fread("Original_MingYi-fungi/PR_ITS_asv_new.txt", sep="\t"), row.names=1)
dim(fungi_table)
fungi_table[1:6, 1:50]

fungi_metadata <- as.data.frame(
  read.delim("Original_MingYi-fungi/Metadata_Transgenic_correct.txt", header=TRUE, row.names=1, sep="\t"))
dim(fungi_metadata)
fungi_metadata

fungi_taxonomy <- as.data.frame(
  read.delim("Original_MingYi-fungi/constax_taxonomy_ITS.txt", header=TRUE, row.names=1, sep="\t"))
dim(fungi_taxonomy)
fungi_taxonomy[1:6,]

fungi_DNAseq <- read_fasta("Original_MingYi-fungi/ITS_ASVs_filtered_new.fasta") %>% 
  dplyr::select("ASV_ID"=2, "DNA_seq"=1)
dim(fungi_DNAseq)

fungi_table_filt <-
  fungi_table %>% 
  rownames_to_column("Sample_ID") %>% 
  mutate(Sample_ID = gsub("ITS", "", Sample_ID)) %>% 
  column_to_rownames("Sample_ID") %>% 
  t(x = .) %>%
  as.data.frame() %>%
  rownames_to_column("ASV_ID") %>% 
  mutate(ASV_ID = gsub("ASV", "", ASV_ID)) %>% 
  mutate(ASV_ID = gsub("^0+", "", ASV_ID)) %>% 
  mutate(ASV_ID = paste("ASV_", ASV_ID, sep = "")) %>% 
  as_tibble()

fungi_table_filt[1:20, 1:10]
dim(fungi_table_filt)

fungi_metadata_filt <-
fungi_metadata %>%
  rownames_to_column("Sample_ID") %>% 
  mutate(Sample_ID = gsub("ITS", "", Sample_ID)) %>% 
  as_tibble()

fungi_metadata_filt

fungi_taxonomy_filt <-
  fungi_taxonomy %>%
  rownames_to_column("ASV_ID") %>%
  mutate(
    Kingdom = gsub("_1", "", Kingdom),
    Phylum = gsub("_1", "", Phylum),
    Class = gsub("_1", "", Class),
    Order = gsub("_1", "", Order),
    Family = gsub("_1", "", Family),
    Genus = gsub("_1", "", Genus),
    Species = gsub("_1", "", Species)
  ) %>%
  as_tibble()

fungi_taxonomy_filt
fungi_taxonomy_filt$ASV_ID


# Checking the High_level_taxonomy ---------------------------------------------
taxonomy_check_fungi <-
  fungi_taxonomy_filt %>%
  dplyr::select(ASV_ID, Kingdom, High_level_taxonomy) %>%
  left_join(x = .,
            y = fungi_table_filt, by = "ASV_ID") %>%
  pivot_longer(
    cols = c(-ASV_ID, -High_level_taxonomy, -Kingdom),
    names_to = "Sample_ID",
    values_to = "Counts"
  ) %>%
  left_join(y = fungi_metadata_filt,
            by = "Sample_ID") %>%
  mutate(
    High_level_taxonomy = ifelse(High_level_taxonomy == "", "Unlcassified", High_level_taxonomy),
    High_level_taxonomy = ifelse(High_level_taxonomy == "Eukaryota_kgd_Incertae_sedis", "Unlcassified", High_level_taxonomy),
    Kingdom = ifelse(Kingdom == "", "Unlcassified", Kingdom)
  )

taxonomy_check_fungi

anyNA(taxonomy_check_fungi)
anyNA(taxonomy_check_fungi$Counts)
levels(as.factor(taxonomy_check_fungi$High_level_taxonomy))
levels(as.factor(taxonomy_check_fungi$Kingdom))
anyNA(taxonomy_check_fungi$High_level_taxonomy)

taxonomy_plot_fungi_high <-
  taxonomy_check_fungi %>%
  subset(!SampleRound%in%c("AtTransplant")) %>% 
    dplyr::select(High_level_taxonomy, Sample_ID, SampleType, SampleRound, Counts) %>% 
    group_by(High_level_taxonomy, SampleType, SampleRound, Sample_ID) %>% 
    summarise(across(c(Counts), list(mean = mean, sum = sum))) %>% 
    ungroup() %>% 
      ggplot(aes(x = Sample_ID, y = Counts_sum, fill = High_level_taxonomy)) +
      geom_bar(stat = "identity") +
      facet_grid( ~ SampleRound ~ SampleType, scales = "free_x", space = "free_x") +
      themePlain() +
      scale_fill_manual(values = c(Alveolata="#dba4a4",
                                   Apusozoa="#521899",
                                   Fungi="#cc1c1c",
                                   Metazoa="#b7ffdb",
                                   Protista="#825121",
                                   Rhizaria="#fcb067",
                                   Stramenopila="#d8d6d4",
                                   Unlcassified="#000000",
                                   Viridiplantae="#014443")) +
    labs(title = "High_level_taxonomy", x=NULL, y="Reads number")

taxonomy_plot_fungi_high

taxonomy_plot_fungi_king <-
taxonomy_check_fungi %>%
    dplyr::select(Kingdom, Sample_ID, SampleType, SampleRound, Counts) %>% 
    group_by(Kingdom, SampleType, SampleRound, Sample_ID) %>% 
    summarise(across(c(Counts), list(mean = mean, sum = sum))) %>% 
    ungroup() %>% 
      ggplot(aes(x = Sample_ID, y = Counts_sum, fill = Kingdom)) +
      geom_bar(stat = "identity") +
      facet_grid( ~ SampleRound ~ SampleType, scales = "free_x", space = "free_x") +
      themePlain() +
      scale_fill_manual(values = c(Alveolata="#dba4a4",
                                   Amoebozoa="#111b77",
                                   Anthophyta="#bfc5ff",
                                   Apusozoa="#521899",
                                   Fungi="#cc1c1c",
                                   Metazoa="#b7ffdb",
                                   Protista="#825121",
                                   Rhizaria="#fcb067",
                                   Rhodoplantae="#ffff9e",
                                   Stramenopila="#d8d6d4",
                                   Unlcassified="#000000",
                                   Viridiplantae="#014443")) +
  labs(title = "Kingdom", x=NULL, y="Reads number")

taxonomy_plot_fungi_king


# Detect untarget taxa to remove, looking at both the taxonomy at Kingdom level
# and the High_level_taxonomy assigned by constax2
anyNA(taxonomy_check_fungi)
anyNA(taxonomy_check_fungi$Counts)
anyNA(taxonomy_check_fungi$High_level_taxonomy)
anyNA(taxonomy_check_fungi$Kingdom)

levels(as.factor(taxonomy_check_fungi$High_level_taxonomy))
levels(as.factor(taxonomy_check_fungi$Kingdom))

table(fungi_taxonomy_filt$High_level_taxonomy) %>% 
  
  as.data.frame() %>% 
  rename("Taxon"=Var1)

table(fungi_taxonomy_filt$Kingdom) %>% 
  as.data.frame() %>% 
  rename("Taxon"=1, "Freq"=2)

table(fungi_taxonomy_filt$Phylum) %>% 
  as.data.frame() %>% 
  rename("Taxon"=1, "Freq"=2)


untarget_fungi <- c("Amoebozoa","Metazoa", "Protista", "Rhizaria", "Viridiplantae", "Stramenopila",
                     "Apusozoa", "Rhodoplantae", "Anthophyta", "Alveolata",
                     "Ichthyosporia", "Choanoflagellozoa", "Eukaryota_kgd_Incertae_sedis")

fungi_taxonomy_new <-
  fungi_taxonomy_filt %>%
  subset(!High_level_taxonomy %in% untarget_fungi &
           !Kingdom %in% untarget_fungi) %>%
  subset(Kingdom %in% c("Fungi") &
           HL_hit_query_cover > 60 |
           HL_hit_percent_id > 60) %>%
  subset(!Kingdom %in% c("")) %>%
  filter(Kingdom %in% "Fungi") %>% 
  subset(!Phylum %in% c("Anthophyta", "Cercozoa"))

fungi_taxonomy_new

levels(as.factor(fungi_taxonomy_new$Kingdom))
levels(as.factor(fungi_taxonomy_new$Phylum))
fungi_taxonomy_new[fungi_taxonomy_new$Phylum=="Cercozoa", ]
fungi_taxonomy_new[fungi_taxonomy_new$Phylum=="Anthophyta", ]

# extarcting ASV_ID for all Glomeromycota --------------------------------------
Glomero_ASVs <-
left_join(fungi_taxonomy_filt,
          fungi_table_filt,
          by = "ASV_ID") %>% 
  filter(Phylum %in% "Glomeromycota", .preserve = FALSE) %>% 
  dplyr::select(ASV_ID, starts_with("TR")) %>% 
  pivot_longer(cols = -ASV_ID, names_to = "Sample_ID", values_to ="Counts") %>% 
  left_join(x=.,
            y= fungi_metadata_filt,
            by = "Sample_ID")

Glomero_ASVs

Glomero_ASVs %>% 
    ggplot(aes(x = Sample_ID, y = Counts, fill = "red")) +
    geom_bar(stat = "identity") +
    facet_grid( ~ SampleRound ~ SampleType, scales = "free_x", space = "free_x") +
    themePlain() +
  theme(legend.position = "none")

# Filter out Glomero in the leaf and inflorescence
fungi_no_glo <-
full_join(
  fungi_table_filt %>% 
    column_to_rownames("ASV_ID") %>% 
    t(x=.) %>% 
    as.data.frame() %>% 
    rownames_to_column("Sample_ID"),
  fungi_metadata_filt,
  by = "Sample_ID") %>% 
  as_tibble() %>% 
  pivot_longer(cols = c(-Sample_ID, -SampleType, 
                        -SampleRound, -Treatment), 
               names_to = "ASV_ID", values_to = "Count") %>% 
subset( ! (SampleType%in%c("Inflorescence") & ASV_ID%in%Glomero_ASVs$ASV_ID)) %>% 
subset( ! (SampleType%in%c("Leaf") & ASV_ID%in%Glomero_ASVs$ASV_ID)) %>% 
mutate_all(., ~replace_na(.,0))

fungi_no_glo

# Test again after removing the Glomero
left_join(fungi_no_glo, 
          fungi_taxonomy_filt, by = "ASV_ID") %>% 
  filter(Phylum %in% "Glomeromycota", .preserve = FALSE) %>% 
  ggplot(aes(x = Sample_ID, y = Count, fill = "red")) +
  geom_bar(stat = "identity") +
  facet_grid( ~ SampleRound ~ SampleType, scales = "free_x", space = "free_x") +
  themePlain()


# Regenerate the dataframe
fungi_table_filt_noGlo <-
fungi_no_glo %>% 
  dplyr::select(ASV_ID, Sample_ID, Count) %>% 
  tidyr::pivot_wider(id_cols = ASV_ID, names_from = Sample_ID, values_from = Count) %>% 
  mutate_all(., ~replace_na(.,0))

fungi_table_filt_noGlo %>% as.data.frame()


# Plot BLAST hits % coverage and % identity before and after the filtering
par(mfrow = c(2, 2))
hist(
  fungi_taxonomy_filt$HL_hit_percent_id,
  breaks = 50,
  xlim = c(0, 100),
  ylim = c(0, 5000),
  xlab = "% Identity",
  main = "Sequence identity - Pre"
)
hist(
  fungi_taxonomy_new$HL_hit_percent_id,
  breaks = 10,
  xlim = c(0, 100),
  ylim = c(0, 5000),
  xlab = "% Identity",
  main = "Sequence identity - Post"
)

hist(
  fungi_taxonomy_filt$HL_hit_query_cover,
  breaks = 50,
  xlim = c(0, 100),
  ylim = c(0, 5000),
  xlab = "% Query covearge",
  main = "Query coverage - Pre"
)
hist(
  fungi_taxonomy_new$HL_hit_query_cover,
  breaks = 10,
  xlim = c(0, 100),
  ylim = c(0, 5000),
  xlab = "% Query covearge",
  main = "Query coverage - Post"
)
dev.off()

# Double check that only fungi are in the dataset
table(fungi_taxonomy_new$Kingdom)
table(fungi_taxonomy_new$High_level_taxonomy)

# Removing specific ASV contaminants of the leaves
fungi_taxonomy_new <-
  fungi_taxonomy_new %>% 
  filter(!ASV_ID%in%c("ASV_1267","ASV_1065","ASV_536"))

# Import bacteria datasets -----------------------------------------------------
bact_table <- data.frame(data.table::fread("Original_Ming-Yi_bacteria/PR_16S_asv.txt", sep=","), row.names=1)
dim(bact_table)
head(bact_table)

bact_metadata <- read.csv("Original_Ming-Yi_bacteria/Metadata_Transgenic_bacteria.csv", header=T, row.names=1)
dim(bact_metadata)

bact_taxonomy <- read.delim("Original_Ming-Yi_bacteria/constax_taxonomy.txt", header=TRUE, row.names=1, sep="\t")
dim(bact_taxonomy)

bact_table_filt <-
  bact_table %>%
  rownames_to_column("Sample_ID") %>% 
  mutate(Sample_ID = gsub("16S", "", Sample_ID)) %>% 
  column_to_rownames("Sample_ID") %>% 
  t(x = .) %>%
  as_tibble() %>%
  rownames_to_column("ASV_ID") %>%
  mutate(ASV_ID = paste("ASV_", ASV_ID, sep = "")) %>% 
  as_tibble()

bact_table_filt[1:20, 1:20]
bact_table_filt

bact_metadata_filt <-
  bact_metadata %>%
  rownames_to_column("Sample_ID") %>% 
  mutate(Sample_ID = gsub("16S", "", Sample_ID)) %>% 
  as_tibble()

bact_metadata_filt

bact_DNAseq <-
  bact_table %>%
  t(x = .) %>%
  as.data.frame() %>%
  rownames_to_column("DNAseq") %>%
  as_tibble() %>%
  dplyr::select(DNAseq) %>%
  rownames_to_column("ASV_ID") %>%
  mutate(ASV_ID = paste("ASV_", ASV_ID, sep = ""))

bact_DNAseq

# removing Rank_8 that does not have any classification
unique(bact_taxonomy$Rank_7)
unique(bact_taxonomy$Rank_8)

bact_taxonomy_filt <-
  bact_taxonomy %>%
  dplyr::select(-Rank_8) %>%
  rownames_to_column("ASV_ID") %>%
  rename(
    Kingdom = Rank_1,
    Phylum = Rank_2,
    Class = Rank_3,
    Order = Rank_4,
    Family = Rank_5,
    Genus = Rank_6,
    Species = Rank_7
  ) %>%
  mutate(
    Kingdom = gsub("_1", "", Kingdom),
    Phylum = gsub("_1", "", Phylum),
    Class = gsub("_1", "", Class),
    Order = gsub("_1", "", Order),
    Family = gsub("_1", "", Family),
    Genus = gsub("_1", "", Genus),
    Species = gsub("_1", "", Species)
  ) %>%
  as_tibble()

bact_taxonomy_filt 

# Checking the High_level_taxonomy ---------------------------------------------
# and Testing for Chloroplast and Mitochondria contamination
taxonomy_check_bact <-
  bact_taxonomy_filt %>%
  dplyr::select(ASV_ID, High_level_taxonomy) %>%
  left_join(bact_table_filt, 
            by = "ASV_ID") %>%
  pivot_longer(
    cols = c(-ASV_ID,-High_level_taxonomy),
    names_to = "Sample_ID",
    values_to = "Counts"
  ) %>%
  left_join(y = bact_metadata_filt,
            by = "Sample_ID") %>%
  mutate(
    High_level_taxonomy = as.factor(High_level_taxonomy),
    High_level_taxonomy = gsub("_.*", "", High_level_taxonomy)
  ) %>%
  mutate(
    High_level_taxonomy = ifelse(High_level_taxonomy == "", "Unlcassified", High_level_taxonomy)
  ) 

taxonomy_check_bact
levels(as.factor(taxonomy_check_bact$High_level_taxonomy))

taxonomy_plot_bact <-
taxonomy_check_bact %>% 
  subset(!SampleRound%in%c("AtTransplant")) %>% 
  dplyr::select(High_level_taxonomy, Sample_ID, SampleType, SampleRound, Counts) %>% 
  group_by(High_level_taxonomy, SampleType, SampleRound, Sample_ID) %>% 
  summarise(across(c(Counts), list(mean = mean, sum = sum))) %>% 
    ggplot(aes(x = Sample_ID, y = Counts_sum, fill = High_level_taxonomy)) +
    geom_bar(stat = "identity") +
    facet_grid( ~ SampleRound ~ SampleType, scales = "free_x", space = "free_x") +
    themePlain() +
  scale_fill_manual(values = paletteCB6)

taxonomy_plot_bact
  
# Detect untarget taxa to remove, looking at both the taxonomy at Kingdom level
# and the High_level_taxonomy assigned by constax2
dim(bact_taxonomy_filt)
table(bact_taxonomy_filt$High_level_taxonomy)
table(bact_taxonomy_filt$Kingdom)

table(bact_taxonomy_filt$High_level_taxonomy) %>% 
  as.data.frame() %>% 
  rename("Taxon"=1, "Freq"=2)

table(bact_taxonomy_filt$Kingdom) %>% 
  as.data.frame() %>% 
  rename("Taxon"=1, "Freq"=2)

table(bact_taxonomy_filt$Phylum) %>% 
  as.data.frame() %>% 
  rename("Taxon"=1, "Freq"=2)

untarget_bact <- c("Mitochondria","Chloroplast")

bact_taxonomy_new <-
  bact_taxonomy_filt %>%
  subset(!Order %in% untarget_bact &
           !Family %in% untarget_bact) %>%
  subset(!High_level_taxonomy %in% untarget_bact &
           !High_level_taxonomy %in% untarget_bact) %>%
  subset(HL_hit_query_cover > 60 |
           HL_hit_percent_id > 60) %>%
  subset(!Kingdom %in% c(""))

bact_taxonomy_new

# Checking the Cyanobacteria
bact_taxonomy_new %>% 
  subset(Phylum%in%c("Cyanobacteria")) %>% 
  group_by(Order) %>% 
  count()

# Plot BLAST hits % coverage and % identity before and after the filtering
par(mfrow = c(2, 2))
hist(
  bact_taxonomy_filt$HL_hit_percent_id,
  breaks = 50,
  xlim = c(0, 100),
  ylim = c(0, 40000),
  xlab = "% Identity",
  main = "Sequence identity - Pre"
)
hist(
  bact_taxonomy_new$HL_hit_percent_id,
  breaks = 10,
  xlim = c(0, 100),
  ylim = c(0, 40000),
  xlab = "% Identity",
  main = "Sequence identity - Post"
)

hist(
  bact_taxonomy_filt$HL_hit_query_cover,
  breaks = 50,
  xlim = c(0, 100),
  ylim = c(0, 40000),
  xlab = "% Query covearge",
  main = "Query coverage - Pre"
)
hist(
  bact_taxonomy_new$HL_hit_query_cover,
  breaks = 10,
  xlim = c(0, 100),
  ylim = c(0, 40000),
  xlab = "% Query covearge",
  main = "Query coverage - Post"
)
dev.off()

# Double check that only bact are in the dataset
table(bact_taxonomy_new$Kingdom)
table(bact_taxonomy_new$High_level_taxonomy) %>% as.data.frame()
any(bact_taxonomy_new$Order == "Chloroplast")
any(bact_taxonomy_new$Family == "Mitochondria")

# ***********************************************-------------------------------
# Combine fungal dataset for GitHub -------------------------------------------
save(fungi_table, fungi_metadata, fungi_taxonomy, fungi_DNAseq, file = "fungi_data.Rdata")
save(bact_table, bact_metadata, bact_taxonomy, bact_DNAseq, file = "bacteria_data.Rdata")

load("fungi_data.Rdata")
load("bacteria_data.Rdata")

# ***********************************************-------------------------------
# REFORMAT METADATA ------------------------------------------------------------
# match fungal and bacterial metadata

# Generate dataset with correct factor order
FixData <- function(dataframe){
    new_dataframe <-
    dataframe %>% 
    rename("Niche"=SampleType, "Status"=SampleRound)  %>% 
    mutate(Niche = recode(Niche, `Senescence Leaf`="Senescent Leaf", BulkSoil="Soil"),
           Niche = fct_relevel(Niche, 
                               "Inflorescence","Leaf", "Senescent Leaf", 
                               "Root", "Rhizosphere","Soil"),
           Status=recode(Status, 
                         AtTransplant="At-Transplant",
                         PreTransplant="Pre-Transplant",
                         PostTransplant="Post-Transplant"),
           Status=fct_relevel(Status, 
                              "At-Transplant","Pre-Transplant","Post-Transplant"),
           Treatment = recode(Treatment, 
                        Trans = "GMO", None = "non-GMO"),
           Treatment = fct_relevel(Treatment, "GMO","non-GMO")) %>% 
    as_tibble()
  return(new_dataframe)
}


# removing ITS flag from samples name and remove samples "At-Transplant". Also 
# add a variable for easy splitting the Niche factor. 
# I also remove the Pre-Transplant Leaf and Senescent Leaf
meta_fungi_fix <-
  FixData(fungi_metadata_filt) %>% 
  mutate(Source = recode(Niche, Inflorescence ="IN",  Leaf= "LF", `Senescent Leaf` = "SL", 
                         Root= "RT", Rhizosphere = "RI", Soil ="SO")) %>% 
  filter(!Status %in% c("At-Transplant")) %>% 
  filter(!Status == "Pre-Transplant" | !Niche == "Leaf") %>% 
  filter(!Niche == "Senescent Leaf") %>% 
  droplevels()

meta_fungi_fix

meta_bact_fix <-
  FixData(bact_metadata_filt %>% 
            mutate(SampleType = recode(SampleType, SenescenceLeaf = "Senescence Leaf"))) %>% 
  mutate(Source = recode(Niche, Inflorescence ="IN",  Leaf= "LF", `Senescent Leaf` = "SL", 
                         Root= "RT", Rhizosphere = "RI", Soil ="SO")) %>% 
  filter(!Status %in% c("At-Transplant")) %>% 
  filter(!Status == "Pre-Transplant" | !Niche == "Leaf") %>% 
  filter(!Niche == "Senescent Leaf") %>% 
  droplevels()

meta_bact_fix

# REFORMAT TAXONOMY ------------------------------------------------------------
ReformatTaxonomy <- function(taxonomy_tab){
  require(tidyverse)
  lastValue <- function(x) tail(x[!is.na(x)], 1)
  last_taxons<- apply(taxonomy_tab[,c(2:8)], 1, lastValue)
  taxonomy_tab$BestMatch <- last_taxons
  taxonomy_res <-
    taxonomy_tab %>% 
    unite(ASV_ID, BestMatch, col=Taxonomy, sep = " ", remove = FALSE)
  return(taxonomy_res)
}

# You taxonomy table has to have a column with OTU_ID as first column
# to make this function work properly. Also, you need to convert any 
#"","NA","na","N/A","n/a","NaN","nan" to NA before running it.

# We use na_if() within mutate to convert annoying annoying value to NA.

fungi_taxonomy_fix <-
  fungi_taxonomy_new[,1:8] %>% 
  mutate(across(where(is.character), ~na_if(., ""))) %>% 
  ReformatTaxonomy() 

fungi_taxonomy_fix

# gsub using regexp .* = any number of any character
dim(bact_taxonomy_new)
bact_taxonomy_fix <-
  bact_taxonomy_new[,1:8] %>% 
  mutate(across(where(is.character),~gsub("uncultured.*", "", x=.)),
         across(where(is.character),~gsub("Uncultured.*", "", x=.))) %>% 
  mutate(across(where(is.character), ~na_if(., ""))) %>% 
  ReformatTaxonomy() 

bact_taxonomy_fix

# REFROMAT DATA ----------------------------------------------------------------
# removing the ITS flag in the names, transpose to columns as ASVs, and remove 
# samples "At-Transplant"

identical(meta_fungi_fix$Sample_ID, meta_bact_fix$Sample_ID)
sample_order <- match(meta_bact_fix$Sample_ID, meta_fungi_fix$Sample_ID)
meta_fungi_fix <- meta_fungi_fix[sample_order, ]
# metadata names are now identical!

fungi_table_mod <-
  fungi_table_filt_noGlo %>% 
  dplyr::select(ASV_ID, meta_fungi_fix$Sample_ID) %>% # match the metadata
  column_to_rownames("ASV_ID") %>% 
  filter(rowSums(x=.) > 0) %>% # remove ASV that sum up to 0
  dplyr::select(where( ~is.numeric(.x) && sum(.x) != 0)) %>% # remove samples that sum up to 0
  rownames_to_column("ASV_ID") %>%
  as_tibble()

colnames(fungi_table_mod)
fungi_table_mod

bact_table_mod <-
  bact_table_filt %>% 
  dplyr::select(ASV_ID, meta_bact_fix$Sample_ID) %>% 
  column_to_rownames("ASV_ID") %>% 
  filter(rowSums(x=.) > 0) %>% 
  dplyr::select(where( ~is.numeric(.x) && sum(.x) != 0)) %>% 
  rownames_to_column("ASV_ID") %>%
  as_tibble()

colnames(bact_table_mod)
bact_table_mod

identical(colnames(fungi_table_mod), colnames(bact_table_mod))
# TRUE

# MATCH TABLE, METADATA, TAXONOMY ----------------------------------------------
dim(fungi_table_mod)
dim(fungi_taxonomy_fix)
dim(meta_fungi_fix)

dim(bact_table_mod)
dim(bact_taxonomy_fix)
dim(meta_bact_fix)

# Remember that the number of ASVs is affected by the filtering the samples 
#"At-Transplant" in the otu_table  
fungi_table_fix <-
fungi_table_mod %>% 
  filter(ASV_ID %in% fungi_taxonomy_fix$ASV_ID) %>% 
  column_to_rownames("ASV_ID") %>% 
  filter(rowSums(x=.) > 0) %>% 
  dplyr::select(where( ~is.numeric(.x) && sum(.x) != 0)) %>% 
  rownames_to_column("ASV_ID") %>%
  as_tibble()

dim(fungi_table_fix)

bact_table_fix <-
  bact_table_mod %>% 
  filter(ASV_ID %in% bact_taxonomy_fix$ASV_ID) %>% 
  column_to_rownames("ASV_ID") %>% 
  filter(rowSums(x=.) > 0) %>% 
  dplyr::select(where( ~is.numeric(.x) && sum(.x) != 0)) %>% 
  rownames_to_column("ASV_ID") %>%
  as_tibble()

dim(bact_table_fix)

# compare fungi and bacteria metadata
identical(colnames(fungi_table_fix), colnames(bact_table_fix))

# TRUE
colnames(fungi_table_fix); dim(fungi_table_fix)
colnames(bact_table_fix); dim(bact_table_fix)

fungi_table_fix %>% 
  column_to_rownames("ASV_ID") %>% colSums()

bact_table_fix %>% 
  column_to_rownames("ASV_ID") %>% colSums() 

# ***********************************************-------------------------------
# FIND RAREFACTION DEPTH -------------------------------------------------------

# Calculating Good's coverage. The fraction of sequences that appear in an OTU
# that have been seen more than once, and allows estimating what percent of the total 
# species is represented in a sample.
# Coverage = 1 - (number of individuals in species / total number of individuals)
# Example: If I have Goods = 0.96 it means that 4% of your reads in that sample are
# from OTUs that appear only once in that sample.

# Calculate long dataframe with stats 
RareStats <- function(dataframe){
# calculate distribution outliers
findoutlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5*IQR(x) | x > quantile(x, 0.75) + 1.5*IQR(x))
}
# generate output dataframe
df_output <-
  dataframe %>% 
  pivot_longer(-ASV_ID, names_to = "Sample_ID", values_to = "Seq_No") %>% 
  group_by(Sample_ID) %>% 
  summarize(Read_No = sum(Seq_No),
            Singlton_No = sum(Seq_No == 1),
            Goods = 100*(1 - Singlton_No / Read_No)) %>% 
  mutate(outlier = ifelse(findoutlier(log10(Read_No)), Read_No, NA)) %>% 
  ungroup()
return(df_output)
}


# Plotting fucntion 
PloRareStats <- function(dataframe){
  require(ggrepel)
  plot_output <-
    ggarrange(
      dataframe %>%
        ggplot(aes(x = Read_No)) +
        geom_histogram(binwidth = 5000, 
                    fill = "firebrick", color = "firebrick") +
        themeStats() +
        labs(title = "Histogram"), 
      dataframe %>%
        ggplot(aes(x = Read_No)) +
        geom_histogram(binwidth = 1000, 
                     fill = "firebrick",color = "firebrick") +
        coord_cartesian(xlim = c(0, 25000)) +
        themeStats() +
        labs(title = "Histogram Zoom"), 
      dataframe %>%
        ggplot(aes(x = Read_No, y = Goods)) +
        geom_point(shape = 1, color = "firebrick") +
        themeStats() +
        labs(title = "Good's Coverage"), 
      dataframe %>%
        ggplot(aes(x = 1, y = Read_No)) +
        geom_jitter(shape=1, color = "firebrick") +
        scale_y_log10() +
        themeStats() +
        labs(title = "Log10 jitter"), 
      dataframe %>%
        ggplot(aes(x = 1, y = Read_No)) +
        geom_boxplot(color = "firebrick") +
        scale_y_log10() +
        geom_text_repel(
          data = filter(dataframe, !is.na(outlier)),
          mapping = aes(x = 1, y = Read_No, label = outlier), 
          max.overlaps = 15, size = 3) +
        themeStats() +
        labs(title = "Log10 boxplot"), 
      dataframe %>%
        arrange(Read_No) %>%
        ggplot(aes(x = 1:nrow(.), y = Read_No)) +
        geom_bar(stat = "identity", color = "firebrick") +
        themeStats() +
        labs(title = "Ranked"), 
    ncol = 3,
    nrow = 2, 
    align = "hv", 
    labels = c("A", "B", "C", "D", "E", "F"))
  return(plot_output)
}

rare_fungi <-
  RareStats(fungi_table_fix) %>% 
  arrange(Read_No) %>% 
  left_join(meta_fungi_fix, by = "Sample_ID")

rare_fungi

rare_fungi_plot <-
  RareStats(fungi_table_fix) %>% PloRareStats()
rare_fungi_plot


rare_bact <-
  RareStats(bact_table_fix) %>% 
  arrange(Read_No) %>% 
  left_join(meta_bact_fix, by = "Sample_ID")

rare_bact

rare_bact_plot <-
  RareStats(bact_table_fix) %>% PloRareStats()
rare_bact_plot

# ***********************************************-------------------------------
# RAREFACTION CURVES -----------------------------------------------------------
# Creating data frames ---------------------------------------------------------
# Generating data frames with samples that match the rarefaction depth
# This is necessary in vegan because samples for samples that have < the
# depth cutoff chosen no rarefaction is applied.

# Selected depth cutoffs
min_depth_fungi = 3090
min_depth_bact = 12126

#minDepth_fungi <-
#  fungi_table_fix %>% 
#  column_to_rownames("ASV_ID") %>% 
#  t(x=.) %>% 
#  as.data.frame()

# Fungi
df_filt_fungi <-
  fungi_table_fix %>% 
  pivot_longer(-ASV_ID, names_to = "Sample_ID", values_to = "Seq_No") %>% 
  group_by(Sample_ID) %>%
  mutate(Read_No = sum(Seq_No)) %>% 
  filter(Read_No >= min_depth_fungi) %>% 
  ungroup() %>%
  dplyr::select(-Read_No) %>% 
  pivot_wider(names_from ="ASV_ID", values_from = "Seq_No", values_fill = 0) %>% 
  column_to_rownames("Sample_ID") 

dim(df_filt_fungi)
df_filt_fungi %>%  as_tibble()

df_filt_fungi %>% 
  rownames_to_column("Sample_ID") %>% 
  left_join(y = meta_fungi_fix, by="Sample_ID") %>% 
  count(Niche)

# Bacteria
df_filt_bact <-
  bact_table_fix %>% 
  pivot_longer(-ASV_ID, names_to = "Sample_ID", values_to = "Seq_No") %>% 
  group_by(Sample_ID) %>%
  mutate(Read_No = sum(Seq_No)) %>% 
  filter(Read_No >= min_depth_bact) %>% 
  ungroup() %>%
  dplyr::select(-Read_No) %>% 
  pivot_wider(names_from ="ASV_ID", values_from = "Seq_No", values_fill = 0) %>% 
  column_to_rownames("Sample_ID")

dim(df_filt_bact)
df_filt_bact %>%  as_tibble()

df_filt_bact %>% 
  rownames_to_column("Sample_ID") %>% 
  left_join(y = meta_bact_fix, by="Sample_ID") %>% 
  count(Niche)

# Calculate probabilities that OTU occur in a rarefied community of size sample.
# Not sure if this is useful and I need to find a way to plot this.

drare_fungi <-
  df_filt_fungi %>% 
  drarefy(sample=min_depth_fungi) %>%
  as_tibble(rownames="Sample_ID") %>%
  pivot_longer(-Sample_ID)

drare_fungi

drare_fungi %>%
  ggplot(aes(x = value)) +
  geom_histogram(binwidth = 0.1, 
                 fill = "firebrick", color = "firebrick") +
  themeStats() +
  labs(title = "Histogram")

drare_bact <-
  df_filt_bact %>% 
  drarefy(sample=min_depth_bact) %>%
  as_tibble(rownames="Sample_ID") %>%
  pivot_longer(-Sample_ID)

drare_bact

# Extract the number of OTUs found only for specific set of OTU number. For example, 
# the the 1, 1001, 2001, 3001, ... sequence 
# The smaller is step the more precise the curve will be.

rarecurve_fungi <- 
  rarecurve(df_filt_fungi,  step=1000)

str(rarecurve_fungi)
as.data.frame(rarecurve_fungi[[1]])
rarecurve_fungi[[2]]

rarecurve_bact <- 
  rarecurve(df_filt_bact,  step=1000)

ExtrRare <- function(rare_obj, dataframe, metadata){
  require(tidyverse)
  # The dataframes and metadata must a column with Sample_ID
  # example on how to check using ifelse
  if ("Sample_ID" %in% colnames(metadata)) {
    cat("\nSample_ID variable is present. Good to go! \n")
  } else {
    cat("\nSample_ID variable is absent. Creating it from rownames. \n")
    metadata$Sample_ID <- rownames(metadata)
  }
  # Using purr::map_dfr and bind rows to generate a dataframe for
  # a list of named numbers and values 
  rare_joint <-
    map_dfr(rare_obj, bind_rows) %>%
    bind_cols(Sample_ID = rownames(dataframe), .) %>%
    pivot_longer(-Sample_ID) %>%
    drop_na() %>%
    mutate(n_seqs = as.numeric(str_replace(name, "N", ""))) %>%
    dplyr::select(-name) %>%
    left_join(x = .,
              y = metadata,
              by = "Sample_ID")
  return(rare_joint)
}

df_rarecurve_fungi <-
  ExtrRare(rarecurve_fungi, df_filt_fungi, meta_fungi_fix)

df_rarecurve_bact <-
  ExtrRare(rarecurve_bact, df_filt_bact, meta_bact_fix)

# Different way to extract the cureves
ExtrRareBase <- function(x){
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

ExtrRareBase(rarecurve_fungi)

# Then plotting rarecurve in ggplot2
PlotRareCurve <- function(rare_joint, depth){
  plot_rare <-
    ggplot(rare_joint, aes(x=n_seqs, y=value, group=Sample_ID, color=Niche)) +  
    geom_line() +
    geom_vline(xintercept = depth, color="black", lineTreatment = "dashed") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 30, size = 8, hjust = 1, vjust = 1, color = "black"),
      axis.text.y = element_text(angle = 0, size = 8, hjust = 1, vjust = 0.5, color = "black"),
      axis.title = element_text(angle = 0, size = 10, hjust = 0.5, vjust = 0.5),
      strip.text = element_text(face="bold"), strip.background = element_blank(),
      legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm")) +
    guides(color = guide_legend(ncol = 3, override.aes = list(shape = 15, size = 3)),
           shape = guide_legend(ncol = 1, override.aes = list(color = "black", size=2.5)))
  return(plot_rare)
}


PlotRareCurve(df_rarecurve_fungi, min_depth_fungi) +
  labs(title="ITS", x= "Number of DNA reads", y= "Number of OTUs")


# Plotting the graphs made with the two methods
plot_curves <-
  ggarrange(
    PlotRareCurve(df_rarecurve_fungi, min_depth_fungi) +
      scale_color_manual(values = paletteCB6) +
      labs(title = "Fungal Microbiome", x = "Number of DNA reads", y = "Number of OTUs"),
    PlotRareCurve(df_rarecurve_bact, min_depth_bact) +
      scale_color_manual(values = paletteCB6) +
      labs(title = "Bacterial Microbiome", x = "Number of DNA reads", y = "Number of OTUs"),
    common.legend = TRUE,
    legend = "bottom"
  )

plot_curves

# ***********************************************-------------------------------
# RAREFY TO EVEN DEPTH ---------------------------------------------------------

# 1) Calculating rarefied distances using avgdist() - as in Mothur.
rare_dist_fungi <-
fungi_table_fix %>% 
  pivot_longer(-ASV_ID, names_to = "Sample_ID", values_to = "Seq_No") %>% 
  mutate(Sample_ID = gsub("ITS", "", Sample_ID)) %>% 
  pivot_wider(names_from ="ASV_ID", values_from = "Seq_No", values_fill = 0) %>% 
  column_to_rownames("Sample_ID") %>% 
  avgdist(x=., dmethod="bray", sample=min_depth_fungi, iterations = 100)

rare_dist_fungi
dim(as.matrix(rare_dist_fungi))

rare_dist_bact <-
  bact_table_fix %>% 
  pivot_longer(-ASV_ID, names_to = "Sample_ID", values_to = "Seq_No") %>% 
  mutate(Sample_ID = gsub("16S", "", Sample_ID)) %>% 
  pivot_wider(names_from ="ASV_ID", values_from = "Seq_No", values_fill = 0) %>% 
  column_to_rownames("Sample_ID") %>% 
  avgdist(x=., dmethod="bray", sample=min_depth_bact, iterations = 100)

rare_dist_bact
dim(as.matrix(rare_dist_bact))

# 2) Function to rarefy data 100 times and generate an average, based on what 
# is implemented in Mothur and avgdist but on the data frame not on the distance matrix
# that will be calculated after. In this way I can subset the data frame and the rerun 
# the analysis on a subset of the data. The standard format in vegan is a data frame with 
# samples as rows and taxa as columns. Warning messages: In mean.default(X[[i]], ...) : 
# argument is not numeric or logical: returning NA is caused by the fact that aggregate tries
# to generate the mean of the Sample_ID column. The new grouping column is by default called 
# Group.1 by the aggreagte() function and renamed again to Sample_ID.

rarefyData <- function(dataframe, depth_level){
  require(tidyverse)
  require(data.table)
  
  #df <-
  #dataframe %>% 
  #  column_to_rownames("ASV_ID") %>% 
  #  t(x=.) %>% 
  #  as.data.frame()
  
  com_iter <- vector(mode = "list", length =  100)
  for (i in seq_along(com_iter)) {
    com_iter[[i]] <- as.data.frame(
      vegan::rrarefy(dataframe, sample = depth_level)) %>% 
      rownames_to_column("Sample_ID")
  }
  
  mean_100 <- do.call(rbind, com_iter)
  mean_100 <- aggregate(mean_100, 
                        list(mean_100$Sample_ID), 
                        mean)
  
  mean_df <-
    mean_100 %>% 
    dplyr::select(-Sample_ID) %>% 
    rename(Sample_ID = "Group.1") %>% 
    column_to_rownames("Sample_ID") %>% 
    as.data.frame() %>% 
    filter(rowSums(x=.) >= depth_level) %>% 
    dplyr::select(where( ~is.numeric(.x) && sum(.x) != 0))
  
  print(mean_df %>% as_tibble())
  return(mean_df)
}


set.seed(010101)

#https://stackoverflow.com/questions/31465415/combine-multiple-data-frames-and-calculate-average
#https://stackoverflow.com/questions/59980978/get-the-mean-across-list-of-dataframes-by-rows

fungi_ev_100itr <-
  rarefyData(df_filt_fungi, min_depth_fungi)
dim(fungi_ev_100itr)
fungi_ev_100itr %>% as_tibble()

fungi_ev_100itr


# check if filtered correctly
any(fungi_ev_100itr %>% colSums() ==0)
any(fungi_ev_100itr %>% rowSums() < min_depth_fungi)

set.seed(010102)
bact_ev_100itr <-
  rarefyData(df_filt_bact, min_depth_bact)
dim(bact_ev_100itr)
bact_ev_100itr %>% as_tibble()

# check if filtered correctly
any(bact_ev_100itr %>% colSums() ==0)
any(bact_ev_100itr %>% rowSums() < min_depth_bact)


# rrarefy
set.seed(122)
nmds_fungi_100 <- metaMDS(vegdist(fungi_ev_100itr, method = "bray"), k = 2, trymax = 200,
                      autotransform = FALSE, engine = c("monoMDS"))
str(nmds_fungi_100)

# avgdist
set.seed(123)
nmds_fungi <- metaMDS(rare_dist_fungi, trymax = 200, 
                      autotransform = FALSE, engine = c("monoMDS"))
str(nmds_fungi)

# Plotting the ordinations 
PlotNMDS <- function(nmds, metadata){
  nmds_plot <-
    nmds$points %>%
    as_tibble(rownames = "Sample_ID") %>%
    inner_join(., metadata,
               by="Sample_ID") %>%
    rename(Axis.1=MDS1, Axis.2=MDS2) %>% 
    ggplot(aes(x=Axis.1, y=Axis.2, color=Niche, shape=Status)) +
    geom_point() +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5, color = "black"),
      axis.text.y = element_text(angle = 0, size = 8, hjust = 1, vjust = 0.5, color = "black"),
      axis.title = element_text(angle = 0, size = 10, hjust = 0.5, vjust = 0.5),
      strip.text = element_text(face="bold"),
      strip.background = element_blank(),
      legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm")) +
    guides(color = guide_legend(ncol = 3, override.aes = list(shape = 15, size = 3.5)),
           shape = guide_legend(ncol = 1, override.aes = list(color = "black", size=2.5)))
  return(nmds_plot)
}


PlotNMDS(nmds_fungi, meta_fungi_fix) +
  labs(title = "Fungal Microbiome")

# Comparing avgdist() and rrarefy with 100 iterations and making an average
# Fig. 2 Betadiv ---------------------------------------------------------------

plot_betadiv_all <-
  ggarrange(
    PlotNMDS(nmds_fungi, meta_fungi_fix) +
      labs(title = "Fungal Microbiome\navgdist()"),
    PlotNMDS(nmds_fungi_100, meta_fungi_fix) +
      labs(title = "Fungal Microbiome\nrrarefy() 100 iter"),
    common.legend = TRUE,
    legend = "bottom"
  )

plot_betadiv_all

# ***********************************************-------------------------------
# METADATA MATCH ---------------------------------------------------------------

meta_fungi_100itr <-
meta_fungi_fix %>% 
  filter(Sample_ID %in% rownames(fungi_ev_100itr)) 

meta_bact_100itr <-
meta_bact_fix %>% 
  filter(Sample_ID %in% rownames(bact_ev_100itr)) 


table(
  meta_fungi_100itr %>% 
    mutate(Group = paste(Status, Niche, sep="-")) %>% 
    pull(Group),
  meta_fungi_100itr$Status)

table(
  meta_fungi_100itr %>% 
    mutate(Group = paste(Status, Niche, sep="-")) %>% 
    pull(Group),
  meta_fungi_100itr$Treatment)


table(
  meta_bact_100itr %>% 
    mutate(Group = paste(Status, Niche, sep="-")) %>% 
    pull(Group),
  meta_bact_100itr$Status)

table(
  meta_bact_100itr %>% 
    mutate(Group = paste(Status, Niche, sep="-")) %>% 
    pull(Group),
  meta_bact_100itr$Treatment)

# ***********************************************-------------------------------
# BETA DIVERSITY ---------------------------------------------------------------
set.seed(010323)
nmds_bact_100 <- metaMDS(vegdist(bact_ev_100itr, method = "bray"), k = 2, trymax = 200, 
                          autotransform = FALSE, engine = c("monoMDS"))
str(nmds_bact_100)

# Plotting all the graphs
# Fig. S4 NMDS ---------------------------------------------------
plot_beta <-
  ggarrange(
    PlotNMDS(nmds_fungi_100, meta_fungi_fix) +
      scale_color_manual(values = paletteCB6) +
      labs(title = "Fungal Microbiome"),
    PlotNMDS(nmds_bact_100, meta_bact_fix) +
      scale_color_manual(values = paletteCB6) +
      labs(title = "Bacterial Microbiome"),
    common.legend = TRUE,
    legend = "bottom"
  )

plot_beta

# Stress plot ------------------------------------------------------------------
# Function stressplot draws a Shepard plot where ordination distances
# are plotted against community dissimilarities, and the fit is shown as a
# monotone step line. The correlation based on stress is R2 = 1−S2.
# The “fit-based R2” is the correlation between the fitted values θ(d) and
# ordination distances ˜d, or between the step line and the points. This
# should be linear even when the fit is strongly curved and is often known
# as the “linear fit”. These two correlations are both based on the residuals
# in the Shepard plot, but they differ in their null models. In linear fit, the
# null model is that all ordination distances are equal, and the fit is a flat
# horizontal line. This sounds sensible, but you need N − 1 dimensions for
# the null model of N points, and this null model is geometrically impossible 
# in the ordination space. The basic stress uses the null model where all
# observations are put in the same point, which is geometrically possible.
# https://www.mooreecology.com/uploads/2/4/2/1/24213970/vegantutor.pdf

stressplot(nmds_fungi_100)
stressplot(nmds_bact_100)

# Non-metric fit and Linear fit R2
1-nmds_fungi$stress^2
cor(stressplot(nmds_fungi)$yf, stressplot(nmds_fungi)$y)^2

# Create a tibble that contains the data from stressplot

# position of the annotation layer 
# Annotation position and hjust, vjust use to modify it
# xpos ypos        annotateText        hjust    vjust
# 1 -Inf -Inf Bottom Left (h0,v0)        0        0
# 2 -Inf  Inf    Top Left (h0,v1)        0        1
# 3  Inf -Inf  Bottom Right h1,v0        1        0
# 4  Inf  Inf     Top Right h1,v1        1        1


PlotStress <- function(nmds){
  df_stress <-tibble(Observed = stressplot(nmds)$x,
                distance = stressplot(nmds)$y,
                fitted = stressplot(nmds)$yf) %>% 
    pivot_longer(cols = c(distance, fitted), names_to = "Groups")
  
  nmds_R2 <- round(1-nmds$stress^2, digits = 3)
  linear_R2 <- round(cor(stressplot(nmds)$yf, 
                         stressplot(nmds)$y)^2, digits = 3)
  stress <- nmds$stress

  stress_plot <-
    df_stress %>% 
    ggplot(aes(x = Observed, y = value)) +
    geom_point(data = df_stress %>% 
                 filter(Groups == "distance"), shape=1, alpha=0.5) +
    geom_step(data = df_stress %>% 
                filter(Groups == "fitted"), col = "red", direction = "vh") +
    labs(x = "Observed Dissimilarity", y = "Ordination Distance") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5, color = "black"),
      axis.text.y = element_text(angle = 0, size = 8, hjust = 1, vjust = 0.5, color = "black"),
      axis.title = element_text(angle = 0, size = 10, hjust = 0.5, vjust = 0.5)) +
    annotate("text", -Inf, Inf, label = paste("`Non-metric fit:`","~italic(R)^2==", nmds_R2), parse=TRUE, 
             color = "black",size = 3, hjust = -0.1, vjust = 1.5) +
    annotate("text", -Inf, Inf, label = paste("`Linear fit:`", "~italic(R)^2==", linear_R2), parse=TRUE,
             color = "black", size = 3, hjust = -0.1, vjust = 3) +
    annotate("text", -Inf, Inf, label = paste("Stress:", stress), parse=TRUE,
           color = "black", size = 3, hjust = -0.1, vjust = 6)
    #annotate("text",  -Inf, Inf, label=paste("Non-metric fit", "~italic(x)~';'~italic(r)^2==", nmds_R2), parse=TRUE, hjust = -0.1, vjust = 1.2)
  return(stress_plot)
}

PlotStress(nmds_fungi_100)
PlotStress(nmds_bact_100) 

# Plotting all the graphs
# Fig. SX Shepard plot ---------------------------------------------------------
stressplot_all <-
  ggarrange(
    PlotStress(nmds_fungi_100) + labs(title = "Shepard diagram\nfungal Microbiome"),
    PlotStress(nmds_bact_100) + labs(title = "Shepard diagram\nbacterial Microbiome"),
    common.legend = TRUE,
    legend = "bottom"
  )
stressplot_all

# Goodness of fit stats --------------------------------------------------------
goodness(nmds_fungi_100)
goodness(nmds_bact_100)

# Plotting the ordinations 
PlotGoodness <- function(nmds, metadata){
  gof_df <-
    nmds$points %>%
    as_tibble(rownames = "Sample_ID") %>%
    inner_join(., metadata,
               by="Sample_ID") %>%
    rename(Axis.1=MDS1, Axis.2=MDS2) %>% 
    mutate(Goodness = goodness(nmds))
  
  print(gof_df)

  goodness_plot <-
    gof_df %>% 
    ggplot(aes(x=Axis.1, y=Axis.2)) +
    geom_point(shape=16) +
    geom_point(shape=1, size=gof_df$Goodness*500, color="red") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5, color = "black"),
      axis.text.y = element_text(angle = 0, size = 8, hjust = 1, vjust = 0.5, color = "black"),
      axis.title = element_text(angle = 0, size = 10, hjust = 0.5, vjust = 0.5))
  return(goodness_plot)
}


PlotGoodness(nmds_fungi_100, meta_fungi_100itr)
PlotGoodness(nmds_bact_100, meta_bact_100itr)


# Plotting all the graphs
# Fig SXX Goodnest of fit in the NMDS ------------------------------------------
good_fit <-
  ggarrange(
    PlotGoodness(nmds_fungi_100, meta_fungi_100itr) +
      labs(title = "Goodness of fit\nfungal Microbiome"),
    PlotGoodness(nmds_bact_100, meta_bact_100itr) +
      labs(title = "Goodness of fot\nacterial Microbiome"),
    common.legend = TRUE,
    legend = "bottom"
  )
good_fit

# ***********************************************-------------------------------
# BETA DIVERSITY ON SPLIT DATASETS ---------------------------------------------

# Plotting splitted data
PlotSplitNMDS <- function(nmds, metadata){
  nmds_plot <-
    nmds$points %>%
    as_tibble(rownames = "Sample_ID") %>%
    inner_join(., metadata,
               by="Sample_ID") %>%
    rename(Axis.1=MDS1, Axis.2=MDS2) %>% 
    ggplot(aes(x=Axis.1, y=Axis.2, color=Treatment, shape=Status)) +
    scale_color_manual(values = paletteCB2) +
    geom_point(size=2) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5, color = "black"),
      axis.text.y = element_text(angle = 0, size = 8, hjust = 1, vjust = 0.5, color = "black"),
      axis.title = element_text(angle = 0, size = 10, hjust = 0.5, vjust = 0.5),
      strip.text = element_text(face="bold"),
      strip.background = element_blank(),
      legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm")) +
    guides(color = guide_legend(ncol = 2, override.aes = list(shape = 15, size = 3.5)),
           shape = guide_legend(ncol = 2, override.aes = list(color = "black", size=2.5)))
  return(nmds_plot)
}



SplitNMDS <- function(dataframe, metadata, Variable){
  require(tidyverse)
  require(vegan)
  
  otu <-
    dataframe %>%
    rownames_to_column("Sample_ID") %>% 
    left_join(metadata, by = "Sample_ID") %>% 
    filter(Niche %in% c(Variable)) %>% 
    column_to_rownames("Sample_ID") %>% 
    dplyr::select(where( ~is.numeric(.x) && sum(.x) != 0)) 

  # nmds
  nmds_out <- metaMDS(vegdist(otu, method = "bray"), k = 2, trymax = 200,
                            autotransform = FALSE, engine = c("monoMDS"))
  
  print(nmds_out$stress)
  return(nmds_out)
}


set.seed(201)
SplitNMDS(fungi_ev_100itr, meta_fungi_100itr, "Soil")

PlotSplitNMDS(SplitNMDS(fungi_ev_100itr, meta_fungi_100itr, "Soil"),meta_fungi_100itr)
  

# Multiplot Fungal Microbiome
set.seed(401)
split_nmds_fungi <-
  ggarrange(
  PlotSplitNMDS(SplitNMDS(fungi_ev_100itr, meta_fungi_100itr, "Soil"),meta_fungi_100itr) +
    labs(title = "Soil"),
  PlotSplitNMDS(SplitNMDS(fungi_ev_100itr, meta_fungi_100itr, "Rhizosphere"),meta_fungi_100itr)+
    labs(title = "Rhizosphere"),
  PlotSplitNMDS(SplitNMDS(fungi_ev_100itr, meta_fungi_100itr, "Root"),meta_fungi_100itr)+
    labs(title = "Root"),
  PlotSplitNMDS(SplitNMDS(fungi_ev_100itr, meta_fungi_100itr, "Leaf"),meta_fungi_100itr)+
    labs(title = "Leaf"),
  PlotSplitNMDS(SplitNMDS(fungi_ev_100itr, meta_fungi_100itr, "Inflorescence"),meta_fungi_100itr)+
    labs(title = "Inflorescence"),
  common.legend = TRUE,
  ncol = 3,
  nrow = 2,
  legend = "bottom"
)

split_nmds_fungi

set.seed(402)
split_nmds_bact <-
ggarrange(
  PlotSplitNMDS(SplitNMDS(bact_ev_100itr, meta_bact_100itr, "Soil"),meta_bact_100itr) +
    labs(title = "Soil"),
  PlotSplitNMDS(SplitNMDS(bact_ev_100itr, meta_bact_100itr, "Rhizosphere"),meta_bact_100itr)+
    labs(title = "Rhizosphere"),
  PlotSplitNMDS(SplitNMDS(bact_ev_100itr, meta_bact_100itr, "Root"),meta_bact_100itr)+
    labs(title = "Root"),
  PlotSplitNMDS(SplitNMDS(bact_ev_100itr, meta_bact_100itr, "Leaf"),meta_bact_100itr)+
    labs(title = "Leaf"),
  PlotSplitNMDS(SplitNMDS(bact_ev_100itr, meta_bact_100itr, "Inflorescence"),meta_bact_100itr)+
    labs(title = "Inflorescence"),
  common.legend = TRUE,
  ncol = 3,
  nrow = 2,
  legend = "bottom"
)

split_nmds_bact

# Check the stress. 
# If stress is nearly zero, the fit is almost perfect. Usually,
# this means that the data set is too small for the requested analysis, and there
# may be several different solutions that are almost as perfect.

ggarrange(
PlotStress(SplitNMDS(fungi_ev_100itr, meta_fungi_100itr, "Soil"))+
  labs(title = "Soil"),
PlotStress(SplitNMDS(fungi_ev_100itr, meta_fungi_100itr, "Rhizosphere"))+
  labs(title = "Rhizosphere"),
PlotStress(SplitNMDS(fungi_ev_100itr, meta_fungi_100itr, "Root"))+
  labs(title = "Root"),
PlotStress(SplitNMDS(fungi_ev_100itr, meta_fungi_100itr, "Leaf"))+
  labs(title = "Leaf"),
PlotStress(SplitNMDS(fungi_ev_100itr, meta_fungi_100itr, "Inflorescence"))+
  labs(title = "Inflorescence"),
common.legend = TRUE,
ncol = 3,
nrow = 2,
legend = "bottom"
)


ggarrange(
  PlotStress(SplitNMDS(bact_ev_100itr, meta_bact_100itr, "Soil"))+
    labs(title = "Soil"),
  PlotStress(SplitNMDS(bact_ev_100itr, meta_bact_100itr, "Rhizosphere"))+
    labs(title = "Rhizosphere"),
  PlotStress(SplitNMDS(bact_ev_100itr, meta_bact_100itr, "Root"))+
    labs(title = "Root"),
  PlotStress(SplitNMDS(bact_ev_100itr, meta_bact_100itr, "Leaf"))+
    labs(title = "Leaf"),
  PlotStress(SplitNMDS(bact_ev_100itr, meta_bact_100itr, "Inflorescence"))+
    labs(title = "Inflorescence"),
  common.legend = TRUE,
  ncol = 3,
  nrow = 2,
  legend = "bottom"
)

# PCOA PLOTS -------------------------------------------------------------------

# Calculating PCoA in vegan ----------------------------------------------------
SplitPCOA <- function(dataframe, metadata, Variable){
  require(tidyverse)
  require(vegan)
  otu <-
    dataframe %>%
    rownames_to_column("Sample_ID") %>% 
    left_join(metadata, by = "Sample_ID") %>% 
    filter(Niche %in% c(Variable)) %>% 
    column_to_rownames("Sample_ID") %>% 
    dplyr::select(where( ~is.numeric(.x) && sum(.x) != 0)) 
  
  pcoa <-
    cmdscale(vegdist(otu, method = "bray"), eig=TRUE)
  pcoa_df <-
    as.data.frame(pcoa$points) %>% 
    rename("Axis.1"=1, "Axis.2"=2) %>% 
    rownames_to_column("Sample_ID") %>% 
    left_join(metadata,
              by= "Sample_ID")

  return(pcoa_df)
}


SplitPCOA(bact_ev_100itr, meta_bact_100itr, "Soil")


PlotSplitPCOA <- function(pcoa_df, Var){
  
  if (Var == "single_shape"){
    shapes = 17
  } else {
    shapes = c(16, 17)
  }
  
pcoa_plot <-
  pcoa_df %>% 
  ggplot(aes(x=Axis.1, y=Axis.2, color=Treatment, shape=Status)) +
  scale_color_manual(values = paletteCB2) +
  scale_shape_manual(values = shapes) +
  geom_point(size=2) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5, color = "black"),
    axis.text.y = element_text(angle = 0, size = 8, hjust = 1, vjust = 0.5, color = "black"),
    axis.title = element_text(angle = 0, size = 10, hjust = 0.5, vjust = 0.5),
    strip.text = element_text(face="bold"),
    strip.background = element_blank(),
    legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm")) +
  guides(color = guide_legend(ncol = 2, override.aes = list(shape = 15, size = 3.5)),
         shape = guide_legend(ncol = 2, override.aes = list(color = "black", size=2.5)))

return(pcoa_plot)
}

PlotSplitPCOA(SplitPCOA(bact_ev_100itr, meta_bact_100itr, "Leaf"), "single_shape")


split_pcoa_fungi <-
  ggarrange(
    PlotSplitPCOA(SplitPCOA(fungi_ev_100itr, meta_fungi_100itr, "Soil"), "double_shape") +
      labs(title = "Soil"),
    PlotSplitPCOA(SplitPCOA(fungi_ev_100itr, meta_fungi_100itr, "Rhizosphere"), "double_shape")+
      labs(title = "Rhizosphere"),
    PlotSplitPCOA(SplitPCOA(fungi_ev_100itr, meta_fungi_100itr, "Root"), "double_shape")+
      labs(title = "Root"),
    PlotSplitPCOA(SplitPCOA(fungi_ev_100itr, meta_fungi_100itr, "Leaf"), "single_shape")+
      labs(title = "Leaf"),
    PlotSplitPCOA(SplitPCOA(fungi_ev_100itr, meta_fungi_100itr, "Inflorescence"), "single_shape")+
      labs(title = "Inflorescence"),
    common.legend = TRUE,
    ncol = 3,
    nrow = 2,
    legend = "bottom"
  )

split_pcoa_fungi


split_pcoa_bact <-
  ggarrange(
    PlotSplitPCOA(SplitPCOA(bact_ev_100itr, meta_bact_100itr, "Soil")) +
      labs(title = "Soil"),
    PlotSplitPCOA(SplitPCOA(bact_ev_100itr, meta_bact_100itr, "Rhizosphere"))+
      labs(title = "Rhizosphere"),
    PlotSplitPCOA(SplitPCOA(bact_ev_100itr, meta_bact_100itr, "Root"))+
      labs(title = "Root"),
    PlotSplitPCOA(SplitPCOA(bact_ev_100itr, meta_bact_100itr, "Leaf"))+
      labs(title = "Leaf"),
    PlotSplitPCOA(SplitPCOA(bact_ev_100itr, meta_bact_100itr, "Inflorescence"))+
      labs(title = "Inflorescence"),
    common.legend = TRUE,
    ncol = 3,
    nrow = 2,
    legend = "bottom"
  )

split_pcoa_bact


# ***********************************************-------------------------------
# PERMANOVAs -------------------------------------------------------------------

# Check if data and metadata samples are in the same order. If not, you need 
# to match one to the other and reorder.
identical(rownames(fungi_ev_100itr), meta_fungi_100itr$Sample_ID)

#identical(rownames(fungi_ev_100itr), meta_fungi_100itr$Sample_ID)
#sample_order <- match(meta_fungi_100itr$Sample_ID, rownames(fungi_ev_100itr))
#fungi_ev_100itr <- fungi_ev_100itr[sample_order, ]
# metadata names are now identical!

identical(meta_bact_100itr$Sample_ID, rownames(bact_ev_100itr))
length(rownames(fungi_ev_100itr))
length(meta_fungi_100itr$Sample_ID)

# Adonis function
Adjadonis2 <- function(formula, metadata, stratum){
  require(tidyverse)
  require(vegan)
  
  df_adonis <-
      adonis2(formula, 
              metadata,
              method = "bray",
              starta = stratum,
              permutations = 999)

  df_adj <-
  cbind(df_adonis,
        as.data.frame(p.adjust(df_adonis$`Pr(>F)`, 
                               method = "BH")) %>% 
          rename("p.adj"=1))
  
  ad_list <-
    list(formula, df_adonis,df_adj)
  
  return(ad_list)
}


set.seed(01)
adonis_test <-
  Adjadonis2(formula(fungi_ev_100itr ~ Status + Niche * Treatment),
             meta_fungi_100itr,
             stratum = NULL)
adonis_test


# Blocking by Status
Adonis2All <- function(dataframe, metadata){
  set.seed(11)
  ad1 <- Adjadonis2(formula(dataframe ~ Status + Niche * Treatment),metadata,stratum = NULL)
  set.seed(12)
  ad2 <- Adjadonis2(formula(dataframe ~ Status + Niche * Treatment),metadata,stratum = "Status")
  set.seed(13)
  ad3 <- Adjadonis2(formula(dataframe ~ Status +  Treatment * Niche),metadata,stratum = NULL)
  set.seed(14)
  ad4 <- Adjadonis2(formula(dataframe ~ Status +  Treatment * Niche),metadata, stratum = "Status")
  # Blocking by Niche
  set.seed(15)
  ad5 <- Adjadonis2(formula(dataframe ~ Niche + Status * Treatment),metadata,stratum = NULL)
  set.seed(16)
  ad6 <- Adjadonis2(formula(dataframe ~ Niche + Status * Treatment),metadata, stratum = "Niche")
  set.seed(17)
  ad7 <- Adjadonis2(formula(dataframe ~ Niche +  Treatment * Status),metadata,stratum = NULL)
  set.seed(18)
  ad8 <- Adjadonis2(formula(dataframe ~ Niche +  Treatment * Status),metadata, stratum = "Niche")
  # Blocking by Treatment
  set.seed(19)
  ad9 <- Adjadonis2(formula(dataframe ~ Treatment + Status * Niche),metadata,stratum = NULL)
  set.seed(20)
  ad10 <- Adjadonis2(formula(dataframe ~ Treatment + Niche * Status),metadata,stratum = "Treatment")
  set.seed(21)
  ad11 <- Adjadonis2(formula(dataframe ~ Treatment + Status * Niche),metadata,stratum = NULL)
  set.seed(22)
  ad12 <- Adjadonis2(formula(dataframe ~ Treatment + Status * Niche),metadata, stratum = "Treatment")
  
  combine_ad <- function(ad) {
    ad_df <-
      cbind(ad[[3]],
            data.frame(rep(as.character(ad[[1]])[3],
                           times = nrow(ad[[3]]))) %>%
              rename(Model = 1))
    return(ad_df)}
  
  return(
    data.frame(
    rbind(
      combine_ad(ad1),combine_ad(ad2),combine_ad(ad3),combine_ad(ad4),
      combine_ad(ad5),combine_ad(ad6),combine_ad(ad7),combine_ad(ad8),
      combine_ad(ad9), combine_ad(ad10), combine_ad(ad11), combine_ad(ad12)),
    Stratum = c(rep("-", times=6),rep("Status", times=6),rep("-", times=6),rep("Status", times=6),
                rep("-", times=6),rep("Niche", times=6),rep("-", times=6),rep("Niche", times=6),
                rep("-", times=6),rep("Treatment", times=6),rep("-", times=6),rep("Treatment", times=6))))
  }


main_adonis2_fungi <- Adonis2All(fungi_ev_100itr, meta_fungi_100itr)
main_adonis2_fungi

main_adonis2_bact <- Adonis2All(bact_ev_100itr, meta_bact_100itr)
main_adonis2_bact


# BETA-DISPERSION ---------------------------------------------------------------
BetadispExtr <- function(dataframe, metadata, Var){
  require(tidyverse)
  require(vegan)
  
  metadata_mod <-
    metadata %>% 
    column_to_rownames("Sample_ID") %>% 
    as.data.frame()
  
  disp <-
    betadisper(
      vegan::vegdist(dataframe, method="bray"),
      metadata_mod[,Var], type = c("centroid"),bias.adjust = FALSE,
      sqrt.dist = FALSE)
  anova_d <-
    anova(disp,
          permutations = how(nperm=9999))
  p_adj <-
    round(p.adjust(anova_d$`Pr(>F)`,"BH"), 6)
  dist_var <-
    vegan::permutest(disp,
                     permutations = 9999,
                     pairwise = T)
  return(list(dist_var, p_adj, disp))
}


betadisper(vegan::vegdist(fungi_ev_100itr, method="bray"), meta_fungi_100itr[,"Status"], type = c("centroid"))

BetadispExtr(fungi_ev_100itr, meta_fungi_100itr,"Status")

# Extract all betadisper models 
BetadispAll <- function(df, meta){
  
  betadisp_all <-
    data.frame(
      rbind(
        BetadispExtr(df, meta, "Status")[[1]]$tab,
        BetadispExtr(df, meta, "Niche")[[1]]$tab,
        BetadispExtr(df, meta, "Treatment")[[1]]$tab
      ),
      Padj = c(
        BetadispExtr(df, meta, "Status")[[2]],
        BetadispExtr(df, meta, "Niche")[[2]],
        BetadispExtr(df, meta, "Treatment")[[2]]
      ),
      Model = c(rep("Status", 2), rep("Niche", 2), rep("Treatment", 2))
    )
  return(betadisp_all)
}


# This does not check for interaction factors - we could do that but we need
# to combine the variables, for example paste(Status, Niche, sep="-")
betadisper_fungi <-
  BetadispAll(fungi_ev_100itr, meta_fungi_100itr)
betadisper_fungi

betadisper_bacteria <-
  BetadispAll(bact_ev_100itr, meta_bact_100itr)
betadisper_bacteria


# SPLIT PERMANOVAs ----------------------------------------------------------
SplitAdjadonis2 <- function(dataframe, metadata, stratum, Variable, Var1, Var2){
  require(tidyverse)
  require(vegan)
  
  otu <-
    dataframe %>%
    rownames_to_column("Sample_ID") %>% 
    left_join(metadata, by = "Sample_ID") %>% 
    filter(Niche %in% c(Variable)) %>% 
    column_to_rownames("Sample_ID") %>% 
    dplyr::select(where( ~is.numeric(.x) && sum(.x) != 0)) 
  
  print(otu %>% dim())
  
  # One example on how to handle formulas within a function
  # useful for model with interactions so you can test different
  # order of factors
  outcome <- "otu"
  variables <- c(Var1, Var2)
  
  model <-
    as.formula(
    paste(outcome,
          paste(variables, collapse = "*"),
          sep = " ~ "))
  map <-
  metadata %>% 
    filter(Niche %in% c(Variable)) 
  
  print(
    table(map$Treatment, map$Status)
    )
  
  df_adonis <-
    adonis2(model, 
            map,
            method = "bray",
            starta = stratum,
            permutations = 999)
  
  df_adj <-
    cbind(df_adonis,
          as.data.frame(p.adjust(df_adonis$`Pr(>F)`, 
                                 method = "BH")) %>% 
            rename("p.adj"=1))
  
  ad_list <-
    list(formula, df_adonis,df_adj)
  
  return(ad_list)
}


set.seed(21)
adonis_split_test <-
SplitAdjadonis2(fungi_ev_100itr, 
                meta_fungi_100itr,
                stratum = NULL,
                "Soil",
                "Treatment", 
                "Status")

adonis_split_test

# Combine split adonis2 
SplitAdjadonis2All <- function(dataframe, metadata, Variable){
    adonis_all <-
      data.frame(
        rbind(SplitAdjadonis2(dataframe, metadata, stratum = NULL, Variable,"Treatment", "Status")[[3]],
              SplitAdjadonis2(dataframe, metadata, stratum = "Treatment", Variable,"Treatment", "Status")[[3]],
              SplitAdjadonis2(dataframe, metadata, stratum = NULL, Variable,"Status", "Treatment")[[3]],
              SplitAdjadonis2(dataframe, metadata, stratum = "Status", Variable,"Status", "Treatment")[[3]]),
        Factor = rep(Variable, times=20),
        Strata = c(rep("-", times=5), rep("Treatment", times=5), rep("-", times=5), rep("Status", times=5)))
    
    return(adonis_all)
  }

# Fungi 
set.seed(100)
SplitAdjadonis2All(fungi_ev_100itr, meta_fungi_100itr, "Soil")
set.seed(101)
SplitAdjadonis2All(fungi_ev_100itr, meta_fungi_100itr, "Root")
set.seed(102)
SplitAdjadonis2All(fungi_ev_100itr, meta_fungi_100itr, "Rhizosphere")



# For Inflorescence and Senescent Leaf we have to run it separately
set.seed(105)
SplitAdjadonis2(fungi_ev_100itr, meta_fungi_100itr, stratum = NULL, "Inflorescence", Var1 = "Treatment", Var2 = NULL)[[3]]
set.seed(106)
SplitAdjadonis2(fungi_ev_100itr, meta_fungi_100itr, stratum = NULL, "Leaf", Var1 = "Treatment", Var2 = NULL)[[3]]


# Bacteria
set.seed(107)
SplitAdjadonis2All(bact_ev_100itr, meta_bact_100itr, "Soil")
set.seed(108)
SplitAdjadonis2All(bact_ev_100itr, meta_bact_100itr, "Root")
set.seed(109)
SplitAdjadonis2All(bact_ev_100itr, meta_bact_100itr, "Rhizosphere")


# For Inflorescence and Senescent Leaf we have to run it separately
set.seed(111)
SplitAdjadonis2(bact_ev_100itr, meta_bact_100itr, stratum = NULL, "Inflorescence", Var1 = "Treatment", Var2 = NULL)[[3]]
set.seed(112)
SplitAdjadonis2(bact_ev_100itr, meta_bact_100itr, stratum = NULL, "Leaf", Var1 = "Treatment", Var2 = NULL)[[3]]


# SPLIT BETA-SISPERSION --------------------------------------------------------
SplitBetadispExtr <- function(dataframe, metadata, Variable, Var1){
  require(tidyverse)
  require(vegan)
  otu <-
    dataframe %>%
    rownames_to_column("Sample_ID") %>% 
    left_join(metadata, by = "Sample_ID") %>% 
    filter(Niche %in% c(Variable)) %>% 
    column_to_rownames("Sample_ID") %>% 
    dplyr::select(where( ~is.numeric(.x) && sum(.x) != 0)) 
  
  map <-
    metadata %>% 
    filter(Niche %in% c(Variable)) 
  
  disp <-
    betadisper(
      vegan::vegdist(otu, method="bray"),
      map[,Var1])
  anova_d <-
    anova(disp,
          permutations = how(nperm=9999))
  p_adj <-
    round(p.adjust(anova_d$`Pr(>F)`,"BH"), 4)
  dist_var <-
    vegan::permutest(disp,
                     permutations = 9999,
                     pairwise = T)
  return(list(dist_var, p_adj, disp))
}



SplitBetadispExtr(fungi_ev_100itr, meta_fungi_100itr,"Soil", "Status")[[1]]$tab
SplitBetadispExtr(fungi_ev_100itr, meta_fungi_100itr,"Soil", "Treatment")
SplitBetadispExtr(fungi_ev_100itr, meta_fungi_100itr,"Soil", "Treatment")

# Extract all betadisper models for split dataframe by Niche 
SplitBetadispAll <- function(df, meta){
  
  betadisp_all <-
    data.frame(
      rbind(
        SplitBetadispExtr(df, meta,"Soil", "Status")[[1]]$tab,
        SplitBetadispExtr(df, meta,"Soil", "Treatment")[[1]]$tab,
        SplitBetadispExtr(df, meta,"Root", "Status")[[1]]$tab,
        SplitBetadispExtr(df, meta,"Root", "Treatment")[[1]]$tab,
        SplitBetadispExtr(df, meta,"Rhizosphere", "Status")[[1]]$tab,
        SplitBetadispExtr(df, meta,"Rhizosphere", "Treatment")[[1]]$tab,
        SplitBetadispExtr(df, meta,"Leaf", "Treatment")[[1]]$tab,
        SplitBetadispExtr(df, meta,"Inflorescence", "Treatment")[[1]]$tab
      ),
      Padj = c(
        SplitBetadispExtr(df, meta,"Soil", "Status")[[2]],
        SplitBetadispExtr(df, meta,"Soil", "Treatment")[[2]],
        SplitBetadispExtr(df, meta,"Root", "Status")[[2]],
        SplitBetadispExtr(df, meta,"Root", "Treatment")[[2]],
        SplitBetadispExtr(df, meta,"Rhizosphere", "Status")[[2]],
        SplitBetadispExtr(df, meta,"Rhizosphere", "Treatment")[[2]],
        SplitBetadispExtr(df, meta,"Leaf", "Treatment")[[2]],
        SplitBetadispExtr(df, meta,"Inflorescence", "Treatment")[[2]]
      ),
      Factor = c(rep("Soil", 4), rep("Root", 4), rep("Rhizosphere", 4), rep("Leaf", 2),
                rep("Inflorescence", 2)),
      Model = c(rep("Status", 2), rep("Treatment", 2),rep("Status", 2), rep("Treatment", 2),
                 rep("Status", 2), rep("Treatment", 2),rep("Treatment", 2),
                rep("Treatment", 2))
    )
  return(betadisp_all)
}


# This does not check for interaction factors - we could do that but we need
# to combine the variables, for example paste(Status, Niche, sep="-")
split_betadisper_fungi <-
  SplitBetadispAll(fungi_ev_100itr, meta_fungi_100itr)

split_betadisper_fungi

split_betadisper_bact <-
  SplitBetadispAll(bact_ev_100itr, meta_bact_100itr)

split_betadisper_bact


# ***********************************************-------------------------------
# ALPHA DIVERSITY --------------------------------------------------------------
# By default vegan rarefy() is using the bootstrapped randomized 
# approach as default. You should filter the samples that are lower
# the rarefaction cutoff you chose otherwise the function will 
# return non-rarefied species richness (and standard error = 0)
# Empirical rarefaction approach: run 100 rarefaction of the data 
#and calculate richness and Shannon

# RICHNESS
GetRich <-function(dataframe, depth_level){
  observed_list <- list(mode = "list")

  set.seed(032)
  
  for (i in 1:99) {
    output  <- vegan::rrarefy(dataframe, sample = depth_level) %>%
      specnumber(x = ., MARGIN = 1)
    observed_list[[length(observed_list) + 1]] <- output
  }
  # Convert the named vectors list into a dataframe. I did not find
  # a good tidtverse solution for this
  alpha_out <-
    as.data.frame(do.call(cbind, observed_list)) %>%
    dplyr::select(-mode) %>% 
    as_tibble(rownames="Sample_ID") %>% 
    pivot_longer(-Sample_ID, names_to = "Round", values_to = "Observed") %>% 
    group_by(Sample_ID) %>% 
    mutate_at(3, as.numeric) %>% 
    summarise(Mean = mean(Observed),
              SD = sd(Observed))
  
  return(alpha_out)
}


rrarefy(df_filt_fungi, sample = min_depth_fungi)

rich_fungi <-
  GetRich(df_filt_fungi, min_depth_fungi)
rich_fungi

rich_bact <-
  GetRich(df_filt_bact, min_depth_bact)
rich_bact

# DIVERSITY INDEX
GetIndex <-function(dataframe, depth_level, measure){
index_list <- list(mode = "list")

set.seed(042)

for (i in 1:99) {
  output  <- vegan::rrarefy(dataframe, sample = depth_level) %>%
    diversity(x=., index = measure, MARGIN=1)
  index_list[[length(index_list) + 1]] <- output
}

# I used a trick to avoid explicitly naming the variables since they depend
# on "measure". just listing the functions in summarize_at worked. Also, remember
# when you group by a column that column is not a columns anymore but a grouping 
# factor, so won't count as column in the numbering.
index_out <-
  as.data.frame(do.call(cbind, index_list)) %>%
  dplyr::select(-mode) %>% 
  as_tibble(rownames="Sample_ID") %>% 
  pivot_longer(-Sample_ID, names_to = "Round", values_to = measure) %>% 
  mutate_at(3, as.numeric) %>% 
  group_by(Sample_ID) %>% 
  dplyr::summarise_at(2, list(Mean = mean, SD = sd),na.rm = TRUE)

  return(index_out)
}


shannon_fungi <-
  GetIndex(df_filt_fungi, min_depth_fungi, "shannon")
shannon_fungi 

shannon_bact <-
  GetIndex(df_filt_bact, min_depth_bact, "shannon")
shannon_bact 

alpha_fungi <-
  full_join(rich_fungi, shannon_fungi, by = "Sample_ID") %>%
  rename(
    Richness = Mean.x,
    Richness_SD = SD.x,
    Shannon = Mean.y,
    Shannon_SD = SD.y) %>%
  mutate(Shannon_Eq = Shannon / log(Richness)) %>%
  left_join(y = meta_fungi_100itr, by = "Sample_ID") %>% 
  mutate(Group = paste(Niche, Treatment, sep="_"))

alpha_fungi

alpha_bact <-
  full_join(rich_bact, shannon_bact, by = "Sample_ID") %>%
  rename(
    Richness = Mean.x,
    Richness_SD = SD.x,
    Shannon = Mean.y,
    Shannon_SD = SD.y) %>%
  mutate(Shannon_Eq = Shannon / log(Richness)) %>%
  left_join(y = meta_bact_fix, by = "Sample_ID" )%>% 
  mutate(Group = paste(Niche, Treatment, sep="_"))

alpha_bact

# Variations in alpha metrics ------------------------
# Plotting just simple histogrham to have a sense of the variation
# across samples

alphaHisto <- function(){
  require(grid)
  require(gridExtra)
  
title1 = text_grob("Rarefied richness", size = 12, face = 2)
title2 = text_grob("Shannon index", size = 12, face = 2)
title3 = text_grob("Standard deviation", size = 12, face = 2)

alpha_histo_plot <-
  ggarrange(
    grid.arrange(rbind(
      alpha_fungi %>%
        dplyr::select(Sample_ID, Richness, Richness_SD, Shannon, Shannon_SD,Shannon_Eq) %>% 
        pivot_longer(-Sample_ID, names_to = "Group", values_to = "Value") %>% 
        mutate(Kingdom = "Fungi"),
      alpha_bact %>%
        dplyr::select(Sample_ID, Richness, Richness_SD, Shannon, Shannon_SD,Shannon_Eq) %>% 
        pivot_longer(-Sample_ID, names_to = "Group", values_to = "Value") %>% 
        mutate(Kingdom = "Bacteria")) %>% 
        filter(Group %in% "Richness") %>%
        ggplot(aes(x = Value)) +
        geom_histogram(bins = 100) +
        facet_wrap(~Kingdom) +
        themeStats(), top=title1),
    grid.arrange(rbind(
      alpha_fungi %>%
        dplyr::select(Sample_ID, Richness, Richness_SD, Shannon, Shannon_SD,Shannon_Eq) %>% 
        pivot_longer(-Sample_ID, names_to = "Group", values_to = "Value") %>% 
        mutate(Kingdom = "Fungi"),
      alpha_bact %>%
        dplyr::select(Sample_ID, Richness, Richness_SD, Shannon, Shannon_SD,Shannon_Eq) %>% 
        pivot_longer(-Sample_ID, names_to = "Group", values_to = "Value") %>% 
        mutate(Kingdom = "Bacteria")) %>% 
        filter(Group %in% "Richness_SD") %>%
        ggplot(aes(x = Value)) +
        geom_histogram(bins = 100) +
        facet_wrap(~Kingdom) +
        themeStats(), top=title3),
    grid.arrange(rbind(
      alpha_fungi %>%
        dplyr::select(Sample_ID, Richness, Richness_SD, Shannon, Shannon_SD,Shannon_Eq) %>% 
        pivot_longer(-Sample_ID, names_to = "Group", values_to = "Value") %>% 
        mutate(Kingdom = "Fungi"),
      alpha_bact %>%
        dplyr::select(Sample_ID, Richness, Richness_SD, Shannon, Shannon_SD,Shannon_Eq) %>% 
        pivot_longer(-Sample_ID, names_to = "Group", values_to = "Value") %>% 
        mutate(Kingdom = "Bacteria")) %>% 
        filter(Group %in% "Shannon") %>%
        ggplot(aes(x = Value)) +
        geom_histogram(bins = 100) +
        facet_wrap(~Kingdom) +
        themeStats(), top=title2),
    grid.arrange(rbind(
      alpha_fungi %>%
        dplyr::select(Sample_ID, Richness, Richness_SD, Shannon, Shannon_SD,Shannon_Eq) %>% 
        pivot_longer(-Sample_ID, names_to = "Group", values_to = "Value") %>% 
        mutate(Kingdom = "Fungi"),
      alpha_bact %>%
        dplyr::select(Sample_ID, Richness, Richness_SD, Shannon, Shannon_SD,Shannon_Eq) %>% 
        pivot_longer(-Sample_ID, names_to = "Group", values_to = "Value") %>% 
        mutate(Kingdom = "Bacteria")) %>% 
        filter(Group %in% "Shannon_SD") %>%
        ggplot(aes(x = Value)) +
        geom_histogram(bins = 100) +
        facet_wrap(~Kingdom) +
        themeStats(), top=title3),
  ncol = 1,
  nrow = 4, 
  align = "hv")

return(alpha_histo_plot)
}

alphaHisto()


# Extract significant differences --------------------------------------------
CompSampl <- function(df, formula){
  require(multcompView)
  require(tidyverse)
  df1 <-
    df %>%
    compare_means(data = ., formula, 
                  method = "wilcox.test", p.adjust.method = "BH") %>% 
    dplyr::select(group1, group2, p.adj)
  df2 <-
    df1 %>% 
    dplyr::select(group2, group1, p.adj)
  #print(df2)
  colnames(df2) <- 
    c("group1", "group2", "p.adj") 
  rbind(df1, df2) %>% 
    as.data.frame() %>% 
    xtabs(p.adj ~ group2 + group1, data=.) %>% 
    as.dist(diag = TRUE) -> dist
  #print(dist)
  res <-
    as.data.frame(multcompLetters(x=dist, reversed = FALSE)['Letters'])
  return(res)
}

alpha_fungi %>% 
  filter(Status %in% "Pre-Transplant") %>% 
  mutate(Group = as.factor(Group),
         Group = fct_relevel(Group, 
                             "Root_GMO","Root_non-GMO",
                             "Rhizosphere_GMO","Rhizosphere_non-GMO",
                             "Soil_GMO","Soil_non-GMO")) %>% 
  arrange(factor(Group)) %>% 
  CompSampl(df = ., formula(Richness ~ Group)) 

# Reformat letters DF ---------------------------------------------
MakeLetters <- function(dataframe, formula, Var){

  df_letters <-
  dataframe %>% 
    filter(Status %in% c(Var)) %>% 
    CompSampl(df = ., formula) %>% 
    rownames_to_column("Group") 
  
  if (Var == "Post-Transplant"){
    fix_letters <-
    df_letters %>% 
      mutate(Group = as.factor(Group),
             Group = fct_relevel(Group, 
                                 "Inflorescence_GMO","Inflorescence_non-GMO",
                                 "Leaf_GMO","Leaf_non-GMO",
                                 "Root_GMO","Root_non-GMO",
                                 "Rhizosphere_GMO","Rhizosphere_non-GMO",
                                 "Soil_GMO","Soil_non-GMO"),
             Letters = as.factor(Letters)) %>% 
      arrange(factor(Group, levels=c("Inflorescence_GMO","Inflorescence_non-GMO",
                                     "Leaf_GMO","Leaf_non-GMO",
                                     "Root_GMO","Root_non-GMO",
                                     "Rhizosphere_GMO","Rhizosphere_non-GMO",
                                     "Soil_GMO","Soil_non-GMO")))
  } else {
    fix_letters <-
    df_letters %>% 
      mutate(Group = as.factor(Group),
             Group = fct_relevel(Group, 
                                 "Root_GMO","Root_non-GMO",
                                 "Rhizosphere_GMO","Rhizosphere_non-GMO",
                                 "Soil_GMO","Soil_non-GMO"),
             Letters = as.factor(Letters)) %>% 
      arrange(factor(Group, levels=c("Root_GMO","Root_non-GMO",
                                     "Rhizosphere_GMO","Rhizosphere_non-GMO",
                                     "Soil_GMO","Soil_non-GMO")))
  }
  return(fix_letters)
}


rich_sig_post_fungi <-
  MakeLetters(alpha_fungi, formula(Richness ~ Group), "Post-Transplant")
rich_sig_pre_fungi <-
  MakeLetters(alpha_fungi, formula(Richness ~ Group), "Pre-Transplant") 
shan_sig_post_fungi <-
  MakeLetters(alpha_fungi, formula(Shannon ~ Group), "Post-Transplant")
shan_sig_pre_fungi <-
  MakeLetters(alpha_fungi, formula(Shannon ~ Group), "Pre-Transplant") 

rich_sig_post_bact <-
  MakeLetters(alpha_bact, formula(Richness ~ Group), "Post-Transplant")
rich_sig_pre_bact <-
  MakeLetters(alpha_bact, formula(Richness ~ Group), "Pre-Transplant") 
shan_sig_post_bact <-
  MakeLetters(alpha_bact, formula(Shannon ~ Group), "Post-Transplant")
shan_sig_pre_bact <-
  MakeLetters(alpha_bact, formula(Shannon ~ Group), "Pre-Transplant") 


# Plotting the data -------------------------------
PlotAlpha <- function(dataframe, Var, my_labels, y_limit, y_labels){
  rich_plot <-
    dataframe %>% 
    ggplot(aes(x = Treatment, y = get(Var), fill = Niche, colour = Niche)) +
    geom_boxplot(outlier.shape = 1, outlier.size = 0.8, outlier.stroke = 1,
                 position = position_dodge(preserve = "single"), alpha=0.6, lwd = 0.5) +
    stat_summary(fun=mean, geom="point", shape=18, size=1.9, color="black", fill="black") +
    stat_summary(geom = "text", label = my_labels$Letters, fun= max, 
                 aes(y = y_labels), size=3, color="black") +
    facet_grid(~ Niche, scales = "free_x", space="free_x") +
    theme_bw() +
    expand_limits(y = c(0, y_limit)) +
    theme(strip.text.x = element_text(angle = 0, size = 9,hjust = 0.5, vjust = 0.5, color = "black"),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5,color = "black"),
          legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
          legend.title =element_text(size = 9, face = "bold", hjust = 0.5,color = "black"),
          legend.text = element_text(angle = 0, size = 8, hjust = 0, vjust = 0.5, color = "black"),
          axis.text.x = element_text(angle = 33, size = 8, hjust = 1, vjust = 1.05, color = "black"),
          axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5, color = "black"),
          axis.title = element_text(angle = 0, size = 8, face = "bold"),
          strip.background = element_blank(), 
          legend.position = "bottom")
  guides(fill = guide_legend(order = 1, ncol = 3, title = "Niche",
                            override.aes = list(shape = 22, size=3)))
  return(rich_plot)
}


alpha_fungi %>% 
  filter(Status %in% c("Post-Transplant")) %>% 
  group_by(Group) %>% 
  summarise(mean = mean(Richness), sd = sd(Richness)) %>% 
  separate(Group, c("Niche", "Treatment"), sep = "_") %>% 
  mutate(
    Niche = as.factor(Niche),
    Niche = fct_relevel(Niche, 
    "Inflorescence","Leaf",
    "Root","Rhizosphere","Soil")) %>% 
  PlotAlpha(rich_sig_post_fungi,  y_limit = 430, y_labels=430) 

alpha_fungi %>% 
  filter(Status %in% c("Post-Transplant")) %>% 
  PlotAlpha("Richness", rich_sig_post_fungi,  y_limit = 430, y_labels=430) 

alpha_fungi %>% 
  filter(Status %in% c("Post-Transplant")) %>% 
  ggplot(aes(x = Treatment, y = Richness, fill = Niche)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1, outlier.stroke = 1,
               position = position_dodge(preserve = "single"), alpha=0.6, lwd = 0.5) +
  stat_summary(fun=mean, geom="point", shape=18, size=1.9, color="black", fill="black") +
  facet_grid(~ Niche, scales = "free_x", space="free_x") 


# Observed richness plot ----------------------------------------------------
Fig_richness <-
ggarrange(
  alpha_fungi %>%
    filter(Status %in% c("Post-Transplant")) %>%
    mutate(
      Niche = as.factor(Niche),
      Niche = fct_relevel(Niche,
                          "Inflorescence", "Leaf",
                          "Root", "Rhizosphere", "Soil")
    ) %>%
    PlotAlpha("Richness", rich_sig_post_fungi,  y_limit = 420, y_labels = 420) +
    scale_fill_manual(values = paletteCB6) +
    scale_color_manual(values = paletteCB6) +
    labs(title = "Fungi Post-Transplant", y = "Richness"),
  
  alpha_fungi %>%
    filter(Status %in% c("Pre-Transplant")) %>%
    mutate(Niche = as.factor(Niche),
           Niche = fct_relevel(Niche,
                               "Root", "Rhizosphere", "Soil")) %>%
    PlotAlpha("Richness",rich_sig_pre_fungi,  y_limit = 420, y_labels = 420) +
    scale_fill_manual(values = paletteCB6[c(3, 4, 5)]) +
    scale_color_manual(values = paletteCB6[c(3, 4, 5)]) +
    labs(title = "Fungi Pre-Transplant", y = "Richness"), 

  alpha_bact %>%
    filter(Status %in% c("Post-Transplant")) %>%
    mutate(
      Niche = as.factor(Niche),
      Niche = fct_relevel(Niche,
                          "Inflorescence", "Leaf",
                          "Root", "Rhizosphere", "Soil")
    ) %>%
    PlotAlpha("Richness",rich_sig_post_bact,  y_limit = 2100, y_labels = 2100) +
    scale_fill_manual(values = paletteCB6) +
    scale_color_manual(values = paletteCB6) +
    labs(title = "Bacteria Post-Transplant", y = "Richness"), 

  alpha_bact %>%
    filter(Status %in% c("Pre-Transplant")) %>%
    mutate(Niche = as.factor(Niche),
           Niche = fct_relevel(Niche,
                               "Root", "Rhizosphere", "Soil")) %>%
    PlotAlpha("Richness",rich_sig_pre_bact,  y_limit = 2100, y_labels = 2100) +
    scale_fill_manual(values = paletteCB6[c(3, 4, 5)]) +
    scale_color_manual(values = paletteCB6[c(3, 4, 5)]) +
    labs(title = "Bacteria Pre-Transplant", y = "Richness"),
  
widths = c(1, 0.67),
ncol = 2, nrow = 2,
labels = c("A", "B", "C", "D"),
common.legend = TRUE, 
legend = "bottom")


Fig_richness

title4=text_grob("Rarefied richness",size=12, face=2)
grid.arrange(Fig_richness, top = title4)

Fig_richness + 
  plot_annotation(
    title = "Fungal rarefied richness", 
    theme = theme(plot.title = element_text(
      size = 12, face="bold", hjust = 0.5))) +
  theme(legend.position = "none")


# Shannon index plot ------------------------------------------------------

EH_fungi <-
  ggarrange(
    alpha_fungi %>%
      filter(Status %in% c("Post-Transplant")) %>%
      mutate(
        Niche = as.factor(Niche),
        Niche = fct_relevel(Niche,
                            "Inflorescence", "Leaf",
                            "Root", "Rhizosphere", "Soil")
      ) %>%
      PlotAlpha("Shannon", shan_sig_post_fungi,y_limit = 7.5, y_labels = 7.5) +
      scale_fill_manual(values = paletteCB6) +
      scale_color_manual(values = paletteCB6) +
      labs(title = "Fungi Post-Transplant", y="Shannon"),
    
    alpha_fungi %>%
      filter(Status %in% c("Pre-Transplant")) %>%
      mutate(Niche = as.factor(Niche),
             Niche = fct_relevel(Niche,
                                 "Root", "Rhizosphere", "Soil")) %>%
      PlotAlpha("Shannon", shan_sig_pre_fungi,  y_limit = 7.5, y_labels = 7.5) +
      scale_fill_manual(values = paletteCB6[c(3, 4, 5)]) +
      scale_color_manual(values = paletteCB6[c(3, 4, 5)]) +
      labs(title = "Fungi Post-Transplant", y="Shannon"),
    
    alpha_bact %>%
      filter(Status %in% c("Post-Transplant")) %>%
      mutate(
        Niche = as.factor(Niche),
        Niche = fct_relevel(Niche,
                            "Inflorescence", "Leaf",
                            "Root", "Rhizosphere", "Soil")
      ) %>%
      PlotAlpha("Shannon", shan_sig_post_bact, y_limit = 7.5, y_labels = 7.5) +
      scale_fill_manual(values = paletteCB6) +
      scale_color_manual(values = paletteCB6) +
      labs(title = "Fungi Post-Transplant", y="Shannon"),
    
    alpha_bact %>%
      filter(Status %in% c("Pre-Transplant")) %>%
      mutate(Niche = as.factor(Niche),
             Niche = fct_relevel(Niche,
                                 "Root", "Rhizosphere", "Soil")) %>%
      PlotAlpha("Shannon", shan_sig_pre_bact,  y_limit = 7.5, y_labels = 7.5) +
      scale_fill_manual(values = paletteCB6[c(3, 4, 5)]) +
      scale_color_manual(values = paletteCB6[c(3, 4, 5)]) +
      labs(title = "Fungi Post-Transplant", y="Shannon"),
    
    widths = c(1, 0.67),
    ncol = 2, nrow = 2,
    labels = c("A", "B", "C", "D"),
    common.legend = TRUE, 
    legend = "bottom")


EH_fungi

title4=text_grob("Shannon Index",size=12, face=2)
grid.arrange(EH_fungi, top = title4)

# ***********************************************-------------------------------
# DIFFERENTIAL ABUNDANCE -------------------------------------------------------

# Split DF as needed
SplitDF <- function(dataframe, metadata, Var1, Var2){
  require(tidyverse)
  require(vegan)
  
  otu <-
    dataframe %>%
    rownames_to_column("Sample_ID") %>% 
    left_join(metadata, by = "Sample_ID") %>% 
    filter(Niche %in% c(Var1) & Status %in% c(Var2)) %>% 
    column_to_rownames("Sample_ID") %>% 
    dplyr::select(where( ~is.numeric(.x) && sum(.x) != 0)) 

  map <-
    metadata %>% 
    filter(Niche %in% c(Var1) & Status %in% c(Var2))
  
  dim(otu) %>%  print()
  dim(map) %>%  print()
  
  return(list(otu, map))
}

SplitDF(fungi_ev_100itr, meta_fungi_100itr, "Leaf", "Post-Transplant")
SplitDF(fungi_ev_100itr, meta_fungi_100itr, "Leaf", "Post-Transplant")[[1]]


ExtractLongDF <- function(dataframe, metadata, taxonomy, Var1, Var2){
  require(tidyverse)
  split_df <-
    SplitDF(dataframe, metadata, Var1, Var2)
  
  long_df <-  
    split_df[[1]] %>% 
    t() %>%
    as.data.frame() %>%
    mutate(Sum = rowSums(x = .)) %>%  
    left_join(x = rownames_to_column(., var = "ASV_ID"),
              y = taxonomy,
              by = "ASV_ID") %>%
    dplyr::select(.data = ., starts_with("TR"), Taxonomy) %>%
    group_by(Taxonomy) %>% 
    pivot_longer(cols=c(-Taxonomy),names_to="Sample_ID") %>% 
    mutate(Abund = value / sum(value)) %>% 
    left_join(x=., 
              y= as.data.frame(split_df[[2]]),
              by="Sample_ID") %>% 
    na.omit()
  
  return(long_df)
}


split_fungi_pre_root <-
  ExtractLongDF(fungi_ev_100itr, meta_fungi_100itr,fungi_taxonomy_fix, "Root", "Pre-Transplant")

split_fungi_pre_rhizo <-
  ExtractLongDF(fungi_ev_100itr, meta_fungi_100itr,fungi_taxonomy_fix,"Rhizosphere", "Pre-Transplant")
split_fungi_pre_soil <-
  ExtractLongDF(fungi_ev_100itr, meta_fungi_100itr, fungi_taxonomy_fix, "Soil", "Pre-Transplant")

split_fungi_post_leaf <-
  ExtractLongDF(fungi_ev_100itr, meta_fungi_100itr,fungi_taxonomy_fix,"Leaf", "Post-Transplant")
split_fungi_post_inf <-
  ExtractLongDF(fungi_ev_100itr, meta_fungi_100itr,fungi_taxonomy_fix,"Inflorescence", "Post-Transplant")
split_fungi_post_root <-
ExtractLongDF(fungi_ev_100itr, meta_fungi_100itr,fungi_taxonomy_fix,"Root", "Post-Transplant")
split_fungi_post_rhizo <-
  ExtractLongDF(fungi_ev_100itr, meta_fungi_100itr,fungi_taxonomy_fix,"Rhizosphere", "Post-Transplant")
split_fungi_post_soil <-
  ExtractLongDF(fungi_ev_100itr, meta_fungi_100itr,fungi_taxonomy_fix,"Soil", "Post-Transplant")

split_bact_pre_root <-
  ExtractLongDF(bact_ev_100itr, meta_bact_100itr,bact_taxonomy_fix, "Root", "Pre-Transplant")
split_bact_pre_rhizo <-
  ExtractLongDF(bact_ev_100itr, meta_bact_100itr,bact_taxonomy_fix,"Rhizosphere", "Pre-Transplant")
split_bact_pre_soil <-
  ExtractLongDF(bact_ev_100itr, meta_bact_100itr, bact_taxonomy_fix, "Soil", "Pre-Transplant")

split_bact_post_leaf <-
  ExtractLongDF(bact_ev_100itr, meta_bact_100itr,bact_taxonomy_fix,"Leaf", "Post-Transplant")
split_bact_post_inf <-
  ExtractLongDF(bact_ev_100itr, meta_bact_100itr,bact_taxonomy_fix,"Inflorescence", "Post-Transplant")
split_bact_post_root <-
  ExtractLongDF(bact_ev_100itr, meta_bact_100itr,bact_taxonomy_fix,"Root", "Post-Transplant")
split_bact_post_rhizo <-
  ExtractLongDF(bact_ev_100itr, meta_bact_100itr,bact_taxonomy_fix,"Rhizosphere", "Post-Transplant")
split_bact_post_soil <-
  ExtractLongDF(bact_ev_100itr, meta_bact_100itr,bact_taxonomy_fix,"Soil", "Post-Transplant")


GetWilcoxPval <- function(df, formula){
  require(tidyverse)
  require(broom)
  res_pval <-
    df %>% 
    nest(data = -Taxonomy) %>%
    mutate(test = map(.x=data, ~wilcox.test(formula, data=.x) %>% tidy)) %>%
    unnest(test) %>%
    mutate(p.adjust = p.adjust(p.value, method="BH")) %>%
    filter(p.adjust < 0.05) %>%
    dplyr::select(Taxonomy, p.adjust)
  
  return(res_pval) 
}

# Plotting ---------------------------------------------------------------------
# I think I prefer median pointrange plots

PlotBoxplot <- function(df){
  require(tidyverse)
  bar_plot <- 
    df %>%
    ggplot(aes(., x = sqrt(value), y = Taxonomy, fill = Treatment, colour=Treatment)) +
    geom_boxplot(outlier.shape = 1, outlier.size = 0.8, outlier.stroke = 1,
                 position = position_dodge(preserve = "single"), alpha=0.6, lwd = 0.5) +
    #stat_summary(fun=mean, geom="point", shape=18, size=1.9, color="black", fill="black") +
    theme_classic() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
          axis.title = element_text(size = 10, hjust = 0.5, face = "bold"),
          strip.background = element_blank()) +
  labs(x=expression(sqrt("Rarefied abundance"))) +
         scale_y_discrete(limits = rev)
  return(bar_plot)
}


rbind(
GetWilcoxPval(split_fungi_pre_soil, formula(value ~ Treatment)) %>% 
  mutate(Niche = "Soil"),
GetWilcoxPval(split_fungi_pre_rhizo, formula(value ~ Treatment)) %>% 
  mutate(Niche = "Rhizosphere"),
GetWilcoxPval(split_fungi_pre_root, formula(value ~ Treatment)) %>% 
  mutate(Niche = "Root"))

# Plotting all them together is a mess, plotting by Niche and then 
# we'll decide how to imporve it.

rbind(
GetWilcoxPval(split_fungi_pre_soil, formula(value ~ Treatment)) %>% 
  left_join(x=., y=split_fungi_pre_soil, 
            multiple = "all", by = "Taxonomy"),
GetWilcoxPval(split_fungi_pre_rhizo, formula(value ~ Treatment)) %>% 
  left_join(x=., y=split_fungi_pre_rhizo, 
            multiple = "all", by = "Taxonomy"),
GetWilcoxPval(split_fungi_pre_root, formula(value ~ Treatment)) %>% 
  left_join(x=., y=split_fungi_pre_root, 
            multiple = "all", by = "Taxonomy")) %>% 
  PlotBoxplot(df = .) +
  facet_grid(~Niche) #scales = "free", space="free"
  

# Plotting all Niches and Status differential abundant ASVs -------------------

# Fungi Pre-Transplant -----------------------------------------------------
GetWilcoxPval(split_fungi_pre_soil, formula(value ~ Treatment)) %>%
    left_join(
      x = .,
      y = split_fungi_pre_soil,
      multiple = "all",
      by = "Taxonomy") %>%
    PlotBoxplot(df = .) + 
  labs(title = "Soil Pre-Transplant")

GetWilcoxPval(split_fungi_pre_rhizo, formula(value ~ Treatment)) %>%
    left_join(
      x = .,
      y = split_fungi_pre_rhizo,
      multiple = "all",
      by = "Taxonomy") %>%
    PlotBoxplot(df = .) +
  labs(title = "Rhizosphere Pre-Transplant")

GetWilcoxPval(split_fungi_pre_root, formula(value ~ Treatment)) %>%
    left_join(
      x = .,
      y = split_fungi_pre_root,
      multiple = "all",
      by = "Taxonomy") %>%
  PlotBoxplot(df = .) +
  labs(title = "Root Pre-Transplant")


# Fungi Post-Transplant ----------------------------------------------
GetWilcoxPval(split_fungi_post_leaf, formula(value ~ Treatment)) %>%
  left_join(
    x = .,
    y = split_fungi_post_leaf,
    multiple = "all",
    by = "Taxonomy") %>%
  PlotBoxplot(df = .) + 
  labs(title = "Leaf Post-Transplant")

GetWilcoxPval(split_fungi_post_inf, formula(value ~ Treatment)) %>%
  left_join(
    x = .,
    y = split_fungi_post_inf,
    multiple = "all",
    by = "Taxonomy") %>%
  PlotBoxplot(df = .) + 
  labs(title = "Inflorescence Post-Transplant")

GetWilcoxPval(split_fungi_post_root, formula(value ~ Treatment)) %>%
  left_join(
    x = .,
    y = split_fungi_post_root,
    multiple = "all",
    by = "Taxonomy") %>%
  PlotBoxplot(df = .) + 
  labs(title = "Root Post-Transplant")

GetWilcoxPval(split_fungi_post_rhizo, formula(value ~ Treatment)) %>%
  left_join(
    x = .,
    y = split_fungi_post_rhizo,
    multiple = "all",
    by = "Taxonomy") %>%
  PlotBoxplot(df = .) + 
  labs(title = "Rhizosphere Post-Transplant")

GetWilcoxPval(split_fungi_post_soil, formula(value ~ Treatment)) %>%
  left_join(
    x = .,
    y = split_fungi_post_soil,
    multiple = "all",
    by = "Taxonomy") %>%
  PlotBoxplot(df = .) + 
  labs(title = "Soil Post-Transplant")
    
# Bacteria Pre-Transplant -------------------------------------------
GetWilcoxPval(split_bact_pre_soil, formula(value ~ Treatment)) %>%
  left_join(
    x = .,
    y = split_bact_pre_soil,
    multiple = "all",
    by = "Taxonomy") %>%
  PlotBoxplot(df = .) + 
  labs(title = "Soil Pre-Transplant")

GetWilcoxPval(split_bact_pre_rhizo, formula(value ~ Treatment)) %>%
  left_join(
    x = .,
    y = split_bact_pre_rhizo,
    multiple = "all",
    by = "Taxonomy") %>%
  PlotBoxplot(df = .) +
  labs(title = "Rhizosphere Pre-Transplant")

GetWilcoxPval(split_bact_pre_root, formula(value ~ Treatment)) %>%
  left_join(
    x = .,
    y = split_bact_pre_root,
    multiple = "all",
    by = "Taxonomy") %>%
  PlotBoxplot(df = .) +
  labs(title = "Root Pre-Transplant")


# Bacteria Post-Transplant -----------------------------------------------------
GetWilcoxPval(split_bact_post_leaf, formula(value ~ Treatment)) %>%
  left_join(
    x = .,
    y = split_bact_post_leaf,
    multiple = "all",
    by = "Taxonomy") %>%
  PlotBoxplot(df = .) + 
  labs(title = "Leaf Post-Transplant")

GetWilcoxPval(split_bact_post_inf, formula(value ~ Treatment)) %>%
  left_join(
    x = .,
    y = split_bact_post_inf,
    multiple = "all",
    by = "Taxonomy") %>%
  PlotBoxplot(df = .) + 
  labs(title = "Inflorescence Post-Transplant")

GetWilcoxPval(split_bact_post_root, formula(value ~ Treatment)) %>%
  left_join(
    x = .,
    y = split_bact_post_root,
    multiple = "all",
    by = "Taxonomy") %>%
  PlotBoxplot(df = .) + 
  labs(title = "Root Post-Transplant")

GetWilcoxPval(split_bact_post_rhizo, formula(value ~ Treatment)) %>%
  left_join(
    x = .,
    y = split_bact_post_rhizo,
    multiple = "all",
    by = "Taxonomy") %>%
  PlotBoxplot(df = .) + 
  labs(title = "Rhizosphere Post-Transplant")

GetWilcoxPval(split_bact_post_soil, formula(value ~ Treatment)) %>%
  left_join(
    x = .,
    y = split_bact_post_soil,
    multiple = "all",
    by = "Taxonomy") %>%
  PlotBoxplot(df = .) + 
  labs(title = "Soil Post-Transplant")

# ***********************************************-------------------------------
# DESEQ2 Differential Abundance ------------------------------------------------
# For other methods https://www.nicholas-ollberding.com/post/identifying-differentially-abundant-features-in-microbiome-data/

library(phyloseq)
library(DESeq2)
library(Biostrings)
library(apeglm)

colnames(physeq_fungi)
rownames(meta_fungi_100itr)

dim(fungi_ev_100itr)
dim(meta_fungi_100itr)

fungi_taxonomy_100itr <-
  fungi_taxonomy_fix %>%
  subset(ASV_ID %in% colnames(fungi_ev_100itr))

dim(fungi_taxonomy_100itr)

physeq_fungi_ev <- phyloseq::phyloseq(
  otu_table(fungi_ev_100itr %>%
              t(x = .) %>%
              as.data.frame(), taxa_are_rows = TRUE),
  sample_data(meta_fungi_100itr %>%
                column_to_rownames("Sample_ID")),
  tax_table(as.matrix(
    fungi_taxonomy_100itr %>%
      column_to_rownames("ASV_ID")
  ))
)

physeq_fungi_ev
head(physeq_fungi_ev@tax_table)
head(physeq_fungi_ev@otu_table)

bact_taxonomy_100itr <-
  bact_taxonomy_fix %>%
  subset(ASV_ID %in% colnames(bact_ev_100itr))


physeq_bact_ev <- phyloseq::phyloseq(
  otu_table(bact_ev_100itr %>%
              t(x = .) %>%
              as.data.frame(), taxa_are_rows = TRUE),
  sample_data(meta_bact_100itr %>%
                column_to_rownames("Sample_ID")),
  tax_table(as.matrix(
    bact_taxonomy_100itr %>%
      column_to_rownames("ASV_ID")
  ))
)

physeq_bact_ev
head(physeq_bact_ev@tax_table)
head(physeq_bact_ev@otu_table)


D2Stats <- function(physeq, des_df){
  require(DESeq2)
  require(tidyverse)
  
  # calculate geometric means prior to estimate size factors
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  
  #geoMeans <- apply(counts(des_df), 1, gm_mean)
  #diagdds <- estimateSizeFactors(des_df, geoMeans = geoMeans)
  
  #diagdds <- DESeq(des_df, fitType="parametric")
  diagdds <- DESeq(des_df, 
                   test = "Wald", 
                   fitType = "local", 
                   sfType = "poscounts")
  
  plotMA(diagdds)
  
  # Calculate results (p < 0.01)
  results <- results(diagdds, 
                     alpha = 0.01, 
                     cooksCutoff = FALSE,
                     contrast = c("Treatment", "GMO","non-GMO"))
  
  # Tabulate results
  sigtab_results <- results[which(results$padj < 0.05), ]
  summary(sigtab_results)
  
  sigtab_results_tax <- cbind(
    as(sigtab_results, "data.frame"), 
    as(tax_table(physeq)[row.names(sigtab_results), ], "matrix"))
  
  return(sigtab_results_tax)
}


# Function to generate all DESeq2 stats at once
D2Stats <- function(physeq, des_df){
  require(DESeq2)
  require(apeglm)
  require(tidyverse)
  
  # calculate geometric means prior to estimate size factors
  #gm_mean = function(x, na.rm=TRUE){
  #  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  #}
  
  #geoMeans = apply(counts(des_df), 1, gm_mean)
  #diagdds = estimateSizeFactors(des_df, geoMeans = geoMeans)
  #diagdds = DESeq(des_df, fitType="parametric")
  
  diagdds <- DESeq(des_df, test = "Wald", fitType = "local", sfType = "poscounts")
  
  plotDispEsts(diagdds)
  plotMA(diagdds)
  
  # using 'apeglm' for LFC shrinkage.
  #res_df <- DESeq2::lfcShrink(diagdds, coef=2, type="normal", quiet=TRUE) 
  
  # Calculate results 
  res_df <- DESeq2::results(diagdds, 
                     alpha = 0.05, 
                     contrast = c("Treatment", "GMO","non-GMO"))
  
  # Tabulate results
  #sigtab_results <-
  #data.frame(res_df) %>%
  #  dplyr::arrange(padj) %>% 
  #  dplyr::filter(padj < 0.05)
  
  sigtab_results <- res_df[which(res_df$padj < 0.05), ]
  print(summary(sigtab_results))
  
  sigtab_results_tax <- cbind(
    as(sigtab_results, "data.frame"), 
    as(tax_table(physeq)[row.names(sigtab_results), ], "matrix"))
  
  return(sigtab_results_tax)
}


D2Plot <- function(DS_dataframe){
  require(tidyverse)
  
  ds_plot <-
    DS_dataframe %>% 
      arrange(desc(log2FoldChange)) %>% 
      ggplot(aes(x=reorder(Taxonomy, -log2FoldChange), 
                 y=log2FoldChange, color=Phylum)) + 
      geom_point(size=4.5) + 
      geom_hline(yintercept=0, linetype="dashed", color = "red") +
      theme(axis.text.x = element_text(angle = 90, size= 10, hjust = 1)) +
      annotate("text", x=Inf, y = Inf, label="GMO", size=6, vjust=1.5, hjust=1.2) +
      annotate("text", x=-Inf, y =-Inf, label="non-GMO", size=6, vjust=-0.4, hjust=-0.1) +
        theme_bw() +
        theme(strip.text.x = element_text(angle = 0, size = 9,hjust = 0.5, vjust = 0.5, color = "black"),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5,color = "black"),
          legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
          legend.title =element_text(size = 9, face = "bold", hjust = 0.5,color = "black"),
          legend.text = element_text(angle = 0, size = 8, hjust = 0, vjust = 0.5, color = "black"),
          axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5, color = "black"),
          axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5, color = "black"),
          axis.title.y = element_text(size = 8, face = "bold"),
          axis.title.x = element_blank(),
          strip.background = element_blank(), 
          legend.position = "bottom")
    
  return(ds_plot)
}

# Fungi
# Pre-Transplant
plot_des_fungi_soil_pre <-
subset_samples(physeq_fungi_ev,
                 Niche %in% c("Soil") & Status %in% c("Pre-Transplant")) %>%
  prune_taxa(taxa_sums(x = .) > 0, .) %>%
  phyloseq_to_deseq2(., ~ Treatment) %>% 
  D2Stats(physeq = physeq_fungi_ev, des_df = .) %>% 
  mutate(Taxonomy = gsub("\\s+$", "", Taxonomy)) %>%
  D2Plot(DS_dataframe = .)

plot_des_fungi_soil_pre


plot_des_fungi_rhizo_pre <-
  subset_samples(physeq_fungi_ev,
                 Niche %in% c("Rhizosphere") & Status %in% c("Pre-Transplant")) %>%
  prune_taxa(taxa_sums(x = .) > 0, .) %>%
  phyloseq_to_deseq2(., ~ Treatment) %>% 
  D2Stats(physeq = physeq_fungi_ev, des_df = .) %>% 
  mutate(Taxonomy = gsub("\\s+$", "", Taxonomy)) %>%
  D2Plot(DS_dataframe = .)

plot_des_fungi_rhizo_pre

plot_des_fungi_root_pre <-
  subset_samples(physeq_fungi_ev,
                 Niche %in% c("Root") & Status %in% c("Pre-Transplant")) %>%
  prune_taxa(taxa_sums(x = .) > 0, .) %>%
  phyloseq_to_deseq2(., ~ Treatment) %>% 
  D2Stats(physeq = physeq_fungi_ev, des_df = .) %>% 
  mutate(Taxonomy = gsub("\\s+$", "", Taxonomy)) %>%
  D2Plot(DS_dataframe = .)

plot_des_fungi_root_pre

# Post-Transplant
plot_des_fungi_soil_post <-
  subset_samples(physeq_fungi_ev,
                 Niche %in% c("Soil") & Status %in% c("Post-Transplant")) %>%
  prune_taxa(taxa_sums(x = .) > 0, .) %>%
  phyloseq_to_deseq2(., ~ Treatment) %>% 
  D2Stats(physeq = physeq_fungi_ev, des_df = .) %>% 
  mutate(Taxonomy = gsub("\\s+$", "", Taxonomy)) %>%
  D2Plot(DS_dataframe = .)

plot_des_fungi_soil_post

plot_des_fungi_rhizo_post <-
  subset_samples(physeq_fungi_ev,
                 Niche %in% c("Rhizosphere") & Status %in% c("Post-Transplant")) %>%
  prune_taxa(taxa_sums(x = .) > 0, .) %>%
  phyloseq_to_deseq2(., ~ Treatment) %>% 
  D2Stats(physeq = physeq_fungi_ev, des_df = .) %>% 
  mutate(Taxonomy = gsub("\\s+$", "", Taxonomy)) %>%
  D2Plot(DS_dataframe = .)

plot_des_fungi_rhizo_post

plot_des_fungi_root_post <-
  subset_samples(physeq_fungi_ev,
                 Niche %in% c("Root") & Status %in% c("Post-Transplant")) %>%
  prune_taxa(taxa_sums(x = .) > 0, .) %>%
  phyloseq_to_deseq2(., ~ Treatment) %>% 
  D2Stats(physeq = physeq_fungi_ev, des_df = .) %>% 
  mutate(Taxonomy = gsub("\\s+$", "", Taxonomy)) %>%
  D2Plot(DS_dataframe = .)

plot_des_fungi_root_post

plot_des_fungi_inflo_post <-
  subset_samples(physeq_fungi_ev,
                 Niche %in% c("Inflorescence") & Status %in% c("Post-Transplant")) %>%
  prune_taxa(taxa_sums(x = .) > 0, .) %>%
  phyloseq_to_deseq2(., ~ Treatment) %>% 
  D2Stats(physeq = physeq_fungi_ev, des_df = .) %>% 
  mutate(Taxonomy = gsub("\\s+$", "", Taxonomy)) %>%
  D2Plot(DS_dataframe = .)

plot_des_fungi_inflo_post

plot_des_fungi_leaf_post <-
  subset_samples(physeq_fungi_ev,
                 Niche %in% c("Leaf") & Status %in% c("Post-Transplant")) %>%
  prune_taxa(taxa_sums(x = .) > 0, .) %>%
  phyloseq_to_deseq2(., ~ Treatment) %>% 
  D2Stats(physeq = physeq_fungi_ev, des_df = .) %>% 
  mutate(Taxonomy = gsub("\\s+$", "", Taxonomy)) %>%
  D2Plot(DS_dataframe = .)

plot_des_fungi_leaf_post


# Bacteria
# Pre-Transplant
plot_des_bact_soil_pre <-
  subset_samples(physeq_bact_ev,
                 Niche %in% c("Soil") & Status %in% c("Pre-Transplant")) %>%
  prune_taxa(taxa_sums(x = .) > 0, .) %>%
  phyloseq_to_deseq2(., ~ Treatment) %>% 
  D2Stats(physeq = physeq_bact_ev, des_df = .) %>% 
  mutate(Taxonomy = gsub("\\s+$", "", Taxonomy)) %>%
  D2Plot(DS_dataframe = .)

plot_des_bact_soil_pre


plot_des_bact_rhizo_pre <-
  subset_samples(physeq_bact_ev,
                 Niche %in% c("Rhizosphere") & Status %in% c("Pre-Transplant")) %>%
  prune_taxa(taxa_sums(x = .) > 0, .) %>%
  phyloseq_to_deseq2(., ~ Treatment) %>% 
  D2Stats(physeq = physeq_bact_ev, des_df = .) %>% 
  mutate(Taxonomy = gsub("\\s+$", "", Taxonomy)) %>%
  D2Plot(DS_dataframe = .)

plot_des_bact_rhizo_pre

plot_des_bact_root_pre <-
  subset_samples(physeq_bact_ev,
                 Niche %in% c("Root") & Status %in% c("Pre-Transplant")) %>%
  prune_taxa(taxa_sums(x = .) > 0, .) %>%
  phyloseq_to_deseq2(., ~ Treatment) %>% 
  D2Stats(physeq = physeq_bact_ev, des_df = .) %>% 
  mutate(Taxonomy = gsub("\\s+$", "", Taxonomy)) %>%
  D2Plot(DS_dataframe = .)

plot_des_bact_root_pre

# Post-Transplant
plot_des_bact_soil_post <-
  subset_samples(physeq_bact_ev,
                 Niche %in% c("Soil") & Status %in% c("Post-Transplant")) %>%
  prune_taxa(taxa_sums(x = .) > 0, .) %>%
  phyloseq_to_deseq2(., ~ Treatment) %>% 
  D2Stats(physeq = physeq_bact_ev, des_df = .) %>% 
  mutate(Taxonomy = gsub("\\s+$", "", Taxonomy)) %>%
  D2Plot(DS_dataframe = .)

plot_des_bact_soil_post

plot_des_bact_rhizo_post <-
  subset_samples(physeq_bact_ev,
                 Niche %in% c("Rhizosphere") & Status %in% c("Post-Transplant")) %>%
  prune_taxa(taxa_sums(x = .) > 0, .) %>%
  phyloseq_to_deseq2(., ~ Treatment) %>% 
  D2Stats(physeq = physeq_bact_ev, des_df = .) %>% 
  mutate(Taxonomy = gsub("\\s+$", "", Taxonomy)) %>%
  D2Plot(DS_dataframe = .)

plot_des_bact_rhizo_post

plot_des_bact_root_post <-
  subset_samples(physeq_bact_ev,
                 Niche %in% c("Root") & Status %in% c("Post-Transplant")) %>%
  prune_taxa(taxa_sums(x = .) > 0, .) %>%
  phyloseq_to_deseq2(., ~ Treatment) %>% 
  D2Stats(physeq = physeq_bact_ev, des_df = .) %>% 
  mutate(Taxonomy = gsub("\\s+$", "", Taxonomy)) %>%
  D2Plot(DS_dataframe = .)

plot_des_bact_root_post

plot_des_bact_inflo_post <-
  subset_samples(physeq_bact_ev,
                 Niche %in% c("Inflorescence") & Status %in% c("Post-Transplant")) %>%
  prune_taxa(taxa_sums(x = .) > 0, .) %>%
  phyloseq_to_deseq2(., ~ Treatment) %>% 
  D2Stats(physeq = physeq_bact_ev, des_df = .) %>% 
  mutate(Taxonomy = gsub("\\s+$", "", Taxonomy)) %>%
  D2Plot(DS_dataframe = .)

plot_des_bact_inflo_post

plot_des_bact_leaf_post <-
  subset_samples(physeq_bact_ev,
                 Niche %in% c("Leaf") & Status %in% c("Post-Transplant")) %>%
  prune_taxa(taxa_sums(x = .) > 0, .) %>%
  phyloseq_to_deseq2(., ~ Treatment) %>% 
  D2Stats(physeq = physeq_bact_ev, des_df = .) %>% 
  mutate(Taxonomy = gsub("\\s+$", "", Taxonomy)) %>%
  D2Plot(DS_dataframe = .)

plot_des_bact_leaf_post

# Testing using corncob
library(corncob)

subset_samples(physeq_fungi_ev,
               Niche %in% c("Root") & Status %in% c("Post-Transplant")) %>%
  prune_taxa(taxa_sums(x = .) > 0, .) %>%
  differentialTest( data = .,
                    formula = ~ Treatment,
                    phi.formula = ~ Treatment,
                    formula_null = ~ 1,
                    phi.formula_null = ~ Treatment,
                    test = "Wald", 
                    boot = FALSE,
                    fdr_cutoff = 0.05) 

plot_des_fungi_soil_post <-
  subset_samples(physeq_fungi_ev,
                 Niche %in% c("Soil") & Status %in% c("Post-Transplant")) %>%
  prune_taxa(taxa_sums(x = .) > 0, .) %>%
  phyloseq_to_deseq2(., ~ Treatment) %>% 
  D2Stats(physeq = physeq_fungi_ev, des_df = .) %>% 
  mutate(Taxonomy = gsub("\\s+$", "", Taxonomy)) %>%
  D2Plot(DS_dataframe = .)

des_fungi_soil_pre <-
  subset_samples(physeq_fungi_ev,
                 Niche %in% c("Soil") & Status %in% c("Pre-Transplant")) %>%
  prune_taxa(taxa_sums(x = .) > 0, .) %>%
  phyloseq_to_deseq2(., ~ Treatment)

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

des_fungi_soil_pre

geoMeans = apply(counts(des_fungi_soil_pre), 1, gm_mean)
diagdds = estimateSizeFactors(des_fungi_soil_pre, geoMeans = geoMeans)
diagdds = DESeq(des_fungi_soil_pre, fitType="parametric")

# Calculate results (p < 0.01)
results <- results(diagdds, alpha = 0.01, 
                   contrast = c("Treatment", "GMO","non-GMO"))

# Tabulate results
sigtab_results <- results[which(results$padj < 0.05), ]
summary(sigtab_results)

sigtab_results_tax <- cbind(
  as(sigtab_results, "data.frame"), 
  as(tax_table(physeq_fungi_ev)[row.names(sigtab_results), ], "matrix"))

sigtab_results_tax %>% 
  mutate(Taxonomy = gsub("\\s+$", "", Taxonomy)) %>% 
  arrange(log2FoldChange) %>% 
  ggplot(aes(x=reorder(Taxonomy, desc(log2FoldChange)), y=log2FoldChange, color=Phylum)) + 
  geom_point(size=4.5) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  themePlain() +
  theme(axis.text.x = element_text(angle = 90, size= 10, hjust = 1)) +
  annotate("text", x=Inf, y = Inf, label="GMO", size=6, vjust=1.5, hjust=1.2) +
  annotate("text", x=-Inf, y =-Inf, label="non-GMO", size=6, vjust=-0.4, hjust=-0.1)

sigtab_results_tax %>% 
  mutate(Taxonomy = gsub("\\s+$", "", Taxonomy))



# position of the annotation layer ---------------------------------------------
# Annotation position and hjust, vjust use to modify it
# xpos ypos        annotateText        hjust    vjust
# 1 -Inf -Inf Bottom Left (h0,v0)        0        0
# 2 -Inf  Inf    Top Left (h0,v1)        0        1
# 3  Inf -Inf  Bottom Right h1,v0        1        0
# 4  Inf  Inf     Top Right h1,v1        1        1


# ***********************************************-------------------------------
# ABUNDANCE TABLES AND STACKED BARS --------------------------------------------
StackedBarsDF <- function(dataframe, taxonomy, metadata, Rank){
  require(tidyverse)
  Var <- enquo(Rank)
  
  # Transpose to rows being the ASV, join metadata and taxonomy.
  # Sum up by Phylum and generate a wide data frame as the facets.
  wide_df <-
    dataframe %>% 
    t() %>%
    as.data.frame() %>%
    left_join(x = rownames_to_column(., var = "ASV_ID"),
              y = taxonomy,
              by = "ASV_ID") %>%
    dplyr::select(.data = ., starts_with("TR"), !!Var) %>%
    group_by(!!Var) %>% 
    pivot_longer(cols=c(-!!Var),names_to="Sample_ID") %>% 
    left_join(x=., 
              y= metadata,
              by="Sample_ID") %>% 
    na.omit() %>% 
    dplyr::select(!!Var, Niche, Status, Treatment, value) %>% 
    group_by(!!Var, Niche, Status, Treatment) %>% 
    summarise(across(c(value), list(sum = sum))) %>% 
    ungroup() %>% 
    unite("Group", Niche, Status, Treatment, sep = "_") %>% 
    pivot_wider(names_from = Group, values_from = value_sum) 
  
  # Use map to divide each cell by the column sum to scale to 0-100
  long_df <-
    wide_df %>% 
    map_if(is.numeric, ~./sum(.)) %>%
    as_tibble() %>% 
    pivot_longer(cols = -!!Var, names_to = "Group", values_to = "RelAbund") %>% 
    separate(Group, c("Niche", "Status","Treatment"),sep = "_") 

  return(long_df)
}

StackedBarsDF(fungi_ev_100itr, fungi_taxonomy_fix, meta_fungi_100itr, Phylum)


PlotBar <- function(dataframe){
  require(tidyverse)
  
  bar_plot <- 
    dataframe %>% 
    ggplot(aes(x=Treatment, y=RelAbund, fill=Phylum)) +
    geom_bar(stat = "identity") +
    facet_wrap(Status~Niche, ncol = 5, nrow = 2) +
    theme_bw() +
    theme(strip.text.x = element_text(angle = 0, size = 9,hjust = 0.5, vjust = 0.5, color = "black"),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5,color = "black"),
          legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
          legend.title =element_text(size = 9, face = "bold", hjust = 0.5,color = "black"),
          legend.text = element_text(angle = 0, size = 8, hjust = 0, vjust = 0.5, color = "black"),
          axis.text.x = element_text(angle = 33, size = 8, hjust = 1, vjust = 1.05, color = "black"),
          axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5, color = "black"),
          axis.title = element_text(angle = 0, size = 8, face = "bold"),
          strip.background = element_blank(), 
          legend.position = "bottom") 
  
  return(bar_plot)
}


fungi_phylum_relab <-
  StackedBarsDF(fungi_ev_100itr, fungi_taxonomy_fix, meta_fungi_100itr, Phylum) 
fungi_phylum_relab

fungi_phylum_relab %>% 
  subset(Phylum == "Glomeromycota")

write.csv(fungi_phylum_relab, "fungi_phylum.csv")

bact_phylum_relab <-
  StackedBarsDF(bact_ev_100itr, bact_taxonomy_fix, meta_bact_100itr, Phylum) 
bact_phylum_relab

write.csv(bact_phylum_relab, "bacteria_phylum.csv")


# Plotting the bargraphs -------------------------------------------------------
fungi_phylum_relab %>% 
  as.data.frame() %>% 
  PlotBar() + 
  scale_fill_manual(values = paletteCB15) 

# Just the first 15 Phyla, the rest is grouped into others.
bacteria_phyla <-
  bact_phylum_relab %>% 
  dplyr::select(Phylum, RelAbund) %>% 
  group_by(Phylum) %>% 
  summarise(across(c(RelAbund), list(Sum = sum))) %>%
  arrange(desc(RelAbund_Sum)) %>% 
  top_n(15) %>% 
  pull(Phylum)

bact_phylum_relab %>% 
  mutate(
    Phylum = ifelse(Phylum%in%bacteria_phyla, Phylum, "Others")) %>% 
  as.data.frame() %>% 
  PlotBar() + 
  scale_fill_manual(values = paletteCB15) 

