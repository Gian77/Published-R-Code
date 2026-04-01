
# *********************************************---------------------------------
# ************ USEARCH ANALYSIS ***************---------------------------------
#
# Manuscript:   The impact of drought on microbiome assembly in Arabidopsis
# Authors:      ...
# Affiliation:  1 Department of Plant Soil and Microbial Sciences, Michigan State University, East Lansing MI 48824
#               2 Great Lakes Bioenergy Research Center, Michigan State University, East Lansing MI 48824
#               ...
#
# Journal:      Nature 
# Date:         Nov 02, 2023
#
# *********************************************---------------------------------

# R options --------------------------------------------------------------------
options(scipen = 9999, pillar.sigfig=6, digits=6, max.print=10000000) 
#rm(list = ls())

# Check the lib paths ----------------------------------------------------------
#.libPaths()
#.libPaths("/home/gian/R/x86_64-pc-linux-gnu-library/4.4")

# Define the target version path you want to use
#desired_lib_path <- "/home/gian/R/x86_64-pc-linux-gnu-library/4.4"

# Check if the current .libPaths() points to 4.5
#current_lib_paths <- .libPaths()
#if (any(grepl("4\\.5", current_lib_paths))) {
#  message("Detected R 4.5 lib path. Resetting to R 4.4 lib path...")
#  .libPaths(desired_lib_path)
#} else {
#  message("No R 4.5 lib path detected. No change made.")
#}

# Optional: print the active library paths
print(.libPaths())

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
  lindia,
  LinDA,
  tm,
  topicmodels,
  ldatuning,
  install=TRUE
)
               

# Solve known conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("intersect", "base")
conflict_prefer("annotate", "ggplot2")
conflict_prefer("phyloseq", "tax_glom")
conflicts_prefer(phyloseq::tax_glom)
conflicts_prefer(base::unname)
conflicts_prefer(ggplot2::annotate)
conflicts_prefer(dplyr::filter)

# Import tables ----------------------------------------------------------------
asv_table <-
  read.csv("dataset/R1/otutable_UNOISE_230bp.txt", row.names = 1, sep = "\t")
head(asv_table)
colnames(asv_table)
dim(asv_table)

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

# check blast to Syncom
Fig_S1_identity_coverage <- 
  ggarrange(
    data.frame(identity = constax_asv$Isolate_percent_id) %>% 
      ggplot(aes(x=identity)) +
      geom_histogram() + 
      theme_bw() +
      labs(title = "ASV SynCom\nsequence BLAST identity",
           x = "Identity %",
           y = NULL),
    data.frame(coverage = constax_asv$HL_hit_query_cover) %>% 
      ggplot(aes(x=coverage)) +
      geom_histogram() +
      theme_bw() +
      labs(title = "ASV SynCom\nsequence BLAST coverage",
           x = "Coverage %",
           y = NULL))

Fig_S1_identity_coverage

# ***** Figure S1 ***** --------------------------------------------------------
ggsave("figures/Fig_S1_blast_to_syncom.pdf", plot = Fig_S1_identity_coverage, device = "pdf")

# Filtering taxonomy -----------------------------------------------------------
ReformatTaxonomy <- function(taxonomy_tab){
  require(tidyverse)
  
  lastValue <- function(x) tail(x[!is.na(x)], 1)
  last_taxons<- apply(taxonomy_tab[,c(2:7)], 1, lastValue)
  taxonomy_tab$BestMatch <- last_taxons
  taxonomy_res <-
    taxonomy_tab %>% 
    unite(ASV_ID, BestMatch, col=Taxonomy, sep = "-", remove = FALSE)
  
  return(taxonomy_res)
}


any(constax_asv$family=="Mitochondria")

constax_filt <-
  constax_asv %>%
  filter(Isolate_percent_id >= 60) %>% 
  filter(Isolate_query_cover >= 60) %>%
  filter(!order %in% "Chloroplast") %>% 
  filter(!family %in% "Mitochondria") %>% 
  rownames_to_column("ASV_ID") %>%
  separate(Isolate, c("Kingdom", "Phylum", "Class", "Order", 
                      "Family", "Genus", "Isolate"), sep = ";") %>% 
  mutate(Isolate = gsub("^_", "", Isolate)) %>% 
  dplyr::select(ASV_ID, domain, phylum, class, order, family, genus, Isolate) %>% 
  mutate(across(c("domain", "phylum", "class","order",
                  "family","genus", "Isolate"), ~na_if(., ""))) %>% 
  ReformatTaxonomy() %>% 
  column_to_rownames("ASV_ID")


head(constax_filt)
dim(constax_filt)

# NOTE. As we discussed with Sarah and Brittni we are using a 99% identity to the 
# Syncom. This will reduce the ASVs to less than 100.

# Importing sequences ----------------------------------------------------------
usearch_asv <- 
  readDNAStringSet("dataset/R1/asv_230bp.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
usearch_asv

# ***********************************************-------------------------------
# GENERATE PHYLOSEQ OBJECT -----------------------------------------------------

# Filter to match the SynCom ---------------------------------------------------
physeq_all <-
  phyloseq(
    otu_table(asv_table, taxa_are_rows = TRUE),
    sample_data(metadata),
    tax_table(as.matrix(constax_filt)),
    usearch_asv) %>% 
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.)

physeq_all

# ***********************************-------------------------------------------
# DECONTAMINATION --------------------------------------------------------------
sort(sample_data(physeq_all)$Genotype)

# Adding control variable
sample_data(physeq_all)$is.neg <-
  sample_data(physeq_all)$Genotype == "pbs" |
  sample_data(physeq_all)$Genotype == "RTSF" |
  sample_data(physeq_all)$Genotype == "kit" 

head(physeq_all@sam_data)

# Detecting contaminant frequencies
contam_asv <-
  decontam::isContaminant(
    physeq_all,
    method = "prevalence",
    neg = "is.neg",
    threshold = 0.1)

contam_asv
table(contam_asv$contaminant) 
contam_asv %>% filter(contaminant == TRUE)

# Check contaminants taxonomy
contam_asv %>% 
  filter(contaminant == TRUE) %>% 
  rownames_to_column("ASV_ID") %>% 
  left_join(constax_asv %>% 
              rownames_to_column("ASV_ID"),
            by="ASV_ID")

# Check contaminants abundance
asv_table %>% 
  mutate(Sum = rowSums(across(is.numeric))) %>%
  rownames_to_column("asvID") %>% 
  right_join(contam_asv %>%
             filter(contaminant == TRUE) %>% 
             rownames_to_column("asvID"), 
             by = "asvID") %>% 
  select(!starts_with("sample"))


# Removing contaminants --------------------------------------------------------
remove_taxa = function(physeq, badTaxa) {
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

physeq_filt <-
  remove_taxa(physeq_all, 
              rownames(subset(contam_asv, contaminant  %in% c("TRUE")))) %>% 
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.)

physeq_filt
head(physeq_filt@sam_data)

# ***********************************************-------------------------------
# DATA NORMALIZATION -----------------------------------------------------------

# I think rarefying the data is not feasible unless we consider a removal of 
# a large number of samples. Alternative approach is to remove all samples below
# 10 reads and then use metagenomeSeq normalization method.

# https://bioconductor.org/packages/release/bioc/vignettes/metagenomeSeq/inst/doc/metagenomeSeq.pdf

CSSNorm <-function(physeq){
  require(metagenomeSeq)
  physeq %>% 
    phyloseq_to_metagenomeSeq() -> physeq_CSS
  #p_biom <-cumNormStatFast(physeq_CSS)
  p_biom <-cumNormStat(physeq_CSS)
  biom_quant <-cumNorm(physeq_CSS, p=p_biom)
  physeq_CSS <- MRcounts(biom_quant, norm=T)
  physeq_mSeq <- physeq
  otu_table(physeq_mSeq) <- otu_table(physeq_CSS, taxa_are_rows=TRUE)
  return(physeq_mSeq)
}

physeq_filt_mSeq <- 
  physeq_filt %>% 
  prune_samples(sample_sums(x=.) > 10, x =.) %>% # removing samples with less than 10 reads
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  CSSNorm()

physeq_filt_mSeq
head(physeq_filt_mSeq@otu_table)

write.csv(
  otu_table(physeq_filt_mSeq) %>% as.data.frame(),
  "tables/physeq_filt_mSeq_asv.csv"
)
write.csv(
  tax_table(physeq_filt_mSeq) %>% as.matrix() %>% as.data.frame(),
  "tables/physeq_filt_mSeq_taxaonomy.csv"
)

# Check the distributions of corrected read-counts per ASV across samples to 
# determine how good was the MetagenomeSeq normalization. I am also adding the 
# simple Hellinger transformation as a comparison. 
# See Ecologically meaningful transformations for ordination of species data
# https://link.springer.com/article/10.1007/s004420100716

# this is the Hellinger transforme data
physeq_sqrt <- 
  transform_sample_counts(physeq_filt %>%
                          prune_samples(sample_sums(x = .) > 10, x = .) %>%
                          prune_taxa(taxa_sums(x = .) > 0, x = .), 
                          function(x) sqrt(x / sum(x)))
head(physeq_sqrt@otu_table)

df_mSeq_filt <- 
  data.frame(
    counts = sort(sample_sums(physeq_filt %>%
                           prune_samples(sample_sums(x = .) > 10, x = .) %>%
                           prune_taxa(taxa_sums(x = .) > 0, x = .))),
    Hellinger = sort(sample_sums(physeq_sqrt)),
    mSeq = sort(sample_sums(physeq_filt_mSeq))) 

head(df_mSeq_filt)

# Calculating variance across normalization methods

# ***** Table S1 ***** ---------------------------------------------------------
Table_S1_compare_norm_methods <-
  df_mSeq_filt %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "method", values_to = "count") %>%
  group_by(method) %>%
  summarise(across(c(count), list(mean = mean, var = var))) %>%
  mutate(count_CV = count_var / count_mean) %>%
  as.data.frame()

write.csv(Table_S1_compare_norm_methods,
          "tables/Table_S1_compare_norm_methods.csv")

# plotting
Fig_S2_sample_normal_compare <-
  ggarrange(
    # Hellinger plot
    df_mSeq_filt %>%
      select(Hellinger) %>%
      rownames_to_column("SampleID") %>%
      pivot_longer(-SampleID) %>%
      ggplot(aes(x = SampleID, y = value)) +
      geom_bar(stat = "identity", position = "dodge", fill = "lightblue") +
      geom_text(
        data = Table_S1_compare_norm_methods %>% filter(method == "Hellinger"),
        aes(x = 10, y = max(df_mSeq_filt$Hellinger) * 0.9, 
          label = paste0("CV = ", round(count_CV, 2))),
        inherit.aes = FALSE, hjust = 0,color = "blue",size = 4) +
      theme_classic() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5)) +
      labs(title = "Hellinger", y = "Counts"),
    
    # mSeq plot
    df_mSeq_filt %>%
      select(mSeq) %>%
      rownames_to_column("SampleID") %>%
      pivot_longer(-SampleID) %>%
      ggplot(aes(x = SampleID, y = value)) +
      geom_bar(stat = "identity", position = "dodge", fill = "orange") +
      geom_text(
        data = Table_S1_compare_norm_methods %>% filter(method == "mSeq"),
        aes(x = 10,y = max(df_mSeq_filt$mSeq) * 0.9,
          label = paste0("CV = ", round(count_CV, 2))),
        inherit.aes = FALSE, hjust = 0,color = "darkorange",size = 4) +
      theme_classic() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5)) +
      labs(title = "MetagenomeSeq", y = "Counts"),
    
    # counts plot
    df_mSeq_filt %>%
      select(counts) %>%
      rownames_to_column("SampleID") %>%
      pivot_longer(-SampleID) %>%
      ggplot(aes(x = SampleID, y = value)) +
      geom_bar(stat = "identity", position = "dodge", fill = "lightgreen") +
      geom_text(
        data = Table_S1_compare_norm_methods %>% filter(method == "counts"),
        aes(x = 10,y = max(df_mSeq_filt$counts) * 0.9,
          label = paste0("CV = ", round(count_CV, 2))),
        inherit.aes = FALSE,hjust = 0,color = "darkgreen",size = 4) +
      theme_classic() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5)) +
      labs(title = "Raw counts", y = "Counts"),
    
    ncol = 1,
    nrow = 3,
    align = "hv"
  )


Fig_S2_sample_normal_compare

# ***** Figure S2 ***** --------------------------------------------------------
ggsave("figures/Fig_S2_sample_normal_compare.pdf",
       plot = Fig_S2_sample_normal_compare,
       device = "pdf")

# SUBSET PHYLOSEQ DATA ---------------------------------------------------------
# subset to just the SynCom data. I will keep going with the just the ASVs.

physeq_sc <-
  physeq_sqrt %>% 
  subset_samples(Inoculation %in% "sc") %>% 
  subset_samples(Genotype %in% c("col", "cp1", "ed4",
                                 "nd1", "pd4", "sd2" )) %>%
  subset_samples(Treatment %in% c("w", "d")) %>% 
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.)

physeq_sc

# usingspeedyseq to add directly to the sam_data with dplyr wrappers
physeq_sc <- 
  physeq_sc %>%
  mutate_sample_data(
    Replicate = as.factor(Replicate),
    Genotype = as.factor(Genotype),
    Treatment = as.factor(Treatment),
    Genotype = fct_recode(Genotype,
                          "Col-0" = "col", "cpr1" = "cp1", "pad4" = "pd4",
                          "sid2" = "sd2", "ndr1" = "nd1", "erd4" = "ed4"),
    GenotypeTreatment = paste(Genotype, Treatment, sep = "_"),
    Treatment = fct_relevel(Treatment, "w","d"),
    Genotype = fct_relevel(Genotype,"Col-0","cpr1","pad4","sid2","ndr1","erd4"),
    GenotypeTreatment = fct_relevel(GenotypeTreatment,
                                    "Col-0_w","Col-0_d","cpr1_w","cpr1_d","pad4_w","pad4_d",
                                    "sid2_w","sid2_d","ndr1_w","ndr1_d","erd4_w","erd4_d")
    )

physeq_sc@sam_data %>%as_tibble()
physeq_sc@sam_data$Treatment
physeq_sc@sam_data$GenotypeTreatment


# Extracting no plant data for Sarah's plots -----------------------------------
physeq_all@sam_data

physeq_no_plant <-
  physeq_all %>% 
  subset_samples(Treatment %in% "p" | Genotype %in% c("inoc", "np" )) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.)

physeq_no_plant


# this is the Hellinger transforme data
physeq_no_plant_sqrt <- 
  transform_sample_counts(physeq_no_plant %>%
                            prune_samples(sample_sums(x = .) > 10, x = .) %>%
                            prune_taxa(taxa_sums(x = .) > 0, x = .), 
                          function(x) sqrt(x / sum(x)))

head(physeq_no_plant_sqrt@otu_table)

physeq_no_plant_sqrt@sam_data
write.csv(as.data.frame(as.matrix(physeq_no_plant_sqrt@sam_data)), 
          "tables/physeq_noSC_samples.csv")

physeq_no_plant_sqrt@otu_table
write.csv(as.data.frame(as.matrix(physeq_no_plant_sqrt@otu_table)), 
          "tables/physeq_noSC_otu_table.csv")


#***********************************************--------------------------------
# BETA DIVERSITY ---------------------------------------------------------------

# Extract dataset --------------------------------------------------------------
asv_sc <- t(as.data.frame(otu_table(physeq_sc)))
metadata_sc <- as.data.frame(as.matrix(sample_data(physeq_sc)))
taxonomy_sc <- as.data.frame(as.matrix(tax_table(physeq_sc)))

write.csv(asv_sc, "tables/normalized_asv_table.csv")
write.csv(pcoa_hell_df, "tables/metadata.csv")
write.csv(taxonomy_sc, "tables/taxonomy.csv")

# Euclidean distance -----------------------------------------------------------
# using Euclidean distance on Hellinger-transformed data is the most 
# straightforward and commonly used approach

pcoa_hell <- 
  cmdscale(vegdist(asv_sc, method = "euclidean"), k = 2, eig=TRUE)  

# Variance explained 
expvar_hell <- pcoa_hell$eig / sum(pcoa_hell$eig)
expvar_hell[1:2]
# [1] 0.170552 0.110547

pcoa_hell_df <- 
  as.data.frame(pcoa_hell$points) %>% 
  rename(Axis.1 = V1 , Axis.2 = V2) %>%
  rownames_to_column("SampleID") %>% 
  full_join(metadata_sc %>% 
              rownames_to_column("SampleID"),
            by="SampleID")

pcoa_hell_df

write.csv(pcoa_hell_df, "tables/pcoa_coordinates.csv")

# plot PCoA ordination ----------------------------------
# Treatment by replicate based on euclidean distance on Hellinger transformed data.
fig_pcoa_hell <-
  pcoa_hell_df %>% 
  ggplot(aes(x=Axis.1,y=Axis.2, color = Treatment, shape = Replicate)) +
  geom_point(size = 2) +
  theme_bw() +
  #stat_ellipse() +
  theme(plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_markdown(size = 10, hjust = 0.5),
        strip.text = element_markdown(size = 10, face = "bold"),
        axis.text.x = element_markdown(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
        axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
        legend.position = "right", 
        legend.title = element_text(hjust= 0.5)) +
  guides(color = guide_legend(ncol=1, override.aes = list(shape = 15, size = 3.5)),
         shape = guide_legend(ncol=1, override.aes = list(color = "black", size=2.5)),
         size = guide_legend(ncol=1)) +
  coord_cartesian(xlim = c(-0.5, 0.4),ylim = c(-0.3, 0.6)) +
  labs(title = "PCoA", 
       subtitle = "Euclidean distance on<br>Hellinger transformed ASV counts",
       x = "Axis.1 [17.06%]",
       y = "Axis.2 [11.05%]") + 
  scale_color_manual(values=c('#f46d43', '#74add1'))

fig_pcoa_hell


# Treatment by Genotype if it make sense 
palette_cb12 <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288","#404040",
                  "#6699CC","#661100", "#999933", "#44AA99","#AA4499", "#888888")
pie(rep(1, 12), col = palette_cb12)

palette_cb12 <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf',
                  '#e0f3f8','#abd9e9','#74add1','#4575b4','#313695','#1A1447')
pie(rep(1, 12), col = palette_cb12)

palette_cb12_reversed <- c('#a50026', '#e0f3f8', '#d73027', '#abd9e9', '#f46d43', '#74add1',
                          '#fdae61', '#1A1447', '#fee090', '#313695', '#ffffbf', '#4575b4')
pie(rep(1, 12), col = palette_cb12_reversed)

fig_pcoa_hell_gen <-
  pcoa_hell_df %>% 
  ggplot(aes(x=Axis.1, y=Axis.2, color = interaction(Genotype, Treatment, sep="_"), shape = Treatment)) +
  geom_point(size = 2) +
  theme_bw() +
  #stat_ellipse() +
  theme(plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_markdown(size = 10, hjust = 0.5),
        strip.text = element_markdown(size = 10, face = "bold"),
        axis.text.x = element_markdown(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
        axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
        legend.position = "right", legend.title.align = 0.5) +
  guides(color = guide_legend(ncol=1, override.aes = list(shape = 15, size = 3.5)),
         shape = guide_legend(ncol=1, override.aes = list(color = "black", size=2.5)),
         size = guide_legend(ncol=1)) +
  coord_cartesian(xlim = c(-0.5, 0.4),ylim = c(-0.3, 0.6)) +
  labs(title = "PCoA", 
       subtitle = "Euclidean distance on<br>Hellinger transformed ASV counts",
       x = "Axis.1 [17.06%]",
       y = "Axis.2 [11.05%]",
       color = NULL, 
       shape = NULL) +
  scale_color_manual(values=palette_cb12)

fig_pcoa_hell_gen

# Jaccard distance -----------------------------------------------------------
pcoa_jacc <- 
  cmdscale(vegdist(asv_sc, method = "jaccard"), k = 2, eig=TRUE)  

# Variance explained 
expvar_jacc <- pcoa_jacc$eig / sum(pcoa_jacc$eig)
expvar_jacc[1:2]
# [1] 0.1358477 0.0782763

pcoa_jacc_df <- 
  as.data.frame(pcoa_jacc$points) %>% 
  rename(V1="Axis.1" , V2="Axis.2") %>%
  rownames_to_column("SampleID") %>% 
  full_join(metadata_sc %>% 
              rownames_to_column("SampleID"),
            by="SampleID") %>% 
  mutate(Genotype = fct_recode(Genotype,
                               "Col-0" = "col", "cpr1" = "cp1", "pad4" = "pd4",
                               "sid2" = "sd2", "ndr1" = "nd1", "erd4" = "ed4"))

pcoa_jacc_df

# plot PCoA ordination Treatment by Replicate
fig_pcoa_jac <-
  pcoa_jacc_df %>% 
  ggplot(aes(x=Axis.1,y=Axis.2, color = Treatment, shape = Replicate)) +
  geom_point(size = 2) +
  theme_bw() +
  #stat_ellipse() +
  theme(plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_markdown(size = 10, hjust = 0.5),
        strip.text = element_markdown(size = 10, face = "bold"),
        axis.text.x = element_markdown(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
        axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
        legend.position = "right", legend.title.align = 0.5) +
  guides(color = guide_legend(ncol=1, override.aes = list(shape = 15, size = 3.5)),
         shape = guide_legend(ncol=1, override.aes = list(color = "black", size=2.5)),
         size = guide_legend(ncol=1)) +
  #coord_cartesian(xlim = c(-0.5, 0.4),ylim = c(-0.3, 0.6)) +
  labs(title = "PCoA", 
       subtitle = "Jaccard distance on<br>Hellinger transformed ASV counts",
       x = "Axis.1 [13.58%]",
       y = "Axis.2 [7.83%]")+ 
  scale_color_manual(values=c('#f46d43', '#74add1'))

fig_pcoa_jac

# # Treatment by Genotype 
fig_pcoa_jac_jen <-
  pcoa_jacc_df %>% 
  ggplot(aes(x=Axis.1, y=Axis.2, color = interaction(Genotype, Treatment, sep="_"), shape = Treatment)) +
  geom_point(size = 2) +
  theme_bw() +
  theme(plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_markdown(size = 10, hjust = 0.5),
        strip.text = element_markdown(size = 10, face = "bold"),
        axis.text.x = element_markdown(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
        axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
        legend.position = "right", legend.title.align = 0.5) +
  guides(color = guide_legend(ncol=1, override.aes = list(shape = 15, size = 3.5)),
         shape = guide_legend(ncol=1, override.aes = list(color = "black", size=2.5)),
         size = guide_legend(ncol=1)) +
  #coord_cartesian(xlim = c(-0.5, 0.4),ylim = c(-0.3, 0.6)) +
  labs(title = "PCoA", 
       subtitle = "Jaccard distance on<br>Hellinger transformed ASV counts",
       x = "Axis.1 [13.58%]",
       y = "Axis.2 [7.83%]",
       color = NULL, 
       shape = NULL) +
  scale_color_manual(values=palette_cb12)

fig_pcoa_jac_jen


# ***** Figure 1 ***** ---------------------------------------------------------
Fig_1_pcoa_all <-
  ggarrange(
  ggarrange(fig_pcoa_hell, fig_pcoa_jac,
            ncol = 1, nrow = 2,
            common.legend = TRUE, 
            legend = "right"),
  ggarrange(fig_pcoa_hell_gen, fig_pcoa_jac_jen,
            ncol = 1, nrow = 2,
            common.legend = TRUE,
            legend = "right"),
  ncol = 2,
  nrow = 1
)

Fig_1_pcoa_all

ggsave("figures/Fig_1_pcoa_all.pdf",
       plot = Fig_1_pcoa_all,
       device = "pdf")


# to label the ellipses
pcoa_jacc_df %>%
  group_by(Treatment) %>%
  summarise(Axis.1 = mean(Axis.1), Axis.2 = mean(Axis.2))

pcoa_jacc_df %>% 
  ggplot(aes(x = Axis.1, y = Axis.2, color = Treatment, shape = Treatment)) +
  geom_point(size = 2) +
  stat_ellipse(aes(group = Treatment, color = Treatment)) +
  geom_text(data = centroids, aes(x = Axis.1, y = Axis.2, label = Treatment),
            size = 4, hjust = 0.5, vjust = -0.8) +
  theme_bw()


# ***********************************************-------------------------------
# PERMANOVA --------------------------------------------------------------------

# These below are the two models than make the make sense to me. Permutations are 
# restricted within the levels of Replicate. This means the test will only compare
# samples within the same replicate to assess the effect of Treatment or Genotype.
# I try inverting the factors but doea not change the results.

set.seed(991)

# Not strictly necessary I think, but I need to check!
perm_for_adonis <- how(nperm = 999)
setBlocks(perm_for_adonis) <- with(metadata_sc, Replicate)

adonis2_sc1 <-
  adonis2(formula = asv_sc ~ Treatment * Genotype, data = metadata_sc, 
          method = "euclidean", strata = metadata_sc$Replicate, by = "terms",
          permutations = perm_for_adonis)

adonis2_sc1

adonis2_sc2<-
  adonis2(formula = asv_sc ~ Genotype * Treatment, data = metadata_sc, 
          method = "euclidean", strata = metadata_sc$Replicate, by = "terms",
          permutations = perm_for_adonis)

adonis2_sc2

# ***** Table S2 ***** ---------------------------------------------------------
main_adonis_table <- 
rbind(
  as.data.frame(adonis2_sc1) %>% 
    mutate(model = "asv_sc ~ Treatment * Genotype",
           block = "Replicate"),
  as.data.frame(adonis2_sc2) %>% 
    mutate(model = "asv_sc ~ Genotype * Treatment",
           block = "Replicate")
  )

write.csv(main_adonis_table, "tables/main_adonis_table.csv")

# Exploring the significant Interaction -------------------------------------

# a. PCoA Plot: If you have already calculated the PCoA coordinates, you can plot
# the interaction between Genotype and Treatment on the first few PCoA axes.

Fig_S3_pcoa_interactions <-
pcoa_hell_df %>% 
  ggplot(aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(fill = interaction(Genotype, Treatment), shape = Treatment), 
             size = 2, color = "black") +
  stat_ellipse(aes(color = interaction(Genotype, Treatment))) +
  theme_bw() +
  facet_grid(~Genotype) +
  theme(plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_markdown(size = 10, hjust = 0.5),
        strip.text = element_markdown(size = 10, face = "bold"),
        axis.text.x = element_markdown(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
        axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        legend.key.height = unit(0.4, "cm"), 
        legend.key.width = unit(0.4, "cm"),
        legend.position = "right", 
        legend.title.align = 0.5) +
  guides(fill = "none", # Remove fill legend if not needed
         color = "none",
         shape = guide_legend(ncol = 1, override.aes = list(color = "black", size = 2.5)),
         size = guide_legend(ncol = 2)) +
  coord_cartesian(xlim = c(-0.5, 0.4), ylim = c(-0.3, 0.6)) +
  labs(title = "PCoA", 
       subtitle = "Euclidean distance<br>on Hellinger transformed ASV counts",
       x = "Axis.1 [17.06%]",
       y = "Axis.2 [11.05%]") +
  scale_fill_manual(values = palette_cb12) +
  scale_color_manual(values = palette_cb12) +
  scale_shape_manual(values = c(21, 24))

Fig_S3_pcoa_interactions

# ***** Figure S3 ***** --------------------------------------------------------
ggsave("figures/Fig_S3_pcoa_interactions.pdf",
       plot = Fig_S3_pcoa_interactions,
       device = "pdf")


# b. Boxplot of PCoA Axis Scores: You can also plot the distribution of PCoA scores 
# for each Genotype-Treatment combination.

Fig_S4_interaction_axis1 <-
pcoa_hell_df %>% 
  ggplot(aes(x = reorder(interaction(Genotype, Treatment, sep = "_")), 
             y = Axis.1, 
             fill = reorder(interaction(Genotype, Treatment, sep = "_")))) +
  geom_boxplot(show.legend = FALSE, outliers = TRUE, outlier.shape = 8) + # Disable legend for boxplot
  geom_jitter(aes(fill = reorder(interaction(Genotype, Treatment, sep = "_")))) + # Dummy layer for fill legend
  theme_bw() +
  theme(
    plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
    strip.text = element_markdown(size = 10, face = "bold"),
    axis.text.x = element_markdown(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
    axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
    legend.position = "right", legend.title.align = 0.5) +
  scale_fill_manual(values = palette_cb12_reversed) +
  scale_color_manual(values = palette_cb12_reversed) +
  guides(fill = guide_legend(override.aes = list(size = 3, shape = 22))) +
  labs(x = NULL, fill = NULL, title = "PCoA Axis 1: Genotype * Treatment Interaction")
               
Fig_S4_interaction_axis1

Fig_S4_interaction_axis2 <-
  pcoa_hell_df %>% 
  ggplot(aes(x = reorder(interaction(Genotype, Treatment, sep = "_")), 
             y = Axis.2, 
             fill = reorder(interaction(Genotype, Treatment, sep = "_")))) +
  geom_boxplot(show.legend = FALSE, outliers = TRUE, outlier.shape = 8) +
  geom_jitter(aes(fill = reorder(interaction(Genotype, Treatment, sep = "_")))) +
  theme_bw() +
  theme(
    plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
    strip.text = element_markdown(size = 10, face = "bold"),
    axis.text.x = element_markdown(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
    axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
    legend.position = "right", legend.title.align = 0.5) +
  scale_fill_manual(values = palette_cb12_reversed) +
  scale_color_manual(values = palette_cb12_reversed) +
  guides(fill = guide_legend(override.aes = list(size = 3, shape = 22))) +
  labs(x = NULL, fill=NULL,title = "PCoA Axis 2: Genotype * Treatment Interaction")

Fig_S4_interaction_axis2
    
# ***** Figure S4 ***** --------------------------------------------------------
Fig_S4_pcoa_boxplot_interactions <-
ggarrange(
  Fig_S3_pcoa_interactions,
  Fig_S4_interaction_axis1,
  Fig_S4_interaction_axis2, 
  ncol = 1,
  nrow = 3,
  align = "hv", 
  labels = c("a", "b", "c"),
  heights = c(1, 0.7, 0.7))

Fig_S4_pcoa_boxplot_interactions

ggsave("figures/Fig_S4_pcoa_boxplot_interactions.pdf",
       plot = Fig_S4_pcoa_boxplot_interactions,
       device = "pdf")


# If the boxes for different groups do not overlap much, it suggests that the groups
# differ significantly along that PCoA axis.

# ***********************************************-------------------------------
# PERMANOVA PLOTS --------------------------------------------------------------
str(adonis2_sc1)

permanova_plot <-
data.frame(group = c("Treatment", "Genotype", "Treatment:Genotype"),
           R2 = adonis2_sc1$R2[1:3],
           betadisp=c("*", NA, "*")) %>% 
  mutate(lls = ifelse(!is.na(betadisp), paste(round(R2,3), betadisp, sep = " "), round(R2,3))) %>% 
  mutate(Factor = fct_relevel(group,"Treatment", "Genotype", "Treatment:Genotype", )) %>% 
  ggplot(aes(x=Factor, y=R2)) +
  geom_bar(stat = "identity", fill = "grey80") +
  geom_text(aes(label = lls), 
            size = 4, 
            hjust= 0.2, vjust = -0.8, 
            angle=25,
            colour = "black") +
  theme_bw() +
  theme(plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
        strip.text = element_markdown(size = 10, face = "bold"),
        axis.text.x = element_markdown(color = "black",angle = 25, size = 10, hjust = 1, vjust = 1),
        axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        legend.title = element_blank(), 
        legend.text = element_markdown(face = "italic", size=8),
        legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.3, "cm"),
        legend.position = "none", legend.title.align = 0.5) +
  ylim(0, 0.05) +
  guides(color = guide_legend(ncol = 2),
         fill = guide_legend(ncol = 2)) +
  labs(title = "PERMANOVA", 
       y=bquote(italic(R^2))) 

permanova_plot

ggsave("figures/Fig_permanova_plot.pdf",
       plot = permanova_plot,
       device = "pdf")


# ***********************************************-------------------------------
# PAIR-WISE PERMANOVA ----------------------------------------------------------
head(metadata_sc)

# Calculate number of comparisons 
PairCalc <- function(n){
  n_comparison = n * (n - 1) / 2
  return(n_comparison)
}

# 6 Genotypes per 2 Treatments is 12 groups
PairCalc(12)

# Custom function for pairwise adonis2 
pairwise_permanova <- function(sp_matrix,
                               metadata,
                               Var, 
                               dist = "euclidean", 
                               strata = NULL,  # Default to NULL if no strata is provided
                               perm = 999) {
  require(vegan)
  require(tidyverse)
  
  # List contrasts: Create pairwise combinations of genotypes
  group_var <- metadata %>% 
    pull(Var) %>% 
    as.character() # Convert to character to avoid factor issues
  
  if (length(unique(group_var)) < 2) {
    stop("Error: Fewer than two unique levels in the grouping variable.")
  }
  
  # Generate all unique combinations of Var (6 genotypes -> 15 comparisons)
  groups <- as.data.frame(t(combn(unique(group_var), m = 2))) 
  
  contrasts <- data.frame(
    group1 = groups$V1, group2 = groups$V2,
    R2 = NA, F_value = NA, df1 = NA, df2 = NA, p_value = NA
  )
  
  # Loop through the unique pairwise comparisons
  for (i in seq(nrow(contrasts))) {
    # Subset the metadata and OTU table for the current pair of genotypes
    sp_subset <- group_var == contrasts$group1[i] | group_var == contrasts$group2[i] 
    contrast_matrix <- sp_matrix[sp_subset, ]  # Subset the OTU matrix
    metadata_subset <- metadata[sp_subset, ]  # Subset the metadata for the same samples
    
    # Conditional to check if strata is NULL
    if (is.null(strata)) {
      # If strata is NULL, just run adonis without strata
      fit <- vegan::adonis2(
        contrast_matrix ~ group_var[sp_subset],  # Model formula
        method = dist,                          # Distance method (e.g., Euclidean, Bray-Curtis)
        perm = perm,                            # Number of permutations
        parallel = 8                             # Parallel processing (optional)
      )
    } else {
      # If strata is provided, use it in the model
      strata_var <- metadata_subset[[strata]]
      fit <- vegan::adonis2(
        contrast_matrix ~ group_var[sp_subset],  # Model formula
        method = dist,                          # Distance method (e.g., Euclidean, Bray-Curtis)
        strata = strata_var,                    # Strata (blocking factor, e.g., Replicate)
        perm = perm,                            # Number of permutations
        parallel = 8                             # Parallel processing (optional)
      )
    }
    
    # Extract results from adonis2 output
    contrasts$R2[i] <- round(fit$R2[1], digits = 3)
    contrasts$F_value[i] <- round(fit[["F"]][1], digits = 3)
    contrasts$df1[i] <- fit$Df[1]
    contrasts$df2[i] <- fit$Df[2]
    contrasts$p_value[i] <- fit$`Pr(>F)`[1]
  }
  
  # Return results (including p-values, R2, F-values)
  return(contrasts)
}


pair_adonis2 <- 
  pairwise_permanova(asv_sc, metadata_sc, "GenotypeTreatment", strata = "Replicate", perm = 999)

# IMPORTANT NOTE! Adjusting the p_value depends on which comparisons we want to 
# perform. i.e. how many independent tests we do.

Table_S3_pairwise_adonis <-
pair_adonis2 %>% 
  mutate(Adj_p_value = p.adjust(p_value, method = "BH"))

Table_S3_pairwise_adonis

# ***** Table S3 ***** ---------------------------------------------------------

write.csv(Table_S3_pairwise_adonis,"tables/Table_S3_pairwise_adonis.csv")

# ***********************************************-------------------------------
# BETA-DISPERSION --------------------------------------------------------------

# Interesting fact! within a function it is better to specify the otu
# generated from a phyloeq object in multiple steps otherwise 
# does not really work if you include it in ().
# e.g. otu <- as.data.frame(as.matrix(t(otu_table(physeq))))

BetadispExtr <- function(sp_matrix,
                         metadata,
                         group_var,
                         distance_mat,
                         perm_num){
  
  require(tidyverse)
  require(vegan)
  
  # Set up the permutation structure
  perm <- how(nperm = perm_num)
  setBlocks(perm) <- with(metadata, Replicate) # keep the blocking
  
  # using centroids since I have normalized using Hellinger,
  # therefore I do not expect huge variances
  disp <- betadisper(vegan::vegdist(sp_matrix, method = distance_mat),
                     metadata[[group_var]], type = "centroid")
  
  # Perform permutation test 
  disp_perm <- vegan::permutest(disp, permutations = perm, pairwise = TRUE)
  
  # Calculate group sizes
  group_sizes <- table(metadata[[group_var]])
  
  # Compute degrees of freedom for pairwise comparisons
  pairwise_comparisons <- combn(names(group_sizes), 2, simplify = FALSE)
  
  # Create a data frame of degrees of freedom for each pairwise comparison
  df_results <- map(pairwise_comparisons, function(pair) {
    n1 <- group_sizes[pair[1]]
    n2 <- group_sizes[pair[2]]
    df_residual <- n1 + n2 - 2
    df_group <- 1
    tibble(group1 = pair[1],
           group2 = pair[2],
           df_group = df_group,
           df_residual = df_residual)
  }) %>% bind_rows()
  
  # Merge pairwise t-values, p-values, and degrees of freedom
  disp_pair <- 
    full_join(
      data.frame(t = abs(disp_perm$statistic[-1])) %>%
        rownames_to_column("name") %>%
        mutate(name = gsub(" \\(t)", "", name)) %>% 
        separate(name, into = c("group1", "group2"), sep = "-", remove = FALSE),
        data.frame(p_value = disp_perm$pairwise$permuted) %>%
        rownames_to_column("name"),
      by = "name"
    ) %>%
    left_join(df_results, by = c("group1", "group2")) # Add df information
  
  # Return the full list
  return(list(disp = disp, disp_perm = disp_perm, disp_pair = disp_pair))
}

# ***** Table S4 ***** ---------------------------------------------------------

Table_S4_betadisper_resutls_main <-
  rbind(
    as.data.frame(anova(BetadispExtr(asv_sc, metadata_sc, "Treatment", "euclidean", 999)$disp)) %>% 
      mutate(group = "Treatment"),
    as.data.frame(anova(BetadispExtr(asv_sc, metadata_sc, "Genotype", "euclidean", 999)$disp))%>% 
      mutate(group = "Genotype"),
    as.data.frame(anova(BetadispExtr(asv_sc, metadata_sc, "GenotypeTreatment", "euclidean", 999)$disp))%>% 
      mutate(group = "GenotypeTreatment")
    )


Table_S4_betadisper_resutls_main

write.csv(Table_S4_betadisper_resutls_main,"tables/Table_S4_betadisper_resutls_main.csv")


# ***** Table S5 ***** ---------------------------------------------------------
# The function fails becasue of the "-" in the genotype name Col-0. Need to 
# fix that so it will run correctly. I will made a fixed metadata for this analysis.

metadata_fix <-
  metadata_sc %>% 
  mutate(Genotype = gsub("-", "", Genotype),
         GenotypeTreatment = paste(Genotype, Treatment, sep = "_"))
  
metadata_fix

Table_S5_betadisper_results_pair <-
  BetadispExtr(asv_sc, metadata_fix, "GenotypeTreatment", "euclidean", 999)$disp_pair %>% 
  mutate(Adj_p_value = p.adjust(p_value, method = "BH"))

Table_S5_betadisper_results_pair

write.csv(Table_S5_betadisper_results_pair,"tables/Table_S5_betadisper_results_pair.csv")

# *********************************************---------------------------------
# ALPHA DIVERSITY --------------------------------------------------------------
# ADDING ALPHA METETRICS --------------------------------------------------------
AlphaMetrics <- function(physeq){
  require(vegan)
  require(tidyverse)
  
  sample_data(physeq)$ReadNo <- sample_sums(physeq)
  sample_data(physeq)$hill_0 <- as.data.frame(as.matrix(t(physeq@otu_table))) %>% renyi(scale=c(0), hill=TRUE)
  sample_data(physeq)$hill_1 <- as.data.frame(as.matrix(t(physeq@otu_table))) %>% renyi(scale=c(1), hill=TRUE)
  sample_data(physeq)$hill_2 <- as.data.frame(as.matrix(t(physeq@otu_table))) %>% renyi(scale=c(2), hill=TRUE)
  sample_data(physeq)$Richness <- specnumber(as.data.frame(otu_table(physeq)), MARGIN = 2)
  sample_data(physeq)$invSimpson <- diversity(as.data.frame(otu_table(physeq)), index="inv", MARGIN = 2)
  sample_data(physeq)$Shannon <- diversity(as.data.frame(otu_table(physeq)), index="shannon", MARGIN = 2)
  sample_data(physeq)$EH <- 1 - sample_data(physeq)$Shannon / log(sample_data(physeq)$Richness)
  
  # Normalize hill 1 to 0-1
  sample_data(physeq)$hill_1_norm <-  sample_data(physeq)$hill_1 / max(sample_data(physeq)$hill_1)
  sample_data(physeq)$hill_2_norm <-  sample_data(physeq)$hill_2 / max(sample_data(physeq)$hill_2)
  
  return(physeq)
}


physeq_sc <- AlphaMetrics(physeq_sc) 
physeq_sc %>% as_tibble()
physeq_sc@sam_data

# Need to reclassify each variable due to the conversion from a sam_data. The
# interesting thing is the hill number have weird class renij 
metadata_sc <- 
  as.data.frame(as.matrix(sample_data(physeq_sc))) %>% 
  mutate(hill_0 = as.numeric(hill_0),
         hill_1 = as.numeric(hill_1), 
         hill_1_norm = as.numeric(hill_1_norm),
         Replicate = as.factor(Replicate),
         Genotype = as.factor(Genotype),
         Treatment = as.factor(Treatment)) 

# TESTING USING WILCOX ---------------------------------------------------------

# Calculate number of multiple pairwise comparisons based on the number of levels
PairCalc <- function(n){
  n_comparison = n * (n - 1) / 2
  return(n_comparison)
}

PairCalc(6)

# Testing for significant differences ------------------------------------------
# Using wilcox.test, for BH you do not need to specify the number of comparisons.

CompSampl <- function(df, formula, adjust_method){
  require(multcompView)
  require(lazyeval)
  
  test_CC <-compare_means(formula,data = df,
                  method = "wilcox.test",
                  p.adjust.method = adjust_method)
  
  # Assign the data frame to the global environment
  assign("test_CC", test_CC, envir = .GlobalEnv)
  print(test_CC)
  
  test_CC <- as.data.frame(test_CC)[,c(2,3,5)] # to change form p to p.adj do 4 to 5
  test_CC2 <- data.frame(test_CC[,2], test_CC[,1], test_CC[,3])
  
  colnames(test_CC2) <- c("group1", "group2", "p.adj") # change p to p.adj
  rbind(test_CC, test_CC2) -> test_all
  print(test_all)
  
  # Assign the data frame to the global environment
  assign("test_all", test_all, envir = .GlobalEnv)
  
  as.dist(xtabs(test_all[, 3] ~ (test_all[, 2] + test_all[, 1])), diag = TRUE) -> dist_CC
  # print(dist_CC)
  data.frame(multcompLetters(dist_CC)['Letters']) -> res_CC
  print(res_CC)
  return(res_CC)
}

# Example usage (replace metadata_sc and formula as needed)
CompSampl(metadata_sc, formula(hill_0 ~ GenotypeTreatment), "BH")
CompSampl(metadata_sc, formula(hill_0 ~ Genotype), "BH")
CompSampl(metadata_sc, formula(hill_0 ~ Replicate), "BH")

# plotting Richness function ---------------------------------------------------

# I changed this function to work form a dataframe instead of a phyloseq object 
# becasue it is hard to mantain the order of the factors when you extract the 
# sam_data from a phyloseq object.

PlotRich <- function(dataframe, X_var, Y_var, label_df=NULL, labels_y=NULL){
  require(phyloseq)
  require(tidyverse)
  require(forcats)
  
  # Ensure Y_var is numeric
  dataframe <- dataframe %>%
    mutate(!!Y_var := as.numeric(!!sym(Y_var)))
  
  # print the dataframe levels 
  print(levels(dataframe[[X_var]]))
  
  # Calculate labels_y based on the maximum value of Y_var
  labels_y <- max(dataframe[[Y_var]]) + 0.1 * max(dataframe[[Y_var]])
  
  # plot
  rich_plot <-
    ggplot(dataframe, aes(x = get(X_var), y = !!sym(Y_var))) +
    geom_jitter(position=position_jitter(0.4), size=2, shape=16, #width = 0.25,
                aes(color=get(X_var))) +
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar",width = 0.8, color = "black") +
    
    # use geom_text() to plot the labels
    geom_text(data = label_df,
              aes(x = group, y = labels_y, label = Letters),
              size = 4, color = "black") +
    
    expand_limits(y = 0) +
    theme_bw() +
    theme(plot.title = element_markdown(size = 12, face = "bold", vjust = 0.5, hjust = 0.5),
          axis.title = element_markdown(),
          axis.text.x = element_markdown(size = 12, angle = 0, hjust = 1, vjust = 1, color = "black"),
          axis.text.y = element_markdown(size = 12, angle = 0, hjust = 0.5,  color = "black"),
          legend.position = "bottom", 
          legend.title=element_blank())
  
  return(rich_plot)
}


PlotRich(dataframe = sample_data(AddFit(physeq_sc, formula(hill_0 ~ Replicate), "2")[[3]]) %>% 
           as.matrix() %>% 
           as.data.frame() %>% 
           mutate(Genotype = fct_relevel(Genotype,"Col-0","cpr1","pad4","sid2","ndr1","erd4")),
         X_var =  "Genotype",
         Y_var =  ".resid",
         label_df = AddFit(physeq_sc, formula(hill_0 ~ Replicate), "2")[[2]] %>% 
           CompSampl(., formula(.resid ~ Genotype), "BH") %>% 
           rownames_to_column(var = "group") %>%
           select(group, Letters)
)




# We can see the effect of replicate clearly form the plots. We can model this
# and use the residuals for out Alpha diversity comparisons.

# ***** Figure S5 ***** --------------------------------------------------------
Fig_S5_Replicate_effect_on_Hill <-
  ggarrange(
    # hill_0
    PlotRich(metadata_sc, "Replicate", "hill_0",
             metadata_sc %>% 
               CompSampl(., formula(hill_0 ~ Replicate), "BH") %>% 
               rownames_to_column(var = "group") %>%
               select(group, Letters)) +
      theme(
        axis.text.x = element_markdown(size = 12, angle = 0, 
                                       hjust = 0.5, vjust = 0.5, color = "black")) +
      labs(title = "Effect of Replicate on hill_0", 
           x="Replicate"),
    
    # hill_1
    PlotRich(metadata_sc, "Replicate", "hill_1",
             metadata_sc %>% 
               CompSampl(., formula(hill_1 ~ Replicate), "BH") %>% 
               rownames_to_column(var = "group") %>%
               select(group, Letters)) +
      theme(
        axis.text.x = element_markdown(size = 12, angle = 0, 
                                       hjust = 0.5, vjust = 0.5, color = "black")) +
      labs(title = "Effect of Replicate on hill_1", 
           x="Replicate"),
    ncol = 2,
    nrow = 1,
    align = "hv", 
    labels = c("a", "b"))

Fig_S5_Replicate_effect_on_Hill


ggsave("figures/Fig_S5_Replicate_effect_on_Hill.pdf",
       plot = Fig_S5_Replicate_effect_on_Hill,
       device = "pdf")

# Replicates are different! This could be expected in a manipulative experiments.
# I will run a linear model to remove the effect of replicate form the data and re
# run the comparisons.

# Modeling alpha diversity -----------------------------------------------------

# Hill_0  ----------------------------------------------------------------------
# NOTE: sample pd43scw is a high leverage point, so I think it is worth removing 
# it for hill_0.

# First visualize the relationship
metadata_sc %>% 
  filter(SampleName != "pd43scw") %>% 
  ggplot(aes(x = hill_0, fill = Replicate, colour = Replicate)) + 
  theme_bw() +
  geom_histogram(binwidth = 5, alpha = 1, position = "identity") +
  labs(title = "Hill_0 distribution in 3 replicates")+
  theme(legend.position = "none", 
        legend.title = element_blank()) +
  facet_grid(~Replicate)


# start with a linear model, numeric vs numeric
contrasts(metadata_sc$Replicate)

fit_hill_0_m1 <-
  metadata_sc %>%
  filter(SampleName != "pd43scw") %>% 
  mutate(Replicate = relevel(Replicate, ref = "2")) %>% 
  lm(hill_0 ~ Replicate, data =.)

summary(fit_hill_0_m1)
anova(fit_hill_0_m1)
confint(fit_hill_0_m1)

# Generate diagnostic plots ----------------------------------------------------
# Using the lindia package 
hill_0_diagnostic <- gg_diagnose(fit_hill_0_m1, plot.all = FALSE)
plot_all(hill_0_diagnostic[c(1,2,3,4,5,6,7,8)])


# Define a better custom function
diagnostic_plots <- function(fit) {
  
  # Extract diagnostic data using broom + direct model functions
  diag_data <- augment(fit) %>%
    mutate(
      .hat = hatvalues(fit),           # Extract leverage directly
      .cooksd = cooks.distance(fit)    # Extract Cook's distance directly
    )
  
  # Histogram of residuals
  p1 <- ggplot(diag_data, aes(.resid), shape =1) +
    geom_histogram(aes(y = ..density..), bins = 30, color = "grey", fill = "grey") +
    stat_function(
      fun = dnorm, 
      args = list(mean = mean(diag_data$.resid), sd = sd(diag_data$.resid)),
      color = "red"
    ) +
    labs(x = "Residuals", y = "Density", title = "Histogram of Residuals") +
    theme_bw()
  
  # Normal Q-Q
  p2 <- ggplot(diag_data, aes(sample = .std.resid)) +
    stat_qq(size = 2, shape =1) +
    stat_qq_line(color = "red") +
    labs(title = "Normal Q-Q") +
    theme_bw()
  
  # Scale-Location
  p3 <- ggplot(diag_data, aes(.fitted, sqrt(abs(.std.resid)))) +
    geom_point(size = 2, shape = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(x = "Fitted Values", y = expression(sqrt("|Standardized Residuals|")), 
         title = "Scale-Location") +
    theme_bw()
  
  # Generate Cook's distance contour data
  leverage_seq <- seq(0, max(diag_data$.hat), length.out = 100)
  cooks_dist <- sqrt(4 / length(diag_data$.cooksd) * leverage_seq * (1 - leverage_seq))
  cooks_data <- data.frame(
    .hat = leverage_seq,
    cooks_dist = cooks_dist
  )
  
  # Cook's distance plot (stem plot style)
  cooks_threshold <- 4 / length(diag_data$.cooksd)  # Influential point threshold
  
  p4 <- ggplot(diag_data, aes(x = seq_along(.cooksd), y = .cooksd)) +
    
    # Vertical line segments (stems)
    geom_segment(aes(xend = seq_along(.cooksd), yend = 0), color = "black") +
    
    # Add points at the top of stems
    geom_point(size = 1.5, shape = 1, color = "black") +
    
    # Add threshold line
    geom_hline(yintercept = cooks_threshold, linetype = "dashed", color = "red") +
    
    # Label influential points
    geom_text(
      aes(label = ifelse(.cooksd > cooks_threshold, seq_along(.cooksd), "")),
      color = "red",
      hjust = -0.2, vjust = -0.5,
      size = 3
    ) +
    
    labs(x = "Leverage", y = "Standardized Residuals", 
         title = "Residuals vs Leverage") +
    theme_bw()
  
  # Combine plots using ggarrange
  combined_plot <- ggarrange(
    p1, p2, p3, p4, 
    ncol = 2, nrow = 2, 
    labels = c("a", "b", "c", "d")
  )
  
  # Return the combined plot
  return(combined_plot)
}


diagnostic_plots(fit_hill_0_m1)

# ***** Figure S6 ***** --------------------------------------------------------
grid.arrange(diagnostic_plots(fit_hill_0_m1),
    top = text_grob("lm(hill_0 ~ Replicate) diagnostic plots", size = 12, face = 2))

ggsave("figures/Fig_S6_diagnostics_hill_0.pdf",
       plot = grid.arrange(diagnostic_plots(fit_hill_0_m1),
                           top = text_grob("lm(hill_0 ~ Replicate) diagnostic plots", 
                                           size = 12, face = 2)),
       device = "pdf")


# Hill_1  ----------------------------------------------------------------------
# Since Rep 2 and Rep 3 are not differente, I will use Rep 1 as the reference group.

fit_hill_1_m1 <-
  metadata_sc %>%
  filter(SampleName != "pd43scw") %>% 
  mutate(Replicate = relevel(Replicate, ref = "2")) %>% 
  lm(hill_1 ~ Replicate, data =.)

summary(fit_hill_1_m1)
anova(fit_hill_1_m1)
confint(fit_hill_0_m1)


# ***** Figure S7 ***** --------------------------------------------------------
diagnostic_plots(fit_hill_1_m1)

grid.arrange(diagnostic_plots(fit_hill_1_m1),
             top = text_grob("lm(hill_1 ~ Replicate) diagnostic plots", size = 12, face = 2))

ggsave("figures/Fig_S7_diagnostics_hill_1.pdf",
       plot = grid.arrange(diagnostic_plots(fit_hill_1_m1),
                           top = text_grob("lm(hill_1 ~ Replicate) diagnostic plots", 
                                           size = 12, face = 2)),
       device = "pdf")


# ADD ALPHA DIVERSITY RESIDUALS ------------------------------------------------
# Using the broom package that has a function for this 

AddFit <- function(physeq, formula, ref_group){
  require(tidyverse)
  require(broom)
  
  dataframe <- 
    physeq@sam_data %>% 
    as.matrix() %>% 
    as.data.frame() %>%
    rownames_to_column("SampleID") %>% 
    filter(SampleID != "pd43scw") %>% 
    select(SampleID, Replicate, 
           Genotype, Treatment, GenotypeTreatment,
           hill_0, hill_1, hill_1_norm) %>% 
    mutate(hill_0 = as.numeric(hill_0),
           hill_1 = as.numeric(hill_1),
           hill_1_norm = as.numeric(hill_1_norm),
           Replicate = as.factor(Replicate)) %>% 
    mutate(Replicate = relevel(Replicate, ref = ref_group))
  
  model_fit <- lm(formula, data = dataframe)
  
  dataframe_new <-
    broom::augment(model_fit, dataframe) 
  
  physeq_fit <-
    phyloseq(
      otu_table(physeq, taxa_are_rows = TRUE),
      sample_data(dataframe_new %>% column_to_rownames("SampleID")),
      tax_table(physeq),
      refseq(physeq)) %>% 
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
    prune_samples(sample_sums(x=.) > 0, x =.)
  
  return(list(anova(model_fit), dataframe_new, physeq_fit))
}

AddFit(physeq_sc, formula(hill_0 ~ Replicate), "2")
AddFit(physeq_sc, formula(hill_1 ~ Replicate), "2")

# Plotting the residuals of Hill_0 and Hill_1 ----------------------------------

# ***** Figure S8 ***** --------------------------------------------------------
Fig_S8_Replicate_effect_on_Hill_resid <-
  ggarrange(
    # hill_0
    PlotRich(dataframe = sample_data(AddFit(physeq_sc, formula(hill_0 ~ Replicate), "2")[[3]]) %>% 
               as.matrix() %>% 
               as.data.frame() %>% 
               mutate(Genotype = fct_relevel(Genotype,"Col-0","cpr1","pad4","sid2","ndr1","erd4")),
             X_var =  "Replicate",
             Y_var =  ".resid",
             label_df = AddFit(physeq_sc, formula(hill_0 ~ Replicate), "2")[[2]] %>% 
               CompSampl(., formula(.resid ~ Replicate), "BH") %>% 
               rownames_to_column(var = "group") %>%
               select(group, Letters)) +
      theme(
        axis.text.x = element_markdown(size = 12, angle = 0, 
                                       hjust = 0.5, vjust = 0.5, color = "black")) +
      labs(title = "Effect of Replicate on Reduals<br>lm(hill_0 ~ Replicate)", 
           x="Replicate"),
    
    # hill_1
    PlotRich(dataframe = sample_data(AddFit(physeq_sc, formula(hill_1 ~ Replicate), "2")[[3]]) %>% 
               as.matrix() %>% 
               as.data.frame() %>% 
               mutate(Genotype = fct_relevel(Genotype,"Col-0","cpr1","pad4","sid2","ndr1","erd4")),
             X_var =  "Replicate",
             Y_var =  ".resid",
             label_df = AddFit(physeq_sc, formula(hill_1 ~ Replicate), "2")[[2]] %>% 
               CompSampl(., formula(.resid ~ Replicate), "BH") %>% 
               rownames_to_column(var = "group") %>%
               select(group, Letters)) +
      theme(
        axis.text.x = element_markdown(size = 12, angle = 0, 
                                       hjust = 0.5, vjust = 0.5, color = "black")) +
      labs(title = "Effect of Replicate on Reduals<br>lm(hill_1 ~ Replicate)", 
           x="Replicate"),
    ncol = 2,
    nrow = 1,
    align = "hv", 
    labels = c("a", "b"))

Fig_S8_Replicate_effect_on_Hill_resid

ggsave("figures/Fig_S8_Replicate_effect_on_Hill_resid.pdf",
       plot = Fig_S8_Replicate_effect_on_Hill_resid,
       device = "pdf")

# Plotting the alpha diversity -------------------------------------------------

# Define new palette for plotting Alpha
palette_cb6 <- c('#a50026','#a50026','#d73027','#d73027',
                 '#f46d43','#f46d43','#fdae61','#fdae61',
                 '#fee090','#fee090','#ffffbf','#ffffbf')
pie(rep(1, 12), col = palette_cb6)


Fig_hill_0_resid <- 
  ggarrange(
    PlotRich(dataframe = sample_data(AddFit(physeq_sc, formula(hill_0 ~ Replicate), "2")[[3]]) %>% 
           as.matrix() %>% 
           as.data.frame() %>% 
           mutate(Treatment = fct_relevel(Treatment, "w","d")),
         X_var =  "Treatment",
         Y_var =  ".resid",
         label_df = data.frame(group=c("d", "w"),
                               Letters=c("b", "a"))) +
      scale_color_manual(values=c('grey60', 'grey60')) +
      labs(title = "Treatment", x = NULL, y = "Residulas"),
PlotRich(dataframe = sample_data(AddFit(physeq_sc, formula(hill_0 ~ Replicate), "2")[[3]]) %>% 
           as.matrix() %>% 
           as.data.frame() %>% 
           mutate(Genotype = fct_relevel(Genotype,"Col-0","cpr1","pad4","sid2","ndr1","erd4")),
         X_var =  "Genotype",
         Y_var =  ".resid",
         label_df = AddFit(physeq_sc, formula(hill_0 ~ Replicate), "2")[[2]] %>% 
           CompSampl(., formula(.resid ~ Genotype), "BH") %>% 
           rownames_to_column(var = "group") %>%
           select(group, Letters)) +
       scale_color_manual(values = palette_cb6) +
       labs(title = "Genotype", x = NULL, y = "Residulas"),
PlotRich(dataframe = sample_data(AddFit(physeq_sc, formula(hill_0 ~ Replicate), "2")[[3]]) %>% 
           as.matrix() %>% 
           as.data.frame() %>% 
           mutate(GenotypeTreatment = fct_relevel(GenotypeTreatment,
                                                  "Col-0_w","Col-0_d","cpr1_w","cpr1_d","pad4_w","pad4_d",
                                                  "sid2_w","sid2_d","ndr1_w","ndr1_d","erd4_w","erd4_d")),
         X_var =  "GenotypeTreatment",
         Y_var =  ".resid",
         label_df = AddFit(physeq_sc, formula(hill_0 ~ Replicate), "2")[[2]] %>% 
           CompSampl(., formula(.resid ~ GenotypeTreatment), "BH") %>% 
           rownames_to_column(var = "group") %>%
           select(group, Letters)) +
       scale_color_manual(values = palette_cb6) +
       labs(title = "Genotype x Treatement", 
       x = NULL,
       y = NULL),
nrow = 1,
ncol = 3,
align = "hv",
widths = c(0.15, 0.27, 0.48),
labels = c("a", "b", "c"))

Fig_hill_0_resid


Fig_hill_1_resid <- 
  ggarrange(
    PlotRich(dataframe = sample_data(AddFit(physeq_sc, formula(hill_1 ~ Replicate), "2")[[3]]) %>% 
               as.matrix() %>% 
               as.data.frame() %>% 
               mutate(Treatment = fct_relevel(Treatment, "w","d")),
             X_var =  "Treatment",
             Y_var =  ".resid",
             label_df = data.frame(group=c("d", "w"),
                                   Letters=c("b", "a"))) +
      scale_color_manual(values=c('grey60', 'grey60')) +
      labs(title = "Treatment", x = NULL, y = "Residulas"),
    PlotRich(dataframe = sample_data(AddFit(physeq_sc, formula(hill_1 ~ Replicate), "2")[[3]]) %>% 
               as.matrix() %>% 
               as.data.frame() %>% 
               mutate(Genotype = fct_relevel(Genotype,"Col-0","cpr1","pad4","sid2","ndr1","erd4")),
             X_var =  "Genotype",
             Y_var =  ".resid",
             label_df = AddFit(physeq_sc, formula(hill_1 ~ Replicate), "2")[[2]] %>% 
               CompSampl(., formula(.resid ~ Genotype), "BH") %>% 
               rownames_to_column(var = "group") %>%
               select(group, Letters)) +
      scale_color_manual(values = palette_cb6) +
      labs(title = "Genotype", x = NULL, y = "Residulas"),
    PlotRich(dataframe = sample_data(AddFit(physeq_sc, formula(hill_1 ~ Replicate), "2")[[3]]) %>% 
               as.matrix() %>% 
               as.data.frame() %>% 
               mutate(GenotypeTreatment = fct_relevel(GenotypeTreatment,
                                                      "Col-0_w","Col-0_d","cpr1_w","cpr1_d","pad4_w","pad4_d",
                                                      "sid2_w","sid2_d","ndr1_w","ndr1_d","erd4_w","erd4_d")),
             X_var = "GenotypeTreatment",
             Y_var = ".resid",
             label_df = AddFit(physeq_sc, formula(hill_1 ~ Replicate), "2")[[2]] %>% 
               CompSampl(., formula(.resid ~ GenotypeTreatment), "BH") %>% 
               rownames_to_column(var = "group") %>%
               select(group, Letters)) +
      scale_color_manual(values = palette_cb6) +
      labs(title = "Genotype x Treatement", 
           x = NULL,
           y = NULL),
    nrow = 1,
    ncol = 3,
    align = "hv",
    widths = c(0.15, 0.27, 0.48),
    labels = c("a", "b", "c"))

Fig_hill_1_resid



Fig_2_alpha_diversity_residuals <-
ggarrange(
  grid.arrange(Fig_hill_0_resid,
               top = text_grob("lm(Hill 0 ~ Replicate)", size = 14, face = 2)),
  grid.arrange(Fig_hill_1_resid,
               top = text_grob("lm(Hill 1 ~ Replicate)", size = 14, face = 2)),
  ncol = 1,
  nrow = 2)

Fig_2_alpha_diversity_residuals

# ***** Figure 2 ***** --------------------------------------------------------

ggsave("figures/Fig_2_alpha_diversity_residuals.pdf", Fig_2_alpha_diversity_residuals,
       device = "pdf")


# Figure to put in the main figures

Fig_alpha_treatment_only <-
ggarrange(
  PlotRich(dataframe = sample_data(AddFit(physeq_sc, formula(hill_0 ~ Replicate), "2")[[3]]) %>% 
             as.matrix() %>% 
             as.data.frame() %>% 
             mutate(Treatment = fct_relevel(Treatment, "w","d")),
           X_var =  "Treatment",
           Y_var =  ".resid",
           label_df = data.frame(group=c("d", "w"),
                                 Letters=c("b", "a"))) +
    scale_color_manual(values=c('grey40', 'grey80')) +
    labs(title = "Hill 0", x = NULL, y = "Residulas"),
  PlotRich(dataframe = sample_data(AddFit(physeq_sc, formula(hill_1 ~ Replicate), "2")[[3]]) %>% 
           as.matrix() %>% 
           as.data.frame() %>% 
           mutate(Treatment = fct_relevel(Treatment, "w","d")),
         X_var =  "Treatment",
         Y_var =  ".resid",
         label_df = data.frame(group=c("d", "w"),
                               Letters=c("b", "a"))) +
  scale_color_manual(values=c('grey40', 'grey80')) +
  labs(title = "Hill 1", x = NULL, y = NULL, fill=NULL),
ncol=2,
align="hv",
legend="none")

Fig_alpha_treatment_only

ggsave("figures/Fig_alpha_treatment_only.pdf", Fig_alpha_treatment_only,
       device = "pdf")

# ***********************************************-------------------------------
# TOPIC MODEL ANALYSIS ---------------------------------------------------------

pacman::p_load(
  LinDA,
  tm,
  topicmodels,
  ldatuning,
  install=TRUE
)


# pooling to SynCom Isolate ----------------------------------------------------
# Format the Isolate column to make sure that tax_glom will work correctly.
physeq_iso <- physeq_sc

tax_table(physeq_iso) <- tax_table(apply(tax_table(physeq_sc), 2, function(x) gsub("_", "", x)))
tax_table(physeq_iso) <- tax_table(apply(tax_table(physeq_sc), 2, function(x) gsub("\\.", "", x)))

physeq_iso@tax_table %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  select(domain, Isolate) %>% 
  head()

physeq_iso <-
  phyloseq(
    otu_table(physeq_iso, taxa_are_rows = TRUE),
    sample_data(physeq_iso),
    tax_table(as.matrix(physeq_iso@tax_table %>% 
                          as.matrix() %>% 
                          as.data.frame() %>% 
                          select(domain, Isolate))),
    refseq(physeq_iso)) %>% 
  tax_glom(taxrank = "Isolate") %>% 
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.) 

physeq_iso
physeq_iso@tax_table
head(physeq_iso@sam_data)

# Extarct taxonomy lables ------------------------------------------------------
# for matching with isolate names 
write.csv(
  tax_table(physeq_iso) %>% as.data.frame(),
  "dataset/physeq_iso_tax_table.csv"
)

# importing correct isolate names ----------------------------------------------
isolate_names <- read.csv("dataset/physeq_iso_tax_table_with_names.csv", row.names=1) 
isolate_names

tax_table(physeq_iso) <- tax_table(as.matrix(isolate_names))
tax_table(physeq_iso)

# IDENTIFYING TOPICS -----------------------------------------------------------
# Generating a count matrix 
count_matrix <- as.data.frame(t(as.matrix(otu_table(physeq_iso))))
str(count_matrix)
head(count_matrix)

#docs <- tm::Corpus(VectorSource(count_matrix))
#dtm <- tm::DocumentTermMatrix(docs)

# Choosing the metric ----------------------------------------------------------

# CaoJuan2009: This metric tends to prefer models with fewer topics. It is based 
# on the similarity of topics, where lower values indicate better performance.
# It penalizes models with highly similar topics, thus favoring models where 
# topics are more distinct.

# Arun2010: This metric is based on the KL divergence between the document-topic 
# distribution and the topic-word distribution. It often favors a higher number 
# of topics compared to other metrics.

# Deveaud2014: This metric is designed to measure the coherence of topics. It 
# tends to balance between too few and too many topics, providing a middle ground.

# I need to modify the tax_table before I can merge by Isolate because there are
# unique values at other columns.

# ***********************************************-------------------------------
# VEM method for topics calling  -----------------------------------------------

# In LDA, a document is viewed as a distribution over topics, while a topic is a 
# distribution over words. To generate a document, LDA firstly samples a 
# document-specific multinomial distribution over topics from a Dirichlet distribution; 
# then repeatedly samples the words in the document from the corresponding 
# multinomial distribution.

# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02655-5

topics_VEM <- 
  FindTopicsNumber(
    round(count_matrix * 1000), # should expect no decimals so we can round
    topics = seq(from = 2, to = 36, by = 1),
    metrics = c("CaoJuan2009", "Arun2010"), #"Deveaud2014"
    method = "VEM",
    control = list(seed = 2025),
    mc.cores = 6,
    verbose = TRUE
  )

topics_VEM
str(topics_VEM)


topics_Gibbs <- 
  FindTopicsNumber(
    round(count_matrix * 1000), # should expect no decimals so we can round
    topics = seq(from = 2, to = 36, by = 1),
    metrics = c("CaoJuan2009", "Arun2010"), #"Deveaud2014"
    method = "Gibbs",
    control = list(seed = 2025),
    mc.cores = 6,
    verbose = TRUE
  )

topics_Gibbs
str(topics_Gibbs)

# NOTEs
# VEM (Variational Expectation-Maximization)
# Uses a variational approximation to compute the posterior distribution of the latent topics.
# It optimizes the variational lower bound using an iterative algorithm.
# Implemented through the topicmodels::LDA() function in R.
# Pros: Faster for small and medium-sized datasets.
# - Works well when the number of topics is not very large.
# - More deterministic (less variation in results between runs).

# Gibbs (Gibbs Sampling)
# Uses a Markov Chain Monte Carlo (MCMC) method to estimate the posterior distribution of the latent topics.
# Samples from the posterior distribution over many iterations until convergence.
# Implemented through the topicmodels::LDA() function in R with method = "Gibbs".
# Pros: More accurate for large and sparse datasets.
# - Tends to provide more stable results over multiple runs.
# - Less sensitive to poor starting conditions.

# Plot topic metrics -----------------------------------------------------------
topics_df <-
rbind(
  topics_VEM %>% mutate(method = "VEM"),
  topics_Gibbs %>% mutate(method = "Gibbs")
  )
  #pivot_longer(cols = c(-topics,-method), names_to = "group", values_to = "value")

topics_df


Fig_S9_topics_metrics <-
  ggarrange(
    topics_df %>% 
      ggplot(aes(x = topics, y = CaoJuan2009)) +
      geom_line() +
      geom_point() +
      facet_grid(~method, scales = "free") +
      theme_bw() +
      theme(plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
            strip.text = element_markdown(size = 10, face = "bold"),
            axis.text.x = element_markdown(angle = 0, size = 7, hjust = 0.5, vjust = 1.05),
            axis.text.y = element_markdown(angle = 0, size = 7, hjust = 1, vjust = 0.5),
            strip.background = element_blank(),
            legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
            legend.position = "right", legend.title.align = 0.5) +
      labs(title = "CaoJuan2009") + 
      scale_x_continuous(breaks = topics_df$topics),
    topics_df %>% 
      ggplot(aes(x = topics, y = Arun2010)) +
      geom_line() +
      geom_point()+
      facet_grid(~method, scales = "free") +
      theme_bw() +
      theme(plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
            strip.text = element_markdown(size = 10, face = "bold"),
            axis.text.x = element_markdown(angle = 0, size = 7, hjust = 0.5, vjust = 1.05),
            axis.text.y = element_markdown(angle = 0, size = 7, hjust = 1, vjust = 0.5),
            strip.background = element_blank(),
            legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
            legend.position = "right", legend.title.align = 0.5) +
      labs(title = "Arun2010") +
      scale_x_continuous(breaks = topics_df$topics),
    ncol = 1, 
    nrow = 2, 
    labels=c("a", "b"))

Fig_S9_topics_metrics

# ***** Figure S9 ***** --------------------------------------------------------
ggsave("figures/Fig_S9_topics_metrics.pdf", 
       grid.arrange(Fig_S9_topics_metrics, 
                   top = text_grob("LDA model score", size = 12, face = 2, hjust=0.5)),
       device = "pdf")

# PICKING TOPICS ---------------------------------------------------------------
# Selecting the number of topic to consider give the consensus

# I decide to go with VEM becasue Gibbs seems as many topics at the syncom members
# There are three strong signals at 8, 17 and 30 I will try all at first!

lda_VEM_k8 <-
  LDA(
    round(count_matrix * 1000),
    k = 8,
    method = "VEM",
    control = list(seed = 2025))

lda_VEM_k8
str(lda_VEM_k8)

lda_VEM_k17 <-
  LDA(
    round(count_matrix * 1000),
    k = 17,
    method = "VEM",
    control = list(seed = 2056)) 

lda_VEM_k17
str(lda_VEM_k17)

lda_VEM_k30 <-
  LDA(
    round(count_matrix * 1000),
    k = 30,
    method = "VEM",
    control = list(seed = 2027))

lda_VEM_k30
str(lda_VEM_k30)

# per-topic-per-word probabilities
beta_lda_VEM_k8 <- as.data.frame(tidy(lda_VEM_k8, matrix = "beta"))
beta_lda_VEM_k8

beta_lda_VEM_k17 <- as.data.frame(tidy(lda_VEM_k17, matrix = "beta"))
beta_lda_VEM_k17

beta_lda_VEM_k30 <- as.data.frame(tidy(lda_VEM_k30, matrix = "beta"))
beta_lda_VEM_k30


#per-document-per-topic probabilities
gamma_lda_VEM_k8 <- 
  as.data.frame(tidy(lda_VEM_k8, matrix = "gamma")) %>%
  arrange(document, topic)

gamma_lda_VEM_k17 <- 
  as.data.frame(tidy(lda_VEM_k17, matrix = "gamma")) %>%
  arrange(document, topic)

gamma_lda_VEM_k30 <- 
  as.data.frame(tidy(lda_VEM_k30, matrix = "gamma")) %>%
  arrange(document, topic)

# Multiply the per-document-per-topic probabilities by sample counts
getCountProb <- function(physeq, gamma_lda){
  
  library(tidyverse)
  
  lib_size <- 
    data.frame(sample_sums(physeq)) %>%
    dplyr::rename("readNo" = 1) %>%
    rownames_to_column(var = "document")
  
  tm_lda <- 
    left_join(lib_size, gamma_lda) %>%
    mutate(topic_count = readNo * gamma,
           topic_count = round(topic_count, 0)) %>%
    dplyr::select(-readNo, -gamma) %>%
    pivot_wider(names_from = topic, values_from = topic_count) %>%
    dplyr::rename_with(~ paste0("Topic_", .), -document) %>%
    column_to_rownames(var = "document") %>%
    t(.) %>%
    data.frame(.)
  
  ps_topic <- 
    phyloseq(
      sample_data(physeq),
      otu_table(tm_lda, taxa_are_rows = TRUE))
  
  return(ps_topic)
  
}

ps_topic_VEM_k8 <- getCountProb(physeq_iso, gamma_lda_VEM_k8)
ps_topic_VEM_k8

ps_topic_VEM_k17 <- getCountProb(physeq_iso, gamma_lda_VEM_k17)
ps_topic_VEM_k17

ps_topic_VEM_k30 <- getCountProb(physeq_iso, gamma_lda_VEM_k30)
ps_topic_VEM_k30

# ***********************************************-------------------------------
# FITTING LDA to TOPICS --------------------------------------------------------
fitLinda <- function(ps_topic){
  library(tidyverse)
  
  data.frame(sample_data(ps_topic)) %>% 
    mutate(Treatment = fct_relevel(Treatment, "w", "d")) %>% 
    head() %>% print()
  
  linda_vem <- linda(otu.tab = data.frame(otu_table(ps_topic)), 
                     meta = data.frame(sample_data(ps_topic)) %>% 
                       mutate(Treatment = fct_relevel(Treatment, "w", "d")), 
                     formula = "~ Treatment + (1|Replicate)", # Replicate as random effect 
                     imputation = FALSE,
                     alpha = 0.05, 
                     n.cores = 12)
  
  return(linda_vem)
}


linda_VEM_k8 <- fitLinda(ps_topic_VEM_k8)
linda_VEM_k8

linda_VEM_k17 <- fitLinda(ps_topic_VEM_k17)
linda_VEM_k17

linda_VEM_k30 <- fitLinda(ps_topic_VEM_k30)
linda_VEM_k30

# NOTE! The 'linda_sc$output$Treatmentd' means that we have positive log2Foldchange
# for the level 'd' in the metadata. If you want to change the direction just relevel
# the factor in the meta, as I did above.

# In out study three topics that have positive log2 fold changes and low FDR p-values 
# (< 0.05) for samples obtained from drought samples. We have 1 that have negative 
# log2 fold values and low FDR p-values highlighting enrichment of this topic in 
# the watered samples.


# Try plotting this in a different way
predicted_topics <-
  rbind(
    linda_VEM_k8$output$Treatmentd %>% 
      as.data.frame() %>% 
      dplyr::arrange(padj) %>%
      rownames_to_column(var = "Topic") %>% 
      mutate(reject = ifelse(padj < 0.05, "Yes", "No")) %>%
      separate(Topic, into = c("t", "t_n"), remove = FALSE, convert = TRUE) %>%
      mutate(Topic = gsub("_", " ", Topic)) %>%
      arrange(desc(t_n)) %>%
      mutate(Condition = rep("k=8", times = 8)), 
    linda_VEM_k17$output$Treatmentd %>% 
      as.data.frame() %>% 
      dplyr::arrange(padj) %>%
      rownames_to_column(var = "Topic") %>% 
      mutate(reject = ifelse(padj < 0.05, "Yes", "No")) %>%
      separate(Topic, into = c("t", "t_n"), remove = FALSE, convert = TRUE) %>%
      mutate(Topic = gsub("_", " ", Topic)) %>%
      arrange(desc(t_n)) %>% 
      mutate(Condition = rep("k=17", times = 17)),
    linda_VEM_k30$output$Treatmentd %>% 
      as.data.frame() %>% 
      dplyr::arrange(padj) %>%
      rownames_to_column(var = "Topic") %>% 
      mutate(reject = ifelse(padj < 0.05, "Yes", "No")) %>%
      separate(Topic, into = c("t", "t_n"), remove = FALSE, convert = TRUE) %>%
      mutate(Topic = gsub("_", " ", Topic)) %>%
      arrange(desc(t_n)) %>% 
      mutate(Condition = rep("k=30", times = 30))
  ) %>% 
  mutate(Condition = fct_relevel(Condition, "k=8", "k=17", "k=30"))

predicted_topics



Fig_S10_predicted_topics <-
predicted_topics %>%
  group_by(Condition) %>%
  complete(Topic = paste("Topic", 1:28), fill = list(Log2_Fold_Change = NA)) %>%
  mutate(Topic = fct_relevel(Topic, paste("Topic", 1:max(as.numeric(gsub("Topic ", "", Topic)))))) %>%
  mutate(Topic = fct_rev(Topic)) %>%  # Reverse the order of topics
  ungroup() %>%
  #ggplot(aes(x = Topic, y = log2FoldChange, fill = reject)) +
  mutate(color = if_else(pvalue <= 0.05, "#ff7f00", "grey"),
         color = if_else(padj <= 0.05, "#1f78b4", color)) %>% 
  ggplot(aes(x = Topic, y = log2FoldChange, fill = color)) +
  geom_col(width = 0.8) + # Set a reasonable width to control thickness of bars
  facet_wrap(~Condition, ncol = 3, scales = "free_y") + # Create the facets for each group
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Dashed line at y = 0
  coord_flip() + # Flip coordinates to have horizontal bars
  theme_bw() +
  annotate("text", x= Inf, y = Inf, label = "bold(drought)", 
           parse = TRUE, hjust=1.4, vjust=1.4, col="black", size = 3) +
  theme(plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_markdown(size = 10, hjust = 0.5),
        strip.text = element_markdown(size = 12, face = "bold"),
        axis.text.x = element_markdown(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
        axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
        legend.position = "none") +
  scale_fill_manual(values = c("#1f78b4", "#ff7f00","grey")) +
  labs(title = "Predicted topics", 
       x = "",
       y = "Log2 Fold-Change")

Fig_S10_predicted_topics


# ***** Figure S10 ***** --------------------------------------------------------
ggsave("figures/Fig_S10_predicted_topics.pdf", 
       Fig_S10_predicted_topics, device = "pdf")


# Logic on how to look into this. I think we first need evaluate the number of topics. 
# We saw that 10, 18 and 28 could be good number as showed by CaoJuan2009 and Arun2010
# metrics. Then for each of the three, we look into significantly enriched topics
# using the linda modeling approach. The next thing to do is to look into the
# probabilities that each single SC member have within the topics. The final step is 
# to fit the linda model for each topic and get a fold change for each specific SC member.
# Makes sense???

ExtractTopicProb <- function(physeq, beta_lda, beta_cutoff, k, topic_n){
  library(tidyverse)
  
  # extract one topic at a time
  topic_prob <-
    beta_lda %>%
    filter(beta >= beta_cutoff) %>%  #filtering beta to reduce noise, e.g. 0.02
    dplyr::filter(topic %in% c(topic_n)) %>% 
    arrange(desc(beta)) %>%
    arrange(desc(term)) %>%
    mutate(
      topic_isolate = paste(topic, term, sep = "_"),
      topic_isolate = factor(topic_isolate, levels = topic_isolate)) %>% 
    left_join(., physeq@tax_table %>% 
                as.matrix() %>% 
                as.data.frame() %>% 
                rownames_to_column("term"),
              by="term") %>% 
    mutate(
      topicN = rep(paste("topic=", topic_n, sep = ""), times = nrow(.)),
      k = paste("k=", k, sep = ""))
  
    #left_join(., isolate_names %>% 
    #            mutate(species = old_name), by = "species")
  
  return(topic_prob)
}

# test the function 
ExtractTopicProb(physeq_iso, beta_lda_VEM_k8, 0.02, 8, 2)

# Extract beta for drought -----------------------------------------------------

Alltopic_beta_df_drought <- rbind(
  # model with 8 topics
  ExtractTopicProb(physeq_iso, beta_lda_VEM_k8, 0.02, 8, 2),
  ExtractTopicProb(physeq_iso, beta_lda_VEM_k8, 0.02, 8, 6),
  
  # model with 17 topics
  ExtractTopicProb(physeq_iso, beta_lda_VEM_k17, 0.02, 17, 7),
  ExtractTopicProb(physeq_iso, beta_lda_VEM_k17, 0.02, 17, 8),
  ExtractTopicProb(physeq_iso, beta_lda_VEM_k17, 0.02, 17, 9),
  ExtractTopicProb(physeq_iso, beta_lda_VEM_k17, 0.02, 17, 15),
  
  # model with 30 topics
  ExtractTopicProb(physeq_iso, beta_lda_VEM_k30, 0.02, 30, 7),
  ExtractTopicProb(physeq_iso, beta_lda_VEM_k30, 0.02, 30, 14),
  ExtractTopicProb(physeq_iso, beta_lda_VEM_k30, 0.02, 30, 15)
  ) %>% 
  mutate(topicN = as.factor(topicN)) 


Alltopic_beta_df_drought
head(Alltopic_beta_df_drought)

# plotting function 
PlotBetaProb <- function(dataframe) {
  require(tidyverse)
  
  plot <- 
    dataframe %>%
    as.data.frame() %>%
    #mutate(k = fct_relevel(k, "k=8", "k=17", "k=30")) %>% 
    ggplot(aes(x = factor(topic), y = markdown_name, fill = beta)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(beta, 2)), size = 3, color = "black") +
    scale_fill_gradient(low = "#f7fbff", high = "#1f78b4") + 
    facet_grid(k ~ topic, scales = "free", space = "free") +
    theme_bw() +
    theme(
      plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_markdown(size = 10, hjust = 0.5),
      strip.text = element_markdown(
        size = 10, face = "bold"),
      axis.text = element_markdown(hjust = 0.5, vjust = 0.5, size = 10),
      strip.background = element_blank(),
      legend.key.height = unit(0.4, "cm"), 
      legend.key.width = unit(0.4, "cm"),
      panel.grid = element_blank()
    ) +
    labs(title = "Beta Probabilities",
         x = "Topic",
         y = NULL,
         fill = "Beta")
  
  return(plot)
}


Fig_S11_beta_drought <-
  Alltopic_beta_df_drought %>% 
  mutate(k = fct_relevel(k, "k=8", "k=17", "k=30")) %>% 
  PlotBetaProb()

Fig_S11_beta_drought <-
  Alltopic_beta_df_drought %>% 
  mutate(k = fct_relevel(k, "k=8", "k=17", "k=30")) %>% 
  PlotBetaProb() +
  scale_fill_gradient(low = "#ffffe5", high = "#d73027") +
  labs(subtitle = "Topics enriched in drought")

Fig_S11_beta_drought

# ***** Figure S11 ***** -------------------------------------------------------
ggsave("figures/Fig_S11_beta_drought.pdf", 
       Fig_S11_beta_drought, device = "pdf")


# Extract beta for well watered ------------------------------------------------

Alltopic_beta_df_watered <- rbind(
  # model with 8 topics
  ExtractTopicProb(physeq_iso, beta_lda_VEM_k8, 0.02, 8, 4),
  ExtractTopicProb(physeq_iso, beta_lda_VEM_k8, 0.02, 8, 5),
  
  # model with 17 topics
  ExtractTopicProb(physeq_iso, beta_lda_VEM_k17, 0.02, 17, 1),
  
  # model with 30 topics
  ExtractTopicProb(physeq_iso, beta_lda_VEM_k30, 0.02, 30, 1),
  ExtractTopicProb(physeq_iso, beta_lda_VEM_k30, 0.02, 30, 10),
  ExtractTopicProb(physeq_iso, beta_lda_VEM_k30, 0.02, 30, 18),
  ExtractTopicProb(physeq_iso, beta_lda_VEM_k30, 0.02, 30, 27)) %>% 
  mutate(topicN = as.factor(topicN)) 

Alltopic_beta_df_watered


Fig_S12_beta_water <-
  Alltopic_beta_df_watered %>% 
  mutate(k = fct_relevel(k, "k=8", "k=17", "k=30")) %>% 
  PlotBetaProb() +
  labs(subtitle = "Topics enriched in well-watered")

#Fig_S12_beta_water <-
#  PlotBetaProb(Alltopic_beta_df_watered) +
#  scale_fill_gradient(low = "#ffffe5", high = "#d73027") +
#  labs(subtitle = "Topics enriched in well watered")

Fig_S12_beta_water

# ***** Figure S12 ***** -------------------------------------------------------
ggsave("figures/Fig_S12_beta_water.pdf", 
       Fig_S12_beta_water, device = "pdf")


# FIT LINDA WITHIN EACH TOPIC --------------------------------------------------
# Fit linda from the LinDA package to specific topic to then visualize log fold 
# change for the SC member within a topic.

RunLinda <- function(physeq, beta_df, k_val, topic_n) {
  library(tidyverse)
  library(LinDA)
  
  # Clean up the beta table to make sure k is numeric
  beta_clean <- beta_df %>%
    mutate(
      k = as.numeric(gsub("k=", "", k)),
      topicN = as.character(topicN)
    ) %>%
    filter(topic == topic_n & k == k_val)
  
  if (nrow(beta_clean) == 0) {
    stop("No matching taxa found for the specified topic and k value")
  }
  
  # Get the taxa names to keep
  taxa_to_keep <- beta_clean$term
  
  if (length(taxa_to_keep) == 0) {
    stop("No matching taxa found after filtering")
  }
  
  # Subset the OTU table
  otu_tab <- as.data.frame(
    otu_table(
      prune_taxa(taxa_names(physeq) %in% taxa_to_keep, physeq)
    ) %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
      prune_samples(sample_sums(x = .) > 0, x = .)
  )
  
  # Subset the metadata table
  meta <- as.data.frame(as.matrix(
    sample_data(
      prune_taxa(taxa_names(physeq) %in% taxa_to_keep, physeq) %>%
        prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
        prune_samples(sample_sums(x = .) > 0, x = .)))) %>%
    mutate(Treatment = fct_relevel(Treatment, "w", "d"))
  
  # Check for missing metadata or low sample count
  if (nrow(meta) == 0) {
    stop("No matching samples found after filtering")
  }
  
  # Fit the `linda` model
  model <- linda(
    otu.tab = otu_tab,
    meta = meta,
    formula = "~ Treatment + (1 | Replicate)",  
    imputation = FALSE,
    p.adj.method = "BH",
    alpha = 0.05,
    n.cores = 6
  )
  
  return(model)

}


RunLinda(physeq_iso, Alltopic_beta_df_drought, 8, 2)


# Extract the Log2Fold change for each topic and each k 
ExtractLindaLog2Fold <- function(linda_df, physeq_obj){
  
  linda_res <-
    linda_df$output$Treatmentd %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "ZOTU") %>% 
    # need to match with the corrected isolate name
    left_join(., physeq_obj@tax_table %>% 
                as.matrix() %>% 
                as.data.frame() %>% 
                rownames_to_column("ZOTU"), 
              by="ZOTU") %>% 
    mutate(signif = if_else(pvalue<=0.05, TRUE, FALSE))
  
  return(linda_res)
}


ExtractLindaLog2Fold(RunLinda(physeq_iso, Alltopic_beta_df_drought, 8, 2), physeq_iso) %>% 
  mutate(k = "k=8", topicN = "topic=2")

Pertopic_Fold2Change_drought <-
  rbind(
    # model with 8 topics
    ExtractLindaLog2Fold(RunLinda(physeq_iso, Alltopic_beta_df_drought, 8, 2), physeq_iso) %>% 
      mutate(k = "k=8", topicN = "2"),
    ExtractLindaLog2Fold(RunLinda(physeq_iso, Alltopic_beta_df_drought, 8, 6), physeq_iso) %>% 
      mutate(k = "k=8", topicN = "6"),
    # model with 17 topics
    ExtractLindaLog2Fold(RunLinda(physeq_iso, Alltopic_beta_df_drought, 17, 7), physeq_iso) %>% 
      mutate(k = "k=17", topicN = "7"),
    ExtractLindaLog2Fold(RunLinda(physeq_iso, Alltopic_beta_df_drought, 17, 8), physeq_iso) %>% 
      mutate(k = "k=17", topicN = "8"),
    ExtractLindaLog2Fold(RunLinda(physeq_iso, Alltopic_beta_df_drought, 17, 9), physeq_iso) %>% 
      mutate(k = "k=17", topicN = "9"),
    ExtractLindaLog2Fold(RunLinda(physeq_iso, Alltopic_beta_df_drought, 17, 15), physeq_iso) %>% 
      mutate(k = "k=17", topicN = "15"),
    # model with 30 topics
    ExtractLindaLog2Fold(RunLinda(physeq_iso, Alltopic_beta_df_drought, 30, 7), physeq_iso) %>% 
      mutate(k = "k=30", topicN = "7"),
    ExtractLindaLog2Fold(RunLinda(physeq_iso, Alltopic_beta_df_drought, 30, 14), physeq_iso) %>% 
      mutate(k = "k=30", topicN = "14"),
    ExtractLindaLog2Fold(RunLinda(physeq_iso, Alltopic_beta_df_drought, 30, 15), physeq_iso) %>% 
      mutate(k = "k=30", topicN = "15")
    )

Pertopic_Fold2Change_drought


# Plot LogFoldChange using heatmap
head(Pertopic_Fold2Change_drought)
range(Pertopic_Fold2Change_drought$log2FoldChange)


Fig_S13_Fold2Change_drought <-
Pertopic_Fold2Change_drought %>%
  as.data.frame() %>%
  mutate(k = fct_relevel(k, "k=8", "k=17", "k=30")) %>% 
  mutate(topicN = as.factor(topicN),
         topicN = fct_relevel(topicN, "2", "6", "7", "8", "18", "9", "14", "15")) %>% 
  mutate(label_fmt = ifelse(signif,
    paste0("**", round(log2FoldChange, 2), "**"),  # bold if TRUE
    as.character(round(log2FoldChange, 2))         # plain if FALSE
  )) %>% 
  ggplot(aes(x = factor(topicN), y = markdown_name, fill = log2FoldChange)) +
  geom_tile(color = "white") +
  #geom_text(aes(label = round(log2FoldChange, 2)), size = 3, color = "black") +
  geom_richtext(
    aes(label = label_fmt),
    fill = NA, label.color = NA, color = "black",  # No background or box border
    size = 3) +
  scale_fill_gradient2(name = "LFC",
                       low = "#313695", 
                       high = "#d73027", 
                       mid = "white", 
                       na.value = "white",
                       limits = c(-1, 1),
                       breaks = c(-1, -0.5, 0, 0.5, 1)) +
  facet_grid(k ~ topicN, scales = "free", space = "free") +
  theme_bw() +
  theme(
    plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_markdown(size = 10, hjust = 0.5),
    strip.text = element_markdown(
      size = 10, face = "bold"),
    axis.text = element_markdown(hjust = 0.5, vjust = 0.5, size = 10),
    strip.background = element_blank(),
    legend.key.height = unit(0.4, "cm"), 
    legend.key.width = unit(0.4, "cm"),
    panel.grid = element_blank()
  ) +
  labs(title = "Differentially abundant syncom members", 
       subtitle = "Topics enriched in drought",
       x = "Topic",
       y = NULL,
       fill = "LFC")


Fig_S13_Fold2Change_drought

# ***** Figure S13 ***** -------------------------------------------------------
ggsave("figures/Fig_S13_Fold2Change_drought.pdf", 
       Fig_S13_Fold2Change_drought, device = "pdf")


Pertopic_Fold2Change_watered <-
  rbind(
    # model with 8 topics
    ExtractLindaLog2Fold(RunLinda(physeq_iso, Alltopic_beta_df_watered, 8, 4), physeq_iso) %>% 
      mutate(k = "k=8", topicN = "4"),
    ExtractLindaLog2Fold(RunLinda(physeq_iso, Alltopic_beta_df_watered, 8, 5), physeq_iso) %>% 
      mutate(k = "k=8", topicN = "5"),
    # model with 17 topics
    ExtractLindaLog2Fold(RunLinda(physeq_iso, Alltopic_beta_df_watered, 17, 1), physeq_iso) %>% 
      mutate(k = "k=17", topicN = "1"),
    # model with 30 topics
    ExtractLindaLog2Fold(RunLinda(physeq_iso, Alltopic_beta_df_watered, 30, 1), physeq_iso) %>% 
      mutate(k = "k=30", topicN = "1"),
    ExtractLindaLog2Fold(RunLinda(physeq_iso, Alltopic_beta_df_watered, 30, 10), physeq_iso) %>% 
      mutate(k = "k=30", topicN = "10"),
    ExtractLindaLog2Fold(RunLinda(physeq_iso, Alltopic_beta_df_watered, 30, 18), physeq_iso) %>% 
      mutate(k = "k=30", topicN = "18"),
    ExtractLindaLog2Fold(RunLinda(physeq_iso, Alltopic_beta_df_watered, 30, 27), physeq_iso) %>% 
      mutate(k = "k=30", topicN = "27")
  )

Pertopic_Fold2Change_watered

Fig_S14_Fold2Change_watered <-
Pertopic_Fold2Change_watered %>%
  as.data.frame() %>%
  mutate(k = fct_relevel(k, "k=8", "k=17", "k=30")) %>% 
  mutate(topicN = as.factor(topicN),
         topicN = fct_relevel(topicN, "1", "4", "5", "10", "18", "27")) %>% 
  mutate(log2FoldChange = log2FoldChange * -1) %>% # Invert the logfold sign
  mutate(label_fmt = ifelse(signif,
                            paste0("**", round(log2FoldChange, 2), "**"),  # bold if TRUE
                            as.character(round(log2FoldChange, 2))         # plain if FALSE
  )) %>% 
  ggplot(aes(x = factor(topicN), y = markdown_name, fill = log2FoldChange)) +
  geom_tile(color = "white") +
  #geom_text(aes(label = round(log2FoldChange, 2)), size = 3, color = "black") +
  geom_richtext(
    aes(label = label_fmt),
    fill = NA, label.color = NA, color = "black",  # No background or box border
    size = 3) +
  scale_fill_gradient2(name = "LFC",
                       low = "#d73027", 
                       high = "#313695", 
                       mid = "white", 
                       na.value = "white",
                       limits = c(-1, 1),
                       breaks = c(-1, -0.5, 0, 0.5, 1)) +
  facet_grid(k ~ topicN, scales = "free", space = "free") +
  theme_bw() +
  theme(
    plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_markdown(size = 10, hjust = 0.5),
    strip.text = element_markdown(
      size = 10, face = "bold"),
    axis.text = element_markdown(hjust = 0.5, vjust = 0.5, size = 10),
    strip.background = element_blank(),
    legend.key.height = unit(0.4, "cm"), 
    legend.key.width = unit(0.4, "cm"),
    panel.grid = element_blank()
  ) +
  labs(title = "Differentially abundant syncom members", 
       subtitle = "Topics enriched in well-watered",
       x = "Topic",
       y = NULL,
       fill = "LFC")

Fig_S14_Fold2Change_watered

# ***** Figure S14 ***** -------------------------------------------------------
ggsave("figures/Fig_S14_Fold2Change_watered.pdf", 
       Fig_S14_Fold2Change_watered, device = "pdf")


# Running Linda between genotypes ----------------------------------------------
physeq_iso@sam_data

physeq_iso_Col0_pad4 <-
physeq_iso %>% 
  subset_samples(Genotype %in% c("Col-0", "pad4") & Treatment %in% c("w")) %>% 
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.) 

physeq_iso_Col0_pad4@sam_data

physeq_iso_Col0_erd4 <-
physeq_iso %>% 
  subset_samples(Genotype %in% c("Col-0", "erd4") & Treatment %in% c("w")) %>% 
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.) 

physeq_iso_Col0_erd4@sam_data


# Generating a count matrix 
count_matrix_Col0_pad4 <- as.data.frame(t(as.matrix(otu_table(physeq_iso_Col0_pad4))))
str(count_matrix_Col0_pad4)

count_matrix_Col0_erd4 <- as.data.frame(t(as.matrix(otu_table(physeq_iso_Col0_erd4))))
str(count_matrix_Col0_erd4)

# I just go with VEM to detect topics becasue it seems it does a better job than Gibbs
# with this dataset.

# Col0-pad4
dim(count_matrix_Col0_pad4)

topics_VEM_Col0_pad4 <- 
  FindTopicsNumber(
    round(count_matrix_Col0_pad4 * 1000), # should expect no decimals so we can round
    topics = seq(from = 2, to = 36, by = 1),
    metrics = c("CaoJuan2009", "Arun2010"), #"Deveaud2014"
    method = "VEM",
    control = list(seed = 2026),
    mc.cores = 6,
    verbose = TRUE
  )

topics_VEM_Col0_pad4

# Col0-pad4
dim(count_matrix_Col0_erd4)

topics_VEM_Col0_erd4 <- 
  FindTopicsNumber(
    round(count_matrix_Col0_erd4 * 1000), # should expect no decimals so we can round
    topics = seq(from = 2, to = 33, by = 1),
    metrics = c("CaoJuan2009", "Arun2010"), #"Deveaud2014"
    method = "VEM",
    control = list(seed = 2026),
    mc.cores = 6,
    verbose = TRUE
  )

topics_VEM_Col0_erd4

# plotting
Fig_S15_topics_metrics_Col0_pad4_erd4 <-
  ggarrange(
    topics_VEM_Col0_pad4 %>% 
      pivot_longer(cols = c(-topics), names_to = "group", values_to = "value") %>% 
      ggplot(aes(x = topics, y = value)) +
      geom_line() +
      geom_point() +
      facet_wrap(~ group, scales = "free_y") +
      theme_bw() +
      theme(plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
            strip.text = element_markdown(size = 10, face = "bold"),
            axis.text.x = element_markdown(angle = 0, size = 7, hjust = 0.5, vjust = 1.05),
            axis.text.y = element_markdown(angle = 0, size = 7, hjust = 1, vjust = 0.5),
            strip.background = element_blank(),
            legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
            legend.position = "right", legend.title.align = 0.5) +
      labs(title = "Col-0 vs. pad4 in drought conditions") + 
      scale_x_continuous(breaks = topics_df$topics),
    topics_VEM_Col0_erd4 %>% 
      pivot_longer(cols = c(-topics), names_to = "group", values_to = "value") %>% 
      ggplot(aes(x = topics, y = value)) +
      geom_line() +
      geom_point()+
      facet_wrap(~ group, scales = "free_y") +
      theme_bw() +
      theme(plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
            strip.text = element_markdown(size = 10, face = "bold"),
            axis.text.x = element_markdown(angle = 0, size = 7, hjust = 0.5, vjust = 1.05),
            axis.text.y = element_markdown(angle = 0, size = 7, hjust = 1, vjust = 0.5),
            strip.background = element_blank(),
            legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
            legend.position = "right", legend.title.align = 0.5) +
      labs(title = "Col-0 vs. erd4 in drought conditions") +
      scale_x_continuous(breaks = topics_df$topics),
    ncol = 1, 
    nrow = 2, 
    labels=c("a", "b"))

Fig_S15_topics_metrics_Col0_pad4_erd4

# ***** Figure S15 ***** --------------------------------------------------------
ggsave("figures/Fig_S15_topics_metrics_Col0_pad4_erd4.pdf", 
       grid.arrange(Fig_S15_topics_metrics_Col0_pad4_erd4, 
                    top = text_grob("LDA model score", size = 12, face = 2, hjust=0.5)),
       device = "pdf")

# ***********************************************-------------------------------
# TOPIC ANALYSIS for Col-0 vs pad4 ---------------------------------------------
# For the Co0-pad4 there are three strong signals at 8, 10, 16 and 27

# per-topic-per-word probabilities, "beta" -------------------------------------
beta_lda_VEM_Col0_pad4_k8 <-
  as.data.frame(tidy(
    LDA(
    round(count_matrix_Col0_pad4 * 1000),
    k = 8,
    method = "VEM",
    control = list(seed = 2012)),
    matrix = "beta"))

beta_lda_VEM_Col0_pad4_k8

beta_lda_VEM_Col0_pad4_k10 <-
  as.data.frame(tidy(
    LDA(
      round(count_matrix_Col0_pad4 * 1000),
      k = 10,
      method = "VEM",
      control = list(seed = 2013)),
    matrix = "beta"))

beta_lda_VEM_Col0_pad4_k10

beta_lda_VEM_Col0_pad4_k16 <-
  as.data.frame(tidy(
    LDA(
      round(count_matrix_Col0_pad4 * 1000),
      k = 16,
      method = "VEM",
      control = list(seed = 2014)),
    matrix = "beta"))

beta_lda_VEM_Col0_pad4_k16

beta_lda_VEM_Col0_pad4_k26 <-
  as.data.frame(tidy(
    LDA(
      round(count_matrix_Col0_pad4 * 1000),
      k = 26,
      method = "VEM",
      control = list(seed = 2015)),
    matrix = "beta"))

beta_lda_VEM_Col0_pad4_k26

# per-document-per-topic probabilities, "gamma" --------------------------------
gamma_lda_VEM_Col0_pad4_k8 <-
  as.data.frame(tidy(
    LDA(
      round(count_matrix_Col0_pad4 * 1000),
      k = 8,
      method = "VEM",
      control = list(seed = 2012)),
    matrix = "gamma")) %>%
  arrange(document, topic)

gamma_lda_VEM_Col0_pad4_k8

gamma_lda_VEM_Col0_pad4_k10 <-
  as.data.frame(tidy(
    LDA(
      round(count_matrix_Col0_pad4 * 1000),
      k = 10,
      method = "VEM",
      control = list(seed = 2013)),
    matrix = "gamma")) %>%
  arrange(document, topic)

gamma_lda_VEM_Col0_pad4_k10

gamma_lda_VEM_Col0_pad4_k16 <-
  as.data.frame(tidy(
    LDA(
      round(count_matrix_Col0_pad4 * 1000),
      k = 16,
      method = "VEM",
      control = list(seed = 2014)),
    matrix = "gamma")) %>%
  arrange(document, topic)

gamma_lda_VEM_Col0_pad4_k16

gamma_lda_VEM_Col0_pad4_k26 <-
  as.data.frame(tidy(
    LDA(
      round(count_matrix_Col0_pad4 * 1000),
      k = 26,
      method = "VEM",
      control = list(seed = 2015)),
    matrix = "gamma")) %>%
  arrange(document, topic)

gamma_lda_VEM_Col0_pad4_k26

# Multiply the per-document-per-topic probabilities by sample counts
ps_topic_VEM_k8_Col0_pad4 <- getCountProb(physeq_iso_Col0_pad4, gamma_lda_VEM_Col0_pad4_k8)
ps_topic_VEM_k10_Col0_pad4 <- getCountProb(physeq_iso_Col0_pad4, gamma_lda_VEM_Col0_pad4_k10)
ps_topic_VEM_k16_Col0_pad4 <- getCountProb(physeq_iso_Col0_pad4, gamma_lda_VEM_Col0_pad4_k16)
ps_topic_VEM_k26_Col0_pad4 <- getCountProb(physeq_iso_Col0_pad4, gamma_lda_VEM_Col0_pad4_k26)

# FITTING LDA to TOPICS --------------------------------------------------------
fitLindaGenotype <- function(ps_topic, genotype_1, genotype_2){
  library(tidyverse)
  
  data.frame(sample_data(ps_topic)) %>% 
    mutate(Genotype = fct_relevel(Genotype, genotype_1, genotype_2)) %>% 
    head() %>% print()
  
  linda_vem <- linda(otu.tab = data.frame(otu_table(ps_topic)), 
                     meta = data.frame(sample_data(ps_topic)) %>% 
                       mutate(Genotype = fct_relevel(Genotype, genotype_1, genotype_2)), 
                     formula = "~ Genotype + (1|Replicate)", # Replicate as random effect 
                     imputation = FALSE,
                     alpha = 0.05, 
                     n.cores = 12)
  
  return(linda_vem)
}

# NOTE! The output$`GenotypeCol-0` means that we have positive log2Foldchange
# for the level 'Col-0' in the metadata. If you want to change the direction just relevel
# the factor in the function below, e.g. "Col-0" first then "pad4".

linda_VEM_k8_Col0_pad4 <- fitLindaGenotype(ps_topic_VEM_k8_Col0_pad4, "pad4", "Col-0")
linda_VEM_k10_Col0_pad4 <- fitLindaGenotype(ps_topic_VEM_k10_Col0_pad4, "pad4", "Col-0")
linda_VEM_k16_Col0_pad4 <- fitLindaGenotype(ps_topic_VEM_k16_Col0_pad4, "pad4", "Col-0")
linda_VEM_k26_Col0_pad4 <- fitLindaGenotype(ps_topic_VEM_k26_Col0_pad4, "pad4", "Col-0")

# Try plotting this in a different way
predicted_topics_Col0_pad4 <-
  rbind(
    linda_VEM_k8_Col0_pad4$output$`GenotypeCol-0`%>% 
      as.data.frame() %>% 
      dplyr::arrange(padj) %>%
      rownames_to_column(var = "Topic") %>% 
      mutate(reject = ifelse(padj < 0.05, "Yes", "No")) %>%
      separate(Topic, into = c("t", "t_n"), remove = FALSE, convert = TRUE) %>%
      mutate(Topic = gsub("_", " ", Topic)) %>%
      arrange(desc(t_n)) %>%
      mutate(Condition = rep("k=8", times = 8)), 
    linda_VEM_k10_Col0_pad4$output$`GenotypeCol-0` %>% 
      as.data.frame() %>% 
      dplyr::arrange(padj) %>%
      rownames_to_column(var = "Topic") %>% 
      mutate(reject = ifelse(padj < 0.05, "Yes", "No")) %>%
      separate(Topic, into = c("t", "t_n"), remove = FALSE, convert = TRUE) %>%
      mutate(Topic = gsub("_", " ", Topic)) %>%
      arrange(desc(t_n)) %>% 
      mutate(Condition = rep("k=10", times = 10)),
    linda_VEM_k16_Col0_pad4$output$`GenotypeCol-0` %>% 
      as.data.frame() %>% 
      dplyr::arrange(padj) %>%
      rownames_to_column(var = "Topic") %>% 
      mutate(reject = ifelse(padj < 0.05, "Yes", "No")) %>%
      separate(Topic, into = c("t", "t_n"), remove = FALSE, convert = TRUE) %>%
      mutate(Topic = gsub("_", " ", Topic)) %>%
      arrange(desc(t_n)) %>% 
      mutate(Condition = rep("k=16", times = 16)),
    linda_VEM_k26_Col0_pad4$output$`GenotypeCol-0` %>% 
      as.data.frame() %>% 
      dplyr::arrange(padj) %>%
      rownames_to_column(var = "Topic") %>% 
      mutate(reject = ifelse(padj < 0.05, "Yes", "No")) %>%
      separate(Topic, into = c("t", "t_n"), remove = FALSE, convert = TRUE) %>%
      mutate(Topic = gsub("_", " ", Topic)) %>%
      arrange(desc(t_n)) %>% 
      mutate(Condition = rep("k=26", times = 26))
  ) %>% 
  mutate(Condition = fct_relevel(Condition, "k=8", "k=10", "k=16", "k=26"))

predicted_topics_Col0_pad4


Fig_S16_predicted_topics_Col0_pad4 <-
  predicted_topics_Col0_pad4 %>%
  group_by(Condition) %>%
  complete(Topic = paste("Topic", 1:26), fill = list(Log2_Fold_Change = NA)) %>%
  mutate(Topic = fct_relevel(Topic, paste("Topic", 1:max(as.numeric(gsub("Topic ", "", Topic)))))) %>%
  mutate(Topic = fct_rev(Topic)) %>%  # Reverse the order of topics
  ungroup() %>%
  #ggplot(aes(x = Topic, y = log2FoldChange, fill = reject)) +
  mutate(color = if_else(pvalue <= 0.05, "#ff7f00", "grey"),
         color = if_else(padj <= 0.05, "#1f78b4", color)) %>% 
  ggplot(aes(x = Topic, y = log2FoldChange, fill = color)) +
  geom_col(width = 0.8) + # Set a reasonable width to control thickness of bars
  facet_wrap(~Condition, ncol = 4, scales = "free_y") + # Create the facets for each group
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Dashed line at y = 0
  coord_flip() + # Flip coordinates to have horizontal bars
  theme_bw() +
  annotate("text",  x= Inf, y = Inf, label = "bold(`Col-0`)", 
           parse = TRUE, hjust=1.4, vjust=1.4, col="black", , size = 3) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
        axis.text.y = element_text(angle = 0, size = 8, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
        legend.position = "none") +
  scale_fill_manual(values = c("#1f78b4", "#ff7f00","grey")) +
  labs(title = "Predicted topics between Col-0 and pad4 genotypes", 
       x = "",
       y = "Log2 Fold-Change")

Fig_S16_predicted_topics_Col0_pad4

# ***** Figure S16 ***** --------------------------------------------------------
ggsave("figures/Fig_S16_predicted_topics_Col0_pad4.pdf", 
       Fig_S16_predicted_topics_Col0_pad4, device = "pdf")


# Extract beta for enriched in Col-0 -------------------------------------------
ExtractTopicProb(physeq_iso_Col0_pad4, beta_lda_VEM_Col0_pad4_k10, 0.02, 10, 2)

Alltopic_beta_df_genotype_Col0_pad4 <- 
  rbind(
  # model with 10 topics
  ExtractTopicProb(physeq_iso_Col0_pad4, beta_lda_VEM_Col0_pad4_k10, 0.02, 10, 2),
  
  # model with 30 topics
  ExtractTopicProb(physeq_iso_Col0_pad4, beta_lda_VEM_Col0_pad4_k26, 0.02, 26, 7),
  ExtractTopicProb(physeq_iso_Col0_pad4, beta_lda_VEM_Col0_pad4_k26, 0.02, 26, 9),
  ExtractTopicProb(physeq_iso_Col0_pad4, beta_lda_VEM_Col0_pad4_k26, 0.02, 26, 14),
  ExtractTopicProb(physeq_iso_Col0_pad4, beta_lda_VEM_Col0_pad4_k26, 0.02, 26, 18),
  ExtractTopicProb(physeq_iso_Col0_pad4, beta_lda_VEM_Col0_pad4_k26, 0.02, 26, 20),
  ExtractTopicProb(physeq_iso_Col0_pad4, beta_lda_VEM_Col0_pad4_k26, 0.02, 26, 22),
  ExtractTopicProb(physeq_iso_Col0_pad4, beta_lda_VEM_Col0_pad4_k26, 0.02, 26, 24)
  ) %>% 
  mutate(topicN = as.factor(topicN)) 


Alltopic_beta_df_genotype_Col0_pad4
head(Alltopic_beta_df_genotype_Col0_pad4)

# plotting 

Fig_S17_beta_genotype_Col0_pad4 <-
  PlotBetaProb(Alltopic_beta_df_genotype_Col0_pad4)

Fig_S17_beta_genotype_Col0_pad <-
  Alltopic_beta_df_genotype_Col0_pad4 %>% 
  mutate(k = fct_relevel(k, "k=10", "k=26")) %>% 
  PlotBetaProb() +
  scale_fill_gradient(low = "#ffffe5", high = "#d73027") +
  labs(subtitle = "Topics enriched in Col-0 compared to pad4 genotypes under drought conditions")

Fig_S17_beta_genotype_Col0_pad

# ***** Figure S17 ***** -------------------------------------------------------
ggsave("figures/Fig_S17_beta_genotype_Col0_pad4.pdf", 
       Fig_S17_beta_genotype_Col0_pad, device = "pdf")


# Fit linda between Col-0 and pad4 ---------------------------------------------

RunLindaGenotype <- function(physeq, beta_df, k_val, topic_n, 
                             genotype_1, genotype_2) {
  library(tidyverse)
  library(LinDA)
  
  # Clean up the beta table to make sure k is numeric
  beta_clean <- beta_df %>%
    mutate(
      k = as.numeric(gsub("k=", "", k)),
      topicN = as.character(topicN)
    ) %>%
    filter(topic == topic_n & k == k_val)
  
  if (nrow(beta_clean) == 0) {
    stop("No matching taxa found for the specified topic and k value")
  }
  
  # Get the taxa names to keep
  taxa_to_keep <- beta_clean$term
  
  if (length(taxa_to_keep) == 0) {
    stop("No matching taxa found after filtering")
  }
  
  # Subset the OTU table
  otu_tab <- as.data.frame(
    otu_table(
      prune_taxa(taxa_names(physeq) %in% taxa_to_keep, physeq)
    ) %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
      prune_samples(sample_sums(x = .) > 0, x = .)
  )
  
  # Subset the metadata table
  meta <- as.data.frame(as.matrix(
    sample_data(
      prune_taxa(taxa_names(physeq) %in% taxa_to_keep, physeq) %>%
        prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
        prune_samples(sample_sums(x = .) > 0, x = .)))) %>%
    mutate(Genotype = fct_relevel(Genotype, genotype_1, genotype_2))
  
  # Check for missing metadata or low sample count
  if (nrow(meta) == 0) {
    stop("No matching samples found after filtering")
  }
  
  # Fit the `linda` model
  model <- linda(
    otu.tab = otu_tab,
    meta = meta,
    formula = "~ Genotype + (1 | Replicate)",  
    imputation = FALSE,
    p.adj.method = "BH",
    alpha = 0.05,
    n.cores = 6
  )
  
  return(model)
  
}


RunLindaGenotype(physeq_iso_Col0_pad4, 
                 Alltopic_beta_df_genotype_Col0_pad4, 10, 2,"pad4", "Col-0")


ExtractLindaLog2FoldGenotype <- function(linda_df, physeq_obj){
  
  linda_res <-
    linda_df$output$`GenotypeCol-0` %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "ZOTU") %>% 
    # need to match with the corrected isolate name
    left_join(., physeq_obj@tax_table %>% 
                as.matrix() %>% 
                as.data.frame() %>% 
                rownames_to_column("ZOTU"), 
              by="ZOTU") %>% 
    mutate(signif = if_else(pvalue <= 0.05, TRUE, FALSE))
  
  return(linda_res)
}



ExtractLindaLog2FoldGenotype(
  RunLindaGenotype(physeq_iso_Col0_pad4, 
                   Alltopic_beta_df_genotype_Col0_pad4,10,2,"pad4","Col-0"),
  physeq_iso_Col0_pad4) 


Pertopic_Fold2Change_Col0_pad4 <-
  rbind(
    # k=10
    ExtractLindaLog2FoldGenotype(
    RunLindaGenotype(physeq_iso_Col0_pad4, 
                     Alltopic_beta_df_genotype_Col0_pad4,10,2,"pad4","Col-0"),
    physeq_iso_Col0_pad4)%>% 
      mutate(k = "k=10", topicN = "2"),
    # k=26
    ExtractLindaLog2FoldGenotype(
      RunLindaGenotype(physeq_iso_Col0_pad4, 
                       Alltopic_beta_df_genotype_Col0_pad4,26,7,"pad4","Col-0"),
      physeq_iso_Col0_pad4)%>% 
      mutate(k = "k=26", topicN = "7"),
    ExtractLindaLog2FoldGenotype(
      RunLindaGenotype(physeq_iso_Col0_pad4, 
                       Alltopic_beta_df_genotype_Col0_pad4,26,9,"pad4","Col-0"),
      physeq_iso_Col0_pad4)%>% 
      mutate(k = "k=26", topicN = "9"),
    ExtractLindaLog2FoldGenotype(
      RunLindaGenotype(physeq_iso_Col0_pad4, 
                       Alltopic_beta_df_genotype_Col0_pad4,26,14,"pad4","Col-0"),
      physeq_iso_Col0_pad4) %>% 
      mutate(k = "k=26", topicN = "14"),
    ExtractLindaLog2FoldGenotype(
      RunLindaGenotype(physeq_iso_Col0_pad4, 
                       Alltopic_beta_df_genotype_Col0_pad4,26,18,"pad4","Col-0"),
      physeq_iso_Col0_pad4) %>% 
      mutate(k = "k=26", topicN = "18"),
    ExtractLindaLog2FoldGenotype(
      RunLindaGenotype(physeq_iso_Col0_pad4, 
                       Alltopic_beta_df_genotype_Col0_pad4,26,20,"pad4","Col-0"),
      physeq_iso_Col0_pad4) %>% 
      mutate(k = "k=26", topicN = "20"),
    ExtractLindaLog2FoldGenotype(
      RunLindaGenotype(physeq_iso_Col0_pad4, 
                       Alltopic_beta_df_genotype_Col0_pad4,26,22,"pad4","Col-0"),
      physeq_iso_Col0_pad4) %>% 
      mutate(k = "k=26", topicN = "22"),
    ExtractLindaLog2FoldGenotype(
      RunLindaGenotype(physeq_iso_Col0_pad4, 
                       Alltopic_beta_df_genotype_Col0_pad4,26,24,"pad4","Col-0"),
      physeq_iso_Col0_pad4) %>% 
      mutate(k = "k=26", topicN = "24")
    )

Pertopic_Fold2Change_Col0_pad4

Fig_S18_Fold2Change_Col0_pad4 <-
  Pertopic_Fold2Change_Col0_pad4 %>%
  mutate(topicN_char = as.character(topicN)) %>%
  mutate(topicN = fct_relevel(topicN_char, "2", "7", "9", "14", "18", "20", "22", "24")) %>%  # Corrected line
  mutate(k = fct_relevel(k, "k=10", "k=26")) %>% 
  mutate(label_fmt = ifelse(signif,
                            paste0("**", round(log2FoldChange, 2), "**"),  # bold if TRUE
                            as.character(round(log2FoldChange, 2))         # plain if FALSE
  )) %>% 
  ggplot(aes(x = factor(topicN), y = markdown_name, fill = log2FoldChange)) +
  geom_tile(color = "white") +
  #geom_text(aes(label = round(log2FoldChange, 2)), size = 3, color = "black") +
  geom_richtext(
    aes(label = label_fmt),
    fill = NA, label.color = NA, color = "black",  # No background or box border
    size = 3) +
  scale_fill_gradient2(name = "LFC",
                       low = "#313695", 
                       high = "#d73027", 
                       mid = "white", 
                       na.value = "white",
                       limits = c(-1, 1),
                       breaks = c(-1, -0.5, 0, 0.5, 1)) +
  facet_grid(k ~ topicN, scales = "free", space = "free") +
  theme_bw() +
  theme(
    plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_markdown(size = 10, hjust = 0.5),
    strip.text = element_markdown(
      size = 10, face = "bold"),
    axis.text = element_markdown(hjust = 0.5, vjust = 0.5, size = 10),
    strip.background = element_blank(),
    legend.key.height = unit(0.4, "cm"), 
    legend.key.width = unit(0.4, "cm"),
    panel.grid = element_blank()
  ) +
  labs(title = "Differentially abundant syncom members", 
       subtitle = "Topics enriched in Col-0 compared to pad4 genotypes under drought conditions",
       x = "Topic",
       y = NULL,
       fill = "LFC")


Fig_S18_Fold2Change_Col0_pad4

# ***** Figure S18 ***** -------------------------------------------------------
ggsave("figures/Fig_S18_Fold2Change_Col0_pad4.pdf", 
       Fig_S18_Fold2Change_Col0_pad4, device = "pdf")

# ***********************************************-------------------------------
# TOPIC ANALYSIS for Col-0 vs erd4 ---------------------------------------------
# For the Col0-erd4 there are three strong signals at 10, 13 and 30

# per-topic-per-word probabilities, "beta" -------------------------------------
beta_lda_VEM_Col0_erd4_k10 <-
  as.data.frame(tidy(
    LDA(
      round(count_matrix_Col0_erd4 * 1000),
      k = 10,
      method = "VEM",
      control = list(seed = 2012)),
    matrix = "beta"))

beta_lda_VEM_Col0_erd4_k10

beta_lda_VEM_Col0_erd4_k13 <-
  as.data.frame(tidy(
    LDA(
      round(count_matrix_Col0_erd4 * 1000),
      k = 13,
      method = "VEM",
      control = list(seed = 2013)),
    matrix = "beta"))

beta_lda_VEM_Col0_erd4_k13

beta_lda_VEM_Col0_erd4_k30 <-
  as.data.frame(tidy(
    LDA(
      round(count_matrix_Col0_erd4 * 1000),
      k = 30,
      method = "VEM",
      control = list(seed = 2015)),
    matrix = "beta"))

beta_lda_VEM_Col0_erd4_k30

# per-document-per-topic probabilities, "gamma" --------------------------------
gamma_lda_VEM_Col0_erd4_k10 <-
  as.data.frame(tidy(
    LDA(
      round(count_matrix_Col0_erd4 * 1000),
      k = 10,
      method = "VEM",
      control = list(seed = 2023)),
    matrix = "gamma")) %>%
  arrange(document, topic)

gamma_lda_VEM_Col0_erd4_k10

gamma_lda_VEM_Col0_erd4_k13 <-
  as.data.frame(tidy(
    LDA(
      round(count_matrix_Col0_erd4 * 1000),
      k = 13,
      method = "VEM",
      control = list(seed = 2024)),
    matrix = "gamma")) %>%
  arrange(document, topic)

gamma_lda_VEM_Col0_erd4_k13

gamma_lda_VEM_Col0_erd4_k30 <-
  as.data.frame(tidy(
    LDA(
      round(count_matrix_Col0_erd4 * 1000),
      k = 30,
      method = "VEM",
      control = list(seed = 2025)),
    matrix = "gamma")) %>%
  arrange(document, topic)

gamma_lda_VEM_Col0_erd4_k30

# Multiply the per-document-per-topic probabilities by sample counts
ps_topic_VEM_k10_Col0_erd4 <- getCountProb(physeq_iso_Col0_erd4, gamma_lda_VEM_Col0_erd4_k10)
ps_topic_VEM_k13_Col0_erd4 <- getCountProb(physeq_iso_Col0_erd4, gamma_lda_VEM_Col0_erd4_k13)
ps_topic_VEM_k30_Col0_erd4 <- getCountProb(physeq_iso_Col0_erd4, gamma_lda_VEM_Col0_erd4_k30)

# FITTING LDA to TOPICS --------------------------------------------------------

# NOTE! The output$`GenotypeCol-0` means that we have positive log2Foldchange
# for the level 'Col-0' in the metadata. If you want to change the direction just relevel
# the factor in the function below, e.g. "Col-0" first then "erd4".

linda_VEM_k10_Col0_erd4 <- fitLindaGenotype(ps_topic_VEM_k10_Col0_erd4, "erd4", "Col-0")
linda_VEM_k13_Col0_erd4 <- fitLindaGenotype(ps_topic_VEM_k13_Col0_erd4, "erd4", "Col-0")
linda_VEM_k30_Col0_erd4 <- fitLindaGenotype(ps_topic_VEM_k30_Col0_erd4, "erd4", "Col-0")

# Try plotting this in a different way
predicted_topics_Col0_erd4 <-
  rbind(
    linda_VEM_k10_Col0_erd4$output$`GenotypeCol-0` %>% 
      as.data.frame() %>% 
      dplyr::arrange(padj) %>%
      rownames_to_column(var = "Topic") %>% 
      mutate(reject = ifelse(padj < 0.05, "Yes", "No")) %>%
      separate(Topic, into = c("t", "t_n"), remove = FALSE, convert = TRUE) %>%
      mutate(Topic = gsub("_", " ", Topic)) %>%
      arrange(desc(t_n)) %>% 
      mutate(Condition = rep("k=10", times = 10)),
    linda_VEM_k13_Col0_erd4$output$`GenotypeCol-0` %>% 
      as.data.frame() %>% 
      dplyr::arrange(padj) %>%
      rownames_to_column(var = "Topic") %>% 
      mutate(reject = ifelse(padj < 0.05, "Yes", "No")) %>%
      separate(Topic, into = c("t", "t_n"), remove = FALSE, convert = TRUE) %>%
      mutate(Topic = gsub("_", " ", Topic)) %>%
      arrange(desc(t_n)) %>% 
      mutate(Condition = rep("k=13", times = 13)),
    linda_VEM_k30_Col0_erd4$output$`GenotypeCol-0` %>% 
      as.data.frame() %>% 
      dplyr::arrange(padj) %>%
      rownames_to_column(var = "Topic") %>% 
      mutate(reject = ifelse(padj < 0.05, "Yes", "No")) %>%
      separate(Topic, into = c("t", "t_n"), remove = FALSE, convert = TRUE) %>%
      mutate(Topic = gsub("_", " ", Topic)) %>%
      arrange(desc(t_n)) %>% 
      mutate(Condition = rep("k=30", times = 30))
  ) %>% 
  mutate(Condition = fct_relevel(Condition, "k=10", "k=13", "k=30"))

predicted_topics_Col0_erd4


Fig_S19_predicted_topics_Col0_erd4 <-
  predicted_topics_Col0_erd4 %>%
  group_by(Condition) %>%
  complete(Topic = paste("Topic", 1:30), fill = list(Log2_Fold_Change = NA)) %>%
  mutate(Topic = fct_relevel(Topic, paste("Topic", 1:max(as.numeric(gsub("Topic ", "", Topic)))))) %>%
  mutate(Topic = fct_rev(Topic)) %>%  # Reverse the order of topics
  ungroup() %>%
  mutate(color = if_else(pvalue <= 0.05, "#ff7f00", "grey"),
         color = if_else(padj <= 0.05, "#1f78b4", color)) %>% 
  ggplot(aes(x = Topic, y = log2FoldChange, fill = color)) +
  geom_col(width = 0.8) + # Set a reasonable width to control thickness of bars
  facet_wrap(~Condition, ncol = 4, scales = "free_y") + # Create the facets for each group
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Dashed line at y = 0
  coord_flip() + # Flip coordinates to have horizontal bars
  theme_bw() +
  annotate("text",  x= Inf, y = Inf, label = "bold(`Col-0`)", 
           parse = TRUE, hjust=1.4, vjust=1.4, col="black", , size = 3) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
        axis.text.y = element_text(angle = 0, size = 8, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
        legend.position = "none") +
  scale_fill_manual(values = c("#1f78b4", "#ff7f00","grey")) +
  labs(title = "Predicted topics between Col-0 and erd4 genotypes", 
       x = "",
       y = "Log2 Fold-Change")

Fig_S19_predicted_topics_Col0_erd4

# ***** Figure S19 ***** --------------------------------------------------------
ggsave("figures/Fig_S19_predicted_topics_Col0_erd4.pdf", 
       Fig_S19_predicted_topics_Col0_erd4, device = "pdf")


# Extract beta for enriched in Col-0 -------------------------------------------
topic_beta_k10_genotype_Col0 <-
  PlotBetaProb(ExtractTopicProb(physeq_iso_Col0_erd4, beta_lda_VEM_Col0_erd4_k10, 0.02, 10, 3))

Pertopic_Fold2Change_k10_Col0 <-
ExtractLindaLog2FoldGenotype(
  RunLindaGenotype(physeq_iso_Col0_erd4, 
                   ExtractTopicProb(physeq_iso_Col0_erd4, beta_lda_VEM_Col0_erd4_k10, 0.02, 10, 3), 
                   10, 3, "erd4","Col-0"),
  physeq_iso_Col0_erd4) %>% 
  mutate(k = "k=10", topicN = "3")

Pertopic_Fold2Change_k10_Col0

Fig_S20_Fold2Change_k10_Col0 <-
ggarrange(
  topic_beta_k10_genotype_Col0 +
    labs(subtitle = "Topics enriched<br>in Col-0 compared to erd4 genotypes<br>under drought conditions"),
  Pertopic_Fold2Change_k10_Col0 %>% 
  mutate(label_fmt = ifelse(signif,
                            paste0("**", round(log2FoldChange, 2), "**"),  # bold if TRUE
                            as.character(round(log2FoldChange, 2))         # plain if FALSE
  )) %>% 
  ggplot(aes(x = factor(topicN), y = markdown_name, fill = log2FoldChange)) +
  geom_tile(color = "white") +
  #geom_text(aes(label = round(log2FoldChange, 2)), size = 3, color = "black") +
  geom_richtext(
    aes(label = label_fmt),
    fill = NA, label.color = NA, color = "black",  # No background or box border
    size = 3) +
  scale_fill_gradient2(name = "LFC",
                       low = "#313695", 
                       high = "#d73027", 
                       mid = "white", 
                       na.value = "white",
                       limits = c(-1, 1),
                       breaks = c(-1, -0.5, 0, 0.5, 1)) +
  facet_grid(k ~ topicN, scales = "free", space = "free") +
  theme_bw() +
  theme(
    plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_markdown(size = 10, hjust = 0.5),
    strip.text = element_markdown(
      size = 10, face = "bold"),
    axis.text = element_markdown(hjust = 0.5, vjust = 0.5, size = 10),
    strip.background = element_blank(),
    legend.key.height = unit(0.4, "cm"), 
    legend.key.width = unit(0.4, "cm"),
    panel.grid = element_blank()
  ) +
  labs(title = "Differential abundance", 
       subtitle = "Topics enriched<br>in Col-0 compared to erd4 genotypes<br>under drought conditions",
       x = "Topic",
       y = NULL,
       fill = "LFC"),
  ncol=2,
  nrow=1, 
  align="hv")

Fig_S20_Fold2Change_k10_Col0

# ***** Figure S20 ***** -------------------------------------------------------
ggsave("figures/Fig_S20_Fold2Change_k10_Col0.pdf", 
       Fig_S20_Fold2Change_k10_Col0, device = "pdf")


# Extract beta for enriched in erd4 --------------------------------------------
Alltopic_beta_df_genotype_Col0_erd4 <- 
  rbind(
    # model with 10 topics
    ExtractTopicProb(physeq_iso_Col0_erd4, beta_lda_VEM_Col0_erd4_k13, 0.02, 13, 4),
    ExtractTopicProb(physeq_iso_Col0_erd4, beta_lda_VEM_Col0_erd4_k13, 0.02, 13, 5),
    ExtractTopicProb(physeq_iso_Col0_erd4, beta_lda_VEM_Col0_erd4_k13, 0.02, 13, 12),
    
    # model with 30 topics
    ExtractTopicProb(physeq_iso_Col0_erd4, beta_lda_VEM_Col0_erd4_k30, 0.02, 30, 18),
    ExtractTopicProb(physeq_iso_Col0_erd4, beta_lda_VEM_Col0_erd4_k30, 0.02, 30, 22),
    ExtractTopicProb(physeq_iso_Col0_erd4, beta_lda_VEM_Col0_erd4_k30, 0.02, 30, 30)
  ) %>% 
  mutate(topicN = as.factor(topicN)) 


Alltopic_beta_df_genotype_Col0_erd4
head(Alltopic_beta_df_genotype_Col0_erd4)

# plotting 

Fig_S21_beta_genotype_Col0_erd4 <-
  PlotBetaProb(Alltopic_beta_df_genotype_Col0_erd4)

Fig_S21_beta_genotype_Col0_erd4 <-
  Alltopic_beta_df_genotype_Col0_erd4 %>% 
  mutate(k = fct_relevel(k, "k=13", "k=30")) %>% 
  PlotBetaProb(.) +
  #scale_fill_gradient(low = "#ffffe5", high = "#d73027") +
  labs(subtitle = "Topics enriched in erd4 compared to Col-0 genotypes under drought conditions")

Fig_S21_beta_genotype_Col0_erd4

# ***** Figure S21 ***** -------------------------------------------------------
ggsave("figures/Fig_S21_beta_genotype_Col0_erd4.pdf", 
       Fig_S21_beta_genotype_Col0_erd4, device = "pdf")


# Fit linda between Col-0 and erd4 ---------------------------------------------

RunLindaGenotype <- function(physeq, beta_df, k_val, topic_n, 
                             genotype_1, genotype_2) {
  library(tidyverse)
  library(LinDA)
  
  # Clean up the beta table to make sure k is numeric
  beta_clean <- beta_df %>%
    mutate(
      k = as.numeric(gsub("k=", "", k)),
      topicN = as.character(topicN)
    ) %>%
    filter(topic == topic_n & k == k_val)
  
  if (nrow(beta_clean) == 0) {
    stop("No matching taxa found for the specified topic and k value")
  }
  
  # Get the taxa names to keep
  taxa_to_keep <- beta_clean$term
  
  if (length(taxa_to_keep) == 0) {
    stop("No matching taxa found after filtering")
  }
  
  # Subset the OTU table
  otu_tab <- as.data.frame(
    otu_table(
      prune_taxa(taxa_names(physeq) %in% taxa_to_keep, physeq)
    ) %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
      prune_samples(sample_sums(x = .) > 0, x = .)
  )
  
  # Subset the metadata table
  meta <- as.data.frame(as.matrix(
    sample_data(
      prune_taxa(taxa_names(physeq) %in% taxa_to_keep, physeq) %>%
        prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
        prune_samples(sample_sums(x = .) > 0, x = .)))) %>%
    mutate(Genotype = fct_relevel(Genotype, genotype_1, genotype_2))
  
  # Check for missing metadata or low sample count
  if (nrow(meta) == 0) {
    stop("No matching samples found after filtering")
  }
  
  # Fit the `linda` model
  model <- linda(
    otu.tab = otu_tab,
    meta = meta,
    formula = "~ Genotype + (1 | Replicate)",  
    imputation = FALSE,
    p.adj.method = "BH",
    alpha = 0.05,
    n.cores = 6
  )
  
  return(model)
  
}


RunLindaGenotype(physeq_iso_Col0_erd4, 
                 Alltopic_beta_df_genotype_Col0_erd4, 13, 4,"erd4", "Col-0")

ExtractLindaLog2FoldGenotype(
  RunLindaGenotype(physeq_iso_Col0_erd4, 
                   Alltopic_beta_df_genotype_Col0_erd4,13, 4,"erd4","Col-0"),
  physeq_iso_Col0_erd4) 


Pertopic_Fold2Change_Col0_erd4 <-
  rbind(
    # k=10
    ExtractLindaLog2FoldGenotype(
      RunLindaGenotype(physeq_iso_Col0_erd4, 
                       Alltopic_beta_df_genotype_Col0_erd4,13,4,"erd4","Col-0"),
      physeq_iso_Col0_erd4)%>% 
      mutate(k = "k=13", topicN = "4"),
    ExtractLindaLog2FoldGenotype(
      RunLindaGenotype(physeq_iso_Col0_erd4, 
                       Alltopic_beta_df_genotype_Col0_erd4,13,5,"erd4","Col-0"),
      physeq_iso_Col0_erd4)%>% 
      mutate(k = "k=13", topicN = "5"),
    ExtractLindaLog2FoldGenotype(
      RunLindaGenotype(physeq_iso_Col0_erd4, 
                       Alltopic_beta_df_genotype_Col0_erd4,13,12,"erd4","Col-0"),
      physeq_iso_Col0_erd4)%>% 
      mutate(k = "k=13", topicN = "12"),
    # k=30
    ExtractLindaLog2FoldGenotype(
      RunLindaGenotype(physeq_iso_Col0_erd4, 
                       Alltopic_beta_df_genotype_Col0_erd4,30,18,"erd4","Col-0"),
      physeq_iso_Col0_erd4) %>% 
      mutate(k = "k=30", topicN = "18"),
    ExtractLindaLog2FoldGenotype(
      RunLindaGenotype(physeq_iso_Col0_erd4, 
                       Alltopic_beta_df_genotype_Col0_erd4,30,22,"erd4","Col-0"),
      physeq_iso_Col0_erd4) %>% 
      mutate(k = "k=30", topicN = "22"),
    ExtractLindaLog2FoldGenotype(
      RunLindaGenotype(physeq_iso_Col0_erd4, 
                       Alltopic_beta_df_genotype_Col0_erd4,30,30,"erd4","Col-0"),
      physeq_iso_Col0_erd4) %>% 
      mutate(k = "k=30", topicN = "30")
  )

Pertopic_Fold2Change_Col0_erd4

Fig_S22_Fold2Change_Col0_erd4 <-
  Pertopic_Fold2Change_Col0_erd4 %>%
  mutate(log2FoldChange = log2FoldChange * -1) %>% # Invert the logfold sign
  mutate(topicN_char = as.character(topicN)) %>%
  mutate(topicN = fct_relevel(topicN_char, "4", "5", "12", "18", "22", "30")) %>%  # Corrected line
  mutate(k = fct_relevel(k, "k=13", "k=30")) %>% 
  mutate(label_fmt = ifelse(signif,
                            paste0("**", round(log2FoldChange, 2), "**"),  # bold if TRUE
                            as.character(round(log2FoldChange, 2))         # plain if FALSE
  )) %>% 
  ggplot(aes(x = factor(topicN), y = markdown_name, fill = log2FoldChange)) +
  geom_tile(color = "white") +
  #geom_text(aes(label = round(log2FoldChange, 2)), size = 3, color = "black") +
  geom_richtext(
    aes(label = label_fmt),
    fill = NA, label.color = NA, color = "black",  # No background or box border
    size = 3) +
  scale_fill_gradient2(name = "LFC",
                       low = "#d73027", 
                       high = "#313695", 
                       mid = "white", 
                       na.value = "white",
                       limits = c(-1, 1),
                       breaks = c(-1, -0.5, 0, 0.5, 1)) +
  facet_grid(k ~ topicN, scales = "free", space = "free") +
  theme_bw() +
  theme(
    plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_markdown(size = 10, hjust = 0.5),
    strip.text = element_markdown(
      size = 10, face = "bold"),
    axis.text = element_markdown(hjust = 0.5, vjust = 0.5, size = 10),
    strip.background = element_blank(),
    legend.key.height = unit(0.4, "cm"), 
    legend.key.width = unit(0.4, "cm"),
    panel.grid = element_blank()
  ) +
  labs(title = "Differentially abundant syncom members", 
       subtitle = "Topics enriched in erd4 compared to Col-0 genotypes under drought conditions",
       x = "Topic",
       y = NULL,
       fill = "LFC")


Fig_S22_Fold2Change_Col0_erd4

# ***** Figure S13 ***** -------------------------------------------------------
ggsave("figures/Fig_S22_Fold2Change_Col0_erd4.pdf", 
       Fig_S22_Fold2Change_Col0_erd4, device = "pdf")


# ***********************************************-------------------------------
# LDA MODEL FOR ALL SYNCOM MEMBERS ---------------------------------------------
# We can also fit the LinDA model to each genera individually to get some idea 
# of their association with drought we do not directly account for the
# other features (for better or worse).

# FIXED SLOPE ------------------------------------------------------------------

linda_iso_drought <- linda(
  otu.tab = data.frame(otu_table(physeq_iso)),
  meta = data.frame(sample_data(physeq_iso)) %>%
    mutate(Treatment = fct_relevel(Treatment, "w", "d")),
  formula = "~ Treatment + ( 1 | Replicate)",
  # Replicate as random effect
  imputation = FALSE,
  p.adj.method = "BH",
  alpha = 0.05,
  n.cores = 6
)

linda_iso_drought


# Plotting
Fig_SXX_linda_symcom <-
linda_iso_drought$output$Treatmentd %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ZOTU") %>% 
  left_join(., physeq_iso@tax_table %>% 
              as.matrix() %>% 
              as.data.frame() %>% 
              rownames_to_column("ZOTU"),
            by="ZOTU") %>% 
  mutate(color = if_else(pvalue <= 0.05, "#ff7f00", "grey"),
         color = if_else(padj <= 0.05, "#1f78b4", color)) %>% 
  ggplot(aes(x = markdown_name, y = log2FoldChange, fill = color)) +
  geom_col(alpha = 0.8) +
  geom_hline(yintercept = 0, linetype="dotted") +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_markdown(size = 10, face = "italic", hjust = 0.5),
        strip.text = element_markdown(size = 10, face = "bold"),
        axis.text.x = element_markdown(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
        axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
        legend.position = "none") +
  ggplot2::annotate("text",  x= Inf, y = Inf, label = "bold(drought)", 
                    parse = TRUE, hjust=1.1, vjust=1, col="black", size=4) +
  scale_fill_manual(values = c("#1f78b4", "#ff7f00", "grey" )) +
  labs(title = "SynCom members",
       y = "Log2 Fold-Change", 
       x = "") 

Fig_SXX_linda_symcom

Fig_SXX_linda_symcom + 
  scale_y_continuous(limits=c(-0.35, +0.35))

###### Custom figure for Sarah ##### -------------------------------------------

Fig_linda_all_syncom_df <-
  linda_iso_drought$output$Treatmentd %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ZOTU") %>% 
  left_join(., physeq_iso@tax_table %>% 
              as.matrix() %>% 
              as.data.frame() %>% 
              rownames_to_column("ZOTU"),
            by="ZOTU") 

Fig_linda_all_syncom_df

syncom_order <- read.csv(file = "dataset/New SynCom order.csv", header = FALSE)
syncom_order


Fig_linda_all_syncom_df$Isolate
syncom_order$V1

#Clean up both isolate name columns to make them comparable
clean_names <- function(x) {
  x %>%
    str_replace_all("_", "") %>%  # remove underscores
    str_replace_all("\\.", "")    # remove dots
}


df1 <- Fig_linda_all_syncom_df %>%
  mutate(Isolate_clean = clean_names(Isolate))

df2 <- syncom_order %>%
  mutate(Isolate_clean = clean_names(V1),
         Order = row_number())

# 2. Join order information into df1
df1_ordered <- df1 %>%
  left_join(df2 %>% select(Isolate_clean, Order), by = "Isolate_clean") %>%
  mutate(Isolate_clean = factor(Isolate_clean, 
                                levels = df2$Isolate_clean[order(df2$Order)])) %>%
  arrange(Isolate_clean) 

df1_ordered$markdown_name <- as.factor(df1_ordered$markdown_name) 
df1_ordered
  
identical(df2$Isolate_clean, as.character(df1_ordered$Isolate_clean))

# 3. Plot with y-axis following the new order
Fig_3_syncom_linda <- 
  df1_ordered %>% 
  mutate(
    color = if_else(pvalue <= 0.05, "#ff7f00", "grey"),
    color = if_else(padj <= 0.05, "#1f78b4", color),
    markdown_name = factor(markdown_name, levels = rev(unique(markdown_name)))
  ) %>% 
  mutate(markdown_name = recode(markdown_name,
                                "*Pseudomonas* sp. MF397" = "*Pseudomonas umsongensis* CL58"
  )) %>% 
  ggplot(aes(x = markdown_name, y = log2FoldChange, fill = color)) +
  geom_col(alpha = 0.8) +
  geom_hline(yintercept = 0, linetype="dotted") +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_markdown(size = 10, face = "italic", hjust = 0.5),
    strip.text = element_markdown(size = 10, face = "bold"),
    axis.text.x = element_markdown(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
    axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.4, "cm"),
    legend.position = "none"
  ) +
  ggplot2::annotate("text", x = Inf, y = Inf, label = "bold(drought)", 
                    parse = TRUE, hjust = 1.1, vjust = 1, col = "black", size = 4) +
  scale_fill_manual(values = c("#1f78b4", "#ff7f00", "grey")) +
  labs(title = "SynCom members",
       y = "Log2 Fold-Change", 
       x = "")

Fig_3_syncom_linda




# ***** Figure SXX1 ***** -------------------------------------------------------
ggsave("figures/Fig_SXX_linda_symcom.pdf", 
       Fig_SXX_linda_symcom, device = "pdf")

# RANDOM SLOPE -----------------------------------------------------------------

# testing random slope for curiosity
linda_iso_drought_randomslope <- linda(
  otu.tab = data.frame(otu_table(physeq_iso)),
  meta = data.frame(sample_data(physeq_iso)) %>%
    mutate(Treatment = fct_relevel(Treatment, "w", "d")),
  formula = "~ Treatment + ( Treatment | Replicate)",
  # Replicate as random effect
  imputation = FALSE,
  p.adj.method = "BH",
  alpha = 0.05,
  n.cores = 6
)

linda_iso_drought_randomslope$output$Treatmentd %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ZOTU") %>% 
  left_join(., physeq_iso@tax_table %>% 
              as.matrix() %>% 
              as.data.frame() %>% 
              rownames_to_column("ZOTU"),
            by="ZOTU") %>% 
  mutate(color = if_else(pvalue <= 0.05, "#ff7f00", "grey"),
         color = if_else(padj <= 0.05, "#1f78b4", color)) %>% 
  ggplot(aes(x = markdown_name, y = log2FoldChange, fill = color)) +
  geom_col(alpha = 0.8) +
  geom_hline(yintercept = 0, linetype="dotted") +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_markdown(size = 10, face = "italic", hjust = 0.5),
        strip.text = element_markdown(size = 10, face = "bold"),
        axis.text.x = element_markdown(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
        axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
        legend.position = "none") +
  ggplot2::annotate("text",  x= Inf, y = Inf, label = "bold(drought)", 
                    parse = TRUE, hjust=1.1, vjust=1, col="black", size=4) +
  scale_fill_manual(values = c("#1f78b4", "#ff7f00", "grey" )) +
  labs(title = "SynCom members",
       y = "Log2 Fold-Change", 
       x = "") 

# Not really worth it!

# RUNNING LINDA on all Syncom members for Col-0 vs. pad4 -----------------------

# Running this analysis on just the two selected geneotypes in drought conditions.

linda_iso_genotype_Col0_pad4 <- linda(
  otu.tab = data.frame(otu_table(physeq_iso_Col0_pad4)),
  meta = data.frame(sample_data(physeq_iso_Col0_pad4)) %>%
    mutate(Genotype = fct_relevel(Genotype, "pad4", "Col-0")),
  formula = "~ Genotype + (1|Replicate)", 
  # Replicate as random effect
  imputation = FALSE,
  p.adj.method = "BH",
  alpha = 0.05,
  n.cores = 6
)

linda_iso_genotype_Col0_pad4

Fig_SXX_linda_syncom_Col0_pad4 <-
linda_iso_genotype_Col0_pad4$output$`GenotypeCol-0` %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ZOTU") %>% 
  left_join(., physeq_iso_Col0_pad4@tax_table %>% 
              as.matrix() %>% 
              as.data.frame() %>% 
              rownames_to_column("ZOTU"),
            by="ZOTU") %>% 
  mutate(color = if_else(pvalue <= 0.05, "#ff7f00", "grey"),
         color = if_else(padj <= 0.05, "#1f78b4", color)) %>% 
  ggplot(aes(x = markdown_name, y = log2FoldChange, fill = color)) +
  geom_col(alpha = 0.8) +
  geom_hline(yintercept = 0, linetype="dotted") +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_markdown(size = 10, face = "italic", hjust = 0.5),
        strip.text = element_markdown(size = 10, face = "bold"),
        axis.text.x = element_markdown(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
        axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
        legend.position = "none") +
  ggplot2::annotate("text",  x= Inf, y = Inf, label = "bold(Col-0)", 
                    parse = TRUE, hjust=1.4, vjust=1.4, col="black", size=3) +
  scale_fill_manual(values = c("#1f78b4", "#ff7f00", "grey" )) +
  labs(title = "SynCom members",
       y = "Log2 Fold-Change", 
       x = "") 


Fig_SXX_linda_syncom_Col0_pad4

# ***** Figure SXX2 ***** -------------------------------------------------------
ggsave("figures/Fig_SXX_linda_symcom_Col0_pad4.pdf", 
       Fig_SXX_linda_syncom_Col0_pad4, device = "pdf")



# RUNNING LINDA on all Syncom members for Col-0 vs. erd4 -----------------------

linda_iso_genotype_Col0_erd4 <- linda(
  otu.tab = data.frame(otu_table(physeq_iso_Col0_erd4)),
  meta = data.frame(sample_data(physeq_iso_Col0_erd4)) %>%
    mutate(Genotype = fct_relevel(Genotype, "erd4", "Col-0")),
  formula = "~ Genotype + (1|Replicate)", 
  # Replicate as random effect
  imputation = FALSE,
  p.adj.method = "BH",
  alpha = 0.05,
  n.cores = 6
)

linda_iso_genotype_Col0_erd4

Fig_SXX_linda_syncom_Col0_erd4 <-
  linda_iso_genotype_Col0_erd4$output$`GenotypeCol-0` %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ZOTU") %>% 
  left_join(., physeq_iso_Col0_erd4@tax_table %>% 
              as.matrix() %>% 
              as.data.frame() %>% 
              rownames_to_column("ZOTU"),
            by="ZOTU") %>% 
  mutate(color = if_else(pvalue <= 0.05, "#ff7f00", "grey"),
         color = if_else(padj <= 0.05, "#1f78b4", color)) %>% 
  ggplot(aes(x = markdown_name, y = log2FoldChange, fill = color)) +
  geom_col(alpha = 0.8) +
  geom_hline(yintercept = 0, linetype="dotted") +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_markdown(size = 10, face = "italic", hjust = 0.5),
        strip.text = element_markdown(size = 10, face = "bold"),
        axis.text.x = element_markdown(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
        axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
        legend.position = "none") +
  ggplot2::annotate("text",  x= Inf, y = Inf, label = "bold(Col-0)", 
                    parse = TRUE, hjust=1.4, vjust=1.4, col="black", size=3) +
  scale_fill_manual(values = c("#1f78b4", "#ff7f00", "grey" )) +
  labs(title = "SynCom members",
       y = "Log2 Fold-Change", 
       x = "") 


Fig_SXX_linda_syncom_Col0_erd4

# ***** Figure SXX2 ***** -------------------------------------------------------
ggsave("figures/Fig_SXX_linda_symcom_Col0_erd4.pdf", 
       Fig_SXX_linda_syncom_Col0_erd4, device = "pdf")


# ***********************************************-------------------------------
# TOPIC MODEL FIGUREs for Manuscript --------------------------------------------

topics_k17 <-
  predicted_topics %>%
  filter(Condition == "k=17") %>% 
  group_by(Condition) %>%
  complete(Topic = paste("Topic", 1:28), fill = list(Log2_Fold_Change = NA)) %>%
  mutate(Topic = fct_relevel(Topic, paste("Topic", 1:max(as.numeric(gsub("Topic ", "", Topic)))))) %>%
  mutate(Topic = fct_rev(Topic)) %>%  # Reverse the order of topics
  ungroup() %>%
  #ggplot(aes(x = Topic, y = log2FoldChange, fill = reject)) +
  mutate(color = if_else(pvalue <= 0.05, "#ff7f00", "grey"),
         color = if_else(padj <= 0.05, "#1f78b4", color)) %>% 
  ggplot(aes(x = Topic, y = log2FoldChange, fill = color)) +
  geom_col(width = 0.8) + # Set a reasonable width to control thickness of bars
  facet_wrap(~Condition, ncol = 3, scales = "free_y") + # Create the facets for each group
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Dashed line at y = 0
  coord_flip() + # Flip coordinates to have horizontal bars
  theme_bw() +
  annotate("text", x= Inf, y = Inf, label = "bold(drought)", 
           parse = TRUE, hjust=1.4, vjust=1.4, col="black", size = 3) +
  theme(plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_markdown(size = 10, hjust = 0.5),
        strip.text = element_markdown(size = 12, face = "bold"),
        axis.text.x = element_markdown(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
        axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
        legend.position = "none") +
  scale_fill_manual(values = c("#1f78b4", "#ff7f00","grey")) +
  labs(title = "Predicted topics", x = "", y = "Log2 Fold-Change")

topics_k17

beta_prob_k17 <-
  Alltopic_beta_df_drought %>% 
  filter(k == "k=17") %>% 
  PlotBetaProb() +
  #scale_fill_gradient(low = "#ffffe5", high = "#DDCC77") +
  scale_fill_gradient(low = "#ffffe5", high = "grey40") +
  #scale_fill_gradient(low = "#ffffe5", high = "#d73027") +
  labs(title = "SynCom members Beta probilities in topics enriched in drought",
       x = NULL)

beta_prob_k17


fold2Change_drought_k17 <-
  Pertopic_Fold2Change_drought %>%
  as.data.frame() %>%
  filter(k == "k=17") %>% 
  mutate(topicN = as.factor(topicN),
         topicN = fct_relevel(topicN, "7", "8", "9", "15")) %>% 
  mutate(label_fmt = ifelse(signif,
                            paste0("**", round(log2FoldChange, 2), "**"),  # bold if TRUE
                            as.character(round(log2FoldChange, 2))         # plain if FALSE
  )) %>% 
  ggplot(aes(x = factor(topicN), y = markdown_name, fill = log2FoldChange)) +
  geom_tile(color = "white") +
  #geom_text(aes(label = round(log2FoldChange, 2)), size = 3, color = "black") +
  geom_richtext(
    aes(label = label_fmt),
    fill = NA, label.color = NA, color = "black",  # No background or box border
    size = 3) +
  scale_fill_gradient2(name = "LFC",
                       low = "#313695", 
                       high = "#d73027", 
                       mid = "white", 
                       na.value = "white",
                       limits = c(-1, 1),
                       breaks = c(-1, -0.5, 0, 0.5, 1)) +
  facet_grid(k ~ topicN, scales = "free", space = "free") +
  theme_bw() +
  theme(
    plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_markdown(size = 10, hjust = 0.5),
    strip.text = element_markdown(
      size = 10, face = "bold"),
    axis.text = element_markdown(hjust = 0.5, vjust = 0.5, size = 10),
    strip.background = element_blank(),
    legend.key.height = unit(0.4, "cm"), 
    legend.key.width = unit(0.4, "cm"),
    panel.grid = element_blank()
  ) +
  labs(title = "Differentially abundant SynCom members",
       x = NULL,
       y = NULL,
       fill = "LFC")


fold2Change_drought_k17

Fig_4_Topic_Models <- 
ggarrange(topics_k17,
    ggarrange(
      beta_prob_k17 +
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()),
      fold2Change_drought_k17 +
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()),
      ncol=1, 
      nrow=2,
      labels=c("B", "C")
    ),
ncol=2,
widths=c(1, 3),
labels="A"
)

Fig_4_Topic_Models

# ***** Figure SXX2 ***** -------------------------------------------------------
ggsave("figures/Fig_4_Topic_Models.pdf", 
       Fig_4_Topic_Models, device = "pdf")









# ***********************************************-------------------------------
# CO-OCCURRENCE ANALYSIS -------------------------------------------------------
library(cooccur)

asv_iso <- as.data.frame(as.matrix(otu_table(physeq_iso))) %>% 
  rownames_to_column("OTU_ID") %>% 
  left_join(physeq_iso@tax_table %>% 
              as.matrix() %>% 
              as.data.frame() %>% 
              rownames_to_column("OTU_ID"),
            by = "OTU_ID") %>% 
  select(starts_with("sample"), new_name) %>% 
  column_to_rownames("new_name")


head(asv_iso)


asv_sc <- as.data.frame(as.matrix(otu_table(physeq_sc)))%>% 
  rownames_to_column("OTU_ID") %>% 
  left_join(physeq_sc@tax_table %>% 
              as.matrix() %>% 
              as.data.frame() %>% 
              rownames_to_column("OTU_ID"),
            by = "OTU_ID") %>%
  mutate(Isolate = gsub(" ", "_", Isolate)) %>%
  mutate(new_name = paste(OTU_ID, Isolate, sep="_")) %>% 
  select(starts_with("sample"), new_name) %>% 
  column_to_rownames("new_name")
 
asv_sc

cooccur_syncom_iso <- cooccur(asv_iso, type = "spp_site",
                              thresh = TRUE, spp_names = TRUE)
cooccur_syncom_iso

cooccur_syncom_iso$

  

pair(mod = cooccur_syncom_iso, spp = "Ralstonia sp. CL21")

summary(cooccur_syncom_iso)

plot(cooccur_syncom_iso)


cooccur_syncom_sc <- cooccur(asv_sc, type = "spp_site",
                              thresh = TRUE, spp_names = TRUE)
cooccur_syncom_sc



