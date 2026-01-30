#*************************************************************************************-------
# Manuscript Title: Site-specific factors rather than nitrogen impact arbuscular 
#                   mycorrhizal fungi diversity in bioenergy switchgrass monocultures
# Authors:          Shuang Liu, Gian Maria Niccolò Benucci, Alden Dirks, 
#                   Lukas Bell-Dereske, Sarah Evans, Gregory Bonito
# Code Developer:   Gian MN Benucci 2025
# Citation:         ...
#                   
# DOI               ...
# PMID:             ...
# ************************************************************************************--------

# Positron options to restore Rstudio projects
load(".Rdata")

# R options --------------------------------------------------------------------
options(scipen = 9999, pillar.sigfig = 6, digits = 6, max.print = 10000000)
#rm(list = ls())

# Check the lib paths ----------------------------------------------------------
.libPaths()

# Load R packages --------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  pak,
  styler,
  magrittr,
  Biostrings,
  ape,
  msa,
  DECIPHER,
  phangorn,
  decontam,
  phyloseq,
  speedyseq,
  tidysq,
  vegan,
  AICcPermanova,
  tidyverse,
  ggtext,
  ggpubr,
  cowplot,
  gridExtra,
  ggrepel,
  scales,
  agricolae,
  install=FALSE
)

# Session Info -----------------------------------------------------------------
sessionInfo()

# **********************************************************************--------
# ***** PATHS ***** ------------------------------------------------------------

data_path <-
  ("/home/gian/Dropbox/6_PROJETCS/2025_amf_community_pacbio")

data_path

# results ----------------------------------------------------------------------

results_path <-
  ("/home/gian/Dropbox/6_PROJETCS/2025_amf_community_pacbio/github")

results_path

# **********************************************************************--------
# ***** IMPORT ***** -----------------------------------------------------------

# Decided to work with the 99% clustered ASVs.

# otu table --------------------------------------------------------------------

otutable_99 <-read.delim(file.path(data_path, "datasets/otutab_99.txt"), 
                         row.names = 1)
head(otutable_99)

# NCBI blasTAX taxonomy --------------------------------------------------------

taxonomy_99 <- extract_blasTAX(
  tax_path = file.path(data_path, "datasets/taxonomy_blast_99.txt"),
  namemap_path = file.path(data_path, "datasets/name_mapping_99.txt")
)

head(taxonomy_99)
dim(taxonomy_99)

# These below were Unclassified  99% OTUs
taxonomy_99 %>% filter(Query %in% "Query_998")
taxonomy_99 %>% filter(Query %in% "Query_1124")

# Adding columns and finalize taxonomy -----------------------------------------
taxonomy_99_fix <-
  FinalizeTaxonomy(
    taxonomy_99 %>%
      dplyr::select(
        "Zotu",
        "Query",
        "Kingdom" ,
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species"
      )
  ) %>%
  full_join(taxonomy_99 %>% dplyr::select(Query, S_score), by = "Query") %>%
  mutate(S_score = as.numeric(S_score)) %>%
  mutate_if(is.character, ~ replace(., is.na(.), "Unclassified")) %>%
  column_to_rownames("Zotu") %>%
  filter(Phylum %in% "Mucoromycota")

head(taxonomy_99_fix)
as.factor(taxonomy_99_fix$Phylum)


# 99% OTUs ---------------------------------------------------------------------

otu_99 <- readDNAStringSet(
  file.path(data_path, "datasets/otus_99.fasta"),
  format = "fasta",
  seek.first.rec = TRUE,
  use.names = TRUE
)
otu_99
class(otu_99)

# metadata  --------------------------------------------------------------------
metadata_99 <-
  read.csv(
    file = file.path(data_path, "datasets/metadata_pacbio.csv"),
    header = TRUE
  ) %>%
  column_to_rownames("SampleID") %>%
  janitor::clean_names()

head(metadata_99)

# Phylogenetic tree ------------------------------------------------------------

# I run this analysis on the HPCC 

# options
cat("Setting up R options.\n\n")
options(scipen = 9999)

# libraries
cat("Loading R libraries.\n\n")
library(tidyverse)
library(msa)
library(phangorn)
library(ape)
library(DECIPHER)
library(future)
library(future.apply)

# path
cat("Load data PATH.\n\n")
data_path = "/mnt/home/benucci/Phylogeny_in_R/otu_alignment/"
results_path = "/mnt/home/benucci/Phylogeny_in_R/results/"

# import
cat("Import otu representative sequences file into R.\n\n")
otu_99 <- readDNAStringSet(
  file.path(data_path, "otus_99.fasta"),
  format = "fasta",
  seek.first.rec = TRUE,
  use.names = TRUE
)

# Set up parallel processing
cat("Set up parallel processing...\n\n")
plan("multicore", workers = 42)
cat("Parallel workers:", future::nbrOfWorkers(), "\n\n")

# align
cat("Align sequences in parallel...\n\n")
results <- future_lapply(list("simple", "custom", "desc"), function(type) {
  if (type == "simple") {
    msa(otu_99, method = "Muscle")
  } else if (type == "custom") {
    msa(
      otu_99,
      method = "Muscle",
      gapOpening = -10,
      gapExtension = -0.5,
      maxiters = 2,
      cluster = "upgmamax",
      order = "input",
      verbose = TRUE
    )
  } else if (type == "decipher") {
    # DECIPHER alignment (directly on DNAStringSet)
    AlignSeqs(otu_99, processors = NULL)  # NULL = all available cores
  }
}, future.seed = TRUE)

aligned_otu_99_simple <- results[[1]]
aligned_otu_99_custom <- results[[2]]
aligned_otu_99_desc <- results[[3]]

# save
cat("All done! Save results.\n\n")
saveRDS(object = aligned_otu_99_simple,
        file = file.path(results_path, "alignment_simple.RDS"))
saveRDS(object = aligned_otu_99_custom,
        file = file.path(results_path, "alignment_custom.RDS"))
saveRDS(object = aligned_otu_99_desc,
        file = file.path(results_path, "alignment_desc.RDS"))

# Import alignments -------------------------------------------------------------
# Made 3 different alingments with different parameters 

# Convert MSA object to phangorn format
align_simple <-
  readRDS(file = file.path(data_path, "sequence_alignments/alignment_simple.RDS")) %>% 
  as.phyDat(type = "DNA") 

str(align_simple)

align_custom <-
  readRDS(file = file.path(data_path, "sequence_alignments/alignment_custom.RDS")) %>%
  as.phyDat(type = "DNA") 

# This is deshipher
align_deshipher <-
  readRDS(file = file.path(data_path, "sequence_alignments/alignment_desc.RDS")) 

#Calculate distance matrix (for starting tree)
dm_simple <- dist.ml(align_simple)

# Create a starting tree (NJ tree)
NJ_tree_simple <- NJ(dm_simple)

# Optimize tree with Maximum Likelihood
ML_tree_simple <- pml(NJ_tree_simple, align_simple)

# Optimize tree parameters
ML_tree_simple_optimized <- optim.pml(ML_tree_simple,
                           optNni = TRUE,      # Optimize topology
                           optBf = TRUE,       # Optimize base frequencies
                           optQ = TRUE,        # Optimize rate matrix
                           optGamma = TRUE,    # Optimize gamma shape
                           model = "GTR")      # Use GTR model

# Extract the ML tree
ML_tree_simple_final <- ML_tree_simple_optimized$tree

# Plot the tree
plot(ML_tree_simple_final, main = "Maximum Likelihood Tree")

ggtree(ML_tree_simple_final) + 
  geom_tiplab(size = 3) +
  geom_nodelab(aes(label = label), size = 2, color = "red") +
  theme_tree2()


# *********************************************---------------------------------
# ***** MAKE PHYLOSEQ OBJECT ***** ---------------------------------------------
dim(otutable_99)
dim(metadata_99)
dim(taxonomy_99_fix)
str(otu_99)
str(raxml_tree_filtered)

rownames(otutable_99)
rownames(taxonomy_99_fix)
names(otu_99)
raxml_tree_filtered$tip.label


physeq_AMF <-
  phyloseq(
    otu_table(otutable_99, taxa_are_rows = TRUE),
    sample_data(metadata_99),
    tax_table(as.matrix(taxonomy_99_fix)),
    DNAStringSet(otu_99), 
    phy_tree(raxml_tree_filtered)) %>% 
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.)

physeq_AMF

physeq_AMF@phy_tree
head(physeq_AMF@otu_table)
head(physeq_AMF@sam_data)

is.na(physeq_AMF@sam_data$collection_date)

physeq_AMF %>%
  subset_samples(!is.na(collection_date))

dim(Zotutable_99)
dim(metadata_99)
dim(blastax_99)

# Intersection across *all* components
taxa_all <- Reduce(
  intersect,
  list(
    rownames(otutable_99),
    rownames(taxonomy_99_fix),
    names(otu_99),
    raxml_tree_filtered$tip.label
  )
)
length(taxa_all)  # this will almost certainly be 243

length(setdiff(rownames(otutable_99), taxa_all))
head(setdiff(rownames(otutable_99), taxa_all))

length(setdiff(rownames(taxonomy_99_fix), taxa_all))
head(setdiff(rownames(taxonomy_99_fix), taxa_all))

length(setdiff(names(otu_99), taxa_all))
head(setdiff(names(otu_99), taxa_all))

length(setdiff(raxml_tree_filtered$tip.label, taxa_all))
head(setdiff(raxml_tree_filtered$tip.label, taxa_all))

# **** BUILD THE PHYLOEGENTY TREE **** -----------------------------------------
library(msa)
library(phangorn)
library(ape)
library(job)

# This code will run in a background RStudio job

job::job({
  aligned_otu_99 <- msa(otu_99, method = "Muscle") 
  print(paste("Background job result:", aligned_otu_99))
})

aligned <- aligned_msa@unmasked   









##################################################################--------------
##################################################################--------------
##################################################################--------------

physeq_glo99 <-
  phyloseq(
    otu_table(Zotutable_99, taxa_are_rows = TRUE),
    sample_data(metadata_99),
    tax_table(as.matrix(blastax_99_fix)),
    Zotu_99,
    raxml_tree_filtered) %>% 
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.)

physeq_glo99

phy_tree(physeq_glo99)



# Filter out control samples ---------------------------------------------------
head(physeq_glo99@sam_data)
sort(sample_sums(physeq_glo99))

# Not sure why the MMPRNTCtrl7BB5 has 15158 reads
physeq_AMF_filt <-
  subset_samples(physeq_glo99, is.na(collectionDate) == FALSE) %>% 
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.)

physeq_AMF_filt

# ***** OPTIMAL RAREFACTION DEPTH ***** ----------------------------------------
# Calculating Good's coverage. The fraction of sequences that appear in an OTU
# that have been seen more than once, and allows estimating what percent of the total 
# species is represented in a sample.
# Coverage = 1 - (number of individuals in species / total number of individuals)
# Example: If I have Goods = 0.96 it means that 4% of your reads in that sample are
# from OTUs that appear only once in that sample.

# Calculate long dataframe with stats 
RareStats <- function(physeq){
  require(tidyverse)
  # calculate distribution outliers
  findoutlier <- function(x) {
    return(x < quantile(x, 0.25) - 1.5*IQR(x) | x > quantile(x, 0.75) + 1.5*IQR(x))
  }
  # generate output dataframe
  df_output <-
    as.data.frame(as.matrix(physeq@otu_table)) %>% 
    rownames_to_column("OTU_ID") %>% 
    pivot_longer(-OTU_ID, names_to = "SampleID", values_to = "Seq_No") %>% 
    group_by(SampleID) %>% 
    summarize(ReadNo = sum(Seq_No),
              Singlton_No = sum(Seq_No == 1),
              Goods = 100*(1 - Singlton_No / ReadNo)) %>% 
    mutate(outlier = ifelse(findoutlier(log10(ReadNo)), ReadNo, NA)) %>% 
    ungroup()
  
  return(df_output)
}


# Plotting fucntion 
PloRareStats <- function(dataframe){
  require(ggrepel)
  plot_output <-
    ggarrange(
      dataframe %>%
        ggplot(aes(x = ReadNo)) +
        geom_histogram(binwidth = 5000, 
                       fill = "firebrick", color = "firebrick") +
        theme_bw() +
        theme(axis.text.x = element_markdown(angle = 33, hjust = 1, vjust = 1)) +
        labs(title = "Histogram"), 
      dataframe %>%
        ggplot(aes(x = ReadNo)) +
        geom_histogram(binwidth = 1000, 
                       fill = "firebrick",color = "firebrick") +
        coord_cartesian(xlim = c(0, 25000)) +
        theme_bw() +
        theme(axis.text.x = element_markdown(angle = 33, hjust = 1, vjust = 1)) +
        labs(title = "Histogram Zoom"), 
      dataframe %>%
        ggplot(aes(x = ReadNo, y = Goods)) +
        geom_point(shape = 1, color = "firebrick") +
        theme_bw() +
        theme(axis.text.x = element_markdown(angle = 33, hjust = 1, vjust = 1)) +
        labs(title = "Good's Coverage"), 
      dataframe %>%
        ggplot(aes(x = 1, y = ReadNo)) +
        geom_jitter(shape=1, color = "firebrick") +
        theme_bw() +
        theme(axis.text.x = element_markdown(angle = 33, hjust = 1, vjust = 1)) +
        scale_y_log10() +
        labs(title = "Log10 jitter",
             x = "Data set"), 
      dataframe %>%
        ggplot(aes(x = 1, y = ReadNo)) +
        geom_boxplot(color = "firebrick") +
        theme_bw() +
        theme(axis.text.x = element_markdown(angle = 33, hjust = 1, vjust = 1)) +
        scale_y_log10() +
        geom_text_repel(
          data = filter(dataframe, !is.na(outlier)),
          mapping = aes(x = 1, y = ReadNo, label = outlier), 
          max.overlaps = 15, size = 3) +
        labs(title = "Log10 boxplot",
             x = "Data set"), 
      dataframe %>%
        arrange(ReadNo) %>%
        ggplot(aes(x = 1:nrow(.), y = ReadNo)) +
        geom_bar(stat = "identity", color = "firebrick") +
        theme_bw() +
        theme(axis.text.x = element_markdown(angle = 33, hjust = 1, vjust = 1)) +
        labs(title = "Ranked",
             x = "Samples"), 
      ncol = 3,
      nrow = 2, 
      align = "hv", 
      labels = c("a", "b", "c", "d", "e", "f"))
  return(plot_output)
}


RareStats(physeq_AMF_filt) %>% dplyr::arrange(ReadNo)

# >>>> FIGURE S1 <<<< ---------------------------------------------------------- 
# rarefaction depth ------------------------------------------------------------
Fig_S1_rarefaction <- 
  PloRareStats(RareStats(physeq_AMF_filt))
Fig_S1_rarefaction

grid.arrange(Fig_S1_rarefaction, top = text_grob(
  "Pre-rarefaction diagnostic metrics", size = 14, face = 2
))

ggsave(
  plot = grid.arrange(Fig_S1_rarefaction, top = text_grob(
    "Pre-rarefaction diagnostic metrics", size = 14, face = 2
  )),
  path = results_path,
  filename = "Fig_S1_rarefaction.pdf"
)

#************* RAREFACTION CURVES ************----------------------------------
rarecurve_glo <- 
  physeq_AMF_filt@otu_table %>% 
  t() %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rarecurve(x=., step=1000, tidy = TRUE)

head(rarecurve_glo)
head(physeq_AMF_filt@sam_data)

# Then plotting rarecurve in ggplot2
PlotRareCurve <- function(rare_tidy, depth=NULL, Col){
  plot_rare <-
    ggplot(rare_tidy, aes(x=Sample, y=Species, group=sample_id, color=get(Col))) +  
    geom_line() +
    geom_vline(xintercept = depth, color="black", linetype = "dashed") +
    theme_bw() +
    theme(
      plot.title = element_markdown(size =12, face = "bold",hjust = 0.5, vjust = 0.5),
      plot.subtitle = element_markdown(size = 10, face = "bold",hjust = 0.5, vjust = 0.5),
      axis.text.x = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
      axis.text.y = element_markdown(size = 8),
      axis.title = element_markdown(size = 10, face = "bold",hjust = 0.5, vjust = 0.5),
      legend.key.height = unit(0.5, "cm"),
      legend.key.width = unit(0.5, "cm"),
      legend.title = element_blank(), 
      legend.text = element_text(size = 8),
      legend.position = "right")
  return(plot_rare)
}

fungi_depth = 9929

# **** Figure S2 - rarefaction curves **** -------------------------------------
Fig_S2_rare_curves <- ggarrange(
  rarecurve_glo %>% 
    rename("sample_id" = Site) %>% 
    left_join(.,
              physeq_AMF_filt@sam_data %>% 
                as.matrix() %>% 
                as.data.frame() %>% 
                rownames_to_column("sample_id"),
              by= "sample_id") %>% 
    PlotRareCurve(., fungi_depth, "site_id") +
    labs(title="Glomeromycotina", x= "Number of DNA reads", y= "Number of OTUs"),
  rarecurve_glo %>% 
    rename("sample_id" = Site) %>% 
    left_join(.,
              physeq_AMF_filt@sam_data %>% 
                as.matrix() %>% 
                as.data.frame() %>% 
                rownames_to_column("sample_id"),
              by= "sample_id") %>% 
    PlotRareCurve(., fungi_depth, "fert_status") +
    labs(title="Glomeromycotina", x= "Number of DNA reads", y= "Number of OTUs"),
  ncol = 2,
  nrow = 1, 
  labels = c("a", "b"),
  common.legend = FALSE,
  align = "hv",
  legend = "bottom")

grid.arrange(Fig_S2_rare_curves, top = text_grob(
  "Rarefaction curves", size = 14, face = 2
))

ggsave(
  plot = grid.arrange(Fig_S2_rare_curves, top = text_grob(
    "Rarefaction curves", size = 14, face = 2
  )),
  path = results_path,
  filename = "Fig_S2_rare_curves.pdf"
)


#************** RAREFACTION ****************------------------------------------
rarefyData <- function(physeq, depth_level){
  require(tidyverse)
  
  dataframe <- 
    as.data.frame(as.matrix(t(physeq@otu_table)))
  
  com_iter <- vector(mode = "list", length =  100)
  
  for (i in seq_along(com_iter)) {
    com_iter[[i]] <- as.data.frame(
      vegan::rrarefy(dataframe, sample = depth_level)
    ) %>% rownames_to_column("SampleID")
  }
  
  mean_100 <- do.call(rbind, com_iter)
  mean_100 <- mean_100 %>% 
    group_by(SampleID) %>% 
    summarise(across(everything(), mean)) %>% 
    filter(rowSums(across(where(is.numeric))) >= depth_level)
  
  print(mean_100 %>% as_tibble())
  return(mean_100)
}

# AMF ------------------------------------------------------------------------

# recreating the phyloseq objects with the rarefied otus
set.seed(2332)

physeq_glo_rare <-
  phyloseq(
    otu_table(rarefyData(physeq_AMF_filt, fungi_depth) %>% 
                column_to_rownames("SampleID") %>% 
                t() %>% 
                as.matrix() %>% 
                as.data.frame(), taxa_are_rows = TRUE),
    physeq_AMF_filt@sam_data,
    physeq_AMF_filt@tax_table,
    physeq_AMF_filt@refseq,
    physeq_AMF_filt@phy_tree) %>% 
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.)

physeq_glo_rare

sample_sums(physeq_glo_rare)
head(physeq_glo_rare@sam_data)
physeq_glo_rare@otu_table[100:110, 40:50]

# subset_samples(physeq_glo_rare, New_Name == "Sample178")@sam_data

# Adding alpha metrics ---------------------------------------------------------
AlphaMetrics <- function(physeq) {
  sample_data(physeq)$ReadNo <- sample_sums(physeq)
  sample_data(physeq)$hill_0 <- as.data.frame(as.matrix(t(physeq@otu_table))) %>% renyi(scales = c(0), hill = TRUE)
  sample_data(physeq)$hill_1 <- as.data.frame(as.matrix(t(physeq@otu_table))) %>% renyi(scales = c(1), hill = TRUE)
  sample_data(physeq)$hill_2 <- as.data.frame(as.matrix(t(physeq@otu_table))) %>% renyi(scales = c(2), hill = TRUE)
  sample_data(physeq)$pielou_j <- log(sample_data(physeq)$hill_1) / log(sample_data(physeq)$hill_0)
  return(physeq)
}

# Where Pileau J is the classic 0–1 evenness measure based on Shannon.

physeq_glo_rare <- 
  AlphaMetrics(physeq_glo_rare)
physeq_glo_rare

head(physeq_glo_rare@sam_data)
AlphaMetrics(physeq_glo_rare) %>% sam_data() %>% pull(pielou_j) %>% range()

# **** ALPHA DIVERSITY **** ----------------------------------------------------
# Testing for significant differences 

library(rstatix)

physeq_glo_rare@sam_data %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(hill_0 = as.numeric(hill_0)) %>%
  as_tibble() %>%
  wilcox_test(formula(hill_0 ~ siteID)) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

physeq_glo_rare@sam_data %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(hill_0 = as.numeric(hill_0)) %>%
  as_tibble() %>%
  mutate(interaction_group = interaction(siteID, fertStatus)) %>%
  wilcox_test(formula(hill_0 ~ interaction_group)) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()



get_pairwise_letters <- function(df, 
                                 formula, 
                                 method = "wilcox.test", 
                                 p.adjust = "BH", 
                                 threshold = 0.05) {
  require(rstatix)
  require(multcompView)
  
  # Extract variable names from formula
  vars <- all.vars(formula)
  response <- vars[1]
  group <- vars[2]
  
  # Pairwise tests
  results <- df %>%
    wilcox_test(formula) %>%
    adjust_pvalue(method = p.adjust) %>%
    add_significance()
  
  # Convert to named vector
  pvals <- setNames(
    results$p.adj,
    paste(results$group1, results$group2, sep = "-")
  )
  
  # Generate letters
  letters <- multcompLetters(pvals, threshold = threshold)
  
  # Return as dataframe
  letters_df <- data.frame(
    group = names(letters$Letters),
    letters = letters$Letters,
    row.names = NULL
  )
  
  colnames(letters_df)[1] <- group
  
  return(list(
    results = results,
    letters = letters_df
  ))
}

alpha_siteID <- 
  get_pairwise_letters(
  df = physeq_glo_rare@sam_data %>%
    as.matrix() %>%
    as.data.frame() %>%
    mutate(hill_0 = as.numeric(hill_0)),
  formula = formula(hill_0 ~ siteID)
)

# Test plotting
# Get summary stats for plotting
alpha_siteID_summary_stats <- 
  physeq_glo_rare@sam_data %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(hill_0 = as.numeric(hill_0)) %>%
  select(hill_0, siteID, fertStatus) %>%
  group_by(siteID) %>%
  summarise(mean = mean(hill_0),
            se = sd(hill_0) / sqrt(n()),
            median = median(hill_0)) %>%
  left_join(alpha_siteID$letters, by = "siteID")

alpha_siteID_summary_stats


PlotRich <- function(physeq, 
                     X_var, 
                     Y_var, 
                     cld_df) {
  # Removed: require(phyloseq)
  # Removed: require(tidyverse)
  
  # Use !!rlang::sym() for unquoting if developing a true package, but 
  # the existing !!sym() works fine when tidyverse is not loaded via library()
  # as long as dplyr/ggplot2 are explicitly referenced.
    
    # 1. Data extraction and preparation
    dataframe <-
      phyloseq::sample_data(physeq) %>%
      as.data.frame() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(!!Y_var := as.numeric(!!rlang::sym(Y_var)))
    
    # DEBUG: Check if X_var exists in the data
    if (!X_var %in% colnames(dataframe)) {
      stop(paste0("Column '", X_var, "' not found in sample_data. Available columns: ", 
                  paste(colnames(dataframe), collapse = ", ")))
    }
    
    # 2. Prepare cld_df with correct column name
    # Get the first column name from cld_df (should be the grouping variable)
    cld_group_col <- colnames(cld_df)[1]
    
    # Rename it to match X_var for the join
    cld_df_renamed <- cld_df %>%
      dplyr::rename(!!X_var := !!rlang::sym(cld_group_col))
    
    # 3. Calculate statistics for label placement
    label_data <- dataframe %>%
      dplyr::group_by(!!rlang::sym(X_var)) %>%
      dplyr::summarise(
        max_y = max(!!rlang::sym(Y_var), na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      # Join with renamed cld_df
      dplyr::left_join(cld_df_renamed, by = X_var) %>%
      dplyr::mutate(
        label_y_position = max_y + 0.05 * max(dataframe[[Y_var]], na.rm = TRUE)
      )
    
    # 4. Plotting
    rich_plot <-
      ggplot2::ggplot(dataframe, ggplot2::aes(x = !!rlang::sym(X_var), y = !!rlang::sym(Y_var))) +
      
      ggplot2::geom_jitter(
        ggplot2::aes(color = !!rlang::sym(X_var)),
        position = ggplot2::position_jitter(width = 0.3, seed = 123), 
        size = 2.5, 
        shape = 16,
        alpha = 0.7
      ) +
      
      ggplot2::stat_summary(
        fun = median, 
        fun.min = median, 
        fun.max = median,
        geom = "crossbar", 
        width = 0.6, 
        linewidth = 0.5, 
        color = "black"
      ) +
      
      # Add IQR bars
      ggplot2::stat_summary(
        fun.min = function(z) quantile(z, 0.25),
        fun.max = function(z) quantile(z, 0.75),
        geom = "errorbar",
        width = 0.3,
        linewidth = 0.3,
        color = "black",
        alpha = 0.7
      ) +
      
      # CLD labels
      ggplot2::geom_text(
        data = label_data, 
        ggplot2::aes(x = !!rlang::sym(X_var), y = label_y_position, label = letters),
        size = 4, 
        color = "black", 
        fontface = "bold",
        inherit.aes = FALSE
      ) +
      
      ggplot2::labs(
        x = X_var,
        y = Y_var,
        color = X_var
      ) +
      
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
      
      ggplot2::theme_classic() + 
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = ggplot2::element_text(size = 10, face = "bold"),
        axis.text.x = ggplot2::element_text(size = 9, color = "black", 
                                            angle = 45, hjust = 1, vjust = 1),
        axis.text.y = ggplot2::element_text(size = 9, color = "black"),
        legend.position = "none",
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )
    
    return(rich_plot)
  }

PlotRich(physeq=physeq_glo_rare,
         X_var="siteID",
         Y_var="hill_0",
         cld_df=alpha_siteID$letters)

# evenness as pielou J index
PlotRich(physeq=physeq_glo_rare,
         X_var="siteID",
         Y_var="pielou_j",
         cld_df=get_pairwise_letters(
           df = physeq_glo_rare@sam_data %>%
             as.matrix() %>%
             as.data.frame() %>%
             mutate(pielou_j = as.numeric(pielou_j)),
           formula = formula(pielou_j ~ siteID)
         )$letters
)





# **** GLMR approach **** ------------------------------------------------------
library(MASS)
library(DHARMa)
library(emmeans)
library(multcomp)
library(car)

# Fit Poisson GLM
model_pois <- glm(formula = hill_0 ~ siteID * fertStatus, 
                  data = data.frame(physeq_glo_rare@sam_data), 
                  family = poisson)

# Check for overdispersion
summary(model_pois)
# If residual deviance >> residual df, you have overdispersion
# For a Poisson regression model, overdispersion is indicated if the Residual 
# Deviance is significantly larger than its Degrees of Freedom.
dispersion_parameter = 1177.6 / 96
dispersion_parameter

# If the ratio is significantly greater than 1: Overdispersion is present.

model_qpois <- glm(formula = hill_0 ~ siteID * fertStatus, 
                   data = data.frame(physeq_glo_rare@sam_data), 
                   family = quasipoisson)

summary(model_qpois)

# Better: use DHARMa for diagnostic plots
simulateResiduals(model_pois, plot = TRUE)

# INTERPRETATION. 
# 1) QQ Plot Residuals (Left Panel)This plot checks if the model residuals follow the 
# expected theoretical distribution (uniform distribution, in the case of DHARMa-scaled
# residuals). Data Points: The observed residuals (y-axis) closely follow the 
# diagonal red line (Expected values, x-axis). This is what you want to see.
# * Kolmogorov-Smirnov (KS) Test:$p = 0.44192$ (n.s. - not significant). 
#   Since the p-value is much greater than your typical significance
#   level (e.g., alpha = 0.05$), we do not reject the null hypothesis that the 
#   residuals are uniformly distributed. This confirms that the residuals are 
#   well-behaved and follow the expected pattern.
# * Dispersion Test: p = 0.384$ (n.s.). This test specifically checks for over
#   underdispersion. Since the p-value is non-significant, the Negative Binomial 
#   model has successfully accounted for the dispersion in your data. 
# * Outlier Test: p = 1 (n.s.). The model does not show evidence of significant 
#   outliers. 
# 2) DHARMa Residual vs. Predicted (Right Panel). This plot checks for violations
# of homogeneity of variance and whether any predictor variables were missed 
# (i.e., patterns in the residuals across the range of predicted values). 
# * Horizontal Dashed Lines: The center dashed line is the mean residual 
#   (should be around 0.5), and the other two show the 0.25 and 0.75 quantiles.
# * Solid Black Line: This is a fitted line (often a smooth trend) through the 
#   median residuals. It is nearly flat, indicating no major mean shifts or trends
#   across the predicted values.
# * Solid Red Curves (Quantile Deviations): These represent the 0.25$ and 0.75 
#   quantiles of the residuals. They are nearly flat and parallel to the mean, 
#   which confirms no strong pattern in variance (homogeneity of variance is fine).
# * Combined adjusted quantile test n.s.: The formal test for the quantile 
#   deviations is not significant, providing statistical confirmation that there
#   are no systematic deviations in the residuals based on the predicted response.

# Summary of the DHARMa ResultsThe DHARMa plots confirm that the Negative Binomial model:
# - Has normally distributed residuals (QQ Plot).
# - Does not suffer from overdispersion (Dispersion test).
# - Has no significant outliers (Outlier test).
# - Has no systematic issues with mean or variance (Residuals vs. Predicted Plot).
# Bottom line. The model is good enough to proceed with confidence!

# If overdispersed, use negative binomial
model_nb <- MASS::glm.nb(formula = hill_0 ~ siteID * fertStatus, 
                         data = data.frame(physeq_glo_rare@sam_data)
                         )

simulateResiduals(model_nb, plot = TRUE)

# Analysis of Deviance (ANOVA-style test)  
Anova(model_nb, type = "II")

# Analysis of Deviance: This is the equivalent of an ANOVA for models based on 
# likelihood (like the Negative Binomial model). Instead of comparing sums of 
# squares, it compares the change in deviance (which relates to the log-likelihood)
# between nested models to test the effect of each term.
# * Type II Tests: This is a method of partitioning the variance (or deviance) 
#   that tests the effect of a term after accounting for all other terms in the model,
#   except for the terms that contain it.

# For example, the effect of siteID is tested after accounting for the main effect
# of fertStatus (but not the siteID:fertStatus interaction). Type II is generally 
# the standard and preferred method when there is an interaction term in the model,
# unless the interaction is significant (which it is not here).

# Post-hoc pairwise comparisons with compact letter display
emm <- emmeans(model_nb, ~ siteID * fertStatus) 
pairs(emm)

# Compact letter display
cld(emm, Letters = letters, adjust = "fdr")

# The emmeans function first calculated the Estimated Marginal Means (EMMs) for 
# every combination of your two factors (siteID and fertStatus). Then, pairs(emm)
# performed Tukey-adjusted pairwise comparisons between every single one of 
# those 10 factor level combinations (5 sites x 2 fertilizer statuses).

# 1. The Scale
# The most important note is at the bottom: Results are given on the log (not the
# response) scale. 
# * Estimate: The values in the estimate column are the estimated difference between
#   the two groups on the logarithm scale. Because the model used 
#   a log-link (as is standard for the Negative Binomial family), these estimates
#   represent the log Ratio of Means. For example, if the estimate is X, the mean
#   of the first group is e^X times the mean of the second group.2. 
# * The P-Value Adjustment. The p-values have been adjusted using the Tukey method.
#   You performed 45 separate comparisons (a "family" of tests). Running 
#   this many tests increases the Family-Wise Error Rate (FWER), meaning you're 
#   more likely to get a false positive (Type I error). The Tukey Adjustment  
#   method controls the FWER, ensuring that the probability of making at least
#   one Type I error across all 45 tests is kept at your chosen significance 
#   level (e.g., alpha=0.05). This is why most of the p$-values are quite high.

# The significant differences are almost exclusively comparisons involving RHN versus
# other sites (ESC and LUX), and a comparison involving two fertilized treatments 
# (RHN FERT vs. ESC UNFERT). Specifically, RHN UNFERT (the unfertilized RHN site) 
# appears to have a higher hill_0 diversity than many other combinations.

# SUGGESTIONS! Simplify Interpretation. The Anova() result showed your interaction 
# term was not significant (p=0.701) and your fertStatus main effect was not 
# significant (p=0.281). Since the interaction and fertStatus are non-significant, 
# the most meaningful way to interpret the highly significant siteID effect is to
# average over the fertStatus factor.

emmeans(model_nb, ~ siteID) %>% pairs(adjust = "tukey")

# This will give you 10 comparisons (instead of 45), which will have more statistical
# power and be easier to interpret, isolating the effect you know is significant. 
# Would you like to run and interpret that simplified test?

# The output includes a note: NOTE: Results may be misleading due to involvement 
# in interactions. While the Anova() test suggested the interaction (siteID:fertStatus)
# was non-significant (p=0.701), this note is a standard cautionary flag from emmeans. 
# Given your large interaction p-value, averaging over fertStatus is statistically 
# sound and the results below are the best way to describe the main effect of siteID.

# The table compares the Estimated Marginal Means (EMMs) of your five sites, 
# averaged over both fertilized and unfertilized treatments. The results are 
# on the log scale, meaning the estimate is the log ratio of the means. 
# We look for p-values <= 0.05 (your significance threshold) to find which sites 
# are significantly different from each other.

# To summarize the patterns based on the estimates:
# * Lowest Diversity (ESC): The ESC site is consistently on the lower end, being
# significantly less diverse than three other sites (HAN, LC, and RHN). 
# * Highest Diversity (RHN): The RHN site is consistently on the higher end. It 
# is significantly more diverse than the ESC and LUX sites, and while its 
# difference from HAN and LC is large (negative estimates of -0.12$ and -0.126),
# the difference is not quite statistically significant after the Tukey adjustment.
# * Intermediate Sites (HAN, LC, LUX): These sites sit in the middle, generally 
# not differing significantly from each other.

# Get estimated marginal means
emm <- emmeans(model_nb, ~  siteID * fertStatus)
emm_df <- as.data.frame(emm)

# Get compact letter display
cld_results <- cld(emm, Letters = letters, adjust = "fdr")
cld_df <- as.data.frame(cld_results)
cld_df




df_model_nb <- as.data.frame(model.frame(model_nb))  
df_model_nb

df_pred$pred <- predict(model_nb, type = "response")  

# 2. Points with confidence intervals (better shows raw data):
ggplot(df_model_nb, aes(x = siteID, y = hill_0, color = fertStatus)) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  geom_pointrange(data = cld_df,
                  aes(x = siteID, y = emmean, 
                      ymin = asymp.LCL, ymax = asymp.UCL),
                  position = position_dodge(0.3),
                  size = 1, linewidth = 1.2) +
  geom_text(data = cld_df,
            aes(label = .group, y = asymp.UCL + 2),
            position = position_dodge(0.3), 
            color = "black", size = 4) +
  labs(y = "Species Richness", x = "Site") +
  theme_classic()
  
# 3. Interaction plot (shows interaction effect clearly):
  r# Interaction plot
emmip(model_nb, fertStatus ~ siteID) +
  theme_classic() +
  labs(y = "Predicted Species Richness", x = "Site",
       color = "Fertilization")

# Or with ggplot
ggplot(cld_df, aes(x = siteID, y = emmean, 
                   color = fertilization, group = fertilization)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  labs(y = "Species Richness", x = "Site") +
  theme_classic() +
  theme(legend.position = "top")

# 4. Boxplot with model predictions overlaid:
  rggplot(df, aes(x = site, y = richness, fill = fertilization)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = fertilization), 
              position = position_jitterdodge(jitter.width = 0.2),
              alpha = 0.6) +
  geom_point(data = cld_df, aes(y = emmean),
             position = position_dodge(0.75),
             size = 4, shape = 18) +
  geom_text(data = cld_df,
            aes(label = .group, y = upper.CL + 2),
            position = position_dodge(0.75)) +
  labs(y = "Species Richness", x = "Site") +
  theme_classic()
  
# 5. Faceted plot (if interaction is significant):
  rggplot(df, aes(x = site, y = richness)) +
  geom_boxplot(aes(fill = site), alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  geom_pointrange(data = cld_df,
                  aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
                  color = "red", size = 1) +
  geom_text(data = cld_df,
            aes(label = .group, y = upper.CL + 2)) +
  facet_wrap(~ fertilization, ncol = 2) +
  labs(y = "Species Richness", x = "Site") +
  theme_bw() +
  theme(legend.position = "none")




# Plot
ggplot(cld_df, aes(x = siteID, y = emmean, fill = fertStatus)) +
    geom_col(position = position_dodge(0.9), color = "black") +
    geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL ),
                  position = position_dodge(0.9), width = 0.25) +
    geom_text(aes(label = .group, y = asymp.UCL  + 1),
              position = position_dodge(0.9), size = 4) +
    labs(y = "Species Richness", x = "Site",
         fill = "Fertilization") +
    theme_classic() +
    theme(legend.position = "top")





































































# Test the plotting function
PlotRich(
  physeq_ITS_rare_root, "Treatment", "hill_0",
  physeq_ITS_rare_root@sam_data %>%
    as.matrix() %>%
    as.data.frame() %>%
    mutate(hill_0 = as.numeric(hill_0)) %>%
    CompSampl(., formula(hill_0 ~ Treatment), comparisons = PairCalc(4)) %>%
    pull(Letters),
  200
) +
  scale_color_manual(values = c("#825121", "#CC2D35", "#FFB400", "#00A6ED"))

# Using wilcox.test
CompSampl <- function(df, formula, comparisons) {
  require(multcompView)
  require(lazyeval)
  
  test_CC <-
    compare_means(formula, data = df, method = "wilcox.test", p.adjust.method = "none")
  
  test_CC$adj.pval <-
    p.adjust(test_CC$p, method = "BH", n = comparisons)
  # print(test_CC)
  test_CC <-
    as.data.frame(test_CC)[, c(2, 3, 5)] # to change form p to p.adj do 4 to 5
  test_CC2 <-
    data.frame(test_CC[, 2], test_CC[, 1], test_CC[, 3])
  colnames(test_CC2) <-
    c("group1", "group2", "p.adj") # change p to p.adj
  rbind(test_CC, test_CC2) -> test_all
  print(test_all)
  
  dist_CC <- as.dist(xtabs(test_all[, 3] ~ (test_all[, 2] + test_all[, 1])), diag = TRUE)
  # print(dist_CC)
  res_CC <-
    data.frame(multcompLetters(dist_CC)["Letters"])
  return(res_CC)
}


# Calculate number of multiple pairwise comparisons based on the number of levels
PairCalc <- function(n) {
  n_comparison <- n * (n - 1) / 2
  return(n_comparison)
}

PairCalc(6)

# Extracting multiple comparisons and p.adj
physeq_ITS_rare_root@sam_data %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(hill_0 = as.numeric(hill_0)) %>%
  as_tibble() %>%
  CompSampl(., formula(hill_0 ~ Treatment), comparisons = PairCalc(4))

# plotting Richness function ---------------------------------------------------
PlotRich <- function(physeq, X_var, Y_var, my_labels, labels_y) {
  require(phyloseq)
  require(tidyverse)
  
  # extract dataframe
  dataframe <-
    physeq@sam_data %>%
    as.matrix() %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(!!Y_var := as.numeric(!!sym(Y_var)))
  
  # print factor order
  
  # Calculate labels_y based on the maximum value of Y_var
  labels_y <- max(dataframe[, Y_var]) + 0.1 * max(dataframe[, Y_var])
  
  # plot
  rich_plot <-
    ggplot(dataframe, aes(x = get(X_var), y = !!sym(Y_var))) +
    geom_jitter(
      position = position_jitter(0.4), size = 1, shape = 16,
      aes(color = get(X_var))
    ) +
    # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
    #             geom="pointrange", color="black", shape=18, size=0.8) +
    stat_summary(
      fun = median, fun.min = median, fun.max = median,
      geom = "crossbar", width = 0.8, linewidth = 0.3, color = "black"
    ) +
    # stat_summary(
    #  geom="pointrange",
    #  fun.min = function(z) { quantile(z, 0.25) },
    #  fun.max = function(z) { quantile(z, 0.75) },
    #  fun = median, color="black", shape=18, size=1,
    #  show.legend = FALSE) +
    stat_summary(
      geom = "text", angle = 0, label = my_labels,
      fun = max, aes(y = labels_y), size = 3, color = "black"
    ) +
    # scale_color_manual(values=c("grey", "blue")) +
    expand_limits(y = 0) +
    theme_bw() +
    theme(
      plot.title = element_markdown(size = 8, face = "bold", vjust = 0.5, hjust = 0.5),
      axis.title = element_markdown(size = 7),
      axis.text.x = element_markdown(size = 7, colour = "black", angle = 33, hjust = 1, vjust = 1),
      axis.text.y = element_markdown(size = 7, angle = 0, hjust = 0.5),
      legend.position = "none"
    )
  
  return(rich_plot)
}

# Test the plotting function
PlotRich(
  physeq_ITS_rare_root, "Treatment", "hill_0",
  physeq_ITS_rare_root@sam_data %>%
    as.matrix() %>%
    as.data.frame() %>%
    mutate(hill_0 = as.numeric(hill_0)) %>%
    CompSampl(., formula(hill_0 ~ Treatment), comparisons = PairCalc(4)) %>%
    pull(Letters),
  200
) +
  scale_color_manual(values = c("#825121", "#CC2D35", "#FFB400", "#00A6ED"))







# ********************************************************************----------
# ********************************************************************----------
# ********************************************************************----------







# Importing metadata -----------------------------------------------------------
metadata_99 <-
  read.csv("data_final_99/metadata_pacbio.csv", header = TRUE) %>% 
  column_to_rownames("SampleID")

head(metadata_99)
dim(metadata_99)

# write the metadata
write.csv(metadata_99, file="metadata_MMNPRNT_pacbio.csv")

# Importing metadata
metadata_HMV <- read.csv("metadata_kbs/metadata_pacbio_HMV.csv")
head(metadata_HMV)
dim(metadata_HMV)

metadata_merge <- 
  left_join(metadata_99 %>% rownames_to_column("SampleID"), 
            metadata_HMV,
            by = "SampleID")

head(metadata_merge)
write_csv(metadata_merge, "metadata_merge_amf.csv")

# Importing sequences ----------------------------------------------------------
Zotu_99 <- 
  readDNAStringSet("data_final_99/asv_filt_nico_glo_99.fasta", format="fasta", 
                   seek.first.rec=TRUE, use.names=TRUE)
Zotu_99

# Importing zotu_table ---------------------------------------------------------
Zotutable_99 <-
  read.delim("data_final_99/zotu_table_nico_glo_99.txt", header = TRUE, sep = "\t") %>% 
  rename("Zotu" = 1) %>% 
  mutate(Zotu = gsub("_", "", Zotu)) %>% 
  column_to_rownames("Zotu")

head(Zotutable_99)
dim(Zotutable_99)


# Importing RAX ML tree --------------------------------------------------------
# Replace the file path with the actual location of your RAxML tree file
raxml_tree <- read.tree("RAXML_results/RAxML_bestTree.result")
raxml_tree <- read.tree("RAXML_results/RAxML_bipartitions.result")
raxml_tree

# 1. Extract the tip labels from the tree
tip_labels <- raxml_tree$tip.label
zotu_labels <- grep("^Zotu", tip_labels, value = TRUE)
raxml_tree_filtered <- drop.tip(raxml_tree, setdiff(tip_labels, zotu_labels))

#Relabel the Zotu tips to keep only the "Zotu" and the number (remove everything after "_")
raxml_tree_filtered$tip.label <- gsub("^(Zotu[0-9]+)_.*", "\\1", raxml_tree_filtered$tip.label)
raxml_tree_filtered

raxml_tree_filtered$tip.label

# *********************************************---------------------------------
# GENERATE PHYLOSEQ OBJECT -----------------------------------------------------

colnames(Zotutable_99)
colnames(Zotutable_99)
rownames(metadata_99)

dim(Zotutable_99)
dim(metadata_99)
dim(blastax_99)


physeq_glo99 <-
  phyloseq(
    otu_table(Zotutable_99, taxa_are_rows = TRUE),
    sample_data(metadata_99),
    tax_table(as.matrix(blastax_99_fix)),
    Zotu_99,
    raxml_tree_filtered) %>% 
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.)

physeq_glo99

















# COLOR PALETTE *******************************---------------------------------

palette_v2 <- c("#825121","#CC2D35","#FFB400","#00A6ED",
                "#7FB800","#2D3142", "#979797","#F0E442")

palette_v2

palette_v3 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E11", "#CC79A7")


# ***********************************************************************-------
# Import datasets for beta-diversity cutoff -------------------------------------

# OTU/ASV tables
otutable_90 <-read.delim(file.path(data_path, "data_final_99/ASV_clusters/otutab_90.txt"), 
                         row.names = 1)

otutable_91 <-read.delim(file.path(data_path, "data_final_99/ASV_clusters/otutab_91.txt"), 
                         row.names = 1)

otutable_92 <-read.delim(file.path(data_path, "data_final_99/ASV_clusters/otutab_92.txt"), 
                         row.names = 1)

otutable_93 <-read.delim(file.path(data_path, "data_final_99/ASV_clusters/otutab_93.txt"), 
                         row.names = 1)

otutable_94 <-read.delim(file.path(data_path, "data_final_99/ASV_clusters/otutab_94.txt"), 
                         row.names = 1)

otutable_95 <-read.delim(file.path(data_path, "data_final_99/ASV_clusters/otutab_95.txt"), 
                         row.names = 1)

otutable_96 <-read.delim(file.path(data_path, "data_final_99/ASV_clusters/otutab_96.txt"), 
                         row.names = 1)

otutable_97 <-read.delim(file.path(data_path, "data_final_99/ASV_clusters/otutab_97.txt"), 
                         row.names = 1)

otutable_98 <-read.delim(file.path(data_path, "data_final_99/ASV_clusters/otutab_98.txt"), 
                         row.names = 1)

otutable_99 <-read.delim(file.path(data_path, "data_final_99/ASV_clusters/otutab_99.txt"), 
                         row.names = 1)

otutable_100 <-read.delim(file.path(data_path, "data_final_99/ASV_clusters/otutab_00.txt"), 
                         row.names = 1)


# Taxonomy 
extract_blasTAX <- function(tax_path,
                            namemap_path){
  
  taxonomy <- read.delim(tax_path, header = TRUE, strip.white = TRUE, sep="\t") 
  name_mapping <- read.delim(namemap_path, header = FALSE, strip.white = TRUE, sep="\t") 

  if (length(setdiff(name_mapping$V2, taxonomy$OTU_ID)) > 0) {
    cli::cli_alert_warning(
      paste("No significant similarity found for:", 
            paste(setdiff(name_mapping$V2, taxonomy$OTU_ID), collapse = ", "))
    )
  }

  taxonomy_clean <-
    taxonomy %>% 
    rename(Query = 1, Query_Score = 2) %>% 
    mutate_all(~ str_trim(., side = "both")) %>% # remove trailing spaces in both sides
    as.matrix() %>% 
    as.data.frame() %>% 
    dplyr::arrange(Query)
  
  name_mapping_clean <-
    name_mapping %>% 
    rename(Zotu = V1, Query = V2) %>% 
    mutate(Zotu = gsub(">", "", Zotu)) %>% 
    mutate_all(~ str_trim(., side = "both")) %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    dplyr::arrange(Query)
  
  blastax <-
    full_join(name_mapping_clean, taxonomy_clean, by="Query")
  
  return(blastax)
  
}



taxonomy_90 <- extract_blasTAX(
  tax_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_90_taxonomy/taxonomy_blast_final.txt"
  ),
  namemap_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_90_taxonomy/name_mapping.txt"
  )
)

taxonomy_90


taxonomy_91 <- extract_blasTAX(
  tax_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_91_taxonomy/taxonomy_blast_final.txt"
  ),
  namemap_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_91_taxonomy/name_mapping.txt"
  )
)

taxonomy_91

taxonomy_92 <- extract_blasTAX(
  tax_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_92_taxonomy/taxonomy_blast_final.txt"
  ),
  namemap_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_92_taxonomy/name_mapping.txt"
  )
)

taxonomy_92


taxonomy_93 <- extract_blasTAX(
  tax_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_93_taxonomy/taxonomy_blast_final.txt"
  ),
  namemap_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_93_taxonomy/name_mapping.txt"
  )
)

taxonomy_93


taxonomy_94 <- extract_blasTAX(
  tax_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_94_taxonomy/taxonomy_blast_final.txt"
  ),
  namemap_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_94_taxonomy/name_mapping.txt"
  )
)

taxonomy_94


taxonomy_95 <- extract_blasTAX(
  tax_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_95_taxonomy/taxonomy_blast_final.txt"
  ),
  namemap_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_95_taxonomy/name_mapping.txt"
  )
)

taxonomy_95


taxonomy_96 <- extract_blasTAX(
  tax_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_96_taxonomy/taxonomy_blast_final.txt"
  ),
  namemap_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_96_taxonomy/name_mapping.txt"
  )
)

taxonomy_96


taxonomy_97 <- extract_blasTAX(
  tax_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_97_taxonomy/taxonomy_blast_final.txt"
  ),
  namemap_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_97_taxonomy/name_mapping.txt"
  )
)

taxonomy_97


taxonomy_98 <- extract_blasTAX(
  tax_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_98_taxonomy/taxonomy_blast_final.txt"
  ),
  namemap_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_98_taxonomy/name_mapping.txt"
  )
)

taxonomy_98


taxonomy_99 <- extract_blasTAX(
  tax_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_99_taxonomy/taxonomy_blast_final.txt"
  ),
  namemap_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_99_taxonomy/name_mapping.txt"
  )
)

taxonomy_99


taxonomy_100 <- extract_blasTAX(
  tax_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_100_taxonomy/taxonomy_blast_final.txt"
  ),
  namemap_path = file.path(
    data_path,
    "data_final_99/ASV_clusters/otus_100_taxonomy/name_mapping.txt"
  )
)

taxonomy_100


# sequences 
Zotu_90 <- readDNAStringSet(file.path(data_path,"data_final_99/ASV_clusters/otus_90.fasta"), format="fasta", 
                   seek.first.rec=TRUE, use.names=TRUE)
Zotu_90

Zotu_91 <- readDNAStringSet(file.path(data_path,"data_final_99/ASV_clusters/otus_91.fasta"), format="fasta", 
                            seek.first.rec=TRUE, use.names=TRUE)
Zotu_91

Zotu_92 <- readDNAStringSet(file.path(data_path,"data_final_99/ASV_clusters/otus_92.fasta"), format="fasta", 
                            seek.first.rec=TRUE, use.names=TRUE)
Zotu_92

Zotu_93 <- readDNAStringSet(file.path(data_path,"data_final_99/ASV_clusters/otus_93.fasta"), format="fasta", 
                            seek.first.rec=TRUE, use.names=TRUE)
Zotu_93

Zotu_94 <- readDNAStringSet(file.path(data_path,"data_final_99/ASV_clusters/otus_94.fasta"), format="fasta", 
                            seek.first.rec=TRUE, use.names=TRUE)
Zotu_94

Zotu_95 <- readDNAStringSet(file.path(data_path,"data_final_99/ASV_clusters/otus_95.fasta"), format="fasta", 
                            seek.first.rec=TRUE, use.names=TRUE)
Zotu_95

Zotu_96 <- readDNAStringSet(file.path(data_path,"data_final_99/ASV_clusters/otus_96.fasta"), format="fasta", 
                            seek.first.rec=TRUE, use.names=TRUE)
Zotu_96

Zotu_97 <- readDNAStringSet(file.path(data_path,"data_final_99/ASV_clusters/otus_97.fasta"), format="fasta", 
                            seek.first.rec=TRUE, use.names=TRUE)
Zotu_97

Zotu_98 <- readDNAStringSet(file.path(data_path,"data_final_99/ASV_clusters/otus_98.fasta"), format="fasta", 
                            seek.first.rec=TRUE, use.names=TRUE)
Zotu_98

Zotu_99 <- readDNAStringSet(file.path(data_path,"data_final_99/ASV_clusters/otus_99.fasta"), format="fasta", 
                            seek.first.rec=TRUE, use.names=TRUE)
Zotu_99

Zotu_100 <- readDNAStringSet(file.path(data_path,"data_final_99/ASV_clusters/asv.fasta"), format="fasta", 
                            seek.first.rec=TRUE, use.names=TRUE)
Zotu_100


# Metadata ---------------------------------------------------------------------
metadata_99 <-
  read.csv(file=file.path(data_path, "data_final_99/metadata_pacbio.csv"), header = TRUE) %>% 
  column_to_rownames("SampleID") %>% 
  janitor::clean_names()

head(metadata_99)
dim(metadata_99)

# create phyloseq objects ------------------------------------------------------

dim(otutable_90);
dim(metadata_99)
dim(taxonomy_90)
dim(Zotu_90)


generate_phyloseq <- function(otu, metadata, taxonomy, sequences) {
  # convert to data.frame
  otu       <- as.data.frame(otu, stringsAsFactors = FALSE)
  taxonomy  <- as.data.frame(taxonomy, stringsAsFactors = FALSE) %>% 
    column_to_rownames("Zotu")
  metadata  <- as.data.frame(metadata, stringsAsFactors = FALSE)
  
  # check alignment
  taxa_missing <- setdiff(rownames(otu), rownames(taxonomy))
  if (length(taxa_missing)) {
    cli::cli_alert_warning("Some OTUs missing in taxonomy: {length(taxa_missing)}")
  }
  sample_missing <- setdiff(colnames(otu), rownames(metadata))
  if (length(sample_missing)) {
    cli::cli_alert_warning("Some samples missing in metadata: {length(sample_missing)}")
  }
  
  # build phyloseq
  ps <- phyloseq(
    otu_table(otu, taxa_are_rows = TRUE),
    sample_data(metadata),
    tax_table(as.matrix(taxonomy)),
    refseq(sequences)
  ) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    prune_samples(sample_sums(.) > 0, .)
  
  cli::cli_alert_success("phyloseq object created with {ntaxa(ps)} taxa and {nsamples(ps)} samples")
  
  ps
}


physeq_90 <- generate_phyloseq(otu=otutable_90,metadata=metadata_99,
                               taxonomy=taxonomy_90,sequences=Zotu_90)
physeq_90
physeq_91 <- generate_phyloseq(otu=otutable_91,metadata=metadata_99,
                               taxonomy=taxonomy_91,sequences=Zotu_91)
physeq_91
physeq_92 <- generate_phyloseq(otu=otutable_92,metadata=metadata_99,
                               taxonomy=taxonomy_92,sequences=Zotu_92)
physeq_92
physeq_93 <- generate_phyloseq(otu=otutable_93,metadata=metadata_99,
                               taxonomy=taxonomy_93,sequences=Zotu_93)
physeq_93
physeq_94 <- generate_phyloseq(otu=otutable_94,metadata=metadata_99,
                               taxonomy=taxonomy_94,sequences=Zotu_94)
physeq_94
physeq_95 <- generate_phyloseq(otu=otutable_95,metadata=metadata_99,
                               taxonomy=taxonomy_95,sequences=Zotu_95)
physeq_95
physeq_96 <- generate_phyloseq(otu=otutable_96,metadata=metadata_99,
                               taxonomy=taxonomy_96,sequences=Zotu_96)
physeq_96
physeq_97 <- generate_phyloseq(otu=otutable_97,metadata=metadata_99,
                               taxonomy=taxonomy_97,sequences=Zotu_97)
physeq_97
physeq_98 <- generate_phyloseq(otu=otutable_98,metadata=metadata_99,
                               taxonomy=taxonomy_98,sequences=Zotu_98)
physeq_98
physeq_99 <- generate_phyloseq(otu=otutable_99,metadata=metadata_99,
                               taxonomy=taxonomy_99,sequences=Zotu_99)
physeq_99
physeq_100 <- generate_phyloseq(otu=otutable_100,metadata=metadata_99,
                               taxonomy=taxonomy_100,sequences=Zotu_100)
physeq_100

# rarefaction 
pak::pak("germs-lab/BRCore")
pak::pak("germs-lab/BRCore@22d4b66")
library(BRCore)

sort(sample_sums(physeq_90))[10]   
hist(sort(sample_sums(physeq_90)))

physeq_90_rare <- multi_rarefy(
  physeq = physeq_90,
  depth_level = 600,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_90)

otu_table(physeq_90)
otu_table(physeq_90_rare)

sort(sample_sums(physeq_91))[10]  
physeq_91_rare <- multi_rarefy(
  physeq = physeq_91,
  depth_level = 1400,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_91)

sort(sample_sums(physeq_92))[10]  
physeq_92_rare <- multi_rarefy(
  physeq = physeq_92,
  depth_level = 3380,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_92)

sort(sample_sums(physeq_93))[10]  
physeq_93_rare <- multi_rarefy(
  physeq = physeq_93,
  depth_level = 2000,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_93)

sort(sample_sums(physeq_94))[10]  
physeq_94_rare <- multi_rarefy(
  physeq = physeq_94,
  depth_level = 4500,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_94)

sort(sample_sums(physeq_95))[10]  
physeq_95_rare <- multi_rarefy(
  physeq = physeq_95,
  depth_level = 5500,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_95)

sort(sample_sums(physeq_96))[10]  
physeq_96_rare <- multi_rarefy(
  physeq = physeq_96,
  depth_level = 8570,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_96)

sort(sample_sums(physeq_97))[10]  
physeq_97_rare <- multi_rarefy(
  physeq = physeq_97,
  depth_level = 10000,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_97)

sort(sample_sums(physeq_98))[10]  
physeq_98_rare <- multi_rarefy(
  physeq = physeq_98,
  depth_level = 10000,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_98)

sort(sample_sums(physeq_99))[10]  
physeq_99_rare <- multi_rarefy(
  physeq = physeq_99,
  depth_level = 10500,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_99)

sort(sample_sums(physeq_100))[10]  
physeq_100_rare <- multi_rarefy(
  physeq = physeq_100,
  depth_level = 10500,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_100)


# I think I should rarefy at the same depth

# plotting ordinations 

plotPCOA <- function(physeq, Col=NULL, She=NULL){
  
  taxa_num <- ncol(otu_table(physeq))
  print(taxa_num)
  
  pcoa <- phyloseq::ordinate(physeq, method = "PCoA", distance = "bray") 
  
  plot_res <-
    plot_ordination(physeq = physeq, 
                    ordination=pcoa,  
                    type="samples", 
                    color=Col, 
                    shape=She) +
    theme_bw() +
    theme(plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
          strip.text = element_markdown(size = 10, face = "bold"),
          axis.text.x = element_markdown(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
          axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
          strip.background = element_blank(),
          legend.key.height = unit(0.4, "cm"), 
          legend.key.width = unit(0.4, "cm"),
          legend.position = "right", 
          legend.title.align = 0.5,
          legend.text = element_markdown(size = 8)) +
    guides(color = guide_legend(nrow=1, override.aes = list(shape = 15, size = 3.5)),
           shape = guide_legend(nrow=1, override.aes = list(color = "black", size=2.5)),
           size = guide_legend(ncol=1))
  
  return(plot_res)
  
}

sample_data(physeq_100_rare)

ggarrange(
plotPCOA(physeq = physeq_100_rare, Col = "site_id", She = "fert_status"),
plotPCOA(physeq = physeq_99_rare, Col = "site_id", She = "fert_status"),
plotPCOA(physeq = physeq_98_rare, Col = "site_id", She = "fert_status"),
plotPCOA(physeq = physeq_97_rare, Col = "site_id", She = "fert_status"),
plotPCOA(physeq = physeq_96_rare, Col = "site_id", She = "fert_status"),
plotPCOA(physeq = physeq_95_rare, Col = "site_id", She = "fert_status"),
plotPCOA(physeq = physeq_94_rare, Col = "site_id", She = "fert_status"),
plotPCOA(physeq = physeq_93_rare, Col = "site_id", She = "fert_status"),
plotPCOA(physeq = physeq_92_rare, Col = "site_id", She = "fert_status"),
plotPCOA(physeq = physeq_91_rare, Col = "site_id", She = "fert_status"),
plotPCOA(physeq = physeq_90_rare, Col = "site_id", She = "fert_status"),
ncol = 3, nrow = 4,
common.legend = TRUE,
legend = "bottom")

# rerun at same depth of 600 reads per sample

physeq_90_rare_s <- multi_rarefy(
  physeq = physeq_90,
  depth_level = 600,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_90)

physeq_91_rare_s <- multi_rarefy(
  physeq = physeq_91,
  depth_level = 600,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_91)

physeq_92_rare_s <- multi_rarefy(
  physeq = physeq_92,
  depth_level = 600,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_92)

physeq_93_rare_s <- multi_rarefy(
  physeq = physeq_93,
  depth_level = 600,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_93)

physeq_94_rare_s <- multi_rarefy(
  physeq = physeq_94,
  depth_level = 600,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_94)

physeq_95_rare_s <- multi_rarefy(
  physeq = physeq_95,
  depth_level = 600,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_95)

physeq_96_rare_s <- multi_rarefy(
  physeq = physeq_96,
  depth_level = 600,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_96)

physeq_97_rare_s <- multi_rarefy(
  physeq = physeq_97,
  depth_level = 600,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_97)

physeq_98_rare_s <- multi_rarefy(
  physeq = physeq_98,
  depth_level = 600,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_98)

physeq_99_rare_s <- multi_rarefy(
  physeq = physeq_99,
  depth_level = 600,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_99)

physeq_100_rare_s <- multi_rarefy(
  physeq = physeq_100,
  depth_level = 600,
  num_iter = 10,
  threads = 8,
  set_seed = 10
) %>%
  do_phyloseq(otu_rare = ., physeq =  physeq_100)

otu_table(physeq_100_rare_s)

# plotting

plot_beta_percent <-
  ggarrange(
plotPCOA(physeq = physeq_100_rare_s, Col = "site_id", She = "fert_status") +
  labs(title="PCoA 100% read similarity"),
plotPCOA(physeq = physeq_99_rare_s, Col = "site_id", She = "fert_status")+
  labs(title="PCoA 99% read similarity"),
plotPCOA(physeq = physeq_98_rare_s, Col = "site_id", She = "fert_status")+
  labs(title="PCoA 98% read similarity"),
plotPCOA(physeq = physeq_97_rare_s, Col = "site_id", She = "fert_status")+
  labs(title="PCoA 97% read similarity"),
plotPCOA(physeq = physeq_96_rare_s, Col = "site_id", She = "fert_status")+
  labs(title="PCoA 96% read similarity"),
plotPCOA(physeq = physeq_95_rare_s, Col = "site_id", She = "fert_status")+
  labs(title="PCoA 95% read similarity"),
plotPCOA(physeq = physeq_94_rare_s, Col = "site_id", She = "fert_status")+
  labs(title="PCoA 94% read similarity"),
plotPCOA(physeq = physeq_93_rare_s, Col = "site_id", She = "fert_status")+
  labs(title="PCoA 93% read similarity"),
plotPCOA(physeq = physeq_92_rare_s, Col = "site_id", She = "fert_status")+
  labs(title="PCoA 92% read similarity"),
plotPCOA(physeq = physeq_91_rare_s, Col = "site_id", She = "fert_status")+
  labs(title="PCoA 91% read similarity"),
plotPCOA(physeq = physeq_90_rare_s, Col = "site_id", She = "fert_status")+
  labs(title="PCoA 90% read similarity"),
ncol=3,
nrow=4,
common.legend=TRUE,
legend="bottom")

plot_beta_percent


plot_beta_percent <-
  ggarrange(
    plotPCOA(physeq = physeq_100_rare_s, Col = "site_id", She = "fert_status") +
      labs(title="100%"),
    plotPCOA(physeq = physeq_99_rare_s, Col = "site_id", She = "fert_status")+
      labs(title="99%"),
    plotPCOA(physeq = physeq_98_rare_s, Col = "site_id", She = "fert_status")+
      labs(title="98%"),
    plotPCOA(physeq = physeq_97_rare_s, Col = "site_id", She = "fert_status")+
      labs(title="97%"),
    plotPCOA(physeq = physeq_96_rare_s, Col = "site_id", She = "fert_status")+
      labs(title="96%"),
    plotPCOA(physeq = physeq_95_rare_s, Col = "site_id", She = "fert_status")+
      labs(title="95%"),
    plotPCOA(physeq = physeq_94_rare_s, Col = "site_id", She = "fert_status")+
      labs(title="94%"),
    plotPCOA(physeq = physeq_93_rare_s, Col = "site_id", She = "fert_status")+
      labs(title="93%"),
    plotPCOA(physeq = physeq_92_rare_s, Col = "site_id", She = "fert_status")+
      labs(title="92%"),
    plotPCOA(physeq = physeq_91_rare_s, Col = "site_id", She = "fert_status")+
      labs(title="91%"),
    plotPCOA(physeq = physeq_90_rare_s, Col = "site_id", She = "fert_status")+
      labs(title="90%"),
    ncol=11,
    nrow=1,
    common.legend=TRUE,
    legend="bottom")

plot_beta_percent


# INTERPRETATION
#Looking at your figure, you're observing that as you move from 100% similarity (ASVs) down to 90% similarity (more aggressive clustering), the variance explained by the first two PCoA axes increases substantially - from about 12.3% total at 100% to 35.7% at 90%.
#This is a common and expected pattern in microbiome ordination analyses. Here's why:
#  The underlying mechanism
#Noise reduction through aggregation: At 100% similarity (ASVs), you're working with the finest possible resolution. This means:

#You have many OTUs, including rare variants that may differ by just 1-2 nucleotides
#These fine-scale differences often represent sequencing errors, intraspecific variation, or ecologically similar strains
#This creates "noise" that obscures the major ecological patterns

#As you cluster more aggressively (97% → 90%), you're:
  
#  Collapsing ecologically similar taxa into single units
#Averaging out stochastic variation and technical noise
#Reducing the dimensionality of your data while retaining biologically meaningful signal

#Concentration of signal: When related sequences are clustered together, the ecological patterns that distinguish your samples (fertilized vs. unfertilized, site differences) become more prominent relative to the background variation. The major axes of variation can now capture more of the total variance because the signal-to-noise ratio has improved.
#What this means for your data
#The clear separation between FERT and UNFERT samples is visible across all clustering levels, but becomes progressively clearer (and the site-level clustering becomes more apparent) as you increase clustering. This suggests the fertilization effect and site differences operate at a broader taxonomic level than single nucleotide variants.
#Practical consideration
#While 90% shows the highest variance explained, 97% similarity is typically recommended for ITS data as it better balances biological resolution with noise reduction and aligns with species-level hypotheses in fungal ecology.


#Here are several statistical approaches to formally test whether the signal strengthens at broader taxonomic levels:
#  1. PERMANOVA with effect size comparison (Most direct)
#Run PERMANOVA at each clustering level and compare the R² values:
#  r# For each clustering level (100%, 99%, ..., 90%)
adonis2(distance_matrix ~ fert_status * site_id, 
        data = metadata, 
        permutations = 999)
#What to look for:
# If R² increases as you move from 100% → 90%, the explanatory variables account for more variation at broader taxonomic levels
#Compare both main effects (fert_status, site_id) and their interaction

#2. ANOSIM R-statistic across clustering levels
#r# For each clustering level
anosim(distance_matrix, grouping = metadata$fert_status, permutations = 999)
#What to look for:
  #ANOSIM R ranges from -1 to 1 (R ≈ 1 means strong separation)
#Plot R-statistic vs. clustering threshold
#If R increases with clustering, group separation strengthens

#3. Multivariate dispersion analysis (PERMDISP)
#Check if the increased variance explained is due to tighter within-group clustering:
#  r# For each clustering level
dist_matrix <- vegdist(otu_table, method = "bray")
mod <- betadisper(dist_matrix, metadata$fert_status)
permutest(mod)
#What to look for:
  
#  If within-group dispersion decreases with clustering, groups become more cohesive
#This supports your noise-reduction hypothesis

#4. Mantel test correlation
#Test whether distance matrices at different clustering levels correlate differently with environmental variables:
#  r# For each clustering level
mantel(distance_matrix, env_distance_matrix, method = "pearson")
#What to look for:
# Increasing Mantel r with broader clustering suggests environmental factors better predict community dissimilarity

#5. Variance partitioning
#Decompose variance to see how much is explained by fert_status vs site_id at each level:
#  r# Using vegan package
varpart(otu_table, ~fert_status, ~site_id, data = metadata)
#What to look for:
  
#  Track unique and shared variance across clustering levels
#Shows which factor's signal strengthens more with clustering

#6. Simple visualization approach
#Create a summary plot:
r# Extract R² from PERMANOVA and variance explained from PCoA
clustering_levels <- c(100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90)
results <- data.frame(
  clustering = clustering_levels,
  permanova_r2_fert = ...,
 
  permanova_r2_site = ...,
  pcoa_var_axis1 = ...,
  pcoa_var_axis2 = ...
)

# Plot trends
plot(clustering, permanova_r2_fert, type = "b")

#My recommendation
#Start with PERMANOVA R² comparison (#1) - it directly tests your hypothesis and is straightforward to interpret. Supplement with ANOSIM (#2) for a non-parametric confirmation. If you want to publish this observation, include the variance partitioning (#5) to show how the relative importance of factors changes with taxonomic resolution.
#The key result would be showing that R² (or other effect sizes) systematically increases as clustering becomes more aggressive, which statistically confirms what you're seeing visually in the PCoA plots.






 # ********************************************************************----------
# ********************************************************************----------
# ********************************************************************----------

taxonomy <-read.delim(file.path(data_path, "data_final_99/ASV_clusters/otus_90_taxonomy/taxonomy_blast_final.txt"), 
                      header = TRUE, strip.white = TRUE, sep="\t") 

dim(taxonomy)
head(taxonomy)
taxonomy[, 1:5]

taxonomy$OTU_ID

taxonomy %>% 
  rename(Query = 1, Query_Score = 2) %>% 
  mutate_all(~ str_trim(., side = "both")) %>% # remove trailing spaces in both sides
  as.matrix() %>% 
  as.data.frame() %>% 
  dplyr::arrange(Query)

name_mapping <-read.delim(file.path(data_path, "data_final_99/ASV_clusters/otus_90_taxonomy/name_mapping.txt"), 
                          header = FALSE, strip.white = TRUE, sep="\t") 

dim(name_mapping)
head(name_mapping)
name_mapping

name_mapping$V2

length(unique(taxonomy$OTU_ID)); length(unique(name_mapping$V2))

# Items present in taxonomy but NOT in name_mapping
setdiff(taxonomy$OTU_ID, name_mapping$V2)
setdiff(name_mapping$V2, taxonomy$OTU_ID)

"Query_92"  "Query_101"

identical(taxonomy$OTU_ID, name_mapping$V2)

name_mapping %>% 
  rename(Zotu = V1, Query = V2) %>% 
  mutate(Zotu = gsub(">", "", Zotu)) %>% 
  mutate_all(~ str_trim(., side = "both")) %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  dplyr::arrange(Query)



# functions
extract_blasTAX <- function(tax_path, namemap_path) {
  library(dplyr)
  library(stringr)
  # Read as character only; strip whitespace
  taxonomy <- read.delim(tax_path, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, colClasses = "character",
                         fileEncoding = "UTF-8")
  name_mapping <- read.delim(namemap_path, header = FALSE, sep = "\t",
                             stringsAsFactors = FALSE, colClasses = "character",
                             fileEncoding = "UTF-8")
  
  # Clean helpers
  strip_bom <- function(x) sub("^\ufeff", "", x)        # remove BOM if present
  trim_all  <- function(df) mutate(across(everything(), ~ str_trim(strip_bom(.x))))
  
  taxonomy_clean <- taxonomy %>%
    rename(Query = 1, Query_Score = 2) %>%
    trim_all() %>%
    arrange(Query)
  
  name_mapping_clean <- name_mapping %>%
    rename(Zotu = V1, Query = V2) %>%
    mutate(Zotu = gsub(">", "", Zotu, fixed = TRUE)) %>%
    trim_all() %>%
    arrange(Query)
  
  # Strict but type-safe comparison
  same_order <- identical(taxonomy_clean$Query, name_mapping_clean$Query)
  same_chars <- identical(as.character(taxonomy_clean$Query),
                          as.character(name_mapping_clean$Query))
  
  if (!same_order) {
    # Helpful diagnostics
    cli::cli_alert_warning("Query vectors are not identical (order/type/values).")
    # What differs?
    cli::cli_alert_info("Length: tax={.val {length(taxonomy_clean$Query)}}, map={.val {length(name_mapping_clean$Query)}}")
    cli::cli_alert_info("Type:   tax={.val {class(taxonomy_clean$Query)}}, map={.val {class(name_mapping_clean$Query)}}")
    # Show any set mismatches
    only_in_tax <- setdiff(taxonomy_clean$Query, name_mapping_clean$Query)
    only_in_map <- setdiff(name_mapping_clean$Query, taxonomy_clean$Query)
    if (length(only_in_tax)) cli::cli_alert_danger("Only in taxonomy: {head(only_in_tax, 5)}")
    if (length(only_in_map)) cli::cli_alert_danger("Only in name map: {head(only_in_map, 5)}")
    # If only types/order differ but values match:
    if (same_chars) cli::cli_alert_success("Values match after character coercion.")
  } else {
    cli::cli_alert_success("Query vectors are identical.")
  }
  
  blastax <- dplyr::full_join(name_mapping_clean, taxonomy_clean, by = "Query")
  return(blastax)
}













#************* BETA DIVERSITY  *************------------------------------------
# Custom PCoA plot 
plot_pcoa <- function(ord, phy, Col=NULL, She=NULL){
  
  require(phyloseq)
  
  plot_res <-
    plot_ordination(physeq = phy, 
                    ordination=ord,  
                    type="samples", 
                    color=Col, 
                    shape=She) +
    theme_bw() +
    theme(plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
          strip.text = element_markdown(size = 10, face = "bold"),
          axis.text.x = element_markdown(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
          axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
          strip.background = element_blank(),
          legend.key.height = unit(0.4, "cm"), 
          legend.key.width = unit(0.4, "cm"),
          legend.position = "right", 
          legend.title.align = 0.5,
          legend.text = element_markdown(size = 8)) +
    guides(color = guide_legend(ncol=3, override.aes = list(shape = 15, size = 3.5)),
           shape = guide_legend(ncol=1, override.aes = list(color = "black", size=2.5)),
           size = guide_legend(ncol=1))
  
  return(plot_res)
  
}




ordinate(physeq_glo_rare, method = "PCoA", distance = "bray") %>% 
  plot_pcoa(phy = physeq_glo_rare, ord=., Col="siteID", NULL) +
  geom_point(size = 2) +
  scale_color_manual(values = palette_v2) +
  theme_bw() +
  theme(plot.title = element_markdown(size = 8, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title = element_markdown(size = 7),
        axis.text.x = element_markdown(size = 7, colour = "black", hjust = 1, vjust = 1),
        axis.text.y = element_markdown(size = 7, angle = 0, hjust = 0.5),
        legend.position = "bottom") +
  labs(title="Soil AMF")


# RAXML tree -------------------------------------------------------------------
# Replace the file path with the actual location of your RAxML tree file
raxml_tree <- read.tree(file.path(data_path, "RAXML_results/RAxML_bestTree.result"))
raxml_tree <- read.tree(file.path(data_path, "RAXML_results/RAxML_bipartitions.result"))
raxml_tree

# Extract the tip labels from the tree
tip_labels <- raxml_tree$tip.label
zotu_labels <- grep("^Zotu", tip_labels, value = TRUE)
raxml_tree_filtered <- drop.tip(raxml_tree, setdiff(tip_labels, zotu_labels))

# Relabel the Zotu tips to keep only the "Zotu" and the number (remove everything after "_")
raxml_tree_filtered$tip.label <- gsub("^(Zotu[0-9]+)_.*", "\\1", raxml_tree_filtered$tip.label)
raxml_tree_filtered
class(raxml_tree_filtered)


