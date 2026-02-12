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
# load(".Rdata")

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
  tidytree,
  vegan,
  AICcPermanova,
  tidyverse,
  ggtext,
  ggpubr,
  ggtree,
  cowplot,
  gridExtra,
  ggrepel,
  scales,
  agricolae,
  BRCore,
  install=FALSE
)



# Session Info -----------------------------------------------------------------
sessionInfo()

# **********************************************************************--------
# ***** PATHS ***** ------------------------------------------------------------

data_path <-
  ("/home/gian/Dropbox/6_PROJETCS/Published-R-Code/Liu_etal_2026_AMF_FertUnfertSwitchgrass")

data_path

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

# I run this analysis on the HPCC. Then I import the phylogenetic tree in here.

# NOTE. I am importing 2 trees. One genrated with RAXML, the other generated with 
# iqtree2 that has bayesian posterior probabilities.
tree_raxml <- 
  read.tree("phylogeny/otus99_mafft_trim_spprt.raxml.support") 

str(tree_raxml)
ggtree::ggtree(tree_raxml)

tree_iqtree2 <- read.tree("phylogeny/otus99_mafft_trim_iq2.treefile")

str(tree_iqtree2)
ggtree::ggtree(tree_iqtree2)

# **********************************************************************--------
# ***** MAKE PHYLOSEQ OBJECT ***** ---------------------------------------------
dim(otutable_99)
dim(metadata_99)
dim(taxonomy_99_fix)
str(otu_99)
str(tree_raxml)

# Final phyloseq object --------------------------------------------------------
physeq_AMF <-
  phyloseq(
    otu_table(otutable_99, taxa_are_rows = TRUE),
    sample_data(metadata_99),
    tax_table(as.matrix(taxonomy_99_fix)),
    DNAStringSet(otu_99), 
    phy_tree(tree_raxml)) %>% 
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.)

physeq_AMF

physeq_AMF@phy_tree
head(physeq_AMF@otu_table)
head(physeq_AMF@sam_data)

# **********************************************************************--------
# ***** DECONTAMINATION ***** --------------------------------------------------
# decontamination from a phyloseq object ---------------------------------------
sample_data(physeq_AMF) %>% as.matrix()

sample_data(physeq_AMF)$is.neg <-
  sample_data(physeq_AMF)$description == "Control"

contam_AMF <-
  decontam::isContaminant(physeq_AMF,
                          method = "prevalence",
                          neg = "is.neg",
                          threshold = 0.1
  )

contam_AMF
table(contam_AMF$contaminant)
contam_AMF %>% filter(contaminant == TRUE)

# Check contaminants taxonomy
left_join(
  contam_AMF %>% 
    filter(contaminant == TRUE) %>% 
    rownames_to_column("OTU_ID"), 
  taxonomy_99_fix%>% 
    rownames_to_column("OTU_ID"),
  by = "OTU_ID") %>% 
  dplyr::select( OTU_ID, freq, BestMatch) %>% 
  left_join(
    otutable_99 %>% 
      rownames_to_column("OTU_ID") %>% 
      filter(OTU_ID %in% c(contam_AMF %>% 
                             filter(contaminant == TRUE) %>% 
                             rownames_to_column("OTU_ID") %>% 
                             pull(OTU_ID))) %>% 
      mutate(Abund = rowSums(across(where(is.numeric)))) %>% 
      dplyr::select(OTU_ID, Abund),
    by = "OTU_ID")

# NOTE! These looks like all real taxa to me. I am not going to remove any!

# function to remove taxa by OTU name
remove_taxa <- function(physeq, badTaxa) {
  allTaxa <- taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

subset(contam_AMF, contaminant %in% c("TRUE"))

# Filtering the phyloseq object
physeq_AMF_clean <-
  remove_taxa(
    physeq_AMF,
    rownames(subset(contam_AMF, contaminant %in% c("TRUE")))
  ) %>%
  subset_samples(!is.neg %in% TRUE) %>% # remove the control samples
  prune_taxa(taxa_sums(x = .) > 0, x = .) # make sure there aren't OTUs that are 0

# ***** FINAL PHYLOSEQ OBJECTS ***** -------------------------------------------
physeq_AMF_clean
physeq_AMF_clean@sam_data

# ***** RAREFACTION ***** ------------------------------------------------------

# Add rarefaction metrics
physeq_AMF_clean <- add_rarefaction_metrics(data = physeq_AMF_clean)
physeq_AMF_clean@sam_data %>% as.matrix()

rarefaction_plot <- plot_rarefaction_metrics(physeq_AMF_clean)
print(rarefaction_plot)

# Identify best rarefaction depth cutoff 
physeq_AMF_clean@sam_data %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  arrange(read_num)
  
rare_depth_cutoff = 4568

# Perform mutiple rarefaction
AMF_otutable_rarefied <-
  multi_rarefy(
    physeq = physeq_AMF_clean,
    depth_level = rare_depth_cutoff,
    num_iter = 100,
    threads = 8,
    set_seed = 1026
  )

rowSums(AMF_otutable_rarefied)
dim(AMF_otutable_rarefied)

# Update otu_table
physeq_AMF_rare <-
  do_phyloseq(physeq = physeq_AMF_clean, otu_rare = AMF_otutable_rarefied)

physeq_AMF_rare
sample_sums(physeq_AMF_rare)

# Update metadata 
sample_data(physeq_AMF_rare)$site_id

sample_data(physeq_AMF_rare) <- sample_data(
  as.data.frame(as.matrix(physeq_AMF_rare@sam_data)) %>%
    mutate(
      site = site_id,
      site = recode(
        site,
        LUX = "Lux Harbor",
        LC = "Lake City",
        HAN = "Hancock",
        RHN = "Rhinelander",
        ESC = "Escanaba"
      )
    ) %>%
    select(
      collection_date,
      site_id,
      site,
      fert_status,
      plot_rep,
      pseudo_no,
      x_cord,
      y_cord
    )
)

sample_data(physeq_AMF_rare) 

# save the metadata
write.csv(
  x = as.data.frame(as.matrix(physeq_AMF_rare@sam_data)) %>%
    mutate(
      site = site_id,
      site = recode(
        site,
        LUX = "Lux Harbor",
        LC = "Lake City",
        HAN = "Hancock",
        RHN = "Rhinelander",
        ESC = "Escanaba"
      )
    ) %>%
    select(
      collection_date,
      site_id,
      site,
      fert_status,
      plot_rep,
      pseudo_no,
      x_cord,
      y_cord
    ),
  file = file.path(data_path, "datasets/medatata_filtered.csv")
)






# **** ALPHA DIVERSITY **** ----------------------------------------------------
# Adding alpha metrics ---------------------------------------------------------
AlphaMetrics <- function(physeq) {
  sample_data(physeq)$ReadNo <- sample_sums(physeq)
  sample_data(physeq)$hill_0 <- as.data.frame(as.matrix(t(physeq@otu_table))) %>% 
    renyi(scales = c(0), hill = TRUE)
  sample_data(physeq)$hill_1 <- as.data.frame(as.matrix(t(physeq@otu_table))) %>% 
    renyi(scales = c(1), hill = TRUE)
  sample_data(physeq)$hill_2 <- as.data.frame(as.matrix(t(physeq@otu_table))) %>% 
    renyi(scales = c(2), hill = TRUE)
  sample_data(physeq)$pielou_j <- log(sample_data(physeq)$hill_1) / log(sample_data(physeq)$hill_0)
  return(physeq)
}

# NOTE. Pileau J is the classic 0–1 evenness measure based on Shannon.

physeq_AMF_rare <- 
  AlphaMetrics(physeq_AMF_rare)
physeq_AMF_rare

head(physeq_AMF_rare@sam_data)
AlphaMetrics(physeq_AMF_rare) %>% sam_data() %>% pull(pielou_j) %>% range()










