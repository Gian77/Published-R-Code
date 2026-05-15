#************************************************************************-------
# Manuscript Title: Site-specific factors rather than nitrogen impact arbuscular 
#                   mycorrhizal fungi diversity in bioenergy switchgrass monocultures
# Authors:          Shuang Liu, Gian Maria Niccolò Benucci, Alden Dirks, 
#                   Lukas Bell-Dereske, Sarah Evans, Gregory Bonito
# Code Developer:   Gian MN Benucci 2025
# Citation:         ...
#                   
# DOI               ...
# PMID:             ...
# **********************************************************************--------

# >>>>>>>>> PREFILTERING <<<<<<<<<< --------------------------------------------

# Load R packages. Define your package list as a character vector
pkgs <- c(
  # --- Project & Package Management ---
  "renv",            # Reproducible environments (lockfile management)
  "pak",             # Fast, modern package installation dependency engine
  "cli",             # Attractive and informative command-line interfaces  
  "styler",          # Automated R code formatting
  "janitor",         # Data cleaning
  "magrittr",        # Pipe operators (%>%)
  
  # --- Sequence Analysis & Phylogenetics ---
  "Biostrings",      # DNA/RNA/AA sequence containers
  "ape",             # Core tree handling
  "msa",             # Multiple Sequence Alignment
  "DECIPHER",        # Sequence alignment & chimera detection
  "phangorn",        # Phylogenetic analysis (MP, ML)
  "tidysq",          # Tidy processing of sequences
  "tidytree",        # Tidy manipulation of trees
  "ggtree",          # Visualization of trees
  
  # --- Microbial Ecology & Diversity ---
  "decontam",        # Contaminant OTU removal
  "phyloseq",        # Integration of OTU tables & taxonomy
  "speedyseq",       # High-performance phyloseq
  "vegan",           # Community ecology (Ordination)
  "AICcPermanova",   # Model selection for PERMANOVA
  
  # --- Data Science & Visualization ---
  "tidyverse",       # Core suite (ggplot2, dplyr, etc.)
  "ggtext",          # Markdown labels
  "ggpubr",          # Publication-ready themes
  "cowplot",         # Plot arrangement
  "gridExtra",       # Miscellaneous graphics
  "ggrepel",         # Non-overlapping labels
  "scales",          # Axis transformations
  "ggh4x",           # Extended facets
  "ggstar",          # Geometric shapes
  "ggtreeExtra",     # Annotation visualization
  "ggnewscale",      # Multiple scales
  
  # --- Statistics & Modeling ---
  "agricolae",       # Agricultural stats
  "BRCore",          # Internal core functions
  "lme4",            # Mixed-Effects Models (GLMM)
  "glmmTMB",         # Fast GLMMs
  "robustlmm",       # Outlier resistance
  "DHARMa",          # Residual diagnostics
  "parallel",        # Multi-core support
  "ggeffects",       # Marginal effects
  "sjPlot",          # Visualizing models
  "broom.mixed",     # Tidy summaries
  "merTools",        # Analyzing lme4 objects
  "multcompView",    # Significance letters
  "gllvm",           # Latent Variable Models
  "Maaslin2",        # Microbiome stats
  
  # --- Networks & Hierarchical Viz ---
  "ggdendro",        # Dendrogram data
  "igraph"            # Network visualization
)

# Load them all silently
invisible(lapply(pkgs, library, character.only = TRUE))

# Tracking package versions with renv ------------------------------------------
# renv::init()      # initializes renv in your project
# renv::restore()   # installs all packages from the lockfile
renv::snapshot()   # updates the lockfile
renv::status()
renv::install("")

# Then commit and push the updated renv.lock file!

# NOTE. Positron options to restore Rstudio projects load(".Rdata")

# R options --------------------------------------------------------------------
options(scipen = 9999, pillar.sigfig = 6, digits = 6, max.print = 10000000)
#rm(list = ls())

# Check the lib paths ----------------------------------------------------------
.libPaths()

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

as.factor(taxonomy_99$Phylum)
as.factor(taxonomy_99$Class)

# NCBI blasTAX taxonomy --------------------------------------------------------
taxonomy_99 <-
  extract_blasTAX(
    tax_path = file.path(data_path, "datasets/taxonomy_blast_99.txt"),
    namemap_path = file.path(data_path, "datasets/name_mapping_99.txt")
  ) 

head(taxonomy_99)
dim(taxonomy_99)

# Adding columns and finalize taxonomy -----------------------------------------
taxonomy_99_fix <-
  taxonomy_99 %>%
  dplyr::select(
    "Zotu", "Kingdom" ,"Phylum", "Class","Order","Family","Genus","Species") %>%
  filter(Class %in% "Glomeromycetes") %>% 
  FinalizeTaxonomy()

# These below were Unclassified 99% OTUs
taxonomy_99_fix %>% filter(Query %in% "Query_998")
taxonomy_99_fix %>% filter(Query %in% "Query_1124")

taxonomy_99_fix %>% subset(Genus %in% c("Mortierella", "Jimgerdemannia"))
taxonomy_99_fix %>% subset(Class %in% c("Mortierellomycetes", "Endogonomycetes"))

# Check taxonomy table consistency across ranks --------------------------------

CheckTaxonomyConsistency(taxonomy_99_fix, 
                         return_long=TRUE)

taxonomy_99_fix %>%
  mutate(Family = if_else(Family == "Archaeosporaceae", "Ambisporaceae", Family)) %>% 
  filter(Genus == "Ambispora")

# Fix the taxonomy table
taxonomy_99 <- 
  taxonomy_99 %>%
  mutate(Family = if_else(Family == "Archaeosporaceae", "Ambisporaceae", Family))

taxonomy_99_fix %>% subset(Class %in% c("Mortierellomycetes", "Endogonomycetes"))
taxonomy_99_fix %>%  subset(Genus %in% c("Mortierella", "Jimgerdemannia"))

head(taxonomy_99_fix)
as.factor(taxonomy_99_fix$Class)

# 99% OTUs ---------------------------------------------------------------------

otu_99 <- readDNAStringSet(
  file.path(data_path, "datasets/otus_99.fasta"),
  format = "fasta",
  seek.first.rec = TRUE,
  use.names = TRUE
)
otu_99
class(otu_99)

# WARNING. From the first round of phyloegentic tree we were able to pull out some
# chmeric non-fungal Zotus, so we are removing them before bulding the new tree.

zotus_to_remove <- 
  c("Zotu63","Zotu11850","Zotu9564","Zotu8216","Zotu4487",
    "Zotu7296","Zotu11164","Zotu8988","Zotu11484","Zotu9933",
    "Zotu9315","Zotu11229","Zotu7349","Zotu11714","Zotu10456",
    "Zotu9208","Zotu6781","Zotu8728","Zotu9297","Zotu6810",
    "Zotu6455","Zotu10620","Zotu11493","Zotu9262","Zotu6276",
    "Zotu5784","Zotu10971","Zotu6429","Zotu7095","Zotu6503",
    "Zotu10503","Zotu10263","Zotu9058","Zotu9607","Zotu621",
    "Zotu9068","Zotu7634","Zotu9266","Zotu7636","Zotu8526",
    "Zotu11324","Zotu5785","Zotu9920","Zotu9611","Zotu8450",
    "Zotu8265","Zotu945","Zotu11542","Zotu8271")

zotus_to_remove

otu_99_filtered <- otu_99[ !names(otu_99) %in% zotus_to_remove ]

# Sanity check
length(otu_99)
length(otu_99_filtered)
length(zotus_to_remove)
sum(names(otu_99) %in% zotus_to_remove)

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

# metadata  --------------------------------------------------------------------
metadata_99 <-
  read.csv(
    file = file.path(data_path, "datasets/metadata_pacbio.csv"),
    header = TRUE
  ) 

metadata_99 <- 
  metadata_99 %>%
  column_to_rownames("SampleID") %>%
  janitor::clean_names() %>% 
  mutate(
    site = site_id,
    site = factor(recode(
      site, LUX = "Lux Arbor",LC = "Lake City",HAN = "Hancock",
      RHN = "Rhinelander",ESC = "Escanaba")),
    site = fct_relevel(site, "Lux Arbor", "Lake City", "Escanaba", "Rhinelander","Hancock" ),
    fert_status = factor(recode(
      fert_status, "FERT" = "Fertilized", "UNFERT" = "Control")),
    fert_status = fct_relevel(fert_status, "Fertilized", "Control"),
    plot_rep = as.factor(plot_rep),
    ) %>%
  dplyr::select(
    site, fert_status, plot_rep, pseudo_no, collection_date)

head(metadata_99)
str(metadata_99)

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
head(physeq_AMF@tax_table)

any(physeq_AMF@tax_table %>% 
  as.data.frame() %>% 
  pull(Class) == "Mortierellomycetes")

tax_table(physeq_AMF) %>%
  as.data.frame() %>%
  pull(Class) %>%
  unique() %>%
  grep("Mortierellomycetes|Endogonomycetes", ., value = TRUE)

physeq_AMF %>% subset_taxa(Genus %in% c("Glomus"))

physeq_AMF@refseq

writeXStringSet(
  physeq_AMF@refseq,
  filepath = file.path(data_path, "datasets/otus_99_filtered.fasta"),
  format = "fasta"
)

# **********************************************************************--------
# ***** DECONTAMINATION ***** --------------------------------------------------
# decontamination from a phyloseq object ---------------------------------------

sample_data(physeq_AMF) %>% as.matrix()

sample_data(physeq_AMF)$is.neg <-
  sample_data(physeq_AMF)$site == "Control"

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

# INTERPRETATION. These looks like all real taxa to me. I am not going to 
# remove any taxa! I will just remove the control sample(s).

# Filtering the phyloseq object
physeq_AMF_clean <-
physeq_AMF %>%
  subset_samples(!is.neg %in% TRUE) %>% # remove the control samples
  prune_taxa(taxa_sums(x = .) > 0, x = .) # make sure there aren't OTUs that are 0

physeq_AMF_clean
physeq_AMF_clean@sam_data

# ***** FINAL PHYLOSEQ OBJECTS ***** -------------------------------------------
# Final phyloseq object for downstream analysis ---------------------
physeq_AMF_clean
physeq_AMF_clean@sam_data

# save the metadata
write.csv(
  x = as.data.frame(as.matrix(physeq_AMF_clean@sam_data)),
  file = file.path(data_path, "datasets/medatata_amf.csv")
)

# **********************************************************************--------