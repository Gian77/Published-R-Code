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
physeq_AMF_clean
physeq_AMF_clean@sam_data

# save the metadata
write.csv(
  x = as.data.frame(as.matrix(physeq_AMF_clean@sam_data)),
  file = file.path(data_path, "datasets/medatata_amf.csv")
)

# **********************************************************************--------
# **** 1. ALPHA DIVERSITY **** -------------------------------------------------

# NOTE. I am not attaching the metrics as these metrics are not the final one
# calculated over iterations.

add_rarefaction_metrics(physeq_AMF_clean) %>% 
  sample_data() %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  arrange(read_num)

# ***** FIGURE S1 - identify depth ***** ---------------------------------------
rarefaction_plot <- 
  add_rarefaction_metrics(physeq_AMF_clean) %>% 
  plot_rarefaction_metrics()

rarefaction_plot 

ggsave(plot = rarefaction_plot, 
       path = data_path,
       filename = "results/Fig_S1_identify_raredepth_plots.pdf")


# **********************************************************************--------
# ***** RAREFACTION CURVES ***** -----------------------------------------------

rare_depth_cutoff = 6512

rarecurve_AMF <-
  physeq_AMF_clean@otu_table %>%
  t() %>%
  as.matrix() %>%
  as.data.frame() %>%
  rarecurve(x = ., step = 20, sample = rare_depth_cutoff, tidy = TRUE)

head(rarecurve_AMF)
dim(rarecurve_AMF)

Fig_S2_rarecurve <-
rarecurve_AMF %>%
  dplyr::rename(sample_id = Site) %>%
  left_join(.,
            physeq_AMF_clean@sam_data %>%
              as.matrix() %>%
              as.data.frame() %>%
              rownames_to_column("sample_id"),
            by = "sample_id") %>%
  ggplot(aes(x = Sample, y = Species, group = sample_id, color = fert_status)) +
  geom_line() +
  geom_vline(xintercept = rare_depth_cutoff, color = "black", linetype = "dashed") +
  theme_bw() +
  theme(
    plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5, vjust = 0.5),
    plot.subtitle = element_markdown(size = 10, face = "bold", hjust = 0.5, vjust = 0.5),
    axis.text.x = element_markdown(angle = 0, size = 10, hjust = 0.5, vjust = 0.5),
    axis.text.y = element_markdown(size = 10),
    axis.title = element_markdown(size = 10, face = "bold", hjust = 0.5, vjust = 0.5),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = c(0.98, 0.02),
    legend.justification = c("right", "bottom")) +
  labs(title = "Rarefaction curves Arbuscular Mycorrhizal Fungi (AMF)", 
       x = "Number of DNA reads", 
       y = "Number of 99% OTUs") 


Fig_S2_rarecurve

ggsave(plot = Fig_S2_rarecurve, 
       path = data_path,
       filename = "results/Fig_S2_rarefaction_curves.pdf")


# Extarct data.frame from phyloseq object
meta_AMF <- 
  as.data.frame(as.matrix(physeq_AMF_clean@sam_data)) %>% 
  mutate(
    site = as.factor(site),
    site = fct_relevel(site, "Lux Arbor", "Lake City", "Escanaba", "Rhinelander","Hancock" ),
    site_plot = factor(paste(site, plot_rep, sep=":")),
    fert_status = as.factor(fert_status),
    fert_status = fct_relevel(fert_status, "Fertilized", "Control"),
  )

str(meta_AMF)

otutable_AMF <- as.data.frame(t(as.matrix(physeq_AMF_clean@otu_table)))

str(otutable_AMF)

# NOTE. E[D(X)] instead of  D(E[X]) is the right way to calculate alpha metrics. 
# First you calculate the diversity metric for each rarefied iteration, then you
# average those diversity metrics across all iterations to get the expected value
# of the diversity metric after rarefaction.

# Calculate hill numbers
hill_all <- map_dfr(c(0, 1, 2), function(q) {
  calc_hill_rarefied(
    otu_mat = otutable_AMF,
    depth   = rare_depth_cutoff,
    n_iter  = 100,
    q       = q
  ) %>%
    dplyr::select(sample_id, hill_mean) %>%
    mutate(q = paste0("q", q))
})  %>%
  pivot_wider(
    names_from  = q,
    values_from = hill_mean,
    names_prefix = "hill_"
    )

hill_all

alpha_df <-
  meta_AMF %>%
  rownames_to_column("sample_id") %>% 
  inner_join(hill_all, by = "sample_id") %>% 
  mutate(across(.cols = 1:4, .fns = as.factor))

str(alpha_df)
head(alpha_df)


# hill_q0 -----------------------------------------------------------------------

alpha_df %>% 
  ggplot2::ggplot(aes(x = hill_q0)) +
  geom_histogram(binwidth = 20) +
  labs(title = "hill_q0") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(hill_q0)), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)

alpha_df %>% 
  ggplot(aes(x = log(hill_q0))) +
  geom_histogram(binwidth = 0.1) +
  labs(title = "hill_q0") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(log(hill_q0))), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)

# Model Comparison -------------------------------------------------------------
# Does fertilizer affect richness? H0: mean richness (FERT) = 0
summary(run_lmem(alpha_df, hill_q0, "baseline"))
summary(run_lmem(alpha_df, hill_q0, "fixslope"))
summary(run_lmem(alpha_df, hill_q0, "randomslope"))

# Does fertilizer improve fit?
anova( run_lmem(alpha_df, hill_q0, "baseline"), 
       run_lmem(alpha_df, hill_q0, "fixslope") )

# Does slope vary among sites?
anova( run_lmem(alpha_df, hill_q0, "fixslope"), 
       run_lmem(alpha_df, hill_q0, "randomslope") )

# INTERPRETATION. So the statistically supported model is the baseline model with 
# no fertilizer effect and no slope variation among sites. Richness varies strongly
# among plots and moderately among sites, but not due to fertilizer.

# hill_q1 -----------------------------------------------------------------------

alpha_df %>% 
  ggplot2::ggplot(aes(x = hill_q1)) +
  geom_histogram(binwidth = 15) +
  labs(title = "hill_q1") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(hill_q1)), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)

alpha_df %>% 
  ggplot(aes(x = log(hill_q1))) +
  geom_histogram(binwidth = 0.1) +
  labs(title = "hill_q1") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(log(hill_q1))), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)

# Model Comparison -------------------------------------------------------------
# Does fertilizer affect Shannon? H0: mean Shannon (FERT) = 0
summary(run_lmem(alpha_df, hill_q1, "baseline"))
summary(run_lmem(alpha_df, hill_q1, "fixslope"))
summary(run_lmem(alpha_df, hill_q1, "randomslope"))

# Does fertilizer improve fit?
anova( run_lmem(alpha_df, hill_q1, "baseline"), 
       run_lmem(alpha_df, hill_q1, "fixslope") )

# Does slope vary among sites?
anova( run_lmem(alpha_df, hill_q1, "fixslope"), 
       run_lmem(alpha_df, hill_q1, "randomslope") )

# The "Zero Variance" Culprit
# Look closely at the Random Effects section for all three models. You will see 
# this line: site (Intercept)  0.000  0.000
# A "singular fit" occurs when one of your random effects has a variance of exactly
# (or very near) zero.

# hill_q2 -----------------------------------------------------------------------

alpha_df %>% 
  ggplot2::ggplot(aes(x = hill_q2)) +
  geom_histogram(binwidth = 10) +
  labs(title = "hill_q2") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(hill_q2)), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)

alpha_df %>% 
  ggplot(aes(x = log(hill_q2))) +
  geom_histogram(binwidth = 0.1) +
  labs(title = "hill_q2") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(log(hill_q2))), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)

# Model Comparison -------------------------------------------------------------
# Does fertilizer affect Shannon? H0: mean Shannon (FERT) = 0
summary(run_lmem(alpha_df, hill_q2, "baseline"))
summary(run_lmem(alpha_df, hill_q2, "fixslope"))
summary(run_lmem(alpha_df, hill_q2, "randomslope"))

# Does fertilizer improve fit?
anova( run_lmem(alpha_df, hill_q2, "baseline"), 
       run_lmem(alpha_df, hill_q2, "fixslope") )

# Does slope vary among sites?
anova( run_lmem(alpha_df, hill_q2, "fixslope"), 
       run_lmem(alpha_df, hill_q2, "randomslope") )

# ***** MODEL DIAGNOSTICS ****** -----------------------------------------------

# DHARMA diagnostics -----------------------------------------------------------
# Generate simulated residuals (n=1000 is standard) and plot.

diagnostics_dharma(
  model = run_lmem(alpha_df, hill_q0, "baseline"),
  group_var1 = alpha_df$site,
  group_var2 = alpha_df$plot_rep
)

diagnostics_dharma(
  model = run_lmem(alpha_df, hill_q0, "fixslope"),
  group_var1 = alpha_df$site,
  group_var2 = alpha_df$plot_rep
)

diagnostics_dharma(
  model = run_lmem(alpha_df, hill_q1, "baseline"),
  group_var1 = alpha_df$site,
  group_var2 = alpha_df$plot_rep
)

diagnostics_dharma(
  model = run_lmem(alpha_df, hill_q2, "fixslope"),
  group_var1 = alpha_df$site,
  group_var2 = alpha_df$plot_rep
)


# ***** VISUALIZING MIX MODELS ***** ------------------------------------------

# 1. Basic Boxplot by Fertilizer Status ----------------------------------------
ggplot(alpha_df, 
       aes(x = fert_status, y = hill_q0, fill = fert_status)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(title = "Hill Number by Fertilizer Status",
       x = "Fertilizer Status",
       y = "Hill Number (Diversity)") +
  theme_minimal() +
  theme(legend.position = "none")

# 2. Show Random Effects (Sites) -----------------------------------------------
# Color by site to show between-site variability
ggplot(alpha_df, 
       aes(x = fert_status, y = hill_q0, color = site)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5),
             alpha = 0.6, size = 2) +
  stat_summary(aes(group = site), fun = mean, geom = "line", 
               position = position_dodge(0.5)) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", 
               color = "black", size = 4, shape = 18) +
  labs(title = "Hill Number by Fertilizer Status",
       subtitle = "Colors = different sites, black diamonds = fixed effect means",
       x = "Fertilizer Status",
       y = "Hill Number",
       color = "Site") +
  theme_minimal()

# 3. Facet by Site (Best for Seeing Nested Structure) --------------------------

# Show each site separately
ggplot(alpha_df, 
       aes(x = fert_status, y = hill_q0, fill = fert_status)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.5, aes(color = plot_rep)) +
  stat_summary(aes(group = site), fun = mean, geom = "line", 
               position = position_dodge(0.5)) +
  facet_wrap(~ site, ncol = 5) +
  labs(title = "Hill Number by Fertilizer Status Across Sites",
       subtitle = "Each panel = one site, shows within-site variation",
       x = "Fertilizer Status",
       y = "Hill Number") +
  theme_minimal() +
  theme(legend.position = "bottom")

# 4. Model Predictions with Random Effects (ggeffects) -------------------------

# Site-level predicted means
ggeffects::ggpredict(run_lmem(alpha_df, hill_q0, "baseline"), 
                     terms = "site", type = "random") %>% plot()

ggeffects::ggpredict(run_lmem(alpha_df, hill_q0, "fixslope"), 
            terms = c("fert_status", "site")) %>% plot() 


fit_hill_q0_base <- lmer(
  hill_q0 ~ 1 + (1 | site) + (1 | site_plot),
  data = alpha_df,
  REML = FALSE
)

pred_hill_q0_base <-
  ggpredict(fit_hill_q0_base, terms = "site_plot", type = "random")

plot(pred_hill_q0_base) +
  labs(title = "Model Predictions: Hill Number by Fertilizer Status",
       subtitle = "Lines = site-specific predictions (random effects)",
       x = "Fertilizer Status",
       y = "Predicted Hill Number") 

# 5. Show Fixed + Random Effects (merTools) ------------------------------------

# Extract fixed and random effects
cbind(
  alpha_df,
  run_lmem(alpha_df, hill_q0, "baseline") %>% 
  merTools::predictInterval(newdata = alpha_df, 
                            which = "full", 
                            level = 0.95)) %>% 
  rename(fitted = fit) %>%
  ggplot(aes(x = fert_status, y = hill_q0)) +
  geom_jitter(aes(color = site), width = 0.2, alpha = 0.5) +
  geom_point(aes(y = fitted, group = site, color = site), 
             size = 3, shape = 17) +
  geom_errorbar(aes(ymin = lwr, ymax = upr, group = site, color = site),
                width = 0.2, alpha = 0.6) +
  labs(title = "Observed vs Predicted Hill Numbers",
       subtitle = "Triangles = model predictions per site",
       x = "Fertilizer Status",
       y = "Hill Number",
       color = "Site") +
  theme_minimal()

# 6. Plot Random Effects Deviations (broom.mixed) ------------------------------

broom.mixed::tidy(
  run_lmem(alpha_df, hill_q0, "baseline"), effects = "ran_vals") %>% 
  ggplot(aes(x = estimate, y = level, color = group)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = estimate - 1.96*std.error, 
                    xmax = estimate + 1.96*std.error),
                 height = 0.2, 
                orientation = "y") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Random Effect Deviations by Site and Plot",
       x = "Deviation from Fixed Effect",
       y = "Group Level",
       color = "Random Effect") +
  theme_minimal()

broom.mixed::tidy(
  run_lmem(alpha_df, hill_q0, "fixslope"), effects = "ran_vals") %>% 
  ggplot(aes(x = estimate, y = level, color = group)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = estimate - 1.96*std.error, 
                    xmax = estimate + 1.96*std.error),
                height = 0.2, 
                orientation = "y") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Random Effect Deviations by Site and Plot",
       x = "Deviation from Fixed Effect",
       y = "Group Level",
       color = "Random Effect") +
  theme_minimal()

# INTERPRETATION. These are Random effects plots so they do not change between
# baseline and fixslope models.

# 9. Full Diagnostic Plot Set (sjPlot) -----------------------------------------

sjPlot::plot_model(run_lmem(alpha_df, hill_q0, "baseline"), type = "re")
sjPlot::plot_model(run_lmem(alpha_df, hill_q0, "fixslope"), type = "re")

# Plot fixed effects with confidence intervals (can't do on the baseline as 
# there is no fixed effects)
plot_model(run_lmem(alpha_df, hill_q0, "fixslope"), type = "est", 
           show.values = TRUE, 
           value.offset = 0.3) +
  labs(title = "Fixed Effect Estimates",
       subtitle = "Effect of fertilizer status on Hill number") +
  theme_minimal()

# Comprehensive model diagnostics and predictions
tab_model(run_lmem(alpha_df, hill_q0, "baseline"), show.stat = TRUE)
tab_model(run_lmem(alpha_df, hill_q0, "fixslope"), show.stat = TRUE) 

tab_model(run_lmem(alpha_df, hill_q1, "baseline"), show.stat = TRUE)
tab_model(run_lmem(alpha_df, hill_q1, "fixslope"), show.stat = TRUE) 


# 10. Raw Data + Model Predictions (emmeand) -----------------------------------

# Panel A: Raw data
ggplot(alpha_df, aes(x = fert_status, y = hill_q0, fill = fert_status)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(title = "A. Raw Data",
       x = "Fertilizer Status", y = "Hill Number") +
  theme_minimal() +
  theme(legend.position = "none")

# Panel B: Model predictions
emmeans::emmeans(run_lmem(alpha_df, hill_q0, "fixslope"), ~ fert_status) %>% 
  as.data.frame() %>% 
  ggplot(aes(x = fert_status, y = emmean, fill = fert_status)) +
  geom_col(alpha = 0.7) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  labs(title = "B. Model-Predicted Means",
       subtitle = "Error bars = 95% CI",
       x = "Fertilizer Status", y = "Predicted Hill Number") +
  theme_minimal() +
  theme(legend.position = "none")

# Selecting BEST VISUALIZATION APPROACH ----------------------------------------

palette_fert <- c("#CC2D35", "#2D3142")

palette_site <- c("#009E73", "#0072B2", "#825121", "#E69F00", "#CC79A7")

palette_taxa <-c("#D21E2C","#058ED9","#ae09ea","#dba4a4","#117744",
                 "#F7F7C5","#283dff","#521899","#82807f","#014443",
                 "#FDDB8E","#bfc5ff","#111b77","#d8d6d4","#b7ffdb",
                 "#fcb067","#ffb7ef","#560d0d","#60ffaf", "#A5A518", "#000000") 

palette_taxa <-c("#D21E2C","#058ED9","#ae09ea","#dba4a4","#117744",
                 "#F7F7C5","#283dff","#521899","#82807f","#014443",
                 "#FDDB8E","#bfc5ff","#111b77","#d8d6d4","#b7ffdb",
                 "#fcb067","#ffb7ef","#560d0d","#60ffaf", "#A5A518", "white") 
                 
c("#ea7f17","#fa7efc","#a35151","#825121",)

# Hill 0 and Hill 2 ------------------------------------------------------------

# INTERPRETATION. Hill numbers are not independent indices. They are the same diversity 
# continuum at different q values. So including all three is often redundant unless 
# you want to demonstrate consistency across evenness weighting.
# hill_q0 (q = 0) Pure richness, sensitive to rare taxa, influenced by sampling noise
# hill_q2 (q = 2) Strongly dominated by abundant taxa, less sensitive to rare taxa, 
# more stable across samples.
# hill_q0 + hill_q2 They bracket the diversity spectrum "Fertilization did not affect 
# richness (q=0) nor diversity weighted toward dominant taxa (q=2)."

# If hill_q0 = 150 → there are 150 observed taxa
# If hill_q2 = 24, it means: The community has the same dominance structure as a 
# community with 24 equally abundant species.


# calculate global limits across all strains
x_limits_hill_q0 <- 
  alpha_df %>%
  group_by(site) %>% 
  group_map(~broom.mixed::tidy(run_lmem(alpha_df, hill_q0, "baseline"), 
                               effects = "ran_vals")) %>%
  bind_rows() %>%
  summarise(
    min = min(estimate - 1.96*std.error),
    max = max(estimate + 1.96*std.error)
  )

x_limits_hill_q0

intercept_hill_q0 <- fixef(run_lmem(alpha_df, hill_q0, "baseline"))["(Intercept)"]

x_limits_hill_q2 <- 
  alpha_df %>%
  group_by(site) %>% 
  group_map(~broom.mixed::tidy(run_lmem(alpha_df, hill_q2, "baseline"), 
                               effects = "ran_vals")) %>%
  bind_rows() %>%
  summarise(
    min = min(estimate - 1.96*std.error),
    max = max(estimate + 1.96*std.error)
  )

x_limits_hill_q2

intercept_hill_q2 <- fixef(run_lmem(alpha_df, hill_q2, "baseline"))["(Intercept)"]

# plotting
Figure_1_alpha <-
  ggarrange(
    ggarrange(
# hill_q0: Raw data with fixed effect means  
ggplot(alpha_df, 
       aes(x = fert_status, y = hill_q0, color = site)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5),
             alpha = 0.8, size = 2) +
  stat_summary(aes(group = site), fun = mean, geom = "line", 
               position = position_dodge(0.5)) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", 
               color = "black", size = 4, shape = 18) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", 
               color = "black", size = 1) +
  labs(title = "Hill 0 (richness)",
       subtitle = "Black line connects fixed effect means",
       x = "Treatment",
       y = "OTUs",
       color = NULL) +
  scale_color_manual(values = palette_site) +
  theme_classic() +
  theme(
    plot.title = element_markdown(size =12, face = "bold",hjust = 0.5, vjust = 0.5),
    plot.subtitle = element_markdown(size = 10,hjust = 0.5, vjust = 0.5),
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 8),
    legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"),
    legend.title = element_blank(), legend.text = element_text(size = 8)
    #legend.margin=ggplot2::margin(0,5,0,0),
    #legend.box.margin=ggplot2::margin(0,5,0,0)
    )+
  guides(color = guide_legend(ncol=1)),

# hill_q2: Raw data with fixed effect means  
ggplot(alpha_df, 
       aes(x = fert_status, y = hill_q2, color = site)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5),
             alpha = 0.8, size = 2) +
  stat_summary(aes(group = site), fun = mean, geom = "line", 
               position = position_dodge(0.5)) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", 
               color = "black", size = 4, shape = 18) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", 
               color = "black", size = 1) +
  labs(title = "Hill 2 (inverse Simpson) ",
       subtitle = "Black line connects fixed effect means",
       x = "Treatment",
       y = "Effective OTU number",
       color = NULL) +
  scale_color_manual(values = palette_site) +
  theme_classic() +
  theme(
    plot.title = element_markdown(size =12, face = "bold",hjust = 0.5, vjust = 0.5),
    plot.subtitle = element_markdown(size = 10, hjust = 0.5, vjust = 0.5),
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 8),
    legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"),
    legend.title = element_blank(), legend.text = element_text(size = 8)
    #legend.margin=ggplot2::margin(0,5,0,0),
    #legend.box.margin=ggplot2::margin(0,5,0,0)
  )+
  guides(color = guide_legend(ncol=1)),

ncol = 2, 
nrow = 1,
labels = c("A", "C"),
legend = "right",
common.legend = TRUE),

ggarrange(
# hill_q0: Random effect deviations
broom.mixed::tidy(
  run_lmem(alpha_df, hill_q0, "baseline"), effects = "ran_vals") %>% 
  as.data.frame() %>%
  separate(level, into = c("plot_rep", "site"), sep = ":", 
           extra = "merge", fill = "left") %>% 
  mutate(level = ifelse(is.na(plot_rep), site, paste(site, plot_rep, sep = ":"))) %>% 
  mutate(group = recode(group, "site" = "Site", 
                        "plot_rep:site" = "Site:Plot")) %>%
  arrange(level) %>% 
  ggplot(aes(x = estimate, y = level, color = group)) +
  geom_point(size = 4, shape = 18) +
  geom_errorbar(aes(xmin = estimate - 1.96*std.error, 
                    xmax = estimate + 1.96*std.error),
                height = 0, 
                orientation = "y") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Hill 0 Random Effects",
       subtitle = NULL,
       x = "Deviations from fixed-effect means",
       y = "Plot:site",
       color = "Plot:site") +
  scale_color_manual(values = palette_fert)+
  theme_classic()+ 
  theme(
    plot.title = element_markdown(size =12, face = "bold",hjust = 0.5, vjust = 0.5),
    plot.subtitle = element_markdown(size = 10, hjust = 0.5, vjust = 0.5),
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 8),
    legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"),
    legend.title = element_blank(), legend.text = element_text(size = 8)
    #legend.margin=ggplot2::margin(0,5,0,0),
    #legend.box.margin=ggplot2::margin(0,5,0,0)
    ) +
  guides(color = guide_legend(ncol = 1)) +
  scale_x_continuous(
    limits = c(x_limits_hill_q0$min, x_limits_hill_q0$max),
    sec.axis = sec_axis(
      transform = ~ . + intercept_hill_q0,
      name ="OTUs")
    ),

# hill_q2: Random effect deviations
broom.mixed::tidy(
  run_lmem(alpha_df, hill_q2, "baseline"), effects = "ran_vals") %>% 
  as.data.frame() %>%
  separate(level, into = c("plot_rep", "site"), sep = ":", 
           extra = "merge", fill = "left") %>% 
  mutate(level = ifelse(is.na(plot_rep), site, paste(site, plot_rep, sep = ":"))) %>% 
  mutate(group = recode(group, "site" = "Site", 
                        "plot_rep:site" = "Site:Plot")) %>%
  arrange(level) %>% 
  ggplot(aes(x = estimate, y = level, color = group)) +
  geom_point(size = 4, shape = 18) +
  geom_errorbar(aes(xmin = estimate - 1.96*std.error, 
                    xmax = estimate + 1.96*std.error),
                height = 0, 
                orientation = "y") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Hill 2 Random Effects",
       subtitle = NULL,
       x = "Deviations from fixed-effect means",
       y = "Plot:site",
       color = "Plot:site") +
  scale_color_manual(values = palette_fert)+
  theme_classic()+ 
  theme(
    plot.title = element_markdown(size =12, face = "bold",hjust = 0.5, vjust = 0.5),
    plot.subtitle = element_markdown(size = 10, hjust = 0.5, vjust = 0.5),
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 8),
    legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"),
    legend.title = element_blank(), legend.text = element_text(size = 8)
    #legend.margin=ggplot2::margin(0,5,0,0),
    #legend.box.margin=ggplot2::margin(0,5,0,0)
    ) +
  guides(color = guide_legend(ncol = 1)) +
  scale_x_continuous(
    limits = c(x_limits_hill_q2$min, x_limits_hill_q2$max),
    sec.axis = sec_axis(
      transform = ~ . + intercept_hill_q2,
      name = "inverse Simpson")
    ),

ncol = 2, 
nrow = 1,
labels = c("B","D"),
legend = "right",
common.legend = TRUE),
ncol = 1, 
nrow = 2, heights=c(0.6, 0.8)
)

Figure_1_alpha

# ***** FIGURE 1 - Alpha-diversity ***** ---------------------------------------
ggsave(
  file.path(data_path, "results/Fig_1_alphadiv.pdf"),
  plot = ggpubr::annotate_figure(
    Figure_1_alpha,
    top = text_grob("ALPHA DIVERSITY", size = 12, face = "bold")
  ),
  device = "pdf"
)


# ***** FUTHER COMPARISONS ***** -----------------------------------------------

# Additional question we can have: Which site is richer in term of diversity
# metrics ?

fit_hill0_site <- lmer(hill_q0 ~ site + (1 | site:plot_rep),
               data = alpha_df, REML = FALSE)

anova(fit_hill0_site)   # overall test: does richness differ among sites?
summary(fit_hill0_site)

# Pairwise comparisons: “which sites differ?”
fit_hill0_site_emm <- emmeans::emmeans(fit_hill0_site, ~ site)

fit_hill0_site_emm                       # site means (model-adjusted)
pairs(fit_hill0_site_emm, adjust = "tukey")  # Tukey-corrected pairwise tests
confint(fit_hill0_site_emm)              # 95% CIs for each site mean

as.data.frame(fit_hill0_site_emm) |>
  dplyr::arrange(dplyr::desc(emmean))

# if you also want to control for fertilizer while comparing sites

m_site_adj <- lmer(hill_q0 ~ site + fert_status + (1 | site:plot_rep),
                   data = alpha_df, REML = FALSE)

emm_site_adj <- emmeans::emmeans(m_site_adj, ~ site)
pairs(emm_site_adj, adjust = "tukey")

# This answers: “Which sites are richer after adjusting for fertilizer?”

# NOTE. The (1 | site:plot_rep) still belong to the model even when testing for
# differences between sites, because plots are still your experimental units, 
# and you have 3 pseudoreplicates per plot. Even when site is fixed, the 
# hierarchical structure does not disappear. 


# **********************************************************************--------
# **** 2. BETA DIVERSITY **** --------------------------------------------------

set.seed(26031)

# Generating rarefied distance matrix ------------------------------------------
rare_depth_cutoff = 6512

dist_avg_AMF <- avgdist(
  x = otutable_AMF,
  sample = rare_depth_cutoff,
  iterations = 100,
  dmethod = "bray"
)

dist_avg_AMF
labels(dist_avg_AMF)

# Filter out the samples that were dropped sure to rarefaction from the metadata file
meta_AMF_rare <- meta_AMF[ labels(dist_avg_AMF), ]
dim(meta_AMF_rare)

# Verify alignment (should be TRUE)
identical(labels(dist_avg_AMF), rownames(meta_AMF[ labels(dist_avg_AMF), ]))

# PERMANOVA ------------------------------------------------------------------- 

# INTERPRETATION. The adonis2 function in vegan can be run using:

# >>>> 1) by = "terms" (sequential/Type I) -------------------------------------
# Tests each term after only the terms that come before it in the formula. This 
# means the results are order-dependent — putting site before fert_status gives 
# site "first dibs" on explaining variance, and fert_status is tested on whatever
# is left over. This is appropriate when you have a strong a priori reason to 
# believe terms should be entered in a specific hierarchical order.

# >>>> 2) by = "margin" (marginal/Type III) ------------------------------------
# Tests each term after all other terms are already in the model, regardless of 
# order. This is what you've been using, and for your question — "does fertilization
# matter after accounting for everything else?" — it's the right choice.

# Since your primary question is whether fertilization has an effect independent 
# of site and plot structure, by = "margin" is the correct and more conservative 
# test. Using by = "terms" would actually be more likely to find a spurious 
# fertilization effect if you put fert_status early in the formula before site
# absorbs its share of variance.

# 1) Does fertilization affect community composition? --------------------------
# treatment is applied at the plot level, the appropriate unit of permutation 
# within site and among plots, not withing plots. This is a correct test for 
# fertilizer while accounting for site structure.

adonis2(
  dist_avg_AMF ~ fert_status,
  data = meta_AMF_rare,
  method = "bray",
  permutations = how(blocks = meta_AMF_rare$site, nperm = 999),
  by = "margin"
)


# NOTE that using permutations = how(blocks = meta_rare$site) is the same as 
# adding the strata = site was used to restrict permutations within groups i.e.
# treatment labels are shuffled only within each site.

# 2) Does site impact community composition? -----------------------------------

# This model tests site effect, fertilization effect within site, and keeps 
# permutations restricted within site. Gives marginal tests (Type III-like logic)

adonis2(
  dist_avg_AMF ~ site,
  data = meta_AMF_rare,
  method = "bray",
  permutations = 999,
  by = "term"
)

adonis2(
  dist_avg_AMF ~ fert_status + site,
  data = meta_AMF_rare,
  method = "bray",
  permutations = 999,
  by = "term"
)

adonis2(
  dist_avg_AMF ~ fert_status + site,
  data = meta_AMF_rare,
  method = "bray",
  permutations = how(blocks = meta_AMF_rare$site, nperm = 999),
  by = "margin"
)

adonis2(
  dist_avg_AMF ~ fert_status + site,
  data = meta_AMF_rare,
  method = "bray",
  permutations = how(blocks = meta_AMF_rare$site, nperm = 999),
  by = "terms"
)

# 3) Does fertilizer effect differ among sites? --------------------------------
# NOTE. Should we check for an interaction between site and fert_status? Usually, 
# interactions are justified when you hypothesized context-dependent effects. For
# example that fert_status impact site A but not site B. If you have no reason to
# expect that, then the interaction is not justified.

adonis2(
  dist_avg_AMF ~ fert_status * site,
  data = meta_AMF_rare,
  method = "bray",
  permutations = 999,
  by = "margin"
)

adonis2(
  dist_avg_AMF ~ fert_status * site,
  data = meta_AMF_rare,
  method = "bray",
  permutations = how(blocks = meta_AMF_rare$site, nperm = 999),
  by = "margin"
)

adonis2(
  dist_avg_AMF ~ fert_status * site,
  data = meta_AMF_rare,
  method = "bray",
  permutations = how(blocks = meta_AMF_rare$site, nperm = 999),
  by = "terms"
)

# 4) Does plot-level (beyond site) impact beta diversity? ----------------------
adonis2(
  dist_avg_AMF ~ site_plot,
  data = meta_AMF_rare,
  method = "bray",
  permutations = how(blocks = meta_AMF_rare$site, nperm = 999),
  by = "margin"
)

adonis2(
  dist_avg_AMF ~ site + site_plot,
  data = meta_AMF_rare,
  method = "bray",
  permutations = how(blocks = meta_AMF_rare$site, nperm = 999),
  by = "terms"
)

# 5) Does fertilizer effect vary among plots within sites? ---------------------
adonis2(
  dist_avg_AMF ~ site_plot * fert_status,
  data = meta_AMF_rare,
  method = "bray",
  permutations = how(blocks = meta_AMF_rare$site, nperm = 999),
  by = "margin"
)

# INTERPRETATION. This appears significant, we only have 2 samples per plot 
# (treatment/control). In this situation the interaction is essentially capturing
# plot-specific treatment differences. Given n = 2 per plot, this term absorbs 
# nearly all treatment contrast variability (No independent estimate of treatment 
# variance). So it is mathematically allowed, but biologically unstable. 
# I would not emphasize this interaction given the paired design and limited 
# replication. The interaction term is essentially capturing: Differences in 
# pairwise Bray–Curtis distances among plots. Not a replicated random slope.
# Not a properly estimated treatment heterogeneity.

# BETADISPER -------------------------------------------------------------------
# Check homogeneity of variances -----------------------------------------------

# The logic is identical to PERMANOVA — the permutation scheme must match the 
# scale at which the grouping variable exists:
# _Between sites → free permutations
# _Within sites (fert_status, plot) → restricted within sites

bd_fert <- betadisper(dist_avg_AMF, meta_AMF_rare$fert_status)
permutest(bd_fert, permutations = how(blocks = meta_AMF_rare$site, nperm = 999))
TukeyHSD(bd_fert)

bd_site <- betadisper(dist_avg_AMF, meta_AMF_rare$site)
permutest(bd_site, permutations = 999)
TukeyHSD(bd_site)

bd_plot <- betadisper(dist_avg_AMF, meta_AMF_rare$site_plot)
permutest(bd_plot, permutations = how(blocks = meta_AMF_rare$site, nperm = 999))
TukeyHSD(bd_plot)

# NOTE. anova() is ok, but for publication you want permutest() because:
# _It doesn't assume normality of distances.
# _It keeps your permutation scheme consistent with your PERMANOVA — a reviewer 
#    will notice if you restricted permutations in adonis2 but then used a free 
#    parametric test for betadisper.
# _It's the approach vegan's own documentation recommends for this kind of data.


# PLOTTING ---------------------------------------------------------------------
# community composition --------------------------------------------------------

# Run NMDS
set.seed(270226)

amf_nmds <- metaMDS(dist_avg_AMF, k = 2, trymax = 200)
amf_pcoa <- cmdscale(dist_avg_AMF, k = 2, eig = TRUE)

plot_ordination(
  ord = amf_nmds,
  meta = meta_AMF_rare,
  col_var = "site",
  shape_var = "fert_status"
)

plot_ordination(
  ord = amf_pcoa,
  meta = meta_AMF_rare,
  col_var = "site",
  shape_var = "fert_status", 
  ellipse = FALSE
)

plot_ordination(
  ord = amf_pcoa,
  meta = meta_AMF_rare,
  col_var = "site",
  shape_var = "fert_status", 
  ellipse = TRUE
)


# dispersion -------------------------------------------------------------------
plot_betadisper(dist_matrix = dist_avg_AMF, 
                grouping = meta_AMF_rare$fert_status)

plot_betadisper(dist_matrix = dist_avg_AMF, 
                grouping = meta_AMF_rare$site)

plot_betadisper(dist_matrix = dist_avg_AMF, 
                grouping = meta_AMF_rare$site_plot)


# Plotting results together

Figure_2_beta_long <-
ggarrange(
  plot_ordination(
    ord = amf_pcoa,
    meta = meta_AMF_rare,
    col_var = "site",
    shape_var = "fert_status",
    ellipse = FALSE,
    legend_inside = TRUE
  ) +
    labs(title = "PCoA") +
    scale_color_manual(values = palette_site),
  plot_betadisper(
        dist_matrix = dist_avg_AMF,
        grouping = meta_AMF_rare$fert_status,
        signif_label = "Treatment, F=1.006, p=0.243"
        ) +
        labs(title = "Within-treatment variance",
             y = NULL),
  plot_betadisper(
    dist_matrix = dist_avg_AMF,
    grouping = meta_AMF_rare$site,
    signif_label = "Site, F=6.88, p=0.001"
  ) +
    labs(title = "Within-site variance",
         y = NULL),
  #plot_betadisper(
  #  dist_matrix = dist_avg_AMF,
  #  grouping = meta_AMF_rare$site_plot,
  #  signif_label = "Plot, F=4.777, p=0.001"
  #) +
  #  labs(title = "Within-plot variance",
  #       y = NULL),
  ncol = 1,
  nrow = 3,
  align = "hv",
  heights =  c(1.4, 0.6, 0.8),
  labels = c("A", "B", "C")
)

# Using ggpubr -----------------------------------------------------------------
# ***** FIGURE 2 - Beta-diversity ***** ----------------------------------------
Figure_2_beta_long

ggsave(
  file.path(data_path, "results/Fig_2_betadiv.pdf"),
  plot = ggpubr::annotate_figure(
    Figure_2_beta_long,
    top = text_grob("BETA DIVERSITY", size = 12, face = "bold")
  ),
  device = "pdf"
)


ggsave(
  file.path(data_path, "results/Fig_S3_betadiv.pdf"),
  plot = ggpubr::annotate_figure(
    plot_betadisper(
        dist_matrix = dist_avg_AMF,
        grouping = meta_AMF_rare$site_plot,
        signif_label = "Plot, F=4.777, p=0.001"
      ) +
        labs(title = "Within- and between-plot variance",
             y = NULL),
    top = text_grob("BETA DIVERSITY", size = 12, face = "bold")
  ),
  device = "pdf"
)



# INTERPREATTION. NMDS space is different than PCoA space, that is why the spread
# (variances) of the samples within sites is different — the axes and distances are
# not the same thing.
# _ betadisper() calculates distances to centroids in PCoA space:
# _ Takes the Bray-Curtis distance matrix
# _ Runs a PCoA (cmdscale()) on it to convert distances into Euclidean coordinates\
#    in multivariate space
# _ Finds the centroid of each group in that PCoA space
# _ Calculates the Euclidean distance from each sample to its group centroid (or 
#   medians) in PCoA space. Those distances are what get plotted and tested.

# _ NMDS plot shows dispersion visually in NMDS space, which is a non-metric 
# rank-based reduction 

# Using patchwork --------------------------------------------------------------
pcoa_plot <- plot_ordination(
  ord = amf_pcoa,
  meta = meta_rare,
  col_var = "site",
  shape_var = "fert_status",
  ellipse = FALSE, 
  legend_inside = TRUE
) +
  labs(title = "PCoA")

btadisper_plot_fert <- 
  plot_betadisper(
    dist_matrix = bray_dist,
    grouping = meta_rare$fert_status,
    signif_label = "Treatment, F=1.006, p=0.243"
  ) +
  labs(title = "Within-treatment variance",
       y = NULL)

btadisper_plot_site <- 
  plot_betadisper(
    dist_matrix = bray_dist,
    grouping = meta_rare$site,
    signif_label = "Site, F=6.88, p=0.001"
  ) +
  labs(title = "Within-site variance",
       y = NULL)

btadisper_plot_rep <- 
  plot_betadisper(
    dist_matrix = bray_dist,
    grouping = meta_rare$site_plot,
    signif_label = "Plot, F=4.777, p=0.001"
  ) +
  labs(title = "Within-plot variance",
       y = NULL)

Figure_2_beta <-
  (pcoa_plot | (btadisper_plot_fert / btadisper_plot_site) | btadisper_plot_rep) +
  plot_layout(widths = c(1.5, 1, 1.3)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

Figure_2_beta

# ***** FIGURE 2 - Beta-diversity ***** -----------------------------------------

# Need to modify this if want to add the title!
ggsave(
  file.path(data_path, "results/Fig_2_betadiv_wide.pdf"),
  plot = grid.arrange(
    Figure_2_beta,
    top = text_grob("BETA DIVERSITY", size = 12, face = "bold")
  ),
  device = "pdf"
)



# **********************************************************************--------
 # **** 3. COMPOSITION AND STRUCTURE **** --------------------------

physeq_AMF_rare <-
  physeq_AMF_clean %>%
  rarefy_even_depth(rngseed = 26031, sample.size = rare_depth_cutoff) %>% 
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.) 

physeq_AMF_rare

# Condense at genus level
physeq_AMF_rare_Gen <-
  physeq_AMF_rare %>%
  transform_sample_counts(function(x) 100 * x / sum(x)) %>% 
  #subset_taxa(!Genus %in% c("Mortierella", "Jimgerdemannia", "Unclassified")) %>% 
  speedyseq::select_tax_table(Phylum,Class,Order,Family,Genus) %>% 
  phyloseq::tax_glom(., taxrank="Genus")

physeq_AMF_rare_Gen

as.matrix(physeq_AMF_rare_Gen@tax_table) %>% as.data.frame()
as.matrix(physeq_AMF_rare_Gen@otu_table) %>% as.data.frame()
as.matrix(physeq_AMF_rare_Gen@sam_data) %>% as.data.frame()

# Create a long format data.frame ----------------------------------------------
physeq_AMF_rare_Gen %>%
  psmelt() %>%
  arrange(Genus) %>%
  head()

# Plotting 
bar_charts_gen <-
  physeq_AMF_rare_Gen %>%
  psmelt() %>%
  mutate(Genus = fct_relevel(Genus, "Unclassified", after = Inf)) %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(color = "black", linewidth = 0.05) + # width=1 removes gaps between bars
  ggh4x::facet_nested(. ~ site + fert_status, 
               scales = "free_x", 
               space = "free_x", 
               switch = "x",
               nest_line = element_line(color = "black")) + # Optional line to join groups
  theme_classic() +
  theme(
    plot.title=element_text(size = 12, face="bold", hjust=0.5),
    plot.subtitle=element_text(size = 10),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    strip.placement = "outside", 
    strip.background = element_blank(),
    strip.text = element_text(size = 9, face = "bold"),
    panel.spacing = unit(0.3, "lines"), 
    legend.key.height = unit(0.3, "cm"),
    legend.key.width  = unit(0.3, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(face = "italic", size = 8)
  ) +
  scale_fill_manual(values = palette_taxa) +
  labs(title = "Genera relative abundance across sites",
    y = "DNA sequences (%)") +
  guides(
    fill = guide_legend(ncol = 1, override.aes = list(shape = 15)))


bar_charts_gen

# ***** FIGURE 3 - Genus Composition ***** -------------------------------------

ggsave(
  file.path(data_path, "results/Fig_3_bar_charts_gen_relab.pdf"),
  plot = grid.arrange(
    bar_charts_gen,
    top = text_grob("ARBUSCULAR MYCORRHIZAL FUNGI COMMUNITY COMPOSITION", 
                    size = 12, face = "bold")
  ),
  device = "pdf"
)

# Condense at family level
physeq_AMF_rare_Fam <-
  physeq_AMF_rare %>%
  transform_sample_counts(function(x) 100 * x / sum(x)) %>% 
  #subset_taxa(!Genus %in% c("Mortierella", "Jimgerdemannia", "Unclassified")) %>% 
  speedyseq::select_tax_table(Phylum,Class,Order,Family) %>% 
  phyloseq::tax_glom(., taxrank="Family")

physeq_AMF_rare_Fam %>%
  psmelt() %>%
  arrange(Family) %>%
  head()


bar_charts_fam <-
  physeq_AMF_rare_Fam %>%
  psmelt() %>%
  mutate(Family = fct_relevel(Family, "Unclassified", after = Inf)) %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", width = 1) + # width=1 removes gaps between bars
  ggh4x::facet_nested(. ~ site + fert_status, 
                      scales = "free_x", 
                      space = "free_x", 
                      switch = "x",
                      nest_line = element_line(color = "black")) + # Optional line to join groups
  theme_classic() +
  theme(
    plot.title=element_text(size = 12, face="bold", hjust=0.5),
    plot.subtitle=element_text(size = 10),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    strip.placement = "outside", 
    strip.background = element_blank(),
    strip.text = element_text(size = 9, face = "bold"),
    panel.spacing = unit(0.3, "lines"), 
    legend.key.height = unit(0.3, "cm"),
    legend.key.width  = unit(0.3, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(face = "italic", size = 8)
  ) +
  scale_fill_manual(values = c(palette_taxa[1:9], "black")) +
  labs(y = "DNA reads (%)") +
  guides(
    fill = guide_legend(ncol = 1, 
                        override.aes = list(shape = 15)))

bar_charts_fam

# ***** FIGURE 4 - Family Composition ***** -------------------------------------

ggsave(
  file.path(data_path, "results/Fig_3_bar_charts_fam_relab.pdf"),
  plot = grid.arrange(
    bar_charts_fam,
    top = text_grob("ARBUSCULAR MYCORRHIZAL FUNGI COMMUNITY COMPOSITION", 
                    size = 12, face = "bold")
  ),
  device = "pdf"
)


# **********************************************************************--------
# **** 5. DIFFERENTIAL ABUNDANCE **** ------------------------------------------

str(physeq_AMF_rare)

metadata_AMF_rare <- 
  as.data.frame(as.matrix(physeq_AMF_rare@sam_data)) %>% 
  mutate(across(.cols = 1:2, .fns = as.factor))

str(metadata_AMF_rare)

otutable_AMF_rare <- 
  as.data.frame(t(as.matrix(physeq_AMF_rare@otu_table)))

str(otutable_AMF_rare)

# Identify differentially abundant taxa  ---------------------------------------

# OTU level --------------------------------------------------------------------

table(metadata_AMF_rare$site, metadata_AMF_rare$fert_status)

# Test differential abundance between Fertilized vs Control, controlling for site

mlin2_fert_status <- 
  Maaslin2(
    cores = 4,
    output = "demo_output",
    input_data = otutable_AMF_rare,
    input_metadata = metadata_AMF_rare %>% dplyr::select(fert_status, site),
    fixed_effects = "fert_status",
    random_effects = "site",
    reference = "fert_status,Control", # no space between variable and level 
    normalization = "NONE",
    standardize = FALSE,
    plot_scatter = FALSE
  )

mlin2_fert_status$results %>%
  dplyr::filter(metadata == "fert_status", qval <= 0.05)

mlin2_fert_status$results %>%
  dplyr::filter(metadata == "site", qval <= 0.05)

# INTERPRETATION. Multiple testing penalty is brutal here. When testing ~2870 
# features BH correction inflates q-values. You need very strong, consistent 
# effects to pass q ≤ 0.05

# Look at raw p-values (not q-values)
mlin2_fert_status$results %>%
  dplyr::filter(metadata == "fert_status") %>%
  dplyr::arrange(pval) %>%
  head(20)

# Check effect sizes
mlin2_fert_status$results %>%
  dplyr::filter(metadata == "fert_status") %>%
  dplyr::arrange(desc(abs(coef))) %>%
  head(20)

# Genus level ------------------------------------------------------------------

# NOTE! I remove all Unclassified OTUs clusters at genus as they will inflate the 
# p-value correction and they are not informative. 

metadata_AMF_rare_gen <- 
  sample_data(
    physeq_AMF_rare_Gen %>% 
      subset_taxa(Genus != "Unclassified") %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.) 
  )%>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  mutate(across(.cols = 1:2, .fns = as.factor)) 

otutable_AMF_rare_gen <- 
  otu_table(
    physeq_AMF_rare_Gen %>% 
      subset_taxa(Genus != "Unclassified") %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.)
  ) %>% 
  t() %>% 
  as.matrix() %>% 
  as.data.frame() 

mlin2_fert_status_gen <- 
  Maaslin2(
    cores = 4,
    output = "demo_output",
    input_data = otutable_AMF_rare_gen,
    input_metadata = metadata_AMF_rare_gen %>% dplyr::select(fert_status, site),
    fixed_effects = "fert_status",
    random_effects = "site",
    reference = "fert_status,Control", # no space between variable and level 
    normalization = "NONE",
    standardize = FALSE,
    plot_scatter = FALSE
  )

mlin2_fert_status_gen$results %>%
  dplyr::filter(metadata == "fert_status", qval <= 0.05)

mlin2_fert_status_gen$results %>%
  dplyr::filter(metadata == "site", qval <= 0.05)

physeq_AMF_rare_Gen@tax_table %>% 
  as.data.frame() %>% 
  rownames_to_column("otu") %>%
  filter(otu %in% "Zotu3633")

# What about at family level?
metadata_AMF_rare_fam <- 
  sample_data(
    physeq_AMF_rare_Fam %>% 
      subset_taxa(Family != "Unclassified") %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.) 
  )%>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  mutate(across(.cols = 1:2, .fns = as.factor)) 

otutable_AMF_rare_fam <- 
  otu_table(
    physeq_AMF_rare_Fam %>% 
      subset_taxa(Family != "Unclassified") %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.)
  ) %>% 
  t() %>% 
  as.matrix() %>% 
  as.data.frame() 

mlin2_fert_status_fam <- 
  Maaslin2(
    cores = 4,
    output = "demo_output",
    input_data = otutable_AMF_rare_fam,
    input_metadata = metadata_AMF_rare_fam %>% dplyr::select(fert_status, site),
    fixed_effects = "fert_status",
    random_effects = "site",
    reference = "fert_status,Control", # no space between variable and level 
    normalization = "NONE",
    standardize = FALSE,
    plot_scatter = FALSE
  )

mlin2_fert_status_fam$results %>%
  dplyr::filter(metadata == "fert_status", qval <= 0.05)

mlin2_fert_status_fam$results %>%
  dplyr::filter(metadata == "site", qval <= 0.05)

physeq_AMF_rare_Fam@tax_table %>% 
  as.data.frame() %>% 
  rownames_to_column("otu") %>%
  filter(otu %in% "Zotu3633")

# RESUTLS. No significant differences between fertilized and control samples at 
# OTU and Family levels but Rhizophagus (Zotu3633) is higher in fertilized plots.

# Now, I can compare with a fixed-effect site model, no random effects! We can report 
# this if we explain what we did. For direct compariosn I can use the dataset that
# was agglomerated by Genus.

mlin2_site <- 
  Maaslin2(
    cores = 4,
    output = "demo_output",
    input_data = otutable_AMF_rare_gen,
    input_metadata = metadata_AMF_rare_gen %>% dplyr::select(fert_status, site),
    fixed_effects = c("fert_status", "site"),
    reference = "site,Escanaba",
    correction = "BH",
    normalization = "NONE",
    standardize = FALSE,
    plot_scatter = FALSE
  )

mlin2_site$results %>%
  dplyr::filter(metadata == "fert_status", qval <= 0.05)

# RESULTS. This confirm the result from above, the mixed effect model.

mlin2_site$results %>%
  dplyr::filter(metadata == "site", qval <= 0.05)

# RESULTS. This is consistent with the PERMANOVA result, site has a strong effect
# on community composition. These difference are related to the Escanaba site as 
# reference level. We should run it pairwise to know which sites differ from each 
# other and then use a stringent p value correction method.


# pairwise maaslin2 ------------------------------------------------------------
run_maaslin2_all_pairs <- function(physeq, 
                                   factor_variable, 
                                   qval_threshold = 0.05, 
                                   norm_method, 
                                   trans_method,
                                   ...) {
  
  # extract datasets
  input_metadata <- as.data.frame(as.matrix(physeq@sam_data)) %>% 
    dplyr::select(.data[[factor_variable]])
  
  input_data <- as.data.frame(t(as.matrix(physeq@otu_table)))
  
  # Get unique levels of the factor variable
  factor_levels <- unique(input_metadata[[factor_variable]])
  
  # Store results
  significant_results <- list()
  
  # Loop through each level and set it as the reference
  for (ref in factor_levels) {
    cat("\nRunning Maaslin2 with reference:", ref, "\n")
    
    result <- Maaslin2(
      cores = 6,
      input_data = input_data,
      input_metadata = input_metadata,
      output = tempfile(),  # Avoid creating permanent folders
      reference = c(factor_variable, ref),
      min_abundance = 0.001,
      min_prevalence = 0.1,
      correction = "BH",
      normalization = norm_method,
      transform = trans_method,
      standardize = FALSE,
      plot_scatter = FALSE,
      ...
    )
    
    # Extract significant results and add reference group
    sig_results <- result$results %>%
      filter(qval <= qval_threshold) %>%
      mutate(reference_group = ref)
    
    # Store results
    significant_results[[ref]] <- sig_results
  }
  
  # Combine all significant results into a single dataframe
  final_results <- bind_rows(significant_results)
  
  return(final_results)
}

physeq_AMF_rare_Gen@tax_table

mlin2_site <- 
  run_maaslin2_all_pairs(
    physeq = 
      physeq_AMF_rare_Gen %>% 
      subset_taxa(Genus != "Unclassified") %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.),
    factor_variable = "site",
    qval_threshold = 0.05, 
    norm_method = "NONE", 
    trans_method = "LOG") 

dim(mlin2_site)
head(mlin2_site)


mlin2_fert_status <- 
  run_maaslin2_all_pairs(
    physeq = 
      physeq_AMF_rare_Gen %>% 
      subset_taxa(Genus != "Unclassified") %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.),
    factor_variable = "fert_status",
    qval_threshold = 0.05, 
    norm_method = "NONE", 
    trans_method = "LOG") 

dim(mlin2_fert_status)
head(mlin2_fert_status)



# Plotting Maaslin2 ------------------------------------------------------------

# INTERPRETATION. The issue is the sign interpretation. Since MaAsLin2 coefficient
# s are always coef = value - reference = Escanaba - Lux Arbor

PlotMaaslin2 <- function(maaslin_results, 
                         physeq_object,
                         value_type = c("ln", "log2")) {
  
  value_type <- match.arg(value_type)
  
  # Process the data
  cleaned_results <- 
    maaslin_results %>%
    rowwise() %>%
    mutate(pair_raw = paste(value, reference_group, sep = " vs "),
           pair = paste(sort(c(value, reference_group)), collapse = "-")) %>%
    ungroup() %>%
    distinct(feature, pair, .keep_all = TRUE) %>%
    # Direction-aware label
    mutate(pair_label = paste0(value, " vs ", reference_group)) %>%
    left_join(as.data.frame(as.matrix(physeq_object@tax_table)) %>%
                rownames_to_column("feature"), by = "feature") %>%
    mutate(Genus_makrdown = str_c("*", str_trim(Genus), "*"),
           coef_log2   = coef / log(2),
           stderr_log2 = stderr / log(2))
  
  # Choose what to plot
  if (value_type == "ln") {
    cleaned_results <- cleaned_results %>%
      mutate(plot_value = coef,
             label_value = round(coef, 2)) %>% 
      na.omit()
    
    fill_label <- "ln fold change"
    
    print(range(cleaned_results$coef))
    
    fill_scale <- scale_fill_gradient2(
      name = "lFC",
      low = "blue", high = "red", mid = "#FFFFCC",
      na.value = "white",
      limits = c(-8, 8),
      breaks = c(-8, -4, 0, 4, 8)
    )
    
  } else if (value_type == "log2") {
    
    cleaned_results <- cleaned_results %>%
      mutate(plot_value = coef_log2,
             label_value = round(coef_log2, 2)) %>% 
      na.omit()
    
    fill_label <- "log2 fold change"
    
    print(range(cleaned_results$coef_log2))
    
    fill_scale <- scale_fill_gradient2(
      name = "l2FC",
      low = "blue", high = "red", mid = "#FFFFCC",
      na.value = "white",
      limits = c(-12.5, 12.5),
      breaks = c(-12.5, -6.25, 0, 6.25, 12.5)
    )
  }
  
  # Plot
  p <- 
    ggplot(cleaned_results, aes(x = pair_label, y = Genus_makrdown, fill = plot_value)) + 
    geom_tile() +
    geom_text(aes(label = label_value), size = 4, show.legend = FALSE) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "white", color = "white"),
      plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_markdown(size = 10, hjust = 0.5),
      axis.text.x = element_text(angle = 33, hjust = 1, vjust = 1, size = 10),
      axis.text.y = element_markdown(size = 10),
      axis.title = element_markdown(size = 10), 
      legend.title = element_text(size = 10, face = "bold"), 
      legend.text = element_text(size = 9),
      legend.key.height = unit(0.4, "cm"),  
      legend.key.width = unit(0.4, "cm"),
      legend.position = "right"
    ) +
    fill_scale
  
  return(p)
}


Fig_4_diff_abund_heat <-
  PlotMaaslin2(maaslin_results = mlin2_site %>% filter(qval <= 0.05),
               physeq_object = physeq_AMF_rare_Gen,
               value_type = "log2") +
  labs(
    title = "Differentially abundant genera",
    x = "Pairwise site comparison (test site vs reference site)",
    y = "Genus"
  ) +
  theme(
    plot.subtitle = ggtext::element_markdown(lineheight = 1.2)
  )

ggarrange(
  bar_charts_gen,
  Fig_4_diff_abund_heat,
  ncol = 1, 
  nrow = 2,
  heights = c(0.5, 1), 
  labels = c("A", "B")
  )


ggsave(
  file.path(data_path, "results/Fig_3_Compos_DiffAbund.pdf"),
  plot = grid.arrange(
    ggarrange(
      bar_charts_gen,
      Fig_4_diff_abund_heat,
      ncol = 1, 
      nrow = 2,
      heights = c(0.5, 1), 
      labels = c("A", "B")
    ),
    top = text_grob("ARBUSCULAR MYCORRHIZAL FUNGI COMMUNITY COMPOSITION", 
                    size = 12, face = "bold")
  ),
  device = "pdf"
)



# **********************************************************************--------
# **** 4. PHYLOGENETIC TREE **** ------------------------------------------

# https://yulab-smu.top/treedata-book/chapter10.html

palette_bestmatch <-
  c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77",
    "#283dff","#636bb7","#bfc5ff","#014443","#195637",
    "#117744","#60ffaf","#b7ffdb","#825121","#ea7f17",
    "#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d",
    "#5b5b19","#fcfc00","#ffff9e","#ffb7ef","#fa7efc",
    "#ae09ea","#521899","#1e0047")

# A. Generating a 99% OTUs tree ---------------------------------------------------
#Working with the RAXML tree for now
tree_raxml_AMF <- phy_tree(physeq_AMF_rare)

# Generating metadata for plotting ---------------------------------------------
melted_AMF <- 
  physeq_AMF_rare %>%
  psmelt() %>%
  arrange(Genus)

head(melted_AMF)
dim(melted_AMF)

dat1 <- melted_AMF %>%
  group_by(OTU, BestMatch, Genus, Phylum, Family) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus)) %>% 
  rename(label = OTU)

range(dat1$MeanAbundance)
range(sqrt(dat1$MeanAbundance))

dat2 <- melted_AMF %>%
  group_by(OTU, site) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  rename(ID = OTU, Sites = site)

dat3 <- melted_AMF %>%
  group_by(OTU, site) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  rename(ID = OTU, Sites = site)

# Build original plot ---------------------------------------------------------------
tree_AMF <- 
  tree_raxml_AMF %>% 
  ggtree(layout = "circular", branch.length = "none", size = 0.1, open.angle = 5)

tree_AMF <-
  tree_raxml_AMF %>% 
  ggtree(branch.length = "none", size = 0.15, open.angle = 5)

tree_AMF <- 
  tree_raxml_AMF %>% 
  ggtree(layout = "fan", size = 0.15, open.angle = 5)


# Layer 1: tip points — color = Genus, size = mean abundance
tree_AMF <- 
  tree_AMF %<+% dat1 +
  geom_tippoint(
    aes(fill = BestMatch, size = sqrt(MeanAbundance)),
    shape = 21, stroke = 0.1, alpha = 0.85) +
  scale_fill_manual(
    values = palette_bestmatch,
    guide  = guide_legend(
      keywidth = 0.5, keyheight = 0.5, order = 1,
      override.aes = list(shape = 22, size = 3)
    ),
    na.translate = FALSE
  ) +
  scale_size_continuous(
    range = c(0.1, 4),
    name  = "sqrt(Abundance)",
    guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 2,
                         override.aes = list(shape = 21))
  )

tree_AMF <- 
  tree_AMF +
  geom_tiplab2(
    aes(label = BestMatch),
    size = 1,
    align = TRUE,
    offset = 0.02
  )

tree_AMF 

# Layer 2: site abundance heatmap ring
tree_AMF <- 
  tree_AMF +
  new_scale_fill() +
  geom_fruit(
    data    = dat2,
    geom    = geom_tile,
    mapping = aes(y = ID, x = Sites, alpha = Abundance, fill = Sites),
    color   = "grey90", offset = 0.04, size = 0.02
  ) +
  scale_alpha_continuous(
    range = c(0, 1),
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 4)
  ) +
  scale_fill_manual(
    values = palette_site,
    guide  = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 3)
  )

# Layer 3: total abundance bar ring
tree_AMF <- 
  tree_AMF +
  new_scale_fill() +
  geom_fruit(
    data        = dat3,
    geom        = geom_bar,
    mapping     = aes(y = ID, x = TotalAbundance, fill = Sites),
    pwidth      = 0.38,
    orientation = "y",
    stat        = "identity"
  ) +
  scale_fill_manual(
    values = palette_site,
    guide  = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 3)
  )

# Theme 
tree_AMF <- 
  tree_AMF +
  geom_treescale(fontsize = 2, linesize = 0.3) +
  theme(
    legend.position   = c(0.93, 0.5),
    legend.background = element_rect(fill = NA),
    legend.title      = element_text(size = 6.5),
    legend.text       = element_text(size = 4.5),
    legend.spacing.y  = unit(0.02, "cm")
  )

tree_AMF

tree_AMF + layout_circular() + 
  theme(legend.position=c(0.05, 0.7)) 

tree_AMF + layout_rectangular() + 
  theme(legend.position=c(0.7, 0.2)) 

# B. Reducing tree complexity my agglomerating taxa ----------------------------

sample_data(physeq_AMF_rare)$fert_status_site <- 
  paste(sample_data(physeq_AMF_rare)$fert_status, 
        sample_data(physeq_AMF_rare)$site,
        sep="-")

physeq_AMF_rare@sam_data

merged_AMF <- merge_samples(physeq_AMF_rare, "fert_status_site")
merged_AMF@sam_data

merged_AMF <- 
  tax_glom(merged_AMF, taxrank = "BestMatch") %>% 
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.)

# Checks
merged_AMF@tax_table %>% as.data.frame()
merged_AMF@tax_table %>% as.data.frame() %>% subset(Genus == "Diversispora")
merged_AMF

# Melting 
melt_merged_AMF <- 
  psmelt(merged_AMF) %>% 
  separate(col = Sample, into = c("fert_status", "site"), sep = "-") %>% 
  mutate(val = sqrt(Abundance))

#melt_merged_AMF$Species
#separate(col = Sample, into = c("fert_status", "site"), sep = "-")
melt_merged_AMF %>% subset(Genus == "Diversispora")
head(melt_merged_AMF)
dim(melt_merged_AMF)

# NOTE. There are different Zotus classified with the same BestMatch name.
tree_raxml_AMF_merged <- phy_tree(merged_AMF)
tree_raxml_AMF_merged

# Generating metadata for plotting ---------------------------------------------
dat1_merged <- 
  melt_merged_AMF %>%
  group_by(OTU, BestMatch, Genus, Phylum, Family) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus)) %>% 
  rename(label = OTU)

range(dat1_merged$MeanAbundance)
range(sqrt(dat1_merged$MeanAbundance))

dat2_merged <-
  melt_merged_AMF %>%
  group_by(OTU, site) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  rename(ID = OTU, Sites = site)

dat2_merged

dat3_merged <- 
  melt_merged_AMF %>%
  group_by(OTU, site) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  rename(ID = OTU, Sites = site)

dat3_merged

dat4_merged <- 
  melt_merged_AMF %>%
  group_by(OTU, fert_status) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  rename(ID = OTU, N_addition = fert_status)

dat4_merged


# 1) Build plot version 1 ---------------------------------------------------------
tree_AMF_merged <- 
  tree_raxml_AMF_merged %>% 
  ggtree(layout = "circular", branch.length = "none", size = 0.15, open.angle = 5)

tree_AMF_merged <- 
  tree_raxml_AMF_merged %>% 
  ggtree(layout = "fan", size = 0.15, open.angle = 5)

tree_AMF_merged

# Layer 1: tip points — color = Genus, size = mean abundance
tree_AMF_merged <- 
  tree_AMF_merged %<+% dat1_merged +
  geom_tippoint(
    aes(fill = BestMatch, size = sqrt(MeanAbundance)),
    shape = 21, stroke = 0.1, alpha = 0.85) +
  scale_fill_manual(
    values = palette_bestmatch,
    guide  = guide_legend(
      keywidth = 0.5, keyheight = 0.5, order = 1,
      override.aes = list(shape = 21, size = 3)
    ),
    na.translate = FALSE
  ) +
  scale_size_continuous(
    range = c(0.1, 4),
    name  = "sqrt(Abundance)",
    guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 2,
                         override.aes = list(shape = 21))
  )

tree_AMF_merged <- 
  tree_AMF_merged +
  geom_tiplab2(
    aes(label = BestMatch),
    size = 3,
    align = TRUE,
    offset = 0.2
  )

tree_AMF_merged

# Layer 2: site abundance heatmap ring
tree_AMF_merged <- 
  tree_AMF_merged +
  new_scale_fill() +
  geom_fruit(
    data    = dat2_merged,
    geom    = geom_tile,
    mapping = aes(y = ID, x = Sites, alpha = Abundance, fill = Sites),
    color   = "grey90", offset = 0.04, size = 0.02
  ) +
  scale_alpha_continuous(
    range = c(0, 1),
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 4)
  ) +
  scale_fill_manual(
    values = palette_site,
    guide  = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 3)
  )

# Layer 3: total abundance for fert_status
tree_AMF_merged <- 
  tree_AMF_merged +
  new_scale_fill() +
  geom_fruit(
    data    = dat4_merged,
    geom    = geom_tile,
    mapping = aes(y = ID, x = N_addition, alpha = Abundance, fill = N_addition),
    color   = "grey90", offset = 0.04, size = 0.02
  ) +
  scale_alpha_continuous(
    range = c(0, 1),
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 4)
  ) +
  scale_fill_manual(
    values = palette_site,
    guide  = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 3)
  )

# Layer 4: total abundance for site
tree_AMF_merged <- 
  tree_AMF_merged +
  new_scale_fill() +
  geom_fruit(
    data        = dat3_merged,
    geom        = geom_bar,
    mapping     = aes(y = ID, x = TotalAbundance, fill = Sites),
    pwidth      = 0.38,
    orientation = "y",
    stat        = "identity"
  ) +
  scale_fill_manual(
    values = palette_site,
    guide  = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 3)
  )

# Theme 
tree_AMF_merged <- 
  tree_AMF_merged +
  geom_treescale(fontsize = 2, linesize = 0.3) +
  theme(
    legend.position   = c(0.93, 0.5),
    legend.background = element_rect(fill = NA),
    legend.title      = element_text(size = 6.5),
    legend.text       = element_text(size = 4.5),
    legend.spacing.y  = unit(0.02, "cm")
  )

tree_AMF_merged


tree_AMF_merged <- tree_AMF_merged + 
  geom_tiplab(
    aes(label = BestMatch), 
    size = 2.5,           # Keep it small for 1,700 tips!
    offset = 0.5,         # Push it slightly past the tippoint
    linetype = "dotted",   # "33" or "dotted"
    linewidth = 0.2,       # Keep it subtle
    color = "grey50", 
    geom = "text",        # Use 'text' to avoid the dotted line logic here
    align = TRUE         # Set to FALSE so it stays next to the point
  )


tree_AMF_merged

tree_AMF_merged + layout_circular() + 
  theme(legend.position=c(0.05, 0.7)) 

tree_AMF_merged + layout_rectangular() + 
  theme(legend.position=c(0.05, 0.7)) 

# Imporvements -----------------------------------------------------------------
# https://www.earlham.ac.uk/articles/plotting-phylogenetic-trees-r-alternating-clade-highlights

tree_AMF <- ggtree(tree_raxml_AMF) %<+% dat1

genera_nodes <- 
  dat1 %>%
  group_by(BestMatch) %>%
  filter(n() > 1) %>% # Only collapse if there's more than one tip
  summarize(node = ggtree::MRCA(tree_raxml_AMF, label))

genera_nodes

for(n in genera_nodes$node) {
  tree_AMF <- collapse(tree_AMF, node = n, mode = 'mixed', alpha = 0.4)
}

tree_AMF + 
  geom_tippoint(aes(fill = BestMatch, size = sqrt(MeanAbundance)), 
                shape = 21, stroke = 0.1) +
  geom_tiplab(aes(subset = !is.na(Genus)), offset = 0.2) + # Adjust as needed
  scale_fill_manual(values = palette_bestmatch)

tree_AMF + 
  geom_tippoint(aes(fill = BestMatch, size = sqrt(MeanAbundance)), 
                shape = 21, stroke = 0.1) +
  geom_tiplab(aes(subset = !is.na(Genus)), offset = 0.2) + # Adjust as needed
  scale_fill_manual(values = palette_bestmatch)

# This is Goood! 
# The Better Approach: geom_strip()

genera_nodes <- 
  dat1 %>%
  filter(label %in% tree_raxml_AMF$tip.label) %>%
  group_by(BestMatch) %>%
  summarize(
    t1 = first(label),
    t2 = last(label),
    .groups = "drop"
  ) %>%
  filter(t1 != t2) # Only draw strips for groups with more than 1 tip

genera_nodes

tree_AMF <- 
  ggtree(tree_raxml_AMF) %<+% dat1 +
  geom_tippoint(aes(fill = BestMatch, size = sqrt(MeanAbundance)), 
                shape = 21, stroke = 0.1, alpha = 0.8) +
  scale_fill_manual(values = palette_bestmatch)

for(i in 1:nrow(genera_nodes)) {
  tree_AMF <- 
    tree_AMF + 
    geom_strip(
    taxa1 = genera_nodes$t1[i], 
    taxa2 = genera_nodes$t2[i], 
    label = genera_nodes$BestMatch[i],
    offset = 0.1,       # Adjust this value based on your tree's branch lengths
    barsize = 1.5, 
    extend = 0.2, 
    fontsize = 3,
    offset.text = 0.1
  )
}

tree_AMF

tree_AMF + scale_fill_manual(
  values = palette_bestmatch,
  guide = guide_legend(
    title = "Taxa",
    ncol = 1,               # Force single column
    byrow = TRUE,
    override.aes = list(
      shape = 22,           # 22 is a square with a border
      size = 4,             # Smaller square size
      stroke = 0.2          # Thin border for the square
    )
  ),
  na.translate = FALSE
) +
  theme(
    legend.key.height = unit(0.4, "cm"), # Reduces vertical space between items
    legend.text = element_text(size = 8), # Smaller font size
    legend.title = element_text(size = 9, face = "bold")
  ) +
  xlim(0, 4) 

# Ttsting another option

tree_AMF <- 
  ggtree(tree_raxml_AMF, linewidth=0.5) %<+% dat1 +
  geom_tiplab(aes(label=label), size=2)

tree_AMF

#Make dataframe for clade nodes
genera_nodes <- data.frame(
  clade=unique(dat1$BestMatch),
  node=NA
)

genera_nodes

#Find the most recent common ancestor for each clade
for (i in 1:length(genera_nodes$clade)) {
  
  genera_nodes$node[i] <- MRCA(
    tree_raxml_AMF,
    dat1$label[dat1$BestMatch == genera_nodes$clade[i]]
  )
  
}

tree_raxml_AMF

tree_AMF <- 
  ggtree(tree_raxml_AMF, linetype=NA) %<+% dat1 +
  geom_highlight(data=genera_nodes, 
                 aes(node=node, fill=clade),
                 alpha=1,
                 align="right",
                 extend=0.1,
                 show.legend=FALSE) +
  geom_tree(linewidth=0.5) +
  geom_tiplab(aes(label=BestMatch), size=2)

tree_AMF



  
# 2) Build plot version 2 ------------------------------------------------------

tree_AMF_merged <- 
  merged_AMF %>% 
  ggtree(layout="fan", open.angle=1) + 
  geom_tippoint(mapping=aes(color = BestMatch, size = Abundance),
                show.legend=FALSE) +
  scale_color_manual(
    values = palette_bestmatch)

tree_AMF_merged

#tree_AMF_merged <- rotate_tree(tree_AMF_merged, -90)

tree_AMF_merged <- 
  tree_AMF_merged +
  geom_fruit(
    data=melt_merged_AMF %>% dplyr::select(OTU, val),
    geom=geom_boxplot,
    mapping = aes(
      y=OTU,
      x=val,
      group=OTU,
      fill=BestMatch,
    ),
    size=.2,
    outlier.size=0.5,
    outlier.stroke=0.08,
    outlier.shape=21,
    axis.params=list(
      axis       = "x",
      text.size  = 1.8,
      hjust      = 1,
      vjust      = 0.5,
      nbreak     = 3,
    ),
    grid.params=list()
  ) 


tree_AMF_merged <- 
  tree_AMF_merged +
  scale_fill_manual(
    values = palette_bestmatch,
    guide  = guide_legend(
      keywidth = 0.5, 
      keyheight = 0.5, 
      ncol = 1,
      override.aes = list(shape = 21, size = 6)
    ),
    na.translate = FALSE
  ) +
  theme(
    legend.title=element_text(size=9), 
    legend.text=element_text(size=7) 
  )


tree_AMF_merged <- 
  tree_AMF_merged + 
  geom_tiplab(
    aes(label = BestMatch), 
    size = 2.5,           # Keep it small for 1,700 tips!
    offset = 0.5,         # Push it slightly past the tippoint
    linetype = "dotted",   # "33" or "dotted"
    linewidth = 0.2,       # Keep it subtle
    color = "grey50", 
    geom = "text",        # Use 'text' to avoid the dotted line logic here
    align = TRUE         # Set to FALSE so it stays next to the point
  )

tree_AMF_merged

# Adding side abudances
tree_AMF_merged <- 
  tree_AMF_merged +
  new_scale_fill() +
  geom_fruit(
    data    = dat3_merged,
    geom    = geom_tile,
    mapping = aes(y = ID, x = Sites, alpha = val, fill = Sites),
    color   = "grey90", offset = 0.04, size = 0.02
  ) +
  scale_alpha_continuous(
    range = c(0, 1),
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 4)
  ) +
  scale_fill_manual(
    values = palette_site,
    guide  = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 3)
  )

tree_AMF_merged

# **********************************************************************--------
# B. Generating a 97% OTUs tree ------------------------------------------------

tree_raxml_97 <- 
  read.tree("phylogeny/otus_97_filtered_mafft_trim_spprt.raxml.support") 
str(tree_raxml_97)
ggtree::ggtree(tree_raxml_97)

# Rooting the tree to the outgroup ---------------------------------------------
tree_raxml_97_rooted <- root(tree_raxml_97, outgroup = "Zotu6503", resolve.root = TRUE)
ggtree::ggtree(tree_raxml_97_rooted)

tree_iqtree2_97 <- read.tree("phylogeny/otus_97_filtered_mafft_trim_iq2.treefile")
str(tree_iqtree2_97)
ggtree::ggtree(tree_iqtree2_97)

# Generate the phyloseq object -------------------------------------------------
# otutable
otutable_97 <- 
  read.delim(file.path(data_path, "datasets/otutab_97.txt"), row.names = 1)

# Taxonomy
taxonomy_97 <- 
  extract_blasTAX(tax_path = file.path(data_path, "datasets/taxonomy_blast_97.txt"),
                namemap_path = file.path(data_path, "datasets/name_mapping_97.txt")) %>% 
  dplyr::select(
    "Zotu", "Kingdom" ,"Phylum", "Class","Order","Family","Genus","Species") %>%
  filter(Phylum %in% "Mucoromycota") %>% 
  FinalizeTaxonomy()

# I will use Mortierellomycetes as outgroup to root the tree?
taxonomy_97 %>% subset(Class %in% c("Mortierellomycetes", "Endogonomycetes"))

# Fixing classificaiton inconsistencies
CheckTaxonomyConsistency(taxonomy_97, 
                         return_long=TRUE)

CheckTaxonomyConsistency(taxonomy_97, 
                         return_long=FALSE)

taxonomy_97

taxonomy_97 <- 
  taxonomy_97 %>%
  mutate(Family = if_else(Family == "Archaeosporaceae", "Ambisporaceae", Family))

# OTUs
zotu_97 <- readDNAStringSet(file.path(data_path,"datasets/otus_97.fasta"), 
                            format="fasta", seek.first.rec=TRUE, use.names=TRUE)

zotu_97_filtered <- zotu_97[ !names(zotu_97) %in% zotus_to_remove ]


# Removing only Jimgerdemannia
zotu_97_filtered <- zotu_97[ !names(zotu_97) %in% "Zotu11542" ]

# Phyloseq
physeq_97 <- generate_phyloseq(otu=otutable_97,
                               metadata=metadata_99,
                               taxonomy=taxonomy_97 %>% column_to_rownames("Zotu"),
                               sequences=zotu_97_filtered) # I want to keep the outgroups
physeq_97
physeq_97@sam_data

# Adding phylotree
physeq_97 <-
  phyloseq(
    physeq_97@otu_table,
    physeq_97@sam_data,
    physeq_97@tax_table,
    physeq_97@refseq,
    phy_tree(tree_raxml_97_rooted)
  ) %>% 
  phyloseq::prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  phyloseq::prune_samples(sample_sums(x=.) > 0, x =.)

physeq_97
print(sample_data(physeq_97), n=113)
sort(sample_sums(physeq_97))

# Remove control samples and Rarefaction (single round in this case) 
physeq_97_rare <-
  physeq_97 %>%
  subset_samples(!site %in% "Control") %>% 
  rarefy_even_depth(rngseed = 260414, sample.size = 6434) %>% 
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.) 

physeq_97_rare
print(physeq_97_rare@sam_data, n=109)

# Working with the RAXML tree for now
tree_raxml_AMF_97 <- phy_tree(physeq_97_rare)
ggtree::ggtree(tree_raxml_AMF_97)

phyloseq::plot_tree(physeq_97_rare, "treeonly", label.tips = "taxa_names", text.size=2) + 
  coord_polar(theta = "y") 

# Generating metadata for the phyloegentic tree --------------------------------
melted_AMF_97 <- 
  physeq_97_rare %>%
  psmelt() %>%
  arrange(Genus)

head(melted_AMF_97)
dim(melted_AMF_97)

dat1_97 <- 
  melted_AMF_97 %>%
  group_by(OTU, BestMatch, Genus, Phylum, Family) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus)) %>% 
  rename(label = OTU) %>% 
  mutate(OTU = paste(label, BestMatch, sep="-"))

range(dat1_97$MeanAbundance)
range(sqrt(dat1_97$MeanAbundance))

write.csv(x= dat1_97, file= file.path(data_path, "datasets/tree_metadata_97.csv"))
dat1_97_mod <- read.csv(file=file.path(data_path, "datasets/tree_metadata_97_mod_SL.csv"))

dat2_97 <- 
  melted_AMF_97 %>% 
  group_by(OTU, site) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  rename(ID = OTU, Sites = site)

dat3_97 <- 
  melted_AMF_97 %>%
  group_by(OTU, site) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  rename(ID = OTU, Sites = site)

# Grouping to clades -----------------------------------------------------------

# Make dataframe for clade nodes
genera_nodes <- data.frame(
  clade=unique(dat1_97$BestMatch),
  node=NA
)

# Find the most recent common ancestor for each clade
for (i in 1:length(genera_nodes$clade)) {
  
  genera_nodes$node[i] <- MRCA(
    tree_raxml_AMF_97,
    dat1_97$label[dat1_97$BestMatch == genera_nodes$clade[i]]
  )
  
}

# Add column with alternating binary value. This is based on the ggtree data
genera_nodes <-
  genera_nodes[match(
    tree_AMF_97$data %>%
      filter(isTip == "TRUE") %>%
      arrange(y) %>%
      pull(BestMatch) %>%
      unique(),
    genera_nodes$clade
  ), ] %>% 
  dplyr::mutate(highlight = rep(c(0,1),
                         length.out=length(genera_nodes$clade)))

genera_nodes

# Color palette for BestMatch --------------------------------------------------
palette_bestmatch <- c(
  # --- GLOMERALES (Deep Purples to Sky Blues) ---
  # Core Glomus & Rhizophagus (Purples)
  "Glomus"                = "#781156",
  "Glomus macrocarpum"     = "#A51876",
  "Glomus tetrastratosum"  = "#D21E96",
  "Rhizophagus"           = "#E43FAD",
  "Funneliformis"         = "#EA6CC0",
  "Septoglomus"           = "#F098D3",
  "uncultured Glomus"     = "#C2629E", # Moved to Glomerales
  "Glomeraceae"           = "#4D1136", # Darker version for higher rank
  
  # Dominikia & Silvaspora Clade (Deep Blues)
  "Dominikia"             = "#114578",
  "Microkamienskia"       = "#185EA5",
  "Microdominikia litorea"= "#1E78D2",
  "Silvaspora"            = "#3F91E4",
  "Oehlia"                = "#6CABEA",
  
  # --- DIVERSISPORALES (Teals and Seafoams) ---
  "Diversispora"            = "#117878",
  "Diversispora versiformis"= "#18A5A5",
  "Diversisporaceae"        = "#083B3B", # Higher rank
  "Acaulospora"             = "#3FE4E4",
  "Acaulospora brasiliensis" = "#6CEAEA",
  "Dentiscutata heterogama" = "#98F0F0",
  "Gigaspora"               = "#117845", # Branching into green-teals
  "Gigaspora margarita"     = "#18A55E",
  "Gigaspora rosea"         = "#1ED278",
  "Scutellospora"           = "#3FE491",
  "Cetraspora gilmorei"     = "#6CEAAB",
  
  # --- ARCHAEOSPORALES (Earth Greens) ---
  "Archaeospora"            = "#4B5D16",
  "Archaeospora trappei"    = "#787811",
  "Ambispora"               = "#A5A518",
  "Archaeosporaceae"        = "#2E330D", # Higher rank
  "Archaeosporales"         = "#D2D21E", # Higher rank
  "uncultured Archaeosporales" = "#E4E43F", # Moved here
  
  # --- PARAGLOMERALES (Warm Oranges & Browns) ---
  "Paraglomus"              = "#A55E18",
  "Paraglomus laccatum"     = "#D2781E",
  "Paraglomus brasilianum"  = "#E4913F",
  "Paraglomerales"          = "#784511", # Higher rank
  "uncultured Paraglomus"   = "#EAAB6C", # Moved here
  
  # --- ENTROPHOSPORALES & OTHERS (Reds and Pinks) ---
  "Entrophospora"           = "#781122",
  "Entrophospora claroidea" = "#A5182F",
  "Entrophospora drummondii" = "#D21E2C",
  "Complexispora"           = "#E43F5B",
  "Glomeromycetes"          = "#EA6C81", # Broadest group in Red
  "uncultured Glomeromycotina" = "#F0C498"  # Pale beige for the unknown
)

# "#000000","#242424","#484848","#6D6D6D","#919191","#B6B6B6"

# Build tree plot ---------------------------------------------------------------

tree_AMF_97 <- 
  tree_raxml_97_rooted %>% 
  ggtree(layout = "circular", branch.length = "none", size = 0.1, open.angle = 5) %<+% dat1_97 

tree_AMF_97 +
  geom_tippoint(
    aes(fill = BestMatch, size = sqrt(MeanAbundance)),
    shape = 21, stroke = 0.1, alpha = 0.85) +
  geom_tiplab(aes(label=OTU), size=2) +
  scale_fill_manual(values = bestmatch_colors_97, na.translate = FALSE) +
  geom_text2(aes(subset = !isTip, label = node), size = 2, color = "black") +
  theme(legend.position = "none")

tree_AMF_97 <- 
  tree_AMF_97 +
  geom_tippoint(
    aes(fill = BestMatch, size = sqrt(MeanAbundance)),
    shape = 21, stroke = 0.1, alpha = 0.85) +
  scale_fill_manual(values = bestmatch_colors_97, na.translate = FALSE) +
  scale_size_continuous(range = c(0.331801/5, 22.456278/5),
                        name = "sqrt(Read Abundance)") +
  geom_highlight(data=genera_nodes, 
                 aes(node=node, fill=clade),
                 alpha=0.2,
                 align="right",
                 extend=0.1,
                 show.legend = FALSE) + # try with TRUE as well
  geom_tree(linewidth=0.3) +
  geom_text2(aes(subset = !isTip, label = node), size = 2, color = "black") +
  geom_tiplab(aes(label=OTU), size=2) +
  geom_cladelab(data=genera_nodes,
                mapping=aes(node=node, label=clade, color=clade),
                parse = TRUE,
                fontsize=2,
                align="TRUE",
                angle="auto",
                offset=6,
                offset.text=0.05, 
                show.legend = FALSE) +
    guides(fill = guide_legend(ncol=1, keywidth = 0.2, keyheight = 0.2, order = 1,
                             override.aes = list(shape = 22, size = 3)),
         shape = guide_legend(ncol=1, keywidth = 0.1, keyheight = 0.1, order = 2,
                             override.aes = list(shape = 21, size = 3))) 

tree_AMF_97


# Plotting all together
ggtree(tree_raxml_97_rooted, layout="circular", linetype=NA) %<+% dat1_97 +
  geom_highlight(data=genera_nodes, 
                 aes(node=node, fill=as.factor(highlight)),
                 alpha=1,
                 align="right",
                 extend=0.04,
                 show.legend=FALSE) +
  geom_cladelab(data=genera_nodes,
                mapping=aes(node=node, label=clade),
                fontsize=2,
                align="TRUE",
                angle="auto",
                offset=0.04,
                offset.text=0.01) +
  geom_tree(linewidth=0.3) +
  #geom_tippoint() +
  geom_tippoint(aes(size = sqrt(MeanAbundance)), 
                shape = 16, color = "darkred", stroke = 0.2, alpha = 0.85, 
                show.legend=FALSE) +
  #xlim(0, 0.35) +
  scale_fill_manual(values=c("#F5F5F5", "#ECECEC")) +
  scale_size_continuous(range = c(0.331801/5, 22.456278/5),
                        name = "sqrt(Read Abundance)") 


# Plotting all together
ggtree(tree_raxml_97_rooted, layout="circular", linetype=NA) %<+% dat1_97 +
  geom_highlight(data=genera_nodes, 
                 aes(node=node, fill=clade),
                 alpha=1,
                 align="right",
                 extend=0.04,
                 show.legend=FALSE) +
  geom_cladelab(data=genera_nodes,
                mapping=aes(node=node, label=clade),
                fontsize=2,
                align="TRUE",
                angle="auto",
                offset=0.04,
                offset.text=0.01) +
  geom_tree(linewidth=0.3) +
  #geom_tippoint() +
  geom_tiplab(aes(fill=BestMatch), size=2, show.legend = FALSE) 



# Need to modufy this if want to add the title!
ggsave(
  file.path(data_path, "results/Fig_2_betadiv_wide.pdf"),
  plot = grid.arrange(
    Figure_2_beta,
    top = text_grob("BETA DIVERSITY", size = 12, face = "bold")
  ),
  device = "pdf"
)




tree_AMF_97 <- 
  tree_raxml_AMF_97 %>% 
  ggtree(layout = "circular", branch.length = "none", size = 0.1, open.angle = 5) 

tree_AMF_97

tree_AMF_97 +
  geom_text2(aes(subset = !isTip, label = node), size = 2, color = "black")

tree_AMF_97 <- 
  tree_AMF_97 %<+% dat1_97 +
  geom_tippoint(
    aes(fill = BestMatch, size = sqrt(MeanAbundance)),
    shape = 21, stroke = 0.1, alpha = 0.85) +
  scale_fill_manual(
    values = bestmatch_colors_97,
    guide  = guide_legend(
      keywidth = 0.5, keyheight = 0.5, order = 1,
      override.aes = list(shape = 21, size = 3)
    ),
    na.translate = FALSE
  ) +
  scale_size_continuous(
    range = c(0.1, 4),
    name  = "sqrt(Read Abundance)",
    guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 2,
                         override.aes = list(shape = 21))
  )

tree_AMF_97 <- 
  tree_AMF_97 +
  geom_tiplab2(
    aes(label = BestMatch),
    size = 2,
    align = TRUE,
    offset = 0.02
  )

tree_AMF_97 + theme(legend.position="none")


# Test tree 1 (warnings aren't harmful)














# Improvements -----------------------------------------------------------------

tree_AMF_97 <- 
  ggtree(tree_raxml_97_rooted, linewidth=0.5) %<+% dat1_97 +
  geom_tiplab(aes(label=label), size=2)

tree_AMF_97

#Make dataframe for clade nodes
genera_nodes <- data.frame(
  clade=unique(dat1_97$BestMatch),
  node=NA
)

genera_nodes

#Find the most recent common ancestor for each clade
for (i in 1:length(genera_nodes$clade)) {
  
  genera_nodes$node[i] <- MRCA(
    tree_raxml_97_rooted,
    dat1_97$label[dat1_97$BestMatch == genera_nodes$clade[i]]
  )
  
}

tree_raxml_97_rooted

tree_AMF_97 <- 
  ggtree(tree_raxml_97_rooted, linetype=NA) %<+% dat1_97 +
  geom_highlight(data=genera_nodes, 
                 aes(node=node, fill=clade),
                 alpha=1,
                 align="right",
                 extend=0.1,
                 show.legend=FALSE) + # try with TRUE as well
  geom_tree(linewidth=0.5) +
  geom_tiplab(aes(label=BestMatch), size=2)

tree_AMF_97 + xlim(0, 3) 


tree_AMF_97 + layout_circular() + 
  theme(legend.position=c(0.05, 0.7)) 

# Alternating highlight colors

genera_nodes <- genera_nodes[match(tree_AMF_97$data %>%
                               filter(isTip == "TRUE") %>%
                               arrange(y) %>%
                               pull(BestMatch) %>%
                               unique(),
                               genera_nodes$clade),]
#Add column with alternating binary value
genera_nodes$highlight <- rep(c(0,1),
                           length.out=length(genera_nodes$clade))
genera_nodes


#Add highlights
tree_AMF_97 <- 
  ggtree(tree_raxml_97_rooted, linetype=NA) %<+% dat1_97 +
  geom_highlight(data=genera_nodes, 
                 aes(node=node, fill=as.factor(highlight)),
                 alpha=1,
                 align="right",
                 extend=0.1,
                 show.legend=FALSE) +
  geom_tree(linewidth=0.5) +
  geom_tiplab(aes(label=BestMatch), size=2) +
  scale_fill_manual(values=c("#F5F5F5", "#ECECEC"))

tree_AMF_97

#Add clade labels
tree_AMF_97 +
  geom_cladelab(data=genera_nodes,
                mapping=aes(node=node, label=clade),
                fontsize=2,
                align=TRUE,
                offset=0.1,
                offset.text=0.02)


tree_raxml_97_rooted %>% 
  ggtree(layout = "circular", branch.length = "none", size = 0.1, open.angle = 5)































tree_AMF_merged <- tree_AMF_merged + 
  geom_tiplab(
    aes(label = BestMatch),       # We don't want text, just the line
    size = 2,
    linetype = "dotted",   # "33" or "dotted"
    linewidth = 0.2,       # Keep it subtle
    color = "grey50", 
    align = TRUE,          # This draws the line from tip to the alignment margin
    offset = -1             # Adjust if there is a gap
  )



# Adding geneus labels 
genus_nodes <- dat1 %>%
  filter(!is.na(Taxa)) %>%                    # genus-level only
  group_by(Taxa) %>%
  summarise(tips = list(label), .groups = "drop") %>%
  mutate(
    ntips = lengths(tips),
    node  = map_int(tips, function(t) {
      t <- unlist(t)
      if (length(t) > 1) ape::getMRCA(tree_raxml_AMF, t)
      else which(tree_raxml_AMF$tip.label == t)
    })
  )

# diagnostic — check for genera with MRCA near the root
root_node <- length(tree_raxml_AMF$tip.label) + 1
genus_nodes %>% dplyr::select(Taxa, ntips, node) %>% 
  arrange(node) %>% 
  print(n = Inf)

# optional: drop genera whose MRCA is suspiciously close to root
# adjust threshold after inspecting the diagnostic above
genus_nodes_plot <- genus_nodes %>% filter(node > root_node + 10)

  geom_cladelab(
  data       = genus_nodes_plot,
  mapping    = aes(node = node, label = Taxa),
  hjust      = 0.5,   angle      = "auto",
  barsize    = NA,    horizontal = FALSE,
  fontsize   = 1.4,   fontface   = "italic"
) +



# Get the MRCA node for each genus group
genus_nodes <- dat1 %>%
  group_by(Genus) %>%
  summarise(tips = list(label), .groups = "drop") %>%
  filter(Genus != "Unclassified") %>%
  rowwise() %>%
  mutate(
    node = ifelse(
      length(tips) > 1,
      ape::getMRCA(tree_raxml_AMF, unlist(tips)),  # use ape:: explicitly + unlist()
      which(tree_raxml_AMF$tip.label == unlist(tips)[[1]])
    )
  ) %>%
  ungroup()

genus_nodes

nodedf <- data.frame(node = genus_nodes$node)
labdf  <- data.frame(
  node  = genus_nodes$node,
  label = genus_nodes$Genus
)


physeq_AMF_clean
physeq_AMF_rare
physeq_AMF_rare_Gen

tree_raxml_AMF <- phy_tree(physeq_AMF_rare)
str(tree_raxml_AMF)

# Generating metadata for tree 
# Generating abundance, taxonomy, metadata table for the tree plot 

melted_AMF <- 
  physeq_AMF_rare %>%
  psmelt() %>%
  arrange(Genus)

head(melted_AMF)
dim(melted_AMF)

any(physeq_AMF_rare@tax_table == "Mortierellaceae")
is.na(physeq_AMF_rare@tax_table) 

# Create the base taxonomy mapping (1 row per OTU)
taxonomy_4_tree <- 
  melted_AMF %>%
  dplyr::select(OTU, Family, Genus, BestMatch) %>%
  distinct()

unique(taxonomy_4_tree$BestMatch)

# Calculate Total Abundance and Frequencies
tree_metadata <- 
  melted_AMF %>%
  mutate(present = ifelse(Abundance > 0, 1, 0)) %>%
  group_by(OTU) %>%
  dplyr::summarise(
    Abundance = sum(Abundance),
    Fertilized = sum(present[fert_status == "Fertilized"]),
    Control    = sum(present[fert_status == "Control"]),
    `Lux Arbor` = sum(present[site == "Lux Arbor"]),
    `Lake City`   = sum(present[site == "Lake City"]),
    Escanaba   = sum(present[site == "Escanaba"]),
    Rhinelander = sum(present[site == "Rhinelander"]),
    Hancock    = sum(present[site == "Hancock"])
  ) %>%
  left_join(taxonomy_4_tree, by = "OTU") %>%
  dplyr::select(OTU, Family, Genus, BestMatch, Abundance, everything())

tree_metadata

unique(tree_metadata$Genus)
unique(tree_metadata$BestMatch)

palette_bestmatch <-
  c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77",
    "#283dff","#636bb7","#bfc5ff","#014443","#195637",
    "#117744","#60ffaf","#b7ffdb","#825121","#ea7f17",
    "#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d",
    "#5b5b19","#fcfc00","#ffff9e","#ffb7ef","#fa7efc",
    "#ae09ea","#521899","#1e0047")


tree_AMF <- 
  ggtree(tree_raxml_AMF, layout="rectangular", branch.length="branch.length") %<+% tree_metadata +
  geom_tippoint(aes(color = BestMatch, size = Abundance), alpha = 0.8)  +
  scale_color_manual(values = palette_bestmatch) +
  guides(color = guide_legend(ncol = 1, override.aes = list(shape = 15, size=5))) 


tree_AMF + 
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = OTU, fill = Fertilized), 
    pwidth = 0.1) +
  scale_fill_gradient2() + 
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = OTU, fill = Control), 
    pwidth = 0.1) +
  scale_fill_gradient2()
  

abund_4_tree <- 
  melted_AMF %>%
  group_by(OTU, site, Family, Genus, BestMatch) %>% # no fert_status
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  rename(label = OTU, abundance = Abundance) %>%
  mutate(Family = as.character(Family))

head(abund_4_tree)
dim(abund_4_tree)
any(is.na(abund_4_tree))
levels(as.factor(abund_4_tree$Family))

# Tree to tibble 
tree_data <- as_tibble(tree_raxml_AMF)
tree_data

tree_metadata <- 
  abund_4_tree %>%
  tidyr::pivot_wider(
    id_cols = c(label, Family, Genus, BestMatch),
    names_from = "site",
    values_from = abundance) %>% 
  mutate(Genus = ifelse(is.na(Genus), "Unclassified", Genus))

dim(tree_metadata)
tree_metadata

# create list of OTUs per Genus
genus_list <- 
  tree_metadata %>%
  mutate(Genus = ifelse(is.na(Genus), "Unclassified", Genus)) %>%
  group_by(Genus) %>%
  summarise(tips = list(label), .groups = "drop") %>%
  { setNames(.$tips, .$Genus) }

tree_grouped <- groupOTU(tree_raxml_AMF, genus_list)

genus_list <- setNames(genus_list$tips, genus_list$Genus)

# group tree
tree_grouped <- groupOTU(tree_raxml_AMF, genus_list)
tree_grouped

tree_grouped@data$group[tree_grouped@data$group == "0"] <- NA

# 2) Build plot version 2 -
ggtree(tree_raxml_AMF) %<+% tree_metadata +
  geom_tippoint(aes(color = Genus, size = abundance), alpha = 0.8)  +
  #geom_tiplab(aes(color = Genus), size = 2) +
  theme(legend.position="right") +
  scale_color_manual(values = palette_taxa,na.translate = FALSE) +
  guides(color = guide_legend(ncol = 1, override.aes = list(shape = 15, size=3))) 


tree_AMF <- 
  ggtree(tree_grouped, aes(color = group), size = 0.3) +
  scale_color_discrete(na.translate = FALSE)

tree_AMF

tree_AMF <- ggtree(as.treedata(tree_tbl2), layout="rectangular", branch.length="none")

tree_AMF <- tree_AMF + geom_tree(aes(color = Genus)) +
  theme(legend.position="right") +
  scale_color_manual(values = palette_taxa,na.translate = FALSE) +
  guides(color = guide_legend(ncol = 1, override.aes = list(shape = 15, size=3))) 

tree_AMF



genus_list <- 
  tree_metadata %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus) %>%
  summarise(tips = list(label), .groups = "drop")


tree_tbl2 <- tree_data %>%
  left_join(tree_metadata %>% distinct(label, .keep_all = TRUE),
            by = "label")

as.treedata(tree_tbl2)

tree_grouped@data$group[tree_grouped@data$group == "0"] <- NA

ggtree(tree_grouped, aes(color = group), size = 0.3, layout="rectangular") +
  scale_color_discrete(na.translate = FALSE)

setdiff(tree_grouped$tip.label, tree_metadata$label)
setdiff(tree_metadata$label, tree_grouped$tip.label)


ggtree(tree_grouped, aes(color = group), size = 0.3) +
  scale_color_discrete(breaks = function(x) x[x != "0"])

tree_AMF <- 
  ggtree(tree_raxml_AMF, layout = "rectangular", branch.length = "none") %<+% abund_4_tree +
  geom_tippoint(aes(color = Genus), size = 1) +
  theme(legend.position = "right")

tree_AMF <- ggtree(tree_raxml_AMF, layout = "rectangular", branch.length = "none") %<+% abund_4_tree +
  geom_tiplab(aes(color = Family), size = 2) +
  theme(legend.position = "right")

tree_AMF



# Add the Prevalence Heatmap Ring
tree_AMF <- tree_AMF + geom_fruit(
  geom = geom_tile,
  mapping = aes(y=label, x=site, fill=abundance),
  pwidth = 0.1,
  axis.params = list(axis="x", text.size=2),
  grid.params = list()
) +
  scale_fill_viridis_c()

tree_AMF
taxonomy_4_tree
# Add the Abundance Bar Ring
tree_AMF <- tree_AMF + geom_fruit(
  geom = geom_bar,
  mapping = aes(y=label, x=abundance, fill=site),
  stat="identity",
  orientation="y",
  offset=0.15,
  pwidth = 0.1
)

tree_AMF

tree_AMF <- tree_AMF + 
  geom_star(
    aes(fill=Family, starshape=fert_status, size=Abundance), 
    starstroke=0.3,
    na.rm = TRUE
  ) + 
  scale_starshape_manual(values=c(15, 1)) # 15 is a square, 1 is a circle

tree_AMF


taxonomy_4_tree <- 
  as.data.frame(tax_table(physeq_AMF_rare)) %>% 
  dplyr::select(-Query, -Kingdom, -S_score) %>% 
  rownames_to_column("tip.label")

head(taxonomy_4_tree)
dim(taxonomy_4_tree)

all(tree_raxml_AMF$tip.label %in% taxonomy_4_tree$tip.label)
all(taxonomy_4_tree$tip.label %in% tree_raxml_AMF$tip.label)

setdiff(tree_raxml_AMF$tip.label, taxonomy_4_tree$tip.label)
setdiff(taxonomy_4_tree$tip.label, tree_raxml_AMF$tip.label)



metadata_4_tree <- 
  metadata_99 %>% rownames_to_column("label")

str(metadata_4_tree)
head(metadata_4_tree)

taxonomy_4_tree <- as.data.frame(tax_table(physeq_AMF_rare)) %>% 
  dplyr::select(-Query, -Kingdom, -S_score) %>% 
  rownames_to_column("label")

head(taxonomy_4_tree)
str(taxonomy_4_tree)



# Start over 
# https://yulab-smu.top/treedata-book/chapter10.html

library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(ggnewscale)
library(dplyr)
library(phyloseq)

tree_raxml_AMF <- phy_tree(physeq_AMF_rare)

dat1 <- melted_AMF %>%
  group_by(OTU, Genus, Phylum, Family) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus)) %>% 
  rename(label = OTU)

dat2 <- melted_AMF %>%
  group_by(OTU, site) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  rename(ID = OTU, Sites = site)


dat3 <- melted_AMF %>%
  group_by(OTU, site) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  rename(ID = OTU, Sites = site)

# ── 5. Color palette for genera ───────────────────────────────────────────────
genus_colors <- c(
  "Acaulospora"    = "#E41A1C", "Ambispora"      = "#4575B4",
  "Archaeospora"   = "#762A83", "Cetraspora"     = "#F4A582",
  "Complexispora"  = "#1B7837", "Dentiscutata"   = "#FFFFBF",
  "Diversispora"   = "#2166AC", "Dominikia"      = "#000080",
  "Entrophospora"  = "#808080", "Funneliformis"  = "#00441B",
  "Gigaspora"      = "#FDB863", "Glomus"         = "#9ECAE1",
  "Microdominikia" = "#1A1A2E", "Microkamienskia"= "#BDBDBD",
  "Oehlia"         = "#66C2A5", "Paraglomus"     = "#F46D43",
  "Rhizophagus"    = "#FF69B4", "Scutellospora"  = "#4D0000",
  "Septoglomus"    = "#74C476", "Silvaspora"     = "#8B8B00",
  "Unclassified"   = "#000000"
)

site_colors <- c(
  "Lux Arbor"  = "#0000FF",
  "Escanaba"   = "#FFA500",
  "Maple River"= "#FF0000"   # add/rename to match your actual site names
)

# Build the plot
tree_AMF <- ggtree(tree_raxml_AMF, layout = "fan", size = 0.15, open.angle = 5)
tree_AMF

# Layer 1: tip points — color = Genus, size = mean abundance
tree_AMF <- tree_AMF %<+% dat1 +
  geom_tippoint(
    aes(fill = Genus, size = log1p(MeanAbundance)),
    shape = 21, stroke = 0.1, alpha = 0.85) +
  scale_fill_manual(
    values = palette_bestmatch,
    guide  = guide_legend(
      keywidth = 0.5, keyheight = 0.5, order = 1,
      override.aes = list(shape = 21, size = 3)
    ),
    na.translate = FALSE
  ) +
  scale_size_continuous(
    range = c(0.1, 4),
    name  = "log(Abundance + 1)",
    guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 2,
                         override.aes = list(shape = 21))
  )

# Layer 2: site abundance heatmap ring
tree_AMF <- 
  tree_AMF +
  new_scale_fill() +
  geom_fruit(
    data    = dat2,
    geom    = geom_tile,
    mapping = aes(y = ID, x = Sites, alpha = Abundance, fill = Sites),
    color   = "grey90", offset = 0.04, size = 0.02
  ) +
  scale_alpha_continuous(
    range = c(0, 1),
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 4)
  ) +
  scale_fill_manual(
    values = palette_site,
    guide  = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 3)
  )

# Layer 3: total abundance bar ring
tree_AMF <- 
  tree_AMF +
  new_scale_fill() +
  geom_fruit(
    data        = dat3,
    geom        = geom_bar,
    mapping     = aes(y = ID, x = TotalAbundance, fill = Sites),
    pwidth      = 0.38,
    orientation = "y",
    stat        = "identity"
  ) +
  scale_fill_manual(
    values = palette_site,
    guide  = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 3)
  )

# ── 7. Theme ──────────────────────────────────────────────────────────────────
tree_AMF <- 
  tree_AMF +
  geom_treescale(fontsize = 2, linesize = 0.3) +
  theme(
    legend.position   = c(0.93, 0.5),
    legend.background = element_rect(fill = NA),
    legend.title      = element_text(size = 6.5),
    legend.text       = element_text(size = 4.5),
    legend.spacing.y  = unit(0.02, "cm")
  )

tree_AMF

tree_AMF + layout_rectangular() + 
  theme(legend.position=c(.05, .7))





tree <- tree_hmptree
dat1 <- df_tippoint
dat2 <- df_ring_heatmap
dat3 <- df_barplot_attr

# adjust the order
dat2$Sites <- factor(dat2$Sites, 
                     levels=c("Stool (prevalence)", "Cheek (prevalence)",
                              "Plaque (prevalence)","Tongue (prevalence)",
                              "Nose (prevalence)", "Vagina (prevalence)",
                              "Skin (prevalence)"))
dat3$Sites <- factor(dat3$Sites, 
                     levels=c("Stool (prevalence)", "Cheek (prevalence)",
                              "Plaque (prevalence)", "Tongue (prevalence)",
                              "Nose (prevalence)", "Vagina (prevalence)",
                              "Skin (prevalence)"))
# extract the clade label information. Because some nodes of tree are
# annotated to genera, which can be displayed with high light using ggtree.
nodeids <- nodeid(tree, tree$node.label[nchar(tree$node.label)>4])
nodedf <- data.frame(node=nodeids)
nodelab <- gsub("[\\.0-9]", "", tree$node.label[nchar(tree$node.label)>4])
# The layers of clade and hightlight
poslist <- c(1.6, 1.4, 1.6, 0.8, 0.1, 0.25, 1.6, 1.6, 1.2, 0.4,
             1.2, 1.8, 0.3, 0.8, 0.4, 0.3, 0.4, 0.4, 0.4, 0.6,
             0.3, 0.4, 0.3)
labdf <- data.frame(node=nodeids, label=nodelab, pos=poslist)

# The circular layout tree.
p <- ggtree(tree, layout="fan", size=0.15, open.angle=5) +
  geom_hilight(data=nodedf, mapping=aes(node=node),
               extendto=6.8, alpha=0.3, fill="grey", color="grey50",
               size=0.05) +
  geom_cladelab(data=labdf, 
                mapping=aes(node=node, 
                            label=label,
                            offset.text=pos),
                hjust=0.5,
                angle="auto",
                barsize=NA,
                horizontal=FALSE, 
                fontsize=1.4,
                fontface="italic"
  )

p <- p %<+% dat1 + geom_star(
  mapping=aes(fill=Phylum, starshape=Type, size=Size),
  position="identity",starstroke=0.1) +
  scale_fill_manual(values=c("#FFC125","#87CEFA","#7B68EE","#808080",
                             "#800080", "#9ACD32","#D15FEE","#FFC0CB",
                             "#EE6A50","#8DEEEE", "#006400","#800000",
                             "#B0171F","#191970"),
                    guide=guide_legend(keywidth = 0.5, 
                                       keyheight = 0.5, order=1,
                                       override.aes=list(starshape=15)),
                    na.translate=FALSE)+
  scale_starshape_manual(values=c(15, 1),
                         guide=guide_legend(keywidth = 0.5, 
                                            keyheight = 0.5, order=2),
                         na.translate=FALSE)+
  scale_size_continuous(range = c(1, 2.5),
                        guide = guide_legend(keywidth = 0.5, 
                                             keyheight = 0.5, order=3,
                                             override.aes=list(starshape=15)))

p <- p + new_scale_fill() +
  geom_fruit(data=dat2, geom=geom_tile,
             mapping=aes(y=ID, x=Sites, alpha=Abundance, fill=Sites),
             color = "grey50", offset = 0.04,size = 0.02)+
  scale_alpha_continuous(range=c(0, 1),
                         guide=guide_legend(keywidth = 0.3, 
                                            keyheight = 0.3, order=5)) +
  geom_fruit(data=dat3, geom=geom_bar,
             mapping=aes(y=ID, x=HigherAbundance, fill=Sites),
             pwidth=0.38, 
             orientation="y", 
             stat="identity",
  ) +
  scale_fill_manual(values=c("#0000FF","#FFA500","#FF0000",
                             "#800000", "#006400","#800080","#696969"),
                    guide=guide_legend(keywidth = 0.3, 
                                       keyheight = 0.3, order=4))+
  geom_treescale(fontsize=2, linesize=0.3, x=4.9, y=0.1) +
  theme(legend.position=c(0.93, 0.5),
        legend.background=element_rect(fill=NA),
        legend.title=element_text(size=6.5),
        legend.text=element_text(size=4.5),
        legend.spacing.y = unit(0.02, "cm"),
  )
p





# ***********************************************************************-------
# **** 6. CHEMISTRY DATA **** --------------------------------------------------

# Import the chemistry and yield data ------------------------------------------
chem_data <- 
  read.csv(file = file.path(data_path, "datasets/soil_metadata_yield.csv")) %>%
  janitor::clean_names()

head(chem_data)

# Average the AMF counts by site, fert_status, and plot_rep --------------------
sapply(sample_data(physeq_AMF_rare), class)

mean_physeq_AMF_rare <-
  physeq_AMF_rare %>% 
  speedyseq::select_sample_data(site, fert_status, plot_rep) %>% 
  speedyseq::mutate_sample_data(
    site_name = gsub(" ", "", site),
    sample_id = paste(site_name, fert_status, plot_rep , sep="_")
  ) %>% 
  merge_samples("sample_id", fun = mean) 

mean_physeq_AMF_rare@sam_data
mean_physeq_AMF_rare@otu_table

mean_otutable_AMF_rare <- 
  as.data.frame(otu_table(mean_physeq_AMF_rare)) 
head(mean_otutable_AMF_rare)

#mean_otutable_AMF_rare$hill_q0 <- specnumber(mean_otutable_AMF_rare)
#mean_otutable_AMF_rare$hill_q1 <- exp(diversity(mean_otutable_AMF_rare, index = "shannon"))
#mean_otutable_AMF_rare$hill_q2 <- diversity(mean_otutable_AMF_rare, index = "invsimpson")
mean_otutable_AMF_rare$hill_q0 <- vegan::renyi(x = mean_otutable_AMF_rare, scales = c(0), hill = TRUE)
mean_otutable_AMF_rare$hill_q1 <- vegan::renyi(x = mean_otutable_AMF_rare, scales = c(1), hill = TRUE)
mean_otutable_AMF_rare$hill_q2 <- vegan::renyi(x= mean_otutable_AMF_rare, scales = c(2), hill = TRUE)

mean_otutable_AMF_rare <- 
  mean_otutable_AMF_rare %>% 
  mutate(sample_id = rownames(.)) %>% 
  separate(sample_id, into = c("site", "fert_status", "plot_rep"), 
           sep = "_", remove = FALSE)

head(mean_otutable_AMF_rare)
setdiff(mean_otutable_AMF_rare$sample_id, chem_data$sample_id)

# OPTIONAL. We already have too many predictors and not enoug observations. 
# Adding community composition metrics, bray and jaccard.

pcoa_meanAMF_bray <- 
  ordinate(mean_physeq_AMF_rare, method="PCoA", distance="bray")

pcoa_meanAMF_jacc <- 
  ordinate(mean_physeq_AMF_rare, method="PCoA", distance="jaccard")

# combined dataset
otu_chem_data <-
  mean_otutable_AMF_rare %>% 
  left_join(chem_data, by = "sample_id") %>% 
  left_join(
    as.data.frame(pcoa_meanAMF_bray$vectors) %>% 
      dplyr::select(Axis.1, Axis.2) %>% 
      rownames_to_column("sample_id"), 
    by= "sample_id") %>% 
  rename(bray.1 = Axis.1, bray.2 = Axis.2) %>% 
  left_join(
    as.data.frame(pcoa_meanAMF_jacc$vectors) %>% 
      dplyr::select(Axis.1, Axis.2) %>% 
      rownames_to_column("sample_id"),
    by= "sample_id") %>% 
  rename(jacc.1 = Axis.1, jacc.2 = Axis.2)

head(otu_chem_data)

ggplot(otu_chem_data, aes(x=bray.1, y=bray.2, 
                          col = site.x, shape = fert_status)) +
  geom_point()

ggplot(otu_chem_data, aes(x=jacc.1, y=jacc.2, 
                          col = site.x, shape = fert_status)) +
  geom_point()

adonis2(
  otu_chem_data ~ fert_status,
  data = meta_rare,
  method = "bray",
  permutations = how(blocks = meta_rare$site, nperm = 999),
  by = "margin"
)


# Testing effect of N and P on biomass and hill 0 ------------------------------
lmer_otu_chem_data <-
  otu_chem_data %>% 
  dplyr::select(sample_id, hill_q0, hill_q2, dry_matter_yield_mg_ha, p_ppm, site.x, fert_status) %>% 
  column_to_rownames("sample_id") %>% 
  rename(site = site.x, dry_matter = dry_matter_yield_mg_ha) %>% 
  mutate(across(.cols = c(1:4), .fns = as.numeric)) %>% 
  mutate(site = as.factor(site), 
         fert_status = as.factor(fert_status))

lmer_otu_chem_data
str(lmer_otu_chem_data)


# MIXED MODELS -----------------------------------------------------------------

# Distribution of outcome variable
ggarrange(
  lmer_otu_chem_data %>% 
    ggplot2::ggplot(aes(x = hill_q0)) +
    geom_histogram(binwidth = 20) +
    labs(title = "hill_q0") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
    geom_vline(aes(
      xintercept = mean(hill_q0)), 
      color = "darkred",
      linetype = "dashed", linewidth = 2, show.legend = FALSE),
  lmer_otu_chem_data %>% 
    ggplot2::ggplot(aes(x = hill_q2)) +
    geom_histogram(binwidth = 5) +
    labs(title = "hill_q2") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
    geom_vline(aes(
      xintercept = mean(hill_q2)), 
      color = "darkred",
      linetype = "dashed", linewidth = 2, show.legend = FALSE),
  lmer_otu_chem_data %>% 
    ggplot2::ggplot(aes(x = dry_matter)) +
    geom_histogram(binwidth = 2) +
    labs(title = "Yield") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
    geom_vline(aes(
      xintercept = mean(dry_matter)), 
      color = "darkred",
      linetype = "dashed", linewidth = 2, show.legend = FALSE),
  ncol = 3,
  nrow = 1)

# NOTE. Running the model with fixed slope as the standard in this case. 

# Model for hill_q0
fit_hill_q0_fixslope <- lmer(
  hill_q0 ~ fert_status + p_ppm + fert_status:p_ppm + (1 | site),
  data = lmer_otu_chem_data, REML = FALSE
)

summary(fit_hill_q0_fixslope)
anova(fit_hill_q0_fixslope)
plot(fit_hill_q0_fixslope)

diagnostics_dharma(
  model     = fit_hill_q0_fixslope,
  group_var1 = lmer_otu_chem_data$site,
  group_var2 = NULL
)

# If I try running glmmTMB still not significant
glmmTMB(hill_q0 ~ fert_status + p_ppm + fert_status:p_ppm + (1 | site),
        dispformula = ~site, data = lmer_otu_chem_data, 
        REML = FALSE) %>% summary()

# Model for diversity
fit_hill_q2_fixslope <- lmer(
  hill_q2 ~ fert_status + p_ppm + fert_status:p_ppm + (1 | site),
  data = lmer_otu_chem_data, REML = FALSE
)

summary(fit_hill_q2_fixslope)
plot(fit_hill_q2_fixslope)

diagnostics_dharma(
  model     = fit_hill_q2_fixslope,
  group_var1 = lmer_otu_chem_data$site,
  group_var2 = NULL
)

# Model for yield
fit_yield_fixslope <- lmer(
  dry_matter ~ fert_status + p_ppm + fert_status:p_ppm + (1 | site),
  data = lmer_otu_chem_data, REML = FALSE
)

summary(fit_yield_fixslope)
anova(fit_yield_fixslope)

diagnostics_dharma(
  model     = fit_yield_fixslope,
  group_var1 = lmer_otu_chem_data$site,
  group_var2 = NULL
)

# INTERPRETATION. The DHARMa diagnostic is checking whether your model residuals 
# are uniformly distributed between 0 and 1 (i.e., whether the model fits well), 
# not whether your raw data are normal. 
# dry_mass is normally distributed meaning the marginal distribution looks okay,
# However, the categorical test is flagging that after accounting for your predictors, 
# the residuals within Escanaba and Rhinelander deviate significantly from uniformity
# (shown in red). 

# The real problem: your model is fitting the overall data well on average, but is
# systematically mis-fitting Escanaba and Rhinelander specifically. 
# This suggests:
# _ Those sites may have different variance (heteroscedasticity by site)!
# _ There could be a site-level effect not captured in your fixed structure, aka 
#    dry_matter ~ fert_status + p_ppm + fert_status:p_ppm 
#    And the interaction term (e.g., treatment × site) which is already present. 
# _ Consider adding site as a random effect if not already included, or allowing 
#   random slopes.

# If Escanaba and Rhinelander are no longer red, the random effect successfully 
# absorbed the site-level variance. I may need random slopes, as it allows 
# fert_status effect to vary by site

fit_yield_randslope <-
  lmer(dry_matter ~ fert_status + p_ppm + fert_status:p_ppm + (1 + fert_status | site),
       data = lmer_otu_chem_data, REML = FALSE)

summary(fit_yield_randslope)
plot(fit_yield_randslope)

diagnostics_dharma(
  model      = fit_yield_randslope,
  group_var1 = lmer_otu_chem_data$site,
  group_var2 = NULL
)

anova(fit_yield_fixslope, fit_yield_randslope)

# INTERPRETATION. 
# _ The singular fit warning is the key issue here, and it's telling me the random
#   slopes model is too complex for your data. Look at this in the random
#   effects: fert_statusFertilized  0.887 and 0.942 and 1.00  ← correlation = 1.00
#   (perfect!)
# _ The categorical DHARMa test is still red after that, it may simply reflect 
#   genuine heterogeneity across sites that you cannot model away with only 5 
#   levels. In that case you have a few options:
#   _ Accept it and note it as a limitation — the overall fit is excellent (KS p = 0.996,
#     dispersion fine, no outliers), and the site-level issue may be unavoidable.
#   _ Use a site-stratified analysis as a sensitivity check.
#   _ Switch to glmmTMB (Fit a generalized linear mixed model (GLMM) using Template
#     Model Builder (TMB).) which allows more flexible variance structures.

# glm model for dry_matter

fitGLM_yield_fixslope <- 
  glmmTMB(dry_matter ~ fert_status * p_ppm + (1|site),
                    dispformula = ~site, data = lmer_otu_chem_data, 
          REML = FALSE)

fitGLM_yield_fixslope_no_p  <- 
  glmmTMB(dry_matter ~ fert_status + p_ppm + (1|site),
                    dispformula = ~site, data = lmer_otu_chem_data, 
          REML = FALSE)

anova(fitGLM_yield_fixslope_no_p, fitGLM_yield_fixslope)

# Adding the interaction doe snot imporve the model fit!

# NOTE. Use ML (REML = FALSE) for comparing fixed-effect structures. If LRT confirms
# no interaction, refit the additive model with REML for reporting.

# INTERPRETATION. Think of “family” as the distribution of Y and of “link” as 
# the transformation of the mean of Y that makes the relationship with 
# predictors linear.

# _ Family: which probability distribution you assume for the response (Gaussian, 
#   Binomial, Poisson, Gamma, etc.). It controls the mean–variance relationship 
#   and the likelihood.
# _ Link: the function that maps the mean onto the real line so a linear predictor
#   makes sense, and it determines how you interpret coefficients (differences, 
#   ratios, odds ratios).

# What about p_ppm? It is not significant, but it is not causing the model to fail.
# Stepwise removal based on p-values inflates Type I error on the remaining terms, 
# biases their estimates away from zero, and produces overconfident standard errors. 

# Is p_ppm there for a scientific reason? If your hypothesis or experimental design
# involves soil P as a relevant covariate — which seems likely given you're 
# studying fertilization effects — then it belongs in the model regardless of 
# p-value. You're controlling for baseline soil P differences across plots and 
# sites, which is exactly what a covariate should do.

# Does it confound the fertilization effect? Compare the fert_status coefficient
# with and without p_ppm:

fitGLM_yield_fixslope_with_p <- glmmTMB(dry_matter ~ fert_status + p_ppm + (1|site),
                      dispformula = ~site, data = lmer_otu_chem_data)
fitGLM_yield_fixslope_no_p   <- glmmTMB(dry_matter ~ fert_status + (1|site),
                      dispformula = ~site, data = lmer_otu_chem_data)

summary(fitGLM_yield_fixslope_with_p)$coefficients$cond
summary(fitGLM_yield_fixslope_no_p)$coefficients$cond

# If the fertilization estimate changes meaningfully when you remove p_ppm, 
# then p_ppm is doing real adjustment work even though its own coefficient 
# isn't significant. Keep it.

# The fertilization effect changes by 0.028 — roughly 1% of its magnitude, and 
# well within a fraction of one SE. p_ppm is doing essentially no confounding 
# adjustment work. The fertilization effect stands on its own.

AIC(fitGLM_yield_fixslope_with_p, fitGLM_yield_fixslope_no_p)
183.802-182.145

# ΔAIC = 1.66 confirms the two additive models are essentially equivalent in fit,
# so parsimony favors the simpler one.

# NOTE. Final model. I will refit with no interaction and no p_ppm!
fitGLM_yield_fixslope <- 
  glmmTMB(dry_matter ~ fert_status + (1 | site),
          dispformula = ~site,   # allows variance to differ by site
          data = lmer_otu_chem_data, 
          REML = TRUE)

summary(fitGLM_yield_fixslope)


# FIGURE S3 glmmTMB diagnostics ------------------------------------------------
set.seed(3241)
plot_dharma_glm_dia <- 
  diagnostics_dharma(
  model      = fitGLM_yield_fixslope,
  group_var1 = lmer_otu_chem_data$site,
  group_var2 = NULL
)

plot_dharma_glm_dia

ggsave(
  file.path(data_path, "results/Fig_XX_dharma_glm_dia.pdf"),
    plot_dharma_glm_dia,
  device = "pdf"
)


# Visualize --------------------------------------------------------------------

fitGLM_yield_fixslope_pred_fixed <- 
  ggpredict(fitGLM_yield_fixslope, terms = "fert_status")

fitGLM_yield_fixslope_pred_site  <- 
  ggpredict(fitGLM_yield_fixslope, terms = c("fert_status", "site"), type = "random")

# glmmTMB model ----------------------------------------------------------------
glm_model_plot <- 
  ggplot() +
  # Raw data, jittered
  geom_jitter(data = lmer_otu_chem_data,
              aes(x = fert_status, y = dry_matter, color = site),
              width = 0.15, alpha = 0.8, size = 2, shape = 16) +
  # Site-level means connected by lines (shows random intercepts)
  geom_line(data = fitGLM_yield_fixslope_pred_site,
            aes(x = x, y = predicted, group = group, color = group),
            alpha = 0.5, linewidth = 0.6, show.legend = FALSE) +
  #geom_point(data = fitGLM_yield_fixslope_pred_site,
  #           aes(x = x, y = predicted, color = group),
  #           size = 2.5) +
  # Population-level mean with CI
  #geom_errorbar(data = fitGLM_yield_fixslope_pred_fixed,
  #              aes(x = x, ymin = conf.low, ymax = conf.high),
  #              width = 0.1, linewidth = 0.8) +
  geom_point(data = fitGLM_yield_fixslope_pred_fixed,
             aes(x = x, y = predicted),
             size = 4, shape = 18) +
  geom_line(data = fitGLM_yield_fixslope_pred_fixed,
            aes(x = x, y = predicted, group = 1),
            linewidth = 1.5) +
  labs(title = "Switchgrass yield ",
       x = "Treatment",
       y = "Dry matter yield (mg/ha)") +
  scale_color_manual(values = palette_site)+
  theme_classic()+ 
  theme(
    plot.title = element_markdown(size =12, face = "bold",hjust = 0.5, vjust = 0.5),
    plot.subtitle = element_markdown(size = 10, hjust = 0.5, vjust = 0.5),
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 8),
    legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"),
    legend.title = element_blank(), legend.text = element_text(size = 8)
    #legend.margin=ggplot2::margin(0,5,0,0),
    #legend.box.margin=ggplot2::margin(0,5,0,0)
  ) +
  scale_y_continuous(n.breaks = 10)

glm_model_plot

# Raw data and dispersion ------------------------------------------------------

x_limits_yield <- 
  lmer_otu_chem_data %>%
  group_by(site) %>% 
  group_map(~ broom.mixed::tidy(fitGLM_yield_fixslope, effects = "ran_vals") ) %>%
  bind_rows() %>%
  summarise(
    min = min(estimate - 1.96*std.error),
    max = max(estimate + 1.96*std.error)
  )

x_limits_yield

intercept_yield <- fixef(fitGLM_yield_fixslope)[[1]]["fert_statusFertilized"]
intercept_yield

plot_raw_chem_data <- 
  ggplot(lmer_otu_chem_data, 
       aes(x = fert_status, y =dry_matter  , color = site)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5),
             alpha = 0.8, size = 2, shape = 16) +
  stat_summary(aes(group = site), fun = mean, geom = "line", 
               position = position_dodge(0.5)) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", 
               color = "black", size = 4, shape = 18) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", 
               color = "black", size = 1) +
  labs(title = "Switchgrass yield ",
       #subtitle = "Black line connects fixed effect means",
       x = "Treatment",
       y = "Dry matter yield (mg/ha)") +
  scale_color_manual(values = palette_site) +
  theme_classic() +
  theme(
    plot.title = element_markdown(size =12, face = "bold",hjust = 0.5, vjust = 0.5),
    plot.subtitle = element_markdown(size = 10,hjust = 0.5, vjust = 0.5),
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 8),
    legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"),
    legend.title = element_blank(), legend.text = element_text(size = 8)
    #legend.margin=ggplot2::margin(0,5,0,0),
    #legend.box.margin=ggplot2::margin(0,5,0,0)
  )+
  guides(color = guide_legend(ncol=1))

plot_raw_chem_data

lmer_otu_chem_data %>% dplyr::select(dry_matter, site) %>% 
  mutate(biomass = dry_matter) %>% 
  ggplot(aes(x = site, y = biomass)) + 
  geom_boxplot() +
  coord_flip()

lmer_otu_chem_data %>% 
  group_by(site, fert_status) %>% 
  summarise(sd_within = sd(dry_matter), n = n(), .groups = "drop")


plot_deviation_from_fix <- 
  broom.mixed::tidy(fitGLM_yield_fixslope, effects = "ran_vals") %>% 
  as.data.frame() %>%
  dplyr::rename( site = "level") %>%  
  arrange(site) %>% 
  ggplot(aes(x = estimate, y = site, color = site)) +
  geom_point(size = 4, shape = 18) +
  geom_errorbar(aes(xmin = estimate - 1.96*std.error, 
                    xmax = estimate + 1.96*std.error),
                height = 0, 
                orientation = "y") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Deviation from Fixed Effect",
       x = "Deviation from fixed-effect means",
       y = NULL) +
  scale_color_manual(values = palette_site)+
  theme_classic()+ 
  theme(
    plot.title = element_markdown(size =12, face = "bold",hjust = 0.5, vjust = 0.5),
    plot.subtitle = element_markdown(size = 10, hjust = 0.5, vjust = 0.5),
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 8),
    legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"),
    legend.title = element_blank(), legend.text = element_text(size = 8)
    #legend.margin=ggplot2::margin(0,5,0,0),
    #legend.box.margin=ggplot2::margin(0,5,0,0)
  ) +
  guides(color = guide_legend(ncol = 1)) +
  scale_x_continuous(
    limits = c(x_limits_yield$min, x_limits_yield$max),
    sec.axis = sec_axis(
      transform = ~ . + intercept_yield,
      name = "Dry matter yield (mg/ha)")
  )

plot_deviation_from_fix

# If you top show disperions then mayeb consider the code below.

# To get a truly linear mg/ha axis, transform the data itself before plotting 
# and let the asymmetric CIs show as asymmetric error bars. Symmetric Wald CIs 
# live naturally on the log scale, and showing them on a log-spaced axis preserves
# that symmetry.

# Pull dispersion coefs
# and SEs straight from the model
disp_est <- fixef(fitGLM_yield_fixslope)$disp
disp_se  <- sqrt(diag(vcov(fitGLM_yield_fixslope)$disp))
log_intercept <- disp_est[1]

plot_deviation_from_fix <- 
  data.frame(
    term      = names(disp_est),
    estimate  = as.numeric(disp_est),
    std.error = as.numeric(disp_se)
  ) %>% 
  mutate(
    site = ifelse(term == "(Intercept)", 
                  "Escanaba", 
                  gsub("^site", "", term)),
    log_sd    = ifelse(term == "(Intercept)", 
                       estimate, 
                       estimate + log_intercept),
    # Back-transform point estimate and CI bounds to mg/ha
    sd_mg_ha  = exp(log_sd),
    sd_lo     = exp(log_sd - 1.96 * std.error),
    sd_hi     = exp(log_sd + 1.96 * std.error)
  ) %>% 
  arrange(site) %>% 
  ggplot(aes(x = sd_mg_ha, y = site, color = site)) +
  geom_point(size = 4, shape = 18) +
  geom_errorbar(aes(xmin = sd_lo, xmax = sd_hi),
                height = 0, 
                orientation = "y") +
  geom_vline(xintercept = exp(log_intercept), 
             linetype = "dashed", color = "black") +
  labs(title = "Within-site residual variability",
       x = "Residual SD of dry matter yield (mg/ha)",
       y = NULL) +
  scale_color_manual(values = palette_site) +
  theme_classic() +
  theme(
    plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
    axis.text  = element_markdown(size = 8),
    legend.key.height = unit(0.5, "cm"), 
    legend.key.width  = unit(0.5, "cm"),
    legend.title = element_blank(), 
    legend.text  = element_text(size = 8)
  ) +
  guides(color = guide_legend(ncol = 1))

plot_deviation_from_fix

# **** FIGURE 4 - glmmTMB ------------------------------------------------------

Figure_X_yield <-
  ggarrange(
    plot_raw_chem_data,
    glm_model_plot,
    plot_deviation_from_fix,
    ncol = 3, nrow = 1,
    labels = c("A","B", "C"),
    #align = "hv",
    legend = "right", 
    common.legend = TRUE)

Figure_X_yield

ggsave(
  file.path(data_path, "results/Fig_XX_yield.pdf"),
  plot = ggpubr::annotate_figure(
    Figure_X_yield,
    top = text_grob("EFFECT OF NITROGEN ON SWITCHGRASS YEILD", size = 12, face = "bold")
  ),
  device = "pdf"
)

# Envfit analysis --------------------------------------------------------------

otumat_mean_AMF_rare <- as(otu_table(mean_physeq_AMF_rare), "matrix")
head(otumat_mean_AMF_rare)

metadata_mean_AMF_rare <- data.frame(sample_data(mean_physeq_AMF_rare))
metadata_mean_AMF_rare

# Extract dataset
env_fit_chem_data <- 
  chem_data %>% 
  column_to_rownames("sample_id") %>% 
  mutate(
    treatment = as.numeric(factor(treatment)) - 1
  ) %>% 
  mutate(across(where(is.numeric) & !treatment, 
                ~ (. - min(.)) / (max(.) - min(.))))

env_fit_chem_data
dim(env_fit_chem_data)

# Reorder env_fit_chem_data to match the order of otumat_mean_AMF_rare
env_fit_chem_data <- 
  env_fit_chem_data[match(rownames(otumat_mean_AMF_rare), rownames(env_fit_chem_data)), ]

# Verify again
identical(rownames(otumat_mean_AMF_rare), rownames(env_fit_chem_data))

dist_mean_AMF <- vegdist(otumat_mean_AMF_rare, method = "bray")
dist_mean_AMF

pcoa_mean_AMF <- cmdscale(dist_mean_AMF, k = 2, eig = TRUE)
pcoa_mean_AMF

# Fit the data
env_fit_AMF <- envfit(pcoa_mean_AMF, env_fit_chem_data[, 2:11], 
                      permutations = how(blocks = env_fit_chem_data$site, nperm = 999),
                      na.rm = FALSE)

env_fit_AMF


env_fit_AMF_plot <-
plot_ordination_envft(ord = pcoa_mean_AMF, 
                      meta = env_fit_chem_data, 
                      col_var = "site", 
                      shape_var = "treatment", 
                      env = env_fit_AMF, 
                      arrow_mul = 0.2,
                      p_threshold = 0.05) +
  labs(title = "Principal Coordinate Analysis") +
  scale_color_manual(values = palette_site) +
  scale_shape_manual(values = c(16, 17),
                     labels = c("0" = "Control", "1" = "Fertilized")) +
  theme_classic() + 
  theme(
    plot.title = element_markdown(size =12, face = "bold",hjust = 0.5, vjust = 0.5),
    plot.subtitle = element_markdown(size = 10, hjust = 0.5, vjust = 0.5),
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 8),
    legend.key.height = unit(0.6, "cm"), 
    legend.key.width = unit(0.6, "cm"),
    legend.title = element_blank(), 
    legend.text = element_text(size = 8)
    #legend.margin=ggplot2::margin(0,5,0,0),
    #legend.box.margin=ggplot2::margin(0,5,0,0)
  ) +
  guides(
    color = guide_legend(ncol = 1, override.aes = list(shape = 15)))



env_fit_AMF_plot


# Create a text string of the significant results
sig_results <- env_fit_AMF$vectors$r[env_fit_AMF$vectors$pvals <= 0.05]
r2_text <- paste(paste0(names(sig_results), ": R²=", 
                        round(sig_results, 2)), collapse = "\n")

# Add to your ggplot
env_fit_AMF_plot + annotate("text", x = Inf, y = -Inf, label = r2_text, 
             hjust = 1.1, vjust = -0.5, size = 3, fontface = "italic")


Figure_X_yield <-
  ggarrange(
    glm_model_plot,
    plot_deviation_from_fix,
    env_fit_AMF_plot +
      annotate("text", x = Inf, y = -Inf, label = r2_text, 
               hjust = 1.1, vjust = -0.5, size = 3, fontface = "italic"),
    ncol = 3, nrow = 1,
    labels = c("A","B", "C"),
    widths = c(0.5, 0.6, 0.7),
    align = "h",
    legend = "bottom", 
    common.legend = TRUE)

Figure_X_yield

ggsave(
  file.path(data_path, "results/Fig_5_yield.pdf"),
  plot = ggpubr::annotate_figure(
    Figure_X_yield,
    top = text_grob("EFFECT OF NITROGEN ON SWITCHGRASS YEILD", size = 12, face = "bold")
  ),
  device = "pdf"
)


# **** END OF ANALYSIS **************************************************-------
# ***********************************************************************-------
# ***********************************************************************-------

# **********************************************************************--------
# BETA DIVERSITY 90% OTUS to ASVs ----------------------------------------------

# Import otu_tables ------------------------------------------------------------

otutable_90 <-read.delim(file.path(data_path, "datasets/otutab_90.txt"), row.names = 1)
otutable_91 <-read.delim(file.path(data_path, "datasets/otutab_91.txt"), row.names = 1)
otutable_92 <-read.delim(file.path(data_path, "datasets/otutab_92.txt"), row.names = 1)
otutable_93 <-read.delim(file.path(data_path, "datasets/otutab_93.txt"), row.names = 1)
otutable_94 <-read.delim(file.path(data_path, "datasets/otutab_94.txt"), row.names = 1)
otutable_95 <-read.delim(file.path(data_path, "datasets/otutab_95.txt"), row.names = 1)
otutable_96 <-read.delim(file.path(data_path, "datasets/otutab_96.txt"), row.names = 1)
otutable_97 <-read.delim(file.path(data_path, "datasets/otutab_97.txt"), row.names = 1)
otutable_98 <-read.delim(file.path(data_path, "datasets/otutab_98.txt"), row.names = 1)
otutable_99 <-read.delim(file.path(data_path, "datasets/otutab_99.txt"), row.names = 1)
otutable_100 <-read.delim(file.path(data_path, "datasets/otutab_100_asv.txt"), row.names = 1)

# Import taxonomy --------------------------------------------------------------
taxonomy_90 <- extract_blasTAX(
  tax_path = file.path(data_path, "datasets/taxonomy_blast_90.txt"),
  namemap_path = file.path(data_path, "datasets/name_mapping_90.txt")) %>% 
  extract_Glomero()

taxonomy_91 <- extract_blasTAX(
  tax_path = file.path(data_path, "datasets/taxonomy_blast_91.txt"),
  namemap_path = file.path(data_path, "datasets/name_mapping_91.txt")) %>% 
  extract_Glomero()

taxonomy_92 <- extract_blasTAX(
  tax_path = file.path(data_path, "datasets/taxonomy_blast_92.txt"),
  namemap_path = file.path(data_path, "datasets/name_mapping_92.txt")) %>% 
  extract_Glomero()

taxonomy_93 <- extract_blasTAX(
  tax_path = file.path(data_path, "datasets/taxonomy_blast_93.txt"),
  namemap_path = file.path(data_path, "datasets/name_mapping_93.txt")) %>% 
  extract_Glomero()

taxonomy_94 <- extract_blasTAX(
  tax_path = file.path(data_path, "datasets/taxonomy_blast_94.txt"),
  namemap_path = file.path(data_path, "datasets/name_mapping_94.txt")) %>% 
  extract_Glomero()

taxonomy_95 <- extract_blasTAX(
  tax_path = file.path(data_path, "datasets/taxonomy_blast_95.txt"),
  namemap_path = file.path(data_path, "datasets/name_mapping_95.txt")) %>% 
  extract_Glomero()

taxonomy_96 <- extract_blasTAX(
  tax_path = file.path(data_path, "datasets/taxonomy_blast_96.txt"),
  namemap_path = file.path(data_path, "datasets/name_mapping_96.txt")) %>% 
  extract_Glomero()

taxonomy_97 <- extract_blasTAX(
  tax_path = file.path(data_path, "datasets/taxonomy_blast_97.txt"),
  namemap_path = file.path(data_path, "datasets/name_mapping_97.txt")) %>% 
  extract_Glomero()

taxonomy_98 <- extract_blasTAX(
  tax_path = file.path(data_path, "datasets/taxonomy_blast_98.txt"),
  namemap_path = file.path(data_path, "datasets/name_mapping_98.txt")) %>% 
  extract_Glomero()

taxonomy_99 <- extract_blasTAX(
  tax_path = file.path(data_path, "datasets/taxonomy_blast_99.txt"),
  namemap_path = file.path(data_path, "datasets/name_mapping_99.txt")) %>% 
  extract_Glomero()

taxonomy_100 <- extract_blasTAX(
  tax_path = file.path(data_path, "datasets/taxonomy_blast_100_asv.txt"),
  namemap_path = file.path(data_path, "datasets/name_mapping_100_asv.txt")) %>% 
  extract_Glomero()

# Import sequences -------------------------------------------------------------
zotu_90 <- readDNAStringSet(file.path(data_path,"datasets/otus_90.fasta"), 
                            format="fasta", seek.first.rec=TRUE, use.names=TRUE)

zotu_91 <- readDNAStringSet(file.path(data_path,"datasets/otus_91.fasta"), 
                            format="fasta", seek.first.rec=TRUE, use.names=TRUE)

zotu_92 <- readDNAStringSet(file.path(data_path,"datasets/otus_92.fasta"), 
                            format="fasta", seek.first.rec=TRUE, use.names=TRUE)

zotu_93 <- readDNAStringSet(file.path(data_path,"datasets/otus_93.fasta"), 
                            format="fasta", seek.first.rec=TRUE, use.names=TRUE)

zotu_94 <- readDNAStringSet(file.path(data_path,"datasets/otus_94.fasta"), 
                            format="fasta", seek.first.rec=TRUE, use.names=TRUE)

zotu_95 <- readDNAStringSet(file.path(data_path,"datasets/otus_95.fasta"), 
                            format="fasta", seek.first.rec=TRUE, use.names=TRUE)

zotu_96 <- readDNAStringSet(file.path(data_path,"datasets/otus_96.fasta"), 
                            format="fasta", seek.first.rec=TRUE, use.names=TRUE)

zotu_97 <- readDNAStringSet(file.path(data_path,"datasets/otus_97.fasta"), 
                            format="fasta", seek.first.rec=TRUE, use.names=TRUE)

zotu_98 <- readDNAStringSet(file.path(data_path,"datasets/otus_98.fasta"), 
                            format="fasta", seek.first.rec=TRUE, use.names=TRUE)

zotu_99 <- readDNAStringSet(file.path(data_path,"datasets/otus_99.fasta"), 
                            format="fasta", seek.first.rec=TRUE, use.names=TRUE)

zotu_100 <- readDNAStringSet(file.path(data_path,"datasets/otu_100_asv.fasta"), 
                             format="fasta", seek.first.rec=TRUE, use.names=TRUE)

# Metadata ---------------------------------------------------------------------
metadata_99

# Create individual phyloseq objects -------------------------------------------
physeq_90 <- generate_phyloseq(otu=otutable_90,metadata=metadata_99,
                               taxonomy=taxonomy_90,sequences=zotu_90)
physeq_90
physeq_91 <- generate_phyloseq(otu=otutable_91,metadata=metadata_99,
                               taxonomy=taxonomy_91,sequences=zotu_91)
physeq_91
physeq_92 <- generate_phyloseq(otu=otutable_92,metadata=metadata_99,
                               taxonomy=taxonomy_92,sequences=zotu_92)
physeq_92
physeq_93 <- generate_phyloseq(otu=otutable_93,metadata=metadata_99,
                               taxonomy=taxonomy_93,sequences=zotu_93)
physeq_93
physeq_94 <- generate_phyloseq(otu=otutable_94,metadata=metadata_99,
                               taxonomy=taxonomy_94,sequences=zotu_94)
physeq_94
physeq_95 <- generate_phyloseq(otu=otutable_95,metadata=metadata_99,
                               taxonomy=taxonomy_95,sequences=zotu_95)
physeq_95
physeq_96 <- generate_phyloseq(otu=otutable_96,metadata=metadata_99,
                               taxonomy=taxonomy_96,sequences=zotu_96)
physeq_96
physeq_97 <- generate_phyloseq(otu=otutable_97,metadata=metadata_99,
                               taxonomy=taxonomy_97,sequences=zotu_97)
physeq_97

writeXStringSet(
  physeq_97@refseq,
  filepath = file.path(data_path, "datasets/otus_97_filtered.fasta"),
  format = "fasta"
)

physeq_98 <- generate_phyloseq(otu=otutable_98,metadata=metadata_99,
                               taxonomy=taxonomy_98,sequences=zotu_98)
physeq_98

# This is the selected dataset and analyzed in details
physeq_99 <- generate_phyloseq(otu=otutable_99,metadata=metadata_99,
                               taxonomy=taxonomy_99,sequences=zotu_99)
physeq_99

physeq_100 <- generate_phyloseq(otu=otutable_100,metadata=metadata_99,
                                taxonomy=taxonomy_100,sequences=zotu_100)
physeq_100

# Multiple rarefaction ---------------------------------------------------------

# NOTE. Because the goal is to compare the effect of clustering resolution (90–100%)
# on β-diversity and the nitrogen effect, the key requirement is keeping the 
# sampling process identical across datasets. Otherwise differences you observe 
# may come from sequencing depth artifacts rather than clustering resolution. For
# these reasons the most robust strategy is to rarefy all datasets to the exact 
# same depth using the same subset of samples.

# Rarefaction through our R package BRCore
pak::pak("germs-lab/BRCore")
#pak::pak("germs-lab/BRCore@22d4b66")
library(BRCore)

sort(sample_sums(physeq_90))
hist(sort(sample_sums(physeq_90)))

# Select same samples and equal depth!

merge_otu_tables <- function(){
  
  # 1. Put your dataframes into a named list
  # (Make sure the names are strings like "t90", "t91")
  otu_tab_list <- list(
    "t90"  = otutable_90,  "t91"  = otutable_91,  "t92"  = otutable_92,
    "t93"  = otutable_93,  "t94"  = otutable_94,  "t95"  = otutable_95,
    "t96"  = otutable_96,  "t97"  = otutable_97,  "t98"  = otutable_98,
    "t99"  = otutable_99,  "t100" = otutable_100
  )
  
  # 2. Extract sums and join
  master_depths <- names(otu_tab_list) %>%
    map(function(name) {
      df <- otu_tab_list[[name]]
      
      # Calculate column sums (Total reads per sample)
      data.frame(reads = colSums(df)) %>%
        # Rename the 'reads' column to the threshold name (e.g., t90)
        rename(!!name := reads) %>%
        rownames_to_column("sample_id")
    }) %>%
    reduce(full_join, by = "sample_id") %>%
    # Removing the control sample
    filter(sample_id != "MMPRNTCtrl7BB5") %>% 
    column_to_rownames("sample_id") %>% 
    arrange(t90)

  return(master_depths)
}

dim(merge_otu_tables())
merge_otu_tables()


# Test 1. Filter the data to same samples and rarefy at 600 sequences per sample.
sample_to_keep <- 
  merge_otu_tables() %>% filter(t90 > 600) %>% rownames()
  
# sample removed 
merge_otu_tables() %>% filter(t90 < 600) %>% rownames() 

prune_samples(samples = sample_to_keep, physeq_90)

physeq_90_rare <- multi_rarefy(
  physeq = prune_samples(samples = sample_to_keep, physeq_90),
  depth_level = 600,
  num_iter = 100,
  threads = 8,
  set_seed = 270226
) %>%
  update_otu_table(
    otu_rare = .,
    physeq =  prune_samples(samples = sample_to_keep, physeq_90)
  )

otu_table(physeq_90)
otu_table(physeq_90_rare)

physeq_91_rare <- multi_rarefy(
  physeq = prune_samples(samples = sample_to_keep, physeq_91),
  depth_level = 600,
  num_iter = 100,
  threads = 8,
  set_seed = 270226
) %>%
  update_otu_table(
    otu_rare = .,
    physeq =  prune_samples(samples = sample_to_keep, physeq_91)
  )

otu_table(physeq_91)
otu_table(physeq_91_rare)

physeq_92_rare <- multi_rarefy(
  physeq = prune_samples(samples = sample_to_keep, physeq_92),
  depth_level = 600,
  num_iter = 100,
  threads = 8,
  set_seed = 270226
) %>%
  update_otu_table(
    otu_rare = .,
    physeq =  prune_samples(samples = sample_to_keep, physeq_92)
  )

otu_table(physeq_92)
otu_table(physeq_92_rare)

physeq_93_rare <- multi_rarefy(
  physeq = prune_samples(samples = sample_to_keep, physeq_93),
  depth_level = 600,
  num_iter = 100,
  threads = 8,
  set_seed = 270226
) %>%
  update_otu_table(
    otu_rare = .,
    physeq =  prune_samples(samples = sample_to_keep, physeq_93)
  )

otu_table(physeq_93)
otu_table(physeq_93_rare)

physeq_94_rare <- multi_rarefy(
  physeq = prune_samples(samples = sample_to_keep, physeq_94),
  depth_level = 600,
  num_iter = 100,
  threads = 8,
  set_seed = 270226
) %>%
  update_otu_table(
    otu_rare = .,
    physeq =  prune_samples(samples = sample_to_keep, physeq_94)
  )

otu_table(physeq_94)
otu_table(physeq_94_rare)


physeq_95_rare <- multi_rarefy(
  physeq = prune_samples(samples = sample_to_keep, physeq_95),
  depth_level = 600,
  num_iter = 100,
  threads = 8,
  set_seed = 270226
) %>%
  update_otu_table(
    otu_rare = .,
    physeq =  prune_samples(samples = sample_to_keep, physeq_95)
  )

otu_table(physeq_95)
otu_table(physeq_95_rare)

physeq_96_rare <- multi_rarefy(
  physeq = prune_samples(samples = sample_to_keep, physeq_96),
  depth_level = 600,
  num_iter = 100,
  threads = 8,
  set_seed = 270226
) %>%
  update_otu_table(
    otu_rare = .,
    physeq =  prune_samples(samples = sample_to_keep, physeq_96)
  )

otu_table(physeq_96)
otu_table(physeq_96_rare)

physeq_97_rare <- multi_rarefy(
  physeq = prune_samples(samples = sample_to_keep, physeq_97),
  depth_level = 600,
  num_iter = 100,
  threads = 8,
  set_seed = 270226
) %>%
  update_otu_table(
    otu_rare = .,
    physeq =  prune_samples(samples = sample_to_keep, physeq_97)
  )

otu_table(physeq_97)
otu_table(physeq_97_rare)

physeq_98_rare <- multi_rarefy(
  physeq = prune_samples(samples = sample_to_keep, physeq_98),
  depth_level = 600,
  num_iter = 100,
  threads = 8,
  set_seed = 270226
) %>%
  update_otu_table(
    otu_rare = .,
    physeq =  prune_samples(samples = sample_to_keep, physeq_98)
  )

otu_table(physeq_98)
otu_table(physeq_98_rare)


physeq_99_rare <- multi_rarefy(
  physeq = prune_samples(samples = sample_to_keep, physeq_99),
  depth_level = 600,
  num_iter = 100,
  threads = 8,
  set_seed = 270226
) %>%
  update_otu_table(
    otu_rare = .,
    physeq =  prune_samples(samples = sample_to_keep, physeq_99)
  )

otu_table(physeq_99)
otu_table(physeq_99_rare)

physeq_100_rare <- multi_rarefy(
  physeq = prune_samples(samples = sample_to_keep, physeq_100),
  depth_level = 600,
  num_iter = 100,
  threads = 8,
  set_seed = 270226
) %>%
  update_otu_table(
    otu_rare = .,
    physeq =  prune_samples(samples = sample_to_keep, physeq_100)
  )

otu_table(physeq_100)
otu_table(physeq_100_rare)


# Adjusting the metadata ------------------------------------------------------- 

metadata_filt <- 
  metadata_99 %>% 
  dplyr::select(site_id, fert_status, plot_rep) %>% 
  rownames_to_column("sample_id") %>% 
  filter(sample_id %in% sample_to_keep) %>%
  mutate(
    site = site_id, 
    site = factor(recode(
      site,LUX = "Lux Harbor",LC = "Lake City",HAN = "Hancock",
      RHN = "Rhinelander",ESC = "Escanaba")),
    fert_status = factor(recode(
      fert_status, "FERT" = "Fertilized" ,  "UNFERT" ="Control"))
    ) 

metadata_filt

# Generating a legend for site to plot -----------------------------------------
pcoa_legend <-
  as_ggplot(get_legend(
    plot_ordination(
      ord = generate_ordination(ps = physeq_90_rare, method = "PCOA"),
      meta = metadata_filt,
      col_var = "site",
      shape_var = "fert_status",
      legend_inside = FALSE,
      ellipse = FALSE
    ) +
      scale_color_manual(values = palette_site) + 
      theme(legend.text = element_markdown(size = 10))
  ))

pcoa_legend

# Plotting 90% to ASV ordinations ----------------------------------------------

plot_beta_percent_sites <-
  ggarrange(
    plot_ordination(
      ord = generate_ordination(ps = physeq_90_rare, method = "PCOA"),
      meta = metadata_filt,col_var = "site", shape_var = "fert_status",
      legend_inside = FALSE,ellipse = FALSE) +
      scale_color_manual(values = palette_site) +
      labs(title="90% OTUs"),
    plot_ordination(
      ord = generate_ordination(ps = physeq_91_rare, method = "PCOA"),
      meta = metadata_filt,col_var = "site", shape_var = "fert_status",
      legend_inside = FALSE,ellipse = FALSE) +
      scale_color_manual(values = palette_site) +
      labs(title="91% OTUs"),
    plot_ordination(
      ord = generate_ordination(ps = physeq_92_rare, method = "PCOA"),
      meta = metadata_filt,col_var = "site", shape_var = "fert_status",
      legend_inside = FALSE,ellipse = FALSE) +
      scale_color_manual(values = palette_site) +
      labs(title="92% OTUs"),
    plot_ordination(
      ord = generate_ordination(ps = physeq_93_rare, method = "PCOA"),
      meta = metadata_filt,col_var = "site", shape_var = "fert_status",
      legend_inside = FALSE,ellipse = FALSE) +
      scale_color_manual(values = palette_site) +
      labs(title="93% OTUs"),
    plot_ordination(
      ord = generate_ordination(ps = physeq_94_rare, method = "PCOA"),
      meta = metadata_filt,col_var = "site", shape_var = "fert_status",
      legend_inside = FALSE,ellipse = FALSE) +
      scale_color_manual(values = palette_site) +
      labs(title="94% OTUs"),
    plot_ordination(
      ord = generate_ordination(ps = physeq_95_rare, method = "PCOA"),
      meta = metadata_filt,col_var = "site", shape_var = "fert_status",
      legend_inside = FALSE,ellipse = FALSE) +
      scale_color_manual(values = palette_site) +
      labs(title="95% OTUs"),
    plot_ordination(
      ord = generate_ordination(ps = physeq_96_rare, method = "PCOA"),
      meta = metadata_filt,col_var = "site", shape_var = "fert_status",
      legend_inside = FALSE,ellipse = FALSE) +
      scale_color_manual(values = palette_site) +
      labs(title="96% OTUs"),
    plot_ordination(
      ord = generate_ordination(ps = physeq_97_rare, method = "PCOA"),
      meta = metadata_filt,col_var = "site", shape_var = "fert_status",
      legend_inside = FALSE,ellipse = FALSE) +
      scale_color_manual(values = palette_site) +
      labs(title="97% OTUs"),
    plot_ordination(
      ord = generate_ordination(ps = physeq_98_rare, method = "PCOA"),
      meta = metadata_filt,col_var = "site", shape_var = "fert_status",
      legend_inside = FALSE,ellipse = FALSE) +
      scale_color_manual(values = palette_site) +
      labs(title="98% OTUs"),
    plot_ordination(
      ord = generate_ordination(ps = physeq_99_rare, method = "PCOA"),
      meta = metadata_filt,col_var = "site", shape_var = "fert_status",
      legend_inside = FALSE,ellipse = FALSE) +
      scale_color_manual(values = palette_site) +
      labs(title="99% OTUs"),
    plot_ordination(
      ord = generate_ordination(ps = physeq_100_rare, method = "PCOA"),
      meta = metadata_filt,col_var = "site", shape_var = "fert_status",
      legend_inside = FALSE, ellipse = FALSE) +
      scale_color_manual(values = palette_site) +
      labs(title="ASVs"),
    pcoa_legend,
    ncol=6,
    nrow=2,
    legend= "none")

plot_beta_percent_sites

# ***** FIGURE 3 - Beta-diversity ***** ----------------------------------------
ggsave(
  file.path(data_path, "results/Fig_3_plot_beta_percent_sites.pdf"),
  plot = ggpubr::annotate_figure(
    plot_beta_percent_sites,
    top = text_grob("BETA DIVERSITY AT 90% to ASVs SEQUENCE SIMILARITY", size = 12, face = "bold")
  ),
  device = "pdf"
)


# This is from a phyloseq object
plot_beta_percent_sites <-
  ggarrange(
    plot_ordination_from_phyloseq(ps = physeq_100_rare, 
                  col_var = "site_id", shape_var = "fert_status") +
      scale_color_manual(values = palette_site)+
      labs(title="ASVs"),
    plot_ordination_from_phyloseq(ps = physeq_99_rare, 
                  col_var = "site_id", shape_var = "fert_status")+
      scale_color_manual(values = palette_site)+
      labs(title="99% OTUs"),
    plot_ordination_from_phyloseq(ps = physeq_98_rare, 
                   col_var = "site_id", shape_var = "fert_status")+
      scale_color_manual(values = palette_site)+
      labs(title="OTUs 98%"),
    plot_ordination_from_phyloseq(ps = physeq_97_rare,
                   col_var = "site_id", shape_var = "fert_status")+
      scale_color_manual(values = palette_site)+
      labs(title="OTUs 97%"),
    plot_ordination_from_phyloseq(ps = physeq_96_rare,
                   col_var = "site_id", shape_var = "fert_status")+
      scale_color_manual(values = palette_site)+
      labs(title="OTUs 96%"),
    plot_ordination_from_phyloseq(ps = physeq_95_rare,
                   col_var = "site_id", shape_var = "fert_status")+
      scale_color_manual(values = palette_site)+
      labs(title="OTUs 95%"),
    plot_ordination_from_phyloseq(ps = physeq_94_rare,
                   col_var = "site_id", shape_var = "fert_status")+
      scale_color_manual(values = palette_site)+
      labs(title="OTUs 94%"),
    plot_ordination_from_phyloseq(ps = physeq_93_rare,
                   col_var = "site_id", shape_var = "fert_status")+
      scale_color_manual(values = palette_site)+
      labs(title="OTUs 93%"),
    plot_ordination_from_phyloseq(ps = physeq_92_rare,
                  col_var = "site_id", shape_var = "fert_status")+
      scale_color_manual(values = palette_site)+
      labs(title="OTUs 92%"),
    plot_ordination_from_phyloseq(ps = physeq_91_rare,
                  col_var = "site_id", shape_var = "fert_status")+
      scale_color_manual(values = palette_site)+
      labs(title="OTUs 91%"),
    plot_ordination_from_phyloseq(ps = physeq_90_rare,
                   col_var = "site_id", shape_var = "fert_status")+
      scale_color_manual(values = palette_site)+
      labs(title="OTUs 90%"),
    pcoa_legend,
    ncol=3,
    nrow=4,
    common.legend=FALSE,
    legend="none")

plot_beta_percent_sites

# **********************************************************************--------
# Extracting PERMANOVA results -------------------------------------------------

adonis_betadisp_90to100 <- rbind(
  extract_adonis(ps = physeq_90_rare),
  extract_adonis(ps = physeq_91_rare),
  extract_adonis(ps = physeq_92_rare),
  extract_adonis(ps = physeq_93_rare),
  extract_adonis(ps = physeq_94_rare),
  extract_adonis(ps = physeq_95_rare),
  extract_adonis(ps = physeq_96_rare),
  extract_adonis(ps = physeq_97_rare),
  extract_adonis(ps = physeq_98_rare),
  extract_adonis(ps = physeq_99_rare),
  extract_adonis(ps = physeq_100_rare)
)

adonis_betadisp_90to100


# **********************************************************************--------
# ASSESS CLUSTERING THRESHOLD --------------------------------------------------

# Looking for species with S_score of 1 (or lower if we want to test other taxa). 
# That means all the BLAST hits were the same species. 

as.data.frame(as.matrix(physeq_100_rare@tax_table)) %>% 
  filter(S_score > 0.9999) %>% 
  count(Species)

# Selecting Paraglomus brasilianum across datasets and counting
assess_clustering_threshold(
  tax_rank = "Species",
  taxon_name =  "Paraglomus brasilianum",
  ps_list = list("90" = physeq_90_rare,"91" = physeq_91_rare,"92" = physeq_92_rare,
                 "93" = physeq_93_rare,"94" = physeq_94_rare,"95" = physeq_95_rare,
                 "96" = physeq_96_rare,"97" = physeq_97_rare,"98" = physeq_98_rare,
                 "99" = physeq_99_rare,"100" = physeq_100_rare),
  score_col = "S_score",
  min_score = 0.99999999)

assess_taxon_sites(
  taxon_name =  "Paraglomus brasilianum",
  ps_list = list(
    "90" = physeq_90_rare, "91" = physeq_91_rare, "92" = physeq_92_rare,
    "93" = physeq_93_rare, "94" = physeq_94_rare, "95" = physeq_95_rare,
    "96" = physeq_96_rare, "97" = physeq_97_rare, "98" = physeq_98_rare,
    "99" = physeq_99_rare, "100" = physeq_100_rare
  ),
  site_var = "site_id"
)

# WARNING. Let's now check at genus level, as I believe that clustering may generate 
# sequences that are not necessarily classified at the same rank!

assess_clustering_threshold(
  tax_rank = "Genus",
  taxon_name =  "Paraglomus",
  ps_list = list("90" = physeq_90_rare,"91" = physeq_91_rare,"92" = physeq_92_rare,
                 "93" = physeq_93_rare,"94" = physeq_94_rare,"95" = physeq_95_rare,
                 "96" = physeq_96_rare,"97" = physeq_97_rare,"98" = physeq_98_rare,
                 "99" = physeq_99_rare,"100" = physeq_100_rare),
  score_col = "S_score",
  min_score = 0.1)

assess_taxon_sites(
  tax_rank = "Genus",
  taxon_name =  "Paraglomus",
  ps_list = list(
    "90" = physeq_90_rare, "91" = physeq_91_rare, "92" = physeq_92_rare,
    "93" = physeq_93_rare, "94" = physeq_94_rare, "95" = physeq_95_rare,
    "96" = physeq_96_rare, "97" = physeq_97_rare, "98" = physeq_98_rare,
    "99" = physeq_99_rare, "100" = physeq_100_rare
  ),
  site_var = "site_id"
)


# Selecting Cetraspora gilmorei
assess_clustering_threshold(
  tax_rank = "Species",
  taxon_name =  "Cetraspora gilmorei",
  ps_list = list("90" = physeq_90_rare,"91" = physeq_91_rare,"92" = physeq_92_rare,
                 "93" = physeq_93_rare,"94" = physeq_94_rare,"95" = physeq_95_rare,
                 "96" = physeq_96_rare,"97" = physeq_97_rare,"98" = physeq_98_rare,
                 "99" = physeq_99_rare,"100" = physeq_100_rare),
  score_col = "S_score",
  min_score = 0.1)

as.data.frame(as.matrix(physeq_100_rare@tax_table)) %>% 
  filter(S_score > 0.8) %>% 
  count(Genus)

assess_clustering_threshold(
  tax_rank = "Genus",
  taxon_name =  "Cetraspora",
  ps_list = list("90" = physeq_90_rare,"91" = physeq_91_rare,"92" = physeq_92_rare,
                 "93" = physeq_93_rare,"94" = physeq_94_rare,"95" = physeq_95_rare,
                 "96" = physeq_96_rare,"97" = physeq_97_rare,"98" = physeq_98_rare,
                 "99" = physeq_99_rare,"100" = physeq_100_rare),
  score_col = "S_score",
  min_score = 0.1)


# Selecting Paraglomus laccatum
assess_clustering_threshold(
  tax_rank = "Species",
  taxon_name =  "Paraglomus laccatum",
  ps_list = list("90" = physeq_90_rare,"91" = physeq_91_rare,"92" = physeq_92_rare,
                 "93" = physeq_93_rare,"94" = physeq_94_rare,"95" = physeq_95_rare,
                 "96" = physeq_96_rare,"97" = physeq_97_rare,"98" = physeq_98_rare,
                 "99" = physeq_99_rare,"100" = physeq_100_rare),
  score_col = "S_score",
  min_score = 0.1)

assess_taxon_sites(
  taxon_name =  "Paraglomus laccatum",
  ps_list = list(
    "90" = physeq_90_rare, "91" = physeq_91_rare, "92" = physeq_92_rare,
    "93" = physeq_93_rare, "94" = physeq_94_rare, "95" = physeq_95_rare,
    "96" = physeq_96_rare, "97" = physeq_97_rare, "98" = physeq_98_rare,
    "99" = physeq_99_rare, "100" = physeq_100_rare
  ),
  site_var = "site_id"
)

# Selecting Glomus macrocarpum
assess_clustering_threshold(
  tax_rank = "Species",
  taxon_name =  "Glomus macrocarpum",
  ps_list = list("90" = physeq_90_rare,"91" = physeq_91_rare,"92" = physeq_92_rare,
                 "93" = physeq_93_rare,"94" = physeq_94_rare,"95" = physeq_95_rare,
                 "96" = physeq_96_rare,"97" = physeq_97_rare,"98" = physeq_98_rare,
                 "99" = physeq_99_rare,"100" = physeq_100_rare),
  score_col = "S_score",
  min_score = 0.1)

# INTERPRETATION. T
# Things to consider for interpretation
# _ High ASV counts at 100% are normal — don’t interpret as hundreds of species.
# _ OTU counts at 95–96% often approximate species-level richness.
# _ Compare across taxa — some taxa are more diverse (many ASVs) than others.
# _ Check site distribution — some OTUs may appear only in certain sites after
#   clustering, which affects beta-diversity and PERMANOVA results.












# plotting

plot_beta_percent_sites <-
  ggarrange(
    plot_pcoa(physeq = physeq_100_rare, Col = "site_id", She = "fert_status") +
      labs(title="ASVs"),
    plot_pcoa(physeq = physeq_99_rare, Col = "site_id", She = "fert_status")+
      labs(title="99% OTUs"),
    plot_pcoa(physeq = physeq_98_rare, Col = "site_id", She = "fert_status")+
      labs(title="OTUs 98%"),
    plot_pcoa(physeq = physeq_97_rare, Col = "site_id", She = "fert_status")+
      labs(title="OTUs 97%"),
    plot_pcoa(physeq = physeq_96_rare, Col = "site_id", She = "fert_status")+
      labs(title="OTUs 96%"),
    plot_pcoa(physeq = physeq_95_rare, Col = "site_id", She = "fert_status")+
      labs(title="OTUs 95%"),
    plot_pcoa(physeq = physeq_94_rare, Col = "site_id", She = "fert_status")+
      labs(title="OTUs 94%"),
    plot_pcoa(physeq = physeq_93_rare, Col = "site_id", She = "fert_status")+
      labs(title="OTUs 93%"),
    plot_pcoa(physeq = physeq_92_rare, Col = "site_id", She = "fert_status")+
      labs(title="OTUs 92%"),
    plot_pcoa(physeq = physeq_91_rare, Col = "site_id", She = "fert_status")+
      labs(title="OTUs 91%"),
    plot_pcoa(physeq = physeq_90_rare, Col = "site_id", She = "fert_status")+
      labs(title="OTUs 90%"),
    ncol=3,
    nrow=4,
    common.legend=TRUE,
    legend="bottom")

plot_beta_percent_sites


plot_beta_percent_species <-
  ggarrange(
    plot_pcoa(physeq = physeq_100_rare, Col = "site_id", She = "fert_status",
             type = "species") +
      labs(title="ASVs"),
    plot_pcoa(physeq = physeq_99_rare, Col = "site_id", She = "fert_status",
             type = "species")+
      labs(title="OTUs 99%"),
    plot_pcoa(physeq = physeq_98_rare, Col = "site_id", She = "fert_status",
             type = "species")+
      labs(title="OTUs 98%"),
    plot_pcoa(physeq = physeq_97_rare, Col = "site_id", She = "fert_status",
             type = "species")+
      labs(title="OTUs 97%"),
    plot_pcoa(physeq = physeq_96_rare, Col = "site_id", She = "fert_status",
             type = "species")+
      labs(title="OTUs 96%"),
    plot_pcoa(physeq = physeq_95_rare, Col = "site_id", She = "fert_status",
             type = "species")+
      labs(title="OTUs 95%"),
    plot_pcoa(physeq = physeq_94_rare, Col = "site_id", She = "fert_status",
             type = "species")+
      labs(title="OTUs 94%"),
    plot_pcoa(physeq = physeq_93_rare, Col = "site_id", She = "fert_status",
             type = "species")+
      labs(title="OTUs 93%"),
    plot_pcoa(physeq = physeq_92_rare, Col = "site_id", She = "fert_status",
             type = "species")+
      labs(title="OTUs 92%"),
    plot_pcoa(physeq = physeq_91_rare, Col = "site_id", She = "fert_status",
             type = "species")+
      labs(title="OTUs 91%"),
    plot_pcoa(physeq = physeq_90_rare, Col = "site_id", She = "fert_status",
             type = "species")+
      labs(title="OTUs 90%"),
    ncol=3,
    nrow=4,
    common.legend=TRUE,
    legend="bottom")

plot_beta_percent_species

# INTERPRETATION. 
# PCoA axes report proportions of variance, not absolute ecological signal. When 
# you change clustering resolution, you change the structure of the distance 
# matrix, which directly affects those proportions. If the total variation shrinks
# or redistributes, the percent explained by early axes can increase, even if 
# the biological signal did not strengthen.




# With this analysis we identified a resolution that preserves ecological signal 
# without over-fragmenting taxa.


# https://maarjam.ut.ee/



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



metadata_AMF_rare <- 
  as.data.frame(as.matrix(sample_data(physeq_AMF_rare))) %>% 
  mutate(site_name = gsub(" ", "", site)) %>% 
  mutate(sample_id = paste(site_name, fert_status, plot_rep , sep="_")) %>% 
  mutate(across(.cols = c(3:7, 9:13), .fns = as.numeric))

head(metadata_AMF_rare)




mean_metadata_AMF_rare <- 
  as(sample_data(mean_physeq_AMF_rare), "matrix") %>% 
  as.data.frame() %>% 
  dplyr::select(site, fert_status, plot_rep) %>%
  mutate(sample_id = rownames(.)) %>% 
  #rownames_to_column("sample_id") %>% 
  separate(sample_id, into = c("site", "fert_status", "plot_rep"), 
           sep = "_", remove = FALSE)

head(mean_metadata_AMF_rare)


mean_physeq_AMF_rare@tax_table %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  pull(Class) %>% 
  unique()

mean_physeq_AMF_rare@tax_table %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  pull(Genus) %>% 
  unique()

mean_physeq_AMF_rare@tax_table %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  subset(Genus %in% c("Mortierella", "Jimgerdemannia"))


mlin2_fert_status <- 
  Maaslin2(
    cores = 4,
    output = "demo_output", 
    min_abundance = 0.001,
    min_prevalence = 0.1,
    min_variance = 0.001,
    input_data = otutable_AMF_rare, 
    input_metadata = metadata_AMF_rare %>% dplyr::select(fert_status, site),  
    #fixed_effects  = c("fert_status"),  # fert_status as fixed 
    reference = c("fert_status,Control"),
    random_effects = c("site"),         # site as random
    correction = "BH",
    normalization = "NONE",
    standardize = FALSE, 
    plot_scatter = FALSE)

str(mlin2_fert_status)

mlin2_fert_status$results %>% filter(qval <= 0.05)

mlin2_fert_status$results %>% 
  filter(metadata == "site")%>% 
  filter(qval <= 0.05)

# INETERPRETATION. When site is included as random effect, its variation is absorbed
# into the random effect structure. The model accounts for between-site variance 
# when estimating coefficients. This reduces residual variance, making it harder 
# to detect site effects as fixed effects. The 16 OTUs that survive are those with
# very strong site signals that persist even after the random effect soaks up site
# variance.

mlin2_site <- 
  Maaslin2(
    cores = 4,
    output = "demo_output_site", 
    min_abundance = 0.001,
    min_prevalence = 0.1,
    min_variance = 0.001,
    input_data = otutable_AMF_rare, 
    input_metadata = metadata_AMF_rare %>% dplyr::select(site),
    fixed_effects  = c("site"),   # site as fixed only
    # no random_effects
    reference = c("site,Escanaba"),  # no space
    correction = "BH", # correction = "bonferroni",
    normalization = 'NONE',
    standardize = FALSE, 
    plot_scatter = FALSE)

str(mlin2_site)

mlin2_site$results %>% 
  filter(metadata == "site") %>% 
  filter(qval <= 0.05)


# Check effect sizes
mlin2_site$results %>%
  dplyr::filter(metadata == "site") %>%
  dplyr::arrange(desc(abs(coef)), qval) %>%
  head(40)






# Plotting 
mlin2_site$results %>%
  dplyr::filter(metadata == "site", qval <= 0.05) %>%
  left_join(as.data.frame(as.matrix(physeq_AMF_rare@tax_table)) %>%
              rownames_to_column("feature"), 
            by = "feature") %>% 
  ggplot(aes(x = value, y = BestMatch, fill = coef)) + 
  geom_tile() +
  geom_text(aes(label = round(coef, 2)), size = 4, show.legend = FALSE) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "white", color = "white"),
    plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5, vjust = 0.5),
    plot.subtitle = element_markdown(size = 10, face = "bold", hjust = 0.5, vjust = 0.5),
    axis.text.x = element_text(angle = 33, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_markdown(size = 10),
    axis.title = element_markdown(size = 10), 
    legend.title = element_text(size = 10, face = "bold"), 
    legend.text = element_text(size = 9),
    legend.key.height = unit(0.4, "cm"),  
    legend.key.width = unit(0.4, "cm"),
    legend.position = "right"
  ) +
  scale_fill_gradient2(name = "lfc",
                       low = "blue", high = "red", mid = "#FFFFCC", na.value = "white",
                       limits = c(-8, 8),
                       breaks = c(-8, -4, 0, 4, 8)) 



fitGLM_yield_fixslope <- 
  glmmTMB(dry_matter ~ fert_status + p_ppm + fert_status:p_ppm + (1 | site),
          dispformula = ~site,   # allows variance to differ by site
          data = lmer_otu_chem_data)

summary(fitGLM_yield_fixslope)

diagnostics_dharma(
  model      = fitGLM_yield_fixslope,
  group_var1 = lmer_otu_chem_data$site,
  group_var2 = NULL
)

# Trying gamma function
fitGLM_yield_fixslope_gamma <- 
  glmmTMB(dry_matter ~ fert_status + p_ppm + fert_status:p_ppm + (1 | site),
          dispformula = ~site,   # allows variance to differ by site
          family = Gamma(link = "log"),
          data = lmer_otu_chem_data, 
          REML = TRUE)

summary(fitGLM_yield_fixslope)

diagnostics_dharma(
  model      = fitGLM_yield_fixslope,
  group_var1 = lmer_otu_chem_data$site,
  group_var2 = NULL
)

anova(fitGLM_yield_fixslope, fitGLM_yield_fixslope_gamma)

# INTERPRETATION. Think of “family” as the distribution of Y and of “link” as 
# the transformation of the mean of Y that makes the relationship with 
# predictors linear.

# _ Family: which probability distribution you assume for the response (Gaussian, 
#   Binomial, Poisson, Gamma, etc.). It controls the mean–variance relationship 
#   and the likelihood.
# _ Link: the function that maps the mean onto the real line so a linear predictor
#   makes sense, and it determines how you interpret coefficients (differences, 
#   ratios, odds ratios).





# INTERPRETATION. A Gaussian mixed-effects model was used to analyze dry matter 
# yield, including fertilization status and soil phosphorus as fixed effects, 
# and a random intercept for site to account for baseline differences among 
# locations. Additionally, the model allowed residual variance to differ by site 
# to accommodate heteroskedasticity. Fertilization had a strong positive effect 
# on yield, increasing dry matter by approximately 2.77 units (SE = 0.53, p < 0.001), 
# whereas soil phosphorus concentration was not significantly associated with 
# yield (p = 0.52). The random intercept for site indicated substantial between-site 
# variation (SD = 1.46). Site-specific residual variance was significantly different
# among some locations, reflecting variation in within-site variability. 
# DHARMa diagnostics showed no evidence of over- or under-dispersion, outliers, 
# or severe deviations from model assumptions, indicating that the model provides
# a reliable description of the data.

# glm model for hill_q0
fitGLM_hill_q0_fixslope <- 
  glmmTMB(hill_q0 ~ fert_status + p_ppm + fert_status:p_ppm + (1 | site),
          dispformula = ~site,   # allows variance to differ by site
          data = lmer_otu_chem_data)

summary(fitGLM_hill_q0_fixslope)

diagnostics_dharma(
  model      = fitGLM_hill_q0_fixslope,
  group_var1 = lmer_otu_chem_data$site)


# Impact of richness on yield --------------------------------------------------

# Model for diversity
fit_yield_hill_fixslope <- lmer(
  dry_matter ~ hill_q0 + p_ppm + hill_q0:p_ppm + (1 | site),
  data = lmer_otu_chem_data,
  REML = FALSE
)

summary(fit_yield_hill_fixslope)
anova(fit_yield_hill_fixslope)

lmer_otu_chem_data$hill_q0_sc <- scale(lmer_otu_chem_data$hill_q0)
lmer_otu_chem_data$p_ppm_sc  <- scale(lmer_otu_chem_data$p_ppm)

fit_yield_hill_fixslope <- lmer(
  dry_matter ~ hill_q0_sc + p_ppm_sc + hill_q0_sc:p_ppm_sc + (1 | site),
  data = lmer_otu_chem_data,
  REML = FALSE
)


diagnostics_dharma(
  model      = fit_yield_hill_fixslope,
  group_var1 = lmer_otu_chem_data$site)

# VISUALIZING THE MODELS -------------------------------------------------------

predict_response(fitGLM_yield_fixslope, terms = "fert_status")

predict_response(fitGLM_yield_fixslope, terms = "fert_status") %>% 
  plot() + 
  labs(title = "Effect of Fertilization on Dry Matter",
       y = "Predicted Dry Matter", 
       x = "Fertilization Status") +
  theme_minimal()

ggpredict(fitGLM_yield_fixslope, terms = "fert_status") %>% 
  ggplot(aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) +
  geom_jitter(data = lmer_otu_chem_data,
              aes(x = fert_status, y = dry_matter),
              width = 0.1, alpha = 0.5) +
  labs(x = "Fertilization", y = "Dry matter") +
  theme_bw()

plot_model(fitGLM_yield_fixslope, type = "re")

ggplot(lmer_otu_chem_data,
       aes(x = site, y = dry_matter)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5) +
  theme_bw()




# Load R packages --------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  # --- Project & Package Management ---
  renv,           # Reproducible environments (lockfile management)
  pak,            # Fast, modern package installation dependency engine
  cli,            # Attractive and informative command-line interfaces  
  styler,         # Automated R code formatting to tidyverse standards
  janitor,        # Data cleaning (heavily used for cleaning column names)
  magrittr,       # Pipe operators (%>%) and functional sequence aliases
  
  # --- Sequence Analysis & Phylogenetics ---
  Biostrings,     # Memory-efficient containers for DNA/RNA/AA sequences
  ape,            # Analysis of Phylogenetics and Evolution (core tree handling)
  msa,            # Multiple Sequence Alignment (Interface to Clustal, Muscle)
  DECIPHER,       # Tools for sequence alignment, chimera detection, and probes
  phangorn,       # Phylogenetic analysis (MP, ML, and distance-based methods)
  tidysq,         # Tidy processing of biological sequences
  tidytree,       # Tidyverse-style manipulation of phylogenetic tree objects
  ggtree,         # Visualization and annotation of phylogenetic trees
  
  # --- Microbial Ecology & Diversity ---
  decontam,       # Statistical identification/removal of contaminant OTUs
  phyloseq,       # Integration of OTU tables, taxonomy, and metadata
  speedyseq,      # High-performance optimizations for phyloseq functions
  vegan,          # Community ecology (Ordination, Alpha/Beta diversity)
  AICcPermanova,  # Model selection for PERMANOVA using AICc
  
  # --- Data Science & Visualization ---
  tidyverse,      # Core suite (ggplot2, dplyr, tidyr, purrr, etc.)
  ggtext,         # Markdown/HTML rendering for ggplot2 labels
  ggpubr,         # Publication-ready ggplot2 themes and figure arrangements
  cowplot,        # Powerful plot arrangement and theme adjustments
  gridExtra,      # Miscellaneous functions for "grid" graphics (tables/plots)
  ggrepel,        # Non-overlapping text labels for ggplot2
  scales,         # Internal scaling and axis transformations for plots
  ggh4x,          # Extensions for ggplot2 facets and axes (nested facets, etc.) 
  ggstar,         # Adding geometric shapes to ggplot
  ggtreeExtra,    # Annotation Phylogenetics Software Visualization
  ggnewscale,
  
  # --- Statistics & Modeling ---
  agricolae,      # Agricultural research stats (Tukey's, LSD, experimental design)
  BRCore,         # Custom/Internal core functions for standardized workflows
  lme4,           # Linear and Generalized Linear Mixed-Effects Models (GLMM)
  glmmTMB,        # Fast GLMMs with support for Zero-Inflation/Overdispersion
  robustlmm,      # Robust Linear Mixed-Effects Models (outlier resistance)
  DHARMa,         # Residual diagnostics for hierarchical (multi-level) models
  parallel,       # Parallel computing support for multi-core processing
  ggeffects,      # Estimated Marginal Means and marginal effects for models
  sjPlot,         # Visualizing Mixed-Effects models and diagnostic tables
  broom.mixed,    # Tidy summary outputs (tidiers) for mixed-effect models
  merTools,       # Tools for analyzing lme4/merMod objects (uncertainty/prediction)
  multcompView,   # Visualizing multiple comparison significance (compact letter display)
  gllvm,          # Generalized Linear Latent Variable Models (multivariate stats)
  Maaslin2,       # Multivariate Association with Linear Models (microbiome stats) 
  maaslin3,       # Multivariate Association with Linear Models (microbiome stats)  
  
  # --- Networks & Hierarchical Viz ---
  ggdendro,       # Extract/Plot dendrogram data using ggplot2
  ggraph,         # Visualization of networks, graphs, and complex trees
  
  install = FALSE
)

PlotMaaslin2 <- function(maaslin_results, 
                         physeq_object) {
  
  as.data.frame(as.matrix(physeq_object@tax_table)) %>%
    rownames_to_column("feature") %>% print()
  
  # Process the data: Remove duplicate comparisons
  cleaned_results <- 
    maaslin_results %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(value, reference_group)), collapse = "-")) %>%
    ungroup() %>%
    distinct(feature, pair, .keep_all = TRUE) %>%
    dplyr::mutate(pair = as.factor(pair)) %>% 
    left_join(as.data.frame(as.matrix(physeq_object@tax_table)) %>%
                rownames_to_column("feature"), by = "feature") %>%
    mutate(Genus_makrdown = str_c("*", str_trim(Genus), "*")) %>% 
    mutate(coef_log2   = coef / log(2),
           stderr_log2 = stderr / log(2))
  #mutate(Taxonomy = paste(Genus_makrdown, str_c("(", str_trim(feature), ")"))) %>% 
  #mutate(foldChange = exp(coef))
  #mutate(BestMatch_makrdown = str_c("*", str_trim(BestMatch), "*")) %>%
  #mutate(Taxonomy = paste(BestMatch_makrdown, str_c("(", str_trim(feature), ")"))) 
  #print(cleaned_results, n=19)
  print(range(cleaned_results$coef_log2))
  print(cleaned_results$Genus_makrdown)
  print(cleaned_results$Genus)
  
  # Generate the heatmap plot
  p <- 
    ggplot(cleaned_results, aes(x = pair, y = Genus_makrdown, fill = coef)) + 
    geom_tile() +
    geom_text(aes(label = round(coef, 2)), size = 4, show.legend = FALSE) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "white", color = "white"),
      plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5, vjust = 0.5),
      plot.subtitle = element_markdown(size = 10, hjust = 0.5, vjust = 0.5),
      axis.text.x = element_text(angle = 33, hjust = 1, vjust = 1, size = 10),
      axis.text.y = element_markdown(size = 10),
      axis.title = element_markdown(size = 10), 
      legend.title = element_text(size = 10, face = "bold"), 
      legend.text = element_text(size = 9),
      legend.key.height = unit(0.4, "cm"),  
      legend.key.width = unit(0.4, "cm"),
      legend.position = "right"
    ) +
    #facet_grid(~pair) +
    scale_fill_gradient2(name = "lfc",
                         low = "blue", high = "red", mid = "#FFFFCC", na.value = "white",
                         limits = c(-12.5, 12.5),
                         breaks = c(-12.5, -6.5, 0, 6.5, 12.5))
  
  return(p)
}

head(mlin2_site)

# Are there NAs introduced by the conversion?
mlin2_site %>%
  mutate(coef_log2 = coef / log(2)) %>%
  summarise(
    n_total    = n(),
    n_na_coef  = sum(is.na(coef)),
    n_na_log2  = sum(is.na(coef / log(2)))
  )

# Check if missing comparisons were present before filtering
mlin2_site %>%
  filter(qval <= 0.05) %>%
  distinct(value) %>%
  pull(value)

# Plotting
PlotMaaslin2(maaslin_results = mlin2_site %>% dplyr::filter(qval <= 0.05),
             physeq_object =  physeq_AMF_rare_Gen) +
  labs(title = "Differential abundance of genera among sites",
       subtitle = "lfc = log fold change, fold change = exp(lfc)<br>
       Negative lfc = more abundant in the reference site<br>
       Positive lfc = more abundant in the value site",
       x = "Pairwise site comparison (site vs reference site)",
       y = "Genus")




PlotMaaslin2(maaslin_results = mlin2_fert_status_gen$results,
             physeq_object = physeq_AMF_rare_Gen,
             value_type = "log2")



PlotMaaslin2(maaslin_results = mlin2_site %>% filter(qval <= 0.05),
             physeq_object =  physeq_AMF_rare_Gen,
             value_type= "log2") +
  labs(title = "Differential abundance of genera among sites",
       subtitle = "l2FC = log2 fold change, fold change = 2<sup>l2FC</sup><br>
                Positive l2FC = More abundant in the test site<br>
                Negative l2FC = More abundant in the reference site (i.e., reference)",
       x = "Pairwise site comparison",
       y = "Genus")


PlotMaaslin2(maaslin_results = mlin2_site 
             %>% filter(qval <= 0.05),
             physeq_object =  physeq_AMF_rare_Gen,
             value_type= "log2") +
  labs(title = "Differential abundance of genera among sites",
       subtitle = "l2FC = log2 fold change, fold change = 2 exp(l2FC)<br>
                   Positive l2FC = More abundant in the fist site<br>
                   Negative l2FC = more abundant in the second site (i.e. reference)",
       x = "Pairwise site comparison (site vs reference site)",
       y = "Genus")

PlotMaaslin2(maaslin_results = mlin2_site 
             %>% filter(qval <= 0.05),
             physeq_object =  physeq_AMF_rare_Gen,
             value_type= "ln") +
  labs(title = "Differential abundance of genera among sites",
       subtitle = "l2FC = log2 fold change, fold change = 2 exp(l2FC)<br>
       Negative lfc = more abundant in the reference site<br>
       Positive lfc = more abundant in the value site",
       x = "Pairwise site comparison (site vs reference site)",
       y = "Genus")


labs(
  title = "Differential abundance of genera among sites",
  subtitle = "Values represent log fold change (ln scale). Direction: first site relative to second (A/B)",
  x = "Pairwise site comparison",
  y = "Genus"
) 


# Plot check
mlin2_site %>% filter(value == "Escanaba", reference_group == "Hancock") %>% 
  filter(feature == "Zotu885") %>% pull(coef)
mlin2_site %>% filter(value == "Escanaba", reference_group == "Hancock") %>% 
  filter(feature == "Zotu9149") %>% pull(coef)



mlin2_fert_status_gen$results %>%
  dplyr::filter(metadata == "fert_status", qval <= 0.05)


# INTERPRETATION. Differential abundance analysis revealed strong site-specific 
# structuring of AMF genera, with several taxa showing large log fold changes 
# across pairwise site comparisons. Genera such as Septoglomus, Diversispora, 
# and Microkamienskia exhibited pronounced enrichment in specific sites, with 
# log fold changes exceeding 6 in some contrasts, indicating substantial shifts 
# in relative abundance. In contrast, genera such as Paraglomus and Oehlia tended
# to be consistently depleted in certain sites across multiple comparisons. 
# Overall, the magnitude and consistency of these patterns indicate that site-level
# environmental factors strongly influence AMF community composition.

# Figure X. Differential abundance of AMF genera across sites.
# Heatmap of log fold changes (lfc) from Maaslin2 analysis for pairwise site 
# comparisons. Values represent the natural log-transformed fold change in 
# genus-level abundance between sites (first site relative to second). Positive 
# values (red) indicate higher abundance in the first site, while negative values
# (blue) indicate higher abundance in the second site. Only significant 
# associations (q ≤ 0.05) are shown. The large magnitude of several log fold 
# changes indicates strong site-specific structuring of AMF communities.

# What if I want to use log2foldchange instead 


mlin2_site <- 
  run_maaslin2_all_pairs(
    physeq = 
      physeq_AMF_rare_Gen %>% 
      subset_taxa(Genus != "Unclassified") %>%
      prune_taxa(taxa_sums(.) > 0, .) %>% 
      prune_samples(sample_sums(.) > 0, .) %>% 
      # Add a pseudocount (e.g., 1) so log2(0+1) = 0
      transform_sample_counts(fun = function(x) log2(x + 1)),
    factor_variable = "site",
    qval_threshold = 0.05, 
    norm_method = "NONE", 
    trans_method = "NONE"
  )


PlotMaaslin2(maaslin_results = mlin2_site 
             %>% filter(qval <= 0.05),
             physeq_object =  physeq_AMF_rare_Gen) +
  labs(title = "Differential abundance of genera among sites",
       subtitle = "lfc = log fold change, fold change = exp(lfc)",
       x = "Pairwise site comparison",
       y = "Genus")



