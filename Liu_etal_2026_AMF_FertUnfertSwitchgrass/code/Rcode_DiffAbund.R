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

# >>>>>>>>> DIFFERENTIAL ABUNDANCE <<<<<<<<<< ----------------------------------

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




