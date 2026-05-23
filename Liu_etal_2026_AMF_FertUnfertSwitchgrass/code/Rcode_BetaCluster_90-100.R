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

# >>>>>>>>> PHYLOGENETIC TREE <<<<<<<<<< ---------------------------------------

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
  dplyr::select(site, fert_status, plot_rep) %>% 
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

# ***** FIGURE 3 - Beta-diversity ***** 
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
#Start with PERMANOVA R² comparison (#1) - it directly tests your hypothesis and 
# is straightforward to interpret. Supplement with ANOSIM (#2) for a non-parametric 
# confirmation. If you want to publish this observation, include the variance 
# partitioning (#5) to show how the relative importance of factors changes with 
# taxonomic resolution.
#The key result would be showing that R² (or other effect sizes) systematically 
# increases as clustering becomes more aggressive, which statistically confirms
# what you're seeing visually in the PCoA plots.

