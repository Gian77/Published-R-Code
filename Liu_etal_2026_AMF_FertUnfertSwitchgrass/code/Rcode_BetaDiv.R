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

# >>>>>>>>> BETA DIVERSITY <<<<<<<<<< ------------------------------------------

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
  file.path(data_path, "figures/Fig_2_betadiv.pdf"),
  plot = ggpubr::annotate_figure(
    Figure_2_beta_long,
    top = text_grob("BETA DIVERSITY", size = 12, face = "bold")
  ),
  device = "pdf"
)

# ***** FIGURE S4 - betadisper plot ***** --------------------------------------
ggsave(
  file.path(data_path, "figures/Fig_SX_betadiv_plot.pdf"),
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
  file.path(data_path, "figures/Fig_3_bar_charts_gen_relab.pdf"),
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





