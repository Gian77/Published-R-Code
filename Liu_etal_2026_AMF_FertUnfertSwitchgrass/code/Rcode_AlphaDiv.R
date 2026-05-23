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

# >>>>>>>>> ALPHA DIVERSITY <<<<<<<<<< -----------------------------------------

# Load R packages. Define your package list as a character vector
pkgs <- c(
  # --- Project & Package Management ---
  "renv",            # Reproducible environments (lockfile management)
  "pak",             # Fast, modern package installation dependency engine
  "cli",             # Attractive and informative command-line interfaces  
  "styler",          # Automated R code formatting
  "janitor",         # Data cleaning
  "magrittr",        # Pipe operators (%>%)
  "parallel",        # {arallel processes}
  
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
# **** 1. ALPHA DIVERSITY **** -------------------------------------------------

# Attaching rarefaction metrics ------------------------------------------------

add_rarefaction_metrics(physeq_AMF_clean) %>% 
  sample_data() %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  arrange(read_num)


rarefaction_plot <- 
  add_rarefaction_metrics(physeq_AMF_clean) %>% 
  plot_rarefaction_metrics()

rarefaction_plot 

# ***** FIGURE S1 - identify depth ***** ---------------------------------------
ggsave(plot = rarefaction_plot, 
       path = data_path,
       filename = "figures/Fig_SX_identify_raredepth_plots.pdf")

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

# ***** FIGURE S2 - rarefaction curves ***** -----------------------------------
ggsave(plot = Fig_S2_rarecurve, 
       path = data_path,
       filename = "figures/Fig_SX_rarecurves.pdf")

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

# ***** FIGURE S3 - hill q=0 histograms ***** ----------------------------------
ggarrange(
alpha_df %>% 
  ggplot2::ggplot(aes(x = hill_q0)) +
  geom_histogram(binwidth = 15) +
  labs(title = "hill_q0") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(hill_q0)), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE),

alpha_df %>% 
  ggplot(aes(x = log(hill_q0))) +
  geom_histogram(binwidth = 0.1) +
  labs(title = "hill_q0") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(log(hill_q0))), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE),
ncol = 2)

# Model Comparison -------------------------------------------------------------
# Does fertilizer affect richness? H0: mean richness (FERT) = 0
summary(run_lmem(alpha_df, hill_q0, "baseline", reml = TRUE))
summary(run_lmem(alpha_df, hill_q0, "fixslope", reml = TRUE))
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

# ***** FIGURE S4 - hill q=1 histograms ***** ----------------------------------
ggarrange(
alpha_df %>% 
  ggplot2::ggplot(aes(x = hill_q1)) +
  geom_histogram(binwidth = 5) +
  labs(title = "hill_q1") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(hill_q1)), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE),

alpha_df %>% 
  ggplot(aes(x = log(hill_q1))) +
  geom_histogram(binwidth = 0.2) +
  labs(title = "hill_q1") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(log(hill_q1))), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE),
ncol = 2) 

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
# ***** FIGURE S5 - hill q=2 histograms ***** ----------------------------------

ggarrange(
alpha_df %>% 
  ggplot2::ggplot(aes(x = hill_q2)) +
  geom_histogram(binwidth = 5) +
  labs(title = "hill_q2") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(hill_q2)), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE),

alpha_df %>% 
  ggplot(aes(x = log(hill_q2))) +
  geom_histogram(binwidth = 0.3) +
  labs(title = "hill_q2") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(log(hill_q2))), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE),
ncol = 2) 

# Model Comparison -------------------------------------------------------------
# Does fertilizer affect Shannon? H0: mean Shannon (FERT) = 0
summary(run_lmem(alpha_df, hill_q2, "baseline", reml = TRUE))
summary(run_lmem(alpha_df, hill_q2, "fixslope", reml = TRUE))
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

# ***** FIGURE S6 - diagnostics model hill q=0 ***** ---------------------------
diagnostics_dharma(
  model = run_lmem(alpha_df, hill_q0, "baseline", reml = TRUE),
  group_var1 = alpha_df$site,
  group_var2 = alpha_df$plot_rep
)

diagnostics_dharma(
  model = run_lmem(alpha_df, hill_q0, "fixslope"),
  group_var1 = alpha_df$site,
  group_var2 = alpha_df$plot_rep
)

# ***** FIGURE S7 - diagnostics model hill q=2 ***** ---------------------------
diagnostics_dharma(
  model = run_lmem(alpha_df, hill_q2, "baseline", reml = TRUE),
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

# Color palettes ---------------------------------------------------------------

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
                     color = "black", linewidth = 4, shape = 18) +
        stat_summary(aes(group = 1), fun = mean, geom = "line", 
                     color = "black", linewidth = 1) +
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
                     color = "black", linewidth = 4, shape = 18) +
        stat_summary(aes(group = 1), fun = mean, geom = "line", 
                     color = "black", linewidth = 1) +
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
  file.path(data_path, "figures/Fig_1_alphadiv.pdf"),
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

