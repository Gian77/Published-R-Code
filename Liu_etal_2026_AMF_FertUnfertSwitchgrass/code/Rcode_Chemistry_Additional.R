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
    top = text_grob("EFFECT OF NITROGEN ON SWITCHGRASS YIELD", size = 12, face = "bold")
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


