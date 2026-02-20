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
  lme4,
  glmmTMB,
  robustlmm,
  DHARMa,
  parallel,
  ggeffects,
  sjPlot,
  broom.mixed,
  merTools,
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

# Extract dataset
alpha_df <-
  physeq_AMF_rare@sam_data %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(across(.cols = c(5:13), .fns = as.numeric)) %>%
  mutate(
    site = factor(site),
    plot_rep = factor(paste("R", plot_rep, sep = "")),
    fert_status = factor(recode(
      fert_status, "FERT" = "Fertilized" ,  "UNFERT" ="Control"
    ))
  )

str(alpha_df)
head(alpha_df)

# hill_0 -----------------------------------------------------------------------

alpha_df %>% 
  ggplot2::ggplot(aes(x = hill_0)) +
  geom_histogram(binwidth = 20) +
  labs(title = "hill_0") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(hill_0)), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)

alpha_df %>% 
  ggplot(aes(x = log(hill_0))) +
  geom_histogram(binwidth = 0.1) +
  labs(title = "hill_0") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(log(hill_0))), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)

# Model Comparison -------------------------------------------------------------
# Does fertilizer affect richness? H0: mean richness (FERT) = 0
summary(run_lmem(alpha_df, hill_0, "baseline"))
summary(run_lmem(alpha_df, hill_0, "fixslope"))
summary(run_lmem(alpha_df, hill_0, "randomslope"))

# Does fertilizer improve fit?
anova( run_lmem(alpha_df, hill_0, "baseline"), 
       run_lmem(alpha_df, hill_0, "randomslope") )

# Does slope vary among sites?
anova( run_lmem(alpha_df, hill_0, "fixslope"), 
       run_lmem(alpha_df, hill_0, "randomslope") )

# INTERPRETATION. So the statistically supported model is the baseline model with 
# no fertilizer effect and no slope variation among sites. Richness varies strongly
# among plots and moderately among sites, but not due to fertilizer.

# hill_1 -----------------------------------------------------------------------

alpha_df %>% 
  ggplot2::ggplot(aes(x = hill_1)) +
  geom_histogram(binwidth = 15) +
  labs(title = "hill_1") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(hill_1)), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)

alpha_df %>% 
  ggplot(aes(x = log(hill_1))) +
  geom_histogram(binwidth = 0.1) +
  labs(title = "hill_1") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(log(hill_1))), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)

# Model Comparison -------------------------------------------------------------
# Does fertilizer affect Shannon? H0: mean Shannon (FERT) = 0
summary(run_lmem(alpha_df, hill_1, "baseline"))
summary(run_lmem(alpha_df, hill_1, "fixslope"))
summary(run_lmem(alpha_df, hill_1, "randomslope"))

# Does fertilizer improve fit?
anova( run_lmem(alpha_df, hill_1, "baseline"), 
       run_lmem(alpha_df, hill_1, "randomslope") )

# Does slope vary among sites?
anova( run_lmem(alpha_df, hill_1, "fixslope"), 
       run_lmem(alpha_df, hill_1, "randomslope") )

# The "Zero Variance" Culprit
# Look closely at the Random Effects section for all three models. You will see 
# this line: site (Intercept)  0.000  0.000
# A "singular fit" occurs when one of your random effects has a variance of exactly
# (or very near) zero.

# hill_2 -----------------------------------------------------------------------

alpha_df %>% 
  ggplot2::ggplot(aes(x = hill_2)) +
  geom_histogram(binwidth = 10) +
  labs(title = "hill_2") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(hill_2)), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)

alpha_df %>% 
  ggplot(aes(x = log(hill_2))) +
  geom_histogram(binwidth = 0.1) +
  labs(title = "hill_2") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(log(hill_2))), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)

# Model Comparison -------------------------------------------------------------
# Does fertilizer affect Shannon? H0: mean Shannon (FERT) = 0
summary(run_lmem(alpha_df, hill_2, "baseline"))
summary(run_lmem(alpha_df, hill_2, "fixslope"))
summary(run_lmem(alpha_df, hill_2, "randomslope"))

# Does fertilizer improve fit?
anova( run_lmem(alpha_df, hill_2, "baseline"), 
       run_lmem(alpha_df, hill_2, "randomslope") )

# Does slope vary among sites?
anova( run_lmem(alpha_df, hill_2, "fixslope"), 
       run_lmem(alpha_df, hill_2, "randomslope") )

# pielou_j ----------------------------------------------------------------------

alpha_df %>% 
  ggplot2::ggplot(aes(x = pielou_j)) +
  geom_histogram(binwidth = 0.05) +
  labs(title = "pielou j") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(pielou_j)), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)

alpha_df %>% 
  ggplot(aes(x = log(pielou_j))) +
  geom_histogram(binwidth = 0.1) +
  labs(title = "log(pielou j)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(log(pielou_j))), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)

alpha_df %>% 
  rownames_to_column("sample_id") %>% 
  ggplot(aes(x = "All Samples", y = pielou_j)) +
  geom_boxplot(outlier.color = "darkred") +
  geom_text_repel(aes(label = sample_id), color = "darkred") +
  labs(title = "Pielou's Evenness (pielou_j)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


# NOTE. Pileau J does not follow a perfect normal distribution, in particular
# there is/are outliers.

# Option A) Use lmer but check residuals ---------------------------------------
# Try the Same Baseline Model 
fit_pj_base<- run_lmem(alpha_df, pielou_j, "baseline")
summary(fit_pj_base)

# ***** MODEL DIAGNOSTICS ****** -----------------------------------------------

# Custom diagnostics -----------------------------------------------------------
diagnostic_plots(run_lmem(alpha_df, hill_0, "baseline"))
diagnostic_plots(run_lmem(alpha_df, hill_0, "fixslope"))

# DHARMA diagnostics -----------------------------------------------------------
# Generate simulated residuals (n=1000 is standard) and plot.
diagnostics_dharma(
  model = fit_pj_base,
  group_var1 = alpha_df$site,
  group_var2 = alpha_df$plot_rep
)

# INTERPRETASTION 
# The KS test is significant. The residual distribution is significantly different 
# from what the fitted Gaussian model expects. In other words, the residuals are 
# not normally distributed, and the model is not a good fit for the data.
# The dispersion test is n.s. so No over- or under-dispersion is detected.
# No strong outlier problem.
# Within-group Deviations detected, that means residual distribution differs 
# systematically across sites and plots. This is consistent with: Bounded data,
# Skew, Non-Gaussian residual shape.

# Additional option: Logit Transform + LMM 
alpha_df$logit_pj <- qlogis(alpha_df$pielou_j)

fit_logit_pj_fixslipe <- lmer(
  logit_pj ~ fert_status + (1 | site/plot_rep),
  data = alpha_df
)

summary(fit_logit_pj_fixslipe)

diagnostics_dharma(
  model     = fit_logit_pj_fixslipe,
  group_var1 = alpha_df$site,
  group_var2 = alpha_df$plot_rep
)

# option B) Beta Mixed Model (BEST OPTION) -------------------------------------
# Because Pielou J ∈ (0,1), the statistically principled model is

pielou_j_betaMM_base <- glmmTMB(
  pielou_j ~1 + (1 | site/plot_rep),
  family = beta_family(),
  data = alpha_df
)

summary(pielou_j_betaMM_base)

pielou_j_betaMM_fixslipe <- glmmTMB(
  pielou_j ~ fert_status + (1 | site/plot_rep),
  family = beta_family(),
  data = alpha_df
)

summary(pielou_j_betaMM_fixslipe)

anova(pielou_j_betaMM_base, pielou_j_betaMM_fixslipe) 

# option C) robust mixed effect models -----------------------------------------

# A robust mixed model downweights extreme observations so that outliers have less
# influence on parameter estimates, instead of deleting them. robustlmm::rlmer) 
# replace the normal likelihood with an M-estimation framework. Instead of minimizing 
# ∑ri2 they minimize ∑ρ(ri).
# Practically: Large residuals are downweighted, extreme observations influence 
# the model less, fixed effects become more stable. rlmer keep the outlier(s),
# estimates a weight for it, reduces its leverage automatically. Use robust LMM when
# most data are roughly Gaussian, on only a few extreme points exist, you believe 
# those points are real but atypical. 

# NOTE. Robust mixed models are not the appropriate primary solution for Pielou’s
# J simply because it is bounded between 0 and 1. They do not enforce predictions 
# within (0,1), model mean–variance relationships typical of proportions, correct
# skewness caused by boundary constraints, and ddress heteroscedasticity inherent 
# to bounded data.

summary(run_lmem_robust(alpha_df, pielou_j, "baseline"))
summary(run_lmem_robust(alpha_df, pielou_j, "fixslope"))
summary(run_lmem_robust(alpha_df, pielou_j, "randomslope"))


# **********************************************************************--------
# ***** VISUALIZING MIX MODELS ***** ------------------------------------------

# 1. Basic Boxplot by Fertilizer Status ----------------------------------------
ggplot(alpha_df, 
       aes(x = fert_status, y = hill_0, fill = fert_status)) +
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
       aes(x = fert_status, y = hill_0, color = site)) +
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
       aes(x = fert_status, y = hill_0, fill = fert_status)) +
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
ggeffects::ggpredict(run_lmem(alpha_df, hill_0, "baseline"), 
                     terms = "site", type = "random") %>% plot()

ggeffects::ggpredict(run_lmem(alpha_df, hill_0, "fixslope"), 
            terms = c("fert_status", "site")) %>% plot() 


# Plot-level predicted means
alpha_df$site_plot <- 
  interaction(alpha_df$site, alpha_df$plot_rep, drop = TRUE)

fit_hill_0_base <- lmer(
  hill_0 ~ 1 + (1 | site) + (1 | site_plot),
  data = alpha_df,
  REML = FALSE
)

pred_hill_0_base <-
  ggpredict(fit_hill_0_base, terms = "site_plot", type = "random")

plot(pred_hill_0_base) +
  labs(title = "Model Predictions: Hill Number by Fertilizer Status",
       subtitle = "Lines = site-specific predictions (random effects)",
       x = "Fertilizer Status",
       y = "Predicted Hill Number") 

# 5. Show Fixed + Random Effects (merTools) ------------------------------------

# Extract fixed and random effects
cbind(
  alpha_df,
  run_lmem(alpha_df, hill_0, "baseline") %>% 
  merTools::predictInterval(newdata = alpha_df, 
                            which = "full", 
                            level = 0.95)) %>% 
  rename(fitted = fit) %>%
  ggplot(aes(x = fert_status, y = hill_0)) +
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
  run_lmem(alpha_df, hill_0, "baseline"), effects = "ran_vals") %>% 
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
  run_lmem(alpha_df, hill_0, "fixslope"), effects = "ran_vals") %>% 
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

# 9. Full Diagnostic Plot Set (sjPlot) -----------------------------------------

sjPlot::plot_model(run_lmem(alpha_df, hill_0, "baseline"), type = "re")
sjPlot::plot_model(run_lmem(alpha_df, hill_0, "fixslope"), type = "re")

# Plot fixed effects with confidence intervals (can't do on the baseline as 
# there is no fixed effects)
plot_model(run_lmem(alpha_df, hill_0, "fixslope"), type = "est", 
           show.values = TRUE, 
           value.offset = 0.3) +
  labs(title = "Fixed Effect Estimates",
       subtitle = "Effect of fertilizer status on Hill number") +
  theme_minimal()

# Comprehensive model diagnostics and predictions
tab_model(run_lmem(alpha_df, hill_0, "baseline"), show.stat = TRUE)
tab_model(run_lmem(alpha_df, hill_0, "fixslope"), show.stat = TRUE) 


# 10. Raw Data + Model Predictions (emmeand) -----------------------------------

# Panel A: Raw data
ggplot(alpha_df, aes(x = fert_status, y = hill_0, fill = fert_status)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(title = "A. Raw Data",
       x = "Fertilizer Status", y = "Hill Number") +
  theme_minimal() +
  theme(legend.position = "none")

# Panel B: Model predictions
emmeans::emmeans(run_lmem(alpha_df, hill_0, "fixslope"), ~ fert_status) %>% 
  as.data.frame() %>% 
  ggplot(aes(x = fert_status, y = emmean, fill = fert_status)) +
  geom_col(alpha = 0.7) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  labs(title = "B. Model-Predicted Means",
       subtitle = "Error bars = 95% CI",
       x = "Fertilizer Status", y = "Predicted Hill Number") +
  theme_minimal() +
  theme(legend.position = "none")

# **********************************************************************--------
# Selecting BEST VISUALIZATION APPROACH ----------------------------------------

palette_fert <- c("#CC2D35", "#2D3142")

palette_site <- c("#009E73", "#0072B2", "#825121", "#E69F00", "#CC79A7")

# Hill 0 and Hill 2 ------------------------------------------------------------

# INTERPRETATION. Hill numbers are not independent indices. They are the same diversity 
# continuum at different q values. So including all three is often redundant unless you want
# to demonstrate consistency across evenness weighting.
# Hill_0 (q = 0) Pure richness, sensitive to rare taxa, influenced by sampling noise
# Hill_2 (q = 2) Strongly dominated by abundant taxa, less sensitive to rare taxa, 
# more stable across samples.
# Hill_0 + Hill_2 They bracket the diversity spectrum "Fertilization did not affect 
# richness (q=0) nor diversity weighted toward dominant taxa (q=2)."

# If Hill_0 = 150 → there are 150 observed taxa
# If Hill_2 = 24, it means: The community has the same dominance structure as a 
# community with 24 equally abundant species.

Figure_1_alpha <-
  ggarrange(
    ggarrange(
# Hill_0: Raw data with fixed effect means  
ggplot(alpha_df, 
       aes(x = fert_status, y = hill_0, color = site)) +
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

# Hill_2: Raw data with fixed effect means  
ggplot(alpha_df, 
       aes(x = fert_status, y = hill_2, color = site)) +
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
# Hill_0: Random effect deviations
broom.mixed::tidy(
  run_lmem(alpha_df, hill_0, "baseline"), effects = "ran_vals") %>% 
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
       subtitle = "deviations by site & plot",
       x = "Deviation from Fixed Effect",
       y = "Site:plot",
       color = NULL) +
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
  guides(color = guide_legend(ncol = 1)),

# Hill_2: Random effect deviations
broom.mixed::tidy(
  run_lmem(alpha_df, hill_2, "baseline"), effects = "ran_vals") %>% 
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
       subtitle = "deviations by site & plot",
       x = "Deviation from Fixed Effect",
       y = "Site:plot",
       color = NULL) +
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
  guides(color = guide_legend(ncol = 1)),

ncol = 2, nrow = 1,
labels = c("B","D"),
legend = "right",
common.legend = TRUE),

ncol = 1, 
nrow = 2
)

Figure_1_alpha

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

fit_hill0_site <- lmer(hill_0 ~ site + (1 | site:plot_rep),
               data = alpha_df, REML = FALSE)

anova(fit_hill0_site)   # overall test: does richness differ among sites?
summary(fit_hill0_site)

# Pairwise comparisons: “which sites differ?”
fit_hill0_site_emm <- emmeans(fit_hill0_site, ~ site)

fit_hill0_site_emm                       # site means (model-adjusted)
pairs(fit_hill0_site_emm, adjust = "tukey")  # Tukey-corrected pairwise tests
confint(fit_hill0_site_emm)              # 95% CIs for each site mean

as.data.frame(fit_hill0_site_emm) |>
  dplyr::arrange(dplyr::desc(emmean))

# if you also want to control for fertilizer while comparing sites

m_site_adj <- lmer(hill_0 ~ site + fert_status + (1 | site:plot_rep),
                   data = alpha_df, REML = FALSE)

emm_site_adj <- emmeans(m_site_adj, ~ site)
pairs(emm_site_adj, adjust = "tukey")

# This answers: “Which sites are richer after adjusting for fertilizer?”

# NOTE. The (1 | site:plot_rep) still belong to the model even when testing for
# differences between sites, because plots are still your experimental units, 
# and you have 3 pseudoreplicates per plot. Even when site is fixed, the 
# hierarchical structure does not disappear. 












