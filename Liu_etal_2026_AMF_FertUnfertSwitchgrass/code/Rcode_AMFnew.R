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
  renv,
  pak,
  styler,
  janitor,
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
  multcompView,
  install=FALSE
)


# Tracking package versions with renv ------------------------------------------
renv::init()      # initializes renv in your project
renv::restore()   # installs all packages from the lockfile
renv::snapshot()  # updates the lockfile

# Then commit and push the updated renv.lock file!

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

# save
writeXStringSet(
  otu_99_filtered,
  filepath = file.path(data_path, "datasets/otu_99_filtered.fasta"),
  format = "fasta"
)

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
  ) %>%
  column_to_rownames("SampleID") %>%
  janitor::clean_names()

head(metadata_99)

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

# NOTE. These looks like all real taxa to me. I am not going to remove any!

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

# **********************************************************************--------
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


# **********************************************************************--------
# **** 1. ALPHA DIVERSITY **** -------------------------------------------------
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
       run_lmem(alpha_df, hill_0, "fixslope") )

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
       run_lmem(alpha_df, hill_1, "fixslope") )

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
       run_lmem(alpha_df, hill_2, "fixslope") )

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

# Selecting BEST VISUALIZATION APPROACH ----------------------------------------

palette_fert <- c("#CC2D35", "#2D3142")

palette_site <- c("#009E73", "#0072B2", "#825121", "#E69F00", "#CC79A7")

# Hill 0 and Hill 2 ------------------------------------------------------------

# INTERPRETATION. Hill numbers are not independent indices. They are the same diversity 
# continuum at different q values. So including all three is often redundant unless 
# you want to demonstrate consistency across evenness weighting.
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

# ***** FIGURE 1 - Alpha-diversity ***** ----------------------------------------
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

# **********************************************************************--------
# **** 2. BETA DIVERSITY **** --------------------------------------------------

# Extracting metadata and otutable ---------------------------------------------
otutable_rare <- 
  physeq_AMF_rare %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>% 
  otu_table() %>%  
  as.matrix() %>% 
  t() %>% 
  as.data.frame()

otutable_rare

meta_rare <- 
  physeq_AMF_rare %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>% 
  sample_data() %>%   # Extract metadata table 
  as.matrix() %>% 
  as.data.frame() %>% 
  mutate(fert_status = factor(recode(
    fert_status, "FERT" = "Fertilized" ,  "UNFERT" ="Control"
  )))

meta_rare$site_plot <- 
  as.factor(
    paste(meta_rare$plot_rep, meta_rare$site, sep=":"))
  
meta_rare$site_plot <- with(
  meta_rare,
  factor(site_plot,
         levels = unique(site_plot[order(site, plot_rep)])))

head(meta_rare)


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
  otutable_rare ~ fert_status,
  data = meta_rare,
  method = "bray",
  permutations = how(blocks = meta_rare$site, nperm = 999),
  by = "margin"
)


# NOTE that using permutations = how(blocks = meta_rare$site) is the same as 
# adding the strata = site was used to restrict permutations within groups i.e.
# treatment labels are shuffled only within each site.

# 2) Does site impact community composition? -----------------------------------

# This model tests site effect, fertilization effect within site, and keeps 
# permutations restricted within site. Gives marginal tests (Type III-like logic)

adonis2(
  otutable_rare ~ site,
  data = meta_rare,
  method = "bray",
  permutations = 999,
  by = "term"
)

adonis2(
  otutable_rare ~ fert_status + site,
  data = meta_rare,
  method = "bray",
  permutations = 999,
  by = "term"
)

adonis2(
  otutable_rare ~ fert_status + site,
  data = meta_rare,
  method = "bray",
  permutations = how(blocks = meta_rare$site, nperm = 999),
  by = "margin"
)

adonis2(
  otutable_rare ~ fert_status + site,
  data = meta_rare,
  method = "bray",
  permutations = how(blocks = meta_rare$site, nperm = 999),
  by = "terms"
)

# 3) Does fertilizer effect differ among sites? --------------------------------
# NOTE. Should we check for an interaction between site and fert_status? Usually, 
# interactions are justified when you hypothesized context-dependent effects. For
# example that fert_status impact site A but not site B. If you have no reason to
# expect that, then the interaction is not justified.

adonis2(
  otutable_rare ~ fert_status * site,
  data = meta_rare,
  method = "bray",
  permutations = 999,
  by = "margin"
)

adonis2(
  otutable_rare ~ fert_status * site,
  data = meta_rare,
  method = "bray",
  permutations = how(blocks = meta_rare$site, nperm = 999),
  by = "margin"
)

adonis2(
  otutable_rare ~ fert_status * site,
  data = meta_rare,
  method = "bray",
  permutations = how(blocks = meta_rare$site, nperm = 999),
  by = "terms"
)

# 4) Does plot-level (beyond site) impact beta diversity? ----------------------
adonis2(
  otutable_rare ~ site_plot,
  data = meta_rare,
  method = "bray",
  permutations = how(blocks = meta_rare$site, nperm = 999),
  by = "margin"
)

adonis2(
  otutable_rare ~ site + site_plot,
  data = meta_rare,
  method = "bray",
  permutations = how(blocks = meta_rare$site, nperm = 999),
  by = "terms"
)

# 5) Does fertilizer effect vary among plots within sites? ---------------------
adonis2(
  otutable_rare ~ site_plot * fert_status,
  data = meta_rare,
  method = "bray",
  permutations = how(blocks = meta_rare$site, nperm = 999),
  by = "margin"
)

# INTERPRETATION. This appears significant, we only have 2 samples per plot 
# (treatment/control). In this situation the interaction is essentially capturing
# plot-specific treatment differences. Given n = 2 per plot, this term absorbs 
# nearly all treatment contrast variability (No independent estimate of treatment 
# variance). So it is mathematically allowed, but biologically unstable. 
# I would not emphasize this interaction given the paired 
# design and limited replication. The interaction term is essentially capturing:
# Differences in pairwise Bray–Curtis distances among plots.
# Not a replicated random slope.
# Not a properly estimated treatment heterogeneity.

# BETADISPER -------------------------------------------------------------------
# Check homogeneity of variances -----------------------------------------------
bray_dist <- vegdist(otutable_rare, method = "bray")

# The logic is identical to PERMANOVA — the permutation scheme must match the 
# scale at which the grouping variable exists:
# _Between sites → free permutations
# _Within sites (fert_status, plot) → restricted within sites

bd_fert <- betadisper(bray_dist, meta_rare$fert_status)
permutest(bd_fert, permutations = how(blocks = meta_rare$site, nperm = 999))
TukeyHSD(bd_fert)

bd_site <- betadisper(bray_dist, meta_rare$site)
permutest(bd_site, permutations = 999)
TukeyHSD(bd_site)

bd_plot <- betadisper(bray_dist, meta_rare$site_plot)
permutest(bd_plot, permutations = how(blocks = meta_rare$site, nperm = 999))
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

amf_nmds <- metaMDS(otutable_rare, distance = "bray", k = 2, trymax = 100)
amf_pcoa <- cmdscale(bray_dist, k = 2, eig = TRUE)

plot_ordination(
  ord = amf_nmds,
  meta = meta_rare,
  col_var = "site",
  shape_var = "fert_status"
)

plot_ordination(
  ord = amf_pcoa,
  meta = meta_rare,
  col_var = "site",
  shape_var = "fert_status", 
  ellipse = FALSE
)


# dispersion -------------------------------------------------------------------
plot_betadisper(dist_matrix = bray_dist, 
                grouping = meta_rare$fert_status)

plot_betadisper(dist_matrix = bray_dist, 
                grouping = meta_rare$site)

plot_betadisper(dist_matrix = bray_dist, 
                grouping = meta_rare$site_plot)


# ***** FIGURE 2 - Beta-diversity ***** -----------------------------------------

# using ggarrange --------------------------------------------------------------

Figure_2_beta_long <-
ggarrange(
  plot_ordination(
    ord = amf_pcoa,
    meta = meta_rare,
    col_var = "site",
    shape_var = "fert_status",
    ellipse = FALSE,
    legend_inside = TRUE
  ) +
    labs(title = "PCoA") +
    scale_color_manual(values = palette_site),
  plot_betadisper(
        dist_matrix = bray_dist,
        grouping = meta_rare$fert_status,
        signif_label = "Treatment, F=1.006, p=0.243"
        ) +
        labs(title = "Within-treatment variance",
             y = NULL),
  plot_betadisper(
    dist_matrix = bray_dist,
    grouping = meta_rare$site,
    signif_label = "Site, F=6.88, p=0.001"
  ) +
    labs(title = "Within-site variance",
         y = NULL),
  plot_betadisper(
    dist_matrix = bray_dist,
    grouping = meta_rare$site_plot,
    signif_label = "Plot, F=4.777, p=0.001"
  ) +
    labs(title = "Within-plot variance",
         y = NULL),
  ncol = 1,
  nrow = 4,
  align = "h",
  heights =  c(1.4, 0.6, 0.6, 1.2),
  labels = c("A", "B", "C", "D")
)

Figure_2_beta_long

ggsave(
  file.path(data_path, "results/Fig_2_betadiv.pdf"),
  plot = ggpubr::annotate_figure(
    Figure_2_beta_long,
    top = text_grob("BETA DIVERSITY", size = 12, face = "bold")
  ),
  device = "pdf"
)



ggarrange(
  plot_ordination(
    ord = amf_nmds,
    meta = meta_rare,
    col_var = "site",
    shape_var = "fert_status"
  ) +
    labs(title = "NMDS"),
  plot_betadisper(
    dist_matrix = bray_dist,
    grouping = meta_rare$site,
    signif_label = "Site, F=6.88, p=001"
  ) +
    labs(title = "Group dispersion"),
  ncol = 2,
  nrow = 1,
  align = "h",
  legend = "right",
  widths = c(1.4, 1),
  labels = c("A", "B")
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
                 "99" = physeq_99_rare,"100" = physeq_100_rare))

# Selecting Cetraspora gilmorei
assess_clustering_threshold(
  tax_rank = "Species",
  taxon_name =  "Cetraspora gilmorei",
  ps_list = list("90" = physeq_90_rare,"91" = physeq_91_rare,"92" = physeq_92_rare,
                 "93" = physeq_93_rare,"94" = physeq_94_rare,"95" = physeq_95_rare,
                 "96" = physeq_96_rare,"97" = physeq_97_rare,"98" = physeq_98_rare,
                 "99" = physeq_99_rare,"100" = physeq_100_rare))









# I think I should rarefy at the same depth

# plotting ordinations 

plot_pcoa <- function(physeq, 
                      Col=NULL, 
                      She=NULL, 
                      type = "samples"){
  
  taxa_num <- ncol(otu_table(physeq))
  print(taxa_num)
  
  pcoa <- phyloseq::ordinate(physeq, method = "PCoA", distance = "bray") 
  pcoa
  
  plot_res <-
    phyloseq::plot_ordination(physeq = physeq, 
                    ordination=pcoa,  
                    type=type, 
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





plot_pcoa(physeq = physeq_100_rare, Col = "site_id", She = "fert_status")
plot_pcoa(physeq = physeq_100_rare, Col = "site_id", She = "fert_status", type = "species")

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









