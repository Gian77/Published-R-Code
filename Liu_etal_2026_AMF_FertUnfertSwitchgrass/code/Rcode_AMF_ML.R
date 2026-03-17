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
  styler,
  magrittr,
  job,
  vegan,
  AICcPermanova,
  tidyverse,
  ggtext,
  ggpubr,
  gridExtra,
  ggdendro,
  ggraph,
  patchwork,
  GGally, # not needed 
  reshape2, # not needed 
  caret,
  mlr3,
  mlr3verse,
  mlr3filters,
  mlr3pipelines,
  mlr3learners,
  mlr3extralearners,
  mlr3tuning,
  paradox,
  mlr3misc,
  mlr3viz,
  iml,
  shapviz,
  yaImpute,
  precrec,
  data.table,
  pls,
  fastDummies,
  janitor,
  knitr,
  install = TRUE
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
# ***** MACHINE LEARNING ***** -------------------------------------------------

# Prepare the data -------------------------------------------------------------
model_data <- 
  otu_chem_data %>%
  dplyr::select(
    dry_matter_yield_mg_ha,
    site.x, fert_status,
    hill_0, hill_1, hill_2,
    p_ppm, k_ppm, ca_ppm, mg_ppm, cec, p_h,
    moisture_at_harvest
  ) %>%
dplyr::rename(site = site.x) %>% 
  mutate(
    fert_status = as.factor(fert_status),
    site        = as.factor(site)
  ) %>% 
mutate(across(.cols = c(1, 4:13), .fns = as.numeric))

head(model_data)
dim(model_data)


# ******************************************************************************

# BENCHMARK: 3 LEARNERS × 2 RESAMPLING STRATEGIES


# ── 1. TASK
model_data_encoded <- 
  model_data %>%
  mutate(fert_fertilized = as.integer(fert_status == "Fertilized")) %>% 
  select(-fert_status)

task <- as_task_regr(
  model_data_encoded,
  target = "dry_matter_yield_mg_ha",
  id     = "biomass"
)

# stratify to site 
task$set_col_roles("site", roles = "stratum")
task

# ── 2. LEARNERS

# Elastic Net (alpha = 0.5, L1 + L2)
learner_enet <- as_learner(
  po("encode", method = "one-hot") %>>%
    lrn("regr.cv_glmnet", 
        alpha = 0.5, 
        standardize = TRUE, 
        nfolds = 4,
        id = "enet")
)
learner_enet$id <- "Elastic_Net"

# Ridge (alpha = 0, L2 only)
learner_ridge <- as_learner(
  po("encode", method = "one-hot") %>>%
    lrn("regr.cv_glmnet", 
        alpha = 0, 
        standardize = TRUE, 
        nfolds = 4,
        id = "ridge")
)
learner_ridge$id <- "Ridge"

# Random Forest
learner_rf <- lrn("regr.ranger",
                  importance  = "permutation",
                  num.trees   = 500,
                  respect.unordered.factors = "order",
                  id = "RF"
)

learner_rf$id <- "Random_Forest"

learners <- list(learner_enet, learner_ridge, learner_rf)

# ── 3. RESAMPLING STRATEGIES 

# resampling 
resampling_blocked = rsmp("cv", folds = 5, repeats = 10)
resampling_blocked <- rsmp("repeated_cv", folds = 4, repeats = 20)

# ── 4. BENCHMARK GRID
design <- benchmark_grid(
  tasks      = task,
  learners   = learners,
  resamplings = list(resampling_blocked)
)

# ── 5. RUN BENCHMARK 
set.seed(31626)
bmr <- benchmark(design, store_models = TRUE)

# ── 6. RESULTS TABLE
measures <- list(
  msr("regr.rsq"),
  msr("regr.rmse"),
  msr("regr.mae")
)

results <- bmr$aggregate(measures) %>%
  as.data.frame() %>%
  select(learner_id, resampling_id, 
         regr.rsq, regr.rmse, regr.mae) %>%
  arrange(resampling_id, desc(regr.rsq))

results

# ── 7. VISUALIZATION 

# --- 7a. R² comparison across learners and resampling
bmr$aggregate(measures) %>%
  as.data.frame() %>%
  ggplot(aes(x = learner_id, y = regr.rsq, fill = resampling_id)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey40") +
  scale_fill_manual(values = c("Site_LOO" = "#185FA5", 
                               "Standard_LOO" = "#A32D2D")) +
  labs(
    title    = "Benchmark: R² by learner",
    x        = NULL, 
    y = "R²", 
    fill = "Resampling"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

# --- 7b. Per-fold performance (shows variance across sites)
bmr$score(measures) %>%
  as.data.frame() %>%
  ggplot(aes(x = learner_id, y = regr.rsq, col = learner_id)) +
  geom_jitter(width = 0.1, size = 3) +
  stat_summary(fun = mean, geom = "crossbar", 
               width = 0.4, col = "black") +
  labs(
    title    = "Per-fold R²",
    subtitle = "Each point = one fold | Bar = mean",
    x        = NULL, y = "R²"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

# --- 7c. RMSE heatmap
bmr$aggregate(measures) %>%
  as.data.frame() %>%
  ggplot(aes(x = learner_id, y = resampling_id, fill = regr.rmse)) +
  geom_tile(col = "white", linewidth = 1) +
  geom_text(aes(label = round(regr.rmse, 2)),
            size = 5, fontface = "bold",
            col = "white") +
  scale_fill_gradient(low = "#185FA5", high = "#A32D2D") +
  labs(
    title = "RMSE heatmap — lower is better",
    x = NULL, y = NULL, fill = "RMSE"
  ) +
  theme_bw(base_size = 13)

# ── 8. EXTRACT BEST LEARNER 
best_learner_id <- 
  results %>%
  slice_max(regr.rsq, n = 1) %>%
  pull(learner_id)

best_learner_id

# ── 9. VARIABLE IMPORTANCE FROM BEST LEARNER 

# Retrain best learner on full data
best_learner <- learners[[which(sapply(learners, `[[`, "id") == best_learner_id)]]
best_learner$train(task)

rr = resample(task, best_learner, resampling_blocked, store_models = TRUE)
print(rr$aggregate(list(msr("regr.rsq"), msr("regr.rmse"), msr("regr.mae"))))

# Extract importance depending on learner type
if (best_learner_id == "Random_Forest") {
  imp_df <- data.frame(
    Variable   = names(best_learner$importance()),
    Importance = as.numeric(best_learner$importance())
  ) %>%
    filter(Importance > 0) %>%
    arrange(desc(Importance)) %>%
    mutate(Variable = factor(Variable, levels = Variable))
  
  ggplot(imp_df, aes(x = Importance, y = Variable)) +
    geom_col(fill = "#185FA5") +
    labs(title = "RF permutation importance",
         x = "Permutation importance", y = NULL) +
    theme_bw(base_size = 13)
  
} else {
  # Elastic Net or Ridge — extract coefficients
  glmnet_model <- best_learner$model$enet$model  # adjust id if needed
  coef_mat <- coef(glmnet_model, s = "lambda.min")
  
  coef_df <- data.frame(
    Variable    = rownames(coef_mat)[-1],
    Coefficient = as.numeric(coef_mat)[-1]
  ) %>%
    filter(Coefficient != 0) %>%
    arrange(desc(abs(Coefficient))) %>%
    mutate(
      Direction = ifelse(Coefficient > 0, "Positive", "Negative"),
      Variable  = factor(Variable, levels = Variable)
    )
  
  ggplot(coef_df, aes(x = Coefficient, y = Variable, fill = Direction)) +
    geom_col() +
    scale_fill_manual(values = c("Positive" = "#185FA5", 
                                 "Negative" = "#A32D2D")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(title = paste(best_learner_id, "— surviving predictors"),
         x = "Standardized coefficient", y = NULL) +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom")
}


# ── 10 NULL MODEL
run_permuted_iteration <- function(task, learner, resampling) {
  
  # 1. Gather all required column names
  # We need the features, the target, AND the grouping column (site)
  all_cols <- c(task$feature_names, task$target_names, "site")
  
  # Get the data for just those columns
  perm_data <- task$data(cols = all_cols)
  
  # Shuffle the target
  perm_data[[task$target_names]] <- sample(perm_data[[task$target_names]])
  
  # 2. Create the new task
  task_perm <- TaskRegr$new(
    id = "perm_task",
    backend = perm_data,
    target = task$target_names
  )
  
  # 3. Re-assign the site grouping
  task_perm$set_col_roles("site", roles = "group")
  
  # 4. Run Resampling (cloning to ensure a clean slate)
  rr_perm <- resample(task_perm, learner$clone(), resampling)
  
  return(rr_perm$aggregate(msr("regr.rsq"))[[1]])
}

# Run 100 permutations (This may take a minute)
null_rsq_distribution <-
  replicate(10,
            run_permuted_iteration(task, best_learner, resampling_blocked))

data.frame(rsq = null_rsq_distribution) %>% 
  ggplot(aes(x = rsq)) +
  geom_histogram(bins = 30, fill = "gray", color = "white") +
  geom_vline(xintercept = -0.32, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Permutation Test (Null Model)",
       subtitle = "Red line is your actual Model R² (-0.32)",
       x = "Cross-Validated R²",
       y = "Frequency") +
  theme_classic()



