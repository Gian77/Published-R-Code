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
  mlr3,              # Core package for machine learning in R using the mlr3 framework
  mlr3verse,         # Meta-package loading common mlr3 extensions for a full ML workflow
  mlr3filters,       # Feature selection and filter methods for mlr3
  mlr3pipelines,     # Pipeline operators for building complex preprocessing and modeling workflows
  mlr3learners,      # Collection of standard machine learning learners for mlr3
  mlr3extralearners, # Community-contributed set of additional learners for mlr3
  mlr3tuning,        # Hyperparameter tuning and optimization tools for mlr3
  paradox,           # Parameter space definition and sampling framework used by mlr3
  mlr3misc,          # Utility functions supporting other mlr3 packages
  mlr3viz,           # Visualization tools for mlr3 objects (learners, benchmarks, etc.)
  iml,               # Interpretable Machine Learning – model-agnostic interpretation tools
  shapviz,           # Fast SHAP value computation and visualization
  styler,            # Code formatter for enforcing tidyverse-style R code layout
  magrittr,          # Provides the pipe operator (%>%) used for functional chaining
  job,               # Lightweight parallel computing framework for R jobs
  AICcPermanova,     # PERMANOVA implementation with AICc model selection for ecological data
  tidyverse,         # Core collection of tidy data science packages (ggplot2, dplyr, tidyr, etc.)
  ggtext,            # Enhanced text rendering in ggplot2 (markdown, HTML formatting)
  ggpubr,            # Publication-ready ggplot2 visualizations with simple syntax helpers
  gridExtra,         # Functions to arrange multiple grid-based plots in a layout
  patchwork,         # Combine separate ggplots into a unified layout using + and /
  data.table,        # High-performance data manipulation and aggregation
  pls,               # Partial Least Squares and Principal Component Regression
  fastDummies,       # Fast creation of dummy (one-hot encoded) variables
  janitor,           # Data cleaning and exploration helper functions
  knitr,             # Dynamic report generation integrating R code and text
  install = FALSE
)


# GGally, # not needed 
# reshape2, # not needed 
# caret,
# yaImpute,
# precrec,

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

# (A) CHEMSITRY AND AMF METRICS ------------------------------------------------

# Prepare the data -------------------------------------------------------------
model_data <- 
  otu_chem_data %>%
  dplyr::select(
    dry_matter_yield_mg_ha,
    site.x, fert_status,
    hill_0, hill_1, hill_2,
    bray.1,bray.2,jacc.1,jacc.2,
    p_ppm, k_ppm, ca_ppm, mg_ppm, cec, p_h,
    moisture_at_harvest
  ) %>%
dplyr::rename(site = site.x) %>% 
  mutate(
    fert_status = as.factor(fert_status),
    site        = as.factor(site)
  ) %>% 
mutate(across(.cols = c(1, 4:17), .fns = as.numeric))

head(model_data)
dim(model_data)

# Check for 
plot(model_data$hill_2, model_data$hill_0)
plot(model_data$p_ppm, model_data$p_h)
plot(model_data$ca_ppm, model_data$mg_ppm)
plot(model_data$ca_ppm, model_data$cec)

# ******************************************************************************

# 1) recode fert_status to 0-1 and remove site ---------------------------------
model_data_encoded <- 
  model_data %>%
  mutate(Status_fertilized = as.integer(fert_status == "Fertilized")) %>% 
  dplyr::select(-fert_status, -moisture_at_harvest) %>% 
  dplyr::select(-hill_1, -jacc.1, -jacc.2, -bray.2, -mg_ppm, -cec)

head(model_data_encoded)
str(model_data_encoded)

# 2) Create the task, gorup by site, remove site -------------------------------
task <- as_task_regr(
  model_data_encoded,
  target = "dry_matter_yield_mg_ha",
  id     = "biomass"
)

task

# 3) Learners for benchmarking -------------------------------------------------

# WARNING. Predicotrs are on different scale so it is necessary to rescale.
learner_null <- lrn("regr.featureless")
learner_null$id <- "Featureless"

learner_enet <- as_learner(
  po("encode", method = "one-hot") %>>%
    po("scale") %>>%          # scaling it here 
    lrn("regr.cv_glmnet",
        alpha = 0.5,
        standardize = FALSE,  # IMPORTANT: avoid double scaling
        nfolds = 4)
)
learner_enet$id <- "Elastic_Net"

learner_ridge <- as_learner(
  po("encode", method = "one-hot") %>>%
    lrn("regr.cv_glmnet", 
        alpha = 0, 
        standardize = FALSE,  # IMPORTANT: avoid double scaling
        nfolds = 4)
)
learner_ridge$id <- "Ridge"

learner_lasso <- as_learner(
  po("encode", method = "one-hot") %>>%
    po("scale") %>>%
    lrn("regr.cv_glmnet",
        alpha = 1,
        standardize = FALSE, # IMPORTANT: avoid double scaling
        nfolds = 4)
)
learner_lasso$id <- "Lasso"

learner_kriging <- as_learner(
  po("encode", method = "one-hot") %>>%
    po("scale") %>>%
    lrn("regr.km")
)
learner_kriging$id <- "Kriging_GPR"

learner_tree <- as_learner(
  po("encode", method = "one-hot") %>>%
    po("scale") %>>%
    lrn("regr.rpart", cp = 0.01)
)
learner_tree$id <- "Decision_Tree"

learner_svm <- as_learner(
  po("encode", method = "one-hot") %>>%
    po("scale") %>>%
    lrn("regr.ksvm",
        type = "eps-svr",
        kernel = "rbfdot")
)
learner_svm$id <- "SVM"

learner_rf <- as_learner(
  po("encode", method = "one-hot") %>>%
    po("scale") %>>%
    lrn("regr.ranger",
        importance = "permutation",
        num.trees = 500,
        respect.unordered.factors = "order")
)
learner_rf$id <- "Random_Forest"

learner_xgb <- as_learner(
  po("encode", method = "one-hot") %>>%
    po("scale") %>>%
    lrn("regr.xgboost",
        nrounds = 100,
        eta = 0.05,
        max_depth = 2,
        lambda = 1,
        alpha = 1)
)
learner_xgb$id <- "XGBoost"

learners <- list(learner_null,
                 learner_enet, learner_ridge, learner_lasso,learner_kriging,
                 learner_tree, learner_svm, learner_rf, learner_xgb)

# Reduce the learners to a more simple set
learners <- list(learner_null, 
                 learner_ridge, learner_enet, learner_lasso, 
                 learner_svm, learner_rf)

learners <- list(learner_null, learner_enet, learner_tree, 
                 learner_svm, learner_rf)


# 4) Group resampling by site --------------------------------------------------

# Group by site and remove site
task$set_col_roles("site", roles = "group")

# Remove the site column
task$select(setdiff(task$feature_names, "site"))
task

# Resampling strategy with site as Random effect
resampling_outer <- rsmp("cv", folds = 5)
resampling_outer$instantiate(task)
resampling_outer

# NOTE. To avoidn model instability (see below) inclueding repeats may help,
# especially when dealing with small datasets. This below, runs the entire 
# benchmark 10 times and averages the results. It is much harder for a "lucky" 
# split to crown a fake winner if that model has to perform well across 50 
# different test sets (10 repeats × 5 folds).

resampling_outer <- rsmp("repeated_cv", folds = 5, repeats = 3)
resampling_outer$instantiate(task)
resampling_outer

# 5) Benckmark -----------------------------------------------------------------
set.seed(64387628)

bmr <- mlr3::benchmark(
  benchmark_grid(
    tasks = task,
    learners = learners,
    resamplings = resampling_outer),
  store_models = TRUE
)

bmr

autoplot(bmr, measure = msr("regr.rmse"), type = "boxplot")

# WARNING. Small data set can cause model instability, that is why if I ran this 
# multiple time I almost always have a new best_learner (See the NOTE above). 

# 6) Model performace summary --------------------------------------------------
results <- 
  bmr$aggregate(list(msr("regr.rsq"),msr("regr.rmse"),msr("regr.mae")))

results %>% arrange(desc(regr.rsq))

# Alternative way
results_df <- 
  bmr$score(list(msr("regr.rsq"),msr("regr.rmse"),msr("regr.mae"))) %>%
  as.data.frame()

results_df %>%
  group_by(learner_id) %>%
  summarise(
    rsq_mean  = mean(regr.rsq, na.rm = TRUE),
    rsq_sd    = sd(regr.rsq, na.rm = TRUE),
    rmse_mean = mean(regr.rmse, na.rm = TRUE),
    rmse_sd   = sd(regr.rmse, na.rm = TRUE),
    mae_mean  = mean(regr.mae, na.rm = TRUE),
    mae_sd    = sd(regr.mae, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(rmse_sd)


# 7) Select best learner and etrain best learner on full data

# Selection based on RMSE
best_learner_id <- 
  results[order(regr.rmse)][1]$learner_id 
best_learner_id

# Selection based on R2, decreasing = TRUE becasue is negative 
best_learner_id <- 
  results[order(regr.rsq, decreasing = TRUE)][1]$learner_id 
best_learner_id

# or 
best_learner_id <- "Random_Forest"

best_learner <- 
  learners[[which(sapply(learners, function(x) x$id) == best_learner_id)]]$clone()
best_learner

# 7) Run a null model to verify validity of the model --------------------------

# INTERPRETATION. For a Null Model, global shuffling is the standard approach.
# Basically I am trying to answer this question: "If yield was totally random and
#    had nothing to do with soil chemistry OR site location, what would 
#    the R2 look like?"

# In the function while the resampling (the rsmp("cv") part) correctly keeps sites 
# grouped during the training/testing phase, the shuffling step happens globally 
# across the entire dataset before the sites are even defined in the new task.
# Because sample() is called on the entire column of dry_matter_yield_mg_ha, a 
# yield value that originally belonged to Hancock could be reassigned to a 
# row in Escanaba.

# In this situation, the null model is "dumber." It doesn't even know that some 
# sites naturally have higher yields than others.

global_permuted_iteration <- 
  function(task, 
           learner,
           resampling) {
  
  # We need the features, the target, AND the grouping column (site)
  all_cols <- c(task$feature_names, task$target_names, "site")
  
  # Get the data for just those columns
  perm_data <- task$data(cols = all_cols)
  
  # Shuffle the target
  perm_data[[task$target_names]] <- sample(perm_data[[task$target_names]])
  
  # 2. Create the new task
  task_perm <- mlr3::TaskRegr$new(
    id = "perm_task",
    backend = perm_data,
    target = task$target_names
  )
  
  # 3. Re-assign the site grouping
  task_perm$set_col_roles("site", roles = "group")
  
  # 4. Run Resampling (cloning to ensure a clean slate)
  rr_perm <- mlr3::resample(task_perm, learner$clone(), resampling)
  
  return(rr_perm$aggregate(msr("regr.rsq"))[[1]])
}

# Run 100 permutations (This may take a minute)
null_rsq_distribution <-
  replicate(
    n = 100,
    expr = global_permuted_iteration(
      task = task,
      learner =  best_learner,
      resampling =  resampling_outer
    )
  )

head(null_rsq_distribution)

# NOTE. A "Constrained" Null Model (to see if your features matter beyond just 
# knowing which site you are in), you should shuffle within each site. This keeps
# the "Site Mean" intact but breaks the link between the microbes/chemistry 
# and the yield. In this situation the null model "knows" the site means but 
# doesn't know the chemistry. This null distribution will likely shift closer 
# to your actual model's R2..

constrained_permuted_iteration <- 
  function(task, 
           learner, 
           group_var = "site",
           resampling) {

    all_cols <- c(task$feature_names, task$target_names, group_var)
    perm_data <- 
      task$data(cols = all_cols) %>%
      group_by(!!sym(group_var)) %>% 
      mutate(!!sym(task$target_names) := sample(!!sym(task$target_names))) %>%
      ungroup()
  
  task_perm <- mlr3::TaskRegr$new(id = "perm_task", 
                            backend = perm_data, 
                            target = task$target_names)
  task_perm$set_col_roles(group_var, roles = "group")
  
  rr_perm <- mlr3::resample(task_perm, learner$clone(), resampling)
  return(rr_perm$aggregate(msr("regr.rsq"))[[1]])
}


# Run 100 permutations (This may take a minute)
null_rsq_distribution_constr <-
  replicate(
    n=100,
    expr = constrained_permuted_iteration(
      task = task,
      learner =  best_learner,
      group_var = "site",
      resampling =  resampling_outer
    )
  )

head(null_rsq_distribution_constr)

# **** FIGURE XX - Null Models **** --------------------------------------------

Fig_XX_null_models <- 
  ggarrange(
  data.frame(rsq = null_rsq_distribution) %>% 
  ggplot(aes(x = rsq)) +
  geom_histogram(bins = 30, fill = "gray", color = "white") +
  # The Red Line
  geom_vline(xintercept = results[order(regr.rmse)][1]$regr.rsq, 
             color = "darkred", 
             linetype = "dashed",
             linewidth = 1) +
  # ADDING THE LABEL ON THE LINE
  annotate("text", x = -Inf, y = Inf,
           label = paste("R² =", round(results[order(regr.rmse)][1]$regr.rsq, 3)),
           color = "darkred", hjust = -0.1, vjust = 1.1, fontface = "bold") +
  labs(title = "Null Model (Global shuffle)",
       # Using paste() here makes the variable render correctly
       #subtitle = paste("Actual Model R² =", round(results[order(regr.rmse)][1]$regr.rsq, 2)),
       x = "Cross-Validated R²",
       y = "Frequency") +
  theme_classic()+ 
    theme(
      plot.title = element_markdown(size =12, face = "bold",hjust = 0.5, vjust = 0.5),
      plot.subtitle = element_markdown(size = 10, hjust = 0.5, vjust = 0.5),
      axis.text.x = element_markdown(size = 8),
      axis.text.y = element_markdown(size = 8),
      legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"),
      legend.title = element_blank(), legend.text = element_text(size = 8)) ,
  #ylim(0, 20) ,
data.frame(rsq = null_rsq_distribution_constr) %>% 
  ggplot(aes(x = rsq)) +
  geom_histogram(bins = 30, fill = "gray", color = "white") +
  # The Red Line
  geom_vline(xintercept = results[order(regr.rmse)][1]$regr.rsq, 
             color = "darkred", 
             linetype = "dashed",
             linewidth = 1) +
  # ADDING THE LABEL ON THE LINE
  #annotate("text", x = -Inf, y = Inf,
  #         label = paste("Model R² =", round(results[order(regr.rmse)][1]$regr.rsq, 3)),
  #         color = "darkred", hjust = -0.1, vjust = 1.1, fontface = "bold") +
  labs(title = "Null Model (Within-site shuffle)",
       # Using paste() here makes the variable render correctly
       #subtitle = paste("Actual Model R² =", round(results[order(regr.rmse)][1]$regr.rsq, 2)),
       x = "Cross-Validated R²",
       y = "Frequency") +
  theme_classic()+ 
  theme(
    plot.title = element_markdown(size =12, face = "bold",hjust = 0.5, vjust = 0.5),
    plot.subtitle = element_markdown(size = 10, hjust = 0.5, vjust = 0.5),
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 8),
    legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"),
    legend.title = element_blank(), legend.text = element_text(size = 8))  ,
  #ylim(0, 30),
labels = c("A", "B"),
ncol = 2,
nrow = 1)

Fig_XX_null_models


# 8) Re-running with the best learner and auto-tuning -----------------------------

# Define tuning space for SVM
search_space_svm <- ps(
  regr.ksvm.C = p_dbl(lower = 0.1, upper = 10, logscale = TRUE),
  regr.ksvm.sigma = p_dbl(lower = 0.01, upper = 10, logscale = TRUE)
)

search_space_rf <- ps(
  regr.ranger.mtry = p_int(lower = 2,upper = length(task$feature_names)),
  regr.ranger.min.node.size = p_int(lower = 2,upper = 10),
  regr.ranger.sample.fraction = p_dbl(lower = 0.6, upper = 1)
)

# IMPORTANT. For inner resampling within tuning (that must also be grouped)
# I use 4 folds because in each outer fold, the training set only has 4 sites left,
# the other one is kept for testing, you do not want the model to see it and 
# learn from it.
resampling_inner <- rsmp("cv", folds = 4)

# Or to stabilize th model try 
resampling_inner <- rsmp("repeated_cv", folds = 4, repeats = 10)

# AutoTuner within in the learner
at_svm <- auto_tuner(
  learner = best_learner,
  resampling = resampling_inner,
  measure = msr("regr.rmse"),
  search_space = search_space_svm,
  terminator = trm("evals", n_evals = 30),
  tuner = mlr3verse::tnr("random_search")
)

at_svm

at_rf <- auto_tuner(
  learner = best_learner,
  resampling = resampling_inner,
  measure = msr("regr.rmse"),
  search_space = search_space_rf,
  terminator = trm("evals", n_evals = 40),  # slightly higher than SVM
  tuner = mlr3verse::tnr("random_search")
)

at_rf

# Resample with tuning (nested cv) for st_rf
set.seed(657564)

rr_tuned <- 
  mlr3::resample(
  learner =   at_rf,
  task = task,
  resampling = resampling_outer,
  store_models = TRUE
)

rr_tuned$aggregate(list(msr("regr.rsq"),msr("regr.rmse"),msr("regr.mae")))

autoplot(rr_tuned, measure = msr("regr.rmse"), type = "boxplot")

# INTERPRETATION. The R2 drop from the benchmark RF to the tuned RF, why?
# Well, this is a very well-known phenomenon called Tuning Bias (or "overfitting
# the validation set"). This is a classic "Small Data" trap. With 38 samples when
# you tune the model, the AutoTuner is looking at only 4 sites in the inner 
# loop (about 30 samples). This is what's happening:
# _ Untuned RF: Uses the default mtry (usually $\sqrt{p}$) and nodesize. These 
#   defaults are designed to be "safe" and conservative.
# _ Tuned RF: The tuner tries 
#   40 different combinations of settings. By pure chance, it might find a specific
#   mtry and nodesize that happens to work perfectly for the 4 training sites. 
#   It thinks, "Aha! I found the secret!" but it actually just found a lucky 
#   coincidence. When it tries to apply that "lucky" setting to the 5th site 
#   (the outer loop), it fails miserably because that coincidence doesn't 
#   exist there.
# Also, the Inner Resampling (where the tuning happens) is telling the model what 
# works within the training sites. But your Outer Resampling is asking if it 
# works on a new site. The model is essentially "over-optimizing" for the training
# environment.

# Increasing combinations to mre thna 40 will  likely just increasing the risk 
# of finding a "False Prophet". So better to stick with the benchmark here and say:
# Due to the limited sample size (n=38) and high inter-site variability, 
# hyperparameter tuning led to overfitting of the training sites. Therefore, 
# the more conservative default parameters were utilized to ensure model stability
# across locations."

# Inspect chosen hyperparameters
rr_tuned$learners[[1]]$model$tuning_instance

# Or across folds
lapply(rr_tuned$learners, function(x) x$model$tuning_instance)

# 9) Fold wise permuation importance (no leakage) ------------------------------

get_fold_importance <- function(rr, task) {
  
  imps <- lapply(seq_len(rr$iters), function(i) {
    
    model <- rr$learners[[i]]
    
    test_idx <- rr$resampling$test_set(i)
    X_test <- as.data.frame(task$data(rows = test_idx, cols = task$feature_names))
    y_test <- task$data(rows = test_idx, cols = task$target_names)[[1]]
    
    pred_fun <- function(model, newdata) {
      newdata[[task$target_names]] <- NA_real_
      new_task <- TaskRegr$new(
        id = "newdata",
        backend = newdata,
        target = task$target_names
      )
      model$predict(new_task)$response
    }
    
    predictor <- Predictor$new(
      model = model,
      data  = X_test,
      y     = y_test,
      predict.function = pred_fun
    )
    
    #imp <- FeatureImp$new(predictor, loss = "rmse")
    # This is calculates as:
    # importance = permuted_error - original_error
    imp <- FeatureImp$new(
      predictor,
      loss = "rmse",
      compare = "difference"
    )
    
    dt <- as.data.table(imp$results)
    dt[, fold := i]
    
    dt
  })
  
  rbindlist(imps)
}

# 1) Using the tuned model --------------------------------------------------------
imp_df <-
  get_fold_importance(rr = rr_tuned, task = task)

str(imp_df)

# Filter the folds based on model performance
rr_tuned$score(msr("regr.rsq"))$regr.rsq 

fold_sites <- lapply(seq_len(rr_tuned$iters), function(i) {
  test_idx <- rr_tuned$resampling$test_set(i)
  unique(task$data(rows = test_idx, cols = "site")[[1]])
})

fold_sites

# INTERPRETATION. Fold 2 is not perfoming as other folds. 

# Compare response distribution
task$data(cols = c("site", task$target_names)) %>%
  as.data.frame() %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(
    mean = mean(dry_matter_yield_mg_ha),
    sd   = sd(dry_matter_yield_mg_ha)
  )

# 2) Auto-tunign ineffective - Retaining benchmark model----------------------------

# Extracting the model form the benchmark instead of using the tuned version.

rr_svm_from_bmr <- bmr$resample_result(learner_id = "SVM")

lapply(seq_len(rr_svm_from_bmr$iters), function(i) {
  test_idx <- rr_svm_from_bmr$resampling$test_set(i)
  unique(task$data(rows = test_idx, cols = "site")[[1]])
})

imp_df_svm <-
  get_fold_importance(rr = rr_svm_from_bmr, task = task)

imp_df_svm

# Aggregate importance
imp_summary <-
  imp_df_svm %>%
  as.data.frame() %>%
  group_by(feature) %>% 
  dplyr::summarise(
    mean_importance = mean(importance),
    sd_importance   = sd(importance)
  ) %>% 
  arrange(mean_importance)

imp_summary

# Rank stability importance across folds
imp_df_svm[, rank := rank(importance), by = fold]

rank_summary <- imp_df_svm[, .(
  mean_rank = mean(rank),
  sd_rank   = sd(rank)
), by = feature][order(mean_rank)]

rank_summary

rank_summary[, importance_score := (max(mean_rank) + 1) - mean_rank]


# Random Forest
rr_rf_from_bmr <- bmr$resample_result(learner_id = "Random_Forest")

lapply(seq_len(rr_rf_from_bmr$iters), function(i) {
  test_idx <- rr_rf_from_bmr$resampling$test_set(i)
  unique(task$data(rows = test_idx, cols = "site")[[1]])
})

imp_df_rf <-
  get_fold_importance(rr = rr_rf_from_bmr, task = task)

str(imp_df_rf)

# Aggregate importance
# Aggregate importance
imp_summary <-
  imp_df_rf %>%
  as.data.frame() %>%
  group_by(feature) %>% 
  dplyr::summarise(
    mean_importance = mean(importance),
    sd_importance   = sd(importance)
  ) %>% 
  arrange(mean_importance)

imp_summary


imp_summary <- imp_df[, .(
  mean_imp = mean(importance),
  sd_imp   = sd(importance)
), by = feature][order(-mean_imp)]

imp_summary

# Rank stability importance across folds
imp_df[, rank := rank(-importance), by = fold]

rank_summary <- imp_df[, .(
  mean_rank = mean(rank),
  sd_rank   = sd(rank)
), by = feature][order(mean_rank)]

print(rank_summary)

rank_summary[, importance_score := (max(mean_rank) + 1) - mean_rank]

# **** FIGURE XXX - Ranked features per fold **** ------------------------------

Fig_XX_feature_rank <- 
  rank_summary %>% 
ggplot(aes(x = reorder(feature, -mean_rank), y = mean_rank)) +
  #geom_point(size = 4, shape = 18) +
  geom_col(fill = "grey70") +
  geom_errorbar(aes(ymin = mean_rank - sd_rank,
                    ymax = mean_rank + sd_rank),
                    width = 0) +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Predictor",
    y = "Mean Rank (lower = more important)",
    title = "Variable Importance",
    subtitle = "Random Forest (Grouped CV)") +
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
  )

Fig_XX_feature_rank

# INTERPRETATION. To assess the stability of feature importance across different 
# geographic contexts, we performed a Rank Stability Analysis within our 
# Leave-One-Site-Out (LOSO) cross-validation framework (1 site used for testing 
# and 4 for training). For each of the five folds (representing the five experimental
# sites), predictors were ranked from 1 (highest contribution to model accuracy) 
# to 14 (lowest contribution) based on their permutation importance (increase in 
# RMSE). We then calculated the mean rank and standard deviation of ranks for 
# each feature across all folds. This dual-metric approach allows for the 
# differentiation between 'Global Drivers' (low mean rank, low SD) and 
# weak contributors and 'Context-Dependent Drivers' (high mean rank, high SD).
# Remember that "Rank" and "Importance Value" move in opposite directions.

# **** FIGURE XXX - Importance per fold **** -----------------------------------

Fig_XX_imortance_per_fold <- 
  imp_df %>%
  group_by(fold) %>%
  slice_max(order_by = importance, n = 3) %>%
  ungroup() %>% 
  mutate(site = recode_values(fold,
                    1 ~ "Hancock",2 ~ "Escanaba",3 ~ "Lake City",
                    4 ~ "Lux Arbor",5 ~ "Rhinelander")) %>% 
  mutate(feature_key = factor(paste0(fold, "_", feature)),
         feature_key = gsub("\\.", "_", feature_key),
         feature_key = reorder(feature_key, importance)) %>%
  ggplot(aes(x = importance, 
             y = feature_key, 
             fill = as.factor(site))) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~site, scales = "free_y") +
  # 4. Use gsub to strip the "1_" prefix from the axis labels
  #scale_y_discrete(labels = function(x) gsub("^.*_", "", x)) +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    title = "Site-Specific Top 3 Predictors (Features)",
    subtitle = "Ranked by Permutation Importance (RMSE Increase)",
    x = "Importance (Contribution to Model Error)",
    y = "Predictor"
  ) +
  theme_classic()+ 
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_markdown(size =12, face = "bold",hjust = 0.5, vjust = 0.5),
    plot.subtitle = element_markdown(size = 10, hjust = 0.5, vjust = 0.5),
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 8),
    legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"),
    legend.title = element_blank(), legend.text = element_text(size = 8)
  ) +
  scale_fill_manual(values = palette_site)

Fig_XX_imortance_per_fold

# INTERPRETATION. Feature importance was quantified using permutation importance 
# (increase in RMSE) calculated on the held-out site for each cross-validation 
# fold. This approach identifies predictors that are not only influential during 
# training but are also critical for the model's ability to generalize to new 
# geographic locations.
# Raw Importance (RMSE): Here, a higher number is better. 

# **** FIGURE XXX - Stabilty plot **** -----------------------------------------

Fig_XX_stability_plot <- 
rank_summary %>%
  as.data.frame() %>% 
  mutate(Type = case_when(
    # Low Mean Rank (< 6) and Low SD (< 4.5)
    mean_rank < 6 & sd_rank < 4.5 ~ "Global Driver", 
    # Low Mean Rank (< 6) but High SD (> 4.5)
    mean_rank < 6 & sd_rank >= 4.5 ~ "Context-Dependent",
    TRUE ~ "Low Importance"
  )) %>% 
  ggplot(aes(x = mean_rank, y = sd_rank, color = Type)) +
  # 95% CI shaded region for Global Drivers (optional aesthetic)
  annotate("rect", xmin = -Inf, xmax = 6, ymin = -Inf, ymax = 4.5, 
           fill = "seagreen", alpha = 0.1) +
  geom_point(aes(size = 10-mean_rank), show.legend = FALSE) +
  # Use ggrepel for clean labeling
  geom_text_repel(aes(label = feature), 
                  box.padding = 0.5, point.padding = 0.5, 
                  segment.color = 'black', color = "black", size = 4) +
  # Invert X-axis because Lower Rank = Better
  #scale_x_reverse(limits = c(12, 1), breaks = seq(1, 12, 2)) +
  # Color palette
  scale_color_manual(values = c("goldenrod3", "seagreen4", "grey60")) +
  # Set size based on importance
  scale_size_continuous(range = c(3, 7)) +
  labs(
    title = "Predictor Rank Stability",
    subtitle = "Analysis of LOSO Folds (n=5)",
    x = "Mean Rank (lower = more important)",
    y = expression("Rank Instability (" * SD[rank] * ")")
  ) +
  theme_classic()+ 
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_markdown(size =12, face = "bold",hjust = 0.5, vjust = 0.5),
    plot.subtitle = element_markdown(size = 10, hjust = 0.5, vjust = 0.5),
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 8),
    legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"),
    legend.title = element_blank(), legend.text = element_text(size = 8),
    legend.position = "bottom",
  )

Fig_XX_stability_plot


ggarrange(
  Fig_XX_null_models,
  Fig_XX_feature_rank,
  Fig_XX_imortance_per_fold,
  Fig_XX_stability_plot,
  ncol = 2, nrow = 2
)


# ***********************************************************************-------
# (A) 99% OTUS to GENUS LEVEL --------------------------------------------------




# ***********************************************************************-------
# ***********************************************************************-------
# ***********************************************************************-------

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
resampling_blocked = rsmp("cv", folds = 5)
resampling_blocked <- rsmp("repeated_cv", folds = 4, repeats = 20)

# ── 4. BENCHMARK GRID
design <- benchmark_grid(
  tasks      = task,
  learners   = learners,
  resamplings = list(resampling_blocked)
)

# ── 5. RUN BENCHMARK 
set.seed(31726)
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



