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

# R options --------------------------------------------------------------------
options(scipen = 9999, pillar.sigfig = 6, digits = 6, max.print = 10000000)
#rm(list = ls())

# Check the lib paths ----------------------------------------------------------
.libPaths()

# ********** FUNCTIONS ********** ----------------------------------------------

# Stadardize NCBI blasTAX taxonomy ---------------------------------------------
extract_blasTAX <- function(tax_path,
                            namemap_path){
  
  taxonomy <- read.delim(tax_path, header = TRUE, strip.white = TRUE, sep="\t") 
  name_mapping <- read.delim(namemap_path, header = FALSE, strip.white = TRUE, sep="\t") 
  
  if (length(setdiff(name_mapping$V2, taxonomy$OTU_ID)) > 0) {
    cli::cli_alert_warning(
      paste("No significant similarity found for:", 
            paste(setdiff(name_mapping$V2, taxonomy$OTU_ID), collapse = ", "))
    )
  }
  
  taxonomy_clean <-
    taxonomy %>% 
    rename(Query = 1, Query_Score = 2) %>% 
    mutate_all(~ str_trim(., side = "both")) %>% # remove trailing spaces in both sides
    as.matrix() %>% 
    as.data.frame() %>% 
    dplyr::arrange(Query)
  
  name_mapping_clean <-
    name_mapping %>% 
    rename(Zotu = V1, Query = V2) %>% 
    mutate(Zotu = gsub(">", "", Zotu)) %>% 
    mutate_all(~ str_trim(., side = "both")) %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    dplyr::arrange(Query)
  
  blastax <-
    full_join(name_mapping_clean, taxonomy_clean, by="Query")
  
  return(blastax)
  
}


extract_blasTAX(
  tax_path = file.path(
    data_path,
    "datasets/taxonomy_blast_99.txt"),
  namemap_path = file.path(
    data_path,
    "datasets/name_mapping_99.txt"
  )
)



# Finalize and modify taxonomy table -------------------------------------------

# Check taxonomy table consistency across ranks --------------------------------
CheckTaxonomyConsistency <- function(tax_table, return_long = FALSE) {
  
  # Load safely (no masking issues)
  requireNamespace("dplyr")
  requireNamespace("tibble")
  
  # Helper: valid classification
  is_classified <- function(x) {
    !is.na(x) & x != "" & x != "Unclassified"
  }
  
  # Prepare table
  tax_df <- tax_table %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature")
  
  # Internal function to check conflicts between two ranks
  check_rank_conflict <- function(df, group_rank, test_rank) {
    
    summary_tbl <- df %>%
      dplyr::group_by(.data[[group_rank]]) %>%
      dplyr::filter(
        is_classified(.data[[group_rank]]) &
          dplyr::n_distinct(.data[[test_rank]][is_classified(.data[[test_rank]])]) > 1
      ) %>%
      dplyr::summarise(
        n_values = dplyr::n_distinct(.data[[test_rank]][is_classified(.data[[test_rank]])]),
        values = paste(unique(.data[[test_rank]][is_classified(.data[[test_rank]])]), collapse = "; "),
        OTUs = paste(feature, collapse = "; "),
        .groups = "drop"
      )
    
    long_tbl <- df %>%
      dplyr::group_by(.data[[group_rank]]) %>%
      dplyr::filter(
        is_classified(.data[[group_rank]]) &
          dplyr::n_distinct(.data[[test_rank]][is_classified(.data[[test_rank]])]) > 1
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(feature, dplyr::all_of(group_rank), dplyr::all_of(test_rank))
    
    return(list(summary = summary_tbl, long = long_tbl))
  }
  
  # Run checks
  genus_family   <- check_rank_conflict(tax_df, "Genus", "Family")
  family_order   <- check_rank_conflict(tax_df, "Family", "Order")
  order_class    <- check_rank_conflict(tax_df, "Order", "Class")
  class_phylum   <- check_rank_conflict(tax_df, "Class", "Phylum")
  
  # Return
  if (return_long) {
    return(list(
      genus_family_conflicts = genus_family$long,
      family_order_conflicts = family_order$long,
      order_class_conflicts  = order_class$long,
      class_phylum_conflicts = class_phylum$long
    ))
  } else {
    return(list(
      genus_family_conflicts = genus_family$summary,
      family_order_conflicts = family_order$summary,
      order_class_conflicts  = order_class$summary,
      class_phylum_conflicts = class_phylum$summary
    ))
  }
}



CheckTaxonomyConsistency(taxonomy_99, 
                         return_long=FALSE)

CheckTaxonomyConsistency(taxonomy_99, 
                         return_long=TRUE)

CheckTaxonomyConsistency(taxonomy_99, 
                         return_long=TRUE)$genus_family_conflicts %>% 
  filter(Genus == "Ambispora")


# Replace blanks ---------------------------------------------------------------
blank2na = function(x, na.strings=c('','.','NA','na','N/A','n/a','NaN','nan')) {
  if (is.factor(x)) {
    lab = attr(x, 'label', exact = T)
    labs1 <- attr(x, 'labels', exact = T)
    labs2 <- attr(x, 'value.labels', exact = T)
    # trimws will convert factor to character
    x = trimws(x,'both')
    if (! is.null(lab)) lab = trimws(lab,'both')
    if (! is.null(labs1)) labs1 = trimws(labs1,'both')
    if (! is.null(labs2)) labs2 = trimws(labs2,'both')
    if (!is.null(na.strings)) {
      # convert to NA
      x[x %in% na.strings] = NA
      # also remember to remove na.strings from value labels 
      labs1 = labs1[! labs1 %in% na.strings]
      labs2 = labs2[! labs2 %in% na.strings]
    }
    # the levels will be reset here
    x = factor(x)
    if (! is.null(lab)) attr(x, 'label') <- lab
    if (! is.null(labs1)) attr(x, 'labels') <- labs1
    if (! is.null(labs2)) attr(x, 'value.labels') <- labs2
  } else if (is.character(x)) {
    lab = attr(x, 'label', exact = T)
    labs1 <- attr(x, 'labels', exact = T)
    labs2 <- attr(x, 'value.labels', exact = T)
    # trimws will convert factor to character
    x = trimws(x,'both')
    if (! is.null(lab)) lab = trimws(lab,'both')
    if (! is.null(labs1)) labs1 = trimws(labs1,'both')
    if (! is.null(labs2)) labs2 = trimws(labs2,'both')
    if (!is.null(na.strings)) {
      # convert to NA
      x[x %in% na.strings] = NA
      # also remember to remove na.strings from value labels 
      labs1 = labs1[! labs1 %in% na.strings]
      labs2 = labs2[! labs2 %in% na.strings]
    }
    if (! is.null(lab)) attr(x, 'label') <- lab
    if (! is.null(labs1)) attr(x, 'labels') <- labs1
    if (! is.null(labs2)) attr(x, 'value.labels') <- labs2
  } else {
    x = x
  }
  return(x)
}


# Finalize OTU/ASVs taxonomy ---------------------------------------------------
FinalizeTaxonomy <- function(taxonomy){
  taxonomy$Species <- 
    gsub(" sp ", "", taxonomy$Species)
  taxonomy[] = lapply(taxonomy, blank2na, na.strings=c('','NA','na','N/A','n/a','NaN','nan'))
  lastValue <- function(x) tail(x[!is.na(x)], 1)
  taxonomy$Genus <- as.character(taxonomy$Genus)
  taxonomy[which(is.na(taxonomy$Genus) == FALSE),]$Genus <-
    paste(taxonomy$Genus[is.na(taxonomy$Genus) == FALSE], "sp.", sep = " ")
  last_taxons<- apply(taxonomy[,c(2:8)], 1, lastValue)
  taxonomy$BestMatch <- last_taxons
  taxonomy$BestMatch <-
    gsub("_", " ", taxonomy$BestMatch)
  taxonomy$Taxonomy <-
    paste(taxonomy$Zotu, taxonomy$BestMatch, sep = "-")
  taxonomy$BestMatch <- 
    gsub(" sp.", "", taxonomy$BestMatch)
  taxonomy$Genus <- 
    gsub(" sp.", "", taxonomy$Genus)
  return(taxonomy)
}


FinalizeTaxonomy(taxonomy_99)

# Calculate hill numbers across interations (rarefied hill) --------------------

calc_hill_rarefied <- function(otu_mat, 
                               depth, 
                               n_iter = 100, 
                               q = 0,
                               return_iter = FALSE, 
                               seed = 2026) {
  
  # Filter samples by depth (CRITICAL)
  keep_samples <- rowSums(otu_mat) >= depth
  
  if (any(!keep_samples)) {
    message(sprintf(
      "Filtering samples: keeping %d / %d (removed %d below depth %d)",
      sum(keep_samples), length(keep_samples), sum(!keep_samples), depth
    ))
  }
  
  otu_mat <- otu_mat[keep_samples, , drop = FALSE]
  
  # Check
  if (!q %in% c(0, 1, 2)) {
    stop("q must be 0, 1, or 2")
  }
  
  set.seed(seed)
  
  n_samples <- nrow(otu_mat)
  samples <- rownames(otu_mat)
  
  # helper: compute Hill number for one matrix
  hill_fun <- function(mat) {
    
    if (q == 0) {
      return(vegan::specnumber(mat))
    }
    
    p <- mat / rowSums(mat)
    
    if (q == 1) {
      return(exp(-rowSums(p * log(p + 1e-12))))
    }
    
    if (q == 2) {
      return(1 / rowSums(p^2))
    }
  }
  
  # storage
  hill_mat <- matrix(NA, nrow = n_samples, ncol = n_iter,
                     dimnames = list(samples, paste0("iter_", seq_len(n_iter))))
  
  # iterate rarefaction
  for (i in seq_len(n_iter)) {
    rare_mat <- vegan::rrarefy(otu_mat, sample = depth)
    hill_mat[, i] <- hill_fun(rare_mat)
  }
  
  # return per-iteration if requested
  if (return_iter) {
    return(as.data.frame(hill_mat))
  }
  
  # otherwise return mean + sd
  data.frame(
    sample_id = samples,
    hill_mean = rowMeans(hill_mat),
    hill_sd   = apply(hill_mat, 1, sd),
    q = q,
    depth = depth,
    n_iter = n_iter
  )
}


# Mean + SD (recommended) across iterations
calc_hill_rarefied(otu_mat = otutable_AMF, depth = 5000, n_iter = 100, q = 0)

# Get all iterations (long-format later if needed)
calc_hill_rarefied(otu_mat = otutable_AMF,, depth = 1000, 
                   n_iter = 100, q = 1, return_iter = TRUE)

# Mixed effect models ----------------------------------------------------------

run_lmem <- function(df,
                     response,
                     model_type = c("baseline", "fixslope", "randomslope"),
                     reml = FALSE) {
  
  stopifnot(requireNamespace("lme4", quietly = TRUE))
  
  model_type <- match.arg(model_type)
  
  # capture response as a string, allowing unquoted column names
  resp <- rlang::as_name(rlang::ensym(response))
  
  # build formulas
  if (model_type == "baseline") {
    # Baseline (Random intercept only, no fixed predictors)
    form <- stats::reformulate(
      termlabels = c("1", "(1 | site/plot_rep)"),
      response   = resp
    )
  } else if (model_type == "fixslope") {
    # Fixed slope (Random Intercepts Only)
    form <- stats::reformulate(
      termlabels = c("fert_status", "(1 | site/plot_rep)"),
      response   = resp
    )
  } else if (model_type == "randomslope") {
    # Random Intercept and Random Slope
    form <- stats::reformulate(
      termlabels = c("fert_status", "(fert_status | site)", "(1 | site:plot_rep)"),
      response   = resp
    )
  }
  
  lme4::lmer(form, data = df, REML = reml)
}

run_lmem(alpha_df, hill_0, "fixslope")

# Robust mixed effect models ---------------------------------------------------

run_lmem_robust <- function(df,
                     response,
                     model_type = c("baseline", "fixslope", "randomslope"),
                     reml = FALSE) {
  
  stopifnot(requireNamespace("lme4", quietly = TRUE))
  
  model_type <- match.arg(model_type)
  
  # capture response as a string, allowing unquoted column names
  resp <- rlang::as_name(rlang::ensym(response))
  
  # build formulas
  if (model_type == "baseline") {
    # Baseline (Random intercept only, no fixed predictors)
    form <- stats::reformulate(
      termlabels = c("1", "(1 | site/plot_rep)"),
      response   = resp
    )
  } else if (model_type == "fixslope") {
    # Fixed slope (Random Intercepts Only)
    form <- stats::reformulate(
      termlabels = c("fert_status", "(1 | site/plot_rep)"),
      response   = resp
    )
  } else if (model_type == "randomslope") {
    # Random Intercept and Random Slope
    form <- stats::reformulate(
      termlabels = c("fert_status", "(fert_status | site)", "(1 | site:plot_rep)"),
      response   = resp
    )
  }
  
  robustlmm::rlmer(form, data = df, REML = reml)
}

summary(run_lmem_robust(alpha_df, hill_0, "fixslope"))




# DHARMA diagnostics plots -----------------------------------------------------

# DHARMa provides simulation-based residual diagnostics that are often more 
# informative than raw residual plots, especially in hierarchical models.

diagnostics_dharma <- function(model, 
                               group_var1, 
                               group_var2 = NULL,
                               n = 1000) {
  
  # Generate simulated residuals
  resid <- simulateResiduals(fittedModel = model, n = n)
  
  # Set plot layout and ensure it resets on exit (even if function errors)
  n_plots <- if (!is.null(group_var2)) 6 else 5
  old_par <- par(mfrow = c(n_plots, 1))
  on.exit(par(old_par), add = TRUE)   # always resets, even on error
  
  # 1) KS test for correct distribution of residuals
  testUniformity(resid)
  
  # 2) Dispersion test (under and overdispersion)
  testDispersion(resid)
  
  # 3) Outlier test
  testOutliers(resid, type = "bootstrap")
  
  # 4) KS test against categorical predictor
  testCategorical(resid, group_var1)
  
  # 5) KS test by interaction (only if group_var2 provided)
  if (!is.null(group_var2)) {
    testCategorical(resid, interaction(group_var1, group_var2, drop = TRUE))
  }
  
  # 6) Quantile deviations
  plotResiduals(resid, group = lmer_otu_chem_data[,group_var1])
  
  invisible(resid)
}

diagnostics_dharma(
  model     = pielou_j_betaMM_base,
  group_var1 = alpha_df$site,
  group_var2 = alpha_df$plot_rep
)


# Generate vegan NMDS or PCoA from a phyloseq object -----------------------------------------
generate_ordination <- function(ps, 
                                method = "NMDS",
                                dist_method = "bray",
                                k = 2) {
  
  # 1. Extract OTU table and ensure samples are ROWS for vegan
  otu_tab <- as.data.frame(otu_table(ps))
  if (taxa_are_rows(ps)) {
    otu_tab <- t(otu_tab)
  }
  
  # 2. Run Engine
  if (toupper(method) == "NMDS") {
    # Returns a 'monoMDS' / 'metaMDS' object
    ord <- metaMDS(otu_tab, distance = dist_method, k = k, trymax = 100, trace = FALSE)
    
  } else if (toupper(method) == "PCOA") {
    # Returns a list with points, eigenvalues, etc.
    d <- vegdist(otu_tab, method = dist_method)
    ord <- cmdscale(d, k = k, eig = TRUE)
    
  } else {
    stop("Method must be 'NMDS' or 'PCoA'")
  }
  
  return(ord)
}

generate_ordination(ps = physeq_90_rare, method = "PCOA")


# Plot ordination from vegan objects -------------------------------------------
plot_ordination <- function(ord,
                            meta,
                            col_var,
                            shape_var,
                            ellipse = TRUE,
                            ellipse_level = 0.95,
                            point_size = 3,
                            point_alpha = 0.8,
                            ellipse_type = "t", 
                            legend_inside = FALSE) {
  
  # Ensure factors
  meta[[col_var]]   <- as.factor(meta[[col_var]])
  meta[[shape_var]] <- as.factor(meta[[shape_var]])
  
  # Extract site scores + axis labels
  if (inherits(ord, "metaMDS")) {
    
    scores_df <- as.data.frame(scores(ord, display = "sites"))
    x_col <- colnames(scores_df)[1]
    y_col <- colnames(scores_df)[2]
    x_lab <- "NMDS1"
    y_lab <- "NMDS2"
    stress_label <- paste("Stress =", round(ord$stress, 3))
    
  } else if (is.list(ord) && !is.null(ord$points)) {
    scores_df <- as.data.frame(ord$points)
    colnames(scores_df)[1:2] <- c("Axis.1", "Axis.2")
    x_col <- "Axis.1"
    y_col <- "Axis.2"
    
    if (!is.null(ord$eig)) {
      eig_pos <- ord$eig[ord$eig > 0]
      eig_pct <- round(eig_pos / sum(eig_pos) * 100, 1)
      x_lab <- paste0("PCoA1 (", eig_pct[1], "%)")
      y_lab <- paste0("PCoA2 (", eig_pct[2], "%)")
    } else {
      x_lab <- "Axis.1"
      y_lab <- "Axis.2"
    }
    
    stress_label <- NULL
    
  } else {
    stop("`ord` must be metaMDS or a cmdscale-like object with $points.")
  }
  
  # Attach metadata
  scores_df <- cbind(scores_df, meta)
  
  # Build plot
  p <- ggplot(
    scores_df,
    aes(x = .data[[x_col]],y = .data[[y_col]],
      color = .data[[col_var]],shape = .data[[shape_var]] )) +
    geom_point(size = point_size, alpha = point_alpha) +
    #stat_ellipse(
    #  aes(group = .data[[col_var]], color = .data[[col_var]]),
    #  type = ellipse_type,
    #  level = ellipse_level,
    #  linewidth = 0.6, 
    #  show.legend = FALSE) +
    labs(x = x_lab,
         y = y_lab,
         color = col_var,
         shape = shape_var) +
    theme_classic() +
    theme(
      plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
      strip.text = element_markdown(size = 10, face = "bold"),
      axis.text.x = element_markdown(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
      axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
      strip.background = element_blank(),
      legend.key.height = unit(0.4, "cm"),
      legend.key.width  = unit(0.4, "cm"),
      legend.title = element_blank(),
      legend.text  = element_markdown(size = 8)
    ) +
    guides(
      color = guide_legend(ncol = 1, override.aes = list(shape = 15, size = 3.5)),
      shape = guide_legend(
        ncol = 1,
        override.aes = list(color = "black", size = 2.5)
      )
    )
  
  # Add ellipses
  if (ellipse) {
    p <- p + stat_ellipse(
      aes(group = .data[[col_var]]),
      type = ellipse_type,
      level = ellipse_level,
      linewidth = 0.6,
      show.legend = FALSE
    )
  }
  
  # Add stress if NMDS
  if (!is.null(stress_label)) {
    p <- p + annotate(
      "text",
      x = -Inf,
      y = -Inf,
      label = stress_label,
      hjust = -0.1,   # push slightly right
      vjust = -1.0,   # push slightly up
      size = 3.5
    )
  }
  
  if (legend_inside) {
    p <- p +
      theme(legend.position = c(1, 0.04),        # x, y from bottom-left (0,0) to top-right (1,1)
            legend.justification = c(1, 0),           # anchor point of the legend box
            legend.background = element_blank()) 
    }
  
  return(p)
}

plot_ordination(
  ord = amf_nmds,
  meta = meta_AMF_rare,
  col_var = "site",
  shape_var = "fert_status",
  legend_inside = TRUE
)


# plot_ordination with env fit fitter vectors ----------------------------------

plot_ordination_envft <- function(ord,
                                  meta,
                                  col_var,
                                  shape_var,
                                  env = NULL,
                                  p_threshold = 0.05,
                                  arrow_mul = 1,
                                  ellipse = FALSE,
                                  ellipse_level = 0.95,
                                  point_size = 3,
                                  point_alpha = 0.8,
                                  ellipse_type = "t", 
                                  legend_inside = FALSE) {
  
  # 1. Extract Scores from Ordination Object
  if (inherits(ord, "metaMDS")) {
    scores_df <- as.data.frame(vegan::scores(ord, display = "sites"))
    x_lab <- "NMDS1"
    y_lab <- "NMDS2"
    stress_label <- paste("Stress =", round(ord$stress, 3))
  } else if (is.list(ord) && !is.null(ord$points)) {
    # Handles cmdscale/PCoA
    scores_df <- as.data.frame(ord$points)
    colnames(scores_df)[1:2] <- c("Axis.1", "Axis.2")
    
    if (!is.null(ord$eig)) {
      eig_pos <- ord$eig[ord$eig > 0]
      eig_pct <- round(eig_pos / sum(eig_pos) * 100, 1)
      x_lab <- paste0("PCoA1 (", eig_pct[1], "%)")
      y_lab <- paste0("PCoA2 (", eig_pct[2], "%)")
    } else {
      x_lab <- "Axis.1"
      y_lab <- "Axis.2"
    }
    stress_label <- NULL
  } else {
    stop("`ord` must be metaMDS or a cmdscale-like object with $points.")
  }
  
  # 2. Attach Metadata and ensure factors
  scores_df <- cbind(scores_df, meta)
  scores_df[[col_var]]   <- as.factor(scores_df[[col_var]])
  scores_df[[shape_var]] <- as.factor(scores_df[[shape_var]])
  
  # Set column names for mapping
  x_col <- colnames(scores_df)[1]
  y_col <- colnames(scores_df)[2]
  
  # 3. Build Base Plot
  p <- ggplot(scores_df, aes(x = .data[[x_col]], y = .data[[y_col]], 
                             color = .data[[col_var]], shape = .data[[shape_var]])) +
    geom_point(size = point_size, alpha = point_alpha) +
    labs(x = x_lab, y = y_lab, color = col_var, shape = shape_var) +
    theme_classic() +
    theme(legend.text = element_markdown(),
          axis.text = element_markdown())
  
  # 4. Add Envfit Vectors if provided
  if (!is.null(env)) {
    # Calculate appropriate scaling
    multiplier <- vegan::ordiArrowMul(env) * arrow_mul
    vec_df <- as.data.frame(vegan::scores(env, display = "vectors")) * multiplier
    
    # Add p-values and filter
    vec_df$p_val <- env$vectors$pvals
    vec_df$var   <- rownames(vec_df)
    vec_df <- vec_df[vec_df$p_val <= p_threshold, ]
    
    if (nrow(vec_df) > 0) {
      p <- p +
        geom_segment(data = vec_df,
                     aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
                     arrow = arrow(length = unit(0.2, "cm")),
                     color = "black", linewidth = 0.7, inherit.aes = FALSE) +
        geom_text(data = vec_df,
                  aes(x = Dim1, y = Dim2, label = var),
                  size = 3.5, vjust = -0.7, color = "black", 
                  fontface = "bold", inherit.aes = FALSE)
    }
  }
  
  # 5. Add Ellipses
  if (ellipse) {
    p <- p + stat_ellipse(aes(group = .data[[col_var]]), 
                          type = ellipse_type, level = ellipse_level, 
                          linewidth = 0.6, show.legend = FALSE)
  }
  
  # 6. Formatting and Legend
  if (legend_inside) {
    p <- p + theme(legend.position = c(0.98, 0.02), 
                   legend.justification = c(1, 0),
                   legend.background = element_blank())
  }
  
  if (!is.null(stress_label)) {
    p <- p + annotate("text", x = -Inf, y = -Inf, label = stress_label, 
                      hjust = -0.1, vjust = -1, size = 3.5)
  }
  
  return(p)
}

plot_ordination_envft(ord = pcoa_mean_AMF, 
                     meta = env_fit_chem_data, 
                     col_var = "treatment", 
                     shape_var = "site", 
                     env = env_fit_AMF, 
                     arrow_mul = 4,
                     p_threshold = 0.05)


# plot ordination using phyloseq -----------------------------------------------

plot_ordination_from_phyloseq <- function(ps, 
                                          method = "PCOA", 
                                          dist_method = "bray", 
                                          col_var = NULL, 
                                          shape_var = NULL, 
                                          legend_inside = FALSE) {
  
  # 1. Extract OTU table and ensure samples are ROWS for vegan
  otu_tab <- as.data.frame(otu_table(ps))
  if (taxa_are_rows(ps)) { otu_tab <- t(otu_tab) }
  
  # 2. Run Engine
  if (toupper(method) == "NMDS") {
    ord <- meta_MDS(otu_tab, distance = dist_method, k = 2, trymax = 100, trace = FALSE)
    df_coords <- as.data.frame(scores(ord, display = "sites"))
    axis_labs <- c("NMDS1", "NMDS2")
    subtitle <- paste("Stress:", round(ord$stress, 3))
    
  } else if (toupper(method) == "PCOA") {
    d <- vegdist(otu_tab, method = dist_method)
    ord <- cmdscale(d, k = 2, eig = TRUE)
    df_coords <- as.data.frame(ord$points)
    # Calculate % variation for axes
    var_exp <- round(100 * (ord$eig / sum(ord$eig)), 1)
    axis_labs <- c(paste0("PCoA1 (", var_exp[1], "%)"), 
                   paste0("PCoA2 (", var_exp[2], "%)"))
    
  } else {
    stop("Method must be 'NMDS' or 'PCoA'")
  }
  
  # 3. Combine Coordinates with Metadata
  colnames(df_coords) <- c("Axis1", "Axis2")
  df_coords$ID <- rownames(df_coords)
  
  meta_df <- data.frame(sample_data(ps)) %>% rownames_to_column("ID")
  plot_data <- inner_join(df_coords, meta_df, by = "ID")
  
  # 4. Check if grouping variables exist
  if (!is.null(col_var) && !col_var %in% colnames(plot_data)) {
    stop(paste("Error: Column", col_var, "not found in phyloseq sample_data!"))
  }
  
  # 5. Build Plot
  p <- ggplot(plot_data, aes(x = Axis1, y = Axis2)) +
    geom_point(aes(color = .data[[col_var]], shape = .data[[shape_var]]), 
               size = 3, alpha = 0.7) +
    labs(x = axis_labs[1], 
         y = axis_labs[2], 
         color = col_var, 
         shape = shape_var) +
    theme_classic() +
    guides(
      color = guide_legend(ncol = 1, override.aes = list(shape = 15, size = 3.5)),
      shape = guide_legend(
        ncol = 1,
        override.aes = list(color = "black", size = 2.5)
      ))
  
  if (legend_inside) {
    p <- p + theme(legend.position = "inside", 
                   legend.position.inside = c(0.02, 0.02),
                   legend.justification = c(0, 0), 
                   legend.title = element_blank())
  }
  
  return(p)
}

plot_ordination_from_phyloseq(
  ps = physeq_90_rare, 
  method = "PCoA", 
  col_var = "site_id", 
  shape_var = "fert_status",
  legend_inside = TRUE)


# plot beta-dispersion ---------------------------------------------------------
plot_betadisper <- function(dist_matrix, 
                            grouping, 
                            offset = 0.02, 
                            jitter_width = 0.2, 
                            alpha = 0.5, 
                            point_size = 1.5, 
                            signif_label = NULL, 
                            scale_axis=NULL) {
  
  
  # Run betadisper and Tukey HSD
  bd <- betadisper(dist_matrix, grouping, type = "median")
  tukey <- TukeyHSD(bd)
  print(tukey)
  
  # Extract distances
  disp_df <- data.frame(
    distances = bd$distances,
    group = bd$group
  )
  
  # Calculate means and max for label placement
  means_df <- disp_df %>%
    group_by(group) %>%
    summarise(
      mean_dist = mean(distances),
      max_dist  = max(distances),
      .groups   = "drop"
    )
  
  # Get Tukey letters only if more than 2 groups
  if (nlevels(factor(grouping)) > 2) {
    letters_raw <- multcompLetters(tukey$group[, "p adj"])$Letters
    letters_df <- data.frame(
      group = names(letters_raw),
      letter = letters_raw
    )
    letters_df <- left_join(letters_df, means_df, by = "group")
  } else {
    # For 2 groups, just use the p-value as annotation
    p_val <- tukey$group[, "p adj"]
    letters_df <- means_df %>%
      mutate(letter = c("a", ifelse(p_val < 0.05, "b", "a")))
  }
  
  
  # Plot
  p <- ggplot(disp_df, aes(y = group, x = distances)) +
    geom_jitter(height = 0.2, width = 0, alpha = 0.5, size = 1.5, color = "grey50") +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.7) +
    #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, color = "black") +
    geom_text(data = letters_df,
              aes(y = group, x = max_dist + 0.02, label = letter),
              size = 4) + #, fontface = "bold") +
    theme_classic() +
    theme(
      plot.title = element_markdown(size =12, face = "bold",hjust = 0.5, vjust = 0.5),
      plot.subtitle = element_markdown(size = 10,hjust = 0.5, vjust = 0.5),
      axis.text.x = element_markdown(size = 8, angle = 0, hjust = 0.5),
      axis.text.y = element_markdown(size = 8),
      legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"),
      legend.title = element_blank(), legend.text = element_text(size = 8)
      #legend.margin=ggplot2::margin(0,5,0,0),
      #legend.box.margin=ggplot2::margin(0,5,0,0)
    )+
    guides(color = guide_legend(ncol=1))
  
  # Return both plot and objects invisibly for further use
  #return(invisible(list(
  #  plot    = p,
  #  betadisper = bd,
  #  tukey   = tukey,
  #  data    = disp_df
  #)))
  
  # Add stress if NMDS
  if (!is.null(signif_label)) {
    p <- p + annotate(
      "text",
      x = -Inf,
      y = -Inf,
      label = signif_label,
      hjust = -0.1,   # push slightly right
      vjust = -0.8,   # push slightly up
      size = 2.5
    )
  }
  
  # Scale axis
  if (!is.null(scale_axis)) {
    p <- p +
      scale_x_continuous(limits = c(0, max(disp_df$distances) + offset + 0.05),
                         expand = c(0, 0))
  }
  
  return(p)
}

# Basic usage 
plot_betadisper(dist_matrix = bray_dist, 
                grouping = meta_rare$site, 
                signif_label = "Site, F=6.88, p=001", 
                scale_axis = TRUE)


# Generate a phyloseq object ---------------------------------------------------

generate_phyloseq <- function(otu, metadata, taxonomy, sequences) {
  # convert to data.frame
  otu       <- as.data.frame(otu, stringsAsFactors = FALSE)
  taxonomy  <- as.data.frame(taxonomy, stringsAsFactors = FALSE)
  metadata  <- as.data.frame(metadata, stringsAsFactors = FALSE)
  
  # check alignment
  taxa_missing <- setdiff(rownames(otu), rownames(taxonomy))
  if (length(taxa_missing)) {
    cli::cli_alert_warning("Some OTUs missing in taxonomy: {length(taxa_missing)}")
  }
  sample_missing <- setdiff(colnames(otu), rownames(metadata))
  if (length(sample_missing)) {
    cli::cli_alert_warning("Some samples missing in metadata: {length(sample_missing)}")
  }
  
  # build refseq object first
  refseq_obj <- refseq(sequences)
  
  # build phyloseq
  ps <- phyloseq(
    otu_table(otu, taxa_are_rows = TRUE),
    sample_data(metadata),
    tax_table(as.matrix(taxonomy)),
    refseq_obj
  ) %>%
    phyloseq::prune_taxa(taxa_sums(.) > 0, .) %>%
    phyloseq::prune_samples(sample_sums(.) > 0, .)
  
  cli::cli_alert_success("phyloseq object created with {ntaxa(ps)} taxa and {nsamples(ps)} samples")
  
  ps
}

generate_phyloseq(
  otu = otutable_90,
  metadata = metadata_99,
  taxonomy = taxonomy_90,
  sequences = zotu_90
)


# Multi PEROMANOVA for R2, F, and p --------------------------------------------
extract_adonis <- function(ps) {
  
  # OTU matrix and metadata
  otu <- as(t(otu_table(ps)), "matrix")
  otu <- as.data.frame(otu)
  metadata <- as.data.frame(as(sample_data(ps), "matrix"))
  
  metadata <- 
    metadata %>% 
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
  
  # run PERMANOVA
  ad1 <- adonis2(
    otu ~ fert_status,
    data = metadata,
    method = "bray",
    permutations = how(blocks = metadata$site, nperm = 999),
    by = "margin"
  )
  
  ad2 <- adonis2(
    otu ~ fert_status + site,
    data = metadata,
    method = "bray",
    permutations = 999,
    by = "margin"
  )
  
  bray_dist <- vegdist(otu, method = "bray")
  bd <- betadisper(bray_dist, metadata$site)
  bd_perm <- permutest(bd, permutations = 999)[[1]]

  # dataset name
  ps_name <- deparse(substitute(ps))
  
  res1 <- data.frame(
    dataset = ps_name,
    predictor = "fert_status",
    R2 = ad1$R2[1],
    F = ad1$F[1],
    p = ad1$`Pr(>F)`[1],
    stringsAsFactors = FALSE
  )
  
  res2 <- data.frame(
    dataset = ps_name,
    predictor = "site",
    R2 = ad2$R2[2],
    F = ad2$F[2],
    p = ad2$`Pr(>F)`[2],
    stringsAsFactors = FALSE
  )
  
  bd_summary <- data.frame(
    predictor = "site",
    Df = bd_perm["Groups", "Df"],
    F = bd_perm["Groups", "F"],
    p = bd_perm["Groups", "Pr(>F)"],
    stringsAsFactors = FALSE
  )
  
  return(cbind(res1, res2, bd_summary))
}

extract_adonis(ps = physeq_90_rare)

# Maaslin plotting function - Differential analysis ----------------------------

PlotMaaslin2 <- function(
    maaslin_results,
    physeq_object,
    group_levels   = c("Fertilized", "Control"),         # customizable factor order
    lfc_limits     = NULL,                               # auto-detect if NULL
    show_ns        = TRUE,                               # show non-significant taxa?
    tile_border    = 0.3                                 # geom_tile border size
) {
  
  require(tidyverse)
  require(ggtext)
  require(phyloseq)
  
  # ── 1. Validate inputs ────────────────────────────────────────────────────
  stopifnot(
    is.data.frame(maaslin_results),
    inherits(physeq_object, "phyloseq"),
    !is.null(tax_table(physeq_object))
  )
  
  required_cols <- c("feature", "value", "coef", "qval")
  missing_cols  <- setdiff(required_cols, names(maaslin_results))
  if (length(missing_cols) > 0)
    stop("Missing columns in maaslin_results: ", paste(missing_cols, collapse = ", "))
  
  # ── 2. Extract & join taxonomy ────────────────────────────────────────────
  tax_df <- as.data.frame(as.matrix(tax_table(physeq_object))) %>%
    rownames_to_column("feature")
  
  # ── 3. Clean & process results ────────────────────────────────────────────
  cleaned_results <- maaslin_results %>%
    # Remove duplicate pairwise comparisons
    rowwise() %>%
    mutate(pair = as.character(value)) %>%
    ungroup() %>%
    distinct(feature, pair, .keep_all = TRUE) %>%
    
    # Significance flag
    mutate(
      sig_label = case_when(
        qval < 0.001 ~ "***",
        qval < 0.01  ~ "**",
        qval < 0.05  ~ "*",
        TRUE ~ ""
      ),
      coef_label = if_else(sig_label != "",
                           paste0(round(coef, 2), "\n", sig_label),
                           as.character(round(coef, 2)))
    ) %>%
    
    # Filter non-significant if requested
    { if (!show_ns) filter(., qval <= 0.05) else . } %>%
    
    # Factor ordering for group facets (only levels that exist in data)
    { if ("group" %in% names(.))
      mutate(., group = fct_relevel(group,
                                    intersect(group_levels, unique(group))))
      else . } %>%
    
    # Join taxonomy
    left_join(tax_df, by = "feature") %>%
    
    # Build italic markdown labels — fall back to feature ID if BestMatch is absent
    mutate(
      BestMatch_clean = if_else(
        !is.na(BestMatch) & str_trim(BestMatch) != "",
        str_trim(BestMatch),
        feature
      ),
      Taxonomy = glue::glue("*{BestMatch_clean}* ({str_trim(feature)})")
    ) %>%
    
    # Order taxa by mean absolute coefficient for a cleaner plot
    mutate(Taxonomy = fct_reorder(Taxonomy, abs(coef), .fun = mean))
  
  # ── 4. Auto-scale fill limits ─────────────────────────────────────────────
  if (is.null(lfc_limits)) {
    max_abs <- max(abs(cleaned_results$coef), na.rm = TRUE)
    max_abs <- ceiling(max_abs)          # round up to nearest integer
    lfc_limits <- c(-max_abs, max_abs)
  }
  lfc_breaks <- pretty(lfc_limits, n = 5)
  
  # ── 5. Build plot ─────────────────────────────────────────────────────────
  p <- ggplot(cleaned_results,
              aes(x = pair, y = Taxonomy, fill = coef)) +
    
    geom_tile(color = "grey85", linewidth = tile_border) +
    
    geom_text(
      aes(label = coef_label),
      size      = 3.2,
      lineheight = 0.85,
      show.legend = FALSE
    ) +
    
    scale_fill_gradient2(
      name     = "Log-fold\nchange",
      low      = "#2166AC",   # colorblind-friendlier blue
      mid      = "#F7F7F7",   # neutral white-grey midpoint
      high     = "#B2182B",   # colorblind-friendlier red
      midpoint = 0,
      limits   = lfc_limits,
      breaks   = lfc_breaks,
      oob      = scales::squish   # squish values outside limits rather than NA
    ) +
    
    #labs(
    #  x = NULL,
    #  y = NULL,
    #  caption = "† q < 0.25 | * q < 0.05 | ** q < 0.01 | *** q < 0.001"
    #) +
    
    { if ("group" %in% names(cleaned_results))
      facet_grid(~ group, scales = "free_x", space = "free_x")
      else list() } +
    
    theme_bw(base_size = 11) +
    theme(
      # Facet strips
      strip.text       = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "grey96", color = "grey70"),
      
      # Axes
      axis.text.x      = element_text(angle = 35, hjust = 1, vjust = 1, size = 9),
      axis.text.y      = element_markdown(size = 9),
      axis.ticks       = element_line(linewidth = 0.3),
      
      # Legend
      legend.title     = element_text(size = 9, face = "bold"),
      legend.text      = element_text(size = 8),
      legend.key.height = unit(0.5, "cm"),
      legend.key.width  = unit(0.35, "cm"),
      legend.position  = "right",
      
      # Grid & panel
      panel.grid       = element_blank(),
      panel.border     = element_rect(color = "grey70"),
      
      # Caption
      plot.caption     = element_text(size = 7.5, hjust = 0, color = "grey40")
    )
  
  return(p)
}


str(mlin2_fert_status)

PlotMaaslin2(maaslin_results=mlin2_fert_status$results,
             physeq_object=physeq_AMF_rare)



# Assess clustering threshold through taxon counts -----------------------------

# NOTE. The min_score is based on the S_score that this BLAST generated
# taxonomy has, it won't be present if the classifciation is generated with a 
# different classifier than blasTAX.
assess_clustering_threshold <- function(taxon_name, 
                                        ps_list,
                                        tax_rank = "Species",
                                        score_col = "S_score",
                                        min_score = 0.9999) {
  
  counts <- sapply(ps_list, function(ps) {
    tax_mat <- as(tax_table(ps), "matrix")
    if (!tax_rank %in% colnames(tax_mat)) stop(paste("Tax rank", tax_rank, "not found"))
    if (!score_col %in% colnames(tax_mat)) stop(paste("Score column", score_col, "not found"))
    
    # Filter by minimum score
    idx <- which(tax_mat[, tax_rank] == taxon_name & tax_mat[, score_col] >= min_score)
    length(idx)
  })
  
  data.frame(
    Threshold = names(ps_list),
    OTU_Count = counts,
    Taxon = taxon_name,
    stringsAsFactors = FALSE
  )
}


assess_clustering_threshold(
  tax_rank = "Species",
  taxon_name =  "Paraglomus brasilianum",
  ps_list = list("90" = physeq_90_rare,"91" = physeq_91_rare,"92" = physeq_92_rare,
                 "93" = physeq_93_rare,"94" = physeq_94_rare,"95" = physeq_95_rare,
                 "96" = physeq_96_rare,"97" = physeq_97_rare,"98" = physeq_98_rare,
                 "99" = physeq_99_rare,"100" = physeq_100_rare),
  score_col = "S_score",
  min_score = 0.99999999)


# Calculate how many taxon ASV or OTUs are across site_var ----------------------
assess_taxon_sites <- function(taxon_name, 
                               ps_list, 
                               tax_rank = "Species", 
                               site_var = "site") {
  
  results <- lapply(names(ps_list), function(thresh) {
    ps <- ps_list[[thresh]]
    
    # Taxonomy table
    tax_mat <- as(tax_table(ps), "matrix")
    
    # Identify OTUs for this taxon
    hits <- rownames(tax_mat)[tax_mat[, tax_rank] == taxon_name]
    
    if (length(hits) == 0) {
      return(data.frame(
        Threshold = thresh,
        Total_OTUs = 0,
        stringsAsFactors = FALSE
      ))
    }
    
    # OTU table
    otu_mat <- as(otu_table(ps), "matrix")
    if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
    otu_sub <- otu_mat[hits, , drop = FALSE]
    
    # Metadata
    meta <- as(sample_data(ps), "data.frame")
    
    site_levels <- unique(meta[[site_var]])
    site_assign <- rep(NA_character_, length(hits))
    
    # Assign each OTU to the first site where it occurs
    for (i in seq_along(hits)) {
      otu <- hits[i]
      samples_with_otu <- which(otu_sub[otu, ] > 0)
      sites_with_otu <- unique(meta[[site_var]][samples_with_otu])
      site_assign[i] <- sites_with_otu[1]  # assign to first site
    }
    
    # Count OTUs per site
    site_counts <- table(factor(site_assign, levels = site_levels))
    
    df <- data.frame(
      Threshold = thresh,
      Total_OTUs = length(hits),
      t(as.matrix(site_counts)),
      stringsAsFactors = FALSE
    )
    
    colnames(df)[3:ncol(df)] <- site_levels
    df
  })
  
  do.call(rbind, results)
}

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





# ****************************************************************--------------


# plotting contaminant OTUs ----------------------------------------------------
PlotContam <- function(df, contam){
  # Make phyloseq object of presence-absence in negative controls and true samples
  physeq_pa <- transform_sample_counts(df, function(abund) 1*(abund>0))
  physeq_pa_neg <- subset_samples(physeq_pa, is.neg%in%c("TRUE"))
  physeq_pa_pos <- subset_samples(physeq_pa, is.neg%in%c("FALSE"))
  # Make data.frame of prevalence in positive and negative samples
  df_contam <- data.frame(pa.pos=taxa_sums(physeq_pa_pos), 
                          pa.neg=taxa_sums(physeq_pa_neg),
                          contaminant=contam$contaminant, 
                          Pvalue=contam$p)
  head(df_contam) %T>% print()
  # plotting 
  ggplot(data=df_contam, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
    geom_point(size=2, alpha=0.7) +
    labs(x="Prevalence in negative controls", y="Prevalence in true samples") +
    theme_classic() +
    scale_colour_manual("Contaminant OTUs", values = c("grey", "red")) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.x = element_text(angle =0, size = 8, hjust = 0.5, vjust = 1), 
          axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8)) -> plot_cont
  return(plot_cont)
}

PlotContam(physeq_ITS, contam_ITS)
PlotContam(physeq_16s, contam_16s) 
PlotContam(physeq_18s, contam_18s) 

sample_data(physeq_ITS)$Index <- as.numeric(as.character(seq(nrow(sample_data(physeq_ITS)))))
sample_data(physeq_16s)$Index <- as.numeric(as.character(seq(nrow(sample_data(physeq_16s)))))
sample_data(physeq_18s)$Index <- as.numeric(as.character(seq(nrow(sample_data(physeq_18s)))))

# Function to plot sample depth ------------------------------------------------
PlotDepth <- function(physeq){
  df <- as(sample_data(physeq), "matrix")
  df <- as.data.frame(df)
  # reconvert to numeric
  df$LibSize <- as.numeric(as.character(df$LibSize))
  df$Index <- as.numeric(as.character(df$Index))
  # order
  df <- df[order(df$LibSize), ]
  df$Index <- seq(nrow(df))
  # inspect
  str(df) %T>% print()
  head(df) %T>% print()
  ggplot(data=df, aes(x=Index, y=LibSize, color=is.neg)) +
    geom_point(alpha =0.7, size=2) +
    theme_classic() +
    scale_colour_manual("Negative control", values = c("grey", "red")) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.x = element_text(angle =0, size = 8, hjust = 0.5, vjust = 1), 
          axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8)) -> plot_dist
  return(plot_dist)  
}

PlotDepth(physeq_ITS)



# *** FIGURE S1 - dcontam ------------------------------------------------------
ggarrange(
  ggarrange(PlotDepth(physeq_ITS) +
              labs(title="ITS", 
                   subtitle = "Samples read depth", 
                   x="Sample index", 
                   y="Read number"),
            PlotDepth(physeq_16s) +
              labs(title="16S", 
                   subtitle = "Samples read depth",
                   x="Sample index",
                   y="Read number"),
            PlotDepth(physeq_18s) +
              labs(title="18S", 
                   subtitle = "Samples read depth",
                   x="Sample index",
                   y="Read number"),
            labels = c("A","B","C"),
            widths = c(1,1,1),
            align = "hv" ,
            ncol = 3, 
            nrow = 1, 
            common.legend = TRUE, 
            legend = c("bottom")),
  ggarrange(
    PlotContam(physeq_ITS, contam_ITS) +
      labs(subtitle="Contaminants n. 13"),
    PlotContam(physeq_16s, contam_16s) +
      labs(subtitle="Contaminants n. 38"),
    PlotContam(physeq_18s, contam_18s) +
      labs(subtitle="Contaminants n. 40"),
    widths = c(1,1,1),
    labels = c("C","D","E"),
    align = "hv" ,
    ncol = 3, 
    nrow = 1,
    common.legend = TRUE,
    legend = c("bottom")),
  widths =  c(1, 1.2),
  ncol = 1, 
  nrow = 2) -> Fig_S1

Fig_S1

# BARPLOTS ---------------------------------------------------------------------
palette_30 <-c("#a35151","#2b8c8a","#dba4a4","#111b77","#283dff","#636bb7","#bfc5ff","#195637","#117744","#60ffaf",
               "#b7ffdb","#cc1c1c","#ff0000","#fcb067","#ffe8d3","#d8d6d4","#82807f","#3f3e3d","#560d0d","#825121",
               "#5b5b19","#fcfc00","#ffff9e","#521899","#ae09ea","#fa7efc","#ffb7ef","#a0fffc","#14fffa","#e0f8fc",
               "#FF8000","#000000")






ExtractBar <- function(physeq, last_gen, Rank){
  #tax_table(physeq)[is.na(tax_table(physeq))]<-"Unclassified"
  print(head(tax_table(physeq)))
  top30 <- 
    names(sort(taxa_sums(physeq), TRUE)[1:last_gen]) # get 30 top genera
  physeq_gen <- 
    prune_taxa(top30, physeq)
  print(head(physeq_gen@sam_data))
  df_bar <- 
    physeq_gen %>%
    tax_glom(taxrank = Rank) %>%                     
    transform_sample_counts(function(x) {x/sum(x)} ) %>% 
    psmelt() %>%                                         
    #filter(Abundance > 0.01) %>%                         
    arrange(get(Rank))                                     
  print(levels(df_bar[,Rank]))
  return(df_bar)
}

bar_df_ITS <-
  ExtractBar(ppro_pcoa_ITS, 50, "Genus")
bar_df_ITS$Genus


PlotBar <- function(df, Rank){
  barplot <-
    ggplot(df, aes(x = Description, 
                   y = Abundance, 
                   fill = get(Rank))) + 
    geom_bar(stat = "identity") +
    theme_classic() +
    scale_fill_manual(values =palette_30) +
    #facet_grid(~Niche, scales = "free_x", space="free_x") +
    #facet_wrap(~ Niche ~ Crop, scales = "free_x", ncol=10) + 
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8)) +
    guides(fill = guide_legend(ncol = 1, title = "Taxa"))
  return(barplot)
}

physeq_ITS_gen <-
  ppro_pcoa_ITS %>%
  subset_samples(Niche%in%c("Soil", "Leaf", "Root"))
physeq_ITS_gen <-
  prune_taxa(taxa_sums(physeq_ITS_gen) > 0, physeq_ITS_gen) 
PlotBar(ExtractBar(physeq_ITS_gen, 50, "Genus"), "Genus")

# **********************************************************************--------

# Diagnostics plots ------------------------------------------------------ 
diagnostic_plots <- function(fit) {
  
  # Extract diagnostic data using broom + direct model functions
  diag_data <- 
    data.frame(
      # This grabs the actual response variable used in the model (e.g., hill_0)
      # Change [[2]] to the name of your predictor if you want to color by it specifically
      .treatment_group = model.frame(fit)[, 2],
      .cooksd = cooks.distance(fit),
      .hat = hatvalues(fit),
      .fitted = fitted(fit),
      .resid = residuals(fit),
      .std.resid = residuals(fit, type = "pearson")
      # NOTE. In LMMs, standardized residuals can be computed using the Pearson 
      # method, which accounts for the variance structure of the model. This is 
      # often more appropriate than simple standardization for linear models.
    )
  
  print(head(diag_data))
  
  # 2. Histogram
  p1 <- ggplot(diag_data, aes(x = .resid)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20, color = "grey", fill = "grey") +
    stat_function(
      fun = dnorm, 
      args = list(mean = mean(diag_data$.resid), sd = sd(diag_data$.resid)),
      color = "red"
    ) +
    labs(x = "Residuals", y = "Density", title = "Histogram of Residuals") +
    theme_bw()
  
  # 3. Normal Q-Q
  p2 <- ggplot(diag_data, aes(sample = .std.resid)) +
    stat_qq(shape = 1) + # Fixed: shape goes here
    stat_qq_line(color = "red") +
    labs(title = "Normal Q-Q") +
    theme_bw()
  
  # 4. Scale-Location
  p3 <- ggplot(diag_data, aes(x = .fitted, 
                              y = sqrt(abs(.std.resid)), 
                              color = .treatment_group)) +
    geom_point(shape = 1) + # Fixed: shape goes here
    geom_smooth(method = "lm", formula = 'y ~ x', se = FALSE, color = "red") +
    labs(x = "Fitted Values", y = expression(sqrt("|Pearson Resid|")), 
         title = "Scale-Location") +
    theme_bw() +
    theme(legend.title = element_blank()) +
    scale_color_manual(values = c("grey", "red"))
  
  # 5. Influence (Cook's D)
  cooks_threshold <- 4 / nrow(diag_data)
  
  p4 <- ggplot(diag_data, aes(x = seq_along(.cooksd), y = .cooksd)) +
    geom_segment(aes(xend = seq_along(.cooksd), yend = 0)) +
    geom_point(shape = 1) +
    geom_hline(yintercept = cooks_threshold, linetype = "dashed", color = "red") +
    geom_text(
      aes(label = ifelse(.cooksd > cooks_threshold, seq_along(.cooksd), "")),
      color = "red",
      hjust = -0.2, vjust = -0.5,
      size = 3)+
    labs(x = "Observation Index", y = "Cook's Distance", title = "Influence Check") +
    theme_bw()
  
  # Combine
  return(ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, labels = "AUTO"))
}


diagnostic_plots(run_lmem(alpha_df, hill_0, "fixslope"))


