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



# FInalize and modify taxonomy table -------------------------------------------

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


# DHARMA diagnostics plots -----------------------------------------------------

# DHARMa provides simulation-based residual diagnostics that are often more 
# informative than raw residual plots, especially in hierarchical models.

diagnostics_dharma <- function(model, 
                               group_var1, 
                               group_var2 = NULL,
                               n = 1000) {
  
  # Generate simulated residuals
  resid <- simulateResiduals(fittedModel = model, n = n)
  
  par(mfrow = c(5, 1))
  
  # goodness-of-fit tests
  # 1) KS test for correct distribution of residuals - if the overall distribution
  # conforms to expectations.
  testUniformity(resid)
  
  # 2) Dispersion test (under and overdispersion). Tests if the simulated dispersion
  # is equal to the observed dispersion.
  testDispersion(resid)
  
  # 3) Outlier test (observations outside simulation envelope). Tests if there are 
  # more simulation outliers than expected.
  testOutliers(resid, type = "bootstrap")
  
  # 4) KS test for correct distribution within and between groups - tests residuals 
  # against a categorical predictor.
  testCategorical(resid, group_var1)
  
  # KS test by interaction(group_var1, group_var2) (e.g., site:plot)
  if (!is.null(group_var2)) {
    testCategorical(resid, interaction(group_var1, group_var2, drop = TRUE))
  }
  
  par(mfrow = c(1, 1))
  
  # Return residuals invisibly in case user wants to do further checks
  invisible(resid)
}

diagnostics_dharma(
  model     = pielou_j_base,
  group_var1 = alpha_df$site,
  group_var2 = alpha_df$plot_rep
)


# Plot beta-dversity -----------------------------------------------------------

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
  
  if (!is.null(legend_inside)) {
    p <- p +
      theme(legend.position = c(0.02, 0.04),        # x, y from bottom-left (0,0) to top-right (1,1)
            legend.justification = c(0, 0),           # anchor point of the legend box
            legend.background = element_blank()) 
    }
  
  return(p)
}

plot_ordination(
  ord = amf_nmds,
  meta = meta_rare,
  col_var = "site",
  shape_var = "fert_status",
  legend_inside = TRUE
)


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


# ****************************************************************--------------







