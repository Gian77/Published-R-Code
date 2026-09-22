# ==============================================================================
# 00_setup.R  -- options, packages, config, and shared helper functions
# Geese vineyard soil microbiome (16S + ITS) pipeline.
# Sourced by every numbered script in code/. Run order: 01 -> 02 -> ... -> 14.
# ==============================================================================

options(scipen = 9999)        # avoid scientific notation
options(max.print = 1e6)

suppressPackageStartupMessages({
  library(phyloseq)
  library(Biostrings)
  library(tidyverse)
  library(vegan)
  library(ggpubr)
  library(ggrepel)
  library(ggtext)
  library(multcompView)
  library(gridExtra)
  library(grid)
})

# ------------------------------------------------------------------ config ----
# Every path is derived from the paper directory (the folder holding code/,
# functions/ and datasets/), so a fresh clone runs with no editing:
#
#     Rscript code/04_alpha.R        # from anywhere -- ROOT comes from --file=
#
# Each script in code/ defines ROOT before sourcing this file. When this file is
# sourced on its own (line-by-line in an IDE, say), set GEESE_PROJ_DIR or make
# the paper directory the working directory.
PROJ_DIR <- local({
  if (exists("ROOT") && is.character(ROOT) && nzchar(ROOT[1]))
    return(normalizePath(ROOT[1], mustWork = FALSE))
  env <- Sys.getenv("GEESE_PROJ_DIR", "")
  if (nzchar(env)) return(normalizePath(env, mustWork = FALSE))
  normalizePath(getwd(), mustWork = FALSE)
})
if (!dir.exists(file.path(PROJ_DIR, "datasets")))
  stop("Cannot find the paper directory: ", PROJ_DIR,
       "\nRun the scripts as 'Rscript code/<script>.R', set the working directory",
       "\nto the folder containing code/ and datasets/, or export GEESE_PROJ_DIR=/path.")

DATA_DIR <- file.path(PROJ_DIR, "datasets")         # everything shipped with the code
PROC_DIR <- file.path(DATA_DIR, "phyloseq")         # saved phyloseq objects (.rds)
TAB_DIR  <- file.path(PROJ_DIR, "results")          # result CSVs        (generated)
FIG_DIR  <- file.path(PROJ_DIR, "figures")          # figure PDFs        (generated)
OUT_DIR  <- file.path(PROJ_DIR, "results")          # misc outputs (maaslin3, etc.)
for (d in c(PROC_DIR, TAB_DIR, FIG_DIR, OUT_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# raw inputs (read by 01_import.R)
ITS_DIR      <- file.path(DATA_DIR, "ITS")
BAC_DIR      <- file.path(DATA_DIR, "16S")
MAP_ITS      <- file.path(DATA_DIR, "mapping_ITS.txt")
MAP_16S      <- file.path(DATA_DIR, "mapping_16s.txt")

# MaAsLin3 results for a marker: a fresh run of 11_maaslin3.R writes them under
# results/, otherwise fall back on the copies shipped in datasets/ -- so 13 can
# rebuild Figure 4 without refitting every model.
maaslin_dir <- function(marker) {
  fresh <- file.path(OUT_DIR, paste0("maaslin3_", marker))
  if (file.exists(file.path(fresh, "all_results.tsv"))) fresh
  else file.path(DATA_DIR, paste0("maaslin3_", marker))
}

# Geese grazing levels, fixed order + consistent colour/shape across all figures
GEESE_LEVELS <- c("none", "low", "high")
GEESE_COLS   <- c(none = "#56B4E9", low = "#009E73", high = "#D55E00")
GEESE_SHAPES <- c(none = 16, low = 17, high = 15)   # solid circle / triangle / square

# convenience scales (apply to any ggplot mapped on Geese)
scale_geese_col   <- function(...) scale_color_manual("Geese", values = GEESE_COLS, ...)
scale_geese_shape <- function(...) scale_shape_manual("Geese", values = GEESE_SHAPES, ...)

# Marker (dataset) colours -- deliberately outside the Geese palette so a bar
# coloured by marker is never confused with a point coloured by grazing level.
MARKER_LEVELS <- c("fungi", "prokaryote")
MARKER_LABELS <- c(fungi = "Fungi", prokaryote = "Prokaryotes")
MARKER_COLS   <- c(Fungi = "#CC79A7", Prokaryotes = "#0072B2")
MARKER_GREYS  <- c(Fungi = "grey75",  Prokaryotes = "grey35")   # greyscale-safe version

# x-axis padding for the none/low/high categories in jitter/pointrange panels.
# LARGER = the 3 columns are pulled closer together (more margin at the panel
# edges, less white space *between* none/low/high). Applied to Fig 1 (alpha),
# Fig 2 dispersion, and Fig 3 (per-phylum richness).
GEESE_X_EXPAND <- 0.85

# Wrap a ggarrange/grob panel row (or column) with a bold group label -- above
# it (side = "top", e.g. "Fungi"/"Prokaryotes" over a row of panels) or to its
# left (side = "left", e.g. the same labels beside a row when panel columns are
# reserved for something else, such as ordination type) -- for combining
# marker-level blocks into one multi-panel figure. ----------------------------
title_block <- function(grob, title, size = 12, side = c("top", "left")) {
  side <- match.arg(side)
  lbl  <- if (side == "top") grid::textGrob(title, gp = grid::gpar(fontsize = size, fontface = "bold"))
          else grid::textGrob(title, rot = 90, gp = grid::gpar(fontsize = size, fontface = "bold"))
  if (side == "top") gridExtra::arrangeGrob(grob, top  = lbl)
  else                gridExtra::arrangeGrob(grob, left = lbl)
}

set.seed(2024)

# =========================================================== helpers ==========

# Append a "Taxonomy" column = "<OTU_ID>-<deepest classified rank>", where the
# deepest rank is the last non-empty value across whatever rank columns exist
# (Kingdom..Species), regardless of how deep the classification goes. ----------
RANK_COLS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
ReformatTaxonomy <- function(taxonomy_tab) {
  ranks     <- intersect(RANK_COLS, names(taxonomy_tab))
  lastValue <- function(x) { x <- x[!is.na(x) & x != ""]; if (length(x)) tail(x, 1) else NA_character_ }
  taxonomy_tab$BestMatch <- apply(taxonomy_tab[, ranks, drop = FALSE], 1, lastValue)
  taxonomy_tab %>%
    unite(OTU_ID, BestMatch, col = Taxonomy, sep = "-", remove = FALSE)
}

# Named vector OTU_ID -> "<deepest classified rank> (OTU_ID)", for relabeling
# OTU IDs with readable taxonomy on axes/plots without touching modeling code. -
best_taxon_label <- function(physeq) {
  tt <- as.data.frame(as.matrix(tax_table(physeq)), stringsAsFactors = FALSE)
  ranks <- intersect(RANK_COLS, names(tt))
  last_value <- function(x) {
    x <- x[!is.na(x) & x != "" & x != "Unclassified"]
    if (length(x)) tail(x, 1) else "Unclassified"
  }
  best <- apply(tt[, ranks, drop = FALSE], 1, last_value)
  setNames(paste0(best, " (", rownames(tt), ")"), rownames(tt))
}

# Apply ReformatTaxonomy to a phyloseq object's tax_table in place. ------------
reformat_physeq_tax <- function(physeq) {
  tt <- physeq@tax_table %>%
    as.matrix() %>% as.data.frame() %>%
    rownames_to_column("OTU_ID") %>%
    mutate(across(everything(), ~ na_if(.x, "Unclassified"))) %>%
    ReformatTaxonomy()
  m <- as.matrix(tt)
  rownames(m) <- tt$OTU_ID          # keep OTU IDs as taxa names (also kept as a column)
  physeq@tax_table <- tax_table(m)
  physeq
}

# Remove a set of taxa from a phyloseq object. ---------------------------------
remove_taxa <- function(physeq, badTaxa) {
  keep <- setdiff(taxa_names(physeq), badTaxa)
  prune_taxa(keep, physeq)
}

# samples x taxa community matrix from a phyloseq object. ----------------------
comm_matrix <- function(physeq) {
  otu <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) otu <- t(otu)
  otu
}

# Single rarefaction to even depth (integer counts) for count-based analyses
# (unique OTUs, indicator species, heat trees, random forest). -----------------
rarefy_once <- function(physeq, depth, seed = 2024) {
  suppressMessages(
    rarefy_even_depth(physeq, sample.size = depth, rngseed = seed,
                      replace = FALSE, trimOTUs = TRUE, verbose = FALSE))
}

# Hill-number alpha diversity (vegan::renyi, hill = TRUE) averaged over `iters`
# rarefactions -- average the METRIC, NOT the counts. Returns per-sample
# effective OTU numbers: q0 = richness, q1 = exp(Shannon), q2 = inverse Simpson.
HILL_Q <- c(0, 1, 2)
alpha_rarefied <- function(physeq, depth, iters = 100, seed = 2024, qs = HILL_Q) {
  comm <- comm_matrix(physeq)                       # samples x taxa
  set.seed(seed)
  acc <- matrix(0, nrow = nrow(comm), ncol = length(qs),
                dimnames = list(rownames(comm), paste0("q", qs)))
  for (i in seq_len(iters)) {
    r <- vegan::rrarefy(comm, sample = depth)
    acc <- acc + as.matrix(vegan::renyi(r, scales = qs, hill = TRUE))
  }
  acc <- acc / iters
  data.frame(SampleID = rownames(comm), acc, check.names = FALSE, row.names = NULL)
}

# Bray-Curtis distance averaged over `iters` rarefactions (vegan::avgdist). -----
avgdist_rarefied <- function(physeq, depth, iters = 100, method = "bray", seed = 2024) {
  set.seed(seed)
  vegan::avgdist(comm_matrix(physeq), sample = depth, iterations = iters,
                 dmethod = method, meanfun = mean)
}

# Pairwise Wilcoxon (BH) -> compact-letter-display groups. --------------------
CompSampl <- function(df, formula) {
  test_CC <- ggpubr::compare_means(formula, data = df, method = "wilcox.test",
                                   p.adjust.method = "BH")
  test_CC  <- as.data.frame(test_CC)[, c(2, 3, 5)]            # group1, group2, p.adj
  test_CC2 <- data.frame(test_CC[, 2], test_CC[, 1], test_CC[, 3])
  colnames(test_CC2) <- c("group1", "group2", "p.adj")
  test_all <- rbind(test_CC, test_CC2)
  dist_CC  <- as.dist(xtabs(test_all[, 3] ~ (test_all[, 2] + test_all[, 1])), diag = TRUE)
  data.frame(multcompLetters(dist_CC)["Letters"])
}

# Compute CLD letters for an alpha metric directly from a phyloseq object. -----
alpha_letters <- function(physeq, metric, group = "Geese") {
  physeq@sam_data %>%
    as.matrix() %>% as.data.frame() %>%
    mutate(!!metric := as.numeric(!!sym(metric))) %>%
    CompSampl(formula(paste(metric, "~", group))) %>%
    pull(Letters)
}

# Jitter + median/IQR plot of an alpha metric, with CLD letters. --------------
PlotRich <- function(physeq, X_var, Y_var, my_labels, labels_y) {
  dataframe <- physeq@sam_data %>%
    as.matrix() %>% as.data.frame() %>% as_tibble() %>%
    mutate(!!Y_var := as.numeric(!!sym(Y_var)),
           !!X_var := factor(!!sym(X_var), levels = GEESE_LEVELS))
  ggplot(dataframe, aes(x = get(X_var), y = !!sym(Y_var))) +
    geom_jitter(position = position_jitter(0.12), size = 2,
                aes(color = get(X_var), shape = get(X_var))) +
    stat_summary(geom = "pointrange",
                 fun.min = function(z) quantile(z, 0.25),
                 fun.max = function(z) quantile(z, 0.75),
                 fun = median, color = "black", shape = 5, size = 0.5,
                 show.legend = FALSE) +
    stat_summary(geom = "text", angle = 0, label = my_labels,
                 fun = max, aes(y = labels_y), size = 3, color = "black") +
    expand_limits(y = 0) +
    scale_x_discrete(expand = expansion(add = GEESE_X_EXPAND)) +   # tighter columns
    theme_bw() +
    theme(plot.title  = element_markdown(size = 10, face = "bold", vjust = 0.5, hjust = 0.5),
          axis.title  = element_markdown(),
          axis.text.x = element_markdown(angle = 33, hjust = 1, vjust = 1),
          axis.text.y = element_markdown(angle = 0, hjust = 0.5),
          plot.margin = margin(4, 4, 4, 4),
          legend.position = "none")
}

# adonis2 (PERMANOVA) from a precomputed distance. Omnibus test -- NO multiple-
# testing correction here; BH is applied only to the pairwise contrasts. -------
adonis_from_dist <- function(dist, meta, group = "Geese", perm = 999) {
  f <- as.formula(paste("dist ~", group))
  adonis2(f, data = meta, permutations = perm, parallel = 8)
}

# Pairwise PERMANOVA across levels of `group_vec` from a precomputed distance. -
# modified from https://gist.github.com/mcgoodman/58c9d1257fd1625954a4ffa1c3301939
pairwise_adonis_from_dist <- function(dist, group_vec, adj = "BH", perm = 999) {
  D      <- as.matrix(dist)
  groups <- as.data.frame(t(combn(unique(group_vec), m = 2)))
  contrasts <- data.frame(group1 = groups$V1, group2 = groups$V2,
                          R2 = NA, F_value = NA, df1 = NA, df2 = NA, p_value = NA)
  for (i in seq(nrow(contrasts))) {
    idx <- group_vec %in% c(contrasts$group1[i], contrasts$group2[i])
    fit <- vegan::adonis2(as.dist(D[idx, idx]) ~ group_vec[idx],
                          perm = perm, parallel = 8)
    contrasts$R2[i]      <- round(fit$R2[1], 3)
    contrasts$F_value[i] <- round(fit[["F"]][1], 3)
    contrasts$df1[i]     <- fit$Df[1]
    contrasts$df2[i]     <- fit$Df[2]
    contrasts$p_value[i] <- fit$`Pr(>F)`[1]
  }
  contrasts$p_value <- round(p.adjust(contrasts$p_value, method = adj), 3)
  list(contrasts = contrasts, "p-value adjustment" = adj, permutations = perm)
}

# betadisper + permutest from a precomputed distance. The omnibus permutest p is
# left uncorrected; BH is applied only to the pairwise contrasts (in 05_beta and
# the dispersion plot below). -------------------------------------------------
betadisp_from_dist <- function(dist, group_vec) {
  disp     <- betadisper(dist, group_vec)
  dist_var <- vegan::permutest(disp, permutations = 999, pairwise = TRUE)
  list(dist_var, NA, disp)
}

# PCoA scatter from a precomputed distance (cmdscale), coloured by Geese. ------
# (no per-sample labels -- the 3 Geese shapes/colours are enough to read the plot)
plot_pcoa_dist <- function(dist, meta, title) {
  cs   <- cmdscale(dist, k = 2, eig = TRUE)
  varp <- round(100 * cs$eig / sum(cs$eig[cs$eig > 0]), 1)
  data.frame(X = cs$points[, 1], Y = cs$points[, 2],
             Geese = factor(meta$Geese, levels = GEESE_LEVELS)) %>%
    ggplot(aes(X, Y, color = Geese, shape = Geese)) +
    geom_point(size = 3.5) +
    theme_bw() + scale_geese_col() + scale_geese_shape() +
    labs(title = title, x = paste0("PCo1 - ", varp[1], "%"),
         y = paste0("PCo2 - ", varp[2], "%")) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 8), axis.text = element_text(size = 7),
          plot.margin = margin(2, 2, 2, 2), legend.position = "right")
}

# Constrained ordination (CAP / dbRDA) of a distance on Geese. ----------------
plot_cap_dist <- function(dist, meta, title) {
  cap <- vegan::capscale(dist ~ Geese, data = meta)
  sc  <- as.data.frame(vegan::scores(cap, display = "sites", choices = 1:2))
  colnames(sc) <- c("X", "Y")
  sc$Geese <- factor(meta$Geese, levels = GEESE_LEVELS)
  list(plot = ggplot(sc, aes(X, Y, color = Geese, shape = Geese)) +
         geom_point(size = 3.5) +
         theme_bw() + scale_geese_col() + scale_geese_shape() +
         labs(title = title, x = "CAP1", y = "CAP2") +
         theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
               axis.title = element_text(size = 8), axis.text = element_text(size = 7),
               plot.margin = margin(2, 2, 2, 2), legend.position = "right"),
       anova = anova(cap, permutations = 999))
}

# Boxplot-style dispersion plot with CLD letters, from a precomputed distance. -
PlotBetadisper <- function(dist, group_vec, y_axis) {
  bd <- betadisp_from_dist(dist, group_vec)
  signif_beta <- data.frame(multcompLetters(
    p.adjust(bd[[1]]$pairwise$observed, method = "BH"))["Letters"])
  data.frame(value    = bd[[1]]$groups,
             distance = bd[[3]]$distances) %>%
    mutate(value = factor(value, levels = GEESE_LEVELS)) %>%
    ggplot(aes(x = value, y = distance)) +
    geom_jitter(position = position_jitter(0.12), size = 2.5,
                aes(color = value, shape = value)) +
    stat_summary(geom = "pointrange",
                 fun.min = function(z) quantile(z, 0.25),
                 fun.max = function(z) quantile(z, 0.75),
                 fun = median, color = "black", shape = 5, size = 0.5,
                 show.legend = FALSE) +
    stat_summary(geom = "text", angle = 0, label = signif_beta$Letters,
                 fun = max, aes(y = y_axis), size = 3, color = "black") +
    scale_x_discrete(expand = expansion(add = GEESE_X_EXPAND)) +   # tighter columns
    theme_bw() +
    theme(plot.title  = element_markdown(size = 9, face = "bold", vjust = 0.5, hjust = 0.5),
          axis.title  = element_markdown(size = 8, hjust = 0.5, vjust = 0.5),
          axis.text.x = element_markdown(size = 7),
          axis.text.y = element_markdown(size = 7),
          plot.margin = margin(2, 2, 2, 2),
          legend.position = "none")
}

# Significance stars for a p-value vector (ns / * / ** / ***). ----------------
p_stars <- function(p) {
  cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("***", "**", "*", "ns"), right = TRUE) %>% as.character()
}

# "none"/"high" -> "none vs high", always ordered along GEESE_LEVELS so the same
# contrast gets the same label whatever order adonis reported it in. ----------
contrast_label <- function(g1, g2) {
  mapply(function(a, b) {
    o <- c(a, b)[order(match(c(a, b), GEESE_LEVELS))]
    paste(o[1], "vs", o[2])
  }, g1, g2, USE.NAMES = FALSE)
}
CONTRAST_LEVELS <- c("overall", "none vs low", "none vs high", "low vs high")

# PERMANOVA effect sizes as a ggplot: R2 (% of variation in Bray-Curtis distance
# explained by grazing level) for the omnibus test and for each pairwise
# contrast, fungi vs prokaryotes side by side, labelled with R2 + significance.
# Reads the CSVs written by 05_beta.R so the figure can never drift from the
# tables. ---------------------------------------------------------------------
adonis_r2_df <- function(adonis_csv = file.path(TAB_DIR, "adonis_results.csv"),
                         pairwise_csv = file.path(TAB_DIR, "pairwise_adonis.csv")) {
  om <- read.csv(adonis_csv, check.names = FALSE) %>%
    filter(Term == "Model") %>%
    transmute(Dataset, Contrast = "overall", R2, p = `Pr(>F)`, Test = "Omnibus")
  pw <- read.csv(pairwise_csv, check.names = FALSE) %>%
    transmute(Dataset, Contrast = contrast_label(group1, group2),
              R2, p = p_value, Test = "Pairwise (BH)")
  bind_rows(om, pw) %>%
    mutate(Contrast = factor(Contrast, levels = CONTRAST_LEVELS),
           Marker   = factor(MARKER_LABELS[Dataset], levels = MARKER_LABELS),
           stars    = p_stars(p),
           label    = sprintf("%.0f%%%s", 100 * R2, ifelse(stars == "ns", "", stars))) %>%
    arrange(Contrast, Marker)
}

# `flip = TRUE` turns it into horizontal bars with "overall" on top; `fill_values`
# takes MARKER_GREYS for a neutral (greyscale-safe) version; the text sizes are
# arguments so the panel can be scaled up when it sits next to ordinations.
# Pass a single-marker `df` (plus a one-value `fill_values` and a shared `ylim`)
# to get one panel per marker on a common scale, e.g. as the third column of a
# per-marker figure row.
plot_adonis_r2 <- function(df = adonis_r2_df(), title = "PERMANOVA (Bray-Curtis)",
                           flip = FALSE, fill_values = MARKER_COLS,
                           label_size = 2.6, text_size = 8, title_size = 10,
                           ylim = NULL, show_legend = TRUE,
                           caption = paste("overall = omnibus test; pairwise p BH-adjusted.",
                                           "* p<0.05  ** p<0.01  *** p<0.001")) {
  sep <- if (flip) length(CONTRAST_LEVELS) - 0.5 else 1.5   # omnibus | pairwise divider
  p <- ggplot(df, aes(x = Contrast, y = 100 * R2, fill = Marker)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65,
             colour = "grey25", linewidth = 0.25) +
    geom_text(aes(label = label), position = position_dodge(width = 0.75),
              size = label_size,
              hjust = if (flip) -0.12 else 0.5,
              vjust = if (flip) 0.5   else -0.35) +
    geom_vline(xintercept = sep, linetype = "dashed", colour = "grey60") +
    scale_fill_manual("", values = fill_values) +
    scale_y_continuous(limits = ylim,
                       expand = expansion(mult = c(0, if (flip) 0.18 else 0.14))) +
    labs(title = title, x = NULL, y = expression("Variation explained, " * R^2 * " (%)"),
         caption = caption) +
    theme_bw() +
    theme(plot.title   = element_text(size = title_size, face = "bold", hjust = 0.5),
          plot.caption = element_text(size = 6, colour = "grey30", hjust = 0.5),
          axis.title   = element_text(size = text_size + 1),
          axis.text.x  = element_text(size = text_size),
          axis.text.y  = element_text(size = text_size),
          legend.position = if (show_legend) "bottom" else "none",
          legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.35, "cm"),
          legend.text  = element_text(size = text_size),
          plot.margin  = margin(4, 4, 4, 4))
  # blank the grid lines that run ALONG the bars (they add nothing); which theme
  # element that is swaps when the panel is flipped
  p <- p + if (flip) theme(panel.grid.major.y = element_blank())
           else      theme(panel.grid.major.x = element_blank())
  if (flip) p <- p + scale_x_discrete(limits = rev(CONTRAST_LEVELS)) + coord_flip()
  p
}

# Per-sample number of OTUs observed in each phylum (presence/absence within a
# sample). "Unclassified" is dropped -- it is a taxonomy-DB gap, not a phylum.
# Used by 07_unique_otus.R (Fig 3) and 14_fig3_single_row.R. -------------------
phylum_richness <- function(physeq) {
  otu <- comm_matrix(physeq)                          # samples x taxa
  tax <- as.data.frame(as.matrix(tax_table(physeq)))
  phy <- tax$Phylum[match(colnames(otu), tax$OTU_ID)]
  phy[is.na(phy) | phy == ""] <- "Unclassified"
  meta <- as(sample_data(physeq), "data.frame")
  pres <- (otu > 0)
  # sum presence per phylum, per sample
  agg <- t(rowsum(t(pres) * 1, group = phy))          # samples x phyla
  as.data.frame(agg) %>% rownames_to_column("SampleID") %>%
    pivot_longer(-SampleID, names_to = "Phylum", values_to = "n_otus") %>%
    mutate(Geese = factor(meta[SampleID, "Geese"], levels = GEESE_LEVELS)) %>%
    filter(Phylum != "Unclassified")
}

# OTUs present in one group but absent in the other (presence/absence). -------
find_unique_otus <- function(physeq, Group, level1, level2) {
  otu  <- as(otu_table(physeq), "matrix")
  meta <- data.frame(sample_data(physeq))
  l1 <- rownames(meta[meta[, Group] == level1, ])
  l2 <- rownames(meta[meta[, Group] == level2, ])
  in1 <- rownames(otu[, l1, drop = FALSE])[apply(otu[, l1, drop = FALSE], 1, function(x) any(x > 0))]
  in2 <- rownames(otu[, l2, drop = FALSE])[apply(otu[, l2, drop = FALSE], 1, function(x) any(x > 0))]
  list(level1 = setdiff(in1, in2), level2 = setdiff(in2, in1))
}

# Per-phylum count of unique OTUs for a pair of groups, as a grouped barplot. --
UniqueDFplot <- function(physeq, Fact, Lev_1, Lev_2, max_y) {
  uni_list <- find_unique_otus(physeq, Fact, Lev_1, Lev_2)
  tax_df <- physeq@tax_table %>% as.matrix() %>% as.data.frame()
  hi_no <- full_join(
    tax_df %>% filter(OTU_ID %in% uni_list$level1) %>% dplyr::count(Phylum),
    tax_df %>% filter(OTU_ID %in% uni_list$level2) %>% dplyr::count(Phylum),
    by = "Phylum") %>%
    rename(Lev_1 = n.x, Lev_2 = n.y) %>%
    mutate(Phylum = if_else(is.na(Phylum), "Unclassified", Phylum))
  hi_no %>%
    pivot_longer(-Phylum, names_to = "Geese", values_to = "Count") %>%
    mutate(Geese = factor(Geese, levels = c("Lev_1", "Lev_2"), labels = c(Lev_1, Lev_2))) %>%
    ggplot(aes(x = Phylum, y = Count, color = Geese, fill = Geese)) +
    geom_bar(stat = "identity",
             position = position_dodge(preserve = "single", width = 0.8), width = 0.6) +
    geom_text(aes(label = Count), position = position_dodge(width = 0.6),
              vjust = -0.5, size = 2) +
    theme_bw() +
    theme(plot.title  = element_markdown(size = 10, face = "bold", vjust = 0.5, hjust = 0.5),
          axis.title  = element_markdown(size = 8),
          legend.title = element_blank(),
          axis.text.x = element_markdown(angle = 33, size = 7, hjust = 1, vjust = 1),
          axis.text.y = element_markdown(angle = 0, size = 7, hjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"),
          legend.position = c(0.08, 0.93),
          legend.justification = c("left", "top"),
          legend.box = element_blank(),
          legend.margin = margin(-10, -10, -10, -10)) +
    scale_y_continuous(limits = c(0, max_y)) +
    scale_fill_manual(values = c("gray28", "gray70")) +
    scale_color_manual(values = c("gray28", "gray70"))
}

# Rarefaction curve as a ggplot (OTUs vs reads), one line per sample coloured
# by Geese, with a dashed line at the rarefaction depth. ----------------------
rarecurve_ggplot <- function(physeq, depth, step, title) {
  comm <- comm_matrix(physeq)                       # samples x taxa (integer)
  rc   <- vegan::rarecurve(comm, step = step, tidy = TRUE)   # cols: Site, Sample, Species
  meta <- as(sample_data(physeq), "data.frame") %>% rownames_to_column("Site")
  rc %>%
    left_join(meta[, c("Site", "Geese")], by = "Site") %>%
    mutate(Geese = factor(Geese, levels = GEESE_LEVELS)) %>%
    ggplot(aes(x = Sample, y = Species, group = Site, color = Geese)) +
    geom_line(linewidth = 0.4) +
    geom_vline(xintercept = depth, linetype = "dashed", colour = "grey40") +
    scale_geese_col() +
    labs(title = title, x = "Number of reads", y = "Number of OTUs") +
    theme_bw() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          legend.position = "bottom")
}

# Save a phyloseq object to datasets/phyloseq/<name>.rds and echo a one-line summary.
save_physeq <- function(physeq, name) {
  path <- file.path(PROC_DIR, paste0(name, ".rds"))
  saveRDS(physeq, path)
  message(sprintf("[saved] %s : %d taxa x %d samples", path, ntaxa(physeq), nsamples(physeq)))
  invisible(path)
}

message("00_setup.R loaded: helpers + config ready.")
