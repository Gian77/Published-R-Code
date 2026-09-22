# ==============================================================================
# 07_unique_otus.R  -- per-sample OTU richness WITHIN each phylum, compared
#   across grazing levels. Only phyla with a SIGNIFICANT difference (Kruskal-
#   Wallis, BH across phyla) are plotted, as median/IQR pointrange (alpha style)
#   with pairwise-Wilcoxon compact-letter groups. "Unclassified" is dropped --
#   it is a taxonomy-DB gap, not a biological phylum. Fungi + prokaryotes are
#   combined into a single figure (one panel per significant phylum).
# Inputs : datasets/phyloseq/physeq_{fungi,prokaryote}_rar.rds
# Outputs: results/phylum_richness_tests_<marker>.csv
#          figures/Fig_3_phylum_richness.pdf
# ==============================================================================

# Locate the paper directory: this file is <ROOT>/code/<script>.R
ROOT <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f)) dirname(dirname(normalizePath(f[1]))) else getwd()
})
source(file.path(ROOT, "functions", "00_setup.R"))

# phylum_richness() now lives in 00_setup.R -- 14_fig3_single_row.R needs it too.

analyse_marker <- function(marker) {
  physeq <- readRDS(file.path(PROC_DIR, sprintf("physeq_%s_rar.rds", marker)))
  pr <- phylum_richness(physeq)

  # Kruskal-Wallis per phylum across grazing levels, BH-adjusted across phyla
  # (Unclassified already excluded, so it no longer dilutes the correction)
  tests <- pr %>% group_by(Phylum) %>%
    summarise(p = tryCatch(kruskal.test(n_otus ~ Geese)$p.value, error = function(e) NA_real_),
              .groups = "drop") %>%
    mutate(padj = p.adjust(p, "BH"))
  write.csv(tests, file.path(TAB_DIR, sprintf("phylum_richness_tests_%s.csv", marker)), row.names = FALSE)

  sig <- tests %>% filter(!is.na(padj) & padj <= 0.05) %>% arrange(padj) %>% pull(Phylum)
  message(sprintf("%s: %d phyla (excl. Unclassified) with significant OTU-richness difference", marker, length(sig)))
  if (length(sig) == 0) return(NULL)

  dat <- pr %>% filter(Phylum %in% sig) %>% mutate(Marker = tools::toTitleCase(marker))

  # pairwise-Wilcoxon compact letters + a y position, per significant phylum
  letters_df <- dat %>% group_by(Phylum) %>% group_modify(~ {
    lt <- CompSampl(.x, formula(n_otus ~ Geese))
    tibble(Geese = factor(rownames(lt), levels = GEESE_LEVELS),
           Letters = lt$Letters, y = max(.x$n_otus) * 1.08)
  }) %>% ungroup() %>% mutate(Marker = tools::toTitleCase(marker))

  list(data = dat, letters = letters_df, sig = sig)
}

fungi <- analyse_marker("fungi")
prok  <- analyse_marker("prokaryote")

# ------------------------------------------------------- Figure 3 (combined) --
# One figure, one panel per significant phylum; fungi phyla first (by padj),
# then prokaryote phyla (by padj); "Marker: Phylum" as the strip label.
blocks <- Filter(Negate(is.null), list(fungi, prok))

if (length(blocks) == 0) {
  message("No phylum reached padj <= 0.05 in either marker -- Fig 3 skipped.")
} else {
  dat_all     <- bind_rows(lapply(blocks, `[[`, "data"))
  letters_all <- bind_rows(lapply(blocks, `[[`, "letters"))
  facet_levels <- unlist(lapply(blocks, function(b) paste(unique(b$data$Marker), b$sig, sep = ": ")))
  dat_all     <- dat_all     %>% mutate(Facet = factor(paste(Marker, Phylum, sep = ": "), levels = facet_levels))
  letters_all <- letters_all %>% mutate(Facet = factor(paste(Marker, Phylum, sep = ": "), levels = facet_levels))

  n_total <- length(facet_levels)
  ncol_use <- min(6, n_total)

  Fig3 <- ggplot(dat_all, aes(Geese, n_otus)) +
    geom_jitter(position = position_jitter(0.12), size = 2, aes(color = Geese, shape = Geese)) +
    stat_summary(geom = "pointrange",
                 fun.min = function(z) quantile(z, 0.25),
                 fun.max = function(z) quantile(z, 0.75),
                 fun = median, color = "black", shape = 5, size = 0.4, show.legend = FALSE) +
    geom_text(data = letters_all, aes(x = Geese, y = y, label = Letters), size = 3, inherit.aes = FALSE) +
    facet_wrap(~ Facet, scales = "free_y", ncol = ncol_use) +
    scale_geese_col() + scale_geese_shape() +
    scale_x_discrete(expand = expansion(add = GEESE_X_EXPAND)) +
    labs(x = "Geese", y = "Number of OTUs in phylum") +
    theme_bw() +
    theme(strip.text = element_text(size = 7, face = "bold"),
          axis.text = element_text(size = 7), axis.title = element_text(size = 8),
          legend.position = "bottom", plot.margin = margin(2, 2, 2, 2))

  ggsave(file.path(FIG_DIR, "Fig_3_phylum_richness.pdf"), Fig3,
         width = 2.1 * ncol_use, height = 2.1 * ceiling(n_total / ncol_use) + 0.7)
  message("[fig] figures/Fig_3_phylum_richness.pdf")
}

# drop the old separate-marker files from earlier runs, now superseded by the
# single combined Fig_3_phylum_richness.pdf
old_files <- file.path(FIG_DIR, c("Fig_3_phylum_fungi.pdf", "Fig_4_phylum_prokaryotes.pdf"))
invisible(file.remove(old_files[file.exists(old_files)]))
