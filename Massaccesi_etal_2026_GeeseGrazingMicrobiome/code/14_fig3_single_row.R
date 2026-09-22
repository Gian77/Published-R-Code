# ==============================================================================
# 14_fig3_single_row.R -- single-row variant of Fig 3 (per-phylum OTU richness):
#   every significant phylum in ONE row instead of the 6-per-row wrap, with the
#   "Marker: Phylum" strip labels broken over several lines so they fit above the
#   narrower panels.
#
# This is an ADDITIONAL figure: figures/Fig_3_phylum_richness.pdf and the
# Figures_final/Fig_3 slot it feeds are left exactly as they are.
#
# Same numbers as 07_unique_otus.R -- Kruskal-Wallis (BH across phyla) for which
# phyla are shown, pairwise Wilcoxon (BH) for the compact-letter groups -- read
# back from the tables that script wrote rather than recomputed, so the two
# figures can never disagree about which phyla are significant.
#
# Inputs : datasets/phyloseq/physeq_{fungi,prokaryote}_rar.rds
#          results/phylum_richness_tests_{fungi,prokaryote}.csv   (07_unique_otus.R)
# Outputs: figures/Fig_3_phylum_richness_singlerow.pdf
# ==============================================================================

# Locate the paper directory: this file is <ROOT>/code/<script>.R
ROOT <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f)) dirname(dirname(normalizePath(f[1]))) else getwd()
})
source(file.path(ROOT, "functions", "00_setup.R"))

marker_block <- function(marker) {
  tests <- read.csv(file.path(TAB_DIR, sprintf("phylum_richness_tests_%s.csv", marker)))
  sig   <- tests %>% filter(!is.na(padj) & padj <= 0.05) %>% arrange(padj) %>% pull(Phylum)
  if (length(sig) == 0) return(NULL)

  physeq <- readRDS(file.path(PROC_DIR, sprintf("physeq_%s_rar.rds", marker)))
  dat <- phylum_richness(physeq) %>%
    filter(Phylum %in% sig) %>%
    mutate(Marker = tools::toTitleCase(marker))

  letters_df <- dat %>% group_by(Phylum) %>% group_modify(~ {
    lt <- CompSampl(.x, formula(n_otus ~ Geese))
    tibble(Geese = factor(rownames(lt), levels = GEESE_LEVELS),
           Letters = lt$Letters, y = max(.x$n_otus) * 1.08)
  }) %>% ungroup() %>% mutate(Marker = tools::toTitleCase(marker))

  list(data = dat, letters = letters_df, sig = sig)
}

blocks <- Filter(Negate(is.null), list(marker_block("fungi"), marker_block("prokaryote")))
if (length(blocks) == 0) stop("No phylum reached padj <= 0.05 in either marker.")

dat_all      <- bind_rows(lapply(blocks, `[[`, "data"))
letters_all  <- bind_rows(lapply(blocks, `[[`, "letters"))
facet_levels <- unlist(lapply(blocks, function(b) paste(unique(b$data$Marker), b$sig, sep = ": ")))
dat_all      <- dat_all     %>% mutate(Facet = factor(paste(Marker, Phylum, sep = ": "), levels = facet_levels))
letters_all  <- letters_all %>% mutate(Facet = factor(paste(Marker, Phylum, sep = ": "), levels = facet_levels))

# "Prokaryote: candidate division WPS-1" is far wider than a single panel, so the
# marker always gets its own line and the phylum is then wrapped to the panel
# width. str_wrap() is applied to the phylum alone -- run over the whole label it
# would treat the line break after the marker as ordinary whitespace and undo it.
wrap_strip <- function(x) {
  vapply(x, function(s) {
    parts <- strsplit(s, ": ", fixed = TRUE)[[1]]
    paste0(parts[1], ":\n", str_wrap(paste(parts[-1], collapse = ": "), width = 15))
  }, character(1), USE.NAMES = FALSE)
}

n_total <- length(facet_levels)

Fig3row <- ggplot(dat_all, aes(Geese, n_otus)) +
  geom_jitter(position = position_jitter(0.12), size = 2, aes(color = Geese, shape = Geese)) +
  stat_summary(geom = "pointrange",
               fun.min = function(z) quantile(z, 0.25),
               fun.max = function(z) quantile(z, 0.75),
               fun = median, color = "black", shape = 5, size = 0.4, show.legend = FALSE) +
  geom_text(data = letters_all, aes(x = Geese, y = y, label = Letters), size = 3, inherit.aes = FALSE) +
  facet_wrap(~ Facet, scales = "free_y", nrow = 1,
             labeller = labeller(Facet = wrap_strip)) +
  scale_geese_col() + scale_geese_shape() +
  scale_x_discrete(expand = expansion(add = GEESE_X_EXPAND)) +
  # headroom for the compact-letter groups, which sit at 1.08 * max and otherwise
  # touch the top of the shorter panels
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = "Geese", y = "Number of OTUs in phylum") +
  theme_bw() +
  theme(strip.text = element_text(size = 7, face = "bold", lineheight = 1.05),
        axis.text = element_text(size = 7), axis.title = element_text(size = 8),
        legend.position = "bottom", plot.margin = margin(2, 2, 2, 2))

# Same per-panel width as the wrapped version (2.1 in), so the three grazing
# levels stay as readable as they are in Fig 3; the row just gets long.
ggsave(file.path(FIG_DIR, "Fig_3_phylum_richness_singlerow.pdf"), Fig3row,
       width = 2.1 * n_total, height = 4.6, limitsize = FALSE)
message(sprintf("[fig] figures/Fig_3_phylum_richness_singlerow.pdf  (%d panels in one row)", n_total))
