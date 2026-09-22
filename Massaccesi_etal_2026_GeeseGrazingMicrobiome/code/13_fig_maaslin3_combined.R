# ==============================================================================
# 13_fig_maaslin3_combined.R -- MaAsLin3 results as ONE two-panel figure:
#   A = per-marker coefficient summary plots (the content of Fig_10_maaslin3_summary)
#   B = counts of significant associations by direction (the content of Fig_9_maaslin3)
# A is given the larger share of the page; B sits centred underneath it.
#
# This is an ADDITIONAL figure: Fig_9_maaslin3.pdf / Fig_10_maaslin3_summary.pdf and
# the Figures_final/Fig_4, Fig_5 slots they feed are left exactly as they are.
#
# Nothing is refitted here -- the panels are rebuilt from what 11_maaslin3.R already
# wrote to results/maaslin3_<marker>/, so this script is cheap to re-run.
#
# Two colour versions of panel B are written, identical in every other respect:
#   Fig_4_5_maaslin3_combined.pdf     bars in the Geese palette (blue/orange)
#   Fig_4_5_maaslin3_combined_v2.pdf  bars in panel A's own P_FDR colours -- the
#     abundance facets in MaAsLin3's magenta, the prevalence facets in its teal,
#     with the lighter tint of each for "depleted" and full strength for "enriched"
#
# Inputs : results/maaslin3_{fungi,prokaryote}/{all,significant}_results.tsv (11_maaslin3.R)
#          datasets/phyloseq/physeq_{fungi,prokaryote}_clean.rds  (taxonomy labels)
# Outputs: figures/Fig_4_5_maaslin3_combined.pdf
#          figures/Fig_4_5_maaslin3_combined_v2.pdf
# ==============================================================================

# Locate the paper directory: this file is <ROOT>/code/<script>.R
ROOT <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f)) dirname(dirname(normalizePath(f[1]))) else getwd()
})
source(file.path(ROOT, "functions", "00_setup.R"))
suppressPackageStartupMessages(library(maaslin3))

# Same text scale as the final Fig 2 (see code/12_fig_alpha_beta.R), so axis text and
# axis titles match across the manuscript figures.
AX_TEXT   <- 10   # tick labels
AX_TITLE  <- 11   # axis titles
STRIP     <- 11   # facet strip labels (bold)
PAN_TITLE <- 12   # panel titles
TAG       <- list(size = 18, face = "bold")   # A / B panel tags

harmonise <- function(p) {
  p + theme(axis.text     = element_text(size = AX_TEXT),
            axis.title    = element_text(size = AX_TITLE),
            strip.text    = element_text(size = STRIP, face = "bold"),
            legend.text   = element_text(size = AX_TEXT),
            legend.title  = element_text(size = AX_TITLE),
            plot.title    = element_text(size = PAN_TITLE, face = "bold", hjust = 0.5))
}

# ------------------------------------------------------- A: summary plots -----
# Reuses MaAsLin3's own plotting code on the unmodified results on disk, with OTU
# IDs relabelled to "<deepest classified rank> (OTU_ID)" -- same as 11_maaslin3.R.
summary_panel <- function(marker, physeq) {
  out <- maaslin_dir(marker)
  res <- read.delim(file.path(out, "all_results.tsv"))
  res$model[res$model == "abundance"]  <- "linear"
  res$model[res$model == "prevalence"] <- "logistic"
  lbl <- best_taxon_label(physeq)
  res$feature <- ifelse(res$feature %in% names(lbl), lbl[res$feature], res$feature)
  # maaslin3_summary_plot() insists on writing its own PDF here, so the directory
  # has to exist even when only the *_results.tsv files were carried over
  figs_dir <- file.path(OUT_DIR, paste0("maaslin3_", marker), "figures")
  dir.create(figs_dir, showWarnings = FALSE, recursive = TRUE)
  plt <- maaslin3:::maaslin3_summary_plot(
    res, file.path(figs_dir, "summary_plot_taxonomy.pdf"), figs_dir,
    first_n = 25, max_significance = 0.1)
  harmonise(plt$final)
}

p_fungi <- summary_panel("fungi",      readRDS(file.path(PROC_DIR, "physeq_fungi_clean.rds")))
p_prok  <- summary_panel("prokaryote", readRDS(file.path(PROC_DIR, "physeq_prokaryote_clean.rds")))

panel_A <- ggarrange(title_block(p_fungi, "Fungi"),
                     title_block(p_prok,  "Prokaryotes"),
                     ncol = 2, nrow = 1)

# -------------------------------------------------- B: association counts -----
# Both the abundance and the prevalence row of every jointly-significant feature
# are counted (MaAsLin3's intended reporting unit); rows with an unestimable
# coefficient carry no direction and are dropped. Marker names come from
# MARKER_LABELS so they read the same as in Figs 2-3.
read_sig <- function(marker) {
  f <- file.path(maaslin_dir(marker), "significant_results.tsv")
  if (!file.exists(f)) return(NULL)
  read.delim(f) %>% mutate(Marker = marker)
}
sig <- bind_rows(read_sig("fungi"), read_sig("prokaryote")) %>%
  filter(metadata == "Geese", !is.na(coef)) %>%
  mutate(value     = factor(value, levels = GEESE_LEVELS),
         Direction = ifelse(coef > 0, "enriched", "depleted"),
         Marker    = factor(MARKER_LABELS[Marker], levels = MARKER_LABELS))

counts_panel <- function(fill_scale, fill_aes) {
  harmonise(
    ggplot(sig, aes(x = value, fill = {{ fill_aes }})) +
      geom_bar(position = "dodge") +
      facet_grid(Marker ~ model, scales = "free_y") +
      fill_scale +
      labs(title = "MaAsLin3 significant associations (joint FDR q <= 0.1, vs. ungrazed)",
           x = "Geese density", y = "Number of OTUs", fill = NULL) +
      theme_bw() +
      theme(legend.position = "bottom"))
}

# v1: the original Geese palette.
panel_B <- counts_panel(
  scale_fill_manual(values = c(enriched = "#D55E00", depleted = "#56B4E9")), Direction)

# v2: panel A's own colours. MaAsLin3 draws the Abundance P_FDR legend as a ramp
# from #8B008B to white and the Prevalence one from #008B8B to white (see
# maaslin3:::make_coef_plot), so taking a tint off each ramp keeps B on exactly the
# palette A already uses. Bars have to be keyed on model AND direction to carry
# both, since a single discrete scale cannot vary per facet.
tint <- function(col, pct) colorRampPalette(c("white", col))(101)[pct + 1]
ABUND_COL <- "#8B008B"   # magenta = abundance, as in A's Abundance P_FDR legend
PREV_COL  <- "#008B8B"   # teal    = prevalence, as in A's Prevalence P_FDR legend
FILL_KEYS <- c("abundance depleted", "abundance enriched",
               "prevalence depleted", "prevalence enriched")
sig$FillKey <- factor(paste(sig$model, sig$Direction), levels = FILL_KEYS)

panel_B_v2 <- counts_panel(
  scale_fill_manual(values = setNames(c(tint(ABUND_COL, 42), ABUND_COL,
                                        tint(PREV_COL, 42),  PREV_COL), FILL_KEYS),
                    drop = FALSE),
  FillKey)

# B is a 2 x 2 facet grid; stretched across the full width of the summary plots it
# would be all whitespace, so it is padded to the middle ~60% of the page.
row_A <- ggarrange(panel_A, ncol = 1, nrow = 1, labels = "A", font.label = TAG)
assemble <- function(pB) {
  row_B <- ggarrange(NULL, pB, NULL, ncol = 3, widths = c(0.2, 0.6, 0.2),
                     labels = c("", "B", ""), font.label = TAG)
  gridExtra::arrangeGrob(row_A, row_B, ncol = 1,
                         heights = unit.c(unit(1.65, "null"), unit(1, "null")))
}

Fig    <- assemble(panel_B)
Fig_v2 <- assemble(panel_B_v2)
# 18 in rather than the 22 in of the standalone summary figure: the point sizes
# above are absolute, so the wider the page the smaller this figure's text ends up
# relative to the others once it is placed in the manuscript. Below ~16 in the
# beta-coefficient tick labels of the four forest facets start to collide, so 18 is
# about as tight as this content goes.
ggsave(file.path(FIG_DIR, "Fig_4_5_maaslin3_combined.pdf"), Fig,
       width = 18, height = 13, limitsize = FALSE)
message("[fig] figures/Fig_4_5_maaslin3_combined.pdf")

ggsave(file.path(FIG_DIR, "Fig_4_5_maaslin3_combined_v2.pdf"), Fig_v2,
       width = 18, height = 13, limitsize = FALSE)
message("[fig] figures/Fig_4_5_maaslin3_combined_v2.pdf")
