# ==============================================================================
# 12_fig_alpha_beta.R  -- assemble the merged diversity figure (old Fig 1 alpha
#   + old Fig 2 beta in a single 12-panel figure) and the PERMANOVA effect-size
#   figure, both from artefacts already produced by 04_alpha.R / 05_beta.R
#   (no rarefaction or permutation is repeated here).
#
#   Layout of the merged figure: 3 columns x 4 rows, one marker per pair of rows
#     Fungi        a b c = Hill q0/q1/q2      d e f = PCoA / CAP / dispersion
#     Prokaryotes  g h i = Hill q0/q1/q2      j k l = PCoA / CAP / dispersion
#   with ONE shared Geese legend at the bottom instead of four.
#
# Inputs : datasets/phyloseq/panels_{alpha,beta}.rds  (04_alpha.R, 05_beta.R)
#          results/{adonis_results,pairwise_adonis}.csv          (05_beta.R)
# Outputs: figures/Fig_2_beta_permanova.pdf        <- FINAL Fig 2: PCoA + dispersion
#                                                     + flipped grey PERMANOVA (a-e)
#          figures/Fig_1_alpha_beta.pdf            (12-panel alpha+beta merge, a-l)
#          figures/Fig_2_permanova.pdf             (standalone R2 barplot)
#          figures/Fig_1_alpha_beta_permanova.pdf  (13-panel variant, a-m)
# The final set uses figures/Fig_1_alpha.pdf (from 04_alpha.R) as Fig 1 and
# Fig_2_beta_permanova.pdf as Fig 2; the other three files here are alternates.
# ==============================================================================

# Locate the paper directory: this file is <ROOT>/code/<script>.R
ROOT <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f)) dirname(dirname(normalizePath(f[1]))) else getwd()
})
source(file.path(ROOT, "functions", "00_setup.R"))

panels_alpha <- readRDS(file.path(PROC_DIR, "panels_alpha.rds"))
panels_beta  <- readRDS(file.path(PROC_DIR, "panels_beta.rds"))

LAB <- list(size = 13, face = "bold")

# In the stacked layout the per-panel "Geese" x-axis title of the alpha panels
# would run into the titles of the ordination row below it, and the rotated tick
# labels waste vertical space -- none/low/high are short enough to sit flat.
tidy_alpha <- function(p) {
  p + labs(x = NULL) +
    theme(legend.position = "none",
          axis.text.x = element_markdown(angle = 0, hjust = 0.5, vjust = 1),
          plot.margin = margin(2, 3, 2, 3))
}
tidy_beta <- function(p) p + theme(legend.position = "none",
                                   plot.margin = margin(6, 3, 2, 3))

# a row of three panels, legends suppressed (one shared legend is added at the end)
panel_row <- function(panels, labs, tidy) {
  ggarrange(plotlist = lapply(panels, tidy),
            ncol = 3, nrow = 1, align = "hv", labels = labs, font.label = LAB,
            legend = "none")
}

# the single Geese legend, taken from one of the ordination panels
geese_legend <- get_legend(
  panels_beta$fungi[[1]] +
    theme(legend.position = "bottom", legend.box.margin = margin(0, 0, 0, 0)))

marker_block <- function(marker, title, labs_alpha, labs_beta) {
  ggarrange(panel_row(panels_alpha[[marker]], labs_alpha, tidy_alpha),
            panel_row(panels_beta[[marker]],  labs_beta,  tidy_beta),
            ncol = 1, nrow = 2) %>%
    title_block(title, side = "left")
}

block_fungi <- marker_block("fungi",      "Fungi",       letters[1:3], letters[4:6])
block_prok  <- marker_block("prokaryote", "Prokaryotes", letters[7:9], letters[10:12])

# NOTE: assembled with arrangeGrob, not ggarrange(legend.grob=) -- ggarrange
# silently drops already-assembled grobs (like the title_block'd rows) when a
# legend.grob is supplied, which yields a page with only the legend on it.
Fig1 <- gridExtra::arrangeGrob(block_fungi, block_prok, geese_legend,
                               ncol = 1, heights = unit.c(unit(1, "null"), unit(1, "null"),
                                                          unit(0.9, "cm")))
ggsave(file.path(FIG_DIR, "Fig_1_alpha_beta.pdf"), Fig1, width = 9.5, height = 11)
message("[fig] figures/Fig_1_alpha_beta.pdf")

# ------------------------------------------------- PERMANOVA effect sizes -----
r2_df <- adonis_r2_df()
write.csv(r2_df, file.path(TAB_DIR, "adonis_effect_sizes.csv"), row.names = FALSE)
print(r2_df)

p_perm <- plot_adonis_r2(r2_df)
ggsave(file.path(FIG_DIR, "Fig_2_permanova.pdf"), p_perm, width = 5.5, height = 3.6)
message("[fig] figures/Fig_2_permanova.pdf")

# ================================================== FINAL Figure 2 (beta) =====
# Three columns -- PCoA | dispersion | PERMANOVA bars -- x two rows (Fungi,
# Prokaryotes). The PERMANOVA panel is split per marker so each row is
# self-contained; both share one R2 scale so the bars stay comparable between
# markers (grazing explains far more prokaryotic than fungal variation).
#
# All six panels go into ONE 3x2 grid with align = "hv": assembling the two rows
# separately (as the 12-panel merge above still does) lets the plotting regions
# drift apart wherever the y-axis tick labels differ in width, so the panel edges
# no longer line up between rows. The Fungi/Prokaryotes row labels are therefore
# added as a separate left-hand column rather than per row with title_block().
PERM_YLIM <- c(0, 1.18 * max(100 * r2_df$R2))

# One text scale for the whole figure. The ordination/dispersion panels come out
# of 05_beta.R at 7-8 pt, which is too small once the figure is placed in the
# manuscript, so they are bumped here to match the PERMANOVA panel.
AX_TEXT   <- 10   # tick labels
AX_TITLE  <- 11   # axis titles
PAN_TITLE <- 12   # per-panel titles (PCoA / Dispersion / PERMANOVA)
# extra bottom margin so the "PCo1 - xx%" axis title cannot be clipped by the
# shared legend / the row below it
PAN_MARGIN <- margin(2, 3, 9, 3)

# `md = TRUE` for panels whose theme was built with ggtext elements (the
# dispersion panels from PlotBetadisper): ggplot2 refuses to merge an
# element_text over an element_markdown, so the replacement has to be the same
# class as what the panel already carries.
# The sizes are set on axis.text.x/.y and axis.title.x/.y, not on the axis.text /
# axis.title parents: PlotBetadisper() sets the children, and a child always wins
# over the parent however late the parent is added.
restyle <- function(p, md = FALSE) {
  el <- if (md) element_markdown else element_text
  p + theme(plot.title   = el(size = PAN_TITLE, face = "bold", hjust = 0.5),
            axis.title.x = el(size = AX_TITLE),
            axis.title.y = el(size = AX_TITLE),
            axis.text.x  = el(size = AX_TEXT),
            axis.text.y  = el(size = AX_TEXT),
            legend.position = "none",
            plot.margin = PAN_MARGIN)
}

# none/low/high sit at fixed x = 1,2,3, so the only way to pull them together is
# to pad the panel edges harder than the figure-wide GEESE_X_EXPAND does.
tighten_x <- function(p)
  suppressMessages(p + scale_x_discrete(expand = expansion(add = 1.3)))

perm_panel <- function(marker_label, grey) {
  plot_adonis_r2(dplyr::filter(r2_df, Marker == marker_label),
                 title = "PERMANOVA", flip = TRUE,
                 fill_values = setNames(grey, marker_label),
                 label_size = 3.4, text_size = AX_TEXT, title_size = PAN_TITLE,
                 ylim = PERM_YLIM, show_legend = FALSE, caption = NULL) +
    theme(plot.margin = PAN_MARGIN)
}

# legend without the "Geese" title, and matching the enlarged panel text
geese_legend_fig2 <- get_legend(
  restyle(panels_beta$fungi[[1]]) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = AX_TITLE),
          legend.box.margin = margin(0, 0, 0, 0)))

beta_panels <- list(
  restyle(panels_beta$fungi[[1]]),                          # a  PCoA
  tighten_x(restyle(panels_beta$fungi[[3]], md = TRUE)),    # b  dispersion
  perm_panel("Fungi", MARKER_GREYS[["Fungi"]]),             # c  PERMANOVA
  restyle(panels_beta$prokaryote[[1]]),                     # d
  tighten_x(restyle(panels_beta$prokaryote[[3]], md = TRUE)), # e
  perm_panel("Prokaryotes", MARKER_GREYS[["Prokaryotes"]])) # f

beta_grid <- ggarrange(plotlist = beta_panels, ncol = 3, nrow = 2, align = "hv",
                       widths = c(1.2, 0.95, 1.5),
                       labels = letters[1:6], font.label = LAB, legend = "none")

row_titles <- gridExtra::arrangeGrob(
  grid::textGrob("Fungi",       rot = 90, gp = grid::gpar(fontsize = 12, fontface = "bold")),
  grid::textGrob("Prokaryotes", rot = 90, gp = grid::gpar(fontsize = 12, fontface = "bold")),
  ncol = 1)

Fig2 <- gridExtra::arrangeGrob(
  gridExtra::arrangeGrob(row_titles, beta_grid, ncol = 2,
                         widths = unit.c(unit(0.55, "cm"), unit(1, "null"))),
  geese_legend_fig2,
  ncol = 1, heights = unit.c(unit(1, "null"), unit(0.8, "cm")))
ggsave(file.path(FIG_DIR, "Fig_2_beta_permanova.pdf"), Fig2, width = 11, height = 6.4)
message("[fig] figures/Fig_2_beta_permanova.pdf")

# ------------------------------- variant: the same 12 panels + PERMANOVA (m) ---
perm_row <- ggarrange(NULL, p_perm, NULL, ncol = 3, widths = c(0.5, 1, 0.5),
                      labels = c("", "m", ""), font.label = LAB)
Fig1b <- gridExtra::arrangeGrob(Fig1, perm_row, ncol = 1,
                                heights = unit.c(unit(1, "null"), unit(0.33, "null")))
ggsave(file.path(FIG_DIR, "Fig_1_alpha_beta_permanova.pdf"), Fig1b,
       width = 9.5, height = 14)
message("[fig] figures/Fig_1_alpha_beta_permanova.pdf")
