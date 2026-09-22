# ==============================================================================
# 11_maaslin3.R  -- MaAsLin3 differential abundance + prevalence modelling of
#                   Geese grazing level, for fungi and prokaryotes.
# MaAsLin3 models both the abundance (of present features) and the prevalence
# (presence/absence) associations, which DESeq2 does not.
# Inputs : datasets/phyloseq/physeq_{fungi,prokaryote}_clean.rds  (raw counts; TSS inside)
# Outputs: results/maaslin3_<marker>/  (all_results.tsv, significant_results.tsv, plots)
#          figures/Fig_9_maaslin3.pdf  (individually-significant OTU counts by direction)
#          figures/Fig_10_maaslin3_summary.pdf (combined, taxonomy-labelled summary plots)
# ==============================================================================

# Locate the paper directory: this file is <ROOT>/code/<script>.R
ROOT <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f)) dirname(dirname(normalizePath(f[1]))) else getwd()
})
source(file.path(ROOT, "functions", "00_setup.R"))
suppressPackageStartupMessages(library(maaslin3))

run_maaslin3 <- function(physeq, marker) {
  feat <- as.data.frame(t(as.data.frame(otu_table(physeq))))     # samples x OTUs
  meta <- as(sample_data(physeq), "data.frame")
  meta$Geese <- factor(meta$Geese, levels = GEESE_LEVELS)        # reference = "none"

  out <- file.path(OUT_DIR, paste0("maaslin3_", marker))
  dir.create(out, showWarnings = FALSE, recursive = TRUE)

  maaslin3(
    input_data      = feat,
    input_metadata  = meta,
    output          = out,
    formula         = ~ Geese,
    normalization   = "TSS",
    transform       = "LOG",
    augment         = TRUE,
    standardize     = TRUE,
    max_significance = 0.1,
    cores           = 4,
    plot_summary_plot = TRUE,
    plot_associations = TRUE)

  sig <- file.path(out, "significant_results.tsv")
  if (file.exists(sig)) {
    n <- nrow(read.delim(sig))
    message(sprintf("%s: %d significant associations -> %s", marker, n, sig))
  }
}

run_maaslin3(readRDS(file.path(PROC_DIR, "physeq_fungi_clean.rds")),      "fungi")
run_maaslin3(readRDS(file.path(PROC_DIR, "physeq_prokaryote_clean.rds")), "prokaryote")

# ------------------------------------------------------------- Figure 9 -------
# Summarise the significant associations: number of features per grazing level
# and per model (abundance vs prevalence), and the direction of the effect.
# MaAsLin3 reports a feature as significant on its JOINT FDR (Fisher-combined
# abundance + prevalence q-value <= 0.1); significant_results.tsv therefore keeps
# both the abundance and prevalence row for each such feature, and we count both
# models here (this is MaAsLin3's intended reporting unit). Rows with an
# unestimable coefficient (NA -- e.g. logistic separation) carry no direction and
# are dropped.
read_sig <- function(marker) {
  f <- file.path(OUT_DIR, paste0("maaslin3_", marker), "significant_results.tsv")
  if (!file.exists(f)) return(NULL)
  read.delim(f) %>% mutate(Marker = marker)
}
sig <- bind_rows(read_sig("fungi"), read_sig("prokaryote")) %>%
  filter(metadata == "Geese", !is.na(coef)) %>%
  mutate(value     = factor(value, levels = GEESE_LEVELS),
         Direction = ifelse(coef > 0, "enriched", "depleted"),
         Marker    = tools::toTitleCase(Marker))

Fig9 <- ggplot(sig, aes(x = value, fill = Direction)) +
  geom_bar(position = "dodge") +
  facet_grid(Marker ~ model, scales = "free_y") +
  scale_fill_manual(values = c(enriched = "#D55E00", depleted = "#56B4E9")) +
  labs(title = "MaAsLin3 significant associations (joint FDR q <= 0.1, vs. ungrazed)",
       x = "Geese density", y = "Number of OTUs", fill = NULL) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        legend.position = "bottom")
ggsave(file.path(FIG_DIR, "Fig_9_maaslin3.pdf"), Fig9, width = 8, height = 6)
message("[fig] figures/Fig_9_maaslin3.pdf")

# ------------------------------------------------------------- Figure 10 ------
# Combine the per-marker MaAsLin3 summary plots (coefficient forest plots of the
# top significant Geese associations) into one panel, relabeling OTU IDs to
# "<deepest classified rank> (OTU_ID)" so the taxonomy is readable directly off
# the plot. Reuses MaAsLin3's own internal plotting code on the unmodified
# results/models on disk -- same content as results/maaslin3_<marker>/figures/
# summary_plot.png, just combined and taxonomy-labelled.
summary_panel <- function(marker, physeq) {
  out <- file.path(OUT_DIR, paste0("maaslin3_", marker))
  res <- read.delim(file.path(out, "all_results.tsv"))
  res$model[res$model == "abundance"]  <- "linear"
  res$model[res$model == "prevalence"] <- "logistic"
  lbl <- best_taxon_label(physeq)
  res$feature <- ifelse(res$feature %in% names(lbl), lbl[res$feature], res$feature)
  figs_dir <- file.path(out, "figures")
  plt <- maaslin3:::maaslin3_summary_plot(
    res, file.path(figs_dir, "summary_plot_taxonomy.pdf"), figs_dir,
    first_n = 25, max_significance = 0.1)
  plt$final
}

p_fungi_summary <- summary_panel("fungi",      readRDS(file.path(PROC_DIR, "physeq_fungi_clean.rds")))
p_prok_summary  <- summary_panel("prokaryote", readRDS(file.path(PROC_DIR, "physeq_prokaryote_clean.rds")))

block_fungi10 <- title_block(p_fungi_summary, "Fungi")
block_prok10  <- title_block(p_prok_summary,  "Prokaryotes")
# Wider-than-tall layout with A/B panel tags; the extra width keeps the taxonomy
# feature names on the left readable while the reduced height makes a more
# manuscript-friendly (landscape) figure.
Fig10 <- ggarrange(block_fungi10, block_prok10, ncol = 2, nrow = 1,
                   labels = c("A", "B"), font.label = list(size = 18, face = "bold"))
ggsave(file.path(FIG_DIR, "Fig_10_maaslin3_summary.pdf"), Fig10, width = 22, height = 8.5, limitsize = FALSE)
message("[fig] figures/Fig_10_maaslin3_summary.pdf")
