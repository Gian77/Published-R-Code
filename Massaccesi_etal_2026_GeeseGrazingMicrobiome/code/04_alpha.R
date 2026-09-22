# ==============================================================================
# 04_alpha.R  -- Hill-number alpha diversity (vegan::renyi, q = 0, 1, 2),
#                averaged over 100 rarefactions (the METRIC is averaged, not the
#                counts) + pairwise Wilcoxon/BH letters -> Figure 1.
#   q0 = OTU richness | q1 = exp(Shannon) | q2 = inverse Simpson  (effective #OTUs)
# Inputs : datasets/phyloseq/physeq_{fungi,prokaryote}_clean.rds  (raw counts)
#          datasets/phyloseq/physeq_{fungi,prokaryote}_rar.rds    (for sample_data)
# Outputs: results/alpha_diversity.csv; re-saves *_rar with Hill metrics;
#          figures/Fig_1_alpha.pdf
# ==============================================================================

# Locate the paper directory: this file is <ROOT>/code/<script>.R
ROOT <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f)) dirname(dirname(normalizePath(f[1]))) else getwd()
})
source(file.path(ROOT, "functions", "00_setup.R"))

depths <- c(fungi      = min(sample_sums(readRDS(file.path(PROC_DIR, "physeq_fungi_clean.rds")))),
            prokaryote = min(sample_sums(readRDS(file.path(PROC_DIR, "physeq_prokaryote_clean.rds")))))

HILL_COLS <- paste0("q", HILL_Q)  # "q0" "q1" "q2"

# Compute Hill numbers from the CLEAN counts, attach to the rarefied object's
# sample_data (same 18 samples) for plotting.
attach_alpha <- function(marker) {
  clean <- readRDS(file.path(PROC_DIR, sprintf("physeq_%s_clean.rds", marker)))
  rar   <- readRDS(file.path(PROC_DIR, sprintf("physeq_%s_rar.rds",   marker)))
  a <- alpha_rarefied(clean, depths[[marker]], iters = 100)
  rownames(a) <- a$SampleID
  sd <- as(sample_data(rar), "data.frame")
  a  <- a[rownames(sd), ]
  for (q in HILL_COLS) sample_data(rar)[[q]] <- a[[q]]
  save_physeq(rar, sprintf("physeq_%s_rar", marker))
  list(rar = rar, alpha = a %>% mutate(Marker = marker))
}

fungi <- attach_alpha("fungi")
prok  <- attach_alpha("prokaryote")
physeq_fungi_rar      <- fungi$rar
physeq_prokaryote_rar <- prok$rar

# tidy alpha table (effective OTU numbers) for the manuscript / supplementary
geese_lu <- bind_rows(
  as(sample_data(physeq_fungi_rar), "data.frame") %>% rownames_to_column("SampleID"),
  as(sample_data(physeq_prokaryote_rar), "data.frame") %>% rownames_to_column("SampleID")
) %>% dplyr::select(SampleID, Geese) %>% distinct()
alpha_tab <- bind_rows(fungi$alpha, prok$alpha) %>% left_join(geese_lu, by = "SampleID")
write.csv(alpha_tab, file.path(TAB_DIR, "alpha_diversity.csv"), row.names = FALSE)

for (m in c("fungi", "prokaryote")) {
  p <- get(sprintf("physeq_%s_rar", m))
  for (q in HILL_COLS) { message(sprintf("%s %s letters:", m, q)); print(alpha_letters(p, q)) }
}

# ---------------------------------------------------------------- Figure 1 ----
# automatically place the CLD letters just above the largest value in each panel
panel <- function(physeq, metric, title, ylab) {
  vals <- as.numeric(as(sample_data(physeq), "data.frame")[[metric]])
  PlotRich(physeq, "Geese", metric, alpha_letters(physeq, metric), max(vals) * 1.08) +
    scale_geese_col() + scale_geese_shape() +
    labs(title = title, y = ylab, x = "Geese") +
    theme(plot.margin = margin(2, 2, 2, 2))
}

hill_titles <- c(q0 = "q = 0 (richness)", q1 = "q = 1 (exp Shannon)", q2 = "q = 2 (inv Simpson)")
panels_fungi <- lapply(HILL_COLS, function(q) panel(physeq_fungi_rar,      q, hill_titles[q], "Effective #OTUs"))
panels_prok  <- lapply(HILL_COLS, function(q) panel(physeq_prokaryote_rar, q, hill_titles[q], "Effective #OTUs"))

# All 6 panels in a single row, split into two titled blocks (Fungi / Prokaryotes)
# so the marker name is shown once instead of repeated in every panel title.
row_fungi <- ggarrange(plotlist = panels_fungi, ncol = 3, nrow = 1, align = "hv",
                        labels = letters[1:3], common.legend = TRUE, legend = "bottom")
row_prok  <- ggarrange(plotlist = panels_prok,  ncol = 3, nrow = 1, align = "hv",
                        labels = letters[4:6], common.legend = TRUE, legend = "bottom")

block_fungi <- title_block(row_fungi, "Fungi")
block_prok  <- title_block(row_prok,  "Prokaryotes")

Fig1 <- ggarrange(block_fungi, block_prok, ncol = 2, nrow = 1)
ggsave(file.path(FIG_DIR, "Fig_1_alpha.pdf"), Fig1, width = 11, height = 3.2)
message("[fig] figures/Fig_1_alpha.pdf")

# keep the bare panels so code/12_fig_alpha_beta.R can re-assemble them together
# with the beta panels without re-running the 100 rarefactions
saveRDS(list(fungi = panels_fungi, prokaryote = panels_prok),
        file.path(PROC_DIR, "panels_alpha.rds"))
