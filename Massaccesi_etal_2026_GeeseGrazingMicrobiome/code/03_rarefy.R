# ==============================================================================
# 03_rarefy.R  -- rarefaction curves + a single even-depth rarefaction used by
#                 the count-based downstream analyses (unique OTUs, indicator
#                 species, heat trees, random forest).
# NOTE on method: alpha diversity (04) and beta diversity (05) do NOT use a
# rarefied count table -- they average the metric / the distance over 100
# rarefactions (alpha_rarefied / avgdist). A single rarefaction here is only a
# convenience integer table for analyses that need fixed, comparable depth.
# Inputs : datasets/phyloseq/physeq_{fungi,prokaryote}_clean.rds
# Outputs: figures/FigS1_rarecurve_fungi.pdf, FigS2_rarecurve_prokaryotes.pdf
#          datasets/phyloseq/physeq_{fungi,prokaryote}_rar.rds
# ==============================================================================

# Locate the paper directory: this file is <ROOT>/code/<script>.R
ROOT <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f)) dirname(dirname(normalizePath(f[1]))) else getwd()
})
source(file.path(ROOT, "functions", "00_setup.R"))

physeq_fungi_clean      <- readRDS(file.path(PROC_DIR, "physeq_fungi_clean.rds"))
physeq_prokaryote_clean <- readRDS(file.path(PROC_DIR, "physeq_prokaryote_clean.rds"))

# --------------------------------------------------------- rarefaction depths
fungi_depth <- min(sample_sums(physeq_fungi_clean))
bact_depth  <- min(sample_sums(physeq_prokaryote_clean))
message(sprintf("rarefaction depth -- fungi: %d | prokaryotes: %d", fungi_depth, bact_depth))

# ------------------------------------ rarefaction curves with ggplot2 (Fig S1/S2)
ggsave(file.path(FIG_DIR, "FigS1_rarecurve_fungi.pdf"),
       rarecurve_ggplot(physeq_fungi_clean, fungi_depth, step = 200, "Fungi (ITS)"),
       width = 6, height = 4.5)
ggsave(file.path(FIG_DIR, "FigS2_rarecurve_prokaryotes.pdf"),
       rarecurve_ggplot(physeq_prokaryote_clean, bact_depth, step = 500, "Prokaryotes (16S)"),
       width = 6, height = 4.5)
message("[fig] figures/FigS1_rarecurve_fungi.pdf, FigS2_rarecurve_prokaryotes.pdf")

# --------------------------------------------- single even-depth rarefaction --

physeq_fungi_rar      <- reformat_physeq_tax(rarefy_once(physeq_fungi_clean,      fungi_depth))
physeq_prokaryote_rar <- reformat_physeq_tax(rarefy_once(physeq_prokaryote_clean, bact_depth))

save_physeq(physeq_fungi_rar,      "physeq_fungi_rar")
save_physeq(physeq_prokaryote_rar, "physeq_prokaryote_rar")
