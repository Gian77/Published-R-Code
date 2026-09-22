# ==============================================================================
# run_all.R -- run the whole pipeline in order, from a clone of this directory.
#
#     Rscript code/run_all.R
#
# Every script resolves its paths from the paper directory (this file's parent's
# parent), reads its inputs from datasets/ and writes to results/ and figures/,
# so each one can also be run on its own, in any order after 01-03:
#
#     Rscript code/05_beta.R
#
# Steps 05 (permutation tests) and 11 (MaAsLin3 model fitting) are the slow ones.
# Everything else runs in seconds to minutes.
# ==============================================================================

ROOT <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f)) dirname(dirname(normalizePath(f[1]))) else getwd()
})
CODE_DIR <- file.path(ROOT, "code")

scripts <- file.path(CODE_DIR, c(
  # ---- build the phyloseq objects ----
  "01_import.R",                  # raw OTU tables + taxonomy -> phyloseq (fungi, prokaryotes)
  "02_decontam.R",                # decontam + drop controls          -> *_clean
  "03_rarefy.R",                  # rarefaction curves + even depth   -> *_rar, Figs S1-S2
  # ---- core analyses ----
  "04_alpha.R",                   # Hill-number alpha diversity       -> Fig 1
  "05_beta.R",                    # PCoA, dispersion, PERMANOVA       -> Fig 2 panels + tables
  "07_unique_otus.R",             # per-phylum OTU richness           -> Fig 3
  "11_maaslin3.R",                # MaAsLin3 abundance + prevalence   -> Fig 4 halves
  # ---- figure assembly (no new computation) ----
  "12_fig_alpha_beta.R",          # beta panels + PERMANOVA effect sizes -> final Fig 2
  "13_fig_maaslin3_combined.R",   # MaAsLin3 as one two-panel figure     -> final Fig 4
  "14_fig3_single_row.R"          # Fig 3 re-laid out as a single row
))

for (s in scripts) {
  if (!file.exists(s)) { message(sprintf("[skip] %s (not present)", s)); next }
  message(sprintf("\n========== running %s ==========", basename(s)))
  t0 <- Sys.time()
  source(s, local = new.env())
  message(sprintf("---------- done %s (%.0fs) ----------",
                  basename(s), as.numeric(Sys.time() - t0, units = "secs")))
}
message("\nALL DONE")
