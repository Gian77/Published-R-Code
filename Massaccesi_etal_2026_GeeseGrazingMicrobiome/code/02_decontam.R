# ==============================================================================
# 02_decontam.R  -- identify & remove contaminants (decontam, prevalence method),
#                   then drop the negative-control samples.
# Inputs : datasets/phyloseq/physeq_fungi_filt.rds, physeq_prokaryote_filt.rds
# Outputs: datasets/phyloseq/physeq_fungi_clean.rds, physeq_prokaryote_clean.rds
# ==============================================================================

# Locate the paper directory: this file is <ROOT>/code/<script>.R
ROOT <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f)) dirname(dirname(normalizePath(f[1]))) else getwd()
})
source(file.path(ROOT, "functions", "00_setup.R"))
suppressPackageStartupMessages(library(decontam))

physeq_fungi_filt      <- readRDS(file.path(PROC_DIR, "physeq_fungi_filt.rds"))
physeq_prokaryote_filt <- readRDS(file.path(PROC_DIR, "physeq_prokaryote_filt.rds"))

REAL_LEVELS <- c("none", "low", "high")   # non-control sample groups

# Run decontam (prevalence) using control samples as negatives, drop the
# flagged contaminants, then keep only real samples.  `threshold` differs by
# marker, matching the original analysis (0.2 fungi, 0.5 prokaryotes).
decontaminate <- function(physeq, threshold) {
  sample_data(physeq)$is.neg <- sample_data(physeq)$Geese == "control"
  contam <- isContaminant(physeq, method = "prevalence", neg = "is.neg",
                          threshold = threshold)
  message(sprintf("  contaminants flagged: %d", sum(contam$contaminant)))
  bad   <- rownames(subset(contam, contaminant == TRUE))
  clean <- remove_taxa(physeq, bad)
  # FIX: original used `Geese == c("none","low","high")` (recycled ==, wrong);
  # use %in% so every real sample is retained. prune_samples (not subset_samples)
  # avoids NSE scoping issues when this script is sourced into a fresh env.
  clean <- prune_samples(sample_data(clean)$Geese %in% REAL_LEVELS, clean)
  otu_table(clean) <- otu_table(clean)[rowSums(otu_table(clean)) > 0, ]
  clean
}

message("Fungi:")
physeq_fungi_clean <- decontaminate(physeq_fungi_filt, threshold = 0.2)
message("Prokaryotes:")
physeq_prokaryote_clean <- decontaminate(physeq_prokaryote_filt, threshold = 0.5)

# sanity: no empty taxa/samples
stopifnot(!any(taxa_sums(physeq_fungi_clean) == 0),
          !any(sample_sums(physeq_fungi_clean) == 0),
          !any(taxa_sums(physeq_prokaryote_clean) == 0),
          !any(sample_sums(physeq_prokaryote_clean) == 0))

save_physeq(physeq_fungi_clean,      "physeq_fungi_clean")
save_physeq(physeq_prokaryote_clean, "physeq_prokaryote_clean")
