# Install all R packages needed for the geese microbiome pipeline into the `R` conda env.
# Run via: Rscript code/install_packages.R   (long-running)

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(Ncpus = 8)

ensure <- function(pkgs, installer) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message(sprintf("[install] %s ...", p))
      tryCatch(installer(p),
               error = function(e) message(sprintf("[FAIL] %s: %s", p, conditionMessage(e))))
    } else {
      message(sprintf("[ok] %s already installed", p))
    }
  }
}

# --- CRAN packages ---
cran <- c("ggpubr", "multcompView", "ggrepel", "ggtext",
          "lazyeval", "BiocManager")
ensure(cran, function(p) install.packages(p))

# --- Bioconductor packages ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc <- c("decontam", "maaslin3")
ensure(bioc, function(p) BiocManager::install(p, update = FALSE, ask = FALSE))

# --- maaslin3 fallback from GitHub if not on this Bioc release ---
if (!requireNamespace("maaslin3", quietly = TRUE)) {
  message("[install] maaslin3 from GitHub (biobakery/maaslin3) ...")
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  tryCatch(remotes::install_github("biobakery/maaslin3", upgrade = "never"),
           error = function(e) message(sprintf("[FAIL] maaslin3 github: %s", conditionMessage(e))))
}

# --- report ---
all_pkgs <- c(cran, bioc, "phyloseq", "vegan", "tidyverse", "Biostrings")
message("\n==== FINAL STATUS ====")
for (p in unique(all_pkgs)) {
  message(sprintf("%-14s %s", p, requireNamespace(p, quietly = TRUE)))
}
message("DONE")
