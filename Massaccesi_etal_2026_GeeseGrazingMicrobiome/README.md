# 📘 Code for Massaccesi et al.

## **Title:** "Soil microbial microbiome in a vineyard grazed with geese"

## 📖 Citation (unpublished)

**Authors:** Massaccesi, Benucci, Massacci
<!-- TODO before publishing: complete author list with affiliations, and the corresponding author -->

Corresponding author: benucci@msu.edu

## **Journal:** _TODO — journal, status and DOI once submitted/accepted._

**Status:** Unpublished. To use this code, please contact the corresponding author.

---

## 🔬 What this is

Soil fungal (ITS) and prokaryotic (16S rRNA) amplicon analysis of an organic vineyard grazed by
geese at three intensities (none / low / high; 6 replicates each, plus extraction and mock
controls). The scripts here are the ones behind the four manuscript figures: they go from the
UPARSE OTU tables to the finished PDFs via decontamination, rarefaction, Hill-number alpha
diversity, PCoA / dispersion / PERMANOVA, per-phylum OTU richness, and MaAsLin3 differential
abundance and prevalence modelling.

> **Study design caveat:** each grazing level corresponds to a *single* vineyard area, with the six
> replicates being subsamples within that area. Grazing is therefore confounded with area, and all
> results are reported as associations, not causal effects of grazing.

---

## 📂 Repository Structure

```
├── code
│   ├── 01_import.R                  raw OTU tables + taxonomy + fasta -> phyloseq objects
│   ├── 02_decontam.R                decontam (prevalence) + drop controls        -> *_clean
│   ├── 03_rarefy.R                  rarefaction curves + even-depth rarefaction  -> *_rar
│   ├── 04_alpha.R                   Hill-number alpha diversity (q = 0, 1, 2)
│   ├── 05_beta.R                    PCoA / CAP, betadisper, PERMANOVA
│   ├── 07_unique_otus.R             per-phylum OTU richness (Kruskal-Wallis + Wilcoxon)
│   ├── 11_maaslin3.R                MaAsLin3 abundance + prevalence models
│   ├── 12_fig_alpha_beta.R          beta panels + PERMANOVA effect sizes (figure assembly)
│   ├── 13_fig_maaslin3_combined.R   MaAsLin3 as one two-panel figure (figure assembly)
│   ├── 14_fig3_single_row.R         per-phylum richness re-laid out in one row
│   ├── install_packages.R           CRAN + Bioconductor installs
│   └── run_all.R                    runs every step above, in order
├── datasets
│   ├── 16S
│   │   ├── otu_table_16s_UPARSE.txt     OTU table (UPARSE)
│   │   ├── taxonomy_16s_RDP.txt         RDP classifier assignments + per-rank confidence
│   │   └── otus.fasta                   OTU representative sequences
│   ├── ITS
│   │   ├── otu_table_ITS_UPARSE_R1.txt  OTU table (UPARSE, forward read)
│   │   ├── otus_R1.fasta                OTU representative sequences
│   │   ├── constax_taxonomy.txt         CONSTAX2 consensus ranks (Kingdom..Species)
│   │   └── rdp_taxonomy.txt             RDP classifier kingdom call (Fungi / mock flag)
│   ├── mapping_16s.txt                  sample metadata, prokaryotes
│   ├── mapping_ITS.txt                  sample metadata, fungi
│   ├── phyloseq
│   │   ├── physeq_{fungi,prokaryote}_filt.rds    imported + prefiltered (01)
│   │   ├── physeq_{fungi,prokaryote}_clean.rds   decontaminated (02)
│   │   ├── physeq_{fungi,prokaryote}_rar.rds     rarefied (03)
│   │   └── panels_{alpha,beta}.rds              saved ggplot panels (04, 05)
│   └── maaslin3_{fungi,prokaryote}
│       ├── all_results.tsv              every MaAsLin3 abundance/prevalence row
│       └── significant_results.tsv      jointly significant features
├── functions
│   └── 00_setup.R                   paths, palettes, plotting theme, shared helpers
├── misc
├── LICENSE
└── README.md
```

Only the files the pipeline actually reads are included — the sequencing outputs also hold UNOISE
and reverse-read (R2) variants, which this analysis does not use. `results/` and `figures/` are
created on the fly and are not tracked.

The gaps in the script numbering (06, 08, 09, 10) are analysis steps of the wider project that none
of the four figures depend on; the numbering is kept as-is so it still matches the order the
pipeline was written and run in.

---

## 🚀 Clone and Use

```bash
# clone the main repo
git clone https://github.com/Gian77/Published-R-Code.git

cd Published-R-Code/Massaccesi_etal_2026_GeeseGrazingMicrobiome

Rscript code/install_packages.R    # first time only
Rscript code/run_all.R             # the whole pipeline
Rscript code/04_alpha.R            # or any single step
```

Paths are resolved from the location of the script, so there is nothing to edit and no working
directory to set — `Rscript code/05_beta.R` works from anywhere. Running a chunk at a time inside
an IDE works too, as long as this folder is the working directory (or `GEESE_PROJ_DIR` points at
it). Results are written to `results/` (CSV, MaAsLin3 output) and `figures/` (PDF).

**Run order.** Scripts `01 → 02 → 03` build the phyloseq objects; `04`, `05`, `07` and `11` are
independent of each other and can be run in any order afterwards; `12`, `13` and `14` only
re-assemble figures. The phyloseq objects and the MaAsLin3 results are shipped in `datasets/`, so
`04`, `05`, `07`, `13` and `14` run straight from a fresh clone. Two exceptions:

- `12_fig_alpha_beta.R` reads `results/adonis_results.csv` and `results/pairwise_adonis.csv`, so
  **`05_beta.R` has to be run first**.
- `11_maaslin3.R` rewrites the MaAsLin3 output under `results/`; `13` uses it when present and
  falls back on the copies in `datasets/` otherwise.

`05_beta.R` (permutation tests over 100 rarefactions) and `11_maaslin3.R` (model fitting, 4 cores)
are the slow steps — minutes, not seconds. Everything else is fast. Note that `Rscript` leaves a
stray `Rplots.pdf` in whatever directory it is launched from.

---

## 🖼️ Figures

| Manuscript figure | Built by | File in `figures/` |
|---|---|---|
| Fig 1 — Hill-number alpha diversity (q = 0, 1, 2), fungi and prokaryotes | `04_alpha.R` | `Fig_1_alpha.pdf` |
| Fig 2 — PCoA + dispersion per marker, with PERMANOVA effect sizes | `05_beta.R` → `12_fig_alpha_beta.R` | `Fig_2_beta_permanova.pdf` |
| Fig 3 — per-phylum OTU richness across grazing levels | `07_unique_otus.R` | `Fig_3_phylum_richness.pdf` |
| Fig 4 — MaAsLin3: coefficient summary (A) and enriched/depleted association counts (B) | `11_maaslin3.R` → `13_fig_maaslin3_combined.R` | `Fig_4_5_maaslin3_combined.pdf` |
| Figs S1–S2 — rarefaction curves | `03_rarefy.R` | `FigS1_rarecurve_fungi.pdf`, `FigS2_rarecurve_prokaryotes.pdf` |

Alternative layouts of the same results, produced alongside the originals:

- `Fig_4_5_maaslin3_combined_v2.pdf` — Figure 4 with panel B recoloured onto panel A's P_FDR
  palette (`13_fig_maaslin3_combined.R`).
- `Fig_9_maaslin3.pdf` and `Fig_10_maaslin3_summary.pdf` — the two halves of Figure 4 as
  standalone plots (`11_maaslin3.R`).
- `Fig_3_phylum_richness_singlerow.pdf` — Figure 3 with all phyla in a single row
  (`14_fig3_single_row.R`).
- `Fig_1_alpha_beta.pdf`, `Fig_2_permanova.pdf`, `Fig_1_alpha_beta_permanova.pdf` —
  alpha+beta merges and the standalone R² barplot (`12_fig_alpha_beta.R`).

---

## 📊 Statistics behind each figure

- **Alpha diversity:** Hill numbers via `vegan::renyi(hill = TRUE)`, averaged over 100 rarefactions
  to the minimum library size — the *metric* is averaged, not the counts. Points are per-sample, the
  open diamond is the median, whiskers are the 25th–75th percentiles, letters are pairwise Wilcoxon
  (BH-adjusted).
- **Beta diversity:** Bray–Curtis averaged over 100 rarefactions (`vegan::avgdist`), PCoA by
  `cmdscale`, dispersion by `betadisper` + `permutest`, PERMANOVA by `adonis2` with 999
  permutations. Omnibus p-values are uncorrected; pairwise contrasts are BH-adjusted.
- **Per-phylum richness:** per-sample OTU counts within each phylum, Kruskal–Wallis per phylum with
  BH across phyla (only p ≤ 0.05 phyla are plotted), then pairwise Wilcoxon (BH) compact letters.
  "Unclassified" is excluded as a database gap.
- **MaAsLin3:** `~ Geese` with TSS normalization, LOG transform, augmentation and standardization,
  `max_significance = 0.1`. A feature is significant on its joint (Fisher-combined abundance +
  prevalence) FDR, and both of its rows are counted in the association-count panel.

---

## 💻 Environment

Reproduced on the MSU ICER HPCC under a conda environment with R 4.5.2:

| package | version |
|---|---|
| R | 4.5.2 |
| phyloseq | 1.54.2 |
| vegan | 2.7.5 |
| ggplot2 | 4.0.3 |
| ggpubr | 0.6.3 |
| ggtext | 0.1.2 |
| maaslin3 | 1.2.0 |
| decontam | 1.30.0 |
| multcompView | 0.1.11 |
| gridExtra | 2.3 |
| Biostrings | 2.78.0 |
| tidyverse | 2.0.0 |

All randomness is seeded (`set.seed(2024)` in `functions/00_setup.R`, plus explicit seeds inside the
rarefaction helpers), so reruns reproduce the shipped results exactly.
