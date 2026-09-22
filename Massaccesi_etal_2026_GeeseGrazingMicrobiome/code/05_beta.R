# ==============================================================================
# 05_beta.R  -- beta diversity from Bray-Curtis distances AVERAGED over 100
#               rarefactions (vegan::avgdist). Single Figure 2 with PCoA (a,b),
#               CAP (c,d) and betadispersion (e,f); plus PERMANOVA + pairwise.
# Inputs : datasets/phyloseq/physeq_{fungi,prokaryote}_clean.rds  (raw counts)
# Outputs: figures/Fig_2_beta.pdf
#          results/{adonis_results,pairwise_adonis,betadisper_results,pairwise_betadisper}.csv
# ==============================================================================

# Locate the paper directory: this file is <ROOT>/code/<script>.R
ROOT <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f)) dirname(dirname(normalizePath(f[1]))) else getwd()
})
source(file.path(ROOT, "functions", "00_setup.R"))

physeq_fungi_clean      <- readRDS(file.path(PROC_DIR, "physeq_fungi_clean.rds"))
physeq_prokaryote_clean <- readRDS(file.path(PROC_DIR, "physeq_prokaryote_clean.rds"))

fungi_depth <- min(sample_sums(physeq_fungi_clean))
bact_depth  <- min(sample_sums(physeq_prokaryote_clean))

# --- averaged Bray-Curtis distances (100 rarefactions each) -------------------
d_fungi <- avgdist_rarefied(physeq_fungi_clean,      fungi_depth, iters = 100)
d_prok  <- avgdist_rarefied(physeq_prokaryote_clean, bact_depth,  iters = 100)

meta_fungi <- as(sample_data(physeq_fungi_clean),      "data.frame")[labels(d_fungi), ]
meta_prok  <- as(sample_data(physeq_prokaryote_clean), "data.frame")[labels(d_prok), ]
meta_fungi$Geese <- factor(meta_fungi$Geese, levels = GEESE_LEVELS)
meta_prok$Geese  <- factor(meta_prok$Geese,  levels = GEESE_LEVELS)

# --------------------------------------------------- ordinations (Figure 2) ---
cap_fungi <- plot_cap_dist(d_fungi, meta_fungi, "CAP")
cap_prok  <- plot_cap_dist(d_prok,  meta_prok,  "CAP")
message("Fungi CAP anova:");      print(cap_fungi$anova)
message("Prokaryote CAP anova:"); print(cap_prok$anova)

# ordination panels (the merged Figure 2 is assembled at the end, together with
# the dispersion panels)
p_pcoa_fungi <- plot_pcoa_dist(d_fungi, meta_fungi, "PCoA")
p_pcoa_prok  <- plot_pcoa_dist(d_prok,  meta_prok,  "PCoA")
p_cap_fungi  <- cap_fungi$plot
p_cap_prok   <- cap_prok$plot

# ------------------------------------------------------- PERMANOVA (tables) ---
ad_fungi <- adonis_from_dist(d_fungi, meta_fungi, "Geese")
ad_prok  <- adonis_from_dist(d_prok,  meta_prok,  "Geese")
adonis_results <- rbind(
  as.data.frame(ad_fungi) %>% rownames_to_column("Term") %>% mutate(Dataset = "fungi"),
  as.data.frame(ad_prok)  %>% rownames_to_column("Term") %>% mutate(Dataset = "prokaryote"))
write.csv(adonis_results, file.path(TAB_DIR, "adonis_results.csv"), row.names = FALSE)
print(adonis_results)

pairwise_adonis <- rbind(
  pairwise_adonis_from_dist(d_fungi, as.character(meta_fungi$Geese))$contrasts %>% mutate(Dataset = "fungi"),
  pairwise_adonis_from_dist(d_prok,  as.character(meta_prok$Geese))$contrasts  %>% mutate(Dataset = "prokaryote"))
write.csv(pairwise_adonis, file.path(TAB_DIR, "pairwise_adonis.csv"), row.names = FALSE)
print(pairwise_adonis)

# ----------------------------------------------------- betadispersion ---------
bd_fungi <- betadisp_from_dist(d_fungi, as.character(meta_fungi$Geese))
bd_prok  <- betadisp_from_dist(d_prok,  as.character(meta_prok$Geese))
betadisper_results <- rbind(
  bd_fungi[[1]]$tab %>% mutate(Dataset = "fungi",      Model = "Geese"),
  bd_prok[[1]]$tab  %>% mutate(Dataset = "prokaryote", Model = "Geese"))
write.csv(betadisper_results, file.path(TAB_DIR, "betadisper_results.csv"))

pairwise_betadisper <- rbind(
  bd_fungi[[1]]$pairwise$permuted %>% as.data.frame() %>% rownames_to_column("Group") %>%
    separate(Group, c("group1", "group2"), sep = "-") %>% dplyr::rename(pval = 3) %>%
    mutate(padj = round(p.adjust(pval, "BH"), 4), Dataset = "fungi",      Model = "Geese"),
  bd_prok[[1]]$pairwise$permuted  %>% as.data.frame() %>% rownames_to_column("Group") %>%
    separate(Group, c("group1", "group2"), sep = "-") %>% dplyr::rename(pval = 3) %>%
    mutate(padj = round(p.adjust(pval, "BH"), 4), Dataset = "prokaryote", Model = "Geese"))
write.csv(pairwise_betadisper, file.path(TAB_DIR, "pairwise_betadisper.csv"))

# dispersion panels
p_disp_fungi <- PlotBetadisper(d_fungi, as.character(meta_fungi$Geese), 0.65) +
  scale_geese_col() + scale_geese_shape() +
  labs(title = "Dispersion", y = "Distance to centroid", x = NULL)
p_disp_prok  <- PlotBetadisper(d_prok, as.character(meta_prok$Geese), 0.70) +
  scale_geese_col() + scale_geese_shape() +
  labs(title = "Dispersion", y = "Distance to centroid", x = NULL)

# ------------------------------------------------------------- Figure 2 -------
# One figure, 3 columns (PCoA | CAP | Dispersion) x 2 rows (Fungi | Prokaryotes);
# the marker name is shown once as a row label instead of in every panel title.
row_fungi <- ggarrange(p_pcoa_fungi, p_cap_fungi, p_disp_fungi,
                       labels = letters[1:3], align = "hv", ncol = 3, nrow = 1,
                       common.legend = TRUE, legend = "bottom")
row_prok  <- ggarrange(p_pcoa_prok, p_cap_prok, p_disp_prok,
                       labels = letters[4:6], align = "hv", ncol = 3, nrow = 1,
                       common.legend = TRUE, legend = "bottom")

block_fungi <- title_block(row_fungi, "Fungi",       side = "left")
block_prok  <- title_block(row_prok,  "Prokaryotes", side = "left")

Fig2 <- ggarrange(block_fungi, block_prok, ncol = 1, nrow = 2)
ggsave(file.path(FIG_DIR, "Fig_2_beta.pdf"), Fig2, width = 9, height = 6.5)
message("[fig] figures/Fig_2_beta.pdf")

# keep the bare panels for the merged alpha+beta figure (see code/12_fig_alpha_beta.R)
saveRDS(list(fungi      = list(p_pcoa_fungi, p_cap_fungi, p_disp_fungi),
             prokaryote = list(p_pcoa_prok,  p_cap_prok,  p_disp_prok)),
        file.path(PROC_DIR, "panels_beta.rds"))
