#*********************** Bellucci et al 2025 **********************-------------
# Manuscript:     "soprene-emitting Transgenic Tobacco Shapes Root Microbiome and
#                  Enhances Growth of Co-cultivated Non-emitting Plants"
# Authors:         Bellucci M, Mostofa M, Benucci GMN, Kabir A, Khan I, Lombardi M,
#                  Locato V, Bonito G, Loreto F, Sharkey T
# Code Developer:  Gian MN Benucci
# Affiliation:     Michigan State University, ...
# Journal:         ...
# Date:            August 25, 2025
# **********************************************************************--------

options(scipen = 9999, pillar.sigfig = 6, digits = 6, max.print = 10000000)

rm(list = ls())

# Check the lib paths ----------------------------------------------------------
.libPaths()

# **********************************************************************--------
# ***** SETUP ***** ------------------------------------------------------------

# Load R packages --------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  styler,
  phyloseq,
  speedyseq,
  Biostrings,
  ape,
  vegan,
  AICcPermanova,
  ggtext,
  tidyverse,
  ggpubr,
  gridExtra,
  ggrepel,
  scales,
  magrittr,
  cowplot,
  agricolae
)


# Color palettes ---------------------------------------------------------------
palette1 <- c(
  "#FFDB6D", "#C4961A", "#F4EDCA", "#D16103",
  "#C3D7A4", "#52854C", "#4E84C4", "#293352"
)

palette2 <- c(
  "#825121", "#CC2D35", "#FFB400", "#00A6ED",
  "#7FB800", "#2D3142", "#979797", "#F0E442"
)


palette2 <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E11", "#CC79A7"
)

# **********************************************************************--------
# ***** PATHS ***** ------------------------------------------------------------
# datasets ---------------------------------------------------------------------

data_path <- 
  ("/home/gian/Dropbox/6_PROJETCS/2025_isoprene_tobaccomicrobiome_sharkey/datasets")

# results ----------------------------------------------------------------------

results_path <- 
  ("/home/gian/Dropbox/6_PROJETCS/2025_isoprene_tobaccomicrobiome_sharkey/github/")

# **********************************************************************--------
# ***** IMPORT ***** -----------------------------------------------------------

# Import datasets --------------------------------------------------------------

otutable_ITS <-
  read.delim(file.path(data_path, "/data_ITS/otutable_UPARSE_225bp.txt"),
    row.names = 1
  )

head(otutable_ITS)

otutable_16S <-
  read.delim(file.path(data_path, "/data_16S/otutable_UPARSE_225bp.txt"),
    row.names = 1
  )

head(otutable_16S)

# sequences --------------------------------------------------------------------
sequences_ITS <-
  Biostrings::readDNAStringSet(
    file.path(data_path, "/data_ITS/otus_225bp.fasta"),
    format = "fasta",
    seek.first.rec = TRUE,
    use.names = TRUE
  )
sequences_ITS

sequences_16S <-
  Biostrings::readDNAStringSet(
    file.path(data_path, "/data_16S/otus_225bp.fasta"),
    format = "fasta",
    seek.first.rec = TRUE,
    use.names = TRUE
  )
sequences_16S

# taxonomy ---------------------------------------------------------------------
taxonomy_ITS <-
  read.delim(
    file.path(data_path, "/data_ITS/constax_taxonomy.txt"),
    header = TRUE,
    row.names = 1,
    sep = "\t"
  )
head(taxonomy_ITS)

taxonomy_16S <-
  read.delim(
    file.path(data_path, "/data_16S/constax_taxonomy.txt"),
    header = TRUE,
    row.names = 1,
    sep = "\t"
  )
head(taxonomy_16S)

# metadata ---------------------------------------------------------------------

# 7 replication for each genotype and treatment
# soil control means pot with only soil
# fumigated soil mean: pot with only soil fumigated with isoprene
# IE+NE mix means that the soil belongs to the pot where was planted IE and NE together
# IE-mix or NE mix means that root or rhizosphere belong to IE or NE planted in
#    the mix pot (where was planted IE and NE togheter)

metadata_ITS <-
  read.delim(file.path(data_path, "/data_ITS/metadata_ITS.txt"), header = TRUE, sep = "\t") %>%
  column_to_rownames("SampleID")
head(metadata_ITS)

metadata_16S <-
  read.delim(file.path(data_path, "/data_16S//metadata_16S.txt"), header = TRUE, sep = "\t") %>%
  column_to_rownames("SampleID")
head(metadata_16S)

# Filter taxonomy --------------------------------------------------------------

# ITS taxonomy -----------------------------------------------------------------
head(taxonomy_ITS)
table(taxonomy_ITS$High_level_taxonomy)
table(taxonomy_ITS$Kingdom)

# Based on the High_level_taxonomy I can set the untarget
untarget_taxa <- c(
  "Amoebozoa", "Choanoflagellozoa", "Eukaryota_kgd_Incertae_sedis",
  "Ichthyosporia", "Metazoa", "Protista", "Rhizaria", "Viridiplantae", ""
)

unique(
  c(
    unique(taxonomy_ITS$High_level_taxonomy),
    unique(taxonomy_ITS$Kingdom)
  )
)


apply(taxonomy_ITS, 2, function(x) which(x %in% untarget_taxa))

# Remove non-target taxa and non-classified taxa that did not have any hit
# when blasted against UNITE at 60% conevrage and identity. ALso, remove
# fungal OTUs that were classified as Fungi bu had low hit coverage < 60%
# and identity < 60%

bad_otu <-
  subset(
    taxonomy_ITS,
    High_level_taxonomy %in% untarget_taxa |
      Kingdom %in% untarget_taxa |
      High_level_taxonomy %in% c("") & Kingdom %in% c("") |
      Kingdom %in% c("Fungi") & HL_hit_query_cover < 60 |
      Kingdom %in% c("Fungi") & HL_hit_percent_id < 60
  ) %>%
  filter(!(Genus == "" & Kingdom == "Fungi_1" & Phylum != "" |
    Genus == "" & Kingdom == "Fungi_1" & Phylum != "" |
    Genus == "" & Kingdom == "Fungi_1" & Phylum != "" |
    Genus == "" & Kingdom == "Fungi_1" & Phylum != "")) %>%
  filter(!(Kingdom == "Fungi_1" & Phylum != "" & Class != "" &
    Order != "" & Family != "" & Genus != "")) %>% # save some potentially real fungi
  as.data.frame()

dim(bad_otu)

write.csv(x = bad_otu, 
          file = file.path(results_path, "uparse_bad_taxa.csv"))

# Re-format taxonomy table
ReformatTaxonomy <- function(taxonomy_tab, range_col) {
  require(tidyverse)

  taxonomy_tab <-
    taxonomy_tab %>% mutate_all(na_if, "")
  lastValue <- function(x) tail(x[!is.na(x)], 1)
  last_taxons <- apply(taxonomy_tab[, 1:range_col], 1, lastValue)
  taxonomy_tab$BestMatch <- last_taxons
  taxonomy_res <-
    taxonomy_tab %>%
    unite(OTU_ID, BestMatch, col = Taxonomy, sep = " ", remove = FALSE)
  return(taxonomy_res)
}

taxonomy_ITS_filt <-
  taxonomy_ITS %>%
  rownames_to_column("OTU_ID") %>%
  filter(!OTU_ID %in% rownames(bad_otu)) %>%
  mutate(across(c(
    "Kingdom", "Phylum", "Class", "Order",
    "Family", "Genus", "Species"
  ), ~ gsub("_1", "", x = .))) %>%
  dplyr::select(
    "OTU_ID", "Kingdom", "Phylum", "Class", "Order",
    "Family", "Genus", "Species"
  ) %>%
  as.data.frame() %>%
  ReformatTaxonomy(taxonomy_tab = ., range_col = 8) %>%
  column_to_rownames("OTU_ID")

dim(taxonomy_ITS_filt)
head(taxonomy_ITS_filt)


# 16S taxonomy -----------------------------------------------------------------
unwanted_amp <- c(
  "Mitochondria", "Chloroplast", "mitochondria", "chloroplast",
  "Chloroplast_1", "Mitochondria_1"
)
apply(taxonomy_16S, 2, function(x) which(x %in% unwanted_amp))

head(taxonomy_16S)

taxonomy_16S_filt <-
  taxonomy_16S %>%
  rename(
    Kingdom = Rank_1, Phylum = Rank_2, Class = Rank_3,
    Order = Rank_4, Family = Rank_5, Genus = Rank_6, Species = Rank_7
  ) %>%
  mutate(across(c(
    "Kingdom", "Phylum", "Class", "Order",
    "Family", "Genus", "Species"
  ), ~ gsub("_1", "", x = .))) %>%
  filter(!Order %in% "Chloroplast") %>%
  filter(!Family %in% "Mitochondria") %>%
  rownames_to_column("OTU_ID") %>%
  dplyr::select(
    "OTU_ID", "Kingdom", "Phylum", "Class", "Order",
    "Family", "Genus", "Species"
  ) %>%
  as.data.frame() %>%
  ReformatTaxonomy(taxonomy_tab = ., range_col = 8) %>%
  column_to_rownames("OTU_ID") %>%
  dplyr::select(-BestMatch) %>%
  separate(Taxonomy, into = c("OTU_ID", "BestMatch"), sep = " ", remove = FALSE)

head(taxonomy_16S_filt)


# Create the phyloseq objects --------------------------------------------------

# ITS phyloseq object
physeq_ITS <-
  phyloseq::phyloseq(
    otu_table(otutable_ITS, taxa_are_rows = TRUE),
    sample_data(metadata_ITS),
    tax_table(as.matrix(taxonomy_ITS_filt)),
    sequences_ITS
  )

physeq_ITS

table(sample_data(physeq_ITS)$Treatment)
table(sample_data(physeq_ITS)$Compartment)

# 16S phyloseq object
physeq_16S <-
  phyloseq::phyloseq(
    otu_table(otutable_16S, taxa_are_rows = TRUE),
    sample_data(metadata_16S),
    tax_table(as.matrix(taxonomy_16S_filt)),
    sequences_16S
  )

physeq_16S
table(sample_data(physeq_16S)$Compartment)


# **********************************************************************--------
# ***** DECONTAMINATION ***** --------------------------------------------------
# decontamination from a phyloseq object ---------------------------------------
physeq_ITS@sam_data

sample_data(physeq_ITS)$is.neg <-
  sample_data(physeq_ITS)$SampleNum == "negative_control"

contam_ITS <-
  decontam::isContaminant(physeq_ITS,
    method = "prevalence",
    neg = "is.neg",
    threshold = 0.5
  )

contam_ITS
table(contam_ITS$contaminant)
contam_ITS %>% filter(contaminant == TRUE)

# function to remove taxa by OTU name
remove_taxa <- function(physeq, badTaxa) {
  allTaxa <- taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

subset(contam_ITS, contaminant %in% c("TRUE"))

# filtering the phyloseq object
# Fungi
physeq_ITS_clean <-
  remove_taxa(
    physeq_ITS,
    rownames(subset(contam_ITS, contaminant %in% c("TRUE")))
  ) %>%
  subset_samples(!is.neg %in% TRUE) %>% # remove the control samples
  prune_taxa(taxa_sums(x = .) > 0, x = .) # make sure there aren't OTUs that are 0

physeq_ITS_clean
physeq_ITS_clean@sam_data

# Bacteria
physeq_16S@sam_data

physeq_16S_clean <-
  physeq_16S %>%
  subset_samples(!SampleNum == "Mock") %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .)

physeq_16S_clean

# **********************************************************************--------
# ***** FINAL PHYLOSEQ OBJECTS ***** -------------------------------------------
physeq_ITS_clean
physeq_16S_clean

# Sequencing results -----------------------------------------------------------
sum(sample_sums(physeq_ITS_clean))
mean(sample_sums(physeq_ITS_clean))
sd(sample_sums(physeq_ITS_clean))

sum(sample_sums(physeq_16S_clean))
mean(sample_sums(physeq_16S_clean))
sd(sample_sums(physeq_16S_clean))

# Github datasets --------------------------------------------------------------

saveRDS(object = physeq_ITS_clean, file = file.path(results_path, "phyloseq_ITS.rds"))
saveRDS(object = physeq_16S_clean, file = file.path(results_path, "phyloseq_16S.rds"))

# **********************************************************************--------
# ***** DATA MINING and PRELIMINARY EXPLORATION ***** --------------------------

# Adding alpha metrics ---------------------------------------------------------

AlphaMetrics <- function(physeq) {
  require(vegan)
  require(tidyverse)

  sample_data(physeq)$ReadNo <- sample_sums(physeq)
  sample_data(physeq)$hill_0 <- as.data.frame(as.matrix(t(physeq@otu_table))) %>% renyi(scale = c(0), hill = TRUE)
  sample_data(physeq)$hill_1 <- as.data.frame(as.matrix(t(physeq@otu_table))) %>% renyi(scale = c(1), hill = TRUE)
  sample_data(physeq)$hill_2 <- as.data.frame(as.matrix(t(physeq@otu_table))) %>% renyi(scale = c(2), hill = TRUE)
  sample_data(physeq)$Richness <- specnumber(as.data.frame(otu_table(physeq)), MARGIN = 2)
  sample_data(physeq)$invSimpson <- diversity(as.data.frame(otu_table(physeq)), index = "inv", MARGIN = 2)
  sample_data(physeq)$Shannon <- diversity(as.data.frame(otu_table(physeq)), index = "shannon", MARGIN = 2)
  sample_data(physeq)$EH <- 1 - sample_data(physeq)$Shannon / log(sample_data(physeq)$Richness)

  return(physeq)
}


physeq_ITS_clean <-
  AlphaMetrics(physeq_ITS_clean)
physeq_ITS_clean

physeq_ITS_clean@sam_data

physeq_16S_clean <-
  AlphaMetrics(physeq_16S_clean)
physeq_16S_clean

physeq_16S_clean@sam_data

# Visualizing read and richness distribution -----------------------------------

df_histo_prelim <-
  rbind(
    physeq_ITS_clean@sam_data %>%
      as.matrix() %>%
      as.data.frame() %>%
      dplyr::select(ReadNo, Richness, Compartment, Treatment) %>%
      rownames_to_column("SampleID") %>%
      mutate(
        SampleID = paste(SampleID, "F", sep = ""),
        Kingdom = "Fungi",
        Kingdom = as.factor(Kingdom)
      ),
    physeq_16S_clean@sam_data %>%
      as.matrix() %>%
      as.data.frame() %>%
      dplyr::select(ReadNo, Richness, Compartment, Treatment) %>%
      rownames_to_column("SampleID") %>%
      mutate(
        SampleID = paste(SampleID, "B", sep = ""),
        Kingdom = "Bacteria",
        Kingdom = as.factor(Kingdom)
      )
  ) %>%
  as_tibble() %>%
  mutate(
    ReadNo = as.numeric(ReadNo),
    Richness = as.numeric(Richness)
  )

plot_readNO_histo <-
  df_histo_prelim %>%
  ggplot(aes(x = ReadNo, fill = Compartment)) +
  geom_histogram(binwidth = 2000, boundary = 0, closed = "left", alpha = 0.6) +
  labs(title = "Distribution of the samples depth", y = "Samples") +
  theme_bw() +
  facet_wrap(~Kingdom, ncol = 2, scales = "free") +
  theme(
    plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5, vjust = 0.5),
    plot.subtitle = element_markdown(size = 10, face = "bold", hjust = 0.5, vjust = 0.5),
    axis.text.x = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
    axis.text.y = element_markdown(size = 8),
    axis.title = element_markdown(size = 10, face = "bold", hjust = 0.5, vjust = 0.5),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.position = "bottom"
  )

plot_Richness_histo <-
  df_histo_prelim %>%
  ggplot(aes(x = Richness, fill = Compartment)) +
  geom_histogram(binwidth = 50, boundary = 0, closed = "left", alpha = 0.6) +
  labs(title = "Distribution of the samples richness", y = "Samples") +
  theme_bw() +
  facet_wrap(~Kingdom, ncol = 2, scales = "free") +
  theme(
    plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5, vjust = 0.5),
    plot.subtitle = element_markdown(size = 10, face = "bold", hjust = 0.5, vjust = 0.5),
    axis.text.x = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
    axis.text.y = element_markdown(size = 8),
    axis.title = element_markdown(size = 10, face = "bold", hjust = 0.5, vjust = 0.5),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.position = "bottom"
  )

# >>>>> FIGURE S2 <<<<< --------------------------------------------------------
# histogram of Read No and Richness across Compartments
Fig_S2 <-
  ggarrange(plot_readNO_histo,
    plot_Richness_histo,
    ncol = 1,
    nrow = 2,
    labels = c("a", "b"),
    common.legend = TRUE,
    legend = "bottom"
  )

Fig_S2

ggsave(plot = Fig_S2, 
       path = results_path,
       filename = "Fig_S2.pdf")

# **********************************************************************--------
# ***** RAREFACTION CURVES ***** -----------------------------------------------
sort(sample_sums(physeq_ITS_clean))

rarecurve_fungi <-
  physeq_ITS_clean@otu_table %>%
  t() %>%
  as.matrix() %>%
  as.data.frame() %>%
  rarecurve(x = ., step = 1000, tidy = TRUE)

head(rarecurve_fungi)

sort(sample_sums(physeq_16S_clean))

rarecurve_bacteria <-
  physeq_16S_clean@otu_table %>%
  t() %>%
  as.matrix() %>%
  as.data.frame() %>%
  rarecurve(x = ., step = 1000, tidy = TRUE)

head(rarecurve_bacteria)


# Then plotting rarecurve in ggplot2
PlotRareCurve <- function(rare_tidy, depth = NULL, Col) {
  plot_rare <-
    ggplot(rare_tidy, aes(x = Sample, y = Species, group = SampleID, color = get(Col))) +
    geom_line() +
    geom_vline(xintercept = depth, color = "black", linetype = "dashed") +
    theme_bw() +
    theme(
      plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5, vjust = 0.5),
      plot.subtitle = element_markdown(size = 10, face = "bold", hjust = 0.5, vjust = 0.5),
      axis.text.x = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
      axis.text.y = element_markdown(size = 8),
      axis.title = element_markdown(size = 10, face = "bold", hjust = 0.5, vjust = 0.5),
      legend.key.height = unit(0.5, "cm"),
      legend.key.width = unit(0.5, "cm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.position = "right"
    )
  return(plot_rare)
}

# Setting up desired depth -----------------------------------------------------
fungi_depth <- 19849
bacteria_depth <- 31244

# >>>>> FIGURE S3 <<<<< --------------------------------------------------------

Fig_S3_rarecurve <-
  ggarrange(
    rarecurve_fungi %>%
      rename("SampleID" = Site) %>%
      left_join(.,
        physeq_ITS_clean@sam_data %>%
          as.matrix() %>%
          as.data.frame() %>%
          rownames_to_column("SampleID"),
        by = "SampleID"
      ) %>%
      PlotRareCurve(., fungi_depth, "Compartment") +
      labs(title = "ITS", x = "Number of DNA reads", y = "Number of OTUs") +
      ggplot2::annotate("text",
        x = 13000, y = 460, label = "19849",
        fontface = "bold", hjust = 0.5, vjust = 0.5,
        col = "black"
      ),
    rarecurve_bacteria %>%
      rename("SampleID" = Site) %>%
      left_join(.,
        physeq_16S_clean@sam_data %>%
          as.matrix() %>%
          as.data.frame() %>%
          rownames_to_column("SampleID"),
        by = "SampleID"
      ) %>%
      PlotRareCurve(., bacteria_depth, "Compartment") +
      labs(title = "16S", x = "Number of DNA reads", y = "Number of OTUs") +
      ggplot2::annotate("text",
        x = 24000, y = 3800, label = "31244",
        fontface = "bold", hjust = 0.5, vjust = 0.5,
        col = "black"
      ),
    ncol = 1,
    nrow = 2,
    labels = c("a", "b"),
    common.legend = TRUE,
    align = "hv",
    legend = "bottom"
  )

Fig_S3_rarecurve

ggsave(plot = Fig_S3_rarecurve, 
       path = results_path,
       filename = "Fig_S3_rarefaction_curves.pdf")

# **********************************************************************--------
# ***** OPTIMAL RAREFACTION DEPTH ***** ----------------------------------------

# Calculating Good's coverage. The fraction of sequences that appear in an OTU
# that have been seen more than once, and allows estimating what percent of the total
# species is represented in a sample.
# Coverage = 1 - (number of individuals in species / total number of individuals)
# Example: If I have Goods = 0.96 it means that 4% of your reads in that sample are
# from OTUs that appear only once in that sample.

# Calculate long dataframe with stats
RareStats <- function(physeq) {
  require(tidyverse)

  # calculate distribution outliers
  findoutlier <- function(x) {
    return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
  }

  # generate output dataframe
  df_output <-
    as.data.frame(as.matrix(physeq@otu_table)) %>%
    rownames_to_column("OTU_ID") %>%
    pivot_longer(-OTU_ID, names_to = "SampleID", values_to = "Seq_No") %>%
    group_by(SampleID) %>%
    summarize(
      ReadNo = sum(Seq_No),
      Singlton_No = sum(Seq_No == 1),
      Goods = 100 * (1 - Singlton_No / ReadNo)
    ) %>%
    mutate(outlier = ifelse(findoutlier(log10(ReadNo)), ReadNo, NA)) %>%
    ungroup()

  return(df_output)
}


# Plotting function
PloRareStats <- function(dataframe) {
  require(ggrepel)
  plot_output <-
    ggarrange(
      dataframe %>%
        ggplot(aes(x = ReadNo)) +
        geom_histogram(
          binwidth = 5000,
          fill = "firebrick", color = "firebrick"
        ) +
        theme_bw() +
        theme(axis.text.x = element_markdown(angle = 33, hjust = 1, vjust = 1)) +
        labs(title = "Histogram"),
      dataframe %>%
        ggplot(aes(x = ReadNo)) +
        geom_histogram(
          binwidth = 1000,
          fill = "firebrick", color = "firebrick"
        ) +
        coord_cartesian(xlim = c(0, 25000)) +
        theme_bw() +
        theme(axis.text.x = element_markdown(angle = 33, hjust = 1, vjust = 1)) +
        labs(title = "Histogram Zoom"),
      dataframe %>%
        ggplot(aes(x = ReadNo, y = Goods)) +
        geom_point(shape = 1, color = "firebrick") +
        theme_bw() +
        theme(axis.text.x = element_markdown(angle = 33, hjust = 1, vjust = 1)) +
        labs(title = "Good's Coverage"),
      dataframe %>%
        ggplot(aes(x = 1, y = ReadNo)) +
        geom_jitter(shape = 1, color = "firebrick") +
        theme_bw() +
        theme(axis.text.x = element_markdown(angle = 33, hjust = 1, vjust = 1)) +
        scale_y_log10() +
        labs(title = "Log10 jitter"),
      dataframe %>%
        ggplot(aes(x = 1, y = ReadNo)) +
        geom_boxplot(color = "firebrick") +
        theme_bw() +
        theme(axis.text.x = element_markdown(angle = 33, hjust = 1, vjust = 1)) +
        scale_y_log10() +
        geom_text_repel(
          data = filter(dataframe, !is.na(outlier)),
          mapping = aes(x = 1, y = ReadNo, label = outlier),
          max.overlaps = 15, size = 3
        ) +
        labs(title = "Log10 boxplot"),
      dataframe %>%
        arrange(ReadNo) %>%
        ggplot(aes(x = 1:nrow(.), y = ReadNo)) +
        geom_bar(stat = "identity", color = "firebrick") +
        theme_bw() +
        theme(axis.text.x = element_markdown(angle = 33, hjust = 1, vjust = 1)) +
        labs(title = "Ranked"),
      ncol = 3,
      nrow = 2,
      align = "hv",
      labels = c("a", "b", "c", "d", "e", "f")
    )
  return(plot_output)
}


RareStats(physeq_ITS_clean)

rare_fungi <-
  RareStats(physeq_ITS_clean) %>%
  arrange(ReadNo) %>%
  left_join(
    physeq_ITS_clean@sam_data %>%
      as.matrix() %>%
      as.data.frame() %>%
      rownames_to_column("SampleID") %>%
      dplyr::select(
        SampleID, Compartment, Treatment,
        hill_0, hill_1, hill_2,
        Richness, invSimpson, Shannon
      ),
    by = "SampleID"
  )

rare_fungi

rare_bact <-
  RareStats(physeq_16S_clean) %>%
  arrange(ReadNo) %>%
  left_join(
    physeq_16S_clean@sam_data %>%
      as.matrix() %>%
      as.data.frame() %>%
      rownames_to_column("SampleID") %>%
      dplyr::select(
        SampleID, Compartment, Treatment,
        hill_0, hill_1, hill_2,
        Richness, invSimpson, Shannon
      ),
    by = "SampleID"
  )
rare_bact

# plotting rarestats -----------------------------------------------------------
# >>>>> FIGURE S4 <<<<< --------------------------------------------------------
title1 <- text_grob("Fungi", size = 12, face = 2)

Fig_S4_rare_fungi <- rare_fungi %>% PloRareStats()
grid.arrange(Fig_S4_rare_fungi, top = title1)

ggsave(plot = grid.arrange(Fig_S4_rare_fungi, top = title1), 
       path = results_path,
       filename = "Fig_S4_rarefaction_depth_fungi.pdf")

# >>>>> FIGURE S4 <<<<< --------------------------------------------------------
title2 <- text_grob("Bacteria", size = 12, face = 2)

Fig_S5_rare_bact <- rare_bact %>% PloRareStats()
grid.arrange(Fig_S5_rare_bact, top = title2)

ggsave(plot = grid.arrange(Fig_S5_rare_bact, top = title1), 
       path = results_path,
       filename = "Fig_S5_rarefaction_depth_bacteria.pdf")

# **********************************************************************--------
# ***** RAREFACTION ***** ------------------------------------------------------
rarefyData <- function(physeq, depth_level) {
  require(tidyverse)

  dataframe <-
    as.data.frame(as.matrix(t(physeq@otu_table)))

  com_iter <- vector(mode = "list", length = 100)

  for (i in seq_along(com_iter)) {
    com_iter[[i]] <- as.data.frame(
      vegan::rrarefy(dataframe, sample = depth_level)
    ) %>% rownames_to_column("SampleID")
  }

  mean_100 <- do.call(rbind, com_iter)
  mean_100 <- mean_100 %>%
    group_by(SampleID) %>%
    summarise(across(everything(), mean)) %>%
    filter(rowSums(across(where(is.numeric))) >= depth_level)

  print(mean_100 %>% as_tibble())
  return(mean_100)
}


# Fungi ------------------------------------------------------------------------

# recreating the phyloseq objects with the rarefied otus
set.seed(2332)

# test
rarefyData(physeq_ITS_clean, fungi_depth)

# add to the fungal phyloseq object
physeq_ITS_rare <-
  phyloseq(
    otu_table(rarefyData(physeq_ITS_clean, fungi_depth) %>%
      column_to_rownames("SampleID") %>%
      t() %>%
      as.matrix() %>%
      as.data.frame(), taxa_are_rows = TRUE),
    physeq_ITS_clean@sam_data,
    physeq_ITS_clean@tax_table,
    physeq_ITS_clean@refseq
  ) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_ITS_rare


# Bacteria ---------------------------------------------------------------------

# NOTE This require rarefying for each niche separately which is fine since we do
# not aim to compare between niches. Please see M&M for details.

physeq_16S_roots <-
  physeq_16S_clean %>%
  subset_samples(Compartment %in% c("Root")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_16S_roots
physeq_16S_roots@sam_data

# Soil and rhizosphere together
physeq_16S_soilrhizo <-
  physeq_16S_clean %>%
  subset_samples(Compartment %in% c("soil", "Rhizosphere")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_16S_soilrhizo
physeq_16S_soilrhizo@sam_data

# Rarefying Root separately form Soil and Rhizosphere --------------------------

# Quick stats before deciding depths -------------------------------------------
RareStats(physeq_16S_roots) %>%
  arrange(ReadNo) %>%
  PloRareStats()

RareStats(physeq_16S_soilrhizo) %>%
  arrange(ReadNo) %>%
  PloRareStats()

bacteria_depth_root <- 440
bacteria_depth_soilrhizo <- 31244

# rarefying --------------------------------------------------------------------
physeq_16S_roots_rare <-
  phyloseq(
    otu_table(rarefyData(physeq_16S_roots, bacteria_depth_root) %>%
      column_to_rownames("SampleID") %>%
      t() %>%
      as.matrix() %>%
      as.data.frame(), taxa_are_rows = TRUE),
    physeq_16S_roots@sam_data,
    physeq_16S_roots@tax_table,
    physeq_16S_roots@refseq
  ) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_16S_roots_rare


physeq_16S_soilrhizo_rare <-
  phyloseq(
    otu_table(rarefyData(physeq_16S_soilrhizo, bacteria_depth_soilrhizo) %>%
      column_to_rownames("SampleID") %>%
      t() %>%
      as.matrix() %>%
      as.data.frame(), taxa_are_rows = TRUE),
    physeq_16S_soilrhizo@sam_data,
    physeq_16S_soilrhizo@tax_table,
    physeq_16S_soilrhizo@refseq
  ) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_16S_soilrhizo_rare

# Recreating the phyloseq object (for convenience) with roots, rhizosphere and
# soil but with niches rarefied at different depth.

# merging tables first
combined_otutab_rare <-
  full_join(
    as.data.frame(physeq_16S_roots_rare@otu_table) %>%
      rownames_to_column("otuID"),
    as.data.frame(physeq_16S_soilrhizo_rare@otu_table) %>%
      rownames_to_column("otuID"),
    by = "otuID"
  ) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames("otuID")

dim(combined_otutab_rare)
colSums(combined_otutab_rare)

# visual check!
write.csv(x = combined_otutab_rare, file = file.path(results_path, "combined_otutab_rare.csv"))

# merging
physeq_16S_rare <-
  phyloseq(
    otu_table(combined_otutab_rare, taxa_are_rows = TRUE),
    physeq_16S_clean@sam_data,
    physeq_16S_clean@tax_table,
    physeq_16S_clean@refseq
  )

physeq_16S_rare

# Quick check after merging
sample_sums(physeq_16S_rare)

RareStats(physeq_16S_rare) %>%
  arrange(ReadNo) %>%
  PloRareStats()

# **********************************************************************--------
# ***** BETA DIVERSITY ***** ---------------------------------------------------

sum(sample_sums(physeq_ITS_rare))
mean(sample_sums(physeq_ITS_rare))
sd(sample_sums(physeq_ITS_rare))

sum(sample_sums(physeq_16S_clean))
mean(sample_sums(physeq_16S_clean))
sd(sample_sums(physeq_16S_clean)) # confirms different rarefaction depths!

# Plotting Fungi ---------------------------------------------------------------
table(physeq_ITS_rare@sam_data$Compartment, physeq_ITS_rare@sam_data$Treatment)

# Separating dataset according to research questions and experimental setup.

physeq_ITS_rare_root <-
  subset_samples(physeq_ITS_rare, Compartment %in% c("Root")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_ITS_rare_rhizo <-
  subset_samples(physeq_ITS_rare, Compartment %in% c("Rhizosphere")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_ITS_rare_soil <-
  subset_samples(physeq_ITS_rare, Compartment %in% c("soil") &
    Treatment %in% c("IE", "NE", "IE_NE_MIX")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_ITS_rare_soil@sam_data

physeq_ITS_rare_control <-
  subset_samples(physeq_ITS_rare, Treatment %in% c("control", "fumigated")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_ITS_rare_control@sam_data

# function for plotting beta diversity
plot_pcoa <- function(ord, phy, Col=NULL, She=NULL){
  
  require(phyloseq)
  
  plot_res <-
    plot_ordination(physeq = phy, 
                    ordination=ord,  
                    type="samples", 
                    color=Col, 
                    shape=She) +
    theme_bw() +
    theme(plot.title = element_markdown(size = 12, face = "bold", hjust = 0.5),
          strip.text = element_markdown(size = 10, face = "bold"),
          axis.text.x = element_markdown(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
          axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
          strip.background = element_blank(),
          legend.key.height = unit(0.4, "cm"), 
          legend.key.width = unit(0.4, "cm"),
          legend.position = "right", 
          legend.title.align = 0.5,
          legend.text = element_markdown(size = 8)) +
    guides(color = guide_legend(ncol=3, override.aes = list(shape = 15, size = 3.5)),
           shape = guide_legend(ncol=1, override.aes = list(color = "black", size=2.5)),
           size = guide_legend(ncol=1))
  
  return(plot_res)
  
}

# plotting
title3 <- text_grob("Bray-Curtis PCoA Fungi", size = 12, face = 2)

fungi_beta <-
  ggarrange(
    ordinate(physeq_ITS_rare_root, method = "PCoA", distance = "bray") %>%
      plot_pcoa(phy = physeq_ITS_rare_root, ord = ., Col = "Treatment", NULL) +
      geom_point(size = 1) +
      scale_color_manual(values = palette2) +
      theme_bw() +
      theme(
        plot.title = element_markdown(size = 8, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title = element_markdown(size = 7),
        axis.text.x = element_markdown(size = 7, colour = "black", hjust = 1, vjust = 1),
        axis.text.y = element_markdown(size = 7, angle = 0, hjust = 0.5),
        legend.position = "none"
      ) +
      labs(title = "Root"),
    ordinate(physeq_ITS_rare_rhizo, method = "PCoA", distance = "bray") %>%
      plot_pcoa(phy = physeq_ITS_rare_rhizo, ord = ., Col = "Treatment", NULL) +
      geom_point(size = 1) +
      scale_color_manual(values = palette2) +
      theme_bw() +
      theme(
        plot.title = element_markdown(size = 8, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title = element_markdown(size = 7),
        axis.text.x = element_markdown(size = 7, colour = "black", hjust = 1, vjust = 1),
        axis.text.y = element_markdown(size = 7, angle = 0, hjust = 0.5),
        legend.position = "none"
      ) +
      labs(title = "Rhizosphere"),
    ordinate(physeq_ITS_rare_soil, method = "PCoA", distance = "bray") %>%
      plot_pcoa(phy = physeq_ITS_rare_soil, ord = ., Col = "Treatment", NULL) +
      geom_point(size = 1) +
      scale_color_manual(values = c("#825121", "#CC79A7", "#FFB400")) +
      theme_bw() +
      theme(
        plot.title = element_markdown(size = 8, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title = element_markdown(size = 7),
        axis.text.x = element_markdown(size = 7, colour = "black", hjust = 1, vjust = 1),
        axis.text.y = element_markdown(size = 7, angle = 0, hjust = 0.5),
        legend.position = "none"
      ) +
      labs(title = "Soil"),
    ordinate(physeq_ITS_rare_control, method = "PCoA", distance = "bray") %>%
      plot_pcoa(phy = physeq_ITS_rare_control, ord = ., Col = "Treatment", NULL) +
      geom_point(size = 1) +
      scale_color_manual(values = c("#7FB800", "#2D3142")) +
      theme_bw() +
      theme(
        plot.title = element_markdown(size = 8, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title = element_markdown(size = 7),
        axis.text.x = element_markdown(size = 7, colour = "black", hjust = 1, vjust = 1),
        axis.text.y = element_markdown(size = 7, angle = 0, hjust = 0.5),
        legend.position = "none"
      ) +
      labs(title = "Control"),
    ncol = 4,
    nrow = 1,
    align = "hv",
    labels = c("A", "B", "C", "D"),
    common.legend = FALSE,
    legend = "none"
  )

fungi_beta

# Plotting Bacteria ------------------------------------------------------------
table(physeq_16S_rare@sam_data$Compartment, physeq_16S_rare@sam_data$Treatment)

# Fix the rarefied taxonomy
physeq_16S_rare

modify_taxonomy <- function(physeq) {
  tax_table <-
    tax_table(physeq) %>%
    as.matrix() %>%
    as.data.frame() %>%
    mutate(
      OTU_ID = word(Taxonomy, 1), # Extract the first word as OTU ID
      BestMatch = word(Taxonomy, 2, -1) %>%
        str_replace_all("_", " ") # Remove underscores in classification name
    ) %>%
    dplyr::select(OTU_ID, BestMatch, everything()) %>% # Reorder columns
    mutate(Taxonomy = paste(OTU_ID, BestMatch, sep = " ")) %>%
    dplyr::select(-OTU_ID)

  return(tax_table)
}

# check
tax_table(as.matrix(modify_taxonomy(physeq_16S_rare))) %>% head()

tax_table(physeq_16S_rare) <-
  tax_table(as.matrix(modify_taxonomy(physeq_16S_rare)))


physeq_16S_rare_root <-
  subset_samples(physeq_16S_rare, Compartment %in% c("Root")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_16S_rare_rhizo <-
  subset_samples(physeq_16S_rare, Compartment %in% c("Rhizosphere")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_16S_rare_soil <-
  subset_samples(physeq_16S_rare, Compartment %in% c("soil") &
    Treatment %in% c("IE", "NE", "IE_NE MIX")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_16S_rare_soil@sam_data

physeq_16S_rare_control <-
  subset_samples(physeq_16S_rare, Treatment %in% c("control", "fumigated")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)


title4 <- text_grob("Bray-Curtis PCoA Bacteria", size = 12, face = 2)

bacteria_beta <-
  ggarrange(
    ordinate(physeq_16S_rare_root, method = "PCoA", distance = "bray") %>%
      plot_pcoa(phy = physeq_16S_rare_root, ord = ., Col = "Treatment", NULL) +
      geom_point(size = 1) +
      scale_color_manual(values = palette2) +
      theme_bw() +
      theme(
        plot.title = element_markdown(size = 8, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title = element_markdown(size = 7),
        axis.text.x = element_markdown(size = 7, colour = "black", hjust = 1, vjust = 1),
        axis.text.y = element_markdown(size = 7, angle = 0, hjust = 0.5),
        legend.position = "none"
      ) +
      labs(title = "Root"),
    ordinate(physeq_16S_rare_rhizo, method = "PCoA", distance = "bray") %>%
      plot_pcoa(phy = physeq_16S_rare_rhizo, ord = ., Col = "Treatment", NULL) +
      geom_point(size = 1) +
      scale_color_manual(values = palette2) +
      theme_bw() +
      theme(
        plot.title = element_markdown(size = 8, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title = element_markdown(size = 7),
        axis.text.x = element_markdown(size = 7, colour = "black", hjust = 1, vjust = 1),
        axis.text.y = element_markdown(size = 7, angle = 0, hjust = 0.5),
        legend.position = "none"
      ) +
      labs(title = "Rhizosphere"),
    ordinate(physeq_16S_rare_soil, method = "PCoA", distance = "bray") %>%
      plot_pcoa(phy = physeq_16S_rare_soil, ord = ., Col = "Treatment", NULL) +
      geom_point(size = 1) +
      scale_color_manual(values = c("#825121", "#CC79A7", "#FFB400")) +
      theme_bw() +
      theme(
        plot.title = element_markdown(size = 8, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title = element_markdown(size = 7),
        axis.text.x = element_markdown(size = 7, colour = "black", hjust = 1, vjust = 1),
        axis.text.y = element_markdown(size = 7, angle = 0, hjust = 0.5),
        legend.position = "none"
      ) +
      labs(title = "Soil"),
    ordinate(physeq_16S_rare_control, method = "PCoA", distance = "bray") %>%
      plot_pcoa(phy = physeq_16S_rare_control, ord = ., Col = "Treatment", NULL) +
      geom_point(size = 1) +
      scale_color_manual(values = c("#7FB800", "#2D3142")) +
      theme_bw() +
      theme(
        plot.title = element_markdown(size = 8, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title = element_markdown(size = 7),
        axis.text.x = element_markdown(size = 7, colour = "black", hjust = 1, vjust = 1),
        axis.text.y = element_markdown(size = 7, angle = 0, hjust = 0.5),
        legend.position = "none"
      ) +
      labs(title = "Control"),
    ncol = 4,
    nrow = 1,
    align = "hv",
    labels = c("E", "F", "G", "H"),
    common.legend = FALSE,
    legend = "none"
  )


bacteria_beta

# Make the legend first
MakeLegend <- function() {
  data_legend <-
    data.frame(
      label = c(
        "IE", "NE", "IE_MIX", "NE_MIX",
        "IE_NE_MIX", "Fumigated", "Control"
      ),
      color = c("#825121", "#FFB400", "#CC2D35", "#00A6ED", "#CC79A7", "#2D3142", "#7FB800")
    ) %>%
    mutate(label = fct_relevel(
      label,
      "IE", "NE", "IE_MIX", "NE_MIX",
      "IE_NE_MIX", "Fumigated", "Control"
    ))

  legend_obj <-
    ggplot(data_legend, aes(x = 1, y = label, color = label)) +
    geom_point() +
    scale_color_manual(values = data_legend$color) +
    theme(
      legend.position = "right",
      legend.key.height = unit(0.2, "cm"),
      legend.key.width = unit(0.2, "cm")
    ) +
    guides(color = guide_legend(
      order = 1, ncol = 7, title = NULL,
      override.aes = list(shape = 15, size = 3)
    ))

  return(get_legend(legend_obj))
}

MakeLegend() %>% as_ggplot()


#  >>>>> FIGURE 3 <<<<< --------------------------------------------------------
Fig_3_betadiv <- ggarrange(
  ggarrange(
    grid.arrange(fungi_beta, top = text_grob("Bray-Curtis PCoA Fungi", size = 12, face = 2)),
    grid.arrange(bacteria_beta, top = text_grob("Bray-Curtis PCoA Bacteria", size = 11, face = 2)),
    ncol = 1,
    nrow = 2
  ),
  MakeLegend() %>% as_ggplot(),
  ncol = 1,
  nrow = 2,
  heights = c(0.9, 0.1)
)

Fig_3_betadiv

ggsave(plot = Fig_3_betadiv, 
       path = results_path,
       filename = "Fig_3_beta_diversity.pdf")

# **********************************************************************--------
# ***** PERMANOVA ***** --------------------------------------------------------

Adonis2All <- function(physeq, method, strata = NULL, perm = 999) {
  require(tidyverse)
  require(vegan)

  # Extract OTU and metadata
  otu <- physeq@otu_table %>%
    t() %>%
    as.matrix() %>%
    as.data.frame()
  metadata <- physeq@sam_data %>%
    as.matrix() %>%
    as.data.frame()

  # Run analysis
  df_adonis <- adonis2(
    formula(otu ~ Treatment),
    metadata,
    method = method,
    strata = NULL,
    permutations = perm,
    parallel = 8
  )

  # Adjust p-values
  df_adj <- cbind(
    df_adonis,
    as.data.frame(p.adjust(df_adonis$`Pr(>F)`, method = "BH")) %>%
      rename("p.adj" = 1)
  )

  # Return results
  return(list(
    formula = formula(otu ~ Treatment),
    df_adonis = df_adonis,
    df_adj = df_adj
  ))
}

# checks
Adonis2All(physeq_ITS_rare_root, "bray")[[3]]
Adonis2All(physeq_ITS_rare_root, "bray")

# running PRMANOVA -------------------------------------------------------------
set.seed(030225)

adonis_result <-
  data.frame(
    rbind(
      Adonis2All(physeq_ITS_rare_root, "bray")[[3]],
      Adonis2All(physeq_ITS_rare_rhizo, "bray")[[3]],
      Adonis2All(physeq_ITS_rare_soil, "bray")[[3]],
      Adonis2All(physeq_ITS_rare_control, "bray")[[3]],
      Adonis2All(physeq_16S_rare_root, "bray")[[3]],
      Adonis2All(physeq_16S_rare_rhizo, "bray")[[3]],
      Adonis2All(physeq_16S_rare_soil, "bray")[[3]],
      Adonis2All(physeq_16S_rare_control, "bray")[[3]]
    ),
    Kingdom = c(rep("Fungi", 12), rep("Bacteria", 12)),
    Dataset = c(Dataset = c(
      rep("Root", 3),
      rep("Rhizosphere", 3),
      rep("Soil", 3),
      rep("Control", 3),
      rep("Root", 3),
      rep("Rhizosphere", 3),
      rep("Soil", 3),
      rep("Control", 3)
    ))
  ) %>%
  rownames_to_column("Group") %>%
  mutate(Group = rep(c("Treatment", "Residuals", "Total"), 8))

adonis_result

#  >>>>> TABLE 1 <<<<< --------------------------------------------------------
write.csv(x = adonis_result, 
          file = file.path(results_path,"Table1_adonis_result.csv"))

# **********************************************************************--------
# ***** BETA-DISPERSION ***** --------------------------------------------------
BetadispExtr <- function(physeq, method, Var) {
  require(tidyverse)
  require(vegan)

  # Interesting fact! within a function it is better to specify the otu
  # generated from a phyloeq object in multiple steps otherwise
  # does not really work if you include it in ().
  # e.g. otu <- as.data.frame(as.matrix(t(otu_table(physeq))))
  otu <-
    physeq@otu_table %>%
    t() %>%
    as.matrix() %>%
    as.data.frame()

  metadata <-
    physeq@sam_data %>%
    as.matrix() %>%
    as.data.frame()

  disp <-
    betadisper(
      vegan::vegdist(otu, method = method),
      metadata[, Var]
    )
  anova_d <-
    anova(disp,
      permutations = how(nperm = 999)
    )
  p_adj <-
    round(p.adjust(
      anova_d$`Pr(>F)`,
      "BH"
    ), 4)
  dist_var <-
    vegan::permutest(disp,
      permutations = 999,
      pairwise = T
    )

  # Extract pairwise comparisons into a data frame, if present
  pairwise_df <- NULL
  pw <- dist_var$pairwise

  if (!is.null(pw)) {
    # Often $observed and $permuted are named numeric vectors
    if (is.numeric(pw$observed) & is.numeric(pw$permuted)) {
      pairwise_df <- data.frame(
        Comparison = names(pw$observed),
        Observed   = as.numeric(pw$observed),
        Permuted   = as.numeric(pw$permuted)
      ) %>%
        # If you want to split "IE-NE_MIX" into two columns "Group1" and "Group2":
        tidyr::separate(col = "Comparison", into = c("group1", "group2"), sep = "-") %>%
        mutate(Permuted_adj = round(p.adjust(Permuted, method = "BH"), digits = 3))
    }
  }

  # 6) Return a list with all relevant objects
  return(list(
    dist_var       = dist_var, # The permutest object
    p_adj_anova    = p_adj, # Adjusted p-value from the global ANOVA
    betadisp_model = disp, # The betadisper object
    pairwise_df    = pairwise_df # Tidy data frame of pairwise comparisons
  ))
}


# checks
res <- BetadispExtr(physeq_ITS_rare_root, "bray", "Treatment")
str(res)

res$pairwise_df
res$dist_var

BetadispExtr(physeq_ITS_rare_root, "bray", "Treatment")[[1]]$tab

# running betadisper ------------------------------------------------------------
multipleBetadisper <- function() {
  require(vegan)

  betadisp_res <-
    data.frame(
      rbind(
        BetadispExtr(physeq_ITS_rare_root, "bray", "Treatment")[[1]]$tab,
        BetadispExtr(physeq_ITS_rare_rhizo, "bray", "Treatment")[[1]]$tab,
        BetadispExtr(physeq_ITS_rare_soil, "bray", "Treatment")[[1]]$tab,
        BetadispExtr(physeq_ITS_rare_control, "bray", "Treatment")[[1]]$tab,
        BetadispExtr(physeq_16S_rare_root, "bray", "Treatment")[[1]]$tab,
        BetadispExtr(physeq_16S_rare_rhizo, "bray", "Treatment")[[1]]$tab,
        BetadispExtr(physeq_16S_rare_soil, "bray", "Treatment")[[1]]$tab,
        BetadispExtr(physeq_16S_rare_control, "bray", "Treatment")[[1]]$tab
      ),
      Padj = c(
        BetadispExtr(physeq_ITS_rare_root, "bray", "Treatment")[[2]],
        BetadispExtr(physeq_ITS_rare_rhizo, "bray", "Treatment")[[2]],
        BetadispExtr(physeq_ITS_rare_soil, "bray", "Treatment")[[2]],
        BetadispExtr(physeq_ITS_rare_soil, "bray", "Treatment")[[2]],
        BetadispExtr(physeq_16S_rare_control, "bray", "Treatment")[[2]],
        BetadispExtr(physeq_16S_rare_rhizo, "bray", "Treatment")[[2]],
        BetadispExtr(physeq_16S_rare_soil, "bray", "Treatment")[[2]],
        BetadispExtr(physeq_16S_rare_control, "bray", "Treatment")[[2]]
      ),
      Kigdom = c(rep("Fungi", 8), rep("Bacteria", 8)),
      Dataset = c(
        "Root", "Root",
        "Rhizosphere", "Rhizosphere",
        "Soil", "Soil",
        "Control", "Control",
        "Root", "Root",
        "Rhizosphere", "Rhizosphere",
        "Soil", "Soil",
        "Control", "Control"
      )
    )

  return(betadisp_res)
}

betadisper_result <-
  multipleBetadisper() %>%
  rownames_to_column("Group") %>%
  mutate(Group = rep(c("Treatment", "Residuals"), 8))

betadisper_result

write.csv(x = betadisper_result,
          file = file.path(results_path, "Table1_betadisper_result.csv"))

# **********************************************************************--------
# ***** PAIRWISE PERMANOVA ***** -----------------------------------------------
pairwise_permanova <- function(physeq,
                               Var,
                               dist = "bray",
                               adj = "BH",
                               perm = 999) {
  require(vegan)
  require(tidyverse)

  sp_matrix <- as.data.frame(t(as.matrix(otu_table(physeq,taxa_are_rows = TRUE))))
  metadata <- as.data.frame(as.matrix(sample_data(physeq)))

  ## list contrasts
  group_var <- metadata %>% pull(Var)
  groups <- as.data.frame(t(combn(unique(group_var), m = 2)))

  contrasts <- data.frame(
    group1 = groups$V1, group2 = groups$V2,
    R2 = NA, F_value = NA, df1 = NA, df2 = NA, p_value = NA
  )

  for (i in seq(nrow(contrasts))) {
    sp_subset <- group_var == contrasts$group1[i] | group_var == contrasts$group2[i]
    contrast_matrix <- sp_matrix[sp_subset, ]

    ## fit contrast using adonis
    fit <- vegan::adonis2(
      contrast_matrix ~ group_var[sp_subset],
      method = dist,
      perm = perm,
      parallel = 8
    )

    contrasts$R2[i] <- round(fit$R2[1], digits = 3)
    contrasts$F_value[i] <- round(fit[["F"]][1], digits = 3)
    contrasts$df1[i] <- fit$Df[1]
    contrasts$df2[i] <- fit$Df[2]
    contrasts$p_value[i] <- fit$`Pr(>F)`[1]
  }

  ## adjust p-values for multiple comparisons
  contrasts$p_value <- round(p.adjust(contrasts$p_value, method = adj), digits = 3)

  return(list(
    contrasts = contrasts,
    "p-value adjustment" = adj,
    permutations = perm
  ))
}

# This should not highlight differences, and indeed roots we don't have differences.
pairwise_permanova(physeq_ITS_rare_root, "Treatment")

pairwise_adonis_table <-
  rbind(
    pairwise_permanova(physeq_ITS_rare_root, "Treatment")$contrasts %>%
      mutate(Kingdom = "Fungi", Dataset = "Root"),
    pairwise_permanova(physeq_ITS_rare_rhizo, "Treatment")$contrasts %>%
      mutate(Kingdom = "Fungi", Dataset = "Rhizosphere"),
    pairwise_permanova(physeq_ITS_rare_soil, "Treatment")$contrasts %>%
      mutate(Kingdom = "Fungi", Dataset = "Soil"),
    pairwise_permanova(physeq_ITS_rare_control, "Treatment")$contrasts %>%
      mutate(Kingdom = "Fungi", Dataset = "Control"),
    pairwise_permanova(physeq_16S_rare_root, "Treatment")$contrasts %>%
      mutate(Kingdom = "Bacteria", Dataset = "Root"),
    pairwise_permanova(physeq_16S_rare_rhizo, "Treatment")$contrasts %>%
      mutate(Kingdom = "Bacteria", Dataset = "Rhizosphere"),
    pairwise_permanova(physeq_16S_rare_soil, "Treatment")$contrasts %>%
      mutate(Kingdom = "Bacteria", Dataset = "Soil"),
    pairwise_permanova(physeq_16S_rare_control, "Treatment")$contrasts %>%
      mutate(Kingdom = "Bacteria", Dataset = "Control")
  )

pairwise_adonis_table

#  >>>>> TABLE 2 <<<<< ---------------------------------------------------------
write.csv(x = pairwise_adonis_table,
          file = file.path(results_path,"Table2_pairwise_adonis_table.csv"))

# **********************************************************************--------
# *****  PAIRWISE BETADISPER ***** ---------------------------------------------

BetadispExtr(physeq_ITS_rare_root, "bray", "Treatment")$pairwise_df

pairwise_betadisper_table <-
  rbind(
    BetadispExtr(physeq_ITS_rare_root, "bray", "Treatment")$pairwise_df %>%
      mutate(Kingdom = "Fungi", Dataset = "Root"),
    BetadispExtr(physeq_ITS_rare_rhizo, "bray", "Treatment")$pairwise_df %>%
      mutate(Kingdom = "Fungi", Dataset = "Rhizosphere"),
    BetadispExtr(physeq_ITS_rare_soil, "bray", "Treatment")$pairwise_df %>%
      mutate(Kingdom = "Fungi", Dataset = "Soil"),
    BetadispExtr(physeq_ITS_rare_control, "bray", "Treatment")$pairwise_df %>%
      mutate(Kingdom = "Fungi", Dataset = "Control"),
    BetadispExtr(physeq_16S_rare_root, "bray", "Treatment")$pairwise_df %>%
      mutate(Kingdom = "Bacteria", Dataset = "Root"),
    BetadispExtr(physeq_16S_rare_rhizo, "bray", "Treatment")$pairwise_df %>%
      mutate(Kingdom = "Bacteria", Dataset = "Rhizosphere"),
    BetadispExtr(physeq_16S_rare_soil, "bray", "Treatment")$pairwise_df %>%
      mutate(Kingdom = "Bacteria", Dataset = "Soil"),
    BetadispExtr(physeq_16S_rare_control, "bray", "Treatment")$pairwise_df %>%
      mutate(Kingdom = "Bacteria", Dataset = "Control")
  )

pairwise_betadisper_table

write.csv(x = pairwise_betadisper_table, 
          file = file.path(results_path, "Table2_pairwise_betadisper_table.csv"))

# Combine pairwise adonis and betadisper results -------------------------------
# to match groups comparisons

betadisper_ordered <- pairwise_betadisper_table %>%
  mutate(
    group1 = as.character(group1),
    group2 = as.character(group2),
    Group_A = pmin(group1, group2),
    Group_B = pmax(group1, group2)
  ) %>%
  arrange(Group_A, Group_B)


adonis_ordered <- pairwise_adonis_table %>%
  mutate(
    group1 = as.character(group1),
    group2 = as.character(group2),
    Group_A = pmin(group1, group2),
    Group_B = pmax(group1, group2)
  ) %>%
  arrange(Group_A, Group_B)

identical(betadisper_ordered$Group_A, adonis_ordered$Group_A)
identical(betadisper_ordered$Group_B, adonis_ordered$Group_B)

pairwise_combined <- left_join(
  betadisper_ordered,
  adonis_ordered,
  by = c("Group_A", "Group_B", "Kingdom", "Dataset"),
  suffix = c("_betadisper", "_adonis")
)

pairwise_combined

write.csv(x = pairwise_combined,
          file = file.path(results_path,"Table2_pairwise_combined_table.csv"))

# **********************************************************************--------
# ***** ALPHA DIVERSITY ***** --------------------------------------------------

AlphaMetrics <- function(physeq) {
  require(vegan)
  require(tidyverse)

  sample_data(physeq)$ReadNo <- sample_sums(physeq)
  sample_data(physeq)$hill_0 <- as.data.frame(as.matrix(t(physeq@otu_table))) %>% renyi(scale = c(0), hill = TRUE)
  sample_data(physeq)$hill_1 <- as.data.frame(as.matrix(t(physeq@otu_table))) %>% renyi(scale = c(1), hill = TRUE)
  sample_data(physeq)$hill_2 <- as.data.frame(as.matrix(t(physeq@otu_table))) %>% renyi(scale = c(2), hill = TRUE)
  sample_data(physeq)$Richness <- specnumber(as.data.frame(otu_table(physeq)), MARGIN = 2)
  sample_data(physeq)$invSimpson <- diversity(as.data.frame(otu_table(physeq)), index = "inv", MARGIN = 2)
  sample_data(physeq)$Shannon <- diversity(as.data.frame(otu_table(physeq)), index = "shannon", MARGIN = 2)
  sample_data(physeq)$EH <- 1 - sample_data(physeq)$Shannon / log(sample_data(physeq)$Richness)

  return(physeq)
}

# Running this will the metrics to the @sam_data table
AlphaMetrics(physeq_ITS_rare_root) %>% sample_data()
AlphaMetrics(physeq_ITS_rare_rhizo) %>% sample_data()
AlphaMetrics(physeq_ITS_rare_soil) %>% sample_data()
AlphaMetrics(physeq_ITS_rare_control) %>% sample_data()

AlphaMetrics(physeq_16S_rare_root) %>% sample_data()
AlphaMetrics(physeq_16S_rare_rhizo) %>% sample_data()
AlphaMetrics(physeq_16S_rare_soil) %>% sample_data()
AlphaMetrics(physeq_16S_rare_control) %>% sample_data()


# Testing for significant differences ------------------------------------------

# Testing for significant differences ------------------------------------------
# Using wilcox.test
CompSampl <- function(df, formula, comparisons) {
  require(multcompView)
  require(lazyeval)

  test_CC <-
    compare_means(formula, data = df, method = "wilcox.test", p.adjust.method = "none")

  test_CC$adj.pval <-
    p.adjust(test_CC$p, method = "BH", n = comparisons)

  # print(test_CC)

  test_CC <-
    as.data.frame(test_CC)[, c(2, 3, 5)] # to change form p to p.adj do 4 to 5
  test_CC2 <-
    data.frame(test_CC[, 2], test_CC[, 1], test_CC[, 3])
  colnames(test_CC2) <-
    c("group1", "group2", "p.adj") # change p to p.adj
  rbind(test_CC, test_CC2) -> test_all

  print(test_all)

  dist_CC <- as.dist(xtabs(test_all[, 3] ~ (test_all[, 2] + test_all[, 1])), diag = TRUE)
  # print(dist_CC)
  res_CC <-
    data.frame(multcompLetters(dist_CC)["Letters"])
  return(res_CC)
}


# Calculate number of multiple pairwise comparisons based on the number of levels
PairCalc <- function(n) {
  n_comparison <- n * (n - 1) / 2
  return(n_comparison)
}

PairCalc(6)

# Extracting multiple comparisons and p.adj
physeq_ITS_rare_root@sam_data %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(hill_0 = as.numeric(hill_0)) %>%
  as_tibble() %>%
  CompSampl(., formula(hill_0 ~ Treatment), comparisons = PairCalc(4))

# plotting Richness function ---------------------------------------------------
PlotRich <- function(physeq, X_var, Y_var, my_labels, labels_y) {
  require(phyloseq)
  require(tidyverse)

  # extract dataframe
  dataframe <-
    physeq@sam_data %>%
    as.matrix() %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(!!Y_var := as.numeric(!!sym(Y_var)))

  # print factor order

  # Calculate labels_y based on the maximum value of Y_var
  labels_y <- max(dataframe[, Y_var]) + 0.1 * max(dataframe[, Y_var])

  # plot
  rich_plot <-
    ggplot(dataframe, aes(x = get(X_var), y = !!sym(Y_var))) +
    geom_jitter(
      position = position_jitter(0.4), size = 1, shape = 16,
      aes(color = get(X_var))
    ) +
    # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
    #             geom="pointrange", color="black", shape=18, size=0.8) +
    stat_summary(
      fun = median, fun.min = median, fun.max = median,
      geom = "crossbar", width = 0.8, linewidth = 0.3, color = "black"
    ) +
    # stat_summary(
    #  geom="pointrange",
    #  fun.min = function(z) { quantile(z, 0.25) },
    #  fun.max = function(z) { quantile(z, 0.75) },
    #  fun = median, color="black", shape=18, size=1,
    #  show.legend = FALSE) +
    stat_summary(
      geom = "text", angle = 0, label = my_labels,
      fun = max, aes(y = labels_y), size = 3, color = "black"
    ) +
    # scale_color_manual(values=c("grey", "blue")) +
    expand_limits(y = 0) +
    theme_bw() +
    theme(
      plot.title = element_markdown(size = 8, face = "bold", vjust = 0.5, hjust = 0.5),
      axis.title = element_markdown(size = 7),
      axis.text.x = element_markdown(size = 7, colour = "black", angle = 33, hjust = 1, vjust = 1),
      axis.text.y = element_markdown(size = 7, angle = 0, hjust = 0.5),
      legend.position = "none"
    )

  return(rich_plot)
}

# Test the plotting function
PlotRich(
  physeq_ITS_rare_root, "Treatment", "hill_0",
  physeq_ITS_rare_root@sam_data %>%
    as.matrix() %>%
    as.data.frame() %>%
    mutate(hill_0 = as.numeric(hill_0)) %>%
    CompSampl(., formula(hill_0 ~ Treatment), comparisons = PairCalc(4)) %>%
    pull(Letters),
  200
) +
  scale_color_manual(values = c("#825121", "#CC2D35", "#FFB400", "#00A6ED"))

# Multiplot Fungi --------------------------------------------------------------

# Hill_0 -----------------------------------------------------------------------
plot_hill_0 <-
  ggarrange(
    ggarrange(
      grid.arrange(
        ggarrange(
          PlotRich(
            physeq_ITS_rare_root, "Treatment", "hill_0",
            physeq_ITS_rare_root@sam_data %>%
              as.matrix() %>%
              as.data.frame() %>%
              mutate(hill_0 = as.numeric(hill_0)) %>%
              CompSampl(., formula(hill_0 ~ Treatment), comparisons = PairCalc(4)) %>%
              pull(Letters),
            200
          ) +
            scale_color_manual(values = c("#825121", "#CC2D35", "#FFB400", "#00A6ED")) +
            labs(title = "Root", x = NULL),
          PlotRich(
            physeq_ITS_rare_rhizo, "Treatment", "hill_0",
            physeq_ITS_rare_rhizo@sam_data %>%
              as.matrix() %>%
              as.data.frame() %>%
              mutate(hill_0 = as.numeric(hill_0)) %>%
              CompSampl(., formula(hill_0 ~ Treatment), comparisons = PairCalc(4)) %>%
              pull(Letters),
            600
          ) +
            scale_color_manual(values = c("#825121", "#CC2D35", "#FFB400", "#00A6ED")) +
            labs(title = "Rhizosphere", x = NULL),
          PlotRich(
            physeq_ITS_rare_soil, "Treatment", "hill_0",
            physeq_ITS_rare_soil@sam_data %>%
              as.matrix() %>%
              as.data.frame() %>%
              mutate(hill_0 = as.numeric(hill_0)) %>%
              CompSampl(., formula(hill_0 ~ Treatment), comparisons = PairCalc(3)) %>%
              pull(Letters),
            600
          ) +
            scale_color_manual(values = c("#825121", "#CC79A7", "#FFB400")) +
            labs(title = "Soil", x = NULL),
          PlotRich(
            physeq_ITS_rare_control, "Treatment", "hill_0",
            physeq_ITS_rare_control@sam_data %>%
              as.matrix() %>%
              as.data.frame() %>%
              mutate(hill_0 = as.numeric(hill_0)) %>%
              CompSampl(., formula(hill_0 ~ Treatment), comparisons = PairCalc(2)) %>%
              pull(Letters),
            600
          ) +
            scale_color_manual(values = c("#7FB800", "#2D3142")) +
            labs(title = "Control", x = NULL),
          nrow = 4,
          ncol = 1,
          # align = "hv",
          labels = c("A", "B", "C", "D")
        ),
        top = title1
      ),
      grid.arrange(
        ggarrange(
          PlotRich(
            physeq_16S_rare_root, "Treatment", "hill_0",
            physeq_16S_rare_root@sam_data %>%
              as.matrix() %>%
              as.data.frame() %>%
              mutate(hill_0 = as.numeric(hill_0)) %>%
              CompSampl(., formula(hill_0 ~ Treatment), comparisons = PairCalc(4)) %>%
              pull(Letters),
            500
          ) +
            scale_color_manual(values = c("#825121", "#CC2D35", "#FFB400", "#00A6ED")) +
            labs(title = "Root", x = NULL),
          PlotRich(
            physeq_16S_rare_rhizo, "Treatment", "hill_0",
            physeq_16S_rare_rhizo@sam_data %>%
              as.matrix() %>%
              as.data.frame() %>%
              mutate(hill_0 = as.numeric(hill_0)) %>%
              CompSampl(., formula(hill_0 ~ Treatment), comparisons = PairCalc(4)) %>%
              pull(Letters),
            5000
          ) +
            scale_color_manual(values = c("#825121", "#CC2D35", "#FFB400", "#00A6ED")) +
            labs(title = "Rhizosphere", x = NULL),
          PlotRich(
            physeq_16S_rare_soil, "Treatment", "hill_0",
            physeq_16S_rare_soil@sam_data %>%
              as.matrix() %>%
              as.data.frame() %>%
              mutate(hill_0 = as.numeric(hill_0)) %>%
              CompSampl(., formula(hill_0 ~ Treatment), comparisons = PairCalc(3)) %>%
              pull(Letters),
            5000
          ) +
            scale_color_manual(values = c("#825121", "#CC79A7", "#FFB400")) +
            labs(title = "Soil", x = NULL),
          PlotRich(
            physeq_16S_rare_control, "Treatment", "hill_0",
            physeq_16S_rare_control@sam_data %>%
              as.matrix() %>%
              as.data.frame() %>%
              mutate(hill_0 = as.numeric(hill_0)) %>%
              CompSampl(., formula(hill_0 ~ Treatment), comparisons = PairCalc(2)) %>%
              pull(Letters),
            5000
          ) +
            scale_color_manual(values = c("#7FB800", "#2D3142")) +
            labs(title = "Control", x = NULL),
          nrow = 4,
          ncol = 1,
          # align = "hv",
          labels = c("E", "F", "G", "H")
        ),
        top = title2
      ),
      nrow = 1,
      ncol = 2
    ),
    # MakeLegend() %>% as_ggplot(),
    # nrow = 2,
    # ncol = 1,
    heights = c(.95, 0.07)
  )


plot_hill_0

# Hill 1 -----------------------------------------------------------------------
plot_hill_1 <-
  ggarrange(
    ggarrange(
      grid.arrange(
        ggarrange(
          PlotRich(
            physeq_ITS_rare_root, "Treatment", "hill_1",
            physeq_ITS_rare_root@sam_data %>%
              as.matrix() %>%
              as.data.frame() %>%
              mutate(hill_1 = as.numeric(hill_1)) %>%
              CompSampl(., formula(hill_1 ~ Treatment), comparisons = PairCalc(4)) %>%
              pull(Letters),
            200
          ) +
            scale_color_manual(values = c("#825121", "#CC2D35", "#FFB400", "#00A6ED")) +
            labs(title = "Root", x = NULL),
          PlotRich(
            physeq_ITS_rare_rhizo, "Treatment", "hill_1",
            physeq_ITS_rare_rhizo@sam_data %>%
              as.matrix() %>%
              as.data.frame() %>%
              mutate(hill_1 = as.numeric(hill_1)) %>%
              CompSampl(., formula(hill_1 ~ Treatment), comparisons = PairCalc(4)) %>%
              pull(Letters),
            600
          ) +
            scale_color_manual(values = c("#825121", "#CC2D35", "#FFB400", "#00A6ED")) +
            labs(title = "Rhizosphere", x = NULL),
          PlotRich(
            physeq_ITS_rare_soil, "Treatment", "hill_1",
            physeq_ITS_rare_soil@sam_data %>%
              as.matrix() %>%
              as.data.frame() %>%
              mutate(hill_1 = as.numeric(hill_1)) %>%
              CompSampl(., formula(hill_1 ~ Treatment), comparisons = PairCalc(3)) %>%
              pull(Letters),
            600
          ) +
            scale_color_manual(values = c("#825121", "#CC79A7", "#FFB400")) +
            labs(title = "Soil", x = NULL),
          PlotRich(
            physeq_ITS_rare_control, "Treatment", "hill_1",
            physeq_ITS_rare_control@sam_data %>%
              as.matrix() %>%
              as.data.frame() %>%
              mutate(hill_1 = as.numeric(hill_1)) %>%
              CompSampl(., formula(hill_1 ~ Treatment), comparisons = PairCalc(2)) %>%
              pull(Letters),
            600
          ) +
            scale_color_manual(values = c("#7FB800", "#2D3142")) +
            labs(title = "Control", x = NULL),
          nrow = 4,
          ncol = 1,
          # align = "hv",
          labels = c("A", "B", "C", "D")
        ),
        top = title1
      ),
      grid.arrange(
        ggarrange(
          PlotRich(
            physeq_16S_rare_root, "Treatment", "hill_1",
            physeq_16S_rare_root@sam_data %>%
              as.matrix() %>%
              as.data.frame() %>%
              mutate(hill_1 = as.numeric(hill_1)) %>%
              CompSampl(., formula(hill_1 ~ Treatment), comparisons = PairCalc(6)) %>%
              pull(Letters),
            500
          ) +
            scale_color_manual(values = c("#825121", "#CC2D35", "#FFB400", "#00A6ED")) +
            labs(title = "Root", x = NULL),
          PlotRich(
            physeq_16S_rare_rhizo, "Treatment", "hill_1",
            physeq_16S_rare_rhizo@sam_data %>%
              as.matrix() %>%
              as.data.frame() %>%
              mutate(hill_1 = as.numeric(hill_1)) %>%
              CompSampl(., formula(hill_1 ~ Treatment), comparisons = PairCalc(6)) %>%
              pull(Letters),
            5000
          ) +
            scale_color_manual(values = c("#825121", "#CC2D35", "#FFB400", "#00A6ED")) +
            labs(title = "Rhizosphere", x = NULL),
          PlotRich(
            physeq_16S_rare_soil, "Treatment", "hill_1",
            physeq_16S_rare_soil@sam_data %>%
              as.matrix() %>%
              as.data.frame() %>%
              mutate(hill_1 = as.numeric(hill_1)) %>%
              CompSampl(., formula(hill_1 ~ Treatment), comparisons = PairCalc(3)) %>%
              pull(Letters),
            5000
          ) +
            scale_color_manual(values = c("#825121", "#CC79A7", "#FFB400")) +
            labs(title = "Soil", x = NULL),
          PlotRich(
            physeq_16S_rare_control, "Treatment", "hill_1",
            physeq_16S_rare_control@sam_data %>%
              as.matrix() %>%
              as.data.frame() %>%
              mutate(hill_1 = as.numeric(hill_1)) %>%
              CompSampl(., formula(hill_1 ~ Treatment), comparisons = PairCalc(2)) %>%
              pull(Letters),
            5000
          ) +
            scale_color_manual(values = c("#7FB800", "#2D3142")) +
            labs(title = "Control", x = NULL),
          nrow = 4,
          ncol = 1,
          # align = "hv",
          labels = c("E", "F", "G", "H")
        ),
        top = title2
      ),
      nrow = 1,
      ncol = 2
    ),
    # MakeLegend() %>% as_ggplot(),
    # nrow = 2,
    # ncol = 1,
    heights = c(.95, 0.07)
  )


plot_hill_1

#  >>>>> FIGURE 4 <<<<< --------------------------------------------------------

# Multiplot hill_0 and hill_1

Fig_4_alphadiv <- ggarrange(
  ggarrange(
    grid.arrange(plot_hill_0, top = text_grob("Hill 0", size = 14, face = 2)),
    grid.arrange(plot_hill_1, top = text_grob("Hill 1", size = 14, face = 2)),
    ncol = 2,
    nrow = 1
  ),
  MakeLegend() %>% as_ggplot(),
  nrow = 2,
  heights = c(.95, 0.07)
)

Fig_4_alphadiv

ggsave(plot = Fig_4_alphadiv, 
       path = results_path,
       filename = "Fig_4_alpha_diversity.pdf")

# **********************************************************************--------
# DIFFERENTIAL ABUNDANCE -------------------------------------------------------

# NOTE Due to the difficulties in running this analysis at OTU level or higher, I
# used the BestMatch column in the data set as my rank for this comparative analysis.
# Please see the manuscript M&M for details and rationale.
# Additionally, this analysis was run only on the pairwise comparisons that were
# found to be significant in the pairwise PERMANOVA for consistency and logic.

# Splitting datasets to compared groups.

# fungi
physeq_ITS_rhizo_IE_NE <-
  subset_samples(
    physeq_ITS_rare,
    Compartment %in% c("Rhizosphere") & Treatment %in% c("IE", "NE")
  ) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_ITS_rhizo_IE_NE@sam_data
head(physeq_ITS_rhizo_IE_NE@tax_table)

physeq_ITS_rhizo_NE_NE_MIX <-
  subset_samples(
    physeq_ITS_rare,
    Compartment %in% c("Rhizosphere") & Treatment %in% c("NE", "NE_MIX")
  ) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)


physeq_ITS_rhizo_NE_NE_MIX
head(physeq_ITS_rhizo_NE_NE_MIX@tax_table)

# bacteria
physeq_16S_root_IE_IE_MIX <-
  subset_samples(
    physeq_16S_rare_root,
    Compartment %in% c("Root") & Treatment %in% c("IE", "IE_MIX")
  ) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_16S_root_IE_NE_MIX <-
  subset_samples(
    physeq_16S_rare_root,
    Compartment %in% c("Root") & Treatment %in% c("IE", "NE_MIX")
  ) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_16S_root_NE_IE_MIX <-
  subset_samples(
    physeq_16S_rare_root,
    Compartment %in% c("Root") & Treatment %in% c("NE", "IE_MIX")
  ) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_16S_root_NE_NE_MIX <-
  subset_samples(
    physeq_16S_rare_root,
    Compartment %in% c("Root") & Treatment %in% c("NE", "NE_MIX")
  ) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

# rhizo and control
physeq_16S_rhizo_IE_NE_MIX <-
  subset_samples(
    physeq_16S_rare_rhizo,
    Compartment %in% c("Rhizosphere") & Treatment %in% c("IE", "NE_MIX")
  ) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_16S_rhizo_IE_MIX__NE_MIX <-
  subset_samples(
    physeq_16S_rare_rhizo,
    Compartment %in% c("Rhizosphere") & Treatment %in% c("IE_MIX", "NE_MIX")
  ) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

# control
physeq_16S_soil_fumigated_control <-
  subset_samples(
    physeq_16S_rare_control,
    Compartment %in% c("soil") & Treatment %in% c("fumigated", "control")
  ) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)


# Extract a long dataframe physeq object ---------------------------------------
ExtractLongDF <- function(physeq) {
  require(phyloseq)
  require(tidyverse)

  otu <- physeq@otu_table %>%
    as.matrix() %>%
    as.data.frame()
  taxonomy <- physeq@tax_table %>%
    as.matrix() %>%
    as.data.frame()
  metadata <- physeq@sam_data %>%
    as.matrix() %>%
    as.data.frame()

  long_df <-
    otu %>%
    rownames_to_column("OTU_ID") %>%
    left_join(., taxonomy %>%
      tidyr::separate(Taxonomy,
        into = c("OTU_ID", "BestMatch"),
        sep = " ", extra = "merge", fill = "right", remove = FALSE
      ) %>%
      dplyr::select(OTU_ID, Taxonomy, BestMatch), by = "OTU_ID") %>%
    dplyr::select(.data = ., starts_with("sam"), OTU_ID, Taxonomy, BestMatch) %>%
    group_by(BestMatch, Taxonomy, OTU_ID) %>%
    pivot_longer(cols = c(-OTU_ID, -Taxonomy, -BestMatch), names_to = "Sample_ID", values_to = "Count") %>%
    mutate(Abund = Count / sum(Count)) %>%
    inner_join(
      x = .,
      y = metadata %>% rownames_to_column("Sample_ID"),
      by = "Sample_ID"
    ) %>%
    ungroup()

  return(long_df)
}


# test the function
df_long_ITS_NE_NE_MIX <- ExtractLongDF(physeq_ITS_rhizo_NE_NE_MIX)
df_long_ITS_NE_NE_MIX

# Test using wilcox.test -------------------------------------------------------
GetWilcoxPval <- function(df) {
  require(tidyverse)
  require(broom)

  summary_df <- df %>%
    group_by(BestMatch, Treatment) %>%
    summarise(
      Sum = sum(Count, na.rm = TRUE),
      median_abund = median(Count, na.rm = TRUE),
      mean_abund = mean(Count, na.rm = TRUE),
      sd_abund = sd(Count, na.rm = TRUE),
      group_size = n(),
      present = sum(Count > 0),
      # total_present = sum(present),
      ratio = present / group_size,
      # prevalence = sum(Count > 0, na.rm = TRUE) / group_size,
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Treatment,
      values_from = c(
        Sum, median_abund, mean_abund, sd_abund,
        group_size, present, ratio
      ),
      names_glue = "{.value}_{Treatment}",
      values_fn = mean # or first, median, etc.
    )

  res_pval <- df %>%
    group_by(BestMatch) %>%
    nest(data = -BestMatch) %>%
    mutate(test = map(data, ~ wilcox.test(Count ~ Treatment, data = .x, exact = FALSE) %>%
      tidy())) %>%
    unnest(test) %>%
    mutate(p.adjust = round(p.adjust(p.value, method = "BH"), 4)) %>%
    filter(p.adjust <= 0.05)

  # Step 3: Join with summary stats
  final_result <- res_pval %>%
    left_join(summary_df, by = "BestMatch") %>%
    as.data.frame()

  return(final_result)
}


# fungi
diff_ITS_NE_NE_MIX <-
  ExtractLongDF(physeq_ITS_rhizo_NE_NE_MIX) %>%
  GetWilcoxPval(.) %>%
  filter(if_all(starts_with("Sum_"), ~ . > 0)) %>% # removing groups with 0 in one treatment
  mutate(Group = "NE-NE_MIX")

diff_ITS_NE_NE_MIX

# NOTE some taxa were excluded fo simplicity and because the taxonomic classification
# was poor. With no good/reliable taxonomic classification any biological conclusion \
# can be misinterpreted.

diff_ITS_NE_IE <-
  ExtractLongDF(physeq_ITS_rhizo_IE_NE) %>%
  GetWilcoxPval(.) %>%
  filter(if_all(starts_with("Sum_"), ~ . > 0)) %>%
  mutate(Group = "NE-IE") %>%
  filter(!BestMatch %in% c("Fungi", "Rozellomycota", "Endogonomycetes", "Agaricomycetes"))

diff_ITS_NE_IE

# bacteria
diff_16S_root_IE_IE_MIX <-
  ExtractLongDF(physeq_16S_root_IE_IE_MIX) %>%
  GetWilcoxPval(.) %>%
  filter(if_all(starts_with("Sum_"), ~ . > 0)) %>%
  mutate(Group = "IE-IE_MIX")

diff_16S_root_IE_IE_MIX

diff_16S_root_IE_IE_MIX$data[[4]]$Taxonomy
diff_16S_root_IE_IE_MIX$data[[6]]$Taxonomy
diff_16S_root_IE_IE_MIX$data[[7]]$Taxonomy

diff_16S_root_IE_IE_MIX <-
  diff_16S_root_IE_IE_MIX %>%
  filter(!BestMatch %in% c(
    "Bacteria", "uncultured bacterium 946",
    "metagenome69", "Uncultured bacterium 500"
  ))

diff_16S_root_IE_NE_MIX <-
  ExtractLongDF(physeq_16S_root_IE_NE_MIX) %>%
  GetWilcoxPval(.) %>%
  filter(if_all(starts_with("Sum_"), ~ . > 0)) %>%
  mutate(Group = "IE-NE_MIX")

diff_16S_root_IE_NE_MIX

diff_16S_root_IE_NE_MIX$data[[2]]$Taxonomy
diff_16S_root_IE_NE_MIX$data[[6]]$Taxonomy

diff_16S_root_IE_NE_MIX <-
  diff_16S_root_IE_NE_MIX %>%
  filter(!BestMatch %in% c("uncultured bacterium121", "Uncultured bacterium166"))


diff_16S_root_NE_IE_MIX <-
  ExtractLongDF(physeq_16S_root_NE_IE_MIX) %>%
  GetWilcoxPval(.) %>%
  filter(if_all(starts_with("Sum_"), ~ . > 0)) %>%
  mutate(Group = "NE-IE_MIX")

diff_16S_root_NE_IE_MIX

diff_16S_root_NE_IE_MIX <-
  diff_16S_root_NE_IE_MIX %>%
  filter(!BestMatch %in% c("uncultured bacterium 2722", "D05-2"))


diff_16S_root_NE_NE_MIX <-
  ExtractLongDF(physeq_16S_root_NE_NE_MIX) %>%
  GetWilcoxPval(.) %>%
  filter(if_all(starts_with("Sum_"), ~ . > 0)) %>%
  mutate(Group = "NE-NE_MIX")

diff_16S_root_NE_NE_MIX

diff_16S_root_NE_NE_MIX <-
  diff_16S_root_NE_NE_MIX %>%
  filter(!BestMatch %in% c("uncultured bacterium 689", "Uncultured 43", "D05-2"))

# rhizo and soil controls
diff_16S_rhizo_IE_NE_MIX <-
  ExtractLongDF(physeq_16S_rhizo_IE_NE_MIX) %>%
  GetWilcoxPval(.) %>%
  filter(if_all(starts_with("Sum_"), ~ . > 0)) %>%
  mutate(Group = "IE-NE_MIX")

diff_16S_rhizo_IE_NE_MIX

diff_16S_rhizo_IE_NE_MIX <-
  diff_16S_rhizo_IE_NE_MIX %>%
  filter(!str_detect(BestMatch, regex("^uncultured", ignore_case = TRUE))) %>%
  filter(!str_detect(BestMatch, regex("^metagenome", ignore_case = TRUE))) %>%
  filter(!BestMatch %in% c(
    "R7c24", "Sm2d12", "S0134 terrestrial group",
    "Bacterium enrichment culture clone auto9 4w",
    "Jtb23", "Ns11-12 marine group"
  ))


diff_16S_rhizo_IE_MIX__NE_MIX <-
  ExtractLongDF(physeq_16S_rhizo_IE_MIX__NE_MIX) %>%
  GetWilcoxPval(.) %>%
  filter(if_all(starts_with("Sum_"), ~ . > 0)) %>%
  mutate(Group = "IE_MIX-NE_MIX")

diff_16S_rhizo_IE_MIX__NE_MIX

diff_16S_rhizo_IE_MIX__NE_MIX <-
  diff_16S_rhizo_IE_MIX__NE_MIX %>%
  filter(!str_detect(BestMatch, regex("^uncultured", ignore_case = TRUE))) %>%
  filter(!str_detect(BestMatch, regex("^metagenome", ignore_case = TRUE))) %>%
  filter(!str_detect(BestMatch, regex("-", ignore_case = TRUE))) %>%
  filter(!BestMatch %in% c(
    "Olb13", "cidobacteria bacterium cb 286306",
    "Subgroup7", "Wd2101 soil group",
    "Acidobacteria bacterium cb 286306",
    "Uncultivated soil bacterium clone c028",
    "Mbnt15", "Dev008", "Sbr1031"
  )) %>%
  filter(!BestMatch %in% c(
    "Candidatus Protochlamydia", "Proteobacteria",
    "Bdellovibrio"
  ))

# soil
diff_16S_soil_fumigated_control <-
  ExtractLongDF(physeq_16S_soil_fumigated_control) %>%
  GetWilcoxPval(.) %>%
  filter(if_all(starts_with("Sum_"), ~ . > 0)) %>%
  mutate(Group = "fumigated-control")

diff_16S_soil_fumigated_control

diff_16S_soil_fumigated_control <-
  diff_16S_soil_fumigated_control %>%
  filter(!str_detect(BestMatch, regex("^uncultured", ignore_case = TRUE))) %>%
  filter(!str_detect(BestMatch, regex("^metagenome", ignore_case = TRUE))) %>%
  filter(!str_detect(BestMatch, regex("-", ignore_case = TRUE))) %>%
  filter(!BestMatch %in% c(
    "Birii41", "Mnd1", "Ellin516", "Clostridium sensu stricto 6",
    "Swb02", "A4b", "Lineage iv", "Pla4 lineage",
    "Akyh767", "Wastewater metagenome 21", "Ellin516",
    "Bacteroidetes bacterium 20", "Env.ops7",
    "Akyh767", "Possible genus 04", "Sm2d12"
  )) %>%
  filter(!BestMatch %in% c(
    "Woesearchaeales", "Parachlamydiaceae", "Coxiella",
    "Bdellovibrio"
  ))

# plotting BestMatch -----------------------------------------------------------

PlotBestMatch <- function(physeq, diff_wilcox) {
  require(phyloseq)
  require(tidyverse)

  plot_diff <-
    psmelt(physeq) %>%
    filter(BestMatch %in% diff_wilcox$BestMatch) %>%
    mutate(BestMatch = str_trim(BestMatch)) %>%
    mutate(BestMatch = if_else(
      str_ends(BestMatch, "ceae") |
        str_ends(BestMatch, "ales") |
        str_ends(BestMatch, "etes") |
        str_ends(BestMatch, "ota"),
      BestMatch,
      str_replace(BestMatch, "(.*)", "*\\1*")
    )) %>%
    group_by(BestMatch, Sample, Treatment) %>%
    ggplot(aes(x = Abundance + 1e-4, y = BestMatch, colour = Treatment, group = interaction(BestMatch, Treatment))) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), shape = 16, size = 1.5) +
    geom_boxplot(outlier.shape = NA, fatten = 3, alpha = 0.4, color = "black") +
    # stat_summary(fun.data = median_hilow, fun.args = list(conf.int=0.5),
    #             geom="pointrange", shape = 5, size = 0.5, linewidth = 0.4,stroke = 0.4,
    #             position = position_dodge(width=0.8),
    #             color="black", show.legend = FALSE) +
    # stat_summary(fun.data = mean, geom="point", shape = 1, color = "blue", show.legend = FALSE) +
    scale_x_log10() +
    theme_bw() +
    theme(
      axis.text.y = element_markdown(size = 8, angle = 0, hjust = 1, vjust = 0.5),
      axis.text.x = element_markdown(size = 8, hjust = 1, vjust = 0),
      plot.title = element_markdown(hjust = 1),
      legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
      legend.title = element_blank(), legend.background = element_blank(),
      legend.text = element_markdown(size = 8),
      legend.position = "right"
    ) +
    guides(colour = guide_legend(ncol = 1, override.aes = list(shape = 15, size = 3))) +
    labs(y = NULL, x = "log10(Abundance)")

  return(plot_diff)
}


# fungi
plot_ITS_NE_NE_MIX_no_zero <- PlotBestMatch(physeq_ITS_rhizo_NE_NE_MIX, diff_ITS_NE_NE_MIX)
plot_ITS_NE_NE_MIX_no_zero

plot_ITS_NE_IE_no_zero <- PlotBestMatch(physeq_ITS_rhizo_IE_NE, diff_ITS_NE_IE)
plot_ITS_NE_IE_no_zero

# bacteria root
plot_16S_IE_IE_MIX_no_zero <- PlotBestMatch(physeq_16S_root_IE_IE_MIX, diff_16S_root_IE_IE_MIX)
plot_16S_IE_IE_MIX_no_zero

plot_16S_IE_NE_MIX_no_zero <- PlotBestMatch(physeq_16S_root_IE_NE_MIX, diff_16S_root_IE_NE_MIX)
plot_16S_IE_NE_MIX_no_zero

plot_16S_NE_IE_MIX_no_zero <- PlotBestMatch(physeq_16S_root_NE_IE_MIX, diff_16S_root_NE_IE_MIX)
plot_16S_NE_IE_MIX_no_zero

plot_16S_NE_NE_MIX_no_zero <- PlotBestMatch(physeq_16S_root_NE_NE_MIX, diff_16S_root_NE_NE_MIX)
plot_16S_NE_NE_MIX_no_zero

# Rhizosphere and soil
plot_16S_rhizo_IE_NE_MIX <- PlotBestMatch(physeq_16S_rhizo_IE_NE_MIX, diff_16S_rhizo_IE_NE_MIX)
plot_16S_rhizo_IE_NE_MIX

plot_16S_rhizo_IE_MIX__NE_MIX <- PlotBestMatch(physeq_16S_rhizo_IE_MIX__NE_MIX, diff_16S_rhizo_IE_MIX__NE_MIX)
plot_16S_rhizo_IE_MIX__NE_MIX

plot_16S_soil_fumigated_control <- PlotBestMatch(physeq_16S_soil_fumigated_control, diff_16S_soil_fumigated_control)
plot_16S_soil_fumigated_control

# checks
dplyr::intersect(diff_16S_soil_fumigated_control$BestMatch, diff_16S_rhizo_IE_NE_MIX$BestMatch)
dplyr::intersect(diff_16S_soil_fumigated_control$BestMatch, diff_16S_rhizo_IE_NE_MIX %>%
  separate(BestMatch, c("Genus", "Species"), sep = " ", remove = FALSE) %>%
  pull(Genus))

dplyr::intersect(diff_16S_soil_fumigated_control$BestMatch, diff_16S_rhizo_IE_MIX__NE_MIX$BestMatch)
dplyr::intersect(diff_16S_soil_fumigated_control$BestMatch, diff_16S_rhizo_IE_MIX__NE_MIX %>%
  separate(BestMatch, c("Genus", "Species"), sep = " ", remove = FALSE) %>%
  pull(Genus))

# Composite plot --------------------
data.frame(
  label = c("IE", "NE", "IE_MIX", "NE_MIX", "IE_NE_MIX", "Fumigated", "Control"),
  color = c("#825121", "#FFB400", "#CC2D35", "#00A6ED", "#CC79A7", "#2D3142", "#7FB800")
) %>%
  ggplot(aes(x = label, y = 1, fill = color)) +
  geom_col() +
  scale_fill_identity() +
  theme(
    axis.text.x = element_text(angle = 90, size = 20, hjust = 1, vjust = 0.5),
    legend.position = "none"
  )

# fungi

fungi_wilcox <-
  ggarrange(
    plot_ITS_NE_NE_MIX_no_zero +
      scale_color_manual(values = c("#FFB400", "#00A6ED")),
    plot_ITS_NE_IE_no_zero +
      scale_color_manual(values = c("#825121", "#FFB400")),
    ncol = 1,
    nrow = 2, labels = c("A", "B"),
    align = "v",
    heights = c(12, 10)
  )

fungi_wilcox

#  >>>>> FIGURE 5 <<<<< --------------------------------------------------------
ggsave(
  path = results_path,
  filename = "Fig_5_diff_abund_fungi_rhizosphere.pdf",
  plot =
    grid.arrange(fungi_wilcox,
      top = text_grob("Rhizosphere differentially abundant fungi",
        size = 12, face = 2
      )
    ),
  device = "pdf"
)

# Saving 5.02 x 5.55 in image

bacteria_root_wilcox <-
  ggarrange(
    plot_16S_NE_NE_MIX_no_zero +
      scale_color_manual(values = c("#FFB400", "#00A6ED")),
    plot_16S_IE_NE_MIX_no_zero +
      scale_color_manual(values = c("#825121", "#00A6ED")),
    plot_16S_NE_IE_MIX_no_zero +
      scale_color_manual(values = c("#CC2D35", "#FFB400")),
    plot_16S_IE_IE_MIX_no_zero +
      scale_color_manual(values = c("#825121", "#CC2D35")),
    ncol = 1,
    nrow = 4,
    align = "v",
    labels = c("A", "B", "C", "D"),
    heights = c(7, 6, 7, 7)
  )

bacteria_root_wilcox

#  >>>>> FIGURE 6 <<<<< --------------------------------------------------------
ggsave(
  path = results_path,
  filename = "Fig_6_diff_abund_bacteria_root_wilcox.pdf",
  plot =
    grid.arrange(bacteria_root_wilcox,
      top = text_grob("Root differentially abundant bacteria",
        size = 12, face = 2
      )
    ),
  device = "pdf"
)

# Saving 5.02 x 8.76 in image

bacteria_rhizo_soil_wilcox <-
  ggarrange(
    plot_16S_rhizo_IE_NE_MIX +
      scale_color_manual(values = c("#825121", "#00A6ED")) +
      theme(
        legend.position = "bottom",
        plot.subtitle = element_text(vjust = 0.5, hjust = 0.5)
      ) +
      guides(colour = guide_legend(
        ncol = 2, override.aes = list(shape = 15, size = 3)
      )) +
      labs(subtitle = "Rhizosphere"),
    plot_16S_rhizo_IE_MIX__NE_MIX +
      scale_color_manual(values = c("#CC2D35", "#00A6ED")) +
      theme(
        legend.position = "bottom",
        plot.subtitle = element_text(vjust = 0.5, hjust = 0.5)
      ) +
      guides(colour = guide_legend(
        ncol = 2, override.aes = list(shape = 15, size = 3)
      )) +
      labs(subtitle = "Rhizosphere"),
    plot_16S_soil_fumigated_control +
      scale_color_manual(values = c("#2D3142", "#7FB800")) +
      theme(
        legend.position = "bottom",
        plot.subtitle = element_text(vjust = 0.5, hjust = 0.5)
      ) +
      guides(colour = guide_legend(
        ncol = 2, override.aes = list(shape = 15, size = 3)
      )) +
      labs(subtitle = "Soil"),
    ncol = 3,
    nrow = 1,
    align = "v",
    labels = c("A", "B", "C")
  )

bacteria_rhizo_soil_wilcox

#  >>>>> FIGURE 7 <<<<< --------------------------------------------------------
ggsave(
  path = results_path,
  filename = "Fig_7_diff_abund_bacteria_rhizo_soil.pdf",
  plot = grid.arrange(bacteria_rhizo_soil_wilcox,
    top = text_grob("Differentially abundant bacteria",
      size = 12, face = 2
    )
  ),
  device = "pdf"
)

# Saving 10.5 x 7.65 in image

# **********************************************************************--------
# *****  PHENOPTYPIC ANALYSIS ***** --------------------------------------------

# Import datasets --------------------------------------------------------------

fig1 <- read.csv(
  file = file.path(data_path, "FIG1.csv"), 
  header = TRUE)

fig1

leaf_area <- read.csv(
  file = file.path(data_path, "/AREA.csv"), 
  header = TRUE)

leaf_area

#plant_length <- read.csv(file = "datasets/LENGTH.csv", header = TRUE)
#plant_length

plant_height <- read.csv(
  file = file.path(data_path, "HEIGTH.csv"), 
  header = TRUE)
plant_height

dry_mass <- read.csv(
  file = file.path(data_path, "DW.csv"),
  header = TRUE)

dry_mass

#  >>>>> FIGURE 1 <<<<< -------------------------------------------------------

# >>>>> FIGURE 1B <<<<< --------------------------------------------------------
# leaf area --------------------------------------------------------------------

fig1 %>%
  filter(group == "FW (g)") %>%
  pivot_longer(-group, names_to = "plant", values_to = "value") %>%
  mutate(plant = as.factor(plant)) %>%
  t.test(value ~ plant, data = .)


Fig1B_Freshweight <-
  fig1 %>%
  filter(group == "FW (g)") %>%
  pivot_longer(-group, names_to = "plant", values_to = "value") %>%
  mutate(plant = as.factor(plant)) %>%
  ggplot(aes(x = plant, y = value, colour = plant)) +
  geom_jitter(position = position_jitter(0.4), size = 2, shape = 16) +
  stat_summary(
    fun = mean, fun.min = mean, fun.max = mean,
    geom = "crossbar", width = 0.8, color = "black"
  ) +
  stat_summary(
    geom = "text", angle = 0,
    label = c("a", "b"),
    fun = max, aes(y = 2.1), size = 4, color = "black"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 12, face = "bold", vjust = 0.5, hjust = 0.5),
    axis.title.y = element_text(size = 10, colour = "black", hjust = 0.5, vjust = 0.5),
    axis.text.x = element_text(size = 10, angle = 0, colour = "black", hjust = 0.5, vjust = 0.5),
    axis.text.y = element_text(size = 10, angle = 0, hjust = 0.5),
    legend.position = "none"
  ) +
  scale_y_continuous(n.breaks = 10, limits = c(0, 2.1)) +
  labs(title = "Fresh biomass (g)", x = NULL, y = NULL) +
  scale_color_manual(values = c("#825121", "#FFB400"))

Fig1B_Freshweight

# >>>>> FIGURE 1C <<<<< --------------------------------------------------------
# dry mass ---------------------------------------------------------------------

fig1 %>%
  filter(group == "DW (g)") %>%
  pivot_longer(-group, names_to = "plant", values_to = "value") %>%
  mutate(plant = as.factor(plant)) %>%
  t.test(value ~ plant, data = .)


Fig1C_DryMass <-
  fig1 %>%
  filter(group == "DW (g)") %>%
  pivot_longer(-group, names_to = "plant", values_to = "value") %>%
  mutate(plant = as.factor(plant)) %>%
  ggplot(aes(x = plant, y = value, colour = plant)) +
  geom_jitter(position = position_jitter(0.4), size = 2, shape = 16) +
  stat_summary(
    fun = mean, fun.min = mean, fun.max = mean,
    geom = "crossbar", width = 0.8, color = "black"
  ) +
  stat_summary(
    geom = "text", angle = 0,
    label = c("a", "b"),
    fun = max, aes(y = 0.06), size = 4, color = "black"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 12, face = "bold", vjust = 0.5, hjust = 0.5),
    axis.title.y = element_text(size = 10, colour = "black", hjust = 0.5, vjust = 0.5),
    axis.text.x = element_text(size = 10, angle = 0, colour = "black", hjust = 0.5, vjust = 0.5),
    axis.text.y = element_text(size = 10, angle = 0, hjust = 0.5),
    legend.position = "none"
  ) +
  scale_y_continuous(n.breaks = 6, limits = c(0, 0.06)) +
  labs(title = "Dry biomass (g)", x = NULL, y = NULL) +
  scale_color_manual(values = c("#825121", "#FFB400"))

Fig1C_DryMass

# >>>>> FIGURE 1D <<<<< --------------------------------------------------------
# root length ------------------------------------------------------------------

fig1 %>%
  filter(group == "Total root length (mm)") %>%
  pivot_longer(-group, names_to = "plant", values_to = "value") %>%
  mutate(plant = as.factor(plant)) %>%
  t.test(value ~ plant, data = .)


Fig1D_RootLength <-
  fig1 %>%
  filter(group == "Total root length (mm)") %>%
  pivot_longer(-group, names_to = "plant", values_to = "value") %>%
  mutate(plant = as.factor(plant)) %>%
  ggplot(aes(x = plant, y = value / 100, colour = plant)) +
  geom_jitter(position = position_jitter(0.4), size = 2, shape = 16) +
  stat_summary(
    fun = mean, fun.min = mean, fun.max = mean,
    geom = "crossbar", width = 0.8, color = "black"
  ) +
  stat_summary(
    geom = "text", angle = 0,
    label = c("a", "b"),
    fun = max, aes(y = 175), size = 4, color = "black"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 12, face = "bold", vjust = 0.5, hjust = 0.5),
    axis.title.y = element_text(size = 10, colour = "black", hjust = 0.5, vjust = 0.5),
    axis.text.x = element_text(size = 10, angle = 0, colour = "black", hjust = 0.5, vjust = 0.5),
    axis.text.y = element_text(size = 10, angle = 0, hjust = 0.5),
    legend.position = "none"
  ) +
  scale_y_continuous(n.breaks = 8, limits = c(0, 180)) +
  labs(title = "Total root length (cm)", x = NULL, y = NULL) +
  scale_color_manual(values = c("#825121", "#FFB400"))

Fig1D_RootLength

# >>>>> FIGURE 1B <<<<< --------------------------------------------------------
# root depth --------------------------------------------------------------------

fig1 %>%
  filter(group == "Depth (mm)") %>%
  pivot_longer(-group, names_to = "plant", values_to = "value") %>%
  mutate(plant = as.factor(plant)) %>%
  t.test(value ~ plant, data = .)


Fig1E_RootDepth <-
  fig1 %>%
  filter(group == "Depth (mm)") %>%
  pivot_longer(-group, names_to = "plant", values_to = "value") %>%
  mutate(plant = as.factor(plant)) %>%
  ggplot(aes(x = plant, y = value / 10, colour = plant)) +
  geom_jitter(position = position_jitter(0.4), size = 2, shape = 16) +
  stat_summary(
    fun = mean, fun.min = mean, fun.max = mean,
    geom = "crossbar", width = 0.8, color = "black"
  ) +
  stat_summary(
    geom = "text", angle = 0,
    label = c("a", "b"),
    fun = max, aes(y = 32), size = 4, color = "black"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 12, face = "bold", vjust = 0.5, hjust = 0.5),
    axis.title.y = element_text(size = 10, colour = "black", hjust = 0.5, vjust = 0.5),
    axis.text.x = element_text(size = 10, angle = 0, colour = "black", hjust = 0.5, vjust = 0.5),
    axis.text.y = element_text(size = 10, angle = 0, hjust = 0.5),
    legend.position = "none"
  ) +
  scale_y_continuous(n.breaks = 8, limits = c(0, 32)) +
  labs(title = "Root depth (cm)", x = NULL, y = NULL) +
  scale_color_manual(values = c("#825121", "#FFB400"))

Fig1E_RootDepth

# Multiplot --------------------------------------------------------------------

Fig1BCDE_multi_pheno <-
  ggarrange(
    Fig1B_Freshweight,
    Fig1C_DryMass,
    Fig1D_RootLength,
    Fig1E_RootDepth,
    ncol = 2,
    nrow = 2,
    align = "hv",
    labels = c("A", "B", "C", "D")
  )

Fig1BCDE_multi_pheno

ggsave(path = results_path,
       filename = "Fig1_BCDE_for_manuel.pdf",
       plot = Fig1BCDE_multi_pheno)

#  >>>>> FIGURE 2 <<<<< --------------------------------------------------------

# kurskal test
leaf_area %>%
  mutate(SampleID = paste("Sample", 1:10, sep = "_")) %>%
  pivot_longer(-SampleID, names_to = "plant", values_to = "area") %>%
  kruskal.test(area ~ plant, data = .)

leaf_area %>%
  mutate(SampleID = paste("Sample", 1:10, sep = "_")) %>%
  pivot_longer(-SampleID, names_to = "plant", values_to = "area") %>%
  aov(area ~ plant, data = .) %>%
  summary()

# Tukey test
leaf_area_tukey <-
  leaf_area %>%
  mutate(SampleID = paste("Sample", 1:10, sep = "_")) %>%
  pivot_longer(-SampleID, names_to = "plant", values_to = "area") %>%
  aov(area ~ plant, data = .) %>%
  HSD.test(trt = "plant",
           group = TRUE,
           console = TRUE)

leaf_area_tukey$groups %>%
  rownames_to_column("plant") %>%
  arrange(plant)

# >>>>> FIGURE 2D <<<<< --------------------------------------------------------

Fig2D_leaf_area_plot <-
  leaf_area %>%
  mutate(SampleID = paste("Sample", 1:10, sep = "_")) %>%
  pivot_longer(-SampleID, names_to = "plant", values_to = "area") %>%
  ggplot(aes(x = plant, y = area, colour = plant)) +
  geom_jitter(position = position_jitter(0.4), size = 2, shape = 16) +
  stat_summary(
    fun = mean, fun.min = mean, fun.max = mean,
    geom = "crossbar", width = 0.8, color = "black"
  ) +
  stat_summary(
    geom = "text", angle = 0,
    label = c("a", "b", "ab", "c"),
    fun = max, aes(y = 250), size = 4, color = "black"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 12, face = "bold", vjust = 0.5, hjust = 0.5),
    axis.title.y = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(size = 10, angle = 33, colour = "black", hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10, angle = 0, hjust = 0.5),
    legend.position = "none"
  ) +
  scale_y_continuous(n.breaks = 10, limits = c(0, 250)) +
  labs(title = expression(bold(paste("Leaf area (", cm^2, ")"))), y = NULL, x = NULL) +
  scale_color_manual(values = c("#825121", "#CC2D35", "#FFB400", "#00A6ED"))

Fig2D_leaf_area_plot

# >>>>> FIGURE 2E <<<<< --------------------------------------------------------

# plant_height -----------------------------------------------------------------
plant_height %>%
  mutate(SampleID = paste("Sample", 1:10, sep = "_")) %>%
  pivot_longer(-SampleID, names_to = "plant", values_to = "height") %>%
  kruskal.test(height ~ plant, data = .)

plant_height %>%
  mutate(SampleID = paste("Sample", 1:10, sep = "_")) %>%
  pivot_longer(-SampleID, names_to = "plant", values_to = "height") %>%
  aov(height ~ plant, data = .) %>%
  summary()

# Tukey test
leaf_height_tukey <-
  plant_height %>%
  mutate(SampleID = paste("Sample", 1:10, sep = "_")) %>%
  pivot_longer(-SampleID, names_to = "plant", values_to = "height") %>%
  aov(height ~ plant, data = .) %>%
  HSD.test(trt = "plant", 
           group = TRUE,
           console = TRUE)

leaf_height_tukey$groups %>%
  rownames_to_column("plant") %>%
  arrange(plant)

Fig2E_leaf_height_plot <-
  plant_height %>%
  mutate(SampleID = paste("Sample", 1:10, sep = "_")) %>%
  pivot_longer(-SampleID, names_to = "plant", values_to = "height") %>%
  ggplot(aes(x = plant, y = height, colour = plant)) +
  geom_jitter(position = position_jitter(0.4), size = 2, shape = 16) +
  stat_summary(
    fun = mean, fun.min = mean, fun.max = mean,
    geom = "crossbar", width = 0.8, color = "black"
  ) +
  stat_summary(
    geom = "text", angle = 0,
    label = c("a", "b", "c", "d"),
    fun = max, aes(y = 30), size = 4, color = "black"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 12, face = "bold", vjust = 0.5, hjust = 0.5),
    axis.title.y = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(size = 10, angle = 33, colour = "black", hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10, angle = 0, hjust = 0.5),
    legend.position = "none"
  ) +
  scale_y_continuous(n.breaks = 10, limits = c(0, 30)) +
  labs(
    title = "Plant height (cm)",
    y = NULL, x = NULL
  ) +
  scale_color_manual(values = c("#825121", "#CC2D35", "#FFB400", "#00A6ED"))


Fig2E_leaf_height_plot

# >>>>> FIGURE 2C <<<<< --------------------------------------------------------

# Shoot dry mass ---------------------------------------------------------------
dry_mass %>%
  mutate(SampleID = paste("Sample", 1:10, sep = "_")) %>%
  pivot_longer(-SampleID, names_to = "plant", values_to = "dry_mass") %>%
  kruskal.test(dry_mass ~ plant, data = .)

dry_mass %>%
  mutate(SampleID = paste("Sample", 1:10, sep = "_")) %>%
  pivot_longer(-SampleID, names_to = "plant", values_to = "dry_mass") %>%
  aov(dry_mass ~ plant, data = .) %>%
  summary()

# Tukey test
leaf_mass_tukey <-
  dry_mass %>%
  mutate(SampleID = paste("Sample", 1:10, sep = "_")) %>%
  pivot_longer(-SampleID, names_to = "plant", values_to = "dry_mass") %>%
  aov(dry_mass ~ plant, data = .) %>%
  HSD.test(trt = "plant", 
           group = TRUE,
           console = TRUE)

leaf_mass_tukey$groups %>%
  rownames_to_column("plant") %>%
  arrange(plant)

Fig2F_plant_mass_plot <-
  dry_mass %>%
  mutate(SampleID = paste("Sample", 1:10, sep = "_")) %>%
  pivot_longer(-SampleID, names_to = "plant", values_to = "dry_mass") %>%
  ggplot(aes(x = plant, y = dry_mass, colour = plant)) +
  geom_jitter(position = position_jitter(0.4), size = 2, shape = 16) +
  stat_summary(
    fun = mean, fun.min = mean, fun.max = mean,
    geom = "crossbar", width = 0.8, color = "black"
  ) +
  stat_summary(
    geom = "text", angle = 0,
    label = c("a", "b", "b", "c"),
    fun = max, aes(y = 6.5), size = 4, color = "black"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 12, face = "bold", vjust = 0.5, hjust = 0.5),
    axis.title.y = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(size = 10, angle = 33, colour = "black", hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10, angle = 0, hjust = 0.5),
    legend.position = "none"
  ) +
  scale_y_continuous(n.breaks = 10, limits = c(0, 6.5)) +
  labs(title = "Plant dry mass (g)", y = NULL, x = NULL) +
  scale_color_manual(values = c("#825121", "#CC2D35", "#FFB400", "#00A6ED"))

Fig2F_plant_mass_plot


# Multiplot --------------------------------------------------------------------

Fig2DEF_pheno_all <-
  ggarrange(
    Fig2D_leaf_area_plot,
    Fig2E_leaf_height_plot,
    Fig2F_plant_mass_plot,
    ncol = 3,
    nrow = 1,
    align = "hv",
    labels = c("A", "B", "C")
  )

Fig2DEF_pheno_all

ggsave(path = results_path,
       filename = "Fig2_DEF_for_manuel.pdf", 
       plot = Fig2DEF_pheno_all)
