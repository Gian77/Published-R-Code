# ************************************************ -----------------------------
# Manuscript:     "Mycorrhizal competition release and microbial dynamics in native 
#                  and non-native Tuber melanosporum habitats."
# Authors:         Gian Maria Niccolo Benucci, Sergi Garcia-Barreda, Sergio Sánchez,
#                  Pedro Marco, Ana De Miguel, Francois Le Tacon, Giorgio Marozzi,
#                  Leonardo Baciarelli Failini, Harry Eslick, Todd Elliott, Aurelie
#                  Deveau, Claude Murat, Domizia Donnini, and Gregory Bonito
# Code Developer:  Gian MN Benucci
# Affiliation:     Michigan State University, ...
# Journal:         Applied and Environmental Microbiology
#
# Citation         Benucci et al. 2025
#
# Date:            October 15, 2025
# ************************************************ -----------------------------

# ***** ENVIRONEMNT SETUP ***** ------------------------------------------------

# Options ----------------------------------------------------------------------
options(scipen = 9999, pillar.sigfig = 6, digits = 6, max.print = 20000000)
# rm(list = ls())

# Check the lib paths ----------------------------------------------------------
.libPaths()

# Session Info -----------------------------------------------------------------
sessionInfo()

# libraries --------------------------------------------------------------------
pacman::p_load(
  conflicted,
  styler,
  phyloseq,
  ape,
  vegan,
  AICcPermanova,
  mvabund,
  tidyverse,
  data.table,
  MASS,
  ggpubr,
  ggtext,
  magrittr,
  ggrepel,
  scales,
  cowplot,
  gridExtra,
  formattable,
  fungaltraits,
  ggcorrplot,
  indicspecies,
  Boruta,
  ggvenn,
  ppcor,
  randomForest,
  RColorBrewer,
  mediation,
  lmPerm,
  rstatix,
  DHARMa,
  lmerTest,
  glmmTMB,
  emmeans,
  multcomp,
  knitr,
  BRCore,
  install = FALSE
)


# Solve known conflicts
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::combine)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(base::intersect)
conflicts_prefer(ggplot2::annotate)
conflicts_prefer(phyloseq::tax_glom)
conflicts_prefer(dplyr::combine)
conflicts_prefer(dplyr::last)
conflicts_prefer(ggpubr::get_legend)
conflicts_prefer(lmerTest::lmer)
conflicts_prefer(dplyr::slice)

# ************************************************ -----------------------------
# ***** PATHS ***** ------------------------------------------------------------
# datasets ---------------------------------------------------------------------

project_dir <-
  ("/home/gian/Dropbox/6_PROJETCS/2025_CompetitionRelease_Tuber_EuropeAustralia_AEMgithub/")

project_dir

# ************************************************ -----------------------------
# ***** IMPORT ***** -----------------------------------------------------------

# Import datasets --------------------------------------------------------------
physeq_ITS <- import_biom(file.path(project_dir, "datasets/otu_table_ITS_uparse.biom"))
sample_data(physeq_ITS) <- import_qiime_sample_data(file.path(project_dir, "datasets/mapping_ITS_GPS_soil.txt"))
colnames(tax_table(physeq_ITS)) <- c(
  k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species"
)
otus_ITS <- Biostrings::readDNAStringSet(
  file.path(project_dir, "datasets/otus_ITS_uparse.fasta"),
  format = "fasta",
  seek.first.rec = TRUE,
  use.names = TRUE
)
physeq_ITS <- merge_phyloseq(physeq_ITS, otus_ITS)
physeq_ITS

str(physeq_ITS)
sample_sums(physeq_ITS)
head(physeq_ITS@tax_table)

physeq_16S <- import_biom(file.path(project_dir, "datasets/otu_table_16s_uparse.biom"))
sample_data(physeq_16S) <- import_qiime_sample_data(file.path(project_dir, "datasets/mapping_16S_GPS_soil.txt"))
colnames(tax_table(physeq_16S)) <- c(
  k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species"
)
otus_16S <- Biostrings::readDNAStringSet(
  file.path(project_dir, "datasets/otus_16s_uparse.fasta"),
  format = "fasta",
  seek.first.rec = TRUE,
  use.names = TRUE
)
physeq_16S <- merge_phyloseq(physeq_16S, otus_16S)
physeq_16S

str(physeq_16S)
sample_sums(physeq_16S)
head(physeq_16S@tax_table)

# cleaning ITS taxonomy from extra characters ----------------------------------
tax_table(physeq_ITS)[, "Kingdom"] <- gsub("PMI", "", tax_table(physeq_ITS)[, "Kingdom"])
tax_table(physeq_ITS)[, "Kingdom"] <- gsub("NVP", "", tax_table(physeq_ITS)[, "Kingdom"])
tax_table(physeq_ITS)[, "Kingdom"] <- gsub("k__", "", tax_table(physeq_ITS)[, "Kingdom"])
tax_table(physeq_ITS)[, "Phylum"] <- gsub("p__", "", tax_table(physeq_ITS)[, "Phylum"])
tax_table(physeq_ITS)[, "Class"] <- gsub("c__", "", tax_table(physeq_ITS)[, "Class"])
tax_table(physeq_ITS)[, "Order"] <- gsub("o__", "", tax_table(physeq_ITS)[, "Order"])
tax_table(physeq_ITS)[, "Family"] <- gsub("f__", "", tax_table(physeq_ITS)[, "Family"])
tax_table(physeq_ITS)[, "Genus"] <- gsub("g__", "", tax_table(physeq_ITS)[, "Genus"])
tax_table(physeq_ITS)[, "Species"] <- gsub("s__", "", tax_table(physeq_ITS)[, "Species"])
head(tax_table(physeq_ITS))

tax_table(physeq_16S)[, "Kingdom"] <- gsub("k__", "", tax_table(physeq_16S)[, "Kingdom"])
tax_table(physeq_16S)[, "Phylum"] <- gsub("p__", "", tax_table(physeq_16S)[, "Phylum"])
tax_table(physeq_16S)[, "Class"] <- gsub("c__", "", tax_table(physeq_16S)[, "Class"])
tax_table(physeq_16S)[, "Order"] <- gsub("o__", "", tax_table(physeq_16S)[, "Order"])
tax_table(physeq_16S)[, "Family"] <- gsub("f__", "", tax_table(physeq_16S)[, "Family"])
tax_table(physeq_16S)[, "Genus"] <- gsub("g__", "", tax_table(physeq_16S)[, "Genus"])
tax_table(physeq_16S)[, "Species"] <- gsub("s__", "", tax_table(physeq_16S)[, "Species"])
head(tax_table(physeq_16S))

# filtering out non-fungal OTUs ------------------------------------------------
# and Mitochondria/Chloroplast sequences.
physeq_ITS <- subset_taxa(physeq_ITS, Kingdom == "Fungi")
physeq_ITS

# pay attention here if you want to keep Archaea  ------------------------------
physeq_16S <- subset_taxa(physeq_16S, Kingdom == "Bacteria")
physeq_16S <- subset_taxa(physeq_16S, Phylum!="Chloroplast")
physeq_16S <- subset_taxa(physeq_16S, Class!="Chloroplast")
physeq_16S <- subset_taxa(physeq_16S, Order!="Chloroplast")
physeq_16S <- subset_taxa(physeq_16S, Family!="Chloroplast")
physeq_16S <- subset_taxa(physeq_16S, Genus!="Chloroplast")
physeq_16S <- subset_taxa(physeq_16S, Phylum!="chloroplast")
physeq_16S <- subset_taxa(physeq_16S, Class!="chloroplast")
physeq_16S <- subset_taxa(physeq_16S, Order!="chloroplast")
physeq_16S <- subset_taxa(physeq_16S, Family!="chloroplast")
physeq_16S <- subset_taxa(physeq_16S, Genus!="chloroplast")
physeq_16S <- subset_taxa(physeq_16S, Phylum!="Mitochondria")
physeq_16S <- subset_taxa(physeq_16S, Class!="Mitochondria")
physeq_16S <- subset_taxa(physeq_16S, Order!="Mitochondria")
physeq_16S <- subset_taxa(physeq_16S, Family!="Mitochondria")
physeq_16S <- subset_taxa(physeq_16S, Genus!="Mitochondria")
physeq_16S <- subset_taxa(physeq_16S, Phylum!="mitochondria")
physeq_16S <- subset_taxa(physeq_16S, Class!="mitochondria")
physeq_16S <- subset_taxa(physeq_16S, Order!="mitochondria")
physeq_16S <- subset_taxa(physeq_16S, Family!="mitochondria")
physeq_16S <- subset_taxa(physeq_16S, Genus!="mitochondria")
physeq_16S

# Adding "unclassified" to empty fields ----------------------------------------
tax_table(physeq_ITS)[tax_table(physeq_ITS)==""]<- NA
tax_table(physeq_ITS)[tax_table(physeq_ITS)=="unidentified"]<- NA
tax_table(physeq_ITS)[tax_table(physeq_ITS)=="Unidentified"]<- NA
tax_table(physeq_ITS)[is.na(tax_table(physeq_ITS))]<-"Unclassified"

tax_table(physeq_16S)[tax_table(physeq_16S)==""]<- NA
tax_table(physeq_16S)[tax_table(physeq_16S)=="unidentified"]<- NA
tax_table(physeq_16S)[tax_table(physeq_16S)=="Unidentified"]<- NA
tax_table(physeq_16S)[is.na(tax_table(physeq_16S))]<-"Unclassified"

# Change OTUs names ------------------------------------------------------------
otu_names_ITS <- 1:nrow(as.data.frame(otu_table(physeq_ITS)))
otu_names_ITS <- paste("FOTU", otu_names_ITS, sep = "_")
otu_names_ITS
taxa_names(physeq_ITS) <- paste(otu_names_ITS)
head(otu_table(physeq_ITS))
head(tax_table(physeq_ITS))
sample_data(physeq_ITS)

otu_names_16s <- 1:nrow(as.data.frame(otu_table(physeq_16S)))
otu_names_16s <- paste("BOTU", otu_names_16s, sep = "_")
otu_names_16s
taxa_names(physeq_16S) <- paste(otu_names_16s)
head(otu_table(physeq_16S))
head(tax_table(physeq_16S))
sample_data(physeq_16S)

# Save final clean objects -----------------------------------------------------
saveRDS(object = physeq_ITS, file = file.path(project_dir, "github/phyloseq_ITS.RDS"))
saveRDS(object = physeq_16S, file = file.path(project_dir, "github/phyloseq_16S.RDS"))

# ************************************************ -----------------------------
# ***** DECONTAMINATION ***** --------------------------------------------------
# decontamination from a phyloseq object ---------------------------------------

# Fungi ------------------------------------------------------------------------
table(physeq_ITS@sam_data$brule)

sample_data(physeq_ITS)$is.neg <-
  sample_data(physeq_ITS)$brule == "Control"

contam_ITS <-
  decontam::isContaminant(physeq_ITS,
                          method = "prevalence",
                          neg = "is.neg",
                          threshold = 0.5)

table(contam_ITS$contaminant)

contam_ITS %>% 
  filter(contaminant == TRUE) %>%
  rownames_to_column("otu_id") %>%
  left_join(
    tax_table(physeq_ITS) %>%
      as.matrix() %>%
      as.data.frame() %>%
      rownames_to_column("otu_id"),
    by = "otu_id") 

# remove taxa by OTU name
remove_taxa <- function(physeq, badTaxa) {
  allTaxa <- taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}


# Removing contaminants with no Phylum classification
physeq_ITS_clean <-
  remove_taxa(physeq_ITS, rownames(subset(
    contam_ITS %>% 
      filter(contaminant == TRUE) %>%
      rownames_to_column("otu_id") %>%
      left_join(
        tax_table(physeq_ITS) %>%
          as.matrix() %>%
          as.data.frame() %>%
          rownames_to_column("otu_id"),
        by = "otu_id") %>%
      filter(Order %in% "Unclassified"), # Unclassified at Order,
    contaminant %in% c("TRUE")
  ))) %>%
  subset_samples(!is.neg %in% TRUE) %>% # remove the control samples
  prune_taxa(taxa_sums(x = .) > 0, x = .) # make sure there aren't OTUs that are 0

physeq_ITS_clean
physeq_ITS_clean@sam_data

# Bacteria ---------------------------------------------------------------------
table(physeq_16S@sam_data$brule)

sample_data(physeq_16S)$is.neg <-
  sample_data(physeq_16S)$brule == "Control"

contam_16S <-
  decontam::isContaminant(physeq_16S,
                          method = "prevalence",
                          neg = "is.neg",
                          threshold = 0.5)

table(contam_16S$contaminant)

contam_16S %>% 
  filter(contaminant == TRUE) %>%
  rownames_to_column("otu_id") %>%
  left_join(
    tax_table(physeq_16S) %>%
      as.matrix() %>%
      as.data.frame() %>%
      rownames_to_column("otu_id"),
    by = "otu_id") 


# Removing contaminants with no Phylum classification
physeq_16S_clean <-
  remove_taxa(physeq_16S, rownames(subset(
    contam_16S %>% 
      filter(contaminant == TRUE) %>%
      rownames_to_column("otu_id") %>%
      left_join(
        tax_table(physeq_16S) %>%
          as.matrix() %>%
          as.data.frame() %>%
          rownames_to_column("otu_id"),
        by = "otu_id") %>%
      filter(Order %in% "Unclassified"), # Unclassified at Order,
    contaminant %in% c("TRUE")
  ))) %>%
  subset_samples(!is.neg %in% TRUE) %>% # remove the control samples
  prune_taxa(taxa_sums(x = .) > 0, x = .) # make sure there aren't OTUs that are 0

physeq_16S_clean


# ************************************************ -----------------------------
# ***** CLEAN PHYLOSEQ OBJECTS ***** -------------------------------------------
physeq_ITS_clean
physeq_16S_clean

# ***** SEPARATE DATASETS ***** ------------------------------------------------

# SOIL -------------------------------------------------------------------------
physeq_fungi <- 
  phyloseq(physeq_ITS_clean@otu_table, 
           sample_data(
             read.delim(file.path(project_dir, "datasets/mapping_ITS_GPS_soil.txt"), 
                        header = TRUE, sep = "\t", row.names=1) %>% 
             filter(even_sample %in% "tosample") %>% 
             select("continent","state","site", "site_code", "brule", "management", "latitude", "longitude",
                    "pH", "Olsen.P", "K", "Ca", "Mg", "Zn", "Mn", "Cu", "Fe", "NO3", "NH4") %>% 
               rename(P = Olsen.P)),
           physeq_ITS_clean@tax_table,
           physeq_ITS_clean@refseq) %>% 
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.)
           
physeq_fungi

physeq_bact <- 
  phyloseq(physeq_16S_clean@otu_table, 
           sample_data(
             read.delim(file.path(project_dir, "datasets/mapping_16S_GPS_soil.txt"), 
                        header = TRUE, sep = "\t", row.names=1) %>% 
               filter(even_sample %in% "tosample") %>% 
               select("continent","country", "site", "site_code", "brule", "management", "latitude", "longitude",
                      "pH", "Olsen.P", "K", "Ca", "Mg", "Zn", "Mn", "Cu", "Fe", "NO3", "NH4") %>% 
               rename(P = Olsen.P, state = country)),
           physeq_16S_clean@tax_table,
           physeq_16S_clean@refseq) %>% 
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.)

physeq_bact

# Rarefaction ------------------------------------------------------------------
library(BRCore)

add_rarefaction_metrics(data = physeq_fungi) %>% 
  plot_rarefaction_metrics()

physeq_fungi_rare <-
  do_phyloseq(
    physeq = physeq_fungi,
    otu_rare = multi_rarefy(
      physeq = physeq_fungi,
      depth_level = 3818,
      num_iter = 100,
      threads = 8,
      set_seed = 2635
    )
  ) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)


add_rarefaction_metrics(data = physeq_bact) %>% 
  plot_rarefaction_metrics()

physeq_bact_rare <-
  do_phyloseq(
    physeq = physeq_bact,
    otu_rare = multi_rarefy(
      physeq = physeq_bact,
      depth_level = 4437,
      num_iter = 100,
      threads = 8,
      set_seed = 9873
    )
  ) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_bact_rare

# Mycorrhizal fungi ------------------------------------------------------------
# Extracting ECM fungi from the rarefied dataset
fungal_traits()
as.matrix(tax_table(physeq_fungi_rare)) %>% head()

ecm_guilds <- c("Ectomycorrhizal",
                "Ectomycorrhizal-Fungal Parasite",
                "Ectomycorrhizal-Undefined Saprotroph",
                "Ectomycorrhizal-Endophyte-Ericoid Mycorrhizal-Litter Saprotroph-Orchid Mycorrhizal",
                "Ectomycorrhizal-Orchid Mycorrhizal-Root Associated Biotroph",
                "Ectomycorrhizal-Wood Saprotroph",
                "Dung Saprotroph-Ectomycorrhizal")

ecm_taxonomy <-
  left_join(
    x = as.data.frame(as.matrix(tax_table(physeq_fungi_rare))) %>%
      rownames_to_column("OTU_ID"),
    y = as.data.frame(fungal_traits()),
    by = "Genus",
    unmatched = "drop",
    relationship = "many-to-many"
  ) %>%
  filter(guild_fg %in% ecm_guilds) %>%
  dplyr::select(
    OTU_ID,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    confidence_fg,
    guild_fg,
    trophic_mode_fg,
    growth_form_fg,
    ifungorum_number
  ) %>%
  distinct(OTU_ID, .keep_all = TRUE) %>%
  as_tibble() %>% #filter to only "Highly Probable" and "Probable"
  filter(confidence_fg %in% c("Probable", "Highly Probable")) %>% 
  mutate(
    ifungorum_number = ifelse(is.na(ifungorum_number), "missing", ifungorum_number),
    growth_form_fg = ifelse(is.null(growth_form_fg), "missing", growth_form_fg)
  )

unique(ecm_taxonomy$guild_fg)
unique(ecm_taxonomy$OTU_ID)

# Generate phyloseq object
physeq_ecm <-
  merge_phyloseq(
    otu_table(physeq_fungi_rare, taxa_are_rows = TRUE),
    sample_data(physeq_fungi_rare),
    tax_table(as.matrix(
      ecm_taxonomy %>% column_to_rownames("OTU_ID")
    )),
    refseq(physeq_fungi_rare)
  ) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_ecm

any(rowSums(otu_table(physeq_ecm)) == 0)
any(colSums(otu_table(physeq_ecm)) == 0)
sort(rowSums(otu_table(physeq_ecm)))
sort(colSums(otu_table(physeq_ecm)))

# GLEBA ------------------------------------------------------------------------

sample_data(physeq_16S_clean)$sample_name <- 
  sample_names(physeq_16S_clean)

physeq_gleba <-
  physeq_16S_clean %>%
  subset_samples(brule %in% "gleba") %>% # Remove NOT T. melanosporum.
  subset_samples(!site %in% c("Yarra", "Needles", "Launceston")) %>%
  subset_samples(
    !sample_name %in% c(
      "EU_16S_Tmel_12","EU_16S_Tmel_36","EU_16S_Tmel_24",
      "EU_16S_Tmel_40","EU_16S_Tmel_48")
    ) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_gleba

# rarefaction 
add_rarefaction_metrics(data = physeq_gleba) %>% 
  plot_rarefaction_metrics()

physeq_gleba_rare <-
  do_phyloseq(
    physeq = physeq_gleba,
    otu_rare = multi_rarefy(
      physeq = physeq_gleba,
      depth_level = 2000,
      num_iter = 100,
      threads = 8,
      set_seed = 8256
    )
  ) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_gleba_rare

# OUTGROUP ---------------------------------------------------------------------

physeq_bact_out <-
  physeq_16S_clean %>%
  subset_samples(brule %in% c("outgroup")) %>%
  subset_samples(!place %in% "Nancy_Claude") %>% # We did not include this orchard
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_bact_out
physeq_bact_out@sam_data

# rarefaction 
add_rarefaction_metrics(data = physeq_bact_out) %>% 
  plot_rarefaction_metrics()

sample_sums(physeq_bact_out) %>% sort()

physeq_bact_out_rare <-
  do_phyloseq(
    physeq = physeq_bact_out,
    otu_rare = multi_rarefy(
      physeq = physeq_bact_out,
      depth_level = 19333,
      num_iter = 100,
      threads = 8,
      set_seed = 4299
    )
  ) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
  prune_samples(sample_sums(x = .) > 0, x = .)

physeq_bact_out_rare





























