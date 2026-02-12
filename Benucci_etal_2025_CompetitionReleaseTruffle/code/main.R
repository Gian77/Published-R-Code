#************************************************************************-------
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
  #AICcPermanova,
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
  #BRCore,
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

# load data for plotting
permanova_summary <- 
  readRDS(file = file.path(project_dir, "github/permanova_summary.RDS"))

# ************************************************ -----------------------------
# ***** BETA DIVERSITY  ***** --------------------------------------------------

# Fungi ------------------------------------------------------------------------
physeq_fungi_all <- AlphaMetrics(physeq_fungi_rare)

sample_data(physeq_fungi_all)$sample_id <- 
  sample_names(physeq_fungi_all)

# rearrange factors
sample_data(physeq_fungi_all) <-
  sample_data(
    as.data.frame(as.matrix(physeq_fungi_all@sam_data)) %>% 
      mutate(site = recode(site,
                           S.Demetrio = "San Demetrio",
                           Uzes = "Nimes",
                           Manjimup2 = "Jardee", 
                           Chateauvert = "Grignan",
                           Angouleme = "Cognac",
                           `Espira-de-l-Agly` = "Perpignan",
                           Parnans = "Romans-Sur-Isere")) %>% 
      mutate(TreeID = paste0(site_code, str_extract(sample_id, "\\d+"))) %>% 
      mutate(continent_brule = interaction(continent, brule, sep = "_")) %>% 
      mutate(site_brule = interaction(site, brule, sep = "_")) %>% 
      select(sample_id, TreeID, 
             continent, continent_brule, 
             site, site_brule, 
             brule, management,
             hill_0, hill_1, hill_2) %>% 
      mutate(hill_0 = as.numeric(hill_0), 
             hill_1 = as.numeric(hill_1),
             hill_2 = as.numeric(hill_2)) %>%
      mutate(TreeID = recode(TreeID, EP219 = "EP290", EP288 = "EP298")) %>% 
      janitor::clean_names() %>% 
      mutate(site = fct_relevel(site, 
                                "Yarra", "Wattles", "Launceston", "Needles", "Mole",
                                "Camberra", "Warri", "Braidwood", "Manjimup",
                                "Pemberton", "Jardee", 
                                "Cuneo", "San Demetrio", "Spoleto", "Norcia",
                                "Cognac", "Grignan", "Perpignan", "Nimes", "Romans-Sur-Isere",
                                "Albentosa", "Moncayo", "Zuniga", "Acedo")) %>% 
      mutate(site_brule = paste(site, brule, sep = " ")) %>% 
      mutate(site_brule = recode_factor(
        site_brule,
        `Yarra inside` = "Ya-In", `Yarra outside` = "Ya-Out",
        `Wattles inside` = "Wa-In", `Wattles outside` = "Wa-Out",
        `Launceston inside` = "La-In", `Launceston outside` = "La-Out",
        `Needles inside` = "Ne-In", `Needles outside` = "Ne-Out",
        `Mole inside` = "Mol-In", `Mole outside` = "Mol-Out",
        `Camberra inside` = "Ca-In", `Camberra outside` = "Ca-Out",
        `Warri inside` = "Wi-In", `Warri outside` = "Wi-Out",
        `Braidwood inside` = "Br-In", `Braidwood outside` = "Br-Out",
        `Manjimup inside` = "Ma-In", `Manjimup outside` = "Ma-Out",
        `Pemberton inside` = "Pem-In", `Pemberton outside` = "Pem-Out",
        `Jardee inside` = "Ja-In", `Jardee outside` = "Ja-Out",
        `Cuneo inside` = "Cu-In", `Cuneo outside` = "Cu-Out",
        `San Demetrio inside` = "SD-In", `San Demetrio outside` = "SD-Out",
        `Spoleto inside` = "Sp-In", `Spoleto outside` = "Sp-Out",
        `Norcia inside` = "No-In", `Norcia outside` = "No-Out",
        `Cognac inside` = "Co-In", `Cognac outside` = "Co-Out",
        `Grignan inside` = "Gr-In", `Grignan outside` = "Gr-Out",
        `Perpignan inside` = "Per-In", `Perpignan outside` = "Per-Out",
        `Nimes inside` = "Ni-In", `Nimes outside` = "Ni-Out",
        `Romans-Sur-Isere inside` = "Ro-In", `Romans-Sur-Isere outside` = "Ro-Out",
        `Albentosa inside` = "Al-In", `Albentosa outside` = "Al-Out",
        `Moncayo inside` = "Mon-In", `Moncayo outside` = "Mon-Out",
        `Zuniga inside` = "Zu-In", `Zuniga outside` = "Zu-Out",
        `Acedo inside` = "Ac-In", `Acedo outside` = "Ac-Out")) %>% 
      mutate_at(c("hill_0", "hill_1", "hill_2"), as.numeric) %>% 
      filter(management %in% c("cultivated")) %>%
      droplevels()
  ) 

head(physeq_fungi_all@sam_data)

# Ectomycorrhizal fungi --------------------------------------------------------
physeq_ecm_all_all <- AlphaMetrics(physeq_ecm_all)

sample_data(physeq_ecm_all_all)$sample_id <- 
  sample_names(physeq_ecm_all_all)

# rearrange factors
sample_data(physeq_ecm_all_all) <-
  sample_data(
    as.data.frame(as.matrix(physeq_ecm_all_all@sam_data)) %>% 
      mutate(site = recode(site,
                           S.Demetrio = "San Demetrio",
                           Uzes = "Nimes",
                           Manjimup2 = "Jardee", 
                           Chateauvert = "Grignan",
                           Angouleme = "Cognac",
                           `Espira-de-l-Agly` = "Perpignan",
                           Parnans = "Romans-Sur-Isere")) %>% 
      mutate(TreeID = paste0(site_code, str_extract(sample_id, "\\d+"))) %>% 
      mutate(continent_brule = interaction(continent, brule, sep = "_")) %>% 
      mutate(site_brule = interaction(site, brule, sep = "_")) %>% 
      select(sample_id, TreeID, 
             continent, continent_brule, 
             site, site_brule, 
             brule, management,
             hill_0, hill_1, hill_2) %>% 
      mutate(hill_0 = as.numeric(hill_0), 
             hill_1 = as.numeric(hill_1),
             hill_2 = as.numeric(hill_2)) %>%
      mutate(TreeID = recode(TreeID, EP219 = "EP290", EP288 = "EP298")) %>% 
      janitor::clean_names() %>% 
      mutate(site = fct_relevel(site, 
                                "Yarra", "Wattles", "Launceston", "Needles", "Mole",
                                "Camberra", "Warri", "Braidwood", "Manjimup",
                                "Pemberton", "Jardee", 
                                "Cuneo", "San Demetrio", "Spoleto", "Norcia",
                                "Cognac", "Grignan", "Perpignan", "Nimes", "Romans-Sur-Isere",
                                "Albentosa", "Moncayo", "Zuniga", "Acedo")) %>% 
      mutate(site_brule = paste(site, brule, sep = " ")) %>% 
      mutate(site_brule = recode_factor(
        site_brule,
        `Yarra inside` = "Ya-In", `Yarra outside` = "Ya-Out",
        `Wattles inside` = "Wa-In", `Wattles outside` = "Wa-Out",
        `Launceston inside` = "La-In", `Launceston outside` = "La-Out",
        `Needles inside` = "Ne-In", `Needles outside` = "Ne-Out",
        `Mole inside` = "Mol-In", `Mole outside` = "Mol-Out",
        `Camberra inside` = "Ca-In", `Camberra outside` = "Ca-Out",
        `Warri inside` = "Wi-In", `Warri outside` = "Wi-Out",
        `Braidwood inside` = "Br-In", `Braidwood outside` = "Br-Out",
        `Manjimup inside` = "Ma-In", `Manjimup outside` = "Ma-Out",
        `Pemberton inside` = "Pem-In", `Pemberton outside` = "Pem-Out",
        `Jardee inside` = "Ja-In", `Jardee outside` = "Ja-Out",
        `Cuneo inside` = "Cu-In", `Cuneo outside` = "Cu-Out",
        `San Demetrio inside` = "SD-In", `San Demetrio outside` = "SD-Out",
        `Spoleto inside` = "Sp-In", `Spoleto outside` = "Sp-Out",
        `Norcia inside` = "No-In", `Norcia outside` = "No-Out",
        `Cognac inside` = "Co-In", `Cognac outside` = "Co-Out",
        `Grignan inside` = "Gr-In", `Grignan outside` = "Gr-Out",
        `Perpignan inside` = "Per-In", `Perpignan outside` = "Per-Out",
        `Nimes inside` = "Ni-In", `Nimes outside` = "Ni-Out",
        `Romans-Sur-Isere inside` = "Ro-In", `Romans-Sur-Isere outside` = "Ro-Out",
        `Albentosa inside` = "Al-In", `Albentosa outside` = "Al-Out",
        `Moncayo inside` = "Mon-In", `Moncayo outside` = "Mon-Out",
        `Zuniga inside` = "Zu-In", `Zuniga outside` = "Zu-Out",
        `Acedo inside` = "Ac-In", `Acedo outside` = "Ac-Out")) %>% 
      mutate_at(c("hill_0", "hill_1", "hill_2"), as.numeric) %>% 
      filter(management %in% c("cultivated")) %>%
      droplevels()
  ) 

head(physeq_ecm_all_all@sam_data)

# Bacteria ---------------------------------------------------------------------
physeq_bact_all <- AlphaMetrics(physeq_bact_rare)

physeq_bact_all@sam_data %>% head()

sample_data(physeq_bact_all)$sample_id <- 
  sample_names(physeq_bact_all)

# rearrange factors
sample_data(physeq_bact_all) <-
  sample_data(
    as.data.frame(as.matrix(physeq_bact_all@sam_data)) %>% 
      mutate(site = recode(site,
                           S.Demetrio = "San Demetrio",
                           Uzes = "Nimes",
                           Manjimup2 = "Jardee", 
                           Chateauvert = "Grignan",
                           Angouleme = "Cognac",
                           `Espira-de-l-Agly` = "Perpignan",
                           Parnans = "Romans-Sur-Isere")) %>% 
      mutate(
        # extract the tree number right before N/F/B/HB at end of sample_id
        tree_num = str_extract(sample_id, "\\d+(?=(N|F|B|HB)$)"),
        tree_id   = paste0(site_code, tree_num)
      ) %>% 
      mutate(continent_brule = interaction(continent, brule, sep = "_")) %>% 
      mutate(site_brule = interaction(site, brule, sep = "_")) %>% 
      select(sample_id, tree_id, 
             continent, continent_brule, 
             site, site_brule, 
             brule, management,
             hill_0, hill_1, hill_2) %>% 
      mutate(hill_0 = as.numeric(hill_0), 
             hill_1 = as.numeric(hill_1),
             hill_2 = as.numeric(hill_2)) %>%
      mutate(tree_id = recode(tree_id, EP219 = "EP290", EP288 = "EP298"),
             tree_id = recode(tree_id, EH5 = "EH1")) %>% 
      janitor::clean_names() %>% 
      mutate(site = fct_relevel(site, 
                                "Yarra", "Wattles", "Launceston", "Needles", "Mole",
                                "Camberra", "Warri", "Braidwood", "Manjimup",
                                "Pemberton", "Jardee", 
                                "Cuneo", "San Demetrio", "Spoleto", "Norcia",
                                "Cognac", "Grignan", "Perpignan", "Nimes", "Romans-Sur-Isere",
                                "Albentosa", "Moncayo", "Zuniga", "Acedo")) %>% 
      mutate(site_brule = paste(site, brule, sep = " ")) %>% 
      mutate(site_brule = recode_factor(
        site_brule,
        `Yarra inside` = "Ya-In", `Yarra outside` = "Ya-Out",
        `Wattles inside` = "Wa-In", `Wattles outside` = "Wa-Out",
        `Launceston inside` = "La-In", `Launceston outside` = "La-Out",
        `Needles inside` = "Ne-In", `Needles outside` = "Ne-Out",
        `Mole inside` = "Mol-In", `Mole outside` = "Mol-Out",
        `Camberra inside` = "Ca-In", `Camberra outside` = "Ca-Out",
        `Warri inside` = "Wi-In", `Warri outside` = "Wi-Out",
        `Braidwood inside` = "Br-In", `Braidwood outside` = "Br-Out",
        `Manjimup inside` = "Ma-In", `Manjimup outside` = "Ma-Out",
        `Pemberton inside` = "Pem-In", `Pemberton outside` = "Pem-Out",
        `Jardee inside` = "Ja-In", `Jardee outside` = "Ja-Out",
        `Cuneo inside` = "Cu-In", `Cuneo outside` = "Cu-Out",
        `San Demetrio inside` = "SD-In", `San Demetrio outside` = "SD-Out",
        `Spoleto inside` = "Sp-In", `Spoleto outside` = "Sp-Out",
        `Norcia inside` = "No-In", `Norcia outside` = "No-Out",
        `Cognac inside` = "Co-In", `Cognac outside` = "Co-Out",
        `Grignan inside` = "Gr-In", `Grignan outside` = "Gr-Out",
        `Perpignan inside` = "Per-In", `Perpignan outside` = "Per-Out",
        `Nimes inside` = "Ni-In", `Nimes outside` = "Ni-Out",
        `Romans-Sur-Isere inside` = "Ro-In", `Romans-Sur-Isere outside` = "Ro-Out",
        `Albentosa inside` = "Al-In", `Albentosa outside` = "Al-Out",
        `Moncayo inside` = "Mon-In", `Moncayo outside` = "Mon-Out",
        `Zuniga inside` = "Zu-In", `Zuniga outside` = "Zu-Out",
        `Acedo inside` = "Ac-In", `Acedo outside` = "Ac-Out")) %>% 
      mutate_at(c("hill_0", "hill_1", "hill_2"), as.numeric) %>% 
      filter(management %in% c("cultivated")) %>%
      droplevels()
  ) 

head(physeq_bact_all@sam_data)

# ************************************************ -----------------------------
# BRAY and JACCARD DISTANCE - HEATMAPS -----------------------------------------

# Generate heatmaps for sites --------------------------------------------------
generateHeat(physeq_fungi_all)

# **** FIGURE 1 **** -----------------------------------------------------------
# Heatmaps global diversity
figure_a_betadif <-
  ggarrange(
    generateHeat(physeq_fungi_all) + labs(title = "Fungi"),
    generateHeat(physeq_ecm_all_all) + labs(title = "Ectomycorrhizal fungi"),
    generateHeat(physeq_bact_all) + labs(title = "Bacteria"),
    ncol = 3,
    labels = c("A", "B", "C"),
    common.legend = TRUE,
    legend = "right"
  )


grid.arrange(
  figure_a_betadif,
  top = text_grob("GLOBAL SCALE BETA DIVERSITY OF TRUFFLE ORCHARDS MICROBIOMES",
                  size = 12, face = 2)
)

# Compare bray and jaccard medians ---------------------------------------------

# objects with sites needed for the function to work!
australian_sites <- c("Yarra","Wattles","Launceston","Needles","Mole",
                      "Camberra","Warri","Braidwood","Manjimup",
                      "Pemberton","Jardee")

european_sites <-c("Cuneo","San Demetrio","Spoleto","Norcia",
                   "Cognac","Grignan","Perpignan","Nimes","Romans-Sur-Isere",
                   "Albentosa", "Moncayo", "Zuniga","Acedo")

palette_cont_heat <- c("#C7DAD7", "#5991EE", "#075AFF",
                       "#FFC69E", "#FF715A", "#FF0000")

# **** FIGURE S3 **** -----------------------------------------------------------

figure_S3_median_comp1 <-
  ggarrange(
  grid.arrange(
    compareBetadiv(physeq_fungi_all),
    top =title1),
  grid.arrange(
    compareBetadiv(physeq_ecm_all_all),
    top=title2),
  grid.arrange(
    compareBetadiv(physeq_bact_all),
    top=title3),
  nrow = 1,
  ncol = 3)

figure_S3_median_comp1

# Generate heatmaps for sites x brule ------------------------------------------
generateHeatSplit(physeq_fungi_all, "Australia")
generateHeatSplit(physeq_bact_all, "Europe") 

median(c( 
  (0.9750336 - 0.5610267)/2 +  0.5610267,
  (1.000000 - 0.234778)/2 + 0.234778,
  (0.7686971 - 0.3811134 )/2 + 0.3811134,
  
  (0.9750336 - 0.5610267)/2 +  0.5610267,
  (1.000000 - 0.2957093)/2 + 0.2957093,
  (0.8476821 - 0.4485012 )/2 + 0.4485012
))

# **** FIGURE 2 **** -----------------------------------------------------------

figure_a_betadif_EU <-
  ggarrange(
    physeq_fungi_all %>% 
      subset_samples(continent %in% c("Europe")) %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.) %>% 
      generateHeatSplit(., "Europe") + 
      labs(title = "Fungi", subtitle = "Europe"),
    physeq_ecm_all_all %>% 
      subset_samples(continent %in% c("Europe")) %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.) %>% 
      generateHeatSplit(.,"Europe") + 
      labs(title = "Ectomycorrhizal fungi", subtitle = "Europe"),
    physeq_bact_all %>% 
      subset_samples(continent %in% c("Europe")) %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.) %>% 
      generateHeatSplit(.,"Europe") + 
      labs(title = "Bacteria", subtitle = "Europe"),
    labels = c("A","B","C"),
    ncol = 3, 
    nrow = 1, 
    align = "hv",
    legend = "none")

figure_a_betadif_EU

figure_a_betadif_AU <-
  ggarrange(
    physeq_fungi_all %>% 
      subset_samples(continent %in% c("Australia")) %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.) %>% 
      generateHeatSplit(.,"Australia" ) + 
      labs(title = "Fungi", subtitle = "Australia"),
    physeq_ecm_all_all %>% 
      subset_samples(continent %in% c("Australia")) %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.) %>% 
      generateHeatSplit(.,"Australia") + 
      labs(title = "Ectomycorrhizal fungi", subtitle = "Australia"), 
    physeq_bact_all %>% 
      subset_samples(continent %in% c("Australia")) %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.) %>% 
      generateHeatSplit(.,"Australia") + 
      labs(title = "Bacteria", subtitle = "Australia"),
    labels = c("D","E","F"),
    nrow = 1,
    ncol = 3,
    align = "hv",
    common.legend = TRUE,
    legend = "bottom")

figure_a_betadif_AU


heat_split_multiplot <- 
  ggarrange(
    figure_a_betadif_EU,
    figure_a_betadif_AU + theme(plot.margin = unit(c(0,32,0,62), "pt")),
    nrow = 2)

heat_split_multiplot

# Generating a legend for site to plot -----------------------------------------
legend_data <-
  data.frame(Groups = c("Ya = Yarra",
                        "Wa = Wattles",
                        "La = Launceston",
                        "Ne = Needles",
                        "Mol = Mole",
                        "Ca = Camberra",
                        "Wi = Warri",
                        "Br = Braidwood",
                        "Ma = Manjimup", 
                        "Pem = Pemberton", 
                        "Ja = Jardee", 
                        "Cu = Cuneo",
                        "SD = San Demetrio",
                        "Sp = Spoleto",
                        "No = Norcia",
                        "Co = Cognac",
                        "Gr = Grignan",
                        "Per = Perpignan",
                        "Ni = Nimes",
                        "Ro = Romans-Sur-Isere", 
                        "Al = Albentosa", 
                        "Mon = Moncayo", 
                        "Zu = Zuniga",
                        "Ac = Acedo",
                        "",
                        "In = Inside brûlé",
                        "Out = Outside brûlé"))

ggplot_table_all <-
  ggtexttable(
    legend_data,
    rows = NULL,
    cols = NULL,
    theme = ttheme(
      padding = unit(c(1, 1), "mm"),
      tbody.style = tbody_style(
        fill = "white",
        color = "black", #grey30
        hjust = 0,
        vjust = 0.5,
        size = 9,
        linewidth = 0.05,
        x = 0)))

ggplot_table_all
ggplot_table_all + theme(plot.margin=ggplot2::margin(0,0,0,0))


Fig_heat_split_all <-
  ggarrange(
    heat_split_multiplot,
    ggplot_table_all + theme(plot.margin=ggplot2::margin(0,0,0,0)),
    ncol = 2,
    nrow = 1,
    widths = c(1, 0.13)
  )


grid.arrange(
  Fig_heat_split_all,
  top = text_grob("GLOBAL SCALE BETA DIVERSITY OF TRUFFLE ORCHARDS MICROBIOMES
                  INSIDE AND OUTSIDE THE BRÛLÉ", 
                  size = 12, face = 2))


# Compare bray and jaccard medians sites x brule -------------------------------

# **** FIGURE S4 **** -----------------------------------------------------------
figure_S4_median_comp_EU <-
ggarrange(
  grid.arrange(
    physeq_fungi_all %>% 
      subset_samples(continent %in% c("Europe")) %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.) %>% 
      compareBetadivSplit(.),
    top =title1),
  grid.arrange(
    physeq_ecm_all_all %>% 
      subset_samples(continent %in% c("Europe")) %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.) %>% 
      compareBetadivSplit(.),
    top=title2),
  grid.arrange(
    physeq_bact_all %>% 
      subset_samples(continent %in% c("Europe")) %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.) %>% 
      compareBetadivSplit(.),
    top=title3),
  nrow = 1,
  ncol = 3)

figure_S4_median_comp_EU

# **** FIGURE S5 **** -----------------------------------------------------------
figure_S4_median_comp_AU <-
ggarrange(
  grid.arrange(
    physeq_fungi_all %>% 
      subset_samples(continent %in% c("Australia")) %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.) %>% 
      compareBetadivSplit(.),
    top =title1),
  grid.arrange(
    physeq_ecm_all_all %>% 
      subset_samples(continent %in% c("Australia")) %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.) %>% 
      compareBetadivSplit(.),
    top=title2),
  grid.arrange(
    physeq_bact_all %>% 
      subset_samples(continent %in% c("Australia")) %>%
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.) %>% 
      compareBetadivSplit(.),
    top=title3),
  nrow = 1,
  ncol = 3)


figure_S4_median_comp_AU


# ************************************************ -----------------------------
# PRINCIPAL COORDINATE ANALYSIS ------------------------------------------------
# Fungi ------------------------------------------------------------------------

palette_site <-c("#cc1c1c","#560d0d","#a35151", "#dba4a4",
                 "#ea7f17","#fcb067","#FDDB8E", "#ffff9e",
                 "#111b77","#283dff","#058ED9","#bfc5ff",
                 "#ffb7ef","#fa7efc","#ae09ea","#521899",
                 "#014443","#117744","#60ffaf","#b7ffdb",
                 "#d8d6d4","#82807f","#825121","#000000")


# whole
 plot_pcoa_analysis(physeq_fungi_all,
                    title = "Fungi", 
                    method = "jaccard",
                    show_legend = "right")
 
 plot_pcoa_analysis(physeq_fungi_all,
                    title = "Fungi", 
                    method = "bray",
                    show_legend = "right")
 
# split
 plot_pcoa_analysis(physeq_fungi_all,
                    title = "Fungi", 
                    subset_continent = "Australia",
                    method = "jaccard",
                    show_legend = "right")
 
# Plotting ---------------------------------------------------------------------

figure_S1_pcoa_main <-
  ggarrange(
    plot_pcoa_analysis(physeq_fungi_all,
                       title = "Fungi", 
                       method = "jaccard",
                       show_legend = "none"), 

    plot_pcoa_analysis(physeq_ecm_all,
                       title="Ectomycorrhizal fungi",
                       method = "jaccard",
                       show_legend = "none"),

    plot_pcoa_analysis(physeq_bact_all,
                       title = "Bacteria", 
                       method = "jaccard",
                       show_legend = "none"),
    ncol = 3, 
    nrow = 1, 
    align = "hv",
    labels = c("A", "B", "C"))

 figure_S1_pcoa_main 

figure_S1_pcoa_split <-
  ggarrange(
    subset_samples(physeq_fungi_all, continent%in%c("Europe")) %>% 
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.) %>% 
      plot_pcoa_analysis(title = "Fungi", 
                         subset_continent = "Europe",
                         method = "jaccard",
                         show_legend = "none"),
    subset_samples(physeq_fungi_all, continent%in%c("Australia")) %>% 
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.) %>% 
      plot_pcoa_analysis(title = "Fungi", 
                         subset_continent = "Australia",
                         method = "jaccard",
                         show_legend = "none"),
    
    subset_samples(physeq_ecm_all, continent%in%c("Europe")) %>% 
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.) %>% 
      plot_pcoa_analysis(title = "Ectomycorrhizal fungi", 
                         subset_continent = "Europe",
                         method = "jaccard",
                         show_legend = "none"),
    subset_samples(physeq_ecm_all, continent%in%c("Australia")) %>% 
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.) %>% 
      plot_pcoa_analysis(title = "Ectomycorrhizal fungi", 
                         subset_continent = "Australia",
                         method = "jaccard",
                         show_legend = "none"),
    
    subset_samples(physeq_bact_all, continent%in%c("Europe")) %>% 
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.) %>% 
      plot_pcoa_analysis(title = "Bacteria", 
                         subset_continent = "Europe",
                         method = "jaccard",
                         show_legend = "none"),
    subset_samples(physeq_bact_all, continent%in%c("Australia")) %>% 
      prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
      prune_samples(sample_sums(x=.) > 0, x =.) %>% 
      plot_pcoa_analysis(title = "Bacteria", 
                         subset_continent = "Australia",
                         method = "jaccard",
                         show_legend = "none"),
    
    ncol = 6, 
    nrow = 1,
    align = "hv",
    legend = "none")

figure_S1_pcoa_split

# Create legend ----------------------------------------------------------------
c(rep("#00843D",11), rep("#003399", 13))

pcoa_legend <- as_ggplot(
  get_legend(
    plot_pcoa_analysis(physeq_fungi_all,
                       title = "Fungi", 
                       method = "jaccard",
                       show_legend = "right") +
      themePlot() +
      theme(legend.title =element_blank(), 
            legend.text = element_markdown(size = 8),
            legend.spacing.y = unit(-0.1, "cm")) +
      guides(color = guide_legend(ncol = 1, 
                                  override.aes = list(shape = 15, size = 3.5)),
             shape = guide_legend(ncol=1, 
                                  override.aes = list(color = "black", size=2.5))) +
      scale_color_manual(values = palette_site) +
      scale_shape_manual(values = c(16,17)))
)

pcoa_legend

# **** FIGURE S6 ---------------------------------------------------------------
# need to fix axis tick labels, not sure why they do not align correctly.
Fig2_all <-
  ggarrange(
    ggarrange(
      figure_S1_pcoa_main,
      figure_S1_pcoa_split,
      nrow = 2,
      ncol = 1,
      heights = c(0.65,0.35)
    ),
    pcoa_legend,
    ncol = 2,
    nrow = 1,
    widths = c(0.7, 0.1), 
    align = "hv"
  )

Fig2_all

# ************************************************ -----------------------------
# PERMANOVA --------------------------------------------------------------------

# permanova bray curtis --------------------------------------------------------
physeq_fungi_all %>%
  subset_samples(management %in% c("cultivated")) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  Adonis2All(method = "bray")

physeq_ecm_all %>%
  subset_samples(management %in% c("cultivated")) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  Adonis2All(method = "bray")

physeq_bact_all %>%
  subset_samples(management %in% c("cultivated")) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  Adonis2All(method = "bray")

# permanova jaccard ------------------------------------------------------------
physeq_fungi_all %>%
  subset_samples(management %in% c("jaccard")) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  Adonis2All(method = "bray")

physeq_ecm_all %>%
  subset_samples(management %in% c("jaccard")) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  Adonis2All(method = "bray")

physeq_bact_all %>%
  subset_samples(management %in% c("jaccard")) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  Adonis2All(method = "bray")

# BETA-DISPERSION ---------------------------------------------------------------

# betadisper bray curtis --------------------------------------------------------
physeq_fungi_all %>%
  subset_samples(management %in% c("cultivated")) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  multipleBetadisper(method = "bray")

physeq_ecm_all %>%
  subset_samples(management %in% c("cultivated")) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  multipleBetadisper(method = "bray")

physeq_bact_all %>%
  subset_samples(management %in% c("cultivated")) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  multipleBetadisper(method = "bray")

# betadisper jaccard ------------------------------------------------------------
physeq_fungi_all %>%
  subset_samples(management %in% c("jaccard")) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  multipleBetadisper(method = "bray")

physeq_ecm_all %>%
  subset_samples(management %in% c("jaccard")) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  multipleBetadisper(method = "bray")

physeq_bact_all %>%
  subset_samples(management %in% c("jaccard")) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  Adonis2All(method = "bray")

# Split permanova on Australia and Europe separately ---------------------------

# Split permanova bray curtis --------------------------------------------------

# Fungi

rbund(
  rbind(
    physeq_fungi_all %>%
      subset_samples(management %in% c("cultivated")) %>%
      prune_taxa(taxa_sums(.) > 0, .) %>%
      prune_samples(sample_sums(.) > 0, .) %>%
      Adonis2Split(method = "bray", subset_continent = c("Europe")) %>%
      mutate(Continent = "Europe"),
    physeq_fungi_all %>%
      subset_samples(management %in% c("cultivated")) %>%
      prune_taxa(taxa_sums(.) > 0, .) %>%
      prune_samples(sample_sums(.) > 0, .) %>%
      Adonis2Split(method = "bray", subset_continent = c("Australia")) %>%
      mutate(Continent = "Australia")
  ) %>% mutate(Kingdom = "Fungi"),

  # Ectomycorrhial fungi
  rbind(
    physeq_ecm_all %>%
      subset_samples(management %in% c("cultivated")) %>%
      prune_taxa(taxa_sums(.) > 0, .) %>%
      prune_samples(sample_sums(.) > 0, .) %>%
      Adonis2Split(method = "bray", subset_continent = c("Europe")) %>%
      mutate(Continent = "Europe"),
    physeq_ecm_all %>%
      subset_samples(management %in% c("cultivated")) %>%
      prune_taxa(taxa_sums(.) > 0, .) %>%
      prune_samples(sample_sums(.) > 0, .) %>%
      Adonis2Split(method = "bray", subset_continent = c("Australia")) %>%
      mutate(Continent = "Australia")
  ) %>% mutate(Kingdom = "Ectomycorrhizal fungi"),

  # Bacteria
  rbind(
    physeq_bact_all %>%
      subset_samples(management %in% c("cultivated")) %>%
      prune_taxa(taxa_sums(.) > 0, .) %>%
      prune_samples(sample_sums(.) > 0, .) %>%
      Adonis2Split(method = "bray", subset_continent = c("Europe")) %>%
      mutate(Continent = "Europe"),
    physeq_bact_all %>%
      subset_samples(management %in% c("cultivated")) %>%
      prune_taxa(taxa_sums(.) > 0, .) %>%
      prune_samples(sample_sums(.) > 0, .) %>%
      Adonis2Split(method = "bray", subset_continent = c("Australia")) %>%
      mutate(Continent = "Australia")
  ) %>% mutate(Kingdom = "Bacteria")
)

# Split permanova jaccard-------------------------------------------------------

rbund(
  rbind(
    physeq_fungi_all %>%
      subset_samples(management %in% c("cultivated")) %>%
      prune_taxa(taxa_sums(.) > 0, .) %>%
      prune_samples(sample_sums(.) > 0, .) %>% 
      Adonis2Split(method="jaccard", subset_continent = c("Europe")) %>% 
      mutate(Continent = "Europe"),
    
    physeq_fungi_all %>%
      subset_samples(management %in% c("cultivated")) %>%
      prune_taxa(taxa_sums(.) > 0, .) %>%
      prune_samples(sample_sums(.) > 0, .) %>% 
      Adonis2Split(method="jaccard", subset_continent = c("Australia")) %>% 
      mutate(Continent = "Australia")
  ) %>% mutate(Kingdom = "Fungi"),
  
  # Ectomycorrhial fungi
  rbind(
    physeq_ecm_all %>%
      subset_samples(management %in% c("cultivated")) %>%
      prune_taxa(taxa_sums(.) > 0, .) %>%
      prune_samples(sample_sums(.) > 0, .) %>% 
      Adonis2Split(method="jaccard", subset_continent = c("Europe")) %>% 
      mutate(Continent = "Europe"),
    
    physeq_ecm_all %>%
      subset_samples(management %in% c("cultivated")) %>%
      prune_taxa(taxa_sums(.) > 0, .) %>%
      prune_samples(sample_sums(.) > 0, .) %>% 
      Adonis2Split(method="jaccard", subset_continent = c("Australia")) %>% 
      mutate(Continent = "Australia")
  )%>% mutate(Kingdom = "Ectomycorrhizal fungi"),
  
  # Bacteria 
  rbind(
    physeq_bact_all %>%
      subset_samples(management %in% c("cultivated")) %>%
      prune_taxa(taxa_sums(.) > 0, .) %>%
      prune_samples(sample_sums(.) > 0, .) %>% 
      Adonis2Split(method="jaccard", subset_continent = c("Europe")) %>% 
      mutate(Continent = "Europe"),
    
    physeq_bact_all %>%
      subset_samples(management %in% c("cultivated")) %>%
      prune_taxa(taxa_sums(.) > 0, .) %>%
      prune_samples(sample_sums(.) > 0, .) %>% 
      Adonis2Split(method="jaccard", subset_continent = c("Australia")) %>% 
      mutate(Continent = "Australia")
  ) %>% mutate(Kingdom = "Bacteria")
)

# SPLIT BETADISPER -------------------------------------------------------------

# # Split betadisper on Australia datastes
betadisp_fungi_AU_bray <-
  spliteBetadisper(physeq_fungi_all, "bray", "Australia")
betadisp_fungi_AU_bray
write.csv(betadisp_fungi_AU_bray, "betadisp_fungi_AU_bray.csv")

betadisp_bact_AU_bray <-
  spliteBetadisper(physeq_bact_all, "bray", "Australia")
write.csv(betadisp_bact_AU_bray, "betadisp_bact_AU_bray.csv")

betadisp_ecm_AU_bray <-
  spliteBetadisper(physeq_ecm, "bray", "Australia")
write.csv(betadisp_ecm_AU_bray, "betadisp_ecm_AU_bray.csv")

betadisp_fungi_AU_jac <-
  spliteBetadisper(physeq_fungi_all, "jac", "Australia")
betadisp_fungi_AU_jac
write.csv(betadisp_fungi_AU_jac, "betadisp_fungi_AU_jac.csv")

betadisp_bact_AU_jac <-
  spliteBetadisper(physeq_bact_all, "jac", "Australia")
write.csv(betadisp_bact_AU_jac, "betadisp_bact_AU_jac.csv")

betadisp_ecm_AU_jac <-
  spliteBetadisper(physeq_ecm, "jac", "Australia")
write.csv(betadisp_ecm_AU_jac, "betadisp_ecm_AU_jac.csv")


# Split betadisper on Europe datastes
betadisp_fungi_EU_bray <-
  spliteBetadisper(physeq_fungi_cult, "bray", "Europe")
betadisp_fungi_EU_bray
write.csv(betadisp_fungi_EU_bray, "betadisp_fungi_EU_bray.csv")

betadisp_bact_EU_bray <-
  spliteBetadisper(physeq_bact_cult, "bray", "Europe")
write.csv(betadisp_bact_EU_bray, "betadisp_bact_EU_bray.csv")

betadisp_ecm_EU_bray <-
  spliteBetadisper(physeq_ecm_cult, "bray", "Europe")
write.csv(betadisp_ecm_EU_bray, "betadisp_ecm_EU_bray.csv")

betadisp_fungi_EU_jac <-
  spliteBetadisper(physeq_fungi_cult, "jac", "Europe")
betadisp_fungi_EU_jac
write.csv(betadisp_fungi_EU_jac, "betadisp_fungi_EU_jac.csv")

betadisp_bact_EU_jac <-
  spliteBetadisper(physeq_bact_cult, "jac", "Europe")
write.csv(betadisp_bact_EU_jac, "betadisp_bact_EU_jac.csv")

betadisp_ecm_EU_jac <-
  spliteBetadisper(physeq_ecm_cult, "jac", "Europe")
write.csv(betadisp_ecm_EU_jac, "betadisp_ecm_EU_jac.csv")

# **** Table S1 - adonis and betadisper -----------------------------------------

# PERMANOVA summary -----------------------------------------------------------

# Permanova and betadisper run and analysis generate a data frame for Figure 3 

permanova_summary <-
  rbind(
    read.csv(file.path(project_dir, "results/csv/permanova_all.csv")) %>% 
      mutate(Continent = "Entire dataset") %>% 
      dplyr::select(Dataset,Method,Continent,Factor,R2, betadisper),
    read.csv(file.path(project_dir, "results/csv/permanova_split.csv"))
    ) %>% 
      mutate(stars = ifelse(betadisper<=0.05, paste("*"), NA)) %>% 
      mutate(lls = ifelse(!is.na(stars), paste(R2, stars, sep = " "), R2))


permanova_summary 
permanova_summary 

saveRDS(file = file.path(project_dir, "github//permanova_summary.rds"),
        permanova_summary)

# **** FIGURE 2 permanova summary stats ----------------------------------------

# plotting the R2 in easy interpretable way

figure3_permanova_summary <-
  permanova_summary %>% 
  mutate(Dataset = fct_relevel(Dataset, "Fungi", "Ectomycorrhizal fungi", "Bacteria"),
         Factor =fct_relevel(Factor, "continent","brule","site", "continent:brule","site:brule"),
         Continent =fct_relevel(Continent, "Entire dataset","Europe","Australia")) %>% 
  ggplot(aes(x=Dataset, y=R2, fill=Factor, color=Factor)) +
  geom_bar(stat = "identity") +
  #geom_text_repel(aes(label = round(R2, 3)), colour = "white") +
  geom_text(aes(label = lls), 
            size = 3, 
            position=position_stack(vjust=0.6), 
            colour = "black") +
  facet_wrap(~Method~Continent, nrow=1) +
  themePlot() +
  theme(
    plot.title= element_markdown(size = 12),
    axis.title.x = element_blank(),
    axis.text.x = element_markdown(size=10, angle=25, hjust = 1, vjust = 1),
    axis.text.y = element_markdown(size =10),
    strip.text.x = element_markdown(size= 10),
    legend.text = element_markdown(size=10, face = "italic"),
    legend.position = "bottom",
    legend.title = element_blank()) +
  labs(title = "Distance-based Permutational Multivariate Analysis of Variance", 
       y=bquote(R^2)) +
  scale_fill_manual(values = c("#CC2D35","#599861","#FF934F","#848FA2","#825121")) +
  scale_color_manual(values = c("#CC2D35","#599861","#FF934F","#848FA2","#825121"))


figure3_permanova_summary


# ************************************************ -----------------------------
# PAIRWISE PERMANOVA -----------------------------------------------------------
# modified from https://gist.github.com/mcgoodman/58c9d1257fd1625954a4ffa1c3301939

Pairwise_adonis2(physeq_fungi_all, "site")

# Saving the tables
pairwise_permanova(physeq_fungi_all, "site") %>% 
  write.csv(., "pairwPerm_fungi_site.csv")
pairwise_permanova(physeq_ecm_all, "Site") %>% 
  write.csv(., "pairwPerm_ecm_site.csv")
pairwise_permanova(physeq_bact_all, "Site") %>% 
  write.csv(., "pairwPerm_bact_site.csv")

pairwise_permanova(physeq_fungi_all, "Site", dist = "jaccard") %>% 
  write.csv(., "pairwPerm_fungi_site_jac.csv")
pairwise_permanova(physeq_ecm_all, "Site", dist = "jaccard") %>% 
  write.csv(., "pairwPerm_ecm_site_jac.csv")
pairwise_permanova(physeq_bact_all, "Site", dist = "jaccard") %>% 
  write.csv(., "pairwPerm_bact_site_jac.csv")


# **** FIGURE S7 **** ---------------------------------------------------------- 

# Adonis Pairwise comparisons bray ---------------------------------------------
heatmap_plot_bray  <- 
  ggarrange(
    getPairCompare(physeq_fungi_all, distance = "bray") %>% 
      PairPermCorrplot(.) +
      scale_fill_gradient(expression(R^2), low="red", high="yellow", na.value = "white") +
      labs(title = "Fungi"),
    getPairCompare(physeq_ecm_all, distance = "bray") %>% 
      PairPermCorrplot(.) +
      scale_fill_gradient(expression(R^2), low="red", high="yellow", na.value = "white") +
      labs(title = "Ectomycorrhizal fungi"),
    getPairCompare(physeq_bact_all, distance = "bray") %>% 
      PairPermCorrplot(.) +
      scale_fill_gradient(expression(R^2), low="red", high="yellow", na.value = "white") +
      labs(title = "Bacteria"),
    ncol = 3,
    nrow = 1,
    labels = c("A", "B", "C")
  )

heatmap_plot_bray

# donis Pairwise comparisons jaccard -------------------------------------------
heatmap_plot_jac  <- 
  ggarrange(
    getPairCompare(physeq_fungi_all, distance = "jaccard") %>% 
      PairPermCorrplot(.) +
      scale_fill_gradient(expression(R^2), low="red", high="yellow", na.value = "white") +
      labs(title = "Fungi"),
    getPairCompare(physeq_ecm_all, distance = "jaccard") %>% 
      PairPermCorrplot(.) +
      scale_fill_gradient(expression(R^2), low="red", high="yellow", na.value = "white") +
      labs(title = "Ectomycorrhizal fungi"),
    getPairCompare(physeq_bact_all, distance = "jaccard") %>% 
      PairPermCorrplot(.) +
      scale_fill_gradient(expression(R^2), low="red", high="yellow", na.value = "white") +
      labs(title = "Bacteria"),
    ncol = 3,
    nrow = 1,
    labels = c("D", "E", "F")
  )

heatmap_plot_jac

# ************************************************ -----------------------------
# MVABUND - modeling microbiome using glm --------------------------------------

mvabund_fungi_cult <- DFmvabund(physeq_fungi_cult, NULL)
mvabund_ecm_cult <- DFmvabund(physeq_ecm_cult, NULL)
mvabund_bact_cult <- DFmvabund(physeq_bact_cult, NULL)
mvabund_fungi_AU <- DFmvabund(physeq_fungi_cult, "Australia")
mvabund_fungi_EU <- DFmvabund(physeq_fungi_cult, "Europe")
mvabund_ecm_AU <- DFmvabund(physeq_ecm_cult, "Australia")
mvabund_ecm_EU <- DFmvabund(physeq_ecm_cult, "Europe")
mvabund_bact_AU <- DFmvabund(physeq_bact_cult, "Australia")
mvabund_bact_EU <- DFmvabund(physeq_bact_cult, "Europe")

mvabund_bact_EU[[2]]$Site
mvabund_fungi_cult[[1]] 

saveRDS(mvabund_fungi_cult, file = "mvabund_dsets/mvabund_fungi_cult.rds")
saveRDS(mvabund_ecm_cult, file = "mvabund_dsets/mvabund_ecm_cult.rds")
saveRDS(mvabund_bact_cult, file = "mvabund_dsets/mvabund_bact_cult.rds")

saveRDS(mvabund_fungi_AU, file = "mvabund_dsets/mvabund_fungi_AU.rds")
saveRDS(mvabund_fungi_EU, file = "mvabund_dsets/mvabund_fungi_EU.rds")
saveRDS(mvabund_ecm_AU, file = "mvabund_dsets/mvabund_ecm_AU.rds")
saveRDS(mvabund_ecm_EU, file = "mvabund_dsets/mvabund_ecm_EU.rds")
saveRDS(mvabund_bact_AU, file = "mvabund_dsets/mvabund_bact_AU.rds")
saveRDS(mvabund_bact_EU, file = "mvabund_dsets/mvabund_bact_EU.rds")

# Run this on the HPCC since it takes forever....
# re-importing results

# **** Table 2 - many glm ------------------------------------------------------
readRDS("mvabund_dsets/mvabund_fungi_AU_result.rds")$table
readRDS("mvabund_dsets/mvabund_fungi_EU_result.rds")$table
readRDS("mvabund_dsets/mvabund_ecm_AU_result.rds")$table
readRDS("mvabund_dsets/mvabund_ecm_EU_result.rds")$table
readRDS("mvabund_dsets/mvabund_bact_AU_result.rds")$table
readRDS("mvabund_dsets/mvabund_bact_EU_result.rds")$table


write.csv(readRDS("mvabund_dsets/mvabund_fungi_EU_result.rds")$table, "mvabund_fungi_EU_result.csv")
write.csv(readRDS("mvabund_dsets/mvabund_fungi_AU_result.rds")$table, "mvabund_fungi_AU_result.csv")
write.csv(readRDS("mvabund_dsets/mvabund_ecm_EU_result.rds")$table, "mvabund_ecm_EU_result.csv")
write.csv(readRDS("mvabund_dsets/mvabund_ecm_AU_result.rds")$table, "mvabund_ecm_AU_result.csv")
write.csv(readRDS("mvabund_dsets/mvabund_bact_EU_result.rds")$table, "mvabund_bact_EU_result.csv")
write.csv(readRDS("mvabund_dsets/mvabund_bact_AU_result.rds")$table, "mvabund_bact_AU_result.csv")


# ******************************************************************------------
# ALPHA DIVERSITY --------------------------------------------------------------

# NOTE! Older code was replaced by that in the revisions.R 

# Continent x Brule x Site -----------------------------------------------------

# Using dataframne generated for revisions!

df_alpha_fungi_all_new
df_alpha_ecm_new
df_alpha_bact_all_new

title4 = text_grob("Richness (hill_0)", size = 12, face = 2)
title5 = text_grob("invSimpson (hill_2)", size = 12, face = 2)

# Australia --------------------------------------------------------------------

PlotSiteRich_AU(df_alpha_fungi_all_new,
                formula(hill_0 ~ site_brule),
                "Richness",
                350,
                355,
                "Fungi")

# Australia Rarefied Richness -------------------------------------------------

alpha_plot_AU_sites <-
  grid.arrange(
    ggarrange(
      PlotSiteRich_AU(
        df_alpha_fungi_all_new,
        formula(hill_0 ~ site_brule),
        "Richness (hill_0)",
        350,
        350,
        "Fungi"
      ),
      PlotSiteRich_AU(
        df_alpha_ecm_new,
        formula(hill_0 ~ site_brule),
        "Richness (hill_0)",
        20,
        20,
        "Ectomycorrhizal fungi"
      ),
      PlotSiteRich_AU(
        df_alpha_bact_all_new,
        formula(hill_0 ~ site_brule),
        "Richness (hill_0)",
        1600,
        1600,
        "Bacteria"
      ),
      ncol = 1,
      nrow = 3,
      labels = c("A", "B", "C")
    ),
    top = title4
  )

alpha_plot_AU_sites


# Australia InverseSimpson ---------------------------------------------------------------

invSimp_plot_AU_sites <-
  grid.arrange(
    ggarrange(
      PlotSiteRich_AU(
        df_alpha_fungi_all_new,
        formula(hill_2 ~ site_brule),
        "invSimpson (hill_2)",
        50,
        50,
        "Fungi"
      ),
      PlotSiteRich_AU(
        df_alpha_ecm_new,
        formula(hill_2 ~ site_brule),
        "invSimpson (hill_2)",
        5,
        5,
        "Ectomycorrhizal fungi"
      ),
      PlotSiteRich_AU(
        df_alpha_bact_all_new,
        formula(hill_2 ~ site_brule),
        "invSimpson (hill_2)",
        375,
        375,
        "Bacteria"
      ),
      ncol = 1,
      nrow = 3,
      labels = c("A", "B", "C")
    ),
    top = title5
  )

invSimp_plot_AU_sites

# Adding abbreviation legend
legend_data_Au <-
  data.frame(Groups = c("Ya = Yarra",
                        "Wa = Wattles",
                        "La = Launceston",
                        "Ne = Needles",
                        "Mol = Mole",
                        "Ca = Camberra",
                        "Wi = Warri",
                        "Br = Braidwood",
                        "Ma = Manjimup", 
                        "Pem = Pemberton", 
                        "Ja = Jardee", 
                        "",
                        "In = Inside brûlé",
                        "Out = Outside brûlé"))

legend_data_Au

ggplot_table_Au <-
  ggtexttable(
    legend_data_Au,
    rows = NULL,
    cols = NULL,
    theme = ttheme(
      padding = unit(c(1, 1.8), "mm"),
      tbody.style = tbody_style(
        fill = "white",
        color = "black",
        hjust = 0,
        vjust = 0,
        size = 7,
        linewidth = 0.5,
        x = 0.02)))

ggplot_table_Au


# **** FIGURE S7 - richness x site x brule -----------------------------------
# plotting the combined graph and legend table
fig_rich_brule_AU <-
  ggarrange(
    alpha_plot_AU_sites,
    ggplot_table_Au,
    ncol = 2,
    nrow = 1,
    widths = c(1, 0.15)
  )

fig_rich_brule_AU

# **** FIGURE S8 - invSimpson x site x brule -----------------------------------
fig_invSimp_brule_AU <-
  ggarrange(
    invSimp_plot_AU_sites,
    ggplot_table_Au,
    ncol = 2,
    nrow = 1,
    widths = c(1, 0.15)
  )

fig_invSimp_brule_AU

# Europe -----------------------------------------------------------------------


alpha_plot_EU_sites <-
  grid.arrange(
    ggarrange(
      PlotSiteRich_EU(
        df_alpha_fungi_all_new,
        formula(hill_0 ~ site_brule),
        "Richness (hill_0)",
        500,
        530,
        "Fungi"
      ),
      PlotSiteRich_EU(
        df_alpha_ecm_new,
        formula(hill_0 ~ site_brule),
        "Richness (hill_0)",
        40,
        42,
        "Ectomycorrhizal fungi"
      ),
      PlotSiteRich_EU(
        df_alpha_bact_all_new,
        formula(hill_0 ~ site_brule),
        "Richness (hill_0)",
        1600,
        1670,
        "Bacteria"
      ),
      ncol = 1,
      nrow = 3,
      labels = c("A", "B", "C")
    ),
    top = title4
  )

alpha_plot_EU_sites


invSimp_plot_EU_sites <-
  grid.arrange(
    ggarrange(
      PlotSiteRich_EU(
        df_alpha_fungi_all_new,
        formula(hill_2 ~ site_brule),
        "invSimpson (hill_2)",
        75,
        80,
        "Fungi"
      ),
      PlotSiteRich_EU(
        df_alpha_ecm_new,
        formula(hill_2 ~ site_brule),
        "invSimpson (hill_2)",
        8,
        8.5,
        "Ectomycorrhizal fungi"
      ),
      PlotSiteRich_EU(
        df_alpha_bact_all_new,
        formula(hill_2 ~ site_brule),
        "invSimpson (hill_2)",
        400,
        450,
        "Bacteria"
      ),
      ncol = 1,
      nrow = 3,
      labels = c("A", "B", "C")
    ),
    top = title5
  )

invSimp_plot_EU_sites


# Adding abbreviation legend
legend_data_Eu <-
  data.frame(Groups = c("Cu = Cuneo",
                        "SD = San Demetrio",
                        "Sp = Spoleto",
                        "No = Norcia",
                        "Co = Cognac",
                        "Gr = Grignan",
                        "Per = Perpignan",
                        "Ni = Nimes",
                        "Ro = Romans-Sur-Isere", 
                        "Al = Albentosa", 
                        "Mon = Moncayo", 
                        "Zu = Zuniga",
                        "Ac = Acedo",
                        "",
                        "In = Inside brûlé",
                        "Out = Outside brûlé"))


legend_data_Eu

ggplot_table <-
  ggtexttable(
    legend_data_Eu,
    rows = NULL,
    cols = NULL,
    theme = ttheme(
      padding = unit(c(1, 1.8), "mm"),
      tbody.style = tbody_style(
        fill = "white",
        color = "black",
        hjust = 0,
        vjust = 0,
        size = 7,
        linewidth = 0.5,
        x = 0.02)))

ggplot_table

# **** FIGURE S5 - richness x site x brule -----------------------------------
# plotting the combined graph and legend table
fig_rich_brule_fungi <-
  ggarrange(
    alpha_plot_EU_sites,
    ggplot_table,
    ncol = 2,
    nrow = 1,
    widths = c(1, 0.2)
  )

fig_rich_brule_fungi

# **** FIGURE S6 - invSimpson x site x brule -----------------------------------
fig_invSimp_brule_fungi <-
  ggarrange(
    invSimp_plot_EU_sites,
    ggplot_table,
    ncol = 2,
    nrow = 1,
    widths = c(1, 0.2)
  )

fig_invSimp_brule_fungi

# Discretizing T. melanosporum abudnance ---------------------------------------

# **** FIGURE S12 **** ----------------------------------------------------------

# NOTE! See the revision.R for this. 

# ************************************************ -----------------------------

# **** FIGURE S13 **** ---------------------------------------------------------

# All Fungi, ECM and Bacteira --------------------------------------------------

# Fix taxonomies
head(tax_table(physeq_fungi_all))
head(tax_table(physeq_ecm_all))
head(tax_table(physeq_bact_all))

# Extrat dataframes for plotting top abundanc OTUs ----------------------------- 

ExtractTop(physeq_ecm_all, "Genus", 15)
ExtractTop(physeq_fungi_all, "Genus", 25) 

physeq_fungi_all %>% 
  subset_taxa(Genus=="Tuber") %>% 
  ExtractTop("Species", 20) 

ExtractTop(physeq_bact_all, "Class", 20) 

#How many genera?
physeq_fungi_all@tax_table %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  dplyr::select(Genus) %>% 
  summarise(unique_level_count = n_distinct(Genus))

all_genera <- 
  subset_taxa(physeq_fungi_all, Genus %in% c(
    ExtractTop(physeq_fungi_all, "Genus", 1000))) %>% 
  tax_glom(., taxrank="Genus") %>% 
  psmelt() %>% 
  ExtractAbund(Genus)

all_genera
write.csv(all_genera, "all_fungi_genera.csv")

subset_taxa(physeq_fungi_all, Genus %in% c(
  ExtractTop(physeq_fungi_all, "Genus", 31))) %>% 
  tax_glom(., taxrank="Genus") %>% 
  psmelt() %>% 
  ExtractAbund(Genus)

sum(sample_sums(physeq_fungi_all))

# just ectomycorrhizal
all_ecto <-
  subset_taxa(physeq_ecm, Genus %in% c(
    ExtractTop(physeq_ecm, "Genus", 1000))) %>% 
  tax_glom(., taxrank="Genus") %>% 
  psmelt() %>% 
  ExtractAbund(.,physeq_ecm, Genus)

all_ecto
write.csv(all_ecto, "all_ecto.csv")

all_species <-
  subset_taxa(physeq_ecm, Species %in% c(
    ExtractTop(physeq_ecm, "Species", 1000))) %>% 
  tax_glom(., taxrank="Species") %>% 
  psmelt() %>% 
  ExtractAbund(Species) %>% 
  arrange(Species)

all_species
write.csv(all_species, "all_species.csv")

# bacteria
first100_bact<- 
  subset_taxa(physeq_bact_all, Genus %in% c(
    ExtractTop(physeq_bact_all, "Genus", 100))) %>% 
  tax_glom(., taxrank="Genus") %>% 
  psmelt() %>% 
  ExtractAbund(.,physeq_bact_all, Genus)

first100_bact

subset_taxa(physeq_bact_all, Class %in% c(
  ExtractTop(physeq_bact_all, "Class", 25))) %>% 
  tax_glom(., taxrank="Class") %>% 
  psmelt() %>% 
  ExtractAbund(Class)


# Plot abundance barplots ------------------------------------------------------
colors_bar_ecm <- c("#D21E2C","#014443","#b7ffdb","#F7F7C5","#ea7f17",
                    "#bfc5ff","#111b77","#ae09ea","#FDDB8E","#283dff",
                    "#a35151","#117744","#dba4a4","#058ED9","#82807f",
                    "#ffb7ef","#ffff9e","#d8d6d4","#000000","#560d0d",
                    "#521899","#fcb067","#60ffaf","#825121","#fa7efc")

# ECM 
Fig_ECM_Abund <-
  subset_taxa(physeq_ecm_all, Genus %in% c(
    ExtractTop(physeq_ecm_all, "Genus", 20) 
  )) %>% 
  tax_glom(., taxrank="Genus") %>% 
  psmelt() %>% 
  arrange(Genus) %>% 
  filter(Genus != "Unclassified" & Genus != "unclassified") %>% 
  mutate(site_brule = recode_factor(
    site_brule,
    `Yarra inside` = "Ya-In" ,`Yarra outside` = "Ya-Out",
    `Wattles inside`="Wa-In",`Wattles outside`="Wa-Out",
    `Launceston inside` = "La-In" ,`Launceston outside` = "La-Out",
    `Needles inside`="Ne-In",`Needles outside`="Ne-Out",
    `Mole inside`="Mol-In",`Mole outside`="Mol-Out",
    `Camberra inside` = "Ca-In" ,`Camberra outside` = "Ca-Out",
    `Warri inside`="Wi-In",`Warri outside`="Wi-Out",
    `Braidwood inside`="Br-In",`Braidwood outside`="Br-Out",
    `Manjimup inside` = "Ma-In" ,`Manjimup outside` = "Ma-Out",
    `Pemberton inside`="Pem-In",`Pemberton outside`="Pem-Out",
    `Jardee inside`="Ja-In",`Jardee outside`="Ja-Out",
    `Cuneo inside` = "Cu-In" ,`Cuneo outside` = "Cu-Out",
    `San Demetrio inside`="SD-In",`San Demetrio outside`="SD-Out",
    `Spoleto inside`="Sp-In",`Spoleto outside`="Sp-Out",
    `Norcia inside`= "No-In",`Norcia outside`="No-Out",
    `Cognac inside` = "Co-In" ,`Cognac outside` = "Co-Out",
    `Grignan inside`="Gr-In",`Grignan outside`="Gr-Out",
    `Perpignan inside`="Per-In",`Perpignan outside`="Per-Out",
    `Nimes inside`="Ni-In",`Nimes outside`="Ni-Out",
    `Romans-Sur-Isere inside`= "Ro-In",`Romans-Sur-Isere outside`="Ro-Out",
    `Albentosa inside` = "Al-In" ,`Albentosa outside` = "Al-Out",
    `Moncayo inside`="Mon-In",`Moncayo outside`="Mon-Out",
    `Zuniga inside`="Zu-In",`Zuniga outside`="Zu-Out",
    `Acedo inside`= "Ac-In",`Acedo outside`="Ac-Out")) %>%
  filter(Genus != "Unlcassified" & Genus != "unclassified") %>% 
  mutate(Genus = as.factor(Genus),
         Genus = fct_relevel(
           Genus, 
           "Tuber","Scleroderma","Inocybe","Hebeloma","Tomentella",
           "Trichophaea","Tarzetta","Russula","Descomyces","Entoloma", 
           "Geopora","Wilcoxina","Cortinarius","Astraeus","Picoa",
           "Hydnobolites","Melanogaster","Clavulina","Genea","Boletus")) %>% 
  ggplot(aes(x = site_brule, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity") +
  themePlot() +
  theme(axis.title.x = element_blank(),
        axis.text.x =  element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        legend.text = element_text(face = "italic", size=7.5),
        legend.title = element_blank()) +
  guides(fill= guide_legend(ncol=1)) +
  scale_fill_manual(values = colors_bar_ecm)+
  scale_y_continuous(limits = c(-50,20000), expand = c(0, 0))

Fig_ECM_Abund

# Fungi
colors_bar_fungi <-c("#D21E2C","#ae09ea","#F7F7C5","#283dff","#014443",
                     "#b7ffdb","#FDDB8E","#a35151","#117744","#ea7f17",
                     "#dba4a4","#058ED9","#82807f","#bfc5ff", "#ffb7ef", 
                     "#ffb7ef","#ffff9e","#d8d6d4","#000000","#111b77")


Fig_fungi_Abund <-
  subset_taxa(physeq_fungi_all, Genus %in% c(
    ExtractTop(physeq_fungi_all, "Genus", 32) 
  )) %>% 
  tax_glom(., taxrank="Genus") %>% 
  psmelt() %>% 
  arrange(Genus) %>% 
  mutate(Genus = recode(Genus, "Gibberella" = "Fusarium" , 
                        "Cylindrocarpon" = "Ilyonectria")) %>% 
  filter(Genus != "Unclassified" & Genus != "unclassified") %>% 
  mutate(site_brule = recode_factor(
    site_brule,
    `Yarra inside` = "Ya-In" ,`Yarra outside` = "Ya-Out",
    `Wattles inside`="Wa-In",`Wattles outside`="Wa-Out",
    `Launceston inside` = "La-In" ,`Launceston outside` = "La-Out",
    `Needles inside`="Ne-In",`Needles outside`="Ne-Out",
    `Mole inside`="Mol-In",`Mole outside`="Mol-Out",
    `Camberra inside` = "Ca-In" ,`Camberra outside` = "Ca-Out",
    `Warri inside`="Wi-In",`Warri outside`="Wi-Out",
    `Braidwood inside`="Br-In",`Braidwood outside`="Br-Out",
    `Manjimup inside` = "Ma-In" ,`Manjimup outside` = "Ma-Out",
    `Pemberton inside`="Pem-In",`Pemberton outside`="Pem-Out",
    `Jardee inside`="Ja-In",`Jardee outside`="Ja-Out",
    `Cuneo inside` = "Cu-In" ,`Cuneo outside` = "Cu-Out",
    `San Demetrio inside`="SD-In",`San Demetrio outside`="SD-Out",
    `Spoleto inside`="Sp-In",`Spoleto outside`="Sp-Out",
    `Norcia inside`= "No-In",`Norcia outside`="No-Out",
    `Cognac inside` = "Co-In" ,`Cognac outside` = "Co-Out",
    `Grignan inside`="Gr-In",`Grignan outside`="Gr-Out",
    `Perpignan inside`="Per-In",`Perpignan outside`="Per-Out",
    `Nimes inside`="Ni-In",`Nimes outside`="Ni-Out",
    `Romans-Sur-Isere inside`= "Ro-In",`Romans-Sur-Isere outside`="Ro-Out",
    `Albentosa inside` = "Al-In" ,`Albentosa outside` = "Al-Out",
    `Moncayo inside`="Mon-In",`Moncayo outside`="Mon-Out",
    `Zuniga inside`="Zu-In",`Zuniga outside`="Zu-Out",
    `Acedo inside`= "Ac-In",`Acedo outside`="Ac-Out")) %>%
  filter(Genus != "Unlcassified" & Genus != "unclassified") %>% 
  mutate(Genus = as.factor(Genus),
         Genus = fct_relevel(
           Genus, 
           "Tuber","Fusarium","Mortierella","Cryptococcus","Scleroderma",
           "Inocybe","Hebeloma","Ilyonectria","Tomentella",
           "Metarhizium","Hygrocybe","Tetracladium","Trichophaea","Exophiala", 
           "Humicola", "Geminibasidium", "Archaeorhizomyces", "Tarzetta")) %>% 
  ggplot(aes(x = site_brule, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity") +
  themePlot() +
  theme(axis.title.x = element_blank(),
        axis.text.x =  element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        legend.text = element_text(face = "italic", size=7.5),
        legend.title = element_blank()) +
  guides(fill= guide_legend(ncol=1)) +
  scale_fill_manual(values = colors_bar_fungi) +
  scale_y_continuous(limits = c(-100,25000), expand = c(0, 0))

Fig_fungi_Abund


# Bacteria 
colors_bar_bact <- c("#D21E2C","#014443","#b7ffdb","#F7F7C5","#ea7f17",
                     "#bfc5ff","#111b77","#ae09ea","#FDDB8E","#283dff",
                     "#a35151","#117744","#dba4a4","#058ED9","#82807f",
                     "#ffb7ef","#ffff9e","#d8d6d4","#000000","#560d0d",
                     "#521899","#fcb067","#60ffaf","#825121","#fa7efc")

bact_dataset_bars <-
  subset_taxa(physeq_bact_all, Class %in% c(
    ExtractTop(physeq_bact_all, "Class", 25) 
  )) %>% 
  tax_glom(., taxrank="Class") %>% 
  psmelt() %>% 
  arrange(Class) %>% 
  mutate(site_brule = paste(site, brule, sep = " ")) %>% 
  mutate(site_brule = recode_factor(
    site_brule,
    `Yarra inside` = "Ya-In" ,`Yarra outside` = "Ya-Out",
    `Wattles inside`="Wa-In",`Wattles outside`="Wa-Out",
    `Launceston inside` = "La-In" ,`Launceston outside` = "La-Out",
    `Needles inside`="Ne-In",`Needles outside`="Ne-Out",
    `Mole inside`="Mol-In",`Mole outside`="Mol-Out",
    `Camberra inside` = "Ca-In" ,`Camberra outside` = "Ca-Out",
    `Warri inside`="Wi-In",`Warri outside`="Wi-Out",
    `Braidwood inside`="Br-In",`Braidwood outside`="Br-Out",
    `Manjimup inside` = "Ma-In" ,`Manjimup outside` = "Ma-Out",
    `Pemberton inside`="Pem-In",`Pemberton outside`="Pem-Out",
    `Jardee inside`="Ja-In",`Jardee outside`="Ja-Out",
    `Cuneo inside` = "Cu-In" ,`Cuneo outside` = "Cu-Out",
    `San Demetrio inside`="SD-In",`San Demetrio outside`="SD-Out",
    `Spoleto inside`="Sp-In",`Spoleto outside`="Sp-Out",
    `Norcia inside`= "No-In",`Norcia outside`="No-Out",
    `Cognac inside` = "Co-In" ,`Cognac outside` = "Co-Out",
    `Grignan inside`="Gr-In",`Grignan outside`="Gr-Out",
    `Perpignan inside`="Per-In",`Perpignan outside`="Per-Out",
    `Nimes inside`="Ni-In",`Nimes outside`="Ni-Out",
    `Romans-Sur-Isere inside`= "Ro-In",`Romans-Sur-Isere outside`="Ro-Out",
    `Albentosa inside` = "Al-In" ,`Albentosa outside` = "Al-Out",
    `Moncayo inside`="Mon-In",`Moncayo outside`="Mon-Out",
    `Zuniga inside`="Zu-In",`Zuniga outside`="Zu-Out",
    `Acedo inside`= "Ac-In",`Acedo outside`="Ac-Out")) %>%
  filter(Class != "Unlcassified" & Class != "unclassified") %>% 
  mutate(Class = as.factor(Class),
         Class = fct_relevel(
           Class, 
           "Alphaproteobacteria", "Acidobacteria-6", "Thermoleophilia", 
           "Actinobacteria", "Betaproteobacteria",
           "Ellin6529","Deltaproteobacteria", "[Chloracidobacteria]", 
           "Planctomycetia","Thaumarchaeota", 
           "Gammaproteobacteria", "Rubrobacteria", "Acidimicrobiia", 
           "Anaerolineae","[Saprospirae]",
           "Cytophagia","Bacilli","TK10","iii1-8","MB-A2-108",
           "Phycisphaerae", "[Spartobacteria]", "Gemmatimonadetes",
           "Thermomicrobia", "Nitrospira"))

bact_dataset_bars %>% 
  select(Abundance, site_brule, Class) %>% 
  group_by(site_brule) %>% 
  summarise(Abund = sum(Abundance)) %>% 
  as.data.frame()

Fig_bacteria_Abund <-
  bact_dataset_bars %>% 
  ggplot(aes(x = site_brule, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity") +
  themePlot() +
  theme(axis.title.x = element_blank(),
        axis.text.x =  element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        legend.text = element_text(face = "italic", size = 7.5),
        legend.title = element_blank(),
        legend.key.height = unit(0.3, "cm"), 
        legend.key.width = unit(0.3, "cm")) +
  guides(fill= guide_legend(ncol=1)) +
  scale_fill_manual(values = colors_bar_bact) +
  scale_y_continuous(limits = c(-100,40000), expand = c(0, 0)) 

Fig_bacteria_Abund

# Composite figure Fungi
Fig_barplot_all_fungi <-
  ggarrange(Fig_fungi_Abund + theme(legend.margin=ggplot2::margin(0,0,0,0),
                                    legend.justification = "left",
                                    legend.key.height = unit(0.25, "cm"), 
                                    legend.key.width = unit(0.25, "cm")) +
              annotate("text",  x=Inf, y = Inf, label = "Fungi", 
                       vjust=1.2, hjust=1.2),
            Fig_ECM_Abund + theme(legend.margin=ggplot2::margin(0,0,0,0),
                                  legend.justification = "left",
                                  legend.key.height = unit(0.25, "cm"), 
                                  legend.key.width = unit(0.25, "cm")) +
              annotate("text",  x=Inf, y = Inf, label = "Ectomycorrhizal fungi", 
                       vjust=1.2, hjust=1.05),
            Fig_Tabund + theme(legend.margin=ggplot2::margin(0,0,0,0),
                               legend.justification = "left",
                               legend.key.height = unit(0.25, "cm"), 
                               legend.key.width = unit(0.25, "cm"))+
              annotate("text",  x=Inf, y = Inf, label = "Genus~italic( Tuber)", 
                       parse=TRUE,
                       vjust=1.2, hjust=1.1),
            ncol = 1, 
            nrow = 3, 
            align = "hv", 
            labels = c("A", "B", "C"))

Fig_barplot_all_fungi


grid.arrange(
  Fig_barplot_all_fungi,
  top = text_grob(
    "FUNGAL AND BACTERIAL MICROBIOME COMPOSITION",
    size = 12,
    face = 2
  )
)

