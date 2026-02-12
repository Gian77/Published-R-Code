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
# Custom ggplot2 theme ---------------------------------------------------------
themePlot <- function(){ 
  require(ggtext)
  theme_bw() +
    theme(plot.title = element_markdown(size = 10, face = "bold", hjust = 0.5),
          strip.text = element_markdown(size = 10, face = "bold"),
          axis.text.x = element_markdown(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
          axis.text.y = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
          strip.background = element_blank(),
          legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm")) 
}

# Color palette ----------------------------------------------------------------
palette_cont <-
  c("#00843D", "#003399")

palette_cont_brule <-
  c("#77BC97", "#00843D", "#6685c2", "#003399")

tuber_colors <- c(
  "Tuber melanosporum" = "#D21E2C",
  "Tuber rufum f. nitidum" = "#560d0d",
  "Tuber melosporum" =  "#a35151",
  "Tuber oligospermum" =  "#ea7f17",
  "Tuber spp." =   "#82807f",
  "Tuber brumale" = "#111b77",
  "Tuber rufum" = "#283dff",
  "Tuber rufum f lucidum" =   "#058ED9",
  "Tuber rufum f ferrugineum" =   "#fa7efc",
  "Tuber mexiusanum" =   "#521899",
  "Tuber borchii" =   "#117744",
  "Tuber gennadii" = "#60ffaf",
  "Tuber lyonii" = "#825121",
  "Tuber aestivum" = "#000000"
)

fungi_colors <- c(
  "Tuber" = "#D21E2C",
  "Fusarium" =  "#560d0d",
  "Mortierella" = "#a35151",
  "Cryptococcus" = "#dba4a4",
  "Scleroderma" =   "#ffb7ef",
  "Inocybe" =  "#fa7efc",
  "Hebeloma" =  "#ae09ea",
  "Ilyonectria" = "#521899",
  "Tomentella" =   "#111b77",
  "Metarhizium" =  "#058ED9",
  "Hygrocybe" = "#114477",
  "Tetracladium" =  "#014443",
  "Trichophaea" = "#283dff",
  "Exophiala" =  "#ea7f17",
  "Humicola" = "#FDDB8E",
  "Geminibasidium" = "#ffff9e",
  "Archaeorhizomyces" = "#825121",
  "Tarzetta" = "#bfc5ff",
  "Cladosporium" = "#cccccc",
  "Haematonectria" = "#330033"
)

ecm_colors <- c(
  "Tuber" = "#D21E2C",
  "Scleroderma" = "#ffb7ef",
  "Inocybe" =  "#fa7efc",
  "Hebeloma" =  "#ae09ea",
  "Tomentella" =   "#111b77",
  "Trichophaea" =  "#283dff",
  "Tarzetta" =   "#bfc5ff",
  "Russula" =  "#117744",
  "Descolea" = "#60ffaf",
  "Entoloma" =  "#b7ffdb",
  "Geopora" = "#AAAA44",
  "Wilcoxina" = "#fcb067",
  "Cortinarius" =  "#F7F7C5",
  "Astraeus" =    "#82807f",
  "Picoa" = "#d8d6d4",
  "Hydnobolites" = "#333333",
  "Melanogaster" = "#000000",
  "Clavulina" = "#500050",
  "Genea" = "#662700",
  "Boletus" = "#331900"
)

colors_bar_bact <- c(
  "#D21E2C",
  "#014443",
  "#b7ffdb",
  "#F7F7C5",
  "#ea7f17",
  "#bfc5ff",
  "#111b77",
  "#ae09ea",
  "#FDDB8E",
  "#283dff",
  "#a35151",
  "#117744",
  "#dba4a4",
  "#058ED9",
  "#82807f",
  "#ffb7ef",
  "#ffff9e",
  "#d8d6d4",
  "#000000",
  "#560d0d",
  "#521899",
  "#fcb067",
  "#60ffaf",
  "#825121",
  "#fa7efc"
)

palette_site_gleba <- c(
  "#ea7f17",
  "#111b77",
  "#bfc5ff",
  "#ffb7ef",
  "#fa7efc",
  "#ae09ea",
  "#521899",
  "#014443",
  "#60ffaf",
  "#b7ffdb",
  "#d8d6d4",
  "#82807f",
  "#825121",
  "#000000"
)


colors_bar_gleba <-c("#D21E2C","#058ED9","#ae09ea","#dba4a4","#117744",
                     "#F7F7C5","#283dff","#521899","#82807f","#014443",
                     "#FDDB8E","#bfc5ff","#111b77","#d8d6d4","#b7ffdb",
                     "#fcb067","#ffb7ef","#000000","#560d0d","#60ffaf",
                     "#ea7f17","#fa7efc","#a35151","#825121","#A5A518")

# **************************************************----------------------------
# IMPORT -----------------------------------------------------------------------

load(file = "github/phyloseq_objects.RData")
ls()

# PATH -------------------------------------------------------------------------

project_dir <- 
  ("/home/gian/Dropbox/6_PROJETCS/2025_CompetitionRelease_Tuber_EuropeAustralia_AEMgithub/")

project_dir

# original
physeq_ITS <- readRDS(file=file.path(project_dir, "github/phyloseq_ITS.RDS"))
physeq_16s <- readRDS(file=file.path(project_dir, "github/phyloseq_16S.RDS"))

# rarefied
physeq_fungi_all <- readRDS(file = file.path(project_dir, "github/physeq_fungi_rare.RDS"))
physeq_ecm_all <- readRDS(file = file.path(project_dir, "github/physeq_ecm_rare.RDS"))
physeq_bact_all <- readRDS(file = file.path(project_dir, "github/physeq_bact_rare.RDS"))
physeq_gleba_ev <- readRDS(file = file.path(project_dir, "github/physeq_bact_truffle.RDS"))
physeq_bact_out_ev <- readRDS(file = file.path(project_dir, "github/physeq_bact_outgroup.RDS"))

# soil data
soil_data <- readRDS(file = file.path(project_dir, "github/soil_data.RDS"))

#df_alpha_fungi_all_new <- readRDS(file = file.path(project_dir, "github/df_alpha_fungi_all.RDS"))
#df_alpha_ecm_new <- readRDS(file = file.path(project_dir, "github/df_alpha_ecm_all.RDS"))
#df_alpha_bact_all_new <- readRDS(file = file.path(project_dir, "github/df_alpha_bact_all.RDS"))

# **************************************************----------------------------
# FIGURES and TABLES -----------------------------------------------------------

# **************************************************----------------------------
# Generate dataframes ----------------------------------------------------------

# Fungi ------------------------------------------------------------------------
head(physeq_fungi_all@sam_data)
table(physeq_fungi_all@sam_data$site, physeq_fungi_all@sam_data$brule)
table(physeq_fungi_all@sam_data$site_code,  physeq_fungi_all@sam_data$brule)

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
                           Parnans = "Romans-Sur-Isere")))

df_alpha_fungi_all_new <- 
  as.data.frame(as.matrix(sample_data(physeq_fungi_all))) %>% 
  mutate(TreeID = paste0(site_code, str_extract(SampleID, "\\d+"))) %>% 
  mutate(continent_brule = interaction(continent, brule, sep = "_")) %>% 
  mutate(site_brule = interaction(site, brule, sep = "_")) %>% 
  select(SampleID, TreeID, 
         continent, continent_brule, 
         site, site_brule, 
         brule, management,
         hill_0, hill_1, hill_2) %>% 
  mutate(hill_0= as.numeric(hill_0), 
         hill_1= as.numeric(hill_1),
         hill_2= as.numeric(hill_2)) %>%
  mutate(TreeID = recode(TreeID, EP219 = "EP290", EP288 = "EP298" )) %>% 
  janitor::clean_names() %>% 
  mutate(site =  fct_relevel(site, 
                             "Yarra","Wattles","Launceston","Needles","Mole",
                             "Camberra","Warri","Braidwood","Manjimup",
                             "Pemberton","Jardee", 
                             "Cuneo","San Demetrio","Spoleto","Norcia",
                             "Cognac","Grignan","Perpignan","Nimes","Romans-Sur-Isere",
                             "Albentosa", "Moncayo", "Zuniga","Acedo")) %>% 
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
  mutate_at(c("hill_0", "hill_1", "hill_2"), as.numeric) %>% 
  filter(management %in% c("cultivated")) %>%
  droplevels() 

unique(df_alpha_fungi_all_new$site)
head(df_alpha_fungi_all_new)
table(df_alpha_fungi_all_new$tree_id)
table(df_alpha_fungi_all_new$brule, df_alpha_fungi_all_new$tree_id)
table(df_alpha_fungi_all_new$site, df_alpha_fungi_all_new$brule)


# Ectomycorrhizal fungi --------------------------------------------------------
head(physeq_ecm_all@sam_data)
table(physeq_ecm_all@sam_data$site, physeq_ecm_all@sam_data$brule)
table(physeq_ecm_all@sam_data$site_code,  physeq_ecm_all@sam_data$brule)

sample_data(physeq_ecm_all) <-
  sample_data(
    as.data.frame(as.matrix(physeq_ecm_all@sam_data)) %>% 
      mutate(site = recode(site,
                           S.Demetrio = "San Demetrio",
                           Uzes = "Nimes",
                           Manjimup2 = "Jardee",
                           Chateauvert = "Grignan",
                           Angouleme = "Cognac",
                           `Espira-de-l-Agly` = "Perpignan",
                           Parnans = "Romans-Sur-Isere")))

df_alpha_ecm_new <- 
  as.data.frame(as.matrix(sample_data(physeq_ecm_all))) %>% 
  mutate(tree_id = paste0(site_code, str_extract(SampleID, "\\d+"))) %>% 
  mutate(continent_brule = interaction(continent, brule, sep = "_")) %>% 
  mutate(site_brule = interaction(site, brule, sep = "_")) %>% 
  select(SampleID, tree_id, 
         continent, continent_brule, 
         site, site_brule, 
         brule, management,
         hill_0, hill_1, hill_2) %>% 
  mutate(hill_0= as.numeric(hill_0), 
         hill_1= as.numeric(hill_1),
         hill_2= as.numeric(hill_2)) %>%
  mutate(tree_id = recode(tree_id, EP219 = "EP290", EP288 = "EP298" )) %>% 
  janitor::clean_names() %>% 
  mutate(site =  fct_relevel(site, 
                             "Yarra","Wattles","Launceston","Needles","Mole",
                             "Camberra","Warri","Braidwood","Manjimup",
                             "Pemberton","Jardee", 
                             "Cuneo","San Demetrio","Spoleto","Norcia",
                             "Cognac","Grignan","Perpignan","Nimes","Romans-Sur-Isere",
                             "Albentosa", "Moncayo", "Zuniga","Acedo")) %>% 
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
  mutate_at(c("hill_0", "hill_1", "hill_2"), as.numeric) %>% 
  filter(management %in% c("cultivated")) %>%
  droplevels() 

unique(df_alpha_ecm_new$site)
head(df_alpha_ecm_new)
table(df_alpha_ecm_new$tree_id)
table(df_alpha_ecm_new$brule, df_alpha_ecm_new$tree_id)
table(df_alpha_ecm_new$site, df_alpha_ecm_new$brule)


# Bacteria ---------------------------------------------------------------------
head(physeq_bact_all@sam_data)
table(physeq_bact_all@sam_data$site, physeq_bact_all@sam_data$brule)
table(physeq_bact_all@sam_data$site_code,  physeq_bact_all@sam_data$brule)

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
                           Parnans = "Romans-Sur-Isere")))

df_alpha_bact_all_new <- 
  as.data.frame(as.matrix(sample_data(physeq_bact_all))) %>% 
  mutate(
    # extract the tree number right before N/F/B/HB at end of SampleID
    tree_num = str_extract(SampleID, "\\d+(?=(N|F|B|HB)$)"),
    tree_id   = paste0(site_code, tree_num)
  ) %>% 
  mutate(continent_brule = interaction(continent, brule, sep = "_")) %>% 
  mutate(site_brule = interaction(site, brule, sep = "_")) %>% 
  select(SampleID, tree_id, 
         continent, continent_brule, 
         site, site_brule, 
         brule, management,
         hill_0, hill_1, hill_2) %>% 
  mutate(hill_0= as.numeric(hill_0), 
         hill_1= as.numeric(hill_1),
         hill_2= as.numeric(hill_2)) %>%
  mutate(tree_id = recode(tree_id, EP219 = "EP290", EP288 = "EP298"),
         tree_id = recode(tree_id, EH5 = "EH1")) %>% 
  janitor::clean_names() %>% 
  mutate(site =  fct_relevel(site, 
                             "Yarra","Wattles","Launceston","Needles","Mole",
                             "Camberra","Warri","Braidwood","Manjimup",
                             "Pemberton","Jardee", 
                             "Cuneo","San Demetrio","Spoleto","Norcia",
                             "Cognac","Grignan","Perpignan","Nimes","Romans-Sur-Isere",
                             "Albentosa", "Moncayo", "Zuniga","Acedo")) %>% 
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
    `Acedo inside`= "Ac-In",`Acedo outside`="Ac-Out"))%>% 
  mutate_at(c("hill_0", "hill_1", "hill_2"), as.numeric) %>% 
  filter(management %in% "cultivated") %>% 
  droplevels() 


head(df_alpha_bact_all_new)
table(df_alpha_bact_all_new$tree_id) 
table(df_alpha_bact_all_new$brule, df_alpha_bact_all_new$tree_id)
table(df_alpha_bact_all_new$site, df_alpha_bact_all_new$brule)

# **************************************************----------------------------
# Plot titles ------------------------------------------------------------------
title1 = text_grob("Fungi", size = 12, face = 2)
title2 = text_grob("Ectomycorrhizal fungi", size = 12, face = 2)
title3 = text_grob("Bacteria", size = 12, face = 2)
title4 = text_grob("Tuber melanosporum", size = 12, face = 4)

# Plot labels ------------------------------------------------------------------

# Continent labels
continent_labels <- c(
  Australia = "<span style='color:#00843D'>Australia</span>",
  Europe = "<span style='color:#003399'>Europe</span>")


# Continent x bule labels
cont_brule_labels <- c(
  AI = "<span style='color:#00843D'>AI</span>",
  AO = "<span style='color:#00843D'>AO</span>",
  EI = "<span style='color:#003399'>EI</span>",
  EO = "<span style='color:#003399'>EO</span>"
)

#cont_brule_labels <- c(
#  `Australia<br>Inside` = "<span style='color:#00843D'>Australia<br>Inside</span>",
#  `Australia<br>Outside` = "<span style='color:#00843D'>Australia<br>Outside</span>",
#  `Europe<br>Inside` = "<span style='color:#003399'>Europe<br>Inside</span>",
#  `Europe<br>Outside` = "<span style='color:#003399'>Europe<br>Outside</span>"
#)

# Site labels ------------------------------------------------------------------
make_site_labels <- function(df = df_alpha_fungi_all_new){
  site_levels <- levels(df$site_brule)
  
  idx_break <- match("Ja-Out", site_levels)
  
  site_colors <- ifelse(seq_along(site_levels) <= idx_break,
                        "#00843D",  # green
                        "#003399")  # blue
  names(site_colors) <- site_levels
  
  gray_levels <- c("Ac-In","Ac-Out","Ro-In","Ro-Out","No-In","No-Out")
  site_colors[intersect(gray_levels, names(site_colors))] <- "#808080"
  
  site_labels <- purrr::map_chr(
    site_levels,
    ~ sprintf("<span style='color:%s'>%s</span>", site_colors[.x], .x)
  )
  names(site_labels) <- site_levels
  
  site_labels  # <- return this
}

site_labels <- make_site_labels()
site_labels


# **************************************************----------------------------

# **** FIGURE 4 **** -----------------------------------------------------------

# 1) First Approach - averages -------------------------------------------------
# Calculating within site averages to overcome pseudoreplication 

alpha_continent_brule_mean <-
  ggarrange(
    ggarrange(df_alpha_fungi_all_new %>% 
                group_by(site, brule, continent) %>%
                summarise(mean = mean(hill_0), .groups = "drop") %>% 
                mutate(continent_brule = paste(brule, continent, sep="_")) %>% 
                mutate(continent_brule = as.factor(continent_brule),
                       continent_brule = recode(continent_brule,
                                                "outside_Australia"= "AO", 
                                                "inside_Australia" = "AI", 
                                                "inside_Europe" = "EI", 
                                                "outside_Europe" = "EO"), 
                       continent_brule = fct_relevel(continent_brule, "AI","AO", "EI", "EO")) %>% 
                PlotRich(., "continent_brule", "mean",
                         CompSampl(.,  formula(mean ~ continent_brule)) %>% 
                           pull(Letters) %>% as.character(), 400) + 
                labs(title = "Richness (hill_0)", y=NULL) +
                scale_color_manual(values = palette_cont_brule)+ 
                scale_shape_manual(values = c(1,0,2,5)) +
                theme(
                  axis.title.y = element_markdown(size = 10),
                  axis.text.x = element_markdown(size = 10, hjust = 0.5, vjust = 0.5)) +
                scale_x_discrete(labels = cont_brule_labels),
              
              df_alpha_fungi_all_new %>% 
                group_by(site, brule, continent) %>%
                summarise(mean = mean(hill_2), .groups = "drop") %>% 
                mutate(continent_brule = paste(brule, continent, sep="_")) %>% 
                mutate(continent_brule = as.factor(continent_brule),
                       continent_brule = recode(continent_brule,
                                                "outside_Australia"= "AO", 
                                                "inside_Australia" = "AI", 
                                                "inside_Europe" = "EI", 
                                                "outside_Europe" = "EO"), 
                       continent_brule = fct_relevel(continent_brule, "AI","AO", "EI", "EO")) %>% 
                PlotRich(., "continent_brule", "mean",
                         CompSampl(.,  formula(mean ~ continent_brule)) %>% 
                           pull(Letters) %>% as.character(), 50) + 
                labs(title = "invSimpson (hill_2)", y=NULL) +
                scale_color_manual(values = palette_cont_brule)+ 
                scale_shape_manual(values = c(1,0,2,5)) +
                theme( axis.title.y = element_markdown(size = 10),
                       axis.text.x = element_markdown(size = 10, hjust = 0.5, vjust = 0.5))+
                scale_x_discrete(labels = cont_brule_labels),
              nrow = 1,
              ncol = 2, 
              labels = c("D")),
    ggarrange(df_alpha_ecm_new %>% 
                group_by(site, brule, continent) %>%
                summarise(mean = mean(hill_0), .groups = "drop") %>% 
                mutate(continent_brule = paste(brule, continent, sep="_")) %>% 
                mutate(continent_brule = as.factor(continent_brule),
                       continent_brule = recode(continent_brule,
                                                "outside_Australia"= "AO", 
                                                "inside_Australia" = "AI", 
                                                "inside_Europe" = "EI", 
                                                "outside_Europe" = "EO"), 
                       continent_brule = fct_relevel(continent_brule, "AI","AO", "EI", "EO")) %>% 
                PlotRich(., "continent_brule", "mean",
                         CompSampl(.,  formula(mean ~ continent_brule)) %>% 
                           pull(Letters) %>% as.character(), 24)+ 
                labs(title = "Richness (hill_0)", y=NULL) +
                scale_color_manual(values = palette_cont_brule)+ 
                scale_shape_manual(values = c(1,0,2,5)) +
                theme( axis.title.y = element_markdown(size = 10),
                       axis.text.x = element_markdown(size = 10, hjust = 0.5, vjust = 0.5))+
                scale_x_discrete(labels = cont_brule_labels),
              df_alpha_ecm_new %>% 
                group_by(site, brule, continent) %>%
                summarise(mean = mean(hill_2), .groups = "drop") %>% 
                mutate(continent_brule = paste(brule, continent, sep="_")) %>% 
                mutate(continent_brule = as.factor(continent_brule),
                       continent_brule = recode(continent_brule,
                                                "outside_Australia"= "AO", 
                                                "inside_Australia" = "AI", 
                                                "inside_Europe" = "EI", 
                                                "outside_Europe" = "EO"), 
                       continent_brule = fct_relevel(continent_brule, "AI","AO", "EI", "EO")) %>% 
                PlotRich(., "continent_brule", "mean",
                         CompSampl(.,  formula(mean ~ continent_brule)) %>% 
                           pull(Letters) %>% as.character(), 5) + 
                labs(title = "invSimpson (hill_2)", y=NULL) +
                scale_color_manual(values = palette_cont_brule)+ 
                scale_shape_manual(values = c(1,0,2,5)) +
                theme(axis.title.y = element_markdown(size = 10),
                      axis.text.x = element_markdown(size = 10, hjust = 0.5, vjust = 0.5))+
                scale_x_discrete(labels = cont_brule_labels),
              nrow = 1,
              ncol = 2, 
              labels = c("E")),
    ggarrange(df_alpha_bact_all_new %>% 
                group_by(site, brule, continent) %>%
                summarise(mean = mean(hill_0), .groups = "drop") %>% 
                mutate(continent_brule = paste(brule, continent, sep="_")) %>% 
                mutate(continent_brule = as.factor(continent_brule),
                       continent_brule = recode(continent_brule,
                                                "outside_Australia"= "AO", 
                                                "inside_Australia" = "AI", 
                                                "inside_Europe" = "EI", 
                                                "outside_Europe" = "EO"), 
                       continent_brule = fct_relevel(continent_brule, "AI","AO", "EI", "EO")) %>% 
                PlotRich(., "continent_brule", "mean",
                         CompSampl(.,  formula(mean ~ continent_brule)) %>% 
                           pull(Letters) %>% as.character(), 1500) + 
                labs(title = "Richness (hill_0)", y=NULL) +
                scale_color_manual(values = palette_cont_brule)+ 
                scale_shape_manual(values = c(1,0,2,5)) +
                theme(axis.title.y = element_markdown(size = 10),
                      axis.text.x = element_markdown(size = 10, hjust = 0.5, vjust = 0.5))+
                scale_x_discrete(labels = cont_brule_labels),
              df_alpha_bact_all_new %>% 
                group_by(site, brule, continent) %>%
                summarise(mean = mean(hill_2), .groups = "drop") %>% 
                mutate(continent_brule = paste(brule, continent, sep="_")) %>% 
                mutate(continent_brule = as.factor(continent_brule),
                       continent_brule = recode(continent_brule,
                                                "outside_Australia"= "AO", 
                                                "inside_Australia" = "AI", 
                                                "inside_Europe" = "EI", 
                                                "outside_Europe" = "EO"), 
                       continent_brule = fct_relevel(continent_brule, "AI","AO", "EI", "EO")) %>% 
                PlotRich(., "continent_brule", "mean",
                         CompSampl(.,  formula(mean ~ continent_brule)) %>% 
                           pull(Letters) %>% as.character(), 300) + 
                labs(title = "invSimpson (hill_2)", y=NULL) +
                scale_color_manual(values = palette_cont_brule) + 
                scale_shape_manual(values = c(1,0,2,5)) +
                theme(axis.title.y = element_markdown(size = 10),
                      axis.text.x = element_markdown(size = 10,hjust = 0.5, vjust = 0.5))+
                scale_x_discrete(labels = cont_brule_labels),
              nrow = 1,
              ncol = 2, 
              labels = c("F")),
    ncol = 3, 
    nrow = 1)

alpha_continent_brule_mean

# 2) Second Approach - linear mixed effect model -------------------------------
# hill_0 models ----------------------------------------------------------------

# NOTE. To obtain p-values and significance levels directly in your ANOVA table for a 
# linear mixed model, you must use the lmerTest package instead of base lme4

# All fungi --------------------------------------------------------------------

# checking hill_0 distribution

# NOTE! I will transform to log scale to better achieve normality!

df_alpha_fungi_all_new %>% 
  ggplot(aes(x = log(hill_0))) +
  geom_histogram(binwidth = 0.05) +
  labs(title = "ITS hill_0") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(log(hill_0))), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)

# Outlier on the low end, nothing major

# Build linear mixed-effects models
# Dependent variable: log-transformed richness
# Fixed effects:      continent, brûlé, and their interaction
# Random effect:      site

model_rich_fungi <- lmer(log(hill_0) ~ continent * brule + (1 | site), 
                   data = df_alpha_fungi_all_new)
summary(model_rich_fungi)

# Trying using a random slope but does not work... It's ok to assume identical
# effect of site across sites.
lmer(hill_0 ~ continent * brule + (1 + brule | site/tree_id),
     data = df_alpha_fungi_all_new)

model_rich_fungi_2 <- lmer(log(hill_0) ~ continent * brule + (1 | site/tree_id), 
                     data = df_alpha_fungi_all_new)
summary(model_rich_fungi_2)

anova(model_rich_fungi, model_rich_fungi_2)

anova(
  update(model_rich_fungi, REML = FALSE),
  update(model_rich_fungi_2, REML = FALSE)
)

anova(model_rich_fungi, type = "1")
anova(model_rich_fungi, type = "2")
anova(model_rich_fungi, type = "3")

# Diagnostics 
sim_fungi <- simulateResiduals(model_rich_fungi)

par(mfrow = c(2, 2))
# KS test for correct distribution of residuals
testUniformity(sim_fungi)
# KS test for correct distribution within and between groups
testCategorical(sim_fungi, df_alpha_fungi_all_new$continent_brule)
# Dispersion test - for details see ?testDispersion
testDispersion(sim_fungi) # tests under and overdispersion
# Outlier test (number of observations outside simulation envelope)
testOutliers(sim_fungi, type = "bootstrap")
par(mfrow = c(1, 1)) 

# ECM fungi --------------------------------------------------------------------
df_alpha_ecm_new %>% 
  ggplot(aes(x = log(hill_0))) +
  geom_histogram(binwidth = 0.5) +
  labs(title = "ITS hill_0") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(log(hill_0))), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)

# Perfect!

model_rich_ecm <- lmer(log(hill_0) ~ continent * brule + (1 | site),
                       data = df_alpha_ecm_new)
summary(model_rich_ecm)

model_rich_ecm_2 <- lmer(log(hill_0) ~ continent * brule + (1 | site/tree_id), 
                     data = df_alpha_ecm_new)
summary(model_rich_ecm_2)

anova(model_rich_ecm, model_rich_ecm_2)

anova(
  update(model_rich_ecm, REML = FALSE),
  update(model_rich_ecm_2, REML = FALSE)
)

# Diagnostics 
sim_ecm <- simulateResiduals(model_rich_ecm_2)

par(mfrow = c(2, 2))
# KS test for correct distribution of residuals
testUniformity(sim_ecm)
# KS test for correct distribution within and between groups
testCategorical(sim_ecm, df_alpha_ecm_new$continent_brule)
# Dispersion test - for details see ?testDispersion
testDispersion(sim_ecm) # tests under and overdispersion
# Outlier test (number of observations outside simulation envelope)
testOutliers(sim_ecm, type = "bootstrap")
par(mfrow = c(1, 1)) 

# Bacteria ---------------------------------------------------------------------
df_alpha_bact_all_new %>% 
  ggplot(aes(x = log(hill_0))) +
  geom_histogram(binwidth = 0.05) +
  labs(title = "Bacteria hill_0") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(log(hill_0))), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)

# Seems there is a few outliears but not major 

model_rich_bact <- lmer(log(hill_0) ~ continent * brule + (1 | site), 
                        data = df_alpha_bact_all_new)
summary(model_rich_bact)

model_rich_bact_2 <- lmer(log(hill_0) ~ continent * brule + (1 | site/tree_id), 
                     data = df_alpha_bact_all_new)
summary(model_rich_bact_2)

anova(model_rich_bact, model_rich_bact_2)

anova(
  update(model_rich_bact, REML = FALSE),
  update(model_rich_bact_2, REML = FALSE)
)

# Diagnostics 
sim_bact <- simulateResiduals(model_rich_bact)

par(mfrow = c(2, 2))
# KS test for correct distribution of residuals
testUniformity(sim_bact)
# KS test for correct distribution within and between groups
testCategorical(sim_bact, df_alpha_bact_all_new$continent_brule)
# Dispersion test - for details see ?testDispersion
testDispersion(sim_bact) # tests under and overdispersion
# Outlier test (number of observations outside simulation envelope)
testOutliers(sim_bact, type = "bootstrap")
par(mfrow = c(1, 1)) 

# hill_2 models ----------------------------------------------------------------

# All fungi --------------------------------------------------------------------
df_alpha_fungi_all_new %>% 
  ggplot(aes(x = log(hill_1))) +
  geom_histogram(binwidth = 0.3) +
  labs(title = "ITS hill_2") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(log(hill_1))), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)

model_Shan_fungi <- lmer(log(hill_1) ~ continent * brule + (1 | site),
                       data = df_alpha_fungi_all_new)

summary(model_Shan_fungi)

model_Shan_fungi_2 <- lmer(log(hill_1) ~ continent * brule + (1 | site/tree_id), 
                            data = df_alpha_fungi_all_new)
summary(model_Shan_fungi_2)

anova(
  update(model_Shan_fungi, REML = FALSE),
  update(model_Shan_fungi_2, REML = FALSE)
)

# Diagnostics 
sim_invSim_fungi <- simulateResiduals(model_Shan_fungi)

par(mfrow = c(2, 2))
# KS test for correct distribution of residuals
testUniformity(sim_invSim_fungi)
# KS test for correct distribution within and between groups
testCategorical(sim_invSim_fungi, catPred = df_alpha_fungi_all_new$continent_brule)
# Dispersion test - for details see ?testDispersion
testDispersion(sim_invSim_fungi) # tests under and overdispersion
# Outlier test (number of observations outside simulation envelope)
testOutliers(sim_invSim_fungi, type = "bootstrap")
par(mfrow = c(1, 1)) 

# ECM fungi --------------------------------------------------------------------

df_alpha_ecm_new %>% 
  ggplot(aes(x = log(hill_1))) +
  geom_histogram(bins = 20) +
  labs(title = "ITS invSimpson") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(log(hill_1))), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)


model_Shan_ecm <- lmer(log(hill_1) ~ continent * brule + (1 | site),
                       data = df_alpha_ecm_new)
summary(model_Shan_ecm)


model_Shan_ecm_2 <- lmer(log(hill_1) ~ continent * brule + (1 | site/tree_id), 
                     data = df_alpha_ecm_new)
summary(model_Shan_ecm_2)

anova(
  update(model_Shan_ecm, REML = FALSE),
  update(model_Shan_ecm_2, REML = FALSE)
)

# Diagnostics 
sim_invSim_ecm <- simulateResiduals(model_Shan_ecm)

par(mfrow = c(2, 2))
# KS test for correct distribution of residuals
testUniformity(model_Shan_ecm)
# KS test for correct distribution within and between groups
testCategorical(model_Shan_ecm, df_alpha_ecm_new$continent_brule)
# Dispersion test - for details see ?testDispersion
testDispersion(model_Shan_ecm) # tests under and overdispersion
# Outlier test (number of observations outside simulation envelope)
testOutliers(model_Shan_ecm, type = "bootstrap")
par(mfrow = c(1, 1)) 

# Bacteria ---------------------------------------------------------------------

df_alpha_bact_all_new %>% 
  ggplot(aes(x = log(hill_1))) +
  geom_histogram(bins = 20) +
  labs(title = "ITS hill_1") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(log(hill_1))), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)

model_Shan_bact <- lmer(log(hill_1) ~ continent * brule + (1 | site),
                          data = df_alpha_bact_all_new)
summary(model_Shan_bact)

model_Shan_bact_2 <- lmer(log(hill_1) ~ continent * brule + (1 | site/tree_id), 
                            data = df_alpha_bact_all_new)
summary(model_Shan_bact_2)

anova(
  update(model_Shan_bact, REML = FALSE),
  update(model_Shan_bact_2, REML = FALSE)
)

# Diagnostics 
sim_invSim_bact <- simulateResiduals(model_Shan_bact)

par(mfrow = c(2, 2))
# KS test for correct distribution of residuals
testUniformity(sim_invSim_bact)
# KS test for correct distribution within and between groups
testCategorical(sim_invSim_bact, df_alpha_bact_all_new$continent_brule)
# Dispersion test - for details see ?testDispersion
testDispersion(sim_invSim_bact) # tests under and overdispersion
# Outlier test (number of observations outside simulation envelope)
testOutliers(model_Shan_ecm, type = "bootstrap")
par(mfrow = c(1, 1))


# Figure S15 -------------------------------------------------------------------
par(mfrow = c(3, 8))
# fungi
testUniformity(sim_fungi) +
  mtext("A", side = 3, line = 1, adj = 0, font = 2, cex = 1.2)
testCategorical(sim_fungi, df_alpha_fungi_all_new$continent_brule)
testDispersion(sim_fungi)+ 
testOutliers(sim_fungi, type = "bootstrap")
testUniformity(sim_invSim_fungi)
  mtext("B", side = 3, line = 1, adj = 0, font = 2, cex = 1.2)
testCategorical(sim_invSim_fungi, df_alpha_fungi_all_new$continent_brule)
testDispersion(sim_invSim_fungi) 
testOutliers(sim_invSim_fungi, type = "bootstrap")
# ecm
testUniformity(sim_ecm)+
  mtext("C", side = 3, line = 1, adj = 0, font = 2, cex = 1.2)
testCategorical(sim_ecm, df_alpha_ecm_new$continent_brule)
testDispersion(sim_ecm) 
testOutliers(sim_ecm, type = "bootstrap")
testUniformity(sim_invSim_ecm)+
  mtext("D", side = 3, line = 1, adj = 0, font = 2, cex = 1.2)
testCategorical(sim_invSim_ecm, df_alpha_ecm_new$continent_brule)
testDispersion(sim_invSim_ecm) 
testOutliers(sim_invSim_ecm, type = "bootstrap")
# bacteria
testUniformity(sim_bact)+
  mtext("E", side = 3, line = 1, adj = 0, font = 2, cex = 1.2)
testCategorical(sim_bact, df_alpha_bact_all_new$continent_brule)
testDispersion(sim_bact) 
testOutliers(sim_bact, type = "bootstrap")
testUniformity(sim_invSim_bact) +
  mtext("F", side = 3, line = 1, adj = 0, font = 2, cex = 1.2)
testCategorical(sim_invSim_bact, df_alpha_bact_all_new$continent_brule)
testDispersion(sim_invSim_bact) 
testOutliers(sim_invSim_bact, type = "bootstrap")

par(mfrow = c(1, 1))


# EULER DIAGRAMS ---------------------------------------------------------------

plot_euler_diag <- ggarrange(
  PlotVenn(physeq_fungi_all, "continent"),
  PlotVenn(physeq_ecm_all, "continent"),
  PlotVenn(physeq_bact_all, "continent"),
  ncol = 3, 
  nrow = 1,
  labels=c("G", "H", "I")
)

plot_euler_diag

# Final figure -----------------------------------------------------------------

alpha_continent_all <- 
  ggarrange(
    ggarrange(
    emm_plot_sig_h(model_rich_fungi, title = "All fungi", type = "response", adjust = "none")+
      theme(axis.title.x = element_blank()),
    emm_plot_sig_h(model_rich_ecm, title = "Ectomycorrhizal fungi", type = "response", adjust = "none")+
      theme(axis.title.x = element_blank()),
    emm_plot_sig_h(model_rich_bact, title = "Bacteria", type = "response", adjust = "none")+
      theme(axis.title.x = element_blank()),
    emm_plot_sig_h(model_Shan_fungi, title = "All fungi", type = "response", adjust = "none") +
      theme(plot.title = element_blank()),
    emm_plot_sig_h(model_Shan_ecm, title = "Ectomycorrhizal fungi", type = "response", adjust = "none")+
      theme(plot.title = element_blank()),
    emm_plot_sig_h(model_rich_bact, title = "Bacteria", type = "response", adjust = "none")+
      theme(plot.title = element_blank()),
    ncol = 3,
    nrow = 2,
    labels = c("A", "B", "C", "D", "E", "F")),
    plot_euler_diag,
ncol = 1,
nrow = 2
)

alpha_continent_all

ggpubr::annotate_figure(
  alpha_continent_all,
  top = text_grob("ALPHA DIVERSITY", 
                  size = 12, 
                  face = "bold"))

ggsave(
  "alphadiv.pdf",
  plot = ggpubr::annotate_figure(
    alpha_continent_all,
    top = text_grob("ALPHA DIVERSITY", size = 12, face = "bold")
  ),
  device = "pdf"
)

# **************************************************----------------------------
# Table S5 ---------------------------------------------------------------------
# Model Summary ----------------------------------------------------------------

summary(model_rich_fungi)
pairs(emmeans(model_rich_fungi, ~ continent | brule, type ="response"))
cld(emmeans(model_rich_fungi, ~ continent | brule, type = "response"), 
    by = "brule", Letters = letters, adjust = "BH")

summary(model_rich_ecm)
pairs(emmeans(model_rich_ecm, ~ continent | brule, type ="response"))
cld(emmeans(model_rich_ecm, ~ continent | brule, type = "response"), 
    by = "brule", Letters = letters, adjust = "BH")

summary(model_rich_bact)
pairs(emmeans(model_rich_bact, ~ continent | brule, type ="response"))
cld(emmeans(model_rich_bact, ~ continent | brule, type = "response"), 
    by = "brule", Letters = letters, adjust = "BH")

summary(model_Shan_fungi)
pairs(emmeans(model_Shan_fungi, ~ continent | brule, type ="response"))
cld(emmeans(model_Shan_fungi, ~ continent | brule, type = "response"), 
    by = "brule", Letters = letters, adjust = "BH")

summary(model_Shan_ecm)
pairs(emmeans(model_Shan_ecm, ~ continent | brule, type ="response"))
cld(emmeans(model_Shan_ecm, ~ continent | brule, type = "response"), 
    by = "brule", Letters = letters, adjust = "BH")

summary(model_Shan_bact)
pairs(emmeans(model_Shan_bact, ~ continent | brule, type ="response"))
cld(emmeans(model_Shan_bact, ~ continent | brule, type = "response"), 
    by = "brule", Letters = letters, adjust = "BH")


# Should consider also the brule
cld(emmeans(model_rich_ecm, ~ continent | brule, type ="response"), 
    Letters = letters, adjust = "BH")
cld(emmeans(model_rich_ecm, ~ brule | continent, type ="response"), 
    Letters = letters, adjust = "BH")

cld(emmeans(model_Shan_ecm, ~ continent | brule, type ="response"), 
    Letters = letters, adjust = "BH")
cld(emmeans(model_Shan_ecm, ~ brule | continent, type ="response"),
    Letters = letters, adjust = "BH")






# **************************************************----------------------------

# **** FIGURE 5 **** -----------------------------------------------------------

melt_Tuber_abund <-
  physeq_fungi_all %>% 
  subset_taxa(Genus=="Tuber") %>% 
  psmelt() %>%                                         
  arrange(Genus) %>% 
  mutate(Species = recode_factor(Species, 
                                 `Tuber sp` = "Tuber spp.",
                                 Unlcassified = "Tuber spp.")) %>% 
  mutate(Taxonomy = paste(OTU, Species, sep = " "))

head(melt_Tuber_abund)

# Recreate the Site_Brule factor exactly as in your pipeline
melt_Tuber_abund2 <- 
  melt_Tuber_abund %>% 
  mutate(
    site_brule =
      fct_relevel(
        site,
        "Yarra","Wattles","Launceston","Needles",
        "Mole","Camberra","Warri","Braidwood","Manjimup",
        "Pemberton","Jardee",
        "Cuneo","San Demetrio","Spoleto","Norcia","Cognac",
        "Grignan","Perpignan","Nimes","Romans-Sur-Isere",
        "Albentosa","Moncayo","Zuniga","Acedo")) %>% 
  mutate(site_brule = paste(site, brule, sep = " ")) %>% 
  mutate(site_brule = recode_factor(
    site_brule,
    `Yarra inside` = "Ya-In" ,`Yarra outside` = "Ya-Out",
    `Wattles inside`="Wa-In",`Wattles outside`="Wa-Out",
    `Launceston inside` = "La-In" ,`Launceston outside`="La-Out",
    `Needles inside`="Ne-In",`Needles outside`="Ne-Out",
    `Mole inside`="Mol-In",`Mole outside`="Mol-Out",
    `Camberra inside` = "Ca-In" ,`Camberra outside`="Ca-Out",
    `Warri inside`="Wi-In",`Warri outside`="Wi-Out",
    `Braidwood inside`="Br-In",`Braidwood outside`="Br-Out",
    `Manjimup inside` = "Ma-In" ,`Manjimup outside`="Ma-Out",
    `Pemberton inside`="Pem-In",`Pemberton outside`="Pem-Out",
    `Jardee inside`="Ja-In",`Jardee outside`="Ja-Out",
    `Cuneo inside` = "Cu-In" ,`Cuneo outside`="Cu-Out",
    `San Demetrio inside`="SD-In",`San Demetrio outside`="SD-Out",
    `Spoleto inside`="Sp-In",`Spoleto outside`="Sp-Out",
    `Norcia inside`= "No-In",`Norcia outside`="No-Out",
    `Cognac inside` = "Co-In" ,`Cognac outside`="Co-Out",
    `Grignan inside`="Gr-In",`Grignan outside`="Gr-Out",
    `Perpignan inside`="Per-In",`Perpignan outside`="Per-Out",
    `Nimes inside`="Ni-In",`Nimes outside`="Ni-Out",
    `Romans-Sur-Isere inside`= "Ro-In",`Romans-Sur-Isere outside`="Ro-Out",
    `Albentosa inside` = "Al-In" ,`Albentosa outside`="Al-Out",
    `Moncayo inside`="Mon-In",`Moncayo outside`="Mon-Out",
    `Zuniga inside`="Zu-In",`Zuniga outside`="Zu-Out",
    `Acedo inside`= "Ac-In",`Acedo outside`="Ac-Out")) %>% 
  select(OTU, SampleID, Abundance,
         continent, site, site_brule, brule, management,
         Kingdom, Phylum, Class, Order, Family, Genus, Species,Taxonomy) %>% 
  janitor::clean_names()

head(melt_Tuber_abund2)


# C) Truffle species stacked bars -------------------------------------------------
Fig_Tabund_revised <-
  melt_Tuber_abund2 %>% 
  mutate(species = as.factor(species),
         species = fct_relevel(species, 
                               "Tuber melanosporum","Tuber rufum f. nitidum",
                               "Tuber melosporum","Tuber oligospermum", 
                               "Tuber spp.","Tuber brumale","Tuber rufum",
                               "Tuber rufum f lucidum", "Tuber rufum f ferrugineum",
                               "Tuber mexiusanum","Tuber borchii","Tuber gennadii",
                               "Tuber lyonii", "Tuber aestivum")) %>% 
  ggplot(aes(x = site_brule, y = abundance, fill = species)) + 
  geom_bar(stat = "identity") +
  themePlot() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_markdown(size = 10,angle = 90, hjust = 1, vjust = 0.5),
    legend.text  = element_markdown(size = 9)
  ) +
  scale_x_discrete(labels = site_labels) +
  scale_fill_manual(
    "Species",
    values = tuber_colors, 
    labels = c(
      "*Tuber melanosporum*",
      "*Tuber rufum* f. *nitidum*",
      "*Tuber melosporum*",
      "*Tuber oligospermum*", 
      "*Tuber* spp.",
      "*Tuber brumale*",
      "*Tuber rufum*",
      "*Tuber rufum* f. *lucidum*", 
      "*Tuber rufum* f. *ferrugineum*",
      "*Tuber mexiusanum*",
      "*Tuber borchii*",
      "*Tuber gennadii*",
      "*Tuber lyonii*", 
      "*Tuber aestivum*"
    )
  ) +
  scale_y_continuous(limits = c(-20, 11000), expand = c(0, 0))

Fig_Tabund_revised

# D) All Fungi, E) Ectomycorrihzal fungi and F) Bacteria stached bars ----------

# Extract top 
ExtractTop(physeq_fungi_all, "Genus", 32)

# Fungi
test_cont_brule <-
  subset_taxa(physeq_fungi_all, Genus %in% c(
    ExtractTop(physeq_fungi_all, "Genus", 32) 
  )) %>% 
  tax_glom(., taxrank="Genus") %>% 
  psmelt() %>% 
  arrange(Genus) %>% 
  mutate(Genus = recode(Genus, "Gibberella" = "Fusarium" , 
                        "Cylindrocarpon" = "Ilyonectria")) %>% 
  mutate(brule_continent = paste(brule, continent, sep="_")) %>% 
  mutate(brule_continent = as.factor(brule_continent),
         brule_continent = recode(brule_continent,
                                  "outside_Australia"= "AO", 
                                  "inside_Australia" = "AI", 
                                  "inside_Europe" = "EI", 
                                  "outside_Europe" = "EO"), 
         brule_continent = fct_relevel(brule_continent, 
                                       "AI","AO", "EI", "EO")) %>% 
  janitor::clean_names()

head(test_cont_brule)

AUEU_plot_bar_fungi <-
  test_cont_brule %>% 
  filter(!genus == "unclassified") %>% 
  select(abundance, brule_continent, genus) %>% 
  group_by(genus, brule_continent) %>% 
  summarise(abund = sum(abundance)) %>% 
  arrange(abund) %>% 
  filter(genus != "Unlcassified" & genus != "unclassified") %>% 
  mutate(genus = as.factor(genus),
         genus = fct_relevel(
           genus, 
           "Tuber","Fusarium","Mortierella","Cryptococcus","Scleroderma",
           "Inocybe","Hebeloma","Ilyonectria","Tomentella",
           "Metarhizium","Hygrocybe","Tetracladium","Trichophaea","Exophiala", 
           "Humicola", "Geminibasidium", "Archaeorhizomyces", "Tarzetta", 
           "Cladosporium", "Haematonectria")) %>% 
  ggplot(aes(x = brule_continent, y = abund, fill = genus)) + 
  geom_bar(stat = "identity") +
  themePlot() +
  theme(axis.title.x = element_blank(),
        axis.text.x =  element_markdown(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
        legend.text = element_text(face = "italic", size=9)) +
  guides(fill= guide_legend(ncol=1)) +
  scale_fill_manual(values = fungi_colors) +
  scale_x_discrete(labels = cont_brule_labels) +
  scale_y_continuous(limits = c(-10, 250000), expand = c(0, 0)) +
  labs(title="Fungi", y="Abundance")

AUEU_plot_bar_fungi


# Ectomycorrhizal 
test_cont_brule_ecm <-
  subset_taxa(physeq_ecm_all, Genus %in% c(
    ExtractTop(physeq_ecm_all, "Genus", 20))) %>% 
  tax_glom(., taxrank="Genus") %>% 
  psmelt() %>% 
  arrange(Genus) %>% 
  mutate(Genus = recode(Genus, "Descomyces" = "Descolea")) %>% 
  mutate(brule_continent = paste(brule, continent, sep="_")) %>% 
  mutate(brule_continent = as.factor(brule_continent),
         brule_continent = recode(brule_continent,
                                  "outside_Australia"= "AO", 
                                  "inside_Australia" = "AI", 
                                  "inside_Europe" = "EI", 
                                  "outside_Europe" = "EO"), 
         brule_continent = fct_relevel(brule_continent,
                                       "AI","AO", "EI", "EO")) %>% 
  janitor::clean_names()

AUEU_plot_bar_ecm <-
  test_cont_brule_ecm %>% 
  filter(!genus == "unclassified") %>% 
  select(abundance, brule_continent, genus) %>% 
  group_by(genus, brule_continent) %>% 
  summarise(abund = sum(abundance)) %>% 
  arrange(abund) %>% 
  filter(genus != "Unlcassified" & genus != "unclassified") %>% 
  mutate(genus = as.factor(genus),
         genus = fct_relevel(
           genus, 
           "Tuber","Scleroderma","Inocybe","Hebeloma","Tomentella",
           "Trichophaea","Tarzetta","Russula","Descolea","Entoloma", 
           "Geopora","Wilcoxina","Cortinarius","Astraeus","Picoa",
           "Hydnobolites","Melanogaster","Clavulina","Genea","Boletus")) %>% 
  ggplot(aes(x = brule_continent, y = abund, fill = genus)) + 
  geom_bar(stat = "identity") +
  themePlot() +
  theme(axis.title.x = element_blank(),
        axis.text.x =  element_markdown(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
        legend.text = element_text(face = "italic", size=9)) +
  guides(fill= guide_legend(ncol=1)) +
  scale_fill_manual(values = ecm_colors) +
  scale_x_discrete(labels = cont_brule_labels) +
  scale_y_continuous(limits = c(-10, 125000), expand = c(0, 0)) +
  labs(title="Ectomycorrhizal fungi", y="Abundance")

AUEU_plot_bar_ecm

# Bacteria 
test_cont_brule_bact <-
  subset_taxa(physeq_bact_all, Class %in% c(
    ExtractTop(physeq_bact_all, "Class", 20))) %>% 
  tax_glom(., taxrank="Class") %>% 
  psmelt() %>% 
  arrange(Class) %>% 
  mutate(brule_continent = paste(brule, continent, sep="_")) %>% 
  mutate(brule_continent = as.factor(brule_continent),
         brule_continent = recode(brule_continent,
                                  "outside_Australia"= "AO", 
                                  "inside_Australia" = "AI", 
                                  "inside_Europe" = "EI", 
                                  "outside_Europe" = "EO"), 
         brule_continent = fct_relevel(brule_continent,
                                       "AI","AO", "EI", "EO")) %>% 
  janitor::clean_names()

AUEU_plot_bar_bact <-
  test_cont_brule_bact %>% 
  select(abundance, brule_continent, class) %>% 
  group_by(class, brule_continent) %>% 
  summarise(abund = sum(abundance)) %>% 
  arrange(abund) %>% 
  filter(class != "Unlcassified" & class != "unclassified") %>%
  mutate(class = as.factor(class),
         class = fct_relevel(
         class, 
         "Alphaproteobacteria", "Acidobacteria-6", "Thermoleophilia", "Actinobacteria", "Betaproteobacteria",
         "Ellin6529","Deltaproteobacteria", "[Chloracidobacteria]", "Planctomycetia","Thaumarchaeota", 
         "Gammaproteobacteria", "Rubrobacteria", "Acidimicrobiia", "Anaerolineae","[Saprospirae]",
         "Cytophagia","Bacilli","TK10","iii1-8","MB-A2-108")) %>% 
  ggplot(aes(x = brule_continent, y = abund, fill = class)) + 
  geom_bar(stat = "identity") +
  themePlot() +
  theme(axis.title.x = element_blank(),
        axis.text.x =  element_markdown(size = 10, angle = 0, hjust = 0.5, vjust = 0.5),
        legend.text = element_text(face = "italic", size=9)) +
  guides(fill= guide_legend(ncol=1)) +
  scale_fill_manual(values = colors_bar_bact) +
  scale_x_discrete(labels = cont_brule_labels) +
  scale_y_continuous(limits = c(-10, 150000), expand = c(0, 0)) +
  labs(title="Bacteria", y="Abundance")

AUEU_plot_bar_bact

# 1) First Approach - averages -------------------------------------------------
# Extract T. melanosporum abundance for each sample -----------------------------
Tmel_readNo <-
  physeq_fungi_all %>% 
  subset_taxa(Species == "Tuber melanosporum") %>% 
  otu_table() %>% 
  as.matrix() %>% 
  t(.) %>% 
  as.data.frame() %>% 
  dplyr::select(FOTU_4) %>% 
  rownames_to_column("sample_id") 

Tmel_readNo

df_alpha_fungi_all_new_Tmel <-
  inner_join(df_alpha_fungi_all_new, Tmel_readNo, by = "sample_id")

head(df_alpha_fungi_all_new_Tmel)

# checking distribution 
df_alpha_fungi_all_new_Tmel %>% 
  ggplot(aes(x = FOTU_4)) +
  geom_histogram(binwidth = 500) +
  labs(title = "Tuber melanosporum Read abundance") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(FOTU_4)), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)

df_alpha_fungi_all_new_Tmel %>% 
  ggplot(aes(x = log(FOTU_4))) +
  geom_histogram(binwidth = 0.5) +
  labs(title = "Tuber melanosporum Read abundance") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(
    xintercept = mean(log(FOTU_4))), 
    color = "darkred",
    linetype = "dashed", linewidth = 2, show.legend = FALSE)



# Average per site difference 
fig_cont_brule_FOTU4_abund <- 
  df_alpha_fungi_all_new_Tmel %>% 
  group_by(site, brule, continent) %>%
  summarise(mean = mean(FOTU_4), .groups = "drop") %>% 
  mutate(brule_continent = paste(brule, continent, sep="_")) %>% 
  mutate(brule_continent = as.factor(brule_continent),
         brule_continent = recode(brule_continent,
                                  "outside_Australia"= "AO", 
                                  "inside_Australia" = "AI", 
                                  "inside_Europe" = "EI", 
                                  "outside_Europe" = "EO"), 
         brule_continent = fct_relevel(brule_continent, "AI","AO", "EI", "EO")) %>% 
  PloMedCat(dataframe=., 
            X_var="brule_continent", 
            Y_var="mean",
            my_labels=CompSampl(df_alpha_fungi_all_new_Tmel %>% 
                        group_by(site, brule, continent) %>%
                        summarise(mean = mean(FOTU_4), .groups = "drop") %>% 
                        mutate(brule_continent = paste(brule, continent, sep="_")) %>% 
                        mutate(brule_continent = as.factor(brule_continent),
                               brule_continent = recode(brule_continent,
                                                        "outside_Australia"= "AO", 
                                                        "inside_Australia" = "AI", 
                                                        "inside_Europe" = "EI", 
                                                        "outside_Europe" = "EO"), 
                               brule_continent = fct_relevel(brule_continent, "AI","AO", "EI", "EO")),  
                      formula(mean ~ brule_continent)) %>% 
              pull(Letters) %>% as.character(), 
            labels_y=1300) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_markdown(size = 10),
        axis.text.x =  element_markdown(size = 10, angle = 0, hjust = 0.5, vjust = 0.5)) +
  #scale_y_log10() +
  labs(title = "*T. melanosporum*<br>abundance",
       y = "Average abundance per site ",
       x = NULL) +
  scale_x_discrete(labels = cont_brule_labels) +
  scale_color_manual(values = c("#D21E2C","#D21E2C","#D21E2C","#D21E2C"))+ 
  scale_shape_manual(values = c(1,0,2,5))

fig_cont_brule_FOTU4_abund

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

# 2) Second Approach - mixed effect model --------------------------------------
# We cannot use linear mixed effect models since T.mel abundance it is really zero
# inflated... always close to 0 for outside the brule.

# I will use a negative binomial model with a zero inflated component!

library(glmmTMB)

# Step 1 — compare random effects without zero inflation
m_site <- glmmTMB(
  FOTU_4 ~ continent * brule + (1 | site),
  family = nbinom2,
  data = df_alpha_fungi_all_new_Tmel
)

m_tree <- glmmTMB(
  FOTU_4 ~ continent * brule + (1 | site/tree_id),
  family = nbinom2,
  data = df_alpha_fungi_all_new_Tmel
)

anova(m_site, m_tree)
# Adding tree_id does not improve the model, 
# Tree-level variance for T. melanosporum abundance is essentially zero
# All the extra heterogeneity is already captured at the site level

# Final model 
m_site <- glmmTMB(
  FOTU_4 ~ continent * brule + (1 | site),
  family = nbinom2,
  data = df_alpha_fungi_all_new_Tmel
)

# and add/test zero inflation
m_Tmel_zi1 <- glmmTMB(
  FOTU_4 ~ continent * brule + (1 | site),
  ziformula = ~ continent,
  family = nbinom2,
  data = df_alpha_fungi_all_new_Tmel
)

m_Tmel_zi2 <- glmmTMB(
  FOTU_4 ~ continent * brule + (1 | site),
  ziformula = ~ brule,
  family = nbinom2,
  data = df_alpha_fungi_all_new_Tmel
)

m_Tmel_zi3 <- glmmTMB(
  FOTU_4 ~ continent * brule + (1 | site),
  ziformula = ~ continent * brule,
  family = nbinom2,
  data = df_alpha_fungi_all_new_Tmel
)


m_Tmel_zi0 <- glmmTMB(
  FOTU_4 ~ continent * brule + (1 | site),
  ziformula = ~ 1,
  family = nbinom2,
  data = df_alpha_fungi_all_new_Tmel
)

AIC(m_site, m_Tmel_zi0, m_Tmel_zi1, m_Tmel_zi2, m_Tmel_zi3)
anova(m_site, m_Tmel_zi0, m_Tmel_zi1, m_Tmel_zi2, m_Tmel_zi3)

# The zero inflated model us supported for a better prediciton!

# final model
m_Tmel_zi2 <- glmmTMB(
  FOTU_4 ~ continent * brule + (1 | site),
  ziformula = ~ brule,
  family = nbinom2,
  data = df_alpha_fungi_all_new_Tmel
)

summary(m_Tmel_zi2)

# Best model diagnostics (Figure S16) ------------------------------------------

mm_Tmel_negbinom_sim <- simulateResiduals(m_Tmel_zi2)

par(mfrow = c(2, 2))
# KS test for correct distribution of residuals
testUniformity(mm_Tmel_negbinom_sim)
# KS test for correct distribution within and between groups
testCategorical(mm_Tmel_negbinom_sim, df_alpha_fungi_all_new_Tmel$continent_brule)
# Dispersion test - for details see ?testDispersion
testDispersion(mm_Tmel_negbinom_sim) # tests under and overdispersion
# Outlier test (number of observations outside simulation envelope)
testOutliers(mm_Tmel_negbinom_sim, type = "bootstrap")
par(mfrow = c(1, 1)) 

# The issue I can see is the Dispersion test
# However:
# - The deviation is small
# - The fitted value (red line) is well inside the simulated distribution
# - This is not extreme overdispersion
# It is common for ZINB models with:
# > many zeros
# > heterogeneous ecological processes
# > large sample sizes
# > DHARMa is again being conservative, not flagging a catastrophic problem.

plot_grid(
  emm_plot_sig_h_zinb_cond(m_Tmel_zi2, 
                           title ="*T. melanosporum* (conditional abundance)"),
  plot_zi_presence(m_Tmel_zi2, 
                   title = "*T. melanosporum* (presence probability)"),
  ncol = 1, 
  nrow = 2, 
  align = "hv", 
  axis = "r",
  rel_heights = c(0.55, 0.35),
  labels = c("A", "B"))

emmeans(m_Tmel_zi2, ~ continent | brule, component = "cond",type = "response")
emmeans(m_Tmel_zi2, ~ brule, component = "zi", type = "response") 


# Abundance of Tuber melanosporum was analyzed using a zero-inflated negative 
# binomial mixed-effects model to account for strong zero inflation and overdispersion, 
# with site included as a random effect. The probability of absence was substantially
# higher outside the brûlé, independent of continent. Conditional on presence, 
# T. melanosporum abundance did not differ between Europe and Australia inside 
# the brûlé. In contrast, outside the brûlé, conditional abundance was
# significantly higher in Australia than in Europe (Australia/Europe ratio = 3.38,
# p = 0.0005). These results indicate that continental differences in abundance 
# are expressed primarily outside the brûlé, whereas inside the brûlé abundance 
# converges across continents.

Fig_5_bars_abund <-
  ggarrange(
    ggarrange(
      plot_grid(
        plot_grid(
          emm_plot_sig_h_zinb_cond(m_Tmel_zi2, 
                                   title ="*T. melanosporum* (conditional abundance)"),
          plot_zi_presence(m_Tmel_zi2,
                           title = "*T. melanosporum* (presence probability)"),
          ncol = 1, 
          nrow = 2, 
          align = "hv", 
          axis = "r",
          rel_heights = c(0.55, 0.35),
          labels = c("A", "B")),
        Fig_Tabund_revised +
          theme(axis.text.x = element_markdown(size = 9, angle = 45, hjust = 1, vjust = 1.1)) +
          labs(title = "Genus *Tuber* abundance"),
        ncol=2,
        nrow=1,
        #align="hv",
        #axis = "tb",
        rel_widths=c(0.25, 0.75)),
      plot_grid(
          AUEU_plot_bar_fungi,
          AUEU_plot_bar_ecm,
          AUEU_plot_bar_bact,
          ncol=3,
          nrow=1,
          labels = c("D", "E", "F")),
    ncol = 1, 
    nrow=2),
    ggplot_table_all + 
      theme(plot.margin=ggplot2::margin(0,0,0,0)),
    ncol = 2,
    nrow = 1, 
    widths = c(0.9, 0.1)
  )

Fig_5_bars_abund

# Final figure -----------------------------------------------------------------
grid.arrange(
  Fig_5_bars_abund,
  top = text_grob("MICROBIOME COMPOSITION AND ABUNDANCE", size = 12, face = 2))

# Caption. Model-based predictions from a zero-inflated negative binomial mixed 
# model. Panel A shows the predicted probability of presence of Tuber melanosporum,
# which is high inside the brûlé and markedly reduced outside in both continents. 
# Panel B shows predicted abundance conditional on presence, indicating similarly 
# high abundance inside the brûlé in Europe and Australia, but significantly lower
# abundance outside the brûlé in Europe compared to Australia. Error bars
# represent 95% confidence intervals.

# Zero-inflated negative binomial mixed-model predictions for T. melanosporum. 
# (A) Presence probability: points show observed prevalence per orchard (Site), 
# and red points with 95% CIs show model-predicted probability of presence. 
# (B) Conditional abundance: points show observed abundances among non-zero 
# samples, and red points with 95% CIs show the model-predicted mean abundance 
# given presence. Models included continent, brûlé, and their interaction as 
# fixed effects and Site as a random intercept.

# **************************************************----------------------------

# **** FIGURE 6 **** -----------------------------------------------------------

# 1) Log10 read abundance plots per genus --------------------------------------

all_abund_fungi_genus <-
  left_join(
    subset_taxa(physeq_fungi_all, Genus %in% c(
      ExtractTop(physeq_fungi_all, "Genus", 1000))) %>% 
      tax_glom(., taxrank="Genus") %>%
      otu_table() %>% 
      as.matrix() %>% 
      as.data.frame() %>% 
      rownames_to_column("OTU_ID"),
    subset_taxa(physeq_fungi_all, Genus %in% c(
      ExtractTop(physeq_fungi_all, "Genus", 1000))) %>% 
      tax_glom(., taxrank="Genus") %>% 
      tax_table() %>% 
      as.matrix() %>% 
      as.data.frame() %>% 
      rownames_to_column("OTU_ID"),
    by = "OTU_ID") %>% 
  dplyr::select(contains("ITS"), Genus) %>% 
  pivot_longer(-Genus, values_to = "ReadNo", names_to = "SampleID") %>% 
  left_join(., physeq_fungi_all@sam_data %>% 
              as.matrix() %>% 
              as.data.frame() %>% 
              dplyr::select(SampleID, continent, brule, site), 
            by="SampleID") %>% 
  mutate(continent_brule = interaction(continent, brule, sep = "_")) %>% 
  janitor::clean_names()

all_abund_fungi_genus
unique(all_abund_fungi_genus$genus)

# Manually curated guild selection ---------------------------------------------

# Ectomycorrhizal genera -------------------------------------------------------
ecm_genera <- c("Tuber", "Hebeloma", "Wilcoxina", "Scleroderma", "Tomentella", 
                "Descomyces", "Tarzetta", "Inocybe", "Cortinarius", "Amanita", 
                "Tylospora", "Cenococcum", "Laccaria", "Paxillus", "Lactarius", 
                "Russula", "Tricholoma", "Boletus", "Pisolithus", "Rhizopogon", 
                "Genea", "Terfezia", "Hydnobolites", "Melanogaster", "Picoa", 
                "Geopora", "Balsamia", "Ruhlandiella", "Amphinema", "Gyroporus", 
                "Sarcosphaera", "Alnicola", "Clavulina", "Elaphocordyceps", 
                "Hygrocybe", "Suillus", "Hymenogaster", "Entoloma", "Astraeus")

# Arbuscular mycorrhizal genera
am_genera <- c("Diversispora", "Glomus", "Claroideoglomus", "Funneliformis", 
               "Ambispora", "Rhizophagus", "Paraglomus", "Septoglomus", 
               "Archaeospora", "Acaulospora", "Scutellospora")

# Orchid mycorrhizal genera
orchid_mycorrhizal <- c("Ceratobasidium", "Tulasnella", "Sebacina", "Serendipita")

# Ericoid mycorrhizal genera
ericoid_mycorrhizal <- c("Oidiodendron", "Pseudotomentella")

# All mycorrhizal genera combined
all_mycorrhizal <- c(ecm_genera)

# plotting significant differences of the top 15

subset_taxa(physeq_fungi_all, Genus %in% ecm_genera) %>% 
  ExtractTop("Genus", 20)

top_ecm <- c("Tuber","Scleroderma","Inocybe","Hebeloma","Tomentella",
             "Tarzetta","Hymenogaster","Russula","Descomyces","Geopora",
             "Wilcoxina","Cortinarius", "Picoa", "Hydnobolites", 
             "Melanogaster")

ecm_plot <-
  all_abund_fungi_genus %>%
  filter(!genus %in% c("Unlcassified", "unlcassified")) %>%
  group_by(genus, site, brule, continent) %>%
  summarise(read_no = mean(read_no), .groups = "drop") %>%
  mutate(continent_brule = interaction(continent, brule, sep = "_")) %>%
  PlotMultiDiff_sig_upg(
    genera_set = top_ecm,
    filter_sig = FALSE
  ) +
  scale_shape_manual(
    values = c(1, 0, 2, 5),
    labels = c(
      "AI=Australia inside", "AO=Australia outside",
      "EI=Europe inside", "EO=Europe outside"
    )
  ) +
  scale_x_discrete(labels = cont_brule_labels) +
  labs(title = "Top ectomycorrhizal genera")

ecm_plot

# All fungal genera ------------------------------------------------------------
# NOTE - not sure I will add this to the final plot!

# remove the ecto from the genera_set
unique(all_abund_fungi_genus$Genus)

non_ecto_fungi <- 
  all_abund_fungi_genus %>% 
  filter(!genus %in% all_mycorrhizal) %>% 
  filter(!genus %in% c("Unclassified", "unclassified")) %>% 
  pull(genus) %>% 
  unique()

non_ecto_fungi 

fungi_plot <- 
  all_abund_fungi_genus %>% 
  filter(!genus %in% c("Unclassified", "unclassified")) %>% 
  group_by(genus, site, brule, continent) %>%
  summarise(read_no = mean(read_no), .groups = "drop") %>% 
  mutate(continent_brule = paste(continent, brule, sep = "_")) %>% 
  PlotMultiDiff_sig(num_rows = NULL, 
                    top_n = 12,
                    genera_set = non_ecto_fungi) +
  scale_shape_manual(values = c(1,0,2,5),
                     labels = c("AI=Australia inside","AO=Australia outside",
                                "EI=Europe inside","EO=Europe outside")) +
  scale_x_discrete(labels = cont_brule_labels) +
  labs(title = "Non-symbiotic fungal genera")

fungi_plot

# Bacteria genera --------------------------------------------------------------

all_abund_bact_genus <-
  left_join(
    subset_taxa(physeq_bact_all, Genus %in% c(
      ExtractTop(physeq_bact_all, "Genus", 100))) %>% 
      tax_glom(., taxrank="Genus") %>% 
      otu_table() %>% 
      as.matrix() %>% 
      as.data.frame() %>% 
      rownames_to_column("OTU_ID"),
    subset_taxa(physeq_bact_all, Genus %in% c(
      ExtractTop(physeq_bact_all, "Genus", 100))) %>% 
      tax_glom(., taxrank="Genus") %>% 
      tax_table() %>% 
      as.matrix() %>% 
      as.data.frame() %>% 
      rownames_to_column("OTU_ID"),
    by = "OTU_ID") %>% 
  dplyr::select(contains("16S"), Genus) %>% 
  pivot_longer(-Genus, values_to = "ReadNo", names_to = "SampleID") %>% 
  left_join(., physeq_bact_all@sam_data %>% 
              as.matrix() %>% 
              as.data.frame() %>% 
              dplyr::select(SampleID, continent, brule, site), 
            by="SampleID") %>% 
  mutate(continent_brule = interaction(continent, brule, sep = "_")) %>% 
  janitor::clean_names()

head(all_abund_bact_genus)

# Find top abundance bacteria 
physeq_bact_all %>% 
subset_taxa(!Genus %in% c("Unlcassified", "unlcassified")) %>% 
  ExtractTop("Genus", 20)

top_bact <- c("*Ca.* Nitrososphaera", "Rhodoplanes", "Rubrobacter", 
              "Kaistobacter","Sphingomonas", "Mycobacterium", "Balneimonas",
              "Skermanella","Bacillus", "Gemmata", "Arthrobacter", 
              "Pseudonocardia","Nitrospira", "*Ca.* Udaeobacter", 
              "Bradyrhizobium")

bact_plot <-
all_abund_bact_genus %>%
  filter(!genus %in% c("Unlcassified", "unlcassified")) %>%
  mutate(genus = ifelse(genus %in% "DA101",
    paste("*Ca.* Udaeobacter"), paste(genus)
  )) %>%
  mutate(genus = ifelse(genus %in% "Candidatus Nitrososphaera",
    paste("*Ca.* Nitrososphaera"), paste(genus)
  )) %>%
  group_by(genus, site, brule, continent) %>%
  summarise(read_no = mean(read_no), .groups = "drop") %>%
  mutate(continent_brule = paste(continent, brule, sep = "_")) %>%
  PlotMultiDiff_sig_upg(
    genera_set = top_bact,
    filter_sig = FALSE) +
  scale_shape_manual(
    values = c(1, 0, 2, 5),
    labels = c(
      "AI=Australia inside", "AO=Australia outside",
      "EI=Europe inside", "EO=Europe outside"
    )
  ) +
  scale_x_discrete(labels = cont_brule_labels) +
  labs(title = "Top bacteria genera")

bact_plot

# Indicator analysis -----------------------------------------------------------

# Discretizing T. melanosporum abudnance ---------------------------------------
# Covert T.mel abundance into high, mid (low), low (0)
head(Tmel_readNo)

# kmeans for discrete groups ---------------------------------------------------
Tmel_readNo$clusters <-
  kmeans(Tmel_readNo$FOTU_4, 3)$cluster

Tmel_readNo <-
  Tmel_readNo %>%
  mutate(
    abundance = ifelse(FOTU_4 >= 1000, "High", FOTU_4),
    abundance = ifelse(FOTU_4 < 1000 & FOTU_4 != 0, "Low", abundance),
    abundance = ifelse(FOTU_4 == 0, "Zero", abundance)
  ) %>%
  mutate(
    clusters = as.factor(clusters),
    clusters = recode_factor(clusters,
      "1" = "Mid", "2" = "High", "3" = "Low"),
    clusters = fct_relevel(
      clusters,
      "Low", "Mid", "High"
    )
  )


table(Tmel_readNo$abundance)

# Histogram kmeans
discrete_Tmel_abund <-
  ggarrange(
    Tmel_readNo %>% 
      ggplot(aes(x = FOTU_4)) +
      geom_histogram(binwidth = 200, 
                     aes(fill = clusters)) +
      themePlot() +
      labs(x="Sequence No.",
           y="Number of samples",
           title = "Sequence abundance distribution") +
      theme(legend.title = element_blank())+
      scale_fill_manual(values=c("High"="#999999", "Mid" = "#E69F00", "Low" = "#56B4E9")),
    ggplot(Tmel_readNo, 
           aes(y = FOTU_4, 
               fill = clusters, color=clusters)) +
      geom_boxplot(outlier.shape = 1, outlier.size = 0.8, outlier.alpha = 0.8) +
      themePlot() +
      labs(x="Clusters",
           y="Sequence No.",
           title = "Kmeans clusters") +
      theme(legend.title = element_blank()) +
      scale_fill_manual(values=c("High"="#999999", "Mid" = "#E69F00", "Low" = "#56B4E9")) +
      scale_color_manual(values=c("High"="#999999", "Mid" = "#E69F00", "Low" = "#56B4E9")),
    ncol = 2,
    nrow = 1,
    widths = c(0.7, 0.3),
    common.legend = TRUE,
    legend = "bottom")

discrete_Tmel_abund

grid.arrange(discrete_Tmel_abund, top=title4)

# 2) Indicator analysis --------------------------------------------------------
# T. melanosporum abundance indicators 

# Make sure both fungi and bacteria have the same sample names 
identical(physeq_fungi_all@sam_data$SampleID, physeq_bact_all@sam_data$SampleID)

physeq_fungi_all@sam_data$sample_id <- 
  gsub("_ITS", "", physeq_fungi_all@sam_data$SampleID)
head(physeq_fungi_all@sam_data)

physeq_bact_all@sam_data$sample_id <- 
  gsub("_16S", "", physeq_bact_all@sam_data$SampleID)
head(physeq_bact_all@sam_data)

identical(physeq_fungi_all@sam_data$sample_id, physeq_bact_all@sam_data$sample_id)

# Not present in first but present in second object
setdiff(physeq_fungi_all@sam_data$sample_id, physeq_bact_all@sam_data$sample_id) 
setdiff(physeq_bact_all@sam_data$sample_id, physeq_fungi_all@sam_data$sample_id) 

physeq_fungi_all %>% 
  subset_samples(Site %in% "Cognac") %>% 
  sample_data()

physeq_bact_all %>% 
  subset_samples(Site %in% "Cognac") %>% 
  sample_data()

sample_names(subset_samples(physeq_bact_all, sample_id %in% "EU_EH1F")) # bacteria
sample_names(subset_samples(physeq_bact_all, sample_id %in% "EU_EH5F"))

sample_names(subset_samples(physeq_fungi_all, sample_id %in% "EU_EH1F"))
sample_names(subset_samples(physeq_fungi_all, sample_id %in% "EU_EH5F")) # fungi

# Renaming "EU_EH1F" and "EU_EH5F"
physeq_bact_all@sam_data["EU_16S_EH1F", "sample_id"] <- "EU_EH5F"

identical(physeq_fungi_all@sam_data$sample_id, physeq_bact_all@sam_data$sample_id)
sample_order <- match(physeq_bact_all@sam_data$sample_id, physeq_fungi_all@sam_data$sample_id)
sample_order

physeq_fungi_all@sam_data <- physeq_fungi_all@sam_data[sample_order, ]
identical(physeq_fungi_all@sam_data$sample_id, physeq_bact_all@sam_data$sample_id)

identical(physeq_fungi_all@sam_data$SampleID, Tmel_readNo$SampleID)
Tmel_readNo <- Tmel_readNo[sample_order, ]

Tmel_readNo$sample_id <- gsub("_ITS", "", Tmel_readNo$SampleID)

# Tuber melanosporum
identical(physeq_fungi_all@sam_data$sample_id, Tmel_readNo$sample_id)
physeq_fungi_all@sam_data$FOTU_4 <- Tmel_readNo$FOTU_4
physeq_fungi_all@sam_data$clusters <- Tmel_readNo$clusters
physeq_fungi_all@sam_data$Tmel_abund <- Tmel_readNo$abundance

ind_fungi_Tmel_rg <-
  physeq_fungi_all %>% 
  subset_taxa(Species!="Tuber melanosporum") %>% 
  GetIndicators(., "clusters", "r.g")

ind_fungi_Tmel_rg

ind_fungi_Tmel_IndValg <-
  physeq_fungi_all %>% 
  subset_taxa(Species!="Tuber melanosporum") %>% 
  GetIndicators(., "clusters", "IndVal.g")

ind_fungi_Tmel_IndValg

# Brule'
table(physeq_fungi_all@sam_data$brule)
table(physeq_fungi_all@sam_data$site)

ind_fungi_brule_rg <-
  physeq_fungi_all %>% 
  subset_taxa(Species!="Tuber melanosporum") %>% 
  GetIndicators(., "brule", "r")

ind_fungi_brule_rg %>% head()
ind_fungi_brule_rg %>% 
  filter(p.adjust <= 0.05)

# Bacteria
identical(physeq_bact_all@sam_data$sample_id, Tmel_readNo$sample_id)
physeq_bact_all@sam_data$FOTU_4 <- Tmel_readNo$FOTU_4
physeq_bact_all@sam_data$clusters <- Tmel_readNo$clusters
physeq_bact_all@sam_data$Tmel_abund <- Tmel_readNo$abundance

ind_bact_Tmel_rg <-
  physeq_bact_all %>% 
  GetIndicators(., "clusters", "r.g")
ind_bact_Tmel_rg %>%  head()

ind_bact_Tmel_IndValg <-
  physeq_bact_all %>% 
  GetIndicators(., "clusters", "IndVal.g")
ind_bact_Tmel_IndValg

# Brule'
table(physeq_bact_all@sam_data$brule)

ind_bact_brule_rg <-
  physeq_bact_all %>% 
  GetIndicators(., "brule", "r")

ind_bact_brule_rg %>% head()
ind_bact_brule_rg %>% 
  filter(p.adjust <= 0.05)

# the "stat" value in the multipatt() output represents the square root of the 
# Indicator Value, which is a measure of the fidelity and specificity of a species
# to a particular group in the Indicator Species Analysis. A higher IV indicates 
# a stronger association between the species and the group.

# adjusting the p-value for multiple testing is not necessary, according to the 
# original authors of the method, De Caceres and Legendre (2009). The authors 
# recommend focusing on the raw indicator values (A and B) and their combined 
# indicator value (IV), rather than relying solely on the p-values. These indicator 
# values provide a measure of the fidelity and specificity of a species to a 
# particular group, which can be more informative than the p-values alone.

# 3) Random Forest Feature selection -------------------------------------------
# Extracting features with Random forest

# Using  RF classification -----------------------------------------------------
RF_res_fungi <-
  df4Boruta(physeq_fungi_all, "clusters") %>% 
  mutate(clusters = as.factor(clusters)) %>% 
  selectFeat(df = ., formula(clusters ~ .)) %>% 
  unique(x=.)

RF_res_fungi

df_RF_plus_ind_fungi <-
  df4Plotting(physeq_fungi_all, RF_res_fungi, ind_fungi_Tmel_rg)

df_RF_plus_ind_fungi

RF_res_bact <-
  df4Boruta_bact(physeq_bact_all, "clusters") %>% 
  mutate(clusters = as.factor(clusters)) %>% 
  selectFeat(df = ., formula(clusters ~ .)) %>% 
  unique(x=.)

RF_res_bact

# Using RF regression ----------------------------------------------------------
RF_reg_res_fungi <-
  df4Boruta(physeq_fungi_all, "FOTU_4") %>% 
  mutate(Tmel = as.numeric(FOTU_4)) %>% 
  selectFeat(df = ., formula(Tmel ~ .)) %>% 
  unique(x=.)

RF_reg_res_fungi

RF_reg_res_bact <-
  df4Boruta_bact(physeq_bact_all, "FOTU_4") %>% 
  mutate(Tmel = as.numeric(FOTU_4)) %>% 
  selectFeat(df = ., formula(Tmel ~ .)) %>% 
  unique(x=.)

RF_reg_res_bact

# Putting data sets together ---------------------------------------------------
# Fungi ------------------------------------------------------------------------
fungi_n_hi_in <-
  full_join(
    ind_fungi_brule_rg %>% 
      filter(s.inside == 1) %>% 
      #filter(p.value <= 0.01) %>% 
      arrange(desc(stat)) %>% 
      dplyr::select(Family) %>% 
      dplyr::count(Family) %>% 
      filter(Family != "Unlcassified" & Family != "Incertae sedis" & Family != "unclassified") %>% 
      rename("Inside" = n) %>% 
      arrange(desc(Inside)) %>% 
      head(10),
    ind_fungi_Tmel_rg %>% 
      filter(s.High == 1) %>% 
      #filter(p.value <= 0.01) %>% 
      arrange(desc(stat)) %>% 
      dplyr::select(Family) %>% 
      dplyr::count(Family) %>% 
      filter(Family != "Unlcassified" & Family != "Incertae sedis" & Family != "unclassified") %>% 
      rename("High_g" = n) %>% 
      arrange(desc(High_g)) %>% 
      head(10),
    by="Family") %>% 
  full_join(.,
            ind_fungi_Tmel_IndValg %>% 
              filter(s.High == 1) %>% 
              #filter(p.value <= 0.01) %>% 
              arrange(desc(stat)) %>% 
              dplyr::select(Family) %>% 
              dplyr::count(Family) %>% 
              filter(Family != "Unlcassified" & Family != "Incertae sedis" & Family != "unclassified") %>% 
              rename("High_v" = n) %>% 
              arrange(desc(High_v)) %>% 
              head(10),
            by="Family") %>% 
  full_join(
    physeq_fungi_all@tax_table %>% 
      as.matrix() %>% 
      as.data.frame() %>% 
      rownames_to_column("OTU_ID") %>% 
      filter(OTU_ID %in% RF_res_fungi) %>% 
      dplyr::select(Family) %>% 
      dplyr::count(Family) %>% 
      filter(Family != "Unlcassified" & Family != "Incertae sedis" & Family != "unclassified") %>% 
      rename("RF_c" = n) %>% 
      arrange(desc(RF_c)) %>% 
      head(10), 
    by="Family") %>% 
  full_join(
    physeq_fungi_all@tax_table %>% 
      as.matrix() %>% 
      as.data.frame() %>% 
      rownames_to_column("OTU_ID") %>% 
      filter(OTU_ID %in% RF_reg_res_fungi) %>% 
      dplyr::select(Family) %>% 
      dplyr::count(Family) %>% 
      filter(Family != "Unlcassified" & Family != "Incertae sedis" & Family != "unclassified") %>% 
      rename("RF_r" = n) %>% 
      arrange(desc(RF_r)) %>% 
      head(10),
    by="Family") %>% 
  pivot_longer(-Family) 

fungi_n_hi_in

fungi_n_low_out <-
  full_join(
    ind_fungi_brule_rg %>% 
      filter(s.outside == 1) %>% 
      #filter(p.value <= 0.01) %>% 
      arrange(desc(stat)) %>% 
      dplyr::select(Family) %>% 
      dplyr::count(Family) %>% 
      filter(Family != "Unlcassified" & Family != "Incertae sedis" & Family != "unclassified") %>% 
      rename("Outside" = n) %>% 
      arrange(desc(Outside)) %>% 
      head(10),
    ind_fungi_Tmel_rg %>% 
      filter(s.Low == 1) %>% 
      #filter(p.value <= 0.01) %>% 
      arrange(desc(stat)) %>% 
      dplyr::select(Family) %>% 
      dplyr::count(Family) %>% 
      filter(Family != "Unlcassified" & Family != "Incertae sedis" & Family != "unclassified") %>% 
      rename("Low_g" = n) %>% 
      arrange(desc(Low_g)) %>% 
      head(10),
    by="Family") %>% 
  full_join(.,
            ind_fungi_Tmel_IndValg %>% 
              filter(s.Low == 1) %>% 
              #filter(p.value <= 0.01) %>% 
              arrange(desc(stat)) %>% 
              dplyr::select(Family) %>% 
              dplyr::count(Family) %>% 
              filter(Family != "Unlcassified" & Family != "Incertae sedis" & Family != "unclassified") %>% 
              rename("Low_v" = n) %>% 
              arrange(desc(Low_v)) %>% 
              head(10),
            by="Family") %>% 
  pivot_longer(-Family) 


rbind(
  fungi_n_hi_in,
  fungi_n_low_out) %>% 
  as.data.frame() %>% 
  ggplot(aes(x=name, y=Family, size=value)) +
  geom_point() +
  theme_bw() +
  theme(
    plot.title = element_text(size=12),
    axis.text.x = element_text(size=8, angle = 20, hjust = 1, vjust = 1),
    axis.text.y = element_text(size=8),
    legend.key.height = unit(0.6, "cm"), legend.key.width = unit(0.4, "cm"),
    legend.position = "right")


# Bacteria ---------------------------------------------------------------------
bact_n_hi_in <-
  full_join(
    ind_bact_brule_rg %>% 
      filter(s.inside == 1) %>% 
      #filter(p.value <= 0.01) %>% 
      arrange(desc(stat)) %>% 
      dplyr::select(Family) %>% 
      dplyr::count(Family) %>% 
      filter(Family != "Unlcassified" & Family != "Incertae sedis" & Family != "unclassified") %>% 
      rename("Inside" = n) %>% 
      arrange(desc(Inside)) %>% 
      head(10),
    ind_bact_Tmel_rg %>% 
      filter(s.High == 1) %>% 
      #filter(p.value <= 0.01) %>% 
      arrange(desc(stat)) %>% 
      dplyr::select(Family) %>% 
      dplyr::count(Family) %>% 
      filter(Family != "Unlcassified" & Family != "Incertae sedis" & Family != "unclassified") %>% 
      rename("High_g" = n) %>% 
      arrange(desc(High_g)) %>% 
      head(10),
    by="Family") %>% 
  full_join(.,
            ind_bact_Tmel_IndValg %>% 
              filter(s.High == 1) %>% 
              #filter(p.value <= 0.01) %>% 
              arrange(desc(stat)) %>% 
              dplyr::select(Family) %>% 
              dplyr::count(Family) %>% 
              filter(Family != "Unlcassified" & Family != "Incertae sedis" & Family != "unclassified") %>% 
              rename("High_v" = n) %>% 
              arrange(desc(High_v)) %>% 
              head(10),
            by="Family") %>% 
  full_join(
    physeq_bact_all@tax_table %>% 
      as.matrix() %>% 
      as.data.frame() %>% 
      rownames_to_column("OTU_ID") %>% 
      filter(OTU_ID %in% RF_res_bact) %>% 
      dplyr::select(Family) %>% 
      dplyr::count(Family) %>% 
      filter(Family != "Unlcassified" & Family != "Incertae sedis" & Family != "unclassified") %>% 
      rename("RF_c" = n) %>% 
      arrange(desc(RF_c)) %>% 
      head(10),  
    by = "Family"
  ) %>% 
  full_join(
    physeq_bact_all@tax_table %>% 
      as.matrix() %>% 
      as.data.frame() %>% 
      rownames_to_column("OTU_ID") %>% 
      filter(OTU_ID %in% RF_reg_res_bact) %>% 
      dplyr::select(Family) %>% 
      dplyr::count(Family) %>% 
      filter(Family != "Unlcassified" & Family != "Incertae sedis" & Family != "unclassified") %>% 
      rename("RF_r" = n) %>% 
      arrange(desc(RF_r)) %>% 
      head(10),
    by="Family") %>% 
  pivot_longer(-Family)

bact_n_hi_in

bact_n_low_out <-
  full_join(
    ind_bact_brule_rg %>% 
      filter(s.outside == 1) %>% 
      #filter(p.value <= 0.01) %>% 
      arrange(desc(stat)) %>% 
      dplyr::select(Family) %>% 
      dplyr::count(Family) %>% 
      filter(Family != "Unlcassified" & Family != "Incertae sedis" & Family != "unclassified") %>% 
      rename("Outside" = n) %>% 
      arrange(desc(Outside)) %>% 
      head(10),
    ind_bact_Tmel_rg %>% 
      filter(s.Low == 1) %>% 
      #filter(p.value <= 0.01) %>% 
      arrange(desc(stat)) %>% 
      dplyr::select(Family) %>% 
      dplyr::count(Family) %>% 
      filter(Family != "Unlcassified" & Family != "Incertae sedis" & Family != "unclassified") %>% 
      rename("Low_g" = n) %>% 
      arrange(desc(Low_g)) %>% 
      head(10),
    by="Family") %>% 
  full_join(.,
            ind_bact_Tmel_IndValg %>% 
              filter(s.Low == 1) %>% 
              #filter(p.value <= 0.01) %>% 
              arrange(desc(stat)) %>% 
              dplyr::select(Family) %>% 
              dplyr::count(Family) %>% 
              filter(Family != "Unlcassified" & Family != "Incertae sedis" & Family != "unclassified") %>% 
              rename("Low_v" = n) %>% 
              arrange(desc(Low_v)) %>% 
              head(10),
            by="Family") %>% 
  pivot_longer(-Family) 

bact_n_low_out

rbind(
  bact_n_hi_in,
  bact_n_low_out) %>% 
  as.data.frame() %>% 
  ggplot(aes(x=name, y=Family, size=value)) +
  geom_point() +
  theme_bw() +
  theme(
    plot.title = element_text(size=12),
    axis.text.x = element_text(size=8, angle = 20, hjust = 1, vjust = 1),
    axis.text.y = element_text(size=8),
    legend.key.height = unit(0.6, "cm"), legend.key.width = unit(0.4, "cm"),
    legend.position = "right")


# Final figure -----------------------------------------------------------------
plot_isa_and_RF <-
  ggarrange(
    rbind(
      fungi_n_hi_in,
      fungi_n_low_out) %>% 
      as.data.frame() %>% 
      ggplot(aes(x=name, y=Family, size=value)) +
      geom_point(color="grey60") +
      theme_bw() +
      theme(
        plot.title = element_text(size=12, face = "bold", hjust = 0.5, vjust = 0.5),
        axis.text.x = element_text(size=8, angle = 33, hjust = 1, vjust = 1),
        axis.text.y = element_text(size=8),
        legend.key.height = unit(0.6, "cm"), legend.key.width = unit(0.4, "cm"),
        legend.position = "right") +
      scale_size_continuous("Count",
                            breaks = c(1, 10, 30, 60, 90), # Points to label in the legend
                            range = c(1, 9), # Control the size range in the plot
                            labels = c("1","10", "30","60","90")) +
      labs(title = "Fungi top predictors", 
           x=NULL, y=NULL)+
      scale_x_discrete(labels = site_labels),
    rbind(
      bact_n_hi_in,
      bact_n_low_out) %>% 
      as.data.frame() %>% 
      ggplot(aes(x=name, y=Family, size=value)) +
      geom_point(color="grey60") +
      theme_bw() +
      theme(
        plot.title = element_text(size=12, face = "bold", hjust = 0.5, vjust = 0.5),
        axis.text.x = element_text(size=8, angle = 33, hjust = 1, vjust = 1),
        axis.text.y = element_text(size=8),
        legend.key.height = unit(0.6, "cm"), legend.key.width = unit(0.4, "cm"),
        legend.position = "right") +
      scale_size_continuous("Count",
                            breaks = c(2, 10, 20, 30), # Points to label in the legend
                            range = c(1.5, 4), # Control the size range in the plot
                            labels = c("2","10", "20", "30")) +
      labs(title = "Bacteria top predictors", 
           x=NULL, y=NULL)+
      scale_x_discrete(labels = site_labels),
    ncol = 2,
    nrow = 1,
    align = "hv",
    labels = c("C", "D"))


Multiplot_ISA <-
  ggarrange(ecm_plot,
            bact_plot,
            plot_isa_and_RF,
            ncol = 1,
            nrow = 3, 
            align = "hv",
            labels = c("A", "B"),
            heights = c(0.22, 0.22, 0.4),
            common.legend = TRUE, 
            legend = "bottom")

Multiplot_ISA

grid.arrange(
  Multiplot_ISA,
  top = text_grob(
    "TOP ABUNDANT AND PREDICTORS OF TRUFFLE ABUNDANCE TAXA",
    size = 12,
    face = 2
  )
)


ggsave(
  "Multiplot_ISA.pdf",
  plot = grid.arrange(
    Multiplot_ISA,
    top = text_grob(
      "TOP ABUNDANT AND PREDICTORS OF TRUFFLE ABUNDANCE TAXA",
      size = 12,
      face = 2
    )
  ),
  device = "pdf",width = 9.5, height = 12.5
)

# **************************************************----------------------------

# **** FIGURE 7 **** -----------------------------------------------------------

# 1) Truffle PCoA -------------------------------------------------------------- 

# Bray Curtis ------------------------------------------------------------------

pcoa_gleba_bray <-
  phyloseq::ordinate(physeq_gleba_ev,
    method = "PCoA",
    distance = "bray", binary = TRUE
  )

# variance explained for the first two axes
var_expl <- round(100 * pcoa_gleba_bray$values$Relative_eig[1:2], 1)

x_lab <- sprintf("Axis.1 [%.1f%%]", var_expl[1])
y_lab <- sprintf("Axis.2 [%.1f%%]", var_expl[2])

# Bacteria PCOA 
plot_pcoa_gleba_bray <-
  plot_ordination(physeq_gleba_ev, pcoa_gleba_bray, 
                  type="samples", color="Site", shape="management", 
                  title="PCoA") +
  geom_point(size=2) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 3)),
         shape = guide_legend(ncol=1, override.aes = list(color = "black", size=2))) +
  theme(plot.title = element_text(size = 12, hjust = 0.5, vjust = 0.5, face="bold"),
        axis.text  = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.title =element_blank(), 
        legend.text  = element_text(size = 9), 
        legend.key.height = unit(0.2, "cm"), 
        legend.key.width  = unit(0.2, "cm"),
        legend.spacing.y  = unit(-0.1, "cm")) +
  scale_color_manual(values = palette_site_gleba) +
  scale_shape_manual(values = c(1,2)) +
  labs(x = x_lab, y = y_lab)

plot_pcoa_gleba_bray

# Jaccard ----------------------------------------------------------------------
pcoa_gleba_jac <- 
  phyloseq::ordinate(physeq_gleba_ev, method = "PCoA", distance = "jaccard", binary=TRUE)

# variance explained for the first two axes
var_expl <- round(100 * pcoa_gleba_jac$values$Relative_eig[1:2], 1)

x_lab <- sprintf("Axis.1 [%.1f%%]", var_expl[1])
y_lab <- sprintf("Axis.2 [%.1f%%]", var_expl[2])

plot_pcoa_gleba_jac <-
  plot_ordination(physeq_gleba_ev, pcoa_gleba_jac, 
                  type="samples", color="Site", shape="management", 
                  title="PCoA") +
  geom_point(size=2) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 3)),
         shape = guide_legend(ncol=1, override.aes = list(color = "black", size=2))) +
  theme(plot.title = element_text(size = 12, hjust = 0.5, vjust = 0.5, face="bold"),
        axis.text  = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.title =element_blank(), 
        legend.text  = element_text(size = 9), 
        legend.key.height = unit(0.2, "cm"), 
        legend.key.width  = unit(0.2, "cm"),
        legend.spacing.y  = unit(-0.1, "cm")) +
  scale_color_manual(values = palette_site_gleba) +
  scale_shape_manual(values = c(1,2)) +
  labs(x = x_lab, y = y_lab)

plot_pcoa_gleba_jac

# Bray has a very little difference from jaccard.

# 2) Outgroup Soil PCoA -------------------------------------------------------- 
physeq_bact_out_ev

# Match the same sites that are in the gleba dataset
physeq_bact_out_ev_filt <-
  physeq_bact_out_ev %>% 
  subset_samples(Site %in% c("Braidwood","Manjimup","Spoleto","Cuneo",
                             "San Demetrio","Nimes","Grignan","Cognac",
                             "Zuniga","Albentosa","Moncayo","Norcia",
                             "Romans-Sur-Isere","Acedo"))
otu_table(physeq_bact_out_ev_filt) <- 
  otu_table(physeq_bact_out_ev_filt)[which(rowSums(otu_table(physeq_bact_out_ev_filt)) > 0),] 
physeq_bact_out_ev_filt

# Jaccard PCoA
physeq_bact_out_ev_jac <- 
  phyloseq::ordinate(physeq_bact_out_ev_filt, method = "PCoA", distance = "jaccard", binary=TRUE)

var_expl_out <- round(100 * physeq_bact_out_ev_jac$values$Relative_eig[1:2], 1)

x_lab_out <- sprintf("Axis.1 [%.1f%%]", var_expl_out[1])
y_lab_out <- sprintf("Axis.2 [%.1f%%]", var_expl_out[2])

plot_bact_out_ev_jac <-
  plot_ordination(physeq_bact_out_ev_filt, physeq_bact_out_ev_jac, 
                  type="samples", color="Site", shape="management", 
                  title="PCoA") +
  geom_point(size=2) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 3)),
         shape = guide_legend(ncol=1, override.aes = list(color = "black", size=2))) +
  theme(plot.title = element_text(size = 12, hjust = 0.5, vjust = 0.5, face="bold"),
        axis.text  = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.title =element_blank(), 
        legend.text  = element_text(size = 9), 
        legend.key.height = unit(0.2, "cm"), 
        legend.key.width  = unit(0.2, "cm"),
        legend.spacing.y  = unit(-0.1, "cm")) +
  scale_color_manual(values = palette_site_gleba) +
  scale_shape_manual(values = c(1,2)) +
  labs(x = x_lab_out, y = y_lab_out)

plot_bact_out_ev_jac


# 3) Bacteria genera in truffle gleba -------------------------------------------
Fig_bacteria_Abund_gleba <-
  subset_taxa(physeq_gleba_ev, Genus %in% c(
    ExtractTop(physeq_gleba_ev, "Genus", 32) 
  )) %>% 
  tax_glom(., taxrank="Genus") %>% 
  psmelt() %>% 
  arrange(Genus) %>% 
  filter(Genus != "Unlcassified" & Genus != "unclassified",
         Genus != "Halomonas") %>% 
  mutate(Genus = as.factor(Genus),
         Genus = recode_factor(
           Genus, "Candidatus Xiphinematobacter" = "C. Xiphinematobacter"),
         Genus = fct_relevel(
           Genus, 
           "Bradyrhizobium", "Pedobacter", "Pseudomonas", "Polaromonas",
           "Bacillus","Phyllobacterium","Nesterenkonia", "Paucibacter", 
           "Agrobacterium", "Chitinophaga", "Phaeospirillum","Mycobacterium",
           "Sphingomonas", "Sphingopyxis", "Janthinobacterium", "Paenibacillus", 
           "Brochothrix", "Virgibacillus", "Pelosinus", "Dyadobacter",
           "Sphingobacterium","Pseudoxanthomonas", "Lysobacter","Stenotrophomonas", 
           "Flavobacterium")) %>% 
  # filtering sites and samples 
  #filter(!Site %in% c("Yarra", "Needles", "Launceston"),
  #       !SampleID %in% c("EU_16S_Tmel_12","EU_16S_Tmel_24","EU_16S_Tmel_32","EU_16S_Tmel_40")) %>% 
  mutate(
    Site_code = Site,
    Site_code = recode_factor(
      Site_code,
      "Braidwood" = "Br","Manjimup" = "Ma", "Cuneo" = "Cu","San Demetrio" = "SD",
      "Spoleto" = "Sp","Norcia" = "No","Grignan" = "Gr","Cognac" = "Co",
      "Nimes" = "Ni","Romans-Sur-Isere" = "Ro","Albentosa" = "Al","Moncayo" = "Mo",
      "Zuniga" = "Zu", "Acedo" = "Ac")) %>% 
  ggplot(aes(x = SampleID, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", width = 0.9) +
  themePlot() +
  facet_grid(~Site_code, scales = "free", space = "free") +
  theme(plot.title = element_markdown(size = 12, hjust = 0.5, vjust = 0.5),
        axis.text = element_markdown(),
        axis.title.y = element_markdown(size = 10),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_markdown(size = 12),
        legend.text = element_markdown(size = 10, face = "italic"),
        legend.title = element_blank(), 
        legend.key.height = unit(0.4, "cm"), 
        legend.key.width = unit(0.4, "cm"), 
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 5)) +
  scale_fill_manual(values = colors_bar_gleba) +
  labs(title = "TOP ABUNDANT BACTERIAL TAXA IN TRUFFLES",
       y = "Sequence abundance")

Fig_bacteria_Abund_gleba 

# Colored legend for sites -----------------------------------------------------
gleba_legend <- 
  ggplot() +
  geom_richtext(
    aes(
      x = 0.5, y = 0.5,   # center of the panel
      label = paste0(
        "<span style='color:#003399'>Br=Braidwood, Ma=Manjimup</span>, ",
        "<span style='color:#00843D'>Cu=Cuneo, SD=San Demetrio, Sp=Spoleto, No=Norcia, Gr=Grignan</span><br>",
        "<span style='color:#00843D'>Co=Cognac, Ni=Nimes, Ro=Romans-Sur-Isere, Al=Albentosa, Mo=Moncayo, Zu=Zuniga, Ac=Acedo</span>"
      )
    ),
    fill        = NA,
    label.color = NA,
    hjust = 0.5,          # center horizontally
    vjust = 0.5,          # center vertically
    size  = 4
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void()

gleba_legend

# barplots and legend 
Fig_barplot_gleba_bact <-
  ggarrange(
    Fig_bacteria_Abund_gleba,
    gleba_legend,
    ncol = 1,
    nrow = 2,
    heights =  c(0.9, 0.1)
  )

Fig_barplot_gleba_bact

# 4) Bradyrhizobium abundance analysis --------------------------------------------

# with rarefied data and with orchards removed...
Brady_gleba <-
  full_join(
    physeq_gleba_ev@otu_table %>%
      as.data.frame() %>% 
      rownames_to_column("OTU_ID"), 
    physeq_gleba_ev@tax_table %>%
      as.matrix() %>% 
      as.data.frame() %>% 
      rownames_to_column("OTU_ID"),
    by="OTU_ID") %>% 
  filter(Genus == "Bradyrhizobium") %>% 
  dplyr::select(OTU_ID, contains("_Tmel_")) %>% 
  pivot_longer(cols = -OTU_ID, names_to = "SampleID", values_to = "ReadNo") %>% 
  group_by(SampleID) %>% 
  summarise(Count = sum(ReadNo)) %>% 
  as.data.frame() %>% 
  left_join(., 
            physeq_gleba_ev@sam_data %>%
              as.matrix() %>% 
              as.data.frame() %>% 
              dplyr::select(SampleID, Site, continent, management),
            by="SampleID"
  )


Brady_gleba %>% 
  group_by(continent, Site, management) %>% 
  summarise(ReadNo = sum(Count)) 

# Brady in Outgroup trees
Brady_out <-
  full_join(
    physeq_bact_out_ev@otu_table %>%
      as.data.frame() %>% 
      rownames_to_column("OTU_ID"), 
    physeq_bact_out_ev@tax_table %>%
      as.matrix() %>% 
      as.data.frame() %>% 
      rownames_to_column("OTU_ID"),
    by="OTU_ID") %>% 
  filter(Genus == "Bradyrhizobium") %>% 
  dplyr::select(OTU_ID, contains("_16S_")) %>% 
  pivot_longer(cols = -OTU_ID, names_to = "SampleID", values_to = "ReadNo") %>% 
  group_by(SampleID) %>% 
  summarise(Count = sum(ReadNo)) %>% 
  as.data.frame() %>% 
  left_join(., 
            physeq_bact_out_ev@sam_data %>%
              as.matrix() %>% 
              as.data.frame() %>% 
              dplyr::select(SampleID, Site, continent, management),
            by="SampleID"
  )

Brady_out

Brady_out %>% 
  group_by(continent, Site, management) %>% 
  summarise(ReadNo = sum(Count)) 


# Brady in orchard soils
Brady_soil <-
  full_join(
    physeq_bact_all@otu_table %>%
      as.data.frame() %>% 
      rownames_to_column("OTU_ID"), 
    physeq_bact_all@tax_table %>%
      as.matrix() %>% 
      as.data.frame() %>% 
      rownames_to_column("OTU_ID"),
    by="OTU_ID") %>% 
  filter(Genus == "Bradyrhizobium") %>% 
  dplyr::select(OTU_ID, contains("_16S_")) %>% 
  pivot_longer(cols = -OTU_ID, names_to = "SampleID", values_to = "ReadNo") %>% 
  group_by(SampleID) %>% 
  summarise(Count = sum(ReadNo)) %>% 
  as.data.frame() %>% 
  left_join(., 
            physeq_bact_all@sam_data %>%
              as.matrix() %>% 
              as.data.frame() %>% 
              dplyr::select(SampleID, Site, continent, brule, management),
            by="SampleID"
  )

Brady_soil %>% head()

# Identifying Bradyrhizobium OTUs ---------------------------------------------
as.matrix(physeq_gleba_ev@tax_table) %>% 
  as.data.frame() %>% 
  filter(Genus =="Bradyrhizobium")

as.matrix(physeq_bact_out_ev@tax_table) %>% 
  as.data.frame() %>% 
  filter(Genus =="Bradyrhizobium")

as.matrix(physeq_bact_all@tax_table) %>% 
  as.data.frame() %>% 
  filter(Genus =="Bradyrhizobium")


physeq_gleba_ev@otu_table %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("OTU_ID") %>% 
  filter(OTU_ID %in% c("BOTU_41", "BOTU_21404", "BOTU_21998", "BOTU_22017")) %>% 
  mutate(Sum = rowSums(across(where(is.numeric)))) %>% 
  dplyr::select(OTU_ID, Sum)

physeq_bact_all@otu_table %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("OTU_ID") %>% 
  filter(OTU_ID %in% c("BOTU_41", "BOTU_1678")) %>% 
  mutate(Sum = rowSums(across(where(is.numeric)))) %>% 
  dplyr::select(OTU_ID, Sum)

physeq_bact_out_ev@otu_table %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("OTU_ID") %>% 
  filter(OTU_ID %in% c("BOTU_41", "BOTU_1678")) %>% 
  mutate(Sum = rowSums(across(where(is.numeric)))) %>% 
  dplyr::select(OTU_ID, Sum)


Brady_soil
Brady_gleba
Brady_out


df_Brady_mean <-
  full_join(
    full_join(
      Brady_soil %>% 
        group_by(Site, continent, management, brule) %>% 
        summarise(ReadNo = mean(Count)) %>% 
        rename("Brady_So" = ReadNo) %>% 
        dplyr::select(-continent) %>% 
        pivot_wider(id_cols = c(Site,continent,management), 
                    values_from = Brady_So, names_from = brule) %>% 
        rename("Brady_in_brule" = inside, "Brady_out_brule" = outside) %>% 
        mutate(Group = paste(Site, continent, management, sep = "_")) %>% 
        ungroup() %>% 
        dplyr::select(Brady_in_brule, Brady_out_brule, Group),
      Brady_gleba %>% 
        group_by(Site, continent, management) %>% 
        summarise(ReadNo = mean(Count)) %>%
        rename("Brady_Gl" = ReadNo) %>% 
        mutate(Group = paste(Site, continent, management, sep = "_")) %>% 
        ungroup() %>% 
        dplyr::select(Brady_Gl, Group), 
      by="Group"),
    Brady_out %>% 
      group_by(Site, continent,management) %>% 
      summarise(ReadNo = mean(Count)) %>% 
      rename("Brady_Ou" = ReadNo) %>% 
      mutate(Group = paste(Site, continent, management, sep = "_")) %>% 
      ungroup() %>% 
      dplyr::select(Brady_Ou, Group),
    by="Group") %>% 
  filter(!is.na(Brady_Gl)) %>% 
  pivot_longer(cols = -Group, names_to = "Niche", values_to = "Count") %>% 
  separate(Group, c("Site", "Continent", "Management"), sep = "_")  %>% 
  mutate(Cont_Niche = paste(Continent, Niche, sep = "_"))

df_Brady_mean

df_Brady_mean %>% 
  CompSampl(.,  formula(Count ~ Cont_Niche)) %>% 
  pull(Letters)

# plotting
bray_mean_plot <-
  df_Brady_mean %>% 
  ggplot(., aes(Cont_Niche, Count)) +
  geom_jitter(position=position_jitter(0.5), size=1.8, shape=1, color="darkgrey") +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               color="#D21E2C", shape=18, size=0.9) +
  stat_summary(geom = 'text',
               label = c( "ab","a","a","a","b","a","a","a"), 
               fun= max, aes(y = 150000), size=4, color="black") +
  themePlot() + 
  facet_grid(~Continent, space = "free", scales = "free") +
  scale_x_discrete(labels=c("Gleba", "Inside<br>brûlé", "Outside<br>brûlé", "Outgroup",
                            "Gleba", "Inside<br>brûlé", "Outside<br>brûlé", "Outgroup")) +
  theme(strip.text.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 8),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 8, color = "black", 
                                       angle = 0, hjust = 0.5),
        axis.text.y = element_markdown(size = 9, angle = 0, hjust = 0.5),
        legend.position = "none") +
  scale_y_log10(limits = c(0.1, 150000)) +
  labs(y = "*Bradyrhizobium* spp. log<sub>10</sub>(abundance)")

bray_mean_plot


# Final figure ----------------------------------------------------------------- 
Fig_Comb_gleba_bact <-
  ggarrange(
    Fig_barplot_gleba_bact,
    ggarrange(
      bray_mean_plot,
      ggarrange(
        plot_pcoa_gleba_jac + labs(title = "PCoA truffles"),
        plot_bact_out_ev_jac + labs(title = "PCoA outgroup soil"),
        common.legend = TRUE,
        legend = "right",
        labels = c("C", "D")),
      ncol = 2,
      nrow = 1,
      align = "hv",
      labels = "B",
      widths =  c(0.5, 0.65)),
    ncol = 1,
    nrow = 2,
    labels = "A",
    heights =  c(0.65, 0.36)
  )

Fig_Comb_gleba_bact

ggsave(
  "Reivsions_AEM/Fig_7_gleba_bacteria_revised_new.pdf",
  Fig_Comb_gleba_bact,
  device = "pdf"
)


# ****************************************************************--------------
# FIGURE S13 -------------------------------------------------------------------

soil_comparison <-
  soil_data %>% 
  dplyr::select("continent","Brule", "Site", 
                "pH","Olsen.P","K","Ca","Mg","Zn","Mn","Cu","Fe","NO3","NH4") %>% 
  rename(brule = "Brule", site = "Site", P = "Olsen.P") %>% 
  pivot_longer(cols = -c(continent, site, brule), names_to = "parameter", values_to = "value") %>% 
  mutate(group = interaction(continent, brule, sep = "_")) %>% 
  janitor::clean_names()

soil_comparison


# Kruskal-Wallis- No significance!----------------------------------------------
# For each soil parameter, does the metric differ between continents within each
# brulé category (inside vs outside)?
kw_results <- 
  soil_comparison %>%
  filter(!site %in% c("Romans-Sur-Isere", "Norcia", "Acedo")) %>% # remove the natural plots
  group_by(parameter, brule) %>%
  rstatix::kruskal_test(value ~ continent) %>%
  rstatix::adjust_pvalue(method = "BH")

print(kw_results, n=22)

kw_results_cont <- 
  soil_comparison %>%
  filter(!site %in% c("Romans-Sur-Isere", "Norcia", "Acedo")) %>%
  group_by(parameter, brule) %>%
  summarise(test = list(kruskal.test(value ~ continent)), 
            .groups = "drop") %>%
  mutate(tidy_test = map(test, broom::tidy)) %>%
  unnest(tidy_test, names_sep = "_")

print(kw_results_cont, n = 22)


kw_results_cont_brule <- 
  soil_comparison %>%
  filter(!site %in% c("Romans-Sur-Isere", "Norcia", "Acedo")) %>%
  mutate(group = interaction(continent, brule, sep = "_")) %>%
  group_by(parameter) %>%
  summarise(test = list(kruskal.test(value ~ group)), 
            .groups = "drop") %>%
  mutate(tidy_test = map(test, broom::tidy)) %>%
  unnest(tidy_test, names_sep = "_") %>% 
  ungroup() %>%
  mutate(p_adj = p.adjust(tidy_test_p.value , method = "BH"))

print(kw_results_cont_brule, n=11)

# Dunn test - No significance!--------------------------------------------------
dunn_results <- 
  soil_comparison %>%
  filter(!site %in% c("Romans-Sur-Isere", "Norcia", "Acedo")) %>% 
  group_by(parameter) %>%
  dunn_test(value ~ group, p.adjust.method = "BH")

print(dunn_results, n = 66)

# Wilcoxon - No significance!---------------------------------------------------
soil_wilcox <- 
  soil_comparison %>%
  filter(!site %in% c("Romans-Sur-Isere", "Norcia", "Acedo")) %>% 
  group_by(parameter, brule) %>%
  summarize(
    test = list(wilcox.test(value ~ continent, exact = FALSE)),
    tidy = list(broom::tidy(test[[1]])),
    samples = n(),
    .groups = "drop"
  ) %>%
  unnest(tidy) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"),
         p_label = sprintf("%.3g", p_adj)) %>% 
  as.data.frame()

soil_wilcox

# Labels
param_labs <- c(
  "pH"      = "pH",
  "P" = "P (mg L\u207B\u00B9)",
  "K"       = "K (mg kg\u207B\u00B9)",
  "Ca"      = "Ca (mg kg\u207B\u00B9)",
  "Mg"      = "Mg (mg kg\u207B\u00B9)",
  "Zn"      = "Zn (mg kg\u207B\u00B9)",
  "Mn"      = "Mn (mg kg\u207B\u00B9)",
  "Cu"      = "Cu (mg kg\u207B\u00B9)",
  "Fe"      = "Fe (mg kg\u207B\u00B9)",
  "NO3"     = "NO\u2083\u207B (mg kg\u207B\u00B9)",
  "NH4"     = "NH\u2084\u207A (mg kg\u207B\u00B9)"
)


plot_soil_comparison <- 
  soil_comparison %>% 
  filter(!site %in% c("Romans-Sur-Isere", "Norcia", "Acedo")) %>% 
  ggplot(aes(x = brule, y = value, colour = continent, shape = continent)) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.6),
    size = 1.2) +
  stat_summary(
    fun.data = mean_sdl, fun.args = list(mult = 1),
    geom = "pointrange",
    colour = "black", shape = 18, size = 0.5,
    position = position_dodge(width = 0.6),
    aes(group = continent)
  ) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 6, 
             labeller = labeller(parameter = param_labs)) +
  themePlot() +
  theme(strip.text = element_markdown(size = 10),
        axis.title.y = element_markdown(size = 8),
        axis.title.x = element_blank(),
        axis.text.x  = element_markdown(size = 8, color = "black"),
        axis.text.y  = element_markdown(size = 8),
        legend.position = "right") +
  labs(y = "Soil data") +
  scale_color_manual(values = palette_cont) +
  scale_shape_manual(values = c(1, 2)) +
guides(colour = guide_legend(override.aes = list(size = 2, stroke = 1.5)),
       shape  = guide_legend(override.aes = list(size = 2, stroke = 1.5)))

plot_soil_comparison

# Genrate a table for comparisons ----------------------------------------------
tab_df <- 
  kw_results_cont_brule %>%
  select(parameter, tidy_test_statistic, tidy_test_p.value, p_adj) %>%
  rename(stat = tidy_test_statistic, p = tidy_test_p.value, p.adj = p_adj) %>%
  mutate(across(c(stat, p, p.adj), ~ signif(.x, 3)))

tab_grob <- tableGrob(
  tab_df,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 7,  # <- smaller text
    padding = unit(c(1, 1), "mm")
  )
)

# Final figure -----------------------------------------------------------------
plot_soil_comparison + 
  patchwork::wrap_elements(full = tab_grob) +
  patchwork::plot_layout(widths = c(5, 0.7))


# ****************************************************************--------------
# GRAPHICS RESPONDING REVIWERS -------------------------------------------------

# mean pH and pH distribution --------------------------------------------------

soil_comparison %>% 
  filter(parameter == "pH") %>% 
  ggplot(aes(x = value)) +
  geom_histogram(bins = 6) +
  facet_grid(~brule) + 
  labs(title = "pH") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  geom_vline(data = soil_comparison %>% 
               filter(parameter == "pH") %>% 
               group_by(brule) %>% 
               summarise(mean = mean(value)),
             aes(xintercept = mean), color = "darkred",
             linetype = "dashed", linewidth = 1, show.legend = FALSE)

soil_comparison %>% 
  filter(parameter == "pH") %>% 
  group_by(brule) %>% 
  summarise(mean = mean(value))


soil_comparison %>% 
  filter(parameter == "pH") %>% 
  ggplot(aes(x = value)) +
  geom_histogram(bins = 6) +
  facet_grid(~continent + brule) + 
  labs(title = "pH") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(data = soil_comparison %>% 
               filter(parameter == "pH") %>% 
               group_by(continent, brule) %>% 
               summarise(mean = mean(value)),
             aes(xintercept = mean), color = "darkred",
             linetype = "dashed", linewidth = 1, show.legend = FALSE)


# Precise Number of truffle sampled in each site -------------------------------
# samples per site per continent
df_counts <- 
  physeq_gleba_ev@sam_data %>%
  as.data.frame() %>%
  group_by(continent, site) %>%
  summarise(n_truffles = n(), .groups = "drop")

df_counts

# Mean number of truffles per site within each continent -----------------------
df_counts %>%
  group_by(continent) %>%
  summarise(mean_truffles_per_site = mean(n_truffles),
            sd_truffles_per_site   = sd(n_truffles),
            n_sites                = n())


# Sequencing results -----------------------------------------------------------
sort(sample_sums(physeq_ITS))
sort(sample_sums(physeq_16S))

mean_ITS  <- physeq_ITS@otu_table  %>% sample_sums() %>% mean()
mean_16S  <- physeq_16S@otu_table  %>% sample_sums() %>% mean()

ggarrange(
data.frame(depth = physeq_ITS@otu_table %>% sample_sums()) %>% 
ggplot(aes(x = depth)) +
  geom_histogram(binwidth = 1000) +
  labs(title = "ITS") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(xintercept = physeq_ITS@otu_table %>%
                   sample_sums() %>% 
                   mean()), color = "darkred",
             linetype = "dashed", linewidth = 1, show.legend = FALSE) +
  annotate(
    "text",
    x = mean_ITS,
    y = Inf,
    label = paste0("Mean = ", round(mean_ITS, 0)),
    vjust = 1.5,
    hjust = -0.1,
    color = "darkred",
    size = 3.5
  ),

data.frame(depth = physeq_16S@otu_table %>% sample_sums()) %>% 
  ggplot(aes(x = depth)) +
  geom_histogram(binwidth = 1000) +
  labs(title = "16S") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + 
  geom_vline(aes(xintercept = physeq_16S@otu_table %>%
                   sample_sums() %>% 
                   mean()), color = "darkred",
             linetype = "dashed", linewidth = 1, show.legend = FALSE) +
  annotate(
    "text",
    x = mean_16S,
    y = Inf,
    label = paste0("Mean = ", round(mean_16S, 0)),
    vjust = 1.5,
    hjust = -0.1,
    color = "darkred",
    size = 3.5
  ),
ncol = 2)




