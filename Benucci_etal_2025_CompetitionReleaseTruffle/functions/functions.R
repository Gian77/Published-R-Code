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
# **********************************************************************--------

# ***** FUNCTIONS ***** --------------------------------------------------------

# ***** MAIN ******-------------------------------------------------------------
AlphaMetrics <- function(physeq){
  require(vegan)
  require(tidyverse)
  
  #sample_data(physeq)$read_no <- sample_sums(physeq)
  
  otu <- as.data.frame(as.matrix(t(physeq@otu_table)))
  hill <- vegan::renyi(otu, hill = TRUE, scales = c(0, 1, 2))
  
  if (!identical(rownames(hill), sample_names(physeq))) {
    stop("Sample ordering mismatch between physeq and renyi() output!")
  }
  
  sample_data(physeq)$hill_0 <- hill[, "0", drop = TRUE]
  sample_data(physeq)$hill_1 <- hill[, "1", drop = TRUE]
  sample_data(physeq)$hill_2 <- hill[, "2", drop = TRUE]

  return(physeq)
}


AlphaMetrics(physeq_fungi_all) %>% sample_data() %>% head()

# **********************************************************************--------
# **** FIGURE 1 **** -----------------------------------------------------------

# Generate betadiversty heatmaps -----------------------------------------------
generateHeat <- function(physeq){
  
  require(tidyverse)
  
  b <- 
    as.data.frame(physeq@otu_table) %>% 
    t(.) %>% 
    as.data.frame() %>% 
    vegdist(method = "bray") %>% 
    as.matrix() %>% 
    as_tibble(rownames = "A") %>%
    pivot_longer(-A, names_to="B", values_to="distances") %>% 
    mutate(sample_id = A) %>% 
    left_join(., 
              as.data.frame(as.matrix(physeq@sam_data)) %>% 
                dplyr::select(sample_id, site, brule, continent), 
              by="sample_id") %>% 
    mutate(site_A = site, brule_A = brule, Continent_A = continent) %>% 
    mutate(sample_id = B) %>% 
    dplyr::select(-site, -brule, -continent) %>% 
    left_join(., 
              as.data.frame(as.matrix(physeq@sam_data)) %>% 
                dplyr::select(sample_id, site, brule, continent), 
              by="sample_id") %>% 
    mutate(site_B = site, brule_B = brule, Continent_B = continent) %>% 
    dplyr::select(-site, -brule, -sample_id, -continent) %>% 
    mutate(Group_A = paste(site_A, brule_A, sep = "_"),
           Group_B = paste(site_B, brule_B, sep = "_")) 
  
  
  j <-
    as.data.frame(physeq@otu_table) %>% 
    t(.) %>% 
    as.data.frame() %>% 
    vegdist(method = "jaccard") %>% 
    as.matrix() %>% 
    as_tibble(rownames = "A") %>%
    pivot_longer(-A, names_to="B", values_to="distances") %>% 
    mutate(sample_id = A) %>% 
    left_join(., 
              as.data.frame(as.matrix(physeq@sam_data)) %>% 
                dplyr::select(sample_id, site, brule, continent), 
              by="sample_id") %>% 
    mutate(site_A = site, brule_A = brule, Continent_A = continent) %>% 
    mutate(sample_id = B) %>% 
    dplyr::select(-site, -brule, -continent) %>% 
    left_join(., 
              as.data.frame(as.matrix(physeq@sam_data)) %>% 
                dplyr::select(sample_id, site, brule, continent), 
              by="sample_id") %>% 
    mutate(site_B = site, brule_B = brule, Continent_B = continent) %>% 
    dplyr::select(-site, -brule, -sample_id, -continent) %>% 
    mutate(Group_A = paste(site_A, brule_A, sep = "_"),
           Group_B = paste(site_B, brule_B, sep = "_")) 
  
  
  labels <- tibble(
    x=c(18, 6),
    y=c(2, 22),
    label=c("Jaccard","Bray-Curtis")
  )
  
  
  heat_df <-
    inner_join(
      b %>% 
        mutate(
          site_A = as.factor(site_A),
          site_B = as.factor(site_B)) %>% 
        mutate(site_A = factor(site_A, levels = unique(site_A[order(Continent_A)])),
               site_B = factor(site_B, levels = unique(site_B[order(Continent_B)]))
        ) %>%
        arrange(site_A, site_B) %>% 
        dplyr::select(distances, site_A, site_B) %>% 
        group_by(site_A, site_B) %>% 
        summarise(median = median(distances)) %>% 
        ungroup() %>% 
        filter(!(site_A == site_B)),
      j %>% 
        mutate(
          site_A = as.factor(site_A),
          site_B = as.factor(site_B)) %>% 
        mutate(site_A = factor(site_A, levels = unique(site_A[order(Continent_A)])),
               site_B = factor(site_B, levels = unique(site_B[order(Continent_B)]))
        ) %>%
        arrange(site_A, site_B) %>% 
        dplyr::select(distances, site_A, site_B) %>% 
        group_by(site_A, site_B) %>% 
        summarise(median = median(distances)) %>% 
        ungroup() %>% 
        filter(!(site_A == site_B)),
      by=c("site_A", "site_B")) %>% 
    #dplyr::select(site_A, site_B, bray=median.x, jaccard=median.y) %>% 
    #mutate(distances = if_else(as.numeric(as.factor(site_A)) <= 
    #                             as.numeric(as.factor(site_B)), bray, jaccard))
    dplyr::select(site_A, site_B, bray = median.x, jaccard = median.y) %>%
    mutate(
      distances = case_when(
        as.numeric(as.factor(site_A)) <= as.numeric(as.factor(site_B)) ~ bray,
        TRUE ~ jaccard)) %>% 
    mutate(
      Method = case_when(
        as.numeric(as.factor(site_A)) <= as.numeric(as.factor(site_B)) ~ "bray",
        TRUE ~ "jaccard")) %>% 
    mutate(Continent_A = if_else(site_A %in%c(australian_sites), "Australia", "Europe"),
           Continent_B = if_else(site_B %in%c(australian_sites), "Australia", "Europe"),
           Group = paste(Continent_A, Continent_B, Method, sep="_" )) 
  
  assign("heat_object_main", heat_df, envir = .GlobalEnv)
  
  print(range(heat_df$distances))
  
  heat_df$Method %>% table() %>% print()
  heat_df$Group %>% table() %>% print()
  
  #range_value <-
  #(range(heat_df$distances)[2] - range(heat_df$distances)[1])/2 + range(heat_df$distances)[1]
  
  heat_plot <-
    heat_df %>% 
    ggplot(aes(x=site_A, y=site_B, fill=distances)) +
    geom_tile() +
    theme_classic() +
    scale_fill_gradient2(NULL, low = "#075AFF",mid = "#FFFFCC",high = "#FF0000", 
                         midpoint=0.7633176,
                         limits = c(0.4497408, 1)) +
    geom_text(data=labels, 
              aes(x=x, y=y, label=label), 
              inherit.aes=FALSE, 
              size=5) +
    theme(
      plot.title = element_markdown(size =12, face = "bold",hjust = 0.5, vjust = 0.5),
      plot.subtitle = element_markdown(size = 10, face = "bold",hjust = 0.5, vjust = 0.5),
      axis.text.x = element_markdown(size = 8, angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_markdown(size = 8),
      axis.title = element_blank(),
      legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"),
      legend.title = element_blank(), legend.text = element_text(size = 8),
      legend.position = "bottom",
      legend.margin=ggplot2::margin(0,5,0,0),
      legend.box.margin=ggplot2::margin(0,5,0,0)) 
  
  return(heat_plot)
  
}

generateHeat(physeq_fungi_all)
generateHeat(physeq_ecm_all)
generateHeat(physeq_bact_all)

# Function to mulitplot split -------------------------------------------------- 
generateHeatSplit <- function(physeq, Cont){
  
  # Note: Var can be EU or AU
  require(tidyverse)
  
  b <-
    as.data.frame(physeq@otu_table) %>% 
    t(.) %>% 
    as.data.frame() %>% 
    vegdist(method = "bray") %>% 
    as.matrix() %>% 
    as_tibble(rownames = "A") %>%
    pivot_longer(-A, names_to="B", values_to="distances") %>% 
    mutate(sample_id = A) %>% 
    left_join(., 
              as.data.frame(as.matrix(physeq@sam_data)) %>% 
                dplyr::select(sample_id, site, brule, continent), 
              by="sample_id") %>% 
    mutate(site_A = site, brule_A = brule, Continent_A = continent) %>% 
    mutate(sample_id = B) %>% 
    dplyr::select(-site, -brule, -continent) %>% 
    left_join(., 
              as.data.frame(as.matrix(physeq@sam_data)) %>% 
                dplyr::select(sample_id, site, brule, continent), 
              by="sample_id") %>% 
    mutate(site_B = site, brule_B = brule, Continent_B = continent) %>% 
    dplyr::select(-site, -brule, -sample_id, -continent) %>% 
    mutate(Group_A = paste(site_A, brule_A, sep = "_"),
           Group_B = paste(site_B, brule_B, sep = "_")) 
  
  j <-
    as.data.frame(physeq@otu_table) %>% 
    t(.) %>% 
    as.data.frame() %>% 
    vegdist(method = "jaccard") %>% 
    as.matrix() %>% 
    as_tibble(rownames = "A") %>%
    pivot_longer(-A, names_to="B", values_to="distances") %>% 
    mutate(sample_id = A) %>% 
    left_join(., 
              as.data.frame(as.matrix(physeq@sam_data)) %>% 
                dplyr::select(sample_id, site, brule, continent), 
              by="sample_id") %>% 
    mutate(site_A = site, brule_A = brule, Continent_A = continent) %>% 
    mutate(sample_id = B) %>% 
    dplyr::select(-site, -brule, -continent) %>% 
    left_join(., 
              as.data.frame(as.matrix(physeq@sam_data)) %>% 
                dplyr::select(sample_id, site, brule, continent), 
              by="sample_id") %>% 
    mutate(site_B = site, brule_B = brule, Continent_B = continent) %>% 
    dplyr::select(-site, -brule, -sample_id, -continent) %>% 
    mutate(Group_A = paste(site_A, brule_A, sep = "_"),
           Group_B = paste(site_B, brule_B, sep = "_"))
  
  heat_df <-
    inner_join(
      b %>% 
        dplyr::select(distances, Group_A, Group_B) %>% 
        group_by(Group_A, Group_B) %>% 
        summarise(median = median(distances)) %>% 
        ungroup() %>% 
        separate(Group_A, c("site_A", "Brule_A"), sep = "_", remove = FALSE) %>% 
        separate(Group_B, c("site_B", "Brule_B"), sep = "_", remove = FALSE) %>% 
        mutate(Brule_A = as.factor(Brule_A),
               Brule_B = as.factor(Brule_B)) %>% 
        mutate(Group_A = factor(Group_A, levels = unique(Group_A[order(Brule_A)])),
               Group_B = factor(Group_B, levels = unique(Group_B[order(Brule_B)]))) %>% 
        arrange(Group_A, Group_B), #filter(!(Group_A == Group_B)),
      j %>% 
        dplyr::select(distances, Group_A, Group_B) %>% 
        group_by(Group_A, Group_B) %>% 
        summarise(median = median(distances)) %>% 
        ungroup() %>% 
        separate(Group_A, c("site_A", "Brule_A"), sep = "_", remove = FALSE) %>% 
        separate(Group_B, c("site_B", "Brule_B"), sep = "_", remove = FALSE) %>% 
        mutate(Brule_A = as.factor(Brule_A),
               Brule_B = as.factor(Brule_B)) %>% 
        mutate(Group_A = factor(Group_A, levels = unique(Group_A[order(Brule_A)])),
               Group_B = factor(Group_B, levels = unique(Group_B[order(Brule_B)]))) %>% 
        arrange(Group_A, Group_B), #filter(!(Group_A == Group_B)),
      by=c("Group_A", "Group_B")) %>% 
    dplyr::select(Group_A, Group_B, bray=median.x, jaccard=median.y) %>% 
    mutate(
      distances = case_when(
        as.numeric(as.factor(Group_A)) <= as.numeric(as.factor(Group_B)) ~ bray,
        TRUE ~ jaccard)) %>% 
    mutate(
      Method = case_when(
        as.numeric(as.factor(Group_A)) <= as.numeric(as.factor(Group_B)) ~ "bray",
        TRUE ~ "jaccard")) %>%
    separate(Group_A, c("site_A", "Brule_A"), sep = "_", remove = FALSE) %>% 
    separate(Group_B, c("site_B", "Brule_B"), sep = "_", remove = FALSE) %>% 
    mutate(Group = paste(Brule_A, Brule_B, Method, sep="_" )) %>% 
    mutate(Brule_A = as.factor(Brule_A),
           Brule_B = as.factor(Brule_B)) %>% 
    mutate(Group_A = factor(Group_A, levels = unique(Group_A[order(Brule_A)])),
           Group_B = factor(Group_B, levels = unique(Group_B[order(Brule_B)]))) %>% 
    mutate(Group_A = recode_factor(
      Group_A,
      `Yarra_inside` = "Ya-In" ,`Yarra_outside` = "Ya-Out",
      `Wattles_inside`="Wa-In",`Wattles_outside`="Wa-Out",
      `Launceston_inside` = "La-In" ,`Launceston_outside` = "La-Out",
      `Needles_inside`="Ne-In",`Needles_outside`="Ne-Out",
      `Mole_inside`="Mol-In",`Mole_outside`="Mol-Out",
      `Camberra_inside` = "Ca-In" ,`Camberra_outside` = "Ca-Out",
      `Warri_inside`="Wi-In",`Warri_outside`="Wi-Out",
      `Braidwood_inside`="Br-In",`Braidwood_outside`="Br-Out",
      `Manjimup_inside` = "Ma-In" ,`Manjimup_outside` = "Ma-Out",
      `Pemberton_inside`="Pem-In",`Pemberton_outside`="Pem-Out",
      `Jardee_inside`="Ja-In",`Jardee_outside`="Ja-Out",
      `Cuneo_inside` = "Cu-In" ,`Cuneo_outside` = "Cu-Out",
      `San Demetrio_inside`="SD-In",`San Demetrio_outside`="SD-Out",
      `Spoleto_inside`="Sp-In",`Spoleto_outside`="Sp-Out",
      `Norcia_inside`= "No-In",`Norcia_outside`="No-Out",
      `Cognac_inside` = "Co-In" ,`Cognac_outside` = "Co-Out",
      `Grignan_inside`="Gr-In",`Grignan_outside`="Gr-Out",
      `Perpignan_inside`="Per-In",`Perpignan_outside`="Per-Out",
      `Nimes_inside`="Ni-In",`Nimes_outside`="Ni-Out",
      `Romans-Sur-Isere_inside`= "Ro-In",`Romans-Sur-Isere_outside`="Ro-Out",
      `Albentosa_inside` = "Al-In" ,`Albentosa_outside` = "Al-Out",
      `Moncayo_inside`="Mon-In",`Moncayo_outside`="Mon-Out",
      `Zuniga_inside`="Zu-In",`Zuniga_outside`="Zu-Out",
      `Acedo_inside`= "Ac-In",`Acedo_outside`="Ac-Out")) %>%
    mutate(Group_B = recode_factor(
      Group_B,
      `Yarra_inside` = "Ya-In" ,`Yarra_outside` = "Ya-Out",
      `Wattles_inside`="Wa-In",`Wattles_outside`="Wa-Out",
      `Launceston_inside` = "La-In" ,`Launceston_outside` = "La-Out",
      `Needles_inside`="Ne-In",`Needles_outside`="Ne-Out",
      `Mole_inside`="Mol-In",`Mole_outside`="Mol-Out",
      `Camberra_inside` = "Ca-In" ,`Camberra_outside` = "Ca-Out",
      `Warri_inside`="Wi-In",`Warri_outside`="Wi-Out",
      `Braidwood_inside`="Br-In",`Braidwood_outside`="Br-Out",
      `Manjimup_inside` = "Ma-In" ,`Manjimup_outside` = "Ma-Out",
      `Pemberton_inside`="Pem-In",`Pemberton_outside`="Pem-Out",
      `Jardee_inside`="Ja-In",`Jardee_outside`="Ja-Out",
      `Cuneo_inside` = "Cu-In" ,`Cuneo_outside` = "Cu-Out",
      `San Demetrio_inside`="SD-In",`San Demetrio_outside`="SD-Out",
      `Spoleto_inside`="Sp-In",`Spoleto_outside`="Sp-Out",
      `Norcia_inside`= "No-In",`Norcia_outside`="No-Out",
      `Cognac_inside` = "Co-In" ,`Cognac_outside` = "Co-Out",
      `Grignan_inside`="Gr-In",`Grignan_outside`="Gr-Out",
      `Perpignan_inside`="Per-In",`Perpignan_outside`="Per-Out",
      `Nimes_inside`="Ni-In",`Nimes_outside`="Ni-Out",
      `Romans-Sur-Isere_inside`= "Ro-In",`Romans-Sur-Isere_outside`="Ro-Out",
      `Albentosa_inside` = "Al-In" ,`Albentosa_outside` = "Al-Out",
      `Moncayo_inside`="Mon-In",`Moncayo_outside`="Mon-Out",
      `Zuniga_inside`="Zu-In",`Zuniga_outside`="Zu-Out",
      `Acedo_inside`= "Ac-In",`Acedo_outside`="Ac-Out")) %>% 
    mutate(Group_A = factor(Group_A, levels = unique(Group_A[order(Brule_A)])),
           Group_B = factor(Group_B, levels = unique(Group_B[order(Brule_B)]))) %>% 
    filter(!(Group_A == Group_B)) 
  
  print(heat_df$Group_A)
  print(range(heat_df$distances))
  print(table(heat_df$Method))
  print(table(heat_df$Group))
  
  assign("heat_object_split", heat_df, envir = .GlobalEnv)
  
  if (Cont == "Australia"){
    labels_df <- tibble(
      x=c(18, 6),
      y=c(2, 21),
      label=c("Jaccard", "Bray-Curtis"))
    
  } else {
    labels_df <- tibble(
      x=c(22, 6),
      y=c(2, 25),
      label=c("Jaccard", "Bray-Curtis"))
  }
  
  heat_plot <-
    heat_df %>% 
    ggplot(aes(x=Group_A, y=Group_B, fill=distances)) +
    geom_tile() +
    theme_classic() +
    scale_fill_gradient2(NULL, low = "#075AFF",mid = "#FFFFCC",high = "#FF0000", 
                         midpoint=0.6479731,
                         limits = c(0.234778, 1)) +
    geom_text(data=labels_df, 
              aes(x=x, y=y, label=label), 
              inherit.aes=FALSE, 
              size=5) +
    theme(
      plot.title = element_markdown(size =12, face = "bold", hjust = 0.5 , vjust = 0.5),
      plot.subtitle = element_markdown(size = 10, face = "bold", hjust = 0.5, vjust = 0.5),
      axis.text.x = element_markdown(size = 8, angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_markdown(size = 8),
      axis.title = element_blank(),
      legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"),
      legend.title = element_blank(), legend.text = element_text(size = 8),
      legend.position = "bottom",
      legend.margin=ggplot2::margin(0,5,0,0),
      legend.box.margin=ggplot2::margin(0,5,0,0)) 
  
  
  return(heat_plot)
  
}


physeq_bact_all %>% 
  subset_samples(continent %in% c("Europe")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.) %>% 
  generateHeatSplit(., "Europe") 


# **********************************************************************--------
# **** FIGURE S3-S6 **** -------------------------------------------------------

# Compare beta diversity across groups -----------------------------------------
compareBetadiv <- function(physeq){
  
  require(tidyverse)
  
  b <- 
    as.data.frame(physeq@otu_table) %>% 
    t(.) %>% 
    as.data.frame() %>% 
    vegdist(method = "bray") %>% 
    as.matrix() %>% 
    as_tibble(rownames = "A") %>%
    pivot_longer(-A, names_to="B", values_to="distances") %>% 
    mutate(sample_id = A) %>% 
    left_join(., 
              as.data.frame(as.matrix(physeq@sam_data)) %>% 
                dplyr::select(sample_id, site, brule, continent), 
              by="sample_id") %>% 
    mutate(site_A = site, brule_A = brule, Continent_A = continent) %>% 
    mutate(sample_id = B) %>% 
    dplyr::select(-site, -brule, -continent) %>% 
    left_join(., 
              as.data.frame(as.matrix(physeq@sam_data)) %>% 
                dplyr::select(sample_id, site, brule, continent), 
              by="sample_id") %>% 
    mutate(site_B = site, brule_B = brule, Continent_B = continent) %>% 
    dplyr::select(-site, -brule, -sample_id, -continent) %>% 
    mutate(Group_A = paste(site_A, brule_A, sep = "_"),
           Group_B = paste(site_B, brule_B, sep = "_")) 
  
  
  j <-
    as.data.frame(physeq@otu_table) %>% 
    t(.) %>% 
    as.data.frame() %>% 
    vegdist(method = "jaccard") %>% 
    as.matrix() %>% 
    as_tibble(rownames = "A") %>%
    pivot_longer(-A, names_to="B", values_to="distances") %>% 
    mutate(sample_id = A) %>% 
    left_join(., 
              as.data.frame(as.matrix(physeq@sam_data)) %>% 
                dplyr::select(sample_id, site, brule, continent), 
              by="sample_id") %>% 
    mutate(site_A = site, brule_A = brule, Continent_A = continent) %>% 
    mutate(sample_id = B) %>% 
    dplyr::select(-site, -brule, -continent) %>% 
    left_join(., 
              as.data.frame(as.matrix(physeq@sam_data)) %>% 
                dplyr::select(sample_id, site, brule, continent), 
              by="sample_id") %>% 
    mutate(site_B = site, brule_B = brule, Continent_B = continent) %>% 
    dplyr::select(-site, -brule, -sample_id, -continent) %>% 
    mutate(Group_A = paste(site_A, brule_A, sep = "_"),
           Group_B = paste(site_B, brule_B, sep = "_")) 
  
  
  labels <- tibble(
    x=c(18, 6),
    y=c(2, 22),
    label=c("Jaccard","Bray-Curtis")
  )
  
  australian_sites <- c("Yarra","Wattles","Launceston","Needles","Mole",
                        "Camberra","Warri","Braidwood","Manjimup",
                        "Pemberton","Jardee")  
  
  heat_df <-
    inner_join(
      b %>% 
        mutate(
          site_A = as.factor(site_A),
          site_B = as.factor(site_B)) %>% 
        mutate(site_A = factor(site_A, levels = unique(site_A[order(Continent_A)])),
               site_B = factor(site_B, levels = unique(site_B[order(Continent_B)]))
        ) %>%
        arrange(site_A, site_B) %>% 
        dplyr::select(distances, site_A, site_B) %>% 
        group_by(site_A, site_B) %>% 
        summarise(median = median(distances)) %>% 
        ungroup() %>% 
        filter(!(site_A == site_B)),
      j %>% 
        mutate(
          site_A = as.factor(site_A),
          site_B = as.factor(site_B)) %>% 
        mutate(site_A = factor(site_A, levels = unique(site_A[order(Continent_A)])),
               site_B = factor(site_B, levels = unique(site_B[order(Continent_B)]))
        ) %>%
        arrange(site_A, site_B) %>% 
        dplyr::select(distances, site_A, site_B) %>% 
        group_by(site_A, site_B) %>% 
        summarise(median = median(distances)) %>% 
        ungroup() %>% 
        filter(!(site_A == site_B)),
      by=c("site_A", "site_B")) %>% 
    dplyr::select(site_A, site_B, bray = median.x, jaccard = median.y) %>%
    mutate(
      distances = case_when(
        as.numeric(as.factor(site_A)) <= as.numeric(as.factor(site_B)) ~ bray,
        TRUE ~ jaccard)) %>% 
    mutate(
      Method = case_when(
        as.numeric(as.factor(site_A)) <= as.numeric(as.factor(site_B)) ~ "bray",
        TRUE ~ "jaccard")) %>% 
    mutate(Continent_A = if_else(site_A %in%c(australian_sites), "Australia", "Europe"),
           Continent_B = if_else(site_B %in%c(australian_sites), "Australia", "Europe"),
           Group = paste(Continent_A, Continent_B, Method, sep="_" )) 
  
  assign("heat_object_test", heat_df, envir = .GlobalEnv)
  
  
  diff_df <-
    heat_df %>% 
    mutate(Group = recode(Group, 
                          Australia_Australia_bray="AU-AU-Bray", 
                          Australia_Australia_jaccard= "AU-AU-Jaccard",
                          Europe_Europe_bray = "EU-EU-Bray",
                          Europe_Europe_jaccard= "EU-EU-Jaccard",
                          Australia_Europe_bray = "AU-EU-bray",
                          Europe_Australia_jaccard = "AU-EU-jaccard"),
           Group = as.factor(Group),
           Group = fct_relevel(Group, 
                               "EU-EU-Bray", "AU-AU-Bray", "AU-EU-bray",
                               "EU-EU-Jaccard","AU-AU-Jaccard","AU-EU-jaccard")) %>% 
    dplyr::select(Group, distances) %>% 
    CompSampl(df = .,  formula(distances ~ Group)) %>% 
    rownames_to_column("ID") %>% 
    mutate(ID = as.factor(ID),
           ID = fct_relevel(ID, "EU-EU-Bray", "AU-AU-Bray", "AU-EU-bray",
                            "EU-EU-Jaccard","AU-AU-Jaccard","AU-EU-jaccard")) %>% 
    arrange(ID)
  
  range_plot <-
    heat_df %>% 
    mutate(Group = recode(Group, 
                          Australia_Australia_bray="AU-AU-Bray", 
                          Australia_Australia_jaccard= "AU-AU-Jaccard",
                          Europe_Europe_bray = "EU-EU-Bray",
                          Europe_Europe_jaccard= "EU-EU-Jaccard",
                          Australia_Europe_bray = "AU-EU-bray",
                          Europe_Australia_jaccard = "AU-EU-jaccard"),
           Group = as.factor(Group),
           Group = fct_relevel(Group, 
                               "EU-EU-Bray", "AU-AU-Bray", "AU-EU-bray",
                               "EU-EU-Jaccard","AU-AU-Jaccard","AU-EU-jaccard")) %>% 
    ggplot(aes(x = Group, y = distances, color = Group, group=Group)) +
    geom_jitter(aes(color = Group), alpha = 0.8, shape=1, size=1.5) +
    stat_summary(
      geom="pointrange",
      fun.min = function(z) { quantile(z,0.25) },
      fun.max = function(z) { quantile(z,0.75) },
      fun = median, color="black", shape=18, size=1,
      show.legend = FALSE) +
    labs(title = "Median pointrange") +
    stat_summary(geom = "text", label = diff_df$Letters , fun= max, 
                 aes(y = 1.05), size=4, color="black") +
    theme_classic() +
    theme(
      plot.title = element_markdown(size =10, face = "bold",hjust = 0.5, vjust = 0.5),
      plot.subtitle = element_markdown(size = 10, face = "bold",hjust = 0.5, vjust = 0.5),
      axis.text.x = element_markdown(angle = 20, size = 8, hjust = 0.6, vjust = 0.7),
      axis.text.y = element_markdown(size = 8),
      axis.title = element_blank(),
      legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"),
      legend.title = element_blank(), legend.text = element_text(size = 8),
      legend.position = "top",
      legend.margin=ggplot2::margin(0,0,0,0),
      legend.box.margin=ggplot2::margin(0, 0, 0, 0)) +
    scale_color_manual(values = palette_cont_heat) +
    guides(color = guide_legend(override.aes = list(shape = 15, size = 3.5)))
  
  
  histo_plot <-
    heat_df %>% 
    mutate(Group = recode(Group, 
                          Australia_Australia_bray="AU-AU-Bray", 
                          Australia_Australia_jaccard= "AU-AU-Jaccard",
                          Europe_Europe_bray = "EU-EU-Bray",
                          Europe_Europe_jaccard= "EU-EU-Jaccard",
                          Australia_Europe_bray = "AU-EU-bray",
                          Europe_Australia_jaccard = "AU-EU-jaccard"),
           Group = as.factor(Group),
           Group = fct_relevel(Group, 
                               "EU-EU-Bray", "AU-AU-Bray", "AU-EU-bray",
                               "EU-EU-Jaccard","AU-AU-Jaccard","AU-EU-jaccard")) %>% 
    ggplot(aes(x = distances, fill=Group, group=Group)) +
    geom_histogram(binwidth = 0.01, boundary = 0, closed="left", alpha = 0.4) +
    labs(title = "Median histogram") +
    theme_classic() +
    theme(
      plot.title = element_markdown(size =10, face = "bold",hjust = 0.5, vjust = 0.5),
      plot.subtitle = element_markdown(size = 10, face = "bold",hjust = 0.5, vjust = 0.5),
      axis.text.x = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
      axis.text.y = element_markdown(size = 8),
      axis.title = element_blank(),
      legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"),
      legend.title = element_blank(), legend.text = element_text(size = 8),
      legend.position = "right",
      legend.margin=ggplot2::margin(0,0,0,0),
      legend.box.margin=ggplot2::margin(0, 0, 0, 0)) +
    geom_vline(data = heat_df %>% 
                 mutate(Group = recode(Group, 
                                       Australia_Australia_bray="AU-AU-Bray", 
                                       Australia_Australia_jaccard= "AU-AU-Jaccard",
                                       Europe_Europe_bray = "EU-EU-Bray",
                                       Europe_Europe_jaccard= "EU-EU-Jaccard",
                                       Australia_Europe_bray = "AU-EU-bray",
                                       Europe_Australia_jaccard = "AU-EU-jaccard"),
                        Group = as.factor(Group),
                        Group = fct_relevel(Group, 
                                            "EU-EU-Bray", "AU-AU-Bray", "AU-EU-bray",
                                            "EU-EU-Jaccard","AU-AU-Jaccard","AU-EU-jaccard")) %>% 
                 group_by(Group) %>% 
                 summarise(across(c(distances), list(mean = mean, median = median))),
               aes(xintercept = distances_median, color = Group),
               linetype = "dashed", linewidth = 1, show.legend = FALSE) +
    scale_color_manual(values = palette_cont_heat) +
    scale_fill_manual(values = palette_cont_heat)
  
  plot_all <-
    ggarrange(
      range_plot,
      histo_plot,
      ncol = 1, 
      nrow = 2,
      heights = c(1, 0.7),
      common.legend = TRUE,
      legend = "bottom")
  
  return(plot_all)
  
}


compareBetadiv(physeq_fungi_all)




# Compare betadiversity across groups, sites x brule ---------------------------
# Function to plot point range for split data 
compareBetadivSplit <- function(physeq){
  
  require(tidyverse)
  
  b <- 
    as.data.frame(physeq@otu_table) %>% 
    t(.) %>% 
    as.data.frame() %>% 
    vegdist(method = "bray") %>% 
    as.matrix() %>% 
    as_tibble(rownames = "A") %>%
    pivot_longer(-A, names_to="B", values_to="distances") %>% 
    mutate(sample_id = A) %>% 
    left_join(., 
              as.data.frame(as.matrix(physeq@sam_data)) %>% 
                dplyr::select(sample_id, site, brule, continent), 
              by="sample_id") %>% 
    mutate(site_A = site, brule_A = brule, Continent_A = continent) %>% 
    mutate(sample_id = B) %>% 
    dplyr::select(-site, -brule, -continent) %>% 
    left_join(., 
              as.data.frame(as.matrix(physeq@sam_data)) %>% 
                dplyr::select(sample_id, site, brule, continent), 
              by="sample_id") %>% 
    mutate(site_B = site, brule_B = brule, Continent_B = continent) %>% 
    dplyr::select(-site, -brule, -sample_id, -continent) %>% 
    mutate(Group_A = paste(site_A, brule_A, sep = "_"),
           Group_B = paste(site_B, brule_B, sep = "_")) 
  
  
  j <-
    as.data.frame(physeq@otu_table) %>% 
    t(.) %>% 
    as.data.frame() %>% 
    vegdist(method = "jaccard") %>% 
    as.matrix() %>% 
    as_tibble(rownames = "A") %>%
    pivot_longer(-A, names_to="B", values_to="distances") %>% 
    mutate(sample_id = A) %>% 
    left_join(., 
              as.data.frame(as.matrix(physeq@sam_data)) %>% 
                dplyr::select(sample_id, site, brule, continent), 
              by="sample_id") %>% 
    mutate(site_A = site, brule_A = brule, Continent_A = continent) %>% 
    mutate(sample_id = B) %>% 
    dplyr::select(-site, -brule, -continent) %>% 
    left_join(., 
              as.data.frame(as.matrix(physeq@sam_data)) %>% 
                dplyr::select(sample_id, site, brule, continent), 
              by="sample_id") %>% 
    mutate(site_B = site, brule_B = brule, Continent_B = continent) %>% 
    dplyr::select(-site, -brule, -sample_id, -continent) %>% 
    mutate(Group_A = paste(site_A, brule_A, sep = "_"),
           Group_B = paste(site_B, brule_B, sep = "_")) 
  
  
  heat_df <-
    inner_join(
      b %>% 
        dplyr::select(distances, Group_A, Group_B) %>% 
        group_by(Group_A, Group_B) %>% 
        summarise(median = median(distances)) %>% 
        ungroup() %>% 
        separate(Group_A, c("site_A", "Brule_A"), sep = "_", remove = FALSE) %>% 
        separate(Group_B, c("site_B", "Brule_B"), sep = "_", remove = FALSE) %>% 
        mutate(Brule_A = as.factor(Brule_A),
               Brule_B = as.factor(Brule_B)) %>% 
        mutate(Group_A = factor(Group_A, levels = unique(Group_A[order(Brule_A)])),
               Group_B = factor(Group_B, levels = unique(Group_B[order(Brule_B)]))) %>% 
        arrange(Group_A, Group_B),
      #filter(!(Group_A == Group_B)),
      j %>% 
        dplyr::select(distances, Group_A, Group_B) %>% 
        group_by(Group_A, Group_B) %>% 
        summarise(median = median(distances)) %>% 
        ungroup() %>% 
        separate(Group_A, c("site_A", "Brule_A"), sep = "_", remove = FALSE) %>% 
        separate(Group_B, c("site_B", "Brule_B"), sep = "_", remove = FALSE) %>% 
        mutate(Brule_A = as.factor(Brule_A),
               Brule_B = as.factor(Brule_B)) %>% 
        mutate(Group_A = factor(Group_A, levels = unique(Group_A[order(Brule_A)])),
               Group_B = factor(Group_B, levels = unique(Group_B[order(Brule_B)]))) %>% 
        arrange(Group_A, Group_B), 
      #filter(!(Group_A == Group_B)),
      by=c("Group_A", "Group_B")) %>% 
    dplyr::select(Group_A, Group_B, bray=median.x, jaccard=median.y) %>% 
    #mutate(distances = if_else(as.numeric(as.factor(Group_A)) <= 
    #                             as.numeric(as.factor(Group_B)), bray, jaccard)) %>% 
    mutate(
      distances = case_when(
        as.numeric(as.factor(Group_A)) <= as.numeric(as.factor(Group_B)) ~ bray,
        TRUE ~ jaccard)) %>% 
    mutate(
      Method = case_when(
        as.numeric(as.factor(Group_A)) <= as.numeric(as.factor(Group_B)) ~ "bray",
        TRUE ~ "jaccard")) %>% 
    separate(Group_A, c("site_A", "Brule_A"), sep = "_", remove = FALSE) %>% 
    separate(Group_B, c("site_B", "Brule_B"), sep = "_", remove = FALSE) %>% 
    mutate(Brule_A = as.factor(Brule_A),
           Brule_B = as.factor(Brule_B)) %>% 
    mutate(Group = paste(Brule_A, Brule_B, Method, sep="_" )) %>% 
    mutate(Group_A = factor(Group_A, levels = unique(Group_A[order(Brule_A)])),
           Group_B = factor(Group_B, levels = unique(Group_B[order(Brule_B)]))) %>% 
    filter(!(Group_A == Group_B)) 
  
  assign("heat_object_test", heat_df, envir = .GlobalEnv)
  
  print(range(heat_df$distances))
  print(table(heat_df$Method))
  print(table(heat_df$Group))
  
  diff_df <-
    heat_df %>% 
    mutate(Group = recode(Group, 
                          inside_inside_bray="In-In-Bray", 
                          inside_inside_jaccard= "In-In-Jaccard",
                          inside_outside_bray = "In-Out-Bray",
                          outside_inside_jaccard= "In-Out-Jaccard",
                          outside_outside_bray = "Out-Out-bray",
                          outside_outside_jaccard = "Out-Out-jaccard"),
           Group = as.factor(Group),
           Group = fct_relevel(Group, 
                               "In-In-Bray", "Out-Out-bray", "In-Out-Bray",
                               "In-In-Jaccard","Out-Out-jaccard","In-Out-Jaccard")) %>% 
    dplyr::select(Group, distances) %>% 
    CompSampl(df = .,  formula(distances ~ Group)) %>% 
    rownames_to_column("ID") %>% 
    mutate(ID = as.factor(ID),
           ID = fct_relevel(ID,"In-In-Bray", "Out-Out-bray", "In-Out-Bray",
                            "In-In-Jaccard","Out-Out-jaccard","In-Out-Jaccard")) %>% 
    arrange(ID)
  
  
  range_plot <-
    heat_df %>% 
    mutate(Group = recode(Group, 
                          inside_inside_bray="In-In-Bray", 
                          inside_inside_jaccard= "In-In-Jaccard",
                          inside_outside_bray = "In-Out-Bray",
                          outside_inside_jaccard= "In-Out-Jaccard",
                          outside_outside_bray = "Out-Out-bray",
                          outside_outside_jaccard = "Out-Out-jaccard"),
           Group = as.factor(Group),
           Group = fct_relevel(Group, 
                               "In-In-Bray", "Out-Out-bray", "In-Out-Bray",
                               "In-In-Jaccard","Out-Out-jaccard","In-Out-Jaccard")) %>% 
    ggplot(aes(x = Group, y = distances, color = Group, group=Group)) +
    geom_jitter(aes(color = Group), alpha = 0.8, shape=1, size=1.5) +
    stat_summary(
      geom="pointrange",
      fun.min = function(z) { quantile(z,0.25) },
      fun.max = function(z) { quantile(z,0.75) },
      fun = median, color="black", shape=18, size=1,
      show.legend = FALSE) +
    labs(title = "Median pointrange") +
    stat_summary(geom = "text", label = diff_df$Letters , fun= max, 
                 aes(y = 1.05), size=4, color="black") +
    theme_classic() +
    theme(
      plot.title = element_markdown(size =10, face = "bold",hjust = 0.5, vjust = 0.5),
      plot.subtitle = element_markdown(size = 10, face = "bold",hjust = 0.5, vjust = 0.5),
      axis.text.x = element_markdown(angle = 20, size = 8, hjust = 0.5, vjust = 0.6),
      axis.text.y = element_markdown(size = 8),
      axis.title = element_blank(),
      legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"),
      legend.title = element_blank(), legend.text = element_text(size = 8),
      legend.position = "top",
      legend.margin=ggplot2::margin(0,0,0,0),
      legend.box.margin=ggplot2::margin(0, 0, 0, 0)) +
    scale_color_manual(values = palette_cont_heat) +
    guides(color = guide_legend(override.aes = list(shape = 15, size = 3.5)))
  
  
  histo_plot <-
    heat_df %>% 
    mutate(Group = recode(Group, 
                          inside_inside_bray="In-In-Bray", 
                          inside_inside_jaccard= "In-In-Jaccard",
                          inside_outside_bray = "In-Out-Bray",
                          outside_inside_jaccard= "In-Out-Jaccard",
                          outside_outside_bray = "Out-Out-bray",
                          outside_outside_jaccard = "Out-Out-jaccard"),
           Group = as.factor(Group),
           Group = fct_relevel(Group, 
                               "In-In-Bray", "Out-Out-bray", "In-Out-Bray",
                               "In-In-Jaccard","Out-Out-jaccard","In-Out-Jaccard")) %>% 
    ggplot(aes(x = distances, fill=Group, group=Group)) +
    geom_histogram(binwidth = 0.01, boundary = 0, closed="left", alpha = 0.4) +
    labs(title = "Median histogram") +
    theme_classic() +
    theme(
      plot.title = element_markdown(size =10, face = "bold",hjust = 0.5, vjust = 0.5),
      plot.subtitle = element_markdown(size = 10, face = "bold",hjust = 0.5, vjust = 0.5),
      axis.text.x = element_markdown(angle = 0, size = 8, hjust = 1, vjust = 0.5),
      axis.text.y = element_markdown(size = 8),
      axis.title = element_blank(),
      legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"),
      legend.title = element_blank(), legend.text = element_text(size = 8),
      legend.position = "right",
      legend.margin=ggplot2::margin(0,0,0,0),
      legend.box.margin=ggplot2::margin(0, 0, 0, 0)) +
    geom_vline(data = heat_df %>% 
                 mutate(Group = recode(Group, 
                                       inside_inside_bray="In-In-Bray", 
                                       inside_inside_jaccard= "In-In-Jaccard",
                                       inside_outside_bray = "In-Out-Bray",
                                       outside_inside_jaccard= "In-Out-Jaccard",
                                       outside_outside_bray = "Out-Out-bray",
                                       outside_outside_jaccard = "Out-Out-jaccard"),
                        Group = as.factor(Group),
                        Group = fct_relevel(Group, 
                                            "In-In-Bray", "Out-Out-bray", "In-Out-Bray",
                                            "In-In-Jaccard","Out-Out-jaccard","In-Out-Jaccard")) %>% 
                 group_by(Group) %>% 
                 summarise(across(c(distances), list(mean = mean, median = median))),
               aes(xintercept = distances_median, color = Group),
               linetype = "dashed", linewidth = 1, show.legend = FALSE) +
    scale_color_manual(values = palette_cont_heat) +
    scale_fill_manual(values = palette_cont_heat)
  
  plot_all <-
    ggarrange(
      range_plot,
      histo_plot,
      ncol = 1, 
      nrow = 2,
      heights = c(1, 0.7),
      common.legend = TRUE,
      legend = "bottom")
  
  return(plot_all)
  
}

# **********************************************************************--------
# **** FIGURE S6 **** ----------------------------------------------------------

# Plot PCOA ordinations in phyloseq -------------------------------------------- 
physeq_fungi_all %>% 
  subset_samples(continent %in% c("Europe")) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.) %>% 
  compareBetadivSplit(.)



plot_pcoa_analysis <- function(physeq_obj,
                               title = "PCoA Analysis",
                               subset_continent = NULL,
                               method = "jaccard",
                               show_legend = "none") {
  
  if (!is.null(subset_continent)) {
    md <- as.data.frame(sample_data(physeq_obj))
    keep_ids <- rownames(md)[md$continent %in% subset_continent]
    physeq_obj <- prune_samples(keep_ids, physeq_obj)
    
    physeq_obj <- prune_taxa(taxa_sums(physeq_obj) > 0, physeq_obj)
    physeq_obj <- prune_samples(sample_sums(physeq_obj) > 0, physeq_obj)
    
    title <- paste(title, "-", paste(subset_continent, collapse = ", "))
  }
  
  binary_flag <- (method == "jaccard")
  
  pcoa_result <- ordinate(physeq_obj,
                          method = "PCoA",
                          distance = method,
                          binary = binary_flag)
  
  plot_ordination(physeq_obj, pcoa_result,
                  type = "samples",
                  color = "site",
                  shape = "brule",
                  title = title) +
    ggplot2::geom_point(size = 2) +
    themePlot() +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(shape = 15, size = 3.5)),
      shape = ggplot2::guide_legend(ncol = 1, override.aes = list(color = "black", size = 2.5))
    ) +
    ggplot2::scale_color_manual(values = palette_site) +
    ggplot2::scale_shape_manual(values = c(16, 17)) +
    ggplot2::theme(
      legend.position = show_legend,
      plot.title = ggtext::element_markdown(size = 12)
    )
}

plot_pcoa_analysis(physeq_fungi_all,
                   title = "Fungi", 
                   method = "bray",
                   show_legend = "right")

plot_pcoa_analysis(physeq_fungi_all,
                   title = "Fungi", 
                   subset_continent = "Australia",
                   method = "jaccard",
                   show_legend = "right")


# **** FIGURE 3 **** -----------------------------------------------------------

# Automate PERMANOVA for TABLES ------------------------------------------------

Adjadonis2 <- function(formula, metadata, method, stratum, perm){
  require(tidyverse)
  require(vegan)
  df_adonis <-
    adonis2(formula, 
            metadata,
            method = method,
            starta = stratum,
            permutations = perm, 
            parallel = 8)
  df_adj <-
    cbind(df_adonis,
          as.data.frame(p.adjust(df_adonis$`Pr(>F)`, 
                                 method = "BH")) %>% 
            rename("p.adj"=1))
  ad_list <-
    list(formula, df_adonis,df_adj)
  
  return(ad_list)
}


# AIC adonis2 - model comparisons ----------------------------------------------
# Function to calculate AICc for PERMANOVA. Requires input from adonis or adonis2 {vegan}

AICc.PERMANOVA2 <- function(adonis2.model) {
  
  # check to see if object is an adonis2 model...
  
  if (is.na(adonis2.model$SumOfSqs[1]))
    stop("object not output of adonis2 {vegan} ")
  
  # Ok, now extract appropriate terms from the adonis model Calculating AICc
  # using residual sum of squares (RSS or SSE) since I don't think that adonis
  # returns something I can use as a likelihood function... maximum likelihood
  # and MSE estimates are the same when distribution is gaussian See e.g.
  # https://www.jessicayung.com/mse-as-maximum-likelihood/;
  # https://towardsdatascience.com/probability-concepts-explained-maximum-likelihood-estimation-c7b4342fdbb1
  # So using RSS or MSE estimates is fine as long as the residuals are
  # Gaussian https://robjhyndman.com/hyndsight/aic/ If models have different
  # conditional likelihoods then AIC is not valid. However, comparing models
  # with different error distributions is ok (above link).
  
  
  RSS <- adonis2.model$SumOfSqs[ length(adonis2.model$SumOfSqs) - 1 ]
  MSE <- RSS / adonis2.model$Df[ length(adonis2.model$Df) - 1 ]
  
  nn <- adonis2.model$Df[ length(adonis2.model$Df) ] + 1
  
  k <- nn - adonis2.model$Df[ length(adonis2.model$Df) - 1 ]
  
  
  # AIC : 2*k + n*ln(RSS/n)
  # AICc: AIC + [2k(k+1)]/(n-k-1)
  
  # based on 
  # https://en.wikipedia.org/wiki/Akaike_information_criterion;
  # https://www.statisticshowto.datasciencecentral.com/akaikes-information-criterion/ ;
  # https://www.researchgate.net/post/What_is_the_AIC_formula;
  # http://avesbiodiv.mncn.csic.es/estadistica/ejemploaic.pdf;
  # https://medium.com/better-programming/data-science-modeling-how-to-use-linear-regression-with-python-fdf6ca5481be 
  
  # AIC.g is generalized version of AIC = 2k + n [Ln( 2(pi) RSS/n ) + 1]
  # AIC.pi = k + n [Ln( 2(pi) RSS/(n-k) ) +1],
  
  AIC <- 2*k + nn*log(RSS/nn)
  AIC.g <- 2*k + nn * (1 + log( 2 * pi * RSS / nn))
  AIC.MSE <- 2*k + nn * log(MSE)
  AIC.pi <- k + nn*(1 + log( 2*pi*RSS/(nn-k) )   )
  AICc <- AIC + (2*k*(k + 1))/(nn - k - 1)
  AICc.MSE <- AIC.MSE + (2*k*(k + 1))/(nn - k - 1)
  AICc.pi <- AIC.pi + (2*k*(k + 1))/(nn - k - 1)
  
  output <- list("AIC" = AIC, "AICc" = AICc, "AIC.g" = AIC.g, 
                 "AIC.MSE" = AIC.MSE, "AICc.MSE" = AICc.MSE,
                 "AIC.pi" = AIC.pi, "AICc.pi" = AICc.pi, "k" = k, "N" = nn)
  
  return(output)   
  
}


otu_fungi_cult <- 
  physeq_fungi_all %>%
  subset_samples(management %in% c("cultivated")) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>% 
  otu_table() %>%   # Extract OTU table
  as.matrix() %>% 
  t() %>% 
  as.data.frame()

meta_fungi_cult <- 
  physeq_fungi_all %>%
  subset_samples(management %in% c("cultivated")) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>% 
  sample_data() %>%   # Extract metadata table 
  as.matrix() %>% 
  as.data.frame()

# This is my custom function above which I will use
Adjadonis2(
  formula = formula(otu_fungi_cult ~ continent),
  meta_fungi_cult,
  method = "jaccard",
  stratum = brule,
  perm = how(blocks = meta_fungi_cult$brule, nperm = 999, observed = TRUE)
)[[2]] %>%
  AICc.PERMANOVA2()

# This is form the AICcPermanova R package
Adjadonis2(
  formula = formula(otu_fungi_cult ~ continent),
  meta_fungi_cult,
  method = "jaccard",
  stratum = brule,
  perm = how(blocks = meta_fungi_cult$brule, nperm = 999, observed = TRUE)
)[[2]] %>%
  AICc_permanova2()


# Full model adonis2 -----------------------------------------------------------

Adonis2All <- function(physeq, method = "bray", nperm = 999) {
  
  stopifnot(inherits(physeq, "phyloseq"))
  
  # --- metadata (force plain base types) ---
  meta_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
  meta_df <- as.data.frame(meta_df)  # drop any lingering classes
  
  if (!all(c("continent", "brule") %in% names(meta_df))) {
    stop("sample_data must contain columns: continent, brule")
  }
  
  meta_df$continent <- factor(meta_df$continent)
  meta_df$brule     <- factor(meta_df$brule)
  
  # adonis2/terms.formula is happiest with a list or environment
  meta <- as.list(meta_df)
  
  # --- OTU (samples x taxa) ---
  otu <- t(as.matrix(phyloseq::otu_table(physeq)))
  
  # --- distance ---
  if (identical(method, "jaccard")) {
    otu_use <- (otu > 0) * 1
    dist_obj <- vegan::vegdist(otu_use, method = "jaccard")
  } else {
    dist_obj <- vegan::vegdist(otu, method = method)
  }
  
  # --- blocked permutations ---
  perm_brule <- permute::how(blocks = meta_df$brule,     nperm = nperm, observed = TRUE)
  perm_cont  <- permute::how(blocks = meta_df$continent, nperm = nperm, observed = TRUE)
  
  # --- helpers ---
  pad_aic <- function(aic_df, n_rows) {
    if (is.null(aic_df) || nrow(aic_df) == 0) {
      aic_df <- data.frame(matrix(nrow = 1, ncol = 9))
      names(aic_df) <- c("AIC","AICc","AIC.g","AIC.MSE","AICc.MSE","AIC.pi","AICc.pi","k","N")
    }
    aic_df <- aic_df[1, , drop = FALSE]
    
    out <- aic_df[rep(1, n_rows), , drop = FALSE]
    if (n_rows > 1) out[2:n_rows, ] <- NA
    out
  }
  
  run_one <- function(rhs, model, stratum, perm) {
    
    f <- stats::as.formula(paste0("dist_obj ~ ", rhs))
    environment(f) <- environment()  # ensure dist_obj is visible
    
    res <- vegan::adonis2(f, data = meta, permutations = perm)
    
    ad_df <- as.data.frame(res)
    ad_df$Term <- rownames(ad_df)
    rownames(ad_df) <- NULL
    ad_df$Model <- model
    ad_df$Stratum <- stratum
    ad_df$Permutations <- if (is.numeric(perm)) as.character(perm) else "blocked"
    
    aic1 <- tryCatch(as.data.frame(AICc.PERMANOVA2(res)), error = function(e) NULL)
    cbind(ad_df, pad_aic(aic1, nrow(ad_df)))
  }
  
  dplyr::bind_rows(
    run_one("continent",         "continent",        "-",                   nperm),
    run_one("continent",         "continent",        "brule (blocked)",     perm_brule),
    
    run_one("continent * brule", "continent*brule",  "-",                   nperm),
    run_one("continent * brule", "continent*brule",  "brule (blocked)",     perm_brule),
    run_one("brule * continent", "brule*continent",  "brule (blocked)",     perm_brule),
    
    run_one("brule",             "brule",            "-",                   nperm),
    run_one("continent * brule", "continent*brule",  "continent (blocked)", perm_cont)
  )
}


# permanova using bray-curtis on just cultivates sites
physeq_fungi_all %>%
  subset_samples(management %in% c("cultivated")) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  Adonis2All(method = "bray")

# Split model adonis2 ----------------------------------------------------------

# Permanova on split contiennt x brule -----------------------------------------
Adonis2Split <- function(physeq,
                         method = "bray",
                         subset_continent = NULL,
                         nperm = 999) {
  
  stopifnot(inherits(physeq, "phyloseq"))
  stopifnot(is.character(method), length(method) == 1)
  
  # --- metadata (plain base types) ---
  meta_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
  meta_df <- as.data.frame(meta_df)
  
  req <- c("continent", "brule", "site")
  miss <- setdiff(req, names(meta_df))
  if (length(miss) > 0) {
    stop("sample_data must contain columns: ", paste(miss, collapse = ", "))
  }
  
  # --- optional continent subset (no NSE) ---
  if (!is.null(subset_continent)) {
    keep_ids <- rownames(meta_df)[meta_df$continent %in% subset_continent]
    physeq <- phyloseq::prune_samples(keep_ids, physeq)
    physeq <- phyloseq::prune_taxa(phyloseq::taxa_sums(physeq) > 0, physeq)
    physeq <- phyloseq::prune_samples(phyloseq::sample_sums(physeq) > 0, physeq)
    
    meta_df <- as.data.frame(phyloseq::sample_data(physeq), stringsAsFactors = FALSE)
    meta_df <- as.data.frame(meta_df)
  }
  
  # factors
  meta_df$continent <- factor(meta_df$continent)
  meta_df$brule     <- factor(meta_df$brule)
  meta_df$site      <- factor(meta_df$site)
  
  # adonis2/terms.formula is happiest with a list
  meta <- as.list(meta_df)
  
  # --- OTU (samples x taxa) ---
  otu <- t(as.matrix(phyloseq::otu_table(physeq)))
  
  # --- distance ---
  dist_obj <- if (identical(method, "jaccard")) {
    vegan::vegdist((otu > 0) * 1, method = "jaccard")
  } else {
    vegan::vegdist(otu, method = method)
  }
  
  # --- permutation designs (blocked) ---
  perm_brule <- permute::how(blocks = meta_df$brule, nperm = nperm, observed = TRUE)
  perm_site  <- permute::how(blocks = meta_df$site,  nperm = nperm, observed = TRUE)
  
  # --- helpers ---
  pad_aic <- function(aic_df, n_rows) {
    if (is.null(aic_df) || nrow(aic_df) == 0) {
      aic_df <- data.frame(matrix(nrow = 1, ncol = 9))
      names(aic_df) <- c("AIC","AICc","AIC.g","AIC.MSE","AICc.MSE","AIC.pi","AICc.pi","k","N")
    }
    aic_df <- aic_df[1, , drop = FALSE]
    out <- aic_df[rep(1, n_rows), , drop = FALSE]
    if (n_rows > 1) out[2:n_rows, ] <- NA
    out
  }
  
  run_one <- function(rhs, model, stratum, perm, perm_desc) {
    f <- stats::as.formula(paste0("dist_obj ~ ", rhs))
    environment(f) <- environment()  # make dist_obj visible for adonis2
    
    res <- vegan::adonis2(f, data = meta, permutations = perm)
    
    ad_df <- as.data.frame(res)
    ad_df$Term <- rownames(ad_df)
    rownames(ad_df) <- NULL
    
    ad_df$Model <- model
    ad_df$Stratum <- stratum
    ad_df$Permutations <- perm_desc
    
    aic1 <- tryCatch(as.data.frame(AICc.PERMANOVA2(res)), error = function(e) NULL)
    cbind(ad_df, pad_aic(aic1, nrow(ad_df)))
  }
  
  # --- model suite (mirrors your intent, but without Adjadonis2/NSE) ---
  out <- dplyr::bind_rows(
    # Unrestricted interaction tests (two orderings; adonis2 treats them the same)
    run_one("brule * site", "brule*site", "-",           nperm, as.character(nperm)),
    run_one("site * brule", "site*brule", "-",           nperm, as.character(nperm)),
    
    # Main effects with a stratum (blocked permutations)
    run_one("site",         "site",       "brule",       perm_brule, "blocked by brule"),
    run_one("brule",        "brule",      "site",        perm_site,  "blocked by site"),
    
    # Interaction with blocked permutations (two strata)
    run_one("brule * site", "brule*site", "brule",       perm_brule, "blocked by brule"),
    run_one("brule * site", "brule*site", "site",        perm_site,  "blocked by site")
  )
  
  out
}


physeq_fungi_all %>%
  subset_samples(management %in% "cultivated") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  Adonis2Split(method="bray", subset_continent = c("Europe"))


# Calcuate beta-dispersion across gorups ---------------------------------------
BetadispExtr <- function(physeq, method, Var){
  otu <- as.data.frame(otu_table(physeq))
  metadata <- as.matrix(sample_data(physeq))
  metadata <- as.data.frame(metadata)
  disp <-
    betadisper(
      vegan::vegdist(t(otu), method=method),
      metadata[,Var])
  anova_d <-
    anova(disp,
          permutations = how(nperm=999))
  p_adj <-
    round(p.adjust(anova_d$`Pr(>F)`,
                   "BH"), 4)
  dist_var <-
    vegan::permutest(disp,
                     permutations = 999,
                     pairwise = T)
  return(list(dist_var, p_adj, disp))
}

head(physeq_fungi_all@sam_data)

physeq_fungi_all %>%
  subset_samples(management %in% c("cultivated")) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  BetadispExtr("bray", "site")


multipleBetadisper <- function(physeq, method){
  require(vegan)
  
  betadisp_res <-
    data.frame(
      rbind(
        BetadispExtr(physeq, method, "site")[[1]]$tab,
        BetadispExtr(physeq, method, "brule")[[1]]$tab,
        BetadispExtr(physeq, method,"continent")[[1]]$tab),
      Padj = c(
        BetadispExtr(physeq, method,"site")[[2]],
        BetadispExtr(physeq, method, "brule")[[2]],
        BetadispExtr(physeq, method, "continent")[[2]]),
      Variable = c(rep("Site",2),rep("brule",2),rep("continent",2))
    )
  
  return(betadisp_res)
}


physeq_fungi_all %>%
  subset_samples(management %in% c("cultivated")) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  multipleBetadisper(method = "bray") 

# Calculating betadisper on splitted dataset by continent x brule --------------
spliteBetadisper <- function(physeq, method, Var){
  require(phyloseq)
  require(vegan)
  
  if (Var == "Australia"){
    physeq_sub <- 
      subset_samples(physeq, continent%in%c("Australia"))
    otu_table(physeq_sub) <- 
      otu_table(physeq_sub)[which(rowSums(otu_table(physeq_sub)) > 0),] 
  } else {
    physeq_sub <- 
      subset_samples(physeq, continent%in%c("Europe"))
    otu_table(physeq_sub) <- 
      otu_table(physeq_sub)[which(rowSums(otu_table(physeq_sub)) > 0),] 
  }
  
  betadisp_res <-
    data.frame(
      rbind(
        BetadispExtr(physeq_sub, method, "Site")[[1]]$tab,
        BetadispExtr(physeq_sub, method, "brule")[[1]]$tab),
      Padj = c(
        BetadispExtr(physeq_sub, method, "Site")[[2]],
        BetadispExtr(physeq_sub, method, "brule")[[2]]),
      Variable = c(rep("Site",2), rep("Brule",2)))
  
  return(betadisp_res)
  
}

physeq_fungi_all %>%
  subset_samples(management %in% c("cultivated")) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  spliteBetadisper("bray", "Australia")

# **** Figure S7 **** ----------------------------------------------------------

# Pairways adonis --------------------------------------------------------------
Pairwise_adonis2 <- function(physeq,
                             Var, 
                             dist = "bray", 
                             adj = "BH", 
                             perm = 999) {
  require(vegan)
  require(tidyverse)
  
  sp_matrix <- as.data.frame(t(as.matrix(otu_table(physeq))))
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
    contrast_matrix <- sp_matrix[sp_subset,]
    
    ## fit contrast using adonis
    fit <- vegan::adonis2(
      contrast_matrix ~ group_var[sp_subset],
      method = dist, 
      perm = perm,
      parallel = 8)
    
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

Pairwise_adonis2(physeq_fungi_all, "site")

# Generate pairwise comparisons from phyloseq object
getPairCompare <- function(physeq, distance){
  require(phyloseq)
  require(tidyverse)
  
  data_pairs <-
    Pairwise_adonis2(physeq, "site", dist = distance)$contrast %>% 
    dplyr::select(group1, group2, R2, p_value) %>% 
    mutate(group1 = as.factor(group1),
           group2 = as.factor(group2)) %>% 
    as_tibble() %>% 
    mutate(r2 = if_else(p_value>=0.05, NA, R2)) %>% 
    mutate(group1 = fct_relevel(group1, 
                                "Yarra","Launceston","Needles","Mole","Manjimup",
                                "Pemberton","Jardee","Camberra","Warri","Braidwood",
                                "Wattles","Perpignan","Norcia","Spoleto","Cuneo",
                                "San Demetrio","Nimes","Grignan","Romans-Sur-Isere","Cognac",
                                "Zuniga","Acedo","Albentosa"),
           group2 = fct_relevel(group2,
                                "Launceston","Needles","Mole","Manjimup","Pemberton",
                                "Jardee","Camberra","Warri","Braidwood","Wattles",
                                "Perpignan","Norcia","Spoleto","Cuneo","San Demetrio",
                                "Nimes","Grignan","Romans-Sur-Isere","Cognac","Zuniga",
                                "Acedo","Albentosa","Moncayo"))
  return(data_pairs)
}


adonPairs_physeq_ecm <- 
  getPairCompare(physeq_ecm_all, distance = "bray" )

# Diagonal corrplot ------------------------------------------------------------
PairPermCorrplot <- function(dataframe){
  require(tidyverse)
  
  heat_plot <- dataframe %>% 
    ggplot(aes(x=group1, y=group2, fill=r2)) +
    geom_tile() +
    theme_classic() +
    scale_x_discrete(position="top") +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5), 
      axis.text.x = element_text(size= 8, angle = 90, hjust = 0),
      axis.text.x.top = element_text(vjust = 0.5),
      axis.text.y = element_text(size = 8),
      axis.title = element_blank(),
      plot.background = element_blank(),
      axis.ticks = element_blank(),
      axis.line.x = element_blank(), 
      axis.line.y = element_blank(),
      legend.position = c(0.8,0.3)) 
  
  return(heat_plot)
}

adonPairs_physeq_ecm %>% 
  PairPermCorrplot(.) +
  scale_fill_gradient(expression(R^2), low="red", high="yellow", na.value = "white") +
  labs(title = "Fungi") +
  scale_x_discrete(position="top")


# Analysis of Deviance for Multivariate Generalized Linear Model ---------------
# Fits for Abundance Data

# Extract mvabund dataframe ----------------------------------------------------
DFmvabund <- function(physeq, Var){
  require(mvabund)
  require(tidyverse)
  require(phyloseq)
  
  if (is.null(Var)){
    
    metadata <- 
      physeq@sam_data %>% 
      as.matrix() %>% 
      as.data.frame()
    
    mvabund_DF <-
      physeq@otu_table %>% 
      t() %>% 
      as.data.frame() %>% 
      mvabund()
    
  } else {
    
    sub_set <- subset(sample_data(physeq), continent %in% Var)
    physeq_filt <- merge_phyloseq(otu_table(physeq),
                                  tax_table(physeq),
                                  refseq(physeq),
                                  sub_set)
    otu_table(physeq_filt) <-
      otu_table(physeq_filt)[which(rowSums(otu_table(physeq_filt)) > 0),]
    
    metadata <- 
      physeq_filt@sam_data %>% 
      as.matrix() %>% 
      as.data.frame()
    
    mvabund_DF <-
      physeq_filt@otu_table %>% 
      t() %>% 
      as.data.frame() %>% 
      mvabund()
  }
  return(list(mvabund_DF, metadata))
}

# NOE. Must run this on an HPC cluster or will take forever! No way to 
# parallelize it

# **** FIGURE S8-11 **** -------------------------------------------------------

PlotSiteRich_AU <- function(dataframe, formula, Var, letters_y, limit_y, title){
  require(tidyverse)
  
  PlotRich <- function(dataframe, X_var, Y_var, my_labels, labels_y){
    rich_plot <-
      ggplot(dataframe, aes(x = get(X_var), y = get(Y_var))) +
      geom_jitter(position=position_jitter(0.4), size=1,
                  aes(color=get(X_var),fill=get(X_var), shape = get(X_var))) +
      stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
                   geom="pointrange", color="black", shape=18, size=0.8) +
      stat_summary(geom = 'text', angle=0, label = my_labels, fun= max, 
                   aes(y = labels_y), size=3, color="black") +
      themePlot() +
      expand_limits(y = 0) +
      theme(axis.text.y = element_markdown(angle = 90, hjust = 0.5),
            axis.title.x= element_blank(),
            legend.position = "none")
    return(rich_plot)
  }
  
  Vic_fungi_alpha <-
    dataframe %>% 
    subset(state %in% c("Victoria")) %>%
    PlotRich(.,"Site_Brule",Var,
             CompSampl(., formula) %>% 
               pull(Letters) %>% as.character(), letters_y) +
    labs(title = "Victoria", y = title) +
    theme(
      axis.title.y = element_markdown(size = 10, face = "bold"),
      axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_continuous(n.breaks = 6, limits = c(0, limit_y)) +
    scale_color_manual(values = c("grey60","grey60","grey60","grey60")) 
  
  print(Vic_fungi_alpha)
  
  Tas_fungi_alpha <-
    dataframe %>% 
    subset(state %in% c("Tasmania")) %>%
    PlotRich(.,"Site_Brule",Var,
             CompSampl(., formula) %>% 
               pull(Letters) %>% as.character(), letters_y) +
    labs(title = "Tasmania", y = NULL) +
    theme(
      axis.title.y = element_markdown(size = 10, face = "bold"),
      axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_continuous(n.breaks = 6,limits = c(0, limit_y)) +
    scale_color_manual(values = c("grey60","grey60","grey60","grey60","grey60","grey60"))
  
  print(Tas_fungi_alpha)
  
  NSW_fungi_alpha <-
    dataframe %>% 
    subset(state %in% c("New South Wales")) %>%
    PlotRich(.,"Site_Brule",Var,
             CompSampl(., formula) %>% 
               pull(Letters) %>% as.character(), letters_y) +
    labs(title = "New South Wales", y = NULL) +
    theme(
      axis.title.y = element_markdown(size = 10, face = "bold"),
      axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_continuous(n.breaks = 6, limits = c(0, limit_y)) +
    scale_color_manual(values = c("grey60","grey60","grey60","grey60"))
  
  print(NSW_fungi_alpha)
  
  WA_fungi_alpha <-
    dataframe %>% 
    subset(state %in% c("Western Australia")) %>%
    PlotRich(.,"Site_Brule",Var,
             CompSampl(., formula) %>% 
               pull(Letters) %>% as.character(), letters_y) +
    labs(title = "Western Australia", y = NULL) +
    theme(
      axis.title.y = element_markdown(size = 10, face = "bold"),
      axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_continuous(n.breaks = 6, limits = c(0, limit_y)) +
    scale_color_manual(values = c("grey60","grey60","grey60","grey60","grey60","grey60"))
  
  print(WA_fungi_alpha)
  
  rich_inout <-
    ggarrange(Vic_fungi_alpha,
              Tas_fungi_alpha,
              NSW_fungi_alpha,
              WA_fungi_alpha,
              nrow = 1,
              align = "h", 
              widths = c(0.20, 0.25, 0.25, 0.25))
  
  rich_inout
  
  return(rich_inout)
}


PlotSiteRich_AU(
  df_alpha_fungi_all_new,
  formula(hill_0 ~ site_brule),
  "Richness (hill_0)",
  350,
  355,
  "Fungi"
)


PlotSiteRich_EU <- function(dataframe, formula, Var, letters_y, limit_y, title){
  require(tidyverse)
  
  PlotRich <- function(dataframe, X_var, Y_var, my_labels, labels_y){
    rich_plot <-
      ggplot(dataframe, aes(x = get(X_var), y = get(Y_var))) +
      geom_jitter(position=position_jitter(0.4), size=1,
                  aes(color=get(X_var),fill=get(X_var), shape = get(X_var))) +
      stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
                   geom="pointrange", color="black", shape=18, size=0.8) +
      stat_summary(geom = 'text', angle=90, label = my_labels, fun= max, 
                   aes(y = labels_y), size=3, color="black") +
      themePlot() +
      expand_limits(y = 0) +
      theme(axis.text.y = element_markdown(angle = 90, hjust = 0.5),
            axis.title.x= element_blank(),
            legend.position = "none")
    return(rich_plot)
  }
  
  Italy_fungi_alpha <-
    dataframe %>%
    subset(state %in% c("Italy")) %>%
    PlotRich(.,"Site_Brule",Var,
             CompSampl(.,  formula) %>%
               pull(Letters) %>% as.character(), letters_y) +
    labs(title = "Italy", y = title) +
    theme(
      axis.title.y = element_markdown(size = 10, face = "bold"),
      axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_continuous(n.breaks = 6, limits = c(0, limit_y)) +
    scale_color_manual(values = c("grey60","grey60","grey60","grey60", 
                                  "grey60","grey60","grey60","grey60")) +
    scale_shape_manual(values = c(16, 17, 16, 3, 12, 8, 6, 9))
  
  print(Italy_fungi_alpha)
  
  France_fungi_alpha <-
    dataframe %>%
    subset(state %in% c("France")) %>%
    PlotRich(.,"Site_Brule",Var,
             CompSampl(.,  formula) %>%
               pull(Letters) %>% as.character(), letters_y) +
    labs(title = "France", y = NULL) +
    theme(
      axis.title.y = element_markdown(size = 10, face = "bold"),
      axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_continuous(n.breaks = 6, limits = c(0, limit_y)) +
    scale_color_manual(values = c("grey60","grey60","grey60","grey60","grey60","grey60", 
                                  "grey60","grey60","grey60","grey60")) +
    scale_shape_manual(values = c(16, 17, 16, 3, 12, 8, 6, 9, 10, 4))
  
  Spain_fungi_alpha <-
    dataframe %>%
    subset(state %in% c("Spain")) %>%
    PlotRich(.,"Site_Brule",Var,
             CompSampl(.,  formula) %>%
               pull(Letters) %>% as.character(), letters_y) +
    labs(title = "Spain", y = NULL) +
    scale_y_continuous(n.breaks = 6, limits = c(0, limit_y)) +
    theme(
      axis.title.y = element_markdown(size = 10, face = "bold"),
      axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_color_manual(values = c("grey60","grey60","grey60","grey60", 
                                  "grey60","grey60", "grey60","grey60")) +
    scale_shape_manual(values = c(16, 17, 16, 3, 12, 8, 6, 9))
  
  
  rich_inout_fungi <-
    ggarrange(Italy_fungi_alpha +
                stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
                             geom="pointrange", 
                             color=c("black","black","black","black","black","black","red","red"), 
                             shape=18, size=0.8),
              France_fungi_alpha +  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
                                                 geom="pointrange", 
                                                 color=c("black","black","black","black","black","black","black","black","red","red"), 
                                                 shape=18, size=0.8),
              Spain_fungi_alpha +
                stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
                             geom="pointrange", 
                             color=c("black","black","black","black","black","black","red","red"), 
                             shape=18, size=0.8),
              ncol = 3,
              nrow = 1,
              align = "h", 
              widths = c(0.33, 0.37, 0.30))
  
  rich_inout_fungi
  
  return(rich_inout_fungi)
}


PlotSiteRich_EU(
  df_alpha_fungi_all_new,
  formula(hill_0 ~ site_brule),
  "Richness (hill_0)",
  350,
  355,
  "Fungi"
)


# **********************************************************************--------
# ***** REVISIONS ******--------------------------------------------------------

# **********************************************************************--------

# **** FIGURE 4 **** -----------------------------------------------------------

# Pairwise wilcox comparisons --------------------------------------------------
CompSampl <- function(df, formula){
  require(multcompView)
  require(tidyverse)
  
  df1 <-
    df %>%
    compare_means(data = ., formula, 
                  method = "wilcox.test", p.adjust.method = "BH") %>% 
    dplyr::select(group1, group2, p.adj)
  
  #print(df1)
  
  df2 <-
    df1 %>% 
    dplyr::select(group2, group1, p.adj)
  colnames(df2) <- 
    c("group1", "group2", "p.adj") 
  rbind(df1, df2) %>% 
    as.data.frame() %>% 
    xtabs(p.adj ~ group2 + group1, data=.) %>% 
    as.dist(diag = TRUE) -> dist
  res <-
    as.data.frame(multcompLetters(x=dist, reversed = FALSE)['Letters'])
  return(res)
}

CompSampl(df_alpha_fungi_all_new, formula(hill_0 ~ continent)) 


# Plot alpha diversity metrics as point-range plots ----------------------------
PlotRich <- function(dataframe,
                     X_var,
                     Y_var,
                     my_labels,
                     labels_y) {
  rich_plot <-
    ggplot(dataframe, aes(x = get(X_var), y = get(Y_var))) +
    geom_jitter(
      position = position_jitter(0.4),
      size = 1,
      aes(
        color = get(X_var),
        fill = get(X_var),
        shape = get(X_var)
      )
    ) +
    stat_summary(
      fun.data = mean_sdl,
      fun.args = list(mult = 1),
      geom = "pointrange",
      color = "black",
      shape = 18,
      size = 0.8
    ) +
    stat_summary(
      geom = 'text',
      angle = 0,
      label = my_labels,
      fun = max,
      aes(y = labels_y),
      size = 3,
      color = "black"
    ) +
    themePlot() +
    expand_limits(y = 0) +
    theme(
      axis.text.y = element_markdown(angle = 90, hjust = 0.5),
      axis.title.x = element_blank(),
      legend.position = "none"
    )
  return(rich_plot)
}


PlotRich(
  df_alpha_fungi_all_new,
  "continent",
  "hill_0",
  CompSampl(df_alpha_fungi_all_new, formula(hill_0 ~ continent)) %>%
    pull(Letters) %>% as.character(),
  600
)

# Plot venn diagrams -----------------------------------------------------------
PlotVenn <- function(physeq, Var){
  library(ggvenn)
  
  merge <- merge_samples(physeq, Var) %>% 
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
    prune_samples(sample_sums(x=.) > 0, x =.)
  
  otu <- as.data.frame(t(otu_table(merge))) %>% 
    mutate_all(~ ifelse(. > 0, TRUE, FALSE)) 
  
  ggvenn(
    otu,
    c("Australia", "Europe"),
    fill_color   = palette_cont,
    fill_alpha   = 0.5,
    stroke_size  = 0.05,
    text_size    = 3.5,
    set_name_size = 5,
    auto_scale   = TRUE
  )
}

PlotVenn(physeq_ecm_all, "continent")

# Plot mixed model means plus 95%CI horizontally (preferred) -------------------
emm_plot_sig_h <- function(model, title,
                           type = c("response","link"),
                           adjust = "none") {
  
  type <- match.arg(type)
  
  emm <- emmeans(model, ~ continent | brule, type = type)
  print(pairs(emm, by = "brule") )
  
  cld_obj <- cld(emm, by = "brule", Letters = letters, adjust = adjust)
  
  df <- as.data.frame(cld_obj)
  yvar <- if ("response" %in% names(df)) "response" else "emmean"
  df$.group <- gsub("\\s+", "", df$.group)
  
  # compute an offset so letters are to the RIGHT of the CI
  x_range <- range(c(df$lower.CL, df$upper.CL), na.rm = TRUE)
  x_pad   <- diff(x_range) * 0.05   # 5% of range; increase if needed
  
  df <- df %>%
    mutate(
      label_x = upper.CL + x_pad
    ) %>% 
    mutate(brule = recode(brule, inside = "Inside", outside = "Outside"))
  
  print(df)
  
  ggplot(df, aes(y = continent, x = .data[[yvar]])) +
    facet_wrap(~ brule, ncol = 2) +
    # CI bars (horizontal)
    geom_errorbar(
      aes(xmin = lower.CL, xmax = upper.CL),
      height = 0,
      linewidth = 1.4,
      color = "grey"
    ) +
    # Mean as diamond
    geom_point(shape = 18, size = 2.5) +
    # Letters to the right of the CI
    geom_text(aes(x = label_x, label = .group),
              vjust = 0.35, hjust = 0, size = 3) +
    labs(title = title, x = "Estimated marginal mean (95% CI)", y = NULL) +
    coord_cartesian(xlim = c(x_range[1], x_range[2] + 3*x_pad)) +
    themePlot() +
    theme(
      plot.title = element_markdown(size = 12),
      axis.title.x = element_markdown(size = 10),
      axis.title.y = element_markdown(size = 10),
      axis.text.y = element_markdown(size = 10),
      strip.background = element_rect(fill = NA, color = NA),
      strip.text = element_markdown(face = "bold", size = 10),
      panel.spacing = unit(0.9, "lines")) +
    scale_y_discrete(labels = continent_labels)
}

emm_plot_sig_h(model_rich_fungi, title = "All fungi", type = "response")

# Plot mixed model means plus 95%CI Vertically ---------------------------------
emm_plot_sig_v <- function(model, title,
                           type = c("response","link"),
                           adjust = "none") {
  
  type <- match.arg(type)
  
  emm <- emmeans(model, ~ continent | brule, type = type)
  print(pairs(emm, by = "brule"))
  
  cld_obj <- cld(emm, by = "brule", Letters = letters, adjust = adjust)
  
  df <- as.data.frame(cld_obj)
  yvar <- if ("response" %in% names(df)) "response" else "emmean"
  df$.group <- gsub("\\s+", "", df$.group)
  
  # vertical padding so letters sit above the CI
  y_range <- range(c(df$lower.CL, df$upper.CL), na.rm = TRUE)
  y_pad   <- diff(y_range) * 0.06   # increase/decrease if needed
  
  df <- df %>%
    mutate(
      label_y = upper.CL + y_pad
    ) %>%
    mutate(brule = recode(brule, inside = "Inside", outside = "Outside"))
  
  ggplot(df, aes(x = continent, y = .data[[yvar]])) +
    facet_wrap(~ brule, ncol = 2) +
    geom_errorbar(
      aes(ymin = lower.CL, ymax = upper.CL),
      width = 0, linewidth = 1.5, color = "grey") +
    # Mean as diamond
    geom_point(shape = 16, size = 2.5) +
    # Letters above the CI
    geom_text(aes(y = label_y, label = .group),
              vjust = 0, size = 3) +
    labs(title = title, y = "Estimated marginal mean (95% CI)", x = NULL) +
    coord_cartesian(ylim = c(y_range[1], y_range[2] + 3*y_pad)) +
    themePlot() +
    theme(
      axis.title = element_markdown(size = 12),
      axis.text.y = element_markdown(size = 12),
      strip.background = element_rect(fill = NA, color = NA),
      strip.text = element_markdown(face = "bold", size = 10),
      panel.spacing = unit(0.9, "lines")) +
    scale_y_discrete(labels = continent_labels)
}


emm_plot_sig_h(model_rich_ecm_2, title = "All fungi", 
               type = "response", adjust = "none")

# **********************************************************************--------

# **** FIGURE 5 **** -----------------------------------------------------------
# Extract top taxa -------------------------------------------------------------
ExtractTop <- function(physeq, rank, top){
  require(phyloseq)
  require(tidyverse)
  
  top_rank <-
    cbind(
      physeq %>%
        tax_glom(taxrank = rank) %>% 
        tax_table() %>% 
        as.matrix() %>% 
        as.data.frame() %>% 
        rownames_to_column("OTU_label") %>% 
        dplyr::select(OTU_label, Phylum, Class,
                      Order, Family, Genus, Species),
      physeq %>%
        tax_glom(taxrank = rank) %>% 
        otu_table() %>% 
        as.matrix() %>% 
        as.data.frame() %>% 
        rownames_to_column("OTU_ID") %>%
        mutate(Sum = rowSums(across(where(is.numeric)))) %>% 
        dplyr::select(OTU_ID, Sum)) %>% 
    arrange(desc(Sum)) %>% 
    slice(1:top) %>% 
    pull(rank)  
  
  return(top_rank)
}

ExtractTop(physeq_ecm_all, "Genus", 15)

# Extract clean dataframe form phyloseq object ---------------------------------

# Plot mean abudnance continent x brule ----------------------------------------
PloMedCat <- function(dataframe, X_var, Y_var, my_labels, labels_y){
  
  rich_plot <-
    ggplot(dataframe, aes(x = get(X_var), y = get(Y_var))) +
    geom_jitter(position=position_jitter(0.4), size=1, shape=1,
                aes(color=get(X_var))) +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
                 geom="pointrange", color="red", shape=18, size=0.8) +
    stat_summary(geom = 'text', angle=0, label = my_labels, fun= max, 
                 aes(y = labels_y), size=3, color="black") +
    themePlot() +
    expand_limits(y = 0) +
    theme_bw() +
    theme(plot.title = element_markdown(size = 10, face = "bold", hjust = 0.5),
          strip.text = element_markdown(size = 10, face = "bold"),
          axis.text.x = element_markdown(angle = 0, size = 8, hjust = 0.5, vjust = 1.05),
          axis.text.y = element_markdown(angle = 90, size = 8, hjust = 0.5, vjust = 0.5),
          strip.background = element_blank(),
          legend.position = "none")
  return(rich_plot)
}


# Plot negative beinomial conditinal model -------------------------------------
emm_plot_sig_h_zinb_cond <- function(model, 
                                     title = "Negative binomial conditional model",
                                     adjust = "none",
                                     type = c("response", "link"),
                                     level = 0.95) {
  type <- match.arg(type)
  require(dplyr); require(ggplot2); require(emmeans); require(multcomp)
  
  emm <- emmeans(model, ~ continent | brule,
                 component = "cond",
                 type = type)
  
  # add CIs explicitly (glmmTMB often uses asymp.* names)
  emm_ci <- as.data.frame(confint(emm, level = level))
  
  # add letters based on pairwise comparisons within brule
  cld_obj <- multcomp::cld(emm, by = "brule", Letters = letters, adjust = adjust)
  cld_df  <- as.data.frame(cld_obj) %>%
    select(continent, brule, .group)
  
  df <- emm_ci %>%
    left_join(cld_df, by = c("continent", "brule")) %>%
    mutate(.group = gsub("\\s+", "", .group))
  
  # pick x variable and CI variable names
  xvar <- if (type == "response") "response" else "emmean"
  
  lcl <- dplyr::coalesce(df$lower.CL, df$asymp.LCL)
  ucl <- dplyr::coalesce(df$upper.CL, df$asymp.UCL)
  
  df <- df %>%
    mutate(lower = lcl,
           upper = ucl)
  
  # guard against missing CIs
  if (all(is.na(df$lower)) || all(is.na(df$upper))) {
    stop("Could not find CI columns (lower/upper). Try type='link' or ensure confint() works on the emmeans object.")
  }
  
  x_range <- range(c(df$lower, df$upper), na.rm = TRUE)
  x_pad   <- diff(x_range) * 0.05
  
  df <- df %>%
    mutate(label_x = upper + x_pad,
           brule = recode(brule, inside = "Inside", outside = "Outside"))
  
  print(df)
  
  ggplot(df, aes(y = continent, x = .data[[xvar]])) +
    facet_wrap(~ brule, ncol = 1, strip.position = "right") +
    geom_errorbarh(aes(xmin = lower, xmax = upper),
                   height = 0, linewidth = 1.4, color = "grey") +
    geom_point(shape = 18, size = 2) +
    geom_text(aes(x = label_x, label = .group),
              vjust = 0.35, hjust = 0, size = 3) +
    labs(
      title = title,
      x = if (type == "response")
        "Mean abundance (given presence)"
      else
        "Conditional linear predictor (95% CI)",
      y = NULL
    ) +
    coord_cartesian(xlim = c(x_range[1], x_range[2] + 3 * x_pad)) +
    themePlot() +
    themePlot() +
    theme(
      axis.title.x = element_markdown( size = 10), 
      strip.background = element_rect(fill = NA, color = NA),
      strip.text = element_markdown(face = "bold", size = 10),
      panel.spacing = unit(0.9, "lines")
    ) +
    scale_y_discrete(labels = continent_labels)
}

emm_plot_sig_h_zinb_cond(m_Tmel_zi2)

# Plot negative binomial presence probability model ----------------------------
plot_zi_presence <- function(model,
                             title = "Negative binomial presence probability",
                             level = 0.95,
                             adjust = "none",
                             Letters = letters) {
  
  library(dplyr)
  library(ggplot2)
  library(emmeans)
  library(multcomp)
  
  # 1) ZI EMMs on response scale, then convert to presence prob
  df <- emmeans(model, ~ brule, component = "zi", type = "response") %>%
    confint(level = level) %>%
    as.data.frame() %>%
    transmute(
      SE,
      brule,
      prob  = 1 - response,
      lower = 1 - asymp.UCL,
      upper = 1 - asymp.LCL
    )
  
  # 2) Letters computed on ZI link scale (recommended inferential scale)
  emm_link <- emmeans(model, ~ brule, component = "zi", type = "link")
  cld_df <- as.data.frame(multcomp::cld(emm_link, Letters = Letters, adjust = adjust)) %>%
    select(brule, .group) %>%
    mutate(.group = gsub("\\s+", "", .group))
  
  df <- df %>%
    left_join(cld_df, by = "brule") %>%
    mutate(brule = as.factor(brule),
           brule = recode_factor(brule, inside = "Inside", outside = "Outside"), 
           brule = fct_relevel(brule, "Inside", "Outside"))
  
  # 3) place letters to the right of CI
  x_range <- range(c(df$lower, df$upper), na.rm = TRUE)
  x_pad   <- diff(x_range) * 0.06
  
  df <- df %>%
    mutate(label_x = pmin(1, upper + x_pad)) %>% 
    mutate(brule = fct_relevel(brule,  "Outside", "Inside"))
  
  print(df$brule)
  print(df)
  
  ggplot(df, aes(y = brule, x = prob)) +
    geom_errorbarh(aes(xmin = lower, xmax = upper),
                   height = 0, linewidth = 1.4, color = "grey") +
    geom_point(shape = 18, size = 2) +
    geom_text(aes(x = label_x, label = .group),
              hjust = 0, vjust = 0.35, size = 3) +
    labs(title = title,
         x = "Presence probability (1 - zero inflation)",
         y = NULL) +
    coord_cartesian(xlim = c(0, 1)) +
    themePlot() +
    theme(axis.title.x = element_markdown( size = 10))
}

plot_zi_presence(m_Tmel_zi2)

# **********************************************************************--------

# **** FIGURE 6 **** -----------------------------------------------------------
DiffAbundSet <- function(dataframe) {
  
  safe_CompSampl <- purrr::safely(
    function(df) {
      CompSampl(
        df      = df,
        formula = read_no ~ continent_brule
      )
    })
  
  out <- dataframe %>%
    tidyr::nest(data = -genus) %>%
    dplyr::mutate(
      res = purrr::map(data, safe_CompSampl),
      error = purrr::map_chr(
        res,
        function(x) {
          if (is.null(x$error)) NA_character_ else as.character(x$error)
        }
      ),
      Sig = purrr::map(res, "result")
    ) %>%
    {
      bad <- dplyr::filter(., !is.na(error))$genus
      if (length(bad) > 0) {
        message(
          "Skipping genera with failed CompSampl: ",
          paste(bad, collapse = ", ")
        )
      }
      dplyr::filter(., is.na(error))
    } %>%
    dplyr::mutate(
      Group = purrr::map(
        Sig,
        function(x) rownames(x)
      )
    ) %>%
    tidyr::unnest(cols = c(Sig, Group)) %>%
    dplyr::mutate(
      brule_cont = dplyr::recode(
        Group,
        "Australia_inside"  = "AI",
        "Australia_outside" = "AO",
        "Europe_inside"     = "EI",
        "Europe_outside"    = "EO"
      )
    )
  
  return(out)
}

all_abund_fungi_genus %>% 
  filter(!genus %in% c("Unlcassified", "unclassified")) %>% 
  group_by(genus, site, brule, continent) %>%
  summarise(read_no = mean(read_no), .groups = "drop") %>% 
  mutate(continent_brule = interaction(continent, brule, sep = "_")) %>% 
  DiffAbundSet()

# Plot mean log10(abundances) plus SD ------------------------------------------
# Only significant different genera are plotted and top_n will mark how many of those
PlotMultiDiff_sig <- function(dataframe, 
                              num_rows, 
                              top_n = NULL,
                              genera_set = NULL){
  
  # Extract factor levels with zero sum
  zero_sum <-
    dataframe %>%
    group_by(genus, continent_brule) %>%
    summarize(Sum = sum(read_no), .groups = 'drop') %>%
    filter(Sum == 0)
  
  # remove the groups for the data
  dataframe_filt <-
    dataframe %>%
    anti_join(zero_sum, by = c("genus", "continent_brule"))
  
  # set the genera I want to plot
  if (!is.null(genera_set)) {
    dataframe_filt <- dataframe_filt %>%
      filter(genus %in% genera_set)
  }
  
  # NEW: keep only genera with at least one significant comparison
  # (pairwise Wilcoxon across continent_brule; BH-adjusted)
  sig_genera <- dataframe_filt %>%
    group_by(genus) %>%
    summarise(
      pmin = tryCatch({
        pw <- pairwise.wilcox.test(read_no, continent_brule,
                                   p.adjust.method = "BH", exact = FALSE)
        min(pw$p.value, na.rm = TRUE)
      }, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    filter(is.finite(pmin), pmin < 0.05) %>%
    pull(genus)
  
  dataframe_filt <- dataframe_filt %>%
    filter(genus %in% sig_genera)
  
  # If nothing passes, return an empty plot with a message
  if (nrow(dataframe_filt) == 0) {
    return(
      ggplot() +
        theme_void() +
        labs(title = "No genera with significant pairwise Wilcoxon differences (BH < 0.05).")
    )
  }
  
  # change x axis order based on high to low abundant in total sequences
  x_order <-
    dataframe_filt %>%
    group_by(genus) %>%
    summarise(Sum = sum(read_no), .groups = "drop") %>%
    arrange(desc(Sum)) %>%
    pull(genus)
  
  print(x_order)
  
  if (!is.null(top_n)) {
    x_order <- head(x_order, top_n)
    dataframe_filt <- dataframe_filt %>%
      filter(genus %in% x_order)
  }
  
  # extract significance letters (your existing workflow)
  sig_labels <-
    dataframe_filt %>%
    group_by(genus) %>%
    filter(dplyr::n_distinct(continent_brule) > 1) %>%
    ungroup() %>%
    DiffAbundSet() %>%
    arrange(genus) %>%
    dplyr::select(genus, Group, Letters) %>%
    rename(continent_brule = Group) %>%
    mutate(
      genus = as.factor(genus),
      genus = fct_relevel(genus, x_order),
      continent_brule = as.factor(continent_brule),
      continent_brule = fct_relevel(
        continent_brule,
        "Australia_inside", "Australia_outside",
        "Europe_inside", "Europe_outside"
      ),
      brule_cont = recode(
        continent_brule,
        "Australia_inside"  = "AI",
        "Australia_outside" = "AO",
        "Europe_inside"     = "EI",
        "Europe_outside"    = "EO"
      )
    ) %>%
    arrange(genus, continent_brule)
  
  # plotting
  multidiffplot <-
    dataframe_filt %>%
    mutate(
      genus = as.factor(genus),
      genus = fct_relevel(genus, x_order),
      continent_brule = as.factor(continent_brule),
      continent_brule = fct_relevel(
        continent_brule,
        "Australia_inside", "Australia_outside",
        "Europe_inside", "Europe_outside"
      ),
      brule_cont = recode(
        continent_brule,
        "Australia_inside"  = "AI",
        "Australia_outside" = "AO",
        "Europe_inside"     = "EI",
        "Europe_outside"    = "EO"
      ),
      brule_cont = as.factor(brule_cont)
    ) %>%
    ggplot(aes(x = brule_cont, y = read_no)) +
    geom_jitter(position=position_jitter(0.4), size=1, alpha=0.5,
                aes(color = brule_cont, shape=brule_cont)) +
    stat_summary(aes(group=brule_cont), color="black", 
                 position = position_dodge(0.6),
                 fun.data=mean_sdl, fun.args = list(mult=1),
                 geom="pointrange", shape=18, size=0.7) +
    geom_text(data = sig_labels,
              aes(x = brule_cont, label = Letters, y = 5000),
              inherit.aes = FALSE,size = 3,color = "black") +
    scale_y_log10() +
    theme_bw() +
    facet_grid(~genus, scales = "free", space = "free", switch = "y") +
    #facet_wrap(~ genus,nrow = num_rows, scales = "free", space = "free_x", strip.position = "top") +
    theme(
      strip.background = element_blank(),
      plot.background = element_blank(),
      panel.background = element_blank(),
      plot.title = element_markdown(size = 12, hjust = 0.5, vjust = 0.5),
      axis.text.x = element_markdown(size = 10, angle = 90, hjust = 0.5, vjust = 0.5),
      #axis.text.y  = element_blank(),
      #axis.ticks.y = element_blank(),
      #strip.text.x = element_markdown(size = 10, angle = 0, face = "italic"),
      strip.text.x = element_markdown(size  = 10, angle  = 25, face   = "italic", hjust  = 0, vjust  = 0),
      #  margin = ggplot2::margin(l = 0, r = 0, t = 0, b = 0)), strip.placement =  "outside", strip.switch.pad.grid = unit(0, "pt"),
      legend.text = element_markdown(size = 8),
      legend.title = element_blank(),
      legend.key.height = unit(0.3, "cm"),
      legend.key.width = unit(0.3, "cm"),
      legend.position = "bottom"
    ) +
    scale_color_manual(
      values = palette_cont_brule,
      labels = c(
        "AI=Australia inside", "AO=Australia outside",
        "EI=Europe inside", "EO=Europe outside"
      )
    ) +
    #guides(
    #  color = guide_legend(order = 1, nrow = 1, title = NULL,
    #                       override.aes = list(shape = 15, size = 3))
    #) +
    labs(
      title = "Differentially abundant taxa",
      x = NULL,
      y = expression(log10(Abundance))
    ) +
    guides(shape = guide_legend(override.aes = list(size = 2, stroke = 2)))
  
  return(multidiffplot)
}

all_abund_fungi_genus %>% 
  group_by(genus, site, brule, continent) %>%
  summarise(read_no = mean(read_no), .groups = "drop") %>% 
  mutate(continent_brule = interaction(continent, brule, sep = "_")) %>% 
  PlotMultiDiff_sig(num_rows = 2,
                    genera_set = top_ecm) +
  scale_shape_manual(values = c(1,0,2,5),
                     labels = c("AI=Australia inside","AO=Australia outside",
                                "EI=Europe inside","EO=Europe outside")) +
  scale_x_discrete(labels = cont_brule_labels)



# Optional charging

PlotMultiDiff_sig_upg <- function(dataframe,
                              num_rows,
                              top_n = NULL,
                              genera_set = NULL,
                              filter_sig = TRUE,          # NEW
                              alpha = 0.05,
                              p_adjust_method = "BH"){
  
  # Extract factor levels with zero sum
  zero_sum <-
    dataframe %>%
    dplyr::group_by(genus, continent_brule) %>%
    dplyr::summarize(Sum = sum(read_no), .groups = "drop") %>%
    dplyr::filter(Sum == 0)
  
  # remove the groups for the data
  dataframe_filt <-
    dataframe %>%
    dplyr::anti_join(zero_sum, by = c("genus", "continent_brule"))
  
  # set the genera I want to plot
  if (!is.null(genera_set)) {
    dataframe_filt <- dataframe_filt %>%
      dplyr::filter(genus %in% genera_set)
  }
  
  # OPTIONAL: keep only genera with at least one significant comparison
  if (isTRUE(filter_sig)) {
    
    sig_genera <- dataframe_filt %>%
      dplyr::group_by(genus) %>%
      dplyr::summarise(
        pmin = tryCatch({
          pw <- pairwise.wilcox.test(
            x = read_no,
            g = continent_brule,
            p.adjust.method = p_adjust_method,
            exact = FALSE
          )
          suppressWarnings(min(pw$p.value, na.rm = TRUE))
        }, error = function(e) NA_real_),
        .groups = "drop"
      ) %>%
      dplyr::filter(is.finite(pmin), pmin < alpha) %>%
      dplyr::pull(genus)
    
    dataframe_filt <- dataframe_filt %>%
      dplyr::filter(genus %in% sig_genera)
    
    # If nothing passes, return an empty plot with a message
    if (nrow(dataframe_filt) == 0) {
      return(
        ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::labs(
            title = paste0(
              "No genera with significant pairwise Wilcoxon differences (",
              p_adjust_method, " < ", alpha, ")."
            )
          )
      )
    }
  }
  
  # change x axis order based on high to low abundant in total sequences
  x_order <-
    dataframe_filt %>%
    dplyr::group_by(genus) %>%
    dplyr::summarise(Sum = sum(read_no), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(Sum)) %>%
    dplyr::pull(genus)
  
  if (!is.null(top_n)) {
    x_order <- head(x_order, top_n)
    dataframe_filt <- dataframe_filt %>%
      dplyr::filter(genus %in% x_order)
  }
  
  # significance letters ALWAYS computed for the remaining genera
  sig_labels <-
    dataframe_filt %>%
    dplyr::group_by(genus) %>%
    dplyr::filter(dplyr::n_distinct(continent_brule) > 1) %>%
    dplyr::ungroup() %>%
    DiffAbundSet() %>%
    dplyr::arrange(genus) %>%
    dplyr::select(genus, Group, Letters) %>%
    dplyr::rename(continent_brule = Group) %>%
    dplyr::mutate(
      genus = as.factor(genus),
      genus = forcats::fct_relevel(genus, x_order),
      continent_brule = as.factor(continent_brule),
      continent_brule = forcats::fct_relevel(
        continent_brule,
        "Australia_inside", "Australia_outside",
        "Europe_inside", "Europe_outside"
      ),
      brule_cont = dplyr::recode(
        continent_brule,
        "Australia_inside"  = "AI",
        "Australia_outside" = "AO",
        "Europe_inside"     = "EI",
        "Europe_outside"    = "EO"
      )
    ) %>%
    dplyr::arrange(genus, continent_brule)
  
  # plotting
  multidiffplot <-
    dataframe_filt %>%
    dplyr::mutate(
      genus = as.factor(genus),
      genus = forcats::fct_relevel(genus, x_order),
      continent_brule = as.factor(continent_brule),
      continent_brule = forcats::fct_relevel(
        continent_brule,
        "Australia_inside", "Australia_outside",
        "Europe_inside", "Europe_outside"
      ),
      brule_cont = dplyr::recode(
        continent_brule,
        "Australia_inside"  = "AI",
        "Australia_outside" = "AO",
        "Europe_inside"     = "EI",
        "Europe_outside"    = "EO"
      ),
      brule_cont = as.factor(brule_cont)
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = brule_cont, y = read_no)) +
    ggplot2::geom_jitter(position = ggplot2::position_jitter(0.4),
      size = 1, alpha = 0.5,
      ggplot2::aes(color = brule_cont, shape = brule_cont)
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(group = brule_cont),
      color = "black",
      position = ggplot2::position_dodge(0.6),
      fun.data = mean_sdl, fun.args = list(mult = 1),
      geom = "pointrange", shape = 18, size = 0.7
    ) +
    ggplot2::geom_text(
      data = sig_labels,
      ggplot2::aes(x = brule_cont, label = Letters, y = 5000),
      inherit.aes = FALSE, size = 3, color = "black"
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(~genus, scales = "free", space = "free", switch = "y") +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.title = ggtext::element_markdown(size = 12, face = "bold",
                                            hjust = 0.5, vjust = 0.5),
      axis.text.x = ggtext::element_markdown(size = 10, angle = 90, hjust = 0.5, vjust = 0.5),
      strip.text.x = ggtext::element_markdown(size = 10, angle = 25, face = "italic", hjust = 0, vjust = 0),
      legend.text = ggtext::element_markdown(size = 8),
      legend.title = ggplot2::element_blank(),
      legend.key.height = grid::unit(0.3, "cm"),
      legend.key.width = grid::unit(0.3, "cm"),
      legend.position = "bottom"
    ) +
    ggplot2::scale_color_manual(
      values = palette_cont_brule,
      labels = c(
        "AI=Australia inside", "AO=Australia outside",
        "EI=Europe inside", "EO=Europe outside"
      )
    ) +
    ggplot2::labs(
      title = "Differentially abundant taxa",
      x = NULL,
      y = expression(log10(Abundance))
    ) +
    ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 2, stroke = 2)))
  
  return(multidiffplot)
}



all_abund_fungi_genus %>%
  group_by(genus, site, brule, continent) %>%
  summarise(read_no = mean(read_no), .groups = "drop") %>%
  mutate(continent_brule = interaction(continent, brule, sep = "_")) %>%
  PlotMultiDiff_sig_upg(
    num_rows = 2,
    genera_set = top_ecm,
    filter_sig = FALSE
  )



# Identify indicators ----------------------------------------------------------
GetIndicators <-function(physeq, var, func){
  require(phyloseq)
  require(indicspecies)
  require(tidyverse)
  
  otu <- as.data.frame(otu_table(physeq))
  metadata = as(sample_data(physeq), "data.frame")
  taxa <- as.data.frame(as.matrix(tax_table(physeq)))
  
  perm_isa <- 
    how(blocks = metadata$site, nperm = 999, observed = TRUE)

  multipatt <- multipatt(t(otu), 
                         metadata[,var], 
                         func = func,
                         control = perm_isa,
                         duleg=TRUE)
  isa_list <-
    multipatt$sign %>% 
    mutate(p.adjust = p.adjust(p.value, "fdr")) %>% 
    rownames_to_column("OTU_ID") %>% 
    full_join(x=., 
              y= taxa %>% rownames_to_column("OTU_ID"),
              by = "OTU_ID") %>% 
    filter(p.value <= 0.05)
  
  return(isa_list)
}

physeq_fungi_all %>% 
  subset_taxa(Species!="Tuber melanosporum") %>% 
  GetIndicators(., "Clusters", "r.g")


# Selecting on subset data sets ------------------------------------------------
selectFeat <- function(df, formula){
  require(Boruta)
  
  sel_otu <- vector(mode = "character")
  sel_attr <- list()
  
  for(i in 1:3) {
    sel_attr[[i]] <-Boruta(formula, df, pValue = 0.05, 
                           mcAdj = TRUE, maxRuns=100, doTrace = 3)
    print(sel_attr[[i]])
    sel_otu <- append(sel_otu, getSelectedAttributes(sel_attr[[i]], withTentative = TRUE))
  }
  return(sel_otu)
}


# Built data frame forfeature selection in Boruta ------------------------------
df4Boruta <- function(physeq, Var){
  
  require(tidyverse)
  
  df_for_boruta <-
    otu_table(
      physeq %>% 
        subset_taxa(Species!="Tuber melanosporum"))  %>%
    t(x = .) %>%
    as.data.frame(x = .) %>% 
    rownames_to_column("SampleID") %>% 
    left_join(x=., 
              y=physeq@sam_data %>% 
                as.matrix() %>% 
                as.data.frame(),
              by = "SampleID") %>% 
    dplyr::select(contains("FOTU"), c(Var))
  
  return(df_for_boruta)
}

df4Boruta(physeq_fungi_all, "clusters")


df4Boruta_bact <- function(physeq, Var){
  
  require(tidyverse)
  
  df_for_boruta <-
    otu_table(physeq) %>%
    t(.) %>%
    as.data.frame() %>% 
    rownames_to_column("SampleID") %>% 
    left_join(x=., 
              y=physeq@sam_data %>% 
                as.matrix() %>% 
                as.data.frame(),
              by = "SampleID") %>% 
    dplyr::select(contains("BOTU"), c(Var))
  
  return(df_for_boruta)
}


df4Plotting <- function(physeq, RF_otu, ind_rg){
  require(phyloseq)
  require(tidyverse)
  
  df_for_plot <-
    left_join(
      physeq@otu_table %>% 
        as.data.frame() %>% 
        rownames_to_column("OTU_ID") %>% 
        filter(OTU_ID %in% RF_otu) %>% 
        pivot_longer(cols = -OTU_ID, 
                     names_to = "Sample_ID", 
                     values_to = "Count"),
      physeq@tax_table %>% 
        as.matrix() %>% 
        as.data.frame() %>% 
        rownames_to_column("OTU_ID") %>% 
        mutate(across(
          c("Kingdom", "Phylum", "Class","Order",
            "Family","Genus", "Species"), 
          ~na_if(., "Unlcassified"))) %>% 
        ReformatTaxonomy(.),
      by="OTU_ID") %>% 
    left_join(., ind_rg %>% 
                dplyr::select(OTU_ID, index) %>% 
                rename("r.g" = index) %>%
                mutate(r.g = 
                         recode_factor(r.g, 
                                       "1" = "High",
                                       "2" = "Mid",
                                       "3" = "Low")),
              by="OTU_ID")
  
  return(df_for_plot)
}

df4Plotting(physeq_fungi_all, RF_res_fungi, ind_fungi_Tmel_rg)

# **********************************************************************--------