# ************ DATA ANALYSIS ************ -----------------------------------------------------
# Project name: Fungal and Bacterial microbiome of swithgrass  
# Manuscript:   Fungal and bacterial community variation in switchgrass soils and roots
# Authors:      Lukas, ..., Sarah
# Affiliation:  GLBRC - Michigan State University
# Journal:      
# Date:         February 7, 2022
# ****************************************** --------------------------------------------------

# WORKING ENVIRONMENT SETUP -------------------------------------------------------------------
options(scipen = 9999)
options(max.print=100000000)
#rm(list = ls())

# loading required packages -------------------------------------------------------------------
library(phyloseq); packageVersion("phyloseq")
library(Biostrings)
library(ggplot2)
library(ape)
library(dplyr)
library(ggpubr)
library(vegan)
library(magrittr)
library(Misfunk)

# Import datasets ------------------------------------------------------------------------------
load("GLBRC018_OTU_bact_MMPRNT_mock_bact_rar.RData")
load("GLBRC018_OTU_fung_MMPRNT_mock_fung_rar.RData")

load("GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root.RData")
load("GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil.RData")
load("GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root.RData")
load("GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil.RData")

load("GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5.RData")
load("GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil.RData")
ls()

# MATCH SAMPLES -------------------------------------------------------------------------------
# Match samples to the same ones for fungi and bacteria ---------------------------------------

# All sites
# root
identical(
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root@sam_data$sampleID_long,
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root@sam_data$sampleID_long
)


intersect(
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root@sam_data$sampleID_long,
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root@sam_data$sampleID_long
) -> good_root
good_root

physeq_fungi_root <-
  subset_samples(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_root,
                 sampleID_long %in% good_root)
physeq_fungi_root@otu_table <-
  physeq_fungi_root@otu_table[which(rowSums(physeq_fungi_root@otu_table) > 0), ]
physeq_fungi_root

sample_names(physeq_fungi_root) <-
  physeq_fungi_root@sam_data$sampleID_long
physeq_fungi_root@sam_data
head(physeq_fungi_root@otu_table)

physeq_bact_root <-
  subset_samples(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_root,
                 sampleID_long %in% good_root)
physeq_bact_root@otu_table <-
  physeq_bact_root@otu_table[which(rowSums(physeq_bact_root@otu_table) > 0), ]
physeq_bact_root

sample_names(physeq_bact_root) <-
  physeq_bact_root@sam_data$sampleID_long
physeq_bact_root@sam_data

# soil
identical(
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil@sam_data$sampleID_long,
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil@sam_data$sampleID_long
)


intersect(
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil@sam_data$sampleID_long,
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil@sam_data$sampleID_long
) -> good_soil

GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil@sam_data$sampleID_long %in%
  GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil@sam_data$sampleID_long

GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil@sam_data$sampleID_long %in%
  GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil@sam_data$sampleID_long

good_soil
length(good_soil)

physeq_fungi_soil <-
  subset_samples(GLBRC018_OTU_fung_MMPRNT_All_sites_G5_soil,
                 sampleID_long %in% good_soil)
physeq_fungi_soil@otu_table <-
  physeq_fungi_soil@otu_table[which(rowSums(physeq_fungi_soil@otu_table) > 0), ]
physeq_fungi_soil

sample_names(physeq_fungi_soil) <-
  physeq_fungi_soil@sam_data$sampleID_long
physeq_fungi_soil@sam_data

physeq_bact_soil <-
  subset_samples(GLBRC018_OTU_bact_MMPRNT_All_sites_G5_soil,
                 sampleID_long %in% good_soil)
physeq_bact_soil@otu_table <-
  physeq_bact_soil@otu_table[which(rowSums(physeq_bact_soil@otu_table) > 0),]
physeq_bact_soil

sample_names(physeq_bact_soil) <-
  physeq_bact_soil@sam_data$sampleID_long
physeq_bact_soil@sam_data


# Lux Arbor time series
head(GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil@sam_data)
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil@sam_data$sampType
GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil@sam_data$sampleID_long

head(GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5@sam_data)
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5@sam_data$Root_soil
GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5@sam_data$sampleID_long


physeq_bact_LUX_root <-
  GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil %>%
  subset_samples(sampType %in% "Root") %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .)

physeq_fungi_LUX_root <-
  GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5 %>%
  subset_samples(Root_soil %in% "Root") %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .)


physeq_bact_LUX_soil <-
  GLBRC018_OTU_bact_MMPRNT_rar_LUX_G5_fil %>%
  subset_samples(sampType %in% "Soil") %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .)

physeq_fungi_LUX_soil <-
  GLBRC018_OTU_fung_MMPRNT_rar_LUX_G5 %>%
  subset_samples(Root_soil %in% "Soil") %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .)


intersect(
  physeq_bact_LUX_root@sam_data$sampleID_long,
  physeq_fungi_LUX_root@sam_data$sampleID_long
) -> good_time_root
length(good_time_root)

intersect(
  physeq_bact_LUX_soil@sam_data$sampleID_long,
  physeq_fungi_LUX_soil@sam_data$sampleID_long
) -> good_time_soil
length(good_time_soil)


physeq_bact_LUX_root_filt <-
  subset_samples(physeq_bact_LUX_root, sampleID_long %in% good_time_root) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .)
physeq_bact_LUX_root_filt

physeq_fungi_LUX_root_filt <-
  subset_samples(physeq_fungi_LUX_root, sampleID_long %in% good_time_root) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .)
physeq_fungi_LUX_root_filt

physeq_bact_LUX_soil_filt <-
  subset_samples(physeq_bact_LUX_soil, sampleID_long %in% good_time_soil) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .)
physeq_bact_LUX_soil_filt

physeq_fungi_LUX_soil_filt <-
  subset_samples(physeq_fungi_LUX_soil, sampleID_long %in% good_time_soil) %>%
  prune_taxa(taxa_sums(x = .) > 0, x = .)
physeq_fungi_LUX_soil_filt

sum(duplicated(sample_data(physeq_bact_LUX_root_filt)$sampleID_long))
sum(duplicated(sample_data(physeq_fungi_LUX_root_filt)$sampleID_long))
sum(duplicated(sample_data(physeq_bact_LUX_soil_filt)$sampleID_long))
sum(duplicated(sample_data(physeq_fungi_LUX_soil_filt)$sampleID_long))

sample_names(physeq_bact_LUX_root_filt) <-
  physeq_bact_LUX_root_filt@sam_data$sampleID_long
physeq_bact_LUX_root_filt@sam_data

sample_names(physeq_fungi_LUX_root_filt) <-
  physeq_fungi_LUX_root_filt@sam_data$sampleID_long
physeq_fungi_LUX_root_filt@sam_data

sample_names(physeq_bact_LUX_soil_filt) <-
  physeq_bact_LUX_soil_filt@sam_data$sampleID_long
physeq_bact_LUX_soil_filt@sam_data

sample_names(physeq_fungi_LUX_soil_filt) <-
  physeq_fungi_LUX_soil_filt@sam_data$sampleID_long
physeq_fungi_LUX_soil_filt@sam_data

# ********************************-------------------------------------------------------------
# PROCRUSTES ROTATION ANALYSIS ----------------------------------------------------------------
library(ade4)
library(vegan)

# calculating NMDS ----------------------------------------------------------------------------
set.seed(2021)
nmds_its_root = ordinate(physeq_fungi_root, method ="NMDS", distance="bray")
nmds_its_root
stressplot(nmds_its_root)

nmds_16s_root = ordinate(physeq_bact_root, method ="NMDS", distance="bray")
nmds_16s_root

nmds_its_soil = ordinate(physeq_fungi_soil, method ="NMDS", distance="bray")
nmds_its_soil

nmds_16s_soil = ordinate(physeq_bact_soil, method ="NMDS", distance="bray")
nmds_16s_soil

# LUX arbor Bacteria and Fungi
nmds_16s_soil_time <- ordinate(
  physeq_bact_LUX_soil_filt,
  method = "NMDS",
  distance = "bray",
  autotransform = TRUE,
  trymax = 100
)
nmds_16s_soil_time
stressplot(nmds_16s_soil_time)
plot(nmds_16s_soil_time, "sites")

nmds_16s_root_time <- ordinate(
  physeq_bact_LUX_root_filt,
  method = "NMDS",
  distance = "bray",
  autotransform = TRUE,
  trymax = 100
)
nmds_16s_root_time
stressplot(nmds_16s_root_time)

nmds_its_soil_time <- ordinate(
  physeq_fungi_LUX_soil_filt,
  method = "NMDS",
  distance = "bray",
  autotransform = TRUE,
  trymax = 100
)
nmds_its_soil_time
plot(nmds_its_soil_time, "sites")
stressplot(nmds_its_soil_time)

nmds_its_root_time <- ordinate(
  physeq_fungi_LUX_root_filt,
  method = "NMDS",
  distance = "bray",
  autotransform = TRUE,
  trymax = 100
)
nmds_its_root_time
stressplot(nmds_its_root_time)
plot(nmds_its_root_time, "sites")


# Extract vectors and incorporate sample data -------------------------------------------------
df_its_root <- data.frame(nmds_its_root$points)
df_16s_root <- data.frame(nmds_16s_root$points)

df_its_root <-
  left_join(
    tibble::rownames_to_column(
      as(physeq_fungi_root@sam_data, "data.frame") %>%
        dplyr::select(sampleID_long, siteID, FertStatus)
    ),
    tibble::rownames_to_column(df_its_root),
    by = "rowname"
  )

df_16s_root <-
  left_join(
    tibble::rownames_to_column(
      as(physeq_bact_root@sam_data, "data.frame") %>%
        dplyr::select(sampleID_long, siteID, FertStatus)
    ),
    tibble::rownames_to_column(df_16s_root),
    by = "rowname"
  )

rownames(df_its_root) <- df_its_root$rowname
df_its_root

rownames(df_16s_root) <- df_16s_root$rowname
df_16s_root


df_its_soil <- data.frame(nmds_its_soil$points)
df_16s_soil <- data.frame(nmds_16s_soil$points)

df_its_soil <-
  left_join(
    tibble::rownames_to_column(
      as(physeq_fungi_soil@sam_data, "data.frame") %>%
        dplyr::select(sampleID_long, siteID, FertStatus)
    ),
    tibble::rownames_to_column(df_its_soil),
    by = "rowname"
  )

df_16s_soil <-
  left_join(
    tibble::rownames_to_column(
      as(physeq_bact_soil@sam_data, "data.frame") %>%
        dplyr::select(sampleID_long, siteID, FertStatus)
    ),
    tibble::rownames_to_column(df_16s_soil),
    by = "rowname"
  )

rownames(df_its_soil) <- df_its_soil$rowname
df_its_soil

rownames(df_16s_soil) <- df_16s_soil$rowname
df_16s_soil


# Lux Arbor 
df_its_root_time <- data.frame(nmds_its_root_time$points)
df_16s_root_time <- data.frame(nmds_16s_root_time$points)

df_its_root_time <-
  left_join(
    tibble::rownames_to_column(
      as(physeq_fungi_LUX_root_filt@sam_data, "data.frame") %>%
        dplyr::select(sampleID_long, collectionDate, FertStatus)
    ),
    tibble::rownames_to_column(df_its_root_time),
    by = "rowname"
  )

df_16s_root_time <-
  left_join(
    tibble::rownames_to_column(
      as(physeq_bact_LUX_root_filt@sam_data, "data.frame") %>%
        dplyr::select(sampleID_long, collectionDate, FertStatus)
    ),
    tibble::rownames_to_column(df_16s_root_time),
    by = "rowname"
  )

rownames(df_its_root_time) <- df_its_root_time$rowname
df_its_root_time$Time <-
  as.Date(df_its_root_time$collectionDate, format = "%m/%d/%Y")

rownames(df_16s_root_time) <- df_16s_root_time$rowname
df_16s_root_time


df_its_soil_time <- data.frame(nmds_its_soil_time$points)
df_16s_soil_time <- data.frame(nmds_16s_soil_time$points)

df_its_soil_time <-
  left_join(
    tibble::rownames_to_column(
      as(physeq_fungi_LUX_soil_filt@sam_data, "data.frame") %>%
        dplyr::select(sampleID_long, collectionDate, FertStatus)
    ),
    tibble::rownames_to_column(df_its_soil_time),
    by = "rowname"
  )

df_16s_soil_time <-
  left_join(
    tibble::rownames_to_column(
      as(physeq_bact_LUX_soil_filt@sam_data, "data.frame") %>%
        dplyr::select(sampleID_long, collectionDate, FertStatus)
    ),
    tibble::rownames_to_column(df_16s_soil_time),
    by = "rowname"
  )

rownames(df_its_soil_time) <- df_its_soil_time$rowname
df_its_soil_time

rownames(df_16s_soil_time) <- df_16s_soil_time$rowname
df_16s_soil_time


df_its_root_time$Time <-
  as.Date(df_its_root_time$collectionDate, format = "%m/%d/%Y")
df_its_soil_time$Time <-
  as.Date(df_its_soil_time$collectionDate, format = "%m/%d/%Y")
df_16s_root_time$Time <-
  as.Date(df_16s_root_time$collectionDate, format = "%m/%d/%Y")
df_16s_soil_time$Time <-
  as.Date(df_16s_soil_time$collectionDate, format = "%m/%d/%Y")

# rearranging the sample order ----------------------------------------------------------------
identical(rownames(df_its_root), rownames(df_16s_root))
order_root <- match(rownames(df_16s_root), rownames(df_its_root))
df_its_root <- df_its_root[order_root, ]

identical(rownames(df_its_soil), rownames(df_16s_soil))
order_soil <- match(rownames(df_16s_soil), rownames(df_its_soil))
df_its_soil <- df_its_soil[order_soil, ]


identical(rownames(df_its_root_time), rownames(df_16s_root_time))
order_root_time <-
  match(rownames(df_16s_root_time), rownames(df_its_root_time))
df_its_root_time <- df_its_root_time[order_root_time, ]

identical(rownames(df_its_soil_time), rownames(df_16s_soil_time))
order_soil_time <-
  match(rownames(df_16s_soil_time), rownames(df_its_soil_time))
df_its_soil_time <- df_its_soil_time[order_soil_time, ]


# Run procrustes analysis ---------------------------------------------------------------------
pro_root <- protest(
  df_its_root[, c("MDS1", "MDS2")],
  df_16s_root[, c("MDS1", "MDS2")],
  scores = "sites",
  permutations = how(nperm = 9999)
)
pro_root
pro_root$ss

pro_soil <- protest(
  df_its_soil[, c("MDS1", "MDS2")],
  df_16s_soil[, c("MDS1", "MDS2")],
  scores = "sites",
  permutations = how(nperm = 9999)
)
pro_soil
pro_soil$ss

# Lux Arbor 
pro_root_time <- protest(
  df_its_root_time[, c("MDS1", "MDS2")],
  df_16s_root_time[, c("MDS1", "MDS2")],
  scores = "sites",
  permutations = how(nperm = 9999)
)
pro_root_time
pro_root_time$ss

pro_soil_time <- protest(
  df_its_soil_time[, c("MDS1", "MDS2")],
  df_16s_soil_time[, c("MDS1", "MDS2")],
  scores = "sites",
  permutations = how(nperm = 9999)
)
pro_soil_time
pro_soil_time$ss


# Plotting ------------------------------------------------------------------------------------
require(scales)

PlotProcrust <- function(pro_df, ord, Var){
  pro_plot <- data.frame(
    pcoa1 = pro_df$Yrot[, 1],
    pcoa2 = pro_df$Yrot[, 2],
    xpcoa1 = pro_df$X[, 1],
    xpcoa2 = pro_df$X[, 2]
  )
  if (identical(rownames(pro_plot), rownames(ord))==TRUE) {
    pro_plot$FertStatus <- ord$FertStatus
    pro_plot$siteID <- ord$siteID
    pro_plot$Time <- ord$Time
    # time_range <- range(ord$Time)
    # pro_plot$collectionDate <-
    #   as.Date(ord$Time, origin = time_range[1])
    head(pro_plot) %>% print()
    ggplot(pro_plot) +
      theme(legend.key.height = unit(0.35, "cm"), legend.key.width = unit(0.3, "cm")) +
      theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
      theme(strip.text.x = element_text(size = 9, face = "bold")) +
      theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5)) + 
      theme(axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
      theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
      theme(axis.title = element_text(angle = 0, size = 9, face = "bold")) +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.05)) + # to adjust decimals
      scale_y_continuous(labels = scales::number_format(accuracy = 0.05)) +
      scale_shape_manual("", 
                         values = c(19, 2), 
                         labels=c("Fert","Unfert")) +
      grids(linetype = "dashed") +
      theme(panel.background = element_blank(),
            legend.position = "bottom",
            legend.background = element_blank(),
            legend.key = element_blank(),
            axis.line = element_line()) +
      geom_point(aes(x=pcoa1, y=pcoa2, 
                     colour=get(Var), shape=FertStatus), size=1.5, alpha=0.9) +
      geom_point(aes(x=xpcoa1, y=xpcoa2, 
                     colour=get(Var), shape=FertStatus), size=1.5, alpha=0.9) +
      geom_segment(aes(x=xpcoa1,y=xpcoa2,xend=pcoa1,yend=pcoa2,
                       colour=get(Var)), size=0.3,
                   arrow=arrow(length=unit(0.3,"cm")), show.legend=F) +
      geom_hline(yintercept = 0,linetype="dashed") +
      geom_vline(xintercept = 0,linetype="dashed") +
      geom_abline(slope = -1/pro_df$rotation[1,2]) +
      geom_abline(slope = pro_df$rotation[1,2]) -> plot_proc
    return(plot_proc)
  }else{
    stop("Error: pro_plot and ord rownames are not matching!")
  }
}


# Test
PlotProcrust(pro_root, df_its_root, "siteID") + 
  labs(title = "Procrustes Rotation Root\nFungi - Bacteria", x="Axis.1", y="Axis.2") +
  annotate("text", Inf, -Inf, label = "italic(m) ^ 2 == 0.395", parse = TRUE, size = 3, hjust = 1.1, vjust = -0.5) +
  annotate("text", Inf, -Inf, label = "italic(r) ^ 2 == 0.778", parse = TRUE, size = 3, hjust = 1.1, vjust = -1.6) +
  theme(legend.position = "top",
        legend.margin = margin(0,-0.3, 0, 0, unit = "cm")) +
  guides(color = guide_legend(ncol = 1), #title.position="top"
         shape = guide_legend(ncol = 1)) +
  scale_color_brewer(palette = "Set1",
                     name=NULL,
                     labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhineland")) 


PlotProcrust(pro_root_time, df_its_root_time, "Time") +
  labs(title = "Lux Arbor Procrustes Rotation Root\nFungi - Bacteria", x="Axis.1", y="Axis.2") +
  scale_color_gradient(name=NULL, trans = "date", low = "#B2DF8A", high = "#33A02C", 
                       guide = guide_colourbar(direction = "vertical"))  +
  theme(legend.position = "right",
        legend.margin = margin(0,-0.1, 0, 0, unit = "cm")) +
  guides(shape = guide_legend(ncol = 1)) +
  annotate("text", Inf, -Inf, label = "italic(m) ^ 2 == 0.526", parse = TRUE, size = 3, hjust = 1.1, vjust = -0.5)+
  annotate("text", Inf, -Inf, label = "italic(r) ^ 2 == 0.689", parse = TRUE, size = 3, hjust = 1.1, vjust = -1.6)

# Playing with breaks
time_breaks_roots <- c("2018-05-29", "2018-06-25", "2018-07-30","2018-08-20", "2018-09-17", "2018-10-03")

time_breaks_soil <- c( "2018-03-19","2018-04-30", "2018-05-15", "2018-05-29", "2018-06-11", 
                       "2018-06-25", "2018-07-09", "2018-07-30", "2018-08-08", "2018-08-20",
                       "2018-09-04", "2018-09-17", "2018-10-03", "2018-10-15", "2018-11-05")

time_breaks_soil <- c( "2018-03-19","2018-04-30", "", "2018-05-29", "", 
                       "2018-06-25", "", "2018-07-30", "", "2018-08-20",
                       "", "2018-09-17", "", "2018-10-15", "2018-11-05")

time_breaks_soil_labels <- c("Mar", "Apr","" ,"May", "", 
                             "Jun","", "Jul","", "Aug",
                             "", "Sep","", "Oct", "Nov")

# Colors
root_time_color <- colorRampPalette(c("#B2DF8A", "#33A02C"))(6)
soil_time_color <-  colorRampPalette(c("#FFFF99", "#B15928"))(6)

# plots
ggarrange(
  PlotProcrust(pro_root, df_its_root, "siteID") + 
    labs(title = "All sites root communitues", x="Axis.1", y="Axis.2") +
    annotate("text", -Inf, Inf, label = "italic(m) ^ 2 == 0.395", parse = TRUE, size = 3, hjust = -0.1, vjust = 1.5) +
    annotate("text", -Inf, Inf, label = "italic(r) == 0.778", parse = TRUE, size = 3, hjust = -0.1, vjust = 3.5) +
    theme(legend.position = "top",
          legend.margin = margin(0, 0.1, 0, 0, unit = "cm")) +
    guides(color = guide_legend(ncol = 1, order = 1), 
           shape = guide_legend(ncol = 1, order = 2)) +
    scale_color_brewer(palette = "Set1",
                       name=NULL,
                       labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhineland")),
  PlotProcrust(pro_soil, df_its_soil, "siteID") +
    labs(title = "All sites soil communitues", x="Axis.1", y="Axis.2") +
    annotate("text", -Inf, Inf, label = "italic(m) ^ 2 == 0.267", parse = TRUE, size = 3, hjust = -0.1, vjust = 5) +
    annotate("text", -Inf, Inf, label = "italic(r) == 0.856", parse = TRUE, size = 3, hjust = -0.1, vjust = 8) +
    theme(legend.position = "top",
          legend.margin = margin(0, 0.1, 0, 0, unit = "cm")) +
    guides(color = guide_legend(ncol = 1, order = 1), 
           shape = guide_legend(ncol = 1, order = 2)) +
    scale_color_brewer(palette = "Set1",
                       name=NULL,
                       labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhineland")),
  labels = c("a)","b)"),
  widths = c(1, 1),
  align = "hv" ,
  common.legend = TRUE,
  ncol = 1, nrow = 2, 
  legend = c("right")) -> Fig_X_A

Fig_X_A

# test - ask if people like it
PlotProcrust(pro_soil_time, df_its_soil_time, "Time") +
  scale_color_gradient(name=NULL, trans = "date", low = "#FFFF99",high = "#B15928", 
                       breaks = as.Date(time_breaks_soil), 
                       guide = guide_colourbar(direction = "vertical"),
                       labels = time_breaks_soil_labels) +
  theme(legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.3, "cm"),
        legend.position = "right",
        legend.margin = margin(0,-0.1, 0, 0, unit = "cm")) +
  guides(shape = guide_legend(ncol = 1)) +
  annotate("text", -Inf, Inf, label = "italic(m) ^ 2 == 0.537", parse = TRUE, size = 3, hjust = -0.1, vjust = 1.5) +
  annotate("text", -Inf, Inf, label = "italic(r) == 0.681", parse = TRUE, size = 3, hjust = -0.1, vjust = 3.5)


# standard plot
ggarrange(
  PlotProcrust(pro_root_time, df_its_root_time, "Time") +
    labs(title = "Lux Arbor root communitues", x="Axis.1", y="Axis.2") +
    scale_color_gradient(name=NULL, trans = "date", low = "#B2DF8A", high = "#33A02C", 
                         guide = guide_colourbar(direction = "vertical"))  +
    theme(legend.position = "right",
          legend.margin = margin(0,-0.1, 0, 0, unit = "cm")) +
    guides(shape = guide_legend(ncol = 1)) +
    annotate("text", -Inf, Inf, label = "italic(m) ^ 2 == 0.526", parse = TRUE, size = 3, hjust = -0.1, vjust = 1.5) +
    annotate("text", -Inf, Inf, label = "italic(r) == 0.689", parse = TRUE, size = 3, hjust = -0.1, vjust = 3.5),
  PlotProcrust(pro_soil_time, df_its_soil_time, "Time") +
    labs(title = "Lux Arbor soil communitues", x="Axis.1", y="Axis.2") +
    scale_color_gradient(name=NULL, trans = "date", low = "#FFFF99",high = "#B15928", 
                         guide = guide_colourbar(direction = "vertical"))  +
    theme(legend.position = "right",
          legend.margin = margin(0,-0.1, 0, 0, unit = "cm")) +
    guides(shape = guide_legend(ncol = 1)) +
    annotate("text", -Inf, Inf, label = "italic(m) ^ 2 == 0.537", parse = TRUE, size = 3, hjust = -0.1, vjust = 1.5) +
    annotate("text", -Inf, Inf, label = "italic(r) == 0.681", parse = TRUE, size = 3, hjust = -0.1, vjust = 3.5),
  labels = c("c)","d)"),
  widths = c(1, 1),
  align = "hv" ,
  common.legend = FALSE,
  ncol = 1, nrow = 2, 
  legend = c("right")) -> Fig_X_B


Fig_X_B

save(Fig_X_A, Fig_X_B,
     pro_root, pro_soil, 
     df_its_root, df_its_soil,
     pro_root_time, pro_soil_time,
     df_its_root_time, df_its_soil_time, 
     file = "Fig_S1_rotation.Rdata")

# *** FIGURE 1 - Procustes Rotation Ordinations -----------------------------------------------
Fig_procrust <-
  ggarrange(Fig_X_A,
            Fig_X_B,
            heights = c(1,1),
            widths = c(1,1),
            ncol = 2,
            nrow = 1
  )

Fig_procrust

library(gridExtra)
title1 = text_grob("Procrustes Rotations", size = 12, face = 2)
grid.arrange(Fig_procrust, top=title1)


# m2 statistic provides an overall measure of the concordance between two data sets 
# (i.e. how close the two data configuration match),
# Procrustes analysis allows one to determine how much variance in one matrix is 
# attributable to the variance in the other. Further, visual inspection of a 
# Procrustes plot, in which the residuals between points from each matrix are
# mapped, can allow the identification of individual objects that have (relatively)
# unusual concordance.

# PROCRUSTES ERRORS -----------------------------------------------------------------------
# Residual Plot - This allows identification of samples with the worst fit. The horizontal
# lines, from bottom to top, are the 25% (dashed), 50% (solid), and 75% (dashed) quantiles 
# of the residuals. 

# General function for residual plot exploration
ProcrustRes <- function(pro_df, ord){
  as.data.frame(residuals(pro_df)) -> res_df
  colnames(res_df) <- "Res"
  res_df$FertStatus <- ord$FertStatus
  res_df$siteID <- ord$siteID 
  res_df %T>% print()
  res_df$SampleID <- rownames(res_df)
  res_df$siteID <- dplyr::recode(res_df$siteID,
                                 ESC = "Escanaba",
                                 HAN = "Hancock",
                                 LC = "Lake City",
                                 LUX ="Lux Arbor",
                                 RHN = "Rhineland")
  if (identical(rownames(res_df), rownames(ord))==TRUE){
    ggplot(res_df, aes(x=SampleID, y=Res, fill = FertStatus)) + 
      geom_bar(stat="identity")+
      theme_classic() +
      theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
      theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
      theme(strip.text.x = element_text(size = 9, face = "bold")) +
      theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
      theme(axis.text.x = element_text(angle = 90, size = 5.5, hjust = 1, vjust = 0.5)) +
      theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
      theme(axis.title = element_text(angle = 0, size = 9, face = "bold")) +
      ylim(0, 0.18) +
      labs(y="Procrustes residual error", x="Samples") +
      geom_hline(yintercept = quantile(res_df$Res)[2], linetype="dashed", color="grey30") +
      geom_hline(yintercept = quantile(res_df$Res)[3], linetype="solid", color="grey30") +
      geom_hline(yintercept = quantile(res_df$Res)[4], linetype="dashed", color="grey30") +
      theme(legend.position = "bottom") -> re_plot
    return(re_plot)
  }else{
    stop("Error: pro_plot and ord rownames are not matching!")
  }
}

ProcrustRes(pro_root, df_its_root) + 
  scale_fill_grey("") +
  facet_grid(~siteID + FertStatus, scales = "free_x", space="free_x") 

# Plot Procrustes errors for all samples ------------------------------------------------------
ggarrange(ProcrustRes(pro_root, df_its_root) +
            facet_grid(~siteID, scales = "free_x", space="free_x") +
            labs(title = "Root samples procrustes errors", x=NULL)+ 
            scale_fill_grey(""), 
          ProcrustRes(pro_soil, df_its_soil) +
            facet_grid(~siteID, scales = "free_x", space="free_x")+ 
            labs(title = "Soil sampless procrustes errors")+ 
            scale_fill_grey(""),
          labels = c("A","B"),
          widths = c(1, 1),
          align = "hv",
          ncol = 1, 
          nrow = 2, 
          common.legend = TRUE,
          legend = c("bottom")) -> Fig_x_proc_err

Fig_x_proc_err

# # if you use Market as a variable
# ProcrustRes(pro_its_lsu, ord_its) + scale_fill_manual(values=paletteCB4)
# ProcrustRes(pro_its_16s, ord_its) + scale_fill_manual(values=paletteCB4)
# ProcrustRes(pro_16s_lsu, ord_its) + scale_fill_manual(values=paletteCB4)

# ERRORS COMPARISON ANALYSIS -----------------------------------------------------------------

# extracting dataset
ExtrProcDf <- function(procdf, df, Var){
  res_df <- as.data.frame(residuals(procdf))
  colnames(res_df) <- "Res_Err"
  res_df$Sample_ID <- rownames(res_df)
  identical(rownames(res_df), rownames(res_df))
  order_proc_its_lsu <- match(rownames(df), rownames(res_df))
  res_df <- res_df[order_proc_its_lsu,]
  res_df$FertStatus <- df$FertStatus
  if (Var == "site"){
    res_df$siteID <- df$siteID
    res_df$siteID <- dplyr::recode(res_df$siteID,
                                   ESC = "Escanaba",
                                   HAN = "Hancock",
                                   LC = "Lake City",
                                   LUX ="Lux Arbor",
                                   RHN = "Rhineland")
  }else{
    res_df$Time <- as.factor(df$Time)
  }
  print(identical(rownames(res_df), rownames(res_df))) 
  #res_df$Type <- res_df$Type %>% 
  #  recode_factor('in' = "Inside", 'out'="Outside")
  return(res_df)
}

ExtrProcDf(pro_root, df_its_root, "site") -> res_root
ExtrProcDf(pro_soil, df_its_soil, "site") -> res_soil

ExtrProcDf(pro_root_time, df_its_root_time, "time") -> res_root_time
ExtrProcDf(pro_soil_time, df_its_soil_time, "time") -> res_soil_time

# Kruskal Wallis
compare_means(Res_Err ~ FertStatus,
              res_root,
              method = "kruskal.test",
              p.adjust.method = "fdr") #  ns
compare_means(Res_Err ~ siteID,
              res_root,
              method = "kruskal.test",
              p.adjust.method = "fdr") # p.adj = 0.00000000045 ****

compare_means(Res_Err ~ FertStatus,
              res_soil,
              method = "kruskal.test",
              p.adjust.method = "fdr") # ns
compare_means(Res_Err ~ siteID,
              res_soil,
              method = "kruskal.test",
              p.adjust.method = "fdr") # p.adj = 0.0000002 ****

# Lux Arobor
compare_means(Res_Err ~ FertStatus,
              res_root_time,
              method = "kruskal.test",
              p.adjust.method = "fdr") #  p.adj = 0.003 **
compare_means(Res_Err ~ Time,
              res_root_time,
              method = "kruskal.test",
              p.adjust.method = "fdr") # p.adj = 0.054 ns

# Soil has much more samplig dates than root, so I test same dates as roots and all the dates
compare_means(Res_Err ~ FertStatus,
              res_soil_time,
              method = "kruskal.test",
              p.adjust.method = "fdr") #  p.adj = 0.025 *
compare_means(Res_Err ~ Time,
              res_soil_time,
              method = "kruskal.test",
              p.adjust.method = "fdr") # p.adj = 0.18 ns

res_soil_time %>%
  subset(
    Time %in% c(
      "2018-05-29",
      "2018-06-25",
      "2018-07-30",
      "2018-08-20",
      "2018-09-17",
      "2018-10-03"
    )
  ) %>%
  compare_means(
    data = .,
    Res_Err ~ FertStatus,
    method = "kruskal.test",
    p.adjust.method = "fdr"
  ) # p.adj = 0.046 *

res_soil_time %>%
  subset(
    Time %in% c(
      "2018-05-29",
      "2018-06-25",
      "2018-07-30",
      "2018-08-20",
      "2018-09-17",
      "2018-10-03"
    )
  ) %>%
  compare_means(
    data = .,
    Res_Err ~ Time,
    method = "kruskal.test",
    p.adjust.method = "fdr"
  ) # p.adj = 0.056 ns

str(res_soil_time)
as.factor(res_soil_time$Time)

# Now separating Fertilized vs Unfertilized plots ---------------------------------------------
res_root %>%
  subset(FertStatus %in% "Fert") %>%
  compare_means(
    data = .,
    Res_Err ~ siteID,
    method = "kruskal.test",
    p.adjust.method = "fdr"
  ) # p.adj =  0.00027 ***

res_root %>%
  subset(FertStatus %in% "Unfert") %>%
  compare_means(
    data = .,
    Res_Err ~ siteID,
    method = "kruskal.test",
    p.adjust.method = "fdr"
  ) # p.adj = 0.000009 ***


res_soil %>%
  subset(FertStatus %in% "Fert") %>%
  compare_means(
    data = .,
    Res_Err ~ siteID,
    method = "kruskal.test",
    p.adjust.method = "fdr"
  ) # p.adj =  0.0027 ***

res_soil %>%
  subset(FertStatus %in% "Unfert") %>%
  compare_means(
    data = .,
    Res_Err ~ siteID,
    method = "kruskal.test",
    p.adjust.method = "fdr"
  ) # p.adj = 0.000042 ***

# Lux Arbor
res_root_time %>%
  subset(FertStatus %in% "Fert") %>%
  compare_means(
    data = .,
    Res_Err ~ Time,
    method = "kruskal.test",
    p.adjust.method = "fdr"
  ) # 0.0011 **

res_root_time %>%
  subset(FertStatus %in% "Unfert") %>%
  compare_means(
    data = .,
    Res_Err ~ Time,
    method = "kruskal.test",
    p.adjust.method = "fdr"
  ) # 0.51


res_soil_time %>%
  subset(FertStatus %in% "Fert") %>%
  compare_means(
    data = .,
    Res_Err ~ Time,
    method = "kruskal.test",
    p.adjust.method = "fdr"
  ) # 0.49

res_soil_time %>%
  subset(FertStatus %in% "Unfert") %>%
  compare_means(
    data = .,
    Res_Err ~ Time,
    method = "kruskal.test",
    p.adjust.method = "fdr"
  ) # 0.51

# Using same dates as root samples
res_soil_time %>% 
  subset(FertStatus%in%"Fert") %>%
  subset(Time%in%c("2018-05-29","2018-06-25","2018-07-30","2018-08-20","2018-09-17","2018-10-03")) %>%
  compare_means(data = ., Res_Err ~ Time, method = "kruskal.test", p.adjust.method = "fdr") # p.adj = 0.51

res_soil_time %>% 
  subset(FertStatus%in%"Unfert") %>%
  subset(Time%in%c("2018-05-29","2018-06-25","2018-07-30","2018-08-20","2018-09-17","2018-10-03")) %>%
  compare_means(data = ., Res_Err ~ Time, method = "kruskal.test", p.adjust.method = "fdr") # p.adj = 0.77 ns

# calculating residuals means/median ----------------------------------------------------------
ExtProcLab <- function(dataframe, formula, method, Var){
  require(multcompView)
  compare_means(formula, data = dataframe, method = method, p.adjust.method = "fdr") -> test_CC 
  test_CC <- as.data.frame(test_CC)[,c(2,3,5)] # to change form p to p.adj do 4 to 5
  test_CC2 <- data.frame(test_CC[,2], test_CC[,1], test_CC[,3])
  colnames(test_CC2) <- c("group1", "group2", "p.adj") # change p to p.adj
  rbind(test_CC, test_CC2) -> test_all
  as.dist(xtabs(test_all[, 3] ~ (test_all[, 2] + test_all[, 1])), diag = TRUE) -> dist_CC
  data.frame(multcompLetters(dist_CC)['Letters']) -> res_CC
  res_CC$sample <- rownames(res_CC)
  if (Var == "site") {
    res_CC %>% slice(match(
      c("Escanaba", 
        "Hancock",
        "Lake City",
        "Lux Arbor",
        "Rhineland"),
      sample
    )) -> res_CC
    # } else {
    #   res_CC %>% slice(match(
    #     c(
    #       "2018-05-29",
    #       "2018-06-25",
    #       "2018-07-30",
    #       "2018-08-20",
    #       "2018-09-17",
    #       "2018-10-03"
    #     ),
    #     sample
    #   )) -> res_CC
  }
  as.character(res_CC$Letters) -> res_CC
  return(res_CC)
}


ExtProcLab(res_root, formula(Res_Err ~ siteID), "wilcox.test", "site") -> labels_res_root
labels_res_root

ExtProcLab(res_soil, formula(Res_Err ~ siteID), "wilcox.test", "site") -> labels_res_soil
labels_res_soil

# Fertilized vs unfertilized - only the significant groups! 
# All sites 
res_root %>%
  subset(FertStatus %in% "Fert") %>%
  ExtProcLab(dataframe = ., formula(Res_Err ~ siteID), "wilcox.test", "site") -> labels_root_fert
labels_root_fert

res_root %>%
  subset(FertStatus %in% "Unfert") %>%
  ExtProcLab(dataframe = ., formula(Res_Err ~ siteID), "wilcox.test", "site") -> labels_root_unfert
labels_root_unfert

res_soil %>%
  subset(FertStatus %in% "Fert") %>%
  ExtProcLab(dataframe = ., formula(Res_Err ~ siteID), "wilcox.test", "site") -> labels_soil_fert
labels_soil_fert

res_soil %>%
  subset(FertStatus %in% "Unfert") %>%
  ExtProcLab(dataframe = ., formula(Res_Err ~ siteID), "wilcox.test", "site") -> labels_soil_unfert
labels_soil_unfert

# Lux Arbor
res_root_time %>%
  subset(FertStatus %in% "Fert") %>%
  ExtProcLab(dataframe = , formula(Res_Err ~ Time), "wilcox.test", "time") -> labels_root_fert_time
labels_root_fert_time

res_root_time %>%
  subset(FertStatus %in% "Unfert") %>%
  ExtProcLab(dataframe = , formula(Res_Err ~ Time), "wilcox.test", "time") -> labels_root_unfert_time
labels_root_unfert_time


# TEST PROCRUSTES ERRORS within sites ---------------------------------------------------------
# test significance for Fertilization within site
ExtProcSite <- function(dataframe, formula, method, Var, set){
  require(multcompView)
  if (set == "site"){
    dataframe %>% 
      subset(siteID%in%Var) -> df_new
  }else{
    dataframe %>% 
      subset(Time%in%Var) -> df_new
  }
  print(dim(df_new))
  test_CC <- 
    compare_means(formula, data = df_new, method = method, p.adjust.method = "fdr")
  test_CC <- 
    as.data.frame(test_CC)[,c(2,3,5)] # to change form p to p.adj do 4 to 5
  test_CC2 <- 
    data.frame(test_CC[,2], test_CC[,1], test_CC[,3])
  colnames(test_CC2) <- c("group1", "group2", "p.adj") # change p to p.adj
  test_all <- 
    rbind(test_CC, test_CC2)
  print(test_all)
  dist_CC <-
    as.dist(xtabs(test_all[, 3] ~ (test_all[, 2] + test_all[, 1])), diag = TRUE)
  res_CC <-
    data.frame(multcompLetters(dist_CC)['Letters'])
  res_CC$sample <- rownames(res_CC)
  res_CC %>% 
    slice(match(c("Fert", "Unfert"), sample)) -> res_CC
  as.character(res_CC$Letters) -> res_CC
  return(res_CC)
}


# All sites
ExtProcSite(res_root, formula(Res_Err ~ FertStatus), "wilcox.test", "Escanaba", "site")
ExtProcSite(res_root, formula(Res_Err ~ FertStatus), "wilcox.test", "Hancock", "site")
ExtProcSite(res_root, formula(Res_Err ~ FertStatus), "wilcox.test", "Lake City", "site")
ExtProcSite(res_root, formula(Res_Err ~ FertStatus), "wilcox.test", "Lux Arbor", "site")
ExtProcSite(res_root, formula(Res_Err ~ FertStatus), "wilcox.test", "Rhineland", "site") # a, b

label_root <- c("ns", "ns", "ns", "ns", "ns", "ns", "ns", "ns", "ns", "**")
label_root <- c("a", "a", "a", "a", "a", "a", "a", "a", "a", "b")
letters_root_fert <- label_root[1:5]
letters_root_unfert <- label_root[6:10]

ExtProcSite(res_soil, formula(Res_Err ~ FertStatus), "wilcox.test", "Escanaba", "site")
ExtProcSite(res_soil, formula(Res_Err ~ FertStatus), "wilcox.test", "Hancock", "site")
ExtProcSite(res_soil, formula(Res_Err ~ FertStatus), "wilcox.test", "Lake City", "site")
ExtProcSite(res_soil, formula(Res_Err ~ FertStatus), "wilcox.test", "Lux Arbor", "site")
ExtProcSite(res_soil, formula(Res_Err ~ FertStatus), "wilcox.test", "Rhineland", "site")

label_soil <- c("a", "a", "a", "a", "a", "a", "a", "a", "a", "a")
letters_soil_fert <- label_soil[1:5]
letters_soil_unfert <- label_soil[6:10]


# Lux Arbor
ExtProcSite(res_root_time, formula(Res_Err ~ FertStatus), "wilcox.test", "2018-05-29", "time")
ExtProcSite(res_root_time, formula(Res_Err ~ FertStatus), "wilcox.test", "2018-06-25", "time")
ExtProcSite(res_root_time, formula(Res_Err ~ FertStatus), "wilcox.test", "2018-07-30", "time")
ExtProcSite(res_root_time, formula(Res_Err ~ FertStatus), "wilcox.test", "2018-08-20", "time")
ExtProcSite(res_root_time, formula(Res_Err ~ FertStatus), "wilcox.test", "2018-09-17", "time") 
ExtProcSite(res_root_time, formula(Res_Err ~ FertStatus), "wilcox.test", "2018-10-03", "time") 

label_root_time <- c("a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "b")
letters_root_time_fert <- label_root_time[1:6]
letters_root_time_unfert <- label_root_time[7:12]


# plotting functions for Residuals ------------------------------------------------------------
ProcrustResAver <- function(dataframe, my_labels=NULL, Var){
  #calculating where to put the letters
  max(dataframe[,"Res_Err"] + 0.1* max(dataframe[,"Res_Err"])) -> labels_y
  # plot
  ggplot(dataframe, aes(x = get(Var), 
                        y = Res_Err, 
                        color =  get(Var))) +
    geom_boxplot(outlier.colour="black", outlier.size=1.2, outlier.shape = 8, lwd=0.5, width = 0.7) +
    stat_summary(fun=mean, geom="point", shape=18, size=2, color="red", fill="red") +
    stat_summary(geom = 'text', label = my_labels, fun= max, aes(y = labels_y), size=3.5, color="black") +
    geom_jitter(size=1.5, alpha=0.5, aes(shape=FertStatus)) +
    #facet_grid(~Type, scales = "free_x")+
    theme_classic() +
    expand_limits(y = 0) +
    scale_shape_manual(values = c(19, 2), labels=c("Fert","Unfert")) +
    theme(strip.text.x = element_text(size = 8)) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 8)) +
    theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1, vjust = 1)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
    grids(linetype = "dashed") +
    theme(legend.position="none") -> rich_plot 
  return(rich_plot)
}


ProcrustResAver(res_root, labels_res_root, "siteID") + 
  scale_color_brewer(palette = "Set1",name=NULL,
                     labels=c("Escanaba","Hancock","Lake City","Lux Arbor","Rhineland")) +
  labs(title = "Root sample procrustes errors", x=NULL, y="Procrustes residual error") 


# Plotting and adding letters for siteID/FertStatus comparisons
# Generating letters

labels_root_fert
labels_root_unfert
labels_soil_fert
labels_soil_unfert
labels_root_fert_time
labels_root_unfert_time

letters_root_fert
letters_root_unfert
letters_soil_fert
letters_soil_unfert
letters_root_time_fert 
letters_root_time_unfert

com_letters_root_fert <- paste(labels_root_fert, letters_root_fert, sep = "/")
com_letters_root_unfert <- paste(labels_root_unfert, letters_root_unfert, sep = "/")

com_letters_soil_fert <- paste(labels_soil_fert, letters_soil_fert, sep = "/")
com_letters_soil_unfert <- paste(labels_soil_unfert, letters_soil_unfert, sep = "/")

com_letters_root_time_fert <- paste(labels_root_fert_time, letters_root_time_fert, sep = "/")
com_letters_root_time_unfert <- paste(labels_root_unfert_time, letters_root_time_unfert, sep = "/")

# All sites
res_root %>%
  subset(FertStatus %in% "Fert") %>%
  ProcrustResAver(dataframe = ., com_letters_root_fert, "siteID") +
  scale_color_brewer(
    palette = "Set1",
    name = NULL,
    labels = c("Escanaba", "Hancock", "Lake City", "Lux Arbor", "Rhineland")
  ) +
  labs(title = "Root sample procrustes errors", x = NULL, y = "Procrustes residual error")

res_root %>%
  subset(FertStatus %in% "Unfert") %>%
  ProcrustResAver(dataframe = ., com_letters_root_unfert, "siteID") +
  scale_color_brewer(
    palette = "Set1",
    name = NULL,
    labels = c("Escanaba", "Hancock", "Lake City", "Lux Arbor", "Rhineland")
  ) +
  labs(title = "Root sample procrustes errors", x = NULL, y = "Procrustes residual error")


res_soil %>%
  subset(FertStatus %in% "Fert") %>%
  ProcrustResAver(dataframe = ., com_letters_soil_fert, "siteID") +
  scale_color_brewer(
    palette = "Set1",
    name = NULL,
    labels = c("Escanaba", "Hancock", "Lake City", "Lux Arbor", "Rhineland")
  ) +
  labs(title = "Root sample procrustes errors", x = NULL, y = "Procrustes residual error")

res_soil %>%
  subset(FertStatus %in% "Unfert") %>%
  ProcrustResAver(dataframe = ., com_letters_soil_unfert, "siteID") +
  scale_color_brewer(
    palette = "Set1",
    name = NULL,
    labels = c("Escanaba", "Hancock", "Lake City", "Lux Arbor", "Rhineland")
  ) +
  labs(title = "Root sample procrustes errors", x = NULL, y = "Procrustes residual error")

# Plotting and adding letters for Time/FertStatus comparisons
# Lux Arbor 
res_root_time %>%
  subset(FertStatus %in% "Fert") %>%
  ProcrustResAver(dataframe = ., com_letters_root_time_fert, "Time") +
  scale_color_manual(values = root_time_color) +
  labs(title = "Lux Arbor\nRoot sample procrustes errors", x = NULL, y =
         "Residual error")

res_root_time %>%
  subset(FertStatus %in% "Unfert") %>%
  ProcrustResAver(dataframe = ., com_letters_root_time_unfert, "Time") +
  scale_color_manual(values = root_time_color) +
  labs(title = "Lux Arbor\nRoot sample procrustes errors", x = NULL, y =
         "Residual error")

# *** FIGURE 2 Procrustes Errors --------------------------------------------------------------
ggarrange(
  res_root %>%
    subset(FertStatus %in% "Fert") %>%
    ProcrustResAver(dataframe = ., com_letters_root_fert, "siteID") +
    scale_color_brewer(
      palette = "Set1",
      name = NULL,
      labels = c("Escanaba", "Hancock", "Lake City", "Lux Arbor", "Rhineland")
    ) +
    labs(title = "All sites\nfertilized root samples", x = NULL, y =
           "Residual error"),
  res_soil %>%
    subset(FertStatus %in% "Fert") %>%
    ProcrustResAver(dataframe = ., com_letters_soil_fert, "siteID") +
    scale_color_brewer(
      palette = "Set1",
      name = NULL,
      labels = c("Escanaba", "Hancock", "Lake City", "Lux Arbor", "Rhineland")
    ) +
    labs(title = "All sites\nfertilized soil samples", x = NULL, y =
           "Residual error"),
  res_root_time %>%
    subset(FertStatus %in% "Fert") %>%
    ProcrustResAver(dataframe = ., com_letters_root_time_fert, "Time") +
    scale_color_manual(values = root_time_color) +
    labs(title = "Lux Arbor\nfertilized root samples", x = NULL, y =
           "Residual error"),
  
  res_root %>%
    subset(FertStatus %in% "Unfert") %>%
    ProcrustResAver(dataframe = ., com_letters_root_unfert, "siteID") +
    scale_color_brewer(
      palette = "Set1",
      name = NULL,
      labels = c("Escanaba", "Hancock", "Lake City", "Lux Arbor", "Rhineland")
    ) +
    labs(title = "All sites\nunfertilized root samples", x = NULL, y =
           "Residual error"),
  res_soil %>%
    subset(FertStatus %in% "Unfert") %>%
    ProcrustResAver(dataframe = ., com_letters_soil_unfert, "siteID") +
    scale_color_brewer(
      palette = "Set1",
      name = NULL,
      labels = c("Escanaba", "Hancock", "Lake City", "Lux Arbor", "Rhineland")
    ) +
    labs(title = "All sites\nunfertilized soil samples", x = NULL, y =
           "Residual error"),
  res_root_time %>%
    subset(FertStatus %in% "Unfert") %>%
    ProcrustResAver(dataframe = ., com_letters_root_time_unfert, "Time") +
    scale_color_manual(values = root_time_color) +
    labs(title = "Lux Arbor\nunfertilized root samples", x = NULL, y =
           "Residual error"),
  labels = c("A", "B", "C", "D", "E", "F"),
  align = "hv",
  ncol = 3,
  nrow = 2
) -> Fig_X_procrust_sign

Fig_X_procrust_sign

library(gridExtra)
title = text_grob("Procrustes Errors", size = 12, face = 2)
grid.arrange(Fig_X_procrust_sign, top=title)


save(Fig_X_procrust_sign, 
     res_root, res_soil, res_root_time,
     com_letters_root_fert, com_letters_root_unfert,
     com_letters_soil_fert, com_letters_soil_unfert,
     com_letters_root_time_fert, com_letters_root_time_unfert,
     file = "Fig_S2_residuals.Rdata")



