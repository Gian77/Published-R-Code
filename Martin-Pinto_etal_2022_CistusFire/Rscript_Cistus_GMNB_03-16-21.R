# ************ DATA ANALYSIS ***************************************** --------------------------------------------
# Project name: Fungi and Bacteria on Cistus
# Manuscript:   
# Authors:      Martin Pinto P, Benucci GM ...,
# Affiliation:  Universidad de Valladolid, Michigan State University, ...
# Journal:      
# Date:         March 16, 2021
#
# R Code:       Benucci GMN
# Please contact <gian.benucci@gmail.com> before using this code. 
# ******************************************************************** --------------------------------------------

# starting packages -----------------------------------------------------------------------------------------------
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(ape)
library(dplyr)
library(ggpubr)

# PALETTE ----------------------------------------------------------------------------------------------------------
palette3 = c("#058ED9", "#599861", "#FF934F")
pie(rep(1, length(palette3)),
    labels = sprintf("%d (%s)",
                     seq_along(palette3), palette3),
    col = palette3)

# IMPORTING DATASETS ----------------------------------------------------------------------------------------------
dir()

# A) ITS UPARSE ---------------------------------------------------------------------------------------------------
read.csv(
  file = "otu_table_tax_fun.csv",
  sep = ",",
  header = TRUE,
  row.names = 1
) -> df_fungi
str(df_fungi)
head(df_fungi)
dim(df_fungi)
rownames(df_fungi) <- paste("F", rownames(df_fungi), sep = "")

# B) Bacteria -----------------------------------------------------------------------------------------------------
read.csv(
  file = "otu_table_tax_bact.csv",
  sep = ",",
  header = TRUE,
  row.names = 1
) -> df_bac
str(df_bac)
head(df_bac)
rownames(df_bac) <- paste("B", rownames(df_bac), sep = "")

read.csv(
  file = "metadata.csv",
  sep = ",",
  header = TRUE,
  row.names = 1
) -> metadata
str(metadata)
head(metadata)

# IMPORTING FunGuild AND Mycelium DATA ----------------------------------------------------------------------------
read.csv(
  "mycelium_quantity.csv",
  header = TRUE,
  row.names = 1,
  sep = ","
) -> myceliumDNA
myceliumDNA

read.csv(
  "funguild_input.csv",
  header = TRUE,
  row.names = 1,
  sep = ","
) -> funGuild
rownames(funGuild) <- paste("F", rownames(funGuild), sep = "")
funGuild_otu <- rownames(funGuild)
funGuild_otu

# IMPORTING NEW TAXONOMIES ----------------------------------------------------------------------------------------
library(tidyr)

# A) Fungi --------------------------------------------------------------------------------------------------------
# Importing taxonomies at 0.6 confidence --------------------------------------------------------------------------
taxonomy_ITS06 <-
  read.delim(
    "taxonomy_assignments_euk_06/consensus_taxonomy.txt",
    header = TRUE,
    row.names = 1,
    sep = "\t"
  )

head(taxonomy_ITS06)
taxonomy_ITS06[1:100,]

# Importing taxonomies at 0.8 confidence --------------------------------------------------------------------------
taxonomy_ITS08 <-
  read.delim(
    "taxonomy_assignments_euk_08/consensus_taxonomy.txt",
    header = TRUE,
    row.names = 1,
    sep = "\t"
  )

head(taxonomy_ITS08)
taxonomy_ITS08[1:100,]

# check for identical ordering
identical(rownames(taxonomy_ITS06), rownames(taxonomy_ITS08))
taxonomy_ITS08$Kingdom_06 <- taxonomy_ITS06$Kingdom
head(taxonomy_ITS08)
taxonomy_ITS08[1:100,]
dim(taxonomy_ITS08)
levels(taxonomy_ITS08$Kingdom)

# how many unclassified OTUs in the two taxonomies?
nrow(taxonomy_ITS06[taxonomy_ITS06$Kingdom != "Fungi", ])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom != "Fungi", ])

# Non-target taxa
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom == "Alveolata", ])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom == "Glaucocystoplantae", ])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom == "Heterolobosa", ])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom == "Rhizaria", ])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom == "Metazoa", ])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom == "Viridiplantae", ])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom == "Stramenopila", ])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom == "", ])

nrow(taxonomy_ITS08[taxonomy_ITS06$Kingdom == "Fungi", ])

# removing non-fungal OTUs
subset(taxonomy_ITS08, taxonomy_ITS08$Kingdom_06 == "Fungi") -> taxonomy_ITS08_filt
dim(taxonomy_ITS08_filt)
str(taxonomy_ITS08_filt)

# checking filtering results
nrow(taxonomy_ITS08_filt[taxonomy_ITS08_filt$Kingdom != "Fungi",])
nrow(taxonomy_ITS08_filt[taxonomy_ITS08_filt$Kingdom_06 != "Fungi",])

# classify all OTUs at 60% identity to Fungi
taxonomy_ITS08_filt$Kingdom <- taxonomy_ITS08_filt$Kingdom_06

# formatting ranks
head(taxonomy_ITS08_filt)
taxonomy_ITS08_filt <- taxonomy_ITS08_filt[, 1:7]
taxonomy_ITS08_filt$OTU_ID <- rownames(taxonomy_ITS08_filt)

taxonomy_ITS08_filt <- taxonomy_ITS08_filt[, c("Kingdom",
                                               "Phylum",
                                               "Class",
                                               "Order",
                                               "Family",
                                               "Genus",
                                               "Species",
                                               "OTU_ID")]

rownames(taxonomy_ITS08_filt) <-
  paste("F", rownames(taxonomy_ITS08_filt), sep = "")


# MATCHING TAXONOMY TO OTUS ---------------------------------------------------------------------------------------
dim(df_fungi)
dim(taxonomy_ITS08_filt)

# match taxonomy to OTUs
taxonomy_ITS08_filt[rownames(df_fungi), ] -> tax_fungi
dim(tax_fungi)
tax_fungi %>% drop_na() -> tax_fungi

# 2 taxa are not in the dataframe so I am removing them
df_fungi[rownames(tax_fungi), ] -> df_fungi
dim(df_fungi)
sort(colSums(df_fungi[, 1:27]), decreasing = TRUE)

# EXTRACTING ECM OTU TABLE ----------------------------------------------------------------------------------------
length(funGuild_otu)
funGuild_otu

#match otus to ecm taxa
df_fungi[funGuild_otu, ] -> df_ecm
df_ecm

df_ecm %>% drop_na() -> df_ecm
sort(colSums(df_ecm[, 1:27]), decreasing = TRUE)
dim(df_ecm)

# getting ecm taxa
taxonomy_ITS08_filt[rownames(df_ecm), ] -> tax_ecm
dim(tax_ecm)
tax_ecm

# Reformat taxonomy -----------------------------------------------------------------------------------------------
source("../R_functions/ReformatTaxonomy.R")

tax_ecm$Species <-
  ifelse(grepl("sp", tax_ecm$Species) == FALSE,
         paste(tax_ecm$Species),"")

tax_ecm$OTU_ID <- as.factor(tax_ecm$OTU_ID)
tax_ecm$Species <- as.factor(tax_ecm$Species)

ReformatTaxonomy <- function(taxa) {
  taxa[taxa == ""] <- NA
  taxa$Genus <- as.character(taxa$Genus)
  taxa[which(is.na(taxa$Genus) == FALSE),]$Genus <- paste(taxa$Genus[is.na(taxa$Genus) ==
                                                                       FALSE], "sp.", sep = " ")
  taxa <- taxa[c(8, 1, 2, 3, 4, 5, 6, 7)]
  taxa[] = lapply(taxa,
                  blank2na,
                  na.strings = c('', 'NA', 'na', 'N/A', 'n/a', 'NaN', 'nan'))
  lastValue <- function(x)
    tail(x[!is.na(x)], 1)
  last_taxons <- apply(taxa[, 1:8], 1, lastValue)
  taxa$BestMatch <- last_taxons
  taxa[, "BestMatch"] <- gsub("_", " ", taxa[, "BestMatch"])
  taxa$Taxonomy <- paste(taxa$OTU_ID, taxa$BestMatch, sep = "-")
  taxa[, "Genus"] <- gsub(" sp.", "", taxa[, "Genus"])
  return(taxa)
}


ReformatTaxonomy(as.data.frame(as.matrix(tax_ecm))) -> tax_ecm_filt
head(tax_ecm_filt)

ReformatTaxonomy(as.data.frame(as.matrix(tax_fungi))) -> tax_fungi
head(tax_fungi)


# LIBRAY SIZES and ------------------------------------------------------------------------------------------------
library(vegan)
library(ggpubr)

colSums(df_fungi[,1:27]) -> metadata$Reads_fun
colSums(df_bac[,1:27]) -> metadata$Reads_bac

PlotDist <- function(df, Size) {
  ggplot(df, aes(x = get(Size))) +
    labs(x = "Read number", y = "Sample number") +
    geom_histogram(binwidth = 1000,
                   colour = "grey80",
                   fill = "grey80") +
    geom_vline(
      aes(xintercept = mean(get(Size), na.rm = T)),
      color = "red",
      linetype = "dashed",
      size = 0.8
    ) +
    expand_limits(x = 0, y = 0) +
        theme_classic() +
        theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
        theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) + 
        theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
        theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
        theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
        theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
        theme(legend.position="bottom") -> plot_dist
      return(plot_dist)  
}

PlotDist(metadata, "Reads_fun") + labs(title="Fungi\nDistribution of sample libraries")
PlotDist(metadata, "Reads_bac") + labs(title="Bacteria\nDistribution of sample libraries")

# *** FIGURE S1 ----------------------------------------------------------------------------------------------------
ggarrange(
  PlotDist(metadata, "Reads_fun") + labs(title = "Fungi\nDistribution of sample libraries"),
  PlotDist(metadata, "Reads_bac") + labs(title = "Bacteria\nDistribution of sample libraries"),
  labels = c("A", "B"),
  widths = c(1, 1),
  align = "hv" ,
  common.legend = TRUE,
  ncol = 2,
  nrow = 1,
  legend = c("bottom")
) -> Fig_1

Fig_1


# RAREFACTION CURVES ----------------------------------------------------------------------------------------------
min(colSums(df_fungi[, 1:27])) -> min_fungi
min(colSums(df_bac[, 1:27])) -> min_bact

# creating and empty character
curve_colors_Site <- vector(mode = "character", length = 27)
curve_colors_Site[metadata$Site == "burnt"] <- palette3[1]
curve_colors_Site[metadata$Site == "clearing"] <- palette3[2]
curve_colors_Site[metadata$Site == "old-growth"] <- palette3[3]
curve_colors_Site

# *** FIGURE S2 ---------------------------------------------------------------------------------------------------
par(
  mar = c(3.5, 3.5, 2, 0) + 0.1,
  mgp = c(2.4, 0.8, 0),
  las = 0,
  mfrow = c(1, 2)
)

rarecurve(
  t(df_fungi[, 1:27]),
  step = 100,
  col = curve_colors_Site,
  label = FALSE,
  cex.axis=0.8,
  cex.lab=0.8,
  cex.main=1,
  main = "Fungi",
  ylab = "Number of OTUs",
  xlab = "Number of DNA reads")
abline(v = min_fungi, col="red", lwd=3, lty=2)
rarecurve(
  t(df_bac[, 1:27]),
  step = 100,
  col = curve_colors_Site,
  label = FALSE,
  cex.axis=0.8,
  cex.lab=0.8,
  cex.main=1,
  main = "Bacteria",
  ylab = "",
  xlab = "Number of DNA reads")
abline(v = min_bact, col="red", lwd=3, lty=2)
legend(
  "bottomright",
  legend = c("burnt", "clearing", "old-growth"),
  col = palette3,
  lty = 2,
  cex = 0.8,
  box.lty = 0,
  lwd = 3,
  bty = "n"
) 
dev.off()

# Phyloseq objects ---------------------------------------------------------------------------------------
head(df_fungi)

physeq_fungi <-
  phyloseq(
    otu_table(df_fungi[1:27], taxa_are_rows = TRUE),
    sample_data(metadata),
    tax_table(as.matrix(df_fungi %>% select(Taxon))))

physeq_fungi
head(physeq_fungi@tax_table)
physeq_fungi@sam_data
physeq_fungi@otu_table

head(df_bac)

physeq_bac <-
  phyloseq(
    otu_table(df_bac[1:27], taxa_are_rows = TRUE),
    sample_data(metadata),
    tax_table(as.matrix(df_bac %>% select(Taxon))))

physeq_bac
head(physeq_bac@tax_table)
physeq_bac@sam_data
physeq_bac@otu_table

# Rarefy even depth ----------------------------------------------------------------------------------------
#fungi
set.seed(2020)
sort(colSums(otu_table(physeq_fungi)), decreasing = TRUE)
physeq_fungi_ev <-
  rarefy_even_depth(physeq_fungi, sample.size = 19678, 
                                    rngseed = FALSE, replace = TRUE, 
                                    trimOTUs = TRUE, verbose = TRUE) 
otu_table(physeq_fungi_ev) <- 
  otu_table(physeq_fungi_ev)[which(rowSums(otu_table(physeq_fungi_ev)) >= 1),]
physeq_fungi_ev
colSums(otu_table(physeq_fungi_ev))
any(taxa_sums(physeq_fungi_ev) == 0)

sort(colSums(otu_table(physeq_bac)), decreasing = TRUE)
physeq_bac_ev <-
  rarefy_even_depth(physeq_bac, sample.size = 2337, 
                    rngseed = FALSE, replace = TRUE, 
                    trimOTUs = TRUE, verbose = TRUE) 
otu_table(physeq_bac_ev) <- 
  otu_table(physeq_bac_ev)[which(rowSums(otu_table(physeq_bac_ev)) >= 1),]
physeq_bac_ev
colSums(otu_table(physeq_bac_ev))
any(taxa_sums(physeq_bac_ev) == 0)


otu_fungi_ev <- as.data.frame(otu_table(physeq_fungi_ev))
otu_bact_ev <- as.data.frame(otu_table(physeq_bac_ev))

# ************************************************************************************-----------------------------
# BETA DIVERSITY --------------------------------------------------------------------------------------------------
# Non-metric Multidimensional Scaling
fungi_NMDS = metaMDS(
  t(df_fungi[, 1:27]),
  k = 2,
  trymax = 200,
  distance = "bray",
  weakties = TRUE
)
fungi_NMDS
stressplot(fungi_NMDS)

identical(rownames(as.data.frame(fungi_NMDS$points)), rownames(metadata))
cbind(metadata, as.data.frame(fungi_NMDS$points)) -> metadata_fungi
metadata_fungi


bact_NMDS = metaMDS(
  t(df_bac[, 1:27]),
  k = 2,
  trymax = 200,
  distance = "bray",
  weakties = TRUE
)
bact_NMDS
stressplot(bact_NMDS)

identical(rownames(as.data.frame(bact_NMDS$points)), rownames(metadata))
cbind(metadata, as.data.frame(bact_NMDS$points)) -> metadata_bact
metadata_bact

# Rarefied data 
fungi_NMDS_ev = metaMDS(
  t(otu_fungi_ev),
  k = 2,
  trymax = 200,
  distance = "bray",
  weakties = TRUE
)
fungi_NMDS_ev
stressplot(fungi_NMDS_ev)

identical(rownames(as.data.frame(fungi_NMDS_ev$points)), rownames(metadata))
cbind(metadata_fungi, 
      as.data.frame(fungi_NMDS_ev$points) %>% 
        rename(axis1="MDS1", axis2="MDS2")) -> metadata_fungi
metadata_fungi

bact_NMDS_ev = metaMDS(
  t(otu_bact_ev),
  k = 2,
  trymax = 200,
  distance = "bray",
  weakties = TRUE
)
bact_NMDS_ev
stressplot(bact_NMDS_ev)

identical(rownames(as.data.frame(bact_NMDS_ev$points)), rownames(metadata))
cbind(metadata_bact, 
      as.data.frame(bact_NMDS_ev$points) %>% 
        rename(axis1="MDS1", axis2="MDS2")) -> metadata_bact
metadata_bact


# plot ordination 
PlotOrdin <- function(dataframe, nmds) {
  ggplot(dataframe, aes(NMDS1, NMDS2, color = Site, shape = Treatment)) +
    geom_point(size = 1.5, stroke = 1) +
    scale_colour_manual("Site", values = palette3) +
    scale_shape_manual(name = "Treatment", values = c(0, 1, 2, 5)) +
    theme_classic() +
          theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
          theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
          theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
          theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
          theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
          theme(legend.title = element_text(size = 8, face = "bold"), 
                legend.text = element_text(size = 7)) +
          theme(legend.position="bottom") +
    annotate("text", x=-Inf, y = Inf,
             label=paste("italic(stress) ==", round(nmds$stress, 3)), parse = TRUE, size=2.5, vjust=1.5, hjust=-0.1) +
    guides(
      color = guide_legend(nrow = 1),
      shape = guide_legend(nrow = 1)) +
    grids(linetype = "dashed") -> plot_ord
  return(plot_ord)
}


PlotOrdin(metadata_fungi, fungi_NMDS) + labs(title = "Fungi")
PlotOrdin(metadata_bact, bact_NMDS) + labs(title = "Prokaryote")


# plot ordination 
PlotOrdin <- function(dataframe, nmds) {
  ggplot(dataframe, aes(axis1, axis2, color = Site, shape = Treatment)) +
    geom_point(size = 1.5, stroke = 1) +
    scale_colour_manual("Site", values = palette3) +
    scale_shape_manual(name = "Treatment", values = c(0, 1, 2, 5)) +
    theme_classic() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
    theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 8, face = "bold"), 
          legend.text = element_text(size = 7)) +
    theme(legend.position="bottom") +
    annotate("text", x=-Inf, y = Inf,
             label=paste("italic(stress) ==", round(nmds$stress, 3)), parse = TRUE, size=2.5, vjust=1.5, hjust=-0.1) +
    guides(
      color = guide_legend(nrow = 1),
      shape = guide_legend(nrow = 1)) +
    grids(linetype = "dashed") -> plot_ord
  return(plot_ord)
}


PlotOrdin(metadata_fungi, fungi_NMDS_ev) + labs(title = "Fungi")
PlotOrdin(metadata_bact, bact_NMDS_ev) + labs(title = "Prokaryote")


# Check for differences in rarefied vs non-rarefied
pro_fungi_fungiEV <- protest(fungi_NMDS,
                             fungi_NMDS_ev,
                             scores = "sites",
                             permutations = how(nperm = 9999))
pro_fungi_fungiEV
pro_fungi_fungiEV$ss



# PERMANOVA -------------------------------------------------------------------------------------------------------
adonis(t(df_fungi[, 1:27]) ~ Site * Treatment,
       metadata,
       method = "bray",
       permutations = 9999) -> adonis_fungi
adonis_fungi$aov.tab
# Site            2    1.3825 0.69124  4.5786 0.27988 0.0001 ***
round(p.adjust(adonis_fungi$aov.tab$`Pr(>F)`, "bonferroni"), 4)

adonis(t(df_fungi[, 1:27]) ~ Treatment * Site,
       metadata,
       method = "bray",
       permutations = 9999) -> adonis_fungi
adonis_fungi$aov.tab
#Site            2    1.3262 0.66312  4.3923 0.26850 0.0001 ***
round(p.adjust(adonis_fungi$aov.tab$`Pr(>F)`, "bonferroni"), 4)

adonis2(t(df_fungi[, 1:27]) ~ Treatment,
       metadata,
       method = "bray",
       permutations = 9999, 
       strata = metadata$Site) -> adonis_fungi_st
adonis_fungi_st
 


adonis(t(df_bac[, 1:27]) ~ Site * Treatment,
       metadata,
       method = "bray",
       permutations = 9999) -> adonis_bact
adonis_bact$aov.tab
#Site            2   0.55186 0.275928  3.5904 0.23631 0.0001 ***
round(p.adjust(adonis_bact$aov.tab$`Pr(>F)`, "bonferroni"), 4)

adonis(t(df_bac[, 1:27]) ~  Treatment * Site,
       metadata,
       method = "bray",
       permutations = 9999) -> adonis_bact
adonis_bact$aov.tab
#Site            2   0.48731 0.243657  3.1705 0.20867 0.0001 ***
round(p.adjust(adonis_bact$aov.tab$`Pr(>F)`, "bonferroni"), 4)

adonis(t(df_bac[, 1:27]) ~  Treatment,
       metadata,
       method = "bray",
       permutations = 9999, 
       strata = metadata$Site) -> adonis_bact_st
adonis_bact_st$aov.tab


# rarefied data 
adonis(t(otu_fungi_ev) ~ Site * Treatment,
       metadata,
       method = "bray",
       permutations = 9999) -> adonis_fungi_ev
adonis_fungi_ev$aov.tab

adonis(t(otu_fungi_ev) ~ Treatment * Site,
       metadata,
       method = "bray",
       permutations = 9999) -> adonis_fungi_ev
adonis_fungi_ev$aov.tab

adonis2(t(otu_fungi_ev) ~ Treatment,
        metadata,
        method = "bray",
        permutations = 9999, 
        strata = metadata$Site) -> adonis_fungi_st
adonis_fungi_st


adonis(t(otu_bact_ev) ~ Site * Treatment,
       metadata,
       method = "bray",
       permutations = 9999) -> adonis_bact_ev
adonis_bact_ev$aov.tab

adonis(t(otu_bact_ev) ~ Treatment * Site,
       metadata,
       method = "bray",
       permutations = 9999) -> adonis_bact_ev
adonis_bact_ev$aov.tab

adonis(t(otu_bact_ev) ~  Treatment,
       metadata,
       method = "bray",
       permutations = 9999, 
       strata = metadata$Site) -> adonis_bact_st
adonis_bact_st$aov.tab

# *** FIGURE 1 ordination -----------------------------------------------------------------------------------------
ggarrange(PlotOrdin(metadata_fungi, fungi_NMDS) + 
            labs(title = "Fungi") +
          annotate("text", -Inf, Inf, 
                   label = expression(paste("Site ", italic(R) ^ 2,"= 28.0%***"),
                                      parse=TRUE), size = 2.5, hjust = -0.1, vjust = 2.5),
          PlotOrdin(metadata_bact, bact_NMDS) + 
            labs(title = "Prokaryote") +
            annotate("text", -Inf, Inf, 
                     label = expression(paste("Site ", italic(R) ^ 2,"= 23.6%***"),
                                        parse=TRUE), size = 2.5, hjust = -0.1, vjust = 2.5),
          labels = c("A","B"),
          widths = c(1, 1),
          align = "hv" ,
          common.legend = TRUE,
          ncol = 2, nrow = 1, 
          legend = c("bottom")) -> Fig_1_beta

Fig_1_beta


# MULTIVARIATE DISPERSION -----------------------------------------------------------------------------------------

# Test the homogenity of group variances
table(metadata_fungi$Site)
permdisp_fungi_site <-
  betadisper(vegdist(t(df_fungi[, 1:27]), method = "bray"), metadata_fungi$Site)
anova(permdisp_fungi_site, permutations = 9999)
permutest(permdisp_fungi_site,
          permutations = 9999,
          pairwise = T) -> dist_fungi_site
round(p.adjust(
  anova(permdisp_fungi_site, permutations = 9999)$`Pr(>F)`,
  "bonferroni"
), 4)

boxplot(permdisp_fungi_site)

table(metadata_fungi$Treatment)
permdisp_fungi_treat <-
  betadisper(vegdist(t(df_fungi[, 1:27]), method = "bray"), metadata_fungi$Treatment)
anova(permdisp_fungi_treat, permutations = 9999)
permutest(permdisp_fungi_treat,
          permutations = 9999,
          pairwise = T) -> dist_fungi_treat
round(p.adjust(
  anova(permdisp_fungi_treat, permutations = 999)$`Pr(>F)`,
  "bonferroni"
), 4)

boxplot(permdisp_fungi_treat)

permdisp_bact_site <-
  betadisper(vegdist(t(df_bac[, 1:27]), method = "bray"), metadata_bact$Site)
anova(permdisp_bact_site, permutations = 9999)
permutest(permdisp_bact_site,
          permutations = 9999,
          pairwise = T) -> dist_bact_site
round(p.adjust(
  anova(permdisp_bact_site, permutations = 9999)$`Pr(>F)`,
  "bonferroni"
), 4) # p = 0.0031

boxplot(permdisp_bact_site)

permdisp_bact_treat <-
  betadisper(vegdist(t(df_bac[, 1:27]), method = "bray"), metadata_bact$Treatment)
anova(permdisp_bact_treat, permutations = 9999)
permutest(permdisp_bact_treat,
          permutations = 9999,
          pairwise = T) -> dist_bact_treat
round(p.adjust(
  anova(permdisp_bact_treat, permutations = 999)$`Pr(>F)`,
  "bonferroni"
), 4)

boxplot(permdisp_bact_treat)

library(multcompView)

# Significant differences - letters -------------------------------------------------------------------------------
PermDispLabel <- function(meta){
meta$Site %>% 
  recode_factor("burnt" = "burnt", "clearing"="clearing", "old-growth"="old_growth") -> meta$Site
permdisp_bact_site <- betadisper(vegdist(t(df_bac[,1:27]), method="bray"), meta$Site) 
anova(permdisp_bact_site, permutations = 9999)
permutest(permdisp_bact_site, permutations = 9999, pairwise = T) -> dist_bact_site
round(p.adjust(anova(permdisp_bact_site, permutations = 9999)$`Pr(>F)`, "bonferroni"), 4) # p = 0.0031 

data.frame(multcompLetters(p.adjust(dist_bact_site$pairwise$observed
                                    , method="bonferroni"))['Letters']) -> label_bact_site
return(label_bact_site)
}

PermDispLabel(metadata_bact) -> label_bact_site
label_bact_site


# plot multivariate dispersion 
PlotBetadisper <- function(betadisp, value, palette, my_labels, metadata){
  # creating a label and dataframe
  max(betadisp$distances + 0.1 * max(betadisp$distances)) -> labels_y
  #labels_y = 0.9
  data.frame(betadisp$group, betadisp$distances) -> df
  colnames(df) <- c("value", "distance")
  if (identical(rownames(df), rownames(metadata)) == TRUE) {
    df$Site <- metadata$Site
    df$Treatment <- metadata$Treatment
    # plotting
    ggplot(df, aes(x=value, y=distance, color=value)) +
      geom_jitter(size=1.2, alpha=0.5, aes(shape=Treatment), colour="black") +
      geom_boxplot(color="black", outlier.colour="black", outlier.shape = 8, alpha=0.6, lwd = 0.7) +
      stat_summary(fun=mean, geom="point", shape=18, size=1.9, color="red", fill="red") +
      stat_summary(geom = 'text', label = my_labels, fun= max, aes(y = labels_y), size=3, color="black") +
        scale_shape_manual("Treatment", values = c(0,1,2,5)) +
          theme_classic() +
          theme(strip.text.x = element_text(size = 9, face = "bold")) +
          theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
          theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1, vjust = 1)) +
          theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
          theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
          theme(axis.title.y =  element_text(angle = 90, size = 8, face = "bold")) + 
          #theme(axis.title.x = element_blank()) +
      grids(linetype = "dashed") +
          theme(legend.position="none") -> betaplot
    return(betaplot)
  }else{
    stop("Error: dataframe and metadata rownames are not matching!")
  }
}

# Non-significant permanova not included in the plot
PlotBetadisper(
  permdisp_bact_site,
  "Site",
  palette3,
  as.character(label_bact_site$Letters),
  metadata_bact) +
  labs(title = "ITS\n(Site)", y = "Distance to centroids", x = "Site")

PlotBetadisper(permdisp_fungi_treat,
               "Treatment",
               palette3,
               my_labels = "",
               metadata_bact) +
  labs(title = "ITS\n(Treatment)", y = "Distance to centroids", x = "Site")


# *** FIGURE 1 dispersion -----------------------------------------------------------------------------------------
ggarrange(PlotBetadisper(permdisp_fungi_site, "Site", palette3, my_labels = "", metadata_fungi) + 
            labs(title = "Fungi\n(Site)", y ="Distance to centroids", x= "Site"),
          PlotBetadisper(permdisp_fungi_treat, "Treatment", palette3, my_labels = "", metadata_fungi) + 
            labs(title = "Fungi\n(Treatment)", y =NULL, x= "Treatment"),
          PlotBetadisper(permdisp_bact_site, "Site", palette3, as.character(label_bact_site$Letters), metadata_bact) + 
            labs(title = "Bacbteria\n(Site)", y =NULL, x= "Site"),
          PlotBetadisper(permdisp_fungi_treat, "Treatment", palette3, my_labels = "", metadata_bact) + 
            labs(title = "Bacteria\n(Treatment)", y =NULL, x= "Treatment"),
          widths = c(1,1,1,1),
          labels = c("D","B","C","D"),
          align = "hv",
          ncol = 4, 
          nrow = 1) -> betadisp_plot

betadisp_plot

# *** FIGURE 1 COMPLETE -------------------------------------------------------------------------------------------
ggarrange(PlotOrdin(metadata_fungi, fungi_NMDS) + 
            labs(title = "NMDS") +
            annotate("text", -Inf, Inf, 
                     label = expression(paste("Site ", italic(R) ^ 2,"= 28.0%***"),
                                        parse=TRUE), size = 2.5, hjust = -0.1, vjust = 2.5),
          PlotBetadisper(permdisp_fungi_site, "Site", palette3, my_labels = "", metadata_fungi) + 
             labs(title = "Dispersion", y ="Distance to centroids", x= "Site"),
          PlotOrdin(metadata_bact, bact_NMDS) + 
            labs(title = "NMDS") +
            annotate("text", -Inf, Inf, 
                     label = expression(paste("Site ", italic(R) ^ 2,"= 23.6%***"),
                                        parse=TRUE), size = 2.5, hjust = -0.1, vjust = 2.5),
          PlotBetadisper(permdisp_bact_site, "Site", palette3, as.character(label_bact_site$Letters), metadata_bact) + 
            labs(title = "Dispersion", y ="Distance to centroids", x= "Site"),
          labels = c("A","B", "C","D"),
          widths = c(1, 0.5, 1, 0.5),
          align = "v" ,
          common.legend = TRUE,
          ncol = 4, 
          nrow = 1, 
          legend = c("bottom")) -> Fig_1_all

Fig_1_all

# PROCRUSTES ROTATION ANALYSIS ------------------------------------------------------------------------------------
require(ade4)
require(vegan)

pro_fungi_bact <- protest(fungi_NMDS,
                          bact_NMDS,
                          scores = "sites",
                          permutations = how(nperm = 9999))
pro_fungi_bact
pro_fungi_bact$ss

# Plotting --------------------------------------------------------------------------------------------------------
pro_plot <- data.frame(
  nmds1 = pro_fungi_bact$Yrot[, 1],
  nmds2 = pro_fungi_bact$Yrot[, 2],
  xnmds1 = pro_fungi_bact$X[, 1],
  xnmds2 = pro_fungi_bact$X[, 2]
)

pro_plot

identical(rownames(pro_plot), rownames(metadata_fungi))
pro_plot$Site <- metadata_fungi$Site
pro_plot$Treatment <- metadata_fungi$Treatment

# *** FIGURE 2 - Procustes Rotation Ordinations -------------------------------------------------------------------
ggplot(pro_plot) +
      theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
      theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
      theme(strip.text.x = element_text(size = 9, face = "bold")) +
      theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5)) + 
      theme(axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
      theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
      theme(axis.title = element_text(angle = 0, size = 9, face = "bold")) +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.05)) + # to adjust decimals
      scale_y_continuous(labels = scales::number_format(accuracy = 0.05)) +
      grids(linetype = "dashed") +
      theme(panel.background = element_blank(),
            legend.position = "bottom",
            legend.background = element_blank(),
            legend.key = element_blank(),
            axis.line = element_line()) +
      geom_point(aes(x=nmds1, y=nmds2, colour=Site, shape=Treatment), size=1.5, alpha=0.9) +
      geom_point(aes(x=xnmds1, y=xnmds2, colour=Site, shape=Treatment), size=1.5, alpha=0.9) +
      geom_segment(aes(x=xnmds1,y=xnmds2,xend=nmds1,yend=nmds2,colour=Site), size=0.5,
                   arrow=arrow(length=unit(0.3,"cm")), show.legend=F) +
      geom_hline(yintercept = 0,linetype="dashed") +
      geom_vline(xintercept = 0,linetype="dashed") +
      geom_abline(slope = -1/pro_fungi_bact$rotation[1,2]) +
      geom_abline(slope = pro_fungi_bact$rotation[1,2]) +
      scale_colour_manual(values=palette3, name = "Market") +
      scale_shape_manual(values = c(0,1,2,5)) +
  theme(legend.position="right") -> Fig_3_plot_proc

Fig_3_plot_proc

Fig_3_plot_proc +
annotate("text", -Inf, Inf, label = "italic(m) ^ 2 == 0.3812", parse = TRUE, size = 3, hjust = -0.2, vjust = 1.5) +
annotate("text", -Inf, Inf, label = "italic(r) ^ 2 == 0.7866", parse = TRUE, size = 3, hjust = -0.2, vjust = 3.0) +
annotate("text", -Inf, Inf, label = "italic(p) == 0.0001", parse = TRUE, size = 3, hjust = -0.2, vjust = 5) 

# *****************************************************************************************------------------------
# RANDOM FOREST MODELS --------------------------------------------------------------------------------------------
library(Boruta)
library(randomForest)

metadata_fungi$richness <-
  specnumber(t(df_fungi[, 1:27]), MARGIN = 1)
metadata_fungi$ecm_richness <-
  specnumber(t(df_ecm[, 1:27]), MARGIN = 1)
metadata_fungi$shannon <-
  diversity(t(df_fungi[, 1:27]), index = "shannon", MARGIN = 1)
metadata_fungi$ecm_shannon <-
  diversity(t(df_ecm[, 1:27]), index = "shannon", MARGIN = 1)
identical(rownames(metadata_fungi), rownames(myceliumDNA))
metadata_fungi$ngDNA <- myceliumDNA$ng
head(metadata_fungi)

metadata_bact$richness <- specnumber(t(df_bac[, 1:27]), MARGIN = 1)
metadata_bact$shannon <-
  diversity(t(df_bac[, 1:27]), index = "shannon", MARGIN = 1)
head(metadata_bact)

# >>> Model1: ECM Fungi -> Bacterial richness ---------------------------------------------------------------------
# transformation for microbiome sequencing data - used for compositional 
library(compositions)

# best tranfromation is always the z-scores 
df_fungi_tran <- scale(t(df_ecm[,1:27]), center = TRUE, scale = TRUE) # each value is converted into a Z-score
# df_fungi_tran <- scale(asinh(t(df_fungi[,1:27])), center=TRUE, scale=FALSE)  #inverse hyperbolic sine and then to mean center by sample
# df_fungi_tran <- clr(t(df_fungi[,1:27]))  # centred log-ratio (CLR)

df_fungi_tran <- as.data.frame(df_fungi_tran)
dim(df_fungi_tran)
range(df_fungi_tran)

# adding the classification variable and sample names -------------------------------------------------------------
identical(rownames(metadata_bact), rownames(df_fungi_tran))
df_fungi_tran$Richness <- metadata_bact$richness # Bacterial richness
head(df_fungi_tran)
dim(df_fungi_tran)

set.seed(1)
rfe_ecm_rich <- Boruta(
  Richness ~ .,
  df_fungi_tran,
  pValue = 0.05,
  mcAdj = TRUE,
  maxRuns = ncol(df_fungi_tran),
  doTrace = 3
)

rfe_ecm_rich

# Get significant variables including tentatives
rfe_ecm_rich <-
  getSelectedAttributes(rfe_ecm_rich, withTentative = TRUE)
rfe_ecm_rich

df_fungi_tran[, rfe_ecm_rich] -> df_fungi_tran_sel1
identical(rownames(df_fungi_tran_sel1), rownames(df_fungi_tran))
df_fungi_tran_sel1$Richness <- df_fungi_tran$Richness

# try tuning the model first
round(sqrt(ncol(df_fungi_tran_sel1[, 1:(ncol(df_fungi_tran_sel1) - 1)])))

set.seed(2)
bestmtry_ecm_rich <-
  tuneRF(
    x = df_fungi_tran_sel1[, 1:(ncol(df_fungi_tran_sel1) - 1)],
    y = df_fungi_tran_sel1$Richness,
    mtryStart = 2,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )

bestmtry_ecm_rich

set.seed(3)
RF_ecm_rich1 <-
  randomForest(
    x = df_fungi_tran_sel1[, 1:(ncol(df_fungi_tran_sel1) - 1)],
    y = df_fungi_tran_sel1$Richness,
    ntree = 1001,
    mtry = 1,
    importance = TRUE,
    proximity = TRUE
  )

RF_ecm_rich1

# >>> Model2: Bacteria -> ECM richness ----------------------------------------------------------------------------
df_bact_tran <- scale(t(df_bac[,1:27]), center = TRUE, scale = TRUE)
df_bact_tran <- as.data.frame(df_bact_tran)
dim(df_bact_tran)
range(df_bact_tran)

identical(rownames(metadata_fungi), rownames(df_bact_tran))
df_bact_tran$Richness <- metadata_fungi$ecm_richness # ECM Richness
head(df_bact_tran)
dim(df_bact_tran)

set.seed(5)
rfe_bact_ecm_rich <- Boruta(
  Richness ~ .,
  df_bact_tran,
  pValue = 0.05,
  mcAdj = TRUE,
  maxRuns = ncol(df_bact_tran),
  doTrace = 3
)

rfe_bact_ecm_rich

rfe_bact_ecm_rich <-
  getSelectedAttributes(rfe_bact_ecm_rich, withTentative = TRUE)
rfe_bact_ecm_rich

df_bact_tran[, rfe_bact_ecm_rich] -> df_bact_tran_sel2
identical(rownames(df_bact_tran_sel2), rownames(df_bact_tran))
df_bact_tran_sel2$Richness <- df_bact_tran$Richness

# try tuning the model first
round(sqrt(ncol(df_bact_tran_sel2[,1:(ncol(df_bact_tran_sel2)-1)])))

set.seed(6)
bestmtry_bact_ecm_rich <-
  tuneRF(
    x = df_bact_tran_sel2[, 1:(ncol(df_bact_tran_sel2) - 1)],
    y = df_bact_tran_sel2$Richness,
    mtryStart = 2,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )

bestmtry_bact_ecm_rich

set.seed(6)
RF_bact_ecm_rich2 <-
  randomForest(
    x = df_bact_tran_sel2[, 1:(ncol(df_bact_tran_sel2) - 1)],
    y = df_bact_tran_sel2$Richness,
    ntree = 1001,
    mtry = 1,
    importance = TRUE,
    proximity = TRUE
  )

RF_bact_ecm_rich2
plot(RF_bact_ecm_rich2)

# >>> Model3: ECM Fungi -> Bacterial Shannon ----------------------------------------------------------------------
df_fungi_tran <-
  scale(t(df_ecm[, 1:27]), center = TRUE, scale = TRUE)
df_fungi_tran <- as.data.frame(df_fungi_tran)
dim(df_fungi_tran)
range(df_fungi_tran)

# adding the classification variable and sample names -------------------------------------------------------------
identical(rownames(metadata_bact), rownames(df_fungi_tran))
df_fungi_tran$Shannon <- metadata_bact$shannon # Bacterial Shannon
head(df_fungi_tran)
dim(df_fungi_tran)


set.seed(10)
rfe_ecm_shan <-  Boruta(
  Shannon ~ .,
  df_fungi_tran,
  pValue = 0.05,
  mcAdj = TRUE,
  maxRuns = ncol(df_fungi_tran),
  doTrace = 3
)

rfe_ecm_shan

# Get significant variables including tentatives
rfe_ecm_shan <-
  getSelectedAttributes(rfe_ecm_shan, withTentative = TRUE)
rfe_ecm_shan

df_fungi_tran[, rfe_ecm_shan] -> df_fungi_tran_sel3
identical(rownames(df_fungi_tran_sel3), rownames(df_fungi_tran))
df_fungi_tran_sel3$Shannon <- df_fungi_tran$Shannon
dim(df_fungi_tran_sel3)

# try tuning the model first
round(sqrt(ncol(df_fungi_tran_sel3[, 1:(ncol(df_fungi_tran_sel3) - 1)])))

set.seed(11)
bestmtry_ecm_shan <-
  tuneRF(
    x = df_fungi_tran_sel3[, 1:(ncol(df_fungi_tran_sel3) - 1)],
    y = df_fungi_tran_sel3$Shannon,
    mtryStart = 3,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )

bestmtry_ecm_shan

set.seed(12)
RF_ecm_shan3 <-
  randomForest(
    x = df_fungi_tran_sel3[, 1:(ncol(df_fungi_tran_sel3) - 1)],
    y = df_fungi_tran_sel3$Shannon,
    ntree = 1001,
    mtry = 3,
    importance = TRUE,
    proximity = TRUE
  )

RF_ecm_shan3
plot(RF_ecm_shan3)

# >>> Model4: Bacteria -> ECM Shannon -----------------------------------------------------------------------------
df_bact_tran <- scale(t(df_bac[,1:27]), center = TRUE, scale = TRUE)
df_bact_tran <- as.data.frame(df_bact_tran)
dim(df_bact_tran)
range(df_bact_tran)

identical(rownames(metadata_fungi), rownames(df_bact_tran))
df_bact_tran$Shannon <- metadata_fungi$ecm_shannon # ECM Shannon
head(df_bact_tran)
dim(df_bact_tran)

set.seed(20)
rfe_bact_ecm_shan <- Boruta(
  Shannon ~ .,
  df_fungi_tran,
  pValue = 0.05, 
  mcAdj = TRUE,
  maxRuns = ncol(df_bact_tran),
  doTrace = 3)

rfe_bact_ecm_shan

# Get significant variables including tentatives
rfe_bact_ecm_shan <- getSelectedAttributes(rfe_bact_ecm_shan, withTentative = TRUE)
rfe_bact_ecm_shan

df_bact_tran[,rfe_bact_ecm_shan] -> df_bact_tran_sel4
identical(rownames(df_bact_tran_sel4), rownames(df_bact_tran))
df_bact_tran_sel4$Shannon <- df_bact_tran$Shannon

# try tuning the model first
round(sqrt(ncol(df_bact_tran_sel4[,1:(ncol(df_bact_tran_sel4)-1)])))

set.seed(21)
bestmtry_bact_ecm_shan <-
  tuneRF(
    x = df_bact_tran_sel4[, 1:(ncol(df_bact_tran_sel4) - 1)],
    y = df_bact_tran_sel4$Shannon,
    mtryStart = 2,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )

bestmtry_bact_ecm_shan

set.seed(22)
RF_bact_ecm_shan4 <-
  randomForest(
    x = df_bact_tran_sel4[, 1:(ncol(df_bact_tran_sel4) - 1)],
    y = df_bact_tran_sel4$Shannon,
    ntree = 1001,
    mtry = 1,
    importance = TRUE,
    proximity = TRUE
  )

RF_bact_ecm_shan4
plot(RF_bact_ecm_shan4)

# *********************************************************************************--------------------------------
# >>> Model5: All Fungi -> Bacterial richness ---------------------------------------------------------------------
df_fungi_tran <-
  scale(t(df_fungi[, 1:27]), center = TRUE, scale = TRUE)
df_fungi_tran <- as.data.frame(df_fungi_tran)
dim(df_fungi_tran)
range(df_fungi_tran)

# adding the classification variable and sample names -------------------------------------------------------------
identical(rownames(metadata_bact), rownames(df_fungi_tran))
df_fungi_tran$Richness <-
  metadata_bact$richness # Bacterial richness
head(df_fungi_tran)
dim(df_fungi_tran)

set.seed(30)
rfe_fungi_rich <- Boruta(
  Richness ~ .,
  df_fungi_tran,
  pValue = 0.05,
  mcAdj = TRUE,
  maxRuns = ncol(df_fungi_tran),
  doTrace = 3
)

rfe_fungi_rich

# Get significant variables including tentatives
rfe_fungi_rich <-
  getSelectedAttributes(rfe_fungi_rich, withTentative = TRUE)
rfe_fungi_rich

df_fungi_tran[,rfe_fungi_rich] -> df_fungi_tran_sel5
identical(rownames(df_fungi_tran_sel5), rownames(df_fungi_tran))
df_fungi_tran_sel5$Richness <- df_fungi_tran$Richness
dim(df_fungi_tran_sel5)

# try tuning the model first
round(sqrt(ncol(df_fungi_tran_sel5[, 1:(ncol(df_fungi_tran_sel5) - 1)])))

set.seed(31)
bestmtry_fungi_rich <-
  tuneRF(
    x = df_fungi_tran_sel5[, 1:(ncol(df_fungi_tran_sel5) - 1)],
    y = df_fungi_tran_sel5$Richness,
    mtryStart = 3,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )

bestmtry_fungi_rich

set.seed(32)
RF_fungi_rich5 <-
  randomForest(
    x = df_fungi_tran_sel5[, 1:(ncol(df_fungi_tran_sel5) - 1)],
    y = df_fungi_tran_sel5$Richness,
    ntree = 1001,
    mtry = 1,
    importance = TRUE,
    proximity = TRUE
  )

RF_fungi_rich5
plot(RF_fungi_rich5)

# >>> Model6: All Bacteria -> Fungal richness ---------------------------------------------------------------------
df_bact_tran <- scale(t(df_bac[, 1:27]), center = TRUE, scale = TRUE)
df_bact_tran <- as.data.frame(df_bact_tran)
dim(df_bact_tran)
range(df_bact_tran)

identical(rownames(metadata_fungi), rownames(df_bact_tran))
df_bact_tran$Richness <- metadata_fungi$richness # Fungal Richness
head(df_bact_tran)
dim(df_bact_tran)

set.seed(40)
rfe_bact_rich <- Boruta(
  Richness ~ .,
  df_fungi_tran,
  pValue = 0.05,
  mcAdj = TRUE,
  maxRuns = ncol(df_bact_tran),
  doTrace = 3
)


rfe_bact_rich

# Get significant variables including tentatives
rfe_bact_rich <-
  getSelectedAttributes(rfe_bact_rich, withTentative = TRUE)
rfe_bact_rich

df_bact_tran[, rfe_bact_rich] -> df_bact_tran_sel6
identical(rownames(df_bact_tran_sel6), rownames(df_bact_tran))
df_bact_tran_sel6$Richness <- df_bact_tran$Richness

# try tuning the model first
round(sqrt(ncol(df_bact_tran_sel6[, 1:(ncol(df_bact_tran_sel6) - 1)])))

set.seed(41)
bestmtry_bact_rich <-
  tuneRF(
    x = df_bact_tran_sel6[, 1:(ncol(df_bact_tran_sel6) - 1)],
    y = df_bact_tran_sel6$Richness,
    mtryStart = 3,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )

bestmtry_bact_rich

set.seed(42)
RF_bact_rich6 <-
  randomForest(
    x = df_bact_tran_sel6[, 1:(ncol(df_bact_tran_sel6) - 1)],
    y = df_bact_tran_sel6$Richness,
    ntree = 1001,
    mtry = 1,
    importance = TRUE,
    proximity = TRUE
  )

RF_bact_rich6
plot(RF_bact_rich6)

# >>> Model7: All Fungi -> Bacterial Shannon ----------------------------------------------------------------------
df_fungi_tran <-
  scale(t(df_fungi[, 1:27]), center = TRUE, scale = TRUE)
df_fungi_tran <- as.data.frame(df_fungi_tran)
dim(df_fungi_tran)
range(df_fungi_tran)

# adding the classification variable and sample names -------------------------------------------------------------
identical(rownames(metadata_bact), rownames(df_fungi_tran))
df_fungi_tran$Shannon <- metadata_bact$shannon # Bacterial shannon
head(df_fungi_tran)
dim(df_fungi_tran)

set.seed(50)
rfe_fungi_shan <- Boruta(
  Shannon ~ .,
  df_fungi_tran,
  pValue = 0.05, 
  mcAdj = TRUE,
  maxRuns = ncol(df_fungi_tran),
  doTrace = 3)

rfe_fungi_shan

# Get significant variables including tentatives
rfe_fungi_shan <-
  getSelectedAttributes(rfe_fungi_shan, withTentative = TRUE)
rfe_fungi_shan

df_fungi_tran[, rfe_fungi_shan] -> df_fungi_tran_sel7
identical(rownames(df_fungi_tran_sel7), rownames(df_fungi_tran))
df_fungi_tran_sel7$Shannon <- df_fungi_tran$Shannon
dim(df_fungi_tran_sel7)

# try tuning the model first
round(sqrt(ncol(df_fungi_tran_sel7[, 1:(ncol(df_fungi_tran_sel7) - 1)])))

set.seed(51)
bestmtry_fungi_shan <-
  tuneRF(
    x = df_fungi_tran_sel7[, 1:(ncol(df_fungi_tran_sel7) - 1)],
    y = df_fungi_tran_sel7$Shannon,
    mtryStart = 4,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )

bestmtry_fungi_shan

set.seed(52)
RF_fungi_shan7 <-
  randomForest(
    x = df_fungi_tran_sel7[, 1:(ncol(df_fungi_tran_sel7) - 1)],
    y = df_fungi_tran_sel7$Shannon,
    ntree = 1001,
    mtry = 1,
    importance = TRUE,
    proximity = TRUE
  )

RF_fungi_shan7
plot(RF_fungi_shan7)

# >>> Model8: Bacteria -> Fungal shannon --------------------------------------------------------------------------
df_bact_tran <- scale(t(df_bac[, 1:27]), center = TRUE, scale = TRUE)
df_bact_tran <- as.data.frame(df_bact_tran)
dim(df_bact_tran)
range(df_bact_tran)

identical(rownames(metadata_fungi), rownames(df_bact_tran))
df_bact_tran$Shannon <- metadata_fungi$shannon # Fungal Shannon
head(df_bact_tran)
dim(df_bact_tran)

set.seed(60)
rfe_bact_shan <- Boruta(
  Shannon ~ .,
  df_bact_tran,
  pValue = 0.05,
  mcAdj = TRUE,
  maxRuns = ncol(df_bact_tran),
  doTrace = 3
)

rfe_bact_shan

# Get significant variables including tentatives
rfe_bact_shan <-
  getSelectedAttributes(rfe_bact_shan, withTentative = TRUE)
rfe_bact_shan

df_bact_tran[, rfe_bact_rich] -> df_bact_tran_sel8
identical(rownames(df_bact_tran_sel8), rownames(df_bact_tran))
df_bact_tran_sel8$Shannon <- df_bact_tran$Shannon

# try tuning the model first
round(sqrt(ncol(df_bact_tran_sel8[, 1:(ncol(df_bact_tran_sel8) - 1)])))

set.seed(61)
bestmtry_bact_shan <-
  tuneRF(
    x = df_bact_tran_sel8[, 1:(ncol(df_bact_tran_sel8) - 1)],
    y = df_bact_tran_sel8$Shannon,
    mtryStart = 2,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )

bestmtry_bact_shan

set.seed(62)
RF_bact_shan8 <-
  randomForest(
    x = df_bact_tran_sel8[, 1:(ncol(df_bact_tran_sel8) - 1)],
    y = df_bact_tran_sel8$Shannon,
    ntree = 1001,
    mtry = 1,
    importance = TRUE,
    proximity = TRUE
  )

RF_bact_shan8
plot(RF_bact_shan8)

# *********************************************************************************--------------------------------
# Modelying mycelium DNA

# >>> Model9: Bacteria -> Boletus DNA -----------------------------------------------------------------------------
df_bact_tran <- scale(t(df_bac[, 1:27]), center = TRUE, scale = TRUE)
df_bact_tran <- as.data.frame(df_bact_tran)
dim(df_bact_tran)
range(df_bact_tran)

identical(rownames(metadata_fungi), rownames(df_bact_tran))
df_bact_tran$DNA <- metadata_fungi$ngDNA # ECM Richness
head(df_bact_tran)
dim(df_bact_tran)

set.seed(70)
rfe_bact_DNA <- Boruta(
  DNA ~ .,
  df_bact_tran,
  pValue = 0.05,
  mcAdj = TRUE,
  maxRuns = ncol(df_bact_tran),
  doTrace = 3
)

rfe_bact_DNA

# Get significant variables including tentatives
rfe_bact_DNA <-
  getSelectedAttributes(rfe_bact_DNA, withTentative = TRUE)
rfe_bact_DNA

df_bact_tran[,rfe_bact_DNA] -> df_bact_tran_sel9
identical(rownames(df_bact_tran_sel9), rownames(df_bact_tran))
df_bact_tran_sel9$DNA <- df_bact_tran$DNA

# try tuning the model first
round(sqrt(ncol(df_bact_tran_sel9[, 1:(ncol(df_bact_tran_sel9) - 1)])))

set.seed(71)
bestmtry_bact_DNA <-
  tuneRF(
    x = df_bact_tran_sel9[, 1:(ncol(df_bact_tran_sel9) - 1)],
    y = df_bact_tran_sel9$DNA,
    mtryStart = 2,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )

bestmtry_bact_DNA

set.seed(72)
RF_bact_DNA9 <-
  randomForest(
    x = df_bact_tran_sel9[, 1:(ncol(df_bact_tran_sel9) - 1)],
    y = df_bact_tran_sel9$DNA,
    ntree = 1001,
    mtry = 1,
    importance = TRUE,
    proximity = TRUE
  )

RF_bact_DNA9
plot(RF_bact_DNA9)

# >>> Model10: Fungi -> Boletus DNA -------------------------------------------------------------------------------
identical(rownames(df_fungi), rownames(tax_fungi))
tax_fungi[tax_fungi$Genus == "Boletus",]
df_fungi[!(rownames(df_fungi) %in% "FOTU_390"), ] -> df_fungi_filt
dim(df_fungi_filt)

df_fungi_filt_tran <-
  scale(t(df_fungi_filt[, 1:27]), center = TRUE, scale = TRUE)
df_fungi_filt_tran <- as.data.frame(df_fungi_filt_tran)
dim(df_fungi_filt_tran)
range(df_fungi_filt_tran)

identical(rownames(metadata_fungi), rownames(df_fungi_filt_tran))
df_fungi_filt_tran$DNA <- metadata_fungi$ngDNA # ECM Richness
head(df_fungi_filt_tran)
dim(df_fungi_filt_tran)

set.seed(80)
rfe_fungi_filt_DNA <- Boruta(
  DNA ~ .,
  df_fungi_filt_tran,
  pValue = 0.05,
  mcAdj = TRUE,
  maxRuns = ncol(df_fungi_filt_tran),
  doTrace = 3
)

rfe_fungi_filt_DNA

# Get significant variables including tentatives
rfe_fungi_filt_DNA <-
  getSelectedAttributes(rfe_fungi_filt_DNA, withTentative = TRUE)
rfe_fungi_filt_DNA

df_fungi_filt_tran[, rfe_fungi_filt_DNA] -> df_fungi_filt_tran_sel10
identical(rownames(df_fungi_filt_tran_sel10),
          rownames(df_fungi_filt_tran))
df_fungi_filt_tran_sel10$DNA <- df_fungi_filt_tran$DNA
dim(df_fungi_filt_tran)

# try tuning the model first
round(sqrt(ncol(df_fungi_filt_tran_sel10[, 1:(ncol(df_fungi_filt_tran_sel10) -
                                                1)])))

set.seed(81)
bestmtry_fungi_filt_DNA <-
  tuneRF(
    x = df_fungi_filt_tran_sel10[, 1:(ncol(df_fungi_filt_tran_sel10) - 1)],
    y = df_fungi_filt_tran_sel10$DNA,
    mtryStart = 3,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )

bestmtry_fungi_filt_DNA

set.seed(82)
RF_fungi_filt_DNA10 <-
  randomForest(
    x = df_fungi_filt_tran_sel10[, 1:(ncol(df_fungi_filt_tran_sel10) - 1)],
    y = df_fungi_filt_tran_sel10$DNA,
    ntree = 1001,
    mtry = 3,
    importance = TRUE,
    proximity = TRUE
  )

RF_fungi_filt_DNA10
plot(RF_fungi_filt_DNA10)

# >>> Model11: ECM Fungi -> Boletus DNA ---------------------------------------------------------------------------
identical(rownames(df_ecm), rownames(tax_ecm))
tax_ecm[tax_ecm$Genus == "Boletus", ]
df_ecm[!(rownames(df_ecm) %in% "FOTU_390"),] -> df_ecm_filt
dim(df_ecm_filt)

df_ecm_filt_tran <-
  scale(t(df_ecm_filt[, 1:27]), center = TRUE, scale = TRUE)
df_ecm_filt_tran <- as.data.frame(df_ecm_filt_tran)
dim(df_ecm_filt_tran)
range(df_ecm_filt_tran)

identical(rownames(metadata_fungi), rownames(df_ecm_filt_tran))
df_ecm_filt_tran$DNA <- metadata_fungi$ngDNA # ECM Richness
head(df_ecm_filt_tran)
dim(df_ecm_filt_tran)

set.seed(90)
rfe_ecm_filt_DNA <- Boruta(
  DNA ~ .,
  df_ecm_filt_tran,
  pValue = 0.05,
  mcAdj = TRUE,
  maxRuns = ncol(df_ecm_filt_tran),
  doTrace = 3
)

rfe_ecm_filt_DNA

# Get significant variables including tentatives
rfe_ecm_filt_DNA <-
  getSelectedAttributes(rfe_ecm_filt_DNA, withTentative = TRUE)
rfe_ecm_filt_DNA

df_ecm_filt_tran[, rfe_ecm_filt_DNA[2:length(rfe_ecm_filt_DNA)]] -> df_ecm_filt_tran_sel11
identical(rownames(df_ecm_filt_tran_sel11),
          rownames(df_ecm_filt_tran))
df_ecm_filt_tran_sel11$DNA <- df_ecm_filt_tran$DNA
dim(df_ecm_filt_tran)

# try tuning the model first
round(sqrt(ncol(df_ecm_filt_tran_sel11[, 1:(ncol(df_ecm_filt_tran_sel11) -
                                              1)])))

set.seed(91)
bestmtry_ecm_filt_DNA <-
  tuneRF(
    x = df_ecm_filt_tran_sel11[, 1:(ncol(df_ecm_filt_tran_sel11) - 1)],
    y = df_ecm_filt_tran_sel11$DNA,
    mtryStart = 3,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )

bestmtry_ecm_filt_DNA


set.seed(92)
RF_ecm_filt_DNA11 <-
  randomForest(
    x = df_ecm_filt_tran_sel11[, 1:(ncol(df_ecm_filt_tran_sel11) - 1)],
    y = df_ecm_filt_tran_sel11$DNA,
    ntree = 1001,
    mtry = 3,
    importance = TRUE,
    proximity = TRUE
  )

RF_ecm_filt_DNA11
plot(RF_ecm_filt_DNA11)

# MODELS SIGNIFICANCE using permitations --------------------------------------------------------------------------
PermRFTest <- function(RF, dataframe, tryRF) {
  require(rfUtilities)
  set.seed(100)
  perm_RF <-
    rf.significance(
      x = RF,
      xdata = dataframe[, 1:(ncol(dataframe) - 1)],
      nperm = 999,
      nmtry = tryRF,
      ntree = 1001
    )
  return(perm_RF)
}

PermRFTest(RF_ecm_rich1, df_fungi_tran_sel1, 1) #p-value:  0.001
PermRFTest(RF_bact_ecm_rich2, df_bact_tran_sel2, 1) #p-value:  0.001
PermRFTest(RF_ecm_shan3, df_fungi_tran_sel3, 3) #p-value:  0.001
PermRFTest(RF_bact_ecm_shan4, df_bact_tran_sel4, 1) #p-value:  0.003

PermRFTest(RF_fungi_rich5, df_fungi_tran_sel5, 1) #p-value:  0.001
PermRFTest(RF_bact_rich6, df_bact_tran_sel6, 1) #p-value:  0.001
PermRFTest(RF_fungi_shan7, df_fungi_tran_sel7, 1) #p-value:  0.001
PermRFTest(RF_bact_shan8, df_bact_tran_sel8, 1) #p-value:  0.001

PermRFTest(RF_bact_DNA9, df_bact_tran_sel9, 1) #p-value:  0.001
PermRFTest(RF_fungi_filt_DNA10, df_fungi_filt_tran_sel10, 3) #p-value:  0.001
PermRFTest(RF_ecm_filt_DNA11, df_ecm_filt_tran_sel11, 1) #p-value:  0.001


# PLOTTING BEST MODELS ---------------------------------------------------------------------------------------------

# 1 - PLOTTING LINE ------------------------------------------------------------------------------------------------
PlotLine <- function(rf_model, metadata) {
  require(ggpmisc)
  require(ggpubr)
  df_model <- data.frame(actual = rf_model$y, pred = rf_model$predicted)
  if (identical(rownames(metadata), rownames(df_model)) == TRUE) {
    df_model$Site <- metadata$Site
    ggplot(data = df_model, aes(x = actual, y = pred, color = Site)) +
      geom_point() +
      geom_smooth(
        method = "lm",
        formula = "y ~ x",
        se = FALSE,
        color = "black"
      ) +
      scale_colour_manual("Site", values = palette3) +
      theme_classic() +
      theme(legend.key.height = unit(0.1, "cm"), legend.key.width = unit(0.1, "cm")) +
      theme(legend.title = element_blank(), legend.text = element_text(size = 7)) +
      theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
      theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
      theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
      theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) -> line_plot
    return(line_plot)
  }else{
    stop("Error: dataframe and metadata rownames are not matching!")
  }
}

PlotLine(RF_ecm_rich1, metadata_fungi) + theme(
  legend.justification = c(1, 0),
  legend.position = c(1, 0),
  legend.background = element_blank()
) +
  labs(title = "Observed vs. Predicted", y = "Predicted richness", x = "Observed richness") 


# 2 - PLOTTING ERROR ----------------------------------------------------------------------------------------------
PlotError <- function(rf_model) {
  model_df <- (data.frame(Trees = 1:1001, Error = rf_model$mse))
  ggplot(data = model_df, aes(x = Trees, y = Error)) +
    labs(title = "Model Errors", y = "Error", x = "Tree") +
    theme_classic() +
    geom_line(color = "red", size = 1.2) + 
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size =8)) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) -> error_plot
  return(error_plot)
}


PlotError(RF_ecm_rich1) +
  annotate("text", x=Inf, y = Inf,
           label=paste("Mean squared error:", round(last(RF_ecm_rich1$mse), 2)), size=2.5, vjust=1, hjust=1) +
  annotate("text", x=Inf, y = Inf, 
           label= paste("% Var explained:", round(last(RF_ecm_rich1$rsq*100),2)), size=2.5, vjust=3, hjust=1) +
  annotate("text", x=Inf, y = Inf, label = "italic(p) == 0.001", parse = TRUE, size=2.5, vjust=4, hjust=1)

# R squared calculation
1 - (sum((RF_ecm_rich1$y-RF_ecm_rich1$predicted)^2)/sum((RF_ecm_rich1$y-mean(RF_ecm_rich1$y))^2))

# %IncMSE indicates the increase of the Mean Squared Error when 
# given variable is randomly permuted. %IncMSE - It is computed from permuting 
# test data: For each tree, the prediction error on test is recorded 
#(Mean Squared Error - MSE ). Then the same is done after permuting each
#predictor variable. The difference between the two are then averaged over all
#trees, and normalized by the standard deviation of the differences. If the
#standard deviation of the differences is equal to 0 for a variable, the
#division is not done (but the average is almost always equal to 0 in that case).
#Higher the difference is, more important the variable. MSE = mean((actual_y - predicted_y)^2)

# 3- PLOTTING MOST IMPORTANT FEATURES ------------------------------------------------------------------------------
PlotFeature <- function(rf_model, taxa){
  imp_RF <- as.data.frame(rf_model$importance)
  imp_RF$features <- rownames(imp_RF)
  imp_RF <- arrange(imp_RF, desc(`%IncMSE`))
  rownames(imp_RF) <- imp_RF$features
  # adding marker taxon info
  taxa[rownames(taxa) %in% rownames(imp_RF),] -> taxa_RF
  identical(rownames(taxa_RF), rownames(imp_RF))
  order_taxa <- match(rownames(taxa_RF), rownames(imp_RF))
  imp_RF <- imp_RF[order_taxa, ]
  imp_RF$Taxonomy <- taxa_RF$Taxon
  imp_RF <- imp_RF[order(imp_RF$`%IncMSE`, decreasing = TRUE), ] 
  ggplot(data=imp_RF) + 
    geom_bar(aes(x= reorder(Taxonomy, -`%IncMSE`),
                 y= `%IncMSE`), color="grey80", fill="grey80",stat="identity") +
    labs(x= "OTU", y= "% Increase in Mean Squared Error") +
    coord_flip() +
    #ylim(0, 0.03) +
    theme_classic() +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0)) + 
    theme(axis.text.x = element_text(angle = 0, size = 7 ,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 1, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) -> plot_importance
  return(plot_importance)
}

PlotFeature(RF_ecm_rich1, tax_ecm_filt) + labs(title = "Selected OTUs") 



# FINDING THE BEST MODELS -----------------------------------------------------------------------------------------
library(grid)
library(gridExtra)

# ECM community
title1=text_grob("ECM Fungi -> Bacterial Richness",size=12,face=2)
plots1 <- list(PlotError(RF_ecm_rich1) +
                 annotate("text", x=Inf, y = Inf,
                          label=paste("Mean squared error:", round(last(RF_ecm_rich1$mse), 2)), size=2.5, vjust=1, hjust=1) +
                 annotate("text", x=Inf, y = Inf, 
                          label= paste("% Var explained:", round(last(RF_ecm_rich1$rsq*100),2)), size=2.5, vjust=3, hjust=1) +
                 annotate("text", x=Inf, y = Inf, label = "italic(p) == 0.001", parse = TRUE, size=2.5, vjust=4, hjust=1),
               PlotLine(RF_ecm_rich1, metadata_fungi) +  
                 theme(legend.justification = c(1, 0),
                       legend.position = c(1, 0),
                       legend.background = element_blank()) +
                 labs(title = "Observed vs. Predicted", y="Predicted richness", x="Observed richness"),
               PlotFeature(RF_ecm_rich1, tax_ecm_filt) + labs(title = "Selected OTUs"))

grid.arrange(grobs = plots1, top = title1, ncol=3,  widths=c(2,2,3))

title2=text_grob("Bacteria -> ECM richness",size=12,face=2)
plots2 <- list(PlotError(RF_bact_ecm_rich2) +
                 annotate("text", x=Inf, y = Inf,
                          label=paste("Mean squared error:", round(last(RF_bact_ecm_rich2$mse), 2)), size=2.5, vjust=1, hjust=1) +
                 annotate("text", x=Inf, y = Inf, 
                          label= paste("% Var explained:", round(last(RF_bact_ecm_rich2$rsq*100),2)), size=2.5, vjust=3, hjust=1) +
                 annotate("text", x=Inf, y = Inf, label = "italic(p) == 0.001", parse = TRUE, size=2.5, vjust=4, hjust=1),
               PlotLine(RF_bact_ecm_rich2, metadata_bact) +  
                 theme(legend.justification = c(1, 0),
                       legend.position = c(1, 0),
                       legend.background = element_blank()) +
                 labs(title = "Observed vs. Predicted", y="Predicted richness", x="Observed richness"),
               PlotFeature(RF_bact_ecm_rich2, df_bac) + labs(title = "Selected OTUs"))

grid.arrange(grobs = plots2, top = title2, ncol=3,  widths=c(2,2,3))

title3=text_grob("ECM Fungi -> Bacterial Shannon",size=12,face=2)
plots3 <- list(PlotError(RF_ecm_shan3) +
                 annotate("text", x=Inf, y = Inf,
                          label=paste("Mean squared error:", round(last(RF_ecm_shan3$mse), 2)), size=2.5, vjust=1, hjust=1) +
                 annotate("text", x=Inf, y = Inf, 
                          label= paste("% Var explained:", round(last(RF_ecm_shan3$rsq*100),2)), size=2.5, vjust=3, hjust=1) +
                 annotate("text", x=Inf, y = Inf, label = "italic(p) == 0.001", parse = TRUE, size=2.5, vjust=4, hjust=1),
               PlotLine(RF_ecm_shan3, metadata_fungi) +  
                 theme(legend.justification = c(1, 0),
                       legend.position = c(1, 0),
                       legend.background = element_blank()) +
                 labs(title = "Observed vs. Predicted", y="Predicted Shannon", x="Observed Shannon"),
               PlotFeature(RF_ecm_shan3, tax_ecm_filt) + labs(title = "Selected OTUs"))

grid.arrange(grobs = plots3, top = title3, ncol=3,  widths=c(2,2,3))


title4=text_grob("Bacteria -> ECM Shannon",size=12,face=2)
plots4 <- list(PlotError(RF_bact_ecm_shan4) +
                 annotate("text", x=Inf, y = Inf,
                          label=paste("Mean squared error:", round(last(RF_bact_ecm_shan4$mse), 2)), size=2.5, vjust=1, hjust=1) +
                 annotate("text", x=Inf, y = Inf, 
                          label= paste("% Var explained:", round(last(RF_bact_ecm_shan4$rsq*100),2)), size=2.5, vjust=3, hjust=1) +
                 annotate("text", x=Inf, y = Inf, label = "italic(p) == 0.003", parse = TRUE, size=2.5, vjust=4, hjust=1),
               PlotLine(RF_bact_ecm_shan4, metadata_bact) +  
                 theme(legend.justification = c(1, 0),
                       legend.position = c(1, 0),
                       legend.background = element_blank()) +
                 labs(title = "Observed vs. Predicted", y="Predicted Shannon", x="Observed Shannon"),
               PlotFeature(RF_bact_ecm_shan4, df_bac) + labs(title = "Selected OTUs"))

grid.arrange(grobs = plots4, top = title4, ncol=3,  widths=c(2,2,3))

# WHOLE COMMUNITY -------------------------------------------------------------------------------------------------
title5=text_grob("All Fungi -> Bacterial richness",size=12,face=2)
plots5 <- list(PlotError(RF_fungi_rich5) +
                 annotate("text", x=Inf, y = Inf,
                          label=paste("Mean squared error:", round(last(RF_fungi_rich5$mse), 2)), size=2.5, vjust=1, hjust=1) +
                 annotate("text", x=Inf, y = Inf, 
                          label= paste("% Var explained:", round(last(RF_fungi_rich5$rsq*100),2)), size=2.5, vjust=3, hjust=1) +
                 annotate("text", x=Inf, y = Inf, label = "italic(p) == 0.003", parse = TRUE, size=2.5, vjust=4, hjust=1),
               PlotLine(RF_fungi_rich5, metadata_fungi) +  
                 theme(legend.justification = c(1, 0),
                       legend.position = c(1, 0),
                       legend.background = element_blank()) +
                 labs(title = "Observed vs. Predicted", y="Predicted richness", x="Observed richness"),
               PlotFeature(RF_fungi_rich5, tax_fungi) + labs(title = "Selected OTUs"))

grid.arrange(grobs = plots5, top = title5, ncol=3,  widths=c(2,2,3))


title6=text_grob("All Bacteria -> Fungal richness",size=12,face=2)
plots6 <- list(PlotError(RF_bact_rich6) +
                 annotate("text", x=Inf, y = Inf,
                          label=paste("Mean squared error:", round(last(RF_bact_rich6$mse), 2)), size=2.5, vjust=1, hjust=1) +
                 annotate("text", x=Inf, y = Inf, 
                          label= paste("% Var explained:", round(last(RF_bact_rich6$rsq*100),2)), size=2.5, vjust=3, hjust=1) +
                 annotate("text", x=Inf, y = Inf, label = "italic(p) == 0.001", parse = TRUE, size=2.5, vjust=4, hjust=1),
               PlotLine(RF_bact_rich6, metadata_bact) +  
                 theme(legend.justification = c(1, 0),
                       legend.position = c(1, 0),
                       legend.background = element_blank()) +
                 labs(title = "Observed vs. Predicted", y="Predicted richness", x="Observed richness"),
               PlotFeature(RF_bact_rich6, df_bac) + labs(title = "Selected OTUs"))

grid.arrange(grobs = plots6, top = title6, ncol=3,  widths=c(2,2,3))


title7=text_grob("All Fungi -> Bacterial Shannon",size=12,face=2)
plots7 <- list(PlotError(RF_fungi_shan7) +
                 annotate("text", x=Inf, y = Inf,
                          label=paste("Mean squared error:", round(last(RF_fungi_shan7$mse), 2)), size=2.5, vjust=1, hjust=1) +
                 annotate("text", x=Inf, y = Inf, 
                          label= paste("% Var explained:", round(last(RF_fungi_shan7$rsq*100),2)), size=2.5, vjust=3, hjust=1) +
                 annotate("text", x=Inf, y = Inf, label = "italic(p) == 0.001", parse = TRUE, size=2.5, vjust=4, hjust=1),
               PlotLine(RF_fungi_shan7, metadata_fungi) +  
                 theme(legend.justification = c(1, 0),
                       legend.position = c(1, 0),
                       legend.background = element_blank()) +
                 labs(title = "Observed vs. Predicted", y="Predicted Shannon", x="Observed Shannon"),
               PlotFeature(RF_fungi_shan7, tax_fungi) + labs(title = "Selected OTUs"))

grid.arrange(grobs = plots7, top = title7, ncol=3,  widths=c(2,2,3))


title8=text_grob("All Bacteria -> Fungal Shannon",size=12,face=2)
plots8 <- list(PlotError(RF_bact_shan8) +
                 annotate("text", x=Inf, y = Inf,
                          label=paste("Mean squared error:", round(last(RF_bact_shan8$mse), 2)), size=2.5, vjust=1, hjust=1) +
                 annotate("text", x=Inf, y = Inf, 
                          label= paste("% Var explained:", round(last(RF_bact_shan8$rsq*100),2)), size=2.5, vjust=3, hjust=1) +
                 annotate("text", x=Inf, y = Inf, label = "italic(p) == 0.001", parse = TRUE, size=2.5, vjust=4, hjust=1),
               PlotLine(RF_bact_shan8, metadata_bact) +  
                 theme(legend.justification = c(1, 0),
                       legend.position = c(1, 0),
                       legend.background = element_blank()) +
                 labs(title = "Observed vs. Predicted", y="Predicted Shannon", x="Observed Shannon"),
               PlotFeature(RF_bact_shan8, df_bac) + labs(title = "Selected OTUs"))

grid.arrange(grobs = plots8, top = title8, ncol=3,  widths=c(2,2,3))

# BOLETUS DNA ------------------------------------------------------------------------------------------------------

title9=text_grob("Bacteria -> Boletus DNA",size=12,face=2)
plots9 <- list(PlotError(RF_bact_DNA9) +
                 annotate("text", x=Inf, y = Inf,
                          label=paste("Mean squared error:", round(last(RF_bact_DNA9$mse), 2)), size=2.5, vjust=1, hjust=1) +
                 annotate("text", x=Inf, y = Inf, 
                          label= paste("% Var explained:", round(last(RF_bact_DNA9$rsq*100),2)), size=2.5, vjust=3, hjust=1) +
                 annotate("text", x=Inf, y = Inf, label = "italic(p) == 0.001", parse = TRUE, size=2.5, vjust=4, hjust=1),
               PlotLine(RF_bact_DNA9, metadata_bact) +  
                 theme(legend.justification = c(1, 0),
                       legend.position = c(1, 0),
                       legend.background = element_blank()) +
                 labs(title = "Observed vs. Predicted", y="Predicted DNA (ng)", x="Observed DNA (ng)"),
               PlotFeature(RF_bact_DNA9, df_bac) + labs(title = "Selected OTUs"))

grid.arrange(grobs = plots9, top = title9, ncol=3,  widths=c(2,2,3))

title10=text_grob("All Fungi -> Boletus DNA (no Boletus)",size=12,face=2)
plots10 <- list(PlotError(RF_fungi_filt_DNA10) +
                 annotate("text", x=Inf, y = Inf,
                          label=paste("Mean squared error:", round(last(RF_fungi_filt_DNA10$mse), 2)), size=2.5, vjust=1, hjust=1) +
                 annotate("text", x=Inf, y = Inf, 
                          label= paste("% Var explained:", round(last(RF_fungi_filt_DNA10$rsq*100),2)), size=2.5, vjust=3, hjust=1) +
                 annotate("text", x=Inf, y = Inf, label = "italic(p) == 0.001", parse = TRUE, size=2.5, vjust=4, hjust=1),
               PlotLine(RF_fungi_filt_DNA10, metadata_fungi) +  
                 theme(legend.justification = c(1, 0),
                       legend.position = c(1, 0),
                       legend.background = element_blank()) +
                 labs(title = "Observed vs. Predicted", y="Predicted DNA (ng)", x="Observed DNA (ng)"),
               PlotFeature(RF_fungi_filt_DNA10, tax_fungi) + labs(title = "Selected OTUs"))

grid.arrange(grobs = plots10, top = title10, ncol=3,  widths=c(2,2,3))


title11=text_grob("ECM Fungi -> Boletus DNA (no Boletus)",size=12,face=2)
plots11 <- list(PlotError(RF_ecm_filt_DNA11) +
                  annotate("text", x=Inf, y = Inf,
                           label=paste("Mean squared error:", round(last(RF_ecm_filt_DNA11$mse), 2)), size=2.5, vjust=1, hjust=1) +
                  annotate("text", x=Inf, y = Inf, 
                           label= paste("% Var explained:", round(last(RF_ecm_filt_DNA11$rsq*100),2)), size=2.5, vjust=3, hjust=1) +
                  annotate("text", x=Inf, y = Inf, label = "italic(p) == 0.001", parse = TRUE, size=2.5, vjust=4, hjust=1),
                PlotLine(RF_ecm_filt_DNA11, metadata_fungi) +  
                  theme(legend.justification = c(1, 0),
                        legend.position = c(1, 0),
                        legend.background = element_blank()) +
                  labs(title = "Observed vs. Predicted", y="Predicted DNA (ng)", x="Observed DNA (ng)"),
                PlotFeature(RF_ecm_filt_DNA11, tax_fungi) + labs(title = "Selected OTUs"))

grid.arrange(grobs = plots11, top = title11, ncol=3,  widths=c(2,2,3))

# selected plots 
grid.arrange(grobs = plots3, top = title3, ncol=3,  widths=c(2,2,3))
grid.arrange(grobs = plots4, top = title4, ncol=3,  widths=c(2,2,3))
grid.arrange(grobs = plots7, top = title7, ncol=3,  widths=c(2,2,3))
grid.arrange(grobs = plots8, top = title8, ncol=3,  widths=c(2,2,3))

grid.arrange(grobs = plots9, top = title9, ncol=3,  widths=c(2,2,3))
grid.arrange(grobs = plots10, top = title10, ncol=3,  widths=c(2,2,3)) # no Boletus otu
grid.arrange(grobs = plots11, top = title11, ncol=3,  widths=c(2,2,3))  # no Boletus otu

#**** FIGURE RF models ---------------------------------------------------------------------------------------------
ggarrange(grid.arrange(grobs = plots3, top = title3, ncol=3,  widths=c(2,2,3)),
          grid.arrange(grobs = plots4, top = title4, ncol=3,  widths=c(2,2,3)),
          grid.arrange(grobs = plots7, top = title7, ncol=3,  widths=c(2,2,3)),
          grid.arrange(grobs = plots8, top = title8, ncol=3,  widths=c(2,2,3)),
          labels = c("A","B","C","D"),
          align = "hv" ,
          ncol = 1, 
          nrow =4) -> Fig_RF_div

Fig_RF_div

ggarrange(grid.arrange(grobs = plots1, top = title1, ncol=3,  widths=c(2,2,3)),
          grid.arrange(grobs = plots2, top = title2, ncol=3,  widths=c(2,2,3)),
          grid.arrange(grobs = plots5, top = title5, ncol=3,  widths=c(2,2,3)),
          grid.arrange(grobs = plots6, top = title6, ncol=3,  widths=c(2,2,3)),
          labels = c("A","B","C","D"),
          align = "hv" ,
          ncol = 1, 
          nrow =4) -> Fig_RF_rich

Fig_RF_rich

ggarrange(grid.arrange(grobs = plots9, top = title9, ncol=3,  widths=c(2,2,3)),
          grid.arrange(grobs = plots10, top = title10, ncol=3,  widths=c(2,2,3)),
          grid.arrange(grobs = plots11, top = title11, ncol=3,  widths=c(2,2,3)),
          labels = c("A","B","C"),
          align = "hv" ,
          ncol = 1, 
          nrow = 3) -> Fig_RF_DNA

Fig_RF_DNA



