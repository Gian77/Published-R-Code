# ************ DATA ANALYSIS ***************************************** -------------------------------------
# Project name: Switchgrass phyllosphere microbiome 
# Manuscript:   Host genetic control of succession in the switchgrass leaf fungal microbiome
# Authors:      VanWallendael A, Benucci GMN, ..., Lowry DB.
# Affiliation:  Michigan State University
# Journal:      Nature Comm.
# Date:         March 12, 2021
#
# R Code:       Benucci GMN
# Please contact <gian.benucci@gmail.com> if you intend using this code. 
# ******************************************************************** --------------------------------------

# loading required packages ---------------------------------------------------------------------------------
library(phyloseq); packageVersion("phyloseq")
library(Biostrings)
library(ggplot2)
library(ape)
library(dplyr)
library(magrittr)
library(styler) 
library(ggpubr)
library(vegan)
library(SpiecEasi)
library(tidygraph)
library(igraph)
library(ggraph)
library(gridExtra)


# Importing Acer data ---------------------------------------------------------------------------------------
load("decont_all_lib_jan31_last.rda")
otu_table(decont_all_lib) <-
  otu_table(decont_all_lib)[which(rowSums(otu_table(decont_all_lib)) > 0), ]

load("KBS_fungi_nonnorm_jan31_last.rda")
otu_table(KBS_fungi_nonnorm) <-
  otu_table(KBS_fungi_nonnorm)[which(rowSums(otu_table(KBS_fungi_nonnorm)) > 0), ]

# Generate working objects ----------------------------------------------------------------------------------
decont_all_lib -> physeq_all
physeq_all

KBS_fungi_nonnorm -> physeq_kel
physeq_kel

# Add taxonomy ----------------------------------------------------------------------------------------------
# load Taxonomy function
source("../R_functions/ReformatTaxonomy.R")

taxonomy_ITS07 <-
  read.delim(
    "constax_taxonomy_07.txt",
    header = TRUE,
    row.names = 1,
    sep = "\t"
  )

taxonomy_ITS07 = as.data.frame(taxonomy_ITS07)
taxonomy_ITS07 <-
  data.frame(OTU_ID = rownames(taxonomy_ITS07), taxonomy_ITS07)
head(taxonomy_ITS07)
str(taxonomy_ITS07)

# add NA
taxonomy_ITS07[taxonomy_ITS07==""] <- NA
any(taxonomy_ITS07=="")

taxonomy_ITS07$Species <- gsub(" sp ", "", taxonomy_ITS07$Species)
any(grepl(" sp ", taxonomy_ITS07$Species)==TRUE)

# change factor to characters
taxonomy_ITS07 <- 
  taxonomy_ITS07 %>% mutate_if(is.factor, as.character) 

taxonomy_ITS07[which(is.na(taxonomy_ITS07$Genus) == FALSE),]$Genus <-
  paste(taxonomy_ITS07$Genus[is.na(taxonomy_ITS07$Genus) == FALSE], "sp.", sep = " ")

rownames(taxonomy_ITS07) <- taxonomy_ITS07$OTU_ID
taxonomy_ITS07[1:100,]

# Generate Taxonomy label -----------------------------------------------------------------------------------
FinalizeTaxonomy <- function(taxa_table) {
  taxa_table[] = lapply(taxa_table,
                        blank2na,
                        na.strings = c('', 'NA', 'na', 'N/A', 'n/a', 'NaN', 'nan'))
  lastValue <- function(x)
    tail(x[!is.na(x)], 1)
  last_taxons <- apply(taxa_table[, 1:8], 1, lastValue)
  taxa_table$BestMatch <- last_taxons
  taxa_table[, "BestMatch"] <-
    gsub("_", " ", taxa_table[, "BestMatch"])
  taxa_table$Taxonomy <-
    paste(taxa_table$OTU_ID, taxa_table$BestMatch, sep = "-")
  taxa_table[, "Genus"] <- gsub(" sp.", "", taxa_table[, "Genus"])
  return(taxa_table)
}

FinalizeTaxonomy(taxonomy_ITS07) -> constax_taxonomy_07
dim(constax_taxonomy_07)
head(constax_taxonomy_07)
constax_taxonomy_07

# save the file for Acer
save(constax_taxonomy_07, file = "constax_taxonomy_07.RData")

# Filter out non-target taxa --------------------------------------------------------------------------------
str(constax_taxonomy_07)
table(constax_taxonomy_07$High_level_taxonomy)
#Eukaryota_kgd_Incertae_sedis  Fungi  Metazoa  MockK   Protista  Viridiplantae 
#             1                4924      3       13        3           608 

constax_taxonomy_07_filt <-
  subset(
    constax_taxonomy_07,
    High_level_taxonomy == "Fungi" &
      HL_hit_percent_id >= 60 &
      HL_hit_query_cover >= 85
  )
dim(constax_taxonomy_07_filt)
dim(constax_taxonomy_07)

dim(constax_taxonomy_07[constax_taxonomy_07$High_level_taxonomy == "Fungi" &
                          constax_taxonomy_07$HL_hit_query_cover >= 85 &
                          constax_taxonomy_07$HL_hit_percent_id >= 60, ])

dim(subset(constax_taxonomy_07, High_level_taxonomy == "Fungi"))
dim(subset(
  constax_taxonomy_07,
  High_level_taxonomy == "Fungi" &
    HL_hit_percent_id < 60
))
dim(subset(
  constax_taxonomy_07,
  High_level_taxonomy == "Fungi" &
    HL_hit_query_cover < 85
))

# Import otu sequences 
otus_seq <-
  readDNAStringSet(
    "otus_all_lib_R1.fasta",
    format = "fasta",
    seek.first.rec = TRUE,
    use.names = TRUE)

# Add to phyloseq object 
physeq_all <- phyloseq(
  otu_table(physeq_all, taxa_are_rows = TRUE),
  sample_data(physeq_all),
  tax_table(as.matrix(constax_taxonomy_07_filt)),
  otus_seq
) 

physeq_all
head(sample_data(physeq_all))
head(tax_table(physeq_all))
str(physeq_all)

physeq_kel <- phyloseq(
  otu_table(physeq_kel, taxa_are_rows = TRUE),
  sample_data(physeq_kel),
  tax_table(as.matrix(constax_taxonomy_07_filt)),
  otus_seq
)

physeq_kel
head(sample_data(physeq_kel))
head(tax_table(physeq_kel))
str(physeq_kel)


# Quantifying Viridiplantae
constax_taxonomy_07_plants <-
  subset(constax_taxonomy_07,
       High_level_taxonomy == "Viridiplantae")


physeq_plants <- phyloseq(
  otu_table(physeq_all, taxa_are_rows = TRUE),
  sample_data(physeq_all),
  tax_table(as.matrix(constax_taxonomy_07_plants)),
  otus_seq
) 

write.dna(file="switchgrass_reads.fasta",
          refseq(physeq_plants), colsep="",
          format = "fasta")

# **********************************************************************-------------------------------------
# KBS DATASET -----------------------------------------------------------------------------------------------

head(sample_data(physeq_kel))
levels(as.factor(sample_data(physeq_kel)$Date))
levels(as.factor(sample_data(physeq_kel)$Subpopulations))
levels(as.factor(sample_data(physeq_kel)$Eco))

count(as.data.frame(as.matrix(sample_data(physeq_kel))), vars = Date)
count(as.data.frame(as.matrix(sample_data(physeq_kel))), vars = Subpopulations)
count(as.data.frame(as.matrix(sample_data(physeq_kel))), vars = Eco)
count(as.data.frame(as.matrix(sample_data(physeq_kel))), vars = site)
count(as.data.frame(as.matrix(sample_data(physeq_kel))), vars = inf)
count(as.data.frame(as.matrix(sample_data(physeq_kel))), vars = ident)
count(as.data.frame(as.matrix(sample_data(physeq_kel))), vars = PLANT_ID)

# look for NA
sample_data(physeq_kel)[is.na(sample_data(physeq_kel))] <-
  "Intermediate"
sample_data(physeq_kel)$Subpopulations
any(is.na(sample_data(physeq_kel)))

# ****************************************************************************-------------------------------
# INDICATOR SPECIES ANALYSIS --------------------------------------------------------------------------------
# Focal site KBS --------------------------------------------------------------------------------------------
GetIndicators <-function(dataframe, var){
  require(phyloseq) 
  require(indicspecies)
  require(dplyr)
  otu <- as(otu_table(dataframe), "matrix")
  otu <- as.data.frame(otu)
  metadata = as(sample_data(dataframe), "data.frame")
  taxa <- as(tax_table(dataframe), "matrix")
  taxa <- as.data.frame(taxa)
  multipatt <- multipatt(t(otu), metadata[,var], func = "IndVal.g",
                         control=how(nperm=999), duleg=TRUE)
  multipatt -> multipatt_fdr
  multipatt_fdr$sign$p.value <- p.adjust(multipatt_fdr$sign$p.value, "fdr")
  multipatt_fdr$sign[which(
    multipatt_fdr$sign$p.value <= 0.05), ] -> indicator_taxa
  taxa$OTU <- rownames(taxa)
  data.frame(OTU = as.factor(row.names(indicator_taxa)), indicator_taxa) %>%
    dplyr::left_join(taxa, by="OTU") -> indicator_taxa
  rownames(indicator_taxa) <- indicator_taxa$OTU
  indicator_taxa <- arrange(indicator_taxa, desc(stat))
  return(indicator_taxa)
}

# indicator p.value <=0.05 after fdr correction
GetIndicators(physeq_kel, "Date") -> ind_fungi_date
subset(ind_fungi_date, ind_fungi_date$stat >= 0.75) -> ind_fungi_date
head(ind_fungi_date)
dim(ind_fungi_date)

GetIndicators(physeq_kel, "Subpopulations") -> ind_fungi_sub
subset(ind_fungi_sub, ind_fungi_date$stat >= 0.75) -> ind_fungi_sub
head(ind_fungi_sub)
dim(ind_fungi_sub)

GetIndicators(physeq_kel, "inf") -> ind_fungi_inf
subset(ind_fungi_inf, ind_fungi_inf$stat >= 0.75) -> ind_fungi_inf
head(ind_fungi_inf)
dim(ind_fungi_inf)

GetIndicators(physeq_kel, "Eco") -> ind_fungi_eco
subset(ind_fungi_eco, ind_fungi_date$stat >= 0.75) -> ind_fungi_eco
head(ind_fungi_eco)
dim(ind_fungi_eco)

save(ind_fungi_inf, file = "indicators.RData")

# Site dataset ----------------------------------------------------------------------------------------------
readRDS("phyloseq_shared.rds") -> physeq_shared
sample_shared <-
  rownames(sample_data(physeq_shared))

head(sample_data(physeq_all))
levels(sample_data(physeq_all)$site)

physeq_pruned <-
  prune_samples(sample_shared, physeq_all)


GetIndicators(physeq_pruned, "Date") -> ind_site_fungi_date
subset(ind_site_fungi_date, ind_site_fungi_date$stat >= 0.75) -> ind_site_fungi_date
head(ind_site_fungi_date)
dim(ind_site_fungi_date)

GetIndicators(physeq_pruned, "Subpopulations") -> ind_site_fungi_sub
subset(ind_site_fungi_sub, ind_site_fungi_date$stat >= 0.75) -> ind_site_fungi_sub
head(ind_site_fungi_sub)
dim(ind_site_fungi_sub)

GetIndicators(physeq_pruned, "inf") -> ind_site_fungi_inf
subset(ind_site_fungi_inf, ind_site_fungi_inf$stat >= 0.75) -> ind_site_fungi_inf
head(ind_site_fungi_inf)
dim(ind_site_fungi_inf)

GetIndicators(physeq_pruned, "Eco") -> ind_site_fungi_eco # No indicators
subset(ind_site_fungi_eco, ind_site_fungi_date$stat >= 0.75) -> ind_site_fungi_eco
head(ind_site_fungi_eco)
dim(ind_site_fungi_eco)


# ***********************************************************------------------------------------------------
# EXTRACT CORE ----------------------------------------------------------------------------------------------
source("../R_functions/ExtarctCore.R")

ExtractCore(
  physeq_kel_filt,
  Var = Date,
  method = "increase",
  Group = NULL,
  Level = NULL
) -> otu_core_date

otu_core_date[1]

ExtractCore(
  physeq_kel_filt,
  Var = Subpopulations,
  method = "increase",
  Group = NULL,
  Level = NULL
) -> otu_core_sub

otu_core_sub[1]

union(otu_core_date[[1]], otu_core_sub[[1]]) -> core_all

# Core taxonomy
filterTaxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

taxa_kel_filt <- 
  as(tax_table(filterTaxa(physeq_kel_filt, core_all)), "matrix")
taxa_kel_filt <- 
  as.data.frame(taxa_kel_filt)
taxa_kel_filt

# Genus
taxa_kel_filt<- 
  as.data.frame(table(taxa_kel_filt$Genus))
taxa_kel_filt$Percent <-
  taxa_kel_filt$Freq/sum(taxa_kel_filt$Freq)*100
grid.table(taxa_kel_filt)

# Class
taxa_kel_filt<- 
  as.data.frame(table(taxa_kel_filt$Class))
taxa_kel_filt$Percent <-
  taxa_kel_filt$Freq/sum(taxa_kel_filt$Freq)*100
grid.table(taxa_kel_filt)

# NEUTRAL MODELS---------------------------------------------------------------------------------------------
FitNeutral(otu_core_date) -> nfit_otu_core_date
FitNeutral(otu_core_sub) -> nfit_otu_core_sub

# plot increase
PlotBCincrease(otu_core_date, 150) + labs(title = "Poplar fungi Root")
PlotBCincrease(otu_core_sub, 150) + labs(title = "Poplar fungi Soil")

# plotting neutral models
ggarrange(PlotNeutral(nfit_otu_core_date) + labs(title = "Fungi Root"),
          PlotNeutral(nfit_otu_core_sub) + labs(title = "Fungi Soil "),
          ncol = 4, 
          nrow = 1)


# Extracting core members at each time point ----------------------------------------------------------------
SubsetData <- function(physeq, core, Var){
  sub_set <- subset(sample_data(physeq), Date == Var)
  physeq_filt <- merge_phyloseq(otu_table(physeq),
                                tax_table(physeq),
                                refseq(physeq),
                                sub_set)
  otu_table(physeq_filt) <-
    otu_table(physeq_filt)[which(rowSums(otu_table(physeq_filt)) > 0),]
  physeq_filt <-
    prune_taxa(core, physeq_filt)
  return(physeq_filt)
}

SubsetData(physeq_kel_filt, core_all, "158") -> physeq_kel_158
SubsetData(physeq_kel_filt, core_all, "212") -> physeq_kel_212
SubsetData(physeq_kel_filt, core_all, "233") -> physeq_kel_233
SubsetData(physeq_kel_filt, core_all, "260") -> physeq_kel_260
SubsetData(physeq_kel_filt, core_all, "286") -> physeq_kel_286


# *******************************************************************************----------------------------
# MICROBIAL NETWORKS ----------------------------------------------------------------------------------------
source("../R_functions/MicrobNet.R")

MakeSENet(physeq_kel_158) -> spiec_physeq_kel_158 #100, 1e-1, thresh = 0.05
MakeSENet(physeq_kel_212) -> spiec_physeq_kel_212 
MakeSENet(physeq_kel_233) -> spiec_physeq_kel_233
MakeSENet(physeq_kel_260) -> spiec_physeq_kel_260
MakeSENet(physeq_kel_286) -> spiec_physeq_kel_286

getStability(spiec_physeq_kel_158)
getStability(spiec_physeq_kel_212)
getStability(spiec_physeq_kel_233)
getStability(spiec_physeq_kel_260)
getStability(spiec_physeq_kel_286)

plot(adj2igraph(getRefit(spiec_physeq_kel_158)), 
     vertex.color=c(rep(1, ntaxa(physeq_kel_158))) +1,
     vertex.size=9)

plot(adj2igraph(getRefit(spiec_physeq_kel_212)), 
     vertex.color=c(rep(1,ntaxa(physeq_kel_212))) +1,
     vertex.size=9)

plot(adj2igraph(getRefit(spiec_physeq_kel_233)), 
     vertex.color=c(rep(1,ntaxa(physeq_kel_233))) +1,
     vertex.size=9)

plot(adj2igraph(getRefit(spiec_physeq_kel_260)), 
     vertex.color=c(rep(1,ntaxa(physeq_kel_260))) +1,
     vertex.size=9)

plot(adj2igraph(getRefit(spiec_physeq_kel_286)), 
     vertex.color=c(rep(1,ntaxa(physeq_kel_286))) +1,
     vertex.size=9)


# Make network ----------------------------------------------------------------------------------------------
GetNetwork(physeq_kel_158, NULL, spiec_physeq_kel_158) -> network_kel_158
plot_network(network_kel_158)
GetNetwork(physeq_kel_212, NULL, spiec_physeq_kel_212) -> network_kel_212
plot_network(network_kel_212)
GetNetwork(physeq_kel_233, NULL, spiec_physeq_kel_233) -> network_kel_233
plot_network(network_kel_233)
GetNetwork(physeq_kel_260, NULL, spiec_physeq_kel_260) -> network_kel_260
plot_network(network_kel_260)
GetNetwork(physeq_kel_286, NULL, spiec_physeq_kel_286) -> network_kel_286
plot_network(network_kel_286)


# Test random -----------------------------------------------------------------------------------------------
TesToNULL(network_kel_158)
TesToNULL(network_kel_212)
TesToNULL(network_kel_233)
TesToNULL(network_kel_260)
TesToNULL(network_kel_286)

TesToRandom(network_kel_158, 110)
TesToRandom(network_kel_212, 110)
TesToRandom(network_kel_233, 110)
TesToRandom(network_kel_260, 110)
TesToRandom(network_kel_286, 110)

# Calculate abundances --------------------------------------------------------------------------------------
AbundTotal <- function(physeq) {
  # Extracting abudnances and taxonomy
  physeq %>%
    transform_sample_counts(function(x) {x / sum(x) *100}) %>% 
    taxa_sums() %>% 
    as.data.frame() -> otu
  colnames(otu) <- "TotAbund"
  otu <-
    data.frame(OTU_ID = as.factor(row.names(otu)), otu)
  dim(otu) %>% print()
  # calclating abundances for groups
  # Infection
  physeq %>%
    transform_sample_counts(function(x) {x / sum(x) *100}) -> physeq_ab
  dim(otu_table(physeq_ab)) %T>% print()
  if (any(sample_data(physeq_ab)$inf=="I")){
    otu$InfAbund <- physeq_ab %>%
      subset_samples(inf %in% c("I")) %>%
      taxa_sums()
  } else{
    otu$InfAbund <- 0
  }
  #otu$InfAbund %T>% print()
  otu$HealthyAbund <- physeq_ab %>%
    subset_samples(inf %in% c("U")) %>%
    taxa_sums()
  #Subpopulations
  otu$Atl <- physeq_ab %>%
    subset_samples(Subpopulations %in% c("Atlantic")) %>%
    taxa_sums()
  otu$Gul <- physeq_ab %>%
    subset_samples(Subpopulations %in% c("Gulf")) %>%
    taxa_sums()
  otu$Int <- physeq_ab %>%
    subset_samples(Subpopulations %in% c("Intermediate")) %>%
    taxa_sums()
  otu$Mid <- physeq_ab %>%
    subset_samples(Subpopulations %in% c("Midwest")) %>%
    taxa_sums()
  print(head(otu))
  # # Date
  # otu$Abund158 <- physeq %>%
  #   subset_samples(Date %in% c("158")) %>%
  #   taxa_sums()
  # otu$Abund212 <- physeq %>%
  #   subset_samples(Date %in% c("212")) %>%
  #   taxa_sums()
  # otu$Abund233 <- physeq %>%
  #   subset_samples(Date %in% c("233")) %>%
  #   taxa_sums()
  # otu$Abund260 <- physeq %>%
  #   subset_samples(Date %in% c("260")) %>%
  #   taxa_sums()
  # otu$Abund286 <- physeq %>%
  #   subset_samples(Date %in% c("286")) %>%
  #   taxa_sums()
  # converting in % of the TotAbund
  otu$InfAbund <- otu$InfAbund/otu$TotAbund *100
  otu$HealthyAbund <-  otu$HealthyAbund/otu$TotAbund *100
  otu$Atl <-  otu$Atl/otu$TotAbund*100
  otu$Gul <- otu$Gul/otu$TotAbund*100
  otu$Int <- otu$Int/otu$TotAbund*100
  otu$Mid <- otu$Mid/otu$TotAbund*100
  # otu$Abund158 <- otu$Abund158/otu$TotAbund*100
  # otu$Abund212 <- otu$Abund212/otu$TotAbund*100
  # otu$Abund233 <- otu$Abund233/otu$TotAbund*100
  # otu$Abund260 <- otu$Abund260/otu$TotAbund*100
  # otu$Abund286 <- otu$Abund286/otu$TotAbund*100
  return(otu)
}

AbundTotal(physeq_kel_158) -> abund_kel_all_158
AbundTotal(physeq_kel_212) -> abund_kel_all_212
AbundTotal(physeq_kel_233) -> abund_kel_all_233
AbundTotal(physeq_kel_260) -> abund_kel_all_260
AbundTotal(physeq_kel_286) -> abund_kel_all_286

# calculate Node attribtues ---------------------------------------------------------------------------------
CalcNetAttrib <- function(network, Var) {
  require(igraph)
  require(ggrepel)
  require(microbiome)
  nodes <-
    data.frame(matrix(ncol = 8, nrow = length(V(network)$name)))
  colnames(nodes) <-
    c("Degree",
      "Modules",
      "Betweenness",
      "Date_ind",
      "Infection_ind",
      "Subpopulation_ind",
      "Eco_ind",
      "Hubs")
  rownames(nodes) <- V(network)$name
  nodes <- data.frame(OTU_ID = as.factor(row.names(nodes)), nodes)
  dim(nodes) %T>% print()
  components(network) %T>% print()
  #Degree
  nodes$Degree <- igraph::degree(network, mode = "all")
  # removing weights
  network_new <- delete_edge_attr(network, "weight") # delete weights
  #Moduels
  #modules <- cluster_walktrap(network_new)
  modules <- cluster_fast_greedy(network_new)
  nodes$Modules <- as.factor(igraph::membership(modules))
  # Betweenness
  nodes$Betweenness <- igraph::betweenness(network_new, directed = FALSE)
  # Indicators inf
  ind_fungi_inf[rownames(ind_fungi_inf) %in% rownames(nodes), ] -> net_ind_fungi_inf
  if (any(net_ind_fungi_inf$s.I=="1") ==TRUE){ 
    nodes[rownames(nodes) %in%
            rownames(net_ind_fungi_inf[net_ind_fungi_inf$s.I =="1", ]), ]$Infection_ind <- "Inf"
  }else{}
  if (any(net_ind_fungi_inf$s.U=="1") ==TRUE){ 
    nodes[rownames(nodes) %in%
            rownames(net_ind_fungi_inf[net_ind_fungi_inf$s.U =="1", ]), ]$Infection_ind <- "Un"
  }else{}
  # Indicators date
  ind_fungi_date[rownames(ind_fungi_date) %in% rownames(nodes), ] -> net_ind_fungi_date
  if (any(net_ind_fungi_date$s.158 == "1") == TRUE) {
    nodes[rownames(nodes) %in%
            rownames(net_ind_fungi_date[net_ind_fungi_date$s.158 == "1", ]), ]$Date_ind <- "d158"
  }else{}
  if (any(net_ind_fungi_date$s.212 == "1") == TRUE) {
    nodes[rownames(nodes) %in%
            rownames(net_ind_fungi_date[net_ind_fungi_date$s.212 == "1", ]), ]$Date_ind <- "d212"
  }else{}
  if (any(net_ind_fungi_date$s.233 == "1") == TRUE) {
    nodes[rownames(nodes) %in%
            rownames(net_ind_fungi_date[net_ind_fungi_date$s.233 == "1", ]), ]$Date_ind <- "d233"
  }else{}
  if (any(net_ind_fungi_date$s.260 == "1") == TRUE) {
    nodes[rownames(nodes) %in%
            rownames(net_ind_fungi_date[net_ind_fungi_date$s.260 == "1", ]), ]$Date_ind <- "d260"
  }else{}
  if (any(net_ind_fungi_date$s.286 == "1") == TRUE) {
    nodes[rownames(nodes) %in%
            rownames(net_ind_fungi_date[net_ind_fungi_date$s.286 == "1", ]), ]$Date_ind <- "d286"
  }else{}
  #Subpopulation indicators
  ind_fungi_sub[rownames(ind_fungi_sub) %in% rownames(nodes), ] -> net_ind_fungi_sub
  if (any(net_ind_fungi_sub$s.I=="1") ==TRUE){ 
    nodes[rownames(nodes) %in%
            rownames(net_ind_fungi_sub[net_ind_fungi_sub$s.Atlantic == "1", ]), ]$Subpopulation_ind <- "Atlantic"
  }else{}
  if (any(net_ind_fungi_sub$s.I=="1") ==TRUE){ 
    nodes[rownames(nodes) %in%
            rownames(net_ind_fungi_sub[net_ind_fungi_sub$s.Gulf == "1", ]), ]$Subpopulation_ind <- "Gulf"
  }else{}
  if (any(net_ind_fungi_sub$s.I=="1") ==TRUE){ 
    nodes[rownames(nodes) %in%
            rownames(net_ind_fungi_sub[net_ind_fungi_sub$s.Intermediate == "1", ]), ]$Subpopulation_ind <- "Intermediate"
  }else{}
  if (any(net_ind_fungi_sub$s.I=="1") ==TRUE){ 
    nodes[rownames(nodes) %in%
            rownames(net_ind_fungi_sub[net_ind_fungi_sub$s.Midwest == "1", ]), ]$Subpopulation_ind <- "Midwest"
  }else{}
  # Ecotype indicators
  ind_fungi_sub[rownames(ind_fungi_eco) %in% rownames(nodes), ] -> net_ind_fungi_eco
  if (any(net_ind_fungi_eco$s.eastcoast=="1") ==TRUE){ 
    nodes[rownames(nodes) %in%
            rownames(net_ind_fungi_eco[net_ind_fungi_eco$s.eastcoast == "1", ]), ]$Eco_ind <- "Eastcoats"
  }else{}
  if (any(net_ind_fungi_eco$s.Gulfcoast=="1") ==TRUE){ 
    nodes[rownames(nodes) %in%
            rownames(net_ind_fungi_eco[net_ind_fungi_eco$s.Gulfcoast == "1", ]), ]$Eco_ind <- "Gulfcoast"
  }else{}
  if (any(net_ind_fungi_eco$s.midwest=="1") ==TRUE){ 
    nodes[rownames(nodes) %in%
            rownames(net_ind_fungi_eco[net_ind_fungi_eco$s.midwest == "1", ]), ]$Eco_ind <- "Midwest"
  }else{}
  if (any(net_ind_fungi_eco$s.Texas=="1") ==TRUE){ 
    nodes[rownames(nodes) %in%
            rownames(net_ind_fungi_eco[net_ind_fungi_eco$s.Texas == "1", ]), ]$Eco_ind <- "Texas"
  }else{}
  #Hubs
  df_zipi <- ZiPi(network_new, modules=nodes$Modules)
  df_zipi %T>% assign(paste("df_zipi", Var, sep = "_"),., envir = .GlobalEnv) # saving the df
  df_zipi$Key  <- ifelse(df_zipi$P>=0.62, "Connectors", NA)
  df_zipi$Key  <- ifelse(df_zipi$Z>=2.5, "Module hubs", paste(df_zipi$Key))
  df_zipi$Key  <- ifelse(df_zipi$P>=0.62 & df_zipi$Z>=2.5, "Network hubs", paste(df_zipi$Key))
  df_zipi$Key  <- ifelse(df_zipi$P<0.62 & df_zipi$Z<2.5, "Peripherals", paste(df_zipi$Key))
  df_zipi$Key  <- ifelse(is.na(df_zipi$Key), "Peripherals", paste(df_zipi$Key))
  dim(df_zipi) %T>% print()
  if (identical(as.character(nodes$OTU_ID), as.character(df_zipi$names))) {
    print("dataframes order is identical!")
    nodes$Key <- df_zipi$Key
  } else {
    print("dataframes order is not identical! Reordering.")
    sample_order <- match(as.character(df_zipi$names), as.character(nodes$OTU_ID))
    df_zipi <- df_zipi[sample_order, ]
    nodes$Key <- df_zipi$Key
  }
  return(nodes)
}

# Test
CalcNetAttrib(net_kel_158, "158") -> Attr_kel_158
Attr_kel_158
df_zipi_158


# Calculate node table ---------------------------------------------------------------------------------------
CalcNodeTable <-
  function(physeq,abund,network,
           neutral_fit_f = NULL,neutral_fit_b = NULL,
           Var) {
    CalcNetAttrib(network, Var)  -> nodes
    inner_join(abund, nodes, by = "OTU_ID") -> nodes_abund
    # extarcting taxa
    taxa <- as(tax_table(physeq), "matrix")
    taxa <- as.data.frame(taxa)
    inner_join(nodes_abund, taxa, by = "OTU_ID") -> nodes_attr
    #list(as_tibble(nodes), as_tibble(abund)) -> nodes_attr
    rownames(nodes_attr) <- nodes_attr$OTU_ID
    if ((!is.null(neutral_fit_f) | !is.null(neutral_fit_b))) {
      neutral_fit <-
        rbind(neutral_fit_f[[2]], neutral_fit_b[[2]])
      if (identical(neutral_fit$OTU_ID, rownames(nodes_attr))) {
        
      } else{
        sample_order <- match(rownames(nodes_attr), neutral_fit$OTU_ID)
        neutral_fit <- neutral_fit[sample_order, ]
        nodes_attr$fit_class <- neutral_fit$fit_class
      }
    }
    return(nodes_attr)
  }


CalcNodeTable(
  physeq_kel_158,
  abund_kel_all_158,
  network_kel_158,
  nfit_otu_core_date,
  nfit_otu_core_sub,
  "158"
) -> nodes_kel_all_158
nodes_kel_all_158$Class <-
  ifelse(is.na(nodes_kel_all_158$Class),
         "Unclassified",
         paste(nodes_kel_all_158$Class))

CalcNodeTable(
  physeq_kel_212,
  abund_kel_all_212,
  network_kel_212,
  nfit_otu_core_date,
  nfit_otu_core_sub,
  "212"
) -> nodes_kel_all_212
nodes_kel_all_212$Class <-
  ifelse(is.na(nodes_kel_all_212$Class),
         "Unclassified",
         paste(nodes_kel_all_212$Class))

CalcNodeTable(
  physeq_kel_233,
  abund_kel_all_233,
  network_kel_233,
  nfit_otu_core_date,
  nfit_otu_core_sub,
  "233"
) -> nodes_kel_all_233
nodes_kel_all_233$Class <-
  ifelse(is.na(nodes_kel_all_233$Class),
         "Unclassified",
         paste(nodes_kel_all_233$Class))

CalcNodeTable(
  physeq_kel_260,
  abund_kel_all_260,
  network_kel_260,
  nfit_otu_core_date,
  nfit_otu_core_sub,
  "260"
) -> nodes_kel_all_260
nodes_kel_all_260$Class <-
  ifelse(is.na(nodes_kel_all_260$Class),
         "Unclassified",
         paste(nodes_kel_all_260$Class))

CalcNodeTable(
  physeq_kel_286,
  abund_kel_all_286,
  network_kel_286,
  nfit_otu_core_date,
  nfit_otu_core_sub,
  "286"
) -> nodes_kel_all_286
nodes_kel_all_286$Class <-
  ifelse(is.na(nodes_kel_all_286$Class),
         "Unclassified",
         paste(nodes_kel_all_286$Class))


# Edge attributes -------------------------------------------------------------------------------------------
CalcEdges(network_kel_158, nodes_kel_all_158) -> edges_kel_all_158
CalcEdges(network_kel_212, nodes_kel_all_212) -> edges_kel_all_212
CalcEdges(network_kel_233, nodes_kel_all_233) -> edges_kel_all_233
CalcEdges(network_kel_260, nodes_kel_all_260) -> edges_kel_all_260
CalcEdges(network_kel_286, nodes_kel_all_286) -> edges_kel_all_286

# NETWORK STATS SUMMARY -------------------------------------------------------------------------------------
cbind(
  NetStats(nodes_kel_all_158, edges_kel_all_158, network_kel_158, "d158"),
  NetStats(nodes_kel_all_212, edges_kel_all_212, network_kel_212, "d212"),
  NetStats(nodes_kel_all_233, edges_kel_all_233, network_kel_233, "d233"),
  NetStats(nodes_kel_all_260, edges_kel_all_260, network_kel_260, "d260"),
  NetStats(nodes_kel_all_286, edges_kel_all_286, network_kel_286, "d286")
) -> network_kel_all_summary

network_kel_all_summary

#network stats
grid.table(network_kel_all_summary)
grid.table(network_kel_all_summary[c(1:5, 17, 19:21,24),])
grid.table(network_kel_all_summary[c(1,2,5,6,7,11,16),])

write.csv(network_kel_all_summary, file = "network_kel_all_summary.csv")

# PLOTTING THE NETWORK  -------------------------------------------------------------------------------------

MakeTabGraph(nodes_kel_all_158, edges_kel_all_158) -> ggraph_kel_all_158
MakeTabGraph(nodes_kel_all_212, edges_kel_all_212) -> ggraph_kel_all_212
MakeTabGraph(nodes_kel_all_233, edges_kel_all_233) -> ggraph_kel_all_233
MakeTabGraph(nodes_kel_all_260, edges_kel_all_260) -> ggraph_kel_all_260
MakeTabGraph(nodes_kel_all_286, edges_kel_all_286) -> ggraph_kel_all_286

# Classes
levels(unique(as.factor(nodes_kel_all_158$Class),
              as.factor(nodes_kel_all_212$Class),
              as.factor(nodes_kel_all_233$Class),
              as.factor(nodes_kel_all_260$Class), 
              as.factor(nodes_kel_all_286$Class)))

# NETWORK TAXONOMY ------------------------------------------------------------------------------------------
PlotNetw <-function(tab_raph, Var, title){
  set.seed(21922)
  plot <-
    tab_raph %>%
    ggraph(layout = "fr") +
    geom_edge_link(aes(edge_colour = Direction), edge_width = 0.05) +
    geom_node_point(aes(color = Class, fill = Class, shape = Key, size = sqrt(get(Var)))) + # log(abundance)
    # geom_node_text(aes(filter = Date_ind, label = OTU_ID),
    #                repel=TRUE, color="black",  size =2.5) +
    #facet_edges(~InterKing+Direction) +
    scale_edge_color_manual(values = c("red", "black")) + 
    scale_edge_linetype_manual(values = c("solid", "dashed")) +
    scale_edge_width(range = c(0.005, 0.5))+
    scale_size(range = c(0.01, 3)) +
    scale_shape_manual(values = c(Connectors = 18, 
                                  `Module hubs` = 15, 
                                  `Network hubs` = 17, 
                                  Peripherals = 16)) +
    # Note! If the color is set to taxonomy should specify each color 
    # for each taxa or they are going to be different acorss networks.
    scale_color_manual(values = c(Agaricomycetes = "#F0E442", 
                                  Agaricostilbomycetes = "#D55E00", 
                                  Cystobasidiomycetes = "#AD6F3B", 
                                  Dothideomycetes = "#E69F00", 
                                  Exobasidiomycetes = "#fff0d6",
                                  Leotiomycetes = "#009E73",
                                  Microbotryomycetes = "#56B4E9", 
                                  Sordariomycetes = "#CC79A7", 
                                  Taphrinomycetes = "#8569D5",
                                  Tremellomycetes = "#808080", 
                                  Unclassified = "black")) + #na.value = "black") +
    theme_graph() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.margin = unit(c(0.4,0.01,0.4,0.01), "cm"), 
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 8, face = "bold", hjust = 0, vjust = 0.5), 
          legend.text = element_text(size = 7), legend.position = "bottom") +
    labs(title = title) +
    guides(size = FALSE, 
           shape = guide_legend(ncol = 2, order = 3),
           fill = guide_legend(ncol = 5),
           colour = guide_legend(ncol=5),
           edge_colour = guide_legend(ncol = 1, order = 1), 
           edge_width = FALSE)
  return(plot)
}


legend_net_tax <-
  get_legend(PlotNetw(ggraph_kel_all_233, "TotAbund", "d233"))
as_ggplot(legend_net_tax)


# *** NETWORK plot ------------------------------------------------------------------------------------------
net_graph <- 
  ggarrange(
    ggarrange(
      PlotNetw(ggraph_kel_all_158, "TotAbund", "d158") +
        geom_node_point(aes(filter = Date_ind%in%"d158"), 
                        shape = 1, size = 3, show.legend = F) +
        geom_node_text(aes(filter = Date_ind%in%"d158", label = OTU_ID),
                       repel=TRUE, color="black",  size =2),
      PlotNetw(ggraph_kel_all_212, "TotAbund", "d212"),
      PlotNetw(ggraph_kel_all_233, "TotAbund", "d233"),
      PlotNetw(ggraph_kel_all_260, "TotAbund", "d260") +
        geom_node_point(aes(filter = Date_ind%in%"d260"), 
                        shape = 1, size = 3, show.legend = F) +
        geom_node_text(aes(filter = Date_ind%in%"d260", label = OTU_ID),
                       repel=TRUE, color="black",  size =2),
      PlotNetw(ggraph_kel_all_286, "TotAbund", "d286")+
        geom_node_point(aes(filter = Date_ind%in%"d286"), 
                        shape = 1, size = 3, show.legend = F) +
        geom_node_text(aes(filter = Date_ind%in%"d286", label = OTU_ID),
                       repel=TRUE, color="black",  size =2),
      ncol = 5,
      nrow = 1,
      #labels = c("A", "B", "C", "D"),
      widths = c(1, 1, 1, 1, 1),
      common.legend = TRUE,
      legend = "none"),
    as_ggplot(legend_net_tax),
    ncol = 1,
    nrow = 2,
    heights = c(1, 0.15))


net_graph
str(net_graph)

# CLASS BARPLOT ----------------------------------------------------------------------------------------------
palette_class <- c("#F0E442", "#D55E00",  "#AD6F3B", "#E69F00", "#fff0d6",
                   "#009E73", "#56B4E9","#CC79A7", "#8569D5","#808080", "black")

PlotBar <- function(nodes){
  df <- nodes[, c("OTU_ID", "TotAbund", "Class")]
  df$RelAb <-
    df$TotAbund/sum(df$TotAbund)*100
  plot_bars <- 
    ggplot(df, aes(x = Class, y = sqrt(TotAbund), color = Class, fill=Class)) + 
    geom_bar(stat = "identity") +
    #geom_boxplot() +
    ylim(0, 500) +
    theme_classic()+
    scale_fill_manual(values = palette_class) +
    scale_color_manual(values = palette_class) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), 
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title = element_text(angle = 0, size = 7, face = "bold"),
          axis.text.y = element_text(size = 7), 
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 8, face = "bold", hjust = 0, vjust = 0.5), 
          legend.text = element_text(size = 7), legend.position = "bottom") +
    #theme(axis.text.x = element_text(angle = 90, size = 7, vjust = 0.5, hjust = 1)) +
    guides(color=guide_legend(ncol=1)) 
  return(plot_bars)
}

# **** BARPLOT plot -----------------------------------------------------------------------------------------
barplot_time <-
  ggarrange(PlotBar(nodes_kel_all_158) + labs(title=NULL, x=NULL, y = "Cumulative Abundance"),
            PlotBar(nodes_kel_all_212) + labs(title=NULL, x=NULL, y = NULL),
            PlotBar(nodes_kel_all_233) + labs(title=NULL, x=NULL, y = NULL),
            PlotBar(nodes_kel_all_260) + labs(title=NULL, x=NULL, y = NULL),
            PlotBar(nodes_kel_all_286) + labs(title=NULL, x=NULL, y = NULL),
            ncol = 5,
            nrow = 1,
            align = "hv",
            #labels = c("A", "B", "C", "D"),
            #widths = c(1, 1, 1, 1),
            heights = c(1, 1, 1, 1),
            common.legend = TRUE,
            legend = "none")

barplot_time

# **** COMBINED plot -----------------------------------------------------------------------------------------
ggarrange(net_graph,
          barplot_time,
          ncol = 1,
          nrow = 2,
          heights = c(1, 0.5))

net_graph_all <- 
  ggarrange(
    ggarrange(
      PlotNetw(ggraph_kel_all_158, "TotAbund", "d158") +
        geom_node_point(aes(filter = Date_ind%in%"d158"), 
                        shape = 1, size = 3, show.legend = F) +
        geom_node_text(aes(filter = Date_ind%in%"d158", label = OTU_ID),
                       repel=TRUE, color="black",  size =2),
      PlotNetw(ggraph_kel_all_212, "TotAbund", "d212"),
      PlotNetw(ggraph_kel_all_233, "TotAbund", "d233"),
      PlotNetw(ggraph_kel_all_260, "TotAbund", "d260") +
        geom_node_point(aes(filter = Date_ind%in%"d260"), 
                        shape = 1, size = 3, show.legend = F) +
        geom_node_text(aes(filter = Date_ind%in%"d260", label = OTU_ID),
                       repel=TRUE, color="black",  size =2),
      PlotNetw(ggraph_kel_all_286, "TotAbund", "d286")+
        geom_node_point(aes(filter = Date_ind%in%"d286"), 
                        shape = 1, size = 3, show.legend = F) +
        geom_node_text(aes(filter = Date_ind%in%"d286", label = OTU_ID),
                       repel=TRUE, color="black",  size =2),
      ncol = 5,
      nrow = 1,
      #labels = c("A", "B", "C", "D"),
      widths = c(1, 1, 1, 1, 1),
      common.legend = TRUE,
      legend = "none"),
    barplot_time,
    as_ggplot(legend_net_tax),
    ncol = 1,
    nrow = 3,
    heights = c(1, 0.5, 0.2))

net_graph_all


# NETWORK INFECTION ABUNDANCE --------------------------------------------------------------------------------
PlotNetw <-function(tab_raph, Var, title){
  set.seed(21922)
  plot <-
    tab_raph %>%
    ggraph(layout = "fr") +
    geom_edge_link(aes(edge_colour = Direction), edge_width = 0.05) +
    geom_node_point(aes(color = InfAbund, fill = InfAbund, 
                        shape = Key, size = sqrt(get(Var)))) + # log(abundance)
    geom_node_text(aes(filter = Infection_ind%in%"Inf",
                       label = OTU_ID), repel=TRUE, color="black",  size =2) +
    # #facet_edges(~InterKing+Direction) +
    scale_edge_color_manual(values = c("red", "black")) + 
    scale_edge_linetype_manual(values = c("solid", "dashed")) +
    scale_edge_width(range = c(0.005, 0.5))+
    scale_size(range = c(0.01, 3)) +
    scale_shape_manual(values = c(Connectors = 18, 
                                  `Module hubs` = 15, 
                                  `Network hubs` = 17, 
                                  Peripherals = 16)) +
    theme_graph() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.margin = unit(c(0.4,0.01,0.4,0.01), "cm"), 
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 8, face = "bold", hjust = 0, vjust = 0.5), 
          legend.text = element_text(size = 7), legend.position = "bottom") +
    labs(title = title) +
    guides(size = FALSE, 
           shape = guide_legend(ncol = 2, order = 3),
           edge_colour = guide_legend(ncol = 1, order = 1), 
           edge_width = FALSE)
  return(plot)
}

legend_net_inf <-
  get_legend(PlotNetw(ggraph_kel_all_260, "TotAbund", "day 158"))
as_ggplot(legend_net_inf)

# **** NETWORK plot 
network_abund <-
  ggarrange(
    ggarrange(
      PlotNetw(ggraph_kel_all_158, "TotAbund", "day 158") +
        geom_node_point(aes(filter = Infection_ind%in%"Inf"), 
                        shape = 1, size = 3, show.legend = F),
      PlotNetw(ggraph_kel_all_212, "TotAbund", "day 212") +
        geom_node_point(aes(filter = Infection_ind%in%"Inf"), 
                        shape = 1, size = 3, show.legend = F),
      PlotNetw(ggraph_kel_all_233, "TotAbund", "day 233") +
        geom_node_point(aes(filter = Infection_ind%in%"Inf"), 
                        shape = 1, size = 3, show.legend = F),
      PlotNetw(ggraph_kel_all_260, "TotAbund", "day 260") +
        geom_node_point(aes(filter =Infection_ind%in%"Inf"), 
                        shape = 1, size = 3, show.legend = F),
      PlotNetw(ggraph_kel_all_286, "TotAbund", "day 286")+
        geom_node_point(aes(filter = Infection_ind%in%"Inf"), 
                        shape = 1, size = 3, show.legend = F),
      ncol = 5,
      nrow = 1,
      #labels = c("A", "B", "C", "D"),
      widths = c(1, 1, 1, 1, 1),
      common.legend = TRUE,
      legend = "none"),
    as_ggplot(legend_net_inf),
    ncol = 1,
    nrow = 2,
    heights = c(1, 0.1))

network_abund


# INVESTIGATING CLASS-to-CLASS connections -------------------------------------------------------------------
# Class level HEATMAP ----------------------------------------------------------------------------------------
HeatEdge <- function(df_edges, Var, Date){
  df_edges <-
    df_edges %>%
    mutate(ModuleConnectivity = ifelse(V1_Class==V2_Class, "in-class", "out-class"))
  df_edges <- 
    subset(df_edges, Direction%in%Var)
  df_count <-
    dplyr::count(df_edges, 
                 vars= Direction, V1_Class, V2_Class, ModuleConnectivity) 
  print("Subsetted dataframe")
  str(df_count) %T>% print()
  #df_count %T>% print()
  # df_count$vars <-
  #   as.factor(df_count$vars)
  # df_count_filt <-
  #   subset(df_count, vars%in%Var)
  df_count_filt <-
    df_count[, c(2,3,5)]
  print("Subsetted dataframe dimensions")
  dim(df_count_filt) %T>% print()
  # extract Class data
  df_count_filt$V1_Class <- as.character(df_count_filt$V1_Class)
  df_count_filt$V2_Class <- as.character(df_count_filt$V2_Class)
  faclist <- vector("list", nrow(df_count_filt))
  for (i in 1:nrow(df_count_filt)) {
    faclist[[i]] <-
      sort(c(df_count_filt$V1_Class[i], df_count_filt$V2_Class[i]))
    faclist[[i]] <- paste(faclist[[i]][1], faclist[[i]][2], sep = "_")
  }
  df_count_filt$v1v2 <- unlist(faclist)
  # print("Subsetted dataframe")
  # df_count_filt %T>% print()
  df_count_filt %>% 
    group_by(v1v2) %>% summarize(Sum=sum(n)) %>%
    as.data.frame() -> df_new
  # print("Final dataframe dimensions")
  # dim(df_new)  %T>% print()
  df_final <-
    left_join(df_new, df_count_filt[,c(1:2,4)], "v1v2")
  df_final <-
    df_final[!duplicated(df_final$v1v2), ]
  df_final$Date <- rep(Date, times=nrow(df_final))
  return(df_final)
}


HeatEdge(edges_kel_all_158, "positive", "d158")
HeatEdge(edges_kel_all_233, "positive", "d233")
HeatEdge(edges_kel_all_233, "negative", "d233")


df_bind_positive <-
  rbind(
    HeatEdge(edges_kel_all_158, "positive", "d158"),
    HeatEdge(edges_kel_all_212, "positive", "d212"),
    HeatEdge(edges_kel_all_233, "positive", "d233"),
    HeatEdge(edges_kel_all_260, "positive", "d260"),
    HeatEdge(edges_kel_all_286, "positive", "d286")) 

df_bind_negative <-
  rbind(
    HeatEdge(edges_kel_all_158, "negative", "d158"),
    HeatEdge(edges_kel_all_212, "negative", "d212"),
    HeatEdge(edges_kel_all_233, "negative", "d233"),
    HeatEdge(edges_kel_all_260, "negative", "d260"),
    HeatEdge(edges_kel_all_286, "negative", "d286")) 


df_bing_all <-
  rbind(data.frame(Weight = rep("positive", nrow(df_bind_positive)),
                   df_bind_positive),
        data.frame(Weight = rep("negative", nrow(df_bind_negative)), 
                   df_bind_negative)) 

df_bing_all

# Just Tremellomycetes and Dothideomycetes -------------------------------------------------------------------
TremDot_pos <-
  subset(df_bind_positive, V1_Class%in%c("Dothideomycetes", "Tremellomycetes") &
           V2_Class%in%c("Dothideomycetes", "Tremellomycetes"))

TremDot_neg <-
  subset(df_bind_negative, V1_Class%in%c("Dothideomycetes", "Tremellomycetes") &
           V2_Class%in%c("Dothideomycetes", "Tremellomycetes"))

df_TremDot <-
  rbind(data.frame(Weight = rep("positive", nrow(TremDot_pos)),
                   TremDot_pos),
        data.frame(Weight = rep("negative", nrow(TremDot_neg)), 
                   TremDot_neg)) 
df_TremDot <- df_TremDot[, c(1,6,3,4,5)]

df_TremDot$Edges <-
  c(rep(165, 3), rep(264, 3), rep(269, 3), rep(251, 3), rep(247, 3),
    rep(165, 3), rep(264, 3), rep(269, 3), rep(251, 3), rep(247, 3))

df_TremDot$Mean <- 
  round(df_TremDot$Sum/df_TremDot$Edges*100, digits = 1)

df_TremDot

# Calculating Class abundances -------------------------------------------------------------------------------
subset(df_TremDot, 
       V1_Class%in%c("Dothideomycetes", "Tremellomycetes") &
         V2_Class%in%c("Dothideomycetes", "Tremellomycetes")) %>%
  mutate(Connectivity = ifelse(V1_Class==V2_Class, "in-class", "out-class")) %>%
  group_by(Connectivity, Date) %>%
  summarise(Total = sum(Mean)) %>%
  as.data.frame() %>%
  grid.table()

subset(df_TremDot, 
       V1_Class%in%c("Dothideomycetes", "Tremellomycetes") &
         V2_Class%in%c("Dothideomycetes", "Tremellomycetes")) %>%
  mutate(Connectivity = ifelse(V1_Class==V2_Class, "in-class", "out-class")) %>%
  group_by(Connectivity, Weight, Date) %>%
  summarise(Total = sum(Mean)) %>%
  as.data.frame() %>%
  grid.table()


# Plotting Edge Heatmap --------------------------------------------------------------------------------------
PlotHeat <- function(dataframe){
  heat_plot <-
    ggplot(dataframe, aes(x=V1_Class, y=V2_Class, fill=Mean)) +
    geom_tile(aes(height = 0.92, width = 0.92)) + 
    geom_text(aes(label = Mean), size=2.5, color="white") +
    theme_classic() +
    facet_grid(~Weight~Date) +
    labs(title = "Within- and between-Class connections") +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), 
          axis.text.x =  element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 8), 
          axis.title = element_blank(),
          legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"), 
          legend.title = element_text("Edges %", size = 9, face = "bold",  
                                      hjust = 0, vjust = 0.5), 
          legend.text = element_text(size = 9), legend.position = "right") +
    #guides(fill = guide_legend(title = "Edges")) +
    scale_fill_gradient2("Edges %", low="blue", mid="grey", high="red", #colors in the scale
                         midpoint=max(dataframe$Mean)/2,   #same midpoint for plots (mean of the range)
                         limits=c(0, max(dataframe$Mean))) # 70 for positive
  return(heat_plot)
}

# **** HEATMAP plot ------------------------------------------------------------------------------------------
PlotHeat(df_TremDot)


# Genus level HEATMAP ----------------------------------------------------------------------------------------
HeatEdgeGen <- function(df_edges, Var, Date){
  df_edges <-
    df_edges %>%
    mutate(ModuleConnectivity = ifelse(V1_Genus==V2_Genus, "in-class", "out-class"))
  df_edges <- 
    subset(df_edges, Direction%in%Var)
  df_count <-
    dplyr::count(df_edges, 
                 vars= Direction, V1_Genus, V2_Genus, ModuleConnectivity) 
  print("Subsetted dataframe")
  str(df_count) %T>% print()
  #df_count %T>% print()
  # df_count$vars <-
  #   as.factor(df_count$vars)
  # df_count_filt <-
  #   subset(df_count, vars%in%Var)
  df_count_filt <-
    df_count[, c(2,3,5)]
  print("Subsetted dataframe dimensions")
  dim(df_count_filt) %T>% print()
  # extract Class data
  df_count_filt$V1_Genus <- as.character(df_count_filt$V1_Genus)
  df_count_filt$V2_Genus <- as.character(df_count_filt$V2_Genus)
  faclist <- vector("list", nrow(df_count_filt))
  for (i in 1:nrow(df_count_filt)) {
    faclist[[i]] <-
      sort(c(df_count_filt$V1_Genus[i], df_count_filt$V2_Genus[i]))
    faclist[[i]] <- paste(faclist[[i]][1], faclist[[i]][2], sep = "_")
  }
  df_count_filt$v1v2 <- unlist(faclist)
  # print("Subsetted dataframe")
  # df_count_filt %T>% print()
  df_count_filt %>% 
    group_by(v1v2) %>% summarize(Sum=sum(n)) %>%
    as.data.frame() -> df_new
  # print("Final dataframe dimensions")
  # dim(df_new)  %T>% print()
  df_final <-
    left_join(df_new, df_count_filt[,c(1:2,4)], "v1v2")
  df_final <-
    df_final[!duplicated(df_final$v1v2), ]
  df_final$Date <- rep(Date, times=nrow(df_final))
  return(df_final)
}


df_bind_positive_gen <-
  rbind(
    HeatEdgeGen(edges_kel_all_158, "positive", "d158"),
    HeatEdgeGen(edges_kel_all_212, "positive", "d212"),
    HeatEdgeGen(edges_kel_all_233, "positive", "d233"),
    HeatEdgeGen(edges_kel_all_260, "positive", "d260"),
    HeatEdgeGen(edges_kel_all_286, "positive", "d286")) 

df_bind_negative_gen <-
  rbind(
    HeatEdgeGen(edges_kel_all_158, "negative", "d158"),
    HeatEdgeGen(edges_kel_all_212, "negative", "d212"),
    HeatEdgeGen(edges_kel_all_233, "negative", "d233"),
    HeatEdgeGen(edges_kel_all_260, "negative", "d260"),
    HeatEdgeGen(edges_kel_all_286, "negative", "d286")) 


df_bing_all_gen <-
  rbind(data.frame(Weight = rep("positive", nrow(df_bind_positive_gen)),
                   df_bind_positive_gen),
        data.frame(Weight = rep("negative", nrow(df_bind_negative_gen)), 
                   df_bind_negative_gen)) 
df_bing_all_gen



PlotHeatGen <- function(dataframe){
  heat_plot <-
    ggplot(dataframe, aes(x=V1_Genus, y=V2_Genus, fill=Sum)) +
    geom_tile(aes(height = 0.92, width = 0.92)) + 
    geom_text(aes(label = Sum), size=2.5, color="white") +
    theme_classic() +
    facet_grid(~Weight~Date, as.table = TRUE) +
    labs(title = "Within- and between-Class connections") +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), 
          axis.text.x =  element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 8), 
          axis.title = element_blank(),
          legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"), 
          legend.title = element_text("Edges %", size = 9, face = "bold",  
                                      hjust = 0, vjust = 0.5), 
          legend.text = element_text(size = 9), legend.position = "right") +
    #guides(fill = guide_legend(title = "Edges")) +
    scale_fill_gradient2("Edges %", low="blue", mid="grey", high="red", #colors in the scale
                         midpoint=max(dataframe$Sum)/2,   #same midpoint for plots (Sum of the range)
                         limits=c(0, max(dataframe$Sum))) # 70 for positive
  return(heat_plot)
}


PlotHeatGen(subset(df_bing_all_gen,
                   V2_Genus%in%c("Dioszegia")))
PlotHeatGen(subset(df_bing_all_gen,
                   V1_Genus%in%c("Dioszegia")))

PlotHeatGen(subset(df_bing_all_gen,
                   V2_Genus%in%c("Papiliotrema")))
PlotHeatGen(subset(df_bing_all_gen,
                   V1_Genus%in%c("Papiliotrema")))


# HEATMAP graph for OTUs ------------------------------------------------------------------------------------
data_heat <-rbind(
  nodes_kel_all_158[nodes_kel_all_158$Infection_ind%in%"Inf",],
  nodes_kel_all_158[!is.na(nodes_kel_all_158$Date_ind),])

data_heat <-
  data_heat[,c(1,2,5:8,32)]
data_heat$Group <-
  c(rep("Infection", times=6),rep("Date", times=10))

data_heat

DataHeat <- function(df){
  df[,c(1:3,7:8)] -> df1
  colnames(df1)[3] <- "Abund"
  df[,c(1:2,4,7:8)] -> df2
  colnames(df2)[3] <- "Abund"
  df[,c(1:2,6,7:8)] -> df3
  colnames(df3)[3] <- "Abund"
  df[,c(1:2,6,7:8)] -> df4
  colnames(df4)[3] <- "Abund"
  df_all <-
    rbind(
      rbind(df1,df2),
      rbind(df3,df4))
  return(df_all) 
}


DataHeat(data_heat) -> data_melt
rownames(data_melt) <- 1:nrow(data_melt)
data_melt$Sub <- c(
  rep("Atl", times=16), 
  rep("Gul", times=16), 
  rep("Int", times=16), 
  rep("Mid", times=16)
)

data_melt

subplot <-
  ggplot(data_melt, aes(x=Sub, y=Taxonomy, fill=Abund)) +
  geom_tile(aes(height = 0.92, width = 0.92)) + 
  theme_classic() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), 
        axis.text.x =  element_text(size = 7),
        axis.text.y = element_text(size = 7), 
        axis.title = element_blank(),
        legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
        legend.title = element_text(size = 8, face = "bold", hjust = 0, vjust = 0.5), 
        legend.text = element_text(size = 7), legend.position = "none")

inf_plot<-
  ggplot(data_melt, aes(x=Group, y=Taxonomy, fill=Group)) +
  geom_tile(aes(height = 0.92, width = 0.92))+
  theme_classic() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), 
        axis.text.x =  element_text(size = 7),
        axis.text.y = element_blank(), axis.ticks.y=element_blank(), 
        axis.title = element_blank(),
        legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
        legend.title = element_text(size = 8, face = "bold", hjust = 0, vjust = 0.5), 
        legend.text = element_text(size = 7), legend.position = "none") +
  scale_x_discrete(labels=c("D", "I"))

ggarrange(subplot,
          inf_plot,
          ncol = 2,
          nrow = 1,
          widths = c(1, 0.2))


# TESTING GUILDS --------------------------------------------------------------------------------------------
library(fungaltraits)

fungal_traits() -> db_traits
dim(db_traits)
colnames(db_traits)
db_traits_filt <- 
  db_traits[,c("Genus", "species",  "guild_fg",
               "trophic_mode_fg","trait_fg" )]

sort(rownames(nodes_kel_all_158))
nodes_kel_all_158[is.na(nodes_kel_all_158$Genus), ]

table(nodes_kel_all_158$Genus)
table(nodes_kel_all_158$Class)
dplyr::count(nodes_kel_all_158, vars= Class) -> net_otu_by_class
colnames(net_otu_by_class) <- c("Class", "OTU_number")

plot_class <-
  ggplot(net_otu_by_class, aes(x=Class, y=OTU_number)) +
  geom_bar(stat = "identity") +
  theme_classic()+
  ylim(0, 50) +
  theme(plot.title = element_text(size = 10, face = "bold", vjust = 0.5),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), 
        axis.text.x =  element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 7, angle = 0, hjust = 1, vjust = 0.5), 
        axis.title = element_text(size = 8, face = "bold"),
        legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
        legend.title = element_text(size = 8, face = "bold", hjust = 0, vjust = 0.5), 
        legend.text = element_text(size = 7), legend.position = "none")

plot_class


# Saving data for Acer ---------------------------------------------------------------------------------------
save(network_abund, 
     net_graph,
     barplot_time,
     data_melt,
     ind_fungi_date,
     ind_fungi_inf,
     network_kel_all_summary,
     file = "network.RData")

save(plot_class,
     net_otu_by_class,
     guild_match,
     file = "network_additional.RData")


# NETWORK PROPERTIES CLASS AVERAGES -------------------------------------------------------------------------
CalcNetAttribClass <- function(network, physeq) {
  require(igraph)
  require(ggrepel)
  require(microbiome)
  nodes <-
    data.frame(matrix(ncol = 3, nrow = length(V(network)$name)))
  colnames(nodes) <-
    c("Degree",
      "Modules",
      "Betweenness")
  rownames(nodes) <- V(network)$name
  nodes <- data.frame(OTU_ID = as.factor(row.names(nodes)), nodes)
  dim(nodes) %T>% print()
  components(network) %T>% print()
  #Degree
  nodes$Degree <- igraph::degree(network, mode = "all")
  # removing weights
  network_new <- delete_edge_attr(network, "weight") # delete weights
  #Moduels
  #modules <- cluster_walktrap(network_new)
  modules <- cluster_fast_greedy(network_new)
  nodes$Modules <- as.numeric(igraph::membership(modules))
  # Betweenness
  nodes$Betweenness <- igraph::betweenness(network_new, directed = FALSE)
  # extarcting taxa
  taxa <- as(tax_table(physeq), "matrix")
  taxa <- as.data.frame(taxa)
  inner_join(nodes, taxa[, c(1,4)], by = "OTU_ID") -> nodes_attr
  nodes_attr$Class <- 
    ifelse(is.na(nodes_attr$Class), "Unclassified", paste(nodes_attr$Class))
  return(nodes_attr)
}



CalcNetAttribClass(network_kel_158, physeq_kel_158) -> Class_net_attr_158
str(Class_net_attr_158)

CalcNetAttribClass(network_kel_212, physeq_kel_212) -> Class_net_attr_212
CalcNetAttribClass(network_kel_233, physeq_kel_233) -> Class_net_attr_233
CalcNetAttribClass(network_kel_260, physeq_kel_260) -> Class_net_attr_230
CalcNetAttribClass(network_kel_286, physeq_kel_286) -> Class_net_attr_286


merge(merge(merge(merge(
  Class_net_attr_158 %>% 
    group_by(Class) %>%
    summarise(Deg = mean(Degree)),
  Class_net_attr_212 %>% 
    group_by(Class) %>%
    summarise(Deg = mean(Degree)), all = TRUE, by="Class"),
  Class_net_attr_233 %>% 
    group_by(Class) %>%
    summarise(Deg = mean(Degree)), all = TRUE, by="Class"),
  Class_net_attr_230 %>% 
    group_by(Class) %>%
    summarise(Deg = mean(Degree)), all = TRUE, by="Class"),
  Class_net_attr_286 %>% 
    group_by(Class) %>%
    summarise(Deg = mean(Degree)), all = TRUE, by="Class") -> Deg_all

Deg_all
colnames(Deg_all) <- c("Class", paste("Deg", 1:5, sep="_"))
rowMeans(Deg_all[, 2:6])
apply(Deg_all[, 2:6], 1, mean)

Class_net_attr_158 %>% 
  group_by(Class) %>%
  summarise(Mod = unique(Modules))


merge(merge(merge(merge(
  Class_net_attr_158 %>% 
    group_by(Class) %>%
    summarise(Mod = mean(Modules)),
  Class_net_attr_212 %>% 
    group_by(Class) %>%
    summarise(Mod = mean(Modules)), all = TRUE, by="Class"),
  Class_net_attr_233 %>% 
    group_by(Class) %>%
    summarise(Mod = mean(Modules)), all = TRUE, by="Class"),
  Class_net_attr_230 %>% 
    group_by(Class) %>%
    summarise(Mod = mean(Modules)), all = TRUE, by="Class"),
  Class_net_attr_286 %>% 
    group_by(Class) %>%
    summarise(Mod = mean(Modules)), all = TRUE, by="Class") -> Mod_all

Mod_all
apply(Mod_all[, 2:6], 1, mean)


merge(merge(merge(merge(
  Class_net_attr_158 %>% 
    group_by(Class) %>%
    summarise(Btw = mean(Betweenness)),
  Class_net_attr_212 %>% 
    group_by(Class) %>%
    summarise(Btw = mean(Betweenness)), all = TRUE, by="Class"),
  Class_net_attr_233 %>% 
    group_by(Class) %>%
    summarise(Btw = mean(Betweenness)), all = TRUE, by="Class"),
  Class_net_attr_230 %>% 
    group_by(Class) %>%
    summarise(Btw = mean(Betweenness)), all = TRUE, by="Class"),
  Class_net_attr_286 %>% 
    group_by(Class) %>%
    summarise(Btw = mean(Betweenness)), all = TRUE, by="Class") -> Btw_all

Btw_all
apply(Btw_all[, 2:6], 1, mean)

final_stat <-
  data.frame(Class= Deg_all$Class,
             Degree = apply(Deg_all[, 2:6], 1, mean),
             Modularity = apply(Mod_all[, 2:6], 1, mean),
             Betweenness = apply(Btw_all[, 2:6], 1, mean))

final_stat

