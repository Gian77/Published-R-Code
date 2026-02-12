# ************************* DATA ANALYSIS *************************** -------------------------------------------
# Project name: fungicide on corn and soy 
# Manuscript:   
# Authors:      Noel Z, ..., Benucci GMN, ..., Bonito G.
# R Code:       Benucci GMN
# Affiliation:  Michigan State University
# Journal:      New Phytol.
# Date:         March 11, 2021
#
# Please contact <gian.benucci@gmail.com> before using this code. 
# ******************************************************************** -------------------------------------------

# required packages ----------------------------------------------------------------------------------------------
library(phyloseq); packageVersion("phyloseq")
library(Biostrings)
library(ggplot2)
library(ape)
library(dplyr)
library(ggpubr)
library(vegan); packageVersion("vegan")
library(randomForest)
library(rfUtilities)
library(scales)
library(compositions)
library(Boruta)
library(magrittr)
library(gridExtra)
library(tidyr)

# ******************************************************************************----------------------------------
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
source('functions_themes.R')
# PLOTTING FUNCTIONS ---------------------------------------------------------------------------------------------
# 1 - Plotting error ---------------------------------------------------------------------------------------------
PlotError <- function(rf_model){
  model_df <- (data.frame(Trees = 1:1001, Error = rf_model$mse))
  ggplot(data=model_df, aes(x=Trees, y=Error)) +
    labs(title = "Model Errors", y="Error", x="Tree") +
    theme_classic() +
    geom_line(color="red", size=0.8) +
    ylim(0, NA) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size =8)) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) -> error_plot
  return(error_plot)
}

# 2 - Plotting Line ----------------------------------------------------------------------------------------------
PlotLine <- function(rf_model, metadata){
  require(ggpmisc)
  require(ggpubr)
  df_model <- data.frame(actual=rf_model$y, pred=rf_model$predicted)
  if (identical(rownames(metadata), rownames(df_model))==TRUE){
    df_model$Fungicide <- metadata$Fungicide
    df_model$Date <- metadata$DateSampled
    df_model$Date <-
      factor(df_model$Date, levels = c("3-Aug-18","16-Aug-18", "27-Aug-18"))
    # df_model$Date <-
    #   recode(df_model$Date, 
    #          `3-Aug-18` = "0 dpf", 
    #          `16-Aug-18` = "13 dfp", 
    #          `27-Aug-18` = "33 dpf")
    df_model %T>% print()
    ggplot(data=df_model, aes(x=actual, y=pred)) +
      geom_point(aes(shape=Date), size=1, stroke=0.8, color="black") +
      geom_smooth(method = "lm", formula = "y ~ x", se = TRUE, color="black", size=0.8) +
      theme_classic() +
      scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
      scale_shape_manual(values = c(1, 2, 3),
                         labels = c('0 dpf', '13 dpf', '33 dpf')) +
      theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
            legend.background=element_blank()) +
      theme(legend.title = element_text(size = 8, face = "bold"), 
            legend.text = element_text(size =7)) +
      theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
      theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
      theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
      theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
      theme(legend.title = element_blank()) -> line_plot
    return(line_plot)
  }else{
    stop("Error: dataframe and metadata rownames are not matching!")
  }
}

# 3- Plotting features -------------------------------------------------------------------------------------------
PlotFeature <- function(rf_model, taxa){
  imp_RF <- as.data.frame(rf_model$importance)
  imp_RF$features <- rownames(imp_RF)
  imp_RF <- dplyr::arrange(imp_RF, desc(`%IncMSE`))
  rownames(imp_RF) <- imp_RF$features
  # adding marker taxon info
  taxa[rownames(taxa)%in%rownames(imp_RF), ] -> taxa_RF
  identical(rownames(taxa_RF), rownames(imp_RF))
  order_taxa <- match(rownames(taxa_RF), rownames(imp_RF))
  imp_RF <- imp_RF[order_taxa,]
  imp_RF$Taxonomy <- taxa_RF$Taxon
  imp_RF <- imp_RF[order(imp_RF$`%IncMSE`, decreasing = TRUE),] 
  imp_RF %T>% print()
  ggplot(data=imp_RF) + 
    geom_bar(aes(x= reorder(Taxonomy, -`%IncMSE`),
                 y= `%IncMSE`), color="grey80", fill="grey80",stat="identity") +
    labs(x= "OTU", y= "MSE Increase %") +
    coord_flip() +
    #ylim(0, 0.03) +
    theme_classic() +
    #theme(plot.margin=unit(c(7,9,7,7),"pt")) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 0, size = 7 ,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 1, vjust = 0.5)) +
    theme(axis.title.x  = element_text(angle = 0, size = 8, face = "bold", hjust = 0.5),
          axis.title.y = element_blank()) -> plot_importance
  return(plot_importance)
}

# 4- Heatmap of features abundances ------------------------------------------------------------------------------
PlotHeat <- function(dataframe){
  heat_plot <-
    ggplot(dataframe, aes(x=Date, y=OTU, fill=Sum)) +
    geom_tile() + 
    #geom_text(aes(label = ), size=2, color="white") +
    theme_classic() +
    #labs(title = "OTU Read abundance") +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"), 
          axis.text.x =  element_text(size = 7, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_blank(),  
          axis.title = element_blank(),
          axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 8, face = "bold",  
                                      hjust = 0, vjust = 0.5), 
          legend.text = element_text(size = 8), legend.position = "right") +
    scale_x_discrete(labels=c("3-Aug-18" = "0 dpf", "16-Aug-18" = "13 dpf","27-Aug-18" = "33 dpf")) +
    scale_fill_gradient2("Read No.", low="#999999", mid="white", high="#E69F00", #colors in the scale
                         midpoint=max(dataframe$Sum)/2,   #same midpoint for plots (mean of the range)
                         limits=c(0, max(dataframe$Sum))) # 70 for positive
  return(heat_plot)
}

# 5 - Extract feaures seq for BLAST ------------------------------------------------------------------------------
filterTaxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

# Importing datasets ---------------------------------------------------------------------------------------------
dir()

# CSS normalized 
physeq_fungi_CSS <- 
  readRDS("Preprocessing/Fungicide_Fungi_RootsAndLeaves_1000seqMin_CSSNorm_022421.rds")
dim(physeq_fungi_CSS@tax_table)
head(tax_table(physeq_fungi_CSS))
physeq_prok_CSS <- 
  readRDS("Preprocessing/Fungicide_Bacteria_RootsAndLeaves_newtaxa_CSSnorm_010621.rds")
dim(physeq_prok_CSS@tax_table)
head(tax_table(physeq_prok_CSS))

# Non-normalized
physeq_fungi_no_norm <- 
  readRDS("Preprocessing/Fungicide_Fungi_RootsAndLeaves_1000readsmin_nonorm_022421.rds")
dim(physeq_fungi_no_norm@tax_table)
head(tax_table(physeq_fungi_no_norm))
physeq_prok_no_norm <- 
  readRDS("Preprocessing/Fungicide_Bacteria_RootsAndLeaves_newtaxa_nonorm_010621.rds")
dim(physeq_prok_no_norm@tax_table)
head(tax_table(physeq_prok_no_norm))

# get a sense of the data
head(sample_data(physeq_fungi_CSS))
head(tax_table(physeq_fungi_CSS))

dplyr::count(as.data.frame(as.matrix(
  sample_data(physeq_fungi_CSS))),
     Objective,
      Treatment,
      Crop,
      Compartment,
      Fungicide,
      Rep)

# Removing Fungal contaminants -----------------------------------------------------------------------------------
filterTaxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}


physeq_prok_CSS <-
  filterTaxa(physeq_prok_CSS, c("BOTU_14673","BOTU_1143")) 

# >>> Hypothesis <<< ---------------------------------------------------------------------------------------------
# Is there variation in fungal diversity that can be explained by bacterial 
# diversity - Richness or Beta Diversity?

# Crop DateSampled Treatment Fungicide  n
# 1  Soy   16-Aug-18        T2         F 12
# 2  Soy   27-Aug-18        T2         F 11
# 3  Soy    3-Aug-18        T2         F 12

# >>> Random Forest models for Treatment Conventional <<< ---------------------------------------------------------------------------------------------

# Extract Fingicide-treated samples
SubsetData <- function(physeq_1, physeq_2){
  dplyr::intersect(rownames(sample_data(physeq_1)),
                   rownames(sample_data(physeq_2))) -> to_keep
  subset(otu_table(physeq_1), select = to_keep) -> otu_table(physeq_1)
  otu_table(physeq_1) <-
    otu_table(physeq_1)[which(rowSums(otu_table(physeq_1)) > 0),]
  physeq_new <-
    subset_samples(physeq_1, Crop %in% c("Soy") &
                     Compartment %in% c("Leaf") &
                     Treatment %in% c("T1") &
                     Fungicide %in% c("F"))
  sample_data(physeq_new) <-
    sample_data(physeq_new)[, c("Crop",  "DateSampled", "Treatment", "Fungicide")]
  otu_table(physeq_new) <-
    otu_table(physeq_new)[which(rowSums(otu_table(physeq_new)) > 0), ] 
  sample_data(physeq_new)$ReadNo <- sample_sums(physeq_new)
  dplyr::count(as.data.frame(as.matrix(sample_data(physeq_new))),
   Crop, DateSampled, Fungicide, Treatment) %T>% print()
  return(physeq_new)
}


SubsetData(physeq_fungi_CSS, physeq_prok_CSS) -> physeq_fungi_CSS_sub
physeq_fungi_CSS_sub
SubsetData(physeq_fungi_no_norm, physeq_prok_no_norm) -> physeq_fungi_no_norm_sub
physeq_fungi_no_norm_sub

SubsetData(physeq_prok_CSS, physeq_fungi_CSS) -> physeq_prok_CSS_sub
physeq_prok_CSS_sub
SubsetData(physeq_prok_no_norm, physeq_fungi_no_norm) -> physeq_prok_no_norm_sub
physeq_prok_no_norm_sub

# Extract Control Samples 
SubsetControl <- function(physeq_1, physeq_2){
  dplyr::intersect(rownames(sample_data(physeq_1)),
                   rownames(sample_data(physeq_2))) -> to_keep
  subset(otu_table(physeq_1), select = to_keep) -> otu_table(physeq_1)
  otu_table(physeq_1) <-
    otu_table(physeq_1)[which(rowSums(otu_table(physeq_1)) > 0),]
  physeq_new <-
    subset_samples(physeq_1, Crop %in% c("Soy") &
                     Compartment %in% c("Leaf") &
                     Treatment %in% c("T1") &
                     Fungicide %in% c("C"))
  sample_data(physeq_new) <-
    sample_data(physeq_new)[, c("Crop",  "DateSampled", "Treatment", "Fungicide")]
  otu_table(physeq_new) <-
    otu_table(physeq_new)[which(rowSums(otu_table(physeq_new)) > 0), ] 
  sample_data(physeq_new)$ReadNo <- sample_sums(physeq_new)
  dplyr::count(as.data.frame(as.matrix(sample_data(physeq_new))),
        Crop, DateSampled, Fungicide, Treatment) %T>% print()
  return(physeq_new)
}


SubsetControl(physeq_fungi_CSS, physeq_prok_CSS) -> physeq_fungi_CSS_sub_C
physeq_fungi_CSS_sub_C
SubsetControl(physeq_fungi_no_norm, physeq_prok_no_norm) -> physeq_fungi_no_norm_sub_C
physeq_fungi_no_norm_sub_C

SubsetControl(physeq_prok_CSS, physeq_fungi_CSS) -> physeq_prok_CSS_sub_C
physeq_prok_CSS_sub_C
SubsetControl(physeq_prok_no_norm, physeq_fungi_no_norm) -> physeq_prok_no_norm_sub_C
physeq_prok_no_norm_sub_C


# *** Extract Bulleribasidiaceae ---------------------------------------------------------------------------------
FilterFamily <- function(physeq){
  physeq_new <- subset_taxa(physeq, Family=="Bulleribasidiaceae")
  otu_table(physeq_new) <-
    otu_table(physeq_new)[which(rowSums(otu_table(physeq_new)) > 0), ] 
  otu = as(otu_table(physeq_new), "matrix")
  otu <- as.data.frame(t(otu))
  metadata = as(sample_data(physeq_new), "data.frame")
  metadata$mean <- apply(X = otu, MARGIN = 1, FUN = mean)
  metadata$richness <- specnumber(otu, MARGIN = 1)
  metadata$shannon <- vegan::diversity(otu,index = "shannon", MARGIN=1) 
  print(dim(metadata))
  print(dim(otu))
  list <- list(otu, metadata)
  return(list)
}

# Main samples
FilterFamily(physeq_fungi_CSS_sub) -> physeq_fungi_CSS_Bul
physeq_fungi_CSS_Bul[1]
physeq_fungi_CSS_Bul[2]
dplyr::count(physeq_fungi_CSS_Bul[[2]],Crop,DateSampled,Treatment)

FilterFamily(physeq_fungi_no_norm_sub) -> physeq_fungi_no_norm_Bul_T1
physeq_fungi_no_norm_Bul_T1[1]
physeq_fungi_no_norm_Bul_T1[2]
#dplyr::count(physeq_fungi_no_norm_sub[[2]],Crop,DateSampled,Treatment)

# Controls
FilterFamily(physeq_fungi_CSS_sub_C) -> physeq_fungi_CSS_Bul_C
physeq_fungi_CSS_Bul_C[1]
physeq_fungi_CSS_Bul_C[2]

FilterFamily(physeq_fungi_no_norm_sub_C) -> physeq_fungi_no_norm_Bul_C
physeq_fungi_no_norm_Bul_C[1]
physeq_fungi_no_norm_Bul_C[2]

# **************************************--------------------------------------------------------------------------
# RANDOM FOREST MODELS -------------------------------------------------------------------------------------------
# Extract datasets -----------------------------------------------------------------------------------------------
otu_prok_CSS_sub = as(otu_table(physeq_prok_CSS_sub), "matrix")
otu_prok_CSS_sub <- as.data.frame(t(otu_prok_CSS_sub))
taxa_prok_CSS_sub = as(tax_table(physeq_prok_CSS_sub), "matrix")
taxa_prok_CSS_sub = as.data.frame(taxa_prok_CSS_sub)

otu_prok_no_norm_sub = as(otu_table(physeq_prok_no_norm_sub), "matrix")
otu_prok_no_norm_sub <- as.data.frame(t(otu_prok_no_norm_sub))
taxa_prok_no_norm_sub = as(tax_table(physeq_prok_no_norm_sub), "matrix")
taxa_prok_no_norm_sub = as.data.frame(taxa_prok_no_norm_sub)

otu_prok_no_norm_sub_C = as(otu_table(physeq_prok_no_norm_sub_C), "matrix")
otu_prok_no_norm_sub_C <- as.data.frame(t(otu_prok_no_norm_sub_C))
taxa_prok_no_norm_sub_C = as(tax_table(physeq_prok_no_norm_sub_C), "matrix")
taxa_prok_no_norm_sub_C = as.data.frame(taxa_prok_no_norm_sub_C)

# Test different data transformation 
otu_prok_tran <- otu_prok_no_norm_sub # non-normalized
#otu_prok_tran <-
#  as.data.frame(scale(otu_prok_no_norm_sub, center = TRUE, scale = TRUE)) # each value is converted into a Z-score
#otu_prok_tran <- otu_prok_CSS_sub # Just metagenomeSeq
#otu_prok_tran <-
#  scale(asinh(otu_prok_no_norm_sub),
#        center = TRUE,
#        scale = FALSE)  #inverse hyperbolic sine and then to mean center by sample
#otu_prok_tran <-
#  as.data.frame(clr(otu_prok_no_norm_sub))  # centred log-ratio (CLR)

# The best transformation is just using reads count in terms of MSE and % Var. All 
# tranformations were giving comparable results with the Boruta function.
identical(rownames(as.data.frame(physeq_fungi_no_norm_Bul_T1[2])), rownames(otu_prok_tran))
sample_order <-
  match(rownames(as.data.frame(physeq_fungi_no_norm_Bul_T1[2])), rownames(otu_prok_tran))
otu_prok_tran <- otu_prok_tran[sample_order,]
identical(rownames(as.data.frame(physeq_fungi_no_norm_Bul_T1[2])), rownames(otu_prok_tran))
dim(otu_prok_tran)
range(otu_prok_tran)

# >>> A) Model: bacterial otu -> fungal richness Conventional Management -----------------------------------------------------------------
# adding the classification variable and sample names 
otu_prok_tran$Richness <-
  as.data.frame(physeq_fungi_no_norm_Bul_T1[2])$richness
head(otu_prok_tran)
dim(otu_prok_tran)

# Recursive feature selection
set.seed(11021)
rfe_prok_rich_T1 <- Boruta(
  Richness ~ .,
  otu_prok_tran,
  pValue = 0.05,
  mcAdj = TRUE,
  maxRuns = ncol(otu_prok_tran),
  doTrace = 3
)

rfe_prok_rich_T1

# Get significant variables including tentatives
rfe_prok_rich_T1 <-
  getSelectedAttributes(rfe_prok_rich_T1, withTentative = TRUE)
rfe_prok_rich_T1

otu_prok_tran[, c(rfe_prok_rich_T1)] -> otu_prok_tran_sel_T1
identical(rownames(otu_prok_tran_sel_T1), rownames(otu_prok_tran))
otu_prok_tran_sel_T1$Richness <- otu_prok_tran$Richness

# try tuning the model first
round(sqrt(ncol(otu_prok_tran_sel_T1[, 1:(ncol(otu_prok_tran_sel_T1) - 1)])))

set.seed(11022)
bestmtry_rich_T1 <-
  tuneRF(
    x = otu_prok_tran_sel_T1[, 1:(ncol(otu_prok_tran_sel_T1) - 1)],
    y = otu_prok_tran_sel_T1$Richness,
    mtryStart = 3,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )

bestmtry_rich_T1

set.seed(11033)
RF_prok_rich_T1 <-
  randomForest(
    x = otu_prok_tran_sel_T1[, 1:(ncol(otu_prok_tran_sel_T1) - 1)],
    y = otu_prok_tran_sel_T1$Richness,
    ntree = 1001,
    mtry = 1,
    importance = TRUE,
    proximity = TRUE
  )

RF_prok_rich_T1
plot(RF_prok_rich_T1)

# Assessing model significance using permitations
set.seed(110324)
perm_RF_prok_rich_T1 <-
  rf.significance(
    x = RF_prok_rich_T1,
    xdata = otu_prok_tran_sel_T1[, 1:(ncol(otu_prok_tran_sel_T1) - 1)],
    nperm = 999,
    nmtry = 1,
    ntree = 1001
  )

perm_RF_prok_rich_T1 # model significant = 0.002

# *** FIGURE RF model Richness Conventional Management ---------------------------------------------------------------------------------

#df for heatmap
otu_prok_tran_sel_T1$Date <-
  as.data.frame(physeq_fungi_no_norm_Bul_T1[2])$DateSampled

df_melt_T1 <-
  subset(otu_prok_tran_sel_T1, select = -c(Richness)) %>%
  reshape2::melt(value.name = "Count", variable.name = "OTU") %>%
  group_by(Date, OTU) %>%
  dplyr::summarise(Sum = sum(Count)) %>%
  as.data.frame()

rel.abund.bacteria.selected.T1 <-
  physeq_prok_no_norm %>% 
  subset_samples(Compartment == "Leaf" & Crop == "Soy" & Treatment == "T1") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt.fast() %>%
  subset(OTU %in% rfe_prok_rich_T1) %>%
  mutate(Date = as.Date(DateSampled, "%d-%b-%y")) %>%
  ggplot(aes(x = Date, y = Abundance, color = Fungicide, linetype = Fungicide)) +
  stat_summary(fun.y=mean,geom="line") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  theme_classic()+
  facet_wrap(~Taxonomy, scales = "free_y") + 
  scale_y_continuous(labels = scales::percent)

df_melt_T1$Date <-
  factor(df_melt_T1$Date, levels = c("3-Aug-18", "16-Aug-18", "27-Aug-18"))

# for T1
df_melt_T1$OTU <-
  factor(df_melt_T1$OTU, levels = c(rfe_prok_rich_T1))

# Correcting Taxonomy 
str(taxa_prok_no_norm_sub)
taxa_prok_no_norm_sub$Taxonomy <-
  as.character(taxa_prok_no_norm_sub$Taxonomy)

taxa_prok_no_norm_sub[taxa_prok_no_norm_sub == "BOTU_24-Methylobacterium-methylorubrum sp."] <-
  "BOTU_24-Methylobacterium sp."
taxa_prok_no_norm_sub[taxa_prok_no_norm_sub$Taxonomy == "BOTU_24-Methylobacterium sp.",]

ggpubr::ggarrange(PlotError(RF_prok_rich_T1) +
            annotate("text", x=Inf, y = Inf, 
                     label= paste("% Explained var.:", round(last(RF_prok_rich_T1$rsq*100),2)), size=2.5, vjust=1, hjust=1) +
            annotate("text", x=Inf, y = Inf,
                     label=paste("MSE:", round(last(RF_prok_rich_T1$mse), 2)), size=2.5, vjust=3, hjust=1) +
            annotate("text", x=Inf, y = Inf, 
                     label = paste("italic(p) ==", round(perm_RF_prok_rich_T1$pValue, 4)), parse = TRUE, size=2.5, vjust=4, hjust=1),
          PlotLine(RF_prok_rich_T1,  as.data.frame(physeq_fungi_no_norm_Bul_T1[2])) + 
            theme(legend.position = c(0.2, 0.95)) +
            labs(title = "Bulleribacidaceae\nrichness ", y="Predicted", x="Observed"),
          PlotFeature(RF_prok_rich_T1, taxa_prok_no_norm_sub) +  
            labs(title = "Bacterial\nOTU importance"),
          PlotHeat(df_melt_T1),
          labels = c("A","B","C",""),
          widths = c(0.6, 0.6, 1, 0.4),
          align = "h",
          ncol = 4, 
          nrow =1) -> Fig_1_RF_rich_T1

Fig_1_RF_rich_T1

# >>> B) Control Model: bacterial otu -> fungal richness Conventional Management ------------------------------------------------------------
# generate OTU dataframe
otu_prok_tran_CR <- otu_prok_no_norm_sub_C 
identical(rownames(as.data.frame(physeq_fungi_no_norm_Bul_C[2])), rownames(otu_prok_tran_CR))
sample_order_CR <- match(rownames(as.data.frame(physeq_fungi_no_norm_Bul_C[2])), rownames(otu_prok_tran_CR))
otu_prok_tran_CR <- otu_prok_tran_CR[sample_order_CR,]
identical(rownames(as.data.frame(physeq_fungi_no_norm_Bul_C[2])), rownames(otu_prok_tran_CR))
dim(otu_prok_tran_CR)

# adding the classification variable
otu_prok_tran_CR$Richness <- as.data.frame(physeq_fungi_no_norm_Bul_C[2])$richness
head(otu_prok_tran_CR)
dim(otu_prok_tran_CR)

# select OTUs
otu_prok_tran_CR %>% 
  select(rfe_prok_rich) -> otu_prok_tran_sel_CR
identical(rownames(otu_prok_tran_sel_CR), rownames(otu_prok_tran_CR))
otu_prok_tran_sel_CR$Richness <- otu_prok_tran_CR$Richness

# try tuning the model first
round(sqrt(ncol(otu_prok_tran_sel_CR[,1:(ncol(otu_prok_tran_sel_CR)-1)])))

set.seed(110325)
bestmtry_rich_CR <- tuneRF(x = otu_prok_tran_sel_CR[,1:(ncol(otu_prok_tran_sel_CR)-1)], 
                           y = otu_prok_tran_sel_CR$Richness,
                           mtryStart = 3, 
                           ntreeTry= 1001, 
                           improve=0.01,
                           stepFactor=0.5,
                           nodesize=1,
                           doBest = TRUE,
                           trace=TRUE,
                           plot=TRUE)

bestmtry_rich_CR

set.seed(110326)
RF_prok_rich_CR <- randomForest(x=otu_prok_tran_sel_CR[,1:(ncol(otu_prok_tran_sel_CR)-1)],
                                y=otu_prok_tran_sel_CR$Richness,
                                ntree=1001,
                                mtry = 6,
                                importance=TRUE,
                                proximity=TRUE)

RF_prok_rich_CR
plot(RF_prok_rich_CR)

# Significance 
set.seed(110327)
perm_RF_prok_rich_CR <- rf.significance(x=RF_prok_rich_CR, xdata=otu_prok_tran_sel_CR[,1:(ncol(otu_prok_tran_sel_CR)-1)],
                                        nperm=999,
                                        nmtry=1, 
                                        ntree=1001)  
perm_RF_prok_rich_CR # model not significant

# >>> C) Model: bacterial otu -> fungal NMDS1 Conventional Management --------------------------------------------------------------------

otu_fungi_CSS_Bul <- as(physeq_fungi_CSS_Bul[[1]], "matrix")
otu_fungi_CSS_Bul <- as.data.frame(t(otu_fungi_CSS_Bul))

nmds_Bul <- metaMDS(t(otu_fungi_CSS_Bul), k=2, trymax=200, autotransform = FALSE)
df_nmds_Bul <-  as.data.frame(nmds_Bul$points)

otu_prok_tran <- otu_prok_no_norm_sub 
identical(rownames(df_nmds_Bul), rownames(otu_prok_tran))
sample_order <- match(rownames(df_nmds_Bul), rownames(otu_prok_tran))
otu_prok_tran <- otu_prok_tran[sample_order,]
identical(rownames(df_nmds_Bul), rownames(otu_prok_tran))
dim(otu_prok_tran)

# adding the classification variable and sample names 
otu_prok_tran$MDS1 <- df_nmds_Bul$MDS1
head(otu_prok_tran)
dim(otu_prok_tran)

# Recursive feature selectione
set.seed(110328)
rfe_prok_nmds1 <- Boruta(
  MDS1 ~ .,
  otu_prok_tran,
  pValue = 0.05,
  mcAdj = TRUE,
  maxRuns = ncol(otu_prok_tran),
  doTrace = 3
)
rfe_prok_nmds1

# Get significant variables including tentatives
rfe_prok_nmds1 <- getSelectedAttributes(rfe_prok_nmds1, withTentative = TRUE)
rfe_prok_nmds1

otu_prok_tran[, c(rfe_prok_nmds1)] -> otu_prok_tran_sel_nmds1
identical(rownames(otu_prok_tran_sel_nmds1), rownames(otu_prok_tran))
otu_prok_tran_sel_nmds1$MDS1 <- otu_prok_tran$MDS1

# try tuning the model first
round(sqrt(ncol(otu_prok_tran_sel_nmds1[,1:(ncol(otu_prok_tran_sel_nmds1)-1)])))

set.seed(110329)
bestmtry_nmds1 <-
  tuneRF(
    x = otu_prok_tran_sel_nmds1[, 1:(ncol(otu_prok_tran_sel_nmds1) - 1)],
    y = otu_prok_tran_sel_nmds1$MDS1,
    mtryStart = 3,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )

bestmtry_nmds1

set.seed(110330)
RF_prok_nmds1 <-
  randomForest(
    x = otu_prok_tran_sel_nmds1[, 1:(ncol(otu_prok_tran_sel_nmds1) - 1)],
    y = otu_prok_tran_sel_nmds1$MDS1,
    ntree = 1001,
    mtry = 3,
    importance = TRUE,
    proximity = TRUE
  )

RF_prok_nmds1
plot(RF_prok_nmds1)

# Assessing model significance using permitations 
set.seed(110331)
perm_RF_prok_nmds1 <-
  rf.significance(
    x = RF_prok_nmds1,
    xdata = otu_prok_tran_sel_nmds1[, 1:(ncol(otu_prok_tran_sel_nmds1) - 2)],
    nperm = 999,
    nmtry = 1,
    ntree = 1001
  )
perm_RF_prok_nmds1 # model significant = 0.001

# melting the data 
identical(rownames(otu_prok_tran_sel_nmds1), 
          rownames(as.data.frame(physeq_fungi_no_norm_Bul[2])))

otu_prok_tran_sel_nmds1$Date <-
  as.data.frame(physeq_fungi_no_norm_Bul[2])$DateSampled

df_melt_nmds1 <-
  subset(otu_prok_tran_sel_nmds1, select=-c(MDS1)) %>%
  reshape2::melt(value.name = "Count", variable.name = "OTU" ) %>%
  group_by(Date, OTU) %>% 
  dplyr::summarise(Sum = sum(Count)) %>%
  as.data.frame()

df_melt_nmds1$Date <-
  factor(df_melt_nmds1$Date, levels = c("3-Aug-18","16-Aug-18", "27-Aug-18"))

df_melt_nmds1$OTU <-
  factor(df_melt_nmds1$OTU, levels = c("BOTU_483", "BOTU_5895","BOTU_557", 
                                       "BOTU_80", "BOTU_109", "BOTU_5403",
                                       "BOTU_528","BOTU_44","BOTU_129"))

# Correcting Taxonomy 
taxa_prok_no_norm_sub[taxa_prok_no_norm_sub == "BOTU_129-Birii41"] <- "BOTU_129-Myxococcales"
taxa_prok_no_norm_sub[taxa_prok_no_norm_sub$Taxonomy=="BOTU_129-Myxococcales", ]

taxa_prok_no_norm_sub[taxa_prok_no_norm_sub == "BOTU_528-67-14"] <- "BOTU_528-Actinobacteria"
taxa_prok_no_norm_sub[taxa_prok_no_norm_sub$Taxonomy=="BOTU_528-Actinobacteria", ]

taxa_prok_no_norm_sub[taxa_prok_no_norm_sub == "BOTU_483-Microcoleus pcc-7113 sp."] <- "BOTU_483-Microcoleus sp."
taxa_prok_no_norm_sub[taxa_prok_no_norm_sub$Taxonomy=="BOTU_483-Microcoleus sp.", ]

taxa_prok_no_norm_sub[taxa_prok_no_norm_sub == "BOTU_17-Methylobacterium-methylorubrum sp."] <-
  "BOTU_17-Methylobacterium sp."
taxa_prok_no_norm_sub[taxa_prok_no_norm_sub$Taxonomy == "BOTU_17-Methylobacterium sp.",]


# *** FIGURE RF model Beta DiversityConventional Management ---------------------------------------------------------------------------
ggarrange(PlotError(RF_prok_nmds1) +
            annotate("text", x=Inf, y = Inf, 
                     label= paste("% Explained var.:", round(last(RF_prok_nmds1$rsq*100),2)), size=2.5, vjust=1, hjust=1) +
            annotate("text", x=Inf, y = Inf,
                     label=paste("MSE:", round(last(RF_prok_nmds1$mse), 2)), size=2.5, vjust=3, hjust=1) +
            annotate("text", x=Inf, y = Inf, 
                     label = paste("italic(p) ==", round(perm_RF_prok_nmds1$pValue, 4)), parse = TRUE, size=2.5, vjust=4, hjust=1),
          PlotLine(RF_prok_nmds1,  as.data.frame(physeq_fungi_no_norm_Bul[2])) + 
            theme(legend.position = c(0.2, 0.95)) +
            labs(title = "Bulleribacidaceae\nbeta diversity", y="Predicted", x="Observed"),
          PlotFeature(RF_prok_nmds1, taxa_prok_no_norm_sub) +
            labs(title = "Bacterial\nOTU importance"),
          PlotHeat(df_melt_nmds1),
          labels = c("A","B","C",""),
          widths = c(0.6, 0.6, 1, 0.4),
          align = "h",
          ncol = 4, 
          nrow =1) -> Fig_1_RF_nmds1

Fig_1_RF_nmds1


# >>> D) Model Control: bacterial otu -> fungal NMDS1 Conventional Management ------------------------------------------------------------

# generate OTU dataframe
otu_prok_tran_CD <- otu_prok_no_norm_sub_C
identical(rownames(as.data.frame(physeq_fungi_no_norm_Bul_C[2])), rownames(otu_prok_tran_CD))
sample_order_CD <-
  match(rownames(as.data.frame(physeq_fungi_no_norm_Bul_C[2])), rownames(otu_prok_tran_CD))
otu_prok_tran_CD <- otu_prok_tran_CD[sample_order_CD, ]
identical(rownames(as.data.frame(physeq_fungi_no_norm_Bul_C[2])), rownames(otu_prok_tran_CD))
dim(otu_prok_tran_CD)

# Calculate NMDS
otu_fungi_CSS_Bul_CD <- as(physeq_fungi_CSS_Bul_C[[1]], "matrix")
otu_fungi_CSS_Bul_CD <- as.data.frame(t(otu_fungi_CSS_Bul_CD))

nmds_Bul_CD <-
  metaMDS(
    t(otu_fungi_CSS_Bul_CD),
    k = 2,
    trymax = 200,
    autotransform = FALSE)

stressplot(nmds_Bul_CD)
df_nmds_Bul_CD <-  as.data.frame(nmds_Bul_CD$points)

# adding the classification variable
otu_prok_tran_CD$MDS1 <- df_nmds_Bul_CD$MDS1
head(otu_prok_tran_CD)
dim(otu_prok_tran_CD)

# select OTUs - BOTU_528 doesn't exist!
otu_prok_tran_CD %>%
  select(rfe_prok_nmds1[-8]) -> otu_prok_tran_sel_CD
identical(rownames(otu_prok_tran_sel_CD), rownames(otu_prok_tran_CD))
otu_prok_tran_sel_CD$MDS1 <- otu_prok_tran_CD$MDS1

# try tuning the model first
round(sqrt(ncol(otu_prok_tran_sel_CD[,1:(ncol(otu_prok_tran_sel_CD)-1)])))

set.seed(110332)
bestmtry_nmds1_CD <-
  tuneRF(
    x = otu_prok_tran_sel_CD[, 1:(ncol(otu_prok_tran_sel_CD) - 1)],
    y = otu_prok_tran_sel_CD$MDS1,
    mtryStart = 3,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )

bestmtry_nmds1_CD

set.seed(110333)
RF_prok_nmds1_CD <-
  randomForest(
    x = otu_prok_tran_sel_CD[, 1:(ncol(otu_prok_tran_sel_CD) - 1)],
    y = otu_prok_tran_sel_CD$MDS1,
    ntree = 1001,
    mtry = 1,
    importance = TRUE,
    proximity = TRUE
  )

RF_prok_nmds1_CD
plot(RF_prok_nmds1_CD)

# Assessing model significance using permitations  
set.seed(110334)
perm_RF_prok_nmds1_CD <-
  rf.significance(
    x = RF_prok_nmds1_CD,
    xdata = otu_prok_tran_sel_CD[, 1:(ncol(otu_prok_tran_sel_CD) - 1)],
    nperm = 999,
    nmtry = 1,
    ntree = 1001
  )
perm_RF_prok_nmds1_CD # model not significant = 0.79







# >>> Random Forest models for No-till Treatment <<< ---------------------------------------------------------------------------------------------

# Extract Fingicide-treated samples
SubsetData <- function(physeq_1, physeq_2){
  dplyr::intersect(rownames(sample_data(physeq_1)),
                   rownames(sample_data(physeq_2))) -> to_keep
  subset(otu_table(physeq_1), select = to_keep) -> otu_table(physeq_1)
  otu_table(physeq_1) <-
    otu_table(physeq_1)[which(rowSums(otu_table(physeq_1)) > 0),]
  physeq_new <-
    subset_samples(physeq_1, Crop %in% c("Soy") &
                     Compartment %in% c("Leaf") &
                     Treatment %in% c("T2") &
                     Fungicide %in% c("F"))
  sample_data(physeq_new) <-
    sample_data(physeq_new)[, c("Crop",  "DateSampled", "Treatment", "Fungicide")]
  otu_table(physeq_new) <-
    otu_table(physeq_new)[which(rowSums(otu_table(physeq_new)) > 0), ] 
  sample_data(physeq_new)$ReadNo <- sample_sums(physeq_new)
  dplyr::count(as.data.frame(as.matrix(sample_data(physeq_new))),
               Crop, DateSampled, Fungicide, Treatment) %T>% print()
  return(physeq_new)
}


SubsetData(physeq_fungi_CSS, physeq_prok_CSS) -> physeq_fungi_CSS_sub
physeq_fungi_CSS_sub
SubsetData(physeq_fungi_no_norm, physeq_prok_no_norm) -> physeq_fungi_no_norm_sub
physeq_fungi_no_norm_sub

SubsetData(physeq_prok_CSS, physeq_fungi_CSS) -> physeq_prok_CSS_sub
physeq_prok_CSS_sub
SubsetData(physeq_prok_no_norm, physeq_fungi_no_norm) -> physeq_prok_no_norm_sub
physeq_prok_no_norm_sub

# Extract Control Samples 
SubsetControl <- function(physeq_1, physeq_2){
  dplyr::intersect(rownames(sample_data(physeq_1)),
                   rownames(sample_data(physeq_2))) -> to_keep
  subset(otu_table(physeq_1), select = to_keep) -> otu_table(physeq_1)
  otu_table(physeq_1) <-
    otu_table(physeq_1)[which(rowSums(otu_table(physeq_1)) > 0),]
  physeq_new <-
    subset_samples(physeq_1, Crop %in% c("Soy") &
                     Compartment %in% c("Leaf") &
                     Treatment %in% c("T2") &
                     Fungicide %in% c("C"))
  sample_data(physeq_new) <-
    sample_data(physeq_new)[, c("Crop",  "DateSampled", "Treatment", "Fungicide")]
  otu_table(physeq_new) <-
    otu_table(physeq_new)[which(rowSums(otu_table(physeq_new)) > 0), ] 
  sample_data(physeq_new)$ReadNo <- sample_sums(physeq_new)
  dplyr::count(as.data.frame(as.matrix(sample_data(physeq_new))),
               Crop, DateSampled, Fungicide, Treatment) %T>% print()
  return(physeq_new)
}


SubsetControl(physeq_fungi_CSS, physeq_prok_CSS) -> physeq_fungi_CSS_sub_C
physeq_fungi_CSS_sub_C
SubsetControl(physeq_fungi_no_norm, physeq_prok_no_norm) -> physeq_fungi_no_norm_sub_C
physeq_fungi_no_norm_sub_C

SubsetControl(physeq_prok_CSS, physeq_fungi_CSS) -> physeq_prok_CSS_sub_C
physeq_prok_CSS_sub_C
SubsetControl(physeq_prok_no_norm, physeq_fungi_no_norm) -> physeq_prok_no_norm_sub_C
physeq_prok_no_norm_sub_C


# *** Extract Bulleribasidiaceae ---------------------------------------------------------------------------------
FilterFamily <- function(physeq){
  physeq_new <- subset_taxa(physeq, Family=="Bulleribasidiaceae")
  otu_table(physeq_new) <-
    otu_table(physeq_new)[which(rowSums(otu_table(physeq_new)) > 0), ] 
  otu = as(otu_table(physeq_new), "matrix")
  otu <- as.data.frame(t(otu))
  metadata = as(sample_data(physeq_new), "data.frame")
  metadata$mean <- apply(X = otu, MARGIN = 1, FUN = mean)
  metadata$richness <- specnumber(otu, MARGIN = 1)
  metadata$shannon <- vegan::diversity(otu,index = "shannon", MARGIN=1) 
  print(dim(metadata))
  print(dim(otu))
  list <- list(otu, metadata)
  return(list)
}

# Main samples
FilterFamily(physeq_fungi_CSS_sub) -> physeq_fungi_CSS_Bul
physeq_fungi_CSS_Bul[1]
physeq_fungi_CSS_Bul[2]
dplyr::count(physeq_fungi_CSS_Bul[[2]],Crop,DateSampled,Treatment)

FilterFamily(physeq_fungi_no_norm_sub) -> physeq_fungi_no_norm_Bul_T2
physeq_fungi_no_norm_Bul_T2[1]
physeq_fungi_no_norm_Bul_T2[2]
#dplyr::count(physeq_fungi_no_norm_sub[[2]],Crop,DateSampled,Treatment)

# Controls
FilterFamily(physeq_fungi_CSS_sub_C) -> physeq_fungi_CSS_Bul_C
physeq_fungi_CSS_Bul_C[1]
physeq_fungi_CSS_Bul_C[2]

FilterFamily(physeq_fungi_no_norm_sub_C) -> physeq_fungi_no_norm_Bul_C
physeq_fungi_no_norm_Bul_C[1]
physeq_fungi_no_norm_Bul_C[2]

# **************************************--------------------------------------------------------------------------
# RANDOM FOREST MODELS -------------------------------------------------------------------------------------------
# Extract datasets -----------------------------------------------------------------------------------------------
otu_prok_CSS_sub = as(otu_table(physeq_prok_CSS_sub), "matrix")
otu_prok_CSS_sub <- as.data.frame(t(otu_prok_CSS_sub))
taxa_prok_CSS_sub = as(tax_table(physeq_prok_CSS_sub), "matrix")
taxa_prok_CSS_sub = as.data.frame(taxa_prok_CSS_sub)

otu_prok_no_norm_sub = as(otu_table(physeq_prok_no_norm_sub), "matrix")
otu_prok_no_norm_sub <- as.data.frame(t(otu_prok_no_norm_sub))
taxa_prok_no_norm_sub = as(tax_table(physeq_prok_no_norm_sub), "matrix")
taxa_prok_no_norm_sub = as.data.frame(taxa_prok_no_norm_sub)

otu_prok_no_norm_sub_C = as(otu_table(physeq_prok_no_norm_sub_C), "matrix")
otu_prok_no_norm_sub_C <- as.data.frame(t(otu_prok_no_norm_sub_C))
taxa_prok_no_norm_sub_C = as(tax_table(physeq_prok_no_norm_sub_C), "matrix")
taxa_prok_no_norm_sub_C = as.data.frame(taxa_prok_no_norm_sub_C)

# Test different data transformation 
otu_prok_tran <- otu_prok_no_norm_sub # non-normalized
#otu_prok_tran <-
#  as.data.frame(scale(otu_prok_no_norm_sub, center = TRUE, scale = TRUE)) # each value is converted into a Z-score
#otu_prok_tran <- otu_prok_CSS_sub # Just metagenomeSeq
#otu_prok_tran <-
#  scale(asinh(otu_prok_no_norm_sub),
#        center = TRUE,
#        scale = FALSE)  #inverse hyperbolic sine and then to mean center by sample
#otu_prok_tran <-
#  as.data.frame(clr(otu_prok_no_norm_sub))  # centred log-ratio (CLR)

# The best transformation is just using reads count in terms of MSE and % Var. All 
# tranformations were giving comparable results with the Boruta function.
identical(rownames(as.data.frame(physeq_fungi_no_norm_Bul_T2[2])), rownames(otu_prok_tran))
sample_order <-
  match(rownames(as.data.frame(physeq_fungi_no_norm_Bul_T2[2])), rownames(otu_prok_tran))
otu_prok_tran <- otu_prok_tran[sample_order,]
identical(rownames(as.data.frame(physeq_fungi_no_norm_Bul_T2[2])), rownames(otu_prok_tran))
dim(otu_prok_tran)
range(otu_prok_tran)

# >>> E) Model: bacterial otu -> fungal richness No-till Management -----------------------------------------------------------------
# adding the classification variable and sample names 
otu_prok_tran$Richness <-
  as.data.frame(physeq_fungi_no_norm_Bul_T2[2])$richness
head(otu_prok_tran)
dim(otu_prok_tran)

# Recursive feature selection
set.seed(11021)
rfe_prok_rich_T2 <- Boruta(
  Richness ~ .,
  otu_prok_tran,
  pValue = 0.05,
  mcAdj = TRUE,
  maxRuns = ncol(otu_prok_tran),
  doTrace = 3
)
rfe_prok_rich_T2

# Get significant variables including tentatives
rfe_prok_rich_T2 <-
  getSelectedAttributes(rfe_prok_rich_T2, withTentative = TRUE)
rfe_prok_rich_T2

otu_prok_tran[, c(rfe_prok_rich_T2)] -> otu_prok_tran_sel
identical(rownames(otu_prok_tran_sel), rownames(otu_prok_tran))
otu_prok_tran_sel$Richness <- otu_prok_tran$Richness

# try tuning the model first
round(sqrt(ncol(otu_prok_tran_sel[, 1:(ncol(otu_prok_tran_sel) - 1)])))

set.seed(11022)
bestmtry_rich_T2 <-
  tuneRF(
    x = otu_prok_tran_sel[, 1:(ncol(otu_prok_tran_sel) - 1)],
    y = otu_prok_tran_sel$Richness,
    mtryStart = 3,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )

bestmtry_rich_T2

set.seed(11033)
RF_prok_rich_T2 <-
  randomForest(
    x = otu_prok_tran_sel[, 1:(ncol(otu_prok_tran_sel) - 1)],
    y = otu_prok_tran_sel$Richness,
    ntree = 1001,
    mtry = 1,
    importance = TRUE,
    proximity = TRUE
  )

RF_prok_rich_T2
plot(RF_prok_rich_T2)

# Assessing model significance using permitations
set.seed(110324)
perm_RF_prok_rich_T2 <-
  rf.significance(
    x = RF_prok_rich_T2,
    xdata = otu_prok_tran_sel[, 1:(ncol(otu_prok_tran_sel) - 1)],
    nperm = 999,
    nmtry = 1,
    ntree = 1001
  )

perm_RF_prok_rich_T2 # model significant = 0.002

# *** FIGURE RF model Richness No-till Management ---------------------------------------------------------------------------------

#df for heatmap
otu_prok_tran_sel$Date <-
  as.data.frame(physeq_fungi_no_norm_Bul_T2[2])$DateSampled

df_melt_T2 <-
  subset(otu_prok_tran_sel, select = -c(Richness)) %>%
  reshape2::melt(value.name = "Count", variable.name = "OTU") %>%
  group_by(Date, OTU) %>%
  dplyr::summarise(Sum = sum(Count)) %>%
  as.data.frame()

rel.abund.bacteria.selected.T2 <-
  physeq_prok_no_norm %>% 
  subset_samples(Compartment == "Leaf" & Crop == "Soy" & Treatment == "T2") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) %>%
  psmelt.fast() %>%
  subset(OTU %in% rfe_prok_rich_T2) %>%
  mutate(Date = as.Date(DateSampled, "%d-%b-%y")) %>%
  ggplot(aes(x = Date, y = Abundance, color = Fungicide, linetype = Fungicide)) +
  stat_summary(fun.y=mean,geom="line") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
  stat_compare_means(method = "t.test", label = "p.signif") +
  theme_classic() +
  facet_wrap(~Taxonomy, scales = "free_y") + 
  scale_y_continuous(labels = scales::percent)

df_melt_T2$Date <-
  factor(df_melt_T2$Date, levels = c("3-Aug-18", "16-Aug-18", "27-Aug-18"))

# for T2
df_melt_T2$OTU <-
  factor(df_melt_T2$OTU, levels = c(rfe_prok_rich_T2))

# Correcting Taxonomy 
str(taxa_prok_no_norm_sub)
taxa_prok_no_norm_sub$Taxonomy <-
  as.character(taxa_prok_no_norm_sub$Taxonomy)

taxa_prok_no_norm_sub[taxa_prok_no_norm_sub == "BOTU_24-Methylobacterium-methylorubrum sp."] <-
  "BOTU_24-Methylobacterium sp."
taxa_prok_no_norm_sub[taxa_prok_no_norm_sub$Taxonomy == "BOTU_24-Methylobacterium sp.",]

ggpubr::ggarrange(PlotError(RF_prok_rich_T2) +
                    annotate("text", x=Inf, y = Inf, 
                             label= paste("% Explained var.:", round(last(RF_prok_rich_T2$rsq*100),2)), size=2.5, vjust=1, hjust=1) +
                    annotate("text", x=Inf, y = Inf,
                             label=paste("MSE:", round(last(RF_prok_rich_T2$mse), 2)), size=2.5, vjust=3, hjust=1) +
                    annotate("text", x=Inf, y = Inf, 
                             label = paste("italic(p) ==", round(perm_RF_prok_rich_T2$pValue, 4)), parse = TRUE, size=2.5, vjust=4, hjust=1),
                  PlotLine(RF_prok_rich_T2,  as.data.frame(physeq_fungi_no_norm_Bul_T2[2])) + 
                    theme(legend.position = c(0.2, 0.95)) +
                    labs(title = "Bulleribacidaceae\nrichness ", y="Predicted", x="Observed"),
                  PlotFeature(RF_prok_rich_T2, taxa_prok_no_norm_sub) +  
                    labs(title = "Bacterial\nOTU importance"),
                  PlotHeat(df_melt_T2),
                  labels = c("A","B","C",""),
                  widths = c(0.6, 0.6, 1, 0.4),
                  align = "h",
                  ncol = 4, 
                  nrow =1) -> Fig_1_RF_rich_T2

Fig_1_RF_rich_T2

# >>> F) Control Model: bacterial otu -> fungal richness No-till Management ------------------------------------------------------------
# generate OTU dataframe
otu_prok_tran_CR <- otu_prok_no_norm_sub_C 
identical(rownames(as.data.frame(physeq_fungi_no_norm_Bul_C[2])), rownames(otu_prok_tran_CR))
sample_order_CR <- match(rownames(as.data.frame(physeq_fungi_no_norm_Bul_C[2])), rownames(otu_prok_tran_CR))
otu_prok_tran_CR <- otu_prok_tran_CR[sample_order_CR,]
identical(rownames(as.data.frame(physeq_fungi_no_norm_Bul_C[2])), rownames(otu_prok_tran_CR))
dim(otu_prok_tran_CR)

# adding the classification variable
otu_prok_tran_CR$Richness <- as.data.frame(physeq_fungi_no_norm_Bul_C[2])$richness
head(otu_prok_tran_CR)
dim(otu_prok_tran_CR)

# select OTUs
otu_prok_tran_CR %>% 
  select(rfe_prok_rich) -> otu_prok_tran_sel_CR
identical(rownames(otu_prok_tran_sel_CR), rownames(otu_prok_tran_CR))
otu_prok_tran_sel_CR$Richness <- otu_prok_tran_CR$Richness

# try tuning the model first
round(sqrt(ncol(otu_prok_tran_sel_CR[,1:(ncol(otu_prok_tran_sel_CR)-1)])))

set.seed(110325)
bestmtry_rich_CR <- tuneRF(x = otu_prok_tran_sel_CR[,1:(ncol(otu_prok_tran_sel_CR)-1)], 
                           y = otu_prok_tran_sel_CR$Richness,
                           mtryStart = 3, 
                           ntreeTry= 1001, 
                           improve=0.01,
                           stepFactor=0.5,
                           nodesize=1,
                           doBest = TRUE,
                           trace=TRUE,
                           plot=TRUE)

bestmtry_rich_CR

set.seed(110326)
RF_prok_rich_CR <- randomForest(x=otu_prok_tran_sel_CR[,1:(ncol(otu_prok_tran_sel_CR)-1)],
                                y=otu_prok_tran_sel_CR$Richness,
                                ntree=1001,
                                mtry = 6,
                                importance=TRUE,
                                proximity=TRUE)

RF_prok_rich_CR
plot(RF_prok_rich_CR)

# Significance 
set.seed(110327)
perm_RF_prok_rich_CR <- rf.significance(x=RF_prok_rich_CR, xdata=otu_prok_tran_sel_CR[,1:(ncol(otu_prok_tran_sel_CR)-1)],
                                        nperm=999,
                                        nmtry=1, 
                                        ntree=1001)  
perm_RF_prok_rich_CR # model not significant



# ******************************************************************************----------------------------------
# FINAL FIGURES  -------------------------------------------------------------------------------------------------


# *** FIGURE 1 ---------------------------------------------------------------------------------------------------
Fig_X_RF_Bul <-
  ggarrange(PlotLine(RF_prok_rich_T1,  as.data.frame(physeq_fungi_no_norm_Bul_T1[2])) + 
              theme(legend.position = c(0.2, 0.95)) +
              labs(title = "Bulleribasidiaceae\nrichness conv.", y="Predicted", x="Observed"),
            PlotFeature(RF_prok_rich_T1, taxa_prok_no_norm_sub) +  
              labs(title = "Bacterial\nOTU importance"),
            PlotHeat(df_melt_T1),
            PlotLine(RF_prok_rich_T2,  as.data.frame(physeq_fungi_no_norm_Bul_T2[2])) + 
              theme(legend.position = c(0.2, 0.95)) +
              labs(title = "Bulleribasidiaceae\nrichness no-till", y="Predicted", x="Observed"),
            PlotFeature(RF_prok_rich_T2, taxa_prok_no_norm_sub) +
              labs(title = "Bacterial\nOTU importance"),
            PlotHeat(df_melt_T2),
            labels = c("A","B","", "C", "D", ""),
            widths = c(0.6, 1, 0.45, 0.6, 1, 0.45),
            align = "hv",
            ncol = 3,
            nrow =2)

Fig_X_RF_Bul


Fig_X_RF_Bul_Error <-
  ggarrange(PlotError(RF_prok_rich) +
              annotate("text", x=Inf, y = Inf, 
                       label= paste("% Explained var.:", round(last(RF_prok_rich$rsq*100),2)), size=2.5, vjust=1, hjust=1) +
              annotate("text", x=Inf, y = Inf,
                       label=paste("MSE:", round(last(RF_prok_rich$mse), 2)), size=2.5, vjust=3, hjust=1) +
              annotate("text", x=Inf, y = Inf, 
                       label = paste("italic(p) ==", round(perm_RF_prok_rich$pValue, 4)), parse = TRUE, size=2.5, vjust=4, hjust=1) +
              labs(title = "Bulleribasidiaceae\nrichness"),
            
            PlotError(RF_prok_nmds1) +
              annotate("text", x=Inf, y = Inf, 
                       label= paste("% Explained var.:", round(last(RF_prok_nmds1$rsq*100),2)), size=2.5, vjust=1, hjust=1) +
              annotate("text", x=Inf, y = Inf,
                       label=paste("MSE:", round(last(RF_prok_nmds1$mse), 2)), size=2.5, vjust=3, hjust=1) +
              annotate("text", x=Inf, y = Inf, 
                       label = paste("italic(p) ==", round(perm_RF_prok_nmds1$pValue, 4)), parse = TRUE, size=2.5, vjust=4, hjust=1) +
              labs(title = "Bulleribasidiaceae\nbeta diversity"),
            labels = c("A","B"),
            widths = c(1,1),
            align = "hv",
            ncol = 2, 
            nrow = 1)

Fig_X_RF_Bul_Error

Fig_X_RF_Bul_Error_Contr <-
  ggarrange(PlotError(RF_prok_rich_CR) +
              annotate("text", x=Inf, y = Inf, 
                       label= paste("% Explained var.:", round(last(RF_prok_rich_CR$rsq*100),2)), size=2.5, vjust=1, hjust=1) +
              annotate("text", x=Inf, y = Inf,
                       label=paste("MSE:", round(last(RF_prok_rich_CR$mse), 2)), size=2.5, vjust=3, hjust=1) +
              annotate("text", x=Inf, y = Inf, 
                       label = paste("italic(p) ==", round(perm_RF_prok_rich_CR$pValue, 4)), parse = TRUE, size=2.5, vjust=4, hjust=1)+
              labs(title = "Bulleribasidiaceae\nrichness"),
            
            PlotError(RF_prok_nmds1_CD) +
              annotate("text", x=Inf, y = Inf, 
                       label= paste("% Explained var.:", round(last(RF_prok_nmds1_CD$rsq*100),2)), size=2.5, vjust=1, hjust=1) +
              annotate("text", x=Inf, y = Inf,
                       label=paste("MSE:", round(last(RF_prok_nmds1_CD$mse), 2)), size=2.5, vjust=3, hjust=1) +
              annotate("text", x=Inf, y = Inf, 
                       label = paste("italic(p) ==", round(perm_RF_prok_nmds1_CD$pValue, 4)), parse = TRUE, size=2.5, vjust=4, hjust=1) +
              labs(title = "Bulleribasidiaceae\nbeta diversity"),
            labels = c("C","D"),
            widths = c(1,1),
            align = "hv",
            ncol = 2, 
            nrow = 1)

Fig_X_RF_Bul_Error_Contr

title1 = text_grob("Treatment",size=12, face=2)
title2 = text_grob("Control",size=12, face=2)

# *** FIGURE 2 ---------------------------------------------------------------------------------------------------
Fig_2_Errors = ggarrange(grid.arrange(Fig_X_RF_Bul_Error, top = title1),
                  grid.arrange(Fig_X_RF_Bul_Error_Contr, top = title2),
                  widths = c(1,1),
                  align = "hv", 
                  ncol = 1, 
                  nrow = 2)

Fig_2_Errors

PlotError(RF_prok_rich) +
  annotate("text", x=Inf, y = Inf,
           label=paste("Mean squared error:", round(last(RF_prok_rich$mse), 2)), size=2.5, vjust=1, hjust=1) +
  annotate("text", x=Inf, y = Inf, 
           label= paste("% Var explained:", round(last(RF_prok_rich$rsq*100),2)), size=2.5, vjust=3, hjust=1) +
  annotate("text", x=Inf, y = Inf, 
           label = paste("italic(p) ==", round(perm_RF_prok_rich$pValue, 4)), parse = TRUE, size=2.5, vjust=4, hjust=1)

PlotLine(RF_prok_rich, as.data.frame(physeq_fungi_no_norm_Bul[2])) + theme(legend.position = c(0.1, 0.9))

PlotFeature(RF_prok_rich, taxa_prok_no_norm_sub) + labs(title = "Top Bacterial\nOTUs explaining Fungal richness") 

PlotHeat(df_melt) 

write.dna(refseq(filterTaxa(physeq_prok_no_norm, rfe_prok_rich)),
          format="fasta", 
          colsep="", 
          file="bact_rich.fasta")

write.dna(refseq(filterTaxa(physeq_prok_no_norm, rfe_prok_nmds1)),
          format="fasta", 
          colsep="", 
          file="bact_nmds1.fasta")
