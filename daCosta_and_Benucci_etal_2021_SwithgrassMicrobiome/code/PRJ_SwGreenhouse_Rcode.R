# Project: GLBRC greenhouse Pedro 

# Aithors: Chou Ming-Yi, Benucci GMN, Beschoren da Costa P, Bonito G (preliminary order)
# Institution: Michigan State University

# WORKING ENVIRONMENT SETUP --------------------------------------------------------------------------------
options(scipen = 999)
options(max.print = 100000000)

#rm(list = ls(all=TRUE)) # removes all variables in the global environment so you start fresh
Sys.time() # prints out the time and date you ran the code
set.seed(1977) #to make reproduceble results
detach(package:phyloseq, unload = TRUE) #to reorder packages
search() #to search into the environment
#rm(list= ls()[!(ls() %in% c('keepThis','andThis'))]) #remove all but ...
session.info() #to get session information

# Citation resourches
citation()
version$version.string
devtools::session_info()
print(citation("MASS"), style = "text")
library(purrr)
c("vegan", "phyloseq", "ape") %>%
  map(citation) %>%
  print(style = "text")


# NECESSARY PACKAGES ---------------------------------------------------------------------------------------
library(phyloseq)
library(Biostrings)
library(ape)

library(tidyverse)
library(dplyr)
library(tidyr)

library(ggplot2)
library(ggpubr)
library(magrittr)

library(decontam)

library(stringr)
library(purrr)

library(vegan) 
library(jtools)
library(interactions)

library(multcompView)
library(multcomp)
library(sjPlot)
library(car)

library(caret)

library(grid)
library(gridExtra)

library(Boruta)
library(randomForest)
library(rfUtilities)

library(agricolae)
library(ggfortify)

library(rcompanion)

# Palettes -------------------------------------------------------------------------------------------------

Pal_geno <- c("#F2AA91","#F6511D","#FFB400","#00A6ED","#7FB800","#0D2C54")
Pal_soil <- c("#F6511D","#FFB400","#00A6ED","#7FB800")

pie(rep(1, length(Pal_geno)), labels = sprintf("%d (%s)",
   seq_along(Pal_geno),Pal_geno), col = Pal_geno)


# LOAD ALL LIBRARIES ---------------------------------------------------------------------------------------
# Fungi ----------------------------------------------------------------------------------------------------
Lib_all_OTU_ITS <-
  read.delim(
    "/home/gian/Documents/GREENHOSUE_glbrc_project/Data/Joint_Fun_otu_tables/FOTU_table_ITS_UPARSE.txt",
    row.names = 1
  )

Lib_all_mapping_ITS <-
  read.delim(
    "/home/gian/Documents/GREENHOSUE_glbrc_project/Data/Joint_Fun_otu_tables/Lib_all_Fun_mapping.txt",
    row.names = 1
  )

All_otus_ITS <-
  readDNAStringSet(
    "/home/gian/Documents/GREENHOSUE_glbrc_project/Data/Joint_Fun_otu_tables/FOTUs.fasta",
    format = "fasta",
    seek.first.rec = TRUE,
    use.names = TRUE
  )

levels(factor(Lib_all_mapping_ITS$Sample_name))

# Prokaryotes ----------------------------------------------------------------------------------------------
Lib_all_OTU_16S <-
  read.delim(
    "/home/gian/Documents/GREENHOSUE_glbrc_project/Data/Joint_Bacterial_otu_tables/otu_table_16S_UPARSE.txt",
    row.names = 1
  )
Lib_all_mapping_16S <-
  read.delim(
    "/home/gian/Documents/GREENHOSUE_glbrc_project/Data/Joint_Bacterial_otu_tables/Lib_all_mapping.txt",
    row.names = 1
  )
All_otus_16S <-
  readDNAStringSet(
    "/home/gian/Documents/GREENHOSUE_glbrc_project/Data/Joint_Bacterial_otu_tables/otus.fasta",
    format = "fasta",
    seek.first.rec = TRUE,
    use.names = TRUE
  )

levels(as.factor(Lib_all_mapping_16S$Sample_name))

# ***********************************************************************************************-----------
# ADJUST MAPPING FILE --------------------------------------------------------------------------------------
classify_mock_control_sample<- function (Lib_mapping) {
  Control_list<-c("Control1","Control2","Control3","Control4","Control5","Control6",
                  "Control_neg_noseq", "Control_noseq")
  Mock_list<- c("Mock1","Mock2", "Mock_resuspendBar","Mock_1ul", "Mock2_176_SpikeIn",
                "Mock2_376_Line6_G3_SpikeIn","Mock2_MMPRNT-1953_SpikeIn")
  Lib_mapping$band_1_RAW_concentration <-
    if_else(is.na(Lib_mapping$band_1_RAW_concentration), 0, Lib_mapping$band_1_RAW_concentration)
  Lib_mapping$band_2_RAW_concentration <-
    if_else(is.na(Lib_mapping$band_2_RAW_concentration), 0, Lib_mapping$band_2_RAW_concentration)
  Lib_rownames<-(rownames(Lib_mapping))
  Lib_mapping <- 
    mutate(Lib_mapping, 
           Control_mock_sample = if_else(Lib_mapping$Sample_name%in%c(Control_list), "Control",
                                        if_else(Lib_mapping$Sample_name%in%c(Mock_list), "Mock","Sample" ))) %>%
    mutate(AmpliconCon = band_1_RAW_concentration + band_2_RAW_concentration)
  rownames(Lib_mapping)<-Lib_rownames 
  return(Lib_mapping)
}

Lib_all_mapping_ITS_new <-
  classify_mock_control_sample(Lib_all_mapping_ITS)
Lib_all_mapping_16S_new <-
  classify_mock_control_sample(Lib_all_mapping_16S)


head(Lib_all_mapping_ITS_new)

Lib_all_mapping_ITS_new %>%
  select(Sample_name, Control_mock_sample, Control_mock_sample)


# Adjust PPCR concentration --------------------------------------------------------------------------------
# Keep DNA concentration only of bands between 250 and 750bp 
# as amplicon_concentration, remove_non_target_bands.

ModConc <- function(map){
  map_rownames <-rownames(map) 
  map1 <- map[c("band_1_size","band_2_size","band_1_RAW_concentration","band_2_RAW_concentration")]
  #first set size those below 300 to or above 550 zero. if band size is zero, 
  #then concentration changes to zero. Finaly, then  unite() the 2 columns of
  #band size and the 2 columns of band cocnentration
  map1[is.na(map1)] = 0 #turns NA to zero; input for next command cannot be NA
  map1[map1$band_1_size<300 | map1$band_1_size>600,]$band_1_RAW_concentration<-NA  # this turns primer dimer bands into NA values for band
  map1[map1$band_2_size<300 | map1$band_2_size>600,]$band_2_RAW_concentration<-NA  # this turns primer dimer bands into NA values for band
  map1[map1$band_1_size<300 | map1$band_1_size>600,]$band_1_size<-NA # this turns primer dimer bands into NA values for band size
  map1[map1$band_2_size<300 | map1$band_2_size>600,]$band_2_size<-NA  # this turns primer dimer bands into NA values for band size
  map1[is.na(map1)] = 0 #turns NA to  "" to allow unite to mix both columns
  map2 <- map1 %>% 
    mutate(
    amplicon_size = pmax(map1$band_1_size, map1$band_2_size),
    amplicon_concentration =  pmax(map1$band_1_RAW_concentration,  map1$band_2_RAW_concentration))
  #turn empy value to 1, the aproximated detection limit
  map2[map2$amplicon_size==0,]$amplicon_size <- 1 
  #turn empy value to 1, the aproximated detection limit
  map2[map2$amplicon_concentration==0,]$amplicon_concentration <- 1
  head(map2) %T>% print()
  map_last <-
    cbind(map, map2[,c(5,6)])
return(map_last)
}


Lib_all_mapping_ITS_mod <- 
  ModConc(Lib_all_mapping_ITS_new)
head(Lib_all_mapping_ITS_mod)
Lib_all_mapping_16S_mod <- 
  ModConc(Lib_all_mapping_16S_new)
head(Lib_all_mapping_16S_mod)


# GENERATE PHYLOSEQ FILES ----------------------------------------------------------------------------------
# generates fungal phyloseq object
Lib_all_ITS <-
  phyloseq(
    otu_table(Lib_all_OTU_ITS, taxa_are_rows = TRUE),
    sample_data(Lib_all_mapping_ITS_mod),
    All_otus_ITS)

Lib_all_ITS
head(sample_data(Lib_all_ITS))

table(Lib_all_ITS@sam_data$Control_mock_sample)

# generates phyloseq object
Lib_all_16S <-
  phyloseq(
    otu_table(Lib_all_OTU_16S, taxa_are_rows = TRUE),
    sample_data(Lib_all_mapping_16S_mod),
    All_otus_16S)

Lib_all_16S
head(sample_data(Lib_all_16S))

table(Lib_all_16S@sam_data$Control_mock_sample)

# EXTRACT PROJECT PHYLOSEQ OBJETCS ------------------------------------------------------------------------ 
# >>> GREEENHOUSE STUDY ------------------------------------------------------------------------------------

# Fungi
greenhouse_ITS <- subset_samples(Lib_all_ITS, Project_name%in%c("Greenhouse") & 
                                   Library_N%in%c("Lib_6", "Lib_8"))
otu_table(greenhouse_ITS) <- 
  otu_table(greenhouse_ITS)[which(rowSums(otu_table(greenhouse_ITS)) > 0),] 
greenhouse_ITS

count(as.data.frame(as.matrix(sample_data(greenhouse_ITS))), 
      Project_name, Library_N, DNA_plate, Control_mock_sample)

head(greenhouse_ITS@sam_data)
table(greenhouse_ITS@sam_data$Control_mock_sample)


# Extracting controls ands mock samples
greenhouse_ITS_contr_mock <- 
  subset_samples(Lib_all_ITS, 
                 DNA_plate%in%c("PLATE_3_ROOT","PLATE_4_ROOT","PLATE_6_ROOT","PLATE_11_SOIL") &
                   Project_name%in%c("Quality"))
otu_table(greenhouse_ITS_contr_mock) <- 
  otu_table(greenhouse_ITS_contr_mock)[which(rowSums(otu_table(greenhouse_ITS_contr_mock)) > 0),] 

greenhouse_ITS_contr_mock

count(as.data.frame(as.matrix(sample_data(greenhouse_ITS_contr_mock))), 
      Project_name, Library_N, DNA_plate, Control_mock_sample)

greenhouse_ITS_contr_mock@sam_data

# Bacteria
greenhouse_16S <- subset_samples(Lib_all_16S, Project_name%in%c("Greenhouse") &
                                   Library_n%in%c("Lib_2", "Lib_4"))
otu_table(greenhouse_16S) <- 
  otu_table(greenhouse_16S)[which(rowSums(otu_table(greenhouse_16S)) > 0),] 
greenhouse_16S

count(as.data.frame(as.matrix(sample_data(greenhouse_16S))), 
      Project_name, Library_n, DNA_plate, Control_mock_sample)

head(greenhouse_16S@sam_data)


# Extracting controls ands mock samples
greenhouse_16S_contr_mock <- 
  subset_samples(Lib_all_16S, 
                 DNA_plate%in%c("PLATE_3_ROOT","PLATE_4_ROOT","PLATE_6_ROOT","PLATE_11_SOIL") &
                   Project_name%in%c("Quality"))
otu_table(greenhouse_16S_contr_mock) <- 
  otu_table(greenhouse_16S_contr_mock)[which(rowSums(otu_table(greenhouse_16S_contr_mock)) > 0),] 

greenhouse_16S_contr_mock

count(as.data.frame(as.matrix(sample_data(greenhouse_16S_contr_mock))), 
      Project_name, Library_n, DNA_plate, Control_mock_sample)

greenhouse_16S_contr_mock@sam_data


# Remove re-do samples with lowest read number -------------------------------------------------------------
head(greenhouse_ITS@sam_data)
greenhouse_ITS@sam_data$Description
duplicated(greenhouse_ITS@sam_data$Description)

table(duplicated(greenhouse_ITS@sam_data$Description))

head(greenhouse_16S@sam_data)
greenhouse_16S@sam_data$Description
duplicated(greenhouse_16S@sam_data$Description)

# Extracting samples with the highest number of reads
RemoveDuplicates <- function(physeq){
  print("Duplicated samples")
  print(table(duplicated(physeq@sam_data$Description)))
  physeq <- 
    subset_samples(physeq, Description%in%physeq@sam_data[
      duplicated(physeq@sam_data$Description), ]$Description)
  df <- 
    data.frame(sample_sums(physeq), row.names = row.names(physeq@sam_data))
  colnames(df) <- "ReadNo"
  df$SampleID <- row.names(physeq@sam_data)
  df$Description <- physeq@sam_data$Description
  df <- df[order(df$Description),]
  df %T>% print()
  res <- 
    df %>% 
    group_by(Description) %>%
    top_n(1, -abs(ReadNo )) %>% 
    as.data.frame()
  res <- res$SampleID
  return(res)
}


RemoveDuplicates(greenhouse_ITS)

greenhouse_ITS_clean <-
  subset_samples(
    greenhouse_ITS,
    !sample_names(greenhouse_ITS) %in% RemoveDuplicates(greenhouse_ITS))
duplicated(greenhouse_ITS_clean@sam_data$Description)

RemoveDuplicates(greenhouse_16S)

greenhouse_16S_clean <-
  subset_samples(
    greenhouse_16S,
    !sample_names(greenhouse_16S) %in% RemoveDuplicates(greenhouse_16S))
duplicated(greenhouse_16S_clean@sam_data$Description)

# since there are samples replicated multiple times I run it 3 times
greenhouse_16S_clean <-
  subset_samples(
    greenhouse_16S_clean,
    !sample_names(greenhouse_16S_clean) %in% RemoveDuplicates(greenhouse_16S_clean))
duplicated(greenhouse_16S_clean@sam_data$Description)

greenhouse_16S_clean <-
  subset_samples(
    greenhouse_16S_clean,
    !sample_names(greenhouse_16S_clean) %in% RemoveDuplicates(greenhouse_16S_clean))
duplicated(greenhouse_16S_clean@sam_data$Description)


# SEPARATE MOCK form SAMPLES -------------------------------------------------------------------------------
greenhouse_ITS_contr <-
  subset_samples(greenhouse_ITS_contr_mock, Control_mock_sample%in%c("Control"))
greenhouse_ITS_contr
greenhouse_ITS_contr@sam_data
count(as.data.frame(as.matrix(sample_data(greenhouse_ITS_contr))), 
      Project_name, Library_N, DNA_plate, Control_mock_sample)

greenhouse_ITS_mock <-
  subset_samples(greenhouse_ITS_contr_mock, Control_mock_sample%in%c("Mock"))
greenhouse_ITS_mock
greenhouse_ITS_mock@sam_data
count(as.data.frame(as.matrix(sample_data(greenhouse_ITS_mock))), 
      Project_name, Library_N, DNA_plate, Control_mock_sample)


greenhouse_16S_contr <-
  subset_samples(greenhouse_16S_contr_mock, Control_mock_sample%in%c("Control"))
greenhouse_16S_contr
greenhouse_16S_contr@sam_data
count(as.data.frame(as.matrix(sample_data(greenhouse_16S_contr))), 
      Project_name, Library_n, DNA_plate, Control_mock_sample)

greenhouse_16S_mock <-
  subset_samples(greenhouse_16S_contr_mock, Control_mock_sample%in%c("Mock"))
greenhouse_16S_mock
greenhouse_16S_mock@sam_data
count(as.data.frame(as.matrix(sample_data(greenhouse_16S_mock))), 
      Project_name, Library_n, DNA_plate, Control_mock_sample)


# GENERATE PHYLOSEQ OBJECTS FOR GREENHOUSE -----------------------------------------------------------------
physeq_fungi_gs <-
  merge_phyloseq(greenhouse_ITS_clean, greenhouse_ITS_contr)
physeq_fungi_gs

count(as.data.frame(as.matrix(sample_data(physeq_fungi_gs))), 
      Project_name, Library_N, DNA_plate, Control_mock_sample)
physeq_fungi_gs@sam_data

physeq_bact_gs <-
  merge_phyloseq(greenhouse_16S_clean, greenhouse_16S_contr)
physeq_bact_gs

count(as.data.frame(as.matrix(sample_data(physeq_bact_gs))), 
      Project_name, Library_n, DNA_plate, Control_mock_sample)
physeq_bact_gs@sam_data

# PCR BAD ANALYSES -----------------------------------------------------------------------------------------
# GEN BAND plot --------------------------------------------------------------------------------------------
quiaxcel_band_plot<-function(Lib_mapping){
  # true bands must be between 300 and 550 bp; primer dimer concentrtion is less than 15ng/ul
  Band_size_and_cocentration1 <- Lib_mapping%>%
    select(band_1_size, band_1_RAW_concentration, Control_mock_sample)%>%
    rename(band_size=band_1_size, band_concentration=band_1_RAW_concentration)
  Band_size_and_cocentration2 <- Lib_mapping%>%
    select(band_2_size, band_2_RAW_concentration, Control_mock_sample)%>%
    rename(band_size=band_2_size, band_concentration=band_2_RAW_concentration)
  # generate a long format
  Band_size_and_cocentration <-
    bind_rows(Band_size_and_cocentration1, Band_size_and_cocentration2)
  Band_size_and_cocentration[is.na(Band_size_and_cocentration)] <- NA
  head(Band_size_and_cocentration) %T>% print()
  plot_band <-
    ggplot(Band_size_and_cocentration, 
           aes(x=band_concentration, y=band_size, color=Control_mock_sample))+
           geom_point() +
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.x = element_text(angle =0, size = 8, hjust = 0.5, vjust = 1), 
          axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8),
          legend.position="right") +
    guides(color = guide_legend(title = "Sample type", ncol = 1))
return(plot_band)
  }


#view plot of band DNA concentration per band size
quiaxcel_band_plot(Lib_all_mapping_ITS_new) +
  labs(title = "PCR band concentration/size ITS",
       x= "Concentration (ng/ul)", 
       y="Size (bp)")

quiaxcel_band_plot(Lib_all_mapping_16S_new) +
  labs(title = "PCR band concentration/size 16S",
       x= "Concentration (ng/ul)", 
       y="Size (bp)")


# Plotting size histograms ---------------------------------------------------------------------------------
# There are some NA that will prevent the size to be plotted 
PlotHist <- function(df, width){
  ggplot(df, aes(x=band_1_size)) +
    geom_histogram(binwidth=width, colour="grey80", fill="grey80") +
    #geom_vline(aes(xintercept=mean(band_1_size, na.rm=TRUE)), 
    #                color="red", linetype="dashed", size=0.8) +
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.x = element_text(angle =0, size = 8, hjust = 0.5, vjust = 1), 
          axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8),
          legend.position="none") -> plot_dist
  return(plot_dist)  
}


PlotHist(Lib_all_mapping_ITS_new, 1) +
  labs(title="Distribution of band sizes", x="band size (bp)", y="Sample number")

PlotHist(Lib_all_mapping_16S_new, 1) +
  xlim(0, 1000) +
  labs(title="Distribution of band sizes", x="band size (bp)", y="Sample number")


# *** Band plots  ------------------------------------------------------------------------------------------
ggarrange(
  ggarrange(
    quiaxcel_band_plot(Lib_all_mapping_ITS_new) +
      ylim(0, 2500) +
            labs(title = "PCR band ITS",
                 x= "Concentration (ng/ul)", 
                 y="Size (bp)"),
    quiaxcel_band_plot(Lib_all_mapping_16S_new) +
      ylim(0, 2500) +
            labs(title = "PCR band 16S",
                 x= "Concentration (ng/ul)", 
                 y="Size (bp)"),
          labels = c("A","B"),
          align = "hv",
          ncol =2,
          nrow = 1,
          common.legend = TRUE,
          legend=c("right")),
  ggarrange(
    PlotHist(Lib_all_mapping_ITS_new, 1) +
            labs(title="", x="band size (bp)", y="Sample number"),
    PlotHist(Lib_all_mapping_16S_new, 1) +
            xlim(0, 1000) +
            labs(title="", x="band size (bp)", y="Sample number"),
          labels = c("C","D"),
          align = "hv",
          ncol = 2, 
          nrow = 1,
          legend = c("right")),
  ncol = 1,
  nrow = 2) -> band_plot

band_plot

# DECONTAMINATIONS -----------------------------------------------------------------------------------------

# INSPECTING LIBRARY SIZES ---------------------------------------------------------------------------------
sample_data(physeq_fungi_gs)$LibrarySize <-
  sample_sums(physeq_fungi_gs)
sample_data(physeq_fungi_gs)$Index <-
  seq(nrow(sample_data(physeq_fungi_gs)))
sample_data(physeq_fungi_gs)$is.neg <-
  sample_data(physeq_fungi_gs)$Control_mock_sample == "Control"
head(sample_data(physeq_fungi_gs))

sample_data(physeq_bact_gs)$LibrarySize <-
  sample_sums(physeq_bact_gs)
sample_data(physeq_bact_gs)$Index <-
  seq(nrow(sample_data(physeq_bact_gs)))
sample_data(physeq_bact_gs)$is.neg <-
  sample_data(physeq_bact_gs)$Control_mock_sample == "Control"
head(sample_data(physeq_bact_gs))

write.csv(physeq_fungi_gs@sam_data, "mapping_good_ITS.csv")
write.csv(physeq_bact_gs@sam_data, "mapping_good_16s.csv")


# Plotting sample depth ------------------------------------------------------------------------------------
PlotDepth <- function(physeq){
  df <-
    as(sample_data(physeq), "matrix")
  df <- 
    as.data.frame(df)[,c(26:28)]
  # reconvert to numeric
  df$LibrarySize <- as.numeric(as.character(df$LibrarySize))
  df$Index <- as.numeric(as.character(df$Index))
  # order
  df <- df[order(df$LibrarySize), ]
  df$Index <- seq(nrow(df))
  # inspect
  str(df) %T>% print()
  head(df) %T>% print()
  ggplot(data=df, aes(x=Index, y=LibrarySize, color=is.neg)) +
    geom_point(alpha =0.7) +
    theme_classic() +
    scale_colour_manual("Negative control", values = c("grey", "red")) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.x = element_text(angle =0, size = 8, hjust = 0.5, vjust = 1), 
          axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8)) -> plot_dist
  return(plot_dist)  
}

PlotDepth(physeq_fungi_gs) +
  labs(title="Fungi", 
       subtitle = "Samples read depth", 
       x="Sample index", 
       y="Read number")

PlotDepth(physeq_bact_gs) +
  labs(title="Prokaryotes", 
       subtitle = "Samples read depth",
       x="Sample index",
       y="Read number")


# DECONTAMINATION ------------------------------------------------------------------------------------------
head(physeq_fungi_gs@sam_data)
physeq_fungi_gs@sam_data$DNA_extraction_yield # do not know why has many NA
physeq_fungi_gs@sam_data[is.na(physeq_fungi_gs@sam_data$DNA_extraction_yield),]
physeq_fungi_gs@sam_data$AmpliconCon

# Error: conc must be positive numeric.
physeq_fungi_gs@sam_data$AmpliconCon <- physeq_fungi_gs@sam_data$AmpliconCon + 0.000001

contam_fungi_gs <- isContaminant(physeq_fungi_gs,
                            method="either", 
                            neg="is.neg",
                            batch = "DNA_plate",
                            batch.combine = "fisher",
                            conc="AmpliconCon",
                            threshold=c(0.1, 0.5))
                                                    
table(contam_fungi_gs$contaminant)


physeq_bact_gs@sam_data$AmpliconCon <- physeq_bact_gs@sam_data$AmpliconCon + 0.000001

contam_bact_gs <- isContaminant(physeq_bact_gs,
                                 method="either", 
                                 neg="is.neg",
                                 batch = "DNA_plate",
                                 batch.combine = "fisher",
                                 conc="AmpliconCon",
                                 threshold=c(0.1, 0.5))

table(contam_bact_gs$contaminant) 


# plotting contaminant OTUs --------------------------------------------------------------------------------
PlotContam <- function(df, contam){
  # Make phyloseq object of presence-absence in negative controls and true samples
  physeq_pa <- transform_sample_counts(df, function(abund) 1*(abund>0))
  physeq_pa_neg <- subset_samples(physeq_pa, Control_mock_sample%in%c("Control"))
  physeq_pa_pos <- subset_samples(physeq_pa, Control_mock_sample%in%c("Sample"))
  # Make data.frame of prevalence in positive and negative samples
  df_contam <- data.frame(pa.pos=taxa_sums(physeq_pa_pos), 
                          pa.neg=taxa_sums(physeq_pa_neg),
                          contaminant=contam$contaminant, 
                          Pvalue=contam$p)
  head(df_contam) %T>% print()
  # plotting 
  ggplot(data=df_contam, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
    geom_point(size=0.8, alpha=0.7) +
    labs(x="Prevalence in negative controls", y="Prevalence in true samples") +
    theme_classic() +
    scale_colour_manual("Contaminant OTUs", values = c("grey", "red")) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.x = element_text(angle =0, size = 8, hjust = 0.5, vjust = 1), 
          axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8))-> plot_cont
  return(plot_cont)
}

PlotContam(physeq_fungi_gs, contam_fungi_gs)
PlotContam(physeq_bact_gs, contam_bact_gs) 

table(contam_fungi_gs$contaminant)
table(contam_bact_gs$contaminant)

table(contam_fungi_gs$contaminant)/sum(table(contam_fungi_gs$contaminant))*100
table(contam_bact_gs$contaminant)/sum(table(contam_bact_gs$contaminant))*100

# *** FIGURE S1 - dcontam ----------------------------------------------------------------------------------
ggarrange(
  ggarrange(PlotDepth(physeq_fungi_gs) +
              labs(title="Fungi", 
                   subtitle = "Samples read depth", 
                   x="Sample index", 
                   y="Read number"),
            PlotDepth(physeq_bact_gs) +
              labs(title="Prokaryotes", 
                   subtitle = "Samples read depth",
                   x="Sample index",
                   y="Read number"),
            labels = c("A","B"),
            widths = c(1,1,1,1),
            align = "hv" ,
            ncol = 2, 
            nrow = 1, 
            common.legend = TRUE, 
            legend = c("bottom")),
  ggarrange(
    PlotContam(physeq_fungi_gs, contam_fungi_gs) +
      labs(subtitle="Contaminants"),
    PlotContam(physeq_bact_gs, contam_bact_gs) +
      labs(subtitle="Contaminants"),
            widths = c(1,1),
            labels = c("C","D"),
            align = "hv" ,
            ncol = 2, 
            nrow = 1,
            common.legend = TRUE,
            legend = c("bottom")),
  widths =  c(1, 1.2),
  ncol = 1, 
  nrow = 2) -> Fig_S1

Fig_S1


# REMOVING CONTAMINANTS ------------------------------------------------------------------------------------
# function to remove bad taxa
remove_taxa = function(physeq, badTaxa) {
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

physeq_fungi_gs_clean <-
  remove_taxa(physeq_fungi_gs, rownames(subset(
    contam_fungi_gs, contaminant %in% c("TRUE")
  )))

physeq_bact_gs_clean <-
  remove_taxa(physeq_bact_gs, rownames(subset(
    contam_bact_gs, contaminant %in% c("TRUE")
  )))

# MAPPING FILES FOR SRA -----------------------------------------------------------------------------
write.csv(sample_data(physeq_fungi_gs_clean), "physeq_fungi_gs_clean.csv")
write.csv(sample_data(physeq_bact_gs_clean), "physeq_bact_gs_clean.csv")

sample_names(physeq_bact_gs_clean)
as.vector(physeq_bact_gs_clean@otu_table[, c("Amp489")])

physeq_bact_gs_clean@sam_data[c("Amp489"),]


# REMOVING ALL CONTROLS ----------------------------------------------------------------------------------- 
# fungi
physeq_fungi_gs_clean <-
  subset_samples(physeq_fungi_gs_clean, Control_mock_sample %in% c("Sample"))
otu_table(physeq_fungi_gs_clean) <-
  otu_table(physeq_fungi_gs_clean)[which(rowSums(otu_table(physeq_fungi_gs_clean)) > 0),]
physeq_fungi_gs_clean
head(physeq_fungi_gs_clean@otu_table)


# bacteria
physeq_bact_gs_clean <-
  subset_samples(physeq_bact_gs_clean, Control_mock_sample %in% c("Sample"))
otu_table(physeq_bact_gs_clean) <-
  otu_table(physeq_bact_gs_clean)[which(rowSums(otu_table(physeq_bact_gs_clean)) > 0),]
physeq_bact_gs_clean
head(physeq_bact_gs_clean@otu_table)

# EXPORT SEQUENCES for CONSTAX -----------------------------------------------------------------------------
write.dna(refseq(physeq_fungi_gs_clean), format="fasta", file = "PRJ_Greenhouse_seq_ITS.fasta", colsep="")
write.dna(refseq(physeq_bact_gs_clean), format="fasta", file = "PRJ_Greenhouse_seq_16S.fasta", colsep="")

# Importing CONSTAX taxonomies -----------------------------------------------------------------------------

# FUNGI ----------------------------------------------------------------------------------------------------
taxonomy_ITS <-
  read.delim(
    "constax_taxonomy_ITS.txt",
    header = TRUE,
    row.names = 1,
    sep = "\t")

head(taxonomy_ITS)
taxonomy_ITS[1:100, ]
dim(taxonomy_ITS)

table(taxonomy_ITS$High_level_taxonomy)
table(taxonomy_ITS$Kingdom)
# NOTE! There are still some mock among the samples. We have some cross-talks we need to quantify!
untarget_uparse_ITS <- c("Alveolata","Amoebozoa","Choanoflagellozoa","Metazoa","Cercozoa", "Protista",
                     "Rhizaria","Stramenopila", "Viridiplantae", "MockK", "Mockc","Mockk","Mockp")

apply(taxonomy_ITS, 2, function(x) which(x %in% untarget_uparse_ITS))
apply(taxonomy_ITS, 2, function(x) which(x %in% c("MockK")))


# Remove non-target taxa and non-classified taxa that did not have any hit
# when blasted against UNITE at 60% conevrage and identity. ALso, remove 
# fungal OTUs that were classified as Fungi bu had low hit coverage < 60%
# and identity < 60% 

taxonomy_ITS_bad <-
  subset(
    taxonomy_ITS,
    High_level_taxonomy %in% untarget_uparse_ITS |
      Kingdom %in% untarget_uparse_ITS |
      High_level_taxonomy %in% c("") &
      Kingdom %in% c(""))

dim(taxonomy_ITS_bad)

taxonomy_ITS_filt <-
  taxonomy_ITS[!(rownames(taxonomy_ITS) %in% rownames(taxonomy_ITS_bad)), ]
taxonomy_ITS_filt <- 
  data.frame(OTU_ID = rownames(taxonomy_ITS_filt), taxonomy_ITS_filt)
dim(taxonomy_ITS_filt)
head(taxonomy_ITS_filt)
table(taxonomy_ITS_filt$Kingdom)
table(taxonomy_ITS_filt$High_level_taxonomy)

# looking at non-classified OTUs:
# Unclassified OTUs that hit Fungi at > 60%
dim(subset(
  taxonomy_ITS,
  High_level_taxonomy %in% c("Fungi") &
    Kingdom %in% c("")))

# Fungal OTUs that were classified as Fungi bu had low hit
# coverage < 60% and identity < 60% 
dim(subset(
  taxonomy_ITS,
  Kingdom %in% c("Fungi") &
    HL_hit_query_cover < 60 |
    HL_hit_percent_id < 60))

# How many isolates we were able to detect with the MiSeq?
subset(taxonomy_ITS,!Isolate %in% "")


# BACTERIA -------------------------------------------------------------------------------------------------
taxonomy_16S <-
  read.delim(
    "constax_taxonomy_16S.txt",
    header = TRUE,
    row.names = 1, 
    sep = "\t")

head(taxonomy_16S)
taxonomy_16S[1:100, ]
dim(taxonomy_16S)

taxonomy_16S <- 
  taxonomy_16S[-8]
taxonomy_16S <- 
  data.frame(OTU_ID = rownames(taxonomy_16S), taxonomy_16S)

# cleaning taxonomy labels 
taxonomy_16S <- 
  map_dfc(
    1:ncol(taxonomy_16S),
    ~ str_remove_all(taxonomy_16S[,.x],"_1")) %>%
  as.data.frame()

# readding labels
colnames(taxonomy_16S) <-
  c("OTU_ID",
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
    "Isolate",
    "Isolate Isolate_percent_id",
    "Isolate_query_cover",
    "High_level_taxonomy",
    "HL_hit_percent_id",
    "HL_hit_query_cover")

# remember! now there are all characters
str(taxonomy_16S)
rownames(taxonomy_16S) <- taxonomy_16S$OTU_ID

# Looking at non-target
table(taxonomy_16S$High_level_taxonomy)
table(taxonomy_16S$Kingdom)
table(taxonomy_16S$Family)
table(taxonomy_16S$Isolate)
# NOTE! There are still some mock among the samples. We have some cross-talks we need to quantify!
untarget_uparse <- c("Mitochondria","Chloroplast")

apply(taxonomy_16S, 2, function(x) which(x %in% untarget_uparse))
apply(taxonomy_16S, 2, function(x) which(x %in% c("Mock")))

# Remove non-target taxa and non-classified taxa that did not have any hit
# when blasted against UNITE at 60% conevrage and identity. ALso, remove 
# fungal OTUs that were classified as Fungi bu had low hit coverage < 60%
# and identity < 60% 

taxonomy_16S_bad <- 
  subset(
    taxonomy_16S,
    High_level_taxonomy%in%untarget_uparse |
      Kingdom%in%untarget_uparse |
      Order%in%untarget_uparse |
      Family%in%untarget_uparse |
      High_level_taxonomy%in%c("") &
      Kingdom%in%c(""))

dim(taxonomy_16S_bad)

apply(taxonomy_16S_bad, 2, function(x) which(x %in% untarget_uparse))

taxonomy_16S_filt <-
  taxonomy_16S[!(rownames(taxonomy_16S) %in% rownames(taxonomy_16S_bad)), ]
dim(taxonomy_16S_filt)

apply(taxonomy_16S_filt, 2, function(x) which(x %in% untarget_uparse))

head(taxonomy_16S_filt)
table(taxonomy_16S_filt$Kingdom)
table(taxonomy_16S_filt$High_level_taxonomy)
rownames(taxonomy_16S_filt) <- taxonomy_16S_filt$OTU_ID


# Extract Chloroplast and Mitocondria 16S reads ------------------------------------------------------------
taxonomy_16S$HL_hit_percent_id <- as.numeric(taxonomy_16S$HL_hit_percent_id)
rownames(taxonomy_16S) <- taxonomy_16S$OTU_ID

taxonomy_Chloro_Mito <-
    subset(taxonomy_16S, 
      Order%in%untarget_uparse |
      Family%in%untarget_uparse |
      High_level_taxonomy%in%untarget_uparse |
      Kingdom%in%untarget_uparse)

taxonomy_Chloro_Mito

physeq_Chloro_Mito <- phyloseq(otu_table(physeq_bact_gs_clean, taxa_are_rows = TRUE),
                               sample_data(physeq_bact_gs_clean),
                               tax_table(as.matrix(taxonomy_Chloro_Mito)),
                               refseq(physeq_bact_gs_clean)) 

physeq_Chloro_Mito
head(physeq_Chloro_Mito@sam_data)


# Reformat Taxonomy ----------------------------------------------------------------------------------------
source("../R_functions/ReformatTaxonomy.R")

FinalizeTaxonomy <- function(constax){
  constax$Species <- 
    gsub(" sp ", "", constax$Species)
  constax[] = lapply(constax, blank2na, na.strings=c('','NA','na','N/A','n/a','NaN','nan'))
  lastValue <- function(x) tail(x[!is.na(x)], 1)
  constax$Genus <- as.character(constax$Genus)
  constax[which(is.na(constax$Genus) == FALSE),]$Genus <-
    paste(constax$Genus[is.na(constax$Genus) == FALSE], "sp.", sep = " ")
  last_taxons<- apply(constax[,c(1:8)], 1, lastValue)
  constax$BestMatch <- last_taxons
  constax$BestMatch <-
    gsub("_", " ", constax$BestMatch)
  constax$Taxonomy <-
    paste(constax$OTU_ID, constax$BestMatch, sep = "-")
  constax$BestMatch <- 
    gsub(" sp.", "", constax$BestMatch)
  constax$Genus <- 
    gsub(" sp.", "", constax$Genus)
  return(constax)
}


taxonomy_ITS_cor <- 
  FinalizeTaxonomy(taxonomy_ITS_filt)
taxonomy_ITS_cor[1:50, ]

taxonomy_16S_cor <- 
  FinalizeTaxonomy(taxonomy_16S_filt)
taxonomy_16S_cor[1:50, ]


# Extract mock sequences across samples --------------------------------------------------------------------

# Fungi
head(taxonomy_ITS)
levels(taxonomy_ITS$Isolate)

taxonomy_ITS_mock1 <-
  subset(taxonomy_ITS, Isolate %in% c("Mock_sp1", "Mock_sp2","Mock_sp3", "Mock_sp4",
                                    "Mock_sp5", "Mock_sp6","Mock_sp7", "Mock_sp8",
                                    "Mock_sp9", "Mock_sp10","Mock_sp11", "Mock_sp12"))
taxonomy_ITS_mock2 <-
  subset(taxonomy_ITS, High_level_taxonomy %in% "MockK")

taxonomy_ITS_mock <-
  full_join(
data.frame(OTU_ID = rownames(taxonomy_ITS_mock1),
           taxonomy_ITS_mock1),
data.frame(OTU_ID = rownames(taxonomy_ITS_mock2),
           taxonomy_ITS_mock2), by="OTU_ID")

taxonomy_ITS_mock
rownames(taxonomy_ITS_mock) <- taxonomy_ITS_mock$OTU_ID

physeq_fungi_gs_clean_mock <- phyloseq(otu_table(physeq_fungi_gs_clean, taxa_are_rows = TRUE),
                               sample_data(physeq_fungi_gs_clean),
                               tax_table(as.matrix(taxonomy_ITS_mock)),
                               refseq(physeq_fungi_gs_clean)) 

physeq_fungi_gs_clean_mock
physeq_fungi_gs_clean_mock@sam_data
physeq_fungi_gs_clean_mock@tax_table
physeq_fungi_gs_clean_mock@otu_table

# Sample Amp1369 is actually a mock sample, has to be removed.
physeq_fungi_gs_clean_mock@sam_data[rownames(physeq_fungi_gs_clean_mock@sam_data)%in%"Amp1369", ]
as.character(otu_table(subset_samples(physeq_fungi_gs_clean_mock, Sample_name%in%"305")))

physeq_fungi_gs_clean@otu_table <- 
  subset(physeq_fungi_gs_clean@otu_table, select = -c(Amp1369))
physeq_fungi_gs_clean@otu_table <- 
  physeq_fungi_gs_clean@otu_table[which(rowSums(physeq_fungi_gs_clean@otu_table) > 0),] 


# Bacteria 
head(taxonomy_16S)
taxonomy_16S$Kingdom
table(taxonomy_16S$Isolate)

taxonomy_16S_mock <-
  subset(taxonomy_16S, Isolate %in% c("Mock_Bacillus_subtilis6S0","Mock_Escherichia_coli6S_8",
                                      "Mock_Lactobacillus_fermentum6S_2","Mock_Listeria_monocytogenes6S_6",
                                      "Mock_Pseudomonas_aeruginosa6S_4","Mock_Salmonella_enterica6S_6"))
taxonomy_16S_mock

physeq_bact_gs_clean_mock <- phyloseq(otu_table(physeq_bact_gs_clean, taxa_are_rows = TRUE),
                                      sample_data(physeq_bact_gs_clean),
                                      tax_table(as.matrix(taxonomy_16S_mock)),
                                      refseq(physeq_bact_gs_clean)) 

physeq_bact_gs_clean_mock
physeq_bact_gs_clean_mock@sam_data
physeq_bact_gs_clean_mock@tax_table
physeq_bact_gs_clean_mock@otu_table

# MOCK COMMUNITY ANALYSIS ----------------------------------------------------------------------------------
physeq_fungi_gs_clean_mock %>%
  tax_table()

physeq_bact_gs_clean_mock %>%
  tax_table()

# IMPORT PLANT AND SOIL METADATA ---------------------------------------------------------------------------
metadata_gs <- read.csv(
    "/home/gian/Documents/GREENHOSUE_glbrc_project/Data/Greenhouse_data/greenhouse_plants.csv",
    header = TRUE, sep = ",")

colnames(metadata_gs)
head(metadata_gs)
dim(metadata_gs)

# Fungi
meta_fungi <- as(physeq_fungi_gs_clean@sam_data, "data.frame")
meta_fungi$Sample_ID <- rownames(meta_fungi) 
head(meta_fungi)
dim(meta_fungi)
str(meta_fungi)
names(meta_fungi)[names(meta_fungi) == "Sample_name"] <- "Pot"

metadata_fungi_gs <-
  merge(meta_fungi, metadata_gs, by="Pot")
rownames(metadata_fungi_gs) <- 
  metadata_fungi_gs$Sample_ID
dim(metadata_fungi_gs)
head(metadata_fungi_gs)

# Bacteria
meta_bact <- as(physeq_bact_gs_clean@sam_data, "data.frame")
meta_bact$Sample_ID <- rownames(meta_bact) 
head(meta_bact)
dim(meta_bact)
str(meta_bact)
names(meta_bact)[names(meta_bact) == "Sample_name"] <- "Pot"

metadata_bact_gs <-
  merge(meta_bact, metadata_gs, by="Pot")
rownames(metadata_bact_gs) <- 
  metadata_bact_gs$Sample_ID
dim(metadata_bact_gs)
head(metadata_bact_gs)

# Generating Phyloseq objetcs for GREENHOSUE STUDY ---------------------------------------------------------
# Look out for Mock OTUs first

dim(taxonomy_ITS_cor[rownames(taxonomy_ITS_cor)%in%rownames(taxonomy_ITS_mock), ])

physeq_fungi_new <- phyloseq(otu_table(physeq_fungi_gs_clean, taxa_are_rows = TRUE),
                       sample_data(metadata_fungi_gs),
                       tax_table(as.matrix(taxonomy_ITS_cor)),
                       refseq(physeq_fungi_gs_clean)) 
physeq_fungi_new
head(physeq_fungi_new@sam_data)
head(physeq_fungi_new@tax_table)
head(physeq_fungi_new@otu_table)

count(as.data.frame(as.matrix(sample_data(physeq_fungi_new))), 
      Genotype, Soil_location)

taxonomy_16S_cor[rownames(taxonomy_16S_cor)%in%rownames(taxonomy_16S_mock), ]

physeq_bact_new <- phyloseq(otu_table(physeq_bact_gs_clean, taxa_are_rows = TRUE),
                            sample_data(metadata_bact_gs),
                            tax_table(as.matrix(taxonomy_16S_cor)),
                            refseq(physeq_bact_gs_clean)) 
physeq_bact_new

physeq_bact_new <-
  remove_taxa(physeq_bact_new, rownames(taxonomy_16S_mock))
otu_table(physeq_bact_new) <-
  otu_table(physeq_bact_new)[which(rowSums(otu_table(physeq_bact_new)) > 0),]

physeq_bact_new
head(physeq_bact_new@sam_data)
head(physeq_bact_new@tax_table)
head(physeq_bact_new@otu_table)

# Match sample for bateria and fungi 
setdiff(physeq_fungi_new@sam_data$Pot, physeq_bact_new@sam_data$Pot) # present in first not in second object
setdiff(physeq_bact_new@sam_data$Pot, physeq_fungi_new@sam_data$Pot) # present in second not in first object

# use same pots for fungi and bacteria
sample_names(subset_samples(physeq_bact_new, Pot%in%"305"))

physeq_bact_new <-
  subset_samples(physeq_bact_new, Pot!="305")
physeq_bact_new@otu_table <- 
  physeq_bact_new@otu_table[which(rowSums(physeq_bact_new@otu_table) > 0),] 

physeq_bact_new

sort(sample_sums(physeq_bact_new))

# EXPORT DATASETS ------------------------------------------------------------------------------------------
saveRDS(physeq_fungi_new, "fungi_data.rds")
saveRDS(physeq_bact_new, "bacteria_data.rds")

write.table(sample_data(physeq_fungi_new), "filtered_data/its_metadata.txt", sep = "\t",
            quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(otu_table(physeq_fungi_new), "filtered_data/its_otu_table.txt", sep = "\t",
            quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(tax_table(physeq_fungi_new), "filtered_data/its_taxonomy.txt", sep="\t",
            quote=FALSE, row.names = TRUE, col.names = TRUE)
write.dna(refseq(physeq_fungi_new), format="fasta", file = "filtered_data/its_OTUs.fasta", colsep="")

write.table(sample_data(physeq_bact_new), "filtered_data/16s_metadata.txt", sep = "\t",
            quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(otu_table(physeq_bact_new), "filtered_data/16s_otu_table.txt", sep = "\t",
            quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(tax_table(physeq_bact_new), "filtered_data/16s_taxonomy.txt", sep="\t",
            quote=FALSE, row.names = TRUE, col.names = TRUE)
write.dna(refseq(physeq_bact_new), format="fasta", file = "filtered_data/16s_OTUs.fasta", colsep="")


# Chloroplast and Mitocondria 
write.dna(physeq_Chloro_Mito@refseq, format="fasta", 
          file = "PRJ_Greenhouse_Chloro_Mito_16S.fasta", colsep="")
write.table(physeq_Chloro_Mito@tax_table, 
            "PRJ_Greenhouse_taxonomy_Chloro_Mito_16S.txt", sep="\t", 
            quote=FALSE, row.names = TRUE, col.names = TRUE)

write.table(sort(taxa_sums(physeq_Chloro_Mito), decreasing = TRUE), 
            "PRJ_Greenhouse_abundance_Chloro_Mito_16S.txt", sep="\t", 
            quote=FALSE, row.names = TRUE, col.names = TRUE)

# write mock 
write.dna(refseq(physeq_bact_gs_clean_mock), format="fasta", file = "Mock_OTUs_16S.fasta", colsep="")
write.dna(refseq(physeq_fungi_gs_clean_mock), format="fasta", file = "mock_OTUs_ITS.fasta", colsep="")


# ******************************************************************************************----------------
# ALPHA DIVERSITY ------------------------------------------------------------------------------------------

# Remove Soil Samples Fungi
physeq_fungi_new <-
  subset_samples(physeq_fungi_new, Ecotype%in%c("Lowland", "Upland"))
otu_table(physeq_fungi_new) <-
  otu_table(physeq_fungi_new)[which(rowSums(otu_table(physeq_fungi_new)) > 0),]
physeq_fungi_new

count(as.data.frame(as.matrix(sample_data(physeq_fungi_new))), 
      Genotype, Soil_location)

# Remove Soil Samples bact
physeq_bact_new <-
  subset_samples(physeq_bact_new, Ecotype%in%c("Lowland", "Upland"))
otu_table(physeq_bact_new) <-
  otu_table(physeq_bact_new)[which(rowSums(otu_table(physeq_bact_new)) > 0),]
physeq_bact_new

count(as.data.frame(as.matrix(sample_data(physeq_bact_new))), 
      Genotype, Soil_location)


# RICHNESS -------------------------------------------------------------------------------------------------
MakeDf <- function(physeq, var){
  otu <- as(physeq@otu_table, "matrix")
  otu <- t(as.data.frame(otu))
  metadata <- as(physeq@sam_data, "data.frame")
  print("Are metadata sample identical to otu_table?")
  identical(rownames(otu), rownames(metadata)) %T>% print()
  df <- data.frame(richness = specnumber(otu, MARGIN = 1))
  df$rarefied <- rarefy(otu, sample=var, se = FALSE, MARGIN = 1)
  df$shannon <- diversity(otu, index = "shannon", MARGIN=1)#/log(df$richness)
  df$simpson <- diversity(otu, index = "simpson", MARGIN=1)
  df$readNo <- apply(otu,1,sum)
  df$Soil <- metadata$Soil_location
  df$Genotype <- metadata$Genotype
  df$Ecotype <- metadata$Ecotype
  df$Pot <- metadata$Pot
  print(table(df$Genotype, df$Soil))
  return(df)
}

alpha_df_fungi <-
  MakeDf(physeq_fungi_new, 10000) 
dim(alpha_df_fungi)

alpha_df_bact <-
  MakeDf(physeq_bact_new, 5000) 
dim(alpha_df_bact)

identical(alpha_df_fungi$Pot, alpha_df_bact$Pot)
sample_order <- match(alpha_df_bact$Pot, alpha_df_fungi$Pot)
alpha_df_fungi <- alpha_df_fungi[sample_order, ]
identical(alpha_df_fungi$Pot, alpha_df_bact$Pot)

# Removing samples that have same value for richness and rarefied, see ? rarefy
alpha_df_fungi_rare <-
  alpha_df_fungi[!alpha_df_fungi$richness == alpha_df_fungi$rarefied, ]
dim(alpha_df_fungi_rare)

alpha_df_bact_rare <-
  alpha_df_bact[!(as.numeric(alpha_df_bact$richness) == alpha_df_bact$rarefied), ]
dim(alpha_df_bact_rare)


# Visualize the data first ---------------------------------------------------------------------------------
ReadEffect <- function(df, index, var){
  read_plot <-
    df %>%
    ggplot(aes(x = readNo, 
               y = log(get(index)), 
               col = get(var))) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    theme_classic() +
    xlim(0, NA) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.x = element_text(angle =45, size = 8, hjust = 1, vjust = 1), 
          axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8)) 
  return(read_plot)
}


ReadEffect(alpha_df_fungi, "richness", "Soil") +
  labs(title = "Read No vs. Alpha index", x="Read No", y ="log(richness)") +
  guides(color = guide_legend(ncol = 1, title = "Soil"))


# Plotting Richness ----------------------------------------------------------------------------------------
PlotAlpha <- function(df){
  plot_alpha <- 
    ggplot(df, aes(x=Genotype, y=Var, color=Soil)) +
    geom_boxplot(outlier.shape = 1, outlier.size = 1, outlier.stroke = 1,
                 position = position_dodge(preserve = "single"), alpha=0.6, lwd = 0.5) +
    theme_classic() +
    #ylim(NA, 7.8) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.x = element_text(angle =45, size = 8, hjust = 1, vjust = 1), 
          axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8))
  return(plot_alpha)
}


alpha_df_fungi_filt %>%
  dplyr::select(Genotype, Soil, richness, shannon) %>%
  mutate(Var = log(richness)) %>% 
  PlotAlpha()

PlotInt <- function(df){
  plot_alpha <- 
    ggplot(df, aes(x=Genotype, y=Var, color=Soil, group=Soil)) +
    #geom_point() +
    geom_line() +
    theme_classic() +
    # ylim(NA, 7.5) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.x = element_text(angle =45, size = 8, hjust = 1, vjust = 1), 
          axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8))
  return(plot_alpha)
}

alpha_df_fungi %>%
  group_by(Genotype, Soil) %>%
  summarise(Var = mean(log(richness))) %>%
PlotInt() + labs(title = "Interaction plot", x="Genotype", y = "Mean of log(richness)") 

# Model diagnostic plots -----------------------------------------------------------------------------------
DiagPlots <- function(df, fit, var){
  par(mfrow=c(2,3))
  par(mar = c(4,4,1.6,1))
  plot(fit, which=1)
  plot(fit, which=2)
  plot(fit, which=3)
  plot(fit, which=5)
  hist(residuals(fit), breaks = 15,
       main = "Residuals distribution", 
       xlab = "Residuals",
       cex.main=1.5, font.main = 1)
  plot(fit$fitted.values, var, 
       xlab="predicted log(richness)", 
       ylab="actual log(richness)", 
       main = "Predicted vs Actual")
  abline(a=0, b=1, col="red", lwd=1, lty=1)
  return()
}

var = log(alpha_df_fungi$richness)
DiagPlots(alpha_df_fungi, fit_fungi_rich_m4, var)

autoplot(fit_fungi_rich_m4, which = 1:6, label.size = 3) +
  theme(axis.text.x = element_text(angle =45, size = 8, hjust = 1, vjust = 1))


# VISUALIZE readNo vs Richness relationship ----------------------------------------------------------------
# Plotting
title1 = text_grob("Read No vs. Alpha diversity metric", size = 12, face = 2)

grid.arrange(
  ggarrange(
    ggarrange(ReadEffect(alpha_df_fungi, "richness", "Soil") +
                labs(x="Read No", y ="log(richness)") +
                guides(color = guide_legend(nrow = 1, title = "Soil")) + ylim(NA, 6),
              ReadEffect(alpha_df_fungi_rare, "rarefied", "Soil") +
                labs(x="Read No", y ="log(rarefied)") +
                guides(color = guide_legend(nrow = 1, title = "Soil"))+ ylim(NA, 6),
              labels = c("A", "B"),
              common.legend = TRUE, 
              legend = "bottom"),
    ggarrange(ReadEffect(alpha_df_fungi, "richness", "Genotype") +
                labs(x="Read No", y ="log(richness)") +
                guides(color = guide_legend(nrow = 1, title = "Genotype"))+ ylim(NA, 6) +
              scale_shape_manual(values = c(0,1,2,5,3,8), 
                                 labels = c("Alamo", "Blackwell", "Cave-in-rock", "Kanlow", "Shelter","Southlow")),
              ReadEffect(alpha_df_fungi_rare, "rarefied", "Genotype") +
                labs(x="Read No", y ="log(rarefied)") +
                guides(color = guide_legend(nrow = 1, title = "Genotype"))+ ylim(NA, 6) +
                scale_shape_manual(values = c(0,1,2,5,3,8),
                                   labels = c("Alamo", "Blackwell", "Cave-in-rock", "Kanlow", "Shelter","Southlow")),
              labels = c("C","D"),
              common.legend = TRUE, 
              legend = "bottom"),
    ncol = 1, nrow = 2),
  top = title1)


grid.arrange(
  ggarrange(
    ggarrange(ReadEffect(alpha_df_bact[-188,], "richness", "Soil") + #remove 1 outlier for
                labs(x="Read No", y ="log(richness)") +
                guides(color = guide_legend(nrow = 1, title = "Soil"))+ ylim(NA, 7.7),
              ReadEffect(alpha_df_bact_rare, "rarefied", "Soil") +
                labs(x="Read No", y ="log(rarefied)") +
                guides(color = guide_legend(nrow = 1, title = "Soil"))+ ylim(NA, 7.7),
              labels = c("A", "B"),
              common.legend = TRUE, 
              legend = "bottom"),
    ggarrange(ReadEffect(alpha_df_bact[-188,], "richness", "Genotype") +
                labs(x="Read No", y ="log(richness)") +
                guides(color = guide_legend(nrow = 1, title = "Genotype"))+ ylim(NA, 7.7) +
                scale_shape_manual(values = c(0,1,2,5,3,8),
                                   labels = c("Alamo", "Blackwell", "Cave-in-rock", "Kanlow", "Shelter","Southlow")),
              ReadEffect(alpha_df_bact_rare, "rarefied", "Genotype") +
                labs(x="Read No", y ="log(rarefied)") +
                guides(color = guide_legend(nrow = 1, title = "Genotype"))+ ylim(NA, 7.7) +
                scale_shape_manual(values = c(0,1,2,5,3,8),
                                   labels = c("Alamo", "Blackwell", "Cave-in-rock", "Kanlow", "Shelter","Southlow")),
              labels = c("C","D"),
              common.legend = TRUE, 
              legend = "bottom"),
    ncol = 1, nrow = 2),
  top = title1)

# *************************************************************************************************-------------
# ROBUSTNESS - Robust linear models rlm() MASS package
# removing autliers - robustfying data
library(party)
library(flexplot)

hist(log(alpha_df_fungi_filt$richness), breaks = 20)

flexplot(richness ~ 1, data = alpha_df_fungi_filt) # visualize richness distrib

rf_fit <- party::cforest(richness ~ readNo + Soil + Genotype, data = alpha_df_fungi_filt)
estimates(rf_fit)

flexplot(richness ~ 1, data = alpha_df_fungi_filt)

# GLM MODELS ---------------------------------------------------------------------------------------------------

# Filtering outliers defined as a value 1.5 times of interquartile range above
# upper quartile (Q3) or below lower quartile (Q1). Simple boxplot method.

par(mfrow=c(1,3))
car::Boxplot(log(alpha_df_fungi$richness), id.method="y")
car::Boxplot(alpha_df_fungi$shannon, id.method="y")
car::Boxplot(alpha_df_fungi$readNo, id.method="y")
dev.off()

alpha_df_fungi[c(134, 171,32, 98, 30, 53), ]
alpha_df_fungi_rare[rownames(alpha_df_fungi_rare)== "Amp1411",]
alpha_df_fungi_filt <- 
  subset(alpha_df_fungi, !Pot%in%c("184", "109", "89", "87", "303", "198"))

par(mfrow=c(1,3))
car::Boxplot(log(alpha_df_bact$richness), id.method="y")
car::Boxplot(alpha_df_bact$readNo, id.method="y")
car::Boxplot(alpha_df_bact$shannon, id.method="y")
dev.off()

alpha_df_bact[c(150, 185, 115, 188, 190, 188), ]
alpha_df_bact_filt <- 
  subset(alpha_df_bact, !Pot%in%c("326", "93", "284", "117", "256", "117", "128"))


# 1) FUNGI RICHNESS -------------------------------------------------------------------------------------
fit_fungi_rich_m1 = glm(log(richness) ~ readNo, family = "Gamma", data=alpha_df_fungi_filt)
fit_fungi_rich_m2 = glm(log(richness)~ readNo + Soil,  family = "Gamma", data=alpha_df_fungi_filt)
fit_fungi_rich_m3 = glm(log(richness) ~ readNo + Soil + Genotype,  family = "Gamma",data=alpha_df_fungi_filt)
fit_fungi_rich_m4 = glm(log(richness) ~ readNo +  Soil + Genotype + Genotype:Soil, family = "Gamma", data=alpha_df_fungi_filt)

compareGLM(fit_fungi_rich_m1, fit_fungi_rich_m2, fit_fungi_rich_m3, fit_fungi_rich_m4)
anova(fit_fungi_rich_m1, fit_fungi_rich_m2, fit_fungi_rich_m3, fit_fungi_rich_m4, test="Chisq")

Anova(fit_fungi_rich_m4)
Anova(fit_fungi_rich_m3)
Anova(fit_fungi_rich_m4, test.statistic=c("F"), type = 2, 
      contrasts = c("contr.sum","contr.poly"))

 Anova(fit_fungi_rich_m4, test.statistic=c("F"), type = 3, 
       contrasts = c("contr.sum","contr.poly"))

shapiro.test(fit_fungi_rich_m4$residuals)
leveneTest(log(richness) ~ Soil, data=alpha_df_fungi_filt)
leveneTest(log(richness) ~ Genotype, data=alpha_df_fungi_filt)

var = log(alpha_df_fungi_filt$richness)
DiagPlots(alpha_df_fungi_filt, fit_fungi_rich_m4, var)

# 2) FUNGI SHANNON --------------------------------------------------------------------------------------------
fit_fungi_shan_m1 = glm(shannon ~ readNo, family = "Gamma",data=alpha_df_fungi_filt)
fit_fungi_shan_m2 = glm(shannon ~ readNo + Genotype, family = "Gamma",data=alpha_df_fungi_filt)
fit_fungi_shan_m3 = glm(shannon ~ readNo + Genotype + Soil, family = "Gamma",data=alpha_df_fungi_filt)
fit_fungi_shan_m4 = glm(shannon ~ readNo + Genotype + Soil + Genotype:Soil, family = "Gamma", data=alpha_df_fungi_filt)

compareGLM(fit_fungi_shan_m1, fit_fungi_shan_m2, fit_fungi_shan_m3, fit_fungi_shan_m4)
anova(fit_fungi_shan_m1, fit_fungi_shan_m2, fit_fungi_shan_m3, fit_fungi_shan_m4, test="Chisq")

Anova(fit_fungi_shan_m3)
Anova(fit_fungi_shan_m3, test.statistic=c("F"), type = 2, 
      contrasts = c("contr.sum","contr.poly"))

var = alpha_df_fungi_filt$shannon
DiagPlots(alpha_df_fungi_filt, fit_fungi_shan_m3, var)

shapiro.test(fit_fungi_shan_m3$residuals)
leveneTest(shannon ~ Soil, data=alpha_df_fungi_filt)
leveneTest(shannon ~ Genotype, data=alpha_df_fungi_filt)


# 3) BACTERIA RICHNESS ----------------------------------------------------------------------------------------------
fit_bact_rich_m1 = glm(log(richness) ~ readNo, family = "Gamma",data=alpha_df_bact_filt)
fit_bact_rich_m2 = glm(log(richness)~ readNo + Genotype, family = "Gamma",data=alpha_df_bact_filt)
fit_bact_rich_m3 = glm(log(richness) ~ readNo + Genotype + Soil, family = "Gamma",data=alpha_df_bact_filt)
fit_bact_rich_m4 = glm(log(richness) ~ readNo + Genotype + Soil + Genotype:Soil, family = "Gamma", data=alpha_df_bact_filt)

compareGLM(fit_bact_rich_m1, fit_bact_rich_m2, fit_bact_rich_m3, fit_bact_rich_m4)
anova(fit_bact_rich_m1, fit_bact_rich_m2, fit_bact_rich_m3, fit_bact_rich_m4,  test="Chisq")

Anova(fit_bact_rich_m1)
Anova(fit_bact_rich_m3, test.statistic=c("F"), type = 2, 
      contrasts = c("contr.sum","contr.poly"))

shapiro.test(fit_bact_rich_m3$residuals)
leveneTest(log(richness) ~ Soil, data=alpha_df_bact_filt)
leveneTest(log(richness) ~ Genotype, data=alpha_df_bact_filt)

var = log(alpha_df_bact_filt$richness)
DiagPlots(alpha_df_bact_filt, fit_bact_rich_m3, var)


# 4) BACTERIA SHANNON -----------------------------------------------------------------------------------------
fit_bact_shan_m1 = glm(shannon ~ readNo,family = "Gamma",  data=alpha_df_bact_filt)
fit_bact_shan_m2 = glm(shannon ~ readNo + Genotype,family = "Gamma",  data=alpha_df_bact_filt)
fit_bact_shan_m3 = glm(shannon ~ readNo + Genotype + Soil, family = "Gamma", data=alpha_df_bact_filt)
fit_bact_shan_m4 = glm(shannon ~ readNo + Genotype + Soil + Genotype:Soil, family = "Gamma",data=alpha_df_bact_filt)

compareGLM(fit_bact_shan_m1, fit_bact_shan_m2, fit_bact_shan_m3, fit_bact_shan_m4)
anova(fit_bact_shan_m1, fit_bact_shan_m2, fit_bact_shan_m3, fit_bact_shan_m4, test="Chisq")

Anova(fit_bact_shan_m3)
Anova(fit_bact_shan_m3, test.statistic=c("F"), type = 2, 
      contrasts = c("contr.sum","contr.poly"))

var = alpha_df_bact_filt$shannon
DiagPlots(alpha_df_bact_filt, fit_bact_shan_m3, var)

shapiro.test(fit_bact_shan_m3$residuals)
leveneTest(shannon ~ Soil, data=alpha_df_bact_filt)
leveneTest(shannon ~ Genotype, data=alpha_df_bact_filt)


# *** TABLE 1 - glm models ------------------------------------------------------------------

Anova(fit_fungi_rich_m4, test.statistic=c("F"), type = 2, 
      contrasts = c("contr.sum","contr.poly"))

Anova(fit_fungi_shan_m3, test.statistic=c("F"), type = 2, 
      contrasts = c("contr.sum","contr.poly"))

Anova(fit_bact_rich_m3, test.statistic=c("F"), type = 2, 
      contrasts = c("contr.sum","contr.poly"))

Anova(fit_bact_shan_m3, test.statistic=c("F"), type = 2, 
      contrasts = c("contr.sum","contr.poly"))


Anova(fit_fungi_rich_m4)
Anova(fit_fungi_shan_m3)
Anova(fit_bact_rich_m3)
Anova(fit_bact_shan_m3)


# Plotting soil and genotypes  -----------------------------------------------------------------

# Post-hoc Tukey tests among the three experimental treatments with partial residuals,
# after accounting for differential sequencing

residuals(fit_bact_rich_m1, type = c("deviance"))

PlotAlphaDiv <- function(df, fit, Var){
  res <-
    as.data.frame(residuals(fit, type = c("deviance")))
  colnames(res) <- "residuals"
  left_join(tibble::rownames_to_column(df),
            tibble::rownames_to_column(res),
            by="rowname") -> df
  head(df) %>% print()
  df$Genotype <- gsub("-", "_", df$Genotype)
  print(anova(aov(df$residuals ~ df[,Var])))
  tukey <-
    TukeyHSD(aov(df$residuals ~ df[,Var]))
  print(multcompLetters(extract_p(tukey[[1]])))
  plot_div <-
    ggplot(df, aes(x=get(Var), y=residuals)) +
    geom_boxplot(outlier.shape = 1, outlier.size = 1, outlier.stroke = 1,
               position = position_dodge(preserve = "single"), alpha=0.6, lwd = 0.5) +
    stat_summary(fun=mean, geom="point", shape=18, size=1.9, color="red", fill="red") +
    theme_classic() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(angle = 0, size = 10, face = "bold"),
        axis.text.x = element_text(angle =45, size = 8, hjust = 1, vjust = 1), 
        axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
        legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
        legend.title = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 8)) +
    labs(x=NULL)
return(plot_div)
}


PlotAlphaDiv(alpha_df_fungi_filt, fit_fungi_rich_1, "Soil") +
  stat_summary(geom = 'text', label = c("c", "a", "b", "a"), 
               fun = max, aes(y = 0.04), size = 3.5, color = "black") +
  labs(title = "log(richness)", x=NULL)

PlotAlphaDiv(alpha_df_fungi_filt, fit_fungi_rich_m1, "Genotype") +
  ylim(NA, 0.04)



# **** FIGURE 1 - glm models ---------------------------------------------------------------------- 
title5 = text_grob("Log(richness)", size = 12, face = 2)
title6 = text_grob("Shannon index", size = 12, face = 2)


ggarrange(
grid.arrange(
ggarrange(
  PlotAlphaDiv(alpha_df_fungi_filt, fit_fungi_rich_m1, "Soil") +
    stat_summary(geom = 'text', label = c("c", "a", "b", "a") , 
                 fun = max, aes(y = 0.2), size = 3.5, color = "black") +
    labs(title = "Fungi"),
  PlotAlphaDiv(alpha_df_fungi_filt, fit_fungi_rich_m1, "Genotype") +
    ylim(NA, 0.2) +
    labs(title = "Fungi", y=NULL) +
    scale_x_discrete(labels=c(Alamo = "Alamo",
                              Blackwell = "Blackwell",
                              `Cave_in_rock` = "Cave-in-rock",
                              Kanlow ="Kanlow",
                              Shlelter="Shelter",
                              Southlow="Southlow")),
  ncol = 2,
  nrow = 1, 
  widths = c(1,1.2),
  align = "h", labels = "A"), 
top=title5),
grid.arrange(
  ggarrange(
    PlotAlphaDiv(alpha_df_fungi_filt, fit_fungi_shan_m1, "Soil") +
      stat_summary(geom = 'text', label = c("c", "a", "ab", "b"), 
                   fun = max, aes(y = 0.5), size = 3.5, color = "black") +
      labs(title = "Fungi"),
    PlotAlphaDiv(alpha_df_fungi_filt, fit_fungi_shan_m1, "Genotype") +
      ylim(NA, 0.5) +
      labs(title = "Fungi",y=NULL) +
      scale_x_discrete(labels=c(Alamo = "Alamo",
                                Blackwell = "Blackwell",
                                `Cave_in_rock` = "Cave-in-rock",
                                Kanlow ="Kanlow",
                                Shlelter="Shelter",
                                Southlow="Southlow")),
    ncol = 2,
    nrow = 1, 
    widths = c(1,1.2),
    align = "h", labels = "B"), 
  top=title6),
  ggarrange(
    PlotAlphaDiv(alpha_df_bact_filt, fit_bact_rich_m1, "Soil") +
      stat_summary(geom = 'text', label = c("c", "a", "b", "b"), 
                   fun = max, aes(y = 0.2), size = 3.5, color = "black") +
      labs(title = "Bacteria") +
      ylim(-0.1, 0.2),
    PlotAlphaDiv(alpha_df_bact_filt, fit_bact_rich_m1, "Genotype") +
      ylim(-0.1, 0.2) +
      labs(title = "Bacteria",y=NULL)+
      scale_x_discrete(labels=c(Alamo = "Alamo",
                                Blackwell = "Blackwell",
                                `Cave_in_rock` = "Cave-in-rock",
                                Kanlow ="Kanlow",
                                Shlelter="Shelter",
                                Southlow="Southlow")),
    ncol = 2,
    nrow = 1, 
    widths = c(1,1.2),
    align = "h", labels = "C"),
  ggarrange(
    PlotAlphaDiv(alpha_df_bact_filt, fit_bact_shan_m1, "Soil") +
      stat_summary(geom = 'text', label = c("a", "a", "b", "b"), 
                   fun = max, aes(y = 0.3), size = 3.5, color = "black") +
      labs(title = "Bacteria"),
    PlotAlphaDiv(alpha_df_bact_filt, fit_bact_shan_m1, "Genotype") +
      ylim(NA, 0.3) +
      labs(title = "Bacteria",y=NULL)+
      scale_x_discrete(labels=c(Alamo = "Alamo",
                                Blackwell = "Blackwell",
                                `Cave_in_rock` = "Cave-in-rock",
                                Kanlow ="Kanlow",
                                Shlelter="Shelter",
                                Southlow="Southlow")),
    ncol = 2,
    nrow = 1, 
    widths = c(1,1.2),
    align = "h", labels = "D"),
ncol = 2,
nrow = 2)


# *** FIGURE 1 - alpha diversity ---------------------------------------------------------------------------

title2 = text_grob("Alpha diversity", size = 12, face = 2)

grid.arrange(
  ggarrange(alpha_df_fungi_filt %>%
              dplyr::select(Genotype, Soil, richness, shannon) %>%
              mutate(Var = log(richness)) %>% 
              PlotAlpha() + 
              scale_color_manual(values = Pal_soil) +
              labs(title = "Fungi", x=NULL, y = "log(richness)"),
            alpha_df_bact_filt %>%
              dplyr::select(Genotype, Soil, richness, shannon) %>%
              mutate(Var = log(richness)) %>% 
              PlotAlpha() +
              scale_color_manual(values = Pal_soil) +
              labs(title = "Bacteria", x=NULL, y = NULL),
            alpha_df_fungi_filt %>%
              dplyr::select(Genotype, Soil, richness, shannon) %>%
              mutate(Var = shannon) %>% 
              PlotAlpha() +
              scale_color_manual(values = Pal_soil) +
              labs(title = NULL, x="Soil", y = "Shannon index"),
            alpha_df_bact_filt %>%
              dplyr::select(Genotype, Soil, richness, shannon) %>%
              mutate(Var = shannon) %>% 
              PlotAlpha() +
              scale_color_manual(values = Pal_soil) +
              labs(title = NULL, x="Soil", y = NULL),
            labels = c("A","","B",""),
            common.legend = TRUE, 
            legend = "bottom",
            ncol = 2,
            nrow = 2),
  top = title2) 


# **********************************************************************------------------------------------
# BIOMASS MODELS -------------------------------------------------------------------------------------------
# *** Fungi *** --------------------------------------------------------------------------------------------
head(physeq_fungi_new@sam_data)

otu_fungi_new <- as(physeq_fungi_new@otu_table, "matrix")
otu_fungi_new <- t(as.data.frame(otu_fungi_new))
str(otu_fungi_new)

taxa_fungi_new <- as(physeq_fungi_new@tax_table, "matrix")
taxa_fungi_new <- as.data.frame(taxa_fungi_new)
head(taxa_fungi_new)

meta_fungi_new <- 
  as(physeq_fungi_new@sam_data, "data.frame") %>%
  dplyr::select(Pot, N_Collared_Leaves_12wk, Total_leaves_12wks, Dead_leaves_12wks,
       Yellow_leaves_12wks, N_tillers_12wks, Flower_12wk, Spikelets_Emerged_percent_12wks,
       Anters_Emerged_percent_12wks, Stress_12wks, plant_height_12wks_cm, Stage_12wks,
       n_tillers_16wks, Plant_height_16wks_cm, Root_lenght_16wks_cm,
       developmental.estage_16wkE, developmental.stage_16wks_detail, Development_16wks_ERS, 
       Flower_16wk, nodes_LeadingTiller_16wks, aerial_part_dry_weight_16wks_grams, 
       AverageDryWeightPerGroup, PotMinusAGroupDryWeightAverage, 
       pH, P, K, Ca, Mg, OM, NO3)
str(meta_fungi_new)

# Calculating PCoA in vegan --------------------------------------------------------------------------------
cmdscale(vegdist(otu_fungi_new, method = "bray"), eig=TRUE) -> pcoa_its

as.data.frame(pcoa_its$points) -> pcoa_fungi
colnames(pcoa_fungi) <- c("PCoA.1", "PCoA.2")
identical(rownames(pcoa_fungi), rownames(alpha_df_fungi_filt))

alpha_df_fungi_filt <-
  left_join(tibble::rownames_to_column(alpha_df_fungi_filt), #keep this samples only
            tibble::rownames_to_column(pcoa_fungi), by="rowname")
str(alpha_df_fungi_filt)

# nmds_fungi <- metaMDS(otu_fungi_new, k=2, trymax=200, autotransform = TRUE)
# stressplot(nmds_fungi)
# df_nmds_fungi <-  as.data.frame(nmds_fungi$points)
# head(df_nmds_fungi)

# Calculating dispersion -----------------------------------------------------------------------------------
head(otu_fungi_new)
otu_fungi_new_filt <- 
  subset(otu_fungi_new, rownames(otu_fungi_new) %in% alpha_df_fungi_filt$rowname) %>%
  as.data.frame() 

otu_fungi_new_filt <- 
  otu_fungi_new_filt[, colSums(otu_fungi_new_filt)>0]

str(otu_fungi_new_filt)
str(alpha_df_fungi_filt)

# reorder
rownames(alpha_df_fungi_filt) <- alpha_df_fungi_filt$rowname
identical(rownames(alpha_df_fungi_filt), rownames(otu_fungi_new_filt))
order_fungi <- match(rownames(otu_fungi_new_filt), rownames(alpha_df_fungi_filt))
alpha_df_fungi_filt <- alpha_df_fungi_filt[order_fungi,]

permdisp_fungi_soil <- 
  betadisper(vegdist(otu_fungi_new_filt, method = "bray"), alpha_df_fungi_filt$Soil) 
dist_fungi_soil<-
  data.frame(permdisp_fungi_soil$group, permdisp_fungi_soil$distances)
colnames(dist_fungi_soil) <- c("value", "dispSoil")
dist_fungi_soil

permdisp_fungi_genotype <- 
  betadisper(vegdist(otu_fungi_new_filt, method = "bray"), alpha_df_fungi_filt$Genotype) 
dist_fungi_genotype<-
  data.frame(permdisp_fungi_genotype$group, permdisp_fungi_genotype$distances)
colnames(dist_fungi_genotype) <- c("value", "dispGenotype")
dist_fungi_genotype

identical(rownames(dist_fungi_genotype), rownames(dist_fungi_soil))

meta_fungi_merged <-
cbind(dist_fungi_genotype, dist_fungi_soil) %>%
  dplyr::select(dispGenotype, dispSoil) %>%
  tibble::rownames_to_column() %>%
  left_join(alpha_df_fungi_filt, by="rowname") %>%
  left_join(tibble::rownames_to_column(meta_fungi_new[,-1]), by="rowname")
str(meta_fungi_merged)


# transform to z-score -------------------------------------------------------------------------------------
#apply(mod_fungi_3[,2:8], 2, function(x) scale(x, center = TRUE, scale = TRUE))
VarStand <- function(df) {
  for (n in names(df)) { 
    if (class(df[[n]]) == "numeric" | class(df[[n]]) == "integer"){
      var = paste(n,"_z", sep="")
      df[[var]] <- scale(df[[n]], center = TRUE, scale = TRUE)
      df[[n]] = NULL
    }
  }
  return(df)
}

mod_fungi_1 <- meta_fungi_merged[, c("Pot",
                                     "richness", 
                                     "shannon", 
                                     "PCoA.1",
                                     "PCoA.2",
                                     "dispSoil",
                                     "dispGenotype", 
                                     "readNo")]

VarStand(mod_fungi_1)

# mod_fungi_2 <- meta_fungi_merged[, c("Pot",
#                                      "aerial_part_dry_weight_16wks_grams",
#                                      "N_Collared_Leaves_12wk", 
#                                      "Total_leaves_12wks", 
#                                      "Dead_leaves_12wks",
#                                      "Yellow_leaves_12wks", 
#                                      "N_tillers_12wks", 
#                                      "plant_height_12wks_cm",
#                                      "n_tillers_16wks", 
#                                      "Plant_height_16wks_cm", 
#                                      "Root_lenght_16wks_cm",
#                                      "nodes_LeadingTiller_16wks")]

mod_fungi_2 <- meta_fungi_merged[, c("Pot",
                                     "aerial_part_dry_weight_16wks_grams")]

VarStand(mod_fungi_2)

mod_fungi_3 <- meta_fungi_merged[, c("Pot",
                                     "pH", 
                                     "P", 
                                     "K",
                                     "Ca", 
                                     "Mg", 
                                     "OM", 
                                     "NO3")]

VarStand(mod_fungi_3)
mod_fungi_3 %>% mutate_if(is.numeric, as.factor)

mod_fungi_4 <- meta_fungi_merged[, c("Pot",
                                     #"Flower_12wk", 
                                     #"Stress_12wks",
                                     #"Stage_12wks",
                                     #"developmental.estage_16wkE", 
                                     #"developmental.stage_16wks_detail", 
                                     #"Development_16wks_ERS", 
                                     #"Flower_16wk",
                                     "Soil",
                                     "Genotype",
                                     "Ecotype")]


# using chemistry variables as numeric ---------------------------------------------------------------------
# physeq_fungi_new@sam_data$SoilGen <-
#   paste(physeq_fungi_new@sam_data$Soil, physeq_fungi_new@sam_data$Genotype)

meta_fungi_filt_2 <-
  left_join(mod_fungi_2, 
            VarStand(mod_fungi_1), by="Pot") %>%
  left_join(VarStand(mod_fungi_3), by="Pot") %>%
  left_join(mod_fungi_4, by="Pot")

# remove Pot and remove 1 sample with too many NA
meta_fungi_filt_2 <-
  meta_fungi_filt_2[,-1][complete.cases(meta_fungi_filt_2[,-1])==TRUE, ]
str(meta_fungi_filt_2)

#recoding varaiables 
# meta_fungi_filt_2$Stress_12wks <-
#   as.factor(ifelse(meta_fungi_filt_2$Stress_12wks=="", "N", paste(meta_fungi_filt_2$Stress_12wks)))

# rename colnames
colnames(meta_fungi_filt_2) <-c("Biomass",
                                "Richness","Shannon",
                                "PCoA.1","PCoA.2","Disp.Soil","Disp.Genotype",
                                "Read.No",
                                "pH","P","K","Ca","Mg","OM","NO3",
                                #"Stress.12wks",
                                "Soil","Genotype","Ecotype")
head(meta_fungi_filt_2)

# loop to store 100 boruta runs
sign_var <- vector(mode = "character")

for(i in 1:99) {
  sel_attr[i]  <- Boruta(Biomass ~., meta_fungi_filt_2, pValue = 0.05, 
                         mcAdj = TRUE, maxRuns=100, doTrace = 3)
  sign_var <- append(sign_var, getSelectedAttributes(sel_attr, withTentative = TRUE))
}

sign_var
unique(sign_var)

df_fungi_RF_2 <-
  meta_fungi_filt_2[,c(unique(sign_var), "Biomass")]

# try tuning the model first
round(sqrt(ncol(df_fungi_RF_2[, 1:(ncol(df_fungi_RF_2) - 1)])))

set.seed(12345)

bestmtry_fungi_Bmass_2 <-
  tuneRF(
    x = df_fungi_RF_2[, 1:(ncol(df_fungi_RF_2) - 1)],
    y = df_fungi_RF_2$Biomass,
    mtryStart = 4,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )


RF_fungi_Bmass_2 <-
  randomForest(
    x = df_fungi_RF_2[, 1:(ncol(df_fungi_RF_2) - 1)],
    y = df_fungi_RF_2$Biomass,
    ntree = 1001,
    mtry = 2,
    importance = TRUE,
    proximity = TRUE
  )

RF_fungi_Bmass_2
plot(RF_fungi_Bmass_2)

# Assessing model significance using permitations
set.seed(110324)
perm_RF_fungi_Bmass_2 <-
  rf.significance(
    x = RF_fungi_Bmass_2,
    xdata = df_fungi_RF_2[, 1:(ncol(df_fungi_RF_2) - 1)],
    nperm = 999,
    nmtry = 4,
    ntree = 1001
  )

perm_RF_fungi_Bmass_2 # model significant = 0.001 


# 1 - Plotting error ---------------------------------------------------------------------------------------
PlotError <- function(rf_model){
  model_df <- (data.frame(Trees = 1:1001, Error = rf_model$mse))
  ggplot(data=model_df, aes(x=Trees, y=Error)) +
    labs(title = "Model Errors", y="Error", x="Tree") +
    theme_classic() +
    geom_line(color="red", size=0.8) +
    ylim(0, NA) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.x = element_text(angle =0, size = 7, hjust = 0.5, vjust = 1), 
          axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8)) -> error_plot
  return(error_plot)
}


PlotError(RF_fungi_Bmass) +
  annotate("text", x=Inf, y = Inf,
           label=paste("Mean squared error:", round(last(RF_fungi_Bmass$mse), 2)), size=2.5, vjust=1, hjust=1) +
  annotate("text", x=Inf, y = Inf, 
           label= paste("% Var explained:", round(last(RF_fungi_Bmass$rsq*100),2)), size=2.5, vjust=3, hjust=1) +
  annotate("text", x=Inf, y = Inf, 
           label = paste("italic(p) ==", round(perm_RF_fungi_Bmass$pValue, 4)), parse = TRUE, size=2.5, vjust=4, hjust=1)


# 2 - Plotting Line ----------------------------------------------------------------------------------------
PlotLine <- function(rf_model, metadata){
  df_model <- data.frame(actual=rf_model$y, pred=rf_model$predicted)
  if (identical(rownames(metadata), rownames(df_model))==TRUE){
    df_model$Soil <- metadata$Soil
    df_model$Genotype <- metadata$Genotype
    df_model %T>% print()
    ggplot(data=df_model, aes(x=actual, y=pred)) +
      geom_point(aes(shape=Genotype, color=Soil), size=1.5, stroke=0.5) +
      geom_smooth(method = "lm", formula = "y ~ x", se = TRUE, color="black", size=0.5) +
      #geom_smooth(method = "loess", formula = "y ~ x", se = TRUE, color="red", size=0.8) +
      theme_classic() +
      scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
      # # scale_shape_manual(values = c(1, 2, 3),
      #                    labels = c('0 dpf', '13 dpf', '33 dpf')) +
      # guides(color = guide_legend(nrow = 2),
      #        shape = guide_legend(nrow = 3)) +
      theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
            axis.title = element_text(angle = 0, size = 10, face = "bold"),
            axis.text.x = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 1), 
            axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
            legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
            legend.title = element_blank(), legend.background = element_blank(), 
            legend.text = element_text(size = 8)) -> line_plot
    return(line_plot)
  }else{
    stop("Error: dataframe and metadata rownames are not matching!")
  }
}

# 3  Inf -Inf  Bottom Right h1,v0        1        0

PlotLine(RF_fungi_Bmass, df_fungi_RF) + theme(legend.position = c(0.1, 0.8)) +
  annotate("text", x=Inf, y = -Inf,
         label=paste("Mean squared error:", round(last(RF_fungi_Bmass$mse), 2)), size=2.5, vjust=-6, hjust=1) +
  annotate("text", x=Inf, y = -Inf, 
           label= paste("% Var explained:", round(last(RF_fungi_Bmass$rsq*100),2)), size=2.5, vjust=-4, hjust=1) +
  annotate("text", x=Inf, y = -Inf, 
           label = paste("italic(p) ==", round(perm_RF_fungi_Bmass$pValue, 4)), parse = TRUE, size=2.5, vjust=-1.5, hjust=1)



# 3- Plotting features -------------------------------------------------------------------------------------
PlotFeature <- function(rf_model, var){
  imp_RF <- as.data.frame(rf_model$importance)
  imp_RF$features <- rownames(imp_RF)
  colnames(imp_RF) <- c("IncMSE","IncNodePurity","features")
  imp_RF <- arrange(imp_RF, desc(imp_RF[,var]))
  imp_RF %T>% print()
  ggplot(data=imp_RF) + 
    geom_bar(aes(x= reorder(features, -get(var)),
                 y= get(var)), color="grey80", fill="grey80",stat="identity") +
    coord_flip() +
    #ylim(0, 0.03) +
    theme_classic() +
    theme(plot.margin=unit(c(7,9,7,7),"pt")) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.x = element_text(angle =45, size = 7, hjust = 1, vjust = 1), 
          axis.text.y = element_text(angle = 0, size = 7, hjust = 1, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8)) -> plot_importance
  return(plot_importance)
}


PlotFeature(RF_fungi_Bmass, "IncMSE" ) + labs(title = "Predictors") 
PlotFeature(RF_fungi_Bmass, "IncNodePurity" ) + labs(title = "Predictors") 
as.data.frame(RF_fungi_Bmass$importance)


# 4 PCA of Predictors --------------------------------------------------------------------------------------
# Test for colinearity using PCA

# PlotNeutralPCA <- function(df, Var){
#   pca_plot <-
#     autoplot(
#       prcomp(x = df, scale= TRUE, center=TRUE), # calculates principal compenents and pltos with ggplot2
#       data = Var, label = FALSE, shape = "Genotype", colour="Soil", # add metadate, labels of objects
#       loadings = TRUE, loadings.colour = "black", size=1.5,
#       frame = FALSE, frame.colour = "Soil", loadings.label.colour = "black",
#       loadings.label = TRUE, loadings.label.size = 3, loadings.label.repel = TRUE) +
#     labs(title = "PCA") +
#     # scale_colour_manual(values = paletteCB4) +
#     # scale_fill_manual(values = paletteCB4) +
#     # scale_shape_manual(values = c(21,22,24)) +
#     theme_classic() +
#     theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
#           plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
#           axis.title = element_text(angle = 0, size = 10, face = "bold"),
#           axis.text.x = element_text(angle =0, size = 7, hjust = 0.5, vjust = 1), 
#           axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
#           legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
#           legend.title = element_text(size = 10, face = "bold"), 
#           legend.text = element_text(size = 8))
#   # guides(color = guide_legend(ncol = 2), #title.position="top"
#   #         fill = guide_legend(ncol= 2),
#   #         shape = guide_legend(ncol = 1)) +
#   # theme(legend.margin = margin(0,-0.5,0,0, unit="cm")) # reduce space betweem legends
#   return(pca_plot)
# }
# 
# 
# PlotNeutralPCA(meta_fungi_filt_2[,c(1:6, 8:15)], meta_fungi_filt_2)
# 
# prcomp_fungi <- 
#   prcomp(meta_fungi_filt_2[c(1:6, 8:15)], center = TRUE, scale = TRUE)
# 
# autoplot(prcomp_fungi, loadings = TRUE, loadings.label = TRUE)
# pca_points_fungi <- 
#   as_tibble(prcomp_fungi$x) %>% 
#   bind_cols(meta_fungi_filt_2)
# pca_points_fungi
# summary(lm(PC1 ~ NO3, pca_points_fungi))
# 
# prcomp_fungi
# 
# PlotNeutralPCA <- function(df){
#   pca_plot <-
#     ggplot(df, aes(x = PC1, y = PC2)) +
#     geom_point(aes(colour = Soil, shape= Genotype)) +
#     # autoplot(
#     #   prcomp(x = df, scale= TRUE, center=TRUE), # calculates principal compenents and pltos with ggplot2
#     #   data = Var, label = FALSE, shape = "Genotype", colour="Soil", # add metadate, labels of objects
#     #   loadings = TRUE, loadings.colour = "black", size=1.5,
#     #   frame = FALSE, frame.colour = "Soil", loadings.label.colour = "black",
#     #   loadings.label = TRUE, loadings.label.size = 3, loadings.label.repel = TRUE) +
#     labs(title = "PCA") +
#     # scale_colour_manual(values = paletteCB4) +
#     # scale_fill_manual(values = paletteCB4) +
#     # scale_shape_manual(values = c(21,22,24)) +
#     theme_classic() +
#     theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
#           plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
#           axis.title = element_text(angle = 0, size = 10, face = "bold"),
#           axis.text.x = element_text(angle =0, size = 7, hjust = 0.5, vjust = 1), 
#           axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
#           legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
#           legend.title = element_text(size = 10, face = "bold"), 
#           legend.text = element_text(size = 8))
#   # guides(color = guide_legend(ncol = 2), #title.position="top"
#   #         fill = guide_legend(ncol= 2),
#   #         shape = guide_legend(ncol = 1)) +
#   # theme(legend.margin = margin(0,-0.5,0,0, unit="cm")) # reduce space betweem legends
#   return(pca_plot)
# }
# 
# 
# pca_load_fungi <-
#   as_tibble(prcomp_fungi$rotation, rownames = 'variable') %>%
#   mutate(variable = dplyr::recode(variable,
#                                   "Biomass" = "Biomass",
#                                   "Richness" = "Richness",
#                                   "Shannon" = "Shannon",
#                                   "PCoA.1" = "PCoA.1",
#                                   "PCoA.2" = "PCoA.2",
#                                   "Disp.Soil" = "Disp.Soil",
#                                   "Read.No" = "Read.No",
#                                   "pH"="pH",
#                                   "P" = expression(PO[4]^{"-"}),
#                                   "K" = expression(K^{"+"}),
#                                   "Ca" =expression(Ca^{"2+"}),
#                                   "Mg" =expression(Mg^{"2+"}),
#                                   "OM" = "OM",
#                                   "NO3" = expression(NO[3]^{"-"})))
# 
# 
# 
# PlotNeutralPCA(pca_points_fungi) +
#   geom_segment(data = pca_load_fungi, 
#                aes(x = 0, y = 0, 
#                    xend = PC1*6,
#                    yend = PC2*6),
#                arrow = arrow(length = unit(1/2, "picas"))) +
#   annotate("text", x = (pca_load_fungi$PC1*6), y = (pca_load_fungi$PC2*5.2),
#            label = c(
#              "Biomass",
#              "Richness",
#              "Shannon",
#              "PCoA.1",
#              "PCoA.2",
#              "Disp.Soil",
#              "Read.No",
#              "pH",
#              expression(PO[4]^{"-"}),
#              expression(K^{"+"}),
#              expression(Ca^{"2+"}),
#              expression(Mg^{"2+"}),
#              "OM"="OM",
#              expression(NO[3]^{"-"})),
#            size = 3.5) 

# New functions 

PlotNeutralPCAfungi <- function(df){
  #generate pca
  pca_res <-
    df %>%
    select(-Soil, -Genotype, -Ecotype, -Disp.Genotype) %>%
    prcomp(x = ., center = TRUE, scale = TRUE)
  # extract points for the samples
  pca_points <- 
    as_tibble(pca_res$x) %>% 
    bind_cols(df)
  # extract point for the loadings
  pca_load <-
    as_tibble(pca_res$rotation, rownames = 'variable')
  # % variation for each axis
  axis_var <- 
    round(as.vector(summary(pca_res)$importance[2,])*100,1)
  
  pca_plot <-
    ggplot(pca_points, aes(x = PC1, y = PC2)) +
    geom_point(aes(colour = Soil, shape= Genotype)) +
    labs(title = "PCA",
         x=as.expression(paste("PC1 (",axis_var[1],"%)"), parse=TRUE),
         y=as.expression(paste("PC2 (",axis_var[2],"%)"), parse=TRUE)) +
    theme_classic() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.x = element_text(angle =0, size = 7, hjust = 0.5, vjust = 1), 
          axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8))
    pca_plot <-
      pca_plot + 
      #grids(linetype = "dashed") +
      geom_segment(data = pca_load, 
                         aes(x = 0, y = 0, 
                             xend = PC1*7.7,
                             yend = PC2*7.7),
                         arrow = arrow(length = unit(1/2, "picas"))) +
            annotate("text", x = (pca_load$PC1*8.3), y = (pca_load$PC2*8),
                     size = 3,
                     label = c(
                       "Biomass",
                       "Richness",
                       "Shannon",
                       "PCoA.1",
                       "PCoA.2",
                       "Disp.Soil",
                       "Read.No",
                       "pH",
                       expression(PO[4]^{"3-"}),
                       expression(K^{"+"}),
                       expression(Ca^{"2+"}),
                       expression(Mg^{"2+"}),
                       "OM",
                       expression(NO[3]^{"-"})))
  return(pca_plot)
}

PlotNeutralPCAfungi(meta_fungi_filt_2)

PlotNeutralPCAbact <- function(df){
  pca_res <-
    df %>%
    select(-Soil, -Genotype, -Ecotype) %>%
    prcomp(x = ., center = TRUE, scale = TRUE)
  pca_points <- 
    as_tibble(pca_res$x) %>% 
    bind_cols(df)
  pca_load <-
    as_tibble(pca_res$rotation, rownames = 'variable')
  axis_var <- 
    round(as.vector(summary(pca_res)$importance[2,])*100,1)
  pca_plot <-
    ggplot(pca_points, aes(x = PC1, y = PC2)) +
    geom_point(aes(colour = Soil, shape= Genotype)) +
    labs(title = "PCA",
         x=as.expression(paste("PC1 (",axis_var[1],"%)"), parse=TRUE),
         y=as.expression(paste("PC2 (",axis_var[2],"%)"), parse=TRUE)) +
    theme_classic() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.x = element_text(angle =0, size = 7, hjust = 0.5, vjust = 1), 
          axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8))
  pca_plot <-
    pca_plot + 
    #grids(linetype = "dashed") +
    geom_segment(data = pca_load, 
                            aes(x = 0, y = 0, 
                                xend = PC1*7.7,
                                yend = PC2*7.7),
                            arrow = arrow(length = unit(1/2, "picas"))) +
    annotate("text", x = (pca_load$PC1*8.6), y = (pca_load$PC2*8),
             size = 3,
             label = c(
               "Biomass",
               "Richness",
               "Shannon",
               "PCoA.1",
               "PCoA.2",
               "Disp.Soil",
               "Disp.Genotype",
               "Read.No",
               "pH",
               expression(PO[4]^{"3-"}),
               expression(K^{"+"}),
               expression(Ca^{"2+"}),
               expression(Mg^{"2+"}),
               "OM",
               expression(NO[3]^{"-"})))
  return(pca_plot)
}

PlotNeutralPCAbact(meta_bact_filt_2)


# Multicollinearity does not affect the accuracy of predictive models, including regression models. 
# Take the attached image as an example. The features in the x and y axis are clearly correlated; 
# however, you need both of them to create an accurate classifier. If you discard one of them for 
# being highly correlated with the other one, the performance of your model will decrease.
#If you want to remove the collinearity, you can always use PCA to project the data into a new 
# space where the 'new features' will be orthogonal to each other. You can then, train your model
# with the new features, but you will find that the performance is the same. You simply rotated 
# your original decision boundary. Now, where multicollinearity becomes 'an issue' is when you want 
# to 'interpret' the parameters learned by your model. In other words, you cannot say that the 
# feature with the 'biggest weight' is 'the most important' when the features are correlated.
# Note that this is independent on the accuracy of the model, this is only the interpretation part, 
# which in my opinion, you should not be doing anyway. 

# Using chemistry data as factor ---------------------------------------------------------------------------
meta_fungi_filt <-
  left_join(mod_fungi_2, 
            VarStand(mod_fungi_1), by="Pot") %>%
  left_join(mutate_if(mod_fungi_3, is.numeric, as.factor), by="Pot") %>%
  left_join(mod_fungi_4, by="Pot")

# remove Pot and remove 1 sample with too many NA
meta_fungi_filt <-
  meta_fungi_filt[,-1][complete.cases(meta_fungi_filt[,-1])==TRUE, ]
str(meta_fungi_filt)

#recoding varaiables 
meta_fungi_filt$Stress_12wks <-
  as.factor(ifelse(meta_fungi_filt$Stress_12wks=="", "N", paste(meta_fungi_filt$Stress_12wks)))

# loop to store 100 boruta runs
sign_var <- vector(mode = "character")

for(i in 1:99) {
  sel_attr[i]  <- Boruta(aerial_part_dry_weight_16wks_grams ~., meta_fungi_filt, pValue = 0.05, 
                         mcAdj = TRUE, maxRuns=100, doTrace = 3)
  sign_var <- append(sign_var, getSelectedAttributes(sel_attr, withTentative = TRUE))
}

sign_var
unique(sign_var)

df_fungi_RF <-
  meta_fungi_filt[,c(unique(sign_var), "aerial_part_dry_weight_16wks_grams")]

# try tuning the model first
round(sqrt(ncol(df_fungi_RF[, 1:(ncol(df_fungi_RF) - 1)])))

set.seed(12345)

bestmtry_fungi_Bmass <-
  tuneRF(
    x = df_fungi_RF[, 1:(ncol(df_fungi_RF) - 1)],
    y = df_fungi_RF$aerial_part_dry_weight_16wks_grams,
    mtryStart = 4,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )

RF_fungi_Bmass <-
  randomForest(
    x = df_fungi_RF[, 1:(ncol(df_fungi_RF) - 1)],
    y = df_fungi_RF$aerial_part_dry_weight_16wks_grams,
    ntree = 1001,
    mtry = 2,
    importance = TRUE,
    proximity = TRUE
  )

RF_fungi_Bmass
plot(RF_fungi_Bmass)

# Assessing model significance using permitations
set.seed(110324)
perm_RF_fungi_Bmass <-
  rf.significance(
    x = RF_fungi_Bmass,
    xdata = df_fungi_RF[, 1:(ncol(df_fungi_RF) - 1)],
    nperm = 999,
    nmtry = 1,
    ntree = 1001
  )

perm_RF_fungi_Bmass # model significant = 0.001 

# FIGURE SXX - RF MODELS with factor 
ggarrange(PlotLine(RF_fungi_Bmass, df_fungi_RF) + 
              theme(legend.position = c(0.15, 0.8)) +
              scale_color_manual(values = Pal_soil) +
              scale_shape_manual(values = c(15,16,17,18,3,8)) +
              labs(title = "Random Forest ", y="Predicted\nAerial Plant Biomass", x="Observed\nAerial Plant Biomass") +
              annotate("text", x=Inf, y = -Inf,
                       label=paste("Mean squared error:", round(last(RF_fungi_Bmass$mse), 2)), size=2.5, vjust=-6, hjust=1) +
              annotate("text", x=Inf, y = -Inf, 
                       label= paste("% Var explained:", round(last(RF_fungi_Bmass$rsq*100),2)), size=2.5, vjust=-4, hjust=1) +
              annotate("text", x=Inf, y = -Inf, 
                       label = paste("italic(p) ==", round(perm_RF_fungi_Bmass$pValue, 4)), parse = TRUE, size=2.5, vjust=-1.5, hjust=1),
            PlotFeature(RF_fungi_Bmass, "IncMSE") +
              labs(title = "Predictors", x=NULL, y="% Increase MSE "),
            PlotFeature(RF_fungi_Bmass,  "IncNodePurity") +
              labs(title = "Predictors", x=NULL, y="Increase Node Purity"),
            PlotError(RF_fungi_Bmass) +
              annotate("text", x=Inf, y = Inf, 
                       label= paste("% Explained var.:", round(last(RF_fungi_Bmass$rsq*100),2)), size=2.5, vjust=1, hjust=1) +
              annotate("text", x=Inf, y = Inf,
                       label=paste("MSE:", round(last(RF_fungi_Bmass$mse), 2)), size=2.5, vjust=3, hjust=1) +
              annotate("text", x=Inf, y = Inf, 
                       label = paste("italic(p) ==", round(perm_RF_fungi_Bmass$pValue, 4)), parse = TRUE, size=2.5, vjust=4, hjust=1),
          labels = c("A","B","C", "D"),
          widths = c(1.1, 0.6, 0.6, 0.5),
          align = "h",
          ncol = 4, 
          nrow =1)

# *** Bacteria *** -----------------------------------------------------------------------------------------
head(physeq_bact_new@sam_data)

otu_bact_new <- as(physeq_bact_new@otu_table, "matrix")
otu_bact_new <- t(as.data.frame(otu_bact_new))
str(otu_bact_new)

colSums(otu_bact_new)

taxa_bact_new <- as(physeq_bact_new@tax_table, "matrix")
taxa_bact_new <- as.data.frame(taxa_bact_new)
head(taxa_bact_new)

meta_bact_new <- 
  as(physeq_bact_new@sam_data, "data.frame") %>%
  dplyr::select(Pot, N_Collared_Leaves_12wk, Total_leaves_12wks, Dead_leaves_12wks,
                Yellow_leaves_12wks, N_tillers_12wks, Flower_12wk, Spikelets_Emerged_percent_12wks,
                Anters_Emerged_percent_12wks, Stress_12wks, plant_height_12wks_cm, Stage_12wks,
                n_tillers_16wks, Plant_height_16wks_cm, Root_lenght_16wks_cm,
                developmental.estage_16wkE, developmental.stage_16wks_detail, Development_16wks_ERS, 
                Flower_16wk, nodes_LeadingTiller_16wks, aerial_part_dry_weight_16wks_grams, 
                AverageDryWeightPerGroup, PotMinusAGroupDryWeightAverage, 
                pH, P, K, Ca, Mg, OM, NO3)
str(meta_bact_new)


# Calculating PCoA in vegan --------------------------------------------------------------------------------
cmdscale(vegdist(otu_bact_new, method = "bray"), eig=TRUE) -> pcoa_16s

as.data.frame(pcoa_16s$points) -> pcoa_bact
colnames(pcoa_bact) <- c("PCoA.1", "PCoA.2")
identical(rownames(pcoa_bact), rownames(alpha_df_bact_filt))

alpha_df_bact_filt <-
  left_join(tibble::rownames_to_column(alpha_df_bact_filt), #keep this samples only
            tibble::rownames_to_column(pcoa_bact), by="rowname")
str(alpha_df_bact_filt)

nmds_bact <- metaMDS(otu_bact_new, k=2, trymax=200, distance = "bray", autotransform = TRUE, weakties = TRUE)
stressplot(nmds_bact)
nmds_bact <- metaMDS(otu_bact_new, previous.best = nmds_bact)
df_nmds_bact <-  as.data.frame(nmds_bact$points)
head(df_nmds_bact)

# Calculating dispersion
head(otu_bact_new)
otu_bact_new_filt <- 
  subset(otu_bact_new, rownames(otu_bact_new) %in% alpha_df_bact_filt$rowname) %>%
  as.data.frame() 

otu_bact_new_filt <- 
  otu_bact_new_filt[, colSums(otu_bact_new_filt)>0]

str(otu_bact_new_filt)
str(alpha_df_bact_filt)

# reorder
rownames(alpha_df_bact_filt) <- alpha_df_bact_filt$rowname
identical(rownames(alpha_df_bact_filt), rownames(otu_bact_new_filt))
order_bact <- match(rownames(otu_bact_new_filt), rownames(alpha_df_bact_filt))
alpha_df_bact_filt <- alpha_df_bact_filt[order_bact,]

permdisp_bact_soil <- 
  betadisper(vegdist(otu_bact_new_filt, method = "bray"), alpha_df_bact_filt$Soil) 
dist_bact_soil<-
  data.frame(permdisp_bact_soil$group, permdisp_bact_soil$distances)
colnames(dist_bact_soil) <- c("value", "dispSoil")
dist_bact_soil

permdisp_bact_genotype <- 
  betadisper(vegdist(otu_bact_new_filt, method = "bray"), alpha_df_bact_filt$Genotype) 
dist_bact_genotype<-
  data.frame(permdisp_bact_genotype$group, permdisp_bact_genotype$distances)
colnames(dist_bact_genotype) <- c("value", "dispGenotype")
dist_bact_genotype

identical(rownames(dist_bact_genotype), rownames(dist_bact_soil))

meta_bact_merged <-
  cbind(dist_bact_genotype, dist_bact_soil) %>%
  dplyr::select(dispGenotype, dispSoil) %>%
  tibble::rownames_to_column() %>%
  left_join(alpha_df_bact_filt, by="rowname") %>%
  left_join(tibble::rownames_to_column(meta_bact_new[,-1]), by="rowname")
str(meta_bact_merged)


# transform to z-score -------------------------------------------------------------------------------------
mod_bact_1 <- meta_bact_merged[, c("Pot",
                                   "richness", 
                                   "shannon", 
                                   "PCoA.1",
                                   "PCoA.2",
                                   "dispSoil",
                                   "dispGenotype", 
                                   "readNo")]

VarStand(mod_bact_1)

mod_bact_2 <- meta_bact_merged[, c("Pot",
                                   "aerial_part_dry_weight_16wks_grams")]
VarStand(mod_bact_2)

mod_bact_3 <- meta_bact_merged[, c("Pot",
                                   "pH", 
                                   "P", 
                                   "K",
                                   "Ca", 
                                   "Mg", 
                                   "OM", 
                                   "NO3")]

VarStand(mod_bact_3)
mod_bact_3 %>% mutate_if(is.numeric, as.factor)

mod_bact_4 <- meta_bact_merged[, c("Pot",
                                   #"Flower_12wk", 
                                   #"Stress_12wks",
                                   #"Stage_12wks",
                                   #"developmental.estage_16wkE", 
                                   #"developmental.stage_16wks_detail", 
                                   #"Development_16wks_ERS", 
                                   #"Flower_16wk",
                                   "Soil",
                                   "Genotype",
                                   "Ecotype")]

# RF with chemistry variables as numeric -------------------------------------------------------------------
# physeq_bact_new@sam_data$SoilGen <-
#   paste(physeq_bact_new@sam_data$Soil, physeq_bact_new@sam_data$Genotype)

meta_bact_filt_2 <-
  left_join(mod_bact_2, 
            VarStand(mod_bact_1), by="Pot") %>%
  left_join(VarStand(mod_bact_3), by="Pot") %>%
  left_join(mod_bact_4, by="Pot")

# remove Pot and remove 1 sample with too many NA
meta_bact_filt_2 <-
  meta_bact_filt_2[,-1][complete.cases(meta_bact_filt_2[,-1])==TRUE, ]
str(meta_bact_filt_2)

#recoding varaiables 
# meta_bact_filt_2$Stress_12wks <-
#   as.factor(ifelse(meta_bact_filt_2$Stress_12wks=="", "N", paste(meta_bact_filt_2$Stress_12wks)))

# rename colnames
colnames(meta_bact_filt_2) <-c("Biomass",
                               "Richness","Shannon",
                               "PCoA.1","PCoA.2","Disp.Soil","Disp.Genotype",
                               "Read.No",
                               "pH","P","K","Ca","Mg","OM","NO3",
                               #"Stress.12wks",
                               "Soil","Genotype","Ecotype")
head(meta_bact_filt_2)

# loop to store 100 boruta runs
sign_var_bact <- vector(mode = "character")

for(i in 1:99) {
  sel_attr[i]  <- Boruta(Biomass ~., meta_bact_filt_2, pValue = 0.05, 
                         mcAdj = TRUE, maxRuns=100, doTrace = 3)
  sign_var_bact <- append(sign_var_bact, getSelectedAttributes(sel_attr, withTentative = TRUE))
}

sign_var_bact
unique(sign_var_bact)
unique(sign_var)

df_bact_RF_2 <-
  meta_bact_filt_2[,c(unique(sign_var_bact), "Biomass")]

# try tuning the model first
round(sqrt(ncol(df_bact_RF_2[, 1:(ncol(df_bact_RF_2) - 1)])))

set.seed(12345)

bestmtry_bact_Bmass_2 <-
  tuneRF(
    x = df_bact_RF_2[, 1:(ncol(df_bact_RF_2) - 1)],
    y = df_bact_RF_2$Biomass,
    mtryStart = 4,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )


RF_bact_Bmass_2 <-
  randomForest(
    x = df_bact_RF_2[, 1:(ncol(df_bact_RF_2) - 1)],
    y = df_bact_RF_2$Biomass,
    ntree = 1001,
    mtry = 2,
    importance = TRUE,
    proximity = TRUE
  )

RF_bact_Bmass_2
plot(RF_bact_Bmass_2)

# Assessing model significance using permutations
set.seed(110324)
perm_RF_bact_Bmass_2 <-
  rf.significance(
    x = RF_bact_Bmass_2,
    xdata = df_bact_RF_2[, 1:(ncol(df_bact_RF_2) - 1)],
    nperm = 999,
    nmtry = 2,
    ntree = 1001
  )

perm_RF_bact_Bmass_2 # model significant = 0.001 

# Changing tick labels of the bar plot ---------------------------------------------------------------------
sort(rownames(RF_fungi_Bmass_2$importance))
sort(rownames(RF_bact_Bmass_2$importance))

PlotFeature(RF_fungi_Bmass_2, "IncMSE") +
  scale_x_discrete(labels=c(
    "Richness" = "Richness", "Shannon" = "Shannon",
    "PCoA.1" = "PCoA.Axis1", "PCoA.2"="PCoA.Axis2",
    "Disp.Soil"="Disp.Soil", "Read.No"="Read.No",
    "pH"="pH",
    "P" = expression(PO[4]^{"-"}),
    "K" = expression(K^{"+"}),
    "Ca" =expression(Ca^{"2+"}),
    "Mg" =expression(Mg^{"2+"}),
    "OM"="OM",
    "NO3" = expression(NO[3]^{"-"}) ,
    "Soil"="Soil","Genotype"="Genotype", "Ecotype"="Ecotype"))

sort(rownames(RF_bact_Bmass_2$importance))

PlotFeature(RF_bact_Bmass_2, "IncMSE") +
  scale_x_discrete(labels=c(
    "Richness"="Richness", "Shannon" = "Shannon",
    "PCoA.1" = "PCoA.Axis1", "PCoA.2"="PCaA.Axis2",
    "Disp.Soil" = "Disp.Soil" ,"Disp.Genotype" = "Disp.Genotype", 
    "Read.No" = "Read.No", "pH"="pH",
    "P" = expression(PO[4]^{"-"}),
    "K" = expression(K^{"+"}),
    "Ca" =expression(Ca^{"2+"}),
    "Mg" =expression(Mg^{"2+"}),
    "OM" = "OM",
    "NO3" = expression(NO[3]^{"-"}) , 
    "Soil" = "Soil", "Genotype"="Genotype", "Ecotype"="Ecotype"))


# FIGURE 3 - RF MODELS fungi -----------------------------------------------------------------------------
# adjust colnames

pca_fungi <-
  PlotNeutralPCAfungi(meta_fungi_filt_2) +
  scale_color_manual(values = Pal_soil) +
  scale_fill_manual(values = Pal_soil) +
  scale_shape_manual(values = c(0,1,2,5,3,8)) + 
  theme(legend.position = "none") 

# modify loading thickness
pca_fungi$layers[[2]]$aes_params$size <- 0.3
pca_fungi$layers[[2]]$geom_params$arrow$length <- unit(6, units = "points")
pca_fungi

Fig_3_RF_fungi <-
  ggarrange(PlotLine(RF_fungi_Bmass_2, df_fungi_RF_2) + 
              #theme(legend.position = c(0.35, 0.8)) +
              theme(legend.position = "none") +
              scale_color_manual(values = Pal_soil) +
              scale_shape_manual(values = c(0,1,2,5,3,8)) +
              labs(title = "Random Forest ", y="Predicted\nAerial Plant Biomass", x="Observed\nAerial Plant Biomass") +
              annotate("text", x=Inf, y = -Inf, colour="black",
                       label=paste("Mean squared error:", round(last(RF_fungi_Bmass_2$mse), 2)), size=2.5, vjust=-6, hjust=1) +
              annotate("text", x=Inf, y = -Inf, colour="black",
                       label= paste("% Var explained:", round(last(RF_fungi_Bmass_2$rsq*100),2)), size=2.5, vjust=-4, hjust=1) +
              annotate("text", x=Inf, y = -Inf, colour="black",
                       label = paste("italic(p) ==", round(perm_RF_fungi_Bmass_2$pValue, 4)), parse = TRUE, size=2.5, vjust=-1.5, hjust=1),
            PlotFeature(RF_fungi_Bmass_2, "IncMSE") +
              scale_x_discrete(labels=c(
                "Richness" = "Richness", "Shannon" = "Shannon",
                "PCoA.1" = "PCoA.1", "PCoA.2"="PCoA.2",
                "Disp.Soil"="Disp.Soil", "Read.No"="Read.No",
                "pH"="pH",
                "P" = expression(PO[4]^{"3-"}),
                "K" = expression(K^{"+"}),
                "Ca" =expression(Ca^{"2+"}),
                "Mg" =expression(Mg^{"2+"}),
                "OM"="OM",
                "NO3" = expression(NO[3]^{"-"}) ,
                "Soil"="Soil","Genotype"="Genotype", "Ecotype"="Ecotype")) +
              labs(title = "Predictors", x=NULL, y="% Increase MSE ") +
              scale_y_continuous(labels = scales::number_format(accuracy = 0.01)),
            pca_fungi,
            labels = c("A","B","C"),
            widths = c(1, 0.65, 1),
            align = "h",
            ncol = 3, 
            nrow =1)

Fig_3_RF_fungi

# FIGURE SXX - Random Forest supplementary ---------------------------------------------------------------

Fig_SXX_fungi <-
  ggarrange(
    PlotFeature(RF_fungi_Bmass_2,  "IncNodePurity") +
      labs(title = "Predictors", x=NULL, y="Increase Node Purity"),
    PlotError(RF_fungi_Bmass_2) +
      annotate("text", x=Inf, y = Inf, 
               label= paste("% Explained var.:", round(last(RF_fungi_Bmass_2$rsq*100),2)), size=2.5, vjust=1, hjust=1) +
      annotate("text", x=Inf, y = Inf,
               label=paste("MSE:", round(last(RF_fungi_Bmass_2$mse), 2)), size=2.5, vjust=3, hjust=1) +
      annotate("text", x=Inf, y = Inf, 
               label = paste("italic(p) ==", round(perm_RF_fungi_Bmass_2$pValue, 4)), parse = TRUE, size=2.5, vjust=4, hjust=1),
    labels = c("A","B"),
    widths = c(1, 1),
    align = "h",
    ncol = 2, 
    nrow =1)

Fig_SXX_fungi


# FIGURE 3 - RF MODELS bacteria-----------------------------------------------------------------------------

pca_bact <-
  PlotNeutralPCAbact(meta_bact_filt_2) +
  scale_color_manual(values = Pal_soil) +
  scale_fill_manual(values = Pal_soil) +
  scale_shape_manual(values = c(0,1,2,5,3,8)) +
  theme(legend.position = "none")

# modify loading thickness
pca_bact$layers[[2]]$aes_params$size <- 0.3
pca_bact$layers[[2]]$geom_params$arrow$length <- unit(6, units = "points")
pca_bact


Fig_3_RF_bact <-
  ggarrange(PlotLine(RF_bact_Bmass_2, df_bact_RF_2) + 
              theme(legend.position = "none") +
              scale_color_manual(values = Pal_soil) +
              scale_shape_manual(values = c(0,1,2,5,3,8)) +
              labs(title = "Random Forest ", y="Predicted\nAerial Plant Biomass", x="Observed\nAerial Plant Biomass") +
              annotate("text", x=Inf, y = -Inf, color="black",
                       label=paste("Mean squared error:", round(last(RF_bact_Bmass_2$mse), 2)), size=2.5, vjust=-6, hjust=1) +
              annotate("text", x=Inf, y = -Inf, 
                       label= paste("% Var explained:", round(last(RF_bact_Bmass_2$rsq*100),2)), size=2.5, vjust=-4, hjust=1) +
              annotate("text", x=Inf, y = -Inf, 
                       label = paste("italic(p) ==", round(perm_RF_bact_Bmass_2$pValue, 4)), parse = TRUE, size=2.5, vjust=-1.5, hjust=1) +
              theme(legend.position = "none"),
            PlotFeature(RF_bact_Bmass_2, "IncMSE") +
              scale_x_discrete(labels=c(
                "Richness"="Richness", "Shannon" = "Shannon",
                "PCoA.1" = "PCoA.1", "PCoA.2"="PCoA.2",
                "Disp.Soil" = "Disp.Soil" ,"Disp.Genotype" = "Disp.Genotype", 
                "Read.No" = "Read.No", "pH"="pH",
                "P" = expression(PO[4]^{"3-"}),
                "K" = expression(K^{"+"}),
                "Ca" =expression(Ca^{"2+"}),
                "Mg" =expression(Mg^{"2+"}),
                "OM" = "OM",
                "NO3" = expression(NO[3]^{"-"}) , 
                "Soil" = "Soil", "Genotype"="Genotype", "Ecotype"="Ecotype")) +
              labs(title = "Predictors", x=NULL, y="% Increase MSE ") +
              scale_y_continuous(labels = scales::number_format(accuracy = 0.01)),
            pca_bact,
            labels = c("D","E","F"),
            widths = c(1, 0.65, 1),
            align = "h",
            ncol = 3, 
            nrow =1)

Fig_3_RF_bact

# FIGURE SXX - Random Forest supplementary ---------------------------------------------------------------

Fig_SXX_bact <-
  ggarrange(
    PlotFeature(RF_bact_Bmass_2,  "IncNodePurity") +
      labs(title = "Predictors", x=NULL, y="Increase Node Purity"),
    PlotError(RF_bact_Bmass_2) +
      annotate("text", x=Inf, y = Inf, 
               label= paste("% Explained var.:", round(last(RF_bact_Bmass_2$rsq*100),2)), size=2.5, vjust=1, hjust=1) +
      annotate("text", x=Inf, y = Inf,
               label=paste("MSE:", round(last(RF_bact_Bmass_2$mse), 2)), size=2.5, vjust=3, hjust=1) +
      annotate("text", x=Inf, y = Inf, 
               label = paste("italic(p) ==", round(perm_RF_bact_Bmass_2$pValue, 4)), parse = TRUE, size=2.5, vjust=4, hjust=1),
    labels = c("A","B"),
    widths = c(1, 1),
    align = "h",
    ncol = 2, 
    nrow =1)

Fig_SXX_bact


# *** FINAL FIGURE 3 ---------------------------------------------------------------------------------------

# extracting the legend for plotting 
get_legend(
  PlotLine(RF_bact_Bmass_2, df_bact_RF_2) + 
    theme(legend.position = c(0.2, 0.8)) +
    scale_color_manual(values = Pal_soil) +
    guides(color = guide_legend(nrow = 2),
           shape = guide_legend(nrow = 2)) +
    scale_shape_manual(values = c(0,1,2,5,3,8), 
                       labels = c("Alamo"=expression(bold("Alamo")), 
                                  "Blackwell", 
                                  "Cave-in-rock", 
                                  "Kanlow"=expression(bold("Kanlow")), 
                                  "Shelter",
                                  "Southlow")) +
    theme(legend.position = "bottom")) -> legend1

as_ggplot(legend1)


title3 = text_grob("Fungi", size = 12, face = 2)
title4 = text_grob("Bacteria", size = 12, face = 2)


ggarrange(
  grid.arrange(Fig_3_RF_fungi, top=title3),
  grid.arrange(Fig_3_RF_bact, top=title4),
  as_ggplot(legend1),
  ncol = 1,
  nrow = 3,
  heights = c(1,1,0.1))
  

# Treating as a factor -------------------------------------------------------------------------------------
meta_bact_filt <-
  left_join(mod_bact_2, 
            VarStand(mod_bact_1), by="Pot") %>%
  left_join(mutate_if(mod_bact_3, is.numeric, as.factor), by="Pot") %>%
  left_join(mod_bact_4, by="Pot")

# remove Pot and remove 1 sample with too many NA
meta_bact_filt <-
  meta_bact_filt[,-1][complete.cases(meta_bact_filt[,-1])==TRUE, ]
str(meta_bact_filt)

#recoding varaiables 
meta_bact_filt$Stress_12wks <-
  as.factor(ifelse(meta_bact_filt$Stress_12wks=="", "N", paste(meta_bact_filt$Stress_12wks)))


# loop to store 100 boruta runs
sign_var <- vector(mode = "character")

for(i in 1:99) {
  sel_attr[i]  <- Boruta(aerial_part_dry_weight_16wks_grams ~., meta_bact_filt, pValue = 0.05, 
                         mcAdj = TRUE, maxRuns=100, doTrace = 3)
  sign_var <- append(sign_var, getSelectedAttributes(sel_attr, withTentative = TRUE))
}

sign_var
unique(sign_var)

df_bact_RF <-
  meta_bact_filt[,c(unique(sign_var), "aerial_part_dry_weight_16wks_grams")]

# try tuning the model first
round(sqrt(ncol(df_bact_RF[, 1:(ncol(df_bact_RF) - 1)])))

set.seed(12345)

bestmtry_bact_Bmass <-
  tuneRF(
    x = df_bact_RF[, 1:(ncol(df_bact_RF) - 1)],
    y = df_bact_RF$aerial_part_dry_weight_16wks_grams,
    mtryStart = 4,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )

RF_bact_Bmass <-
  randomForest(
    x = df_bact_RF[, 1:(ncol(df_bact_RF) - 1)],
    y = df_bact_RF$aerial_part_dry_weight_16wks_grams,
    ntree = 1001,
    mtry = 2,
    importance = TRUE,
    proximity = TRUE
  )

RF_bact_Bmass
plot(RF_bact_Bmass)

# Assessing model significance using permitations
set.seed(110324)
perm_RF_bact_Bmass <-
  rf.significance(
    x = RF_bact_Bmass,
    xdata = df_bact_RF[, 1:(ncol(df_bact_RF) - 1)],
    nperm = 999,
    nmtry = 1,
    ntree = 1001
  )

perm_RF_bact_Bmass # model significant = 0.001 



# *****************************************************************-----------------------------------------
# INFLUENTIAL OTUs on plant biomass ------------------------------------------------------------------------
# Fungi ----------------------------------------------------------------------------------------------------

str(otu_fungi_new)
str(meta_fungi_new)
identical(rownames(otu_fungi_new), rownames(meta_fungi_new))

df_fungi_RF_otu <-
  data.frame(otu_fungi_new, 
             biomass = meta_fungi_new$aerial_part_dry_weight_16wks_grams) 

str(df_fungi_RF_otu)
head(df_fungi_RF_otu)

# Recursive feature selection
set.seed(110321)

rfe_fungi_biom <- Boruta(
  biomass ~ .,
  df_fungi_RF_otu,
  pValue = 0.05,
  mcAdj = TRUE,
  maxRuns = ncol(df_fungi_RF_otu),
  doTrace = 3
)

rfe_fungi_biom

# Get significant variables including tentatives
impfe_fungi_biom <-
  getSelectedAttributes(rfe_fungi_biom, withTentative = TRUE)
impfe_fungi_biom

#subset to important features only
df_fungi_RF_otu[, c(impfe_fungi_biom)] -> df_fungi_RF_otu_sel
identical(rownames(df_fungi_RF_otu_sel), rownames(df_fungi_RF_otu))
df_fungi_RF_otu_sel$biomass <- df_fungi_RF_otu$biomass

# select optimal mtry
round(sqrt(ncol(df_fungi_RF_otu_sel[, 1:(ncol(df_fungi_RF_otu_sel) - 1)])))

set.seed(110322)
bestmtry_fungi_biom <-
  tuneRF(
    x = df_fungi_RF_otu_sel[, 1:(ncol(df_fungi_RF_otu_sel) - 1)],
    y = df_fungi_RF_otu_sel$biomass,
    mtryStart = 7,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )

bestmtry_fungi_biom

set.seed(110323)
RF_fungi_biom <-
  randomForest(
    x = df_fungi_RF_otu_sel[, 1:(ncol(df_fungi_RF_otu_sel) - 1)],
    y = df_fungi_RF_otu_sel$biomass,
    ntree = 1001,
    mtry = 7,
    importance = TRUE,
    proximity = TRUE
  )

RF_fungi_biom
plot(RF_fungi_biom)

# Assessing model significance using permitations
set.seed(110324)
perm_RF_fungi_biom <-
  rf.significance(
    x = RF_fungi_biom,
    xdata = df_fungi_RF_otu_sel[, 1:(ncol(df_fungi_RF_otu_sel) - 1)],
    nperm = 999,
    nmtry = 1,
    ntree = 1001
  )

perm_RF_fungi_biom # model significant = 0.001 


# Dataframe for PCA 
df_fungi_RF_otu_sel_pca <- 
  left_join(tibble::rownames_to_column(df_fungi_RF_otu_sel), #keep this samples only
            tibble::rownames_to_column(as(physeq_fungi_new@sam_data, "data.frame") %>%
                                         dplyr::select(Genotype, Soil_location) %>%
                                         rename(Soil = Soil_location)), by="rowname")
head(df_fungi_RF_otu_sel_pca)
colnames(df_fungi_RF_otu_sel_pca)[56] <- "Biomass"


# Plotting the model... Need a different function for plotting features.

# New plotting functions -----------------------------------------------------------------------------------

PlotOTU <- function(rf_model, taxa, Var){
  require(tibble)
  imp_RF <- as.data.frame(rf_model$importance)
  imp_RF$features <- rownames(imp_RF)
  imp_RF <- arrange(imp_RF, desc(Var))
  rownames(imp_RF) <- imp_RF$features
  # adding marker taxon info
  taxa[rownames(taxa)%in%rownames(imp_RF), ] -> taxa_RF
  identical(rownames(taxa_RF), rownames(imp_RF))
  order_taxa <- match(rownames(taxa_RF), rownames(imp_RF))
  imp_RF <- imp_RF[order_taxa,]
  imp_RF$Taxonomy <- taxa_RF$Taxon
  imp_RF <- imp_RF[order(imp_RF[,Var], decreasing = TRUE),] 
  imp_RF <- left_join(tibble::rownames_to_column(imp_RF),
                      tibble::rownames_to_column(
                        taxa[rownames(imp_RF), c(9:11, 16)]))
  imp_RF$Taxonomy <-
    gsub("FOTU_", "", imp_RF$Taxonomy)
  imp_RF$Taxonomy <-
    gsub("OTU_", "",imp_RF$Taxonomy)
  imp_RF$glbrc <-
    ifelse(is.na(imp_RF$Isolate), "no", "yes")
  imp_RF %T>% print()
  ggplot(data=imp_RF) + 
    geom_bar(aes(x= reorder(Taxonomy, -get(Var)),
                 y= get(Var),
                 color= glbrc, fill=glbrc), stat="identity") +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0)) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.y = element_text(angle =0, size = 8, hjust = 0.5, vjust = 1), 
          axis.text.x = element_text(angle = 90, size = 7 ,hjust = 1, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8), 
          plot.margin=unit(c(1.5,1.5,1.5,1.5),"pt"),
          legend.position = "none") -> plot_importance
  return(plot_importance)
}


PlotOTU(RF_fungi_biom, taxa_fungi_new, "IncNodePurity") + 
  labs(title = "Fungal OTUs affecting areal plant biomass",
       x= NULL, y= "Node Purity")

PlotOTU(RF_fungi_biom, taxa_fungi_new, "%IncMSE") + 
  labs(title = "Fungal OTUs affecting areal plant biomass",
       x= NULL, y= "% Increase MSE")


PlotLineBiom <- function(rf_model, metadata){
  df_model <- 
    left_join(tibble::rownames_to_column(
      data.frame(actual = rf_model$y,
                 pred = rf_model$predicted)), 
      metadata %>%
        dplyr::select(rowname, Genotype, Soil), by = "rowname")
    df_model %T>% print()
    ggplot(data=df_model, aes(x=actual, y=pred)) +
      geom_point(aes(shape=Genotype, color=Soil), size=1.5, stroke=0.5) +
      geom_smooth(data=df_model, method = "lm", formula = "y ~ x", se = TRUE, color="black", size=0.5) +
      theme_classic() +
      scale_y_continuous(labels = scales::number_format(accuracy = 1)) +
      theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
            axis.title = element_text(angle = 0, size = 10, face = "bold"),
            axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 1), 
            axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
            legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
            legend.title = element_blank(), legend.background = element_blank(), 
            legend.text = element_text(size = 8)) +
      guides(color = guide_legend(nrow = 2),shape = guide_legend(nrow = 2)) -> line_plot
    return(line_plot)
}

PlotLineBiom(RF_fungi_biom, df_fungi_RF_otu_sel_pca) + theme(legend.position = c(0.1, 0.8))

PlotPCAotu <- function(df, Var){
  colnames(df) <-
    gsub("FOTU_", "", colnames(df))
  colnames(df) <-
    gsub("OTU_", "", colnames(df))
  pca_plot <-
    autoplot(
      prcomp(x = df, scale= TRUE, center=TRUE), # calculates principal compenents and pltos with ggplot2
        data = Var, label = FALSE, shape = "Genotype", colour="Soil", # add metadate, labels of objects
        loadings = TRUE, loadings.colour = "black", size=1.5, stroke=0.8,
        max.overlaps = getOption("ggrepel.max.overlaps", default = 0),
      frame = FALSE, frame.colour = "Soil", loadings.label.colour = "black",
      loadings.label = TRUE, loadings.label.size = 1.5, loadings.label.repel = TRUE) +
    labs(title = "PCA") +
    # scale_colour_manual(values = paletteCB4) +
    # scale_fill_manual(values = paletteCB4) +
    # scale_shape_manual(values = c(21,22,24)) +
    theme_classic() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.x = element_text(angle =0, size = 8, hjust = 0.5, vjust = 1), 
          axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8))
  # guides(color = guide_legend(ncol = 2), #title.position="top"
  #         fill = guide_legend(ncol= 2),
  #         shape = guide_legend(ncol = 1)) +
  # theme(legend.margin = margin(0,-0.5,0,0, unit="cm")) # reduce space betweem legends
  return(pca_plot)
}


str(df_fungi_RF_otu_sel_pca)
PlotPCAotu(df_fungi_RF_otu_sel_pca[,2:56], df_fungi_RF_otu_sel_pca[, 57:58])


# Extract feaures seq for BLAST ----------------------------------------------------------------------------
filterTaxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

write.dna(refseq(filterTaxa(physeq_fungi_new, 
                            rownames(taxa_fungi_new[colnames(df_fungi_RF_otu_sel), ])[-55])),
          format="fasta", 
          colsep="", 
          file="fungi_plant_biomass.fasta")


write.dna(refseq(filterTaxa(physeq_bact_new, 
                            rownames(taxa_bact_new[colnames(df_bact_RF_otu_sel), ])[-55])),
          format="fasta", 
          colsep="", 
          file="bact_plant_biomass.fasta")




# Correcting taxonomy tables ------------------------------------------------------------------------------ 

taxa_fungi_new[colnames(df_fungi_RF_otu_sel), ][, c(1,9,16)]
taxa_fungi_new$Taxonomy <- as.character(taxa_fungi_new$Taxonomy)

taxa_fungi_new[taxa_fungi_new == "FOTU_492-Glomeraceae"] <- "FOTU_492-Rhizophagus irregularis"
taxa_fungi_new[taxa_fungi_new == "FOTU_1-Ascomycota"] <- "FOTU_1-Setophoma sp."
taxa_fungi_new[taxa_fungi_new == "FOTU_191-Glomeraceae"] <- "FOTU_191-Glomus sp."
taxa_fungi_new[taxa_fungi_new == "FOTU_126-Chaetothyriales"] <- "FOTU_191-Knufia sp."
taxa_fungi_new[taxa_fungi_new == "FOTU_10071-Glomeraceae"] <- "FOTU_10071-Rhizophagus sp."
taxa_fungi_new[taxa_fungi_new == "FOTU_139-Branch06"] <- "FOTU_139-Sordariomycetes"
taxa_fungi_new[taxa_fungi_new == "FOTU_2313-Glomeraceae"] <- "FOTU_2313-Glomus sp."
taxa_fungi_new[taxa_fungi_new == "FOTU_57-Hypocreales"] <- "FOTU_57-Fusarium sp."
taxa_fungi_new[taxa_fungi_new == "FOTU_690-Fungi"] <- "FOTU_690-Glomus sp."
taxa_fungi_new[taxa_fungi_new == "FOTU_258-Glomus"] <- "FOTU_258-Glomus sp."
taxa_fungi_new[taxa_fungi_new == "FOTU_3189-Glomeraceae"] <- "FOTU_258-Glomus cf. macrocarpum"
taxa_fungi_new[taxa_fungi_new == "FOTU_725-Glomeraceae"] <- "FOTU_725-Glomus sp."
taxa_fungi_new[taxa_fungi_new == "FOTU_1183-Glomeraceae"] <- "FOTU_1183-Glomus sp."
taxa_fungi_new[taxa_fungi_new == "FOTU_772-Glomeraceae"] <- "FOTU_772-Glomus sp."
taxa_fungi_new[taxa_fungi_new == "FOTU_967-Piriformospora sp."] <- "FOTU_967-Serendipita sp."
taxa_fungi_new[taxa_fungi_new == "FOTU_10512-Glomeraceae"] <- "FOTU_10512-Septoglomus sp."
taxa_fungi_new[taxa_fungi_new == "FOTU_10558-Glomeraceae"] <- "FOTU_10558-Glomus sp."
taxa_fungi_new[taxa_fungi_new == "FOTU_4264-Glomeraceae"] <- "FOTU_4264-Glomus sp."


taxa_bact_new[colnames(df_bact_RF_otu_sel), ][, c(1,9,16)]
taxa_bact_new$Taxonomy <- as.character(taxa_bact_new$Taxonomy)

taxa_bact_new[taxa_bact_new == "OTU_64-67-14"] <- "OTU_64-Solirubrobacter sp."
taxa_bact_new[taxa_bact_new == "OTU_516-Ffch7168 sp."] <- "OTU_516-Bacteria"
taxa_bact_new[taxa_bact_new == "OTU_306-Uncultured7 sp."] <- "OTU_306-Rhizobiales"
taxa_bact_new[taxa_bact_new == "OTU_2263-Uncultured 60"] <- "OTU_2263-Acidobacteria"
taxa_bact_new[taxa_bact_new == "OTU_754-Uncultured 76 sp."] <- "OTU_754-Gemmatimonas sp."
taxa_bact_new[taxa_bact_new == "OTU_568-Kd4-96"] <- "OTU_568-Acidobacteria"
taxa_bact_new[taxa_bact_new == "OTU_7284-Uncultured 76 sp."] <- "OTU_7284-Gemmatimonadetes"
taxa_bact_new[taxa_bact_new == "OTU_232-type III"] <- "OTU_232-Cand. Moen. glomeromycotorum"
taxa_bact_new[taxa_bact_new == "OTU_2472-Allorhizobium-neorhizobium-pararhizobium-rhizobium sp."] <- "OTU_2472-Rhizobium sp."
taxa_bact_new[taxa_bact_new == "OTU_427-67-14"] <- "OTU_427-Actinobacteria"
taxa_bact_new[taxa_bact_new == "OTU_451-Plta13"] <- "OTU_451-Xanthomonadales"
taxa_bact_new[taxa_bact_new == "OTU_4189-Uncultured51 sp."] <- "OTU_4189-Rhizobiales"
taxa_bact_new[taxa_bact_new == "OTU_796-Uncultured bacterium 2214 sp."] <- "OTU_796-Methylocystaceae"
taxa_bact_new[taxa_bact_new == "OTU_20144-D05-2"] <- "OTU_20144-Bacteria"
taxa_bact_new[taxa_bact_new == "OTU_4834-Uncultured 77 sp."] <- "OTU_4834-Myxococcales"
taxa_bact_new[taxa_bact_new == "OTU_727-Uncultured mollicutes bacterium 8 sp."] <- "OTU_727-Cand. Moen. glomeromycotorum"
taxa_bact_new[taxa_bact_new == "OTU_962-A4b"] <- "OTU_962-Bacteria"
taxa_bact_new[taxa_bact_new == "OTU_2597-Uncultured 87 sp."] <- "OTU_2597-Rhizobiales"
taxa_bact_new[taxa_bact_new == "OTU_6199-Uncultured31 sp."] <- " OTU_6199-Bacteria"
taxa_bact_new[taxa_bact_new == "OTU_760-Uncultured36 sp."] <- "OTU_760-Burkholderiales"
taxa_bact_new[taxa_bact_new == "OTU_700-Uncultured2 sp."] <- "OTU_700-Blastopirellula sp."
taxa_bact_new[taxa_bact_new == "OTU_6802-Sm2d12"] <- "OTU_6802-Bacteria"
taxa_bact_new[taxa_bact_new == "OTU_699-Uncultured 76 sp."] <- "OTU_699-Gemmatimonas sp."

# NEUTRAL MODELS ------------------------------------------------------------------------------------
library(tyRa)
require(minpack.lm)
require(Hmisc)
require(stats4)

neutral_model_fungi <- 
  tyRa::fit_sncm(spp = t(otu_table(physeq_fungi_new)), pool=NULL, taxon=data.frame(tax_table(physeq_fungi_new)))

plot_sncm_fit(neutral_model_fungi, fill = NULL, title = "Model Fit") +
  theme_classic()

neutral_model_bact <- 
  tyRa::fit_sncm(spp = t(otu_table(physeq_bact_new)), pool=NULL, taxon=data.frame(tax_table(physeq_bact_new)))

plot_sncm_fit(neutral_model_bact, fill = NULL, title = "Model Fit") +
  theme_classic()


# INDICATOR SPECIES ANALYSIS ------------------------------------------------------------------
library(indicspecies)

GetIndicators <-function(dataframe, var){
  require(phyloseq); require(indicspecies); require(dplyr)
  otu <- as.data.frame(otu_table(dataframe))
  metadata = as(sample_data(dataframe), "data.frame")
  taxa <- as.data.frame(as.matrix(tax_table(dataframe)))
  multipatt <- multipatt(t(otu), metadata[,var], func = "r.g",
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


# indicator value >0.5 and p-value <0.05 after fdr correction
head(physeq_fungi_new@sam_data)

ind_ITS_soil <-
  GetIndicators(
    physeq_fungi_new %>%
      subset_samples(Root_soil%in%"Root"),
    "Soil_location")

head(ind_ITS_soil)
dim(ind_ITS_soil)

ind_16s_soil <-
  GetIndicators(
    physeq_bact_new %>%
      subset_samples(Root_soil%in%"Root"),
    "Soil_location")

head(ind_16s_soil)
dim(ind_16s_soil)

# RF AND INDICATORS MATCH ---------------------------------------------------------------------
GetIndTab <- function(df_ind) {
  df_ind$IndVal <- 
    df_ind$index %>%
    recode_factor(
      "1" = "Hancock",
      "2" = "Lake City",
      "3" = "Lux Arbor",
      "4" = "Rhineland")
  return(df_ind[, c(1,7,8,25)])
}

# Match indicators with RF models features accuracy -------------------------------------------
ind_ITS_soil_tab <-
  GetIndTab(ind_ITS_soil)

ind_16s_soil_tab <-
  GetIndTab(ind_16s_soil)


PlotOTUInd <- function(rf_model, taxa, NeutMod, Var, df_ind, df_name){
  require(tibble)
  imp_RF <- as.data.frame(rf_model$importance)
  imp_RF$features <- rownames(imp_RF)
  imp_RF <- arrange(imp_RF, desc(Var))
  rownames(imp_RF) <- imp_RF$features
  # adding marker taxon info
  taxa[rownames(taxa)%in%rownames(imp_RF), ] -> taxa_RF
  identical(rownames(taxa_RF), rownames(imp_RF))
  order_taxa <- match(rownames(taxa_RF), rownames(imp_RF))
  imp_RF <- imp_RF[order_taxa,]
  imp_RF$Taxonomy <- taxa_RF$Taxon
  imp_RF <- imp_RF[order(imp_RF[,Var], decreasing = TRUE),] 
  imp_RF <- left_join(tibble::rownames_to_column(imp_RF),
                      tibble::rownames_to_column(
                        taxa[rownames(imp_RF), c(9:11, 16)]))
  imp_RF$Taxonomy <-
    gsub("FOTU_", "", imp_RF$Taxonomy)
  imp_RF$Taxonomy <-
    gsub("OTU_", "",imp_RF$Taxonomy)
  imp_RF$glbrc <-
    ifelse(is.na(imp_RF$Isolate), "no", "yes")
  # adding neutral model results
  df_neutral <-
    as.data.frame(NeutMod[2])
  df_neutral$rowname <- rownames(df_neutral)
  new_df <-
    df_neutral %>%
    dplyr::select("predictions.fit_class", "rowname")
  new_df <-
    left_join(imp_RF, new_df, by= "rowname")
  # adding match with indicators
  colnames(df_ind)[1] <- "rowname"
  final_df <-
    left_join(new_df, df_ind, by= "rowname")
  head(print(final_df))
  # saving intermediate df to R environemnt, pick the right name
  final_df %T>% 
    #assign(paste(df_name, Var1, Var2, sep = ""),., envir = .GlobalEnv) %>% # saving the plot 
    assign(df_name, ., envir = .GlobalEnv) 
  # plotting
  ggplot(data=final_df) + 
    geom_bar(aes(x= reorder(Taxonomy, -get(Var)),
                 y= get(Var),
                 color= IndVal, fill=IndVal), stat="identity") +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0)) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.y = element_text(angle =0, size = 8, hjust = 0.5, vjust = 1), 
          axis.text.x = element_text(angle = 90, size = 7 ,hjust = 1, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8), 
          plot.margin=unit(c(1.5,1.5,1.5,1.5),"pt"),
          legend.position = "none") -> plot_importance
  return(plot_importance)
}


PlotOTUInd(RF_fungi_biom, taxa_fungi_new, neutral_model_fungi, 
           "%IncMSE", ind_ITS_soil_tab, "rf_taxa_fungi") +
  labs(title = "Most important OTUs for plant biomass",
       x= NULL, y= "% Increase MSE") +
  scale_colour_manual(values = Pal_soil) +
  scale_fill_manual(values = Pal_soil)

PlotOTUInd(RF_bact_biom, taxa_bact_new, neutral_model_bact, 
           "%IncMSE", ind_16s_soil_tab, "rf_taxa_bact") +
  labs(title = "Most important OTUs for plant biomass",
       x= NULL, y= "% Increase MSE") +
  scale_colour_manual(values = Pal_soil) +
  scale_fill_manual(values = Pal_soil)



# *** FIGURE 4 ---------------------------------------------------------------------------------------------
pca_fungi_otu <-
  PlotPCAotu(df_fungi_RF_otu_sel_pca[,2:56], df_fungi_RF_otu_sel_pca[, 57:58]) +
    scale_color_manual(values = Pal_soil) +
    scale_fill_manual(values = Pal_soil) +
    scale_shape_manual(values = c(0,1,2,5,3,8)) +
    theme(legend.position = "none")

# modify loading thickness
pca_fungi_otu$layers[[2]]$aes_params$size <- 0.2
pca_fungi_otu$layers[[2]]$geom_params$arrow$length <- unit(4, units = "points")
pca_fungi_otu

Fig_5_RF_fungi_OTU <-
  ggarrange(
    ggarrange(
      PlotLineBiom(RF_fungi_biom, df_fungi_RF_otu_sel_pca) + 
        theme(legend.position = "bottom") +
      theme(legend.position = c(0.15, 0.8)) +
        scale_color_manual(values = Pal_soil) +
        scale_shape_manual(values = c(0,1,2,5,3,8), 
                           labels = c("Alamo", 
                                      "Blackwell", 
                                      "Cave-in-rock", 
                                      "Kanlow", 
                                      "Shelter",
                                      "Southlow")) +
        labs(title = "Random Forest ", y="Predicted\nAerial Plant Biomass", x="Observed\nAerial Plant Biomass") +
        annotate("text", x=Inf, y = -Inf,
                 label=paste("Mean squared error:", round(last(RF_fungi_biom$mse), 2)), size=2.5, vjust=-6, hjust=1) +
        annotate("text", x=Inf, y = -Inf, 
                 label= paste("% Var explained:", round(last(RF_fungi_biom$rsq*100),2)), size=2.5, vjust=-4, hjust=1) +
        annotate("text", x=Inf, y = -Inf, 
                 label = paste("italic(p) ==", round(perm_RF_fungi_biom$pValue, 4)), parse = TRUE, size=2.5, vjust=-1.5, hjust=1),
      pca_fungi_otu,
      labels = c("A","B"),
            widths = c(1, 1),
            align = "h",
            ncol = 2,
            nrow = 1, 
      common.legend = TRUE,
      legend = "bottom"),
    PlotOTUInd(RF_fungi_biom, taxa_fungi_new, neutral_model_fungi, 
               "%IncMSE", ind_ITS_soil_tab, "rf_taxa_fungi") +
      labs(title = "Most important OTUs for plant biomass",
           x= NULL, y= "% Increase MSE") +
      scale_colour_manual(values = Pal_soil) +
      scale_fill_manual(values = Pal_soil),
    labels = c("", "C"),
    heights = c(1, 1),
    ncol = 1, 
    nrow =2)

Fig_5_RF_fungi_OTU

grid.arrange(Fig_5_RF_fungi_OTU, top=title3)


PlotOTU(RF_fungi_biom, taxa_fungi_new, "%IncMSE") + 
  labs(title = "Most important OTUs for plant biomass",
       x= NULL, y= "% Increase MSE") +
  scale_colour_manual(values = c("grey80", "grey40")) +
  scale_fill_manual(values = c("grey80", "grey40"))



# matching root isolaets -----------------------------------------------------------------------------------
isolates_fungi <-
  taxa_fungi_new[colnames(df_fungi_RF_otu_sel), ] %>%
  dplyr::select(Family, Isolate, Isolate_percent_id, Isolate_query_cover,Taxonomy)

isolates_fungi <-
left_join(
  tibble::rownames_to_column(isolates_fungi),
  tibble::rownames_to_column(as.data.frame(RF_fungi_biom$importance)),
  by="rowname")

isolates_fungi <-
  left_join(
    isolates_fungi,
    data.frame(rowname = names(All_otus_ITS), sequence = paste(All_otus_ITS)), 
    by="rowname")

str(isolates_fungi)
write.csv(arrange(isolates_fungi, desc(`%IncMSE`)), "RF_fungi_otu.csv")

# Plot taxonomic proprotions ------------------------------------------------------------
fungi_bar_plot<-
  as.data.frame(
  table(gsub(".*-","",isolates_fungi$Taxonomy))) %>%
  ggplot(aes(Var1, Freq)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(angle = 0, size = 10, face = "bold"),
        axis.text.y = element_text(angle = 0, size = 7, hjust = 1, vjust = 0.5), 
        axis.text.x = element_text(angle = 0, size = 7 ,hjust = 0.5, vjust = 0.5),
        legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
        legend.title = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 8), 
        plot.margin=unit(c(1.5,1.5,1.5,1.5),"pt")) +
  labs(title = "Fungal OTUs", x="Taxon", y="Number of OTUs")
  
fungi_bar_plot


# Bacteria -------------------------------------------------------------------------------------------------

str(otu_bact_new)
str(meta_bact_new)
identical(rownames(otu_bact_new), rownames(meta_bact_new))

df_bact_RF_otu <-
  data.frame(otu_bact_new, 
             biomass = meta_bact_new$aerial_part_dry_weight_16wks_grams) 

str(df_bact_RF_otu)
head(df_bact_RF_otu)

# Recursive feature selection
set.seed(110321)

rfe_bact_biom <- Boruta(
  biomass ~ .,
  df_bact_RF_otu,
  pValue = 0.05,
  mcAdj = TRUE,
  maxRuns = ncol(df_bact_RF_otu),
  doTrace = 3
)

rfe_bact_biom

# Get significant variables including tentatives
impfe_bact_biom <-
  getSelectedAttributes(rfe_bact_biom, withTentative = TRUE)
impfe_bact_biom

#subset to important features only
df_bact_RF_otu[, c(impfe_bact_biom)] -> df_bact_RF_otu_sel
identical(rownames(df_bact_RF_otu_sel), rownames(df_bact_RF_otu))
df_bact_RF_otu_sel$biomass <- df_bact_RF_otu$biomass

# select optimal mtry
round(sqrt(ncol(df_bact_RF_otu_sel[, 1:(ncol(df_bact_RF_otu_sel) - 1)])))

set.seed(110322)
bestmtry_bact_biom <-
  tuneRF(
    x = df_bact_RF_otu_sel[, 1:(ncol(df_bact_RF_otu_sel) - 1)],
    y = df_bact_RF_otu_sel$biomass,
    mtryStart = 7,
    ntreeTry = 1001,
    improve = 0.01,
    stepFactor = 0.5,
    nodesize = 1,
    doBest = TRUE,
    trace = TRUE,
    plot = TRUE
  )

bestmtry_bact_biom

set.seed(110323)
RF_bact_biom <-
  randomForest(
    x = df_bact_RF_otu_sel[, 1:(ncol(df_bact_RF_otu_sel) - 1)],
    y = df_bact_RF_otu_sel$biomass,
    ntree = 1001,
    mtry = 14,
    importance = TRUE,
    proximity = TRUE
  )

RF_bact_biom
plot(RF_bact_biom)

# Assessing model significance using permitations
set.seed(110324)
perm_RF_bact_biom <-
  rf.significance(
    x = RF_bact_biom,
    xdata = df_bact_RF_otu_sel[, 1:(ncol(df_bact_RF_otu_sel) - 1)],
    nperm = 999,
    nmtry = 1,
    ntree = 1001
  )

perm_RF_bact_biom # model significant = 0.001 


# Dataframe for PCA 
df_bact_RF_otu_sel_pca <- 
  left_join(tibble::rownames_to_column(df_bact_RF_otu_sel), #keep this samples only
            tibble::rownames_to_column(as(physeq_bact_new@sam_data, "data.frame") %>%
                                         dplyr::select(Genotype, Soil_location) %>%
                                         rename(Soil = Soil_location)), by="rowname")
str(df_bact_RF_otu_sel_pca)
head(df_bact_RF_otu_sel_pca)

head(df_bact_RF_otu_sel_pca)
colnames(df_bact_RF_otu_sel_pca)[54] <- "Biomass"


# Plotting the model... 

# *** FIGURE 5 ---------------------------------------------------------------------------------------------

pca_bact_otu <-
  PlotPCAotu(df_bact_RF_otu_sel_pca[,2:54], df_bact_RF_otu_sel_pca[, 55:56]) +
  scale_color_manual(values = Pal_soil) +
  scale_fill_manual(values = Pal_soil) +
  scale_shape_manual(values = c(0,1,2,5,3,8)) +
  theme(legend.position = "none")

# modify loading thickness
pca_bact_otu$layers[[2]]$aes_params$size <- 0.25
pca_bact_otu$layers[[2]]$geom_params$arrow$length <- unit(4, units = "points")
pca_bact_otu

Fig_5_RF_bact_OTU <-
  ggarrange(
    ggarrange(
      PlotLineBiom(RF_bact_biom, df_bact_RF_otu_sel_pca) + 
        theme(legend.position = "bottom") +
        theme(legend.position = c(0.15, 0.8)) +
        scale_color_manual(values = Pal_soil) +
        scale_shape_manual(values = c(0,1,2,5,3,8), 
                           labels = c("Alamo", 
                                      "Blackwell", 
                                      "Cave-in-rock", 
                                      "Kanlow", 
                                      "Shelter",
                                      "Southlow"))+
        labs(title = "Random Forest ", y="Predicted\nAerial Plant Biomass", x="Observed\nAerial Plant Biomass") +
        annotate("text", x=Inf, y = -Inf,
                 label=paste("Mean squared error:", round(last(RF_bact_biom$mse), 2)), size=2.5, vjust=-6, hjust=1) +
        annotate("text", x=Inf, y = -Inf, 
                 label= paste("% Var explained:", round(last(RF_bact_biom$rsq*100),2)), size=2.5, vjust=-4, hjust=1) +
        annotate("text", x=Inf, y = -Inf, 
                 label = paste("italic(p) ==", round(perm_RF_bact_biom$pValue, 4)), parse = TRUE, size=2.5, vjust=-1.5, hjust=1),
      pca_bact_otu,
      labels = c("A","B"),
      widths = c(1, 1),
      align = "h",
      ncol = 2,
      nrow = 1, 
      common.legend = TRUE,
      legend = "bottom"),
    PlotOTUInd(RF_bact_biom, taxa_bact_new, neutral_model_bact, 
               "%IncMSE", ind_16s_soil_tab, "rf_taxa_bact") +
      labs(title = "Most important OTUs for plant biomass",
           x= NULL, y= "% Increase MSE") +
      scale_colour_manual(values = Pal_soil) +
      scale_fill_manual(values = Pal_soil),
    labels = c("", "C"),
    heights = c(1, 1),
    ncol = 1, 
    nrow =2)

Fig_5_RF_bact_OTU

grid.arrange(Fig_5_RF_bact_OTU, top=title4)


PlotOTU(RF_bact_biom, taxa_bact_new, "%IncMSE") + 
  labs(title = "Most important OTUs for plant biomass",
       x= NULL, y= "% Increase MSE") +
  scale_colour_manual(values = c("grey80", "grey40")) +
  scale_fill_manual(values = c("grey80", "grey40"))


# matching root isolaets 
isolates_bact <-
taxa_bact_new[colnames(df_bact_RF_otu_sel), ] %>%
  dplyr::select(Class, Isolate, `Isolate Isolate_percent_id`, Isolate_query_cover, Taxonomy) 

isolates_bact <-
  left_join(
    tibble::rownames_to_column(isolates_bact),
    tibble::rownames_to_column(as.data.frame(RF_bact_biom$importance)),
    by="rowname")

isolates_bact <-
  left_join(
    isolates_bact,
    data.frame(rowname = names(All_otus_16S), sequence = paste(All_otus_16S)), 
    by="rowname")

str(isolates_bact)
write.csv(arrange(isolates_bact, desc(`%IncMSE`)), "RF_bact_otu.csv")


# taxonomic proprotions
bact_bar_plot <-
  as.data.frame(
  table(gsub(".*-","",isolates_bact$Taxonomy))) %>%
  ggplot(aes(Var1, Freq)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title = element_text(angle = 0, size = 10, face = "bold"),
        axis.text.y = element_text(angle = 0, size = 7, hjust = 1, vjust = 0.5), 
        axis.text.x = element_text(angle = 0, size = 7 ,hjust = 0.5, vjust = 0.5),
        legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
        legend.title = element_text(size = 10, face = "bold"), 
        legend.text = element_text(size = 8), 
        plot.margin=unit(c(1.5,1.5,1.5,1.5),"pt")) +
  labs(title = "Bacterial OTUs", x="Taxon", y="Number of OTUs")

bact_bar_plot


# FIGURE SX - RF taxa barplot ------------------------------------------------------------------------------
ggarrange(
  ggarrange(fungi_bar_plot,
            ncol = 1,
            nrow = 2, 
            heights = c(1, 0.3),
            labels = "A"), 
  bact_bar_plot, 
  labels = c("","B"),
  ncol = 2,
  nrow = 1)



# **********************************************************************------------------------------------
# BETA DIVERSITY -------------------------------------------------------------------------------------------
str(meta_fungi_merged)
str(meta_bact_merged)

# match samples to those used for modeling: pot 340 raw 55, 124
meta_fungi_merged[rownames(meta_fungi_merged)==55,]
meta_fungi_merged_filt <-
  meta_fungi_merged[-55, ]
str(meta_fungi_merged_filt)

meta_bact_merged[rownames(meta_bact_merged)==124,]
meta_bact_merged_filt <-
  meta_bact_merged[-124, ]
str(meta_bact_merged_filt)


subset(meta_bact_merged, rowname%in%c("Amp489"))
subset(meta_bact_merged, rowname%in%c("Amp504"))


# Calculating fit of environemtnal variables ---------------------------------------------------------------
GenerateEnvFit <- function(df){
  envfit(df[,c("PCoA.1","PCoA.2")],
         df[, c("pH","P","K","Ca","Mg","OM","NO3")],
         perm = 9999) -> res_envfit
  df_envfit <- as.data.frame(scores(res_envfit, display = "vectors"))
  df_envfit <- cbind(df_envfit, Var = rownames(df_envfit))
  df_envfit$pvals <-res_envfit$vectors$pvals
  df_envfit <- subset(df_envfit, df_envfit$pvals<=0.05)
  df_envfit
  return(df_envfit)
}

set.seed(1)
GenerateEnvFit(meta_fungi_merged_filt) -> res_envfit_fungi
res_envfit_fungi

GenerateEnvFit(meta_bact_merged_filt) -> res_envfit_bact
res_envfit_bact

# plot ordination 
PlotPCOA <- function(df, envfit){
  pcoa <-
    ggplot(df,
         aes(x=PCoA.1, y=PCoA.2, color = Soil, shape = Genotype)) +
    geom_point(size = 1.5) +
    scale_colour_manual("Soil", values = Pal_soil) +
    scale_shape_manual(values = c(0,1,2,5,3,8)) +
    theme_classic() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.x = element_text(angle =0, size = 7, face = "bold", hjust = 0.5, vjust = 1), 
          axis.text.y = element_text(angle = 0, size = 7, face = "bold",hjust = 0.5, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8), legend.position="bottom") +
    guides(color = guide_legend(nrow = 2),shape = guide_legend(nrow = 2)) +
    #grids(linetype = "dashed") +
    geom_segment(data = envfit, inherit.aes = FALSE,
                           mapping = aes(x = 0, 
                                         xend = 1.05 * envfit$PCoA.1/2, 
                                         y = 0, 
                                         yend = 1.05 * envfit$PCoA.2/2), 
                           color = "black",
                           arrow = arrow(length = unit(0.02, "npc"))) +
    annotate("text", x = envfit$PCoA.1*0.555, y = envfit$PCoA.2*0.555,
             size = 3,
             label = c(
               "pH",
               expression(PO[4]^{"3-"}),
               expression(K^{"+"}),
               expression(Ca^{"2+"}),
               expression(Mg^{"2+"}),
               "OM",
               expression(NO[3]^{"-"})))
    # geom_text(data = envfit,inherit.aes = FALSE,
    #           mapping = aes(x = 1.05 * envfit[,1]/2, y = 1.05 * envfit[,2]/2, 
    #                         label = Var),
    #           size = 3) 
return(pcoa)
}
  
PlotPCOA(meta_fungi_merged_filt, res_envfit_fungi)

plot_pcoa_fungi <-
  PlotPCOA(meta_fungi_merged_filt, res_envfit_fungi) + 
  labs(title = "Fungi") 

plot_pcoa_fungi$layers[[2]]$aes_params$size <- 0.3
plot_pcoa_fungi$layers[[2]]$geom_params$arrow$length <- unit(6, units = "points")
plot_pcoa_fungi


plot_pcoa_bact <-
  PlotPCOA(meta_bact_merged_filt, res_envfit_bact) + 
  labs(title = "Bacteria") 

plot_pcoa_bact$layers[[2]]$aes_params$size <- 0.3
plot_pcoa_bact$layers[[2]]$geom_params$arrow$length <- unit(6, units = "points")
plot_pcoa_bact

plot_pcoa_bact + scale_x_reverse() + scale_y_reverse()


# Final FIGURE 2 -------------------------------------------------------------------------------------------
ggarrange(
  plot_pcoa_fungi,
  plot_pcoa_bact + scale_x_reverse() + scale_y_reverse(),
  ncol = 2,
  nrow = 1,
  common.legend = TRUE,
  legend = "none")


## PERMANOVA -----------------------------------------------------------------------------------------------
identical(rownames(alpha_df_fungi_filt), rownames(otu_fungi_new_filt))
str(otu_fungi_new_filt)
identical(rownames(alpha_df_bact_filt), rownames(otu_bact_new_filt))
str(otu_bact_new_filt)

adonis_fungi <- adonis(otu_fungi_new_filt ~ readNo + Genotype * Soil, 
                       alpha_df_fungi_filt, method = "bray", permutations=9999)
adonis_fungi

adonis_fungi2 <- adonis(otu_fungi_new_filt ~ readNo + Soil * Genotype, 
                       alpha_df_fungi_filt, method = "bray", permutations=9999)
adonis_fungi2

adonis_bact <- adonis(otu_bact_new_filt ~ readNo + Genotype * Soil, 
                      alpha_df_bact_filt, method = "bray", permutations=9999)
adonis_bact

adonis_bact2 <- adonis(otu_bact_new_filt ~ readNo + Soil * Genotype, 
                      alpha_df_bact_filt, method = "bray", permutations=9999)
adonis_bact2

# BETA DISPERSION ------------------------------------------------------------------------------------------

anova(permdisp_fungi_soil, permutations = 9999)
data.frame(multcompLetters(p.adjust(
  permutest(
    permdisp_fungi_soil,
    permutations = 9999,
    pairwise = T
  )$pairwise$observed,
  method = "BH"
))['Letters']) -> pair_fungi_soil


anova(permdisp_fungi_genotype, permutations = 9999)

anova(permdisp_bact_soil, permutations = 9999)
data.frame(multcompLetters(p.adjust(
  permutest(
    permdisp_bact_soil,
    permutations = 9999,
    pairwise = T
  )$pairwise$observed,
  method = "BH"
))['Letters']) -> pair_bact_soil


anova(permdisp_bact_genotype, permutations = 9999)


# plot multivariate dispersion 
PlotBetadisper <- function(betadisp, Var, my_labels, metadata){
  # creating a label and dataframe
  max(betadisp$distances + 0.1* max(betadisp$distances)) -> labels_y
  #labels_y = 0.9
  data.frame(betadisp$group, betadisp$distances) -> df
  colnames(df) <- c("Fact", "distance")
  #print(head(df))
  metadata <-
    left_join(tibble::rownames_to_column(df),
              tibble::rownames_to_column(metadata), by="rowname")
  head(metadata) %T>% print()
  # metadata$Fact <- dplyr::recode(metadata$Fact,
  #                                Alamo = "Alamo",
  #                                Blackwell = "Blackwell",
  #                                `Cave-in-rock` = "Cave-in-rock",
  #                                Kanlow ="Kanlow",
  #                                Shleter="Shelter",
  #                                Southlow="Southlow")
  # metadata$Genotype <- dplyr::recode(metadata$Genotype,
  #                                      Alamo = "Alamo",
  #                                      Blackwell = "Blackwell",
  #                                      `Cave-in-rock` = "Cave-in-rock",
  #                                      Kanlow ="Kanlow",
  #                                      Shleter="Shelter",
  #                                      Southlow="Southlow")
  #labels = c("Alamo", "Blackwell", "Cave-in-rock", "Kanlow", "Shelter","Southlow")
    # plotting
    betaplot <-
      ggplot(metadata, aes(x=get(Var), y=distance)) +
      geom_jitter(aes(shape = Genotype, color = Soil), alpha = 0.8, lwd=1) +
      geom_boxplot(outlier.colour="black", outlier.shape = 8, outlier.size = 1,
                   alpha=0.6, lwd = 0.4) +
      stat_summary(fun=mean, geom="point", shape=18, size=1.9, color="red", fill="red") +
      stat_summary(geom = 'text', label = my_labels, fun= max, aes(y = labels_y), size=3, color="black") +
      theme_classic() +
      scale_color_manual(values = Pal_soil) +
      scale_shape_manual(values = c(0,1,2,5,3,8), 
                         labels = c("Alamo"=expression(bold("Alamo")), 
                                    "Blackwell",
                                    "Cave-in-rock",
                                    "Kanlow"=expression(bold("Kanlow")),
                                    "Shelter",
                                    "Southlow")) +
      theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
            axis.title = element_text(angle = 0, size = 10, face = "bold"),
            axis.text.x = element_text(angle =45, size = 8, hjust = 1, vjust = 1), 
            axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5),
            legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
            legend.title = element_text(size = 10, face = "bold"), 
            legend.text = element_text(size = 8), legend.position="bottom")  +
      guides(color = guide_legend(nrow = 2),
             shape = guide_legend(nrow = 2)) 
    return(betaplot)
}



PlotBetadisper(permdisp_fungi_soil, "Soil", as.character(pair_fungi_soil$Letters), alpha_df_fungi_filt)
PlotBetadisper(permdisp_bact_soil, "Soil", as.character(pair_bact_soil$Letters), alpha_df_bact_filt)
PlotBetadisper(permdisp_fungi_genotype, "Genotype", NA, alpha_df_fungi_filt) +
  scale_x_discrete(labels=c(Alamo = "Alamo",
                            Blackwell = "Blackwell",
                            `Cave-in-rock` = "Cave-in-rock",
                            Kanlow ="Kanlow",
                            Shlelter="Shelter",
                            Southlow="Southlow"))

PlotBetadisper(permdisp_bact_genotype, "Genotype", NA, alpha_df_bact_filt)


# *** FIGURE 2B - Beta Dispersion --------------------------------------------------------------------------
ggarrange(PlotBetadisper(permdisp_fungi_soil, "Soil", as.character(pair_fungi_soil$Letters), alpha_df_fungi_filt) + 
            labs(title = "Soil", y ="Distance to centroid", x= NULL),
          PlotBetadisper(permdisp_fungi_genotype, "Genotype", NA, alpha_df_fungi_filt) + 
            labs(title = "Genotype", y ="Distance to centroid", x= NULL)+
            scale_x_discrete(labels=c(Alamo = "Alamo",
                                      Blackwell = "Blackwell",
                                      `Cave-in-rock` = "Cave-in-rock",
                                      Kanlow ="Kanlow",
                                      Shlelter="Shelter",
                                      Southlow="Southlow")),
          PlotBetadisper(permdisp_bact_soil, "Soil", as.character(pair_bact_soil$Letters), alpha_df_bact_filt) + 
            labs(title = "Soil", y ="Distance to centroid", x= NULL), 
          PlotBetadisper(permdisp_bact_genotype, "Genotype", NA, alpha_df_bact_filt) +
            labs(title = "Genotype", y ="Distance to centroid", x= NULL) +
            scale_x_discrete(labels=c(Alamo = "Alamo",
                                      Blackwell = "Blackwell",
                                      `Cave-in-rock` = "Cave-in-rock",
                                      Kanlow ="Kanlow",
                                      Shlelter="Shelter",
                                      Southlow="Southlow")),
          widths = c(1, 1.25, 1, 1.25),
          labels = c("C","D","E","F"),
          align = "hv",
          ncol = 4, 
          nrow = 1,
          common.legend = TRUE, 
          legend = c("bottom")) -> betadisp_plot

betadisp_plot

# Calculating axis % of variation in cmdscale --------------------------------------------------------------
var_axes_fungi <-
  round(pcoa_its$eig*100/sum(pcoa_its$eig),1)
var_axes_bact <-
  round(pcoa_16s$eig*100/sum(pcoa_16s$eig),1)


# *** FIGURE 2 - complete ----------------------------------------------------------------------------------
ggarrange(
  ggarrange(
    plot_pcoa_fungi +
      annotate("text", -Inf, Inf, 
               label = expression(paste("Read No. ", italic(R) ^ 2,"= 7.7%***"),
                                  parse=TRUE), size = 2.5, hjust = -0.05, vjust = 1.1) +
      annotate("text", -Inf, Inf, 
               label = expression(paste("Genotype ", italic(R) ^ 2,"= 2.8%***"),
                                  parse=TRUE), size = 2.5, hjust = -0.05, vjust = 2) +
      annotate("text", -Inf, Inf, 
               label = expression(paste("Soil ", italic(R) ^ 2,"= 43.4%***"),
                                  parse=TRUE), size = 2.5, hjust = -0.06, vjust = 3.2) +
      labs(x=as.expression(paste("PCoA.1 (",var_axes_fungi[1],"%)"), parse=TRUE),
           y=as.expression(paste("PCoA.2 (",var_axes_fungi[2],"%)"), parse=TRUE)),
    plot_pcoa_bact  +
      annotate("text", Inf, -Inf, 
               label = expression(paste("Read No. ", italic(R) ^ 2,"= 12.9%***"),
                                  parse=TRUE), size = 2.5, hjust = -0.05, vjust = 1.1) +
      annotate("text", Inf, -Inf, 
               label = expression(paste("Genotype ", italic(R) ^ 2,"= 3.4%***"),
                                  parse=TRUE), size = 2.5, hjust = -0.05, vjust = 2) +
      annotate("text", Inf, -Inf, 
               label = expression(paste("Soil ", italic(R) ^ 2,"= 32.4%***"),
                                  parse=TRUE), size = 2.5, hjust = -0.06, vjust = 3.2) + 
      labs(x=as.expression(paste("PCoA.1 (",var_axes_bact[1],"%)"), parse=TRUE),
           y=as.expression(paste("PCoA.2 (",var_axes_bact[2],"%)"), parse=TRUE)) +
    scale_x_reverse() + 
    scale_y_reverse(),
    ncol = 2,
    nrow = 1,
    common.legend = TRUE,
    legend = "none",
    labels = c("A","B")),
  betadisp_plot,
  nrow = 2, 
  heights = c(1, 1)) -> Fig2_beta_div

Fig2_beta_div

# *******************************************************************---------------------------------------
# TAXONOMIC COMPOSITION ------------------------------------------------------------------------------------

AbundData <- function(physeq, rank){
  tax_table(physeq)[tax_table(physeq)==""]<- NA
  tax_table(physeq)[tax_table(physeq)=="unidentified"]<- NA
  tax_table(physeq)[is.na(tax_table(physeq))]<-"Unclassified"
  tax_table(physeq) <- tax_table(physeq)[, !(colnames(tax_table(physeq)) %in% c("OTU_ID"))]
  physeq_phylum <- tax_glom(physeq, rank)
  otu_physeq <- taxa_sums(physeq_phylum)/sum(taxa_sums(physeq_phylum))*100
  tax_physeq <- as(tax_table(physeq_phylum), "matrix")
  tax_physeq <- as.data.frame(tax_physeq)
  #tax_physeq <- tax_physeq[c(2)]
  tax_physeq$abundance <- as.vector(otu_physeq)
  tax_physeq <- tax_physeq[order(tax_physeq$abundance, decreasing = TRUE),] 
  return(tax_physeq)
}


AbundData(physeq_fungi_new, "Phylum")
AbundData(physeq_bact_new, "Phylum")

str(physeq_fungi_new@sam_data)
physeq_fungi_new@sam_data
tail(physeq_fungi_new@sam_data)
table(physeq_fungi_new@sam_data$Ecotype)


#*********************************************************************************--------------------------
# Models - From Ming-Yi 
# # Colinear removal + MLR + Stepwise Regression  ###
# library(statsr)
# library(MASS)
# library(dplyr)
# library(ggplot2)
# library(BAS)
# library(olsrr)
# library(corpcor)
# 
# 
# full.model <- lm(Severity16_Raw ~ pH+OM+N+C+P+S+Ca+Mg+K+Al+Cu+Fe+Mn+Mo+Na+Zn , data = soil.sever)
# summary(full.model)
# 
# step.model <- stepAIC(full.model, direction = "backward", trace = FALSE)
# summary(step.model)
# 
# # remove colinear #
# variables=soil.sever[(1:18),(11:27)]
# 
# # run vif_func at: https://gist.github.com/fawda123/4717702 #
# vif_func(in_frame=variables,thresh=10,trace=T)

# rerun the model #
# full.model <- lm(Severity16_Raw ~ pH+OM+K+C+P+S++Cu+Fe+Mn+Mo+Na+Zn , data = soil.sever)
# step.model <- stepAIC(full.model, direction = "backward", trace = FALSE)
# summary(step.model)

# 
# data.frame(richness = c(10, 12, 23, 25, 2, 45,8, 5, 9, 10,23,44), 
#            Calcium = c(4,4,5,5,10,10,12,12,30,30, 23,23),
#            Genotype = c("Alamo", "Alamo", "Cave-in-rock","Cave-in-rock","Kanlow","Kanlow",
#                         "Alamo", "Alamo", "Cave-in-rock","Cave-in-rock","Kanlow","Kanlow"),
#            Soil = c("Hancock", "Hancock", "Hancock", "Hancock", "Hancock", "Hancock",
#                     "Lake city", "Lake city","Lake city", "Lake city","Lake city", "Lake city"))


# Additional analyses --------------------------------------------------------------------------------------

# FUNGI
# now, we separate samples as high-beta (hancock/lake city) and low beta (Lux/rhineland)
physeq_fungi_new_H  <- 
  subset_samples(physeq_fungi_new, Soil_location == "Hancock" | Soil_location == "Lake City")
physeq_fungi_new_H <- 
  prune_taxa(taxa_sums(physeq_fungi_new_H) > 0, physeq_fungi_new_H) 


physeq_fungi_new_L <- 
  subset_samples(physeq_fungi_new, Soil_location == "Rhineland" | Soil_location == "Lux Arbor")
physeq_fungi_new_L <- 
  prune_taxa(taxa_sums(physeq_fungi_new_L) > 0, physeq_fungi_new_L) 


## high diversity samples, simple model ##
metadata_fungi_H <- as(sample_data(physeq_fungi_new_H), "data.frame")
adonis(phyloseq::distance(t(otu_table(physeq_fungi_new_H)), method="bray") ~ LibrarySize + Genotype * Soil_location, permutations=999, data = metadata_fungi_H)
adonis(phyloseq::distance(t(otu_table(physeq_fungi_new_H)), method="bray") ~ LibrarySize + Soil_location * Genotype, permutations=999, data = metadata_fungi_H)
adonis(phyloseq::distance(t(otu_table(physeq_fungi_new_H)), method="bray") ~ LibrarySize + Genotype, 
       strata = metadata_fungi_H$Soil_location, permutations=999, data = metadata_fungi_H)

anova(
  betadisper(phyloseq::distance(t(otu_table(physeq_fungi_new_H)), method="bray"), metadata_fungi_H$Genotype), 
  permutations = 999)


## low diversity soil, simple model ##
metadata_fungi_L <- as(sample_data(physeq_fungi_new_L), "data.frame")
adonis(phyloseq::distance(t(otu_table(physeq_fungi_new_L)), method="bray") ~ LibrarySize + Genotype * Soil_location, permutations=999, data = metadata_fungi_L)
adonis(phyloseq::distance(t(otu_table(physeq_fungi_new_L)), method="bray") ~ LibrarySize + Soil_location * Genotype, permutations=999, data = metadata_fungi_L)
adonis(phyloseq::distance(t(otu_table(physeq_fungi_new_L)), method="bray") ~ LibrarySize + Genotype, 
       strata = metadata_fungi_L$Soil_location, permutations=999, data = metadata_fungi_L)

anova(
  betadisper(phyloseq::distance(t(otu_table(physeq_fungi_new_L)), method="bray"), metadata_fungi_L$Genotype), 
  permutations = 999)



# BACTERIA
# now, we separate samples as high-beta (hancock/lake city) and low beta (Lux/rhineland)
physeq_bact_new_H  <- 
  subset_samples(physeq_bact_new, Soil_location == "Hancock" | Soil_location == "Lake City")
physeq_bact_new_H <- 
  prune_taxa(taxa_sums(physeq_bact_new_H) > 0, physeq_bact_new_H) 


physeq_bact_new_L <- 
  subset_samples(physeq_bact_new, Soil_location == "Rhineland" | Soil_location == "Lux Arbor")
physeq_bact_new_L <- 
  prune_taxa(taxa_sums(physeq_bact_new_L) > 0, physeq_bact_new_L) 


## high diversity samples, simple model ##
metadata_bact_H <- as(sample_data(physeq_bact_new_H), "data.frame")
adonis(phyloseq::distance(t(otu_table(physeq_bact_new_H)), method="bray") ~ LibrarySize + Genotype * Soil_location, permutations=999, data = metadata_bact_H)
adonis(phyloseq::distance(t(otu_table(physeq_bact_new_H)), method="bray") ~ LibrarySize + Soil_location * Genotype, permutations=999, data = metadata_bact_H)
adonis(phyloseq::distance(t(otu_table(physeq_bact_new_H)), method="bray") ~ LibrarySize + Genotype, 
       strata = metadata_bact_H$Soil_location, permutations=999, data = metadata_bact_H)


anova(
  betadisper(phyloseq::distance(t(otu_table(physeq_bact_new_H)), method="bray"), metadata_bact_H$Genotype), 
  permutations = 999)


## low diversity soil, simple model ##
metadata_bact_L <- as(sample_data(physeq_bact_new_L), "data.frame")
adonis(phyloseq::distance(t(otu_table(physeq_bact_new_L)), method="bray") ~ LibrarySize + Genotype * Soil_location, permutations=999, data = metadata_bact_L)
adonis(phyloseq::distance(t(otu_table(physeq_bact_new_L)), method="bray") ~ LibrarySize + Soil_location * Genotype, permutations=999, data = metadata_bact_L)
adonis(phyloseq::distance(t(otu_table(physeq_bact_new_L)), method="bray") ~ LibrarySize + Genotype, 
       strata = metadata_bact_L$Soil_location, permutations=999, data = metadata_bact_L)

anova(
  betadisper(phyloseq::distance(t(otu_table(physeq_bact_new_L)), method="bray"), metadata_bact_L$Genotype), 
  permutations = 999)

# Plot diversity in soils soils ----------------------------------------------------------------------------

# Extract soils

physeq_fungi_new_soil  <- 
  subset_samples(physeq_fungi_new, Genotype == "Control")
physeq_fungi_new_soil <- 
  prune_taxa(taxa_sums(physeq_fungi_new_soil) > 0, physeq_fungi_new_soil) 
physeq_fungi_new_soil

physeq_fungi_new_soil@sam_data


physeq_bact_new_soil  <- 
  subset_samples(physeq_bact_new, Genotype == "Control")
physeq_bact_new_soil <- 
  prune_taxa(taxa_sums(physeq_bact_new_soil) > 0, physeq_bact_new_soil) 
physeq_bact_new_soil

physeq_bact_new_soil@sam_data

# Create dataframes for plotting Alpha diversity in soils
alpha_df_fungi_soil <-
  MakeDf(physeq_fungi_new_soil, 10000) 
alpha_df_fungi_soil

alpha_df_bact_soil <-
  MakeDf(physeq_bact_new_soil, 10000) 
alpha_df_bact_soil


CompSamplSoil <- function(dataframe, formula){
  require(multcompView)
  compare_means(formula, data = dataframe,method = "wilcox.test",
                p.adjust.method = "BH") -> test_CC 
  test_CC <- as.data.frame(test_CC)[,c(2,3,5)] # to change form p to p.adj do 4 to 5
  test_CC2 <- data.frame(test_CC[,2], test_CC[,1], test_CC[,3])
  colnames(test_CC2) <- c("group1", "group2", "p.adj") # change p to p.adj
  rbind(test_CC, test_CC2) -> test_all
  as.dist(xtabs(test_all[, 3] ~ (test_all[, 2] + test_all[, 1])), diag = TRUE) -> dist_CC
  data.frame(multcompLetters(dist_CC)['Letters']) -> res_CC
  res_CC$sample <- rownames(res_CC)
  res_CC %>% slice(match(c("Hancock", "Lake City","Rhineland","Lux Arbor"), sample)) -> res_CC
  return(res_CC)
}


CompSamplSoil(alpha_df_fungi_soil, formula(richness ~ Soil))
CompSamplSoil(alpha_df_bact_soil, formula(richness ~ Soil))


PlotAlphaDivSoil <- function(df, Var){
  plot_div <-
    ggplot(df, aes(x=Soil, y=get(Var))) +
    geom_boxplot(outlier.shape = 1, outlier.size = 1, outlier.stroke = 1,
                 position = position_dodge(preserve = "single"), alpha=0.6, lwd = 0.5) +
    stat_summary(fun=mean, geom="point", shape=18, size=1.9, color="red", fill="red") +
    ylim(0, NA) +
    theme_classic() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.x = element_text(angle =45, size = 8, hjust = 1, vjust = 1), 
          axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8)) +
    labs(x=NULL)
  return(plot_div)
}


PlotAlphaDivSoil(alpha_df_fungi_soil, "richness") +
  stat_summary(geom = 'text', 
               label = CompSamplSoil(alpha_df_fungi_soil, formula(richness ~ Soil))$Letters, 
               fun = max, aes(y = 1000), size = 3.5, color = "black") + 
  labs(title = "Fungi Richness Soil")


ggarrange(
  PlotAlphaDivSoil(alpha_df_fungi_soil, "richness") +
    stat_summary(geom = 'text', 
                 label = CompSamplSoil(alpha_df_fungi_soil, formula(richness ~ Soil))$Letters, 
                 fun = max, aes(y = 1000), size = 3.5, color = "black") + 
    labs(title = "Fungi soil", y="OTU richness"),
  PlotAlphaDivSoil(alpha_df_fungi, "richness") +
    stat_summary(geom = 'text', 
                 label = CompSamplSoil(alpha_df_fungi, formula(richness ~ Soil))$Letters, 
                 fun = max, aes(y = 1000), size = 3.5, color = "black") + 
    labs(title = "Fungi roots",y="OTU richness"),
  PlotAlphaDivSoil(alpha_df_bact_soil, "richness") +
    stat_summary(geom = 'text', 
                 label = CompSamplSoil(alpha_df_bact_soil, formula(richness ~ Soil))$Letters, 
                 fun = max, aes(y = 5000), size = 3.5, color = "black") + 
    labs(title = "Bacteria soil",y="OTU richness"),
  PlotAlphaDivSoil(alpha_df_bact, "richness") +
    stat_summary(geom = 'text', 
                 label = CompSamplSoil(alpha_df_bact, formula(richness ~ Soil))$Letters, 
                 fun = max, aes(y = 5000), size = 3.5, color = "black") + 
    labs(title = "Bacteria roots", y="OTU richness"),
ncol = 2,
nrow = 2)



# NEUTRAL MODELS ------------------------------------------------------------------------------------------
library(tyRa)
require(minpack.lm)
require(Hmisc)
require(stats4)

neutral_model_fungi <- 
  tyRa::fit_sncm(spp = t(otu_table(physeq_fungi_new)), pool=NULL, taxon=data.frame(tax_table(physeq_fungi_new)))

plot_sncm_fit(neutral_model_fungi, fill = NULL, title = "Model Fit") +
  theme_classic()

neutral_model_bact <- 
  tyRa::fit_sncm(spp = t(otu_table(physeq_bact_new)), pool=NULL, taxon=data.frame(tax_table(physeq_bact_new)))

plot_sncm_fit(neutral_model_bact, fill = NULL, title = "Model Fit") +
  theme_classic()


PlotOTU <- function(rf_model, taxa, NeutMod, Var, df_name){
  require(tibble)
  imp_RF <- as.data.frame(rf_model$importance)
  imp_RF$features <- rownames(imp_RF)
  imp_RF <- arrange(imp_RF, desc(Var))
  rownames(imp_RF) <- imp_RF$features
  # adding marker taxon info
  taxa[rownames(taxa)%in%rownames(imp_RF), ] -> taxa_RF
  identical(rownames(taxa_RF), rownames(imp_RF))
  order_taxa <- match(rownames(taxa_RF), rownames(imp_RF))
  imp_RF <- imp_RF[order_taxa,]
  imp_RF$Taxonomy <- taxa_RF$Taxon
  imp_RF <- imp_RF[order(imp_RF[,Var], decreasing = TRUE),] 
  imp_RF <- left_join(tibble::rownames_to_column(imp_RF),
                      tibble::rownames_to_column(
                        taxa[rownames(imp_RF), c(9:11, 16)]))
  imp_RF$Taxonomy <-
    gsub("FOTU_", "", imp_RF$Taxonomy)
  imp_RF$Taxonomy <-
    gsub("OTU_", "",imp_RF$Taxonomy)
  imp_RF$glbrc <-
    ifelse(is.na(imp_RF$Isolate), "no", "yes")
  # adding neutral model results
  df_neutral <-
    as.data.frame(NeutMod[2])
  df_neutral$rowname <- rownames(df_neutral)
  new_df <-
    df_neutral %>%
    dplyr::select("predictions.fit_class", "rowname")
  new_df <-
  left_join(imp_RF, new_df, by= "rowname")
  new_df %T>% print()
 # saving intermediate df to R environemnt, pick the right name
   new_df %T>% 
    print() %T>%
    #assign(paste(df_name, Var1, Var2, sep = ""),., envir = .GlobalEnv) %>% # saving the plot 
    assign(df_name, ., envir = .GlobalEnv) 
  # plotting
  ggplot(data=imp_RF) + 
    geom_bar(aes(x= reorder(Taxonomy, -get(Var)),
                 y= get(Var),
                 color= glbrc, fill=glbrc), stat="identity") +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0)) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title = element_text(angle = 0, size = 10, face = "bold"),
          axis.text.y = element_text(angle =0, size = 8, hjust = 0.5, vjust = 1), 
          axis.text.x = element_text(angle = 90, size = 7 ,hjust = 1, vjust = 0.5),
          legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
          legend.title = element_text(size = 10, face = "bold"), 
          legend.text = element_text(size = 8), 
          plot.margin=unit(c(1.5,1.5,1.5,1.5),"pt"),
          legend.position = "none") -> plot_importance
  return(plot_importance)
}


PlotOTU(RF_fungi_biom, taxa_fungi_new, neutral_model_fungi, "%IncMSE", "rf_taxa_fungi")
head(rf_taxa_fungi)
table(rf_taxa_fungi$predictions.fit_class)

PlotOTU(RF_bact_biom, taxa_bact_new, neutral_model_bact, "%IncMSE", "rf_taxa_bact")
head(rf_taxa_bact)
table(rf_taxa_bact$predictions.fit_class)


write.csv(rf_taxa_fungi, "rf_taxa_fungi.csv")
write.csv(rf_taxa_bact, "rf_taxa_bact.csv")


# Spearmna corraltions -------------------------------------------------------------------------------------
library(corrplot)
head(df_bact_RF_2)
head(df_fungi_RF_2)

soil_data_fungi_norm <-
  data.frame(
    apply(df_fungi_RF_2[, -c(14,15,16)], 2, function(x) scale(x, center = TRUE, scale = TRUE))
  )
soil_data_fungi_norm


cor_fungi_soil <- cor(soil_data_fungi_norm, method = "spearman")
corrplot(cor_fungi_soil, 
         method = 'number', 
         type = 'lower', 
         title="Spearman correlations Fungi",
         diag = FALSE,
         mar=c(0,0,1,0)) # colorful number

soil_data_bact_norm <-
  data.frame(
    apply(df_bact_RF_2[, -c(6,15,16,17)], 2, function(x) scale(x, center = TRUE, scale = TRUE))
  )
soil_data_bact_norm


cor_bact_soil <- cor(soil_data_bact_norm, method = "spearman")
corrplot(cor_bact_soil, 
         method = 'number', 
         type = 'upper', 
         title="Spearman correlations Bacteria",
         diag = FALSE,
         mar=c(0,0,1,0)) # colorful number


soil_data <- meta_fungi_merged_filt[, 37:43]
soil_data$Biomass <- meta_fungi_filt$aerial_part_dry_weight_16wks_grams 
soil_data

soil_data_norm <-
  data.frame(
    apply(soil_data, 2, function(x) scale(x, center = TRUE, scale = TRUE))
  )
soil_data_norm

hist(soil_data_norm$Biomass)
hist(soil_data_norm$P)
hist(soil_data_norm$NO3)


cor_soil <- cor(soil_data_norm, method = "spearman")
corrplot(cor_soil, 
         method = 'number', 
         type = 'upper', 
         title="Spearman correlations",
         diag = FALSE,
         mar=c(0,0,1,0)) # colorful number


# OTU that correlate to plant biomass ----------------------------------------------------------------------
library(corrplot)
features_fungi <- c("FOTU_227","FOTU_57","FOTU_772", "FOTU_12")
features_bact <- c("BOTU_700", "BOTU_291", "BOTU_64", "BOTU_727")

physeq_fungi_new

# function to remove bad taxa
removeSamples = function(physeq, goodSamples) {
    allSamples = sample_names(physeq)
    mySamples <- allSamples[(allSamples %in% goodSamples)]
    return(prune_samples(mySamples, physeq))
  }
  
physeq_fungi_filt <-
  removeSamples(physeq_fungi_new, meta_fungi_merged_filt$rowname)
physeq_fungi_filt

rf_features_fungi <-
  subset_taxa(physeq_fungi_filt, OTU_ID%in%features_fungi)
rf_features_fungi@tax_table
rf_features_fungi@otu_table

df_otu_cor <-
  as.data.frame(t(as(rf_features_fungi@otu_table, "matrix")))
identical(rownames(df_otu_cor),rownames(rf_features_fungi@sam_data))

df_otu_cor$Biomass <- rf_features_fungi@sam_data$aerial_part_dry_weight_16wks_grams
df_otu_cor

df_otu_cor_norm <-
  data.frame(
    apply(df_otu_cor, 2, function(x) scale(x, center = TRUE, scale = TRUE))
  )
df_otu_cor_norm

cor_otu_cor_norm <- cor(df_otu_cor_norm, method = "spearman")
cor_otu_cor_norm

corrplot(cor_otu_cor_norm, 
         method = 'number', 
         type = 'upper', 
         title="Spearman correlations",
         diag = FALSE,
         mar=c(0,0,1,0)) # colorful number

# SEQUENCING RESULTS ----------------------------------------------------------------------------------------
physeq_fungi_new %>% 
  subset_samples(Ecotype%in%c("Lowland", "Upland"))


mean(sample_sums(physeq_fungi_filt))
sd(sample_sums(physeq_fungi_filt))
min(sample_sums(physeq_fungi_filt))
max(sample_sums(physeq_fungi_filt))

mean(sample_sums(physeq_fungi_new_soil))
sd(sample_sums(physeq_fungi_new_soil))
min(sample_sums(physeq_fungi_new_soil))
max(sample_sums(physeq_fungi_new_soil))


mean(sample_sums(physeq_bact_new_soil))
sd(sample_sums(physeq_bact_new_soil))

mean(sample_sums(
  physeq_bact_new %>% 
    subset_samples(Ecotype%in%c("Lowland", "Upland"))
))

sd(sample_sums(
  physeq_bact_new %>% 
    subset_samples(Ecotype%in%c("Lowland", "Upland"))
))


# subsetting data for Ming-Yi -------------------------------------------------------------------------------
rownames(as.data.frame(RF_fungi_biom$importance))

MSE_features_fungi <-
  subset_taxa(physeq_fungi_new, OTU_ID%in%rownames(as.data.frame(RF_fungi_biom$importance)))
MSE_features_fungi
MSE_features_fungi@tax_table

rownames(as.data.frame(RF_bact_biom$importance))

MSE_features_bact <-
  subset_taxa(physeq_bact_new, OTU_ID%in%rownames(as.data.frame(RF_bact_biom$importance)))
MSE_features_bact
MSE_features_bact@tax_table


MSE_features_fungi@sam_data$Ecotype
MSE_features_fungi <-
  subset_samples(MSE_features_fungi, Ecotype%in%c("Upland", "Lowland"))
otu_table(MSE_features_fungi) <-
  otu_table(MSE_features_fungi)[which(rowSums(otu_table(MSE_features_fungi)) > 0),]

MSE_features_bact@sam_data$Ecotype
MSE_features_bact <-
  subset_samples(MSE_features_bact, Ecotype%in%c("Upland", "Lowland"))
otu_table(MSE_features_bact) <-
  otu_table(MSE_features_bact)[which(rowSums(otu_table(MSE_features_bact)) > 0),]

# Save an object to a file
saveRDS(MSE_features_fungi, file = "MSE_features_fungi.rds")
saveRDS(MSE_features_bact, file = "MSE_features_bact.rds")
# Restore the object
readRDS(file = "MSE_features_fungi.rds")
readRDS(file = "MSE_features_bact.rds")


# LINEAR MODELS --------------------------------------------------------------------------------------------

## read data ##
fun <- readRDS("MSE_features_fungi.rds")
bac <- readRDS("MSE_features_bact.rds")


#Fungi
fun.fac <- t(as.data.frame(fun@otu_table))
fun.bio <- as.data.frame(fun@sam_data)

df<-as.data.frame(cbind(fun.fac, fun.bio$aerial_part_dry_weight_16wks_grams))
colnames(df)[55]<- "biomass"

full.model <- lm(biomass ~ . , data = df)
summary(full.model)

step.model <- stepAIC(full.model, direction = "backward", trace = FALSE)
summary(step.model)


#Bacteria
bac.fac <- t(as.data.frame(bac@otu_table))
bac.bio <- as.data.frame(bac@sam_data)

df.bac<-as.data.frame(cbind(bac.fac, bac.bio$aerial_part_dry_weight_16wks_grams))
colnames(df.bac)[53]<- "biomass"

full.model.bac <- lm(biomass ~ . , data = df.bac)
summary(full.model.bac)

step.model.bac <- stepAIC(full.model.bac, direction = "backward", trace = FALSE)
summary(step.model.bac)

#### Mantel Test ####

### Correlation between microbial composition and plant overall phenotype ###
# read file for correlation #
#Mantel 16S
library(fastDummies)

otus.bray.s # 16S microbe matrices
otus.bray.s.f # ITS microbe matrices

d.16s.pheno <- m.summer[,(44:52)] # 16S plant phenotype data
d.ITS.pheno <- m.summer.f[,(44:52)] # ITS plant phenotype data

d.16s.pheno.1 <- d.16s.pheno[,-(4:5)]
d.ITS.pheno.1 <- d.ITS.pheno[,-(4:5)]

d.16s.pheno.1 <- dummy_columns(d.16s.pheno.1)
d.ITS.pheno.1 <- dummy_columns(d.ITS.pheno.1)

d.16s.pheno.2 <- d.16s.pheno.1[,-(4:5)]
d.ITS.pheno.2 <- d.ITS.pheno.1[,-(4:5)]

m.16s.pheno <- vegdist(d.16s.pheno.2,"bray") # 16S plant phenotype matric
m.ITS.pheno <- vegdist(d.ITS.pheno.2,"bray") # ITS plant phenotype matric

set.seed(15110)
Mantel.plant.16S <- mantel(otus.bray.s, m.16s.pheno, method="spear", permutations=999, strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))
Mantel.plant.ITS <- mantel(otus.bray.s.f, m.ITS.pheno , method="spear", permutations=999, strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))
Mantel.Plant.chem<- mantel(m.16s.pheno,m.16s.chem , method="spear", permutations=999, strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))

# Mantel ITS
otus.fac_filt <- subset(otus.fac, !Pot%in%c("326", "93", "284", "117", "256", "117", "128","340"))
otus.summer <- subset(otus.fac_filt, Genotype != "Control") # run if extra filtering
o.summer <-otus.summer[,(1:12072)]
m.summer <-otus.summer[,(12073:12135)]
otus.bray.s <- vegdist(o.summer) # 16S microbe matrices without 340

otus.fac.f_filt <- subset(otus.fac.f, !Pot%in%c("184", "109", "89", "87", "303", "198","340"))
otus.summer.f <- subset(otus.fac.f_filt, Genotype != "Control") # run if extra filtering
o.summer.f <-otus.summer.f[,(1:5837)]
m.summer.f <-otus.summer.f[,(5838:5901)]
otus.bray.s.f <- vegdist(o.summer.f) # ITS microbe matrices without 340

d.16s.chem <- m.summer[,(55:61)] # 16S soil chem data
d.ITS.chem <- m.summer.f[,(55:61)] # ITS soil chem data

m.16s.chem <- vegdist(d.16s.chem,"bray") # 16S soil chem matric
m.ITS.chem <- vegdist(d.ITS.chem,"bray") # ITS soil chem matric

set.seed(15110)
Mantel.chem.16S <- mantel(otus.bray.s, m.16s.chem, method="spear", permutations=999, strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))
Mantel.chem.ITS <- mantel(otus.bray.s.f, m.ITS.chem , method="spear", permutations=999, strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))



