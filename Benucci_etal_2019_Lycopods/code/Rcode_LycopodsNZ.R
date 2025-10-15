# ****************DATA ANALYSIS ********************************************************************************************* -----
# Project name: Fungal and Bacterial Associations in Early Diverging Lycopodiaceae
# Manuscript:   Evidence for Co-Evolutionary History of Early Diverging Lycopodiaceae Plants With Fungi But Not Bacteria
# Authors:      Gian Maria Niccolò Benucci, Delaney Burnard, Lara D. Shepherd, Gregory Bonito, Andrew Munkacsi
# Affiliation:  Michigan State University, Victoria University of Wellington, Museum of New Zealand Te Papa Tongarewa
# Journal:      Forntiers in Microbiology
# Date:         October 3, 2019
# *************************************************************************************************************************** -----

# ___________WORKING ENVIRONMENT SETUP _____________ ---------------------------------------------------
options(scipen = 999) #to use decimals
options(max.print=100000000) # to print more lines on the display
options(verbose=TRUE)
#rm(list = ls()) # use with caution

# >>> IMPORTING DATASETS -------------------------------------------------------------------------------
library(phyloseq); packageVersion("phyloseq")
library(ggplot2)
library(Biostrings)
library(ape)

# IMPORTING DATASETS -----------------------------------------------------------------------------------
getwd(); dir(); ls() 

# >>> COLOR PALETTES --------------------------------------------------------------------------------------

palette_fungi = c("#781156","#A51876","#D21E96","#E43FAD","#EA6CC0","#F098D3","#114578","#185EA5",
                  "#1E78D2","#3F91E4","#6CABEA","#98C4F0","#117878","#18A5A5","#3FE4E4","#6CEAEA",
                  "#98F0F0", "#117845","#18A55E","#1ED278","#3FE491","#6CEAAB","#98F0C4","#787811",
                  "#A5A518","#D2D21E","#E4E43F","#EAEA6C","#F0F098","#F7F7C5","#784511","#A55E18",
                  "#D2781E","#E4913F","#EAAB6C","#F0C498","#781122","#A5182F","#D21E2C")

palette_bact = c("#781156","#A51876","#D21E96","#E43FAD","#EA6CC0","#F098D3","#114578","#185EA5",
                 "#1E78D2","#3F91E4","#6CABEA","#98C4F0","#117878","#18A5A5","#3FE4E4","#6CEAEA",
                 "#98F0F0", "#117845","#18A55E","#1ED278","#3FE491","#6CEAAB","#98F0C4","#787811",
                 "#A5A518","#D2D21E","#E4E43F","#EAEA6C","#D21E2C")


# Import fungal ITS data ------------------------------------------------------------------------------
otus_ITS_uparse_R1 <- read.delim("analysis_ITS/otu_table_ITS_UPARSE_R1.txt",row.names=1) 
otus_phy_ITS_uparse_R1 <-otu_table(otus_ITS_uparse_R1, taxa_are_rows = TRUE)

metadata_ITS_uparse_R1 <-read.delim("analysis_ITS/mapping_ITS.txt", row.names=1, header=TRUE, sep="\t")
metadata_phy_ITS_uparse_R1 <-sample_data(metadata_ITS_uparse_R1)

taxonomy_ITS_uparse_R1_cons <-read.delim("analysis_ITS/outputs_UPARSE_R1/consensus_taxonomy.txt", header=TRUE, row.names=1)
head(taxonomy_ITS_uparse_R1_cons)
taxonomy_ITS_uparse_R1_RDP <-read.delim("analysis_ITS/outputs_UPARSE_R1/otu_taxonomy_rdp_final.txt", header=TRUE, row.names=1)
# Taxonomy_ITS_uparse_R1_RDP <- subset(taxonomy_ITS_uparse_R1_RDP, 
# select=c(Kingdom, Phylum, Class, Order, Family, Genus, Species))
identical(rownames(taxonomy_ITS_uparse_R1_cons), rownames(taxonomy_ITS_uparse_R1_RDP))
otu_order <- match(rownames(taxonomy_ITS_uparse_R1_RDP), rownames(taxonomy_ITS_uparse_R1_cons))
taxonomy_ITS_uparse_R1_cons <- taxonomy_ITS_uparse_R1_cons[otu_order,]
colnames(taxonomy_ITS_uparse_R1_RDP)[2] <- "RDP"
taxonomy_ITS_uparse_R1_cons$RDP <- taxonomy_ITS_uparse_R1_RDP$RDP
head(taxonomy_ITS_uparse_R1_cons)
taxonomy_phy_ITS_uparse_R1_cons <- tax_table(as.matrix(taxonomy_ITS_uparse_R1_cons))

otus_seq_ITS_uparse_R1 <- readDNAStringSet("analysis_ITS/otus_R1.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

physeq_obj_ITS_uparse_R1 <- phyloseq(otus_phy_ITS_uparse_R1, 
                                     metadata_phy_ITS_uparse_R1,
                                     taxonomy_phy_ITS_uparse_R1_cons,
                                     otus_seq_ITS_uparse_R1) 

physeq_obj_ITS_uparse_R1
str(physeq_obj_ITS_uparse_R1)
tax_table(physeq_obj_ITS_uparse_R1)
sample_data(physeq_obj_ITS_uparse_R1)

tax_table(physeq_obj_ITS_uparse_R1)[tax_table(physeq_obj_ITS_uparse_R1)==""]<- NA
tax_table(physeq_obj_ITS_uparse_R1)[is.na(tax_table(physeq_obj_ITS_uparse_R1))]<-"Unclassified"

# removing non fungal taxa -----------------------------------------------------------------------------
sort(unique(as.data.frame(tax_table(physeq_obj_ITS_uparse_R1))$Kingdom))
sort(unique(as.data.frame(tax_table(physeq_obj_ITS_uparse_R1))$Phylum))
sort(unique(as.data.frame(tax_table(physeq_obj_ITS_uparse_R1))$RDP))

physeq_obj_ITS_uparse_R1 <- subset_taxa(physeq_obj_ITS_uparse_R1, RDP!="Rhizaria" & Phylum!="Cercozoa" &
                                          RDP!="Chromista" & RDP!="Metazoa")
physeq_obj_ITS_uparse_R1

# function to remove taxa present in the negative controls --------------------------------------------
remove_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

# non-target taxa to be removed - checkd using BLAST against GenBank
bad_taxa <- c("OTU_56","OTU_4363","OTU_5787","OTU_5599","OTU_4924","OTU_4720","OTU_4382",
              "OTU_3569","OTU_3401","OTU_5811","OTU_5307","OTU_5696","OTU_5650","OTU_5302","OTU_5089",
              "OTU_4503","OTU_4258","OTU_3898","OTU_3248","OTU_3134","OTU_2943","OTU_2937","OTU_2295",
              "OTU_5592","OTU_5324","OTU_5309","OTU_5229","OTU_5158","OTU_4818","OTU_4269","OTU_4217",
              "OTU_4132","OTU_4093","OTU_3914","OTU_3878","OTU_3543","OTU_2607","OTU_2360","OTU_5516",
              "OTU_4618","OTU_4181","OTU_4025","OTU_3755","OTU_2702","OTU_3311","OTU_1837","OTU_5188",
              "OTU_996","OTU_5723","OTU_5689","OTU_5326","OTU_5182","OTU_4945","OTU_4824","OTU_3423",
              "OTU_5852","OTU_5680","OTU_5091","OTU_4779","OTU_3981","OTU_3391","OTU_3252","OTU_2867",
              "OTU_5133","OTU_4384","OTU_4201","OTU_3609","OTU_3367","OTU_2796","OTU_4545","OTU_5367",
              "OTU_3317","OTU_3133","OTU_2723","OTU_2494","OTU_2475","OTU_1543","OTU_5839","OTU_5743",
              "OTU_5572","OTU_4837","OTU_3376","OTU_3143","OTU_2289","OTU_2155","OTU_2037","OTU_1584",
              "OTU_5380","OTU_4623","OTU_4607","OTU_3870","OTU_3869","OTU_3745","OTU_3340","OTU_1820",
              "OTU_2976","OTU_1449","OTU_5781","OTU_4537","OTU_4451","OTU_4159","OTU_3558","OTU_2154",
              "OTU_5187","OTU_4976","OTU_4806","OTU_4361","OTU_3823","OTU_1619","OTU_5123","OTU_4939",
              "OTU_4230","OTU_3908","OTU_3318","OTU_2991","OTU_2541","OTU_2391","OTU_1774","OTU_5477",
              "OTU_4812","OTU_4789","OTU_4709","OTU_3928","OTU_3604","OTU_3076","OTU_2721","OTU_2321",
              "OTU_2090","OTU_5580","OTU_4678","OTU_3915","OTU_3892","OTU_3868","OTU_2705","OTU_3952",
              "OTU_2560","OTU_2336","OTU_1997","OTU_1638","OTU_5768","OTU_5602","OTU_5597","OTU_5223",
              "OTU_4927","OTU_3772","OTU_3409","OTU_3388","OTU_3160","OTU_2755","OTU_2134","OTU_2085",
              "OTU_5185","OTU_5021","OTU_4989","OTU_4130","OTU_2471","OTU_1914","OTU_1377","OTU_1127",
              "OTU_4097","OTU_3860","OTU_3846","OTU_1208","OTU_2941","OTU_5866","OTU_5479","OTU_4758",
              "OTU_3044","OTU_1973","OTU_4317","OTU_3978","OTU_4408","OTU_2026","OTU_2998","OTU_2040",
              "OTU_4988","OTU_2243","OTU_4171","OTU_3316","OTU_4296","OTU_3796","OTU_2869","OTU_2337",
              "OTU_2213","OTU_3758","OTU_5159","OTU_5109","OTU_2891","OTU_3166","OTU_3371","OTU_4975",
              "OTU_5439","OTU_5335","OTU_4952","OTU_5892","OTU_5656","OTU_5629","OTU_4050")
bad_taxa

physeq_obj_ITS_uparse_R1 <- remove_taxa(physeq_obj_ITS_uparse_R1, bad_taxa)
physeq_obj_ITS_uparse_R1
tax_table(physeq_obj_ITS_uparse_R1)
sample_data(physeq_obj_ITS_uparse_R1)

# exporting datasets ----------------------------------------------------------------------------------
write.csv(tax_table(physeq_obj_ITS_uparse_R1), "taxonomy_ITS.csv")
write.csv(otu_table(physeq_obj_ITS_uparse_R1), "otu_table_ITS.csv")
write.dna(refseq(physeq_obj_ITS_uparse_R1), format="fasta", colsep="", file="rep_seq_ITS.fasta")

# extarcting taxa belonging to Endogonales ------------------------------------------------------------
sort(unique(as.data.frame(as.matrix(tax_table(physeq_obj_ITS_uparse_R1)))$Kingdom))
sort(unique(as.data.frame(as.matrix(tax_table(physeq_obj_ITS_uparse_R1)))$Phylum))
sort(unique(as.data.frame(as.matrix(tax_table(physeq_obj_ITS_uparse_R1)))$Class))
sort(unique(as.data.frame(as.matrix(tax_table(physeq_obj_ITS_uparse_R1)))$Order)) # Endogonales 

physeq_ITS_Endogonales <- subset_taxa(physeq_obj_ITS_uparse_R1, Order == "Endogonales")
physeq_ITS_Endogonales
otu_table(physeq_ITS_Endogonales)
tax_table(physeq_ITS_Endogonales)

sum(taxa_sums(physeq_ITS_Endogonales))
sum(taxa_sums(physeq_ITS_Endogonales))*100/sum(taxa_sums(physeq_obj_ITS_uparse_R1))

# >>> INSPECT LIBRARY SIZES ----------------------------------------------------------------------------
# Make a data frame with a column for the read counts of each sample
sample_sum_ITS_uparse_R1 <- data.frame(sum = sample_sums(physeq_obj_ITS_uparse_R1))
sample_sum_ITS_uparse_R1

# Histogram of sample read counts
ggplot(sample_sum_ITS_uparse_R1, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 10000) +
  labs(title="Distribution of sample sequencing depth",
       x="Read counts", y="Number of samples") + 
  xlim(0, 1100000) +
  theme(plot.title = element_text(size = 14, hjust = 0.5))

# >>> FILTERING ----------------------------------------------------------------------------------------
# Oliver et al. 2015, PCR errors and tag switching
# Lindhal et al. 2013, tag switching 

physeq_obj_ITS_uparse_R1 -> physeq_obj_ITS_uparse_R1_filt
otu_table(physeq_obj_ITS_uparse_R1_filt)[otu_table(physeq_obj_ITS_uparse_R1_filt) <= 4] <- 0 ### tag switching
otu_table(physeq_obj_ITS_uparse_R1_filt) <- otu_table(physeq_obj_ITS_uparse_R1_filt)[which(rowSums(otu_table(physeq_obj_ITS_uparse_R1_filt)) >= 10),] ### PCR Errors 
physeq_obj_ITS_uparse_R1_filt

# setting-up variables order --------------------------------------------------------------------------
colnames(sample_data(physeq_obj_ITS_uparse_R1_filt)) <- c("BarcodeSequence","LinkerPrimerSequence","Plant")
sample_data(physeq_obj_ITS_uparse_R1_filt)


# add a avraible with plant names abbreviations -------------------------------------------------------
library("vegan")
make.cepnames(c("Lycopodium volubile","Huperzia australiana","Lycopodiella fastigiatum","Lycopodium deuterodensum",
                "Lycopodiella lateralis","Phlegmariurus varius","Lycopodiella cernua","Lycopodiella diffusa",
                "Lycopodium scariosum","Selaginella kraussiana")) -> plant_cepnames

plant_cepnames
sample_data(physeq_obj_ITS_uparse_R1_filt)$Cepname <- plant_cepnames

sample_data(physeq_obj_ITS_uparse_R1_filt)$Plant <- factor(sample_data(physeq_obj_ITS_uparse_R1_filt)$Plant,
                                                           levels=c("Lycopodium_deuterodensum",
                                                                    "Lycopodium_scariosum",
                                                                    "Lycopodium_volubile",
                                                                    "Lycopodiella_cernua",
                                                                    "Lycopodiella_diffusa",
                                                                    "Lycopodiella_fastigiatum",
                                                                    "Lycopodiella_lateralis",
                                                                    "Huperzia_australiana",
                                                                    "Phlegmariurus_varius",
                                                                    "Selaginella_kraussiana")) 

sample_data(physeq_obj_ITS_uparse_R1_filt)$Cepname <- factor(sample_data(physeq_obj_ITS_uparse_R1_filt)$Cepname,
                                                             levels=c("Lycodeut",
                                                                      "Lycoscar",
                                                                      "Lycovolu",
                                                                      "Lycocern",
                                                                      "Lycodiff",
                                                                      "Lycofast",
                                                                      "Lycolate",
                                                                      "Hupeaust",
                                                                      "Phlevari",
                                                                      "Selakrau"))
# >>> RAREFACTION CURVES ------------------------------------------------------------------------------

# converting to table objetcs -------------------------------------------------------------------------
otu_fungi_ITS_uparse_R1 <- as.data.frame(otu_table(physeq_obj_ITS_uparse_R1_filt))
taxa_fungi_ITS_uparse_R1 <- as.data.frame(as.matrix(tax_table(physeq_obj_ITS_uparse_R1_filt)))
metadata_fungi_ITS_uparse_R1 <- as.data.frame(as.matrix(sample_data(physeq_obj_ITS_uparse_R1_filt)))
dim(otu_fungi_ITS_uparse_R1)
metadata_fungi_ITS_uparse_R1

# *** FIGURE S2 - rarefaction curves -----------------------------------------------------------------
 metadata_fungi_ITS_uparse_R1$color <- c("red","brown","orange","yellow1","green",
                                        "cyan","blue","deeppink","grey","black") 
metadata_fungi_ITS_uparse_R1

rarecurve(t(otu_fungi_ITS_uparse_R1), col = metadata_fungi_ITS_uparse_R1$color, label = FALSE, 
          sample=min(colSums(otu_fungi_ITS_uparse_R1)), step = 1000, 
          main="Fungi - ITS", ylab = "Number of OTUs", xlab = "Number of DNA reads") -> rare_fungi_uparse_R1

legend("bottomright", legend=metadata_fungi_ITS_uparse_R1$Cepname,
       col=metadata_fungi_ITS_uparse_R1$color, lty=1, cex=0.8, box.lty=1) #box.lty=0 remove the legend border

# >>> DATA RAREFACTION -------------------------------------------------------------------------------
set.seed(2019)

min(sample_sums(physeq_obj_ITS_uparse_R1_filt))
data.frame(colSums(otu_table(physeq_obj_ITS_uparse_R1_filt)))

physeq_fungi_uparse_R1_ev = rarefy_even_depth(physeq_obj_ITS_uparse_R1_filt, sample.size = min(sample_sums(physeq_obj_ITS_uparse_R1_filt)), 
                                              rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE) 
otu_table(physeq_fungi_uparse_R1_ev) <- otu_table(physeq_fungi_uparse_R1_ev)[which(rowSums(otu_table(physeq_fungi_uparse_R1_ev)) >= 1),]
physeq_fungi_uparse_R1_ev

colSums(otu_table(physeq_fungi_uparse_R1_ev))
any(taxa_sums(physeq_fungi_uparse_R1_ev) == 0)

# >>> TAXONOMIC COMPOSITION ---------------------------------------------------------------------------
library("dplyr")

otu_fungi_uparse_R1_phylum <- physeq_fungi_uparse_R1_ev %>%
  tax_glom(taxrank = "Class") %>%                     
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>%                                         
  #filter(Abundance > 0.01) %>%                         
  arrange(Class)                                     

otu_fungi_uparse_R1_phylum
head(otu_fungi_uparse_R1_phylum)
unique(otu_fungi_uparse_R1_phylum$Phylum)
unique(otu_fungi_uparse_R1_phylum$Class)

# *** FIGURE 1A - Fungi barplot -----------------------------------------------------------------------
library("ggpubr")

otu_fungi_uparse_R1_phylum
levels(otu_fungi_uparse_R1_phylum$Class)
levels(otu_fungi_uparse_R1_phylum$Class)[16] <- "Unclassified"

otu_fungi_uparse_R1_phylum$Class <- factor(otu_fungi_uparse_R1_phylum$Class,levels=c("Agaricomycetes","Agaricostilbomycetes","Archaeorhizomycetes","Archaeosporomycetes","Atractiellomycetes","Cystobasidiomycetes","Dacrymycetes",
                                                                                     "Dothideomycetes","Endogonomycetes","Eurotiomycetes","Exobasidiomycetes","Geoglossomycetes","Glomeromycetes","Gs11",
                                                                                     "Gs37","Kickxellomycetes","Lecanoromycetes","Leotiomycetes","Malasseziomycetes","Microbotryomycetes",
                                                                                     "Mortierellomycetes","Mucoromycetes","Olpidiomycetes","Orbiliomycetes","Paraglomeromycetes","Pezizomycetes","Pucciniomycetes",
                                                                                     "Rhizophydiomycetes","Saccharomycetes","Sordariomycetes","Taphrinomycetes","Tremellomycetes","Tritirachiomycetes","Umbelopsidomycetes",
                                                                                     "Ustilaginomycetes","Wallemiomycetes","Xylonomycetes","Zoopagomycetes","Unclassified"))

barplot_fungi_uparse_R1 = ggplot(otu_fungi_uparse_R1_phylum, aes(x = Cepname, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  coord_flip() +
  labs(title="Fungi ITS", x="", y = "Relative Abundance") +
  scale_fill_manual(values = palette_fungi) +
  #facet_grid(~Plant, scales = "free_x", space="free_x") +
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
  theme(strip.text.x = element_text(size = 8, face = "bold")) +
  theme(axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 0.5)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  #theme(axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")

barplot_fungi_uparse_R1

# >>> ABUNDANCES --------------------------------------------------------------------------------------
physeq_fungi_phylum <- tax_glom(physeq_fungi_uparse_R1_ev, "Phylum")
otu_fungi_abund_phy <- taxa_sums(physeq_fungi_phylum)/sum(taxa_sums(physeq_fungi_phylum))*100
tax_abund_fungi_phy <- as(tax_table(physeq_fungi_phylum), "matrix")
tax_abund_fungi_phy <- as.data.frame(tax_abund_fungi_phy)
tax_abund_fungi_phy <- tax_abund_fungi_phy[c(2)]
tax_abund_fungi_phy$abundance <- as.vector(otu_fungi_abund_phy)
tax_abund_fungi_phy <- tax_abund_fungi_phy[order(tax_abund_fungi_phy$abundance, decreasing = TRUE),] 
tax_abund_fungi_phy

physeq_fungi_class = tax_glom(physeq_fungi_uparse_R1_ev, "Class")
otu_fungi_abund_class = taxa_sums(physeq_fungi_class)/sum(taxa_sums(physeq_fungi_class))*100
tax_abund_fungi_class <- as(tax_table(physeq_fungi_class), "matrix")
tax_abund_fungi_class <- as.data.frame(tax_abund_fungi_class)
tax_abund_fungi_class <- tax_abund_fungi_class[c(1:3)]
tax_abund_fungi_class$abundance <- as.vector(otu_fungi_abund_class)
tax_abund_fungi_class <- tax_abund_fungi_class[order(tax_abund_fungi_class$abundance, decreasing = TRUE),] 
tax_abund_fungi_class

physeq_fungi_order = tax_glom(physeq_fungi_uparse_R1_ev, "Order")
otu_fungi_abund_ord = taxa_sums(physeq_fungi_order)/sum(taxa_sums(physeq_fungi_order))*100
tax_abund_fungi_ord <- as(tax_table(physeq_fungi_order), "matrix")
tax_abund_fungi_ord <- as.data.frame(tax_abund_fungi_ord)
tax_abund_fungi_ord <- tax_abund_fungi_ord[c(1:4)]
tax_abund_fungi_ord$abundance <- as.vector(otu_fungi_abund_ord)
tax_abund_fungi_ord <- tax_abund_fungi_ord[order(tax_abund_fungi_ord$abundance, decreasing = TRUE),] 
tax_abund_fungi_ord

physeq_fungi_genus = tax_glom(physeq_fungi_uparse_R1_ev, "Genus")
otu_fungi_abund_gen = taxa_sums(physeq_fungi_genus)/sum(taxa_sums(physeq_fungi_genus))*100
tax_abund_fungi_gen <- as(tax_table(physeq_fungi_genus), "matrix")
tax_abund_fungi_gen <- as.data.frame(tax_abund_fungi_gen)
tax_abund_fungi_gen <- tax_abund_fungi_gen[c(1:6)]
tax_abund_fungi_gen$abundance <- as.vector(otu_fungi_abund_gen)
tax_abund_fungi_gen <- tax_abund_fungi_gen[order(tax_abund_fungi_gen$abundance, decreasing = TRUE),] 
tax_abund_fungi_gen

# >>> COMAPRING DISTANCE MATRICES ---------------------------------------------------------------------
library("vegan")
library("dendextend")

otu_fungi_ev <- as.data.frame(otu_table(physeq_fungi_uparse_R1_ev))
taxa_fungi_ev <- as.data.frame(as.matrix(tax_table(physeq_fungi_uparse_R1_ev)))
metadata_fungi_ev <- as.data.frame(as.matrix(sample_data(physeq_fungi_uparse_R1_ev)))

identical(colnames(otu_fungi_ev), rownames(metadata_fungi_ev))

colnames(otu_fungi_ev) <- metadata_fungi_ev$Cepname
head(otu_fungi_ev)

# Generatingbray-curtis fungal distance ---------------------------------------------------------------
vegan::vegdist(t(otu_fungi_ev), method="bray") -> dist_fungi
dist_fungi <- as.dist(dist_fungi)
str(dist_fungi)
dist_fungi

# importin plant genetic distances --------------------------------------------------------------------
read.csv("dist_91.csv", header = TRUE, row.names = 1) -> plant_df_91
dist_91 <- as.dist(plant_df_91)
str(dist_91)
dist_91

dist_91_reordered <- as.dist(as.matrix(dist_91)[attr(dist_fungi, "Labels"),attr(dist_fungi, "Labels")])
all.equal(dist_fungi, dist_91_reordered, check.attributes = FALSE)
dist_91_reordered

# mantel test between the two distances ---------------------------------------------------------------
vegan::mantel(dist_fungi, dist_91_reordered, method="spearman", permutations=9999)


# >>> GENERATING CLUSTERS -----------------------------------------------------------------------------
fungi_clust <-hclust(dist_fungi,"complete")
plot(fungi_clust)

fungi_dendro <- as.dendrogram(fungi_clust)
plot(fungi_dendro)

# Importing the phylogenetic tree form newhick format ------------------------------------------------
read.tree(file = "Lycophyte_rbcL_ALL18Nov_Newick_3.tre") -> plant_tree
is.ultrametric(plant_tree)
chronos(plant_tree) -> plant_tree_dendro
is.ultrametric(plant_tree_dendro)
plant_tree_dendro

plant_tree_dendro$tip.label
metadata_fungi_ev$Cepname

plot(plant_tree_dendro)

# *** FIGURE 4A - fungi ITS tanglegram ----------------------------------------------------------------

dend_list <- dendlist(as.dendrogram(plant_tree_dendro), fungi_dendro)

# untangling the trees before calculating the entanglemnt 
dend_list_unt<- dend_list %>% untangle(method = "random") 
dend_list_unt

dend_list_unt %>% tanglegram(common_subtrees_color_lines=TRUE, 
             common_subtrees_color_branches = TRUE,
             highlight_distinct_edges=FALSE, 
             highlight_branches_lwd=FALSE,
             main_left = "Plant\nPhylogeny",
             main_right = "Community\nDissimilarity",
             main="Fungi",
             sub = paste("entanglement =", round(entanglement(dend_list_unt,
                              L=1, leaves_matching_method = "labels"), 2)),
             hang = FALSE, 
             type = "r",
             margin_inner=4.5,
             lwd=2.5,
             lab.cex = 1,
             cex_main = 1.2,
             cex_sub=1,
             center = FALSE)

# calculating entanglement ----------------------------------------------------------------------------
entanglement(dend_list_unt, L=1, leaves_matching_method = "labels")

# >>> HEATMAP TREES -----------------------------------------------------------------------------------
library("metacoder")
library("tibble")

biom_ITS_gen = tax_glom(physeq_fungi_uparse_R1_ev, taxrank="Genus")
head(otu_table(biom_ITS_gen))
head(tax_table(biom_ITS_gen))

tax_table(biom_ITS_gen)[tax_table(biom_ITS_gen)=="Unclassified"]<- ""
tax_table(biom_ITS_gen) <- tax_table(biom_ITS_gen)[,1:6]

otu_fungi_gen <- as.data.frame(otu_table(biom_ITS_gen))
otu_fungi_gen$sample_id <- rownames(otu_fungi_gen)
head(otu_fungi_gen)

metadata_fungi_gen <- as.data.frame(as.matrix(sample_data(biom_ITS_gen)))
metadata_fungi_gen$sample_id <- rownames(metadata_fungi_gen)
metadata_fungi_gen$Cepname <- as.character(metadata_fungi_gen$Cepname)
metadata_fungi_gen
str(metadata_fungi_gen)

# reformatting taxonomies for the metacoder object ----------------------------------------------------
# colnames(otu_fungi_gen) <- metadata_fungi_gen$Cepname
# rownames(metadata_fungi_gen) <- metadata_fungi_gen$Cepname
taxa_fungi_gen <- as.data.frame(as.matrix(tax_table(biom_ITS_gen)))
taxa_fungi_gen <- taxa_fungi_gen[,c(1:6)]
cols <- c("Kingdom","Phylum","Class","Order","Family","Genus")
taxa_fungi_gen$lineage <- do.call(paste, c(taxa_fungi_gen[cols], sep=";"))
head(taxa_fungi_gen)

otu_fungi_gen$taxonomy <- taxa_fungi_gen$lineage
head(otu_fungi_gen)

# transform a data.frame into a tibble 
tibbl_fungi <- tibble::as_data_frame(otu_fungi_gen)
tibbl_fungi

tibbl_metadata_fungi <- tibble::as_data_frame(metadata_fungi_gen)
tibbl_metadata_fungi

# creating the metacoder object ------------------------------------------------------------------------ 
heat_tree_ITS <- parse_tax_data(tibbl_fungi,
                                class_cols = "taxonomy",
                                class_sep = ";")
print(heat_tree_ITS)
heat_tree_ITS$data$tax_data
#data_test <- tibbl_fungi[1:500,1:14]

# Removing low-abundance counts -----------------------------------------------------------------------
heat_tree_ITS$data$tax_data <- zero_low_counts(heat_tree_ITS, "tax_data", min_count = 5)
no_reads <- rowSums(heat_tree_ITS$data$tax_data[, tibbl_metadata_fungi$sample_id]) == 0
sum(no_reads)

# removing OTUs with no reads left 
heat_tree_ITS <- filter_obs(heat_tree_ITS, data  = "tax_data", ! no_reads, drop_taxa = TRUE)
print(heat_tree_ITS)

# calculate taxa proportions --------------------------------------------------------------------------
heat_tree_ITS$data$tax_data <- calc_obs_props(heat_tree_ITS, "tax_data")

# calculating abundance per-taxon ---------------------------------------------------------------------
heat_tree_ITS$data$tax_abund <- calc_taxon_abund(heat_tree_ITS, "tax_data", cols = tibbl_metadata_fungi$sample_id)

# adding taconomic occurrences ------------------------------------------------------------------------
heat_tree_ITS$data$tax_occ <- calc_n_samples(heat_tree_ITS, 
                                             "tax_abund", 
                                             cols = tibbl_metadata_fungi$sample_id,
                                             groups = tibbl_metadata_fungi$Cepname)

heat_tree_ITS$data$tax_occ <- calc_n_samples(heat_tree_ITS, 
                                             "tax_abund", 
                                             groups = tibbl_metadata_fungi$Cepname)

print(heat_tree_ITS)

# Test plotting the heatmap tree
set.seed(1)
heat_tree(heat_tree_ITS, 
          node_label = heat_tree_ITS$taxon_names(),
          node_size = heat_tree_ITS$n_obs(),
          node_color = heat_tree_ITS$n_obs(), 
          overlap_avoidance = 0.01,
          node_label_size_range = c(0.01, 0.04), 
          tree_label_size = as.numeric(0.6),
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

# Plotting read depth ---------------------------------------------------------------------------------
# To plot read depth, you first need to add up the number of reads per taxon.
# The function `calc_taxon_abund` is good for this. 
heat_tree_ITS$data$taxon_counts <- calc_taxon_abund(heat_tree_ITS, data = "tax_data")
heat_tree_ITS$data$taxon_counts$total <- rowSums(heat_tree_ITS$data$taxon_counts[, -1]) # -1 = taxon_id column

# *** FIGURE XXX - heatmpa tree test plot ------------------------------------------------------------- 
heat_tree(heat_tree_ITS, 
          node_label = taxon_names, 
          node_size = total, 
          node_color = Selakrau,
          overlap_avoidance = 0.01,
          node_color_range = c("#FFFFFF", "darkorange3", "#4e567d", "gold"),
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations)


heat_tree_ITS$data$diff_table <- compare_groups(heat_tree_ITS,
                                                dataset = "tax_abund",
                                                cols = tibbl_metadata_fungi$sample_id, # What columns of sample data to use
                                                groups = tibbl_metadata_fungi$Cepname) # What category each sample is assigned to
print(heat_tree_ITS$data$diff_table)
range(heat_tree_ITS$data$diff_table$wilcox_p_value, finite = TRUE) 

filter(heat_tree_ITS$data$diff_table, treatment_2 == "Selakrau")

set.seed(999)
heat_tree(heat_tree_ITS, 
          node_label = taxon_names,
          node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
          node_color = log2_median_ratio, # A column from `obj$data$diff_table`
          node_color_interval = c(-2, 2), # The range of `log2_median_ratio` to display
          node_color_range = c("cyan", "gray", "tan"), # The color palette used
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Log 2 ratio of median proportions",
          overlap_avoidance = 0.01,
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

# mergins samples to compare all Lycopodiaceae with the outgrpup --------------------------------------
head(otu_fungi_gen)

otu_fungi_gen2 <- otu_fungi_gen
otu_fungi_gen2$Others <- rowSums(otu_fungi_gen2[,c(1:9)]) 
otu_fungi_gen2$Outgroup <- otu_fungi_gen$sample10
otu_fungi_gen2$taxonomy <- otu_fungi_gen$taxonomy
otu_fungi_gen2 <- otu_fungi_gen2[,11:14]
otu_fungi_gen2 <- otu_fungi_gen2[,c(3,4,1,2)]
head(otu_fungi_gen2)

metadata_fungi_gen2 <- metadata_fungi_gen[9:10,]
rownames(metadata_fungi_gen2) <- c("Others","Outgroup")
metadata_fungi_gen2$Cepname <- c("Others","Outgroup")
metadata_fungi_gen2$sample_id <- c("Others","Outgroup")
metadata_fungi_gen2

# transform a data.frame into a tibble ----------------------------------------------------------------
tibbl_fungi2 <- as_tibble(otu_fungi_gen2)
tibbl_fungi2

tibbl_metadata_fungi2 <- as_tibble(metadata_fungi_gen2)
tibbl_metadata_fungi2

# create metacoder object 
heat_tree_ITS_2 <- parse_tax_data(tibbl_fungi2,
                                  class_cols = "taxonomy",
                                  class_sep = ";")
print(heat_tree_ITS_2)
heat_tree_ITS_2$data$tax_data

#data_test <- tibbl_fungi[1:500,1:14]

# Removing low-abundance counts -----------------------------------------------------------------------
heat_tree_ITS_2$data$tax_data <- zero_low_counts(heat_tree_ITS_2, "tax_data", min_count = 5)
no_reads <- rowSums(heat_tree_ITS_2$data$tax_data[, tibbl_metadata_fungi2$sample_id]) == 0
sum(no_reads)

# removing OTUs with no reads left 
heat_tree_ITS_2 <- filter_obs(heat_tree_ITS_2, data  = "tax_data", ! no_reads, drop_taxa = TRUE)
print(heat_tree_ITS_2)

# calculate taxa proportions --------------------------------------------------------------------------
heat_tree_ITS_2$data$tax_data <- calc_obs_props(heat_tree_ITS_2, "tax_data")

# abundance per-taxon ---------------------------------------------------------------------------------
heat_tree_ITS_2$data$tax_abund <- calc_taxon_abund(heat_tree_ITS_2, "tax_data", 
                                                   cols = tibbl_metadata_fungi2$sample_id)

# number of samples have reads for each taxon ---------------------------------------------------------
heat_tree_ITS_2$data$tax_occ <- calc_n_samples(heat_tree_ITS_2, "tax_abund", 
                                               cols = tibbl_metadata_fungi2$sample_id,
                                               groups = tibbl_metadata_fungi2$Cepname)
print(heat_tree_ITS_2)


# comparing the Outgroup with Others
heat_tree_ITS_2$data$diff_table <- compare_groups(heat_tree_ITS_2,
                                                  data = "tax_abund",
                                                  func = NULL,
                                                  cols = tibbl_metadata_fungi2$sample_id, # What columns of sample data to use
                                                  groups = tibbl_metadata_fungi2$Cepname) # What category each sample is assigned to

print(heat_tree_ITS_2$data$diff_table)
heat_tree_ITS_2$data$tax_abund
heat_tree_ITS_2$taxon_names()

# function used to compare the groups,
# p.value cannot be calculated because we have no replicates!

function(abund_1, abund_2) {
  log_ratio <- log2(median(abund_1) / median(abund_2))
  if (is.nan(log_ratio)) {
    log_ratio <- 0
  }
  list(log2_median_ratio = log_ratio,
       median_diff = median(abund_1) - median(abund_2),
       mean_diff = mean(abund_1) - mean(abund_2),
       wilcox_p_value = wilcox.test(abund_1, abund_2)$p.value)
}

range(heat_tree_ITS_2$data$diff_table$median_diff, finite = TRUE) 
range(heat_tree_ITS_2$data$diff_table$mean_diff, finite = TRUE) 


# *** FIGURE 3 - Heatmap tree fungi ITS ---------------------------------------------------------------
set.seed(999)
heat_tree(heat_tree_ITS_2, 
          node_label = heat_tree_ITS_2$taxon_names(),
          node_size = heat_tree_ITS_2$n_obs(), # n_obs is a function that calculates, in this case, the number of OTUs per taxon
          node_color = heat_tree_ITS_2$data$diff_table$log2_median_ratio, # A column from `obj$data$diff_table`
          node_color_interval = c(-5, 5), # The range of `log2_median_ratio` to display
          node_color_range = c("tomato", "gray", "skyblue1"), # The color palette used
          #node_size_axis_label = "OTU count",
          #node_color_axis_label = "Log 2 ratio of median proportions",
          node_label_size_range = c(0.01, 0.042), 
          tree_label_size = as.numeric(0.6),
          #overlap_avoidance = 0.01,
          node_color_axis_label = "Lycopodaceae vs. Outgroup",
          node_size_axis_label = "Number of Genera",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations


# UPSET R ----------------------------------------------------------------------------------------------
library("UpSetR")
library("dplyr")

head(otu_fungi_ev)

otu_fungi_ev -> otu_fungi_ev_pa
otu_fungi_ev_pa[otu_fungi_ev_pa>0] <- 1
head(otu_fungi_ev_pa)

rownames(metadata_fungi_ev) <- metadata_fungi_ev$Cepname
metadata_fungi_ev

sample_order <- match(colnames(otu_fungi_ev_pa), rownames(metadata_fungi_ev))
otu_fungi_ev_pa <- otu_fungi_ev_pa[, sample_order]
identical(colnames(otu_fungi_ev_pa), rownames(metadata_fungi_ev)) # now TRUE

# *** FIGURE S4 - upset plot fungi ITS ----------------------------------------------------------------
upset(otu_fungi_ev_pa, 
      sets=metadata_fungi_ev$Cepname,
      #nintersects = 100, # if NA all intersections are plotted
      main.bar.color = "black", 
      group.by = "degree",
      cutoff = 200,
      number.angles = 0, 
      point.size = 2, 
      line.size = 1,
      mainbar.y.label = "OTU number", 
      sets.x.label = "OTU number",
      text.scale = c(1.3, 1.3, 1, 1, 1, 0.75),
      #text.scale = c(1, 1, 1, 1, 1, 1),
      empty.intersections = "on") -> upset_fungi

upset_fungi

# *******************************************************************************************---
# Import fungal 18S data (generic primers NS31_f3, AML2_f3) ------------------------------------
# hereafter calles SSU dataset for simplicity

otus_SSU_uparse_R1 <- read.delim("analysis_SSU/otu_table_ITS_UPARSE_R1.txt",row.names=1) 
otus_phy_SSU_uparse_R1 <-otu_table(otus_SSU_uparse_R1, taxa_are_rows = TRUE)

metadata_SSU_uparse_R1 <-read.delim("analysis_SSU/mapping_SSU.txt", row.names=1, header=TRUE, sep="\t")
metadata_phy_SSU_uparse_R1 <-sample_data(metadata_SSU_uparse_R1)

taxonomy_SSU_uparse_R1 <-read.csv("analysis_SSU/otu_SSU_R1_SILVA_filt.csv", header=TRUE, row.names=1)
head(taxonomy_SSU_uparse_R1)

# filtering SILVA taxonomy using .07 confidence value at each Rank --------------------------------
ifelse(taxonomy_SSU_uparse_R1$conf>=0.7, paste(taxonomy_SSU_uparse_R1$Kingdom), NA) -> Kingdom
ifelse(taxonomy_SSU_uparse_R1$conf.1>=0.7, paste(taxonomy_SSU_uparse_R1$Phylum), NA) -> Phylum
ifelse(taxonomy_SSU_uparse_R1$conf.2>=0.7, paste(taxonomy_SSU_uparse_R1$Class), NA) -> Class
ifelse(taxonomy_SSU_uparse_R1$conf.3>=0.7, paste(taxonomy_SSU_uparse_R1$Order), NA) -> Order
ifelse(taxonomy_SSU_uparse_R1$conf.4>=0.7, paste(taxonomy_SSU_uparse_R1$Family), NA) -> Family
ifelse(taxonomy_SSU_uparse_R1$conf.4>=0.7, paste(taxonomy_SSU_uparse_R1$Family), NA) -> Genus # trick to make all NA
#ifelse(taxonomy_SSU_uparse_R1$conf.5>=0.7, paste(taxonomy_SSU_uparse_R1$Genus), NA) -> Genus
taxonomy_SSU <- cbind(Kingdom, Phylum, Class, Order, Family, Genus)
rownames(taxonomy_SSU) <- rownames(taxonomy_SSU_uparse_R1)
taxonomy_phy_SSU_uparse_R1 <- tax_table(as.matrix(taxonomy_SSU))
head(taxonomy_phy_SSU_uparse_R1)

otus_seq_SSU_uparse_R1 <- readDNAStringSet("analysis_SSU/otus_R1.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

physeq_obj_SSU_uparse_R1 <- phyloseq(otus_phy_SSU_uparse_R1, 
                                     metadata_phy_SSU_uparse_R1,
                                     taxonomy_phy_SSU_uparse_R1,
                                     otus_seq_SSU_uparse_R1) 

physeq_obj_SSU_uparse_R1
str(physeq_obj_SSU_uparse_R1)
tax_table(physeq_obj_SSU_uparse_R1)
sample_data(physeq_obj_SSU_uparse_R1)

tax_table(physeq_obj_SSU_uparse_R1)[tax_table(physeq_obj_SSU_uparse_R1)==""]<- NA
tax_table(physeq_obj_SSU_uparse_R1)[is.na(tax_table(physeq_obj_SSU_uparse_R1))]<-"Unclassified"

# exporting datasets ----------------------------------------------------------------------------------
write.csv(tax_table(physeq_obj_SSU_uparse_R1), "taxonomy_physeq_SSU.csv")
write.dna(refseq(physeq_obj_SSU_uparse_R1), format="fasta", colsep="", file="sequences_physeq_SSU.fasta")

# > ABUNDANCES -----------------
physeq_fungi_phylum_SSU <- tax_glom(physeq_obj_SSU_uparse_R1, "Phylum")
otu_fungi_abund_phy_SSU <- taxa_sums(physeq_fungi_phylum_SSU)/sum(taxa_sums(physeq_fungi_phylum_SSU))*100
tax_abund_fungi_phy_SSU <- as(tax_table(physeq_fungi_phylum_SSU), "matrix")
tax_abund_fungi_phy_SSU <- as.data.frame(tax_abund_fungi_phy_SSU)
tax_abund_fungi_phy_SSU <- tax_abund_fungi_phy_SSU[c(2)]
tax_abund_fungi_phy_SSU$abundance <- as.vector(otu_fungi_abund_phy_SSU)
tax_abund_fungi_phy_SSU <- tax_abund_fungi_phy_SSU[order(tax_abund_fungi_phy_SSU$abundance, decreasing = TRUE),] 
tax_abund_fungi_phy_SSU

# Proportion of Endogonales in the SSU ----------------------------------------------------------------
taxonomy_ranks <- as.matrix(tax_table(physeq_obj_SSU_uparse_R1))[,c("Kingdom","Phylum","Class","Order","Family","Genus")]
taxonomy_ranks
tax_table(physeq_obj_SSU_uparse_R1) <- taxonomy_ranks

sort(unique(as.data.frame(as.matrix(tax_table(physeq_obj_SSU_uparse_R1)))$Kingdom))
sort(unique(as.data.frame(as.matrix(tax_table(physeq_obj_SSU_uparse_R1)))$Phylum)) # Glomero
sort(unique(as.data.frame(as.matrix(tax_table(physeq_obj_SSU_uparse_R1)))$Class))
sort(unique(as.data.frame(as.matrix(tax_table(physeq_obj_SSU_uparse_R1)))$Order)) 

physeq_SSU_Glomeromycota <- subset_taxa(physeq_obj_SSU_uparse_R1, Phylum == "Glomeromycota")
physeq_SSU_Glomeromycota
otu_table(physeq_SSU_Glomeromycota)
tax_table(physeq_SSU_Glomeromycota)

sum(taxa_sums(physeq_SSU_Glomeromycota))*100/sum(taxa_sums(physeq_obj_SSU_uparse_R1))

# extracting sequences for phylogeny ------------------------------------------------------------------
write.dna(refseq(physeq_SSU_Glomeromycota), format="fasta", colsep="", file="taxa_glomeromycota_all.fasta")

# **************************************************************************************************---
# Import fungal 18S data (Endogone Endo18S-1F, NS6_f4) ------------------------------------------------
# hereafter calles 18S dataset for simplicity

otus_18S_uparse_R1 <- read.delim("analysis_18S/otu_table_ITS_UPARSE_R1.txt",row.names=1) 
otus_phy_18S_uparse_R1 <-otu_table(otus_18S_uparse_R1, taxa_are_rows = TRUE)

metadata_18S_uparse_R1 <-read.delim("analysis_18S/mapping_18S.txt", row.names=1, header=TRUE, sep="\t")
metadata_phy_18S_uparse_R1 <-sample_data(metadata_18S_uparse_R1)

taxonomy_18S_uparse_R1 <-read.csv("analysis_18S/otu_R1_SILVA_filt.csv", header=TRUE, row.names=1)
taxonomy_phy_18S_uparse_R1 <- tax_table(as.matrix(taxonomy_18S_uparse_R1))
taxonomy_phy_18S_uparse_R1

otus_seq_18S_uparse_R1 <- readDNAStringSet("analysis_18S/otus_R1.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
otus_seq_18S_uparse_R1

physeq_obj_18S_uparse_R1 <- phyloseq(otus_phy_18S_uparse_R1, 
                                     metadata_phy_18S_uparse_R1,
                                     taxonomy_phy_18S_uparse_R1,
                                     otus_seq_18S_uparse_R1) 

physeq_obj_18S_uparse_R1
str(physeq_obj_18S_uparse_R1)
tax_table(physeq_obj_18S_uparse_R1)
sample_data(physeq_obj_18S_uparse_R1)

write.csv(tax_table(physeq_obj_18S_uparse_R1), "taxonomy_physeq_18S_uparse_R1.csv")
write.dna(refseq(physeq_obj_18S_uparse_R1), format="fasta", colsep="", file="sequences_physeq_18S_uparse_R1.fasta")

tax_table(physeq_obj_18S_uparse_R1)[tax_table(physeq_obj_18S_uparse_R1)==""]<- NA
tax_table(physeq_obj_18S_uparse_R1)[is.na(tax_table(physeq_obj_18S_uparse_R1))]<-"Unclassified"

# filtering SILVA taxonomy using .07 confidence value at each Rank ------------------------------------
ifelse(taxonomy_18S_uparse_R1$conf>=0.7, paste(taxonomy_18S_uparse_R1$Kingdom), NA) -> Kingdom
ifelse(taxonomy_18S_uparse_R1$conf.1>=0.7, paste(taxonomy_18S_uparse_R1$Phylum), NA) -> Phylum
ifelse(taxonomy_18S_uparse_R1$conf.2>=0.7, paste(taxonomy_18S_uparse_R1$Class), NA) -> Class
ifelse(taxonomy_18S_uparse_R1$conf.3>=0.7, paste(taxonomy_18S_uparse_R1$Order), NA) -> Order
ifelse(taxonomy_18S_uparse_R1$conf.4>=0.7, paste(taxonomy_18S_uparse_R1$Family), NA) -> Family
ifelse(taxonomy_18S_uparse_R1$conf.4>=0.7, paste(taxonomy_18S_uparse_R1$Family), NA) -> Genus # trick to make all NA
#ifelse(taxonomy_18S_uparse_R1$conf.5>=0.7, paste(taxonomy_18S_uparse_R1$Genus), NA) -> Genus
taxonomy_18S <- cbind(Kingdom, Phylum, Class, Order, Family, Genus)
rownames(taxonomy_18S) <- rownames(taxonomy_18S_uparse_R1)
taxonomy_phy_18S_uparse_R1 <- tax_table(as.matrix(taxonomy_18S))
head(taxonomy_phy_18S_uparse_R1)

otus_seq_18S_uparse_R1 <- readDNAStringSet("analysis_18S/otus_R1.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

physeq_obj_18S_uparse_R1 <- phyloseq(otus_phy_18S_uparse_R1, 
                                     metadata_phy_18S_uparse_R1,
                                     taxonomy_phy_18S_uparse_R1,
                                     otus_seq_18S_uparse_R1) 

physeq_obj_18S_uparse_R1
str(physeq_obj_18S_uparse_R1)
tax_table(physeq_obj_18S_uparse_R1)
sample_data(physeq_obj_18S_uparse_R1)

tax_table(physeq_obj_18S_uparse_R1)[tax_table(physeq_obj_18S_uparse_R1)==""]<- NA
tax_table(physeq_obj_18S_uparse_R1)[is.na(tax_table(physeq_obj_18S_uparse_R1))]<-"Unclassified"

# exporting datasets to manually check for Endogoanles ------------------------------------------------
write.csv(tax_table(physeq_obj_18S_uparse_R1), "taxonomy_physeq_18S.csv")
write.dna(refseq(physeq_obj_18S_uparse_R1), format="fasta", colsep="", file="sequences_physeq_18S.fasta")

# extracting interesting OTUS for phylogeny -----------------------------------------------------------
library("reshape2")

selected_otus_18S<- c("OTU_870","OTU_1002","OTU_1319","OTU_407", "OTU_1316","OTU_43","OTU_758")

physeq_otu_18S <- extract_taxa(physeq_obj_18S_uparse_R1, selected_otus_18S) 
sample_data(physeq_otu_18S)
otu_table(physeq_otu_18S)

sum(taxa_sums(physeq_otu_18S))*100/sum(taxa_sums(physeq_obj_18S_uparse_R1))

ordered_ceps <- c("Lycodeut","Lycodiff","Lycolate","Lycovolu","Lycocern",
                  "Selakrau","Hupeaust", "Lycofast","Phlevari","Lycoscar")

sample_data(physeq_otu_18S)$Cepname <- ordered_ceps
sample_data(physeq_otu_18S)

sample_data(physeq_otu_18S)$Cepname <- factor(sample_data(physeq_otu_18S)$Cepname,
                                              levels=c("Lycodeut",
                                                       "Lycoscar",
                                                       "Lycovolu",
                                                       "Lycocern",
                                                       "Lycodiff",
                                                       "Lycofast",
                                                       "Lycolate",
                                                       "Hupeaust",
                                                       "Phlevari",
                                                       "Selakrau"))

tax_table(physeq_otu_18S) <- cbind(tax_table(physeq_otu_18S), OTU_ID=taxa_names(physeq_otu_18S))
tax_table(physeq_otu_18S)

otu_table_18S <- as.data.frame(otu_table(physeq_otu_18S))
colnames(otu_table_18S) <- ordered_ceps
otu_table_18S

otu_table_18S <- otu_table_18S[c("Lycodeut","Lycoscar","Lycovolu","Lycocern","Lycodiff",
                                 "Lycofast","Lycolate","Hupeaust","Phlevari","Selakrau")]

otu_table_18S[otu_table_18S=="0"]<- NA
otu_table_18S

otu_table_18S$OTU_ID <- rownames(otu_table_18S) 

otu_table_18S_melted <-melt(otu_table_18S)
otu_table_18S_melted
str(otu_table_18S_melted)

colnames(otu_table_18S_melted) <- c("OTU_ID", "Species", "Abundance")
otu_table_18S_melted[otu_table_18S_melted=="0"]<- NA

otu_table_18S_melted$OTU_ID <- as.factor(otu_table_18S_melted$OTU_ID)
otu_table_18S_melted$OTU_ID <- factor(otu_table_18S_melted$OTU_ID, levels = c("OTU_758",
                                                                              "OTU_407",
                                                                              "OTU_1002",
                                                                              "OTU_1319",
                                                                              "OTU_1316",
                                                                              "OTU_43",
                                                                              "OTU_870"))

# *** FIGURE 2 - Endogonales abundance baloon plot-----------------------------------------------------
baloon_18S <- ggplot(otu_table_18S_melted, aes(x = Species, y = OTU_ID)) +
  labs(x="", y="") +
  geom_point(aes(size=Abundance),shape=21, colour="firebrick", fill="firebrick") +
  scale_size(range = c(0.002679528, 10)) +
  scale_y_discrete(limits = rev(levels(otu_table_18S_melted$OTU_ID))) +
  #scale_y_discrete(limits=c("OTU_407","OTU_1002","OTU_1319","OTU_1316","OTU_870","OTU_43","OTU_758")) +
  theme(axis.text.y = element_text(size = 10, colour = "black")) +
  theme(axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1)) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  theme(panel.background=element_blank(), panel.border = element_rect(colour = "grey", fill=NA, size=0.5))

baloon_18S

# **************************************************************************************************---
# Importing bacterial 16S data ------------------------------------------------------------------------
otus_16s_uparse <- read.delim("analysis_16S/otu_table_16S_UPARSE.txt",row.names=1) 
otus_phy_16s_uparse <- otu_table(otus_16s_uparse, taxa_are_rows = TRUE)

metadata_16s_uparse <- read.delim("analysis_16S/mapping_16s.txt", row.names=1, header=TRUE, sep="\t")
metadata_phy_16s_uparse <- sample_data(metadata_16s_uparse)

# importing RDP taxonomy ------------------------------------------------------------------------------
taxonomy_16s_uparse_RDP <- read.delim("analysis_16S/taxonomy_assignments/otus_taxonomy_RDP.txt", header = TRUE, row.names = 1)
head(taxonomy_16s_uparse_RDP)

ifelse(taxonomy_16s_uparse_RDP$D_Score>=0.7, paste(taxonomy_16s_uparse_RDP$Domain), NA) -> Kingdom
ifelse(taxonomy_16s_uparse_RDP$P_Score>=0.7, paste(taxonomy_16s_uparse_RDP$Phylum), NA) -> Phylum
ifelse(taxonomy_16s_uparse_RDP$C_Score>=0.7, paste(taxonomy_16s_uparse_RDP$Class), NA) -> Class
ifelse(taxonomy_16s_uparse_RDP$O_Score>=0.7, paste(taxonomy_16s_uparse_RDP$Order), NA) -> Order
ifelse(taxonomy_16s_uparse_RDP$F_Score>=0.7, paste(taxonomy_16s_uparse_RDP$Family), NA) -> Family
ifelse(taxonomy_16s_uparse_RDP$G_Score>=0.7, paste(taxonomy_16s_uparse_RDP$Genus), NA) -> Genus
taxonomy <- cbind(Kingdom, Phylum, Class, Order, Family, Genus)
rownames(taxonomy) <- rownames(taxonomy_16s_uparse_RDP)
taxonomy_phy_16s_uparse_RDP <- tax_table(as.matrix(taxonomy))
head(taxonomy_phy_16s_uparse_RDP)

zotus_seq_16s_uparse <- readDNAStringSet("analysis_16S/otus.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

physeq_obj_16s_uparse <- phyloseq(otus_phy_16s_uparse, 
                                  metadata_phy_16s_uparse,
                                  taxonomy_phy_16s_uparse_RDP,
                                  zotus_seq_16s_uparse) 

physeq_obj_16s_uparse
sum(taxa_sums(physeq_obj_16s_uparse))

head(sample_data(physeq_obj_16s_uparse))
head(tax_table(physeq_obj_16s_uparse))

tax_table(physeq_obj_16s_uparse)[tax_table(physeq_obj_16s_uparse)==""]<- NA
tax_table(physeq_obj_16s_uparse)[is.na(tax_table(physeq_obj_16s_uparse))]<-"Unclassified"

# filtering prokaryote taxonomy -----------------------------------------------------------------------
physeq_obj_16s_uparse -> physeq_prokaryote
physeq_prokaryote

# looking for chloroplast in the RDP taxonomy 
any(tax_table(physeq_obj_16s_uparse) == "Chloroplast")
sort(unique(as.data.frame(tax_table(physeq_obj_16s_uparse))$Phylum))

# removing cloroplast according RDP -------------------------------------------------------------------
physeq_prokaryote <- subset_taxa(physeq_obj_16s_uparse, Phylum != "Cyanobacteria/Chloroplast")
physeq_prokaryote

physeq_prokaryote_plastids <- subset_taxa(physeq_obj_16s_uparse, Phylum == "Cyanobacteria/Chloroplast")
physeq_prokaryote_plastids
otu_table(physeq_prokaryote_plastids)
sum(taxa_sums(physeq_prokaryote_plastids))
sum(taxa_sums(physeq_prokaryote_plastids))/sum(taxa_sums(physeq_obj_16s_uparse))*100


# removing mitochondria (and additional plastid reads) using ------------------------------------------
# SILVA reference database v 123 and SINTAX taxonomy assignemnts 
library("dplyr")
library("tidyr")
library("stringr")

connect_16s_uparse_SILVA <-readLines("analysis_16S/taxonomy_assignments/otu_taxonomy.sintax")
connect_16s_uparse_SILVA <- gsub("\tk:.*\t\\+\tk:", ",", connect_16s_uparse_SILVA) 
head(connect_16s_uparse_SILVA)
str(connect_16s_uparse_SILVA)

max(count.fields(textConnection(connect_16s_uparse_SILVA), sep = ","))

taxonomy_16s_uparse_SILVA <-read.table(textConnection(connect_16s_uparse_SILVA),
                                       header = FALSE, 
                                       sep = ",",
                                       col.names = paste0("V",seq_len(8)),
                                       fill = TRUE)

rownames(taxonomy_16s_uparse_SILVA) <- taxonomy_16s_uparse_SILVA$V1
colnames(taxonomy_16s_uparse_SILVA) <- c("OTU_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
str(taxonomy_16s_uparse_SILVA)
head(taxonomy_16s_uparse_SILVA)

otu_chlo_mito <- rownames(taxonomy_16s_uparse_SILVA[taxonomy_16s_uparse_SILVA$Class=="c:Chloroplast" |
                                                      taxonomy_16s_uparse_SILVA$Family=="f:Mitochondria", ])
otu_chlo_mito

# REMOVING UNWANTED taxa ------------------------------------------------------------------------------
remove_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

physeq_prokaryote <- remove_taxa(physeq_prokaryote, otu_chlo_mito)
physeq_prokaryote

# evaluating cloroplast and mitocondrion abundance according to SILVA db
physeq_prokaryote_mito <- otu_table(physeq_obj_16s_uparse)[otu_chlo_mito, ]
sum(taxa_sums(physeq_prokaryote_mito))
sum(taxa_sums(physeq_prokaryote_mito))/sum(taxa_sums(physeq_obj_16s_uparse))*100

# adding OTU_ID to the taxonomy table, always usefull!
tax_table(physeq_prokaryote) <- cbind(tax_table(physeq_prokaryote), OTU_ID=taxa_names(physeq_prokaryote))
tax_table(physeq_prokaryote) 

# reordering based on the most abundant 
tax_order <- names(sort(taxa_sums(physeq_prokaryote), TRUE))
as.data.frame(as.matrix(tax_table(physeq_prokaryote)))[tax_order,]

# exporting datasets ----------------------------------------------------------------------------------
write.csv(tax_table(physeq_prokaryote), "taxonomy_16S.csv")
write.csv(otu_table(physeq_prokaryote), "otu_table_16S.csv")
write.dna(refseq(physeq_prokaryote), format="fasta", colsep="", file="rep_seq_16S.fasta")

# >>> INSPECT LIBRARY SIZES ---------------------------------------------------------------------------

# Make a data frame with a column for the read counts of each sample
sample_sum_16s_uparse <- data.frame(sum = sample_sums(physeq_prokaryote))
sample_sum_16s_uparse

# Histogram of sample read counts
ggplot(sample_sum_16s_uparse, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 10000) +
  labs(title="Distribution of sample sequencing depth",
       x="Read counts", y="Number of samples") + 
  xlim(0,1000000) +
  theme(plot.title = element_text(size = 14, hjust = 0.5))

# >>> FILTERING ---------------------------------------------------------------------------------------
# Oliver et al. 2015, PCR errors and tag switching
# Lindhal et al. 2013, tag switching - that's a good  one!

physeq_prokaryote -> physeq_obj_16s_uparse_R1_filt
otu_table(physeq_obj_16s_uparse_R1_filt)[otu_table(physeq_obj_16s_uparse_R1_filt) <= 4] <- 0 ### tag switching
otu_table(physeq_obj_16s_uparse_R1_filt) <- otu_table(physeq_obj_16s_uparse_R1_filt)[which(rowSums(otu_table(physeq_obj_16s_uparse_R1_filt)) >= 10),] ### PCR Errors 
physeq_obj_16s_uparse_R1_filt

# setting-up variables order --------------------------------------------------------------------------
colnames(sample_data(physeq_obj_16s_uparse_R1_filt)) <- c("BarcodeSequence","LinkerPrimerSequence","Plant")
sample_data(physeq_obj_16s_uparse_R1_filt)


# add a avraible with plant names abbreviations -------------------------------------------------------
library("vegan")

ordered_ceps_16s <- c("Lycocern", "Lycofast","Phlevari","Lycodiff","Lycolate",
                      "Lycovolu","Lycodeut","Lycoscar", "Selakrau","Hupeaust")

sample_data(physeq_obj_16s_uparse_R1_filt)$Cepname <- ordered_ceps_16s
sample_data(physeq_obj_16s_uparse_R1_filt)

sample_data(physeq_obj_16s_uparse_R1_filt)$Plant <- factor(sample_data(physeq_obj_16s_uparse_R1_filt)$Plant,
                                                           levels=c("Lycopodium_deuterodensum",
                                                                    "Lycopodium_scariosum",
                                                                    "Lycopodium_volubile",
                                                                    "Lycopodiella_cernua",
                                                                    "Lycopodiella_diffusa",
                                                                    "Lycopodiella_fastigiatum",
                                                                    "Lycopodiella_lateralis",
                                                                    "Huperzia_australiana",
                                                                    "Phlegmariurus_varius",
                                                                    "Selaginella_kraussiana")) 

sample_data(physeq_obj_16s_uparse_R1_filt)$Cepname <- factor(sample_data(physeq_obj_16s_uparse_R1_filt)$Cepname,
                                                             levels=c("Lycodeut",
                                                                      "Lycoscar",
                                                                      "Lycovolu",
                                                                      "Lycocern",
                                                                      "Lycodiff",
                                                                      "Lycofast",
                                                                      "Lycolate",
                                                                      "Hupeaust",
                                                                      "Phlevari",
                                                                      "Selakrau"))

# >>> RAREFACTION CURVES ------------------------------------------------------------------------------
# converting to table objetcs -------------------------------------------------------------------------
otu_bact_16s_uparse_R1 <- as.data.frame(otu_table(physeq_obj_16s_uparse_R1_filt))
taxa_bact_16s_uparse_R1 <- as.data.frame(as.matrix(tax_table(physeq_obj_16s_uparse_R1_filt)))
metadata_bact_16s_uparse_R1 <- as.data.frame(as.matrix(sample_data(physeq_obj_16s_uparse_R1_filt)))
dim(otu_bact_16s_uparse_R1)
metadata_bact_16s_uparse_R1

metadata_bact_16s_uparse_R1$color <- c("blue", "orange", "cyan" ,"deeppink", "green",
                                       "red","yellow1","grey","black","brown")
metadata_bact_16s_uparse_R1


# *** FIGURE S2 - rarefaction curves ------------------------------------------------------------------
rarecurve(t(otu_bact_16s_uparse_R1), col = metadata_bact_16s_uparse_R1$color, label = FALSE, 
          sample=min(colSums(otu_bact_16s_uparse_R1)), step = 1000, 
          main="Prokaryotes - 16S", ylab = "Number of OTUs", xlab = "Number of DNA reads") -> rare_bact_uparse_R1

legend("bottomright", legend=metadata_bact_16s_uparse_R1$Cepname,
       col=metadata_bact_16s_uparse_R1$color, lty=1, cex=0.8, box.lty=1) #box.lty=0 remove the legend border

# >>> DATA RAREFACTION --------------------------------------------------------------------------------
set.seed(2018)

min(sample_sums(physeq_obj_16s_uparse_R1_filt))
data.frame(colSums(otu_table(physeq_obj_16s_uparse_R1_filt)))

physeq_bact_uparse_R1_ev = rarefy_even_depth(physeq_obj_16s_uparse_R1_filt, sample.size = min(sample_sums(physeq_obj_16s_uparse_R1_filt)), 
                                             rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE) 
otu_table(physeq_bact_uparse_R1_ev) <- otu_table(physeq_bact_uparse_R1_ev)[which(rowSums(otu_table(physeq_bact_uparse_R1_ev)) >= 1),]
physeq_bact_uparse_R1_ev
sample_data(physeq_bact_uparse_R1_ev)

colSums(otu_table(physeq_bact_uparse_R1_ev))
any(taxa_sums(physeq_bact_uparse_R1_ev) == 0)

# >>> TAXONOMIC COMPOSITION ---------------------------------------------------------------------------
library("dplyr")

otu_bact_uparse_R1_phylum <- physeq_bact_uparse_R1_ev %>%
  tax_glom(taxrank = "Class") %>%                     
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>%                                         
  #filter(Abundance > 0.01) %>%                         
  arrange(Class)                                     

otu_bact_uparse_R1_phylum
head(otu_bact_uparse_R1_phylum)
unique(otu_bact_uparse_R1_phylum$Phylum)
unique(otu_bact_uparse_R1_phylum$Class)

# *** FIGURE 1B - bacteria barplot --------------------------------------------------------------------
library("ggpubr")

levels(otu_bact_uparse_R1_phylum$Class)
levels(otu_bact_uparse_R1_phylum$Class)[29] <- "Unclassified"

barplot_bact_uparse_R1 = ggplot(otu_bact_uparse_R1_phylum, aes(x = Cepname, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  coord_flip() +
  labs(title="Bacteria 16S", x="", y = "Relative Abundance") +
  scale_fill_manual(values = palette_bact) +
  #facet_grid(~Plant, scales = "free_x", space="free_x") +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
  theme(strip.text.x = element_text(size = 8, face = "bold")) +
  theme(axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 0.5)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  #theme(axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")

barplot_bact_uparse_R1

# >>> ABUNDANCES --------------------------------------------------------------------------------------
physeq_bact_phylum <- tax_glom(physeq_bact_uparse_R1_ev, "Phylum")
otu_bact_abund_phy <- taxa_sums(physeq_bact_phylum)/sum(taxa_sums(physeq_bact_phylum))*100
tax_abund_bact_phy <- as(tax_table(physeq_bact_phylum), "matrix")
tax_abund_bact_phy <- as.data.frame(tax_abund_bact_phy)
tax_abund_bact_phy <- tax_abund_bact_phy[c(2)]
tax_abund_bact_phy$abundance <- as.vector(otu_bact_abund_phy)
tax_abund_bact_phy <- tax_abund_bact_phy[order(tax_abund_bact_phy$abundance, decreasing = TRUE),] 
tax_abund_bact_phy

physeq_bact_class = tax_glom(physeq_bact_uparse_R1_ev, "Class")
otu_bact_abund_class = taxa_sums(physeq_bact_class)/sum(taxa_sums(physeq_bact_class))*100
tax_abund_bact_class <- as(tax_table(physeq_bact_class), "matrix")
tax_abund_bact_class <- as.data.frame(tax_abund_bact_class)
tax_abund_bact_class <- tax_abund_bact_class[c(1:3)]
tax_abund_bact_class$abundance <- as.vector(otu_bact_abund_class)
tax_abund_bact_class <- tax_abund_bact_class[order(tax_abund_bact_class$abundance, decreasing = TRUE),] 
tax_abund_bact_class

# >>> COMAPRING DISTANCE MATRICES ---------------------------------------------------------------------
otu_bact_ev <- as.data.frame(otu_table(physeq_bact_uparse_R1_ev))
taxa_bact_ev <- as.data.frame(as.matrix(tax_table(physeq_bact_uparse_R1_ev)))
metadata_bact_ev <- as.data.frame(as.matrix(sample_data(physeq_bact_uparse_R1_ev)))
identical(colnames(otu_bact_ev), rownames(metadata_bact_ev))

colnames(otu_bact_ev) <- metadata_bact_ev$Cepname
head(otu_bact_ev)

# creating bray-curtis distance dissimilarity ---------------------------------------------------------
vegan::vegdist(t(otu_bact_ev), method="bray") -> dist_bact
dist_bact <- as.dist(dist_bact)
str(dist_bact)
dist_bact

# >>> CLUSTER ANALYSIS --------------------------------------------------------------------------------
bact_clust <- hclust(dist_bact,"complete")
plot(bact_clust)

bact_dendro <- as.dendrogram(bact_clust)
plot(bact_dendro)

# mantel test between the two distances ---------------------------------------------------------------
vegan::mantel(dist_bact, dist_91_reordered, method="spearman", permutations=9999)


# *** FIGURE 4B - bacteria tanglegram -----------------------------------------------------------------
dend_list_b <- dendlist(as.dendrogram(plant_tree_dendro), bact_dendro)
dend_list_b

dend_list_bunt<- dend_list_b %>% untangle(method = "random") 
dend_list_bunt

dend_list_bunt %>% tanglegram(common_subtrees_color_lines=TRUE,
             common_subtrees_color_branches = TRUE,
             highlight_distinct_edges=FALSE, 
             highlight_branches_lwd=FALSE,
             main_left = "Plant\nPhylogeny",
             main_right = "Community\nDissimilarity",
             main="Prokaryotes",
             sub = paste("entanglement =", round(entanglement(dend_list_bunt, 
                                L=1, leaves_matching_method = "labels"), 2)),
             hang = FALSE, 
             type = "r",
             margin_inner=4.5,
             lwd=2,
             lab.cex = 1,
             cex_main = 1.2,
             cex_sub=1,
             center = FALSE)

# calculating entanglement ----------------------------------------------------------------------------
entanglement(dend_list_bunt, L=1, leaves_matching_method = "labels")

# UPSET R ---------------------------------------------------------------------------------------------
library("UpSetR")
library("dplyr")

head(otu_bact_ev)

otu_bact_ev -> otu_bact_ev_pa
otu_bact_ev_pa[otu_bact_ev_pa>0] <- 1
head(otu_bact_ev_pa)
dim(otu_bact_ev_pa)

rownames(metadata_bact_ev) <- metadata_bact_ev$Cepname
metadata_bact_ev

sample_order_bact <- match(colnames(otu_bact_ev_pa), rownames(metadata_bact_ev))
otu_bact_ev_pa <- otu_bact_ev_pa[, sample_order_bact]
identical(colnames(otu_bact_ev_pa), rownames(metadata_bact_ev)) # now TRUE

# *** FIGURE S5 - upset plot bacteria ----------------------------------------------------------------
upset(otu_bact_ev_pa, 
      sets=metadata_bact_ev$Cepname,
      main.bar.color = "black", 
      #number.angles = 90, 
      point.size = 2, 
      line.size = 1, 
      mainbar.y.label = "OTU number", 
      sets.x.label = "OTU number", 
      text.scale = c(1.3, 1.3, 1, 1, 1, 0.75),
      empty.intersections = "on") -> upset_bact

upset_bact


#>>> BETA DIVERSITY ------------------------------------------------------------------------------
library("ggrepel")

sample_data(physeq_fungi_uparse_R1_ev)
physeq_fungi_uparse_R1_ev

sample_data(physeq_bact_uparse_R1_ev)
physeq_bact_uparse_R1_ev

pcoa_fungi = phyloseq::ordinate(physeq_fungi_uparse_R1_ev, method ="PCoA", distance="bray")

plot_pcoa_fungi = plot_ordination(physeq_fungi_uparse_R1_ev, pcoa_fungi, shape = "Cepname", 
                                  title = "PCoA Fungi") + 
  theme_classic() +
  geom_point(size=2, alpha=1, color="black") + 
  scale_shape_manual(values=c(0:10)) +
  geom_text_repel(aes(label=sample_data(physeq_fungi_uparse_R1_ev)$Cepname), size = 3) + 
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 0.5)) +
  theme(axis.text.y = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 10, face = "bold")) + 
  theme(legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8)) +
  theme(legend.position="bottom") 

plot_pcoa_fungi


pcoa_bact = phyloseq::ordinate(physeq_bact_uparse_R1_ev, method ="PCoA", distance="bray")

plot_pcoa_bact = plot_ordination(physeq_bact_uparse_R1_ev, pcoa_bact, shape = "Cepname",
                                 title = "PCoA Bacteria") + 
  theme_classic() +
  geom_point(size=2, alpha=1, color="black") + 
  scale_shape_manual(values=c(0:10)) +
  geom_text_repel(aes(label=sample_data(physeq_bact_uparse_R1_ev)$Cepname), size = 3, force = 5) + 
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
  theme(axis.text.x = element_text(angle = 0, size = 10, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 10, face = "bold")) + 
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8)) +
  theme(legend.position="bottom") 

plot_pcoa_bact

# *** FIGURE S5 - ordinations ---------------------------------------------------------------------
library("ggpubr")

ggarrange(plot_pcoa_bact,
          plot_pcoa_fungi,
          labels = c("A", "B"),
          widths = c(1,1),
          align = "none", 
          ncol = 2, nrow = 1,
          common.legend = TRUE,
          legend = c("none"))

# ***********************************************************************************************---
# ADDITIONAL CODE for REVISION ----------------------------------------------------------------------

# >>> TAXONOMIC COMPOSITION FOR EACH MISEQ DATASET --------------------------------------------------

# ITS dataset ---------------------------------------------------------------------------------------
sum(taxa_sums(physeq_obj_ITS_uparse_R1))
physeq_obj_ITS_uparse_R1

sort(unique(as.data.frame(as.matrix(tax_table(physeq_obj_ITS_uparse_R1)))$Phylum))
sort(unique(as.data.frame(as.matrix(tax_table(physeq_obj_ITS_uparse_R1)))$Class))
sort(unique(as.data.frame(as.matrix(tax_table(physeq_obj_ITS_uparse_R1)))$Order)) # Endogonales 

physeq_ITS_Endogonales <- subset_taxa(physeq_obj_ITS_uparse_R1, Order == "Endogonales")
physeq_ITS_Endogonales
sum(taxa_sums(physeq_ITS_Endogonales))
sum(taxa_sums(physeq_ITS_Endogonales))*100/sum(taxa_sums(physeq_obj_ITS_uparse_R1))

physeq_ITS_Basidio <- subset_taxa(physeq_obj_ITS_uparse_R1, Phylum == "Basidiomycota")
physeq_ITS_Basidio
sum(taxa_sums(physeq_ITS_Basidio))

physeq_ITS_Asco <- subset_taxa(physeq_obj_ITS_uparse_R1, Phylum == "Ascomycota")
physeq_ITS_Asco
sum(taxa_sums(physeq_ITS_Asco))

physeq_ITS_Glomero <- subset_taxa(physeq_obj_ITS_uparse_R1, Phylum == "Glomeromycota")
physeq_ITS_Glomero
sum(taxa_sums(physeq_ITS_Glomero))

# SSU dataset ---------------------------------------------------------------------------------------
sum(taxa_sums(physeq_obj_SSU_uparse_R1))
physeq_obj_SSU_uparse_R1

sort(unique(as.data.frame(as.matrix(tax_table(physeq_obj_SSU_uparse_R1)))$Phylum))
sort(unique(as.data.frame(as.matrix(tax_table(physeq_obj_SSU_uparse_R1)))$Class))
sort(unique(as.data.frame(as.matrix(tax_table(physeq_obj_SSU_uparse_R1)))$Order)) # Endogonales 

physeq_SSU_Basidio <- subset_taxa(physeq_obj_SSU_uparse_R1, Phylum == "Basidiomycota")
physeq_SSU_Basidio
sum(taxa_sums(physeq_SSU_Basidio))

physeq_SSU_Asco <- subset_taxa(physeq_obj_SSU_uparse_R1, Phylum == "Ascomycota")
physeq_SSU_Asco
sum(taxa_sums(physeq_SSU_Asco))

physeq_SSU_Glomero <- subset_taxa(physeq_obj_SSU_uparse_R1, Phylum == "Glomeromycota")
physeq_SSU_Glomero
sum(taxa_sums(physeq_SSU_Glomero))


# 18S dataset ---------------------------------------------------------------------------------------
sum(taxa_sums(physeq_obj_18S_uparse_R1))
physeq_obj_18S_uparse_R1

sort(unique(as.data.frame(as.matrix(tax_table(physeq_obj_18S_uparse_R1)))$Phylum))
sort(unique(as.data.frame(as.matrix(tax_table(physeq_obj_18S_uparse_R1)))$Class))
sort(unique(as.data.frame(as.matrix(tax_table(physeq_obj_18S_uparse_R1)))$Order)) # Endogonales 

physeq_18S_Basidio <- subset_taxa(physeq_obj_18S_uparse_R1, Phylum == "Basidiomycota")
physeq_18S_Basidio
sum(taxa_sums(physeq_18S_Basidio))

physeq_18S_Asco <- subset_taxa(physeq_obj_18S_uparse_R1, Phylum == "Ascomycota")
physeq_18S_Asco
sum(taxa_sums(physeq_18S_Asco))

physeq_18S_Glomero <- subset_taxa(physeq_obj_18S_uparse_R1, Phylum == "Glomeromycota")
physeq_18S_Glomero
sum(taxa_sums(physeq_18S_Glomero))

physeq_18S_Endogonales <- subset_taxa(physeq_obj_18S_uparse_R1, Order == "Endogonales")
physeq_18S_Endogonales
sum(taxa_sums(physeq_18S_Endogonales))
sum(taxa_sums(physeq_18S_Endogonales))*100/sum(taxa_sums(physeq_obj_18S_uparse_R1))

# 16S dataset ---------------------------------------------------------------------------------------
sum(taxa_sums(physeq_obj_16s_uparse))
physeq_obj_16s_uparse

otu_plastid_RDP <- subset_taxa(physeq_obj_16s_uparse, Phylum == "Cyanobacteria/Chloroplast")
taxa_names(otu_plastid_RDP)

otu_plastid_SILVA <- rownames(taxonomy_16s_uparse_SILVA[taxonomy_16s_uparse_SILVA$Class=="c:Chloroplast", ])
otu_plastid_SILVA

otu_plastid <- unique(c(taxa_names(otu_plastid_RDP), otu_plastid_SILVA))
otu_plastid

physeq_prokaryote_plastid <- otu_table(physeq_obj_16s_uparse)[otu_plastid, ]
sum(taxa_sums(physeq_prokaryote_plastid))
physeq_prokaryote_plastid
rowSums(physeq_prokaryote_plastid)

otu_mito <- rownames(taxonomy_16s_uparse_SILVA[taxonomy_16s_uparse_SILVA$Family=="f:Mitochondria", ])
otu_mito

physeq_prokaryote_mito <- otu_table(physeq_obj_16s_uparse)[otu_mito, ]
sum(taxa_sums(physeq_prokaryote_mito))
physeq_prokaryote_mito
rowSums(physeq_prokaryote_mito)


otu_plastid_mito <- unique(c(otu_plastid, otu_mito))
physeq_bacteria <- otu_table(physeq_obj_16s_uparse)[!(rownames(otu_table(physeq_obj_16s_uparse))%in%otu_plastid_mito), ]
dim(physeq_bacteria)
physeq_bacteria
rowSums(physeq_bacteria)


# >>> ACCUMILATION CURVES ------------------------------------------------------------------------

library("iNEXT")

# 16S curves -------------------------------------------------------------------------------------

df_bact <- data.frame(matrix(ncol = 1, nrow = 579))
rownames_df_bact <- paste("OTU", 1:579, sep = "_")
colnames(df_bact) <- "empty"
rownames(df_bact) <- rownames_df_bact
head(df_bact)

otu_bacteria <- as.data.frame(rowSums(physeq_bacteria))
colnames(otu_bacteria) <- "Bacteria"
head(otu_bacteria)

df_bact <- merge(df_bact, otu_bacteria, by=0, all=TRUE)
rownames(df_bact) <- df_bact$Row.names
df_bact <- df_bact[, c(2,3)]
head(df_bact)


otu_plastids <- as.data.frame(rowSums(physeq_prokaryote_plastid))
colnames(otu_plastids) <- "Plastids"
head(otu_plastids)

df_bact <- merge(df_bact, otu_plastids, by=0, all=TRUE)
rownames(df_bact) <- df_bact$Row.names
df_bact <- df_bact[, c(2:4)]
head(df_bact)


otu_mitos <- as.data.frame(rowSums(physeq_prokaryote_mito))
colnames(otu_mitos) <- "Mitochondria"
head(otu_mitos)

df_bact <- merge(df_bact, otu_mitos, by=0, all=TRUE)
rownames(df_bact) <- df_bact$Row.names
df_bact <- df_bact[, c(3:5)]
head(df_bact)

df_bact[is.na(df_bact)]<- 0

inext_bact_all <- iNEXT(df_bact, q=0, datatype="abundance", endpoint = 7500000)

plot_inext_16s <- ggiNEXT(inext_bact_all, type = 1, se = FALSE) +
  theme_classic() +
  labs(title="515F/803R", x="Number of DNA reads", y="Number of OTUs") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 10, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 10)) +
  theme(legend.position="right")

plot_inext_16s 

# plot_inext_16s$layers
# plot_inext_16s$layers <- plot_inext_16s$layers[-1]

# ITS curves -----------------------------------------------------------------------------------

df_ITS <- data.frame(matrix(ncol = 1, nrow = 5679))
rownames_df_ITS <- paste("OTU", 1:5679, sep = "_")
colnames(df_ITS) <- "empty"
rownames(df_ITS) <- rownames_df_ITS
head(df_ITS)

otu_ITS <- as.data.frame(rowSums(otu_table(physeq_obj_ITS_uparse_R1)))
colnames(otu_ITS) <- "All Fungi"
head(otu_ITS)

df_ITS <- merge(df_ITS, otu_ITS, by=0, all=TRUE)
rownames(df_ITS) <- df_ITS$Row.names
df_ITS <- df_ITS[, c(2,3)]
head(df_ITS)


otu_ITS_endo <- as.data.frame(rowSums(otu_table(physeq_ITS_Endogonales)))
colnames(otu_ITS_endo) <- "Endogonales"
head(otu_ITS_endo)

df_ITS <- merge(df_ITS, otu_ITS_endo, by=0, all=TRUE)
rownames(df_ITS) <- df_ITS$Row.names
df_ITS <- df_ITS[, c(2:4)]
head(df_ITS)


otu_ITS_basid <- as.data.frame(rowSums(otu_table(physeq_ITS_Basidio)))
colnames(otu_ITS_basid) <- "Basidiomycota"
head(otu_ITS_basid)

df_ITS <- merge(df_ITS, otu_ITS_basid, by=0, all=TRUE)
rownames(df_ITS) <- df_ITS$Row.names
df_ITS <- df_ITS[, c(3:5)]
head(df_ITS)


otu_ITS_asco <- as.data.frame(rowSums(otu_table(physeq_ITS_Asco)))
colnames(otu_ITS_asco) <- "Ascomycota"
head(otu_ITS_asco)

df_ITS <- merge(df_ITS, otu_ITS_asco, by=0, all=TRUE)
rownames(df_ITS) <- df_ITS$Row.names
df_ITS <- df_ITS[, c(2:5)]
head(df_ITS)


otu_ITS_glom <- as.data.frame(rowSums(otu_table(physeq_ITS_Glomero)))
colnames(otu_ITS_glom) <- "Glomeromycotina"
head(otu_ITS_glom)

df_ITS <- merge(df_ITS, otu_ITS_glom, by=0, all=TRUE)
rownames(df_ITS) <- df_ITS$Row.names
df_ITS <- df_ITS[, c(2:6)]
head(df_ITS)

df_ITS[is.na(df_ITS)]<- 0

inext_ITS_all <- iNEXT(df_ITS, q=0, datatype="abundance", endpoint =10000000 )

plot_inext_ITS <- ggiNEXT(inext_ITS_all, type = 1, se = FALSE) +
  theme_classic() +
  scale_y_continuous(limits = c(0,6000), breaks=c(0, 2000, 4000, 6000)) +
  labs(title="ITS1F/ITS4", x="Number of DNA reads", y="Number of OTUs") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 10, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 10)) +
  theme(legend.position="right")

plot_inext_ITS 

# 18S curves Endo18S-1F-NS6_f4 ----------------------------------------------

df_18S <- data.frame(matrix(ncol = 1, nrow = 432))
rownames_df_18S <- paste("OTU", 1:432, sep = "_")
colnames(df_18S) <- "empty"
rownames(df_18S) <- rownames_df_18S
head(df_18S)

otu_18S <- as.data.frame(rowSums(otu_table(physeq_obj_18S_uparse_R1)))
colnames(otu_18S) <- "Endo18S"
head(otu_18S)

df_18S <- merge(df_18S, otu_18S, by=0, all=TRUE)
rownames(df_18S) <- df_18S$Row.names
df_18S <- df_18S[, c(2,3)]
head(df_18S)


otu_18S_endo <- as.data.frame(rowSums(otu_table(physeq_18S_Endogonales)))
colnames(otu_18S_endo) <- "Endogonales"
head(otu_18S_endo)

df_18S <- merge(df_18S, otu_18S_endo, by=0, all=TRUE)
rownames(df_18S) <- df_18S$Row.names
df_18S <- df_18S[, c(2:4)]
head(df_18S)


otu_18S_basid <- as.data.frame(rowSums(otu_table(physeq_18S_Basidio)))
colnames(otu_18S_basid) <- "Basidiomycota"
head(otu_18S_basid)

df_18S <- merge(df_18S, otu_18S_basid, by=0, all=TRUE)
rownames(df_18S) <- df_18S$Row.names
df_18S <- df_18S[, c(3:5)]
head(df_18S)


otu_18S_asco <- as.data.frame(rowSums(otu_table(physeq_18S_Asco)))
colnames(otu_18S_asco) <- "Ascomycota"
head(otu_18S_asco)

df_18S <- merge(df_18S, otu_18S_asco, by=0, all=TRUE)
rownames(df_18S) <- df_18S$Row.names
df_18S <- df_18S[, c(2:5)]
head(df_18S)


otu_18S_glom <- as.data.frame(rowSums(otu_table(physeq_18S_Glomero)))
colnames(otu_18S_glom) <- "Glomeromycotina"
head(otu_18S_glom)

df_18S <- merge(df_18S, otu_18S_glom, by=0, all=TRUE)
rownames(df_18S) <- df_18S$Row.names
df_18S <- df_18S[, c(2:6)]
head(df_18S)

df_18S[is.na(df_18S)]<- 0

inext_18S_all <- iNEXT(df_18S, q=0, datatype="abundance", endpoint = 80000)

plot_inext_18S <- ggiNEXT(inext_18S_all, type = 1, se = FALSE) +
  theme_classic() +
  labs(title="Endo18S-1F/NS6_f4", x="Number of DNA reads", y="Number of OTUs") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 10, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 10)) +
  theme(legend.position="right")

plot_inext_18S 

# SSU curves NS31_f3-AML2_f3 ----------------------------------------------

df_SSU <- data.frame(matrix(ncol = 1, nrow = 248))
rownames_df_SSU <- paste("OTU", 1:248, sep = "_")
colnames(df_SSU) <- "empty"
rownames(df_SSU) <- rownames_df_SSU
head(df_SSU)

otu_SSU <- as.data.frame(rowSums(otu_table(physeq_obj_SSU_uparse_R1)))
colnames(otu_SSU) <- "18S"
head(otu_SSU)

df_SSU <- merge(df_SSU, otu_SSU, by=0, all=TRUE)
rownames(df_SSU) <- df_SSU$Row.names
df_SSU <- df_SSU[, c(2,3)]
head(df_SSU)


otu_SSU_basid <- as.data.frame(rowSums(otu_table(physeq_SSU_Basidio)))
colnames(otu_SSU_basid) <- "Basidiomycota"
head(otu_SSU_basid)

df_SSU <- merge(df_SSU, otu_SSU_basid, by=0, all=TRUE)
rownames(df_SSU) <- df_SSU$Row.names
df_SSU <- df_SSU[, c(3:4)]
head(df_SSU)


otu_SSU_asco <- as.data.frame(rowSums(otu_table(physeq_SSU_Asco)))
colnames(otu_SSU_asco) <- "Ascomycota"
head(otu_SSU_asco)

df_SSU <- merge(df_SSU, otu_SSU_asco, by=0, all=TRUE)
rownames(df_SSU) <- df_SSU$Row.names
df_SSU <- df_SSU[, c(2:4)]
head(df_SSU)


otu_SSU_glom <- as.data.frame(rowSums(otu_table(physeq_SSU_Glomero)))
colnames(otu_SSU_glom) <- "Glomeromycotina"
head(otu_SSU_glom)

df_SSU <- merge(df_SSU, otu_SSU_glom, by=0, all=TRUE)
rownames(df_SSU) <- df_SSU$Row.names
df_SSU <- df_SSU[, c(2:5)]
head(df_SSU)

df_SSU[is.na(df_SSU)]<- 0

inext_SSU_all <- iNEXT(df_SSU, q=0, datatype="abundance", endpoint = 1500)

plot_inext_SSU <- ggiNEXT(inext_SSU_all, type = 1, se = FALSE) +
  theme_classic() +
  scale_y_continuous(limits = c(0,10), breaks=c(0, 2, 4, 6, 8, 10, 12)) +
  labs(title="NS31_f3/AML2_f3", x="Number of DNA reads", y="Number of OTUs") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, size = 10, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 10, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 10)) +
  theme(legend.position="right")

plot_inext_SSU 

# *** FIGURE S2 accumulation curves -------------------------------------------------------
ggarrange(plot_inext_16s,plot_inext_ITS,
          plot_inext_SSU,plot_inext_18S,
          labels = c("A","B","C","D"),
          widths = c(1,1,1,1),
          align = "none", 
          ncol = 2, nrow = 2,
          common.legend = FALSE,
          legend = c("right"))


