
---
  title: "Scripts to Analyse Microbial Community Data"
author: "Gian MN Benucci, Ph.D."
date: "January 11th, 2017"
output:
  pdf_document:
  latex_engine: xelatex
html_document: default
---
  
  
  ```{r, echo=FALSE, warning=FALSE, message=FALSE}
```

'''RANDOM FOREST https://rpubs.com/michberr/randomforestmicrobe
https://www.r-bloggers.com/random-forest-classification-of-mushrooms/
PHYLOSEQ function https://github.com/joey711/phyloseq/issues/442
DESEQ DIFF ABUND http://joey711.github.io/phyloseq-extensions/DESeq2.html
CO-OCCURRENCE https://cran.r-project.org/web/packages/EcoSimR/vignettes/CoOccurrenceVignette.html
RANDOM FOREST https://rpubs.com/michberr/randomforestmicrobe
Plotting READS https://rpubs.com/marschmi/133626
BIOENV and ENVFIT http://menugget.blogspot.com/2011/06/clarke-and-ainsworths-bioenv-and-bvstep.html
https://github.com/joey711/phyloseq/issues/105
http://stackoverflow.com/questions/14711470/plotting-envfit-vectors-vegan-package-in-ggplot2/25425258#25425258
Rmarkdown for RStudio http://rmarkdown.rstudio.com/html_document_format.html
RAM package https://cran.r-project.org/web/packages/RAM/RAM.pdf
http://rpackages.ianhowson.com/cran/RAM/
http://fuzzysim.r-forge.r-project.org/fuzzySim-tutorial.html

'''


### PRELIMINARY SETUP ------------------------------------------------------------------------------------------------------

# setting the work environment 
rm(list = ls(all=TRUE)) # removes all variables in the global environment so you start fresh
Sys.time() # prints out the time and date you ran the code
options(scipen = 999) #to use decimals
options(max.print=100000000) # to print more lines on the display
set.seed(1977) #to make reproduceble results
detach(package:phyloseq, unload=TRUE) #to reorder packages
search() #to search into the environment
rm(list= ls()[!(ls() %in% c('keepThis','andThis'))]) #remove all but ...

session.info() #to get session information

#shiny-phyloseq: web-based phyloseq resource for plotting 
install.packages("shiny")
shiny::runGitHub("shiny-phyloseq","joey711")

# load previous data
load(".RData")

### NECESSARY PACKAGES -----------------------------------------------------------------------------------------
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(plyr)
library(dplyr)
library(ggrepel)
library(vegan)
library(DESeq2)
library(metagenomeSeq)
library(indicspecies)
library(limma)

library(microbiome)
library(RAM)
library(BiodiversityR)
library(biom)
library(gplots)
library(RVAideMemoire)
library(lavdsv)
library(mvabund)
library(boral)

library(RColorBrewer)
library(interactiveDisplay)
library(picante)
library(grid)
library(Hmisc)
library(corrplot)
library(psych)
library(ggdendro)
library(colorRamps)
library(multtest)
library(data.table)



# install metagenomeSeq can be a trouble 
# running biocValid(fix=TRUE) sorted it out


### IMPORTING DATA into R ----------------------------------------------------------------------------------------------------------------------

# 1) imoprt from .biom 
#importing ITS biom data

#names(bs1) <- gsub("\\s.+$", "", names(bs1)) #rename otus sequence header if longer than OTU_###

biom_ITS = import_biom("otu_table_ITS.biom")
map_ITS = import_qiime_sample_data("mapping_ITS.txt")
sample_data(biom_ITS) <- map_ITS
colnames(tax_table(biom_ITS)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
otus_rep_ITS <- readDNAStringSet("otus_ITS.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
biom_ITS <- merge_phyloseq(biom_ITS, otus_rep_ITS)
biom_ITS
sample_data(biom_ITS)
head(tax_table(biom_ITS))
refseq(biom_ITS)

biom_LSU = import_biom("otu_table_LSU.biom")
map_LSU = import_qiime_sample_data("mapping_LSU.txt")
sample_data(biom_LSU) <- map_LSU
colnames(tax_table(biom_LSU)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus")
otus_rep_LSU <- readDNAStringSet("otus_LSU.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
biom_LSU <- merge_phyloseq(biom_LSU, otus_rep_LSU)
biom_LSU
sample_data(biom_LSU)
head(tax_table(biom_LSU))

biom_16s = import_biom("otu_table_16s.biom")
map_16s = import_qiime_sample_data("mapping_16s.txt")
sample_data(biom_16s) <- map_16s
colnames(tax_table(biom_16s)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
otus_rep_16s <- readDNAStringSet("otus_16s.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
biom_16s <- merge_phyloseq(biom_16s, otus_rep_16s)
biom_16s
sample_data(biom_16s)
head(tax_table(biom_16s))



### data PREFILTERING --------------------------------------------------------------------------------------------------------------------

# cleaning ITS taxonomy from extra characters
tax_table(biom_ITS)[, "Kingdom"] <- gsub("PMI", "", tax_table(biom_ITS)[, "Kingdom"])
tax_table(biom_ITS)[, "Kingdom"] <- gsub("NVP", "", tax_table(biom_ITS)[, "Kingdom"])

tax_table(biom_ITS)[, "Kingdom"] <- gsub("k__", "", tax_table(biom_ITS)[, "Kingdom"])
tax_table(biom_ITS)[, "Phylum"] <- gsub("p__", "", tax_table(biom_ITS)[, "Phylum"])
tax_table(biom_ITS)[, "Class"] <- gsub("c__", "", tax_table(biom_ITS)[, "Class"])
tax_table(biom_ITS)[, "Order"] <- gsub("o__", "", tax_table(biom_ITS)[, "Order"])
tax_table(biom_ITS)[, "Family"] <- gsub("f__", "", tax_table(biom_ITS)[, "Family"])
tax_table(biom_ITS)[, "Genus"] <- gsub("g__", "", tax_table(biom_ITS)[, "Genus"])
tax_table(biom_ITS)[, "Species"] <- gsub("s__", "", tax_table(biom_ITS)[, "Species"])

head(otu_table(biom_ITS))
head(tax_table(biom_ITS))
head(sample_data(biom_ITS))

write.csv(otu_table(biom_ITS), "otu_table_ITS_test.csv")
write.csv(tax_table(biom_ITS), "tax_table_ITS_test.csv")


tax_table(biom_LSU)[, "Kingdom"] <- gsub("k__", "", tax_table(biom_LSU)[, "Kingdom"])
tax_table(biom_LSU)[, "Phylum"] <- gsub("p__", "", tax_table(biom_LSU)[, "Phylum"])
tax_table(biom_LSU)[, "Class"] <- gsub("c__", "", tax_table(biom_LSU)[, "Class"])
tax_table(biom_LSU)[, "Order"] <- gsub("o__", "", tax_table(biom_LSU)[, "Order"])
tax_table(biom_LSU)[, "Family"] <- gsub("f__", "", tax_table(biom_LSU)[, "Family"])
tax_table(biom_LSU)[, "Genus"] <- gsub("g__", "", tax_table(biom_LSU)[, "Genus"])
tax_table(biom_LSU)[, "Species"] <- gsub("s__", "", tax_table(biom_LSU)[, "Species"])

head(otu_table(biom_LSU))
head(tax_table(biom_LSU))
head(sample_data(biom_LSU))


tax_table(biom_16s)[, "Kingdom"] <- gsub("k__", "", tax_table(biom_16s)[, "Kingdom"])
tax_table(biom_16s)[, "Phylum"] <- gsub("p__", "", tax_table(biom_16s)[, "Phylum"])
tax_table(biom_16s)[, "Class"] <- gsub("c__", "", tax_table(biom_16s)[, "Class"])
tax_table(biom_16s)[, "Order"] <- gsub("o__", "", tax_table(biom_16s)[, "Order"])
tax_table(biom_16s)[, "Family"] <- gsub("f__", "", tax_table(biom_16s)[, "Family"])
tax_table(biom_16s)[, "Genus"] <- gsub("g__", "", tax_table(biom_16s)[, "Genus"])
tax_table(biom_16s)[, "Species"] <- gsub("s__", "", tax_table(biom_16s)[, "Species"])

head(otu_table(biom_16s))
head(tax_table(biom_16s))
head(sample_data(biom_16s))


# filtering out non-fungal OTUs
biom_ITS <- subset_taxa(biom_ITS, Kingdom == "Fungi")
biom_ITS

biom_LSU <- subset_taxa(biom_LSU, Kingdom == "Fungi")
biom_LSU

# filtering out non-bacterial OTUs
# pay attention here if you want to keep Archaea!!
biom_16s <- subset_taxa(biom_16s, Kingdom == "Bacteria")

biom_16s <- subset_taxa(biom_16s, Phylum!="Chloroplast")
biom_16s <- subset_taxa(biom_16s, Class!="Chloroplast")
biom_16s <- subset_taxa(biom_16s, Order!="Chloroplast")
biom_16s <- subset_taxa(biom_16s, Family!="Chloroplast")
biom_16s <- subset_taxa(biom_16s, Genus!="Chloroplast")
biom_16s <- subset_taxa(biom_16s, Phylum!="chloroplast")
biom_16s <- subset_taxa(biom_16s, Class!="chloroplast")
biom_16s <- subset_taxa(biom_16s, Order!="chloroplast")
biom_16s <- subset_taxa(biom_16s, Family!="chloroplast")
biom_16s <- subset_taxa(biom_16s, Genus!="chloroplast")

biom_16s <- subset_taxa(biom_16s, Phylum!="Mitochondria")
biom_16s <- subset_taxa(biom_16s, Class!="Mitochondria")
biom_16s <- subset_taxa(biom_16s, Order!="Mitochondria")
biom_16s <- subset_taxa(biom_16s, Family!="Mitochondria")
biom_16s <- subset_taxa(biom_16s, Genus!="Mitochondria")
biom_16s <- subset_taxa(biom_16s, Phylum!="mitochondria")
biom_16s <- subset_taxa(biom_16s, Class!="mitochondria")
biom_16s <- subset_taxa(biom_16s, Order!="mitochondria")
biom_16s <- subset_taxa(biom_16s, Family!="mitochondria")
biom_16s <- subset_taxa(biom_16s, Genus!="mitochondria")

biom_16s

# manual correction of the tax_table - OPTIONAL!
write.csv(tax_table(biom_ITS) , file = "tax_tab_biom_ITS.csv")
tax_tab_biom_ITS <- read.csv("tax_tab_biom_ITS.csv", header=T, row.names =1)
tax_table(biom_ITS) <- tax_table(as.matrix(tax_tab_biom_ITS))


### PRELIMINARY DATA VISUALIZATION ----------------------------------------------------------------------------------------------------------------------------------

# check the ITS sequencing depth of each sample 
library(plyr)

sums_biom_ITS <- data.frame(colSums(otu_table(biom_ITS)))
colnames(sums_biom_ITS) <- "Sample_TotalSeqs"
sums_biom_ITS$sample <- row.names(sums_biom_ITS)
sums_biom_ITS$Description <- sample_data(biom_ITS)$Description
sums_biom_ITS$Site2 <- sample_data(biom_ITS)$Site2 # if you want to add another variable to use further
sums_biom_ITS <- arrange(sums_biom_ITS, Sample_TotalSeqs)
#sums_biom_ITS <- arrange(sums_biom_ITS, Description) # to order according another variable in the mapping
sums_biom_ITS

write.csv(sums_biom_ITS, file = "library_ITS.csv")

# create a plot of the number of ITS sequences per sample
library(ggplot2)
sums_biom_ITS <- sums_biom_ITS[order(sums_biom_ITS$Description), ]
sums_biom_ITS$Description <- factor(sums_biom_ITS$Description, levels = sums_biom_ITS$Description[order(sums_biom_ITS$Description)])
sums_biom_ITS$Description

ggplot(sums_biom_ITS, aes(x=reorder(Description, Sample_TotalSeqs), y = Sample_TotalSeqs)) + 
  ylab("Number of Sequences per Sample") +
  geom_bar(stat = "identity", colour="black",fill="cornflowerblue")  + xlab("Sample Name") + 
  ggtitle("Total Number of Sequences per Sample ITS") + 
  theme(axis.text.x = element_text(colour = "black", size=10, angle=90, hjust = 1, vjust = 1)) 

# plot boxplots to see outliers rep according an a priori variable 
ggplot(sums_biom_ITS, aes(x=reorder(Description, Sample_TotalSeqs), y = Sample_TotalSeqs)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=16, outlier.size=4) +
  geom_jitter(shape=2, position=position_jitter(0.2), aes(colour = Description)) + 
  geom_text(aes(label=sample), size = 3) +
  theme(axis.text.x=element_text(angle=90, hjust=1))

# adding samples name and labels
library("ggrepel")

ggplot(sums_biom_ITS, aes(x=reorder(Site2, Sample_TotalSeqs), y = Sample_TotalSeqs)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=16, outlier.size=4) +
  geom_point() + geom_text_repel(aes(label=sample), size = 3) + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  labs(x = "Sample Replicate", y = "Read Number") +
  ggtitle("Sample distribution ITS")


# create a plot of the number of LSU sequences per sample
sums_biom_LSU <- data.frame(colSums(otu_table(biom_LSU)))
colnames(sums_biom_LSU) <- "Sample_TotalSeqs"
sums_biom_LSU$sample <- row.names(sums_biom_LSU)
sums_biom_LSU$Description <- sample_data(biom_LSU)$Description
sums_biom_LSU$Site2 <- sample_data(biom_LSU)$Site2 # if you want to add another variable to use further
sums_biom_LSU <- arrange(sums_biom_LSU, Sample_TotalSeqs)
#sums_biom_LSU <- arrange(sums_biom_LSU, Description) # to order according another variable in the mapping
sums_biom_LSU


library(ggplot2)
sums_biom_LSU <- sums_biom_LSU[order(sums_biom_LSU$Description), ]
sums_biom_LSU$Description <- factor(sums_biom_LSU$Description, levels = sums_biom_LSU$Description[order(sums_biom_LSU$Description)])
sums_biom_LSU$Description

ggplot(sums_biom_LSU, aes(x=reorder(Description, Sample_TotalSeqs), y = Sample_TotalSeqs)) + 
  ylab("Number of Sequences per Sample") +
  geom_bar(stat = "identity", colour="black",fill="cornflowerblue")  + xlab("Sample Name") + 
  ggtitle("Total Number of Sequences per Sample LSU") + 
  theme(axis.text.x = element_text(colour = "black", size=10, angle=90, hjust = 1, vjust = 1)) 

# plot boxplots to see outliers rep according an a priori variable 
ggplot(sums_biom_LSU, aes(x=reorder(Description, Sample_TotalSeqs), y = Sample_TotalSeqs)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=16, outlier.size=4) +
  geom_jitter(shape=2, position=position_jitter(0.2), aes(colour = Description)) + 
  geom_text(aes(label=sample), size = 3) +
  theme(axis.text.x=element_text(angle=90, hjust=1))

# adding samples name and labels
library("ggrepel")

ggplot(sums_biom_LSU, aes(x=reorder(Site2, Sample_TotalSeqs), y = Sample_TotalSeqs)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=16, outlier.size=4) +
  geom_point() + geom_text_repel(aes(label=sample), size = 3) + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  labs(x = "Sample Replicate", y = "Read Number") +
  ggtitle("Sample distribution LSU")



# check the 16S sequencing depth of each sample 
sums_biom_16s <- data.frame(colSums(otu_table(biom_16s)))
colnames(sums_biom_16s) <- "Sample_TotalSeqs"
sums_biom_16s$sample <- row.names(sums_biom_16s)
sums_biom_16s$Description <- sample_data(biom_16s)$Description
sums_biom_16s$Site2 <- sample_data(biom_16s)$Site2
sums_biom_16s <- arrange(sums_biom_16s, Sample_TotalSeqs)
sums_biom_16s

write.csv(sums_biom_16s, file = "library_16s.csv")

# create a plot of the number of 16s sequences per sample
ggplot(sums_biom_16s, aes(x=reorder(Description, Sample_TotalSeqs), y = Sample_TotalSeqs)) +
  ylab("Number of Sequences per Sample") +
  geom_bar(stat = "identity", colour="black",fill="cornflowerblue")  + xlab("Sample Name") + 
  ggtitle("Total Number of Sequences per Sample 16S") + 
  theme(axis.text.x = element_text(colour = "black", size=10, angle=90, hjust = 1, vjust = 1))

# plot boxplots to see outliers rep
ggplot(sums_biom_16s, aes(x=reorder(Description, Sample_TotalSeqs), y = Sample_TotalSeqs)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=16, outlier.size=4) +
  geom_jitter(shape=2, position=position_jitter(0.2), aes(colour = Description)) + 
  geom_text(aes(label=sample), size = 3) +
  theme(axis.text.x = element_text(colour = "black", size=10, angle=90, hjust = 1, vjust = 1))

# adding samples name and labels
ggplot(sums_biom_16s, aes(x=reorder(Site2, Sample_TotalSeqs), y = Sample_TotalSeqs)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=16, outlier.size=4) +
  geom_point() + geom_text_repel(aes(label=sample), size = 3) +
  theme(axis.text.x = element_text(colour = "black", size=10, angle=90, hjust = 1, vjust = 1)) +
  labs(x = "Sample Replicate", y = "Read Number") +
  ggtitle("Sample distribution 16S")  


# check distribution of sampling depth 

ggplot(sums_biom_ITS, aes(x = Sample_TotalSeqs)) + # Histogram of sample read counts
  geom_histogram(color = "black", fill = "indianred", binwidth = 25) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) +
  scale_color_manual(values="green")

ggplot(sums_biom_16s, aes(x = Sample_TotalSeqs)) + # Histogram of sample read counts
  geom_histogram(color = "black", fill = "indianred", binwidth = 25) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) +
  scale_color_manual(values="green")


# check also here for library exploration
# http://joey711.github.io/phyloseq-demo/Restroom-Biogeography.html

# check OTUs included in this dataset that have no counted reads
any(taxa_sums(biom_ITS) == 0)

readsumsdf_ITS = data.frame(nreads = sort(taxa_sums(biom_ITS), TRUE), sorted = 1:ntaxa(biom_ITS), type = "OTUs")
readsumsdf_ITS = rbind(readsumsdf_ITS, data.frame(nreads = sort(sample_sums(biom_ITS), TRUE), sorted = 1:nsamples(biom_ITS), type = "Samples"))

title = "Total number of reads"
head(readsumsdf_ITS)

p = ggplot(readsumsdf_ITS, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

set.seed(2) # set seed 
# rarefy at even depth of 1000 reads - best to find a random number see below!
biom_ITS_ev = rarefy_even_depth(biom_ITS, sample.size = 10000)
# And now an alternative to random subsampling, a simple proportional transformation
biom_ITS_prop = transform_sample_counts(biom_ITS, function(x) 1000 * x/sum(x))

# let's replot the sample sums of each of these new data objects,
# to convince ourselves that all of the samples now sum to 500.
par(mfrow = c(1, 2)) # don't forget to close device at the end of plotting
title = "Sum of reads for each sample, biom_ITS even depth"
plot(sort(sample_sums(biom_ITS_ev), TRUE), type = "h", main = title, ylab = "reads", 
     ylim = c(0, 1000))
title = "Sum of reads for each sample, biom_ITS proportional abundances"
plot(sort(sample_sums(biom_ITS_prop), TRUE), type = "h", main = title, ylab = "reads", 
     ylim = c(0, 1000))

dev.off() # to set par() to original default values





### SUBSETTING TO SPECIFIC SAMPLES, OTUS, or TAXA --------------------------------------------------------------------------------------
# http://joey711.github.io/phyloseq/preprocess
# https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-basics.html

# filtering samples with low sequencing depth form the otu_table
#df[, !(colnames(df) %in% c("x","bar","foo"))]

sample_data(biom_ITS)
sample_data(biom_16s)
sample_data(biom_LSU)

biom_ITS_filt <- subset_samples(biom_ITS, Keep%in%c("keep"))
otu_table(biom_ITS_filt) <- otu_table(biom_ITS_filt)[which(rowSums(otu_table(biom_ITS_filt)) > 1),] 
biom_ITS_filt

biom_LSU_filt <- subset_samples(biom_LSU, Keep%in%c("keep"))
otu_table(biom_LSU_filt) <- otu_table(biom_LSU_filt)[which(rowSums(otu_table(biom_LSU_filt)) > 1),] 
biom_LSU_filt

biom_16s_filt <- subset_samples(biom_16s, Keep%in%c("keep"))
otu_table(biom_16s_filt) <- otu_table(biom_16s_filt)[which(rowSums(otu_table(biom_16s_filt)) > 1),] 
biom_16s_filt

sample_data(biom_ITS_filt)
sample_data(biom_LSU_filt)
sample_data(biom_16s_filt)



### FILTERING AND NORMALIZING DATA ----------------------------------------------------------------------------------------------

library(dplyr)

# One approach for Fungi is to remove OTUs with less than 10 reads and less than 5 in at 
# least one sample. Please check Peter Kennedy tougths in google drive doc and 
# see Oliver et al. (2015, Fungal Ecology) and Lindahl et al. (2013, New Phyt)

#See functions all() and any(). The apply() function can be used to run functions over 
#rows or columns. (MARGIN = 1 is rows, MARGIN = 2 is columns, etc). Use df[, -1] to 
#ignore the id variable when doing the comparisons.

any(taxa_sums(biom_ITS_filt) == 0)
any(taxa_sums(biom_LSU_filt) == 0)
any(taxa_sums(biom_16s_filt) == 0)

biom_ITS_cv <- biom_ITS_filt
otu_table(biom_ITS_cv) <- otu_table(biom_ITS_cv)[rowSums(otu_table(biom_ITS_cv) > 0) >= 5, ]
otu_table(biom_ITS_cv) <- otu_table(biom_ITS_cv)[which(rowSums(otu_table(biom_ITS_cv)) >= 10),] 
biom_ITS_cv
sample_data(biom_ITS_cv)

#another apprach below
#apply(otu_table(biom_ITS_cv), 1, function(x) sum(x>0)) >= 5
#otu_table(biom_ITS_cv) <- otu_table(biom_ITS_cv)[apply(otu_table(biom_ITS_cv), 1, function(x) sum(x>0)) >= 5, ]


biom_LSU_cv <- biom_LSU_filt
otu_table(biom_LSU_cv) <- otu_table(biom_LSU_cv)[rowSums(otu_table(biom_LSU_cv) > 0) >= 5, ]
otu_table(biom_LSU_cv) <- otu_table(biom_LSU_cv)[which(rowSums(otu_table(biom_LSU_cv)) >= 10),] 
biom_LSU_cv


#apply(otu_table(biom_ITS_cv), 1, function(x) sum(x!=0))
#apply(otu_table(biom_ITS_cv), 1, function(x) sum(x>0))

biom_16s_cv <- biom_16s_filt
otu_table(biom_16s_cv) <- otu_table(biom_16s_cv)[rowSums(otu_table(biom_16s_cv) > 0) >= 5, ]
otu_table(biom_16s_cv) <- otu_table(biom_16s_cv)[which(rowSums(otu_table(biom_16s_cv)) >= 10),] 
biom_16s_cv
biom_ITS_cv

write.csv(otu_table(biom_ITS_cv), "1259otu_table_ITS.csv")
write.csv(tax_table(biom_ITS_cv), "1259tax_table_ITS.csv")
write.csv(sample_data(biom_ITS_cv), "map54_ITS.csv")
write.table(refseq(biom_ITS_cv), "1259repseq_ITS.fasta")

write.csv(otu_table(biom_LSU_cv), "1061otu_table_LSU.csv")
write.csv(tax_table(biom_LSU_cv), "1061tax_table_LSU.csv")
write.csv(sample_data(biom_LSU_cv), "map51_LSU.csv")
write.table(refseq(biom_LSU_cv), "1061repseq_LSU.fasta")

write.csv(otu_table(biom_16s_cv), "5791otu_table_16s.csv")
write.csv(tax_table(biom_16s_cv), "5791tax_table_16s.csv")
write.csv(sample_data(biom_16s_cv), "map54_16s.csv")
write.table(refseq(biom_16s_cv), "5791repseq_16s.fasta")


# Hellinger standardization
#Legendre, P. & Gallagher, E.D. (2001) Ecologically meaningful
#transformations for ordination of species data. Oecologia 129; 271–280. 

library(vegan)

biom_ITS_cv -> biom_ITS_hell
otu_table(biom_ITS_hell) <- otu_table(decostand(otu_table(biom_ITS_cv), method = "hellinger"), taxa_are_rows=TRUE)
otu_table(biom_ITS_hell)

biom_LSU_cv -> biom_LSU_hell
otu_table(biom_LSU_hell) <- otu_table(decostand(otu_table(biom_LSU_cv), method = "hellinger"), taxa_are_rows=TRUE)
otu_table(biom_LSU_hell)

biom_16s_cv -> biom_16s_hell
otu_table(biom_16s_hell) <- otu_table(decostand(otu_table(biom_16s_cv), method = "hellinger"), taxa_are_rows=TRUE)
otu_table(biom_16s_hell)

biom_ITS_ecm -> biom_ITS_ecm_hell
otu_table(biom_ITS_ecm_hell) <- otu_table(decostand(otu_table(biom_ITS_ecm), method = "hellinger"), taxa_are_rows=TRUE)
otu_table(biom_ITS_ecm_hell)


# normalizing data counts into proportional relative abundances or OTUs frequencies
biom_ITS_prop = transform_sample_counts(biom_ITS_cv, function(OTU) OTU/sum(OTU))
biom_LSU_prop = transform_sample_counts(biom_LSU_cv, function(OTU) OTU/sum(OTU))
biom_16s_prop = transform_sample_counts(biom_16s_cv, function(OTU) OTU/sum(OTU))

biom_ITS_prop -> biom_ITS_hell2
otu_table(biom_ITS_hell2) <- otu_table(decostand(otu_table(biom_ITS_prop), method = "hellinger"), taxa_are_rows=TRUE)
otu_table(biom_ITS_hell2)

biom_LSU_prop -> biom_LSU_hell2
otu_table(biom_LSU_hell2) <- otu_table(decostand(otu_table(biom_LSU_prop), method = "hellinger"), taxa_are_rows=TRUE)
otu_table(biom_LSU_hell2)






### DATA COUNTS in COMMUNITY DATA -------------------------------------------------------------------------------------------

#### >> Extracting relative abundances ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
biom_ITS_cv_phy = tax_glom(biom_ITS_cv, "Phylum")
otu_table(biom_ITS_cv_phy)
tax_table(biom_ITS_cv_phy)

otu_ITS_cv_phy = taxa_sums(biom_ITS_cv_phy)/sum(taxa_sums(biom_ITS_cv_phy))*100
as.data.frame(otu_ITS_cv_phy)
as.vector(otu_ITS_cv_phy)
tax_ITS_cv_phy <- data.frame(tax_table(biom_ITS_cv_phy))
tax_ITS_cv_phy <- tax_ITS_cv_phy[c(1:2)]
tax_ITS_cv_phy
tax_ITS_cv_phy$abundance <- as.vector(otu_ITS_cv_phy)
tax_ITS_cv_phy <- tax_ITS_cv_phy[order(tax_ITS_cv_phy$abundance, decreasing = TRUE),] 
tax_ITS_cv_phy



biom_LSU_cv_phy = tax_glom(biom_LSU_cv, "Phylum")
otu_table(biom_LSU_cv_phy)
tax_table(biom_LSU_cv_phy)

otu_LSU_cv_phy = taxa_sums(biom_LSU_cv_phy)/sum(taxa_sums(biom_LSU_cv_phy))*100
as.data.frame(otu_LSU_cv_phy)
as.vector(otu_LSU_cv_phy)
tax_LSU_cv_phy <- data.frame(tax_table(biom_LSU_cv_phy))
tax_LSU_cv_phy <- tax_LSU_cv_phy[c(1:2)]
tax_LSU_cv_phy
tax_LSU_cv_phy$abundance <- as.vector(otu_LSU_cv_phy)
tax_LSU_cv_phy <- tax_LSU_cv_phy[order(tax_LSU_cv_phy$abundance, decreasing = TRUE),] 
tax_LSU_cv_phy



biom_16s_cv_phy = tax_glom(biom_16s_cv, "Phylum")
otu_table(biom_16s_cv_phy)
tax_table(biom_16s_cv_phy)

otu_16s_cv_phy = taxa_sums(biom_16s_cv_phy)/sum(taxa_sums(biom_16s_cv_phy))*100
as.data.frame(otu_16s_cv_phy)
as.vector(otu_16s_cv_phy)

tax_16s_cv_phy <- data.frame(tax_table(biom_16s_cv_phy))
tax_16s_cv_phy <- tax_16s_cv_phy[c(1:2)]
tax_16s_cv_phy
tax_16s_cv_phy$abundance <- as.vector(otu_16s_cv_phy)
tax_16s_cv_phy <- tax_16s_cv_phy[order(tax_16s_cv_phy$abundance, decreasing = TRUE),] 
tax_16s_cv_phy



biom_16s_cv_gen = tax_glom(biom_16s_cv, "Genus")
otu_16s_cv_gen = taxa_sums(biom_16s_cv_gen)/sum(taxa_sums(biom_16s_cv_gen))*100

tax_16s_cv_gen <- data.frame(tax_table(biom_16s_cv_gen))
tax_16s_cv_gen <- tax_16s_cv_gen[c(5:6)]
tax_16s_cv_gen
tax_16s_cv_gen$abundance <- as.vector(otu_16s_cv_gen)
tax_16s_cv_gen <- tax_16s_cv_gen[order(tax_16s_cv_gen$abundance, decreasing = TRUE),] 
tax_16s_cv_gen

#####


biom_ITS_cv_gen = tax_glom(biom_ITS_cv, "Genus")
biom_ITS_cv_gen


otu_ITS_aus_gen = taxa_sums(biom_ITS_AU)/sum(taxa_sums(biom_ITS_AU))*100
otu_ITS_aus_gen
tax_ITS_aus_gen <- data.frame(tax_table(biom_ITS_aus_gen))
tax_ITS_aus_gen <- tax_ITS_aus_gen[c(1:7)]
tax_ITS_aus_gen
tax_ITS_aus_gen$abundance <- as.vector(otu_ITS_aus_gen)
tax_ITS_aus_gen <- tax_ITS_aus_gen[order(tax_ITS_aus_gen$abundance, decreasing = TRUE),] 
tax_ITS_aus_gen

## OTUs abundances and data counts fungi
biom_ITS_ab <- transform_sample_counts(biom_ITS_cv, function(x) x/sum(x))
biom_ITS_ab

otu_table(biom_ITS_ab)
taxa_sums(fungi_top100)/54 # abundance of each of the first 50 top OTUs
sum(taxa_sums(fungi_top100)/54)  # total abundance of the first 50 top OTUs
sum(taxa_sums(biom_ITS_ab)/54) # sum of total abundance must be 1
sample_sums(biom_ITS_ab) # total abundace of each sample must be 1

otus_ITS <- data.frame(tax_table(biom_ITS_ab))
otus_ITS <- otus_ITS[c(2:7)]
otus_ITS$abundance <- as.vector(taxa_sums(biom_ITS_ab)/54)
otus_ITS <- otus_ITS[order(otus_ITS$abundance, decreasing = TRUE),] 
otus_ITS[1:50,] # select first 50 rows
otus_ITS

#tax_Gen <- arrange(tax_Gen, abundance, decreasing=TRUE)
#top100ITS <- names(sort(taxa_sums(biom_ITS_ab), TRUE)[1:50]) # sampling first 50 taxa
#fungi_top100 <- prune_taxa(top100ITS, biom_ITS_ab)
#fungi_top100

write.csv(taxa_sums(fungi_top100) , file = "fungi_taxonomy.csv")
write.csv(tax_table(fungi_top100)/54, file = "fungi_top100.csv")

biom_ITS_Phy = tax_glom(biom_ITS_cv, "Phylum")
biom_ITS_Phy_ab = transform_sample_counts(biom_ITS_Phy, function(x) x/sum(x))
biom_ITS_Phy_ab
sum(taxa_sums(biom_ITS_Phy_ab)/54)
sample_sums(biom_ITS_Phy_ab)
tax_Phy <- data.frame(tax_table(biom_ITS_Phy))
tax_Phy <- tax_Phy[c(2)]
tax_Phy$abundance <- as.vector(taxa_sums(biom_ITS_Phy_ab)/54)
tax_Phy <- tax_Phy[order(tax_Phy$abundance, decreasing = TRUE),] 
#tax_Phy <- arrange(tax_Phy, "abundance", decreasing=TRUE)
tax_Phy

arrange(tax_Phy, tax_Phy$abundance)

#you can use now the custom function I have made 
source("../abund_Phy.R")
abund_Phy(biom_ITS, 54)
abund_Phy(biom_16s, 54)

biom_ITS_Phy = tax_glom(biom_ITS_cv, "Phylum")
write.table(otu_table(biom_ITS_cv), file = "ITS_phy.txt", sep="\t", eol="\n")



source("../rel_abund.R")
rel_abund(biom_ITS, 54, "Phylum")



biom_ITS_Ord = tax_glom(biom_ITS_cv, "Order")
tax_table(biom_ITS_Ord)


biom_ITS_Ord_ab = transform_sample_counts(biom_ITS_Ord, function(x) x/sum(x))
taxa_sums(biom_ITS_Ord_ab)/54
sum(taxa_sums(biom_ITS_Ord_ab)/54)
sample_sums(biom_ITS_Ord_ab)
tax_table(biom_ITS_Ord_ab)
tax_Ord <- data.frame(tax_table(biom_ITS_Ord))
tax_Ord <- tax_Ord[c(4)]
tax_Ord
tax_Ord$abundance <- as.vector(taxa_sums(biom_ITS_Ord_ab)/54)
#tax_Ord <- arrange(tax_Ord, abundance, decreasing=TRUE)
tax_Ord <- tax_Ord[order(tax_Ord$abundance, decreasing = TRUE),] 


biom_ITS_Gen = tax_glom(biom_ITS_cv, "Genus")
biom_ITS_Gen_ab = transform_sample_counts(biom_ITS_Gen, function(x) x/sum(x))
tax_table(biom_ITS_Gen_ab)
taxa_sums(biom_ITS_Gen_ab)/54
sum(taxa_sums(biom_ITS_Gen_ab)/54) 
sample_sums(biom_ITS_Gen_ab)
tax_Gen <- data.frame(tax_table(biom_ITS_Gen))
tax_Gen <- tax_Gen[c(6)]
tax_Gen$abundance <- as.vector(taxa_sums(biom_ITS_Gen_ab)/54)
#tax_Gen <- arrange(tax_Gen, abundance, decreasing=TRUE)
tax_Gen <- tax_Gen[order(tax_Gen$abundance, decreasing = TRUE),] 
tax_Gen


# OTUs abundances and data counts fungal LSU

biom_LSU_ab <- transform_sample_counts(biom_LSU_cv, function(x) x/sum(x))
taxa_sums(biom_LSU_ab)/51 # abundance of each of the first 100 top OTUs
sum(taxa_sums(biom_LSU_ab)/51) # sum of total abundance must be 1
sample_sums(biom_LSU_ab) # total abundace of each sample must be 1

otus_LSU <- data.frame(tax_table(biom_LSU_ab))
otus_LSU <- otus_LSU[c(1:6)]
otus_LSU$abundance <- as.vector(taxa_sums(biom_LSU_ab )/51)
otus_LSU <- otus_LSU[order(otus_LSU$abundance, decreasing = TRUE),] 
otus_LSU[1:50,]


biom_LSU_Phy = tax_glom(biom_LSU_cv, "Phylum")
biom_LSU_Phy_ab = transform_sample_counts(biom_LSU_Phy, function(x) x/sum(x))
sum(taxa_sums(biom_LSU_Phy_ab)/51)
sample_sums(biom_LSU_Phy_ab)
tax_Phy <- data.frame(tax_table(biom_LSU_Phy))
tax_Phy <- tax_Phy[c(2)]
tax_Phy$abundance <- as.vector(taxa_sums(biom_LSU_Phy_ab)/51)
tax_Phy <- tax_Phy[order(tax_Phy$abundance, decreasing = TRUE),] 
tax_Phy


biom_LSU_Ord = tax_glom(biom_LSU_cv, "Order")
biom_LSU_Ord_ab = transform_sample_counts(biom_LSU_Ord, function(x) x/sum(x))
taxa_sums(biom_LSU_Ord_ab)/51
sum(taxa_sums(biom_LSU_Ord_ab)/51)
sample_sums(biom_LSU_Ord_ab)
tax_table(biom_LSU_Ord_ab)
tax_Ord <- data.frame(tax_table(biom_LSU_Ord))
tax_Ord <- tax_Ord[c(4)]
tax_Ord
tax_Ord$abundance <- as.vector(taxa_sums(biom_LSU_Ord_ab)/51)
tax_Ord <- tax_Ord[order(tax_Ord$abundance, decreasing = TRUE),] 
tax_Ord 

biom_LSU_Gen = tax_glom(biom_LSU_cv, "Genus")
biom_LSU_Gen_ab = transform_sample_counts(biom_LSU_Gen, function(x) x/sum(x))
tax_table(biom_LSU_Gen_ab)
taxa_sums(biom_LSU_Gen_ab)/51
sum(taxa_sums(biom_LSU_Gen_ab)/51) 
sample_sums(biom_LSU_Gen_ab)
tax_Gen <- data.frame(tax_table(biom_LSU_Gen))
tax_Gen <- tax_Gen[c(6)]
tax_Gen$abundance <- as.vector(taxa_sums(biom_LSU_Gen_ab)/51)
tax_Gen <- tax_Gen[order(tax_Gen$abundance, decreasing = TRUE),] 
tax_Gen





## OTUs abundances and data counts bacteria
biom_16s_ab <- transform_sample_counts(biom_16s_cv, function(x) x/sum(x))
taxa_sums(biom_16s_ab)/54 # abundance of each of the first OTUs
sum(taxa_sums(biom_16s_ab)/54) # sum of total abundance must be 1
sample_sums(biom_16s_ab) # total abundace of each sample must be 1

otus_bact <- data.frame(tax_table(biom_16s_ab))
otus_bact <- otus_bact[c(1:7)]
otus_bact$abundance <- as.vector(taxa_sums(biom_16s_ab)/54)
otus_bact <- otus_bact[order(otus_bact$abundance, decreasing = TRUE),] 
otus_bact[1:50,]

write.csv(taxa_sums(fungi_top100)/54 , file = "fungi_test.csv")
write.csv(tax_table(fungi_top100), file = "fungi_test.csv")

biom_16s_phyla = tax_glom(biom_16s_cv, "Phylum")
biom_16s_phyla_ab = transform_sample_counts(biom_16s_phyla, function(x) x/sum(x))
taxa_sums(biom_16s_phyla_ab)/54
sum(taxa_sums(biom_16s_phyla_ab)/54)
sample_sums(biom_16s_phyla_ab)
tax_Phy <- data.frame(tax_table(biom_16s_phyla_ab))
tax_Phy <- tax_Phy[c(2)]
tax_Phy$abundance <- as.vector(taxa_sums(biom_16s_phyla_ab)/54)
tax_Phy <- tax_Phy[order(tax_Phy$abundance, decreasing = TRUE),] 
tax_Phy

biom_16s_Ord = tax_glom(biom_16s_cv, "Order")
biom_16s_Ord_ab = transform_sample_counts(biom_16s_Ord, function(x) x/sum(x))
taxa_sums(biom_16s_Ord_ab)/54
sum(taxa_sums(biom_16s_Ord_ab)/54)
sample_sums(biom_16s_Ord_ab)
tax_Ord <- data.frame(tax_table(biom_16s_Ord_ab))
tax_Ord <- tax_Ord[c(4)]
tax_Ord$abundance <- as.vector(taxa_sums(biom_16s_Ord_ab)/51)
tax_Ord <- tax_Ord[order(tax_Ord$abundance, decreasing = TRUE),] 
tax_Ord 


biom_16s_gen = tax_glom(biom_16s_cv, "Genus")
biom_16s_gen_ab = transform_sample_counts(biom_16s_gen, function(x) x/sum(x))
tax_table(biom_16s_gen_ab)
taxa_sums(biom_16s_gen_ab)/54
sum(taxa_sums(biom_16s_gen_ab)/54) 
sample_sums(biom_16s_gen_ab)
tax_Gen <- data.frame(tax_table(biom_16s_gen_ab))
tax_Gen <- tax_Gen[c(6)]
tax_Gen$abundance <- as.vector(taxa_sums(biom_16s_gen_ab)/51)
tax_Gen <- tax_Gen[order(tax_Gen$abundance, decreasing = TRUE),] 
tax_Gen



### DATA EXPLORATION USING BARPLOTS ----------------------------------------------------------------------------------------------

# basic exploratory barplots fungi
plot_bar(biom_ITS_ab, x="Description", fill="Class")

top50ITS <- names(sort(taxa_sums(biom_ITS_ab), TRUE)[1:50]) # first 50 top OTUs 
biom_ITS_top50 <- prune_taxa(top50ITS, biom_ITS_ab)
plot_bar(biom_ITS_top50, x="Site2", fill="Class")

# plotting merged samples from mapping categories
biom_ITS_des = merge_samples(biom_ITS, "Description") # merging samples
sample_data(biom_ITS_des)$Description <- levels(sample_data(biom_ITS)$Description) # relabelling 
biom_ITS_des = transform_sample_counts(biom_ITS_des, function(x) 100 * x/sum(x)) # calculating abundances
plot_bar(biom_ITS_des, x="Description", fill="Class")




# PLOT BARS using a custom function -------------------------------------------------------------------------------------------------------------------------------------------------
#https://github.com/joey711/phyloseq/issues/442

# dump("add2", file="plot_ordered_bar.R")
# plot_ordered_bar(physeq, x = "Sample", y = "Abundance", fill = NULL, leg_size = 0.5, title = NULL, facet_grid = NULL)
source("../plot_ordered_bar.R")

head(sample_data(biom_ITS_cv))
tax_table(biom_ITS_cv)

write.csv(tax_table(biom_ITS_cv) , file = "biom_ITS_cv.csv") #correct the tax_table
taxonomy_ITS_cv <- read.csv("biom_ITS_cv_corrected.csv", header=T, row.names =1)
tax_table(biom_ITS_cv) <- tax_table(as.matrix(taxonomy_ITS_cv))


biom_ITS_tax = tax_glom(biom_ITS_cv, "Order")
tax_table(biom_ITS_tax)

biom_ITS_tax_mer = merge_samples(biom_ITS_tax, "Site1")
sample_data(biom_ITS_tax_mer)$Site1 <- levels(sample_data(biom_ITS_tax)$Site1)

biom_ITS_tax_mer_ab = transform_sample_counts(biom_ITS_tax_mer, function(x) 100 * x/sum(x))
biom_ITS_tax_mer_ab
tax_table(biom_ITS_tax_mer_ab)
ITS_tax_mer_top20 <- names(sort(taxa_sums(biom_ITS_tax_mer_ab), TRUE)[1:20]) 
biom_ITS_tax_mer_ab_top20 <- prune_taxa(ITS_tax_mer_top20, biom_ITS_tax_mer_ab)
biom_ITS_tax_mer_ab_top20

write.csv(tax_table(biom_ITS_tax_mer_ab_top20) , file = "biom_ITS_tax_mer_ab_top20.csv") #correct the tax_table
taxonomy_ITS_top20 <- read.csv("biom_ITS_tax_mer_ab_top20_corrected.csv", header=T, row.names =1)
tax_table(biom_ITS_tax_mer_ab_top20) <- tax_table(as.matrix(taxonomy_ITS_top20))

#fungi_top100 <- subset_taxa(fungi_top100, Genus!="NA")
#fungi_top100 <- subset_taxa(fungi_top100, Genus!="g__unidentified")

# Changing order of levels
#sample_data(biom_ITS)$Level<-factor(sample_data(biom_ITS)$Level, levels=c("0","10","50","100"))
#levels(sample_data(biom_ITS)$Level)
#sample_data(biom_ITS)
# to extract a vector from the previously created dataframe
#x_axis <- sums_biom_ITS[["Description"]]
#str(x_axis)
#x_axis
#sample_data(biom_ITS)$Descritption<-factor(sample_data(biom_ITS)$Description, levels=x_axis)


facets <- c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Glomeromycota", "Rozellomycota","ZUnlassified","Zygomycota")
p_ITS <- plot_ordered_bar(biom_ITS_tax_mer_ab[tax_table(biom_ITS_tax_mer_ab)$Phylum %in% facets,],
                          x = "Site1", fill="Order", leg_size = 0.4, title="ITS")

p_ITS <- plot_ordered_bar(biom_ITS_tax_mer_ab, x = "Site1", fill="Order", leg_size = 0.4, title="ITS")


p_ITS + ylim(0, 100) +
scale_fill_manual(values=palette3) + 
guides(fill=guide_legend(reverse = TRUE, ncol=2, keywidth = 0.7, keyheight = 0.7)) +
labs(title="ITS", x="Site", y="Relative abundance") +
theme(panel.background = element_blank()) + 
theme(axis.text.x = element_text(vjust=0.5, size=10)) +
theme(axis.text.y = element_text(hjust=0.5, size=10)) +
theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
#theme(axis.line = element_line(colour = "black",size = 1, linetype = "solid")) +
theme(legend.position="right")

str(p_ITS)

ITS_bar = p_ITS + scale_fill_manual(values = c("Agaricales"="#781156",
                                             "Archaeorhizomycetales"="#A51876",
                                             "Atractiellales"="#D21E96",
                                             "Auriculariales"="#E43FAD",
                                             "Boletales"="#EA6CC0",
                                             "Botryosphaeriales"="#F098D3",
                                             "Cantharellales"="#114578",
                                             "Capnodiales"="#185EA5",
                                             "Chaetosphaeriales"="#1E78D2",
                                             "Chaetothyriales"="#3F91E4",
                                             "Chytridiales"="#6CABEA",
                                             "Coniochaetales"="#98C4F0",
                                             "Corticiales"="#117878",
                                             "Diaporthales"="#18A5A5",
                                             "Diversisporales"="#3FE4E4",
                                             "Eurotiales"="#6CEAEA",
                                             "Geminibasidiales"="#98F0F0",
                                             "Glomerales"="#117845",
                                             "Glomerellales"="#18A55E",
                                             "Helotiales"="#1ED278",
                                             "Hymenochaetales"="#3FE491",
                                             "Hypocreales"="#6CEAAB",
                                             "Jahnulales"="#98F0C4",
                                             "Kickxellales"="#787811",
                                             "Leotiales"="#A5A518",
                                             "Magnaporthales"="#D2D21E",
                                             "Microascales"="#E4E43F",
                                             "Mortierellales"="#EAEA6C",
                                             "Onygenales"="#F0F098",
                                             "Ophiostomatales"="#F7F7C5",
                                             "Orbiliales"="#784511",
                                             "Pezizales"="#A55E18",
                                             "Phallales"="#D2781E",
                                             "Pleosporales"="#E4913F",
                                             "Polyporales"="#EAAB6C",
                                             "Rhizophlyctidales"="#F0C498",
                                             "Rhizophydiales"="#781122",
                                             "Russulales"="#A5182F",
                                             "Saccharomycetales"="#D21E2C",
                                             "Sordariales"="#E43F5B",
                                             "Spizellomycetales"="#EA6C81",
                                             "Sporidiobolales"="#F098A7",
                                             "Thelephorales"="#242424",
                                             "Trechisporales"="#484848",
                                             "Tremellales"="#6D6D6D",
                                             "Trichosphaeriales"="#919191",
                                             "Trichosporonales"="#B6B6B6",
                                             "Xylariales"="#DADADA",
                                             "ZUnclassified"="#FFFFFF")) + 
  #theme(panel.background = element_blank()) + 
  guides(fill=guide_legend(reverse = TRUE, ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  #facet_wrap(~Phylum, ncol=4) +
  #facet_grid(~Phylum, nrow=2) +
  ylim(0, 100) +
  labs(title="ITS", x="", y="Relative abundance") +
  theme(panel.background = element_blank()) + 
  theme(axis.text.x = element_text(vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(hjust=0.5, size=10)) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  theme(axis.line = element_line(colour = "black",size = 1, linetype = "solid")) +
  theme(legend.position="right")



#scale_x_discrete(expand = c(0.2, 0)) +
#scale_x_discrete(limits = c("Beach","Cassava","Casuarina","Grassland","Mahogany",
 #                           "Forest","Pine","Stream bank","Sugarcane"), expand=c(0.2, 0)) +
 

# Change the line type and color of axis lines
#p + theme( axis.line = element_line(colour = "black",size = 1, linetype = "solid"))
  

  
 # ylim(0, 700) +
  #theme_minimal() +
  #geom_text_repel(aes(label= site), size=3, force = 25, col="black") +
  #facet_wrap(~site+groups) +
  #stat_smooth(method = "lm") +
  #geom_smooth(method = "lm") +
  #scale_colour_manual(values=c("#FF0000","#D03300","#A16600","#729900","#43CC00","#15FF00")) + 
  #scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000","#000000")) + 
  #scale_colour_manual(values=palette_CB6) + 
  #theme(axis.text.x = element_text(angle = 90, hjust=0.5, size=8)) +
  #theme(axis.text.y = element_text(hjust=0.5, size=8)) +
  #theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
#  labs(title="", x="Year", y="Observed richness") +
 # theme(axis.title = element_text(size = 8)) +  
  #expand_limits(x = 0, y = 0) +
  #theme(legend.position="none")


#p1 + scale_fill_manual(values=custom_col21) + ylim(0, 100) #+ coord_flip() #+ theme_bw()
#p1 + scale_fill_manual(values=custom_col21) + ylim(0, 50) + facet_wrap(~Genus) #+ coord_flip() #+ theme_bw()
#p1 + facet_wrap(~Genus) + ggtitle("Faceting is better than stacking...") + ylim(0, 50)
#p1 + geom_bar(aes(color=Genus,fill=Genus), stat="identity", position="stack") + ylim(0, 50)
#p1 + scale_fill_manual(values=custom_col21) + coord_flip()
#p1 + scale_colour_brewer(palette = "Set1") + coord_flip()
#p1 + scale_fill_manual(values=custom_col21) + geom_bar(aes(color=Genus,fill=Genus), stat="identity", position="stack")
#p1 + scale_fill_manual(values=custom_col21) + geom_bar(stat = "identity")


head(sample_data(biom_LSU_cv))
tax_table(biom_LSU_cv)

write.csv(tax_table(biom_LSU_cv) , file = "biom_LSU_cv.csv") #correct the tax_table
taxonomy_LSU_cv <- read.csv("biom_LSU_cv_corrected.csv", header=T, row.names =1)
tax_table(biom_LSU_cv) <- tax_table(as.matrix(taxonomy_LSU_cv))

biom_LSU_tax = tax_glom(biom_LSU_cv, "Order")
biom_LSU_tax_mer = merge_samples(biom_LSU_tax, "Site1")
sample_data(biom_LSU_tax_mer)$Site1 <- levels(sample_data(biom_LSU_tax)$Site1)

biom_LSU_tax_mer_ab = transform_sample_counts(biom_LSU_tax_mer, function(x) 100 * x/sum(x))
biom_LSU_tax_mer_ab
tax_table(biom_LSU_tax_mer_ab)

write.csv(tax_table(biom_LSU_tax_mer_ab) , file = "biom_LSU_tax_mer_ab.csv") #correct the tax_table
taxonomy_LSU_tax_mer_ab <- read.csv("biom_LSU_tax_mer_ab_corrected.csv", header=T, row.names =1)
tax_table(biom_LSU_tax_mer_ab) <- tax_table(as.matrix(taxonomy_LSU_tax_mer_ab))


#ITS_tax_mer_top20 <- names(sort(taxa_sums(biom_LSU_tax_mer_ab), TRUE)[1:20]) 
#biom_LSU_tax_mer_ab_top20 <- prune_taxa(ITS_tax_mer_top20, biom_LSU_tax_mer_ab)
#biom_LSU_tax_mer_ab_top20

p_LSU <- plot_ordered_bar(biom_LSU_tax_mer_ab, x = "Site1", fill="Order", leg_size = 0.4, title="LSU")

p_LSU + scale_fill_manual(values=palette3) + ylim(0, 100) #+ coord_flip() #+ theme_bw()

p_LSU  + ylim(0, 100) +
  scale_fill_manual(values=palette3) + 
  guides(fill=guide_legend(reverse = TRUE, ncol=2, keywidth = 0.7, keyheight = 0.7)) +
  labs(title="LSU", x="Site", y="Relative abundance") +
  theme(panel.background = element_blank()) + 
  theme(axis.text.x = element_text(vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(hjust=0.5, size=10)) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  #theme(axis.line = element_line(colour = "black",size = 1, linetype = "solid")) +
  theme(legend.position="right")



LSU_bar = p_LSU + scale_fill_manual(values = c("Agaricales"="#781156",
                                     "Archaeorhizomycetales"="#A51876",
                                     "Atractiellales"="#D21E96",
                                     "Auriculariales"="#E43FAD",
                                     "Blastocladiales"="#117878",
                                     "Boletales"="#EA6CC0",
                                     "Boliniales"="#98F0F0",
                                     "Botryosphaeriales"="#F098D3",
                                     "Cantharellales"="#114578",
                                     "Capnodiales"="#185EA5",
                                     "Chaetosphaeriales"="#1E78D2",
                                     "Chaetothyriales"="#3F91E4",
                                     "Chytridiales"="#6CABEA",
                                     "Cladochytriales"="#18A55E",
                                     "Classiculales"="#787811",
                                     "Coniochaetales"="#98C4F0",
                                     "Coronophorales"="#A5A518",
                                     "Cystobasidiomycetes"="#EAEA6C",
                                     "Diaporthales"="#18A5A5",
                                     "Diversisporales"="#3FE4E4",
                                     "Dothideomycetes"="#F0F098",
                                     "Erythrobasidiales"="#784511",
                                     "Eurotiales"="#6CEAEA",
                                     "Geastrales"="#F0C498",
                                     "Glomerales"="#117845",
                                     "Gomphales"="#A5182F",
                                     "Helicobasidiales"="#F098A7",
                                     "Helotiales"="#1ED278",
                                     "Hymenochaetales"="#3FE491",
                                     "Hypocreales"="#6CEAAB",
                                     "Hysteriales"="#484848",
                                     "Jahnulales"="#98F0C4",
                                     "Leotiomycetes"="#919191",
                                     "Lobulomycetales"="#B6B6B6",
                                     "Magnaporthales"="#D2D21E",
                                     "Melanosporales"="#bbfa0f",
                                     "Microascales"="#E4E43F",
                                     "Ophiostomatales"="#F7F7C5",
                                     "Paraglomerales"="#f95a04",
                                     "Pezizales"="#A55E18",
                                     "Phallales"="#D2781E",
                                     "Phyllachorales"="#a99249",
                                     "Pleosporales"="#E4913F",
                                     "Polyporales"="#EAAB6C",
                                     "Rhizophydiales"="#781122",
                                     "Saccharomycetales"="#D21E2C",
                                     "Sebacinales"="#E18E10",
                                     "Septobasidiales"="#D8A916",
                                     "Sordariales"="#E43F5B",
                                     "Spizellomycetales"="#EA6C81",
                                     "Thelephorales"="#242424",
                                     "Tremellales"="#6D6D6D",
                                     "Xylariales"="#DADADA",
                                     "ZUnclassified"="#FFFFFF")) + 
  #theme(panel.background = element_blank()) + 
  guides(fill=guide_legend(reverse = TRUE, ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  #facet_wrap(~Phylum, ncol=4) +
  #facet_grid(~Phylum, nrow=2) +
  ylim(0, 100) +
  labs(title="LSU", x="Site", y="") +
  theme(panel.background = element_blank()) + 
  theme(axis.text.x = element_text(vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(hjust=0.5, size=10)) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  theme(axis.line = element_line(colour = "black",size = 1, linetype = "solid")) +
  theme(legend.position="right")


#-#-#
head(sample_data(biom_16s_cv))
tax_table(biom_16s_cv)

write.csv(tax_table(biom_16s_cv) , file = "biom_16s_cv.csv") #correct the tax_table
taxonomy_16s_cv <- read.csv("biom_16s_cv_corrected.csv", header=T, row.names =1)
tax_table(biom_16s_cv) <- tax_table(as.matrix(taxonomy_16s_cv))


biom_16s_tax = tax_glom(biom_16s_cv, "Order")
tax_table(biom_16s_tax)

biom_16s_tax_mer = merge_samples(biom_16s_tax, "Site1")
sample_data(biom_16s_tax_mer)$Site1 <- levels(sample_data(biom_16s_tax)$Site1)

biom_16s_tax_mer_ab = transform_sample_counts(biom_16s_tax_mer, function(x) 100 * x/sum(x))
biom_16s_tax_mer_ab
tax_table(biom_16s_tax_mer_ab)
tax_16s_mer_top20 <- names(sort(taxa_sums(biom_16s_tax_mer_ab), TRUE)[1:50]) 
biom_16s_tax_mer_ab_top20 <- prune_taxa(tax_16s_mer_top20, biom_16s_tax_mer_ab)
biom_16s_tax_mer_ab_top20

write.csv(tax_table(biom_16s_tax_mer_ab_top20) , file = "biom_16s_tax_mer_ab_top20.csv") #correct the tax_table
taxonomy_16s_top20 <- read.csv("biom_16s_tax_mer_ab_top20_corrected.csv", header=T, row.names =1)
tax_table(biom_16s_tax_mer_ab_top20) <- tax_table(as.matrix(taxonomy_16s_top20))

#fungi_top100 <- subset_taxa(fungi_top100, Genus!="NA")
#fungi_top100 <- subset_taxa(fungi_top100, Genus!="g__unidentified")

# Changing order of levels
#sample_data(biom_ITS)$Level<-factor(sample_data(biom_ITS)$Level, levels=c("0","10","50","100"))
#levels(sample_data(biom_ITS)$Level)
#sample_data(biom_ITS)
# to extract a vector from the previously created dataframe
#x_axis <- sums_biom_ITS[["Description"]]
#str(x_axis)
#x_axis
#sample_data(biom_ITS)$Descritption<-factor(sample_data(biom_ITS)$Description, levels=x_axis)



p_16s <- plot_ordered_bar(biom_16s_tax_mer_ab_top20, x = "Site1", fill="Order", leg_size = 0.4, title="ITS")


p_16s + ylim(0, 100) +
  scale_fill_manual(values=palette3) + 
  guides(fill=guide_legend(reverse = TRUE, ncol=2, keywidth = 0.7, keyheight = 0.7)) +
  labs(title="ITS", x="Site", y="Relative abundance") +
  theme(panel.background = element_blank()) + 
  theme(axis.text.x = element_text(vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(hjust=0.5, size=10)) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  #theme(axis.line = element_line(colour = "black",size = 1, linetype = "solid")) +
  theme(legend.position="right")


bar_16s = p_16s + scale_fill_manual(values = c("0319-7L14"="#781156",
                                               "32-20"="#A51876",
                                               "Acidimicrobiales"="#D21E96",
                                               "Acidobacteriales"="#E43FAD",
                                               "Actinomycetales"="#EA6CC0",
                                               "AKYG885"="#F098D3",
                                               "B07_WMSP1"="#114578",
                                               "Bacillales"="#185EA5",
                                               "Burkholderiales"="#1E78D2",
                                               "Caulobacterales"="#3F91E4",
                                               "CCU21"="#6CABEA",
                                               "[Chthoniobacterales]"="#98C4F0",
                                               "Cytophagales"="#117878",
                                               "Ellin329"="#18A5A5",
                                               "Ellin6513"="#3FE4E4",
                                               "[Entotheonellales]"="#6CEAEA",
                                               "Gaiellales"="#98F0F0",
                                               "Gemmatales"="#117845",
                                               "iii1-15"="#18A55E",
                                               "Ktedonobacterales"="#1ED278",
                                               "MND1"="#3FE491",
                                               "Myxococcales"="#6CEAAB",
                                               "NB1-j"="#98F0C4",
                                               "Nitrososphaerales"="#787811",
                                               "Nitrospirales"="#A5A518",
                                               "[Pedosphaerales]"="#D2D21E",
                                               "Phycisphaerales"="#E4E43F",
                                               "Pirellulales"="#EAEA6C",
                                               "PK29"="#F0F098",
                                               "RB41"="#F7F7C5",
                                               "Rhizobiales"="#784511",
                                               "Rhodospirillales"="#A55E18",
                                               "Rubrobacterales"="#D2781E",
                                               "[Saprospirales]"="#E4913F",
                                               "SBR1031"="#EAAB6C",
                                               "Sediment-1"="#F0C498",
                                               "Solibacterales"="#781122",
                                               "Solirubrobacterales"="#A5182F",
                                               "Sphingomonadales"="#D21E2C",
                                               "Sva0725"="#E43F5B",
                                               "Syntrophobacterales"="#EA6C81",
                                               "Thiotrichales"="#F098A7",
                                               "WD2101"="#242424",
                                               "Xanthomonadales"="#6D6D6D",
                                               "ZUnclassified"="#FFFFFF")) + 
  #theme(panel.background = element_blank()) + 
  guides(fill=guide_legend(reverse = TRUE, ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  #facet_wrap(~Phylum, ncol=4) +
  #facet_grid(~Phylum, nrow=2) +
  ylim(0, 100) +
  labs(title="16S", x="", y="") +
  theme(panel.background = element_blank()) + 
  theme(axis.text.x = element_text(vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(hjust=0.5, size=10)) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  theme(axis.line = element_line(colour = "black",size = 1, linetype = "solid")) +
  theme(legend.position="right")


bar_16s



library("ggpubr")
ggarrange(ITS_bar,LSU_bar,bar_16s,
          labels = c("A", "B","C"),
          widths = c(1, 1,1),
          align = "h", ncol = 3, nrow = 1,
          legend="right",
          common.legend = FALSE)









#-#-#
head(sample_data(biom_16s_cv))

write.csv(tax_table(biom_16s_cv) , file = "biom_16s_cv.csv") #correct the tax_table
taxonomy_16s_cv <- read.csv("biom_16s_cv_corrected.csv", header=T, row.names =1)
tax_table(biom_16s_cv) <- tax_table(as.matrix(taxonomy_16s_cv))


biom_16s_tax = tax_glom(biom_16s_cv, "Class")
biom_16s_tax_mer = merge_samples(biom_16s_tax, "Site1")
sample_data(biom_16s_tax_mer)$Site1 <- levels(sample_data(biom_16s_tax)$Site1)

biom_16s_tax_mer_ab = transform_sample_counts(biom_16s_tax_mer, function(x) 100 * x/sum(x))
biom_16s_tax_mer_ab
ITS_tax_mer_top20 <- names(sort(taxa_sums(biom_16s_tax_mer_ab), TRUE)[1:50]) 
biom_16s_tax_mer_ab_top20 <- prune_taxa(ITS_tax_mer_top20, biom_16s_tax_mer_ab)
biom_16s_tax_mer_ab_top20

write.csv(tax_table(biom_16s_tax_mer_ab_top20) , file = "biom_16s_tax_mer_ab_top20.csv") #correct the tax_table
taxonomy_16s_top20 <- read.csv("biom_16s_tax_mer_ab_top20_corrected.csv", header=T, row.names =1)
tax_table(biom_16s_tax_mer_ab_top20) <- tax_table(as.matrix(taxonomy_16s_top20))

p_16s <- plot_ordered_bar(biom_16s_tax_mer_ab, x = "Site1", fill="Class", leg_size = 0.4, title="16S")

p_16s + scale_fill_manual(values=palette3) + ylim(0, 100) #+ coord_flip() #+ theme_bw()

p_16s  + ylim(0, 100) +
  scale_fill_manual(values=palette3) + 
  guides(fill=guide_legend(reverse = TRUE, ncol=3, keywidth = 0.7, keyheight = 0.7)) +
  labs(title="16S", x="Site", y="Relative abundance") +
  theme(panel.background = element_blank()) + 
  theme(axis.text.x = element_text(vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(hjust=0.5, size=10)) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  #theme(axis.line = element_line(colour = "black",size = 1, linetype = "solid")) +
  theme(legend.position="right")










# check for layers
p1$layers

p2 <- plot_ordered_bar(biom_ITS_gen_mer20, x="Soil", y="Abundance", fill="Family", leg_size=0.5)
p2 + coord_flip() + scale_fill_manual(values=custom_col21)

p3 <- plot_bar(biom_ITS_ev, "Genus", fill="Genus", facet_grid=Level~Soil)
p3 + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") 

p4 <-plot_bar(biom_ITS_top50, x="Soil", fill="Genus", facet_grid=~Soil)
p4 + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")


# trying to plot using ggplot() instead of plot_bar()
physeqdf$Genus <- factor(physeqdf$Genus, levels = c(Isolate)) ## to reorder

physeqdf <- psmelt(fungi_top100)
physeqdf
str(physeqdf)
p <- ggplot(physeqdf, aes(x=Sample, y=Abundance, fill=Genus, order = as.factor(Genus)))
p + geom_bar(stat="identity") + coord_flip() + scale_fill_manual(values=palette20) + theme_bw()

p <- plot_taxa_bar(ent10, "Genus", x = "SeqTech", fill = "TaxaGroup") + facet_wrap(~Enterotype)


##exploratory barplots bacteria
biom_16s_class = tax_glom(biom_16s, "Class")
biom_16s_class_mer= merge_samples(biom_16s_class, "Description")
sample_data(biom_16s_class_mer)$Description <- levels(sample_data(biom_16s_class)$Description)
biom_16s_class_mer = transform_sample_counts(biom_16s_class_mer, function(x) 100 * x/sum(x))
top20_16s <- names(sort(taxa_sums(biom_16s_class_mer), TRUE)[1:20])
biom_16s_class_top20<- prune_taxa(top20_16s, biom_16s_class_mer)
sum(sample_sums(biom_16s_class_top20)/34)


p1 <- plot_ordered_bar(biom_16s_class_top20, x = "Description", fill="Class", leg_size = 0.4, title="Bacteria Barplots") 
p1 + coord_flip() + scale_fill_manual(values=custom_col21) + ylim(0, 100)
p1 + scale_fill_manual(values=palette20) + ylim(0, 100) #+ theme_bw() + coord_flip() 

## there are plotted single otu abundances and grouped as fill= in the legend??
level_16s = merge_samples(biom_16s, "Description")
sample_data(level_16s)$Description <- levels(sample_data(biom_16s)$Description)
level_16s = transform_sample_counts(level_16s, function(x) 100 * x/sum(x))
level_16s
top10016s <- names(sort(taxa_sums(biom_16s), TRUE)[1:100]) 
top10016s
bact_top100 <- prune_taxa(top10016s, level_16s)
bact_top100
otu_table(bact_top100)
tax_table(bact_top100)
#bact_top100 <- subset_taxa(bact_top100, Genus!="o__")
#bact_top100 <- subset_taxa(bact_top100, Genus!="g__unidentified")
bact_top100

p1 <- plot_ordered_bar(bact_top100, x = "Description", fill="Class", leg_size = 0.4, title="Bacteria Barplots") 
p1 + coord_flip() + scale_fill_manual(values=palette5)
p1 + coord_flip() + scale_fill_manual(values=palette) + geom_bar(aes(color=Class,fill=Class), stat="identity", position="stack")

p1 + facet_wrap(~Class) + ggtitle("Faceting is better than stacking...") + ylim(0, 50)



### faceting most abundat otus
TopOTUs <- names(sort(taxa_sums(biom_ITS_ab), TRUE)[1:10])
fungi10 <- prune_taxa(TopOTUs, biom_ITS_ab)
p1 <- plot_ordered_bar(fungi10, x = "Description", fill="Species", leg_size = 0.4, title="Fungi Barplots") 
p1 + scale_fill_manual(values=palette20) + ylim(0, 50) + facet_wrap(~Genus) #+ coord_flip() #+ theme_bw()

plot_bar(fungi10, fill="Species", facet_grid = ~Species) + scale_fill_manual(values=palette20)

TopOTUs <- names(sort(taxa_sums(biom_ITS_ab), TRUE)[1:10])
fungi10 <- prune_taxa(TopOTUs, biom_ITS_ab)
plot_bar(fungi5, fill="Genus", facet_grid = ~Species) + scale_fill_manual(values=rainbow(10))


class_barITS<-plot_bar(fungi10, "Genus", fill = "Genus", facet_grid=~genotype2~Level2)
class_barITS<-class_barITS + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
class_barITS


## faceting instead of staking barplots

biom_ITS_fam = tax_glom(biom_ITS_cv, "Family")
biom_ITS_fam_mer = merge_samples(biom_ITS_fam, "Site2")
sample_data(biom_ITS_fam_mer)$Site2 <- levels(sample_data(biom_ITS_fam)$Site2)
#sample_data(biom_ITS_fam_mer)$Treatment <- levels(sample_data(biom_ITS_fam)$Treatment)
biom_ITS_fam_mer_ab = transform_sample_counts(biom_ITS_fam_mer, function(x) 100 * x/sum(x))
top20ITS <- names(sort(taxa_sums(biom_ITS_fam_mer_ab), TRUE)[1:10]) 
biom_ITS_fam_mer20 <- prune_taxa(top20ITS, biom_ITS_fam_mer_ab)



plot_bar(biom_ITS_fam_mer20, "Site2", fill="Family", facet_grid=~Family)

p = plot_bar(biom_ITS_gen_mer20,"Site1", fill="Phylum", facet_grid=~Site1)
p + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

http://joey711.github.io/phyloseq-demo/Restroom-Biogeography


******************************************************************************************************************************************


### filter and then plot bars for each level of a variable
require(grDevices)
library(graphics)
library(plyr)

biom_ITS_Site2 = merge_samples(biom_ITS_cv, "Site2") 

biom_ITS_Cassava1 <- subset_samples(biom_ITS_Site2, Site2%in%c("3"))
biom_ITS_Cassava1
otu_table(biom_ITS_Cassava1)[apply(otu_table(biom_ITS_Cassava1), MARGIN = 1, function(x) any(x >= 5)), ]
otu_table(biom_ITS_Cassava1) <- otu_table(biom_ITS_Cassava1)[which(rowSums(otu_table(biom_ITS_Cassava1)) > 10),] 
biom_ITS_Cassava1
tax_table(biom_ITS_Cassava1)
biom_ITS_Cassava1_prop <- transform_sample_counts(biom_ITS_Cassava1, function(x) x/sum(x))
top10_ITS_Cassava1 <- names(sort(taxa_sums(biom_ITS_Cassava1_prop), TRUE)[1:10]) # sampling first 10 taxa
biom_ITS_Cassava1 <- prune_taxa(top10_ITS_Cassava1, biom_ITS_Cassava1_prop)
tax_table(biom_ITS_Cassava1)

fungi_ITS_Cassava1 <- data.frame(t(as.matrix(otu_table(biom_ITS_Cassava1))))
fungi_ITS_Cassava1
fungi_ITS_Cassava1_tax <- data.frame(as.matrix(tax_table(biom_ITS_Cassava1)))
fungi_ITS_Cassava1_tax

fungi_ITS_Cassava1 <- fungi_ITS_Cassava1[order(fungi_ITS_Cassava1$Cassava1, decreasing = TRUE),] 

fungi_ITS_Cassava1$OTU<- row.names(fungi_ITS_Cassava1_tax)
fungi_ITS_Cassava1$Species <- sample_data(fungi_ITS_Cassava1_tax)$Species
fungi_ITS_Cassava1$Names = paste(fungi_ITS_Cassava1$OTU, fungi_ITS_Cassava1$Species, sep=" ")
fungi_ITS_Cassava1$Names = paste(fungi_ITS_Cassava1$Names, fungi_ITS_Cassava1$Cassava1, sep=" ")


#try this
https://rpubs.com/michberr/randomforestmicrobe



pie(fungi_ITS_Cassava1$Cassava1, labels=fungi_ITS_Cassava1$Names, 
    col=custom_col21, cex=0.8, main="Cassava1")

plot_bar(biom_ITS_Cassava1, x="Description", fill="Genus")


biom_ITS_Site2 = merge_samples(biom_ITS_cv, "Site2") 
biom_ITS_Cassava1_prop <- transform_sample_counts(biom_ITS_Cassava1, function(x) x/sum(x))


## for filtering look script
filetering_Site2.R


# work on this 
p = plot_bar(biom_ITS_Site2, "Genus", fill="Genus", facet_grid=~Site2)
p + geom_point(aes(x=Genus, y=Abundance), color="black", position="jitter", size=3)








#labels=y$Names, col=terrain.colors(nrow(y)) , main=c(y$country[1]))
#colnames(fungi_ITS_Cassava1) <- "Taxonomy"

sums_biom_ITS <- data.frame(colSums(otu_table(biom_ITS)))
colnames(sums_biom_ITS) <- "Sample_TotalSeqs"
sums_biom_ITS$sample <- row.names(sums_biom_ITS)
sums_biom_ITS$Description <- sample_data(biom_ITS)$Description
sums_biom_ITS$Site2 <- sample_data(biom_ITS)$Site2 # if you want to add another variable to use further
sums_biom_ITS <- arrange(sums_biom_ITS, Sample_TotalSeqs)
#sums_biom_ITS <- arrange(sums_biom_ITS, Description) # to order according another variable in the mapping
str(sums_biom_ITS)

pie <- ggplot(fungi_ITS_Cassava1, aes(x = factor(1), fill = factor(Cassava1))) + geom_bar(width = 1)
pie + coord_polar(theta = "y")


### pie charts
fungi_ITS_cv <- as.data.frame(t(as.matrix(otu_table(biom_ITS_cv))))
fungi_ITS_cv



pie + coord_polar(theta = "y") + facet_wrap(~factor(gear))






### INDICATOR SPECIES ANALYSIS ------------------------------------------------------------------------------------------------------------

#indicator species analysis (isa) Fungi
library(indicspecies)

biom_ITS_symb
biom_ITS_endo
biom_ITS_pat

otu_ITS_symb <- as.data.frame(t(otu_table(biom_ITS_symb)))
metadata_ITS_symb <- as.data.frame(as.matrix(sample_data(biom_ITS_symb)[,4:9]))
tax_table(biom_ITS_symb)
otu_ITS_symb
metadata_ITS_symb

isa_fungi_its<-multipatt(otu_ITS_symb, metadata_ITS_symb$Site1, control=how(nperm=9999), duleg=TRUE, func="r")
isa_fungi_its$sign$p.value<-p.adjust(isa_fungi_its$sign$p.value, "fdr")
summary(isa_fungi_its)
isa_fungi_its

otu_ITS_endo <- as.data.frame(t(otu_table(biom_ITS_endo)))
metadata_ITS_endo <- as.data.frame(as.matrix(sample_data(biom_ITS_endo)[,4:9]))
tax_table(biom_ITS_endo)
otu_ITS_endo
metadata_ITS_endo

isa_fungi_its<-multipatt(otu_ITS_endo, metadata_ITS_endo$Site1, control=how(nperm=9999), duleg=TRUE, func="r")
isa_fungi_its$sign$p.value<-p.adjust(isa_fungi_its$sign$p.value, "fdr")
summary(isa_fungi_its)



otu_ITS_pat <- as.data.frame(t(otu_table(biom_ITS_pat)))
metadata_ITS_pat <- as.data.frame(as.matrix(sample_data(biom_ITS_pat)[,4:9]))
tax_table(biom_ITS_pat)
otu_ITS_pat
metadata_ITS_pat

isa_fungi_its<-multipatt(otu_ITS_pat, metadata_ITS_pat$Site1, control=how(nperm=9999), duleg=TRUE, func="r")
isa_fungi_its$sign$p.value<-p.adjust(isa_fungi_its$sign$p.value, "fdr")
summary(isa_fungi_its)




# inspecting the indicator value components
#(1) Component A is the probability that the surveyed site  belongs  to
#the  target  site  group - specificity or positive
#predictive value. (2) Component B is the probability of  finding
#the species in sites belonging to the site group - fidelity or sensitivity.
summary(isa_fungi, indvalcomp=TRUE)
summary(isa_bact, indvalcomp=TRUE)

isa_bact_out <- capture.output(summary(isa_bact, indvalcomp=TRUE))
write.csv(as.data.frame(isa_bact_out), file = "isa_bact.csv")
write.table(as.data.frame(isa_bact_out), file = "isa_bact.txt", sep="\t", eol = "\n")
#cat("indicator species analysis bacteria", isa_bact_out, file="isa_bact.txt", sep="\t", append=TRUE)

# display the result of the indicator species analysis for all species
summary(isa_fungi, alpha=1)

# to list also species that occur in sites belonging to all groups
isa_fungi$sign

# Analyzing species ecological preferences with correlation indices
# using Pearson's phi coefficient of association. the func="r.g." is for
# unequal number of sites inside the groups
fungi_cv -> fungi_pa
fungi_pa[fungi_pa > 0] <- 1
head(fungi_pa)
fungi_phi = multipatt(as.data.frame(fungi_pa), metadata_ITS$Level, func = "r.g", control = how(nperm=999))
summary(fungi_phi)

# the phi and point biserial coefficients can take negative values expressing the fact that a 
# species tends to 'avoid' particular environmental conditions (i.e., group pf samples)
round(head(fungi_phi$str),8)
round(head(isa_fungi$str),8)

# other species to groups association functions see ?stressoc for details
fungi_prefstat = strassoc(fungi_cv, cluster=metadata_ITS$Level, func="A.g")
fungi_prefstat = strassoc(as.data.frame(fungi), cluster=metadata_ITS$Treatment, func="A.g", nboot=99)
fungi_prefstat
round(head(fungi_prefstat),8)

# sometimes useful to know the  proportion  of sites of a given
# site group where one or another indicator is found
coverage(as.data.frame(fungi_cv), isa_fungi)
plotcoverage(as.data.frame(fungi), isa_fungi, group="0", lty=1)
plotcoverage(as.data.frame(fungi), isa_fungi, group="TNPo3", lty=2, col="blue", add=TRUE)
plotcoverage(as.data.frame(fungi), isa_fungi, group="TNPo5", lty=3, col="red", add=TRUE)
legend(x = 0.01, y=20,legend=c("TNPo1","TNPo3", "TNPo5"), lty=c(1,2,3), col=c("black","blue","red"), bty="n")

# Species combinations as indicators of site groups as two or three species,
# when found together, bear more ecological information than a single one
fungi_comb = combinespecies(as.data.frame(fungi), max.order = 2)$XC
dim(fungi_comb)
isa_fungi_comb = multipatt(fungi_soil_comb, metadata_ITS$Treatment, duleg = TRUE, control = how(nperm=999))
isa_fungi_comb$sign$p.value<-p.adjust(isa_fungi_comb$sign$p.value, "fdr")
summary(isa_fungi_comb, indvalcomp = TRUE)


# trying http://ecology.msu.montana.edu/labdsv/R/labs/lab1/lab1.html
library(labdsv)

dis.bc <- dsvdis(t(fungi_vd),'bray/curtis') # returns a dissimilarity matrix
clust <- sample(1:4,nrow(t(fungi_vd)),replace=TRUE)
indval(t(fungi_vd),clust)


## SIMPER
sim <-simper(as.data.frame(t(r1_ITS_veg)), map$Level)
sim
summary(sim, perm=9999)




### EXTRACT DATA FROM ISA 
library(dplyr)

summary(isa_fungi)
head(isa_fungi$sign)
str(isa_fungi$sign)

result_isa_fungi <- subset(isa_fungi$sign, p.value <= 0.05, select = c(s.0, s.10, s.100, s.50,p.value))
result_isa_fungi <- subset(result_isa_fungi, rowSums(result_isa_fungi) <= 2, select = c(p.value))
result_isa_fungi


#extract significant elements
isa_fungi_sig <- subset(cbind(isa_fungi$str, isa_fungi$sign$p.value), isa_fungi$sign$p.value < 0.05)
isa_fungi_sig


#subset your table accordingly, here a phyloseq object, as an example
biom_multipatt <- biom
otu_table(biom_multipatt) <-otu_table(biom)[row.names(multipatt_sig), ]







##### perform mantel test on community data

http://rfunctions.blogspot.com/2016/10/simple-and-partial-mantel-tests.html
http://stats.idre.ucla.edu/r/faq/how-can-i-perform-a-mantel-test-in-r/
http://www.unh.edu/halelab/BIOL933/labs/lab6.pdf  # data transofrmation
  
  
library(ade4)

bacteria_cv
fungi_ITS_cv
fungi_LSU_cv

geo_fiji_16s<- read.table("geo_fiji_16s.csv", header=T, row.names=1, sep=",")
geo_fiji_ITS<- read.table("geo_fiji_ITS.csv", header=T, row.names=1, sep=",")
geo_fiji_LSU<- read.table("geo_fiji_LSU.csv", header=T, row.names=1, sep=",")
geo_fiji_ITS

commdist_ITS <- vegdist(t(fungi_ITS_cv), method="bray")
commdist_ITS
geodists_ITS <- vegdist(cbind(geo_fiji_ITS$Lon, fiji_geo$Lat), method="euclidean")
geodists_ITS

geo_fiji_ITS$Lat <- log10(geo_fiji_ITS$Lat)
geo_fiji_ITS
geo_fiji_ITS$Lon <- log10(geo_fiji_ITS$Lon)

mantel.rtest(commdist_ITS, geodists_ITS, nrepet = 9999)

plot(commdist_ITS, geodists_ITS,pch=16,cex=0.5,col="black",bty="l",xlab="Spatial distance",ylab="Species composition dissimilarity")
mantel(commdist_ITS, geodists_ITS, method="pearson", permutations=999)

mantel_correg_ITS = mantel.correlog(commdist_ITS, geodists_ITS, nperm=999)
summary(mantel_correg_ITS)
mantel_correg_ITS  
plot(mantel_correg_ITS)






# BUBBLE CHARTS IN R
http://flowingdata.com/2010/11/23/how-to-make-bubble-charts/

symbols(result_isa_fungi$p.value, result_isa_fungi$p.value, circles=result_isa_fungi$p.value)

library(reshape2)
library(ggplot2)

fungi_melt <- melt(fungi_cv, id.vars = "Species", variable.name="Samples", value.name = "Size")
str(fungi_melt)
head(fungi_melt)


ggplot(fungi_melt, aes(x = Var2, y = Size)) +
  geom_point(aes(size = Size)) + 
  scale_size(range = range(fungi_melt$Size)) +
  theme_bw()






### DIFFERENTIAL ABUNDANCE ANALYSIS USING DESeq2
# http://joey711.github.io/phyloseq-extensions/DESeq2.html
# http://joey711.github.io/phyloseq-extensions/DESeq2.html

library(DESeq2)

# subsetting to specific sample groups
biom_ITS_nogenus <- subset_taxa(biom_ITS, Genus!="g__unidentified")
biom_ITS_nogenus

biom_ITS_des2 = phyloseq_to_deseq2(biom_ITS_nogenus, ~ Level2)
biom_ITS_des2
biom_ITS_des2 = DESeq(biom_ITS_des2, test="Wald", fitType="parametric")
res = results(biom_ITS_des2, cooksCutoff = FALSE)
res

res.corrected = p.adjust(res$pvalue, "fdr")
res.corrected

alpha = res.corrected
alpha = 0.01

sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(biom_ITS)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)
sigtab


biom_ITS_des2 = phyloseq_to_deseq2(biom_16s, ~ Level2)
biom_ITS_des2
biom_ITS_des2 = DESeq(biom_ITS_des2, test="Wald", fitType="parametric")
res = results(biom_ITS_des2, cooksCutoff = FALSE)
res
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(biom_16s)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)
sigtab

write.csv(sigtab, file = "sigtab.csv") #correct the DESeq2 result table
sigtab <- read.csv("sigtab.csv", header=T, row.names =1)


# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Species, function(x) max(x))
x = sort(x, TRUE)
sigtab$Species = factor(as.character(sigtab$Species), levels=names(x))

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

dev.off()

ggplot(sigtab, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + geom_point(size=4) + 
  theme(plot.title = element_text("Differential Abundances", size=16), 
        axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=12)) #+ scale_color_manual(values=palette20)


ggplot(sigtab, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + geom_boxplot() + #+ geom_point(size=4)
  theme(plot.title = element_text("Differential Abundances"), 
        axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=12)) #+ scale_color_manual(values=palette20)

# + theme(axis.text.x = element_text(colour="grey20",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),
#      axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
#     axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
#    axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"))


## summerize results in a table
posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ] #select otus that increase according treatment
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus", "Species")]
posigtab

sort(posigtab$log2FoldChange, TRUE)





### HEATMAPS 

biom_ITS_Ord = tax_glom(biom_ITS_cv, "Order")
biom_ITS_Ord_ab = transform_sample_counts(biom_ITS_Ord, function(x) x/sum(x))

biom_ITS_ab <- transform_sample_counts(biom_ITS_cv, function(x) x/sum(x))
heat_ITS <- prune_taxa(names(sort(taxa_sums(biom_ITS_ab), TRUE)[1:30]), biom_ITS_cv)

write.csv(tax_table(heat_ITS) , file = "taxtab_ITS_top30.csv") #correct the tax_table
taxtab_ITS_top30 <- read.csv("taxtab_ITS_top30.csv", header=T, row.names =1)
tax_table(heat_ITS) <- tax_table(as.matrix(taxtab_ITS_top30))

plot_heat_ITS <- plot_heatmap(heat_ITS, sample.label = "Site2", sample.order="Site2")
plot_heat_ITS + theme (axis.text.y = element_text(size=6))
plot_heat_ITS + theme(axis.text.x = element_text(size=9, angle=90, hjust=1), axis.text.y = element_text(size=9))

library("scales")
plot_heat_ITS <- plot_heatmap(heat_ITS, "PCoA", "bray", "Site2", "Order", sample.order="Site2", title="Heatmap ITS", low = "#000033",
                              high = "#66CCFF", trans=identity_trans())

plot_heat_ITS + theme(axis.text.y = element_text(size=6))


plot_heatmap(TOP20KEGG, method="PCoA", distance="bray", taxa.label="Pathway", title=title, trans = identity_trans())

plot_heat_ITS <- plot_heatmap(heat_ITS, "NMDS", "bray", "Site2", "Genus", sample.order="Site2")

## Heatmap examples - selecting only fungal OTUs

heat <- subset_taxa(fungi, Kingdom == "k__Fungi")
heat
plot_heatmap(heat, sample.label = "Site")
heat_plot1 <- plot_heatmap(heat, sample.label = "Site", sample.order="Site")
heat_plot1 + theme (axis.text.y = element_text(size=3))

heat_ab <- transform_sample_counts(heat, function(x) x/sum(x))
heat_ab

heat_ab <- prune_taxa(names(sort(taxa_sums(heat), TRUE)[1:50]), fungi)
plot_heatmap(heat_ab, sample.label = "Site", sample.order="Site")
heat_plot2 <- plot_heatmap(heat_ab, sample.label = "Site", sample.order="Site")
heat_plot2 + theme (axis.text.y = element_text(size=6))

heat_plot2 <- plot_heatmap(heat_ab, sample.label = "Site", sample.order="Site", na.value="white")
heat_plot2 + theme (axis.text.y = element_text(size=6))

heat_plot2 + theme(axis.text.x = element_text(size = 11, angle = 90, hjust = 1),
                   axis.text.y = element_text(size = 11))

## Re-label HEATMAP by a sample variable and taxonomic family
fungi_ab <- transform_sample_counts(fungi, function(x) x/sum(x))
p <- plot_heatmap(fungi, "NMDS", "bray", "Site", "Species", sample.order="Site")
p + theme (axis.text.y = element_text(size=4))

p2 <- plot_heatmap(heat_ab, "NMDS", "bray", "Site", "Species", sample.order="Site")
p2 + theme (axis.text.y = element_text(size=6))

p3 <- plot_heatmap(heat_ab, "NMDS", "bray", "Site", sample.order="Site")
p3 + theme (axis.text.y = element_text(size=6))





## PERFORM CORRELATION ANALYSIS AMONG OTUs

# corr all var in a matrix and plot
table<-read.table("table.txt", skip=0, comment.char = "", check.names = FALSE, header=TRUE, row.names=1, sep="\t") #import metadata
r_cor<-cor(table, use="pairwise.complete.obs") # deal with NA values
corrplot(r_cor, order = "hclust", method="square")

p_cor<-rcorr(as.matrix(table))$P
corrplot(r_cor, p.mat = p_cor, sig.level = 0.05, order="hclust") #corrplot with significance test

cex.before <- par("cex")
par(cex = 0.5)
png(height=1200, width=1500, pointsize=15, file="cor.png")
corrplot.mixed(r_cor, upper='number',lower='square', p.mat=p_cor, insig='blank', order="hclust", tl.cex = 1/par("cex"), cl.cex = 1/par("cex")) #mixed corrplot with r and p-values 
dev.off()


#corr var columns of a matrix with var columns in a second matrix 
biom_class<-t(read.table("class.txt", skip=1, comment.char = "", check.names = FALSE, header=TRUE, row.names=1, sep="\t")) #import QIIME class table
class_table_corr.test<-corr.test(as.matrix(biom_class), as.matrix(table), adjust="fdr")
corrplot(class_table_corr.test$r, p.mat = class_table_corr.test$p, sig.level = 0.05) #hybrid corrplot with significance test


#corr one var with another
corr_plot<-ggplot(data=table,aes(x=as.data.frame(biom_class)$`k__Bacteria;p__Bacteroidetes;c__Bacteroidia`, y=table$WHtR)) + theme_bw() + geom_point() + geom_smooth() + ylab("Bacteroidia") + xlab("WHtR")
cor.test(as.data.frame(biom_class)$`k__Bacteria;p__Bacteroidetes;c__Bacteroidia`, table$WHtR)


#corr 2 matrices via Mantel test
biom_table_ab<-biom_table/rowSums(biom_table)
biom_table_log<-decostand(biom_table_ab, method="log")
biom_table_bray<-vegdist(biom_table_log, method="bray")
biom_class_log<-decostand(biom_class, method="log") #QIIME class table was already in % abundance
biom_class_bray<-vegdist(biom_class_log, method="bray")
mantel(biom_table_bray, biom_class_bray, method="pearson", permutations=9999, parallel=4)

-------
  #links
  -------
  http://joey711.github.io/phyloseq/tutorials-index
https://web.stanford.edu/class/bios221/labs/phyloseq/lab_phyloseq.html
http://joey711.github.io/phyloseq/preprocess
https://github.com/microbiome/microbiome/blob/master/vignettes/vignette.md











### CORE MICROBIOME ANALYSIS
# https://github.com/microbiome/microbiome/tree/master/vignettes

library(microbiome)

#abundance boxplot of specific OTUs in relation to max 3 variables at a time
bp <-boxplot_abundance(biom_ITS, x = "Level", y = "OTU_117", line = "Level", color = "Level", log10 = TRUE)
bp

#core heatmap
biom_core_heatmapITS <- plot_core(Zygo, plot.type = "heatmap", palette = "spectral") #try also with biom_prop or biom_log or subset biom

subset_ITS_Level <- subset_samples(biom_ITS, Level%in%Soil2) 
subset_ITS_loam <- subset_samples(biom_ITS, Soil2=="loam")

biom_core_heatmapITS <- plot_core(subset_ITS_Level, plot.type = "heatmap", palette = "spectral") #try also with biom_prop or biom_log or subset biom


#abundance boxplot of specific OTUs in relation to max 3 variables at a time
boxplot_otu101<-boxplot_abundance(biom_ITS, x = "Level", y = "OTU_101", line = "Plant", color = "Soil2", log10 = TRUE)
boxplot_otu101

#core heatmap
biom_core_heatmap<-plot_core(Zygo, plot.type = "heatmap", palette = "spectral") #try also with biom_prop or biom_log or subset biom


#list of core OTUs with a given detection threshold
otus_prevalence<-prevalence(biom_ITS, detection.threshold = 10, sort = TRUE)
otus_prevalence

#list the OTUs present over in a given fraction of the samples
otu_core_50<-prevalent_taxa(biom, detection.threshold = 10, prevalence.threshold = 0.5)


df$taxa <- c(rep("Mortierella", times=3), rep("Ilionectrya", times=3), rep("control", times=3))

test








### ALPHA DIVERSITY -------------------------------------------------------------------------------------------------------------

# Changing order of levels and transform into a character


otu_ITS_cv <- as.data.frame(otu_table(biom_ITS_cv))
metadata_ITS_cv <- as.data.frame(as.matrix(sample_data(biom_ITS_cv)[,4:9]))

sums_ITS_cv <- data.frame(colSums(otu_table(biom_ITS_cv)))
colnames(sums_ITS_cv) <- "readNO"
sums_ITS_cv$richness <- specnumber(otu_ITS_cv, MARGIN = 2)
sums_ITS_cv$rarefied <- rarefy(otu_ITS_cv,2189, se = FALSE, MARGIN = 2)
sums_ITS_cv$site <- metadata_ITS_cv$Site1
sums_ITS_cv$Description <- row.names(sums_ITS_cv)
#sums_ITS_cv <- arrange(sums_ITS_cv, sample)
#sums_ITS_cv$group <- c(rep("Forest", times=34),rep("Poplar", times=27))
sums_ITS_cv
str(sums_ITS_cv)

p = ggplot(sums_ITS_cv, aes(x=site, y=rarefied), color=site) 

p + geom_boxplot(outlier.colour="black", outlier.shape=8, outlier.size=4) + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=4, col="red") +
  theme(axis.text.x = element_text(angle = -90, hjust=0, size=10)) +
  scale_colour_manual(values=palette_CB9) + 
labs(title="OTU", x="Habitat", y="Rarefied richness")






otu_16s_cv <- as.data.frame(otu_table(biom_16s_cv))
metadata_16s_cv <- as.data.frame(as.matrix(sample_data(biom_16s_cv)[,4:9]))

sums_16s_cv <- data.frame(colSums(otu_table(biom_16s_cv)))
colnames(sums_16s_cv) <- "readNO"
sums_16s_cv$richness <- specnumber(otu_16s_cv, MARGIN = 2)
sums_16s_cv$rarefied <- rarefy(otu_16s_cv,17479, se = FALSE, MARGIN = 2)
sums_16s_cv$site <- metadata_16s_cv$Site1
sums_16s_cv$Description <- row.names(sums_16s_cv)
#sums_16s_cv <- arrange(sums_16s_cv, sample)
#sums_16s_cv$group <- c(rep("Forest", times=34),rep("Poplar", times=27))
sums_16s_cv
str(sums_16s_cv)












plot_alpha_its





sample_data(biom_ITS_cv)
str(biom_ITS_cv)

plot_alpha_its <- plot_richness(biom_ITS_cv, x = "Site1", color = "Site1", shape= "Site1", measures = c("Observed","Shannon"))
plot_alpha_its <- plot_richness(biom_ITS_cv, x = "Site2", color = "Site2", measures = c("Observed"))

plot_alpha_its + 
geom_boxplot(outlier.colour="black", outlier.fill = "black") + 
geom_point(size=2, alpha = 0.9) + 
scale_colour_manual(values = c("Beach_Intertidal1" = "red",
                               "Beach_Intertidal2" = "red",
                               "Cassava1" = "red",
                               "Cassava2" = "blue",
                               "Cassuarina1"="black",
                               "Cassuarina2"="red",
                               "Grassland2"="black",
                               "Mohogany1"="red",
                               "Mohogany2"="black",
                               "Mohogany3"="blue",
                               "Native_Forest1"="red",
                               "Native_Forest3"="blue",
                               "Pinus_carribeana"="black",
                               "Stream_Bank1"="blue",
                               "Stream_Bank2"="blue",
                               "Stream_Bank3"="blue",
                               "Sugarcane1"="black",
                               "Sugarcane2"="black")) +
scale_shape_manual(values=0:19) +
ggtitle("A. ITS") +
theme(axis.text.x = element_text(angle = 90, hjust=0.5, size=10)) +
theme(legend.position="none")

  
plot_alpha_its <- plot_richness(biom_ITS_cv, x = "Site1", color = "Site1", measures = c("Observed","Shannon"), title="(A)")
plot_alpha_lsu <- plot_richness(biom_LSU_cv, x = "Site1", color = "Site1", measures = c("Observed","Shannon"), title="(B)")
plot_alpha_16s <- plot_richness(biom_16s_cv, x = "Site1", color = "Site1", measures = c("Observed","Shannon"), title="(C)")

rich_its = plot_alpha_its + geom_boxplot(outlier.colour="black", outlier.shape = 8) +
  scale_colour_manual(values=palette_CB9) + 
  #scale_colour_manual(values=c("#0000FF","#0000FF","#0000FF","#0000FF","#0000FF","#0000FF",
   #                          "#0000FF","#0000FF","#0000FF")) +
  #stat_summary(fun.y=mean, geom="point", shape=18, size=2, col="red") + 
  theme(legend.position="none")

rich_lsu = plot_alpha_lsu + geom_boxplot(outlier.colour="black", outlier.shape = 8) +
  scale_colour_manual(values=palette_CB9) + 
  #scale_colour_manual(values=c("#0000FF","#0000FF","#0000FF","#0000FF","#0000FF","#0000FF",
   #                            "#0000FF","#0000FF","#0000FF")) +
  #stat_summary(fun.y=mean, geom="point", shape=18, size=2, col="red") + 
  theme(legend.position="none")

rich_16s = plot_alpha_16s + geom_boxplot(outlier.colour="black", outlier.shape = 8) +
  scale_colour_manual(values=palette_CB9) + 
  #scale_colour_manual(values=c("#0000FF","#0000FF","#0000FF","#0000FF","#0000FF","#0000FF",
  #                             "#0000FF","#0000FF","#0000FF")) +
  #stat_summary(fun.y=mean, geom="point", shape=18, size=2, col="red") +
  theme(legend.position="none")

# multiplot
source("../multiplot.R")
multiplot(rich_its, rich_lsu, rich_16s, cols=3)
dev.off()




### trying to add p-value to alpha div
data("soilrep")
erich <- estimate_richness(soilrep, measures = c("Observed", "Shannon", "Simpson", "Fisher"))
ttest <- t(sapply(erich, function(x) unlist(t.test(x~sample_data(soilrep)$warmed)[c("estimate","p.value","statistic","conf.int")])))
ttest


#http://evomics.org/wp-content/uploads/2016/01/phyloseq-Lab-01-Answers.html






# to explore the resutled onject 
str(plot_alpha_its)
plot_alpha_its$data
plot_alpha_its$layers
write.csv(plot_alpha_its$data, file = "observed_shannon_ITS.csv")


plot_alpha_lsu <- plot_richness(biom_LSU_cv, x = "Site1", color = "Site1", measures = c("Observed","Shannon"))

plot_alpha_lsu <- plot_richness(biom_LSU_cv, x = "Site2", color = "Site1", shape= "Site2", measures = c("Observed","Simpson"))
plot_alpha_lsu + geom_boxplot() + 
  scale_colour_manual(values=palette_CB9) + 
  geom_point(size=2, alpha = 0.7) + 
  scale_shape_manual(values=0:19) + 
  ggtitle("(B) Alpha Diversity LSU") +
  theme(axis.text.x = element_text(angle = 90, hjust=0.5, size=10)) +
  theme(legend.position="none")

write.csv(plot_alpha_lsu$data, file = "observed_shannon_LSU.csv")

plot_alpha_16s <- plot_richness(biom_16s_cv, x = "Site1", color = "Site1", measures = c("Observed","Shannon"))

plot_alpha_16s <- plot_richness(biom_16s_cv, x = "Site2", color = "Site1", shape= "Site2", measures = c("Observed","Simpson"))
plot_alpha_16s + geom_boxplot() + scale_colour_manual(values=palette_CB9) + geom_point(size=2, alpha = 0.7) + 
  scale_shape_manual(values=0:26) + ggtitle("(C) Alpha Diversity 16S") +
  theme(axis.text.x = element_text(angle = 90, hjust=0.5, size=10)) + 
  theme(legend.position="none")


write.csv(plot_alpha_16s$data, file = "observed_shannon_16S.csv")



#to plot more than 6 shapes
scale_shape_manual(values=1:nlevels(sample_data(biom_ITS_ev)$chemical))


pr <- plot_richness(biom_ITS_cv, x = "Level", color = "Level", measures = c("Observed","Shannon"))
pr + geom_boxplot(outlier.colour = "red") + scale_colour_manual(values=palette_CB4)

pr <- plot_richness(biom_ITS_cv, x = "Level", color = "Level", shape="Soil", measures = c("Observed", "Shannon" ))
pr + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size=12)) + scale_colour_manual(values=palette_CB4)


pr <- plot_richness(biom_16s_cv, x = "Level", color = "Level", shape= "Soil", measures = c("Observed","Shannon"))
pr + geom_boxplot() + scale_colour_manual(values=palette_CB4) + geom_point(size=2, alpha = 0.7) + # coord_flip() + ggtitle("Alpha Diversity Bacteria")
  theme(axis.text.x = element_text(angle = 90, hjust=0.5, size=10))

pr <- plot_richness(biom_16s_cv, x = "Level", color = "Level", measures = c("Observed","Shannon"))
pr + geom_boxplot(outlier.colour = "red") + scale_colour_manual(values=palette_CB4) +
  theme(axis.text.x = element_text(angle = 90, hjust=0.5, size=10))



# transform back to numeric to use as continuous variable 
sample_data(biom_ITS_cv)$Level <- as.numeric(as.character(sample_data(biom_ITS)$Level, levels=c("0","10","50","100")))
sample_data(biom_16s_cv)$Level <- as.numeric(as.character(sample_data(biom_16s)$Level, levels=c("0","10","50","100")))

pr <- plot_richness(biom_ITS_cv, x = "Level", color = "Level", shape="Soil", measures = c("Observed", "Shannon" ))
pr + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size=10)) + stat_smooth(method = lm) #+ geom_boxplot()

pr <- plot_richness(biom_ITS_cv, x = "Level", shape= "Soil", measures = c("Observed", "Shannon" ))
pr + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size=10)) + stat_smooth(method = lm)

pr <- plot_richness(biom_16s_cv, x = "Level", color = "Level", shape="Soil", measures = c("Observed", "Shannon" ))
pr + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size=10)) + stat_smooth(method = lm) #+ geom_boxplot()


pr <- plot_richness(biom_16s_cv, x = "Level", shape= "Soil", measures = c("Observed", "Shannon" ))
pr + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size=10)) + stat_smooth(method = lm)






# **********---------
# http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html



### BETA DIVERSITY in Phyloseq -------------------------------------------------------------------------------------------
write.csv(sample_data(biom_ITS_cv) , file = "map_chem.csv")
map_chem_new <- read.csv("map_chem_new.csv", header=T, row.names =1)
map_chem_new

# fungi ------

sample_data(biom_ITS_cv) <- map_chem_new
head(sample_data(biom_ITS_cv))

biom_ITS_ev = rarefy_even_depth(biom_ITS_cv, sample.size = 1000)


biom_ITS_cv -> biom_ITS_hell
otu_table(biom_ITS_hell) <- otu_table(decostand(otu_table(biom_ITS_cv), method = "hellinger"), taxa_are_rows=TRUE)
otu_table(biom_ITS_hell)
sample_data(biom_ITS_hell)

otu_ITS_hell <- as.data.frame(otu_table(biom_ITS_cv))
taxa_ITS_hell <- as.data.frame(as.matrix(tax_table(biom_ITS_hell)))
metadata_ITS_hell <- as.data.frame(as.matrix(sample_data(biom_ITS_hell)))

identical(colnames(otu_ITS_hell2), rownames(metadata_ITS_hell))


# transform factor in numeric -------

#as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
#indx <- sapply(metadata_ITS_hell, is.factor)
#breast[indx] <- lapply(breast[indx], function(x) as.numeric(as.character(x)))

metadata_ITS_hell$P <- as.numeric(as.character(metadata_ITS_hell$P))
metadata_ITS_hell$K <- as.numeric(as.character(metadata_ITS_hell$K))
metadata_ITS_hell$Ca <- as.numeric(as.character(metadata_ITS_hell$Ca))
metadata_ITS_hell$Mg <- as.numeric(as.character(metadata_ITS_hell$Mg))
metadata_ITS_hell$Zn <- as.numeric(as.character(metadata_ITS_hell$Zn))
metadata_ITS_hell$Mn <- as.numeric(as.character(metadata_ITS_hell$Mn))
metadata_ITS_hell$Fe <- as.numeric(as.character(metadata_ITS_hell$Fe))
metadata_ITS_hell$NO3 <- as.numeric(as.character(metadata_ITS_hell$NO3))
metadata_ITS_hell$NH4 <- as.numeric(as.character(metadata_ITS_hell$NH4))

str(metadata_ITS_hell)


# perform bioenv function --------------
bioenv_ITS <- vegan::bioenv(t(otu_ITS_hell) ~ P + K + Ca + Mg + Zn + Mn + Fe + NO3 + NH4, 
                     metadata_ITS_hell, method = "pearson", index = "bray")
summary(bioenv_ITS)
bioenv_ITS


CAP_ITS <- ordinate(biom_ITS_hell, method = "CAP", distance = "bray", ~  K + Ca + Mn + NH4)
CAP_ITS
anova(CAP_ITS)

plot_cap_ITS = plot_ordination(biom_ITS_ev, CAP_ITS, 
                       color = "Site2", shape="Site2", title = "A.") + 
  geom_point(size=3.2) + 
  geom_point(aes(fill=Site2), stroke = 1.5) +
  scale_shape_manual(values=setNames(c(0,1,2,6,4,3,7,8,9,10,11,12,13,14,15,18,16,20), sort(unique(sample_data(biom_ITS_hell)$Site2)))) + 
  theme_bw() + 
  #geom_text_repel(aes(label=X.SampleID), size = 3) +
  scale_colour_manual(values = c("Beach_Intertidal1" = "#771155", "Beach_Intertidal2" = "#771155","Cassava1" = "#CC99BB", "Cassava2" = "#CC99BB",
                                 "Cassuarina1"="#114477","Cassuarina2"="#114477", "Grassland2"="#4477AA","Mohogany1"="#117744",
                                 "Mohogany2"="#117744","Mohogany3"="#117744","Native_Forest1"="#777711","Native_Forest3"="#777711",
                                 "Pinus_carribeana"="#AA7744","Stream_Bank1"="#DDDD77","Stream_Bank2"="#DDDD77","Stream_Bank3"="#DDDD77",
                                 "Sugarcane1"="#E69F00","Sugarcane2"="#E69F00"))

plot_cap_ITS


# Now add the environmental variables as arrows ----------
arrowmat <- vegan::scores(CAP_ITS, display = "bp")
arrowmat

# Add labels, make a data.frame --------------------------
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping ---------------------
arrow_map <- aes(xend = CAP1, yend = CAP2,x = 0,y = 0,shape = NULL, color = "red", label = labels)
arrow_map
label_map <- aes(x = 1.3 * CAP1, y = 1.3 * CAP2, shape = NULL, color = NULL, label = labels)
label_map
arrowhead = arrow(length = unit(0.02, "npc"))
arrowhead

# Make a new graphic -------------------------------------
plot_cap_ITS_arrow = plot_cap_ITS + 
  geom_segment(mapping = arrow_map, size = .6, data = arrowdf,color = "gray", arrow = arrowhead) + 
  geom_text(mapping = label_map, size = 5, data = arrowdf, show.legend = FALSE)

plot_cap_ITS_arrow



# bacteria -------
identical(sample_data(biom_ITS_cv)$Site2, sample_data(biom_16s_cv)$Site2) 

write.csv(sample_data(biom_16s_cv) , file = "map_chem_16s.csv")
map_chem_new <- read.csv("map_chem_new_16s.csv", header=T, row.names =1)
map_chem_new


sample_data(biom_16s_cv) <- map_chem_new
head(sample_data(biom_16s_cv))

biom_16s_ev = rarefy_even_depth(biom_16s_cv, sample.size = 1000)


biom_16s_cv -> biom_16s_hell
otu_table(biom_16s_hell) <- otu_table(vegan::decostand(otu_table(biom_16s_cv), method = "hellinger"), taxa_are_rows=TRUE)
otu_table(biom_16s_hell)
sample_data(biom_16s_hell)

otu_16s_hell <- as.data.frame(otu_table(biom_16s_cv))
taxa_16s_hell <- as.data.frame(as.matrix(tax_table(biom_16s_hell)))
metadata_16s_hell <- as.data.frame(as.matrix(sample_data(biom_16s_hell)))

# transform factor in numeric -------

#as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
#indx <- sapply(metadata_16s_hell, is.factor)
#breast[indx] <- lapply(breast[indx], function(x) as.numeric(as.character(x)))

metadata_16s_hell$P <- as.numeric(as.character(metadata_16s_hell$P))
metadata_16s_hell$K <- as.numeric(as.character(metadata_16s_hell$K))
metadata_16s_hell$Ca <- as.numeric(as.character(metadata_16s_hell$Ca))
metadata_16s_hell$Mg <- as.numeric(as.character(metadata_16s_hell$Mg))
metadata_16s_hell$Zn <- as.numeric(as.character(metadata_16s_hell$Zn))
metadata_16s_hell$Mn <- as.numeric(as.character(metadata_16s_hell$Mn))
metadata_16s_hell$Fe <- as.numeric(as.character(metadata_16s_hell$Fe))
metadata_16s_hell$NO3 <- as.numeric(as.character(metadata_16s_hell$NO3))
metadata_16s_hell$NH4 <- as.numeric(as.character(metadata_16s_hell$NH4))

str(metadata_16s_hell)


# perform bioenv function --------------
bioenv_16s <- vegan::bioenv(t(otu_16s_hell) ~ P + K + Ca + Mg + Zn + Mn + Fe + NO3 + NH4, 
                            metadata_16s_hell, method = "pearson", index = "bray")
summary(bioenv_16s)
bioenv_16s

CAP_16s <- ordinate(biom_16s_hell, method = "CAP", distance = "bray", ~  P + K + Ca + Mn + Fe + NH4)
CAP_16s
anova(CAP_16s)


plot_cap_16s = plot_ordination(biom_16s_hell, CAP_16s, color = "Site2", shape="Site2", title = "B.") +
  geom_point(size=3.2) + 
  geom_point(aes(fill=Site2), stroke = 1.5) +
  scale_shape_manual(values=setNames(c(0,1,2,6,4,3,7,8,9,10,11,12,13,14,15,18,16,20), sort(unique(sample_data(biom_16s_hell)$Site2)))) + 
  theme_bw() + 
  #geom_text_repel(aes(label=Description), size = 3) +
  scale_colour_manual(values = c("Beach_Intertidal1" = "#771155",
                                 "Beach_Intertidal2" = "#771155","Cassava1" = "#CC99BB", "Cassava2" = "#CC99BB",
                                 "Cassuarina1"="#114477","Cassuarina2"="#114477", "Grassland2"="#4477AA","Mohogany1"="#117744",
                                 "Mohogany2"="#117744","Mohogany3"="#117744","Native_Forest1"="#777711","Native_Forest3"="#777711",
                                 "Pinus_carribeana"="#AA7744","Stream_Bank1"="#DDDD77","Stream_Bank2"="#DDDD77","Stream_Bank3"="#DDDD77",
                                 "Sugarcane1"="#E69F00","Sugarcane2"="#E69F00"))
plot_cap_16s


# Now add the environmental variables as arrows ----------
arrowmat <- vegan::scores(CAP_16s, display = "bp")
arrowmat

# Add labels, make a data.frame --------------------------
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping ---------------------
arrow_map <- aes(xend = CAP1,yend = CAP2,x = 0,y = 0,shape = NULL, color = "red", label = labels)
arrow_map
label_map <- aes(x = 1.3 * CAP1, y = 1.3 * CAP2, shape = NULL, color = NULL, label = labels)
label_map
arrowhead = arrow(length = unit(0.02, "npc"))
arrowhead

# Make a new graphic -------------------------------------
plot_cap_16s_arrow = plot_cap_16s + 
  geom_segment(mapping = arrow_map, size = .6, data = arrowdf,color = "gray", arrow = arrowhead) + 
  geom_text(mapping = label_map, size = 5, data = arrowdf, show.legend = FALSE)


plot_cap_16s_arrow


library("ggpubr")
ggarrange(plot_cap_ITS_arrow, plot_cap_16s_arrow,
          widths = c(1,1),
          align = "hv", ncol = 2, nrow = 1,
          legend="bottom", common.legend = TRUE)




# >>> CLUSTERING  ------------
chemical_data <- as.data.frame(as.matrix(sample_data(biom_ITS_hell)[c(1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52),c(6,10:18)]))                
chemical_data        

chemical_data <- as.data.frame(as.matrix(sample_data(biom_ITS_hell)[,c(6,10:18)]))                
chemical_data  

chemical_data$Site1 <- sample_data(biom_ITS_hell)$Site1


distance_soil <- dist(as.matrix(chemical_data))   # find distance matrix 
hc <- hclust(distance_soil, method = "ward.D")  # apply hirarchical clustering 
plot(hc)                       # plot the dendrogram     
                
plot(hc, cex = 0.8, hang = -1)


write.csv(chemical_data, "Table_soil.csv")







# fungi ITS
#NMDS_ITS_ev <- ordinate(biom_ITS_ev, method ="NMDS", distance="bray", try=100)

#nmds_ev = plot_ordination(biom_ITS_ev, NMDS_ITS_ev, type="sites", color="Site2", shape="Site2", title=" (A) NMDS Fungi ITS - stress 0.169")

#nmds_ev + geom_point(size=3.2) + 
 # geom_point(aes(fill=Site2), stroke = 1.5) +
#  scale_shape_manual(values=setNames(c(0,1,2,6,4,3,7,8,9,10,11,12,13,14,15,18,16,20), sort(unique(sample_data(biom_ITS_hell)$Site2)))) + 
 # theme_bw() + geom_text_repel(aes(label=X.SampleID), size = 3) +
#  scale_colour_manual(values = c("Beach_Intertidal1" = "#771155", "Beach_Intertidal2" = "#771155","Cassava1" = "#CC99BB", "Cassava2" = "#CC99BB",
#                                 "Cassuarina1"="#114477","Cassuarina2"="#114477", "Grassland2"="#4477AA","Mohogany1"="#117744",
 #                                "Mohogany2"="#117744","Mohogany3"="#117744","Native_Forest1"="#777711","Native_Forest3"="#777711",
  #                               "Pinus_carribeana"="#AA7744","Stream_Bank1"="#DDDD77","Stream_Bank2"="#DDDD77","Stream_Bank3"="#DDDD77",
  #                               "Sugarcane1"="#E69F00","Sugarcane2"="#E69F00"))

library("ggrepel")

NMDS_ITS_hell <- ordinate(biom_ITS_hell, method ="NMDS", distance="bray", try=100)

sample_data(biom_ITS_hell)$Site2  <- with(biom_ITS_hell, reorder(sample_data(biom_ITS_hell)$Site2)) # to reorder before plotting
sample_data(biom_ITS_hell)

NMDS_ITS_hell
nmds_hell = plot_ordination(biom_ITS_hell, NMDS_ITS_hell, type="sites", color="Site2", shape="Site2", title=" (A) NMDS Fungi ITS - stress 0.169")

nmds_hell + geom_point(size=3.2) + 
geom_point(aes(fill=Site2), stroke = 1.5) +
scale_shape_manual(values=setNames(c(0,1,2,6,4,3,7,8,9,10,11,12,13,14,15,18,16,20), sort(unique(sample_data(biom_ITS_hell)$Site2)))) + 
theme_bw() + geom_text_repel(aes(label=X.SampleID), size = 3) +
scale_colour_manual(values = c("Beach_Intertidal1" = "#771155", "Beach_Intertidal2" = "#771155","Cassava1" = "#CC99BB", "Cassava2" = "#CC99BB",
                                 "Cassuarina1"="#114477","Cassuarina2"="#114477", "Grassland2"="#4477AA","Mohogany1"="#117744",
                                 "Mohogany2"="#117744","Mohogany3"="#117744","Native_Forest1"="#777711","Native_Forest3"="#777711",
                                 "Pinus_carribeana"="#AA7744","Stream_Bank1"="#DDDD77","Stream_Bank2"="#DDDD77","Stream_Bank3"="#DDDD77",
                                 "Sugarcane1"="#E69F00","Sugarcane2"="#E69F00"))
  


#scale_shape_manual(values=0:25) + 
#+ theme(plot.title = element_text(hjust = 0.5)) # to center the title
#+ scale_shape_manual(values=setNames(c(0,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,18), sort(unique(sample_data(biom_ITS_hell)$Site2)))) # Assign value pairs for scale_color_manual and scale_shape_manual in ggplot
#scale_color_manual(values=custom_col21) +

pcoa_ITS <- ordinate(biom_ITS_hell, method ="PCoA", distance="bray")
pcoa_ITS <- ordinate(biom_ITS_hell, method ="PCoA", distance(biom_ITS_hell, method = "jaccard", binary = TRUE))

p_pcoa_ITS = plot_ordination(biom_ITS_hell, pcoa_ITS, type="sites", color="Site2", shape="Site2", title=" (A) PCoA Fungi ITS")
p_pcoa_ITS + geom_point(size=3) + geom_point(aes(fill=Site2), stroke = 1.5) + scale_shape_manual(values=0:25) + 
  scale_color_manual(values=custom_col21) + theme_bw() + geom_text_repel(aes(label=X.SampleID), size = 3)


###
NMDS_LSU_hell <- ordinate(biom_LSU_hell, method ="NMDS", distance="bray", try=100)
NMDS_LSU_hell <- ordinate(biom_LSU_hell, method ="NMDS", distance(biom_LSU_hell, method = "jaccard", binary = TRUE), try=100)

NMDS_LSU_hell
nmds_LSU = plot_ordination(biom_LSU_hell, NMDS_LSU_hell, type="sites", color="Site2", shape="Site2", title=" (B) NMDS Fungi LSU - stress 0.159")

nmds_LSU + geom_point(size=3.2) + 
geom_point(aes(fill=Site2), stroke = 1.5) + 
scale_shape_manual(values=setNames(c(0,1,2,6,4,7,8,9,10,11,12,13,14,15,18,16,20), sort(unique(sample_data(biom_LSU_hell)$Site2)))) + 
theme_bw() + geom_text_repel(aes(label=X.SampleID), size = 3) +
scale_colour_manual(values = c("Beach_Intertidal1" = "#771155","Beach_Intertidal2" = "#771155","Cassava1" = "#CC99BB", "Cassava2" = "#CC99BB",
                                 "Cassuarina1"="#114477", "Grassland2"="#4477AA","Mohogany1"="#117744",
                                 "Mohogany2"="#117744","Mohogany3"="#117744","Native_Forest1"="#777711","Native_Forest3"="#777711",
                                 "Pinus_carribeana"="#AA7744","Stream_Bank1"="#DDDD77","Stream_Bank2"="#DDDD77","Stream_Bank3"="#DDDD77",
                                 "Sugarcane1"="#E69F00","Sugarcane2"="#E69F00"))


nmds_lsu = nmds_LSU + geom_point(size=3) + geom_point(aes(fill=Site2), stroke = 1.5) + #scale_shape_manual(values=0:25) + 
  scale_shape_manual(values=setNames(c(0,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17), sort(unique(sample_data(biom_LSU_hell)$Site2)))) + #scale_color_manual(values=custom_col21) + 
  scale_color_manual(values=c("#771155","#AA4488", "#CC99BB", "#114477", "#4477AA", "#117777", "#44AAAA", "#77CCCC", "#117744", 
                              "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455","#DD7788")) +
  theme_bw() + geom_text_repel(aes(label=X.SampleID), size = 3) 





pcoa_LSU <- ordinate(biom_LSU_hell, method ="PCoA", distance="bray")
pcoa_LSU <- ordinate(biom_LSU_hell, method ="PCoA", distance(biom_LSU_hell, method = "jaccard", binary = TRUE))

p_prop_LSU = plot_ordination(biom_LSU_hell, pcoa_LSU, type="sites", color="Site2", shape="Site2", title=" (B) PCoA Fungi LSU")
p_prop_LSU + geom_point(size=3) + geom_point(aes(fill=Site2), stroke = 1.5) + scale_shape_manual(values=0:25) + 
  scale_color_manual(values=custom_col21) + theme_bw() + geom_text_repel(aes(label=X.SampleID), size = 3)


###
NMDS_16s_hell <- ordinate(biom_16s_hell, method ="NMDS", distance="bray", try=200)
NMDS_16s_hell <- ordinate(biom_16s_hell, method ="NMDS", distance(biom_16s_hell, method = "jaccard", binary = TRUE), try=100)

NMDS_16s_hell
bact_hell = plot_ordination(biom_16s_hell, NMDS_16s_hell, type="sites", color="Site2", shape="Site2", title=" (C) NMDS Bacteria 16S - stress 0.08")

nmds_16s = bact_hell + geom_point(size=3.2) + 
geom_point(aes(fill=Site2), stroke = 1.5) + 
scale_shape_manual(values=setNames(c(0,1,2,6,4,3,7,8,9,10,11,12,13,14,15,18,16,20), sort(unique(sample_data(biom_16s_hell)$Site2)))) + 
theme_bw() + geom_text_repel(aes(label=Description), size = 3) +
scale_colour_manual(values = c("Beach_Intertidal1" = "#771155",
                                   "Beach_Intertidal2" = "#771155","Cassava1" = "#CC99BB", "Cassava2" = "#CC99BB",
                                   "Cassuarina1"="#114477","Cassuarina2"="#114477", "Grassland2"="#4477AA","Mohogany1"="#117744",
                                   "Mohogany2"="#117744","Mohogany3"="#117744","Native_Forest1"="#777711","Native_Forest3"="#777711",
                                   "Pinus_carribeana"="#AA7744","Stream_Bank1"="#DDDD77","Stream_Bank2"="#DDDD77","Stream_Bank3"="#DDDD77",
                                   "Sugarcane1"="#E69F00","Sugarcane2"="#E69F00"))

nmds_16s

pcoa_16s <- ordinate(biom_16s_hell, method ="PCoA", distance="bray")
pcoa_16s <- ordinate(biom_16s_hell, method ="PCoA", distance(biom_16s_hell, method = "jaccard", binary = TRUE))

p_prop_16s = plot_ordination(biom_16s_hell, pcoa_16s, type="sites", color="Site2", shape="Site2", title="(C) PCoA Bacteria 16S")
p_prop_16s + geom_point(size=3) + geom_point(aes(fill=Site2), stroke = 1.5) + scale_shape_manual(values=0:25) + 
  scale_color_manual(values=custom_col21) + theme_bw() + geom_text_repel(aes(label=Description), size = 3)


# multiplot
source("../multiplot.R")
multiplot(nmds_its, nmds_lsu, nmds_16s, cols=3)
dev.off()



# cap analysis
CAP_ITS <- ordinate(biom_ITS_hell, method ="CAP", distance="bray", ~ Site2) 
CAP_ITS = plot_ordination(biom_ITS_hell, CAP_ITS, type="sites", color="Site2", shape="Site2", title=" CAP Fungi ITS")
CAP_ITS + geom_point(size=3) + geom_point(aes(fill=Site2), stroke = 1.5) + 
  scale_shape_manual(values=0:25) + scale_color_manual(values=custom_col21) + theme_bw()

# trying to do plot the same ordinations with Phyloseq
RDA_ITS <- ordinate(biom_ITS_prop, method = "RDA", ~ Site2)
RDA_ITS = plot_ordination(biom_ITS_prop, RDA_ITS, type="sites", color="Site2", shape="Site2", title=" RDA Fungi ITS")
RDA_ITS + geom_point(size=3) + geom_point(aes(fill=Site2), stroke = 1.5) + scale_shape_manual(values=0:25) + 
  scale_color_manual(values=custom_col21) + theme_bw() + geom_text_repel(aes(label=X.SampleID), size = 3)

RDA_LSU <- ordinate(biom_LSU_hell, method = "RDA", ~ Site2)
RDA_LSU = plot_ordination(biom_LSU_hell, RDA_LSU, type="sites", color="Site2", shape="Site2", title=" RDA Fungi LSU")
RDA_LSU + geom_point(size=3) + geom_point(aes(fill=Site2), stroke = 1.5) + scale_shape_manual(values=0:25) + 
  scale_color_manual(values=custom_col21) + theme_bw() + geom_text_repel(aes(label=X.SampleID), size = 3)




# CALCULATING NETWEORK network with the merged fasta
# still working on MERGING FUNGI AND BACTERIA

biom_ITS_net <- make_network(biom_ITS_cv, dist.fun = "jaccard", max.dist = 0.90, line_weight = 0.2, label=NULL)
biom_ITS_net_plot <-plot_network(biom_ITS_net, biom_ITS_cv, color="Site2", shape="Site2")
biom_ITS_net_plot + scale_colour_manual(values=custom_col21) + scale_shape_manual(values=0:25) + 
  theme(plot.title = element_text(hjust = 0.5)) + ggtitle("(C) Network Fungi ITS")

biom_LSU_net <- make_network(biom_LSU_cv, dist.fun = "jaccard", max.dist = 0.80, line_weight = 0.2, label=NULL)
biom_LSU_net_plot <-plot_network(biom_LSU_net, biom_LSU_cv, color="Site2", shape="Site2")
biom_LSU_net_plot + scale_colour_manual(values=custom_col21) + scale_shape_manual(values=0:25) + 
  theme(plot.title = element_text(hjust = 0.5)) + ggtitle("(C) Network Fungi LSU")



## new version
biom_ITS_net_plot2 <- plot_net(biom_ITS_cv, distance = "jaccard", type = "samples", maxdist = 0.90, color = "Site2", shape = "Site2", rescale=TRUE)

biom_ITS_net_plot2 + scale_colour_manual(values=custom_col21) + scale_shape_manual(values=0:25) + 
  theme(plot.title = element_text(hjust = 0.5)) + ggtitle("(C) Network Fungi ITS")




## barplots using merged data
merged_ITS
otu_table(merged_ITS)
tax_table(merged_ITS)
sample_data(merged_ITS)


plot_bar(merged_ITS, x="Site1", fill="Class")

biom_ITS_des = merge_samples(merged_ITS, "Site1") # merging samples
sample_data(biom_ITS_des)$Site1 <- levels(sample_data(merged_ITS)$Site1) # relabelling 
biom_ITS_des = transform_sample_counts(biom_ITS_des, function(x) 100 * x/sum(x)) # calculating abundances
plot_bar(biom_ITS_des, x="Site1", fill="Class")

top20ITS <- names(sort(taxa_sums(biom_ITS_des), TRUE)[1:20]) 
biom_ITS_des <- prune_taxa(top20ITS, biom_ITS_des)
biom_ITS_des

p1 <- plot_ordered_bar(biom_ITS_des, x = "Site1", fill="Genus", leg_size = 0.4, title="Fungi ITS Barplots")
p1 + scale_fill_manual(values=custom_col21) + ylim(0, 90)

tax_table(biom_ITS_des)

##################################################################

##################################################################

##################################################################
### Analysing mycorrhizal fungi only - extracted with fungild ####




biom_ITS_ecm = import_biom("otu_table_ITS_LAST_guilds_json.biom")
map_ITS = import_qiime_sample_data("mapping_ITS.txt")
sample_data(biom_ITS_ecm) <- map_ITS
colnames(tax_table(biom_ITS_ecm)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
otus_rep_ITS_ecm <- readDNAStringSet("otus_ITS.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
biom_ITS_ecm <- merge_phyloseq(biom_ITS_ecm, otus_rep_ITS_ecm)
biom_ITS_ecm
sample_data(biom_ITS_ecm)
head(tax_table(biom_ITS_ecm))
refseq(biom_ITS_ecm)


biom_ITS_pat = import_biom("otu_table_ITS_pathogen_json.biom")
map_ITS = import_qiime_sample_data("mapping_ITS.txt")
sample_data(biom_ITS_pat) <- map_ITS
colnames(tax_table(biom_ITS_pat)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
otus_rep_ITS_pat <- readDNAStringSet("otus_ITS.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
biom_ITS_pat <- merge_phyloseq(biom_ITS_pat, otus_rep_ITS_pat)
biom_ITS_pat
sample_data(biom_ITS_pat)
head(tax_table(biom_ITS_pat))
refseq(biom_ITS_pat)




#tax_table(biom_ITS_ecm)
#tax_table(biom_ITS_ecm) <- tax_table(biom_ITS_ecm)[-which(row.names(tax_table(biom_ITS_ecm)) == "OTU_1162"), ] 
#otu_table(biom_ITS_ecm) <- otu_table(biom_ITS_ecm)[which(rowSums(otu_table(biom_ITS_ecm)) > 2),] 

biom_ITS_ecm_gen = tax_glom(biom_ITS_ecm, "Genus")
biom_ITS_ecm_gen
as.data.frame(otu_table(biom_ITS_ecm_gen))
tax_table(biom_ITS_ecm_gen)

biom_ITS_pat_gen = tax_glom(biom_ITS_pat, "Genus")
biom_ITS_pat_gen
as.data.frame(otu_table(biom_ITS_pat_gen))
tax_table(biom_ITS_pat_gen)



write.csv(as.data.frame(otu_table(biom_ITS_ecm_gen)), file = "ecm_noabund.csv")
write.csv(as.data.frame(tax_table(biom_ITS_ecm_gen)), file = "ecm_noabund_tax.csv")


read.table("ecm_otus.csv", header=TRUE, row.names = 1, sep = ",") -> ecm_otus
read.table("endophyte_otus.csv", header=TRUE, row.names = 1, sep = ",") -> endo_otus
ecm_otus
endo_otus

as.matrix(ecm_otus)
otu_table(biom_ITS_ecm_gen)

otu_table(biom_ITS_ecm_gen) <- otu_table(as.matrix(endo_otus), taxa_are_rows = TRUE)
otu_table(biom_ITS_ecm_gen)
tax_table(biom_ITS_ecm_gen)
biom_ITS_ecm_gen

# that's good for OTUs abundances 
biom_ITS_ecm_gen_ab = transform_sample_counts(biom_ITS_ecm_gen, function(x) 100* x/sum(x))
biom_ITS_ecm_gen_ab

# if I want phylum abundances?
#colSums(as.data.frame(otu_table(biom_ITS_ecm_gen)))

biom_ITS_ecm_gen = tax_glom(biom_ITS_ecm, "Genus")
biom_ITS_ecm_gen

otu_ITS_ecm_gen = taxa_sums(biom_ITS_ecm_gen)/sum(taxa_sums(biom_ITS_ecm_gen))*100
otu_ITS_ecm_gen

tax_ITS_ecm_gen <- data.frame(tax_table(biom_ITS_ecm_gen))
tax_ITS_ecm_gen <- tax_ITS_ecm_gen[c(1:6)]
tax_ITS_ecm_gen

tax_ITS_ecm_gen$abundance <- as.vector(otu_ITS_ecm_gen)
tax_ITS_ecm_gen <- tax_ITS_ecm_gen[order(tax_ITS_ecm_gen$abundance, decreasing = TRUE),] 
tax_ITS_ecm_gen



#### pathogens
biom_ITS_pat_gen = tax_glom(biom_ITS_pat, "Genus")
biom_ITS_pat_gen

otu_ITS_pat_gen = taxa_sums(biom_ITS_pat_gen)/sum(taxa_sums(biom_ITS_pat_gen))*100
otu_ITS_pat_gen

tax_ITS_pat_gen <- data.frame(tax_table(biom_ITS_pat_gen))
tax_ITS_pat_gen <- tax_ITS_pat_gen[c(1:6)]
tax_ITS_pat_gen

tax_ITS_pat_gen$abundance <- as.vector(otu_ITS_pat_gen)
tax_ITS_pat_gen <- tax_ITS_pat_gen[order(tax_ITS_pat_gen$abundance, decreasing = TRUE),] 
tax_ITS_pat_gen


###
library(vegan)
data(dune)
dune
rowSums(sweep(dune, 1, rowSums(dune), "/"))
###



sample_data(biom_ITS_ecm)
head(otu_table(biom_ITS_ecm))

#biom_ITS_ecm_gen = tax_glom(biom_ITS_ecm, "Genus")
biom_ITS_ecm_gen
biom_ITS_ecm_gen_merg = merge_samples(biom_ITS_ecm_gen, "Site1")
biom_ITS_ecm_gen_merg
sample_data(biom_ITS_ecm_gen_merg)
otu_table(biom_ITS_ecm_gen_merg)

sample_data(biom_ITS_ecm_gen_merg)$Site1 <- levels(sample_data(biom_ITS_ecm_gen)$Site1)
biom_ITS_ecm_gen_merg_ab = transform_sample_counts(biom_ITS_ecm_gen_merg, function(x) 100 * x/sum(x))
biom_ITS_ecm_gen_merg_ab
otu_table(biom_ITS_ecm_gen_merg_ab)
tax_table(biom_ITS_ecm_gen_merg_ab)
sum(otu_table(biom_ITS_ecm_gen_merg_ab))


write.csv(tax_table(biom_ITS_ecm_gen_merg_ab), file = "biom_ITS_ecm_gen_merg_ab.csv")
tax_tab_ecm_gen_merg_ab <- read.csv("biom_ITS_ecm_gen_merg_ab.csv", header=T, row.names =1)
tax_tab_ecm_gen_merg_ab
tax_table(biom_ITS_ecm_gen_merg_ab) <- tax_table(as.matrix(tax_tab_ecm_gen_merg_ab))


source("../plot_ordered_bar_byGian.R")
p1 <- plot_ordered_bar(biom_ITS_ecm_gen_merg_ab, x = "Site1", fill="Genus", leg_size = 0.4, title="B.")
p1


p1 + scale_fill_manual(values=custom_col21) + ylim(0, 100) #+ coord_flip() #+ theme_bw()
p1 + scale_fill_manual(values=custom_col42) + ylim(0, 75) + facet_wrap(~Genus) #+ coord_flip() #+ theme_bw()


p1 + scale_color_manual(values=custom_col21) + ylim(0, 100) #+ coord_flip() #+ theme_bw()
p1 + facet_wrap(~Site2) + ggtitle("Faceting is better than stacking...") + ylim(0, 100)
p1 + geom_bar(aes(color=Genus,fill=Genus), stat="identity", position="stack") + ylim(0, 50)
p1 + scale_fill_manual(values=custom_col21) + coord_flip()
p1 + scale_colour_brewer(palette = "Set1") + coord_flip()
p1 + scale_fill_manual(values=custom_col21) + geom_bar(aes(color=Genus,fill=Genus), stat="identity", position="stack")
p1 + scale_fill_manual(values=custom_col21) + geom_bar(stat = "identity")

p1 + theme(legend.text = element_text(face="italic")) +
  scale_fill_manual(values = c("Acaulospora" = "#560d0d",
                               "Amanita" = "#a35151",
                               "Ambispora" = "#dba4a4", 
                               "Cenococcum" = "#cc1c1c",
                               "Ceratobasidium"="#111b77",
                               "Chloridium"="#283dff", 
                               "Chromelosporium"="#636bb7",
                               "Claroideoglomus"="#bfc5ff",
                               "Corymbiglomus"="#014443",
                               "Dentiscutata"="#195637",
                               "Diversispora"="#117744",
                               "Entoloma"="#60ffaf",
                               "Gigaspora"="#b7ffdb",
                               "Glomus"="#825121",
                               "Inocybe"="#ea7f17",
                               "Oidiodendron"="#fcb067",
                               "Paraglomus"="#ffe8d3",
                               "Peziza"="#d8d6d4",
                               "Racocetra"="#82807f",
                         "Ramaria" = "#3f3e3d",
                         "Rhizophagus" = "#000000",
                         "Rhizopogon" = "#5b5b19", 
                         "Scleroderma" = "#fcfc00",
                         "Scutellospora"="#ffff9e",
                         "Septoglomus"="#ffb7ef", 
                         "Serendipita"="#fa7efc",
                         "Stephanospora"="#ae09ea",
                         "Tomentella"="#521899",
                         "Wilcoxina"="#1e0047",
                         "unidentified1"="white",
                         "unidentified2"="white",
                         "unidentified2"="white")) + 
                        theme(panel.background = element_blank()) +
  guides(fill=guide_legend(ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  theme(legend.position="right")


p1 + theme(legend.text = element_text(face="italic")) +
  scale_fill_manual(values = c("Annulohypoxylon" = "#560d0d",
                               "Anthostomella" = "#a35151",
                               "Camarosporium" = "#dba4a4", 
                               "Capronia" = "#cc1c1c",
                               "Chalara"="#111b77",
                               "Cladosporium"="#283dff", 
                               "Colletotrichum"="#636bb7",
                               "Cytospora"="#bfc5ff",
                               "Diaporthe"="#014443",
                               "Kabatiella"="#195637",
                               "Lecythophora"="#117744",
                               "Microdochium"="#60ffaf",
                               "Paecilomyces"="#b7ffdb",
                               "Periconia"="#825121",
                               "Phialophora"="#ea7f17",
                               "Phomopsis"="#fcb067",
                               "Torula"="#ffe8d3",
                               "Trichoderma"="#d8d6d4",
                               "Williopsis"="#82807f",
                               "Xylaria" = "#3f3e3d")) + 
                              theme(panel.background = element_blank())

##### plot pathogen fungi

biom_ITS_pat_gen_merg = merge_samples(biom_ITS_pat_gen, "Site1")
biom_ITS_pat_gen_merg_ab = transform_sample_counts(biom_ITS_pat_gen_merg, function(x) 100 * x/sum(x))

badTaxa = c("OTU_1494","OTU_2750","OTU_2559","OTU_3578","OTU_3405","OTU_3221","OTU_5205","OTU_3915","OTU_2624","OTU_3451"
             ,"OTU_3118","OTU_6904","OTU_4023","OTU_5200","OTU_5862","OTU_4387","OTU_5089","OTU_4523","OTU_4615","OTU_4945")

pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

biom_ITS_pat_gen_merg_ab_prun = pop_taxa(biom_ITS_pat_gen_merg_ab, badTaxa)
biom_ITS_pat_gen_merg_ab_prun

#top26 <- names(sort(taxa_sums(biom_ITS_pat_gen_merg_ab), TRUE)[1:26]) # sampling first 100 taxa
#biom_ITS_pat_gen_merg_ab_26 <- prune_taxa(top26, biom_ITS_pat_gen_merg_ab)

p1 <- plot_ordered_bar(biom_ITS_pat_gen_merg_ab_prun, x = "Site1", fill="Genus", leg_size = 0.4, title="C.")
p1 + scale_fill_manual(values=custom_col42) + ylim(0, 100) #+ coord_flip() #+ theme_bw()

p1 + theme(legend.text = element_text(face="italic")) +
  scale_fill_manual(values = c("Aspergillus" = "#560d0d",
                               "Calonectria" = "#a35151",
                               "Ceratocystis" = "#dba4a4", 
                               "Clonostachys" = "#cc1c1c",
                               "Collophora"="#111b77",
                               "Coniothyrium"="#283dff", 
                               "Curvularia"="#636bb7",
                               "Cylindrocarpon"="#bfc5ff",
                               "Cylindrocladiella"="#014443",
                               "Endomelanconiopsis"="#195637",
                               "Fusarium"="#117744",
                               "Ganoderma"="#60ffaf",
                               "Gibberella"="#b7ffdb",
                               "Gibellulopsis"="#825121",
                               "Glomerella"="#ea7f17",
                               "Haematonectria"="#fcb067",
                               "Lasiodiplodia"="#ffe8d3",
                               "Magnaporthe"="#d8d6d4",
                               "Monographella"="#82807f",
                               "Mycoleptodiscus" = "#3f3e3d",
                               "Pestalotiopsis" = "#000000",
                               "Phaeoacremonium" = "#5b5b19", 
                               "Phoma" = "#fcfc00",
                               "Pyricularia"="#ffff9e",
                               "Rhizophydium"="#ffb7ef", 
                               "Sclerotinia"="#fa7efc",
                               "Spizellomyces"="#ae09ea",
                               "Thanatephorus"="#521899",
                               "Trichothecium"="#a0fffc",
                               "Volutella"="#1e0047")) + 
                      theme(panel.background = element_blank())

                      
                        



### other plot
p1 + theme(legend.text = element_text(face="italic")) +
  scale_fill_manual(values = c("Amanita" = "#560d0d",
                               "Entoloma" = "#a35151",
                               "Hygrocybe" = "#cc1c1c", 
                               "Inocybe" = "#dba4a4",
                               "Rhizopogon"="#111b77",
                               "Scleroderma"="#283dff", 
                               "Stephanospora"="#636bb7",
                               "Ramaria"="#bfc5ff",
                               "Tomentella"="#195637",
                               "Ceratobasidium"="#117744",
                               "Serendipita"="#60ffaf",
                               "Cenococcum"="#b7ffdb",
                               "Peziza"="#825121",
                               "Chromelosporium"="#ea7f17",
                               "Oidiodendron"="#fcb067",
                               "Chloridium"="#ffe8d3",
                               "Phialophora"="#d8d6d4",
                               "Capronia"="#82807f",
                               "Lecythophora" = "#3f3e3d",
                               "Glomus" = "#000000",
                               "Rhizophagus" = "#5b5b19", 
                               "Septoglomus" = "#fcfc00",
                               "Claroideoglomus"="#ffff9e",
                               "Paraglomus"="#521899", 
                               "Corymbiglomus"="#ae09ea",
                               "Diversispora"="#fa7efc",
                               "Gigaspora"="#ffb7ef",
                               "Racocetra"="#014443",
                               "Scutellospora"="#2b8c8a",
                               "Dentiscutata"="#14fffa",
                               "Acaulospora"="#a0fffc",
                               "Ambispora"="#e0f8fc",
                               "unidentified"="white"))

####
#tax_table(biom_ITS_ecm_gen_merg_ab) <- factor(tax_table(biom_ITS_ecm_gen_merg_ab), levels = c("OTU_5856","OTU_22","OTU_1129","OTU_575","OTU_1173","OTU_2360","OTU_428",
 #                                                                                             "OTU_4070","OTU_1533","OTU_1194","OTU_2637","OTU_5315","OTU_148","OTU_2028",
  #                                                                                            "OTU_1093","OTU_219","OTU_1233","OTU_2424","OTU_55","OTU_308","OTU_39",
   #                                                                                           "OTU_475","OTU_1822","OTU_531","OTU_2974","OTU_4177","OTU_373","OTU_2945","OTU_3995",
     #                                                                                         "OTU_125","OTU_53","OTU_2642","OTU_1215","OTU_24","OTU_3041"))

#tax_table(biom_ITS_ecm_gen_merg_ab)$Genus <- factor(tax_table(biom_ITS_ecm_gen_merg_ab)$Genus, 
           #                                         c("Amanita","Entoloma","Hygrocybe","Inocybe","Rhizopogon",
              #                                        "Scleroderma","Stephanospora","Ramaria","Tomentella","Ceratobasidium",
             #                                         "Serendipita","Cenococcum","Peziza","Chromelosporium","Oidiodendron","Chloridium",
               #                                       "Phialophora","Capronia","Lecythophora","Glomus","Rhizophagus","Septoglomus",
              #                                        "Claroideoglomus","Paraglomus","Corymbiglomus","Diversispora","Gigaspora",
                 #                                     "Racocetra","Scutellospora","Dentiscutata","Acaulospora","Ambispora","unidentified",
               #                                       "unidentified","unidentified"))

p1 + scale_fill_manual(values=colours[1:(36+1)])
plot_bar(ent10, "Genus", fill = "Genus", facet_grid = SeqTech ~ Enterotype)


merged_ITS
merged_ITS_gen = tax_glom(merged_ITS, "Genus")
merged_ITS_gen
merged_ITS_gen_merg = merge_samples(merged_ITS_gen, "Site1")
merged_ITS_gen_merg
sample_data(merged_ITS_gen_merg)
sample_data(merged_ITS_gen_merg)$Site1 <- levels(sample_data(merged_ITS_gen)$Site1)
merged_ITS_gen_merg_ab = transform_sample_counts(merged_ITS_gen_merg, function(x) 100 * x/sum(x))
merged_ITS_gen_merg_ab
tax_table(biom_ITS_ecm_gen_merg_ab)






################################################################################
################################################################################
### BETA DIVERSITY in Vegan ----------------------------------------------------------------------------------
sample_data(biom_ITS_hell)
sample_data(biom_LSU_hell)

biom_ITS_hell2 <- subset_samples(biom_ITS_hell, Site2%in%c("Beach_Intertidal1","Beach_Intertidal2","Cassava1","Cassava2",
                                                           "Cassuarina1","Grassland2","Mohogany1","Mohogany2",
                                                           "Mohogany3","Native_Forest1","Native_Forest3","Pinus_carribeana",
                                                           "Stream_Bank1","Stream_Bank2","Stream_Bank3","Sugarcane1","Sugarcane2"))

biom_ITS_hell2
sample_data(biom_ITS_hell2)

otu_ITS_hell2 <- as.data.frame(otu_table(biom_ITS_hell2))
taxa_ITS_hell2 <- as.data.frame(as.matrix(tax_table(biom_ITS_hell2)))
metadata_ITS_hell2 <- as.data.frame(as.matrix(sample_data(biom_ITS_hell2)[,2:9]))
identical(colnames(otu_ITS_hell2), rownames(metadata_ITS_hell2))



fungi_ITS_hl <- as.matrix(otu_table(biom_ITS_hell))
fungi_ITS_hl
fungi_LSU_hl <- as.matrix(otu_table(biom_LSU_hell))
fungi_LSU_hl
bacteria_hl <- as.matrix(otu_table(biom_16s_hell))
bacteria_hl

# import, filter reorder samples before analyses in vegan
metadata_ITS = read.csv("mapping_ITS.csv", header=T, row.names=1, sep=",")
metadata_LSU = read.csv("mapping_LSU.csv", header=T, row.names=1, sep=",")
metadata_16s = read.csv("mapping_16s.csv", header=T, row.names=1, sep=",")

metadata_ITS <- subset(metadata_ITS, Keep%in%c("keep"))
metadata_LSU <- subset(metadata_LSU, Keep%in%c("keep"))
metadata_16s <- subset(metadata_16s, Keep%in%c("keep"))

identical(colnames(fungi_ITS_hl), rownames(metadata_ITS))
identical(colnames(fungi_LSU_hl), rownames(metadata_LSU))
identical(colnames(bacteria_hl), rownames(metadata_16s))

# Non-metric Multidimensional Scaling 
otu_ITS_hell2

fungi_NMDS_ITS = metaMDS(t(otu_ITS_hell2), k=2, trymax=200, distance="bray", weakties = TRUE)
fungi_NMDS_ITS_2 = metaMDS(t(otu_ITS_hell2), k=2, trymax=200, distance="raup", weakties = TRUE)
fungi_NMDS_ITS_ja = metaMDS(t(otu_ITS_hell2), k=2, trymax=200, distance="jaccard", weakties = TRUE)

plot(fungi_NMDS_ITS_2, display = c("sites", "species"))
fungi_NMDS_ITS_2
stressplot(fungi_NMDS_ITS_2, main="Stressplot NMDS ITS")


fungi_NMDS_ITS = metaMDS(t(fungi_ITS_hl), k=2, trymax=200, distance="bray", weakties = TRUE)
fungi_NMDS_ITS
stressplot(fungi_NMDS_ITS, main="Stressplot NMDS ITS")

library("BiodiversityR")
plot(fungi_NMDS_ITS, display = c("sites", "species"))
fig <- ordiplot(fungi_NMDS_ITS, type = "none")
text(fig, "species", col="blue", cex=0.5)
text(fig, "sites", col="red", bg="red", cex=1)


fungi_NMDS_LSU = metaMDS(t(fungi_LSU_hl), k=2, trymax=200, distance="bray", weakties = TRUE)

fungi_NMDS_LSU
plot(fungi_NMDS_LSU, display = c("sites", "species"))
stressplot(fungi_NMDS_LSU, main="Stressplot NMDS LSU")

fungi_NMDS_16s = metaMDS(t(bacteria_hl), k=2, trymax=200, distance="bray", weakties = TRUE)
fungi_NMDS_16s
stressplot(fungi_NMDS_16s, main="Stressplot NMDS 16S")


par(mfrow = c(1, 3))
stressplot(fungi_NMDS_ITS, main="Stressplot NMDS ITS")
stressplot(fungi_NMDS_LSU, main="Stressplot NMDS LSU")
stressplot(fungi_NMDS_16s, main="Stressplot NMDS 16S")
dev.off()


#>>> procustes analysis ------------------------------------------------------------------------------------------------------------------------------------------
# http://www.flutterbys.com.au/stats/tut/tut15.1.html

library(vegan)

fungi_NMDS_LSU
fungi_NMDS_ITS

protest_NMDS <- protest(fungi_NMDS_LSU, fungi_NMDS_ITS, permutations = 9999)

protest_NMDS_ja <- protest(fungi_NMDS_ITS_2, fungi_NMDS_ITS, permutations = 999)
#protest_NMDS_16s <- protest(fungi_NMDS_16s, fungi_NMDS_ITS_2, permutations = 9999)

procrustes_NMDS <- procrustes(fungi_NMDS_LSU, fungi_NMDS_ITS_2)
plot(procustes_NMDS)
procustes_NMDS


plot(protest_NMDS, kind = "1", main = "Procrustes \nsuperimposition plot")
plot(protest_NMDS, kind = "2", main = "Procrustes \nrotation residuals plot")


plot(protest_NMDS)
protest_NMDS

plot(protest_NMDS_2,  kind = "2")
protest_NMDS_2

plot(protest_NMDS_ja)
protest_NMDS_ja

par(mfrow = c(1,2))
plot(protest_NMDS, kind = "1", main = "Procrustes \nsuperimposition plot")
plot(protest_NMDS, kind = "2", main = "Procrustes \nrotation residuals plot")
par(mfrow = c(1,1))






## perform CAP analysis --------------------------------------------
cap_its <- capscale(t(fungi_ITS_hl) ~ Site2, metadata_ITS, distance = "bray")
cap_its
anova(cap_its, by = "margin")

cap_lsu <- capscale(t(fungi_LSU_hl) ~ Site2, metadata_LSU, distance = "bray")
cap_lsu
anova(cap_lsu, by = "margin")

cap_bact <- capscale(t(bacteria_hl) ~ Site2, metadata_16s, distance = "bray")
cap_bact
anova(cap_bact, by = "margin")

RsquareAdj(cap_fungi)$adj.r.squared
RsquareAdj(cap_lsu)$adj.r.squared
RsquareAdj(cap_bact)$adj.r.squared



## perform RDA - Redundancy Analysis
rda_ITS <- rda(formula = t(fungi_ITS_hl) ~ Site2, data = metadata_ITS)
rda_ITS
anova(rda_ITS)








### MULTIVARIATE ANALYSIS OR VARIANCE - PERMANOVA


# Testing for beta dispersion - you need vegan objects
vegan::vegdist(t(fungi_ITS_hl), method="bray") -> dist_ITS_hl
vegan::vegdist(t(fungi_LSU_hl), method="bray") -> dist_LSU_hl
vegan::vegdist(t(bacteria_hl), method="bray") -> dist_16s_hl

permdisp_ITS <- betadisper(dist_ITS_hl, metadata_ITS$Site2)
permdisp_LSU <- betadisper(dist_LSU_hl, metadata_LSU$Site2)
permdisp_16s <- betadisper(dist_16s_hl, metadata_16s$Site2)

anova(permdisp_ITS)
anova(permdisp_LSU)
anova(permdisp_16s)

permutest(permdisp_ITS, permutations = 9999)
permutest(permdisp_LSU, permutations = 99999)
permutest(permdisp_16s, permutations = 9999)

par(mar=c(9, 4.1, 4.1, 2.1))
boxplot(permdisp_16s,las=3, cex=1, main="permdisp 16s, ANOVA p=0.399")
dev.off()

plot(permdisp_ITS)
plot(permdisp_LSU)
plot(permdisp_16s)

par(mfrow = c(1, 3))
plot(permdisp_ITS)
plot(permdisp_LSU)
plot(permdisp_16s)
dev.off()

par(mfrow = c(3,1), mar=c(8, 4.1, 4.1, 2.1))
boxplot(permdisp_ITS,las=3, cex=1, main="permdisp ITS, ANOVA p=0.078, perm.=9999")
boxplot(permdisp_LSU,las=3, cex=1, main="permdisp LSU, ANOVA p=0.053, perm.=9999")
boxplot(permdisp_16s,las=3, cex=1, main="permdisp 16s, ANOVA p=0.390, perm.=9999")
dev.off()




########################
# bacteria

vegan::vegdist(t(bact_Gen_hl), method="bray") -> dist_bact_Gen_hl

permdisp_bact1 <- betadisper(dist_bact_Gen_hl, metadata_16s$Plant)
permdisp_bact2 <- betadisper(dist_bact_Gen_hl, metadata_16s$Level)
permdisp_bact3 <- betadisper(dist_bact_Gen_hl, metadata_16s$Soil)

anova(permdisp_bact1)
anova(permdisp_bact2)
anova(permdisp_bact3)
permutest(permdisp_bact2, permutations = 999)









# Analysis of variance using distance matrices — for partitioning distance matrices 
# among sources of variation and fitting linear models (e.g., factors, polynomial regression)
# to distance matrices; uses a permutation test with pseudo-F ratios. 

adonis(t(fungi_ITS_hl) ~ Site2, data=metadata_ITS, method="bray", permutations=9999)
adonis(t(fungi_LSU_hl) ~ Site2, data=metadata_LSU, method="bray", permutations=9999)
adonis(t(bacteria_hl) ~ Site2, data=metadata_16s, method="bray", permutations=9999)
       
adonis(FormulaParsimoniousModel, data=metadata, permutations=9999, by="margin")

# Try adonis using type II test instead of tyype I
library("RVAideMemoire")

adonis.II(t(dist_ITS_hl) ~ Site2, data=metadata_ITS, permutations=9999, by="margin")
adonis.II(t(dist_LSU_hl) ~ Site2, data=metadata_LSU, permutations=9999, by="margin")
adonis.II(t(dist_LSU_hl) ~ Site2, data=metadata_LSU, permutations=9999, by="margin")


##more code to explore for adonis
#https://rpubs.com/michberr/randomforestmicrobe

adonis.site <- phyloseq_to_adonis(
  physeq = sites.scale, 
  dist = "bray", 
  formula = "Station"
)

#permutation ANOVA for alpha values and pairwise comparisons
obs_Level2_permanova<-perm.anova(bact~Level2, data=map_diversity16s, nperm=9999) #example with observed OTU richness values
obs_var1_pairwise<-pairwise.perm.t.test(map_diversity$obs, map_diversity$var1, nperm=9999, p.method="fdr")








################################################################################


### USING RAM PACKAGE
#https://rdrr.io/cran/RAM/
#https://cran.r-project.org/web/packages/RAM/index.html

# importing mappings files from txt
metadata_16s = read.table("mapping_16s.txt", header=TRUE, row.names=1, sep="\t")
metadata_16s

metadata_ITS = read.table("mapping_ITS.txt", header=TRUE, row.names=1, sep="\t")
metadata_ITS


# filtering to equal reps x group
# metadata_ITS <- metadata_ITS[-c(10, 16, 17, 33, 37, 43, 44),]
# metadata_16s <- metadata_16s[-c(13, 17, 18, 35, 45, 46),]

# importing otu_tabes form biom.txt files with taxonomy
# Remember to remove the first line starts with # and 
# convert the taxonomy column into a character
otutab_16s <- read.table("otu_table_16s.txt", header=T, row.names=1, sep="\t")
head(otutab_16s)
str(otutab_16s)
otutab_16s$taxonomy <- as.character(otutab_16s$taxonomy)

otutab_ITS <- read.table("otu_table_ITS_.txt", header=T, row.names=1, sep="\t")
head(otutab_ITS)
str(otutab_ITS)
otutab_ITS$taxonomy <- as.character(otutab_ITS$taxonomy)


# changing order of levels. It should happen you need
# to change variable level order for a better display
sample_data(biom_R1)$Time<-factor(sample_data(biom_R1)$Time, levels=c("day1","day5","day20"))
levels(sample_data(biom_R1)$Time)
sample_data(biom_R1)


# plots the top 10 OTUs (by default) at five ranks
data=list(bacteria=otutab_16s)
data

group.top.number(data=list(bacteria=otutab_16s), ranks=c("p","f","g"), top=10, drop.unclassified=TRUE, 
                 main="Relative Abundance of the Top 10 Taxa", cex.x=10)

group.top.percent(data=list(bacteria=otutab_16s), top=10, drop.unclassified=TRUE)


group.top.number(data=list(fungi=otutab_ITS), ranks=c("p","f","g"), top=10, drop.unclassified=TRUE, 
                 main="Relative Abundance of the Top 10 Taxa", cex.x=10)

group.top.percent(data=list(fungi=otutab_ITS), top=10, drop.unclassified=TRUE)

# stacked barplot of abundances
factor.abundance(data=list(bacteria=otutab_16s), rank="genus", meta=metadata_16s,
                 meta.factor=c("Description"), top=20,
                 drop.unclassified=FALSE)

factor.abundance(data=list(fungi=otutab_ITS), rank="genus", meta=metadata_ITS,
                 meta.factor=c("Description"), top=20,
                 drop.unclassified=FALSE)

# PERMANOVA using adonis
# test OTUs
data <- list(fungi=otutab_ITS)
assist.ado(data=data, is.OTU=TRUE, meta=metadata_ITS, ranks=NULL,
           data.trans="log", dist=NULL)

# test taxa at different ranks
ranks <- c("p", "c", "o", "f", "g")
ado <- assist.ado(data=data, is.OTU=TRUE, 
                  meta=meta, ranks=ranks,
                  data.trans="log", dist="bray" )
# test genera
g1 <- tax.abund(otu1=ITS1, rank="g", drop.unclassified=TRUE)
data <- list(g1=g1)
assist.ado(data=data, is.OTU=FALSE, 
           meta=meta, ranks=NULL,
           data.trans="log", dist="bray" )


## Ordinations
data(ITS1, meta)
its1 <- filter.OTU(data=list(fungi=otutab_ITS), percent=0.001)[[1]]
factors=c("Soil", "Genotype")
factors

# plot sites
ord1 <- OTU.ord(its1, meta=metadata_ITS)
ord1 <- OTU.ord(its1,  meta=metadata_gen_ITS, data.trans="total",
                factors=factors, mode="rda", biplot.sig=0.1,
                taxa=20, biplot.scale=1.5, cex.point=5,
                plot.species=FALSE, rank="f", plot.scaling=3,
                group=c(Soil="Soil", Genotype="Genotype"))


# plot taxa
ord2 <- OTU.ord(its1,  meta=metadata_gen_ITS, data.trans="total", 
                plot.scaling=-1,
                factors=factors, mode="cca", biplot.sig=0.1,
                taxa=20, biplot.scale=3, cex.point=3,
                plot.species=TRUE, rank="g")


# isa analysis 
group.indicators(data=list(bacteria=otutab_16s), is.OTU=TRUE, thresholds = c(A = 0.85,
      B = 0.8,stat = 0.8,p.value = 0.05), metadata_16s,factor = c(Description="Genotype"),rank="f")


group.indicators(data=list(Bacteria=otutab_16s), is.OTU=TRUE, thresholds = c(A = 0.85, 
      B = 0.8,stat = 0.8,p.value = 0.05), metadata_16s,factor = c(Description="Description"),rank="g")


# plot heatmap
group.heatmap.simple(otutab_ITS, is.OTU=TRUE, meta=metadata_ITS,
                     row.factor=c(Soil="Soil"), dendro="row",
                     rank="g", top=10, drop.unclassified=TRUE,
                     leg.x=-0.06)















##############################################################################################################
# useful code for shapes and lines - ggplot2 shapes
d=data.frame(p=c(0:25,32:127))
ggplot() +
  scale_y_continuous(name="") +
  scale_x_continuous(name="") +
  scale_shape_identity() +
  geom_point(data=d, mapping=aes(x=p%%16, y=p%/%16, shape=p), size=5, fill="red") +
  geom_text(data=d, mapping=aes(x=p%%16, y=p%/%16+0.25, label=p), size=3)






# plotting library 

sums_biom_ITS <- data.frame(colSums(otu_table(biom_ITS_cv)))
colnames(sums_biom_ITS) <- "Sample_TotalSeqs"
sums_biom_ITS$sample <- row.names(sums_biom_ITS)
sums_biom_ITS$Description <- sample_data(biom_ITS_cv)$Description
sums_biom_ITS$Site2 <- sample_data(biom_ITS_cv)$Site2 # if you want to add another variable to use further
sums_biom_ITS <- arrange(sums_biom_ITS, Sample_TotalSeqs)
#sums_biom_ITS <- arrange(sums_biom_ITS, Description) # to order according another variable in the mapping
sums_biom_ITS

ggplot(sums_biom_ITS, aes(x=reorder(Site2, Sample_TotalSeqs), y = Sample_TotalSeqs)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=16, outlier.size=4) +
  geom_point() + geom_text_repel(aes(label=sample), size = 3) + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  labs(x = "Sample Replicate", y = "Read Number") +
  ggtitle("Sample distribution ITS")



######

## New script for barplots and relative abundances


biom_ITS_symb = import_biom("mycorrhizal_genus.biom")
map_ITS = import_qiime_sample_data("mapping_ITS.txt")
sample_data(biom_ITS_symb) <- map_ITS
colnames(tax_table(biom_ITS_symb)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
otus_rep_ITS_symb <- readDNAStringSet("otus_ITS.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
biom_ITS_symb <- merge_phyloseq(biom_ITS_symb, otus_rep_ITS_symb)
biom_ITS_symb
sample_data(biom_ITS_symb)
tax_table(biom_ITS_symb)
refseq(biom_ITS_symb)


biom_ITS_symb_gen = tax_glom(biom_ITS_symb, "Genus")
biom_ITS_symb_gen
biom_ITS_symb_gen_merg = merge_samples(biom_ITS_symb_gen, "Site1")
biom_ITS_symb_gen_merg
sample_data(biom_ITS_symb_gen_merg)
otu_table(biom_ITS_symb_gen_merg)

sample_data(biom_ITS_symb_gen_merg)$Site1 <- levels(sample_data(biom_ITS_symb_gen)$Site1)
biom_ITS_symb_gen_merg_ab = transform_sample_counts(biom_ITS_symb_gen_merg, function(x) 100 * x/sum(x))
biom_ITS_symb_gen_merg_ab
otu_table(biom_ITS_symb_gen_merg_ab)
tax_table(biom_ITS_symb_gen_merg_ab)
sum(otu_table(biom_ITS_symb_gen_merg_ab))

source("../plot_ordered_bar_byGian.R")
p1 <- plot_ordered_bar(biom_ITS_symb_gen_merg_ab, x = "Site1", fill="Genus", leg_size = 0.4, title="(A)")
p1

myco_bar = p1 + theme(legend.text = element_text(face="italic")) +
  scale_fill_manual(values = c("Acaulospora" = "#560d0d",
                               "Amanita" = "#a35151",
                               "Ambispora" = "#dba4a4", 
                               "Cenococcum" = "#cc1c1c",
                               "Ceratobasidium"="#111b77",
                               "Chloridium"="#283dff", 
                               "Chromelosporium"="#636bb7",
                               "Claroideoglomus"="#bfc5ff",
                               "Corymbiglomus"="#014443",
                               "Dentiscutata"="#195637",
                               "Diversispora"="#117744",
                               "Entoloma"="#60ffaf",
                               "Gigaspora"="#b7ffdb",
                               "Glomus"="#825121",
                               "Inocybe"="#ea7f17",
                               "Oidiodendron"="#fcb067",
                               "Paraglomus"="#ffe8d3",
                               "Peziza"="#d8d6d4",
                               "Racocetra"="#82807f",
                               "Ramaria" = "#3f3e3d",
                               "Rhizophagus" = "#000000",
                               "Rhizopogon" = "#5b5b19", 
                               "Scleroderma" = "#fcfc00",
                               "Scutellospora"="#ffff9e",
                               "Septoglomus"="#ffb7ef", 
                               "Serendipita"="#fa7efc",
                               "Stephanospora"="#ae09ea",
                               "Tomentella"="#521899",
                               "Wilcoxina"="#1e0047",
                               "unidentified1"="white",
                               "unidentified2"="white",
                               "unidentified2"="white")) + 
  theme(panel.background = element_blank()) + 
  guides(fill=guide_legend(ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  theme(legend.position="right") +
  guides(fill=guide_legend(reverse = TRUE, ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  #facet_wrap(~Phylum, ncol=4) +
  #facet_grid(~Phylum, nrow=2) +
  ylim(0, 100) +
  labs(title="Mycorrhizal", x="", y="Relative abundance") +
  theme(panel.background = element_blank()) + 
  theme(axis.text.x = element_text(vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(hjust=0.5, size=10)) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  theme(axis.line = element_line(colour = "black",size = 1, linetype = "solid")) +
  theme(legend.position="right")


#theme(legend.position="none")
  
  #theme(legend.text=element_text(size=10), legend.key.size = unit(1, "mm"))
  #theme(legend.position = "bottom") +
  #theme(legend.direction = "horizontal")

## now endophytic fungi

biom_ITS_endo = import_biom("endophyte_genus.biom")
map_ITS = import_qiime_sample_data("mapping_ITS.txt")
sample_data(biom_ITS_endo) <- map_ITS
colnames(tax_table(biom_ITS_endo)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
otus_rep_ITS_endo <- readDNAStringSet("otus_ITS.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
biom_ITS_endo <- merge_phyloseq(biom_ITS_endo, otus_rep_ITS_endo)
biom_ITS_endo
sample_data(biom_ITS_endo)
tax_table(biom_ITS_endo)
refseq(biom_ITS_endo)

biom_ITS_endo_gen = tax_glom(biom_ITS_endo, "Genus")
biom_ITS_endo_gen_merg = merge_samples(biom_ITS_endo_gen, "Site1")
sample_data(biom_ITS_endo_gen_merg)$Site1 <- levels(sample_data(biom_ITS_endo_gen)$Site1)
biom_ITS_endo_gen_merg_ab = transform_sample_counts(biom_ITS_endo_gen_merg, function(x) 100 * x/sum(x))

p2 <- plot_ordered_bar(biom_ITS_endo_gen_merg_ab, x = "Site1", fill="Genus", leg_size = 0.4, title="(B)")
p2

endo_bar = p2 + theme(legend.text = element_text(face="italic")) + 
  scale_fill_manual(values = c("Annulohypoxylon" = "#560d0d",
                               "Anthostomella" = "#a35151",
                               "Camarosporium" = "#dba4a4", 
                               "Capronia" = "#cc1c1c",
                               "Chalara"="#111b77",
                               "Cladosporium"="#283dff", 
                               "Colletotrichum"="#636bb7",
                               "Cytospora"="#bfc5ff",
                               "Diaporthe"="#014443",
                               "Kabatiella"="#195637",
                               "Lecythophora"="#117744",
                               "Microdochium"="#60ffaf",
                               "Paecilomyces"="#b7ffdb",
                               "Periconia"="#825121",
                               "Phialophora"="#ea7f17",
                               "Phomopsis"="#fcb067",
                               "Torula"="#ffe8d3",
                               "Trichoderma"="#d8d6d4",
                               "Williopsis"="#82807f",
                               "Xylaria" = "#3f3e3d")) + 
  theme(panel.background = element_blank()) + 
  guides(fill=guide_legend(ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  theme(legend.position="right") + guides(fill=guide_legend(reverse = TRUE, ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  #facet_wrap(~Phylum, ncol=4) +
  #facet_grid(~Phylum, nrow=2) +
  ylim(0, 100) +
  labs(title="Endophytic", x="Habitat", y="") +
  theme(panel.background = element_blank()) + 
  theme(axis.text.x = element_text(vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(hjust=0.5, size=10)) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  theme(axis.line = element_line(colour = "black",size = 1, linetype = "solid")) +
  theme(legend.position="right")

## now pathogens

biom_ITS_pat = import_biom("pathogen_genus.biom")
map_ITS = import_qiime_sample_data("mapping_ITS.txt")
sample_data(biom_ITS_pat) <- map_ITS
colnames(tax_table(biom_ITS_pat)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
otus_rep_ITS_pat <- readDNAStringSet("otus_ITS.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
biom_ITS_pat <- merge_phyloseq(biom_ITS_pat, otus_rep_ITS_pat)
biom_ITS_pat
sample_data(biom_ITS_pat)
tax_table(biom_ITS_pat)
refseq(biom_ITS_pat)

biom_ITS_pat_gen = tax_glom(biom_ITS_pat, "Genus")
biom_ITS_pat_gen_merg = merge_samples(biom_ITS_pat_gen, "Site1")
biom_ITS_pat_gen_merg
tax_table(biom_ITS_pat_gen_merg)
otu_table(biom_ITS_pat_gen_merg)

sample_data(biom_ITS_pat_gen_merg)$Site1 <- levels(sample_data(biom_ITS_pat_gen)$Site1)
biom_ITS_pat_gen_merg_ab = transform_sample_counts(biom_ITS_pat_gen_merg, function(x) 100 * x/sum(x))


badTaxa = c("OTU_1494","OTU_2750","OTU_2559","OTU_3578","OTU_3405","OTU_3221","OTU_5205","OTU_3915","OTU_2624","OTU_3451"
            ,"OTU_3118","OTU_6904","OTU_4023","OTU_5200","OTU_5862","OTU_4387","OTU_5089","OTU_4523","OTU_4615","OTU_4945")

pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

biom_ITS_pat_gen_merg_ab_prun = pop_taxa(biom_ITS_pat_gen_merg_ab, badTaxa)
biom_ITS_pat_gen_merg_ab_prun

#top26 <- names(sort(taxa_sums(biom_ITS_pat_gen_merg_ab), TRUE)[1:26]) # sampling first 100 taxa
#biom_ITS_pat_gen_merg_ab_26 <- prune_taxa(top26, biom_ITS_pat_gen_merg_ab)

p3 <- plot_ordered_bar(biom_ITS_pat_gen_merg_ab_prun, x = "Site1", fill="Genus", leg_size = 0.4, title="(C)")
p3

pat_bar = p3 + theme(legend.text = element_text(face="italic")) +
  scale_fill_manual(values = c("Aspergillus" = "#560d0d",
                               "Calonectria" = "#a35151",
                               "Ceratocystis" = "#dba4a4", 
                               "Clonostachys" = "#cc1c1c",
                               "Collophora"="#111b77",
                               "Coniothyrium"="#283dff", 
                               "Curvularia"="#636bb7",
                               "Cylindrocarpon"="#bfc5ff",
                               "Cylindrocladiella"="#014443",
                               "Endomelanconiopsis"="#195637",
                               "Fusarium"="#117744",
                               "Ganoderma"="#60ffaf",
                               "Gibberella"="#b7ffdb",
                               "Gibellulopsis"="#825121",
                               "Glomerella"="#ea7f17",
                               "Haematonectria"="#fcb067",
                               "Lasiodiplodia"="#ffe8d3",
                               "Magnaporthe"="#d8d6d4",
                               "Monographella"="#82807f",
                               "Mycoleptodiscus" = "#3f3e3d",
                               "Pestalotiopsis" = "#000000",
                               "Phaeoacremonium" = "#5b5b19", 
                               "Phoma" = "#c1c15d",
                               "Pyricularia"="#fcfc00",
                               "Rhizophydium"="#ffff9e", 
                               "Sclerotinia"="#ffb7ef",
                               "Spizellomyces"="#fa7efc",
                               "Thanatephorus"="#ae09ea",
                               "Trichothecium"="#521899",
                               "Volutella"="#1e0047")) + 
  theme(panel.background = element_blank()) + 
  guides(fill=guide_legend(ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  theme(legend.position="right") +
  guides(fill=guide_legend(ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  theme(legend.position="right") + guides(fill=guide_legend(reverse = TRUE, ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  #facet_wrap(~Phylum, ncol=4) +
  #facet_grid(~Phylum, nrow=2) +
  ylim(0, 100) +
  labs(title="Pathogenic", x="", y="") +
  theme(panel.background = element_blank()) + 
  theme(axis.text.x = element_text(vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(hjust=0.5, size=10)) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  theme(axis.line = element_line(colour = "black",size = 1, linetype = "solid")) +
  theme(legend.position="right")


library("ggpubr")
ggarrange(myco_bar,endo_bar,pat_bar,
          labels = c("A", "B","C"),
          widths = c(1,1,1.05),
          align = "h", ncol = 3, nrow = 1,
          legend="right",
          common.legend = FALSE)



biom_16S_gen = tax_glom(biom_16s_cv, "Genus")
biom_16S_gen
biom_16S_gen_merg = merge_samples(biom_16S_gen, "Site1")
biom_16S_gen_merg
sample_data(biom_16S_gen_merg)
otu_table(biom_16S_gen_merg)

sample_data(biom_16S_gen_merg)$Site1 <- levels(sample_data(biom_16S_gen)$Site1)
biom_16S_gen_merg_ab = transform_sample_counts(biom_16S_gen_merg, function(x) 100 * x/sum(x))
biom_16S_gen_merg_ab
otu_table(biom_16S_gen_merg_ab)
tax_table(biom_16S_gen_merg_ab)
sum(otu_table(biom_16S_gen_merg_ab))


badTaxa = c("OTU_30")

pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

biom_16s_gen_merg_ab_prun = pop_taxa(biom_16S_gen_merg_ab, badTaxa)
biom_16s_gen_merg_ab_prun


top_30_bact <- names(sort(taxa_sums(biom_16s_gen_merg_ab_prun), TRUE)[1:30]) # first 30 top OTUs 
top_30_bact
biom_16s_top30 <- prune_taxa(top_30_bact, biom_16s_gen_merg_ab_prun)
biom_16s_top30
tax_table(biom_16s_top30)

# get sequnces for BLAST against GenBank
write.csv(refseq(biom_16S_gen) , file = "refseq_16s_gen.csv")
# manual correction of the tax_table - OPTIONAL!
write.csv(tax_table(biom_16s_top30) , file = "tax_16s_top30.csv")
tax_biom_16s_top30 <- read.csv("tax_16s_top30.csv", header=T, row.names =1)
tax_table(biom_16s_top30) <- tax_table(as.matrix(tax_biom_16s_top30))
tax_table(biom_16s_top30)


source("../plot_ordered_bar_byGian.R")

p_bact <- plot_ordered_bar(biom_16s_top30, x = "Site1", fill="Genus", leg_size = 0.4, title="(D)")
p_bact

p_bact + theme(legend.text = element_text(face="italic")) +
  scale_fill_manual(values = c("Actinoallomurus" = "#560d0d",
                               "Aeribacillus" = "#a35151",
                               "Bacillus" = "#dba4a4", 
                               "Bradyrhizobium" = "#cc1c1c",
                               "Burkholderia"="#111b77",
                               "Candidatus Koribacter"="#283dff", 
                               "Candidatus Nitrososphaera"="#636bb7",
                               "Candidatus Solibacter"="#bfc5ff",
                               "Candidatus Xiphinematobacter"="#014443",
                               "Flavobacterium"="#195637",
                               "Gemmata"="#117744",
                               "Geobacillus"="#60ffaf",
                               "Hyphomicrobium"="#b7ffdb",
                               "Kaistobacter"="#825121",
                               "Ktedonobacter"="#ea7f17",
                               "Mesorhizobium"="#fcb067",
                               "Mycobacterium"="#ffe8d3",
                               "Nitrospira"="#d8d6d4",
                               "Nocardioides"="#82807f",
                               "Pedomicrobium" = "#3f3e3d",
                               "Phenylobacterium" = "#000000",
                               "Pirellula" = "#5b5b19", 
                               "Planctomyces" = "#c1c15d",
                               "Pseudonocardia"="#fcfc00",
                               "Rhizobium"="#ffff9e", 
                               "Rhodoplanes"="#ffb7ef",
                               "[Spartobacteria] DA101"="#fa7efc",
                               "Sphingomonas"="#ae09ea",
                               "Steroidobacter"="#521899",
                               "Streptomyces"="#1e0047")) + 
  theme(panel.background = element_blank()) + 
  guides(fill=guide_legend(ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  theme(legend.position="right") +
  guides(fill=guide_legend(ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  theme(legend.position="right") + guides(fill=guide_legend(reverse = TRUE, ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  #facet_wrap(~Phylum, ncol=4) +
  #facet_grid(~Phylum, nrow=2) +
  ylim(0, 100) +
  labs(title="Bacteria and Archaea", x="Habitat", y="Relative abundance") +
  theme(panel.background = element_blank()) + 
  theme(axis.text.x = element_text(vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(hjust=0.5, size=10)) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  theme(axis.line = element_line(colour = "black",size = 1, linetype = "solid")) +
  theme(legend.position="right")





#################################################################################################
#################################################################################################

## Calculate Relative Abundacnes

biom_ITS_cv
biom_ITS_phy = tax_glom(biom_ITS_cv, "Phylum")
otu_ITS_phy = taxa_sums(biom_ITS_phy)/sum(taxa_sums(biom_ITS_phy))*100
tax_ITS_phy <- data.frame(tax_table(biom_ITS_phy))
tax_ITS_phy <- tax_ITS_phy[c(1:2)]
tax_ITS_phy$abundance <- as.vector(otu_ITS_phy)
tax_ITS_phy <- tax_ITS_phy[order(tax_ITS_phy$abundance, decreasing = TRUE),] 
tax_ITS_phy

biom_LSU_cv
biom_LSU_phy = tax_glom(biom_LSU_cv, "Phylum")
otu_LSU_phy = taxa_sums(biom_LSU_phy)/sum(taxa_sums(biom_LSU_phy))*100
tax_LSU_phy <- data.frame(tax_table(biom_LSU_phy))
tax_LSU_phy <- tax_LSU_phy[c(1:2)]
tax_LSU_phy$abundance <- as.vector(otu_LSU_phy)
tax_LSU_phy <- tax_LSU_phy[order(tax_LSU_phy$abundance, decreasing = TRUE),] 
tax_LSU_phy

biom_16s_cv
biom_16s_phy = tax_glom(biom_16s_cv, "Phylum")
otu_16s_phy = taxa_sums(biom_16s_phy)/sum(taxa_sums(biom_16s_phy))*100
tax_16s_phy <- data.frame(tax_table(biom_16s_phy))
tax_16s_phy <- tax_16s_phy[c(1:2)]
tax_16s_phy$abundance <- as.vector(otu_16s_phy)
tax_16s_phy <- tax_16s_phy[order(tax_16s_phy$abundance, decreasing = TRUE),] 
tax_16s_phy




# different guilds

biom_ITS_symb
tax_table(biom_ITS_symb)
otu_table(biom_ITS_symb)

otu_ITS_symb_gen = taxa_sums(biom_ITS_symb)/sum(taxa_sums(biom_ITS_symb))*100
tax_ITS_symb_gen <- data.frame(tax_table(biom_ITS_symb))
tax_ITS_symb_gen <- tax_ITS_symb_gen[c(1:6)]
tax_ITS_symb_gen$abundance <- as.vector(otu_ITS_symb_gen)
tax_ITS_symb_gen <- tax_ITS_symb_gen[order(tax_ITS_symb_gen$abundance, decreasing = TRUE),] 
tax_ITS_symb_gen

otu_ITS_endo_gen = taxa_sums(biom_ITS_endo)/sum(taxa_sums(biom_ITS_endo))*100
tax_ITS_endo_gen <- data.frame(tax_table(biom_ITS_endo))
tax_ITS_endo_gen <- tax_ITS_endo_gen[c(1:6)]
tax_ITS_endo_gen$abundance <- as.vector(otu_ITS_endo_gen)
tax_ITS_endo_gen <- tax_ITS_endo_gen[order(tax_ITS_endo_gen$abundance, decreasing = TRUE),] 
tax_ITS_endo_gen

biom_ITS_pat_gen = tax_glom(biom_ITS_pat, "Genus")
otu_ITS_pat_gen = taxa_sums(biom_ITS_pat_gen)/sum(taxa_sums(biom_ITS_pat_gen))*100
tax_ITS_pat_gen <- data.frame(tax_table(biom_ITS_pat_gen))
tax_ITS_pat_gen <- tax_ITS_pat_gen[c(1:6)]
tax_ITS_pat_gen$abundance <- as.vector(otu_ITS_pat_gen)
tax_ITS_pat_gen <- tax_ITS_pat_gen[order(tax_ITS_pat_gen$abundance, decreasing = TRUE),] 
tax_ITS_pat_gen




## alpha diversity

plot_alpha_its <- plot_richness(biom_ITS_cv, x = "Site1", color = "Site1", measures = c("Observed","Simpson"))
plot_alpha_lsu <- plot_richness(biom_LSU_cv, x = "Site1", color = "Site1", measures = c("Observed","Shannon"))
plot_alpha_16s <- plot_richness(biom_16s_cv, x = "Site1", color = "Site1", measures = c("Observed","Shannon"))

rich_its = plot_alpha_its + geom_boxplot(outlier.colour="black", outlier.shape = 8) +
scale_colour_manual(values=palette_CB9) +
theme(legend.position="none")

rich_lsu = plot_alpha_lsu + geom_boxplot(outlier.colour="black", outlier.shape = 8) +
scale_colour_manual(values=palette_CB9) +
theme(legend.position="none")

rich_16s = plot_alpha_16s + geom_boxplot(outlier.colour="black", outlier.shape = 8) +
scale_colour_manual(values=palette_CB9) +
theme(legend.position="none")
rich_16s

# multiplot
source("../multiplot.R")
multiplot(rich_its, rich_lsu, rich_16s, cols=3)
dev.off()


## plotting heatmaps
library("superheat")
library("vegan")
data(dune)
dune

superheat(t(dune), left.label.size = 0.3, bottom.label.size = 0.1,
          pretty.order.cols = TRUE)

superheat(t(dune), left.label.size = 0.3, bottom.label.size = 0.1,
          pretty.order.rows = TRUE, 
          heat.pal = c("#b35806", "white", "#542788"))


superheat(t(dune), left.label.size = 0.3, bottom.label.size = 0.1,
          heat.pal = c("#b35806", "white", "#542788"), row.dendrogram = TRUE)
          
          
          
         























