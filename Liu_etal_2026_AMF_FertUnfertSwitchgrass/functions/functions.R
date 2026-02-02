#************************************************************************-------
# Manuscript Title: Site-specific factors rather than nitrogen impact arbuscular 
#                   mycorrhizal fungi diversity in bioenergy switchgrass monocultures
# Authors:          Shuang Liu, Gian Maria Niccolò Benucci, Alden Dirks, 
#                   Lukas Bell-Dereske, Sarah Evans, Gregory Bonito
# Code Developer:   Gian MN Benucci 2025
# Citation:         ...
#                   
# DOI               ...
# PMID:             ...
# **********************************************************************--------

# R options --------------------------------------------------------------------
options(scipen = 9999, pillar.sigfig = 6, digits = 6, max.print = 10000000)
#rm(list = ls())

# Check the lib paths ----------------------------------------------------------
.libPaths()

# ********** FUNCTIONS ********** ----------------------------------------------

# Stadardize NCBI blasTAX taxonomy ---------------------------------------------
extract_blasTAX <- function(tax_path,
                            namemap_path){
  
  taxonomy <- read.delim(tax_path, header = TRUE, strip.white = TRUE, sep="\t") 
  name_mapping <- read.delim(namemap_path, header = FALSE, strip.white = TRUE, sep="\t") 
  
  if (length(setdiff(name_mapping$V2, taxonomy$OTU_ID)) > 0) {
    cli::cli_alert_warning(
      paste("No significant similarity found for:", 
            paste(setdiff(name_mapping$V2, taxonomy$OTU_ID), collapse = ", "))
    )
  }
  
  taxonomy_clean <-
    taxonomy %>% 
    rename(Query = 1, Query_Score = 2) %>% 
    mutate_all(~ str_trim(., side = "both")) %>% # remove trailing spaces in both sides
    as.matrix() %>% 
    as.data.frame() %>% 
    dplyr::arrange(Query)
  
  name_mapping_clean <-
    name_mapping %>% 
    rename(Zotu = V1, Query = V2) %>% 
    mutate(Zotu = gsub(">", "", Zotu)) %>% 
    mutate_all(~ str_trim(., side = "both")) %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    dplyr::arrange(Query)
  
  blastax <-
    full_join(name_mapping_clean, taxonomy_clean, by="Query")
  
  return(blastax)
  
}


extract_blasTAX(
  tax_path = file.path(
    data_path,
    "datasets/taxonomy_blast_99.txt"),
  namemap_path = file.path(
    data_path,
    "datasets/name_mapping_99.txt"
  )
)



# FInalize and modify taxonomy table -------------------------------------------

# Replace blanks ---------------------------------------------------------------
blank2na = function(x, na.strings=c('','.','NA','na','N/A','n/a','NaN','nan')) {
  if (is.factor(x)) {
    lab = attr(x, 'label', exact = T)
    labs1 <- attr(x, 'labels', exact = T)
    labs2 <- attr(x, 'value.labels', exact = T)
    # trimws will convert factor to character
    x = trimws(x,'both')
    if (! is.null(lab)) lab = trimws(lab,'both')
    if (! is.null(labs1)) labs1 = trimws(labs1,'both')
    if (! is.null(labs2)) labs2 = trimws(labs2,'both')
    if (!is.null(na.strings)) {
      # convert to NA
      x[x %in% na.strings] = NA
      # also remember to remove na.strings from value labels 
      labs1 = labs1[! labs1 %in% na.strings]
      labs2 = labs2[! labs2 %in% na.strings]
    }
    # the levels will be reset here
    x = factor(x)
    if (! is.null(lab)) attr(x, 'label') <- lab
    if (! is.null(labs1)) attr(x, 'labels') <- labs1
    if (! is.null(labs2)) attr(x, 'value.labels') <- labs2
  } else if (is.character(x)) {
    lab = attr(x, 'label', exact = T)
    labs1 <- attr(x, 'labels', exact = T)
    labs2 <- attr(x, 'value.labels', exact = T)
    # trimws will convert factor to character
    x = trimws(x,'both')
    if (! is.null(lab)) lab = trimws(lab,'both')
    if (! is.null(labs1)) labs1 = trimws(labs1,'both')
    if (! is.null(labs2)) labs2 = trimws(labs2,'both')
    if (!is.null(na.strings)) {
      # convert to NA
      x[x %in% na.strings] = NA
      # also remember to remove na.strings from value labels 
      labs1 = labs1[! labs1 %in% na.strings]
      labs2 = labs2[! labs2 %in% na.strings]
    }
    if (! is.null(lab)) attr(x, 'label') <- lab
    if (! is.null(labs1)) attr(x, 'labels') <- labs1
    if (! is.null(labs2)) attr(x, 'value.labels') <- labs2
  } else {
    x = x
  }
  return(x)
}


# Finalize ---------------------------------------------------------------------
FinalizeTaxonomy <- function(taxonomy){
  taxonomy$Species <- 
    gsub(" sp ", "", taxonomy$Species)
  taxonomy[] = lapply(taxonomy, blank2na, na.strings=c('','NA','na','N/A','n/a','NaN','nan'))
  lastValue <- function(x) tail(x[!is.na(x)], 1)
  taxonomy$Genus <- as.character(taxonomy$Genus)
  taxonomy[which(is.na(taxonomy$Genus) == FALSE),]$Genus <-
    paste(taxonomy$Genus[is.na(taxonomy$Genus) == FALSE], "sp.", sep = " ")
  last_taxons<- apply(taxonomy[,c(2:8)], 1, lastValue)
  taxonomy$BestMatch <- last_taxons
  taxonomy$BestMatch <-
    gsub("_", " ", taxonomy$BestMatch)
  taxonomy$Taxonomy <-
    paste(taxonomy$Zotu, taxonomy$BestMatch, sep = "-")
  taxonomy$BestMatch <- 
    gsub(" sp.", "", taxonomy$BestMatch)
  taxonomy$Genus <- 
    gsub(" sp.", "", taxonomy$Genus)
  return(taxonomy)
}


FinalizeTaxonomy(taxonomy_99)













