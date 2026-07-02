# >>> EXTRACT LAS TAXONOMIC LEVEL ------------------------------------------------------------------
# thanks to https://rdrr.io/github/jerryzhujian9/ezR/src/R/basic.R

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

# In the tax_table add a column naming the highest resolution taxonomy 
# achieved for each OTU, remove _ and add sp. to genera
ReformatTaxonomy <- function(dataframe){
  taxa_table = as(tax_table(dataframe), "matrix")
  taxa_table = as.data.frame(taxa_table)
  # remember to do run this function only once on your dataframe
  taxa_table$Genus <- as.character(taxa_table$Genus)
  taxa_table[taxa_table=="Unclassified"]<- NA
  taxa_table[taxa_table=="Unidentified"]<- NA
  taxa_table[taxa_table==""]<- NA
  #to remove species that aready have sp
  taxa_table$Species <- 
    gsub(" sp ", "", taxa_table$Species)
  # add sp. to species with just the Genus
  taxa_table[which(is.na(taxa_table$Genus) == FALSE),]$Genus <-
    paste(taxa_table$Genus[is.na(taxa_table$Genus) == FALSE], "sp.", sep = " ")
    taxa_table <- taxa_table[c(8,1,2,3,4,5,6,7,9,10)]
  taxa_table[] = lapply(taxa_table, blank2na, na.strings=c('','NA','na','N/A','n/a','NaN','nan'))
  lastValue <- function(x) tail(x[!is.na(x)], 1)
  last_taxons<- apply(taxa_table[,1:8], 1, lastValue)
  taxa_table$BestMatch <- last_taxons
  taxa_table[, "BestMatch"] <-
    gsub("_", " ", taxa_table[, "BestMatch"])
  taxa_table$Taxonomy <-
    paste(taxa_table$OTU_ID, taxa_table$BestMatch, sep = "-")
  taxa_table[, "Genus"] <- gsub(" sp.", "", taxa_table[, "Genus"])
  tax_table(dataframe) <- tax_table(as.matrix(taxa_table))
  return(dataframe)
}


# # In the tax_table add a column naming the highest resolution taxonomy 
# # achieved for each OTU, remove _ and add sp. to genera
# ReformatTaxonomy <- function(dataframe, taxa){
#   taxa_table <- as.data.frame(as.matrix(tax_table(dataframe)))
#   # remember to do run this function only once on your dataframe
#   taxa_table$Genus <- as.character(taxa_table$Genus)
#   taxa_table[taxa_table=="Unclassified"]<- NA
#   taxa_table[taxa_table==""]<- NA
#   taxa_table[which(is.na(taxa_table$Genus) == FALSE), ]$Genus <- paste(
#     taxa_table$Genus[is.na(taxa_table$Genus)==FALSE], "sp.", sep = " ")
#   taxa_table$OTU <- row.names(taxa_table)
#   if (taxa == "ITS"){
#     taxa_table <- taxa_table[c(8,1,2,3,4,5,6,7)]
#   }else{
#     taxa_table <- taxa_table[c(7,1,2,3,4,5,6)]
#   }
#   taxa_table[] = lapply(taxa_table, blank2na, na.strings=c('','NA','na','N/A','n/a','NaN','nan'))
#   lastValue <- function(x) tail(x[!is.na(x)], 1)
#   last_taxons<- apply(taxa_table, 1, lastValue)
#   taxa_table$BestMatch <- last_taxons
#   taxa_table[, "BestMatch"] <- gsub("_", " ", taxa_table[, "BestMatch"])
#   taxa_table$Taxonomy <- paste(taxa_table$OTU, taxa_table$BestMatch, sep="-")
#   taxa_table[, "Genus"] <- gsub(" sp.", "", taxa_table[, "Genus"])
#   tax_table(dataframe) <- tax_table(as.matrix(taxa_table))
#   return(dataframe)
# }


