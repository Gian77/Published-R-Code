# ExtractCore function ------------------------------------------------------------------------------------
# Modofied form Stopnisek and Shade 2019 - Current Opinion in Microbiology
# https://www.sciencedirect.com/science/article/pii/S1369527419300426?via%3Dihub



# Extract random columns fucntion  ------------------------------------------------------------------------

# Example: SplitSampleDf(physeq)

# Thsi fucntion randomly extract (seed is set) the max number of samples (columns)
# shared in each in each group (variable) present in a phyloseq object 

SplitSampleDf <- function(physeq){
  set.seed(2182023)
  dplyr::count(as.data.frame(as.matrix(sample_data(physeq))), vars = Filt) -> sam_group
  # min shared samples
  min(sam_group$n) -> n
  #metadata
  meta <- as(phyloseq::sample_data(physeq), Class = "data.frame")
  meta <- data.frame(
    SID = phyloseq::sample_names(physeq),
    meta,
    stringsAsFactors = FALSE)
  # generate a list and sample n samples each element
  result <- plyr::dlply(.data = meta, .variables = "Filt", .fun = function(z){ z$SID })
  result <- plyr::llply(.data = result, .fun = function(z){sample(z, n)})
  # unlist and filter the original phyloseq
  unlist(result) -> sample_names
  otu_table(physeq) <- 
    subset(otu_table(physeq), select = c(as.vector(sample_names)))
  otu_table(physeq) <- 
    otu_table(physeq)[which(rowSums(otu_table(physeq)) > 0),] 
  dplyr::count(as.data.frame(as.matrix(sample_data(physeq))), vars = Filt) %T>% print()
  cat("\n")
  physeq %T>% print()
  return(physeq)
}



# Extract Core OTUs function ------------------------------------------------------------------------------

# Example: ExtractCore(phyloseq_object, Date", "elbow", Treatment", "Fertilized")

# This function will filter the dataset by Fertilized samples only and calculate
# the core according to the different sampling dates.

ExtractCore <- function(physeq, Var, method, Group=NULL, Level=NULL){
  require(magrittr)
  require(vegan)
  require(dplyr)
  require(tidyverse)
  library(lazyeval)
  #set.seed(100)
  # input dataset needs to be rarified and minimum depth included
  if (min(sample_sums(physeq)) == max(sample_sums(physeq))) {
    print("Dataset is rarefied at a depth of:")
    nReads = min(sample_sums(physeq))
    nReads %T>% print()
    rare <- physeq
    rare %T>% print()
    taxon <- as(tax_table(rare), "matrix")
    taxon <- as.data.frame(taxon)
    dim(taxon)  %T>% print()
  } else {
    print("Performing Rarefaction. Minimum depth:")
    nReads = min(sample_sums(physeq))
    nReads %T>% print()
    rare <-
      rarefy_even_depth(
        physeq, sample.size = nReads,
        trimOTUs = TRUE,
        rngseed = TRUE,
        replace = TRUE)
    taxon <- as(tax_table(rare), "matrix")
    taxon <- as.data.frame(taxon)
  }
  # choosing a subset or using the whole phyloseq object as is
  if (is.null(Group)) {
    print("Dataset not subsetted")
    otu <- rare@otu_table %>% as("matrix")
    map <- rare@sam_data %>% as("data.frame")
  } else{
    print("Dataset subsetted using the Group variable")
    sub_group <- substitute(Group)
    sub_set <- subset(sample_data(rare), eval(parse(text=sub_group)) %in% Level)
    physeq1 <- merge_phyloseq(otu_table(rare),
                              tax_table(rare),
                              refseq(rare),
                              sub_set)
    otu_table(physeq1) <- otu_table(physeq1)[which(rowSums(otu_table(physeq1)) > 0), ]
    otu <- physeq1@otu_table %>% as("matrix")
    map <- physeq1@sam_data %>% as("data.frame")
    print("Grouping Factor")
    map[,Group] %T>% print()
    taxon <- as(tax_table(physeq1), "matrix")
    taxon <- as.data.frame(taxon)
  }
  map$SampleID <- rownames(map)
  print("Check: dimension of datframe and metadata")
  dim(otu) %T>% print() # funcitons form magrittr package
  dim(map) %T>% print() 
  # calculating occupancy and abundance
  otu_PA <-
    1 * ((otu > 0) == 1)                                             # presence-absence data
  otu_occ <-
    rowSums(otu_PA) / ncol(otu_PA)                                   # occupancy calculation
  otu_rel <-
    apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)     # mean relative abundance
  occ_abun <-
    add_rownames(as.data.frame(cbind(otu_occ, otu_rel)), "otu")     # combining occupancy and abundance data frame
  # NOTE! add_rownames is deprecated and generates a warning, a bug of tidyverse, 
  # alternative you can use:
  # occ_abun <- tibble::rownames_to_column(as.data.frame(cbind(otu_occ, otu_rel)),"otu")
  # Ranking OTUs based on their occupancy
  # For caluclating raking index we included following conditions:
  # - time-specific occupancy (sumF) = frequency of detection within time point (genotype or site)
  # - replication consistency (sumG) = has occupancy of 1 in at least one time point (genotype or site) (1 if occupancy 1, else 0)
  Var <- enquo(Var) # lazy evaluation
  PresenceSum <-
    data.frame(otu = as.factor(row.names(otu)), otu) %>%
    gather(SampleID, abun,-otu) %>%
    left_join(map, by = "SampleID") %>%
    group_by(otu, !!Var) %>%
    dplyr::summarise(
      time_freq = sum(abun > 0) / length(abun), # frequency of detection between time points
      coreTime = ifelse(time_freq == 1, 1, 0)) %>% # 1 only if occupancy 1 with specific time, 0 if not
    group_by(otu) %>%
    dplyr::summarise(
      sumF = sum(time_freq),
      sumG = sum(coreTime),
      nS = length(!!Var)* 2,
      Index = (sumF + sumG) / nS) # calculating weighting Index based on number of points detected
  # PresenceSum %T>% print()
  # ranking otus
  otu_ranked <- occ_abun %>%
    left_join(PresenceSum, by = "otu") %>%
    transmute(otu = otu,
              rank = Index) %>%
    arrange(desc(rank))
  #otu_ranked %T>% print()
  # calculating BC dissimilarity based on the 1st ranked OTU
  otu_start = otu_ranked$otu[1]
  start_matrix <- as.matrix(otu[otu_start, ])
  start_matrix <- t(start_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x)
    sum(abs(start_matrix[, x[1]] - start_matrix[, x[2]])) / (2 * nReads))
  x_names <-
    apply(combn(ncol(start_matrix), 2), 2, function(x)
      paste(colnames(start_matrix)[x], collapse = "-"))
  # creating a data.frame and adding first OTU name as 1
  df_s <- data.frame(x_names,x)
  names(df_s)[2] <- 1 
  # Initialize your data structures:  calculating the contribution of ranked OTUs to the BC similarity
  BCaddition <- NULL
  BCaddition <- rbind(BCaddition,df_s)
  # calculating BC dissimilarity based on additon of ranked OTUs from 2nd to 500th. 
  # Can be set to the entire length of OTUs in the dataset.
  # it might take some time if more than 5000 OTUs are included.
  for(i in 2:1000){ #nrow(otu_ranked)
    otu_add=otu_ranked$otu[i]
    add_matrix <- as.matrix(otu[otu_add,])
    add_matrix <- t(add_matrix)
    start_matrix <- rbind(start_matrix, add_matrix)
    y <-
      apply(combn(ncol(start_matrix), 2), 2, function(y)
        sum(abs(start_matrix[, y[1]] - start_matrix[, y[2]])) / (2 * nReads))
    df_a <- data.frame(x_names, y)
    names(df_a)[2] <- i 
    BCaddition <- left_join(BCaddition, df_a, by=c('x_names'))
  }
  # Calculating the BC dissimilarity of the whole dataset (not needed if the second loop 
  # is already including all OTUs)
  z <-
    apply(combn(ncol(otu), 2), 2, function(z)
      sum(abs(otu[, z[1]] - otu[, z[2]])) / (2 * nReads))
  # overwrite the names here
  x_names <-
    apply(combn(ncol(otu), 2), 2, function(x)
      paste(colnames(otu)[x], collapse = "-"))
  df_full <- data.frame(x_names, z)
  names(df_full)[2] <- length(rownames(otu))
  BCfull <- left_join(BCaddition, df_full, by='x_names')
  # ranking the obtained BC 
  rownames(BCfull) <- BCfull$x_names
  temp_BC <- BCfull
  temp_BC$x_names <- NULL
  temp_BC_matrix <- as.matrix(temp_BC)
  BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
    gather(comparison, BC, -rank) %>%
    group_by(rank) %>%
    # Calculate mean Bray-Curtis dissimilarity
    summarise(MeanBC=mean(BC)) %>%            
    arrange(desc(-MeanBC)) %>%
    # Calculate proportion of the dissimilarity explained by the n number of ranked OTUs 
    mutate(proportionBC=MeanBC/max(MeanBC))
  #BC_ranked %T>% print()
  # Calculating the increase BC
  Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
  increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
  BC_ranked <- left_join(BC_ranked, increaseDF)
  BC_ranked <- BC_ranked[-nrow(BC_ranked),]
  if (method=="elbow"){
    fo_difference <- function(pos){
      left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
      right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
      return(left - right)
    }
    BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)
    elbow <- which.max(BC_ranked$fo_diffs)
    core_otus <- otu_ranked$otu[1:elbow]
    core_otus %T>% print()
  }
  # Creating threshold for core inclusion - last call method using a
  # final increase in BC similarity of equal or greater than 2%
  else{
    lastCall <-
      as.numeric(as.character(dplyr::last(
        subset(BC_ranked, IncreaseBC >= 1.02)$rank)))
    # lastCall <- last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))
    core_otus <- otu_ranked$otu[1:lastCall]
    core_otus %T>% print()
  }
  # Adding Core otus for reating occupancy abundance plot
  occ_abun$fill <- 'no'
  occ_abun$fill[occ_abun$otu %in% core_otus] <- "core"
  return_list <-
    list(core_otus, BC_ranked, otu_ranked, occ_abun, otu, map, taxon)
  return(return_list)
}


# PlotBCincrease funciton ---------------------------------------------------------------------------------

# Example: PlotBCincrease(core_otus, 400)

# This function will plot the last 2% BC increase in the dataset

PlotBCincrease <- function(BC_ranked, max){
  BC_ranked <- as.data.frame(BC_ranked[[2]])
  dim(BC_ranked) %T>% print()
  BC_ranked$rank <- factor(BC_ranked$rank, levels=BC_ranked$rank)
  BC_ranked <- BC_ranked[complete.cases(BC_ranked), ]
  plot <- ggplot(BC_ranked[1:max,], 
                 aes(x=rank[1:max], y=proportionBC)) +
    geom_point(pch=21, fill='white', alpha=0.25, size=1.5) +
    theme_classic() + 
    theme(strip.background = element_blank(),axis.text.x = element_text(size=6, angle=90)) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
    theme(axis.text.x = element_text(angle = 45, size = 7,hjust = 1, vjust = 1.15)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    scale_x_discrete(limits=BC_ranked$rank[1:max], # making the x axis more readable
                     breaks=BC_ranked$rank[1:max][seq(1,length(BC_ranked$rank[1:max]),by=10)]) +
    geom_vline(xintercept=last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.03)])), 
               lty=4, col="red", cex=0.5) +
    labs(x="ranked OTUs",y="Bray-Curtis similarity") +
    annotate(geom="text", 
             x=last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.03)])), 
             y=0.2, 
             label=paste("Last 2% increase\n(",last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=
                                                                                       1.03)]))," OTUs)", 
                         sep=""), color="red", size=3, )
  return(plot)
}

# FitNeutral function ------------------------------------------------------------

# Example: FitNeutral(core_otus) 

# Fitting neutral models through the sncm.fit.R function modified form 
# Burns et al. - 2016 ISME JOURNAL

FitNeutral <- function(otu_core) {
  source("../R_functions/sncm.fit.R")
  #generating df
  taxon <- otu_core[[7]]
  spp <- t(otu_core[[5]])
  occ_abun <- otu_core[[4]]
  names(occ_abun)[names(occ_abun) == "otu"] <- "OTU_ID"
  occ_abun %T>% print()
  # source community pool
  meta <- otu_core[[6]]
  # using root as the pool gave weird results
  # soil_pool <- spp[ rownames(subset(meta, Origin=="Soil")), ]
  # root_pool <- spp[ rownames(subset(meta, Origin=="Root")), ]
  #fitting model
  obs.np <- sncm.fit(spp, taxon, stats = FALSE, pool = NULL)
  sta.np <- sncm.fit(spp, taxon, stats = TRUE, pool = NULL)
  #
  obs.np$fit_class <- "As predicted"
  obs.np[which(obs.np$freq < obs.np$pred.lwr), "fit_class"] <-
    "Below prediction"
  obs.np[which(obs.np$freq > obs.np$pred.upr), "fit_class"] <-
    "Above prediction"
  obs.np[which(is.na(obs.np$freq)), "fit_class"] <- "NA"
  as.data.frame(left_join(occ_abun, obs.np, by = 'OTU_ID')) -> fit_table
  #
  sta.np$above.pred <-
    sum(obs.np$freq > (obs.np$pred.upr), na.rm = TRUE) / sta.np$Richness
  sta.np$below.pred <-
    sum(obs.np$freq < (obs.np$pred.lwr), na.rm = TRUE) / sta.np$Richness
  fit_res <- as.data.frame(sta.np)
  rownames(fit_res) <- "Value"
  fit_res
  list_tab <- list(fit_res, fit_table)
  return(list_tab)
}


# PlotNeutral function ----------------------------------------------------------------

#Example: PlotNeutral(core_otus)

# This function will plot the result of a neutral model fit, and will color
# the core otus that are driven by the host (bove the model) or dispersal limited
# (below the model). The neutral OTUs are plotted as grey circles.

PlotNeutral <- function(obs){
  obs1 <- 
    as.data.frame(obs[2])
  obs1 <- obs1[!is.na(obs1$p), ]
  obs2 <- 
    as.data.frame(obs[1])
  plotn <- ggplot(obs1, aes(x=log10(otu_rel), y=otu_occ)) +
    # geom_point(data=subset(obs1, fill%in%"no" | fit_class%in%"As predicted"),
    #            aes(x=log10(otu_rel), y=otu_occ), pch=21, fill="white", alpha=0.25, size=1.5) +
    geom_point(data=subset(obs1, fill%in%"no" | fit_class%in%"As predicted"),
               aes(x=log10(otu_rel), y=otu_occ), pch=21, fill="white", alpha=0.25, size=0.8) +
    geom_point(data=subset(obs1, fill%in%"core" & fit_class%in%"As predicted"),
               aes(x=log10(otu_rel), y=otu_occ), pch=21, fill="gray60", color="gray60", size=0.8) +
    geom_point(data=subset(obs1, fill%in%"core" & fit_class%in%"Above prediction"), 
               aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='#CB54D6', color= '#CB54D6',size=0.8) +
    geom_point(data=subset(obs1, fill%in%"core" & fit_class%in%"Below prediction"), 
               aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='#09104d', color= '#09104d', size=0.8) +
    geom_line(color='red', data=obs1, size=0.7, aes(y=obs1$freq.pred, x=log10(obs1$p)), alpha=0.25) +
    geom_line(color='black', lty='twodash', size=0.7, data=obs1, aes(y=obs1$pred.upr, x=log10(obs1$p)), alpha=0.25)+
    geom_line(color='black', lty='twodash', size=0.7, data=obs1, aes(y=obs1$pred.lwr, x=log10(obs1$p)), alpha=0.25)+
    labs(x="Log10(mean abundance)", y="Occupancy") + #mean relative abundance
    theme_classic() +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) +
    theme(axis.text.x = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 7)) +
    theme(legend.position = "right") +
    annotate("text", -Inf, Inf, label = paste("italic(R)^2 ==", round(obs2$Rsqr, 3)),
             parse = TRUE, size = 2.2, hjust = -0.2, vjust = 1.2) +
    annotate("text", -Inf, Inf, label = paste("italic(m) ==", round(obs2$m, 3)),
             parse = TRUE, size = 2.2, hjust = -0.2, vjust = 3.2)
  return(plotn)
}







