# Microbial Netwoorks ------------------------------------------------------------------------------------
# Set of functions to generate microbial networks.

# Note, you will have to start with the already filtered datasets,
# containing the OTUs you plan to be present in the networks.

# Example: MakeSENet(physeq_object1, physeq_object2) 

# This function will generate a network starting form two dataset, one
# form bactera and one form fungi using the SpiecEasi package. If the 
# physeq_object2 is not present the function will not use it in the 
# claculations


MakeSENet <- function(fungi_core, prok_core) {
  require(SpiecEasi)
  param <- list(rep.num=100, seed=10010, ncores=12, thresh=0.03) # increase thresh for denser
  if (!missing(prok_core)) {
  #subset to equal samples
  dplyr::intersect(rownames(sample_data(fungi_core)),
                   rownames(sample_data(prok_core))) -> keep
  subset(otu_table(fungi_core), select = keep) -> otu_table(fungi_core)
  otu_table(fungi_core) <-
    otu_table(fungi_core)[which(rowSums(otu_table(fungi_core)) > 0),]
  subset(otu_table(prok_core), select = keep) -> otu_table(prok_core)
  otu_table(prok_core) <-
    otu_table(prok_core)[which(rowSums(otu_table(prok_core)) > 0),]
  fungi_core %T>% print()
  prok_core %T>% print()
  # make network
  spiec_net <- SpiecEasi::spiec.easi(list(fungi_core,prok_core),
                                     method="mb",
                                     lambda.min.ratio=1e-1, # the higher the more edges, higher density (5e-1 to 1e-1, 1e-2, 1e-3)
                                     nlambda=100, # number of slices between minimum amount of edges to maximum
                                     sel.criterion ="stars", # selects the lambda and modules the most stable connections
                                     pulsar.select = TRUE,
                                     pulsar.params=param) #iteration of starts, number of times sample the 80% of the data within each lambda
  print("Optimal Lambda")
  getOptLambda(spiec_net) %T>% print()
  print("Network stability - 0.05 optimum")
  getStability(spiec_net) %T>% print()
  } else {
  spiec_net <- SpiecEasi::spiec.easi(fungi_core,
                                     method="mb",
                                     lambda.min.ratio=1e-1, # the higher the more edges, higher density (5e-1 to 1e-1, 1e-2, 1e-3)
                                     nlambda=100, # number of slices between minimum amount of edges to maximum
                                     sel.criterion ="stars", # selects the lambda and modules the most stable connections
                                     pulsar.select = TRUE,
                                     pulsar.params=param) #iteration of starts, number of times sample the 80% of the data within each lambda
  print("Optimal Lambda")
  getOptLambda(spiec_net) %T>% print()
  print("Network stability - 0.05 optimum")
  getStability(spiec_net) %T>% print()
  }
  return(spiec_net)
}


# GetNetwork -----------------------------------------------------------------------------

# Example: GetNetwork(physeq_object1, physeq_object2, MakeSENet_object) 

# this function will genrate the network object and add edge weigths


GetNetwork <- function(physeq_f, physeq_b, SpiecEasi_object){
  require(igraph)
  if (!missing(physeq_b)) {
  names_taxa <- c(taxa_names(physeq_f), taxa_names(physeq_b))
  network <- adj2igraph(Matrix::drop0(getRefit(SpiecEasi_object)),
                        vertex.attr = list(name=names_taxa),
                        rmEmptyNodes = FALSE)
  symbeta_net = symBeta(getOptBeta(SpiecEasi_object), mode='maxabs')
  weight <- Matrix::summary(t(symbeta_net))[,3]
  E(network)$weight <- weight
  } else {
  names_taxa <- c(taxa_names(physeq_f))
  network <- adj2igraph(Matrix::drop0(getRefit(SpiecEasi_object)),
                          vertex.attr = list(name=names_taxa),
                          rmEmptyNodes = FALSE)
  symbeta_net = symBeta(getOptBeta(SpiecEasi_object), mode='maxabs')
  weight <- Matrix::summary(t(symbeta_net))[,3]
  E(network)$weight <- weight
  }
  return(network)
}



# TesNULL - compare to a Random network ----------------------------------------------------------

#Example: TesToNULL(network)

# This funcion will compare the network degree distribution to that of a random network
# with the same vertes and edges number using the erdos.renyi algorithm. The comparison 
# is performed 100 times and difference tested with Kolmogorov–Smirnov, ks.test() 
# Degree distribution of the two networks is also plotted.

# NOTE, this may not benecessary because SpiecEasi networks are themselves based on a
# repeated subsampling procedure, resulting in edge stability values. It's doing 
# subsampling internally!

TesToNULL <- function(network){
  # network properties
  v <- gorder(network)
  e <- gsize(network)
  random <- erdos.renyi.game(v, e, type = "gnm", directed = FALSE)
  x <- degree_distribution(random) # this gets the degree distribution of the random graph 
  y <- degree_distribution(network) # get degree distribution of actual network
  print("Plotting Degree Distribution")
  par(mfrow=c(2, 2), mar=c(2.5,2.5,2,1))
  plot(density(degree_distribution(random)), main="Random") #, xlim=range(-0.1, 0.4), ylim=range(0, 14))
  plot(density(degree_distribution(network)), main="Real") #, xlim=range(-0.1, 0.4), ylim=range(0, 14))
  plot(degree(random), xlab="Vertices", ylab = "Degree")
  plot(degree(network), xlab="Vertices", ylab = "Degree")
  print(ks.test(x,y))
  df <- NULL
  n <- 100
  for(i in 1:n){
    #newmatrix <- sample_pa(n = v, start.graph =network,out.pref = TRUE)
    newmatrix <- erdos.renyi.game(v, e, type = "gnm", directed = FALSE)
    dist <- degree_distribution(newmatrix)
    df<-rbind(df, dist)
  }
  # test against 100 random networks
  ks.test_res <- NULL
  for (i in 1:nrow(df)){
    ks.test_new <- ks.test(df[i,],y)
    ks.test_res$D[i] <- ks.test_new$statistic
    ks.test_res$p[i] <- ks.test_new$p.value
  }
  ks.test_res$p -> pval
  length(which(pval<=0.05))/length(pval)*100 -> Prop
  ks.test_res <- c(ks.test_res, Prop = list(Prop))
  return(ks.test_res)
}



TesToRandom <- function(network, nodes) {
  y <-
    degree_distribution(network) # get degree distribution of actual network
  degamat_con_below <- NULL
  n <- 100
  for (i in 1:n) {
    newmatrix <- sample_pa(n = nodes, start.graph = network, out.pref = TRUE, directed=FALSE)
    degmat_con_below <- degree_distribution(newmatrix)
    degamat_con_below <- rbind(degamat_con_below, degmat_con_below)
  }
  # test against 100 random networks
  ks.test_con_below <- NULL
  for (i in 1:nrow(degamat_con_below)) {
    ks.test_new <- ks.test(degamat_con_below[i, ], y, alternative = "greater")
    ks.test_con_below$D[i] <- ks.test_new$statistic
    ks.test_con_below$p[i] <- ks.test_new$p.value
  }
  return(ks.test_con_below)
}



# REMOVING WEAK EDGES --------------------------------------------------------

#Example: FilterEdge(network) 

# This function will remove weak edges based on the standard deviaiton form the 
# mean edge weights, re-centere to 0.

FilterEdge <- function(network){
  mean <- 
    mean(E(network)$weight)
  sd <-
    sd(E(network)$weight)
  network_new <- 
    delete_edges(network, which(E(network)$weight <= 1 * sd &
                                  E(network)$weight >= 1 *-sd))
  par(mfrow=c(1, 2))
  # Before
  hist(E(network)$weight, breaks = 100, 
       col="lightblue", border="lightblue", 
       xlab="Weight", ylab="Frequency", 
       main = "Original\nWeight Distribution")
  abline(v=mean(E(network)$weight), col="red", lty = 1, lwd=2)
  abline(v=sd(E(network)$weight), col="blue", lty = 2, lwd=2)
  abline(v=-sd(E(network)$weight), col="blue", lty = 2, lwd=2)
  # After
  hist(E(network_new)$weight, breaks = 100, 
       col="lightblue", border="lightblue", 
       xlab="Weight", ylab="Frequency",
       main = "Filtered\nWeight Distribution")
  abline(v=mean(E(network_new)$weight), col="red", lty = 1, lwd=2)
  # abline(v=sd(E(network_new)$weight), col="blue", lty = 2, lwd=2)
  # abline(v=-sd(E(network_new)$weight), col="blue", lty = 2, lwd=2)
  return(network_new)
}



# Plot Pi and Zi function ---------------------------------------------------------------------------

# Example: PlotZiPi(df_zipi)

# this function will generate a classic Pi-Zi plot based on the definiton 

PlotZiPi <- function(df){
  require(ggrepel)
  plot_between <- ggplot() +
    geom_point(data = df, aes(x = P, y = Z, size=Degree), shape=16, alpha = 0.5) +
    geom_point(data= df,
               aes(x = P, y = Z, size=Degree), shape=16, alpha = 0.5) +
    #scale_size_continuous(range = c(min(nodes$Degree)/10, max(nodes$Degree)/10)) +
    geom_text_repel(data = subset(df, df$Hubs=="key"),
                    aes(x = Betweenness, y = Degree,
                        label = OTU_ID), size=2.5, force = 10, color="red") +
    # geom_vline(xintercept = percent_95_betw, linetype = "dashed", colour = "blue") +
    # geom_hline(yintercept = percent_95_degree, linetype = "dashed", colour = "blue") +
    #geom_text_repel(data= subset(df, !Key%in%c("Peripherals") & Degree>10),
    #                aes(x = P, y = Z, label = BestMatch), size=2.5, force = 10, color="blue") +
    #scale_color_manual(values = palette1) +
    geom_vline(xintercept = 0.62, linetype = "dashed", colour = "red") +
    geom_hline(yintercept = 2.5, linetype = "dashed", colour = "red") +
    theme_classic() +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) +
    theme(axis.text.x = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm")) +
    theme(legend.position="right") +
    guides(color = guide_legend(override.aes = list(size = 4)),
           size = FALSE) +
    labs(x="Among module connectivity (Pi)", 
         y="Within module connectivity (Zi)", 
         title=paste("Pi-Zi Plot"))
  return(plot_between)
}

# ************************************************************************------------------
# EXTRAXTING NETWORK ATTRIBUTES -------------------------------------------------------------

#Extract Taxonomy function --------------------------------------------------------------------

#Example: ExtarctTaxa(physeq_fungi, physeq_prok) 

# Generating Taxonomy

ExtarctTaxa <- function(physeq_f, physeq_b){
  dplyr::intersect(rownames(sample_data(physeq_f)),
                   rownames(sample_data(physeq_b))) -> to_keep
  subset(otu_table(physeq_f), select = to_keep) -> otu_table(physeq_f)
  otu_table(physeq_f) <-
    otu_table(physeq_f)[which(rowSums(otu_table(physeq_f)) > 0),]
  subset(otu_table(physeq_b), select = to_keep) -> otu_table(physeq_b)
  otu_table(physeq_b) <-
    otu_table(physeq_b)[which(rowSums(otu_table(physeq_b)) > 0),]
  rbind(tax_table(physeq_f), 
        tax_table(physeq_b))-> taxa
  taxa <- as.data.frame(taxa)
  taxa <- taxa[, c(1:3, 7, 9:10)]
  taxa$Phylum <-
    ifelse(is.na(taxa$Phylum),"Unclassified",
           paste(taxa$Phylum))
  taxa$Genus <-
    ifelse(is.na(taxa$Genus),"Unclassified",
           paste(taxa$Genus))
  return(taxa)
}


# Genrate composite phyloseq objects  -------------------------------------------------------------

# Example: GeneratePhyseqNet(physeq_fungi, physeq_prok, "Group", "Level") 

GeneratePhyseqNet <- function(physeq_f, physeq_b, Group, Level) {
  require(lazyeval)
  physeq_f %>% transform_sample_counts(function(x) {
    x / sum(x) *100}) -> physeq_f
  physeq_b %>% transform_sample_counts(function(x) {
    x / sum(x) *100}) -> physeq_b
  dplyr::intersect(rownames(sample_data(physeq_f)),
                   rownames(sample_data(physeq_b))) -> to_keep
  print("Number of shared samples")
  length(to_keep) %T>% print()
  subset(otu_table(physeq_f), select = to_keep) -> otu_table(physeq_f)
  subset(otu_table(physeq_b), select = to_keep) -> otu_table(physeq_b)
  merge_phyloseq(physeq_f, physeq_b) -> physeq_fb
  sub_group <- substitute(Group)
  sub_set <- subset(sample_data(physeq_fb), eval(parse(text=sub_group)) %in% Level)
  physeq_filt <- merge_phyloseq(otu_table(physeq_fb),
                                tax_table(physeq_fb),
                                refseq(physeq_fb),
                                sub_set)
  otu_table(physeq_filt) <-
    otu_table(physeq_filt)[which(rowSums(otu_table(physeq_filt)) > 0), ]
  physeq_filt %T>% print()
  return(physeq_filt)
}

# Calculate Abundances ----------------------------------------------------------

#Example: ExtractAbund(physeq_fb_pop, "Poplar", "0-10cm")

# This funciton generates relative abundance for each group

ExtractAbund <- function(physeq_crop, Crop, Var){
  AbundCrop <- as.data.frame(taxa_sums(physeq_crop))
  AbundCrop <-
    data.frame(OTU_ID = as.factor(row.names(AbundCrop)), AbundCrop)
  colnames(AbundCrop)[2] <- "AbundCrop"
  #head(AbundCrop) %T>% print() 
  # subset to Depth
  sub_set <- subset(sample_data(physeq_crop), Section %in% Var)
  physeq_filt <- merge_phyloseq(otu_table(physeq_crop),
                                tax_table(physeq_crop),
                                refseq(physeq_crop),
                                sub_set)
  otu_table(physeq_filt) <-
    otu_table(physeq_filt)[which(rowSums(otu_table(physeq_filt)) > 0), ]
  AbundDepth <- as.data.frame(taxa_sums(physeq_filt))
  AbundDepth <-
    data.frame(OTU_ID = as.factor(row.names( AbundDepth)),  AbundDepth)
  colnames(AbundDepth)[2] <- "AbundDepth"
  # Niche
  AbundDepth$RootAbund <- physeq_filt %>%
    subset_samples(Origin %in% c("Root")) %>%
    taxa_sums()
  AbundDepth$SoilAbund <- physeq_filt %>%
    subset_samples(Origin %in% c("Soil")) %>%
    taxa_sums()
  # join the df
  abund_total <-
    dplyr::inner_join(AbundCrop, AbundDepth, by = "OTU_ID")
  rownames(abund_total) <- abund_total$OTU_ID
  head(abund_total) %T>% print() 
  # converting in % of the AbundCrop
  abund_total$RootAbund <- abund_total$RootAbund * 100 /abund_total$AbundDepth
  abund_total$SoilAbund <- abund_total$SoilAbund * 100 / abund_total$AbundDepth
  abund_total$AbundDepth <- abund_total$AbundDepth * 100 / abund_total$AbundCrop
  abund_total$Depth <- Var
  colnames(abund_total)[6] <- Crop
  #colnames(abund_total)[3] <- Var
  #abund_total$Crop <- rep(Crop, nrow(abund_total))
  return(abund_total)
}


# Calculating zipi for finding keystone microbes ----------------------------------------------------

#Example ZiPi()

ZiPi <- function(netw="network",modules="modules"){
  names= V(netw)$name
  total_connections=c()
  module_connections=c()
  number_of_modules=c()
  meanZ=c()
  sdZ=c()
  Z=c()
  P=c()
  for(i in 1:length(names)){
    total_connections[i]=sum(netw[i]>0)
    module_connections[i]=sum(netw[i][which(modules==modules[i])]>0)
    KitKi=c()
    for(j in 1:length(unique(modules))){
      KitKi[j]=((sum(netw[i][which(modules==j)]))/total_connections[i])^2
      
    }
    P[i]=1-sum(KitKi)
  }
  for(i in 1:length(unique(modules))){
    meanZ[i]=mean(module_connections[which(modules==i)])
    sdZ[i]=sd(module_connections[which(modules==i)])
  }
  for(i in 1:length(names)){
    Z[i]=(module_connections[i]-meanZ[modules[i]])/sdZ[modules[i]]
    
  }
  return(cbind.data.frame(names,"module"=modules,module_connections,total_connections,Z,P))
}



# Calculate vertex attribtues -----------------------------------------

#Example: head(CalcNetAttr(network_Pop_10, "10cm")) 

# This Funciton caluclated natwork attribuets. Also generate and save a 
# df_PiZi dataframe for further plotting.

CalcNetAttr <- function(network, Var1, Var2){
  require(igraph)
  require(microbiome)
  nodes <-
    data.frame(matrix(ncol = 6, nrow = length(V(network)$name)))
  colnames(nodes) <-
    c("Degree","Modules","Betweenness",
      "Niche_ind","Crop_ind","Depth_ind")
  rownames(nodes) <- V(network)$name
  nodes <- data.frame(OTU_ID = as.factor(row.names(nodes)), nodes)
  #nodes <- data.frame(Network = rep(Var, nrow(nodes)), nodes)
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
  # Indicators Niche
  rbind(ind_fungi_niche, ind_prok_niche) -> ind_fb_niche
  ind_fb_niche[rownames(ind_fb_niche) %in% rownames(nodes), ] -> net_ind_df_niche
  nodes[rownames(nodes) %in%
          rownames(net_ind_df_niche[net_ind_df_niche$s.Root =="1", ]), ]$Niche_ind <- "Root"
  nodes[rownames(nodes) %in%
          rownames(net_ind_df_niche[net_ind_df_niche$s.Soil == "1", ]), ]$Niche_ind <-"Soil"
  # Indicators Crop
  rbind(ind_fungi_crop, ind_prok_crop) -> ind_fb_crop
  ind_fb_crop[rownames(ind_fb_crop) %in% rownames(nodes), ] -> net_ind_df_crop
  if (any(net_ind_df_crop$s.Poplar=="1") ==TRUE){ # this becasue there are no indicators for Poplar
  nodes[rownames(nodes) %in%
          rownames(net_ind_df_crop[net_ind_df_crop$s.Poplar == "1", ]), ]$Crop_ind <- "Poplar"
  }else{}
  if (any(net_ind_df_crop$s.Prairie=="1") ==TRUE){ # this becasue there are no indicators for Poplar
  nodes[rownames(nodes) %in%
          rownames(net_ind_df_crop[net_ind_df_crop$s.Prairie == "1", ]), ]$Crop_ind <- "Prairie"
  }else{}
  if (any(net_ind_df_crop$s.Switchgras=="1") ==TRUE){ # this becasue there are no indicators for Poplar
  nodes[rownames(nodes) %in%
          rownames(net_ind_df_crop[net_ind_df_crop$s.Switchgrass == "1", ]), ]$Crop_ind <- "Switchgrass"
  }else{}
  # Indicators Depth
  rbind(ind_fungi_depth, ind_prok_depth) -> ind_fb_depth
  ind_fb_depth[rownames(ind_fb_depth) %in% rownames(nodes), ] -> net_ind_df_depth
  if (any(net_ind_df_depth$s.10=="1") ==TRUE){ # this becasue there are no indicators at 0-10cm
    nodes[rownames(nodes) %in%
            rownames(net_ind_df_depth[net_ind_df_depth$s.10 == "1", ]), ]$Depth_ind <-"0-10cm"
  }else{}
  if (any(net_ind_df_depth$s.25=="1") ==TRUE){ # this becasue there are no indicators at 10-25cm
    nodes[rownames(nodes) %in%
            rownames(net_ind_df_depth[net_ind_df_depth$s.25 == "1", ]), ]$Depth_ind <-"10-25cm"
  }else{}
  if (any(net_ind_df_depth$s.50=="1") ==TRUE){ # this becasue there are no indicators at 25-50cm
    nodes[rownames(nodes) %in%
            rownames(net_ind_df_depth[net_ind_df_depth$s.50 == "1", ]), ]$Depth_ind <-"25-50cm"
  }else{}
  if (any(net_ind_df_depth$s.100=="1") ==TRUE){ # this becasue there are no indicators at 50-100cm
    nodes[rownames(nodes) %in%
            rownames(net_ind_df_depth[net_ind_df_depth$s.100 == "1", ]), ]$Depth_ind <-"50-100cm"
  }else{}
  #Hubs
  df_zipi <- ZiPi(network_new, modules=nodes$Modules)
  df_zipi$Key  <- ifelse(df_zipi$P>=0.62, "Connectors", NA)
  df_zipi$Key  <- ifelse(df_zipi$Z>=2.5, "Module hubs", paste(df_zipi$Key))
  df_zipi$Key  <- ifelse(df_zipi$P>=0.62 & df_zipi$Z>=2.5, "Network hubs", paste(df_zipi$Key))
  df_zipi$Key  <- ifelse(df_zipi$P<0.62 & df_zipi$Z<2.5, "Peripherals", paste(df_zipi$Key))
  df_zipi$Key  <- ifelse(is.na(df_zipi$Key), "Peripherals", paste(df_zipi$Key))
  head(df_zipi) %T>% print()
  if (identical(as.character(nodes$OTU_ID), as.character(df_zipi$names))) {
    print("dataframes order is identical!")
    nodes$Key <- df_zipi$Key
  } else {
    print("dataframes order is not identical! Reordering.")
    sample_order <- match(as.character(df_zipi$names), as.character(nodes$OTU_ID))
    df_zipi <- df_zipi[sample_order, ]
    nodes$Key <- df_zipi$Key
  }
  df_zipi <-
    left_join(df_zipi, nodes)
  df_zipi %T>% assign(paste("df_zipi", Var1, Var2, sep = "_"),., envir = .GlobalEnv) # saving the df
  return(nodes)
}


# NODE TABLE ----------------------------------------------------------------------

#Example: CalcNodes(abund_Pop_10, network_Pop_10, nfit_fungi_Pop_10, nfit_prok_Pop_10, "0-10cm") 

# This function will generate the node table

CalcNodes <- function(abund, network, neutral_fit_f=NULL, neutral_fit_b=NULL, Var1,Var2) {
  CalcNetAttr(network, Var1,Var2) -> nodes
  inner_join(abund, nodes, by = "OTU_ID") -> nodes_abund
  inner_join(nodes_abund, taxa_all, by = "OTU_ID") -> nodes_attr
  rownames(nodes_attr) <- nodes_attr$OTU_ID
  if ((!is.null(neutral_fit_f) | !is.null(neutral_fit_b))){
    rbind(neutral_fit_f[[2]], neutral_fit_b[[2]]) -> neutral_fit
    if (identical(neutral_fit$OTU_ID, rownames(nodes_attr))){
    }else{
      sample_order <- match(rownames(nodes_attr), neutral_fit$OTU_ID)
      neutral_fit <- neutral_fit[sample_order, ]
      nodes_attr$fit_class <- neutral_fit$fit_class
    }
  }
  #list(as_tibble(nodes), as_tibble(abund)) -> nodes_attr
  return(nodes_attr)
}



# EDGE TABLE --------------------------------------------------------------

#Example: CalcEdges(network_Pop_10, nodes_Pop_10cm,"0-10cm") 

CalcEdges <- function(network, node, Var){
  edge <- as.data.frame(ends(network, E(network)))
  edge$Depth <- rep(Var, gsize(network))
  edge$Weight <- E(network)$weight
  edge$Direction <- ifelse(E(network)$weight>0,"positive", "negative")
  edge$V1_taxa <- node$Kingdom[match(edge$V1, rownames(node))]
  edge$V2_taxa <- node$Kingdom[match(edge$V2, rownames(node))]
  edge$V1_Phylum <- node$Phylum[match(edge$V1, rownames(node))]
  edge$V2_Phylum <- node$Phylum[match(edge$V2, rownames(node))]
  edge$InterKing <- ifelse(edge$V1_taxa == edge$V2_taxa,
                           "Intrakingdom", "Interkingdom")
  edge$InterKingTaxa <- paste(edge$V1_taxa, edge$V2_taxa, sep = "-")
  edge$InterKingTaxa[edge$InterKingTaxa == "Bacteria-Archaea"] <-
    "Archaea-Bacteria"
  edge$AbundV1_Soil <- node$SoilAbund[match(edge$V1, rownames(node))]
  edge$AbundV2_Soil <- node$SoilAbund[match(edge$V2, rownames(node))]
  edge$AbundV1_Root <- node$RootAbund[match(edge$V1, rownames(node))]
  edge$AbundV2_Root <- node$RootAbund[match(edge$V2, rownames(node))]
  edge$SoilV1V2_Sum <- edge$AbundV1_Soil + edge$AbundV2_Soil
  edge$RootV1V2_Sum <- edge$AbundV1_Root + edge$AbundV2_Root
  return(edge)
}


# NETWORK STATS SUMMARY -------------------------------------------------------------------

# Example: NetStats(nodes_Pop_10cm, edges_Pop_10cm, network_Pop_10, "0-10cm")

# This funciton generates most common stats to define the characterisitcs 
# of a network

NetStats <- function(node, edge, network, neutral_fit_f, neutral_fit_b, Var){
  require(igraph)
  # creating dataframe
  net_summary <- data.frame(matrix(ncol=1, nrow =41))
  colnames(net_summary) <- Var
  rownames(net_summary) <- c("Nodes", "Edges", "Positive_Edges", "Negative_Edges",
                             "Ratio",
                             "Above", "Neutral", "Below",
                             "Above_all", "Neutral_all", "Below_all",
                             "Transitivity","Modularity","Modules", "Average_Module_size", 
                             "Module_SD", "Average_Degree","Degree_SD", "Average_Betweenes", 
                             "Betweenness_SD", "Path_length", "Edge_density", "Network_diameter",
                             "Module_hubs", "Network_hubs", "Conncetors", "Peripherals",
                             "Fungal_OTUs","Bacterial_OTUs","Archaeal_OTUs",
                             "Interkind_positive","Interkind_negative", 
                             "Root_indicators","Soil_indicators", 
                             "Pop_ind", "Pra_ind", "Switch_ind",
                             "10cm_ind","25cm_ind","50cm_ind","100cm_ind")
  # adding metrics
  net_summary[,Var][1] <- length(rownames(node))
  net_summary[,Var][2] <- length(rownames(edge))
  net_summary[,Var][3] <- nrow(edge[edge$Weight>0,])
  net_summary[,Var][4] <- nrow(edge[edge$Weight<0,])
  net_summary[,Var][5] <- net_summary[3,]/net_summary[4,]
  net_summary[,Var][6] <- table(node$fit_class)[[1]]/nrow(node)*100 # this is just the core
  net_summary[,Var][7] <- table(node$fit_class)[[2]]/nrow(node)*100
  net_summary[,Var][8] <- table(node$fit_class)[[3]]/nrow(node)*100
  if (is.null(neutral_fit_b)){
    neutral_fit_b[[2]] -> neutral_fit
  }else{
    rbind(neutral_fit_f[[2]], neutral_fit_b[[2]]) -> neutral_fit
  }
  net_summary[,Var][9] <- table(neutral_fit$fit_class)[[1]]/nrow(neutral_fit)*100 # this is for all taxa
  net_summary[,Var][10] <- table(neutral_fit$fit_class)[[2]]/nrow(neutral_fit)*100
  net_summary[,Var][11] <- table(neutral_fit$fit_class)[[3]]/nrow(neutral_fit)*100
  net_summary[,Var][12] <- transitivity(network, type = "average") # clustering coefficient
  net_summary[,Var][13] <- modularity(network, node$Modules) # modularity with no weights
  net_summary[,Var][14] <- length(unique(node$Modules))
  net_summary[,Var][15] <- mean(table(node$Modules))
  net_summary[,Var][16] <- sd(table(node$Modules))
  net_summary[,Var][17] <- mean(node$Degree)
  net_summary[,Var][18] <- sd(node$Degree)
  net_summary[,Var][19] <- mean(node$Betweenness)
  net_summary[,Var][20] <- sd(node$Betweenness)
  net_summary[,Var][21] <- mean_distance(network) # Path Length
  net_summary[,Var][22] <- edge_density(network, loops = FALSE) # Density
  # removing weights
  network_new <- delete_edge_attr(network, "weight") # delete weights
  net_summary[,Var][23] <- diameter(network_new, directed = TRUE, unconnected = TRUE) #diameter
  net_summary[,Var][24] <- length(rownames(subset(node, Key%in%"Module hubs")))
  net_summary[,Var][25] <- length(rownames(subset(node, Key%in%"Network hubs")))
  net_summary[,Var][26] <- length(rownames(subset(node, Key%in%"Connectors")))
  net_summary[,Var][27] <- length(rownames(subset(node, Key%in%"Peripherals")))
  net_summary[,Var][28] <- length(rownames(node[node$Kingdom=="Fungi",]))
  net_summary[,Var][29] <- length(rownames(node[node$Kingdom=="Bacteria",]))
  net_summary[,Var][30] <- length(rownames(node[node$Kingdom=="Archaea",]))
  net_summary[,Var][31] <- length(rownames(subset(edge[edge$InterKingTaxa=="Fungi-Bacteria" | 
                                                         edge$InterKingTaxa=="Archaea-Bacteria"|
                                                         edge$InterKingTaxa=="Fungi-Archaea", ], Weight>0)))
  net_summary[,Var][32] <- length(rownames(subset(edge[edge$InterKingTaxa=="Fungi-Bacteria" | 
                                                         edge$InterKingTaxa=="Archaea-Bacteria"|
                                                         edge$InterKingTaxa=="Fungi-Archaea", ], Weight<0)))
  net_summary[,Var][33] <- length(rownames(subset(node, Niche_ind%in%"Root")))
  net_summary[,Var][34] <- length(rownames(subset(node, Niche_ind%in%"Soil")))
  net_summary[,Var][35] <- length(rownames(subset(node, Crop_ind%in%"Poplar")))
  net_summary[,Var][36] <- length(rownames(subset(node, Crop_ind%in%"Prairie")))
  net_summary[,Var][37] <- length(rownames(subset(node, Crop_ind%in%"Switchgrass")))
  net_summary[,Var][38] <- length(rownames(subset(node, Depth_ind%in%"0-10cm")))
  net_summary[,Var][39] <- length(rownames(subset(node, Depth_ind%in%"10-25cm")))
  net_summary[,Var][40] <- length(rownames(subset(node, Depth_ind%in%"25-50cm")))
  net_summary[,Var][41] <- length(rownames(subset(node, Depth_ind%in%"50-100cm")))
  # reduce decimals
  net_summary <-
    apply(net_summary, 2, function(x) round(x, digits = 3))
  return(net_summary)
}



# JUST CORE STATS ---------------------------
NetStatsCore <- function(node, edge, network, neutral_fit_f, neutral_fit_b, Var1, Var2){
  require(igraph)
  # creating dataframe
  network_new <- 
    igraph::delete.vertices(network, which(degree(network)==0)) # delete single nodes
  CalcNetAttr(network_new, Var1, Var2) -> node_new
  vcount(network_new) %T>% print()
  dim(node_new) %T>% print()
  net_summary <- data.frame(matrix(ncol=1, nrow =29))
  colnames(net_summary) <- Var1
  rownames(net_summary) <- c("Nodes", "Edges", "Positive_Edges", "Negative_Edges",
                             "Ratio",
                             "Above", "Neutral", "Below",
                             "Above_all", "Neutral_all", "Below_all",
                             "Transitivity","Modularity","Modules", "Average_Module_size", 
                             "Module_SD", "Average_Degree","Degree_SD", "Average_Betweenes", 
                             "Betweenness_SD", "Path_length", "Edge_density", "Network_diameter",
                             "Module_hubs", "Network_hubs", "Conncetors", "Peripherals",
                             "Interkind_positive","Interkind_negative")
  # adding metrics
  net_summary[,Var1][1] <- length(rownames(node_new))
  net_summary[,Var1][2] <- length(rownames(edge))
  net_summary[,Var1][3] <- nrow(edge[edge$Weight>0,])
  net_summary[,Var1][4] <- nrow(edge[edge$Weight<0,])
  net_summary[,Var1][5] <- net_summary[3,]/net_summary[4,]
  net_summary[,Var1][6] <- table(node$fit_class)[[1]]/nrow(node)*100 # this is just the core
  net_summary[,Var1][7] <- table(node$fit_class)[[2]]/nrow(node)*100
  net_summary[,Var1][8] <- table(node$fit_class)[[3]]/nrow(node)*100
  if (is.null(neutral_fit_b)){
    neutral_fit_b[[2]] -> neutral_fit
  }else{
    rbind(neutral_fit_f[[2]], neutral_fit_b[[2]]) -> neutral_fit
  }
  net_summary[,Var1][9] <- table(neutral_fit$fit_class)[[1]]/nrow(neutral_fit)*100 # this is for all taxa
  net_summary[,Var1][10] <- table(neutral_fit$fit_class)[[2]]/nrow(neutral_fit)*100
  net_summary[,Var1][11] <- table(neutral_fit$fit_class)[[3]]/nrow(neutral_fit)*100
  net_summary[,Var1][12] <- transitivity(network_new, type = "average") # clustering coefficient
  
  net_summary[,Var1][13] <- modularity(network_new, node_new$Modules) # modularity with no weights
 
   net_summary[,Var1][14] <- length(unique(node_new$Modules))
  net_summary[,Var1][15] <- mean(table(node_new$Modules))
  net_summary[,Var1][16] <- sd(table(node_new$Modules))
  net_summary[,Var1][17] <- mean(node_new$Degree)
  net_summary[,Var1][18] <- sd(node_new$Degree)
  net_summary[,Var1][19] <- mean(node_new$Betweenness)
  net_summary[,Var1][20] <- sd(node_new$Betweenness)
  net_summary[,Var1][21] <- mean_distance(network_new) # Path Length
  net_summary[,Var1][22] <- edge_density(network_new, loops = FALSE) # Density
  # removing weights
  network_no_weight <-
    delete_edge_attr(network_new, "weight") # delete weights
  net_summary[,Var1][23] <- diameter(  network_no_weight, directed = TRUE, unconnected = TRUE) #diameter
  net_summary[,Var1][24] <- length(rownames(subset(node_new, Key%in%"Module hubs")))
  net_summary[,Var1][25] <- length(rownames(subset(node_new, Key%in%"Network hubs")))
  net_summary[,Var1][26] <- length(rownames(subset(node_new, Key%in%"Connectors")))
  net_summary[,Var1][27] <- length(rownames(subset(node_new, Key%in%"Peripherals")))
  net_summary[,Var1][28] <- length(rownames(subset(edge[edge$InterKingTaxa=="Fungi-Bacteria" | 
                                                         edge$InterKingTaxa=="Archaea-Bacteria"|
                                                         edge$InterKingTaxa=="Fungi-Archaea", ], Weight>0)))
  net_summary[,Var1][29] <- length(rownames(subset(edge[edge$InterKingTaxa=="Fungi-Bacteria" | 
                                                         edge$InterKingTaxa=="Archaea-Bacteria"|
                                                         edge$InterKingTaxa=="Fungi-Archaea", ], Weight<0)))
  # reduce decimals
  net_summary <-
    apply(net_summary, 2, function(x) round(x, digits = 3))
  dim(node_new[-c(unique(node_new$Modules)), ]) %T>% print()
  return(net_summary)
}


# EXREACT NEUTRAL OTUS ---------------------------------------------------------------------
# StatsNeutral 
StatsNeutral <- function(node, network, neutral_fit_f, neutral_fit_b, Var1, Var2){
  require(igraph)
  # creating dataframe
  network_new <- 
    igraph::delete.vertices(network, which(degree(network)==0)) # delete single nodes
  CalcNetAttr(network_new, Var1, Var2) -> node_new
  vcount(network_new) %T>% print()
  dim(node_new) %T>% print()
  net_summary <- data.frame(matrix(ncol=1, nrow =6))
  colnames(net_summary) <- "Value"
  rownames(net_summary) <- c("Above_fun", "Neutral_fun", "Below_fun",
                             "Above_prok", "Neutral_prok", "Below_prok")
  # adding metrics
  if (is.null(neutral_fit_b)){
    neutral_fit_b[[2]] -> neutral_fit
  }else{
    rbind(neutral_fit_f[[2]], neutral_fit_b[[2]]) -> neutral_fit
  }
  neutral_fit_fun <-
    subset(neutral_fit, Kingdom%in%c("Fungi") & fill%in%c("core"))
  net_summary[,"Value"][1] <- table(neutral_fit_fun$fit_class)[[1]]/nrow(neutral_fit_fun)*100 # this is for all taxa
  net_summary[,"Value"][2] <- table(neutral_fit_fun$fit_class)[[2]]/nrow(neutral_fit_fun)*100
  net_summary[,"Value"][3] <- table(neutral_fit_fun$fit_class)[[3]]/nrow(neutral_fit_fun)*100
  neutral_fit_prok <-
    subset(neutral_fit, Kingdom%in%c("Bacteria","Archaea") & fill%in%c("core"))
  neutral_fit_prok %T>% print()
  net_summary[,"Value"][4] <- table(neutral_fit_prok$fit_class)[[1]]/nrow(neutral_fit_prok)*100 # this is for all taxa
  net_summary[,"Value"][5] <- table(neutral_fit_prok$fit_class)[[2]]/nrow(neutral_fit_prok)*100
  net_summary[,"Value"][6] <- table(neutral_fit_prok$fit_class)[[3]]/nrow(neutral_fit_prok)*100
  net_summary$Kingdom <- c(rep("Fungi", 3), rep("Prokaryote", 3))
  net_summary$Depth <- rep(Var1, 6)
  net_summary$Crop <- rep(Var2, 6)
  return(net_summary)
}



# StatsNeutralFun -  

# Example StatsNeutralFun(nfit_fungi_Pop_10, "0-10cm","Poplar")


StatsNeutralFun <- function(neutral_fit, Var1, Var2){
   neutral_fit[[2]] -> neutral_fit
  # calculate
  neutral_fit_fun <-
    subset(neutral_fit, Kingdom%in%c("Fungi"))
  dim(neutral_fit_fun) %T>% print()
  net_summary <-
    melt(table(neutral_fit_fun$fit_class)/nrow(neutral_fit_fun)*100)
  net_summary$Kingdom <- rep("Fungi", nrow(net_summary))
  net_summary$Depth <- rep(Var1, nrow(net_summary))
  net_summary$Crop <- rep(Var2, nrow(net_summary))
  net_summary <- as.data.frame(net_summary)
  return(net_summary)
}


StatsNeutralProk <- function(neutral_fit, Var1, Var2){
  neutral_fit[[2]] -> neutral_fit
  # calculate
  neutral_fit_fun <-
    subset(neutral_fit, Kingdom%in%c("Bacteria","Archaea"))
  dim(neutral_fit_fun) %T>% print()
  net_summary <-
    melt(table(neutral_fit_fun$fit_class)/nrow(neutral_fit_fun)*100)
  net_summary$Kingdom <- rep("Prokaryote", nrow(net_summary))
  net_summary$Depth <- rep(Var1, nrow(net_summary))
  net_summary$Crop <- rep(Var2, nrow(net_summary))
  net_summary <- as.data.frame(net_summary)
  return(net_summary)
}



StatsNeutralCoreFun <- function(neutral_fit, Var1, Var2){
  neutral_fit[[2]] -> neutral_fit
  # calculate
  neutral_fit_fun <-
    subset(neutral_fit, Kingdom%in%c("Fungi") & fill%in%c("core"))
  dim(neutral_fit_fun) %T>% print()
  net_summary <-
    melt(table(neutral_fit_fun$fit_class)/nrow(neutral_fit_fun)*100)
  net_summary$Kingdom <- rep("Fungi", nrow(net_summary))
  net_summary$Depth <- rep(Var1, nrow(net_summary))
  net_summary$Crop <- rep(Var2, nrow(net_summary))
  net_summary <- as.data.frame(net_summary)
  return(net_summary)
}


StatsNeutralCoreProk <- function(neutral_fit, Var1, Var2){
  neutral_fit[[2]] -> neutral_fit
  # calculate
  neutral_fit_fun <-
    subset(neutral_fit, Kingdom%in%c("Bacteria","Archaea")  & fill%in%c("core"))
  dim(neutral_fit_fun) %T>% print()
  net_summary <-
    melt(table(neutral_fit_fun$fit_class)/nrow(neutral_fit_fun)*100)
  net_summary$Kingdom <- rep("Prokaryote", nrow(net_summary))
  net_summary$Depth <- rep(Var1, nrow(net_summary))
  net_summary$Crop <- rep(Var2, nrow(net_summary))
  net_summary <- as.data.frame(net_summary)
  return(net_summary)
}










# Network attributes -----------------------------------------------------------------

# Example: NetStats(nodes_Pop_10cm, edges_Pop_10cm, network_Pop_10, "0-10cm")

# This fucntion will calculate several netwrok assciated attribuets could be 
# interesting to explore.


NetStats <- function(node, edge, network, Var){
  require(igraph)
  # creating dataframe
  net_summary <- data.frame(matrix(ncol=1, nrow =35))
  colnames(net_summary) <- Var
  rownames(net_summary) <- c("Noes", "Edges", "Positive_Edges", "Negative_Edges",
                             "Ratio",
                             "Transitivity","Modularity","Modules", "Average_Module_size", 
                             "Module_SD", "Average_Degree","Degree_SD", "Average_Betweenes", 
                             "Betweenness_SD", "Path_length", "Edge_density", "Network diameter",
                             "Module_hubs", "Network_hubs", "Conncetors", "Peripherals",
                             "Fungal_OTUs","Bacterial_OTUs","Archaeal_OTUs",
                             "Interkind_positive","Interkind_negative", 
                             "Root_indicators","Soil_indicators", 
                             "Pop_ind", "Pra_ind", "Switch_ind",
                             "10cm_ind","25cm_ind","50cm_ind","100cm_ind")
  # adding metrics
  net_summary[,Var][1] <- length(rownames(node))
  net_summary[,Var][2] <- length(rownames(edge))
  net_summary[,Var][3] <- nrow(edge[edge$Weight>0,])
  net_summary[,Var][4] <- nrow(edge[edge$Weight<0,])
  net_summary[,Var][5] <- net_summary[3,]/net_summary[4,]
  net_summary[,Var][6] <- transitivity(network, type = "average") # clustering coefficient
  net_summary[,Var][7] <- modularity(network, node$Modules) # modularity with no weights
  net_summary[,Var][8] <- length(unique(node$Modules))
  net_summary[,Var][9] <- mean(table(node$Modules))
  net_summary[,Var][10] <- sd(table(node$Modules))
  net_summary[,Var][11] <- mean(node$Degree)
  net_summary[,Var][12] <- sd(node$Degree)
  net_summary[,Var][13] <- mean(node$Betweenness)
  net_summary[,Var][14] <- sd(node$Betweenness)
  net_summary[,Var][15] <- mean_distance(network) # Path Length
  net_summary[,Var][16] <- edge_density(network, loops = FALSE) # Density
  # removing weights
  network_new <- delete_edge_attr(network, "weight") # delete weights
  net_summary[,Var][17] <- diameter(network_new, directed = TRUE, unconnected = TRUE) #diameter
  net_summary[,Var][18] <- length(rownames(subset(node, Key%in%"Module hubs")))
  net_summary[,Var][19] <- length(rownames(subset(node, Key%in%"Network hubs")))
  net_summary[,Var][20] <- length(rownames(subset(node, Key%in%"Connectors")))
  net_summary[,Var][21] <- length(rownames(subset(node, Key%in%"Peripherals")))
  net_summary[,Var][22] <- length(rownames(node[node$Kingdom=="Fungi",]))
  net_summary[,Var][23] <- length(rownames(node[node$Kingdom=="Bacteria",]))
  net_summary[,Var][24] <- length(rownames(node[node$Kingdom=="Archaea",]))
  net_summary[,Var][25] <- length(rownames(subset(edge[edge$InterKingTaxa=="Fungi-Bacteria" | 
                                                         edge$InterKingTaxa=="Archaea-Bacteria"|
                                                         edge$InterKingTaxa=="Fungi-Archaea", ], Weight>0)))
  net_summary[,Var][26] <- length(rownames(subset(edge[edge$InterKingTaxa=="Fungi-Bacteria" | 
                                                         edge$InterKingTaxa=="Archaea-Bacteria"|
                                                         edge$InterKingTaxa=="Fungi-Archaea", ], Weight<0)))
  net_summary[,Var][27] <- length(rownames(subset(node, Niche_ind%in%"Root")))
  net_summary[,Var][28] <- length(rownames(subset(node, Niche_ind%in%"Soil")))
  net_summary[,Var][29] <- length(rownames(subset(node, Crop_ind%in%"Poplar")))
  net_summary[,Var][30] <- length(rownames(subset(node, Crop_ind%in%"Prairie")))
  net_summary[,Var][31] <- length(rownames(subset(node, Crop_ind%in%"Switchgrass")))
  net_summary[,Var][32] <- length(rownames(subset(node, Depth_ind%in%"0-10cm")))
  net_summary[,Var][33] <- length(rownames(subset(node, Depth_ind%in%"10-25cm")))
  net_summary[,Var][34] <- length(rownames(subset(node, Depth_ind%in%"25-50cm")))
  net_summary[,Var][35] <- length(rownames(subset(node, Depth_ind%in%"50-100cm")))
  # reduce decimals
  net_summary <-
    apply(net_summary, 2, function(x) round(x, digits = 3))
  return(net_summary)
}





# PLOTTING THE NETWORK -----------------------------------------------------------

#Example: MakeTabGraph(nodes_Pop_10cm, edges_Pop_10cm)

# Make a table-rgaph with ggraph

MakeTabGraph <- function(nodes_df, edges_df){
  graph <- as_tbl_graph(edges_df, node_key = "OTU_ID", directed = FALSE)
  #graph %T>% print()
  graph_new <- graph %>% 
    activate(nodes) %>% 
    dplyr::left_join(nodes_df, by = c(name = "OTU_ID")) %>% 
    dplyr::rename(OTU_ID = name) %>% 
    full_join(as.data.frame(nodes_df[nodes_df$Degree==0,])) # to include disconnected nodes
  return(graph_new)
}





















