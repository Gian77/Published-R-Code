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

# >>>>>>>>> PHYLOGENETIC TREE <<<<<<<<<< ---------------------------------------

# Load R packages. Define your package list as a character vector
pkgs <- c(
  # --- Project & Package Management ---
  "renv",            # Reproducible environments (lockfile management)
  "pak",             # Fast, modern package installation dependency engine
  "cli",             # Attractive and informative command-line interfaces  
  "styler",          # Automated R code formatting
  "janitor",         # Data cleaning
  "magrittr",        # Pipe operators (%>%)
  
  # --- Sequence Analysis & Phylogenetics ---
  "Biostrings",      # DNA/RNA/AA sequence containers
  "ape",             # Core tree handling
  "msa",             # Multiple Sequence Alignment
  "DECIPHER",        # Sequence alignment & chimera detection
  "phangorn",        # Phylogenetic analysis (MP, ML)
  "tidysq",          # Tidy processing of sequences
  "tidytree",        # Tidy manipulation of trees
  "ggtree",          # Visualization of trees
  
  # --- Microbial Ecology & Diversity ---
  "decontam",        # Contaminant OTU removal
  "phyloseq",        # Integration of OTU tables & taxonomy
  "speedyseq",       # High-performance phyloseq
  "vegan",           # Community ecology (Ordination)
  "AICcPermanova",   # Model selection for PERMANOVA
  
  # --- Data Science & Visualization ---
  "tidyverse",       # Core suite (ggplot2, dplyr, etc.)
  "ggtext",          # Markdown labels
  "ggpubr",          # Publication-ready themes
  "cowplot",         # Plot arrangement
  "gridExtra",       # Miscellaneous graphics
  "ggrepel",         # Non-overlapping labels
  "scales",          # Axis transformations
  "ggh4x",           # Extended facets
  "ggstar",          # Geometric shapes
  "ggtreeExtra",     # Annotation visualization
  "ggnewscale",      # Multiple scales
  
  # --- Statistics & Modeling ---
  "agricolae",       # Agricultural stats
  "BRCore",          # Internal core functions
  "lme4",            # Mixed-Effects Models (GLMM)
  "glmmTMB",         # Fast GLMMs
  "robustlmm",       # Outlier resistance
  "DHARMa",          # Residual diagnostics
  "parallel",        # Multi-core support
  "ggeffects",       # Marginal effects
  "sjPlot",          # Visualizing models
  "broom.mixed",     # Tidy summaries
  "merTools",        # Analyzing lme4 objects
  "multcompView",    # Significance letters
  "gllvm",           # Latent Variable Models
  "Maaslin2",        # Microbiome stats
  
  # --- Networks & Hierarchical Viz ---
  "ggdendro",        # Dendrogram data
  "igraph"            # Network visualization
)

# Load them all silently
invisible(lapply(pkgs, library, character.only = TRUE))

# Tracking package versions with renv ------------------------------------------
# renv::init()      # initializes renv in your project
# renv::restore()   # installs all packages from the lockfile
renv::snapshot()   # updates the lockfile
renv::status()
renv::install("")

# Then commit and push the updated renv.lock file!

# NOTE. Positron options to restore Rstudio projects load(".Rdata")

# R options --------------------------------------------------------------------
options(scipen = 9999, pillar.sigfig = 6, digits = 6, max.print = 10000000)
#rm(list = ls())

# Check the lib paths ----------------------------------------------------------
.libPaths()

# Session Info -----------------------------------------------------------------
sessionInfo()

# **********************************************************************--------
# ***** PATHS ***** ------------------------------------------------------------

data_path <-
  ("/home/gian/Dropbox/6_PROJETCS/Published-R-Code/Liu_etal_2026_AMF_FertUnfertSwitchgrass")

data_path

# **********************************************************************--------
# **** 4. PHYLOGENETIC TREE **** ------------------------------------------

# https://yulab-smu.top/treedata-book/chapter10.html

palette_bestmatch <-
  c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77",
    "#283dff","#636bb7","#bfc5ff","#014443","#195637",
    "#117744","#60ffaf","#b7ffdb","#825121","#ea7f17",
    "#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d",
    "#5b5b19","#fcfc00","#ffff9e","#ffb7ef","#fa7efc",
    "#ae09ea","#521899","#1e0047")

# A. Generating a 99% OTUs tree ---------------------------------------------------
#Working with the RAXML tree for now
tree_raxml_AMF <- phy_tree(physeq_AMF_rare)

# Generating metadata for plotting ---------------------------------------------
melted_AMF <- 
  physeq_AMF_rare %>%
  psmelt() %>%
  arrange(Genus)

head(melted_AMF)
dim(melted_AMF)

dat1 <- melted_AMF %>%
  group_by(OTU, BestMatch, Genus, Phylum, Family) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus)) %>% 
  rename(label = OTU)

range(dat1$MeanAbundance)
range(sqrt(dat1$MeanAbundance))

dat2 <- melted_AMF %>%
  group_by(OTU, site) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  rename(ID = OTU, Sites = site)

dat3 <- melted_AMF %>%
  group_by(OTU, site) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  rename(ID = OTU, Sites = site)

# Build original plot ---------------------------------------------------------------
tree_AMF <- 
  tree_raxml_AMF %>% 
  ggtree(layout = "circular", branch.length = "none", size = 0.1, open.angle = 5)

tree_AMF <-
  tree_raxml_AMF %>% 
  ggtree(branch.length = "none", size = 0.15, open.angle = 5)

tree_AMF <- 
  tree_raxml_AMF %>% 
  ggtree(layout = "fan", size = 0.15, open.angle = 5)


# Layer 1: tip points — color = Genus, size = mean abundance
tree_AMF <- 
  tree_AMF %<+% dat1 +
  geom_tippoint(
    aes(fill = BestMatch, size = sqrt(MeanAbundance)),
    shape = 21, stroke = 0.1, alpha = 0.85) +
  scale_fill_manual(
    values = palette_bestmatch,
    guide  = guide_legend(
      keywidth = 0.5, keyheight = 0.5, order = 1,
      override.aes = list(shape = 22, size = 3)
    ),
    na.translate = FALSE
  ) +
  scale_size_continuous(
    range = c(0.1, 4),
    name  = "sqrt(Abundance)",
    guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 2,
                         override.aes = list(shape = 21))
  )

tree_AMF <- 
  tree_AMF +
  geom_tiplab2(
    aes(label = BestMatch),
    size = 1,
    align = TRUE,
    offset = 0.02
  )

tree_AMF 

# Layer 2: site abundance heatmap ring
tree_AMF <- 
  tree_AMF +
  new_scale_fill() +
  geom_fruit(
    data    = dat2,
    geom    = geom_tile,
    mapping = aes(y = ID, x = Sites, alpha = Abundance, fill = Sites),
    color   = "grey90", offset = 0.04, size = 0.02
  ) +
  scale_alpha_continuous(
    range = c(0, 1),
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 4)
  ) +
  scale_fill_manual(
    values = palette_site,
    guide  = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 3)
  )

# Layer 3: total abundance bar ring
tree_AMF <- 
  tree_AMF +
  new_scale_fill() +
  geom_fruit(
    data        = dat3,
    geom        = geom_bar,
    mapping     = aes(y = ID, x = TotalAbundance, fill = Sites),
    pwidth      = 0.38,
    orientation = "y",
    stat        = "identity"
  ) +
  scale_fill_manual(
    values = palette_site,
    guide  = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 3)
  )

# Theme 
tree_AMF <- 
  tree_AMF +
  geom_treescale(fontsize = 2, linesize = 0.3) +
  theme(
    legend.position   = c(0.93, 0.5),
    legend.background = element_rect(fill = NA),
    legend.title      = element_text(size = 6.5),
    legend.text       = element_text(size = 4.5),
    legend.spacing.y  = unit(0.02, "cm")
  )

tree_AMF

tree_AMF + layout_circular() + 
  theme(legend.position=c(0.05, 0.7)) 

tree_AMF + layout_rectangular() + 
  theme(legend.position=c(0.7, 0.2)) 

# B. Reducing tree complexity my agglomerating taxa ----------------------------

sample_data(physeq_AMF_rare)$fert_status_site <- 
  paste(sample_data(physeq_AMF_rare)$fert_status, 
        sample_data(physeq_AMF_rare)$site,
        sep="-")

physeq_AMF_rare@sam_data

merged_AMF <- merge_samples(physeq_AMF_rare, "fert_status_site")
merged_AMF@sam_data

merged_AMF <- 
  tax_glom(merged_AMF, taxrank = "BestMatch") %>% 
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.)

# Checks
merged_AMF@tax_table %>% as.data.frame()
merged_AMF@tax_table %>% as.data.frame() %>% subset(Genus == "Diversispora")
merged_AMF

# Melting 
melt_merged_AMF <- 
  psmelt(merged_AMF) %>% 
  separate(col = Sample, into = c("fert_status", "site"), sep = "-") %>% 
  mutate(val = sqrt(Abundance))

#melt_merged_AMF$Species
#separate(col = Sample, into = c("fert_status", "site"), sep = "-")
melt_merged_AMF %>% subset(Genus == "Diversispora")
head(melt_merged_AMF)
dim(melt_merged_AMF)

# NOTE. There are different Zotus classified with the same BestMatch name.
tree_raxml_AMF_merged <- phy_tree(merged_AMF)
tree_raxml_AMF_merged

# Generating metadata for plotting ---------------------------------------------
dat1_merged <- 
  melt_merged_AMF %>%
  group_by(OTU, BestMatch, Genus, Phylum, Family) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus)) %>% 
  rename(label = OTU)

range(dat1_merged$MeanAbundance)
range(sqrt(dat1_merged$MeanAbundance))

dat2_merged <-
  melt_merged_AMF %>%
  group_by(OTU, site) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  rename(ID = OTU, Sites = site)

dat2_merged

dat3_merged <- 
  melt_merged_AMF %>%
  group_by(OTU, site) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  rename(ID = OTU, Sites = site)

dat3_merged

dat4_merged <- 
  melt_merged_AMF %>%
  group_by(OTU, fert_status) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  rename(ID = OTU, N_addition = fert_status)

dat4_merged


# 1) Build plot version 1 ---------------------------------------------------------
tree_AMF_merged <- 
  tree_raxml_AMF_merged %>% 
  ggtree(layout = "circular", branch.length = "none", size = 0.15, open.angle = 5)

tree_AMF_merged <- 
  tree_raxml_AMF_merged %>% 
  ggtree(layout = "fan", size = 0.15, open.angle = 5)

tree_AMF_merged

# Layer 1: tip points — color = Genus, size = mean abundance
tree_AMF_merged <- 
  tree_AMF_merged %<+% dat1_merged +
  geom_tippoint(
    aes(fill = BestMatch, size = sqrt(MeanAbundance)),
    shape = 21, stroke = 0.1, alpha = 0.85) +
  scale_fill_manual(
    values = palette_bestmatch,
    guide  = guide_legend(
      keywidth = 0.5, keyheight = 0.5, order = 1,
      override.aes = list(shape = 21, size = 3)
    ),
    na.translate = FALSE
  ) +
  scale_size_continuous(
    range = c(0.1, 4),
    name  = "sqrt(Abundance)",
    guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 2,
                         override.aes = list(shape = 21))
  )

tree_AMF_merged <- 
  tree_AMF_merged +
  geom_tiplab2(
    aes(label = BestMatch),
    size = 3,
    align = TRUE,
    offset = 0.2
  )

tree_AMF_merged

# Layer 2: site abundance heatmap ring
tree_AMF_merged <- 
  tree_AMF_merged +
  new_scale_fill() +
  geom_fruit(
    data    = dat2_merged,
    geom    = geom_tile,
    mapping = aes(y = ID, x = Sites, alpha = Abundance, fill = Sites),
    color   = "grey90", offset = 0.04, size = 0.02
  ) +
  scale_alpha_continuous(
    range = c(0, 1),
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 4)
  ) +
  scale_fill_manual(
    values = palette_site,
    guide  = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 3)
  )

# Layer 3: total abundance for fert_status
tree_AMF_merged <- 
  tree_AMF_merged +
  new_scale_fill() +
  geom_fruit(
    data    = dat4_merged,
    geom    = geom_tile,
    mapping = aes(y = ID, x = N_addition, alpha = Abundance, fill = N_addition),
    color   = "grey90", offset = 0.04, size = 0.02
  ) +
  scale_alpha_continuous(
    range = c(0, 1),
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 4)
  ) +
  scale_fill_manual(
    values = palette_site,
    guide  = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 3)
  )

# Layer 4: total abundance for site
tree_AMF_merged <- 
  tree_AMF_merged +
  new_scale_fill() +
  geom_fruit(
    data        = dat3_merged,
    geom        = geom_bar,
    mapping     = aes(y = ID, x = TotalAbundance, fill = Sites),
    pwidth      = 0.38,
    orientation = "y",
    stat        = "identity"
  ) +
  scale_fill_manual(
    values = palette_site,
    guide  = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 3)
  )

# Theme 
tree_AMF_merged <- 
  tree_AMF_merged +
  geom_treescale(fontsize = 2, linesize = 0.3) +
  theme(
    legend.position   = c(0.93, 0.5),
    legend.background = element_rect(fill = NA),
    legend.title      = element_text(size = 6.5),
    legend.text       = element_text(size = 4.5),
    legend.spacing.y  = unit(0.02, "cm")
  )

tree_AMF_merged


tree_AMF_merged <- tree_AMF_merged + 
  geom_tiplab(
    aes(label = BestMatch), 
    size = 2.5,           # Keep it small for 1,700 tips!
    offset = 0.5,         # Push it slightly past the tippoint
    linetype = "dotted",   # "33" or "dotted"
    linewidth = 0.2,       # Keep it subtle
    color = "grey50", 
    geom = "text",        # Use 'text' to avoid the dotted line logic here
    align = TRUE         # Set to FALSE so it stays next to the point
  )


tree_AMF_merged

tree_AMF_merged + layout_circular() + 
  theme(legend.position=c(0.05, 0.7)) 

tree_AMF_merged + layout_rectangular() + 
  theme(legend.position=c(0.05, 0.7)) 

# Imporvements -----------------------------------------------------------------
# https://www.earlham.ac.uk/articles/plotting-phylogenetic-trees-r-alternating-clade-highlights

tree_AMF <- ggtree(tree_raxml_AMF) %<+% dat1

genera_nodes <- 
  dat1 %>%
  group_by(BestMatch) %>%
  filter(n() > 1) %>% # Only collapse if there's more than one tip
  summarize(node = ggtree::MRCA(tree_raxml_AMF, label))

genera_nodes

for(n in genera_nodes$node) {
  tree_AMF <- collapse(tree_AMF, node = n, mode = 'mixed', alpha = 0.4)
}

tree_AMF + 
  geom_tippoint(aes(fill = BestMatch, size = sqrt(MeanAbundance)), 
                shape = 21, stroke = 0.1) +
  geom_tiplab(aes(subset = !is.na(Genus)), offset = 0.2) + # Adjust as needed
  scale_fill_manual(values = palette_bestmatch)

tree_AMF + 
  geom_tippoint(aes(fill = BestMatch, size = sqrt(MeanAbundance)), 
                shape = 21, stroke = 0.1) +
  geom_tiplab(aes(subset = !is.na(Genus)), offset = 0.2) + # Adjust as needed
  scale_fill_manual(values = palette_bestmatch)

# This is Goood! 
# The Better Approach: geom_strip()

genera_nodes <- 
  dat1 %>%
  filter(label %in% tree_raxml_AMF$tip.label) %>%
  group_by(BestMatch) %>%
  summarize(
    t1 = first(label),
    t2 = last(label),
    .groups = "drop"
  ) %>%
  filter(t1 != t2) # Only draw strips for groups with more than 1 tip

genera_nodes

tree_AMF <- 
  ggtree(tree_raxml_AMF) %<+% dat1 +
  geom_tippoint(aes(fill = BestMatch, size = sqrt(MeanAbundance)), 
                shape = 21, stroke = 0.1, alpha = 0.8) +
  scale_fill_manual(values = palette_bestmatch)

for(i in 1:nrow(genera_nodes)) {
  tree_AMF <- 
    tree_AMF + 
    geom_strip(
      taxa1 = genera_nodes$t1[i], 
      taxa2 = genera_nodes$t2[i], 
      label = genera_nodes$BestMatch[i],
      offset = 0.1,       # Adjust this value based on your tree's branch lengths
      barsize = 1.5, 
      extend = 0.2, 
      fontsize = 3,
      offset.text = 0.1
    )
}

tree_AMF

tree_AMF + scale_fill_manual(
  values = palette_bestmatch,
  guide = guide_legend(
    title = "Taxa",
    ncol = 1,               # Force single column
    byrow = TRUE,
    override.aes = list(
      shape = 22,           # 22 is a square with a border
      size = 4,             # Smaller square size
      stroke = 0.2          # Thin border for the square
    )
  ),
  na.translate = FALSE
) +
  theme(
    legend.key.height = unit(0.4, "cm"), # Reduces vertical space between items
    legend.text = element_text(size = 8), # Smaller font size
    legend.title = element_text(size = 9, face = "bold")
  ) +
  xlim(0, 4) 

# Ttsting another option

tree_AMF <- 
  ggtree(tree_raxml_AMF, linewidth=0.5) %<+% dat1 +
  geom_tiplab(aes(label=label), size=2)

tree_AMF

#Make dataframe for clade nodes
genera_nodes <- data.frame(
  clade=unique(dat1$BestMatch),
  node=NA
)

genera_nodes

#Find the most recent common ancestor for each clade
for (i in 1:length(genera_nodes$clade)) {
  
  genera_nodes$node[i] <- MRCA(
    tree_raxml_AMF,
    dat1$label[dat1$BestMatch == genera_nodes$clade[i]]
  )
  
}

tree_raxml_AMF

tree_AMF <- 
  ggtree(tree_raxml_AMF, linetype=NA) %<+% dat1 +
  geom_highlight(data=genera_nodes, 
                 aes(node=node, fill=clade),
                 alpha=1,
                 align="right",
                 extend=0.1,
                 show.legend=FALSE) +
  geom_tree(linewidth=0.5) +
  geom_tiplab(aes(label=BestMatch), size=2)

tree_AMF




# 2) Build plot version 2 ------------------------------------------------------

tree_AMF_merged <- 
  merged_AMF %>% 
  ggtree(layout="fan", open.angle=1) + 
  geom_tippoint(mapping=aes(color = BestMatch, size = Abundance),
                show.legend=FALSE) +
  scale_color_manual(
    values = palette_bestmatch)

tree_AMF_merged

#tree_AMF_merged <- rotate_tree(tree_AMF_merged, -90)

tree_AMF_merged <- 
  tree_AMF_merged +
  geom_fruit(
    data=melt_merged_AMF %>% dplyr::select(OTU, val),
    geom=geom_boxplot,
    mapping = aes(
      y=OTU,
      x=val,
      group=OTU,
      fill=BestMatch,
    ),
    size=.2,
    outlier.size=0.5,
    outlier.stroke=0.08,
    outlier.shape=21,
    axis.params=list(
      axis       = "x",
      text.size  = 1.8,
      hjust      = 1,
      vjust      = 0.5,
      nbreak     = 3,
    ),
    grid.params=list()
  ) 


tree_AMF_merged <- 
  tree_AMF_merged +
  scale_fill_manual(
    values = palette_bestmatch,
    guide  = guide_legend(
      keywidth = 0.5, 
      keyheight = 0.5, 
      ncol = 1,
      override.aes = list(shape = 21, size = 6)
    ),
    na.translate = FALSE
  ) +
  theme(
    legend.title=element_text(size=9), 
    legend.text=element_text(size=7) 
  )


tree_AMF_merged <- 
  tree_AMF_merged + 
  geom_tiplab(
    aes(label = BestMatch), 
    size = 2.5,           # Keep it small for 1,700 tips!
    offset = 0.5,         # Push it slightly past the tippoint
    linetype = "dotted",   # "33" or "dotted"
    linewidth = 0.2,       # Keep it subtle
    color = "grey50", 
    geom = "text",        # Use 'text' to avoid the dotted line logic here
    align = TRUE         # Set to FALSE so it stays next to the point
  )

tree_AMF_merged

# Adding side abudances
tree_AMF_merged <- 
  tree_AMF_merged +
  new_scale_fill() +
  geom_fruit(
    data    = dat3_merged,
    geom    = geom_tile,
    mapping = aes(y = ID, x = Sites, alpha = val, fill = Sites),
    color   = "grey90", offset = 0.04, size = 0.02
  ) +
  scale_alpha_continuous(
    range = c(0, 1),
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 4)
  ) +
  scale_fill_manual(
    values = palette_site,
    guide  = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 3)
  )

tree_AMF_merged

# **********************************************************************--------
# B. Generating a 97% OTUs tree ------------------------------------------------

tree_raxml_97 <- 
  read.tree("phylogeny/otus_97_filtered_mafft_trim_spprt.raxml.support") 
str(tree_raxml_97)
ggtree::ggtree(tree_raxml_97)

# Rooting the tree to the outgroup ---------------------------------------------
tree_raxml_97_rooted <- root(tree_raxml_97, outgroup = "Zotu6503", resolve.root = TRUE)
ggtree::ggtree(tree_raxml_97_rooted)

tree_iqtree2_97 <- read.tree("phylogeny/otus_97_filtered_mafft_trim_iq2.treefile")
str(tree_iqtree2_97)
ggtree::ggtree(tree_iqtree2_97)

# Generate the phyloseq object -------------------------------------------------
# otutable
otutable_97 <- 
  read.delim(file.path(data_path, "datasets/otutab_97.txt"), row.names = 1)

# Taxonomy
taxonomy_97 <- 
  extract_blasTAX(tax_path = file.path(data_path, "datasets/taxonomy_blast_97.txt"),
                  namemap_path = file.path(data_path, "datasets/name_mapping_97.txt")) %>% 
  dplyr::select(
    "Zotu", "Kingdom" ,"Phylum", "Class","Order","Family","Genus","Species") %>%
  filter(Phylum %in% "Mucoromycota") %>% 
  FinalizeTaxonomy()

# I will use Mortierellomycetes as outgroup to root the tree?
taxonomy_97 %>% subset(Class %in% c("Mortierellomycetes", "Endogonomycetes"))

# Fixing classificaiton inconsistencies
CheckTaxonomyConsistency(taxonomy_97, 
                         return_long=TRUE)

CheckTaxonomyConsistency(taxonomy_97, 
                         return_long=FALSE)

taxonomy_97

taxonomy_97 <- 
  taxonomy_97 %>%
  mutate(Family = if_else(Family == "Archaeosporaceae", "Ambisporaceae", Family))

# OTUs
zotu_97 <- readDNAStringSet(file.path(data_path,"datasets/otus_97.fasta"), 
                            format="fasta", seek.first.rec=TRUE, use.names=TRUE)

zotu_97_filtered <- zotu_97[ !names(zotu_97) %in% zotus_to_remove ]


# Removing only Jimgerdemannia
zotu_97_filtered <- zotu_97[ !names(zotu_97) %in% "Zotu11542" ]

# Phyloseq
physeq_97 <- generate_phyloseq(otu=otutable_97,
                               metadata=metadata_99,
                               taxonomy=taxonomy_97 %>% column_to_rownames("Zotu"),
                               sequences=zotu_97_filtered) # I want to keep the outgroups
physeq_97
physeq_97@sam_data

# Adding phylotree
physeq_97 <-
  phyloseq(
    physeq_97@otu_table,
    physeq_97@sam_data,
    physeq_97@tax_table,
    physeq_97@refseq,
    phy_tree(tree_raxml_97_rooted)
  ) %>% 
  phyloseq::prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  phyloseq::prune_samples(sample_sums(x=.) > 0, x =.)

physeq_97
print(sample_data(physeq_97), n=113)
sort(sample_sums(physeq_97))

# Remove control samples and Rarefaction (single round in this case) 
physeq_97_rare <-
  physeq_97 %>%
  subset_samples(!site %in% "Control") %>% 
  rarefy_even_depth(rngseed = 260414, sample.size = 6434) %>% 
  prune_taxa(taxa_sums(x = .) > 0, x = .) %>% 
  prune_samples(sample_sums(x=.) > 0, x =.) 

physeq_97_rare
print(physeq_97_rare@sam_data, n=109)

# Working with the RAXML tree for now
tree_raxml_AMF_97 <- phy_tree(physeq_97_rare)
ggtree::ggtree(tree_raxml_AMF_97)

phyloseq::plot_tree(physeq_97_rare, "treeonly", label.tips = "Taxonomy", text.size=2) + 
  coord_polar(theta = "y") 

# Generating metadata for the phyloegentic tree --------------------------------
melted_AMF_97 <- 
  physeq_97_rare %>%
  psmelt() %>%
  arrange(Genus)

head(melted_AMF_97)
dim(melted_AMF_97)

dat1_97 <- 
  melted_AMF_97 %>%
  group_by(OTU, BestMatch, Genus, Phylum, Family) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus)) %>% 
  rename(label = OTU) %>% 
  mutate(OTU = paste(label, BestMatch, sep="-"))

range(dat1_97$MeanAbundance)
range(sqrt(dat1_97$MeanAbundance))

# Modifying the taxonomy adding clades
write.csv(x= dat1_97, file= file.path(data_path, "datasets/tree_metadata_97.csv"))
dat1_97_mod <- read.csv(file=file.path(data_path, "datasets/tree_metadata_97_final.csv"), 
                        row.names = 1, header = TRUE)

head(dat1_97_mod)
levels(as.factor(dat1_97_mod$Taxonomy))

dat2_97 <- 
  melted_AMF_97 %>% 
  group_by(OTU, site) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  rename(ID = OTU, Sites = site)

dat3_97 <- 
  melted_AMF_97 %>%
  group_by(OTU, site) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  rename(ID = OTU, Sites = site)

# Grouping to clades -----------------------------------------------------------

# Make dataframe for clade nodes
genera_nodes <- data.frame(
  clade=unique(dat1_97_mod$Taxonomy),
  node=NA
)

# Find the most recent common ancestor for each clade
for (i in 1:length(genera_nodes$clade)) {
  
  genera_nodes$node[i] <- MRCA(
    tree_raxml_AMF_97,
    dat1_97_mod$label[dat1_97_mod$Taxonomy == genera_nodes$clade[i]]
  )
  
}

tree_AMF_97 <- 
  tree_raxml_AMF_97 %>% 
  ggtree(layout = "circular", branch.length = "none", size = 0.1, open.angle = 5) %<+% dat1_97_mod


# Add column with alternating binary value. This is based on the ggtree data
genera_nodes <-
  genera_nodes[match(
    tree_AMF_97$data %>%
      filter(isTip == "TRUE") %>%
      arrange(y) %>%
      pull(Taxonomy) %>%
      unique(),
    genera_nodes$clade
  ), ] %>% 
  dplyr::mutate(highlight = rep(c(0,1),
                                length.out=length(genera_nodes$clade)))

genera_nodes

# Modifying the taxonomy adding clades
write.csv(x= genera_nodes, file= file.path(data_path, "datasets/genera_nodes_97.csv"))
genera_nodes_mod <- read.csv(file=file.path(data_path, "datasets/genera_nodes_97_final.csv"), 
                             row.names = 1, header = TRUE)

head(genera_nodes_mod)
levels(as.factor(genera_nodes_mod$Taxonomy))

# Color palette for Taxonomy --------------------------------------------------


palette_bestmatch <- c(
  # --- ARCHAEOSPORALES (Earth Greens) ---
  "Ambispora"               = "#787811",
  "Archaeosporales IV"         = "#D2D21E",
  "Archaeosporales III"         = "#D2D21E",
  "Archaeosporales II"         = "#D2D21E",
  "Archaeosporales I"         = "#D2D21E",
  "Archaeosporales" = "#D2D21E",
  "Paraglomus I"              = "#784511", 
  "Paraglomus II"              = "#784511", 
  "Paraglomus laccatum"     = "#EAAB6C",
  "Paraglomus brasilianum"  = "#E4913F",
  
  # --- DIVERSISPORALES (Teals and Seafoams) ---
  "Gigaspora"               = "#083B3B", 
  "Dentiscutata"            = "#117845",
  "Cetraspora"              = "#1ED278",
  "Scutellospora"           = "#6CEAAB",
  "Acaulospora"             = "#117878",
  "Diversispora"            = "#3FE4E4",
  
  # --- ENTROPHOSPORALES & OTHERS (Reds and Pinks) ---
  "Entrophospora"           = "#D21E2C",
  
  # Outgroup
  "Mortierella"             = "black",
  
  # Dominikia & Silvaspora Clade (Deep Blues)
  "Dominikia I"               = "#114578",
  "Dominikia II"                = "#103456",
  "Glomeraceae I"           = "#185EA5",
  "Microdominikia litorea I"  = "#1E78D2",
  "Microdominikia litorea II"  = "#1E78D2",
  
  # --- GLOMERALES ---
  # Group 1 (Warm Magenta/Pink Gradient - Vivid & Deep)
  "Rhizophagus"             = "#E43FAD", # Your original pink
  "Glomeraceae II"          = "#FF00FF", # Pure Magenta
  "Oehlia"                  = "#C71585", # Medium Violet Red
  "Kamienskia"              = "#8B008B", # Dark Magenta
  "Microkamienskia"         = "#550055", # Deep Plum
  
  # Group 2 (Cold Lavender/Purple Gradient - Blue-ish Purples)
  "Glomus"               = "#4B0082", # Indigo
  "Complexispora"           = "#6A5ACD", # Slate Blue
  "Septoglomus"             = "#9370DB", # Medium Purple
  "Funneliformis"           = "#B19CD9"  # Light Pastel Purple
)

# Build tree plot ---------------------------------------------------------------

tree_AMF_97 <- 
  tree_raxml_AMF_97 %>% 
  ggtree(layout = "circular", branch.length = "none", size = 0.1, open.angle = 5) %<+% dat1_97_mod 

#max_depth <- max(ggtree(tree_raxml_AMF_97)$data$x, na.rm = TRUE)
max_depth <- max(tree_AMF_97$data$x, na.rm = TRUE)

# coerce to numeric; empty strings become NA
tree_AMF_97$data$bootstrap <- suppressWarnings(
  as.numeric(tree_AMF_97$data$label)
)

tree_AMF_97 +
  geom_tippoint(
    aes(fill = Taxonomy, size = sqrt(MeanAbundance)),
    shape = 21, stroke = 0.1, alpha = 0.85) +
  scale_size_continuous(range = c(0.331801/5, 22.456278/5),
                        name = "Abundance") +
  geom_cladelab(data = genera_nodes_mod %>% filter(clade != "Mortierella"),
                mapping = aes(node = node, label = clade, color = clade),
                parse = TRUE,
                fontsize = 5,
                align = TRUE,                          # all bars on one ring
                angle = "auto", 
                #horizontal = FALSE, 
                offset = max_depth * 0.02,             # small offset
                offset.text = max_depth * 0.02,
                barsize = 1.5,
                show.legend = FALSE) +
  geom_text2(aes(subset = !isTip, label = node), size = 2, color = "black")

tree_AMF_97 +
  geom_tippoint(
    aes(fill = Taxonomy, size = sqrt(MeanAbundance)),
    shape = 21, stroke = 0.1, alpha = 0.85) +
  scale_size_continuous(range = c(0.331801/5, 22.456278/5),
                        name = "Abundance") +
  geom_cladelab(data = genera_nodes_mod %>% filter(clade != "Mortierella"),
                mapping = aes(node = node, label = clade, color = clade),
                parse = TRUE,
                fontsize = 5,
                align = TRUE,                          # all bars on one ring
                angle = "auto", 
                #horizontal = FALSE, 
                offset = max_depth * 0.02,             # small offset
                offset.text = max_depth * 0.02,
                barsize = 1.5,
                show.legend = FALSE) +
  xlim(NA, max_depth * 1.2) +                          # room for labels, no more
  guides(fill = guide_legend(ncol = 1, keywidth = 0.2, keyheight = 0.2, order = 1,
                             override.aes = list(shape = 22, size = 4))) +
  geom_text2(aes(subset = !isTip & !is.na(bootstrap) & bootstrap >= 70, label = label),
             size = 2, hjust = -0.2, color = "grey30") +
  scale_fill_manual(values = palette_bestmatch, na.translate = FALSE) +
  scale_color_manual(values = palette_bestmatch, na.translate = FALSE) +
  theme(legend.position = "none")

new_scale_fill() +geom_fruit(
  axis.params = list(title ="Abundance in site"),
  data    = dat2_97,
  geom    = geom_tile,
  mapping = aes(y = ID, x = Sites, alpha = Abundance, fill = Sites),
  color   = "grey90", 
  offset = 0.5, 
  size = 0.02) +
  scale_alpha_continuous(
    range = c(0.3, 22.5), 
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 4)
  ) +
  scale_fill_manual(
    values = palette_site,
    guide  = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 3)
  )




# Sites factor levels — already set in your data, but make sure dat3 matches
dat3_97$Sites <- factor(dat3_97$Sites, levels = levels(dat2_97$Sites))

# Site palette — pick 5 distinct colors for your 5 sites
site_palette <- c("Lux Arbor"   = "#0000FF",
                  "Lake City"   = "#FFA500",
                  "Escanaba"    = "#FF0000",
                  "Rhinelander" = "#006400",
                  "Hancock"     = "#800080")

max_depth <- max(tree_AMF_97$data$x, na.rm = TRUE)

p <- tree_AMF_97 +
  # tip points (your existing layer, kept as-is)
  geom_tippoint(
    aes(fill = Taxonomy, size = sqrt(MeanAbundance)),
    shape = 21, stroke = 0.1, alpha = 0.85) +
  scale_size_continuous(range = c(0.331801/5, 22.456278/5), name = "Abundance") +
  scale_fill_manual(values = palette_bestmatch, na.translate = FALSE,
                    guide = guide_legend(ncol = 1, keywidth = 0.3, keyheight = 0.3,
                                         order = 1, override.aes = list(shape = 22, size = 4))) +
  scale_color_manual(values = palette_bestmatch, na.translate = FALSE) +
  
  # genus clade labels
  geom_cladelab(data = genera_nodes %>% filter(clade != "Mortierella"),
                mapping = aes(node = node, label = clade, color = clade),
                parse = TRUE, fontsize = 2.5, align = TRUE, angle = "auto",
                offset = max_depth * 0.02, offset.text = max_depth * 0.02,
                barsize = 0.4, show.legend = FALSE) +
  
  # bootstrap labels
  geom_text2(aes(subset = !isTip & !is.na(bootstrap) & bootstrap >= 70,
                 label = label),
             size = 2, hjust = -0.2, color = "grey30") +
  
  # Ring 1: per-site abundance heatmap
  new_scale_fill() +
  geom_fruit(
    data    = dat2_97,
    geom    = geom_tile,
    mapping = aes(y = ID, x = Sites, alpha = Abundance, fill = Sites),
    color   = "grey90", 
    offset = 0.3, 
    size = 0.02) +
  scale_alpha_continuous(
    range = c(0, 22.5), 
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 4)
  )  +
  
  # Ring 2: per-site total abundance bars
  geom_fruit(
    data = dat3_97,
    geom = geom_bar,
    mapping = aes(y = ID, x = TotalAbundance, fill = Sites),
    pwidth = 0.38,
    orientation = "y",
    stat = "identity"
  ) +
  scale_fill_manual(values = palette_site,
                    guide = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 4)) +
  
  # tree scale bar
  geom_treescale(fontsize = 2, linesize = 0.3) +
  
  xlim(NA, max_depth * 2.2) +
  
  theme(legend.position = c(0.93, 0.5),
        legend.background = element_rect(fill = NA),
        legend.title = element_text(size = 6.5),
        legend.text = element_text(size = 4.5),
        legend.spacing.y = unit(0.02, "cm"))

p

























tree_AMF_97 +
  geom_tippoint(
    aes(fill = Taxonomy, size = sqrt(MeanAbundance)),
    shape = 21, stroke = 0.1, alpha = 0.85) +
  scale_size_continuous(range = c(0.331801/5, 22.456278/5),
                        name = "Abundance") +
  scale_fill_manual(values = palette_bestmatch, na.translate = FALSE) +
  scale_color_manual(values = palette_bestmatch, na.translate = FALSE) +
  geom_cladelab(data = genera_nodes %>% filter(clade != "Mortierella"),
                mapping = aes(node = node, label = clade, color = clade),
                parse = TRUE,
                fontsize = 2.5,
                align = TRUE,                          # all bars on one ring
                angle = "auto",
                #horizontal = FALSE, 
                offset = max_depth * 0.02,             # small offset
                offset.text = max_depth * 0.02,
                barsize = 0.4,
                show.legend = FALSE) +
  xlim(NA, max_depth * 1.2) +                          # room for labels, no more
  guides(fill = guide_legend(ncol = 1, keywidth = 0.2, keyheight = 0.2, order = 1,
                             override.aes = list(shape = 22, size = 4))) +
  geom_text2(aes(subset = !isTip & !is.na(bootstrap) & bootstrap >= 70, label = label),
             size = 2, hjust = -0.2, color = "grey30") +
  new_scale_fill() +
  geom_fruit(
    data        = dat3_97,
    geom        = geom_bar,
    mapping     =  aes(y = ID, x = Sites, fill = Sites, alpha = TotalAbundance),
    pwidth      = 0.38,
    orientation = "y",
    stat        = "identity") +
  geom_fruit(
    data    = dat2_97,
    geom    = geom_tile,
    mapping = aes(y = ID, x = Sites, alpha = Abundance, fill = Sites),
    color   = "grey90", offset = 0.04, size = 0.02
  ) +
  scale_alpha_continuous(
    range = c(0, 1),
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 4)
  )  


#theme(legend.position = "none")








tree_AMF_97 +
  geom_tippoint(
    aes(fill = Taxonomy, size = sqrt(MeanAbundance)),
    shape = 21, stroke = 0.1, alpha = 0.85) +
  scale_size_continuous(range = c(0.331801/5, 22.456278/5),
                        name = "Abundance") +
  scale_fill_manual(values = palette_bestmatch, na.translate = FALSE) +
  scale_color_manual(values = palette_bestmatch, na.translate = FALSE) +
  geom_cladelab(data = genera_nodes %>% filter(clade != "Mortierella"),
                mapping = aes(node = node, label = clade, color = clade),
                parse = TRUE, fontsize = 2.5, align = TRUE, angle = "auto",
                offset = max_depth * 0.02, offset.text = max_depth * 0.02,
                barsize = 0.4, show.legend = FALSE) +
  geom_text2(aes(subset = !isTip & !is.na(bootstrap) & bootstrap >= 70,
                 label = label),
             size = 2, hjust = -0.2, color = "grey30") +
  
  # Ring 1: site presence (categorical heatmap)
  new_scale_fill() +
  geom_fruit(
    data = dat2_97,
    geom = geom_tile,
    mapping = aes(y = ID, x = Sites, fill = Sites, alpha = Abundance),
    offset = 0.08,
    pwidth = 0.15,
    color = "white", size = 0.05
  ) +
  scale_fill_brewer(palette = "Set2", name = "Site") +
  scale_alpha_continuous(range = c(0, 1), name = "Per-site\nabundance") +
  
  # Ring 2: total abundance per site (bars or another tile layer)
  new_scale_fill() +
  geom_fruit(
    data = dat3_97,
    geom = geom_tile,
    mapping = aes(y = ID, x = Sites, fill = TotalAbundance),
    offset = 0.08,
    pwidth = 0.15,
    color = "white", size = 0.05
  ) +
  scale_fill_viridis_c(name = "Total\nabundance", trans = "sqrt") +
  
  #xlim(NA, max_depth * 1.6) +
  guides(fill = guide_legend(ncol = 1, keywidth = 0.2, keyheight = 0.2, order = 1,
                             override.aes = list(shape = 22, size = 4)))




tree_AMF_97 +
  geom_highlight(data=genera_nodes %>% filter(!clade == "Mortierella"), 
                 aes(node=node, fill=as.factor(highlight)),
                 alpha=1,
                 align="right",
                 extend=0.04,
                 show.legend=FALSE) +
  geom_cladelab(data=genera_nodes %>% filter(!clade == "Mortierella"),
                mapping=aes(node=node, label=clade),
                fontsize=3,
                align="TRUE",
                angle="auto",
                horizontal = FALSE, 
                offset=0.04,
                offset.text=0.01) +
  geom_tree(linewidth=0.3) +
  #geom_tippoint() +
  geom_tippoint(aes(size = sqrt(MeanAbundance)), 
                shape = 16, color = "darkred", stroke = 0.2, alpha = 0.85, 
                show.legend=FALSE) +
  #xlim(0, 0.35) +
  scale_fill_manual(values=c("#F5F5F5", "#ECECEC")) +
  scale_size_continuous(range = c(0.331801/5, 22.456278/5),
                        name = "sqrt(Read Abundance)") 


# order genera by their angular position around the tree
gn <- genera_nodes %>% filter(clade != "Mortierella")

# get angle for each clade's MRCA from the ggtree data

gn <- gn %>%
  left_join(tree_AMF_97$data %>% 
              dplyr::select(node, angle), by = "node") %>%
  arrange(angle) %>%
  mutate(ring = rep(c(1, 2), length.out = n()))

ring1 <- gn %>% filter(ring == 1)
ring2 <- gn %>% filter(ring == 2)


tree_AMF_97 +
  geom_tippoint(
    aes(fill = Taxonomy, size = sqrt(MeanAbundance)),
    shape = 21, stroke = 0.1, alpha = 0.85) +
  scale_size_continuous(range = c(0.331801/5, 22.456278/5),
                        name = "Abundance") +
  scale_fill_manual(values = palette_bestmatch, na.translate = FALSE) +
  scale_color_manual(values = palette_bestmatch, na.translate = FALSE) +
  geom_cladelab(data = ring2,
                mapping = aes(node = node, label = clade, color = clade),
                fontsize = 2.2, align = TRUE, angle = "auto",
                horizontal = FALSE,
                offset = max_depth * 0.05,
                offset.text = max_depth * 0.02,
                barsize = 0.5, show.legend = FALSE) +
  geom_cladelab(data = ring1,
                mapping = aes(node = node, label = clade, color = clade),
                fontsize = 2.2, align = TRUE, angle = "auto",
                horizontal = FALSE,
                offset = max_depth * 0.18,
                offset.text = max_depth * 0.02,
                barsize = 0.5, show.legend = FALSE) +
  xlim(NA, max_depth * 1.4)


















tree_AMF_97 <- 
  tree_raxml_AMF_97 %>% 
  ggtree(layout = "rectangular", branch.length = "branch.length", size = 0.1, open.angle = 5) %<+% dat1_97_mod 

scaled_tree_raxml_AMF_97 <- tree_raxml_AMF_97
scaled_tree_raxml_AMF_97$edge.length <- scaled_tree_raxml_AMF_97$edge.length + 0.05




tree_AMF_97 +
  geom_tippoint(
    aes(fill = Taxonomy, size = sqrt(MeanAbundance)),
    shape = 21, stroke = 0.1, alpha = 0.85) +
  geom_tiplab(aes(label=label), size=2) +
  scale_fill_manual(values = palette_bestmatch, na.translate = FALSE) +
  geom_text2(aes(subset = !isTip, label = node), size = 2, color = "black") +
  theme(legend.position = "none")

tree_AMF_97 +
  geom_tippoint(
    aes(fill = Taxonomy, size = sqrt(MeanAbundance)),
    shape = 21, stroke = 0.1, alpha = 0.85) +
  scale_size_continuous(range = c(0.331801/5, 22.456278/5),
                        name = "sqrt(Read Abundance)") +
  scale_fill_manual(values = palette_bestmatch, na.translate = FALSE) +
  scale_color_manual(values = palette_bestmatch, na.translate = FALSE) +
  geom_tree(linewidth=0.1) +
  geom_cladelab(data=genera_nodes %>% filter(!clade == "Mortierella"),
                mapping=aes(node=node, label=clade, color=clade),
                parse = TRUE,
                fontsize=2,
                align="FALSE",
                angle="auto",
                offset=max_depth * 0.05,
                offset.text=0.05, 
                show.legend = FALSE) +
  guides(fill = guide_legend(ncol=1, keywidth = 0.2, keyheight = 0.2, order = 1,
                             override.aes = list(shape = 22, size = 3)),
         shape = guide_legend(ncol=1, keywidth = 0.1, keyheight = 0.1, order = 2,
                              override.aes = list(shape = 21, size = 3))) 



tree_AMF_97 +
  geom_tippoint(
    aes(fill = Taxonomy, size = sqrt(MeanAbundance)),
    shape = 21, stroke = 0.1, alpha = 0.85) +
  scale_fill_manual(values = palette_bestmatch, na.translate = FALSE) +
  scale_size_continuous(range = c(0.331801/5, 22.456278/5),
                        name = "sqrt(Read Abundance)") +
  geom_highlight(data=genera_nodes %>% filter(!clade == "Mortierella"), 
                 aes(node=node, fill=clade),
                 alpha=0.2,
                 align="right",
                 extend=0.1,
                 show.legend = FALSE) + # try with TRUE as well
  geom_tree(linewidth=0.3) +
  geom_text2(aes(subset = !isTip, label = node), size = 2, color = "black") +
  geom_tiplab(aes(label=label), size=2) +
  geom_cladelab(data=genera_nodes,
                mapping=aes(node=node, label=clade, color=clade),
                parse = TRUE,
                fontsize=2,
                align="TRUE",
                angle="auto",
                offset=6,
                offset.text=0.05, 
                show.legend = FALSE) +
  guides(fill = guide_legend(ncol=1, keywidth = 0.2, keyheight = 0.2, order = 1,
                             override.aes = list(shape = 22, size = 3)),
         shape = guide_legend(ncol=1, keywidth = 0.1, keyheight = 0.1, order = 2,
                              override.aes = list(shape = 21, size = 3))) 




# Plotting all together
ggtree(tree_raxml_AMF_97, layout="circular", size = 0.1) %<+% dat1_97_mod +
  geom_highlight(data=genera_nodes %>% filter(!clade == "Mortierella"), 
                 aes(node=node, fill=as.factor(highlight)),
                 alpha=1,
                 align="right",
                 extend=0.04,
                 show.legend=FALSE) +
  geom_cladelab(data=genera_nodes %>% filter(!clade == "Mortierella"),
                mapping=aes(node=node, label=clade),
                fontsize=2,
                align="TRUE",
                angle="auto",
                offset=0.04,
                offset.text=0.01) +
  geom_tree(linewidth=0.3) +
  #geom_tippoint() +
  geom_tippoint(aes(size = sqrt(MeanAbundance)), 
                shape = 16, color = "darkred", stroke = 0.2, alpha = 0.85, 
                show.legend=FALSE) +
  #xlim(0, 0.35) +
  scale_fill_manual(values=c("#F5F5F5", "#ECECEC")) +
  scale_size_continuous(range = c(0.331801/5, 22.456278/5),
                        name = "sqrt(Read Abundance)") 


# Plotting all together
ggtree(tree_raxml_AMF_97, layout="circular", size = 0.1) %<+% dat1_97_mod +
  geom_highlight(data=genera_nodes, 
                 aes(node=node, fill=clade),
                 alpha=1,
                 align="right",
                 extend=0.04,
                 show.legend=FALSE) +
  geom_cladelab(data=genera_nodes,
                mapping=aes(node=node, label=clade),
                fontsize=2,
                align="TRUE",
                angle="auto",
                offset=0.04,
                offset.text=0.01) +
  geom_tree(linewidth=0.3) +
  geom_tiplab(aes(fill=Taxonomy), size=2, show.legend = FALSE) +
  scale_fill_manual(values = palette_bestmatch, na.translate = FALSE)
#geom_tippoint() +




# Need to modufy this if want to add the title!
ggsave(
  file.path(data_path, "results/Fig_2_betadiv_wide.pdf"),
  plot = grid.arrange(
    Figure_2_beta,
    top = text_grob("BETA DIVERSITY", size = 12, face = "bold")
  ),
  device = "pdf"
)




tree_AMF_97 <- 
  tree_raxml_AMF_97 %>% 
  ggtree(layout = "circular", branch.length = "none", size = 0.1, open.angle = 5) 

tree_AMF_97

tree_AMF_97 +
  geom_text2(aes(subset = !isTip, label = node), size = 2, color = "black")

tree_AMF_97 <- 
  tree_AMF_97 %<+% dat1_97 +
  geom_tippoint(
    aes(fill = BestMatch, size = sqrt(MeanAbundance)),
    shape = 21, stroke = 0.1, alpha = 0.85) +
  scale_fill_manual(
    values = bestmatch_colors_97,
    guide  = guide_legend(
      keywidth = 0.5, keyheight = 0.5, order = 1,
      override.aes = list(shape = 21, size = 3)
    ),
    na.translate = FALSE
  ) +
  scale_size_continuous(
    range = c(0.1, 4),
    name  = "sqrt(Read Abundance)",
    guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 2,
                         override.aes = list(shape = 21))
  )

tree_AMF_97 <- 
  tree_AMF_97 +
  geom_tiplab2(
    aes(label = BestMatch),
    size = 2,
    align = TRUE,
    offset = 0.02
  )

tree_AMF_97 + theme(legend.position="none")


# Test tree 1 (warnings aren't harmful)














# Improvements -----------------------------------------------------------------

tree_AMF_97 <- 
  ggtree(tree_raxml_97_rooted, linewidth=0.5) %<+% dat1_97 +
  geom_tiplab(aes(label=label), size=2)

tree_AMF_97

#Make dataframe for clade nodes
genera_nodes <- data.frame(
  clade=unique(dat1_97$BestMatch),
  node=NA
)

genera_nodes

#Find the most recent common ancestor for each clade
for (i in 1:length(genera_nodes$clade)) {
  
  genera_nodes$node[i] <- MRCA(
    tree_raxml_97_rooted,
    dat1_97$label[dat1_97$BestMatch == genera_nodes$clade[i]]
  )
  
}

tree_raxml_97_rooted

tree_AMF_97 <- 
  ggtree(tree_raxml_97_rooted, linetype=NA) %<+% dat1_97 +
  geom_highlight(data=genera_nodes, 
                 aes(node=node, fill=clade),
                 alpha=1,
                 align="right",
                 extend=0.1,
                 show.legend=FALSE) + # try with TRUE as well
  geom_tree(linewidth=0.5) +
  geom_tiplab(aes(label=BestMatch), size=2)

tree_AMF_97 + xlim(0, 3) 


tree_AMF_97 + layout_circular() + 
  theme(legend.position=c(0.05, 0.7)) 

# Alternating highlight colors

genera_nodes <- genera_nodes[match(tree_AMF_97$data %>%
                                     filter(isTip == "TRUE") %>%
                                     arrange(y) %>%
                                     pull(BestMatch) %>%
                                     unique(),
                                   genera_nodes$clade),]
#Add column with alternating binary value
genera_nodes$highlight <- rep(c(0,1),
                              length.out=length(genera_nodes$clade))
genera_nodes


#Add highlights
tree_AMF_97 <- 
  ggtree(tree_raxml_97_rooted, linetype=NA) %<+% dat1_97 +
  geom_highlight(data=genera_nodes, 
                 aes(node=node, fill=as.factor(highlight)),
                 alpha=1,
                 align="right",
                 extend=0.1,
                 show.legend=FALSE) +
  geom_tree(linewidth=0.5) +
  geom_tiplab(aes(label=BestMatch), size=2) +
  scale_fill_manual(values=c("#F5F5F5", "#ECECEC"))

tree_AMF_97

#Add clade labels
tree_AMF_97 +
  geom_cladelab(data=genera_nodes,
                mapping=aes(node=node, label=clade),
                fontsize=2,
                align=TRUE,
                offset=0.1,
                offset.text=0.02)


tree_raxml_97_rooted %>% 
  ggtree(layout = "circular", branch.length = "none", size = 0.1, open.angle = 5)































tree_AMF_merged <- tree_AMF_merged + 
  geom_tiplab(
    aes(label = BestMatch),       # We don't want text, just the line
    size = 2,
    linetype = "dotted",   # "33" or "dotted"
    linewidth = 0.2,       # Keep it subtle
    color = "grey50", 
    align = TRUE,          # This draws the line from tip to the alignment margin
    offset = -1             # Adjust if there is a gap
  )



# Adding geneus labels 
genus_nodes <- dat1 %>%
  filter(!is.na(Taxa)) %>%                    # genus-level only
  group_by(Taxa) %>%
  summarise(tips = list(label), .groups = "drop") %>%
  mutate(
    ntips = lengths(tips),
    node  = map_int(tips, function(t) {
      t <- unlist(t)
      if (length(t) > 1) ape::getMRCA(tree_raxml_AMF, t)
      else which(tree_raxml_AMF$tip.label == t)
    })
  )

# diagnostic — check for genera with MRCA near the root
root_node <- length(tree_raxml_AMF$tip.label) + 1
genus_nodes %>% dplyr::select(Taxa, ntips, node) %>% 
  arrange(node) %>% 
  print(n = Inf)

# optional: drop genera whose MRCA is suspiciously close to root
# adjust threshold after inspecting the diagnostic above
genus_nodes_plot <- genus_nodes %>% filter(node > root_node + 10)

geom_cladelab(
  data       = genus_nodes_plot,
  mapping    = aes(node = node, label = Taxa),
  hjust      = 0.5,   angle      = "auto",
  barsize    = NA,    horizontal = FALSE,
  fontsize   = 1.4,   fontface   = "italic"
) +
  
  
  
  # Get the MRCA node for each genus group
  genus_nodes <- dat1 %>%
  group_by(Genus) %>%
  summarise(tips = list(label), .groups = "drop") %>%
  filter(Genus != "Unclassified") %>%
  rowwise() %>%
  mutate(
    node = ifelse(
      length(tips) > 1,
      ape::getMRCA(tree_raxml_AMF, unlist(tips)),  # use ape:: explicitly + unlist()
      which(tree_raxml_AMF$tip.label == unlist(tips)[[1]])
    )
  ) %>%
  ungroup()

genus_nodes

nodedf <- data.frame(node = genus_nodes$node)
labdf  <- data.frame(
  node  = genus_nodes$node,
  label = genus_nodes$Genus
)


physeq_AMF_clean
physeq_AMF_rare
physeq_AMF_rare_Gen

tree_raxml_AMF <- phy_tree(physeq_AMF_rare)
str(tree_raxml_AMF)

# Generating metadata for tree 
# Generating abundance, taxonomy, metadata table for the tree plot 

melted_AMF <- 
  physeq_AMF_rare %>%
  psmelt() %>%
  arrange(Genus)

head(melted_AMF)
dim(melted_AMF)

any(physeq_AMF_rare@tax_table == "Mortierellaceae")
is.na(physeq_AMF_rare@tax_table) 

# Create the base taxonomy mapping (1 row per OTU)
taxonomy_4_tree <- 
  melted_AMF %>%
  dplyr::select(OTU, Family, Genus, BestMatch) %>%
  distinct()

unique(taxonomy_4_tree$BestMatch)

# Calculate Total Abundance and Frequencies
tree_metadata <- 
  melted_AMF %>%
  mutate(present = ifelse(Abundance > 0, 1, 0)) %>%
  group_by(OTU) %>%
  dplyr::summarise(
    Abundance = sum(Abundance),
    Fertilized = sum(present[fert_status == "Fertilized"]),
    Control    = sum(present[fert_status == "Control"]),
    `Lux Arbor` = sum(present[site == "Lux Arbor"]),
    `Lake City`   = sum(present[site == "Lake City"]),
    Escanaba   = sum(present[site == "Escanaba"]),
    Rhinelander = sum(present[site == "Rhinelander"]),
    Hancock    = sum(present[site == "Hancock"])
  ) %>%
  left_join(taxonomy_4_tree, by = "OTU") %>%
  dplyr::select(OTU, Family, Genus, BestMatch, Abundance, everything())

tree_metadata

unique(tree_metadata$Genus)
unique(tree_metadata$BestMatch)

palette_bestmatch <-
  c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77",
    "#283dff","#636bb7","#bfc5ff","#014443","#195637",
    "#117744","#60ffaf","#b7ffdb","#825121","#ea7f17",
    "#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d",
    "#5b5b19","#fcfc00","#ffff9e","#ffb7ef","#fa7efc",
    "#ae09ea","#521899","#1e0047")


tree_AMF <- 
  ggtree(tree_raxml_AMF, layout="rectangular", branch.length="branch.length") %<+% tree_metadata +
  geom_tippoint(aes(color = BestMatch, size = Abundance), alpha = 0.8)  +
  scale_color_manual(values = palette_bestmatch) +
  guides(color = guide_legend(ncol = 1, override.aes = list(shape = 15, size=5))) 


tree_AMF + 
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = OTU, fill = Fertilized), 
    pwidth = 0.1) +
  scale_fill_gradient2() + 
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = OTU, fill = Control), 
    pwidth = 0.1) +
  scale_fill_gradient2()


abund_4_tree <- 
  melted_AMF %>%
  group_by(OTU, site, Family, Genus, BestMatch) %>% # no fert_status
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  rename(label = OTU, abundance = Abundance) %>%
  mutate(Family = as.character(Family))

head(abund_4_tree)
dim(abund_4_tree)
any(is.na(abund_4_tree))
levels(as.factor(abund_4_tree$Family))

# Tree to tibble 
tree_data <- as_tibble(tree_raxml_AMF)
tree_data

tree_metadata <- 
  abund_4_tree %>%
  tidyr::pivot_wider(
    id_cols = c(label, Family, Genus, BestMatch),
    names_from = "site",
    values_from = abundance) %>% 
  mutate(Genus = ifelse(is.na(Genus), "Unclassified", Genus))

dim(tree_metadata)
tree_metadata

# create list of OTUs per Genus
genus_list <- 
  tree_metadata %>%
  mutate(Genus = ifelse(is.na(Genus), "Unclassified", Genus)) %>%
  group_by(Genus) %>%
  summarise(tips = list(label), .groups = "drop") %>%
  { setNames(.$tips, .$Genus) }

tree_grouped <- groupOTU(tree_raxml_AMF, genus_list)

genus_list <- setNames(genus_list$tips, genus_list$Genus)

# group tree
tree_grouped <- groupOTU(tree_raxml_AMF, genus_list)
tree_grouped

tree_grouped@data$group[tree_grouped@data$group == "0"] <- NA

# 2) Build plot version 2 -
ggtree(tree_raxml_AMF) %<+% tree_metadata +
  geom_tippoint(aes(color = Genus, size = abundance), alpha = 0.8)  +
  #geom_tiplab(aes(color = Genus), size = 2) +
  theme(legend.position="right") +
  scale_color_manual(values = palette_taxa,na.translate = FALSE) +
  guides(color = guide_legend(ncol = 1, override.aes = list(shape = 15, size=3))) 


tree_AMF <- 
  ggtree(tree_grouped, aes(color = group), size = 0.3) +
  scale_color_discrete(na.translate = FALSE)

tree_AMF

tree_AMF <- ggtree(as.treedata(tree_tbl2), layout="rectangular", branch.length="none")

tree_AMF <- tree_AMF + geom_tree(aes(color = Genus)) +
  theme(legend.position="right") +
  scale_color_manual(values = palette_taxa,na.translate = FALSE) +
  guides(color = guide_legend(ncol = 1, override.aes = list(shape = 15, size=3))) 

tree_AMF



genus_list <- 
  tree_metadata %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus) %>%
  summarise(tips = list(label), .groups = "drop")


tree_tbl2 <- tree_data %>%
  left_join(tree_metadata %>% distinct(label, .keep_all = TRUE),
            by = "label")

as.treedata(tree_tbl2)

tree_grouped@data$group[tree_grouped@data$group == "0"] <- NA

ggtree(tree_grouped, aes(color = group), size = 0.3, layout="rectangular") +
  scale_color_discrete(na.translate = FALSE)

setdiff(tree_grouped$tip.label, tree_metadata$label)
setdiff(tree_metadata$label, tree_grouped$tip.label)


ggtree(tree_grouped, aes(color = group), size = 0.3) +
  scale_color_discrete(breaks = function(x) x[x != "0"])

tree_AMF <- 
  ggtree(tree_raxml_AMF, layout = "rectangular", branch.length = "none") %<+% abund_4_tree +
  geom_tippoint(aes(color = Genus), size = 1) +
  theme(legend.position = "right")

tree_AMF <- ggtree(tree_raxml_AMF, layout = "rectangular", branch.length = "none") %<+% abund_4_tree +
  geom_tiplab(aes(color = Family), size = 2) +
  theme(legend.position = "right")

tree_AMF



# Add the Prevalence Heatmap Ring
tree_AMF <- tree_AMF + geom_fruit(
  geom = geom_tile,
  mapping = aes(y=label, x=site, fill=abundance),
  pwidth = 0.1,
  axis.params = list(axis="x", text.size=2),
  grid.params = list()
) +
  scale_fill_viridis_c()

tree_AMF
taxonomy_4_tree
# Add the Abundance Bar Ring
tree_AMF <- tree_AMF + geom_fruit(
  geom = geom_bar,
  mapping = aes(y=label, x=abundance, fill=site),
  stat="identity",
  orientation="y",
  offset=0.15,
  pwidth = 0.1
)

tree_AMF

tree_AMF <- tree_AMF + 
  geom_star(
    aes(fill=Family, starshape=fert_status, size=Abundance), 
    starstroke=0.3,
    na.rm = TRUE
  ) + 
  scale_starshape_manual(values=c(15, 1)) # 15 is a square, 1 is a circle

tree_AMF


taxonomy_4_tree <- 
  as.data.frame(tax_table(physeq_AMF_rare)) %>% 
  dplyr::select(-Query, -Kingdom, -S_score) %>% 
  rownames_to_column("tip.label")

head(taxonomy_4_tree)
dim(taxonomy_4_tree)

all(tree_raxml_AMF$tip.label %in% taxonomy_4_tree$tip.label)
all(taxonomy_4_tree$tip.label %in% tree_raxml_AMF$tip.label)

setdiff(tree_raxml_AMF$tip.label, taxonomy_4_tree$tip.label)
setdiff(taxonomy_4_tree$tip.label, tree_raxml_AMF$tip.label)



metadata_4_tree <- 
  metadata_99 %>% rownames_to_column("label")

str(metadata_4_tree)
head(metadata_4_tree)

taxonomy_4_tree <- as.data.frame(tax_table(physeq_AMF_rare)) %>% 
  dplyr::select(-Query, -Kingdom, -S_score) %>% 
  rownames_to_column("label")

head(taxonomy_4_tree)
str(taxonomy_4_tree)



# Start over 
# https://yulab-smu.top/treedata-book/chapter10.html

library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(ggnewscale)
library(dplyr)
library(phyloseq)

tree_raxml_AMF <- phy_tree(physeq_AMF_rare)

dat1 <- melted_AMF %>%
  group_by(OTU, Genus, Phylum, Family) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus)) %>% 
  rename(label = OTU)

dat2 <- melted_AMF %>%
  group_by(OTU, site) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  rename(ID = OTU, Sites = site)


dat3 <- melted_AMF %>%
  group_by(OTU, site) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  rename(ID = OTU, Sites = site)

# ── 5. Color palette for genera ───────────────────────────────────────────────
genus_colors <- c(
  "Acaulospora"    = "#E41A1C", "Ambispora"      = "#4575B4",
  "Archaeospora"   = "#762A83", "Cetraspora"     = "#F4A582",
  "Complexispora"  = "#1B7837", "Dentiscutata"   = "#FFFFBF",
  "Diversispora"   = "#2166AC", "Dominikia"      = "#000080",
  "Entrophospora"  = "#808080", "Funneliformis"  = "#00441B",
  "Gigaspora"      = "#FDB863", "Glomus"         = "#9ECAE1",
  "Microdominikia" = "#1A1A2E", "Microkamienskia"= "#BDBDBD",
  "Oehlia"         = "#66C2A5", "Paraglomus"     = "#F46D43",
  "Rhizophagus"    = "#FF69B4", "Scutellospora"  = "#4D0000",
  "Septoglomus"    = "#74C476", "Silvaspora"     = "#8B8B00",
  "Unclassified"   = "#000000"
)

site_colors <- c(
  "Lux Arbor"  = "#0000FF",
  "Escanaba"   = "#FFA500",
  "Maple River"= "#FF0000"   # add/rename to match your actual site names
)

# Build the plot
tree_AMF <- ggtree(tree_raxml_AMF, layout = "fan", size = 0.15, open.angle = 5)
tree_AMF

# Layer 1: tip points — color = Genus, size = mean abundance
tree_AMF <- tree_AMF %<+% dat1 +
  geom_tippoint(
    aes(fill = Genus, size = log1p(MeanAbundance)),
    shape = 21, stroke = 0.1, alpha = 0.85) +
  scale_fill_manual(
    values = palette_bestmatch,
    guide  = guide_legend(
      keywidth = 0.5, keyheight = 0.5, order = 1,
      override.aes = list(shape = 21, size = 3)
    ),
    na.translate = FALSE
  ) +
  scale_size_continuous(
    range = c(0.1, 4),
    name  = "log(Abundance + 1)",
    guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 2,
                         override.aes = list(shape = 21))
  )

# Layer 2: site abundance heatmap ring
tree_AMF <- 
  tree_AMF +
  new_scale_fill() +
  geom_fruit(
    data    = dat2,
    geom    = geom_tile,
    mapping = aes(y = ID, x = Sites, alpha = Abundance, fill = Sites),
    color   = "grey90", offset = 0.04, size = 0.02
  ) +
  scale_alpha_continuous(
    range = c(0, 1),
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 4)
  ) +
  scale_fill_manual(
    values = palette_site,
    guide  = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 3)
  )

# Layer 3: total abundance bar ring
tree_AMF <- 
  tree_AMF +
  new_scale_fill() +
  geom_fruit(
    data        = dat3,
    geom        = geom_bar,
    mapping     = aes(y = ID, x = TotalAbundance, fill = Sites),
    pwidth      = 0.38,
    orientation = "y",
    stat        = "identity"
  ) +
  scale_fill_manual(
    values = palette_site,
    guide  = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 3)
  )

# ── 7. Theme ──────────────────────────────────────────────────────────────────
tree_AMF <- 
  tree_AMF +
  geom_treescale(fontsize = 2, linesize = 0.3) +
  theme(
    legend.position   = c(0.93, 0.5),
    legend.background = element_rect(fill = NA),
    legend.title      = element_text(size = 6.5),
    legend.text       = element_text(size = 4.5),
    legend.spacing.y  = unit(0.02, "cm")
  )

tree_AMF

tree_AMF + layout_rectangular() + 
  theme(legend.position=c(.05, .7))





tree <- tree_hmptree
dat1 <- df_tippoint
dat2 <- df_ring_heatmap
dat3 <- df_barplot_attr

# adjust the order
dat2$Sites <- factor(dat2$Sites, 
                     levels=c("Stool (prevalence)", "Cheek (prevalence)",
                              "Plaque (prevalence)","Tongue (prevalence)",
                              "Nose (prevalence)", "Vagina (prevalence)",
                              "Skin (prevalence)"))
dat3$Sites <- factor(dat3$Sites, 
                     levels=c("Stool (prevalence)", "Cheek (prevalence)",
                              "Plaque (prevalence)", "Tongue (prevalence)",
                              "Nose (prevalence)", "Vagina (prevalence)",
                              "Skin (prevalence)"))
# extract the clade label information. Because some nodes of tree are
# annotated to genera, which can be displayed with high light using ggtree.
nodeids <- nodeid(tree, tree$node.label[nchar(tree$node.label)>4])
nodedf <- data.frame(node=nodeids)
nodelab <- gsub("[\\.0-9]", "", tree$node.label[nchar(tree$node.label)>4])
# The layers of clade and hightlight
poslist <- c(1.6, 1.4, 1.6, 0.8, 0.1, 0.25, 1.6, 1.6, 1.2, 0.4,
             1.2, 1.8, 0.3, 0.8, 0.4, 0.3, 0.4, 0.4, 0.4, 0.6,
             0.3, 0.4, 0.3)
labdf <- data.frame(node=nodeids, label=nodelab, pos=poslist)

# The circular layout tree.
p <- ggtree(tree, layout="fan", size=0.15, open.angle=5) +
  geom_hilight(data=nodedf, mapping=aes(node=node),
               extendto=6.8, alpha=0.3, fill="grey", color="grey50",
               size=0.05) +
  geom_cladelab(data=labdf, 
                mapping=aes(node=node, 
                            label=label,
                            offset.text=pos),
                hjust=0.5,
                angle="auto",
                barsize=NA,
                horizontal=FALSE, 
                fontsize=1.4,
                fontface="italic"
  )

p <- p %<+% dat1 + geom_star(
  mapping=aes(fill=Phylum, starshape=Type, size=Size),
  position="identity",starstroke=0.1) +
  scale_fill_manual(values=c("#FFC125","#87CEFA","#7B68EE","#808080",
                             "#800080", "#9ACD32","#D15FEE","#FFC0CB",
                             "#EE6A50","#8DEEEE", "#006400","#800000",
                             "#B0171F","#191970"),
                    guide=guide_legend(keywidth = 0.5, 
                                       keyheight = 0.5, order=1,
                                       override.aes=list(starshape=15)),
                    na.translate=FALSE)+
  scale_starshape_manual(values=c(15, 1),
                         guide=guide_legend(keywidth = 0.5, 
                                            keyheight = 0.5, order=2),
                         na.translate=FALSE)+
  scale_size_continuous(range = c(1, 2.5),
                        guide = guide_legend(keywidth = 0.5, 
                                             keyheight = 0.5, order=3,
                                             override.aes=list(starshape=15)))

p <- p + new_scale_fill() +
  geom_fruit(data=dat2, geom=geom_tile,
             mapping=aes(y=ID, x=Sites, alpha=Abundance, fill=Sites),
             color = "grey50", offset = 0.04,size = 0.02)+
  scale_alpha_continuous(range=c(0, 1),
                         guide=guide_legend(keywidth = 0.3, 
                                            keyheight = 0.3, order=5)) +
  geom_fruit(data=dat3, geom=geom_bar,
             mapping=aes(y=ID, x=HigherAbundance, fill=Sites),
             pwidth=0.38, 
             orientation="y", 
             stat="identity",
  ) +
  scale_fill_manual(values=c("#0000FF","#FFA500","#FF0000",
                             "#800000", "#006400","#800080","#696969"),
                    guide=guide_legend(keywidth = 0.3, 
                                       keyheight = 0.3, order=4))+
  geom_treescale(fontsize=2, linesize=0.3, x=4.9, y=0.1) +
  theme(legend.position=c(0.93, 0.5),
        legend.background=element_rect(fill=NA),
        legend.title=element_text(size=6.5),
        legend.text=element_text(size=4.5),
        legend.spacing.y = unit(0.02, "cm"),
  )
p

