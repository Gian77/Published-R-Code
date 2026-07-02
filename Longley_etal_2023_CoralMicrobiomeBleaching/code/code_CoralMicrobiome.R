# ************ DATA ANALYSIS ***************************************--------------
# Project name: 
# Manuscript:   Coral Microbiome - bleached vs unbleached 
# Authors:      Longley et al. 
# Affiliation:  Michigan State University - GLBRC
# Journal:      Plos Biology
# Date:         Jan 10, 2023
# *********************************************************************------------

options(scipen = 9999)
options(max.print=100000000)

# *********************************************************************------------
# LOADING PACKAGES ----------------------------------------------------------------
library(phyloseq); packageVersion("phyloseq")
library(Biostrings)
library(ggplot2)
library(ggpubr)
library(ape)
library(dplyr)
library(randomForest)
library(rfUtilities) 
library(compositions)
library(mlbench)
library(caret)
library(labgeo)
library(compositions)
library(evoPalette)
library(gridExtra)

# ********************-------------------------------------------------------------
# **** CORAL COMMUNITY COMPOSITION **** --------------------------------------------

# Coral community composition ------------------------------------------------------
coral_com <- read.csv("Victor/community.csv", header = TRUE)
coral_com$Genus <- 
  factor(coral_com$Genus, levels = c("Acropora", "Montipora", "Pocillopora", "Porites", "Others"))
coral_com


plot_comm <-
  ggplot(coral_com, aes(x=Genus, y=abund, fill=Genus)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=abund-stderr, ymax=abund+stderr), width=0.2,position=position_dodge(0.9)) +
  facet_wrap(~Site) +
  theme_classic() +
  scale_fill_manual(values = c(paletteCB4, "grey70")) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
        legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8, hjust = 0, face = "italic"), 
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, face = "italic", size = 8,hjust = 1.1, vjust = 1.1),
        axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5),
        axis.title = element_text(angle = 0, size = 8, face = "bold"),
        legend.position = "none") +
  labs(title = "Coral abundance", y= "Abundance (%)", x=NULL)

plot_comm


coral_res <- read.csv("Victor/extra_figure.csv", header = TRUE)

coral_red_df <-
  coral_res %>% 
  rename(Genus=X) %>% 
  mutate(Group=as.factor(Group)) %>% 
  as_tibble() %>% 
  filter(.data = ., Genus != "Others") 


coral_red_df$Group <-
  factor(coral_red_df$Group, levels = c("Bleached","Unbleached","Dead")) 

plot_res <-
  coral_red_df %>% 
  ggplot(data = ., aes(x=Genus, y=Value, fill=Genus)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Value-SD, ymax=Value+SD), width=0.2,position=position_dodge(0.9)) +
  facet_wrap(~Group, ncol = 1) +
  theme_classic() +
  scale_fill_manual(values = c(paletteCB4, "grey70")) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
        legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8, hjust = 0, face = "italic"), 
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, face = "italic", size = 8,hjust = 1.1, vjust = 1.1),
        axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5),
        axis.title = element_text(angle = 0, size = 8, face = "bold"),
        legend.position = "none") +
  labs(title = "Overall coral Status", y= "Abundance (%)", x=NULL)

plot_res

plot_res <-
  plot_res + stat_summary(
    geom = 'text',
    label = coral_red_df$Letter,
    fun = max,
    aes(y = 98),
    size = 3)



coral_bent <- read.csv("Victor/benthic.csv", header = TRUE)
coral_bent

plot_bent <-
  ggplot(coral_bent, aes(x=Category, y=cover)) + 
  geom_bar(stat="identity", fill= "grey70", position=position_dodge()) +
  geom_errorbar(aes(ymin=cover-stderr, ymax=cover+stderr), width=0.2,position=position_dodge(0.9)) +
  facet_wrap(~Site) +
  theme_classic() +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
        legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8, hjust = 0), 
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, size = 8,hjust = 1.1, vjust = 1.1),
        axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5),
        axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  labs(title = "Benthic cover", y= "Cover (%)", x=NULL)

plot_bent



coral_status <- read.csv("Victor/status.csv", header = TRUE)
coral_status$Genus <- 
  factor(coral_status$Genus, levels = c("All","Acropora", "Montipora", "Pocillopora", "Porites", "Others"))
coral_status

coral_status <-
  left_join(
    tibble::rownames_to_column(
      coral_status %>%
        select(Site, Genus, Healthy, Bleached, Dead) %>%
        reshape2::melt()),
    tibble::rownames_to_column(
      coral_status %>%
        select(Site, Genus, H_stderr, B_stdrerr, D_stderr) %>%
        reshape2::melt()), by="rowname") %>%
  select(Site.x,Genus.x,variable.x,value.x, value.y)

colnames(coral_status) <-
  c("Site", "Genus", "Status", "abund", "stderr")
coral_status$Status <-
  factor(coral_status$Status, levels = c("Healthy", "Bleached", "Dead"))


plot_status <-
  ggplot(coral_status, aes(x=Status, y=abund, fill= Genus)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=abund-stderr, ymax=abund+stderr), width=0.2,position=position_dodge(0.9)) +
  facet_grid(~ Site ~ Genus, labeller = ) +
  theme_classic() +
  scale_fill_manual(values = c("black", paletteCB4, "grey70")) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
        legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8, hjust = 0), 
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, size = 8,hjust = 1.1, vjust = 1.1),
        axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5),
        axis.title = element_text(angle = 0, size = 8, face = "bold"),
        legend.position = "bottom") +
  guides(fill=guide_legend(nrow=1)) +
  labs(title = "Coral status at reef study locations", y= "Abundance (%)", x=NULL)

plot_status


label_letters=c("b","a","b","","b","a","","","","","","","b","a","b","a","a","b",
                "b","a","a","","a","a","a","b","b","","","","","a","b","","","")

plot_status <-
  plot_status + stat_summary(
    geom = 'text',
    label = label_letters,
    fun = max,
    aes(y = 98),
    size = 3)


Fig_1_Coral <-
  ggarrange(
    ggarrange(
      plot_comm, 
      plot_bent, 
      widths = c(1, 1.4),
      align = "hv", 
      labels = c("E", "F")),
    plot_status,
    ncol = 1,
    nrow = 2,
    labels = c("","G"),
    heights = c(1, 1.5),
    align = "hv",
    legend = "bottom")

# **** FIGURE 1 - community composition --------------------------------------------

Fig_1_Coral_final <-
  ggarrange(Fig_1_Coral, plot_res,
            ncol = 2,
            nrow = 1,
            align = "hv",
            widths = c(0.8, 0.17), 
            labels = c("", "H"))

Fig_1_Coral_final

# Temperature graphs ---------------------------------------------------------------

coral_temp <- read.csv("Victor/temperature.csv", header = TRUE)
coral_temp$rows <- as.numeric(rownames(coral_temp))
coral_temp$Date <- 
  as.factor(paste(coral_temp$Year, coral_temp$Month2, coral_temp$Day, sep = "-"))
coral_temp$Date2 <- 
  as.Date(coral_temp$Date, "%Y-%b-%d")
head(coral_temp)
tail(coral_temp)

plot_temp <-
  ggplot(coral_temp, aes(x=Date2, y=SST)) + 
  geom_line() +
  theme_classic() +
  scale_x_date(date_breaks = "2 years" , date_labels = "%Y") +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
        legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8, hjust = 0), 
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5),
        axis.title = element_text(angle = 0, size = 8, face = "bold"),
        legend.position = "bottom") +
  labs(title = NULL, y= "SST (°C)", x=NULL)


plot_temp


plot_dhw <-
  ggplot(coral_temp, aes(x=Date2, y=DHW)) + 
  geom_line() +
  theme_classic() +
  scale_x_date(date_breaks = "2 years" , date_labels = "%Y") +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
        legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8, hjust = 0, face = "italic"), 
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, face = "italic", size = 8,hjust = 1.1, vjust = 1.1),
        axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5),
        axis.title = element_text(angle = 0, size = 8, face = "bold"),
        legend.position = "bottom") +
  labs(title = NULL, y= "DHW", x=NULL)

plot_dhw +
  annotate("rect", xmin = as.Date("1999-12-01"), xmax = as.Date("2000-05-31"),
           alpha = 0.3, fill = "red", ymin = -Inf, ymax = Inf) +
  annotate("rect", xmin = as.Date("2000-12-01"), xmax = as.Date("2001-05-31"),
           alpha = 0.3, fill = "yellow", ymin = -Inf, ymax = Inf) +
  annotate("rect", xmin = as.Date("2001-12-01"), xmax = as.Date("2002-05-31"),
           alpha = 0.3, fill = "red", ymin = -Inf, ymax = Inf) +
  annotate("rect", xmin = as.Date("2005-12-01"), xmax = as.Date("2006-05-31"),
           alpha = 0.3, fill = "yellow", ymin = -Inf, ymax = Inf) +
  annotate("rect", xmin = as.Date("2008-12-01"), xmax = as.Date("2009-05-31"),
           alpha = 0.3, fill = "yellow", ymin = -Inf, ymax = Inf) +
  annotate("rect", xmin = as.Date("2013-12-01"), xmax = as.Date("2014-05-31"),
           alpha = 0.3, fill = "yellow", ymin = -Inf, ymax = Inf) +
  annotate("rect", xmin = as.Date("2014-12-01"), xmax = as.Date("2015-05-31"),
           alpha = 0.3, fill = "yellow", ymin = -Inf, ymax = Inf) +
  annotate("rect", xmin = as.Date("2015-12-01"), xmax = as.Date("2016-05-31"),
           alpha = 0.3, fill = "red", ymin = -Inf, ymax = Inf) -> plot_dhw

plot_dhw


ggarrange(
  plot_temp,
  plot_dhw,
  ncol = 1,
  nrow = 2, 
  labels = c("A", "B")
)

# mean daily temperature
coral_mean_t <- read.csv("Victor/mean_daily_temp.csv", header = TRUE)
coral_mean_t$Date2 <- 
  as.Date(coral_mean_t$Date, "%d-%m-%Y")
head(coral_mean_t)
tail(coral_mean_t)


plot_coral_mean <-
  data.frame(
    coral_mean_t %>%
      select(Vatuoalali, Votua, Satellite) %>%
      reshape2::melt(),
    Date = rep(coral_mean_t$Date2, 3)) %>%
  ggplot(aes(x=Date, y=value, color=variable)) + 
  geom_line() +
  theme_classic() +
  scale_x_date(date_breaks = "weeks" , date_labels = "%Y-%m-%d",
               limits = c(as.Date("2015-12-01"), as.Date("2016-05-31"))) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
        legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8, hjust = 0), 
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5),
        axis.title = element_text(angle = 0, size = 8, face = "bold"),
        legend.position = "bottom") +
  labs(title = NULL, y= "Mean daily T (°C)", x=NULL)

plot_coral_mean +
  annotate("rect", xmin = as.Date("2016-02-18"), xmax = as.Date("2016-03-07"),
           alpha = 0.3, fill = "red", ymin = -Inf, ymax = Inf) +
  scale_colour_manual(values = c("red", "green", "black"), 
                      labels = c("Vatuolalai","Votua", "Satellite")) -> plot_coral_mean

plot_coral_mean

# project period DHW ---------------------------------------------------------------
coral_hot_s <- read.csv("Victor/hot_season_project_period.csv", header = TRUE)
coral_hot_s$Date2 <- as.Date(coral_hot_s$Date, "%d-%m-%Y")
head(coral_hot_s)
tail(coral_hot_s)

plot_coral_hot <-
  data.frame(
    coral_hot_s %>%
      select(Vatuoalali, Votua, Satellite) %>%
      reshape2::melt(),
    Date = rep(coral_hot_s$Date2, 3)) %>%
  ggplot(aes(x=Date, y=value, color=variable)) + 
  geom_line() +
  theme_classic() +
  scale_x_date(date_breaks = "weeks" , date_labels = "%Y-%m-%d",
               limits = c(as.Date("2015-12-01"), as.Date("2016-05-31"))) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
        legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8, hjust = 0), 
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, face = "italic", size = 8,hjust = 1.1, vjust = 1.1),
        axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5),
        axis.title = element_text(angle = 0, size = 8, face = "bold"),
        legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
  labs(title = NULL, y= "DHW", x=NULL)

plot_coral_hot


plot_coral_hot +
  annotate("rect", xmin = as.Date("2016-02-18"), xmax = as.Date("2016-03-07"),
           alpha = 0.3, fill = "red", ymin = -Inf, ymax = Inf) +
  scale_colour_manual(values = c("red", "green", "black"),
                      labels = c("Vatuolalai","Votua", "Satellite")) -> plot_coral_hot

plot_coral_hot

# Max daily temp -------------------------------------------------------------------
coral_max <- read.csv("Victor/max_daily_temp.csv", header = TRUE)
coral_max$Date2 <- as.Date(coral_max$Date, "%d-%m-%Y")
head(coral_max)
tail(coral_max)


plot_coral_max <-
  data.frame(
    coral_max %>%
      select(Vatuoalali, Votua) %>%
      reshape2::melt(),
    Date = rep(coral_max$Date2, 2)) %>%
  ggplot(aes(x=Date, y=value, color=variable)) + 
  geom_line() +
  theme_classic() +
  scale_x_date(date_breaks = "weeks" , date_labels = "%Y-%m-%d",
               limits = c(as.Date("2015-12-01"), as.Date("2016-05-31"))) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm"), 
        legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8, hjust = 0), 
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5),
        axis.title = element_text(angle = 0, size = 8, face = "bold"),
        legend.position = "bottom") +
  labs(title = NULL, y= "Max daily T (°C)", x=NULL)


plot_coral_max +
  annotate("rect", xmin = as.Date("2016-02-18"), xmax = as.Date("2016-03-07"),
           alpha = 0.3, fill = "red", ymin = -Inf, ymax = Inf) +
  scale_colour_manual(values = c("red", "green"),
                      labels = c("Vatuolalai","Votua")) -> plot_coral_max
plot_coral_max


ggarrange(
  plot_coral_hot,
  plot_coral_max,
  plot_coral_mean,
  ncol = 1,
  nrow = 3, 
  labels = c("C", "D", "E"), 
  common.legend = TRUE, 
  legend = "bottom",
  align = "hv"
)

# FIGURE 2 - weather ---------------------------------------------------------------

ggarrange(
  ggarrange(
    plot_temp,
    plot_dhw,
    ncol = 1,
    nrow = 2, 
    labels = c("A", "B")
  ),
  ggarrange(
    plot_coral_mean,
    plot_coral_max,
    plot_coral_hot,
    ncol = 1,
    nrow = 3, 
    labels = c("C", "D", "E"), 
    common.legend = TRUE, 
    legend = "bottom"), 
  ncol = 1,
  nrow = 2,
  align = "hv",
  heights = c(1, 1.5))


ggarrange(
  plot_temp,
  plot_dhw,
  plot_coral_mean,
  plot_coral_max,
  plot_coral_hot,
  ncol = 1,
  nrow = 5, 
  labels = c("A", "B", "C", "D", "E"), 
  common.legend = TRUE, 
  legend = "bottom",
  align = "hv") -> Figure_2_weather

Figure_2_weather

# ****  FIGURE 2 - weather data ---------------------------------------------------

library(gridExtra)
title1=text_grob("Sea Temperatures", size=10, face=2)
grid.arrange(Figure_2_weather, top=title1)


# **** MICROBIOME ANALYSIS **** ---------------------------------------------------

# Color Palettes ------------------------------------------------------------------
launch_evo_palette()
palette_box()

palette_box()$broccoli
palette_box()$fresh

pie(rep(1, length(palette_box()$broccoli)), labels = sprintf("%d (%s)",
        seq_along(palette_box()$broccoli),palette_box()$broccoli), col = palette_box()$broccoli)

pie(rep(1, length(palette_box()$fresh)), labels = sprintf("%d (%s)",
     seq_along(palette_box()$fresh),palette_box()$fresh), col = palette_box()$fresh)

palette_TX=c("#CBD588","#4AC6B1","#434080","#CB4A21", "#111111")
pie(rep(1, length(palette_TX)), labels = sprintf("%d (%s)",
    seq_along(palette_TX),palette_TX), col = palette_TX)


paletteCB4 = c("#E69F00", "#D55E00", "#CBD588", "#0072B2")
pie(rep(1, length(paletteCB4)), labels = sprintf("%d (%s)",
    seq_along(paletteCB4),paletteCB4), col = paletteCB4)


paletteCB2 = c("grey70", "grey17")
pie(rep(1, length(paletteCB2)), labels = sprintf("%d (%s)",
    seq_along(paletteCB2),paletteCB2), col = paletteCB2)

# ********************-------------------------------------------------------------
#> IMPORTING DATASETS -------------------------------------------------------------

# fungi
readRDS("datasets/fungi_final.rds") -> physeq_fun
physeq_fun

sample_data(physeq_fun) <- sample_data(physeq_fun)[,c("Species", "Status", "Description")]
sample_data(physeq_fun)

sample_data(physeq_fun)[sample_data(physeq_fun)$Description=="FJVB42A",]
sample_data(physeq_fun)[rownames(sample_data(physeq_fun))=="FJVB42A", 
                        colnames(sample_data(physeq_fun))=="Species"] <- "Porites_cylindrica"

# prkaryotes
readRDS("datasets/prok_final.rds") -> physeq_prok
physeq_prok

sample_data(physeq_prok) <- sample_data(physeq_prok)[,c("Species", "Status", "Description")]
sample_data(physeq_prok)
sample_data(physeq_prok)[sample_data(physeq_prok)$Description=="FJVB42A",]

# Symbiodinium
readRDS("datasets/symbio_final.rds") -> physeq_symb
physeq_symb

sample_data(physeq_symb) <- sample_data(physeq_symb)[,c("Species", "Status", "Description")]
sample_data(physeq_symb)
sample_data(physeq_symb)[sample_data(physeq_symb)$Description=="FJVB42A",]


# METADATA
write.csv(physeq_symb@sam_data, "map_symbio.csv")
saveRDS(object, file = "my_data.rds")


physeq_symb@sam_data
physeq_fun@sam_data
physeq_bact@sam_data

intersect(rownames(physeq_symb@sam_data),
          rownames(physeq_fun@sam_data),
          rownames(physeq_bact@sam_data)) 

# Refromat taxonomy ---------------------------------------------------------------
source("../R_functions/ReformatTaxonomy.R")

#Fungi
ReformatTaxonomy(physeq_fun) -> physeq_fun_filt
physeq_fun_filt
head(tax_table(physeq_fun_filt))

table(sample_data(physeq_fun_filt)$Status)

# Prok
taxa_prok <- as.data.frame(as.matrix(tax_table(physeq_prok)))
taxa_prok$OTU_ID <- rownames(taxa_prok)
taxa_prok <- taxa_prok[,c(1:7,10,9,8)] # modify this for different column number
identical(rownames(tax_table(physeq_prok)), rownames(taxa_prok))
tax_table(physeq_prok) <- NULL
tax_table(physeq_prok) <- as.matrix(taxa_prok)
head(tax_table(physeq_prok))

ReformatTaxonomy(physeq_prok) -> physeq_prok_filt
head(tax_table(physeq_prok_filt))

# Symb
taxa_symb <- as.data.frame(as.matrix(tax_table(physeq_symb)))
taxa_symb$OTU_ID <- rownames(taxa_symb)
taxa_symb <- taxa_symb[,c(1:7,10,9,8)] # modify this for different column number
taxa_symb[, "Subclade"] <- gsub("_Unclassified", "", taxa_symb[, "Subclade"])

identical(rownames(tax_table(physeq_symb)), rownames(taxa_symb))
tax_table(physeq_symb) <- NULL
tax_table(physeq_symb) <- as.matrix(taxa_symb)
head(tax_table(physeq_symb))

ReformatTaxonomy(physeq_symb) -> physeq_symb_filt
head(tax_table(physeq_symb_filt))

table(sample_data(physeq_symb_filt)$Species)
table(sample_data(physeq_symb_filt)$Status)

# ********************-------------------------------------------------------------
# #> RANDOM FOREST ----------------------------------------------------------------

# problems of random forest
# Ubalanced sample gorups
# Multicollinearity - OTUs that shows same abdundnace trends
# across samples. e.g. contaminants or real OTUs but governed by drift.
# Random forests are biased towards the categorical variable
# having multiple levels (categories)

#https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/#:~:text=A%20popular%20automatic%20method%20for,iteration%20to%20evaluate%20the%20model.
# Feature selection ---------------------------------------------------------------

# **** SYMBIODINIUM *****----------------------------------------------------------
otus_symb <- as.data.frame(otu_table(physeq_symb_filt))
metadata_symb = as(sample_data(physeq_symb_filt), "data.frame")
taxa_symb <- as.data.frame(as.matrix(tax_table(physeq_symb_filt)))
head(taxa_symb)

# centred log-ratio (CLR) transformation 
otus_symb_clr <- clr(otus_symb)
range(otus_symb_clr)

# Adding metadata
otus_symb_clr_Status <- data.frame(t(otus_symb_clr))
otus_symb_clr_Status$Sample <- rownames(otus_symb_clr_Status)
otus_symb_clr_Status$Status <-
  metadata_symb[rownames(otus_symb_clr_Status), "Status"]
head(otus_symb_clr_Status)

hist(
  colSums(otus_symb_clr),
  breaks = 100,
  col = "grey",
  xlab = "Read number",
  ylab = "Frequency",
  main = "Central Log Ratio"
)
dev.off()

# Recursive feature elimination 
rfe_symb_clr_Status <-
  recursive_feature_elimination(
    df = data.frame(otus_symb_clr_Status$Status, otus_symb_clr_Status[, 1:(ncol(otus_symb_clr_Status) - 2)]),
    sizes = seq(1, floor(nrow(otus_symb_clr)), by = 5),
    nfolds = 10,
    repeats = 5,
    fun = caret::rfFuncs,
    method_rfe = "cv",
    cpu_cores = 8)

rfe_symb_clr_Status
rfe_symb_clr_Status$optVariables
rfe_symb_clr_Status$resample

# plot 
ggplot(rfe_symb_clr_Status) + theme_bw() +
  geom_vline(xintercept = 126, linetype="dotted", color = "blue", size=1)

otus_symb_clr_Status[,rfe_symb_clr_Status$optVariables] -> otus_symb_clr_Status_sel
dim(otus_symb_clr_Status_sel)
head(otus_symb_clr_Status_sel)

identical(rownames(otus_symb_clr_Status_sel), rownames(otus_symb_clr_Status))
otus_symb_clr_Status_sel$Sample <- otus_symb_clr_Status$Sample
otus_symb_clr_Status_sel$Status <- otus_symb_clr_Status$Status

# try tuning the model first ------------------------------------------------------
floor(sqrt(ncol(otus_symb_clr_Status_sel) - 2))

set.seed(1)

bestmtry_symb_St <-
  tuneRF(
    x = otus_symb_clr_Status_sel[, 1:(ncol(otus_symb_clr_Status_sel) - 2)],
    y = otus_symb_clr_Status_sel$Status,
    stepFactor = 0.5,
    mtryStart = 11,
    improve = 0.0001,
    ntreeTry = 1001,
    doBest = FALSE,
    importance = TRUE
  )

bestmtry_symb_St

RF_symb_St <-
  randomForest(
    x = otus_symb_clr_Status_sel[, 1:(ncol(otus_symb_clr_Status_sel) - 2)],
    y = otus_symb_clr_Status_sel$Status,
    ntree = 1001,
    mtry = 11,
    nodesize = 1,
    importance = TRUE,
    proximity = TRUE,
    do.trace=200
  )

RF_symb_St
plot(RF_symb_St)

RF_symb_St$confusion
RF_symb_St$predicted

# removing unim-ortant variables using MeanDecreaseGini ---------------------------
head(importance(RF_symb_St))

##PlotMean Decrease Gini and Mean Decrease Accuracy 
varImpPlot(RF_symb_St) 

# ASSESSING MODEL FIT -------------------------------------------------------------
# Now check to see if the random forest is actually big enough... Up to a point, the more 
# trees in the forest, the better. You can tell when you've made enough when the OOB
# no longer improves.

oob_error_symb_st <- data.frame(
  Trees = rep(1:nrow(RF_symb_St$err.rate), times = 3),
  Status = rep(
    c("OOB", "bleached", "unbleached"),
    each = nrow(RF_symb_St$err.rate)),
  Error = c(
    RF_symb_St$err.rate[, "OOB"],
    RF_symb_St$err.rate[, "bleached"],
    RF_symb_St$err.rate[, "unbleached"]
  )
)

oob_error_symb_st$Status <- factor(oob_error_symb_st$Status,
                                   levels = c("bleached", "unbleached", "OOB"))

p_oob_error_symb_st <- ggplot(data=oob_error_symb_st, aes(x=Trees, y=Error)) +
  ggtitle("Model Errors") +
  theme_classic() +
  annotate("text", x=250, y=0.9, label="OOB estimate: 26.56%", size=3) +
  geom_line(aes(color=Status)) + 
  scale_color_manual(values = c(paletteCB2, "red")) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size =8)) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="bottom")

p_oob_error_symb_st

# ASSESS MODEL PERFORMANCE WITH PERMUTATIONS --------------------------------------
set.seed(2)

perm_RF_symb_St <-
  rf.significance(
    x = RF_symb_St,
    xdata = otus_symb_clr_Status_sel[, 1:(ncol(otus_symb_clr_Status_sel) - 2)],
    nperm = 999,
    nmtry = 11,
    ntree = 1001
  )

perm_RF_symb_St

# PLOTTING THE RESULTS ------------------------------------------------------------
# Now let's create an MDS-plot to show how the samples are related to each
# other. Start by converting the proximity matrix into a distance matrix.
dist_symb_st <- as.dist(1 - RF_symb_St$proximity)
mds_symb_st <- cmdscale(dist_symb_st, eig = TRUE, x.ret = TRUE)

## calculate the percentage of variation that each MDS axis accounts for...
var_symb_st <-
  round(mds_symb_st$eig / sum(mds_symb_st$eig) * 100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
values_symb_st <- mds_symb_st$points
mds_symb_st <- data.frame(
  Sample = rownames(values_symb_st),
  X = values_symb_st[, 1],
  Y = values_symb_st[, 2],
  Status = otus_symb_clr_Status_sel$Status
)

# Annotation position and hjust, vjust use to modify it
# xpos ypos        annotateText        hjust    vjust
# 1 -Inf -Inf Bottom Left (h0,v0)        0        0
# 2 -Inf  Inf    Top Left (h0,v1)        0        1
# 3  Inf -Inf  Bottom Right h1,v0        1        0
# 4  Inf  Inf     Top Right h1,v1        1        1

head(mds_symb_st)

predicted_symb_st <- as.data.frame(RF_symb_St$predicted)
colnames(predicted_symb_st) <- "Predicted"
head(predicted_symb_st)
identical(rownames(mds_symb_st),rownames(predicted_symb_st))
mds_symb_st$Predicted <- predicted_symb_st$Predicted

p_nmds_symb_st = ggplot(data=mds_symb_st, aes(x=X, y=Y, label=Sample, color=Status, shape=Predicted)) + 
  geom_point(size=1.2) +
  #geom_text_repel(label=rownames(mds_data), size = 2, min.segment.length = 0.5) +
  theme_classic() +
  xlab(paste("MDS1 - ", var_symb_st[1], "%", sep="")) +
  ylab(paste("MDS2 - ", var_symb_st[2], "%", sep="")) +
  ggtitle("Symbiodiniaceae - Status") +
  scale_color_manual(values = paletteCB2) +
  scale_shape_manual(values = c(16, 17)) +
  annotate("text", -Inf, -Inf,  hjust = - 0.2, vjust = - 1.2, 
           label="OOB = 21.88%", size=2.5) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) + # to adjust decimals
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8)) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=1, order = 1), 
         shape=guide_legend(ncol = 1, order = 2)) +
  theme(legend.position="bottom")

p_nmds_symb_st

# Identifying Important Features --------------------------------------------------
imp_symb_st <- as.data.frame(RF_symb_St$importance)
imp_symb_st$features <- rownames(imp_symb_st)
imp_symb_st <- arrange(imp_symb_st, desc(MeanDecreaseAccuracy))
head(imp_symb_st)

taxa_symb[rownames(imp_symb_st), ] -> taxa_symb_filt
head(taxa_symb_filt)

identical(rownames(taxa_symb_filt), rownames(imp_symb_st))
order_symb_st <-
  match(rownames(taxa_symb_filt), rownames(imp_symb_st))
imp_symb_st <- imp_symb_st[order_symb_st, ]

imp_symb_st$Taxonomy <- taxa_symb_filt$Taxonomy
imp_symb_st <-
  imp_symb_st[order(imp_symb_st$MeanDecreaseAccuracy, decreasing = TRUE),]
head(imp_symb_st)

p_bar_symb_st = ggplot(data=imp_symb_st[1:15,]) + 
  geom_bar(aes(x= reorder(Taxonomy, -MeanDecreaseAccuracy),
               y= MeanDecreaseAccuracy), color="lightblue", fill="lightblue",stat="identity") +
  labs(title = "Top Symbiodinium OTUs\nto classify coral Status", 
       x= "OTU", y= "Mean Decrease\nin Accuracy") +
  theme_classic() +
  coord_flip() +
  #ylim(0, 0.03) +
  scale_y_continuous(breaks=c(0.00, 0.01, 0.02, 0.03), limits = c(0, 0.03)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, size = 7 ,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 5, hjust = 1, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="none") #theme(legend.position=c(0.8, 0.8))

p_bar_symb_st

# ********************-------------------------------------------------------------
# > CASSIFYING BY SPECIES ---------------------------------------------------------

# Adding metadata
otus_symb_clr_species <- data.frame(t(otus_symb_clr))
otus_symb_clr_species$Sample <- rownames(otus_symb_clr_species)

metadata_symb_filt$Species <- metadata_symb_filt$Species %>%
  recode_factor(
    "Acropora_millepora" = "Acropora millepora",
    "Montipora_digitata" = "Montipora digitata",
    "Pocillopora_damicornis" = "Pocillopora damicornis",
    "Porites_cylindrica" = "Porites cylindrica"
  )

otus_symb_clr_species$Species <-
  metadata_symb_filt[rownames(otus_symb_clr_species), "Species"]
head(otus_symb_clr_species)

# Recursive feature elimination 
rfe_symb_clr_Species <-
  recursive_feature_elimination(
    df = data.frame(otus_symb_clr_species$Species, otus_symb_clr_species[, 1:(ncol(otus_symb_clr_species) - 2)]),
    sizes = seq(1, floor(nrow(otus_symb_clr)), by = 5),
    nfolds = 10,
    repeats = 5,
    fun = caret::rfFuncs,
    method_rfe = "cv",
    cpu_cores = 8)

rfe_symb_clr_Species
rfe_symb_clr_Species$optVariables

# plot 
ggplot(rfe_symb_clr_Species) + theme_bw() +
  geom_vline(xintercept = 26, linetype="dotted", color = "blue", size=1)

otus_symb_clr_species[,rfe_symb_clr_Species$optVariables] -> otus_symb_clr_species_sel
dim(otus_symb_clr_species_sel)
head(otus_symb_clr_species_sel)

identical(rownames(otus_symb_clr_species_sel), rownames(otus_symb_clr_species))
otus_symb_clr_species_sel$Sample <- otus_symb_clr_species$Sample
otus_symb_clr_species_sel$Species <- otus_symb_clr_species$Species

# try tuning the model first ------------------------------------------------------
floor(sqrt(ncol(otus_symb_clr_species_sel) - 2))

set.seed(4)

bestmtry_symb_sp <-
  tuneRF(
    x = otus_symb_clr_species_sel[, 1:(ncol(otus_symb_clr_species_sel) - 2)],
    y = otus_symb_clr_species_sel$Species,
    stepFactor = 1.5,
    mtryStart = 5,
    improve = 0.0001,
    ntreeTry = 1001,
    doBest = FALSE,
    importance = TRUE
  )

bestmtry_symb_sp

RF_symb_sp <-
  randomForest(
    x = otus_symb_clr_species_sel[, 1:(ncol(otus_symb_clr_species_sel) - 2)],
    y = otus_symb_clr_species_sel$Species,
    ntree = 1001,
    mtry = 5,
    nodesize = 1,
    importance = TRUE,
    proximity = TRUE,
    do.trace=100
  )

RF_symb_sp
RF_symb_sp$confusion
plot(RF_symb_sp)

# removing unim-ortant variables using MeanDecreaseGini ---------------------------
head(importance(RF_symb_sp))

##PlotMean Decrease Gini and Mean Decrease Accuracy 
varImpPlot(RF_symb_sp) 

# ASSESSING MODEL FIT -------------------------------------------------------------
oob_error_symb_sp <- data.frame(
  Trees = rep(1:nrow(RF_symb_sp$err.rate), times = 5),
  Species = rep(
    c(
      "OOB",
      "Acropora millepora",
      "Montipora digitata",
      "Pocillopora damicornis",
      "Porites cylindrica"
    ),
    each = nrow(RF_symb_sp$err.rate)
  ),
  Error = c(
    RF_symb_sp$err.rate[, "OOB"],
    RF_symb_sp$err.rate[, "Acropora millepora"],
    RF_symb_sp$err.rate[, "Montipora digitata"],
    RF_symb_sp$err.rate[, "Pocillopora damicornis"],
    RF_symb_sp$err.rate[, "Porites cylindrica"]
  )
)

oob_error_symb_sp$Species <- factor(
  oob_error_symb_sp$Species,
  levels = c(
    "Acropora millepora",
    "Montipora digitata",
    "Pocillopora damicornis",
    "Porites cylindrica",
    "OOB"
  )
)

p_oob_error_symb_sp <- ggplot(data=oob_error_symb_sp, aes(x=Trees, y=Error)) +
  ggtitle("Model Errors") +
  theme_classic() +
  annotate("text", x=250, y=0.9, label="OOB estimate: 1.56%", size=3) +
  geom_line(aes(color=Species)) + 
  scale_color_manual(values = c(paletteCB4, "black"), # to make legend text italic by level
      labels = c(expression(italic("Acropora millepora")),
                 expression(italic("Montipora digitata")),
                 expression(italic("Pocillopora damicornis")),
                 expression(italic("Porites cylindrica")),"OOB")) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), 
        legend.text = element_text(size = 7.5, hjust = 0)) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="bottom")

p_oob_error_symb_sp

# ASSESS MODEL PERFORMANCE WITH PERMUTATIONS --------------------------------------
set.seed(7)

perm_RF_symb_sp <-
  rf.significance(
    x = RF_symb_sp,
    xdata = otus_symb_clr_species_sel[, 1:(ncol(otus_symb_clr_species_sel) - 2)],
    nperm = 999,
    nmtry = 5,
    ntree = 1001
  )

perm_RF_symb_sp

# PLOTTING THE RESULTS ------------------------------------------------------------
# Now let's create an MDS-plot to show how the samples are related to each
# other. Start by converting the proximity matrix into a distance matrix.
dist_symb_sp <- as.dist(1 - RF_symb_sp$proximity)
mds_symb_sp <- cmdscale(dist_symb_sp, eig = TRUE, x.ret = TRUE)

#  calculate the percentage of variation that each MDS axis accounts for
var_symb_sp <- round(mds_symb_sp$eig / sum(mds_symb_sp$eig) * 100, 1)

# now make a fancy looking plot that shows the MDS axes and the variation:
values_symb_sp <- mds_symb_sp$points
mds_symb_sp <- data.frame(
  Sample = rownames(values_symb_sp),
  X = values_symb_sp[, 1],
  Y = values_symb_sp[, 2],
  Species = otus_symb_clr_species_sel$Species
)

# Annotation position and hjust, vjust use to modify it
# xpos ypos        annotateText        hjust    vjust
# 1 -Inf -Inf Bottom Left (h0,v0)        0        0
# 2 -Inf  Inf    Top Left (h0,v1)        0        1
# 3  Inf -Inf  Bottom Right h1,v0        1        0
# 4  Inf  Inf     Top Right h1,v1        1        1

head(mds_symb_sp)

predicted_symb_sp <- as.data.frame(RF_symb_sp$predicted)
colnames(predicted_symb_sp) <- "Predicted"
head(predicted_symb_sp)
identical(rownames(mds_symb_sp),rownames(predicted_symb_sp))
mds_symb_sp$Predicted <- predicted_symb_sp$Predicted


p_nmds_symb_sp = ggplot(data=mds_symb_sp, aes(x=X, y=Y, label=Sample, color=Species, shape=Predicted)) + 
  geom_point(size=1.2) +
  theme_classic() +
  xlab(paste("MDS1 - ", var_symb_sp[1], "%", sep="")) +
  ylab(paste("MDS2 - ", var_symb_sp[2], "%", sep="")) +
  ggtitle("Symbiodiniaceae - Species") +
  scale_color_manual(values = paletteCB4) +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  annotate("text", -Inf, -Inf,  hjust = - 0.2, vjust = - 1.2, 
           label="OOB = 0%", size=2.5) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) + # to adjust decimals
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 7.5, hjust = 0, face = "italic")) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=1, order = 1), 
         shape=guide_legend(ncol = 1, order = 2)) +
  theme(legend.position="bottom")

p_nmds_symb_sp

# Identifying Important Features --------------------------------------------------
imp_symb_sp <- as.data.frame(RF_symb_sp$importance)
imp_symb_sp$features <- rownames(imp_symb_sp)
imp_symb_sp <- arrange(imp_symb_sp, desc(MeanDecreaseAccuracy))
head(imp_symb_sp)

taxa_symb[rownames(imp_symb_sp), ] -> taxa_symb_filt
head(taxa_symb_filt)

identical(rownames(taxa_symb_filt), rownames(imp_symb_sp))
order_symb_sp <-
  match(rownames(taxa_symb_filt), rownames(imp_symb_sp))
imp_symb_sp <- imp_symb_sp[order_symb_sp,]

imp_symb_sp$Taxonomy <- taxa_symb_filt$Taxonomy
imp_symb_sp <-
  imp_symb_sp[order(imp_symb_sp$MeanDecreaseAccuracy, decreasing = TRUE),]
head(imp_symb_sp)

p_bar_symb_sp = ggplot(data=imp_symb_sp[1:15,]) + 
  geom_bar(aes(x= reorder(Taxonomy, -MeanDecreaseAccuracy),
               y= MeanDecreaseAccuracy), color="grey80", fill="grey80", stat="identity") +
  labs(title = "Top Symbiodinium OTUs\nto classify coral Species", 
       x= "OTU", y= "Mean Decrease\nin Accuracy") +
  theme_classic() +
  coord_flip() +
  #ylim(0, 0.06) +
  #scale_y_continuous(breaks=c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05), limits = c(0, 0.2)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, size = 7 ,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 5, hjust = 1, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="none") #theme(legend.position=c(0.8, 0.8))

p_bar_symb_sp


# ********************-------------------------------------------------------------
# **** FUNGI *****-----------------------------------------------------------------
otus_fungi <- as.data.frame(otu_table(physeq_fun_filt))
metadata_fungi = as(sample_data(physeq_fun_filt), "data.frame")
taxa_fungi <- as.data.frame(as.matrix(tax_table(physeq_fun_filt)))
head(taxa_fungi)

# centred log-ratio (CLR) transformation 
otus_fungi_clr <- clr(otus_fungi)
range(otus_fungi_clr)

# Adding metadata
otus_fungi_clr_Status <- data.frame(t(otus_fungi_clr))
otus_fungi_clr_Status$Sample <- rownames(otus_fungi_clr_Status)
otus_fungi_clr_Status$Status <-
  metadata_fungi[rownames(otus_fungi_clr_Status), "Status"]
head(otus_fungi_clr_Status)

hist(
  colSums(otus_symb_clr),
  breaks = 100,
  col = "grey",
  xlab = "Read number",
  ylab = "Frequency",
  main = "Central Log Ratio"
)
dev.off()

# Recursive feature elimination 
set.seed(10)

rfe_fungi_clr_Status <-
  recursive_feature_elimination(
    df = data.frame(otus_fungi_clr_Status$Status, otus_fungi_clr_Status[, 1:(ncol(otus_fungi_clr_Status) - 2)]),
    sizes = seq(1, floor(nrow(otus_fungi_clr)), by = 5),
    nfolds = 10,
    repeats = 5,
    fun = caret::rfFuncs,
    method_rfe = "cv", #cross-validation
    cpu_cores = 8)

rfe_fungi_clr_Status
rfe_fungi_clr_Status$optVariables

# plot 
ggplot(rfe_fungi_clr_Status) + theme_bw() +
  geom_vline(xintercept = 26, linetype="dotted", color = "blue", size=1)

otus_fungi_clr_Status[,rfe_fungi_clr_Status$optVariables] -> otus_fungi_clr_Status_sel
dim(otus_fungi_clr_Status_sel)
head(otus_fungi_clr_Status_sel)

identical(rownames(otus_fungi_clr_Status_sel), rownames(otus_fungi_clr_Status))
otus_fungi_clr_Status_sel$Sample <- otus_fungi_clr_Status$Sample
otus_fungi_clr_Status_sel$Status <- otus_fungi_clr_Status$Status

# > CASSIFYING BY STATUS ----------------------------------------------------------
# try tuning the model first ------------------------------------------------------
floor(sqrt(ncol(otus_fungi_clr_Status_sel) - 2))

set.seed(30)

bestmtry_fungi_St <-
  tuneRF(
    x = otus_fungi_clr_Status_sel[, 1:(ncol(otus_fungi_clr_Status_sel) - 2)],
    y = otus_fungi_clr_Status_sel$Status,
    stepFactor = 0.5,
    mtryStart = 4,
    improve = 0.0001,
    ntreeTry = 1001,
    doBest = FALSE,
    importance = TRUE
  )

bestmtry_fungi_St

set.seed(31)

RF_fungi_St <-
  randomForest(
    x = otus_fungi_clr_Status_sel[, 1:(ncol(otus_fungi_clr_Status_sel) - 2)],
    y = otus_fungi_clr_Status_sel$Status,
    ntree = 1001,
    mtry = 2,
    nodesize = 1,
    importance = TRUE,
    proximity = TRUE,
    do.trace=200
  )

RF_fungi_St
RF_fungi_St$confusion
plot(RF_fungi_St)

# removing unim-ortant variables using MeanDecreaseGini ---------------------------
head(importance(RF_fungi_St))

##PlotMean Decrease Gini and Mean Decrease Accuracy 
varImpPlot(RF_fungi_St) 

# ASSESSING MODEL FIT -------------------------------------------------------------
# Now check to see if the random forest is actually big enough... Up to a point, the more 
# trees in the forest, the better. You can tell when you've made enough when the OOB
# no longer improves.

oob_error_fungi_st <- data.frame(
  Trees = rep(1:nrow(RF_fungi_St$err.rate), times = 3),
  Status = rep(
    c("OOB", "bleached", "unbleached"),
    each = nrow(RF_fungi_St$err.rate)
  ),
  Error = c(
    RF_fungi_St$err.rate[, "OOB"],
    RF_fungi_St$err.rate[, "bleached"],
    RF_fungi_St$err.rate[, "unbleached"]
  )
)

oob_error_fungi_st$Status <- factor(oob_error_fungi_st$Status,
                                    levels = c("bleached", "unbleached", "OOB"))

p_oob_error_fungi_st <- ggplot(data=oob_error_fungi_st, aes(x=Trees, y=Error)) +
  ggtitle("Model Errors") +
  theme_classic() +
  annotate("text", x=250, y=0.9, label="OOB estimate: 20.0%", size=3) +
  geom_line(aes(color=Status)) + 
  scale_color_manual(values = c(paletteCB2, "red")) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size =8)) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="bottom")

p_oob_error_fungi_st

# ASSESS MODEL PERFORMANCE WITH PERMUTATIONS --------------------------------------
set.seed(33)
perm_RF_fungi_St <-
  rf.significance(
    x = RF_fungi_St,
    xdata = otus_fungi_clr_Status_sel[, 1:(ncol(otus_fungi_clr_Status_sel) - 2)],
    nperm = 999,
    nmtry = 2,
    ntree = 1001
  )

perm_RF_fungi_St

# PLOTTING THE RESULTS ------------------------------------------------------------
# Now let's create an MDS-plot to show how the samples are related to each 
# other. Start by converting the proximity matrix into a distance matrix.
dist_fungi_st <- as.dist(1 - RF_fungi_St$proximity)
mds_fungi_st <- cmdscale(dist_fungi_st, eig = TRUE, x.ret = TRUE)

## calculate the percentage of variation that each MDS axis accounts for...
var_fungi_st <- round(mds_fungi_st$eig / sum(mds_fungi_st$eig) * 100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
values_fungi_st <- mds_fungi_st$points

head(mds_fungi_st)

mds_fungi_st <- data.frame(
  Sample = rownames(values_fungi_st),
  X = values_fungi_st[, 1],
  Y = values_fungi_st[, 2],
  Status = otus_fungi_clr_Status_sel$Status
)

head(mds_fungi_st)

predicted_fungi_st <- as.data.frame(RF_fungi_St$predicted)
colnames(predicted_fungi_st) <- "Predicted"
head(predicted_fungi_st)
identical(rownames(mds_fungi_st),rownames(predicted_fungi_st))
mds_fungi_st$Predicted <- predicted_fungi_st$Predicted

# Annotation position and hjust, vjust use to modify it
# xpos ypos        annotateText        hjust    vjust
# 1 -Inf -Inf Bottom Left (h0,v0)        0        0
# 2 -Inf  Inf    Top Left (h0,v1)        0        1
# 3  Inf -Inf  Bottom Right h1,v0        1        0
# 4  Inf  Inf     Top Right h1,v1        1        1

p_nmds_fungi_st = ggplot(data=mds_fungi_st, aes(x=X, y=Y, label=Sample, color=Predicted, shape=Status)) + 
  geom_point(size=1.2) +
  #geom_text_repel(label=rownames(mds_data), size = 2, min.segment.length = 0.5) +
  theme_classic() +
  xlab(paste("MDS1 - ", var_fungi_st[1], "%", sep="")) +
  ylab(paste("MDS2 - ", var_fungi_st[2], "%", sep="")) +
  ggtitle("Fungi - Status") +
  scale_color_manual(values = paletteCB2) +
  scale_shape_manual(values = c(16, 17)) +
  annotate("text", -Inf, -Inf,  hjust = - 0.2, vjust = - 1.2, 
           label="OOB = 20.0%", size=2.5) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) + # to adjust decimals
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8)) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=1, order = 1), 
         shape=guide_legend(ncol = 1, order = 2)) +
  theme(legend.position="bottom")

p_nmds_fungi_st

# Identifying Important Features --------------------------------------------------
imp_fungi_st <- as.data.frame(RF_fungi_St$importance)
imp_fungi_st$features <- rownames(imp_fungi_st)
imp_fungi_st <- arrange(imp_fungi_st, desc(MeanDecreaseAccuracy))
head(imp_fungi_st)
dim(imp_fungi_st)

taxa_fungi[rownames(imp_fungi_st), ] -> taxa_fungi_filt
dim(taxa_fungi_filt)
head(taxa_fungi_filt)

identical(rownames(taxa_fungi_filt), rownames(imp_fungi_st))
order_fungi_st <-
  match(rownames(taxa_fungi_filt), rownames(imp_fungi_st))
imp_fungi_st <- imp_fungi_st[order_fungi_st,]

imp_fungi_st$Taxonomy <- taxa_fungi_filt$Taxonomy
imp_fungi_st <-
  imp_fungi_st[order(imp_fungi_st$MeanDecreaseAccuracy, decreasing = TRUE),]
head(imp_fungi_st)

p_bar_fungi_st = ggplot(data=imp_fungi_st) + 
  geom_bar(aes(x= reorder(Taxonomy, -MeanDecreaseAccuracy),
               y= MeanDecreaseAccuracy), color="indianred", fill="indianred",stat="identity") +
  labs(title = "Top Fungi OTUs to\nclassify coral Status", 
       x= "OTU", y= "Mean Decrease\nin Accuracy") +
  theme_classic() +
  coord_flip() +
  #ylim(0, 0.03) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, size = 7 ,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 5, hjust = 1, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="none") #theme(legend.position=c(0.8, 0.8))

p_bar_fungi_st


# ********************-------------------------------------------------------------
# > CASSIFYING BY SPECIES ---------------------------------------------------------
otus_fungi_clr_species <- data.frame(t(otus_fungi_clr))
otus_fungi_clr_species$Sample <- rownames(otus_fungi_clr_species)

metadata_fungi$Species <- metadata_fungi$Species %>%
  recode_factor(
    "Acropora_millepora" = "Acropora millepora",
    "Montipora_digitata" = "Montipora digitata",
    "Pocillopora_damicornis" = "Pocillopora damicornis",
    "Porites_cylindrica" = "Porites cylindrica"
  )

otus_fungi_clr_species$Species <-
  metadata_fungi[rownames(otus_fungi_clr_species), "Species"]
head(otus_fungi_clr_species)
dim(otus_fungi_clr_species)

# Recursive feature elimination 
set.seed(101)

rfe_fungi_clr_Species <-
  recursive_feature_elimination(
    df = data.frame(otus_fungi_clr_species$Species, otus_fungi_clr_species[, 1:(ncol(otus_fungi_clr_species) - 2)]),
    sizes = seq(1, floor(nrow(otus_fungi_clr)), by = 5),
    nfolds = 10,
    repeats = 5,
    fun = caret::rfFuncs,
    method_rfe = "cv",
    cpu_cores = 8)

rfe_fungi_clr_Species
rfe_fungi_clr_Species$optVariables

# plot 
ggplot(rfe_fungi_clr_Species) + theme_bw() +
  geom_vline(xintercept = 26, linetype="dotted", color = "blue", size=1)

otus_fungi_clr_species[,rfe_fungi_clr_Species$optVariables] -> otus_fungi_clr_species_sel
dim(otus_fungi_clr_species_sel)
head(otus_fungi_clr_species_sel)

identical(rownames(otus_fungi_clr_species_sel), rownames(otus_fungi_clr_species))
otus_fungi_clr_species_sel$Sample <- otus_fungi_clr_species$Sample
otus_fungi_clr_species_sel$Species <- otus_fungi_clr_species$Species

# try tuning the model first ------------------------------------------------------
floor(sqrt(ncol(otus_fungi_clr_species_sel) - 2))

set.seed(40)

bestmtry_fungi_sp <-
  tuneRF(
    x = otus_fungi_clr_species_sel[, 1:(ncol(otus_fungi_clr_species_sel) - 2)],
    y = otus_fungi_clr_species_sel$Species,
    stepFactor = 1.5,
    mtryStart = 3,
    improve = 0.0001,
    ntree = 3001,
    nodesize = 1,
    doBest = TRUE,
    importance = TRUE
  )

bestmtry_fungi_sp

set.seed(41)

RF_fungi_sp <-
  randomForest(
    x = otus_fungi_clr_species_sel[, 1:(ncol(otus_fungi_clr_species_sel) - 2)],
    y = otus_fungi_clr_species_sel$Species,
    ntree = 3001,
    mtry = 3,
    nodesize = 1,
    importance = TRUE,
    proximity = TRUE,
    do.trace=200
  )

RF_fungi_sp
plot(RF_fungi_sp)

# removing unim-ortant variables using MeanDecreaseGini ---------------------------
head(importance(RF_fungi_sp))

##PlotMean Decrease Gini and Mean Decrease Accuracy 
varImpPlot(RF_fungi_sp) 

# ASSESSING MODEL FIT -------------------------------------------------------------
oob_error_fungi_sp <- data.frame(
  Trees = rep(1:nrow(RF_fungi_sp$err.rate), times = 5),
  Species = rep(
    c(
      "OOB",
      "Acropora millepora",
      "Montipora digitata",
      "Pocillopora damicornis",
      "Porites cylindrica"
    ),
    each = nrow(RF_fungi_sp$err.rate)
  ),
  Error = c(
    RF_fungi_sp$err.rate[, "OOB"],
    RF_fungi_sp$err.rate[, "Acropora millepora"],
    RF_fungi_sp$err.rate[, "Montipora digitata"],
    RF_fungi_sp$err.rate[, "Pocillopora damicornis"],
    RF_fungi_sp$err.rate[, "Porites cylindrica"]
  )
)

oob_error_fungi_sp$Species <- factor(
  oob_error_fungi_sp$Species,
  levels = c(
    "Acropora millepora",
    "Montipora digitata",
    "Pocillopora damicornis",
    "Porites cylindrica",
    "OOB"
  )
)

p_oob_error_fungi_sp <- ggplot(data=oob_error_fungi_sp, aes(x=Trees, y=Error)) +
  ggtitle("Model Errors") +
  theme_classic() +
  annotate("text", x=250, y=0.9, label="OOB estimate: 43.08%", size=3) +
  geom_line(aes(color=Species)) + 
  scale_color_manual(values = c(paletteCB4, "black"), # to make legend text italic by level
                     labels = c(expression(italic("Acropora millepora")),
                                expression(italic("Montipora digitata")),
                                expression(italic("Pocillopora damicornis")),
                                expression(italic("Porites cylindrica")),"OOB")) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), 
        legend.text = element_text(size = 7.5, hjust = 0)) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="bottom")

p_oob_error_fungi_sp

# ASSESS MODEL PERFORMANCE WITH PERMUTATIONS --------------------------------------
set.seed(53)

perm_RF_fungi_sp <-
  rf.significance(
    x = RF_fungi_sp,
    xdata = otus_fungi_clr_species_sel[, 1:(ncol(otus_fungi_clr_species_sel) - 2)],
    nperm = 999,
    nmtry = 3,
    ntree = 3001
  )

perm_RF_fungi_sp

# PLOTTING THE RESULTS------------------------------------------------------------
# Now let's create an MDS-plot to show how the samples are related to each
# other. Start by converting the proximity matrix into a distance matrix.
dist_fungi_sp <- as.dist(1 - RF_fungi_sp$proximity)
mds_fungi_sp <- cmdscale(dist_fungi_sp, eig = TRUE, x.ret = TRUE)

## calculate the percentage of variation that each MDS axis accounts for...
var_fungi_sp <-
  round(mds_fungi_sp$eig / sum(mds_fungi_sp$eig) * 100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
values_fungi_sp <- mds_fungi_sp$points

mds_fungi_sp <- data.frame(
  Sample = rownames(values_fungi_sp),
  X = values_fungi_sp[, 1],
  Y = values_fungi_sp[, 2],
  Species = otus_fungi_clr_species_sel$Species
)

# Annotation position and hjust, vjust use to modify it
# xpos ypos        annotateText        hjust    vjust
# 1 -Inf -Inf Bottom Left (h0,v0)        0        0
# 2 -Inf  Inf    Top Left (h0,v1)        0        1
# 3  Inf -Inf  Bottom Right h1,v0        1        0
# 4  Inf  Inf     Top Right h1,v1        1        1

head(mds_fungi_sp)

predicted_fungi_sp <- as.data.frame(RF_fungi_sp$predicted)
colnames(predicted_fungi_sp) <- "Predicted"
head(predicted_fungi_sp)
identical(rownames(mds_fungi_sp),rownames(predicted_fungi_sp))
mds_fungi_sp$Predicted <- predicted_fungi_sp$Predicted


p_nmds_fungi_sp = ggplot(data=mds_fungi_sp, aes(x=X, y=Y, label=Sample, color=Species, shape = Predicted)) + 
  geom_point(size=1.2) +
  theme_classic() +
  xlab(paste("MDS1 - ", var_fungi_sp[1], "%", sep="")) +
  ylab(paste("MDS2 - ", var_fungi_sp[2], "%", sep="")) +
  ggtitle("Fungi - Species") +
  scale_color_manual(values = paletteCB4) +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  annotate("text", -Inf, -Inf,  hjust = - 0.2, vjust = - 1.2, 
           label="OOB = 43.08%", size=2.5) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) + # to adjust decimals
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 7.5, hjust = 0, face = "italic")) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=1, order = 1), 
         shape=guide_legend(ncol = 1, order = 2)) +
  theme(legend.position="bottom")

p_nmds_fungi_sp

# Identifying Important Features --------------------------------------------------
imp_fungi_sp <- as.data.frame(RF_fungi_sp$importance)
imp_fungi_sp$features <- rownames(imp_fungi_sp)
imp_fungi_sp <- arrange(imp_fungi_sp, desc(MeanDecreaseAccuracy))
head(imp_fungi_sp)

taxa_fungi[rownames(imp_fungi_sp), ] -> taxa_fungi_filt
dim(taxa_fungi_filt)
head(taxa_fungi_filt)

identical(rownames(taxa_fungi_filt), rownames(imp_fungi_sp))
order_fungi_sp <-
  match(rownames(taxa_fungi_filt), rownames(imp_fungi_sp))
imp_fungi_sp <- imp_fungi_sp[order_fungi_sp, ]

imp_fungi_sp$Taxonomy <- taxa_fungi_filt$Taxonomy
imp_fungi_sp <-
  imp_fungi_sp[order(imp_fungi_sp$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(imp_fungi_sp)

p_bar_fungi_sp = ggplot(data=imp_fungi_sp) + 
  geom_bar(aes(x= reorder(Taxonomy, -MeanDecreaseAccuracy),
               y= MeanDecreaseAccuracy), color="indianred", fill="indianred",stat="identity") +
  labs(title = "Top Fungi OTUs to\nclassify coral Species", 
       x= "OTU", y= "Mean Decrease\nin Accuracy") +
  theme_classic() +
  coord_flip() +
  #ylim(0, 0.06) +
  scale_y_continuous(breaks=c(0.00, 0.01, 0.02, 0.03), limits = c(0, 0.03)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, size = 7 ,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 5, hjust = 1, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="none") #theme(legend.position=c(0.8, 0.8))

p_bar_fungi_sp

# ********************-------------------------------------------------------------
# **** BACTERIA *****--------------------------------------------------------------
otus_bact <- as.data.frame(otu_table(physeq_prok_filt))
metadata_bact = as(sample_data(physeq_prok_filt), "data.frame")
taxa_bact <- as.data.frame(as.matrix(tax_table(physeq_prok_filt)))
head(taxa_bact)

# centred log-ratio (CLR) transformation 
otus_bact_clr <- clr(otus_bact)
range(otus_bact_clr)

# > CASSIFYING BY STATUS ----------------------------------------------------------
# Adding metadata
otus_bact_clr_Status <- data.frame(t(otus_bact_clr))
otus_bact_clr_Status$Sample <- rownames(otus_bact_clr_Status)
otus_bact_clr_Status$Status <-
  metadata_bact[rownames(otus_bact_clr_Status), "Status"]
head(otus_bact_clr_Status)

hist(
  colSums(otus_bact_clr),
  breaks = 100,
  col = "grey",
  xlab = "Read number",
  ylab = "Frequency",
  main = "Central Log Ratio"
)
dev.off()

# Recursive feature elimination
set.seed(201)

rfe_bact_clr_Status <-
  recursive_feature_elimination(
    df = data.frame(otus_bact_clr_Status$Status, otus_bact_clr_Status[, 1:(ncol(otus_bact_clr_Status) - 2)]),
    sizes = seq(1, floor(nrow(otus_bact_clr)), by = 5),
    nfolds = 10,
    repeats = 5,
    fun = caret::rfFuncs,
    method_rfe = "cv",
    cpu_cores = 8)

rfe_bact_clr_Status
rfe_bact_clr_Status$optVariables
rfe_bact_clr_Status$resample

# plot 
ggplot(rfe_bact_clr_Status) + theme_bw()
  geom_vline(xintercept = 571, linetype="dotted", color = "blue", size=1)

otus_bact_clr_Status[,rfe_bact_clr_Status$optVariables] -> otus_bact_clr_Status_sel
dim(otus_bact_clr_Status_sel)
head(otus_bact_clr_Status_sel)

identical(rownames(otus_bact_clr_Status_sel), rownames(otus_bact_clr_Status))
otus_bact_clr_Status_sel$Sample <- otus_bact_clr_Status$Sample
otus_bact_clr_Status_sel$Status <- otus_bact_clr_Status$Status

# try tuning the model first ------------------------------------------------------
floor(sqrt(ncol(otus_bact_clr_Status_sel) - 2))

set.seed(105)

bestmtry_bact_St <-
  tuneRF(
    x = otus_bact_clr_Status[, 1:(ncol(otus_bact_clr_Status) - 2)],
    y = otus_bact_clr_Status$Status,
    stepFactor = 0.5,
    mtryStart = 23,
    improve = 0.0001,
    ntreeTry = 3001,
    doBest = FALSE,
    importance = TRUE
  )

bestmtry_bact_St

set.seed(51)

RF_bact_St <-
  randomForest(
    x = otus_bact_clr_Status[, 1:(ncol(otus_bact_clr_Status) - 2)],
    y = otus_bact_clr_Status$Status,
    ntree = 1001,
    mtry = 5,
    nodesize = 1,
    importance = TRUE,
    proximity = TRUE,
    do.trace=200
  )

RF_bact_St
plot(RF_bact_St)

# removing unim-ortant variables using MeanDecreaseGini ---------------------------
head(importance(RF_bact_St))

##PlotMean Decrease Gini and Mean Decrease Accuracy
varImpPlot(RF_bact_St)

# ASSESSING MODEL FIT -------------------------------------------------------------
# Now check to see if the random forest is actually big enough... Up to a point, the more
# trees in the forest, the better. You can tell when you've made enough when the OOB
# no longer improves.

oob_error_bact_st <- data.frame(
  Trees = rep(1:nrow(RF_bact_St$err.rate), times = 3),
  Status = rep(
    c("OOB", "bleached", "unbleached"),
    each = nrow(RF_bact_St$err.rate)
  ),
  Error = c(
    RF_bact_St$err.rate[, "OOB"],
    RF_bact_St$err.rate[, "bleached"],
    RF_bact_St$err.rate[, "unbleached"]
  )
)

oob_error_bact_st$Status <- factor(oob_error_bact_st$Status,
                                   levels = c("bleached", "unbleached", "OOB"))

p_oob_error_bact_st <- ggplot(data=oob_error_bact_st, aes(x=Trees, y=Error)) +
  ggtitle("Model Errors") +
  theme_classic() +
  annotate("text", x=250, y=0.9, label="OOB estimate: 43.75%", size=3) +
  geom_line(aes(color=Status)) + 
  scale_color_manual(values = c(paletteCB2, "red")) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size =8)) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="bottom")

p_oob_error_bact_st

# ASSESS MODEL PERFORMANCE WITH PERMUTATIONS --------------------------------------
set.seed(63)

perm_RF_bact_St <-
  rf.significance(
    x = RF_bact_St,
    xdata = otus_bact_clr_Status[, 1:(ncol(otus_bact_clr_Status) - 2)],
    nperm = 999,
    nmtry = 5,
    ntree = 1001
  )

perm_RF_bact_St # Model not significant at p=0.124


# PLOTTING THE RESULTS ------------------------------------------------------------
# Now let's create an MDS-plot to show how the samples are related to each
# other. Start by converting the proximity matrix into a distance matrix.
dist_bact_st <- as.dist(1 - RF_bact_St$proximity)
mds_bact_st <- cmdscale(dist_bact_st, eig = TRUE, x.ret = TRUE)

## calculate the percentage of variation that each MDS axis accounts for...
var_bact_st <- round(mds_bact_st$eig / sum(mds_bact_st$eig) * 100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
values_bact_st <- mds_bact_st$points
mds_bact_st <- data.frame(
  Sample = rownames(values_bact_st),
  X = values_bact_st[, 1],
  Y = values_bact_st[, 2],
  Status = otus_bact_clr_Status$Status
)

# Annotation position and hjust, vjust use to modify it
# xpos ypos        annotateText        hjust    vjust
# 1 -Inf -Inf Bottom Left (h0,v0)        0        0
# 2 -Inf  Inf    Top Left (h0,v1)        0        1
# 3  Inf -Inf  Bottom Right h1,v0        1        0
# 4  Inf  Inf     Top Right h1,v1        1        1

head(mds_bact_st)

predicted_bact_st <- as.data.frame(RF_bact_St$predicted)
colnames(predicted_bact_st) <- "Predicted"
head(predicted_bact_st)
identical(rownames(mds_bact_st),rownames(predicted_bact_st))
mds_bact_st$Predicted <- predicted_bact_st$Predicted

p_nmds_bact_st = ggplot(data=mds_bact_st, aes(x=X, y=Y, label=Sample, color=Status, shape=Predicted)) + 
  geom_point(size=1.2) +
  #geom_text_repel(label=rownames(mds_data), size = 2, min.segment.length = 0.5) +
  theme_classic() +
  xlab(paste("MDS1 - ", var_bact_st[1], "%", sep="")) +
  ylab(paste("MDS2 - ", var_bact_st[2], "%", sep="")) +
  ggtitle("Bacteria - Status") +
  scale_color_manual(values = paletteCB2) +
  scale_shape_manual(values = c(16, 17)) +
  annotate("text", -Inf, -Inf,  hjust = - 0.2, vjust = - 1.2, 
           label="OOB = 43.75%", size=2.5) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) + # to adjust decimals
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8)) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=1, order = 1), 
         shape=guide_legend(ncol =1, order = 2)) +
  theme(legend.position="bottom")

p_nmds_bact_st

# Identifying Important Features --------------------------------------------------
imp_bact_st <- as.data.frame(RF_bact_St$importance)
imp_bact_st$features <- rownames(imp_bact_st)
imp_bact_st <- arrange(imp_bact_st, desc(MeanDecreaseAccuracy))
head(imp_bact_st)
dim(imp_bact_st)

taxa_bact[rownames(imp_bact_st),] -> taxa_bact_filt
dim(taxa_bact_filt)
head(taxa_bact_filt)

identical(rownames(taxa_bact_filt), rownames(imp_bact_st))
order_bact_st <-
  match(rownames(taxa_bact_filt), rownames(imp_bact_st))
imp_bact_st <- imp_bact_st[order_bact_st, ]

imp_bact_st$Taxonomy <- taxa_bact_filt$Taxonomy
imp_bact_st <-
  imp_bact_st[order(imp_bact_st$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(imp_bact_st)

p_bar_bact_st = ggplot(data=imp_bact_st[1:20,]) + 
  geom_bar(aes(x= reorder(Taxonomy, -MeanDecreaseAccuracy),
               y= MeanDecreaseAccuracy), color="indianred", fill="indianred",stat="identity") +
  labs(title = "Top Bacteria OTUs to\nclassify coral Status", 
       x= "OTU", y= "Mean Decrease\nin Accuracy") +
  theme_classic() +
  coord_flip() +
  #ylim(0, 0.03) +
  scale_y_continuous(breaks=c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05), limits = c(0, 0.05)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, size = 7 ,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 5, hjust = 1, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="none") #theme(legend.position=c(0.8, 0.8))

p_bar_bact_st


# ********************-------------------------------------------------------------
# > CASSIFYING BY SPECIES ---------------------------------------------------------
otus_bact_clr_species <- data.frame(t(otus_bact_clr))
otus_bact_clr_species$Sample <- rownames(otus_bact_clr_species)

metadata_bact$Species <- metadata_bact$Species %>%
  recode_factor(
    "Acropora_millepora" = "Acropora millepora",
    "Montipora_digitata" = "Montipora digitata",
    "Pocillopora_damicornis" = "Pocillopora damicornis",
    "Porites_cylindrica" = "Porites cylindrica"
  )

otus_bact_clr_species$Species <-
  metadata_bact[rownames(otus_bact_clr_species), "Species"]
head(otus_bact_clr_species)
dim(otus_bact_clr_species)

# Recursive feature elimination
set.seed(301)

rfe_bact_clr_Species <-
  recursive_feature_elimination(
    df = data.frame(otus_bact_clr_species$Species, otus_bact_clr_species[, 1:(ncol(otus_bact_clr_species) - 2)]),
    sizes = seq(1, floor(nrow(otus_bact_clr)), by = 5),
    nfolds = 10,
    repeats = 5,
    fun = caret::rfFuncs,
    method_rfe = "cv",
    cpu_cores = 8)

rfe_bact_clr_Species
rfe_bact_clr_Species$optVariables
rfe_bact_clr_Species$resample

# plot 
ggplot(rfe_bact_clr_Species) + theme_bw() +
geom_vline(xintercept = 571, linetype="dotted", color = "blue", size=1)

otus_bact_clr_species[,rfe_bact_clr_Species$optVariables] -> otus_bact_clr_species_sel
dim(otus_bact_clr_species_sel)
head(otus_bact_clr_species_sel)

identical(rownames(otus_bact_clr_species_sel), rownames(otus_bact_clr_species))
otus_bact_clr_species_sel$Sample <- otus_bact_clr_species$Sample
otus_bact_clr_species_sel$Species <- otus_bact_clr_species$Species

# try tuning the model first ------------------------------------------------------
floor(sqrt(ncol(otus_bact_clr_species_sel) - 2))

set.seed(401)

bestmtry_bact_sp <-
  tuneRF(
    x = otus_bact_clr_species_sel[, 1:(ncol(otus_bact_clr_species_sel) - 2)],
    y = otus_bact_clr_species_sel$Species,
    stepFactor = 0.5,
    mtryStart = 23,
    improve = 0.0001,
    ntreeTry = 3001,
    doBest = FALSE,
    importance = TRUE
  )

bestmtry_bact_sp

set.seed(402)

RF_bact_sp <-
  randomForest(
    x = otus_bact_clr_species_sel[, 1:(ncol(otus_bact_clr_species_sel) - 2)],
    y = otus_bact_clr_species_sel$Species,
    ntree = 3001,
    mtry = 11,
    nodesize = 1,
    importance = TRUE,
    proximity = TRUE,
    do.trace=200
  )

RF_bact_sp
plot(RF_bact_sp)

# removing unim-ortant variables using MeanDecreaseGini ---------------------------
head(importance(RF_bact_sp))

##PlotMean Decrease Gini and Mean Decrease Accuracy
varImpPlot(RF_bact_sp)

# ASSESSING MODEL FIT -------------------------------------------------------------
oob_error_bact_sp <- data.frame(
  Trees = rep(1:nrow(RF_bact_sp$err.rate), times = 5),
  Species = rep(
    c(
      "OOB",
      "Acropora millepora",
      "Montipora digitata",
      "Pocillopora damicornis",
      "Porites cylindrica"
    ),
    each = nrow(RF_bact_sp$err.rate)
  ),
  Error = c(
    RF_bact_sp$err.rate[, "OOB"],
    RF_bact_sp$err.rate[, "Acropora millepora"],
    RF_bact_sp$err.rate[, "Montipora digitata"],
    RF_bact_sp$err.rate[, "Pocillopora damicornis"],
    RF_bact_sp$err.rate[, "Porites cylindrica"]
  )
)

oob_error_bact_sp$Species <- factor(
  oob_error_bact_sp$Species,
  levels = c(
    "Acropora millepora",
    "Montipora digitata",
    "Pocillopora damicornis",
    "Porites cylindrica",
    "OOB"
  )
)

p_oob_error_bact_sp <- ggplot(data=oob_error_bact_sp, aes(x=Trees, y=Error)) +
  ggtitle("Model Errors") +
  theme_classic() +
  annotate("text", x=250, y=0.9, label="OOB estimate:  21.05%", size=3) +
  geom_line(aes(color=Species)) + 
  scale_color_manual(values = c(paletteCB4, "black"), # to make legend text italic by level
                     labels = c(expression(italic("Acropora millepora")),
                                expression(italic("Montipora digitata")),
                                expression(italic("Pocillopora damicornis")),
                                expression(italic("Porites cylindrica")),"OOB")) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), 
        legend.text = element_text(size = 7.5, hjust = 0)) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="bottom")

p_oob_error_bact_sp

# ASSESS MODEL PERFORMANCE WITH PERMUTATIONS --------------------------------------
set.seed(405)

perm_RF_bact_sp <-
  rf.significance(
    x = RF_bact_sp,
    xdata = otus_bact_clr_species_sel[, 1:(ncol(otus_bact_clr_species_sel) - 2)],
    nperm = 999,
    nmtry = 11,
    ntree = 3001
  )

perm_RF_bact_sp

# PLOTTING THE RESULTS ------------------------------------------------------------
# Now let's create an MDS-plot to show how the samples are related to each
# other. Start by converting the proximity matrix into a distance matrix.
dist_bact_sp <- as.dist(1 - RF_bact_sp$proximity)
mds_bact_sp <- cmdscale(dist_bact_sp, eig = TRUE, x.ret = TRUE)

## calculate the percentage of variation that each MDS axis accounts for...
var_bact_sp <- round(mds_bact_sp$eig / sum(mds_bact_sp$eig) * 100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
values_bact_sp <- mds_bact_sp$points

mds_bact_sp <- data.frame(
  Sample = rownames(values_bact_sp),
  X = values_bact_sp[, 1],
  Y = values_bact_sp[, 2],
  Species = otus_bact_clr_species_sel$Species
)

head(mds_bact_sp)

predicted_fungi_sp <- as.data.frame(RF_bact_sp$predicted)
colnames(predicted_fungi_sp) <- "Predicted"
head(predicted_fungi_sp)
identical(rownames(mds_bact_sp),rownames(predicted_fungi_sp))
mds_bact_sp$Predicted <- predicted_fungi_sp$Predicted


# Annotation position and hjust, vjust use to modify it
# xpos ypos        annotateText        hjust    vjust
# 1 -Inf -Inf Bottom Left (h0,v0)        0        0
# 2 -Inf  Inf    Top Left (h0,v1)        0        1
# 3  Inf -Inf  Bottom Right h1,v0        1        0
# 4  Inf  Inf     Top Right h1,v1        1        1

p_nmds_bact_sp = ggplot(data=mds_bact_sp, aes(x=X, y=Y, label=Sample, color=Species, shape=Predicted)) + 
  geom_point(size=1.2) +
  theme_classic() +
  xlab(paste("MDS1 - ", var_bact_sp[1], "%", sep="")) +
  ylab(paste("MDS2 - ", var_bact_sp[2], "%", sep="")) +
  ggtitle("Bacteria - Species") +
  scale_color_manual(values = paletteCB4) +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  annotate("text", -Inf, -Inf,  hjust = - 0.2, vjust = - 1.2, 
           label="OOB = 21.05%", size=2.5) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) + # to adjust decimals
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 7.5, hjust = 0, face = "italic")) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=1, order = 1), 
         shape=guide_legend(ncol = 1, order = 2)) +
  theme(legend.position="bottom")

p_nmds_bact_sp

# Identifying Important Features --------------------------------------------------
imp_bact_sp <- as.data.frame(RF_bact_sp$importance)
imp_bact_sp$features <- rownames(imp_bact_sp)
imp_bact_sp <- arrange(imp_bact_sp, desc(MeanDecreaseAccuracy))
head(imp_bact_sp)

taxa_bact[rownames(imp_bact_sp),] -> taxa_bact_filt
dim(taxa_bact_filt)
head(taxa_bact_filt)

identical(rownames(taxa_bact_filt), rownames(imp_bact_sp))
order_bact_sp <-
  match(rownames(taxa_bact_filt), rownames(imp_bact_sp))
imp_bact_sp <- imp_bact_sp[order_bact_sp, ]

imp_bact_sp$Taxonomy <- taxa_bact_filt$Taxonomy
imp_bact_sp <-
  imp_bact_sp[order(imp_bact_sp$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(imp_bact_sp)

p_bar_bact_sp = ggplot(data=imp_bact_sp[1:15,]) + 
  #geom_bar(aes(x= reorder(Taxonomy, -MeanDecreaseAccuracy),
  #             y= MeanDecreaseAccuracy), color="indianred", fill="indianred",stat="identity") +
  geom_col(aes(x= reorder(Taxonomy, -MeanDecreaseAccuracy),
               y= MeanDecreaseAccuracy), fill = "grey80", position = position_dodge2(width = 0.5, preserve = "single")) +
  labs(title = "Top Bacteria OTUs to\nclassify coral Species", 
       x= "OTU", y= "Mean Decrease\nin Accuracy") +
  theme_classic() +
  coord_flip() +
  #ylim(0, 0.06) +
  scale_y_continuous(breaks=c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05), limits = c(0, 0.05)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, size = 7 ,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 5, hjust = 1, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="none") #theme(legend.position=c(0.8, 0.8))

p_bar_bact_sp

# FIGURE X - A MDS plots proximity matrix -----------------------------------------
ggarrange(
  ggarrange(p_nmds_symb_sp, p_nmds_bact_sp,p_nmds_fungi_sp,
            labels = c("A", "B", "C"),
            widths = c(1, 1, 1),
            heights = c(1, 1, 1),
            align = "hv" ,
            ncol = 3,
            nrow = 1,
            common.legend = TRUE,
            legend = c("right")),
  ggarrange(p_nmds_symb_st, p_nmds_fungi_st,
            labels = c("D", "E"),
            widths = c(1, 1),
            heights = c(1, 1),
            align = "hv" ,
            ncol = 2,
            nrow = 1,
            common.legend = TRUE,
            legend = c("right")),
  widths = c(1, 0.66),
  align = "h",
  ncol = 2,
  nrow = 1,
  legend = c("right")) -> Fig_RF_mds

Fig_RF_mds


# INDICATOR SPECIES ANALYSIS ------------------------------------------------------
library(indicspecies)

GetIndicators <-function(dataframe, var){
  require(phyloseq); require(indicspecies); require(dplyr)
  otu <- as.data.frame(otu_table(dataframe))
  metadata = as(sample_data(dataframe), "data.frame")
  taxa <- as.data.frame(as.matrix(tax_table(dataframe)))
  multipatt <- multipatt(t(otu), metadata[,var], func = "IndVal.g",
                         control=how(nperm=9999), duleg=TRUE)
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


# indicator value >0.5 and p-value < 0.05 after fdr correction
GetIndicators(physeq_fun_filt, "Status") -> ind_fungi_st
ind_fungi_st
GetIndicators(physeq_fun_filt, "Species") -> ind_fungi_sp
ind_fungi_sp

GetIndicators(physeq_symb_filt, "Status") -> ind_symb_st
ind_symb_st
GetIndicators(physeq_symb_filt, "Species") -> ind_symb_sp
ind_symb_sp

GetIndicators(physeq_prok_filt, "Status") -> ind_bact_st
ind_bact_st
GetIndicators(physeq_prok_filt, "Species") -> ind_bact_sp
ind_bact_sp

# RF AND INDICATORS MATCH ---------------------------------------------------------
InditoRF <- function(df_imp, df_r) {
  #subset(ind_ITS_mark, ind_ITS_mark$stat >= 0.5) -> ind_ITS_mark_filt
  df_imp[1:15, ] -> df_imp20
  df_imp20$index <- df_r[rownames(df_imp20),]$index
  df_imp20$stat <- df_r[rownames(df_imp20),]$stat
  df_imp20$index <- df_imp20$index %>%
      recode_factor(
        "1" = "Acropora millepora",
        "2" = "Montipora digitata",
        "3" = "Pocillopora damicornis",
        "4" = "Porites cylindrica")
  return(df_imp20)
}

# Match indicators with RF models features accuracy -------------------------------
InditoRF(df_imp = imp_symb_sp, df_r = ind_symb_sp) -> ind_RF_symb_sp
ind_RF_symb_sp

InditoRF(df_imp = imp_bact_sp, df_r = ind_bact_sp) -> ind_RF_bact_sp
ind_RF_bact_sp

PlotIndtoRF <- function(df){
  require(purrr)
  pal_market = list(values=paletteCB4, na.value="grey80")
  ggplot(data=df, aes(x= reorder(Taxonomy,-MeanDecreaseAccuracy),
                      y= MeanDecreaseAccuracy)) -> plot 
    if (any(colnames(df) == "index")==TRUE){
      plot + geom_bar(aes(color=index, fill=index), stat="identity") +
        geom_text(aes(label=round(stat, digits = 2)), fontface  = "bold", angle = 0,
                  position = position_dodge(1), hjust= -0.2, color="black", size=1.8) +
        invoke(scale_colour_manual, pal_market) + invoke(scale_fill_manual, pal_market) +
        guides(color=guide_legend(title="IndVal.g. value",ncol=1), 
               fill=guide_legend(title="IndVal.g. value", ncol=1)) -> plot1
    }else{
      plot + geom_bar(color="grey80", fill="grey80", stat="identity") -> plot1
    }
    plot1 + theme_classic() +
    coord_flip() +
    #ylim(0, 0.2) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)) +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 45, size = 6 ,hjust = 1, vjust = 1)) +
    theme(axis.text.y = element_text(angle = 0, size = 6, hjust = 1, vjust = 0.5)) +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_text(angle = 90, size = 8, face = "bold")) +
    theme(legend.position="right") -> plot_IND
  return(plot_IND)
}


PlotIndtoRF(ind_RF_symb_sp[1:15,]) + labs(title = "Top ITS OTUs to\nclassify markets",  x= "OTU", y= "Mean Decrease\nin Accuracy") +ylim(0, 0.19)
PlotIndtoRF(ind_RF_bact_sp[1:15,]) + labs(title = "Top ITS OTUs to\nclassify niches",  x= "OTU", y= "Mean Decrease\nin Accuracy") + 
  ylim(0, 0.027) + 
  theme(axis.title.y = element_blank())
PlotIndtoRF(imp_fungi_sp[1:11,]) + labs(title = "Top ITS OTUs to\nclassify niches",  x= "OTU", y= "Mean Decrease\nin Accuracy")
PlotIndtoRF(imp_symb_st[1:15,]) + labs(title = "Top ITS OTUs to\nclassify niches",  x= "OTU", y= "Mean Decrease\nin Accuracy")
PlotIndtoRF(imp_fungi_st[1:15,]) + labs(title = "Top ITS OTUs to\nclassify niches",  x= "OTU", y= "Mean Decrease\nin Accuracy")

  # scale_y_discrete(labels = c("Zotu39-Aspergillus sp.",
  #                             "Zotu2-Eurotiales",
  #                             "Zotu1625-Candida sp.",
  #                             "Zotu44-Trichoderma sp.",
  #                             "Zotu2107-Cladosporium sp.",
  #                             "Zotu1293-Dothideales",          
  #                             "Zotu748-Rhodosporidiobolus sp.",
  #                             "Zotu85-Parengyodontium sp.",
  #                             "Zotu54-Parengyodontium sp.",
  #                             "Zotu13-Eurotiales",
  #                             "Zotu3882-Candida sp.",
  #                             "",
  #                             "",
  #                             "",
  #                             ""))


# FIGURE X - B Top OTUs for classification ---------------------------------
ggarrange(PlotIndtoRF(ind_RF_symb_sp[1:15, ]) + labs(title = NULL,  x = NULL, y = "Mean Decrease in Accuracy")+ 
            ylim(0, 0.19),
          PlotIndtoRF(ind_RF_bact_sp[1:15, ]) + labs(title = NULL,  x = NULL, y = NULL) + 
            ylim(0, 0.027) + 
            theme(axis.title.y = element_blank()),
          PlotIndtoRF(imp_fungi_sp[1:11, ]) + labs(title = NULL,  x = NULL, y = NULL) + 
            theme(axis.title.y = element_blank()),
          PlotIndtoRF(imp_symb_st[1:15, ]) + labs(title = NULL,  x = NULL, y = NULL) + 
            theme(axis.title.y = element_blank()),
          PlotIndtoRF(imp_fungi_st[1:15, ]) + labs(title = NULL,  x = NULL, y = NULL) + 
            theme(axis.title.y = element_blank()), 
          labels = c("F","G","H","I","J"), 
          widths = c(1, 1, 1, 1, 1),
           # heights = c(1, 1, 1, 1, 1),
            align = "hv" ,
            ncol = 5,
            nrow = 1,
            legend = c("none")) -> Fig_RF_bars

Fig_RF_bars


ggarrange(Fig_RF_mds, 
          Fig_RF_bars,
          widths = c(1, 1),
          heights = c(1, 1),
          align = "hv" ,
          ncol = 1,
          nrow = 2) -> Fig_RF_all

Fig_RF_all

# ********************-------------------------------------------------------------
# *** FIGURE 6 - Ranfom Forest Model ----------------------------------------------

ggarrange(ggarrange(ggarrange(p_nmds_symb_st,
                              p_nmds_symb_sp,
                              labels = c("A","B"),
                              align = "hv" ,
                              widths = c(1,1),
                              ncol = 2, 
                              nrow = 1, 
                              common.legend = TRUE,
                              legend = c("bottom")),
                    ggarrange(p_nmds_fungi_sp,
                              p_nmds_bact_sp,
                              labels = c("C", "D"), 
                              align = "hv" ,
                              widths = c(1,1),
                              ncol = 2, 
                              nrow = 1, 
                              common.legend = TRUE,
                              legend = c("bottom")),
                    align = "hv",
                    ncol = 1, 
                    nrow = 2),
          ggarrange(p_bar_symb_st,
                    p_bar_symb_sp,
                    p_bar_fungi_sp,
                    p_bar_bact_sp,
                    labels = c("E","F","G","H"),
                    widths = c(1,1,1,1),
                    align = "hv" ,
                    ncol = 2, 
                    nrow = 2),
          widths = c(0.7,1),
          ncol = 2, 
          nrow = 1) -> Fig_RF_all

Fig_RF_all


# ********************-------------------------------------------------------------
# Importing cleaned datasets for subsequent analysis ------------------------------

#FUNGI ----------------------------------------------------------------------------
# FILTERING TAXONOMY --------------------------------------------------------------
# Importing taxonomies at 0.6 confidence ------------------------------------------
taxonomy_fungi06 <-read.delim("taxonomy_assignments_coral06_i/consensus_taxonomy.txt", header=TRUE, row.names=1, sep = "\t")
head(taxonomy_fungi06)

# Importing taxonomies at 0.8 confidence ------------------------------------------
#taxonomy_fungi08 <-read.delim("taxonomy_assignments_coral08_i/consensus_taxonomy.txt", header=TRUE, row.names=1)
taxonomy_fungi08 <- read.csv("taxonomy_assignments_coral08_i/consensus_taxonomy.csv", sep = ",", row.names = 1) 
head(taxonomy_fungi08)

identical(rownames(taxonomy_fungi06), rownames(taxonomy_fungi08))
# sample_order <- match(rownames(taxonomy_fungi08), rownames(taxonomy_fungi06))
# taxonomy_fungi06 <- taxonomy_fungi06[,sample_order]

taxonomy_fungi08$Kingdom_06 <- taxonomy_fungi06$Kingdom
head(taxonomy_fungi08)

# removing non fungal OTUs 
subset(taxonomy_fungi08, taxonomy_fungi08$Kingdom_06=="Fungi") -> taxonomy_fungi08_filt
str(taxonomy_fungi08_filt)

# removing lab contamnants
taxonomy_fungi08_filt$Isolate_percent_id[is.na(taxonomy_fungi08_filt$Isolate_percent_id)] <- 0
subset(taxonomy_fungi08_filt, taxonomy_fungi08_filt$Isolate_percent_id>=99) -> contam_otus
contam_otus
dim(contam_otus)

subset(taxonomy_fungi08_filt, !(taxonomy_fungi08_filt$Isolate_percent_id>=99)) -> taxonomy_fungi08_filt
head(taxonomy_fungi08_filt)
dim(taxonomy_fungi08_filt)

dim(taxonomy_fungi08_filt[taxonomy_fungi08_filt$Genus=="Mortierella",])
# 47 Mortierellas left

# remove identified contaminants across all samples
otu_to_remove <- c( "OTU_11","OTU_16","OTU_24","OTU_35","OTU_48","OTU_56","OTU_60",
                    "OTU_63","OTU_67","OTU_70","OTU_86","OTU_91","OTU_104","OTU_110",
                    "OTU_111","OTU_124","OTU_126","OTU_171","OTU_173","OTU_217",
                    "OTU_226","OTU_239","OTU_240","OTU_259","OTU_274","OTU_293",
                    "OTU_314","OTU_359","OTU_379","OTU_387","OTU_396","OTU_405",
                    "OTU_407","OTU_461","OTU_490","OTU_540","OTU_578","OTU_603",
                    "OTU_650","OTU_655","OTU_685","OTU_814","OTU_834","OTU_1040",
                    "OTU_1046","OTU_1048","OTU_1063","OTU_1332","OTU_1465",
                    "OTU_1733","OTU_2042","OTU_3600","OTU_4059","OTU_4345",
                    "OTU_4891","OTU_5017","OTU_5405","OTU_5415","OTU_6684")

taxonomy_fungi08_filt <- taxonomy_fungi08_filt[!(rownames(taxonomy_fungi08_filt) %in% otu_to_remove), ]
dim(taxonomy_fungi08_filt)

# filtering out non-fungal otus
taxonomy_fungi08_filt$Kingdom_06
dim(taxonomy_fungi08_filt[taxonomy_fungi08_filt$Kingdom_06=="Fungi",])

# importing otus, metadata, and sequences
otus_fungi <- read.csv("Fungi/pre_normalization_otus_fungi.csv", header = TRUE, row.names = 1) 
map_fungi <- read.csv("Fungi/pre_normalization_map_fungi.csv", header = TRUE, row.names = 1)
map_fungi <- map_fungi[,10:11]

seq_fungi <- readDNAStringSet("Fungi/pre_normalization_fungi_sequences.fasta",
                              format="fasta", seek.first.rec=TRUE, use.names=TRUE)
taxa_fungi <- taxonomy_fungi08_filt

head(otus_fungi)
head(map_fungi)
head(seq_fungi)

head(taxa_fungi)
taxa_fungi <- taxa_fungi[,1:7]

physeq_fungi <- phyloseq(otu_table(otus_fungi, taxa_are_rows = TRUE),
                         sample_data(map_fungi),
                         tax_table(as.matrix(taxa_fungi)),
                         seq_fungi) 

physeq_fungi

identical(rownames(otu_table(physeq_fungi_new)), rownames(taxa_fungi))

# modyfing fungal taxonomy --------------------------------------------------------
source("../R_functions/ReformatTaxonomy.R")

head(tax_table(physeq_fungi))

ReformatTaxonomy(physeq_fungi, "ITS") -> physeq_fungi_new
head(tax_table(physeq_fungi_new))
rownames(tax_table(physeq_fungi_new))

taxa_fungi <- as.data.frame(as.matrix(tax_table(physeq_fungi_new)))
otus_fungi <- as.data.frame(otu_table(physeq_fungi_new))
metadata_fungi <- as.data.frame(as.matrix(sample_data(physeq_fungi_new)))

# BACTERIA ------------------------------------------------------------------------
# read bacterial data
readRDS("Reid_prok_new_constax.rds") -> physeq_bact
physeq_bact

head(tax_table(physeq_bact))

ReformatTaxonomy(physeq_bact, "ITS") -> physeq_bact_new
head(tax_table(physeq_bact_new))
rownames(tax_table(physeq_bact_new))

sample_data(physeq_bact)[, 10:11]
sample_data(physeq_bact)[sample_data(physeq_bact)$Description=="FJVB42A",]
table(sample_data(physeq_bact)[, 10:11])

otu_bact <- as.data.frame(otu_table(physeq_bact_new))
taxa_bact <- as.data.frame(as.matrix(tax_table(physeq_bact_new)))
metadata_bact <- as.data.frame(as.matrix(sample_data(physeq_bact_new)[, 10:11]))

# APICOMPLEXA ---------------------------------------------------------------------
otus_api <- read.csv("Apicomplexa/pre_normalization_otu_apicomplexa.csv", header = TRUE, row.names = 1) 
map_api <- read.csv("Apicomplexa/pre_normalization_map_apicomplexa.csv", header = TRUE, row.names = 1) 
map_api <- map_api[,10:11]

taxa_api <- read.csv("Apicomplexa/pre_normalization_tax_apicomplexa.csv", header = TRUE, row.names = 1) 
seq_api <- readDNAStringSet("Apicomplexa/pre_normalization_apicomplexa_sequences.fasta",
                            format="fasta", seek.first.rec=TRUE, use.names=TRUE)

map_api[rownames(map_api)=="FJVB42A", colnames(map_api)=="Species"] <- "Porites_cylindrica"

identical(rownames(otus_api), rownames(taxa_api))
identical(rownames(otus_api), names(seq_api))

head(otus_api)
head(map_api)
head(taxa_api)
head(seq_api)

physeq_api <- phyloseq(otu_table(otus_api, taxa_are_rows = TRUE),
                       sample_data(map_api),
                       tax_table(as.matrix(taxa_api)),
                       seq_api) 

physeq_api

# SYMBIODINIUM --------------------------------------------------------------------
otus_symb_mseq <- read.csv("Symbiodinium/otu_table_normalized_symbiodinium.csv", header = TRUE, row.names = 1) 
otus_symb <- read.csv("Symbiodinium/pre_normalization_otu_symbiodinium.csv", header = TRUE, row.names = 1) 
metadata_symb <- read.csv("Symbiodinium/pre_normalization_metadata_symbiodinium.csv", header = TRUE, row.names = 1)
metadata_symb <- metadata_symb[,10:11]

taxa_symb <- read.csv("Symbiodinium/pre_normalization_tax_symbiodinium.csv", header = TRUE, row.names = 1) 
seq_symb <- readDNAStringSet("Symbiodinium//pre_normalization_symbiodinium_sequences.fasta",
                             format="fasta", seek.first.rec=TRUE, use.names=TRUE)

head(otus_symb)
head(metadata_symb)
head(taxa_symb)
head(seq_symb)

metadata_symb[rownames(metadata_symb)=="FJVB42A", colnames(metadata_symb)=="Species"] <- "Porites_cylindrica"

physeq_symb <- phyloseq(otu_table(otus_symb, taxa_are_rows = TRUE),
                        sample_data(metadata_symb),
                        tax_table(as.matrix(taxa_symb)),
                        seq_symb) 

physeq_symb

# saving datasets in .RDS format ---------------------------------------------------
saveRDS(physeq_fungi, "physeq_fungi_2226.rds")
saveRDS(physeq_bact, "physeq_bact.rds")
saveRDS(physeq_api, "physeq_api.rds")
saveRDS(physeq_symb, "physeq_symb.rds")

# ********************-------------------------------------------------------------
# Feature selection ---------------------------------------------------------------
library(mlbench)
library(caret)
library(labgeo)

# define the control using a random forest selection function
control <- rfeControl(functions=rfFuncs, method="cv", number=100)
# run the RFE algorithm
rfe_symb_clr_Status <-
  rfe(
    otus_symb_clr_Status_sel[, 1:(ncol(otus_symb_clr_Status_sel) - 2)],
    otus_symb_clr_Status_sel$Status,
    sizes = c(1:floor(ncol(
      otus_symb_clr_Status_sel) / 2)),
    rfeControl = control)

# This takse forever!!! try a different apprach

rfe_symb_clr_Status <-
  recursive_feature_elimination(
    df = data.frame(otus_symb_clr_Status_sel$Status, otus_symb_clr_Status_sel[, 1:(ncol(otus_symb_clr_Status_sel) - 2)]),
    sizes = c(1:floor(ncol(otus_symb_clr_Status_sel) / 2)),
    nfolds = 10,
    repeats = 1,
    fun = caret::rfFuncs,
    method_rfe = "cv",
    cpu_cores = 8)

rfe_symb_clr_Status
rfe_symb_clr_Status$optVariables

otus_symb_clr_Status_sel[,rfe_symb_clr_Status$optVariables] -> otus_symb_clr_Status_sel

# Identify predictors of Symbiodium presence in corals ------------------------------------------

# Note Random forests are built on decision trees, and decision trees are sensitive to class imbalance.
# Each tree is built on a bag, and each bag is a uniform random sample from the data (with replacement).
# Therefore each tree will be biased in the same direction and magnitude (on average) by class imbalance.

read.csv("map_symbio_clades.csv", header = TRUE, row.names = 1) -> metadata_Symb_clades
metadata_Symb_clades

table(metadata_Symb_clades$Symbio)
dim(metadata_Symb_clades)
dim(metadata_bact)

# match to bacteria samples
setdiff(rownames(metadata_Symb_clades), rownames(metadata_bact)) # present in first not in secont object
setdiff(rownames(metadata_bact), rownames(metadata_Symb_clades)) 
to_remove <- c("FJVB21", "FJVB27","FJVB30B",
               "FJVB36B","FJVB42C","FJVB3B",
               "FJVB6B","FJVB45", "FJVB47B", 
               "FJVB16", "FJVB23")

metadata_Symb_clades_new <-
  metadata_Symb_clades[!rownames(metadata_Symb_clades)%in%(to_remove), ] 
dim(metadata_Symb_clades_new)

metadata_bact_new <-
  metadata_bact[!rownames(metadata_bact)%in%(to_remove), ] 
dim(metadata_bact_new)

otus_bact_new <-
  otus_bact[, !colnames(otus_bact)%in%(to_remove)] 
otus_bact_new <-
  otus_bact_new[rowSums(otus_bact_new)>0, ]
dim(otus_bact_new)

taxa_bact_new <-
  taxa_bact[rownames(otus_bact_new), ]
dim(taxa_bact_new)

table(metadata_Symb_clades_new$Symbio)
colnames(metadata_Symb_clades_new)[3] <- "Clade"
metadata_Symb_clades_new

table(metadata_Symb_clades_new$Clade)

# centred log-ratio (CLR) transformation 
otus_bact_new_clr <- clr(otus_bact_new)
range(otus_bact_new_clr)

# BACTERIA ---------------------------------------------------------------------------------------------
# > CASSIFYING SYMBIO CLADES with BACTERIAL OTUs --------------------------------------------------
# Adding metadata
otus_bact_new_clr_Clade <- data.frame(t(otus_bact_new_clr))

identical(rownames(metadata_Symb_clades_new), rownames(otus_bact_new_clr_Clade))
fungi_clade <- match(rownames(otus_bact_new_clr_Clade), rownames(metadata_Symb_clades_new))
metadata_Symb_clades_new <- metadata_Symb_clades_new[fungi_clade,]

otus_bact_new_clr_Clade$Clade <- metadata_Symb_clades_new$Clade
head(otus_bact_new_clr_Clade)

hist(
  colSums(otus_bact_new_clr_Clade[, -2333]),
  breaks = 100,
  col = "grey",
  xlab = "Read number",
  ylab = "Frequency",
  main = "Central Log Ratio"
)
dev.off()

# Recursive feature elimination 
rfe_bact_new_clr_Clade <-
  recursive_feature_elimination(
    df = data.frame(otus_bact_new_clr_Clade$Clade, otus_bact_new_clr_Clade[, -2333]),
    sizes = seq(1, floor(nrow(otus_bact_new_clr)), by = 5),
    nfolds = 10,
    repeats = 5,
    fun = caret::rfFuncs,
    method_rfe = "cv",
    cpu_cores = 8)

rfe_bact_new_clr_Clade
rfe_bact_new_clr_Clade$optVariables
rfe_bact_new_clr_Clade$resample


# plot 
ggplot(rfe_bact_new_clr_Clade) + theme_bw() +
geom_vline(xintercept = 6, linetype="dotted", color = "blue", size=1)

otus_bact_new_clr_Clade[,rfe_bact_new_clr_Clade$optVariables] -> otus_bact_new_clr_Clade_sel
dim(otus_bact_new_clr_Clade_sel)
head(otus_bact_new_clr_Clade_sel)

identical(rownames(otus_bact_new_clr_Clade_sel), rownames(otus_bact_new_clr_Clade))
otus_bact_new_clr_Clade_sel$Clade <- otus_bact_new_clr_Clade$Clade

# try tuning the model first ------------------------------------------------------
floor(sqrt(ncol(otus_bact_new_clr_Clade_sel) - 1))

set.seed(51)

RF_bact_Cl <-
  randomForest(
    x = otus_bact_new_clr_Clade_sel[, -17],
    y = otus_bact_new_clr_Clade_sel$Clade,
    na.action = na.roughfix,
    ntree = 1001,
    mtry = 4,
    nodesize = 1,
    importance = TRUE,
    proximity = TRUE,
    do.trace=200
  )

RF_bact_Cl
plot(RF_bact_Cl)

# removing unim-ortant variables using MeanDecreaseGini ---------------------------
head(importance(RF_bact_Cl))

##PlotMean Decrease Gini and Mean Decrease Accuracy
varImpPlot(RF_bact_Cl)

# ASSESSING MODEL FIT -------------------------------------------------------------
# Now check to see if the random forest is actually big enough... Up to a point, the more
# trees in the forest, the better. You can tell when you've made enough when the OOB
# no longer improves.

oob_error_bact_Cl <- data.frame(
  Trees = rep(1:nrow(RF_bact_Cl$err.rate), times = 5),
  Clade = rep(
    c("OOB", "Cladocopium_C1", "Cladocopium_C15","Cladocopium_C3", "Durusdinium_D1"),
    each = nrow(RF_bact_Cl$err.rate)
  ),
  Error = c(
    RF_bact_Cl$err.rate[, "OOB"],
    RF_bact_Cl$err.rate[, "Cladocopium_C1"],
    RF_bact_Cl$err.rate[, "Cladocopium_C15"],
    RF_bact_Cl$err.rate[, "Cladocopium_C3"],
    RF_bact_Cl$err.rate[, "Durusdinium_D1"]
  )
)

#Cladocopium_C1 Cladocopium_C15  Cladocopium_C3  Durusdinium_D1

oob_error_bact_Cl$Clade <- factor(oob_error_bact_Cl$Clade,
                            levels =  c("Cladocopium_C1", "Cladocopium_C3","Cladocopium_C15", "Durusdinium_D1", "OOB"))

p_oob_error_bact_Cl <- ggplot(data=oob_error_bact_Cl, aes(x=Trees, y=Error)) +
  ggtitle("Bacteria - Errors") +
  theme_classic() +
  annotate("text", x=500, y=1, label="OOB estimate: 27.78%", size=2.5) +
  geom_line(aes(color=Clade)) + 
  scale_color_manual(values = palette_TX)  +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 7.5)) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=1)) +
  theme(legend.position="right")

p_oob_error_bact_Cl


# ASSESS MODEL PERFORMANCE WITH PERMUTATIONS --------------------------------------
set.seed(63)

perm_RF_bact_Cl <-
  rf.significance(
    x = RF_bact_Cl,
    xdata = otus_bact_new_clr_Clade[, -2333],
    nperm = 999,
    nmtry = 4,
    ntree = 1001
  )

perm_RF_bact_Cl # Model significant at p=0.001


#PLOTTING NMDS ------------------------------------------------------------
dist_bact_Cl <- as.dist(1 - RF_bact_Cl$proximity)
mds_bact_Cl <- cmdscale(dist_bact_Cl, eig = TRUE, x.ret = TRUE)
var_bact_Cl <- round(mds_bact_Cl$eig / sum(mds_bact_Cl$eig) * 100, 1)

head(otus_bact_new_clr_Clade)

values_bact_Cl <- mds_bact_Cl$points
values_bact_Cl

identical(rownames(otus_bact_new_clr_Clade), rownames(values_bact_Cl))

mds_bact_Cl <- data.frame(
  Sample = rownames(values_bact_Cl),
  X = values_bact_Cl[, 1],
  Y = values_bact_Cl[, 2],
  Clade = otus_bact_new_clr_Clade$Clade
)

head(mds_bact_Cl)

predicted_bact_Cl <- as.data.frame(RF_bact_Cl$predicted)
colnames(predicted_bact_Cl) <- "Predicted"
head(predicted_bact_Cl)
identical(rownames(mds_bact_Cl),rownames(predicted_bact_Cl))
mds_bact_Cl$Predicted <- predicted_bact_Cl$Predicted

mds_bact_Cl$Clade <-
  factor(mds_bact_Cl$Clade, levels = c("Cladocopium_C1", "Cladocopium_C3",
                                     "Cladocopium_C15", "Durusdinium_D1"))
mds_bact_Cl$Predicted <-
  factor(mds_bact_Cl$Predicted, levels = c("Cladocopium_C1", "Cladocopium_C3",
                                       "Cladocopium_C15", "Durusdinium_D1"))

p_nmds_bact_Cl = ggplot(data=mds_bact_Cl, aes(x=X, y=Y, label=Sample, color=Clade, shape=Predicted)) + 
  geom_point(size=1.2) +
  #geom_text_repel(label=rownames(mds_data), size = 2, min.segment.length = 0.5) +
  theme_classic() +
  xlab(paste("MDS1 - ", var_bact_Cl[1], "%", sep="")) +
  ylab(paste("MDS2 - ", var_bact_Cl[2], "%", sep="")) +
  ggtitle("Bacteria - Clade") +
  scale_color_manual(values = palette_TX) +
  scale_shape_manual(values = c(3,4,8,13)) +
  #annotate("text", -Inf, -Inf,  hjust = - 0.2, vjust = - 1.2, 
  #         label="OOB = 27.78%", size=2.5) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) + # to adjust decimals
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 7.5)) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=1, order = 1), 
         shape=guide_legend(ncol =1, order = 2)) +
  theme(legend.position="right")

p_nmds_bact_Cl


# Identifying Important Features --------------------------------------------------
imp_bact_Cl <- as.data.frame(RF_bact_Cl$importance)
imp_bact_Cl$features <- rownames(imp_bact_Cl)
imp_bact_Cl <- dplyr::arrange(imp_bact_Cl, desc(MeanDecreaseAccuracy))
head(imp_bact_Cl)
dim(imp_bact_Cl)

taxa_bact[rownames(imp_bact_Cl),] -> taxa_bact_filt_Cl
dim(taxa_bact_filt_Cl)
head(taxa_bact_filt_Cl)

identical(rownames(taxa_bact_filt_Cl), rownames(imp_bact_Cl))
order_bact_Cl <-
  match(rownames(taxa_bact_filt_Cl), rownames(imp_bact_Cl))
imp_bact_Cl <- imp_bact_Cl[order_bact_Cl, ]

imp_bact_Cl$Taxonomy <- taxa_bact_filt_Cl$Taxonomy
imp_bact_Cl <-
  imp_bact_Cl[order(imp_bact_Cl$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(imp_bact_Cl)



# Correcting Taxonomy 
imp_bact_Cl$Taxonomy <- as.character(imp_bact_Cl$Taxonomy)

imp_bact_Cl[imp_bact_Cl == "Zotu10491-Alphaproteobacteria Sar11 clade"] <- "Zotu10491-Pelagibacteraceae"
imp_bact_Cl[imp_bact_Cl$Taxonomy=="Zotu10491-Pelagibacteraceae", ]

imp_bact_Cl[imp_bact_Cl == "Zotu644-Alphaproteobacteria"] <- "Zotu644-Rhodospirillaceae"
imp_bact_Cl[imp_bact_Cl$Taxonomy=="Zotu644-Rhodospirillaceae", ]


p_bar_bact_Cl = ggplot(data=imp_bact_Cl[1:15,]) + 
  geom_bar(aes(x= reorder(Taxonomy, -MeanDecreaseAccuracy),
               y= MeanDecreaseAccuracy), color="grey80", fill="grey80",stat="identity") +
  labs(title = "Top Bacteria OTUs", 
       x= "OTU", y= "Mean Decrease\nin Accuracy") +
  theme_classic() +
  coord_flip() +
  #ylim(0, 0.03) +
  #scale_y_continuous(breaks=c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05), limits = c(0, 0.05)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, size = 7 ,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 5, hjust = 1, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="none") #theme(legend.position=c(0.8, 0.8))

p_bar_bact_Cl


# Fig bacteria test
ggarrange(p_oob_error_bact_Cl,
          p_nmds_bact_Cl,
          p_bar_bact_Cl,
          ncol = 3,
          nrow = 1, 
          align = "h", 
          widths = c(1,1,0.65))



# FUNGI ------------------------------------------------------------------------------------------
# match to bacteria samples
dim(metadata_fungi)
dim(metadata_Symb_clades)

setdiff(rownames(metadata_Symb_clades), rownames(metadata_fungi)) # present in first not in secont object
setdiff(rownames(metadata_fungi), rownames(metadata_Symb_clades)) 
to_remove <- c("FJVB49B", "FJVB16", "FJVB23")

metadata_Symb_clades_new <-
  metadata_Symb_clades[!rownames(metadata_Symb_clades)%in%(to_remove), ] 
dim(metadata_Symb_clades_new)

metadata_fungi_new <-
  metadata_fungi[rownames(metadata_fungi)%in%rownames(metadata_Symb_clades_new), ] 
dim(metadata_fungi_new)

otus_fungi_new <-
  otus_fungi[,colnames(otus_fungi)%in%rownames(metadata_Symb_clades_new)] 
otus_fungi_new <-
  otus_fungi_new[rowSums(otus_fungi_new)>0, ]
dim(otus_fungi_new)

taxa_fungi_new <-
  taxa_fungi[rownames(otus_fungi_new), ]
dim(taxa_fungi_new)

dim(metadata_Symb_clades_new)

# centred log-ratio (CLR) transformation 
otus_fungi_new_clr <- clr(otus_fungi_new)
range(otus_fungi_new_clr)

# > CASSIFYING SYMBIO CLADES with FUGNAL OTUs --------------------------------------------------
# Adding metadata
otus_fungi_new_clr_Clade <- data.frame(t(otus_fungi_new_clr))
dim(otus_fungi_new_clr_Clade)

identical(rownames(metadata_Symb_clades_new), rownames(otus_fungi_new_clr_Clade))
fungi_clade <- match(rownames(otus_fungi_new_clr_Clade), rownames(metadata_Symb_clades_new))
metadata_Symb_clades_new <- metadata_Symb_clades_new[fungi_clade, ]

colnames(metadata_Symb_clades_new)[3] <- "Clade"
metadata_Symb_clades_new

otus_fungi_new_clr_Clade$Clade <- metadata_Symb_clades_new$Clade
head(otus_fungi_new_clr_Clade)
dim(otus_fungi_new_clr_Clade)

hist(
  colSums(otus_fungi_new_clr_Clade[, -2585]),
  breaks = 100,
  col = "grey",
  xlab = "Read number",
  ylab = "Frequency",
  main = "Central Log Ratio"
)
dev.off()

# Recursive feature elimination 
rfe_fungi_new_clr_Clade <-
  recursive_feature_elimination(
    df = data.frame(otus_fungi_new_clr_Clade$Clade, otus_fungi_new_clr_Clade[, -2585]),
    sizes = seq(1, floor(nrow(otus_fungi_new_clr)), by = 5),
    nfolds = 10,
    repeats = 5,
    fun = caret::rfFuncs,
    method_rfe = "cv",
    cpu_cores = 8)

rfe_fungi_new_clr_Clade
rfe_fungi_new_clr_Clade$optVariables
rfe_fungi_new_clr_Clade$resample


# plot 
ggplot(rfe_fungi_new_clr_Clade) + theme_bw() +
geom_vline(xintercept = 26, linetype="dotted", color = "blue", size=1)

otus_fungi_new_clr_Clade[,rfe_fungi_new_clr_Clade$optVariables] -> otus_fungi_new_clr_Clade_sel
dim(otus_fungi_new_clr_Clade_sel)
head(otus_fungi_new_clr_Clade_sel)

identical(rownames(otus_fungi_new_clr_Clade_sel), rownames(otus_fungi_new_clr_Clade))
otus_fungi_new_clr_Clade_sel$Clade <- otus_fungi_new_clr_Clade$Clade

# try tuning the model first ------------------------------------------------------
floor(sqrt(ncol(otus_fungi_new_clr_Clade_sel) - 1))
dim(otus_fungi_new_clr_Clade_sel)

set.seed(51)

RF_fungi_Cl <-
  randomForest(
    x = otus_fungi_new_clr_Clade_sel[, -27],
    y = otus_fungi_new_clr_Clade_sel$Clade,
    na.action = na.roughfix,
    ntree = 1001,
    mtry = 5,
    nodesize = 1,
    importance = TRUE,
    proximity = TRUE,
    do.trace=200
  )

RF_fungi_Cl
plot(RF_fungi_Cl)


# ASSESS MODEL PERFORMANCE WITH PERMUTATIONS --------------------------------------
set.seed(63)

perm_RF_fungi_Cl <-
  rf.significance(
    x = RF_fungi_Cl,
    xdata = otus_fungi_new_clr_Clade[, -2585],
    nperm = 999,
    nmtry = 5,
    ntree = 1001
  )

perm_RF_fungi_Cl # Model significant at p=0.001


# ASSESSING MODEL FIT -------------------------------------------------------------
# Now check to see if the random forest is actually big enough... Up to a point, the more
# trees in the forest, the better. You can tell when you've made enough when the OOB
# no longer improves.

oob_error_fungi_Cl <- data.frame(
  Trees = rep(1:nrow(RF_fungi_Cl$err.rate), times = 5),
  Clade = rep(
    c("OOB", "Cladocopium_C1", "Cladocopium_C15","Cladocopium_C3", "Durusdinium_D1"),
    each = nrow(RF_fungi_Cl$err.rate)
  ),
  Error = c(
    RF_fungi_Cl$err.rate[, "OOB"],
    RF_fungi_Cl$err.rate[, "Cladocopium_C1"],
    RF_fungi_Cl$err.rate[, "Cladocopium_C15"],
    RF_fungi_Cl$err.rate[, "Cladocopium_C3"],
    RF_fungi_Cl$err.rate[, "Durusdinium_D1"]
  )
)

#Cladocopium_C1 Cladocopium_C15  Cladocopium_C3  Durusdinium_D1

oob_error_fungi_Cl$Clade <- factor(oob_error_fungi_Cl$Clade,
                                   levels =  c("Cladocopium_C1", "Cladocopium_C3","Cladocopium_C15", "Durusdinium_D1", "OOB"))

p_oob_error_fungi_Cl <- ggplot(data=oob_error_fungi_Cl, aes(x=Trees, y=Error)) +
  ggtitle("Fungi - Errors") +
  theme_classic() +
  annotate("text", x=500, y=0.25, label="OOB estimate: 45.9%", size=2.5) +
  geom_line(aes(color=Clade)) + 
  scale_color_manual(values = palette_TX)  +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 7.5)) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=1)) +
  theme(legend.position="right")

p_oob_error_fungi_Cl


#PLOTTING NMDS ------------------------------------------------------------
dist_fungi_Cl <- as.dist(1 - RF_fungi_Cl$proximity)
mds_fungi_Cl <- cmdscale(dist_fungi_Cl, eig = TRUE, x.ret = TRUE)
var_fungi_Cl <- round(mds_fungi_Cl$eig / sum(mds_fungi_Cl$eig) * 100, 1)

values_fungi_Cl <- mds_fungi_Cl$points
mds_fungi_Cl <- data.frame(
  Sample = rownames(values_fungi_Cl),
  X = values_fungi_Cl[, 1],
  Y = values_fungi_Cl[, 2],
  Clade = otus_fungi_new_clr_Clade$Clade
)

head(mds_fungi_Cl)

predicted_fungi_Cl <- as.data.frame(RF_fungi_Cl$predicted)
colnames(predicted_fungi_Cl) <- "Predicted"
head(predicted_fungi_Cl)
identical(rownames(mds_fungi_Cl),rownames(predicted_fungi_Cl))
mds_fungi_Cl$Predicted <- predicted_fungi_Cl$Predicted

mds_fungi_Cl$Clade <-
  factor(mds_fungi_Cl$Clade, levels = c("Cladocopium_C1", "Cladocopium_C3",
                                        "Cladocopium_C15", "Durusdinium_D1"))
mds_fungi_Cl$Predicted <-
  factor(mds_fungi_Cl$Predicted, levels = c("Cladocopium_C1", "Cladocopium_C3",
                                            "Cladocopium_C15", "Durusdinium_D1"))

p_nmds_fungi_Cl = ggplot(data=mds_fungi_Cl, aes(x=X, y=Y, label=Sample, color=Clade, shape=Predicted)) + 
  geom_point(size=1.2) +
  #geom_text_repel(label=rownames(mds_data), size = 2, min.segment.length = 0.5) +
  theme_classic() +
  xlab(paste("MDS1 - ", var_fungi_Cl[1], "%", sep="")) +
  ylab(paste("MDS2 - ", var_fungi_Cl[2], "%", sep="")) +
  ggtitle("Fungi - Clade") +
  scale_color_manual(values = palette_TX) +
  scale_shape_manual(values = c(3,4,8,13)) +
  #annotate("text", -Inf, -Inf,  hjust = - 0.2, vjust = - 1.2, 
  #         label="OOB = 27.78%", size=2.5) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) + # to adjust decimals
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 7.5)) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=1, order = 1), 
         shape=guide_legend(ncol =1, order = 2)) +
  theme(legend.position="right")

p_nmds_fungi_Cl


# Identifying Important Features --------------------------------------------------
imp_fungi_Cl <- as.data.frame(RF_fungi_Cl$importance)
imp_fungi_Cl$features <- rownames(imp_fungi_Cl)
imp_fungi_Cl <- dplyr::arrange(imp_fungi_Cl, desc(MeanDecreaseAccuracy))
head(imp_fungi_Cl)
dim(imp_fungi_Cl)

taxa_fungi[rownames(imp_fungi_Cl),] -> taxa_fungi_filt_Cl
dim(taxa_fungi_filt_Cl)
head(taxa_fungi_filt_Cl)

identical(rownames(taxa_fungi_filt_Cl), rownames(imp_fungi_Cl))
order_fungi_Cl <-
  match(rownames(taxa_fungi_filt_Cl), rownames(imp_fungi_Cl))
imp_fungi_Cl <- imp_fungi_Cl[order_fungi_Cl, ]

imp_fungi_Cl$Taxonomy <- taxa_fungi_filt_Cl$Taxonomy
imp_fungi_Cl <-
  imp_fungi_Cl[order(imp_fungi_Cl$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(imp_fungi_Cl)


# Correcting Taxonomy 
imp_fungi_Cl$Taxonomy <- as.character(imp_fungi_Cl$Taxonomy)

imp_fungi_Cl[imp_fungi_Cl == "Zotu15-Fungi"] <- "Zotu15-Ganoderma sp."
imp_fungi_Cl[imp_fungi_Cl$Taxonomy=="Zotu15-Ganoderma sp.", ]

imp_fungi_Cl[imp_fungi_Cl == "Zotu2-Eurotiales"] <- "Zotu2-Aspergillus sp."
imp_fungi_Cl[imp_fungi_Cl$Taxonomy=="Zotu2-Aspergillus sp.", ]


imp_fungi_Cl[imp_fungi_Cl == "Zotu2932-Eurotiales"] <- "Zotu2932-Aspergillus sp."
imp_fungi_Cl[imp_fungi_Cl$Taxonomy=="Zotu2932-Aspergillus sp.", ]

imp_fungi_Cl[imp_fungi_Cl == "Zotu2080-Eurotiales"] <- "Zotu2080-Aspergillus restrictus"
imp_fungi_Cl[imp_fungi_Cl$Taxonomy=="Zotu2080-Aspergillus restrictus", ]

imp_fungi_Cl[imp_fungi_Cl == "Zotu1293-Dothideales"] <- "Zotu1293-Hortaea sp."
imp_fungi_Cl[imp_fungi_Cl$Taxonomy=="Zotu1293-Hortaea sp.", ]



p_bar_fungi_Cl = ggplot(data=imp_fungi_Cl[1:15,]) + 
  geom_bar(aes(x= reorder(Taxonomy, -MeanDecreaseAccuracy),
               y= MeanDecreaseAccuracy), color="grey80", fill="grey80",stat="identity") +
  labs(title = "Top Fungi OTUs", 
       x= "OTU", y= "Mean Decrease\nin Accuracy") +
  theme_classic() +
  coord_flip() +
  #ylim(0, 0.03) +
  #scale_y_continuous(breaks=c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05), limits = c(0, 0.05)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, size = 7 ,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 5, hjust = 1, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="none") #theme(legend.position=c(0.8, 0.8))

p_bar_fungi_Cl


# Fig fungieria test
ggarrange(p_oob_error_fungi_Cl,
          p_nmds_fungi_Cl,
          p_bar_fungi_Cl,
          ncol = 3,
          nrow = 1, 
          align = "h", 
          widths = c(1,1,0.65))



# **** FIGURE 7 - SYMBIO-OTU CLADE MODELS -----------------------------------------------------------
title_a = text_grob("Bacteria",size=10, face=2)
title_b = text_grob("Fungi",size=10, face=2)

Fig_all_otu_clade <-
ggarrange(
  ggarrange(p_oob_error_bact_Cl,
            p_oob_error_fungi_Cl,
            ncol = 2,
            nrow = 1,
            labels = c("A","B"),
            align = "h",
            common.legend = TRUE,
            legend = "right"),
  ggarrange(p_nmds_bact_Cl,
          p_nmds_fungi_Cl,
          ncol = 2,
          nrow = 1,
          labels = c("C","D"),
          align = "h",
          common.legend = TRUE,
          legend = "right"),
  ggarrange(p_bar_bact_Cl,
          p_bar_fungi_Cl,
          ncol = 2,
          nrow = 1,
          labels = c("E","F"),
          align = "h",
          common.legend = TRUE,
          legend = "right"),
  ncol = 1,
  nrow = 3,
  heights = c(1,1,1.1))

Fig_all_otu_clade

