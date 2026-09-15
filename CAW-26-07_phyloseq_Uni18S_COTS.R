rm(list= ls())

getwd()

setwd("C:/Users/joelle.lousberg/OneDrive - Cawthron/Documents/Chapter_2/R_step3/Uni18S")

####read in blast taxonomy
####read in blast results
library(ggplot2)
library(stringr) # not strictly required but handy
library(readr)
library(seqinr)
library(data.table)
library(plyr)
library(dplyr)
library(knitr)
library(data.table)
library(reshape2)
library(grid)
library(ape)
library(gtable)
library(car)
library(tibble)
library(cowplot)
library(vegan)
library(phyloseq)
library(ggpubr)
library(CawthronColours)
library(tidyr)
library(scales)

ps_run <- readRDS("ps_run.rds")

###get sample data
sample.data <-fread("C:/Users/joelle.lousberg/OneDrive - Cawthron/Documents/Chapter_2/R_step3/Uni18S/CAW-26-07_sample_sheet.csv")
names(sample.data)

sample.data$Lake <- as.factor(sample.data$Lake)
levels(sample.data$Lake)

sample.data$Sample_type<- as.factor(sample.data$Sample_type)
levels(sample.data$Sample_type)

sample.data <- as.data.frame(sample.data)
rownames(sample.data) <- sample.data$Sample_ID 

sample_data(ps_run) <- sample.data
ps_run
#rename it so it works for the code
phyloseq1 <- ps_run


#_______________________________________

##set a plot theme
my_theme = theme_bw(base_size = 12) + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.text.x = element_text(face = "bold.italic", size = 10),
  axis.title.y = element_text(vjust= 2.2),
  axis.text = element_text(colour = "black"),
  axis.text.x = element_text(size = 10, angle = 45, hjust =1),
  legend.position="top",
  legend.text = element_text(size = 12),  
  legend.box = "horizontal",
  legend.key.size = unit(1,"line")) + 
  (theme(plot.margin = unit(c(.65,.65,.65,.65), "cm")))

#Choose your colour palette
#colour50 <- c(
#  "#E60026", "#FF6F00", "#FFD300", "#A4DD00", "#00B200",
#  "#00C896", "#00BFFF", "#0068FF", "#4900FF", "#7F00FF",
#  "#B200FF", "#FF00E6", "#FF007F", "#FF0033", "#FF4C4C",
#  "#FF944C", "#FFCC4C", "#D4FF4C", "#94FF4C", "#4CFF4C",
#  "#4CFF94", "#4CFFD4", "#4CFFFF", "#4CC9FF", "#4C7FFF",
#  "#FF4C94", "#FF7373", "#FFA573", "#FFD173", "#E6FF73",
#  "#B3FF73", "#73FF73", "#73FFB3", "#73FFE6", "#73FFFF",
#  "#73D9FF", "#7399FF", "#7373FF", "#B373FF", "#E673FF",
#  "#FF73E6", "#FF73B3", "#FF7399", "#FF9999", "#FFCC99"
#)

colourCAW <- c(
  "#d45f74", "#fbb15d", "#7f1424", "#00549e","#6860a0", "#4ba791", "#ffcc28","#9b99cd", "#c3ce58", 
  "#702365", "#4CFFFF", "#acd58e", "#7ac4d3", "#f8d3ca", "#2f4926", "#3398d2", "#4a7637", "#175c7d", 
  "#c3d3c2", "#67bd45", "#c45b28", "#e8b5d4", "#2f725e", "#85243f", "#76c044", "#75a54a", "#bc3635", 
  "#205128", "#35356d", "#083631", "#fbd872", "#c22c43")


###Look at data
Phylum1 <- prune_samples(sample_sums(phyloseq1) > 0,phyloseq1)
#Phylum <- transform_sample_counts(Phylum1, function(x) x/sum(x))
glom <- tax_glom(Phylum1, taxrank = 'phylum', NArm = FALSE)
Phylumdat <- psmelt(glom)
Phylumdat$phylum <- as.character(Phylumdat$phylum)
medians <- ddply(Phylumdat, ~phylum, function(x) c(median=median(x$Abundance)))
remainder <- medians[medians$median <= 0.001,]$phylum
Phylumdat[Phylumdat$phylum %in% remainder,]$phylum <- 'zz_Other'
Phylum1

Phylum2 <- ggplot(Phylumdat, aes(x = Sample_Name, y = Abundance, fill = phylum))  +
  geom_bar(stat = "identity", position = "stack") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.~Lake, scales = "free", space = "free") + 
  xlab("Sample") +
  scale_fill_manual(values= get_pal("caw_cat_1"))+
  scale_y_continuous(expand = c(0,0)) + my_theme
Phylum2


phyloseq_merged <- phyloseq1


saveRDS(phyloseq_merged, file = "phyloseq_Uni18S_26_07_clean.rds")


##_________________________

