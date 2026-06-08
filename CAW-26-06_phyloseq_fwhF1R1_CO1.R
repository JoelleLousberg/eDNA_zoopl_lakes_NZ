rm(list= ls())

getwd()

setwd("C:/Users/joelle.lousberg/OneDrive - Cawthron/Documents/Chapter_2/R_step3/fwhF1R1_CO1")


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

#Get nonchim values out of rds file
# Replace with your file path
seqtab.nochim <- readRDS("fwhF1R1.seqtab.nochim.rds")
nonchim_counts <- rowSums(seqtab.nochim)
nonchim_table <- data.frame(
  Sample = rownames(seqtab.nochim),
  NonChimera_Reads = nonchim_counts
)
write.csv(nonchim_table,
          file = "nonchimera_reads_per_sample.csv",
          row.names = FALSE)


###Use MEGAN COMMUNITY EDITION - to do LCA (lowest common ancestor analysis)
##We need to sort the phytoplankton since it is misassigned in BLAST
blast.tax <- fread("fwhF1R1_CO1.nochim_ASVs_blast_BOLD-ex.txt", header = FALSE)
head(blast.tax)
blast.tax2 = colsplit(blast.tax$V2, ';', c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species")) #since taxonomy table only consist of 2 columns #(feature ID and taxa) the 7 hierarchies of taxa need to be put in different columns.  string: Taxon (header of column that needs to be split up); pattern ;
row.names(blast.tax2) = blast.tax$`V1`
blast.tax2 <- blast.tax2 %>% mutate_all(na_if,"")
blast.tax2$V1 <- rownames(blast.tax2)

move_taxonomy <- function(df) {
  df %>%
    mutate(
      domain= across(everything(), ~ ifelse(grepl("^d__", .), ., NA_character_)) %>%
        tidyr::unite(division, everything(), na.rm = TRUE, sep = ""),
      kingdom = across(everything(), ~ ifelse(grepl("^k__", .), ., NA_character_)) %>%
        tidyr::unite(kingdom, everything(), na.rm = TRUE, sep = ""),
      phylum = across(everything(), ~ ifelse(grepl("^p__", .), ., NA_character_)) %>%
        tidyr::unite(phylum, everything(), na.rm = TRUE, sep = ""),
      class = across(everything(), ~ ifelse(grepl("^c__", .), ., NA_character_)) %>%
        tidyr::unite(class, everything(), na.rm = TRUE, sep = ""),
      order = across(everything(), ~ ifelse(grepl("^o__", .), ., NA_character_)) %>%
        tidyr::unite(order, everything(), na.rm = TRUE, sep = ""),
      family = across(everything(), ~ ifelse(grepl("^f__", .), ., NA_character_)) %>%
        tidyr::unite(family, everything(), na.rm = TRUE, sep = ""),
      genus = across(everything(), ~ ifelse(grepl("^g__", .), ., NA_character_)) %>%
        tidyr::unite(genus, everything(), na.rm = TRUE, sep = ""),
      species = across(everything(), ~ ifelse(grepl("^s__", .), ., NA_character_)) %>%
        tidyr::unite(species, everything(), na.rm = TRUE, sep = "")
    )
}


blast.tax2.sorted <- move_taxonomy(blast.tax2)
blast.tax2.sorted <- blast.tax2.sorted %>% select(V1, domain, kingdom, phylum, class, order, family, genus, species)
blast.tax2.sorted <- as.data.frame(blast.tax2.sorted)
colnames(blast.tax2.sorted) <- c("V1", "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
blast.tax2.sorted <- as.data.table(blast.tax2.sorted)
blast.tax2.sorted[Class == "c__Dinophyceae", Phylum := "p__Myzozoa"]
blast.tax2.sorted[Class == "c__Cryptophyceae", Phylum := "p__Cryptophyta"]
blast.tax2.sorted[Class == "c__Synurophyceae", Phylum := "p__Ochrophyta"]
blast.tax2.sorted[Class == "c__Phaeophyceae", Phylum := "p__Ochrophyta"]
blast.tax2.sorted[Class == "c__Eustigmatophyceae", Phylum := "p__chrophyta"]
blast.tax2.sorted[Class == "c__Chrysophyceae", Phylum := "p__Ochrophyta"]
blast.tax2.sorted <- data.frame(blast.tax2.sorted)
rownames(blast.tax2.sorted) <- blast.tax2.sorted$V1
blast.tax2.sorted <- blast.tax2.sorted %>% select(-V1)

write.csv(blast.tax2.sorted, "MEGANCO1.blast.csv")


CO1.blast.tax.edit<- blast.tax2.sorted
CO1.blast.tax.edit$Row.names <- rownames(CO1.blast.tax.edit)
head(CO1.blast.tax.edit)

##Ignore this when you want to do the part with "primer X found 23 species in 6lakes"
#filled_tax_df <- CO1.blast.tax.edit %>%
#  mutate(
#    Kingdom = ifelse(is.na(Kingdom), Domain, Kingdom),
#    Phylum = ifelse(is.na(Phylum), Kingdom, Phylum),
#    Class = ifelse(is.na(Class), Phylum, Class),
#    Order = ifelse(is.na(Order), Class, Order),
#    Family = ifelse(is.na(Family), Order, Family),
#    Genus = ifelse(is.na(Genus), Family, Genus),
#    Species = ifelse(is.na(Species), Genus, Species)
#  )
#____________________________________________________

ASVCO1 <- fread("fwhF1R1.nochim_ASVs_counts.tsv", header = T)
head(ASVCO1)
ASVCO1 <- as.data.frame(ASVCO1)
rownames(ASVCO1) <- ASVCO1$V1
head(ASVCO1)
ASVCO1 <- subset(ASVCO1, select = -c(`V1`) )
head(ASVCO1)


##________________________________

###get back into phyloseq
ps_otu_table = otu_table(ASVCO1, taxa_are_rows = T)
head(ps_otu_table)

ps_tax = tax_table(as.matrix(CO1.blast.tax.edit))

#ps_tax = tax_table(as.matrix(filled_tax_df))
#head(ps_tax)
#rownames(ps_tax)=rownames(filled_tax_df)
#colnames(ps_tax)=colnames(filled_tax_df)
#head(ps_tax)

####
sample.data <-fread("C:/Users/joelle.lousberg/OneDrive - Cawthron/Documents/Chapter_2/R_step3/fwhF1R1_CO1/CAW-26-06_sample_sheet.csv")
names(sample.data)
sample.data$Sample.Name <- as.factor(sample.data$Sample_Name)
ps_map = sample_data(sample.data)
head(ps_map)
rownames(ps_map)=sample.data$'Sample_ID'
head(ps_map)


##phyloseq object to create plots

phyloseq_fwhF1R1 <- phyloseq(ps_otu_table, ps_map, ps_tax)
phyloseq_fwhF1R1
ntaxa(phyloseq_fwhF1R1)

# Look at the rds file to get numbers after chimera check
#nonchimera_data <- readRDS("fwhF1R1.seqtab.nochim_all.rds")
#sample_sums_df <- data.frame(
#  Sample = rownames(nonchimera_data),
#  NonChimera_Sum = rowSums(nonchimera_data)
#)
#write.csv(sample_sums_df, 
#          file = "nonchimera_sample_sums.csv", 
#          row.names = FALSE)

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


##First look at data
phyloseq1 <- phyloseq_fwhF1R1
###
Phylum1 <- prune_samples(sample_sums(phyloseq1) > 0,phyloseq1)
#Phylum <- transform_sample_counts(Phylum1, function(x) x/sum(x))
glom <- tax_glom(Phylum1, taxrank = 'Phylum', NArm = FALSE)
Phylumdat <- psmelt(glom)
Phylumdat$Phylum <- as.character(Phylumdat$Phylum)
medians <- ddply(Phylumdat, ~Phylum, function(x) c(median=median(x$Abundance)))
remainder <- medians[medians$median <= 0.001,]$Phylum
Phylumdat[Phylumdat$Phylum %in% remainder,]$Phylum <- 'zz_Other'
Phylum1

Phylum2 <- ggplot(Phylumdat, aes(x = Sample_Name, y = Abundance, fill = Phylum))  +
  geom_bar(stat = "identity", position = "stack") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.~Lake, scales = "free", space = "free") + 
  xlab("Sample") +
  scale_fill_manual(values= get_pal("caw_cat_1"))+
  scale_y_continuous(expand = c(0,0)) + my_theme
Phylum2


####Subtraction of contamination from lab blanks for CAW-25-36 tows
phyloseq_25_36_tow <- subset_samples(
  phyloseq1,
  Run %in% c("fwh_25_36_tow", "fwh_25_36_Blank")
)

Controls = subset_samples(phyloseq_25_36_tow, Run %in% c("fwh_25_36_Blank"))
if (sum(sample_sums(Controls)) > 0) {
  Controls <- filter_taxa(Controls, function(x) sum(x) > 0, TRUE)
} else {
  warning("No reads in Controls — skipping filtering step.")
}

sample_sums(Controls)
taxa_sums(Controls)
#You can change the taxon level here
plot_bar(Controls, fill = "Genus")

#Extract controls from samples
Extraction_neg <- subset_samples(phyloseq_25_36_tow, Run %in% c("fwh_25_36_Blank"))
Extraction_neg
Extraction_neg_max <- apply(as.matrix(otu_table(Extraction_neg)), 1, max)
Extraction_neg_max_vec <- as.vector(Extraction_neg_max)
names(Extraction_neg_max_vec) <- taxa_names(Extraction_neg)  # ensure names stay aligned
Extraction <- as(otu_table(phyloseq_25_36_tow), "matrix")

# Ensure taxa are rows
if (!taxa_are_rows(phyloseq_25_36_tow)) {
  Extraction <- t(Extraction)
}
Extractiondf <- as.data.frame(Extraction)
Extractiondf <- Extractiondf[names(Extraction_neg_max_vec), , drop = FALSE]
Extractiondf <- sweep(Extractiondf, 1, Extraction_neg_max_vec, "-")
Extractiondf[Extractiondf < 0] <- 0
rownames(Extractiondf) <- taxa_names(phyloseq_25_36_tow)
colnames(Extractiondf) <- sample_names(phyloseq_25_36_tow)
stopifnot(identical(rownames(Extractiondf), taxa_names(phyloseq_25_36_tow)))
stopifnot(identical(colnames(Extractiondf), sample_names(phyloseq_25_36_tow)))
Extraction_OTU <- otu_table(as.matrix(Extractiondf), taxa_are_rows = TRUE)
Tutorial.subtractextract.ps <- phyloseq(
  Extraction_OTU,
  sample_data(phyloseq_25_36_tow),
  tax_table(phyloseq_25_36_tow)
)

org.ss <- as.data.frame(sample_sums(phyloseq_25_36_tow))
org.ss$Names <- rownames(org.ss)
new.ss <- as.data.frame(sample_sums(Tutorial.subtractextract.ps))
new.ss$Names <- rownames(new.ss)
control.rem.maxsub <- dplyr::left_join(org.ss, new.ss, by="Names")
colnames(control.rem.maxsub) <- c("Original", "Names", "New")
control.rem.maxsub.final <- control.rem.maxsub %>%
  mutate(perc = New/Original*100)
control.rem.maxsub.final

phyloseq1.5 = subset_samples(phyloseq_25_36_tow, Sample_type != "Blank")

##remove anything not at phylum level
phyloseq2 <- subset_taxa(phyloseq1.5, 
                         !Phylum %in% c("d__Eukaryota", "k__Metazoa", "") & 
                           !is.na(Phylum))
phyloseq_fwhF1R1_25_36_clean <- phyloseq2
saveRDS(phyloseq_fwhF1R1_25_36_clean, file = "phyloseq_fwhF1R1_25_36_clean.rds")


####Subtraction of contamination from lab blanks for CAW-26-06 watersamples
phyloseq_26_06_water <- subset_samples(
  phyloseq1,
  Run %in% c("fwh_26_06_water", "fwh_26_06_Blank")
)

Controls = subset_samples(phyloseq_26_06_water, Run %in% c("fwh_26_06_Blank"))
if (sum(sample_sums(Controls)) > 0) {
  Controls <- filter_taxa(Controls, function(x) sum(x) > 0, TRUE)
} else {
  warning("No reads in Controls — skipping filtering step.")
}

sample_sums(Controls)
taxa_sums(Controls)
#You can change the taxon level here
plot_bar(Controls, fill = "Genus")

#Extract controls from samples
Extraction_neg <- subset_samples(phyloseq_26_06_water, Run %in% c("fwh_26_06_Blank"))
Extraction_neg
Extraction_neg_max <- apply(as.matrix(otu_table(Extraction_neg)), 1, max)
Extraction_neg_max_vec <- as.vector(Extraction_neg_max)
names(Extraction_neg_max_vec) <- taxa_names(Extraction_neg)  # ensure names stay aligned
Extraction <- as(otu_table(phyloseq_26_06_water), "matrix")

# Ensure taxa are rows
if (!taxa_are_rows(phyloseq_26_06_water)) {
  Extraction <- t(Extraction)
}
Extractiondf <- as.data.frame(Extraction)
Extractiondf <- Extractiondf[names(Extraction_neg_max_vec), , drop = FALSE]
Extractiondf <- sweep(Extractiondf, 1, Extraction_neg_max_vec, "-")
Extractiondf[Extractiondf < 0] <- 0
rownames(Extractiondf) <- taxa_names(phyloseq_26_06_water)
colnames(Extractiondf) <- sample_names(phyloseq_26_06_water)
stopifnot(identical(rownames(Extractiondf), taxa_names(phyloseq_26_06_water)))
stopifnot(identical(colnames(Extractiondf), sample_names(phyloseq_26_06_water)))
Extraction_OTU <- otu_table(as.matrix(Extractiondf), taxa_are_rows = TRUE)
Tutorial.subtractextract.ps <- phyloseq(
  Extraction_OTU,
  sample_data(phyloseq_26_06_water),
  tax_table(phyloseq_26_06_water)
)

org.ss <- as.data.frame(sample_sums(phyloseq_26_06_water))
org.ss$Names <- rownames(org.ss)
new.ss <- as.data.frame(sample_sums(Tutorial.subtractextract.ps))
new.ss$Names <- rownames(new.ss)
control.rem.maxsub <- dplyr::left_join(org.ss, new.ss, by="Names")
colnames(control.rem.maxsub) <- c("Original", "Names", "New")
control.rem.maxsub.final <- control.rem.maxsub %>%
  mutate(perc = New/Original*100)
control.rem.maxsub.final

phyloseq1.5 = subset_samples(phyloseq_26_06_water, Run != "fwh_26_06_Blank")

##remove anything not at phylum level
phyloseq2 <- subset_taxa(phyloseq1.5, 
                         !Phylum %in% c("d__Eukaryota", "k__Metazoa", "") & 
                           !is.na(Phylum))
phyloseq_fwhF1R1_26_06_water_clean <- phyloseq2
saveRDS(phyloseq_fwhF1R1_26_06_water_clean, file = "phyloseq_fwhF1R1_26_06_water_clean.rds")



####Subtraction of contamination from lab blanks for CAW-26-06 & CAW-25-36 tows & freezer
phyloseq_26_06_tow <- subset_samples(
  phyloseq1,
  Run == "fwh_26_06_tow" |
  Sample_ID %in% c("CAW-26-06-097", "CAW-25-36-009")
)

Controls = subset_samples(phyloseq_26_06_tow, Sample_type %in% c("Blank"))
if (sum(sample_sums(Controls)) > 0) {
  Controls <- filter_taxa(Controls, function(x) sum(x) > 0, TRUE)
} else {
  warning("No reads in Controls — skipping filtering step.")
}

sample_sums(Controls)
taxa_sums(Controls)
#You can change the taxon level here
plot_bar(Controls, fill = "Genus")

#Extract controls from samples
Extraction_neg <- subset_samples(phyloseq_26_06_tow, Sample_type %in% c("Blank"))
Extraction_neg
Extraction_neg_max <- apply(as.matrix(otu_table(Extraction_neg)), 1, max)
Extraction_neg_max_vec <- as.vector(Extraction_neg_max)
names(Extraction_neg_max_vec) <- taxa_names(Extraction_neg)  # ensure names stay aligned
Extraction <- as(otu_table(phyloseq_26_06_tow), "matrix")

# Ensure taxa are rows
if (!taxa_are_rows(phyloseq_26_06_tow)) {
  Extraction <- t(Extraction)
}
Extractiondf <- as.data.frame(Extraction)
Extractiondf <- Extractiondf[names(Extraction_neg_max_vec), , drop = FALSE]
Extractiondf <- sweep(Extractiondf, 1, Extraction_neg_max_vec, "-")
Extractiondf[Extractiondf < 0] <- 0
rownames(Extractiondf) <- taxa_names(phyloseq_26_06_tow)
colnames(Extractiondf) <- sample_names(phyloseq_26_06_tow)
stopifnot(identical(rownames(Extractiondf), taxa_names(phyloseq_26_06_tow)))
stopifnot(identical(colnames(Extractiondf), sample_names(phyloseq_26_06_tow)))
Extraction_OTU <- otu_table(as.matrix(Extractiondf), taxa_are_rows = TRUE)
Tutorial.subtractextract.ps <- phyloseq(
  Extraction_OTU,
  sample_data(phyloseq_26_06_tow),
  tax_table(phyloseq_26_06_tow)
)

org.ss <- as.data.frame(sample_sums(phyloseq_26_06_tow))
org.ss$Names <- rownames(org.ss)
new.ss <- as.data.frame(sample_sums(Tutorial.subtractextract.ps))
new.ss$Names <- rownames(new.ss)
control.rem.maxsub <- dplyr::left_join(org.ss, new.ss, by="Names")
colnames(control.rem.maxsub) <- c("Original", "Names", "New")
control.rem.maxsub.final <- control.rem.maxsub %>%
  mutate(perc = New/Original*100)
control.rem.maxsub.final

phyloseq1.5 = subset_samples(phyloseq_26_06_tow, Run != "fwh_26_06_Blank")

##remove anything not at phylum level
phyloseq2 <- subset_taxa(phyloseq1.5, 
                         !Phylum %in% c("d__Eukaryota", "k__Metazoa", "") & 
                           !is.na(Phylum))
phyloseq_fwhF1R1_26_06_tow_clean <- phyloseq2
saveRDS(phyloseq_fwhF1R1_26_06_tow_clean, file = "phyloseq_fwhF1R1_26_06_tow_clean.rds")

phyloseq_merged <- merge_phyloseq(phyloseq_fwhF1R1_26_06_tow_clean, phyloseq_fwhF1R1_26_06_water_clean, phyloseq_fwhF1R1_25_36_clean)
saveRDS(phyloseq_merged, file = "phyloseq_fwhF1R1_26_06_clean.rds")



