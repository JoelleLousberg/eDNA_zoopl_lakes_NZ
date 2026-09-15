rm(list= ls())

getwd()

setwd("C:/Users/joelle.lousberg/OneDrive - Cawthron/Documents/Chapter_2/R_step2/SSUF04R22_COTS/DADA2_18S_results")

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
library(tidyr)

ps_run <- readRDS("ps_run.rds")

###get sample data
sample.data <-fread("C:/Users/joelle.lousberg/OneDrive - Cawthron/Documents/Chapter_2/R_step2/SSUF04R22_COTS/DADA2_18S_results/CAW-25-36_sample_data.csv")
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
  "#f48521", "#fbd872", "#acd58e", "#76c044", "#4a7637", "#2f4926", "#4ba791",
  "#2f725e", "#004f52", "#175c7d", "#029cbd", "#7ac4d3", "#9b99cd", "#6860a0",
  "#702365", "#e8b5d4", "#85243f", "#c22c43", "#d45f74", "#f8d3ca", "#7f1424",
  "#bc3635", "#ffcc28", "#fbb15d", "#ffe0ae", "#c45b28", "#9ecce5", "#3398d2",
  "#00549e", "#d383b7", "#a184bd", "#35356d", "#c3d3c2", "#75a54a", "#205128",
  "#67bd45", "#c3ce58", "#737f7e", "#083631", "#4CFFFF"
)


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


####Subtraction of contamination from all lab blanks
Controls = subset_samples(phyloseq1, Sample_type %in% c("Blank"))
Controls = filter_taxa(Controls, function(x) sum(x) > 0, TRUE)
sample_sums(Controls)
taxa_sums(Controls)
#You can change the taxon level here
plot_bar(Controls, fill = "genus")

##If error appears, use this code instead:
# Check which taxa have non-zero counts
nonzero_taxa <- taxa_names(Controls)[rowSums(otu_table(Controls)) > 0]

# Only filter/prune if there are any non-zero taxa
if(length(nonzero_taxa) > 0){
  Controls_filtered <- prune_taxa(nonzero_taxa, Controls)
  message("Number of taxa kept: ", length(nonzero_taxa))
} else {
  Controls_filtered <- Controls
  warning("No non-zero taxa found. Keeping original Controls object.")
}
##

###removing ASV contamination - subtraction
Extraction_neg = subset_samples(phyloseq1, Sample_type %in% c("Blank"))
Extraction_neg
Extraction_neg_max <- apply(as.data.frame(as.matrix(t(otu_table(Extraction_neg)))), 1, max)
Extraction_neg_max_vec <- as.vector(Extraction_neg_max)
Extraction_neg_max_vec
Extraction_neg_sums <- colSums(otu_table(Extraction_neg))
Extraction_neg_sums

Extraction = as(otu_table(phyloseq1), "matrix")
Extraction
Extractiondf = as.data.frame(Extraction)
Extractiondf[,1:length(Extractiondf)] <- sweep(Extractiondf[,1:length(Extractiondf)],2,Extraction_neg_max_vec)
Extractiondf <- replace(Extractiondf, Extractiondf < 0, 0)
str(Extractiondf)
Tutorial.subtractextract.ps <- phyloseq(otu_table(Extractiondf, taxa_are_rows=F),
                                        sample_data(phyloseq1), 
                                        tax_table(phyloseq1))


org.ss <- as.data.frame(sample_sums(phyloseq1))
org.ss$Names <- rownames(org.ss)
new.ss <- as.data.frame(sample_sums(Tutorial.subtractextract.ps))
new.ss$Names <- rownames(new.ss)
control.rem.maxsub <- dplyr::left_join(org.ss, new.ss, by="Names")
colnames(control.rem.maxsub) <- c("Original", "Names", "New")
control.rem.maxsub.final <- control.rem.maxsub %>%
  mutate(perc = New/Original*100)
control.rem.maxsub.final

phyloseq1.5 = subset_samples(phyloseq1, Sample_type != "Blank")

##remove anything not at phylum level
phyloseq2 <- subset_taxa(phyloseq1.5, 
                         !phylum %in% c("Eukaryota", "Metazoa") & 
                           !is.na(phylum))
phyloseq_SSUF04R22_clean <- phyloseq2
saveRDS(phyloseq_SSUF04R22_clean, file = "phyloseq_SSUF04R22_clean.rds")


## Load from here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
phyloseq2 <- readRDS("phyloseq_SSUF04R22_clean.rds")
ntaxa(phyloseq2)

##Subset samples for zoopl only 
library(openxlsx)
phyloseq_zoopl <- subset_taxa(phyloseq2, class == "Branchiopoda" | class == "Hexanauplia" | class == "Eurotatoria" )
phyloseq_zoopl_filtered <- subset_samples(phyloseq_zoopl, Lake != "Pooled lakes A")
phyloseq_zoopl_filtered <- subset_samples(phyloseq_zoopl_filtered, Lake != "Pooled lakes B")
phyloseq_zoopl_filtered <- prune_taxa(taxa_sums(phyloseq_zoopl_filtered) >0, phyloseq_zoopl_filtered)

sam_data_df <- as.data.frame(sample_data(phyloseq_zoopl_filtered)) %>%
  tibble::rownames_to_column("Sample")

## Create a raw data table with all ASVs corresponding to zoopl across all 6 lakes incl. rel. abundances
physeq_present <- prune_taxa(
  taxa_sums(phyloseq_zoopl_filtered) > 0,
  phyloseq_zoopl_filtered
)

tax_df <- tax_table(physeq_present) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ASV")
tax_table_merged <- tax_df %>%
  select(ASV, family, genus, species) %>%
  filter(!(is.na(family) & is.na(genus) & is.na(species))) %>%
  distinct()
otu_df <- otu_table(phyloseq_zoopl_filtered)
if (!taxa_are_rows(phyloseq_zoopl_filtered)) {
  otu_df <- t(otu_df)
}
otu_df <- as.data.frame(otu_df)
otu_long <- otu_df %>%
  tibble::rownames_to_column("ASV") %>%
  pivot_longer(
    -ASV,
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  group_by(Sample) %>%
  mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  filter(Abundance > 0)
taxalist_by_sample <- otu_long %>%
  left_join(
    tax_df %>% select(ASV, family, genus, species),
    by = "ASV"
  ) %>%
  left_join(
    sam_data_df,
    by = "Sample"
  ) %>%
  arrange(Sample_Name, desc(Abundance))

write.xlsx(
  taxalist_by_sample,
  "zoopl_by_lake_incl_abundances_per_ASV.xlsx",
  rowNames = FALSE
)


##Create a summary list of what species, genera, family found across 6lakes with their rel. abundances
taxalist_by_sample <- taxalist_by_sample |>
  mutate(
    Taxon = coalesce(species, genus, family),
    RankUsed = case_when(
      !is.na(species) ~ "species",
      is.na(species) & !is.na(genus) ~ "genus",
      is.na(species) & is.na(genus) & !is.na(family) ~ "family",
      TRUE ~ "Unknown"
    )
  )
taxon_summary_total <- taxalist_by_sample |>
  mutate(Abundance = ifelse(is.na(Abundance), 0, Abundance)) |>  # ensure no NAs
  group_by(Taxon, RankUsed) |>
  summarise(
    TotalAbundance = sum(Abundance),
    .groups = "drop"
  ) |>
  # Recalculate relative abundance across the entire dataset
  mutate(TotalRelativeAbundance = TotalAbundance / sum(TotalAbundance)) |>
  arrange(desc(TotalAbundance))
wb <- createWorkbook()
addWorksheet(wb, "zoopl_across6lakes_total")
writeData(wb, "zoopl_across6lakes_total", taxon_summary_total)
saveWorkbook(wb, "zoopl_across6lakes_incl_abundances.xlsx", overwrite = TRUE)


##Create a zoopl list per lake for future heatmaps
taxa_list_by_sample <- taxalist_by_sample |>
  select(Sample, Taxon) |>
  distinct() |>
  arrange(Sample, Taxon) |>
  left_join(sam_data_df, by = "Sample")
wb <- createWorkbook()
addWorksheet(wb, "zoopl_list_by_lake")
writeData(wb, "zoopl_list_by_lake", taxa_list_by_sample)
saveWorkbook(wb, "zoopl_list_by_lake.xlsx", overwrite = TRUE)

##Test for alpha diversity of zoopl ASVs per sample
alpha_div <- estimate_richness(
  phyloseq_zoopl_filtered,
  measures = c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher")
) %>%
  tibble::rownames_to_column("Sample") %>%
  left_join(
    as.data.frame(sample_data(phyloseq_zoopl_filtered)) %>%
      tibble::rownames_to_column("Sample"),
    by = "Sample"
  )
alpha_div


##Test for alpha diversity of zoopl ASVs across 6lakes
test <- merge_samples(phyloseq_zoopl_filtered, "Sample_type")
alpha_div_6lakes <- estimate_richness(test, measures = c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher"))
alpha_div_6lakes


#Number of species, genera and family assignments in rarefied zoopl subset
phyloseq_zoopl_genus <- tax_glom(phyloseq_zoopl_filtered, "genus")
phyloseq_zoopl_genus
#unique(tax_table(phyloseq_zoopl_genus)[, "Genus"])  to see the list directly in R
phyloseq_zoopl_species <- tax_glom(phyloseq_zoopl_filtered, "species")
phyloseq_zoopl_species
phyloseq_zoopl_family <- tax_glom(phyloseq_zoopl_filtered, "family")
sum(sample_sums(phyloseq_zoopl_family))/sum(sample_sums(phyloseq_zoopl_filtered))*100


#_____________________________________________________________________________________
phyloseq2 <- readRDS("phyloseq_SSUF04R22_clean.rds")
ntaxa(phyloseq2)

###Rarefy step
library(vegan)
sample_sums(phyloseq2)
#min(sample_sums(phyloseq3))
sort(sample_sums(phyloseq2))[1:10]
set.seed(100)
phyloseq2_rare = rarefy_even_depth(phyloseq2, sample.size = 11000, 
                                   replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
phyloseq2 <- phyloseq2_rare
ntaxa(phyloseq2)
#_________________________________________________________________________________________

##Subset samples for zoopl only 
library(openxlsx)
phyloseq_zoopl_filtered_rare <- subset_taxa(phyloseq2, class == "Branchiopoda" | class == "Hexanauplia" | class == "Eurotatoria" )
phyloseq_zoopl_filtered_rare <- subset_samples(phyloseq_zoopl_filtered_rare, Lake != "Pooled lakes A")
phyloseq_zoopl_filtered_rare <- subset_samples(phyloseq_zoopl_filtered_rare, Lake != "Pooled lakes B")
phyloseq_zoopl_filtered_rare <- prune_taxa(taxa_sums(phyloseq_zoopl_filtered_rare) >0, phyloseq_zoopl_filtered_rare)

#include the lakes names in excel sheets
sam_data_df <- as.data.frame(sample_data(phyloseq_zoopl_filtered_rare)) %>%
  tibble::rownames_to_column("Sample")

## Create a raw data table with all ASVs corresponding to zoopl across all 6 lakes incl. rel. abundances
physeq_present <- prune_taxa(
  taxa_sums(phyloseq_zoopl_filtered_rare) > 0,
  phyloseq_zoopl_filtered_rare
)

tax_df <- tax_table(physeq_present) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ASV")
tax_table_merged <- tax_df %>%
  select(ASV, family, genus, species) %>%
  filter(!(is.na(family) & is.na(genus) & is.na(species))) %>%
  distinct()
otu_df <- otu_table(phyloseq_zoopl_filtered_rare)
if (!taxa_are_rows(phyloseq_zoopl_filtered_rare)) {
  otu_df <- t(otu_df)
}
otu_df <- as.data.frame(otu_df)
otu_long <- otu_df %>%
  tibble::rownames_to_column("ASV") %>%
  pivot_longer(
    -ASV,
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  group_by(Sample) %>%
  mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  filter(Abundance > 0)
taxalist_by_sample <- otu_long %>%
  left_join(
    tax_df %>% select(ASV, family, genus, species),
    by = "ASV"
  ) %>%
  left_join(
    sam_data_df,
    by = "Sample"
  ) %>%
  arrange(Sample_Name, desc(Abundance))

write.xlsx(taxalist_by_sample, "zoopl_by_lake_incl_abundances_per_ASV_rare.xlsx")



##Create a summary list of what species, genera, family found across 6lakes with their rel. abundances on rarefied data
taxalist_by_sample <- taxalist_by_sample |>
  mutate(
    Taxon = coalesce(species, genus, family),
    RankUsed = case_when(
      !is.na(species) ~ "species",
      is.na(species) & !is.na(genus) ~ "genus",
      is.na(species) & is.na(genus) & !is.na(family) ~ "family",
      TRUE ~ "Unknown"
    )
  )
taxon_summary_total <- taxalist_by_sample |>
  mutate(Abundance = ifelse(is.na(Abundance), 0, Abundance)) |>  
  group_by(Taxon, RankUsed) |>
  summarise(
    TotalAbundance = sum(Abundance),
    .groups = "drop"
  ) |>
  # Recalculate relative abundance across the entire dataset
  mutate(TotalRelativeAbundance = TotalAbundance / sum(TotalAbundance)) |>
  arrange(desc(TotalAbundance))
wb <- createWorkbook()
addWorksheet(wb, "zoopl_across6lakes_total_rare")
writeData(wb, "zoopl_across6lakes_total_rare", taxon_summary_total)
saveWorkbook(wb, "zoopl_across6lakes_incl_abundances_rare.xlsx", overwrite = TRUE)


##Create a zoopl list per lake for future heatmaps
taxa_list_by_sample <- taxalist_by_sample |>
  select(Sample, Taxon) |>
  distinct() |>
  arrange(Sample, Taxon) |>
  left_join(sam_data_df, by = "Sample")
wb <- createWorkbook()
addWorksheet(wb, "zoopl_list_by_lake_rare")
writeData(wb, "zoopl_list_by_lake_rare", taxa_list_by_sample)
saveWorkbook(wb, "zoopl_list_by_lake_rare.xlsx", overwrite = TRUE)

##Test for alpha diversity of zoopl ASVs per sample
alpha_div <- estimate_richness(
  phyloseq_zoopl_filtered_rare,
  measures = c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher")
) %>%
  tibble::rownames_to_column("Sample") %>%
  left_join(
    as.data.frame(sample_data(phyloseq_zoopl_filtered_rare)) %>%
      tibble::rownames_to_column("Sample"),
    by = "Sample"
  )
alpha_div


##Test for alpha diversity of zoopl ASVs across 6lakes on rarefied data
test <- merge_samples(phyloseq_zoopl_filtered_rare, "Sample_type")
alpha_div_6lakes <- estimate_richness(test, measures = c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher"))
alpha_div_6lakes

#Number of species, genera and family assignments in rarefied zoopl subset
phyloseq_zoopl_genus_rare <- tax_glom(phyloseq_zoopl_filtered_rare, "genus")
phyloseq_zoopl_genus_rare
phyloseq_zoopl_species_rare <- tax_glom(phyloseq_zoopl_filtered_rare, "species")
phyloseq_zoopl_species_rare
phyloseq_zoopl_family_rare <- tax_glom(phyloseq_zoopl_filtered_rare, "family")
sum(sample_sums(phyloseq_zoopl_family_rare))/sum(sample_sums(phyloseq_zoopl_filtered_rare))*100






# OLD CODE _______________________________________________________________________________________

### Write a table with all counts as percentages relative abundance but grouped into zooplankton etc (NOT per Species!)
if (!taxa_are_rows(phyloseq2)) {
  phyloseq2 <- t(phyloseq2)  # transpose OTU table
}
phyloseq_relabund <- transform_sample_counts(phyloseq2, function(x) x / sum(x))
otu_relabund <- as.data.frame(otu_table(phyloseq_relabund))
otu_relabund$Taxon <- rownames(otu_relabund)
tax_df <- as.data.frame(tax_table(phyloseq2), stringsAsFactors = FALSE)
tax_df$Taxon <- rownames(tax_df)
relabund_tax <- merge(otu_relabund, tax_df, by = "Taxon", all.x = TRUE)

relabund_tax <- relabund_tax %>%
  mutate(Group = case_when(
    # Zooplankton
    class %in% c("Hexanauplia", "Branchiopoda", "Eurotatoria") ~ "Zooplankton",
    # Phytoplankton / algae
    phylum %in% c("Chlorophyta", "Streptophyta") ~ "Macro- & Greenalgae",
    # Fungi / Microbes
    phylum %in% c("Ascomycota", "Chytridiomycota", "Cryptomycota", "Microsporidia") ~ "Fungi",
    # Protozoa
    phylum %in% c("Ciliophora", "Cercozoa", "Apicomplexa") ~ "Protists",
    # Hydrozoa (Cnidaria)
    phylum == "Cnidaria" & class == "Hydrozoa" ~ "Hydrozoa",# Water mites (Acari, Hydrachnidia)
    phylum == "Arthropoda" & order %in% c("Hydrachnidia", "Trombidiformes") ~ "Water Mites",
    # Other Arthropods (exclude zooplankton and water mites)
    phylum == "Arthropoda" & !class %in% c("Hexanauplia", "Branchiopoda") & !order %in% c("Hydrachnidia") ~ "Other Arthropoda",
    
    # Annelids
    phylum == "Annelida" ~ "Annelids",
    # Molluscs
    phylum == "Mollusca" ~ "Molluscs",
    # Porifera
    phylum == "Porifera" ~ "Sponges",
    
    # Catch-all remaining
    TRUE ~ "Other"
  ))

# Check that grouping worked
table(relabund_tax$Group)

# Create counts_tax from phyloseq2 (rarefied counts)
otu_counts <- as.data.frame(otu_table(phyloseq2))
otu_counts$Taxon <- rownames(otu_counts)
tax_df <- as.data.frame(tax_table(phyloseq2), stringsAsFactors = FALSE)
tax_df$Taxon <- rownames(tax_df)
counts_tax <- merge(otu_counts, tax_df, by = "Taxon", all.x = TRUE)

# Summarize read counts by group
group_counts <- counts_tax %>%
  left_join(select(relabund_tax, Taxon, Group), by = "Taxon") %>%
  group_by(Group) %>%
  summarise_if(is.numeric, sum, na.rm = TRUE) %>%
  ungroup()

# Summarize relative abundance by group
group_relabund <- relabund_tax %>%
  group_by(Group) %>%
  summarise_if(is.numeric, sum, na.rm = TRUE) %>%
  ungroup() %>%
  mutate(
    MeanRelativeAbundance = rowMeans(select(., where(is.numeric)), na.rm = TRUE),
    Percent = MeanRelativeAbundance * 100
  )

sample_cols <- names(group_relabund)[sapply(group_relabund, is.numeric)]
sample_cols <- setdiff(sample_cols, c("MeanRelativeAbundance", "Percent"))
relabund_long <- group_relabund %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "SampleID",
    values_to = "RelAbundance"
  )

sample_meta <- as.data.frame(sample_data(phyloseq2))
sample_meta$SampleID <- rownames(sample_meta)
relabund_long_meta <- relabund_long %>%
  left_join(sample_meta, by = "SampleID")

desired_order <- c("Rototekoiti", "Manapouri", "Chalice", "Sheppard", "Poerua", "Heaton", "Pooled lakes A", "Pooled lakes B")

relabund_long_meta$Lake <- factor(relabund_long_meta$Lake, levels = desired_order)
relabund_long_meta$SampleID <- factor(relabund_long_meta$SampleID,
                                      levels = unique(relabund_long_meta$SampleID))

# Barplot of groups found in each lake (zoopl, phytopl,...)
colourCAW <- c(
  "Arachnids" = "#c3ce58",         
  "Annelids" = "#f48521",            
  "Clitellata" = "#c45b28",        
  "Crayfish" = "#d383b7",           
  "Fish" = "#ffe0ae",               
  "Fungi" = "#a184bd",  
  "Hydrozoa" = "#fbb15d",
  "Insects" = "#c3d3c2",          
  "Malacostraca" = "#7f1424",      
  "Macro- & Greenalgae" = "#67bd45",
  "Molluscs" = "#9b99cd",          
  "Phytoplankton" = "#029cbd",       
  "Protists" = "#3398d2",            
  "Sponges" = "#00549e",             
  "Water mites" = "#4CFFFF",         
  "Unclassified Arthropoda" = "#76c044",
  "Other" = "#A0A0A0",
  "Zooplankton" = "#4ba791"
)
relabund_long_meta$Group <- factor(
  relabund_long_meta$Group,
  levels = names(colourCAW)
)

ggplot(relabund_long_meta, aes(x = Lake, y = RelAbundance*100, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_fill_manual(values = colourCAW) +
  ylab("Relative abundance [%]") +
  xlab("Lake") +
  ggtitle("Rarefied & grouped taxa SSUF04/R22") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

write.csv(group_counts, "Rarefied_Read_counts_per_group_AllTaxa.csv", row.names = FALSE)
write.csv(group_relabund, "Rarefied_Rel_abundances_per_group_AllTaxa.csv", row.names = FALSE)

### If wanted, here is a table what taxa are in Others
#sample_cols <- names(relabund_tax)[sapply(relabund_tax, is.numeric)]
#relabund_long_taxa <- relabund_tax %>%
# pivot_longer(
#   cols = all_of(sample_cols),
#   names_to = "SampleID",
#   values_to = "RelAbundance"
# ) %>%
# left_join(sample_meta, by = "SampleID")
#other_taxa <- relabund_long_taxa %>%
# filter(Group == "Other")

##_____________________________________________________________________________________________


#Write a table with rarefied read counts for all taxa to lowest resolved rank
tax <- as.data.frame(tax_table(phyloseq2), stringsAsFactors = FALSE)
tax[is.na(tax)] <- ""
rank_order <- rev(rank_names(phyloseq2))  # specific → general
unclassified_patterns <- c("uncultured", "unidentified", "unassigned",
                           "metazoa", "environmental", "sp\\.?$", "NA", "unknown")
assign_lowest <- function(row) {
  vals <- row[rank_order]
  is_valid <- sapply(vals, function(v) {
    if (v == "" || is.na(v)) return(FALSE)
    # check against "junk" patterns
    any(grepl(paste(unclassified_patterns, collapse="|"), v, ignore.case=TRUE)) == FALSE
  })
  
  idx <- which(is_valid)
  if (length(idx) == 0) {
    return(c(Rank="Unclassified", Label="Unclassified_all"))
  } else {
    best <- idx[1]
    rankname <- rank_order[best]
    label <- paste0(rankname, ":", vals[best])
    return(c(Rank=rankname, Label=label))
  }
}

assigned <- t(apply(tax[, rank_order, drop=FALSE], 1, assign_lowest))
assigned_df <- as.data.frame(assigned, stringsAsFactors=FALSE)
rownames(assigned_df) <- rownames(tax)
otu_mat <- as.data.frame(otu_table(phyloseq2))
if (!taxa_are_rows(phyloseq2)) otu_mat <- t(otu_mat)
otu_mat <- otu_mat[rownames(assigned_df), , drop=FALSE]
otu_by_label <- rowsum(otu_mat, group=assigned_df$Label)
agg_df <- as.data.frame(t(otu_by_label))
agg_df <- tibble::rownames_to_column(agg_df, var="Sample")

#Verify totals
row_totals <- rowSums(agg_df[,-1])
stopifnot(all(row_totals == 11000))

write.csv(agg_df, "Rarefied_counts_lowest_resolved.csv", row.names=FALSE)
write.csv(assigned_df, "Taxon_assigned_rank_mapping.csv", row.names = TRUE)


##Write a table only for zooplankton rarefied counts to lowest resolved rank
df <- read.csv("Rarefied_counts_lowest_resolved.csv", check.names = FALSE)
assigned_df <- read.csv("Taxon_assigned_rank_mapping.csv", row.names = 1, check.names = FALSE)
tax <- as.data.frame(tax_table(phyloseq2), stringsAsFactors = FALSE)
tax[is.na(tax)] <- ""
get_higher_taxa <- function(label, assigned_df, tax) {
  otus <- rownames(assigned_df)[assigned_df$Label == label]
  tax_sub <- tax[otus, , drop = FALSE]
  return(tax_sub)
}
zoop_groups <- c("Hexanauplia", "Eurotatoria", "Branchiopoda")
zoop_labels <- sapply(colnames(df)[-1], function(label) {
  tax_sub <- get_higher_taxa(label, assigned_df, tax)
  any(zoop_groups %in% unlist(tax_sub))
})
zoop_labels <- names(zoop_labels[zoop_labels])
df_zoop <- df %>%
  select(Sample, all_of(zoop_labels))
df_zoop$Total_Zooplankton <- rowSums(df_zoop[,-1])

write.csv(df_zoop, "Rarefied_counts_Zooplankton.csv", row.names = FALSE)

##__________________________________________________________________________________________




###You are using the rarefied data now with phyloseq2
Phylum1 <- prune_samples(sample_sums(phyloseq2) > 0,phyloseq2)
#Phylum <- transform_sample_counts(Phylum1, function(x) x/sum(x))
glom <- tax_glom(Phylum1, taxrank = 'phylum', NArm = FALSE)
Phylumdat <- psmelt(glom)
Phylumdat$phylum <- as.character(Phylumdat$phylum)
medians <- ddply(Phylumdat, ~phylum, function(x) c(median=median(x$Abundance)))
remainder <- medians[medians$median <= 0.001,]$phylum
Phylumdat[Phylumdat$phylum %in% remainder,]$phylum <- 'zz_Other'

Phylum2 <- ggplot(Phylumdat, aes(x = Sample_Name, y = Abundance, fill = phylum))  +
  geom_bar(stat = "identity", position = "stack") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.~Lake, scales = "free", space = "free") + 
  xlab("Sample") +
  scale_fill_manual(values= get_pal("caw_cat_1"))+
  scale_y_continuous(expand = c(0,0)) + my_theme

Phylum2

#Phylum2 is now a Phylum level plot without the blank contamination from curated phyloseq2 object



##______________________________________________________________


###Plot different Taxa levels

Phylumfig <- plot_bar(subset_taxa(phyloseq2, phylum %in% c("Arthropoda", "Rotifera", "Cnidaria")), 
                      x = "Sample_Name", fill = "phylum") + 
  facet_grid(.~Lake, scales = "free") +
  geom_bar(aes(fill=phylum), stat="identity", position="stack") + 
  scale_fill_manual(values= colourCAW)+
  scale_y_continuous(expand = c(0,0)) + my_theme 
Phylumfig

##Class
phyloseq_rel <- transform_sample_counts(phyloseq2, function(x) x / sum(x))
Classfig <- plot_bar(subset_taxa(phyloseq_rel, class %in% c("Branchiopoda", "Hexanauplia", "Eurotatoria")), 
                     x = "Sample_Name", fill = "class") + 
  facet_grid(.~Lake, scales = "free") +
  geom_bar(aes(fill=class), stat="identity", position="stack") + 
  scale_fill_manual(values= colourCAW)+
  scale_y_continuous(expand = c(0,0)) + my_theme  +
  labs(y = "Relative Abundance to total taxa detected [%]", x = "")
Classfig

##Genus
phyloseq_rel <- transform_sample_counts(phyloseq2, function(x) x / sum(x))
Genusfig <- plot_bar(subset_taxa(phyloseq_rel, class %in% c("Branchiopoda", "Hexanauplia", "Eurotatoria")), 
                     x = "Sample_Name", fill = "genus") + 
  facet_grid(phylum~Lake, scales = "free") +
  geom_bar(aes(fill=genus), stat="identity", position="stack") + 
  scale_fill_manual(values= colourCAW)+
  scale_y_continuous(expand = c(0,0)) + my_theme +
  labs(y = "Relative Abundance to total taxa detected [%]", x = "") 
Genusfig

##Species (if available)
phyloseq_rel <- transform_sample_counts(phyloseq2, function(x) x / sum(x))
Speciesfig <- plot_bar(subset_taxa(phyloseq_rel, class %in% c("Branchiopoda", "Hexanauplia", "Eurotatoria")), 
                       x = "Sample_Name", fill = "species") + 
  facet_grid(phylum~Lake, scales = "free") +
  geom_bar(aes(fill=species), stat="identity", position="stack") + 
  scale_fill_manual(values= colourCAW)+
  scale_y_continuous(expand = c(0,0)) + my_theme + 
  labs(y = "Relative Abundance to total taxa detected [%]", x = "")
Speciesfig


##_______________________________________________________________________________________________________________
### Make a Genus plot with Rotifera, Cladocera and Copepoda separated
target_classes <- c("Branchiopoda", "Hexanauplia", "Eurotatoria")
phyloseq_rel <- transform_sample_counts(phyloseq2, function(x) x / sum(x))
ps_subset_rel <- subset_taxa(phyloseq_rel, class %in% target_classes)
n_genera <- length(unique(tax_table(ps_subset_rel)[, "genus"]))
desired_lake_order <- c("Rototekoiti", "Manapouri", "Chalice", "Sheppard", "Poerua", "Heaton", "Pooled lakes A", "Pooled lakes B")
sample_data(ps_subset_rel)$Lake <- factor(sample_data(ps_subset_rel)$Lake, 
                                          levels = desired_lake_order)
sample_data(ps_subset_rel)$Sample_Name <- factor(sample_data(ps_subset_rel)$Sample_Name, 
                                                 levels = unique(sample_data(ps_subset_rel)$Sample_Name))

Genusfig_relzoopl <- plot_bar(
  ps_subset_rel,
  x = "Sample_Name",
  fill = "genus"
) +
  facet_grid(class ~ Lake, scales = "free") +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colourCAW) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    title = "Relative Abundance of Genera in Branchiopoda, Hexanauplia, and Eurotatoria SSUF04R22",
    y = "Relative Abundance [%]",
    x = "Sample"
  )

print(Genusfig_relzoopl)

#______________________________________________________________________________________________


### Make a Species (and Genus where not available) plot with Rotifera, Cladocera and Copepoda separated
tax <- as.data.frame(as.matrix(tax_table(phyloseq2)))
gen_col <- grep("genus", colnames(tax), ignore.case = TRUE, value = TRUE)[1]
sp_col  <- grep("species", colnames(tax), ignore.case = TRUE, value = TRUE)[1]

tax$Label <- ifelse(
  is.na(tax[[sp_col]]) | tax[[sp_col]] == "",
  tax[[gen_col]],
  paste(tax[[gen_col]], tax[[sp_col]], sep = " ")
)

tax_table(phyloseq2) <- as.matrix(tax)
target_classes <- c("Branchiopoda", "Hexanauplia", "Eurotatoria")
phyloseq_rel <- transform_sample_counts(phyloseq2, function(x) x / sum(x))
ps_subset_rel <- subset_taxa(phyloseq_rel, class %in% target_classes)
desired_lake_order <- c("Rototekoiti", "Manapouri", "Chalice", "Sheppard",
                        "Poerua", "Heaton", "Pooled lakes A", "Pooled lakes B")

sample_data(ps_subset_rel)$Lake <- factor(
  sample_data(ps_subset_rel)$Lake,
  levels = desired_lake_order
)

sample_data(ps_subset_rel)$Sample_Name <- factor(
  sample_data(ps_subset_rel)$Sample_Name,
  levels = unique(sample_data(ps_subset_rel)$Sample_Name)
)

SpeciesGenusFig_relzoopl <- plot_bar(
  ps_subset_rel,
  x = "Sample_Name",
  fill = "Label"
) +
  facet_grid(class ~ Lake, scales = "free") +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colourCAW) +
  scale_y_continuous(
    expand = c(0, 0),
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    title = "Relative Abundance of Zooplankton SSUF04R22 (Best Available Taxonomy: Genus/Species)",
    y = "Relative Abundance [%]",
    x = "Sample"
  )

print(SpeciesGenusFig_relzoopl)








# Write a csv with relative abundances (not rarefied) for species and genus levels for 3 zoopl classes

tax <- as.data.frame(as.matrix(tax_table(phyloseq2)))
gen_col <- grep("genus", colnames(tax), ignore.case = TRUE, value = TRUE)[1]
sp_col  <- grep("species", colnames(tax), ignore.case = TRUE, value = TRUE)[1]

tax$Label <- ifelse(
  is.na(tax[[sp_col]]) | tax[[sp_col]] == "",
  tax[[gen_col]],
  paste(tax[[gen_col]], tax[[sp_col]], sep = " ")
)

tax_table(phyloseq2) <- as.matrix(tax)
target_classes <- c("Branchiopoda", "Eurotatoria", "Hexanauplia")
physeq_filtered <- subset_taxa(phyloseq2, class %in% target_classes)
physeq_rel <- transform_sample_counts(physeq_filtered, function(x) x / sum(x))
physeq_label <- tax_glom(physeq_rel, taxrank = "Label")

#replace OTU names with Label names
tax_labels <- as.data.frame(tax_table(physeq_label))$Label
tax_labels[is.na(tax_labels) | tax_labels == ""] <- "Unclassified"
tax_labels <- make.unique(tax_labels)  # avoid duplicate column names

# Assign Label names as OTU names
taxa_names(physeq_label) <- tax_labels
abund <- as.data.frame(otu_table(physeq_label))
if (taxa_are_rows(physeq_label)) {
  abund <- t(abund)
}
meta <- as.data.frame(sample_data(physeq_label))
abund_df <- cbind(meta, abund)

write.csv(abund_df, "RelAbundance_Genus_Species_SelectedClasses.csv", row.names = TRUE)


