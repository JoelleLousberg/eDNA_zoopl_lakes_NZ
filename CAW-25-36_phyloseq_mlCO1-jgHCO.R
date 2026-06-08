rm(list= ls())

getwd()

setwd("C:/Users/joelle.lousberg/OneDrive - Cawthron/Documents/Chapter_2/R_step2/mlCO1_all")


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

###Use MEGAN COMMUNITY EDITION - to do LCA (lowest common ancestor analysis)

blast.tax <- fread("mlCO1_all.nochim_ASVs_blast_BOLD-ex2.txt", header = FALSE)
head(blast.tax)
blast.tax2 = colsplit(blast.tax$V2, ';', c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species")) #since taxonomy table only consist of 2 columns #(feature ID and taxa) the 7 hierarchies of taxa need to be put in different columns.  string: Taxon (header of column that needs to be split up); pattern ;
row.names(blast.tax2) = blast.tax$`V1`
blast.tax2 <- blast.tax2 %>% mutate_all(na_if,"")
write.csv(blast.tax2, "MEGANCO1.blast.csv")


CO1.blast.tax.edit<- blast.tax2
CO1.blast.tax.edit$Row.names <- rownames(CO1.blast.tax.edit)
head(CO1.blast.tax.edit)

filled_tax_df <- CO1.blast.tax.edit %>%
  mutate(
    Kingdom = ifelse(is.na(Kingdom), Domain, Kingdom),
    Phylum = ifelse(is.na(Phylum), Kingdom, Phylum),
    Class = ifelse(is.na(Class), Phylum, Class),
    Order = ifelse(is.na(Order), Class, Order),
    Family = ifelse(is.na(Family), Order, Family),
    Genus = ifelse(is.na(Genus), Family, Genus),
    Species = ifelse(is.na(Species), Genus, Species)
  )

head(filled_tax_df)

ASVCO1 <- fread("mlCO1_all.nochim_ASVs_counts.tsv", header = T)
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

ps_tax = tax_table(as.matrix(filled_tax_df))
head(ps_tax)
#adding ID's back
rownames(ps_tax)=rownames(filled_tax_df)
colnames(ps_tax)=colnames(filled_tax_df)
head(ps_tax)

####
sample.data <-fread("C:/Users/joelle.lousberg/OneDrive - Cawthron/Documents/Chapter_2/R_step2/mlCO1_all/CAW-25-36_sample_data.csv")
names(sample.data)

sample.data$Sample.Name <- as.factor(sample.data$Sample_Name)

ps_map = sample_data(sample.data)
head(ps_map)
rownames(ps_map)=sample.data$'Sample_ID'
head(ps_map)




##new phyloseq object to create plots

phyloseq_mlCO1 <- phyloseq(ps_otu_table, ps_map, ps_tax)
phyloseq_mlCO1



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



##First look at data
phyloseq1 <- phyloseq_mlCO1
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


####Subtraction of contamination from all lab blanks
Controls = subset_samples(phyloseq1, Sample_type %in% c("Blank"))
Controls = filter_taxa(Controls, function(x) sum(x) > 0, TRUE)
sample_sums(Controls)
taxa_sums(Controls)
#You can change the taxon level here
plot_bar(Controls, fill = "Genus")

#substract NC
Extraction_neg <- subset_samples(phyloseq1, Sample_type %in% c("Blank"))
Extraction_neg
Extraction_neg_max <- apply(as.matrix(otu_table(Extraction_neg)), 1, max)
Extraction_neg_max_vec <- as.vector(Extraction_neg_max)
names(Extraction_neg_max_vec) <- taxa_names(Extraction_neg)  # ensure names stay aligned
Extraction <- as(otu_table(phyloseq1), "matrix")

# Ensure taxa are rows
if (!taxa_are_rows(phyloseq1)) {
  Extraction <- t(Extraction)
}

Extractiondf <- as.data.frame(Extraction)
Extractiondf <- Extractiondf[names(Extraction_neg_max_vec), , drop = FALSE]

Extractiondf <- sweep(Extractiondf, 1, Extraction_neg_max_vec, "-")
Extractiondf[Extractiondf < 0] <- 0

rownames(Extractiondf) <- taxa_names(phyloseq1)
colnames(Extractiondf) <- sample_names(phyloseq1)

stopifnot(identical(rownames(Extractiondf), taxa_names(phyloseq1)))
stopifnot(identical(colnames(Extractiondf), sample_names(phyloseq1)))

Extraction_OTU <- otu_table(as.matrix(Extractiondf), taxa_are_rows = TRUE)

Tutorial.subtractextract.ps <- phyloseq(
  Extraction_OTU,
  sample_data(phyloseq1),
  tax_table(phyloseq1)
)



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
                         !Phylum %in% c("d__Eukaryota", "k__Metazoa") & 
                           !is.na(Phylum))


phyloseq_mlCO1_clean <- phyloseq2
saveRDS(phyloseq_mlCO1_clean, file = "phyloseq_mlCO1_clean.rds")

###Load from here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
phyloseq2 <- readRDS("phyloseq_mlCO1_clean.rds")


###Rarefaction step
library(vegan)
sample_sums(phyloseq2)
#min(sample_sums(phyloseq3))
sort(sample_sums(phyloseq2))[1:10]
set.seed(100)
phyloseq2_rare = rarefy_even_depth(phyloseq2, sample.size = 11000, 
                                   replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
phyloseq2 <- phyloseq2_rare
#_________________________________________________________________________________________



### Write a table with all counts as percentages relative abundance but grouped into zooplankton etc (NOT per Species!)
if (!taxa_are_rows(phyloseq2)) {
  phyloseq2 <- transform_sample_counts(phyloseq2, identity)
}
otu_counts <- as.data.frame(otu_table(phyloseq2))
otu_counts$Taxon <- rownames(otu_counts)
tax_df <- as.data.frame(tax_table(phyloseq2))
tax_df$Taxon <- rownames(tax_df)
counts_tax <- merge(otu_counts, tax_df, by = "Taxon")

phyloseq_relabund <- transform_sample_counts(phyloseq2, function(x) x / sum(x))
otu_relabund <- as.data.frame(otu_table(phyloseq_relabund))
otu_relabund$Taxon <- rownames(otu_relabund)
relabund_tax <- merge(otu_relabund, tax_df, by = "Taxon")

relabund_tax <- relabund_tax %>%
  mutate(across(c(Phylum, Class, Order, Family, Genus), ~trimws(.)),
         Phylum = coalesce(Phylum, ""),
         Class  = coalesce(Class, ""),
         Order  = coalesce(Order, ""),
         Family = coalesce(Family, ""),
         Genus  = coalesce(Genus, ""))

relabund_tax <- relabund_tax %>%
  mutate(Group = case_when(
    Genus %in% c("g__Dinobryon", "g__Trachydiscus", "g__Lindavia", "g__Stephanopyxis") ~ "Phytoplankton",
    Genus %in% c("g__Phytophthora", "g__Pythium", "g__Aphanomyces") ~ "Fungi",
    Genus %in% c("g__Pugetia", "g__Padina", "g__Stictyosiphon", "g__Macrocystis", "g__Dictyota") ~ "Macro- & Greenalgae",
    Genus == "g__Acanthamoeba" ~ "Protists",
    Genus == "g__Alcelaphus" ~ "Fish",
    Genus == "g__Naegleria" ~ "Protists",
    
    Family %in% c("f__Saprolegniaceae", "f__Pythiaceae", "f__Acinetosporaceae") ~ "Fungi",
    Family %in% c("f__Durvillaeaceae", "f__Dictyotaceae", "f__Scytosiphonaceae", "f__Chordariaceae") ~ "Macro- & Greenalgae",
    Family %in% c("f__Gymnodiniaceae", "f__Suessiaceae", "f__Goniochloridaceae", "f__Dinobryaceae", "f__Mallomonadaceae") ~ "Phytoplankton",
    
    Order %in% c("o__Melosirales", "o__Coscinodiscales", "o__Naviculales",
                 "o__Lithodesmiales", "o__Rhizosoleniales", "o__Bacillariales",
                 "o__Fragilariales", "o__Stephanodiscales", "o__Cymbellales") ~ "Phytoplankton",
    Order %in% c("o__Ceramiales", "o__Gigartinales", "o__Ectocarpales", "o__Laminariales") ~ "Macro- & Greenalgae",
    Order %in% c("o__Decapoda") ~ "Crayfish",
    Order %in% c("o__Trombidiformes") ~ "Water mites",
    Order %in% c("o__Haplotaxida") ~ "Annelids",
    Order %in% c("o__Spongillida") ~ "Sponges",
    Order %in% c("o__Dactylopodida") ~ "Protists",
    
    Class %in% c("c__Actinopteri") ~ "Fish",
    Class %in% c("c__Hexanauplia", "c__Branchiopoda", "c__Eurotatoria") ~ "Zooplankton",
    Class %in% c("c__Bivalvia", "c__Gastropoda") ~ "Molluscs",
    Class %in% c("c__Florideophyceae", "c__Phaeophyceae", "c__Chlorophyceae",
                 "c__Mamiellophyceae", "c__Trebouxiophyceae", "c__Chloropicophyceae") ~ "Macro- & Greenalgae",
    Class %in% c("c__Hydrozoa", "c__Scyphozoa", "c__Anthozoa") ~ "Hydrozoa",
    Class %in% c("c__Sordariomycetes", "c__Eurotiomycetes", "c__Leotiomycetes",
                 "c__Microbotryomycetes", "c__Exobasidiomycetes", "c__Arthoniomycetes",
                 "c__Pichiomycetes") ~ "Fungi",
    Class %in% c("c__Clitellata", "c__Catenulida") ~ "Annelids",
    Class %in% c("c__Malacostraca") ~ "Malacostraca",
    Class %in% c("c__Insecta") ~ "Insects",
    Class %in% c("c__Demospongiae") ~ "Sponges",
    Class %in% c("c__Arachnida") ~ "Arachnids",
    
    Phylum == "p__Arthropoda" & Class == "" ~ "Unclassified Arthropoda",
    Phylum == "p__Bacillariophyta" ~ "Phytoplankton",
    Phylum == "p__Oomycota" ~ "Fungi",
    
    # Catch-all
    TRUE ~ "Other"
  ))

# Check that grouping worked
table(relabund_tax$Group)

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
  ggtitle("Rarefied & grouped taxa mlCO1/jgHCO") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

write.csv(group_counts, "Rarefied_ReadCounts_grouped_AllTaxa.csv", row.names = FALSE)
write.csv(group_relabund, "Rarefied_Rel_abundances_grouped_AllTaxa.csv", row.names = FALSE)

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





#Write a table with rarefied read counts for all taxa to lowest resolved rank (CO1 only, change the 11,000 if needed)
if (!exists("sample_depth")) sample_depth <- 11000

# 2. Convert tax_table to plain data.frame and clean whitespace / NAs
tax_df <- as.data.frame(tax_table(phyloseq2), stringsAsFactors = FALSE)
tax_df[] <- lapply(tax_df, function(x) {
  x <- as.character(x)
  x <- gsub("^\\s+|\\s+$", "", x)
  x[x == "NA"] <- ""      # some DBs use literal "NA"
  x[is.na(x)] <- ""
  x
})

# 3. Detect & drop columns that duplicate the rownames (common accidentally added column)
rn <- rownames(tax_df)
dup_cols <- sapply(colnames(tax_df), function(col) all(as.character(tax_df[[col]]) == rn))
if (any(dup_cols)) {
  cat("Dropping tax columns that duplicate rownames:\n")
  print(colnames(tax_df)[dup_cols])
  tax_df <- tax_df[, !dup_cols, drop = FALSE]
}
# drop obvious nonsense names if present
bad_names <- c("Row.names", "row.names", "rowname", "X", "X.1")
drop2 <- intersect(colnames(tax_df), bad_names)
if (length(drop2) > 0) {
  cat("Dropping columns named: ", paste(drop2, collapse = ", "), "\n")
  tax_df <- tax_df[, !(colnames(tax_df) %in% drop2), drop = FALSE]
}

# 4. Prepare rank_order (most specific -> more general) and confirm columns present
all_ranks <- rank_names(phyloseq2)
rank_order <- rev(all_ranks)                 # e.g. species, genus, family, order, class, phylum, kingdom, domain
rank_order <- rank_order[rank_order %in% colnames(tax_df)]   # keep only existing columns
cat("Using rank_order (most specific -> general):\n"); print(rank_order)

# 5. Prepare helper functions: clean tax strings and prefix -> rank mapping
prefix_map <- c(
  "s" = "species", "sp" = "species",
  "g" = "genus",
  "f" = "family",
  "o" = "order",
  "c" = "class",
  "p" = "phylum",
  "k" = "kingdom",
  "d" = "domain"
)
unclassified_patterns <- c("uncultured","unidentified","unassigned","metazoa",
                           "environmental","sp\\.?$","^NA$","^unknown$")

clean_tax_string <- function(x) {
  if (is.na(x) || x == "") return("")
  # remove prefix (g__, s__, etc.), trailing semicolons, and trim whitespace
  y <- sub("^[A-Za-z]+__", "", x)
  y <- sub(";.*$", "", y)
  y <- gsub("^\\s+|\\s+$", "", y)
  y
}

# 6. Build cleaned tax table (no prefixes) for lookup
clean_tax_df <- tax_df
for (col in colnames(tax_df)) {
  clean_tax_df[[col]] <- sapply(tax_df[[col]], clean_tax_string, USE.NAMES = FALSE)
}

# 7. Assignment function: choose first valid entry (most specific → general),
#    infer true rank from prefix if present; strip prefix in label
assign_lowest_with_prefix <- function(otu) {
  vals_raw <- unlist(tax_df[otu, rank_order], use.names = FALSE)
  valid <- sapply(vals_raw, function(v) {
    if (is.null(v) || is.na(v) || v == "") return(FALSE)
    !any(grepl(paste(unclassified_patterns, collapse = "|"), v, ignore.case = TRUE))
  })
  if (!any(valid)) {
    return(data.frame(Rank = "ASV", Label = otu, stringsAsFactors = FALSE))
  }
  best_idx <- which(valid)[1]
  rawval <- vals_raw[best_idx]
  prefix_match <- regmatches(rawval, regexpr("^([A-Za-z]+)__", rawval))
  if (length(prefix_match) && nchar(prefix_match) > 0) {
    prefix_code <- sub("__", "", prefix_match)
    rankname_from_prefix <- prefix_map[[tolower(prefix_code)]]
    rankname <- if (!is.null(rankname_from_prefix)) rankname_from_prefix else rank_order[best_idx]
  } else {
    rankname <- rank_order[best_idx]
  }
  taxclean <- clean_tax_string(rawval)
  label <- paste0(otu, "|", rankname, ":", taxclean)
  return(data.frame(Rank = rankname, Label = label, stringsAsFactors = FALSE))
}

# 8. Apply assignment to all OTUs (ASVs)
otu_names <- rownames(tax_df)
assigned_list <- lapply(otu_names, assign_lowest_with_prefix)
assigned_df <- do.call(rbind, assigned_list)
rownames(assigned_df) <- otu_names

# diagnostic: check for any leftover "Row.names" in labels
bad_row_names_count <- sum(grepl("Row.names", assigned_df$Label, fixed = TRUE))
cat("Remaining labels containing 'Row.names':", bad_row_names_count, "\n")
if (bad_row_names_count > 0) {
  cat("Example problematic mappings (first 20):\n")
  print(head(assigned_df[grepl("Row.names", assigned_df$Label, fixed = TRUE), , drop = FALSE], 20))
  warning("Some labels contain 'Row.names' — inspect tax_df and column names.")
}

cat("Assigned rank counts:\n"); print(sort(table(assigned_df$Rank), decreasing = TRUE))

# 9. Aggregate OTU counts by the assigned Label
otu_mat <- as.data.frame(otu_table(phyloseq2))
if (!taxa_are_rows(phyloseq2)) otu_mat <- t(otu_mat)
otu_mat <- otu_mat[rownames(assigned_df), , drop = FALSE]

otu_by_label <- rowsum(otu_mat, group = assigned_df$Label)
agg_df <- as.data.frame(t(otu_by_label))
agg_df <- tibble::rownames_to_column(agg_df, var = "Sample")

# verify totals (if rarefied)
if (!is.na(sample_depth)) {
  totals_ok <- all(rowSums(agg_df[,-1, drop = FALSE]) == sample_depth)
  cat("Totals per sample equal rarefied depth?", totals_ok, "\n")
  if (!totals_ok) warning("Per-sample totals after aggregation are not all equal to sample_depth.")
}

# 10. Save mapping + aggregated table
write.csv(assigned_df, "Taxon_assigned_rank_mapping_CO1.csv", row.names = TRUE)
write.csv(agg_df, "Rarefied_counts_CO1_AllTaxa.csv", row.names = FALSE)

# 11. Extract zooplankton labels (matches at any rank using cleaned taxonomy)
label_names <- colnames(agg_df)[-1]
label_to_otus <- lapply(label_names, function(lbl) {
  rownames(assigned_df)[assigned_df$Label == lbl]
})
names(label_to_otus) <- label_names

zoop_groups <- c("Hexanauplia","Branchiopoda","Eurotatoria")  # adjust case if needed
zoop_labels_found <- label_names[sapply(label_names, function(lbl) {
  otus <- label_to_otus[[lbl]]
  if (length(otus) == 0) return(FALSE)
  any(tolower(as.matrix(clean_tax_df[otus, , drop = FALSE])) %in% tolower(zoop_groups))
})]

cat("Number of labels identified as zooplankton:", length(zoop_labels_found), "\n")
cat("Example zooplankton labels (first 30):\n"); print(head(zoop_labels_found, 30))

# 12. Build zooplankton-only table
df_zoop <- agg_df %>% select(Sample, all_of(zoop_labels_found))
orig_zoop_cols <- colnames(df_zoop)[-1]

# clean headers for zooplankton
new_zoop_names <- sapply(orig_zoop_cols, function(lbl) {
  if (grepl("\\|", lbl)) {
    taxpart <- sub(".*\\|", "", lbl)
    taxonly <- sub("^[^:]*:", "", taxpart)
    taxonly <- gsub("^\\s+|\\s+$", "", taxonly)
    if (nzchar(taxonly)) return(taxonly) else return(lbl)
  } else {
    return(lbl)
  }
}, USE.NAMES = FALSE)
new_zoop_names_unique <- make.unique(new_zoop_names)
colnames(df_zoop) <- c("Sample", new_zoop_names_unique)
df_zoop$Total_Zooplankton <- rowSums(df_zoop[ , -1, drop = FALSE])

# Save zooplankton tables and mapping
write.csv(df_zoop, "Rarefied_counts_Zooplankton_CO1.csv", row.names = FALSE)
zoop_label_map <- data.frame(
  OriginalLabel = orig_zoop_cols,
  NewHeader = new_zoop_names_unique,
  stringsAsFactors = FALSE
)
write.csv(zoop_label_map, "Zoop_label_mapping.csv", row.names = FALSE)
write.csv(df_zoop, "Rarefied_counts_Zooplankton_CO1.csv", row.names = FALSE)

# 13. Replace ASV headers by cleaned tax names for all samples
new_col_names <- sapply(colnames(agg_df)[-1], function(lbl) {
  if (grepl("\\|", lbl)) {
    taxpart <- sub(".*\\|", "", lbl)
    taxonly <- sub("^[^:]*:", "", taxpart)
    taxonly <- gsub("^\\s+|\\s+$", "", taxonly)
    if (nzchar(taxonly)) return(taxonly) else return(lbl)
  } else {
    return(lbl)
  }
}, USE.NAMES = FALSE)
new_col_names_unique <- make.unique(new_col_names)
agg_df_taxheaders <- agg_df
colnames(agg_df_taxheaders) <- c("Sample", new_col_names_unique)
write.csv(agg_df_taxheaders, "Rarefied_counts_CO1_AllTaxa.csv", row.names = FALSE)

# ------------- done -------------






### Basic plotting of either rarefied or non rarefied data (skip rarefaction step above if not needed!!!)

#First Phylum level plot without blank contamination from curated phyloseq2 object
Phylum1 <- prune_samples(sample_sums(phyloseq2) > 0,phyloseq2)
#Phylum <- transform_sample_counts(Phylum1, function(x) x/sum(x))
glom <- tax_glom(Phylum1, taxrank = 'Phylum', NArm = FALSE)
Phylumdat <- psmelt(glom)
Phylumdat$Phylum <- as.character(Phylumdat$Phylum)
medians <- ddply(Phylumdat, ~Phylum, function(x) c(median=median(x$Abundance)))
remainder <- medians[medians$median <= 0.001,]$Phylum
Phylumdat[Phylumdat$Phylum %in% remainder,]$Phylum <- 'zz_Other'

Phylum2 <- ggplot(Phylumdat, aes(x = Sample_Name, y = Abundance, fill = Phylum))  +
  geom_bar(stat = "identity", position = "stack") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.~Lake, scales = "free", space = "free") + 
  xlab("Sample") +
  scale_fill_manual(values= get_pal("caw_cat_1"))+
  scale_y_continuous(expand = c(0,0)) + my_theme

Phylum2

##______________________________________________________________


##Basic Plot for different Taxa levels

Phylumfig <- plot_bar(subset_taxa(phyloseq2, Phylum %in% c("p__Arthropoda", "p__Rotifera", "p__Cnidaria")), 
                      x = "Sample_Name", fill = "Phylum") + 
  facet_grid(.~Lake, scales = "free") +
  geom_bar(aes(fill=Phylum), stat="identity", position="stack") + 
  scale_fill_manual(values= colourCAW)+
  scale_y_continuous(expand = c(0,0)) + my_theme 
Phylumfig

##Class
phyloseq_rel <- transform_sample_counts(phyloseq2, function(x) x / sum(x))
Classfig <- plot_bar(subset_taxa(phyloseq_rel, Class %in% c("c__Branchiopoda", "c__Hexanauplia", "c__Eurotatoria")), 
                     x = "Sample_Name", fill = "Class") + 
  facet_grid(.~Lake, scales = "free") +
  geom_bar(aes(fill=Class), stat="identity", position="stack") + 
  scale_fill_manual(values= colourCAW)+
  scale_y_continuous(expand = c(0,0)) + my_theme  +
  labs(y = "Relative Abundance to total taxa detected [%]", x = "")
Classfig

##Genus
phyloseq_rel <- transform_sample_counts(phyloseq2, function(x) x / sum(x))
Genusfig <- plot_bar(subset_taxa(phyloseq_rel, Class %in% c("c__Branchiopoda", "c__Hexanauplia", "c__Eurotatoria")), 
                     x = "Sample_Name", fill = "Genus") + 
  facet_grid(Phylum~Lake, scales = "free") +
  geom_bar(aes(fill=Genus), stat="identity", position="stack") + 
  scale_fill_manual(values= colourCAW)+
  scale_y_continuous(expand = c(0,0)) + my_theme +
  labs(y = "Relative Abundance to total taxa detected [%]", x = "") 
Genusfig

##Species (if available)
phyloseq_rel <- transform_sample_counts(phyloseq2, function(x) x / sum(x))
Speciesfig <- plot_bar(subset_taxa(phyloseq_rel, Class %in% c("c__Branchiopoda", "c__Hexanauplia", "c__Eurotatoria")), 
                       x = "Sample_Name", fill = "Species") + 
  facet_grid(Phylum~Lake, scales = "free") +
  geom_bar(aes(fill=Species), stat="identity", position="stack") + 
  scale_fill_manual(values= colourCAW)+
  scale_y_continuous(expand = c(0,0)) + my_theme + 
  labs(y = "Relative Abundance to total taxa detected [%]", x = "")
Speciesfig

##_______________________________________________________________________________________________________________



### Genus plot with Rotifera, Cladocera and Copepoda separated
target_classes <- c("c__Branchiopoda", "c__Hexanauplia", "c__Eurotatoria")
phyloseq_rel <- transform_sample_counts(phyloseq2, function(x) x / sum(x))
ps_subset_rel <- subset_taxa(phyloseq_rel, Class %in% target_classes)
n_genera <- length(unique(tax_table(ps_subset_rel)[, "Genus"]))
desired_lake_order <- c("Rototekoiti", "Manapouri", "Chalice", "Sheppard", "Poerua", "Heaton", "Pooled lakes A", "Pooled lakes B")
sample_data(ps_subset_rel)$Lake <- factor(sample_data(ps_subset_rel)$Lake, 
                                          levels = desired_lake_order)
sample_data(ps_subset_rel)$Sample_Name <- factor(sample_data(ps_subset_rel)$Sample_Name, 
                                                 levels = unique(sample_data(ps_subset_rel)$Sample_Name))

Genusfig_relzoopl_mlCO1 <- plot_bar(
  ps_subset_rel,
  x = "Sample_Name",
  fill = "Genus"
) +
  facet_grid(Class ~ Lake, scales = "free") +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colourCAW) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    title = "Relative Abundance of Genera in Branchiopoda, Hexanauplia, and Eurotatoria mlCO1",
    y = "Relative Abundance [%]",
    x = "Sample"
  )

print(Genusfig_relzoopl_mlCO1)

#______________________________________________________________________________________________


### Species plot with Rotifera, Cladocera and Copepoda separated
target_classes <- c("c__Branchiopoda", "c__Hexanauplia", "c__Eurotatoria")
phyloseq_rel <- transform_sample_counts(phyloseq2, function(x) x / sum(x))
ps_subset_rel <- subset_taxa(phyloseq_rel, Class %in% target_classes)
n_genera <- length(unique(tax_table(ps_subset_rel)[, "Species"]))
desired_lake_order <- c("Rototekoiti", "Manapouri", "Chalice", "Sheppard", "Poerua", "Heaton", "Pooled lakes A", "Pooled lakes B")
sample_data(ps_subset_rel)$Lake <- factor(sample_data(ps_subset_rel)$Lake, 
                                          levels = desired_lake_order)
sample_data(ps_subset_rel)$Sample_Name <- factor(sample_data(ps_subset_rel)$Sample_Name, 
                                                 levels = unique(sample_data(ps_subset_rel)$Sample_Name))

Speciesfig_relzoopl_mlCO1 <- plot_bar(
  ps_subset_rel,
  x = "Sample_Name",
  fill = "Species"
) +
  facet_grid(Class ~ Lake, scales = "free") +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colourCAW) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    title = "Relative Abundance of Species in Branchiopoda, Hexanauplia, and Eurotatoria mlCO1",
    y = "Relative Abundance [%]",
    x = "Sample"
  )

print(Speciesfig_relzoopl_mlCO1)

#______________________________________________________________________________________________



# CSV with relative abundances for species and genus levels only for 3 zoopl classes

target_classes <- c("c__Branchiopoda", "c__Eurotatoria", "c__Hexanauplia")
phyloseq_filtered <- subset_taxa(phyloseq2, Class %in% target_classes)
phyloseq_rel <- transform_sample_counts(phyloseq_filtered, function(x) x / sum(x))

phyloseq_genus <- tax_glom(phyloseq_rel, taxrank = "Genus")
genus_abundance <- as.data.frame(t(otu_table(phyloseq_genus)))
tax_tab_genus <- tax_table(phyloseq_genus)
colnames(genus_abundance) <- as.vector(tax_tab_genus[, "Genus"])

phyloseq_species <- tax_glom(phyloseq_rel, taxrank = "Species")
species_abundance <- as.data.frame(t(otu_table(phyloseq_species)))
tax_tab_species <- tax_table(phyloseq_species)
colnames(species_abundance) <- as.vector(tax_tab_species[, "Species"])

metadata <- as.data.frame(sample_data(phyloseq_rel))
metadata$SampleID <- rownames(metadata)

# Combine into one data frame
combined_df <- cbind(metadata,
                     genus_abundance[, colSums(genus_abundance) > 0, drop = FALSE],
                     species_abundance[, colSums(species_abundance) > 0, drop = FALSE])

write.csv(combined_df, "RelAbundance_Genus_Species_SelectedClasses.csv", row.names = FALSE)

