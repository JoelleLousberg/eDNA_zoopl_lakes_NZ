#Diversity analyses

rm(list= ls())

getwd()

setwd("C:/Users/joelle.lousberg/OneDrive - Cawthron/Documents/Chapter_2/R_step3/Diversity")

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggthemes)      # for clean themes
library(ggsci)         # optional, for color palettes
library(ggpubr)        # for statistical annotation
library(multcompView)
library(phyloseq)
library(vegan)
library(tidyverse)
library(lme4)
library(broom)
library(openxlsx)

my_theme = theme_bw(base_size = 12) + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.text.x = element_text(face = "bold", size = 10),
  axis.title.y = element_text(vjust= 2.2),
  axis.text = element_text(colour = "black"),
  axis.text.x = element_text(size = 10, angle = 45, hjust =1),
  legend.position="top",
  legend.text = element_text(size = 12),  
  legend.box = "horizontal",
  legend.key.size = unit(1,"line")) + 
  (theme(plot.margin = unit(c(.65,.65,.65,.65), "cm")))

colourCAW <- c(
  "#FF914C", "#1EA6AC", "#7f1424", "#00549e", "#6860a0", "#4ba791", "#ffcc28","#9b99cd", "#c3ce58", 
  "#702365", "#4CFFFF", "#acd58e", "#7ac4d3", "#f8d3ca", "#2f4926", "#3398d2", "#4a7637", "#175c7d", 
  "#c3d3c2", "#FF6F00", "#c45b28", "#e8b5d4", "#2f725e", "#85243f", "#75a54a", "#bc3635", 
  "#205128", "#35356d", "#083631", "#fbd872", "#c22c43")
method_cols <- c(
  "morphology" = "#76c044",
  "townet"     = "#FF914C",
  "freezer"    = "#1EA6AC",
  "water"      = "#7f1424")


#Rarefy and subset the three phyloseq objects 3 primers

ps_fwh <- readRDS("phyloseq_fwhF1R1_26_06_clean.rds")
# remove extraction blank
ps_fwh <- subset_samples(ps_fwh, Sample_Name != "Pooled_extraction _blanks")

# remove taxa with zero counts after sample removal
ps_fwh <- prune_taxa(taxa_sums(ps_fwh) > 0, ps_fwh)
ps_fwh <- subset_taxa(ps_fwh, Class == "c__Branchiopoda" | Class == "c__Hexanauplia" | Class == "c__Eurotatoria")
set.seed(100)
ps_fwh <- rarefy_even_depth(ps_fwh, sample.size = 15000, 
                            replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
ntaxa(ps_fwh)


ps_mlCO1 <- readRDS("phyloseq_mlCO1_26_07_clean.rds")
ps_mlCO1 <- subset_taxa(ps_mlCO1, Class == "c__Branchiopoda" | Class == "c__Hexanauplia" | Class == "c__Eurotatoria")
set.seed(100)
ps_mlCO1 <- rarefy_even_depth(ps_mlCO1, sample.size = 15000,
                              replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
ntaxa(ps_mlCO1)


ps_Uni18S <- readRDS("phyloseq_Uni18S_26_07_clean.rds")
ps_Uni18S <- subset_taxa(ps_Uni18S, class == "Branchiopoda" | class == "Hexanauplia" | class == "Eurotatoria")
set.seed(100)
ps_Uni18S <- rarefy_even_depth(ps_Uni18S, sample.size = 15000, 
                               replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
ntaxa(ps_Uni18S)


#Convert read counts to relative abundances
to_relabund <- function(ps) {
  transform_sample_counts(ps, function(x) x / sum(x))
}

ps_fwh_rel <- to_relabund(ps_fwh)
ps_mlCO1_rel <- to_relabund(ps_mlCO1)
ps_Uni18S_rel <- to_relabund(ps_Uni18S)


#Create a function that takes the mean of the three water replicates for diversity analyses
average_water_reps <- function(ps) {

  meta <- as(sample_data(ps), "data.frame")
  meta$SampleID <- rownames(meta)
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) {
    otu <- t(otu)}
  otu_df <- as.data.frame(otu)
  otu_df$SampleID <- rownames(otu_df)
  merged <- dplyr::left_join(meta, otu_df, by = "SampleID")
  # rename water replicates
  merged$Collection_Method <- ifelse(
    merged$Collection_Method %in% c("water1","water2","water3"),
    "water",
    merged$Collection_Method)
  asv_cols <- colnames(otu_df)[colnames(otu_df) != "SampleID"]
  # average relative abundances
  averaged <- merged %>%
    dplyr::group_by(Lake, Collection_Method) %>%
    dplyr::summarise(
      dplyr::across(all_of(asv_cols), mean),
      .groups = "drop")
  # new sample IDs
  averaged$NewSampleID <- paste(
    averaged$Lake,
    averaged$Collection_Method,
    sep = "_")
  
  # OTU matrix
  otu_new <- averaged[, asv_cols]
  otu_new <- as.matrix(otu_new)
  rownames(otu_new) <- averaged$NewSampleID
  otu_new <- otu_table(otu_new, taxa_are_rows = FALSE)
  # metadata dataframe
  meta_new <- averaged[, c("NewSampleID", "Lake", "Collection_Method")]
  meta_new <- as.data.frame(meta_new)
  rownames(meta_new) <- meta_new$NewSampleID
  meta_new$NewSampleID <- NULL
  meta_new <- sample_data(meta_new)
  # confirm names match
  identical(
    sample_names(otu_new),
    rownames(meta_new))
  # rebuild phyloseq
  phyloseq(otu_new, meta_new, tax_table(ps))
}

#Aplly this on all 3 primers:
ps_fwh_final <- average_water_reps(ps_fwh_rel)
ps_mlCO1_final <- average_water_reps(ps_mlCO1_rel)
ps_Uni18S_final <- average_water_reps(ps_Uni18S_rel)

ps_list_rarefied <- list(
  fwh = ps_fwh_final,
  mlCOI = ps_mlCO1_final,
  Uni18S = ps_Uni18S_final)


#Alpha diversity:
library(phyloseq)
library(dplyr)
library(lmerTest)
library(openxlsx)

alpha_all <- lapply(names(ps_list_rarefied), function(p){
  ps <- ps_list_rarefied[[p]]
  # OTU matrix
  otu <- as(otu_table(ps), "matrix")
  # ensure samples are rows
  if (taxa_are_rows(ps)) {otu <- t(otu)}
  # replace NA with 0
  otu[is.na(otu)] <- 0
  # calculate diversity
  observed <- rowSums(otu > 0)
  shannon <- vegan::diversity(otu, index = "shannon")
  simpson <- vegan::diversity(otu, index = "simpson")
  alpha <- data.frame(SampleID = rownames(otu),Observed = observed,Shannon = shannon,Simpson = simpson)
  alpha$Hill_Shannon <- exp(alpha$Shannon)
  # Hill-Simpson (q = 2)
  alpha$Hill_Simpson <- 1 / (1 - alpha$Simpson)
  meta <- as.data.frame(sample_data(ps))
  df <- cbind(alpha, meta)
  df$Primer <- p
  df
}) %>% bind_rows()

write.xlsx(alpha_all,
           "alpha_diversities_allprimers_alllakes_waterpooled.xlsx",
           rowNames = FALSE)

# Exclude Waihau and Waihola
alpha_all <- alpha_all %>%
  dplyr::filter(!Lake %in% c("Waihau", "Waihola"))
alpha_wide <- alpha_all %>%
  select(Lake, Collection_Method, Primer, Shannon, Simpson, Hill_Shannon) %>%
  pivot_wider(
    names_from = Primer,
    values_from = c(Shannon, Simpson, Hill_Shannon),
    names_glue = "{Primer}_{.value}") %>%
  select(
    Lake, Collection_Method,
    starts_with("fwh_"),
    starts_with("mlCOI_"),
    starts_with("Uni18S_"))
write.csv(alpha_wide, "alpha_diversities_wide_allprimers_alllakes_ASVs_waterpooled.csv", row.names = FALSE)


alpha_model <- lmer(
  Hill_Shannon ~ Collection_Method * Primer + (1|Lake),
  data = alpha_all
)
summary(alpha_model)
#ANOVA tells me if after accounting for collection method and lake, is there a main effect of primer on Hill-Shannon diversity?
anova(alpha_model, type = 3)

library(emmeans)

emm <- emmeans(alpha_model, ~ Collection_Method | Primer)
pairs(emm)

#Posthoc comparison that tells me which primer and for which collection method there is a difference in Hill_Shannon diversity 
emm_primer <- emmeans(alpha_model, ~ Primer)
pairs(emm_primer, adjust = "BH")
emm_primer_by_method <- emmeans(alpha_model, ~ Primer | Collection_Method)
pairs(emm_primer_by_method, adjust = "BH")


ggplot(alpha_all,
       aes(Collection_Method, Hill_Shannon, fill = Collection_Method)) +
  geom_boxplot() +
  facet_wrap(~Primer, scales = "free_y",
             labeller = as_labeller(c(
               "fwh"    = "fwhF1/fwhR1",
               "mlCOI"  = "mlCO1intF/jgHCO2198",
               "Uni18S" = "Uni18S/Uni18SR"
             ))) +
  labs(
    title = "Total Zooplankton Hill-Shannon richness across 18 lakes",
    x = " ", y = "Hill-Shannon diversity") +
  scale_fill_manual(values = method_cols) +
  scale_x_discrete(labels = c(
    "freezer" = "Tow net (-20°C)",
    "townet"  = "Tow net (room T)",
    "water"  = "Water sample" )) +
  theme(axis.text.x = element_text(size = 8)) +
  my_theme 
ggsave("Hill-Shannon_diversity_per_method_and_primer_waterpooled.svg",
       width = 8,
       height = 5,
       units = "in")



#Shows the observed ASV richness per collection method and primer to see whether a primer performs best 
alpha_summary <- alpha_all %>%
  dplyr::group_by(Collection_Method, Primer) %>%
  dplyr::summarise(
    mean_richness = mean(Hill_Shannon, na.rm = TRUE),
    se = sd(Hill_Shannon)/sqrt(dplyr::n()),
    .groups = "drop"
  )
ggplot(alpha_summary,
       aes(Collection_Method, mean_richness,
           color = Primer,
           group = Primer)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_richness-se,
                    ymax = mean_richness+se),
                width = 0.1) +
  scale_color_manual(values = colourCAW) +
  my_theme +
  labs(
    y = "Mean zooplankton Hill-Shannon diversity across 18 lakes",
    x = "Collection method"
  )

#Per lake in one big figure
lake_order <- c(
  "Rototekoiti", "Wakatipu", "Ianthe", "Manapouri", "Moeraki","McGregor", "Chalice", "Mapourika", "Sheppard", "Kaweka1",
  "Horseshoe26277", "Poerua", "Kohangapiripiri", "Waihola", "Wiritoa", "34599", "Karangata", "Waihau", "Horseshoe35422", "Kohangatera")
# Convert Lake to factor with specified order
alpha_all$Lake <- factor(alpha_all$Lake, levels = lake_order)
ggplot(alpha_all,
       aes(Collection_Method, Observed,
           color = Primer)) +
  geom_point(position = position_jitter(width = 0.2)) +
  facet_wrap(~Lake) +  # facets will now follow lake_order
  labs(
    y = "Zooplankton observed ASV richness",
    x = "") +
  scale_color_manual(values = colourCAW) +
  my_theme


#Use GLMM with ASVs and Hill_Shannon across Trophic state and surface 
df <- alpha_all
df <- df %>%
  dplyr::filter(!Lake %in% c("Waihau", "Waihola"))

df <- df %>%
  dplyr::mutate(
    Trophic_Group = dplyr::case_when(
      Lake %in% c("Rototekoiti","Wakatipu","Ianthe","Manapouri","Moeraki") ~ "oligotrophic",
      Lake %in% c("McGregor","Chalice","Mapourika","Sheppard","Kaweka1") ~ "mesotrophic",
      Lake %in% c("Horseshoe26277","Poerua","Kohangapiripiri","Wiritoa") ~ "eutrophic",
      Lake %in% c("34599","Karangata","Horseshoe35422","Kohangatera") ~ "supertrophic"
    ) )






##BETA DIVERSITY ASV level

library(phyloseq)
library(vegan)
library(dplyr)

# List of primers
ps_list <- list(
  fwh = ps_fwh_final,
  mlCOI = ps_mlCO1_final,
  Uni18S = ps_Uni18S_final)

# Exclude lakes
exclude_lakes <- c("Waihau", "Waihola")

# NMDS results table
nmds_results <- data.frame(
  Primer = character(),
  Samples = integer(),
  ASVs = integer(),
  Stress = numeric(),
  stringsAsFactors = FALSE
)

# PERMANOVA results table
permanova_results <- data.frame(
  Primer = character(),
  Term = character(),
  F_value = numeric(),
  R2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

bc_list <- list()
nmds_list <- list()

library(pairwiseAdonis)
pairwise_results <- data.frame(
  Primer = character(),
  Comparison = character(),
  F_value = numeric(),
  R2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

#See if water is the driver for differences in NMDS
no_water_results <- data.frame(
  Primer = character(),
  F_value = numeric(),
  R2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)
#Test whether distances btw collection methods are bigger than within a method
within_lake_results <- data.frame(
  Primer = character(),
  mean_same = numeric(),
  mean_diff = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)




# Loop over primers
for(p in names(ps_list)){
  ps <- ps_list[[p]]
  ps <- subset_samples(ps, !Lake %in% exclude_lakes)
  asv <- as.data.frame(otu_table(ps))
  if(taxa_are_rows(ps)) asv <- t(asv)
  # Remove all-zero samples
  row_sums <- rowSums(asv)
  asv <- asv[row_sums > 0, ]
  if(nrow(asv) == 0) next
  meta <- data.frame(sample_data(ps))
  meta <- meta[row_sums > 0, , drop = FALSE]
  
  # Bray-Curtis distance
  bc_dist <- vegdist(asv, method = "bray")
  # Store for Mantel
  bc_list[[p]] <- bc_dist
  
  # Check for dispersion
  disp <- betadisper(bc_dist, meta$Collection_Method)
  print(anova(disp))
  disp_lake <- betadisper(bc_dist, meta$Lake)
  print(anova(disp_lake))
  
  # PERMANOVA
  # PERMANOVA (blocked by Lake)
  set.seed(42)
  perm <- adonis2(
    bc_dist ~ Lake + Collection_Method,
    data = meta,
    permutations = 999,
    by = "terms",)
  
  print(paste("PERMANOVA results for", p))
  print(perm)
  
  perm_df <- as.data.frame(perm)
  perm_df$Term <- rownames(perm_df)
  perm_df <- perm_df[perm_df$Term %in% c("Lake", "Collection_Method"), ]
  
  permanova_results <- rbind(
    permanova_results,
    data.frame(
      Primer = p,
      Term = perm_df$Term,
      F_value = perm_df$F,
      R2 = perm_df$R2,
      p_value = perm_df$`Pr(>F)`
    )
  )
  
  
  # New test btw primers
  dist_mat <- as.matrix(bc_dist)
  dist_df <- as.data.frame(as.table(dist_mat))
  colnames(dist_df) <- c("Sample1", "Sample2", "Distance")
  dist_df <- dist_df[dist_df$Sample1 != dist_df$Sample2, ]
  dist_df <- dist_df[as.character(dist_df$Sample1) < as.character(dist_df$Sample2), ]
  meta$SampleID <- rownames(meta)
  dist_df <- dist_df %>%
    left_join(meta, by = c("Sample1" = "SampleID")) %>%
    rename(Method1 = Collection_Method, Lake1 = Lake) %>%
    left_join(meta, by = c("Sample2" = "SampleID")) %>%
    rename(Method2 = Collection_Method, Lake2 = Lake)
  dist_df <- dist_df %>% filter(Lake1 == Lake2)
  dist_df <- dist_df %>%
    mutate(
      Comparison = ifelse(Method1 == Method2, "same", "different")
    )
  print(table(dist_df$Comparison))
  mean_same <- mean(dist_df$Distance[dist_df$Comparison == "same"])
  mean_diff <- mean(dist_df$Distance[dist_df$Comparison == "different"])
  
  # Statistical test (paired structure approximated via grouping by lake)
  lake_summary <- dist_df %>%
    group_by(Lake1, Comparison) %>%
    summarise(mean_dist = mean(Distance), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = Comparison,
      values_from = mean_dist,
      values_fill = NA
    )
  # Ensure both columns exist
  if(!("same" %in% colnames(lake_summary))) lake_summary$same <- NA
  if(!("different" %in% colnames(lake_summary))) lake_summary$different <- NA
  lake_summary <- lake_summary %>%
    filter(!is.na(same) & !is.na(different))
  if(nrow(lake_summary) > 2){
    test <- wilcox.test(lake_summary$same, lake_summary$different, paired = TRUE)
    
    within_lake_results <- rbind(
      within_lake_results,
      data.frame(
        Primer = p,
        mean_same = mean_same,
        mean_diff = mean_diff,
        p_value = test$p.value
      )
    )
  } else {
    print(paste("Skipping within-lake test for", p, "- insufficient paired data"))
  }
  
  
  # Remove water samples
  subset_idx <- meta$Collection_Method %in% c("freezer", "townet")
  meta_nowater <- meta[subset_idx, ]
  asv_nowater  <- asv[subset_idx, ]
  # Only run if both groups present
  if(length(unique(meta_nowater$Collection_Method)) > 1){
    bc_nowater <- vegdist(asv_nowater, method = "bray")
    perm_nowater <- adonis2(
      bc_nowater ~ Lake + Collection_Method,
      data = meta_nowater,
      permutations = 999,
      by = "terms"
    )
    print(paste("PERMANOVA without water for", p))
    print(perm_nowater)
    cm_row <- as.data.frame(perm_nowater)
    cm_row$Term <- rownames(cm_row)
    cm_row <- cm_row[cm_row$Term == "Collection_Method", ]
    no_water_results <- rbind(
      no_water_results,
      data.frame(
        Primer = p,
        F_value = cm_row$F,
        R2 = cm_row$R2,
        p_value = cm_row$`Pr(>F)` ) )}
  
  
  # NMDS
  nmds <- metaMDS(bc_dist, k = 2, trymax = 100)
  nmds_results <- rbind(
    nmds_results,
    data.frame(
      Primer = p,
      Samples = nrow(asv),
      ASVs = ncol(asv),
      Stress = nmds$stress
    )
  )
  
  nmds <- metaMDS(bc_dist, k = 2, trymax = 100)
  nmds_list[[p]] <- nmds
  
  # NMDS plotting
  nmds_points <- as.data.frame(scores(nmds))
  nmds_points <- cbind(nmds_points, meta)
  nmds_points$Primer <- p
  method_labels <- c(
    "freezer" = "tow net -20°C",
    "townet"  = "tow net room T°",
    "water"  = "water sample")
  primer_labels <- c(
    "fwh"    = "fwhF1/fwhR1",
    "mlCOI"  = "mlCO1intF/jgHCO2198",
    "Uni18S" = "Uni18S/Uni18SR")
  p_plot <- ggplot(
    nmds_points,
    aes(NMDS1, NMDS2, color = Collection_Method, shape = Collection_Method)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = method_cols, labels = method_labels) +
    scale_shape_manual(
      values = c(
        "freezer" = 16,  # triangle
        "townet"  = 17,  # circle
        "water"  = 15  # square
      ),
      labels = method_labels) +
    labs(
      title = paste("NMDS (Bray-Curtis) -", primer_labels[p]),
      color = "Collection Method",shape = "Collection Method") +
    guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
    my_theme
  
  print(p_plot)
  
  ggsave(
    filename = paste0("NMDS_waterpooled_", p, ".svg"),
    plot = p_plot,
    width = 7,
    height = 6
  )
}
#Loop end

nmds_results
permanova_results
no_water_results #PERMANOVA results townet vs freezer
pairwise_results #PERMANOVA results for water vs townet & water vs freezer

# Mantel test between primers
# ---------------------------
primer_pairs <- combn(names(bc_list), 2, simplify = FALSE)

mantel_results <- data.frame(
  Primer1 = character(),
  Primer2 = character(),
  Mantel_r = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for(pair in primer_pairs){
  
  p1 <- pair[1]
  p2 <- pair[2]
  
  # Shared samples
  shared_samples <- intersect(
    attr(bc_list[[p1]], "Labels"),
    attr(bc_list[[p2]], "Labels")
  )
  
  # Skip if too few samples
  if(length(shared_samples) < 3){
    cat("Skipping", p1, "vs", p2, "- not enough shared samples\n")
    next
  }
  
  # Subset matrices
  bc1 <- as.dist(as.matrix(bc_list[[p1]])[shared_samples, shared_samples])
  bc2 <- as.dist(as.matrix(bc_list[[p2]])[shared_samples, shared_samples])
  
  # Mantel test
  mantel_test <- mantel(bc1, bc2, method = "pearson", permutations = 999)
  
  mantel_results <- rbind(
    mantel_results,
    data.frame(
      Primer1 = p1,
      Primer2 = p2,
      Mantel_r = mantel_test$statistic,
      p_value = mantel_test$signif
    )
  )
  
  cat(
    "Mantel test:", p1, "vs", p2,
    "- r =", mantel_test$statistic,
    ", p =", mantel_test$signif, "\n"
  )
}

#PROCRUSTES
primer_pairs <- combn(names(nmds_list), 2, simplify = FALSE)

procrustes_results <- data.frame(
  Primer1 = character(),
  Primer2 = character(),
  Correlation = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for(pair in primer_pairs){
  
  p1 <- pair[1]
  p2 <- pair[2]
  
  coords1 <- scores(nmds_list[[p1]], display = "sites")
  coords2 <- scores(nmds_list[[p2]], display = "sites")
  
  shared_samples <- intersect(rownames(coords1), rownames(coords2))
  
  if(length(shared_samples) < 3){
    cat("Skipping", p1, "vs", p2, "- not enough shared samples\n")
    next
  }
  
  coords1 <- coords1[shared_samples, ]
  coords2 <- coords2[shared_samples, ]
  
  proc <- procrustes(coords1, coords2)
  protest_test <- protest(coords1, coords2, permutations = 999)
  
  procrustes_results <- rbind(
    procrustes_results,
    data.frame(
      Primer1 = p1,
      Primer2 = p2,
      Correlation = protest_test$t0,
      p_value = protest_test$signif
    )
  )
  
  cat(
    "Procrustes:", p1, "vs", p2,
    "- correlation =", protest_test$t0,
    ", p =", protest_test$signif, "\n"
  )
}

plot(proc)

