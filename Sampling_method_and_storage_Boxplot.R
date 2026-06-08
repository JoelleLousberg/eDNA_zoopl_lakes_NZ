rm(list= ls())

getwd()

setwd("C:/Users/joelle.lousberg/OneDrive - Cawthron/Documents/Chapter_2/R_step3/Boxplot")

df  <- read_csv("All_Rare_zooplspecies_permethod_andLake Backup.csv")


method_cols <- c(
  "morphology" = "#76c044",
  "townet"     = "#7f1424",
  "freezer"    = "#73FFFF",
  "water"      = "#FF944C"
)

# Columns with species
primers <- c("fwh", "mlCO1", "Uni18S", "morphology")

# 1) Split species by ';;' and remove duplicates
for (primer in primers) {
  df[[primer]] <- str_split(df[[primer]], pattern = "\\s*;;\\s*|\\s*;\\s*|\\s*,\\s*") %>%
    map(~ .x[.x != "" & !is.na(.x)] %>% unique())
}

# 2) Count unique species per primer
for (primer in primers) {
  df[[paste0(primer, "_count")]] <- map_int(df[[primer]], length)
}

# 3) Count total unique species across all primers per row
df$total_unique_species <- pmap_int(df[primers], ~ length(unique(c(...))))

# 4) View result
result <- df %>%
  select(Lake, Collection_Method, ends_with("_count"), total_unique_species)

print(result)

output_file <- "lake_species_counts.csv"

# Write the result to CSV
write.csv(result, file = output_file, row.names = FALSE)

# Confirmation message (optional)
cat("CSV file saved as", output_file, "\n")


#BOXPLOT 
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggpubr)

df <- read.csv("All_Rare_zooplspecies_permethod_andlake_watersummed.csv", stringsAsFactors = FALSE)
df <- df %>%
  filter(!Lake %in% c("Waihola", "Waihau"))
df <- df %>%
  mutate(Lake = factor(Lake), Collection_Method = factor(Collection_Method))

# Order methods by median richness
df$Collection_Method <- forcats::fct_reorder(
  df$Collection_Method, df$total_species, .fun = median)
y_max <- max(df$total_species, na.rm = TRUE)


#MIXED-EFFECTS MODEL
mod <- lmer(total_species ~ Collection_Method + (1 | Lake),data = df)
anova(mod)

# ESTIMATED MARGINAL MEANS (MODEL RICHNESS)
emm <- emmeans(mod, ~ Collection_Method)
emm_df <- as.data.frame(emm)

#PAIRWISE COMPARISONS + EFFECT SIZES
pairs_emm <- pairs(emm, adjust = "BH")
pairs_df <- as.data.frame(pairs_emm)


pairs_df <- pairs_df %>%
  separate(contrast,into = c("group1", "group2"),sep = " - ")
pairs_df <- pairs_df %>%
  mutate(
    significance = case_when(
      p.value <= 0.001 ~ "***",
      p.value <= 0.01  ~ "**",
      p.value <= 0.05  ~ "*",
      TRUE             ~ "ns"))
    #effect_label = paste0(significance, "\nΔ = ", round(estimate, 1)))
pairs_df$y.position <- y_max + 2 + seq_len(nrow(pairs_df))
pairs_df

#SAMPLE SIZE LABELS
n_df <- df %>%
  dplyr::group_by(Collection_Method) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop")

p <- ggplot(df, aes(Collection_Method, total_species)) +
  # Boxplots
  geom_boxplot(width = 0.6, fill = "grey90", color = "grey70", outlier.shape = NA) +
  # Raw observations
  geom_jitter(width = 0.15, alpha = 0.35, size = 2) +
  # Model estimated means
  geom_point(data = emm_df, aes(y = emmean), color = "darkred", size = 4) +
  # Confidence intervals
  geom_errorbar(data = emm_df,
                aes(ymin = lower.CL,
                    ymax = upper.CL,
                    y = emmean),
                color = "darkred",
                width = 0.15,
                linewidth = 1) +
  # Effect sizes + stars
  stat_pvalue_manual(pairs_df, label = "significance", tip.length = 0.01, size = 4) +
  #labels
  scale_x_discrete(labels = c(
    "morphology" = "Morphology\n(Tow net room T)",
    "townet"     = "eDNA\n(Tow net room T)",
    "freezer"    = "eDNA\n(Tow net -20°C)",
    "water"      = "eDNA\n(Water sample)"
  ))+
  
  labs(x = " ",
    y = "Species richness") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 30, hjust = 1, size = 12),
    axis.title.y = element_text(size = 16))
p
ggsave("Boxplot_morph_eDNA.svg", plot = p, width = 8, height = 6, units = "in")


## Boxplot per trophic category. Dropping Heaton (no eDNA), Waihau (no morph) and Waihola (bad sample)
## =====================================
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggpubr)

#table has no Heaton or Waihau
df <- read.csv("All_Rare_zooplspecies_permethod_andlake_watersummed.csv", stringsAsFactors = FALSE)

df <- df %>%
  filter(!Lake %in% c("Waihola", "Waihau"))

## TROPHIC GROUPS
df <- df %>%
  mutate(
    Trophic_Group = case_when(
      Lake %in% c("Rototekoiti","Wakatipu","Ianthe","Manapouri","Moeraki")
      ~ "oligotrophic",
      Lake %in% c("McGregor","Chalice","Mapourika","Sheppard","Kaweka1")
      ~ "mesotrophic",
      Lake %in% c("Horseshoe26277","Poerua","Kohangapiripiri","Wiritoa")
      ~ "eutrophic",
      Lake %in% c("34599","Karangata","Horseshoe35422","Kohangatera")
      ~ "supertrophic") )

df <- df %>%
  mutate(
    Lake = factor(Lake),
    Collection_Method = factor(Collection_Method),
    Trophic_Group = factor(Trophic_Group, levels = c("oligotrophic","mesotrophic","eutrophic","supertrophic")
    ) )

## INTERACTION MIXED MODEL
mod <- lmer(total_species ~ Collection_Method * Trophic_Group + (1 | Lake), data = df)
anova(mod)

## MODEL PREDICTIONS
emm <- emmeans(mod, ~ Collection_Method | Trophic_Group)
emm_df <- as.data.frame(emm)

## SAMPLE SIZES
n_df <- df %>%
  group_by(Trophic_Group, Collection_Method) %>%
  summarise(n = n(), .groups = "drop")
p <- ggplot(df, aes(Collection_Method,total_species)) +
  ## raw distributions
  geom_boxplot(width = 0.6, fill = "grey92", color = "grey60", outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 2) +
  ## MODEL PREDICTIONS
  geom_point(data = emm_df, aes(y = emmean), color = "darkred", size = 4) +
  geom_errorbar(data = emm_df, aes(ymin = lower.CL, ymax = upper.CL, y = emmean),
                color = "darkred", width = 0.15, linewidth = 1) +
  ## facet by trophic state
  facet_wrap(~ Trophic_Group, nrow = 1) +
  ## sample size labels
  geom_text(data = n_df,aes(y = min(df$total_species) - 1, label = paste0("n=", n)), size = 3) +
  
  labs(
    title = "Collection-method performance across trophic states",
    subtitle = "Red points = mixed-model predictions (±95% CI)",
    x = "Collection Method",
    y = "Total Species Richness") +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold"))
p

#Plot the relative efficience for each method per trophic state relative to the best performing method
emm_eff <- emmeans(mod,~ Collection_Method | Trophic_Group)
eff_df <- as.data.frame(emm_eff)
eff_df <- eff_df %>%
  group_by(Trophic_Group) %>%
  mutate(
    rel_efficiency = 100 * emmean / max(emmean),
    rel_lower = 100 * lower.CL / max(emmean),
    rel_upper = 100 * upper.CL / max(emmean)) %>%
  ungroup()

eff_plot <- ggplot(
  eff_df,
  aes(x = Collection_Method, y = rel_efficiency, colour = Collection_Method )) +
  
  geom_point(size = 4) +
  geom_errorbar(
    aes(ymin = rel_lower, ymax = rel_upper), width = 0.15,linewidth = 0.7) +
  
  facet_wrap(~ Trophic_Group, nrow = 1) +
  scale_colour_manual(values = method_cols) +
  labs(
    title = "Relative sampling efficiency across trophic states",
    subtitle = "Efficiency scaled to best-performing method within each trophic group",
    x = "Collection Method",
    y = "Relative Efficiency (%)") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold")
  )
eff_plot


##BOXPLOT along lake surface area
#table has no Heaton or Waihau
df <- read.csv("All_Rare_zooplspecies_permethod_andlake_watersummed.csv", stringsAsFactors = FALSE)

df <- df %>%
  filter(!Lake %in% c("Waihola", "Waihau"))

## TROPHIC GROUPS
df <- df %>%
  mutate(
    Trophic_Group = case_when(
      Lake %in% c("Kaweka1","Horseshoe26277","Kohangapiripiri","34599","Karangata")
      ~ "<20m^2",
      Lake %in% c("McGregor","Kohangatera","Wiritoa","Rototekoiti","Horseshoe35422")
      ~ "20-40m^2",
      Lake %in% c("Poerua","Sheppard","Chalice","Moeraki")
      ~ "40-400m^2",
      Lake %in% c("Wakatipu","Manapouri","Ianthe","Mapourika")
      ~ ">400m^2") )

df <- df %>%
  mutate(
    Lake = factor(Lake),
    Collection_Method = factor(Collection_Method),
    Trophic_Group = factor(Trophic_Group, levels = c("<20m^2","20-40m^2","40-400m^2",">400m^2")
    ) )

## INTERACTION MIXED MODEL
mod <- lmer(total_species ~ Collection_Method * Trophic_Group + (1 | Lake), data = df)
anova(mod)

## MODEL PREDICTIONS
emm <- emmeans(mod, ~ Collection_Method | Trophic_Group)
emm_df <- as.data.frame(emm)

## SAMPLE SIZES
n_df <- df %>%
  group_by(Trophic_Group, Collection_Method) %>%
  summarise(n = n(), .groups = "drop")
p <- ggplot(df, aes(Collection_Method,total_species)) +
  ## raw distributions
  geom_boxplot(width = 0.6, fill = "grey92", color = "grey60", outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 2) +
  ## MODEL PREDICTIONS
  geom_point(data = emm_df, aes(y = emmean), color = "darkred", size = 4) +
  geom_errorbar(data = emm_df, aes(ymin = lower.CL, ymax = upper.CL, y = emmean),
                color = "darkred", width = 0.15, linewidth = 1) +
  ## facet by trophic state
  facet_wrap(~ Trophic_Group, nrow = 1) +
  ## sample size labels
  geom_text(data = n_df,aes(y = min(df$total_species) - 1, label = paste0("n=", n)), size = 3) +
  
  labs(
    title = "Collection-method performance surface area dependent",
    subtitle = "Red points = mixed-model predictions (±95% CI)",
    x = "Collection Method",
    y = "Total Species Richness") +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold"))
p

#Plot the relative efficience for each method per trophic state relative to the best performing method
emm_eff <- emmeans(mod,~ Collection_Method | Trophic_Group)
eff_df <- as.data.frame(emm_eff)
eff_df <- eff_df %>%
  group_by(Trophic_Group) %>%
  mutate(
    rel_efficiency = 100 * emmean / max(emmean),
    rel_lower = 100 * lower.CL / max(emmean),
    rel_upper = 100 * upper.CL / max(emmean)) %>%
  ungroup()

eff_plot <- ggplot(
  eff_df,
  aes(x = Collection_Method, y = rel_efficiency, colour = Collection_Method)) +
  
  geom_point(size = 4) +
  geom_errorbar(
    aes(ymin = rel_lower, ymax = rel_upper), width = 0.15,linewidth = 0.7) +
  
  facet_wrap(~ Trophic_Group, nrow = 1) +
  scale_colour_manual(values = method_cols) +
  labs(
    title = "Relative sampling efficiency surface area dependent",
    subtitle = "Efficiency scaled to best-performing method within each group",
    x = "Collection Method",
    y = "Relative Efficiency (%)") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold")
  )
eff_plot














library(dplyr)
library(dplyr)
library(readr)
library(stringr)
library(purrr)

# bulletproof pipeline: explicit key ensures Lake + Collection_Method are always present
library(dplyr)
library(readr)
library(stringr)
library(purrr)

# --- read ---
fwh    <- read_csv("fwh_Rare_zooplspecies_relabund_and_reads_permethod_andlake.csv")
mlco1  <- read_csv("mlCO1_Rare_zooplspecies_relabund_and_reads_permethod_andlake.csv")
uni18s <- read_csv("Uni18S_Rare_zooplspecies_relabund_and_reads_permethod_andlake.csv")

# --- cleaning ---
clean_primer <- function(df) {
  df %>%
    select(any_of(c("Lake", "Collection_Method_Sum", "Species"))) %>%
    mutate(
      Lake = as.character(Lake) %>% str_trim(),
      Collection_Method_Sum = as.character(Collection_Method_Sum) %>% str_trim(),
      Species = as.character(Species) %>% str_trim() %>% str_remove("^s__")
    ) %>%
    filter(!is.na(Lake) & !is.na(Collection_Method_Sum)) %>%
    filter(!str_detect(tolower(Species), "unclassified"))
}

fwh_c    <- clean_primer(fwh)
mlco1_c  <- clean_primer(mlco1)
uni18s_c <- clean_primer(uni18s)

# --- robust summary: build explicit key, aggregate, then restore keys ---
make_summary <- function(df, prefix) {
  df %>%
    distinct(Lake, Collection_Method_Sum, Species) %>%
    group_by(Lake, Collection_Method_Sum) %>%
    summarise(
      species_list = list(sort(unique(Species))),
      species_n    = n_distinct(Species),
      .groups = "drop"
    ) %>%
    # Use rename(), not :=
    rename(
      !!paste0(prefix, "_species_list") := species_list,
      !!paste0(prefix, "_species_n")    := species_n
    )
}

fwh_sum    <- make_summary(fwh_c,   "fwh")
mlco1_sum  <- make_summary(mlco1_c, "mlco1")
uni18s_sum <- make_summary(uni18s_c,"uni18s")

# quick check (should show Lake and Collection_Method_Sum plus species columns)
print(head(fwh_sum, 3))
print(head(mlco1_sum, 3))
print(head(uni18s_sum, 3))

# --- merge safely ---
merged <- list(fwh_sum, mlco1_sum, uni18s_sum) %>%
  reduce(full_join, by = c("Lake", "Collection_Method_Sum"))

# --- normalize missing list columns and counts ---
merged <- merged %>%
  mutate(
    fwh_species_list    = map(fwh_species_list, ~ if (is.null(.x) || all(is.na(.x))) character() else .x),
    mlco1_species_list  = map(mlco1_species_list, ~ if (is.null(.x) || all(is.na(.x))) character() else .x),
    uni18s_species_list = map(uni18s_species_list, ~ if (is.null(.x) || all(is.na(.x))) character() else .x),
    fwh_species_n    = replace_na(fwh_species_n, 0),
    mlco1_species_n  = replace_na(mlco1_species_n, 0),
    uni18s_species_n = replace_na(uni18s_species_n, 0)
  )

merged <- merged %>%
       rowwise() %>%
       mutate(total_species = length(unique(c(unlist(fwh),unlist(mlco1),unlist(uni18s_species_list) ) ) ) ) %>%
       ungroup()

# save and preview
saveRDS(merged, "All_Rare_zooplspecies_permethod_andlake_onlyeDNA.rds")
print(head(merged, 10))

merged_export <- merged %>%
  mutate(
    fwh    = sapply(fwh,    paste, collapse = "; "),
    mlCO1  = sapply(mlco1,  paste, collapse = "; "),
    Uni18S = sapply(uni18s, paste, collapse = "; "))
write.csv(merged, "All_Rare_zooplspecies_permethod_andLake_onlyeDNA_watersummed.csv", row.names = FALSE)
