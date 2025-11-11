#### Load libraries ####
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(vegan)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(stringr)
})

#### Inputs & setup ####
stopifnot(exists("fish_long"), exists("output_dir"))
dir.create(file.path(output_dir, "community_composition"), showWarnings = FALSE, recursive = TRUE)
comm_dir <- file.path(output_dir, "community_composition")

fish_long <- fish_long %>%
  mutate(
    Type = factor(Type, levels = c("Dived","Undived")),
    Functional_Group = factor(Functional_Group,
                              levels = c("Grazer","Invertivore","Mesopredator","HTLP"),
                              ordered = TRUE)
  )

##### 1. Functional-group composition between site types #####
fg_comm <- fish_long %>%
  group_by(Site, Type, Functional_Group) %>%
  summarise(Abundance = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Functional_Group, values_from = Abundance, values_fill = 0)

# Hellinger transform
fg_matrix <- fg_comm %>%
  select(Grazer:HTLP) %>%
  decostand(method = "hellinger")

# Bray–Curtis PERMANOVA
set.seed(42)
fg_dist <- vegdist(fg_matrix, method = "bray")
fg_adonis <- adonis2(fg_dist ~ Type, data = fg_comm, permutations = 999)
write.csv(as.data.frame(fg_adonis), file.path(comm_dir, "fg_permanova.csv"))

# NMDS
fg_nmds <- metaMDS(fg_matrix, k = 2, trymax = 200)
saveRDS(fg_nmds, file.path(comm_dir, "fg_nmds.rds"))
fg_scores <- as.data.frame(scores(fg_nmds, display = "sites")) %>%
  bind_cols(fg_comm %>% select(Site, Type))

p_fg <- ggplot(fg_scores, aes(NMDS1, NMDS2, color = Type)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Type), type = "t", linetype = 2) +
  theme_minimal() +
  labs(title = "NMDS of functional-group composition",
       subtitle = "Bray–Curtis dissimilarity, Hellinger-transformed",
       color = "Site type")
p_fg
ggsave(file.path(comm_dir, "plot_fg_nmds.png"), p_fg, width = 7, height = 5, dpi = 300)

# SIMPER: which functional groups drive the difference
fg_sim <- simper(fg_matrix, group = fg_comm$Type)
fg_sim_df <- as.data.frame(summary(fg_sim)$Dived_Undived) %>%
  tibble::rownames_to_column("Functional_Group") %>%
  arrange(desc(average))
write.csv(fg_sim_df, file.path(comm_dir, "fg_simper.csv"), row.names = FALSE)

##### 2. Species-level composition between site types #####
spp_comm <- fish_long %>%
  group_by(Site, Type, Species) %>%
  summarise(Abundance = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0)

# Hellinger transform
spp_matrix <- spp_comm %>%
  select(-Site, -Type) %>%
  decostand(method = "hellinger")

set.seed(42)
spp_dist <- vegdist(spp_matrix, method = "bray")
spp_adonis <- adonis2(spp_dist ~ Type, data = spp_comm, permutations = 999)
write.csv(as.data.frame(spp_adonis), file.path(comm_dir, "spp_permanova.csv"))

# NMDS
spp_nmds <- metaMDS(spp_matrix, k = 2, trymax = 200)
saveRDS(spp_nmds, file.path(comm_dir, "spp_nmds.rds"))
spp_scores <- as.data.frame(scores(spp_nmds, display = "sites")) %>%
  bind_cols(spp_comm %>% select(Site, Type))

p_spp <- ggplot(spp_scores, aes(NMDS1, NMDS2, color = Type)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Type), type = "t", linetype = 2) +
  theme_minimal() +
  labs(title = "NMDS of species composition",
       subtitle = "Bray–Curtis dissimilarity, Hellinger-transformed",
       color = "Site type")
p_spp
ggsave(file.path(comm_dir, "plot_spp_nmds.png"), p_spp, width = 7, height = 5, dpi = 300)

# SIMPER: which species drive the difference #####
spp_sim <- simper(spp_matrix, group = spp_comm$Type)
spp_sim_df <- as.data.frame(summary(spp_sim)$Dived_Undived) %>%
  tibble::rownames_to_column("Species") %>%
  arrange(desc(average))
write.csv(spp_sim_df, file.path(comm_dir, "spp_simper.csv"), row.names = FALSE)


# Read SIMPER table and merge scientific names
spp_sim_df <- readr::read_csv(file.path(comm_dir, "spp_simper.csv"), show_col_types = FALSE) %>%
  left_join(map_spp, by = "Species") %>%
  mutate(
    sci_name = ifelse(is.na(sci_name), Species, sci_name),
    sci_name = factor(sci_name)
  )

# Identify which group has higher mean abundance
spp_mean_by_type <- fish_long %>%
  group_by(Type, Species) %>%
  summarise(mean_abund = mean(Count, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Type, values_from = mean_abund, values_fill = 0) %>%
  mutate(Direction = case_when(Dived > Undived ~ "Dived", TRUE ~ "Undived")) %>%
  select(Species, Direction)

# Merge with SIMPER data
sim_df <- spp_sim_df %>%
  left_join(spp_mean_by_type, by = "Species") %>%
  mutate(
    Significant = ifelse(p <= 0.05, "p ≤ 0.05", "n.s."),
    cum_perc = cumsum(average / sum(average) * 100)
  ) %>%
  filter(cum_perc <= 70) %>%
  mutate(Species_label = forcats::fct_reorder(sci_name, average))

# Reef palette for Dived vs Undived 
reef_cols <- c("Dived" = "#66BFA6", "Undived" = "#007A87")

# Plot
p_simper <- ggplot(sim_df,
                   aes(x = average * 100, y = Species_label,
                       fill = Direction, alpha = Significant)) +
  geom_col(color = "black", linewidth = 0.2) +
  scale_fill_manual(values = reef_cols) +
  scale_alpha_manual(values = c("p ≤ 0.05" = 1, "n.s." = 0.5)) +
  scale_y_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) +
  labs(
    title = "Top species driving reef-type dissimilarity",
    subtitle = "Bray–Curtis SIMPER (Dived vs Undived)",
    x = "Average contribution to dissimilarity (%)",
    y = NULL,
    fill = "Higher on",
    alpha = "Significance"
  ) +
  theme_clean +
  coord_cartesian(xlim = c(0, max(sim_df$average * 100) * 1.25))

ggsave(file.path(comm_dir, "plot_spp_simper_sig.png"),
       p_simper, width = 7.5, height = 6, dpi = 300)


p_simper





##### 3. Ordination stress summary #####
ord_summary <- tibble(
  Level = c("Functional group", "Species"),
  Stress = c(fg_nmds$stress, spp_nmds$stress)
)
write_csv(ord_summary, file.path(comm_dir, "nmds_stress_summary.csv"))

message("✓ Community composition analysis complete. Outputs in: ", normalizePath(comm_dir))
