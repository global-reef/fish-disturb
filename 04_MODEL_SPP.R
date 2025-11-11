#### Load libraries ####
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(glmmTMB); library(DHARMa)
  library(emmeans); library(ggplot2); library(readr); library(broom.mixed)
  library(purrr); library(stringr); library(tibble); library(patchwork)
})

#### Inputs ####
stopifnot(exists("fish_long"), exists("output_dir"))
spp_dir <- file.path(output_dir, "species_models")
dir.create(spp_dir, showWarnings = FALSE, recursive = TRUE)

#### 1) Species totals per transect ####
# fish_long expected columns: Species, Count, Type, Site, Date, TransectOrder, survey_pair, Boats
fish_long <- fish_long %>%
  mutate(
    Type = factor(Type, levels = c("Dived","Undived")),
    TransectOrder = factor(TransectOrder, levels = c("A","B"))
  )

boats_pair <- fish_long %>%
  distinct(survey_pair, TransectOrder, Boats) %>%
  group_by(survey_pair) %>%
  summarise(Boats = suppressWarnings(max(Boats, na.rm = TRUE)), .groups = "drop") %>%
  mutate(Boats = ifelse(is.infinite(Boats), NA, Boats))

totals_spp <- fish_long %>%
  group_by(Type, Site, Date, TransectOrder, survey_pair, Species) %>%
  summarise(SppTotal = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  left_join(boats_pair, by = "survey_pair")

#### 2) Helper for emmeans CI names ####
normalize_emm_cis <- function(df) {
  if (!("response" %in% names(df))) {
    if ("emmean" %in% names(df)) df <- dplyr::rename(df, response = emmean)
    else if ("Estimate" %in% names(df)) df <- dplyr::rename(df, response = Estimate)
    else if ("rate" %in% names(df)) df <- dplyr::rename(df, response = rate)
    else if ("prob" %in% names(df)) df <- dplyr::rename(df, response = prob)
  }
  if (!("lower.CL" %in% names(df))) {
    if ("asymp.LCL" %in% names(df)) df <- dplyr::rename(df, lower.CL = asymp.LCL)
    else if ("lower.HPD" %in% names(df)) df <- dplyr::rename(df, lower.CL = lower.HPD)
  }
  if (!("upper.CL" %in% names(df))) {
    if ("asymp.UCL" %in% names(df)) df <- dplyr::rename(df, upper.CL = asymp.UCL)
    else if ("upper.HPD" %in% names(df)) df <- dplyr::rename(df, upper.CL = upper.HPD)
  }
  df
}

#### 3) Per-species model ####
fit_species_model <- function(sp) {
  message(">>> Fitting species: ", sp)
  df <- totals_spp %>% filter(Species == sp)
  
  has_boats <- "Boats" %in% names(df) && any(!is.na(df$Boats))
  df_use <- if (has_boats) filter(df, !is.na(Boats)) else df
  
  f_base  <- SppTotal ~ Type * TransectOrder + (1|Site) + (1|survey_pair)
  f_boats <- SppTotal ~ Type * TransectOrder + scale(Boats) + (1|Site) + (1|survey_pair)
  
  m_nb    <- glmmTMB(f_base,  family = nbinom2, data = df_use)
  m_nb_zi <- glmmTMB(f_base,  family = nbinom2, ziformula = ~1, data = df_use)
  
  if (has_boats) {
    m_nb_boats    <- glmmTMB(f_boats, family = nbinom2, data = df_use)
    m_nb_boats_zi <- glmmTMB(f_boats, family = nbinom2, ziformula = ~1, data = df_use)
    aic_tab <- AIC(m_nb, m_nb_zi, m_nb_boats, m_nb_boats_zi)
  } else {
    aic_tab <- AIC(m_nb, m_nb_zi)
  }
  
  aic_tab <- as.data.frame(aic_tab) |> tibble::rownames_to_column("model") |> arrange(AIC)
  write_csv(aic_tab, file.path(spp_dir, paste0("AIC_", sp, ".csv")))
  final_name  <- aic_tab$model[1]
  final_model <- get(final_name, inherits = TRUE)
  
  set.seed(42)
  sim <- DHARMa::simulateResiduals(final_model, n = 2000)
  png(file.path(spp_dir, paste0("DHARMa_", sp, ".png")), 900, 700); plot(sim); dev.off()
  invisible(try(DHARMa::testDispersion(sim), silent = TRUE))
  invisible(try(DHARMa::testZeroInflation(sim), silent = TRUE))
  
  has_boats_in_model <- grepl("Boats", paste(deparse(formula(final_model)), collapse = " "))
  if (has_boats_in_model) {
    emm <- emmeans::emmeans(final_model, ~ Type * TransectOrder,
                            at = list(Boats = median(df_use$Boats, na.rm = TRUE)),
                            type = "response")
  } else {
    emm <- emmeans::emmeans(final_model, ~ Type * TransectOrder, type = "response")
  }
  
  emm_slopes <- emmeans::contrast(emm, method = "revpairwise", by = "Type")
  levs <- with(as.data.frame(emm), interaction(Type, TransectOrder, sep = ":"))
  w <- setNames(rep(0, length(levs)), levs); w["Undived:B"] <- 1; w["Undived:A"] <- -1; w["Dived:B"] <- -1; w["Dived:A"] <- 1
  emm_inter <- emmeans::contrast(emm, list(diff_in_diff = w))
  
  write_csv(as.data.frame(summary(emm)),        file.path(spp_dir, paste0("emm_", sp, ".csv")))
  write_csv(as.data.frame(summary(emm_slopes)), file.path(spp_dir, paste0("emm_slopes_", sp, ".csv")))
  write_csv(as.data.frame(summary(emm_inter)),  file.path(spp_dir, paste0("emm_diff_in_diff_", sp, ".csv")))
  
  emm_df <- as.data.frame(summary(emm)) %>% normalize_emm_cis()
  p <- ggplot(emm_df, aes(TransectOrder, response, color = Type, group = Type)) +
    geom_line(size = 1) + geom_point(size = 2) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.05) +
    labs(title = paste0(sp, ": estimated means"),
         x = "Transect order", y = "Expected count") +
    theme_clean
  ggsave(file.path(spp_dir, paste0("plot_emm_", sp, ".png")), p, width = 7, height = 5, dpi = 300)
  
  coefs <- broom.mixed::tidy(final_model, effects = "fixed", conf.int = TRUE)
  write_csv(coefs, file.path(spp_dir, paste0("coefs_", sp, ".csv")))
  
  inter <- filter(coefs, term == "TypeUndived:TransectOrderB")
  tibble(
    Species = sp, final_model = final_name, AIC = aic_tab$AIC[1],
    interaction_est = inter$estimate, interaction_SE = inter$std.error,
    interaction_z = inter$statistic, interaction_p = inter$p.value
  )
}

#### 4) Choose species and fit ####
spp_keep <- totals_spp %>%
  group_by(Species) %>%
  summarise(n_nonzero = sum(SppTotal > 0), .groups = "drop") %>%
  filter(n_nonzero >= 10) %>%
  pull(Species)

spp_summary <- map_dfr(spp_keep, fit_species_model) %>% arrange(AIC)
write_csv(spp_summary, file.path(spp_dir, "species_models_summary.csv"))
print(spp_summary)
message("✓ Species modelling complete. Outputs in: ", normalizePath(spp_dir))

#### 5) Species lookup and scientific names ####
spp_lookup <- tibble::tribble(
  ~Species,            ~Functional_Group, ~Genus,                ~Species_epithet, ~sci_name,
  "Parrotfish",        "Grazer",          "Scarus",              "spp.",           "Scarus spp.",
  "Rabbitfish",        "Grazer",          "Siganus",             "spp.",           "Siganus spp.",
  "Butterflyfish",     "Grazer",          "Chaetodon",           "spp.",           "Chaetodon spp.",
  "Angelfish",         "Invertivore",     "Pomacanthus",         "spp.",           "Pomacanthus spp.",
  "Cleaner_Wrasse",    "Invertivore",     "Labroides",           "dimidiatus",     "Labroides dimidiatus",
  "Batfish",           "Invertivore",     "Ephippidae",          "spp.",           "Ephippidae spp.",
  "Thicklip",          "Invertivore",     "Hemigymnus",          "melapterus",     "Hemigymnus melapterus",
  "Red_Breast",        "Invertivore",     "Cheilinus",           "fasciatus",      "Cheilinus fasciatus",
  "Slingjaw",          "Invertivore",     "Epibulus",            "insidiator",     "Epibulus insidiator",
  "Sweetlips",         "Invertivore",     "Diagramma/Plectorhinchus","spp.",        "Diagramma/ Plectorhinchus spp.",
  "Squirrel.Soldier",  "Invertivore",     "Holocentridae",       "spp.",           "Holocentridae spp.",
  "Triggerfish",       "Invertivore",     "Balistidae",          "spp.",           "Balistidae spp.",
  "Porcupine.Puffer",  "Invertivore",     "Diodon/Tetraodon",    "spp.",           "Diodon/ Tetraodon spp.",
  "Ray",               "Mesopredator",    "Taeniura/Neotrygon",  "spp.",           "Taeniura/ Neotrygon spp.",
  "sml_Snapper",       "Mesopredator",    "Lutjanus",            "spp.",           "Lutjanus (<30cm) spp.",
  "lrg_Snapper",       "HTLP",            "Lutjanus",            "spp.",           "Lutjanus (>30cm) spp.",
  "Eel",               "Mesopredator",    "Gymnothorax",         "spp.",           "Gymnothorax spp.",
  "Trevally",          "HTLP",            "Caranx",              "spp.",           "Caranx spp.",
  "Emperorfish",       "Mesopredator",    "Lethrinus",           "spp.",           "Lethrinus spp.",
  "sml_Grouper",       "Mesopredator",    "Cephalopholis/Epinephelus","spp.",       "Cephalopholis/ Epinephelus spp.",
  "lrg_Grouper",       "HTLP",            "Epinephelus",         "spp.",           "Epinephelus (>30cm)/ Plectropomus spp.",
  "Barracuda",         "HTLP",            "Sphyraena",           "spp.",           "Sphyraena spp."
) %>%
  mutate(Species = as.character(Species), sci_name = as.character(sci_name))

map_spp <- spp_lookup %>% select(Species, sci_name) %>% distinct()

# Ensure types and column names align
spp_summary <- spp_summary %>% mutate(Species = as.character(Species))
if (!"sci_name" %in% names(spp_summary)) spp_summary <- mutate(spp_summary, sci_name = NA_character_)
spp_summary <- spp_summary %>%
  select(-sci_name) %>%
  left_join(map_spp, by = "Species") %>%
  mutate(label = ifelse(is.na(sci_name), Species, sci_name))

#### 6) Forest of interaction effects (scientific names) ####
spp_forest <- spp_summary %>%
  mutate(
    ci_low  = interaction_est - 1.96 * interaction_SE,
    ci_high = interaction_est + 1.96 * interaction_SE
  ) %>%
  arrange(interaction_est) %>%
  mutate(label = factor(label, levels = label))

write_csv(spp_forest, file.path(spp_dir, "species_interaction_forest.csv"))

p_forest <- ggplot(spp_forest, aes(y = label, x = interaction_est)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.12) +
  geom_point(size = 1.8) +
  labs(
    title = "Interaction Undived × B by species",
    x = "Log-count difference in A to B change (Undived vs Dived)", y = NULL
  ) +
  theme_clean
ggsave(file.path(spp_dir, "plot_species_interaction_forest_sciname.png"),
       p_forest, width = 7.5, height = 10, dpi = 300)

#### 7) Emmeans facets and heatmap (scientific names) ####
emm_all <- map_dfr(spp_summary$Species, function(sp) {
  fp <- file.path(spp_dir, paste0("emm_", sp, ".csv"))
  if (!file.exists(fp)) return(NULL)
  read_csv(fp, show_col_types = FALSE) %>%
    normalize_emm_cis() %>%
    mutate(Species = sp)
}) %>%
  left_join(map_spp, by = "Species") %>%
  mutate(
    label = ifelse(is.na(sci_name), Species, sci_name),
    Type = factor(Type, levels = c("Dived","Undived")),
    TransectOrder = factor(TransectOrder, levels = c("A","B"))
  )

p_emm <- ggplot(emm_all, aes(TransectOrder, response, color = Type, group = Type)) +
  geom_line(size = 0.6) + geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.05) +
  facet_wrap(~ label, scales = "free_y", ncol = 5) +
  labs(y = "Expected count", x = "Transect order", title = "Species specific responses") +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(face = "bold"))
ggsave(file.path(spp_dir, "plot_emm_species_facets_sciname.png"),
       p_emm, width = 12, height = 9, dpi = 300)

emm_heat <- emm_all %>%
  unite(cond, Type, TransectOrder, sep = "_") %>%
  select(label, cond, response) %>%
  pivot_wider(names_from = cond, values_from = response) %>%
  pivot_longer(-label, names_to = "Condition", values_to = "response") %>%
  mutate(label = factor(label))

p_heat <- ggplot(emm_heat, aes(Condition, label, fill = response)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  labs(x = "Condition", y = NULL, fill = "Predicted count",
       title = "Species level expected counts") +
  theme_minimal(base_size = 9)
ggsave(file.path(spp_dir, "plot_emm_species_heatmap_sciname.png"),
       p_heat, width = 6.5, height = 10, dpi = 300)

#### 8) Ranked effects and functional composition (uses lookup FG) ####
spp_rank <- spp_summary %>%
  mutate(
    ci_low  = interaction_est - 1.96 * interaction_SE,
    ci_high = interaction_est + 1.96 * interaction_SE,
    abs_est = abs(interaction_est),
    sig = interaction_p < 0.05
  ) %>%
  arrange(desc(abs_est)) %>%
  mutate(label = factor(label, levels = rev(label)))

p_rank <- ggplot(spp_rank, aes(x = interaction_est, y = label, color = sig)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("gray40", "firebrick")) +
  labs(x = "Interaction effect Undived × B", y = NULL,
       title = "Ranked species level interaction effects", color = "p < 0.05") +
  theme_minimal(base_size = 10)
ggsave(file.path(spp_dir, "plot_species_ranked_effects_sciname.png"),
       p_rank, width = 7.5, height = 10, dpi = 300)

spp_effects_fg <- spp_summary %>%
  select(Species, interaction_est) %>%
  left_join(spp_lookup %>% select(Species, Functional_Group), by = "Species") %>%
  mutate(direction = ifelse(interaction_est > 0, "Positive", "Negative"))

p_fg_bar <- ggplot(spp_effects_fg, aes(Functional_Group, fill = direction)) +
  geom_bar(position = "fill") +
  labs(y = "Proportion of species", x = NULL,
       title = "Directional species responses by functional group") +
  theme_minimal(base_size = 10)
ggsave(file.path(spp_dir, "plot_species_direction_by_functional_group.png"),
       p_fg_bar, width = 6.5, height = 4.5, dpi = 300)

#### 9) Quick summary panel ####
suppressWarnings({
  (p_forest + p_fg_bar) + plot_layout(ncol = 2)
})
