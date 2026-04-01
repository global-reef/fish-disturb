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

### building nice spp plots (called-out species only) #####
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(purrr)
  library(patchwork)
})

### expects ####
# theme_clean, reef_cols
# spp_dir contains: species_models_summary.csv and emm_<Species>.csv
# spp_lookup contains: Species, sci_name

normalize_emm_cis <- function(df) {
  # standardise mean column
  if (!("response" %in% names(df))) {
    if ("emmean" %in% names(df)) df <- dplyr::rename(df, response = emmean)
    else if ("Estimate" %in% names(df)) df <- dplyr::rename(df, response = Estimate)
  }
  # standardise CI columns (handle either CL or HPD naming)
  if (!("lower.CL" %in% names(df))) {
    if (all(c("asymp.LCL","asymp.UCL") %in% names(df))) {
      df <- dplyr::rename(df, lower.CL = asymp.LCL, upper.CL = asymp.UCL)
    } else if (all(c("lower.HPD","upper.HPD") %in% names(df))) {
      df <- dplyr::rename(df, lower.CL = lower.HPD, upper.CL = upper.HPD)
    }
  }
  df
}

### 0) species set (as in Results text) ####
spp_called <- c(
  "Angelfish",      # Pomacanthus spp.
  "Thicklip",       # Hemigymnus melapterus
  "Red_Breast",     # Cheilinus fasciatus
  "Rabbitfish",     # Siganus spp.
  "Parrotfish",     # Scarus spp.
  "Butterflyfish"   # Chaetodon spp.
)

### 1) read species summary for sig flag + estimates ####
spp_sum <- readr::read_csv(file.path(spp_dir, "species_models_summary.csv"),
                           show_col_types = FALSE) %>%
  mutate(
    Species = as.character(Species),
    sig = interaction_p < 0.05
  ) %>%
  filter(Species %in% spp_called) %>%
  mutate(Species = factor(Species, levels = spp_called))  # lock order

stopifnot(nrow(spp_sum) == length(spp_called))

### 2) map to scientific labels ####
map_spp <- spp_lookup %>%
  transmute(Species = as.character(Species),
            label   = as.character(sci_name)) %>%
  distinct()

spp_sum <- spp_sum %>%
  left_join(map_spp, by = "Species") %>%
  mutate(label = ifelse(is.na(label) | label == "", as.character(Species), label))

### 3) read ONLY the emmeans tables for these species ####
emm_files <- file.path(spp_dir, paste0("emm_", as.character(spp_sum$Species), ".csv"))
miss <- emm_files[!file.exists(emm_files)]
if (length(miss) > 0) {
  stop("Missing emm files for called-out species:\n", paste(basename(miss), collapse = "\n"))
}

emm_all <- purrr::map_dfr(emm_files, function(fp) {
  sp <- stringr::str_match(basename(fp), "^emm_(.*)\\.csv$")[, 2]
  readr::read_csv(fp, show_col_types = FALSE) %>%
    normalize_emm_cis() %>%
    mutate(Species = sp)
}) %>%
  filter(Species %in% spp_called) %>%   # hard guarantee
  left_join(spp_sum %>% transmute(Species = as.character(Species), label, interaction_est, sig),
            by = "Species") %>%
  mutate(
    response  = as.numeric(response),
    lower.CL  = as.numeric(lower.CL),
    upper.CL  = as.numeric(upper.CL),
    Type = factor(Type, levels = c("Dived", "Undived")),
    TransectOrder = factor(TransectOrder, levels = c("A", "B")),
    x = as.numeric(TransectOrder),
    Species = factor(Species, levels = spp_called),
    label = factor(label, levels = spp_sum$label)
  )

### 4) summarise to one row per Species x Type x TransectOrder (should be 24 rows) ####
emm_all_u <- emm_all %>%
  group_by(Species, label, Type, TransectOrder) %>%
  summarise(
    response = mean(response, na.rm = TRUE),
    lower.CL = mean(lower.CL, na.rm = TRUE),
    upper.CL = mean(upper.CL, na.rm = TRUE),
    x        = first(x),
    sig      = first(sig),
    .groups  = "drop"
  )

# sanity checks
stopifnot(nrow(emm_all_u) == 24)
stopifnot(all(levels(droplevels(emm_all_u$Species)) %in% spp_called))

# ensure every Species x Type has both A and B
chk_AB <- emm_all_u %>%
  count(Species, Type, TransectOrder) %>%
  tidyr::pivot_wider(names_from = TransectOrder, values_from = n, values_fill = 0)

if (any(chk_AB$A == 0 | chk_AB$B == 0)) {
  stop("Some Species x Type are missing A or B rows:\n",
       paste0(capture.output(print(chk_AB %>% filter(A == 0 | B == 0))), collapse = "\n"))
}

### 5) percent-change labels within each Species x Type (should be 12 rows) ####
delta_lab <- emm_all_u %>%
  select(Species, label, Type, TransectOrder, response) %>%   # <-- NO x here
  pivot_wider(
    id_cols    = c(Species, label, Type),
    names_from = TransectOrder,
    values_from = response,
    values_fn  = mean
  ) %>%
  mutate(
    A = as.numeric(A),
    B = as.numeric(B)
  ) %>%
  filter(is.finite(A), is.finite(B), A > 0) %>%
  mutate(
    pct = 100 * (B - A) / A,
    lab = paste0(ifelse(pct >= 0, "+", ""), sprintf("%.0f%%", pct))
  ) %>%
  left_join(
    emm_all_u %>%
      group_by(Species, label, Type) %>%
      summarise(
        x_mid = mean(x),
        y_mid = mean(response),
        .groups = "drop"
      ),
    by = c("Species", "label", "Type")
  ) %>%
  left_join(
    emm_all_u %>%
      group_by(Species, label) %>%
      summarise(y_max = max(response, na.rm = TRUE), .groups = "drop"),
    by = c("Species", "label")
  ) %>%
  mutate(y_pos = y_mid + 0.08 * y_max)

stopifnot(nrow(delta_lab) == 12)
delta_lab %>% count(Species, Type)

### 6) one ggplot per species + patchwork ####
make_one <- function(sp) {
  
  d  <- emm_all_u %>% filter(Species == sp)
  dl <- delta_lab %>% filter(Species == sp)
  
  # hard guard: if we somehow have no rows, return NULL (will be filtered out)
  if (nrow(d) == 0) return(NULL)
  
  title_txt <- unique(as.character(d$label))
  title_txt <- title_txt[!is.na(title_txt)][1]
  
  # use subtitle for the star so the title stays “pure italics”
  sub_txt <- if (any(d$sig, na.rm = TRUE)) "*" else NULL
  
  p <- ggplot(d, aes(x = x, y = response, color = Type, group = Type)) +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = Type),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    geom_text(
      data = dl,
      aes(x = x_mid, y = y_pos, label = lab, color = Type),
      size = 3.2,
      fontface = "bold",
      show.legend = FALSE
    ) +
    scale_x_continuous(breaks = c(1, 2), labels = c("A", "B")) +
    scale_color_manual(values = reef_cols) +
    scale_fill_manual(values  = reef_cols, guide = "none") +
    labs(title = title_txt, subtitle = sub_txt, x = NULL, y = NULL) +
    theme_clean +
    theme(
      plot.title = element_text(face = "italic", size = 12, hjust = 0.5),
      plot.subtitle = element_text(face = "plain", size = 12, hjust = 0.5),
      legend.position = "top",
      legend.title = element_blank()
    )
  
  return(p)
}

plots <- purrr::map(spp_called, make_one) |>   # use spp_called to preserve your intended order
  purrr::keep(~ inherits(.x, "ggplot"))

stopifnot(length(plots) == 6)  # should be exactly your 6 species

p_panel <- patchwork::wrap_plots(plots, ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

p_panel <- p_panel +
  patchwork::plot_annotation(
    caption = "* indicates Type × Transect interaction p < 0.05",
    theme = theme(plot.caption = element_text(hjust = 0, size = 9))
  ) &
  labs(x = "Transect order", y = "Expected abundance")

ggsave(
  file.path(spp_dir, "fig_emm_species_panel_calledout_clean.png"),
  p_panel, width = 10, height = 6.0, dpi = 300, bg = "white"
)

p_panel



library(dplyr)
library(readr)
library(purrr)

# get all species coefficient files
coef_files <- list.files(spp_dir, pattern = "^coefs_.*\\.csv$", full.names = TRUE)

s5_table <- map_dfr(coef_files, function(fp) {
  
  sp <- gsub("coefs_|\\.csv", "", basename(fp))
  
  read_csv(fp, show_col_types = FALSE) %>%
    filter(effect == "fixed") %>%
    transmute(
      Species = sp,
      Term = term,
      Estimate = estimate,
      SE = std.error,
      z = statistic,
      p = p.value,
      CI_lower = conf.low,
      CI_upper = conf.high
    )
})
s5_table <- s5_table %>%
  mutate(
    Term = recode(Term,
                  "(Intercept)" = "Intercept",
                  "TypeUndived" = "Type (Undived vs Dived)",
                  "TransectOrderB" = "Transect order (B vs A)",
                  "TypeUndived:TransectOrderB" = "Type × Transect order",
                  "scale(Boats)" = "Boats"
    )
  )

s5_table <- s5_table %>%
  filter(Term == "Type × Transect order") %>%
  mutate(
    Direction = case_when(
      Estimate > 0 ~ "Undived decline weaker",
      Estimate < 0 ~ "Undived decline stronger",
      TRUE ~ "No difference"
    )
  ) %>%
  select(
    Species,
    Estimate,
    SE,
    z,
    p,
    CI_lower,
    CI_upper,
    Direction
  ) %>%
  arrange(desc(Estimate))
print(s5_table, n=Inf)


##### S6: Model seleciton tables ##### 

library(dplyr)
library(readr)
library(purrr)
library(stringr)

### TOTAL MODEL ####
aic_total <- AIC(m_nb, m_nb_zi, m_nb_boats, m_nb_boats_zi) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Model") %>%
  mutate(Level = "Total", Group = "All")

### FUNCTIONAL GROUPS ####
fg_files <- list.files(file.path(output_dir, "functional_groups"),
                       pattern = "^AIC_.*\\.csv$", full.names = TRUE)

aic_fg <- map_dfr(fg_files, function(fp) {
  g <- str_remove(basename(fp), "AIC_|\\.csv")
  read_csv(fp, show_col_types = FALSE) %>%
    mutate(Level = "Functional group", Group = g)
})

### SPECIES ####
spp_files <- list.files(file.path(output_dir, "species_models"),
                        pattern = "^AIC_.*\\.csv$", full.names = TRUE)

aic_spp <- map_dfr(spp_files, function(fp) {
  sp <- str_remove(basename(fp), "AIC_|\\.csv")
  read_csv(fp, show_col_types = FALSE) %>%
    mutate(Level = "Species", Group = sp)
})

### COMBINE ####
table_s6 <- bind_rows(aic_total, aic_fg, aic_spp) %>%
  arrange(Level, Group, AIC)

table_s6
