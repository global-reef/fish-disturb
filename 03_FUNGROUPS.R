##### Load libraries #####
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(glmmTMB)
  library(DHARMa)
  library(emmeans)
  library(ggplot2)
  library(readr)
  library(broom.mixed)
  library(purrr)
  library(stringr)
})

##### Inputs & assumptions #####
# Assumes you already have: fish_long, output_dir
stopifnot(exists("fish_long"), exists("output_dir"))
dir.create(file.path(output_dir, "functional_groups"), showWarnings = FALSE, recursive = TRUE)
fg_dir <- file.path(output_dir, "functional_groups")

##### 1. Aggregate to functional-group totals per transect #####
# Ensure expected factor ordering (matches earlier scripts)
fish_long <- fish_long %>%
  mutate(
    Functional_Group = factor(Functional_Group,
                              levels = c("Grazer", "Invertivore", "Mesopredator", "HTLP"),
                              ordered = TRUE
    ),
    Type = factor(Type, levels = c("Dived", "Undived")),            # Dived as baseline
    TransectOrder = factor(TransectOrder, levels = c("A", "B"))
  )

# Pair-level Boats (single row per survey_pair); use max as safe default
boats_pair <- fish_long %>%
  distinct(survey_pair, TransectOrder, Boats) %>%
  group_by(survey_pair) %>%
  summarise(Boats = suppressWarnings(max(Boats, na.rm = TRUE)), .groups = "drop") %>%
  mutate(Boats = ifelse(is.infinite(Boats), NA, Boats))

# Group totals
totals_group <- fish_long %>%
  group_by(Type, Site, Date, TransectOrder, survey_pair, Functional_Group) %>%
  summarise(GroupTotal = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  left_join(boats_pair, by = "survey_pair")

##### 2. Helper for emmeans CI column names #####
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

##### 3. Per-group modelling function #####
fit_group_model <- function(g) {
  message(">>> Fitting functional group: ", g)
  df <- totals_group %>% filter(Functional_Group == g)
  
  # Boats present in data?
  has_boats_data <- "Boats" %in% names(df) && any(!is.na(df$Boats))
  df_use <- if (has_boats_data) dplyr::filter(df, !is.na(Boats)) else df
  
  # Base formulas
  f_base  <- GroupTotal ~ Type * TransectOrder + (1|Site) + (1|survey_pair)
  f_boats <- GroupTotal ~ Type * TransectOrder + scale(Boats) + (1|Site) + (1|survey_pair)
  
  # Fit candidates
  m_nb    <- glmmTMB(f_base,  family = nbinom2, data = df_use)
  m_nb_zi <- glmmTMB(f_base,  family = nbinom2, ziformula = ~1, data = df_use)
  
  if (has_boats_data) {
    m_nb_boats    <- glmmTMB(f_boats, family = nbinom2, data = df_use)
    m_nb_boats_zi <- glmmTMB(f_boats, family = nbinom2, ziformula = ~1, data = df_use)
    aic_tab <- AIC(m_nb, m_nb_zi, m_nb_boats, m_nb_boats_zi)
  } else {
    aic_tab <- AIC(m_nb, m_nb_zi)
  }
  
  aic_tab <- as.data.frame(aic_tab) |> tibble::rownames_to_column("model") |> dplyr::arrange(AIC)
  readr::write_csv(aic_tab, file.path(fg_dir, paste0("AIC_", g, ".csv")))
  
  final_name  <- aic_tab$model[1]
  final_model <- get(final_name, inherits = TRUE)
  
  # ---- Diagnostics ----
  set.seed(42)
  sim <- DHARMa::simulateResiduals(final_model, n = 2000)
  png(file.path(fg_dir, paste0("DHARMa_", g, ".png")), 900, 700); plot(sim); dev.off()
  disp <- DHARMa::testDispersion(sim)
  zi   <- tryCatch(DHARMa::testZeroInflation(sim), error = function(e) NULL)
  
  # ---- emmeans ----
  # Does the chosen model actually include Boats?
  formula_text <- paste(deparse(formula(final_model)), collapse = " ")
  has_boats_in_model <- grepl("Boats", formula_text)
  
  if (has_boats_in_model) {
    med_boats <- median(df_use$Boats, na.rm = TRUE)
    emm <- emmeans::emmeans(final_model, ~ Type * TransectOrder,
                            at = list(Boats = med_boats), type = "response")
  } else {
    emm <- emmeans::emmeans(final_model, ~ Type * TransectOrder, type = "response")
  }
  
  # Slopes & diff-in-diff
  emm_slopes <- emmeans::contrast(emm, method = "revpairwise", by = "Type")
  levs <- with(as.data.frame(emm), interaction(Type, TransectOrder, sep = ":"))
  w <- setNames(rep(0, length(levs)), levs)
  w["Undived:B"] <-  1; w["Undived:A"] <- -1; w["Dived:B"] <- -1; w["Dived:A"] <-  1
  emm_inter <- emmeans::contrast(emm, list(diff_in_diff = w))
  
  # Exports
  readr::write_csv(as.data.frame(summary(emm)),          file.path(fg_dir, paste0("emm_", g, ".csv")))
  readr::write_csv(as.data.frame(summary(emm_slopes)),   file.path(fg_dir, paste0("emm_slopes_", g, ".csv")))
  readr::write_csv(as.data.frame(summary(emm_inter)),    file.path(fg_dir, paste0("emm_diff_in_diff_", g, ".csv")))
  
  # Plot
  emm_df <- as.data.frame(summary(emm)) %>% normalize_emm_cis()
  p <- ggplot(emm_df, aes(TransectOrder, response, color = Type, group = Type)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.05) +
    labs(title = paste0(g, " — Estimated means"), x = "Transect order", y = "Expected count") +
    theme_clean
  ggsave(file.path(fg_dir, paste0("plot_emm_", g, ".png")), p, width = 7, height = 5, dpi = 300)
  
  # Coefs
  coefs <- broom.mixed::tidy(final_model, effects = "fixed", conf.int = TRUE)
  readr::write_csv(coefs, file.path(fg_dir, paste0("coefs_", g, ".csv")))
  
  inter_row <- dplyr::filter(coefs, term == "TypeUndived:TransectOrderB")
  tibble::tibble(
    Functional_Group = g,
    final_model = final_name,
    AIC = aic_tab$AIC[1],
    interaction_est = inter_row$estimate,
    interaction_SE  = inter_row$std.error,
    interaction_z   = inter_row$statistic,
    interaction_p   = inter_row$p.value
  )
}


##### 4. Fit models for all four groups #####
groups_to_run <- levels(fish_long$Functional_Group)
groups_to_run <- groups_to_run[groups_to_run %in% c("Grazer","Invertivore","Mesopredator","HTLP")]

fg_summary <- map_dfr(groups_to_run, fit_group_model) %>%
  arrange(AIC)

write_csv(fg_summary, file.path(fg_dir, "functional_groups_summary.csv"))
print(fg_summary)

message("✓ Functional-group modelling complete. Outputs in: ", normalizePath(fg_dir))


##### 5. Combined forest plot of interaction across functional groups #####

# Option A (quick): use fg_summary (SE-based Wald 95% CI)
fg_forest <- fg_summary %>%
  mutate(
    ci_low  = interaction_est - 1.96 * interaction_SE,
    ci_high = interaction_est + 1.96 * interaction_SE
  ) %>%
  # order by estimate (or choose a fixed order you like)
  mutate(Functional_Group = factor(Functional_Group,
                                   levels = fg_summary$Functional_Group[order(interaction_est)]))

readr::write_csv(fg_forest, file.path(fg_dir, "functional_groups_interaction_forest.csv"))

p_forest <- ggplot(fg_forest,
                   aes(y = Functional_Group, x = interaction_est)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.15) +
  geom_point(size = 2) +
  labs(
    title = "Interaction (Undived × B) by functional group",
    x = "Log-count difference in A→B change (Undived vs Dived)",
    y = NULL
  ) +
  theme_clean
p_forest
ggsave(file.path(fg_dir, "plot_interaction_forest.png"),
       p_forest, width = 7.5, height = 4.8, dpi = 300)


##### 6. Optional: p-value adjustment summary (for transparency) #####
fg_forest_padj <- fg_summary %>%
  mutate(
    p_adj_fdr = p.adjust(interaction_p, method = "fdr"),
    p_adj_bonf = p.adjust(interaction_p, method = "bonferroni")
  )
readr::write_csv(fg_forest_padj, file.path(fg_dir, "functional_groups_summary_padj.csv"))
print(fg_forest_padj)


##### Combined emmeans panel across functional groups #####
# Requires: files 'emm_<group>.csv' in fg_dir created by the per-group models

stopifnot(exists("fg_dir"))
groups <- c("Grazer","Invertivore","Mesopredator","HTLP")

# Read, normalize, and combine
emm_all <- purrr::map_dfr(groups, function(g) {
  fp <- file.path(fg_dir, paste0("emm_", g, ".csv"))
  df <- readr::read_csv(fp, show_col_types = FALSE)
  df <- normalize_emm_cis(df)
  df$Functional_Group <- g
  df %>%
    dplyr::mutate(
      Type = factor(Type, levels = c("Dived","Undived")),
      TransectOrder = factor(TransectOrder, levels = c("A","B")),
      Functional_Group = factor(Functional_Group,
                                levels = c("Grazer","Invertivore","Mesopredator","HTLP"))
    )
})

# Faceted emmeans plot (shared y-axis for comparability; set scales="free_y" if needed)
p_emm_panel <- ggplot(emm_all,
                      aes(x = TransectOrder, y = response, color = Type, group = Type)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.05) +
  facet_wrap(~ Functional_Group, nrow = 1) +
  labs(
    title = "Estimated means (emmeans) by functional group",
    x = "Transect order",
    y = "Expected count"
  ) +
  theme_clean +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(fg_dir, "plot_emm_functional_groups_panel.png"),
       p_emm_panel, width = 12, height = 3.8, dpi = 300)

# If counts differ wildly across groups, use free y-scale:
p_emm_panel_free <- p_emm_panel + facet_wrap(~ Functional_Group, nrow = 2, scales = "free_y")
ggsave(file.path(fg_dir, "plot_emm_functional_groups_panel_freeY.png"), p_emm_panel_free, width = 12, height = 3.8, dpi = 300)

