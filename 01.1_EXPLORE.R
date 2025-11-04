# exploratory analysis based on reccommendations from Zuur et al. (2010)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lubridate); library(stringr)
  library(ggplot2); library(GGally); library(forcats); library(readr)
  library(car)
})

# -------- Step 1: Outliers in Y and X --------
# Y boxplots
p1 <- ggplot(totals_transect, aes(x = Type, y = Total, fill = Type)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
  facet_wrap(~ TransectOrder, nrow = 1) +
  guides(fill = "none") +
  labs(title = "Total fish by Type and Transect", x = "Type", y = "Total fish")
save_plot(p1, "s1_box_total_by_type_transect")

# Cleveland dotplot for totals (order index vs value)
totals_transect <- totals_transect %>% arrange(Total) %>% mutate(row_idx = row_number())
p2 <- ggplot(totals_transect, aes(x = Total, y = row_idx, color = Type)) +
  geom_point(alpha = 0.7) +
  labs(title = "Cleveland dotplot of Total fish", x = "Total fish", y = "Row index")
save_plot(p2, "s1_cleveland_total")

# X boxplots for Depth, Vis, Boats
xcols <- c("Depth","Vis","Boats")
present_x <- intersect(xcols, names(fish_long))
if (length(present_x) > 0) {
  xdat <- fish_long %>%
    distinct(Type, Site, Date, TransectOrder, survey_pair, across(all_of(present_x))) %>%
    pivot_longer(cols = all_of(present_x), names_to = "Xvar", values_to = "Xval")
  p3 <- ggplot(xdat, aes(x = Type, y = Xval, fill = Type)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
    facet_wrap(~ Xvar, scales = "free_y") +
    guides(fill = "none") +
    labs(title = "Covariates by Type", x = "Type", y = "Value")
  save_plot(p3, "s1_x_boxplots_by_type")
}

# -------- Step 2: Homogeneity of variance --------
# Quick residual check using a simple quasi-Poisson (exploration only)
m0 <- glm(Total ~ Type * TransectOrder, data = totals_transect, family = quasipoisson)
res_df <- data.frame(fitted = fitted(m0), resid = residuals(m0, type = "pearson"),
                     Type = totals_transect$Type, TransectOrder = totals_transect$TransectOrder)

p4 <- ggplot(res_df, aes(x = fitted, y = resid, color = Type)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~ TransectOrder) +
  labs(title = "Residuals vs Fitted (quasi-Poisson)", x = "Fitted", y = "Pearson residuals")
save_plot(p4, "s2_resid_vs_fitted")

# -------- Step 3: Normality (of residuals) --------
p5 <- ggplot(res_df, aes(x = resid)) + geom_histogram(bins = 30, fill = "grey70", color = "white") +
  labs(title = "Residual histogram", x = "Pearson residuals", y = "Count")
save_plot(p5, "s3_resid_hist")

png(file.path(output_dir, paste0("s3_resid_qq_", analysis_date, ".png")), width=800, height=600)
qqnorm(res_df$resid); qqline(res_df$resid); dev.off()

# -------- Step 4: Zeros and double zeros --------
# Total-level zeros
zero_tab <- totals_transect %>%
  summarise(
    zeros_total = mean(Total == 0),
    zeros_dived = mean(Total == 0 & Type == "Dived")/mean(Type == "Dived"),
    zeros_undiv = mean(Total == 0 & Type == "Undived")/mean(Type == "Undived")
  )
write_csv(zero_tab, file.path(output_dir, paste0("s4_zero_props_", analysis_date, ".csv")))

# Species zero proportions
species_zero <- fish_long %>%
  group_by(Species) %>%
  summarise(zero_prop = mean(Count == 0), .groups = "drop") %>%
  arrange(desc(zero_prop))
write_csv(species_zero, file.path(output_dir, paste0("s4_species_zero_props_", analysis_date, ".csv")))

p6 <- ggplot(species_zero, aes(x = reorder(Species, zero_prop), y = zero_prop)) +
  geom_col() + coord_flip() + ylim(0,1) +
  labs(title = "Zero proportions by species", x = "Species", y = "Proportion zero")
save_plot(p6, "s4_species_zero_bar", w=7, h=9)

# Double-zero heatmap for top species set
top_spp <- species_zero %>% filter(zero_prop < 0.95) %>% slice_min(zero_prop, n = 15) %>% pull(Species)
wide_counts <- fish_long %>%
  filter(Species %in% top_spp) %>%
  group_by(survey_pair, Species) %>%
  summarise(n = sum(Count > 0), .groups = "drop") %>%
  mutate(pres = as.integer(n > 0)) %>%
  select(-n) %>%
  pivot_wider(names_from = Species, values_from = pres, values_fill = 0)

if (nrow(wide_counts) > 0) {
  mat <- as.matrix(wide_counts[,-1, drop=FALSE])
  # pairwise double zeros proportion
  spp <- colnames(mat)
  dz <- matrix(NA_real_, ncol = length(spp), nrow = length(spp), dimnames = list(spp, spp))
  for (i in seq_along(spp)) for (j in seq_along(spp)) {
    both_zero <- mean(mat[,i] == 0 & mat[,j] == 0)
    dz[i,j] <- both_zero
  }
  dz_long <- as.data.frame(as.table(dz))
  names(dz_long) <- c("Spp1","Spp2","prop_dzero")
  p7 <- ggplot(dz_long, aes(Spp1, Spp2, fill = prop_dzero)) +
    geom_tile() + scale_fill_viridis_c(limits = c(0,1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Double-zero proportion heatmap (top species)", x = "", y = "", fill = "Prop DZ")
  save_plot(p7, "s4_double_zero_heatmap", w=8, h=7)
}

# -------- Step 5: Collinearity among covariates --------
covars <- fish_long %>%
  distinct(survey_pair, Depth, Vis, Boats) %>%
  select(Depth, Vis, Boats) %>%
  mutate(across(everything(), as.numeric))
covars <- covars[, colSums(!is.na(covars)) > 0, drop = FALSE]

if (ncol(covars) >= 2) {
  # VIF from a dummy linear model
  suppressWarnings({
    lm_dummy <- lm(scale(Depth, center=TRUE, scale=TRUE) ~ ., data = as.data.frame(scale(covars, center=TRUE, scale=TRUE)))
  })
  vifs <- tryCatch(car::vif(lm_dummy), error = function(e) NA)
  if (!all(is.na(vifs))) {
    vif_tab <- tibble(variable = names(vifs), VIF = as.numeric(vifs))
    write_csv(vif_tab, file.path(output_dir, paste0("s5_vif_", analysis_date, ".csv")))
  }
  # pairs plot
  png(file.path(output_dir, paste0("s5_pairs_covariates_", analysis_date, ".png")), width=900, height=900, res=120)
  print(GGally::ggpairs(as.data.frame(covars)))
  dev.off()
}

# -------- Step 6: Relationships Y vs X --------
if (length(present_x) > 0) {
  yx <- totals_transect %>%
    left_join(fish_long %>% distinct(survey_pair, across(all_of(present_x))), by = "survey_pair") %>%
    pivot_longer(cols = all_of(present_x), names_to = "Xvar", values_to = "Xval")
  p8 <- ggplot(yx, aes(x = Xval, y = Total, color = Type)) +
    geom_point(alpha = 0.6) +
    geom_smooth(se = FALSE, method = "loess") +
    facet_wrap(~ Xvar, scales = "free_x") +
    labs(title = "Total vs covariates", x = "Covariate", y = "Total")
  save_plot(p8, "s6_total_vs_covariates", w=9, h=6)
}

# -------- Step 7: Interactions (Type x TransectOrder) --------
# Spaghetti by pair, colored by Type
p9 <- ggplot(totals_transect, aes(x = TransectOrder, y = Total, group = survey_pair, color = Type)) +
  geom_line(alpha = 0.35) + geom_point(alpha = 0.75) +
  labs(title = "Paired A→B change by Type", x = "Transect order", y = "Total")
save_plot(p9, "s7_pair_lines")

# Simple slopes by Type
p10 <- ggplot(totals_transect, aes(x = as.numeric(TransectOrder), y = Total, color = Type)) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.05, height = 0)) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_continuous(breaks = c(1,2), labels = c("A","B")) +
  labs(title = "Slopes A→B by Type", x = "Transect order", y = "Total")
save_plot(p10, "s7_slopes_by_type")

# Δ = B − A per pair
delta_tab <- totals_transect %>%
  select(Type, survey_pair, TransectOrder, Total) %>%
  pivot_wider(names_from = TransectOrder, values_from = Total) %>%
  mutate(delta_B_minus_A = B - A)
write_csv(delta_tab, file.path(output_dir, paste0("s7_delta_by_pair_", analysis_date, ".csv")))

p11 <- ggplot(delta_tab, aes(x = Type, y = delta_B_minus_A, fill = Type)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 2) +
  guides(fill = "none") +
  labs(title = "Delta (B − A) by Type", x = "Type", y = "B − A")
save_plot(p11, "s7_delta_box")

# -------- Step 8: Independence --------
# Time series per site to see within-site dependence
p12 <- ggplot(totals_transect, aes(x = Date, y = Total, color = TransectOrder)) +
  geom_line(alpha = 0.5) + geom_point(alpha = 0.7) +
  facet_wrap(~ Site, scales = "free_y") +
  labs(title = "Within-site totals over time", x = "Date", y = "Total")
save_plot(p12, "s8_time_by_site", w=10, h=8)

# Lag-1 correlation within site (exploratory)
lag1 <- totals_transect %>%
  arrange(Site, Date, TransectOrder) %>%
  group_by(Site) %>%
  summarise(lag1_cor = suppressWarnings(cor(Total[-n()], Total[-1], use = "complete.obs")),
            n = n(), .groups = "drop")
write_csv(lag1, file.path(output_dir, paste0("s8_lag1_by_site_", analysis_date, ".csv")))

# -------- Done --------
message("✓ Exploration complete. Files saved to: ", normalizePath(output_dir))


#### basic summaries ##### 
# ---- totals & by-site ----
stopifnot(exists("fish_long"), exists("output_dir"))

# Ensure transect-level totals exist
if (!exists("totals_transect")) {
  totals_transect <- fish_long %>%
    group_by(Type, Site, Date, TransectOrder, survey_pair) %>%
    summarise(Total = sum(Count, na.rm = TRUE), .groups = "drop")
}

# Overall counts
summary_overall <- tibble::tibble(
  n_transects      = nrow(totals_transect),
  n_pairs          = dplyr::n_distinct(totals_transect$survey_pair),
  n_sites          = dplyr::n_distinct(totals_transect$Site),
  n_survey_days    = dplyr::n_distinct(totals_transect$Date),
  first_date       = min(totals_transect$Date, na.rm = TRUE),
  last_date        = max(totals_transect$Date, na.rm = TRUE),
  total_fish_count = sum(totals_transect$Total, na.rm = TRUE),
  mean_per_transect= mean(totals_transect$Total, na.rm = TRUE),
  median_per_transect = median(totals_transect$Total, na.rm = TRUE)
)

readr::write_csv(summary_overall, file.path(output_dir, "summary_overall.csv"))

# By site
by_site <- totals_transect %>%
  group_by(Site) %>%
  summarise(
    n_transects   = dplyr::n(),
    n_pairs       = dplyr::n_distinct(survey_pair),
    n_days        = dplyr::n_distinct(Date),
    first_date    = min(Date, na.rm = TRUE),
    last_date     = max(Date, na.rm = TRUE),
    days_span     = as.integer(last_date - first_date),
    mean_total    = mean(Total, na.rm = TRUE),
    median_total  = median(Total, na.rm = TRUE),
    sd_total      = sd(Total, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(by_site, file.path(output_dir, "summary_by_site.csv"))

# By site and Type
by_site_type <- totals_transect %>%
  group_by(Site, Type) %>%
  summarise(
    n_transects = dplyr::n(),
    n_pairs     = dplyr::n_distinct(survey_pair),
    mean_total  = mean(Total, na.rm = TRUE),
    median_total= median(Total, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(by_site_type, file.path(output_dir, "summary_by_site_type.csv"))

# ---- fish counted & density ----
# total fish counted (from long)
total_fish_long <- fish_long %>%
  summarise(total_fish_count = sum(Count, na.rm = TRUE))

readr::write_csv(total_fish_long, file.path(output_dir, "total_fish_from_long.csv"))

# Optional density if area is available
# Expect a column named Area_m2 at transect level (one value per row in totals_transect)
density_available <- "Area_m2" %in% names(totals_transect)
if (density_available) {
  dens_tab <- totals_transect %>%
    mutate(density_m2 = Total / Area_m2) %>%
    group_by(Site) %>%
    summarise(
      mean_density_m2   = mean(density_m2, na.rm = TRUE),
      median_density_m2 = median(density_m2, na.rm = TRUE),
      .groups = "drop"
    )
  readr::write_csv(dens_tab, file.path(output_dir, "density_by_site.csv"))
}

# ---- survey effort ----
present_x <- intersect(c("Depth","Vis","Boats","Duration"), names(fish_long))

effort_base <- totals_transect %>%
  # one row is one transect
  summarise(
    n_transects = dplyr::n(),
    n_pairs     = dplyr::n_distinct(survey_pair),
    n_sites     = dplyr::n_distinct(Site),
    n_days      = dplyr::n_distinct(Date),
    first_date  = min(Date, na.rm = TRUE),
    last_date   = max(Date, na.rm = TRUE),
    .groups = "drop"
  )

# effort by Type and TransectOrder
effort_type_transect <- totals_transect %>%
  count(Type, TransectOrder, name = "n_transects")

readr::write_csv(effort_base,          file.path(output_dir, "effort_overall.csv"))
readr::write_csv(effort_type_transect, file.path(output_dir, "effort_by_type_transect.csv"))

# covariate summaries at survey-pair level to avoid double counting A and B
effort_covars <- NULL
if (length(present_x) > 0) {
  effort_covars <- fish_long %>%
    distinct(survey_pair, Site, Date, across(all_of(present_x))) %>%
    summarise(
      across(all_of(present_x),
             list(min = ~min(.x, na.rm = TRUE),
                  q25 = ~stats::quantile(.x, 0.25, na.rm = TRUE),
                  median = ~stats::median(.x, na.rm = TRUE),
                  mean = ~mean(.x, na.rm = TRUE),
                  q75 = ~stats::quantile(.x, 0.75, na.rm = TRUE),
                  max = ~max(.x, na.rm = TRUE)),
             .names = "{.col}_{.fn}")
    )
  readr::write_csv(effort_covars, file.path(output_dir, "effort_covariates_summary.csv"))
}


# quick console prints
print(summary_overall)
print(by_site)
print(effort_type_transect)
if (density_available) message("Density by site saved (Area_m2 detected).")

fish_long %>%
  distinct(Site) %>%
  arrange(Site) %>%
  print(n = Inf)
