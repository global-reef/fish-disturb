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
})

##### 1. Prepare totals dataset #####
if (!exists("totals_transect")) {
  totals_transect <- fish_long %>%
    group_by(Type, Site, Date, TransectOrder, survey_pair) %>%
    summarise(Total = sum(Count, na.rm = TRUE), .groups = "drop")
}

totals_transect <- totals_transect %>%
  mutate(
    Type = factor(Type, levels = c("Dived", "Undived")),
    TransectOrder = factor(TransectOrder, levels = c("A", "B"))
  )

##### 2. Fit baseline models #####
# Poisson vs NB
m_pois <- glmmTMB(Total ~ Type * TransectOrder + (1|Site) + (1|survey_pair),
                  family = poisson, data = totals_transect)

m_nb <- glmmTMB(Total ~ Type * TransectOrder + (1|Site) + (1|survey_pair),
                family = nbinom2, data = totals_transect)

AIC(m_nb, m_pois)

##### 3. Extend models (zero inflation and Boats) #####
m_nb_zi <- glmmTMB(Total ~ Type * TransectOrder + (1|Site) + (1|survey_pair),
                   ziformula = ~1, family = nbinom2, data = totals_transect)

m_nb_boats <- glmmTMB(Total ~ Type * TransectOrder + scale(Boats) + (1|Site) + (1|survey_pair),
                      family = nbinom2, data = totals_transect)

m_nb_boats_zi <- glmmTMB(Total ~ Type * TransectOrder + scale(Boats) + (1|Site) + (1|survey_pair),
                         ziformula = ~1, family = nbinom2, data = totals_transect)

AIC(m_nb, m_nb_zi, m_nb_boats, m_nb_boats_zi)

##### 4. Select final model and run diagnostics #####
final_model <- m_nb_boats_zi
set.seed(42)
sim_zi <- simulateResiduals(final_model, n = 2000)
png(file.path(output_dir, "final_DHARMa_residuals.png"), 900, 700)
plot(sim_zi)
dev.off()
disp_test <- testDispersion(sim_zi)
zi_test <- testZeroInflation(sim_zi)

capture.output({
  cat("Model: m_nb_boats_zi\n\n")
  print(summary(final_model))
  cat("\nDHARMa dispersion test:\n"); print(disp_test)
  cat("\nDHARMa zero inflation test:\n"); print(zi_test)
}, file = file.path(output_dir, "final_model_summary.txt"))

##### 5. Estimated means and contrasts (median Boats) #####
med_boats <- median(totals_transect$Boats, na.rm = TRUE)

emm <- emmeans(final_model,
               ~ Type * TransectOrder,
               at = list(Boats = med_boats),
               type = "response")

# Within-Type slopes (B - A)
emm_slopes <- contrast(emm, method = "revpairwise", by = "Type")

# Diff-in-diff: (Undived:B - Undived:A) - (Dived:B - Dived:A)
levs <- with(as.data.frame(emm), interaction(Type, TransectOrder, sep = ":"))
coeff <- setNames(rep(0, length(levs)), levs)
coeff["Undived:B"] <-  1; coeff["Undived:A"] <- -1
coeff["Dived:B"]   <- -1; coeff["Dived:A"]   <-  1
emm_interaction <- contrast(emm, list(diff_in_diff = coeff))

##### 6. Export summary tables #####
write_csv(as.data.frame(summary(emm)),
          file.path(output_dir, "emm_TypeXTransect_medianBoats.csv"))
write_csv(as.data.frame(summary(emm_slopes)),
          file.path(output_dir, "emm_withinType_slopes.csv"))
write_csv(as.data.frame(summary(emm_interaction)),
          file.path(output_dir, "emm_diff_in_diff.csv"))

##### 7. Helper for CI columns #####
normalize_emm_cis <- function(df) {
  if (!("lower.CL" %in% names(df)) && "asymp.LCL" %in% names(df)) {
    df <- dplyr::rename(df, lower.CL = asymp.LCL, upper.CL = asymp.UCL)
  }
  if (!("response" %in% names(df)) && "emmean" %in% names(df)) {
    df <- dplyr::rename(df, response = emmean)
  }
  df
}

##### 8. Plot: Type × Transect (median Boats) #####
emm_df <- as.data.frame(summary(emm)) %>% normalize_emm_cis()
emm_plot <- emm_df %>%
  mutate(
    Type = factor(Type, levels = c("Dived","Undived")),
    TransectOrder = factor(TransectOrder, levels = c("A","B")),
    x = as.numeric(TransectOrder)
  )

delta_lab <- emm_plot %>%
  select(Type, TransectOrder, response) %>%
  pivot_wider(names_from = TransectOrder, values_from = response) %>%
  mutate(pct = 100 * (B - A) / A,
         label = paste0(ifelse(pct >= 0, "+", ""), sprintf("%.0f%%", pct))) %>%
  left_join(
    emm_plot %>% group_by(Type) %>% summarise(x_mid = mean(x), y_mid = mean(response), .groups = "drop"),
    by = "Type"
  )

p_eff <- ggplot(emm_plot, aes(x = x, y = response, color = Type, group = Type)) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = Type),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  geom_text(
    data = delta_lab,
    aes(x = x_mid, y = y_mid + 0.05 * max(emm_plot$response),
        label = label, color = Type),
    size = 3.2,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_x_continuous(breaks = c(1, 2), labels = c("A", "B")) +
  scale_color_manual(values = reef_cols) +
  scale_fill_manual(values  = reef_cols, guide = "none") +
  labs(x = "Transect order", y = "Expected total fish") +
  theme_clean +
  theme(legend.position = "top", legend.title = element_blank())

ggsave(file.path(output_dir, "figures", "fig_total_emm_clean.png"),
       p_eff, width = 7, height = 5, dpi = 300, bg = "white")

p_eff 
##### 9. Plot: Type × Transect across Boats levels #####
qs <- quantile(totals_transect$Boats, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
names(qs) <- c("Low", "Median", "High")

emm_boats <- emmeans(final_model,
                     ~ Type * TransectOrder | Boats,
                     at = list(Boats = as.numeric(qs)),
                     type = "response")

emm_boats_df <- as.data.frame(summary(emm_boats)) %>%
  normalize_emm_cis() %>%
  mutate(Boats_level = factor(rep(names(qs), each = nrow(.) / 3),
                              levels = c("Low", "Median", "High")))

write_csv(emm_boats_df, file.path(output_dir, "emm_TypeXTransect_byBoats.csv"))

p_boats <- ggplot(emm_boats_df,
                  aes(x = TransectOrder, y = response, color = Type, group = Type)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.05) +
  facet_wrap(~ Boats_level) +
  labs(title = "Estimated means across Boats levels",
       x = "Transect order", y = "Expected total fish") +
  theme_clean

ggsave(file.path(output_dir, "plot_emm_TypeXTransect_byBoatsLevels.png"),
       p_boats, width = 9, height = 5.5, dpi = 300)
p_boats

##### 10. Export coefficients #####
fixef_tab <- tidy(final_model, effects = "fixed", conf.int = TRUE)
write_csv(fixef_tab, file.path(output_dir, "fixed_effects_logscale.csv"))

zi_coefs <- summary(final_model)$coefficients$zi
write_csv(as.data.frame(zi_coefs),
          file.path(output_dir, "zero_inflation_coefs.csv"))

message("✓ Modelling complete. Outputs saved in: ", normalizePath(output_dir))
summary(final_model)


