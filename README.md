# Fish Disturbance Project

Analysis of fish assemblage responses to diver presence across reef sites in the Gulf of Thailand.  
---

## Repository Structure

| File / Folder | Description |
|----------------|--------------|
| `00_RUN.R` | Master script that runs the full workflow in sequence. |
| `01_CLEAN.R` | Loads and tidies raw data; standardizes columns, factors, and metadata. |
| `01.1_EXPLORE.R` | Basic summaries and exploratory plots for initial data inspection. |
| `02_MODEL_ALLFISH.R` | Fits generalized linear mixed models (GLMMs) for total fish abundance. |
| `03_FUNGROUPS.R` | Fits GLMMs for each functional group (Grazers, Invertivores, Mesopredators, High-Trophic Predators). |
| `04_MODEL_SPP.R` | Runs species-level models for taxa with sufficient occurrences and outputs plots + CSVs. |
| `2025.11.04_fish-disturb-data.csv` | Clean input dataset used for all analyses. |
| `Analysis_2025_11_04/` | Output directory containing model summaries, figures, and CSV results. |
| `fish-disturbance.Rproj` | RStudio project file for environment management. |

---

## Quick Start

Open `fish-disturbance.Rproj` in RStudio and run:

```r
source("00_RUN.R")
```

This executes the workflow in the correct order and saves outputs to `Analysis_YYYY_MM_DD/`.

---

## Dependencies

Install required R packages:

```r
install.packages(c(
  "dplyr", "tidyr", "ggplot2", "glmmTMB", "emmeans",
  "DHARMa", "broom.mixed", "purrr", "stringr", "readr", "tibble"
))
```

---

## Analysis Summary

1. **Data Cleaning**  
   - Reads fish count data.  
   - Formats factors (`Type`, `TransectOrder`) and creates paired transect IDs.

2. **Exploration**  
   - Generates histograms, richness plots, and site summaries.

3. **Functional Group Models**  
   - GLMMs per group (`glmmTMB`, negative binomial or zero-inflated).  
   - Summaries of fixed effects (Type × TransectOrder).

4. **Species Models**  
   - Runs models for species with sufficient data.  
   - Outputs CSV summaries and effect plots.

---

## Notes

- Designed for use with `Global Reef` diver disturbance surveys.  

---

© 2025 Global Reef  
Maintainer: **Scarlett R. Taylor**
scarlett@global-reef.com
