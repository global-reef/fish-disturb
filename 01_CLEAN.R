library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(forcats)
library(ggplot2)

clean_data_disturb <- function(file_path) {
  df <- read.csv(file_path, stringsAsFactors = TRUE, strip.white = TRUE)
  
  # Drop empty rows/cols
  df[df == ""] <- NA
  df <- df[, colSums(!is.na(df)) > 0]
  df <- df[rowSums(!is.na(df)) > 0, ]
  
  # Standardize key columns from this raw structure
  df <- df %>%
    rename(
      DU              = D.U,
      survey_code     = SURVEY.CODE,
      path            = PATH,        # DEEP or SHALLOW if used
      transect        = n,           # A or B
      Site            = SITE,
      Date_raw        = DATE..mm.dd.yyyy.,
      Depth           = AVG.DEPTH..m.,
      Time_raw        = TIME..24hr.,
      Duration        = DURATION..m.,
      Vis             = VIS,
      Weather         = WEATHER..1.5.,
      Current         = CURRENT,
      Boats           = BOATS,
      Researcher      = RESEARCHERS
    ) %>%
    mutate(
      Type = recode(as.character(DU), "D" = "Dived", "U" = "Undived"),
      Type = factor(Type, levels = c("Dived", "Undived")),
      TransectOrder = factor(as.character(transect), levels = c("A","B")),
      Date = as.Date(as.character(Date_raw), format = "%m/%d/%Y"),
      Time = format(as.POSIXct(as.character(Time_raw), format = "%H:%M"), "%H:%M")
    )
  
  # Fix 0025 years to 2025
  df <- df %>%
    mutate(Date = if_else(year(Date) == 25,
                          as.Date(make_date(2025, month(Date), mday(Date))),
                          Date))
  
  # Consistent species column names and renames
  df <- df %>%
    rename(
      sml_Grouper = Grouper.30,
      lrg_Grouper = Grouper.30.1,
      sml_Snapper = Snapper.30,
      lrg_Snapper = Snapper.30.1
    )
  
  # Species columns present in this file
  species_cols <- c(
    "Parrotfish","Rabbitfish","Butterflyfish","Angelfish","Cleaner_Wrasse",
    "Batfish","Thicklip","Red_Breast","Slingjaw","Sweetlips","Squirrel.Soldier",
    "Triggerfish","Porcupine.Puffer","Ray",
    "sml_Snapper","lrg_Snapper",
    "Eel","Trevally","Emperorfish",
    "sml_Grouper","lrg_Grouper","Barracuda","Surgeonfish"
  )
  species_cols <- intersect(species_cols, colnames(df))
  df[species_cols] <- lapply(df[species_cols], function(x) as.numeric(as.character(x)))
  
  # Long format
  fish_long <- df %>%
    pivot_longer(
      cols = all_of(species_cols),
      names_to = "Species",
      values_to = "Count"
    ) %>%
    mutate(
      Count = replace_na(Count, 0),
      # round up fractional entries if any
      Count = ceiling(Count)
    )
  
  # Functional groups as per pelagic pipeline, no wrasse merge
  functional_groups <- tibble::tribble(
    ~Species,               ~Functional_Group,
    "Parrotfish",           "Grazer",
    "Rabbitfish",           "Grazer",
    "Butterflyfish",        "Grazer",
    "Angelfish",            "Invertivore",
    "Cleaner_Wrasse",       "Invertivore",
    "Batfish",              "Invertivore",
    "Thicklip",             "Invertivore",
    "Red_Breast",           "Invertivore",
    "Slingjaw",             "Invertivore",
    "Sweetlips",            "Invertivore",
    "Squirrel.Soldier",     "Invertivore",
    "Triggerfish",          "Invertivore",
    "Porcupine.Puffer",     "Mesopredator",
    "Ray",                  "Mesopredator",
    "sml_Snapper",          "Mesopredator",
    "lrg_Snapper",          "HTLP",
    "Eel",                  "Mesopredator",
    "Trevally",             "HTLP",
    "Emperorfish",          "Mesopredator",
    "sml_Grouper",          "Mesopredator",
    "lrg_Grouper",          "HTLP",
    "Barracuda",            "HTLP",
    "Surgeonfish",          "Grazer"
  )
  
  fish_long <- fish_long %>%
    left_join(functional_groups, by = "Species") %>%
    mutate(
      Functional_Group = factor(Functional_Group,
                                levels = c("Grazer","Invertivore","Mesopredator","HTLP"),
                                ordered = TRUE)
    )
  
  # Create IDs useful for pairing and random effects
  fish_long <- fish_long %>%
    mutate(
      survey_pair = paste(survey_code, path, Date, sep = "_"),
      survey_id   = paste(survey_pair, TransectOrder, sep = "_")
    )
  
  
  fish_long
}

# Run cleaner
fish_long <- clean_data_disturb(file_path)


# Sanity check plot: A vs B by Dived vs Undived
ggplot(totals_transect,
       aes(x = TransectOrder, y = Total, group = survey_pair, color = Type)) +
  geom_line(alpha = 0.4) +
  geom_point(alpha = 0.7) +
  labs(x = "Transect order", y = "Total fish", color = "Site type")



# Basic totals
totals_transect <- fish_long %>%
  group_by(Type, Site, Date, TransectOrder, survey_pair) %>%
  summarise(Total = sum(Count, na.rm = TRUE), .groups = "drop")

# Functional group totals
totals_group <- fish_long %>%
  group_by(Type, Site, Date, TransectOrder, Functional_Group, survey_pair) %>%
  summarise(GroupTotal = sum(Count, na.rm = TRUE), .groups = "drop")

# 1) apply standardizations
fish_long <- fish_long %>%
  mutate(
    Site = as.character(Site),
    Site = stringr::str_squish(Site),
    Site = dplyr::recode(Site,
                         "Japanese Garden"  = "Japanese Gardens",
                         "Tanote"           = "Tanote Bay",
                         "Grapeview 003"    = "Grapeview",
                         "N Sai Nuan 0012"  = "N Sai Nuan 012",
                         "S Nangyuan 001"   = "Nang Yuan S 001",
                         "Buddha Wall 010"  = "Buddha Rock S 010"
    ),
    Site = factor(Site)
  )

# 2) drop Laem Thian S 007 from analyses
drop_sites <- c("Laem Thian S 007")
fish_long <- fish_long %>% filter(!(Site %in% drop_sites))

# Save cleaned
write_csv(fish_long, file.path(output_dir, paste0("fish_long_", analysis_date, ".csv")))
write_csv(totals_transect, file.path(output_dir, paste0("totals_transect_", analysis_date, ".csv")))
write_csv(totals_group, file.path(output_dir, paste0("totals_group_", analysis_date, ".csv")))

# Helper for saving plots
save_plot <- function(p, name, w=8, h=6) ggsave(file.path(output_dir, paste0(name, "_", analysis_date, ".png")), p, width=w, height=h, dpi=300)
