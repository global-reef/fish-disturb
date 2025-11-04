### run all scripts to update analysis

######### 0. Set Analysis Date & Create Output Folder #######

# Enter the date for this analysis
analysis_date <- "2025_11_04"  # Update these  for each analysis run
# file path (adjust date for correct date)
file_path <- "~/Documents/1_GLOBAL REEF/0_PROJECTS/FishDisturb/fish-disturbance/2025.11.04_fish-disturb-data.csv"
# --------------------------
raw_fish <- read.csv(file_path, stringsAsFactors=TRUE, strip.white=TRUE) 


# Create a folder named with the date inside the working directory
output_dir <- file.path(getwd(), paste0("Analysis_", analysis_date))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# packages
# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(lubridate)
library(stringr)
library(forcats)
library(tidyverse)



### running analysis 
# 01 clean 
source("~/Documents/1_GLOBAL REEF/0_PROJECTS/FishDisturb/fish-disturbance/01_CLEAN.R")

# 01.1 explore (based on Zuur et al. 2010)
source("~/Documents/1_GLOBAL REEF/0_PROJECTS/FishDisturb/fish-disturbance/01.1_EXPLORE.R")

# 02 total fish model 
source("~/Documents/1_GLOBAL REEF/0_PROJECTS/FishDisturb/fish-disturbance/02_MODEL_ALLFISH.R")

# 03 functional groups 
source("~/Documents/1_GLOBAL REEF/0_PROJECTS/FishDisturb/fish-disturbance/03_FUNGROUPS.R")

# 04 species specific models 
source("~/Documents/1_GLOBAL REEF/0_PROJECTS/FishDisturb/fish-disturbance/04_MODEL_SPP.R")