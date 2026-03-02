# Please always source this script first to set up the environment needed for 
# this project.


## ----- Loading necessary packages -----
packages <- c("tidyr",
              "dplyr",
              "stringr",
              "ggplot2",
              "readxl",
              "readr",
              "ggthemes",
              "readxl",
              "scales",
              "igraph",
              "ggraph",
              "purrr",
              "here",
              "janitor",
              "brms",
              "gt",
              "modelsummary",
              "tidybayes",
              "gridExtra")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

## ----- Loading data sets----
og_individual_metadata <- read_excel(here::here("Data/raw/Statistics project LMU_individuals metadata_Irene Hoegner.xlsx"))
og_tomb_parameter <- read_excel(here::here("Data/raw/Statistics project LMU_tomb parameters_Irene Hoegner_new.xlsx"), range = "A2:N12")
og_kinship_result <- read_excel(here::here("Data/raw/Statistics project LMU_READv2 kinship results_Irene Hoegner_corr.xlsx"))


## ----- cleaning names of variables -----
og_individual_metadata <- og_individual_metadata %>%
  janitor::clean_names() %>%
  rename(
    `1240k_snps`= x1240k_sn_ps
  )

og_tomb_parameter <- og_tomb_parameter %>%
  janitor::clean_names() %>%
  rename(
    twist_samples_used = twist_capture_samples_included_in_kinship_analysis,
    number_of_analysed_individuals = after_merging_samples_from_same_individual,
    tomb = `x1`
  )

og_kinship_result <- og_kinship_result %>%
  janitor::clean_names() %>%
  rename(
    `1st_type` = x1st_type,
    overlap_nsnps = overlap_nsn_ps,
    nsnps_x_norm = nsn_ps_x_norm
  )


