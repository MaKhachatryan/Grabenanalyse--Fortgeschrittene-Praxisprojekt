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
              "here")

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
