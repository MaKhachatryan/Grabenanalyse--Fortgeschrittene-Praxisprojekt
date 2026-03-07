# This script cleans and processes the raw individual, tomb, and kinship data.
# It creates pairwise datasets, derives analysis variables, and prepares
# filtered datasets that are used as input for the statistical models.

# Load project environment
source("environment_setup.R")

##Individual Metadata
##- Filtering out the Lapoutsi
##- adding the tomb 46, 56, 62 as a group
##- merge the subadult and subadult?

individual_metadata <- og_individual_metadata |> 
  filter(site != "Lapoutsi", tomb %in% c("Amfissa tholos", "Tomb 31", "Tomb 36", "Tomb 46", "Tomb 50", "Tomb 56", "Tomb 62", "Tomb 67")) |> 
  mutate(
    age_estimation = case_when(
      age_estimation %in% c("subadult", "subadult?") ~ "subadult",
      age_estimation %in% c("adolescent or young adult", "adult or undefined") ~ "adult or undefined",
      TRUE ~ age_estimation   
    ),
    age_estimation = factor(
      age_estimation,
      levels = c("subadult", "adult or undefined")
    )
  )

tomb_group_rows <- individual_metadata |>
  filter(tomb %in% c("Tomb 46", "Tomb 56", "Tomb 62")) |>
  mutate(tomb = "Tomb 46 56 62")

individual_metadata <- bind_rows(individual_metadata, tomb_group_rows)



##Tomb Parameter 

combined_tombs <- og_tomb_parameter |> slice(rep(10, 4)) |> mutate(
  "mni_minimum_number_of_individuals" = c(250, 300, 350, 400),
  "percentage_of_successful_samples" = twist_samples_used / max_number_of_individuals,
  "tomb" = paste("Elateia T 46,56,62 (MNI ", c(250, 300, 350, 400), ")", sep = ""))
tomb_parameter <- rbind(og_tomb_parameter, combined_tombs)
tomb_parameter <- tomb_parameter |> rename(Tomb = 1, MNI = 2, `analysed individuals` = 12) |> slice(-c(5, 10)) |> mutate(MNI = as.integer(MNI), length_of_use_of_tomb = readr::parse_number(length_of_use_of_tomb),
                                                                                                                         percentage_of_successful_samples = round(percentage_of_successful_samples, 3)
)

tomb_parameter$MNI[8] <- 20
tomb_parameter$MNI[4] <- 70
tomb_parameter$MNI[6] <- 70

# Absolute plot
df_long <- pivot_longer(tomb_parameter, cols = c(MNI, `analysed individuals`),
                        names_to = "Type", values_to = "Count")
ggplot(df_long, aes(x = Tomb, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip()

##Kinship Result

kinship_result <- og_kinship_result |>
  mutate(rel_original = rel, # keeping the original rel for global threshold
  
  # cleaning rel based on original threshold
  rel = case_when(
    low_data == "yes" ~ "Uncertain",
    TRUE ~ rel
  ),
  
  # Binary relatedness based on the thresholded "rel" column:
  # - 1st/2nd/3rd degree = 1  (related)
  # - Unrelated          = 0  (unrelated)
  # - Uncertain          = NA (cannot classify reliably)
  related_binary = case_when(
    rel %in% c("First Degree", "Second Degree", "Third Degree") ~ 1,
    rel == "Unrelated"                                          ~ 0,
    rel == "Uncertain"                                          ~ NA_real_,
    TRUE                                                        ~ NA_real_
  ),
  
  # Degree outcome under global threshold 
  global_degree = case_when(
    overlap_nsnps < 15000 ~ "Uncertain",
    TRUE ~ rel_original      
  ),
  
  "individual1_id" = substr(ind1, 1, 6),
  "individual2_id" = substr(ind2, 1, 6)) |>
  
  # adding ind1 information
  left_join(individual_metadata |>
              rename(
                tomb1 = tomb,
                sex1 = sex,
                age1 = age_estimation,
                skeletal_element1 = skeletal_element
              ), by = c("individual1_id" = "individual_id")) |>
  
  # adding ind2 information
  left_join(individual_metadata |>
              rename(
                tomb2 = tomb,
                sex2 = sex,
                age2 = age_estimation,
                skeletal_element2 = skeletal_element
              ), by = c("individual2_id" = "individual_id")) |>
  
  dplyr::select(names(og_kinship_result), related_binary, global_degree, tomb1, 
         tomb2, sex1, sex2, age1, age2, skeletal_element1, skeletal_element2) |>
  mutate(
    "individual1_id" = substr(ind1, 1, 6),
    "individual2_id" = substr(ind2, 1, 6)
  ) |>
  filter(!startsWith(`ind1`, "LPS") & !startsWith(`ind2`, "LPS"),
         tomb1 %in% c("Amfissa tholos", "Tomb 31", "Tomb 36", "Tomb 46", "Tomb 50", "Tomb 56", "Tomb 62", "Tomb 67"),
         tomb2 %in% c("Amfissa tholos", "Tomb 31", "Tomb 36", "Tomb 46", "Tomb 50", "Tomb 56", "Tomb 62", "Tomb 67")) 

grouped_tombs <- c("Tomb 46", "Tomb 56", "Tomb 62")
group_name <- "Tomb 46 56 62"

kinship_group_pairs <- kinship_result |>
  filter(tomb1 %in% grouped_tombs | tomb2 %in% grouped_tombs) |>
  
  # Rename the Tomb columns in the DUPLICATE rows to reflect the new group
  mutate(
    tomb1 = ifelse(tomb1 %in% grouped_tombs, group_name, tomb1),
    tomb2 = ifelse(tomb2 %in% grouped_tombs, group_name, tomb2)
  )

kinship_result <- bind_rows(kinship_result, kinship_group_pairs)

tomb_lookup <- tomb_parameter |>
  group_by(Tomb) |>
  summarise(
    sample_success       = first(percentage_of_successful_samples),
    analysed_individuals = first(`analysed individuals`),   
    length_of_use        = first(length_of_use_of_tomb),
    .groups = "drop"
  )

tomb_lookup <- tomb_lookup |>
  mutate(
    Tomb = str_remove(Tomb, "\\s*\\(MNI\\s*\\d+\\)"),
    Tomb = str_replace(Tomb, "T\\s+(\\d)", "T\\1"),
    Tomb = str_replace_all(Tomb, ",", " ")
  ) |>
  
  group_by(Tomb) |>
  summarise(
    sample_success                    = first(sample_success),
    analysed_individuals              = first(analysed_individuals),
    length_of_use                     = first(length_of_use),
    .groups = "drop"
  )


#changing tomb names to match tomb parameter data set
kinship_result <- kinship_result |>
  mutate(
    tomb1 = case_when(
      tomb1 == "Amfissa tholos" ~ "Amfissa tholos",
      str_detect(tomb1, "^Tomb ") ~ paste0("Elateia T", str_remove(tomb1, "Tomb ")),
      TRUE ~ tomb1
    ),
    tomb2 = case_when(
      tomb2 == "Amfissa tholos" ~ "Amfissa tholos",
      str_detect(tomb2, "^Tomb ") ~ paste0("Elateia T", str_remove(tomb2, "Tomb ")),
      TRUE ~ tomb2
    )
  )

#joining kinship data and tomb data for the sample success
kinship_result <- kinship_result |>
  # join for individual 1's tomb
  left_join(
    tomb_lookup |>
      rename(
        tomb1            = Tomb,
        sample_success1  = sample_success,
        analysed_ind1    = analysed_individuals,
        length_of_use1   = length_of_use
      ),
    by = "tomb1"
  ) |>
  
  # join for individual 2's tomb
  left_join(
    tomb_lookup |>
      rename(
        tomb2            = Tomb,
        sample_success2  = sample_success,
        analysed_ind2    = analysed_individuals,
        length_of_use2   = length_of_use
      ),
    by = "tomb2"
  )


# Helper to combine skeletal elements from two individuals into one
combine_skeletal_elements <- function(sk1, sk2) {
  if (is.na(sk1) && is.na(sk2)) return(NA_character_)
  
  # Split by '+' and trim spaces
  split_and_clean <- function(x) {
    if (is.na(x)) return(character(0))
    parts <- unlist(strsplit(x, "\\+"))
    trimws(parts)
  }
  
  els1 <- split_and_clean(sk1)
  els2 <- split_and_clean(sk2)
  
  all_els <- unique(c(els1, els2))
  
  if (length(all_els) == 0) return(NA_character_)
  
  # Define a canonical order for elements
  element_order <- c("petrous bone",
                     "talus",
                     "phalanx",
                     "tooth",
                     "other")
  
  all_els <- all_els[order(match(all_els, element_order))]
  

  paste(all_els, collapse = "+")
}


##clean pairs data

pairs_all <- kinship_result %>%
  mutate(
    sex_pair = case_when(
      sex1 == "XX" & sex2 == "XX" ~ "XX-XX",
      sex1 == "XY" & sex2 == "XY" ~ "XY-XY",
      TRUE                        ~ "XX-XY"
    ),
    
    age_pair = case_when(
      age1 == "subadult" & age2 == "subadult" ~ "both subadult",
      age1 == "subadult" | age2 == "subadult" ~ "mixed",
      TRUE                                    ~ "both adult/undefined"
    ),
    
    skeletal_element_pair = map2_chr(
      skeletal_element1,
      skeletal_element2,
      combine_skeletal_elements
    ),
    
    same_tomb = tomb1 == tomb2
  ) 

pairs_ind_og <- pairs_all |> filter(tomb1 != "Elateia T46 56 62", tomb2 != "Elateia T46 56 62")

pairs_group_og <- pairs_all |> filter(!tomb1 %in% c("Elateia T46", "Elateia T56", "Elateia T62"),
                                   !tomb2 %in% c("Elateia T46", "Elateia T56", "Elateia T62"))

pairs_ind_og <- pairs_ind_og %>%
  mutate(
    bone_quality = case_when(
      grepl("petrous", skeletal_element1) & grepl("petrous", skeletal_element2) ~ "high-high",
      grepl("petrous", skeletal_element1) | grepl("petrous", skeletal_element2) ~ "high-low",
      TRUE ~ "low-low"
    )
  )
pairs_group_og <- pairs_group_og |>
  mutate(
    bone_quality = case_when(
      grepl("petrous", skeletal_element1) & grepl("petrous", skeletal_element2) ~ "high-high",
      grepl("petrous", skeletal_element1) | grepl("petrous", skeletal_element2) ~ "high-low",
      TRUE ~ "low-low"
    )
  )

#pairs that only has over 15k nsnps
pairs_ind <- pairs_ind_og[pairs_ind_og$overlap_nsnps >= 15000, ]
pairs_group <- pairs_group_og[pairs_group_og$overlap_nsnps >= 15000, ]

#pairs filter out by low data
pairs_ind_all <- pairs_ind_og[is.na(pairs_ind_og$low_data), ]
pairs_group_all <- pairs_group_og[is.na(pairs_group_og$low_data), ]

#under 15k overlap NSNPs and treat 3rd degree as a unrelated 
pairs_ind_unrel3rd <- pairs_ind_og %>%
  filter(overlap_nsnps >= 15000) %>%
  mutate(
    # overwrite categorical relationship
    rel = case_when(
      rel == "Third Degree" ~ "Unrelated",
      TRUE                  ~ rel
    ),
    
    # overwrite binary outcome accordingly
    related_binary = case_when(
      rel %in% c("First Degree", "Second Degree") ~ 1,
      rel == "Unrelated"                          ~ 0,
      TRUE                                       ~ NA_real_
    ),
    
    # overwrite global degree as well
    global_degree = case_when(
      global_degree == "Third Degree" ~ "Unrelated",
      TRUE                            ~ global_degree
    )
  )

pairs_group_unrel3rd <- pairs_group_og %>%
  filter(overlap_nsnps >= 15000) %>%
  mutate(
    # overwrite categorical relationship
    rel = case_when(
      rel == "Third Degree" ~ "Unrelated",
      TRUE                  ~ rel
    ),
    
    # overwrite binary outcome accordingly
    related_binary = case_when(
      rel %in% c("First Degree", "Second Degree") ~ 1,
      rel == "Unrelated"                          ~ 0,
      TRUE                                       ~ NA_real_
    ),
    
    # overwrite global degree as well
    global_degree = case_when(
      global_degree == "Third Degree" ~ "Unrelated",
      TRUE                            ~ global_degree
    )
  )

##data sets for within tomb only
pairs_ind_withintomb <- pairs_ind_og %>%
  filter(tomb1 == tomb2, overlap_nsnps >= 15000)

pairs_group_withintomb <- pairs_group_og %>%
  filter(tomb1 == tomb2, overlap_nsnps >= 15000)

###### Data preparation for modeling SNP overlap and MNAR considerations

kinship_result_m <- og_kinship_result |>
  mutate(related_binary = case_when(
    rel %in% c("First Degree", "Second Degree", "Third Degree",
               "Identical twin / same person") ~ 1,
    rel %in% c("Unrelated", "Uncertain",
               "Unrelated/consistent with third degree") ~ 0,
    TRUE ~ 0
  ),
  low_data_fac = case_when(
    low_data == "yes" ~ "yes",
    TRUE              ~ "no"     # all blanks/NA become "no"
  ) %>% 
    factor(levels = c("no", "yes")),
  
  "individual1_id" = substr(ind1, 1, 6),
  "individual2_id" = substr(ind2, 1, 6)) |>
  
  # adding ind1 information
  left_join(individual_metadata |>
              rename(
                tomb1 = tomb,
                sex1 = sex,
                age1 = age_estimation,
                skeletal_element1 = skeletal_element
              ), by = c("individual1_id" = "individual_id")) |>
  
  # adding ind2 information
  left_join(individual_metadata |>
              rename(
                tomb2 = tomb,
                sex2 = sex,
                age2 = age_estimation,
                skeletal_element2 = skeletal_element
              ), by = c("individual2_id" = "individual_id")) |>
  
  dplyr::select(names(og_kinship_result), related_binary, tomb1, related_binary, low_data_fac,
                tomb2, sex1, sex2, age1, age2, skeletal_element1, skeletal_element2) |>
  mutate(
    "individual1_id" = substr(ind1, 1, 6),
    "individual2_id" = substr(ind2, 1, 6)
  ) |>
  filter(!startsWith(`ind1`, "LPS") & !startsWith(`ind2`, "LPS"),
         tomb1 %in% c("Amfissa tholos", "Tomb 31", "Tomb 36", "Tomb 46", "Tomb 50", "Tomb 56", "Tomb 62", "Tomb 67"),
         tomb2 %in% c("Amfissa tholos", "Tomb 31", "Tomb 36", "Tomb 46", "Tomb 50", "Tomb 56", "Tomb 62", "Tomb 67")) 





#changing tomb names to match tomb parameter data set
kinship_result_m <- kinship_result_m |>
  mutate(
    tomb1 = case_when(
      tomb1 == "Amfissa tholos" ~ "Amfissa tholos",
      str_detect(tomb1, "^Tomb ") ~ paste0("Elateia T", str_remove(tomb1, "Tomb ")),
      TRUE ~ tomb1
    ),
    tomb2 = case_when(
      tomb2 == "Amfissa tholos" ~ "Amfissa tholos",
      str_detect(tomb2, "^Tomb ") ~ paste0("Elateia T", str_remove(tomb2, "Tomb ")),
      TRUE ~ tomb2
    )
  )

tomb_param <- tomb_parameter %>%
  rename(
    tomb = Tomb,                         # standard name
    sample_success = percentage_of_successful_samples,
    length_of_use = length_of_use_of_tomb
  ) %>%
  select(tomb, length_of_use, sample_success)



#joining kinship data and tomb data for the sample success and length of use
kinship_result_m <- kinship_result_m |>
  # join for individual 1's tomb
  left_join(
    tomb_param |>
      rename(
        tomb1 = tomb,
        sample_success1 = sample_success,
        length_of_use1 = length_of_use
      ),
    by = "tomb1"
  ) |>
  
  # join for individual 2's tomb
  left_join(
    tomb_param |>
      rename(
        tomb2 = tomb,
        sample_success2 = sample_success,
        length_of_use2 = length_of_use
      ),
    by = "tomb2"
  )



#making pair wise parameters
kinship_result_m <- kinship_result_m %>%
  mutate(
    sex_pair = case_when(
      sex1 == "XX" & sex2 == "XX" ~ "XX-XX",
      sex1 == "XY" & sex2 == "XY" ~ "XY-XY",
      TRUE                        ~ "XX-XY"
    ),
    
    age_pair = case_when(
      age1 == "subadult" & age2 == "subadult" ~ "both subadult",
      age1 == "subadult" | age2 == "subadult" ~ "mixed",
      TRUE                                    ~ "both adult/undefined"
    ),
    
    bone_quality = case_when(
      grepl("petrous", skeletal_element1) & grepl("petrous", skeletal_element2) ~ "high-high",
      grepl("petrous", skeletal_element1) | grepl("petrous", skeletal_element2) ~ "high-low",
      TRUE ~ "low-low"
    ),
    
    same_tomb = tomb1 == tomb2,
    
    sample_success_min = pmin(sample_success1, sample_success2, na.rm = TRUE),
    length_shared = pmin(length_of_use1, length_of_use2, na.rm = TRUE)
  )


#scaling overlap
kinship_result_m$overlap_sc <- scale(kinship_result_m$overlap_nsnps, center = TRUE, scale = TRUE)

## Make sure parameters are factors/binaries as desired:
kinship_result_m <- kinship_result_m %>%
  mutate(
    related_binary = as.integer(related_binary),
    age_pair              = factor(age_pair),
    sex_pair              = factor(sex_pair),
    bone_quality = factor(bone_quality),
    same_tomb             = factor(same_tomb)  # TRUE/FALSE -> factor
  ) %>%
  filter(
    !is.na(related_binary),
    !is.na(overlap_sc),
    !is.na(low_data_fac),
    !is.na(age_pair),
    !is.na(sex_pair),
    !is.na(bone_quality),
    !is.na(same_tomb)
  ) %>%
  droplevels()
