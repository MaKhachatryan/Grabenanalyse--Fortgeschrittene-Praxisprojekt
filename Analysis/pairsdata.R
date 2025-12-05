source("environmentSetUp.R")

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

##Kinship Result

kinship_result <- og_kinship_result |>
  mutate(rel_original = rel, # keeping the original rel for global threshold
  
  # cleaning rel based on original threshold
  rel = case_when(
    rel == "First Degree"  & overlap_nsnps < 500   ~ "Uncertain",
    rel == "Second Degree" & overlap_nsnps < 2000  ~ "Uncertain",
    rel == "Third Degree"  & overlap_nsnps < 15000 ~ "Uncertain",
    TRUE ~ rel
  ),
  
  # Binary relatedness under global threshold (related/unrelated) → (0/1) 
  # - If a pair has < 15000 overlapping SNPs → 0 (unrelated)
  # - If rel = "Unrelated" → 0 (even if SNP coverage is high)
  # - If rel = 1st/2nd/3rd degree AND has ≥15000 SNPs  → 1 (related).
  related_binary = case_when(
    overlap_nsnps < 15000              ~ 0,
    rel_original == "Unrelated"        ~ 0,
    rel_original %in% c("First Degree", "Second Degree", "Third Degree") &
      overlap_nsnps >= 15000           ~ 1,
    TRUE                               ~ NA_real_
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

tomb_success <- tomb_parameter |>
  group_by(Tomb) |>
  summarise(
    sample_success = first(percentage_of_successful_samples),
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
    tomb_success |>
      rename(
        tomb1 = Tomb,
        sample_success1 = sample_success
      ),
    by = "tomb1"
  ) |>
  
  # join for individual 2's tomb
  left_join(
    tomb_success |>
      rename(
        tomb2 = Tomb,
        sample_success2 = sample_success
      ),
    by = "tomb2"
  )


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
    
    same_tomb = tomb1 == tomb2,
    
    sample_success_min = pmin(sample_success1, sample_success2, na.rm = TRUE)
  )

pairs_ind <- pairs_all |> filter(tomb1 != "Elateia T46 56 62", tomb2 != "Elateia T46 56 62")

pairs_group <- pairs_all |> filter(!tomb1 %in% c("Elateia T46", "Elateia T56", "Elateia T62"),
                                   !tomb2 %in% c("Elateia T46", "Elateia T56", "Elateia T62"))
