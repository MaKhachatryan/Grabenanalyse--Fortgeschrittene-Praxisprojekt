source("environmentSetUp.R")
source("Analysis/pairsdata.R")


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

#################################################################################
#-------------------------MODELS-------------------#
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


#------------------1. Overlap GLM-----------------#
#Gaussian GLM with identity link
m_overlap_glm <- glm(
  overlap_sc ~ same_tomb + sex_pair + age_pair +
    bone_quality + sample_success_min +
    length_shared,
  data   = kinship_result_m,
  family = gaussian(link = "identity")
)

summary(m_overlap_glm)
saveRDS(m_overlap_glm, "unofficial work/models save/m_overlap_glm.rds")


#-------2.Relatedness Bayesian Regression Modeling (BRM)------#
priors <- c(
     prior(normal(0, 1), class = "b"),
     prior(normal(0, 1), class = "Intercept"),
     prior(exponential(2), class = "sd")
  )

m_related <- brm(related_binary_na ~ overlap_sc + low_data_fac + age_pair + 
                                     sex_pair + skeletal_element_pair + same_tomb +
                                     (1 | mm(individual1_id, individual2_id)) +
                                     (1 | mm(tomb1, tomb2)),
                                    data    = kinship_result_m,
                                    family  = bernoulli(),
                                   
                                      # INCREASE ITERATIONS
                                     # 6000 total = 3000 warmup (thrown away) + 3000 sampling (kept)
                                     iter = 6000, 
                                     warmup = 3000, 
                                     chains = 4, 
                                     cores = 4, 
                                     seed = 12,
                                     control = list(adapt_delta = 0.99,
                                                    max_treedepth = 15),
                                     prior   = priors
                  )


summary(m_related)
saveRDS(m_related, "unofficial work/models save/m_related1.rds")











