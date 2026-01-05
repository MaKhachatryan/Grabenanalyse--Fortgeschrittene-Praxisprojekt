source("environmentSetUp.R")
source("Analysis/pairsdata.R")

# Try fit and save model

# **Remark:**
  
#- Only model the pairs over 15k in Overlapping SNPs
#- Create new variable to collapse skeletal elements into good and bad bones 
# only (high-quality: petrous, low-quality: others)
#- Use gaussian distrbution for numeric target (kinship coefficient) eventhough
# it's skewed

# Process data
pairs_ind <- pairs_ind[pairs_ind$overlap_nsnps >= 15000, ]
pairs_group <- pairs_group[pairs_group$overlap_nsnps >= 15000, ]


### Separated tomb
#### Binary separated (with categorical bones)
binary_separate_cat <- brm(related_binary ~ sex_pair + age_pair + bone_quality + same_tomb +
                             (1 | mm(individual1_id, individual2_id)) +
                             (1 | mm(tomb1, tomb2)),
                           data = pairs_ind,
                           family = bernoulli(),
                           seed = 123)

summary(binary_separate_cat)

##----------- Improve ------
##### Improved version
imp_binary_separate_cat <- brm(related_binary ~ sex_pair + age_pair + bone_quality + same_tomb +
                                 (1 | mm(individual1_id, individual2_id)) +
                                 (1 | mm(tomb1, tomb2)),
                               data = pairs_ind,
                               family = bernoulli(),
                               control = list(
                                 adapt_delta = 0.999,   # The "Baby Steps" (Default 0.8)
                                 max_treedepth = 15     # Allow longer trajectories (Default 10)
                               ),
                               
                               # INCREASE ITERATIONS
                               # 6000 total = 3000 warmup (thrown away) + 3000 sampling (kept)
                               iter = 6000, 
                               warmup = 3000, 
                               
                               chains = 4, 
                               cores = 4, 
                               seed = 123)

summary(imp_binary_separate_cat)
saveRDS(imp_binary_separate_cat, "unofficial work/models save/imp_binary_separate_cat.rds")


##### Try out the improved model of group version
imp_binary_group_cat <- brm(related_binary ~ sex_pair + age_pair + bone_quality + same_tomb +
                              (1 | mm(individual1_id, individual2_id)) +
                              (1 | mm(tomb1, tomb2)),
                            data = pairs_group,
                            family = bernoulli(),
                            control = list(
                              adapt_delta = 0.999,   # The "Baby Steps" (Default 0.8)
                              max_treedepth = 15     # Allow longer trajectories (Default 10)
                            ),
                            
                            # INCREASE ITERATIONS
                            # 6000 total = 3000 warmup (thrown away) + 3000 sampling (kept)
                            iter = 6000, 
                            warmup = 3000, 
                            
                            chains = 4, 
                            cores = 4, 
                            seed = 123)

summary(imp_binary_group_cat)
saveRDS(imp_binary_group_cat, "unofficial work/models save/imp_binary_group_cat.rds")

### try out numeric target
imp_numeric_sep_cat <- brm(
  kinship_coefficient ~ sex_pair + age_pair + bone_quality + same_tomb +
    (1 | mm(individual1_id, individual2_id)) +
    (1 | mm(tomb1, tomb2)),
  data = pairs_ind,
  
  family = student(),
  
  control = list(
    adapt_delta = 0.999,
    max_treedepth = 15
  ),
  iter = 6000,
  warmup = 3000,
  chains = 4,
  cores = 4,
  seed = 123
)

summary(imp_numeric_sep_cat)
saveRDS(imp_numeric_sep_cat, "unofficial work/models save/imp_numeric_sep_cat.rds")
##----------------------

#### Binary separated (with Overlapping SNPs instead of bones)
binary_separate_ov <- brm(related_binary ~ sex_pair + age_pair + overlap_nsnps + same_tomb +
                            (1 | mm(individual1_id, individual2_id)) +
                            (1 | mm(tomb1, tomb2)),
                          data = pairs_ind,
                          family = bernoulli(),
                          seed = 123)

summary(binary_separate_ov)



#### Numeric separated (with categorical bones)
numeric_separate_cat <- brm(kinship_coefficient ~ sex_pair + age_pair + bone_quality + same_tomb +
                              (1 | mm(individual1_id, individual2_id)) +
                              (1 | mm(tomb1, tomb2)),
                            data = pairs_ind,
                            family = gaussian(),
                            seed = 123)

summary(numeric_separate_cat)


#### Numeric separated (with Overlapping SNPs instead of bones)
numeric_separate_ov <- brm(kinship_coefficient ~ sex_pair + age_pair + overlap_nsnps + same_tomb +
                             (1 | mm(individual1_id, individual2_id)) +
                             (1 | mm(tomb1, tomb2)),
                           data = pairs_ind,
                           family = gaussian(),
                           seed = 123)

summary(numeric_separate_ov)




### Grouped tomb
#### Binary grouped (with categorical bones)
binary_group_cat <- brm(related_binary ~ sex_pair + age_pair + bone_quality + same_tomb +
                          (1 | mm(individual1_id, individual2_id)) +
                          (1 | mm(tomb1, tomb2)),
                        data = pairs_group,
                        family = bernoulli(),
                        seed = 123)

summary(binary_group_cat)






#### Binary grouped (with Overlapping SNPs instead of bones)
binary_group_ov <- brm(related_binary ~ sex_pair + age_pair + overlap_nsnps + same_tomb +
                         (1 | mm(individual1_id, individual2_id)) +
                         (1 | mm(tomb1, tomb2)),
                       data = pairs_group,
                       family = bernoulli(),
                       seed = 123)

summary(binary_group_ov)



#### Numeric grouped (with categorical bones)
numeric_group_cat <- brm(kinship_coefficient ~ sex_pair + age_pair + bone_quality + same_tomb +
                           (1 | mm(individual1_id, individual2_id)) +
                           (1 | mm(tomb1, tomb2)),
                         data = pairs_group,
                         family = gaussian(),
                         seed = 123)

summary(numeric_group_cat)


#### Numeric grouped (with Overlapping SNPs instead of bones)
numeric_group_ov <- brm(kinship_coefficient ~ sex_pair + age_pair + overlap_nsnps + same_tomb +
                          (1 | mm(individual1_id, individual2_id)) +
                          (1 | mm(tomb1, tomb2)),
                        data = pairs_group,
                        family = gaussian(),
                        seed = 123)

summary(numeric_group_ov)


## Save
model_names <- c(
  "binary_separate_cat", "binary_separate_ov",
  "numeric_separate_cat", "numeric_separate_ov"
  #"binary_group_cat", "binary_group_ov",
  #"numeric_group_cat", "numeric_group_ov"
)

for (m in model_names) {
  saveRDS(get(m),
          file = file.path("unofficial work/models save", paste0(m, ".rds"))
  )
}


### NEW MODEL SENSITIVITY ANALYSIS ------
# Now we try to see if the model works without the strict (15k only) threshold
source("Analysis/pairsdata.R")
pairs_ind_all <- pairs_ind[is.na(pairs_ind$low_data), ]
pairs_group_all <- pairs_group[is.na(pairs_group$low_data), ]


all_binary_separate_cat <- brm(related_binary ~ sex_pair + age_pair + bone_quality + same_tomb +
                                 (1 | mm(individual1_id, individual2_id)) +
                                 (1 | mm(tomb1, tomb2)),
                               data = pairs_ind_all,
                               family = bernoulli(),
                               control = list(
                                 adapt_delta = 0.999,   # The "Baby Steps" (Default 0.8)
                                 max_treedepth = 15     # Allow longer trajectories (Default 10)
                               ),
                               
                               # INCREASE ITERATIONS
                               # 6000 total = 3000 warmup (thrown away) + 3000 sampling (kept)
                               iter = 6000, 
                               warmup = 3000, 
                               
                               chains = 4, 
                               cores = 4, 
                               seed = 123)

summary(all_binary_separate_cat)
saveRDS(all_binary_separate_cat, "unofficial work/models save/all_binary_separate_cat.rds")


##### Group version
all_binary_group_cat <- brm(related_binary ~ sex_pair + age_pair + bone_quality + same_tomb +
                              (1 | mm(individual1_id, individual2_id)) +
                              (1 | mm(tomb1, tomb2)),
                            data = pairs_group_all,
                            family = bernoulli(),
                            control = list(
                              adapt_delta = 0.999,   # The "Baby Steps" (Default 0.8)
                              max_treedepth = 15     # Allow longer trajectories (Default 10)
                            ),
                            
                            # INCREASE ITERATIONS
                            # 6000 total = 3000 warmup (thrown away) + 3000 sampling (kept)
                            iter = 6000, 
                            warmup = 3000, 
                            
                            chains = 4, 
                            cores = 4, 
                            seed = 123)

summary(all_binary_group_cat)
saveRDS(all_binary_group_cat, "unofficial work/models save/all_binary_group_cat.rds")

## Compare
ogthreshold_vs_15k <- modelsummary(list("Model with 15k threshold" = imp_binary_separate_cat,
                  "Model with original threshold" = all_binary_separate_cat))
saveRDS(ogthreshold_vs_15k, "unofficial work/models save/ogthreshold_vs_15k")

### Delete bone quality as proxy ------
imp_binary_separate_nobone <- brm(related_binary ~ sex_pair + age_pair + same_tomb +
                                 (1 | mm(individual1_id, individual2_id)) +
                                 (1 | mm(tomb1, tomb2)),
                               data = pairs_ind,
                               family = bernoulli(),
                               control = list(
                                 adapt_delta = 0.999,   # The "Baby Steps" (Default 0.8)
                                 max_treedepth = 15     # Allow longer trajectories (Default 10)
                               ),
                               
                               # INCREASE ITERATIONS
                               # 6000 total = 3000 warmup (thrown away) + 3000 sampling (kept)
                               iter = 6000, 
                               warmup = 3000, 
                               
                               chains = 4, 
                               cores = 4, 
                               seed = 123)

summary(imp_binary_separate_nobone)
saveRDS(imp_binary_separate_nobone, "unofficial work/models save/imp_binary_separate_nobone.rds")


##### Group version
imp_binary_group_nobone <- brm(related_binary ~ sex_pair + age_pair + same_tomb +
                              (1 | mm(individual1_id, individual2_id)) +
                              (1 | mm(tomb1, tomb2)),
                            data = pairs_group,
                            family = bernoulli(),
                            control = list(
                              adapt_delta = 0.999,   # The "Baby Steps" (Default 0.8)
                              max_treedepth = 15     # Allow longer trajectories (Default 10)
                            ),
                            
                            # INCREASE ITERATIONS
                            # 6000 total = 3000 warmup (thrown away) + 3000 sampling (kept)
                            iter = 6000, 
                            warmup = 3000, 
                            
                            chains = 4, 
                            cores = 4, 
                            seed = 123)

summary(imp_binary_group_nobone)
saveRDS(imp_binary_group_nobone, "unofficial work/models save/imp_binary_group_nobone.rds")


## Compare
bone_vs_nobone <- modelsummary(list("Model with bone quality" = imp_binary_separate_cat,
                                        "Model without bone quality" = imp_binary_separate_nobone))
saveRDS(bone_vs_nobone, "unofficial work/models save/bone_vs_nobone.rds")


### Only use same_tomb parameter ------
pairs_ind <- pairs_ind[pairs_ind$overlap_nsnps >= 15000, ]
pairs_group <- pairs_group[pairs_group$overlap_nsnps >= 15000, ]

imp_binary_separate_sametomb <- brm(related_binary ~ same_tomb +
                                    (1 | mm(individual1_id, individual2_id)) +
                                    (1 | mm(tomb1, tomb2)),
                                  data = pairs_ind,
                                  family = bernoulli(),
                                  control = list(
                                    adapt_delta = 0.999,   # The "Baby Steps" (Default 0.8)
                                    max_treedepth = 15     # Allow longer trajectories (Default 10)
                                  ),
                                  
                                  # INCREASE ITERATIONS
                                  # 6000 total = 3000 warmup (thrown away) + 3000 sampling (kept)
                                  iter = 6000, 
                                  warmup = 3000, 
                                  
                                  chains = 4, 
                                  cores = 4, 
                                  seed = 123)

summary(imp_binary_separate_sametomb)
saveRDS(imp_binary_separate_sametomb, "unofficial work/models save/imp_binary_separate_sametomb.rds")


##### Group version
imp_binary_group_sametomb <- brm(related_binary ~ same_tomb +
                                 (1 | mm(individual1_id, individual2_id)) +
                                 (1 | mm(tomb1, tomb2)),
                               data = pairs_group,
                               family = bernoulli(),
                               control = list(
                                 adapt_delta = 0.999,   # The "Baby Steps" (Default 0.8)
                                 max_treedepth = 15     # Allow longer trajectories (Default 10)
                               ),
                               
                               # INCREASE ITERATIONS
                               # 6000 total = 3000 warmup (thrown away) + 3000 sampling (kept)
                               iter = 6000, 
                               warmup = 3000, 
                               
                               chains = 4, 
                               cores = 4, 
                               seed = 123)

summary(imp_binary_group_sametomb)
saveRDS(imp_binary_group_sametomb, "unofficial work/models save/imp_binary_group_sametomb.rds")


## Compare
normal_vs_only_sametomb <- modelsummary(list("Normal (no bone)" = imp_binary_separate_nobone,
                                    "Model with only same_tomb" = imp_binary_separate_sametomb))
saveRDS(normal_vs_only_sametomb, "unofficial work/models save/normal_vs_only_sametomb.rds")








