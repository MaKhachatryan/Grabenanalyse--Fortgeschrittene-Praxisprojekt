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


