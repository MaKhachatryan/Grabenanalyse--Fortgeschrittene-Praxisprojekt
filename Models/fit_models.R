# Load environment and packages
source("environmentSetUp.R")

# Load processed pair data
source("Analysis/pairsdata.R")

# Keep only pairs with 15000 SNP overlap
pairs_ind <- pairs_ind[pairs_ind$overlap_nsnps >= 15000, ]
pairs_group <- pairs_group[pairs_group$overlap_nsnps >= 15000, ]


# ------------------------------------------------------------
# Model 1: Binary relatedness – separate category pairs
# ------------------------------------------------------------

imp_binary_separate_cat <- brm(
  related_binary ~ sex_pair + age_pair + bone_quality + same_tomb +
    (1 | mm(individual1_id, individual2_id)) +
    (1 | mm(tomb1, tomb2)),
  
  data = pairs_ind,
  family = bernoulli(),
  
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

# Inspect model
summary(imp_binary_separate_cat)

# Save fitted model
saveRDS(
  imp_binary_separate_cat,
  file = "Models/saved_models/imp_binary_separate_cat.rds"
)


# ------------------------------------------------------------
# Model 2: Binary relatedness – group version
# ------------------------------------------------------------

imp_binary_group_cat <- brm(
  related_binary ~ sex_pair + age_pair + bone_quality + same_tomb +
    (1 | mm(individual1_id, individual2_id)) +
    (1 | mm(tomb1, tomb2)),
  
  data = pairs_group,
  family = bernoulli(),
  
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

# Inspect model
summary(imp_binary_group_cat)

# Save fitted model
saveRDS(
  imp_binary_group_cat,
  file = "Models/saved_models/imp_binary_group_cat.rds"
)

