# Load environment and packages
source("environmentSetUp.R")

# Load processed data
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


# ---------------------------------------------
# Model 3: Overlap GLM
# ---------------------------------------------

#Gaussian GLM with identity link
m_overlap_glm <- glm(
  overlap_sc ~ same_tomb + sex_pair + age_pair +
    bone_quality + sample_success_min +
    length_shared,
  data   = kinship_result_m,
  family = gaussian(link = "identity")
)

# Inspect model
summary(m_overlap_glm)

# Save fitted model
saveRDS(
  m_overlap_glm,
  file = "Models/saved_models/m_overlap_glm.rds"
)


# -------------------------------------------------------------
# Model 4: Relatedness Bayesian Regression Modeling (BRM)
# -------------------------------------------------------------

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

# Inspect model
summary(m_related)

# Save fitted model
saveRDS(m_related,
        file = "Models/saved_models/m_related1.rds"
)






