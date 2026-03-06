# This script loads the prepared pairwise datasets and fits all models used in
# the analysis, including the main models and several robustness checks.
# The fitted models are saved as .rds files in the folder Analysis/Saved models.

# Load project environment and processed data
source("environment_setup.R")
source("Data/pairs_data.R")

# Keep only pairs with 15000 SNP overlap
pairs_ind <- pairs_ind_og[pairs_ind_og$overlap_nsnps >= 15000, ]
pairs_group <- pairs_group_og[pairs_group_og$overlap_nsnps >= 15000, ]



# 1. Main model (separated) ----
imp_binary_separate_nobone <- brm(related_binary ~ sex_pair + age_pair + same_tomb +
                                    (1 | mm(individual1_id, individual2_id)) +
                                    (1 | mm(tomb1, tomb2)),
                                  data = pairs_ind,
                                  family = bernoulli(),
                                  control = list(
                                    adapt_delta = 0.999,   # Smaller steps (Default 0.8)
                                    max_treedepth = 15     # Allow longer trajectories (Default 10)
                                  ),
                                  
                                  # INCREASE ITERATIONS
                                  # 6000 total = 3000 warmup (thrown away) + 3000 sampling (kept)
                                  iter = 6000, 
                                  warmup = 3000, 
                                  
                                  chains = 4, 
                                  cores = 4, 
                                  seed = 123)

# Inspect model
summary(imp_binary_separate_nobone)

# Save fitted model
saveRDS(
  imp_binary_separate_nobone,
  file = "Analysis/Saved models/imp_binary_separate_nobone.rds"
)


# 2. Main model (grouped) ----
imp_binary_group_nobone <- update(imp_binary_separate_nobone,
                                  newdata = pairs_group,
                                  seed = 123)

# Inspect model
summary(imp_binary_group_nobone)

# Save fitted model
saveRDS(
  imp_binary_group_nobone,
  file = "Analysis/Saved models/imp_binary_group_nobone.rds"
)

# 3. Robust check model (separated): only use "same_tomb" as predictor ----
imp_binary_separate_sametomb <- update(imp_binary_separate_nobone, 
                                       formula. = related_binary ~ same_tomb +
                                         (1 | mm(individual1_id, individual2_id)) +
                                         (1 | mm(tomb1, tomb2)),
                                       seed = 123)

# Inspect model
summary(imp_binary_separate_sametomb)

# Save fitted model
saveRDS(
  imp_binary_separate_sametomb,
  file = "Analysis/Saved models/imp_binary_separate_sametomb.rds"
)

# 4. Robust check model (grouped): only use "same_tomb" as predictor ----
imp_binary_group_sametomb <- update(imp_binary_separate_sametomb,
                                    newdata = pairs_group,
                                    seed = 123)

# Inspect model
summary(imp_binary_group_sametomb)

# Save fitted model
saveRDS(
  imp_binary_group_sametomb,
  file = "Analysis/Saved models/imp_binary_group_sametomb.rds"
)


# 5. Robust check model (separated): set 3rd degree kinship as unrelated ----
imp_binary_separate_unrel3rd <- update(imp_binary_separate_nobone, 
                                       newdata = pairs_ind_unrel3rd,
                                       seed = 123)

# Inspect model
summary(imp_binary_separate_unrel3rd)

# Save fitted model
saveRDS(
  imp_binary_separate_unrel3rd,
  file = "Analysis/Saved models/imp_binary_separate_unrel3rd.rds"
)

# 6. Robust check model (grouped): set 3rd degree kinship as unrelated ----
imp_binary_group_unrel3rd <- update(imp_binary_separate_unrel3rd,
                                    newdata = pairs_group_unrel3rd,
                                    seed = 123)

# Inspect model
summary(imp_binary_group_unrel3rd)

# Save fitted model
saveRDS(
  imp_binary_group_unrel3rd,
  file = "Analysis/Saved models/imp_binary_group_unrel3rd.rds"
)

# 7. Robust check model (separated): set prior N(0,1) ----
weak_prior <- set_prior("normal(0, 1)", class = "b")

imp_binary_separate_nobone_prior <- update(imp_binary_separate_nobone,
                                           prior = weak_prior,
                                           seed = 123)

# Inspect model
summary(imp_binary_separate_nobone_prior)

# Save fitted model
saveRDS(
  imp_binary_separate_nobone_prior,
  file = "Analysis/Saved models/imp_binary_separate_nobone_prior.rds"
)

# 8. Robust check model (grouped): set prior N(0,1) ----
imp_binary_group_nobone_prior <- update(imp_binary_group_nobone,
                                        prior = weak_prior,
                                        seed = 123)

# Inspect model
summary(imp_binary_group_nobone_prior)

# Save fitted model
saveRDS(
  imp_binary_group_nobone_prior,
  file = "Analysis/Saved models/imp_binary_group_nobone_prior.rds"
)

# 9. Robust check model (separated): within tomb pairs only ---- 
imp_binary_separate_withintomb <- update(imp_binary_separate_nobone, 
                                         formula. = related_binary ~ sex_pair + age_pair + 
                                           (1 | mm(individual1_id, individual2_id)) + (1 | tomb1),
                                         newdata = pairs_ind_withintomb,
                                         seed = 123)

# Inspect model
summary(imp_binary_separate_withintomb)

# Save fitted model
saveRDS(
  imp_binary_separate_withintomb,
  file = "Analysis/Saved models/imp_binary_separate_withintomb.rds"
)

# 10. Robust check model (grouped): within tomb pairs only ----
imp_binary_group_withintomb <- update(imp_binary_separate_withintomb, newdata = pairs_group_withintomb,
                                      seed = 123)

# Inspect model
summary(imp_binary_group_withintomb)

# Save fitted model
saveRDS(
  imp_binary_group_withintomb,
  file = "Analysis/Saved models/imp_binary_separate_withintomb.rds"
)

# 11. Robust check model (separated): within tomb pairs only and 3rd degree as unrelated ----
pairs_ind_withintomb_unrel3rd <- pairs_ind_unrel3rd[pairs_ind_unrel3rd$same_tomb == TRUE, ]
pairs_group_withintomb_unrel3rd <- pairs_group_unrel3rd[pairs_group_unrel3rd$same_tomb == TRUE, ]

imp_binary_separate_withintomb_unrel3rd <- update(imp_binary_separate_withintomb,
                                                  newdata = pairs_ind_withintomb_unrel3rd,
                                                  seed = 123)
# Inspect model
summary(imp_binary_separate_withintomb_unrel3rd)

# Save fitted model
saveRDS(
  imp_binary_separate_withintomb_unrel3rd,
  file = "Analysis/Saved models/imp_binary_separate_withintomb_unrel3rd.rds"
)

# 12. Robust check model (grouped): within tomb pairs only and 3rd degree as unrelated ----
imp_binary_group_withintomb_unrel3rd <- update(imp_binary_group_withintomb,
                                               newdata = pairs_group_withintomb_unrel3rd,
                                               seed = 123)
# Inspect model
summary(imp_binary_group_withintomb_unrel3rd)

# Save fitted model
saveRDS(
  imp_binary_group_withintomb_unrel3rd,
  file = "Analysis/Saved models/imp_binary_group_withintomb_unrel3rd.rds"
)

# 13. Robust check model (separated): standard threshold with ske.elem ----
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
# Inspect model
summary(all_binary_separate_cat)

# Save fitted model
saveRDS(
  all_binary_separate_cat,
  file = "Analysis/Saved models/all_binary_separate_cat.rds"
)

# 14. Robust check model (grouped): standard threshold with ske.elem ----
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
# Inspect model
summary(all_binary_group_cat)

# Save fitted model
saveRDS(
  all_binary_group_cat,
  file = "Analysis/Saved models/all_binary_group_cat.rds"
)


# 15. Robust check model (separated): standard threshold with main formula ----
all_binary_separate_nobone <- update(imp_binary_separate_nobone,
                                     newdata = pairs_ind_all,
                                     seed = 123)
# Inspect model
summary(all_binary_separate_nobone)

# Save fitted model
saveRDS(
  all_binary_separate_nobone,
  file = "Analysis/Saved models/all_binary_separate_nobone.rds"
)

# 16. Robust check model (grouped): standard threshold with main formula ----
all_binary_group_nobone <- update(imp_binary_group_nobone,
                                  newdata = pairs_group_all,
                                  seed = 123)
# Inspect model
summary(all_binary_group_nobone)

# Save fitted model
saveRDS(
  all_binary_group_nobone,
  file = "Analysis/Saved models/all_binary_group_nobone.rds"
)

# . Model Overlap GLM
# Gaussian GLM with identity link
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
  file = "Analysis/Saved models/m_overlap_glm.rds"
)

# . Relatedness Bayesian Regression Modeling (BRM)
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
        file = "Analysis/Saved models/m_related1.rds"
)

