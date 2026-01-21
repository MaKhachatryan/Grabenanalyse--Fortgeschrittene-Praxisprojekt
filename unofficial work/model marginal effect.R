source("environmentSetUp.R")
source("Analysis/pairsdata.R")
library(marginaleffects)
library(collapse)

imp_binary_separate_nobone <- readRDS("unofficial work/models save/imp_binary_separate_nobone.rds")
imp_binary_group_nobone <- readRDS("unofficial work/models save/imp_binary_group_nobone.rds")

##marginal effects for separate tomb

###by sex group

me_separate_sex <- avg_predictions(
  model = imp_binary_separate_nobone,
  by = "sex_pair",
  type = "response",# probability scale
  newdata = pairs_ind, 
  re_formula = NA, # population-level (no random effects)
  allow_new_levels = TRUE
)

me_separate_sex

###by age group

me_separate_age <- avg_predictions(
  model = imp_binary_separate_nobone,
  by = "age_pair",
  type = "response",        
  newdata = pairs_ind,
  re_formula = NA,          
  allow_new_levels = TRUE
)

me_separate_age

###by same_tomb
me_separate_tomb <- avg_predictions(
  model = imp_binary_separate_nobone,
  by = "same_tomb",
  type = "response",        
  newdata = pairs_ind,
  re_formula = NA,          
  allow_new_levels = TRUE
)

me_separate_tomb

##for plotting
me_sex_sep  <- me_separate_sex  |> 
  rename(term = sex_pair) |>
  mutate(group = "Sex pair")
me_age_sep  <- me_separate_age  |> 
  rename(term = age_pair) |>
  mutate(group = "Age pair")
me_tomb_sep <- me_separate_tomb |> 
  rename(term = same_tomb) |>
  mutate(
    group = "Same tomb",
    term = if_else(term, "Same tomb", "Different tomb")
    )

me_all_sep <- bind_rows(me_sex_sep, me_age_sep, me_tomb_sep)

ggplot(me_all_sep,
       aes(x = estimate, y = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  facet_wrap(~ group, scales = "free_y") +
  labs(x = "Predicted probability", y = NULL) +
  theme_minimal()


##marginal effect plots 
p_sex_sep <- ggplot(me_sex_sep, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.2
  ) +
  labs(
    title = "Marginal effects by sex pair (separate tombs)",
    x = "Predicted probability",
    y = NULL
  ) +
  theme_minimal(base_size = 14)


p_age_sep <- ggplot(me_age_sep, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.2
  ) +
  labs(
    title = "Marginal effects by age pair (separate tombs)",
    x = "Predicted probability",
    y = NULL
  ) +
  theme_minimal(base_size = 14)


p_tomb_sep <- ggplot(me_tomb_sep, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.2
  ) +
  labs(
    title = "Marginal effects by burial context (separate tombs)",
    x = "Predicted probability",
    y = NULL
  ) +
  theme_minimal(base_size = 14)


##save the plots 


##marginal effects for grouped tombs (Elateia T46, 56, 62)

###by sex group

me_group_sex <- avg_predictions(
  model = imp_binary_group_nobone,
  by = "sex_pair",
  type = "response",# probability scale
  newdata = pairs_group, 
  re_formula = NA, # population-level (no random effects)
  allow_new_levels = TRUE
)

me_group_sex

###by age group

me_group_age <- avg_predictions(
  model = imp_binary_group_nobone,
  by = "age_pair",
  type = "response",        
  newdata = pairs_group,
  re_formula = NA,          
  allow_new_levels = TRUE
)

me_group_age

###by same_tomb
me_group_tomb <- avg_predictions(
  model = imp_binary_group_nobone,
  by = "same_tomb",
  type = "response",        
  newdata = pairs_group,
  re_formula = NA,          
  allow_new_levels = TRUE
)

me_group_tomb

##for plotting

me_sex_gr  <- me_group_sex  |> 
  rename(term = sex_pair) |>
  mutate(group = "Sex pair")
me_age_gr  <- me_group_age  |> 
  rename(term = age_pair) |>
  mutate(group = "Age pair")
me_tomb_gr <- me_group_tomb |> 
  rename(term = same_tomb) |>
  mutate(
    group = "Same tomb",
    term = if_else(term, "Same tomb", "Different tomb")
  )

me_all_gr <- bind_rows(me_sex_gr, me_age_gr, me_tomb_gr)

ggplot(me_all_gr,
       aes(x = estimate, y = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  facet_wrap(~ group, scales = "free_y") +
  labs(x = "Predicted probability", y = NULL) +
  theme_minimal()

##marginal effect plots 
p_sex_gr <- ggplot(me_sex_gr, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.2
  ) +
  labs(
    title = "Marginal effects by sex pair (grouped tombs)",
    x = "Predicted probability",
    y = NULL
  ) +
  theme_minimal(base_size = 14)


p_age_gr <- ggplot(me_age_gr, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.2
  ) +
  labs(
    title = "Marginal effects by age pair (grouped tombs)",
    x = "Predicted probability",
    y = NULL
  ) +
  theme_minimal(base_size = 14)


p_tomb_gr <- ggplot(me_tomb_gr, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.2
  ) +
  labs(
    title = "Marginal effects by burial context (grouped tombs)",
    x = "Predicted probability",
    y = NULL
  ) +
  theme_minimal(base_size = 14)


##save the plots 


##Making a combined table for visual
me_all_sep_df <- me_all_sep |> 
  as.data.frame()

me_all_gr_df <- me_all_gr |> 
  as.data.frame() 


me_combined <- me_all_sep_df |>
  select(
    category = term,
    estimate_separate = estimate,
    low_sep = conf.low,
    high_sep = conf.high
  ) |>
  mutate(
    estimate_separate = sprintf("%.5f", estimate_separate),
    CI_separate = paste0(
      "[", signif(low_sep, 3), ", ", signif(high_sep, 3), "]"
    )
  ) |>
  select(category, estimate_separate, CI_separate) |>
  left_join(
    me_all_gr_df |>
      select(
        category = term,
        estimate_group = estimate,
        low_gr = conf.low,
        high_gr = conf.high
      ) |>
      mutate(
        estimate_group = sprintf("%.5f", estimate_group),
        CI_group = paste0(
          "[", signif(low_gr, 3), ", ", signif(high_gr, 3), "]"
        )
      ) |>
      select(category, estimate_group, CI_group),
    by = "category"
  )


me_combined

#save the table 
saveRDS(me_combined, "unofficial work/me_tables/me_combined.rds")


