# This script computes and visualizes marginal effects 
# (average predicted probabilities) from the binary kinship
# models for both the separate-tomb and grouped-tomb analyses.

# Notes:
# - Predictions are calculated on the probability scale (type = "response").
# - Population-level predictions are used (re_formula = NA),
#   meaning that random effects are not included in the marginal
#   effect estimates.

# Load project environment
source("environment_setup.R")
source("Data/pairs_data.R")

imp_binary_separate_nobone <- readRDS("Analysis/Saved models/imp_binary_separate_nobone.rds")
imp_binary_group_nobone <- readRDS("Analysis/Saved models/imp_binary_group_nobone.rds")

# Marginal effects (Average Predicted Probability) for separate tomb

# By sex group
me_separate_sex <- avg_predictions(
  model = imp_binary_separate_nobone,
  by = "sex_pair",
  type = "response",# probability scale
  newdata = pairs_ind, 
  re_formula = NA, # population-level (no random effects)
  allow_new_levels = TRUE
)

me_separate_sex

# By age group
me_separate_age <- avg_predictions(
  model = imp_binary_separate_nobone,
  by = "age_pair",
  type = "response",        
  newdata = pairs_ind,
  re_formula = NA,          
  allow_new_levels = TRUE
)

me_separate_age

# By same_tomb
me_separate_tomb <- avg_predictions(
  model = imp_binary_separate_nobone,
  by = "same_tomb",
  type = "response",        
  newdata = pairs_ind,
  re_formula = NA,          
  allow_new_levels = TRUE
)

me_separate_tomb

# Preparation plotting
me_sex_sep  <- me_separate_sex  |> 
  rename(term = sex_pair) |>
  mutate(group = "Sex pair")
me_age_sep  <- me_separate_age  |> 
  rename(term = age_pair) |>
  mutate(
    group = "Age pair",
    term = case_when(
      term == "both subadult" ~ "Both Subadults",
      term == "mixed" ~ "Mixed",
      term == "both adult/undefined" ~ "Both Adults/Undefined"
    )
  )
me_tomb_sep <- me_separate_tomb |> 
  rename(term = same_tomb) |>
  mutate(
    group = "Same tomb",
    term = if_else(term, "Same Tomb", "Different Tomb")
    )

me_all_sep <- bind_rows(me_sex_sep, me_age_sep, me_tomb_sep)

# Plot the marginal effect for separated tombs for visualization of all parameters
ggplot(me_all_sep,
       aes(x = estimate, y = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  facet_wrap(~ group, scales = "free_y") +
  labs(x = "Predicted probability", y = NULL) +
  theme_minimal()


# Marginal effect plots for separated tombs
# By sex group
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

# By age group
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

# By same tomb
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



# Marginal effects for grouped tombs (Elateia T46, 56, 62)

# By sex group
me_group_sex <- avg_predictions(
  model = imp_binary_group_nobone,
  by = "sex_pair",
  type = "response",# probability scale
  newdata = pairs_group, 
  re_formula = NA, # population-level (no random effects)
  allow_new_levels = TRUE
)

me_group_sex

# By age group
me_group_age <- avg_predictions(
  model = imp_binary_group_nobone,
  by = "age_pair",
  type = "response",        
  newdata = pairs_group,
  re_formula = NA,          
  allow_new_levels = TRUE
)

me_group_age

# By same_tomb
me_group_tomb <- avg_predictions(
  model = imp_binary_group_nobone,
  by = "same_tomb",
  type = "response",        
  newdata = pairs_group,
  re_formula = NA,          
  allow_new_levels = TRUE
)

me_group_tomb

# Preparation for plotting

me_sex_gr  <- me_group_sex  |> 
  rename(term = sex_pair) |>
  mutate(group = "Sex pair")
me_age_gr  <- me_group_age  |> 
  rename(term = age_pair) |>
  mutate(
    group = "Age pair",
    term = case_when(
      term == "both subadult" ~ "Both Subadults",
      term == "mixed" ~ "Mixed",
      term == "both adult/undefined" ~ "Both Adults/Undefined"
    )
  )
me_tomb_gr <- me_group_tomb |> 
  rename(term = same_tomb) |>
  mutate(
    group = "Same tomb",
    term = if_else(term, "Same Tomb", "Different Tomb")
  )

me_all_gr <- bind_rows(me_sex_gr, me_age_gr, me_tomb_gr)

# Plot the marginal effect for grouped tombs for visualization of all parameters
ggplot(me_all_gr,
       aes(x = estimate, y = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  facet_wrap(~ group, scales = "free_y") +
  labs(x = "Predicted probability", y = NULL) +
  theme_minimal()

# Marginal effect plots for grouped tombs
# By sex group
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

# By age group
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

# By same tomb
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



# Making a combined table for visual
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

# Save the table 
saveRDS(me_combined, "Analysis/Saved models/me_combined.rds")

# Combine marginal effects from both models
me_plot <- bind_rows(
  me_all_sep %>% mutate(model = "Separated Model"),
  me_all_gr  %>% mutate(model = "Grouped Model")
) %>%
  group_by(group) %>%
  mutate(term = factor(term, levels = rev(unique(term)))) %>%
  ungroup()

# One combined plot (two models in one figure)
dodge <- position_dodge(width = 0.55)

p_one <- ggplot(me_plot, aes(x = estimate, y = term, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6, alpha = 0.6) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2, position = dodge, linewidth = 0.8) +
  geom_point(position = dodge, size = 3) +
  facet_wrap(~ group, scales = "free_y", ncol = 1) +
  labs(
    x = "Predicted probability",
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

p_one

# Save the combined marginal effect plot
ggsave(filename = "Plots/marginal_effects_combined.png", 
       plot = p_one, 
       width = 10, 
       height = 8, 
       dpi = 300 )
