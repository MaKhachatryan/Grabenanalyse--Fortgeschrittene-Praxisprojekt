# This script extracts, organizes, and visualizes the estimated coefficients 
# from the binary kinship models for both the separate-tomb and grouped-tomb 
# analyses.
# Coefficients are plotted on the log-odds scale. The output is intended to 
# provide a direct visual comparison of parameter estimates between the two 
# model specifications.

# Load project environment
source("environment_setup.R")

# Load the models (separated and grouped)
imp_binary_separate_nobone <- readRDS("Analysis/Saved models/imp_binary_separate_nobone.rds")
imp_binary_group_nobone <- readRDS("Analysis/Saved models/imp_binary_group_nobone.rds")

# Extract the coefficient from the separated tombs model
coef_sep <- broom.mixed::tidy(
  imp_binary_separate_nobone,
  effects = c("fixed", "ran_pars"),
  conf.int = TRUE
) |>
  mutate(model = "Separated Model")

# Extract the coefficient from the grouped tombs model
coef_gr <- broom.mixed::tidy(
  imp_binary_group_nobone,
  effects = c("fixed", "ran_pars"),
  conf.int = TRUE
) |>
  mutate(model = "Grouped Model")

# Combined the extracted coefficients from both models
coef_all <- bind_rows(coef_sep, coef_gr)

# Plots for visualization for all parameters
p_coef <- ggplot(coef_all,
                 aes(x = estimate,
                     y = reorder(term, estimate),
                     color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    position = position_dodge(width = 0.5),
    height = 0.2
  ) +
  labs(
    x = "Coefficient (log-odds)",
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 14)

p_coef

# Make a neater plot to visualize
# 1) format labels + set facet order
coef_fmt <- coef_all %>%
  mutate(
    # Pretty labels for fixed effects
    term_clean = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "same_tombTRUE" ~ "Same Tomb",
      term == "age_pairmixed" ~ "Mixed",
      term == "age_pairbothsubadult" ~ "Both Subadults",
      term == "sex_pairXXMXY" ~ "XX-XY",
      term == "sex_pairXYMXY" ~ "XY-XY",
      TRUE ~ term
    ),
    
    # Facet group
    facet = case_when(
      str_detect(term, "^age_pair") ~ "Age pair",
      str_detect(term, "^same_tomb") ~ "Same tomb",
      str_detect(term, "^sex_pair") ~ "Sex pair",
      term == "(Intercept)" ~ "Intercept",
      str_detect(term, "^sd__\\(") ~ "Random effects (SD)",
      TRUE ~ "Other"
    ),
    
    term_plot = case_when(
      str_detect(term, "^sd__\\(") ~ group,
      TRUE ~ term_clean
    )
  ) %>%
  filter(facet %in% c("Age pair", "Same tomb", "Sex pair", "Random effects (SD)")) %>%
  mutate(
    facet = factor(
      facet,
      levels = c("Age pair", "Same tomb", "Sex pair", "Random effects (SD)")
    )
  ) %>%
  select(model, facet, term_plot, estimate, conf.low, conf.high)

# Rename the random-effect group labels to something readable
coef_fmt <- coef_fmt %>%
  mutate(
    term_plot = case_when(
      term_plot == "mmindividual1_idindividual2_id" ~ "Individual-pair (ID1–ID2)",
      term_plot == "mmtomb1tomb2" ~ "Tomb-pair (Tomb1–Tomb2)",
      TRUE ~ term_plot
    )
  )

# 2) order rows within each facet nicely
coef_fmt <- coef_fmt %>%
  group_by(facet) %>%
  mutate(term_plot = factor(term_plot, levels = rev(unique(term_plot)))) %>%
  ungroup()

# 3) plot
dodge <- position_dodge(width = 0.55)

p_coef_clean <- ggplot(coef_fmt, aes(x = estimate, y = term_plot, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6, alpha = 0.6) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2, position = dodge, linewidth = 0.8) +
  geom_point(position = dodge, size = 3) +
  facet_wrap(~ facet, scales = "free_y", ncol = 1) +
  labs(x = "Coefficient (log-odds)", y = NULL, color = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

p_coef_clean

# Save the model output plot
ggsave(
  "Plots/model_coefficients_all.png",
  p_coef_clean,
  width = 10,
  height = 10,
  dpi = 300
)
