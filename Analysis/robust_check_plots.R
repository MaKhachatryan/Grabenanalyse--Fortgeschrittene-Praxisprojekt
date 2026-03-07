# This script performs robustness checks for the binary kinship models by 
# comparing the main model specification with several alternative model variants 
# for both the separated-tomb and grouped-tomb analyses.
# Coefficients are plotted on the log-odds scale. The main aim is to evaluate 
# whether the conclusions of the main models remain stable across different 
# assumptions, sample restrictions, and prior specifications.

# Load project environment
source("environment_setup.R")

## -- load all the separate models --
imp_binary_separate_nobone <- readRDS("Analysis/Saved models/imp_binary_separate_nobone.rds")
imp_binary_separate_withintomb <- readRDS("Analysis/Saved models/imp_binary_separate_withintomb.rds")
imp_binary_separate_unrel3rd <- readRDS("Analysis/Saved models/imp_binary_separate_unrel3rd.rds")
imp_binary_separate_sametomb <- readRDS("Analysis/Saved models/imp_binary_separate_sametomb.rds")
imp_binary_separate_nobone_prior <- readRDS("Analysis/Saved models/imp_binary_separate_nobone_prior.rds")

## -- load all the grouped models -- 
imp_binary_group_nobone <- readRDS("Analysis/Saved models/imp_binary_group_nobone.rds")
imp_binary_group_withintomb <- readRDS("Analysis/Saved models/imp_binary_group_withintomb.rds")
imp_binary_group_unrel3rd <- readRDS("Analysis/Saved models/imp_binary_group_unrel3rd.rds")
imp_binary_group_sametomb <- readRDS("Analysis/Saved models/imp_binary_group_sametomb.rds")
imp_binary_group_nobone_prior <- readRDS("Analysis/Saved models/imp_binary_group_nobone_prior.rds")


extract_brms_estimates <- function(model, model_name, probs = c(.025, .975)) {
  
  # ---- fixed effects ----
  fx <- as.data.frame(fixef(model, probs = probs)) |>
    rownames_to_column("term") |>
    transmute(
      model = model_name,
      component = "Fixed effects",
      term = term,
      estimate = Estimate,
      conf.low = Q2.5,
      conf.high = Q97.5
    )
  
  # ---- random effects SDs (robust) ----
  vc <- VarCorr(model, probs = probs)
  
  sd_df <- imap_dfr(vc, function(grp, grp_name) {
    
    # some groups might not have sd (rare, but safe)
    if (!("sd" %in% names(grp))) return(tibble())
    
    sd_obj <- grp$sd
    
    # Case 1: matrix/array with rownames = coef, colnames = Estimate/Q2.5/Q97.5
    if (!is.null(dim(sd_obj))) {
      sd_tbl <- as.data.frame(sd_obj) |>
        rownames_to_column("coef")
      
      # brms usually uses these column names; if not, this still shows you where it fails
      return(sd_tbl |>
               transmute(
                 model = model_name,
                 component = "Random effects (SD)",
                 term = paste0("SD: ", grp_name, " (", coef, ")"),
                 estimate = Estimate,
                 conf.low = Q2.5,
                 conf.high = Q97.5
               ))
    }
    
    # Case 2: named numeric vector (often Intercept only)
    # names typically: Estimate, Q2.5, Q97.5
    if (is.numeric(sd_obj) && !is.null(names(sd_obj))) {
      return(tibble(
        model = model_name,
        component = "Random effects (SD)",
        term = paste0("SD: ", grp_name, " (Intercept)"),
        estimate = unname(sd_obj["Estimate"]),
        conf.low = unname(sd_obj["Q2.5"]),
        conf.high = unname(sd_obj["Q97.5"])
      ))
    }
    
    # Fallback (in case structure is unusual)
    tibble()
  })
  
  bind_rows(fx, sd_df)
}

# Clean the label
clean_term_labels <- function(df) {
  df |>
    mutate(
      term = case_when(
        term == "Intercept" | term == "(Intercept)" ~ "Intercept",
        term == "age_pairmixed" ~ "Mixed",
        term == "age_pairbothsubadult" ~ "Both Subadults",
        term == "sex_pairXYMXY" ~ "XY-XY",
        term == "sex_pairXXMXY" ~ "XX-XY",
        term == "same_tombTRUE" ~ "Same Tomb",
        TRUE ~ term
      ),
      
      term = case_when(
        str_detect(term, "^SD:\\s*mmindividual1_idindividual2_id") ~ "SD (Individual pair intercept)",
        str_detect(term, "^SD:\\s*mmtomb1tomb2") ~ "SD (Tomb pair intercept)",
        TRUE ~ term
      ),
      
      group = case_when(
        group == "Random effects (SD)" ~ "Random effects (SD)",
        group == "Age pair" ~ "Age pair",
        group == "Sex pair" ~ "Sex pair",
        group == "Same tomb" ~ "Same tomb",
        group == "Intercept" ~ "Intercept",
        TRUE ~ group
      )
    )
}

robustcheck_models_separate <- list(
  "Main model" = imp_binary_separate_nobone,
  "Only same_tomb" = imp_binary_separate_sametomb,
  "Within tomb" = imp_binary_separate_withintomb,
  "3rd as unrelated" = imp_binary_separate_unrel3rd
)

robust_df_separate <- imap_dfr(
  robustcheck_models_separate,
  ~ extract_brms_estimates(.x, .y)
)

robust_df_separate <- robust_df_separate |>
  mutate(
    group = case_when(
      component == "Random effects (SD)"        ~ "Random effects (SD)",
      term == "Intercept"                       ~ "Intercept",
      str_detect(term, "^age_pair")             ~ "Age pair",
      str_detect(term, "^sex_pair")             ~ "Sex pair",
      str_detect(term, "^same_tomb")            ~ "Same tomb",
      TRUE                                      ~ "Other"
    ),
    
    model = factor(model, levels = names(robustcheck_models_separate))
  )

robust_df_separate <- clean_term_labels(robust_df_separate)

## -- plot for all of the parameters
pd <- position_dodge(width = 0.6)

ggplot(
  robust_df_separate,
  aes(x = estimate, y = term, color = model)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(position = pd, size = 1.6) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.18, position = pd, linewidth = 0.5) +
  facet_wrap(~ group, scales = "free_y", ncol = 1) +
  labs(
    x = "Estimate (log-odds scale)",
    y = NULL,
    color = "Model",
    title = "Model Robustness Check (Separated tombs)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 14),
    strip.text = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.spacing = unit(0.9, "lines"),
    plot.margin = margin(10, 20, 10, 10)  # extra right margin
  )


##-- plot for each of the group --

model_colors <- c(
  "Main model"        = "#D55E00",  # red/orange
  "Only same_tomb"    = "#009E73",  # green
  "Within tomb"      = "#0072B2",  # blue
  "3rd as unrelated" = "#CC79A7"   # purple
)


plot_one_separate <- function(g) {
  d <- robust_df_separate |> filter(group == g)
  
  ggplot(d, aes(x = estimate, y = term, color = model)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_point(position = pd, size = 2) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                   height = 0.18, position = pd, linewidth = 0.5) +
    scale_color_manual(values = model_colors, drop = FALSE) +
    labs(
      x = "Estimate (log-odds scale)",
      y = NULL,
      color = "Model",
      title = paste("Model Robustness Check (Separated):", g)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 14),
      axis.text.y = element_text(size = 10)
    )
}

##- -- plot separate tombs --

p_age_separate   <- plot_one_separate("Age pair")
p_sex_separate   <- plot_one_separate("Sex pair")
p_tomb_separate  <- plot_one_separate("Same tomb")
p_int_separate   <- plot_one_separate("Intercept")
p_re_separate    <- plot_one_separate("Random effects (SD)")

p_age_separate
p_sex_separate
p_tomb_separate
p_int_separate
p_re_separate


## -- Grouped tombs --
robustcheck_models_group <- list(
  "Main model" = imp_binary_group_nobone,
  "Only same_tomb" = imp_binary_group_sametomb,
  "Within tomb" = imp_binary_group_withintomb,
  "3rd as unrelated" = imp_binary_group_unrel3rd
)

robust_df_group <- imap_dfr(
  robustcheck_models_group,
  ~ extract_brms_estimates(.x, .y)
)

robust_df_group <- robust_df_group |>
  mutate(
    group = case_when(
      component == "Random effects (SD)"        ~ "Random effects (SD)",
      term == "Intercept"                       ~ "Intercept",
      str_detect(term, "^age_pair")             ~ "Age pair",
      str_detect(term, "^sex_pair")             ~ "Sex pair",
      str_detect(term, "^same_tomb")            ~ "Same tomb",
      TRUE                                      ~ "Other"
    ),
    
    model = factor(model, levels = names(robustcheck_models_separate))
  )

robust_df_group <- clean_term_labels(robust_df_group)

plot_one_group <- function(g) {
  d <- robust_df_group |> filter(group == g)
  
  ggplot(d, aes(x = estimate, y = term, color = model)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_point(position = pd, size = 2) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                   height = 0.18, position = pd, linewidth = 0.5) +
    scale_color_manual(values = model_colors, drop = FALSE) +
    labs(
      x = "Estimate (log-odds scale)",
      y = NULL,
      color = "Model",
      title = paste("Model Robustness Check (Grouped):", g)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 14),
      axis.text.y = element_text(size = 10)
    )
}

p_age_group   <- plot_one_group("Age pair")
p_sex_group   <- plot_one_group("Sex pair")
p_tomb_group  <- plot_one_group("Same tomb")
p_int_group   <- plot_one_group("Intercept")
p_re_group    <- plot_one_group("Random effects (SD)")

p_age_group
p_sex_group
p_tomb_group
p_int_group
p_re_group



# ---- choose only the two separate models (Main and Within tomb) ----
models_sep_main_vs_withintomb <- list(
  "Main model"   = imp_binary_separate_nobone,
  "Within tomb"  = imp_binary_separate_withintomb
)

# ---- extract estimates ----
robust_sep_main_vs_withintomb <- imap_dfr(models_sep_main_vs_withintomb, ~ extract_brms_estimates(.x, .y))

# ---- keep only age + sex fixed effects ----
robust_sep_main_vs_withintomb <- robust_sep_main_vs_withintomb |>
  filter(component == "Fixed effects") |>
  mutate(
    group = case_when(
      str_detect(term, "^age_pair") ~ "Age pair",
      str_detect(term, "^sex_pair") ~ "Sex pair",
      TRUE                          ~ NA_character_
    )
  ) |>
  filter(!is.na(group)) |>
  mutate(
    model = factor(model, levels = names(models_sep_main_vs_withintomb))
  )


pd <- position_dodge(width = 0.6)

p_sep_main_vs_withintomb <- ggplot(robust_sep_main_vs_withintomb, aes(x = estimate, y = term, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(position = pd, size = 2) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.18, position = pd, linewidth = 0.5
  ) +
  facet_wrap(~ group, scales = "free_y", ncol = 1) +
  labs(
    x = "Estimate (log-odds scale)",
    y = NULL,
    color = "Model",
    title = "Main vs Within tomb (Separate)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 11),
    axis.text.y = element_text(size = 10)
  )

p_sep_main_vs_withintomb

ggsave("Plots/robust_main_withintomb.png",   p_sep_main_vs_withintomb,  width = 6, height = 5, dpi = 400)

##For the grouped tombs
# ---- choose only the two grouped models (Main and Within tomb) ----
models_grp_main_vs_withintomb <- list(
  "Main model"   = imp_binary_group_nobone,
  "Within tomb"  = imp_binary_group_withintomb
)

# ---- extract estimates ----
robust_grp_main_vs_withintomb <- imap_dfr(
  models_grp_main_vs_withintomb,
  ~ extract_brms_estimates(.x, .y)
)

# ---- keep only age + sex fixed effects ----
robust_grp_main_vs_withintomb <- robust_grp_main_vs_withintomb |>
  filter(component == "Fixed effects") |>
  mutate(
    group = case_when(
      str_detect(term, "^age_pair") ~ "Age pair",
      str_detect(term, "^sex_pair") ~ "Sex pair",
      TRUE                          ~ NA_character_
    ),
    model = factor(model, levels = names(models_grp_main_vs_withintomb))
  ) |>
  filter(!is.na(group))

p_grp_main_vs_withintomb <- ggplot(
  robust_grp_main_vs_withintomb,
  aes(x = estimate, y = term, color = model)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(position = pd, size = 2) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.18, position = pd, linewidth = 0.5
  ) +
  facet_wrap(~ group, scales = "free_y", ncol = 1) +
  labs(
    x = "Estimate (log-odds scale)",
    y = NULL,
    color = "Model",
    title = "Main vs Within tomb (Grouped)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 11),
    axis.text.y = element_text(size = 10)
  )

p_grp_main_vs_withintomb

ggsave(
  "Plots/robust_group_main_withintomb.png",
  p_grp_main_vs_withintomb,
  width = 6, height = 5, dpi = 400
)

## --- compare the main model with the 3 as unrelated ---

models_sep_main_vs_unrel3rd <- list(
  "Main model"        = imp_binary_separate_nobone,
  "3rd as unrelated"  = imp_binary_separate_unrel3rd
)

robust_sep_main_vs_unrel3rd <- imap_dfr(
  models_sep_main_vs_unrel3rd,
  ~ extract_brms_estimates(.x, .y)
)

robust_sep_main_vs_unrel3rd <- robust_sep_main_vs_unrel3rd |>
  filter(component == "Fixed effects") |>
  mutate(
    group = case_when(
      str_detect(term, "^age_pair")  ~ "Age pair",
      str_detect(term, "^sex_pair")  ~ "Sex pair",
      str_detect(term, "^same_tomb") ~ "Same tomb",
      TRUE                           ~ NA_character_
    ),
    model = factor(model, levels = names(models_sep_main_vs_unrel3rd))
  ) |>
  filter(!is.na(group))

robust_sep_main_vs_unrel3rd <- robust_sep_main_vs_unrel3rd |>
  mutate(
    group = factor(
      group,
      levels = c("Age pair", "Sex pair", "Same tomb")
    )
  )

p_sep_main_vs_unrel3rd <- ggplot(
  robust_sep_main_vs_unrel3rd,
  aes(x = estimate, y = term, color = model)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(position = pd, size = 2) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.18, position = pd, linewidth = 0.5
  ) +
  facet_wrap(~ group, scales = "free_y", ncol = 1) +
  labs(
    x = "Estimate (log-odds scale)",
    y = NULL,
    color = "Model",
    title = "Main vs 3rd-as-unrelated (Separate)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 10)
  )

p_sep_main_vs_unrel3rd

ggsave("Plots/robust_main_unrel3rd.png",   p_sep_main_vs_unrel3rd,  
       width = 6, height = 5, dpi = 400)

##Grouped tombs
models_grp_main_vs_unrel3rd <- list(
  "Main model"        = imp_binary_group_nobone,
  "3rd as unrelated"  = imp_binary_group_unrel3rd
)

# ---- extract estimates ----
robust_grp_main_vs_unrel3rd <- imap_dfr(
  models_grp_main_vs_unrel3rd,
  ~ extract_brms_estimates(.x, .y)
)

# ---- keep age + sex + same_tomb fixed effects ----
robust_grp_main_vs_unrel3rd <- robust_grp_main_vs_unrel3rd |>
  filter(component == "Fixed effects") |>
  mutate(
    group = case_when(
      str_detect(term, "^age_pair")  ~ "Age pair",
      str_detect(term, "^sex_pair")  ~ "Sex pair",
      str_detect(term, "^same_tomb") ~ "Same tomb",
      TRUE                           ~ NA_character_
    ),
    model = factor(model, levels = names(models_grp_main_vs_unrel3rd))
  ) |>
  filter(!is.na(group)) |>
  mutate(
    group = factor(group, levels = c("Age pair", "Sex pair", "Same tomb"))
  )

p_grp_main_vs_unrel3rd <- ggplot(
  robust_grp_main_vs_unrel3rd,
  aes(x = estimate, y = term, color = model)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(position = pd, size = 2) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.18, position = pd, linewidth = 0.5
  ) +
  facet_wrap(~ group, scales = "free_y", ncol = 1) +
  labs(
    x = "Estimate (log-odds scale)",
    y = NULL,
    color = "Model",
    title = "Main vs 3rd-as-unrelated (Grouped)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 10)
  )

p_grp_main_vs_unrel3rd

ggsave(
  "Plots/robust_group_main_unrel3rd.png",
  p_grp_main_vs_unrel3rd,
  width = 9, height = 5, dpi = 400
)

## --- compare the main model (flat prior) with weakly informative prior N(0,1) ---

models_sep_flat_vs_wip <- list(
  "Flat prior (main)"                 = imp_binary_separate_nobone,
  "Weakly informative prior N(0,1)"   = imp_binary_separate_nobone_prior
)

# ---- extract estimates ----
robust_sep_flat_vs_wip <- imap_dfr(
  models_sep_flat_vs_wip,
  ~ extract_brms_estimates(.x, .y)
)

# ---- keep age + sex + same_tomb fixed effects ----
robust_sep_flat_vs_wip <- robust_sep_flat_vs_wip |>
  filter(component == "Fixed effects") |>
  mutate(
    group = case_when(
      str_detect(term, "^age_pair")  ~ "Age pair",
      str_detect(term, "^sex_pair")  ~ "Sex pair",
      str_detect(term, "^same_tomb") ~ "Same tomb",
      TRUE                           ~ NA_character_
    ),
    model = factor(model, levels = names(models_sep_flat_vs_wip))
  ) |>
  filter(!is.na(group)) |>
  mutate(
    group = factor(group, levels = c("Age pair", "Sex pair", "Same tomb"))
  )

pd <- position_dodge(width = 0.6)

p_sep_flat_vs_wip <- ggplot(
  robust_sep_flat_vs_wip,
  aes(x = estimate, y = term, color = model)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(position = pd, size = 2) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.18, position = pd, linewidth = 0.5
  ) +
  facet_wrap(~ group, scales = "free_y", ncol = 1) +
  labs(
    x = "Estimate (log-odds scale)",
    y = NULL,
    color = "Model",
    title = "Flat prior vs Weakly informative prior (Separate)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 10)
  )

p_sep_flat_vs_wip

## --- Grouped tombs ---

models_grp_flat_vs_wip <- list(
  "Flat prior (main)"                 = imp_binary_group_nobone,
  "Weakly informative prior N(0,1)"   = imp_binary_group_nobone_prior
)

# ---- extract estimates ----
robust_grp_flat_vs_wip <- imap_dfr(
  models_grp_flat_vs_wip,
  ~ extract_brms_estimates(.x, .y)
)

# ---- keep age + sex + same_tomb fixed effects ----
robust_grp_flat_vs_wip <- robust_grp_flat_vs_wip |>
  filter(component == "Fixed effects") |>
  mutate(
    group = case_when(
      str_detect(term, "^age_pair")  ~ "Age pair",
      str_detect(term, "^sex_pair")  ~ "Sex pair",
      str_detect(term, "^same_tomb") ~ "Same tomb",
      TRUE                           ~ NA_character_
    ),
    model = factor(model, levels = names(models_grp_flat_vs_wip))
  ) |>
  filter(!is.na(group)) |>
  mutate(
    group = factor(group, levels = c("Age pair", "Sex pair", "Same tomb"))
  )

pd <- position_dodge(width = 0.6)

p_grp_flat_vs_wip <- ggplot(
  robust_grp_flat_vs_wip,
  aes(x = estimate, y = term, color = model)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(position = pd, size = 2) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.18, position = pd, linewidth = 0.5
  ) +
  facet_wrap(~ group, scales = "free_y", ncol = 1) +
  labs(
    x = "Estimate (log-odds scale)",
    y = NULL,
    color = "Model",
    title = "Flat prior vs Weakly informative prior (Grouped)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 10)
  )

p_grp_flat_vs_wip

## Fixed cleaned label 

clean_fixed_term_labels <- function(df) {
  df |>
    mutate(
      term = case_when(
        term == "age_pairmixed" ~ "Mixed",
        term == "age_pairbothsubadult" ~ "Both Subadults",
        term == "sex_pairXYMXY" ~ "XY-XY",
        term == "sex_pairXXMXY" ~ "XX-XY",
        term == "same_tombTRUE" ~ "Same Tomb",
        TRUE ~ term
      ),
      # optional: nicer facet titles too (if you ever have them as factors)
      group = case_when(
        group == "Age pair" ~ "Age pair",
        group == "Sex pair" ~ "Sex pair",
        group == "Same tomb" ~ "Same tomb",
        TRUE ~ group
      )
    )
}

### 1. Separate Main vs Within Tomb
robust_sep_main_vs_withintomb <- clean_fixed_term_labels(robust_sep_main_vs_withintomb)

pd <- position_dodge(width = 0.6)

p_sep_main_vs_withintomb_clean <- ggplot(robust_sep_main_vs_withintomb, aes(x = estimate, y = term, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(position = pd, size = 2) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.18, position = pd, linewidth = 0.5
  ) +
  facet_wrap(~ group, scales = "free_y", ncol = 1) +
  labs(
    x = "Estimate (log-odds scale)",
    y = NULL,
    color = "Model",
    title = "Main vs Within tomb (Separate)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 11),
    axis.text.y = element_text(size = 10)
  )

p_sep_main_vs_withintomb_clean 

### 2. Grouped Main vs Within Tomb
robust_grp_main_vs_withintomb <- clean_fixed_term_labels(robust_grp_main_vs_withintomb)

p_grp_main_vs_withintomb_clean <- ggplot(robust_grp_main_vs_withintomb, aes(x = estimate, y = term, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(position = pd, size = 2) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.18, position = pd, linewidth = 0.5
  ) +
  facet_wrap(~ group, scales = "free_y", ncol = 1) +
  labs(
    x = "Estimate (log-odds scale)",
    y = NULL,
    color = "Model",
    title = "Main vs Within tomb (Grouped)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 11),
    axis.text.y = element_text(size = 10)
  )

p_grp_main_vs_withintomb_clean

### 3. Separate Main vs 3rd-as-unrelated
robust_sep_main_vs_unrel3rd <- clean_fixed_term_labels(robust_sep_main_vs_unrel3rd)

p_sep_main_vs_unrel3rd_clean <- ggplot(robust_sep_main_vs_unrel3rd, aes(x = estimate, y = term, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(position = pd, size = 2) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.18, position = pd, linewidth = 0.5
  ) +
  facet_wrap(~ group, scales = "free_y", ncol = 1) +
  labs(
    x = "Estimate (log-odds scale)",
    y = NULL,
    color = "Model",
    title = "Main vs 3rd-as-unrelated (Separate)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 10)
  )

p_sep_main_vs_unrel3rd_clean

### 4. Grouped Main vs 3rd-as-unrelated
robust_grp_main_vs_unrel3rd <- clean_fixed_term_labels(robust_grp_main_vs_unrel3rd)

p_grp_main_vs_unrel3rd_clean <- ggplot(robust_grp_main_vs_unrel3rd, aes(x = estimate, y = term, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(position = pd, size = 2) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.18, position = pd, linewidth = 0.5
  ) +
  facet_wrap(~ group, scales = "free_y", ncol = 1) +
  labs(
    x = "Estimate (log-odds scale)",
    y = NULL,
    color = "Model",
    title = "Main vs 3rd-as-unrelated (Grouped)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 10)
  )

p_grp_main_vs_unrel3rd_clean

### 5. Separate Flat prior vs Weak prior
robust_sep_flat_vs_wip <- clean_fixed_term_labels(robust_sep_flat_vs_wip)

p_sep_flat_vs_wip_clean <- ggplot(robust_sep_flat_vs_wip, aes(x = estimate, y = term, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(position = pd, size = 2) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.18, position = pd, linewidth = 0.5
  ) +
  facet_wrap(~ group, scales = "free_y", ncol = 1) +
  labs(
    x = "Estimate (log-odds scale)",
    y = NULL,
    color = "Model",
    title = "Flat prior vs Weakly informative prior (Separate)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 10)
  )

p_sep_flat_vs_wip_clean

### 6. Grouped Flat prior vs Weak prior
robust_grp_flat_vs_wip <- clean_fixed_term_labels(robust_grp_flat_vs_wip)

p_grp_flat_vs_wip_clean <- ggplot(robust_grp_flat_vs_wip, aes(x = estimate, y = term, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(position = pd, size = 2) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.18, position = pd, linewidth = 0.5
  ) +
  facet_wrap(~ group, scales = "free_y", ncol = 1) +
  labs(
    x = "Estimate (log-odds scale)",
    y = NULL,
    color = "Model",
    title = "Flat prior vs Weakly informative prior (Grouped)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 10)
  )

p_grp_flat_vs_wip_clean

# Save all of the clean plots 
ggsave(
  "Plots/robust_main_withintomb_sep_clean.png",
  p_sep_main_vs_withintomb_clean,
  width = 6, height = 5, dpi = 400
)
ggsave(
  "Plots/robust_group_main_withintomb_clean.png",
  p_grp_main_vs_withintomb_clean,
  width = 6, height = 5, dpi = 400
)
ggsave(
  "Plots/robust_main_unrel3rd_sep_clean.png",
  p_sep_main_vs_unrel3rd_clean,
  width = 6, height = 5, dpi = 400
)
ggsave(
  "Plots/robust_group_main_unrel3rd_clean.png",
  p_grp_main_vs_unrel3rd_clean,
  width = 9, height = 5, dpi = 400
)
ggsave(
  "Plots/robust_flat_vs_weakprior_sep_clean.png",
  p_sep_flat_vs_wip_clean,
  width = 6, height = 5, dpi = 400
)
ggsave(
  "Plots/robust_flat_vs_weakprior_group_clean.png",
  p_grp_flat_vs_wip_clean,
  width = 6, height = 5, dpi = 400
)
