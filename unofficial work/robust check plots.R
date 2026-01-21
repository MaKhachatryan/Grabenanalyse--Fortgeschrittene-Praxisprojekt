source("environmentSetUp.R")
library(tibble)

## -- load all the separate models --
imp_binary_separate_nobone <- readRDS("unofficial work/models save/imp_binary_separate_nobone.rds")
imp_binary_separate_withintomb <- readRDS("unofficial work/models save/imp_binary_separate_withintomb.rds")
imp_binary_separate_unrel3rd <- readRDS("unofficial work/models save/imp_binary_separate_unrel3rd.rds")
imp_binary_separate_sametomb <- readRDS("unofficial work/models save/imp_binary_separate_sametomb.rds")

## -- load all the grouped models -- 
imp_binary_group_nobone <- readRDS("unofficial work/models save/imp_binary_group_nobone.rds")
imp_binary_group_withintomb <- readRDS("unofficial work/models save/imp_binary_group_withintomb.rds")
imp_binary_group_unrel3rd <- readRDS("unofficial work/models save/imp_binary_group_unrel3rd.rds")
imp_binary_group_sametomb <- readRDS("unofficial work/models save/imp_binary_group_sametomb.rds")



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

##save the plots 
ggsave("Plots/robust_age_sep.png",   p_age_separate,  width = 8, height = 6, dpi = 400)
ggsave("Plots/robust_sex_sep.png",   p_sex_separate,  width = 8, height = 6, dpi = 400)
ggsave("Plots/robust_tomb_sep.png",  p_tomb_separate, width = 8, height = 4, dpi = 400)
ggsave("Plots/robust_intercept_sep.png", p_int_separate, width = 10, height = 3.5, dpi = 400)
ggsave("Plots/robust_random_sd_sep.png", p_re_separate, width = 10, height = 4.5, dpi = 400)


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

##save the plots 
ggsave("Plots/robust_age_gr.png",   p_age_group,  width = 8, height = 6, dpi = 400)
ggsave("Plots/robust_sex_gr.png",   p_sex_group,  width = 8, height = 6, dpi = 400)
ggsave("Plots/robust_tomb_gr.png",  p_tomb_group, width = 8, height = 4, dpi = 400)
ggsave("Plots/robust_intercept_gr.png", p_int_group, width = 10, height = 3.5, dpi = 400)
ggsave("Plots/robust_random_sd_gr.png", p_re_group, width = 10, height = 4.5, dpi = 400)

