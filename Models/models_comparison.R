source("environmentSetUp.R")

# Load models ----
## Separated tomb version ----
imp_binary_separate_cat <- readRDS("unofficial work/models save/imp_binary_separate_cat.rds")
all_binary_separate_cat <- readRDS("unofficial work/models save/all_binary_separate_cat.rds")
imp_binary_separate_nobone <- readRDS("unofficial work/models save/imp_binary_separate_nobone.rds")
imp_binary_separate_newtombvar <- readRDS("unofficial work/models save/imp_binary_separate_newtombvar.rds")
imp_binary_separate_withintomb <- readRDS("unofficial work/models save/imp_binary_separate_withintomb.rds")
imp_binary_separate_unrel3rd <- readRDS("unofficial work/models save/imp_binary_separate_unrel3rd.rds")
imp_binary_separate_sametomb <- readRDS("unofficial work/models save/imp_binary_separate_sametomb.rds")
imp_binary_separate_nobone_prior <- readRDS("unofficial work/models save/imp_binary_separate_nobone_prior.rds")
imp_binary_separate_withintomb_unrel3rd <- readRDS("unofficial work/models save/imp_binary_separate_withintomb_unrel3rd.rds")
all_binary_separate_nobone <-  readRDS("unofficial work/models save/all_binary_separate_nobone.rds")

## Grouped tomb version ----
imp_binary_group_cat <- readRDS("unofficial work/models save/imp_binary_group_cat.rds")
all_binary_group_cat <- readRDS("unofficial work/models save/all_binary_group_cat.rds")
imp_binary_group_nobone <- readRDS("unofficial work/models save/imp_binary_group_nobone.rds")
imp_binary_group_newtombvar <- readRDS("unofficial work/models save/imp_binary_group_newtombvar.rds")
imp_binary_group_withintomb <- readRDS("unofficial work/models save/imp_binary_group_withintomb.rds")
imp_binary_group_unrel3rd <- readRDS("unofficial work/models save/imp_binary_group_unrel3rd.rds")
imp_binary_group_sametomb <- readRDS("unofficial work/models save/imp_binary_group_sametomb.rds")
imp_binary_group_nobone_prior <- readRDS("unofficial work/models save/imp_binary_group_nobone_prior.rds")
imp_binary_group_withintomb_unrel3rd <- readRDS("unofficial work/models save/imp_binary_group_withintomb_unrel3rd.rds")
all_binary_group_nobone <- readRDS("unofficial work/models save/all_binary_group_nobone.rds")

# Define compare function ----
compare_model <- function(models, ci_width = 0.95, digits = 2, output = c("gt", "data"),
                          title = "Model comparison"
) {
  output <- match.arg(output)
  
  extract_one <- function(fit, model_name) {
    tidy_draws(fit) |>
      pivot_longer(
        cols = matches("^(b_|sd_)"),
        names_to = "term",
        values_to = "value"
      ) |>
      group_by(term) |>
      median_qi(value, .width = ci_width) |>
      transmute(
        term,
        model = model_name,
        estimate = value,
        ci_low = .lower,
        ci_high = .upper
      )
  }
  
  long <- imap_dfr(models, extract_one)
  
  wide <- long |>
    mutate(
      term = case_when(
        term == "b_Intercept" ~ "Intercept",
        grepl("^b_", term) ~ gsub("^b_", "", term),
        grepl("^sd_", term) ~ gsub("^sd_", "SD(RE): ", term),
        TRUE ~ term
      ),
      CI = sprintf(
        paste0("[%.", digits, "f, %.", digits, "f]"),
        ci_low, ci_high
      ),
      row_order = case_when(
        term == "Intercept" ~ 1,
        grepl("^SD\\(RE\\):", term) ~ 3,
        TRUE ~ 2
      )
    ) |>
    select(term, model, estimate, CI, row_order) |>
    pivot_wider(
      names_from = model,
      values_from = c(estimate, CI),
      names_glue = "{model}_{.value}"
    )
  
  col_order <- c(
    "term",
    unlist(
      lapply(
        names(models),
        function(m) c(paste0(m, "_estimate"), paste0(m, "_CI"))
      )
    )
  )
  
  tab <- wide |>
    arrange(row_order, term) |>
    select(all_of(col_order))
  
  # ---- default: return gt table ----
  if (output == "gt") {
    gt_tab <- gt(tab) |>
      tab_header(
        title = title
      ) |>
      cols_align(align = "left", term) |>
      cols_align(align = "center", -term)
    
    # automatic column spanners per model
    for (m in names(models)) {
      gt_tab <- gt_tab |>
        tab_spanner(
          label = m,
          columns = c(paste0(m, "_estimate"), paste0(m, "_CI"))
        )
    }
    
    return(gt_tab)
  }
  
  return(tab)
}

# Compare model ----
## Compare between 15k threshold and original threhold ----
### Separated version ----
stdthreshold_vs_15k <- list("15k threshold" = imp_binary_separate_cat,
                            "standard threshold" = all_binary_separate_cat)
compare_model(stdthreshold_vs_15k)


## Robust check main model ----
### Separated version ----
robustcheck_models_separate <- list(
  "Main model" = imp_binary_separate_nobone,
  "Only same_tomb" = imp_binary_separate_sametomb,
  "Within tomb" = imp_binary_separate_withintomb,
  "3rd as unrelated" = imp_binary_separate_unrel3rd
)

compare_model(robustcheck_models_separate, title = "Model Robust Check (Separated)")

### Grouped version ----
robustcheck_models_group <- list(
  "Main model" = imp_binary_group_nobone,
  "Only same_tomb" = imp_binary_group_sametomb,
  "Within tomb" = imp_binary_group_withintomb,
  "3rd as unrelated" = imp_binary_group_unrel3rd
)

compare_model(robustcheck_models_group, title = "Model Robust Check (Grouped)")









