# This script checks RE normality assumptions and returns Q-Q plots.
# Input: a fitted model object
# Output: Q-Q plots

# Load project environment
source("environment_setup.R")

imp_binary_separate_nobone <- readRDS("Models/Saved models/imp_binary_separate_nobone.rds")
imp_binary_group_nobone <- readRDS("Models/Saved models/imp_binary_group_nobone.rds")

# This script is to check the assumption of RE normality 

plot_qq_both_effects <- function(model_object, subtitle = "Separated tomb") {
  
  # Extract random effects
  re_list <- ranef(model_object)
  
  plot_list <- list()
  
  for (group_name in names(re_list)) {
    re_df <- as.data.frame(re_list[[group_name]])
    
    col_idx <- grep("Estimate", names(re_df))[1]
    names(re_df)[col_idx] <- "Estimate"
    
    re_df$z_score <- (re_df$Estimate - mean(re_df$Estimate)) / sd(re_df$Estimate)
    
    # Q-Q Plot
    p <- ggplot(re_df, aes(sample = z_score)) +
      stat_qq(size = 2, color = "skyblue4", alpha = 0.8) +
      stat_qq_line(color = "red", linewidth = 1) +
      labs(title = paste("Q-Q Plot:", group_name),
           subtitle = subtitle,
           x = "Theoretical Quantiles",
           y = "Sample Quantiles (Standardized)") +
      theme_minimal() +
      theme(plot.title = element_text(size = 11, face = "bold"))
    
    plot_list[[group_name]] <- p
  }
  
  grid.arrange(grobs = plot_list, ncol = 2)
  
  return(invisible(plot_list))
}

plot_qq_both_effects(imp_binary_separate_nobone, "Separated tombs")
plot_qq_both_effects(imp_binary_group_nobone, "Grouped tombs")
