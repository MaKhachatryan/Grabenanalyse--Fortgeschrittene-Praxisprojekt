source("environmentSetUp.R")
# Check assumption of normality RE

library(dplyr)
library(ggplot2)

library(brms)
library(ggplot2)
library(dplyr)

library(ggplot2)
library(dplyr)

## -------------------------------
## 1. Extract random effects
## -------------------------------

r_effs <- ranef(imp_binary_separate_nobone)

# Individual MM
re_ind <- r_effs$mmindividual1_idindividual2_id
random_ind <- data.frame(
  id = rownames(re_ind[, , "Intercept"]),
  Intercept = re_ind[, "Estimate", "Intercept"]
)

# Tomb MM
re_tomb <- r_effs$mmtomb1tomb2
random_tomb <- data.frame(
  tomb = rownames(re_tomb[, , "Intercept"]),
  Intercept = re_tomb[, "Estimate", "Intercept"]
)

## -------------------------------
## 2. Extract SDs from VarCorr
## -------------------------------

vc <- VarCorr(imp_binary_separate_nobone)

sd_ind <- as.numeric(vc$mmindividual1_idindividual2_id$sd[1])
sd_tomb <- as.numeric(vc$mmtomb1tomb2$sd[1])

## -------------------------------
## 3. Bin random effects for bar charts
## -------------------------------

# Individual MM
bins <- seq(
  floor(min(random_ind$Intercept) - 0.5),
  ceiling(max(random_ind$Intercept) + 0.5),
  by = 0.2
)
random_ind$bin <- cut(random_ind$Intercept, breaks = bins, include.lowest = TRUE)

counts_ind <- random_ind %>%
  group_by(bin) %>%
  summarise(count = n())

# Tomb MM
bins_tomb <- seq(
  floor(min(random_tomb$Intercept) - 0.5),
  ceiling(max(random_tomb$Intercept) + 0.5),
  by = 0.2
)
random_tomb$bin <- cut(random_tomb$Intercept, breaks = bins_tomb, include.lowest = TRUE)

counts_tomb <- random_tomb %>%
  group_by(bin) %>%
  summarise(count = n())

## -------------------------------
## 4. Normal distribution for overlay
## -------------------------------

# Individuals
xseq_ind <- seq(min(random_ind$Intercept) - 0.5,
                max(random_ind$Intercept) + 0.5,
                by = 0.01)
normal_ind <- data.frame(
  x = xseq_ind,
  y = dnorm(xseq_ind, mean = 0, sd = sd_ind) * nrow(random_ind) * 0.2 # scale to match bar height
)

# Tombs
xseq_tomb <- seq(min(random_tomb$Intercept) - 0.5,
                 max(random_tomb$Intercept) + 0.5,
                 by = 0.01)
normal_tomb <- data.frame(
  x = xseq_tomb,
  y = dnorm(xseq_tomb, mean = 0, sd = sd_tomb) * nrow(random_tomb) * 0.2
)

## -------------------------------
## 5. Plot bar charts with normal overlay
## -------------------------------

# Individual MM
p_ind <- ggplot() +
  geom_bar(data = counts_ind, aes(x = bin, y = count), stat = "identity", fill = "grey80", color = "black") +
  geom_line(data = normal_ind, aes(x = x, y = y), color = "red", size = 1) +
  theme_minimal() +
  labs(x = "Random Effect (Individual MM)", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Tomb MM
p_tomb <- ggplot() +
  geom_bar(data = counts_tomb, aes(x = bin, y = count), stat = "identity", fill = "grey80", color = "black") +
  geom_line(data = normal_tomb, aes(x = x, y = y), color = "red", size = 1) +
  theme_minimal() +
  labs(x = "Random Effect (Tomb MM)", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display plots
p_ind
p_tomb

