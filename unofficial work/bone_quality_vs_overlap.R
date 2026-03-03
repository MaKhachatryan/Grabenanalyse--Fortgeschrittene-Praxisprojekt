# Load environment and data
suppressMessages(
  suppressWarnings(
    source(here::here("environment_setup.R"))
  )
)

suppressMessages(
  suppressWarnings(
    source(here::here("Analysis", "pairs_data.R"))
  )
)

# --------------------------------------------
# Plot overlap_nsnps vs bone_quality
# --------------------------------------------

# Boxplot
ggplot(pairs_ind, aes(x = bone_quality, y = overlap_nsnps, fill = bone_quality)) +
  geom_boxplot(alpha = 0.6) +
  labs(
    title = "SNP Overlap by Bone Quality",
    x = "Bone Quality",
    y = "Number of Overlapping SNPs"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


# Jitter plot
ggplot(pairs_ind, aes(x = bone_quality, y = overlap_nsnps, color = bone_quality)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 1.5) +
  labs(
    title = "SNP Overlap by Bone Quality",
    x = "Bone Quality",
    y = "Number of Overlapping SNPs"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none"
  )


# --------------------------------------------
# Plot overlap_nsnps vs skeletal_element_pair
# --------------------------------------------

# Boxplot
ggplot(pairs_ind, aes(x = skeletal_element_pair, y = overlap_nsnps, fill = skeletal_element_pair)) +
  geom_boxplot(alpha = 0.6) +
  labs(
    title = "SNP Overlap by Skeletal Element Pair",
    x = "Skeletal Element Pair",
    y = "Number of Overlapping SNPs"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1) # rotate labels for readability
  )


# Jitter plot
ggplot(pairs_ind, aes(x = skeletal_element_pair, y = overlap_nsnps, color = skeletal_element_pair)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 1.5) +
  labs(
    title = "SNP Overlap by Skeletal Element Pair",
    x = "Skeletal Element Pair",
    y = "Number of Overlapping SNPs"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)  # rotate labels for readability
  )


#### Test

# --------------------------------------------
# The Kruskal–Wallis test
# --------------------------------------------


# “Do the distributions of a numeric variable 
# differ across two or more groups?”

kruskal.test(overlap_nsnps ~ bone_quality, data = pairs_ind)

### output
# data:  overlap_nsnps by bone_quality
# Kruskal-Wallis chi-squared = 2694.4, df = 2, p-value < 2.2e-16
 

### explanation
# In Kruskal-Wallis, the chi-squared value measures
#
# how different the ranks of overlap_nsnps are across the bone_quality groups.
#
# Larger values → greater separation between groups.
# 
#  2694.4 is huge, so the distributions are very different.

# df = number of groups − 1

# p-value < 2.2e-16
# Extremely small.
# 
# Statistically, this rejects the null hypothesis
#
# that the distributions are identical.
#
# Important: With 27,730 rows, almost any difference
#
# becomes “statistically significant”.
#
# So we don’t focus on the p-value.
#
# Focus on effect size and patterns.


### conclusion
# the boxplot and Kruskal-Wallis test show that
#
# overlap_nsnps and bone_quality are associated
 


# --------------------------------------------
# Spearman correlation
# --------------------------------------------

# 1 → worst quality, 3 → best quality
pairs_ind_ord <- pairs_ind |>
  mutate(
    bone_quality_ord = case_when(
      bone_quality == "low-low"  ~ 1,
      bone_quality == "high-low" ~ 2,
      bone_quality == "high-high"~ 3
    )
  )

# corr
cor.test(
  pairs_ind_ord$bone_quality_ord,
  pairs_ind_ord$overlap_nsnps,
  method = "spearman"
)

# output
#
# rho = 0.3103595
#
# p-value < 2.2e-16


# explanation
# ρ ≈ 0.31 → moderate positive association
# 
# Higher bone quality tends to correspond to higher overlap_nsnps
# 
# p-value < 2.2e-16 → extremely small
# 
# As expected with 27,730 rows
# 
# Significance is guaranteed with huge sample; focus on ρ, not p


# conclusion
# Moderate positive correlation supports
#
# the idea that overlap_nsnps and bone_quality are associated













