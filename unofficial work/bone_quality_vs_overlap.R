# Load environment and data
source(here::here("environmentSetUp.R"))
source(here::here("pairsdata.R"))

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



