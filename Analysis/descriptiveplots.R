# Load environment and data
source(here::here("environmentSetUp.R"))
source(here::here("analysis", "pairsdata.R"))

######################################################
### Absolute relation of MNI and analysed individuals
######################################################

# Tombs to exclude from plot (keep in dataset but hide in plot)
tombs_to_exclude <- c(
  "Elateia T 46,56,62 (MNI 300)",
  "Elateia T 46,56,62 (MNI 350)",
  "Elateia T 46,56,62 (MNI 400)"
)

# Filter tombs for plotting
df_long_filtered <- df_long %>% 
  filter(!Tomb %in% tombs_to_exclude)

# Plot stacked bar chart of MNI vs analysed individuals per tomb
plot_mni <- ggplot(df_long_filtered, aes(x = Tomb, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  scale_fill_manual(values = c(
    "MNI" = "#1f78b4",
    "analysed individuals" = "#ff7f00"
  )) +
  scale_x_discrete(labels = c(
    "Elateia T 46,56,62 (MNI 250)" = "Elateia T 46,56,62"  # clean tomb label
  )) +
  theme_minimal(base_size = 12)

plot_mni

#################################################
#### Relative number of skeletal elements by Sex
#################################################

# Count skeletal elements per tomb and sex
# Ensure all combinations exist, fill missing with 0
# Compute relative proportion per tomb
heat_sex <- individual_metadata %>%
  count(tomb, sex, skeletal_element) %>%                        
  complete(tomb, sex, skeletal_element, fill = list(n = 0)) %>% 
  group_by(tomb, sex) %>%
  mutate(prop = n / sum(n)) %>%                                 
  ungroup()

# Plot heatmap of relative proportions by sex
plot_sex <- ggplot(heat_sex, aes(x = skeletal_element, y = tomb, fill = prop)) +
  geom_tile(color = "black") +
  facet_wrap(~ sex, labeller = labeller(
    sex = c("XX" = "Female", "XY" = "Male")
  )) +
  scale_fill_gradient(low = "white", high = "darkblue", labels = scales::percent_format(1)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()) +
  labs(x = "Skeletal Element", y = "Tomb")


################################################
### Relative number of skeletal elements by Age
################################################

# Count skeletal elements per tomb and age group
# Ensure all combinations exist, fill missing with 0
# Compute relative proportion per tomb
heat_age <- individual_metadata %>%
  count(tomb, age_estimation, skeletal_element) %>%
  complete(tomb, age_estimation, skeletal_element, fill = list(n = 0)) %>%
  group_by(tomb, age_estimation) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Plot heatmap of relative proportions by age
plot_age <- ggplot(heat_age, aes(x = skeletal_element, y = tomb, fill = prop)) +
  geom_tile(color = "black") +
  facet_wrap(~ age_estimation, labeller = labeller(
    age_estimation = c("adult or undefined" = "Adult or Undefined",
                       "subadult" = "Subadult")
  )) +
  scale_fill_gradient(low = "white", high = "darkblue", labels = scales::percent_format(1)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()) +
  labs(x = "Skeletal Element", y = "Tomb")

###########################################################
### Absolute and Relative Number of Analysed Pairs Per Tomb
###########################################################

# Step 1: Prepare pairwise counts including cross-tomb pairs
tomb_pairwise <- kinship_result %>%
  pivot_longer(
    cols = c(tomb1, tomb2),
    names_to = "tomb_type",
    values_to = "tomb"
  ) %>%
  mutate(
    Rel_group = case_when(
      related_binary == 1 ~ "Related",
      related_binary == 0 ~ "Unrelated",
      is.na(related_binary) ~ "Uncertain"
    )
  ) %>%
  group_by(tomb, Rel_group) %>%
  summarise(num_pairs = n(), .groups = "drop") %>%
  group_by(tomb) %>%
  mutate(prop = num_pairs / sum(num_pairs)) %>%
  ungroup()

# Step 2: Relative proportions plot
p_relative <- ggplot(tomb_pairwise, aes(x = Rel_group, y = prop, fill = Rel_group)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ tomb) +
  scale_y_continuous(labels = percent_format()) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Proportion of analysed pairs") +
  xlab("Relationship")

# Step 3: Absolute counts plot
p_absolute <- ggplot(tomb_pairwise, aes(x = Rel_group, y = num_pairs, fill = Rel_group)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ tomb) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Number of analysed pairs") +
  xlab("Relationship")

# Step 4: Print both plots
p_relative
p_absolute

#######################
# Save plots
#######################

if (!file.exists("Plots/plot_mni.png")) {
  ggsave("Plots/plot_mni.png", plot = plot_mni, width = 8, height = 5, dpi = 300)
}

if (!file.exists("Plots/heatmap_sex.png")) {
  ggsave("Plots/heatmap_sex.png", plot = plot_sex, width = 8, height = 6, dpi = 300)
}

if (!file.exists("Plots/heatmap_age.png")) {
  ggsave("Plots/heatmap_age.png", plot = plot_age, width = 8, height = 6, dpi = 300)
}

if (!file.exists("Plots/plot_pairs_rel.png")) {
  ggsave("Plots/plot_pairs_rel.png", plot = p_relative, width = 8, height = 5, dpi = 300)
}

if (!file.exists("Plots/plot_pairs_abs.png")) {
  ggsave("Plots/plot_pairs_abs.png", plot = p_absolute, width = 8, height = 5, dpi = 300)
}
