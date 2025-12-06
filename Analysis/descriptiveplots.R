# Load environment and data
source(here::here("environmentSetUp.R"))
source(here::here("pairsdata.R"))

### Absolute relation of MNI and analysed individuals

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


#### Relative number of skeletal elements by Sex

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


### Relative number of skeletal elements by Age

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

# Save plots
ggsave("Plots/plot_mni.png", plot = plot_mni, width = 8, height = 5, dpi = 300)
ggsave("Plots/heatmap_sex.png", plot = plot_sex, width = 8, height = 6, dpi = 300)
ggsave("Plots/heatmap_age.png", plot = plot_age, width = 8, height = 6, dpi = 300)
