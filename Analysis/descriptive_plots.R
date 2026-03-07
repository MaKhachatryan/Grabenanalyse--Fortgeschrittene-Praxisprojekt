# This script creates descriptive plots based on the cleaned project datasets.
# It generates visualizations for tomb characteristics, skeletal element
# distributions, pairwise relatedness counts, SNP overlap quality measures,
# and p0 distributions, and saves all plots as .png files in the Plots folder.

# Load project environment and processed data
source("environment_setup.R")
source("Data/pairs_data.R")

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


# Step 2: Relative proportions plot (Modified)
p_relative_rel <- ggplot(tomb_pairwise, aes(x = Rel_group, y = prop, fill = Rel_group)) +
  geom_bar(stat = "identity") +
  # Add percentage labels
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1)), 
            vjust = -0.5,     # Moves text slightly above the bar
            color = "black",  # Basic black writing
            size = 3) +       # Adjust text size if needed
  facet_wrap(~ tomb) +
  # 'expand' adds 15% padding at the top so labels aren't cut off
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.15))) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Proportion of analysed pairs") +
  xlab("Relationship")


# Step 3: Absolute counts plot (Modified)
p_absolute_rel <- ggplot(tomb_pairwise, aes(x = Rel_group, y = num_pairs, fill = Rel_group)) +
  geom_bar(stat = "identity") +
  # Add count labels
  geom_text(aes(label = num_pairs), 
            vjust = -0.5,     # Moves text slightly above the bar
            color = "black", 
            size = 3) +
  facet_wrap(~ tomb) +
  # 'expand' adds 15% padding at the top
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Number of analysed pairs") +
  xlab("Relationship")





##############################################################
### Number and proportion of within-tomb pairs by relationship
##############################################################


# --------------------------
# Step 1: Filter within-tomb pairs
# --------------------------
within_tomb <- kinship_result %>%
  filter(tomb1 == tomb2) %>%
  mutate(
    Rel_group = case_when(
      related_binary == 1 ~ "Related",
      related_binary == 0 ~ "Unrelated",
      is.na(related_binary) ~ "Uncertain"
    ),
    # Set order of bars
    Rel_group = factor(Rel_group, levels = c("Related", "Unrelated", "Uncertain"))
  )

# --------------------------
# Step 2: Summarise counts per tomb and relationship
# --------------------------
tomb_within_summary <- within_tomb %>%
  group_by(tomb = tomb1, Rel_group) %>%
  summarise(num_pairs = n(), .groups = "drop") %>%
  group_by(tomb) %>%
  mutate(prop = num_pairs / sum(num_pairs)) %>%
  ungroup()

# --------------------------
# Step 3: Absolute counts plot
# --------------------------
p_within_abs <- ggplot(tomb_within_summary, aes(x = Rel_group, y = num_pairs, fill = Rel_group)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = num_pairs), vjust = -0.5, size = 3) +
  facet_wrap(~ tomb) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +  # extra space for labels
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ylab("Number Of Within-tomb Pairs") +
  xlab("Relationship")

# --------------------------
# Step 4: Relative proportions plot
# --------------------------
p_within_rel <- ggplot(tomb_within_summary, aes(x = Rel_group, y = prop, fill = Rel_group)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1)), vjust = -0.5, size = 3) +
  facet_wrap(~ tomb) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.2))) +  # same top padding
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ylab("Proportion Of Within-tomb Pairs") +
  xlab("Relationship")




###########################################################
### Plot Error vs. SNP Count
### Plot CV vs SNP Count
###########################################################

# 1. Calculate the CV for each pair
kinship_p0_cv <- og_kinship_result %>%
  mutate(CV = nonnormalized_p0_serr / nonnormalized_p0)

custom_breaks <- c(15000, 200000, 400000, 600000, 800000)


# 2. Plot CV vs SNP Count (Relative Noise)
plot_p0_overlap <- ggplot(kinship_p0_cv, aes(x = overlap_nsnps, y = CV)) +
  geom_point(alpha = 0.3, color = "darkgreen") +
  # Adds a bold baseline at 0
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  geom_vline(xintercept = 15000, linetype = "dashed", color = "red") +
  # Log scale spreads out the critical 0-50k SNP range
  scale_x_log10(labels = label_comma(), 
                breaks = c(100, 1000, 15000, 100000, 800000)) +
  ylim(0, 0.5) +
  labs(x = "Number of Overlapping SNPs (Log Scale)",
       y = "Coefficient of Variation (Relative Error of P0)") +
  theme_minimal() 




###########################################################
### P0 vs degree bar plot
###########################################################


# 1. Filter data to include only the degrees of interest
df_p0_plot <- og_kinship_result %>%
  filter(!rel %in% c("IdenticalTwins/SameIndividual", 
                     "Unrelated/Consistent with Third Degree"))

# 2. Create a Binned Bar Plot
# We divide p0_mean into 0.05 increments to show distribution clearly
plot_p0_reldegree <- ggplot(df_p0_plot, aes(x = p0_mean, fill = rel)) +
  geom_histogram(binwidth = 0.01, color = "white") +
  # Separate plots for each degree to handle the scale difference
  facet_wrap(~rel, scales = "free_y") + 
  scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    # This makes the degree names bold and larger
    strip.text = element_text(face = "bold", size = 12), 
    axis.title = element_text(face = "bold")
  ) + 
  labs(
    x = "p0 mean (Normalized Mismatch Rate)",
    y = "Count of Individual Pairs"
  )



#######################################################
### Network plots of genetic relatedness
#######################################################

plot_dir <- here::here("Plots")

set.seed(123)

# -------------------------
# Tombs to plot
# -------------------------
tombs <- c(
  "Amfissa tholos",
  "Elateia T31",
  "Elateia T36",
  "Elateia T46",
  "Elateia T50",
  "Elateia T56",
  "Elateia T62",
  "Elateia T67"
)

# -------------------------
# Filename map
# -------------------------
file_map <- c(
  "Amfissa tholos" = "Amfissa_Tholos_network.png",
  "Elateia T31"    = "Elateia_T31_network.png",
  "Elateia T36"    = "Elateia_T36_network.png",
  "Elateia T46"    = "Elateia_T46_network.png",
  "Elateia T50"    = "Elateia_T50_network.png",
  "Elateia T56"    = "Elateia_T56_network.png",
  "Elateia T62"    = "Elateia_T62_network.png",
  "Elateia T67"    = "Elateia_T67_network.png"
)

# -------------------------
# Individual tomb networks
# -------------------------
for (t in tombs) {
  
  edges <- pairs_ind_og |>
    dplyr::filter(
      rel %in% c("First Degree", "Second Degree", "Third Degree"),
      tomb1 == t | tomb2 == t
    ) |>
    dplyr::mutate(
      tomb_relation = ifelse(tomb1 == tomb2,
                             "Within tomb",
                             "Cross tomb"),
      rel = factor(
        rel,
        levels = c("First Degree",
                   "Second Degree",
                   "Third Degree")
      )
    )
  
  if (nrow(edges) == 0) {
    message("No relations for ", t, " — skipped")
    next
  }
  
  nodes <- dplyr::bind_rows(
    edges |>
      dplyr::select(individual_id = individual1_id,
                    tomb = tomb1,
                    sex = sex1),
    edges |>
      dplyr::select(individual_id = individual2_id,
                    tomb = tomb2,
                    sex = sex2)
  ) |>
    dplyr::distinct()
  
  g <- igraph::graph_from_data_frame(
    d = edges |>
      dplyr::select(from = individual1_id,
                    to   = individual2_id,
                    rel,
                    tomb_relation),
    vertices = nodes |>
      dplyr::select(name = individual_id,
                    tomb,
                    sex),
    directed = FALSE
  )
  
  p <- ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(
      aes(width = rel, color = tomb_relation),
      alpha = 0.8
    ) +
    ggraph::geom_node_point(
      aes(shape = sex, color = tomb),
      size = 4
    ) +
    ggraph::scale_edge_width_manual(
      values = c(
        "First Degree"  = 2,
        "Second Degree" = 1.2,
        "Third Degree"  = 0.6
      ),
      name = "Kinship Degree"
    ) +
    ggraph::scale_edge_color_manual(
      values = c(
        "Within tomb" = "black",
        "Cross tomb"  = "red"
      ),
      name = "Tomb Relation"
    ) +
    ggplot2::scale_color_discrete(name = "Tomb") +  
    ggplot2::scale_shape_manual(
      name = "Sex",
      values = c("XX" = 16, "XY" = 17),
      labels = c("XX" = "Female", "XY" = "Male")
    ) +
    ggplot2::theme_void() +
    ggplot2::labs(title = t)
  
  out_file <- file.path(plot_dir, file_map[t])
  
  if (!file.exists(out_file)) {
    ggplot2::ggsave(
      filename = out_file,
      plot = p,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
}

# -------------------------
# All tombs together
# -------------------------

all_edges <- pairs_ind_og |>
  dplyr::filter(
    rel %in% c("First Degree",
               "Second Degree",
               "Third Degree")
  ) |>
  dplyr::mutate(
    tomb_relation = ifelse(tomb1 == tomb2,
                           "Within tomb",
                           "Cross tomb"),
    rel = factor(rel,
                 levels = c("First Degree",
                            "Second Degree",
                            "Third Degree"))
  )

all_nodes <- dplyr::bind_rows(
  all_edges |>
    dplyr::select(individual_id = individual1_id,
                  tomb = tomb1,
                  sex = sex1),
  all_edges |>
    dplyr::select(individual_id = individual2_id,
                  tomb = tomb2,
                  sex = sex2)
) |>
  dplyr::distinct()

g_all <- igraph::graph_from_data_frame(
  d = all_edges |>
    dplyr::select(from = individual1_id,
                  to   = individual2_id,
                  rel,
                  tomb_relation),
  vertices = all_nodes |>
    dplyr::select(name = individual_id,
                  tomb,
                  sex),
  directed = FALSE
)

p_all <- ggraph::ggraph(g_all, layout = "fr") +
  ggraph::geom_edge_link(
    aes(width = rel, color = tomb_relation),
    alpha = 0.8
  ) +
  ggraph::geom_node_point(
    aes(shape = sex, color = tomb),
    size = 4
  ) +
  ggraph::scale_edge_width_manual(
    values = c(
      "First Degree"  = 2,
      "Second Degree" = 1.2,
      "Third Degree"  = 0.6
    )
  ) +
  ggraph::scale_edge_color_manual(
    values = c(
      "Within tomb" = "black",
      "Cross tomb"  = "red"
    )
  ) +
  ggplot2::scale_shape_manual(
    values = c("XX" = 16, "XY" = 17)
  ) +
  ggplot2::theme_void() +
  ggplot2::labs(title = "All Tombs Together")

out_file <- file.path(plot_dir, "All_Tombs_network.png")

if (!file.exists(out_file)) {
  ggplot2::ggsave(
    filename = out_file,
    plot = p_all,
    width = 8,
    height = 6,
    dpi = 300
  )
}


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
  ggsave("Plots/plot_pairs_rel.png", plot = p_relative_rel, width = 8, height = 5, dpi = 300)
}

if (!file.exists("Plots/plot_pairs_abs.png")) {
  ggsave("Plots/plot_pairs_abs.png", plot = p_absolute_rel, width = 8, height = 5, dpi = 300)
}

if (!file.exists("Plots/plot_p0_overlap.png")) {
  ggsave("Plots/plot_p0_overlap.png", plot = plot_p0_overlap, width = 8, height = 5, dpi = 300)
}

if (!file.exists("Plots/plot_p0_reldegree.png")) {
  ggsave("Plots/plot_p0_reldegree.png", plot = plot_p0_reldegree, width = 8, height = 5, dpi = 300)
}


if (!file.exists("Plots/p_within_abs.png")) {
  ggsave("Plots/p_within_abs.png", plot = p_within_abs, width = 8, height = 5, dpi = 300)
}

if (!file.exists("Plots/p_within_rel.png")) {
  ggsave("Plots/p_within_rel.png", plot = p_within_rel, width = 8, height = 5, dpi = 300)
}
