# Load environment and data
source(here::here("environmentSetUp.R"))
source(here::here("analysis", "pairsdata.R"))


###################### Amfissa Tholos
set.seed(123)

# -------------------------
# 1. Filter edges for Amfissa tholos
# -------------------------
amfissa_edges <- pairs_ind_og |>
  filter(tomb1 == "Amfissa tholos" | tomb2 == "Amfissa tholos") |>
  filter(rel %in% c("First Degree", "Second Degree", "Third Degree")) |>
  mutate(
    tomb_relation = ifelse(tomb1 == tomb2, "Within tomb", "Cross tomb"),
    rel = factor(rel, levels = c("First Degree", "Second Degree", "Third Degree"))
  )

# -------------------------
# 2. Create node table
# -------------------------
amfissa_nodes <- bind_rows(
  amfissa_edges |>
    select(individual_id = individual1_id,
           tomb = tomb1,
           sex = sex1),
  amfissa_edges |>
    select(individual_id = individual2_id,
           tomb = tomb2,
           sex = sex2)
) |>
  distinct()

# -------------------------
# 3. Build graph
# -------------------------
g_amfissa <- graph_from_data_frame(
  d = amfissa_edges |>
    select(from = individual1_id,
           to   = individual2_id,
           rel,
           tomb_relation),
  vertices = amfissa_nodes |>
    select(name = individual_id, tomb, sex),
  directed = FALSE
)

# -------------------------
# 4. Plot network with legend
# -------------------------
p_amfissa <- ggraph(g_amfissa, layout = "fr") +
  
  geom_edge_link(
    aes(width = rel, color = tomb_relation),
    alpha = 0.8
  ) +
  
  geom_node_point(
    aes(shape = sex, color = tomb),
    size = 4
  ) +
  
  scale_edge_width_manual(
    values = c(
      "First Degree"  = 2,
      "Second Degree" = 1.2,
      "Third Degree"  = 0.6
    ),
    name = "Kinship Degree"
  ) +
  
  scale_edge_color_manual(
    values = c(
      "Within tomb" = "black",
      "Cross tomb"  = "red"
    ),
    name = "Tomb Relation"
  ) +
  
  scale_color_discrete(name = "Tomb") +
  
  scale_shape_manual(
    name = "Sex",
    values = c("XX" = 16, "XY" = 17), 
    labels = c("XX" = "Female", "XY" = "Male")
  ) +
  
  theme_void() +
  labs(title = "Amfissa Tholos")


##################### Elateia 62

set.seed(123)

# -------------------------
# 1. Filter edges for Elateia T62 only
# -------------------------
elateia62_edges <- pairs_ind_og |>
  filter(rel %in% c("First Degree", "Second Degree", "Third Degree")) |>
  filter(tomb1 == "Elateia T62" | tomb2 == "Elateia T62") |>
  mutate(
    tomb_relation = ifelse(tomb1 == tomb2, "Within tomb", "Cross tomb"),
    rel = factor(rel, levels = c("First Degree", "Second Degree", "Third Degree"))
  )

# -------------------------
# 2. Nodes: only individuals in these edges
# -------------------------
elateia62_nodes <- bind_rows(
  elateia62_edges |>
    select(individual_id = individual1_id,
           tomb = tomb1,
           sex = sex1),
  elateia62_edges |>
    select(individual_id = individual2_id,
           tomb = tomb2,
           sex = sex2)
) |>
  distinct()

# -------------------------
# 3. Build graph
# -------------------------
g_elateia62 <- graph_from_data_frame(
  d = elateia62_edges |>
    select(from = individual1_id,
           to   = individual2_id,
           rel,
           tomb_relation),
  vertices = elateia62_nodes |>
    select(name = individual_id, tomb, sex),
  directed = FALSE
)

# -------------------------
# 4. Plot network with legend
# -------------------------
p_elateia62 <- ggraph(g_elateia62, layout = "fr") +
  
  geom_edge_link(
    aes(width = rel, color = tomb_relation),
    alpha = 0.8
  ) +
  
  geom_node_point(
    aes(shape = sex, color = tomb),
    size = 4
  ) +
  
  scale_edge_width_manual(
    values = c(
      "First Degree"  = 2,
      "Second Degree" = 1.2,
      "Third Degree"  = 0.6
    ),
    name = "Kinship Degree"
  ) +
  
  scale_edge_color_manual(
    values = c(
      "Within tomb" = "black",
      "Cross tomb"  = "red"
    ),
    name = "Tomb Relation"
  ) +
  
  scale_color_discrete(name = "Tomb") +
  
  scale_shape_manual(
    name = "Sex",
    values = c("XX" = 16, "XY" = 17), 
    labels = c("XX" = "Female", "XY" = "Male")
  ) +
  
  theme_void() +
  labs(title = "Elateia T62")

################ All Tombs Together

set.seed(123)

# -------------------------
# 1. Filter edges: all tombs
# -------------------------
all_edges <- pairs_ind_og |>
  filter(rel %in% c("First Degree", "Second Degree", "Third Degree")) |>
  mutate(
    tomb_relation = ifelse(tomb1 == tomb2, "Within tomb", "Cross tomb"),
    rel = factor(rel, levels = c("First Degree", "Second Degree", "Third Degree"))
  )

# -------------------------
# 2. Nodes: all individuals in these edges
# -------------------------
all_nodes <- bind_rows(
  all_edges |>
    select(individual_id = individual1_id,
           tomb = tomb1,
           sex = sex1),
  all_edges |>
    select(individual_id = individual2_id,
           tomb = tomb2,
           sex = sex2)
) |>
  distinct()

# -------------------------
# 3. Build graph
# -------------------------
g_all <- graph_from_data_frame(
  d = all_edges |>
    select(from = individual1_id,
           to   = individual2_id,
           rel,
           tomb_relation),
  vertices = all_nodes |>
    select(name = individual_id, tomb, sex),
  directed = FALSE
)

# -------------------------
# 4. Plot network with legend
# -------------------------
p_all <- ggraph(g_all, layout = "fr") +
  
  geom_edge_link(
    aes(width = rel, color = tomb_relation),
    alpha = 0.8
  ) +
  
  geom_node_point(
    aes(shape = sex, color = tomb),
    size = 4
  ) +
  
  scale_edge_width_manual(
    values = c(
      "First Degree"  = 2,
      "Second Degree" = 1.2,
      "Third Degree"  = 0.6
    ),
    name = "Kinship Degree"
  ) +
  
  scale_edge_color_manual(
    values = c(
      "Within tomb" = "black",
      "Cross tomb"  = "red"
    ),
    name = "Tomb Relation"
  ) +
  
  scale_color_discrete(name = "Tomb") +
  
  scale_shape_manual(
    name = "Sex",
    values = c("XX" = 16, "XY" = 17), 
    labels = c("XX" = "Female", "XY" = "Male")
  ) +
  
  theme_void() +
  labs(title = "All Tombs Together")


################### Grouped Elateia T46 56 62

set.seed(123)

# -------------------------
# 1. Filter edges for grouped Elateia T46 56 62
# -------------------------
elateia_edges <- pairs_group_og |>
  filter(
    (tomb1 == "Elateia T46 56 62" | tomb2 == "Elateia T46 56 62") &
      rel %in% c("First Degree", "Second Degree", "Third Degree")  # only main relations
  ) |>
  mutate(
    tomb_relation = ifelse(tomb1 == tomb2, "Within tomb", "Cross tomb"),
    rel = factor(rel, levels = c("First Degree", "Second Degree", "Third Degree"))
  )

# -------------------------
# 2. Create node table
# -------------------------
elateia_nodes <- bind_rows(
  elateia_edges |>
    select(individual_id = individual1_id,
           tomb = tomb1,
           sex = sex1),
  elateia_edges |>
    select(individual_id = individual2_id,
           tomb = tomb2,
           sex = sex2)
) |>
  distinct()

# -------------------------
# 3. Build graph
# -------------------------
g_elateia <- graph_from_data_frame(
  d = elateia_edges |>
    select(from = individual1_id,
           to   = individual2_id,
           rel,
           tomb_relation),
  vertices = elateia_nodes |>
    select(name = individual_id, tomb, sex),
  directed = FALSE
)

# -------------------------
# 4. Plot network with legend
# -------------------------
p_elateia_group <- ggraph(g_elateia, layout = "fr") +
  
  geom_edge_link(
    aes(width = rel, color = tomb_relation),
    alpha = 0.8
  ) +
  
  geom_node_point(
    aes(shape = sex, color = tomb),
    size = 4
  ) +
  
  scale_edge_width_manual(
    values = c(
      "First Degree"  = 2,
      "Second Degree" = 1.2,
      "Third Degree"  = 0.6
    ),
    name = "Kinship Degree"
  ) +
  
  scale_edge_color_manual(
    values = c(
      "Within tomb" = "black",
      "Cross tomb"  = "red"
    ),
    name = "Tomb Relation"
  ) +
  
  scale_color_discrete(name = "Tomb") +
  
  scale_shape_manual(
    name = "Sex",
    values = c("XX" = 16, "XY" = 17), 
    labels = c("XX" = "Female", "XY" = "Male")
  ) +
  
  theme_void() +
  labs(title = "Group Elateia T46 56 62")


###saving plots

# plot_list <- list(
#   p_amfissa,
#   p_elateia62,
#   p_all,
#   p_elateia_group
# )
# 
# file_names <- c(
#   "amfissa_tholos_network.png",
#   "elateia_T62_network.png",
#   "all_tombs_network.png",
#   "elateia_T46_56_62_grouped_network.png"
# )
# 
# for (i in seq_along(plot_list)) {
#   ggsave(
#     filename = paste0(
#       "unofficial work/network_plots/",
#       file_names[i]
#     ),
#     plot = plot_list[[i]],
#     width = 10,
#     height = 8,
#     dpi = 300
#   )
# }
