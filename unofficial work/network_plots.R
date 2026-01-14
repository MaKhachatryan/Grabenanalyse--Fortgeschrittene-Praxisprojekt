# Load environment and data
source(here::here("environmentSetUp.R"))
source(here::here("analysis", "pairsdata.R"))

############### Amfissa and Elateia


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
# Manual filename map (fixed casing)
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
# Loop over tombs
# -------------------------
for (t in tombs) {
  
  # 1. Filter edges
  edges <- pairs_ind_og |>
    filter(
      rel %in% c("First Degree", "Second Degree", "Third Degree"),
      tomb1 == t | tomb2 == t
    ) |>
    mutate(
      tomb_relation = ifelse(tomb1 == tomb2, "Within tomb", "Cross tomb"),
      rel = factor(rel, levels = c("First Degree", "Second Degree", "Third Degree"))
    )
  
  # Skip if no relations
  if (nrow(edges) == 0) {
    message("No relations for ", t, " — skipped")
    next
  }
  
  # 2. Nodes
  nodes <- bind_rows(
    edges |>
      select(individual_id = individual1_id, tomb = tomb1, sex = sex1),
    edges |>
      select(individual_id = individual2_id, tomb = tomb2, sex = sex2)
  ) |>
    distinct()
  
  # 3. Graph
  g <- graph_from_data_frame(
    d = edges |>
      select(from = individual1_id, to = individual2_id, rel, tomb_relation),
    vertices = nodes |>
      select(name = individual_id, tomb, sex),
    directed = FALSE
  )
  
  # 4. Plot
  p <- ggraph(g, layout = "fr") +
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
    labs(title = t)
  
  # 5. Save (only if file does NOT already exist)
  out_file <- file.path(
    "unofficial work",
    "network_plots",
    file_map[t]
  )
  
  if (!file.exists(out_file)) {
    ggsave(
      filename = out_file,
      plot  = p,
      width = 8,
      height = 6,
      dpi   = 300
    )
  } else {
    message("File exists, skipped: ", file_map[t])
  }
  
}


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

# -------------------------
# 5. Save 
# -------------------------
out_file <- file.path(
  "unofficial work",
  "network_plots",
  "All_Tombs_network.png"
)

if (!file.exists(out_file)) {
  ggsave(
    filename = out_file,
    plot  = p_all,
    width = 8,
    height = 6,
    dpi   = 300
  )
} else {
  message("File exists, skipped: All_Tombs_network.png")
}


################### Grouped Elateia T46 56 62

set.seed(123)

# -------------------------
# 1. Filter edges for grouped Elateia T46 56 62
# -------------------------
elateia_edges <- pairs_group_og |>
  filter(
    (tomb1 == "Elateia T46 56 62" | tomb2 == "Elateia T46 56 62") &
      rel %in% c("First Degree", "Second Degree", "Third Degree")
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

# -------------------------
# 5. Save
# -------------------------
out_file <- file.path(
  "unofficial work",
  "network_plots",
  "Elateia_T46_56_62_network.png"
)

if (!file.exists(out_file)) {
  ggsave(
    filename = out_file,
    plot  = p_elateia_group,
    width = 8,
    height = 6,
    dpi   = 300
  )
} else {
  message("File exists, skipped: Elateia_T46_56_62_network.png")
}
