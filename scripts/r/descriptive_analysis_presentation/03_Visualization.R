args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_all[grep(file_arg, args_all)])

if (length(script_path) == 0) {
  script_dir <- getwd()
} else {
  script_dir <- dirname(normalizePath(script_path))
}

base_dir <- normalizePath(file.path(script_dir, "."))
cluster_dir <- normalizePath(file.path(base_dir, "."))
repo_dir <- normalizePath(file.path(cluster_dir, "."))

output_dir <- file.path(base_dir, "Cluster_visulization_clade2", "descriptive_analysis_presentation", "outputs")
desc_dir <- file.path(output_dir, "descriptive_statistics")
viz_dir <- file.path(output_dir, "visualization")

metadata_path <- file.path(desc_dir, "metadata_case_unformed.csv")
tree_path <- file.path(repo_dir, "2025-07-01_cognac", "iTOL", "cognac_tree.nwk")

if (!dir.exists(viz_dir)) {
  dir.create(viz_dir, recursive = TRUE)
}

required_pkgs <- c("ape", "ggtree", "ggplot2", "tidyverse", "tibble")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(tidyverse)
})

metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE) %>%
  mutate(
    genome_id = as.character(genome_id),
    ST = as.character(ST),
    clade_group = as.character(clade_group),
    ST = if_else(is.na(ST) | ST == "" | ST == "NA", "ST NA", paste0("ST ", ST)),
    clade_group = if_else(is.na(clade_group) | clade_group == "", "Unknown", clade_group)
  ) %>%
  distinct(genome_id, .keep_all = TRUE)

tree <- read.tree(tree_path)

metadata_ids <- metadata$genome_id
tree_ids <- tree$tip.label
matched_ids <- intersect(tree_ids, metadata_ids)
metadata_missing_from_tree <- setdiff(metadata_ids, tree_ids)
tree_extra_tips <- setdiff(tree_ids, metadata_ids)

if (length(matched_ids) == 0) {
  stop("No overlapping genome IDs found between metadata_case_unformed.csv and cognac_tree.nwk")
}

filtered_tree <- keep.tip(tree, matched_ids)
filtered_metadata <- metadata %>%
  filter(genome_id %in% matched_ids)

plot_height <- max(18, nrow(filtered_metadata) * 0.04)

p <- ggtree(filtered_tree, branch.length = "none", linewidth = 0.2)

tip_positions <- p$data %>%
  filter(isTip) %>%
  select(label, x, y)

annotation_data <- filtered_metadata %>%
  transmute(
    label = genome_id,
    clade_label = clade_group,
    ST_label = ST
  ) %>%
  left_join(tip_positions, by = "label")

clade_colors <- c(
  "Clade 1" = "#367db8",
  "Clade 2" = "#fecbe3",
  "Clade 3" = "#fcb359",
  "Clade 4" = "#be7fbc",
  "Clade 5" = "#87d5c6",
  "Clade cryptic" = "#49b243",
  "Unknown" = "#bdbdbd"
)

st_levels <- sort(unique(annotation_data$ST_label))
st_colors <- setNames(grDevices::hcl.colors(length(st_levels), palette = "Dark 3"), st_levels)
if ("ST NA" %in% names(st_colors)) {
  st_colors["ST NA"] <- "#7f7f7f"
}

x_tree_max <- max(tip_positions$x, na.rm = TRUE)
x_clade <- x_tree_max + 2.5
x_st <- x_tree_max + 5.5
y_header <- max(tip_positions$y, na.rm = TRUE) + 20

p_annotated <- p +
  geom_tiplab(size = 1.1, align = TRUE, linetype = NA, offset = 0.3) +
  geom_tile(
    data = annotation_data,
    aes(x = x_clade, y = y, fill = clade_label),
    width = 0.8,
    height = 0.8,
    inherit.aes = FALSE
  ) +
  geom_tile(
    data = annotation_data,
    aes(x = x_st, y = y, fill = ST_label),
    width = 0.8,
    height = 0.8,
    inherit.aes = FALSE
  ) +
  annotate("text", x = x_clade, y = y_header, label = "Clade", fontface = 2, size = 3) +
  annotate("text", x = x_st, y = y_header, label = "ST", fontface = 2, size = 3, hjust = 0) +
  scale_fill_manual(values = clade_colors, drop = FALSE, name = "Clade") +
  scale_color_manual(values = st_colors, drop = FALSE, guide = "none") +
  coord_cartesian(xlim = c(0, x_st + 3), clip = "off") +
  theme(
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 80, 10, 10)
  )

ggsave(
  filename = file.path(viz_dir, "case_unformed_cognac_tree_ST_clade.pdf"),
  plot = p_annotated,
  device = "pdf",
  width = 25,
  height = plot_height,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = file.path(viz_dir, "case_unformed_cognac_tree_ST_clade.png"),
  plot = p_annotated,
  device = "png",
  width = 25,
  height = plot_height,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

cat("Done. Visualization outputs written to:\n", viz_dir, "\n", sep = "")
