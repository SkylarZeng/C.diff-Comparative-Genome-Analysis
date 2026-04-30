args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_all[grep(file_arg, args_all)])

if (length(script_path) == 0) {
  script_dir <- getwd()
} else {
  script_dir <- dirname(normalizePath(script_path))
}

base_dir <- normalizePath(file.path(script_dir, ".."))
output_dir <- file.path(base_dir, "outputs")
desc_dir <- file.path(output_dir, "descriptive_statistics")
viz_dir <- file.path(output_dir, "visualization")

locations_path <- file.path(desc_dir, "locations_filtered_case_unformed_clade1_clade2.tsv")
metadata_path <- file.path(desc_dir, "metadata_case_unformed_clade1_clade2.csv")

if (!dir.exists(viz_dir)) {
  dir.create(viz_dir, recursive = TRUE)
}

required_pkgs <- c("ggplot2", "dplyr", "readr", "tibble", "stringr", "tidyr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(tidyr)
})

save_plot_pdf <- function(plot_obj, filename, width, height) {
  ggsave(
    filename = filename,
    plot = plot_obj,
    device = "pdf",
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
}

is_present <- function(x) {
  !is.na(x) & trimws(as.character(x)) != "" & trimws(as.character(x)) != "NA"
}

locations <- read.delim(locations_path, check.names = FALSE, stringsAsFactors = FALSE) %>%
  mutate(
    PatientID = trimws(as.character(PatientID)),
    clade_group = trimws(as.character(clade_group)),
    FacilityCode = trimws(as.character(FacilityCode)),
    Unit = trimws(as.character(Unit))
  )

metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE) %>%
  mutate(
    patient_id = trimws(as.character(patient_id)),
    clade_group = trimws(as.character(clade_group))
  )

patient_totals <- metadata %>%
  filter(clade_group %in% c("Clade 1", "Clade 2"), !is.na(patient_id), patient_id != "") %>%
  distinct(patient_id, clade_group) %>%
  count(clade_group, name = "total_patients")

build_hierarchy <- function(clade_label) {
  total_patients <- patient_totals %>%
    filter(clade_group == clade_label) %>%
    pull(total_patients)

  loc_clade <- locations %>%
    filter(clade_group == clade_label)

  facility_summary <- loc_clade %>%
    filter(is_present(FacilityCode)) %>%
    distinct(PatientID, FacilityCode) %>%
    count(FacilityCode, name = "patient_n") %>%
    arrange(desc(patient_n), FacilityCode)

  unit_summary <- loc_clade %>%
    filter(is_present(FacilityCode), is_present(Unit)) %>%
    distinct(PatientID, FacilityCode, Unit) %>%
    count(FacilityCode, Unit, name = "patient_n") %>%
    arrange(FacilityCode, desc(patient_n), Unit)

  facility_nodes <- facility_summary %>%
    mutate(
      facility_id = row_number(),
      facility_label = paste0(FacilityCode, "\n", "n=", patient_n)
    )

  if (nrow(facility_nodes) == 0) {
    stop("No facility data found for ", clade_label)
  }

  unit_nodes <- unit_summary %>%
    inner_join(facility_nodes %>% select(FacilityCode, facility_id), by = "FacilityCode") %>%
    group_by(FacilityCode) %>%
    arrange(desc(patient_n), Unit, .by_group = TRUE) %>%
    mutate(unit_order = row_number()) %>%
    ungroup() %>%
    mutate(
      unit_label = paste0(Unit, "\n", "n=", patient_n)
    )

  gap_between_facilities <- 1
  unit_y_map <- numeric(nrow(unit_nodes))
  current_y <- 1

  for (i in seq_len(nrow(facility_nodes))) {
    facility_code <- facility_nodes$FacilityCode[i]
    idx <- which(unit_nodes$FacilityCode == facility_code)
    if (length(idx) == 0) {
      facility_nodes$y[i] <- current_y
      current_y <- current_y + gap_between_facilities + 1
    } else {
      y_vals <- seq(current_y, by = 1, length.out = length(idx))
      unit_y_map[idx] <- y_vals
      facility_nodes$y[i] <- mean(y_vals)
      current_y <- max(y_vals) + gap_between_facilities + 1
    }
  }

  unit_nodes$y <- unit_y_map
  facility_nodes$x <- 2
  unit_nodes$x <- 3.6

  total_node <- tibble(
    node_type = "Total",
    node_name = clade_label,
    patient_n = total_patients,
    label = paste0(clade_label, "\n", "n=", total_patients),
    x = 0.4,
    y = mean(facility_nodes$y)
  )

  facility_plot_nodes <- facility_nodes %>%
    transmute(
      node_type = "Facility",
      node_name = FacilityCode,
      patient_n = patient_n,
      label = facility_label,
      x = x,
      y = y
    )

  unit_plot_nodes <- unit_nodes %>%
    transmute(
      node_type = "Unit",
      node_name = Unit,
      patient_n = patient_n,
      label = unit_label,
      x = x,
      y = y,
      FacilityCode = FacilityCode
    )

  total_to_facility_edges <- facility_nodes %>%
    transmute(
      edge_type = "total_to_facility",
      x = total_node$x + 0.28,
      y = total_node$y,
      xend = x - 0.28,
      yend = y
    )

  facility_to_unit_edges <- unit_nodes %>%
    inner_join(facility_nodes %>% select(FacilityCode, facility_y = y), by = "FacilityCode") %>%
    transmute(
      edge_type = "facility_to_unit",
      x = 2.28,
      y = facility_y,
      xend = 3.32,
      yend = y
    )

  list(
    nodes = bind_rows(total_node, facility_plot_nodes, unit_plot_nodes),
    edges = bind_rows(total_to_facility_edges, facility_to_unit_edges),
    facility_summary = facility_summary,
    unit_summary = unit_summary
  )
}

make_branch_plot <- function(clade_label, fill_total, fill_facility, fill_unit) {
  hierarchy <- build_hierarchy(clade_label)
  nodes <- hierarchy$nodes
  edges <- hierarchy$edges

  box_width_lookup <- c("Total" = 0.55, "Facility" = 0.52, "Unit" = 0.48)
  box_height_lookup <- c("Total" = 0.9, "Facility" = 0.82, "Unit" = 0.74)
  fill_lookup <- c("Total" = fill_total, "Facility" = fill_facility, "Unit" = fill_unit)
  text_size_lookup <- c("Total" = 4.2, "Facility" = 3.4, "Unit" = 2.8)

  nodes <- nodes %>%
    mutate(
      xmin = x - unname(box_width_lookup[node_type]) / 2,
      xmax = x + unname(box_width_lookup[node_type]) / 2,
      ymin = y - unname(box_height_lookup[node_type]) / 2,
      ymax = y + unname(box_height_lookup[node_type]) / 2,
      fill_color = unname(fill_lookup[node_type]),
      text_size = unname(text_size_lookup[node_type])
    )

  plot_height <- max(8, nrow(nodes %>% filter(node_type == "Unit")) * 0.28 + 3)

  p <- ggplot() +
    geom_segment(
      data = edges,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.45,
      color = "#94a3b8"
    ) +
    geom_rect(
      data = nodes,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = nodes$fill_color,
      color = NA
    ) +
    geom_text(
      data = nodes,
      aes(x = x, y = y, label = label, size = text_size),
      color = "white",
      fontface = "bold",
      lineheight = 0.95,
      show.legend = FALSE
    ) +
    scale_size_identity() +
    annotate("text", x = 0.4, y = max(nodes$y) + 1.6, label = "Clade total", fontface = 2, size = 4) +
    annotate("text", x = 2, y = max(nodes$y) + 1.6, label = "Facilities", fontface = 2, size = 4) +
    annotate("text", x = 3.6, y = max(nodes$y) + 1.6, label = "Units", fontface = 2, size = 4) +
    coord_cartesian(
      xlim = c(-0.1, 4.2),
      ylim = c(0, max(nodes$y) + 2.2),
      clip = "off"
    ) +
    labs(
      title = paste0(clade_label, ": patient distribution across facilities and units"),
      subtitle = "Node labels show distinct patient counts"
    ) +
    theme_void(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      plot.margin = margin(20, 20, 20, 20)
    )

  list(plot = p, height = plot_height, hierarchy = hierarchy)
}

clade1 <- make_branch_plot("Clade 1", "#1d4e89", "#4E79A7", "#9ecae1")
clade2 <- make_branch_plot("Clade 2", "#8f2d21", "#E15759", "#f4a6a8")

hierarchy_summary <- bind_rows(
  clade1$hierarchy$facility_summary %>% mutate(clade_group = "Clade 1", level = "Facility", parent = "Clade 1", name = FacilityCode) %>% select(clade_group, level, parent, name, patient_n),
  clade1$hierarchy$unit_summary %>% mutate(clade_group = "Clade 1", level = "Unit", parent = FacilityCode, name = Unit) %>% select(clade_group, level, parent, name, patient_n),
  clade2$hierarchy$facility_summary %>% mutate(clade_group = "Clade 2", level = "Facility", parent = "Clade 2", name = FacilityCode) %>% select(clade_group, level, parent, name, patient_n),
  clade2$hierarchy$unit_summary %>% mutate(clade_group = "Clade 2", level = "Unit", parent = FacilityCode, name = Unit) %>% select(clade_group, level, parent, name, patient_n)
)

write.csv(
  hierarchy_summary,
  file.path(desc_dir, "location_facility_unit_hierarchy_counts.csv"),
  row.names = FALSE
)

ggsave(
  filename = file.path(viz_dir, "clade1_location_facility_unit_branch.png"),
  plot = clade1$plot,
  width = 16,
  height = clade1$height,
  units = "in",
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)
save_plot_pdf(
  plot_obj = clade1$plot,
  filename = file.path(viz_dir, "clade1_location_facility_unit_branch.pdf"),
  width = 16,
  height = clade1$height
)

ggsave(
  filename = file.path(viz_dir, "clade2_location_facility_unit_branch.png"),
  plot = clade2$plot,
  width = 16,
  height = clade2$height,
  units = "in",
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)
save_plot_pdf(
  plot_obj = clade2$plot,
  filename = file.path(viz_dir, "clade2_location_facility_unit_branch.pdf"),
  width = 16,
  height = clade2$height
)

cat("Done. Facility-unit branch diagrams written to:\n", viz_dir, "\n", sep = "")
