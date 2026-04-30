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

metadata_path <- file.path(desc_dir, "metadata_case_unformed_clade1_clade2.csv")
locations_path <- file.path(desc_dir, "locations_filtered_case_unformed_clade1_clade2.tsv")

if (!dir.exists(viz_dir)) {
  dir.create(viz_dir, recursive = TRUE)
}

required_pkgs <- c("ggplot2", "tidyverse")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
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
  !is.na(x) & trimws(as.character(x)) != ""
}

metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
names(metadata) <- trimws(names(metadata))

locations <- read.delim(locations_path, check.names = FALSE, stringsAsFactors = FALSE)

patient_clade <- metadata %>%
  transmute(
    patient_id = trimws(as.character(patient_id)),
    clade_group = trimws(as.character(clade_group))
  ) %>%
  filter(!is.na(patient_id), patient_id != "", clade_group %in% c("Clade 1", "Clade 2")) %>%
  distinct(patient_id, clade_group)

locations <- locations %>%
  mutate(PatientID = trimws(as.character(PatientID)))

make_availability_counts <- function(clade_label) {
  total_patients <- patient_clade %>%
    filter(clade_group == clade_label) %>%
    summarise(n = n_distinct(patient_id)) %>%
    pull(n)

  loc_clade <- locations %>%
    filter(clade_group == clade_label)

  facility_patients <- loc_clade %>%
    filter(is_present(FacilityCode)) %>%
    summarise(n = n_distinct(PatientID)) %>%
    pull(n)

  unit_patients <- loc_clade %>%
    filter(is_present(FacilityCode), is_present(Unit)) %>%
    summarise(n = n_distinct(PatientID)) %>%
    pull(n)

  room_patients <- loc_clade %>%
    filter(is_present(FacilityCode), is_present(Unit), is_present(Room)) %>%
    summarise(n = n_distinct(PatientID)) %>%
    pull(n)

  bed_patients <- loc_clade %>%
    filter(is_present(FacilityCode), is_present(Unit), is_present(Room), is_present(Bed)) %>%
    summarise(n = n_distinct(PatientID)) %>%
    pull(n)

  tibble(
    stage = c("Total patients", "Facility data", "Unit data", "Room data", "Bed data"),
    n = c(total_patients, facility_patients, unit_patients, room_patients, bed_patients),
    pct_total = c(total_patients, facility_patients, unit_patients, room_patients, bed_patients) / total_patients
  )
}

make_diagram_plot <- function(clade_label, fill_color) {
  counts <- make_availability_counts(clade_label) %>%
    mutate(
      x = c(1, 3, 5, 7, 9),
      y = c(1, 1, 1, 1, 1),
      label = paste0(stage, "\n", "n=", n, " (", sprintf("%.1f%%", 100 * pct_total), ")")
    )

  arrows <- tibble(
    x = counts$x[-nrow(counts)] + 0.55,
    xend = counts$x[-1] - 0.55,
    y = 1,
    yend = 1
  )

  ggplot() +
    geom_segment(
      data = arrows,
      aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = 1,
      arrow = grid::arrow(length = grid::unit(0.18, "inches"), type = "closed"),
      color = "#6b7280"
    ) +
    geom_label(
      data = counts,
      aes(x = x, y = y, label = label),
      fill = fill_color,
      color = "white",
      fontface = "bold",
      size = 5,
      label.size = 0,
      label.padding = grid::unit(0.35, "lines")
    ) +
    coord_cartesian(xlim = c(0.2, 9.8), ylim = c(0.5, 1.5), clip = "off") +
    labs(
      title = paste0(clade_label, ": Location data availability"),
      subtitle = "Patient cascade from total cohort to Facility -> Unit -> Room -> Bed"
    ) +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      plot.margin = margin(20, 20, 20, 20)
    )
}

clade1_plot <- make_diagram_plot("Clade 1", "#4E79A7")
clade2_plot <- make_diagram_plot("Clade 2", "#E15759")

clade1_counts <- make_availability_counts("Clade 1") %>% mutate(clade_group = "Clade 1")
clade2_counts <- make_availability_counts("Clade 2") %>% mutate(clade_group = "Clade 2")

write.csv(
  bind_rows(clade1_counts, clade2_counts) %>% select(clade_group, everything()),
  file.path(desc_dir, "location_data_availability_diagram_counts.csv"),
  row.names = FALSE
)

ggsave(
  filename = file.path(viz_dir, "clade1_location_data_availability_diagram.png"),
  plot = clade1_plot,
  width = 16,
  height = 5,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = clade1_plot,
  filename = file.path(viz_dir, "clade1_location_data_availability_diagram.pdf"),
  width = 16,
  height = 5
)

ggsave(
  filename = file.path(viz_dir, "clade2_location_data_availability_diagram.png"),
  plot = clade2_plot,
  width = 16,
  height = 5,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = clade2_plot,
  filename = file.path(viz_dir, "clade2_location_data_availability_diagram.pdf"),
  width = 16,
  height = 5
)

cat("Done. Location availability diagrams written to:\n", viz_dir, "\n", sep = "")
