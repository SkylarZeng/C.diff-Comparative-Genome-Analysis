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

clade1_pairs_path <- file.path(desc_dir, "clade1_post_noncore_snv_long_with_location_overlap.csv")
clade2_pairs_path <- file.path(desc_dir, "clade2_post_noncore_snv_long_with_location_overlap.csv")
tests_path <- file.path(desc_dir, "location_overlap_and_sequential_snv_wilcox_tests.csv")

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

sequential_class_label <- "Sequential_14d"
sequential_file_stub <- "14d"

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

prepare_pairs <- function(path, clade_label) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE) %>%
    mutate(
      snv = suppressWarnings(as.numeric(snv)),
      clade_group = clade_label,
      same_patient = trimws(as.character(same_patient))
    ) %>%
    filter(same_patient == "Between patient", !is.na(snv), snv > 0)
}

make_star_label <- function(p) {
  case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

clade1_pairs <- prepare_pairs(clade1_pairs_path, "Clade 1")
clade2_pairs <- prepare_pairs(clade2_pairs_path, "Clade 2")
all_pairs <- bind_rows(clade1_pairs, clade2_pairs)

tests <- read.csv(tests_path, check.names = FALSE, stringsAsFactors = FALSE) %>%
  mutate(
    clade_group = as.character(clade_group),
    exposure_class = as.character(exposure_class),
    exposure_type = factor(as.character(exposure_type), levels = c("Facility", "Floor", "Unit", "Room", "Bed")),
    star_label = make_star_label(as.numeric(wilcox_p_value))
  ) %>%
  filter(clade_group %in% c("Clade 1", "Clade 2"))

build_plot_data <- function(dat, exposure_mode = c("Overlap", "Sequential_14d")) {
  exposure_mode <- match.arg(exposure_mode)

  selected_cols <- if (exposure_mode == "Overlap") {
    c("facility_overlap", "floor_overlap", "unit_overlap", "room_overlap")
  } else {
    c("facility_sequential_14d", "floor_sequential_14d", "unit_sequential_14d", "room_sequential_14d")
  }

  dat %>%
    select(clade_group, snv, all_of(selected_cols)) %>%
    pivot_longer(
      cols = all_of(selected_cols),
      names_to = "exposure_type",
      values_to = "exposure_flag"
    ) %>%
    mutate(
      exposure_type = recode(
        exposure_type,
        facility_overlap = "Facility",
        floor_overlap = "Floor",
        unit_overlap = "Unit",
        room_overlap = "Room",
        facility_sequential_14d = "Facility",
        floor_sequential_14d = "Floor",
        unit_sequential_14d = "Unit",
        room_sequential_14d = "Room"
      ),
      exposure_type = factor(exposure_type, levels = c("Facility", "Floor", "Unit", "Room", "Bed")),
      exposure_status = if (exposure_mode == "Overlap") {
        if_else(exposure_flag, "Yes", "No")
      } else {
        if_else(exposure_flag, "Yes", "No")
      }
    )
}

make_boxplot <- function(plot_dat, tests_dat, exposure_mode = c("Overlap", "Sequential_14d")) {
  exposure_mode <- match.arg(exposure_mode)

  subtitle_text <- if (exposure_mode == "Overlap") {
    "Between-patient SNV distances by direct location overlap"
  } else {
    "Between-patient SNV distances by sequential same-location exposure within 14 days"
  }

  ann_dat <- plot_dat %>%
    group_by(clade_group, exposure_type) %>%
    summarise(y_pos = max(snv, na.rm = TRUE) * 1.12, .groups = "drop") %>%
    left_join(
      tests_dat %>%
        filter(exposure_class == exposure_mode) %>%
        select(clade_group, exposure_type, star_label),
      by = c("clade_group", "exposure_type")
    )

  ggplot(plot_dat, aes(x = exposure_type, y = snv, fill = exposure_status)) +
    geom_boxplot(
      position = position_dodge(width = 0.75),
      width = 0.65,
      outlier.alpha = 0.08,
      outlier.size = 0.4
    ) +
    geom_text(
      data = ann_dat,
      aes(x = exposure_type, y = y_pos, label = star_label),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 6
    ) +
    facet_wrap(~clade_group, ncol = 2, scales = "free_y") +
    scale_fill_manual(
      values = c("Yes" = "#e15759", "No" = "#4e79a7"),
      name = NULL
    ) +
    scale_y_log10() +
    labs(
      title = if (exposure_mode == "Overlap") {
        "SNV distance by location overlap"
      } else {
        "SNV distance by sequential 14-day location exposure"
      },
      subtitle = subtitle_text,
      x = NULL,
      y = "SNV distance (log scale)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

overlap_plot_dat <- build_plot_data(all_pairs, "Overlap")
sequential_plot_dat <- build_plot_data(all_pairs, "Sequential_14d")

overlap_plot <- make_boxplot(overlap_plot_dat, tests, "Overlap")
sequential_plot <- make_boxplot(sequential_plot_dat, tests, "Sequential_14d")

ggsave(
  filename = file.path(viz_dir, "location_overlap_snv_boxplot_with_significance.png"),
  plot = overlap_plot,
  width = 14,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = overlap_plot,
  filename = file.path(viz_dir, "location_overlap_snv_boxplot_with_significance.pdf"),
  width = 14,
  height = 8
)

ggsave(
  filename = file.path(viz_dir, "location_sequential_14d_snv_boxplot_with_significance.png"),
  plot = sequential_plot,
  width = 14,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = sequential_plot,
  filename = file.path(viz_dir, "location_sequential_14d_snv_boxplot_with_significance.pdf"),
  width = 14,
  height = 8
)

cat("Done. Location overlap/sequential boxplots written to:\n", viz_dir, "\n", sep = "")
