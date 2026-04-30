args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_all[grep(file_arg, args_all)])

if (length(script_path) == 0) {
  script_dir <- getwd()
} else {
  script_dir <- dirname(normalizePath(script_path))
}

base_dir <- normalizePath(file.path(script_dir, ".."))
cluster_dir <- normalizePath(file.path(base_dir, ".."))
output_dir <- file.path(base_dir, "outputs")
desc_dir <- file.path(output_dir, "descriptive_statistics")
viz_dir <- file.path(output_dir, "visualization")

clade1_snv_path <- file.path(cluster_dir, "SNV", "clade1_post_noncore_snv_long.csv")
clade2_snv_path <- file.path(cluster_dir, "SNV", "clade2_post_noncore_snv_long.csv")
metadata_path <- file.path(desc_dir, "metadata_case_unformed_clade1_clade2.csv")

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

plot_width <- 8
plot_height <- 8

metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
names(metadata) <- trimws(names(metadata))

patient_lookup <- metadata %>%
  transmute(
    genome_id = as.character(genome_id),
    patient_id = trimws(as.character(patient_id))
  ) %>%
  mutate(patient_id = na_if(patient_id, "")) %>%
  distinct(genome_id, .keep_all = TRUE)

prepare_snv_pairs <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE) %>%
    rename(snv = snvs) %>%
    mutate(
      id1 = as.character(id1),
      id2 = as.character(id2),
      snv = suppressWarnings(as.numeric(snv))
    ) %>%
    left_join(patient_lookup, by = c("id1" = "genome_id")) %>%
    rename(patient_id_1 = patient_id) %>%
    left_join(patient_lookup, by = c("id2" = "genome_id")) %>%
    rename(patient_id_2 = patient_id) %>%
    mutate(
      same_patient = case_when(
        !is.na(patient_id_1) & !is.na(patient_id_2) & patient_id_1 == patient_id_2 ~ "Within patient",
        !is.na(patient_id_1) & !is.na(patient_id_2) & patient_id_1 != patient_id_2 ~ "Between patient",
        TRUE ~ NA_character_
      )
    )
}

make_patient_histograms <- function(dat, clade_label) {
  dat_base <- dat %>%
    filter(!is.na(same_patient))

  plot_log <- dat_base %>%
    filter(!is.na(snv), snv > 0) %>%
    ggplot(aes(x = snv, fill = same_patient)) +
    geom_histogram(
      bins = 100,
      position = "identity",
      alpha = 0.8
    ) +
    scale_x_log10() +
    scale_fill_manual(
      values = c("Within patient" = "#ECA0B2", "Between patient" = "#54BFB7"),
      name = NULL
    ) +
    labs(
      title = paste0(clade_label, " - post-noncore - Within vs Between patient"),
      x = "SNV distance (log scale)",
      y = "Pairs"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    facet_wrap(~same_patient, ncol = 1, scales = "free_y")

  plot_zoom <- dat_base %>%
    filter(!is.na(snv), snv >= 0, snv <= 10) %>%
    ggplot(aes(x = snv, fill = same_patient)) +
    geom_histogram(
      breaks = seq(-0.5, 10.5, by = 1),
      closed = "right",
      position = "identity",
      alpha = 0.4
    ) +
    scale_fill_manual(
      values = c("Within patient" = "#ECA0B2", "Between patient" = "#54BFB7"),
      name = NULL
    ) +
    scale_x_continuous(
      breaks = 1:10,
      limits = c(0, 10)
    ) +
    labs(
      title = paste0(clade_label, " - post-noncore - Within vs Between patient (0-10 SNVs)"),
      x = "SNV distance",
      y = "Pairs"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  list(log_plot = plot_log, zoom_plot = plot_zoom)
}

clade1_pairs <- prepare_snv_pairs(clade1_snv_path)
clade2_pairs <- prepare_snv_pairs(clade2_snv_path)

clade1_plots <- make_patient_histograms(clade1_pairs, "Clade 1")
clade2_plots <- make_patient_histograms(clade2_pairs, "Clade 2")

write.csv(
  clade1_pairs,
  file.path(desc_dir, "clade1_post_noncore_snv_long_with_patient_labels.csv"),
  row.names = FALSE
)
write.csv(
  clade2_pairs,
  file.path(desc_dir, "clade2_post_noncore_snv_long_with_patient_labels.csv"),
  row.names = FALSE
)

ggsave(
  filename = file.path(viz_dir, "clade1_within_between_patient_snv_log.png"),
  plot = clade1_plots$log_plot,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = clade1_plots$log_plot,
  filename = file.path(viz_dir, "clade1_within_between_patient_snv_log.pdf"),
  width = plot_width,
  height = plot_height
)

ggsave(
  filename = file.path(viz_dir, "clade1_within_between_patient_snv_0_10.png"),
  plot = clade1_plots$zoom_plot,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = clade1_plots$zoom_plot,
  filename = file.path(viz_dir, "clade1_within_between_patient_snv_0_10.pdf"),
  width = plot_width,
  height = plot_height
)

ggsave(
  filename = file.path(viz_dir, "clade2_within_between_patient_snv_log.png"),
  plot = clade2_plots$log_plot,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = clade2_plots$log_plot,
  filename = file.path(viz_dir, "clade2_within_between_patient_snv_log.pdf"),
  width = plot_width,
  height = plot_height
)

ggsave(
  filename = file.path(viz_dir, "clade2_within_between_patient_snv_0_10.png"),
  plot = clade2_plots$zoom_plot,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = clade2_plots$zoom_plot,
  filename = file.path(viz_dir, "clade2_within_between_patient_snv_0_10.pdf"),
  width = plot_width,
  height = plot_height
)

cat("Done. Within-vs-between SNV histograms written to:\n", viz_dir, "\n", sep = "")
