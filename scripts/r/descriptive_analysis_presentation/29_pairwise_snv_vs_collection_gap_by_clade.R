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

required_pkgs <- c("ggplot2", "dplyr", "readr", "tibble")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
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

parse_date_safely <- function(x) {
  as.Date(sub("T.*$", "", as.character(x)))
}

metadata <- read_csv(metadata_path, show_col_types = FALSE) %>%
  transmute(
    genome_id = as.character(genome_id),
    collection_date = parse_date_safely(collection_date),
    clade_group = as.character(clade_group)
  ) %>%
  filter(!is.na(genome_id), !is.na(collection_date)) %>%
  distinct(genome_id, .keep_all = TRUE)

prepare_pairs <- function(path, clade_label) {
  read_csv(path, show_col_types = FALSE) %>%
    rename(snv = snvs) %>%
    mutate(
      id1 = as.character(id1),
      id2 = as.character(id2),
      snv = suppressWarnings(as.numeric(snv))
    ) %>%
    inner_join(metadata, by = c("id1" = "genome_id")) %>%
    rename(collection_date_1 = collection_date, clade_group_1 = clade_group) %>%
    inner_join(metadata, by = c("id2" = "genome_id")) %>%
    rename(collection_date_2 = collection_date, clade_group_2 = clade_group) %>%
    filter(
      clade_group_1 == clade_label,
      clade_group_2 == clade_label,
      !is.na(snv)
    ) %>%
    transmute(
      clade_group = clade_label,
      id1,
      id2,
      snv,
      days_between_collection = abs(as.numeric(collection_date_2 - collection_date_1))
    )
}

clade1_pairs <- prepare_pairs(clade1_snv_path, "Clade 1")
clade2_pairs <- prepare_pairs(clade2_snv_path, "Clade 2")
plot_dat <- bind_rows(clade1_pairs, clade2_pairs) %>%
  mutate(clade_group = factor(clade_group, levels = c("Clade 1", "Clade 2"))) %>%
  filter(snv <= 20)

p <- ggplot(plot_dat, aes(x = days_between_collection, y = snv)) +
  geom_point(alpha = 0.12, size = 1.3, color = "#4E79A7") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1, color = "#E15759") +
  facet_wrap(~clade_group, ncol = 2, scales = "free_y") +
  labs(
    title = "Pairwise SNV distance versus collection gap by clade",
    subtitle = "All pairwise isolate comparisons within each cleaned clade, restricted to SNV <= 20",
    x = "Days between collection dates",
    y = "SNV distance"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(viz_dir, "pairwise_snv_vs_collection_gap_by_clade.png"),
  plot = p,
  width = 11,
  height = 6.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

save_plot_pdf(
  plot_obj = p,
  filename = file.path(viz_dir, "pairwise_snv_vs_collection_gap_by_clade.pdf"),
  width = 11,
  height = 6.5
)

write_csv(
  plot_dat,
  file.path(desc_dir, "pairwise_snv_vs_collection_gap_by_clade_data.csv")
)

cat("Done. Pairwise SNV vs collection-gap outputs written to:\n", viz_dir, "\n", sep = "")
