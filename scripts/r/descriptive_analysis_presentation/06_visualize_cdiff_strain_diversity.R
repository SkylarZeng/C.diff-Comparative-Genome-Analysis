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

if (!dir.exists(viz_dir)) {
  dir.create(viz_dir, recursive = TRUE)
}

required_pkgs <- c("ggplot2", "tidyverse", "scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(scales)
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

plot_width <- 24
plot_height <- 9

metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
names(metadata) <- trimws(names(metadata))

metadata <- metadata %>%
  mutate(
    ST = as.character(ST),
    ST = stringr::str_remove_all(ST, "\\*"),
    ST = trimws(ST),
    ST = na_if(ST, ""),
    clade_group = trimws(as.character(clade_group)),
    collection_date = as.Date(sub("T.*$", "", as.character(collection_date)))
  ) %>%
  filter(!is.na(ST), !is.na(collection_date), clade_group %in% c("Clade 1", "Clade 2")) %>%
  mutate(ST_label = paste0("ST ", ST))

make_st_palette <- function(dat_clade, palette_values) {
  st_levels <- dat_clade %>%
    count(ST_label, sort = TRUE) %>%
    pull(ST_label)

  setNames(rep_len(palette_values, length(st_levels)), st_levels)
}

make_tableau_palette <- function(dat_clade) {
  st_levels <- dat_clade %>%
    count(ST_label, sort = TRUE) %>%
    pull(ST_label)

  tableau_base <- c(
    "#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F",
    "#8CD17D", "#B6992D", "#F1CE63", "#499894", "#86BCB6",
    "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", "#D37295",
    "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6"
  )
  if (length(st_levels) <= length(tableau_base)) {
    colors_out <- tableau_base[seq_along(st_levels)]
  } else {
    colors_out <- grDevices::colorRampPalette(tableau_base)(length(st_levels))
  }

  setNames(colors_out, st_levels)
}

make_original_clade2_palette <- function(dat_clade) {
  st_levels <- dat_clade %>%
    count(ST_label, sort = TRUE) %>%
    pull(ST_label)

  clade2_base <- c(
    "#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD",
    "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA",
    "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77",
    "#771122", "#AA4455", "#DD7788"
  )

  colors_out <- grDevices::colorRampPalette(clade2_base)(length(st_levels))
  setNames(colors_out, st_levels)
}

make_diversity_plot <- function(dat, clade_label, legend_rows, palette_mode = c("tableau", "clade2_original")) {
  palette_mode <- match.arg(palette_mode)
  dat_clade <- dat %>%
    filter(clade_group == clade_label)

  if (nrow(dat_clade) == 0) {
    stop("No rows found for ", clade_label)
  }

  colors_gg <- if (palette_mode == "clade2_original") {
    make_original_clade2_palette(dat_clade)
  } else {
    make_tableau_palette(dat_clade)
  }
  legend_breaks <- names(colors_gg)
  legend_labels <- stringr::str_remove(legend_breaks, "^ST\\s+")

  plot_obj <- ggplot(dat_clade, aes(x = collection_date, fill = factor(ST_label, levels = names(colors_gg)))) +
    geom_density(
      aes(y = after_stat(count)),
      position = "fill",
      adjust = 0.35,
      alpha = 0.95
    ) +
    scale_fill_manual(
      values = colors_gg,
      breaks = legend_breaks,
      labels = legend_labels,
      name = "ST",
      guide = guide_legend(nrow = legend_rows, byrow = TRUE)
    ) +
    scale_y_continuous(labels = percent_format(scale = 100)) +
    labs(
      title = paste0(clade_label, ": C. difficile strain diversity over time"),
      x = "Collection Date",
      y = "Prevalence"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      text = element_text(size = 22),
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(15, 20, 60, 15),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.key.width = unit(1.2, "cm"),
      legend.key.height = unit(0.6, "cm"),
      legend.box.margin = margin(5, 5, 5, 5)
    )

  list(
    plot = plot_obj,
    color_key = tibble(ST_label = names(colors_gg), color = unname(colors_gg))
  )
}

clade1_res <- make_diversity_plot(metadata, "Clade 1", legend_rows = 3, palette_mode = "tableau")
clade2_res <- make_diversity_plot(metadata, "Clade 2", legend_rows = 2, palette_mode = "clade2_original")

write.csv(
  clade1_res$color_key,
  file.path(desc_dir, "clade1_strain_diversity_color_key.csv"),
  row.names = FALSE
)

write.csv(
  clade2_res$color_key,
  file.path(desc_dir, "clade2_strain_diversity_color_key.csv"),
  row.names = FALSE
)

ggsave(
  filename = file.path(viz_dir, "clade1_cdiff_strain_diversity.png"),
  plot = clade1_res$plot,
  device = "png",
  height = plot_height,
  width = plot_width,
  units = "in",
  dpi = 300,
  bg = "white"
)

save_plot_pdf(
  plot_obj = clade1_res$plot,
  filename = file.path(viz_dir, "clade1_cdiff_strain_diversity.pdf"),
  width = plot_width,
  height = plot_height
)

ggsave(
  filename = file.path(viz_dir, "clade2_cdiff_strain_diversity.png"),
  plot = clade2_res$plot,
  device = "png",
  height = plot_height,
  width = plot_width,
  units = "in",
  dpi = 300,
  bg = "white"
)

save_plot_pdf(
  plot_obj = clade2_res$plot,
  filename = file.path(viz_dir, "clade2_cdiff_strain_diversity.pdf"),
  width = plot_width,
  height = plot_height
)

cat("Done. Strain diversity plots written to:\n", viz_dir, "\n", sep = "")
