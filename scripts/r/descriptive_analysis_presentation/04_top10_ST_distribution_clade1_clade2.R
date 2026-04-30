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
metadata_path <- file.path(desc_dir, "metadata_case_unformed.csv")

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

clade2_colors <- c(
  "#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD",
  "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA",
  "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77",
  "#771122", "#AA4455", "#DD7788"
)

clade1_colors <- c(
  "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c",
  "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
  "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f",
  "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"
)

metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
names(metadata) <- trimws(names(metadata))

metadata <- metadata %>%
  mutate(
    ST = as.character(ST),
    ST = stringr::str_remove_all(ST, "\\*"),
    ST = trimws(ST),
    ST = na_if(ST, ""),
    clade_group = trimws(as.character(clade_group))
  )

make_top10_st_plot <- function(dat, clade_label, palette_values) {
  dat_clade <- dat %>%
    filter(clade_group == clade_label, !is.na(ST))

  if (nrow(dat_clade) == 0) {
    stop("No rows found for ", clade_label)
  }

  st_counts <- dat_clade %>%
    count(ST, name = "ST_n", sort = TRUE) %>%
    slice_head(n = 10) %>%
    mutate(
      binned_ST = paste0("ST ", ST),
      ST_color = palette_values[seq_len(n())]
    )

  dat_binned <- dat_clade %>%
    left_join(st_counts %>% select(ST, binned_ST, ST_color), by = "ST") %>%
    filter(!is.na(binned_ST))

  summary_binned <- dat_binned %>%
    count(binned_ST, ST_color, name = "n") %>%
    mutate(total_n = n) %>%
    ungroup()

  st_levels <- summary_binned %>%
    distinct(binned_ST, total_n) %>%
    arrange(total_n, binned_ST) %>%
    pull(binned_ST)

  st_palette <- st_counts %>%
    select(binned_ST, ST_color) %>%
    deframe()

  plot_obj <- ggplot(
    summary_binned,
    aes(y = factor(binned_ST, levels = st_levels), x = n, fill = binned_ST)
  ) +
    geom_col(width = 0.75) +
    geom_text(
      aes(label = paste0("n=", n)),
      hjust = -0.1,
      size = 3.8
    ) +
    scale_fill_manual(values = st_palette, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = paste0(clade_label, ": Top 10 ST Distribution"),
      x = "Number of isolates",
      y = "Sequence type"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  list(
    plot = plot_obj,
    summary = summary_binned %>% arrange(desc(total_n), binned_ST),
    top10 = st_counts
  )
}

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

plot_width <- 6
plot_height <- 4

clade1_res <- make_top10_st_plot(metadata, "Clade 1", clade1_colors)
clade2_res <- make_top10_st_plot(metadata, "Clade 2", clade2_colors)

write.csv(
  clade1_res$summary,
  file.path(desc_dir, "clade1_top10_ST_distribution_summary.csv"),
  row.names = FALSE
)
write.csv(
  clade2_res$summary,
  file.path(desc_dir, "clade2_top10_ST_distribution_summary.csv"),
  row.names = FALSE
)

ggsave(
  filename = file.path(viz_dir, "clade1_top10_ST_distribution.png"),
  plot = clade1_res$plot,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 300,
  bg = "white"
)

save_plot_pdf(
  plot_obj = clade1_res$plot,
  filename = file.path(viz_dir, "clade1_top10_ST_distribution.pdf"),
  width = plot_width,
  height = plot_height
)

ggsave(
  filename = file.path(viz_dir, "clade2_top10_ST_distribution.png"),
  plot = clade2_res$plot,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 300,
  bg = "white"
)

save_plot_pdf(
  plot_obj = clade2_res$plot,
  filename = file.path(viz_dir, "clade2_top10_ST_distribution.pdf"),
  width = plot_width,
  height = plot_height
)

cat("Done. Top 10 ST distribution figures written to:\n", viz_dir, "\n", sep = "")
