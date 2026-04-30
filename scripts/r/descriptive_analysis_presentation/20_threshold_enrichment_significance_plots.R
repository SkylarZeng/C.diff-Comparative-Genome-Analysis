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

enrichment_path <- file.path(desc_dir, "secondary_pairwise_snv_threshold_enrichment_fisher_tests.csv")

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

make_star_label <- function(p_label) {
  case_when(
    is.na(p_label) ~ "NA",
    p_label == "<0.001" ~ "***",
    suppressWarnings(as.numeric(p_label)) < 0.01 ~ "**",
    suppressWarnings(as.numeric(p_label)) < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

enrichment <- read.csv(enrichment_path, check.names = FALSE, stringsAsFactors = FALSE) %>%
  filter(
    analysis %in% c("Clade 1", "Clade 2"),
    (exposure_class == "Overlap" & exposure_type %in% c("Facility", "Floor", "Unit", "Room")) |
      (exposure_class == "Sequential_14d" & exposure_type %in% c("Facility", "Floor", "Unit", "Room", "Bed"))
  ) %>%
  mutate(
    analysis = factor(analysis, levels = c("Clade 1", "Clade 2")),
    exposure_class = factor(exposure_class, levels = c("Overlap", "Sequential_14d")),
    exposure_type = factor(exposure_type, levels = c("Facility", "Floor", "Unit", "Room", "Bed")),
    threshold_rule = factor(threshold_rule, levels = c("SNV <= 2", "SNV <= 5", "SNV <= 10")),
    star_label = make_star_label(p_value_label)
  )

plot_dat <- bind_rows(
  enrichment %>%
    transmute(
      analysis,
      exposure_class,
      exposure_type,
      threshold_rule,
      group = "Yes",
      pct = yes_low_pct,
      p_value_label,
      star_label
    ),
  enrichment %>%
    transmute(
      analysis,
      exposure_class,
      exposure_type,
      threshold_rule,
      group = "No",
      pct = no_low_pct,
      p_value_label,
      star_label
    )
) %>%
  mutate(group = factor(group, levels = c("Yes", "No")))

ann_dat <- plot_dat %>%
  group_by(analysis, exposure_class, exposure_type, threshold_rule) %>%
  summarise(
    y_pos = max(pct, na.rm = TRUE) * 1.15 + 0.2,
    star_label = first(star_label),
    .groups = "drop"
  )

make_plot <- function(exposure_class_value, title_text, file_stub) {
  dat <- plot_dat %>% filter(exposure_class == exposure_class_value)
  ann <- ann_dat %>% filter(exposure_class == exposure_class_value)

  p <- ggplot(dat, aes(x = threshold_rule, y = pct, fill = group)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    geom_text(
      aes(label = sprintf("%.1f", pct)),
      position = position_dodge(width = 0.75),
      vjust = -0.2,
      size = 3
    ) +
    geom_text(
      data = ann,
      aes(x = threshold_rule, y = y_pos, label = star_label),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 5
    ) +
    facet_grid(analysis ~ exposure_type, scales = "free_y") +
    scale_fill_manual(values = c("Yes" = "#008ECE", "No" = "#ECA0B2"), name = NULL) +
    labs(
      title = title_text,
      subtitle = "Bars show percent of pairwise SNV rows below threshold; stars mark Fisher test significance",
      x = NULL,
      y = "Percent below threshold"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(viz_dir, paste0(file_stub, ".png")),
    plot = p,
    width = 16,
    height = 10,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  save_plot_pdf(
    plot_obj = p,
    filename = file.path(viz_dir, paste0(file_stub, ".pdf")),
    width = 16,
    height = 10
  )
}

make_plot(
  exposure_class_value = "Overlap",
  title_text = "Threshold enrichment for low SNV values: direct overlap",
  file_stub = "threshold_enrichment_overlap_significance"
)

make_plot(
  exposure_class_value = "Sequential_14d",
  title_text = "Threshold enrichment for low SNV values: sequential 14-day exposure",
  file_stub = "threshold_enrichment_sequential_14d_significance"
)

cat("Done. Threshold enrichment significance plots written to:\n", viz_dir, "\n", sep = "")
