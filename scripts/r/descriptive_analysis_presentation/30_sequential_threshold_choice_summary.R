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

wilcox_path <- file.path(desc_dir, "sequential_threshold_comparison_wilcox_tests.csv")
enrichment_path <- file.path(desc_dir, "sequential_threshold_comparison_enrichment_tests.csv")

if (!dir.exists(viz_dir)) {
  dir.create(viz_dir, recursive = TRUE)
}

required_pkgs <- c("ggplot2", "dplyr", "readr", "tidyr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
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

wilcox_dat <- read_csv(wilcox_path, show_col_types = FALSE) %>%
  mutate(
    threshold_days = as.numeric(as.character(threshold_days)),
    exposure_type = factor(exposure_type, levels = c("Facility", "Floor", "Unit", "Room")),
    clade_group = factor(clade_group, levels = c("Clade 1", "Clade 2"))
  )

enrichment_dat <- read_csv(enrichment_path, show_col_types = FALSE) %>%
  filter(snv_rule %in% c("Min SNV <= 2", "Min SNV <= 5")) %>%
  mutate(
    threshold_days = as.numeric(as.character(threshold_days)),
    exposure_type = factor(exposure_type, levels = c("Facility", "Floor", "Unit", "Room")),
    clade_group = factor(clade_group, levels = c("Clade 1", "Clade 2")),
    pct_diff = yes_low_pct - no_low_pct
  )

choice_summary <- enrichment_dat %>%
  select(clade_group, exposure_type, threshold_days, snv_rule, yes_low_pct, no_low_pct, pct_diff, p_value_label) %>%
  left_join(
    wilcox_dat %>% select(clade_group, exposure_type, threshold_days, yes_n, yes_median_snv, no_median_snv, median_diff),
    by = c("clade_group", "exposure_type", "threshold_days")
  ) %>%
  arrange(clade_group, exposure_type, snv_rule, threshold_days)

write_csv(
  choice_summary,
  file.path(desc_dir, "sequential_threshold_choice_summary.csv")
)

highlight_df <- tibble(
  xmin = 12,
  xmax = 16,
  ymin = -Inf,
  ymax = Inf
)

p1_dat <- enrichment_dat %>% filter(snv_rule == "Min SNV <= 2")
p2_dat <- enrichment_dat %>% filter(snv_rule == "Min SNV <= 5")
p2_counts <- wilcox_dat

p1 <- ggplot(p1_dat, aes(x = threshold_days, y = pct_diff, color = exposure_type)) +
  geom_rect(
    data = highlight_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "#FCECC0",
    alpha = 0.5
  ) +
  geom_line(linewidth = 1) +
  geom_point(aes(size = yes_low_pct), alpha = 0.9) +
  facet_wrap(~clade_group, ncol = 2) +
  scale_color_manual(values = c("Facility" = "#4E79A7", "Floor" = "#59A14F", "Unit" = "#F28E2B", "Room" = "#E15759")) +
  scale_size_continuous(name = "Sequential yes\n% Min SNV <= 2", range = c(2, 6)) +
  scale_x_continuous(breaks = c(7, 14, 30, 60)) +
  labs(
    title = "Why 14 days is a reasonable sequential threshold",
    subtitle = "Panel A shows enrichment difference for Min SNV <= 2 (Sequential yes minus Sequential no); shaded band marks 14 days",
    x = "Sequential threshold (days)",
    y = "Percent difference"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

p2 <- ggplot(p2_counts, aes(x = threshold_days, y = yes_n, color = exposure_type)) +
  geom_rect(
    data = highlight_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "#FCECC0",
    alpha = 0.5
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.8) +
  facet_wrap(~clade_group, ncol = 2, scales = "free_y") +
  scale_color_manual(values = c("Facility" = "#4E79A7", "Floor" = "#59A14F", "Unit" = "#F28E2B", "Room" = "#E15759")) +
  scale_x_continuous(breaks = c(7, 14, 30, 60)) +
  labs(
    title = "Sequential sample size grows quickly beyond 14 days",
    subtitle = "Panel B shows the number of sequential yes patient-pairs; larger thresholds add more pairs but do not clearly improve the signal",
    x = "Sequential threshold (days)",
    y = "Sequential yes pair count",
    color = "Location level"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(viz_dir, "sequential_threshold_choice_summary_panelA.png"),
  plot = p1,
  width = 12,
  height = 6.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

save_plot_pdf(
  plot_obj = p1,
  filename = file.path(viz_dir, "sequential_threshold_choice_summary_panelA.pdf"),
  width = 12,
  height = 6.5
)

ggsave(
  filename = file.path(viz_dir, "sequential_threshold_choice_summary_panelB.png"),
  plot = p2,
  width = 12,
  height = 6.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

save_plot_pdf(
  plot_obj = p2,
  filename = file.path(viz_dir, "sequential_threshold_choice_summary_panelB.pdf"),
  width = 12,
  height = 6.5
)

writeLines(
  c(
    "Interpretation of sequential-threshold choice:",
    "Clade 1 does not show a strong low-SNV enrichment trend across thresholds, so threshold choice is not driven by Clade 1.",
    "Clade 2 shows the clearest and most interpretable signal at the Unit level around 14 days.",
    "At Unit, Min SNV <= 2 enrichment is strongest and statistically significant at 14 days (1.52% vs 0.49%, p = 0.028).",
    "At Unit, Min SNV <= 5 is also strongest and significant at 14 days (1.82% vs 0.74%, p = 0.041).",
    "A 7-day threshold is more restrictive but gives fewer sequential pairs.",
    "Thirty and 60-day thresholds add many more sequential pairs but dilute the low-SNV enrichment signal.",
    "For that reason, 14 days is a practical compromise between sensitivity and specificity."
  ),
  con = file.path(desc_dir, "sequential_threshold_choice_summary_notes.txt")
)

cat("Done. Sequential threshold choice summary written to:\n", viz_dir, "\n", sep = "")
