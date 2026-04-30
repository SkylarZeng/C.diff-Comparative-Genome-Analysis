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

input_path <- file.path(desc_dir, "within_patient_pairs_st_and_time_for_snv_tests.csv")

if (!dir.exists(viz_dir)) {
  dir.create(viz_dir, recursive = TRUE)
}

required_pkgs <- c("ggplot2", "dplyr", "readr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
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

plot_dat <- read_csv(input_path, show_col_types = FALSE) %>%
  mutate(
    clade = factor(clade, levels = c("Clade 1", "Clade 2")),
    snv = suppressWarnings(as.numeric(snv)),
    days_between_collection = suppressWarnings(as.numeric(days_between_collection)),
    st_same = factor(st_same, levels = c("ST same", "ST different"))
  ) %>%
  filter(!is.na(clade), !is.na(snv), !is.na(days_between_collection), !is.na(st_same))

p <- ggplot(plot_dat, aes(x = days_between_collection, y = snv, color = st_same)) +
  geom_point(alpha = 0.7, size = 2.6) +
  geom_hline(yintercept = 10, linetype = "dashed", linewidth = 0.8, color = "gray40") +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  facet_wrap(~clade, ncol = 2, scales = "free_y") +
  scale_color_manual(
    values = c("ST same" = "#4E79A7", "ST different" = "#E15759"),
    name = NULL
  ) +
  labs(
    title = "Within-patient SNV versus days between collections",
    subtitle = "Dashed line marks SNV = 10",
    x = "Days between collections",
    y = "SNV distance"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(viz_dir, "within_patient_snv_vs_collection_days.png"),
  plot = p,
  width = 11,
  height = 6.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

save_plot_pdf(
  plot_obj = p,
  filename = file.path(viz_dir, "within_patient_snv_vs_collection_days.pdf"),
  width = 11,
  height = 6.5
)

cat("Done. Within-patient SNV vs collection-days plot written to:\n", viz_dir, "\n", sep = "")
