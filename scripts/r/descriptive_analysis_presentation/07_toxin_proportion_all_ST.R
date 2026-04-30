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

metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
names(metadata) <- trimws(names(metadata))

metadata_clean <- metadata %>%
  mutate(
    ST = as.character(ST),
    ST = stringr::str_remove_all(ST, "\\*"),
    ST = trimws(ST),
    ST = na_if(ST, ""),
    toxin_status_updated = case_when(
      as.character(has_any_toxin) %in% c("TRUE", "True", "true", "1") ~ "Toxigenic",
      as.character(has_any_toxin) %in% c("FALSE", "False", "false", "0") ~ "Non-Toxigenic",
      TRUE ~ "Missing"
    )
  ) %>%
  filter(!is.na(ST), toxin_status_updated %in% c("Toxigenic", "Non-Toxigenic"))

st_counts <- metadata_clean %>%
  count(ST, name = "ST_n")

metadata_binned <- metadata_clean %>%
  left_join(st_counts, by = "ST") %>%
  mutate(
    binned_ST = if_else(ST_n > 1, paste0("ST ", ST), "Singletons")
  )

summary_counts <- metadata_binned %>%
  count(binned_ST, toxin_status_updated, name = "n") %>%
  complete(
    binned_ST,
    toxin_status_updated = c("Toxigenic", "Non-Toxigenic"),
    fill = list(n = 0)
  ) %>%
  group_by(binned_ST) %>%
  mutate(
    total_n = sum(n),
    proportion = if_else(total_n > 0, n / total_n, 0)
  ) %>%
  ungroup()

st_levels <- summary_counts %>%
  distinct(binned_ST, total_n) %>%
  arrange(total_n, binned_ST) %>%
  pull(binned_ST)

plot_obj <- ggplot(
  summary_counts,
  aes(
    y = factor(binned_ST, levels = st_levels),
    x = n,
    fill = toxin_status_updated
  )
) +
  geom_col(width = 0.75) +
  geom_text(
    data = summary_counts %>% distinct(binned_ST, total_n),
    aes(
      y = factor(binned_ST, levels = st_levels),
      x = total_n,
      label = paste0("n=", total_n)
    ),
    inherit.aes = FALSE,
    hjust = -0.1,
    size = 3.8
  ) +
  scale_fill_manual(
    values = c("Toxigenic" = "#e15759", "Non-Toxigenic" = "#4e79a7"),
    name = "Toxin status"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.14))
  ) +
  labs(
    title = "Toxigenic vs Non-Toxigenic Distribution Across All ST",
    subtitle = "Singleton STs are combined into one category",
    x = "Number of isolates",
    y = "Sequence type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

write.csv(
  summary_counts %>% arrange(desc(total_n), binned_ST, toxin_status_updated),
  file.path(desc_dir, "all_ST_toxin_proportion_summary.csv"),
  row.names = FALSE
)

ggsave(
  filename = file.path(viz_dir, "all_ST_toxin_proportion.png"),
  plot = plot_obj,
  width = 12,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)

save_plot_pdf(
  plot_obj = plot_obj,
  filename = file.path(viz_dir, "all_ST_toxin_proportion.pdf"),
  width = 12,
  height = 10
)

cat("Done. Toxin proportion plot written to:\n", viz_dir, "\n", sep = "")
