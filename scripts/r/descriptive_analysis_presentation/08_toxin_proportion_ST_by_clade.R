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

plot_width <- 6
plot_height <- 6

metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
names(metadata) <- trimws(names(metadata))

metadata_clean <- metadata %>%
  mutate(
    ST = as.character(ST),
    ST = stringr::str_remove_all(ST, "\\*"),
    ST = trimws(ST),
    ST = na_if(ST, ""),
    clade_group = trimws(as.character(clade_group)),
    toxin_status_updated = case_when(
      as.character(has_any_toxin) %in% c("TRUE", "True", "true", "1") ~ "Toxigenic",
      as.character(has_any_toxin) %in% c("FALSE", "False", "false", "0") ~ "Non-Toxigenic",
      TRUE ~ "Missing"
    )
  ) %>%
  filter(
    !is.na(ST),
    clade_group %in% c("Clade 1", "Clade 2"),
    toxin_status_updated %in% c("Toxigenic", "Non-Toxigenic")
  )

make_st_toxin_plot <- function(dat, clade_label) {
  dat_clade <- dat %>%
    filter(clade_group == clade_label)

  if (nrow(dat_clade) == 0) {
    stop("No rows found for ", clade_label)
  }

  st_counts <- dat_clade %>%
    count(ST, name = "ST_n")

  dat_binned <- dat_clade %>%
    left_join(st_counts, by = "ST") %>%
    mutate(binned_ST = if_else(ST_n > 1, paste0("ST ", ST), "Singletons"))

  summary_binned <- dat_binned %>%
    count(binned_ST, toxin_status_updated, name = "n") %>%
    complete(
      binned_ST,
      toxin_status_updated = c("Toxigenic", "Non-Toxigenic"),
      fill = list(n = 0)
    ) %>%
    group_by(binned_ST) %>%
    mutate(total_n = sum(n)) %>%
    ungroup()

  st_levels <- summary_binned %>%
    distinct(binned_ST, total_n) %>%
    arrange(total_n, binned_ST) %>%
    pull(binned_ST)

  plot_obj <- ggplot(
    summary_binned,
    aes(
      y = factor(binned_ST, levels = st_levels),
      x = n,
      fill = toxin_status_updated
    )
  ) +
    geom_col(width = 0.75) +
    geom_text(
      data = summary_binned %>% distinct(binned_ST, total_n),
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
      values = c("Toxigenic" = "#ECA0B2", "Non-Toxigenic" = "#54BFB7"),
      name = "Toxin status"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = paste0(clade_label, ": ST Distribution by Toxin Status"),
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

  list(
    plot = plot_obj,
    summary = summary_binned %>% arrange(desc(total_n), binned_ST, toxin_status_updated)
  )
}

clade1_res <- make_st_toxin_plot(metadata_clean, "Clade 1")
clade2_res <- make_st_toxin_plot(metadata_clean, "Clade 2")

write.csv(
  clade1_res$summary,
  file.path(desc_dir, "clade1_ST_toxin_proportion_summary.csv"),
  row.names = FALSE
)
write.csv(
  clade2_res$summary,
  file.path(desc_dir, "clade2_ST_toxin_proportion_summary.csv"),
  row.names = FALSE
)

ggsave(
  filename = file.path(viz_dir, "clade1_ST_toxin_proportion.png"),
  plot = clade1_res$plot,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 300,
  bg = "white"
)

save_plot_pdf(
  plot_obj = clade1_res$plot,
  filename = file.path(viz_dir, "clade1_ST_toxin_proportion.pdf"),
  width = plot_width,
  height = plot_height
)

ggsave(
  filename = file.path(viz_dir, "clade2_ST_toxin_proportion.png"),
  plot = clade2_res$plot,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 300,
  bg = "white"
)

save_plot_pdf(
  plot_obj = clade2_res$plot,
  filename = file.path(viz_dir, "clade2_ST_toxin_proportion.pdf"),
  width = plot_width,
  height = plot_height
)

cat("Done. Clade-specific ST toxin proportion plots written to:\n", viz_dir, "\n", sep = "")
