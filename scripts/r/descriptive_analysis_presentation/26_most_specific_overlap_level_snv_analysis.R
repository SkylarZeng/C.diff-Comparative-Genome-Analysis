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

format_p <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x < 0.001 ~ "<0.001",
    TRUE ~ formatC(x, format = "f", digits = 3)
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

safe_pairwise_wilcox <- function(dat) {
  levels_present <- unique(as.character(dat$most_specific_overlap_level))
  if (length(levels_present) < 2) {
    return(tibble())
  }

  combn(sort(levels_present), 2, simplify = FALSE) %>%
    map_dfr(function(x) {
      dat_sub <- dat %>% filter(most_specific_overlap_level %in% x)
      vals_1 <- dat_sub %>% filter(most_specific_overlap_level == x[1]) %>% pull(snv)
      vals_2 <- dat_sub %>% filter(most_specific_overlap_level == x[2]) %>% pull(snv)
      p_val <- suppressWarnings(wilcox.test(vals_1, vals_2)$p.value)
      tibble(
        level_1 = x[1],
        level_2 = x[2],
        n_1 = length(vals_1),
        n_2 = length(vals_2),
        median_1 = median(vals_1),
        median_2 = median(vals_2),
        wilcox_p_value = p_val,
        p_value_label = format_p(p_val)
      )
    })
}

safe_fisher <- function(yes_low, yes_high, no_low, no_high) {
  mat <- matrix(c(yes_low, yes_high, no_low, no_high), nrow = 2, byrow = TRUE)
  if (any(rowSums(mat) == 0) || any(colSums(mat) == 0)) {
    return(NA_real_)
  }
  fisher.test(mat)$p.value
}

prepare_pairs <- function(path, clade_label) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE) %>%
    mutate(
      snv = suppressWarnings(as.numeric(snv)),
      clade_group = clade_label,
      same_patient = trimws(as.character(same_patient))
    ) %>%
    filter(same_patient == "Between patient", !is.na(snv)) %>%
    mutate(
      most_specific_overlap_level = case_when(
        room_overlap %in% TRUE ~ "Room",
        unit_overlap %in% TRUE ~ "Unit",
        floor_overlap %in% TRUE ~ "Floor",
        facility_overlap %in% TRUE ~ "Facility",
        TRUE ~ NA_character_
      ),
      most_specific_overlap_level = factor(
        most_specific_overlap_level,
        levels = c("Facility", "Floor", "Unit", "Room")
      )
    ) %>%
    filter(!is.na(most_specific_overlap_level))
}

clade1_pairs <- prepare_pairs(clade1_pairs_path, "Clade 1")
clade2_pairs <- prepare_pairs(clade2_pairs_path, "Clade 2")
all_pairs <- bind_rows(clade1_pairs, clade2_pairs)

overlap_level_summary <- all_pairs %>%
  group_by(clade_group, most_specific_overlap_level) %>%
  summarise(
    pairwise_snv_n = n(),
    median_snv = median(snv),
    iqr_snv = IQR(snv),
    min_snv = min(snv),
    max_snv = max(snv),
    .groups = "drop"
  )

kruskal_results <- bind_rows(
  all_pairs %>% mutate(analysis = "Overall"),
  clade1_pairs %>% mutate(analysis = "Clade 1"),
  clade2_pairs %>% mutate(analysis = "Clade 2")
) %>%
  group_by(analysis) %>%
  group_modify(~ {
    if (n_distinct(.x$most_specific_overlap_level) < 2) {
      tibble(kruskal_p_value = NA_real_, p_value_label = NA_character_)
    } else {
      p_val <- kruskal.test(snv ~ most_specific_overlap_level, data = .x)$p.value
      tibble(kruskal_p_value = p_val, p_value_label = format_p(p_val))
    }
  }) %>%
  ungroup()

pairwise_results <- bind_rows(
  safe_pairwise_wilcox(all_pairs) %>% mutate(analysis = "Overall"),
  safe_pairwise_wilcox(clade1_pairs) %>% mutate(analysis = "Clade 1"),
  safe_pairwise_wilcox(clade2_pairs) %>% mutate(analysis = "Clade 2")
) %>%
  select(analysis, everything())

thresholds <- c(2, 5, 10)
threshold_results <- crossing(
  clade_group = c("Clade 1", "Clade 2"),
  threshold = thresholds,
  level_1 = c("Facility", "Floor", "Unit"),
  level_2 = c("Floor", "Unit", "Room")
) %>%
  filter(
    match(level_1, c("Facility", "Floor", "Unit", "Room")) <
      match(level_2, c("Facility", "Floor", "Unit", "Room"))
  ) %>%
  rowwise() %>%
  mutate(
    dat = list(all_pairs %>% filter(clade_group == .data$clade_group)),
    vals_1 = list(dat %>% filter(most_specific_overlap_level == level_1) %>% pull(snv)),
    vals_2 = list(dat %>% filter(most_specific_overlap_level == level_2) %>% pull(snv)),
    level_1_low_n = sum(unlist(vals_1) <= threshold, na.rm = TRUE),
    level_1_high_n = sum(unlist(vals_1) > threshold, na.rm = TRUE),
    level_2_low_n = sum(unlist(vals_2) <= threshold, na.rm = TRUE),
    level_2_high_n = sum(unlist(vals_2) > threshold, na.rm = TRUE),
    level_1_low_pct = ifelse((level_1_low_n + level_1_high_n) > 0, 100 * level_1_low_n / (level_1_low_n + level_1_high_n), NA_real_),
    level_2_low_pct = ifelse((level_2_low_n + level_2_high_n) > 0, 100 * level_2_low_n / (level_2_low_n + level_2_high_n), NA_real_),
    fisher_p_value = safe_fisher(level_2_low_n, level_2_high_n, level_1_low_n, level_1_high_n),
    p_value_label = format_p(fisher_p_value)
  ) %>%
  ungroup() %>%
  transmute(
    clade_group,
    level_1,
    level_2,
    threshold_rule = paste0("SNV <= ", threshold),
    level_1_low_n,
    level_1_high_n,
    level_2_low_n,
    level_2_high_n,
    level_1_low_pct,
    level_2_low_pct,
    fisher_p_value,
    p_value_label
  )

write.csv(
  all_pairs,
  file.path(desc_dir, "most_specific_overlap_level_pairwise_snv_rows.csv"),
  row.names = FALSE
)
write.csv(
  overlap_level_summary,
  file.path(desc_dir, "most_specific_overlap_level_snv_summary.csv"),
  row.names = FALSE
)
write.csv(
  kruskal_results,
  file.path(desc_dir, "most_specific_overlap_level_kruskal_tests.csv"),
  row.names = FALSE
)
write.csv(
  pairwise_results,
  file.path(desc_dir, "most_specific_overlap_level_pairwise_wilcox_tests.csv"),
  row.names = FALSE
)
write.csv(
  threshold_results,
  file.path(desc_dir, "most_specific_overlap_level_threshold_comparisons.csv"),
  row.names = FALSE
)

ann_dat <- overlap_level_summary %>%
  group_by(clade_group) %>%
  summarise(y_pos = max(max_snv, na.rm = TRUE) * 1.15, .groups = "drop") %>%
  left_join(
    kruskal_results %>%
      filter(analysis %in% c("Clade 1", "Clade 2")) %>%
      transmute(clade_group = analysis, star_label = make_star_label(p_value_label)),
    by = "clade_group"
  )

boxplot_obj <- ggplot(
  all_pairs,
  aes(x = most_specific_overlap_level, y = snv, fill = most_specific_overlap_level)
) +
  geom_boxplot(
    width = 0.65,
    outlier.alpha = 0.08,
    outlier.size = 0.4
  ) +
  geom_text(
    data = ann_dat,
    aes(x = 2.5, y = y_pos, label = paste0("Kruskal-Wallis: ", star_label)),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 5
  ) +
  facet_wrap(~clade_group, ncol = 2, scales = "free_y") +
  scale_fill_manual(
    values = c(
      "Facility" = "#4E79A7",
      "Floor" = "#59A14F",
      "Unit" = "#F28E2B",
      "Room" = "#E15759"
    ),
    guide = "none"
  ) +
  scale_y_log10() +
  labs(
    title = "SNV by most specific overlap level",
    subtitle = "Between-patient pairwise SNV rows with direct overlap only",
    x = "Most specific overlap level",
    y = "SNV distance (log scale)"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(viz_dir, "most_specific_overlap_level_snv_boxplot.png"),
  plot = boxplot_obj,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = boxplot_obj,
  filename = file.path(viz_dir, "most_specific_overlap_level_snv_boxplot.pdf"),
  width = 12,
  height = 8
)

threshold_plot_dat <- threshold_results %>%
  filter(level_2 == "Room") %>%
  bind_rows(
    threshold_results %>%
      filter(level_2 == "Unit", level_1 %in% c("Facility", "Floor"))
  ) %>%
  mutate(
    comparison = paste(level_2, "vs", level_1),
    comparison = factor(comparison, levels = c("Room vs Facility", "Room vs Floor", "Room vs Unit", "Unit vs Facility", "Unit vs Floor")),
    threshold_rule = factor(threshold_rule, levels = c("SNV <= 2", "SNV <= 5", "SNV <= 10")),
    star_label = make_star_label(p_value_label)
  )

threshold_ann <- threshold_plot_dat %>%
  group_by(clade_group, comparison, threshold_rule) %>%
  summarise(
    y_pos = max(level_1_low_pct, level_2_low_pct, na.rm = TRUE) * 1.15 + 0.2,
    star_label = first(star_label),
    .groups = "drop"
  )

threshold_plot_obj <- bind_rows(
  threshold_plot_dat %>%
    transmute(clade_group, comparison, threshold_rule, group = "More specific", pct = level_2_low_pct),
  threshold_plot_dat %>%
    transmute(clade_group, comparison, threshold_rule, group = "Less specific", pct = level_1_low_pct)
) %>%
  ggplot(aes(x = threshold_rule, y = pct, fill = group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(
    aes(label = sprintf("%.1f", pct)),
    position = position_dodge(width = 0.75),
    vjust = -0.2,
    size = 3
  ) +
  geom_text(
    data = threshold_ann,
    aes(x = threshold_rule, y = y_pos, label = star_label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4.8
  ) +
  facet_grid(clade_group ~ comparison, scales = "free_y") +
  scale_fill_manual(values = c("More specific" = "#E15759", "Less specific" = "#4E79A7"), name = NULL) +
  labs(
    title = "Low-SNV enrichment by most specific overlap level",
    subtitle = "Room overlap compared with broader overlap levels",
    x = NULL,
    y = "Percent below threshold"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(viz_dir, "most_specific_overlap_level_threshold_plot.png"),
  plot = threshold_plot_obj,
  width = 15,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = threshold_plot_obj,
  filename = file.path(viz_dir, "most_specific_overlap_level_threshold_plot.pdf"),
  width = 15,
  height = 8
)

writeLines(
  c(
    "This analysis uses between-patient pairwise SNV rows with direct overlap only.",
    "Each row is assigned to the most specific direct overlap level in a nested hierarchy: Room > Unit > Floor > Facility.",
    "This avoids double-counting the same pair across multiple overlap levels in the main comparison.",
    "Kruskal-Wallis tests compare SNV distributions across overlap levels within each clade.",
    "Pairwise Wilcoxon tests provide level-by-level comparisons.",
    "Threshold comparisons show whether more specific overlap levels are enriched for low-SNV pairs."
  ),
  con = file.path(desc_dir, "most_specific_overlap_level_snv_analysis_notes.txt")
)

cat("Done. Most-specific overlap level SNV analysis outputs written to:\n", desc_dir, "\n", sep = "")
