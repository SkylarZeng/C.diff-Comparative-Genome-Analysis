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

clade1_pairs_path <- file.path(desc_dir, "clade1_post_noncore_snv_long_with_location_overlap.csv")
clade2_pairs_path <- file.path(desc_dir, "clade2_post_noncore_snv_long_with_location_overlap.csv")

if (!dir.exists(desc_dir)) {
  dir.create(desc_dir, recursive = TRUE)
}

required_pkgs <- c("tidyverse")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(tidyverse)
})

safe_wilcox <- function(x_yes, x_no) {
  if (length(x_yes) == 0 || length(x_no) == 0) {
    return(NA_real_)
  }
  suppressWarnings(wilcox.test(x_yes, x_no)$p.value)
}

safe_fisher <- function(yes_low, yes_high, no_low, no_high) {
  mat <- matrix(c(yes_low, yes_high, no_low, no_high), nrow = 2, byrow = TRUE)
  if (any(rowSums(mat) == 0) || any(colSums(mat) == 0)) {
    return(NA_real_)
  }
  fisher.test(mat)$p.value
}

format_p <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x < 0.001 ~ "<0.001",
    TRUE ~ formatC(x, format = "f", digits = 3)
  )
}

exposure_vars <- tribble(
  ~exposure_class, ~exposure_type, ~flag_col,
  "Overlap", "Facility", "facility_overlap",
  "Overlap", "Floor", "floor_overlap",
  "Overlap", "Unit", "unit_overlap",
  "Overlap", "Room", "room_overlap",
  "Sequential_14d", "Facility", "facility_sequential_14d",
  "Sequential_14d", "Floor", "floor_sequential_14d",
  "Sequential_14d", "Unit", "unit_sequential_14d",
  "Sequential_14d", "Room", "room_sequential_14d",
  "Sequential_14d", "Bed", "bed_sequential_14d"
)

prepare_pairs <- function(path, clade_label) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE) %>%
    mutate(
      snv = suppressWarnings(as.numeric(snv)),
      clade_group = clade_label,
      same_patient = trimws(as.character(same_patient))
    ) %>%
    filter(
      same_patient == "Between patient",
      !is.na(snv),
      !is.na(patient_a),
      patient_a != "",
      !is.na(patient_b),
      patient_b != ""
    )
}

collapse_to_min_pair <- function(dat) {
  dat %>%
    group_by(clade_group, patient_a, patient_b) %>%
    summarise(
      min_snv = min(snv, na.rm = TRUE),
      genome_pair_n = n(),
      facility_overlap = any(facility_overlap %in% TRUE),
      floor_overlap = any(floor_overlap %in% TRUE),
      unit_overlap = any(unit_overlap %in% TRUE),
      room_overlap = any(room_overlap %in% TRUE),
      bed_overlap = any(bed_overlap %in% TRUE),
      facility_sequential_14d = any(facility_sequential_14d %in% TRUE),
      floor_sequential_14d = any(floor_sequential_14d %in% TRUE),
      unit_sequential_14d = any(unit_sequential_14d %in% TRUE),
      room_sequential_14d = any(room_sequential_14d %in% TRUE),
      bed_sequential_14d = any(bed_sequential_14d %in% TRUE),
      .groups = "drop"
    )
}

run_continuous_tests <- function(dat, analysis_label, snv_filter_label) {
  map_dfr(seq_len(nrow(exposure_vars)), function(i) {
    flag_col <- exposure_vars$flag_col[i]
    yes_vals <- dat %>% filter(.data[[flag_col]] %in% TRUE) %>% pull(min_snv)
    no_vals <- dat %>% filter(.data[[flag_col]] %in% FALSE) %>% pull(min_snv)
    p_val <- safe_wilcox(yes_vals, no_vals)

    tibble(
      analysis = analysis_label,
      snv_filter = snv_filter_label,
      exposure_class = exposure_vars$exposure_class[i],
      exposure_type = exposure_vars$exposure_type[i],
      yes_n = length(yes_vals),
      no_n = length(no_vals),
      yes_median_snv = ifelse(length(yes_vals) > 0, median(yes_vals), NA_real_),
      no_median_snv = ifelse(length(no_vals) > 0, median(no_vals), NA_real_),
      yes_iqr_snv = ifelse(length(yes_vals) > 0, IQR(yes_vals), NA_real_),
      no_iqr_snv = ifelse(length(no_vals) > 0, IQR(no_vals), NA_real_),
      wilcox_p_value = p_val,
      p_value_label = format_p(p_val)
    )
  })
}

run_continuous_tests_pairwise <- function(dat, analysis_label, snv_filter_label) {
  map_dfr(seq_len(nrow(exposure_vars)), function(i) {
    flag_col <- exposure_vars$flag_col[i]
    yes_vals <- dat %>% filter(.data[[flag_col]] %in% TRUE) %>% pull(snv)
    no_vals <- dat %>% filter(.data[[flag_col]] %in% FALSE) %>% pull(snv)
    p_val <- safe_wilcox(yes_vals, no_vals)

    tibble(
      analysis = analysis_label,
      snv_filter = snv_filter_label,
      exposure_class = exposure_vars$exposure_class[i],
      exposure_type = exposure_vars$exposure_type[i],
      yes_n = length(yes_vals),
      no_n = length(no_vals),
      yes_median_snv = ifelse(length(yes_vals) > 0, median(yes_vals), NA_real_),
      no_median_snv = ifelse(length(no_vals) > 0, median(no_vals), NA_real_),
      yes_iqr_snv = ifelse(length(yes_vals) > 0, IQR(yes_vals), NA_real_),
      no_iqr_snv = ifelse(length(no_vals) > 0, IQR(no_vals), NA_real_),
      wilcox_p_value = p_val,
      p_value_label = format_p(p_val)
    )
  })
}

run_enrichment_tests_pairwise <- function(dat, analysis_label) {
  thresholds <- c(2, 5, 10)

  crossing(exposure_vars, threshold = thresholds) %>%
    mutate(row_id = row_number()) %>%
    split(.$row_id) %>%
    map_dfr(function(x) {
      flag_col <- x$flag_col[[1]]
      threshold <- x$threshold[[1]]

      yes_dat <- dat %>% filter(.data[[flag_col]] %in% TRUE)
      no_dat <- dat %>% filter(.data[[flag_col]] %in% FALSE)

      yes_low <- sum(yes_dat$snv <= threshold, na.rm = TRUE)
      yes_high <- sum(yes_dat$snv > threshold, na.rm = TRUE)
      no_low <- sum(no_dat$snv <= threshold, na.rm = TRUE)
      no_high <- sum(no_dat$snv > threshold, na.rm = TRUE)
      p_val <- safe_fisher(yes_low, yes_high, no_low, no_high)

      tibble(
        analysis = analysis_label,
        exposure_class = x$exposure_class[[1]],
        exposure_type = x$exposure_type[[1]],
        threshold_rule = paste0("SNV <= ", threshold),
        yes_low_n = yes_low,
        yes_high_n = yes_high,
        no_low_n = no_low,
        no_high_n = no_high,
        yes_low_pct = ifelse((yes_low + yes_high) > 0, 100 * yes_low / (yes_low + yes_high), NA_real_),
        no_low_pct = ifelse((no_low + no_high) > 0, 100 * no_low / (no_low + no_high), NA_real_),
        fisher_p_value = p_val,
        p_value_label = format_p(p_val)
      )
    })
}

clade1_raw <- prepare_pairs(clade1_pairs_path, "Clade 1")
clade2_raw <- prepare_pairs(clade2_pairs_path, "Clade 2")
all_raw <- bind_rows(clade1_raw, clade2_raw)

clade1_min <- collapse_to_min_pair(clade1_raw)
clade2_min <- collapse_to_min_pair(clade2_raw)
all_min <- bind_rows(clade1_min, clade2_min)

write.csv(
  clade1_min,
  file.path(desc_dir, "clade1_between_patient_min_snv_by_patient_pair.csv"),
  row.names = FALSE
)
write.csv(
  clade2_min,
  file.path(desc_dir, "clade2_between_patient_min_snv_by_patient_pair.csv"),
  row.names = FALSE
)
write.csv(
  all_min,
  file.path(desc_dir, "all_between_patient_min_snv_by_patient_pair.csv"),
  row.names = FALSE
)

primary_continuous <- bind_rows(
  run_continuous_tests(all_min, "Overall", "No SNV restriction"),
  run_continuous_tests(clade1_min, "Clade 1", "No SNV restriction"),
  run_continuous_tests(clade2_min, "Clade 2", "No SNV restriction")
)

sensitivity_continuous <- bind_rows(
  run_continuous_tests_pairwise(all_raw %>% filter(snv <= 50), "Overall", "SNV <= 50"),
  run_continuous_tests_pairwise(clade1_raw %>% filter(snv <= 50), "Clade 1", "SNV <= 50"),
  run_continuous_tests_pairwise(clade2_raw %>% filter(snv <= 50), "Clade 2", "SNV <= 50")
)

secondary_enrichment <- bind_rows(
  run_enrichment_tests_pairwise(clade1_raw, "Clade 1"),
  run_enrichment_tests_pairwise(clade2_raw, "Clade 2")
)

write.csv(
  primary_continuous,
  file.path(desc_dir, "primary_min_snv_patient_pair_wilcox_tests.csv"),
  row.names = FALSE
)
write.csv(
  sensitivity_continuous,
  file.path(desc_dir, "sensitivity_pairwise_snv_le_50_wilcox_tests.csv"),
  row.names = FALSE
)
write.csv(
  secondary_enrichment,
  file.path(desc_dir, "secondary_pairwise_snv_threshold_enrichment_fisher_tests.csv"),
  row.names = FALSE
)

writeLines(
  c(
    "Primary analysis: minimum SNV per patient pair, comparing Yes vs No with Wilcoxon rank-sum tests.",
    "Secondary analysis: enrichment of pairwise SNV rows with SNV <= 2, <= 5, and <= 10 using Fisher's exact tests.",
    "Sensitivity analysis: rerun the continuous pairwise SNV comparison after restricting to SNV <= 50.",
    "All analyses use between-patient pairs only.",
    "Exposure flags are inherited from the location overlap/sequential pair annotations.",
    "",
    "Interpretation:",
    "The secondary threshold-enrichment analysis is the most interpretable primary result because it directly asks whether exposed pairs are enriched for very low SNV values that are most relevant to recent transmission.",
    "The primary minimum-SNV-per-patient-pair analysis is useful as a patient-level summary, but it does not consistently show that exposed pairs have lower SNV values overall.",
    "The pairwise SNV <= 50 sensitivity analysis is best treated as supportive rather than primary, because it conditions on a restricted SNV range and answers a narrower question.",
    "In these data, the clearest enrichment signal appears in Clade 2, where several overlap categories show higher percentages of low-SNV pairs at thresholds <=2 and <=5."
  ),
  con = file.path(desc_dir, "min_snv_patient_pair_overlap_analyses_notes.txt")
)

cat("Done. Minimum-SNV patient-pair analysis outputs written to:\n", desc_dir, "\n", sep = "")
