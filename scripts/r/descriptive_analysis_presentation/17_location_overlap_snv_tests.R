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

sequential_suffix <- "14d"
sequential_class_label <- "Sequential_14d"

safe_wilcox <- function(x_yes, x_no) {
  if (length(x_yes) == 0 || length(x_no) == 0) {
    return(NA_real_)
  }

  suppressWarnings(
    wilcox.test(x_yes, x_no)$p.value
  )
}

prepare_pairs <- function(path, clade_label) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE) %>%
    mutate(
      snv = suppressWarnings(as.numeric(snv)),
      clade_group = clade_label,
      same_patient = trimws(as.character(same_patient))
    ) %>%
    filter(same_patient == "Between patient", !is.na(snv))
}

run_tests_for_flag <- function(dat, flag_col, exposure_class, exposure_type) {
  yes_vals <- dat %>%
    filter(.data[[flag_col]] %in% TRUE) %>%
    pull(snv)

  no_vals <- dat %>%
    filter(.data[[flag_col]] %in% FALSE) %>%
    pull(snv)

  tibble(
    exposure_class = exposure_class,
    exposure_type = exposure_type,
    yes_n = length(yes_vals),
    no_n = length(no_vals),
    yes_median_snv = ifelse(length(yes_vals) > 0, median(yes_vals), NA_real_),
    no_median_snv = ifelse(length(no_vals) > 0, median(no_vals), NA_real_),
    yes_iqr_snv = ifelse(length(yes_vals) > 0, IQR(yes_vals), NA_real_),
    no_iqr_snv = ifelse(length(no_vals) > 0, IQR(no_vals), NA_real_),
    wilcox_p_value = safe_wilcox(yes_vals, no_vals)
  )
}

run_all_tests <- function(dat, clade_label) {
  bind_rows(
    run_tests_for_flag(dat, "facility_overlap", "Overlap", "Facility"),
    run_tests_for_flag(dat, "floor_overlap", "Overlap", "Floor"),
    run_tests_for_flag(dat, "unit_overlap", "Overlap", "Unit"),
    run_tests_for_flag(dat, "room_overlap", "Overlap", "Room"),
    run_tests_for_flag(dat, "facility_sequential_14d", sequential_class_label, "Facility"),
    run_tests_for_flag(dat, "floor_sequential_14d", sequential_class_label, "Floor"),
    run_tests_for_flag(dat, "unit_sequential_14d", sequential_class_label, "Unit"),
    run_tests_for_flag(dat, "room_sequential_14d", sequential_class_label, "Room"),
    run_tests_for_flag(dat, "bed_sequential_14d", sequential_class_label, "Bed")
  ) %>%
    mutate(clade_group = clade_label) %>%
    select(clade_group, everything())
}

clade1_pairs <- prepare_pairs(clade1_pairs_path, "Clade 1")
clade2_pairs <- prepare_pairs(clade2_pairs_path, "Clade 2")
all_pairs <- bind_rows(clade1_pairs, clade2_pairs)

test_results <- bind_rows(
  run_all_tests(all_pairs, "Overall"),
  run_all_tests(clade1_pairs, "Clade 1"),
  run_all_tests(clade2_pairs, "Clade 2")
) %>%
  mutate(
    p_value_label = case_when(
      is.na(wilcox_p_value) ~ NA_character_,
      wilcox_p_value < 0.001 ~ "<0.001",
      TRUE ~ formatC(wilcox_p_value, format = "f", digits = 3)
    )
  )

write.csv(
  test_results,
  file.path(desc_dir, "location_overlap_and_sequential_snv_wilcox_tests.csv"),
  row.names = FALSE
)

writeLines(
  c(
    "Tests use between-patient SNV pairs only.",
    "For each location level, SNV distances were compared between Yes vs No using the Wilcoxon rank-sum test.",
    "Overlap = direct temporal overlap in the same location.",
    "Sequential_14d = same location, no direct overlap, and later stay begins within 14 days after earlier stay ends.",
    "Special units were excluded upstream in the annotated pair tables."
  ),
  con = file.path(desc_dir, "location_overlap_and_sequential_snv_wilcox_tests_notes.txt")
)

cat("Done. Location overlap/sequential SNV test outputs written to:\n", desc_dir, "\n", sep = "")
