args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_all[grep(file_arg, args_all)])

if (length(script_path) == 0) {
  script_dir <- getwd()
} else {
  script_dir <- dirname(normalizePath(script_path))
}

base_dir <- normalizePath(file.path(script_dir, ".."))
cluster_dir <- normalizePath(file.path(base_dir, ".."))
output_dir <- file.path(base_dir, "outputs")
desc_dir <- file.path(output_dir, "descriptive_statistics")

clade1_snv_path <- file.path(cluster_dir, "SNV", "clade1_post_noncore_snv_long.csv")
clade2_snv_path <- file.path(cluster_dir, "SNV", "clade2_post_noncore_snv_long.csv")
metadata_path <- file.path(desc_dir, "metadata_case_unformed_clade1_clade2.csv")

required_pkgs <- c("tidyverse")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(tidyverse)
})

parse_date_safely <- function(x) {
  as.Date(sub("T.*$", "", as.character(x)))
}

metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
names(metadata) <- trimws(names(metadata))

lookup <- metadata %>%
  transmute(
    genome_id = as.character(genome_id),
    patient_id = trimws(as.character(patient_id)),
    collection_date = parse_date_safely(collection_date)
  ) %>%
  mutate(patient_id = na_if(patient_id, "")) %>%
  distinct(genome_id, .keep_all = TRUE)

prepare_high_snv_table <- function(path, clade_label) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE) %>%
    rename(snv = snvs) %>%
    mutate(
      id1 = as.character(id1),
      id2 = as.character(id2),
      snv = suppressWarnings(as.numeric(snv))
    ) %>%
    left_join(lookup, by = c("id1" = "genome_id")) %>%
    rename(
      patient_id_1 = patient_id,
      collection_date_1 = collection_date
    ) %>%
    left_join(lookup, by = c("id2" = "genome_id")) %>%
    rename(
      patient_id_2 = patient_id,
      collection_date_2 = collection_date
    ) %>%
    mutate(
      same_patient = case_when(
        !is.na(patient_id_1) & !is.na(patient_id_2) & patient_id_1 == patient_id_2 ~ "Within patient",
        !is.na(patient_id_1) & !is.na(patient_id_2) & patient_id_1 != patient_id_2 ~ "Between patient",
        TRUE ~ NA_character_
      ),
      patient_id = if_else(same_patient == "Within patient", patient_id_1, NA_character_),
      days_between_collection = abs(as.numeric(collection_date_2 - collection_date_1)),
      clade = clade_label
    ) %>%
    filter(same_patient == "Within patient", !is.na(snv), snv > 10) %>%
    arrange(desc(snv), desc(days_between_collection))
}

clade1_high <- prepare_high_snv_table(clade1_snv_path, "Clade 1")
clade2_high <- prepare_high_snv_table(clade2_snv_path, "Clade 2")

high_snv_all <- bind_rows(clade1_high, clade2_high) %>%
  select(
    clade,
    patient_id,
    id1,
    id2,
    collection_date_1,
    collection_date_2,
    days_between_collection,
    snv
  )

gap_summary <- high_snv_all %>%
  mutate(
    time_gap_group = case_when(
      is.na(days_between_collection) ~ "Missing collection date",
      days_between_collection <= 30 ~ "<=30 days",
      days_between_collection <= 90 ~ "31-90 days",
      days_between_collection <= 180 ~ "91-180 days",
      days_between_collection <= 365 ~ "181-365 days",
      TRUE ~ ">365 days"
    )
  ) %>%
  group_by(clade, time_gap_group) %>%
  summarise(
    pair_n = n(),
    median_snv = median(snv, na.rm = TRUE),
    median_days_between = median(days_between_collection, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(clade, time_gap_group)

overall_summary <- high_snv_all %>%
  group_by(clade) %>%
  summarise(
    pair_n = n(),
    median_snv = median(snv, na.rm = TRUE),
    median_days_between = median(days_between_collection, na.rm = TRUE),
    min_days_between = min(days_between_collection, na.rm = TRUE),
    max_days_between = max(days_between_collection, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  high_snv_all,
  file.path(desc_dir, "within_patient_pairs_snv_gt10_with_collection_gap.csv"),
  row.names = FALSE
)

write.csv(
  gap_summary,
  file.path(desc_dir, "within_patient_pairs_snv_gt10_gap_summary.csv"),
  row.names = FALSE
)

write.csv(
  overall_summary,
  file.path(desc_dir, "within_patient_pairs_snv_gt10_overall_summary.csv"),
  row.names = FALSE
)

writeLines(
  c(
    "Suggested interpretation:",
    "Within-patient pairs with SNV > 10 may reflect longer time between collections, allowing additional within-host evolution or reinfection by a distinct strain.",
    "Use the pair-level table and the gap summary to show whether high-SNV within-patient pairs are enriched at longer collection intervals."
  ),
  con = file.path(desc_dir, "within_patient_pairs_snv_gt10_notes.txt")
)

cat("Done. High-SNV within-patient explanation tables written to:\n", desc_dir, "\n", sep = "")
