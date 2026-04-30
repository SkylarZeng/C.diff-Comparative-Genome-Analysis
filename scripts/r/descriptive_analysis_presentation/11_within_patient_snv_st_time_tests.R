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

format_p_value <- function(p) {
  if (is.na(p)) {
    return(NA_character_)
  }
  if (p < 0.001) {
    return("<0.001")
  }
  sprintf("%.3f", p)
}

safe_fisher <- function(tab) {
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(NA_real_)
  }
  tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)
}

safe_wilcox <- function(df) {
  df2 <- df %>% filter(!is.na(days_between_collection), !is.na(snv_group))
  if (n_distinct(df2$snv_group) != 2) {
    return(NA_real_)
  }
  counts <- df2 %>% count(snv_group, name = "n")
  if (any(counts$n == 0)) {
    return(NA_real_)
  }
  tryCatch(wilcox.test(days_between_collection ~ snv_group, data = df2)$p.value, error = function(e) NA_real_)
}

metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
names(metadata) <- trimws(names(metadata))

lookup <- metadata %>%
  transmute(
    genome_id = as.character(genome_id),
    patient_id = trimws(as.character(patient_id)),
    ST = trimws(as.character(ST)),
    collection_date = parse_date_safely(collection_date),
    clade_group = trimws(as.character(clade_group))
  ) %>%
  mutate(
    patient_id = na_if(patient_id, ""),
    ST = na_if(ST, "")
  ) %>%
  distinct(genome_id, .keep_all = TRUE)

prepare_pairs <- function(path, clade_label) {
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
      ST_1 = ST,
      collection_date_1 = collection_date,
      clade_group_1 = clade_group
    ) %>%
    left_join(lookup, by = c("id2" = "genome_id")) %>%
    rename(
      patient_id_2 = patient_id,
      ST_2 = ST,
      collection_date_2 = collection_date,
      clade_group_2 = clade_group
    ) %>%
    mutate(
      same_patient = case_when(
        !is.na(patient_id_1) & !is.na(patient_id_2) & patient_id_1 == patient_id_2 ~ "Within patient",
        !is.na(patient_id_1) & !is.na(patient_id_2) & patient_id_1 != patient_id_2 ~ "Between patient",
        TRUE ~ NA_character_
      ),
      st_same = case_when(
        !is.na(ST_1) & !is.na(ST_2) & ST_1 == ST_2 ~ "ST same",
        !is.na(ST_1) & !is.na(ST_2) & ST_1 != ST_2 ~ "ST different",
        TRUE ~ NA_character_
      ),
      days_between_collection = abs(as.numeric(collection_date_2 - collection_date_1)),
      snv_group = case_when(
        !is.na(snv) & snv <= 10 ~ "SNV <= 10",
        !is.na(snv) & snv > 10 ~ "SNV > 10",
        TRUE ~ NA_character_
      ),
      clade = clade_label
    ) %>%
    filter(
      same_patient == "Within patient",
      !is.na(snv_group),
      !is.na(ST_1),
      !is.na(ST_2),
      !is.na(collection_date_1),
      !is.na(collection_date_2),
      !is.na(days_between_collection)
    )
}

run_st_same_test <- function(df, label) {
  df2 <- df %>%
    filter(!is.na(st_same))

  contingency <- table(df2$snv_group, df2$st_same)
  p_val <- safe_fisher(contingency)

  counts_df <- as.data.frame(contingency, stringsAsFactors = FALSE) %>%
    rename(
      snv_group = Var1,
      st_same = Var2,
      pair_n = Freq
    ) %>%
    mutate(
      analysis = label,
      fisher_p_value = format_p_value(p_val)
    ) %>%
    select(analysis, snv_group, st_same, pair_n, fisher_p_value)

  counts_df
}

run_time_test <- function(df, label) {
  df2 <- df

  p_val <- safe_wilcox(df2)

  summary_df <- df2 %>%
    group_by(snv_group) %>%
    summarise(
      pair_n = n(),
      median_days_between = median(days_between_collection, na.rm = TRUE),
      IQR_days_between = IQR(days_between_collection, na.rm = TRUE),
      min_days_between = min(days_between_collection, na.rm = TRUE),
      max_days_between = max(days_between_collection, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      analysis = label,
      wilcox_p_value = format_p_value(p_val)
    ) %>%
    select(analysis, snv_group, pair_n, median_days_between, IQR_days_between, min_days_between, max_days_between, wilcox_p_value)

  summary_df
}

clade1_pairs <- prepare_pairs(clade1_snv_path, "Clade 1")
clade2_pairs <- prepare_pairs(clade2_snv_path, "Clade 2")
all_pairs <- bind_rows(clade1_pairs, clade2_pairs)

st_test_results <- bind_rows(
  run_st_same_test(all_pairs, "Overall"),
  run_st_same_test(clade1_pairs, "Clade 1"),
  run_st_same_test(clade2_pairs, "Clade 2")
)

time_test_results <- bind_rows(
  run_time_test(all_pairs, "Overall"),
  run_time_test(clade1_pairs, "Clade 1"),
  run_time_test(clade2_pairs, "Clade 2")
)

pair_table <- all_pairs %>%
  select(
    clade,
    id1,
    id2,
    patient_id_1,
    ST_1,
    ST_2,
    st_same,
    collection_date_1,
    collection_date_2,
    days_between_collection,
    snv,
    snv_group
  ) %>%
  arrange(clade, snv_group, desc(snv), desc(days_between_collection))

write.csv(
  pair_table,
  file.path(desc_dir, "within_patient_pairs_st_and_time_for_snv_tests.csv"),
  row.names = FALSE
)

write.csv(
  st_test_results,
  file.path(desc_dir, "within_patient_st_same_vs_snv_group_fisher_tests.csv"),
  row.names = FALSE
)

write.csv(
  time_test_results,
  file.path(desc_dir, "within_patient_days_between_by_snv_group_wilcox_tests.csv"),
  row.names = FALSE
)

writeLines(
  c(
    "Definitions used:",
    "Pairs restricted to within-patient pairs from the filtered metadata_case_unformed_clade1_clade2 cohort.",
    "SNV <= 10 and SNV > 10 were compared.",
    "Both tests use the same complete-case subset requiring non-missing ST for both genomes and non-missing collection dates for both genomes.",
    "ST same vs ST different was tested with Fisher's exact test.",
    "Days between collection dates for SNV <= 10 vs SNV > 10 was tested with Wilcoxon rank-sum test."
  ),
  con = file.path(desc_dir, "within_patient_snv_st_time_tests_notes.txt")
)

cat("Done. Within-patient ST and collection-gap test outputs written to:\n", desc_dir, "\n", sep = "")
