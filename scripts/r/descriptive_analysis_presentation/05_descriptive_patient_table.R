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
metadata_path <- file.path(desc_dir, "metadata_case_unformed_clade1_clade2.csv")

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

pick_first_non_missing <- function(x) {
  x <- x[!is.na(x) & trimws(as.character(x)) != ""]
  if (length(x) == 0) {
    return(NA_character_)
  }
  as.character(x[[1]])
}

format_mean_sd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return("NA")
  }
  sprintf("%.1f (%.1f)", mean(x), stats::sd(x))
}

format_n_pct <- function(n, denom) {
  if (is.na(denom) || denom == 0) {
    return("0 (0.0%)")
  }
  sprintf("%d (%.1f%%)", n, 100 * n / denom)
}

format_p_value <- function(p) {
  if (is.na(p)) {
    return("")
  }
  if (p < 0.001) {
    return("<0.001")
  }
  sprintf("%.3f", p)
}

safe_t_test <- function(x, group) {
  dat <- tibble(x = suppressWarnings(as.numeric(x)), group = group) %>%
    filter(!is.na(x), !is.na(group))

  if (n_distinct(dat$group) != 2) {
    return(NA_real_)
  }

  group_sizes <- dat %>%
    count(group, name = "n")

  if (any(group_sizes$n < 2)) {
    return(NA_real_)
  }

  tryCatch(stats::t.test(x ~ group, data = dat)$p.value, error = function(e) NA_real_)
}

safe_cat_test <- function(x, group) {
  dat <- tibble(x = x, group = group) %>%
    filter(!is.na(x), !is.na(group))

  if (n_distinct(dat$group) != 2 || n_distinct(dat$x) < 2) {
    return(NA_real_)
  }

  tab <- table(dat$x, dat$group)
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(NA_real_)
  }

  tryCatch({
    chi_fit <- suppressWarnings(chisq.test(tab))
    if (any(chi_fit$expected < 5)) {
      suppressWarnings(fisher.test(tab, simulate.p.value = TRUE, B = 10000)$p.value)
    } else {
      chi_fit$p.value
    }
  }, error = function(e) NA_real_)
}

make_numeric_row <- function(label, var, total_df, clade1_df, clade2_df) {
  p_val <- safe_t_test(
    x = bind_rows(
      clade1_df %>% transmute(value = .data[[var]], group = "Clade 1"),
      clade2_df %>% transmute(value = .data[[var]], group = "Clade 2")
    )$value,
    group = bind_rows(
      clade1_df %>% transmute(value = .data[[var]], group = "Clade 1"),
      clade2_df %>% transmute(value = .data[[var]], group = "Clade 2")
    )$group
  )

  tibble(
    variable = label,
    total = format_mean_sd(total_df[[var]]),
    clade_1 = format_mean_sd(clade1_df[[var]]),
    clade_2 = format_mean_sd(clade2_df[[var]]),
    p_value = format_p_value(p_val)
  )
}

make_categorical_rows <- function(label, var, total_df, clade1_df, clade2_df, levels = NULL) {
  total_vals <- total_df[[var]]
  clade1_vals <- clade1_df[[var]]
  clade2_vals <- clade2_df[[var]]

  if (is.null(levels)) {
    levels <- c(
      sort(unique(c(
        as.character(total_vals[!is.na(total_vals)]),
        as.character(clade1_vals[!is.na(clade1_vals)]),
        as.character(clade2_vals[!is.na(clade2_vals)])
      )))
    )
  }

  p_val <- safe_cat_test(
    x = bind_rows(
      clade1_df %>% transmute(value = .data[[var]], group = "Clade 1"),
      clade2_df %>% transmute(value = .data[[var]], group = "Clade 2")
    )$value,
    group = bind_rows(
      clade1_df %>% transmute(value = .data[[var]], group = "Clade 1"),
      clade2_df %>% transmute(value = .data[[var]], group = "Clade 2")
    )$group
  )

  bind_rows(
    tibble(
      variable = label,
      total = "",
      clade_1 = "",
      clade_2 = "",
      p_value = format_p_value(p_val)
    ),
    bind_rows(lapply(unname(levels), function(lvl) {
      tibble(
        variable = paste0("  ", lvl),
        total = format_n_pct(sum(total_vals == lvl, na.rm = TRUE), nrow(total_df)),
        clade_1 = format_n_pct(sum(clade1_vals == lvl, na.rm = TRUE), nrow(clade1_df)),
        clade_2 = format_n_pct(sum(clade2_vals == lvl, na.rm = TRUE), nrow(clade2_df)),
        p_value = ""
      )
    }))
  )
}

metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
names(metadata) <- trimws(names(metadata))

if (!"patient_id" %in% names(metadata)) {
  stop("Expected patient_id column in metadata_case_unformed_clade1_clade2.csv")
}

metadata_clean <- metadata %>%
  mutate(
    patient_id = trimws(as.character(patient_id)),
    clade_group = trimws(as.character(clade_group)),
    collection_date = suppressWarnings(as.Date(collection_date)),
    age = suppressWarnings(as.numeric(age)),
    BMI = suppressWarnings(as.numeric(BMI)),
    LOS = suppressWarnings(as.numeric(LOS)),
    GenderCode = case_when(
      trimws(as.character(GenderCode)) %in% c("F", "Female") ~ "Female",
      trimws(as.character(GenderCode)) %in% c("M", "Male") ~ "Male",
      TRUE ~ "Missing"
    ),
    RaceName = case_when(
      is.na(RaceName) | trimws(as.character(RaceName)) == "" ~ "Unknown/Missing",
      trimws(as.character(RaceName)) %in% c("Unknown", "Missing") ~ "Unknown/Missing",
      TRUE ~ trimws(as.character(RaceName))
    ),
    inpatient = case_when(
      as.character(inpatient) %in% c("1", "TRUE", "True", "true", "Yes", "Y") ~ "Yes",
      as.character(inpatient) %in% c("0", "FALSE", "False", "false", "No", "N") ~ "No",
      TRUE ~ "Missing"
    ),
    toxin_status_updated = case_when(
      as.character(has_any_toxin) %in% c("TRUE", "True", "true", "1") ~ "Toxigenic",
      as.character(has_any_toxin) %in% c("FALSE", "False", "false", "0") ~ "Non-Toxigenic",
      TRUE ~ "Missing"
    ),
    hospital_onset_updated = case_when(
      trimws(as.character(hospital_onset_updated)) %in% c("HO", "CO", "rCDI", "unknown/NA") ~ trimws(as.character(hospital_onset_updated)),
      TRUE ~ "unknown/NA"
    )
  ) %>%
  filter(!is.na(patient_id), patient_id != "", clade_group %in% c("Clade 1", "Clade 2")) %>%
  arrange(patient_id, clade_group, collection_date)

patient_total <- metadata_clean %>%
  group_by(patient_id) %>%
  summarise(
    age = pick_first_non_missing(age),
    BMI = pick_first_non_missing(BMI),
    GenderCode = pick_first_non_missing(GenderCode),
    RaceName = pick_first_non_missing(RaceName),
    LOS = pick_first_non_missing(LOS),
    inpatient = pick_first_non_missing(inpatient),
    toxin_status_updated = pick_first_non_missing(toxin_status_updated),
    hospital_onset_updated = pick_first_non_missing(hospital_onset_updated),
    .groups = "drop"
  ) %>%
  mutate(
    age = suppressWarnings(as.numeric(age)),
    BMI = suppressWarnings(as.numeric(BMI)),
    LOS = suppressWarnings(as.numeric(LOS))
  )

patient_by_clade <- metadata_clean %>%
  group_by(clade_group, patient_id) %>%
  summarise(
    age = pick_first_non_missing(age),
    BMI = pick_first_non_missing(BMI),
    GenderCode = pick_first_non_missing(GenderCode),
    RaceName = pick_first_non_missing(RaceName),
    LOS = pick_first_non_missing(LOS),
    inpatient = pick_first_non_missing(inpatient),
    toxin_status_updated = pick_first_non_missing(toxin_status_updated),
    hospital_onset_updated = pick_first_non_missing(hospital_onset_updated),
    .groups = "drop"
  ) %>%
  mutate(
    age = suppressWarnings(as.numeric(age)),
    BMI = suppressWarnings(as.numeric(BMI)),
    LOS = suppressWarnings(as.numeric(LOS))
  )

patient_clade1 <- patient_by_clade %>% filter(clade_group == "Clade 1")
patient_clade2 <- patient_by_clade %>% filter(clade_group == "Clade 2")

total_n <- nrow(patient_total)
clade1_n <- nrow(patient_clade1)
clade2_n <- nrow(patient_clade2)

race_levels <- c(
  "Caucasian",
  "African American",
  "Asian",
  "American Indian or Alaska Native",
  "Other",
  "Unknown/Missing"
)
race_levels <- race_levels[race_levels %in% unique(c(patient_total$RaceName, patient_clade1$RaceName, patient_clade2$RaceName))]

table1_df <- bind_rows(
  make_numeric_row("Age (years)", "age", patient_total, patient_clade1, patient_clade2),
  make_numeric_row("BMI", "BMI", patient_total, patient_clade1, patient_clade2),
  make_categorical_rows("Sex", "GenderCode", patient_total, patient_clade1, patient_clade2, levels = c("Female", "Male", "Missing")),
  make_categorical_rows("Race", "RaceName", patient_total, patient_clade1, patient_clade2, levels = race_levels),
  make_numeric_row("Length of stay (days)", "LOS", patient_total, patient_clade1, patient_clade2),
  make_categorical_rows("Inpatient", "inpatient", patient_total, patient_clade1, patient_clade2, levels = c("Yes", "No", "Missing")),
  make_categorical_rows("Toxin status", "toxin_status_updated", patient_total, patient_clade1, patient_clade2, levels = c("Toxigenic", "Non-Toxigenic", "Missing")),
  make_categorical_rows("Hospital onset updated", "hospital_onset_updated", patient_total, patient_clade1, patient_clade2, levels = c("HO", "CO", "rCDI", "unknown/NA"))
) %>%
  rename(
    `Variable` = variable,
    !!paste0("Total (n=", total_n, ")") := total,
    !!paste0("Clade 1 (n=", clade1_n, ")") := clade_1,
    !!paste0("Clade 2 (n=", clade2_n, ")") := clade_2,
    `P value` = p_value
  )

write.csv(
  table1_df,
  file.path(desc_dir, "05_descriptive_patient_table_case_unformed_clade1_clade2.csv"),
  row.names = FALSE,
  na = ""
)

writeLines(
  c(
    paste0("Total unique patients: ", total_n),
    paste0("Clade 1 patients: ", clade1_n),
    paste0("Clade 2 patients: ", clade2_n),
    "Patients can contribute to both clade columns if they have isolates in both clades.",
    "Continuous variables are shown as mean (SD).",
    "Categorical variables are shown as n (%).",
    "P values compare Clade 1 vs Clade 2 patients.",
    "Hospital onset updated definitions:",
    "HO: positive stool specimen collected more than 3 days after admission.",
    "CO: positive stool specimen collected on days 1 to 3 after admission.",
    "rCDI: recurrent CDI episode occurring within 56 days of a previous CDI episode.",
    "Positive tests collected within 14 days of a prior positive test are treated as part of the same episode, not a new episode.",
    "unknown/NA: insufficient timing/history information to classify as HO, CO, or rCDI."
  ),
  con = file.path(desc_dir, "05_descriptive_patient_table_notes.txt")
)

cat("Done. Descriptive patient table written to:\n", desc_dir, "\n", sep = "")
