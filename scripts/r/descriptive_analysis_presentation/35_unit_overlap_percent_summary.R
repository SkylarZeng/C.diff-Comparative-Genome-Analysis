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

locations_path <- file.path(desc_dir, "locations_filtered_case_unformed_clade1_clade2.tsv")

required_pkgs <- c("tidyverse")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(tidyverse)
})

special_units <- c("ESA", "OR", "DVU", "MPUA", "UADM", "USSS", "ORC", "VHOR", "IRU", "CATH", "MHOR", "INU")
sequential_gap_days <- 14

is_present <- function(x) {
  !is.na(x) & trimws(as.character(x)) != "" & trimws(as.character(x)) != "NA"
}

parse_date_safe <- function(x) {
  suppressWarnings(as.Date(substr(trimws(as.character(x)), 1, 10)))
}

parse_datetime_safe <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | x == "NA"] <- NA_character_
  parsed <- suppressWarnings(as.POSIXct(x, tz = "UTC", format = "%Y-%m-%dT%H:%M:%OSZ"))
  fallback_idx <- is.na(parsed) & !is.na(x)
  if (any(fallback_idx)) {
    parsed[fallback_idx] <- suppressWarnings(as.POSIXct(x[fallback_idx], tz = "UTC"))
  }
  parsed
}

collapse_patient_intervals <- function(df) {
  if (nrow(df) == 0) {
    return(tibble())
  }

  df <- df[order(df$StartDate_dt, df$EndDate_dt), , drop = FALSE]
  out_start <- as.POSIXct(character(), tz = "UTC")
  out_end <- as.POSIXct(character(), tz = "UTC")
  current_start <- df$StartDate_dt[1]
  current_end <- df$EndDate_dt[1]

  if (nrow(df) > 1) {
    for (i in 2:nrow(df)) {
      next_start <- df$StartDate_dt[i]
      next_end <- df$EndDate_dt[i]
      if (next_start <= current_end) {
        current_end <- max(current_end, next_end)
      } else {
        out_start <- c(out_start, current_start)
        out_end <- c(out_end, current_end)
        current_start <- next_start
        current_end <- next_end
      }
    }
  }

  out_start <- c(out_start, current_start)
  out_end <- c(out_end, current_end)

  tibble(StartDate_dt = out_start, EndDate_dt = out_end)
}

find_direct_pairs_by_unit <- function(loc_df) {
  dat <- loc_df %>%
    filter(is_present(Unit)) %>%
    select(PatientID, Unit, StartDate_dt, EndDate_dt) %>%
    distinct() %>%
    group_by(PatientID, Unit) %>%
    group_modify(~ collapse_patient_intervals(.x)) %>%
    ungroup()

  if (nrow(dat) == 0) {
    return(tibble(Unit = character(), patient_a = character(), patient_b = character()))
  }

  by_unit <- split(dat, dat$Unit)
  out <- list()
  idx <- 1L

  for (unit_name in names(by_unit)) {
    group_df <- by_unit[[unit_name]]
    group_df <- group_df[order(group_df$StartDate_dt, group_df$EndDate_dt), , drop = FALSE]
    n <- nrow(group_df)
    if (n < 2) next

    pair_keys_env <- new.env(parent = emptyenv(), hash = TRUE)
    for (i in seq_len(n - 1)) {
      j <- i + 1
      while (j <= n && group_df$StartDate_dt[j] < group_df$EndDate_dt[i]) {
        if (
          group_df$PatientID[i] != group_df$PatientID[j] &&
          group_df$StartDate_dt[i] < group_df$EndDate_dt[j] &&
          group_df$StartDate_dt[j] < group_df$EndDate_dt[i]
        ) {
          pair_a <- min(group_df$PatientID[i], group_df$PatientID[j])
          pair_b <- max(group_df$PatientID[i], group_df$PatientID[j])
          pair_keys_env[[paste(pair_a, pair_b, sep = "||")]] <- TRUE
        }
        j <- j + 1
      }
    }

    pair_keys <- ls(pair_keys_env, all.names = TRUE)
    if (length(pair_keys) == 0) next

    tmp <- tibble(pair_key = pair_keys, Unit = unit_name) %>%
      separate(pair_key, into = c("patient_a", "patient_b"), sep = "\\|\\|", remove = TRUE)
    out[[idx]] <- tmp
    idx <- idx + 1L
  }

  bind_rows(out)
}

find_sequential_pairs_by_unit <- function(loc_df, max_gap_days = 14) {
  dat <- loc_df %>%
    filter(is_present(Unit)) %>%
    select(PatientID, Unit, StartDate_dt, EndDate_dt) %>%
    distinct() %>%
    group_by(PatientID, Unit) %>%
    group_modify(~ collapse_patient_intervals(.x)) %>%
    ungroup()

  if (nrow(dat) == 0) {
    return(tibble(Unit = character(), patient_a = character(), patient_b = character()))
  }

  by_unit <- split(dat, dat$Unit)
  out <- list()
  idx <- 1L

  for (unit_name in names(by_unit)) {
    group_df <- by_unit[[unit_name]]
    group_df <- group_df[order(group_df$StartDate_dt, group_df$EndDate_dt), , drop = FALSE]
    n <- nrow(group_df)
    if (n < 2) next

    pair_keys_env <- new.env(parent = emptyenv(), hash = TRUE)
    for (i in seq_len(n - 1)) {
      j <- i + 1
      while (j <= n && group_df$StartDate_dt[j] <= group_df$EndDate_dt[i] + as.difftime(max_gap_days, units = "days")) {
        gap_days <- as.numeric(difftime(group_df$StartDate_dt[j], group_df$EndDate_dt[i], units = "days"))
        if (
          group_df$PatientID[i] != group_df$PatientID[j] &&
          is.finite(gap_days) &&
          gap_days > 0 &&
          gap_days <= max_gap_days
        ) {
          pair_a <- min(group_df$PatientID[i], group_df$PatientID[j])
          pair_b <- max(group_df$PatientID[i], group_df$PatientID[j])
          pair_keys_env[[paste(pair_a, pair_b, sep = "||")]] <- TRUE
        }
        j <- j + 1
      }
    }

    pair_keys <- ls(pair_keys_env, all.names = TRUE)
    if (length(pair_keys) == 0) next

    tmp <- tibble(pair_key = pair_keys, Unit = unit_name) %>%
      separate(pair_key, into = c("patient_a", "patient_b"), sep = "\\|\\|", remove = TRUE)
    out[[idx]] <- tmp
    idx <- idx + 1L
  }

  bind_rows(out)
}

locations <- read.delim(locations_path, check.names = FALSE, stringsAsFactors = FALSE) %>%
  mutate(
    PatientID = trimws(as.character(PatientID)),
    clade_group = trimws(as.character(clade_group)),
    Unit = trimws(as.character(Unit)),
    Date_date = parse_date_safe(Date),
    StartDate_dt = coalesce(parse_datetime_safe(StartDate), as.POSIXct(Date_date, tz = "UTC")),
    EndDate_dt = coalesce(parse_datetime_safe(EndDate), parse_datetime_safe(StartDate), as.POSIXct(Date_date, tz = "UTC"))
  ) %>%
  filter(
    !is.na(PatientID),
    PatientID != "",
    !is.na(StartDate_dt),
    !is.na(EndDate_dt),
    is_present(Unit),
    !(Unit %in% special_units)
  )

make_summary_for_clade <- function(clade_label) {
  loc_clade <- locations %>% filter(clade_group == clade_label)

  unit_patient_counts <- loc_clade %>%
    group_by(Unit) %>%
    summarise(patient_n = n_distinct(PatientID), .groups = "drop")

  direct_pairs <- find_direct_pairs_by_unit(loc_clade)
  sequential_pairs <- find_sequential_pairs_by_unit(loc_clade, max_gap_days = sequential_gap_days)

  direct_summary <- direct_pairs %>%
    group_by(Unit) %>%
    summarise(direct_pair_n = n(), .groups = "drop")

  sequential_summary <- sequential_pairs %>%
    group_by(Unit) %>%
    summarise(sequential_14d_pair_n = n(), .groups = "drop")

  total_direct_pairs <- sum(direct_summary$direct_pair_n, na.rm = TRUE)
  total_sequential_pairs <- sum(sequential_summary$sequential_14d_pair_n, na.rm = TRUE)

  direct_pct <- if (total_direct_pairs > 0) {
    100 * direct_summary$direct_pair_n / total_direct_pairs
  } else {
    numeric(nrow(direct_summary))
  }

  sequential_pct <- if (total_sequential_pairs > 0) {
    100 * sequential_summary$sequential_14d_pair_n / total_sequential_pairs
  } else {
    numeric(nrow(sequential_summary))
  }

  direct_summary <- direct_summary %>%
    mutate(direct_pair_pct_within_clade = direct_pct)

  sequential_summary <- sequential_summary %>%
    mutate(sequential_14d_pair_pct_within_clade = sequential_pct)

  full_join(direct_summary, sequential_summary, by = "Unit") %>%
    full_join(unit_patient_counts, by = "Unit") %>%
    mutate(
      clade_group = clade_label,
      direct_pair_n = replace_na(direct_pair_n, 0L),
      sequential_14d_pair_n = replace_na(sequential_14d_pair_n, 0L),
      patient_n = replace_na(patient_n, 0L),
      direct_pair_pct_within_clade = replace_na(direct_pair_pct_within_clade, 0),
      sequential_14d_pair_pct_within_clade = replace_na(sequential_14d_pair_pct_within_clade, 0)
    ) %>%
    arrange(desc(direct_pair_n), desc(sequential_14d_pair_n), Unit) %>%
    select(
      clade_group,
      Unit,
      patient_n,
      direct_pair_n,
      direct_pair_pct_within_clade,
      sequential_14d_pair_n,
      sequential_14d_pair_pct_within_clade
    )
}

summary_table <- bind_rows(
  make_summary_for_clade("Clade 1"),
  make_summary_for_clade("Clade 2")
)

write.csv(
  summary_table,
  file.path(desc_dir, "unit_direct_and_sequential_overlap_percent_by_clade.csv"),
  row.names = FALSE
)

writeLines(
  c(
    "This table summarizes unit-level overlap after excluding special units.",
    paste0("Sequential overlap uses the current ", sequential_gap_days, "-day threshold."),
    "direct_pair_n = number of distinct patient-pair overlaps directly overlapping in time within that unit.",
    paste0("sequential_", sequential_gap_days, "d_pair_n = number of distinct patient-pair exposures in that unit with no direct overlap and a 1-", sequential_gap_days, "-day gap."),
    "Percent columns are calculated within each clade separately, using the sum of unit-level direct or sequential patient-pair counts as the denominator.",
    "A patient pair can contribute to more than one unit if they overlap in multiple units."
  ),
  con = file.path(desc_dir, "unit_direct_and_sequential_overlap_percent_by_clade_notes.txt")
)

cat("Done. Unit overlap summary written to:\n", desc_dir, "\n", sep = "")
