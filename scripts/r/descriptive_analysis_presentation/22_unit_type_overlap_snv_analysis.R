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
viz_dir <- file.path(output_dir, "visualization")

clade1_snv_path <- file.path(cluster_dir, "SNV", "clade1_post_noncore_snv_long.csv")
clade2_snv_path <- file.path(cluster_dir, "SNV", "clade2_post_noncore_snv_long.csv")
metadata_path <- file.path(desc_dir, "metadata_case_unformed_clade1_clade2.csv")
locations_path <- file.path(desc_dir, "locations_filtered_case_unformed_clade1_clade2.tsv")

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

safe_wilcox <- function(x_yes, x_no) {
  if (length(x_yes) == 0 || length(x_no) == 0) {
    return(NA_real_)
  }
  suppressWarnings(wilcox.test(x_yes, x_no)$p.value)
}

metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
names(metadata) <- trimws(names(metadata))

patient_lookup <- metadata %>%
  transmute(
    genome_id = as.character(genome_id),
    patient_id = trimws(as.character(patient_id))
  ) %>%
  mutate(patient_id = na_if(patient_id, "")) %>%
  distinct(genome_id, .keep_all = TRUE)

special_units <- c("ESA", "OR", "DVU", "MPUA", "UADM", "USSS", "ORC", "VHOR", "IRU", "CATH", "MHOR", "INU")

locations <- read.delim(locations_path, check.names = FALSE, stringsAsFactors = FALSE) %>%
  mutate(
    PatientID = trimws(as.character(PatientID)),
    clade_group = trimws(as.character(clade_group)),
    FacilityCode = trimws(as.character(FacilityCode)),
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
    is_present(FacilityCode)
  )

prepare_snv_pairs <- function(path, clade_label) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE) %>%
    rename(snv = snvs) %>%
    mutate(
      id1 = as.character(id1),
      id2 = as.character(id2),
      snv = suppressWarnings(as.numeric(snv)),
      clade_group = clade_label
    ) %>%
    left_join(patient_lookup, by = c("id1" = "genome_id")) %>%
    rename(patient_id_1 = patient_id) %>%
    left_join(patient_lookup, by = c("id2" = "genome_id")) %>%
    rename(patient_id_2 = patient_id) %>%
    mutate(
      same_patient = case_when(
        !is.na(patient_id_1) & !is.na(patient_id_2) & patient_id_1 == patient_id_2 ~ "Within patient",
        !is.na(patient_id_1) & !is.na(patient_id_2) & patient_id_1 != patient_id_2 ~ "Between patient",
        TRUE ~ NA_character_
      ),
      patient_a = if_else(!is.na(patient_id_1) & !is.na(patient_id_2), pmin(patient_id_1, patient_id_2), NA_character_),
      patient_b = if_else(!is.na(patient_id_1) & !is.na(patient_id_2), pmax(patient_id_1, patient_id_2), NA_character_)
    )
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

find_unit_overlap_pairs <- function(loc_df) {
  dat <- loc_df %>%
    select(PatientID, StartDate_dt, EndDate_dt, FacilityCode, Unit) %>%
    unite("location_value", FacilityCode, Unit, sep = "||", remove = TRUE) %>%
    distinct() %>%
    group_by(PatientID, location_value) %>%
    group_modify(~ collapse_patient_intervals(.x)) %>%
    ungroup()

  if (nrow(dat) == 0) {
    return(tibble(patient_a = character(), patient_b = character()))
  }

  by_location <- split(dat, dat$location_value)
  pair_keys_env <- new.env(parent = emptyenv(), hash = TRUE)

  for (group_df in by_location) {
    group_df <- group_df[order(group_df$StartDate_dt, group_df$EndDate_dt), , drop = FALSE]
    n <- nrow(group_df)
    if (n < 2) {
      next
    }

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
  }

  pair_keys <- ls(pair_keys_env, all.names = TRUE)
  if (length(pair_keys) == 0) {
    return(tibble(patient_a = character(), patient_b = character()))
  }

  tibble(pair_key = pair_keys) %>%
    separate(pair_key, into = c("patient_a", "patient_b"), sep = "\\|\\|", remove = TRUE)
}

annotate_unit_type_overlap <- function(snv_pairs, clade_label, unit_type = c("normal", "special")) {
  unit_type <- match.arg(unit_type)

  loc_clade <- locations %>%
    filter(clade_group == clade_label) %>%
    filter(if (unit_type == "special") Unit %in% special_units else !(Unit %in% special_units))

  pair_flags <- snv_pairs %>%
    filter(same_patient == "Between patient", !is.na(patient_a), !is.na(patient_b)) %>%
    distinct(patient_a, patient_b)

  overlap_pairs <- find_unit_overlap_pairs(loc_clade) %>%
    mutate(unit_overlap = TRUE)

  pair_flags %>%
    left_join(overlap_pairs, by = c("patient_a", "patient_b")) %>%
    mutate(unit_overlap = replace_na(unit_overlap, FALSE))
}

make_distribution_plot <- function(dat, unit_type_label) {
  ggplot(dat, aes(x = snv, fill = overlap_status)) +
    geom_histogram(
      bins = 40,
      position = "identity",
      alpha = 0.5
    ) +
    scale_fill_manual(values = c("Overlap yes" = "#e15759", "Overlap no" = "#4e79a7"), name = NULL) +
    scale_x_continuous(limits = c(0, 75)) +
    facet_wrap(~clade_group, scales = "free_y", ncol = 2) +
    labs(
      title = paste0("Between-patient SNV by ", unit_type_label, " unit overlap (0-75 SNVs)"),
      subtitle = "Direct overlap uses full timestamps",
      x = "SNV distance",
      y = "Pairs"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

make_summary <- function(dat, unit_type_label) {
  bind_rows(
    dat %>% filter(clade_group == "Clade 1"),
    dat %>% filter(clade_group == "Clade 2")
  ) %>%
    group_by(clade_group, overlap_status) %>%
    summarise(
      pair_n = n(),
      median_snv = median(snv),
      iqr_snv = IQR(snv),
      .groups = "drop"
    ) %>%
    mutate(unit_type = unit_type_label)
}

run_test_summary <- function(dat, unit_type_label) {
  bind_rows(
    tibble(clade_group = "Clade 1", dat = list(dat %>% filter(clade_group == "Clade 1"))),
    tibble(clade_group = "Clade 2", dat = list(dat %>% filter(clade_group == "Clade 2")))
  ) %>%
    rowwise() %>%
    mutate(
      yes_n = sum(dat$overlap_status == "Overlap yes"),
      no_n = sum(dat$overlap_status == "Overlap no"),
      yes_median_snv = ifelse(yes_n > 0, median(dat$snv[dat$overlap_status == "Overlap yes"]), NA_real_),
      no_median_snv = ifelse(no_n > 0, median(dat$snv[dat$overlap_status == "Overlap no"]), NA_real_),
      wilcox_p_value = safe_wilcox(
        dat$snv[dat$overlap_status == "Overlap yes"],
        dat$snv[dat$overlap_status == "Overlap no"]
      )
    ) %>%
    ungroup() %>%
    transmute(
      unit_type = unit_type_label,
      clade_group,
      yes_n,
      no_n,
      yes_median_snv,
      no_median_snv,
      wilcox_p_value,
      p_value_label = case_when(
        is.na(wilcox_p_value) ~ NA_character_,
        wilcox_p_value < 0.001 ~ "<0.001",
        TRUE ~ formatC(wilcox_p_value, format = "f", digits = 3)
      )
    )
}

clade1_pairs <- prepare_snv_pairs(clade1_snv_path, "Clade 1")
clade2_pairs <- prepare_snv_pairs(clade2_snv_path, "Clade 2")

clade1_normal <- annotate_unit_type_overlap(clade1_pairs, "Clade 1", "normal")
clade2_normal <- annotate_unit_type_overlap(clade2_pairs, "Clade 2", "normal")
clade1_special <- annotate_unit_type_overlap(clade1_pairs, "Clade 1", "special")
clade2_special <- annotate_unit_type_overlap(clade2_pairs, "Clade 2", "special")

normal_dat <- bind_rows(
  clade1_pairs %>% left_join(clade1_normal, by = c("patient_a", "patient_b")),
  clade2_pairs %>% left_join(clade2_normal, by = c("patient_a", "patient_b"))
) %>%
  filter(same_patient == "Between patient", !is.na(snv), snv >= 0, snv <= 75) %>%
  mutate(overlap_status = if_else(replace_na(unit_overlap, FALSE), "Overlap yes", "Overlap no"))

special_dat <- bind_rows(
  clade1_pairs %>% left_join(clade1_special, by = c("patient_a", "patient_b")),
  clade2_pairs %>% left_join(clade2_special, by = c("patient_a", "patient_b"))
) %>%
  filter(same_patient == "Between patient", !is.na(snv), snv >= 0, snv <= 75) %>%
  mutate(overlap_status = if_else(replace_na(unit_overlap, FALSE), "Overlap yes", "Overlap no"))

normal_summary <- make_summary(normal_dat, "Normal")
special_summary <- make_summary(special_dat, "Special")
test_summary <- bind_rows(
  run_test_summary(normal_dat, "Normal"),
  run_test_summary(special_dat, "Special")
)

write.csv(normal_summary, file.path(desc_dir, "normal_unit_overlap_snv_distribution_summary.csv"), row.names = FALSE)
write.csv(special_summary, file.path(desc_dir, "special_unit_overlap_snv_distribution_summary.csv"), row.names = FALSE)
write.csv(test_summary, file.path(desc_dir, "unit_type_overlap_snv_wilcox_tests.csv"), row.names = FALSE)

writeLines(
  c(
    "Unit-type overlap analysis uses between-patient SNV pairs only.",
    paste0("Special units are defined as: ", paste(special_units, collapse = ", "), "."),
    "Direct overlap uses full timestamps, so same-day but non-overlapping stays are not counted as direct overlap.",
    "Plots are limited to SNV 0-75 for readability."
  ),
  con = file.path(desc_dir, "unit_type_overlap_snv_analysis_notes.txt")
)

normal_plot <- make_distribution_plot(normal_dat, "normal")
special_plot <- make_distribution_plot(special_dat, "special")

ggsave(
  filename = file.path(viz_dir, "normal_unit_overlap_snv_distribution.png"),
  plot = normal_plot,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = normal_plot,
  filename = file.path(viz_dir, "normal_unit_overlap_snv_distribution.pdf"),
  width = 12,
  height = 8
)

ggsave(
  filename = file.path(viz_dir, "special_unit_overlap_snv_distribution.png"),
  plot = special_plot,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = special_plot,
  filename = file.path(viz_dir, "special_unit_overlap_snv_distribution.pdf"),
  width = 12,
  height = 8
)

cat("Done. Unit-type overlap SNV analysis outputs written to:\n", desc_dir, "\n", sep = "")
