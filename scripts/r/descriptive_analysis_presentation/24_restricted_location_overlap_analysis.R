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

make_star_label <- function(p_label) {
  case_when(
    is.na(p_label) ~ "NA",
    p_label == "<0.001" ~ "***",
    suppressWarnings(as.numeric(p_label)) < 0.01 ~ "**",
    suppressWarnings(as.numeric(p_label)) < 0.05 ~ "*",
    TRUE ~ "ns"
  )
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
    floor = trimws(as.character(floor)),
    Unit = trimws(as.character(Unit)),
    Room = trimws(as.character(Room)),
    StartDate_dt = coalesce(parse_datetime_safe(StartDate), as.POSIXct(parse_date_safe(Date), tz = "UTC")),
    EndDate_dt = coalesce(parse_datetime_safe(EndDate), parse_datetime_safe(StartDate), as.POSIXct(parse_date_safe(Date), tz = "UTC"))
  ) %>%
  filter(
    !is.na(PatientID),
    PatientID != "",
    !is.na(StartDate_dt),
    !is.na(EndDate_dt),
    !(Unit %in% special_units)
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
    ) %>%
    filter(same_patient == "Between patient", !is.na(snv), !is.na(patient_a), !is.na(patient_b))
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

find_overlap_pairs_for_level <- function(loc_df, key_fields) {
  dat <- loc_df %>%
    filter(if_all(all_of(key_fields), is_present)) %>%
    select(PatientID, StartDate_dt, EndDate_dt, all_of(key_fields)) %>%
    unite("location_value", all_of(key_fields), sep = "||", remove = TRUE) %>%
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

compute_restricted_flags <- function(snv_pairs, clade_label) {
  loc_clade <- locations %>%
    filter(clade_group == clade_label)

  pair_flags <- snv_pairs %>%
    distinct(patient_a, patient_b)

  hospital_time_overlap <- find_overlap_pairs_for_level(
    loc_clade %>% transmute(PatientID, StartDate_dt, EndDate_dt, hospital_key = "hospital"),
    c("hospital_key")
  ) %>%
    mutate(hospital_time_overlap = TRUE)

  facility_overlap <- find_overlap_pairs_for_level(loc_clade, c("FacilityCode")) %>%
    mutate(facility_overlap = TRUE)

  floor_overlap <- find_overlap_pairs_for_level(loc_clade, c("floor")) %>%
    mutate(floor_overlap = TRUE)

  unit_overlap <- find_overlap_pairs_for_level(loc_clade, c("FacilityCode", "Unit")) %>%
    mutate(unit_overlap = TRUE)

  room_overlap <- find_overlap_pairs_for_level(loc_clade, c("FacilityCode", "Room")) %>%
    mutate(room_overlap = TRUE)

  pair_flags %>%
    left_join(hospital_time_overlap, by = c("patient_a", "patient_b")) %>%
    left_join(facility_overlap, by = c("patient_a", "patient_b")) %>%
    left_join(floor_overlap, by = c("patient_a", "patient_b")) %>%
    left_join(unit_overlap, by = c("patient_a", "patient_b")) %>%
    left_join(room_overlap, by = c("patient_a", "patient_b")) %>%
    mutate(across(where(is.logical), ~ replace_na(.x, FALSE)))
}

build_restricted_plot_data <- function(dat) {
  bind_rows(
    dat %>%
      filter(hospital_time_overlap) %>%
      transmute(clade_group, snv, exposure_type = "Facility", exposure_status = if_else(facility_overlap, "Overlap yes", "Overlap no")),
    dat %>%
      filter(hospital_time_overlap) %>%
      transmute(clade_group, snv, exposure_type = "Floor", exposure_status = if_else(floor_overlap, "Overlap yes", "Overlap no")),
    dat %>%
      filter(hospital_time_overlap) %>%
      transmute(clade_group, snv, exposure_type = "Unit", exposure_status = if_else(unit_overlap, "Overlap yes", "Overlap no")),
    dat %>%
      filter(hospital_time_overlap) %>%
      transmute(clade_group, snv, exposure_type = "Room", exposure_status = if_else(room_overlap, "Overlap yes", "Overlap no"))
  ) %>%
    mutate(exposure_type = factor(exposure_type, levels = c("Facility", "Floor", "Unit", "Room")))
}

run_threshold_tests <- function(plot_dat) {
  thresholds <- c(2, 5, 10)

  crossing(
    clade_group = c("Clade 1", "Clade 2"),
    exposure_type = c("Facility", "Floor", "Unit", "Room"),
    threshold = thresholds
  ) %>%
    rowwise() %>%
    mutate(
      dat = list(plot_dat %>% filter(clade_group == .data$clade_group, exposure_type == .data$exposure_type)),
      yes_low_n = sum(dat$snv[dat$exposure_status == "Overlap yes"] <= threshold, na.rm = TRUE),
      yes_high_n = sum(dat$snv[dat$exposure_status == "Overlap yes"] > threshold, na.rm = TRUE),
      no_low_n = sum(dat$snv[dat$exposure_status == "Overlap no"] <= threshold, na.rm = TRUE),
      no_high_n = sum(dat$snv[dat$exposure_status == "Overlap no"] > threshold, na.rm = TRUE),
      yes_low_pct = ifelse((yes_low_n + yes_high_n) > 0, 100 * yes_low_n / (yes_low_n + yes_high_n), NA_real_),
      no_low_pct = ifelse((no_low_n + no_high_n) > 0, 100 * no_low_n / (no_low_n + no_high_n), NA_real_),
      fisher_p_value = safe_fisher(yes_low_n, yes_high_n, no_low_n, no_high_n),
      p_value_label = format_p(fisher_p_value)
    ) %>%
    ungroup() %>%
    transmute(
      clade_group,
      exposure_type,
      threshold_rule = paste0("SNV <= ", threshold),
      yes_low_n,
      yes_high_n,
      no_low_n,
      no_high_n,
      yes_low_pct,
      no_low_pct,
      fisher_p_value,
      p_value_label
    )
}

make_distribution_plot <- function(plot_dat) {
  p <- plot_dat %>%
    filter(snv >= 0, snv <= 75) %>%
    ggplot(aes(x = snv, fill = exposure_status)) +
    geom_histogram(
      bins = 40,
      position = "identity",
      alpha = 0.5
    ) +
    scale_fill_manual(values = c("Overlap yes" = "#E0607E", "Overlap no" = "#54BFB7"), name = NULL) +
    scale_x_continuous(limits = c(0, 75)) +
    facet_grid(clade_group ~ exposure_type, scales = "free_y") +
    labs(
      title = "Restricted between-patient SNV by location overlap",
      subtitle = "Facility compares only pairs with concurrent hospitalization; Floor/Unit/Room compare only pairs with concurrent same-facility time",
      x = "SNV distance",
      y = "Pairs"
    ) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(viz_dir, "restricted_location_overlap_snv_0_75.png"),
    plot = p,
    width = 14,
    height = 8,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  save_plot_pdf(
    plot_obj = p,
    filename = file.path(viz_dir, "restricted_location_overlap_snv_0_75.pdf"),
    width = 14,
    height = 8
  )
}

make_threshold_plot <- function(threshold_dat) {
  plot_dat <- bind_rows(
    threshold_dat %>% transmute(clade_group, exposure_type, threshold_rule, group = "Yes", pct = yes_low_pct, p_value_label),
    threshold_dat %>% transmute(clade_group, exposure_type, threshold_rule, group = "No", pct = no_low_pct, p_value_label)
  ) %>%
    mutate(
      group = factor(group, levels = c("Yes", "No")),
      threshold_rule = factor(threshold_rule, levels = c("SNV <= 2", "SNV <= 5", "SNV <= 10")),
      exposure_type = factor(exposure_type, levels = c("Facility", "Floor", "Unit", "Room")),
      star_label = make_star_label(p_value_label)
    )

  ann_dat <- plot_dat %>%
    group_by(clade_group, exposure_type, threshold_rule) %>%
    summarise(
      y_pos = max(pct, na.rm = TRUE) * 1.15 + 0.2,
      star_label = first(star_label),
      .groups = "drop"
    )

  p <- ggplot(plot_dat, aes(x = threshold_rule, y = pct, fill = group)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    geom_text(
      aes(label = sprintf("%.1f", pct)),
      position = position_dodge(width = 0.75),
      vjust = -0.2,
      size = 3.2
    ) +
    geom_text(
      data = ann_dat,
      aes(x = threshold_rule, y = y_pos, label = star_label),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 5
    ) +
    facet_grid(clade_group ~ exposure_type, scales = "free_y") +
    scale_fill_manual(values = c("Yes" = "#e15759", "No" = "#4e79a7"), name = NULL) +
    labs(
      title = "Restricted threshold enrichment for direct overlap",
      subtitle = "Restricted comparison set removes pairs with no plausible epidemiologic opportunity",
      x = NULL,
      y = "Percent below threshold"
    ) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(viz_dir, "restricted_location_overlap_threshold_enrichment.png"),
    plot = p,
    width = 14,
    height = 8,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  save_plot_pdf(
    plot_obj = p,
    filename = file.path(viz_dir, "restricted_location_overlap_threshold_enrichment.pdf"),
    width = 14,
    height = 8
  )
}

clade1_pairs <- prepare_snv_pairs(clade1_snv_path, "Clade 1")
clade2_pairs <- prepare_snv_pairs(clade2_snv_path, "Clade 2")

clade1_flags <- compute_restricted_flags(clade1_pairs, "Clade 1")
clade2_flags <- compute_restricted_flags(clade2_pairs, "Clade 2")

annotated_pairs <- bind_rows(
  clade1_pairs %>% left_join(clade1_flags, by = c("patient_a", "patient_b")),
  clade2_pairs %>% left_join(clade2_flags, by = c("patient_a", "patient_b"))
) %>%
  mutate(across(c(hospital_time_overlap, facility_overlap, floor_overlap, unit_overlap, room_overlap), ~ replace_na(.x, FALSE)))

plot_dat <- build_restricted_plot_data(annotated_pairs)
threshold_dat <- run_threshold_tests(plot_dat)

write.csv(
  threshold_dat,
  file.path(desc_dir, "restricted_location_overlap_threshold_enrichment_fisher_tests.csv"),
  row.names = FALSE
)
write.csv(
  annotated_pairs,
  file.path(desc_dir, "restricted_location_overlap_annotated_pairs.csv"),
  row.names = FALSE
)

writeLines(
  c(
    "Restricted comparison set analysis uses between-patient pairs only.",
    paste0("Special units excluded before building the restriction set: ", paste(special_units, collapse = ", "), "."),
    "All location-level comparisons are restricted to pairs with any overlapping admission stay anywhere in the hospital.",
    "Direct overlap uses full timestamps."
  ),
  con = file.path(desc_dir, "restricted_location_overlap_analysis_notes.txt")
)

make_distribution_plot(plot_dat)
make_threshold_plot(threshold_dat)

cat("Done. Restricted location overlap outputs written to:\n", desc_dir, "\n", sep = "")
