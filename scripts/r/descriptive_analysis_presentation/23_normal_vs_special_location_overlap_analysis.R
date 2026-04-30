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
    is_present(Unit)
  ) %>%
  mutate(unit_type_group = if_else(Unit %in% special_units, "Special", "Normal"))

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

compute_group_flags <- function(snv_pairs, clade_label, unit_group_label) {
  loc_clade <- locations %>%
    filter(clade_group == clade_label, unit_type_group == unit_group_label)

  pair_flags <- snv_pairs %>%
    filter(same_patient == "Between patient", !is.na(patient_a), !is.na(patient_b)) %>%
    distinct(patient_a, patient_b)

  facility_pairs <- find_overlap_pairs_for_level(loc_clade, c("FacilityCode")) %>% mutate(facility_overlap = TRUE)
  floor_pairs <- find_overlap_pairs_for_level(loc_clade, c("floor")) %>% mutate(floor_overlap = TRUE)
  unit_pairs <- find_overlap_pairs_for_level(loc_clade, c("FacilityCode", "Unit")) %>% mutate(unit_overlap = TRUE)
  room_pairs <- find_overlap_pairs_for_level(loc_clade, c("FacilityCode", "Room")) %>% mutate(room_overlap = TRUE)

  pair_flags %>%
    left_join(facility_pairs, by = c("patient_a", "patient_b")) %>%
    left_join(floor_pairs, by = c("patient_a", "patient_b")) %>%
    left_join(unit_pairs, by = c("patient_a", "patient_b")) %>%
    left_join(room_pairs, by = c("patient_a", "patient_b")) %>%
    mutate(across(ends_with("_overlap"), ~ replace_na(.x, FALSE)))
}

build_plot_long <- function(dat, group_label) {
  dat %>%
    filter(same_patient == "Between patient", !is.na(snv), snv >= 0, snv <= 75) %>%
    select(clade_group, snv, facility_overlap, floor_overlap, unit_overlap, room_overlap) %>%
    pivot_longer(
      cols = c(facility_overlap, floor_overlap, unit_overlap, room_overlap),
      names_to = "exposure_type",
      values_to = "exposure_flag"
    ) %>%
    mutate(
      unit_type_group = group_label,
      exposure_type = recode(
        exposure_type,
        facility_overlap = "Facility",
        floor_overlap = "Floor",
        unit_overlap = "Unit",
        room_overlap = "Room"
      ),
      exposure_type = factor(exposure_type, levels = c("Facility", "Floor", "Unit", "Room")),
      exposure_status = if_else(exposure_flag, "Overlap yes", "Overlap no")
    )
}

run_threshold_tests <- function(dat, group_label) {
  thresholds <- c(2, 5, 10)

  crossing(
    clade_group = c("Clade 1", "Clade 2"),
    exposure_type = c("Facility", "Floor", "Unit", "Room"),
    threshold = thresholds
  ) %>%
    rowwise() %>%
    mutate(
      group_dat = list(dat %>% filter(clade_group == .data$clade_group, exposure_type == .data$exposure_type)),
      yes_low_n = sum(group_dat$snv[group_dat$exposure_status == "Overlap yes"] <= threshold, na.rm = TRUE),
      yes_high_n = sum(group_dat$snv[group_dat$exposure_status == "Overlap yes"] > threshold, na.rm = TRUE),
      no_low_n = sum(group_dat$snv[group_dat$exposure_status == "Overlap no"] <= threshold, na.rm = TRUE),
      no_high_n = sum(group_dat$snv[group_dat$exposure_status == "Overlap no"] > threshold, na.rm = TRUE),
      yes_low_pct = ifelse((yes_low_n + yes_high_n) > 0, 100 * yes_low_n / (yes_low_n + yes_high_n), NA_real_),
      no_low_pct = ifelse((no_low_n + no_high_n) > 0, 100 * no_low_n / (no_low_n + no_high_n), NA_real_),
      fisher_p_value = safe_fisher(yes_low_n, yes_high_n, no_low_n, no_high_n),
      p_value_label = format_p(fisher_p_value)
    ) %>%
    ungroup() %>%
    transmute(
      unit_type_group = group_label,
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

make_distribution_plot <- function(plot_dat, group_label, file_stub) {
  p <- ggplot(plot_dat, aes(x = snv, fill = exposure_status)) +
    geom_histogram(
      bins = 40,
      position = "identity",
      alpha = 0.5
    ) +
    scale_fill_manual(values = c("Overlap yes" = "#e15759", "Overlap no" = "#4e79a7"), name = NULL) +
    scale_x_continuous(limits = c(0, 75)) +
    facet_grid(clade_group ~ exposure_type, scales = "free_y") +
    labs(
      title = paste0("Between-patient SNV by location overlap: ", group_label, " unit rows"),
      subtitle = "Direct overlap uses full timestamps",
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
    filename = file.path(viz_dir, paste0(file_stub, ".png")),
    plot = p,
    width = 14,
    height = 8,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  save_plot_pdf(
    plot_obj = p,
    filename = file.path(viz_dir, paste0(file_stub, ".pdf")),
    width = 14,
    height = 8
  )
}

make_threshold_plot <- function(threshold_dat, group_label, file_stub) {
  plot_dat <- bind_rows(
    threshold_dat %>%
      transmute(clade_group, exposure_type, threshold_rule, group = "Yes", pct = yes_low_pct, p_value_label),
    threshold_dat %>%
      transmute(clade_group, exposure_type, threshold_rule, group = "No", pct = no_low_pct, p_value_label)
  ) %>%
    mutate(
      group = factor(group, levels = c("Yes", "No")),
      threshold_rule = factor(threshold_rule, levels = c("SNV <= 2", "SNV <= 5", "SNV <= 10")),
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
      title = paste0("Threshold enrichment for direct overlap: ", group_label, " unit rows"),
      subtitle = "Bars show percent of pairwise SNV rows below threshold",
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
    filename = file.path(viz_dir, paste0(file_stub, ".png")),
    plot = p,
    width = 14,
    height = 8,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  save_plot_pdf(
    plot_obj = p,
    filename = file.path(viz_dir, paste0(file_stub, ".pdf")),
    width = 14,
    height = 8
  )
}

clade1_pairs <- prepare_snv_pairs(clade1_snv_path, "Clade 1")
clade2_pairs <- prepare_snv_pairs(clade2_snv_path, "Clade 2")

clade1_normal_flags <- compute_group_flags(clade1_pairs, "Clade 1", "Normal")
clade2_normal_flags <- compute_group_flags(clade2_pairs, "Clade 2", "Normal")
clade1_special_flags <- compute_group_flags(clade1_pairs, "Clade 1", "Special")
clade2_special_flags <- compute_group_flags(clade2_pairs, "Clade 2", "Special")

normal_pairs <- bind_rows(
  clade1_pairs %>% left_join(clade1_normal_flags, by = c("patient_a", "patient_b")),
  clade2_pairs %>% left_join(clade2_normal_flags, by = c("patient_a", "patient_b"))
) %>%
  mutate(across(ends_with("_overlap"), ~ replace_na(.x, FALSE)))

special_pairs <- bind_rows(
  clade1_pairs %>% left_join(clade1_special_flags, by = c("patient_a", "patient_b")),
  clade2_pairs %>% left_join(clade2_special_flags, by = c("patient_a", "patient_b"))
) %>%
  mutate(across(ends_with("_overlap"), ~ replace_na(.x, FALSE)))

normal_plot_dat <- build_plot_long(normal_pairs, "Normal")
special_plot_dat <- build_plot_long(special_pairs, "Special")

normal_threshold <- run_threshold_tests(normal_plot_dat, "Normal")
special_threshold <- run_threshold_tests(special_plot_dat, "Special")

write.csv(normal_threshold, file.path(desc_dir, "normal_unit_rows_threshold_enrichment_fisher_tests.csv"), row.names = FALSE)
write.csv(special_threshold, file.path(desc_dir, "special_unit_rows_threshold_enrichment_fisher_tests.csv"), row.names = FALSE)

writeLines(
  c(
    "This analysis splits the location-overlap workflow into rows from normal units versus rows from special units.",
    paste0("Special units are defined as: ", paste(special_units, collapse = ", "), "."),
    "Within each subset, direct overlap is recalculated separately for Facility, Floor, Unit, and Room using full timestamps.",
    "Plots are restricted to SNV 0-75 for readability."
  ),
  con = file.path(desc_dir, "normal_vs_special_location_overlap_analysis_notes.txt")
)

make_distribution_plot(normal_plot_dat, "Normal", "normal_unit_rows_location_overlap_snv_0_75")
make_distribution_plot(special_plot_dat, "Special", "special_unit_rows_location_overlap_snv_0_75")
make_threshold_plot(normal_threshold, "Normal", "normal_unit_rows_threshold_enrichment_overlap")
make_threshold_plot(special_threshold, "Special", "special_unit_rows_threshold_enrichment_overlap")

cat("Done. Normal-vs-special location overlap outputs written to:\n", desc_dir, "\n", sep = "")
