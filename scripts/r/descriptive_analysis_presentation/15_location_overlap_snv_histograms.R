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

sequential_gap_days <- 14
sequential_suffix <- "14d"

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
    Bed = trimws(as.character(Bed)),
    Date_date = parse_date_safe(Date),
    StartDate_dt = coalesce(parse_datetime_safe(StartDate), as.POSIXct(Date_date, tz = "UTC")),
    EndDate_dt = coalesce(parse_datetime_safe(EndDate), parse_datetime_safe(StartDate), as.POSIXct(Date_date, tz = "UTC"))
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
      patient_a = if_else(
        !is.na(patient_id_1) & !is.na(patient_id_2),
        pmin(patient_id_1, patient_id_2),
        NA_character_
      ),
      patient_b = if_else(
        !is.na(patient_id_1) & !is.na(patient_id_2),
        pmax(patient_id_1, patient_id_2),
        NA_character_
      )
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
    select(
      PatientID = PatientID,
      StartDate_dt = StartDate_dt,
      EndDate_dt = EndDate_dt,
      all_of(key_fields)
    ) %>%
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

find_sequential_pairs_for_level <- function(loc_df, key_fields, max_gap_days = 60) {
  dat <- loc_df %>%
    filter(if_all(all_of(key_fields), is_present)) %>%
    select(
      PatientID = PatientID,
      StartDate_dt = StartDate_dt,
      EndDate_dt = EndDate_dt,
      all_of(key_fields)
    ) %>%
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

      while (j <= n && group_df$StartDate_dt[j] <= group_df$EndDate_dt[i] + as.difftime(max_gap_days, units = "days")) {
        gap_days <- as.numeric(difftime(group_df$StartDate_dt[j], group_df$EndDate_dt[i], units = "days"))

        if (
          group_df$PatientID[i] != group_df$PatientID[j] &&
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
  }

  pair_keys <- ls(pair_keys_env, all.names = TRUE)

  if (length(pair_keys) == 0) {
    return(tibble(patient_a = character(), patient_b = character()))
  }

  tibble(pair_key = pair_keys) %>%
    separate(pair_key, into = c("patient_a", "patient_b"), sep = "\\|\\|", remove = TRUE)
}

compute_overlap_flags <- function(snv_pairs, clade_label) {
  loc_clade <- locations %>%
    filter(clade_group == clade_label)

  pair_flags <- snv_pairs %>%
    filter(same_patient == "Between patient", !is.na(patient_a), !is.na(patient_b)) %>%
    distinct(patient_a, patient_b)

  if (nrow(pair_flags) == 0) {
    return(tibble())
  }

  facility_pairs <- find_overlap_pairs_for_level(loc_clade, c("FacilityCode")) %>%
    mutate(facility_overlap = TRUE)
  floor_pairs <- find_overlap_pairs_for_level(loc_clade, c("floor")) %>%
    mutate(floor_overlap = TRUE)
  unit_pairs <- find_overlap_pairs_for_level(loc_clade, c("FacilityCode", "Unit")) %>%
    mutate(unit_overlap = TRUE)
  room_pairs <- find_overlap_pairs_for_level(loc_clade, c("FacilityCode", "Room")) %>%
    mutate(room_overlap = TRUE)
  bed_pairs <- find_overlap_pairs_for_level(loc_clade, c("FacilityCode", "Room", "Bed")) %>%
    mutate(bed_overlap = TRUE)
  facility_seq_pairs <- find_sequential_pairs_for_level(loc_clade, c("FacilityCode"), max_gap_days = sequential_gap_days) %>%
    mutate(facility_sequential_14d = TRUE)
  floor_seq_pairs <- find_sequential_pairs_for_level(loc_clade, c("floor"), max_gap_days = sequential_gap_days) %>%
    mutate(floor_sequential_14d = TRUE)
  unit_seq_pairs <- find_sequential_pairs_for_level(loc_clade, c("FacilityCode", "Unit"), max_gap_days = sequential_gap_days) %>%
    mutate(unit_sequential_14d = TRUE)
  room_seq_pairs <- find_sequential_pairs_for_level(loc_clade, c("FacilityCode", "Room"), max_gap_days = sequential_gap_days) %>%
    mutate(room_sequential_14d = TRUE)
  bed_seq_pairs <- find_sequential_pairs_for_level(loc_clade, c("FacilityCode", "Room", "Bed"), max_gap_days = sequential_gap_days) %>%
    mutate(bed_sequential_14d = TRUE)

  pair_flags %>%
    left_join(facility_pairs, by = c("patient_a", "patient_b")) %>%
    left_join(floor_pairs, by = c("patient_a", "patient_b")) %>%
    left_join(unit_pairs, by = c("patient_a", "patient_b")) %>%
    left_join(room_pairs, by = c("patient_a", "patient_b")) %>%
    left_join(bed_pairs, by = c("patient_a", "patient_b")) %>%
    left_join(facility_seq_pairs, by = c("patient_a", "patient_b")) %>%
    left_join(floor_seq_pairs, by = c("patient_a", "patient_b")) %>%
    left_join(unit_seq_pairs, by = c("patient_a", "patient_b")) %>%
    left_join(room_seq_pairs, by = c("patient_a", "patient_b")) %>%
    left_join(bed_seq_pairs, by = c("patient_a", "patient_b")) %>%
    mutate(
      across(matches("(_overlap|_sequential_14d)$"), ~ replace_na(.x, FALSE))
    )
}

make_exposure_plots <- function(dat, clade_label, exposure_mode = c("overlap", "sequential")) {
  exposure_mode <- match.arg(exposure_mode)
  selected_cols <- if (exposure_mode == "overlap") {
    c("facility_overlap", "floor_overlap", "unit_overlap", "room_overlap")
  } else {
    c("facility_sequential_14d", "floor_sequential_14d", "unit_sequential_14d", "room_sequential_14d")
  }

  plot_long <- dat %>%
    filter(same_patient == "Between patient") %>%
    select(snv, all_of(selected_cols)) %>%
    pivot_longer(
      cols = all_of(selected_cols),
      names_to = "exposure_type",
      values_to = "exposure_flag"
    ) %>%
    mutate(
      exposure_type = recode(
        exposure_type,
        facility_overlap = "Facility",
        floor_overlap = "Floor",
        unit_overlap = "Unit",
        room_overlap = "Room",
        facility_sequential_14d = "Facility",
        floor_sequential_14d = "Floor",
        unit_sequential_14d = "Unit",
        room_sequential_14d = "Room"
      ),
      exposure_type = factor(exposure_type, levels = c("Facility", "Floor", "Unit", "Room", "Bed")),
      exposure_status = if (exposure_mode == "overlap") {
        if_else(exposure_flag, "Overlap yes", "Overlap no")
      } else {
        if_else(exposure_flag, "Sequential yes", "Sequential no")
      }
    )

  title_stub <- if (exposure_mode == "overlap") "location overlap" else paste0("location sequential exposure within ", sequential_gap_days, " days")
  subtitle_stub <- if (exposure_mode == "overlap") {
    paste0("Special units excluded: ", paste(special_units, collapse = ", "))
  } else {
    paste0("Sequential = same location, no direct overlap, gap 1-", sequential_gap_days, " days; special units excluded: ", paste(special_units, collapse = ", "))
  }
  fill_values <- if (exposure_mode == "overlap") {
    c("Overlap yes" = "#e15759", "Overlap no" = "#4e79a7")
  } else {
    c("Sequential yes" = "#f28e2b", "Sequential no" = "#4e79a7")
  }

  plot_log <- plot_long %>%
    filter(!is.na(snv), snv > 0) %>%
    ggplot(aes(x = snv, fill = exposure_status)) +
    geom_histogram(
      bins = 100,
      position = "identity",
      alpha = 0.5
    ) +
    scale_x_log10() +
    scale_fill_manual(
      values = fill_values,
      name = NULL
    ) +
    facet_wrap(~exposure_type, scales = "free_y", ncol = 2) +
    labs(
      title = paste0(clade_label, " - between-patient SNV by ", title_stub),
      subtitle = subtitle_stub,
      x = "SNV distance (log scale)",
      y = "Pairs"
    ) +
    theme_minimal(base_size = 32) +
    theme(
      plot.title = element_text(face = "bold", size = 38),
      plot.subtitle = element_text(size = 24),
      strip.text = element_text(size = 26, face = "bold"),
      axis.text = element_text(size = 22),
      axis.title = element_text(size = 24),
      legend.text = element_text(size = 22),
      legend.position = "bottom"
    )

  plot_zoom <- plot_long %>%
    filter(!is.na(snv), snv >= 0, snv <= 75) %>%
    ggplot(aes(x = snv, fill = exposure_status)) +
    geom_histogram(
      bins = 40,
      closed = "right",
      position = "identity",
      alpha = 0.5
    ) +
    scale_fill_manual(
      values = fill_values,
      name = NULL
    ) +
    scale_x_continuous(
      limits = c(0, 75)
    ) +
    facet_wrap(~exposure_type, scales = "free_y", ncol = 2) +
    labs(
      title = paste0(clade_label, " - between-patient SNV by ", title_stub, " (0-75 SNVs)"),
      subtitle = subtitle_stub,
      x = "SNV distance",
      y = "Pairs"
    ) +
    theme_minimal(base_size = 32) +
    theme(
      plot.title = element_text(face = "bold", size = 38),
      plot.subtitle = element_text(size = 24),
      strip.text = element_text(size = 26, face = "bold"),
      axis.text = element_text(size = 22),
      axis.title = element_text(size = 24),
      legend.text = element_text(size = 22),
      legend.position = "bottom"
    )

  list(log_plot = plot_log, zoom_plot = plot_zoom, long_data = plot_long)
}

clade1_pairs <- prepare_snv_pairs(clade1_snv_path, "Clade 1")
clade2_pairs <- prepare_snv_pairs(clade2_snv_path, "Clade 2")

clade1_flags <- compute_overlap_flags(clade1_pairs, "Clade 1")
clade2_flags <- compute_overlap_flags(clade2_pairs, "Clade 2")

clade1_pairs_annotated <- clade1_pairs %>%
  left_join(clade1_flags, by = c("patient_a", "patient_b"))

clade2_pairs_annotated <- clade2_pairs %>%
  left_join(clade2_flags, by = c("patient_a", "patient_b"))

clade1_overlap_plots <- make_exposure_plots(clade1_pairs_annotated, "Clade 1", exposure_mode = "overlap")
clade2_overlap_plots <- make_exposure_plots(clade2_pairs_annotated, "Clade 2", exposure_mode = "overlap")
clade1_sequential_plots <- make_exposure_plots(clade1_pairs_annotated, "Clade 1", exposure_mode = "sequential")
clade2_sequential_plots <- make_exposure_plots(clade2_pairs_annotated, "Clade 2", exposure_mode = "sequential")

overlap_summary <- bind_rows(
  clade1_overlap_plots$long_data %>% mutate(clade_group = "Clade 1", exposure_class = "Overlap"),
  clade2_overlap_plots$long_data %>% mutate(clade_group = "Clade 2", exposure_class = "Overlap"),
  clade1_sequential_plots$long_data %>% mutate(clade_group = "Clade 1", exposure_class = "Sequential_14d"),
  clade2_sequential_plots$long_data %>% mutate(clade_group = "Clade 2", exposure_class = "Sequential_14d")
) %>%
  group_by(clade_group, exposure_class, exposure_type, exposure_status) %>%
  summarise(pair_n = n(), .groups = "drop")

write.csv(
  clade1_pairs_annotated,
  file.path(desc_dir, "clade1_post_noncore_snv_long_with_location_overlap.csv"),
  row.names = FALSE
)
write.csv(
  clade2_pairs_annotated,
  file.path(desc_dir, "clade2_post_noncore_snv_long_with_location_overlap.csv"),
  row.names = FALSE
)
write.csv(
  overlap_summary,
  file.path(desc_dir, "location_overlap_pair_summary_by_clade.csv"),
  row.names = FALSE
)

writeLines(
  c(
    "Overlap analysis uses only the filtered clade 1/clade 2 cohort from metadata_case_unformed_clade1_clade2.csv.",
    "Overlap is evaluated for between-patient SNV pairs only.",
    paste0("Special units excluded before overlap calculations: ", paste(special_units, collapse = ", "), "."),
    "A pair is overlap yes at a location level when both patients have the same location value and their full StartDate-EndDate timestamps overlap.",
    "Same-day but non-overlapping timestamps are treated as sequential rather than direct overlap.",
    "Sequential 14-day exposure means the two patients share the same location value, do not directly overlap in time there, and the later stay begins within 14 days after the earlier stay ends.",
    "floor overlap uses the derived floor label from the filtered location dataset."
  ),
  con = file.path(desc_dir, "location_overlap_analysis_notes.txt")
)

ggsave(
  filename = file.path(viz_dir, "clade1_between_patient_location_overlap_snv_log.png"),
  plot = clade1_overlap_plots$log_plot,
  width = 22,
  height = 14,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = clade1_overlap_plots$log_plot,
  filename = file.path(viz_dir, "clade1_between_patient_location_overlap_snv_log.pdf"),
  width = 22,
  height = 14
)

ggsave(
  filename = file.path(viz_dir, "clade1_between_patient_location_overlap_snv_0_75.png"),
  plot = clade1_overlap_plots$zoom_plot,
  width = 22,
  height = 14,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = clade1_overlap_plots$zoom_plot,
  filename = file.path(viz_dir, "clade1_between_patient_location_overlap_snv_0_75.pdf"),
  width = 22,
  height = 14
)

ggsave(
  filename = file.path(viz_dir, "clade2_between_patient_location_overlap_snv_log.png"),
  plot = clade2_overlap_plots$log_plot,
  width = 22,
  height = 14,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = clade2_overlap_plots$log_plot,
  filename = file.path(viz_dir, "clade2_between_patient_location_overlap_snv_log.pdf"),
  width = 22,
  height = 14
)

ggsave(
  filename = file.path(viz_dir, "clade2_between_patient_location_overlap_snv_0_75.png"),
  plot = clade2_overlap_plots$zoom_plot,
  width = 22,
  height = 14,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = clade2_overlap_plots$zoom_plot,
  filename = file.path(viz_dir, "clade2_between_patient_location_overlap_snv_0_75.pdf"),
  width = 22,
  height = 14
)

ggsave(
  filename = file.path(viz_dir, "clade1_between_patient_location_sequential_14d_snv_log.png"),
  plot = clade1_sequential_plots$log_plot,
  width = 22,
  height = 14,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = clade1_sequential_plots$log_plot,
  filename = file.path(viz_dir, "clade1_between_patient_location_sequential_14d_snv_log.pdf"),
  width = 22,
  height = 14
)

ggsave(
  filename = file.path(viz_dir, "clade1_between_patient_location_sequential_14d_snv_0_75.png"),
  plot = clade1_sequential_plots$zoom_plot,
  width = 22,
  height = 14,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = clade1_sequential_plots$zoom_plot,
  filename = file.path(viz_dir, "clade1_between_patient_location_sequential_14d_snv_0_75.pdf"),
  width = 22,
  height = 14
)

ggsave(
  filename = file.path(viz_dir, "clade2_between_patient_location_sequential_14d_snv_log.png"),
  plot = clade2_sequential_plots$log_plot,
  width = 22,
  height = 14,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = clade2_sequential_plots$log_plot,
  filename = file.path(viz_dir, "clade2_between_patient_location_sequential_14d_snv_log.pdf"),
  width = 22,
  height = 14
)

ggsave(
  filename = file.path(viz_dir, "clade2_between_patient_location_sequential_14d_snv_0_75.png"),
  plot = clade2_sequential_plots$zoom_plot,
  width = 22,
  height = 14,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = clade2_sequential_plots$zoom_plot,
  filename = file.path(viz_dir, "clade2_between_patient_location_sequential_14d_snv_0_75.pdf"),
  width = 22,
  height = 14
)

cat("Done. Location overlap SNV outputs written to:\n", viz_dir, "\n", sep = "")
