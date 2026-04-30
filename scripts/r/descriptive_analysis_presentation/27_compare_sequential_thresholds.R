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
viz_dir <- file.path(output_dir, "visualization")

min_snv_path <- file.path(desc_dir, "all_between_patient_min_snv_by_patient_pair.csv")
locations_path <- file.path(desc_dir, "locations_filtered_case_unformed_clade1_clade2.tsv")

if (!dir.exists(viz_dir)) {
  dir.create(viz_dir, recursive = TRUE)
}

required_pkgs <- c("ggplot2", "dplyr", "tidyr", "readr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
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

safe_fisher <- function(yes_low, yes_high, no_low, no_high) {
  mat <- matrix(c(yes_low, yes_high, no_low, no_high), nrow = 2, byrow = TRUE)
  if (any(rowSums(mat) == 0) || any(colSums(mat) == 0)) {
    return(NA_real_)
  }
  fisher.test(mat)$p.value
}

format_p <- function(x) {
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    x < 0.001 ~ "<0.001",
    TRUE ~ formatC(x, format = "f", digits = 3)
  )
}

make_star_label <- function(p_label) {
  dplyr::case_when(
    is.na(p_label) ~ "ns",
    p_label == "<0.001" ~ "***",
    suppressWarnings(as.numeric(p_label)) < 0.01 ~ "**",
    suppressWarnings(as.numeric(p_label)) < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

collapse_intervals <- function(starts, ends) {
  ord <- order(starts, ends)
  starts <- starts[ord]
  ends <- ends[ord]
  out_starts <- vector("list", 0)
  out_ends <- vector("list", 0)
  cur_start <- starts[1]
  cur_end <- ends[1]

  if (length(starts) > 1) {
    for (i in 2:length(starts)) {
      if (starts[i] <= cur_end) {
        cur_end <- max(cur_end, ends[i])
      } else {
        out_starts[[length(out_starts) + 1]] <- cur_start
        out_ends[[length(out_ends) + 1]] <- cur_end
        cur_start <- starts[i]
        cur_end <- ends[i]
      }
    }
  }

  out_starts[[length(out_starts) + 1]] <- cur_start
  out_ends[[length(out_ends) + 1]] <- cur_end

  data.frame(
    StartDate_dt = as.POSIXct(unlist(out_starts), origin = "1970-01-01", tz = "UTC"),
    EndDate_dt = as.POSIXct(unlist(out_ends), origin = "1970-01-01", tz = "UTC"),
    stringsAsFactors = FALSE
  )
}

compute_min_gap_pairs <- function(loc_df, key_cols, max_gap_days = 60) {
  keep <- rep(TRUE, nrow(loc_df))
  for (col in key_cols) {
    keep <- keep & is_present(loc_df[[col]])
  }
  dat <- loc_df[keep, c("PatientID", "StartDate_dt", "EndDate_dt", key_cols), drop = FALSE]
  if (nrow(dat) == 0) {
    return(data.frame(patient_a = character(), patient_b = character(), min_gap_days = numeric(), stringsAsFactors = FALSE))
  }

  dat$location_key <- do.call(paste, c(dat[key_cols], sep = "||"))
  dat <- unique(dat[, c("PatientID", "StartDate_dt", "EndDate_dt", "location_key"), drop = FALSE])

  combo_key <- paste(dat$PatientID, dat$location_key, sep = "___")
  split_idx <- split(seq_len(nrow(dat)), combo_key)

  collapsed_list <- vector("list", length(split_idx))
  nm <- names(split_idx)
  for (i in seq_along(split_idx)) {
    rows <- split_idx[[i]]
    tmp <- collapse_intervals(dat$StartDate_dt[rows], dat$EndDate_dt[rows])
    parts <- strsplit(nm[i], "___", fixed = TRUE)[[1]]
    tmp$PatientID <- parts[1]
    tmp$location_key <- parts[2]
    collapsed_list[[i]] <- tmp
  }
  collapsed <- do.call(rbind, collapsed_list)
  if (is.null(collapsed) || nrow(collapsed) == 0) {
    return(data.frame(patient_a = character(), patient_b = character(), min_gap_days = numeric(), stringsAsFactors = FALSE))
  }

  loc_split <- split(collapsed, collapsed$location_key)
  pair_gap_env <- new.env(parent = emptyenv(), hash = TRUE)
  max_gap_dt <- as.difftime(max_gap_days, units = "days")

  for (group_df in loc_split) {
    if (nrow(group_df) < 2) {
      next
    }
    ord <- order(group_df$StartDate_dt, group_df$EndDate_dt)
    group_df <- group_df[ord, , drop = FALSE]
    n <- nrow(group_df)

    for (i in seq_len(n - 1)) {
      j <- i + 1
      while (j <= n && group_df$StartDate_dt[j] <= group_df$EndDate_dt[i] + max_gap_dt) {
        gap_days <- as.numeric(difftime(group_df$StartDate_dt[j], group_df$EndDate_dt[i], units = "days"))
        if (
          group_df$PatientID[i] != group_df$PatientID[j] &&
          is.finite(gap_days) &&
          gap_days > 0 &&
          gap_days <= max_gap_days
        ) {
          pair_a <- min(group_df$PatientID[i], group_df$PatientID[j])
          pair_b <- max(group_df$PatientID[i], group_df$PatientID[j])
          pair_key <- paste(pair_a, pair_b, sep = "||")
          if (!exists(pair_key, envir = pair_gap_env, inherits = FALSE) ||
              gap_days < get(pair_key, envir = pair_gap_env, inherits = FALSE)) {
            assign(pair_key, gap_days, envir = pair_gap_env)
          }
        }
        j <- j + 1
      }
    }
  }

  pair_keys <- ls(pair_gap_env, all.names = TRUE)
  if (length(pair_keys) == 0) {
    return(data.frame(patient_a = character(), patient_b = character(), min_gap_days = numeric(), stringsAsFactors = FALSE))
  }

  out <- data.frame(
    pair_key = pair_keys,
    min_gap_days = vapply(pair_keys, function(x) get(x, envir = pair_gap_env, inherits = FALSE), numeric(1)),
    stringsAsFactors = FALSE
  )
  parts <- strsplit(out$pair_key, "\\|\\|")
  out$patient_a <- vapply(parts, `[`, character(1), 1)
  out$patient_b <- vapply(parts, `[`, character(1), 2)
  out[, c("patient_a", "patient_b", "min_gap_days"), drop = FALSE]
}

special_units <- c("ESA", "OR", "DVU", "MPUA", "UADM", "USSS", "ORC", "VHOR", "IRU", "CATH", "MHOR", "INU")
threshold_values <- c(7, 14, 30, 60)
location_levels <- c("Facility", "Floor", "Unit", "Room")
location_colors <- c("Facility" = "#4E79A7", "Floor" = "#59A14F", "Unit" = "#F28E2B", "Room" = "#E15759")

min_snv_pairs <- read_csv(min_snv_path, show_col_types = FALSE) %>%
  transmute(
    clade_group,
    patient_a = as.character(patient_a),
    patient_b = as.character(patient_b),
    snv = as.numeric(min_snv)
  ) %>%
  filter(!is.na(patient_a), !is.na(patient_b), !is.na(snv))

locations <- read_delim(locations_path, delim = "\t", show_col_types = FALSE) %>%
  mutate(
    PatientID = trimws(as.character(PatientID)),
    clade_group = trimws(as.character(clade_group)),
    FacilityCode = trimws(as.character(FacilityCode)),
    floor = trimws(as.character(floor)),
    Unit = trimws(as.character(Unit)),
    Room = trimws(as.character(Room)),
    StartDate_dt = parse_datetime_safe(StartDate),
    EndDate_dt = dplyr::coalesce(parse_datetime_safe(EndDate), StartDate_dt)
  ) %>%
  filter(
    !is.na(PatientID),
    PatientID != "",
    !is.na(StartDate_dt),
    !is.na(EndDate_dt),
    !(Unit %in% special_units)
  )

compute_gap_tables_for_clade <- function(clade_label) {
  cat("Preparing gap tables for", clade_label, "...\n")
  relevant_pairs <- min_snv_pairs[min_snv_pairs$clade_group == clade_label, c("patient_a", "patient_b"), drop = FALSE]
  relevant_patients <- unique(c(relevant_pairs$patient_a, relevant_pairs$patient_b))
  loc_clade <- locations[locations$clade_group == clade_label & locations$PatientID %in% relevant_patients, , drop = FALSE]

  cat("  Facility...\n")
  facility <- compute_min_gap_pairs(loc_clade, c("FacilityCode"), 60)
  cat("  Floor...\n")
  floor <- compute_min_gap_pairs(loc_clade, c("floor"), 60)
  cat("  Unit...\n")
  unit <- compute_min_gap_pairs(loc_clade, c("FacilityCode", "Unit"), 60)
  cat("  Room...\n")
  room <- compute_min_gap_pairs(loc_clade, c("FacilityCode", "Room"), 60)

  list(
    facility = facility,
    floor = floor,
    unit = unit,
    room = room
  )
}

annotate_threshold <- function(pair_df, clade_label, day_threshold, gap_tables) {
  clade_pairs <- pair_df %>%
    filter(clade_group == clade_label) %>%
    distinct(patient_a, patient_b)

  facility_pairs <- subset(gap_tables$facility, min_gap_days <= day_threshold, select = c("patient_a", "patient_b"))
  floor_pairs <- subset(gap_tables$floor, min_gap_days <= day_threshold, select = c("patient_a", "patient_b"))
  unit_pairs <- subset(gap_tables$unit, min_gap_days <= day_threshold, select = c("patient_a", "patient_b"))
  room_pairs <- subset(gap_tables$room, min_gap_days <= day_threshold, select = c("patient_a", "patient_b"))

  names(facility_pairs)[names(facility_pairs) == "patient_a"] <- "patient_a"
  facility_pairs$facility_seq <- TRUE
  floor_pairs$floor_seq <- TRUE
  unit_pairs$unit_seq <- TRUE
  room_pairs$room_seq <- TRUE

  clade_pairs %>%
    left_join(facility_pairs, by = c("patient_a", "patient_b")) %>%
    left_join(floor_pairs, by = c("patient_a", "patient_b")) %>%
    left_join(unit_pairs, by = c("patient_a", "patient_b")) %>%
    left_join(room_pairs, by = c("patient_a", "patient_b")) %>%
    mutate(across(ends_with("_seq"), ~ replace_na(.x, FALSE)))
}

make_long <- function(dat, threshold_days) {
  dat %>%
    transmute(
      clade_group,
      snv,
      threshold_days,
      facility_seq,
      floor_seq,
      unit_seq,
      room_seq
    ) %>%
    pivot_longer(
      cols = c(facility_seq, floor_seq, unit_seq, room_seq),
      names_to = "exposure_type",
      values_to = "seq_flag"
    ) %>%
    mutate(
      exposure_type = recode(
        exposure_type,
        facility_seq = "Facility",
        floor_seq = "Floor",
        unit_seq = "Unit",
        room_seq = "Room"
      ),
      exposure_type = factor(exposure_type, levels = location_levels),
      threshold_days = factor(threshold_days, levels = threshold_values),
      seq_status = if_else(seq_flag, "Sequential yes", "Sequential no")
    )
}

summarize_threshold <- function(dat, threshold_days) {
  test_rows <- list()
  enrich_rows <- list()
  row_i <- 1
  enrich_i <- 1

  for (clade_label in c("Clade 1", "Clade 2")) {
    dat_clade <- dat[dat$clade_group == clade_label, , drop = FALSE]
    for (exposure_type in location_levels) {
      flag_col <- paste0(tolower(exposure_type), "_seq")
      yes_vals <- dat_clade$snv[dat_clade[[flag_col]]]
      no_vals <- dat_clade$snv[!dat_clade[[flag_col]]]
      p_val <- safe_wilcox(yes_vals, no_vals)

      test_rows[[row_i]] <- data.frame(
        clade_group = clade_label,
        exposure_type = exposure_type,
        threshold_days = threshold_days,
        yes_n = length(yes_vals),
        no_n = length(no_vals),
        yes_median_snv = ifelse(length(yes_vals) > 0, median(yes_vals), NA_real_),
        no_median_snv = ifelse(length(no_vals) > 0, median(no_vals), NA_real_),
        median_diff = ifelse(length(yes_vals) > 0 && length(no_vals) > 0, median(yes_vals) - median(no_vals), NA_real_),
        wilcox_p_value = p_val,
        p_value_label = format_p(p_val),
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1

      for (cutoff in c(2, 5, 10)) {
        yes_low_n <- sum(yes_vals <= cutoff, na.rm = TRUE)
        yes_high_n <- sum(yes_vals > cutoff, na.rm = TRUE)
        no_low_n <- sum(no_vals <= cutoff, na.rm = TRUE)
        no_high_n <- sum(no_vals > cutoff, na.rm = TRUE)
        fisher_p <- safe_fisher(yes_low_n, yes_high_n, no_low_n, no_high_n)

        enrich_rows[[enrich_i]] <- data.frame(
          clade_group = clade_label,
          exposure_type = exposure_type,
          threshold_days = threshold_days,
          snv_rule = paste0("Min SNV <= ", cutoff),
          yes_low_n = yes_low_n,
          yes_high_n = yes_high_n,
          no_low_n = no_low_n,
          no_high_n = no_high_n,
          yes_low_pct = ifelse((yes_low_n + yes_high_n) > 0, 100 * yes_low_n / (yes_low_n + yes_high_n), NA_real_),
          no_low_pct = ifelse((no_low_n + no_high_n) > 0, 100 * no_low_n / (no_low_n + no_high_n), NA_real_),
          fisher_p_value = fisher_p,
          p_value_label = format_p(fisher_p),
          stringsAsFactors = FALSE
        )
        enrich_i <- enrich_i + 1
      }
    }
  }

  list(
    test_results = bind_rows(test_rows),
    enrichment_results = bind_rows(enrich_rows)
  )
}

make_median_plot <- function(test_results) {
  p <- ggplot(test_results, aes(x = threshold_days, y = yes_median_snv, group = exposure_type, color = exposure_type)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.8) +
    facet_wrap(~clade_group, ncol = 2, scales = "free_y") +
    scale_color_manual(values = location_colors) +
    labs(
      title = "Sequential yes median minimum SNV by threshold",
      subtitle = "Patient-pair analysis; smaller medians suggest a cleaner sequential signal",
      x = "Sequential threshold (days)",
      y = "Median minimum SNV among sequential yes pairs",
      color = "Location level"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(viz_dir, "sequential_threshold_median_snv_comparison.png"), p, width = 11, height = 6, units = "in", dpi = 300, bg = "white")
  save_plot_pdf(p, file.path(viz_dir, "sequential_threshold_median_snv_comparison.pdf"), 11, 6)
}

make_enrichment_plot <- function(enrichment_results) {
  plot_dat <- enrichment_results %>%
    filter(snv_rule == "Min SNV <= 2") %>%
    mutate(star_label = make_star_label(p_value_label))

  ann_dat <- plot_dat %>%
    group_by(clade_group, exposure_type, threshold_days) %>%
    summarise(
      y_pos = max(yes_low_pct, no_low_pct, na.rm = TRUE) * 1.15 + 0.2,
      star_label = first(star_label),
      .groups = "drop"
    )

  p <- bind_rows(
    plot_dat %>% transmute(clade_group, exposure_type, threshold_days, group = "Sequential yes", pct = yes_low_pct),
    plot_dat %>% transmute(clade_group, exposure_type, threshold_days, group = "Sequential no", pct = no_low_pct)
  ) %>%
    ggplot(aes(x = threshold_days, y = pct, fill = group)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    geom_text(
      aes(label = sprintf("%.1f", pct)),
      position = position_dodge(width = 0.75),
      vjust = -0.2,
      size = 3
    ) +
    geom_text(
      data = ann_dat,
      aes(x = threshold_days, y = y_pos, label = star_label),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 4.5
    ) +
    facet_grid(clade_group ~ exposure_type, scales = "free_y") +
    scale_fill_manual(values = c("Sequential yes" = "#F28E2B", "Sequential no" = "#4E79A7"), name = NULL) +
    labs(
      title = "Sequential threshold comparison for low minimum-SNV enrichment",
      subtitle = "Shown for minimum SNV <= 2; larger yes-vs-no separation suggests a clearer threshold",
      x = "Sequential threshold (days)",
      y = "Percent minimum SNV <= 2"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(viz_dir, "sequential_threshold_low_snv_enrichment.png"), p, width = 13, height = 8, units = "in", dpi = 300, bg = "white")
  save_plot_pdf(p, file.path(viz_dir, "sequential_threshold_low_snv_enrichment.pdf"), 13, 8)
}

clade1_gap_tables <- compute_gap_tables_for_clade("Clade 1")
clade2_gap_tables <- compute_gap_tables_for_clade("Clade 2")
cat("Building threshold comparison summaries...\n")

summary_list <- lapply(threshold_values, function(day_threshold) {
  clade1_flags <- annotate_threshold(min_snv_pairs, "Clade 1", day_threshold, clade1_gap_tables)
  clade2_flags <- annotate_threshold(min_snv_pairs, "Clade 2", day_threshold, clade2_gap_tables)

  threshold_dat <- bind_rows(
    min_snv_pairs %>%
      filter(clade_group == "Clade 1") %>%
      left_join(clade1_flags, by = c("patient_a", "patient_b")),
    min_snv_pairs %>%
      filter(clade_group == "Clade 2") %>%
      left_join(clade2_flags, by = c("patient_a", "patient_b"))
  ) %>%
    mutate(across(ends_with("_seq"), ~ replace_na(.x, FALSE)))

  summarize_threshold(threshold_dat, day_threshold)
})

test_results <- bind_rows(lapply(summary_list, `[[`, "test_results")) %>%
  mutate(
    exposure_type = factor(exposure_type, levels = location_levels),
    threshold_days = factor(threshold_days, levels = threshold_values)
  )
enrichment_results <- bind_rows(lapply(summary_list, `[[`, "enrichment_results")) %>%
  mutate(
    exposure_type = factor(exposure_type, levels = location_levels),
    threshold_days = factor(threshold_days, levels = threshold_values)
  )
cat("Writing outputs and plots...\n")

write_csv(test_results, file.path(desc_dir, "sequential_threshold_comparison_wilcox_tests.csv"))
write_csv(enrichment_results, file.path(desc_dir, "sequential_threshold_comparison_enrichment_tests.csv"))

writeLines(
  c(
    "This analysis compares sequential overlap thresholds of 7, 14, 30, and 60 days.",
    "It uses the minimum SNV per patient pair, rather than all isolate-pair rows, to make threshold choice more interpretable.",
    paste0("Special units excluded before analysis: ", paste(special_units, collapse = ", "), "."),
    "Sequential overlap is defined as same location, no direct timestamp overlap, and later stay begins within the threshold after the earlier stay ends.",
    "A more obvious threshold should show lower median minimum SNV among sequential yes pairs and stronger low-SNV enrichment relative to sequential no pairs."
  ),
  con = file.path(desc_dir, "sequential_threshold_comparison_notes.txt")
)

make_median_plot(test_results)
make_enrichment_plot(enrichment_results)

cat("Done. Sequential threshold comparison outputs written to:\n", desc_dir, "\n", sep = "")
