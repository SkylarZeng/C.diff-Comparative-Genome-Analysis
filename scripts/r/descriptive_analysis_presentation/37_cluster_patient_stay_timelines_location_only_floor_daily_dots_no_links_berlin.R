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
cluster_fig_dir <- file.path(output_dir, "cluster_patient_stay_timelines_location_only_floor_daily_dots_no_links_berlin")

clade1_clusters_path <- file.path(desc_dir, "clade_1_post_noncore_clusters_snv_lt2.csv")
clade2_clusters_path <- file.path(desc_dir, "clade_2_post_noncore_clusters_snv_lt2.csv")
locations_path <- file.path(desc_dir, "locations_filtered_case_unformed_clade1_clade2.tsv")

start_window <- as.Date("2016-02-02")
end_window <- as.Date("2017-12-31")

if (!dir.exists(cluster_fig_dir)) {
  dir.create(cluster_fig_dir, recursive = TRUE)
}

required_pkgs <- c("ggplot2", "tidyverse", "unikn")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(unikn)
})

save_plot_pdf <- function(plot_obj, filename, width, height) {
  ggsave(
    filename = filename,
    plot = plot_obj,
    device = "pdf",
    width = width,
    height = height,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
}

is_present <- function(x) {
  !is.na(x) & trimws(as.character(x)) != "" & trimws(as.character(x)) != "NA"
}

parse_date_safe <- function(x) {
  suppressWarnings(as.Date(substr(trimws(as.character(x)), 1, 10)))
}

make_floor_palette <- function(floor_values) {
  floor_values <- sort(unique(floor_values))
  if (length(floor_values) == 0) {
    return(character())
  }
  base_cols <- unikn::pal_unikn_light
  cols <- if (length(floor_values) <= length(base_cols)) {
    base_cols[seq_along(floor_values)]
  } else {
    grDevices::colorRampPalette(base_cols)(length(floor_values))
  }
  setNames(cols, floor_values)
}

find_collection_floor <- function(cluster_samples, stay_df) {
  if (nrow(cluster_samples) == 0 || nrow(stay_df) == 0) {
    return(cluster_samples %>% mutate(collection_floor = NA_character_))
  }

  map_dfr(seq_len(nrow(cluster_samples)), function(i) {
    row_i <- cluster_samples[i, , drop = FALSE]
    if (is.na(row_i$collection_date)) {
      matches <- stay_df[0, , drop = FALSE]
    } else {
      matches <- stay_df %>%
        filter(
          patient_id == row_i$patient_id,
          StartDate_date <= row_i$collection_date,
          EndDate_date >= row_i$collection_date
        ) %>%
        arrange(StartDate_date, EndDate_date)
    }

    row_i %>%
      mutate(collection_floor = ifelse(nrow(matches) > 0, matches$floor[1], NA_character_))
  })
}

expand_stay_days <- function(stay_df) {
  if (nrow(stay_df) == 0) {
    return(tibble())
  }

  stay_df %>%
    rowwise() %>%
    mutate(stay_day = list(seq(StartDate_date, EndDate_date, by = "day"))) %>%
    ungroup() %>%
    unnest(stay_day)
}

find_same_floor_same_day_links <- function(stay_df, patient_levels) {
  if (nrow(stay_df) < 2) {
    return(tibble(
      floor = character(),
      patient_a = character(),
      patient_b = character(),
      link_date = as.Date(character()),
      patient_plot_id_a = factor(levels = rev(patient_levels$patient_plot_id)),
      patient_plot_id_b = factor(levels = rev(patient_levels$patient_plot_id))
    ))
  }

  split(stay_df, stay_df$floor) %>%
    purrr::imap_dfr(function(floor_df, floor_name) {
      if (!is_present(floor_name) || nrow(floor_df) < 2) {
        return(tibble())
      }

      floor_df <- floor_df %>% arrange(StartDate_date, EndDate_date)
      links <- list()
      idx <- 1L

      for (i in seq_len(nrow(floor_df) - 1)) {
        for (j in seq((i + 1), nrow(floor_df))) {
          if (floor_df$patient_id[i] == floor_df$patient_id[j]) next

          overlap_start <- max(floor_df$StartDate_date[i], floor_df$StartDate_date[j])
          overlap_end <- min(floor_df$EndDate_date[i], floor_df$EndDate_date[j])

          if (!is.na(overlap_start) && !is.na(overlap_end) && overlap_start <= overlap_end) {
            midpoint_num <- as.numeric(overlap_start) + floor(as.numeric(overlap_end - overlap_start) / 2)
            pair_a <- min(floor_df$patient_id[i], floor_df$patient_id[j])
            pair_b <- max(floor_df$patient_id[i], floor_df$patient_id[j])

            links[[idx]] <- tibble(
              floor = floor_name,
              patient_a = pair_a,
              patient_b = pair_b,
              link_date = as.Date(midpoint_num, origin = "1970-01-01")
            )
            idx <- idx + 1L
          }
        }
      }

      if (length(links) == 0) {
        return(tibble())
      }

      bind_rows(links) %>%
        group_by(floor, patient_a, patient_b) %>%
        summarise(link_date = min(link_date), .groups = "drop") %>%
        mutate(
          patient_plot_id_a = factor(patient_levels$patient_plot_id[match(patient_a, patient_levels$patient_id)],
                                     levels = rev(patient_levels$patient_plot_id)),
          patient_plot_id_b = factor(patient_levels$patient_plot_id[match(patient_b, patient_levels$patient_id)],
                                     levels = rev(patient_levels$patient_plot_id))
        )
    })
}

clusters_all <- bind_rows(
  read.csv(clade1_clusters_path, check.names = FALSE, stringsAsFactors = FALSE),
  read.csv(clade2_clusters_path, check.names = FALSE, stringsAsFactors = FALSE)
) %>%
  mutate(
    genome_id = trimws(as.character(genome_id)),
    patient_id = trimws(as.character(patient_id)),
    collection_date = parse_date_safe(collection_date),
    ST = trimws(as.character(ST)),
    clade_group = trimws(as.character(clade_group)),
    cluster_id = as.integer(cluster_id)
  ) %>%
  filter(!is.na(genome_id), genome_id != "", !is.na(patient_id), patient_id != "")

clusters_keep <- clusters_all %>%
  group_by(clade_group, cluster_id) %>%
  summarise(
    genome_n = n(),
    patient_n = n_distinct(patient_id),
    .groups = "drop"
  ) %>%
  filter(genome_n > 2, patient_n > 1)

clusters_plot <- clusters_all %>%
  inner_join(clusters_keep, by = c("clade_group", "cluster_id"))

locations <- read.delim(locations_path, check.names = FALSE, stringsAsFactors = FALSE) %>%
  mutate(
    PatientID = trimws(as.character(PatientID)),
    clade_group = trimws(as.character(clade_group)),
    StartDate_date = parse_date_safe(StartDate),
    EndDate_date = coalesce(parse_date_safe(EndDate), parse_date_safe(StartDate), parse_date_safe(Date)),
    Unit = trimws(as.character(Unit)),
    floor = trimws(as.character(floor))
  ) %>%
  filter(
    !is.na(PatientID),
    PatientID != "",
    !is.na(StartDate_date),
    !is.na(EndDate_date),
    is_present(Unit),
    StartDate_date >= start_window,
    StartDate_date <= end_window
  )

global_floor_values <- unique(locations$floor)
global_floor_values <- global_floor_values[!is.na(global_floor_values) & global_floor_values != ""]
global_floor_palette <- make_floor_palette(global_floor_values)
global_floor_palette <- c(global_floor_palette, "Special Units" = "#4d4d4d")

patient_label_tables <- list()
link_summary_tables <- list()
cluster_keys <- clusters_keep %>%
  arrange(clade_group, cluster_id)

for (i in seq_len(nrow(cluster_keys))) {
  clade_label <- cluster_keys$clade_group[i]
  cluster_id_value <- cluster_keys$cluster_id[i]

  cluster_samples_all <- clusters_plot %>%
    filter(clade_group == clade_label, .data$cluster_id == cluster_id_value) %>%
    arrange(collection_date, patient_id, genome_id)

  stay_source <- locations %>%
    filter(clade_group == clade_label, PatientID %in% cluster_samples_all$patient_id) %>%
    transmute(
      patient_id = PatientID,
      StartDate_date,
      EndDate_date,
      floor
    )

  mapped_patients <- stay_source %>%
    distinct(patient_id)

  if (nrow(mapped_patients) == 0) {
    next
  }

  cluster_samples <- cluster_samples_all %>%
    semi_join(mapped_patients, by = "patient_id")

  if (nrow(cluster_samples) == 0) {
    next
  }

  patient_order <- cluster_samples %>%
    group_by(patient_id) %>%
    summarise(first_collection = min(collection_date, na.rm = TRUE), .groups = "drop") %>%
    arrange(first_collection, patient_id) %>%
    mutate(patient_plot_id = paste0("pt", row_number()))

  stay_df <- stay_source %>%
    semi_join(patient_order, by = "patient_id") %>%
    left_join(patient_order, by = "patient_id") %>%
    mutate(
      patient_plot_id = factor(patient_plot_id, levels = rev(patient_order$patient_plot_id)),
      floor = if_else(is.na(floor) | floor == "", "Special Units", floor)
    )

  stay_day_df <- expand_stay_days(stay_df) %>%
    mutate(patient_plot_id = factor(patient_plot_id, levels = rev(patient_order$patient_plot_id)))

  same_floor_links <- find_same_floor_same_day_links(stay_df, patient_order)

  cluster_samples_plot <- find_collection_floor(cluster_samples, stay_df) %>%
    left_join(patient_order, by = "patient_id") %>%
    mutate(
      patient_plot_id = factor(patient_plot_id, levels = rev(patient_order$patient_plot_id)),
      collection_floor = if_else(is.na(collection_floor) | collection_floor == "", "Special Units", collection_floor)
    )

  patient_label_tables[[length(patient_label_tables) + 1]] <- patient_order %>%
    mutate(clade_group = clade_label, cluster_id = cluster_id_value)

  if (nrow(same_floor_links) > 0) {
    link_summary_tables[[length(link_summary_tables) + 1]] <- same_floor_links %>%
      mutate(clade_group = clade_label, cluster_id = cluster_id_value)
  }

  genome_n <- nrow(cluster_samples)
  patient_n <- n_distinct(cluster_samples$patient_id)
  removed_patient_n <- n_distinct(cluster_samples_all$patient_id) - patient_n

  plot_title <- paste0(clade_label, " cluster ", cluster_id_value, ": patient stays and collection dates")
  plot_subtitle <- paste0(
    "n=", genome_n, " genomes across ", patient_n, " mapped patients",
    if (removed_patient_n > 0) paste0(" (", removed_patient_n, " patients without location data removed)") else "",
    if (nrow(same_floor_links) > 0) paste0("; ", nrow(same_floor_links), " same-floor same-day links shown") else ""
  )

  p <- ggplot()

  if (nrow(same_floor_links) > 0) {
    p <- p +
      geom_segment(
        data = same_floor_links,
        aes(x = link_date, xend = link_date, y = patient_plot_id_a, yend = patient_plot_id_b, color = floor),
        linetype = "dashed",
        linewidth = 0.8,
        alpha = 0.8
      )
  }

  p <- p +
    geom_point(
      data = stay_day_df,
      aes(x = stay_day, y = patient_plot_id, color = floor, shape = "Patient stay"),
      size = 3.4,
      stroke = 1.1,
      fill = "white",
      alpha = 0.6
    ) +
    geom_point(
      data = cluster_samples_plot,
      aes(x = collection_date, y = patient_plot_id, shape = "Sample collection"),
      size = 3.6,
      alpha = 0.95,
      color = "black"
    ) +
    scale_color_manual(values = global_floor_palette, name = "Floor") +
    scale_shape_manual(
      values = c("Patient stay" = 16, "Sample collection" = 17),
      name = NULL
    ) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Hospital stay timeline",
      y = "Patient"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 11),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(size = 4, alpha = 1)),
      shape = guide_legend(order = 2, override.aes = list(color = "#4d4d4d", fill = c("white", "#4d4d4d"), size = c(3.4, 4.2)))
    )

  clade_stub <- tolower(gsub(" ", "_", clade_label))
  file_stub <- paste0(clade_stub, "_cluster_", cluster_id_value, "_patient_stay_timeline_location_only_floor_daily_dots_no_links_berlin")
  plot_height <- max(4.5, patient_n * 0.7 + 2.5)

  ggsave(
    filename = file.path(cluster_fig_dir, paste0(file_stub, ".png")),
    plot = p,
    width = 12,
    height = plot_height,
    units = "in",
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
  save_plot_pdf(
    plot_obj = p,
    filename = file.path(cluster_fig_dir, paste0(file_stub, ".pdf")),
    width = 12,
    height = plot_height
  )
}

patient_label_lookup <- bind_rows(patient_label_tables) %>%
  select(clade_group, cluster_id, patient_plot_id, patient_id, first_collection)

write.csv(
  patient_label_lookup,
  file.path(desc_dir, "cluster_patient_timeline_patient_label_lookup_location_only_floor_daily_dots_no_links_berlin.csv"),
  row.names = FALSE
)

write.csv(
  bind_rows(link_summary_tables),
  file.path(desc_dir, "cluster_patient_same_floor_same_day_links_daily_dots_berlin.csv"),
  row.names = FALSE
)

writeLines(
  c(
    "Figures are generated for post-noncore SNV<2 clusters with more than 2 genomes and more than 1 patient.",
    "Only stay rows with StartDate between 2016-02-02 and 2017-12-31 are included.",
    "This location-only version removes patients without any location rows in the filtered location dataset within that window.",
    "Dashed vertical lines connect patients who share the same floor on overlapping calendar days.",
    "Point colors and the legend represent floor, using unikn::pal_unikn_light.",
    "Patient-stay dots are shown for every day during the stay interval, not only the start and end dates.",
    "Sample collection dates are shown as larger black triangles.",
    "Special Units represents records without a mapped floor value."
  ),
  con = file.path(desc_dir, "cluster_patient_stay_timelines_location_only_floor_daily_dots_no_links_berlin_notes.txt")
)

cat("Done. Cluster patient stay timeline figures with daily stay dots and no links written to:\n", cluster_fig_dir, "\n", sep = "")
