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

location_path <- file.path(cluster_dir, "Metadata", "locations_from_wiens_group_filtered.tsv")
unmapped_location_path <- file.path(cluster_dir, "Metadata", "unmapped_samples_hosp_onset_noPHI.csv")
metadata_path <- file.path(desc_dir, "metadata_case_unformed_clade1_clade2.csv")

if (!dir.exists(desc_dir)) {
  dir.create(desc_dir, recursive = TRUE)
}
if (!dir.exists(viz_dir)) {
  dir.create(viz_dir, recursive = TRUE)
}

required_pkgs <- c("tidyverse", "unikn")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
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
    bg = "white"
  )
}

is_present <- function(x) {
  !is.na(x) & trimws(as.character(x)) != "" & trimws(as.character(x)) != "NA"
}

derive_floor_label <- function(facility_code, unit_value) {
  facility_code <- trimws(as.character(facility_code))
  unit_value <- trimws(as.character(unit_value))

  has_floor <- is_present(facility_code) & is_present(unit_value) & str_detect(unit_value, "^\\d+")
  floor_num <- str_extract(unit_value, "^\\d+")

  if_else(has_floor, paste0(facility_code, "-", floor_num), NA_character_)
}

locations <- read.delim(location_path, check.names = FALSE, stringsAsFactors = FALSE)
unmapped_locations <- read.csv(unmapped_location_path, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
names(metadata) <- trimws(names(metadata))

patient_clade_lookup <- metadata %>%
  transmute(
    patient_id = trimws(as.character(patient_id)),
    sample_id = trimws(as.character(sample_id)),
    clade_group = trimws(as.character(clade_group))
  ) %>%
  filter(
    !is.na(patient_id),
    patient_id != "",
    clade_group %in% c("Clade 1", "Clade 2")
  ) %>%
  distinct(patient_id, sample_id, clade_group)

unmapped_locations_formatted <- unmapped_locations %>%
  mutate(
    sample_id = trimws(as.character(sample_id)),
    collection_date = as.Date(collection_date),
    FacilityCode = recode(trimws(as.character(facility)), "MH" = "UMH", "UH" = "UUH", .default = trimws(as.character(facility))),
    Unit = trimws(as.character(unit)),
    Room = trimws(as.character(room_num)),
    Bed = trimws(as.character(bed_num))
  ) %>%
  inner_join(
    patient_clade_lookup %>% select(patient_id, sample_id, clade_group),
    by = "sample_id",
    relationship = "many-to-many"
  ) %>%
  transmute(
    PatientID = patient_id,
    Date = as.character(collection_date),
    StartDate = as.character(collection_date),
    EndDate = as.character(collection_date),
    Bed = Bed,
    FacilityCode = FacilityCode,
    Room = Room,
    Unit = Unit
  )

locations_combined <- bind_rows(
  locations,
  unmapped_locations_formatted
) %>%
  distinct()

locations_filtered <- locations_combined %>%
  mutate(
    PatientID = trimws(as.character(PatientID)),
    FacilityCode = trimws(as.character(FacilityCode)),
    Unit = trimws(as.character(Unit)),
    Room = trimws(as.character(Room)),
    Bed = trimws(as.character(Bed))
  ) %>%
  inner_join(patient_clade_lookup, by = c("PatientID" = "patient_id"), relationship = "many-to-many") %>%
  filter(FacilityCode != "UBW01") %>%
  mutate(floor = derive_floor_label(FacilityCode, Unit)) %>%
  relocate(clade_group, .after = PatientID)

availability_summary <- locations_filtered %>%
  group_by(clade_group) %>%
  summarise(
    location_rows = n(),
    unique_patients = n_distinct(PatientID),
    facility_rows_available = sum(is_present(FacilityCode)),
    facility_unique_types = n_distinct(FacilityCode[is_present(FacilityCode)]),
    room_rows_available = sum(is_present(Room)),
    room_unique_types = n_distinct(Room[is_present(Room)]),
    unit_rows_available = sum(is_present(Unit)),
    unit_unique_types = n_distinct(Unit[is_present(Unit)]),
    bed_rows_available = sum(is_present(Bed)),
    bed_unique_types = n_distinct(Bed[is_present(Bed)]),
    .groups = "drop"
  )

availability_long <- bind_rows(
  locations_filtered %>%
    group_by(clade_group) %>%
    summarise(location_type = "Facility", rows_available = sum(is_present(FacilityCode)), unique_types = n_distinct(FacilityCode[is_present(FacilityCode)]), .groups = "drop"),
  locations_filtered %>%
    group_by(clade_group) %>%
    summarise(location_type = "Floor", rows_available = sum(is_present(floor)), unique_types = n_distinct(floor[is_present(floor)]), .groups = "drop"),
  locations_filtered %>%
    group_by(clade_group) %>%
    summarise(location_type = "Room", rows_available = sum(is_present(Room)), unique_types = n_distinct(Room[is_present(Room)]), .groups = "drop"),
  locations_filtered %>%
    group_by(clade_group) %>%
    summarise(location_type = "Unit", rows_available = sum(is_present(Unit)), unique_types = n_distinct(Unit[is_present(Unit)]), .groups = "drop"),
  locations_filtered %>%
    group_by(clade_group) %>%
    summarise(location_type = "Bed", rows_available = sum(is_present(Bed)), unique_types = n_distinct(Bed[is_present(Bed)]), .groups = "drop")
) %>%
  arrange(clade_group, match(location_type, c("Facility", "Floor", "Room", "Unit", "Bed")))

facility_keep <- locations_filtered %>%
  filter(is_present(FacilityCode)) %>%
  group_by(clade_group, FacilityCode) %>%
  summarise(patient_n = n_distinct(PatientID), .groups = "drop") %>%
  filter(patient_n > 1) %>%
  select(clade_group, FacilityCode)

locations_plot_filtered <- locations_filtered %>%
  inner_join(facility_keep, by = c("clade_group", "FacilityCode"))

location_patient_summary <- bind_rows(
  locations_plot_filtered %>%
    filter(is_present(FacilityCode)) %>%
    group_by(clade_group, location_type = "Facility", location_value = FacilityCode) %>%
    summarise(patient_n = n_distinct(PatientID), row_n = n(), .groups = "drop"),
  locations_plot_filtered %>%
    filter(is_present(floor)) %>%
    group_by(clade_group, location_type = "Floor", location_value = floor) %>%
    summarise(patient_n = n_distinct(PatientID), row_n = n(), .groups = "drop"),
  locations_plot_filtered %>%
    filter(is_present(Unit)) %>%
    group_by(clade_group, location_type = "Unit", location_value = Unit) %>%
    summarise(patient_n = n_distinct(PatientID), row_n = n(), .groups = "drop")
) %>%
  arrange(clade_group, match(location_type, c("Facility", "Floor", "Unit")), desc(patient_n), location_value)

make_location_plot <- function(df, clade_label, top_n = 15) {
  plot_df <- df %>%
    filter(clade_group == clade_label) %>%
    mutate(location_type = factor(location_type, levels = c("Facility", "Floor", "Unit"))) %>%
    group_by(location_type) %>%
    slice_max(order_by = patient_n, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      panel_label = location_value,
      bar_group = case_when(
        location_type == "Unit" & location_value == "ESA" ~ "ESA",
        location_type == "Unit" & location_value == "OR" ~ "OR",
        location_type == "Unit" & location_value == "DVU" ~ "DVU",
        location_type == "Unit" & location_value == "MPUA" ~ "MPUA",
        location_type == "Unit" & location_value == "UADM" ~ "UADM",
        location_type == "Unit" & location_value == "USSS" ~ "USSS",
        TRUE ~ "Other"
      )
    )

  ggplot(
    plot_df,
    aes(x = patient_n, y = reorder(panel_label, patient_n), fill = bar_group)
  ) +
    geom_col(width = 0.75) +
    geom_text(aes(label = patient_n), hjust = -0.1, size = 5.8) +
    facet_wrap(~location_type, scales = "free_y", ncol = 3) +
    scale_fill_manual(
      values = c(
        "ESA" = unname(unikn::pal_unikn_pref[1]),
        "OR" = unname(unikn::pal_unikn_pref[2]),
        "DVU" = unname(unikn::pal_unikn_pref[3]),
        "MPUA" = unname(unikn::pal_unikn_pref[4]),
        "UADM" = unname(unikn::pal_unikn_pref[5]),
        "USSS" = unname(unikn::pal_unikn_pref[6]),
        "Other" = "#1F77B4"
      ),
      guide = "none"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = paste0(clade_label, ": patient counts by location"),
      subtitle = "Facility, derived floor, and unit counts after removing facilities with only 1 patient",
      x = "Number of patients",
      y = NULL
    ) +
    theme_minimal(base_size = 20) +
    theme(
      plot.title = element_text(face = "bold", size = 24),
      plot.subtitle = element_text(size = 18),
      strip.text = element_text(face = "bold", size = 19),
      axis.text = element_text(size = 17),
      axis.title.x = element_text(size = 18),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}

write.table(
  locations_filtered,
  file.path(desc_dir, "locations_filtered_case_unformed_clade1_clade2.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  availability_summary,
  file.path(desc_dir, "locations_availability_summary_clade1_clade2.csv"),
  row.names = FALSE
)

write.csv(
  availability_long,
  file.path(desc_dir, "locations_availability_by_type_clade1_clade2.csv"),
  row.names = FALSE
)

write.csv(
  location_patient_summary,
  file.path(desc_dir, "locations_patient_counts_by_value_clade1_clade2.csv"),
  row.names = FALSE
)

clade1_location_plot <- make_location_plot(location_patient_summary, "Clade 1")
clade2_location_plot <- make_location_plot(location_patient_summary, "Clade 2")

ggsave(
  filename = file.path(viz_dir, "clade1_location_patient_counts.png"),
  plot = clade1_location_plot,
  width = 16,
  height = 12,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = clade1_location_plot,
  filename = file.path(viz_dir, "clade1_location_patient_counts.pdf"),
  width = 16,
  height = 12
)

ggsave(
  filename = file.path(viz_dir, "clade2_location_patient_counts.png"),
  plot = clade2_location_plot,
  width = 16,
  height = 12,
  units = "in",
  dpi = 300,
  bg = "white"
)
save_plot_pdf(
  plot_obj = clade2_location_plot,
  filename = file.path(viz_dir, "clade2_location_patient_counts.pdf"),
  width = 16,
  height = 12
)

writeLines(
  c(
    "Filtered location file includes only patients present in metadata_case_unformed_clade1_clade2.csv.",
    "Additional one-day location rows are appended from unmapped_samples_hosp_onset_noPHI.csv by sample_id when the sample matches the cleaned clade 1/2 metadata.",
    "For appended rows from unmapped_samples_hosp_onset_noPHI.csv, facility codes are harmonized as MH -> UMH and UH -> UUH.",
    "FacilityCode UBW01 is excluded from the filtered location dataset, summaries, and figures.",
    "If a patient appears in both clades, that patient's location rows appear once for each clade assignment.",
    "rows_available counts non-missing entries for each location field.",
    "unique_types counts distinct non-missing values observed for each location field.",
    "floor is derived when Unit begins with a number and is labeled as FacilityCode-floor, for example UUH-8.",
    "locations_patient_counts_by_value_clade1_clade2.csv reports, for each clade and each Facility/Floor/Unit value used in the figure, how many distinct patients appear in that location.",
    "The clade1_location_patient_counts and clade2_location_patient_counts figures show the top location values within Facility, Floor, and Unit after removing facilities with only 1 patient."
  ),
  con = file.path(desc_dir, "locations_availability_clade1_clade2_notes.txt")
)

cat("Done. Filtered location outputs written to:\n", desc_dir, "\n", sep = "")
