library(data.table)
library(tidyverse)
library(here)

# Step 8: connect cluster assignments to patient-location overlap.

patient_data <- read_csv(here("data", "metadata", "hospital_onset_data.csv"), show_col_types = FALSE)
mlst_data <- read_tsv(here("data", "metadata", "mlst_replaced.tsv"), show_col_types = FALSE) %>%
  mutate(ST = stringr::str_remove(ST, "\\*"))
locations <- read_tsv(here("data", "metadata", "locations_from_wiens_group_filtered.tsv"), show_col_types = FALSE)
pre_clusters <- fread(here("2026_analysis", "clade2_pre_core_clusters.csv"))

patient_mlst_cluster <- inner_join(patient_data, mlst_data, by = "genome_id") %>%
  left_join(pre_clusters, by = c("genome_id" = "genome")) %>%
  filter(clade == "2")

locations_with_cluster <- locations %>%
  left_join(patient_mlst_cluster, by = c("PatientID" = "RDW_patient_id")) %>%
  filter(!is.na(cluster_pre), EndDate > StartDate)

dt <- as.data.table(
  locations_with_cluster %>%
    distinct(PatientID, FacilityCode, Unit, StartDate, EndDate, genome_id, cluster_pre)
)

setorder(dt, FacilityCode, Unit, StartDate)

sequential_unit <- dt[, .SD[order(StartDate)], by = .(FacilityCode, Unit)][
  , `:=`(
    next_PatientID = shift(PatientID, type = "lead"),
    next_StartDate = shift(StartDate, type = "lead"),
    next_cluster = shift(cluster_pre, type = "lead"),
    next_genome = shift(genome_id, type = "lead")
  ),
  by = .(FacilityCode, Unit)
][
  !is.na(next_PatientID) & PatientID != next_PatientID & EndDate <= next_StartDate
][
  , time_between_days := as.integer(difftime(next_StartDate, EndDate, units = "days"))
]

same_cluster_unit <- sequential_unit[cluster_pre == next_cluster]

print(same_cluster_unit[, .(
  n_pairs = .N,
  min_days = min(time_between_days),
  median_days = median(time_between_days),
  max_days = max(time_between_days)
)])
