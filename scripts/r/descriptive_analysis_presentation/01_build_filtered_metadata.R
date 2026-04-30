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
metadata_dir <- file.path(cluster_dir, "Metadata")
output_dir <- file.path(base_dir, "outputs")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

hospital_path <- file.path(metadata_dir, "hospital_onset_data.csv")
tcdA_path <- file.path(metadata_dir, "tcdA_distinct.csv")
tcdB_path <- file.path(metadata_dir, "tcdB_distinct.csv")
mlst_path <- file.path(metadata_dir, "mlst_replaced.tsv")
wiens_path <- file.path(metadata_dir, "main_from_wiens_group_filtered_cdiff_cases.tsv")

hospital <- read.csv(hospital_path, check.names = FALSE, stringsAsFactors = FALSE)
tcdA <- read.csv(tcdA_path, check.names = FALSE, stringsAsFactors = FALSE)
tcdB <- read.csv(tcdB_path, check.names = FALSE, stringsAsFactors = FALSE)
mlst <- read.delim(mlst_path, check.names = FALSE, stringsAsFactors = FALSE)
wiens <- read.delim(wiens_path, check.names = FALSE, stringsAsFactors = FALSE)

if ("ST" %in% colnames(mlst)) {
  mlst$ST <- gsub("\\*", "", mlst$ST)
}

parse_date_safely <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr[x_chr == ""] <- NA_character_
  as.Date(sub("T.*$", "", x_chr))
}

parse_numeric_safely <- function(x) {
  suppressWarnings(as.numeric(trimws(as.character(x))))
}

classify_ho_co <- function(days_from_admit) {
  if (is.na(days_from_admit)) {
    return("unknown/NA")
  }
  if (days_from_admit > 3) {
    return("HO")
  }
  if (days_from_admit >= 0) {
    return("CO")
  }
  "unknown/NA"
}

derive_hospital_onset_updated <- function(df) {
  patient_id <- trimws(as.character(df$RDW_patient_id))
  patient_id[patient_id == ""] <- NA_character_

  onset_date <- parse_date_safely(df$CollectionDateCdiff)
  collection_date <- parse_date_safely(df$collection_date)
  onset_date[is.na(onset_date)] <- collection_date[is.na(onset_date)]

  days_from_admit <- parse_numeric_safely(df$days_btwn_admit_and_onset)
  cdiff14 <- parse_numeric_safely(df$Cdiff14)
  cdiff90 <- parse_numeric_safely(df$Cdiff90)
  cdiff1yr <- parse_numeric_safely(df$Cdiff1yr)

  classifications <- rep("unknown/NA", nrow(df))

  valid_patients <- unique(patient_id[!is.na(patient_id)])

  for (pid in valid_patients) {
    patient_rows <- which(patient_id == pid)
    patient_onset <- onset_date[patient_rows]

    dated_rows <- patient_rows[!is.na(patient_onset)]
    undated_rows <- setdiff(patient_rows, dated_rows)

    if (length(dated_rows) > 0) {
        raw_dates <- as.Date(sort(unique(as.character(onset_date[dated_rows]))))
        grouped_episode_dates <- as.Date(character())

        if (length(raw_dates) > 0) {
          current_anchor <- raw_dates[[1]]
          grouped_episode_dates <- c(grouped_episode_dates, current_anchor)

          if (length(raw_dates) > 1) {
            for (date_idx in 2:length(raw_dates)) {
              candidate_date <- raw_dates[[date_idx]]
              if (as.numeric(candidate_date - current_anchor) > 14) {
                grouped_episode_dates <- c(grouped_episode_dates, candidate_date)
                current_anchor <- candidate_date
              }
            }
          }
        }

        previous_episode_date <- as.Date(NA)

        for (episode_date in grouped_episode_dates) {
          episode_group_rows <- patient_rows[
            !is.na(onset_date[patient_rows]) &
              as.numeric(onset_date[patient_rows] - episode_date) >= 0 &
              as.numeric(onset_date[patient_rows] - episode_date) <= 14
          ]

          if (length(grouped_episode_dates) > 1) {
            next_episode_idx <- match(episode_date, grouped_episode_dates) + 1
            if (!is.na(next_episode_idx) && next_episode_idx <= length(grouped_episode_dates)) {
              next_episode_date <- grouped_episode_dates[[next_episode_idx]]
              episode_group_rows <- episode_group_rows[onset_date[episode_group_rows] < next_episode_date]
            }
          }

          episode_rows <- episode_group_rows

          episode_has_cdiff14 <- any(cdiff14[episode_rows] == 1, na.rm = TRUE)
          episode_has_cdiff90 <- any(cdiff90[episode_rows] == 1, na.rm = TRUE)
          episode_has_cdiff1yr <- any(cdiff1yr[episode_rows] == 1, na.rm = TRUE)
        episode_days <- days_from_admit[episode_rows]
        episode_days <- episode_days[!is.na(episode_days)]
        episode_days_value <- if (length(episode_days) > 0) episode_days[[1]] else NA_real_

        if (!is.na(previous_episode_date) &&
            as.numeric(episode_date - previous_episode_date) <= 56) {
          episode_class <- "rCDI"
        } else if (is.na(previous_episode_date) && episode_has_cdiff14) {
          episode_class <- "rCDI"
        } else if (is.na(previous_episode_date) && episode_has_cdiff90 && !episode_has_cdiff14) {
          episode_class <- "unknown/NA"
        } else if (is.na(previous_episode_date) && episode_has_cdiff1yr && !episode_has_cdiff90 && !episode_has_cdiff14) {
          episode_class <- classify_ho_co(episode_days_value)
        } else {
          episode_class <- classify_ho_co(episode_days_value)
        }

        classifications[episode_rows] <- episode_class
        previous_episode_date <- episode_date
      }
    }

    if (length(undated_rows) > 0) {
      for (row_idx in undated_rows) {
        if (!is.na(cdiff14[row_idx]) && cdiff14[row_idx] == 1) {
          classifications[row_idx] <- "rCDI"
        } else if (!is.na(cdiff90[row_idx]) && cdiff90[row_idx] == 1) {
          classifications[row_idx] <- "unknown/NA"
        } else if (!is.na(cdiff1yr[row_idx]) && cdiff1yr[row_idx] == 1) {
          classifications[row_idx] <- classify_ho_co(days_from_admit[row_idx])
        } else {
          classifications[row_idx] <- classify_ho_co(days_from_admit[row_idx])
        }
      }
    }
  }

  classifications[is.na(patient_id)] <- vapply(
    which(is.na(patient_id)),
    function(i) {
      if (!is.na(cdiff14[i]) && cdiff14[i] == 1) {
        "rCDI"
      } else if (!is.na(cdiff90[i]) && cdiff90[i] == 1) {
        "unknown/NA"
      } else if (!is.na(cdiff1yr[i]) && cdiff1yr[i] == 1) {
        classify_ho_co(days_from_admit[i])
      } else {
        classify_ho_co(days_from_admit[i])
      }
    },
    character(1)
  )

  classifications
}

count_entities <- function(df) {
  isolate_count <- length(unique(df$sample_id[!is.na(df$sample_id) & df$sample_id != ""]))
  patient_count <- length(unique(df$subject_id[!is.na(df$subject_id) & df$subject_id != ""]))
  c(isolates = isolate_count, patients = patient_count)
}

record_step <- function(step_name, rule, df_current, df_previous = NULL) {
  current <- count_entities(df_current)
  if (is.null(df_previous)) {
    excluded <- c(isolates = 0, patients = 0)
  } else {
    prev <- count_entities(df_previous)
    excluded <- prev - current
  }

  data.frame(
    step = step_name,
    rule = rule,
    isolates_kept = unname(current["isolates"]),
    isolates_excluded = unname(excluded["isolates"]),
    patients_kept = unname(current["patients"]),
    patients_excluded = unname(excluded["patients"]),
    stringsAsFactors = FALSE
  )
}

step0 <- hospital
step1 <- step0[!is.na(step0$genome_id) & step0$genome_id != "", ]

tcdA_unique <- tcdA[!duplicated(tcdA$genomeid), ]
tcdB_unique <- tcdB[!duplicated(tcdB$genomeid), ]

tcdA_ids <- unique(tcdA_unique$genomeid[!is.na(tcdA_unique$genomeid) & tcdA_unique$genomeid != ""])
tcdB_ids <- unique(tcdB_unique$genomeid[!is.na(tcdB_unique$genomeid) & tcdB_unique$genomeid != ""])

tcdA_cols <- setdiff(colnames(tcdA_unique), "genomeid")
colnames(tcdA_unique)[match(tcdA_cols, colnames(tcdA_unique))] <- paste0("tcdA_", tcdA_cols)

tcdB_cols <- setdiff(colnames(tcdB_unique), "genomeid")
colnames(tcdB_unique)[match(tcdB_cols, colnames(tcdB_unique))] <- paste0("tcdB_", tcdB_cols)

joined_tcdA <- merge(
  step1,
  tcdA_unique,
  by.x = "genome_id",
  by.y = "genomeid",
  all.x = TRUE,
  sort = FALSE
)
step2 <- joined_tcdA

joined_tcdAB <- merge(
  step2,
  tcdB_unique,
  by.x = "genome_id",
  by.y = "genomeid",
  all.x = TRUE,
  sort = FALSE
)
step3 <- joined_tcdAB
step3$has_tcdA <- step3$genome_id %in% tcdA_ids
step3$has_tcdB <- step3$genome_id %in% tcdB_ids
step3$has_any_toxin <- step3$has_tcdA | step3$has_tcdB

step3$toxin_status <- ifelse(
  step3$has_tcdA & step3$has_tcdB,
  "Both",
  ifelse(
    step3$has_tcdA & !step3$has_tcdB,
    "Only tcdA",
    ifelse(!step3$has_tcdA & step3$has_tcdB, "Only tcdB", "None")
  )
)

cdiff_case_chr <- as.character(step3$cdiff_case)
stool_chr <- as.character(step3$stoolconsistency)
is_formed <- grepl("formed", step3$genome_id, ignore.case = TRUE)
is_cdif <- grepl("CDIF", step3$genome_id, ignore.case = TRUE)

step3$sample_type <- ifelse(
  is_formed & (stool_chr == "0" | is.na(step3$stoolconsistency)),
  "control_formed",
  ifelse(
    is_cdif & cdiff_case_chr == "0" & stool_chr == "0",
    "control_formed",
    ifelse(
      is_cdif & cdiff_case_chr == "0" & stool_chr == "1",
      "control_unformed",
      ifelse(
        is_cdif & cdiff_case_chr == "1" & stool_chr == "1",
        "case_unformed",
        ifelse(
          is_cdif & cdiff_case_chr == "1" & stool_chr == "0",
          "case_formed",
          ifelse(is_cdif & (is.na(step3$cdiff_case) | stool_chr == "2"), "NA", NA_character_)
        )
      )
    )
  )
)

# Keep only metadata columns plus toxin flag variables in final output.
drop_cols <- grepl("^tcdA_|^tcdB_", colnames(step3))
step3_metadata_only <- step3[, !drop_cols, drop = FALSE]

mlst_keep_cols <- intersect(c("genome_id", "ST", "clade"), colnames(mlst))
mlst_small <- mlst[, mlst_keep_cols, drop = FALSE]
if ("genome_id" %in% colnames(mlst_small)) {
  mlst_small <- mlst_small[!duplicated(mlst_small$genome_id), , drop = FALSE]
} else {
  mlst_small <- unique(mlst_small)
}
step3_metadata_only <- merge(
  step3_metadata_only,
  mlst_small,
  by = "genome_id",
  all.x = TRUE,
  sort = FALSE
)

if ("PatientID" %in% colnames(wiens)) {
  wiens_unique <- wiens[!duplicated(wiens$PatientID), , drop = FALSE]
  step3_metadata_only <- merge(
    step3_metadata_only,
    wiens_unique,
    by.x = "RDW_patient_id",
    by.y = "PatientID",
    all.x = TRUE,
    sort = FALSE
  )
}

step3_metadata_only$hospital_onset_updated <- derive_hospital_onset_updated(step3_metadata_only)

steps_table <- do.call(
  rbind,
  list(
    record_step("step0_hospital_onset_data", "Starting dataset", step0, NULL),
    record_step("step1_after_genome_id_present", "Keep rows with non-missing genome_id in hospital_onset_data", step1, step0),
    record_step("step2_joined_tcdA", "Left-join tcdA_distinct by genome_id", step2, step1),
    record_step("step3_joined_tcdB_add_flags", "Left-join tcdB_distinct and add has_tcdA/has_tcdB", step3, step2)
  )
)

write.csv(step1, file.path(output_dir, "step1_hospital_with_genome_id.csv"), row.names = FALSE)
write.csv(joined_tcdA, file.path(output_dir, "step2_joined_with_tcdA_all_rows.csv"), row.names = FALSE)
write.csv(joined_tcdAB, file.path(output_dir, "step3_joined_with_tcdA_tcdB_all_rows.csv"), row.names = FALSE)
write.csv(step3_metadata_only, file.path(output_dir, "metadata_hospital_tcdA_tcdB_with_flags.csv"), row.names = FALSE)
write.csv(steps_table, file.path(output_dir, "filtering_step_counts.csv"), row.names = FALSE)

summary_lines <- c(
  "# Metadata Filtering Summary",
  "",
  "This file is auto-generated by scripts/01_build_filtered_metadata.R.",
  "",
  sprintf("- Start (hospital_onset_data): %d isolates, %d patients", steps_table$isolates_kept[1], steps_table$patients_kept[1]),
  sprintf("- After requiring genome_id: %d isolates, %d patients (excluded: %d isolates, %d patients)",
    steps_table$isolates_kept[2], steps_table$patients_kept[2],
    steps_table$isolates_excluded[2], steps_table$patients_excluded[2]),
  sprintf("- After tcdA left join: %d isolates, %d patients (excluded: %d isolates, %d patients)",
    steps_table$isolates_kept[3], steps_table$patients_kept[3],
    steps_table$isolates_excluded[3], steps_table$patients_excluded[3]),
  sprintf("- Final after tcdB left join + flags: %d isolates, %d patients (excluded: %d isolates, %d patients)",
    steps_table$isolates_kept[4], steps_table$patients_kept[4],
    steps_table$isolates_excluded[4], steps_table$patients_excluded[4]),
  "",
  sprintf("- has_tcdA = TRUE: %d isolates", sum(step3_metadata_only$has_tcdA, na.rm = TRUE)),
  sprintf("- has_tcdB = TRUE: %d isolates", sum(step3_metadata_only$has_tcdB, na.rm = TRUE)),
  "",
  "- hospital_onset_updated categories:",
  paste(sprintf("  - %s: %d", names(table(step3_metadata_only$hospital_onset_updated, useNA = 'ifany')), as.integer(table(step3_metadata_only$hospital_onset_updated, useNA = 'ifany'))), collapse = "\n")
)
writeLines(summary_lines, file.path(output_dir, "filtering_step_summary.md"))

flow_lines <- c(
  "flowchart TD",
  sprintf("A[Hospital onset input\\nIsolates: %d\\nPatients: %d]", steps_table$isolates_kept[1], steps_table$patients_kept[1]),
  sprintf("B[After requiring genome_id\\nIsolates kept: %d\\nPatients kept: %d]", steps_table$isolates_kept[2], steps_table$patients_kept[2]),
  sprintf("C[After tcdA left join\\nIsolates kept: %d\\nPatients kept: %d]", steps_table$isolates_kept[3], steps_table$patients_kept[3]),
  sprintf("D[After tcdB left join + flags\\nIsolates kept: %d\\nPatients kept: %d]", steps_table$isolates_kept[4], steps_table$patients_kept[4]),
  sprintf("X1[Excluded for missing genome_id\\nIsolates: %d\\nPatients: %d]", steps_table$isolates_excluded[2], steps_table$patients_excluded[2]),
  sprintf("X2[Excluded at tcdA step\\nIsolates: %d\\nPatients: %d]", steps_table$isolates_excluded[3], steps_table$patients_excluded[3]),
  sprintf("X3[Excluded at tcdB step\\nIsolates: %d\\nPatients: %d]", steps_table$isolates_excluded[4], steps_table$patients_excluded[4]),
  sprintf("T1[has_tcdA = TRUE\\nIsolates: %d]", sum(step3_metadata_only$has_tcdA, na.rm = TRUE)),
  sprintf("T2[has_tcdB = TRUE\\nIsolates: %d]", sum(step3_metadata_only$has_tcdB, na.rm = TRUE)),
  "A --> B",
  "B --> C",
  "C --> D",
  "B -.-> X1",
  "C -.-> X2",
  "D -.-> X3",
  "D --> T1",
  "D --> T2"
)

writeLines(flow_lines, file.path(output_dir, "filtering_flowchart.mmd"))

cat("Done. Outputs written to:\n", output_dir, "\n", sep = "")
