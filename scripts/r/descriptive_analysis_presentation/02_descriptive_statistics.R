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
metadata_path <- file.path(output_dir, "metadata_hospital_tcdA_tcdB_with_flags.csv")

if (!dir.exists(desc_dir)) {
  dir.create(desc_dir, recursive = TRUE)
}

if (!requireNamespace("flowchart", quietly = TRUE)) {
  stop("Package 'flowchart' is required. Install it with install.packages('flowchart').")
}

suppressPackageStartupMessages({
  library(flowchart)
  library(tidyverse)
})

metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE) %>%
  mutate(
    isolate_id = genome_id,
    patient_id = na_if(trimws(RDW_patient_id), ""),
    cdiff_case_num = suppressWarnings(as.numeric(as.character(cdiff_case))),
    stoolconsistency_num = suppressWarnings(as.numeric(as.character(stoolconsistency))),
    clade_group = if_else(is.na(clade) | trimws(as.character(clade)) == "", "Unknown", paste0("Clade ", as.character(clade)))
  )

count_entities <- function(df) {
  tibble::tibble(
    isolates = dplyr::n_distinct(df$isolate_id[!is.na(df$isolate_id) & df$isolate_id != ""]),
    patients = dplyr::n_distinct(df$patient_id[!is.na(df$patient_id) & df$patient_id != ""])
  )
}

make_box_text <- function(label, df) {
  counts <- count_entities(df)
  paste0(label, "\n(Isolates: ", counts$isolates, "; Patients: ", counts$patients, ")")
}

make_step_row <- function(step, rule, kept_df, excluded_df) {
  kept_counts <- count_entities(kept_df)
  excluded_counts <- count_entities(excluded_df)

  tibble::tibble(
    step = step,
    rule = rule,
    isolates_kept = kept_counts$isolates,
    patients_kept = kept_counts$patients,
    isolates_excluded = excluded_counts$isolates,
    patients_excluded = excluded_counts$patients
  )
}

relabel_last_filter <- function(fc_obj, kept_df, excluded_df, keep_label, excluded_label) {
  keep_row <- nrow(fc_obj$fc) - 1
  excluded_row <- nrow(fc_obj$fc)

  fc_obj$fc$label[keep_row] <- keep_label
  fc_obj$fc$text[keep_row] <- make_box_text(keep_label, kept_df)
  fc_obj$fc$label[excluded_row] <- excluded_label
  fc_obj$fc$text[excluded_row] <- make_box_text(excluded_label, excluded_df)
  fc_obj
}

relabel_last_split <- function(fc_obj, split_df, split_var) {
  split_values <- split_df[[split_var]]
  split_levels <- unique(split_values[!is.na(split_values) & split_values != ""])

  if (length(split_levels) == 0) {
    return(fc_obj)
  }

  split_rows <- seq.int(nrow(fc_obj$fc) - length(split_levels) + 1, nrow(fc_obj$fc))

  for (index in seq_along(split_levels)) {
    split_level <- split_levels[index]
    split_subset <- split_df[split_values == split_level, , drop = FALSE]
    fc_obj$fc$label[split_rows[index]] <- split_level
    fc_obj$fc$text[split_rows[index]] <- make_box_text(split_level, split_subset)
  }

  fc_obj
}

step0 <- metadata

step1 <- step0 %>%
  filter(!is.na(patient_id))
excluded1 <- anti_join(step0, step1, by = "isolate_id")

step2 <- step1 %>%
  filter(isolate == "a")
excluded2 <- anti_join(step1, step2, by = "isolate_id")

step3 <- step2 %>%
  filter(cdiff_case_num == 1)
excluded3 <- anti_join(step2, step3, by = "isolate_id")

step4 <- step3 %>%
  filter(stoolconsistency_num == 1)
excluded4 <- anti_join(step3, step4, by = "isolate_id")

step4_clade12 <- step4 %>%
  filter(clade_group %in% c("Clade 1", "Clade 2"))

clade_counts <- step4 %>%
  group_by(clade_group) %>%
  summarise(
    isolates = n_distinct(isolate_id),
    patients = n_distinct(patient_id[!is.na(patient_id) & patient_id != ""]),
    .groups = "drop"
  ) %>%
  arrange(clade_group)

flow_steps <- bind_rows(
  make_step_row("step0_input_metadata", "Starting merged metadata", step0, step0[0, ]),
  make_step_row("step1_keep_patient_id", "Exclude isolates with missing patient ID", step1, excluded1),
  make_step_row("step2_keep_primary_isolates", "Keep primary isolates only using isolate == 'a'; exclude other isolate labels", step2, excluded2),
  make_step_row("step3_keep_cases", "Keep cases only using cdiff_case == 1; exclude controls/non-cases", step3, excluded3),
  make_step_row("step4_keep_unformed", "Keep unformed stool only using stoolconsistency == 1; exclude formed/other", step4, excluded4)
)

fc_obj <- step0 %>%
  as_fc(label = "Merged metadata input")

fc_obj$fc$text[1] <- make_box_text("Merged metadata input", step0)
fc_obj$fc$text_fs[1] <- 6.5
fc_obj$fc$text_padding[1] <- 0.7

fc_obj <- fc_obj %>%
  fc_filter(!is.na(patient_id), label = "Keep non-missing patient ID", show_exc = TRUE, label_exc = "Excluded missing patient ID")
fc_obj <- relabel_last_filter(fc_obj, step1, excluded1, "Keep non-missing patient ID", "Excluded missing patient ID")

fc_obj <- fc_obj %>%
  fc_filter(isolate == "a", label = "Primary isolates only", show_exc = TRUE, label_exc = "Excluded non-primary isolates")
fc_obj <- relabel_last_filter(fc_obj, step2, excluded2, "Primary isolates only", "Excluded non-primary isolates")

fc_obj <- fc_obj %>%
  fc_filter(cdiff_case_num == 1, label = "Cases only", show_exc = TRUE, label_exc = "Excluded controls/non-cases")
fc_obj <- relabel_last_filter(fc_obj, step3, excluded3, "Cases only", "Excluded controls/non-cases")

fc_obj <- fc_obj %>%
  fc_filter(stoolconsistency_num == 1, label = "Unformed stool only", show_exc = TRUE, label_exc = "Excluded formed/other stool")
fc_obj <- relabel_last_filter(fc_obj, step4, excluded4, "Unformed stool only", "Excluded formed/other stool")

fc_obj <- fc_obj %>%
  fc_split(clade_group)
fc_obj <- relabel_last_split(fc_obj, step4, "clade_group")

fc_obj$fc$text_fs <- ifelse(fc_obj$fc$type == "excl", 5.2, 6.3)
fc_obj$fc$text_padding <- ifelse(fc_obj$fc$type == "excl", 0.55, 0.65)

write.csv(flow_steps, file.path(desc_dir, "flowchart_step_counts_case_unformed.csv"), row.names = FALSE)
write.csv(step4, file.path(desc_dir, "metadata_case_unformed.csv"), row.names = FALSE)
write.csv(step4_clade12, file.path(desc_dir, "metadata_case_unformed_clade1_clade2.csv"), row.names = FALSE)
write.csv(clade_counts, file.path(desc_dir, "clade_counts_case_unformed.csv"), row.names = FALSE)
write.csv(fc_obj$fc, file.path(desc_dir, "flowchart_box_layout_case_unformed.csv"), row.names = FALSE)

png(
  filename = file.path(desc_dir, "case_unformed_flowchart.png"),
  width = 3800,
  height = 1400,
  res = 220
)
fc_draw(fc_obj, title = "Case/Control and Stool Consistency Filtering")
dev.off()

pdf(
  file = file.path(desc_dir, "case_unformed_flowchart.pdf"),
  width = 18,
  height = 7.5
)
fc_draw(fc_obj, title = "Case/Control and Stool Consistency Filtering")
dev.off()

cat("Done. Descriptive statistics outputs written to:\n", desc_dir, "\n", sep = "")