library(tidyverse)
library(readxl)
library(here)

# Step 4: join patient, toxin, and MLST metadata for descriptive context.

ids <- read_csv(here("data", "metadata", "cleaned_passed_genomeid.csv"), col_names = FALSE)$X1

patient_data <- read_csv(here("data", "metadata", "hospital_onset_data.csv")) %>%
  filter(genome_id %in% ids)

patient_data_tx <- read_excel(here("data", "metadata", "summary_toxin_metadata.xlsx"))

mlst_data <- read_tsv(here("data", "metadata", "mlst_replaced.tsv")) %>%
  mutate(ST = str_remove(ST, "\\*"))

patient_and_mlst <- inner_join(patient_data, mlst_data, by = "genome_id")

st_summary <- patient_and_mlst %>%
  count(ST, sort = TRUE)

onset_summary <- patient_data_tx %>%
  count(hospital_onset, sort = TRUE)

print(head(st_summary, 10))
print(onset_summary)
