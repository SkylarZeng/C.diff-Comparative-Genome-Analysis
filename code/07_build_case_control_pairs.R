library(tidyverse)
library(readr)
library(here)

# Step 7: convert SNV distances into within-ST case/control style pair tables.

snv_case_threshold <- 2
snv_control_threshold <- 10

pre_core <- read_csv(here("data", "SNV_matrices", "clade2_pre_core.csv"), show_col_types = FALSE) %>%
  rename(snv_pre = snvs)

post_noncore <- read_csv(here("data", "SNV_matrices", "clade2_post_noncore_snv_long.csv"), show_col_types = FALSE) %>%
  rename(snv_post_noncore = snvs)

post_softcore <- read_csv(here("data", "SNV_matrices", "clade2_post_softcore_snv_long.csv"), show_col_types = FALSE) %>%
  rename(snv_post_softcore = snvs)

mlst_all <- read_tsv(here("data", "metadata", "mlst_replaced.tsv"), show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  select(genome_id, st, clade) %>%
  filter(clade == 2, !is.na(st))

snv_with_st <- pre_core %>%
  full_join(post_noncore, by = c("id1", "id2")) %>%
  full_join(post_softcore, by = c("id1", "id2")) %>%
  left_join(mlst_all %>% rename(st1 = st), by = c("id1" = "genome_id")) %>%
  left_join(mlst_all %>% rename(st2 = st), by = c("id2" = "genome_id")) %>%
  filter(!is.na(st1), st1 == st2) %>%
  mutate(st = st1)

pair_table <- snv_with_st %>%
  transmute(id1, id2, st, snv = snv_post_noncore) %>%
  filter(!is.na(snv)) %>%
  mutate(
    case_control = case_when(
      snv <= snv_case_threshold ~ "case",
      snv > snv_control_threshold ~ "control_gt10",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(case_control))

print(count(pair_table, case_control))
