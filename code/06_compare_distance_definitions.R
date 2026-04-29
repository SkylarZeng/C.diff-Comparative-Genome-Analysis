library(tidyverse)
library(readr)
library(here)

# Step 6: compare strict core, soft-core, and non-core SNV distances.

source(here("scripts", "r", "snv_pipeline_skeleton.R"))

aln <- read_alignment(
  here(
    "variant_calling", "clade2_R20291", "phylo_analysis", "gubbins",
    "clade2_R20291_gubbins_masked_var_sites.fa"
  )
)

secondary_isolates <- read_tsv(
  here("data", "metadata", "secondary_isolate_ids.txt"),
  col_names = FALSE,
  show_col_types = FALSE
)$X1

valid_mat <- valid_position_matrix(aln$dna)

core_long <- distance_matrix_to_long(
  compute_snv_distances(aln$dnabin, identify_core_sites(valid_mat, 1.0))
)

softcore_long <- distance_matrix_to_long(
  compute_snv_distances(aln$dnabin, identify_core_sites(valid_mat, 0.95))
)

noncore_long <- distance_matrix_to_long(
  compute_snv_distances(aln$dnabin, core_sites = NULL)
)

drop_secondary <- function(df) {
  df %>% filter(!(id1 %in% secondary_isolates), !(id2 %in% secondary_isolates))
}

snv_comparison <- drop_secondary(core_long) %>%
  rename(snv_core = snvs) %>%
  inner_join(drop_secondary(softcore_long) %>% rename(snv_softcore = snvs), by = c("id1", "id2")) %>%
  inner_join(drop_secondary(noncore_long) %>% rename(snv_noncore = snvs), by = c("id1", "id2")) %>%
  mutate(
    delta_soft_vs_core = snv_softcore - snv_core,
    delta_noncore_vs_soft = snv_noncore - snv_softcore
  )

print(summary(snv_comparison$delta_soft_vs_core))
print(summary(snv_comparison$delta_noncore_vs_soft))
