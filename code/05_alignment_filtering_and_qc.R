library(here)

# Step 5: alignment QC, missingness review, and genome filtering.
# This uses helper functions from the modular SNV pipeline.

source(here("scripts", "r", "snv_pipeline_skeleton.R"))

aln <- read_alignment(
  here(
    "variant_calling", "clade2_R20291", "phylo_analysis", "gubbins",
    "clade2_R20291_gubbins_masked.fa"
  )
)

valid_mat <- valid_position_matrix(aln$dna)
summary_stats <- summarise_alignment(aln$dna)
genome_qc <- genome_missingness(valid_mat)

filtered <- filter_genomes(aln, valid_mat, threshold = 0.10)
filtered_valid <- valid_position_matrix(filtered$dna)

retention_tbl <- soft_core_retention(
  filtered_valid,
  thresholds = c(1.0, 0.95, 0.90, 0.85)
)

print(summary_stats)
print(head(genome_qc))
print(retention_tbl)
