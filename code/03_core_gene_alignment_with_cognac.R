library(tidyverse)
library(cognac)

# Step 3: core-gene alignment and distance generation with cognac.
# Adapted from the project-specific 2025-07-01 cognac workflow.

ids <- read_csv("data/metadata/passed_ids.csv", show_col_types = FALSE)
isolate_ids <- ids$isolate_id
out_dir <- "results/cognac"

fasta_files <- paste0(
  "/nfs/turbo/umms-esnitkin/Project_Cdiff/Sequence_data/assembly/illumina/prokka/",
  isolate_ids, "/", isolate_ids, ".fna"
)

gff_files <- paste0(
  "/nfs/turbo/umms-esnitkin/Project_Cdiff/Sequence_data/assembly/illumina/prokka/",
  isolate_ids, "/", isolate_ids, ".gff"
)

algn_env <- cognac(
  fastaFiles = fasta_files,
  featureFiles = gff_files,
  keepTempFiles = TRUE,
  outDir = out_dir,
  minGeneNum = 500,
  maxMissGenes = 0.05,
  threadVal = 18,
  njTree = TRUE
)

dist_mat <- CreateAlgnDistMat(algn_env$aaAlgnPath, "raw")
write.table(dist_mat, file = file.path(out_dir, "cdiff_distMat.tsv"))
