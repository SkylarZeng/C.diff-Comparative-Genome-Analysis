# Force R to use the correct Conda library path
# .libPaths("/home/yigezeng/.conda/envs/rstats/lib/R/library")

# Cognac Run
# Start Date: 07/01/25
# Organism: C.diff

library(tidyverse)
library(cognac)

Sys.info()
sessionInfo()

setwd("/nfs/turbo/umms-esnitkin/Project_Cdiff/Analysis/2025-sysbio-UM-transmission-Skylar/2025-07-01_cognac")

outDir = "/nfs/turbo/umms-esnitkin/Project_Cdiff/Analysis/2025-sysbio-UM-transmission-Skylar/2025-07-01_cognac/results/"

ids = read_csv("/nfs/turbo/umms-esnitkin/Project_Cdiff/Analysis/2025-sysbio-UM-transmission-Skylar/data/metadata/passed_ids.csv")
isolate_ids = ids$isolate_id

fasta_files = paste0("/nfs/turbo/umms-esnitkin/Project_Cdiff/Sequence_data/assembly/illumina/prokka/", isolate_ids, "/", isolate_ids, ".fna")
gff_files   = paste0("/nfs/turbo/umms-esnitkin/Project_Cdiff/Sequence_data/assembly/illumina/prokka/", isolate_ids, "/", isolate_ids, ".gff")

algnEnv = cognac(
  fastaFiles = fasta_files,
  featureFiles = gff_files,
  keepTempFiles = TRUE,
  outDir = outDir,
  minGeneNum = 500,
  maxMissGenes = 0.05,
  threadVal = 18,
  njTree = TRUE
)

system(
  paste(
    "FastTree <",
    algnEnv$aaAlgnPath,
    ">",
    file.path(outDir, "cdiff_fastTree.tree")
  )
)

algnEnv$distMat = CreateAlgnDistMat(
  "/nfs/turbo/umms-esnitkin/Project_Cdiff/Analysis/2025-sysbio-UM-transmission-Skylar/2025-07-01_cognac/results/concatenated_gene_aa_alignment.fasta",
  "raw"
)

write.table(
  algnEnv$distMat,
  file = "/nfs/turbo/umms-esnitkin/Project_Cdiff/Analysis/2025-sysbio-UM-transmission-Skylar/2025-07-01_cognac/results/cdiff_distMat.tsv"
)
