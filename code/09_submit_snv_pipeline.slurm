#!/bin/bash
#SBATCH --job-name=snv_pipeline
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

set -euo pipefail

# Step 9: submit the modular SNV pipeline to an HPC scheduler.

PRE="variant_calling/clade2_R20291/phylo_analysis/genome_alignment/pre_alignment.fa"
POST="variant_calling/clade2_R20291/phylo_analysis/gubbins/post_alignment.fa"
GUBBINS_MASK="variant_calling/clade2_R20291/phylo_analysis/gubbins/masked_positions.txt"
SECONDARY="data/metadata/secondary_isolate_ids.txt"
OUTDIR="snv_pipeline_output"

Rscript scripts/r/run_snv_pipeline.R \
  --pre "$PRE" \
  --post "$POST" \
  --gubbins_masked "$GUBBINS_MASK" \
  --ref_id cdiff_R20291_ref_genome \
  --secondary_isolates "$SECONDARY" \
  --outdir "$OUTDIR" \
  --genome_missing_threshold 0.10 \
  --soft_core_threshold 0.90 \
  --soft_core_thresholds 1.0,0.95,0.9,0.85,0.8 \
  --mask_post_based_on_pre
