#!/bin/bash
#SBATCH --job-name=snv_pipeline
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=you@example.com

set -euo pipefail

# Load R (adjust module name for your cluster) or activate conda env
module purge
module load R || true
# OR: source /path/to/conda.sh && conda activate my_r_env

# Paths (adjust if necessary)
PRE="/Users/skylar/Desktop/Snitkin Lab/DataAnalysis/2025-sysbio-UM-transmission/variant_calling/clade2_R20291/phylo_analysis/genome_alignment/2025_07_16_09_49_34_Cdiff_O27_R20291_genome_aln_w_alt_allele_unmapped.fa"
POST="/Users/skylar/Desktop/Snitkin Lab/DataAnalysis/2025-sysbio-UM-transmission/variant_calling/clade2_R20291/phylo_analysis/gubbins/clade2_R20291_gubbins_masked.fa"
GUBBINS_MASK="/Users/skylar/Desktop/Snitkin Lab/DataAnalysis/2025-sysbio-UM-transmission/variant_calling/clade2_R20291/phylo_analysis/gubbins/clade2_R20291_masked_recomb_positions.txt"
SECONDARY="/Users/skylar/Desktop/Snitkin Lab/DataAnalysis/2025-sysbio-UM-transmission/Metadata/secondary_isolate_ids.txt"
OUTDIR="/Users/skylar/Desktop/Snitkin Lab/DataAnalysis/2025-sysbio-UM-transmission/snv_pipeline_output"

# Rscript command
Rscript 2026_analysis/run_snv_pipeline.R \
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

# End of script
