#!/bin/bash
set -euo pipefail

# Step 2: reference-based variant calling with snpkit.
# These two commands capture the call -> parse flow used before custom SNV work.

python snpkit/snpkit.py \
  -type PE \
  -readsdir /Path-To-Your/test_readsdir/ \
  -outdir /Path/test_output_core/ \
  -analysis output_prefix \
  -index KPNIH1 \
  -steps call \
  -cluster cluster \
  -scheduler SLURM \
  -clean

python snpkit/snpkit.py \
  -type PE \
  -readsdir /Path-To-Your/test_readsdir/ \
  -outdir /Path/test_output_core/ \
  -analysis output_prefix \
  -index reference.fasta \
  -steps parse \
  -cluster cluster \
  -gubbins yes \
  -scheduler SLURM
