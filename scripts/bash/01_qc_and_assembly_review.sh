#!/bin/bash
set -euo pipefail

# Step 1: QC and assembly review before downstream analysis.
# This condenses the QC stage that relied on QCD/pubQCD-style workflows.

# Dry-run a download/QC workflow to verify expected inputs and outputs.
snakemake -s workflow/download_genomes.smk --dryrun -p

# Example cluster submission for a pubQCD-style run.
sbatch bash_script_to_download_and_run_raw_reads_pubQCD.sbat

# Typical review rules used before keeping genomes for later analyses:
# - contig count between 10 and 500
# - average coverage greater than 20X
# - assembly size no greater than 7 Mb
