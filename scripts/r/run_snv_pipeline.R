#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# run_snv_pipeline.R
#
# Top‑level script to compute core, soft‑core and non‑core SNV distances from
# whole‑genome alignments.  This script uses helper functions defined in
# `snv_pipeline_skeleton.R` and wraps them with a command‑line interface.
# It produces QC summaries, distance matrices in both matrix and long format,
# and a suite of diagnostic plots.  The script is designed to be run on an
# HPC cluster (e.g. via SLURM) but can also be executed locally.
#
# Usage: Rscript run_snv_pipeline.R --pre <pre_alignment.fa> --post <post_alignment.fa> \
#        --outdir <output_directory> [options]
#
# Example:
#   Rscript run_snv_pipeline.R \
#       --pre data/pre_alignment.fa \
#       --post data/post_alignment.fa \
#       --gubbins_masked positions.txt \
#       --secondary_isolates secondary_isolates.txt \
#       --ref_id cdiff_W0022a_ref_genome \
#       --outdir results/analysis1
#
# Required packages: optparse, ggplot2, readr, dplyr, tibble, ape, Biostrings.
#
# -----------------------------------------------------------------------------

suppression <- suppressPackageStartupMessages
suppression({
  # Load only the necessary tidyverse packages.  optparse is intentionally not
  # loaded here to avoid dependency issues on systems where it is not installed.
  library(ggplot2)
  library(tibble)
  library(dplyr)
  library(readr)
})

# -----------------------------------------------------------------------------
# Source the helper functions
#
# We expect that `snv_pipeline_skeleton.R` is located in the same directory from
# which this script is invoked.  Rscript does not reliably set `sys.frame(1)$ofile`,
# so attempting to use it leads to errors.  Instead, we use the current
# working directory.  If the skeleton file is located elsewhere, you can edit
# this line to provide an absolute or relative path.
script_dir <- getwd()
source(file.path(script_dir, "snv_pipeline_skeleton.R"))

## ----------------------------------------------------------------------------
## Manual command‑line argument parsing
##
## This script avoids relying on the optparse package, which may not be
## installed on all systems.  Arguments are provided in the form
##   --argument_name value
## Boolean flags can be set by including them without a following value.

args <- commandArgs(trailingOnly = TRUE)

# Initialise defaults
opt <- list(
  pre = NULL,
  post = NULL,
  gubbins_masked = NULL,
  ref_id = NULL,
  secondary_isolates = NULL,
  missing_chars = "N,n,-,?,.,X,x",
  genome_missing_threshold = 0.10,
  soft_core_threshold = 0.90,
  soft_core_thresholds = "1.0,0.95,0.9,0.85,0.8",
  mask_post_based_on_pre = FALSE,
  outdir = "snv_pipeline_output",
  keep_invariant = FALSE,
  no_plots = FALSE
)

# Simple parser: loop over args, detect --key and assign next value if present
i <- 1
while (i <= length(args)) {
  arg <- args[[i]]
  if (startsWith(arg, "--")) {
    key <- sub("^--", "", arg)
    # Boolean flags: if next arg is another option or missing, set TRUE
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      opt[[key]] <- TRUE
      i <- i + 1
    } else {
      val <- args[[i + 1]]
      # Convert numeric arguments where appropriate
      if (key %in% c("genome_missing_threshold", "soft_core_threshold")) {
        opt[[key]] <- as.numeric(val)
      } else if (key %in% c("pre", "post", "gubbins_masked", "ref_id",
                            "secondary_isolates", "missing_chars", "soft_core_thresholds", "outdir")) {
        opt[[key]] <- val
      } else if (key %in% c("mask_post_based_on_pre", "keep_invariant", "no_plots")) {
        opt[[key]] <- tolower(val) %in% c("true", "t", "yes", "1")
      } else {
        # Unknown argument, store as character
        opt[[key]] <- val
      }
      i <- i + 2
    }
  } else {
    i <- i + 1
  }
}

# Validate required argument --pre
if (is.null(opt$pre)) {
  stop("--pre argument is required.  Please provide the path to the pre‑Gubbins alignment using --pre")
}

# Split missing characters into vector
missing_chars <- strsplit(opt$missing_chars, ",")[[1]]

# Create output directory
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# Helper function to save plots
save_plot <- function(plot_obj, filename) {
  ggsave(filename, plot = plot_obj, device = "pdf", width = 8, height = 5)
}

# Write summary tables to CSV
write_table <- function(df, filename) {
  readr::write_csv(df, filename)
}

# -----------------------------------------------------------------------------
# Main analysis
# -----------------------------------------------------------------------------

run_analysis <- function(aln_path, alignment_label, post_aln_path = NULL) {
  # Create subdirectory for this alignment
  prefix_dir <- file.path(opt$outdir, alignment_label)
  dir.create(prefix_dir, recursive = TRUE, showWarnings = FALSE)
  
  # -------------------------------------------------------------------------
  # Read alignment and perform QC
  aln_data <- read_alignment(aln_path)
  valid_mat <- valid_position_matrix(aln_data$dna, missing_chars = missing_chars)
  summary_stats <- summarise_alignment(aln_data$dna, missing_chars = missing_chars)
  # Save summary to CSV
  write_table(as_tibble(summary_stats), file.path(prefix_dir, paste0("alignment_summary_", alignment_label, ".csv")))
  
  # Genome missingness
  gm <- genome_missingness(valid_mat)
  write_table(gm, file.path(prefix_dir, paste0("genome_missingness_", alignment_label, ".csv")))
  if (!opt$no_plots) {
    p_gm <- plot_genome_missingness(gm, threshold = opt$genome_missing_threshold)
    save_plot(p_gm, file.path(prefix_dir, paste0("genome_missingness_", alignment_label, ".pdf")))
  }
  
  # Filter genomes on missingness
  filtered <- filter_genomes(aln_data, valid_mat, threshold = opt$genome_missing_threshold)
  # Optionally drop reference genome
  if (!is.null(opt$ref_id)) {
    if (opt$ref_id %in% rownames(filtered$aln)) {
      filtered$aln <- filtered$aln[rownames(filtered$aln) != opt$ref_id, , drop = FALSE]
      filtered$dna <- filtered$dna[names(filtered$dna) != opt$ref_id]
      filtered$removed_ids <- c(filtered$removed_ids, opt$ref_id)
    }
  }
  # Optionally drop secondary isolates before computing distances
  if (!is.null(opt$secondary_isolates)) {
    sec_ids <- readr::read_tsv(opt$secondary_isolates, col_names = FALSE, show_col_types = FALSE)$X1
    sec_ids <- unique(sec_ids)
    drop_ids <- intersect(rownames(filtered$aln), sec_ids)
    keep_ids <- setdiff(rownames(filtered$aln), sec_ids)
    filtered$aln <- filtered$aln[keep_ids, , drop = FALSE]
    filtered$dna <- filtered$dna[keep_ids]
    filtered$removed_ids <- c(filtered$removed_ids, drop_ids)
  }
  # Update kept_ids after additional filtering
  filtered$kept_ids <- rownames(filtered$aln)
  # Update valid matrix and stats after filtering
  valid_mat_filt <- valid_position_matrix(filtered$dna, missing_chars = missing_chars)
  gm_filt <- genome_missingness(valid_mat_filt)
  write_table(gm_filt, file.path(prefix_dir, paste0("genome_missingness_filtered_", alignment_label, ".csv")))
  if (!opt$no_plots) {
    p_gmf <- plot_genome_missingness(gm_filt, threshold = opt$genome_missing_threshold)
    save_plot(p_gmf, file.path(prefix_dir, paste0("genome_missingness_filtered_", alignment_label, ".pdf")))
  }
  sm <- site_missingness(valid_mat_filt)
  write_table(sm, file.path(prefix_dir, paste0("site_missingness_", alignment_label, ".csv")))
  if (!opt$no_plots) {
    p_sm <- plot_site_missingness(sm)
    save_plot(p_sm, file.path(prefix_dir, paste0("site_missingness_", alignment_label, ".pdf")))
  }
  # Soft‑core retention summary
  thresholds <- as.numeric(strsplit(opt$soft_core_thresholds, ",")[[1]])
  retention_tbl <- soft_core_retention(valid_mat_filt, thresholds)
  write_table(retention_tbl, file.path(prefix_dir, paste0("soft_core_retention_", alignment_label, ".csv")))
  if (!opt$no_plots) {
    p_ret <- plot_soft_core_retention(retention_tbl)
    save_plot(p_ret, file.path(prefix_dir, paste0("soft_core_retention_", alignment_label, ".pdf")))
  }
  # Identify core, soft‑core and compute distances
  results <- list()
  # Strict core
  core_sites <- identify_core_sites(valid_mat_filt, threshold = 1.0)
  dist_core <- compute_snv_distances(filtered$aln, core_sites)
  dist_core_long <- distance_matrix_to_long(dist_core)
  results$core <- list(matrix = dist_core, long = dist_core_long,
                       summary = summarise_distances(dist_core_long))
  # Soft‑core (default)
  soft_core_sites <- identify_core_sites(valid_mat_filt, threshold = opt$soft_core_threshold)
  dist_soft_core <- compute_snv_distances(filtered$aln, soft_core_sites)
  dist_soft_core_long <- distance_matrix_to_long(dist_soft_core)
  results$soft_core <- list(matrix = dist_soft_core, long = dist_soft_core_long,
                            summary = summarise_distances(dist_soft_core_long))
  # Non‑core
  dist_noncore <- compute_snv_distances(filtered$aln, core_sites = NULL)
  dist_noncore_long <- distance_matrix_to_long(dist_noncore)
  results$noncore <- list(matrix = dist_noncore, long = dist_noncore_long,
                          summary = summarise_distances(dist_noncore_long))
  # Write distance matrices and summaries
  for (method in names(results)) {
    method_dir <- file.path(prefix_dir, method)
    dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
    # Matrix as TSV
    write.table(results[[method]]$matrix,
                file = file.path(method_dir, paste0("distance_matrix_", method, ".tsv")),
                sep = "\t", quote = FALSE, col.names = NA)
    # Long format as CSV
    write_table(results[[method]]$long,
                file.path(method_dir, paste0("distance_long_", method, ".csv")))
    # Summary statistics
    write_table(results[[method]]$summary,
                file.path(method_dir, paste0("distance_summary_", method, ".csv")))
  }
  # Plot distance distributions
  if (!opt$no_plots) {
    p_dist <- plot_distance_distributions(list(core = results$core$long,
                                               soft_core = results$soft_core$long,
                                               noncore = results$noncore$long))
    save_plot(p_dist, file.path(prefix_dir, paste0("distance_distributions_", alignment_label, ".pdf")))
  }
  # Save kept and removed genomes lists
  write_table(tibble(id = filtered$kept_ids), file.path(prefix_dir, paste0("kept_genomes_", alignment_label, ".csv")))
  write_table(tibble(id = filtered$removed_ids), file.path(prefix_dir, paste0("removed_genomes_", alignment_label, ".csv")))
  # Return results for cross‑step comparison
  list(results = results, filtered = filtered, alignment_label = alignment_label,
       summary_stats = summary_stats, gm = gm, gm_filt = gm_filt, retention = retention_tbl)
}

# -----------------------------------------------------------------------------
# Pre‑alignment analysis
# -----------------------------------------------------------------------------
pre_result <- run_analysis(opt$pre, alignment_label = "pre")

# -----------------------------------------------------------------------------
# Gubbins masking summary (if provided)
# -----------------------------------------------------------------------------
if (!is.null(opt$gubbins_masked)) {
  # Summarise masking
  pre_len <- pre_result$summary_stats$n_sites
  mask_summary <- summarise_masking(opt$gubbins_masked, genome_length = pre_len)
  write_table(as_tibble(mask_summary), file.path(opt$outdir, "gubbins_masking_summary.csv"))
  # If block lengths are available, save them separately
  if (!is.null(mask_summary$block_lengths)) {
    write_table(tibble(block_length = mask_summary$block_lengths),
                file.path(opt$outdir, "gubbins_masking_block_lengths.csv"))
  }
}

# -----------------------------------------------------------------------------
# Post‑alignment analysis (if provided)
# -----------------------------------------------------------------------------
if (!is.null(opt$post)) {
  post_aln_path <- opt$post
  # Optionally mask post alignment based on pre missingness
  if (opt$mask_post_based_on_pre) {
    # Read post‑Gubbins alignment with Biostrings, drop known problematic
    # genome and convert to an ape DNAbin matrix for downstream functions.
    aln_post <- Biostrings::readDNAStringSet(post_aln_path)
    # Drop a specific genome that may cause issues (if present)
    if ("CDIF_SB_1448" %in% names(aln_post)) {
      aln_post <- aln_post[names(aln_post) != "CDIF_SB_1448"]
    }
    cat("Loaded", length(aln_post), "genomes; alignment length =", width(aln_post)[1], "\n")
    # Convert Biostrings -> ape DNAbin MATRIX (important for matrix subsetting)
    aln_post_ape <- ape::as.matrix.DNAbin(ape::as.DNAbin(aln_post))
    post_data <- list(dna = aln_post, dnabin = aln_post_ape)
    masked <- mask_post_based_on_pre(pre_result$filtered$dna, post_data$dna, missing_chars = missing_chars)
    # Write masked alignment to a temporary file
    masked_fa <- file.path(opt$outdir, "masked_post_alignment.fa")
    Biostrings::writeXStringSet(masked$masked_post, masked_fa)
    post_aln_path <- masked_fa
  }
  post_result <- run_analysis(post_aln_path, alignment_label = "post")
  # Cross‑step comparison summary
  comparison <- tibble(
    metric = c("n_genomes", "n_sites", "prop_missing_mean", "core_site_count", "soft_core_site_count", "noncore_pairs", "core_pairs", "soft_core_pairs"),
    pre = c(
      length(pre_result$filtered$kept_ids),
      pre_result$summary_stats$n_sites,
      mean(pre_result$gm_filt$prop_missing),
      sum(identify_core_sites(valid_position_matrix(pre_result$filtered$dna, missing_chars = missing_chars), threshold = 1.0)),
      sum(identify_core_sites(valid_position_matrix(pre_result$filtered$dna, missing_chars = missing_chars), threshold = opt$soft_core_threshold)),
      nrow(pre_result$results$noncore$long),
      nrow(pre_result$results$core$long),
      nrow(pre_result$results$soft_core$long)
    ),
    post = c(
      length(post_result$filtered$kept_ids),
      post_result$summary_stats$n_sites,
      mean(post_result$gm_filt$prop_missing),
      sum(identify_core_sites(valid_position_matrix(post_result$filtered$dna, missing_chars = missing_chars), threshold = 1.0)),
      sum(identify_core_sites(valid_position_matrix(post_result$filtered$dna, missing_chars = missing_chars), threshold = opt$soft_core_threshold)),
      nrow(post_result$results$noncore$long),
      nrow(post_result$results$core$long),
      nrow(post_result$results$soft_core$long)
    )
  )
  write_table(comparison, file.path(opt$outdir, "cross_step_comparison.csv"))
}

message("SNV distance pipeline completed. Outputs are in: ", normalizePath(opt$outdir))