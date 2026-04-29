# SNV Distance Calculation Pipeline (skeleton)

# This script defines a set of functions to perform QC, filtering, and SNV distance
# calculations on whole‑genome alignments.  It is designed to be modular and
# parameterised so that it can be reused across different datasets.  The script
# can be sourced from another R script or executed directly with command line
# arguments (to be added in a wrapper).  The functions below do not produce
# plots or write files by themselves; companion functions later in the script
# generate summaries and figures.  Additional code should be added to parse
# command line options, set up output directories, and call these functions in
# the desired order.

## Load required packages
suppressPackageStartupMessages({
  library(ape)         # for dist.dna and DNAbin handling
  library(Biostrings)  # for reading FASTA and missingness masking
  library(tidyverse)   # for data manipulation and plotting
  library(rlang)       # for error handling
})

# --------------------------------------------------------------------------------
# Helper functions
# --------------------------------------------------------------------------------

#' Read an alignment from a FASTA file
#'
#' @param aln_path Path to the FASTA alignment file
#' @return A DNAStringSet and DNAbin representation of the alignment
read_alignment <- function(aln_path) {
  message("Reading alignment: ", aln_path)
  # Use Biostrings to read the FASTA.  This returns a DNAStringSet, which we
  # keep for position‑level operations.  Convert to DNAbin for distance
  # calculations.
  dna_strings <- Biostrings::readDNAStringSet(aln_path)
  if (length(dna_strings) == 0) {
    abort(paste("No sequences were found in", aln_path))
  }
  # Convert to DNAbin for ape
  dnabin <- ape::as.DNAbin(dna_strings)
  list(dna = dna_strings, dnabin = dnabin)
}

#' Identify valid (non‑missing) positions for each genome
#'
#' Missing characters are defined as N/n, gap ("-"), ambiguous bases ("X"),
#' question mark ("?"), or dot (".").  The default set can be overridden.
#'
#' @param dna_strings A DNAStringSet
#' @param missing_chars Character vector of characters considered missing
#' @return A logical matrix [n_positions x n_genomes] indicating non‑missing data
valid_position_matrix <- function(dna_strings,
                                  missing_chars = c("N", "n", "-", "?", ".", "X", "x")) {
  # Build a consensus matrix that counts characters at each position across
  # genomes.  Biostrings::consensusMatrix returns a matrix with rows
  # corresponding to characters and columns to positions.
  cm <- Biostrings::consensusMatrix(dna_strings, baseOnly = FALSE)
  # Identify rows corresponding to missing characters
  miss_rows <- intersect(rownames(cm), missing_chars)
  # Number of genomes
  n_genomes <- length(dna_strings)
  # A position is valid for a genome if the genome's base at that position is
  # not in missing_chars.  We build this matrix by iterating through genomes.
  n_pos <- width(dna_strings)[1]
  valid <- matrix(TRUE, nrow = n_pos, ncol = n_genomes)
  seqs <- as.character(dna_strings)
  for (i in seq_len(n_genomes)) {
    chars <- unlist(strsplit(seqs[i], split = ""))
    valid[, i] <- !(chars %in% missing_chars)
  }
  colnames(valid) <- names(dna_strings)
  valid
}

#' Summarise alignment at a high level
#'
#' Computes the number of genomes, alignment length, number of invariant sites,
#' number of variable sites, and overall missingness.  Invariant sites are
#' positions where there is at most one unique unambiguous base.  The definition
#' follows guidelines from the Core‑SNP‑filter tool, which treats sites as
#' invariant if the number of unique unambiguous bases is zero or one【289074610468718†L336-L351】.
#'
#' @param dna_strings A DNAStringSet
#' @param missing_chars Characters considered missing
#' @return A list with summary statistics
summarise_alignment <- function(dna_strings,
                                missing_chars = c("N", "n", "-", "?", ".", "X", "x")) {
  n_genomes <- length(dna_strings)
  n_sites <- width(dna_strings)[1]
  # Build consensus matrix to count unambiguous bases at each position
  cm <- Biostrings::consensusMatrix(dna_strings, baseOnly = FALSE)
  # Identify unambiguous base rows (A,C,G,T and lower case)
  unambig_bases <- c("A", "C", "G", "T", "a", "c", "g", "t")
  # For each site, count the number of distinct unambiguous bases
  unambig_counts <- colSums(cm[intersect(rownames(cm), unambig_bases), ] > 0)
  invariant_sites <- sum(unambig_counts <= 1)
  variable_sites <- n_sites - invariant_sites
  # Compute missingness per site
  miss_rows <- intersect(rownames(cm), missing_chars)
  total_missing <- sum(cm[miss_rows, ])
  total_bases <- n_genomes * n_sites
  prop_missing <- total_missing / total_bases
  list(
    n_genomes = n_genomes,
    n_sites = n_sites,
    invariant_sites = invariant_sites,
    variable_sites = variable_sites,
    prop_missing = prop_missing
  )
}

#' Compute genome‑level missingness
#'
#' Calculates the number and proportion of missing positions for each genome.
#'
#' @param valid_matrix Logical matrix from valid_position_matrix
#' @return A data frame with columns id, n_missing, prop_missing
genome_missingness <- function(valid_matrix) {
  n_sites <- nrow(valid_matrix)
  n_missing_per_genome <- colSums(!valid_matrix)
  tibble(
    id = colnames(valid_matrix),
    n_missing = n_missing_per_genome,
    prop_missing = n_missing_per_genome / n_sites
  )
}

#' Compute position‑level missingness
#'
#' Calculates the number and proportion of genomes missing at each site.
#'
#' @param valid_matrix Logical matrix from valid_position_matrix
#' @return A data frame with columns position (1‑based), n_missing, prop_missing
site_missingness <- function(valid_matrix) {
  n_genomes <- ncol(valid_matrix)
  n_missing_per_site <- rowSums(!valid_matrix)
  tibble(
    position = seq_len(nrow(valid_matrix)),
    n_missing = n_missing_per_site,
    prop_missing = n_missing_per_site / n_genomes
  )
}

#' Filter genomes based on missingness threshold
#'
#' Removes genomes whose proportion of missing sites exceeds a threshold.  The
#' threshold is the maximum allowed missingness, so genomes with prop_missing
#' strictly greater than this value are excluded.  The function returns the
#' filtered DNAbin and DNAStringSet, along with a vector of removed genome IDs.
#' The missingness threshold defaults to 0.10 (10%) based on studies that
#' recommend excluding genomes with more than 10% ambiguous bases【289074610468718†L372-L377】.
#'
#' @param aln_list A list containing dnabin and dna elements as returned by
#'   read_alignment()
#' @param valid_matrix Logical matrix from valid_position_matrix
#' @param threshold Proportion of missingness allowed (default 0.10)
#' @return A list with elements aln (filtered DNAbin), dna (filtered DNAStringSet),
#'   kept_ids, removed_ids
filter_genomes <- function(aln_list, valid_matrix, threshold = 0.10) {
  gm <- genome_missingness(valid_matrix)
  removed <- gm %>% filter(prop_missing > threshold)
  kept <- gm %>% filter(prop_missing <= threshold)
  # Subset the alignments
  # Subset DNAStringSet by kept ids (works by name/index)
  filt_dna <- aln_list$dna[kept$id]
  # Ensure DNAbin is returned in a matrix-like form.  Some conversions yield a
  # list-like DNAbin where two‑dimensional subsetting fails; to be robust we
  # regenerate DNAbin from the filtered DNAStringSet.
  filt_dnabin <- ape::as.DNAbin(filt_dna)
  # Ensure DNAbin is matrix‑like so callers can subset with [rows, cols].
  if (!is.null(filt_dnabin)) {
    if (!is.matrix(filt_dnabin)) {
      # single sequence -> make a 1 x L matrix
      filt_dnabin <- matrix(filt_dnabin, nrow = 1)
    }
    # Ensure rownames match sequence names
    if (is.null(rownames(filt_dnabin))) {
      rownames(filt_dnabin) <- names(filt_dna)
    }
    class(filt_dnabin) <- c("DNAbin", class(filt_dnabin)[!class(filt_dnabin) %in% "DNAbin"]) 
  }

  list(aln = filt_dnabin, dna = filt_dna, kept_ids = kept$id, removed_ids = removed$id)
}

#' Identify core and soft‑core sites
#'
#' Returns a logical vector indicating whether each position is part of the core
#' according to a given threshold.  A site is considered core if the fraction of
#' unambiguous bases exceeds the threshold (e.g. 1.0 for strict core, 0.9 for
#' 90% core).  All non‑A/C/G/T characters are treated equally【289074610468718†L336-L351】.
#'
#' @param valid_matrix Logical matrix from valid_position_matrix
#' @param threshold Minimum fraction of genomes with valid bases to retain a site
#'   (e.g. 1.0 = strict core, 0.9 = allow up to 10% missing)
#' @return A logical vector of length n_positions indicating core sites
identify_core_sites <- function(valid_matrix, threshold = 1.0) {
  n_genomes <- ncol(valid_matrix)
  prop_valid <- rowMeans(valid_matrix)
  prop_valid >= threshold
}

#' Compute pairwise SNV distances
#'
#' This function calculates distances using ape::dist.dna.  For core and
#' soft‑core distances, the alignment is subset to the retained sites and
#' pairwise.deletion is set to FALSE so that only those sites are considered.
#' For non‑core distances, the full alignment is used with pairwise.deletion =
#' TRUE so that only overlapping positions between each pair are counted.  The
#' distance model "N" returns the number of differing sites between sequences
#'【125490742559221†L88-L92】.  The function returns a distance matrix.
#'
#' @param dnabin A DNAbin object (filtered genomes)
#' @param core_sites Optional logical vector indicating which sites to include
#'   (for core or soft‑core).  If NULL, non‑core distances are computed
#'   pairwise with pairwise.deletion = TRUE.
#' @return A numeric matrix of pairwise SNV distances
compute_snv_distances <- function(dnabin, core_sites = NULL) {
  if (!is.null(core_sites)) {
    # Subset the alignment by sites (columns)
    if (length(core_sites) != ncol(dnabin)) {
      abort("Length of core_sites does not match alignment width")
    }
    # Keep only columns that are TRUE
    sub_dnabin <- dnabin[, core_sites, drop = FALSE]
    # Use pairwise.deletion = FALSE because sites with missing data should have
    # been filtered out by core_sites definition
    dist_obj <- ape::dist.dna(sub_dnabin, model = "N", pairwise.deletion = FALSE)
  } else {
    # Non‑core distances: use the full alignment and allow pairwise deletion
    dist_obj <- ape::dist.dna(dnabin, model = "N", pairwise.deletion = TRUE)
  }
  as.matrix(dist_obj)
}

#' Convert a distance matrix to a long format table
#'
#' Returns a tibble with columns id1, id2, snvs.  Pairs where id1 == id2 are
#' dropped, and each pair appears only once (id1 < id2).
#'
#' @param dist_mat Numeric matrix of pairwise distances
#' @return A tibble with pairwise distances
distance_matrix_to_long <- function(dist_mat) {
  rn <- rownames(dist_mat)
  cn <- colnames(dist_mat)
  # Expand grid of indices
  df <- expand.grid(id1 = rn, id2 = cn, stringsAsFactors = FALSE)
  df$snvs <- dist_mat[cbind(match(df$id1, rn), match(df$id2, cn))]
  df %>%
    filter(id1 < id2) %>%
    transmute(id1, id2, snvs)
}

#' Summarise pairwise distances
#'
#' Computes basic statistics (min, max, mean, median, quartiles) for the
#' distribution of pairwise distances.  This can be useful to compare
#' core/non‑core/soft‑core results across steps.
#'
#' @param dist_long A tibble produced by distance_matrix_to_long
#' @return A tibble with summary statistics
summarise_distances <- function(dist_long) {
  tibble(
    n_pairs = nrow(dist_long),
    min = min(dist_long$snvs, na.rm = TRUE),
    q1 = quantile(dist_long$snvs, 0.25, na.rm = TRUE),
    median = median(dist_long$snvs, na.rm = TRUE),
    mean = mean(dist_long$snvs, na.rm = TRUE),
    q3 = quantile(dist_long$snvs, 0.75, na.rm = TRUE),
    max = max(dist_long$snvs, na.rm = TRUE)
  )
}

#' Evaluate soft‑core retention across thresholds
#'
#' Generates a summary table of retained site counts for a range of core
#' thresholds (e.g. 100% down to 80%).  This helps visualise how allowing
#' increasing missingness affects the number of sites retained, as shown in
#' large genomic studies where soft‑core thresholds preserved more sites than
#' strict cores【289074610468718†L224-L236】.  The function returns a tibble with
#' columns threshold and n_sites.
#'
#' @param valid_matrix Logical matrix from valid_position_matrix
#' @param thresholds Numeric vector of thresholds (values between 0 and 1)
#' @return A tibble summarising retained site counts
soft_core_retention <- function(valid_matrix, thresholds = seq(1.0, 0.80, by = -0.05)) {
  map_dfr(thresholds, ~{
    core_sites <- identify_core_sites(valid_matrix, threshold = .x)
    tibble(threshold = .x, n_sites = sum(core_sites))
  }) %>%
    arrange(desc(threshold))
}

#' Mask positions in a post‑Gubbins alignment based on pre‑alignment missingness
#'
#' When recombination has been masked in a post‑Gubbins alignment, you may
#' optionally want to remove positions that had any missingness in the pre‑Gubbins
#' alignment.  This function takes both the pre‑ and post‑alignment DNAStringSets,
#' identifies positions with any missing data in the pre‑alignment, and masks
#' those positions in the post‑alignment by replacing them with "N".  This
#' replicates behaviour used in prior analyses to ensure that post‑Gubbins core
#' sites are a subset of pre‑core sites.  The positions are returned for
#' potential downstream summaries.
#'
#' @param pre_dna A DNAStringSet of the pre‑Gubbins alignment
#' @param post_dna A DNAStringSet of the post‑Gubbins alignment
#' @param missing_chars Characters considered missing in the pre‑alignment
#' @return A list containing masked_post (DNAStringSet) and positions_to_mask
mask_post_based_on_pre <- function(pre_dna, post_dna,
                                   missing_chars = c("N", "n", "-", "?", ".", "X", "x")) {
  cm_pre <- Biostrings::consensusMatrix(pre_dna, baseOnly = FALSE)
  miss_rows <- intersect(rownames(cm_pre), missing_chars)
  # Positions with at least one missing data in pre
  positions_to_mask <- which(colSums(cm_pre[miss_rows, , drop = FALSE]) > 0)
  # Mask these positions in post
  masked_post <- post_dna
  for (i in seq_along(post_dna)) {
    seq_vec <- strsplit(as.character(post_dna[i]), split = "")[[1]]
    seq_vec[positions_to_mask] <- "N"
    masked_post[i] <- Biostrings::DNAStringSet(paste0(seq_vec, collapse = ""))
  }
  list(masked_post = masked_post, positions_to_mask = positions_to_mask)
}

#' Summarise recombination masking from a Gubbins output
#'
#' Parses a file listing masked recombination positions (e.g. from
#' gubbins), returning the number of positions masked, the proportion of
#' the genome masked, and optionally the distribution of block lengths.  If
#' block boundaries can be inferred, block_lengths can be calculated; otherwise
#' only counts are returned.  The input file is expected to contain either a
#' column of positions or ranges.  This function should be adapted to the
#' specific format of your GFF/VCF files.
#'
#' @param mask_file Path to the file listing recombination‑masked positions
#' @param genome_length Length of the alignment (number of sites)
#' @return A list with n_masked, prop_masked and optional block_lengths
summarise_masking <- function(mask_file, genome_length) {
  # Read lines; assume tab or space delimited, ignore comments
  mask_data <- readr::read_tsv(mask_file, comment = "#", col_names = FALSE, show_col_types = FALSE)
  # Attempt to detect whether file contains ranges (start/end) or single positions
  if (ncol(mask_data) >= 2) {
    starts <- mask_data[[1]]
    ends <- mask_data[[2]]
    # If any end < start, swap
    ends <- pmax(starts, ends)
    masked_positions <- unlist(Map(seq, starts, ends))
    block_lengths <- ends - starts + 1
  } else {
    masked_positions <- mask_data[[1]]
    block_lengths <- NULL
  }
  n_masked <- length(unique(masked_positions))
  prop_masked <- n_masked / genome_length
  list(n_masked = n_masked, prop_masked = prop_masked, block_lengths = block_lengths)
}

# --------------------------------------------------------------------------------
# Plotting functions (examples)
# --------------------------------------------------------------------------------

#' Plot genome missingness distribution
#'
#' Generates a histogram or density plot of genome‑level missingness.  The figure
#' helps identify outlier genomes with high proportions of missing data and
#' choose appropriate filtering thresholds.
#'
#' @param gm Data frame from genome_missingness()
#' @param threshold Optional threshold line to display
#' @return A ggplot object
plot_genome_missingness <- function(gm, threshold = NULL) {
  p <- ggplot(gm, aes(x = prop_missing)) +
    geom_histogram(binwidth = 0.01, fill = "steelblue", color = "black") +
    labs(x = "Proportion of missing sites per genome", y = "Count of genomes",
         title = "Genome-level missingness") +
    theme_minimal()
  if (!is.null(threshold)) {
    p <- p + geom_vline(xintercept = threshold, linetype = "dashed", color = "red") +
      annotate("text", x = threshold, y = Inf, label = paste0("Threshold = ", threshold), vjust = -0.5, hjust = 1, color = "red")
  }
  p
}

#' Plot per‑site missingness distribution
#'
#' Generates a histogram of the proportion of genomes missing at each site.  This
#' figure highlights highly missing positions and can be used to determine
#' appropriate soft‑core thresholds.
#'
#' @param sm Data frame from site_missingness()
#' @return A ggplot object
plot_site_missingness <- function(sm) {
  ggplot(sm, aes(x = prop_missing)) +
    geom_histogram(binwidth = 0.01, fill = "darkgreen", color = "black") +
    labs(x = "Proportion of genomes missing per site", y = "Count of sites",
         title = "Per‑site missingness distribution") +
    theme_minimal()
}

#' Plot soft‑core retention curve
#'
#' Plots the number of retained sites as a function of core threshold.  This
#' curve illustrates the trade‑off between missingness tolerance and alignment
#' length【289074610468718†L224-L236】.  Thresholds should be sorted in descending order.
#'
#' @param retention_tbl A tibble from soft_core_retention()
#' @return A ggplot object
plot_soft_core_retention <- function(retention_tbl) {
  ggplot(retention_tbl, aes(x = threshold, y = n_sites)) +
    geom_line() + geom_point() +
    scale_x_reverse(breaks = retention_tbl$threshold) +
    labs(x = "Core threshold (proportion of genomes with valid base)", y = "Number of sites retained",
         title = "Soft‑core retention curve") +
    theme_minimal()
}

#' Plot distribution of pairwise SNV distances
#'
#' Plots a density or histogram of SNV distances.  Useful to compare core,
#' non‑core and soft‑core distance distributions or pre‑ vs post‑Gubbins.  A
#' grouping variable can be supplied to facet or colour the distributions.
#'
#' @param dist_list A named list where each element is a tibble of pairwise
#'   distances (id1, id2, snvs)
#' @return A ggplot object
plot_distance_distributions <- function(dist_list) {
  # Bind all into one data frame with a method column
  df <- bind_rows(lapply(names(dist_list), function(nm) mutate(dist_list[[nm]], method = nm)))
  ggplot(df, aes(x = snvs, colour = method)) +
    geom_density() +
    labs(x = "Pairwise SNV distance", y = "Density",
         title = "Distribution of pairwise SNV distances") +
    theme_minimal()
}

# --------------------------------------------------------------------------------
# Notes
# --------------------------------------------------------------------------------

# The functions above provide building blocks for constructing a complete pipeline.
# A typical analysis might proceed as follows:
#
# 1. Read pre‑Gubbins alignment using read_alignment().
# 2. Compute valid_position_matrix() to determine which sites are non‑missing.
# 3. Summarise alignment with summarise_alignment().
# 4. Compute genome and site missingness and plot distributions.
# 5. Filter genomes using filter_genomes() based on a chosen missingness
#    threshold (default 10% missingness is a common cutoff【289074610468718†L372-L377】).
# 6. For the filtered alignment, identify core sites with threshold = 1.0
#    (strict core), soft‑core sites with threshold = 0.9 (allow 10% missing),
#    etc.
# 7. Compute SNV distances using compute_snv_distances() for core,
#    non‑core (pass NULL for core_sites), and soft‑core.  Convert to long format
#    with distance_matrix_to_long() and summarise with summarise_distances().
# 8. If recombination masking files are provided, call summarise_masking() and
#    visualise masked positions.
# 9. Read post‑Gubbins alignment and repeat steps 2–7.  Optionally mask post
#    alignment at positions with missing data in the pre‑alignment using
#    mask_post_based_on_pre().
# 10. Compare pre‑ and post‑Gubbins summaries and distance distributions to
#     identify potential artifacts.
# 11. Save all tables and figures to organised output directories and write
#     summary logs detailing the parameters used and genomes/sites retained.
#
# The actual implementation will need to handle command line arguments,
# configuration files, file I/O and parallelism (e.g. running under SLURM).  Those
# details are left as an exercise for integration into your computing
# environment.  The functions above should help structure the core analytical
# tasks.
