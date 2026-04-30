args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_all[grep(file_arg, args_all)])

if (length(script_path) == 0) {
  script_dir <- getwd()
} else {
  script_dir <- dirname(normalizePath(script_path))
}

base_dir <- normalizePath(file.path(script_dir, ".."))
cluster_dir <- normalizePath(file.path(base_dir, ".."))
output_dir <- file.path(base_dir, "outputs")
desc_dir <- file.path(output_dir, "descriptive_statistics")
viz_dir <- file.path(output_dir, "visualization")

metadata_path <- file.path(desc_dir, "metadata_case_unformed_clade1_clade2.csv")
clade1_snv_path <- file.path(cluster_dir, "SNV", "clade1_post_noncore_snv_long.csv")
clade2_snv_path <- file.path(cluster_dir, "SNV", "clade2_post_noncore_snv_long.csv")

if (!dir.exists(desc_dir)) {
  dir.create(desc_dir, recursive = TRUE)
}
if (!dir.exists(viz_dir)) {
  dir.create(viz_dir, recursive = TRUE)
}

required_pkgs <- c("tidyverse", "igraph", "ggplot2")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(ggplot2)
})

snv_cutoff <- 2

save_plot_pdf <- function(plot_obj, filename, width, height) {
  ggsave(
    filename = filename,
    plot = plot_obj,
    device = "pdf",
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
}

make_st_palette <- function(st_values) {
  st_values <- sort(unique(st_values))
  if (length(st_values) == 0) {
    return(character())
  }

  hues <- seq(15, 375, length.out = length(st_values) + 1)[seq_along(st_values)]
  cols <- hcl(h = hues, c = 100, l = 65)
  names(cols) <- st_values
  cols
}

get_saved_st_palette <- function(clade_label, st_values) {
  key_path <- if (clade_label == "Clade 1") {
    file.path(desc_dir, "clade1_strain_diversity_color_key.csv")
  } else {
    file.path(desc_dir, "clade2_strain_diversity_color_key.csv")
  }

  if (file.exists(key_path)) {
    key_df <- read.csv(key_path, stringsAsFactors = FALSE) %>%
      transmute(
        ST = trimws(as.character(ST_label)),
        color = trimws(as.character(color))
      ) %>%
      filter(!is.na(ST), ST != "", !is.na(color), color != "") %>%
      distinct(ST, .keep_all = TRUE)

    palette <- setNames(key_df$color, key_df$ST)
    missing_st <- setdiff(st_values, names(palette))

    if (length(missing_st) > 0) {
      fallback <- make_st_palette(missing_st)
      palette <- c(palette, fallback)
    }

    return(palette[unique(c(intersect(names(palette), st_values), missing_st))])
  }

  make_st_palette(st_values)
}

metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
names(metadata) <- trimws(names(metadata))

metadata_clean <- metadata %>%
  transmute(
    genome_id = trimws(as.character(genome_id)),
    patient_id = trimws(as.character(patient_id)),
    collection_date = trimws(as.character(collection_date)),
    ST = trimws(as.character(ST)),
    clade_group = trimws(as.character(clade_group))
  ) %>%
  mutate(
    patient_id = na_if(patient_id, ""),
    ST = if_else(is.na(ST) | ST == "" | ST == "NA", "ST NA", paste0("ST ", ST))
  ) %>%
  filter(
    !is.na(genome_id),
    genome_id != "",
    clade_group %in% c("Clade 1", "Clade 2")
  ) %>%
  distinct(genome_id, .keep_all = TRUE)

read_snv_long <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE) %>%
    transmute(
      id1 = trimws(as.character(id1)),
      id2 = trimws(as.character(id2)),
      snvs = suppressWarnings(as.numeric(snvs))
    ) %>%
    filter(
      !is.na(id1),
      id1 != "",
      !is.na(id2),
      id2 != "",
      !is.na(snvs)
    ) %>%
    mutate(
      a = pmin(id1, id2),
      b = pmax(id1, id2)
    ) %>%
    select(id1 = a, id2 = b, snvs)
}

build_cluster_outputs <- function(clade_label, snv_path) {
  metadata_clade <- metadata_clean %>%
    filter(clade_group == clade_label)

  genome_set <- metadata_clade$genome_id
  snv_dat <- read_snv_long(snv_path) %>%
    filter(id1 %in% genome_set, id2 %in% genome_set)

  vertices <- sort(unique(genome_set))
  edges <- snv_dat %>%
    filter(snvs < snv_cutoff) %>%
    select(id1, id2) %>%
    distinct()

  g <- make_empty_graph(n = length(vertices), directed = FALSE) %>%
    set_vertex_attr("name", value = vertices)

  if (nrow(edges) > 0) {
    g <- add_edges(g, as.vector(t(as.matrix(edges))))
  }

  comps <- components(g)

  cluster_membership <- tibble(
    genome_id = names(comps$membership),
    cluster_id = as.integer(comps$membership)
  ) %>%
    left_join(metadata_clade, by = "genome_id") %>%
    arrange(cluster_id, genome_id)

  cluster_patient_counts_all <- cluster_membership %>%
    group_by(cluster_id) %>%
    summarise(
      cluster_size_all = n(),
      patient_n_all = n_distinct(patient_id[!is.na(patient_id) & patient_id != ""]),
      .groups = "drop"
    )

  retained_cluster_ids <- cluster_patient_counts_all %>%
    filter(patient_n_all > 1) %>%
    pull(cluster_id)

  cluster_membership <- cluster_membership %>%
    left_join(cluster_patient_counts_all, by = "cluster_id") %>%
    filter(cluster_id %in% retained_cluster_ids) %>%
    arrange(cluster_id, genome_id)

  cluster_sizes <- cluster_membership %>%
    count(cluster_id, name = "cluster_size") %>%
    left_join(
      cluster_membership %>%
        group_by(cluster_id) %>%
        summarise(patient_n = n_distinct(patient_id[!is.na(patient_id) & patient_id != ""]), .groups = "drop"),
      by = "cluster_id"
    ) %>%
    arrange(desc(cluster_size), cluster_id)

  cluster_pair_ranges <- snv_dat %>%
    left_join(cluster_membership %>% select(genome_id, cluster_id), by = c("id1" = "genome_id")) %>%
    rename(cluster_1 = cluster_id) %>%
    left_join(cluster_membership %>% select(genome_id, cluster_id), by = c("id2" = "genome_id")) %>%
    rename(cluster_2 = cluster_id) %>%
    filter(!is.na(cluster_1), !is.na(cluster_2), cluster_1 == cluster_2) %>%
    mutate(cluster_id = cluster_1) %>%
    group_by(cluster_id) %>%
    summarise(
      n_pairs = n(),
      min_snv = min(snvs, na.rm = TRUE),
      max_snv = max(snvs, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    right_join(cluster_sizes, by = "cluster_id") %>%
    mutate(
      n_pairs = replace_na(n_pairs, 0),
      min_snv = replace_na(min_snv, 0),
      max_snv = replace_na(max_snv, 0)
    ) %>%
    arrange(desc(cluster_size), cluster_id)

  cluster_summary <- if (nrow(cluster_membership) > 0) {
    tibble(
      clade_group = clade_label,
      snv_rule = "SNV < 2",
      filter_rule = "Cross-patient clusters only",
      n_genomes = nrow(cluster_membership),
      n_clusters = n_distinct(cluster_membership$cluster_id),
      n_edges = nrow(edges %>% filter(id1 %in% cluster_membership$genome_id, id2 %in% cluster_membership$genome_id)),
      n_singletons = sum(cluster_sizes$cluster_size == 1),
      min_cluster_size = min(cluster_sizes$cluster_size),
      median_cluster_size = median(cluster_sizes$cluster_size),
      max_cluster_size = max(cluster_sizes$cluster_size),
      n_clusters_removed_same_patient_only = sum(cluster_patient_counts_all$patient_n_all <= 1),
      n_genomes_removed_same_patient_only = sum(cluster_patient_counts_all$cluster_size_all[cluster_patient_counts_all$patient_n_all <= 1])
    )
  } else {
    tibble(
      clade_group = clade_label,
      snv_rule = "SNV < 2",
      filter_rule = "Cross-patient clusters only",
      n_genomes = 0,
      n_clusters = 0,
      n_edges = 0,
      n_singletons = 0,
      min_cluster_size = NA_real_,
      median_cluster_size = NA_real_,
      max_cluster_size = NA_real_,
      n_clusters_removed_same_patient_only = sum(cluster_patient_counts_all$patient_n_all <= 1),
      n_genomes_removed_same_patient_only = sum(cluster_patient_counts_all$cluster_size_all[cluster_patient_counts_all$patient_n_all <= 1])
    )
  }

  st_palette <- get_saved_st_palette(clade_label, cluster_membership$ST)
  non_singleton_clusters <- cluster_sizes %>%
    filter(cluster_size > 1) %>%
    pull(cluster_id)

  vertex_df <- cluster_membership %>%
    filter(cluster_id %in% non_singleton_clusters) %>%
    mutate(
      degree = degree(g, v = genome_id),
      vertex_color = unname(st_palette[ST]),
      vertex_size = case_when(
        clade_label == "Clade 1" & degree == 0 ~ 1.2,
        clade_label == "Clade 1" ~ 1.7,
        degree == 0 ~ 3,
        TRUE ~ 4.5
      )
    )

  plot_edges <- edges %>%
    filter(id1 %in% vertex_df$genome_id, id2 %in% vertex_df$genome_id)

  g_plot <- graph_from_data_frame(
    d = plot_edges,
    directed = FALSE,
    vertices = vertex_df %>% transmute(name = genome_id, ST = ST)
  )

  if (vcount(g_plot) > 0) {
    set.seed(1001)
    lay <- layout_with_fr(g_plot, niter = 5000)
  } else {
    lay <- matrix(numeric(0), ncol = 2)
  }

  clade_stub <- tolower(gsub(" ", "_", clade_label))

  write.csv(
    cluster_membership,
    file.path(desc_dir, paste0(clade_stub, "_post_noncore_clusters_snv_lt2.csv")),
    row.names = FALSE
  )
  write.csv(
    cluster_sizes,
    file.path(desc_dir, paste0(clade_stub, "_post_noncore_cluster_sizes_snv_lt2.csv")),
    row.names = FALSE
  )
  write.csv(
    cluster_pair_ranges,
    file.path(desc_dir, paste0(clade_stub, "_post_noncore_cluster_snv_ranges_snv_lt2.csv")),
    row.names = FALSE
  )
  write.csv(
    cluster_summary,
    file.path(desc_dir, paste0(clade_stub, "_post_noncore_cluster_summary_snv_lt2.csv")),
    row.names = FALSE
  )

  write.csv(
    tibble(ST = names(st_palette), color = unname(st_palette)),
    file.path(desc_dir, paste0(clade_stub, "_post_noncore_adjacency_st_color_key.csv")),
    row.names = FALSE
  )

  png_path <- file.path(viz_dir, paste0(clade_stub, "_post_noncore_adjacency_graph_snv_lt2.png"))
  pdf_path <- file.path(viz_dir, paste0(clade_stub, "_post_noncore_adjacency_graph_snv_lt2.pdf"))
  adjacency_width <- ifelse(clade_label == "Clade 2", 8, 12)
  adjacency_height <- ifelse(clade_label == "Clade 2", 8, 12)

  if (vcount(g_plot) > 0) {
    vertex_plot_df <- vertex_df %>%
      mutate(
        x = lay[, 1],
        y = lay[, 2]
      )

    edge_plot_df <- as_data_frame(g_plot, what = "edges") %>%
      left_join(
        vertex_plot_df %>% select(genome_id, x, y) %>% rename(x = x, y = y),
        by = c("from" = "genome_id")
      ) %>%
      left_join(
        vertex_plot_df %>% select(genome_id, x, y) %>% rename(xend = x, yend = y),
        by = c("to" = "genome_id")
      )

    adjacency_plot <- ggplot() +
      geom_segment(
        data = edge_plot_df,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = "#b0b7c3",
        linewidth = ifelse(clade_label == "Clade 2", 0.6, 0.25),
        alpha = 0.7
      ) +
      geom_point(
        data = vertex_plot_df,
        aes(x = x, y = y, color = ST, size = vertex_size),
        alpha = 0.95
      ) +
      scale_color_manual(
        values = st_palette,
        breaks = names(st_palette),
        labels = stringr::str_remove(names(st_palette), "^ST\\s+"),
        name = "ST"
      ) +
      scale_size_identity() +
      labs(
      title = paste0(clade_label, ": post-noncore adjacency graph (SNV < 2, non-singletons only)"),
        subtitle = "Node color = ST; same-patient-only clusters removed"
      ) +
      theme_void(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 8)
      ) +
      guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))
  } else {
    adjacency_plot <- ggplot() +
      theme_void() +
      labs(title = paste0(clade_label, ": no non-singleton clusters at SNV < 2"))
  }

  ggsave(
    filename = png_path,
    plot = adjacency_plot,
    width = adjacency_width,
    height = adjacency_height,
    units = "in",
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
  save_plot_pdf(
    plot_obj = adjacency_plot,
    filename = pdf_path,
    width = adjacency_width,
    height = adjacency_height
  )

  bar_plot <- cluster_membership %>%
    count(cluster_id, ST, name = "n") %>%
    left_join(cluster_sizes, by = "cluster_id") %>%
    mutate(cluster_id = factor(cluster_id, levels = cluster_sizes$cluster_id)) %>%
    ggplot(aes(x = cluster_id, y = n, fill = ST)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = st_palette, name = "ST") +
    labs(
      title = paste0(clade_label, ": post-noncore clusters (SNV < 2)"),
      subtitle = "Cross-patient clusters only",
      x = "Cluster",
      y = "Number of genomes"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(
    filename = file.path(viz_dir, paste0(clade_stub, "_post_noncore_cluster_barplot_snv_lt2.png")),
    plot = bar_plot,
    width = 12,
    height = 6,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  save_plot_pdf(
    plot_obj = bar_plot,
    filename = file.path(viz_dir, paste0(clade_stub, "_post_noncore_cluster_barplot_snv_lt2.pdf")),
    width = 12,
    height = 6
  )
}

clade_inputs <- tibble(
  clade_label = c("Clade 1", "Clade 2"),
  snv_path = c(clade1_snv_path, clade2_snv_path)
)

purrr::pwalk(
  clade_inputs,
  ~ build_cluster_outputs(..1, ..2)
)

cat("Done. Post-noncore cluster outputs written to:\n", desc_dir, "\n", sep = "")
