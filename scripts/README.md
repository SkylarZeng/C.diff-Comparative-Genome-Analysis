# Scripts

This folder is the main curated code package for the project, structured in a way that matches the metagenomics portfolio style.

## Subfolders

- `bash/`: upstream QC and variant-calling command workflows
- `slurm/`: HPC submission entrypoints
- `r/`: main R notebooks, pipeline code, and downstream analysis scripts
- `r/descriptive_analysis_presentation/`: stepwise cluster-visualization and epidemiology scripts
- `qmd/`: focused Quarto comparison notebooks
- `lib/`: helper utilities used by the analysis

## Suggested Review Path

1. `bash/01_qc_and_assembly_review.sh`
2. `bash/02_variant_calling_with_snpkit.sh`
3. `r/cdiff_cognac_070125.R`
4. `r/SNV_Pipeline_design.Rmd`
5. `r/run_snv_pipeline.R`
6. `qmd/compare_core_vs_soft-core_non_core_distances_ST1_clade2_R20291.qmd`
7. `r/build_pair_table.Rmd`
8. `r/descriptive_analysis_presentation/`

For the full original working tree, including many additional exploratory files, see `../source_code/`.
