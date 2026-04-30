# Step-by-Step Code Path

This folder is the small, reviewer-friendly code layer for the project.

Each file represents one stage of the workflow:

1. `01_qc_and_assembly_review.sh`
2. `02_variant_calling_with_snpkit.sh`
3. `03_core_gene_alignment_with_cognac.R`
4. `04_patient_metadata_overview.R`
5. `05_alignment_filtering_and_qc.R`
6. `06_compare_distance_definitions.R`
7. `07_build_case_control_pairs.R`
8. `08_cluster_overlap_analysis.R`
9. `09_submit_snv_pipeline.slurm`

These are intentionally smaller than the full project scripts. They show the logic of one step at a time.

If you want the fuller curated code package, use `scripts/`.

If you want the full original working tree, use `source_code/`.
