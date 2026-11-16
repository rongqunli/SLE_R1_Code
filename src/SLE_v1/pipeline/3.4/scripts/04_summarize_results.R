if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
summary_root <- file.path(results_root, "04_summary")
dir.create(summary_root, recursive = TRUE, showWarnings = FALSE)

meta <- data.table::fread(file.path(results_root, "03_cross_trial", "cross_trial_fixed_effect_meta.tsv"), sep = "\t")
group_test <- data.table::fread(file.path(results_root, "02_group_stratified", "three_way_interaction_tests.tsv"), sep = "\t")
direction <- data.table::fread(file.path(results_root, "02_group_stratified", "cross_trial_direction_consistency.tsv"), sep = "\t")
visit_counts <- data.table::fread(file.path(results_root, "00_input_audit", "visit_arm_counts.tsv"), sep = "\t")
score_audit <- data.table::fread(file.path(results_root, "00_input_audit", "score_table_audit.tsv"), sep = "\t")

bcell_meta <- meta[program == "GREENYELLOW_B_CELL_COMPLEMENT"]
key_findings <- data.table::rbindlist(list(
  data.table::data.table(
    finding = "overall_bcell_meta",
    scope = paste(bcell_meta$dose, paste0("week", bcell_meta$week), sep = "_"),
    estimate = bcell_meta$meta_estimate,
    p_value = bcell_meta$p_value,
    fdr = bcell_meta$fdr,
    interpretation = "Overall fixed-effect B-cell/complement pharmacodynamic decrease"
   ),
  group_test[program == "GREENYELLOW_B_CELL_COMPLEMENT", .(
    finding = "group_three_way",
    scope = trial,
    estimate = three_way_chisq,
    p_value = three_way_p_value,
    fdr = three_way_fdr,
    interpretation = "Treatment-by-time-by-group likelihood-ratio test"
   )],
  direction[program == "GREENYELLOW_B_CELL_COMPLEMENT", .(
    finding = "shared_group_direction",
    scope = paste(trial_1, trial_2, sep = "_"),
    estimate = same_direction_fraction,
    p_value = NA_real_,
    fdr = NA_real_,
    interpretation = "Descriptive fraction of eight shared contrasts with the same sign"
   )]
), fill = TRUE)
data.table::fwrite(key_findings, file.path(summary_root, "key_findings.tsv"), sep = "\t")

statistical_boundaries <- data.table::data.table(
  boundary = c(
    "analysis_unit", "missing_visits", "grouping_transport", "low_margin", "illuminate2_b_cell_high", "heterogeneity", "clinical_efficacy"
  ),
  statement = c(
    "Participant-level longitudinal mixed models; repeated samples are not independent observations.",
    "Overall models use all available visits under a likelihood-based missing-at-random assumption; paired sensitivity requires baseline and target visit.",
    "Baseline groupings are transported from frozen 3.3 centroids; they are not independently refit in the trials.",
    "Low-margin assignments are retained in the primary grouping models and treated as geometric uncertainty, not probabilities.",
    "B_CELL_HIGH is not estimable in ILLUMINATE-2 because the placebo arm has seven participants; this is not a negative replication.",
    "Cochran Q and I2 are descriptive with two trials and do not exclude between-trial heterogeneity.",
    "Deposited expression records do not contain harmonized individual clinical-response outcomes; findings indicate molecular target engagement, not clinical benefit or a treatment-selection rule."
  )
)
data.table::fwrite(statistical_boundaries, file.path(summary_root, "statistical_boundaries.tsv"), sep = "\t")

sample_summary <- data.table::rbindlist(list(
  score_audit[, .(trial, subjects, samples, primary_rows, finite_primary_rows)],
  visit_counts[, .(trial, time, treatment, subjects, samples)]
), fill = TRUE)
data.table::fwrite(sample_summary, file.path(summary_root, "sample_summary.tsv"), sep = "\t")

result_files <- list.files(results_root, recursive = TRUE, full.names = TRUE)
result_files <- result_files[file.info(result_files)$isdir %in% FALSE]
manifest <- data.table::data.table(
  relative_path = substring(result_files, nchar(results_root) + 2L),
  size_bytes = as.numeric(file.info(result_files)$size),
  md5 = as.character(tools::md5sum(result_files))
)
data.table::fwrite(manifest, file.path(summary_root, "current_result_manifest.tsv"), sep = "\t")
cat("3.4 summary tables and interpretive boundaries completed.\n")
