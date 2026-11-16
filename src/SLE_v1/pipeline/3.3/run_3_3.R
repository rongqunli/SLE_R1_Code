if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
analysis_root <- if (length(file_arg)) {
  normalizePath(dirname(file_arg[[1]]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), "pipeline", "3.3"), winslash = "/", mustWork = TRUE)
}
project_root <- normalizePath(file.path(analysis_root, "..", ".."), winslash = "/", mustWork = TRUE)
results_root <- file.path(project_root, "results", "3.3")
dir.create(file.path(results_root, "logs"), recursive = TRUE, showWarnings = FALSE)

options(stringsAsFactors = FALSE, width = 140)
Sys.setenv(SLE_PROJECT_ROOT = project_root)
log_file <- file.path(results_root, "logs", "pipeline.log")
sink(log_file, append = FALSE, split = TRUE)
on.exit({
  cat("\nSession information\n")
  print(sessionInfo())
  sink()
}, add = TRUE)

analysis_env <- new.env(parent = globalenv())
for (script in c("00_validate_inputs.R", "01_consensus_grouping.R", "02_summarize_results.R")) {
  cat("\nRunning", script, "at", format(Sys.time()), "\n")
  source(file.path(analysis_root, "scripts", script), local = analysis_env)
}

saveRDS(
  list(
    study = "3.3_program_based_patient_groupings",
    status = "completed_with_interpretive_boundaries",
    completed_at = format(Sys.time(), tz = "Asia/Shanghai"),
    required_outputs = c(
      "00_input_audit/input_manifest.tsv",
      "01_grouping/cluster_number_stability.tsv",
      "01_grouping/baseline_subtype_assignments.tsv",
      "01_grouping/frozen_subtype_centroids.tsv",
      "01_grouping/validation_assignment_audit.tsv",
      "01_grouping/trial_baseline_overlap_audit.tsv",
      "01_grouping/distance_vs_correlation_assignment_audit.tsv",
      "02_summary/current_key_findings.tsv"
    )
  ),
  file.path(results_root, "_SUCCESS.rds")
)
cat("\n3.3 pipeline completed at", format(Sys.time()), "\n")

