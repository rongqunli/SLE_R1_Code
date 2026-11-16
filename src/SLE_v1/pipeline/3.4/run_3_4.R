if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
analysis_root <- if (length(file_arg)) {
  normalizePath(dirname(file_arg[[1]]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), "pipeline", "3.4"), winslash = "/", mustWork = TRUE)
}
project_root <- normalizePath(file.path(analysis_root, "..", ".."), winslash = "/", mustWork = TRUE)
results_root <- file.path(project_root, "results", "3.4")
for (directory in c("00_input_audit", "01_overall", "02_group_stratified", "03_cross_trial", "04_summary", "figures", "logs")) {
  dir.create(file.path(results_root, directory), recursive = TRUE, showWarnings = FALSE)
}

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
scripts <- c(
  "00_validate_inputs.R",
  "01_overall_pharmacodynamics.R",
  "02_group_stratified_pharmacodynamics.R",
  "03_cross_trial_synthesis.R",
  "04_summarize_results.R"
)
for (script in scripts) {
  cat("\nRunning", script, "at", format(Sys.time()), "\n")
  source(file.path(analysis_root, "scripts", script), local = analysis_env)
}

saveRDS(
  list(
    study = "3.4_cross_trial_transcriptional_pharmacodynamics",
    status = "completed_with_pharmacodynamic_boundaries",
    completed_at = format(Sys.time(), tz = "Asia/Shanghai"),
    required_outputs = c(
      "00_input_audit/input_manifest.tsv",
      "01_overall/ILLUMINATE-1/mixed_model_results.tsv",
      "01_overall/ILLUMINATE-2/mixed_model_results.tsv",
      "01_overall/paired_change_summary.tsv",
      "02_group_stratified/three_way_interaction_tests.tsv",
      "02_group_stratified/subtype_specific_difference_in_differences.tsv",
      "03_cross_trial/cross_trial_fixed_effect_meta.tsv",
      "04_summary/key_findings.tsv"
    )
  ),
  file.path(results_root, "_SUCCESS.rds")
)
cat("\n3.4 pipeline completed at", format(Sys.time()), "\n")

