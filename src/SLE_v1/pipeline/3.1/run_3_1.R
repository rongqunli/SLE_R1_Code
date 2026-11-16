if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Run this entry point with Rscript.")
analysis_root <- normalizePath(dirname(sub("^--file=", "", file_arg[[1]])), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(analysis_root, "..", ".."), winslash = "/", mustWork = TRUE)
results_root <- file.path(project_root, "results", "3.1")
dir.create(file.path(results_root, "logs"), recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(results_root, "logs", "pipeline.log")
writeLines(c(
  paste0("3.1 pipeline started: ", format(Sys.time(), tz = "Asia/Shanghai")),
  paste0("Project root: ", project_root)
), log_file)

scripts <- c(
  "00_validate_inputs.R",
  "01_prepare_expression.R",
  "02_score_programs.R",
  "03_longitudinal_models.R",
  "04_visualize_results.R",
  "05_gene_associations.R",
  "06_wgcna.R",
  "07_module_selection_audit.R",
  "08_summarize_results.R",
  "09_make_wgcna_figures.R"
)

rscript <- file.path(R.home("bin"), "Rscript")
forward_args <- args
for (script in scripts) {
  script_path <- file.path(analysis_root, "scripts", script)
  if (!file.exists(script_path)) stop("Missing source script: ", script_path)
  write(paste0("Running ", script, " at ", format(Sys.time(), tz = "Asia/Shanghai")), log_file, append = TRUE)
  status <- system2(rscript, c(shQuote(script_path), forward_args))
  if (!identical(as.integer(status), 0L)) {
    stop("Pipeline step failed: ", script)
  }
}

saveRDS(
  list(
    study = "3.1_pediatric_longitudinal_wgcna_and_module_selection",
    status = "completed",
    completed_at = format(Sys.time(), tz = "Asia/Shanghai"),
    required_outputs = c(
      "00_input_audit/input_audit_summary.tsv",
      "03_longitudinal/longitudinal_model_results.tsv",
      "06_wgcna/module_trait_associations.tsv",
      "06_wgcna/module_stability.tsv",
      "07_module_selection_audit/Table_S2_module_selection_audit.tsv",
      "08_summary/current_key_findings.tsv",
      "figures/Figure1C_module_selection_heatmap.png",
      "figures/Figure1D_module_stability.png",
      "figures/Figure1E_representative_hubs.png",
      "figures/FigureS_module_trait_heatmap.png"
    )
  ),
  file.path(results_root, "_SUCCESS.rds")
)
write(paste0("3.1 pipeline completed at ", format(Sys.time(), tz = "Asia/Shanghai")), log_file, append = TRUE)
message("3.1 source pipeline completed through summary step.")
