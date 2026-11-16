if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
analysis_root <- if (length(file_arg)) {
  normalizePath(dirname(file_arg[[1]]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), "pipeline", "3.2"), winslash = "/", mustWork = TRUE)
}
project_root <- normalizePath(file.path(analysis_root, "..", ".."), winslash = "/", mustWork = TRUE)
results_root <- file.path(project_root, "results", "3.2")
dir.create(file.path(results_root, "logs"), recursive = TRUE, showWarnings = FALSE)

options(stringsAsFactors = FALSE, width = 120)
Sys.setenv(SLE_PROJECT_ROOT = project_root)
log_file <- file.path(results_root, "logs", "pipeline.log")
sink(log_file, append = FALSE, split = TRUE)
on.exit({
  cat("\nSession information\n")
  print(sessionInfo())
  sink()
}, add = TRUE)

source(file.path(analysis_root, "scripts", "00_validate_inputs.R"), local = new.env(parent = globalenv()))
for (script in c(
  "01_prepare_full_expression.R",
  "02_limma_and_gsea.R",
  "03_module_preservation.R",
  "04_summarize_results.R"
)) {
  cat("\nRunning", script, "at", format(Sys.time()), "\n")
  source(file.path(analysis_root, "scripts", script), local = new.env(parent = globalenv()))
}

saveRDS(
  list(
    study = "3.2_adult_bulk_and_module_preservation",
    status = "completed_positive",
    completed_at = format(Sys.time(), tz = "Asia/Shanghai"),
    required_outputs = c(
      "00_input_audit/input_manifest.tsv",
      "01_expression/full_gene_expression.rds",
      "02_deg_gsea/full_transcriptome_limma.tsv",
      "02_deg_gsea/fgsea_all_collections.tsv",
      "03_preservation/module_preservation_statistics.tsv",
      "04_summary/current_key_findings.tsv"
    )
  ),
  file.path(results_root, "_SUCCESS.rds")
)
cat("\n3.2 pipeline completed at", format(Sys.time()), "\n")
