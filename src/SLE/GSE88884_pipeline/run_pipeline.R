if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_arg[[1]])
PIPELINE_ROOT <- normalizePath(dirname(script_file), winslash = "/", mustWork = TRUE)

cfg <- list(
  accession = "GSE88884",
  role = "baseline",
  expression_file = "GSE88884_ILLUMINATE1and2_SLEbaselineVsHealthy_preprocessed.txt.gz",
  expected_samples = 1820L,
  expected_rows = 92343L,
  expected_subjects = 1820L
)

source(file.path(PIPELINE_ROOT, "scripts", "_hta_pipeline.R"), local = FALSE)
