if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_arg[[1]])
PIPELINE_ROOT <- normalizePath(dirname(script_file), winslash = "/", mustWork = TRUE)

cfg <- list(
  accession = "GSE88885",
  role = "longitudinal",
  trial = "ILLUMINATE-1",
  expression_file = "GSE88885_ILLUMINATE1_pd_preprocessed.txt.gz",
  expected_samples = 908L,
  expected_rows = 92869L,
  expected_subjects = 309L
)

source(
  file.path(PIPELINE_ROOT, "..", "GSE88884_pipeline", "scripts", "_hta_pipeline.R"),
  local = FALSE
)
