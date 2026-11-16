if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
wrapper_file <- sys.frame(1)$ofile
if (is.null(wrapper_file) || !nzchar(wrapper_file)) wrapper_file <- "scripts/_common.R"
shared <- normalizePath(
  file.path(dirname(wrapper_file), "..", "..", "GSE135779_pipeline", "scripts", "_common.R"),
  winslash = "/", mustWork = FALSE
)
if (!file.exists(shared)) {
  shared <- normalizePath(file.path(getwd(), "GSE135779_pipeline", "scripts", "_common.R"),
                          winslash = "/", mustWork = TRUE)
}
source(shared, local = FALSE)

