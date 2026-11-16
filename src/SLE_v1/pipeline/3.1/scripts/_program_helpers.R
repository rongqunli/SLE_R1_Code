if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
matched_program_genes <- function(matching, program_name) {
  matching$gene[matching$program == program_name & matching$matched]
}

