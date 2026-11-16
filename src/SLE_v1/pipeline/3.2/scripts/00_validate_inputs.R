if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
project_root <- Sys.getenv("SLE_PROJECT_ROOT")
stopifnot(nzchar(project_root))
analysis_root <- file.path(project_root, "pipeline", "3.2")
results_root <- file.path(project_root, "results", "3.2")
dir.create(file.path(results_root, "00_input_audit"), recursive = TRUE, showWarnings = FALSE)

contract <- data.table::fread(
  file.path(analysis_root, "config", "analysis_contract.tsv"),
  sep = "\t", header = TRUE, data.table = TRUE
)
cfg <- setNames(as.list(contract[["value"]]), contract[["key"]])

resolve <- function(path) {
  if (grepl("^[A-Za-z]:[/\\\\]", path)) normalizePath(path, winslash = "/", mustWork = TRUE)
  else normalizePath(file.path(analysis_root, path), winslash = "/", mustWork = TRUE)
}

project_norm <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
assert_internal <- function(path) {
  path_norm <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (!startsWith(path_norm, paste0(project_norm, "/")) && path_norm != project_norm) {
    stop("Input resolves outside SLE_v1: ", path_norm)
  }
  path_norm
}

input_keys <- c(
  "raw_expression", "sample_audit", "frozen_metadata", "annotation_library",
  "program_definitions", "go_gmt", "kegg_gmt", "reference_expression", "reference_wgcna"
)
records <- lapply(input_keys, function(key) {
  path <- assert_internal(resolve(cfg[[key]]))
  is_dir <- dir.exists(path)
  files <- if (is_dir) list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE) else path
  files <- files[file.info(files)$isdir %in% FALSE]
  sizes <- if (length(files)) file.info(files)$size else 0
  checksum_target <- if (is_dir) file.path(path, "DESCRIPTION") else path
  checksum <- if (file.exists(checksum_target)) as.character(tools::md5sum(checksum_target)) else NA_character_
  data.table::data.table(
    input_key = key,
    path = path,
    exists = file.exists(path) || dir.exists(path),
    type = if (is_dir) "directory" else "file",
    n_files = length(files),
    size_bytes = sum(sizes, na.rm = TRUE),
    md5_of_file_or_description = checksum
  )
})
manifest <- data.table::rbindlist(records)
if (any(!manifest$exists)) stop("Missing input(s): ", paste(manifest$input_key[!manifest$exists], collapse = ", "))
data.table::fwrite(manifest, file.path(results_root, "00_input_audit", "input_manifest.tsv"), sep = "\t")

resolved_flags <- vapply(seq_len(nrow(contract)), function(i) {
  if (!(contract[["key"]][[i]] %in% input_keys)) return(NA)
  tryCatch(
    startsWith(assert_internal(resolve(contract[["value"]][[i]])), paste0(project_norm, "/")),
    error = function(e) FALSE
  )
}, logical(1))
contract_audit <- data.table::data.table(
  contract_key = contract[["key"]],
  value = contract[["value"]],
  is_absolute = grepl("^[A-Za-z]:[/\\\\]", contract[["value"]]),
  resolved_inside_project = resolved_flags
)
data.table::fwrite(contract_audit, file.path(results_root, "00_input_audit", "contract_audit.tsv"), sep = "\t")
if (any(contract_audit$is_absolute)) stop("Absolute path found in the 3.2 contract.")
cat("3.2 input audit passed.\n")
