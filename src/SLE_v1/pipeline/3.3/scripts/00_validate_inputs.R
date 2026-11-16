if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
required <- c("data.table", "ggplot2", "cluster", "ComplexHeatmap", "circlize", "ggalluvial")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "))

project_root <- Sys.getenv("SLE_PROJECT_ROOT")
stopifnot(nzchar(project_root))
analysis_root <- file.path(project_root, "pipeline", "3.3")
results_root <- file.path(project_root, "results", "3.3")
dir.create(file.path(results_root, "00_input_audit"), recursive = TRUE, showWarnings = FALSE)

contract <- data.table::fread(file.path(analysis_root, "config", "analysis_contract.tsv"))
cfg <- setNames(as.list(contract$value), contract$key)
resolve <- function(path) normalizePath(file.path(analysis_root, path), winslash = "/", mustWork = TRUE)

project_norm <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
assert_internal <- function(path) {
  path_norm <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (!startsWith(path_norm, paste0(project_norm, "/")) && path_norm != project_norm) {
    stop("Input resolves outside SLE_v1: ", path_norm)
  }
  path_norm
}

input_keys <- c("gse88884_object", "gse88885_object", "gse88886_object", "frozen_metadata")
manifest <- data.table::rbindlist(lapply(input_keys, function(key) {
  path <- assert_internal(resolve(cfg[[key]]))
  info <- file.info(path)
  data.table::data.table(
    input_key = key,
    path = path,
    exists = file.exists(path),
    size_bytes = as.numeric(info$size),
    md5 = as.character(tools::md5sum(path))
  )
}))
if (any(!manifest$exists)) stop("Missing 3.3 input file")
data.table::fwrite(manifest, file.path(results_root, "00_input_audit", "input_manifest.tsv"), sep = "\t")

contract_audit <- data.table::data.table(
  contract_key = contract$key,
  value = contract$value,
  is_absolute = grepl("^[A-Za-z]:[/\\\\]", contract$value),
  resolved_inside_project = vapply(seq_len(nrow(contract)), function(i) {
    if (!(contract$key[[i]] %in% input_keys)) return(NA)
    startsWith(assert_internal(resolve(contract$value[[i]])), paste0(project_norm, "/"))
  }, logical(1))
)
data.table::fwrite(contract_audit, file.path(results_root, "00_input_audit", "contract_audit.tsv"), sep = "\t")
if (any(contract_audit$is_absolute)) stop("Absolute path found in the 3.3 contract")

object_audit <- data.table::rbindlist(lapply(input_keys[1:3], function(key) {
  obj <- readRDS(resolve(cfg[[key]]))
  if (!all(c("metadata", "program_scores") %in% names(obj))) stop("Invalid object: ", key)
  required_score <- c("sample_id", "subject_id", "condition", "time", "program", "score_z", "analysis_role", "expression_column")
  if (!all(required_score %in% names(obj$program_scores))) stop("Program-score contract failed: ", key)
  data.table::data.table(
    input_key = key,
    metadata_rows = nrow(obj$metadata),
    program_score_rows = nrow(obj$program_scores),
    programs = data.table::uniqueN(obj$program_scores$program),
    subjects = data.table::uniqueN(obj$metadata$subject_id)
  )
}))
data.table::fwrite(object_audit, file.path(results_root, "00_input_audit", "object_structure_audit.tsv"), sep = "\t")
cat("3.3 input audit passed.\n")

