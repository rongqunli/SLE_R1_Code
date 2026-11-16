if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
project_root <- Sys.getenv("SLE_PROJECT_ROOT")
if (!nzchar(project_root)) stop("SLE_PROJECT_ROOT is not set")
analysis_root <- file.path(project_root, "pipeline", "3.4")
results_root <- file.path(project_root, "results", "3.4")
audit_root <- file.path(results_root, "00_input_audit")
dir.create(audit_root, recursive = TRUE, showWarnings = FALSE)

contract_boot <- utils::read.delim(
  file.path(analysis_root, "config", "analysis_contract.tsv"),
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
cfg_boot <- setNames(as.list(contract_boot$value), contract_boot$key)
local_r_library <- file.path(analysis_root, cfg_boot$local_r_library)
if (!dir.exists(local_r_library)) stop("Missing internal 3.4 R library: ", local_r_library)
.libPaths(unique(c(local_r_library, .libPaths())))
required_packages <- c("data.table", "lme4", "emmeans", "ggplot2", "ggalluvial")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required packages: ", paste(missing_packages, collapse = ", "))

contract <- data.table::fread(file.path(analysis_root, "config", "analysis_contract.tsv"))
cfg <- setNames(as.list(contract$value), contract$key)
resolve_input <- function(path) normalizePath(file.path(analysis_root, path), winslash = "/", mustWork = TRUE)

project_norm <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
assert_internal <- function(path) {
  path_norm <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (!startsWith(path_norm, paste0(project_norm, "/")) && path_norm != project_norm) {
    stop("Input resolves outside SLE_v1: ", path_norm)
  }
  path_norm
}

input_keys <- c("trial1_scores", "trial2_scores", "grouping_model")
path_keys <- c(input_keys, "local_r_library")
input_manifest <- data.table::rbindlist(lapply(input_keys, function(key) {
  path <- assert_internal(resolve_input(cfg[[key]]))
  info <- file.info(path)
  data.table::data.table(
    input_key = key,
    path = path,
    size_bytes = as.numeric(info$size),
    md5 = as.character(tools::md5sum(path))
  )
}))
data.table::fwrite(input_manifest, file.path(audit_root, "input_manifest.tsv"), sep = "\t")

contract_audit <- data.table::data.table(
  contract_key = contract$key,
  value = contract$value,
  is_absolute = grepl("^[A-Za-z]:[/\\\\]", contract$value),
  resolved_inside_project = vapply(seq_len(nrow(contract)), function(index) {
    if (!(contract$key[[index]] %in% path_keys)) return(NA)
    startsWith(assert_internal(resolve_input(contract$value[[index]])), paste0(project_norm, "/"))
  }, logical(1))
)
data.table::fwrite(contract_audit, file.path(audit_root, "contract_audit.tsv"), sep = "\t")
if (any(contract_audit$is_absolute)) stop("Absolute path found in the 3.4 contract")

read_trial_scores <- function(path, trial) {
  x <- data.table::fread(resolve_input(path), sep = "\t", na.strings = c("", "NA"))
  required_columns <- c(
    "program", "sample_id", "subject_id", "condition", "time", "treatment",
    "expression_column", "score_z", "analysis_role"
  )
  if (!all(required_columns %in% names(x))) stop(trial, " score table is missing required columns")
  if (anyDuplicated(x[, .(program, sample_id)])) stop(trial, " has duplicate program/sample rows")
  x[, trial := trial]
  x[]
}

trial_scores <- data.table::rbindlist(list(
  read_trial_scores(cfg$trial1_scores, "ILLUMINATE-1"),
  read_trial_scores(cfg$trial2_scores, "ILLUMINATE-2")
), use.names = TRUE, fill = TRUE)

primary_programs <- trimws(strsplit(cfg$primary_programs, ",", fixed = TRUE)[[1]])
time_levels <- trimws(strsplit(cfg$time_levels, ",", fixed = TRUE)[[1]])
if (!all(primary_programs %in% unique(trial_scores$program))) stop("A primary program is absent from trial inputs")

score_table_audit <- trial_scores[, .(
  table_rows = .N,
  samples = data.table::uniqueN(sample_id),
  subjects = data.table::uniqueN(subject_id),
  programs = data.table::uniqueN(program),
  primary_rows = sum(program %in% primary_programs & analysis_role == "primary"),
  finite_primary_rows = sum(program %in% primary_programs & analysis_role == "primary" & is.finite(score_z))
), by = trial]
data.table::fwrite(score_table_audit, file.path(audit_root, "score_table_audit.tsv"), sep = "\t")

sample_rows <- unique(trial_scores[, .(trial, sample_id, subject_id, time, treatment)])
visit_arm_audit <- sample_rows[, .(
  samples = .N,
  subjects = data.table::uniqueN(subject_id)
), by = .(trial, time, treatment)]
visit_arm_audit[, time_order := match(time, time_levels)]
data.table::setorder(visit_arm_audit, trial, time_order, treatment)
visit_arm_audit[, time_order := NULL]
data.table::fwrite(visit_arm_audit, file.path(audit_root, "visit_arm_counts.tsv"), sep = "\t")

grouping_model <- readRDS(resolve_input(cfg$grouping_model))
if (!("assignments" %in% names(grouping_model))) stop("Grouping model lacks assignments")
group_assignments <- data.table::as.data.table(grouping_model$assignments)
required_assignment_columns <- c(
  "subject_id", "cohort_partition", "subtype", "assignment_margin", "low_confidence"
)
if (!all(required_assignment_columns %in% names(group_assignments))) {
  stop("Grouping assignments do not satisfy the 3.4 contract")
}
if (anyDuplicated(group_assignments$subject_id)) stop("Grouping assignments contain duplicate subjects")

trial_subjects <- unique(trial_scores[, .(trial, subject_id)])
assignment_audit <- merge(
  trial_subjects,
  group_assignments[, .(subject_id, baseline_grouping = subtype, assignment_margin, low_confidence)],
  by = "subject_id",
  all.x = TRUE
)[, .(
  longitudinal_subjects = .N,
  assigned_subjects = sum(!is.na(baseline_grouping)),
  unassigned_subjects = sum(is.na(baseline_grouping)),
  low_margin_subjects = sum(low_confidence %in% TRUE, na.rm = TRUE),
  low_margin_fraction_among_assigned = mean(low_confidence %in% TRUE, na.rm = TRUE)
), by = trial]
data.table::fwrite(assignment_audit, file.path(audit_root, "grouping_input_audit.tsv"), sep = "\t")

cat("3.4 input and internal-path audit passed.\n")
