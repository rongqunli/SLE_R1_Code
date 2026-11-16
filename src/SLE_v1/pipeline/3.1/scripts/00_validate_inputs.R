if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
source(file.path(dirname(normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]]
), winslash = "/", mustWork = TRUE)), "_common.R"))

require_packages(c("data.table", "digest"), "GSE65391 input audit")

script_path <- script_file()
expression_file <- file.path(cfg$paths$raw, cfg$files$expression)
platform_file <- file.path(cfg$paths$raw, cfg$files$platform)
manifest_file <- file.path(cfg$paths$metadata, "download_manifest.tsv")
output_dir <- file.path(cfg$paths$results, "00_input_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(expression_file, platform_file, manifest_file, cfg$paths$frozen_meta)
if (!all(file.exists(required_inputs))) {
  stop("Missing required input(s): ", paste(required_inputs[!file.exists(required_inputs)], collapse = ", "))
}

outputs <- file.path(output_dir, c(
  "input_audit_summary.tsv", "sample_metadata_audit.tsv",
  "clinical_field_audit.tsv", "patient_visit_summary.tsv"
))
signature <- build_step_signature(
  input_files = required_inputs,
  code_files = c(script_path, file.path(PIPELINE_ROOT, "scripts", "_common.R"),
                 file.path(PIPELINE_ROOT, "config", "config.R"))
)
if (!has_flag("--force") && step_is_complete(output_dir, signature, outputs)) {
  log_message("Step 00 input audit is already complete.")
  quit(save = "no", status = 0L)
}

read_series_sample_ids <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection))
  sample_ids <- NULL
  repeat {
    line <- readLines(connection, n = 1L, warn = FALSE)
    if (!length(line)) break
    if (startsWith(line, "!Sample_geo_accession")) {
      fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
      sample_ids <- gsub('^"|"$', "", fields[-1L])
    }
    if (identical(line, "!series_matrix_table_begin")) break
  }
  if (is.null(sample_ids) || !length(sample_ids)) {
    stop("Could not find !Sample_geo_accession in series matrix header.")
  }
  sample_ids
}

matrix_sample_ids <- read_series_sample_ids(expression_file)
if (anyDuplicated(matrix_sample_ids)) stop("Duplicate sample IDs in expression matrix header.")

meta <- data.table::fread(cfg$paths$frozen_meta, sep = "\t", na.strings = NULL)
if (anyDuplicated(meta$sample_id)) stop("Duplicate sample IDs in frozen GSE65391 metadata.")

sample_audit <- data.table::copy(meta)
sample_audit[, in_expression := sample_id %in% matrix_sample_ids]
sample_audit[, in_frozen_metadata := TRUE]
expression_only <- setdiff(matrix_sample_ids, meta$sample_id)
if (length(expression_only)) {
  sample_audit <- data.table::rbindlist(list(
    sample_audit,
    data.table::data.table(
      sample_id = expression_only,
      in_expression = TRUE,
      in_frozen_metadata = FALSE
    )
  ), use.names = TRUE, fill = TRUE)
}

included <- meta[analysis_include == "yes"]
included_sle <- included[disease_state == "SLE"]
malformed <- included[!grepl("^V[0-9]+$", visit_id)]

missing_codes <- c("not_available", "not_applicable", "Data Not Available", "unknown", "")
clinical_fields <- c(
  "timepoint", "SLEDAI", "renal_activity", "C3", "C4", "anti_dsDNA",
  "treatment", "age", "sex", "batch"
)
clinical_audit <- data.table::rbindlist(lapply(clinical_fields, function(field) {
  values <- as.character(included_sle[[field]])
  available <- !values %in% missing_codes
  numeric_values <- suppressWarnings(as.numeric(values[available]))
  data.table::data.table(
    field = field,
    total_sle_visits = length(values),
    available_values = sum(available),
    explicit_missing = sum(!available),
    numeric_values = sum(is.finite(numeric_values)),
    non_numeric_available = sum(available) - sum(is.finite(numeric_values)),
    unique_available_values = data.table::uniqueN(values[available])
  )
}))

patient_visits <- included[, .(
  disease_state = paste(sort(unique(disease_state)), collapse = ";"),
  n_visits = .N,
  first_timepoint = suppressWarnings(min(as.numeric(timepoint), na.rm = TRUE)),
  last_timepoint = suppressWarnings(max(as.numeric(timepoint), na.rm = TRUE)),
  has_malformed_visit_id = any(!grepl("^V[0-9]+$", visit_id))
), by = patient_id]

summary <- data.table::data.table(
  metric = c(
    "matrix_arrays", "metadata_rows", "included_arrays", "technical_replicates",
    "included_sle_visits", "sle_patients", "healthy_arrays",
    "malformed_visit_ids", "expression_only_samples", "metadata_only_samples"
  ),
  value = c(
    length(matrix_sample_ids), nrow(meta), nrow(included),
    nrow(meta[analysis_include == "no"]), nrow(included_sle),
    data.table::uniqueN(included_sle$patient_id), nrow(included[disease_state == "Healthy"]),
    nrow(malformed), length(expression_only), sum(!sample_audit$in_expression, na.rm = TRUE)
  )
)

expected <- unlist(cfg$expected[c(
  "total_arrays", "included_arrays", "technical_replicates", "sle_visits",
  "sle_patients", "healthy_arrays", "malformed_visit_ids"
)])
observed <- c(
  total_arrays = length(matrix_sample_ids), included_arrays = nrow(included),
  technical_replicates = nrow(meta[analysis_include == "no"]),
  sle_visits = nrow(included_sle), sle_patients = data.table::uniqueN(included_sle$patient_id),
  healthy_arrays = nrow(included[disease_state == "Healthy"]), malformed_visit_ids = nrow(malformed)
)
if (!identical(as.numeric(observed), as.numeric(expected))) {
  stop("Observed input counts do not match frozen config.")
}
if (length(expression_only) || any(!sample_audit$in_expression, na.rm = TRUE)) {
  stop("Expression and frozen metadata sample IDs do not match one-to-one.")
}

snapshot_path <- file.path(
  cfg$paths$metadata,
  paste0("GSE65391_meta_", cfg$meta_version, ".tsv")
)
if (!file.copy(cfg$paths$frozen_meta, snapshot_path, overwrite = TRUE)) {
  stop("Failed to copy frozen metadata snapshot.")
}
snapshot_hash <- digest::digest(snapshot_path, algo = "sha256", file = TRUE, serialize = FALSE)
safe_write_tsv(
  data.table::data.table(file = basename(snapshot_path), sha256 = snapshot_hash),
  file.path(cfg$paths$metadata, "metadata_snapshot_manifest.tsv")
)

safe_write_tsv(summary, outputs[[1]])
safe_write_tsv(sample_audit, outputs[[2]])
safe_write_tsv(clinical_audit, outputs[[3]])
safe_write_tsv(patient_visits, outputs[[4]])
write_session_info(file.path(cfg$paths$logs, "00_input_audit_sessionInfo.txt"))
mark_step_complete(output_dir, signature, outputs)
log_message("Input audit complete: ", nrow(included), " included arrays from ",
            data.table::uniqueN(included$patient_id), " subjects.")
