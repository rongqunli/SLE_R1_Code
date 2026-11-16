if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
source(file.path(dirname(normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]]
), winslash = "/", mustWork = TRUE)), "_common.R"))
source(file.path(PIPELINE_ROOT, "scripts", "_helpers.R"), local = TRUE)
require_packages(c("data.table", "digest"), "GSE121239 input audit")
suppressPackageStartupMessages(library(data.table))


output_dir <- file.path(cfg$paths$results, "00_input_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- data.table::as.data.table(cfg$downloads)[file %in% c(cfg$files$series_matrix, cfg$files$platform)]
manifest[, path := file.path(cfg$paths$downloads, file)]
manifest[, observed_bytes := as.numeric(file.info(path)$size)]
manifest[, observed_sha256 := toupper(vapply(path, digest::digest, character(1), algo = 'sha256', file = TRUE, serialize = FALSE))]
stopifnot(nrow(manifest) == 2L, all(manifest$observed_bytes == manifest$expected_bytes), all(manifest$observed_sha256 == manifest$expected_sha256))
safe_write_tsv(manifest, file.path(output_dir, 'download_manifest.tsv'))
matrix_path <- file.path(cfg$paths$downloads, cfg$files$series_matrix)
matrix_ids <- read_matrix_sample_ids(matrix_path)
meta <- data.table::fread(cfg$paths$frozen_meta, sep = "\t", na.strings = NULL)
raw_meta <- data.table::fread(cfg$paths$raw_geo_meta, sep = "\t", na.strings = NULL)
stopifnot(!anyDuplicated(meta$sample_id), !anyDuplicated(raw_meta$sample_id),
          setequal(meta$sample_id, raw_meta$sample_id), setequal(meta$sample_id, matrix_ids))
raw_fields <- raw_meta[, .(
  sample_id,
  raw_disease = characteristic_disease_state,
  raw_patient_id = characteristic_patient_id,
  raw_sledai = characteristic_sledai,
  visit_date_raw = characteristic_visit_date,
  tissue_raw = characteristic_tissue,
  imputed_neutrophil_percentage = characteristic_imputed_neutrophil_percentage,
  source_description = sample_description
)]
analysis_meta <- merge(meta, raw_fields, by = "sample_id", all.x = TRUE, sort = FALSE)
analysis_meta <- analysis_meta[match(meta$sample_id, sample_id)]
analysis_meta[
  visit_date_raw %in% c("", "NA", "not_available", "not_applicable"),
  visit_date_raw := NA_character_
]
analysis_meta[, `:=`(
  SLEDAI_numeric = parse_numeric(SLEDAI),
  neutrophil_numeric = parse_numeric(imputed_neutrophil_percentage),
  visit_date = data.table::as.IDate(visit_date_raw, format = "%Y-%m-%d")
)]
analysis_meta[disease == "Healthy control", visit_date := data.table::as.IDate(NA)]
analysis_meta[disease == "SLE", years_since_first := as.numeric(
  visit_date - min(visit_date, na.rm = TRUE)
) / 365.25, by = patient_id]
analysis_meta[disease == "Healthy control", years_since_first := 0]
analysis_meta[, visit_order := data.table::frank(
  visit_date, ties.method = "first", na.last = "keep"
), by = patient_id]
analysis_meta[disease == "Healthy control", visit_order := 1]

sle <- analysis_meta[disease == "SLE"]
if (any(!is.finite(sle$SLEDAI_numeric)) || anyNA(sle$visit_date)) {
  stop("All included SLE visits require numeric SLEDAI and a visit date.")
}
variation <- sle[, .(
  n_visits = .N,
  distinct_sledai = data.table::uniqueN(SLEDAI_numeric),
  sledai_min = min(SLEDAI_numeric),
  sledai_max = max(SLEDAI_numeric),
  neutrophil_visits = sum(is.finite(neutrophil_numeric))
), by = patient_id]


stopifnot(nrow(analysis_meta) == 312L, nrow(sle) == 292L, uniqueN(sle$patient_id) == 65L)
safe_write_tsv(analysis_meta, file.path(output_dir, "analysis_metadata.tsv"))
safe_write_tsv(meta, file.path(output_dir, "sample_audit.tsv"))
safe_write_tsv(variation, file.path(output_dir, "longitudinal_feasibility.tsv"))
