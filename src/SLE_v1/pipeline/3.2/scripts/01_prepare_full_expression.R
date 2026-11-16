if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
project_root <- Sys.getenv("SLE_PROJECT_ROOT")
stopifnot(nzchar(project_root))
analysis_root <- file.path(project_root, "pipeline", "3.2")
results_root <- file.path(project_root, "results", "3.2")

required <- c("data.table", "AnnotationDbi")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "))

read_contract <- function(path) {
  x <- data.table::fread(path, sep = "\t", header = TRUE)
  setNames(as.list(x[["value"]]), x[["key"]])
}
resolve_path <- function(path) {
  if (grepl("^[A-Za-z]:[/\\\\]", path)) return(normalizePath(path, winslash = "/", mustWork = TRUE))
  normalizePath(file.path(analysis_root, path), winslash = "/", mustWork = TRUE)
}
cfg <- read_contract(file.path(analysis_root, "config", "analysis_contract.tsv"))
set.seed(as.integer(cfg$random_seed))

annotation_library <- resolve_path(cfg$annotation_library)
.libPaths(unique(c(annotation_library, .libPaths())))
if (!requireNamespace("hta20transcriptcluster.db", quietly = TRUE)) {
  stop("Local hta20transcriptcluster.db was not found after adding: ", annotation_library)
}

out_dir <- file.path(results_root, "01_expression")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw_file <- resolve_path(cfg$raw_expression)
audit <- data.table::fread(resolve_path(cfg$sample_audit), sep = "\t", na.strings = c("", "NA"))
clinical <- data.table::fread(
  resolve_path(cfg$frozen_metadata), sep = "\t",
  na.strings = c("", "NA", "not_available", "not_applicable", "--")
)

needed_audit <- c("sample_id", "subject_id", "condition", "time", "expression_column")
if (!all(needed_audit %in% names(audit))) stop("sample_audit is missing: ", paste(setdiff(needed_audit, names(audit)), collapse = ", "))
if (!all(c("sample_id", "patient_id", "SLEDAI", "C3", "C4", "anti_dsDNA", "age", "sex", "race", "batch") %in% names(clinical))) {
  stop("Frozen metadata does not satisfy the clinical-field contract.")
}
clinical_keep <- clinical[, .(
  sample_id, patient_id, SLEDAI, C3, C4, anti_dsDNA, age, sex, race, batch, cohort,
  analysis_include, analysis_role
)]
metadata <- merge(audit, clinical_keep, by = "sample_id", all.x = TRUE, sort = FALSE)
metadata <- metadata[match(audit$sample_id, sample_id)]

read_full_matrix <- function(path) {
  if (requireNamespace("readr", quietly = TRUE)) {
    message("Reading compressed matrix with readr.")
    tab <- readr::read_tsv(path, show_col_types = FALSE, progress = TRUE, name_repair = "minimal")
    return(data.table::as.data.table(tab))
  }
  message("readr is unavailable; using base read.delim. This uses more memory.")
  data.table::as.data.table(read.delim(gzfile(path), sep = "\t", quote = "\"", check.names = FALSE))
}

raw <- read_full_matrix(raw_file)
if (ncol(raw) < 3L) stop("Expression input has fewer than two samples.")
probe_col <- names(raw)[[1]]
raw_probe_id <- as.character(raw[[probe_col]])
transcript_cluster_id <- sub(":[0-9]+$", "", raw_probe_id)
raw[[probe_col]] <- NULL
sample_columns <- names(raw)

idx <- match(sample_columns, metadata$expression_column)
if (anyNA(idx)) {
  stop("Expression columns not mapped to sample audit: ", paste(head(sample_columns[is.na(idx)], 10), collapse = ", "))
}
metadata <- metadata[idx]
expr_probe <- as.matrix(raw)
storage.mode(expr_probe) <- "double"
rownames(expr_probe) <- raw_probe_id
colnames(expr_probe) <- metadata$sample_id
rm(raw)
invisible(gc())

annotation <- suppressMessages(AnnotationDbi::select(
  hta20transcriptcluster.db::hta20transcriptcluster.db,
  keys = unique(transcript_cluster_id[startsWith(transcript_cluster_id, "TC")]),
  keytype = "PROBEID",
  columns = c("SYMBOL", "ENTREZID")
))
annotation <- data.table::as.data.table(annotation)
annotation[, SYMBOL := toupper(trimws(SYMBOL))]
annotation <- annotation[!is.na(SYMBOL) & nzchar(SYMBOL)]
annotation[, ambiguous_probe := data.table::uniqueN(SYMBOL), by = PROBEID]
annotation <- annotation[ambiguous_probe == 1L]
annotation <- unique(annotation, by = c("PROBEID", "SYMBOL"))

probe_stats <- data.table::data.table(
  raw_probe_id = rownames(expr_probe),
  PROBEID = transcript_cluster_id,
  expression_iqr = apply(expr_probe, 1, stats::IQR, na.rm = TRUE),
  expression_mean = rowMeans(expr_probe, na.rm = TRUE)
)
mapping <- merge(annotation, probe_stats, by = "PROBEID", all.x = TRUE)
mapping <- mapping[is.finite(expression_iqr)]
data.table::setorder(mapping, SYMBOL, -expression_iqr, -expression_mean, PROBEID)
mapping[, selected_for_gene := seq_len(.N) == 1L, by = SYMBOL]
mapping[, passes_iqr := expression_iqr >= as.numeric(cfg$min_expression_iqr)]

selected <- mapping[selected_for_gene & passes_iqr]
if (nrow(selected) < 5000L) warning("Only ", nrow(selected), " genes passed annotation and IQR filtering.")
gene_expr <- expr_probe[match(selected$raw_probe_id, rownames(expr_probe)), , drop = FALSE]
rownames(gene_expr) <- selected$SYMBOL

audit_summary <- data.table::data.table(
  metric = c(
    "raw_rows", "raw_samples", "mapped_unambiguous_probes", "selected_genes",
    "sle_samples", "healthy_samples", "duplicated_sample_ids", "nonfinite_values"
  ),
  value = c(
    nrow(expr_probe), ncol(expr_probe), data.table::uniqueN(mapping$raw_probe_id), nrow(gene_expr),
    sum(metadata$condition == cfg$case_label), sum(metadata$condition == cfg$control_label),
    anyDuplicated(metadata$sample_id), sum(!is.finite(gene_expr))
  )
)
if (anyDuplicated(rownames(gene_expr)) || anyDuplicated(colnames(gene_expr))) stop("Gene or sample identifiers are not unique.")
if (any(!metadata$condition %in% c(cfg$case_label, cfg$control_label))) stop("Unexpected condition labels in GSE88884.")

data.table::fwrite(mapping, file.path(out_dir, "probe_to_gene_mapping_audit.tsv"), sep = "\t")
data.table::fwrite(metadata, file.path(out_dir, "analysis_metadata.tsv"), sep = "\t", na = "NA")
data.table::fwrite(audit_summary, file.path(out_dir, "expression_audit.tsv"), sep = "\t")
saveRDS(
  list(expression = gene_expr, metadata = metadata, mapping = selected, audit = audit_summary),
  file.path(out_dir, "full_gene_expression.rds"),
  compress = FALSE
)
