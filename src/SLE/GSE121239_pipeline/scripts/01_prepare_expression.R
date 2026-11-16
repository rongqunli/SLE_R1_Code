if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
source(file.path(dirname(normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]]
), winslash = "/", mustWork = TRUE)), "_common.R"))
source(file.path(PIPELINE_ROOT, "scripts", "_helpers.R"), local = TRUE)
require_packages(c("data.table", "matrixStats", "ggplot2"), "expression preparation")

script_path <- script_file()
matrix_file <- file.path(cfg$paths$downloads, cfg$files$series_matrix)
platform_file <- file.path(cfg$paths$downloads, cfg$files$platform)
analysis_meta_file <- file.path(cfg$paths$results, "00_input_audit", "analysis_metadata.tsv")
output_dir <- file.path(cfg$paths$results, "01_expression_qc")
object_file <- file.path(cfg$paths$processed, "01_gene_expression.rds")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
outputs <- c(object_file, file.path(output_dir, c(
  "expression_dimension_audit.tsv", "probe_gene_mapping.tsv",
  "probe_selection_audit.tsv", "sample_qc_metrics.tsv",
  "pca_coordinates.tsv", "expression_qc.pdf"
)))
signature <- build_step_signature(
  input_files = c(matrix_file, platform_file, analysis_meta_file),
  code_files = c(
    script_path, file.path(PIPELINE_ROOT, "scripts", "_common.R"),
    file.path(PIPELINE_ROOT, "scripts", "_helpers.R"),
    file.path(PIPELINE_ROOT, "config", "config.R")
  )
)
if (!has_flag("--force") && step_is_complete(output_dir, signature, outputs)) {
  log_message("Step 01 expression preparation is already complete.")
  quit(save = "no", status = 0L)
}

meta <- data.table::fread(analysis_meta_file, sep = "\t", na.strings = NULL)
log_message("Reading official GEO RMA-normalized series matrix.")
dt <- data.table::fread(matrix_file, skip = "ID_REF", quote = "\"", showProgress = TRUE)
data.table::setnames(dt, names(dt)[1L], "probe_id")
dt <- dt[!grepl("^!", probe_id)]
if (nrow(dt) != cfg$expected$probes || anyDuplicated(dt$probe_id)) {
  stop("Unexpected probe matrix structure.")
}
if (!setequal(names(dt)[-1L], meta$sample_id)) stop("Expression columns and metadata differ.")
data.table::setcolorder(dt, c("probe_id", meta$sample_id))
probe_ids <- dt$probe_id
expression <- as.matrix(dt[, -"probe_id"])
storage.mode(expression) <- "double"
rownames(expression) <- probe_ids
if (any(!is.finite(expression))) stop("Non-finite expression values detected.")
rm(dt); gc(FALSE)

mapping <- read_platform_mapping(platform_file)
safe_write_tsv(mapping, file.path(output_dir, "probe_gene_mapping.tsv"))
mapped <- merge(
  data.table::data.table(probe_id = rownames(expression), row = seq_len(nrow(expression))),
  mapping[mapping_status == "unique_symbol", .(probe_id, gene_symbol)],
  by = "probe_id"
)
mapped[, probe_variance := matrixStats::rowVars(expression[row, , drop = FALSE])]
data.table::setorder(mapped, gene_symbol, -probe_variance, probe_id)
mapped[, selected := seq_len(.N) == 1L, by = gene_symbol]
chosen <- mapped[selected == TRUE]
data.table::setorder(chosen, gene_symbol)
gene_expression <- expression[chosen$row, , drop = FALSE]
rownames(gene_expression) <- chosen$gene_symbol
colnames(gene_expression) <- meta$sample_id
safe_write_tsv(
  mapped[, .(probe_id, gene_symbol, probe_variance, selected)],
  file.path(output_dir, "probe_selection_audit.tsv")
)

top_rows <- order(
  matrixStats::rowVars(gene_expression), decreasing = TRUE
)[seq_len(min(3000L, nrow(gene_expression)))]
qc_matrix <- gene_expression[top_rows, , drop = FALSE]
pca <- stats::prcomp(t(qc_matrix), center = TRUE, scale. = FALSE, rank. = 10L)
correlation <- stats::cor(qc_matrix)
diag(correlation) <- NA_real_
median_cor <- matrixStats::colMedians(correlation, na.rm = TRUE)
pc_distance <- rowSums(scale(pca$x[, 1:5, drop = FALSE])^2)
cor_cut <- stats::median(median_cor) -
  cfg$analysis$outlier_mad_multiplier * stats::mad(median_cor)
pc_cut <- stats::median(pc_distance) +
  cfg$analysis$outlier_mad_multiplier * stats::mad(pc_distance)
sample_qc <- data.table::data.table(
  sample_id = colnames(gene_expression),
  expression_median = matrixStats::colMedians(gene_expression),
  expression_iqr = matrixStats::colIQRs(gene_expression),
  median_correlation = median_cor,
  robust_pc_distance = pc_distance,
  low_correlation_flag = median_cor < cor_cut,
  pca_outlier_flag = pc_distance > pc_cut
)
sample_qc[, outlier_flag := low_correlation_flag | pca_outlier_flag]
sample_qc <- merge(sample_qc, meta, by = "sample_id", sort = FALSE)
pca_dt <- data.table::as.data.table(pca$x, keep.rownames = "sample_id")
pca_dt <- merge(pca_dt, meta, by = "sample_id", sort = FALSE)
safe_write_tsv(sample_qc, file.path(output_dir, "sample_qc_metrics.tsv"))
safe_write_tsv(pca_dt, file.path(output_dir, "pca_coordinates.tsv"))
safe_write_tsv(data.table::data.table(
  metric = c("probes", "mapped_unique_probes", "selected_genes", "samples", "qc_flagged_samples"),
  value = c(
    nrow(expression), nrow(mapped), nrow(gene_expression),
    ncol(gene_expression), sum(sample_qc$outlier_flag)
  )
), file.path(output_dir, "expression_dimension_audit.tsv"))

grDevices::pdf(file.path(output_dir, "expression_qc.pdf"), 11, 8.5)
print(
  ggplot2::ggplot(pca_dt, ggplot2::aes(PC1, PC2, color = disease)) +
    ggplot2::geom_point(alpha = 0.75) + ggplot2::theme_bw() +
    ggplot2::labs(title = "GSE121239 official RMA matrix PCA")
)
print(
  ggplot2::ggplot(sample_qc, ggplot2::aes(
    expression_median, expression_iqr, color = outlier_flag
  )) +
    ggplot2::geom_point(alpha = 0.75) + ggplot2::theme_bw() +
    ggplot2::labs(title = "Predefined expression QC flags")
)
grDevices::dev.off()

atomic_save_rds(list(
  expression = gene_expression,
  metadata = as.data.frame(meta),
  selected_probes = as.data.frame(chosen),
  sample_qc = as.data.frame(sample_qc),
  normalization = "GEO deposited RMA-normalized series matrix"
), object_file)
write_session_info(file.path(cfg$paths$logs, "01_expression_qc_sessionInfo.txt"))
mark_step_complete(output_dir, signature, outputs)
log_message(
  "Expression preparation complete: ", nrow(gene_expression),
  " genes x ", ncol(gene_expression), " samples."
)

