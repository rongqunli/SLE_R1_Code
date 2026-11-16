if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
source(file.path(dirname(normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]]
), winslash = "/", mustWork = TRUE)), "_common.R"))

require_packages(c("data.table", "matrixStats", "ggplot2"), "expression preparation and QC")
source(file.path(PIPELINE_ROOT, "scripts", "_expression_helpers.R"), local = TRUE)

script_path <- script_file()
expression_file <- file.path(cfg$paths$raw, cfg$files$expression)
platform_file <- file.path(cfg$paths$raw, cfg$files$platform)
metadata_file <- cfg$paths$frozen_meta
output_dir <- file.path(cfg$paths$results, "01_expression_qc")
object_file <- file.path(cfg$paths$processed, "01_gene_expression.rds")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

outputs <- c(
  object_file,
  file.path(output_dir, c(
    "probe_gene_mapping.tsv", "probe_selection_audit.tsv", "sample_qc_metrics.tsv",
    "pca_coordinates.tsv", "expression_qc.pdf"
  ))
)
signature <- build_step_signature(
  input_files = c(expression_file, platform_file, metadata_file),
  code_files = c(script_path, file.path(PIPELINE_ROOT, "scripts", "_common.R"),
                 file.path(PIPELINE_ROOT, "scripts", "_expression_helpers.R"),
                 file.path(PIPELINE_ROOT, "config", "config.R"))
)
if (!has_flag("--force") && step_is_complete(output_dir, signature, outputs)) {
  log_message("Step 01 expression preparation is already complete.")
  quit(save = "no", status = 0L)
}

log_message("Reading frozen metadata and 972 included expression columns.")
metadata <- data.table::fread(metadata_file, sep = "\t", na.strings = NULL)
metadata <- metadata[analysis_include == "yes"]
if (nrow(metadata) != cfg$expected$included_arrays) stop("Unexpected included metadata rows.")

selected_columns <- c("ID_REF", metadata$sample_id)
expression_dt <- data.table::fread(
  expression_file,
  skip = "ID_REF",
  select = selected_columns,
  quote = "\"",
  showProgress = TRUE
)
data.table::setnames(expression_dt, names(expression_dt)[[1]], "probe_id")
expression_dt <- expression_dt[!grepl("^!", probe_id)]
if (anyDuplicated(expression_dt$probe_id)) stop("Duplicate probe IDs in series matrix.")
if (!all(metadata$sample_id %in% names(expression_dt))) stop("Missing included expression columns.")

probe_ids <- expression_dt$probe_id
expression <- select_expression_matrix(expression_dt, metadata$sample_id)
rownames(expression) <- probe_ids
if (any(!is.finite(expression))) stop("Non-finite values in normalized expression matrix.")
rm(expression_dt); gc(FALSE)

log_message("Reading GPL10558 annotation and resolving unique gene symbols.")
mapping <- read_platform_gene_symbols(platform_file)
mapping[, mapping_status := data.table::fcase(
  !nzchar(gene_symbol_raw) | gene_symbol_raw == "---", "unmapped",
  grepl(" /// ", gene_symbol_raw, fixed = TRUE), "multi_symbol",
  default = "unique_symbol"
)]
mapping[, gene_symbol := ifelse(mapping_status == "unique_symbol", toupper(gene_symbol_raw), NA_character_)]
mapping <- unique(mapping, by = "probe_id")
safe_write_tsv(mapping, file.path(output_dir, "probe_gene_mapping.tsv"))

mapped <- merge(
  data.table::data.table(probe_id = rownames(expression), matrix_row = seq_len(nrow(expression))),
  mapping[mapping_status == "unique_symbol", .(probe_id, gene_symbol)],
  by = "probe_id"
)
if (!nrow(mapped)) stop("No uniquely mapped probes overlap the expression matrix.")

probe_variance <- matrixStats::rowVars(expression[mapped$matrix_row, , drop = FALSE])
selection <- mapped[, .(probe_id, gene_symbol, matrix_row)]
selection[, probe_variance := probe_variance]
selection[, selected := FALSE]
data.table::setorder(selection, gene_symbol, -probe_variance, probe_id)
selection[, selected := seq_len(.N) == 1L, by = gene_symbol]
selection[!is.finite(probe_variance) | probe_variance < cfg$expression$min_gene_variance, selected := FALSE]

chosen <- selection[selected == TRUE]
data.table::setorder(chosen, gene_symbol)
gene_expression <- expression[chosen$matrix_row, , drop = FALSE]
rownames(gene_expression) <- chosen$gene_symbol
colnames(gene_expression) <- metadata$sample_id
rm(expression); gc(FALSE)

safe_write_tsv(selection[, .(probe_id, gene_symbol, probe_variance, selected)],
               file.path(output_dir, "probe_selection_audit.tsv"))

log_message("Calculating sample QC metrics and PCA.")
sample_median <- matrixStats::colMedians(gene_expression)
sample_iqr <- matrixStats::colIQRs(gene_expression)
sample_sd <- matrixStats::colSds(gene_expression)
gene_variance <- matrixStats::rowVars(gene_expression)
top_n <- min(2000L, nrow(gene_expression))
top_rows <- order(gene_variance, decreasing = TRUE)[seq_len(top_n)]
qc_expression <- gene_expression[top_rows, , drop = FALSE]
correlation <- stats::cor(qc_expression, method = "pearson")
diag(correlation) <- NA_real_
median_correlation <- matrixStats::colMedians(correlation, na.rm = TRUE)

pca_fit <- stats::prcomp(t(qc_expression), center = TRUE, scale. = FALSE, rank. = cfg$expression$pca_components)
pca_coordinates <- data.table::as.data.table(pca_fit$x, keep.rownames = "sample_id")
pca_coordinates <- merge(pca_coordinates, metadata, by = "sample_id", sort = FALSE)
pca_coordinates <- pca_coordinates[match(metadata$sample_id, sample_id)]

pc_count <- min(5L, ncol(pca_fit$x))
robust_pc_distance <- rowSums(scale(pca_fit$x[, seq_len(pc_count), drop = FALSE])^2)
distance_cutoff <- stats::median(robust_pc_distance) +
  cfg$expression$outlier_mad_multiplier * stats::mad(robust_pc_distance)
correlation_cutoff <- stats::median(median_correlation) -
  cfg$expression$outlier_mad_multiplier * stats::mad(median_correlation)

sample_qc <- data.table::data.table(
  sample_id = metadata$sample_id,
  expression_median = sample_median,
  expression_iqr = sample_iqr,
  expression_sd = sample_sd,
  missing_fraction = 0,
  median_correlation = median_correlation,
  robust_pc_distance = robust_pc_distance,
  low_correlation_flag = median_correlation < correlation_cutoff,
  pca_outlier_flag = robust_pc_distance > distance_cutoff
)
sample_qc[, outlier_flag := low_correlation_flag | pca_outlier_flag]
sample_qc <- merge(sample_qc, metadata, by = "sample_id", sort = FALSE)
sample_qc <- sample_qc[match(metadata$sample_id, sample_id)]

safe_write_tsv(sample_qc, file.path(output_dir, "sample_qc_metrics.tsv"))
safe_write_tsv(pca_coordinates, file.path(output_dir, "pca_coordinates.tsv"))

sample_qc_plot <- prepare_plot_metadata(sample_qc)
pca_plot <- prepare_plot_metadata(pca_coordinates)
pdf(file.path(output_dir, "expression_qc.pdf"), width = 11, height = 8.5, onefile = TRUE)
print(ggplot2::ggplot(sample_qc_plot, ggplot2::aes(expression_median, expression_iqr, color = disease_state)) +
  ggplot2::geom_point(alpha = 0.7, size = 1.7) +
  ggplot2::facet_wrap(~ batch) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "GSE65391 normalized-expression distributions", color = "State"))
print(ggplot2::ggplot(pca_plot, ggplot2::aes(PC1, PC2, color = disease_state, shape = batch)) +
  ggplot2::geom_point(alpha = 0.75, size = 2) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "GSE65391 PCA", color = "State", shape = "Batch"))
print(ggplot2::ggplot(sample_qc_plot, ggplot2::aes(median_correlation, robust_pc_distance, color = outlier_flag)) +
  ggplot2::geom_point(alpha = 0.75, size = 1.8) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "Predefined expression-QC outlier flags", color = "Flag"))
grDevices::dev.off()

object <- list(
  expression = gene_expression,
  metadata = as.data.frame(metadata),
  selected_probes = as.data.frame(chosen[, .(probe_id, gene_symbol, probe_variance)]),
  sample_qc = as.data.frame(sample_qc),
  pca_variance = pca_fit$sdev^2 / sum(pca_fit$sdev^2)
)
atomic_save_rds(object, object_file, compress = TRUE)
write_session_info(file.path(cfg$paths$logs, "01_expression_qc_sessionInfo.txt"))
mark_step_complete(output_dir, signature, outputs)
log_message("Expression QC complete: ", nrow(gene_expression), " genes x ",
            ncol(gene_expression), " included arrays; ", sum(sample_qc$outlier_flag), " flagged samples.")
