if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript generate_single_cell_qc_summary.R <pipeline_results_root> <audit_output_root>")
}

results_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
audit_root <- normalizePath(args[[2]], winslash = "/", mustWork = FALSE)
figure_dir <- file.path(audit_root, "figures")
table_dir <- file.path(audit_root, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

read_tsv <- function(path) {
  utils::read.delim(path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

qc_dir <- file.path(results_root, "03_qc_cluster")
parameter_dir <- file.path(results_root, "03_parameter_evaluation")
annotation_dir <- file.path(results_root, "04_annotation")

thresholds <- read_tsv(file.path(qc_dir, "sample_specific_qc_thresholds.tsv"))
doublets <- read_tsv(file.path(qc_dir, "scDblFinder_by_sample.tsv"))
qc_summary <- read_tsv(file.path(qc_dir, "qc_summary_by_sample.tsv"))
pc_summary <- read_tsv(file.path(parameter_dir, "pc_selection_summary.tsv"))
resolution_metrics <- read_tsv(file.path(parameter_dir, "resolution_metrics.tsv"))
resolution_stability <- read_tsv(file.path(parameter_dir, "resolution_transition_stability.tsv"))
cluster15_validation <- read_tsv(file.path(annotation_dir, "cluster15_split_validation.tsv"))

before_metadata <- qc_summary[qc_summary$stage == "before_qc", c("sample_id", "cohort", "condition")]
sample_audit <- merge(thresholds, doublets, by = "sample_id", all.x = TRUE, sort = FALSE)
sample_audit <- merge(sample_audit, before_metadata, by = "sample_id", all.x = TRUE, sort = FALSE)
sample_audit$cohort_condition <- paste(sample_audit$cohort, sample_audit$condition, sep = "_")
sample_audit$basic_retention <- sample_audit$cells_after_basic_qc / sample_audit$cells_before
sample_audit$singlet_retention <- sample_audit$singlets / sample_audit$cells_before

raw_cells <- sum(sample_audit$cells_before)
basic_cells <- sum(sample_audit$cells_after_basic_qc)
singlet_cells <- sum(sample_audit$singlets)
removed_low_confidence <- as.numeric(cluster15_validation$value[
  cluster15_validation$metric == "removed_low_confidence_cells"
])
if (length(removed_low_confidence) != 1L || !is.finite(removed_low_confidence)) {
  stop("Could not obtain the cluster-15 low-confidence removal count.")
}
annotated_cells <- singlet_cells - removed_low_confidence

cell_flow <- data.frame(
  stage = factor(
    c("Raw matrix", "After basic QC", "After scDblFinder", "Final annotated atlas"),
    levels = c("Raw matrix", "After basic QC", "After scDblFinder", "Final annotated atlas")
  ),
  cells = c(raw_cells, basic_cells, singlet_cells, annotated_cells),
  removed_from_previous = c(NA, raw_cells - basic_cells, basic_cells - singlet_cells, singlet_cells - annotated_cells),
  stringsAsFactors = FALSE
)

retention_long <- rbind(
  data.frame(
    sample_id = sample_audit$sample_id,
    cohort_condition = sample_audit$cohort_condition,
    stage = "After basic QC",
    retention = sample_audit$basic_retention
  ),
  data.frame(
    sample_id = sample_audit$sample_id,
    cohort_condition = sample_audit$cohort_condition,
    stage = "After scDblFinder",
    retention = sample_audit$singlet_retention
  )
)
sample_order <- sample_audit$sample_id[order(sample_audit$singlet_retention)]
retention_long$sample_id <- factor(retention_long$sample_id, levels = sample_order)
retention_long$stage <- factor(retention_long$stage, levels = c("After basic QC", "After scDblFinder"))

threshold_long <- rbind(
  data.frame(sample_id = sample_audit$sample_id, metric = "Maximum detected genes", value = sample_audit$nFeature_max),
  data.frame(sample_id = sample_audit$sample_id, metric = "Maximum UMI count", value = sample_audit$nCount_max)
)
threshold_long$sample_id <- factor(threshold_long$sample_id, levels = sample_order)

palette <- c(
  Adult_HD = "#5B8FF9",
  Adult_SLE = "#5AD8A6",
  Child_HD = "#F6BD16",
  Child_SLE = "#E8684A"
)

p_flow <- ggplot2::ggplot(cell_flow, ggplot2::aes(stage, cells, fill = stage)) +
  ggplot2::geom_col(width = 0.72, show.legend = FALSE) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::comma(cells)),
    vjust = -0.35,
    size = 3.5
  ) +
  ggplot2::scale_fill_manual(values = c("#6C757D", "#4C78A8", "#59A14F", "#B279A2")) +
  ggplot2::scale_y_continuous(labels = scales::label_comma(), expand = ggplot2::expansion(mult = c(0, 0.12))) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 22, hjust = 1)) +
  ggplot2::labs(title = "A  Cell retention across preprocessing", x = NULL, y = "Cells")

p_retention <- ggplot2::ggplot(
  retention_long,
  ggplot2::aes(sample_id, retention, color = cohort_condition, shape = stage, group = stage)
) +
  ggplot2::geom_point(size = 1.7, alpha = 0.9) +
  ggplot2::scale_color_manual(values = palette, name = NULL) +
  ggplot2::scale_shape_manual(values = c(16, 1), name = NULL) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0.65, 1.01)) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) +
  ggplot2::labs(
    title = "B  Per-sample retention",
    subtitle = "Samples ordered by final singlet retention",
    x = "56 samples",
    y = "Retained fraction"
  )

p_doublet <- ggplot2::ggplot(
  sample_audit,
  ggplot2::aes(cohort_condition, doublet_fraction, fill = cohort_condition)
) +
  ggplot2::geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.78, show.legend = FALSE) +
  ggplot2::geom_jitter(width = 0.12, size = 1.3, alpha = 0.75, show.legend = FALSE) +
  ggplot2::scale_fill_manual(values = palette) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 22, hjust = 1)) +
  ggplot2::labs(
    title = "C  Sample-level doublet calls",
    subtitle = paste0(scales::comma(sum(sample_audit$doublets)), " cells removed"),
    x = NULL,
    y = "Doublet fraction"
  )

p_threshold <- ggplot2::ggplot(
  threshold_long,
  ggplot2::aes(sample_id, value, group = 1)
) +
  ggplot2::geom_line(color = "#4C78A8", linewidth = 0.45) +
  ggplot2::geom_point(color = "#4C78A8", size = 0.9) +
  ggplot2::facet_wrap(~metric, ncol = 1, scales = "free_y") +
  ggplot2::scale_y_continuous(labels = scales::label_comma()) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey94", color = "grey70")
  ) +
  ggplot2::labs(
    title = "D  Sample-adaptive upper thresholds",
    subtitle = "Median + 4 MAD; detected genes capped at 7,500; mitochondrial RNA capped at 15%",
    x = "Samples in the same order as panel B",
    y = NULL
  )

qc_plot <- (p_flow | p_retention) / (p_doublet | p_threshold) +
  patchwork::plot_annotation(
    title = "GSE135779 single-cell quality-control audit",
    subtitle = paste0(
      scales::comma(raw_cells), " raw cells -> ", scales::comma(annotated_cells),
      " final cells across 56 donors"
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(size = 10)
    )
  )

pdf_path <- file.path(figure_dir, "01_QC_cell_retention_summary.pdf")
png_path <- file.path(figure_dir, "01_QC_cell_retention_summary.png")
ggplot2::ggsave(pdf_path, qc_plot, width = 14, height = 9.5, units = "in")
ggplot2::ggsave(png_path, qc_plot, width = 14, height = 9.5, units = "in", dpi = 320)

parameter_summary <- data.frame(
  parameter = c(
    "raw_cells", "after_basic_qc", "after_scdblfinder", "final_annotated_cells",
    "minimum_features", "maximum_percent_mitochondrial", "upper_threshold_rule",
    "variable_features", "computed_pcs", "elbow_pc", "review_recommended_pc",
    "selected_pcs", "cumulative_90_percent_pc", "harmony_group",
    "resolution_grid", "selected_resolution", "selected_resolution_clusters",
    "selected_resolution_smallest_cluster", "selected_resolution_mean_silhouette"
  ),
  value = c(
    raw_cells, basic_cells, singlet_cells, annotated_cells,
    200, 15, "per sample: median + 4 MAD; nFeature_RNA capped at 7500",
    3000, pc_summary$available_pcs, pc_summary$elbow_pc, pc_summary$recommended_pc_for_review,
    paste0("PC1-", pc_summary$pc_used_in_object), pc_summary$cumulative_variance_pc, "batch",
    "0.2;0.4;0.6;0.8;1.0", pc_summary$resolution_used_in_object,
    resolution_metrics$n_clusters[resolution_metrics$selected_in_object],
    resolution_metrics$min_cluster_cells[resolution_metrics$selected_in_object],
    signif(resolution_metrics$mean_silhouette[resolution_metrics$selected_in_object], 4)
  ),
  stringsAsFactors = FALSE
)

write_tsv(cell_flow, file.path(table_dir, "qc_cell_flow.tsv"))
write_tsv(sample_audit, file.path(table_dir, "qc_sample_audit.tsv"))
write_tsv(parameter_summary, file.path(table_dir, "single_cell_preprocessing_parameter_summary.tsv"))
write_tsv(resolution_stability, file.path(table_dir, "resolution_transition_stability_copy.tsv"))

cat("Created:\n", pdf_path, "\n", png_path, "\n", sep = "")
