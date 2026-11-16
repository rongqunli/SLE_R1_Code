if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
source(file.path(dirname(normalizePath(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[[1]]), winslash = "/", mustWork = TRUE)), "_common.R"))

require_packages(
  c(
    "Seurat", "SeuratObject", "data.table", "ggplot2", "patchwork",
    "clustree", "ggraph", "cluster"
  ),
  "PCA and clustering parameter evaluation"
)

input_file <- file.path(cfg$paths$processed, "01_qc_clustered_seurat.rds")
output_dir <- file.path(cfg$paths$results, "03_parameter_evaluation")
checkpoint_dir <- file.path(
  cfg$paths$processed,
  "checkpoints",
  "03_parameter_evaluation"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("QC/clustering object not found: ", input_file, ". Run step 03_qc first.")
}

evaluation_cfg <- cfg$parameter_evaluation
step_signature <- build_step_signature(
  input_files = input_file,
  extra = list(
    step = "03_parameter_evaluation",
    parameter_evaluation = evaluation_cfg
  )
)

critical_outputs <- c(
  file.path(output_dir, "pca_variance.tsv"),
  file.path(output_dir, "pc_selection_summary.tsv"),
  file.path(output_dir, "resolution_metrics.tsv"),
  file.path(output_dir, "resolution_transition_stability.tsv"),
  file.path(output_dir, "pca_elbow_and_variance.pdf"),
  file.path(output_dir, "resolution_clustree.pdf"),
  file.path(output_dir, "resolution_metric_comparison.pdf"),
  file.path(output_dir, "resolution_umap_comparison.pdf")
)

if (step_is_complete(output_dir, step_signature, critical_outputs)) {
  log_message("Step 03_params is already complete and its signature matches. Use --force to rebuild.")
  quit(save = "no", status = 0L)
}

initialize_step(
  step_signature,
  checkpoint_dir,
  owned_paths = output_dir
)

detect_elbow_pc <- function(stdev_values) {
  values <- as.numeric(stdev_values)
  values <- values[is.finite(values)]
  n_dimensions <- length(values)
  if (n_dimensions <= 2L) return(max(1L, n_dimensions))

  x_scaled <- (seq_len(n_dimensions) - 1) / (n_dimensions - 1)
  value_range <- range(values)
  if (!all(is.finite(value_range)) || diff(value_range) <= 0) {
    return(max(2L, min(n_dimensions, 5L)))
  }
  y_scaled <- (values - value_range[[1]]) / diff(value_range)
  chord <- y_scaled[[1]] + (y_scaled[[n_dimensions]] - y_scaled[[1]]) * x_scaled
  as.integer(which.max(chord - y_scaled))
}

adjusted_rand_index <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  x <- as.character(x[keep])
  y <- as.character(y[keep])
  n <- length(x)
  if (n < 2L) return(NA_real_)
  contingency <- table(x, y)
  choose_two <- function(values) sum(values * (values - 1) / 2)
  index <- choose_two(contingency)
  row_pairs <- choose_two(rowSums(contingency))
  column_pairs <- choose_two(colSums(contingency))
  total_pairs <- n * (n - 1) / 2
  expected <- row_pairs * column_pairs / total_pairs
  maximum <- (row_pairs + column_pairs) / 2
  denominator <- maximum - expected
  if (!is.finite(denominator) || abs(denominator) < .Machine$double.eps) {
    return(1)
  }
  (index - expected) / denominator
}

dominant_fraction <- function(parent, child, condition_on = c("child", "parent")) {
  condition_on <- match.arg(condition_on)
  contingency <- table(parent = parent, child = child)
  if (condition_on == "child") {
    sum(apply(contingency, 2, max)) / sum(contingency)
  } else {
    sum(apply(contingency, 1, max)) / sum(contingency)
  }
}

stratified_sample_cells <- function(cells, groups, maximum_cells, seed) {
  cells <- as.character(cells)
  groups <- as.character(groups)
  if (length(cells) <= maximum_cells) return(cells)
  set.seed(seed)
  grouped <- split(cells, groups)
  quota <- max(1L, floor(maximum_cells / length(grouped)))
  selected <- unlist(lapply(grouped, function(group_cells) {
    sample(group_cells, min(length(group_cells), quota))
  }), use.names = FALSE)
  remaining_slots <- maximum_cells - length(selected)
  if (remaining_slots > 0L) {
    remaining <- setdiff(cells, selected)
    selected <- c(selected, sample(remaining, min(remaining_slots, length(remaining))))
  }
  selected
}

object <- readRDS(input_file)
if (!"pca" %in% names(object@reductions)) {
  stop("The step 03 object does not contain a PCA reduction.")
}

pca_stdev <- as.numeric(object[["pca"]]@stdev)
pca_stdev <- pca_stdev[is.finite(pca_stdev)]
if (length(pca_stdev) < 2L) stop("Too few valid PCA dimensions for evaluation.")

pca_variance <- pca_stdev^2
pca_variance_ratio <- pca_variance / sum(pca_variance)
pca_cumulative <- cumsum(pca_variance_ratio)
elbow_pc <- detect_elbow_pc(pca_stdev)
cumulative_pc <- which(
  pca_cumulative >= evaluation_cfg$cumulative_variance_threshold
)[[1]]
recommended_pc <- max(
  elbow_pc,
  min(cumulative_pc, elbow_pc + evaluation_cfg$pc_elbow_buffer)
)
recommended_pc <- min(
  length(pca_stdev),
  max(evaluation_cfg$pc_auto_min, recommended_pc)
)

object_cfg <- object@misc$analysis_config %||% list()
object_dimensions <- object_cfg$dimension_reduction$dims_use %||%
  cfg$dimension_reduction$dims_use
object_pc <- max(as.integer(object_dimensions))
configured_pc <- max(as.integer(cfg$dimension_reduction$dims_use))
object_resolution <- object_cfg$dimension_reduction$resolution %||%
  cfg$dimension_reduction$resolution

pca_table <- data.frame(
  pc = seq_along(pca_stdev),
  standard_deviation = pca_stdev,
  variance = pca_variance,
  variance_fraction = pca_variance_ratio,
  cumulative_variance = pca_cumulative,
  elbow_pc = seq_along(pca_stdev) == elbow_pc,
  recommended_pc = seq_along(pca_stdev) == recommended_pc,
  configured_pc = seq_along(pca_stdev) == configured_pc,
  stringsAsFactors = FALSE
)
safe_write_tsv(pca_table, file.path(output_dir, "pca_variance.tsv"))

selection_summary <- data.frame(
  available_pcs = length(pca_stdev),
  elbow_pc = elbow_pc,
  cumulative_variance_threshold = evaluation_cfg$cumulative_variance_threshold,
  cumulative_variance_pc = cumulative_pc,
  elbow_buffer = evaluation_cfg$pc_elbow_buffer,
  recommended_pc_for_review = recommended_pc,
  pc_used_in_object = object_pc,
  pc_currently_configured = configured_pc,
  resolution_used_in_object = as.numeric(object_resolution),
  resolution_currently_configured = cfg$dimension_reduction$resolution,
  object_and_config_pc_match = identical(object_pc, configured_pc),
  note = paste0(
    "Review the elbow plot before editing config.R. The recommended PC is the elbow, ",
    "extended by at most the configured buffer without exceeding the cumulative-variance PC."
  ),
  stringsAsFactors = FALSE
)
safe_write_tsv(
  selection_summary,
  file.path(output_dir, "pc_selection_summary.tsv")
)

elbow_plot <- ggplot2::ggplot(
  pca_table,
  ggplot2::aes(x = pc, y = standard_deviation)
) +
  ggplot2::geom_line(linewidth = 0.6, color = "#2C3E50") +
  ggplot2::geom_point(size = 1.2, color = "#2C3E50") +
  ggplot2::geom_vline(
    xintercept = elbow_pc,
    color = "#D55E00",
    linewidth = 0.7,
    linetype = 2
  ) +
  ggplot2::geom_vline(
    xintercept = recommended_pc,
    color = "#0072B2",
    linewidth = 0.7,
    linetype = 3
  ) +
  ggplot2::geom_vline(
    xintercept = configured_pc,
    color = "#009E73",
    linewidth = 0.7,
    linetype = 4
  ) +
  ggplot2::annotate(
    "label",
    x = elbow_pc,
    y = max(pca_stdev),
    label = paste0("Elbow = ", elbow_pc),
    hjust = -0.05,
    vjust = 1,
    size = 3.2,
    color = "#D55E00"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::labs(
    title = "PCA elbow",
    subtitle = paste0(
      "Blue: review recommendation = ", recommended_pc,
      "; green: configured = ", configured_pc
    ),
    x = "Principal component",
    y = "Standard deviation"
  )

cumulative_plot <- ggplot2::ggplot(
  pca_table,
  ggplot2::aes(x = pc, y = cumulative_variance)
) +
  ggplot2::geom_line(linewidth = 0.6, color = "#2C3E50") +
  ggplot2::geom_point(size = 1.1, color = "#2C3E50") +
  ggplot2::geom_hline(
    yintercept = evaluation_cfg$cumulative_variance_threshold,
    color = "#CC79A7",
    linewidth = 0.7,
    linetype = 2
  ) +
  ggplot2::geom_vline(
    xintercept = cumulative_pc,
    color = "#CC79A7",
    linewidth = 0.7,
    linetype = 2
  ) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::labs(
    title = "Cumulative variance among computed PCs",
    subtitle = paste0(
      evaluation_cfg$cumulative_variance_threshold * 100,
      "% is reached at PC ", cumulative_pc
    ),
    x = "Principal component",
    y = "Cumulative variance"
  )

ggplot2::ggsave(
  file.path(output_dir, "pca_elbow_and_variance.pdf"),
  elbow_plot / cumulative_plot,
  width = 9,
  height = 9,
  units = "in"
)

metadata <- object@meta.data
resolution_columns_all <- grep("_snn_res\\.", colnames(metadata), value = TRUE)
if (length(resolution_columns_all) < 2L) {
  stop("At least two *_snn_res.* metadata columns are required. Re-run step 03_qc.")
}
resolution_prefixes <- sub(
  "[0-9]+(?:\\.[0-9]+)?$",
  "",
  resolution_columns_all,
  perl = TRUE
)
resolution_prefix <- names(sort(table(resolution_prefixes), decreasing = TRUE))[[1]]
resolution_columns <- resolution_columns_all[
  startsWith(resolution_columns_all, resolution_prefix)
]
resolution_values <- suppressWarnings(as.numeric(sub(
  paste0("^", gsub("([.])", "\\\\\\1", resolution_prefix)),
  "",
  resolution_columns
)))
valid_resolution <- is.finite(resolution_values)
resolution_columns <- resolution_columns[valid_resolution]
resolution_values <- resolution_values[valid_resolution]
resolution_order <- order(resolution_values)
resolution_columns <- resolution_columns[resolution_order]
resolution_values <- resolution_values[resolution_order]
if (length(resolution_columns) < 2L) {
  stop("Could not parse at least two clustering resolutions from metadata.")
}

downstream_reduction <- object@misc$downstream_reduction %||%
  if ("harmony" %in% names(object@reductions)) "harmony" else "pca"
if (!downstream_reduction %in% names(object@reductions)) {
  stop("Downstream reduction is absent from the object: ", downstream_reduction)
}
embedding <- SeuratObject::Embeddings(object, reduction = downstream_reduction)
evaluation_dimensions <- seq_len(min(object_pc, ncol(embedding)))
embedding <- embedding[, evaluation_dimensions, drop = FALSE]

highest_clusters <- metadata[[resolution_columns[[length(resolution_columns)]]]]
names(highest_clusters) <- rownames(metadata)
silhouette_cells <- stratified_sample_cells(
  rownames(metadata),
  highest_clusters,
  evaluation_cfg$silhouette_max_cells,
  cfg$seed
)
embedding_distance <- stats::dist(embedding[silhouette_cells, , drop = FALSE])

donor_column <- c("donor_id", "patient_id", "sample_id")
donor_column <- donor_column[donor_column %in% colnames(metadata)][[1]]

resolution_metrics <- lapply(seq_along(resolution_columns), function(index) {
  column <- resolution_columns[[index]]
  resolution <- resolution_values[[index]]
  clusters <- as.character(metadata[[column]])
  cluster_counts <- table(clusters)
  donor_counts <- tapply(
    as.character(metadata[[donor_column]]),
    clusters,
    function(values) length(unique(values[!is.na(values) & nzchar(values)]))
  )
  sampled_clusters <- as.integer(factor(clusters[match(silhouette_cells, rownames(metadata))]))
  silhouette_values <- cluster::silhouette(sampled_clusters, embedding_distance)[, "sil_width"]
  data.frame(
    resolution = resolution,
    metadata_column = column,
    selected_in_object = isTRUE(all.equal(resolution, as.numeric(object_resolution))),
    n_clusters = length(cluster_counts),
    min_cluster_cells = min(cluster_counts),
    median_cluster_cells = stats::median(cluster_counts),
    max_cluster_cells = max(cluster_counts),
    clusters_below_warning_size = sum(
      cluster_counts < evaluation_cfg$small_cluster_cell_warning
    ),
    min_donors_per_cluster = min(donor_counts),
    median_donors_per_cluster = stats::median(donor_counts),
    silhouette_cells = length(silhouette_values),
    mean_silhouette = mean(silhouette_values),
    median_silhouette = stats::median(silhouette_values),
    negative_silhouette_fraction = mean(silhouette_values < 0),
    stringsAsFactors = FALSE
  )
})
resolution_metrics <- do.call(rbind, resolution_metrics)
safe_write_tsv(
  resolution_metrics,
  file.path(output_dir, "resolution_metrics.tsv")
)

transition_rows <- lapply(seq_len(length(resolution_columns) - 1L), function(index) {
  parent <- metadata[[resolution_columns[[index]]]]
  child <- metadata[[resolution_columns[[index + 1L]]]]
  data.frame(
    resolution_from = resolution_values[[index]],
    resolution_to = resolution_values[[index + 1L]],
    adjusted_rand_index = adjusted_rand_index(parent, child),
    child_from_single_parent_fraction = dominant_fraction(parent, child, "child"),
    parent_retained_in_single_child_fraction = dominant_fraction(parent, child, "parent"),
    stringsAsFactors = FALSE
  )
})
transition_table <- do.call(rbind, transition_rows)
safe_write_tsv(
  transition_table,
  file.path(output_dir, "resolution_transition_stability.tsv")
)

clustree_metadata <- metadata[, resolution_columns, drop = FALSE]
clustree_plot <- clustree::clustree(
  clustree_metadata,
  prefix = resolution_prefix
) +
  ggraph::scale_edge_colour_continuous(guide = "none") +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::labs(
    title = "Cluster transitions across resolutions",
    subtitle = paste0(
      "Current object: PC1-", object_pc,
      ", selected resolution ", object_resolution
    )
  )
ggplot2::ggsave(
  file.path(output_dir, "resolution_clustree.pdf"),
  clustree_plot,
  width = 15,
  height = 10,
  units = "in",
  limitsize = FALSE
)

cluster_number_plot <- ggplot2::ggplot(
  resolution_metrics,
  ggplot2::aes(x = resolution, y = n_clusters)
) +
  ggplot2::geom_line(color = "#0072B2", linewidth = 0.7) +
  ggplot2::geom_point(color = "#0072B2", size = 2) +
  ggplot2::geom_vline(
    xintercept = as.numeric(object_resolution),
    linetype = 2,
    color = "#D55E00"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::labs(title = "Cluster count", x = "Resolution", y = "Clusters")

silhouette_plot <- ggplot2::ggplot(
  resolution_metrics,
  ggplot2::aes(x = resolution, y = mean_silhouette)
) +
  ggplot2::geom_line(color = "#009E73", linewidth = 0.7) +
  ggplot2::geom_point(color = "#009E73", size = 2) +
  ggplot2::geom_vline(
    xintercept = as.numeric(object_resolution),
    linetype = 2,
    color = "#D55E00"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::labs(
    title = "Mean silhouette (stratified sample)",
    x = "Resolution",
    y = "Mean silhouette"
  )

small_cluster_plot <- ggplot2::ggplot(
  resolution_metrics,
  ggplot2::aes(x = resolution, y = min_cluster_cells)
) +
  ggplot2::geom_hline(
    yintercept = evaluation_cfg$small_cluster_cell_warning,
    linetype = 3,
    color = "#CC79A7"
  ) +
  ggplot2::geom_line(color = "#E69F00", linewidth = 0.7) +
  ggplot2::geom_point(color = "#E69F00", size = 2) +
  ggplot2::geom_vline(
    xintercept = as.numeric(object_resolution),
    linetype = 2,
    color = "#D55E00"
  ) +
  ggplot2::scale_y_log10() +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::labs(
    title = "Smallest cluster",
    subtitle = paste0(
      "Dotted horizontal line: warning threshold = ",
      evaluation_cfg$small_cluster_cell_warning
    ),
    x = "Resolution",
    y = "Cells (log10 scale)"
  )

ggplot2::ggsave(
  file.path(output_dir, "resolution_metric_comparison.pdf"),
  cluster_number_plot / silhouette_plot / small_cluster_plot,
  width = 8,
  height = 11,
  units = "in"
)

if (!"umap" %in% names(object@reductions)) {
  stop("The step 03 object does not contain a UMAP reduction.")
}
umap <- as.data.frame(SeuratObject::Embeddings(object, reduction = "umap"))
colnames(umap)[1:2] <- c("UMAP_1", "UMAP_2")
umap$cell_id <- rownames(umap)
umap_cells <- stratified_sample_cells(
  rownames(metadata),
  highest_clusters,
  evaluation_cfg$umap_comparison_max_cells,
  cfg$seed
)
umap <- umap[umap_cells, , drop = FALSE]
resolution_labels <- paste0(
  "resolution ",
  formatC(resolution_values, format = "f", digits = 1)
)
umap_long <- data.table::rbindlist(lapply(seq_along(resolution_columns), function(index) {
  data.frame(
    UMAP_1 = umap$UMAP_1,
    UMAP_2 = umap$UMAP_2,
    resolution = factor(
      resolution_labels[[index]],
      levels = resolution_labels
    ),
    cluster = factor(metadata[umap$cell_id, resolution_columns[[index]]]),
    stringsAsFactors = FALSE
  )
}))
umap_comparison_plot <- ggplot2::ggplot(
  umap_long,
  ggplot2::aes(x = UMAP_1, y = UMAP_2, color = cluster)
) +
  ggplot2::geom_point(size = 0.08, alpha = 0.75, stroke = 0) +
  ggplot2::facet_wrap(~resolution, ncol = 2) +
  ggplot2::coord_equal() +
  ggplot2::theme_void(base_size = 10) +
  ggplot2::theme(
    legend.position = "none",
    strip.text = ggplot2::element_text(face = "bold", size = 10),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
  ) +
  ggplot2::labs(
    title = paste0(
      "Resolution comparison on the existing UMAP (PC1-", object_pc, ")"
    )
  )
ggplot2::ggsave(
  file.path(output_dir, "resolution_umap_comparison.pdf"),
  umap_comparison_plot,
  width = 10,
  height = 12,
  units = "in"
)

write_session_info(
  file.path(cfg$paths$logs, "03_parameter_evaluation_sessionInfo.txt")
)
complete_step_and_cleanup(
  output_dir,
  step_signature,
  critical_outputs,
  cleanup_paths = checkpoint_dir
)

log_message("Parameter evaluation results saved: ", output_dir)
log_message(
  "PCA elbow = ", elbow_pc,
  "; review recommendation = ", recommended_pc,
  "; configured/object PC = ", configured_pc, "/", object_pc, "."
)
