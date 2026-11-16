if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
source(file.path(dirname(normalizePath(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[[1]]), winslash = "/", mustWork = TRUE)), "_common.R"))
annotation_helper_file <- file.path(PIPELINE_ROOT, "scripts", "_annotation_helpers.R")
source(annotation_helper_file)
split_config_file <- file.path(PIPELINE_ROOT, "config", "cluster15_split_config.R")
source(split_config_file)

require_packages(
  c("Seurat", "SeuratObject", "Matrix", "data.table", "ggplot2", "scales"),
  "GSE135779 cell annotation"
)

input_file <- file.path(cfg$paths$processed, "01_qc_clustered_seurat.rds")
output_file <- file.path(cfg$paths$processed, "02_annotated_seurat.rds")
output_dir <- file.path(cfg$paths$results, "04_annotation")
checkpoint_dir <- file.path(cfg$paths$processed, "checkpoints", "04_annotation")
marker_file <- file.path(PIPELINE_ROOT, "config", "pbmc_markers.tsv")
manual_file <- file.path(PIPELINE_ROOT, "config", "manual_cluster_annotations.tsv")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Input object not found: ", input_file, ". Run 03_qc_cluster.R first.")
}
if (!file.exists(marker_file)) {
  stop("Marker configuration not found: ", marker_file)
}

step_signature <- build_step_signature(
  input_files = c(
    input_file, marker_file, manual_file, annotation_helper_file, split_config_file
  ),
  extra = list(
    step = "04_annotation",
    annotation = cfg$annotation,
    cluster15_split = cluster15_split_cfg
  )
)
if (step_is_complete(output_dir, step_signature, output_file)) {
  log_message("Step 04 is already complete and its signature matches. Use --force to rebuild.")
  quit(save = "no", status = 0L)
}
initialize_step(
  step_signature,
  checkpoint_dir,
  owned_paths = c(output_dir, output_file)
)

labels_checkpoint <- file.path(checkpoint_dir, "01_annotation_labels.rds")
annotation_state <- load_checkpoint_rds(labels_checkpoint, step_signature)
sort_cluster_levels <- function(values) {
  levels <- unique(as.character(values))
  numeric_levels <- suppressWarnings(as.numeric(levels))
  if (all(is.finite(numeric_levels))) {
    levels[order(numeric_levels)]
  } else {
    sort(levels)
  }
}
if (is.null(annotation_state)) {
object <- readRDS(input_file)
object$seurat_clusters <- as.character(object$seurat_clusters)
cluster_levels <- sort_cluster_levels(object$seurat_clusters)
marker_sets <- read_marker_table(marker_file)
expression_data <- get_assay_data(object)

marker_sets_matched <- lapply(marker_sets, function(genes) {
  intersect(genes, rownames(expression_data))
})
marker_sets_matched <- marker_sets_matched[lengths(marker_sets_matched) > 0L]

marker_match_summary <- data.frame(
  marker_cell_type = names(marker_sets),
  configured_genes = lengths(marker_sets),
  matched_genes = lengths(marker_sets_matched[names(marker_sets)]),
  stringsAsFactors = FALSE
)
marker_match_summary$matched_genes[is.na(marker_match_summary$matched_genes)] <- 0L

cluster_marker_scores <- list()
for (cluster in cluster_levels) {
  cells <- colnames(object)[object$seurat_clusters == cluster]
  for (marker_name in names(marker_sets_matched)) {
    genes <- marker_sets_matched[[marker_name]]
    values <- expression_data[genes, cells, drop = FALSE]
    cluster_marker_scores[[length(cluster_marker_scores) + 1L]] <- data.frame(
      cluster = cluster,
      marker_cell_type = marker_name,
      mean_expression = mean(Matrix::rowMeans(values)),
      positive_fraction = mean(Matrix::colSums(values > 0) > 0),
      matched_genes = length(genes),
      stringsAsFactors = FALSE
    )
  }
}
cluster_marker_scores <- do.call(rbind, cluster_marker_scores)
cluster_marker_scores$z_score <- ave(
  cluster_marker_scores$mean_expression,
  cluster_marker_scores$marker_cell_type,
  FUN = function(x) {
    if (stats::sd(x) == 0) rep(0, length(x)) else as.numeric(scale(x))
  }
)
marker_suggestion <- do.call(rbind, lapply(split(
  cluster_marker_scores,
  cluster_marker_scores$cluster
), function(x) {
  x <- x[order(-x$z_score, -x$positive_fraction), , drop = FALSE]
  data.frame(
    cluster = x$cluster[[1]],
    marker_suggestion = x$marker_cell_type[[1]],
    marker_z_score = x$z_score[[1]],
    marker_positive_fraction = x$positive_fraction[[1]],
    stringsAsFactors = FALSE
  )
}))

author_labeled <- object@meta.data[
  !is.na(object$author_cell_type_major) & nzchar(object$author_cell_type_major),
  c("seurat_clusters", "author_cell_type_major"),
  drop = FALSE
]

if (nrow(author_labeled) > 0L) {
  author_counts <- as.data.frame(table(
    cluster = author_labeled$seurat_clusters,
    author_cell_type_major = author_labeled$author_cell_type_major
  ), stringsAsFactors = FALSE)
  author_counts <- author_counts[author_counts$Freq > 0, , drop = FALSE]
  author_consensus <- do.call(rbind, lapply(split(author_counts, author_counts$cluster), function(x) {
    x <- x[order(-x$Freq), , drop = FALSE]
    total <- sum(x$Freq)
    data.frame(
      cluster = as.character(x$cluster[[1]]),
      author_consensus = as.character(x$author_cell_type_major[[1]]),
      author_consensus_cells = x$Freq[[1]],
      author_labeled_cells = total,
      author_purity = x$Freq[[1]] / total,
      stringsAsFactors = FALSE
    )
  }))
} else {
  author_consensus <- data.frame(
    cluster = cluster_levels,
    author_consensus = NA_character_,
    author_consensus_cells = 0L,
    author_labeled_cells = 0L,
    author_purity = NA_real_,
    stringsAsFactors = FALSE
  )
}

canonicalize_reference_label <- function(label) {
  label <- as.character(label)
  out <- rep("Unresolved", length(label))
  out[grepl("plasma|plasmablast", label, ignore.case = TRUE)] <- "Plasmablast"
  out[grepl("plasmacytoid", label, ignore.case = TRUE)] <- "pDC"
  out[grepl("dendritic", label, ignore.case = TRUE) & out == "Unresolved"] <- "cDC"
  out[grepl("monocyte|macrophage", label, ignore.case = TRUE)] <- "Monocyte"
  out[grepl("natural killer|^NK", label, ignore.case = TRUE)] <- "NK"
  out[grepl("B cell|B-cell|^B_", label, ignore.case = TRUE)] <- "B"
  out[grepl("T cell|T-cell|^T_", label, ignore.case = TRUE)] <- "T"
  out[grepl("platelet|megakary", label, ignore.case = TRUE)] <- "Platelet"
  out[grepl("eryth", label, ignore.case = TRUE)] <- "Erythroid"
  out
}

marker_to_major <- function(label) {
  label <- as.character(label)
  out <- rep("Unresolved", length(label))
  out[grepl("^CD4_|^Treg$", label)] <- "CD4_T"
  out[grepl("^CD8_", label)] <- "CD8_T"
  out[grepl("^T$", label)] <- "T"
  out[grepl("^NK$", label)] <- "NK"
  out[grepl("^B_", label)] <- "B"
  out[grepl("^Plasmablast$", label)] <- "Plasmablast"
  out[grepl("Monocyte$", label)] <- "Monocyte"
  out[grepl("^cDC$", label)] <- "cDC"
  out[grepl("^pDC$", label)] <- "pDC"
  out[grepl("^Platelet$", label)] <- "Platelet"
  out[grepl("^Erythroid$", label)] <- "Erythroid"
  out[grepl("^Cycling$", label)] <- "Cycling"
  out
}

marker_to_analysis <- function(label) {
  label <- as.character(label)
  out <- rep("Unresolved", length(label))
  detailed_t <- grepl("^CD4_|^CD8_|^Treg$", label)
  out[detailed_t] <- label[detailed_t]
  out[grepl("^NK$", label)] <- "NK"
  out[grepl("^B_", label)] <- label[grepl("^B_", label)]
  out[grepl("^Plasmablast$", label)] <- "Plasmablast"
  out[grepl("^CD14_Monocyte$", label)] <- "CD14_Monocyte"
  out[grepl("^CD16_Monocyte$", label)] <- "CD16_Monocyte"
  out[grepl("^cDC$", label)] <- "cDC"
  out[grepl("^pDC$", label)] <- "pDC"
  out[grepl("^Platelet$", label)] <- "Platelet"
  out[grepl("^Erythroid$", label)] <- "Erythroid"
  out[grepl("^Cycling$", label)] <- "Cycling"
  out
}

singler_table <- data.frame(
  cluster = cluster_levels,
  singler_label = NA_character_,
  singler_major = NA_character_,
  stringsAsFactors = FALSE
)

need_singler <- length(unique(na.omit(object$author_cell_type_major))) == 0L ||
  any(!cluster_levels %in% author_consensus$cluster) ||
  any(
    is.na(author_consensus$author_purity) |
      author_consensus$author_consensus_cells < cfg$annotation$author_cluster_min_cells |
      author_consensus$author_purity < cfg$annotation$author_cluster_min_purity
  )

if (isTRUE(cfg$annotation$run_singler_if_author_labels_missing) && need_singler) {
  if (all(vapply(
    c("SingleR", "celldex", "SingleCellExperiment", "SummarizedExperiment"),
    requireNamespace,
    logical(1),
    quietly = TRUE
  ))) {
    log_message("Running cluster-level SingleR as annotation fallback.")
    cluster_index <- match(object$seurat_clusters, cluster_levels)
    cell_to_cluster <- Matrix::sparseMatrix(
      i = seq_len(ncol(object)),
      j = cluster_index,
      x = 1,
      dims = c(ncol(object), length(cluster_levels))
    )
    cluster_sums <- expression_data %*% cell_to_cluster
    cluster_n <- tabulate(cluster_index, nbins = length(cluster_levels))
    cluster_average <- average_sparse_cluster_columns(cluster_sums, cluster_n)
    cluster_average <- as.matrix(cluster_average)
    colnames(cluster_average) <- cluster_levels

    reference <- celldex::MonacoImmuneData()
    reference_labels <- SummarizedExperiment::colData(reference)$label.main
    prediction <- SingleR::SingleR(
      test = cluster_average,
      ref = reference,
      labels = reference_labels
    )
    singler_table <- data.frame(
      cluster = rownames(prediction),
      singler_label = as.character(prediction$pruned.labels),
      stringsAsFactors = FALSE
    )
    missing_pruned <- is.na(singler_table$singler_label)
    singler_table$singler_label[missing_pruned] <- as.character(prediction$labels[missing_pruned])
    singler_table$singler_major <- canonicalize_reference_label(singler_table$singler_label)
  } else {
    warning("SingleR fallback requested but its packages are unavailable.")
  }
}

review <- data.frame(cluster = cluster_levels, stringsAsFactors = FALSE)
review <- merge(review, author_consensus, by = "cluster", all.x = TRUE)
review <- merge(review, marker_suggestion, by = "cluster", all.x = TRUE)
review <- merge(review, singler_table, by = "cluster", all.x = TRUE)
review$marker_major <- marker_to_major(review$marker_suggestion)

author_high_confidence <- !is.na(review$author_consensus) &
  review$author_consensus_cells >= cfg$annotation$author_cluster_min_cells &
  review$author_purity >= cfg$annotation$author_cluster_min_purity

review$cluster_auto_label <- "Unresolved"
review$cluster_auto_label[author_high_confidence] <- review$author_consensus[author_high_confidence]

use_singler <- review$cluster_auto_label == "Unresolved" &
  !is.na(review$singler_major) & review$singler_major != "Unresolved"
review$cluster_auto_label[use_singler] <- review$singler_major[use_singler]

use_markers <- review$cluster_auto_label == "Unresolved" &
  !is.na(review$marker_major) & review$marker_major != "Unresolved"
review$cluster_auto_label[use_markers] <- review$marker_major[use_markers]

manual <- data.frame(
  cluster = character(0),
  cell_type_major = character(0),
  cell_type_analysis = character(0)
)
if (file.exists(manual_file)) {
  manual <- data.table::fread(
    manual_file,
    sep = "\t",
    header = TRUE,
    data.table = FALSE
  )
  if (!all(c("cluster", "cell_type_major") %in% colnames(manual))) {
    stop("Manual annotation file must contain cluster and cell_type_major columns.")
  }
  if (!"cell_type_analysis" %in% colnames(manual)) {
    manual$cell_type_analysis <- NA_character_
  }
  manual$cluster <- as.character(manual$cluster)
  manual$cell_type_major <- trimws(as.character(manual$cell_type_major))
  manual$cell_type_analysis <- trimws(as.character(manual$cell_type_analysis))
  keep_manual <- (!is.na(manual$cell_type_major) & nzchar(manual$cell_type_major)) |
    (!is.na(manual$cell_type_analysis) & nzchar(manual$cell_type_analysis))
  manual <- manual[keep_manual, , drop = FALSE]
  if (any(!manual$cluster %in% cluster_levels)) {
    stop(
      "Manual annotation contains unknown clusters: ",
      paste(setdiff(manual$cluster, cluster_levels), collapse = ", ")
    )
  }
}

review$manual_label <- NA_character_
if (nrow(manual) > 0L) {
  review$manual_label <- manual$cell_type_major[match(review$cluster, manual$cluster)]
}
review$cluster_final_label <- review$cluster_auto_label
has_manual <- !is.na(review$manual_label) & nzchar(review$manual_label)
review$cluster_final_label[has_manual] <- review$manual_label[has_manual]

analysis_scores <- cluster_marker_scores
analysis_scores$marker_major <- marker_to_major(analysis_scores$marker_cell_type)
analysis_scores$analysis_label <- marker_to_analysis(analysis_scores$marker_cell_type)
analysis_suggestions <- do.call(rbind, lapply(seq_len(nrow(review)), function(i) {
  cluster <- review$cluster[[i]]
  major <- review$cluster_final_label[[i]]
  candidates <- analysis_scores[
    analysis_scores$cluster == cluster &
      analysis_scores$marker_major == major &
      analysis_scores$analysis_label != "Unresolved",
    , drop = FALSE
  ]
  if (nrow(candidates) == 0L) {
    return(data.frame(
      cluster = cluster,
      analysis_auto_label = major,
      analysis_marker_z = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  candidates <- candidates[order(-candidates$z_score, -candidates$positive_fraction), , drop = FALSE]
  data.frame(
    cluster = cluster,
    analysis_auto_label = candidates$analysis_label[[1]],
    analysis_marker_z = candidates$z_score[[1]],
    stringsAsFactors = FALSE
  )
}))
review <- merge(review, analysis_suggestions, by = "cluster", all.x = TRUE, sort = FALSE)
review$manual_analysis_label <- NA_character_
if (nrow(manual) > 0L) {
  review$manual_analysis_label <- manual$cell_type_analysis[
    match(review$cluster, manual$cluster)
  ]
}
review$cluster_analysis_label <- review$analysis_auto_label
has_manual_analysis <- !is.na(review$manual_analysis_label) &
  nzchar(review$manual_analysis_label)
review$cluster_analysis_label[has_manual_analysis] <-
  review$manual_analysis_label[has_manual_analysis]

cluster_final_map <- setNames(review$cluster_final_label, review$cluster)
object$cluster_cell_type_major <- unname(
  cluster_final_map[object$seurat_clusters]
)
object$cell_type_major <- object$cluster_cell_type_major
cluster_analysis_map <- setNames(review$cluster_analysis_label, review$cluster)
object$cluster_cell_type_analysis <- unname(
  cluster_analysis_map[object$seurat_clusters]
)
object$cell_type_analysis <- object$cluster_cell_type_analysis
analysis_major <- marker_to_major(object$cell_type_analysis)
analysis_incompatible <- analysis_major != "Unresolved" &
  object$cell_type_major != "Unresolved" &
  analysis_major != object$cell_type_major
analysis_incompatible[is.na(analysis_incompatible)] <- FALSE
object$cell_type_analysis[analysis_incompatible] <- as.character(
  object$cell_type_major[analysis_incompatible]
)
review$marker_genes_configured <- sum(marker_match_summary$configured_genes)
review$marker_genes_matched <- sum(marker_match_summary$matched_genes)
review$marker_sets_with_matches <- sum(marker_match_summary$matched_genes > 0L)
review$cells_checked_for_analysis_compatibility <- ncol(object)
review$incompatible_analysis_labels_replaced <- sum(analysis_incompatible)
review$incompatible_label_replacement <- "cell_type_major"
object$annotation_source <- rep("cluster_consensus", ncol(object))
if (nrow(manual) > 0L) {
  manual_by_cell <- manual$cell_type_major[match(object$seurat_clusters, manual$cluster)]
  manual_analysis_by_cell <- manual$cell_type_analysis[
    match(object$seurat_clusters, manual$cluster)
  ]
  object$annotation_source[
    (!is.na(manual_by_cell) & nzchar(manual_by_cell)) |
      (!is.na(manual_analysis_by_cell) & nzchar(manual_analysis_by_cell))
  ] <- "manual_cluster_override"
}

split_state <- list(
  assignments = data.frame(),
  removed_cells = data.frame(),
  validation = data.frame()
)
split_cfg <- cluster15_split_cfg
if (isTRUE(split_cfg$enabled)) {
  split_cluster <- as.character(split_cfg$cluster)
  split_cells <- colnames(object)[object$seurat_clusters == split_cluster]
  reference_column <- as.character(split_cfg$reference_column)
  negative_label <- as.character(split_cfg$negative_label)
  positive_label <- as.character(split_cfg$positive_label)
  confidence_threshold <- as.numeric(split_cfg$confidence_threshold)
  if (length(split_cells) == 0L) {
    stop("Configured split cluster is absent: ", split_cluster)
  }
  if (!reference_column %in% colnames(object@meta.data)) {
    stop("Cluster split reference column is absent: ", reference_column)
  }
  if (
    !is.finite(confidence_threshold) ||
      confidence_threshold <= 0.5 || confidence_threshold >= 1
  ) {
    stop("Cluster split confidence_threshold must be between 0.5 and 1.")
  }

  split_scores <- scaled_marker_scores(
    expression_data,
    split_cells,
    list(
      pdc_score = split_cfg$negative_markers,
      plasmablast_score = split_cfg$positive_markers
    )
  )
  reference_labels <- as.character(
    object@meta.data[split_cells, reference_column]
  )
  reference_labels[
    !reference_labels %in% c(negative_label, positive_label)
  ] <- NA_character_
  known_reference <- !is.na(reference_labels)
  reference_counts <- table(factor(
    reference_labels[known_reference],
    levels = c(negative_label, positive_label)
  ))
  if (any(reference_counts < as.integer(split_cfg$k_neighbors) + 1L)) {
    stop(
      "Cluster split reference has too few cells for kNN: ",
      paste(names(reference_counts), reference_counts, sep = "=", collapse = ", ")
    )
  }

  reference_subjects <- as.character(object$subject_id[match(
    split_cells[known_reference], colnames(object)
  )])
  if (any(is.na(reference_subjects) | !nzchar(reference_subjects))) {
    stop("Every cluster split reference cell must have a subject_id.")
  }
  donor_levels <- unique(reference_subjects)
  cv_folds <- min(as.integer(split_cfg$donor_cv_folds), length(donor_levels))
  if (cv_folds < 2L) stop("Cluster split requires at least two donor CV folds.")
  set.seed(cfg$seed + suppressWarnings(as.integer(split_cluster)))
  donor_levels <- sample(donor_levels)
  donor_fold <- setNames(rep(seq_len(cv_folds), length.out = length(donor_levels)), donor_levels)
  reference_fold <- unname(donor_fold[reference_subjects])
  cv_probability <- rep(NA_real_, sum(known_reference))
  reference_scores <- split_scores[known_reference, , drop = FALSE]
  reference_only_labels <- reference_labels[known_reference]
  for (fold in seq_len(cv_folds)) {
    training <- reference_fold != fold
    testing <- reference_fold == fold
    cv_probability[testing] <- marker_knn_probability(
      reference_scores[training, , drop = FALSE],
      reference_only_labels[training],
      reference_scores[testing, , drop = FALSE],
      positive_label = positive_label,
      k = as.integer(split_cfg$k_neighbors),
      exclude_matching_rownames = FALSE
    )
  }
  cv_prediction <- ifelse(cv_probability >= 0.5, positive_label, negative_label)
  cv_confusion <- table(
    reference = factor(reference_only_labels, levels = c(negative_label, positive_label)),
    prediction = factor(cv_prediction, levels = c(negative_label, positive_label))
  )
  cv_accuracy <- mean(cv_prediction == reference_only_labels)
  negative_recall <- cv_confusion[negative_label, negative_label] /
    sum(cv_confusion[negative_label, ])
  positive_recall <- cv_confusion[positive_label, positive_label] /
    sum(cv_confusion[positive_label, ])
  cv_balanced_accuracy <- mean(c(negative_recall, positive_recall))

  split_probability <- marker_knn_probability(
    reference_scores,
    reference_only_labels,
    split_scores,
    positive_label = positive_label,
    k = as.integer(split_cfg$k_neighbors),
    exclude_matching_rownames = TRUE
  )
  lower_threshold <- 1 - confidence_threshold
  threshold_tolerance <- sqrt(.Machine$double.eps)
  split_labels <- ifelse(
    split_probability >= confidence_threshold - threshold_tolerance,
    positive_label,
    ifelse(
      split_probability <= lower_threshold + threshold_tolerance,
      negative_label,
      "Ambiguous"
    )
  )
  split_assignments <- data.frame(
    cell_id = split_cells,
    sample_id = as.character(object$sample_id[match(split_cells, colnames(object))]),
    subject_id = as.character(object$subject_id[match(split_cells, colnames(object))]),
    seurat_cluster = split_cluster,
    pdc_score = split_scores$pdc_score,
    plasmablast_score = split_scores$plasmablast_score,
    probability_plasmablast = split_probability,
    assigned_cell_type = split_labels,
    removed_low_confidence = split_labels == "Ambiguous",
    stringsAsFactors = FALSE
  )
  rownames(split_assignments) <- split_cells

  assigned <- split_labels != "Ambiguous"
  object$cell_type_major[match(split_cells[assigned], colnames(object))] <-
    split_labels[assigned]
  object$cell_type_analysis[match(split_cells[assigned], colnames(object))] <-
    split_labels[assigned]
  object$annotation_source[match(split_cells[assigned], colnames(object))] <-
    "cluster15_marker_knn"
  object$annotation_source[match(split_cells[!assigned], colnames(object))] <-
    "cluster15_low_confidence_removed"

  removed_cells <- split_assignments[split_assignments$removed_low_confidence, , drop = FALSE]
  expected_removed <- as.integer(split_cfg$expected_removed_cells)
  if (nrow(removed_cells) != expected_removed) {
    stop(
      "Cluster split produced ", nrow(removed_cells),
      " low-confidence cells; expected ", expected_removed,
      ". Review markers and thresholds before changing the frozen exclusion set."
    )
  }
  if (isTRUE(split_cfg$remove_ambiguous) && nrow(removed_cells) > 0L) {
    object <- subset(
      object,
      cells = setdiff(colnames(object), removed_cells$cell_id)
    )
  }

  validation <- data.frame(
    metric = c(
      "split_cluster", "cells_before_split", "reference_negative_cells",
      "reference_positive_cells", "donor_cv_folds", "k_neighbors",
      "confidence_threshold", "negative_markers", "positive_markers",
      "expected_removed_cells", "donor_grouped_cv_accuracy",
      "donor_grouped_cv_balanced_accuracy", "negative_label_recall",
      "positive_label_recall", "final_negative_cells", "final_positive_cells",
      "removed_low_confidence_cells"
    ),
    value = as.character(c(
      split_cluster, length(split_cells), reference_counts[[negative_label]],
      reference_counts[[positive_label]], cv_folds, split_cfg$k_neighbors,
      confidence_threshold,
      paste(split_cfg$negative_markers, collapse = ";"),
      paste(split_cfg$positive_markers, collapse = ";"),
      split_cfg$expected_removed_cells, cv_accuracy, cv_balanced_accuracy,
      negative_recall, positive_recall,
      sum(split_labels == negative_label), sum(split_labels == positive_label),
      nrow(removed_cells)
    )),
    stringsAsFactors = FALSE
  )
  split_state <- list(
    assignments = split_assignments,
    removed_cells = removed_cells,
    validation = validation
  )

  review$split_strategy <- NA_character_
  review$split_negative_cells <- NA_integer_
  review$split_positive_cells <- NA_integer_
  review$split_removed_cells <- NA_integer_
  review$split_cv_accuracy <- NA_real_
  review$split_cv_balanced_accuracy <- NA_real_
  split_review_row <- review$cluster == split_cluster
  review$split_strategy[split_review_row] <- "donor_grouped_marker_score_knn"
  review$split_negative_cells[split_review_row] <- sum(split_labels == negative_label)
  review$split_positive_cells[split_review_row] <- sum(split_labels == positive_label)
  review$split_removed_cells[split_review_row] <- nrow(removed_cells)
  review$split_cv_accuracy[split_review_row] <- cv_accuracy
  review$split_cv_balanced_accuracy[split_review_row] <- cv_balanced_accuracy
}

object$cluster_cell_type_major <- object$cell_type_major
object$cluster_cell_type_analysis <- object$cell_type_analysis
object$cluster_annotation_display <- paste0(
  "C", object$seurat_clusters, "\n", object$cell_type_analysis
)

author_annotation_columns <- intersect(
  c("author_subcluster", "author_cell_type_major", "author_label_matched", "cell_type_fine"),
  colnames(object@meta.data)
)
if (length(author_annotation_columns) > 0L) {
  object@meta.data <- object@meta.data[
    , setdiff(colnames(object@meta.data), author_annotation_columns), drop = FALSE
  ]
}

unresolved_fraction <- mean(object$cell_type_major == "Unresolved")
if (unresolved_fraction > cfg$annotation$unresolved_warning_fraction) {
  warning(
    sprintf(
      "%.1f%% of cells remain unresolved. Review cluster_annotation_review.tsv and manual annotations.",
      100 * unresolved_fraction
    )
  )
}

annotation_state <- list(
  object = object,
  review = review,
  marker_sets_matched = marker_sets_matched,
  cluster_marker_scores = cluster_marker_scores,
  split_state = split_state
)
save_checkpoint_rds(
  annotation_state,
  labels_checkpoint,
  step_signature,
  compress = FALSE
)
} else {
  object <- annotation_state$object
  review <- annotation_state$review
  marker_sets_matched <- annotation_state$marker_sets_matched
  cluster_marker_scores <- annotation_state$cluster_marker_scores
  split_state <- annotation_state$split_state
  if (is.null(cluster_marker_scores)) {
    stop(
      "Annotation checkpoint predates cluster marker score persistence. ",
      "Re-run step 04 with --force."
    )
  }
  if (is.null(split_state)) {
    stop("Annotation checkpoint predates cluster split persistence. Re-run step 04 with --force.")
  }
}

safe_write_tsv(review, file.path(output_dir, "cluster_annotation_review.tsv"))
safe_write_tsv(
  split_state$assignments,
  file.path(output_dir, "cluster15_split_cell_assignments.tsv")
)
safe_write_tsv(
  split_state$removed_cells,
  file.path(output_dir, "cluster15_removed_low_confidence_cells.tsv")
)
safe_write_tsv(
  split_state$validation,
  file.path(output_dir, "cluster15_split_validation.tsv")
)

cluster_levels <- sort_cluster_levels(object$seurat_clusters)
cluster_marker_scores$cluster <- factor(
  as.character(cluster_marker_scores$cluster),
  levels = cluster_levels
)
cluster_marker_scores$marker_cell_type <- factor(
  as.character(cluster_marker_scores$marker_cell_type),
  levels = names(marker_sets_matched)
)
cluster_marker_scores <- cluster_marker_scores[
  order(cluster_marker_scores$cluster, cluster_marker_scores$marker_cell_type),
  , drop = FALSE
]
cluster_marker_scores$cluster <- as.character(cluster_marker_scores$cluster)
cluster_marker_scores$marker_cell_type <- as.character(
  cluster_marker_scores$marker_cell_type
)
safe_write_tsv(
  cluster_marker_scores,
  file.path(output_dir, "cluster_marker_scores.tsv")
)

set.seed(cfg$seed)
cells_by_cluster <- split(colnames(object), object$seurat_clusters)
marker_cells <- unlist(lapply(cells_by_cluster, function(cells) {
  if (length(cells) <= cfg$annotation$marker_cells_per_cluster) cells else
    sample(cells, cfg$annotation$marker_cells_per_cluster)
}), use.names = FALSE)
marker_object <- object[, marker_cells, drop = FALSE]
Seurat::Idents(marker_object) <- marker_object$seurat_clusters
cluster_markers <- Seurat::FindAllMarkers(
  marker_object,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = cfg$annotation$marker_min_pct,
  logfc.threshold = cfg$annotation$marker_logfc,
  verbose = FALSE
)
safe_write_tsv(cluster_markers, file.path(output_dir, "cluster_markers_for_annotation.tsv"))

marker_panel_map <- c(
  CD4_T_Naive = "lymphoid",
  CD4_T_Memory = "lymphoid",
  Treg = "lymphoid",
  CD8_T_Naive = "lymphoid",
  CD8_T_GZMK = "lymphoid",
  CD8_T_Cytotoxic = "lymphoid",
  NK = "lymphoid",
  B_Naive = "lymphoid",
  B_Memory = "lymphoid",
  B_ABC = "lymphoid",
  Plasmablast = "lymphoid",
  CD14_Monocyte = "myeloid_dc",
  CD16_Monocyte = "myeloid_dc",
  cDC = "myeloid_dc",
  pDC = "myeloid_dc",
  Platelet = "other",
  Erythroid = "other",
  Cycling = "other"
)
unassigned_marker_types <- setdiff(names(marker_sets_matched), names(marker_panel_map))
if (length(unassigned_marker_types) > 0L) {
  stop(
    "Marker cell types are missing a DotPlot panel assignment: ",
    paste(unassigned_marker_types, collapse = ", ")
  )
}

marker_entry_table <- data.table::rbindlist(lapply(
  seq_along(marker_sets_matched),
  function(marker_type_index) {
    marker_type <- names(marker_sets_matched)[[marker_type_index]]
    genes <- marker_sets_matched[[marker_type]]
    data.frame(
      marker_cell_type = marker_type,
      marker_cell_type_order = marker_type_index,
      gene = genes,
      gene_order_within_type = seq_along(genes),
      plot_panel = unname(marker_panel_map[[marker_type]]),
      stringsAsFactors = FALSE
    )
  }
))
marker_entry_table <- as.data.frame(marker_entry_table)
marker_entry_table$marker_entry_order <- seq_len(nrow(marker_entry_table))
marker_entry_table$marker_key <- paste0(
  sprintf("%03d", marker_entry_table$marker_entry_order),
  "__",
  marker_entry_table$gene
)
marker_entry_table$marker_cell_type_label <- gsub(
  "_",
  "\n",
  marker_entry_table$marker_cell_type,
  fixed = TRUE
)

unique_marker_genes <- unique(marker_entry_table$gene)
cluster_dotplot_data <- Seurat::DotPlot(
  object,
  features = unique_marker_genes,
  group.by = "seurat_clusters",
  assay = "RNA",
  scale = TRUE
)$data
cluster_dotplot_data <- as.data.frame(cluster_dotplot_data)
colnames(cluster_dotplot_data)[
  colnames(cluster_dotplot_data) == "features.plot"
] <- "gene"
cluster_dotplot_data$id <- as.character(cluster_dotplot_data$id)
cluster_dotplot_data <- merge(
  marker_entry_table,
  cluster_dotplot_data,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)
cluster_dotplot_data <- cluster_dotplot_data[
  order(
    cluster_dotplot_data$marker_entry_order,
    match(cluster_dotplot_data$id, cluster_levels)
  ),
  , drop = FALSE
]
cluster_dotplot_data$cluster_label <- factor(
  paste0("Cluster ", cluster_dotplot_data$id),
  levels = rev(paste0("Cluster ", cluster_levels))
)
cluster_dotplot_data$marker_key <- factor(
  cluster_dotplot_data$marker_key,
  levels = marker_entry_table$marker_key
)
marker_type_labels <- unique(marker_entry_table$marker_cell_type_label)
cluster_dotplot_data$marker_cell_type_label <- factor(
  cluster_dotplot_data$marker_cell_type_label,
  levels = marker_type_labels
)
cluster_dotplot_data$avg.exp.scaled[!is.finite(cluster_dotplot_data$avg.exp.scaled)] <- 0
cluster_dotplot_data$pct.exp[!is.finite(cluster_dotplot_data$pct.exp)] <- 0
gene_axis_labels <- stats::setNames(
  marker_entry_table$gene,
  marker_entry_table$marker_key
)

build_cluster_marker_dotplot <- function(plot_data, title) {
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = marker_key,
      y = cluster_label,
      size = pct.exp,
      color = avg.exp.scaled
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::facet_grid(
      cols = ggplot2::vars(marker_cell_type_label),
      scales = "free_x",
      space = "free_x"
    ) +
    ggplot2::scale_x_discrete(labels = gene_axis_labels, drop = TRUE) +
    ggplot2::scale_size_continuous(
      name = "Expressing cells (%)",
      range = c(0.2, 5.2),
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100)
    ) +
    ggplot2::scale_color_gradient2(
      name = "Scaled average\nexpression",
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = scales::squish
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing.x = grid::unit(0.08, "lines"),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 6
      ),
      axis.text.y = ggplot2::element_text(size = 7),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey60"),
      strip.text.x = ggplot2::element_text(size = 7, face = "bold"),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    ) +
    ggplot2::labs(title = title, x = "Marker gene", y = NULL)
}

save_marker_dotplot <- function(plot_data, filename, title) {
  plot <- build_cluster_marker_dotplot(plot_data, title)
  marker_entries <- length(unique(as.character(plot_data$marker_key)))
  ggplot2::ggsave(
    file.path(output_dir, filename),
    plot,
    width = max(12, min(30, marker_entries * 0.19 + 4)),
    height = max(8, length(cluster_levels) * 0.38 + 3),
    units = "in",
    limitsize = FALSE
  )
}

save_marker_dotplot(
  cluster_dotplot_data,
  "cluster_marker_dotplot.pdf",
  "Cluster marker overview"
)
panel_specifications <- list(
  lymphoid = list(
    filename = "cluster_marker_dotplot_lymphoid.pdf",
    title = "Cluster markers - lymphoid"
  ),
  myeloid_dc = list(
    filename = "cluster_marker_dotplot_myeloid_dc.pdf",
    title = "Cluster markers - myeloid and DC"
  ),
  other = list(
    filename = "cluster_marker_dotplot_other.pdf",
    title = "Cluster markers - other populations"
  )
)
for (panel_name in names(panel_specifications)) {
  specification <- panel_specifications[[panel_name]]
  save_marker_dotplot(
    cluster_dotplot_data[cluster_dotplot_data$plot_panel == panel_name, , drop = FALSE],
    specification$filename,
    specification$title
  )
}

score_plot_data <- cluster_marker_scores
score_plot_data$cluster_label <- factor(
  paste0("Cluster ", score_plot_data$cluster),
  levels = rev(paste0("Cluster ", cluster_levels))
)
score_plot_data$marker_cell_type <- factor(
  score_plot_data$marker_cell_type,
  levels = names(marker_sets_matched)
)
score_plot_data$is_cluster_maximum <- ave(
  score_plot_data$z_score,
  score_plot_data$cluster,
  FUN = function(values) as.integer(values == max(values, na.rm = TRUE))
) > 0
score_heatmap <- ggplot2::ggplot(
  score_plot_data,
  ggplot2::aes(x = marker_cell_type, y = cluster_label, fill = z_score)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.25) +
  ggplot2::geom_point(
    data = score_plot_data[which(score_plot_data$is_cluster_maximum), , drop = FALSE],
    shape = 8,
    size = 2,
    color = "black"
  ) +
  ggplot2::scale_fill_gradient2(
    name = "Marker Z-score",
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = ggplot2::element_text(size = 8),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5)
  ) +
  ggplot2::labs(
    title = "Cluster marker-set score heatmap",
    subtitle = "Asterisk marks the highest marker Z-score within each cluster",
    x = "Candidate cell type",
    y = NULL
  )
ggplot2::ggsave(
  file.path(output_dir, "cluster_marker_score_heatmap.pdf"),
  score_heatmap,
  width = 13,
  height = max(8, length(cluster_levels) * 0.4 + 2.5),
  units = "in"
)

umap_celltype <- Seurat::DimPlot(
  object,
  reduction = "umap",
  group.by = "cell_type_major",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) + ggplot2::ggtitle("GSE135779 final major cell types")

umap_analysis <- Seurat::DimPlot(
  object,
  reduction = "umap",
  group.by = "cell_type_analysis",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggplot2::ggtitle("GSE135779 analysis-level cell types") +
  ggplot2::guides(
    colour = ggplot2::guide_legend(
      ncol = 2,
      override.aes = list(size = 3)
    )
  ) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 7))

umap_cluster_analysis <- Seurat::DimPlot(
  object,
  reduction = "umap",
  group.by = "cluster_annotation_display",
  label = TRUE,
  label.size = 3,
  repel = TRUE,
  raster = TRUE
) +
  ggplot2::ggtitle("Cluster-level candidate annotations") +
  Seurat::NoLegend()

cluster15_cells_final <- colnames(object)[object$seurat_clusters == as.character(
  cluster15_split_cfg$cluster
)]
cluster15_umap_coordinates <- SeuratObject::Embeddings(
  object,
  reduction = "umap"
)[cluster15_cells_final, , drop = FALSE]
cluster15_zoom_x <- as.numeric(stats::quantile(
  cluster15_umap_coordinates[, 1], c(0.01, 0.99)
))
cluster15_zoom_y <- as.numeric(stats::quantile(
  cluster15_umap_coordinates[, 2], c(0.01, 0.99)
))
umap_cluster15_split <- Seurat::DimPlot(
  subset(object, cells = cluster15_cells_final),
  reduction = "umap",
  group.by = "cell_type_analysis",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggplot2::coord_fixed(xlim = cluster15_zoom_x, ylim = cluster15_zoom_y) +
  ggplot2::labs(
    title = "Cluster 15: pDC and Plasmablast split",
    subtitle = "Central 98% view; original UMAP coordinates retained"
  )

ggplot2::ggsave(
  file.path(output_dir, "umap_cell_type_major.pdf"),
  umap_celltype,
  width = 10,
  height = 8
)
ggplot2::ggsave(
  file.path(output_dir, "umap_cell_type_analysis.pdf"),
  umap_analysis,
  width = 12,
  height = 9
)
ggplot2::ggsave(
  file.path(output_dir, "umap_cluster_analysis_annotation.pdf"),
  umap_cluster_analysis,
  width = 12,
  height = 9
)
ggplot2::ggsave(
  file.path(output_dir, "umap_cluster15_pdc_plasmablast_split.pdf"),
  umap_cluster15_split,
  width = 8,
  height = 7
)

annotation_counts <- as.data.frame(table(
  sample_id = object$sample_id,
  cohort = object$cohort,
  condition = object$condition,
  cell_type_major = object$cell_type_major,
  cell_type_analysis = object$cell_type_analysis
), stringsAsFactors = FALSE)
safe_write_tsv(annotation_counts, file.path(output_dir, "cell_type_counts_by_sample.tsv"))

object@misc$annotation_review <- review
atomic_save_rds(object, output_file, compress = FALSE)
legacy_dotplot <- file.path(output_dir, "pbmc_marker_dotplot.pdf")
if (file.exists(legacy_dotplot)) unlink(legacy_dotplot, force = TRUE)
write_session_info(file.path(cfg$paths$logs, "04_annotation_sessionInfo.txt"))
critical_outputs <- c(
  output_file,
  file.path(output_dir, "cluster_annotation_review.tsv"),
  file.path(output_dir, "cluster_markers_for_annotation.tsv"),
  file.path(output_dir, "cluster_marker_scores.tsv"),
  file.path(output_dir, "cell_type_counts_by_sample.tsv"),
  file.path(output_dir, "cluster15_split_cell_assignments.tsv"),
  file.path(output_dir, "cluster15_removed_low_confidence_cells.tsv"),
  file.path(output_dir, "cluster15_split_validation.tsv"),
  file.path(output_dir, "cluster_marker_dotplot.pdf"),
  file.path(output_dir, "cluster_marker_dotplot_lymphoid.pdf"),
  file.path(output_dir, "cluster_marker_dotplot_myeloid_dc.pdf"),
  file.path(output_dir, "cluster_marker_dotplot_other.pdf"),
  file.path(output_dir, "cluster_marker_score_heatmap.pdf"),
  file.path(output_dir, "umap_cell_type_major.pdf"),
  file.path(output_dir, "umap_cell_type_analysis.pdf"),
  file.path(output_dir, "umap_cluster_analysis_annotation.pdf"),
  file.path(output_dir, "umap_cluster15_pdc_plasmablast_split.pdf")
)
complete_step_and_cleanup(
  output_dir,
  step_signature,
  critical_outputs,
  cleanup_paths = c(
    checkpoint_dir,
    file.path(output_dir, "umap_cell_type_fine.pdf"),
    file.path(output_dir, "umap_author_subclusters.pdf")
  )
)

log_message("Annotated object saved: ", output_file)
