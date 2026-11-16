if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
average_sparse_cluster_columns <- function(cluster_sums, cluster_n) {
  if (is.null(dim(cluster_sums)) || length(dim(cluster_sums)) != 2L) {
    stop("cluster_sums must be a two-dimensional matrix-like object.")
  }
  cluster_n <- as.numeric(cluster_n)
  if (length(cluster_n) != ncol(cluster_sums)) {
    stop("cluster_n length must equal the number of cluster_sums columns.")
  }
  if (any(!is.finite(cluster_n)) || any(cluster_n <= 0)) {
    stop("Every cluster must contain at least one cell.")
  }
  cluster_sums %*% Matrix::Diagonal(x = 1 / cluster_n)
}

scaled_marker_scores <- function(expression_data, cells, marker_sets) {
  cells <- as.character(cells)
  if (length(cells) == 0L || any(!cells %in% colnames(expression_data))) {
    stop("Every requested cell must be present in expression_data.")
  }
  if (!is.list(marker_sets) || length(marker_sets) < 2L || is.null(names(marker_sets))) {
    stop("marker_sets must be a named list containing at least two marker sets.")
  }
  marker_sets <- lapply(marker_sets, function(genes) {
    unique(intersect(as.character(genes), rownames(expression_data)))
  })
  if (any(lengths(marker_sets) == 0L)) {
    stop(
      "No expression genes matched marker set(s): ",
      paste(names(marker_sets)[lengths(marker_sets) == 0L], collapse = ", ")
    )
  }

  genes <- unique(unlist(marker_sets, use.names = FALSE))
  values <- as.matrix(expression_data[genes, cells, drop = FALSE])
  scaled <- t(scale(t(values)))
  scaled[!is.finite(scaled)] <- 0
  scores <- vapply(marker_sets, function(set_genes) {
    colMeans(scaled[set_genes, , drop = FALSE])
  }, numeric(length(cells)))
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  rownames(scores) <- cells
  scores
}

marker_knn_probability <- function(
    training_scores,
    training_labels,
    query_scores,
    positive_label,
    k = 15L,
    exclude_matching_rownames = TRUE) {
  training_scores <- as.matrix(training_scores)
  query_scores <- as.matrix(query_scores)
  training_labels <- as.character(training_labels)
  if (nrow(training_scores) != length(training_labels)) {
    stop("training_labels length must match training_scores rows.")
  }
  if (ncol(training_scores) != ncol(query_scores)) {
    stop("Training and query score matrices must contain the same features.")
  }
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 1L) {
    stop("k must be one positive integer.")
  }
  k <- as.integer(k)
  if (nrow(training_scores) < k) {
    stop("Training data contain fewer rows than k.")
  }

  center <- colMeans(training_scores)
  spread <- apply(training_scores, 2, stats::sd)
  spread[!is.finite(spread) | spread == 0] <- 1
  training_scaled <- sweep(training_scores, 2, center, "-")
  training_scaled <- sweep(training_scaled, 2, spread, "/")
  query_scaled <- sweep(query_scores, 2, center, "-")
  query_scaled <- sweep(query_scaled, 2, spread, "/")

  vapply(seq_len(nrow(query_scaled)), function(i) {
    distances <- rowSums((training_scaled - matrix(
      query_scaled[i, ],
      nrow = nrow(training_scaled),
      ncol = ncol(training_scaled),
      byrow = TRUE
    ))^2)
    if (
      isTRUE(exclude_matching_rownames) &&
        !is.null(rownames(training_scaled)) &&
        !is.null(rownames(query_scaled))
    ) {
      distances[rownames(training_scaled) == rownames(query_scaled)[[i]]] <- Inf
    }
    available <- sum(is.finite(distances))
    if (available < k) stop("Fewer than k eligible neighbors remain after exclusion.")
    neighbors <- order(distances)[seq_len(k)]
    mean(training_labels[neighbors] == positive_label)
  }, numeric(1))
}
