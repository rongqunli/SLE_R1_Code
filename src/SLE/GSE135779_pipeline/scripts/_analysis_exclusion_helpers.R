if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
mask_excluded_cluster_labels <- function(
    metadata,
    excluded_clusters,
    label_columns,
    exclusion_reason,
    cluster_column = "seurat_clusters") {
  if (!is.data.frame(metadata)) stop("metadata must be a data.frame.")
  required_columns <- unique(c(cluster_column, label_columns))
  missing_columns <- setdiff(required_columns, colnames(metadata))
  if (length(missing_columns) > 0L) {
    stop(
      "Cannot apply cluster exclusions; missing metadata column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }

  excluded_clusters <- unique(as.character(excluded_clusters))
  excluded_clusters <- excluded_clusters[!is.na(excluded_clusters) & nzchar(excluded_clusters)]
  cluster_values <- as.character(metadata[[cluster_column]])
  excluded <- cluster_values %in% excluded_clusters
  excluded[is.na(excluded)] <- FALSE

  original_metadata <- metadata
  metadata$analysis_excluded <- excluded
  metadata$analysis_exclusion_reason <- NA_character_
  metadata$analysis_exclusion_reason[excluded] <- exclusion_reason
  for (column in label_columns) metadata[[column]][excluded] <- NA

  excluded_metadata <- original_metadata[excluded, , drop = FALSE]
  excluded_metadata$analysis_exclusion_reason <- exclusion_reason
  list(
    metadata = metadata,
    excluded_metadata = excluded_metadata,
    excluded = excluded
  )
}
