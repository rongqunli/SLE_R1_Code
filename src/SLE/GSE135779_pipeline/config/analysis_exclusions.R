if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
analysis_exclusions <- list(
  seurat_clusters = c("13"),
  expected_cluster_counts = c(`13` = 1272L),
  reason = "Potential platelet-associated mixed cluster; excluded from primary downstream analyses"
)
