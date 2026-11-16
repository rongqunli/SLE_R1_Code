if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
cluster15_split_cfg <- list(
  enabled = TRUE,
  cluster = "15",
  reference_column = "author_cell_type_major",
  negative_label = "pDC",
  positive_label = "Plasmablast",
  negative_markers = c("GZMB", "CLEC4C", "LILRA4", "TCF4", "IL3RA", "SCT"),
  positive_markers = c(
    "MZB1", "XBP1", "PRDM1", "SDC1", "TNFRSF17", "DERL3", "JCHAIN", "IGJ"
  ),
  k_neighbors = 15L,
  confidence_threshold = 0.80,
  donor_cv_folds = 5L,
  remove_ambiguous = TRUE,
  expected_removed_cells = 8L
)
