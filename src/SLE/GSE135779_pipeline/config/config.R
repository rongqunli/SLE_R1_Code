if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
if (!exists("PIPELINE_ROOT", inherits = FALSE)) {
  stop("PIPELINE_ROOT is not defined. Run scripts with Rscript from this pipeline.")
}

cfg <- list(
  accession = "GSE135779",
  seed = 20260725L,
  meta = list(
    version = "2026-07-25_v1.0",
    frozen_master = normalizePath(
      file.path(
        PIPELINE_ROOT,
        "..",
        "metadata",
        "frozen",
        "2026-07-25_v1.0",
        "by_dataset",
        "GSE135779_meta.tsv"
      ),
      winslash = "/",
      mustWork = FALSE
    )
  ),
  paths = list(
    raw = file.path(PIPELINE_ROOT, "data", "raw"),
    metadata = file.path(PIPELINE_ROOT, "data", "metadata"),
    processed = file.path(PIPELINE_ROOT, "data", "processed"),
    results = file.path(PIPELINE_ROOT, "results"),
    logs = file.path(PIPELINE_ROOT, "results", "logs")
  ),
  urls = list(
    raw_tar = paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135779/suppl/",
      "GSE135779_RAW.tar"
    ),
    genes = paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135779/suppl/",
      "GSE135779_genes.tsv.gz"
    ),
    author_metadata = paste0(
      "https://raw.githubusercontent.com/dnehar/SingleCells_SLE_paper/master/",
      "Meta_caSLE_processed_08092021_small.csv"
    )
  ),
  files = list(
    raw_tar = "GSE135779_RAW.tar",
    genes = "GSE135779_genes.tsv.gz",
    author_metadata = "Meta_caSLE_processed_08092021_small.csv",
    extracted_dir = "GSE135779_RAW"
  ),
  data = list(
    expected_geo_matrices = 56L
  ),
  clinical = list(
    activity_low_max = 3,
    activity_high_min_exclusive = 6
  ),
  inference = list(
    confidence_level = 0.95
  ),
  qc = list(
    min_features = 200L,
    max_features_hard = 7500L,
    max_percent_mt = 15,
    mad_multiplier = 4,
    mad_constant = 1,
    iqr_normalizer = 1.349,
    fallback_upper_quantile = 0.99,
    min_cells_per_sample_after_qc = 100L,
    drop_samples_below_min_cells = TRUE,
    run_scdblfinder = TRUE,
    keep_doublet_class = "singlet"
  ),
  dimension_reduction = list(
    normalization_scale_factor = 10000,
    n_variable_features = 3000L,
    n_pcs = 50L,
    
    
    dims_use = 1:15,
    run_harmony = TRUE,
    harmony_group = "batch",
    
    
    resolution = 0.6,
    resolution_grid = seq(0.2, 1.0, by = 0.2),
    umap_neighbors = 30L,
    umap_min_dist = 0.3
  ),
  parameter_evaluation = list(
    cumulative_variance_threshold = 0.90,
    pc_auto_min = 5L,
    pc_elbow_buffer = 5L,
    silhouette_max_cells = 3000L,
    umap_comparison_max_cells = 20000L,
    small_cluster_cell_warning = 100L
  ),
  annotation = list(
    author_cluster_min_cells = 50L,
    author_cluster_min_purity = 0.60,
    marker_cells_per_cluster = 2000L,
    unresolved_warning_fraction = 0.05,
    marker_min_pct = 0.10,
    marker_logfc = 0.25,
    run_singler_if_author_labels_missing = TRUE
  ),
  pseudobulk = list(
    min_cells_major = 30L,
    min_cells_analysis = 30L,
    min_cells_subcluster = 20L,
    min_samples_per_group = 5L,
    min_samples_continuous = 12L,
    min_expressed_genes = 10L,
    fdr = 0.05,
    logfc = 0.5,
    adult_support_fdr = 0.20,
    abundance_logit_pseudocount = 0.5,
    covariates_child = c("batch", "age", "sex"),
    covariates_adult = c("batch", "age"),
    covariates_activity = c("batch", "age", "sex")
  ),
  programs = list(
    min_matched_genes = 3L,
    min_cells_per_donor = 20L,
    min_samples_per_group = 5L
  )
)
