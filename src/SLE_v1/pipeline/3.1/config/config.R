if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
if (!exists("PIPELINE_ROOT", inherits = FALSE)) {
  stop("PIPELINE_ROOT is not defined. Run scripts through run_pipeline.R or Rscript.")
}

PROJECT_ROOT <- normalizePath(file.path(PIPELINE_ROOT, "..", ".."), winslash = "/", mustWork = TRUE)

cfg <- list(
  accession = "GSE65391",
  seed = 20260726L,
  meta_version = "2026-07-25_v1.0",
  paths = list(
    raw = file.path(PIPELINE_ROOT, "data", "raw"),
    metadata = file.path(PIPELINE_ROOT, "data", "metadata"),
    processed = file.path(PIPELINE_ROOT, "data", "processed"),
    results = file.path(PROJECT_ROOT, "results", "3.1"),
    logs = file.path(PROJECT_ROOT, "results", "3.1", "logs"),
    frozen_meta = normalizePath(
      file.path(PIPELINE_ROOT, "data", "frozen", "2026-07-25_v1.0",
                "by_dataset", "GSE65391_meta.tsv"),
      winslash = "/", mustWork = FALSE
    )
  ),
  files = list(
    expression = "GSE65391_series_matrix.txt.gz",
    platform = "GPL10558.annot.gz",
    programs = file.path(PIPELINE_ROOT, "config", "programs.tsv")
  ),
  expected = list(
    total_arrays = 996L,
    included_arrays = 972L,
    technical_replicates = 24L,
    sle_visits = 924L,
    sle_patients = 158L,
    healthy_arrays = 48L,
    malformed_visit_ids = 1L
  ),
  expression = list(
    min_gene_variance = 1e-8,
    min_program_genes = 3L,
    outlier_mad_multiplier = 6,
    pca_components = 10L
  ),
  longitudinal = list(
    primary_program = "CYTOTOXICITY_SHARED",
    fdr = 0.05,
    primary_formula = paste(
      "score ~ SLEDAI_within + SLEDAI_between + timepoint_scaled +",
      "(1 | patient_id)"
    ),
    composition_formula = paste(
      "score ~ SLEDAI_within + SLEDAI_between + timepoint_scaled +",
      "NK + T_CELL + NEUTROPHIL + ERYTHROID + PLATELET +",
      "(1 | patient_id)"
    )
  )
)
