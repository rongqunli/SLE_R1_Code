if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
if (!exists("PIPELINE_ROOT", inherits = FALSE)) {
  stop("PIPELINE_ROOT is not defined. Run scripts through run_pipeline.R or Rscript.")
}

cfg <- list(
  accession = "GSE121239",
  seed = 20260729L,
  meta_version = "2026-07-25_v1.0",
  paths = list(
    raw = file.path(PIPELINE_ROOT, "data", "raw"),
    downloads = file.path(PIPELINE_ROOT, "data", "raw", "downloads"),
    metadata = file.path(PIPELINE_ROOT, "data", "metadata"),
    processed = file.path(PIPELINE_ROOT, "data", "processed"),
    results = file.path(PIPELINE_ROOT, "results"),
    logs = file.path(PIPELINE_ROOT, "results", "logs"),
    frozen_meta = normalizePath(file.path(
      PIPELINE_ROOT, "..", "metadata", "frozen", "2026-07-25_v1.0",
      "by_dataset", "GSE121239_meta.tsv"
    ), winslash = "/", mustWork = FALSE),
    raw_geo_meta = normalizePath(file.path(
      PIPELINE_ROOT, "..", "metadata", "interim", "geo_sample_metadata",
      "GSE121239_series_matrix.txt.sample_metadata.tsv"
    ), winslash = "/", mustWork = FALSE),
    primary_programs = normalizePath(file.path(
      PIPELINE_ROOT, "..", "GSE200306_pipeline", "results",
      "07_cross_dataset_bridge", "expanded_program_membership.tsv"
    ), winslash = "/", mustWork = FALSE),
    secondary_programs = normalizePath(file.path(
      PIPELINE_ROOT, "..", "GSE65391_pipeline", "config", "programs.tsv"
    ), winslash = "/", mustWork = FALSE),
    frozen_shortlist = normalizePath(file.path(
      PIPELINE_ROOT, "..", "GSE65391_pipeline", "results",
      "08_key_gene_integration", "key_gene_shortlist.tsv"
    ), winslash = "/", mustWork = FALSE),
    gse65391_module_models = normalizePath(file.path(
      PIPELINE_ROOT, "..", "GSE65391_pipeline", "results",
      "06_wgcna", "module_trait_associations.tsv"
    ), winslash = "/", mustWork = FALSE),
    design_contract = file.path(PIPELINE_ROOT, "config", "analysis_contract.tsv")
  ),
  files = list(
    raw_archive = "GSE121239_RAW.tar",
    series_matrix = "GSE121239_series_matrix.txt.gz",
    platform = "GPL13158.annot.gz",
    filelist = "filelist.txt"
  ),
  downloads = data.frame(
    file = c(
      "GSE121239_series_matrix.txt.gz", "GPL13158.annot.gz",
      "GSE121239_RAW.tar", "filelist.txt"
    ),
    url = c(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121239/matrix/GSE121239_series_matrix.txt.gz",
      "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13158/annot/GPL13158.annot.gz",
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121239/suppl/GSE121239_RAW.tar",
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121239/suppl/filelist.txt"
    ),
    expected_bytes = c(92993584, 8313277, 668651520, 24652),
    expected_sha256 = c(
      "37BE0D11D78139158E4A67624C7E45D209D96E12956E0E27E09CECD1057D1F67",
      "D9EDDD653C200C9AC1B282B441079A19E890914929B35792A4E99CE4EBA3B0D0",
      "99D3F5FB3B315846AB0ABAA869BBB37BD69B3693D4B36FE4B5F41E0029B16CEC",
      "BE54D297EF083CEFA155E811072D724ACAFF73F6A66BDEE177767791A2594D60"
    ),
    stringsAsFactors = FALSE
  ),
  expected = list(
    samples = 312L,
    sle_visits = 292L,
    healthy_samples = 20L,
    sle_patients = 65L,
    healthy_donors = 20L,
    cel_files = 312L,
    probes = 54715L,
    neutrophil_visits = 241L,
    patients_with_sledai_variation = 65L,
    top_genes = 100L
  ),
  analysis = list(
    primary_programs = c(
      "CYAN_NEUTROPHIL_IMMATURE",
      "GREENYELLOW_B_CELL_COMPLEMENT",
      "MAGENTA_IFN_ANTIVIRAL"
    ),
    secondary_programs = c("CYTOTOXICITY_SHARED", "CD8_GZMK_STATE"),
    min_program_genes = 5L,
    min_patient_visits = 2L,
    outlier_mad_multiplier = 6,
    fdr = 0.05,
    supportive_fdr = 0.10,
    loo_direction_fraction = 0.90
  )
)


cfg$paths <- cfg$paths[c("raw","downloads","metadata","processed","results","logs","frozen_meta","raw_geo_meta")]
