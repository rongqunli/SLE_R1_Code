if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
study_root <- Sys.getenv("SLE_STUDY_ROOT")
if (!nzchar(study_root)) stop("SLE_STUDY_ROOT is not set. Run through run_pipeline.R.")

for (package in c("SeuratObject", "Matrix", "data.table", "UCell", "ggplot2")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Required package is missing: ", package)
}

source(file.path(study_root, "scripts", "_helpers.R"), local = FALSE)
contract <- read_contract(file.path(study_root, "config", "analysis_contract.tsv"))
input_file <- resolve_contract_path(study_root, contract$input_object)
membership_file <- resolve_contract_path(study_root, contract$module_membership)
pseudobulk_file <- resolve_contract_path(study_root, contract$pseudobulk_results)
results_dir <- file.path(study_root, "results", "01_projection")
figures_dir <- file.path(study_root, "figures")
score_cache_file <- file.path(results_dir, "frozen_program_score_matrices.rds")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

from_aggregates <- identical(Sys.getenv("SLE_FROM_AGGREGATES"), "true")
scores_file <- file.path(results_dir, "donor_celltype_program_scores.tsv")
proportions_file <- file.path(results_dir, "donor_celltype_proportions.tsv")
audit_file <- file.path(results_dir, "frozen_program_gene_matching.tsv")

if (from_aggregates) {
  required <- c(scores_file, proportions_file, audit_file)
  if (any(!file.exists(required))) {
    stop("--from-aggregates requires existing donor score, proportion, and gene audit tables.")
  }
  donor_scores <- data.table::fread(scores_file)
  donor_proportions <- data.table::fread(proportions_file)
  audit <- data.table::fread(audit_file)
  program_state <- list(
    audit = audit,
    programs = stats::setNames(
      lapply(audit$matched_gene_list[audit$passes_scoring_gate], split_values),
      audit$program[audit$passes_scoring_gate]
    )
  )
  n_cells_input <- NA_integer_
  n_cells_primary <- sum(unique(donor_proportions[
    analysis_set == "primary_ucell",
    .(sample_id, n_total_cells)
  ])$n_total_cells)
  n_cells_input <- sum(unique(donor_proportions[
    analysis_set == "sensitivity_include_cluster13",
    .(sample_id, n_total_cells)
  ])$n_total_cells)
} else {
  object <- readRDS(input_file)
  metadata <- object@meta.data
  validate_metadata(metadata)
  counts <- SeuratObject::LayerData(object, assay = SeuratObject::DefaultAssay(object), layer = "counts")
  if (!inherits(counts, "sparseMatrix")) counts <- methods::as(counts, "dgCMatrix")

  module_map <- parse_module_map(contract$module_map)
  program_state <- build_programs(
    utils::read.delim(membership_file, check.names = FALSE),
    module_map,
    rownames(counts),
    as.integer(contract$minimum_matched_genes),
    as.numeric(contract$minimum_match_fraction)
  )
  if (!length(program_state$programs)) stop("No frozen program passed the scoring gate.")
  safe_write_tsv(program_state$audit, audit_file)

  cache_signature <- list(
    object = unname(tools::md5sum(input_file)),
    membership = unname(tools::md5sum(membership_file)),
    programs = program_state$programs
  )
  score_cache <- if (file.exists(score_cache_file)) readRDS(score_cache_file) else NULL
  cache_valid <- is.list(score_cache) && identical(score_cache$signature, cache_signature)
  if (cache_valid) {
    score_ucell <- score_cache$score_ucell
    score_logmean <- score_cache$score_logmean
  } else {
    score_ucell <- UCell::ScoreSignatures_UCell(
      counts,
      features = program_state$programs,
      ncores = 1L
    )
    score_ucell <- as.matrix(score_ucell[colnames(object), , drop = FALSE])
    colnames(score_ucell) <- sub("_UCell$", "", colnames(score_ucell))
    score_logmean <- calculate_log_normalized_mean(counts, program_state$programs)
    saveRDS(
      list(
        signature = cache_signature,
        score_ucell = score_ucell,
        score_logmean = score_logmean
      ),
      score_cache_file,
      compress = TRUE
    )
  }

  excluded_clusters <- split_values(contract$excluded_clusters_primary)
  primary_keep <- !as.character(metadata$seurat_clusters) %in% excluded_clusters
  analysis_sets <- list(
    primary_ucell = list(keep = primary_keep, scores = score_ucell),
    sensitivity_logmean = list(keep = primary_keep, scores = score_logmean),
    sensitivity_include_cluster13 = list(keep = rep(TRUE, nrow(metadata)), scores = score_ucell)
  )

  donor_scores <- data.table::rbindlist(lapply(names(analysis_sets), function(name) {
    state <- analysis_sets[[name]]
    aggregate_donor_scores(
      metadata[state$keep, , drop = FALSE],
      state$scores[state$keep, , drop = FALSE],
      as.integer(contract$minimum_cells_per_donor_celltype),
      name
    )
  }), use.names = TRUE, fill = TRUE)
  donor_proportions <- data.table::rbindlist(lapply(names(analysis_sets), function(name) {
    keep <- analysis_sets[[name]]$keep
    build_donor_proportions(metadata[keep, , drop = FALSE], name)
  }), use.names = TRUE, fill = TRUE)
  safe_write_tsv(donor_scores, scores_file)
  safe_write_tsv(donor_proportions, proportions_file)
  n_cells_input <- ncol(object)
  n_cells_primary <- sum(primary_keep)
}

