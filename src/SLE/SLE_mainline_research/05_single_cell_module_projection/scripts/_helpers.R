if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
read_contract <- function(path) {
  x <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(colnames(x), c("key", "value")) || anyDuplicated(x$key)) {
    stop("Analysis contract must contain unique key/value rows.")
  }
  stats::setNames(as.list(x$value), x$key)
}

split_values <- function(x, separator = ";") {
  values <- trimws(strsplit(x %||% "", separator, fixed = TRUE)[[1]])
  values[nzchar(values)]
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

resolve_contract_path <- function(study_root, path) {
  normalizePath(file.path(study_root, "config", path), winslash = "/", mustWork = TRUE)
}

parse_module_map <- function(x) {
  pieces <- split_values(x)
  pair_list <- strsplit(pieces, "=", fixed = TRUE)
  if (any(lengths(pair_list) != 2L)) stop("module_map entries must use source=program.")
  source <- vapply(pair_list, `[[`, character(1), 1L)
  program <- vapply(pair_list, `[[`, character(1), 2L)
  if (anyDuplicated(source) || anyDuplicated(program)) stop("module_map entries must be unique.")
  stats::setNames(program, source)
}

numeric_or_na <- function(x) suppressWarnings(as.numeric(as.character(x)))

safe_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  data.table::fwrite(x, tmp, sep = "\t", na = "NA", quote = FALSE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not atomically write ", path)
  invisible(path)
}

validate_metadata <- function(metadata) {
  required <- c(
    "sample_id", "subject_id", "cohort", "condition", "age", "sex", "batch",
    "seurat_clusters", "cell_type_analysis"
  )
  missing <- setdiff(required, colnames(metadata))
  if (length(missing)) stop("Object metadata is missing: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

build_programs <- function(membership, module_map, available_genes, min_genes, min_fraction) {
  required <- c("gene", "module")
  if (!all(required %in% colnames(membership))) stop("Module membership lacks gene/module columns.")
  programs <- lapply(names(module_map), function(module) {
    unique(as.character(membership$gene[membership$module == module]))
  })
  names(programs) <- unname(module_map)
  audit <- data.frame(
    source_module = names(module_map),
    program = unname(module_map),
    requested_genes = lengths(programs),
    stringsAsFactors = FALSE
  )
  matched <- lapply(programs, intersect, y = available_genes)
  audit$matched_genes <- lengths(matched)
  audit$match_fraction <- audit$matched_genes / audit$requested_genes
  audit$matched_gene_list <- vapply(matched, paste, character(1), collapse = ";")
  audit$missing_gene_list <- vapply(
    names(programs),
    function(program) paste(setdiff(programs[[program]], matched[[program]]), collapse = ";"),
    character(1)
  )
  audit$passes_scoring_gate <- audit$matched_genes >= min_genes &
    audit$match_fraction >= min_fraction
  list(programs = matched[audit$passes_scoring_gate], audit = audit)
}

calculate_log_normalized_mean <- function(counts, programs, scale_factor = 1e4) {
  library_size <- Matrix::colSums(counts)
  library_size[library_size <= 0] <- 1
  normalized <- counts %*% Matrix::Diagonal(x = scale_factor / library_size)
  normalized@x <- log1p(normalized@x)
  result <- vapply(programs, function(genes) {
    Matrix::colMeans(normalized[genes, , drop = FALSE])
  }, numeric(ncol(counts)))
  if (is.null(dim(result))) result <- matrix(result, ncol = 1L)
  rownames(result) <- colnames(counts)
  colnames(result) <- names(programs)
  result
}

aggregate_donor_scores <- function(metadata, score_matrix, min_cells, analysis_set) {
  dt <- data.table::as.data.table(data.frame(
    metadata,
    score_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  ))
  score_columns <- colnames(score_matrix)
  long <- data.table::melt(
    dt,
    id.vars = setdiff(colnames(dt), score_columns),
    measure.vars = score_columns,
    variable.name = "program",
    value.name = "score"
  )
  output <- long[
    !is.na(cell_type_analysis) & nzchar(cell_type_analysis),
    .(
      n_cells = .N,
      mean_score = mean(score, na.rm = TRUE),
      median_score = stats::median(score, na.rm = TRUE)
    ),
    by = .(
      sample_id, subject_id, cohort, condition, age, sex, batch,
      cell_type_analysis, program
    )
  ]
  output <- output[n_cells >= min_cells]
  output[, analysis_set := analysis_set]
  output[]
}

build_donor_proportions <- function(metadata, analysis_set) {
  dt <- data.table::as.data.table(metadata)
  counts <- dt[
    !is.na(cell_type_analysis) & nzchar(cell_type_analysis),
    .(n_cells = .N),
    by = .(
      sample_id, subject_id, cohort, condition, age, sex, batch,
      cell_type_analysis
    )
  ]
  totals <- counts[, .(n_total_cells = sum(n_cells)), by = sample_id]
  counts <- merge(counts, totals, by = "sample_id", all.x = TRUE, sort = FALSE)
  counts[, cell_proportion := n_cells / n_total_cells]
  counts[, analysis_set := analysis_set]
  counts[]
}

