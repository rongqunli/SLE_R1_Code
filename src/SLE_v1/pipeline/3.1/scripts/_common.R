if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

script_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0L) {
    return(NA_character_)
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
}

pipeline_root_from_script <- function() {
  current_script <- script_file()
  if (is.na(current_script)) {
    candidate <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    if (basename(candidate) == "scripts") {
      return(dirname(candidate))
    }
    if (file.exists(file.path(candidate, "config", "config.R"))) {
      return(candidate)
    }
    stop("Cannot determine pipeline root. Run this file with Rscript.")
  }
  script_dir <- dirname(current_script)
  if (basename(script_dir) == "scripts") {
    return(normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE))
  }
  if (basename(script_dir) == "tests") {
    return(normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE))
  }
  if (file.exists(file.path(script_dir, "config", "config.R"))) {
    return(script_dir)
  }
  stop("Cannot determine pipeline root from: ", current_script)
}

PIPELINE_ROOT <- pipeline_root_from_script()
source(file.path(PIPELINE_ROOT, "config", "config.R"), local = FALSE)

ensure_dirs <- function(paths) {
  for (path in unname(unlist(paths))) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }
  }
  invisible(paths)
}

ensure_dirs(cfg$paths)

require_packages <- function(packages, purpose = NULL) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    detail <- if (is.null(purpose)) "" else paste0(" for ", purpose)
    stop(
      "Missing R package(s)", detail, ": ", paste(missing, collapse = ", "),
      ". Run scripts/00_install_packages.R first."
    )
  }
  invisible(TRUE)
}

has_flag <- function(flag) {
  flag %in% commandArgs(trailingOnly = TRUE)
}

safe_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  data.table::fwrite(x, temporary, sep = "\t", quote = FALSE, na = "NA")
  if (file.exists(path)) unlink(path, force = TRUE)
  if (!file.rename(temporary, path)) {
    stop("Could not atomically move temporary TSV to: ", path)
  }
  invisible(path)
}

safe_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  data.table::fwrite(x, temporary, sep = ",", quote = TRUE, na = "NA")
  if (file.exists(path)) unlink(path, force = TRUE)
  if (!file.rename(temporary, path)) {
    stop("Could not atomically move temporary CSV to: ", path)
  }
  invisible(path)
}

force_requested <- function() {
  has_flag("--force")
}

atomic_save_rds <- function(object, path, compress = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  saveRDS(object, temporary, compress = compress)
  if (file.exists(path)) unlink(path, force = TRUE)
  if (!file.rename(temporary, path)) {
    stop("Could not atomically move temporary RDS to: ", path)
  }
  invisible(path)
}

path_signature <- function(path, hash_contents = FALSE) {
  normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (!file.exists(path)) {
    return(list(path = normalized, exists = FALSE))
  }
  info <- file.info(path)
  signature <- list(
    path = normalized,
    exists = TRUE,
    size = unname(info$size),
    modified = format(info$mtime, "%Y-%m-%dT%H:%M:%OS6%z")
  )
  if (isTRUE(hash_contents)) {
    signature$md5 <- unname(tools::md5sum(path))
  }
  signature
}

build_step_signature <- function(
    input_files = character(0),
    extra = list(),
    code_files = NULL) {
  if (is.null(code_files)) {
    current_script <- script_file()
    code_files <- c(
      current_script,
      file.path(PIPELINE_ROOT, "scripts", "_common.R"),
      file.path(PIPELINE_ROOT, "config", "config.R")
    )
  }
  code_files <- unique(code_files)
  code_files <- code_files[!is.na(code_files) & nzchar(code_files)]
  input_files <- unique(as.character(input_files))
  list(
    checkpoint_schema = 1L,
    code = lapply(code_files, path_signature, hash_contents = TRUE),
    inputs = lapply(input_files, path_signature, hash_contents = FALSE),
    extra = extra
  )
}

signature_matches <- function(path, signature) {
  if (!file.exists(path)) return(FALSE)
  saved <- tryCatch(readRDS(path), error = function(e) NULL)
  identical(success_record_signature(saved), signature)
}

success_record_signature <- function(record) {
  if (is.list(record) && !is.null(record$signature)) record$signature else record
}

normalize_for_comparison <- function(path, must_work = FALSE) {
  tolower(normalizePath(path, winslash = "/", mustWork = must_work))
}

path_is_within <- function(path, parent) {
  child_normalized <- normalize_for_comparison(path, must_work = FALSE)
  parent_normalized <- normalize_for_comparison(parent, must_work = FALSE)
  identical(child_normalized, parent_normalized) ||
    startsWith(child_normalized, paste0(parent_normalized, "/"))
}

protected_cleanup_roots <- function() {
  unique(c(
    PIPELINE_ROOT,
    cfg$paths$raw,
    cfg$paths$metadata,
    file.path(PIPELINE_ROOT, "data", "frozen"),
    file.path(cfg$paths$processed, "checkpoints", "02_build_object"),
    file.path(cfg$paths$results, "02_build_object")
  ))
}

is_protected_cleanup_path <- function(path) {
  path_normalized <- normalize_for_comparison(path, must_work = FALSE)
  protected_file <- normalize_for_comparison(
    file.path(cfg$paths$processed, "00_raw_seurat.rds"),
    must_work = FALSE
  )
  if (identical(path_normalized, protected_file)) return(TRUE)

  roots <- protected_cleanup_roots()
  any(vapply(roots, function(root) {
    root_normalized <- normalize_for_comparison(root, must_work = FALSE)
    if (identical(root_normalized, normalize_for_comparison(PIPELINE_ROOT, FALSE))) {
      identical(path_normalized, root_normalized)
    } else {
      path_is_within(path, root)
    }
  }, logical(1)))
}

validate_cleanup_path <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("Cleanup path must be one non-empty path.")
  }
  if (!path_is_within(path, PIPELINE_ROOT)) {
    stop("Refusing to clean a path outside PIPELINE_ROOT: ", path)
  }
  if (is_protected_cleanup_path(path)) {
    stop("Refusing to clean a protected pipeline path: ", path)
  }
  invisible(normalize_for_comparison(path, must_work = FALSE))
}

path_size_bytes <- function(path) {
  if (file.exists(path) && !dir.exists(path)) {
    return(as.numeric(file.info(path)$size))
  }
  if (!dir.exists(path)) return(0)
  entries <- list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  entries <- entries[file.exists(entries) & !dir.exists(entries)]
  if (length(entries) == 0L) return(0)
  sum(as.numeric(file.info(entries)$size), na.rm = TRUE)
}

cleanup_owned_paths <- function(paths, reason = "final outputs validated") {
  paths <- unique(as.character(paths))
  rows <- lapply(paths, function(path) {
    existed <- file.exists(path) || dir.exists(path)
    bytes <- if (existed) path_size_bytes(path) else 0
    status <- "not_found"
    detail <- "Path did not exist."
    if (existed) {
      validate_cleanup_path(path)
      unlink(path, recursive = TRUE, force = TRUE)
      removed <- !file.exists(path) && !dir.exists(path)
      status <- if (removed) "removed" else "failed"
      detail <- if (removed) reason else "unlink() did not remove the path."
    }
    data.frame(
      path = normalizePath(path, winslash = "/", mustWork = FALSE),
      status = status,
      bytes = bytes,
      reason = detail,
      stringsAsFactors = FALSE
    )
  })
  if (length(rows) == 0L) {
    return(data.frame(
      path = character(0), status = character(0), bytes = numeric(0),
      reason = character(0), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

reset_owned_paths <- function(paths) {
  for (path in unique(paths)) {
    if (!file.exists(path) && !dir.exists(path)) next
    normalized <- normalizePath(path, winslash = "/", mustWork = TRUE)
    pipeline_normalized <- normalizePath(PIPELINE_ROOT, winslash = "/", mustWork = TRUE)
    if (!startsWith(tolower(normalized), paste0(tolower(pipeline_normalized), "/"))) {
      stop("Refusing to reset a path outside PIPELINE_ROOT: ", normalized)
    }
    unlink(path, recursive = TRUE, force = TRUE)
  }
  invisible(paths)
}

initialize_step <- function(
    signature,
    checkpoint_dir,
    owned_paths = character(0),
    force = force_requested()) {
  signature_file <- file.path(checkpoint_dir, "_SIGNATURE.rds")
  stale_partial <- dir.exists(checkpoint_dir) &&
    !signature_matches(signature_file, signature)
  if (isTRUE(force) || stale_partial) {
    reason <- if (isTRUE(force)) "--force requested" else "step signature changed"
    log_message("Resetting step-owned outputs because ", reason, ".")
    reset_owned_paths(c(owned_paths, checkpoint_dir))
  }
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  if (!signature_matches(signature_file, signature)) {
    atomic_save_rds(signature, signature_file, compress = TRUE)
  }
  invisible(signature_file)
}

checkpoint_is_valid <- function(path, signature) {
  file.exists(path) && signature_matches(paste0(path, ".signature.rds"), signature)
}

save_checkpoint_rds <- function(object, path, signature, compress = FALSE) {
  atomic_save_rds(object, path, compress = compress)
  atomic_save_rds(signature, paste0(path, ".signature.rds"), compress = TRUE)
  invisible(path)
}

load_checkpoint_rds <- function(path, signature) {
  if (!checkpoint_is_valid(path, signature)) return(NULL)
  log_message("Reusing checkpoint: ", path)
  readRDS(path)
}

step_success_file <- function(output_dir) {
  file.path(output_dir, "_SUCCESS.rds")
}

step_is_complete <- function(output_dir, signature, required_outputs) {
  success_file <- step_success_file(output_dir)
  if (!file.exists(success_file) || !all(file.exists(required_outputs))) return(FALSE)
  record <- tryCatch(readRDS(success_file), error = function(e) NULL)
  if (!identical(success_record_signature(record), signature)) return(FALSE)
  if (is.list(record) && !is.null(record$outputs)) {
    recorded_paths <- vapply(record$outputs, function(output) {
      as.character(output$path %||% NA_character_)
    }, character(1))
    recorded_paths <- recorded_paths[!is.na(recorded_paths) & nzchar(recorded_paths)]
    if (!all(file.exists(recorded_paths))) return(FALSE)
  }
  TRUE
}

mark_step_complete <- function(output_dir, signature, outputs) {
  outputs <- unique(as.character(outputs))
  missing_outputs <- outputs[!file.exists(outputs)]
  if (length(missing_outputs) > 0L) {
    stop(
      "Cannot mark step complete; required output(s) are missing: ",
      paste(missing_outputs, collapse = ", ")
    )
  }
  record <- list(
    schema_version = 2L,
    signature = signature,
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    outputs = lapply(outputs, path_signature, hash_contents = FALSE)
  )
  atomic_save_rds(record, step_success_file(output_dir), compress = TRUE)
  legacy_record <- file.path(output_dir, "_COMPLETION_RECORD.rds")
  if (file.exists(legacy_record)) unlink(legacy_record, force = TRUE)
  safe_write_tsv(flatten_config(cfg), file.path(output_dir, "_ANALYSIS_CONFIG.tsv"))
  invisible(record)
}

complete_step_and_cleanup <- function(
    output_dir,
    signature,
    outputs,
    cleanup_paths = character(0),
    cleanup_report = NULL) {
  record <- mark_step_complete(output_dir, signature, outputs)
  cleanup <- cleanup_owned_paths(cleanup_paths)
  if (!is.null(cleanup_report)) safe_write_tsv(cleanup, cleanup_report)
  invisible(list(record = record, cleanup = cleanup))
}

combine_result_frames <- function(results) {
  results <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, results)
  if (length(results) == 0L) return(data.frame())
  data.table::rbindlist(results, fill = TRUE, use.names = TRUE) |>
    as.data.frame(stringsAsFactors = FALSE)
}

named_cache_has_fields <- function(cache, expected_names, fields) {
  is.list(cache) &&
    all(expected_names %in% names(cache)) &&
    all(vapply(cache[expected_names], function(item) {
      is.list(item) && all(fields %in% names(item))
    }, logical(1)))
}

numeric_or_na <- function(x) {
  value <- suppressWarnings(as.numeric(as.character(x)))
  value[!is.finite(value)] <- NA_real_
  value
}

derive_sledai_group <- function(x) {
  value <- numeric_or_na(x)
  ifelse(
    is.na(value), NA_character_,
    ifelse(
      value <= cfg$clinical$activity_low_max,
      "Low",
      ifelse(value > cfg$clinical$activity_high_min_exclusive, "High", "Intermediate")
    )
  )
}

flatten_config <- function(x, prefix = "cfg") {
  rows <- list()
  walk <- function(value, path) {
    if (is.list(value)) {
      for (name in names(value)) walk(value[[name]], paste(path, name, sep = "."))
    } else {
      rows[[length(rows) + 1L]] <<- data.frame(
        parameter = path,
        value = paste(as.character(value), collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }
  walk(x, prefix)
  do.call(rbind, rows)
}

get_assay_counts <- function(object, assay = "RNA") {
  tryCatch(
    SeuratObject::GetAssayData(object, assay = assay, layer = "counts"),
    error = function(e) SeuratObject::GetAssayData(object, assay = assay, slot = "counts")
  )
}

get_assay_data <- function(object, assay = "RNA") {
  tryCatch(
    SeuratObject::GetAssayData(object, assay = assay, layer = "data"),
    error = function(e) SeuratObject::GetAssayData(object, assay = assay, slot = "data")
  )
}

read_program_table <- function(path) {
  require_packages("data.table", "reading gene sets")
  x <- data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE)
  if (!all(c("program", "genes") %in% colnames(x))) {
    stop("Gene set file must contain program and genes columns: ", path)
  }
  gene_sets <- lapply(x$genes, function(value) {
    genes <- unique(trimws(unlist(strsplit(value, "[;,|]"))))
    genes[nzchar(genes)]
  })
  names(gene_sets) <- x$program
  gene_sets
}

read_marker_table <- function(path) {
  require_packages("data.table", "reading marker sets")
  x <- data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE)
  if (!all(c("cell_type", "marker_genes") %in% colnames(x))) {
    stop("Marker file must contain cell_type and marker_genes columns: ", path)
  }
  marker_sets <- lapply(x$marker_genes, function(value) {
    genes <- unique(trimws(unlist(strsplit(value, "[;,|]"))))
    genes[nzchar(genes)]
  })
  names(marker_sets) <- x$cell_type
  marker_sets
}

strip_barcode_suffix <- function(x) {
  sub("[-_.][0-9]+$", "", trimws(as.character(x)))
}

make_safe_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
}

write_session_info <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(capture.output(sessionInfo()), con)
  invisible(path)
}

log_message <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
}

set.seed(cfg$seed)
