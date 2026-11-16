if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
source(file.path(dirname(normalizePath(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[[1]]), winslash = "/", mustWork = TRUE)), "_common.R"))

require_packages(c("data.table"), "downloading GSE135779 support files")

raw_tar <- file.path(cfg$paths$raw, cfg$files$raw_tar)
genes_file <- file.path(cfg$paths$raw, cfg$files$genes)
author_metadata <- file.path(cfg$paths$metadata, cfg$files$author_metadata)
extract_dir <- file.path(cfg$paths$raw, cfg$files$extracted_dir)

download_if_missing <- function(url, destination, mode = "wb") {
  if (file.exists(destination) && file.info(destination)$size > 0) {
    log_message("Reuse existing file: ", destination)
    return(invisible(destination))
  }
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  log_message("Downloading: ", url)
  utils::download.file(url, destination, mode = mode, method = "libcurl", quiet = FALSE)
  if (!file.exists(destination) || file.info(destination)$size == 0) {
    stop("Download failed: ", destination)
  }
  invisible(destination)
}

download_if_missing(cfg$urls$genes, genes_file)
download_if_missing(cfg$urls$author_metadata, author_metadata)

if (has_flag("--download-raw")) {
  download_if_missing(cfg$urls$raw_tar, raw_tar)
} else if (!file.exists(raw_tar)) {
  log_message("Raw tar not downloaded. Re-run with --download-raw or place it at: ", raw_tar)
}

if (file.exists(raw_tar)) {
  matrix_files_existing <- if (dir.exists(extract_dir)) {
    list.files(extract_dir, pattern = "_matrix\\.mtx\\.gz$", recursive = TRUE, full.names = TRUE)
  } else {
    character(0)
  }

  if (length(matrix_files_existing) == 0L) {
    dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
    log_message("Extracting ", raw_tar, " to ", extract_dir)
    utils::untar(raw_tar, exdir = extract_dir)
  } else {
    log_message("Reuse extracted directory with ", length(matrix_files_existing), " matrices: ", extract_dir)
  }
}

gene_preview <- data.table::fread(
  genes_file,
  sep = "\t",
  header = FALSE,
  nrows = 5L,
  data.table = FALSE
)
if (ncol(gene_preview) < 2L) {
  stop("Unexpected GSE135779 gene file format: ", genes_file)
}

author_preview <- data.table::fread(
  author_metadata,
  header = TRUE,
  nrows = 5L,
  data.table = FALSE,
  check.names = FALSE
)
required_author_columns <- c("index", "IDs", "Names", "SLEDAI", "subclusters")
if (!all(required_author_columns %in% colnames(author_preview))) {
  stop(
    "Unexpected author metadata format. Missing columns: ",
    paste(setdiff(required_author_columns, colnames(author_preview)), collapse = ", ")
  )
}

download_manifest <- data.frame(
  file = c(raw_tar, genes_file, author_metadata),
  exists = file.exists(c(raw_tar, genes_file, author_metadata)),
  size_bytes = vapply(c(raw_tar, genes_file, author_metadata), function(path) {
    if (file.exists(path)) file.info(path)$size else NA_real_
  }, numeric(1)),
  stringsAsFactors = FALSE
)
safe_write_tsv(download_manifest, file.path(cfg$paths$metadata, "download_manifest.tsv"))

log_message("Download/preparation step completed.")
