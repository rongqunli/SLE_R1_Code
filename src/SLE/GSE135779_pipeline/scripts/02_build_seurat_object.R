if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
source(file.path(dirname(normalizePath(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[[1]]), winslash = "/", mustWork = TRUE)), "_common.R"))

require_packages(
  c("Matrix", "data.table", "Seurat", "SeuratObject"),
  "building the GSE135779 Seurat object"
)

raw_dir <- file.path(cfg$paths$raw, cfg$files$extracted_dir)
genes_file <- file.path(cfg$paths$raw, cfg$files$genes)
author_metadata_file <- file.path(cfg$paths$metadata, cfg$files$author_metadata)
frozen_meta_file <- cfg$meta$frozen_master
output_file <- file.path(cfg$paths$processed, "00_raw_seurat.rds")
output_dir <- file.path(cfg$paths$results, "02_build_object")
checkpoint_dir <- file.path(cfg$paths$processed, "checkpoints", "02_build_object")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(genes_file, author_metadata_file, frozen_meta_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Missing input files: ", paste(missing_files, collapse = "; "),
    ". Run 01_download_data.R first."
  )
}
if (!dir.exists(raw_dir)) {
  stop("Extracted raw directory not found: ", raw_dir)
}

discover_sample_files <- function(path) {
  matrix_files <- sort(list.files(
    path,
    pattern = "^GSM[0-9]+_[^_]+_matrix\\.mtx\\.gz$",
    recursive = TRUE,
    full.names = TRUE
  ))
  if (length(matrix_files) == 0L) {
    stop("No GSM*_matrix.mtx.gz files found under: ", path)
  }

  parsed <- lapply(matrix_files, function(matrix_file) {
    file_name <- basename(matrix_file)
    match_parts <- regexec(
      "^(GSM[0-9]+)_([^_]+)_matrix\\.mtx\\.gz$",
      file_name
    )
    fields <- regmatches(file_name, match_parts)[[1]]
    if (length(fields) != 3L) {
      stop("Unexpected matrix filename: ", file_name)
    }
    prefix <- paste(fields[[2]], fields[[3]], sep = "_")
    barcode_file <- file.path(dirname(matrix_file), paste0(prefix, "_barcodes.tsv.gz"))
    if (!file.exists(barcode_file)) {
      stop("Barcode file missing for ", file_name, ": ", barcode_file)
    }
    data.frame(
      gsm = fields[[2]],
      library_id = fields[[3]],
      matrix_file = matrix_file,
      barcode_file = barcode_file,
      stringsAsFactors = FALSE
    )
  })
  files <- do.call(rbind, parsed)
  if (anyDuplicated(files$gsm) > 0L || anyDuplicated(files$library_id) > 0L) {
    stop("Duplicate GSM or library IDs detected in extracted files.")
  }
  files
}

author_major_type <- function(subcluster) {
  subcluster <- as.character(subcluster)
  out <- rep(NA_character_, length(subcluster))
  out[grepl("^T-", subcluster)] <- "T"
  out[grepl("^Mono-", subcluster)] <- "Monocyte"
  out[grepl("^NK-", subcluster)] <- "NK"
  out[grepl("^B-", subcluster)] <- "B"
  out[grepl("^cDC($|-)", subcluster, ignore.case = TRUE)] <- "cDC"
  out[grepl("^pDC($|-)", subcluster, ignore.case = TRUE)] <- "pDC"
  out[grepl("^PC($|-)", subcluster, ignore.case = TRUE)] <- "Plasmablast"
  out[tolower(subcluster) == "mgk"] <- "Platelet"
  out[tolower(subcluster) == "eryth"] <- "Erythroid"
  out
}

collapse_gene_symbols <- function(count_matrix, gene_symbols) {
  if (nrow(count_matrix) != length(gene_symbols)) {
    stop("Gene annotation length does not match matrix rows.")
  }
  unique_symbols <- unique(gene_symbols)
  symbol_index <- match(gene_symbols, unique_symbols)
  row_mapper <- Matrix::sparseMatrix(
    i = symbol_index,
    j = seq_along(symbol_index),
    x = 1,
    dims = c(length(unique_symbols), length(gene_symbols))
  )
  collapsed <- row_mapper %*% count_matrix
  collapsed <- methods::as(collapsed, "dgCMatrix")
  rownames(collapsed) <- unique_symbols
  collapsed
}

sample_files <- discover_sample_files(raw_dir)
log_message("Detected ", nrow(sample_files), " GEO count matrices.")
if (nrow(sample_files) < cfg$data$expected_geo_matrices) {
  stop(
    "Detected only ", nrow(sample_files), " matrices; expected at least ",
    cfg$data$expected_geo_matrices, ". The raw archive may be incomplete."
  )
}

step_signature <- build_step_signature(
  input_files = c(
    genes_file, author_metadata_file, frozen_meta_file,
    sample_files$matrix_file, sample_files$barcode_file
  ),
  extra = list(step = "02_build_object", expected_matrices = cfg$data$expected_geo_matrices)
)
cache_root_signature <- build_step_signature(
  input_files = c(genes_file, author_metadata_file, frozen_meta_file),
  extra = list(
    step = "02_build_object_sample_cache",
    expected_matrices = cfg$data$expected_geo_matrices
  )
)
initialize_step(
  cache_root_signature,
  checkpoint_dir,
  owned_paths = character(0)
)
initialize_step(
  step_signature,
  file.path(checkpoint_dir, "_final_state"),
  owned_paths = c(output_dir, output_file)
)
if (step_is_complete(output_dir, step_signature, output_file)) {
  log_message("Step 02 is already complete and its signature matches. Use --force to rebuild.")
  quit(save = "no", status = 0L)
}

genes <- data.table::fread(
  genes_file,
  sep = "\t",
  header = FALSE,
  col.names = c("ensembl_id", "gene_symbol"),
  data.table = FALSE
)
if (ncol(genes) < 2L || nrow(genes) == 0L) {
  stop("Invalid gene annotation file: ", genes_file)
}
genes$gene_symbol <- trimws(genes$gene_symbol)
genes$gene_symbol[!nzchar(genes$gene_symbol)] <- genes$ensembl_id[!nzchar(genes$gene_symbol)]

gene_mapping <- data.frame(
  ensembl_id = genes$ensembl_id,
  gene_symbol = genes$gene_symbol,
  collapsed_symbol = genes$gene_symbol,
  duplicated_symbol = duplicated(genes$gene_symbol) |
    duplicated(genes$gene_symbol, fromLast = TRUE),
  stringsAsFactors = FALSE
)
safe_write_tsv(gene_mapping, file.path(output_dir, "gene_annotation_mapping.tsv"))

author_meta <- data.table::fread(
  author_metadata_file,
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)
required_author_columns <- c("index", "IDs", "Names", "SLEDAI", "subclusters")
if (!all(required_author_columns %in% colnames(author_meta))) {
  stop(
    "Author metadata is missing columns: ",
    paste(setdiff(required_author_columns, colnames(author_meta)), collapse = ", ")
  )
}

author_cells <- data.frame(
  library_id = trimws(author_meta$IDs),
  subject_id_author = trimws(author_meta$Names),
  sle_activity = trimws(author_meta$SLEDAI),
  author_subcluster = trimws(author_meta$subclusters),
  barcode_base = strip_barcode_suffix(author_meta$index),
  stringsAsFactors = FALSE
)
author_cells$author_cell_type_major <- author_major_type(author_cells$author_subcluster)
author_cells$join_key <- paste(author_cells$library_id, author_cells$barcode_base, sep = "|")
if (anyDuplicated(author_cells$join_key) > 0L) {
  duplicates <- unique(author_cells$join_key[duplicated(author_cells$join_key)])
  stop("Duplicate library/barcode keys in author metadata: ", paste(head(duplicates), collapse = ", "))
}

author_sample <- unique(author_cells[c(
  "library_id", "subject_id_author", "sle_activity"
)])




frozen_meta <- data.table::fread(
  frozen_meta_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  colClasses = "character",
  check.names = FALSE
)
required_meta_columns <- c(
  "meta_version", "dataset", "sample_id", "patient_id", "visit_id",
  "timepoint", "disease", "disease_state", "disease_activity_group",
  "response", "tissue", "treatment", "treatment_stage", "SLEDAI",
  "renal_activity", "C3", "C4", "anti_dsDNA", "age", "sex", "cohort",
  "batch", "platform", "assay_type", "sc_library", "cell_count",
  "analysis_include"
)
if (!all(required_meta_columns %in% colnames(frozen_meta))) {
  stop(
    "Frozen metadata is missing columns: ",
    paste(setdiff(required_meta_columns, colnames(frozen_meta)), collapse = ", ")
  )
}
frozen_meta <- frozen_meta[frozen_meta$dataset == cfg$accession, , drop = FALSE]
if (nrow(frozen_meta) != cfg$data$expected_geo_matrices) {
  stop(
    "Frozen metadata has ", nrow(frozen_meta), " GSE135779 rows; expected ",
    cfg$data$expected_geo_matrices, "."
  )
}
if (anyDuplicated(frozen_meta$sample_id) > 0L ||
    anyDuplicated(frozen_meta$sc_library) > 0L) {
  stop("Duplicate sample_id or sc_library in frozen GSE135779 metadata.")
}
if (!all(frozen_meta$meta_version == cfg$meta$version)) {
  stop("Frozen metadata version does not match cfg$meta$version.")
}
if (!all(frozen_meta$analysis_include == "yes")) {
  stop("GSE135779 frozen metadata contains analysis_include != yes.")
}

missing_gsm <- setdiff(sample_files$gsm, frozen_meta$sample_id)
extra_gsm <- setdiff(frozen_meta$sample_id, sample_files$gsm)
if (length(missing_gsm) > 0L || length(extra_gsm) > 0L) {
  stop(
    "Count matrices and frozen metadata have different GSM sets. Missing in meta: ",
    paste(missing_gsm, collapse = ", "), "; extra in meta: ",
    paste(extra_gsm, collapse = ", ")
  )
}

meta_for_merge <- frozen_meta
meta_for_merge$gsm <- meta_for_merge$sample_id
meta_for_merge$library_id <- meta_for_merge$sc_library
meta_for_merge$subject_id <- meta_for_merge$patient_id
meta_for_merge$condition <- ifelse(
  meta_for_merge$disease == "Healthy control", "HD", "SLE"
)
meta_for_merge$sle_activity <- meta_for_merge$disease_activity_group

sample_table <- merge(
  sample_files,
  meta_for_merge,
  by = c("gsm", "library_id"),
  all.x = TRUE,
  sort = FALSE
)
if (nrow(sample_table) != nrow(sample_files) ||
    any(is.na(sample_table$patient_id)) ||
    any(is.na(sample_table$condition))) {
  stop("Frozen metadata did not join one-to-one to matrix GSM/library IDs.")
}

sample_table <- sample_table[order(sample_table$cohort, sample_table$condition, sample_table$gsm), ]
rownames(sample_table) <- NULL
safe_write_tsv(sample_table, file.path(cfg$paths$metadata, "sample_metadata.tsv"))

author_library_status <- merge(
  author_sample,
  frozen_meta[c("sc_library", "patient_id", "disease_activity_group")],
  by.x = "library_id",
  by.y = "sc_library",
  all = TRUE
)
author_library_status$present_in_geo_files <- author_library_status$library_id %in%
  sample_table$library_id
author_library_status$subject_matches_frozen <- with(
  author_library_status,
  !is.na(subject_id_author) & !is.na(patient_id) & subject_id_author == patient_id
)
author_library_status$clinical_source_used <- "frozen_meta"
safe_write_tsv(
  author_library_status,
  file.path(output_dir, "author_library_match_status.tsv")
)

matrix_list <- vector("list", nrow(sample_table))
metadata_list <- vector("list", nrow(sample_table))
join_summary <- vector("list", nrow(sample_table))

for (i in seq_len(nrow(sample_table))) {
  row <- sample_table[i, , drop = FALSE]
  sample_cache <- file.path(
    checkpoint_dir,
    paste0("sample_", make_safe_name(row$gsm), "_", make_safe_name(row$library_id), ".rds")
  )
  sample_signature <- build_step_signature(
    input_files = c(
      row$matrix_file, row$barcode_file, genes_file,
      author_metadata_file, frozen_meta_file
    ),
    extra = list(
      step = "02_build_object_sample",
      gsm = row$gsm,
      library_id = row$library_id,
      patient_id = row$patient_id
    )
  )
  cached <- load_checkpoint_rds(sample_cache, sample_signature)
  if (!is.null(cached)) {
    matrix_list[[i]] <- cached$counts
    metadata_list[[i]] <- cached$metadata
    join_summary[[i]] <- cached$join_summary
    next
  }
  log_message(
    "Reading ", row$gsm, " / ", row$library_id, " / ", row$patient_id,
    " (", i, "/", nrow(sample_table), ")"
  )

  barcodes <- data.table::fread(
    row$barcode_file,
    sep = "\t",
    header = FALSE,
    data.table = FALSE
  )[[1]]
  barcodes <- trimws(as.character(barcodes))
  barcode_base <- strip_barcode_suffix(barcodes)

  matrix_con <- gzfile(row$matrix_file, open = "rt")
  raw_matrix <- Matrix::readMM(matrix_con)
  close(matrix_con)
  raw_matrix <- methods::as(raw_matrix, "dgCMatrix")

  if (nrow(raw_matrix) != nrow(genes)) {
    stop("Matrix/gene row mismatch for ", row$gsm)
  }
  if (ncol(raw_matrix) != length(barcodes)) {
    stop("Matrix/barcode column mismatch for ", row$gsm)
  }

  counts <- collapse_gene_symbols(raw_matrix, genes$gene_symbol)
  cell_id <- paste(row$gsm, row$library_id, barcode_base, sep = "_")
  if (anyDuplicated(cell_id) > 0L) {
    stop("Duplicate cell IDs within sample: ", row$gsm)
  }
  colnames(counts) <- cell_id

  cell_join_key <- paste(row$library_id, barcode_base, sep = "|")
  author_index <- match(cell_join_key, author_cells$join_key)
  cell_meta <- data.frame(
    cell_id = cell_id,
    barcode = barcodes,
    barcode_base = barcode_base,
    sample_id = row$sample_id,
    gsm = row$gsm,
    library_id = row$library_id,
    subject_id = row$patient_id,
    patient_id = row$patient_id,
    visit_id = row$visit_id,
    timepoint = row$timepoint,
    condition = row$condition,
    disease = row$disease,
    disease_state = row$disease_state,
    cohort = row$cohort,
    age = row$age,
    sex = row$sex,
    response = row$response,
    tissue = row$tissue,
    treatment = row$treatment,
    treatment_stage = row$treatment_stage,
    SLEDAI = row$SLEDAI,
    sle_activity = row$disease_activity_group,
    renal_activity = row$renal_activity,
    C3 = row$C3,
    C4 = row$C4,
    anti_dsDNA = row$anti_dsDNA,
    batch = row$batch,
    platform = row$platform,
    meta_version = row$meta_version,
    analysis_include = row$analysis_include,
    author_subcluster = author_cells$author_subcluster[author_index],
    author_cell_type_major = author_cells$author_cell_type_major[author_index],
    author_label_matched = !is.na(author_index),
    stringsAsFactors = FALSE,
    row.names = cell_id
  )

  matrix_list[[i]] <- counts
  metadata_list[[i]] <- cell_meta
  join_summary_row <- data.frame(
    sample_id = row$sample_id,
    library_id = row$library_id,
    subject_id = row$patient_id,
    matrix_cells = ncol(counts),
    author_matched_cells = sum(cell_meta$author_label_matched),
    author_match_fraction = mean(cell_meta$author_label_matched),
    stringsAsFactors = FALSE
  )
  join_summary[[i]] <- join_summary_row
  save_checkpoint_rds(
    list(
      counts = matrix_list[[i]],
      metadata = metadata_list[[i]],
      join_summary = join_summary_row
    ),
    sample_cache,
    sample_signature,
    compress = FALSE
  )

  rm(raw_matrix, counts, cell_meta, join_summary_row)
  if (i %% 5L == 0L) gc(FALSE)
}

log_message("Combining sparse matrices.")
combined_counts <- do.call(cbind, matrix_list)
combined_counts <- methods::as(combined_counts, "dgCMatrix")
combined_metadata <- do.call(rbind, metadata_list)
combined_metadata <- combined_metadata[colnames(combined_counts), , drop = FALSE]

if (!identical(rownames(combined_metadata), colnames(combined_counts))) {
  stop("Cell metadata row order does not match count matrix columns.")
}

object <- SeuratObject::CreateSeuratObject(
  counts = combined_counts,
  project = cfg$accession,
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  meta.data = combined_metadata
)
object$condition <- factor(object$condition, levels = c("HD", "SLE"))
object$cohort <- factor(object$cohort, levels = c("Child", "Adult"))

join_summary <- do.call(rbind, join_summary)
safe_write_tsv(join_summary, file.path(output_dir, "author_cell_label_join_summary.tsv"))

object_summary <- data.frame(
  metric = c(
    "genes_original", "gene_symbols_after_collapse", "cells", "samples",
    "subjects", "author_labels_matched"
  ),
  value = c(
    nrow(genes), nrow(object), ncol(object), length(unique(object$sample_id)),
    length(unique(object$subject_id)), sum(object$author_label_matched)
  ),
  stringsAsFactors = FALSE
)
safe_write_tsv(object_summary, file.path(output_dir, "object_summary.tsv"))

atomic_save_rds(object, output_file, compress = FALSE)
write_session_info(file.path(cfg$paths$logs, "02_build_seurat_sessionInfo.txt"))
mark_step_complete(output_dir, step_signature, output_file)

log_message("Raw Seurat object saved: ", output_file)
log_message("Dimensions: ", nrow(object), " genes x ", ncol(object), " cells")
