if (!nzchar(Sys.getenv("SLE_RUN_ROOT"))) stop("Run this analysis through run.py with an external workspace")
select_expression_matrix <- function(expression_table, sample_ids) {
  if (!data.table::is.data.table(expression_table)) {
    stop("expression_table must be a data.table.")
  }
  missing_samples <- setdiff(sample_ids, names(expression_table))
  if (length(missing_samples)) {
    stop("Missing expression sample column(s): ", paste(missing_samples, collapse = ", "))
  }
  matrix <- as.matrix(expression_table[, sample_ids, with = FALSE])
  storage.mode(matrix) <- "double"
  matrix
}

read_platform_gene_symbols <- function(path) {
  annotation <- data.table::fread(
    path,
    skip = 28L,
    header = TRUE,
    sep = "\t",
    quote = "",
    fill = Inf,
    select = c("ID", "Gene symbol"),
    showProgress = FALSE
  )
  data.table::setnames(annotation, c("ID", "Gene symbol"),
                       c("probe_id", "gene_symbol_raw"))
  annotation[, probe_id := as.character(probe_id)]
  annotation[, gene_symbol_raw := trimws(as.character(gene_symbol_raw))]
  annotation
}

prepare_plot_metadata <- function(data) {
  data <- data.table::copy(data.table::as.data.table(data))
  for (field in intersect(c("batch", "disease_state"), names(data))) {
    data[, (field) := factor(as.character(get(field)))]
  }
  data
}
